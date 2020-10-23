use libc::{c_int, c_double, c_char};
use std::ffi::{CString};
use std::path::Path;
use std::{ptr, ptr::null_mut};
use utils::{_string, _last_cpl_err, _last_null_pointer_err};
use raster::{Driver, RasterBand};
use raster::driver::_register_drivers;
use raster::types::GdalType;
use gdal_major_object::MajorObject;
use metadata::Metadata;
use gdal_sys::{self, CPLErr, GDALAccess, GDALDatasetH, GDALDataType, GDALMajorObjectH};

#[cfg(feature = "ndarray")]
use ndarray::Array2;

use errors::*;

pub type GeoTransform = [c_double; 6];

pub struct Dataset {
    c_dataset: GDALDatasetH,
}

impl MajorObject for Dataset {
    unsafe fn gdal_object_ptr(&self) -> GDALMajorObjectH {
        self.c_dataset
    }
}

impl Metadata for Dataset {}

impl Drop for Dataset {
    fn drop(&mut self) {
        unsafe { gdal_sys::GDALClose(self.c_dataset); }
    }
}

impl Dataset {
    pub fn open(path: &Path) -> Result<Dataset> {
        _register_drivers();
        let filename = path.to_string_lossy();
        let c_filename = CString::new(filename.as_ref())?;
        let c_dataset = unsafe { gdal_sys::GDALOpen(c_filename.as_ptr(), GDALAccess::GA_ReadOnly) };
        if c_dataset.is_null() {
            Err(_last_null_pointer_err("GDALOpen"))?;
        }
        Ok(Dataset{c_dataset})
    }

    pub unsafe fn _with_c_ptr(c_dataset: GDALDatasetH) -> Dataset {
        Dataset { c_dataset }
    }

    pub unsafe fn _c_ptr(&self) -> GDALDatasetH {
        self.c_dataset
    }

    pub fn rasterband(&self, band_index: isize) -> Result<RasterBand> {
        unsafe {
            let c_band = gdal_sys::GDALGetRasterBand(self.c_dataset, band_index as c_int);
            if c_band.is_null() {
                Err(_last_null_pointer_err("GDALGetRasterBand"))?;
            }
            Ok(RasterBand::_with_c_ptr(c_band, self))
        }
    }

    pub fn size(&self) -> (usize, usize) {
        let size_x = unsafe { gdal_sys::GDALGetRasterXSize(self.c_dataset) } as usize;
        let size_y = unsafe { gdal_sys::GDALGetRasterYSize(self.c_dataset) } as usize;
        (size_x, size_y)
    }

    /// Get block size from a 'Dataset'.
    /// # Arguments
    /// * band_index - the band_index
    /*
    pub fn size_block(&self, band_index: isize) -> (usize, usize) {
        let band = self.rasterband(band_index)?;
        band.size_block()
    }
    */

    pub fn driver(&self) -> Driver {
        unsafe {
            let c_driver = gdal_sys::GDALGetDatasetDriver(self.c_dataset);
            Driver::_with_c_ptr(c_driver)
        }
    }

    pub fn count(&self) -> isize {
        (unsafe { gdal_sys::GDALGetRasterCount(self.c_dataset) }) as isize
    }

    pub fn projection(&self) -> String {
        let rv = unsafe { gdal_sys::GDALGetProjectionRef(self.c_dataset) };
        _string(rv)
    }

    pub fn set_projection(&self, projection: &str) -> Result<()>{
        let c_projection = CString::new(projection)?;
        unsafe { gdal_sys::GDALSetProjection(self.c_dataset, c_projection.as_ptr()) };
        Ok(())
    }

    /// Affine transformation called geotransformation.
    ///
    /// This is like a linear transformation preserves points, straight lines and planes.
    /// Also, sets of parallel lines remain parallel after an affine transformation.
    /// # Arguments
    /// * transformation - coeficients of transformations
    ///
    /// x-coordinate of the top-left corner pixel (x-offset)
    /// width of a pixel (x-resolution)
    /// row rotation (typically zero)
    /// y-coordinate of the top-left corner pixel
    /// column rotation (typically zero)
    /// height of a pixel (y-resolution, typically negative)
    pub fn set_geo_transform(&self, transformation: &GeoTransform) -> Result<()> {
        assert_eq!(transformation.len(), 6);
        let rv = unsafe {
            gdal_sys::GDALSetGeoTransform(self.c_dataset, transformation.as_ptr() as *mut f64)
        };
        if rv != CPLErr::CE_None {
            Err(_last_cpl_err(rv))?;
        }
        Ok(())
    }

    /// Get affine transformation coefficients.
    ///
    /// x-coordinate of the top-left corner pixel (x-offset)
    /// width of a pixel (x-resolution)
    /// row rotation (typically zero)
    /// y-coordinate of the top-left corner pixel
    /// column rotation (typically zero)
    /// height of a pixel (y-resolution, typically negative)
    pub fn geo_transform(&self) -> Result<GeoTransform> {
        let mut transformation = GeoTransform::default();
        let rv = unsafe {
            gdal_sys::GDALGetGeoTransform(
                self.c_dataset,
                transformation.as_mut_ptr()
            )
        };

        // check if the dataset has a GeoTransform
        if rv != CPLErr::CE_None {
            Err(_last_cpl_err(rv))?;
        }
        Ok(transformation)
    }

    pub fn create_copy(
        &self,
        driver: &Driver,
        filename: &str
    ) -> Result<Dataset> {
        let c_filename = CString::new(filename)?;
        let c_dataset = unsafe { gdal_sys::GDALCreateCopy(
                driver._c_ptr(),
                c_filename.as_ptr(),
                self.c_dataset,
                0,
                null_mut(),
                None,
                null_mut()
            ) };
        if c_dataset.is_null() {
            Err(_last_null_pointer_err("GDALCreateCopy"))?;
        }
        Ok(Dataset{c_dataset})
    }

    pub fn band_type(&self, band_index: isize) -> Result<GDALDataType::Type> {
        self.rasterband(band_index).map(|band| band.band_type())
    }

    /// Read a 'Buffer<u8>' from a 'Dataset'.
    /// # Arguments
    /// * band_index - the band_index
    /// * window - the window position from top left
    /// * window_size - the window size (GDAL will interpolate data if window_size != buffer_size)
    /// * buffer_size - the desired size of the 'Buffer'
    pub fn read_raster(&self,
        band_index: isize,
        window: (isize, isize),
        window_size: (usize, usize),
        size: (usize, usize)
    ) -> Result<ByteBuffer>
    {
        self.read_raster_as::<u8>(
            band_index,
            window,
            window_size,
            size
        )
    }

    /// Read a full 'Dataset' as 'Buffer<T>'.
    /// # Arguments
    /// * band_index - the band_index
    pub fn read_full_raster_as<T: Copy + GdalType>(
        &self,
        band_index: isize,
    ) -> Result<Buffer<T>>
    {
        self.rasterband(band_index)?.read_band_as()
    }

    /// Read a 'Buffer<T>' from a 'Dataset'. T implements 'GdalType'
    /// # Arguments
    /// * band_index - the band_index
    /// * window - the window position from top left
    /// * window_size - the window size (GDAL will interpolate data if window_size != buffer_size)
    /// * buffer_size - the desired size of the 'Buffer'
    pub fn read_raster_as<T: Copy + GdalType>(
        &self,
        band_index: isize,
        window: (isize, isize),
        window_size: (usize, usize),
        size: (usize, usize),
    ) -> Result<Buffer<T>>
    {
        self.rasterband(band_index)?.read_as(window, window_size, size)
    }

    #[cfg(feature = "ndarray")]
    /// Read a 'Array2<T>' from a 'Dataset'. T implements 'GdalType'.
    /// # Arguments
    /// * band_index - the band_index
    /// * window - the window position from top left
    /// * window_size - the window size (GDAL will interpolate data if window_size != array_size)
    /// * array_size - the desired size of the 'Array'
    pub fn read_as_array<T: Copy + GdalType>(
        &self,
        band_index: isize,
        window: (isize, isize),
        window_size: (usize, usize),
        array_size: (usize, usize),
    ) -> Result<Array2<T>>
    {
        self.rasterband(band_index)?.read_as_array(window, window_size, array_size)
    }

    /// Write a 'Buffer<T>' into a 'Dataset'.
    /// # Arguments
    /// * band_index - the band_index
    /// * window - the window position from top left
    /// * window_size - the window size (GDAL will interpolate data if window_size != Buffer.size)
    pub fn write_raster<T: GdalType+Copy>(
        &self,
        band_index: isize,
        window: (isize, isize),
        window_size: (usize, usize),
        buffer: &Buffer<T>
    ) -> Result<()> {
        self.rasterband(band_index)?.write(window, window_size, buffer)
    }

    /// Saves the dataset as a jpeg
    pub fn save_as_png(&self, destination: String, bands: Vec<u8>) -> Result<bool> {
        unsafe {
            println!("destination {:?}", destination);
            let dest_cstring = CString::new(destination).unwrap();
            let dest_c = dest_cstring.as_ptr();
            
            // Select the output format. Starting with GDAL 2.3, if not specified, the format is guessed from the extension (previously was GTiff). Use the short format name.
            let mut args = vec![String::from("-of"), String::from("JPEG")];
            
            // Rescale the input pixels values from the range src_min to src_max to the range dst_min to dst_max. If omitted the output range is 0 to 255. If omitted the input range is automatically computed from the source data. 
            args.push("-scale".to_string());

            // Override the color interpretation of all specified bands. For example -colorinterp red,green,blue,alpha for a 4 band output dataset.
            args.push("-colorinterp".to_string());
            args.push("red,green,blue".to_string());
            
            // Force the output image bands to have a specific data type supported by the driver, which may be one of the following: 
            // Byte, UInt16, Int16, UInt32, Int32, Float32, Float64, CInt16, CInt32, CFloat32 or CFloat64
            args.push("-ot".to_string());
            args.push("Byte".to_string());
            
            // Select an input band band for output. Bands are numbered from 1. Multiple -b switches may be used to select a set of input bands to write to the output file, or to reorder bands. band can also be set to “mask,1” (or just “mask”) to mean the mask band of the first band of the input dataset.
            for band in bands {
                args.push("-b".to_string());
                args.push(band.to_string());
            }

            let mut args_cstring = args
                .iter()
                .map(|x| CString::new(x.clone()).unwrap().into_raw())
                .collect::<Vec<*mut c_char>>();
            args_cstring.push(ptr::null::<i8>() as *mut i8);

            let args_c = args_cstring.as_mut_ptr();
            let options = gdal_sys::GDALTranslateOptionsNew(
                args_c,
                ptr::null::<gdal_sys::GDALTranslateOptionsForBinary>()
                    as *mut gdal_sys::GDALTranslateOptionsForBinary,
            );
            let mut usage_error = 0;
            println!("dest_cstring {:?}", dest_cstring);
            gdal_sys::GDALTranslate(dest_c, self.c_dataset, options, &mut usage_error);
            gdal_sys::GDALTranslateOptionsFree(options);
        }
        Ok(true)
    }

}

pub struct Buffer<T: GdalType> {
    pub size: (usize, usize),
    pub data: Vec<T>,
}

impl<T: GdalType> Buffer<T> {
    pub fn new(size: (usize, usize), data: Vec<T>) -> Buffer<T> {
        Buffer{size, data}
    }
}

pub type ByteBuffer = Buffer<u8>;
