extern crate gdal;

use std::path::Path;
use gdal::raster::{Dataset, RasterBand};
use gdal::metadata::Metadata;

fn main() {

    let path = Path::new("./fixtures/pasture.tiff");
    let dataset = Dataset::open(path).unwrap();
    println!("dataset description: {:?}", dataset.description());

    let rasterband: RasterBand = dataset.rasterband(1).unwrap();
    println!("rasterband description: {:?}", rasterband.description());
    println!("rasterband no_data_value: {:?}", rasterband.no_data_value());
    println!("rasterband type: {:?}", rasterband.band_type());
    println!("rasterband scale: {:?}", rasterband.scale());
    println!("rasterband offset: {:?}", rasterband.offset());
    if let Ok(rv) = rasterband.read_as::<u8>(
        (20, 30),
        (2, 3),
        (2, 3)
    ) {
        println!("{:?}", rv.data);
    }

    dataset.save_as_png("./fixtures/pasture.png".to_string(), vec![4,3,2]).unwrap();
}
