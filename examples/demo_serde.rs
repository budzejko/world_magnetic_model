use world_magnetic_model::GeomagneticField;
use world_magnetic_model::time::Date;
use world_magnetic_model::uom::si::angle::degree;
use world_magnetic_model::uom::si::f32::{Angle, Length};
use world_magnetic_model::uom::si::length::meter;

fn main() {
    let geomagnetic_field = GeomagneticField::new(
        Length::new::<meter>(100.0),
        Angle::new::<degree>(37.03),
        Angle::new::<degree>(-7.91),
        Date::from_ordinal_date(2029, 15).unwrap(),
    )
    .unwrap();

    let json = serde_json::to_string_pretty(&geomagnetic_field).unwrap();
    let from_json: GeomagneticField = serde_json::from_str(&json).unwrap();

    let mut buf = [0u8; 256];
    let bytes = postcard::to_slice(&geomagnetic_field, &mut buf).unwrap();
    let from_postcard: GeomagneticField = postcard::from_bytes(bytes).unwrap();

    println!("JSON ({} bytes):\n{json}", json.len());
    println!("postcard ({} bytes):\n{}", bytes.len(), format_hex(bytes));
    println!(
        "declination JSON {}  postcard {}",
        from_json.declination().get::<degree>(),
        from_postcard.declination().get::<degree>()
    );
}

fn format_hex(bytes: &[u8]) -> String {
    bytes
        .chunks(16)
        .map(|chunk| {
            chunk
                .iter()
                .map(|byte| format!("{byte:02x}"))
                .collect::<Vec<_>>()
                .join(" ")
        })
        .collect::<Vec<_>>()
        .join("\n")
}
