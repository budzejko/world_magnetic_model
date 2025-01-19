mod wmm {
    use world_magnetic_model::time::format_description::well_known::Iso8601;
    use world_magnetic_model::time::Date;
    use world_magnetic_model::uom::si::angle::{degree, mil};
    use world_magnetic_model::uom::si::f32::{Angle, Length};
    use world_magnetic_model::uom::si::length::meter;
    use world_magnetic_model::GeomagneticField;

    pub struct Wmm(GeomagneticField);

    impl Wmm {
        pub fn new(height_in_meters: f32, lat: f32, lon: f32, date: &str) -> Wmm {
            let date = Date::parse(date, &Iso8601::DATE).unwrap();

            Wmm(GeomagneticField::new(
                Length::new::<meter>(height_in_meters),
                Angle::new::<degree>(lat),
                Angle::new::<degree>(lon),
                date,
            )
            .unwrap())
        }

        pub fn declination_deg(&self) -> f32 {
            self.0.declination().get::<degree>()
        }

        pub fn declination_mil(&self) -> f32 {
            self.0.declination().get::<mil>()
        }
    }
}

fn main() {
    let date = "2025-01-09";
    let points = [
        ("Świnoujście", 53.912668, 14.261214),
        ("Jakuszyce", 50.814168, 15.428267),
        ("Sejny", 54.109385, 23.346215),
        ("Wołosate", 49.066563, 22.68012),
    ];

    println!("{}", date);

    for point in points {
        let geomagnetic_field = wmm::Wmm::new(100.0, point.1, point.2, date);

        println!(
            "Declination in {:11} {:.2}° {:.1} mil",
            point.0,
            geomagnetic_field.declination_deg(),
            geomagnetic_field.declination_mil()
        );
    }
}
