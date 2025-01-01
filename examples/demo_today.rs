use time::OffsetDateTime;
use uom::fmt::DisplayStyle::Abbreviation;
use uom::si::angle::degree;
use uom::si::f32::{Angle, Length};
use uom::si::length::meter;
use world_magnetic_model::GeomagneticField;

fn main() {
    let today = OffsetDateTime::now_utc().date();
    let points = [
        ("Świnoujście", 53.912668, 14.261214),
        ("Jakuszyce", 50.814168, 15.428267),
        ("Sejny", 54.109385, 23.346215),
        ("Wołosate", 49.066563, 22.680121),
    ];

    println!("{}", today);

    for point in points {
        let geomagnetic_field = GeomagneticField::new(
            Length::new::<meter>(100.0),
            Angle::new::<degree>(point.1),
            Angle::new::<degree>(point.2),
            today,
        )
        .unwrap();

        println!(
            "Declination in {:11} {:.2} ±{:.2}",
            point.0,
            geomagnetic_field
                .declination()
                .into_format_args(degree, Abbreviation),
            geomagnetic_field
                .declination_uncertainty()
                .into_format_args(degree, Abbreviation)
        );
    }
}
