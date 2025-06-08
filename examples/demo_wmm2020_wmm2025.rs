use world_magnetic_model::time::{Date, Month};
use world_magnetic_model::uom::fmt::DisplayStyle::Abbreviation;
use world_magnetic_model::uom::si::angle::degree;
use world_magnetic_model::uom::si::f32::{Angle, Length};
use world_magnetic_model::uom::si::length::meter;
use world_magnetic_model::GeomagneticField;

fn main() {
    let points = [
        ("Świnoujście", 53.912668, 14.261214),
        ("Jakuszyce", 50.814168, 15.428267),
        ("Sejny", 54.109385, 23.346215),
        ("Wołosate", 49.066563, 22.68012),
    ];

    for point in points {
        for year in 2020..2030 {
            for date in [
                Date::from_ordinal_date(year, 1).unwrap(),
                Date::from_calendar_date(year, Month::December, 31).unwrap(),
            ] {
                let geomagnetic_field = GeomagneticField::new(
                    Length::new::<meter>(100.0),
                    Angle::new::<degree>(point.1),
                    Angle::new::<degree>(point.2),
                    date,
                )
                .unwrap();

                println!(
                    "{} Declination in {:11} {:.2} ±{:.2}",
                    date,
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
        println!();
    }
}
