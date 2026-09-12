extern crate std;

use super::decimal_year_to_date;
use crate::GeomagneticField;
use rstest::rstest;
use uom::si::angle::degree;
use uom::si::f32::{Angle, Length};
use uom::si::length::kilometer;
use uom::si::magnetic_flux_density::nanotesla;

macro_rules! wmm_official_tests {
    ($($file:literal),+ $(,)?) => {
        #[rstest]
        $(
            #[cfg_attr(
                not(has_ncei_test_file = $file),
                ignore = "NCEI test values not present"
            )]
            #[case(concat!("ncei.noaa.gov/", $file))]
        )+
        fn wmm_tests(#[case] test_file: &str) {
            run_wmm_official_tests(test_file);
        }
    };
}

wmm_official_tests!("WMM2020_TestValues.txt", "WMM2025_TestValues.txt");

fn run_wmm_official_tests(test_file: &str) {
    use assert_float_eq::{assert_float_absolute_eq, assert_float_relative_eq};
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::Path;

    let test_path = Path::new(env!("CARGO_MANIFEST_DIR")).join(test_file);
    let test_values_file = File::open(&test_path)
        .unwrap_or_else(|error| panic!("Failed to open {test_file}: {error}"));
    let test_values_lines = BufReader::new(test_values_file).lines();

    for line in test_values_lines {
        match line {
            Ok(content) => {
                if content.trim().starts_with("#") {
                    continue;
                }

                let mut elements = content.split_whitespace();
                let date: f32 = elements
                    .next()
                    .expect("Lack of date!")
                    .parse()
                    .expect("Unable to parse date!");
                let height: f32 = elements
                    .next()
                    .expect("Lack of height!")
                    .parse()
                    .expect("Unable to parse height!");
                let lat: f32 = elements
                    .next()
                    .expect("Lack of lat!")
                    .parse()
                    .expect("Unable to parse lat!");
                let lon: f32 = elements
                    .next()
                    .expect("Lack of lon!")
                    .parse()
                    .expect("Unable to parse lon!");
                let declination: f32 = elements
                    .next()
                    .expect("Lack of declination!")
                    .parse()
                    .expect("Unable to parse declination!");
                let inclination: f32 = elements
                    .next()
                    .expect("Lack of inclination!")
                    .parse()
                    .expect("Unable to parse inclination!");
                let h: f32 = elements
                    .next()
                    .expect("Lack of h!")
                    .parse()
                    .expect("Unable to parse h!");
                let x: f32 = elements
                    .next()
                    .expect("Lack of x!")
                    .parse()
                    .expect("Unable to parse x!");
                let y: f32 = elements
                    .next()
                    .expect("Lack of y!")
                    .parse()
                    .expect("Unable to parse y!");
                let z: f32 = elements
                    .next()
                    .expect("Lack of z!")
                    .parse()
                    .expect("Unable to parse z!");
                let f: f32 = elements
                    .next()
                    .expect("Lack of f!")
                    .parse()
                    .expect("Unable to parse f!");
                let dd_dt: f32 = elements
                    .next()
                    .expect("Lack of dd_dt!")
                    .parse()
                    .expect("Unable to parse dd_dt!");
                let di_dt: f32 = elements
                    .next()
                    .expect("Lack of di_dt!")
                    .parse()
                    .expect("Unable to parse di_dt!");
                let dh_dt: f32 = elements
                    .next()
                    .expect("Lack of dh_dt!")
                    .parse()
                    .expect("Unable to parse dh_dt!");
                let dx_dt: f32 = elements
                    .next()
                    .expect("Lack of dx_dt!")
                    .parse()
                    .expect("Unable to parse dx_dt!");
                let dy_dt: f32 = elements
                    .next()
                    .expect("Lack of dy_dt!")
                    .parse()
                    .expect("Unable to parse dy_dt!");
                let dz_dt: f32 = elements
                    .next()
                    .expect("Lack of dz_dt!")
                    .parse()
                    .expect("Unable to parse dz_dt!");
                let df_dt: f32 = elements
                    .next()
                    .expect("Lack of df_dt!")
                    .parse()
                    .expect("Unable to parse df_dt!");

                let result = GeomagneticField::new(
                    Length::new::<kilometer>(height),
                    Angle::new::<degree>(lat),
                    Angle::new::<degree>(lon),
                    decimal_year_to_date(date),
                )
                .unwrap();

                assert_float_relative_eq!(declination, result.declination().get::<degree>(), 0.04);
                assert_float_absolute_eq!(declination, result.declination().get::<degree>(), 0.006);

                assert_float_relative_eq!(
                    inclination,
                    result.inclination().get::<degree>(),
                    0.0009
                );
                assert_float_absolute_eq!(inclination, result.inclination().get::<degree>(), 0.006);

                assert_float_relative_eq!(h, result.h().get::<nanotesla>(), 0.00005);
                assert_float_absolute_eq!(h, result.h().get::<nanotesla>(), 0.2);

                assert_float_relative_eq!(x, result.x().get::<nanotesla>(), 0.001);
                assert_float_absolute_eq!(x, result.x().get::<nanotesla>(), 0.2);

                assert_float_relative_eq!(y, result.y().get::<nanotesla>(), 0.0013);
                assert_float_absolute_eq!(y, result.y().get::<nanotesla>(), 0.2);

                assert_float_relative_eq!(z, result.z().get::<nanotesla>(), 0.00005);
                assert_float_absolute_eq!(z, result.z().get::<nanotesla>(), 0.3);

                assert_float_relative_eq!(f, result.f().get::<nanotesla>(), 0.000005);
                assert_float_absolute_eq!(f, result.f().get::<nanotesla>(), 0.2);

                assert_float_absolute_eq!(dd_dt, result.declination_rate(), 0.055);

                assert_float_absolute_eq!(di_dt, result.inclination_rate(), 0.055);

                assert_float_absolute_eq!(dh_dt, result.h_rate(), 0.055);

                assert_float_relative_eq!(dx_dt, result.x_rate(), 0.05);
                assert_float_absolute_eq!(dx_dt, result.x_rate(), 0.055);

                assert_float_relative_eq!(dy_dt, result.y_rate(), 0.09);
                assert_float_absolute_eq!(dy_dt, result.y_rate(), 0.055);

                assert_float_relative_eq!(dz_dt, result.z_rate(), 0.012);
                assert_float_absolute_eq!(dz_dt, result.z_rate(), 0.055);

                assert_float_relative_eq!(df_dt, result.f_rate(), 0.02);
                assert_float_absolute_eq!(df_dt, result.f_rate(), 0.055);
            }

            Err(error) => {
                panic!("Error reading lines: {}!", error)
            }
        }
    }
}
