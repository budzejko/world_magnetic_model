extern crate std;

use super::{decimal_year_to_date, years_since_epoch_to_date};
use crate::{WGS84_A, WGS84_E2, WGS84_F, geodetic_to_geocentric, wrap_to_180, years_since_epoch};
use libm::{cosf, sinf, sqrtf};
use rstest::rstest;
use time::Date;
use time::Month::{December, January, July};

#[test]
fn test_grid_variation_wraps_to_plus_minus_180() {
    assert_eq!(wrap_to_180(181.0), -179.0);
    assert_eq!(wrap_to_180(-181.0), 179.0);
    assert_eq!(wrap_to_180(180.0), 180.0);
    assert_eq!(wrap_to_180(-180.0), -180.0);
    assert_eq!(wrap_to_180(0.0), 0.0);
}

#[rstest]
#[case(2020, July, 2, 2020)]
#[case(2021, July, 3, 2020)]
#[case(2023, January, 1, 2020)]
#[case(2023, January, 2, 2020)]
#[case(2024, January, 2, 2020)]
#[case(2024, December, 12, 2020)]
#[case(2024, December, 31, 2020)]
#[case(2025, January, 1, 2025)]
#[case(2029, December, 31, 2025)]
fn test_years_since_epoch(
    #[case] year_int: i32,
    #[case] month: time::Month,
    #[case] day: u8,
    #[case] epoch_year: i32,
) {
    use assert_float_eq::assert_float_absolute_eq;

    let date = Date::from_calendar_date(year_int, month, day).unwrap();
    let days = time::util::days_in_year(year_int) as f64;
    let expected = (year_int - epoch_year) as f64 + f64::from(date.ordinal() - 1) / days;
    let years_since = years_since_epoch(date, epoch_year);
    assert_float_absolute_eq!(f64::from(years_since), expected, 3e-7);
    assert_eq!(years_since_epoch_to_date(epoch_year, years_since), date);
    assert_eq!(
        years_since_epoch(
            years_since_epoch_to_date(epoch_year, years_since),
            epoch_year
        ),
        years_since
    );
}

#[rstest]
#[case(2020.5, 2020, July, 2)]
#[case(2023.0, 2023, January, 1)]
#[case(2025.0, 2025, January, 1)]
#[case(2025.5, 2025, July, 3)]
fn test_decimal_year_to_date(
    #[case] decimal_year: f32,
    #[case] year_int: i32,
    #[case] month: time::Month,
    #[case] day: u8,
) {
    assert_eq!(
        decimal_year_to_date(decimal_year),
        Date::from_calendar_date(year_int, month, day).unwrap()
    );
}

#[test]
fn test_geodetic_to_geocentric_equator() {
    use assert_float_eq::assert_float_absolute_eq;

    let (r, cos_theta, sin_theta, sin_phi, cos_phi) = geodetic_to_geocentric(0.0, 0.0);
    assert_float_absolute_eq!(r, WGS84_A, 0.1);
    assert_float_absolute_eq!(cos_theta, 0.0, 0.0001);
    assert_float_absolute_eq!(sin_theta, 1.0, 0.0001);
    assert_float_absolute_eq!(sin_phi, 0.0, 0.0001);
    assert_float_absolute_eq!(cos_phi, 1.0, 0.0001);
}

#[test]
fn test_geodetic_to_geocentric_north_pole() {
    use assert_float_eq::assert_float_absolute_eq;
    use std::f32::consts::PI;

    let (r, cos_theta, sin_theta, sin_phi, cos_phi) = geodetic_to_geocentric(0.0, PI / 2.0);
    assert_float_absolute_eq!(r, WGS84_A * (1.0 - WGS84_F), 0.1);
    assert_float_absolute_eq!(cos_theta, 1.0, 0.0001);
    assert_float_absolute_eq!(sin_theta, 0.0, 0.0001);
    assert_float_absolute_eq!(sin_phi, 1.0, 0.0001);
    assert_float_absolute_eq!(cos_phi, 0.0, 0.0001);
}

#[test]
fn test_geodetic_to_geocentric_south_pole() {
    use assert_float_eq::assert_float_absolute_eq;
    use std::f32::consts::PI;

    let (r, cos_theta, sin_theta, sin_phi, cos_phi) = geodetic_to_geocentric(0.0, -PI / 2.0);
    assert_float_absolute_eq!(r, WGS84_A * (1.0 - WGS84_F), 0.1);
    assert_float_absolute_eq!(cos_theta, -1.0, 0.0001);
    assert_float_absolute_eq!(sin_theta, 0.0, 0.0001);
    assert_float_absolute_eq!(sin_phi, -1.0, 0.0001);
    assert_float_absolute_eq!(cos_phi, 0.0, 0.0001);
}

#[test]
fn test_geodetic_to_geocentric_45_degrees_latitude() {
    use assert_float_eq::assert_float_absolute_eq;
    use std::f32::consts::PI;

    let (r, cos_theta, sin_theta, sin_phi, cos_phi) = geodetic_to_geocentric(1000.0, PI / 4.0);

    let r_c = WGS84_A / sqrtf(1.0 - WGS84_E2 * sinf(PI / 4.0) * sinf(PI / 4.0));
    let p = (r_c + 1000.0) * cosf(PI / 4.0);
    let z = (r_c * (1.0 - WGS84_E2) + 1000.0) * sinf(PI / 4.0);
    let expected_r = sqrtf(p * p + z * z);

    assert_float_absolute_eq!(r, expected_r, 0.0001);
    assert_float_absolute_eq!(cos_theta, z / expected_r, 0.0001);
    assert_float_absolute_eq!(sin_theta, p / expected_r, 0.0001);
    assert_float_absolute_eq!(sin_phi, sinf(PI / 4.0), 0.0001);
    assert_float_absolute_eq!(cos_phi, cosf(PI / 4.0), 0.0001);
}

#[test]
fn test_geodetic_to_geocentric_below_sea_level() {
    use assert_float_eq::assert_float_absolute_eq;

    let (r, cos_theta, sin_theta, sin_phi, cos_phi) = geodetic_to_geocentric(-1000.0, 0.0);
    assert_float_absolute_eq!(r, WGS84_A - 1000.0, 0.0001);
    assert_float_absolute_eq!(cos_theta, 0.0, 0.0001);
    assert_float_absolute_eq!(sin_theta, 1.0, 0.0001);
    assert_float_absolute_eq!(sin_phi, 0.0, 0.0001);
    assert_float_absolute_eq!(cos_phi, 1.0, 0.0001);
}
