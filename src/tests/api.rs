use super::decimal_year_to_date;
use crate::error::Error::{
    DateOutsideOfValidityRange, HeightOutsideOfValidityRange, InvalidLatitude, InvalidLongitude,
};
use crate::{GeomagneticField, WarningZone, wrap_to_180};
use rstest::rstest;
use time::Date;
use time::Month::{December, January};
use uom::si::angle::{degree, radian};
use uom::si::f32::{Angle, Length, MagneticFluxDensity};
use uom::si::length::{kilometer, meter};
use uom::si::magnetic_flux_density::nanotesla;

#[rstest]
#[case(2023, 15, 0.0069974475, 0.21, 145.0, 128.0, 131.0, 94.0, 157.0)]
#[case(2028, 15, 0.0068380823, 0.20, 138.0, 133.0, 137.0, 89.0, 141.0)]
fn test_uncertainty(
    #[case] year: i32,
    #[case] day: u16,
    #[case] declination_uncertainty_rad: f32,
    #[case] inclination_uncertainty_deg: f32,
    #[case] f_uncertainty_nt: f32,
    #[case] h_uncertainty_nt: f32,
    #[case] x_uncertainty_nt: f32,
    #[case] y_uncertainty_nt: f32,
    #[case] z_uncertainty_nt: f32,
) {
    let result = GeomagneticField::new(
        Length::new::<meter>(100.0),
        Angle::new::<degree>(54.0),
        Angle::new::<degree>(-13.0),
        Date::from_ordinal_date(year, day).unwrap(),
    )
    .unwrap();
    assert_eq!(
        result.declination_uncertainty(),
        Angle::new::<radian>(declination_uncertainty_rad)
    );
    assert_eq!(
        result.inclination_uncertainty(),
        Angle::new::<degree>(inclination_uncertainty_deg)
    );
    assert_eq!(
        result.f_uncertainty(),
        MagneticFluxDensity::new::<nanotesla>(f_uncertainty_nt)
    );
    assert_eq!(
        result.h_uncertainty(),
        MagneticFluxDensity::new::<nanotesla>(h_uncertainty_nt)
    );
    assert_eq!(
        result.x_uncertainty(),
        MagneticFluxDensity::new::<nanotesla>(x_uncertainty_nt)
    );
    assert_eq!(
        result.y_uncertainty(),
        MagneticFluxDensity::new::<nanotesla>(y_uncertainty_nt)
    );
    assert_eq!(
        result.z_uncertainty(),
        MagneticFluxDensity::new::<nanotesla>(z_uncertainty_nt)
    );
}

#[rstest]
#[case(-90.0, 0.42352515)]
#[case(-89.999, 0.423519)]
#[case(-89.99, 0.42346412)]
#[case(-89.9, 0.42293993)]
#[case(-89.0, 0.41797695)]
#[case(-45.0, 0.5829367)]
#[case(0.0, 0.33038762)]
#[case(45.0, 0.35611528)]
#[case(89.0, 2.4856236)]
#[case(89.9, 3.1052732)]
#[case(89.99, 3.1886303)]
#[case(89.999, 3.193964)]
#[case(90.0, 3.194559)]
fn test_declination_uncertainty(#[case] latitude: f32, #[case] expected_uncertainty: f32) {
    let result = GeomagneticField::new(
        Length::new::<meter>(100.0),
        Angle::new::<degree>(latitude),
        Angle::new::<degree>(0.0),
        Date::from_ordinal_date(2023, 15).unwrap(),
    )
    .unwrap();
    assert_eq!(
        result.declination_uncertainty().get::<degree>(),
        expected_uncertainty
    );
}

/// WMM2025.0 magnetic dip poles ([WMM2025 TR](https://doi.org/10.25923/prbc-s316) Table 4): the field
/// is vertical, so \(H \to 0\), \(|I| \to 90^\circ\), \(\sigma_D\) is clamped
/// to \(180^\circ\), and the point is inside the Blackout Zone.
#[rstest]
#[case(85.762, 139.298, 90.0)]
#[case(-63.851, 135.078, -90.0)]
fn test_wmm2025_magnetic_dip_poles(
    #[case] latitude: f32,
    #[case] longitude: f32,
    #[case] expected_inclination: f32,
) {
    use assert_float_eq::assert_float_absolute_eq;

    let result = GeomagneticField::new(
        Length::new::<meter>(0.0),
        Angle::new::<degree>(latitude),
        Angle::new::<degree>(longitude),
        Date::from_ordinal_date(2025, 1).unwrap(),
    )
    .unwrap();

    assert_eq!(result.declination_uncertainty().get::<degree>(), 180.0);
    assert_eq!(
        result.declination_warning(),
        Some(WarningZone::BlackoutZone)
    );
    assert!(result.h().get::<nanotesla>() < 50.0);
    assert_float_absolute_eq!(
        result.inclination().get::<degree>(),
        expected_inclination,
        0.1
    );
}

#[rstest]
#[case(-90.0, None)]
#[case(-89.9, None)]
#[case(-75.0, None)]
#[case(-70.0, Some(WarningZone::CautionZone))]
#[case(-65.0, Some(WarningZone::BlackoutZone))]
#[case(-60.0, Some(WarningZone::CautionZone))]
#[case(-55.0, None)]
#[case(0.0, None)]
#[case(70.0, None)]
#[case(80.0, Some(WarningZone::CautionZone))]
#[case(85.0, Some(WarningZone::BlackoutZone))]
#[case(90.0, Some(WarningZone::BlackoutZone))]
fn test_declination_warning(#[case] latitude: f32, #[case] warning: Option<WarningZone>) {
    let result = GeomagneticField::new(
        Length::new::<meter>(100.0),
        Angle::new::<degree>(latitude),
        Angle::new::<degree>(135.0),
        Date::from_ordinal_date(2023, 15).unwrap(),
    )
    .unwrap();
    assert_eq!(result.declination_warning(), warning);
}

/// Official UPS GV tables (longitude 240° stored as −120°):
/// [WMM2020 TR](https://doi.org/10.25923/ytk1-yx35) Table 6 and WMM2025 TR Table 6.
/// Remaining Table 6 columns are [`test_wmm_tr_table_6_elements`].
#[rstest]
#[case(2020.0, 0.0, 80.0, 0.0, Some(-1.28))]
#[case(2020.0, 0.0, 0.0, 120.0, None)]
#[case(2020.0, 0.0, -80.0, -120.0, Some(-50.64))]
#[case(2020.0, 100.0, 80.0, 0.0, Some(-1.70))]
#[case(2020.0, 100.0, 0.0, 120.0, None)]
#[case(2020.0, 100.0, -80.0, -120.0, Some(-51.22))]
#[case(2022.5, 0.0, 80.0, 0.0, Some(0.01))]
#[case(2022.5, 0.0, 0.0, 120.0, None)]
#[case(2022.5, 0.0, -80.0, -120.0, Some(-50.87))]
#[case(2022.5, 100.0, 80.0, 0.0, Some(-0.41))]
#[case(2022.5, 100.0, 0.0, 120.0, None)]
#[case(2022.5, 100.0, -80.0, -120.0, Some(-51.47))]
#[case(2025.0, 0.0, 80.0, 0.0, Some(1.28))]
#[case(2025.0, 0.0, 0.0, 120.0, None)]
#[case(2025.0, 0.0, -80.0, -120.0, Some(-51.22))]
#[case(2025.0, 100.0, 80.0, 0.0, Some(0.85))]
#[case(2025.0, 100.0, 0.0, 120.0, None)]
#[case(2025.0, 100.0, -80.0, -120.0, Some(-51.79))]
#[case(2027.5, 0.0, 80.0, 0.0, Some(2.59))]
#[case(2027.5, 0.0, 0.0, 120.0, None)]
#[case(2027.5, 0.0, -80.0, -120.0, Some(-51.51))]
#[case(2027.5, 100.0, 80.0, 0.0, Some(2.16))]
#[case(2027.5, 100.0, 0.0, 120.0, None)]
#[case(2027.5, 100.0, -80.0, -120.0, Some(-52.07))]
fn test_wmm_tr_table_6_grid_variation(
    #[case] decimal_year: f32,
    #[case] height_km: f32,
    #[case] lat: f32,
    #[case] lon: f32,
    #[case] expected_gv: Option<f32>,
) {
    use assert_float_eq::assert_float_absolute_eq;

    let result = GeomagneticField::new(
        Length::new::<kilometer>(height_km),
        Angle::new::<degree>(lat),
        Angle::new::<degree>(lon),
        decimal_year_to_date(decimal_year),
    )
    .unwrap();

    match (expected_gv, result.grid_variation()) {
        (None, None) => {}
        (Some(expected), Some(actual)) => {
            assert_float_absolute_eq!(expected, actual.get::<degree>(), 0.03);
        }
        (expected, actual) => {
            panic!("GV mismatch: expected {expected:?}, got {actual:?}");
        }
    }
}

/// Official main-field and secular-variation columns (GV is
/// [`test_wmm_tr_table_6_grid_variation`]):
/// WMM2020 TR Table 6 and WMM2025 TR Table 6.
/// Longitude 240° is stored as −120°.
///
/// TR notes that `f32` synthesis can differ by up to \(0.1\,\mathrm{nT}\);
/// Table 6 SV is printed to \(0.1\,\mathrm{nT/yr}\) and \(0.01^\circ/\mathrm{yr}\).
#[rstest]
#[case(2020.0, 0.0, 80.0, 0.0, [6570.4, -146.3, 54606.0, 6572.0, 55000.1], [83.14, -1.28], [-16.2, 59.0, 42.9, -17.5, 40.5], [0.02, 0.51])]
#[case(2020.0, 0.0, 0.0, 120.0, [39624.3, 109.9, -10932.5, 39624.4, 41104.9], [-15.42, 0.16], [24.2, -60.8, 49.2, 24.0, 10.1], [0.08, -0.09])]
#[case(2020.0, 0.0, -80.0, -120.0, [5940.6, 15772.1, -52480.8, 16853.8, 55120.6], [-72.20, 69.36], [30.4, 1.8, 91.7, 12.4, -83.5], [0.04, -0.10])]
#[case(2020.0, 100.0, 80.0, 0.0, [6261.8, -185.5, 52429.1, 6264.5, 52802.0], [83.19, -1.70], [-15.1, 56.4, 39.2, -16.8, 36.9], [0.02, 0.51])]
#[case(2020.0, 100.0, 0.0, 120.0, [37636.7, 104.9, -10474.8, 37636.9, 39067.3], [-15.55, 0.16], [22.9, -56.1, 45.1, 22.8, 9.8], [0.07, -0.09])]
#[case(2020.0, 100.0, -80.0, -120.0, [5744.9, 14799.5, -49969.4, 15875.4, 52430.6], [-72.37, 68.78], [28.0, 1.4, 85.6, 11.4, -78.1], [0.04, -0.09])]
#[case(2022.5, 0.0, 80.0, 0.0, [6529.9, 1.1, 54713.4, 6529.9, 55101.7], [83.19, 0.01], [-16.2, 59.0, 42.9, -16.2, 40.7], [0.02, 0.52])]
#[case(2022.5, 0.0, 0.0, 120.0, [39684.7, -42.2, -10809.5, 39684.7, 41130.5], [-15.24, -0.06], [24.2, -60.8, 49.2, 24.2, 10.5], [0.08, -0.09])]
#[case(2022.5, 0.0, -80.0, -120.0, [6016.5, 15776.7, -52251.6, 16885.0, 54912.1], [-72.09, 69.13], [30.4, 1.8, 91.7, 12.6, -83.4], [0.04, -0.09])]
#[case(2022.5, 100.0, 80.0, 0.0, [6224.0, -44.5, 52527.0, 6224.2, 52894.5], [83.24, -0.41], [-15.1, 56.4, 39.2, -15.5, 37.1], [0.02, 0.52])]
#[case(2022.5, 100.0, 0.0, 120.0, [37694.0, -35.3, -10362.0, 37694.1, 39092.4], [-15.37, -0.05], [22.9, -56.1, 45.1, 23.0, 10.2], [0.07, -0.09])]
#[case(2022.5, 100.0, -80.0, -120.0, [5815.0, 14803.0, -49755.3, 15904.1, 52235.4], [-72.27, 68.55], [28.0, 1.4, 85.6, 11.6, -78.0], [0.04, -0.09])]
#[case(2025.0, 0.0, 80.0, 0.0, [6521.6, 145.9, 54791.5, 6523.2, 55178.5], [83.21, 1.28], [-8.3, 59.5, 31.1, -7.0, 30.1], [0.01, 0.52])]
#[case(2025.0, 0.0, 0.0, 120.0, [39677.8, -109.6, -10580.2, 39677.9, 41064.3], [-14.93, -0.16], [9.5, -23.1, 79.4, 9.6, -11.2], [0.11, -0.03])]
#[case(2025.0, 0.0, -80.0, -120.0, [6117.5, 15751.9, -52022.5, 16898.1, 54698.2], [-72.00, 68.78], [33.3, -8.6, 95.5, 4.0, -89.6], [0.03, -0.12])]
#[case(2025.0, 100.0, 80.0, 0.0, [6216.0, 92.4, 52598.8, 6216.7, 52964.9], [83.26, 0.85], [-7.7, 56.5, 28.7, -6.9, 27.6], [0.01, 0.52])]
#[case(2025.0, 100.0, 0.0, 120.0, [37688.6, -96.2, -10152.1, 37688.7, 39032.1], [-15.08, -0.15], [9.2, -21.0, 72.9, 9.2, -10.0], [0.11, -0.03])]
#[case(2025.0, 100.0, -80.0, -120.0, [5907.6, 14780.3, -49540.7, 15917.1, 52035.0], [-72.19, 68.21], [30.6, -8.0, 89.2, 3.9, -83.8], [0.03, -0.11])]
#[case(2027.5, 0.0, 80.0, 0.0, [6500.8, 294.5, 54869.4, 6507.5, 55253.9], [83.24, 2.59], [-8.3, 59.5, 31.1, -5.6, 30.3], [0.01, 0.53])]
#[case(2027.5, 0.0, 0.0, 120.0, [39701.6, -167.4, -10381.8, 39702.0, 41036.9], [-14.65, -0.24], [9.5, -23.1, 79.4, 9.6, -10.7], [0.11, -0.03])]
#[case(2027.5, 0.0, -80.0, -120.0, [6200.7, 15730.3, -51783.7, 16908.3, 54474.2], [-71.92, 68.49], [33.3, -8.6, 95.5, 4.2, -89.5], [0.04, -0.12])]
#[case(2027.5, 100.0, 80.0, 0.0, [6196.7, 233.8, 52670.5, 6201.1, 53034.3], [83.29, 2.16], [-7.7, 56.5, 28.7, -5.6, 27.8], [0.01, 0.52])]
#[case(2027.5, 100.0, 0.0, 120.0, [37711.5, -148.7, -9969.8, 37711.8, 39007.4], [-14.81, -0.23], [9.2, -21.0, 72.9, 9.3, -9.7], [0.11, -0.03])]
#[case(2027.5, 100.0, -80.0, -120.0, [5984.0, 14760.1, -49317.7, 15927.0, 51825.7], [-72.10, 67.93], [30.6, -8.0, 89.2, 4.0, -83.7], [0.03, -0.11])]
fn test_wmm_tr_table_6_elements(
    #[case] decimal_year: f32,
    #[case] height_km: f32,
    #[case] lat: f32,
    #[case] lon: f32,
    #[case] xyzhf: [f32; 5],
    #[case] i_d: [f32; 2],
    #[case] sv_nt: [f32; 5],
    #[case] sv_deg: [f32; 2],
) {
    use assert_float_eq::assert_float_absolute_eq;

    let result = GeomagneticField::new(
        Length::new::<kilometer>(height_km),
        Angle::new::<degree>(lat),
        Angle::new::<degree>(lon),
        decimal_year_to_date(decimal_year),
    )
    .unwrap();

    let [x, y, z, h, f] = xyzhf;
    let [inclination, declination] = i_d;
    let [xdot, ydot, zdot, hdot, fdot] = sv_nt;
    let [idot, ddot] = sv_deg;

    const NT: f32 = 0.2;
    const DEG: f32 = 0.03;
    const NT_YR: f32 = 0.15;
    const DEG_YR: f32 = 0.02;

    assert_float_absolute_eq!(x, result.x().get::<nanotesla>(), NT);
    assert_float_absolute_eq!(y, result.y().get::<nanotesla>(), NT);
    assert_float_absolute_eq!(z, result.z().get::<nanotesla>(), NT);
    assert_float_absolute_eq!(h, result.h().get::<nanotesla>(), NT);
    assert_float_absolute_eq!(f, result.f().get::<nanotesla>(), NT);
    assert_float_absolute_eq!(inclination, result.inclination().get::<degree>(), DEG);
    assert_float_absolute_eq!(declination, result.declination().get::<degree>(), DEG);
    assert_float_absolute_eq!(xdot, result.x_rate(), NT_YR);
    assert_float_absolute_eq!(ydot, result.y_rate(), NT_YR);
    assert_float_absolute_eq!(zdot, result.z_rate(), NT_YR);
    assert_float_absolute_eq!(hdot, result.h_rate(), NT_YR);
    assert_float_absolute_eq!(fdot, result.f_rate(), NT_YR);
    assert_float_absolute_eq!(idot, result.inclination_rate(), DEG_YR);
    assert_float_absolute_eq!(ddot, result.declination_rate(), DEG_YR);
}

#[rstest]
#[case(-90.0, true)]
#[case(-54.999, false)]
#[case(-55.0, true)]
#[case(54.999, false)]
#[case(55.0, true)]
#[case(90.0, true)]
fn test_grid_variation_polar_threshold(#[case] latitude: f32, #[case] defined: bool) {
    let result = GeomagneticField::new(
        Length::new::<meter>(0.0),
        Angle::new::<degree>(latitude),
        Angle::new::<degree>(10.0),
        Date::from_ordinal_date(2025, 1).unwrap(),
    )
    .unwrap();
    assert_eq!(result.grid_variation().is_some(), defined);
    if let Some(gv) = result.grid_variation() {
        let d = result.declination().get::<degree>();
        let expected = if latitude >= 0.0 { d - 10.0 } else { d + 10.0 };
        assert_eq!(gv.get::<degree>(), wrap_to_180(expected));
    }
}

#[rstest]
#[case(-91.0)]
#[case(-90.1)]
#[case(-90.00001)]
#[case(90.00001)]
#[case(90.1)]
#[case(91.0)]
#[case(f32::NAN)]
#[case(f32::INFINITY)]
#[case(f32::NEG_INFINITY)]
fn test_invalid_latitude(#[case] latitude: f32) {
    let result = GeomagneticField::new(
        Length::new::<meter>(100.0),
        Angle::new::<degree>(latitude),
        Angle::new::<degree>(0.0),
        Date::from_ordinal_date(2023, 15).unwrap(),
    );
    assert!(result.is_err_and(|e| e == InvalidLatitude));
}

#[rstest]
#[case(-181.0)]
#[case(-180.1)]
#[case(-180.00001)]
#[case(180.00001)]
#[case(180.1)]
#[case(181.0)]
#[case(f32::NAN)]
#[case(f32::INFINITY)]
#[case(f32::NEG_INFINITY)]
fn test_invalid_longitude(#[case] longitude: f32) {
    let result = GeomagneticField::new(
        Length::new::<meter>(100.0),
        Angle::new::<degree>(0.0),
        Angle::new::<degree>(longitude),
        Date::from_ordinal_date(2023, 15).unwrap(),
    );
    assert!(result.is_err_and(|e| e == InvalidLongitude));
}

#[rstest]
#[case(-2.0)]
#[case(-1.1)]
#[case(-1.001)]
#[case(850.001)]
#[case(850.1)]
#[case(851.0)]
#[case(f32::NAN)]
#[case(f32::INFINITY)]
#[case(f32::NEG_INFINITY)]
fn test_height_outside_of_validity_range(#[case] height: f32) {
    let result = GeomagneticField::new(
        Length::new::<kilometer>(height),
        Angle::new::<degree>(0.0),
        Angle::new::<degree>(0.0),
        Date::from_ordinal_date(2023, 15).unwrap(),
    );
    assert!(result.is_err_and(|e| e
        == HeightOutsideOfValidityRange {
            min_height_km: -1.0,
            max_height_km: 850.0
        }));
}

#[rstest]
#[case(-1.0, -90.0, -180.0)]
#[case(850.0, 90.0, 180.0)]
fn test_validity_range(#[case] height: f32, #[case] latitude: f32, #[case] longitude: f32) {
    let result = GeomagneticField::new(
        Length::new::<kilometer>(height),
        Angle::new::<degree>(latitude),
        Angle::new::<degree>(longitude),
        Date::from_ordinal_date(2023, 15).unwrap(),
    );
    assert!(result.is_ok());
    let field = result.unwrap();
    assert!(field.x().get::<nanotesla>().is_finite());
    assert!(field.y().get::<nanotesla>().is_finite());
    assert!(field.z().get::<nanotesla>().is_finite());
}

/// The antimeridian is one meridian: \(+180^\circ\) and \(-180^\circ\) are
/// both legal and must produce the same field. Grid variation is compared
/// as an angle so the closed dateline pair \(\pm 180^\circ\) stays equivalent.
/// Component slack is \(0.01\,\mathrm{nT}\): `f32` \(\sin(\pm\pi)\) is a
/// few ULP from zero and opposite in sign, not a geographic split.
#[rstest]
#[case(0.0)]
#[case(45.0)]
#[case(-45.0)]
#[case(55.0)]
#[case(-55.0)]
#[case(80.0)]
#[case(-80.0)]
#[case(90.0)]
#[case(-90.0)]
fn test_antimeridian_longitude_equivalence(#[case] latitude: f32) {
    use assert_float_eq::assert_float_absolute_eq;

    let date = Date::from_ordinal_date(2025, 1).unwrap();
    let height = Length::new::<meter>(0.0);
    let plus = GeomagneticField::new(
        height,
        Angle::new::<degree>(latitude),
        Angle::new::<degree>(180.0),
        date,
    )
    .unwrap();
    let minus = GeomagneticField::new(
        height,
        Angle::new::<degree>(latitude),
        Angle::new::<degree>(-180.0),
        date,
    )
    .unwrap();

    const NT: f32 = 0.01;
    assert_float_absolute_eq!(
        plus.x().get::<nanotesla>(),
        minus.x().get::<nanotesla>(),
        NT
    );
    assert_float_absolute_eq!(
        plus.y().get::<nanotesla>(),
        minus.y().get::<nanotesla>(),
        NT
    );
    assert_float_absolute_eq!(
        plus.z().get::<nanotesla>(),
        minus.z().get::<nanotesla>(),
        NT
    );
    assert_float_absolute_eq!(
        plus.h().get::<nanotesla>(),
        minus.h().get::<nanotesla>(),
        NT
    );
    assert_float_absolute_eq!(
        plus.f().get::<nanotesla>(),
        minus.f().get::<nanotesla>(),
        NT
    );
    assert_float_absolute_eq!(
        wrap_to_180(plus.declination().get::<degree>() - minus.declination().get::<degree>()),
        0.0,
        1e-4
    );
    assert_float_absolute_eq!(
        plus.inclination().get::<degree>(),
        minus.inclination().get::<degree>(),
        1e-4
    );
    assert_eq!(plus.declination_warning(), minus.declination_warning());

    match (plus.grid_variation(), minus.grid_variation()) {
        (None, None) => {}
        (Some(a), Some(b)) => {
            assert_float_absolute_eq!(
                wrap_to_180(a.get::<degree>() - b.get::<degree>()),
                0.0,
                1e-4
            );
        }
        (a, b) => panic!("grid_variation defined on one side only: {a:?} vs {b:?}"),
    }
}

#[rstest]
#[case(90.0)]
#[case(-90.0)]
fn test_grid_variation_at_geographic_poles(#[case] latitude: f32) {
    use assert_float_eq::assert_float_absolute_eq;

    let date = Date::from_ordinal_date(2025, 1).unwrap();
    for lon in [0.0f32, 90.0, 135.0, 180.0, -120.0] {
        let field = GeomagneticField::new(
            Length::new::<meter>(0.0),
            Angle::new::<degree>(latitude),
            Angle::new::<degree>(lon),
            date,
        )
        .unwrap();
        let gv = field
            .grid_variation()
            .expect("GV is defined at the geographic poles");
        let d = field.declination().get::<degree>();
        let expected = if latitude >= 0.0 { d - lon } else { d + lon };
        assert_float_absolute_eq!(gv.get::<degree>(), wrap_to_180(expected), 1e-5);
    }
}

#[rstest]
#[case(90.0)]
#[case(-90.0)]
fn test_pole_y_finite_and_nonzero_for_some_lon(#[case] latitude: f32) {
    let date = Date::from_ordinal_date(2025, 15).unwrap();
    let mut any_nonzero_y = false;
    for lon in [0.0f32, 90.0, 135.0] {
        let field = GeomagneticField::new(
            Length::new::<meter>(100.0),
            Angle::new::<degree>(latitude),
            Angle::new::<degree>(lon),
            date,
        )
        .unwrap();
        let x = field.x().get::<nanotesla>();
        let y = field.y().get::<nanotesla>();
        let z = field.z().get::<nanotesla>();
        assert!(x.is_finite() && y.is_finite() && z.is_finite());
        any_nonzero_y |= y != 0.0;
    }
    assert!(any_nonzero_y);
}

#[rstest]
#[case(90.0, 89.9)]
#[case(-90.0, -89.9)]
fn test_y_and_h_at_pole_match_nearby_latitude(#[case] pole: f32, #[case] near: f32) {
    use assert_float_eq::assert_float_absolute_eq;
    let date = Date::from_ordinal_date(2025, 15).unwrap();
    for lon in [0.0f32, 90.0, 135.0] {
        let at_pole = GeomagneticField::new(
            Length::new::<meter>(100.0),
            Angle::new::<degree>(pole),
            Angle::new::<degree>(lon),
            date,
        )
        .unwrap();
        let beside = GeomagneticField::new(
            Length::new::<meter>(100.0),
            Angle::new::<degree>(near),
            Angle::new::<degree>(lon),
            date,
        )
        .unwrap();
        assert_float_absolute_eq!(
            at_pole.y().get::<nanotesla>(),
            beside.y().get::<nanotesla>(),
            100.0
        );
        assert_float_absolute_eq!(
            at_pole.h().get::<nanotesla>(),
            beside.h().get::<nanotesla>(),
            100.0
        );
    }
}

#[rstest]
#[case(1918, January, 1)]
#[case(2019, December, 31)]
#[case(2030, January, 1)]
#[case(2229, December, 31)]
fn test_date_outside_validity_range(
    #[case] year: i32,
    #[case] month: time::Month,
    #[case] day: u8,
) {
    let result = GeomagneticField::new(
        Length::new::<kilometer>(0.0),
        Angle::new::<degree>(0.0),
        Angle::new::<degree>(0.0),
        Date::from_calendar_date(year, month, day).unwrap(),
    );
    assert!(result.is_err_and(|e| e
        == DateOutsideOfValidityRange {
            min_date: Date::from_calendar_date(2020, January, 1).unwrap(),
            max_date: Date::from_calendar_date(2029, December, 31).unwrap(),
        }));
}

#[rstest]
#[case(2020, January, 1)]
#[case(2024, December, 31)]
#[case(2025, January, 1)]
#[case(2029, December, 31)]
fn test_date_validity_range(#[case] year: i32, #[case] month: time::Month, #[case] day: u8) {
    let result = GeomagneticField::new(
        Length::new::<kilometer>(0.0),
        Angle::new::<degree>(0.0),
        Angle::new::<degree>(0.0),
        Date::from_calendar_date(year, month, day).unwrap(),
    );
    assert!(result.is_ok());
}
