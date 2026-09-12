extern crate std;

use crate::GeomagneticField;
use crate::WarningZone::{BlackoutZone, CautionZone};
use proptest::prelude::*;
use time::Date;
use uom::si::angle::degree;
use uom::si::f32::{Angle, Length};
use uom::si::length::kilometer;
use uom::si::magnetic_flux_density::nanotesla;

fn almost_equal(a: f32, b: f32, rel: f32, abs: f32) -> bool {
    let diff = (a - b).abs();
    diff <= abs || diff <= rel * a.abs().max(b.abs())
}

/// Inclusive WMM validity window; ordinal is 1..=365 or 1..=366 in leap years.
fn valid_wmm_date() -> impl Strategy<Value = Date> {
    (2020i32..=2029).prop_flat_map(|year| {
        (1u16..=time::util::days_in_year(year))
            .prop_map(move |ordinal| Date::from_ordinal_date(year, ordinal).unwrap())
    })
}

proptest! {
    #[test]
    fn valid_query_invariants(
        height_km in -1.0f32..=850.0,
        latitude_deg in -90.0f32..=90.0,
        longitude_deg in -180.0f32..=180.0,
        date in valid_wmm_date(),
    ) {
        let field = GeomagneticField::new(
            Length::new::<kilometer>(height_km),
            Angle::new::<degree>(latitude_deg),
            Angle::new::<degree>(longitude_deg),
            date,
        )
        .unwrap();

        let x = field.x().get::<nanotesla>();
        let y = field.y().get::<nanotesla>();
        let z = field.z().get::<nanotesla>();
        let h = field.h().get::<nanotesla>();
        let f = field.f().get::<nanotesla>();
        let d = field.declination().get::<degree>();
        let i = field.inclination().get::<degree>();
        let d_unc = field.declination_uncertainty().get::<degree>();

        prop_assert!(x.is_finite() && y.is_finite() && z.is_finite());
        prop_assert!(h.is_finite() && f.is_finite());
        prop_assert!(d.is_finite() && i.is_finite() && d_unc.is_finite());
        prop_assert!(h >= 0.0);
        prop_assert!(f >= h || almost_equal(f, h, 1e-4, 1e-2));
        prop_assert!(almost_equal(h * h, x * x + y * y, 1e-4, 1e-2));
        prop_assert!(almost_equal(f * f, x * x + y * y + z * z, 1e-4, 1e-2));
        prop_assert!((-180.0..=180.0).contains(&d));
        prop_assert!((-90.0..=90.0).contains(&i));
        prop_assert!((0.0..=180.0).contains(&d_unc) && d_unc > 0.0);

        prop_assert_eq!(
            field.declination_warning(),
            if h < 2000.0 {
                Some(BlackoutZone)
            } else if h < 6000.0 {
                Some(CautionZone)
            } else {
                None
            }
        );

        match field.grid_variation() {
            Some(gv) => {
                prop_assert!(latitude_deg.abs() >= 55.0);
                let g = gv.get::<degree>();
                prop_assert!(g.is_finite());
                prop_assert!((-180.0..=180.0).contains(&g));
            }
            None => prop_assert!(latitude_deg.abs() < 55.0),
        }
    }
}
