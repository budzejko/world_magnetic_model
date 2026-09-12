//! Query-plus-results wire format for [`crate::GeomagneticField`].
use crate::{GeomagneticField, WarningZone};
use serde::{Deserialize, Serialize};
use time::Date;
use uom::si::angle::degree;
use uom::si::f32::{Angle, Length};
use uom::si::length::meter;
use uom::si::magnetic_flux_density::nanotesla;

/// Inputs accepted when deserializing; output fields in a payload are ignored.
#[derive(Deserialize)]
pub(crate) struct GeomagneticFieldInput {
    pub(crate) height_m: f32,
    pub(crate) latitude_deg: f32,
    pub(crate) longitude_deg: f32,
    #[serde(with = "date_iso")]
    pub(crate) date: Date,
}

/// Serialized query plus public getters. Units are encoded in field names.
#[derive(Serialize)]
pub(crate) struct GeomagneticFieldView {
    height_m: f32,
    latitude_deg: f32,
    longitude_deg: f32,
    #[serde(with = "date_iso")]
    date: Date,
    #[serde(rename = "x_nT")]
    x_nt: f32,
    #[serde(rename = "y_nT")]
    y_nt: f32,
    #[serde(rename = "z_nT")]
    z_nt: f32,
    #[serde(rename = "h_nT")]
    h_nt: f32,
    #[serde(rename = "f_nT")]
    f_nt: f32,
    declination_deg: f32,
    inclination_deg: f32,
    declination_uncertainty_deg: f32,
    inclination_uncertainty_deg: f32,
    #[serde(rename = "x_uncertainty_nT")]
    x_uncertainty_nt: f32,
    #[serde(rename = "y_uncertainty_nT")]
    y_uncertainty_nt: f32,
    #[serde(rename = "z_uncertainty_nT")]
    z_uncertainty_nt: f32,
    #[serde(rename = "h_uncertainty_nT")]
    h_uncertainty_nt: f32,
    #[serde(rename = "f_uncertainty_nT")]
    f_uncertainty_nt: f32,
    declination_warning: Option<WarningZone>,
    grid_variation_deg: Option<f32>,
}

impl From<GeomagneticField> for GeomagneticFieldView {
    fn from(field: GeomagneticField) -> Self {
        Self {
            height_m: field.height_m,
            latitude_deg: field.latitude_deg,
            longitude_deg: field.longitude_deg,
            date: field.date,
            x_nt: field.x().get::<nanotesla>(),
            y_nt: field.y().get::<nanotesla>(),
            z_nt: field.z().get::<nanotesla>(),
            h_nt: field.h().get::<nanotesla>(),
            f_nt: field.f().get::<nanotesla>(),
            declination_deg: field.declination().get::<degree>(),
            inclination_deg: field.inclination().get::<degree>(),
            declination_uncertainty_deg: field.declination_uncertainty().get::<degree>(),
            inclination_uncertainty_deg: field.inclination_uncertainty().get::<degree>(),
            x_uncertainty_nt: field.x_uncertainty().get::<nanotesla>(),
            y_uncertainty_nt: field.y_uncertainty().get::<nanotesla>(),
            z_uncertainty_nt: field.z_uncertainty().get::<nanotesla>(),
            h_uncertainty_nt: field.h_uncertainty().get::<nanotesla>(),
            f_uncertainty_nt: field.f_uncertainty().get::<nanotesla>(),
            declination_warning: field.declination_warning(),
            grid_variation_deg: field.grid_variation().map(|gv| gv.get::<degree>()),
        }
    }
}

impl TryFrom<GeomagneticFieldInput> for GeomagneticField {
    type Error = crate::Error;

    fn try_from(input: GeomagneticFieldInput) -> Result<Self, Self::Error> {
        GeomagneticField::new(
            Length::new::<meter>(input.height_m),
            Angle::new::<degree>(input.latitude_deg),
            Angle::new::<degree>(input.longitude_deg),
            input.date,
        )
    }
}

/// `YYYY-MM-DD` without `alloc`.
mod date_iso {
    use core::fmt;
    use serde::de::{self, Visitor};
    use serde::{Deserializer, Serializer};
    use time::{Date, Month};

    pub(super) fn serialize<S: Serializer>(date: &Date, serializer: S) -> Result<S::Ok, S::Error> {
        let year = date.year();
        if !(0..=9999).contains(&year) {
            return Err(serde::ser::Error::custom(
                "year out of range for YYYY-MM-DD",
            ));
        }
        let mut buf = [0u8; 10];
        write_digits(&mut buf[0..4], year as u32);
        buf[4] = b'-';
        write_digits(&mut buf[5..7], u8::from(date.month()) as u32);
        buf[7] = b'-';
        write_digits(&mut buf[8..10], u32::from(date.day()));
        serializer.serialize_str(core::str::from_utf8(&buf).map_err(serde::ser::Error::custom)?)
    }

    pub(super) fn deserialize<'de, D: Deserializer<'de>>(
        deserializer: D,
    ) -> Result<Date, D::Error> {
        deserializer.deserialize_str(DateVisitor)
    }

    fn write_digits(buf: &mut [u8], mut value: u32) {
        for slot in buf.iter_mut().rev() {
            *slot = b'0' + (value % 10) as u8;
            value /= 10;
        }
    }

    fn parse_date(s: &str) -> Result<Date, &'static str> {
        let bytes = s.as_bytes();
        if bytes.len() != 10 || bytes[4] != b'-' || bytes[7] != b'-' {
            return Err("date must be YYYY-MM-DD");
        }
        let year: i32 = s[0..4].parse().map_err(|_| "invalid year")?;
        let month_n: u8 = s[5..7].parse().map_err(|_| "invalid month")?;
        let day: u8 = s[8..10].parse().map_err(|_| "invalid day")?;
        let month = Month::try_from(month_n).map_err(|_| "invalid month")?;
        Date::from_calendar_date(year, month, day).map_err(|_| "invalid date")
    }

    struct DateVisitor;

    impl<'de> Visitor<'de> for DateVisitor {
        type Value = Date;

        fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            f.write_str("date in YYYY-MM-DD format")
        }

        fn visit_str<E: de::Error>(self, value: &str) -> Result<Date, E> {
            parse_date(value).map_err(E::custom)
        }
    }
}

#[cfg(test)]
mod tests {
    extern crate std;

    use super::GeomagneticFieldInput;
    use crate::{GeomagneticField, WarningZone};
    use rstest::rstest;
    use std::string::ToString;
    use time::Date;
    use time::Month::January;
    use uom::si::angle::degree;
    use uom::si::f32::{Angle, Length};
    use uom::si::length::meter;
    use uom::si::magnetic_flux_density::nanotesla;

    fn sample_field(year: i32, month: time::Month, day: u8, latitude_deg: f32) -> GeomagneticField {
        GeomagneticField::new(
            Length::new::<meter>(100.0),
            Angle::new::<degree>(latitude_deg),
            Angle::new::<degree>(-7.91),
            Date::from_calendar_date(year, month, day).unwrap(),
        )
        .unwrap()
    }

    fn assert_public_getters_eq(original: &GeomagneticField, restored: &GeomagneticField) {
        assert_eq!(original.x(), restored.x());
        assert_eq!(original.y(), restored.y());
        assert_eq!(original.z(), restored.z());
        assert_eq!(original.h(), restored.h());
        assert_eq!(original.f(), restored.f());
        assert_eq!(original.declination(), restored.declination());
        assert_eq!(original.inclination(), restored.inclination());
        assert_eq!(
            original.declination_uncertainty(),
            restored.declination_uncertainty()
        );
        assert_eq!(
            original.inclination_uncertainty(),
            restored.inclination_uncertainty()
        );
        assert_eq!(original.x_uncertainty(), restored.x_uncertainty());
        assert_eq!(original.y_uncertainty(), restored.y_uncertainty());
        assert_eq!(original.z_uncertainty(), restored.z_uncertainty());
        assert_eq!(original.h_uncertainty(), restored.h_uncertainty());
        assert_eq!(original.f_uncertainty(), restored.f_uncertainty());
        assert_eq!(
            original.declination_warning(),
            restored.declination_warning()
        );
        assert_eq!(original.grid_variation(), restored.grid_variation());
    }

    fn roundtrip_json(original: &GeomagneticField) -> GeomagneticField {
        let json = serde_json::to_string(original).unwrap();
        serde_json::from_str(&json).unwrap()
    }

    fn roundtrip_postcard(original: &GeomagneticField) -> GeomagneticField {
        let mut buf = [0u8; 256];
        let bytes = postcard::to_slice(original, &mut buf).unwrap();
        postcard::from_bytes(bytes).unwrap()
    }

    #[rstest]
    #[case(2024, January, 1, 37.03)]
    #[case(2029, January, 15, 37.03)]
    #[case(2025, January, 1, 80.0)]
    fn roundtrip_getters(
        #[case] year: i32,
        #[case] month: time::Month,
        #[case] day: u8,
        #[case] latitude_deg: f32,
    ) {
        let original = sample_field(year, month, day, latitude_deg);
        assert_public_getters_eq(&original, &roundtrip_json(&original));
        assert_public_getters_eq(&original, &roundtrip_postcard(&original));
        assert_eq!(
            original.x().get::<nanotesla>(),
            roundtrip_postcard(&original).x().get::<nanotesla>()
        );
    }

    #[test]
    fn json_contains_query_and_public_results() {
        let original = sample_field(2029, January, 15, 37.03);
        let value = serde_json::to_value(&original).unwrap();
        let object = value.as_object().unwrap();
        for key in [
            "height_m",
            "latitude_deg",
            "longitude_deg",
            "date",
            "x_nT",
            "y_nT",
            "z_nT",
            "h_nT",
            "f_nT",
            "declination_deg",
            "inclination_deg",
            "declination_uncertainty_deg",
            "inclination_uncertainty_deg",
            "x_uncertainty_nT",
            "y_uncertainty_nT",
            "z_uncertainty_nT",
            "h_uncertainty_nT",
            "f_uncertainty_nT",
            "declination_warning",
            "grid_variation_deg",
        ] {
            assert!(object.contains_key(key), "missing {key}");
        }
        assert!(!object.contains_key("error_model"));
        assert!(!object.contains_key("epoch_year"));
        assert_eq!(object["date"], "2029-01-15");
        assert_eq!(object["height_m"], 100.0);
        assert_eq!(
            object["x_uncertainty_nT"],
            original.x_uncertainty().get::<nanotesla>()
        );
        assert_eq!(
            object["y_uncertainty_nT"],
            original.y_uncertainty().get::<nanotesla>()
        );
        assert_eq!(
            object["z_uncertainty_nT"],
            original.z_uncertainty().get::<nanotesla>()
        );
        assert_eq!(
            object["h_uncertainty_nT"],
            original.h_uncertainty().get::<nanotesla>()
        );
        assert_eq!(
            object["f_uncertainty_nT"],
            original.f_uncertainty().get::<nanotesla>()
        );
        assert_eq!(
            object["declination_uncertainty_deg"],
            original.declination_uncertainty().get::<degree>()
        );
        assert_eq!(
            object["inclination_uncertainty_deg"],
            original.inclination_uncertainty().get::<degree>()
        );
        assert!(object["declination_warning"].is_null());
        assert!(object["grid_variation_deg"].is_null());
    }

    #[test]
    fn json_includes_grid_variation_in_polar_region() {
        let original = sample_field(2025, January, 1, 80.0);
        let value = serde_json::to_value(&original).unwrap();
        assert!(value["grid_variation_deg"].as_f64().is_some());
    }

    #[test]
    fn deserialize_ignores_output_fields() {
        let json = r#"{
                "height_m": 100.0,
                "latitude_deg": 37.03,
                "longitude_deg": -7.91,
                "date": "2029-01-15",
                "declination_deg": 0.0,
                "x_nT": 0.0
            }"#;
        let restored: GeomagneticField = serde_json::from_str(json).unwrap();
        let expected = sample_field(2029, January, 15, 37.03);
        assert_public_getters_eq(&expected, &restored);
    }

    #[test]
    fn warning_zone_roundtrip() {
        for zone in [WarningZone::BlackoutZone, WarningZone::CautionZone] {
            let json = serde_json::to_string(&zone).unwrap();
            let from_json: WarningZone = serde_json::from_str(&json).unwrap();
            assert_eq!(zone, from_json);

            let mut buf = [0u8; 8];
            let bytes = postcard::to_slice(&zone, &mut buf).unwrap();
            let from_postcard: WarningZone = postcard::from_bytes(bytes).unwrap();
            assert_eq!(zone, from_postcard);
        }
    }

    #[test]
    fn reject_nan_latitude() {
        let error = GeomagneticField::try_from(GeomagneticFieldInput {
            height_m: 100.0,
            latitude_deg: f32::NAN,
            longitude_deg: 0.0,
            date: Date::from_calendar_date(2025, January, 1).unwrap(),
        })
        .unwrap_err();
        assert_eq!(error, crate::Error::InvalidLatitude);
    }

    #[test]
    fn reject_invalid_latitude() {
        let error = serde_json::from_str::<GeomagneticField>(
            r#"{"height_m":100.0,"latitude_deg":91.0,"longitude_deg":0.0,"date":"2025-01-01"}"#,
        )
        .unwrap_err();
        assert!(error.to_string().contains("invalid latitude"));
    }

    #[test]
    fn reject_date_outside_validity_range() {
        let error = serde_json::from_str::<GeomagneticField>(
            r#"{"height_m":100.0,"latitude_deg":0.0,"longitude_deg":0.0,"date":"2015-01-01"}"#,
        )
        .unwrap_err();
        assert!(
            error
                .to_string()
                .contains("date outside of model validity range")
        );
    }
}
