//! Errors returned by [`crate::GeomagneticField::new`].
use thiserror::Error;
use time::Date;

/// The error type for [`crate::GeomagneticField::new`].
#[derive(Error, Debug, Copy, Clone, PartialEq)]
pub enum Error {
    /// Date outside the validity range of the bundled WMM epochs.
    #[error("date outside of model validity range ({min_date} to {max_date})")]
    DateOutsideOfValidityRange {
        /// First day of the bundled WMM validity window (inclusive).
        min_date: Date,
        /// Last day of the bundled WMM validity window (inclusive).
        max_date: Date,
    },

    /// Height above the WGS 84 ellipsoid outside the model validity range.
    #[error(
        "height above WGS 84 ellipsoid outside of model validity range ({min_height_km}km to {max_height_km}km)"
    )]
    HeightOutsideOfValidityRange {
        /// Minimum height of validity range in kilometres above WGS 84 ellipsoid.
        min_height_km: f32,
        /// Maximum height of validity range in kilometres above WGS 84 ellipsoid.
        max_height_km: f32,
    },

    /// Latitude outside the valid range \([-90^\circ, 90^\circ]\).
    #[error("invalid latitude (-90° to 90°)")]
    InvalidLatitude,

    /// Longitude outside the valid range \([-180^\circ, 180^\circ]\).
    #[error("invalid longitude (-180° to 180°)")]
    InvalidLongitude,
}

#[cfg(test)]
mod tests {
    extern crate std;

    use super::*;
    use std::string::ToString;
    use time::Month::{December, January};

    #[test]
    fn test_error() {
        fn assert_error<T: core::error::Error>() {}
        assert_error::<Error>();
    }

    #[test]
    fn test_copy() {
        fn assert_copy<T: Copy>() {}
        assert_copy::<Error>();
    }

    #[test]
    fn test_clone() {
        fn assert_clone<T: Clone>() {}
        assert_clone::<Error>();
    }

    #[test]
    fn test_send() {
        fn assert_send<T: Send>() {}
        assert_send::<Error>();
    }

    #[test]
    fn test_sync() {
        fn assert_sync<T: Sync>() {}
        assert_sync::<Error>();
    }

    #[test]
    fn test_display() {
        let date = Error::DateOutsideOfValidityRange {
            min_date: Date::from_calendar_date(2020, January, 1).unwrap(),
            max_date: Date::from_calendar_date(2029, December, 31).unwrap(),
        };
        assert_eq!(
            date.to_string(),
            "date outside of model validity range (2020-01-01 to 2029-12-31)"
        );

        let height = Error::HeightOutsideOfValidityRange {
            min_height_km: -1.0,
            max_height_km: 850.0,
        };
        assert_eq!(
            height.to_string(),
            "height above WGS 84 ellipsoid outside of model validity range (-1km to 850km)"
        );

        assert_eq!(
            Error::InvalidLatitude.to_string(),
            "invalid latitude (-90° to 90°)"
        );
        assert_eq!(
            Error::InvalidLongitude.to_string(),
            "invalid longitude (-180° to 180°)"
        );
    }
}
