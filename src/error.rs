//! Error definitions.
use thiserror::Error;

/// The error type for [crate::GeomagneticField::new] method.
#[derive(Error, Debug, PartialEq)]
pub enum Error {
    /// Error due to the date argument outside of model validity range.
    #[error("date outside of model validity range")]
    DateOutsideOfValidityRange,

    /// Error due to the height above WGS 84 ellipsoid argument outside of model validity range.
    #[error("height above WGS 84 ellipsoid outside of model validity range ({min_height_km}km to {max_height_km}km)")]
    HeightOutsideOfValidityRange {
        min_height_km: f32,
        max_height_km: f32,
    },

    /// Error due to the invalid latitude argument.
    #[error("invalid latitude")]
    InvalidLatitude,

    /// Error due to the invalid longitude argument.
    #[error("invalid longitude")]
    InvalidLongitude,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_error() {
        fn assert_error<T: core::error::Error>() {}
        assert_error::<Error>();
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
}
