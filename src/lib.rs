#![no_std]
#![warn(missing_docs)]
//! This crate is a Rust implementation of the NOAA [World Magnetic Model](https://www.ncei.noaa.gov/products/world-magnetic-model),
//! a mathematical model of the magnetic field produced by the Earth's core and its variation over time.
//!
//! Crate interface is using Units of measurement [uom](https://docs.rs/uom/latest/uom/) crate for physical quantities
//! representation. WMM coeficient files are converted to code constants to avoid file reading in runtime. Implemented models:
//! WMM2020 and WMM2025.
//!
//! # Usage
//! ```rust
//! # use core::error::Error;
//! # use world_magnetic_model::time::Date;
//! # use world_magnetic_model::uom::si::f32::{Angle, Length};
//! # use world_magnetic_model::uom::si::angle::degree;
//! # use world_magnetic_model::uom::si::length::meter;
//! # use world_magnetic_model::GeomagneticField;
//! # fn main() -> Result<(), Box<dyn Error>> {
//! let geomagnetic_field = GeomagneticField::new(
//!     Length::new::<meter>(100.0), // height
//!     Angle::new::<degree>(54.20), // lat
//!     Angle::new::<degree>(18.67), // lon
//!     Date::from_ordinal_date(2023, 15)? // date
//! )?;
//!
//! assert_eq!(
//!     geomagnetic_field.declination().get::<degree>(),
//!     6.439543
//! );
//! assert_eq!(
//!     geomagnetic_field.declination_uncertainty().get::<degree>(),
//!     0.41047537
//! );
//! #     Ok(())
//! # }
//! ```
//!
//! # World Magnetic Model
//!
//! The [World Magnetic Model](https://www.ncei.noaa.gov/products/world-magnetic-model) is the standard model used by
//! the U.S. Department of Defense, the U.K. Ministry of Defence, the North Atlantic Treaty Organization (NATO)
//! and the International Hydrographic Organization (IHO), for navigation, attitude and heading referencing systems
//! using the geomagnetic field. It is also used widely in civilian navigation and heading systems.
//! The model is produced at 5-year intervals, with the current model expiring on December 31, 2029. The current
//! model WMM2025 is produced jointly by the NCEI and the British Geological Survey (BGS). The model, associated
//! software, and documentation are distributed by NCEI on behalf of US National Geospatial-Intelligence Agency
//! and by BGS on behalf of UK Defence Geographic Centre.
//! [\[source\]](https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/gov.noaa.ngdc:WMM2025/html)
//!
//! Please refer to [World Magnetic Model Accuracy, Limitations, and Error Model](https://www.ncei.noaa.gov/products/world-magnetic-model/accuracy-limitations-error-model).
//!
//! # NOAA License Statement
//! The WMM source code is in the public domain and not licensed or under copyright. The information and software
//! may be used freely by the public. As required by 17 U.S.C. 403, third parties producing copyrighted works
//! consisting predominantly of the material produced by U.S. government agencies must provide notice with such
//! work(s) identifying the U.S. Government material incorporated and stating that such material is not subject
//! to copyright protection.
//!
//! The WMM model and associated data files are produced by the U.S. Government and are not subject to copyright.
//!
//! # Credit
//! This work was inspired from [geomag-wmm](https://git.hostux.fr/ConstellationVFR/geomag-wmm).
//!
//! # License
//! Licensed under either of [Apache License, Version 2.0](https://github.com/budzejko/world_magnetic_model/blob/main/LICENSE-APACHE)
//! or [MIT license](https://github.com/budzejko/world_magnetic_model/blob/main/LICENSE-MIT) at your option.

pub use error::Error;
pub use time;
pub use uom;

mod error;
mod math;
mod wmm;
mod wmm_models;

use libm::{asinf, atan2f, cos, cosf, sinf, sqrt, sqrtf, tan};
use time::Date;
#[allow(unused)]
use uom::num_traits::float::FloatCore;
use uom::si::angle::{degree, radian};
use uom::si::f32::{Angle, Length, MagneticFluxDensity};
use uom::si::length::{kilometer, meter};
use uom::si::magnetic_flux_density::nanotesla;

use error::Error::{HeightOutsideOfValidityRange, InvalidLatitude, InvalidLongitude};
use math::{index, schmidt_semi_normalised_associated_legendre};
use wmm_models::WmmErrorModel;

/// Represents the geomagnetic field at given point and date.
///
/// # Examples
/// ```rust
/// # use core::error::Error;
/// # use world_magnetic_model::time::Date;
/// # use world_magnetic_model::uom::si::f32::{Angle, Length, MagneticFluxDensity};
/// # use world_magnetic_model::uom::si::angle::{degree, radian};
/// # use world_magnetic_model::uom::si::length::foot;
/// # use world_magnetic_model::uom::si::magnetic_flux_density::{nanotesla, gauss};
/// # use world_magnetic_model::GeomagneticField;
/// # use world_magnetic_model::WarningZone::CautionZone;
/// # use world_magnetic_model::uom::fmt::DisplayStyle::{Abbreviation, Description};
/// # fn main() -> Result<(), Box<dyn Error>> {
/// // various units can be used on input
/// let geomagnetic_field = GeomagneticField::new(
///     Length::new::<foot>(3000.0), // height
///     Angle::new::<degree>(80.0), // lat
///     Angle::new::<radian>(2.36), // lon
///     Date::from_ordinal_date(2023, 15)? // date
/// )?;
///
/// // declination results can be interpreted as degrees or radians
/// assert_eq!(
///     geomagnetic_field.declination().get::<degree>(),
///     -19.997639
/// );
/// assert_eq!(
///     geomagnetic_field.declination().get::<radian>(),
///     -0.34902462
/// );
/// assert_eq!(
///     geomagnetic_field.declination_uncertainty().get::<degree>(),
///     2.088066
/// );
/// assert_eq!(
///     geomagnetic_field.declination_uncertainty().get::<radian>(),
///     0.03644363
/// );
///
/// // warning conditions can be checked e.g. for user interface warning
/// assert_eq!(
///     geomagnetic_field.declination_warning(),
///     Some(CautionZone)
/// );
///
/// // magnetic flux density results can be interpreted as nT or gauss
/// assert_eq!(
///     geomagnetic_field.f().get::<nanotesla>(),
///     58768.96
/// );
/// assert_eq!(
///     geomagnetic_field.f().get::<gauss>(),
///     0.58768964
/// );
///
/// // results can be formatted with unit abbreviation or full description
/// assert_eq!(
///     format!("{:.2}", geomagnetic_field.f().into_format_args(nanotesla, Abbreviation)),
///     "58768.96 nT"
/// );
/// assert_eq!(
///     format!("{:.2}", geomagnetic_field.f().into_format_args(nanotesla, Description)),
///     "58768.96 nanoteslas"
/// );
/// #     Ok(())
/// # }
/// ```
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct GeomagneticField {
    x: f32, // Northern component of the magnetic field vector
    y: f32, // Eastern component of the magnetic field vector
    z: f32, // Downward component of the magnetic field vector
    #[cfg(test)]
    x_rate: f32, // Yearly rate of change in the northern component
    #[cfg(test)]
    y_rate: f32, // Yearly rate of change in the eastern component
    #[cfg(test)]
    z_rate: f32, // Yearly rate of change in the downward component
    error_model: WmmErrorModel, // Values to calculate uncertainty
}

/// Warning zone type, as defined in
/// [World Magnetic Model Accuracy, Limitations, and Error Model](https://www.ncei.noaa.gov/products/world-magnetic-model/accuracy-limitations-error-model).
#[derive(PartialEq, Debug)]
pub enum WarningZone {
    /// Based on the WMM military specification, “Blackout Zones” (BoZ) are areas around
    /// the north and south magnetic poles where compasses are not accurate and should not be
    /// relied on for navigation. The BoZ are defined as regions around the north and south
    /// magnetic poles where the horizontal intensity of Earth’s magnetic field (H) is less
    /// than 2000 nT. In BoZs, WMM declination values are not accurate and compasses are
    /// unreliable.
    BlackoutZone,
    /// “Caution Zone” (2000 nT <= H < 6000 nT) is area around the perimeter of the BoZs where
    /// compasses should be used with caution because they may not be fully accurate.
    CautionZone,
}

const MIN_HEIGHT_KM: f32 = -1.0;
const MAX_HEIGHT_KM: f32 = 850.0;

// WGS84 ellipsoid parameters
const A: f32 = 6378137.0; // Equation (8)
                          // 6378137.0m Semi Major Axis (a) Equator to Center [m]
const F: f32 = 1.0 / 298.25723; // Equation (8) // Flattening (f) = (a-b)/a
                                // 6356752.3142m Semi Minor Axis (b) Pole to Center
const E_2: f32 = F * (2.0 - F); // Equation (8)
                                // Geomagnetic parameters
const A2: f32 = 6371200.0; // Equation (4) // Geomagnetic Reference Radius

impl GeomagneticField {
    /// Calculates the geomagnetic field at a given location and date.
    ///
    /// # Arguments
    ///
    /// * `height` - the height above World Geodetic System 1984 (WGS 84) ellipsoid (from -1km to +850km).
    /// * `latitude` - the WGS 84 latitude.
    /// * `longitude` - the WGS 84 longitude.
    /// * `date` - the date for the WMM model evaluation.
    ///
    /// # Errors
    /// When [GeomagneticField] cannot be calculated for given argumants, [error::Error] is returned.
    pub fn new(
        height: Length,
        latitude: Angle,
        longitude: Angle,
        date: Date,
    ) -> Result<Self, error::Error> {
        if height < Length::new::<kilometer>(MIN_HEIGHT_KM)
            || height > Length::new::<kilometer>(MAX_HEIGHT_KM)
        {
            return Err(HeightOutsideOfValidityRange {
                min_height_km: MIN_HEIGHT_KM,
                max_height_km: MAX_HEIGHT_KM,
            });
        }

        if latitude < Angle::new::<degree>(-90.0) || latitude > Angle::new::<degree>(90.0) {
            return Err(InvalidLatitude);
        }

        if longitude < Angle::new::<degree>(-180.0) || longitude > Angle::new::<degree>(180.0) {
            return Err(InvalidLongitude);
        }

        let (wmm_model, wmm_error_model) = wmm_models::select_models(date)?;

        let time_delta = date_to_year_decimal(date) - wmm_model.model_version as f32;

        // Translation of coordinates from geodetic to geocentric
        let h = height.get::<meter>();
        let phi = latitude.get::<radian>();
        let lambda = longitude.get::<radian>();
        let r_c = A / sqrtf(1.0 - E_2 * sinf(phi).powi(2)); // Equation (8)
        let p = (r_c + h) * cosf(phi); // Equation (7)
        let z = (r_c * (1.0 - E_2) + h) * sinf(phi); // Equation (7)
        let r = sqrtf(p.powi(2) + z.powi(2)); // Equation (7)
        let phi_prime = asinf(z / r); // Equation (7)

        let psn_sinphi = schmidt_semi_normalised_associated_legendre(sinf(phi_prime));

        let mut g_t = [0.0; 90];
        let mut h_t = [0.0; 90];
        let mut dpsn_sinphi_dphi = [0.0; 90];

        // Prepare matrices
        for n in 1..=12 {
            for m in 0..=n {
                let ix = index(n, m);
                g_t[ix] = wmm_model.g_mfc[ix] + time_delta * wmm_model.g_svc[ix]; // Equation (9): calculating Gauss coefficients for the desired time
                h_t[ix] = wmm_model.h_mfc[ix] + time_delta * wmm_model.h_svc[ix]; // Equation (9): calculating Gauss coefficients for the desired time

                dpsn_sinphi_dphi[ix] = ((n as f64 + 1.0) * tan(phi_prime as f64) * (psn_sinphi[ix] as f64) - // Equation (16): used to calculate the field vector components
                    sqrt((n as f64 + 1.0).powi(2) - (m as f64).powi(2)) / cos(phi_prime as f64) * (psn_sinphi[index(n+1, m)] as f64))
                    as f32;
            }
        }

        // Field vector components
        let mut x_prime = 0.0; // will be computed using Equation (10)
        let mut y_prime = 0.0; // will be computed using Equation (11)
        let mut z_prime = 0.0; // will be computed using Equation (12)
        #[cfg(test)]
        let mut x_dot_prime = 0.0; // will be computed using Equation (13)
        #[cfg(test)]
        let mut y_dot_prime = 0.0; // will be computed using Equation (14)
        #[cfg(test)]
        let mut z_dot_prime = 0.0; // will be computed using Equation (15)

        for n in 1..=12 {
            let mut x_prime_tmp = 0.0;
            let mut y_prime_tmp = 0.0;
            let mut z_prime_tmp = 0.0;
            #[cfg(test)]
            let mut x_dot_prime_tmp = 0.0;
            #[cfg(test)]
            let mut y_dot_prime_tmp = 0.0;
            #[cfg(test)]
            let mut z_dot_prime_tmp = 0.0;

            for m in 0..=n {
                let ix = index(n, m);
                x_prime_tmp += (g_t[ix] * cosf(m as f32 * lambda)
                    + h_t[ix] * sinf(m as f32 * lambda))
                    * dpsn_sinphi_dphi[ix];
                y_prime_tmp += m as f32
                    * (g_t[ix] * sinf(m as f32 * lambda) - h_t[ix] * cosf(m as f32 * lambda))
                    * psn_sinphi[ix];
                z_prime_tmp += (g_t[ix] * cosf(m as f32 * lambda)
                    + h_t[ix] * sinf(m as f32 * lambda))
                    * psn_sinphi[ix];
                #[cfg(test)]
                {
                    x_dot_prime_tmp += (wmm_model.g_svc[ix] * cosf(m as f32 * lambda)
                        + wmm_model.h_svc[ix] * sinf(m as f32 * lambda))
                        * dpsn_sinphi_dphi[ix];
                    y_dot_prime_tmp += m as f32
                        * (wmm_model.g_svc[ix] * sinf(m as f32 * lambda)
                            - wmm_model.h_svc[ix] * cosf(m as f32 * lambda))
                        * psn_sinphi[ix];
                    z_dot_prime_tmp += (wmm_model.g_svc[ix] * cosf(m as f32 * lambda)
                        + wmm_model.h_svc[ix] * sinf(m as f32 * lambda))
                        * psn_sinphi[ix];
                }
            }
            let k_temp = (A2 / r).powi(n as i32 + 2);
            x_prime += k_temp * x_prime_tmp;
            y_prime += k_temp * y_prime_tmp;
            z_prime += (n as f32 + 1.0) * k_temp * z_prime_tmp;

            #[cfg(test)]
            {
                x_dot_prime += k_temp * x_dot_prime_tmp;
                y_dot_prime += k_temp * y_dot_prime_tmp;
                z_dot_prime += (n as f32 + 1.0) * k_temp * z_dot_prime_tmp;
            }
        }

        let x = -x_prime * cosf(phi_prime - phi) + z_prime * sinf(phi_prime - phi);
        let y = y_prime / cosf(phi_prime);
        let z = -x_prime * sinf(phi_prime - phi) - z_prime * cosf(phi_prime - phi);

        #[cfg(test)]
        let x_rate = -x_dot_prime * cosf(phi_prime - phi) + z_dot_prime * sinf(phi_prime - phi);
        #[cfg(test)]
        let y_rate = y_dot_prime / cosf(phi_prime);
        #[cfg(test)]
        let z_rate = -x_dot_prime * sinf(phi_prime - phi) - z_dot_prime * cosf(phi_prime - phi);

        #[cfg(test)]
        let result = Self {
            x,
            y,
            z,
            x_rate,
            y_rate,
            z_rate,
            error_model: wmm_error_model.clone(),
        };

        #[cfg(not(test))]
        let result = Self {
            x,
            y,
            z,
            error_model: wmm_error_model.clone(),
        };

        Ok(result)
    }

    /// Angle between the magnetic flux density (magnetic field) B vector and true north, positive east.
    pub fn declination(&self) -> Angle {
        Angle::new::<radian>(atan2f(self.y, self.x))
    }

    /// Global average error of declination.
    pub fn declination_uncertainty(&self) -> Angle {
        let constant_error = self.error_model.declination_constant_error_factor;
        let variable_error =
            self.error_model.declination_variable_error_factor / self.h().get::<nanotesla>();
        let uncertainty = sqrtf((constant_error).powi(2) + variable_error.powi(2));

        if uncertainty < 180.0 {
            return Angle::new::<degree>(uncertainty);
        }
        Angle::new::<degree>(180.0)
    }

    /// Possible declination warning.
    /// [WarningZone::BlackoutZone] is area around the magnetic pole where
    /// compasses are not accurate and should not be relied on for navigation.
    /// [WarningZone::CautionZone] is area where
    /// compasses should be used with caution because they may not be fully accurate.
    pub fn declination_warning(&self) -> Option<WarningZone> {
        if self.h() < MagneticFluxDensity::new::<nanotesla>(2000.0) {
            Some(WarningZone::BlackoutZone)
        } else if self.h() < MagneticFluxDensity::new::<nanotesla>(6000.0) {
            Some(WarningZone::CautionZone)
        } else {
            None
        }
    }

    /// Angle between the magnetic flux density (magnetic field) B vector and the horizontal plane, positive down.
    pub fn inclination(&self) -> Angle {
        Angle::new::<radian>(atan2f(self.z, self.h().get::<nanotesla>()))
    }

    /// Global average error of inclination.
    pub fn inclination_uncertainty(&self) -> Angle {
        Angle::new::<degree>(self.error_model.inclination_uncertainty)
    }

    /// Horizontal magnetic flux density (magnetic field) B strength.
    pub fn h(&self) -> MagneticFluxDensity {
        MagneticFluxDensity::new::<nanotesla>(sqrtf(self.x.powi(2) + self.y.powi(2)))
    }

    /// Global average error of horizontal magnetic flux density (magnetic field) B strength.
    pub fn h_uncertainty(&self) -> MagneticFluxDensity {
        MagneticFluxDensity::new::<nanotesla>(self.error_model.h_uncertainty)
    }

    /// Geomagnetic flux density (magnetic field) B strength.
    pub fn f(&self) -> MagneticFluxDensity {
        MagneticFluxDensity::new::<nanotesla>(sqrtf(
            self.h().get::<nanotesla>().powi(2) + self.z.powi(2),
        ))
    }

    /// Global average error of geomagnetic flux density (magnetic field) B strength.
    pub fn f_uncertainty(&self) -> MagneticFluxDensity {
        MagneticFluxDensity::new::<nanotesla>(self.error_model.f_uncertainty)
    }

    /// Northern component of the magnetic flux density (magnetic field) B vector.
    pub fn x(&self) -> MagneticFluxDensity {
        MagneticFluxDensity::new::<nanotesla>(self.x)
    }

    /// Global average error of northern component of the magnetic flux density (magnetic field) B vector.
    pub fn x_uncertainty(&self) -> MagneticFluxDensity {
        MagneticFluxDensity::new::<nanotesla>(self.error_model.x_uncertainty)
    }

    /// Eastern component of the magnetic flux density (magnetic field) B vector.
    pub fn y(&self) -> MagneticFluxDensity {
        MagneticFluxDensity::new::<nanotesla>(self.y)
    }

    /// Global average error of eastern component of the magnetic flux density (magnetic field) B vector.
    pub fn y_uncertainty(&self) -> MagneticFluxDensity {
        MagneticFluxDensity::new::<nanotesla>(self.error_model.y_uncertainty)
    }

    /// Downward component of the magnetic flux density (magnetic field) B vector.
    pub fn z(&self) -> MagneticFluxDensity {
        MagneticFluxDensity::new::<nanotesla>(self.z)
    }

    /// Global average error of downward component of the magnetic flux density (magnetic field) B vector.
    pub fn z_uncertainty(&self) -> MagneticFluxDensity {
        MagneticFluxDensity::new::<nanotesla>(self.error_model.z_uncertainty)
    }

    /// Yearly rate of change in declination
    #[cfg(test)]
    fn declination_rate(&self) -> f32 {
        ((self.y_rate * self.x - self.y * self.x_rate) / (self.h().get::<nanotesla>().powi(2)))
            .to_degrees()
    }

    /// Yearly rate of change in inclination
    #[cfg(test)]
    fn inclination_rate(&self) -> f32 {
        ((self.z_rate * self.h().get::<nanotesla>() - self.z * self.h_rate())
            / (self.f().get::<nanotesla>().powi(2)))
        .to_degrees()
    }

    /// Yearly rate of change in horizontal field strength
    #[cfg(test)]
    fn h_rate(&self) -> f32 {
        (self.x * self.x_rate + self.y * self.y_rate) / self.h().get::<nanotesla>()
    }

    /// Yearly rate of change in Magnetic field strength
    #[cfg(test)]
    fn f_rate(&self) -> f32 {
        (self.x * self.x_rate + self.y * self.y_rate + self.z * self.z_rate)
            / self.f().get::<nanotesla>()
    }

    /// Yearly rate of change in the northern component
    #[cfg(test)]
    fn x_rate(&self) -> f32 {
        self.x_rate
    }

    /// Yearly rate of change in the eastern component
    #[cfg(test)]
    fn y_rate(&self) -> f32 {
        self.y_rate
    }

    /// Yearly rate of change in the downward component
    #[cfg(test)]
    fn z_rate(&self) -> f32 {
        self.z_rate
    }
}

/// Function to calculate a decimal year for a given date.
fn date_to_year_decimal(date: Date) -> f32 {
    date.year() as f32
        + (date.ordinal() - 1) as f32 / (time::util::days_in_year(date.year()) as f32)
}

#[cfg(test)]
mod tests {
    extern crate std;
    use super::*;
    use rstest::rstest;
    use time::Month::*;

    #[rstest]
    #[case("ncei.noaa.gov/WMM2020_TestValues.txt")]
    #[case("ncei.noaa.gov/WMM2025_TestValues.txt")]
    fn wmm_tests(#[case] test_file: &str) {
        use assert_float_eq::{assert_float_absolute_eq, assert_float_relative_eq};
        use std::fs::File;
        use std::io::{BufRead, BufReader};
        use std::path::Path;

        let test_values_file =
            File::open(Path::new(test_file)).expect("Test values file not found!");
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
                        year_decimal_to_date(date),
                    )
                    .unwrap();

                    assert_float_relative_eq!(
                        declination,
                        result.declination().get::<degree>(),
                        0.04
                    );
                    assert_float_absolute_eq!(
                        declination,
                        result.declination().get::<degree>(),
                        0.02
                    );

                    assert_float_relative_eq!(
                        inclination,
                        result.inclination().get::<degree>(),
                        0.0009
                    );
                    assert_float_absolute_eq!(
                        inclination,
                        result.inclination().get::<degree>(),
                        0.01
                    );

                    assert_float_relative_eq!(h, result.h().get::<nanotesla>(), 0.0001);
                    assert_float_absolute_eq!(h, result.h().get::<nanotesla>(), 1.0);

                    assert_float_relative_eq!(x, result.x().get::<nanotesla>(), 0.002);
                    assert_float_absolute_eq!(x, result.x().get::<nanotesla>(), 1.0);

                    assert_float_relative_eq!(y, result.y().get::<nanotesla>(), 0.0013);
                    assert_float_absolute_eq!(y, result.y().get::<nanotesla>(), 1.0);

                    assert_float_relative_eq!(z, result.z().get::<nanotesla>(), 0.00005);
                    assert_float_absolute_eq!(z, result.z().get::<nanotesla>(), 1.0);

                    assert_float_relative_eq!(f, result.f().get::<nanotesla>(), 0.04);
                    assert_float_absolute_eq!(f, result.f().get::<nanotesla>(), 1.0);

                    assert_float_absolute_eq!(dd_dt, result.declination_rate(), 0.1);

                    assert_float_absolute_eq!(di_dt, result.inclination_rate(), 0.1);

                    assert_float_absolute_eq!(dh_dt, result.h_rate(), 0.1);

                    assert_float_relative_eq!(dx_dt, result.x_rate(), 0.05);
                    assert_float_absolute_eq!(dx_dt, result.x_rate(), 0.1);

                    assert_float_relative_eq!(dy_dt, result.y_rate(), 0.09);
                    assert_float_absolute_eq!(dy_dt, result.y_rate(), 0.1);

                    assert_float_relative_eq!(dz_dt, result.z_rate(), 0.02);
                    assert_float_absolute_eq!(dz_dt, result.z_rate(), 0.1);

                    assert_float_relative_eq!(df_dt, result.f_rate(), 0.02);
                    assert_float_absolute_eq!(df_dt, result.f_rate(), 0.1);
                }

                Err(error) => {
                    panic!("Error reading lines: {}!", error)
                }
            }
        }
    }

    #[rstest]
    #[case(2023, 15, 0.006997449, 0.21, 145.0, 128.0, 131.0, 94.0, 157.0)]
    #[case(2028, 15, 0.006838081, 0.20, 138.0, 133.0, 137.0, 89.0, 141.0)]
    fn test_uncertainty(
        #[case] year: i32,
        #[case] day: u16,
        #[case] declination_uncertainty: f32,
        #[case] inclination_uncertainty: f32,
        #[case] f_uncertainty: f32,
        #[case] h_uncertainty: f32,
        #[case] x_uncertainty: f32,
        #[case] y_uncertainty: f32,
        #[case] z_uncertainty: f32,
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
            Angle::new::<radian>(declination_uncertainty)
        );
        assert_eq!(
            result.inclination_uncertainty(),
            Angle::new::<degree>(inclination_uncertainty)
        );
        assert_eq!(
            result.f_uncertainty(),
            MagneticFluxDensity::new::<nanotesla>(f_uncertainty)
        );
        assert_eq!(
            result.h_uncertainty(),
            MagneticFluxDensity::new::<nanotesla>(h_uncertainty)
        );
        assert_eq!(
            result.x_uncertainty(),
            MagneticFluxDensity::new::<nanotesla>(x_uncertainty)
        );
        assert_eq!(
            result.y_uncertainty(),
            MagneticFluxDensity::new::<nanotesla>(y_uncertainty)
        );
        assert_eq!(
            result.z_uncertainty(),
            MagneticFluxDensity::new::<nanotesla>(z_uncertainty)
        );
    }

    #[rstest]
    #[case(-90.0, 180.0)]
    #[case(-89.999, 180.0)]
    #[case(-89.99, 180.0)]
    #[case(-89.9, 0.42294145)]
    #[case(-89.0, 0.41796952)]
    #[case(-45.0, 0.5829365)]
    #[case(0.0, 0.33038762)]
    #[case(45.0, 0.35611528)]
    #[case(89.0, 2.4860976)]
    #[case(89.9, 3.1008945)]
    #[case(89.99, 180.0)]
    #[case(89.999, 180.0)]
    #[case(90.0, 180.0)]
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

    #[rstest]
    #[case(-90.0, Some(WarningZone::BlackoutZone))]
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

    #[rstest]
    #[case(-91.0)]
    #[case(-90.1)]
    #[case(-90.00001)]
    #[case(90.00001)]
    #[case(90.1)]
    #[case(91.0)]
    fn test_invalid_latitude(#[case] latitude: f32) {
        let result = GeomagneticField::new(
            Length::new::<meter>(100.0),
            Angle::new::<degree>(latitude),
            Angle::new::<degree>(0.0),
            Date::from_ordinal_date(2023, 15).unwrap(),
        );
        assert_eq!(result.is_err_and(|e| e == InvalidLatitude), true);
    }

    #[rstest]
    #[case(-181.0)]
    #[case(-180.1)]
    #[case(-180.00001)]
    #[case(180.00001)]
    #[case(180.1)]
    #[case(181.0)]
    fn test_invalid_longitude(#[case] longitude: f32) {
        let result = GeomagneticField::new(
            Length::new::<meter>(100.0),
            Angle::new::<degree>(0.0),
            Angle::new::<degree>(longitude),
            Date::from_ordinal_date(2023, 15).unwrap(),
        );
        assert_eq!(result.is_err_and(|e| e == InvalidLongitude), true);
    }

    #[rstest]
    #[case(-101.0)]
    #[case(-100.1)]
    #[case(-100.0001)]
    #[case(850.001)]
    #[case(850.1)]
    #[case(851.0)]
    fn test_height_outside_of_validity_range(#[case] height: f32) {
        let result = GeomagneticField::new(
            Length::new::<kilometer>(height),
            Angle::new::<degree>(0.0),
            Angle::new::<degree>(0.0),
            Date::from_ordinal_date(2023, 15).unwrap(),
        );
        assert_eq!(
            result.is_err_and(|e| e
                == HeightOutsideOfValidityRange {
                    min_height_km: -1.0,
                    max_height_km: 850.0
                }),
            true
        );
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
        use error::Error::DateOutsideOfValidityRange;

        let result = GeomagneticField::new(
            Length::new::<kilometer>(0.0),
            Angle::new::<degree>(0.0),
            Angle::new::<degree>(0.0),
            Date::from_calendar_date(year, month, day).unwrap(),
        );
        assert_eq!(result.is_err_and(|e| e == DateOutsideOfValidityRange), true);
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

    fn year_decimal_to_date(year_decimal: f32) -> Date {
        let year_int = year_decimal as i32;
        let day_ordinal =
            f32::round(year_decimal.fract() * time::util::days_in_year(year_int) as f32);
        Date::from_ordinal_date(year_int, day_ordinal as u16 + 1).unwrap()
    }

    #[rstest]
    #[case(2020, July, 2, 2020.5)]
    #[case(2021, July, 3, 2021.5013698630137)]
    #[case(2023, January, 1, 2023.0)]
    #[case(2023, January, 2, 2023.0027397260274)]
    #[case(2024, January, 2, 2024.0027322404371)]
    #[case(2024, December, 12, 2024.9453551912568)]
    #[case(2024, December, 31, 2024.9972677595629)]
    #[case(2025, January, 1, 2025.0)]
    fn test_year_decimal_conversions(
        #[case] year_int: i32,
        #[case] month: time::Month,
        #[case] day: u8,
        #[case] year_decimal: f32,
    ) {
        assert_eq!(
            date_to_year_decimal(Date::from_calendar_date(year_int, month, day).unwrap()),
            year_decimal
        );
        assert_eq!(
            year_decimal_to_date(year_decimal),
            Date::from_calendar_date(year_int, month, day).unwrap()
        );
        assert_eq!(
            date_to_year_decimal(year_decimal_to_date(year_decimal)),
            year_decimal
        );
        assert_eq!(
            year_decimal_to_date(date_to_year_decimal(
                Date::from_calendar_date(year_int, month, day).unwrap()
            )),
            Date::from_calendar_date(year_int, month, day).unwrap()
        );
    }

    #[test]
    fn test_send() {
        fn assert_send<T: Send>() {}
        assert_send::<GeomagneticField>();
    }

    #[test]
    fn test_sync() {
        fn assert_sync<T: Sync>() {}
        assert_sync::<GeomagneticField>();
    }
}
