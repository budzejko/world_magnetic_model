#![no_std]
#![forbid(unsafe_code)]
#![warn(missing_docs)]
//! This crate is a Rust implementation of the NOAA [World Magnetic Model (WMM)](https://www.ncei.noaa.gov/products/world-magnetic-model),
//! a mathematical representation of the Earth's core magnetic field and its temporal variations.
//!
//! The crate's interface utilizes the [uom (Units of Measurement) crate](https://docs.rs/uom/latest/uom/) to represent physical quantities
//! accurately. WMM coefficient files are converted into code constants to eliminate the need for file reading at runtime.
//! This crate is compatible with `no_std` environments, meaning it does not depend on the Rust standard library and can be used in embedded,
//! bare-metal, or other restricted contexts by relying on the core crate instead.
//!
//! # Bundled models
//!
//! Each release ships two consecutive WMM coefficient sets (currently 2020 and 2025),
//! selected automatically from the evaluation date. A successor model is published while
//! the current epoch is still in force. Keeping both the outgoing and incoming sets in
//! the crate lets embedded firmware reserve a constant footprint for coefficients: an
//! update that replaces the older epoch with the newly published one does not grow the
//! image, and devices already carry the next model when the epoch boundary arrives.
//!
//! Crate versions follow `api.epoch.patch`, where `epoch` is the year of the newest
//! bundled WMM model.
//!
//! # Usage
//!
//! Results are shown as \(D \pm \sigma_D\) (1σ). Field accuracy is the
//! [WMM error model](https://www.ncei.noaa.gov/products/world-magnetic-model/accuracy-limitations-error-model),
//! not the `f32` precision.
//!
//! ```rust
//! # use core::error::Error;
//! # use world_magnetic_model::time::Date;
//! # use world_magnetic_model::uom::si::f32::{Angle, Length};
//! # use world_magnetic_model::uom::si::angle::degree;
//! # use world_magnetic_model::uom::si::length::meter;
//! # use world_magnetic_model::GeomagneticField;
//! # fn main() -> Result<(), Box<dyn Error>> {
//! let geomagnetic_field = GeomagneticField::new(
//!     Length::new::<meter>(100.0), // Height above the WGS 84 ellipsoid
//!     Angle::new::<degree>(37.03), // WGS 84 latitude (negative values for the Southern Hemisphere)
//!     Angle::new::<degree>(-7.91), // WGS 84 longitude (negative values for the Western Hemisphere)
//!     Date::from_ordinal_date(2029, 15)? // Date (15th day of 2029)
//! )?;
//!
//! let d = geomagnetic_field.declination().get::<degree>();
//! let sigma = geomagnetic_field.declination_uncertainty().get::<degree>();
//!
//! assert_eq!(format!("{d:.2} ± {sigma:.2}°"), "-0.17 ± 0.33°");
//! #     Ok(())
//! # }
//! ```
//!
//! # Cargo features
//!
//! Enable **`serde`** for `Serialize`/`Deserialize` on [`GeomagneticField`] and
//! [`WarningZone`]. `GeomagneticField` encodes the constructor query plus public
//! results (declination, intensity, warning, …). Deserialization uses the
//! inputs only and recomputes through [`GeomagneticField::new`]. Enabling
//! `serde` does not pull in `std`. Encode with any serde codec, for example
//! `serde_json` or `postcard`.
//!
//! ```
//! # #[cfg(not(feature = "serde"))]
//! # fn main() {}
//! # #[cfg(feature = "serde")]
//! # fn main() -> Result<(), Box<dyn core::error::Error>> {
//! # use world_magnetic_model::time::Date;
//! # use world_magnetic_model::uom::si::f32::{Angle, Length};
//! # use world_magnetic_model::uom::si::angle::degree;
//! # use world_magnetic_model::uom::si::length::meter;
//! # use world_magnetic_model::GeomagneticField;
//! # let geomagnetic_field = GeomagneticField::new(
//! #     Length::new::<meter>(100.0),
//! #     Angle::new::<degree>(37.03),
//! #     Angle::new::<degree>(-7.91),
//! #     Date::from_ordinal_date(2029, 15)?,
//! # )?;
//! let json = serde_json::to_string(&geomagnetic_field)?;
//! let from_json: GeomagneticField = serde_json::from_str(&json)?;
//!
//! let mut buf = [0u8; 256];
//! let bytes = postcard::to_slice(&geomagnetic_field, &mut buf)?;
//! let from_postcard: GeomagneticField = postcard::from_bytes(bytes)?;
//!
//! assert_eq!(geomagnetic_field.declination(), from_json.declination());
//! assert_eq!(geomagnetic_field.declination(), from_postcard.declination());
//! # Ok(())
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
//! This work was inspired by [geomag-wmm](https://crates.io/crates/geomag-wmm/0.1.0).
//!
//! # License
//! Licensed under either of [Apache License, Version 2.0](https://github.com/budzejko/world_magnetic_model/blob/main/LICENSE-APACHE)
//! or [MIT license](https://github.com/budzejko/world_magnetic_model/blob/main/LICENSE-MIT) at your option.

pub use error::Error;
pub use time;
pub use uom;

mod error;
#[cfg(feature = "serde")]
mod serde_support;
mod synthesis;
mod wmm_data;
mod wmm_models;

use libm::{atan2f, cosf, sinf, sqrtf};
use time::Date;
// `f32::powi` in `no_std` comes from this trait; unused when the std test harness is linked.
#[allow(unused)]
use uom::num_traits::float::FloatCore;
use uom::si::angle::{degree, radian};
use uom::si::f32::{Angle, Length, MagneticFluxDensity};
use uom::si::length::{kilometer, meter};
use uom::si::magnetic_flux_density::nanotesla;

use error::Error::{HeightOutsideOfValidityRange, InvalidLatitude, InvalidLongitude};
use synthesis::{N_COEFF, sin_cos_m_lambda, synthesize_near_pole, synthesize_xyz_prime};
use wmm_models::WmmErrorModel;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

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
///     Length::new::<foot>(3000.0), // Height above the WGS 84 ellipsoid
///     Angle::new::<degree>(80.0), // WGS 84 latitude (negative values for the Southern Hemisphere)
///     Angle::new::<radian>(2.36), // WGS 84 longitude (negative values for the Western Hemisphere)
///     Date::from_ordinal_date(2029, 15)? // Date (15th day of 2029)
/// )?;
///
/// // declination results can be interpreted as degrees or radians
/// assert_eq!(
///     format!("{:.2}", geomagnetic_field.declination().get::<degree>()),
///     "-26.41"
/// );
/// assert_eq!(
///     format!("{:.3}", geomagnetic_field.declination().get::<radian>()),
///     "-0.461"
/// );
/// assert_eq!(
///     format!("{:.2}", geomagnetic_field.declination_uncertainty().get::<degree>()),
///     "2.07"
/// );
/// assert_eq!(
///     format!("{:.3}", geomagnetic_field.declination_uncertainty().get::<radian>()),
///     "0.036"
/// );
///
/// // warning conditions can be checked e.g. for user interface warning
/// assert_eq!(
///     geomagnetic_field.declination_warning(),
///     Some(CautionZone)
/// );
///
/// // grid variation (UPS) is defined at |latitude| ≥ 55°
/// assert_eq!(
///     format!(
///         "{:.2}",
///         geomagnetic_field.grid_variation().unwrap().get::<degree>()
///     ),
///     "-161.63"
/// );
///
/// // magnetic flux density results can be interpreted as nT or gauss
/// assert_eq!(
///     format!("{:.0}", geomagnetic_field.f().get::<nanotesla>()),
///     "59021"
/// );
/// assert_eq!(
///     format!("{:.3}", geomagnetic_field.f().get::<gauss>()),
///     "0.590"
/// );
///
/// // results can be formatted with unit abbreviation or full description
/// assert_eq!(
///     format!("{:.2}", geomagnetic_field.f().into_format_args(nanotesla, Abbreviation)),
///     "59020.67 nT"
/// );
/// assert_eq!(
///     format!("{:.2}", geomagnetic_field.f().into_format_args(nanotesla, Description)),
///     "59020.67 nanoteslas"
/// );
/// assert_eq!(
///     format!("{:.2}", geomagnetic_field.f().into_format_args(gauss, Abbreviation)),
///     "0.59 G"
/// );
/// assert_eq!(
///     format!("{:.2}", geomagnetic_field.f().into_format_args(gauss, Description)),
///     "0.59 gauss"
/// );
/// #     Ok(())
/// # }
/// ```
///
/// # Serialization
///
/// With the `serde` feature, this type serializes the original query together
/// with the public results:
///
/// - `height_m`, `latitude_deg`, `longitude_deg`, `date` (`YYYY-MM-DD`)
/// - `x_nT`, `y_nT`, `z_nT`, `h_nT`, `f_nT`
/// - `declination_deg`, `inclination_deg`
/// - `declination_uncertainty_deg`, `inclination_uncertainty_deg`
/// - `x_uncertainty_nT`, `y_uncertainty_nT`, `z_uncertainty_nT`,
///   `h_uncertainty_nT`, `f_uncertainty_nT`
/// - `declination_warning`, `grid_variation_deg` (`null` when undefined)
///
/// Deserialization reads only the query and calls [`GeomagneticField::new`],
/// so output fields in the payload are ignored. Any serde format works; JSON and
/// `postcard` are covered by the crate tests.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(
    feature = "serde",
    derive(Serialize, Deserialize),
    serde(
        try_from = "serde_support::GeomagneticFieldInput",
        into = "serde_support::GeomagneticFieldView"
    )
)]
pub struct GeomagneticField {
    x: f32,             // Northern component X [nT]
    y: f32,             // Eastern component Y [nT]
    z: f32,             // Downward component Z [nT]
    height_m: f32,      // Height above the WGS 84 ellipsoid [m]
    latitude_deg: f32,  // WGS 84 geodetic latitude [deg], for UPS grid variation
    longitude_deg: f32, // WGS 84 geodetic longitude [deg], for UPS grid variation
    date: Date,         // Evaluation date (selects the WMM epoch)
    #[cfg(test)]
    x_rate: f32, // Yearly rate of change in X [nT/year]
    #[cfg(test)]
    y_rate: f32, // Yearly rate of change in Y [nT/year]
    #[cfg(test)]
    z_rate: f32, // Yearly rate of change in Z [nT/year]
    error_model: WmmErrorModel, // [WMM2020 TR](https://doi.org/10.25923/ytk1-yx35) §3.4; [WMM2025 TR](https://doi.org/10.25923/prbc-s316) §3.4: σ_D is location-dependent; others are global RMS
}

/// Warning zone type, as defined in WMM2025 TR §1.8.
#[derive(Debug, Copy, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum WarningZone {
    /// Based on the WMM military specification, "Blackout Zones" (BoZ) are areas around
    /// the north and south magnetic poles where compasses are not accurate and should not be
    /// relied on for navigation. The BoZ are defined as regions around the north and south
    /// magnetic poles where the horizontal intensity of Earth's magnetic field (H) is less
    /// than 2000 nT. In BoZs, WMM declination values are not accurate and compasses are
    /// unreliable.
    BlackoutZone,
    /// "Caution Zone" (2000 nT <= H < 6000 nT) is an area around the perimeter of the BoZs where
    /// compasses should be used with caution because they may not be fully accurate.
    CautionZone,
}

const MIN_HEIGHT_KM: f32 = -1.0;
const MAX_HEIGHT_KM: f32 = 850.0;
const MIN_LATITUDE_DEG: f32 = -90.0;
const MAX_LATITUDE_DEG: f32 = 90.0;
const MIN_LONGITUDE_DEG: f32 = -180.0;
const MAX_LONGITUDE_DEG: f32 = 180.0;
/// WMM2025 TR §1.8 \(\sigma_D\) clamp; declination cannot exceed \(180^\circ\).
const MAX_DECLINATION_UNCERTAINTY_DEG: f32 = 180.0;
/// Horizontal intensity \(H\) below which compasses are unreliable (WMM2025 TR §1.8 Blackout Zone).
const BLACKOUT_ZONE_H_NT: f32 = 2000.0;
/// Horizontal intensity \(H\) below which compasses should be used with caution (WMM2025 TR §1.8 Caution Zone).
const CAUTION_ZONE_H_NT: f32 = 6000.0;
/// \(\lvert\varphi\rvert\) at and beyond which UPS grid variation is defined
/// ([NOAA WMM C](https://www.ncei.noaa.gov/products/world-magnetic-model) `MAG_PS_*_LAT_DEGREE`; [MIL-PRF-89500B](https://www.ncei.noaa.gov/sites/default/files/2025-09/MIL-PRF-89500B.pdf) uses a strict inequality).
const GRID_VARIATION_MIN_ABS_LATITUDE_DEG: f32 = 55.0;

// WGS 84 ellipsoid (WMM2025 TR §1.2)
const WGS84_A: f32 = 6378137.0; // Semi-major axis a [m]
const WGS84_F: f32 = (1.0_f64 / 298.257223563) as f32;
const WGS84_E2: f32 = WGS84_F * (2.0 - WGS84_F); // First eccentricity squared

// WMM geomagnetic reference radius a (WMM2025 TR §1.2)
const WMM_A: f32 = 6371200.0;

// NOAA WMM C `MAG_SummationSpecial` uses 1e-10 in double; f32 loses accuracy in Y'/sinθ
// and P_n^m ~ sin^m θ well before that, so the polar path is taken sooner.
const NEAR_POLE_SIN_THETA_MIN: f32 = 1e-3;

/// Geodetic WGS 84 \((h, \varphi)\) to geocentric radius and direction sines/cosines.
///
/// # Arguments
///
/// * `height` - Height above the WGS 84 ellipsoid in meters.
/// * `latitude` - WGS 84 geodetic latitude in radians.
///
/// # Returns
///
/// * `r` - Geocentric radius in meters
/// * `cos_theta` - \(\cos\theta = \sin\varphi' = z/r\) (geocentric)
/// * `sin_theta` - \(\sin\theta = \cos\varphi' = p/r\) (geocentric)
/// * `sin_phi` / `cos_phi` - geodetic \(\sin\varphi\), \(\cos\varphi\)
fn geodetic_to_geocentric(height: f32, latitude: f32) -> (f32, f32, f32, f32, f32) {
    let sin_phi = sinf(latitude);
    let cos_phi = cosf(latitude);
    let r_c = WGS84_A / sqrtf(1.0 - WGS84_E2 * sin_phi * sin_phi); // prime vertical radius of curvature
    let p = (r_c + height) * cos_phi; // geocentric cylindrical radius
    let z = (r_c * (1.0 - WGS84_E2) + height) * sin_phi;
    let r = sqrtf(p * p + z * z);
    let cos_theta = z / r;
    let sin_theta = p / r;
    (r, cos_theta, sin_theta, sin_phi, cos_phi)
}

/// East \(Y\) together with φ'-based geocentric \(X'\), \(Z'\) (opposite sign vs
/// WMM2025 TR §1.2 \(X'\), \(Z'\)). Near the geographic poles \(Y\) uses the finite
/// \(P_n^1/\sin\theta\) recurrence instead of \(Y'/\sin\theta\).
#[allow(clippy::too_many_arguments)]
fn synthesize_geocentric_xyz(
    cos_theta: f32,
    sin_theta: f32,
    lambda: f32,
    ratio: f32,
    g: &[f32; N_COEFF],
    h: &[f32; N_COEFF],
    g_dot: &[f32; N_COEFF],
    h_dot: &[f32; N_COEFF],
    time_delta: f32,
) -> (f32, f32, f32) {
    if sin_theta.abs() < NEAR_POLE_SIN_THETA_MIN {
        synthesize_near_pole(
            cos_theta, sin_theta, lambda, ratio, g, h, g_dot, h_dot, time_delta,
        )
    } else {
        let (sin_m, cos_m) = sin_cos_m_lambda(lambda);
        let (x_prime, y_prime, z_prime) = synthesize_xyz_prime(
            cos_theta, sin_theta, &sin_m, &cos_m, ratio, g, h, g_dot, h_dot, time_delta,
        );
        (x_prime, y_prime / sin_theta, z_prime)
    }
}

impl GeomagneticField {
    /// Calculates the geomagnetic field at a given location and date.
    ///
    /// # Arguments
    ///
    /// * `height` - Height above World Geodetic System 1984 (WGS 84) ellipsoid (from -1km to +850km).
    /// * `latitude` - WGS 84 latitude (negative values for the Southern Hemisphere).
    /// * `longitude` - WGS 84 longitude (negative values for the Western Hemisphere).
    /// * `date` - Date for the WMM model evaluation.
    ///
    /// # Errors
    /// When [GeomagneticField] cannot be calculated for given arguments, [error::Error] is returned.
    pub fn new(
        height: Length,
        latitude: Angle,
        longitude: Angle,
        date: Date,
    ) -> Result<Self, error::Error> {
        if !height.get::<kilometer>().is_finite()
            || height < Length::new::<kilometer>(MIN_HEIGHT_KM)
            || height > Length::new::<kilometer>(MAX_HEIGHT_KM)
        {
            return Err(HeightOutsideOfValidityRange {
                min_height_km: MIN_HEIGHT_KM,
                max_height_km: MAX_HEIGHT_KM,
            });
        }

        if !latitude.get::<degree>().is_finite()
            || latitude < Angle::new::<degree>(MIN_LATITUDE_DEG)
            || latitude > Angle::new::<degree>(MAX_LATITUDE_DEG)
        {
            return Err(InvalidLatitude);
        }

        if !longitude.get::<degree>().is_finite()
            || longitude < Angle::new::<degree>(MIN_LONGITUDE_DEG)
            || longitude > Angle::new::<degree>(MAX_LONGITUDE_DEG)
        {
            return Err(InvalidLongitude);
        }

        let (wmm_model, wmm_error_model) = wmm_models::select_models(date)?;

        let time_delta = years_since_epoch(date, wmm_model.epoch_year);

        // Geodetic WGS 84 to geocentric spherical
        let height_m = height.get::<meter>();
        let phi = latitude.get::<radian>();
        let lambda = longitude.get::<radian>();
        let (r, cos_theta, sin_theta, sin_phi, cos_phi) = geodetic_to_geocentric(height_m, phi);

        let ratio = WMM_A / r;
        let (x_prime, y, z_prime) = synthesize_geocentric_xyz(
            cos_theta,
            sin_theta,
            lambda,
            ratio,
            &wmm_model.g,
            &wmm_model.h,
            &wmm_model.g_dot,
            &wmm_model.h_dot,
            time_delta,
        );
        // Secular rates: pass \(\dot g,\dot h\) as the main field with zero epoch offset
        // so the synthesizer returns \(\dot X',\dot Y,\dot Z'\) (Y is already east).
        #[cfg(test)]
        let (x_dot_prime, y_dot, z_dot_prime) = synthesize_geocentric_xyz(
            cos_theta,
            sin_theta,
            lambda,
            ratio,
            &wmm_model.g_dot,
            &wmm_model.h_dot,
            &wmm_model.g_dot,
            &wmm_model.h_dot,
            0.0,
        );

        // Synthesis X', Z' equal -X'_TR, -Z'_TR (φ' vs θ). The minuses recover
        // X = X'_TR cos(φ'−φ) − Z'_TR sin(φ'−φ), Z = X'_TR sin(φ'−φ) + Z'_TR cos(φ'−φ).
        let sin_dphi = cos_theta * cos_phi - sin_theta * sin_phi;
        let cos_dphi = sin_theta * cos_phi + cos_theta * sin_phi;
        let x = -x_prime * cos_dphi + z_prime * sin_dphi;
        let z = -x_prime * sin_dphi - z_prime * cos_dphi;

        #[cfg(test)]
        let x_rate = -x_dot_prime * cos_dphi + z_dot_prime * sin_dphi;
        #[cfg(test)]
        let z_rate = -x_dot_prime * sin_dphi - z_dot_prime * cos_dphi;

        let latitude_deg = latitude.get::<degree>();
        let longitude_deg = longitude.get::<degree>();

        #[cfg(test)]
        let result = Self {
            x,
            y,
            z,
            height_m,
            latitude_deg,
            longitude_deg,
            date,
            x_rate,
            y_rate: y_dot,
            z_rate,
            error_model: *wmm_error_model,
        };

        #[cfg(not(test))]
        let result = Self {
            x,
            y,
            z,
            height_m,
            latitude_deg,
            longitude_deg,
            date,
            error_model: *wmm_error_model,
        };

        Ok(result)
    }

    /// Magnetic declination \(D\): angle from true north to the horizontal
    /// component of B, positive east (WMM2025 TR §1.1.1).
    pub fn declination(&self) -> Angle {
        Angle::new::<radian>(atan2f(self.y, self.x))
    }

    /// Grid variation (grivation) relative to UPS grid north, or `None`
    /// outside the polar regions.
    ///
    /// Near the geographic poles \(D\) changes rapidly with longitude, so WMM
    /// defines an auxiliary angle against Universal Polar Stereographic (UPS)
    /// grid north (MIL-PRF-89500B §A.2.6):
    ///
    /// ```text
    /// GV = D − λ    for φ ≥ 55°
    /// GV = D + λ    for φ ≤ −55°
    /// ```
    ///
    /// \(\varphi\) and \(\lambda\) are geodetic latitude and longitude.
    /// The result is wrapped to \([-180^\circ, 180^\circ]\). Mid-latitudes
    /// return `None`; this crate does not compute UTM grivation (only UPS).
    pub fn grid_variation(&self) -> Option<Angle> {
        let gv = if self.latitude_deg >= GRID_VARIATION_MIN_ABS_LATITUDE_DEG {
            self.declination().get::<degree>() - self.longitude_deg
        } else if self.latitude_deg <= -GRID_VARIATION_MIN_ABS_LATITUDE_DEG {
            self.declination().get::<degree>() + self.longitude_deg
        } else {
            return None;
        };
        Some(Angle::new::<degree>(wrap_to_180(gv)))
    }

    /// Location-dependent 1σ uncertainty of declination (WMM2020 TR §3.4; WMM2025 TR §3.4).
    ///
    /// ```text
    /// σ_D = √(σ_c² + (σ_v / H)²)
    /// ```
    ///
    /// `σ_c` and `σ_v` are epoch constants of the selected WMM;
    /// `H` is the local horizontal intensity. The result is clamped to 180°.
    ///
    /// Unlike the other `*_uncertainty` methods, this is not a global RMS:
    /// it grows as `H` decreases (near the magnetic poles).
    pub fn declination_uncertainty(&self) -> Angle {
        let sigma_c = self.error_model.declination_sigma_c;
        let sigma_v = self.error_model.declination_sigma_v / self.h().get::<nanotesla>();
        let uncertainty = sqrtf(sigma_c.powi(2) + sigma_v.powi(2));

        if uncertainty < MAX_DECLINATION_UNCERTAINTY_DEG {
            return Angle::new::<degree>(uncertainty);
        }
        Angle::new::<degree>(MAX_DECLINATION_UNCERTAINTY_DEG)
    }

    /// Possible declination warning.
    /// [WarningZone::BlackoutZone] is an area around the magnetic pole where
    /// compasses are not accurate and should not be relied on for navigation.
    /// [WarningZone::CautionZone] is an area where
    /// compasses should be used with caution because they may not be fully accurate.
    #[must_use]
    pub fn declination_warning(&self) -> Option<WarningZone> {
        if self.h() < MagneticFluxDensity::new::<nanotesla>(BLACKOUT_ZONE_H_NT) {
            Some(WarningZone::BlackoutZone)
        } else if self.h() < MagneticFluxDensity::new::<nanotesla>(CAUTION_ZONE_H_NT) {
            Some(WarningZone::CautionZone)
        } else {
            None
        }
    }

    /// Magnetic inclination \(I\): angle between the horizontal plane and B,
    /// positive down (WMM2025 TR §1.1.1).
    pub fn inclination(&self) -> Angle {
        Angle::new::<radian>(atan2f(self.z, self.h().get::<nanotesla>()))
    }

    /// Global RMS (1σ) uncertainty of inclination (WMM2020 TR §3.4; WMM2025 TR §3.4).
    ///
    /// Constant for the selected WMM epoch; does not depend on location.
    pub fn inclination_uncertainty(&self) -> Angle {
        Angle::new::<degree>(self.error_model.inclination_uncertainty)
    }

    /// Horizontal intensity \(H\): magnitude of the horizontal component of B
    /// (not the Gauss coefficient \(h_n^m\)).
    pub fn h(&self) -> MagneticFluxDensity {
        MagneticFluxDensity::new::<nanotesla>(sqrtf(self.x.powi(2) + self.y.powi(2)))
    }

    /// Global RMS (1σ) uncertainty of horizontal intensity H (WMM2020 TR §3.4; WMM2025 TR §3.4).
    ///
    /// Constant for the selected WMM epoch; does not depend on location.
    pub fn h_uncertainty(&self) -> MagneticFluxDensity {
        MagneticFluxDensity::new::<nanotesla>(self.error_model.h_uncertainty)
    }

    /// Total intensity \(F\): magnitude of the magnetic flux density (magnetic field) B.
    pub fn f(&self) -> MagneticFluxDensity {
        MagneticFluxDensity::new::<nanotesla>(sqrtf(
            self.h().get::<nanotesla>().powi(2) + self.z.powi(2),
        ))
    }

    /// Global RMS (1σ) uncertainty of total intensity F (WMM2020 TR §3.4; WMM2025 TR §3.4).
    ///
    /// Constant for the selected WMM epoch; does not depend on location.
    pub fn f_uncertainty(&self) -> MagneticFluxDensity {
        MagneticFluxDensity::new::<nanotesla>(self.error_model.f_uncertainty)
    }

    /// Northern component of the magnetic flux density (magnetic field) B vector.
    pub fn x(&self) -> MagneticFluxDensity {
        MagneticFluxDensity::new::<nanotesla>(self.x)
    }

    /// Global RMS (1σ) uncertainty of the northern component X (WMM2020 TR §3.4; WMM2025 TR §3.4).
    ///
    /// Constant for the selected WMM epoch; does not depend on location.
    pub fn x_uncertainty(&self) -> MagneticFluxDensity {
        MagneticFluxDensity::new::<nanotesla>(self.error_model.x_uncertainty)
    }

    /// Eastern component of the magnetic flux density (magnetic field) B vector.
    pub fn y(&self) -> MagneticFluxDensity {
        MagneticFluxDensity::new::<nanotesla>(self.y)
    }

    /// Global RMS (1σ) uncertainty of the eastern component Y (WMM2020 TR §3.4; WMM2025 TR §3.4).
    ///
    /// Constant for the selected WMM epoch; does not depend on location.
    pub fn y_uncertainty(&self) -> MagneticFluxDensity {
        MagneticFluxDensity::new::<nanotesla>(self.error_model.y_uncertainty)
    }

    /// Downward component of the magnetic flux density (magnetic field) B vector.
    pub fn z(&self) -> MagneticFluxDensity {
        MagneticFluxDensity::new::<nanotesla>(self.z)
    }

    /// Global RMS (1σ) uncertainty of the downward component Z (WMM2020 TR §3.4; WMM2025 TR §3.4).
    ///
    /// Constant for the selected WMM epoch; does not depend on location.
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

    /// Yearly rate of change in horizontal intensity \(H\)
    #[cfg(test)]
    fn h_rate(&self) -> f32 {
        (self.x * self.x_rate + self.y * self.y_rate) / self.h().get::<nanotesla>()
    }

    /// Yearly rate of change in total intensity F
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

/// Wrap an angle in degrees to \([-180^\circ, 180^\circ]\) (legacy NOAA WMM C `geomag.c`).
fn wrap_to_180(deg: f32) -> f32 {
    if deg > 180.0 {
        deg - 360.0
    } else if deg < -180.0 {
        deg + 360.0
    } else {
        deg
    }
}

/// Offset \(t-t_0\) in years (WMM2025 TR §1.2).
///
/// Integer year offset plus the day fraction, so the sum stays \(O(1)\) and
/// `f32` keeps ~\(10^{-5}\) day. Forming the calendar year first (e.g. `2024.945`)
/// then subtracting \(t_0\) loses ~0.02 day to the ULP of ~2025.
fn years_since_epoch(date: Date, epoch_year: i32) -> f32 {
    (date.year() - epoch_year) as f32
        + (date.ordinal() - 1) as f32 / (time::util::days_in_year(date.year()) as f32)
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod property_tests;
