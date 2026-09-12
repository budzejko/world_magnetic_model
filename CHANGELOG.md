# Changelog

All notable changes to this project are documented here.

Published crate versions: [crates.io/crates/world_magnetic_model](https://crates.io/crates/world_magnetic_model).

From 1.2025.0 the crate uses `api.epoch.patch`, where `epoch` is the year of
the newest bundled WMM coefficient set.

## [1.2025.0] - 2026-09-12

[crates.io](https://crates.io/crates/world_magnetic_model/1.2025.0)

First 1.x release. Still ships WMM2020 and WMM2025.

### Breaking

- Version scheme is `api.epoch.patch` instead of 0.x (`0.4.0` → `1.2025.0`).
- Rust 2024 edition; MSRV is **1.85**.
- `Error::DateOutsideOfValidityRange` is a struct variant with `min_date` and
  `max_date` (was a unit variant). Display includes that inclusive window
  (`2020-01-01 to 2029-12-31` for the current pair of models).
- `Error::InvalidLatitude` / `InvalidLongitude` display the accepted ranges.
- Re-exported `uom` is **0.38** (was 0.37).
- Optional `serde` no longer uses a default `Serialize`/`Deserialize` derive on
  the internal fields. The wire format is a documented query-plus-results
  object; deserialize reads only the query and calls `GeomagneticField::new`.
  Existing payloads, if any, will not deserialize.
- `GeomagneticField: PartialEq` now includes the constructor query
  (`height`, latitude, longitude, date) as well as \(X\), \(Y\), \(Z\) and the
  error model. Values that compared equal in 0.4.0 can differ if they came
  from different inputs that produced the same field.

### Added

- `GeomagneticField::grid_variation()` — UPS grivation for
  \(\lvert\varphi\rvert \ge 55^\circ\)
  ([NOAA WMM C](https://www.ncei.noaa.gov/products/world-magnetic-model)
  `MAG_PS_*_LAT_DEGREE`). The angle is
  [MIL-PRF-89500B](https://www.ncei.noaa.gov/sites/default/files/2025-09/MIL-PRF-89500B.pdf)
  §A.2.6 (that spec uses a strict inequality). Mid-latitudes return `None`
  (UTM grivation is not computed). The result is wrapped to
  \([-180^\circ, 180^\circ]\).
- `Error` is `Copy` + `Clone`. `WarningZone` is `Copy`.
- Real `serde` support (`default-features = false`, still `no_std`): JSON and
  `postcard` are covered by tests; `demo_serde` example.
- `GeomagneticField::new` rejects non-finite height, latitude, and longitude.
- `demo_sv_extrapolation` example (full WMM vs linear secular-variation
  extrapolation).
- Property tests for constructor invariants (finite outputs, angle ranges).
- `rust-version` in `Cargo.toml`; `#![forbid(unsafe_code)]`.
- docs.rs builds with `all-features`.

### Changed

- Field synthesis rewritten around Schmidt quasi-normalized associated
  Legendre recurrences. Replaces the slower factorial / closed-form path
  in `math.rs`.
- WGS 84 flattening is the standard \(1/298.257223563\) (was \(1/298.25723\)).
- Main-field components and derived angles can differ from 0.4.0 at `f32`
  precision (more near the geographic poles) after the recurrence rewrite,
  polar path, and flattening update.
- Gauss coefficients live in generated `wmm_data.rs`; synthesis is in
  `synthesis.rs`.
- Declination uncertainty documents the
  [WMM2020 TR](https://doi.org/10.25923/ytk1-yx35) §3.4 /
  [WMM2025 TR](https://doi.org/10.25923/prbc-s316) §3.4 model
  \(\sigma_D = \sqrt{\sigma_c^2 + (\sigma_v / H)^2}\) and the \(180^\circ\) clamp.
- Usage docs report \(D \pm \sigma_D\) instead of raw `f32` equality.
- Dependencies: `time` 0.3.55, `thiserror` 2.0.20, `libm` 0.2.16,
  `askama` 0.16.1; dev: `rstest` 0.27, `criterion` 0.8, `assert_float_eq` 1.2.

### Fixed

- Dedicated near-pole summation avoids \(Y'/\sin\theta\) blow-up.
- `+180°` and `-180°` longitude are treated as the same meridian.
- Optional `serde` is usable in `no_std` (previously pulled `std` via default
  features) and no longer depends on `time`/`uom` serde impls. The serialized
  snapshot includes all public 1σ uncertainties (\(\sigma_X,\sigma_Y,\sigma_Z,
  \sigma_H,\sigma_F\) as well as \(\sigma_D\) and \(\sigma_I\)).
- NOAA official test-value files are skipped when absent instead of failing
  the suite.

## [0.4.0] - 2025-06-15

[crates.io](https://crates.io/crates/world_magnetic_model/0.4.0)

### Changed

- Schmidt associated Legendre evaluation uses a precomputed factorial table
  and a closed-form \((n, m)\) index (no runtime `1..=n` fold).

## [0.3.0] - 2025-06-08

[crates.io](https://crates.io/crates/world_magnetic_model/0.3.0)

### Changed

- Re-exported `uom` is **0.37** (was 0.36).
- Other deps: `time` 0.3.41, `thiserror` 2.0.12, `libm` 0.2.14,
  `serde` 1.0.219; build `rinja` → `askama` 0.14; `criterion` 0.6,
  `rstest` 0.25.
- `mise.toml` excluded from the published crate; `Cargo.lock` dropped from
  the repo.
- `serde` imports are gated on the `serde` feature.

## [0.2.0] - 2025-01-19

[crates.io](https://crates.io/crates/world_magnetic_model/0.2.0)

### Added

- `time` and `uom` are re-exported (`world_magnetic_model::time`,
  `world_magnetic_model::uom`) so callers need not depend on them directly.
- Field docs on `Error::HeightOutsideOfValidityRange::{min_height_km, max_height_km}`.
- `#![warn(missing_docs)]`.

### Changed

- `thiserror` 2.0.11.

## [0.1.1] - 2025-01-12

[crates.io](https://crates.io/crates/world_magnetic_model/0.1.1)

### Added

- Examples: `demo_today`, `demo_wmm2020_wmm2025`, `demo_newtype`.

### Changed

- Published crate excludes `.git*`, `benches/`, `build.rs`, `ncei.noaa.gov/`,
  and `templates/`.
- Dropped the unused `uom` `f64` feature (SI `f32` only).
- `time` 0.3.37, `thiserror` 2.0.10, `rstest` 0.24, `serde` 1.0.217.

## [0.1.0] - 2025-01-01

[crates.io](https://crates.io/crates/world_magnetic_model/0.1.0)

Initial release.

- NOAA World Magnetic Model in Rust (`no_std`), WMM2020 and WMM2025
  coefficients compiled in (no runtime file I/O).
- `GeomagneticField::new(height, latitude, longitude, date)` with `uom` SI
  quantities and `time::Date`.
- Outputs: declination \(D\), inclination \(I\), \(X\), \(Y\), \(Z\), \(H\), \(F\),
  1σ uncertainties ([WMM2020 TR](https://doi.org/10.25923/ytk1-yx35) §3.4;
  [WMM2025 TR](https://doi.org/10.25923/prbc-s316) §3.4), blackout / caution
  zones (WMM2025 TR §1.8).
- Optional `serde` feature (derive on the internal struct).
- MIT OR Apache-2.0.
