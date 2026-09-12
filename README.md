# world_magnetic_model

This crate is a Rust implementation of the NOAA [World Magnetic Model (WMM)](https://www.ncei.noaa.gov/products/world-magnetic-model),
a mathematical representation of the Earth's core magnetic field and its temporal variations.

The crate's interface utilizes the [uom (Units of Measurement) crate](https://docs.rs/uom/latest/uom/) to represent physical quantities
accurately. WMM coefficient files are converted into code constants to eliminate the need for file reading at runtime.
This crate is compatible with `no_std` environments, meaning it does not depend on the Rust standard library and can be used in embedded,
bare-metal, or other restricted contexts by relying on the core crate instead.

## Bundled models

Each release ships two consecutive WMM coefficient sets (currently 2020 and 2025),
selected automatically from the evaluation date. A successor model is published while
the current epoch is still in force. Keeping both the outgoing and incoming sets in
the crate lets embedded firmware reserve a constant footprint for coefficients: an
update that replaces the older epoch with the newly published one does not grow the
image, and devices already carry the next model when the epoch boundary arrives.

Crate versions follow `api.epoch.patch`, where `epoch` is the year of the newest
bundled WMM model.

## Usage

Results are shown as \(D \pm \sigma_D\) (1σ). Field accuracy is the
[WMM error model](https://www.ncei.noaa.gov/products/world-magnetic-model/accuracy-limitations-error-model),
not the `f32` precision.

```rust
let geomagnetic_field = GeomagneticField::new(
    Length::new::<meter>(100.0), // Height above the WGS 84 ellipsoid
    Angle::new::<degree>(37.03), // WGS 84 latitude (negative values for the Southern Hemisphere)
    Angle::new::<degree>(-7.91), // WGS 84 longitude (negative values for the Western Hemisphere)
    Date::from_ordinal_date(2029, 15)? // Date (15th day of 2029)
)?;

let d = geomagnetic_field.declination().get::<degree>();
let sigma = geomagnetic_field.declination_uncertainty().get::<degree>();

assert_eq!(format!("{d:.2} ± {sigma:.2}°"), "-0.17 ± 0.33°");
```

## Cargo features

Enable **`serde`** for `Serialize`/`Deserialize` on [`GeomagneticField`] and
[`WarningZone`]. `GeomagneticField` encodes the constructor query plus public
results (declination, intensity, warning, …). Deserialization uses the
inputs only and recomputes through [`GeomagneticField::new`]. Enabling
`serde` does not pull in `std`. Encode with any serde codec, for example
`serde_json` or `postcard`.

```rust
let json = serde_json::to_string(&geomagnetic_field)?;
let from_json: GeomagneticField = serde_json::from_str(&json)?;

let mut buf = [0u8; 256];
let bytes = postcard::to_slice(&geomagnetic_field, &mut buf)?;
let from_postcard: GeomagneticField = postcard::from_bytes(bytes)?;

assert_eq!(geomagnetic_field.declination(), from_json.declination());
assert_eq!(geomagnetic_field.declination(), from_postcard.declination());
```

## World Magnetic Model

The [World Magnetic Model](https://www.ncei.noaa.gov/products/world-magnetic-model) is the standard model used by
the U.S. Department of Defense, the U.K. Ministry of Defence, the North Atlantic Treaty Organization (NATO)
and the International Hydrographic Organization (IHO), for navigation, attitude and heading referencing systems
using the geomagnetic field. It is also used widely in civilian navigation and heading systems.
The model is produced at 5-year intervals, with the current model expiring on December 31, 2029. The current
model WMM2025 is produced jointly by the NCEI and the British Geological Survey (BGS). The model, associated
software, and documentation are distributed by NCEI on behalf of US National Geospatial-Intelligence Agency
and by BGS on behalf of UK Defence Geographic Centre.
[\[source\]](https://www.ncei.noaa.gov/metadata/geoportal/rest/metadata/item/gov.noaa.ngdc:WMM2025/html)

Please refer to [World Magnetic Model Accuracy, Limitations, and Error Model](https://www.ncei.noaa.gov/products/world-magnetic-model/accuracy-limitations-error-model).

## NOAA License Statement
The WMM source code is in the public domain and not licensed or under copyright. The information and software
may be used freely by the public. As required by 17 U.S.C. 403, third parties producing copyrighted works
consisting predominantly of the material produced by U.S. government agencies must provide notice with such
work(s) identifying the U.S. Government material incorporated and stating that such material is not subject
to copyright protection.

The WMM model and associated data files are produced by the U.S. Government and are not subject to copyright.

## Credit
This work was inspired by [geomag-wmm](https://crates.io/crates/geomag-wmm/0.1.0).

## License
Licensed under either of [Apache License, Version 2.0](https://github.com/budzejko/world_magnetic_model/blob/main/LICENSE-APACHE)
or [MIT license](https://github.com/budzejko/world_magnetic_model/blob/main/LICENSE-MIT) at your option.

License: MIT OR Apache-2.0
