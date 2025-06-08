# world_magnetic_model

This crate is a Rust implementation of the NOAA [World Magnetic Model (WMM)](https://www.ncei.noaa.gov/products/world-magnetic-model),
a mathematical representation of the Earth's core magnetic field and its temporal variations.

The crate's interface utilizes the [uom (Units of Measurement) crate](https://docs.rs/uom/latest/uom/) to represent physical quantities
accurately. WMM coefficient files are converted into code constants to eliminate the need for file reading at runtime.
This crate is compatible with no_std environments, meaning it does not depend on the Rust standard library and can be used in embedded,
bare-metal, or other restricted contexts by relying on the core crate instead. The implemented models include WMM2020 and WMM2025.

## Usage
```rust
let geomagnetic_field = GeomagneticField::new(
    Length::new::<meter>(100.0), // Height above the WGS 84 ellipsoid
    Angle::new::<degree>(37.03), // WGS 84 latitude (negative values for the Southern Hemisphere)
    Angle::new::<degree>(-7.91), // WGS 84 longitude (negative values for the Western Hemisphere)
    Date::from_ordinal_date(2029, 15)? // Date (15th day of 2029)
)?;

assert_eq!(
    geomagnetic_field.declination().get::<degree>(),
    -0.17367662
);
assert_eq!(
    geomagnetic_field.declination_uncertainty().get::<degree>(),
    0.32549456
);
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
This work was inspired from [geomag-wmm](https://git.hostux.fr/ConstellationVFR/geomag-wmm).

## License
Licensed under either of [Apache License, Version 2.0](https://github.com/budzejko/world_magnetic_model/blob/main/LICENSE-APACHE)
or [MIT license](https://github.com/budzejko/world_magnetic_model/blob/main/LICENSE-MIT) at your option.

License: MIT OR Apache-2.0
