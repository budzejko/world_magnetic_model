use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use world_magnetic_model::GeomagneticField;
use world_magnetic_model::time::Date;
use world_magnetic_model::uom::si::{
    angle::degree,
    f32::{Angle, Length},
    length::foot,
};

fn geomagnetic_field(latitude_deg: f32) -> Result<GeomagneticField, world_magnetic_model::Error> {
    GeomagneticField::new(
        black_box(Length::new::<foot>(2000.0)),
        black_box(Angle::new::<degree>(latitude_deg)),
        black_box(Angle::new::<degree>(18.0)),
        black_box(Date::from_ordinal_date(2024, 183).unwrap()),
    )
}

fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("wmm");
    group.bench_function("equator", |b| b.iter(|| geomagnetic_field(black_box(0.0))));
    group.bench_function("mid_latitude", |b| {
        b.iter(|| geomagnetic_field(black_box(54.0)))
    });
    group.bench_function("geographic_pole", |b| {
        b.iter(|| geomagnetic_field(black_box(90.0)))
    });
    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
