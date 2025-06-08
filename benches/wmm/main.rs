use criterion::{criterion_group, criterion_main, Criterion};
use std::hint::black_box;
use world_magnetic_model::time::Date;
use world_magnetic_model::uom::si::{
    angle::degree,
    f32::{Angle, Length},
    length::foot,
};
use world_magnetic_model::GeomagneticField;

fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("wwm", |b| {
        b.iter(|| {
            GeomagneticField::new(
                black_box(Length::new::<foot>(2000.0)),
                black_box(Angle::new::<degree>(54.0)),
                black_box(Angle::new::<degree>(18.0)),
                black_box(Date::from_ordinal_date(2024, 183).unwrap()),
            )
        })
    });
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
