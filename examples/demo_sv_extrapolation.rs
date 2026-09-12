//! Compare full WMM declination with linear SV extrapolation.
//!
//! For each `demo_today` location, evaluate declination every half year from
//! 2025-01-01 through 2029-07-01 in three ways:
//!
//! 1. WMM at that date (epoch 2025 after 2025-01-01).
//! 2. Linear extrapolation from D and Ḋ on WMM2025 at 2025-01-01.
//! 3. Linear extrapolation from D and Ḋ on WMM2020 at 2024-12-31.
//!
//! Instantaneous Ḋ uses the official formula
//! `(ẎX − YẊ) / H²` with Ẋ, Ẏ recovered from two same-epoch evaluations
//! one day apart (X and Y are linear in time inside an epoch).

use core::f32::consts::PI;
use world_magnetic_model::GeomagneticField;
use world_magnetic_model::time::{Date, Month, util::days_in_year};
use world_magnetic_model::uom::si::angle::degree;
use world_magnetic_model::uom::si::f32::{Angle, Length};
use world_magnetic_model::uom::si::length::meter;
use world_magnetic_model::uom::si::magnetic_flux_density::nanotesla;

const HEIGHT_M: f32 = 100.0;

const POINTS: [(&str, f32, f32); 4] = [
    ("Świnoujście", 53.912668, 14.261214),
    ("Jakuszyce", 50.814168, 15.428267),
    ("Sejny", 54.109385, 23.346215),
    ("Wołosate", 49.066563, 22.68012),
];

fn main() {
    let height = Length::new::<meter>(HEIGHT_M);
    let ref_2025 = Date::from_calendar_date(2025, Month::January, 1).unwrap();
    let ref_2020 = Date::from_calendar_date(2024, Month::December, 31).unwrap();
    let dates = sample_dates();

    println!("Declination experiment: WMM vs linear SV extrapolation");
    println!("Height {HEIGHT_M} m above WGS 84. Interval 0.5 year from 2025-01-01.");
    println!("Dotted rates are instantaneous Ḋ [°/year] at the reference date.");
    println!();

    for (name, lat_deg, lon_deg) in POINTS {
        let lat = Angle::new::<degree>(lat_deg);
        let lon = Angle::new::<degree>(lon_deg);

        let d0_2025 = declination_deg(height, lat, lon, ref_2025);
        let sv_2025 = declination_rate_deg_per_year(height, lat, lon, ref_2025, next_day(ref_2025));
        let sigma_d = GeomagneticField::new(height, lat, lon, ref_2025)
            .unwrap()
            .declination_uncertainty()
            .get::<degree>();

        let d0_2020 = declination_deg(height, lat, lon, ref_2020);
        let sv_2020 = declination_rate_deg_per_year(height, lat, lon, ref_2020, prev_day(ref_2020));

        println!("{name}  ({lat_deg:.5}°N, {lon_deg:.5}°E)  σ_D(WMM2025) = {sigma_d:.2}°");
        println!("  WMM2025 @ {ref_2025}:  D = {d0_2025:7.4}°   SV = {sv_2025:+.4}°/year");
        println!("  WMM2020 @ {ref_2020}:  D = {d0_2020:7.4}°   SV = {sv_2020:+.4}°/year");
        println!();
        println!(
            "  {:10}  {:>8}  {:>10}  {:>10}  {:>10}  {:>10}",
            "date", "WMM", "ex.2025", "ex.2020", "Δ2025", "Δ2020"
        );

        for date in dates {
            let years_2025 = years_between(ref_2025, date);
            let years_2020 = years_between(ref_2020, date);
            let d_wmm = declination_deg(height, lat, lon, date);
            let d_ex_2025 = d0_2025 + sv_2025 * years_2025;
            let d_ex_2020 = d0_2020 + sv_2020 * years_2020;
            let delta_2025 = d_ex_2025 - d_wmm;
            let delta_2020 = d_ex_2020 - d_wmm;

            println!(
                "  {date}  {d_wmm:8.4}  {d_ex_2025:10.4}  {d_ex_2020:10.4}  {delta_2025:+10.4}  {delta_2020:+10.4}"
            );
        }
        println!();
    }

    println!("Columns in degrees. Δ = extrapolated - WMM.");
}

fn sample_dates() -> [Date; 10] {
    let mut dates = [Date::from_calendar_date(2025, Month::January, 1).unwrap(); 10];
    let mut i = 0;
    for year in 2025..=2029 {
        dates[i] = Date::from_calendar_date(year, Month::January, 1).unwrap();
        dates[i + 1] = Date::from_calendar_date(year, Month::July, 1).unwrap();
        i += 2;
    }
    dates
}

fn declination_deg(height: Length, lat: Angle, lon: Angle, date: Date) -> f32 {
    GeomagneticField::new(height, lat, lon, date)
        .unwrap()
        .declination()
        .get::<degree>()
}

/// Instantaneous Ḋ at `date` from X, Y at `date` and a same-epoch `probe`.
fn declination_rate_deg_per_year(
    height: Length,
    lat: Angle,
    lon: Angle,
    date: Date,
    probe: Date,
) -> f32 {
    let field = GeomagneticField::new(height, lat, lon, date).unwrap();
    let probed = GeomagneticField::new(height, lat, lon, probe).unwrap();
    let dt = years_between(date, probe);
    let x = field.x().get::<nanotesla>();
    let y = field.y().get::<nanotesla>();
    let x_dot = (probed.x().get::<nanotesla>() - x) / dt;
    let y_dot = (probed.y().get::<nanotesla>() - y) / dt;
    (y_dot * x - y * x_dot) / (x * x + y * y) * (180.0 / PI)
}

fn decimal_year(date: Date) -> f32 {
    (date.year() as f32) + (date.ordinal() - 1) as f32 / (days_in_year(date.year()) as f32)
}

fn years_between(from: Date, to: Date) -> f32 {
    decimal_year(to) - decimal_year(from)
}

fn next_day(date: Date) -> Date {
    date.next_day().expect("date has a following day")
}

fn prev_day(date: Date) -> Date {
    date.previous_day().expect("date has a previous day")
}
