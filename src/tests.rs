//! Unit tests that need crate-private items (`*_rate`, geodetic helpers, …).
#![allow(clippy::too_many_arguments)]

extern crate std;

use time::Date;

mod api;
mod geodetic;
mod official;
mod traits;

/// Inverse of [`crate::years_since_epoch`]: keep the year offset and day fraction
/// separate so adding `epoch_year` in `f32` does not undo the precision.
fn years_since_epoch_to_date(epoch_year: i32, years_since: f32) -> Date {
    let year_int = epoch_year + years_since.trunc() as i32;
    let day_ordinal = f32::round(years_since.fract() * time::util::days_in_year(year_int) as f32);
    Date::from_ordinal_date(year_int, day_ordinal as u16 + 1).unwrap()
}

fn decimal_year_to_date(decimal_year: f32) -> Date {
    years_since_epoch_to_date(decimal_year.trunc() as i32, decimal_year.fract())
}
