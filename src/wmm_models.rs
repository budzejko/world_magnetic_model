use crate::synthesis::N_COEFF;
use crate::wmm_data::WMM_MODELS;
use time::Date;

const NUMBER_OF_MODELS: usize = WMM_MODELS.len();
/// WMM release interval and nominal validity span (years).
const WMM_EPOCH_INTERVAL_YEARS: i32 = 5;

#[derive(Debug, PartialEq)]
pub(crate) struct WmmModel {
    /// Epoch year \(t_0\) of this coefficient set.
    pub(crate) epoch_year: i32,
    /// Gauss coefficients \(g_n^m\) [nT].
    pub(crate) g: [f32; N_COEFF],
    /// Gauss coefficients \(h_n^m\) [nT] (not horizontal intensity \(H\)).
    pub(crate) h: [f32; N_COEFF],
    /// Secular variation \(\dot g_n^m\) [nT/year].
    pub(crate) g_dot: [f32; N_COEFF],
    /// Secular variation \(\dot h_n^m\) [nT/year].
    pub(crate) h_dot: [f32; N_COEFF],
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub(crate) struct WmmErrorModel {
    pub(crate) epoch_year: i32,
    /// \(\sigma_c\) (degrees) in the [WMM2020 TR](https://doi.org/10.25923/ytk1-yx35) §3.4; [WMM2025 TR](https://doi.org/10.25923/prbc-s316) §3.4 declination error model.
    pub(crate) declination_sigma_c: f32,
    /// \(\sigma_v\) (nT·deg) in the WMM2020 TR §3.4; WMM2025 TR §3.4 declination error model.
    pub(crate) declination_sigma_v: f32,
    pub(crate) inclination_uncertainty: f32,
    /// Global RMS uncertainty of horizontal intensity \(H\) [nT] (not Gauss \(h_n^m\)).
    pub(crate) h_uncertainty: f32,
    pub(crate) f_uncertainty: f32,
    pub(crate) x_uncertainty: f32,
    pub(crate) y_uncertainty: f32,
    pub(crate) z_uncertainty: f32,
}

const WMM_ERROR_MODELS: [WmmErrorModel; NUMBER_OF_MODELS] = [
    WmmErrorModel {
        epoch_year: 2020,
        declination_sigma_c: 0.26,
        declination_sigma_v: 5625.0,
        inclination_uncertainty: 0.21,
        f_uncertainty: 145.0,
        h_uncertainty: 128.0,
        x_uncertainty: 131.0,
        y_uncertainty: 94.0,
        z_uncertainty: 157.0,
    },
    WmmErrorModel {
        epoch_year: 2025,
        declination_sigma_c: 0.26,
        declination_sigma_v: 5417.0,
        inclination_uncertainty: 0.20,
        f_uncertainty: 138.0,
        h_uncertainty: 133.0,
        x_uncertainty: 137.0,
        y_uncertainty: 89.0,
        z_uncertainty: 141.0,
    },
];

/// Five-year epoch bucket containing `date` (e.g. 2024 → 2020, 2025 → 2025).
/// The bucket may not match a bundled WMM (2019 → 2015, 2030 → 2030).
fn date_to_epoch_year(date: Date) -> i32 {
    date.year() / WMM_EPOCH_INTERVAL_YEARS * WMM_EPOCH_INTERVAL_YEARS
}

const _: () = assert!(NUMBER_OF_MODELS > 0, "at least one WMM model is bundled");

/// Inclusive validity window of the bundled WMM epochs (first day of the
/// earliest epoch through the last day of the latest epoch).
pub(crate) const BUNDLED_MODELS_VALIDITY_RANGE: (Date, Date) = {
    let mut min_epoch = WMM_MODELS[0].epoch_year;
    let mut max_epoch = WMM_MODELS[0].epoch_year;
    let mut i = 1;
    while i < NUMBER_OF_MODELS {
        let year = WMM_MODELS[i].epoch_year;
        if year < min_epoch {
            min_epoch = year;
        }
        if year > max_epoch {
            max_epoch = year;
        }
        i += 1;
    }

    (
        match Date::from_calendar_date(min_epoch, time::Month::January, 1) {
            Ok(date) => date,
            Err(_) => panic!("epoch year has a January 1"),
        },
        match Date::from_calendar_date(
            max_epoch + WMM_EPOCH_INTERVAL_YEARS - 1,
            time::Month::December,
            31,
        ) {
            Ok(date) => date,
            Err(_) => panic!("last validity year has a December 31"),
        },
    )
};

/// WMM2020 TR §3.4; WMM2025 TR §3.4 error model for a bundled epoch year, if present.
pub(crate) fn error_model_by_epoch_year(epoch_year: i32) -> Option<&'static WmmErrorModel> {
    WMM_ERROR_MODELS
        .iter()
        .find(|model| model.epoch_year == epoch_year)
}

/// Main-field and WMM2020 TR §3.4; WMM2025 TR §3.4 error models for `date`, or
/// [`DateOutsideOfValidityRange`](crate::error::Error::DateOutsideOfValidityRange).
pub(crate) fn select_models(
    date: Date,
) -> Result<(&'static WmmModel, &'static WmmErrorModel), crate::error::Error> {
    let epoch_year = date_to_epoch_year(date);

    let wmm_model = WMM_MODELS
        .iter()
        .find(|model| model.epoch_year == epoch_year);
    let wmm_error_model = error_model_by_epoch_year(epoch_year);

    match (wmm_model, wmm_error_model) {
        (Some(wmm_model), Some(wmm_error_model)) => Ok((wmm_model, wmm_error_model)),
        _ => {
            let (min_date, max_date) = BUNDLED_MODELS_VALIDITY_RANGE;
            Err(crate::error::Error::DateOutsideOfValidityRange { min_date, max_date })
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use time::Month::*;

    #[rstest]
    #[case(2019, January, 1, 2015)]
    #[case(2020, January, 1, 2020)]
    #[case(2021, January, 1, 2020)]
    #[case(2022, January, 1, 2020)]
    #[case(2023, January, 1, 2020)]
    #[case(2024, January, 1, 2020)]
    #[case(2024, January, 2, 2020)]
    #[case(2024, December, 12, 2020)]
    #[case(2024, December, 31, 2020)]
    #[case(2025, January, 1, 2025)]
    #[case(2026, January, 1, 2025)]
    #[case(2027, January, 1, 2025)]
    #[case(2028, January, 1, 2025)]
    #[case(2029, January, 1, 2025)]
    #[case(2030, January, 1, 2030)]
    #[case(2031, January, 1, 2030)]
    fn test_date_to_epoch_year(
        #[case] year_int: i32,
        #[case] month: time::Month,
        #[case] day: u8,
        #[case] epoch_year: i32,
    ) {
        assert_eq!(
            date_to_epoch_year(Date::from_calendar_date(year_int, month, day).unwrap()),
            epoch_year
        );
    }

    #[test]
    fn test_bundled_models_validity_range() {
        assert_eq!(
            BUNDLED_MODELS_VALIDITY_RANGE,
            (
                Date::from_calendar_date(2020, January, 1).unwrap(),
                Date::from_calendar_date(2029, December, 31).unwrap(),
            )
        );
    }

    #[test]
    fn test_error_model_by_epoch_year() {
        assert_eq!(
            error_model_by_epoch_year(2020).map(|model| model.epoch_year),
            Some(2020)
        );
        assert_eq!(
            error_model_by_epoch_year(2025).map(|model| model.epoch_year),
            Some(2025)
        );
        assert!(error_model_by_epoch_year(2015).is_none());
        assert!(error_model_by_epoch_year(2030).is_none());
    }

    #[rstest]
    #[case(2020, January, 1, 2020)]
    #[case(2021, January, 1, 2020)]
    #[case(2022, January, 1, 2020)]
    #[case(2023, January, 1, 2020)]
    #[case(2024, January, 1, 2020)]
    #[case(2024, January, 2, 2020)]
    #[case(2024, December, 12, 2020)]
    #[case(2024, December, 31, 2020)]
    #[case(2025, January, 1, 2025)]
    #[case(2026, January, 1, 2025)]
    #[case(2027, January, 1, 2025)]
    #[case(2028, January, 1, 2025)]
    #[case(2029, January, 1, 2025)]
    #[case(2029, December, 31, 2025)]
    fn test_select_models(
        #[case] year_int: i32,
        #[case] month: time::Month,
        #[case] day: u8,
        #[case] epoch_year: i32,
    ) {
        assert!(
            select_models(Date::from_calendar_date(year_int, month, day).unwrap())
                .is_ok_and(|(wmm, wmm_error)| wmm.epoch_year == epoch_year
                    && wmm_error.epoch_year == epoch_year)
        )
    }
}
