use crate::wmm::WMM_MODELS;
use time::Date;

const NUMBER_OF_MODELS: usize = 2;
const WMM_VALIDITY_RANGE_IN_YEARS: i32 = 5;

#[derive(Debug, PartialEq)]
pub(crate) struct WmmModel {
    pub(crate) model_version: i32,
    pub(crate) g_mfc: [f32; 90],
    pub(crate) h_mfc: [f32; 90],
    pub(crate) g_svc: [f32; 90],
    pub(crate) h_svc: [f32; 90],
}

#[derive(Clone, Debug, PartialEq)]
pub(crate) struct WmmErrorModel {
    pub(crate) model_version: i32,
    pub(crate) declination_constant_error_factor: f32,
    pub(crate) declination_variable_error_factor: f32,
    pub(crate) inclination_uncertainty: f32,
    pub(crate) h_uncertainty: f32,
    pub(crate) f_uncertainty: f32,
    pub(crate) x_uncertainty: f32,
    pub(crate) y_uncertainty: f32,
    pub(crate) z_uncertainty: f32,
}

const WMM_ERROR_MODELS: [WmmErrorModel; NUMBER_OF_MODELS] = [
    WmmErrorModel {
        model_version: 2020,
        declination_constant_error_factor: 0.26,
        declination_variable_error_factor: 5625.0,
        inclination_uncertainty: 0.21,
        f_uncertainty: 145.0,
        h_uncertainty: 128.0,
        x_uncertainty: 131.0,
        y_uncertainty: 94.0,
        z_uncertainty: 157.0,
    },
    WmmErrorModel {
        model_version: 2025,
        declination_constant_error_factor: 0.26,
        declination_variable_error_factor: 5417.0,
        inclination_uncertainty: 0.20,
        f_uncertainty: 138.0,
        h_uncertainty: 133.0,
        x_uncertainty: 137.0,
        y_uncertainty: 89.0,
        z_uncertainty: 141.0,
    },
];

/// Function to calculate a model version for a given date.
fn date_to_model_version(date: Date) -> i32 {
    date.year() / WMM_VALIDITY_RANGE_IN_YEARS * WMM_VALIDITY_RANGE_IN_YEARS
}

/// Function that returns WmmModel and WmmErrorModel for given date.
pub(crate) fn select_models(
    date: Date,
) -> Result<(&'static WmmModel, &'static WmmErrorModel), crate::error::Error> {
    let model_version = date_to_model_version(date);

    let wmm_model = WMM_MODELS
        .iter()
        .find(|model| model.model_version == model_version);
    let wmm_error_model = WMM_ERROR_MODELS
        .iter()
        .find(|model| model.model_version == model_version);

    match (wmm_model, wmm_error_model) {
        (Some(wmm_model), Some(wmm_error_model)) => Ok((wmm_model, wmm_error_model)),
        _ => Err(crate::error::Error::DateOutsideOfValidityRange),
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
    fn test_date_to_model_version(
        #[case] year_int: i32,
        #[case] month: time::Month,
        #[case] day: u8,
        #[case] model_version: i32,
    ) {
        assert_eq!(
            date_to_model_version(Date::from_calendar_date(year_int, month, day).unwrap()),
            model_version
        );
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
        #[case] model_version: i32,
    ) {
        assert!(
            select_models(Date::from_calendar_date(year_int, month, day).unwrap())
                .is_ok_and(|(wmm, wmm_error)| wmm.model_version == model_version
                    && wmm_error.model_version == model_version)
        )
    }
}
