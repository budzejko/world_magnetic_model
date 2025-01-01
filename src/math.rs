use libm::{powf, sqrtf};
use uom::num_traits::float::FloatCore;

/// Calculates the factorial of a number.
fn factorial(n: usize) -> f64 {
    (1..=n).fold(1.0, |acc, x| acc * x as f64)
}

/// Calculates the Legendre function (polynomial or associated function).
///
/// # Arguments
/// * `t` - The value to evaluate
/// * `n` - The degree of the polynomial
/// * `m` - The order of the polynomial
///
/// # Returns
/// * The value of the Legendre polynomial
///
/// # Reference
/// * Heiskanen and Moritz, 1967 (https://archive.org/details/HeiskanenMoritz1967PhysicalGeodesy/page/n33/mode/2up): equation 1-62
fn legendre_function(t: f32, n: usize, m: usize) -> f32 {
    let mut p = 0.0;
    for k in 0..=((n - m) / 2) {
        let num = (-1.0_f64).powi(k as i32)
            * factorial(2 * n - 2 * k)
            * (t as f64).powi(n as i32 - m as i32 - 2 * k as i32);
        let den = factorial(k) * factorial(n - k) * factorial(n - m - 2 * k);
        p += (num / den) as f32;
    }
    p * 2.0.powi(-(n as i32)) * powf((1.0 - t.powi(2)), (m as f32 / 2.0))
}

/// Calculates the Schmidt semi-normalised associated Legendre functions
/// for all degrees and orders up to 14.
///
/// # Arguments
/// * `mu` - The value to evaluate
///
/// # Reference
/// * Equation (5) in "The US/UK Workd Magnetic Model for 2020-2025"
pub(crate) fn schmidt_semi_normalised_associated_legendre(mu: f32) -> [f32; 104] {
    let mut psn = [0.0; 104];
    for n in 1..=13 {
        for m in 0..=n {
            let ix = index(n, m);
            if m == 0 {
                psn[ix] = legendre_function(mu, n, m);
            } else {
                psn[ix] = sqrtf((2.0 * factorial(n - m) / factorial(n + m)) as f32)
                    * legendre_function(mu, n, m);
            }
        }
    }
    psn
}

// Calculates index in flat array based on n and m matrix indices.
pub(crate) fn index(n: usize, m: usize) -> usize {
    debug_assert!(n >= 1);
    debug_assert!(n <= 13);
    debug_assert!(m <= n);
    (1..=n).sum::<usize>() + m - 1
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;
    use time::Month::*;

    #[rstest]
    #[case(1, 0, 0)]
    #[case(1, 1, 1)]
    #[case(2, 0, 2)]
    #[case(2, 1, 3)]
    #[case(2, 2, 4)]
    #[case(3, 0, 5)]
    #[case(8, 4, 39)]
    #[case(12, 0, 77)]
    #[case(12, 12, 89)]
    #[case(13, 0, 90)]
    #[case(13, 13, 103)]
    fn test_index(#[case] n: usize, #[case] m: usize, #[case] ix: usize) {
        assert_eq!(index(n, m), ix);
    }
}
