use libm::{powf, sqrtf};
use uom::num_traits::float::FloatCore;

// Precomputed factorials up to 26!
const FACTORIAL: [f64; 27] = [
    1.0,
    1.0,
    2.0,
    6.0,
    24.0,
    120.0,
    720.0,
    5040.0,
    40320.0,
    362880.0,
    3628800.0,
    39916800.0,
    479001600.0,
    6227020800.0,
    87178291200.0,
    1307674368000.0,
    20922789888000.0,
    355687428096000.0,
    6402373705728000.0,
    121645100408832000.0,
    2432902008176640000.0,
    51090942171709440000.0,
    1124000727777607680000.0,
    25852016738884976640000.0,
    620448401733239439360000.0,
    15511210043330985984000000.0,
    403291461126605635584000000.0,
];

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
    debug_assert!(n >= 1);
    debug_assert!(n >= m);
    debug_assert!(n <= 13);
    let mut p = 0.0;
    for k in 0..=((n - m) / 2) {
        let num = (-1.0_f64).powi(k as i32)
            * FACTORIAL[2 * n - 2 * k]
            * (t as f64).powi(n as i32 - m as i32 - 2 * k as i32);
        let den = FACTORIAL[k] * FACTORIAL[n - k] * FACTORIAL[n - m - 2 * k];
        p += (num / den) as f32;
    }
    p * 2.0.powi(-(n as i32)) * powf(1.0 - t.powi(2), m as f32 / 2.0)
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
            let legendre = legendre_function(mu, n, m);
            psn[ix] = if m == 0 {
                legendre
            } else {
                let ratio = FACTORIAL[n - m] / FACTORIAL[n + m];
                sqrtf(2.0 * ratio as f32) * legendre
            };
        }
    }
    psn
}

/// Calculates index in flat array based on n and m matrix indices.
pub(crate) fn index(n: usize, m: usize) -> usize {
    debug_assert!(n >= 1);
    debug_assert!(n >= m);
    debug_assert!(n <= 13);
    n * (n + 1) / 2 + m - 1
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

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
