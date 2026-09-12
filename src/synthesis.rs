#[cfg(test)]
use libm::sqrtf;
use libm::{cosf, sinf};

const N_MAX: usize = 12;
/// Length of associated-Legendre tables indexed \(n(n+1)/2+m\) for
/// \(0 \le m \le n \le n_{\max}\) (includes \(P_0^0\)). Same triangular layout as
/// [NOAA WMM C](https://www.ncei.noaa.gov/products/world-magnetic-model) `MAG_PcupLow`; only two degree rows of \(P\) are kept at runtime.
const N_LEGENDRE: usize = (N_MAX + 1) * (N_MAX + 2) / 2; // 91
/// Number of WMM Gauss coefficients \(n=1..n_{\max}\), \(m=0..n\).
pub(crate) const N_COEFF: usize = N_MAX * (N_MAX + 3) / 2; // 90
/// Array length for azimuthal order \(m = 0..n_{\max}\).
const M_LEN: usize = N_MAX + 1; // 13

/// Flat index \(n(n+1)/2 + m\) of an associated Legendre table that includes
/// \(P_0^0\) at 0 (same layout as NOAA WMM C `MAG_PcupLow`).
#[inline]
fn legendre_index(n: usize, m: usize) -> usize {
    n * (n + 1) / 2 + m
}

/// Schmidt semi-normalization factors \(s_n^m\) relative to Gauss-normalized
/// associated Legendre functions ([WMM2025 TR](https://doi.org/10.25923/prbc-s316) §1.2).
/// NOAA WMM C `MAG_PcupLow` uses the same recurrence and calls it Schmidt
/// quasi-normalization. Independent of \(\theta\). Indexed by `legendre_index`,
/// including \(s_0^0\).
///
/// \(s_n^0 = s_{n-1}^0 (2n-1)/n\),
/// \(s_n^m = s_n^{m-1}\sqrt{(n-m+1)\,(2\text{ if }m=1\text{ else }1)/(n+m)}\).
const SCHMIDT_SEMI_NORM: [f32; N_LEGENDRE] = [
    f32::from_bits(0x3f800000),
    f32::from_bits(0x3f800000),
    f32::from_bits(0x3f800000),
    f32::from_bits(0x3fc00000),
    f32::from_bits(0x3fddb3d7),
    f32::from_bits(0x3f5db3d7),
    f32::from_bits(0x40200000),
    f32::from_bits(0x4043f58d),
    f32::from_bits(0x3ff7def6),
    f32::from_bits(0x3f4a62c3),
    f32::from_bits(0x408c0000),
    f32::from_bits(0x40b1166a),
    f32::from_bits(0x407a708b),
    f32::from_bits(0x4005dd98),
    f32::from_bits(0x3f3d5086),
    f32::from_bits(0x40fc0000),
    f32::from_bits(0x4122aa51),
    f32::from_bits(0x40f5ed43),
    f32::from_bits(0x4096994b),
    f32::from_bits(0x400dfc65),
    f32::from_bits(0x3f33997c),
    f32::from_bits(0x41670000),
    f32::from_bits(0x41973999),
    f32::from_bits(0x416f1b93),
    f32::from_bits(0x411f67b8),
    f32::from_bits(0x40ae9e9e),
    f32::from_bits(0x4014ea85),
    f32::from_bits(0x3f2bf418),
    f32::from_bits(0x41d68000),
    f32::from_bits(0x420de0df),
    f32::from_bits(0x41e7afbc),
    f32::from_bits(0x41a3d3bb),
    f32::from_bits(0x41459538),
    f32::from_bits(0x40c59538),
    f32::from_bits(0x401aff2c),
    f32::from_bits(0x3f25b2d1),
    f32::from_bits(0x42491800),
    f32::from_bits(0x42861000),
    f32::from_bits(0x42605458),
    f32::from_bits(0x4225ada4),
    f32::from_bits(0x41d5e3c6),
    f32::from_bits(0x416d4a14),
    f32::from_bits(0x40dbaff1),
    f32::from_bits(0x40206fd9),
    f32::from_bits(0x3f206fd9),
    f32::from_bits(0x42bdec00),
    f32::from_bits(0x42fece92),
    f32::from_bits(0x42d94cd1),
    f32::from_bits(0x42a5f736),
    f32::from_bits(0x426180c1),
    f32::from_bits(0x4206c387),
    f32::from_bits(0x418b2ef5),
    f32::from_bits(0x40f112a0),
    f32::from_bits(0x40255fe2),
    f32::from_bits(0x3f1beaa7),
    f32::from_bits(0x43346d00),
    f32::from_bits(0x4373493d),
    f32::from_bits(0x4352b122),
    f32::from_bits(0x432547c5),
    f32::from_bits(0x42e9bde1),
    f32::from_bits(0x4293d4cc),
    f32::from_bits(0x422547c5),
    f32::from_bits(0x41a05872),
    f32::from_bits(0x4102ebeb),
    f32::from_bits(0x4029e7ff),
    f32::from_bits(0x3f17f800),
    f32::from_bits(0x43ac3980),
    f32::from_bits(0x43e93177),
    f32::from_bits(0x43cc8624),
    f32::from_bits(0x43a3fbe8),
    f32::from_bits(0x436f8395),
    f32::from_bits(0x431e6c72),
    f32::from_bits(0x42bc3c3d),
    f32::from_bits(0x42466add),
    f32::from_bits(0x41b61490),
    f32::from_bits(0x410d09f0),
    f32::from_bits(0x402e1a24),
    f32::from_bits(0x3f14798b),
    f32::from_bits(0x44250c70),
    f32::from_bits(0x446041c2),
    f32::from_bits(0x4446c84f),
    f32::from_bits(0x44224e22),
    f32::from_bits(0x43f37533),
    f32::from_bits(0x43a702bc),
    f32::from_bits(0x43504c85),
    f32::from_bits(0x42ea1b95),
    f32::from_bits(0x426a1b95),
    f32::from_bits(0x41cc5893),
    f32::from_bits(0x4116eb65),
    f32::from_bits(0x403203d3),
    f32::from_bits(0x3f11593d),
];

/// Recurrence \(k_n^m = ((n-1)^2-m^2)/((2n-1)(2n-3))\). `legendre_index` layout;
/// 0 where unused.
///
/// Used only for \(n>1\) and \(m\le n-2\). Independent of \(\theta\).
const LEGENDRE_K: [f32; N_LEGENDRE] = [
    f32::from_bits(0x00000000),
    f32::from_bits(0x00000000),
    f32::from_bits(0x00000000),
    f32::from_bits(0x3eaaaaab),
    f32::from_bits(0x00000000),
    f32::from_bits(0x00000000),
    f32::from_bits(0x3e888889),
    f32::from_bits(0x3e4ccccd),
    f32::from_bits(0x00000000),
    f32::from_bits(0x00000000),
    f32::from_bits(0x3e83a83b),
    f32::from_bits(0x3e6a0ea1),
    f32::from_bits(0x3e124925),
    f32::from_bits(0x00000000),
    f32::from_bits(0x00000000),
    f32::from_bits(0x3e820821),
    f32::from_bits(0x3e73cf3d),
    f32::from_bits(0x3e430c31),
    f32::from_bits(0x3de38e39),
    f32::from_bits(0x00000000),
    f32::from_bits(0x00000000),
    f32::from_bits(0x3e814afd),
    f32::from_bits(0x3e783e10),
    f32::from_bits(0x3e59364e),
    f32::from_bits(0x3e257eb5),
    f32::from_bits(0x3dba2e8c),
    f32::from_bits(0x00000000),
    f32::from_bits(0x00000000),
    f32::from_bits(0x3e80e526),
    f32::from_bits(0x3e7aa11e),
    f32::from_bits(0x3e652598),
    f32::from_bits(0x3e4157b8),
    f32::from_bits(0x3e0f377f),
    f32::from_bits(0x3d9d89d9),
    f32::from_bits(0x00000000),
    f32::from_bits(0x00000000),
    f32::from_bits(0x3e80a80b),
    f32::from_bits(0x3e7c0fc1),
    f32::from_bits(0x3e6c4ec5),
    f32::from_bits(0x3e520d21),
    f32::from_bits(0x3e2d4ad5),
    f32::from_bits(0x3dfc0fc1),
    f32::from_bits(0x3d888889),
    f32::from_bits(0x00000000),
    f32::from_bits(0x00000000),
    f32::from_bits(0x3e808081),
    f32::from_bits(0x3e7cfcfd),
    f32::from_bits(0x3e70f0f1),
    f32::from_bits(0x3e5cdcdd),
    f32::from_bits(0x3e40c0c1),
    f32::from_bits(0x3e1c9c9d),
    f32::from_bits(0x3de0e0e1),
    f32::from_bits(0x3d70f0f1),
    f32::from_bits(0x00000000),
    f32::from_bits(0x00000000),
    f32::from_bits(0x3e806573),
    f32::from_bits(0x3e7d9f4e),
    f32::from_bits(0x3e741c88),
    f32::from_bits(0x3e644293),
    f32::from_bits(0x3e4e1170),
    f32::from_bits(0x3e31891d),
    f32::from_bits(0x3e0ea99c),
    f32::from_bits(0x3dcae5d8),
    f32::from_bits(0x3d579436),
    f32::from_bits(0x00000000),
    f32::from_bits(0x00000000),
    f32::from_bits(0x3e805220),
    f32::from_bits(0x3e7e1340),
    f32::from_bits(0x3e76603e),
    f32::from_bits(0x3e698b3a),
    f32::from_bits(0x3e579436),
    f32::from_bits(0x3e407b30),
    f32::from_bits(0x3e244029),
    f32::from_bits(0x3e02e321),
    f32::from_bits(0x3db8c82e),
    f32::from_bits(0x3d430c31),
    f32::from_bits(0x00000000),
    f32::from_bits(0x00000000),
    f32::from_bits(0x3e8043d8),
    f32::from_bits(0x3e7e68f2),
    f32::from_bits(0x3e780cb8),
    f32::from_bits(0x3e6d7304),
    f32::from_bits(0x3e5e9bd3),
    f32::from_bits(0x3e4b8728),
    f32::from_bits(0x3e343501),
    f32::from_bits(0x3e18a55e),
    f32::from_bits(0x3df1b07f),
    f32::from_bits(0x3da99b4c),
    f32::from_bits(0x3d321643),
    f32::from_bits(0x00000000),
    f32::from_bits(0x00000000),
];

/// Recurrence coefficient \(k_n^m\) from the `LEGENDRE_K` table.
#[inline]
fn legendre_k(n: usize, m: usize) -> f32 {
    LEGENDRE_K[legendre_index(n, m)]
}

/// After degree `n`: current → prev1, prev1 → prev2, prev2 → scratch.
#[inline]
fn rotate_row_indices(i_prev2: &mut usize, i_prev1: &mut usize, i_cur: &mut usize) {
    let next = *i_prev2;
    *i_prev2 = *i_prev1;
    *i_prev1 = *i_cur;
    *i_cur = next;
}

/// Gauss-normalized \(P_n^m\) and \(\partial P_n^m/\partial\theta\) from rows
/// \(n-1\) and \(n-2\), where \(\theta = \pi/2 - \varphi'\) is geocentric colatitude.
/// Row index is `m`. Seed: \(P_0^0 = 1\) lives in the starting prev1 slot.
///
/// Synthesis applies \(\partial\tilde{P}_n^m/\partial\varphi' = -s\,\partial P_n^m/\partial\theta\).
#[inline]
#[allow(clippy::too_many_arguments)]
fn gauss_norm_p_dp(
    n: usize,
    m: usize,
    cos_theta: f32,
    sin_theta: f32,
    p_prev1: &[f32; M_LEN],
    p_prev2: &[f32; M_LEN],
    dp_prev1: &[f32; M_LEN],
    dp_prev2: &[f32; M_LEN],
) -> (f32, f32) {
    if n == m {
        (
            sin_theta * p_prev1[m - 1],
            sin_theta * dp_prev1[m - 1] + cos_theta * p_prev1[m - 1],
        )
    } else if n == 1 && m == 0 {
        (cos_theta, -sin_theta)
    } else if m > n - 2 {
        (
            cos_theta * p_prev1[m],
            cos_theta * dp_prev1[m] - sin_theta * p_prev1[m],
        )
    } else {
        let k = legendre_k(n, m);
        (
            cos_theta * p_prev1[m] - k * p_prev2[m],
            cos_theta * dp_prev1[m] - sin_theta * p_prev1[m] - k * dp_prev2[m],
        )
    }
}

/// Time-linear Gauss coefficient (WMM2025 TR §1.2):
/// \(g_n^m(t) = g_n^m + \dot g_n^m\,(t-t_0)\).
#[inline]
pub(crate) fn coeff_at_t(main: f32, secular: f32, time_delta: f32) -> f32 {
    main + time_delta * secular
}

/// Two-row associated Legendre recurrence: Schmidt factors are applied as each
/// \(\tilde{P}_n^m\) is formed, then summed with Gauss coefficients at epoch offset.
///
/// Returns \((X'_\varphi,\, Y',\, Z'_\varphi)\). \(Y'\) matches WMM2025 TR \(Y' = Y\sin\theta\).
/// \(X'_\varphi\) and \(Z'_\varphi\) use \(\partial\tilde{P}/\partial\varphi'\) and
/// \(+(n+1)\tilde{P}\), so they equal \(-X'\) and \(-Z'\) of WMM2025 TR §1.2 (which
/// differentiate with respect to colatitude \(\theta\)). The geodetic rotation after synthesis
/// restores the WMM2025 TR NED frame.
///
/// `ratio` is \(a/r\). `sin_m` / `cos_m` are `sin(mλ)`, `cos(mλ)` for `m = 0..=12`.
#[allow(clippy::too_many_arguments)]
pub(crate) fn synthesize_xyz_prime(
    cos_theta: f32,
    sin_theta: f32,
    sin_m: &[f32; M_LEN],
    cos_m: &[f32; M_LEN],
    ratio: f32,
    g: &[f32; N_COEFF],
    h: &[f32; N_COEFF],
    g_dot: &[f32; N_COEFF],
    h_dot: &[f32; N_COEFF],
    time_delta: f32,
) -> (f32, f32, f32) {
    // Slots: prev2, prev1, current. P_0^0 in starting prev1. Rotate, do not copy rows.
    let mut p = [[0.0f32; M_LEN]; 3];
    let mut dp = [[0.0f32; M_LEN]; 3];
    p[1][0] = 1.0;
    let mut i_prev2 = 0usize;
    let mut i_prev1 = 1usize;
    let mut i_cur = 2usize;

    // \((a/r)^{n+2}\); starts at \(n = 1\) so \((a/r)^3\).
    let mut a_over_r_n2 = ratio * ratio * ratio;
    let mut x_prime = 0.0;
    let mut y_prime = 0.0;
    let mut z_prime = 0.0;

    for n in 1..=N_MAX {
        let mut x_tmp = 0.0;
        let mut y_tmp = 0.0;
        let mut z_tmp = 0.0;

        for m in 0..=n {
            let (p_nm, dp_nm) = gauss_norm_p_dp(
                n,
                m,
                cos_theta,
                sin_theta,
                &p[i_prev1],
                &p[i_prev2],
                &dp[i_prev1],
                &dp[i_prev2],
            );
            p[i_cur][m] = p_nm;
            dp[i_cur][m] = dp_nm;

            let s = SCHMIDT_SEMI_NORM[legendre_index(n, m)];
            let psn = s * p_nm;
            let dpsn = -s * dp_nm;
            let ix = coeff_index(n, m);
            let g_t = coeff_at_t(g[ix], g_dot[ix], time_delta);
            let h_t = coeff_at_t(h[ix], h_dot[ix], time_delta);
            let g_c_h_s = g_t * cos_m[m] + h_t * sin_m[m];
            let g_s_h_c = g_t * sin_m[m] - h_t * cos_m[m];
            x_tmp += g_c_h_s * dpsn;
            y_tmp += m as f32 * g_s_h_c * psn;
            z_tmp += g_c_h_s * psn;
        }

        x_prime += a_over_r_n2 * x_tmp;
        y_prime += a_over_r_n2 * y_tmp;
        z_prime += (n as f32 + 1.0) * a_over_r_n2 * z_tmp;
        a_over_r_n2 *= ratio;
        rotate_row_indices(&mut i_prev2, &mut i_prev1, &mut i_cur);
    }
    (x_prime, y_prime, z_prime)
}

/// Near-pole path: φ'-based \(X'\), \(Z'\) (same sign convention as
/// `synthesize_xyz_prime`) and east \(Y\) rather than \(Y'\). Only \(m \in \{0,1\}\).
///
/// Used when \(|\sin\theta|\) is small, not only at \(\sin\theta = 0\). \(Y\) is the east
/// component (WMM2025 TR §1.4 / NOAA WMM C `MAG_SummationSpecial`): the recurrence for
/// \(Q_n = P_n^1 / \sin\theta\) stays finite as \(\sin\theta \to 0\). Terms \(m \ge 2\) vanish
/// as \(\sin^m\theta\) and are omitted (no denormal \(P_n^m\)).
#[allow(clippy::too_many_arguments)]
pub(crate) fn synthesize_near_pole(
    cos_theta: f32,
    sin_theta: f32,
    lambda: f32,
    ratio: f32,
    g: &[f32; N_COEFF],
    h: &[f32; N_COEFF],
    g_dot: &[f32; N_COEFF],
    h_dot: &[f32; N_COEFF],
    time_delta: f32,
) -> (f32, f32, f32) {
    let mut p = [[0.0f32; M_LEN]; 3];
    let mut dp = [[0.0f32; M_LEN]; 3];
    p[1][0] = 1.0;
    let mut i_prev2 = 0usize;
    let mut i_prev1 = 1usize;
    let mut i_cur = 2usize;

    let sin_l = sinf(lambda);
    let cos_l = cosf(lambda);

    let mut a_over_r_n2 = ratio * ratio * ratio;
    let mut x_prime = 0.0;
    let mut y = 0.0;
    let mut z_prime = 0.0;
    // \(Q_0 = 1\), \(Q_1 = P_1^1 / \sin\theta = 1\)
    let mut q_prev2 = 1.0f32;
    let mut q_prev1 = 1.0f32;

    for n in 1..=N_MAX {
        let mut x_tmp = 0.0;
        let mut z_tmp = 0.0;

        let q = if n == 1 {
            1.0
        } else {
            cos_theta * q_prev1 - legendre_k(n, 1) * q_prev2
        };

        for m in 0..=1 {
            let (p_nm, dp_nm) = gauss_norm_p_dp(
                n,
                m,
                cos_theta,
                sin_theta,
                &p[i_prev1],
                &p[i_prev2],
                &dp[i_prev1],
                &dp[i_prev2],
            );
            p[i_cur][m] = p_nm;
            dp[i_cur][m] = dp_nm;

            let s = SCHMIDT_SEMI_NORM[legendre_index(n, m)];
            let psn = s * p_nm;
            let dpsn = -s * dp_nm;
            let ix = coeff_index(n, m);
            let g_t = coeff_at_t(g[ix], g_dot[ix], time_delta);
            let h_t = coeff_at_t(h[ix], h_dot[ix], time_delta);
            let (sin_m, cos_m) = if m == 0 { (0.0, 1.0) } else { (sin_l, cos_l) };
            let g_c_h_s = g_t * cos_m + h_t * sin_m;
            x_tmp += g_c_h_s * dpsn;
            z_tmp += g_c_h_s * psn;
            if m == 1 {
                y += a_over_r_n2 * (g_t * sin_m - h_t * cos_m) * q * s;
            }
        }

        x_prime += a_over_r_n2 * x_tmp;
        z_prime += (n as f32 + 1.0) * a_over_r_n2 * z_tmp;
        a_over_r_n2 *= ratio;
        q_prev2 = q_prev1;
        q_prev1 = q;
        rotate_row_indices(&mut i_prev2, &mut i_prev1, &mut i_cur);
    }
    (x_prime, y, z_prime)
}

/// Two-row Schmidt semi-normalized associated Legendre functions for tests
/// (`n = 1..=12`). Derivatives are \(\partial\tilde{P}_n^m/\partial\varphi'\).
#[cfg(test)]
pub(crate) fn schmidt_semi_normalized_associated_legendre(
    cos_theta: f32,
    sin_theta: f32,
) -> ([f32; N_COEFF], [f32; N_COEFF]) {
    let mut p = [[0.0f32; M_LEN]; 3];
    let mut dp = [[0.0f32; M_LEN]; 3];
    p[1][0] = 1.0;
    let mut i_prev2 = 0usize;
    let mut i_prev1 = 1usize;
    let mut i_cur = 2usize;

    let mut psn = [0.0f32; N_COEFF];
    let mut dpsn = [0.0f32; N_COEFF];

    for n in 1..=N_MAX {
        for m in 0..=n {
            let (p_nm, dp_nm) = gauss_norm_p_dp(
                n,
                m,
                cos_theta,
                sin_theta,
                &p[i_prev1],
                &p[i_prev2],
                &dp[i_prev1],
                &dp[i_prev2],
            );
            p[i_cur][m] = p_nm;
            dp[i_cur][m] = dp_nm;
            let s = SCHMIDT_SEMI_NORM[legendre_index(n, m)];
            let ix = coeff_index(n, m);
            psn[ix] = s * p_nm;
            dpsn[ix] = -s * dp_nm;
        }
        rotate_row_indices(&mut i_prev2, &mut i_prev1, &mut i_cur);
    }
    (psn, dpsn)
}

/// `sin(mλ)`, `cos(mλ)` for `m = 0..=12` by the angle-addition recurrence
/// from `sin λ`, `cos λ`.
///
/// \[
/// \sin((m+1)\lambda) = \sin(m\lambda)\cos\lambda + \cos(m\lambda)\sin\lambda \\
/// \cos((m+1)\lambda) = \cos(m\lambda)\cos\lambda - \sin(m\lambda)\sin\lambda
/// \]
pub(crate) fn sin_cos_m_lambda(lambda: f32) -> ([f32; M_LEN], [f32; M_LEN]) {
    let sin_l = sinf(lambda);
    let cos_l = cosf(lambda);
    let mut sin_m = [0.0f32; M_LEN];
    let mut cos_m = [0.0f32; M_LEN];
    cos_m[0] = 1.0;
    for m in 1..=N_MAX {
        sin_m[m] = sin_m[m - 1] * cos_l + cos_m[m - 1] * sin_l;
        cos_m[m] = cos_m[m - 1] * cos_l - sin_m[m - 1] * sin_l;
    }
    (sin_m, cos_m)
}

/// Flat index of WMM Gauss coefficient \(g_n^m\) / \(h_n^m\) (\(n \ge 1\)):
/// \(n(n+1)/2 + m - 1\).
pub(crate) fn coeff_index(n: usize, m: usize) -> usize {
    debug_assert!(n >= 1);
    debug_assert!(n >= m);
    debug_assert!(n <= N_MAX);
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
    fn test_coeff_index(#[case] n: usize, #[case] m: usize, #[case] ix: usize) {
        assert_eq!(coeff_index(n, m), ix);
    }

    #[test]
    fn test_schmidt_semi_norm_matches_recurrence() {
        let mut s = [0.0f32; N_LEGENDRE];
        s[0] = 1.0;
        for n in 1..=N_MAX {
            let i0 = legendre_index(n, 0);
            s[i0] = s[legendre_index(n - 1, 0)] * (2 * n - 1) as f32 / n as f32;
            for m in 1..=n {
                let two_if_m1 = if m == 1 { 2.0 } else { 1.0 };
                s[i0 + m] = s[i0 + m - 1] * sqrtf((n - m + 1) as f32 * two_if_m1 / (n + m) as f32);
            }
        }
        assert_eq!(s, SCHMIDT_SEMI_NORM);
    }

    #[test]
    fn test_k_legendre_matches_formula() {
        let mut k = [0.0f32; N_LEGENDRE];
        for n in 1..=N_MAX {
            for m in 0..=n {
                if n > 1 && m <= n - 2 {
                    k[legendre_index(n, m)] =
                        (((n - 1) * (n - 1) - m * m) as f32) / (((2 * n - 1) * (2 * n - 3)) as f32);
                }
            }
        }
        assert_eq!(k, LEGENDRE_K);
    }

    fn assert_close(got: f32, expected: f32) {
        let tol = 1e-5 * (1.0 + expected.abs());
        assert!(
            (got - expected).abs() <= tol,
            "{got} !≈ {expected} (tol {tol})"
        );
    }

    #[rstest]
    #[case(0.0)]
    #[case(0.5)]
    #[case(-0.5)]
    #[case(core::f32::consts::FRAC_1_SQRT_2)]
    #[case(-core::f32::consts::FRAC_1_SQRT_2)]
    fn test_schmidt_low_degree(#[case] cos_theta: f32) {
        let sin_theta = sqrtf((1.0 - cos_theta) * (1.0 + cos_theta));
        let (p, dp) = schmidt_semi_normalized_associated_legendre(cos_theta, sin_theta);

        // \(\tilde{P}_1^0 = \cos\theta\), \(\partial/\partial\varphi' = \sin\theta\)
        assert_close(p[coeff_index(1, 0)], cos_theta);
        assert_close(dp[coeff_index(1, 0)], sin_theta);

        // \(\tilde{P}_1^1 = \sin\theta\), \(\partial/\partial\varphi' = -\cos\theta\)
        assert_close(p[coeff_index(1, 1)], sin_theta);
        assert_close(dp[coeff_index(1, 1)], -cos_theta);

        // \(\tilde{P}_2^0 = (3\cos^2\theta-1)/2\), \(\partial/\partial\varphi' = 3\cos\theta\sin\theta\)
        assert_close(
            p[coeff_index(2, 0)],
            (3.0 * cos_theta * cos_theta - 1.0) / 2.0,
        );
        assert_close(dp[coeff_index(2, 0)], 3.0 * cos_theta * sin_theta);

        // \(\tilde{P}_2^1 = \sqrt{3}\,\cos\theta\sin\theta\), \(\partial/\partial\varphi' = \sqrt{3}(\sin^2\theta-\cos^2\theta)\)
        let sqrt3 = libm::sqrtf(3.0);
        assert_close(p[coeff_index(2, 1)], sqrt3 * cos_theta * sin_theta);
        assert_close(
            dp[coeff_index(2, 1)],
            sqrt3 * (sin_theta * sin_theta - cos_theta * cos_theta),
        );

        // \(\tilde{P}_2^2 = (\sqrt{3}/2)\sin^2\theta\), \(\partial/\partial\varphi' = -\sqrt{3}\,\cos\theta\sin\theta\)
        assert_close(p[coeff_index(2, 2)], (sqrt3 / 2.0) * sin_theta * sin_theta);
        assert_close(dp[coeff_index(2, 2)], -sqrt3 * cos_theta * sin_theta);
    }

    #[rstest]
    #[case(0.0)]
    #[case(0.3)]
    #[case(-1.2)]
    #[case(core::f32::consts::PI)]
    fn test_sin_cos_m_lambda(#[case] lambda: f32) {
        let (sin_m, cos_m) = sin_cos_m_lambda(lambda);
        assert_close(sin_m[0], 0.0);
        assert_close(cos_m[0], 1.0);
        for m in 1..=N_MAX {
            assert_close(sin_m[m], sinf(m as f32 * lambda));
            assert_close(cos_m[m], cosf(m as f32 * lambda));
        }
    }

    #[test]
    fn test_synthesize_near_pole_y_finite_when_sin_theta_zero() {
        let mut g = [0.0f32; N_COEFF];
        let mut h = [0.0f32; N_COEFF];
        g[coeff_index(1, 0)] = -29404.0;
        g[coeff_index(1, 1)] = -1450.0;
        h[coeff_index(1, 1)] = 4652.0;
        let zeros = [0.0f32; N_COEFF];
        let (x, y, z) = synthesize_near_pole(1.0, 0.0, 0.5, 1.0, &g, &h, &zeros, &zeros, 0.0);
        assert!(x.is_finite() && y.is_finite() && z.is_finite());
        assert_ne!(y, 0.0);
    }

    #[test]
    fn test_synthesize_near_pole_matches_y_prime_over_sin_theta() {
        let sin_theta = 1e-3f32;
        let cos_theta = sqrtf((1.0 - sin_theta) * (1.0 + sin_theta));
        let lambda = 0.5f32;
        let ratio = 1.01f32;
        let mut g = [0.0f32; N_COEFF];
        let mut h = [0.0f32; N_COEFF];
        g[coeff_index(1, 0)] = -29404.0;
        g[coeff_index(1, 1)] = -1450.0;
        h[coeff_index(1, 1)] = 4652.0;
        g[coeff_index(2, 0)] = -2500.0;
        g[coeff_index(2, 1)] = 2982.0;
        h[coeff_index(2, 1)] = -2991.0;
        let zeros = [0.0f32; N_COEFF];
        let (sin_m, cos_m) = sin_cos_m_lambda(lambda);
        let (x_prime, y_prime, z_prime) = synthesize_xyz_prime(
            cos_theta, sin_theta, &sin_m, &cos_m, ratio, &g, &h, &zeros, &zeros, 0.0,
        );
        let (x, y, z) = synthesize_near_pole(
            cos_theta, sin_theta, lambda, ratio, &g, &h, &zeros, &zeros, 0.0,
        );
        assert_close(x, x_prime);
        assert_close(y, y_prime / sin_theta);
        assert_close(z, z_prime);
    }
}
