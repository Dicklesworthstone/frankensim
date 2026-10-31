//! Frequency-independent (constant loss factor) shell material loss.
//!
//! Wood's measured loss factors are nearly independent of frequency across the audio band. A shell pack whose nine
//! viscosities were written as `V = eta x stiffness / omega_ref` declares that anchor as
//! [`super::CONSTANT_LOSS_REFERENCE_OMEGA`] [rad/s]; its loss stiffness `omega_ref V = eta C` is then constant in
//! frequency. Read as Kelvin-Voigt instead, the same numbers give each mode the loss factor `eta omega / omega_ref`.
//!
//! The time-domain law with that loss stiffness is a generalized Maxwell (Prony) series: [`BRANCHES`] branches whose
//! relaxation rates are log-spaced over [`SPECTRUM_BAND_RAD_S`], weighted by least squares so that the small-loss shape
//! `G(omega) = sum_j g_j x_j / (1 + x_j^2)`, `x_j = omega / omega_j`, is flat at one across that band. A branch of
//! modulus `eta g_j E` on top of `E` is exactly `fs_material::visco::GeneralizedMaxwell` with `E''(omega) = eta G(omega) E`.
//!
//! The modal reduction applies the series at each retained mode's own frequency (the modal strain-energy rule, as
//! `plate_viscosity::hysteretic` does for plates): the Kelvin-Voigt modal matrix `C` of the same viscosities becomes
//! `C_ab s_a s_b` with `s_a = sqrt(omega_ref G(omega_a) / omega_a)`, so mode a's loss factor `C_aa / omega_a` is its
//! strain-energy-weighted `eta` times `G(omega_a)`. A congruence by a positive diagonal keeps `C` positive
//! semidefinite. This legacy modal approximation leaves storage fixed. Selecting
//! [`relaxation::SubstrateRelaxationOptions`] instead installs the SAME positive spectrum as
//! causal material histories, including its storage dispersion and an explicit stiffness reference.
use fs_math::det;

/// Causal substrate storage/loss in the original shell and body-memory owners.
pub mod relaxation;

/// Maxwell branches of the unit loss spectrum.
pub const BRANCHES: usize = 5;
/// The span of the branch relaxation rates and of the least-squares fit [rad/s]: 50 Hz to 10 kHz.
pub const SPECTRUM_BAND_RAD_S: (f64, f64) = (core::f64::consts::TAU * 50.0, core::f64::consts::TAU * 10_000.0);
/// Log-spaced angular frequencies of the least-squares fit.
const FIT_POINTS: usize = 401;

/// The unit (small-loss) constant-loss spectrum: branch relaxation rates [rad/s] and weights.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConstantLossSpectrum {
    rates_rad_s: [f64; BRANCHES],
    weights: [f64; BRANCHES],
}

impl ConstantLossSpectrum {
    /// The least-squares fit of `G = 1` over [`SPECTRUM_BAND_RAD_S`], by Cholesky on the normal equations.
    #[must_use]
    pub fn fit() -> Self {
        Self::fit_over(SPECTRUM_BAND_RAD_S)
    }

    /// The same fit with the branch rates and the fit points spread over `band` [rad/s] instead (for owners whose
    /// retained modes reach past 10 kHz); `fit()` is `fit_over(SPECTRUM_BAND_RAD_S)`.
    #[must_use]
    pub fn fit_over(band: (f64, f64)) -> Self {
        let (lo, hi) = band;
        let span = det::ln(hi / lo);
        let rates: [f64; BRANCHES] =
            core::array::from_fn(|j| lo * det::exp(span * j as f64 / (BRANCHES - 1) as f64));
        let mut normal = [[0.0; BRANCHES]; BRANCHES];
        let mut rhs = [0.0; BRANCHES];
        for k in 0..FIT_POINTS {
            let omega = lo * det::exp(span * k as f64 / (FIT_POINTS - 1) as f64);
            let row: [f64; BRANCHES] = core::array::from_fn(|j| debye(omega / rates[j]));
            for i in 0..BRANCHES {
                rhs[i] += row[i];
                for j in 0..BRANCHES {
                    normal[i][j] += row[i] * row[j];
                }
            }
        }
        // Cholesky: normal = L L^T, then two triangular solves.
        let mut l = [[0.0; BRANCHES]; BRANCHES];
        for i in 0..BRANCHES {
            for j in 0..=i {
                let mut sum = normal[i][j];
                for k in 0..j {
                    sum -= l[i][k] * l[j][k];
                }
                l[i][j] = if i == j { det::sqrt(sum) } else { sum / l[j][j] };
            }
        }
        let mut y = [0.0; BRANCHES];
        for i in 0..BRANCHES {
            let mut sum = rhs[i];
            for k in 0..i {
                sum -= l[i][k] * y[k];
            }
            y[i] = sum / l[i][i];
        }
        let mut weights = [0.0; BRANCHES];
        for i in (0..BRANCHES).rev() {
            let mut sum = y[i];
            for k in i + 1..BRANCHES {
                sum -= l[k][i] * weights[k];
            }
            weights[i] = sum / l[i][i];
        }
        Self { rates_rad_s: rates, weights }
    }

    /// `G(omega)`: the loss factor at `omega` [rad/s] per unit tabulated loss factor.
    #[must_use]
    pub fn shape(&self, omega: f64) -> f64 {
        self.rates_rad_s.iter().zip(&self.weights).map(|(&rate, &g)| g * debye(omega / rate)).sum()
    }

    /// Storage increment per unit loss stiffness, paired with `shape(omega)`.
    /// This is sum g_j (omega tau_j)^2/(1+(omega tau_j)^2), not an
    /// independently fitted frequency correction. Finite frequencies use a
    /// reciprocal form above each relaxation rate to avoid overflow.
    #[must_use]
    pub fn storage_shape(&self, omega: f64) -> f64 {
        self.rates_rad_s.iter().zip(&self.weights)
            .map(|(&rate, &g)| g * debye_storage(omega / rate)).sum()
    }

    /// `(relaxation rate [rad/s], weight)` per branch; the relaxation time is the rate's inverse.
    pub fn branches(&self) -> impl Iterator<Item = (f64, f64)> + '_ {
        self.rates_rad_s.iter().copied().zip(self.weights.iter().copied())
    }
}

/// One Maxwell branch's loss per unit weight, `x / (1 + x^2)` at `x = omega tau`.
fn debye(x: f64) -> f64 {
    // Preserve the original audio-band arithmetic; use the reciprocal only
    // when squaring a very large frequency ratio would lose its finite tail.
    if x.abs() <= 1.0e150 {
        x / (1.0 + x * x)
    } else {
        let inverse = 1.0 / x;
        inverse / (1.0 + inverse * inverse)
    }
}

fn debye_storage(x: f64) -> f64 {
    if x.abs() <= 1.0 {
        let square = x * x;
        square / (1.0 + square)
    } else {
        let inverse = 1.0 / x;
        1.0 / (1.0 + inverse * inverse)
    }
}

/// Rescale a row-major `n x n` Kelvin-Voigt modal viscosity matrix [1/s] (mass-normalized modes at angular frequencies
/// `omega`, viscosities written at `reference_omega`) into the constant-loss law at each mode's own frequency.
pub(super) fn rescale(matrix: &mut [f64], omega: &[f64], reference_omega: f64, spectrum: &ConstantLossSpectrum) {
    let n = omega.len();
    debug_assert_eq!(matrix.len(), n * n);
    let s: Vec<f64> = omega.iter().map(|&w| det::sqrt(reference_omega * spectrum.shape(w) / w)).collect();
    for (a, row) in matrix.chunks_mut(n).enumerate() {
        for (b, value) in row.iter_mut().enumerate() {
            *value *= s[a] * s[b];
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_material::visco::GeneralizedMaxwell;

    fn log_grid(lo_hz: f64, hi_hz: f64, points: usize) -> impl Iterator<Item = f64> {
        (0..points).map(move |k| {
            core::f64::consts::TAU * lo_hz * (hi_hz / lo_hz).powf(k as f64 / (points - 1) as f64)
        })
    }

    // The receipt: every weight is positive, and G stays within 1.5 % of one over the violin body's modal window
    // (100 Hz-2 kHz; measured 0.9903..1.0129 in the numpy probe) and within 4.5 % over the whole fit band
    // (measured 0.9585..1.0143).
    #[test]
    fn the_five_branch_spectrum_is_flat_across_the_body_band() {
        let spectrum = ConstantLossSpectrum::fit();
        assert!(spectrum.branches().all(|(_, g)| g > 0.0), "{spectrum:?}");
        let extreme = |lo, hi| {
            log_grid(lo, hi, 2001).map(|w| spectrum.shape(w)).fold((f64::MAX, f64::MIN), |(a, b), g| (a.min(g), b.max(g)))
        };
        let body = extreme(100.0, 2000.0);
        let band = extreme(50.0, 10_000.0);
        eprintln!("G over 100 Hz-2 kHz {:.4}..{:.4}; over 50 Hz-10 kHz {:.4}..{:.4}; branches {:?}",
            body.0, body.1, band.0, band.1, spectrum.branches().collect::<Vec<_>>());
        assert!(body.0 > 0.985 && body.1 < 1.015, "{body:?}");
        assert!(band.0 > 0.955 && band.1 < 1.045, "{band:?}");
    }

    // The same branches as fs-material's Prony owner: modulus E = 1 plus branches eta g_j at tau_j = 1 / rate_j has
    // E''(omega) = eta G(omega) exactly.
    #[test]
    fn the_spectrum_is_the_loss_modulus_of_a_generalized_maxwell_series() {
        let spectrum = ConstantLossSpectrum::fit();
        for eta in [0.0075, 0.0135, 0.0175] {
            let series = GeneralizedMaxwell::new(1.0, spectrum.branches().map(|(rate, g)| (eta * g, 1.0 / rate)).collect())
                .unwrap();
            for w in log_grid(20.0, 20_000.0, 97) {
                let loss = series.modulus(w).1;
                assert!((loss - eta * spectrum.shape(w)).abs() <= 1e-14 * loss, "{w} {loss}");
            }
        }
    }

    // A Kelvin-Voigt matrix written at omega_ref gives mode a the loss factor eta_a omega_a / omega_ref; after the
    // rescale each mode carries eta_a G(omega_a), and the matrix stays positive definite.
    #[test]
    fn rescaling_gives_each_mode_its_constant_loss_factor() {
        let spectrum = ConstantLossSpectrum::fit();
        let reference = core::f64::consts::TAU * 200.0;
        let omega = [core::f64::consts::TAU * 473.5, core::f64::consts::TAU * 991.0, core::f64::consts::TAU * 1950.0];
        let eta = [0.0135, 0.0075, 0.0162];
        // Strain-energy-weighted loss stiffness H (positive definite) and its Kelvin-Voigt modal matrix H / omega_ref.
        let h = [
            [eta[0] * omega[0] * omega[0], 0.3e3, -0.1e3],
            [0.3e3, eta[1] * omega[1] * omega[1], 0.2e3],
            [-0.1e3, 0.2e3, eta[2] * omega[2] * omega[2]],
        ];
        let mut c: Vec<f64> = h.iter().flatten().map(|v| v / reference).collect();
        for a in 0..3 {
            let kelvin_voigt = c[a * 3 + a] / omega[a];
            assert!((kelvin_voigt - eta[a] * omega[a] / reference).abs() < 1e-12 * kelvin_voigt);
        }
        rescale(&mut c, &omega, reference, &spectrum);
        for a in 0..3 {
            let loss = c[a * 3 + a] / omega[a];
            assert!((loss - eta[a] * spectrum.shape(omega[a])).abs() < 1e-12 * loss, "mode {a}: {loss}");
            for b in 0..3 {
                assert_eq!(c[a * 3 + b], c[b * 3 + a]);
            }
        }
        // Leading minors stay positive.
        let m2 = c[0] * c[4] - c[1] * c[3];
        let m3 = c[0] * (c[4] * c[8] - c[5] * c[7]) - c[1] * (c[3] * c[8] - c[5] * c[6]) + c[2] * (c[3] * c[7] - c[4] * c[6]);
        assert!(c[0] > 0.0 && m2 > 0.0 && m3 > 0.0);
    }
}
