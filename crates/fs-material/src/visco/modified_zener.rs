//! Two-exponent fractional viscoelasticity with an explicit frequency band.
//!
//! `E*(ω) = E0 + (Efit − E0)(iωτ)^α / (1 + (iωτ)^β)` is the
//! five-parameter law used by Lämmlein et al., *Results in Materials* 9
//! (2021), 100137, equations (1)–(3),
//! <https://doi.org/10.1016/j.rinma.2020.100137>.
//! It captures the measured dispersive storage and loss of violin varnishes
//! without forcing the two fitted exponents to agree. At `α = β` it reduces
//! to [`super::FractionalZener`]; at `α = β = 1` it is one Maxwell arm plus
//! an equilibrium spring. With `α > β`, `Efit` is a fit scale, **not** a
//! finite instantaneous modulus. There is no extrapolation or curing model.
//!
//! The admissible region `0 < β ≤ α ≤ 1` gives positive storage and loss at
//! positive frequency. Runtime histories reuse the existing positive Prony
//! owner. Its finite-band approximation checks storage and loss independently;
//! the reported maxima are measured on a dense grid, not continuum enclosures.

use super::{LoweredModel, ScalarLoweringTolerance, ViscoError, lower_scalar_to_prony};
use fs_math::det;

/// An admitted small-strain two-exponent law in pascals and seconds.
/// Coefficients and applicability are immutable after admission.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ModifiedFractionalZener {
    e0_pa: f64,
    e_fit_pa: f64,
    alpha: f64,
    beta: f64,
    tau_s: f64,
    band_hz: (f64, f64),
    band_rad_s: (f64, f64),
}

impl ModifiedFractionalZener {
    /// Admit coefficients and a finite nonzero applicability band in hertz.
    /// A model fitted at a particular temperature or conditioning state must
    /// retain those qualifications in its source-bound consumer.
    ///
    /// # Errors
    /// [`ViscoError::Parameters`] for nonphysical or nonrepresentable input.
    pub fn new(
        e0_pa: f64,
        e_fit_pa: f64,
        alpha: f64,
        beta: f64,
        tau_s: f64,
        band_hz: (f64, f64),
    ) -> Result<Self, ViscoError> {
        if ![e0_pa, e_fit_pa, alpha, beta, tau_s, band_hz.0, band_hz.1]
            .into_iter()
            .all(f64::is_finite)
            || e0_pa <= 0.0
            || e_fit_pa <= e0_pa
            || beta <= 0.0
            || alpha < beta
            || alpha > 1.0
            || tau_s <= 0.0
            || band_hz.0 <= 0.0
            || band_hz.1 <= band_hz.0
        {
            return Err(ViscoError::Parameters {
                what: "modified Zener needs 0 < E0 < Efit, 0 < beta <= alpha <= 1, tau > 0, 0 < f_lo < f_hi",
            });
        }
        let band_rad_s = (
            core::f64::consts::TAU * band_hz.0,
            core::f64::consts::TAU * band_hz.1,
        );
        if !band_rad_s.1.is_finite() {
            return Err(ViscoError::Parameters {
                what: "modified Zener angular-frequency band must remain finite",
            });
        }
        let law = Self {
            e0_pa,
            e_fit_pa,
            alpha,
            beta,
            tau_s,
            band_hz,
            band_rad_s,
        };
        law.modulus(band_rad_s.0)?;
        law.modulus(band_rad_s.1)?;
        Ok(law)
    }

    /// Equilibrium spring modulus [Pa], retained exactly during Prony lowering.
    #[must_use]
    pub const fn equilibrium_modulus_pa(&self) -> f64 {
        self.e0_pa
    }

    /// Source fit scale [Pa]; it is not an infinite-frequency limit if `α > β`.
    #[must_use]
    pub const fn fit_modulus_pa(&self) -> f64 {
        self.e_fit_pa
    }

    /// Numerator fractional exponent.
    #[must_use]
    pub const fn alpha(&self) -> f64 {
        self.alpha
    }

    /// Denominator fractional exponent.
    #[must_use]
    pub const fn beta(&self) -> f64 {
        self.beta
    }

    /// Source relaxation time [s].
    #[must_use]
    pub const fn relaxation_time_s(&self) -> f64 {
        self.tau_s
    }

    /// Admitted applicability band [Hz].
    #[must_use]
    pub const fn band_hz(&self) -> (f64, f64) {
        self.band_hz
    }

    /// Admitted applicability band [rad/s].
    #[must_use]
    pub const fn band_rad_s(&self) -> (f64, f64) {
        self.band_rad_s
    }

    /// Storage and loss modulus `(E′, E″)` [Pa] at `omega` [rad/s].
    ///
    /// # Errors
    /// [`ViscoError::OutOfBand`] outside the supplied source/applicability band;
    /// [`ViscoError::Parameters`] if the requested response is not representable.
    pub fn modulus(&self, omega: f64) -> Result<(f64, f64), ViscoError> {
        if !(omega >= self.band_rad_s.0 && omega <= self.band_rad_s.1) {
            return Err(ViscoError::OutOfBand {
                omega,
                band: self.band_rad_s,
            });
        }
        // Work in log(ωτ), never form a potentially overflowing ωτ. For
        // ωτ > 1 divide numerator and denominator by (ωτ)^(2β); every
        // denominator term is then <= 2 and cancellation is avoided.
        let log_x = det::ln(omega) + det::ln(self.tau_s);
        let half_pi = core::f64::consts::FRAC_PI_2;
        let ca = det::cos(self.alpha * half_pi);
        let sa = det::sin(self.alpha * half_pi);
        let cb = det::cos(self.beta * half_pi);
        let cab = det::cos((self.alpha - self.beta) * half_pi);
        let sab = det::sin((self.alpha - self.beta) * half_pi);
        let (storage, loss) = if log_x <= 0.0 {
            let a = det::exp(self.alpha * log_x);
            let b = det::exp(self.beta * log_x);
            let den = 1.0 + 2.0 * cb * b + b * b;
            (a * (ca + b * cab) / den, a * (sa + b * sab) / den)
        } else {
            let inverse_b = det::exp(-self.beta * log_x);
            let scale = det::exp((self.alpha - self.beta) * log_x);
            let den = 1.0 + 2.0 * cb * inverse_b + inverse_b * inverse_b;
            (
                scale * (ca * inverse_b + cab) / den,
                scale * (sa * inverse_b + sab) / den,
            )
        };
        let delta = self.e_fit_pa - self.e0_pa;
        let ep = self.e0_pa + delta * storage;
        let epp = delta * loss;
        if !(ep > 0.0 && ep.is_finite() && epp >= 0.0 && epp.is_finite()) {
            return Err(ViscoError::Parameters {
                what: "modified Zener response must remain finite and passive on its band",
            });
        }
        Ok((ep, epp))
    }

    /// Material loss factor `η = E″/E′`, without radiation or contact damping.
    pub fn loss_factor(&self, omega: f64) -> Result<f64, ViscoError> {
        let (storage, loss) = self.modulus(omega)?;
        Ok(loss / storage)
    }

    /// Lower on the admitted band to at most `n_terms` positive Maxwell arms.
    /// Both tolerances are relative to their respective component, so a good
    /// stiffness fit cannot conceal an inaccurate decay spectrum. The two-decade
    /// relaxation-time margin includes broad fractional tails beyond this band.
    ///
    /// # Errors
    /// [`ViscoError::Parameters`] for invalid tolerances/term count;
    /// [`ViscoError::FitTolerance`] when either 512-point grid maximum fails.
    pub fn lower_to_prony(
        &self,
        n_terms: usize,
        storage_rel: f64,
        loss_rel: f64,
    ) -> Result<LoweredModel, ViscoError> {
        lower_scalar_to_prony(
            |omega| self.modulus(omega),
            self.e0_pa,
            self.band_hz,
            n_terms,
            100.0,
            ScalarLoweringTolerance {
                complex: None,
                storage: Some(storage_rel),
                loss: Some(loss_rel),
            },
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::visco::{FractionalZener, GeneralizedMaxwell, lower_to_prony};

    // Published Table 2 fits, in order G1/G2/A1/A2/O1/O2. These are modern
    // conditioned film specimens, not measured Stradivari-instrument varnish.
    const FITS: [(f64, f64, f64, f64, f64); 6] = [
        (500e6, 2.08e9, 0.309, 0.288, 0.268),
        (500e6, 2.12e9, 0.346, 0.318, 0.132),
        (70e6, 1.89e9, 0.330, 0.314, 0.125),
        (70e6, 1.83e9, 0.335, 0.313, 0.179),
        (20e6, 1.37e9, 0.391, 0.360, 0.127),
        (20e6, 1.41e9, 0.402, 0.370, 0.075),
    ];

    fn source_law((e0, ef, a, b, tau): (f64, f64, f64, f64, f64)) -> ModifiedFractionalZener {
        ModifiedFractionalZener::new(e0, ef, a, b, tau, (20.0, 20_000.0)).unwrap()
    }

    #[test]
    fn g0_equal_exponents_reduce_to_fractional_zener_and_one_maxwell_arm() {
        for alpha in [0.2, 0.7, 1.0] {
            let law = ModifiedFractionalZener::new(2e9, 5e9, alpha, alpha, 0.003, (20.0, 20_000.0))
                .unwrap();
            let fractional = FractionalZener::new(2e9, 5e9, alpha, 0.003).unwrap();
            let maxwell = GeneralizedMaxwell::new(2e9, vec![(3e9, 0.003)]).unwrap();
            for f in [20.0, 100.0, 1000.0, 20_000.0] {
                let omega = core::f64::consts::TAU * f;
                let (ep, epp) = law.modulus(omega).unwrap();
                let (fp, fpp) = fractional.modulus(omega);
                assert!((ep - fp).abs() / fp < 1e-13);
                assert!((epp - fpp).abs() / fpp < 1e-11);
                if alpha == 1.0 {
                    let (mp, mpp) = maxwell.modulus(omega);
                    assert!((ep - mp).abs() / mp < 1e-13);
                    assert!((epp - mpp).abs() / mpp < 1e-13);
                }
            }
        }
    }

    #[test]
    fn g0_published_spectra_are_passive_dispersive_and_numerically_stable() {
        // Independent complex-arithmetic evaluation of Table 2 at 1 kHz.
        // These are equation checks, not additional measured source samples.
        let at_one_khz = [
            (2.16288588059e9, 1.32446992687e8),
            (2.26157869020e9, 1.68120771267e8),
            (1.88845938660e9, 1.41656174937e8),
            (1.93275048767e9, 1.53507086324e8),
            (1.55486870068e9, 1.44058950841e8),
            (1.56880652665e9, 1.58642242845e8),
        ];
        for (fit, (reference_storage, reference_loss)) in FITS.into_iter().zip(at_one_khz) {
            let law = source_law(fit);
            let (storage, loss) = law.modulus(core::f64::consts::TAU * 1000.0).unwrap();
            assert!((storage - reference_storage).abs() / reference_storage < 1e-10);
            assert!((loss - reference_loss).abs() / reference_loss < 1e-10);
            let mut previous_storage = 0.0;
            let mut previous_eta = f64::INFINITY;
            for i in 0..128 {
                let omega = super::super::lowering_grid_frequency(
                    law.band_rad_s.0,
                    law.band_rad_s.1,
                    i,
                    128,
                );
                let (ep, epp) = law.modulus(omega).unwrap();
                assert!(ep > previous_storage && epp > 0.0);
                assert!(epp / ep < previous_eta);
                previous_storage = ep;
                previous_eta = epp / ep;
            }
        }
        // Product ωτ overflows, although the fractional response is finite.
        let extreme =
            ModifiedFractionalZener::new(1.0, 2.0, 0.4, 0.36, 1e200, (1e-200, 1e200)).unwrap();
        let (ep, epp) = extreme.modulus(extreme.band_rad_s.1).unwrap();
        assert!(ep.is_finite() && epp.is_finite() && epp > 0.0);
    }

    #[test]
    fn g1_all_six_source_fits_lower_with_independent_storage_and_loss_gates() {
        for fit in FITS {
            let law = source_law(fit);
            let lowered = law.lower_to_prony(16, 0.001, 0.001).unwrap();
            assert_eq!(lowered.model.e_inf, law.equilibrium_modulus_pa());
            assert!(lowered.model.terms.len() <= 16);
            assert!(
                lowered
                    .model
                    .terms
                    .iter()
                    .all(|&(e, tau)| e > 0.0 && tau > 0.0)
            );
            assert!(lowered.sup_rel_storage_err <= 0.001);
            assert!(lowered.sup_rel_loss_err <= 0.001);
            // Interleave a denser independent grid with the verification grid;
            // this is an additional measurement, not an analytic enclosure.
            for i in 0..1025 {
                let omega =
                    super::super::lowering_grid_frequency(lowered.band.0, lowered.band.1, i, 1025);
                let (ep, epp) = law.modulus(omega).unwrap();
                let (gp, gpp) = lowered.modulus_checked(omega).unwrap();
                assert!((gp - ep).abs() / ep < 0.001);
                assert!((gpp - epp).abs() / epp < 0.001);
            }
        }
        // A permissive storage gate must not bypass an unattainable loss gate.
        assert!(matches!(
            source_law(FITS[4]).lower_to_prony(4, 1.0, 1e-8),
            Err(ViscoError::FitTolerance { .. })
        ));
    }

    #[test]
    fn g0_band_and_parameter_refusals_are_explicit() {
        let law = source_law(FITS[4]);
        for omega in [
            0.0,
            law.band_rad_s.0 * 0.99,
            law.band_rad_s.1 * 1.01,
            f64::NAN,
            f64::INFINITY,
        ] {
            assert!(matches!(
                law.modulus(omega),
                Err(ViscoError::OutOfBand { .. })
            ));
        }
        for (e0, ef, a, b, tau) in [
            (0.0, 1.0, 0.4, 0.3, 1.0),
            (1.0, 1.0, 0.4, 0.3, 1.0),
            (1.0, 2.0, 0.3, 0.4, 1.0),
            (1.0, 2.0, 1.1, 0.3, 1.0),
            (1.0, 2.0, 0.4, 0.0, 1.0),
            (1.0, 2.0, 0.4, 0.3, f64::NAN),
        ] {
            assert!(ModifiedFractionalZener::new(e0, ef, a, b, tau, (20.0, 20_000.0)).is_err());
        }
        for tol in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(law.lower_to_prony(16, tol, 0.01).is_err());
            assert!(law.lower_to_prony(16, 0.01, tol).is_err());
        }
    }

    #[test]
    fn g3_positive_lowering_is_independent_of_modulus_units() {
        let pa = FractionalZener::new(8e9, 16e9, 0.35, 1e-3).unwrap();
        let mpa = FractionalZener::new(8e3, 16e3, 0.35, 1e-3).unwrap();
        let a = lower_to_prony(&pa, 20.0, 20_000.0, 16, 0.02).unwrap();
        let b = lower_to_prony(&mpa, 20.0, 20_000.0, 16, 0.02).unwrap();
        for f in [20.0, 100.0, 1000.0, 20_000.0] {
            let (ap, app) = a.modulus_checked(core::f64::consts::TAU * f).unwrap();
            let (bp, bpp) = b.modulus_checked(core::f64::consts::TAU * f).unwrap();
            assert!((ap / 1e6 - bp).abs() / bp < 1e-7);
            assert!((app / 1e6 - bpp).abs() / bpp < 1e-7);
        }
    }
}
