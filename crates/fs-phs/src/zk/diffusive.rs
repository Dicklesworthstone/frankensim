//! Passive circular-tube wall memory from the exact Bessel continued fraction.
//!
//! For `t = s rho a²/mu`, the normalized series impedance is
//! `Z/L0 = s + 8/tau + (2/tau) t/(6+t/(8+t/(10+...)))`.
//! Thermal memory uses `F(t) = 2/(2+t/(4+t/(6+...)))`, with
//! `tau_t = Pr tau` and `Y/C0 = s [1+(gamma-1) F(s tau_t)]`.
//! These are DLMF 10.10.1 under the modified-Bessel substitution, the
//! SAME constitutive law as the parent module's frequency-domain oracle.
//!
//! An even truncation of each positive continued fraction is contracted to
//! a positive Jacobi matrix. Its spectral measure gives positive Foster
//! arms `g s/(s+p)`. This retains resistance AND stored energy, without
//! fitting real loss alone, clipping residues, or inventing a tail load.
//! Increasing the pole count is the physical diffusion discretization axis.
//! The declared band is admitted against the existing exact Bessel owner;
//! passivity alone never establishes approximation accuracy.
//!
//! This module supplies coefficients, not a time integrator. In particular,
//! a DG consumer can append these memory states to its existing exponential.

use super::{zk_one_minus_f, zwikker_kosten_f};
use crate::PhsError;
use fs_la::eigen::jacobi_eigh;
use fs_math::{c64::C64, det};

/// Bounded cold eigensolve; a consumer must also admit its augmented state size.
pub const MAX_ZK_DIFFUSIVE_POLES: usize = 64;

/// One positive term `g s/(s+p)` normalized by the inviscid L or C.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ZkDiffusiveTerm {
    /// Normalized resistance or conductance `g` [1/s].
    pub gain_per_s: f64,
    /// Memory relaxation rate `p` [1/s].
    pub rate_per_s: f64,
}

/// Admitted finite-band approximation of the exact cylindrical wall law.
///
/// Multiply the series response by `rho/A`, and the shunt response by
/// `A/(rho*c²)`, for the physical per-length telegraph pair. The returned
/// responses INCLUDE the inviscid `s` term. A consumer which already owns
/// acoustic inertia/compliance adds only the declared drag and memory arms.
#[derive(Clone, Debug, PartialEq)]
pub struct ZkDiffusive {
    band_rad_s: [f64; 2],
    series_drag_per_s: f64,
    series: Vec<ZkDiffusiveTerm>,
    shunt: Vec<ZkDiffusiveTerm>,
    max_complex_error: f64,
    max_real_loss_error: f64,
}

/// Reusable dimensionless circular diffusion spectra at one pole resolution.
///
/// Radius and gas properties scale the SAME two Jacobi spectral measures.
/// A field or radius sweep can therefore prepare these eigensolves once and
/// [`bind`](Self::bind) each local physical state without repeating them.
/// Every binding still performs the original 129-frequency complex/loss
/// accuracy check; a successful state never admits a different state or band.
#[derive(Clone, Debug, PartialEq)]
pub struct ZkDiffusiveBasis {
    series: Vec<(f64, f64)>,
    thermal: Vec<(f64, f64)>,
}

impl ZkDiffusiveBasis {
    /// Prepare the original positive spectral measures, with no physical gas
    /// or frequency band implied. The original pole cap and residual checks
    /// apply. This performs exactly two bounded real symmetric eigensolves.
    pub fn new(poles: usize) -> Result<Self, PhsError> {
        if !(1..=MAX_ZK_DIFFUSIVE_POLES).contains(&poles) {
            return Err(invalid("ZK diffusive physical data, pole count or band"));
        }
        Ok(Self {
            series: measure(6.0, poles)?,
            thermal: measure(2.0, poles)?,
        })
    }

    /// Number of positive viscous and thermal poles in each spectrum.
    #[must_use]
    pub fn poles_per_field(&self) -> usize {
        self.series.len()
    }

    /// Scale and admit one local radius/gas/band through the original owner.
    /// Returns the same coefficients, response and sampled errors as
    /// [`ZkDiffusive::new`] with this pole count, including gamma=1 behavior.
    #[allow(clippy::too_many_arguments)]
    pub fn bind(
        &self,
        radius_m: f64,
        density_kg_m3: f64,
        viscosity_pa_s: f64,
        gamma: f64,
        prandtl: f64,
        band_rad_s: [f64; 2],
        maximum_relative_error: f64,
    ) -> Result<ZkDiffusive, PhsError> {
        ZkDiffusive::build(
            radius_m,
            density_kg_m3,
            viscosity_pa_s,
            gamma,
            prandtl,
            self.series.len(),
            band_rad_s,
            maximum_relative_error,
            Some(self),
        )
    }
}

fn invalid(what: &'static str) -> PhsError {
    PhsError::NotPsd { what }
}

/// Positive Jacobi spectral measure of
/// `a0/(a0+t/(a0+2+t/(a0+4+...))))` with `2*n` denominators.
fn measure(a0: f64, n: usize) -> Result<Vec<(f64, f64)>, PhsError> {
    let c: Vec<f64> = (0..2 * n - 1)
        .map(|i| {
            let a = a0 + 2.0 * i as f64;
            1.0 / (a * (a + 2.0))
        })
        .collect();
    let mut matrix = vec![0.0; n * n];
    matrix[0] = c[0];
    for i in 1..n {
        matrix[i * n + i] = c[2 * i - 1] + c[2 * i];
    }
    for i in 0..n - 1 {
        let value = det::sqrt(c[2 * i] * c[2 * i + 1]);
        matrix[i * n + i + 1] = value;
        matrix[(i + 1) * n + i] = value;
    }
    let (values, vectors) = jacobi_eigh(&matrix, n);
    let scale = matrix.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
    let mut out = Vec::with_capacity(n);
    let mut total_weight = 0.0;
    for k in 0..n {
        let lambda = values[k];
        let weight = vectors[k] * vectors[k];
        if !lambda.is_finite() || lambda <= 0.0 || !weight.is_finite() || weight <= 0.0 {
            return Err(invalid(
                "ZK diffusive spectral measure is not strictly positive",
            ));
        }
        let mut residual = 0.0_f64;
        let mut norm = 0.0;
        for i in 0..n {
            let av: f64 = (0..n).map(|j| matrix[i * n + j] * vectors[j * n + k]).sum();
            residual = residual.max((av - lambda * vectors[i * n + k]).abs());
            norm += vectors[i * n + k] * vectors[i * n + k];
        }
        if residual > 2e-12 * scale || (norm - 1.0).abs() > 2e-12 {
            return Err(invalid("ZK diffusive spectral decomposition is unresolved"));
        }
        out.push((1.0 / lambda, weight));
        total_weight += weight;
    }
    if (total_weight - 1.0).abs() > 2e-12 {
        return Err(invalid(
            "ZK diffusive measure lost its static storage moment",
        ));
    }
    Ok(out)
}

fn response(terms: &[ZkDiffusiveTerm], omega: f64) -> C64 {
    let s = C64::new(0.0, omega);
    terms.iter().fold(C64::ZERO, |sum, term| {
        sum + (s / (s + C64::from_re(term.rate_per_s))).scale(term.gain_per_s)
    })
}

impl ZkDiffusive {
    /// Build positive viscous and thermal memories from physical gas data.
    ///
    /// The fixed 129-point logarithmic band check bounds the sampled relative
    /// error of BOTH complex wall operators AND their real dissipative parts.
    /// It does not claim an interval enclosure or accuracy outside the band.
    /// `maximum_relative_error` must lie strictly between zero and one; increase
    /// the supplied pole count or narrow the physical band when admission fails.
    /// No coefficient or error budget is silently repaired.
    ///
    /// # Errors
    /// Invalid physical data/band/pole count; unresolved positive decomposition;
    /// an exact-Bessel shear refusal; or either sampled error above the budget.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        radius_m: f64,
        density_kg_m3: f64,
        viscosity_pa_s: f64,
        gamma: f64,
        prandtl: f64,
        poles: usize,
        band_rad_s: [f64; 2],
        maximum_relative_error: f64,
    ) -> Result<Self, PhsError> {
        Self::build(
            radius_m,
            density_kg_m3,
            viscosity_pa_s,
            gamma,
            prandtl,
            poles,
            band_rad_s,
            maximum_relative_error,
            None,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn build(
        radius_m: f64,
        density_kg_m3: f64,
        viscosity_pa_s: f64,
        gamma: f64,
        prandtl: f64,
        poles: usize,
        band_rad_s: [f64; 2],
        maximum_relative_error: f64,
        basis: Option<&ZkDiffusiveBasis>,
    ) -> Result<Self, PhsError> {
        if !(1..=MAX_ZK_DIFFUSIVE_POLES).contains(&poles)
            || ![
                radius_m,
                density_kg_m3,
                viscosity_pa_s,
                prandtl,
                band_rad_s[0],
                band_rad_s[1],
            ]
            .iter()
            .all(|v| v.is_finite() && *v > 0.0)
            || !gamma.is_finite()
            || gamma < 1.0
            || band_rad_s[1] <= band_rad_s[0]
            || !maximum_relative_error.is_finite()
            || maximum_relative_error <= 0.0
            || maximum_relative_error >= 1.0
        {
            return Err(invalid("ZK diffusive physical data, pole count or band"));
        }
        let tau_s = density_kg_m3 * radius_m * radius_m / viscosity_pa_s;
        let thermal_tau = tau_s * prandtl;
        if ![
            tau_s,
            thermal_tau,
            8.0 / tau_s,
            band_rad_s[1] / band_rad_s[0],
        ]
        .iter()
        .all(|v| v.is_finite() && *v > 0.0)
        {
            return Err(invalid(
                "ZK diffusive physical scales are not representable",
            ));
        }
        let direct_series;
        let series_measure = if let Some(basis) = basis {
            &basis.series
        } else {
            direct_series = measure(6.0, poles)?;
            &direct_series
        };
        let series = series_measure
            .iter()
            .map(|&(pole, weight)| ZkDiffusiveTerm {
                gain_per_s: weight * pole / (3.0 * tau_s),
                rate_per_s: pole / tau_s,
            })
            .collect::<Vec<_>>();
        let shunt = if gamma == 1.0 {
            Vec::new()
        } else {
            let direct_thermal;
            let thermal_measure = if let Some(basis) = basis {
                &basis.thermal
            } else {
                direct_thermal = measure(2.0, poles)?;
                &direct_thermal
            };
            thermal_measure
                .iter()
                .map(|&(pole, weight)| ZkDiffusiveTerm {
                    gain_per_s: (gamma - 1.0) * weight * pole / thermal_tau,
                    rate_per_s: pole / thermal_tau,
                })
                .collect::<Vec<_>>()
        };
        if series.iter().chain(&shunt).any(|t| {
            !t.gain_per_s.is_finite()
                || t.gain_per_s <= 0.0
                || !t.rate_per_s.is_finite()
                || t.rate_per_s <= 0.0
        }) {
            return Err(invalid(
                "ZK diffusive scaled memory is not positive and finite",
            ));
        }
        let mut model = Self {
            band_rad_s,
            series_drag_per_s: 8.0 / tau_s,
            series,
            shunt,
            max_complex_error: 0.0,
            max_real_loss_error: 0.0,
        };
        for i in 0..=128 {
            let omega =
                band_rad_s[0] * det::pow(band_rad_s[1] / band_rad_s[0], f64::from(i) / 128.0);
            let rv = det::sqrt(omega * tau_s);
            let s = C64::new(0.0, omega);
            let exact_series = s / zk_one_minus_f(rv)? - s;
            let exact_shunt = s * zwikker_kosten_f(rv * det::sqrt(prandtl))?.scale(gamma - 1.0);
            let actual_series =
                C64::from_re(model.series_drag_per_s) + response(&model.series, omega);
            let actual_shunt = response(&model.shunt, omega);
            for (actual, exact) in [(actual_series, exact_series), (actual_shunt, exact_shunt)] {
                if exact == C64::ZERO && actual == C64::ZERO {
                    continue;
                }
                if !exact.re.is_finite() || exact.re <= 0.0 {
                    return Err(invalid(
                        "ZK diffusive exact dissipative reference must be positive",
                    ));
                }
                let complex = (actual - exact).abs() / exact.abs();
                let real = (actual.re - exact.re).abs() / exact.re;
                if !complex.is_finite()
                    || !real.is_finite()
                    || complex > maximum_relative_error
                    || real > maximum_relative_error
                {
                    return Err(invalid(
                        "ZK diffusive pole count does not resolve the declared band",
                    ));
                }
                model.max_complex_error = model.max_complex_error.max(complex);
                model.max_real_loss_error = model.max_real_loss_error.max(real);
            }
        }
        Ok(model)
    }

    /// Poiseuille resistance divided by inviscid inertance [1/s].
    #[must_use]
    pub const fn series_drag_per_s(&self) -> f64 {
        self.series_drag_per_s
    }
    /// Viscous wall arms; each retains positive inertive storage `gain/rate`.
    #[must_use]
    pub fn series_terms(&self) -> &[ZkDiffusiveTerm] {
        &self.series
    }
    /// Thermal wall arms; each retains positive compliant storage `gain/rate`.
    #[must_use]
    pub fn shunt_terms(&self) -> &[ZkDiffusiveTerm] {
        &self.shunt
    }
    /// The admitted physical angular-frequency band [rad/s].
    #[must_use]
    pub const fn band_rad_s(&self) -> [f64; 2] {
        self.band_rad_s
    }
    /// Largest sampled relative error of the complex EXCESS wall operators.
    #[must_use]
    pub const fn max_complex_error(&self) -> f64 {
        self.max_complex_error
    }
    /// Largest sampled relative resistance/conductance error, separately from storage.
    #[must_use]
    pub const fn max_real_loss_error(&self) -> f64 {
        self.max_real_loss_error
    }

    /// Full normalized telegraph pair `(Z/L0,Y/C0)` in the `+i omega` convention.
    ///
    /// # Errors
    /// Frequency outside the declared and sampled admission band.
    pub fn response_at(&self, omega: f64) -> Result<(C64, C64), PhsError> {
        if !omega.is_finite() || !(self.band_rad_s[0]..=self.band_rad_s[1]).contains(&omega) {
            return Err(invalid(
                "ZK diffusive response is outside its admitted band",
            ));
        }
        let s = C64::new(0.0, omega);
        Ok((
            s + C64::from_re(self.series_drag_per_s) + response(&self.series, omega),
            s + response(&self.shunt, omega),
        ))
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
#[path = "diffusive/basis_tests.rs"]
mod basis_tests;
