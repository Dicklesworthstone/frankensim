//! Isothermal tangent-plane searches using the original Peng--Robinson EOS.
//!
//! With `d_i = ln(z_i) + ln(phi_i(z))`, the dimensionless tangent-plane
//! distance is `sum w_i * (ln(w_i) + ln(phi_i(w)) - d_i)`. A negative value
//! exhibits a lower-Gibbs phase split. Two converged nonnegative stationary
//! searches do NOT prove a global minimum over every possible phase.
//!
//! Michelsen, Fluid Phase Equilibria 9 (1982), 1--19,
//! <https://doi.org/10.1016/0378-3812(82)85001-2>.

use super::flash::{PengRobinsonFlashErrorV1 as E, set_composition, wilson_log_k};
use super::{PengRobinsonMixtureV1, PengRobinsonPhaseRootV1, PengRobinsonPhaseV1};
use fs_math::det;
use fs_qty::{Pressure, Temperature};

/// Bounds for each of two incipient-phase searches.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PengRobinsonStabilityOptionsV1 {
    /// Between 1 and 1024 fugacity evaluations per search.
    pub maximum_iterations: usize,
    /// Dimensionless chemical-potential stationarity tolerance.
    pub stationarity_tolerance: f64,
    /// A witness must be more negative than this numerical tolerance.
    pub tangent_plane_tolerance: f64,
}
impl Default for PengRobinsonStabilityOptionsV1 {
    fn default() -> Self {
        Self {
            maximum_iterations: 256,
            stationarity_tolerance: 1.0e-10,
            tangent_plane_tolerance: 1.0e-10,
        }
    }
}

/// One converged stationary search, including its original EOS phase root.
#[derive(Clone, Debug, PartialEq)]
pub struct PengRobinsonStabilityTrialV1 {
    mole_fractions: Vec<f64>,
    root: PengRobinsonPhaseRootV1,
    tangent_plane_distance_over_rt: f64,
    stationarity_residual: f64,
    iterations: usize,
}
impl PengRobinsonStabilityTrialV1 {
    /// Incipient-phase composition in the original component order.
    #[must_use]
    pub fn mole_fractions(&self) -> &[f64] {
        &self.mole_fractions
    }
    /// Minimum-Gibbs mechanical root at this trial composition.
    #[must_use]
    pub const fn root(&self) -> &PengRobinsonPhaseRootV1 {
        &self.root
    }
    /// Dimensionless Gibbs tangent-plane distance relative to the feed root.
    #[must_use]
    pub const fn tangent_plane_distance_over_rt(&self) -> f64 {
        self.tangent_plane_distance_over_rt
    }
    /// Largest departure of a present component from the mean tangent force.
    #[must_use]
    pub const fn stationarity_residual(&self) -> f64 {
        self.stationarity_residual
    }
    /// Number of trial EOS evaluations in this search.
    #[must_use]
    pub const fn iterations(&self) -> usize {
        self.iterations
    }
}

/// Two deterministic Wilson-seeded searches. `instability_witness` is evidence
/// against the feed phase; absence of a witness is explicitly not certification
/// of stability against all phases or compositions.
#[derive(Clone, Debug, PartialEq)]
pub struct PengRobinsonStabilityV1 {
    feed: PengRobinsonPhaseRootV1,
    vapor_seed: PengRobinsonStabilityTrialV1,
    liquid_seed: PengRobinsonStabilityTrialV1,
    tangent_plane_tolerance: f64,
}
impl PengRobinsonStabilityV1 {
    /// Minimum-Gibbs root at the unchanged feed composition.
    #[must_use]
    pub const fn feed(&self) -> &PengRobinsonPhaseRootV1 {
        &self.feed
    }
    /// Stationary search initialized by normalized `z_i K_i`.
    #[must_use]
    pub const fn vapor_seed(&self) -> &PengRobinsonStabilityTrialV1 {
        &self.vapor_seed
    }
    /// Stationary search initialized by normalized `z_i / K_i`.
    #[must_use]
    pub const fn liquid_seed(&self) -> &PengRobinsonStabilityTrialV1 {
        &self.liquid_seed
    }
    /// Most negative searched distance if it exceeds the declared numerical
    /// tolerance. `None` means no witness found, not a global stability proof.
    #[must_use]
    pub fn instability_witness(&self) -> Option<&PengRobinsonStabilityTrialV1> {
        let trial = if self.vapor_seed.tangent_plane_distance_over_rt
            <= self.liquid_seed.tangent_plane_distance_over_rt
        {
            &self.vapor_seed
        } else {
            &self.liquid_seed
        };
        (trial.tangent_plane_distance_over_rt < -self.tangent_plane_tolerance).then_some(trial)
    }
}

impl PengRobinsonMixtureV1 {
    /// Search for negative tangent-plane distance from vapor-like and
    /// liquid-like Wilson seeds. Exact feed absences remain absent; numerical
    /// underflow of a present component refuses rather than deleting species.
    pub fn stability_pt(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        options: PengRobinsonStabilityOptionsV1,
    ) -> Result<PengRobinsonStabilityV1, E> {
        self.stability_pt_cancellable(temperature, pressure, options, || false)
    }

    /// Same source-owned search with a checkpoint before every EOS evaluation.
    pub fn stability_pt_cancellable(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        options: PengRobinsonStabilityOptionsV1,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<PengRobinsonStabilityV1, E> {
        if options.maximum_iterations == 0
            || options.maximum_iterations > 1024
            || !options.stationarity_tolerance.is_finite()
            || options.stationarity_tolerance <= 0.0
            || !options.tangent_plane_tolerance.is_finite()
            || options.tangent_plane_tolerance <= 0.0
        {
            return Err(E::InvalidOptions);
        }
        if cancelled() {
            return Err(E::Cancelled);
        }
        let evaluated = self.evaluate_pt(temperature, pressure)?;
        let feed = evaluated.select(PengRobinsonPhaseV1::Stable);
        let k = wilson_log_k(self, temperature, pressure)?;
        let d: Vec<_> = self
            .mole_fractions()
            .iter()
            .zip(feed.log_fugacity_coefficients())
            .map(|(&z, &phi)| {
                if z > 0.0 {
                    det::ln(z) + phi
                } else {
                    f64::NEG_INFINITY
                }
            })
            .collect();
        let vapor_seed =
            self.stability_trial(temperature, pressure, &d, &k, 1.0, options, &mut cancelled)?;
        let liquid_seed =
            self.stability_trial(temperature, pressure, &d, &k, -1.0, options, &mut cancelled)?;
        if cancelled() {
            return Err(E::Cancelled);
        }
        Ok(PengRobinsonStabilityV1 {
            feed: feed.clone(),
            vapor_seed,
            liquid_seed,
            tangent_plane_tolerance: options.tangent_plane_tolerance,
        })
    }

    fn stability_trial(
        &self,
        temperature: Temperature,
        pressure: Pressure,
        feed_potential: &[f64],
        wilson: &[f64],
        sign: f64,
        options: PengRobinsonStabilityOptionsV1,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<PengRobinsonStabilityTrialV1, E> {
        let mut trial = self.clone();
        let mut log_w: Vec<_> = self
            .mole_fractions()
            .iter()
            .zip(wilson)
            .map(|(&z, &k)| {
                if z > 0.0 {
                    det::ln(z) + sign * k
                } else {
                    f64::NEG_INFINITY
                }
            })
            .collect();
        let mut w = vec![0.0; log_w.len()];
        let mut force = vec![0.0; log_w.len()];
        let mut last = f64::INFINITY;
        for iteration in 0..options.maximum_iterations {
            if cancelled() {
                return Err(E::Cancelled);
            }
            let maximum = log_w.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let mut total = 0.0;
            for i in 0..w.len() {
                w[i] = if self.mole_fractions()[i] == 0.0 {
                    0.0
                } else {
                    det::exp(log_w[i] - maximum)
                };
                if !w[i].is_finite() || (self.mole_fractions()[i] > 0.0 && w[i] <= 0.0) {
                    return Err(E::Unrepresentable("incipient-phase composition"));
                }
                total += w[i];
            }
            let log_total = maximum + det::ln(total);
            if !log_total.is_finite() {
                return Err(E::Unrepresentable("incipient-phase normalization"));
            }
            for i in 0..w.len() {
                w[i] /= total;
                if self.mole_fractions()[i] > 0.0 && w[i] <= 0.0 {
                    return Err(E::Unrepresentable("normalized incipient-phase composition"));
                }
                if w[i] > 0.0 {
                    log_w[i] -= log_total;
                }
            }
            set_composition(&mut trial, &w)?;
            let evaluated = trial.evaluate_pt(temperature, pressure)?;
            let root = evaluated.select(PengRobinsonPhaseV1::Stable);
            let mut tpd = 0.0;
            for i in 0..w.len() {
                if w[i] == 0.0 {
                    continue;
                }
                force[i] = log_w[i] + root.log_fugacity_coefficients()[i] - feed_potential[i];
                tpd += w[i] * force[i];
            }
            let mut residual = 0.0_f64;
            for i in 0..w.len() {
                if w[i] > 0.0 {
                    residual = residual.max((force[i] - tpd).abs());
                }
            }
            if !tpd.is_finite() || !residual.is_finite() {
                return Err(E::Unrepresentable("tangent-plane residual"));
            }
            last = residual;
            if residual <= options.stationarity_tolerance {
                return Ok(PengRobinsonStabilityTrialV1 {
                    mole_fractions: w,
                    root: root.clone(),
                    tangent_plane_distance_over_rt: tpd,
                    stationarity_residual: residual,
                    iterations: iteration + 1,
                });
            }
            // Damped successive substitution on log-amounts. Subtracting the
            // mean force keeps the gauge bounded without changing composition.
            for i in 0..w.len() {
                if w[i] > 0.0 {
                    log_w[i] -= 0.5 * (force[i] - tpd);
                }
            }
        }
        Err(E::StabilityIterationStalled {
            iterations: options.maximum_iterations,
            residual: last,
        })
    }
}
