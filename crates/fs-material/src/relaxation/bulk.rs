//! Positive compression memory from the original molecular heat capacities.
//!
//! For `cv(s) = cv_inf + sum C_j/(1+s*tau_j)`, the matrix inversion lemma gives
//! `K(s) = K_inf - (p R/cv_inf) u^T (s I + D + u u^T)^-1 u`, where
//! `D_jj = 1/tau_j` and `u_j = sqrt(C_j/(cv_inf*tau_j))`.
//! The symmetric positive spectrum therefore realizes the EXACT finite
//! reservoir law as `K_eq + sum delta_K_l s/(s+lambda_l)`. These are physical
//! internal coordinates, not a rational fit, extra bulk viscosity or propagated
//! Stokes--Kirchhoff attenuation. Frozen gas composition and temperature remain
//! the original owner's assumptions. Construction is bounded by 32 reservoirs.

use super::{RelaxationState, debye, positive, state_error};
use crate::MaterialError;
use fs_la::eigen::jacobi_eigh;
use fs_math::{c64::C64, det};

/// One positive Maxwell compression branch in physical pressure/strain units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BulkRelaxationArm {
    /// Increase in instantaneous bulk modulus [Pa].
    pub delta_bulk_modulus_pa: f64,
    /// Relaxation time of compression, after eliminating internal temperatures [s].
    /// This generally differs from an individual caloric reservoir's time.
    pub relaxation_time_s: f64,
}

/// The original molecular law expressed as equilibrium storage plus memory.
#[derive(Clone, Debug, PartialEq)]
pub struct BulkRelaxation {
    source: RelaxationState,
    equilibrium_bulk_modulus_pa: f64,
    frozen_bulk_modulus_pa: f64,
    arms: Vec<BulkRelaxationArm>,
}

impl RelaxationState {
    /// Lower the same evaluated reservoirs to positive time-domain bulk memory.
    ///
    /// Exactly equal relaxation times combine before the small symmetric solve;
    /// zero-strength reservoirs contribute no state. The retained original
    /// source still records every species and process. No spectral window,
    /// fitted pole, guessed rate, or relaxation-strength linearization is used.
    ///
    /// # Errors
    /// Unrepresentable rates, unresolved eigenpairs or lost static/frozen limits.
    /// Callers with cancellation scopes poll before and after this bounded cold
    /// construction; no time integration or runtime state is modified here.
    pub fn bulk_relaxation(&self) -> Result<BulkRelaxation, MaterialError> {
        let cv = self.frozen_cv_j_kg_k();
        let gas = self.gas();
        let scale_bulk = gas.pressure * (gas.specific_gas_constant / cv);
        let equilibrium = self.caloric_response(0.0)?.bulk_modulus_pa.re;
        let frozen = gas.pressure + scale_bulk;
        if !positive(scale_bulk) || !positive(frozen) {
            return Err(state_error("molecular bulk storage is unrepresentable"));
        }
        let mut groups: Vec<(f64, f64)> = Vec::new();
        for process in self.processes() {
            let capacity = process.heat_capacity_j_kg_k;
            if capacity == 0.0 {
                continue;
            }
            if let Some((_, sum)) = groups
                .iter_mut()
                .find(|(tau, _)| *tau == process.relaxation_time_s)
            {
                *sum += capacity;
            } else {
                groups.push((process.relaxation_time_s, capacity));
            }
        }
        let n = groups.len();
        let mut arms = Vec::with_capacity(n);
        if n > 0 {
            let rates: Vec<_> = groups.iter().map(|(tau, _)| tau.recip()).collect();
            if rates.iter().any(|rate| !positive(*rate)) {
                return Err(state_error(
                    "molecular compression rates are unrepresentable",
                ));
            }
            let scale = rates.iter().copied().fold(0.0_f64, f64::max);
            let u: Vec<_> = groups
                .iter()
                .zip(&rates)
                .map(|((_, c), rate)| det::sqrt((c / cv) * (rate / scale)))
                .collect();
            let mut matrix = vec![0.0; n * n];
            for i in 0..n {
                for j in 0..n {
                    matrix[i * n + j] = u[i] * u[j];
                }
                matrix[i * n + i] += rates[i] / scale;
            }
            let norm = matrix.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
            if !positive(norm) || u.iter().any(|v| !positive(*v)) {
                return Err(state_error(
                    "molecular compression spectrum is unrepresentable",
                ));
            }
            let (values, vectors) = jacobi_eigh(&matrix, n);
            let tolerance = 128.0 * f64::EPSILON * n as f64;
            for j in 0..n {
                let lambda = values[j];
                if !positive(lambda) {
                    return Err(state_error("molecular compression spectrum is unresolved"));
                }
                for k in 0..=j {
                    let gram: f64 = (0..n)
                        .map(|i| vectors[i * n + j] * vectors[i * n + k])
                        .sum();
                    if (gram - if j == k { 1.0 } else { 0.0 }).abs() > tolerance {
                        return Err(state_error(
                            "molecular compression eigenbasis is unresolved",
                        ));
                    }
                }
                for i in 0..n {
                    let av: f64 = (0..n).map(|k| matrix[i * n + k] * vectors[k * n + j]).sum();
                    if (av - lambda * vectors[i * n + j]).abs() > tolerance * norm {
                        return Err(state_error("molecular compression eigenpair is unresolved"));
                    }
                }
                let projection: f64 = (0..n).map(|i| u[i] * vectors[i * n + j]).sum();
                if projection == 0.0 {
                    continue;
                }
                let increment = scale_bulk * ((projection * projection) / lambda);
                let rate = lambda * scale;
                let tau = rate.recip();
                if !positive(increment) || !positive(rate) || !positive(tau) {
                    return Err(state_error(
                        "molecular compression branch is unrepresentable",
                    ));
                }
                arms.push(BulkRelaxationArm {
                    delta_bulk_modulus_pa: increment,
                    relaxation_time_s: tau,
                });
            }
            let capacity: f64 = groups.iter().map(|(_, c)| c).sum();
            let target = scale_bulk * (capacity / (cv + capacity));
            let sum: f64 = arms.iter().map(|arm| arm.delta_bulk_modulus_pa).sum();
            if !positive(target) || (sum - target).abs() > 2e-10 * target {
                return Err(state_error(
                    "molecular bulk memory lost its static/frozen limits",
                ));
            }
        }
        Ok(BulkRelaxation {
            source: self.clone(),
            equilibrium_bulk_modulus_pa: equilibrium,
            frozen_bulk_modulus_pa: frozen,
            arms,
        })
    }
}

impl BulkRelaxation {
    /// Unchanged species, rates, heat capacities and equilibrium gas state.
    pub const fn source(&self) -> &RelaxationState {
        &self.source
    }
    /// Locally equilibrated adiabatic modulus, not isothermal wall compression [Pa].
    pub const fn equilibrium_bulk_modulus_pa(&self) -> f64 {
        self.equilibrium_bulk_modulus_pa
    }
    /// Instantaneous modulus with all molecular reservoirs frozen [Pa].
    pub const fn frozen_bulk_modulus_pa(&self) -> f64 {
        self.frozen_bulk_modulus_pa
    }
    /// Positive branches, in ascending decay-rate order.
    pub fn arms(&self) -> &[BulkRelaxationArm] {
        &self.arms
    }

    /// Complex bulk modulus in the original `exp(-i omega t)` convention.
    pub fn response(&self, omega_rad_s: f64) -> Result<C64, MaterialError> {
        if !omega_rad_s.is_finite() || omega_rad_s < 0.0 {
            return Err(state_error(
                "molecular bulk response requires finite nonnegative frequency",
            ));
        }
        let mut value = C64::from_re(self.equilibrium_bulk_modulus_pa);
        for arm in &self.arms {
            // s/(s+1/tau), s=-i*omega. Using reciprocal frequency in the
            // original Debye evaluator retains the small real increment and
            // handles zero/infinite omega*tau without subtractive cancellation.
            let (re, im) = debye((omega_rad_s * arm.relaxation_time_s).recip());
            value = value + C64::new(re, -im).scale(arm.delta_bulk_modulus_pa);
        }
        if !positive(value.re) || !value.im.is_finite() || value.im > 0.0 {
            return Err(state_error("molecular bulk response is unrepresentable"));
        }
        Ok(value)
    }
}

#[cfg(test)]
mod tests;
