//! Positive pressure-driven memories of the original molecular gas law.
//!
//! With Cp(s)=Cp_inf+sum C_j/(1+s tau_j), beta(s)=(1-R/Cp(s))/p0.
//! The matrix inversion lemma gives beta=beta_inf+sum d_beta/(1+s*tau).
//! These pressure-driven times differ from both caloric and compression-driven
//! times. They are an exact finite-reservoir realization, not a fitted law.

use super::{RelaxationState, debye, positive, state_error};
use crate::MaterialError;
use fs_la::eigen::jacobi_eigh;
use fs_math::{c64::C64, det};

/// One positive pressure-memory contribution to compressibility.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CompressibilityArm {
    /// Relaxed increment [1/Pa]. Storage is delta_beta*q^2/2 per volume.
    pub delta_compressibility_pa_inv: f64,
    /// Pressure-memory time [s]: tau*q_dot+q=p.
    pub relaxation_time_s: f64,
}
impl CompressibilityArm {
    /// Original peak exp(-i omega t) response, delta_beta/(1-i omega tau).
    pub fn response(&self, omega: f64) -> Result<C64, MaterialError> {
        if !positive(self.delta_compressibility_pa_inv)
            || !positive(self.relaxation_time_s)
            || !omega.is_finite()
            || omega < 0.0
        {
            return Err(state_error(
                "positive finite compressibility strength/time and nonnegative frequency required",
            ));
        }
        let (re, im) = debye(omega * self.relaxation_time_s);
        Ok(C64::new(re, im).scale(self.delta_compressibility_pa_inv))
    }
}

/// Original molecular source and its exact positive pressure-driven realization.
#[derive(Clone, Debug, PartialEq)]
pub struct CompressibilityRelaxation {
    source: RelaxationState,
    frozen: f64,
    equilibrium: f64,
    arms: Vec<CompressibilityArm>,
}

impl RelaxationState {
    /// Eliminate original internal temperatures at prescribed pressure.
    ///
    /// Only exactly equal caloric times are combined. No branch-frequency
    /// window, attenuation fit or extra Stokes--Kirchhoff loss is introduced.
    /// Up to the existing 32 reservoirs enter one bounded symmetric eigensolve.
    pub fn compressibility_relaxation(&self) -> Result<CompressibilityRelaxation, MaterialError> {
        let gas = self.gas();
        gas.validate_coherence()?;
        let r = gas.specific_gas_constant;
        let cp = self.frozen_cv_j_kg_k() + r;
        let frozen = (self.frozen_cv_j_kg_k() / cp) / gas.pressure;
        let equilibrium = self.caloric_response(0.0)?.bulk_modulus_pa.re.recip();
        if !positive(cp) || !positive(frozen) || !positive(equilibrium) {
            return Err(state_error("molecular pressure storage is unrepresentable"));
        }
        let mut groups: Vec<(f64, f64)> = Vec::new();
        for process in self.processes() {
            if process.heat_capacity_j_kg_k == 0.0 {
                continue;
            }
            if let Some((_, capacity)) = groups
                .iter_mut()
                .find(|(tau, _)| *tau == process.relaxation_time_s)
            {
                *capacity += process.heat_capacity_j_kg_k;
            } else {
                groups.push((process.relaxation_time_s, process.heat_capacity_j_kg_k));
            }
        }
        let n = groups.len();
        let mut arms = Vec::with_capacity(n);
        if n != 0 {
            let rates: Vec<_> = groups.iter().map(|(tau, _)| tau.recip()).collect();
            if rates.iter().any(|x| !positive(*x)) {
                return Err(state_error(
                    "molecular pressure-memory rates are unrepresentable",
                ));
            }
            let scale = rates.iter().copied().fold(0.0_f64, f64::max);
            let u: Vec<_> = groups
                .iter()
                .zip(&rates)
                .map(|((_, c), rate)| det::sqrt((c / cp) * (rate / scale)))
                .collect();
            let mut matrix = vec![0.0; n * n];
            for i in 0..n {
                for j in 0..n {
                    matrix[i * n + j] = u[i] * u[j];
                }
                matrix[i * n + i] += rates[i] / scale;
            }
            let norm = matrix.iter().map(|x| x.abs()).fold(0.0_f64, f64::max);
            if !positive(norm) || u.iter().any(|x| !positive(*x)) {
                return Err(state_error(
                    "molecular pressure-memory spectrum is unrepresentable",
                ));
            }
            let (values, vectors) = jacobi_eigh(&matrix, n);
            let tolerance = 128.0 * f64::EPSILON * n as f64;
            for j in 0..n {
                let lambda = values[j];
                if !positive(lambda) {
                    return Err(state_error(
                        "molecular pressure-memory spectrum is unresolved",
                    ));
                }
                for k in 0..=j {
                    let gram: f64 = (0..n)
                        .map(|i| vectors[i * n + j] * vectors[i * n + k])
                        .sum();
                    if (gram - if j == k { 1.0 } else { 0.0 }).abs() > tolerance {
                        return Err(state_error(
                            "molecular pressure-memory eigenbasis is unresolved",
                        ));
                    }
                }
                for i in 0..n {
                    let av: f64 = (0..n).map(|k| matrix[i * n + k] * vectors[k * n + j]).sum();
                    if (av - lambda * vectors[i * n + j]).abs() > tolerance * norm {
                        return Err(state_error(
                            "molecular pressure-memory eigenpair is unresolved",
                        ));
                    }
                }
                let projection: f64 = (0..n).map(|i| u[i] * vectors[i * n + j]).sum();
                if projection == 0.0 {
                    continue;
                }
                let increment = (r / cp / gas.pressure) * (projection * projection / lambda);
                let time = (lambda * scale).recip();
                if !positive(increment) || !positive(time) {
                    return Err(state_error(
                        "molecular pressure-memory arm is unrepresentable",
                    ));
                }
                arms.push(CompressibilityArm {
                    delta_compressibility_pa_inv: increment,
                    relaxation_time_s: time,
                });
            }
            let capacity: f64 = groups.iter().map(|(_, c)| c).sum();
            let expected = (r / cp / gas.pressure) * (capacity / (cp + capacity));
            let actual: f64 = arms.iter().map(|a| a.delta_compressibility_pa_inv).sum();
            if !positive(expected) || (actual - expected).abs() > 2e-10 * expected {
                return Err(state_error(
                    "molecular pressure memory lost its equilibrium/frozen limits",
                ));
            }
        }
        Ok(CompressibilityRelaxation {
            source: self.clone(),
            frozen,
            equilibrium,
            arms,
        })
    }
}
impl CompressibilityRelaxation {
    pub const fn source(&self) -> &RelaxationState {
        &self.source
    }
    pub const fn frozen_compressibility_pa_inv(&self) -> f64 {
        self.frozen
    }
    pub const fn equilibrium_compressibility_pa_inv(&self) -> f64 {
        self.equilibrium
    }
    pub fn arms(&self) -> &[CompressibilityArm] {
        &self.arms
    }
    /// Physical compressibility in peak exp(-i omega t), including dispersion.
    pub fn response(&self, omega: f64) -> Result<C64, MaterialError> {
        if !omega.is_finite() || omega < 0.0 {
            return Err(state_error(
                "finite nonnegative molecular frequency required",
            ));
        }
        let mut response = C64::from_re(self.frozen);
        for arm in &self.arms {
            response = response + arm.response(omega)?;
        }
        if !positive(response.re) || !response.im.is_finite() || response.im < 0.0 {
            return Err(state_error(
                "molecular compressibility response is unrepresentable",
            ));
        }
        Ok(response)
    }
}

#[cfg(test)]
mod tests;
