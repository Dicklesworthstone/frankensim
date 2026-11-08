//! Matrix-free derivatives through the original anisotropic Maxwell history.
//!
//! The equilibrium, Maxwell and Kelvin tensor factors are fixed. Each branch's
//! logarithmic relaxation time and every imposed strain endpoint may vary. The
//! original six work coordinates per branch survive, including rank-deficient
//! channels; no modal projection, time merging or factor differentiation occurs.
//! A complete experiment propagates its memory directions or cotangents through
//! every step. This is not an adjoint of a surrounding dynamics/thermal solver.

use super::{
    MAX_BRANCHES, Result, TensorProny, TensorPronyHistory, TensorPronyStep, Voigt, apply, dot,
    finite, invalid, mandel, poll, tensor, transpose_add,
};
use crate::visco::PronyRampSensitivity;
use std::sync::Arc;

/// A direction in the original step inputs, with tensors in tensor-shear Voigt
/// order. Memory entries are the original six Gram work coordinates per branch,
/// in the material's branch order, not physical strains or modal amplitudes.
#[derive(Debug, Clone, Copy)]
pub struct TensorPronyStepDirection<'a> {
    pub previous_strain: Voigt,
    pub previous_memory: &'a [Voigt],
    pub strain: Voigt,
    /// One `d ln(tau_j)` per original branch, at fixed physical duration.
    pub log_relaxation_times: &'a [f64],
}

/// Directional derivatives of the physical step observations. The same shape
/// supplies scalar-objective seeds to [`TensorPronyStepLinearization::pull_back`].
/// These are signed derivatives/cotangents, not nonnegative energy ledgers.
///
/// Stress components use tensor-shear Voigt order. A reverse seed pairs with
/// the SIX STORED COMPONENTS by an ordinary dot product. The work/heat outputs
/// already include the correct Mandel shear weights; do not apply them again.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct TensorPronyStepDerivative {
    pub stress: Voigt,
    pub mean_stress: Voigt,
    pub stored_energy: f64,
    pub work: f64,
    pub dissipated: f64,
}

impl TensorPronyStepDerivative {
    fn finite(&self) -> bool {
        finite(&self.stress)
            && finite(&self.mean_stress)
            && [self.stored_energy, self.work, self.dissipated]
                .iter()
                .all(|v| v.is_finite())
    }
}

/// Cotangents for both physical strain endpoints, paired by ordinary six-entry
/// dot products in tensor-shear Voigt order. In a history adjoint, accumulate
/// `previous_strain` into the preceding step's target-strain cotangent.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TensorPronyStrainPullback {
    pub previous_strain: Voigt,
    pub strain: Voigt,
}

#[derive(Debug, Clone, Copy)]
struct BranchLinearization {
    previous: Voigt,
    ramp: PronyRampSensitivity,
}

/// An accepted exact step and its linear-size derivative tape.
///
/// Immutable material tensors are shared with the original history. Each
/// branch retains only its original incoming work coordinates and prepared
/// scalar ramp sensitivity. Preparation reuses the primal exponentials;
/// forward/reverse actions allocate nothing and cost O(number of branches),
/// without constructing a branch-count-squared Jacobian. A caller can retain
/// tapes or reconstruct them from its own history checkpoints.
#[derive(Debug, Clone)]
pub struct TensorPronyStepLinearization {
    model: Arc<TensorProny>,
    previous_strain: Voigt,
    strain: Voigt,
    delta: Voigt,
    dt_s: f64,
    branches: Vec<BranchLinearization>,
    response: TensorPronyStep,
}

impl TensorPronyHistory {
    /// Advance the SAME primal transaction and retain its exact derivative map.
    /// The response is bit-identical to [`Self::advance`] for the same inputs.
    ///
    /// Derivatives hold every admitted tensor factor and physical duration fixed.
    /// They cover original memory, both strain endpoints and branch log times,
    /// including independently integrated heat. They do not differentiate
    /// factorization, source selection, applicability limits or thermal shifts.
    ///
    /// A nonzero Kelvin tensor requires strictly positive duration for this
    /// full gradient map: a zero-time unchanged strain has no unique derivative
    /// with respect to independently varied endpoints. Ordinary `advance`
    /// retains its existing empty-jump behavior. A purely Maxwell/elastic law
    /// retains its exact zero-time elastic-jump derivative.
    ///
    /// Invalid input, cancellation, or failure to reserve the bounded tape
    /// leaves the entire original history unchanged.
    pub fn advance_linearized(
        &mut self,
        strain: Voigt,
        dt_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TensorPronyStepLinearization> {
        poll(&mut cancelled)?;
        if dt_s == 0.0
            && self
                .model
                .kv_pa_s()
                .is_some_and(|c| c.iter().flatten().any(|v| *v != 0.0))
        {
            return Err(invalid(
                "Kelvin endpoint sensitivities require positive duration",
            ));
        }
        let model = Arc::clone(&self.model);
        let previous_strain = self.strain;
        let next_strain = mandel(&strain);
        let delta = core::array::from_fn(|i| next_strain[i] - previous_strain[i]);
        let mut branches = Vec::new();
        branches
            .try_reserve_exact(self.q.len())
            .map_err(|_| invalid("could not reserve the bounded tensor sensitivity tape"))?;
        let response = self.advance_with_ramps(strain, dt_s, &mut cancelled, |r, q| {
            let ramp = PronyRampSensitivity::new(r);
            branches.push(BranchLinearization { previous: *q, ramp });
            *ramp.ramp()
        })?;
        Ok(TensorPronyStepLinearization {
            model,
            previous_strain,
            strain: next_strain,
            delta,
            dt_s,
            branches,
            response,
        })
    }
}

impl TensorPronyStepLinearization {
    /// The original accepted stress, work, storage and physical heat report.
    pub const fn response(&self) -> &TensorPronyStep {
        &self.response
    }

    /// Exact number of six-component memory entries required in each buffer.
    pub fn branch_count(&self) -> usize {
        self.branches.len()
    }

    /// Apply the full step Jacobian to one input direction. On success,
    /// `next_memory` contains the direction to feed into the following step.
    /// Its physical previous-strain direction is this step's `direction.strain`.
    /// Accumulate returned work/heat derivatives to differentiate total ledgers.
    ///
    /// All slices must contain exactly one entry per original branch. Nonfinite
    /// directions, unrepresentable derivatives and cancellation refuse without
    /// changing `next_memory`; signed heat derivatives are never clipped.
    pub fn push_forward(
        &self,
        direction: TensorPronyStepDirection<'_>,
        next_memory: &mut [Voigt],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TensorPronyStepDerivative> {
        poll(&mut cancelled)?;
        let count = self.branch_count();
        if direction.previous_memory.len() != count
            || direction.log_relaxation_times.len() != count
            || next_memory.len() != count
            || !finite(&direction.previous_strain)
            || !finite(&direction.strain)
            || direction.previous_memory.iter().any(|v| !finite(v))
            || direction
                .log_relaxation_times
                .iter()
                .any(|v| !v.is_finite())
        {
            return Err(invalid(
                "complete finite tensor sensitivity directions required",
            ));
        }
        let previous = mandel(&direction.previous_strain);
        let next = mandel(&direction.strain);
        let delta = core::array::from_fn(|i| next[i] - previous[i]);
        let eq = &self.model.equilibrium_gram;
        let old_eq = apply(eq, &self.previous_strain);
        let delta_eq = apply(eq, &self.delta);
        let next_eq = apply(eq, &self.strain);
        let d_old_eq = apply(eq, &previous);
        let d_delta_eq = apply(eq, &delta);
        let d_next_eq = apply(eq, &next);
        let mean_eq = core::array::from_fn(|i| old_eq[i] + 0.5 * delta_eq[i]);
        let d_mean_eq = core::array::from_fn(|i| d_old_eq[i] + 0.5 * d_delta_eq[i]);
        let mut stress = [0.0; 6];
        transpose_add(eq, &d_next_eq, &mut stress);
        let mut mean_stress = [0.0; 6];
        transpose_add(eq, &d_mean_eq, &mut mean_stress);
        let mut stored_energy = dot(&next_eq, &d_next_eq);
        let mut work = dot(&d_mean_eq, &delta_eq) + dot(&mean_eq, &d_delta_eq);
        let mut dissipated = 0.0;
        let mut memory = [[0.0; 6]; MAX_BRANCHES];
        for (j, branch) in self.branches.iter().enumerate() {
            poll(&mut cancelled)?;
            let rows = self.model.branches[j].resolved_gram_rows(&self.model.factors[j]);
            let increment = apply(&rows, &self.delta);
            let d_increment = apply(&rows, &delta);
            let mut mean = [0.0; 6];
            for i in 0..6 {
                let local = branch.ramp.apply(branch.previous[i], increment[i]);
                let pushed = local.push_forward([
                    direction.previous_memory[j][i],
                    d_increment[i],
                    direction.log_relaxation_times[j],
                ]);
                memory[j][i] = pushed[0];
                mean[i] = pushed[1];
                stored_energy += local.response.0 * pushed[0];
                work += pushed[1] * increment[i] + local.response.1 * d_increment[i];
                dissipated += pushed[2];
            }
            transpose_add(&rows, &memory[j], &mut stress);
            transpose_add(&rows, &mean, &mut mean_stress);
        }
        if let Some(kv) = self
            .model
            .kv_pa_s()
            .filter(|c| c.iter().flatten().any(|v| *v != 0.0))
        {
            let rate = self.delta.map(|v| v / self.dt_s);
            let d_rate = delta.map(|v| v / self.dt_s);
            let viscous = apply(kv, &rate);
            let d_viscous = apply(kv, &d_rate);
            let factor = self
                .model
                .kv_factor
                .as_ref()
                .ok_or_else(|| invalid("missing admitted Kelvin sensitivity factor"))?;
            for i in 0..6 {
                stress[i] += d_viscous[i];
                mean_stress[i] += d_viscous[i];
                let channel = (0..6).fold(0.0, |s, k| s + factor[k][i] * rate[k]);
                let d_increment = (0..6).fold(0.0, |s, k| s + factor[k][i] * delta[k]);
                // Differentiate the original independent heat integral. The
                // fixed duration cancels before multiplication, retaining weak loss.
                dissipated += 2.0 * (channel * d_increment);
            }
            work += dot(&d_viscous, &self.delta) + dot(&viscous, &delta);
        }
        let result = TensorPronyStepDerivative {
            stress: tensor(&stress),
            mean_stress: tensor(&mean_stress),
            stored_energy,
            work,
            dissipated,
        };
        if !result.finite() || memory[..count].iter().any(|v| !finite(v)) {
            return Err(invalid(
                "unrepresentable tensor history directional derivative",
            ));
        }
        poll(&mut cancelled)?;
        next_memory.copy_from_slice(&memory[..count]);
        Ok(result)
    }

    /// Apply the transpose of the SAME Jacobian to observation and outgoing
    /// memory seeds. Return both endpoint-strain cotangents and overwrite the
    /// two caller buffers with incoming-memory and branch-log-time cotangents.
    /// For a history, propagate the memory result backward and sum each branch's
    /// log-time contribution over all steps. Cumulative work/heat goals seed
    /// their respective incremental observations at every step.
    ///
    /// All seeds pair with stored components by ordinary dot products. Factors
    /// and tensor-shear conversions are transposed explicitly. Outputs remain
    /// unchanged on shape/data refusal, overflow or cooperative cancellation.
    pub fn pull_back(
        &self,
        seed: TensorPronyStepDerivative,
        next_memory_seed: &[Voigt],
        previous_memory: &mut [Voigt],
        log_relaxation_times: &mut [f64],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TensorPronyStrainPullback> {
        poll(&mut cancelled)?;
        let count = self.branch_count();
        if next_memory_seed.len() != count
            || previous_memory.len() != count
            || log_relaxation_times.len() != count
            || !seed.finite()
            || next_memory_seed.iter().any(|v| !finite(v))
        {
            return Err(invalid("complete finite tensor sensitivity seeds required"));
        }
        // Physical output stress is Mandel stress / WEIGHTS, so this is the
        // transpose conversion for ordinary packed-component cotangents.
        let stress_seed = tensor(&seed.stress);
        let mean_seed = tensor(&seed.mean_stress);
        let eq = &self.model.equilibrium_gram;
        let old_eq = apply(eq, &self.previous_strain);
        let delta_eq = apply(eq, &self.delta);
        let next_eq = apply(eq, &self.strain);
        let eq_stress_seed = apply(eq, &stress_seed);
        let eq_mean_seed = apply(eq, &mean_seed);
        let mut old_eq_seed = [0.0; 6];
        let mut next_eq_seed = [0.0; 6];
        for i in 0..6 {
            let mean = old_eq[i] + 0.5 * delta_eq[i];
            let mean_bar = eq_mean_seed[i] + seed.work * delta_eq[i];
            let delta_bar = seed.work * mean + 0.5 * mean_bar;
            old_eq_seed[i] = mean_bar - delta_bar;
            next_eq_seed[i] = eq_stress_seed[i] + seed.stored_energy * next_eq[i] + delta_bar;
        }
        let mut previous = [0.0; 6];
        let mut next = [0.0; 6];
        transpose_add(eq, &old_eq_seed, &mut previous);
        transpose_add(eq, &next_eq_seed, &mut next);
        let mut memory = [[0.0; 6]; MAX_BRANCHES];
        let mut times = [0.0; MAX_BRANCHES];
        for (j, branch) in self.branches.iter().enumerate() {
            poll(&mut cancelled)?;
            let rows = self.model.branches[j].resolved_gram_rows(&self.model.factors[j]);
            let increment = apply(&rows, &self.delta);
            let stress_bar = apply(&rows, &stress_seed);
            let mean_bar = apply(&rows, &mean_seed);
            let mut increment_bar = [0.0; 6];
            for i in 0..6 {
                let local = branch.ramp.apply(branch.previous[i], increment[i]);
                let pulled = local.pull_back([
                    next_memory_seed[j][i] + stress_bar[i] + seed.stored_energy * local.response.0,
                    mean_bar[i] + seed.work * increment[i],
                    seed.dissipated,
                ]);
                memory[j][i] = pulled[0];
                increment_bar[i] = pulled[1] + seed.work * local.response.1;
                times[j] += pulled[2];
            }
            let mut delta_bar = [0.0; 6];
            transpose_add(&rows, &increment_bar, &mut delta_bar);
            for i in 0..6 {
                next[i] += delta_bar[i];
                previous[i] -= delta_bar[i];
            }
        }
        if let Some(kv) = self
            .model
            .kv_pa_s()
            .filter(|c| c.iter().flatten().any(|v| *v != 0.0))
        {
            let rate = self.delta.map(|v| v / self.dt_s);
            let viscous = apply(kv, &rate);
            let viscous_seed =
                core::array::from_fn(|i| stress_seed[i] + mean_seed[i] + seed.work * self.delta[i]);
            let mut rate_bar = [0.0; 6];
            transpose_add(kv, &viscous_seed, &mut rate_bar);
            let mut delta_bar: Voigt =
                core::array::from_fn(|i| rate_bar[i] / self.dt_s + seed.work * viscous[i]);
            let factor = self
                .model
                .kv_factor
                .as_ref()
                .ok_or_else(|| invalid("missing admitted Kelvin sensitivity factor"))?;
            for i in 0..6 {
                let channel = (0..6).fold(0.0, |s, k| s + factor[k][i] * rate[k]);
                let channel_bar = 2.0 * (seed.dissipated * channel);
                for k in 0..6 {
                    delta_bar[k] += factor[k][i] * channel_bar;
                }
            }
            for i in 0..6 {
                next[i] += delta_bar[i];
                previous[i] -= delta_bar[i];
            }
        }
        // Input strains become Mandel strains by multiplication by WEIGHTS.
        let result = TensorPronyStrainPullback {
            previous_strain: mandel(&previous),
            strain: mandel(&next),
        };
        if !finite(&result.previous_strain)
            || !finite(&result.strain)
            || memory[..count].iter().any(|v| !finite(v))
            || times[..count].iter().any(|v| !v.is_finite())
        {
            return Err(invalid("unrepresentable tensor history cotangent"));
        }
        poll(&mut cancelled)?;
        previous_memory.copy_from_slice(&memory[..count]);
        log_relaxation_times.copy_from_slice(&times[..count]);
        Ok(result)
    }
}
