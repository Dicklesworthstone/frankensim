//! Error-estimated composition of the ORIGINAL conservative chemical step.
//!
//! Compare complete N and 2N step trajectories from one immutable initial gas.
//! Accept the finer actual trajectory, never an extrapolated mixture. This is
//! uniform dyadic refinement of the selected first-order explicit or implicit
//! update, not a second rate law or guarantee of stiffness-independent performance.
//! Chemistry completes once before the existing coupled motion/heat solve.

use super::*;

/// Bounds retained per call, including discarded trial trajectories.
pub const MAX_REACTION_SUBSTEPS: usize = 256;
pub const MAX_REACTION_EVALUATIONS: usize = 2 * MAX_REACTION_SUBSTEPS - 1;

/// Opt-in temporal error estimate and work limits for ONE chemical interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactionSubcycling {
    /// Absolute tolerance for each species amount AND channel extent [mol].
    pub absolute_amount_mol: f64,
    /// Absolute tolerance for the endpoint temperature [K].
    pub absolute_temperature_k: f64,
    /// Absolute tolerance for the endpoint pressure seen by mechanics [Pa].
    pub absolute_pressure_pa: f64,
    /// Nonnegative relative tolerance, applied separately to every observable.
    pub relative_tolerance: f64,
    /// Power of two in 2..=256. The coarse comparison uses half this count.
    pub maximum_substeps: usize,
    /// All attempted original chemical steps, including rejected meshes, 1..=511.
    /// An implicit step additionally retains its own bounded nonlinear work.
    pub maximum_evaluations: usize,
}
impl ReactionSubcycling {
    fn validate(self) -> Result<(), GasReactionError> {
        if !positive(self.absolute_amount_mol) || !positive(self.absolute_temperature_k)
            || !positive(self.absolute_pressure_pa)
            || !self.relative_tolerance.is_finite() || self.relative_tolerance < 0.0
            || !(2..=MAX_REACTION_SUBSTEPS).contains(&self.maximum_substeps)
            || !self.maximum_substeps.is_power_of_two()
            || !(1..=MAX_REACTION_EVALUATIONS).contains(&self.maximum_evaluations)
        {
            return Err(invalid("positive subcycling tolerances and bounded dyadic work required"));
        }
        Ok(())
    }
}

/// Actual accepted subinterval and its original rate/material/stoichiometry report.
/// Times are relative to the caller's chemical interval, not another live clock.
#[derive(Clone, Debug, PartialEq)]
pub struct ReactionSubstep {
    pub start_s: f64,
    pub end_s: f64,
    /// This original single-step report has `subcycling == None`.
    pub reaction: GasReactionStep,
}

/// Only the accepted trajectory is retained; rejected trials change no owner.
#[derive(Clone, Debug, PartialEq)]
pub struct ReactionSubcyclingReport {
    pub duration_s: f64,
    pub evaluations: usize,
    pub meshes_attempted: usize,
    /// Maximum normalized N-vs-2N discrepancy in species, extents, temperature and pressure.
    /// This is an asymptotic first-order error ESTIMATE, not a rigorous bound.
    pub estimated_error_ratio: f64,
    pub steps: Vec<ReactionSubstep>,
}

impl GasReaction {
    /// Select temporal refinement once for all existing parcel, piston, flow,
    /// material-inertia and electrical/thermal reaction consumers. Combine
    /// channels first: there is ONE time policy for the shared inventory.
    /// The default remains the original single explicit chemical step; an
    /// implicit policy, if selected, supplies each original substep instead.
    pub fn with_subcycling(mut self, policy: ReactionSubcycling) -> Result<Self, GasReactionError> {
        policy.validate()?;
        if self.subcycling.is_some() {
            return Err(invalid("reaction already has a subcycling policy"));
        }
        self.subcycling = Some(policy);
        Ok(self)
    }

    #[must_use]
    pub const fn subcycling(&self) -> Option<ReactionSubcycling> { self.subcycling }
}

// Numerical sums only: no independent species or energy storage is introduced.
#[derive(Clone, Copy, Default)]
struct Sum { value: f64, correction: f64 }
impl Sum {
    fn add(&mut self, term: f64) -> Result<(), GasReactionError> {
        let adjusted = term - self.correction;
        let next = self.value + adjusted;
        let correction = (next - self.value) - adjusted;
        if ![term, adjusted, next, correction].iter().all(|v| v.is_finite()) {
            return Err(invalid("chemical trajectory accumulation is unrepresentable"));
        }
        self.value = next;
        self.correction = correction;
        Ok(())
    }
}

struct Trajectory {
    gas: ClosedGasParcel,
    summary: GasReactionStep,
    steps: Vec<ReactionSubstep>,
}

fn limit(evaluations: usize, substeps: usize, ratio: Option<f64>) -> GasReactionError {
    GasReactionError::SubcyclingLimit { evaluations, substeps, last_error_ratio: ratio }
}

/// Replay every candidate mesh from the caller's original state. Only a
/// consumption, entropy-step or implicit-root refusal can request refinement.
/// Source-domain, caloric, arithmetic and cancellation errors remain hard failures.
pub(super) fn candidate(
    law: &GasReaction, initial: &ClosedGasParcel, duration_s: f64,
    controls: GasReactionSolve, policy: ReactionSubcycling,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<GasReactionCandidate, GasReactionError> {
    poll(cancelled)?;
    policy.validate()?;
    if !positive(duration_s) { return Err(invalid("positive finite chemistry interval required")); }
    let mut coarse: Option<Trajectory> = None;
    let mut evaluations = 0;
    let mut meshes = 0;
    let mut count = 1;
    let mut last_ratio = None;
    loop {
        poll(cancelled)?;
        meshes += 1;
        let current = march(law, initial, duration_s, count, controls, policy,
            &mut evaluations, last_ratio, cancelled);
        match current {
            Ok(mut fine) => {
                if let Some(previous) = coarse.as_ref() {
                    let ratio = compare(previous, &fine, policy)?;
                    last_ratio = Some(ratio);
                    if ratio <= 1.0 {
                        let mut step = fine.summary;
                        step.subcycling = Some(ReactionSubcyclingReport {
                            duration_s, evaluations, meshes_attempted: meshes,
                            estimated_error_ratio: ratio, steps: fine.steps,
                        });
                        poll(cancelled)?;
                        return Ok(GasReactionCandidate { gas: fine.gas, step });
                    }
                }
                // Only endpoint/summed observables are needed for comparison.
                // Release discarded source histories before allocating the next mesh.
                fine.steps = Vec::new();
                coarse = Some(fine);
            }
            Err(GasReactionError::ConsumptionLimit { .. } | GasReactionError::EntropyDecrease { .. }
                | GasReactionError::ImplicitOutsideBracket { .. }
                | GasReactionError::ImplicitNotConverged { .. }) => {
                coarse = None;
                last_ratio = None;
            }
            Err(error) => return Err(error),
        }
        if count == policy.maximum_substeps {
            return Err(limit(evaluations, count, last_ratio));
        }
        count *= 2;
    }
}

#[allow(clippy::too_many_arguments)]
fn march(
    law: &GasReaction, initial: &ClosedGasParcel, duration: f64, count: usize,
    controls: GasReactionSolve, policy: ReactionSubcycling, evaluations: &mut usize,
    last_ratio: Option<f64>, cancelled: &mut impl FnMut() -> bool,
) -> Result<Trajectory, GasReactionError> {
    // Split numerical balance allowances, NOT heat, rate, or physical time.
    // This avoids multiplying the admitted energy/species drift by step count.
    let mut inner = controls;
    inner.gas.energy_tolerance_j_kg /= count as f64;
    inner.species_tolerance_mol /= count as f64;
    inner.entropy_tolerance_j_k /= count as f64;
    let mut gas = initial.clone();
    let mut steps = Vec::with_capacity(count);
    let mut start = 0.0;
    for index in 0..count {
        poll(cancelled)?;
        if *evaluations == policy.maximum_evaluations {
            return Err(limit(*evaluations, count, last_ratio));
        }
        // The final boundary is exactly the requested interval. Multiplying a
        // fraction avoids overflow from computing duration * (index + 1) first.
        let end = duration * ((index + 1) as f64 / count as f64);
        let dt = end - start;
        if !positive(dt) { return Err(invalid("chemical subinterval is not representable")); }
        *evaluations += 1;
        let (next, reaction) = law.single_step_candidate(&gas, dt, inner, &mut *cancelled)?.into_parts();
        gas = next;
        steps.push(ReactionSubstep { start_s: start, end_s: end, reaction });
        start = end;
    }
    let summary = summarize(law, initial, &gas, &steps, controls, cancelled)?;
    Ok(Trajectory { gas, summary, steps })
}

fn summarize(
    law: &GasReaction, initial: &ClosedGasParcel, final_gas: &ClosedGasParcel,
    steps: &[ReactionSubstep], controls: GasReactionSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<GasReactionStep, GasReactionError> {
    let mut report = steps[0].reaction.clone();
    // An aggregate extent is not the root of one endpoint equation. Preserve
    // each actual implicit residual only in its corresponding substep.
    report.implicit = None;
    let mut extents = vec![Sum::default(); 1 + report.additional_reactions.len()];
    let mut changes = vec![Sum::default(); report.species.len()];
    let before = law.inventory(initial)?;
    let after = law.inventory(final_gas)?;
    let mut consumed = 0.0_f64;
    for interval in steps {
        poll(cancelled)?;
        for (sum, (_, extent)) in extents.iter_mut().zip(interval.reaction.reaction_extents()) {
            sum.add(extent)?;
        }
        for (sum, species) in changes.iter_mut().zip(&interval.reaction.species) {
            sum.add(species.stoichiometric_change_mol)?;
        }
        consumed = consumed.max(interval.reaction.maximum_consumed_fraction);
    }
    report.extent_mol = extents[0].value;
    for (channel, extent) in report.additional_reactions.iter_mut().zip(&extents[1..]) {
        channel.extent_mol = extent.value;
    }
    let mut maximum = 0.0_f64;
    for (i, species) in report.species.iter_mut().enumerate() {
        poll(cancelled)?;
        species.initial_mol = before[i];
        species.final_mol = after[i];
        species.stoichiometric_change_mol = changes[i].value;
        let residual = (after[i] - before[i]) - changes[i].value;
        if !residual.is_finite() { return Err(invalid("trajectory species balance overflow")); }
        maximum = maximum.max(residual.abs());
    }
    let energy = initial.mass_kg() * (final_gas.state().specific_internal_energy_j_kg()
        - initial.state().specific_internal_energy_j_kg());
    let entropy = initial.mass_kg() * (final_gas.state().specific_entropy_j_kg_k()?
        - initial.state().specific_entropy_j_kg_k()?);
    if final_gas.mass_kg().to_bits() != initial.mass_kg().to_bits()
        || final_gas.volume_m3().to_bits() != initial.volume_m3().to_bits()
        || !energy.is_finite() || energy.abs() > initial.mass_kg() * controls.gas.energy_tolerance_j_kg
        || !entropy.is_finite() || maximum > controls.species_tolerance_mol
    {
        return Err(invalid("full chemical interval failed original mass/species/energy bounds"));
    }
    if entropy < -controls.entropy_tolerance_j_k {
        return Err(GasReactionError::EntropyDecrease {
            change_j_k: entropy, tolerance_j_k: controls.entropy_tolerance_j_k,
        });
    }
    report.maximum_consumed_fraction = consumed;
    report.maximum_species_residual_mol = maximum;
    report.mass_residual_kg = 0.0;
    report.stored_energy_change_j = energy;
    report.entropy_change_j_k = entropy;
    Ok(report)
}

fn scaled_difference(a: f64, b: f64, absolute: f64, relative: f64) -> Result<f64, GasReactionError> {
    let scale = absolute + relative * a.abs().max(b.abs());
    let difference = (a - b).abs();
    let ratio = difference / scale;
    if !positive(scale) || !ratio.is_finite() {
        return Err(invalid("chemical temporal error estimate is unrepresentable"));
    }
    Ok(ratio)
}

fn compare(a: &Trajectory, b: &Trajectory, policy: ReactionSubcycling) -> Result<f64, GasReactionError> {
    let mut error = scaled_difference(a.gas.state().state().temperature, b.gas.state().state().temperature,
        policy.absolute_temperature_k, policy.relative_tolerance)?;
    error = error.max(scaled_difference(a.gas.state().state().pressure, b.gas.state().state().pressure,
        policy.absolute_pressure_pa, policy.relative_tolerance)?);
    for (left, right) in a.summary.species.iter().zip(&b.summary.species) {
        error = error.max(scaled_difference(left.final_mol, right.final_mol,
            policy.absolute_amount_mol, policy.relative_tolerance)?);
    }
    // Cycles can conceal different progress behind nearly identical species.
    for ((_, left), (_, right)) in a.summary.reaction_extents().zip(b.summary.reaction_extents()) {
        error = error.max(scaled_difference(left, right,
            policy.absolute_amount_mol, policy.relative_tolerance)?);
    }
    Ok(error)
}

