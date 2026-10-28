//! Error-controlled constant-pressure chemistry with prescribed external heat.
//!
//! Replay complete N and 2N trajectories through the original coupled implicit
//! endpoint. Accept only the finer physical trajectory, never extrapolated
//! species, temperature or volume. Heat is a constant power over the requested
//! interval; formation energy remains in the original total enthalpy chart.
//! This is dyadic refinement of backward Euler, not a higher-order integrator
//! or a rigorous error bound for arbitrary stiff mechanisms.

/// Adaptive reciprocal exchange with a finite original enthalpy/phase body.
pub mod body;

use super::{ClosedGasParcel, GasReaction, GasReactionError, GasReactionSolve,
    IsobaricGasReactionStep};
use super::super::ImplicitMechanismSolve;
use crate::cavity_phs::piston::reactive::{GasReactionSpecies, mechanism::Sum,
    poll, positive, subcycling::{MAX_REACTION_EVALUATIONS, MAX_REACTION_SUBSTEPS}};
use fs_thermochem::ReactionId;
use std::cell::Cell;

/// Temporal accuracy scales and total work, including discarded trajectories.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsobaricAdvanceSolve {
    /// Applied separately to every species amount and every channel extent.
    pub absolute_amount_mol: f64,
    pub absolute_temperature_k: f64,
    pub absolute_volume_m3: f64,
    pub relative_tolerance: f64,
    /// Power of two in 2..=256. At least two complete meshes must agree.
    pub maximum_substeps: usize,
    /// All attempted original endpoints, including refused trials, in 1..=511.
    pub maximum_endpoint_evaluations: usize,
    /// All original residual/Jacobian evaluations across every attempted mesh,
    /// including failed endpoints, in 1..=1_048_576. Each endpoint also keeps
    /// the original implicit mechanism's (possibly smaller) work limit.
    pub maximum_chemical_evaluations: usize,
}
impl IsobaricAdvanceSolve {
    fn validate(self) -> Result<(), IsobaricAdvanceError> {
        if ![self.absolute_amount_mol, self.absolute_temperature_k,
            self.absolute_volume_m3].into_iter().all(positive)
            || !self.relative_tolerance.is_finite() || self.relative_tolerance < 0.0
            || !(2..=MAX_REACTION_SUBSTEPS).contains(&self.maximum_substeps)
            || !self.maximum_substeps.is_power_of_two()
            || !(1..=MAX_REACTION_EVALUATIONS).contains(&self.maximum_endpoint_evaluations)
            || !(1..=1_048_576).contains(&self.maximum_chemical_evaluations)
        {
            return Err(invalid("positive isobaric accuracy scales and bounded dyadic work required"));
        }
        Ok(())
    }
}

/// No refusal publishes a partial gas state or advances a caller-owned clock.
#[derive(Clone, Debug, PartialEq)]
pub enum IsobaricAdvanceError {
    Invalid(&'static str),
    Reaction(GasReactionError),
    EndpointEvaluationLimit { evaluations: usize, maximum: usize },
    ChemicalEvaluationLimit { evaluations: usize, maximum: usize, endpoint_evaluations: usize },
    RefinementLimit {
        substeps: usize,
        endpoint_evaluations: usize,
        chemical_evaluations: usize,
        last_error_ratio: Option<f64>,
    },
    Cancelled,
}
impl From<GasReactionError> for IsobaricAdvanceError {
    fn from(error: GasReactionError) -> Self {
        if error == GasReactionError::Cancelled { Self::Cancelled }
        else { Self::Reaction(error) }
    }
}
impl core::fmt::Display for IsobaricAdvanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "adaptive constant-pressure chemistry: {self:?}")
    }
}
impl std::error::Error for IsobaricAdvanceError {}
fn invalid(message: &'static str) -> IsobaricAdvanceError {
    IsobaricAdvanceError::Invalid(message)
}

/// One actual accepted interval and its complete original kinetic report.
#[derive(Clone, Debug, PartialEq)]
pub struct IsobaricAdvanceSubstep {
    pub start_s: f64,
    pub end_s: f64,
    pub chemistry: IsobaricGasReactionStep,
}

/// Whole-window balances and accepted history. An aggregate extent is not the
/// root of one backward-Euler equation; those roots remain in `steps` only.
#[derive(Clone, Debug, PartialEq)]
pub struct IsobaricAdvanceReport<I = IsobaricAdvanceSubstep> {
    pub duration_s: f64,
    pub pressure_pa: f64,
    pub endpoint_evaluations: usize,
    pub chemical_evaluations: usize,
    pub meshes_attempted: usize,
    pub refused_meshes: usize,
    /// Maximum normalized N-vs-2N discrepancy over species, channel extents,
    /// temperature and volume. This is an asymptotic first-order estimate.
    pub estimated_error_ratio: f64,
    pub species: Vec<GasReactionSpecies>,
    pub reaction_extents_mol: Vec<(ReactionId, f64)>,
    pub maximum_species_residual_mol: f64,
    pub maximum_consumed_fraction: f64,
    pub mass_residual_kg: f64,
    pub supplied_heat_j: f64,
    /// Compensated sum of accepted original formation-inclusive increments.
    /// No new subtraction of rounded absolute endpoint energies is introduced.
    pub enthalpy_change_j: f64,
    pub stored_energy_change_j: f64,
    /// Original fixed pressure times the actual whole-window volume change.
    pub work_by_gas_j: f64,
    pub enthalpy_residual_j: f64,
    pub energy_residual_j: f64,
    pub entropy_change_j_k: f64,
    /// Sum Q_k/T_k using each accepted endpoint, never total Q/final T.
    pub endpoint_heat_entropy_j_k: f64,
    /// Delta S minus that accepted heat entropy. Includes backward-Euler
    /// discretization entropy; an external reservoir still owns its entropy.
    pub entropy_residual_j_k: f64,
    /// Complete original reports, retained once for the selected physical owner.
    pub steps: Vec<I>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct IsobaricAdvanceCandidate {
    gas: ClosedGasParcel,
    report: IsobaricAdvanceReport,
}
impl IsobaricAdvanceCandidate {
    pub const fn gas(&self) -> &ClosedGasParcel { &self.gas }
    pub const fn report(&self) -> &IsobaricAdvanceReport { &self.report }
    pub fn into_parts(self) -> (ClosedGasParcel, IsobaricAdvanceReport) {
        (self.gas, self.report)
    }
}

#[derive(Default)]
struct Work { endpoints: usize, chemical: usize }

impl GasReaction {
    /// Evolve the original closed gas against its initial constant pressure,
    /// with signed total external heat supplied uniformly over `duration_s`.
    ///
    /// Requires the existing coupled implicit mechanism. Every substep solves
    /// composition, enthalpy, temperature and EOS volume together; its actual
    /// final volume enters every rate equation. No extra chemical heat source,
    /// concentration floor or extrapolated endpoint is introduced. The original
    /// source and inventory checks remain authoritative.
    ///
    /// The finer complete trajectory is returned only after temporal comparison
    /// and whole-window mass/species, H=H0+Q, U+p Delta V=U0+Q and heat-entropy
    /// checks. Those absolute conservation budgets are never multiplied by the
    /// substep count. All rejected work consumes the explicit shared budgets.
    /// Consumption, entropy and nonlinear convergence refusals may refine;
    /// source-domain, arithmetic and cancellation failures are final.
    /// This remains homogeneous ideal-gas evolution under a quasistatic load,
    /// not an open-flow reactor or a guarantee of stiff-independent accuracy.
    #[allow(clippy::too_many_arguments)]
    pub fn isobaric_adaptive_candidate(
        &self, gas: &ClosedGasParcel, duration_s: f64, heat_j: f64,
        controls: GasReactionSolve, accuracy: IsobaricAdvanceSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IsobaricAdvanceCandidate, IsobaricAdvanceError> {
        poll(&mut cancelled)?;
        accuracy.validate()?;
        let numerical = self.isobaric_policy(duration_s, heat_j, controls, &mut cancelled)?;
        let mut work = Work::default();
        let mut coarse: Option<IsobaricAdvanceCandidate> = None;
        let (mut count, mut meshes, mut refused) = (1, 0, 0);
        let mut last_ratio = None;
        loop {
            poll(&mut cancelled)?;
            meshes += 1;
            match march(self, gas, duration_s, heat_j, count, controls, accuracy,
                numerical, &mut work, &mut cancelled)
            {
                Ok(mut fine) => {
                    if let Some(previous) = coarse.as_ref() {
                        let ratio = compare(&previous.gas, &previous.report,
                            &fine.gas, &fine.report, accuracy)?;
                        last_ratio = Some(ratio);
                        if ratio <= 1.0 {
                            fine.report.endpoint_evaluations = work.endpoints;
                            fine.report.chemical_evaluations = work.chemical;
                            fine.report.meshes_attempted = meshes;
                            fine.report.refused_meshes = refused;
                            fine.report.estimated_error_ratio = ratio;
                            poll(&mut cancelled)?;
                            return Ok(fine);
                        }
                    }
                    // Keep only comparison observables from the discarded mesh.
                    fine.report.steps = Vec::new();
                    coarse = Some(fine);
                }
                Err(IsobaricAdvanceError::Reaction(GasReactionError::ConsumptionLimit { .. }
                    | GasReactionError::EntropyDecrease { .. }
                    | GasReactionError::ImplicitOutsideBracket { .. }
                    | GasReactionError::ImplicitNotConverged { .. })) => {
                    coarse = None;
                    last_ratio = None;
                    refused += 1;
                }
                Err(error) => return Err(error),
            }
            if count == accuracy.maximum_substeps {
                return Err(IsobaricAdvanceError::RefinementLimit {
                    substeps: count, endpoint_evaluations: work.endpoints,
                    chemical_evaluations: work.chemical, last_error_ratio: last_ratio,
                });
            }
            count *= 2;
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn march(
    law: &GasReaction, initial: &ClosedGasParcel, duration: f64, heat: f64,
    count: usize, controls: GasReactionSolve, accuracy: IsobaricAdvanceSolve,
    original: &ImplicitMechanismSolve, work: &mut Work,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<IsobaricAdvanceCandidate, IsobaricAdvanceError> {
    let mut numerical = original.clone();
    let mut gas = initial.clone();
    let mut steps = Vec::with_capacity(count);
    let (mut start, mut heat_given) = (0.0, 0.0);
    for index in 0..count {
        poll(cancelled)?;
        if work.endpoints == accuracy.maximum_endpoint_evaluations {
            return Err(IsobaricAdvanceError::EndpointEvaluationLimit {
                evaluations: work.endpoints, maximum: accuracy.maximum_endpoint_evaluations,
            });
        }
        let remaining = accuracy.maximum_chemical_evaluations - work.chemical;
        if remaining == 0 {
            return Err(IsobaricAdvanceError::ChemicalEvaluationLimit {
                evaluations: work.chemical, maximum: accuracy.maximum_chemical_evaluations,
                endpoint_evaluations: work.endpoints,
            });
        }
        let fraction = (index + 1) as f64 / count as f64;
        let (end, heat_end) = (duration * fraction, heat * fraction);
        let (dt, supplied) = (end - start, heat_end - heat_given);
        if !positive(dt) || !supplied.is_finite() || (heat != 0.0 && supplied == 0.0) {
            return Err(invalid("isobaric time interval or prescribed heat is unrepresentable"));
        }
        numerical.maximum_evaluations = original.maximum_evaluations.min(remaining);
        let spent = Cell::new(0);
        work.endpoints += 1;
        let proposal = law.isobaric_candidate_counted(&gas, dt, supplied, controls,
            &numerical, cancelled, &spent);
        // Charge refused roots and source failures before deciding whether any
        // retry is allowed. The original owner checks before each evaluation.
        work.chemical += spent.get();
        let proposal = proposal.map_err(|error| match error {
            GasReactionError::ImplicitEvaluationLimit { .. }
                if remaining <= original.maximum_evaluations =>
                IsobaricAdvanceError::ChemicalEvaluationLimit {
                    evaluations: work.chemical, maximum: accuracy.maximum_chemical_evaluations,
                    endpoint_evaluations: work.endpoints,
                },
            other => other.into(),
        })?;
        let (next, chemistry) = proposal.into_parts();
        gas = next;
        steps.push(IsobaricAdvanceSubstep { start_s: start, end_s: end, chemistry });
        (start, heat_given) = (end, heat_end);
    }
    let (gas, report) = summarize(law, initial, gas, duration, heat, controls, steps,
        |interval| &interval.chemistry, cancelled)?;
    Ok(IsobaricAdvanceCandidate { gas, report })
}

#[allow(clippy::too_many_arguments)]
fn summarize<I>(
    law: &GasReaction, initial: &ClosedGasParcel, gas: ClosedGasParcel,
    duration: f64, heat: f64, controls: GasReactionSolve,
    steps: Vec<I>, chemistry: impl Fn(&I) -> &IsobaricGasReactionStep,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(ClosedGasParcel, IsobaricAdvanceReport<I>), IsobaricAdvanceError> {
    let before = law.inventory(initial)?;
    let after = law.inventory(&gas)?;
    let mut extents = vec![Sum::default(); law.channels().count()];
    let mut changes = vec![Sum::default(); before.len()];
    let (mut supplied, mut heat_entropy, mut work_sum) =
        (Sum::default(), Sum::default(), Sum::default());
    let (mut enthalpy_sum, mut energy_sum, mut entropy_sum) =
        (Sum::default(), Sum::default(), Sum::default());
    let (mut enthalpy_residual, mut energy_residual) = (Sum::default(), Sum::default());
    let mut consumed = 0.0_f64;
    for interval in &steps {
        poll(cancelled)?;
        let interval = chemistry(interval);
        let boundary = &interval.boundary;
        for (sum, (_, extent)) in extents.iter_mut().zip(boundary.reaction.reaction_extents()) {
            sum.add(extent)?;
        }
        for (sum, row) in changes.iter_mut().zip(&boundary.reaction.species) {
            sum.add(row.stoichiometric_change_mol)?;
        }
        supplied.add(boundary.boundary.heat_j)?;
        heat_entropy.add(boundary.endpoint_heat_entropy_j_k)?;
        work_sum.add(boundary.boundary.work_by_gas_j)?;
        enthalpy_sum.add(interval.enthalpy_change_j)?;
        energy_sum.add(boundary.reaction.stored_energy_change_j)?;
        entropy_sum.add(boundary.reaction.entropy_change_j_k)?;
        enthalpy_residual.add(interval.enthalpy_residual_j)?;
        energy_residual.add(boundary.energy_residual_j)?;
        consumed = consumed.max(boundary.reaction.maximum_consumed_fraction);
    }
    let mut maximum = 0.0_f64;
    let mut species = Vec::with_capacity(before.len());
    for (i, component) in law.catalog.thermodynamics().components().iter().enumerate() {
        poll(cancelled)?;
        let residual = (after[i] - before[i]) - changes[i].value;
        if !residual.is_finite() { return Err(invalid("isobaric species balance overflowed")); }
        maximum = maximum.max(residual.abs());
        species.push(GasReactionSpecies { species: component.species().clone(),
            initial_mol: before[i], final_mol: after[i],
            stoichiometric_change_mol: changes[i].value });
    }
    let mass = initial.mass_kg();
    let pressure = initial.state().state().pressure;
    // Every original endpoint already checked its source-resolved storage.
    // These accepted increments telescope; discarded paths contribute nothing.
    let (enthalpy, energy, entropy) = (enthalpy_sum.value, energy_sum.value, entropy_sum.value);
    let work = pressure * (gas.volume_m3() - initial.volume_m3());
    let (h_residual, u_residual) = (enthalpy - heat, energy + work - heat);
    let s_residual = entropy - heat_entropy.value;
    let tolerance = mass * controls.gas.energy_tolerance_j_kg;
    if gas.mass_kg().to_bits() != mass.to_bits()
        || gas.state().state().pressure.to_bits() != pressure.to_bits()
        || maximum > controls.species_tolerance_mol || !positive(tolerance)
        || ![enthalpy, energy, entropy, work, h_residual, u_residual, s_residual]
            .into_iter().all(f64::is_finite)
        || [h_residual, u_residual, enthalpy_residual.value, energy_residual.value,
            supplied.value - heat, work_sum.value - work]
            .into_iter().any(|residual| residual.abs() > tolerance)
    {
        return Err(invalid("whole isobaric interval failed original mass/species/enthalpy/work bounds"));
    }
    if s_residual < -controls.entropy_tolerance_j_k {
        return Err(GasReactionError::EntropyDecrease {
            change_j_k: s_residual, tolerance_j_k: controls.entropy_tolerance_j_k,
        }.into());
    }
    poll(cancelled)?;
    Ok((gas, IsobaricAdvanceReport {
        duration_s: duration, pressure_pa: pressure, endpoint_evaluations: 0,
        chemical_evaluations: 0, meshes_attempted: 0, refused_meshes: 0,
        estimated_error_ratio: 0.0, species,
        reaction_extents_mol: law.channels().zip(extents).map(|((kinetics, _), extent)|
            (kinetics.reverse_rate_closure().equilibrium().reaction().clone(), extent.value)).collect(),
        maximum_species_residual_mol: maximum, maximum_consumed_fraction: consumed,
        mass_residual_kg: 0.0, supplied_heat_j: heat, enthalpy_change_j: enthalpy,
        stored_energy_change_j: energy, work_by_gas_j: work,
        enthalpy_residual_j: h_residual, energy_residual_j: u_residual,
        entropy_change_j_k: entropy, endpoint_heat_entropy_j_k: heat_entropy.value,
        entropy_residual_j_k: s_residual, steps,
    }))
}

fn scaled(a: f64, b: f64, absolute: f64, relative: f64) -> Result<f64, IsobaricAdvanceError> {
    let scale = absolute + relative * a.abs().max(b.abs());
    let ratio = (a - b).abs() / scale;
    if !positive(scale) || !ratio.is_finite() {
        return Err(invalid("isobaric temporal error estimate is unrepresentable"));
    }
    Ok(ratio)
}
fn compare<A, B>(a_gas: &ClosedGasParcel, a: &IsobaricAdvanceReport<A>,
    b_gas: &ClosedGasParcel, b: &IsobaricAdvanceReport<B>,
    policy: IsobaricAdvanceSolve,
) -> Result<f64, IsobaricAdvanceError> {
    let mut error = scaled(a_gas.state().state().temperature, b_gas.state().state().temperature,
        policy.absolute_temperature_k, policy.relative_tolerance)?;
    error = error.max(scaled(a_gas.volume_m3(), b_gas.volume_m3(),
        policy.absolute_volume_m3, policy.relative_tolerance)?);
    for (left, right) in a.species.iter().zip(&b.species) {
        error = error.max(scaled(left.final_mol, right.final_mol,
            policy.absolute_amount_mol, policy.relative_tolerance)?);
    }
    // Reaction cycles may hide substantial progress behind unchanged species.
    for ((_, left), (_, right)) in a.reaction_extents_mol.iter().zip(&b.reaction_extents_mol) {
        error = error.max(scaled(*left, *right, policy.absolute_amount_mol, policy.relative_tolerance)?);
    }
    Ok(error)
}
