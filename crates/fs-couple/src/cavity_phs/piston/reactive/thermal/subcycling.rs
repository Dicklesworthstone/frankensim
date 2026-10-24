//! Error-estimated composition of the original reaction/thermal/motion split.
//!
//! Wall heat and compression change the NEXT chemical substep inside the same
//! requested interval. Compare complete N and 2N trajectories from the immutable
//! initial owners, accepting the finer actual trajectory without extrapolation.
//! This controls the existing first-order split; it is not a monolithic implicit
//! chemical/thermal solve or a guarantee of stiffness-independent performance.

use super::super::subcycling::{MAX_REACTION_EVALUATIONS, MAX_REACTION_SUBSTEPS};
use super::super::{GasReactionSpecies, PistonError, mechanism::Sum, poll, positive};
use super::*;
use crate::port_adoption::{PortTransfer, StepAccount};
use fs_thermochem::ReactionId;

/// Independent physical scales and bounded work for the complete coupled step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReactiveThermalSubcycling {
    /// Absolute tolerance for EACH species amount and EACH channel extent [mol].
    pub absolute_amount_mol: f64,
    pub absolute_temperature_k: f64,
    pub absolute_pressure_pa: f64,
    pub absolute_volume_m3: f64,
    pub absolute_velocity_m_s: f64,
    /// Enthalpy, rather than temperature alone, resolves latent-heat plateaus.
    pub absolute_body_enthalpy_j_kg: f64,
    /// Nonnegative relative tolerance, scaled separately for each observable.
    pub relative_tolerance: f64,
    /// Power of two in 2..=256, as in the original chemical refinement policy.
    pub maximum_substeps: usize,
    /// ALL attempted complete coupled steps, including rejected meshes, 1..=511.
    /// Each original chemistry, heat and motion solve retains its own work cap.
    pub maximum_evaluations: usize,
}
impl ReactiveThermalSubcycling {
    fn validate<E>(self) -> Result<(), ReactiveThermalPistonError<E>> {
        if ![
            self.absolute_amount_mol,
            self.absolute_temperature_k,
            self.absolute_pressure_pa,
            self.absolute_volume_m3,
            self.absolute_velocity_m_s,
            self.absolute_body_enthalpy_j_kg,
        ]
        .into_iter()
        .all(positive)
            || !self.relative_tolerance.is_finite()
            || self.relative_tolerance < 0.0
            || !(2..=MAX_REACTION_SUBSTEPS).contains(&self.maximum_substeps)
            || !self.maximum_substeps.is_power_of_two()
            || !(1..=MAX_REACTION_EVALUATIONS).contains(&self.maximum_evaluations)
        {
            return Err(ReactiveThermalPistonError::Invalid(
                "coupled refinement needs positive absolute scales and bounded dyadic work",
            ));
        }
        Ok(())
    }
}

/// An actual accepted interval, original physical reports and interface source.
#[derive(Clone, Debug, PartialEq)]
pub struct ReactiveThermalSubstep<S> {
    pub start_s: f64,
    pub end_s: f64,
    /// Prescribed external gas heat and constant exterior pressure/force.
    pub input: PistonInput,
    pub body_internal_heat_j: f64,
    /// This step's body endpoint, boundary heat and original source evaluation.
    pub thermal: UniformEnthalpyStep<S>,
    pub step: ReactiveThermalPistonStep,
}

/// Whole-interval accounting, distinct from any one stage or midpoint rate.
#[derive(Clone, Debug, PartialEq)]
pub struct ReactiveThermalSubcyclingReport<S> {
    pub duration_s: f64,
    pub evaluations: usize,
    pub meshes_attempted: usize,
    /// Maximum normalized N-vs-2N discrepancy; an error ESTIMATE, not a bound.
    pub estimated_error_ratio: f64,
    pub species: Vec<GasReactionSpecies>,
    /// Canonical channel identities and their separately accumulated extents.
    pub reaction_extents_mol: Vec<(ReactionId, f64)>,
    pub maximum_species_residual_mol: f64,
    pub gas_energy_change_j: f64,
    pub thermal_energy_change_j: f64,
    pub kinetic_energy_change_j: f64,
    pub exterior_pressure_work_j: f64,
    pub external_force_work_j: f64,
    pub gas_pressure_impulse_n_s: f64,
    pub exterior_pressure_impulse_n_s: f64,
    pub external_force_impulse_n_s: f64,
    pub support_impulse_n_s: f64,
    /// Actual M*(v1-v0) minus accumulated gas/exterior/applied/support impulse.
    pub momentum_residual_n_s: f64,
    pub exchanged_heat_j: f64,
    /// Actual external PLUS wall heat entering the gas over accepted substeps.
    pub gas_boundary_heat_j: f64,
    /// Sum of each actual gas boundary heat divided by its original midpoint T.
    /// Heat reversal need not admit an equivalent single boundary temperature.
    pub gas_boundary_entropy_j_k: f64,
    /// Endpoint body storage plus exchanged heat minus prescribed body heat.
    pub thermal_energy_residual_j: f64,
    /// dU_gas+dH_body+dK+outside_work-force_work-prescribed_gas_heat-body_heat.
    /// Chemical formation energy and reciprocal heat are not extra sources.
    pub total_energy_residual_j: f64,
    /// Only accepted source-bearing steps; rejected trajectories are discarded.
    pub steps: Vec<ReactiveThermalSubstep<S>>,
}

/// Both owners advance together, with a canonical WHOLE-interval body proposal.
#[derive(Clone, Debug, PartialEq)]
pub struct ReactiveThermalSubcyclingCandidate<S, P = GasPiston> {
    piston: P,
    thermal: UniformEnthalpyStep<ReactiveThermalSubcyclingReport<S>>,
}
impl<S, P> ReactiveThermalSubcyclingCandidate<S, P> {
    pub const fn piston(&self) -> &P {
        &self.piston
    }
    /// `external_heat_j` sums ALL accepted wall exchanges, not the final step.
    pub const fn thermal(&self) -> &UniformEnthalpyStep<ReactiveThermalSubcyclingReport<S>> {
        &self.thermal
    }
    pub const fn report(&self) -> &ReactiveThermalSubcyclingReport<S> {
        &self.thermal.report
    }
    pub fn into_parts(self) -> (P, UniformEnthalpyStep<ReactiveThermalSubcyclingReport<S>>) {
        (self.piston, self.thermal)
    }
}
impl<S> ReactiveThermalSubcyclingCandidate<S> {
    pub(in crate::cavity_phs::piston) fn map_piston<P>(
        self,
        bind: impl FnOnce(GasPiston) -> P,
    ) -> ReactiveThermalSubcyclingCandidate<S, P> {
        ReactiveThermalSubcyclingCandidate {
            piston: bind(self.piston),
            thermal: self.thermal,
        }
    }
}

impl GasPiston {
    /// Resolve chemistry, reciprocal finite-body heat and motion repeatedly on
    /// real subintervals, feeding each accepted endpoint into the next reaction.
    ///
    /// Prescribed integrated gas/body heat is distributed at constant power,
    /// with the final remainder retained. Exterior pressure and force are held
    /// constant. Original volume/heat brackets must support every subinterval;
    /// no source law is extrapolated or bracket silently widened. The material
    /// callback must be repeatable and must not debit live external state.
    ///
    /// Species, channel extents, gas T/p/V, wall velocity and body enthalpy each
    /// have an absolute-plus-relative comparison. Macro species and energy
    /// balances retain the ORIGINAL tolerances rather than N times allowances.
    /// An existing reaction integration policy is honored within each coupled
    /// substep; it does not replace this outer splitting-error comparison.
    ///
    /// One macro port account includes actual accepted work, heat and thermal
    /// entropy, and advances the original clock once. Neither discarded meshes
    /// nor source/cancellation refusals publish any owner or boundary transfer.
    #[allow(clippy::too_many_arguments)]
    pub fn reactive_thermal_subcycled_candidate<S, E>(
        &self,
        reaction: &GasReaction,
        input: PistonInput,
        body: UniformEnthalpyStepInput<'_>,
        kinetics: GasReactionSolve,
        controls: PistonSolve,
        exchange: HeatExchangeSolve,
        policy: ReactiveThermalSubcycling,
        mut resolve: impl FnMut(
            &ResolvedThermochemicalGasState,
            EquilibriumPhaseState,
        ) -> Result<ConductanceSample<S>, E>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ReactiveThermalSubcyclingCandidate<S>, ReactiveThermalPistonError<E>> {
        poll(&mut cancelled)?;
        policy.validate()?;
        if !positive(input.duration_s)
            || body.duration_s != input.duration_s
            || !input.heat_j.is_finite()
            || !body.internal_heat_j.is_finite()
        {
            return Err(ReactiveThermalPistonError::Invalid(
                "coupled refinement needs one positive interval and finite integrated heat",
            ));
        }
        let mut coarse = None;
        let mut count = 1;
        let mut evaluations = 0;
        let mut meshes = 0;
        let mut last_ratio = None;
        loop {
            poll(&mut cancelled)?;
            meshes += 1;
            let trial = march(
                self,
                reaction,
                input,
                body,
                kinetics,
                controls,
                exchange,
                policy,
                count,
                &mut evaluations,
                last_ratio,
                &mut resolve,
                &mut cancelled,
            );
            match trial {
                Ok(mut fine) => {
                    if let Some(previous) = coarse.as_ref() {
                        let ratio = compare(previous, &fine, policy)?;
                        last_ratio = Some(ratio);
                        if ratio <= 1.0 {
                            fine.thermal.report.evaluations = evaluations;
                            fine.thermal.report.meshes_attempted = meshes;
                            fine.thermal.report.estimated_error_ratio = ratio;
                            record_account(self, &mut fine, input, exchange.energy_tolerance_j)?;
                            poll(&mut cancelled)?;
                            return Ok(fine);
                        }
                    }
                    // Comparison retains only endpoints and sums, not discarded sources.
                    fine.thermal.report.steps.clear();
                    coarse = Some(fine);
                }
                Err(error) if refinable(&error) => {
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

    /// Use the same material-bound conduction law and enthalpy source at EVERY
    /// accepted/replayed substep, including the evolving gas transport mixture.
    #[allow(clippy::too_many_arguments)]
    pub fn reactive_thermal_layers_subcycled_candidate(
        &self,
        reaction: &GasReaction,
        input: PistonInput,
        body: UniformEnthalpyStepInput<'_>,
        kinetics: GasReactionSolve,
        controls: PistonSolve,
        exchange: HeatExchangeSolve,
        policy: ReactiveThermalSubcycling,
        layers: MaterialConductionLayers<'_>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        ReactiveThermalSubcyclingCandidate<ResolvedMaterialStatePoint>,
        ReactiveThermalPistonError<LayerConductionError>,
    > {
        self.reactive_thermal_subcycled_candidate(
            reaction,
            input,
            body,
            kinetics,
            controls,
            exchange,
            policy,
            |gas, thermal| layers.resolve(gas, thermal),
            cancelled,
        )
    }
}

fn limit<E>(
    evaluations: usize,
    substeps: usize,
    last_error_ratio: Option<f64>,
) -> ReactiveThermalPistonError<E> {
    ReactiveThermalPistonError::SubcyclingLimit {
        evaluations,
        substeps,
        last_error_ratio,
    }
}

// Only timestep/iteration admission can request refinement. Missing source
// support, invalid controls, caloric failures and cancellation remain hard errors.
fn refinable<E>(error: &ReactiveThermalPistonError<E>) -> bool {
    matches!(
        error,
        ReactiveThermalPistonError::Reaction(
            GasReactionError::ConsumptionLimit { .. }
                | GasReactionError::EntropyDecrease { .. }
                | GasReactionError::SubcyclingLimit { .. }
                | GasReactionError::ImplicitOutsideBracket { .. }
                | GasReactionError::ImplicitNotConverged { .. }
                | GasReactionError::ImplicitEvaluationLimit { .. }
        ) | ReactiveThermalPistonError::Thermal(
            ThermalPistonError::OutsideBracket { .. }
                | ThermalPistonError::NotConverged { .. }
                | ThermalPistonError::Piston(
                    PistonError::TimeResolution(
                        fs_phs::linear_exact::ResolvedBandError::Outside { .. }
                    ) | PistonError::OutsideBracket { .. }
                        | PistonError::NotConverged { .. }
                )
        )
    )
}

#[allow(clippy::too_many_arguments)]
fn march<S, E>(
    initial: &GasPiston,
    law: &GasReaction,
    input: PistonInput,
    body: UniformEnthalpyStepInput<'_>,
    kinetics: GasReactionSolve,
    controls: PistonSolve,
    exchange: HeatExchangeSolve,
    policy: ReactiveThermalSubcycling,
    count: usize,
    evaluations: &mut usize,
    last_ratio: Option<f64>,
    resolve: &mut impl FnMut(
        &ResolvedThermochemicalGasState,
        EquilibriumPhaseState,
    ) -> Result<ConductanceSample<S>, E>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ReactiveThermalSubcyclingCandidate<S>, ReactiveThermalPistonError<E>> {
    let mut piston = initial.clone();
    let mut body_state = body.initial;
    let mut steps = Vec::with_capacity(count);
    let (mut start, mut gas_given, mut body_given) = (0.0, 0.0, 0.0);
    // Bound cumulative source residual; other original balances are independently
    // checked at macro endpoints below, without tightening caloric roots to zero.
    let inner_exchange = HeatExchangeSolve {
        source_tolerance_j: exchange.source_tolerance_j / count as f64,
        ..exchange
    };
    for index in 0..count {
        poll(cancelled)?;
        if *evaluations == policy.maximum_evaluations {
            return Err(limit(*evaluations, count, last_ratio));
        }
        let fraction = (index + 1) as f64 / count as f64;
        let end = input.duration_s * fraction;
        let dt = end - start;
        let gas_end = input.heat_j * fraction;
        let body_end = body.internal_heat_j * fraction;
        let subinput = PistonInput {
            duration_s: dt,
            heat_j: gas_end - gas_given,
            ..input
        };
        let subbody = UniformEnthalpyStepInput {
            initial: body_state,
            duration_s: dt,
            internal_heat_j: body_end - body_given,
            ..body
        };
        if !positive(dt) {
            return Err(ReactiveThermalPistonError::Invalid(
                "coupled subinterval is unrepresentable",
            ));
        }
        *evaluations += 1;
        let (next, thermal, step) = piston
            .reactive_thermal_candidate(
                law,
                subinput,
                subbody,
                kinetics,
                controls,
                inner_exchange,
                &mut *resolve,
                &mut *cancelled,
            )?
            .into_parts();
        body_state = thermal.state;
        piston = next;
        steps.push(ReactiveThermalSubstep {
            start_s: start,
            end_s: end,
            input: subinput,
            body_internal_heat_j: subbody.internal_heat_j,
            thermal,
            step,
        });
        (start, gas_given, body_given) = (end, gas_end, body_end);
    }
    summarize(
        initial, piston, law, input, body, kinetics, controls, exchange, steps, cancelled,
    )
}

#[allow(clippy::too_many_arguments)]
fn summarize<S, E>(
    initial: &GasPiston,
    piston: GasPiston,
    law: &GasReaction,
    input: PistonInput,
    body: UniformEnthalpyStepInput<'_>,
    kinetics: GasReactionSolve,
    controls: PistonSolve,
    exchange: HeatExchangeSolve,
    steps: Vec<ReactiveThermalSubstep<S>>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<ReactiveThermalSubcyclingCandidate<S>, ReactiveThermalPistonError<E>> {
    let mut extents = vec![Sum::default(); law.channels().count()];
    let mut species = steps[0].step.reaction.species.clone();
    let mut changes = vec![Sum::default(); species.len()];
    let (
        mut outside,
        mut force,
        mut heat,
        mut gas_heat,
        mut gas_heat_entropy,
        mut chemical_energy,
        mut chemical_entropy,
    ) = (
        Sum::default(),
        Sum::default(),
        Sum::default(),
        Sum::default(),
        Sum::default(),
        Sum::default(),
        Sum::default(),
    );
    let (mut gas_impulse, mut exterior_impulse, mut force_impulse, mut support_impulse) = (
        Sum::default(),
        Sum::default(),
        Sum::default(),
        Sum::default(),
    );
    for interval in &steps {
        poll(cancelled)?;
        let reaction = &interval.step.reaction;
        let thermal = &interval.step.thermal;
        for (sum, (_, extent)) in extents.iter_mut().zip(reaction.reaction_extents()) {
            sum.add(extent)?;
        }
        for (sum, amount) in changes.iter_mut().zip(&reaction.species) {
            sum.add(amount.stoichiometric_change_mol)?;
        }
        outside.add(thermal.piston.exterior_pressure_work_j)?;
        force.add(thermal.piston.external_force_work_j)?;
        let dt = interval.input.duration_s;
        gas_impulse.add(dt * initial.area_m2() * thermal.piston.midpoint.state().pressure)?;
        exterior_impulse.add(dt * initial.area_m2() * interval.input.exterior_pressure_pa)?;
        force_impulse.add(dt * interval.input.external_force_n)?;
        support_impulse.add(thermal.piston.support_impulse_n_s)?;
        heat.add(thermal.exchanged_heat_j)?;
        gas_heat.add(thermal.piston.gas.supplied_heat_j)?;
        gas_heat_entropy.add(
            thermal.piston.gas.supplied_heat_j / thermal.piston.midpoint.state().temperature,
        )?;
        chemical_energy.add(reaction.stored_energy_change_j)?;
        chemical_entropy.add(reaction.entropy_change_j_k)?;
    }
    let before = law.inventory(initial.gas())?;
    let after = law.inventory(piston.gas())?;
    let mut maximum_species_residual = 0.0_f64;
    for (i, amount) in species.iter_mut().enumerate() {
        amount.initial_mol = before[i];
        amount.final_mol = after[i];
        amount.stoichiometric_change_mol = changes[i].value;
        let residual = (after[i] - before[i]) - changes[i].value;
        if !residual.is_finite() {
            return Err(ReactiveThermalPistonError::Invalid(
                "coupled species balance overflow",
            ));
        }
        maximum_species_residual = maximum_species_residual.max(residual.abs());
    }
    let state = steps.last().expect("nonempty admitted mesh").thermal.state;
    let du = initial.gas().mass_kg()
        * (piston.gas().state().specific_internal_energy_j_kg()
            - initial.gas().state().specific_internal_energy_j_kg());
    let dh =
        body.mass_kg * (state.specific_enthalpy_j_kg() - body.initial.specific_enthalpy_j_kg());
    let dk = (0.5 * initial.mass_kg() * (piston.velocity_m_s() - initial.velocity_m_s()))
        * (piston.velocity_m_s() + initial.velocity_m_s());
    let body_residual = (dh + heat.value) - body.internal_heat_j;
    let momentum = initial.mass_kg() * (piston.velocity_m_s() - initial.velocity_m_s())
        - ((gas_impulse.value - exterior_impulse.value) + force_impulse.value)
        - support_impulse.value;
    let total =
        ((du + dh) + dk) + outside.value - force.value - input.heat_j - body.internal_heat_j;
    if ![du, dh, dk, body_residual, momentum, total]
        .into_iter()
        .all(f64::is_finite)
        || piston.gas().mass_kg().to_bits() != initial.gas().mass_kg().to_bits()
        || maximum_species_residual > kinetics.species_tolerance_mol
        || chemical_energy.value.abs()
            > initial.gas().mass_kg() * kinetics.gas.energy_tolerance_j_kg
        || chemical_entropy.value < -kinetics.entropy_tolerance_j_k
        || momentum.abs() > controls.momentum_tolerance_n_s
        || body_residual.abs() > exchange.energy_tolerance_j
        || (gas_heat.value - (input.heat_j + heat.value)).abs() > exchange.energy_tolerance_j
        || total.abs() > exchange.energy_tolerance_j
    {
        return Err(ReactiveThermalPistonError::Invalid(
            "complete coupled interval failed original mass/species/chemical/body/momentum/energy bounds",
        ));
    }
    let reaction_extents_mol = steps[0]
        .step
        .reaction
        .reaction_extents()
        .zip(extents)
        .map(|((id, _), sum)| (id.clone(), sum.value))
        .collect();
    Ok(ReactiveThermalSubcyclingCandidate {
        piston,
        thermal: UniformEnthalpyStep {
            state,
            external_heat_j: -heat.value,
            energy_residual_tolerance_j: exchange.energy_tolerance_j,
            report: ReactiveThermalSubcyclingReport {
                duration_s: input.duration_s,
                evaluations: 0,
                meshes_attempted: 0,
                estimated_error_ratio: 0.0,
                species,
                reaction_extents_mol,
                maximum_species_residual_mol: maximum_species_residual,
                gas_energy_change_j: du,
                thermal_energy_change_j: dh,
                kinetic_energy_change_j: dk,
                exterior_pressure_work_j: outside.value,
                external_force_work_j: force.value,
                gas_pressure_impulse_n_s: gas_impulse.value,
                exterior_pressure_impulse_n_s: exterior_impulse.value,
                external_force_impulse_n_s: force_impulse.value,
                support_impulse_n_s: support_impulse.value,
                momentum_residual_n_s: momentum,
                exchanged_heat_j: heat.value,
                gas_boundary_heat_j: gas_heat.value,
                gas_boundary_entropy_j_k: gas_heat_entropy.value,
                thermal_energy_residual_j: body_residual,
                total_energy_residual_j: total,
                steps,
            },
        },
    })
}

fn compare<S, E>(
    a: &ReactiveThermalSubcyclingCandidate<S>,
    b: &ReactiveThermalSubcyclingCandidate<S>,
    policy: ReactiveThermalSubcycling,
) -> Result<f64, ReactiveThermalPistonError<E>> {
    let mut error = 0.0_f64;
    let mut observe = |left: f64, right: f64, absolute: f64| {
        let scale = absolute + policy.relative_tolerance * left.abs().max(right.abs());
        let ratio = (left - right).abs() / scale;
        if !positive(scale) || !ratio.is_finite() {
            return Err(ReactiveThermalPistonError::Invalid(
                "coupled temporal estimate is unrepresentable",
            ));
        }
        error = error.max(ratio);
        Ok(())
    };
    let ga = a.piston.gas().state().state();
    let gb = b.piston.gas().state().state();
    observe(
        ga.temperature,
        gb.temperature,
        policy.absolute_temperature_k,
    )?;
    observe(ga.pressure, gb.pressure, policy.absolute_pressure_pa)?;
    observe(
        a.piston.gas().volume_m3(),
        b.piston.gas().volume_m3(),
        policy.absolute_volume_m3,
    )?;
    observe(
        a.piston.velocity_m_s(),
        b.piston.velocity_m_s(),
        policy.absolute_velocity_m_s,
    )?;
    observe(
        a.thermal.state.specific_enthalpy_j_kg(),
        b.thermal.state.specific_enthalpy_j_kg(),
        policy.absolute_body_enthalpy_j_kg,
    )?;
    for (left, right) in a.report().species.iter().zip(&b.report().species) {
        observe(left.final_mol, right.final_mol, policy.absolute_amount_mol)?;
    }
    // Reversible cycles can hide different channel progress in equal endpoints.
    for ((_, left), (_, right)) in a
        .report()
        .reaction_extents_mol
        .iter()
        .zip(&b.report().reaction_extents_mol)
    {
        observe(*left, *right, policy.absolute_amount_mol)?;
    }
    Ok(error)
}

fn record_account<S, E>(
    initial: &GasPiston,
    next: &mut ReactiveThermalSubcyclingCandidate<S>,
    input: PistonInput,
    tolerance_j: f64,
) -> Result<(), ReactiveThermalPistonError<E>> {
    let report = next.report();
    let storage = report.gas_energy_change_j + report.kinetic_energy_change_j;
    let entropy = initial.gas().mass_kg()
        * (next
            .piston
            .gas()
            .state()
            .specific_entropy_j_kg_k()
            .map_err(GasReactionError::from)?
            - initial
                .gas()
                .state()
                .specific_entropy_j_kg_k()
                .map_err(GasReactionError::from)?);
    if !storage.is_finite() || !entropy.is_finite() {
        return Err(ReactiveThermalPistonError::Invalid(
            "coupled macro port account is unrepresentable",
        ));
    }
    let transfers = [
        PortTransfer::work(
            report.external_force_work_j - report.exterior_pressure_work_j,
            true,
        ),
        PortTransfer::heat_profile(report.gas_boundary_heat_j, report.gas_boundary_entropy_j_k),
    ];
    // Replace temporary substep accounts; only the original macro clock commits.
    next.piston.ports = initial.ports.clone();
    next.piston.ports.accept(
        StepAccount {
            tick: next.piston.ports.next_tick(),
            duration_s: input.duration_s,
            storage_change_j: Some(storage),
            loss_j: Some(0.0),
            included_work_j: 0.0,
            entropy_change_j_k: entropy,
            tolerance_j: Some(tolerance_j),
            unavailable: None,
        },
        transfers,
    );
    Ok(())
}
