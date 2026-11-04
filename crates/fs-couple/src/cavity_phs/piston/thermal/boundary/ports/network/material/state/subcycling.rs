//! Whole-system refinement of SEALED gases and the finite material assembly.
//!
//! Every real substep advances the original circuit, thermal cells, reactions
//! and independent piston walls together. The next reaction therefore sees the
//! actual wall heat and compression from the previous substep. No nozzle or
//! species-flow stage is implied; this does not replace `state::flow`.
//! N/2N comparison controls the existing first-order split, without introducing
//! a second physical law or claiming monolithic implicit combustion closure.

use super::super::reactive::ReactiveExchangeError;
use super::*;
use crate::cavity_phs::piston::PistonError;
use crate::cavity_phs::piston::reactive::thermal::subcycling::{
    ReactiveThermalSubcycling, scaled_difference,
};
use crate::cavity_phs::piston::reactive::{GasReactionError, GasReactionSpecies, mechanism::Sum};
use crate::cavity_phs::piston::thermal::ThermalPistonError;
use crate::cavity_phs::piston::thermal::boundary::ports::network::NetworkError;
use crate::circuit_enthalpy::CircuitEnthalpyError;
use crate::circuit_enthalpy::material::MaterialNetworkError;
use crate::port_adoption::{PortTransfer, StepAccount};
use fs_thermochem::ReactionId;

/// Reuse the gas/body scales and work bounds, adding electrical storage scales.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialExchangeSubcycling {
    /// `absolute_body_enthalpy_j_kg` applies separately to every finite cell.
    /// Work counts complete attempted assembly steps, not inner Newton trials.
    pub coupled: ReactiveThermalSubcycling,
    pub absolute_charge_c: f64,
    pub absolute_flux_wb: f64,
}
impl MaterialExchangeSubcycling {
    fn validate<E>(self) -> Result<(), MaterialExchangeSubcyclingError<E>> {
        self.coupled.validate::<()>().map_err(|_| {
            MaterialExchangeSubcyclingError::Invalid(
                "positive independent coupled scales and bounded dyadic work required",
            )
        })?;
        if [self.absolute_charge_c, self.absolute_flux_wb]
            .into_iter()
            .any(|v| !v.is_finite() || v <= 0.0)
        {
            return Err(MaterialExchangeSubcyclingError::Invalid(
                "positive finite charge/flux scales required",
            ));
        }
        Ok(())
    }
}

#[derive(Debug)]
pub enum MaterialExchangeSubcyclingError<E = MaterialNetworkError> {
    Invalid(&'static str),
    Step(ReactiveExchangeError<E>),
    Limit {
        evaluations: usize,
        substeps: usize,
        last_error_ratio: Option<f64>,
    },
    Cancelled,
}
impl<E: core::fmt::Debug> core::fmt::Display for MaterialExchangeSubcyclingError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material assembly refinement: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for MaterialExchangeSubcyclingError<E> {}
impl<E> From<ReactiveExchangeError<E>> for MaterialExchangeSubcyclingError<E> {
    fn from(e: ReactiveExchangeError<E>) -> Self {
        match e {
            ReactiveExchangeError::Cancelled => Self::Cancelled,
            other => Self::Step(other),
        }
    }
}
type Error<E> = MaterialExchangeSubcyclingError<E>;

/// Original accepted subinterval accounts, sources and numerical work.
pub struct MaterialExchangeSubstep<I> {
    pub start_s: f64,
    pub end_s: f64,
    pub advance: MaterialExchangeAdvance<I>,
}

/// Whole-interval chamber account; channel extents are never combined across IDs.
#[derive(Clone, Debug, PartialEq)]
pub struct RefinedChamberStep {
    pub chamber: u64,
    /// Empty for a nonreacting chamber; otherwise the original complete catalog.
    pub species: Vec<GasReactionSpecies>,
    pub reaction_extents_mol: Vec<(ReactionId, f64)>,
    pub maximum_species_residual_mol: f64,
    pub gas_energy_change_j: f64,
    pub kinetic_energy_change_j: f64,
    pub external_force_work_j: f64,
    pub exterior_pressure_work_j: f64,
    pub gas_boundary_heat_j: f64,
    pub gas_boundary_entropy_j_k: f64,
    pub gas_pressure_impulse_n_s: f64,
    pub exterior_pressure_impulse_n_s: f64,
    pub external_force_impulse_n_s: f64,
    pub support_impulse_n_s: f64,
    pub momentum_residual_n_s: f64,
    pub energy_residual_j: f64,
}

/// Actual accepted trajectory and independently reconstructed macro balances.
pub struct MaterialExchangeSubcyclingReport<I = ResolvedMaterialStatePoint> {
    pub start_s: f64,
    pub end_s: f64,
    pub evaluations: usize,
    pub meshes_attempted: usize,
    /// Maximum normalized N-vs-2N discrepancy; an estimate, not a rigorous bound.
    pub estimated_error_ratio: f64,
    pub thermal_energy_change_j: f64,
    pub electrical_energy_change_j: f64,
    pub electrical_supply_work_j: f64,
    pub deposited_joule_j: f64,
    /// Boundary/source heat outside the finite assembly, excluding internal transfers.
    pub external_network_heat_j: f64,
    pub prescribed_gas_heat_j: f64,
    /// Includes signed gas-wall source heat, as in the original network account.
    pub non_electrical_heat_j: f64,
    /// Sum of substep maximum cell residuals; bounds every fixed cell's drift.
    pub maximum_cell_residual_bound_j: f64,
    /// Electrical plus cell endpoint storage minus their actual source integrals.
    pub network_energy_residual_j: f64,
    pub total_energy_residual_j: f64,
    pub chambers: Vec<RefinedChamberStep>,
    /// Only actual accepted fine steps; rejected meshes never debit a source.
    pub steps: Vec<MaterialExchangeSubstep<I>>,
}

impl<'m> MaterialExchangeState<'m> {
    /// Refine the complete sealed-chamber circuit/material-heat assembly.
    ///
    /// Electrical port inputs (efforts/flows in the existing circuit convention)
    /// and chamber exterior pressure/applied force are held. Chamber `heat_j`
    /// is integrated external energy, apportioned at constant power with a final
    /// remainder. Thermal `source_w`, convection, contact heat and radiation keep
    /// their original time integration and source owners. `guess_j` is only a
    /// numerical cell-energy guess and is scaled to the actual subinterval.
    ///
    /// Every accepted cell enthalpy, charge, flux, gas inventory and wall momentum
    /// feeds the next physical substep. Compare complete N/2N trajectories using
    /// separate SI scales; accept the finer real endpoint without extrapolation.
    /// Original macro energy, species and momentum tolerances remain in force.
    /// Source/phase refusal or cancellation publishes no subset of these owners.
    ///
    /// Cells, cell masses, circuit topology and sealed chamber identities remain
    /// fixed. This method does not compose `flow` or infer moving material domains.
    /// It refines the existing first-order split, not a monolithic chemistry solve.
    #[allow(clippy::too_many_arguments)]
    pub fn step_subcycled_candidate(
        &self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact],
        controls: NetworkSolve,
        guess_j: &[f64],
        policy: MaterialExchangeSubcycling,
    ) -> Result<(Self, MaterialExchangeSubcyclingReport), Error<MaterialNetworkError>> {
        refine(
            self,
            cx,
            target_s,
            drives,
            reactions,
            controls,
            guess_j,
            policy,
            |state, target, drive, guess| {
                state.step_candidate(
                    prepared,
                    cx,
                    target,
                    ports,
                    drive,
                    reactions,
                    gas_patches,
                    controls,
                    guess,
                )
            },
        )
    }

    /// As `step_subcycled_candidate`, with the existing source-bearing gas law.
    /// The law sees each actual nonlinear material/gas trial and must be pure.
    /// Natural/forced convection or mixed disjoint patches keep the same owners.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn step_subcycled_with_gas_law_candidate<I, E>(
        &self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>],
        surfaces: &[GasSurface],
        controls: NetworkSolve,
        guess_j: &[f64],
        policy: MaterialExchangeSubcycling,
        resolve: impl Fn(
            Contact,
            &DepositionSimulation<'m>,
            &MaterialExchangeSample,
            &ResolvedThermochemicalGasState,
        ) -> Result<ConductanceSample<I>, E>,
    ) -> Result<(Self, MaterialExchangeSubcyclingReport<I>), Error<E>> {
        refine(
            self,
            cx,
            target_s,
            drives,
            reactions,
            controls,
            guess_j,
            policy,
            |state, target, drive, guess| {
                state.step_with_gas_law_candidate(
                    prepared, cx, target, ports, drive, reactions, surfaces, controls, guess,
                    &resolve,
                )
            },
        )
    }

    /// Publish the complete source-bearing refined endpoint after one final gate.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_subcycled(
        &mut self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact],
        controls: NetworkSolve,
        guess_j: &[f64],
        policy: MaterialExchangeSubcycling,
    ) -> Result<MaterialExchangeSubcyclingReport, Error<MaterialNetworkError>> {
        let (next, report) = self.step_subcycled_candidate(
            prepared,
            cx,
            target_s,
            ports,
            drives,
            reactions,
            gas_patches,
            controls,
            guess_j,
            policy,
        )?;
        check(cx)?;
        *self = next;
        Ok(report)
    }

    /// The same all-owner publication with a caller-selected material gas law.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_subcycled_with_gas_law<I, E>(
        &mut self,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>,
        target_s: f64,
        ports: &[f64],
        drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>],
        surfaces: &[GasSurface],
        controls: NetworkSolve,
        guess_j: &[f64],
        policy: MaterialExchangeSubcycling,
        resolve: impl Fn(
            Contact,
            &DepositionSimulation<'m>,
            &MaterialExchangeSample,
            &ResolvedThermochemicalGasState,
        ) -> Result<ConductanceSample<I>, E>,
    ) -> Result<MaterialExchangeSubcyclingReport<I>, Error<E>> {
        let (next, report) = self.step_subcycled_with_gas_law_candidate(
            prepared, cx, target_s, ports, drives, reactions, surfaces, controls, guess_j, policy,
            resolve,
        )?;
        check(cx)?;
        *self = next;
        Ok(report)
    }
}

fn check<E>(cx: &Cx<'_>) -> Result<(), Error<E>> {
    cx.checkpoint().map_err(|_| Error::Cancelled)
}
fn add<E>(sum: &mut Sum, value: f64) -> Result<(), Error<E>> {
    sum.add(value)
        .map_err(|_| Error::Invalid("assembly accumulation is unrepresentable"))
}
fn chemical<E>(chamber: u64, error: GasReactionError) -> Error<E> {
    ReactiveExchangeError::Reaction { chamber, error }.into()
}

struct Mesh<'m, I> {
    state: MaterialExchangeState<'m>,
    report: MaterialExchangeSubcyclingReport<I>,
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn refine<'m, I, E>(
    initial: &MaterialExchangeState<'m>,
    cx: &Cx<'_>,
    target: f64,
    drives: &[ChamberDrive],
    reactions: &[ChamberReaction<'_>],
    controls: NetworkSolve,
    guess: &[f64],
    policy: MaterialExchangeSubcycling,
    evolve: impl Fn(
        &MaterialExchangeState<'m>,
        f64,
        &[ChamberDrive],
        &[f64],
    ) -> Result<
        (MaterialExchangeState<'m>, MaterialExchangeAdvance<I>),
        ReactiveExchangeError<E>,
    >,
) -> Result<
    (
        MaterialExchangeState<'m>,
        MaterialExchangeSubcyclingReport<I>,
    ),
    Error<E>,
> {
    check(cx)?;
    policy.validate()?;
    let duration = target - initial.time_s();
    if !duration.is_finite()
        || duration <= 0.0
        || drives.len() != initial.chambers.len()
        || initial.chambers.len() > controls.maximum_chambers
        || reactions.len() > initial.chambers.len()
        || guess.len() != initial.thermal.cells().len()
        || initial.thermal.cells().len() > controls.network.maximum_cells
        || drives.iter().any(|d| !d.heat_j.is_finite())
    {
        return Err(Error::Invalid(
            "one forward clock and bounded complete assembly inputs required",
        ));
    }
    let mut coarse = None;
    let mut evaluations = 0;
    let mut meshes = 0;
    let mut count = 1;
    let mut last_ratio = None;
    loop {
        check(cx)?;
        meshes += 1;
        let trial = march(
            initial,
            cx,
            target,
            drives,
            reactions,
            controls,
            guess,
            policy,
            count,
            &mut evaluations,
            last_ratio,
            &evolve,
        );
        match trial {
            Ok(mut fine) => {
                if let Some(previous) = &coarse {
                    let ratio = compare(previous, &fine, policy)?;
                    last_ratio = Some(ratio);
                    if ratio <= 1.0 {
                        fine.report.evaluations = evaluations;
                        fine.report.meshes_attempted = meshes;
                        fine.report.estimated_error_ratio = ratio;
                        macro_ports(initial, &mut fine, cx, drives, duration)?;
                        check(cx)?;
                        return Ok((fine.state, fine.report));
                    }
                }
                fine.report.steps.clear();
                coarse = Some(fine);
            }
            Err(Error::Step(e)) if refinable(&e) => {
                coarse = None;
                last_ratio = None;
            }
            Err(e) => return Err(e),
        }
        if count == policy.coupled.maximum_substeps {
            return Err(Error::Limit {
                evaluations,
                substeps: count,
                last_error_ratio: last_ratio,
            });
        }
        count *= 2;
    }
}

fn refinable<E>(e: &ReactiveExchangeError<E>) -> bool {
    matches!(
        e,
        ReactiveExchangeError::Reaction {
            error: GasReactionError::ConsumptionLimit { .. }
                | GasReactionError::EntropyDecrease { .. }
                | GasReactionError::SubcyclingLimit { .. }
                | GasReactionError::ImplicitOutsideBracket { .. }
                | GasReactionError::ImplicitNotConverged { .. }
                | GasReactionError::ImplicitEvaluationLimit { .. },
            ..
        } | ReactiveExchangeError::Exchange(
            CircuitEnthalpyError::NotConverged(_)
                | CircuitEnthalpyError::CircuitBudget
                | CircuitEnthalpyError::Source(NetworkError::Boundary {
                    error: ThermalPistonError::OutsideBracket { .. }
                        | ThermalPistonError::NotConverged { .. }
                        | ThermalPistonError::Piston(
                            PistonError::TimeResolution(
                                fs_phs::linear_exact::ResolvedBandError::Outside { .. }
                            ) | PistonError::OutsideBracket { .. }
                                | PistonError::NotConverged { .. }
                        ),
                    ..
                })
        )
    )
}

#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn march<'m, I, E>(
    initial: &MaterialExchangeState<'m>,
    cx: &Cx<'_>,
    target: f64,
    drives: &[ChamberDrive],
    reactions: &[ChamberReaction<'_>],
    controls: NetworkSolve,
    guess: &[f64],
    policy: MaterialExchangeSubcycling,
    count: usize,
    evaluations: &mut usize,
    last_ratio: Option<f64>,
    evolve: &impl Fn(
        &MaterialExchangeState<'m>,
        f64,
        &[ChamberDrive],
        &[f64],
    ) -> Result<
        (MaterialExchangeState<'m>, MaterialExchangeAdvance<I>),
        ReactiveExchangeError<E>,
    >,
) -> Result<Mesh<'m, I>, Error<E>> {
    let start = initial.time_s();
    let duration = target - start;
    let mut state = MaterialExchangeState::from_parts(
        initial.thermal.clone(),
        initial.circuit.clone(),
        initial.chambers.clone(),
    );
    let mut given = vec![0.0; drives.len()];
    let mut steps = Vec::with_capacity(count);
    for index in 0..count {
        check(cx)?;
        if *evaluations == policy.coupled.maximum_evaluations {
            return Err(Error::Limit {
                evaluations: *evaluations,
                substeps: count,
                last_error_ratio: last_ratio,
            });
        }
        let fraction = (index + 1) as f64 / count as f64;
        let end = if index + 1 == count {
            target
        } else {
            start + duration * fraction
        };
        let before = state.time_s();
        let dt = end - before;
        if !dt.is_finite() || dt <= 0.0 {
            return Err(Error::Invalid("assembly subinterval is unrepresentable"));
        }
        let source_fraction = if index + 1 == count {
            1.0
        } else {
            (end - start) / duration
        };
        let subdrives = drives
            .iter()
            .enumerate()
            .map(|(i, drive)| {
                let cumulative = drive.heat_j * source_fraction;
                let heat_j = cumulative - given[i];
                given[i] = cumulative;
                ChamberDrive { heat_j, ..*drive }
            })
            .collect::<Vec<_>>();
        let subguess = guess
            .iter()
            .map(|v| v * (dt / duration))
            .collect::<Vec<_>>();
        *evaluations += 1;
        let (next, advance) = evolve(&state, end, &subdrives, &subguess)?;
        // This fixed assembly cannot silently create/remove thermal stores.
        if !next
            .thermal
            .cells()
            .keys()
            .eq(initial.thermal.cells().keys())
            || next
                .thermal
                .cells()
                .iter()
                .any(|(id, c)| c.mass_kg.to_bits() != initial.thermal.cells()[id].mass_kg.to_bits())
            || !next.chambers.keys().eq(initial.chambers.keys())
        {
            return Err(Error::Invalid(
                "refined assembly requires fixed cell masses and chamber topology",
            ));
        }
        state = next;
        steps.push(MaterialExchangeSubstep {
            start_s: before,
            end_s: end,
            advance,
        });
    }
    let report = summarize(initial, &state, cx, drives, reactions, controls, steps)?;
    Ok(Mesh { state, report })
}

#[allow(clippy::too_many_arguments)]
fn summarize<I, E>(
    initial: &MaterialExchangeState<'_>,
    state: &MaterialExchangeState<'_>,
    cx: &Cx<'_>,
    drives: &[ChamberDrive],
    reactions: &[ChamberReaction<'_>],
    controls: NetworkSolve,
    steps: Vec<MaterialExchangeSubstep<I>>,
) -> Result<MaterialExchangeSubcyclingReport<I>, Error<E>> {
    let mut supply = Sum::default();
    let mut joule = Sum::default();
    let mut external = Sum::default();
    let mut thermal_input = Sum::default();
    let mut non_electrical = Sum::default();
    let mut cell_residual_bound = Sum::default();
    for interval in &steps {
        check(cx)?;
        let r = &interval.advance.report;
        add(&mut supply, r.electrical.native.supplied)?;
        add(&mut joule, r.deposited_joule_j)?;
        add(&mut non_electrical, r.non_electrical_heat_j)?;
        add(&mut cell_residual_bound, r.maximum_cell_residual_j)?;
        add(&mut external, interval.advance.step.external_network_heat_j)?;
        add(
            &mut thermal_input,
            r.thermal.external_heat_j + r.thermal.incoming_energy_j,
        )?;
    }
    let dh = state.thermal.energy_j() - initial.thermal.energy_j();
    let first = &steps[0].advance.report.electrical;
    let last = &steps
        .last()
        .expect("nonempty mesh")
        .advance
        .report
        .electrical;
    let de = last.stored_after_j - first.stored_before_j;
    let network_residual = (de + dh) - supply.value - non_electrical.value;
    if !network_residual.is_finite()
        || network_residual.abs() > controls.network.total_energy_tolerance_j
        || cell_residual_bound.value > controls.network.cell_energy_tolerance_j
        || !dh.is_finite()
        || !de.is_finite()
        || (dh - thermal_input.value).abs() > controls.network.thermal_energy_tolerance_j
        || (de + joule.value - supply.value).abs() > controls.network.circuit.energy_tolerance_j
    {
        return Err(Error::Invalid(
            "macro electrical or finite-cell energy exceeds the original allowance",
        ));
    }
    let mut gas_heat = Sum::default();
    let mut total = Sum::default();
    for value in [dh, de, -supply.value, -external.value] {
        add(&mut total, value)?;
    }
    let mut chambers = Vec::with_capacity(initial.chambers.len());
    for (&id, before) in &initial.chambers {
        check(cx)?;
        let after = &state.chambers[&id];
        let drive = drives
            .iter()
            .find(|d| d.chamber == id)
            .expect("original drive admission");
        let binding = reactions.iter().find(|r| r.chamber == id);
        let account = chamber_summary(id, before, after, cx, drive, binding, &steps)?;
        for value in [
            account.gas_energy_change_j,
            account.kinetic_energy_change_j,
            account.exterior_pressure_work_j,
            -account.external_force_work_j,
            -drive.heat_j,
        ] {
            add(&mut total, value)?;
        }
        add(&mut gas_heat, drive.heat_j)?;
        chambers.push(account);
    }
    if total.value.abs() > controls.total_energy_tolerance_j {
        return Err(Error::Invalid(
            "macro circuit/cell/gas/kinetic endpoints do not close supplied energy",
        ));
    }
    Ok(MaterialExchangeSubcyclingReport {
        start_s: initial.time_s(),
        end_s: state.time_s(),
        evaluations: 0,
        meshes_attempted: 0,
        estimated_error_ratio: 0.0,
        thermal_energy_change_j: dh,
        electrical_energy_change_j: de,
        electrical_supply_work_j: supply.value,
        deposited_joule_j: joule.value,
        external_network_heat_j: external.value,
        prescribed_gas_heat_j: gas_heat.value,
        non_electrical_heat_j: non_electrical.value,
        maximum_cell_residual_bound_j: cell_residual_bound.value,
        network_energy_residual_j: network_residual,
        total_energy_residual_j: total.value,
        chambers,
        steps,
    })
}

#[allow(clippy::too_many_arguments)]
fn chamber_summary<I, E>(
    id: u64,
    before: &GasPiston,
    after: &GasPiston,
    cx: &Cx<'_>,
    drive: &ChamberDrive,
    binding: Option<&ChamberReaction<'_>>,
    steps: &[MaterialExchangeSubstep<I>],
) -> Result<RefinedChamberStep, Error<E>> {
    let mut sums = [Sum::default(); 10];
    // force work, exterior work, gas heat, Q/T, gas/exterior/applied/support impulses,
    // isolated chemical energy error and isolated chemical entropy change.
    let mut source_residual = Sum::default();
    let mut changes = binding
        .map(|b| vec![Sum::default(); b.reaction.catalog().thermodynamics().components().len()])
        .unwrap_or_default();
    let mut extents = binding
        .map(|b| vec![Sum::default(); b.reaction.channels().count()])
        .unwrap_or_default();
    let mut ids = Vec::new();
    for interval in steps {
        check(cx)?;
        let chamber = interval
            .advance
            .report
            .source
            .chambers
            .iter()
            .find(|c| c.id == id)
            .ok_or(Error::Invalid("accepted substep omitted a chamber account"))?;
        let p = &chamber.boundary.piston;
        let dt = interval.end_s - interval.start_s;
        let values = [
            p.external_force_work_j,
            p.exterior_pressure_work_j,
            p.gas.supplied_heat_j,
            p.gas.supplied_heat_j / p.midpoint.state().temperature,
            dt * before.area_m2() * p.midpoint.state().pressure,
            dt * before.area_m2() * drive.exterior_pressure_pa,
            dt * drive.external_force_n,
            p.support_impulse_n_s,
        ];
        for (sum, value) in sums.iter_mut().zip(values) {
            add(sum, value)?;
        }
        add(&mut source_residual, chamber.boundary.source_residual_j)?;
        if binding.is_some() {
            let reaction = &interval
                .advance
                .report
                .source
                .reactions
                .iter()
                .find(|r| r.chamber == id)
                .ok_or(Error::Invalid(
                    "accepted substep omitted a declared reaction",
                ))?
                .reaction;
            if ids.is_empty() {
                ids = reaction
                    .reaction_extents()
                    .map(|(id, _)| id.clone())
                    .collect();
            }
            for (sum, (_, extent)) in extents.iter_mut().zip(reaction.reaction_extents()) {
                add(sum, extent)?;
            }
            for (sum, species) in changes.iter_mut().zip(&reaction.species) {
                add(sum, species.stoichiometric_change_mol)?;
            }
            add(&mut sums[8], reaction.stored_energy_change_j)?;
            add(&mut sums[9], reaction.entropy_change_j_k)?;
        }
    }
    let mass = before.gas().mass_kg();
    let du = mass
        * (after.gas().state().specific_internal_energy_j_kg()
            - before.gas().state().specific_internal_energy_j_kg());
    let dk = 0.5
        * before.mass_kg()
        * (after.velocity_m_s() - before.velocity_m_s())
        * (after.velocity_m_s() + before.velocity_m_s());
    let energy = (du + dk) + sums[1].value - sums[0].value - sums[2].value;
    let momentum = before.mass_kg() * (after.velocity_m_s() - before.velocity_m_s())
        - ((sums[4].value - sums[5].value) + sums[6].value)
        - sums[7].value;
    if mass.to_bits() != after.gas().mass_kg().to_bits()
        || ![du, dk, energy, momentum].into_iter().all(f64::is_finite)
        || energy.abs() > drive.piston_solve.energy_tolerance_j
        || (energy + source_residual.value).abs() > drive.exchange.energy_tolerance_j
        || momentum.abs() > drive.piston_solve.momentum_tolerance_n_s
    {
        return Err(Error::Invalid(
            "macro chamber mass, energy or momentum exceeds original bounds",
        ));
    }
    let mut species = Vec::new();
    let mut maximum = 0.0_f64;
    if let Some(binding) = binding {
        let old = binding
            .reaction
            .inventory(before.gas())
            .map_err(|e| chemical(id, e))?;
        let new = binding
            .reaction
            .inventory(after.gas())
            .map_err(|e| chemical(id, e))?;
        for (i, component) in binding
            .reaction
            .catalog()
            .thermodynamics()
            .components()
            .iter()
            .enumerate()
        {
            let residual = (new[i] - old[i]) - changes[i].value;
            if !residual.is_finite() {
                return Err(Error::Invalid("macro species balance is unrepresentable"));
            }
            maximum = maximum.max(residual.abs());
            species.push(GasReactionSpecies {
                species: component.species().clone(),
                initial_mol: old[i],
                final_mol: new[i],
                stoichiometric_change_mol: changes[i].value,
            });
        }
        if maximum > binding.controls.species_tolerance_mol
            || sums[8].value.abs() > mass * binding.controls.gas.energy_tolerance_j_kg
            || sums[9].value < -binding.controls.entropy_tolerance_j_k
        {
            return Err(Error::Invalid(
                "macro species/chemical-energy/entropy exceeds original bounds",
            ));
        }
    }
    Ok(RefinedChamberStep {
        chamber: id,
        species,
        reaction_extents_mol: ids
            .into_iter()
            .zip(extents)
            .map(|(id, sum)| (id, sum.value))
            .collect(),
        maximum_species_residual_mol: maximum,
        gas_energy_change_j: du,
        kinetic_energy_change_j: dk,
        external_force_work_j: sums[0].value,
        exterior_pressure_work_j: sums[1].value,
        gas_boundary_heat_j: sums[2].value,
        gas_boundary_entropy_j_k: sums[3].value,
        gas_pressure_impulse_n_s: sums[4].value,
        exterior_pressure_impulse_n_s: sums[5].value,
        external_force_impulse_n_s: sums[6].value,
        support_impulse_n_s: sums[7].value,
        momentum_residual_n_s: momentum,
        energy_residual_j: energy,
    })
}

fn compare<I, E>(
    a: &Mesh<'_, I>,
    b: &Mesh<'_, I>,
    policy: MaterialExchangeSubcycling,
) -> Result<f64, Error<E>> {
    let p = policy.coupled;
    let mut error = 0.0_f64;
    let mut observe = |a, b, absolute| -> Result<(), Error<E>> {
        let ratio = scaled_difference::<()>(a, b, absolute, p.relative_tolerance)
            .map_err(|_| Error::Invalid("assembly temporal estimate is unrepresentable"))?;
        error = error.max(ratio);
        Ok(())
    };
    for (left, right) in a
        .state
        .thermal
        .cells()
        .values()
        .zip(b.state.thermal.cells().values())
    {
        observe(
            left.specific_h_j_kg,
            right.specific_h_j_kg,
            p.absolute_body_enthalpy_j_kg,
        )?;
    }
    for (&left, &right) in a
        .state
        .circuit
        .charges_c()
        .iter()
        .zip(b.state.circuit.charges_c())
    {
        observe(left, right, policy.absolute_charge_c)?;
    }
    for (&left, &right) in a
        .state
        .circuit
        .fluxes_wb()
        .iter()
        .zip(b.state.circuit.fluxes_wb())
    {
        observe(left, right, policy.absolute_flux_wb)?;
    }
    for (id, left) in &a.state.chambers {
        let right = &b.state.chambers[id];
        let x = left.gas().state().state();
        let y = right.gas().state().state();
        observe(x.temperature, y.temperature, p.absolute_temperature_k)?;
        observe(x.pressure, y.pressure, p.absolute_pressure_pa)?;
        observe(
            left.gas().volume_m3(),
            right.gas().volume_m3(),
            p.absolute_volume_m3,
        )?;
        observe(
            left.velocity_m_s(),
            right.velocity_m_s(),
            p.absolute_velocity_m_s,
        )?;
    }
    for (left, right) in a.report.chambers.iter().zip(&b.report.chambers) {
        for (x, y) in left.species.iter().zip(&right.species) {
            observe(x.final_mol, y.final_mol, p.absolute_amount_mol)?;
        }
        for ((_, x), (_, y)) in left
            .reaction_extents_mol
            .iter()
            .zip(&right.reaction_extents_mol)
        {
            observe(*x, *y, p.absolute_amount_mol)?;
        }
    }
    Ok(error)
}

fn macro_ports<I, E>(
    initial: &MaterialExchangeState<'_>,
    mesh: &mut Mesh<'_, I>,
    cx: &Cx<'_>,
    drives: &[ChamberDrive],
    duration: f64,
) -> Result<(), Error<E>> {
    for account in &mesh.report.chambers {
        check(cx)?;
        let id = account.chamber;
        let old = &initial.chambers[&id];
        let next = mesh
            .state
            .chambers
            .get_mut(&id)
            .expect("same accepted chamber identities");
        let entropy = old.gas().mass_kg()
            * (next
                .gas()
                .state()
                .specific_entropy_j_kg_k()
                .map_err(|e| chemical(id, e.into()))?
                - old
                    .gas()
                    .state()
                    .specific_entropy_j_kg_k()
                    .map_err(|e| chemical(id, e.into()))?);
        let storage = account.gas_energy_change_j + account.kinetic_energy_change_j;
        if !storage.is_finite() || !entropy.is_finite() {
            return Err(Error::Invalid(
                "macro gas port storage/entropy is unrepresentable",
            ));
        }
        let drive = drives
            .iter()
            .find(|d| d.chamber == id)
            .expect("admitted original drive");
        next.ports = old.ports.clone();
        next.ports.accept(
            StepAccount {
                tick: next.ports.next_tick(),
                duration_s: duration,
                storage_change_j: Some(storage),
                loss_j: Some(0.0),
                included_work_j: 0.0,
                entropy_change_j_k: entropy,
                tolerance_j: Some(drive.piston_solve.energy_tolerance_j),
                unavailable: None,
            },
            [
                PortTransfer::work(
                    account.external_force_work_j - account.exterior_pressure_work_j,
                    true,
                ),
                PortTransfer::heat_profile(
                    account.gas_boundary_heat_j,
                    account.gas_boundary_entropy_j_k,
                ),
            ],
        );
    }
    Ok(())
}
