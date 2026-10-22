//! Declared reactions share the electrical/contact/radiative wall network.
//! Reuse GasReaction's conservative fixed-volume stage ONCE, followed by the
//! original simultaneous exchange and motion solve. This is the same first-order
//! reaction-then-motion split as GasPiston::reactive_candidate, not an implicit
//! kinetic solve, ignition model or validated combustion mechanism.
use super::exchange::{PreparedMaterialExchangeNetwork, MaterialExchangeSample, MaterialExchangeError};
use super::MaterialContact;
use super::super::{Chamber, NetworkCandidate, NetworkError, NetworkSolve};
use super::exchange::gas_law::GasSurface;
use crate::cavity_phs::piston::thermal::ConductanceSample;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use crate::circuit_enthalpy::{CircuitEnthalpyError, contacts::ContactThermalSample, material::MaterialNetworkError};
use crate::material_circuit::MaterialCircuit;
use fs_exec::Cx;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_thermal::enthalpy::DepositionSimulation;
use crate::cavity_phs::piston::{GasPiston, reactive::{GasReaction, GasReactionError, GasReactionSolve, GasReactionStep}};
use crate::cavity_phs::piston::thermal::boundary::ports::{BoundaryHeatStep, ThermalPortSample};
use crate::cavity_phs::piston::thermal::boundary::ports::network::{Contact, NetworkStep};
use crate::circuit_enthalpy::CircuitEnthalpyReport;
use crate::port_adoption::{PortTransfer, StepAccount};
use std::collections::BTreeMap;

/// At most one declared reaction per existing chamber. Unlisted chambers remain
/// nonreacting and use their original gas model, including its composition.
#[derive(Clone, Copy, Debug)]
pub struct ChamberReaction<'a> {
    pub chamber: u64,
    pub reaction: &'a GasReaction,
    pub controls: GasReactionSolve,
}

/// Complete chemical stage and extensive gas energy over the WHOLE interval.
/// Formation energy is already included in gas internal energy; no heat credit
/// equal to reaction enthalpy is applied to either the gas or its thermal cells.
#[derive(Clone, Debug, PartialEq)]
pub struct ChamberReactionStep {
    pub chamber: u64,
    pub reaction: GasReactionStep,
    pub gas_energy_change_j: f64,
    /// Whole gas/kinetic change and mechanical work minus prescribed gas heat
    /// and the independently evaluated wall heat, including chemical roundoff.
    pub boundary_energy_residual_j: f64,
}

/// Original boundary result with its gas history's one whole-interval port tick.
/// `boundary` retains its post-reaction pressure/work reference; the reaction
/// report above supplies the chemical stage and whole-interval energy change.
#[derive(Clone, Debug, PartialEq)]
pub struct ReactiveExchangeChamber<I = ResolvedMaterialStatePoint> {
    pub id: u64,
    pub piston: GasPiston,
    pub surfaces: Vec<ThermalPortSample<(Contact, I)>>,
    pub boundary: BoundaryHeatStep,
}

#[derive(Debug)]
pub enum ReactiveExchangeError<E = MaterialNetworkError> {
    Invalid(&'static str),
    Reaction { chamber: u64, error: GasReactionError },
    Exchange(CircuitEnthalpyError<NetworkError<MaterialExchangeError, E>>),
    Cancelled,
}
impl<E: core::fmt::Debug> core::fmt::Display for ReactiveExchangeError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "reactive material exchange: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for ReactiveExchangeError<E> {}

/// Sources and physical endpoints from one accepted interval. The circuit and
/// thermal report's energy reference is unchanged; reactions contribute no net
/// supplied heat. Whole-interval balance is in `ReactiveExchangeCandidate::step`.
pub struct ReactiveExchangeSources<I = ResolvedMaterialStatePoint> {
    pub material: ContactThermalSample<MaterialExchangeSample>,
    pub chambers: Vec<ReactiveExchangeChamber<I>>,
    pub reactions: Vec<ChamberReactionStep>,
    pub transferred_heat_j: f64,
}

/// No input is mutated, including on late chemistry/material/work refusal.
pub struct ReactiveExchangeCandidate<'m, I = ResolvedMaterialStatePoint> {
    thermal: DepositionSimulation<'m>,
    circuit: MaterialCircuit,
    report: CircuitEnthalpyReport<ReactiveExchangeSources<I>>,
    step: NetworkStep,
}
impl<'m, I> ReactiveExchangeCandidate<'m, I> {
    pub fn thermal(&self) -> &DepositionSimulation<'m> { &self.thermal }
    pub fn circuit(&self) -> &MaterialCircuit { &self.circuit }
    pub fn chambers(&self) -> &[ReactiveExchangeChamber<I>] { &self.report.source.chambers }
    pub fn reactions(&self) -> &[ChamberReactionStep] { &self.report.source.reactions }
    pub fn report(&self) -> &CircuitEnthalpyReport<ReactiveExchangeSources<I>> { &self.report }
    pub const fn step(&self) -> &NetworkStep { &self.step }
    /// All energy transfers are already included; publish all these histories
    /// together. Reaction extent must not be applied again on consumption.
    pub fn into_parts(self) -> (DepositionSimulation<'m>, MaterialCircuit,
        CircuitEnthalpyReport<ReactiveExchangeSources<I>>, NetworkStep) {
        (self.thermal, self.circuit, self.report, self.step)
    }
}

fn check<E>(cx: &Cx<'_>) -> Result<(), ReactiveExchangeError<E>> {
    cx.checkpoint().map_err(|_| ReactiveExchangeError::Cancelled)
}

impl PreparedMaterialExchangeNetwork<'_, '_> {
    /// Stage reactions once at accepted composition/T/p, then close all existing
    /// electrical, radiative, conductive and pressure/momentum feedback on those
    /// staged gases. Failed Newton trials never consume additional reactants.
    /// Body heat can affect the NEXT reaction step; rates are not reevaluated
    /// at the final temperature in this explicitly first-order split.
    ///
    /// Each original chamber must have a unique ID AND distinct gas owner. The
    /// same checks apply before staging so cloning cannot hide duplicate owners.
    /// Reaction catalogs, consumption limits, species, mass, internal energy and
    /// entropy are admitted by GasReaction. Reaction/source discontinuities and
    /// excessive consumption require shorter caller-selected intervals.
    #[allow(clippy::too_many_arguments)]
    pub fn reactive_step_candidate<'m>(
        &self, cx: &Cx<'_>, thermal: &DepositionSimulation<'m>, target_s: f64,
        circuit: &MaterialCircuit, ports: &[f64], chambers: &[Chamber<'_>],
        reactions: &[ChamberReaction<'_>], gas_patches: &[MaterialContact],
        controls: NetworkSolve, guess_j: &[f64],
    ) -> Result<ReactiveExchangeCandidate<'m>, ReactiveExchangeError> {
        self.reactive_bound(cx, thermal, target_s, chambers, reactions,
            gas_patches.len(), controls, |staged| self.step_candidate(cx, thermal,
                target_s, circuit, ports, staged, gas_patches, controls, guess_j))
    }

    /// Stage the SAME chemistry once, then solve the complete assembly with a
    /// caller-selected gas heat law. The law sees reacted gas energy midpoints
    /// and exact trial cell material samples. Its typed source payload survives
    /// in every surface report. No second chemistry or whole-energy audit.
    ///
    /// See `step_with_gas_law_candidate` for surface, purity and physical limits.
    /// Source errors retain E; cancellation through cx remains cancellation.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn reactive_with_gas_law_candidate<'m, I, E>(
        &self, cx: &Cx<'_>, thermal: &DepositionSimulation<'m>, target_s: f64,
        circuit: &MaterialCircuit, ports: &[f64], chambers: &[Chamber<'_>],
        reactions: &[ChamberReaction<'_>], surfaces: &[GasSurface],
        controls: NetworkSolve, guess_j: &[f64],
        resolve: impl Fn(Contact, &DepositionSimulation<'m>, &MaterialExchangeSample,
            &ResolvedThermochemicalGasState) -> Result<ConductanceSample<I>, E>,
    ) -> Result<ReactiveExchangeCandidate<'m, I>, ReactiveExchangeError<E>> {
        self.reactive_bound(cx, thermal, target_s, chambers, reactions,
            surfaces.len(), controls, |staged| self.step_with_gas_law_candidate(cx,
                thermal, target_s, circuit, ports, staged, surfaces, controls, guess_j, resolve))
    }

    // All contact laws enter the original staging, balance and one-tick account.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn reactive_bound<'m, I, E>(
        &self, cx: &Cx<'_>, thermal: &DepositionSimulation<'m>, target_s: f64,
        chambers: &[Chamber<'_>], reactions: &[ChamberReaction<'_>],
        surface_count: usize, controls: NetworkSolve,
        evolve: impl FnOnce(&[Chamber<'_>]) -> Result<
            NetworkCandidate<'m, ContactThermalSample<MaterialExchangeSample>, I>,
            CircuitEnthalpyError<NetworkError<MaterialExchangeError, E>>>,
    ) -> Result<ReactiveExchangeCandidate<'m, I>, ReactiveExchangeError<E>> {
        use ReactiveExchangeError as Error;
        check::<E>(cx)?;
        let dt = target_s - thermal.time_s();
        if !dt.is_finite() || dt <= 0.0 || chambers.len() > controls.maximum_chambers
            || reactions.len() > chambers.len() || surface_count > controls.maximum_contacts {
            return Err(Error::Invalid("one positive clock and bounded chamber/reaction/surface lists required"));
        }
        let mut originals = BTreeMap::new();
        for (i, chamber) in chambers.iter().enumerate() {
            check::<E>(cx)?;
            if chamber.input.duration_s != dt || originals.insert(chamber.id, chamber).is_some()
                || chambers[..i].iter().any(|c| std::ptr::eq(c.piston, chamber.piston)) {
                return Err(Error::Invalid("unique original chamber IDs, gas owners and matching clocks required"));
            }
        }
        let mut kinetics = BTreeMap::new();
        for reaction in reactions {
            check::<E>(cx)?;
            if !originals.contains_key(&reaction.chamber)
                || kinetics.insert(reaction.chamber, reaction).is_some() {
                return Err(Error::Invalid("exactly one optional reaction per existing chamber required"));
            }
        }
        // Map order makes reaction scheduling independent of caller slice order.
        let mut staged = BTreeMap::new();
        let mut chemical = BTreeMap::new();
        for (&id, original) in &originals {
            check::<E>(cx)?;
            let mut piston = original.piston.clone();
            if let Some(binding) = kinetics.get(&id) {
                let proposal = binding.reaction.parcel_candidate(piston.gas(), dt,
                    binding.controls, || cx.checkpoint().is_err()).map_err(|error| match error {
                        GasReactionError::Cancelled => Error::Cancelled,
                        other => Error::Reaction { chamber: id, error: other },
                    })?;
                let (gas, report) = proposal.into_parts();
                // Within the piston module: preserve inertia, geometry, velocity,
                // source records and port history. This is not a second step.
                piston.gas = gas;
                chemical.insert(id, report);
            }
            staged.insert(id, piston);
        }
        let staged_chambers = originals.iter().map(|(&id, original)| Chamber {
            id, piston: &staged[&id], input: original.input,
            piston_solve: original.piston_solve, exchange: original.exchange,
        }).collect::<Vec<_>>();
        let proposal = evolve(&staged_chambers).map_err(|error| match error {
                CircuitEnthalpyError::Cancelled => Error::Cancelled,
                other => Error::Exchange(other),
            })?;
        let mut reactions = Vec::with_capacity(chemical.len());
        let mut chemical_energy_j = 0.0;
        for chamber in proposal.chambers() {
            check::<E>(cx)?;
            let Some(reaction) = chemical.get(&chamber.id) else { continue; };
            let original = originals[&chamber.id];
            let gas = chamber.boundary.piston().gas();
            let delta = original.piston.gas().mass_kg()
                * (gas.state().specific_internal_energy_j_kg()
                    - original.piston.gas().state().specific_internal_energy_j_kg());
            let p = &chamber.boundary.step().piston;
            let correction = delta - p.gas.stored_energy_change_j;
            let applied_residual = delta + p.kinetic_energy_change_j + p.exterior_pressure_work_j
                - p.external_force_work_j - original.input.heat_j
                - chamber.boundary.step().exchanged_heat_j;
            let boundary_residual = applied_residual + chamber.boundary.step().source_residual_j;
            chemical_energy_j += correction;
            if !delta.is_finite() || !correction.is_finite() || !chemical_energy_j.is_finite()
                || !applied_residual.is_finite() || !boundary_residual.is_finite()
                || (correction - reaction.stored_energy_change_j).abs() > original.piston_solve.energy_tolerance_j
                || applied_residual.abs() > original.piston_solve.energy_tolerance_j
                || boundary_residual.abs() > original.exchange.energy_tolerance_j {
                return Err(Error::Invalid("whole reacting-gas energy and original wall work do not close"));
            }
            reactions.push(ChamberReactionStep { chamber: chamber.id, reaction: reaction.clone(),
                gas_energy_change_j: delta, boundary_energy_residual_j: boundary_residual });
        }
        let mut step = *proposal.step();
        step.total_energy_residual_j += chemical_energy_j;
        if !step.total_energy_residual_j.is_finite()
            || step.total_energy_residual_j.abs() > controls.total_energy_tolerance_j {
            return Err(Error::Invalid("chemistry, circuit, finite bodies and all gas/motion energy do not close"));
        }
        check::<E>(cx)?;
        let (thermal, circuit, report, _) = proposal.into_parts();
        let mut final_chambers = Vec::with_capacity(report.source.chambers.len());
        for chamber in report.source.chambers {
            check::<E>(cx)?;
            let id = chamber.id;
            let (mut piston, surfaces, boundary) = chamber.boundary.into_parts();
            if let Some(reaction) = reactions.iter().find(|r| r.chamber == id) {
                let original = originals[&id];
                // Replace the inner post-reaction tick with ONE account for the
                // same physical interval, exactly as the native reactive piston.
                piston.ports = original.piston.ports.clone();
                piston.ports.accept(StepAccount {
                    tick: piston.ports.next_tick(), duration_s: dt,
                    storage_change_j: Some(reaction.gas_energy_change_j + boundary.piston.kinetic_energy_change_j),
                    loss_j: Some(0.0), included_work_j: 0.0,
                    entropy_change_j_k: boundary.piston.gas.entropy_change_j_k
                        + reaction.reaction.entropy_change_j_k,
                    tolerance_j: Some(original.piston_solve.energy_tolerance_j), unavailable: None,
                }, [
                    PortTransfer::work(boundary.piston.external_force_work_j - boundary.piston.exterior_pressure_work_j, true),
                    PortTransfer::heat(original.input.heat_j + boundary.exchanged_heat_j, boundary.piston.midpoint.state().temperature),
                ]);
            }
            final_chambers.push(ReactiveExchangeChamber { id, piston, surfaces, boundary });
        }
        let source = ReactiveExchangeSources { material: report.source.material,
            chambers: final_chambers, reactions, transferred_heat_j: report.source.transferred_heat_j };
        let report = CircuitEnthalpyReport { source, electrical: report.electrical,
            thermal: report.thermal, nonlinear: report.nonlinear,
            maximum_cell_residual_j: report.maximum_cell_residual_j, circuit_trials: report.circuit_trials,
            deposited_joule_j: report.deposited_joule_j, non_electrical_heat_j: report.non_electrical_heat_j,
            total_energy_residual_j: report.total_energy_residual_j };
        check::<E>(cx)?;
        Ok(ReactiveExchangeCandidate { thermal, circuit, report, step })
    }
}
