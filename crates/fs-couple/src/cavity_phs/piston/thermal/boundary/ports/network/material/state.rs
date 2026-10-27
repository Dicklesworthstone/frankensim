//! One live owner for all circuit, thermal and gas histories in an assembly.
//! Delegate physics to the existing candidate; move every endpoint together.
/// Conservative nozzle connectivity before the shared reacting/thermal advance.
pub mod flow;
/// Bounded temporal refinement of the complete sealed material assembly.
pub mod subcycling;
/// Continuing same-endpoint chemical, pressure, circuit and finite-wall feedback.
pub mod coupled;

use super::exchange::{MaterialExchangeSample, PreparedMaterialExchangeNetwork, gas_law::GasSurface};
use super::reactive::{ChamberReaction, ChamberReactionStep, ReactiveExchangeError};
use super::MaterialContact;
use super::super::{Chamber, Contact, NetworkSolve, NetworkStep};
use crate::cavity_phs::piston::{GasPiston, PistonInput, PistonSolve};
use crate::cavity_phs::piston::thermal::{ConductanceSample, HeatExchangeSolve};
use crate::cavity_phs::piston::thermal::boundary::ports::{BoundaryHeatStep, ThermalPortSample};
use crate::circuit_enthalpy::{CircuitEnthalpyReport, contacts::ContactThermalSample};
use crate::material_circuit::MaterialCircuit;
use fs_exec::Cx;
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_thermal::enthalpy::DepositionSimulation;
use std::collections::BTreeMap;

/// Loads for one existing chamber over the requested interval. Duration comes
/// from the assembly's thermal clock, not an independently supplied gas clock.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ChamberDrive {
    /// Stable chamber ID; exactly one drive is required per live chamber.
    pub chamber: u64,
    /// Prescribed external heat [J], excluding chemistry and internal wall heat.
    pub heat_j: f64,
    /// Held exterior pressure [Pa].
    pub exterior_pressure_pa: f64,
    /// Held externally applied translating-wall force [N].
    pub external_force_n: f64,
    /// Original gas pressure/work and momentum controls.
    pub piston_solve: PistonSolve,
    /// Original signed wall-heat feedback controls.
    pub exchange: HeatExchangeSolve,
}

/// Original chamber account and sources. The actual endpoint has moved into
/// the assembly, rather than being cloned or left for the caller to install.
pub struct ExchangeChamberStep<I = ResolvedMaterialStatePoint> {
    pub id: u64,
    pub surfaces: Vec<ThermalPortSample<(Contact, I)>>,
    /// Post-reaction gas reference, as in the original boundary owner. Whole-
    /// interval reacting-gas energy is retained in the corresponding reaction.
    pub boundary: BoundaryHeatStep,
}
/// All source records used to advance the actual assembled state.
pub struct ExchangeAdvanceSources<I = ResolvedMaterialStatePoint> {
    pub material: ContactThermalSample<MaterialExchangeSample>,
    pub chambers: Vec<ExchangeChamberStep<I>>,
    pub reactions: Vec<ChamberReactionStep>,
    pub transferred_heat_j: f64,
}
/// Existing physical accounts for the interval, not a second history or solver.
pub struct MaterialExchangeAdvance<I = ResolvedMaterialStatePoint> {
    pub report: CircuitEnthalpyReport<ExchangeAdvanceSources<I>>,
    /// Whole-assembly energy, including chemical and radiation roundoff.
    pub step: NetworkStep,
}

/// Co-owned accepted circuit, thermal timeline and all gases/moving walls.
///
/// The constructor assembles caller-supplied states at ONE intended starting
/// instant; circuit and gas objects have no independent absolute timestamp to
/// infer or reconcile an inconsistent initial epoch. Thereafter `advance` uses
/// the thermal clock and publishes all endpoints in one infallible replacement.
/// Constituents are exposed immutably. Explicitly consume `into_parts` to change
/// topology, source catalogs or geometry outside the fixed-assembly evolution.
/// No parallel species inventory, material registry, or numerical law is added.
pub struct MaterialExchangeState<'m> {
    thermal: DepositionSimulation<'m>,
    circuit: MaterialCircuit,
    chambers: BTreeMap<u64, GasPiston>,
}

// A single clock/identity admission for native, staged-flow and simultaneous
// reacting advances. Physics errors remain with each actual boundary owner.
enum StateAdmissionError {
    Invalid(&'static str),
    Cancelled,
}
impl StateAdmissionError {
    fn staged<E>(self) -> ReactiveExchangeError<E> {
        match self {
            Self::Invalid(message) => ReactiveExchangeError::Invalid(message),
            Self::Cancelled => ReactiveExchangeError::Cancelled,
        }
    }
}

impl<'m> MaterialExchangeState<'m> {
    /// Take ownership of already-admitted initial physical states. The map makes
    /// chamber identities unique; source/area bindings are admitted by preparation
    /// and by the existing step. No model or missing thermal state is inferred.
    pub fn from_parts(thermal: DepositionSimulation<'m>, circuit: MaterialCircuit,
        chambers: BTreeMap<u64, GasPiston>) -> Self {
        Self { thermal, circuit, chambers }
    }
    pub fn time_s(&self) -> f64 { self.thermal.time_s() }
    pub const fn thermal(&self) -> &DepositionSimulation<'m> { &self.thermal }
    pub const fn circuit(&self) -> &MaterialCircuit { &self.circuit }
    pub fn chambers(&self) -> &BTreeMap<u64, GasPiston> { &self.chambers }
    /// Explicitly dismantle the co-owned state, without advancing any subsystem.
    pub fn into_parts(self) -> (DepositionSimulation<'m>, MaterialCircuit, BTreeMap<u64, GasPiston>) {
        (self.thermal, self.circuit, self.chambers)
    }

    fn driven_chambers<'a>(
        &self, chamber_states: &'a BTreeMap<u64, GasPiston>, cx: &Cx<'_>,
        target_s: f64, drives: &[ChamberDrive], controls: NetworkSolve,
    ) -> Result<Vec<Chamber<'a>>, StateAdmissionError> {
        use StateAdmissionError as Error;
        cx.checkpoint().map_err(|_| Error::Cancelled)?;
        let dt = target_s - self.time_s();
        if !dt.is_finite() || dt <= 0.0 || self.chambers.len() > controls.maximum_chambers
            || drives.len() != self.chambers.len()
            || !chamber_states.keys().eq(self.chambers.keys()) {
            return Err(Error::Invalid("one forward assembly clock and one drive per live chamber required"));
        }
        let mut ordered = BTreeMap::new();
        for drive in drives {
            cx.checkpoint().map_err(|_| Error::Cancelled)?;
            if !self.chambers.contains_key(&drive.chamber)
                || ordered.insert(drive.chamber, drive).is_some() {
                return Err(Error::Invalid("duplicate drive or unknown live chamber"));
            }
        }
        Ok(ordered.iter().map(|(&id, drive)| Chamber {
            id, piston: &chamber_states[&id],
            input: PistonInput { duration_s: dt, heat_j: drive.heat_j,
                exterior_pressure_pa: drive.exterior_pressure_pa,
                external_force_n: drive.external_force_n },
            piston_solve: drive.piston_solve, exchange: drive.exchange,
        }).collect())
    }

    /// Preview one interval with a caller-selected gas heat law while retaining
    /// this live owner's accepted circuit, thermal, gas, momentum and port
    /// histories. The callback sees each reacted gas energy midpoint and exact
    /// trial cell-material sample; it must be pure and must not advance state.
    ///
    /// This is the persistent-state counterpart of
    /// `PreparedMaterialExchangeNetwork::reactive_with_gas_law_candidate`.
    /// It adds no constitutive law, solver or time discretization.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn step_with_gas_law_candidate<I, E>(
        &self, prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>, target_s: f64, ports: &[f64], drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>], surfaces: &[GasSurface],
        controls: NetworkSolve, guess_j: &[f64],
        resolve: impl Fn(Contact, &DepositionSimulation<'m>, &MaterialExchangeSample,
            &ResolvedThermochemicalGasState) -> Result<ConductanceSample<I>, E>,
    ) -> Result<(Self, MaterialExchangeAdvance<I>), ReactiveExchangeError<E>> {
        self.step_with_gas_law_from_chambers(&self.chambers, prepared, cx,
            target_s, ports, drives, reactions, surfaces, controls, guess_j, resolve)
    }

    /// Publish one caller-selected gas heat-law interval atomically. No
    /// fallible work follows the final checkpoint.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn advance_with_gas_law<I, E>(
        &mut self, prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>, target_s: f64, ports: &[f64], drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>], surfaces: &[GasSurface],
        controls: NetworkSolve, guess_j: &[f64],
        resolve: impl Fn(Contact, &DepositionSimulation<'m>, &MaterialExchangeSample,
            &ResolvedThermochemicalGasState) -> Result<ConductanceSample<I>, E>,
    ) -> Result<MaterialExchangeAdvance<I>, ReactiveExchangeError<E>> {
        let (next, report) = self.step_with_gas_law_candidate(prepared, cx,
            target_s, ports, drives, reactions, surfaces, controls, guess_j, resolve)?;
        cx.checkpoint().map_err(|_| ReactiveExchangeError::Cancelled)?;
        *self = next;
        Ok(report)
    }

    /// Preview one interval without modifying any accepted state. An empty
    /// reaction list uses the original nonreacting path; an empty chamber map
    /// gives electrical/thermal evolution. Every chamber still needs its own
    /// explicit load entry, even at equilibrium or with no gas-surface patches.
    /// All fixed-geometry, homogeneous-gas and first-order splitting limits of
    /// the original reactive/material exchange candidate remain unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate(&self, prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>, target_s: f64, ports: &[f64], drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>], gas_patches: &[MaterialContact],
        controls: NetworkSolve, guess_j: &[f64],
    ) -> Result<(Self, MaterialExchangeAdvance), ReactiveExchangeError> {
        self.step_from_chambers(&self.chambers, prepared, cx, target_s, ports, drives,
            reactions, gas_patches, controls, guess_j)
    }

    // Generic live-state assembly for any admitted gas-side heat law.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub(super) fn step_with_gas_law_from_chambers<I, E>(
        &self, chamber_states: &BTreeMap<u64, GasPiston>,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>,
        target_s: f64, ports: &[f64], drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>], surfaces: &[GasSurface],
        controls: NetworkSolve, guess_j: &[f64],
        resolve: impl Fn(Contact, &DepositionSimulation<'m>, &MaterialExchangeSample,
            &ResolvedThermochemicalGasState) -> Result<ConductanceSample<I>, E>,
    ) -> Result<(Self, MaterialExchangeAdvance<I>), ReactiveExchangeError<E>> {
        use ReactiveExchangeError as Error;
        let chambers = self.driven_chambers(chamber_states, cx, target_s, drives, controls)
            .map_err(StateAdmissionError::staged)?;
        let proposal = prepared.reactive_with_gas_law_candidate(cx, &self.thermal,
            target_s, &self.circuit, ports, &chambers, reactions, surfaces,
            controls, guess_j, resolve)?;
        let (thermal, circuit, report, step) = proposal.into_parts();
        let mut next_chambers = BTreeMap::new();
        let mut accounts = Vec::with_capacity(report.source.chambers.len());
        for chamber in report.source.chambers {
            cx.checkpoint().map_err(|_| Error::Cancelled)?;
            let id = chamber.id;
            if !self.chambers.contains_key(&id)
                || next_chambers.insert(id, chamber.piston).is_some() {
                return Err(Error::Invalid("candidate chamber identities differ from live assembly"));
            }
            accounts.push(ExchangeChamberStep {
                id, surfaces: chamber.surfaces, boundary: chamber.boundary,
            });
        }
        if next_chambers.len() != self.chambers.len() {
            return Err(Error::Invalid("candidate omitted a live chamber"));
        }
        let source = ExchangeAdvanceSources { material: report.source.material,
            chambers: accounts, reactions: report.source.reactions,
            transferred_heat_j: report.source.transferred_heat_j };
        let report = CircuitEnthalpyReport { source, electrical: report.electrical,
            thermal: report.thermal, nonlinear: report.nonlinear,
            maximum_cell_residual_j: report.maximum_cell_residual_j,
            circuit_trials: report.circuit_trials,
            deposited_joule_j: report.deposited_joule_j,
            non_electrical_heat_j: report.non_electrical_heat_j,
            total_energy_residual_j: report.total_energy_residual_j };
        cx.checkpoint().map_err(|_| Error::Cancelled)?;
        Ok((Self { thermal, circuit, chambers: next_chambers },
            MaterialExchangeAdvance { report, step }))
    }

    // Reuse the same state assembly with staged gas inventories, without cloning
    // the accepted thermal timeline or circuit merely to change gas inputs.
    #[allow(clippy::too_many_arguments)]
    fn step_from_chambers(&self, chamber_states: &BTreeMap<u64, GasPiston>,
        prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>, target_s: f64,
        ports: &[f64], drives: &[ChamberDrive], reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact], controls: NetworkSolve, guess_j: &[f64],
    ) -> Result<(Self, MaterialExchangeAdvance), ReactiveExchangeError> {
        use ReactiveExchangeError as Error;
        let chambers = self.driven_chambers(chamber_states, cx, target_s, drives, controls)
            .map_err(StateAdmissionError::staged)?;
        let proposal = prepared.reactive_step_candidate(cx, &self.thermal, target_s,
            &self.circuit, ports, &chambers, reactions, gas_patches, controls, guess_j)?;
        let (thermal, circuit, report, step) = proposal.into_parts();
        let mut next_chambers = BTreeMap::new();
        let mut accounts = Vec::with_capacity(report.source.chambers.len());
        for chamber in report.source.chambers {
            cx.checkpoint().map_err(|_| Error::Cancelled)?;
            let id = chamber.id;
            if !self.chambers.contains_key(&id) || next_chambers.insert(id, chamber.piston).is_some() {
                return Err(Error::Invalid("candidate chamber identities differ from live assembly"));
            }
            accounts.push(ExchangeChamberStep { id, surfaces: chamber.surfaces, boundary: chamber.boundary });
        }
        if next_chambers.len() != self.chambers.len() {
            return Err(Error::Invalid("candidate omitted a live chamber"));
        }
        let source = ExchangeAdvanceSources { material: report.source.material,
            chambers: accounts, reactions: report.source.reactions,
            transferred_heat_j: report.source.transferred_heat_j };
        let report = CircuitEnthalpyReport { source, electrical: report.electrical, thermal: report.thermal,
            nonlinear: report.nonlinear, maximum_cell_residual_j: report.maximum_cell_residual_j,
            circuit_trials: report.circuit_trials, deposited_joule_j: report.deposited_joule_j,
            non_electrical_heat_j: report.non_electrical_heat_j,
            total_energy_residual_j: report.total_energy_residual_j };
        cx.checkpoint().map_err(|_| Error::Cancelled)?;
        Ok((Self { thermal, circuit, chambers: next_chambers }, MaterialExchangeAdvance { report, step }))
    }

    /// Advance ALL owners or none. No fallible work follows the final checkpoint
    /// and no electrical/chemical/wall heat is applied a second time. Failed
    /// iteration, invalid mappings, source refusal and cancellation preserve the
    /// original thermal clock, circuit charge/flux, every gas and every port tick.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(&mut self, prepared: &PreparedMaterialExchangeNetwork<'_, '_>,
        cx: &Cx<'_>, target_s: f64, ports: &[f64], drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>], gas_patches: &[MaterialContact],
        controls: NetworkSolve, guess_j: &[f64],
    ) -> Result<MaterialExchangeAdvance, ReactiveExchangeError> {
        let (next, report) = self.step_candidate(prepared, cx, target_s, ports, drives,
            reactions, gas_patches, controls, guess_j)?;
        cx.checkpoint().map_err(|_| ReactiveExchangeError::Cancelled)?;
        *self = next;
        Ok(report)
    }
}
