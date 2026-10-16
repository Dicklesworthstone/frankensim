//! Finite circuit/enthalpy feedback driving BOTH gases on one shared mechanism.
//! No numerical, caloric, conduction, circuit or pressure law is duplicated.

/// Data-only paired contacts on the original material network.
pub mod material;

use super::*;
use crate::cavity_phs::piston::thermal::boundary::network::{NetworkGasError, NetworkGasSolve, NetworkGasStep};
use crate::circuit_enthalpy::{self, CircuitEnthalpyCandidate, CircuitEnthalpyError,
    CircuitEnthalpyInput, CircuitEnthalpyReport, CircuitThermalCoefficients};
use crate::material_circuit::MaterialCircuit;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::cell::{Cell, RefCell};
use std::collections::BTreeMap;

/// Two finite thermal contact destinations on the existing shared mechanism.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DoubleNetworkDrive {
    /// [Primary, opposing] active cells. The same cell is allowed for separate
    /// non-overlapping patches; its heat store is still counted only once.
    pub cells: [CellId; 2],
    /// Prescribed gas heat and common external mechanical load.
    pub input: DoublePistonInput,
    /// Original coupled pressure/momentum and material work bounds.
    pub piston: DoublePistonSolve,
    /// Both boundary heat brackets and the total paired-proposal work limit.
    pub exchange: DoubleHeatExchangeSolve,
}

/// Exact endpoint source evaluation retained by the original network solver.
pub struct DoubleNetworkSource<S, I> {
    /// Original electrical and thermal material results.
    pub material: S,
    /// Both gas endpoints and both boundary source evaluations.
    pub boundary: DoubleBoundaryCandidate<I>,
    /// Heat actually subtracted from cell source powers, reconstructed AFTER
    /// floating-point insertion and aggregated once per destination [J].
    pub debited_heat_j: f64,
    external_sources_w: BTreeMap<CellId, f64>,
    conductances_w_k: BTreeMap<u64, f64>,
}

/// Existing circuit, thermal cells, two gases and one momentum accepted together.
pub struct DoubleNetworkCandidate<'m, S, I> {
    network: CircuitEnthalpyCandidate<'m, DoubleNetworkSource<S, I>>,
    cells: [CellId; 2],
    step: NetworkGasStep,
}
/// Original per-owner reports and the independently evaluated combined balance.
pub struct DoubleNetworkReport<S, I> {
    /// Network account; its non-electrical heat includes both internal gas debits.
    pub network: CircuitEnthalpyReport<S>,
    /// Original two-port heat and shared pressure/momentum account.
    pub boundary: DoubleBoundaryStep,
    /// Both original interface evaluations, [primary, opposing].
    pub interfaces: [I; 2],
    /// Original finite heat destinations, in the same face order.
    pub cells: [CellId; 2],
    /// Actual total cell-source debit [J], not assumed from requested gas heat.
    pub debited_heat_j: f64,
    /// All-owner balance and cumulative work over all nonlinear probes.
    pub coupled: NetworkGasStep,
}
impl<'m, S, I> DoubleNetworkCandidate<'m, S, I> {
    /// Matching network endpoint with its exact paired boundary proposal.
    pub const fn network(&self) -> &CircuitEnthalpyCandidate<'m, DoubleNetworkSource<S, I>> { &self.network }
    /// Matching gas states and the common mechanical trajectory.
    pub fn boundary(&self) -> &DoubleBoundaryCandidate<I> { &self.network.report().source.boundary }
    /// Combined energy and all discarded/final boundary work.
    pub const fn step(&self) -> &NetworkGasStep { &self.step }
    /// Publish every endpoint together. Do not re-apply the returned Joule or gas heats.
    pub fn into_parts(self) -> (DepositionSimulation<'m>, MaterialCircuit, DoubleActingGasPiston, DoubleNetworkReport<S, I>) {
        let (thermal, circuit, report) = self.network.into_parts();
        let source = report.source;
        let (piston, interfaces, boundary) = source.boundary.into_parts();
        let network = CircuitEnthalpyReport { source: source.material, electrical: report.electrical,
            thermal: report.thermal, nonlinear: report.nonlinear,
            maximum_cell_residual_j: report.maximum_cell_residual_j, circuit_trials: report.circuit_trials,
            deposited_joule_j: report.deposited_joule_j, non_electrical_heat_j: report.non_electrical_heat_j,
            total_energy_residual_j: report.total_energy_residual_j };
        (thermal, circuit, piston, DoubleNetworkReport { network, boundary, interfaces, cells: self.cells,
            debited_heat_j: source.debited_heat_j, coupled: self.step })
    }
}

impl DoubleActingGasPiston {
    /// Close thermal, electrical, two-gas and common mechanical feedback.
    ///
    /// Every original network Newton trial supplies its two cell temperatures.
    /// The paired boundary solve starts from the SAME accepted gases/momentum;
    /// circuit proposals likewise start from the same accepted charges/fluxes.
    /// Both signed gas heats are grouped by cell BEFORE one source-map debit.
    /// Different cells retain their separate heats even when their sum is zero.
    /// No callback may advance or debit external accepted physical history.
    ///
    /// Final external supply is evaluated independently by the original thermal
    /// rate owner with the unmodified external-source map. It is NOT inferred
    /// by adding requested gas heat to an already combined network account.
    /// This makes lost/double-applied transfers visible in the all-owner gate.
    /// The original thermal local residuals and circuit work checks still apply.
    ///
    /// Fixed circuit/thermal topology, gas masses and rigid geometry. Thermal
    /// endpoint rates combined with gas/circuit midpoints are generally FIRST
    /// order. No combustion, valves, conductor-motion work, moving-grid thermal
    /// conduction or automatic acoustic reconfiguration is implied. Custom
    /// callbacks own contact coverage; the material adapter checks it explicitly.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn network_candidate<'m, S, I, E, F>(
        &self, thermal: &DepositionSimulation<'m>, target_s: f64,
        electrical: CircuitEnthalpyInput<'_>, drive: DoubleNetworkDrive, controls: NetworkGasSolve,
        guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'m>) -> Result<CircuitThermalCoefficients<S>, E>,
        resolve_ports: impl Fn(&DepositionSimulation<'m>, &S, &[ResolvedThermochemicalGasState; 2])
            -> Result<[ConductanceSample<I>; 2], F>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DoubleNetworkCandidate<'m, S, I>, CircuitEnthalpyError<NetworkGasError<E, F>>> {
        let cancel = RefCell::new(cancelled);
        let check = || if (*cancel.borrow_mut())() { Err(DepositionError::Cancelled) } else { Ok(()) };
        check()?;
        if drive.cells.iter().any(|id| !thermal.cells().contains_key(id))
            || drive.input.duration_s != target_s-thermal.time_s()
            || controls.maximum_boundary_solves == 0 || !positive(controls.total_energy_tolerance_j) {
            return Err(CircuitEnthalpyError::Invalid("active paired heat cells, one clock and bounded work required"));
        }
        let tolerance = drive.exchange.faces[0].energy_tolerance_j.min(drive.exchange.faces[1].energy_tolerance_j);
        let calls = Cell::new(0usize); let heat_trials = Cell::new(0usize); let gas_trials = Cell::new(0usize);
        let candidate = circuit_enthalpy::step_candidate(thermal, target_s, electrical, controls.network, guess_j, |state| {
            if calls.get() >= controls.maximum_boundary_solves { return Err(NetworkGasError::Budget); }
            let mut coefficients = resolve(state).map_err(NetworkGasError::Material)?;
            let temperatures = [state.temperature_k(drive.cells[0]).map_err(NetworkGasError::Thermal)?,
                state.temperature_k(drive.cells[1]).map_err(NetworkGasError::Thermal)?];
            calls.set(calls.get()+1);
            let boundary = self.thermal_boundaries_candidate(drive.input, drive.piston, drive.exchange, temperatures,
                |gas| resolve_ports(state, &coefficients.source, gas), || (*cancel.borrow_mut())())
                .map_err(NetworkGasError::Boundary)?;
            heat_trials.set(heat_trials.get().checked_add(boundary.step().joint_trials).ok_or(NetworkGasError::Budget)?);
            gas_trials.set(gas_trials.get().checked_add(boundary.step().material_trials).ok_or(NetworkGasError::Budget)?);
            let external_sources_w = coefficients.sources_w.clone();
            let mut debits = BTreeMap::new();
            for (cell, q) in drive.cells.into_iter().zip(boundary.step().exchanged_heat_j) {
                let total = debits.entry(cell).or_insert(0.0);
                *total += q;
                if !total.is_finite() { return Err(NetworkGasError::Transfer); }
            }
            let mut debited_heat_j = 0.0;
            for (cell, q) in debits {
                let watts = q/drive.input.duration_s;
                let source = coefficients.sources_w.entry(cell).or_insert(0.0);
                let previous = *source; let next = previous-watts;
                let actual = drive.input.duration_s*(previous-next);
                if !watts.is_finite() || !next.is_finite() || !actual.is_finite()
                    || (q != 0.0 && (watts == 0.0 || next == previous)) || (actual-q).abs() > tolerance {
                    return Err(NetworkGasError::Transfer);
                }
                *source = next; debited_heat_j += actual;
            }
            if !debited_heat_j.is_finite() { return Err(NetworkGasError::Transfer); }
            let source = DoubleNetworkSource { material: coefficients.source, boundary, debited_heat_j,
                external_sources_w, conductances_w_k: coefficients.conductances_w_k.clone() };
            Ok(CircuitThermalCoefficients { resistors: coefficients.resistors, sources_w: coefficients.sources_w,
                conductances_w_k: coefficients.conductances_w_k, source })
        }, || (*cancel.borrow_mut())());
        let candidate = match candidate {
            Err(CircuitEnthalpyError::Source(NetworkGasError::Boundary(ThermalPistonError::Cancelled))) =>
                return Err(CircuitEnthalpyError::Cancelled),
            other => other?,
        };
        let report = candidate.report(); let source = &report.source;
        let p = &source.boundary.step().piston;
        let external = drive.input.duration_s*candidate.thermal().external_step_power_w(
            &source.external_sources_w, &source.conductances_w_k, check)?;
        let total = (report.electrical.stored_after_j-report.electrical.stored_before_j)
            + report.thermal.stored_energy_change_j + p.gas[0].stored_energy_change_j
            + p.gas[1].stored_energy_change_j + p.kinetic_energy_change_j - p.external_force_work_j
            - report.electrical.native.supplied - external - drive.input.heat_j[0] - drive.input.heat_j[1];
        if !external.is_finite() || !total.is_finite() || total.abs() > controls.total_energy_tolerance_j {
            return Err(CircuitEnthalpyError::Invalid("circuit, finite cells, two gases and one motion do not close external work"));
        }
        check()?;
        Ok(DoubleNetworkCandidate { network: candidate, cells: drive.cells,
            step: NetworkGasStep { external_network_heat_j: external, total_energy_residual_j: total,
                boundary_solves: calls.get(), exchange_trials: heat_trials.get(), material_trials: gas_trials.get() } })
    }
}

#[cfg(test)]
mod tests;
