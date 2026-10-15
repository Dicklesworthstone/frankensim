//! Separate thermal-cell heat debits for surfaces sharing ONE gas cavity.
//! Reuses circuit_enthalpy and the common multi-surface gas boundary solve.
//! Electrical storage, cell enthalpy and piston history never advance in trials.

/// Data-only multiple-surface material binding with combined area admission.
pub mod material;

use super::*;
use crate::cavity_phs::piston::thermal::boundary::ports::{
    ThermalPortsCandidate, BoundaryHeatStep, ThermalPort, ThermalPortSample,
};
use std::collections::BTreeMap;

/// Hard limit on surfaces retained by one coupled network/gas interval.
pub const MAX_NETWORK_SURFACES: usize = 256;

/// Stable physical surface identity and its original finite thermal owner.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct NetworkThermalPort {
    /// Unique surface ID within this cavity, independent of cell or branch order.
    pub id: u64,
    /// Existing active cell. Several distinct non-overlapping patches may name it.
    pub cell: CellId,
}

/// One gas's sources and motion, shared across all its thermal surfaces.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NetworkSurfacesDrive<'a> {
    /// Complete bounded gas-surface binding. Empty is not an implicit insulation.
    pub ports: &'a [NetworkThermalPort],
    /// Prescribed gas heat and mechanical loads for this interval.
    pub input: PistonInput,
    /// Original pressure/momentum solve and its work bounds.
    pub piston: PistonSolve,
    /// NET heat bracket [J]; each individual surface may transfer larger heat.
    pub exchange: HeatExchangeSolve,
}

/// All proposed physical histories and the exact signed per-surface transfer.
pub struct NetworkSurfacesCandidate<'m, S, I> {
    network: CircuitEnthalpyCandidate<'m, (S, ThermalPortsCandidate<I>)>,
    destinations: Vec<NetworkThermalPort>,
    step: NetworkGasStep,
}

/// Original electrical, thermal and gas/momentum accounts for the same interval.
pub struct NetworkSurfacesReport<S, I> {
    /// Network account; non-electrical heat includes the surface debits.
    pub network: CircuitEnthalpyReport<S>,
    /// Original shared gas midpoint, mechanical work and net heat residual.
    pub boundary: BoundaryHeatStep,
    /// Actual signed surface heats and material records, in increasing ID order.
    pub ports: Vec<ThermalPortSample<I>>,
    /// Original heat destinations in the same ID order as `ports`.
    pub destinations: Vec<NetworkThermalPort>,
    /// All-owner energy and total nested boundary work.
    pub coupled: NetworkGasStep,
}
impl<'m, S, I> NetworkSurfacesCandidate<'m, S, I> {
    /// Matching circuit/thermal state with its source and gas candidate.
    pub const fn network(&self) -> &CircuitEnthalpyCandidate<'m, (S, ThermalPortsCandidate<I>)> { &self.network }
    /// All surfaces share this one gas and one mechanical trajectory.
    pub fn boundary(&self) -> &ThermalPortsCandidate<I> { &self.network.report().source.1 }
    /// Original surface-to-cell routing in increasing surface ID order.
    pub fn destinations(&self) -> &[NetworkThermalPort] { &self.destinations }
    /// Independently checked combined energy and nonlinear work.
    pub const fn step(&self) -> &NetworkGasStep { &self.step }
    /// Publish all histories together. Do not re-apply any returned surface heat.
    pub fn into_parts(self) -> (DepositionSimulation<'m>, MaterialCircuit, GasPiston, NetworkSurfacesReport<S, I>) {
        let (thermal, circuit, report) = self.network.into_parts();
        let (source, boundary) = report.source;
        let (gas, ports, boundary) = boundary.into_parts();
        let network = CircuitEnthalpyReport { source, electrical: report.electrical, thermal: report.thermal,
            nonlinear: report.nonlinear, maximum_cell_residual_j: report.maximum_cell_residual_j,
            circuit_trials: report.circuit_trials, deposited_joule_j: report.deposited_joule_j,
            non_electrical_heat_j: report.non_electrical_heat_j, total_energy_residual_j: report.total_energy_residual_j };
        (thermal, circuit, gas, NetworkSurfacesReport { network, boundary, ports,
            destinations: self.destinations, coupled: self.step })
    }
}

impl GasPiston {
    /// Couple multiple finite wall surfaces to ONE uniform gas and moving wall.
    /// The existing network Newton residual supplies every wall's endpoint T.
    /// The common gas boundary solves net heat and momentum together while
    /// retaining Q_i = dt G_i(T_i - T_gas_mid) for EVERY surface, even when the
    /// sum is zero. Each exact evaluated Q_i is debited only from its own cell.
    ///
    /// Every boundary/circuit proposal starts from the accepted histories. A
    /// fresh final network evaluation keeps the matching gas and all material
    /// records. Local thermal gates cannot be bypassed by cancelling cell errors.
    /// The all-owner account includes the difference between applied net gas
    /// heat and evaluated per-surface heats; no artificial split of that residual.
    ///
    /// The custom callback owns surface coverage and model validity. The material
    /// adapter checks summed patch coverage per cell and forbids duplicate ambient
    /// convection. All other physical limits of `network_candidate` remain:
    /// fixed thermal/circuit topology, uniform closed gas, one inertial wall,
    /// generally first-order coupled accuracy, no mass transfer or combustion.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn network_surfaces_candidate<'m, S, I, E, F>(
        &self, thermal: &DepositionSimulation<'m>, target_s: f64,
        electrical: CircuitEnthalpyInput<'_>, drive: NetworkSurfacesDrive<'_>, controls: NetworkGasSolve,
        guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'m>) -> Result<CircuitThermalCoefficients<S>, E>,
        resolve_port: impl Fn(&NetworkThermalPort, &DepositionSimulation<'m>, &S, &ResolvedThermochemicalGasState)
            -> Result<ConductanceSample<I>, F>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<NetworkSurfacesCandidate<'m, S, I>, CircuitEnthalpyError<NetworkGasError<E, F>>> {
        let cancel = RefCell::new(cancelled);
        let check = || if (*cancel.borrow_mut())() { Err(DepositionError::Cancelled) } else { Ok(()) };
        check()?;
        if drive.ports.is_empty() || drive.ports.len() > MAX_NETWORK_SURFACES
            || drive.input.duration_s != target_s-thermal.time_s()
            || controls.maximum_boundary_solves == 0 || !positive(controls.total_energy_tolerance_j) {
            return Err(CircuitEnthalpyError::Invalid("bounded nonempty gas surfaces, one clock and finite work limits required"));
        }
        let mut destinations = drive.ports.to_vec();
        destinations.sort_by_key(|p| p.id);
        for (i, port) in destinations.iter().enumerate() {
            check()?;
            if !thermal.cells().contains_key(&port.cell) || (i > 0 && destinations[i-1].id == port.id) {
                return Err(CircuitEnthalpyError::Invalid("unique surface identities and active receiving cells required"));
            }
        }
        let calls = Cell::new(0usize); let heat_trials = Cell::new(0usize); let gas_trials = Cell::new(0usize);
        let candidate = circuit_enthalpy::step_candidate(thermal, target_s, electrical, controls.network, guess_j, |state| {
            if calls.get() >= controls.maximum_boundary_solves { return Err(NetworkGasError::Budget); }
            let mut coefficients = resolve(state).map_err(NetworkGasError::Material)?;
            let boundaries: Vec<_> = destinations.iter().map(|p| {
                state.temperature_k(p.cell).map(|temperature_k| ThermalPort { id: p.id, temperature_k })
                    .map_err(NetworkGasError::Thermal)
            }).collect::<Result<_, _>>()?;
            calls.set(calls.get()+1);
            let boundary = self.thermal_ports_candidate(drive.input, drive.piston, drive.exchange, &boundaries, MAX_NETWORK_SURFACES,
                |port, gas| {
                    let index = destinations.binary_search_by_key(&port.id, |p| p.id).expect("admitted surface ID");
                    resolve_port(&destinations[index], state, &coefficients.source, gas)
                }, || (*cancel.borrow_mut())()).map_err(NetworkGasError::Boundary)?;
            heat_trials.set(heat_trials.get()+boundary.step().exchange_trials);
            gas_trials.set(gas_trials.get()+boundary.step().material_trials);
            let mut debits = BTreeMap::new();
            // No zero-net early exit: balanced hot/cold walls must both evolve.
            for (destination, port) in destinations.iter().zip(boundary.source()) {
                let heat = debits.entry(destination.cell).or_insert(0.0);
                *heat += port.heat_j;
                if !heat.is_finite() { return Err(NetworkGasError::Transfer); }
            }
            let mut transferred = 0.0;
            for (cell, q) in debits {
                let watts = q/drive.input.duration_s;
                let source = coefficients.sources_w.entry(cell).or_insert(0.0);
                let previous = *source; let next = previous-watts;
                let actual = drive.input.duration_s*watts;
                if !watts.is_finite() || !next.is_finite() || !actual.is_finite()
                    || (q != 0.0 && (watts == 0.0 || next == previous))
                    || (actual-q).abs() > drive.exchange.energy_tolerance_j {
                    return Err(NetworkGasError::Transfer);
                }
                *source = next; transferred += actual;
            }
            if !transferred.is_finite()
                || (transferred-boundary.step().evaluated_heat_j).abs() > drive.exchange.energy_tolerance_j {
                return Err(NetworkGasError::Transfer);
            }
            Ok(CircuitThermalCoefficients { resistors: coefficients.resistors, sources_w: coefficients.sources_w,
                conductances_w_k: coefficients.conductances_w_k, source: (coefficients.source, boundary) })
        }, || (*cancel.borrow_mut())());
        let candidate = match candidate {
            Err(CircuitEnthalpyError::Source(NetworkGasError::Boundary(ThermalPistonError::Cancelled))) =>
                return Err(CircuitEnthalpyError::Cancelled),
            other => other?,
        };
        let report = candidate.report(); let boundary = &report.source.1;
        let p = &boundary.step().piston;
        let external = report.non_electrical_heat_j+boundary.step().evaluated_heat_j;
        let total = (report.electrical.stored_after_j-report.electrical.stored_before_j)
            +report.thermal.stored_energy_change_j+p.gas.stored_energy_change_j+p.kinetic_energy_change_j
            +p.exterior_pressure_work_j-p.external_force_work_j-report.electrical.native.supplied
            -external-drive.input.heat_j;
        if !external.is_finite() || !total.is_finite() || total.abs() > controls.total_energy_tolerance_j {
            return Err(CircuitEnthalpyError::Invalid("circuit, all walls and shared gas/motion do not close external work"));
        }
        check()?;
        Ok(NetworkSurfacesCandidate { network: candidate, destinations,
            step: NetworkGasStep { external_network_heat_j: external, total_energy_residual_j: total,
                boundary_solves: calls.get(), exchange_trials: heat_trials.get(), material_trials: gas_trials.get() } })
    }
}

#[cfg(test)]
mod tests;
