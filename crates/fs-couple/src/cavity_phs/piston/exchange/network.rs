//! Open gas exchange through the existing circuit/enthalpy/thermal-wall solve.
//!
//! One nozzle transfer precedes one coupled network interval. Nonlinear network
//! trials all see the same transferred gas, never another reservoir withdrawal.

/// Same-card electrical, caloric and gas-wall conduction ingress.
pub mod material;

use super::{
    ClosedGasParcel, GasMassTransferStep, GasNozzle, GasNozzleFlow, GasNozzleStepSolve,
    GasPiston, NozzleDirection, NozzleTransfer, PistonError, PistonInput, poll, positive,
};
use crate::cavity_phs::piston::thermal::{ConductanceSample, ThermalPistonError};
use crate::cavity_phs::piston::thermal::boundary::network::{
    NetworkGasCandidate, NetworkGasDrive, NetworkGasError, NetworkGasReport, NetworkGasSolve,
};
use crate::circuit_enthalpy::{CircuitEnthalpyError, CircuitEnthalpyInput, CircuitThermalCoefficients};
use crate::material_circuit::MaterialCircuit;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use fs_thermal::enthalpy::DepositionSimulation;

/// Original stage reports and a balance enclosing every finite energy store.
pub struct NozzleNetworkReport<S, I> {
    /// Initial-pressure-selected upstream side; absent for a closed/no-flow path.
    pub direction: Option<NozzleDirection>,
    /// Source-resolved throat, mass rate and choking state.
    pub flow: GasNozzleFlow,
    /// Species and upstream enthalpy transfer, in original donor order.
    pub transfer: GasMassTransferStep,
    /// Original coupled circuit, thermal network, interface and motion reports.
    /// These stage accounts start AFTER the fixed-volume nozzle transfer.
    pub network: NetworkGasReport<S, I>,
    /// Full chamber internal-energy change, including nozzle enthalpy [J].
    pub chamber_energy_change_j: f64,
    /// Full finite-reservoir internal-energy change [J].
    pub reservoir_energy_change_j: f64,
    /// Independently measured total gas mass discrepancy [kg].
    pub total_mass_residual_kg: f64,
    /// Electrical + cell + both gas + kinetic energy change, plus exterior
    /// pressure work, minus force work, electrical supply and external heat [J].
    /// Nozzle enthalpy, Joule heat and reciprocal wall heat cancel internally.
    pub total_energy_residual_j: f64,
}

/// Four physical owners returned together, with no accepted input modified.
/// A material-bound piston retains its original body through the `P` parameter.
pub struct NozzleNetworkCandidate<'m, S, I, P = GasPiston> {
    thermal: DepositionSimulation<'m>,
    circuit: MaterialCircuit,
    piston: P,
    reservoir: ClosedGasParcel,
    report: NozzleNetworkReport<S, I>,
}
impl<'m, S, I, P> NozzleNetworkCandidate<'m, S, I, P> {
    /// Complete thermal history with material enthalpy and event cursor.
    pub const fn thermal(&self) -> &DepositionSimulation<'m> { &self.thermal }
    /// Proposed original circuit charge and flux stores.
    pub const fn circuit(&self) -> &MaterialCircuit { &self.circuit }
    /// Gas, volume, momentum and any retained source-bound body.
    pub const fn piston(&self) -> &P { &self.piston }
    /// Finite reservoir with its tentative species and caloric endpoint.
    pub const fn reservoir(&self) -> &ClosedGasParcel { &self.reservoir }
    /// Original source records and independently admitted whole-system balance.
    pub const fn report(&self) -> &NozzleNetworkReport<S, I> { &self.report }
    /// Publish all four endpoints together. Every reported internal transfer
    /// has already been applied; none is an additional heat source.
    pub fn into_parts(self) -> (
        DepositionSimulation<'m>, MaterialCircuit, P, ClosedGasParcel, NozzleNetworkReport<S, I>,
    ) {
        (self.thermal, self.circuit, self.piston, self.reservoir, self.report)
    }
}
impl<'m, S, I> NozzleNetworkCandidate<'m, S, I> {
    pub(in crate::cavity_phs::piston) fn map_piston<P>(
        self, bind: impl FnOnce(GasPiston) -> P,
    ) -> NozzleNetworkCandidate<'m, S, I, P> {
        NozzleNetworkCandidate { thermal: self.thermal, circuit: self.circuit,
            piston: bind(self.piston), reservoir: self.reservoir, report: self.report }
    }
}

fn flow_error<E, F>(error: PistonError) -> CircuitEnthalpyError<NetworkGasError<E, F>> {
    if error == PistonError::Cancelled {
        CircuitEnthalpyError::Cancelled
    } else {
        CircuitEnthalpyError::Source(NetworkGasError::Boundary(ThermalPistonError::Piston(error)))
    }
}

impl GasPiston {
    // Both public ingresses use this same transfer and full-endpoint admission.
    // Their existing network solvers continue to own material/source evaluation.
    #[allow(clippy::too_many_arguments)]
    fn nozzle_network_transfer<E, F>(
        &self, reservoir: &ClosedGasParcel, nozzle: GasNozzle,
        thermal: &DepositionSimulation<'_>, target_s: f64, drive: NetworkGasDrive,
        flow: GasNozzleStepSolve, controls: NetworkGasSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NozzleTransfer, CircuitEnthalpyError<NetworkGasError<E, F>>> {
        poll(cancelled).map_err(flow_error)?;
        if !positive(drive.input.duration_s) || !target_s.is_finite()
            || drive.input.duration_s != target_s - thermal.time_s()
            || !positive(controls.total_energy_tolerance_j)
        {
            return Err(CircuitEnthalpyError::Invalid(
                "nozzle and finite network need one positive physical interval and energy tolerance",
            ));
        }
        self.nozzle_transfer(reservoir, nozzle, drive.input.duration_s, flow, cancelled)
            .map_err(flow_error)
    }

    #[allow(clippy::too_many_arguments)]
    fn finish_nozzle_network<'m, S, I, E, F>(
        &self, reservoir: &ClosedGasParcel, mut next: NozzleTransfer,
        network: NetworkGasCandidate<'m, S, I>, drive: NetworkGasDrive,
        flow: GasNozzleStepSolve, controls: NetworkGasSolve,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<NozzleNetworkCandidate<'m, S, I>, CircuitEnthalpyError<NetworkGasError<E, F>>> {
        let (thermal, circuit, piston, stage) = network.into_parts();
        next.piston = piston;
        // The chamber's heat port contains external heat PLUS reciprocal wall
        // heat. This is not the heat supply of the enclosing finite system.
        let gas_input = PistonInput { heat_j: stage.boundary.piston.gas.supplied_heat_j,
            ..drive.input };
        let balance = self.nozzle_balance(reservoir, &next, gas_input, flow,
            controls.total_energy_tolerance_j).map_err(flow_error)?;
        let electrical = &stage.network.electrical;
        let total = ((electrical.stored_after_j - electrical.stored_before_j)
            + stage.network.thermal.stored_energy_change_j)
            + (balance.chamber_energy_j + balance.reservoir_energy_j)
            + balance.kinetic_energy_j + balance.exterior_work_j - balance.force_work_j
            - electrical.native.supplied - stage.coupled.external_network_heat_j
            - drive.input.heat_j;
        if !total.is_finite() || total.abs() > controls.total_energy_tolerance_j {
            return Err(CircuitEnthalpyError::Invalid(
                "nozzle, circuit, finite thermal cells, both gases and motion do not close total energy",
            ));
        }
        self.record_nozzle_account(&mut next, gas_input,
            stage.boundary.piston.midpoint.state().temperature, &balance,
            controls.total_energy_tolerance_j).map_err(flow_error)?;
        let report = NozzleNetworkReport { direction: next.direction, flow: next.flow,
            transfer: next.transfer, network: stage,
            chamber_energy_change_j: balance.chamber_energy_j,
            reservoir_energy_change_j: balance.reservoir_energy_j,
            total_mass_residual_kg: balance.mass_residual_kg, total_energy_residual_j: total };
        poll(cancelled).map_err(flow_error)?;
        Ok(NozzleNetworkCandidate { thermal, circuit, piston: next.piston,
            reservoir: next.reservoir, report })
    }

    /// Exchange gas species and enthalpy, then solve the existing finite
    /// electrical/thermal network, reciprocal wall heat and inertial motion.
    ///
    /// A single nozzle transfer uses accepted initial pressures. Every network
    /// trial then starts from the SAME transferred composition and accepted
    /// circuit/thermal history. Source callbacks see the actual gas energy
    /// midpoint and thermal endpoint, so changed inlet species affect pressure,
    /// conductivity and acoustic properties without a second gas model.
    ///
    /// `drive.input.heat_j` and network sources are external supplies only.
    /// Nozzle enthalpy, circuit Joule heat and reciprocal gas-wall heat are
    /// already exchanged by their original owners. All four endpoints and
    /// source records are returned together after independent mass/energy checks.
    /// Failure or cancellation publishes none of them. For nonzero flow the
    /// chamber's two-port audit explicitly lacks a stream port; use the complete
    /// paired balance in this report, rather than interpreting it as closed gas.
    ///
    /// This is a first-order transfer-then-network split over ONE duration, with
    /// a uniform ideal-gas reservoir, fixed thermal geometry and rigid piston.
    /// Existing donor-fraction, pressure-crossing, source support, event and work
    /// limits remain. Refine time; this does not infer combustion, jet thrust,
    /// variable circuit topology, spatial gas flow or a simultaneous nozzle solve.
    /// Callbacks must be repeatable and must not debit accepted state.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    pub fn nozzle_network_candidate<'m, S, I, E, F>(
        &self, reservoir: &ClosedGasParcel, nozzle: GasNozzle,
        thermal: &DepositionSimulation<'m>, target_s: f64,
        electrical: CircuitEnthalpyInput<'_>, drive: NetworkGasDrive,
        flow: GasNozzleStepSolve, controls: NetworkGasSolve, guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'m>) -> Result<CircuitThermalCoefficients<S>, E>,
        resolve_port: impl Fn(&DepositionSimulation<'m>, &S, &ResolvedThermochemicalGasState)
            -> Result<ConductanceSample<I>, F>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<NozzleNetworkCandidate<'m, S, I>, CircuitEnthalpyError<NetworkGasError<E, F>>> {
        let next = self.nozzle_network_transfer(reservoir, nozzle, thermal, target_s,
            drive, flow, controls, &mut cancelled)?;
        let network = next.piston.network_candidate(thermal, target_s, electrical, drive,
            controls, guess_j, resolve, resolve_port, &mut cancelled)?;
        self.finish_nozzle_network(reservoir, next, network, drive, flow, controls, &mut cancelled)
    }
}

#[cfg(test)]
mod tests;
