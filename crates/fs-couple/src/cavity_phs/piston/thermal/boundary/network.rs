//! Reciprocal gas/momentum coupling to the EXISTING circuit/enthalpy network.
//! Network material states are implicit endpoints; gas force/energy use their
//! original midpoint. The combined scheme is generally first order in time.

/// Same-card conductive boundary with no handwritten material callback.
pub mod material;

use super::*;
use crate::circuit_enthalpy::{self, CircuitEnthalpyCandidate, CircuitEnthalpyError,
    CircuitEnthalpyInput, CircuitEnthalpyReport, CircuitEnthalpySolve, CircuitThermalCoefficients};
use crate::material_circuit::MaterialCircuit;
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::cell::{Cell, RefCell};

/// The one finite thermal cell exchanging heat with this gas cavity.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NetworkGasDrive {
    /// Existing active thermal cell; no duplicate finite body is created.
    pub cell: CellId,
    /// Existing gas heat and mechanical loads, held during this interval.
    pub input: PistonInput,
    /// Existing pressure/momentum work limits.
    pub piston: PistonSolve,
    /// Boundary exchange limits in signed joules INTO the gas.
    pub exchange: HeatExchangeSolve,
}
/// Global work and energy limits in addition to each original owner's controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NetworkGasSolve {
    /// Original Newton, circuit, thermal and local energy admission.
    pub network: CircuitEnthalpySolve,
    /// All boundary solves across numerical probes and final admission.
    /// Each has its own bounded heat, momentum and material iterations.
    pub maximum_boundary_solves: usize,
    /// Absolute all-owner energy discrepancy [J], independent of Newton stopping.
    pub total_energy_tolerance_j: f64,
}
/// Original source failures are not mistaken for a converged thermal residual.
#[derive(Debug)]
pub enum NetworkGasError<E, F> {
    /// Existing material/network callback refused a trial endpoint.
    Material(E),
    /// Original gas, mechanics or thermal boundary callback refused a trial.
    Boundary(ThermalPistonError<F>),
    /// Global boundary-solve budget exhausted, including final checks.
    Budget,
    /// An exchanged heat vanished during transfer to the thermal source map.
    Transfer,
    /// Original thermal-temperature refusal.
    Thermal(DepositionError),
}
impl<E: core::fmt::Debug,F: core::fmt::Debug> core::fmt::Display for NetworkGasError<E,F> {
    fn fmt(&self,f:&mut core::fmt::Formatter<'_>)->core::fmt::Result {write!(f,"network/gas source: {self:?}")}
}
impl<E: std::error::Error+'static,F: std::error::Error+'static> std::error::Error for NetworkGasError<E,F> {}

/// All-owner energy balance; internal electric and gas heat are excluded from supply.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NetworkGasStep {
    /// Source/reservoir heat into the NETWORK excluding circuit and gas exchange [J].
    pub external_network_heat_j: f64,
    /// Actual circuit + cell + gas + kinetic change plus exterior pressure work,
    /// minus force work, external electrical supply and prescribed external heat.
    pub total_energy_residual_j: f64,
    /// All boundary solves, including discarded Newton/Jacobian proposals.
    pub boundary_solves: usize,
    /// All heat-feedback trials within those boundary solves.
    pub exchange_trials: usize,
    /// All nested gas material proposals within those boundary solves.
    pub material_trials: usize,
}
/// Three existing physical owners accepted together, without copying their laws.
pub struct NetworkGasCandidate<'m,S,I> {
    network: CircuitEnthalpyCandidate<'m,(S,PistonBoundaryCandidate<I>)>,
    step: NetworkGasStep,
}
/// Reports from all original owners, with material and boundary sources retained.
pub struct NetworkGasReport<S,I> {
    /// Circuit/thermal report. Its non-electrical heat includes the gas debit.
    pub network: CircuitEnthalpyReport<S>,
    /// Original gas/momentum and boundary-heat account.
    pub boundary: PistonBoundaryStep,
    /// Original interface material evaluation at the accepted states.
    pub interface: I,
    /// All-owner energy and accumulated nonlinear work.
    pub coupled: NetworkGasStep,
}
impl<'m,S,I> NetworkGasCandidate<'m,S,I> {
    /// The unchanged-format circuit and thermal proposal; its source contains
    /// the matching material result and gas proposal, not a second gas history.
    pub const fn network(&self)->&CircuitEnthalpyCandidate<'m,(S,PistonBoundaryCandidate<I>)> {&self.network}
    /// Matching gas, momentum and original boundary source.
    pub fn boundary(&self)->&PistonBoundaryCandidate<I> {&self.network.report().source.1}
    /// Independent combined work and cumulative solver accounting.
    pub const fn step(&self)->&NetworkGasStep {&self.step}
    /// Publish all endpoints together. Both heat transfers are already applied.
    pub fn into_parts(self)->(DepositionSimulation<'m>,MaterialCircuit,GasPiston,NetworkGasReport<S,I>) {
        let (thermal,circuit,report)=self.network.into_parts();
        let (material,boundary)=report.source;
        let (piston,interface,boundary)=boundary.into_parts();
        let network=CircuitEnthalpyReport {source:material,electrical:report.electrical,thermal:report.thermal,
            nonlinear:report.nonlinear,maximum_cell_residual_j:report.maximum_cell_residual_j,
            circuit_trials:report.circuit_trials,deposited_joule_j:report.deposited_joule_j,
            non_electrical_heat_j:report.non_electrical_heat_j,total_energy_residual_j:report.total_energy_residual_j};
        (thermal,circuit,piston,NetworkGasReport {network,boundary,interface,coupled:self.step})
    }
}

impl GasPiston {
    /// Close thermal, electrical and gas/mechanical feedback over ONE interval.
    /// The network's own Newton residual sees its actual endpoint temperatures
    /// and material coefficients. For EACH such trial the existing boundary
    /// solver starts from this SAME gas state, solves motion and heat together,
    /// and subtracts exactly that heat from `drive.cell`. The electrical owner
    /// starts from its SAME accepted charges/fluxes. No iteration drains history.
    ///
    /// `resolve_port` sees the network's original material payload as well as
    /// the actual gas midpoint. The caller can reuse that endpoint evaluation.
    /// Extra network sources, reservoirs and prescribed gas heat are retained.
    /// Callbacks must be repeatable, with no external physical side effects.
    ///
    /// All owners retain their domain, work, event and cancellation checks.
    /// Fixed thermal geometry/C/L/circuit topology; one uniform closed gas and
    /// its existing translating wall. Thermal endpoint rates + gas midpoint
    /// give first-order coupled accuracy in general. This is not combustion,
    /// moving-mesh conduction, shared-piston multicylinder dynamics or an
    /// automatic acoustic solver reconfiguration. Surface coverage for a custom
    /// boundary is the caller's responsibility; the material adapter admits a
    /// non-overlapping exposed patch and refuses concurrent ambient convection.
    #[allow(clippy::too_many_arguments,clippy::type_complexity)]
    pub fn network_candidate<'m,S,I,E,F>(
        &self, thermal:&DepositionSimulation<'m>, target_s:f64,
        electrical:CircuitEnthalpyInput<'_>, drive:NetworkGasDrive, controls:NetworkGasSolve,
        guess_j:&[f64],
        resolve:impl Fn(&DepositionSimulation<'m>)->Result<CircuitThermalCoefficients<S>,E>,
        resolve_port:impl Fn(&DepositionSimulation<'m>,&S,&ResolvedThermochemicalGasState)
            ->Result<ConductanceSample<I>,F>,
        cancelled:impl FnMut()->bool,
    )->Result<NetworkGasCandidate<'m,S,I>,CircuitEnthalpyError<NetworkGasError<E,F>>> {
        let cancel=RefCell::new(cancelled);
        let check=|| if (*cancel.borrow_mut())() {Err(DepositionError::Cancelled)} else {Ok(())};
        check()?;
        if !thermal.cells().contains_key(&drive.cell) || drive.input.duration_s!=target_s-thermal.time_s()
            || controls.maximum_boundary_solves==0 || !positive(controls.total_energy_tolerance_j) {
            return Err(CircuitEnthalpyError::Invalid("active thermal gas boundary, one clock and bounded work required"));
        }
        let calls=Cell::new(0usize);let heat_trials=Cell::new(0usize);let gas_trials=Cell::new(0usize);
        let candidate=circuit_enthalpy::step_candidate(thermal,target_s,electrical,controls.network,guess_j,|state| {
            // Each complete gas-boundary solve is charged BEFORE doing its work.
            if calls.get()>=controls.maximum_boundary_solves {return Err(NetworkGasError::Budget);}
            let mut coefficients=resolve(state).map_err(NetworkGasError::Material)?;
            // The thermal owner has already admitted each endpoint temperature.
            let t=state.temperature_k(drive.cell).map_err(NetworkGasError::Thermal)?;
            calls.set(calls.get()+1);
            let boundary=self.thermal_boundary_candidate(drive.input,drive.piston,drive.exchange,t,
                |gas|resolve_port(state,&coefficients.source,gas),||(*cancel.borrow_mut())())
                .map_err(NetworkGasError::Boundary)?;
            heat_trials.set(heat_trials.get()+boundary.step().exchange_trials);
            gas_trials.set(gas_trials.get()+boundary.step().material_trials);
            let q=boundary.step().exchanged_heat_j;
            let watts=q/drive.input.duration_s;
            let source=coefficients.sources_w.entry(drive.cell).or_insert(0.0);
            let previous=*source;
            let next=previous-watts;
            if !watts.is_finite() || !next.is_finite() || (q!=0.0 && (watts==0.0 || next==previous))
                || (drive.input.duration_s*watts-q).abs()>drive.exchange.energy_tolerance_j {
                return Err(NetworkGasError::Transfer);
            }
            *source=next;
            Ok(CircuitThermalCoefficients {resistors:coefficients.resistors,
                sources_w:coefficients.sources_w,conductances_w_k:coefficients.conductances_w_k,
                source:(coefficients.source,boundary)})
        },||(*cancel.borrow_mut())());
        let candidate=match candidate {
            Err(CircuitEnthalpyError::Source(NetworkGasError::Boundary(ThermalPistonError::Cancelled))) =>
                return Err(CircuitEnthalpyError::Cancelled),
            other=>other?,
        };
        let report=candidate.report();
        let boundary=&report.source.1;
        let p=&boundary.step().piston;
        let external=report.non_electrical_heat_j+boundary.step().exchanged_heat_j;
        let total=(report.electrical.stored_after_j-report.electrical.stored_before_j)
            +report.thermal.stored_energy_change_j+p.gas.stored_energy_change_j+p.kinetic_energy_change_j
            +p.exterior_pressure_work_j-p.external_force_work_j-report.electrical.native.supplied
            -external-drive.input.heat_j;
        if !external.is_finite() || !total.is_finite() || total.abs()>controls.total_energy_tolerance_j {
            return Err(CircuitEnthalpyError::Invalid("circuit, network, gas and motion do not close external work"));
        }
        check()?;
        Ok(NetworkGasCandidate {network:candidate,step:NetworkGasStep {external_network_heat_j:external,
            total_energy_residual_j:total,boundary_solves:calls.get(),
            exchange_trials:heat_trials.get(),material_trials:gas_trials.get()}})
    }
}
