//! Complete material heat transfer on the original current and voltage motors.
//!
//! The native rotor, winding response, source limits and heat deposition are
//! shared with the bare material-network entry points. Only the transport
//! provider changes; no second thermal or mechanical step is introduced.

use super::*;
use super::voltage::{MaterialRotorVoltageDrive, RotorVoltageEnthalpyCandidate};
use crate::circuit_enthalpy::contacts::heat::{
    HeatAssembly, HeatNetworkError, closed_balance,
};
use fs_exec::Cx;
use fs_mbd::rotordynamics::transient::drive::voltage::{
    RotorVoltageControl, RotorVoltageMotor, RotorWindingState,
};

/// Preserve native motor failures and the complete assembly's typed refusal.
#[derive(Debug)]
pub enum MotorHeatError<E> {
    Motor(MotorThermalError),
    Heat(E),
}
impl<E> From<MotorThermalError> for MotorHeatError<E> {
    fn from(error: MotorThermalError) -> Self { Self::Motor(error) }
}
impl<E: core::fmt::Debug> core::fmt::Display for MotorHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "motor heat assembly: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for MotorHeatError<E> {}

type Error<E> = CircuitEnthalpyError<MotorHeatError<E>>;

fn heat_error<E: From<HeatNetworkError>>(error: HeatNetworkError) -> Error<E> {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(MotorHeatError::Heat(other.into())),
    }
}

/// Native motor receivers on one complete, already prepared heat assembly.
///
/// Cells and their original caloric/material bindings come from that assembly,
/// never from a second parallel network. The winding and bridge keep explicit
/// receiving cells; they may coincide only under the native homogenization
/// assumption. Graph resistors, electrical contacts and circuit-basis bindings
/// refuse because this motor has no owner for their work.
///
/// Keep the same assembly across steps to preserve thermal source selections.
/// The accepted motor retains its winding claim and native rotor factorization.
pub struct MotorHeatNetwork<'a, A: HeatAssembly + ?Sized> {
    pub assembly: &'a A,
    pub winding: MotorWinding,
    pub bridge_cell: CellId,
}

impl<A: HeatAssembly + ?Sized> MotorHeatNetwork<'_, A> {
    fn network(&self, cx: &Cx<'_>) -> Result<MotorThermalNetwork<'_>, Error<A::Error>> {
        cx.checkpoint().map_err(|_| CircuitEnthalpyError::Cancelled)?;
        let prepared = self.assembly.prepared_heat_network();
        let bulk = prepared.thermal_material_network().map_err(heat_error)?;
        prepared.admit_field_independent_boundaries(cx, &[]).map_err(heat_error)?;
        Ok(MotorThermalNetwork {
            cells: bulk.cells,
            links: bulk.links,
            winding: self.winding.clone(),
            bridge_cell: self.bridge_cell,
        })
    }

    fn resolve(
        &self,
        cx: &Cx<'_>,
        network: &MotorThermalNetwork<'_>,
        thermal: &DepositionSimulation<'_>,
        original: Option<&ResolvedConductor>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<(CircuitThermalCoefficients<A::Sample>, ResolvedConductor), Error<A::Error>> {
        let coefficients = self.assembly.resolve_heat(cx, thermal)
            .map_err(|error| error.map_source(MotorHeatError::Heat))?;
        // A provider may not return extra electrical ports even when its
        // prepared topology declares no such owner.
        if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
            return Err(heat_error(HeatNetworkError::Invalid(
                "motor heat assembly cannot discard graph resistor or electrical-contact work",
            )));
        }
        let winding = network.resolve_winding(
            thermal, &coefficients.bulk.source.as_ref().material, original, cancel,
        ).map_err(motor_error)?;
        Ok((coefficients.bulk, winding))
    }

    /// Bind the original held-current motor without advancing any history.
    pub fn bind<'r>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        native: &'r RotorMidpoint<'r>,
        state: RotorTransientState,
        parameters: RotorCurrentDrive,
    ) -> Result<MaterialRotorCurrentDrive<'r>, Error<A::Error>> {
        let network = self.network(cx)?;
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        network.bind_resolved(thermal, native, state, parameters,
            || self.resolve(cx, &network, thermal, None, &cancel).map(|(_, winding)| winding),
            &cancel)
    }

    /// Bind the existing finite-inductance motor, retaining current and motion.
    pub fn bind_voltage<'step, 'model>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        native: &RotorVoltageMotor<'step, 'model>,
        state: RotorWindingState,
    ) -> Result<MaterialRotorVoltageDrive<'step, 'model>, Error<A::Error>> {
        let network = self.network(cx)?;
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        network.bind_voltage_resolved(thermal, native, state,
            || self.resolve(cx, &network, thermal, None, &cancel).map(|(_, winding)| winding),
            &cancel)
    }

    /// Solve current-regulated actuation and the COMPLETE thermal assembly on
    /// the same proposed endpoint. Thermal contacts, integrated layers,
    /// radiation and finite coolant laws reach R(T), source curtailment and
    /// both native heat receivers before joint acceptance.
    ///
    /// Full source results remain in report.source.network. Closed-radiation
    /// roundoff is internal: external heat is report.non_electrical_heat_j
    /// minus duration * network.as_ref().radiation_net_power_w(). Its correction
    /// is independently admitted against the original total-energy tolerance.
    /// No correction changes stored energy or deposits another heat source.
    ///
    /// Native current-drive assumptions and first-order thermal accuracy remain.
    /// Geometry, view factors, original caloric laws and mechanical coefficients
    /// stay fixed; no field-dependent thermal boundary or new motor physics is
    /// inferred. Source, work-budget, solver and cancellation refusals publish
    /// neither thermal nor mechanical history.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t, 'r>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        motor: &MaterialRotorCurrentDrive<'r>,
        external_force: &[f64],
        commanded_current_a: f64,
        power: RotorDrivePowerControl,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<MotorEnthalpyCandidate<'t, 'r, A::Sample>, Error<A::Error>> {
        cx.checkpoint().map_err(|_| CircuitEnthalpyError::Cancelled)?;
        admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        let network = self.network(cx)?;
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        let candidate = network.step_resolved_candidate(
            thermal, target_s, motor, external_force, commanded_current_a, power,
            controls, initial_energy_increments_j,
            |state, original| self.resolve(cx, &network, state, original, &cancel),
            &cancel,
        )?;
        closed_balance(
            target_s - thermal.time_s(), candidate.report().total_energy_residual_j,
            candidate.report().source.network.as_ref(), controls.total_energy_tolerance_j,
        ).map_err(heat_error)?;
        cx.checkpoint().map_err(|_| CircuitEnthalpyError::Cancelled)?;
        Ok(candidate)
    }

    /// The same complete transport on the finite-inductance voltage motor.
    /// Native current, magnetic storage, source regeneration and cached rotor
    /// motion share the original electrothermal residual and heat routing.
    /// The source ownership, closed-radiation account and refusal guarantees
    /// of step_candidate also apply.
    #[allow(clippy::too_many_arguments)]
    pub fn step_voltage_candidate<'t, 'step, 'model>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        motor: &MaterialRotorVoltageDrive<'step, 'model>,
        external_force: &[f64],
        commanded_voltage_v: f64,
        power: RotorVoltageControl,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<RotorVoltageEnthalpyCandidate<'t, 'step, 'model, A::Sample>, Error<A::Error>> {
        cx.checkpoint().map_err(|_| CircuitEnthalpyError::Cancelled)?;
        admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        let network = self.network(cx)?;
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        let candidate = network.step_voltage_resolved_candidate(
            thermal, target_s, motor, external_force, commanded_voltage_v, power,
            controls, initial_energy_increments_j,
            |state, original| self.resolve(cx, &network, state, original, &cancel),
            &cancel,
        )?;
        closed_balance(
            target_s - thermal.time_s(), candidate.report().total_energy_residual_j,
            candidate.report().source.network.as_ref(), controls.total_energy_tolerance_j,
        ).map_err(heat_error)?;
        cx.checkpoint().map_err(|_| CircuitEnthalpyError::Cancelled)?;
        Ok(candidate)
    }
}
