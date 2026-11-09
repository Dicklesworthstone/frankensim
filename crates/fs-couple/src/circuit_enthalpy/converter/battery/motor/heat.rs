//! Complete source-owned cooling on the finite battery-to-shaft drive.
//!
//! The original heat assembly supplies contacts, radiation, integrated layers,
//! reservoir exchange and finite coolant links at every proposed endpoint.
//! Native nine-state integration, resistance feedback and all ten physical heat
//! channels are shared with the ordinary material-network entry point.

use super::*;
use crate::circuit_enthalpy::contacts::heat::{HeatAssembly, HeatNetworkError, closed_balance};
use fs_exec::Cx;

/// Preserve native drive/material errors and the complete heat owner's refusal.
#[derive(Debug)]
pub enum BatteryMotorHeatError<E> {
    Drive(ConverterThermalError),
    Heat(E),
}
impl<E> From<ConverterThermalError> for BatteryMotorHeatError<E> {
    fn from(error: ConverterThermalError) -> Self {
        Self::Drive(error)
    }
}
impl<E: core::fmt::Debug> core::fmt::Display for BatteryMotorHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "battery motor heat assembly: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for BatteryMotorHeatError<E> {}

type Error<E> = CircuitEnthalpyError<BatteryMotorHeatError<E>>;
fn heat_error<E: From<HeatNetworkError>>(error: HeatNetworkError) -> Error<E> {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(BatteryMotorHeatError::Heat(other.into())),
    }
}

/// One prepared complete heat assembly and ten explicit drive-loss receivers.
///
/// Nine resistance channels have the original ordering and source semantics of
/// `BatteryMotorThermalNetwork`. Fixed mechanical bearing drag has its own
/// receiver; it may share a body under an explicitly authored homogenization.
/// All cells and topology come from this assembly, including finite coolant
/// inventory. Keep the same assembly and source models across accepted steps.
/// Graph resistors and electrical contacts refuse because this drive owns no
/// additional electrical state or work for them.
pub struct BatteryMotorHeatNetwork<'a, A: HeatAssembly + ?Sized> {
    pub assembly: &'a A,
    pub channels: &'a [ConverterThermalChannel<'a>; BATTERY_MOTOR_RESISTIVE_CHANNELS],
    pub bearing_cell: CellId,
}

impl<A: HeatAssembly + ?Sized> BatteryMotorHeatNetwork<'_, A> {
    fn network(&self, cx: &Cx<'_>) -> Result<BatteryMotorThermalNetwork<'_>, Error<A::Error>> {
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        let prepared = self.assembly.prepared_heat_network();
        let bulk = prepared.thermal_material_network().map_err(heat_error)?;
        prepared
            .admit_field_independent_boundaries(cx, &[])
            .map_err(heat_error)?;
        Ok(BatteryMotorThermalNetwork {
            cells: bulk.cells,
            links: bulk.links,
            channels: self.channels,
            bearing_cell: self.bearing_cell,
        })
    }

    fn resolve(
        &self,
        cx: &Cx<'_>,
        network: &BatteryMotorThermalNetwork<'_>,
        thermal: &DepositionSimulation<'_>,
        fixed: [f64; BATTERY_MOTOR_RESISTIVE_CHANNELS],
        original: Option<&[ConverterChannelMaterial; BATTERY_MOTOR_RESISTIVE_CHANNELS]>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<
        (
            CircuitThermalCoefficients<A::Sample>,
            [ConverterChannelMaterial; BATTERY_MOTOR_RESISTIVE_CHANNELS],
        ),
        Error<A::Error>,
    > {
        let coefficients = self
            .assembly
            .resolve_heat(cx, thermal)
            .map_err(|error| error.map_source(BatteryMotorHeatError::Heat))?;
        if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
            return Err(heat_error(HeatNetworkError::Invalid(
                "battery motor heat assembly cannot discard graph resistor or electrical-contact work",
            )));
        }
        let channels = resolve_channels(
            network.cells,
            self.channels,
            thermal,
            &coefficients.bulk.source.as_ref().material,
            fixed,
            original,
            false,
            cancel,
        )
        .map_err(drive_error)?;
        Ok((coefficients.bulk, channels))
    }

    /// Bind actual drive resistance and complete cooling sources without
    /// advancing any charge, flux, momentum or finite caloric inventory.
    pub fn bind(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        native: &BatteryMotorDrive,
        state: [f64; 9],
    ) -> Result<MaterialBatteryMotorDrive, Error<A::Error>> {
        let network = self.network(cx)?;
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        network.bind_resolved(
            thermal,
            native,
            state,
            |fixed| {
                self.resolve(cx, &network, thermal, fixed, None, &cancel)
                    .map(|(_, channels)| channels)
            },
            &cancel,
        )
    }

    /// Resolve cooling, all nine R(T,phase) sources and the full native drive
    /// in the original finite-cell residual. Finite coolant exchange remains
    /// internal, while prescribed reservoir heat retains its external sign.
    /// Neither is added to the ten electrical/mechanical loss integrals.
    ///
    /// `report.source.network` retains the complete original heat sample.
    /// External heat is `report.non_electrical_heat_j` minus duration times
    /// `network.as_ref().radiation_net_power_w()`: a closed enclosure's numerical
    /// net is internal and receives the existing independent balance check.
    /// No correction changes storage or deposits an additional heat source.
    ///
    /// Fixed storage, bearing drag and motor reciprocity, endpoint resistance
    /// accuracy and all ordinary drive assumptions remain unchanged. Providers
    /// must resolve repeatably without advancing external histories. Any source,
    /// numerical-budget or cancellation refusal publishes no physical history.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        drive: &MaterialBatteryMotorDrive,
        modulation: BridgeModulation,
        shaft_torque_nm: &[f64; 2],
        auxiliary_load: Option<RegulatedPowerLoad>,
        electrical_controls: PowerLoadControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<BatteryMotorEnthalpyCandidate<'t, A::Sample>, Error<A::Error>> {
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        super::super::super::super::admit_driver_unknowns(
            thermal,
            controls,
            initial_energy_increments_j,
        )?;
        let network = self.network(cx)?;
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        let candidate = network.step_resolved_candidate(
            thermal,
            target_s,
            drive,
            modulation,
            controls,
            initial_energy_increments_j,
            1,
            |original, dt, checkpoint| {
                averaged_drive_step(
                    original,
                    &drive.state,
                    dt,
                    shaft_torque_nm,
                    auxiliary_load,
                    electrical_controls,
                    checkpoint,
                )
            },
            |trial, fixed, original| self.resolve(cx, &network, trial, fixed, original, &cancel),
            &cancel,
        )?;
        closed_balance(
            target_s - thermal.time_s(),
            candidate.report().total_energy_residual_j,
            candidate.report().source.network.as_ref(),
            controls.total_energy_tolerance_j,
        )
        .map_err(heat_error)?;
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        Ok(candidate)
    }
    /// Prescribed PWM with the same complete original heat assembly. The native
    /// cycle owns every substep's ten heat/work integrals; original contacts,
    /// radiation and finite coolant exchange remain in the shared caloric solve.
    /// Gate commands are immutable within the candidate and every subdivision
    /// is charged even on exact-resistance trajectory-cache reuse.
    #[allow(clippy::too_many_arguments)]
    pub fn step_switch_candidate<'t>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        drive: &MaterialBatteryMotorDrive,
        intervals: &[SwitchInterval],
        commands: &[BatteryMotorSwitchCommand],
        substeps_per_interval: usize,
        electrical_controls: PowerLoadControls,
        cycle_controls: SwitchCycleControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<
        BatteryMotorEnthalpyCandidate<'t, A::Sample, BatteryMotorSwitchCycleStep>,
        Error<A::Error>,
    > {
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        let integration_steps =
            switched_drive_work(intervals, commands, substeps_per_interval).map_err(drive_error)?;
        super::super::super::super::admit_driver_unknowns(
            thermal,
            controls,
            initial_energy_increments_j,
        )?;
        let network = self.network(cx)?;
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        let candidate = network.step_resolved_candidate(
            thermal,
            target_s,
            drive,
            drive.modulation,
            controls,
            initial_energy_increments_j,
            integration_steps,
            |original, dt, checkpoint| {
                switched_drive_step(
                    original,
                    &drive.state,
                    dt,
                    thermal.time_s(),
                    target_s,
                    intervals,
                    commands,
                    substeps_per_interval,
                    electrical_controls,
                    cycle_controls,
                    checkpoint,
                )
            },
            |trial, fixed, original| self.resolve(cx, &network, trial, fixed, original, &cancel),
            &cancel,
        )?;
        closed_balance(
            target_s - thermal.time_s(),
            candidate.report().total_energy_residual_j,
            candidate.report().source.network.as_ref(),
            controls.total_energy_tolerance_j,
        )
        .map_err(heat_error)?;
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        Ok(candidate)
    }
}
