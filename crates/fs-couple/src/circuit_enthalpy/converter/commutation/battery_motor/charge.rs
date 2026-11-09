//! Source-bound reversible Q(V) device storage and the original finite heat path.
//!
//! All four complete charge laws and their finite voltage support remain fixed
//! along the accepted electrical/thermal history. Actual caloric proposals
//! refresh source applicability and dissipative coefficients. Device energy is
//! the original integral of V dQ, not a sampled C(V) or Q*V/2.
//! Temperature-dependent charge laws need reversible caloric/entropy ownership
//! and are refused here.

use super::super::charge::resolve_charge_curves;
use super::*;
use fs_phs::circuit::capacitance::CapacitorChargeCurve;
use fs_phs::circuit::converter::commutation::capacitance::battery_motor::{
    CapacitiveBatteryMotorCycleStep, PreparedChargeBatteryMotorCycle,
};
use fs_phs::circuit::converter::commutation::capacitance::charge::ChargeCommutationParameters;

type ChargeSources = Arc<[DeviceChargeMaterial; 4]>;

/// Eleven original coordinates, including both independent switching-node charges.
/// Battery chemical/polarization/bus, converter flux/bus/bus/node/node, then
/// motor flux and both shaft momenta. No nine-state promotion resets the nodes.
pub type MaterialChargeBatteryMotorDrive =
    MaterialCommutatedBatteryMotorDrive<11, Arc<[DeviceChargeMaterial; 4]>>;

impl MaterialChargeBatteryMotorDrive {
    pub fn charges(&self) -> &[DeviceChargeMaterial; 4] {
        &self.capacitors
    }
    pub fn charge_parameters(&self, maximum_cells_per_chord: usize) -> ChargeCommutationParameters {
        ChargeCommutationParameters {
            converter: self.parameters,
            output_charge_curves: curves(&self.capacitors),
            maximum_cells_per_chord,
        }
    }
}

/// The original material/caloric source plus four complete reversible charge laws.
#[derive(Clone, Debug)]
pub struct ChargeBatteryMotorThermalSample<S = MaterialNetworkSample> {
    pub drive: CommutatedBatteryMotorThermalSample<S>,
    pub charges: Arc<[DeviceChargeMaterial; 4]>,
    /// Exact chord traversal budget for this proposal; not a physical coefficient.
    pub maximum_cells_per_chord: usize,
}

#[must_use]
pub struct ChargeBatteryMotorEnthalpyCandidate<'t, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'t>,
    drive: MaterialChargeBatteryMotorDrive,
    report:
        CircuitEnthalpyReport<ChargeBatteryMotorThermalSample<S>, CapacitiveBatteryMotorCycleStep>,
}
impl<'t, S> ChargeBatteryMotorEnthalpyCandidate<'t, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn drive(&self) -> &MaterialChargeBatteryMotorDrive {
        &self.drive
    }
    pub fn report(
        &self,
    ) -> &CircuitEnthalpyReport<ChargeBatteryMotorThermalSample<S>, CapacitiveBatteryMotorCycleStep>
    {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        MaterialChargeBatteryMotorDrive,
        CircuitEnthalpyReport<ChargeBatteryMotorThermalSample<S>, CapacitiveBatteryMotorCycleStep>,
    ) {
        (self.thermal, self.drive, self.report)
    }
}

fn curves(sources: &[DeviceChargeMaterial; 4]) -> [CapacitorChargeCurve; 4] {
    core::array::from_fn(|i| sources[i].curve().clone())
}

// The original capacitive receipt already implements the shared seven-owner,
// fourteen-loss and separate device-storage checks. Nonlinear storage uses
// exactly that actual native report; there is no second thermal account.

fn admit_charge_state(
    native: &BatteryMotorDrive,
    state: &[f64; 11],
    authored_output_charge_curves: &[CapacitorChargeCurve; 4],
    sources: &[DeviceChargeMaterial; 4],
    converter: CommutationParameters,
    maximum_cells_per_chord: usize,
) -> Result<(), Error> {
    let resolved = curves(sources);
    if resolved != *authored_output_charge_curves {
        return Err(invalid(
            "incoming electrode charge must retain every knot and bound of its authored original Q(V) laws",
        ));
    }
    let parameters = ChargeCommutationParameters {
        converter,
        output_charge_curves: resolved,
        maximum_cells_per_chord,
    };
    let converter_state = core::array::from_fn(|i| state[3 + i]);
    let converter_energy = parameters
        .stored_energy_j(&converter_state)
        .map_err(model)?;
    // Original constituent storage owners validate their own coordinates. No
    // approximate nine-state model or trajectory is used for the new charge law.
    let battery = BatterySource::new(
        native.electrical().battery_curve().clone(),
        native.electrical().battery_parameters(),
    )
    .map_err(model)?;
    let motor = DcMotor::new(native.motor_parameters()).map_err(model)?;
    finite(
        sum([
            battery.system().hamiltonian(&state[..3]),
            converter_energy,
            motor.system().hamiltonian(&state[8..]),
        ])?
        .to_f64(),
    )?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn advance_charge(
    original: BatteryMotorDrive,
    devices: [SwitchDevice; 4],
    capacitors: &[DeviceChargeMaterial; 4],
    maximum_cells_per_chord: usize,
    state: &[f64; 11],
    start_s: f64,
    target_s: f64,
    intervals: &[CommutationInterval],
    commands: &[BatteryMotorSwitchCommand],
    substeps: usize,
    electrical_controls: PowerLoadControls,
    cycle_controls: SwitchCycleControls,
    cancel: &RefCell<impl FnMut() -> bool>,
) -> Result<([f64; 11], CapacitiveBatteryMotorCycleStep), Error> {
    let mut native = PreparedChargeBatteryMotorCycle::new_charge_controlled(
        original,
        devices,
        curves(capacitors),
        maximum_cells_per_chord,
        intervals,
        substeps,
        electrical_controls,
        cycle_controls,
        || (*cancel.borrow_mut())(),
    )
    .map_err(electrical)?;
    admit_cycle_clock(
        start_s,
        target_s,
        native.period_s(),
        native.integration_steps(),
        intervals
            .len()
            .checked_mul(substeps)
            .ok_or_else(|| invalid("commutation work overflow"))?,
    )?;
    let mut next = [0.0; 11];
    let step = native
        .step_interval_commands_controlled(state, commands, &mut next, || (*cancel.borrow_mut())())
        .map_err(electrical)?;
    Ok((next, step))
}

fn source_sample<S>(
    network: S,
    resistors: Arc<[ConverterChannelMaterial; DRIVE_RESISTORS]>,
    devices: Arc<[DeviceConductionMaterial; 4]>,
    charges: ChargeSources,
    load_heat: Option<RegulatedLoadHeat>,
    motor: DcMotorParameters,
    bearing_cell: CellId,
    maximum_cells_per_chord: usize,
) -> ChargeBatteryMotorThermalSample<S> {
    ChargeBatteryMotorThermalSample {
        drive: CommutatedBatteryMotorThermalSample {
            network,
            resistors,
            devices,
            bearing_cell,
            bearing_drag_nm_s_per_rad: motor.bearing_drag_nm_s_per_rad,
            load_heat,
        },
        charges,
        maximum_cells_per_chord,
    }
}

/// Original finite material network with a complete Q(V) claim for each
/// corresponding semiconductor package. Curves follow input high/low then
/// output high/low, with the native drain-to-source voltage and charge signs.
/// Existing resistance and heat maps apply.
pub struct ChargeBatteryMotorThermalNetwork<'a> {
    pub drive: CommutatedBatteryMotorThermalNetwork<'a>,
    pub charges: &'a [DeviceChargeBinding<'a>; 4],
    /// Bounded exact native traversal; exhaustion refuses without resampling.
    pub maximum_cells_per_chord: usize,
}
impl ChargeBatteryMotorThermalNetwork<'_> {
    fn resolve_curves(
        &self,
        thermal: &DepositionSimulation<'_>,
        sample: &MaterialNetworkSample,
        retained: Option<&MaterialChargeBatteryMotorDrive>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<ChargeSources, Error> {
        if self
            .charges
            .iter()
            .zip(self.drive.devices)
            .any(|(capacitor, device)| capacitor.cell != device.cell)
        {
            return Err(invalid(
                "each output charge law must retain its corresponding original semiconductor package",
            ));
        }
        resolve_charge_curves(
            self.drive.cells,
            self.charges,
            thermal,
            sample,
            retained.map(|old| old.capacitors.as_ref()),
            cancel,
        )
        .map(Arc::new)
    }

    /// Incoming electrode charges must already use every authored Q(V) knot.
    /// The native initializer takes explicit bus/node voltages; neither constant
    /// capacitances nor new equilibrium node values replace these coordinates.
    pub fn bind(
        &self,
        thermal: &DepositionSimulation<'_>,
        native: &BatteryMotorDrive,
        authored_output_charge_curves: &[CapacitorChargeCurve; 4],
        state: [f64; 11],
        cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialChargeBatteryMotorDrive, Error> {
        let cancel = RefCell::new(cancelled);
        self.drive.bind_resolved(
            thermal,
            native,
            state,
            |fixed| {
                let network = self.drive.network();
                let prepared = network
                    .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
                    .map_err(material)?;
                let coefficients = prepared
                    .resolve(thermal, || (*cancel.borrow_mut())())
                    .map_err(material)?;
                let (r, d) = self.drive.resolve_bindings::<11, ChargeSources>(
                    thermal,
                    &coefficients.source,
                    fixed,
                    None,
                    &cancel,
                )?;
                let c = self.resolve_curves(thermal, &coefficients.source, None, &cancel)?;
                Ok((r, d, c))
            },
            |state, c, p| {
                admit_charge_state(
                    native,
                    state,
                    authored_output_charge_curves,
                    c,
                    p,
                    self.maximum_cells_per_chord,
                )
            },
            &cancel,
        )
    }

    /// All eleven coordinates are continuous through every gate and thermal
    /// boundary. Four device capacitor storage changes remain in the receipt;
    /// only fourteen physical loss channels and optional actual load work heat.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        drive: &MaterialChargeBatteryMotorDrive,
        intervals: &[CommutationInterval],
        commands: &[BatteryMotorSwitchCommand],
        substeps_per_interval: usize,
        load_heat_cell: Option<CellId>,
        electrical_controls: PowerLoadControls,
        cycle_controls: SwitchCycleControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<ChargeBatteryMotorEnthalpyCandidate<'t>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        crate::circuit_enthalpy::admit_driver_unknowns(
            thermal,
            controls,
            initial_energy_increments_j,
        )?;
        let network = self.drive.network();
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let (thermal, drive, report) = self.drive.step_resolved_candidate(
            thermal,
            target_s,
            drive,
            intervals,
            commands,
            substeps_per_interval,
            load_heat_cell,
            controls,
            initial_energy_increments_j,
            |trial, fixed, retained| {
                let coefficients = prepared
                    .resolve(trial, || (*cancel.borrow_mut())())
                    .map_err(material)?;
                let (r, d) = self.drive.resolve_bindings(
                    trial,
                    &coefficients.source,
                    fixed,
                    retained,
                    &cancel,
                )?;
                let c = self.resolve_curves(trial, &coefficients.source, retained, &cancel)?;
                Ok((coefficients, r, d, c))
            },
            |original, devices, c, state| {
                advance_charge(
                    original,
                    devices,
                    c,
                    self.maximum_cells_per_chord,
                    state,
                    thermal.time_s(),
                    target_s,
                    intervals,
                    commands,
                    substeps_per_interval,
                    electrical_controls,
                    cycle_controls,
                    &cancel,
                )
            },
            |network, r, d, c, load, motor| {
                source_sample(
                    network,
                    r,
                    d,
                    c,
                    load,
                    motor,
                    self.drive.bearing_cell,
                    self.maximum_cells_per_chord,
                )
            },
            &cancel,
        )?;
        Ok(ChargeBatteryMotorEnthalpyCandidate {
            thermal,
            drive,
            report,
        })
    }
}

/// The same eleven-state source-bound drive inside the complete existing heat
/// assembly, including finite coolant, contact and radiation inventories.
pub struct ChargeBatteryMotorHeatNetwork<'a, A: HeatAssembly + ?Sized> {
    pub drive: CommutatedBatteryMotorHeatNetwork<'a, A>,
    pub charges: &'a [DeviceChargeBinding<'a>; 4],
    /// Bounded exact native traversal; exhaustion refuses without resampling.
    pub maximum_cells_per_chord: usize,
}
impl<A: HeatAssembly + ?Sized> ChargeBatteryMotorHeatNetwork<'_, A> {
    pub fn bind(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        native: &BatteryMotorDrive,
        authored_output_charge_curves: &[CapacitorChargeCurve; 4],
        state: [f64; 11],
    ) -> Result<MaterialChargeBatteryMotorDrive, AssemblyError<A::Error>> {
        let network = ChargeBatteryMotorThermalNetwork {
            drive: self.drive.network(cx)?,
            charges: self.charges,
            maximum_cells_per_chord: self.maximum_cells_per_chord,
        };
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        network.drive.bind_resolved(
            thermal,
            native,
            state,
            |fixed| {
                let (coefficients, r, d) = self.drive.resolve::<11, ChargeSources>(
                    cx,
                    &network.drive,
                    thermal,
                    fixed,
                    None,
                    &cancel,
                )?;
                let c = network
                    .resolve_curves(
                        thermal,
                        &coefficients.source.as_ref().material,
                        None,
                        &cancel,
                    )
                    .map_err(map_drive)?;
                Ok((r, d, c))
            },
            |state, c, p| {
                admit_charge_state(
                    native,
                    state,
                    authored_output_charge_curves,
                    c,
                    p,
                    self.maximum_cells_per_chord,
                )
            },
            &cancel,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        drive: &MaterialChargeBatteryMotorDrive,
        intervals: &[CommutationInterval],
        commands: &[BatteryMotorSwitchCommand],
        substeps_per_interval: usize,
        load_heat_cell: Option<CellId>,
        electrical_controls: PowerLoadControls,
        cycle_controls: SwitchCycleControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<ChargeBatteryMotorEnthalpyCandidate<'t, A::Sample>, AssemblyError<A::Error>> {
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        crate::circuit_enthalpy::admit_driver_unknowns(
            thermal,
            controls,
            initial_energy_increments_j,
        )?;
        let network = ChargeBatteryMotorThermalNetwork {
            drive: self.drive.network(cx)?,
            charges: self.charges,
            maximum_cells_per_chord: self.maximum_cells_per_chord,
        };
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        let duration_s = target_s - thermal.time_s();
        let (thermal, drive, report) = network.drive.step_resolved_candidate(
            thermal,
            target_s,
            drive,
            intervals,
            commands,
            substeps_per_interval,
            load_heat_cell,
            controls,
            initial_energy_increments_j,
            |trial, fixed, retained| {
                let (coefficients, r, d) =
                    self.drive
                        .resolve(cx, &network.drive, trial, fixed, retained, &cancel)?;
                let c = network
                    .resolve_curves(
                        trial,
                        &coefficients.source.as_ref().material,
                        retained,
                        &cancel,
                    )
                    .map_err(map_drive)?;
                Ok((coefficients, r, d, c))
            },
            |original, devices, c, state| {
                advance_charge(
                    original,
                    devices,
                    c,
                    self.maximum_cells_per_chord,
                    state,
                    thermal.time_s(),
                    target_s,
                    intervals,
                    commands,
                    substeps_per_interval,
                    electrical_controls,
                    cycle_controls,
                    &cancel,
                )
            },
            |network, r, d, c, load, motor| {
                source_sample(
                    network,
                    r,
                    d,
                    c,
                    load,
                    motor,
                    self.drive.bearing_cell,
                    self.maximum_cells_per_chord,
                )
            },
            &cancel,
        )?;
        closed_balance(
            duration_s,
            report.total_energy_residual_j,
            report.source.drive.network.as_ref(),
            controls.total_energy_tolerance_j,
        )
        .map_err(heat_error)?;
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        Ok(ChargeBatteryMotorEnthalpyCandidate {
            thermal,
            drive,
            report,
        })
    }
}
