//! Explicit constant device-capacitor storage, finite source/motor and material heat.
//!
//! Every incoming converter charge is the original coupled bus/node electrode
//! charge. Binding requires its authored capacitances to match the four original
//! material claims exactly. Actual caloric proposals refresh claim validity and
//! dissipative laws, while keeping every storage coefficient fixed. Temperature-
//! or voltage-dependent capacitance requires reversible energy/entropy ownership
//! and is refused by this lane.

use super::super::capacitors::resolve_capacitances;
use super::*;
use fs_phs::circuit::converter::commutation::capacitance::CapacitiveCommutationParameters;
use fs_phs::circuit::converter::commutation::capacitance::battery_motor::{
    CapacitiveBatteryMotorCycleStep, PreparedCapacitiveBatteryMotorCycle,
};

type CapacitorSources = Arc<[DeviceCapacitanceMaterial; 4]>;

/// Eleven original coordinates, including both independent switching-node charges.
/// Battery chemical/polarization/bus, converter flux/bus/bus/node/node, then
/// motor flux and both shaft momenta. No nine-state promotion resets the nodes.
pub type MaterialCapacitiveBatteryMotorDrive =
    MaterialCommutatedBatteryMotorDrive<11, Arc<[DeviceCapacitanceMaterial; 4]>>;

impl MaterialCapacitiveBatteryMotorDrive {
    pub fn capacitances(&self) -> &[DeviceCapacitanceMaterial; 4] {
        &self.capacitors
    }
    pub fn capacitive_parameters(&self) -> CapacitiveCommutationParameters {
        CapacitiveCommutationParameters {
            converter: self.parameters,
            output_capacitance_f: values(&self.capacitors),
        }
    }
}

/// The existing material/caloric source plus four original constant charge laws.
#[derive(Clone, Debug)]
pub struct CapacitiveBatteryMotorThermalSample<S = MaterialNetworkSample> {
    pub drive: CommutatedBatteryMotorThermalSample<S>,
    pub capacitances: Arc<[DeviceCapacitanceMaterial; 4]>,
}

#[must_use]
pub struct CapacitiveBatteryMotorEnthalpyCandidate<'t, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'t>,
    drive: MaterialCapacitiveBatteryMotorDrive,
    report: CircuitEnthalpyReport<
        CapacitiveBatteryMotorThermalSample<S>,
        CapacitiveBatteryMotorCycleStep,
    >,
}
impl<'t, S> CapacitiveBatteryMotorEnthalpyCandidate<'t, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn drive(&self) -> &MaterialCapacitiveBatteryMotorDrive {
        &self.drive
    }
    pub fn report(
        &self,
    ) -> &CircuitEnthalpyReport<
        CapacitiveBatteryMotorThermalSample<S>,
        CapacitiveBatteryMotorCycleStep,
    > {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        MaterialCapacitiveBatteryMotorDrive,
        CircuitEnthalpyReport<
            CapacitiveBatteryMotorThermalSample<S>,
            CapacitiveBatteryMotorCycleStep,
        >,
    ) {
        (self.thermal, self.drive, self.report)
    }
}

fn values(sources: &[DeviceCapacitanceMaterial; 4]) -> [f64; 4] {
    core::array::from_fn(|i| sources[i].capacitance_f())
}

impl DriveReceipt for CapacitiveBatteryMotorCycleStep {
    fn account(&self) -> DrivePower<'_> {
        DrivePower {
            native: &self.native,
            storage_change_j: self.storage_change_j,
            heat: self.heat_joules(),
            gate: self.converter.gate.as_ref().map(|gate| GateTransfers {
                source_work_j: gate.source_work_j,
                resistor_heat_j: gate.resistor_heat_j,
            }),
            battery_terminal_work_j: self.battery_terminal_work_j,
            converter_terminal_work_j: self.converter.terminal_work_j,
            motor_electrical_work_j: self.motor_electrical_work_j,
            electromagnetic_work_j: self.electromagnetic_work_j,
            shaft_work_j: self.shaft_work_j,
            load_work_j: self.load_work_j,
        }
    }
    fn admit_storage(&self, allowance_j: f64) -> Result<(), Error> {
        // Four signed capacitor energy changes are storage, never package heat.
        // Independently check the native breakdown before the common seven-owner
        // and fourteen-dissipator physical balances.
        let c = &self.converter;
        let mut components = sum([c.inductor_energy_change_j]
            .into_iter()
            .chain(c.bus_capacitor_energy_change_j)
            .chain(c.device_capacitor_energy_change_j))?;
        if let Some(gate) = &c.gate {
            components = components
                + sum(gate.gate_source_energy_change_j)?
                + sum(gate.gate_drain_energy_change_j)?;
        }
        let stored = Dd::from_f64(finite(c.storage_change_j)?);
        let grouped = Dd::from_f64(finite(self.storage_change_j[3])?);
        for residual in [components - stored, stored - grouped] {
            if !residual.hi.is_finite()
                || !residual.lo.is_finite()
                || residual.abs().to_f64() > allowance_j
            {
                return Err(invalid(
                    "device capacitor storage must remain in the original converter energy account",
                ));
            }
        }
        Ok(())
    }
}

fn admit_capacitive_state(
    native: &BatteryMotorDrive,
    state: &[f64; 11],
    authored_output_capacitance_f: [f64; 4],
    sources: &[DeviceCapacitanceMaterial; 4],
    converter: CommutationParameters,
) -> Result<(), Error> {
    let resolved = values(sources);
    if resolved.map(f64::to_bits) != authored_output_capacitance_f.map(f64::to_bits) {
        return Err(invalid(
            "incoming device charge must retain its authored original capacitances",
        ));
    }
    let parameters = CapacitiveCommutationParameters {
        converter,
        output_capacitance_f: resolved,
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
fn advance_capacitive(
    original: BatteryMotorDrive,
    devices: [SwitchDevice; 4],
    capacitors: &[DeviceCapacitanceMaterial; 4],
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
    let mut native = PreparedCapacitiveBatteryMotorCycle::new_controlled(
        original,
        devices,
        values(capacitors),
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
    capacitances: CapacitorSources,
    load_heat: Option<RegulatedLoadHeat>,
    motor: DcMotorParameters,
    bearing_cell: CellId,
) -> CapacitiveBatteryMotorThermalSample<S> {
    CapacitiveBatteryMotorThermalSample {
        drive: CommutatedBatteryMotorThermalSample {
            network,
            resistors,
            devices,
            bearing_cell,
            bearing_drag_nm_s_per_rad: motor.bearing_drag_nm_s_per_rad,
            load_heat,
        },
        capacitances,
    }
}

/// Original finite material network with an explicit capacitor claim for each
/// corresponding semiconductor package. Existing resistance and heat maps apply.
pub struct CapacitiveBatteryMotorThermalNetwork<'a> {
    pub drive: CommutatedBatteryMotorThermalNetwork<'a>,
    pub capacitances: &'a [DeviceCapacitanceBinding<'a>; 4],
}
impl CapacitiveBatteryMotorThermalNetwork<'_> {
    fn resolve_capacitors(
        &self,
        thermal: &DepositionSimulation<'_>,
        sample: &MaterialNetworkSample,
        retained: Option<&MaterialCapacitiveBatteryMotorDrive>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<CapacitorSources, Error> {
        if self
            .capacitances
            .iter()
            .zip(self.drive.devices)
            .any(|(capacitor, device)| capacitor.cell != device.cell)
        {
            return Err(invalid(
                "each output capacitor must retain its corresponding original semiconductor package",
            ));
        }
        resolve_capacitances(
            self.drive.cells,
            self.capacitances,
            thermal,
            sample,
            retained.map(|old| old.capacitors.as_ref()),
            cancel,
        )
        .map(Arc::new)
    }

    /// Charges must already use the authored capacitance matrix. Initial switch
    /// node voltages are explicit in the native initializer; none are inferred.
    pub fn bind(
        &self,
        thermal: &DepositionSimulation<'_>,
        native: &BatteryMotorDrive,
        authored_output_capacitance_f: [f64; 4],
        state: [f64; 11],
        cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialCapacitiveBatteryMotorDrive, Error> {
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
                let (r, d) = self.drive.resolve_bindings::<11, CapacitorSources>(
                    thermal,
                    &coefficients.source,
                    fixed,
                    None,
                    &cancel,
                )?;
                let c = self.resolve_capacitors(thermal, &coefficients.source, None, &cancel)?;
                Ok((r, d, c))
            },
            |state, c, p| {
                admit_capacitive_state(native, state, authored_output_capacitance_f, c, p)
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
        drive: &MaterialCapacitiveBatteryMotorDrive,
        intervals: &[CommutationInterval],
        commands: &[BatteryMotorSwitchCommand],
        substeps_per_interval: usize,
        load_heat_cell: Option<CellId>,
        electrical_controls: PowerLoadControls,
        cycle_controls: SwitchCycleControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<CapacitiveBatteryMotorEnthalpyCandidate<'t>, Error> {
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
                let c = self.resolve_capacitors(trial, &coefficients.source, retained, &cancel)?;
                Ok((coefficients, r, d, c))
            },
            |original, devices, c, state| {
                advance_capacitive(
                    original,
                    devices,
                    c,
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
                source_sample(network, r, d, c, load, motor, self.drive.bearing_cell)
            },
            &cancel,
        )?;
        Ok(CapacitiveBatteryMotorEnthalpyCandidate {
            thermal,
            drive,
            report,
        })
    }
}

/// The same eleven-state source-bound drive inside the complete existing heat
/// assembly, including finite coolant, contact and radiation inventories.
pub struct CapacitiveBatteryMotorHeatNetwork<'a, A: HeatAssembly + ?Sized> {
    pub drive: CommutatedBatteryMotorHeatNetwork<'a, A>,
    pub capacitances: &'a [DeviceCapacitanceBinding<'a>; 4],
}
impl<A: HeatAssembly + ?Sized> CapacitiveBatteryMotorHeatNetwork<'_, A> {
    pub fn bind(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        native: &BatteryMotorDrive,
        authored_output_capacitance_f: [f64; 4],
        state: [f64; 11],
    ) -> Result<MaterialCapacitiveBatteryMotorDrive, AssemblyError<A::Error>> {
        let network = CapacitiveBatteryMotorThermalNetwork {
            drive: self.drive.network(cx)?,
            capacitances: self.capacitances,
        };
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        network.drive.bind_resolved(
            thermal,
            native,
            state,
            |fixed| {
                let (coefficients, r, d) = self.drive.resolve::<11, CapacitorSources>(
                    cx,
                    &network.drive,
                    thermal,
                    fixed,
                    None,
                    &cancel,
                )?;
                let c = network
                    .resolve_capacitors(
                        thermal,
                        &coefficients.source.as_ref().material,
                        None,
                        &cancel,
                    )
                    .map_err(map_drive)?;
                Ok((r, d, c))
            },
            |state, c, p| {
                admit_capacitive_state(native, state, authored_output_capacitance_f, c, p)
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
        drive: &MaterialCapacitiveBatteryMotorDrive,
        intervals: &[CommutationInterval],
        commands: &[BatteryMotorSwitchCommand],
        substeps_per_interval: usize,
        load_heat_cell: Option<CellId>,
        electrical_controls: PowerLoadControls,
        cycle_controls: SwitchCycleControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<CapacitiveBatteryMotorEnthalpyCandidate<'t, A::Sample>, AssemblyError<A::Error>>
    {
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        crate::circuit_enthalpy::admit_driver_unknowns(
            thermal,
            controls,
            initial_energy_increments_j,
        )?;
        let network = CapacitiveBatteryMotorThermalNetwork {
            drive: self.drive.network(cx)?,
            capacitances: self.capacitances,
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
                    .resolve_capacitors(
                        trial,
                        &coefficients.source.as_ref().material,
                        retained,
                        &cancel,
                    )
                    .map_err(map_drive)?;
                Ok((coefficients, r, d, c))
            },
            |original, devices, c, state| {
                advance_capacitive(
                    original,
                    devices,
                    c,
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
                source_sample(network, r, d, c, load, motor, self.drive.bearing_cell)
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
        Ok(CapacitiveBatteryMotorEnthalpyCandidate {
            thermal,
            drive,
            report,
        })
    }
}
