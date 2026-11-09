//! Finite gate/Miller charge with original finite battery, motor and heat owners.
//!
//! Fifteen authored coordinates retain grouped node and gate-charge histories.
//! Intrinsic Qds, Cgs and Cgd stay fixed; source-bound channel current, driver resistances and
//! other dissipative sources refresh at each actual thermal proposal. Signed
//! gate-port work is external supply/recovery; driver heat has its own authored
//! receiver. No independent gate-supply reservoir is inferred.
//!
//! State order is battery chemical/polarization/bus charges; converter flux,
//! two bus charges, two grouped node charges and four gate charges; then motor
//! flux and rotor/stator momenta. All fifteen coordinates are continuous.

use super::super::gate::resolve_gate_material;
use super::*;
use fs_phs::circuit::converter::commutation::capacitance::battery_motor::{
    CapacitiveBatteryMotorCycleStep, PreparedGateBatteryMotorCycle,
};
use fs_phs::circuit::converter::commutation::capacitance::gate::GateCommutationParameters;
use fs_phs::circuit::converter::commutation::gate::GateChannelResponse;

type GateSources<C = GateConductanceMaterial> = Arc<GateMaterial<C>>;

/// Battery three, converter nine, then motor three. Converter node coordinates
/// include their high-side gate electrode; four gate charges remain independent.
/// Incoming grouped charges must already use the complete intrinsic storage.
pub type MaterialGateBatteryMotorDrive<C = GateConductanceMaterial> =
    MaterialCommutatedBatteryMotorDrive<15, Arc<GateMaterial<C>>>;

impl<C: GateChannelMaterial> MaterialGateBatteryMotorDrive<C> {
    pub fn gate_material(&self) -> &GateMaterial<C> {
        &self.capacitors
    }
    pub fn gate_parameters(
        &self,
        maximum_cells_per_chord: usize,
    ) -> GateCommutationParameters<C::Response> {
        parameters(self.parameters, &self.capacitors, maximum_cells_per_chord)
    }
}

impl<C: GateChannelMaterial> DriveSources for GateSources<C> {
    fn same_electrical_coefficients(&self, other: &Self) -> bool {
        self.as_ref().same_electrical_coefficients(other.as_ref())
    }
    fn gate_driver_cells(&self) -> Option<[CellId; 4]> {
        Some(self.driver_cells())
    }
}

/// Original material/caloric source plus intrinsic storage, gate response and driver sources.
#[derive(Clone, Debug)]
pub struct GateBatteryMotorThermalSample<S = MaterialNetworkSample, C = GateConductanceMaterial> {
    pub drive: CommutatedBatteryMotorThermalSample<S>,
    pub gate: Arc<GateMaterial<C>>,
    /// Exact chord traversal budget for this proposal; not a physical coefficient.
    pub maximum_cells_per_chord: usize,
}

#[must_use]
pub struct GateBatteryMotorEnthalpyCandidate<
    't,
    S = MaterialNetworkSample,
    C = GateConductanceMaterial,
> {
    thermal: DepositionSimulation<'t>,
    drive: MaterialGateBatteryMotorDrive<C>,
    report:
        CircuitEnthalpyReport<GateBatteryMotorThermalSample<S, C>, CapacitiveBatteryMotorCycleStep>,
}
impl<'t, S, C> GateBatteryMotorEnthalpyCandidate<'t, S, C> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn drive(&self) -> &MaterialGateBatteryMotorDrive<C> {
        &self.drive
    }
    pub fn report(
        &self,
    ) -> &CircuitEnthalpyReport<GateBatteryMotorThermalSample<S, C>, CapacitiveBatteryMotorCycleStep>
    {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        MaterialGateBatteryMotorDrive<C>,
        CircuitEnthalpyReport<GateBatteryMotorThermalSample<S, C>, CapacitiveBatteryMotorCycleStep>,
    ) {
        (self.thermal, self.drive, self.report)
    }
}

fn parameters<C: GateChannelMaterial>(
    converter: CommutationParameters,
    material: &GateMaterial<C>,
    maximum_cells_per_chord: usize,
) -> GateCommutationParameters<C::Response> {
    GateCommutationParameters {
        converter,
        intrinsic_drain_charge_curves: core::array::from_fn(|i| {
            material.intrinsic_drain_charge()[i].curve().clone()
        }),
        gate_source_capacitance_f: core::array::from_fn(|i| {
            material.gate_source_capacitance()[i].capacitance_f()
        }),
        gate_drain_capacitance_f: core::array::from_fn(|i| {
            material.gate_drain_capacitance()[i].capacitance_f()
        }),
        channel_conductance: core::array::from_fn(|i| {
            material.channel_conductance()[i].native_response()
        }),
        gate_drivers: *material.drivers(),
        maximum_cells_per_chord,
    }
}

fn admit_gate_state<C: GateChannelMaterial, P: GateChannelResponse + Clone + 'static>(
    native: &BatteryMotorDrive,
    state: &[f64; 15],
    authored: &GateCommutationParameters<P>,
    sources: &GateMaterial<C>,
    converter: CommutationParameters,
    maximum_cells_per_chord: usize,
) -> Result<(), Error> {
    authored.validate().map_err(model)?;
    let resolved = parameters(converter, sources, maximum_cells_per_chord);
    if resolved.intrinsic_drain_charge_curves != authored.intrinsic_drain_charge_curves
        || resolved.gate_source_capacitance_f.map(f64::to_bits)
            != authored.gate_source_capacitance_f.map(f64::to_bits)
        || resolved.gate_drain_capacitance_f.map(f64::to_bits)
            != authored.gate_drain_capacitance_f.map(f64::to_bits)
        || converter.inductance_h.to_bits() != authored.converter.inductance_h.to_bits()
        || converter.bus_capacitance_f.map(f64::to_bits)
            != authored.converter.bus_capacitance_f.map(f64::to_bits)
    {
        return Err(invalid(
            "incoming gate/electrode charge must retain its authored intrinsic Qds, Cgs, Cgd and bus storage",
        ));
    }
    let converter_state = core::array::from_fn(|i| state[3 + i]);
    let converter_energy = resolved.stored_energy_j(&converter_state).map_err(model)?;
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
            motor.system().hamiltonian(&state[12..]),
        ])?
        .to_f64(),
    )?;
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn advance_gate<C: GateChannelMaterial>(
    original: BatteryMotorDrive,
    devices: [SwitchDevice; 4],
    material: &GateMaterial<C>,
    maximum_cells_per_chord: usize,
    state: &[f64; 15],
    start_s: f64,
    target_s: f64,
    intervals: &[CommutationInterval],
    commands: &[BatteryMotorSwitchCommand],
    substeps: usize,
    electrical_controls: PowerLoadControls,
    cycle_controls: SwitchCycleControls,
    cancel: &RefCell<impl FnMut() -> bool>,
) -> Result<([f64; 15], CapacitiveBatteryMotorCycleStep), Error> {
    let p = original.electrical().parameters();
    let gate_parameters = parameters(
        CommutationParameters {
            inductance_h: p.inductance_h,
            bus_capacitance_f: p.bus_capacitance_f,
            winding_resistance_ohm: p.winding_resistance_ohm,
            devices,
        },
        material,
        maximum_cells_per_chord,
    );
    let mut native = PreparedGateBatteryMotorCycle::<C::Response>::new_gate_controlled(
        original,
        gate_parameters,
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
    let mut next = [0.0; 15];
    let step = native
        .step_interval_commands_controlled(state, commands, &mut next, || (*cancel.borrow_mut())())
        .map_err(electrical)?;
    Ok((next, step))
}

fn source_sample<S, C>(
    network: S,
    resistors: Arc<[ConverterChannelMaterial; DRIVE_RESISTORS]>,
    devices: Arc<[DeviceConductionMaterial; 4]>,
    gate: GateSources<C>,
    load_heat: Option<RegulatedLoadHeat>,
    motor: DcMotorParameters,
    bearing_cell: CellId,
    maximum_cells_per_chord: usize,
) -> GateBatteryMotorThermalSample<S, C> {
    GateBatteryMotorThermalSample {
        drive: CommutatedBatteryMotorThermalSample {
            network,
            resistors,
            devices,
            bearing_cell,
            bearing_drag_nm_s_per_rad: motor.bearing_drag_nm_s_per_rad,
            load_heat,
        },
        gate,
        maximum_cells_per_chord,
    }
}

/// Original finite material network with explicit intrinsic gate storage and
/// four independently received driver losses. Authored on/off intervals select
/// driver voltage/resistance branches; actual Vds and Vgs drive the channel law.
pub struct GateBatteryMotorThermalNetwork<'a, B = GateConductanceBinding<'a>> {
    pub drive: CommutatedBatteryMotorThermalNetwork<'a>,
    pub gate: GateMaterialBindings<'a, B>,
    /// Bounded exact native traversal; exhaustion refuses without resampling.
    pub maximum_cells_per_chord: usize,
}
impl<B: GateChannelBinding> GateBatteryMotorThermalNetwork<'_, B> {
    fn resolve_gates(
        &self,
        thermal: &DepositionSimulation<'_>,
        sample: &MaterialNetworkSample,
        retained: Option<&MaterialGateBatteryMotorDrive<B::Material>>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<GateSources<B::Material>, Error> {
        resolve_gate_material(
            self.drive.cells,
            self.drive.devices,
            &self.gate,
            thermal,
            sample,
            retained.map(|old| old.capacitors.as_ref()),
            cancel,
        )
        .map(Arc::new)
    }

    /// Incoming bus, grouped node and gate charges must already use the authored full
    /// intrinsic decomposition. The native initializer takes explicit bus/node/Vgs
    /// voltages. Material binding never resets a gate or reinterprets Coss as Qds.
    pub fn bind<P: GateChannelResponse + Clone + 'static>(
        &self,
        thermal: &DepositionSimulation<'_>,
        native: &BatteryMotorDrive,
        authored: &GateCommutationParameters<P>,
        state: [f64; 15],
        cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialGateBatteryMotorDrive<B::Material>, Error> {
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
                let (r, d) = self
                    .drive
                    .resolve_bindings::<15, GateSources<B::Material>>(
                        thermal,
                        &coefficients.source,
                        fixed,
                        None,
                        &cancel,
                    )?;
                let c = self.resolve_gates(thermal, &coefficients.source, None, &cancel)?;
                Ok((r, d, c))
            },
            |state, c, p| {
                admit_gate_state(native, state, authored, c, p, self.maximum_cells_per_chord)
            },
            &cancel,
        )
    }

    /// All fifteen coordinates remain continuous through commands and thermal
    /// boundaries. Four gate works and four driver heats remain separate from
    /// the original fourteen losses and all reversible capacitor storage.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        drive: &MaterialGateBatteryMotorDrive<B::Material>,
        intervals: &[CommutationInterval],
        commands: &[BatteryMotorSwitchCommand],
        substeps_per_interval: usize,
        load_heat_cell: Option<CellId>,
        electrical_controls: PowerLoadControls,
        cycle_controls: SwitchCycleControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<GateBatteryMotorEnthalpyCandidate<'t, MaterialNetworkSample, B::Material>, Error>
    {
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
                let c = self.resolve_gates(trial, &coefficients.source, retained, &cancel)?;
                Ok((coefficients, r, d, c))
            },
            |original, devices, c, state| {
                advance_gate(
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
        Ok(GateBatteryMotorEnthalpyCandidate {
            thermal,
            drive,
            report,
        })
    }
}

/// The same fifteen-state source-bound drive inside the complete existing heat
/// assembly, including finite coolant, contact and radiation inventories.
pub struct GateBatteryMotorHeatNetwork<'a, A: HeatAssembly + ?Sized, B = GateConductanceBinding<'a>>
{
    pub drive: CommutatedBatteryMotorHeatNetwork<'a, A>,
    pub gate: GateMaterialBindings<'a, B>,
    /// Bounded exact native traversal; exhaustion refuses without resampling.
    pub maximum_cells_per_chord: usize,
}
impl<A: HeatAssembly + ?Sized, B: GateChannelBinding> GateBatteryMotorHeatNetwork<'_, A, B> {
    pub fn bind<P: GateChannelResponse + Clone + 'static>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        native: &BatteryMotorDrive,
        authored: &GateCommutationParameters<P>,
        state: [f64; 15],
    ) -> Result<MaterialGateBatteryMotorDrive<B::Material>, AssemblyError<A::Error>> {
        let network = GateBatteryMotorThermalNetwork {
            drive: self.drive.network(cx)?,
            gate: self.gate,
            maximum_cells_per_chord: self.maximum_cells_per_chord,
        };
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        network.drive.bind_resolved(
            thermal,
            native,
            state,
            |fixed| {
                let (coefficients, r, d) = self.drive.resolve::<15, GateSources<B::Material>>(
                    cx,
                    &network.drive,
                    thermal,
                    fixed,
                    None,
                    &cancel,
                )?;
                let c = network
                    .resolve_gates(
                        thermal,
                        &coefficients.source.as_ref().material,
                        None,
                        &cancel,
                    )
                    .map_err(map_drive)?;
                Ok((r, d, c))
            },
            |state, c, p| {
                admit_gate_state(native, state, authored, c, p, self.maximum_cells_per_chord)
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
        drive: &MaterialGateBatteryMotorDrive<B::Material>,
        intervals: &[CommutationInterval],
        commands: &[BatteryMotorSwitchCommand],
        substeps_per_interval: usize,
        load_heat_cell: Option<CellId>,
        electrical_controls: PowerLoadControls,
        cycle_controls: SwitchCycleControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<
        GateBatteryMotorEnthalpyCandidate<'t, A::Sample, B::Material>,
        AssemblyError<A::Error>,
    > {
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        crate::circuit_enthalpy::admit_driver_unknowns(
            thermal,
            controls,
            initial_energy_increments_j,
        )?;
        let network = GateBatteryMotorThermalNetwork {
            drive: self.drive.network(cx)?,
            gate: self.gate,
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
                    .resolve_gates(
                        trial,
                        &coefficients.source.as_ref().material,
                        retained,
                        &cancel,
                    )
                    .map_err(map_drive)?;
                Ok((coefficients, r, d, c))
            },
            |original, devices, c, state| {
                advance_gate(
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
        Ok(GateBatteryMotorEnthalpyCandidate {
            thermal,
            drive,
            report,
        })
    }
}
