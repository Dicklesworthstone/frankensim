//! Finite chemical, electromagnetic, mechanical and caloric storage together.
//!
//! The original nine-state native drive connects converter voltage to motor
//! current reciprocally. Both rotor and stator remain finite; external shaft
//! work and optional auxiliary-load export stay outside the caloric bodies.
//! Every thermal trial starts from the same accepted charges, fluxes and
//! momenta. Nine electrical resistance sources refresh at that trial's actual
//! material state, while fixed bearing drag heats its explicitly named body.
//!
//! OCV, capacitances, inductances, inertias and reciprocal motor coupling remain
//! fixed. Bearing drag is authored in N m s/rad, never read as an ohmic source.
//! Endpoint resistance feedback has first-order coupled accuracy in general.
//! Prescribed complementary PWM uses the same native drive and ten heat paths.
//! No switching-transition loss, iron loss, reversible battery heat, temperature
//! dependent storage, control or protection law is inferred.

/// Source-owned contacts, radiation, reservoir cooling and finite coolant paths.
pub mod heat;

use super::*;
use crate::circuit_enthalpy::CircuitThermalCoefficients;
use fs_phs::circuit::{
    converter::battery_fed::motor::{BatteryMotorDrive, BatteryMotorStep},
    converter::switching::{
        SwitchCycleControls, SwitchInterval,
        battery_motor::{
            BatteryMotorSwitchCommand, BatteryMotorSwitchCycleStep, PreparedBatteryMotorSwitchCycle,
        },
    },
    motor::{DcMotor, DcMotorParameters},
};

/// Converter five, battery series/polarization, lead, then motor winding.
pub const BATTERY_MOTOR_RESISTIVE_CHANNELS: usize = BATTERY_FED_HEAT_CHANNELS + 1;
/// The nine electrical loss channels followed by internal bearing drag.
pub const BATTERY_MOTOR_HEAT_CHANNELS: usize = BATTERY_MOTOR_RESISTIVE_CHANNELS + 1;

/// Explicit physical loss receivers on the original finite thermal network.
/// Receivers may coincide only when that homogenization is authored here.
pub struct BatteryMotorThermalNetwork<'a> {
    pub cells: &'a [MaterialCell<'a>],
    pub links: &'a [ConductionLink],
    pub channels: &'a [ConverterThermalChannel<'a>; BATTERY_MOTOR_RESISTIVE_CHANNELS],
    /// Receiver of the native fixed rotor-to-stator viscous friction loss.
    pub bearing_cell: CellId,
}

/// Accepted battery, converter and finite rotor/stator state with actual sources.
#[derive(Clone, Debug)]
pub struct MaterialBatteryMotorDrive {
    time_s: f64,
    curve: BatteryOcvCurve,
    battery_parameters: BatteryParameters,
    converter_parameters: BuckBoostParameters,
    modulation: BridgeModulation,
    lead_resistance_ohm: f64,
    motor_parameters: DcMotorParameters,
    state: [f64; 9],
    channels: Arc<[ConverterChannelMaterial; BATTERY_MOTOR_RESISTIVE_CHANNELS]>,
    bearing_cell: CellId,
}
impl MaterialBatteryMotorDrive {
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn battery_curve(&self) -> &BatteryOcvCurve {
        &self.curve
    }
    pub const fn battery_parameters(&self) -> BatteryParameters {
        self.battery_parameters
    }
    pub const fn converter_parameters(&self) -> BuckBoostParameters {
        self.converter_parameters
    }
    pub const fn modulation(&self) -> BridgeModulation {
        self.modulation
    }
    pub const fn lead_resistance_ohm(&self) -> f64 {
        self.lead_resistance_ohm
    }
    pub const fn motor_parameters(&self) -> DcMotorParameters {
        self.motor_parameters
    }
    pub const fn state(&self) -> &[f64; 9] {
        &self.state
    }
    pub fn channels(&self) -> &[ConverterChannelMaterial; BATTERY_MOTOR_RESISTIVE_CHANNELS] {
        &self.channels
    }
    pub const fn bearing_cell(&self) -> CellId {
        self.bearing_cell
    }
}

#[derive(Clone, Debug)]
pub struct BatteryMotorThermalSample<S = MaterialNetworkSample> {
    pub network: S,
    pub channels: Arc<[ConverterChannelMaterial; BATTERY_MOTOR_RESISTIVE_CHANNELS]>,
    pub bearing_cell: CellId,
    /// Explicit fixed native coefficient, with no inferred material authority.
    pub bearing_drag_nm_s_per_rad: f64,
}

/// All nine native coordinates and thermal history publish together.
pub struct BatteryMotorEnthalpyCandidate<'t, S = MaterialNetworkSample, Step = BatteryMotorStep> {
    thermal: DepositionSimulation<'t>,
    drive: MaterialBatteryMotorDrive,
    report: CircuitEnthalpyReport<BatteryMotorThermalSample<S>, Step>,
}
impl<'t, S, Step> BatteryMotorEnthalpyCandidate<'t, S, Step> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn drive(&self) -> &MaterialBatteryMotorDrive {
        &self.drive
    }
    pub fn report(&self) -> &CircuitEnthalpyReport<BatteryMotorThermalSample<S>, Step> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        MaterialBatteryMotorDrive,
        CircuitEnthalpyReport<BatteryMotorThermalSample<S>, Step>,
    ) {
        (self.thermal, self.drive, self.report)
    }
}

fn drive_error<E: From<ConverterThermalError>>(error: Error) -> CircuitEnthalpyError<E> {
    error.map_source(E::from)
}

fn drive_resistances(
    p: BuckBoostParameters,
    b: BatteryParameters,
    lead: f64,
    m: DcMotorParameters,
) -> [f64; BATTERY_MOTOR_RESISTIVE_CHANNELS] {
    let electrical = all_resistances(p, b, lead);
    core::array::from_fn(|i| {
        if i < BATTERY_FED_HEAT_CHANNELS {
            electrical[i]
        } else {
            m.winding_resistance_ohm
        }
    })
}

fn drive_parameters_at(
    mut p: BuckBoostParameters,
    mut b: BatteryParameters,
    mut m: DcMotorParameters,
    channels: &[ConverterChannelMaterial; BATTERY_MOTOR_RESISTIVE_CHANNELS],
) -> (
    BuckBoostParameters,
    BatteryParameters,
    f64,
    DcMotorParameters,
) {
    p.winding_resistance_ohm = channels[0].resistance_ohm;
    p.input_on_resistance_ohm = [channels[1].resistance_ohm, channels[2].resistance_ohm];
    p.output_on_resistance_ohm = [channels[3].resistance_ohm, channels[4].resistance_ohm];
    b.series_resistance_ohm = channels[5].resistance_ohm;
    b.polarization_resistance_ohm = channels[6].resistance_ohm;
    m.winding_resistance_ohm = channels[8].resistance_ohm;
    (p, b, channels[7].resistance_ohm, m)
}

fn drive_native(
    curve: &BatteryOcvCurve,
    p: BuckBoostParameters,
    b: BatteryParameters,
    modulation: BridgeModulation,
    lead: f64,
    m: DcMotorParameters,
) -> Result<BatteryMotorDrive, Error> {
    native(curve, p, b, modulation, lead)?
        .with_dc_motor(DcMotor::new(m).map_err(model)?)
        .map_err(model)
}

// Both original native owners expose the same physical storage/work/heat.
// This private view adds no trajectory, reconstructed loss or caloric solver.
struct DriveEnergy {
    storage: [f64; 7],
    heat: [f64; BATTERY_MOTOR_HEAT_CHANNELS],
    battery_work: f64,
    converter_work: [f64; 2],
    motor_work: f64,
    shaft_work: [f64; 2],
    load_work: f64,
    electrical_external: Option<f64>,
}
trait DriveStep: Clone {
    fn native(&self) -> &fs_phs::PreparedStepRecord;
    fn energy(&self) -> DriveEnergy;
}
impl DriveStep for BatteryMotorStep {
    fn native(&self) -> &fs_phs::PreparedStepRecord {
        &self.native
    }
    fn energy(&self) -> DriveEnergy {
        let e = &self.electrical;
        DriveEnergy {
            storage: [
                e.battery.chemical_energy_change_j,
                e.battery.polarization_energy_change_j,
                e.battery.bus_energy_change_j,
                e.converter_storage_change_j,
                self.magnetic_energy_change_j,
                self.kinetic_energy_change_j[0],
                self.kinetic_energy_change_j[1],
            ],
            heat: self.heat_joules(),
            battery_work: e.battery.terminal_work_j,
            converter_work: e.converter_terminal_work_j,
            motor_work: self.motor_electrical_work_j,
            shaft_work: self.shaft_work_j,
            load_work: e.load_work_j,
            electrical_external: Some(e.external_work_j),
        }
    }
}
impl DriveStep for BatteryMotorSwitchCycleStep {
    fn native(&self) -> &fs_phs::PreparedStepRecord {
        &self.native
    }
    fn energy(&self) -> DriveEnergy {
        DriveEnergy {
            storage: self.storage_change_j,
            heat: self.heat_joules(),
            battery_work: self.battery_terminal_work_j,
            converter_work: self.converter_terminal_work_j,
            motor_work: self.motor_electrical_work_j,
            shaft_work: self.shaft_work_j,
            load_work: self.load_work_j,
            electrical_external: None,
        }
    }
}

struct DriveDriver<Step> {
    drive: MaterialBatteryMotorDrive,
    step: Step,
    heat: [f64; BATTERY_MOTOR_HEAT_CHANNELS],
    total_heat_j: f64,
    net_external_work_j: Dd,
}
impl<Step: DriveStep> ThermalDriver for DriveDriver<Step> {
    type Model = MaterialBatteryMotorDrive;
    type Step = Step;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.total_heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: self.step.native().delta_h,
            supplied_j: self.net_external_work_j.to_f64(),
        }
    }
    fn retained_energy_balance_j(&self) -> Dd {
        Dd::from_f64(self.step.native().delta_h) - self.net_external_work_j
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.drive, self.step)
    }
}

// Audit each physical boundary BEFORE depositing heat. The original averaged
// owner's additional six-state external-port witness remains independently checked.
fn drive_account<Step: DriveStep>(
    step: &Step,
    tolerance_j: f64,
) -> Result<([f64; BATTERY_MOTOR_HEAT_CHANNELS], Dd, Dd), Error> {
    let e = step.energy();
    let heat = e.heat;
    if heat.iter().any(|v| !v.is_finite() || *v < 0.0) {
        return Err(invalid(
            "battery motor drive needs ten nonnegative physical heats",
        ));
    }
    let total = sum(heat)?;
    let exported = Dd::from_f64(finite(e.load_work)?);
    let shaft_work = sum(e.shaft_work)?;
    let external = shaft_work - exported;
    let battery_change = sum(e.storage[..3].iter().copied())?;
    let battery_heat = sum([heat[5], heat[6]])?;
    let converter_change = Dd::from_f64(finite(e.storage[3])?);
    let converter_heat = sum(heat[..5].iter().copied())?;
    let motor_change = sum(e.storage[4..].iter().copied())?;
    let motor_heat = sum([heat[8], heat[9]])?;
    let battery_work = Dd::from_f64(finite(e.battery_work)?);
    let input_work = Dd::from_f64(finite(e.converter_work[0])?);
    let output_work = Dd::from_f64(finite(e.converter_work[1])?);
    let motor_work = Dd::from_f64(finite(e.motor_work)?);
    let delta = Dd::from_f64(finite(step.native().delta_h)?);
    let valid = |r: Dd| r.hi.is_finite() && r.lo.is_finite() && r.abs().to_f64() <= tolerance_j;
    let residuals = [
        battery_change + battery_heat - battery_work,
        converter_change + converter_heat - input_work - output_work,
        battery_work + Dd::from_f64(heat[7]) + input_work,
        output_work + motor_work + exported,
        motor_change + motor_heat - motor_work - shaft_work,
        battery_change + converter_change + motor_change - delta,
        total - Dd::from_f64(finite(step.native().dissipated)?),
        Dd::from_f64(finite(step.native().supplied)?) - external,
        delta + total - external,
    ];
    let original_port_closes = if let Some(work) = e.electrical_external {
        let work = Dd::from_f64(finite(work)?);
        valid(output_work - work + exported) && valid(work + motor_work)
    } else {
        true
    };
    if exported.to_f64() < 0.0 || !original_port_closes || !residuals.into_iter().all(valid) {
        return Err(invalid(
            "battery, converter, motor, shaft and exported load fail independent energy balances",
        ));
    }
    Ok((heat, total, external))
}

// One original averaged native solve; shared by the ordinary and HeatAssembly
// entry points. The returned modulation is accepted command metadata only.
#[allow(clippy::too_many_arguments)]
fn averaged_drive_step(
    original: BatteryMotorDrive,
    state: &[f64; 9],
    dt: f64,
    torques: &[f64; 2],
    load: Option<RegulatedPowerLoad>,
    controls: PowerLoadControls,
    checkpoint: &mut dyn FnMut() -> bool,
) -> Result<([f64; 9], BatteryMotorStep, BridgeModulation), Error> {
    if torques.iter().any(|v| !v.is_finite()) {
        return Err(invalid(
            "battery motor drive needs finite authored shaft torques",
        ));
    }
    let modulation = original.electrical().modulation();
    let mut native = original
        .prepare_controlled(dt, controls, &mut *checkpoint)
        .map_err(electrical)?;
    let mut next = [0.0; 9];
    let step = native
        .step_controlled(state, torques, load, &mut next, &mut *checkpoint)
        .map_err(electrical)?;
    Ok((next, step, modulation))
}

#[allow(clippy::too_many_arguments)]
fn switched_drive_step(
    original: BatteryMotorDrive,
    state: &[f64; 9],
    dt: f64,
    start: f64,
    target: f64,
    intervals: &[SwitchInterval],
    commands: &[BatteryMotorSwitchCommand],
    subdivisions: usize,
    controls: PowerLoadControls,
    cycle: SwitchCycleControls,
    checkpoint: &mut dyn FnMut() -> bool,
) -> Result<([f64; 9], BatteryMotorSwitchCycleStep, BridgeModulation), Error> {
    let mut native = PreparedBatteryMotorSwitchCycle::new_controlled(
        original,
        intervals,
        subdivisions,
        controls,
        cycle,
        &mut *checkpoint,
    )
    .map_err(electrical)?;
    let period = native.period_s();
    let end = finite(start + period)?;
    let clock_roundoff = 0.5 * (clock_spacing(start) + clock_spacing(target));
    let duration_roundoff = 64.0 * f64::EPSILON * period.max(dt);
    if end <= start
        || end.to_bits() != target.to_bits()
        || (period - dt).abs() > clock_roundoff + duration_roundoff
    {
        return Err(invalid(
            "the complete battery motor switch cycle must span the thermal interval",
        ));
    }
    let mut next = [0.0; 9];
    let step = native
        .step_interval_commands_controlled(state, commands, &mut next, &mut *checkpoint)
        .map_err(electrical)?;
    Ok((next, step, native.averaged_modulation()))
}

fn switched_drive_work(
    intervals: &[SwitchInterval],
    commands: &[BatteryMotorSwitchCommand],
    subdivisions: usize,
) -> Result<usize, Error> {
    if commands.len() != intervals.len() {
        return Err(invalid(
            "each battery motor switch interval requires its own shaft/load command",
        ));
    }
    intervals
        .len()
        .checked_mul(subdivisions)
        .filter(|work| *work != 0)
        .ok_or_else(|| invalid("battery motor switching requires representable nonzero work"))
}

impl BatteryMotorThermalNetwork<'_> {
    fn network(&self) -> MaterialNetwork<'_> {
        MaterialNetwork {
            cells: self.cells,
            links: self.links,
            resistors: &[],
        }
    }

    fn admit_bearing(&self, thermal: &DepositionSimulation<'_>) -> Result<(), Error> {
        if !self.cells.iter().any(|c| c.cell == self.bearing_cell)
            || !thermal.cells().contains_key(&self.bearing_cell)
        {
            return Err(invalid(
                "motor bearing heat needs an active original material receiver",
            ));
        }
        Ok(())
    }

    /// Bind actual material resistance without resetting any charge, flux or
    /// momentum. The explicit bearing receiver may share an authored body.
    pub fn bind(
        &self,
        thermal: &DepositionSimulation<'_>,
        native: &BatteryMotorDrive,
        state: [f64; 9],
        cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialBatteryMotorDrive, Error> {
        let cancel = RefCell::new(cancelled);
        self.bind_resolved(
            thermal,
            native,
            state,
            |fixed| {
                let network = self.network();
                let prepared = network
                    .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
                    .map_err(material)?;
                let sample = prepared
                    .resolve(thermal, || (*cancel.borrow_mut())())
                    .map_err(material)?;
                resolve_channels(
                    self.cells,
                    self.channels,
                    thermal,
                    &sample.source,
                    fixed,
                    None,
                    false,
                    &cancel,
                )
            },
            &cancel,
        )
    }

    fn bind_resolved<E: From<ConverterThermalError>>(
        &self,
        thermal: &DepositionSimulation<'_>,
        native: &BatteryMotorDrive,
        state: [f64; 9],
        resolve: impl FnOnce(
            [f64; BATTERY_MOTOR_RESISTIVE_CHANNELS],
        ) -> Result<
            [ConverterChannelMaterial; BATTERY_MOTOR_RESISTIVE_CHANNELS],
            CircuitEnthalpyError<E>,
        >,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<MaterialBatteryMotorDrive, CircuitEnthalpyError<E>> {
        check(cancel).map_err(drive_error)?;
        self.admit_bearing(thermal).map_err(drive_error)?;
        if state.iter().any(|v| !v.is_finite()) || !native.system().hamiltonian(&state).is_finite()
        {
            return Err(drive_error(invalid(
                "battery motor drive needs nine finite original storage coordinates",
            )));
        }
        let electrical = native.electrical();
        let channels = Arc::new(resolve(drive_resistances(
            electrical.parameters(),
            electrical.battery_parameters(),
            electrical.lead_resistance_ohm(),
            native.motor_parameters(),
        ))?);
        let (converter_parameters, battery_parameters, lead, motor_parameters) =
            drive_parameters_at(
                electrical.parameters(),
                electrical.battery_parameters(),
                native.motor_parameters(),
                &channels,
            );
        drive_native(
            electrical.battery_curve(),
            converter_parameters,
            battery_parameters,
            electrical.modulation(),
            lead,
            motor_parameters,
        )
        .map_err(drive_error)?;
        check(cancel).map_err(drive_error)?;
        Ok(MaterialBatteryMotorDrive {
            time_s: thermal.time_s(),
            curve: electrical.battery_curve().clone(),
            battery_parameters,
            converter_parameters,
            modulation: electrical.modulation(),
            lead_resistance_ohm: lead,
            motor_parameters,
            state,
            channels,
            bearing_cell: self.bearing_cell,
        })
    }

    /// Solve finite source, conversion, motor and all heat receivers together.
    /// Both shaft torques enter their absolute-speed ports; positive auxiliary
    /// demand exports work. Neither of those exchanges is bearing/copper heat.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        drive: &MaterialBatteryMotorDrive,
        modulation: BridgeModulation,
        shaft_torque_nm: &[f64; 2],
        auxiliary_load: Option<RegulatedPowerLoad>,
        electrical_controls: PowerLoadControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<BatteryMotorEnthalpyCandidate<'t>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        super::super::super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        let network = self.network();
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        self.step_resolved_candidate(
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
            |trial, fixed, original| {
                let coefficients = prepared
                    .resolve(trial, || (*cancel.borrow_mut())())
                    .map_err(material)?;
                let channels = resolve_channels(
                    self.cells,
                    self.channels,
                    trial,
                    &coefficients.source,
                    fixed,
                    original,
                    false,
                    &cancel,
                )?;
                Ok((coefficients, channels))
            },
            &cancel,
        )
    }

    /// Resolve every authored PWM interval, all nine R(T,phase) sources and all
    /// ten heat receivers in the SAME finite-cell solve. Auxiliary load and
    /// shaft work remain external. Every numerical subdivision consumes one
    /// circuit trial, including exact-coefficient cache reuse. The source/gate
    /// program stays fixed within this transaction; storage laws stay fixed.
    #[allow(clippy::too_many_arguments)]
    pub fn step_switch_candidate<'t>(
        &self,
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
        cancelled: impl FnMut() -> bool,
    ) -> Result<
        BatteryMotorEnthalpyCandidate<'t, MaterialNetworkSample, BatteryMotorSwitchCycleStep>,
        Error,
    > {
        let integration_steps = switched_drive_work(intervals, commands, substeps_per_interval)?;
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        super::super::super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        let network = self.network();
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        self.step_resolved_candidate(
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
            |trial, fixed, original| {
                let coefficients = prepared
                    .resolve(trial, || (*cancel.borrow_mut())())
                    .map_err(material)?;
                let channels = resolve_channels(
                    self.cells,
                    self.channels,
                    trial,
                    &coefficients.source,
                    fixed,
                    original,
                    false,
                    &cancel,
                )?;
                Ok((coefficients, channels))
            },
            &cancel,
        )
    }

    // Only source/transport composition varies. Native integration, heat
    // deposition and acceptance remain identical for every heat provider.
    #[allow(clippy::too_many_arguments)]
    fn step_resolved_candidate<'t, S, E: From<ConverterThermalError>, Step: DriveStep>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        drive: &MaterialBatteryMotorDrive,
        modulation: BridgeModulation,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        integration_steps: usize,
        advance: impl Fn(
            BatteryMotorDrive,
            f64,
            &mut dyn FnMut() -> bool,
        ) -> Result<([f64; 9], Step, BridgeModulation), Error>,
        resolve: impl Fn(
            &DepositionSimulation<'t>,
            [f64; BATTERY_MOTOR_RESISTIVE_CHANNELS],
            Option<&[ConverterChannelMaterial; BATTERY_MOTOR_RESISTIVE_CHANNELS]>,
        ) -> Result<
            (
                CircuitThermalCoefficients<S>,
                [ConverterChannelMaterial; BATTERY_MOTOR_RESISTIVE_CHANNELS],
            ),
            CircuitEnthalpyError<E>,
        >,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<BatteryMotorEnthalpyCandidate<'t, S, Step>, CircuitEnthalpyError<E>> {
        check(cancel).map_err(drive_error)?;
        super::super::super::admit_driver_unknowns(thermal, controls, initial_energy_increments_j)?;
        self.admit_bearing(thermal).map_err(drive_error)?;
        let dt = finite(target_s - thermal.time_s())?;
        if drive.time_s.to_bits() != thermal.time_s().to_bits()
            || self.bearing_cell != drive.bearing_cell
            || dt <= 0.0
            || integration_steps == 0
        {
            return Err(drive_error(invalid(
                "battery motor and heat need the same positive interval, receivers and nonzero native work",
            )));
        }
        let fixed = drive_resistances(
            drive.converter_parameters,
            drive.battery_parameters,
            drive.lead_resistance_ohm,
            drive.motor_parameters,
        );
        let (_, initial) = resolve(thermal, fixed, None)?;
        if initial
            .iter()
            .zip(drive.channels.iter())
            .any(|(a, b)| !same_channel(a, b))
        {
            return Err(drive_error(invalid(
                "battery motor must retain all nine accepted resistance sources",
            )));
        }
        let cache = RefCell::new(
            None::<(
                [u64; BATTERY_MOTOR_RESISTIVE_CHANNELS],
                [f64; 9],
                Step,
                BridgeModulation,
            )>,
        );
        let (thermal, drive, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |trial, remaining| {
                check(cancel).map_err(drive_error)?;
                if remaining < integration_steps {
                    return Err(CircuitEnthalpyError::CircuitBudget);
                }
                let (coefficients, channels) = resolve(trial, fixed, Some(&drive.channels))?;
                let channels = Arc::new(channels);
                let (converter_parameters, battery_parameters, lead, motor_parameters) =
                    drive_parameters_at(
                        drive.converter_parameters,
                        drive.battery_parameters,
                        drive.motor_parameters,
                        &channels,
                    );
                let bits = drive_resistances(
                    converter_parameters,
                    battery_parameters,
                    lead,
                    motor_parameters,
                )
                .map(f64::to_bits);
                let previous = cache
                    .borrow()
                    .as_ref()
                    .filter(|entry| entry.0 == bits)
                    .cloned();
                let reused = previous.is_some();
                let (next, step, accepted_modulation) = if let Some((_, next, step, accepted)) =
                    previous
                {
                    (next, step, accepted)
                } else {
                    let original = drive_native(
                        &drive.curve,
                        converter_parameters,
                        battery_parameters,
                        modulation,
                        lead,
                        motor_parameters,
                    )
                    .map_err(drive_error)?;
                    advance(original, dt, &mut || (*cancel.borrow_mut())()).map_err(drive_error)?
                };
                let (heat, total, net_external_work_j) =
                    drive_account(&step, controls.total_energy_tolerance_j).map_err(drive_error)?;
                check(cancel).map_err(drive_error)?;
                if !reused {
                    *cache.borrow_mut() = Some((bits, next, step.clone(), accepted_modulation));
                }
                Ok(DrivenCoefficients {
                    circuit: DriveDriver {
                        drive: MaterialBatteryMotorDrive {
                            time_s: target_s,
                            curve: drive.curve.clone(),
                            battery_parameters,
                            converter_parameters,
                            modulation: accepted_modulation,
                            lead_resistance_ohm: lead,
                            motor_parameters,
                            state: next,
                            channels: Arc::clone(&channels),
                            bearing_cell: self.bearing_cell,
                        },
                        step,
                        heat,
                        total_heat_j: total.to_f64(),
                        net_external_work_j,
                    },
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: BatteryMotorThermalSample {
                        network: coefficients.source,
                        channels,
                        bearing_cell: self.bearing_cell,
                        bearing_drag_nm_s_per_rad: motor_parameters.bearing_drag_nm_s_per_rad,
                    },
                    circuit_trials: integration_steps,
                })
            },
            |candidate, duration, checkpoint| {
                let mut by_cell: BTreeMap<CellId, Dd> = BTreeMap::new();
                for (cell, heat) in candidate
                    .drive
                    .channels
                    .iter()
                    .map(|c| c.cell)
                    .chain([candidate.drive.bearing_cell])
                    .zip(candidate.heat)
                {
                    checkpoint()?;
                    let total = by_cell.entry(cell).or_insert(Dd::ZERO);
                    *total = *total + Dd::from_f64(heat);
                }
                let mut powers = BTreeMap::new();
                for (cell, heat) in by_cell {
                    checkpoint()?;
                    let heat = finite(heat.to_f64())?;
                    let power = finite(heat / duration)?;
                    if heat > 0.0 && power == 0.0 {
                        return Err(drive_error(invalid(
                            "battery motor heat cannot be represented as a rate",
                        )));
                    }
                    powers.insert(cell, power);
                }
                Ok(powers)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(cancel).map_err(drive_error)?;
        Ok(BatteryMotorEnthalpyCandidate {
            thermal,
            drive,
            report,
        })
    }
}
