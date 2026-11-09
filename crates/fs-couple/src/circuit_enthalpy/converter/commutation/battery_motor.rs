//! Finite battery, device commutation, both motor shafts and material heat.
//!
//! Original battery/motor storage and the native nine-state commutation solve
//! own the electrical/mechanical trajectory. Five ordinary resistance sources
//! and four complete device bundles refresh at each actual thermal proposal.
//! Channel and body-diode heat remain distinct in the native report and enter
//! their authored package body once. No averaged drive receipt is manufactured.
//!
//! OCV, capacitances, inductances, inertias, motor reciprocity and bearing drag
//! stay fixed. Endpoint material feedback needs time refinement. Native sampled
//! device limits do not certify unresolved peaks or semiconductor safe area.
//! The default nine-state lane has no device capacitor states; `capacitance`
//! adds explicit constant output-charge storage. Neither lane models gate-drive
//! energy, reverse recovery, avalanche or a general nonlinear capacitance law.

use super::*;

pub mod capacitance;
use crate::circuit_enthalpy::CircuitThermalCoefficients;
use crate::circuit_enthalpy::contacts::heat::{HeatAssembly, HeatNetworkError, closed_balance};
use crate::circuit_enthalpy::converter::battery::motor::heat::BatteryMotorHeatError;
use fs_exec::Cx;
use fs_phs::circuit::battery::{BatteryOcvCurve, BatteryParameters, BatterySource};
use fs_phs::circuit::converter::battery_fed::motor::BatteryMotorDrive;
use fs_phs::circuit::converter::commutation::battery_motor::{
    CommutatedBatteryMotorCycleStep, PreparedCommutatedBatteryMotorCycle,
};
use fs_phs::circuit::converter::switching::battery_motor::BatteryMotorSwitchCommand;
use fs_phs::circuit::motor::{DcMotor, DcMotorParameters};

/// Converter winding, battery series/polarization, lead, motor winding.
pub const DRIVE_RESISTORS: usize = 5;

/// Original finite bodies and independently owned dissipative channels.
/// Several channels/devices may share a body only through these explicit maps.
pub struct CommutatedBatteryMotorThermalNetwork<'a> {
    pub cells: &'a [MaterialCell<'a>],
    pub links: &'a [ConductionLink],
    pub resistors: &'a [ConverterThermalChannel<'a>; DRIVE_RESISTORS],
    pub devices: &'a [DeviceConductionBinding<'a>; 4],
    pub bearing_cell: CellId,
}

/// Accepted original drive and complete source-bound coefficients.
/// The default nine coordinates retain the existing algebraic-node API; the
/// capacitor child binds eleven coordinates with four immutable charge laws.
#[derive(Clone, Debug)]
pub struct MaterialCommutatedBatteryMotorDrive<const N: usize = 9, C = ()> {
    time_s: f64,
    curve: BatteryOcvCurve,
    battery_parameters: BatteryParameters,
    parameters: CommutationParameters,
    lead_resistance_ohm: f64,
    motor_parameters: DcMotorParameters,
    state: [f64; N],
    capacitors: C,
    resistors: Arc<[ConverterChannelMaterial; DRIVE_RESISTORS]>,
    devices: Arc<[DeviceConductionMaterial; 4]>,
    bearing_cell: CellId,
}
impl<const N: usize, C> MaterialCommutatedBatteryMotorDrive<N, C> {
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn battery_curve(&self) -> &BatteryOcvCurve {
        &self.curve
    }
    pub const fn battery_parameters(&self) -> BatteryParameters {
        self.battery_parameters
    }
    pub const fn parameters(&self) -> CommutationParameters {
        self.parameters
    }
    pub const fn lead_resistance_ohm(&self) -> f64 {
        self.lead_resistance_ohm
    }
    pub const fn motor_parameters(&self) -> DcMotorParameters {
        self.motor_parameters
    }
    pub const fn state(&self) -> &[f64; N] {
        &self.state
    }
    pub fn resistors(&self) -> &[ConverterChannelMaterial; DRIVE_RESISTORS] {
        &self.resistors
    }
    pub fn devices(&self) -> &[DeviceConductionMaterial; 4] {
        &self.devices
    }
    pub const fn bearing_cell(&self) -> CellId {
        self.bearing_cell
    }
}

#[derive(Clone, Debug)]
pub struct CommutatedBatteryMotorThermalSample<S = MaterialNetworkSample> {
    pub network: S,
    pub resistors: Arc<[ConverterChannelMaterial; DRIVE_RESISTORS]>,
    pub devices: Arc<[DeviceConductionMaterial; 4]>,
    pub bearing_cell: CellId,
    /// Authored mechanical coefficient; never interpreted as an ohmic source.
    pub bearing_drag_nm_s_per_rad: f64,
    /// Only actual auxiliary work is routed; None leaves it external.
    pub load_heat: Option<RegulatedLoadHeat>,
}

/// Original electrical, mechanical, source and caloric histories publish together.
#[must_use]
pub struct CommutatedBatteryMotorEnthalpyCandidate<'t, S = MaterialNetworkSample> {
    thermal: DepositionSimulation<'t>,
    drive: MaterialCommutatedBatteryMotorDrive,
    report: CircuitEnthalpyReport<
        CommutatedBatteryMotorThermalSample<S>,
        CommutatedBatteryMotorCycleStep,
    >,
}
impl<'t, S> CommutatedBatteryMotorEnthalpyCandidate<'t, S> {
    pub fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub fn drive(&self) -> &MaterialCommutatedBatteryMotorDrive {
        &self.drive
    }
    pub fn report(
        &self,
    ) -> &CircuitEnthalpyReport<
        CommutatedBatteryMotorThermalSample<S>,
        CommutatedBatteryMotorCycleStep,
    > {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        MaterialCommutatedBatteryMotorDrive,
        CircuitEnthalpyReport<
            CommutatedBatteryMotorThermalSample<S>,
            CommutatedBatteryMotorCycleStep,
        >,
    ) {
        (self.thermal, self.drive, self.report)
    }
}

fn map_drive<E: From<ConverterThermalError>>(error: Error) -> CircuitEnthalpyError<E> {
    error.map_source(E::from)
}
type AssemblyError<E> = CircuitEnthalpyError<BatteryMotorHeatError<E>>;
fn heat_error<E: From<HeatNetworkError>>(error: HeatNetworkError) -> AssemblyError<E> {
    match error {
        HeatNetworkError::Cancelled => CircuitEnthalpyError::Cancelled,
        other => CircuitEnthalpyError::Source(BatteryMotorHeatError::Heat(other.into())),
    }
}

fn resistances(
    p: CommutationParameters,
    b: BatteryParameters,
    lead: f64,
    m: DcMotorParameters,
) -> [f64; DRIVE_RESISTORS] {
    [
        p.winding_resistance_ohm,
        b.series_resistance_ohm,
        b.polarization_resistance_ohm,
        lead,
        m.winding_resistance_ohm,
    ]
}
fn parameters_at(
    mut p: CommutationParameters,
    mut b: BatteryParameters,
    mut m: DcMotorParameters,
    resistors: &[ConverterChannelMaterial; DRIVE_RESISTORS],
    devices: &[DeviceConductionMaterial; 4],
) -> Result<
    (
        CommutationParameters,
        BatteryParameters,
        f64,
        DcMotorParameters,
    ),
    Error,
> {
    p.winding_resistance_ohm = resistors[0].resistance_ohm();
    p.devices = core::array::from_fn(|i| devices[i].parameters());
    b.series_resistance_ohm = resistors[1].resistance_ohm();
    b.polarization_resistance_ohm = resistors[2].resistance_ohm();
    m.winding_resistance_ohm = resistors[4].resistance_ohm();
    p.validate().map_err(model)?;
    Ok((p, b, resistors[3].resistance_ohm(), m))
}

// Reassemble only ORIGINAL component owners. The native commutation constructor
// replaces the assembly's modulation/channel laws; none of their work is reused.
fn original_drive(
    curve: &BatteryOcvCurve,
    p: CommutationParameters,
    b: BatteryParameters,
    lead: f64,
    m: DcMotorParameters,
) -> Result<BatteryMotorDrive, Error> {
    AveragedBuckBoost::new(
        BuckBoostParameters {
            inductance_h: p.inductance_h,
            bus_capacitance_f: p.bus_capacitance_f,
            winding_resistance_ohm: p.winding_resistance_ohm,
            input_on_resistance_ohm: [
                p.devices[0].on_resistance_ohm,
                p.devices[1].on_resistance_ohm,
            ],
            output_on_resistance_ohm: [
                p.devices[2].on_resistance_ohm,
                p.devices[3].on_resistance_ohm,
            ],
        },
        BridgeModulation {
            input_high_fraction: 0.0,
            output_high_fraction: 0.0,
        },
    )
    .map_err(model)?
    .with_battery_source(BatterySource::new(curve.clone(), b).map_err(model)?, lead)
    .map_err(model)?
    .with_dc_motor(DcMotor::new(m).map_err(model)?)
    .map_err(model)
}

fn coefficient_bits(
    p: CommutationParameters,
    b: BatteryParameters,
    lead: f64,
    m: DcMotorParameters,
) -> [u64; 33] {
    // Reuse the existing complete device-key owner, with lead/motor resistance
    // in its two additional resistance slots; append both battery resistances.
    let prefix = super::coefficient_bits(p, [lead, m.winding_resistance_ohm]);
    let mut bits = [0; 33];
    bits[..31].copy_from_slice(&prefix);
    bits[31] = b.series_resistance_ohm.to_bits();
    bits[32] = b.polarization_resistance_ohm.to_bits();
    bits
}

fn admit_original_state(native: &BatteryMotorDrive, state: &[f64; 9]) -> Result<(), Error> {
    if !native.system().hamiltonian(state).is_finite() {
        return Err(invalid(
            "commutated battery drive needs nine supported original coordinates",
        ));
    }
    Ok(())
}

/// A view of actual native physical transfers, never an averaged receipt.
struct DrivePower<'a> {
    native: &'a fs_phs::PreparedStepRecord,
    storage_change_j: [f64; 7],
    heat: [f64; 14],
    battery_terminal_work_j: f64,
    converter_terminal_work_j: [f64; 2],
    motor_electrical_work_j: f64,
    electromagnetic_work_j: f64,
    shaft_work_j: [f64; 2],
    load_work_j: f64,
}
trait DriveReceipt: Clone {
    fn account(&self) -> DrivePower<'_>;
    fn admit_storage(&self, _allowance_j: f64) -> Result<(), Error> {
        Ok(())
    }
}
impl DriveReceipt for CommutatedBatteryMotorCycleStep {
    fn account(&self) -> DrivePower<'_> {
        DrivePower {
            native: &self.native,
            storage_change_j: self.storage_change_j,
            heat: self.heat_joules(),
            battery_terminal_work_j: self.battery_terminal_work_j,
            converter_terminal_work_j: self.converter_terminal_work_j,
            motor_electrical_work_j: self.motor_electrical_work_j,
            electromagnetic_work_j: self.electromagnetic_work_j,
            shaft_work_j: self.shaft_work_j,
            load_work_j: self.load_work_j,
        }
    }
}

fn admit_cycle_clock(
    start_s: f64,
    target_s: f64,
    period_s: f64,
    integration_steps: usize,
    expected_steps: usize,
) -> Result<(), Error> {
    let dt = finite(target_s - start_s)?;
    let represented_end = finite(start_s + period_s)?;
    let allowance = 0.5 * (clock_spacing(start_s) + clock_spacing(target_s))
        + 64.0 * f64::EPSILON * period_s.max(dt);
    if represented_end <= start_s
        || represented_end.to_bits() != target_s.to_bits()
        || (period_s - dt).abs() > allowance
        || integration_steps != expected_steps
    {
        return Err(invalid(
            "the commutated drive cycle must reach the same represented thermal endpoint",
        ));
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn advance_commutated(
    original: BatteryMotorDrive,
    devices: [SwitchDevice; 4],
    state: &[f64; 9],
    start_s: f64,
    target_s: f64,
    intervals: &[CommutationInterval],
    commands: &[BatteryMotorSwitchCommand],
    substeps: usize,
    electrical_controls: PowerLoadControls,
    cycle_controls: SwitchCycleControls,
    cancel: &RefCell<impl FnMut() -> bool>,
) -> Result<([f64; 9], CommutatedBatteryMotorCycleStep), Error> {
    let mut native = PreparedCommutatedBatteryMotorCycle::new_controlled(
        original,
        devices,
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
    let mut next = [0.0; 9];
    let step = native
        .step_interval_commands_controlled(state, commands, &mut next, || (*cancel.borrow_mut())())
        .map_err(electrical)?;
    Ok((next, step))
}

// Independent physical boundaries precede optional auxiliary-work routing.
// No heat is inferred from an energy remainder.
fn physical_account(receipt: &impl DriveReceipt, allowance_j: f64) -> Result<(Dd, Dd), Error> {
    receipt.admit_storage(allowance_j)?;
    let step = receipt.account();
    let heat = step.heat;
    if heat
        .iter()
        .chain([&step.load_work_j])
        .any(|v| !v.is_finite() || *v < 0.0)
    {
        return Err(invalid(
            "commutated drive needs fourteen nonnegative physical heats and actual load work",
        ));
    }
    let total = sum(heat)?;
    let s = step.storage_change_j;
    let battery = sum(s[..3].iter().copied())?;
    let converter = Dd::from_f64(finite(s[3])?);
    let magnetic = Dd::from_f64(finite(s[4])?);
    let kinetic = sum([s[5], s[6]])?;
    let battery_work = Dd::from_f64(finite(step.battery_terminal_work_j)?);
    let input = Dd::from_f64(finite(step.converter_terminal_work_j[0])?);
    let output = Dd::from_f64(finite(step.converter_terminal_work_j[1])?);
    let motor = Dd::from_f64(finite(step.motor_electrical_work_j)?);
    let electromagnetic = Dd::from_f64(finite(step.electromagnetic_work_j)?);
    let shaft = sum(step.shaft_work_j)?;
    let load = Dd::from_f64(step.load_work_j);
    let external = shaft - load;
    let delta = Dd::from_f64(finite(step.native.delta_h)?);
    let converter_heat = sum(heat[..5].iter().chain(&heat[10..]).copied())?;
    let residuals = [
        battery + sum([heat[5], heat[6]])? - battery_work,
        battery_work + Dd::from_f64(heat[7]) + input,
        converter + converter_heat - input - output,
        output + motor + load,
        magnetic + Dd::from_f64(heat[8]) + electromagnetic - motor,
        kinetic + Dd::from_f64(heat[9]) - electromagnetic - shaft,
        battery + converter + magnetic + kinetic - delta,
        total - Dd::from_f64(finite(step.native.dissipated)?),
        Dd::from_f64(finite(step.native.supplied)?) - external,
        delta + total - external,
    ];
    if residuals
        .into_iter()
        .any(|r| !r.hi.is_finite() || !r.lo.is_finite() || r.abs().to_f64() > allowance_j)
    {
        return Err(invalid(
            "finite battery, commutation devices, motor and shafts must close independently",
        ));
    }
    Ok((total, external))
}

struct DriveDriver<const N: usize, C, Step> {
    drive: MaterialCommutatedBatteryMotorDrive<N, C>,
    step: Step,
    load_heat: Option<RegulatedLoadHeat>,
    total_heat_j: f64,
    external_work_j: Dd,
}
impl<const N: usize, C, Step: DriveReceipt> ThermalDriver for DriveDriver<N, C, Step> {
    type Model = MaterialCommutatedBatteryMotorDrive<N, C>;
    type Step = Step;
    fn account(&self) -> DriverAccount<'_> {
        DriverAccount {
            resistors: &[],
            contacts: &[],
            heat_j: self.total_heat_j,
            numerical_dissipation_j: 0.0,
            stored_change_j: self.step.account().native.delta_h,
            supplied_j: self.external_work_j.to_f64(),
        }
    }
    fn retained_energy_balance_j(&self) -> Dd {
        Dd::from_f64(self.step.account().native.delta_h) - self.external_work_j
    }
    fn into_parts(self) -> (Self::Model, Self::Step) {
        (self.drive, self.step)
    }
}

impl CommutatedBatteryMotorThermalNetwork<'_> {
    fn network(&self) -> MaterialNetwork<'_> {
        MaterialNetwork {
            cells: self.cells,
            links: self.links,
            resistors: &[],
        }
    }
    fn admit_receiver(
        &self,
        thermal: &DepositionSimulation<'_>,
        receiver: CellId,
    ) -> Result<(), Error> {
        if !self.cells.iter().any(|c| c.cell == receiver)
            || !thermal.cells().contains_key(&receiver)
        {
            return Err(invalid(
                "commutated drive heat requires an active original material receiver",
            ));
        }
        Ok(())
    }
    fn resolve_bindings<const N: usize, C>(
        &self,
        thermal: &DepositionSimulation<'_>,
        sample: &MaterialNetworkSample,
        fixed: [f64; DRIVE_RESISTORS],
        retained: Option<&MaterialCommutatedBatteryMotorDrive<N, C>>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<
        (
            [ConverterChannelMaterial; DRIVE_RESISTORS],
            [DeviceConductionMaterial; 4],
        ),
        Error,
    > {
        let resistors = resolve_channels(
            self.cells,
            self.resistors,
            thermal,
            sample,
            fixed,
            retained.map(|old| old.resistors.as_ref()),
            false,
            cancel,
        )?;
        let devices = resolve_devices(
            self.cells,
            self.devices,
            thermal,
            sample,
            retained.map(|old| old.devices.as_ref()),
            cancel,
        )?;
        Ok((resistors, devices))
    }

    /// Preserve original OCV, L/C, motor/inertias and all nine coordinates.
    /// Actual complete device bundles replace the original averaged channels;
    /// no gate-dependent operating support is claimed until a program steps.
    pub fn bind(
        &self,
        thermal: &DepositionSimulation<'_>,
        native: &BatteryMotorDrive,
        state: [f64; 9],
        cancelled: impl FnMut() -> bool,
    ) -> Result<MaterialCommutatedBatteryMotorDrive, Error> {
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
                self.resolve_bindings::<9, ()>(thermal, &sample.source, fixed, None, &cancel)
                    .map(|(r, d)| (r, d, ()))
            },
            |state, _, _| admit_original_state(native, state),
            &cancel,
        )
    }

    fn bind_resolved<const N: usize, C, E: From<ConverterThermalError>>(
        &self,
        thermal: &DepositionSimulation<'_>,
        native: &BatteryMotorDrive,
        state: [f64; N],
        resolve: impl FnOnce(
            [f64; DRIVE_RESISTORS],
        ) -> Result<
            (
                [ConverterChannelMaterial; DRIVE_RESISTORS],
                [DeviceConductionMaterial; 4],
                C,
            ),
            CircuitEnthalpyError<E>,
        >,
        admit_state: impl FnOnce(&[f64; N], &C, CommutationParameters) -> Result<(), Error>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<MaterialCommutatedBatteryMotorDrive<N, C>, CircuitEnthalpyError<E>> {
        check(cancel).map_err(map_drive)?;
        self.admit_receiver(thermal, self.bearing_cell)
            .map_err(map_drive)?;
        if state.iter().any(|v| !v.is_finite()) {
            return Err(map_drive(invalid(
                "commutated battery drive needs finite original charge and momentum coordinates",
            )));
        }
        let p = native.electrical().parameters();
        let b = native.electrical().battery_parameters();
        let m = native.motor_parameters();
        let (resistors, devices, capacitors) = resolve([
            p.winding_resistance_ohm,
            b.series_resistance_ohm,
            b.polarization_resistance_ohm,
            native.electrical().lead_resistance_ohm(),
            m.winding_resistance_ohm,
        ])?;
        let p = CommutationParameters {
            inductance_h: p.inductance_h,
            bus_capacitance_f: p.bus_capacitance_f,
            winding_resistance_ohm: p.winding_resistance_ohm,
            devices: core::array::from_fn(|i| devices[i].parameters()),
        };
        let (parameters, battery_parameters, lead, motor_parameters) =
            parameters_at(p, b, m, &resistors, &devices).map_err(map_drive)?;
        original_drive(
            native.electrical().battery_curve(),
            parameters,
            battery_parameters,
            lead,
            motor_parameters,
        )
        .map_err(map_drive)?;
        admit_state(&state, &capacitors, parameters).map_err(map_drive)?;
        check(cancel).map_err(map_drive)?;
        Ok(MaterialCommutatedBatteryMotorDrive {
            time_s: thermal.time_s(),
            curve: native.electrical().battery_curve().clone(),
            battery_parameters,
            parameters,
            lead_resistance_ohm: lead,
            motor_parameters,
            state,
            capacitors,
            resistors: Arc::new(resistors),
            devices: Arc::new(devices),
            bearing_cell: self.bearing_cell,
        })
    }

    /// One joint finite-source/shaft/caloric proposal. Every material evaluation
    /// refreshes all five ordinary resistances and all seven roles of each device.
    /// Fourteen native heat channels remain distinct; an optional wholly
    /// dissipative auxiliary receiver absorbs actual delivered work only.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        drive: &MaterialCommutatedBatteryMotorDrive,
        intervals: &[CommutationInterval],
        commands: &[BatteryMotorSwitchCommand],
        substeps_per_interval: usize,
        load_heat_cell: Option<CellId>,
        electrical_controls: PowerLoadControls,
        cycle_controls: SwitchCycleControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<CommutatedBatteryMotorEnthalpyCandidate<'t>, Error> {
        let cancel = RefCell::new(cancelled);
        check(&cancel)?;
        crate::circuit_enthalpy::admit_driver_unknowns(
            thermal,
            controls,
            initial_energy_increments_j,
        )?;
        let network = self.network();
        let prepared = network
            .prepare_with_bases(thermal, &[], || (*cancel.borrow_mut())())
            .map_err(material)?;
        let (thermal, drive, report) = self.step_resolved_candidate(
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
                let (resistors, devices) =
                    self.resolve_bindings(trial, &coefficients.source, fixed, retained, &cancel)?;
                Ok((coefficients, resistors, devices, ()))
            },
            |original, devices, _, state| {
                advance_commutated(
                    original,
                    devices,
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
            |network, resistors, devices, (), load_heat, motor_parameters| {
                CommutatedBatteryMotorThermalSample {
                    network,
                    resistors,
                    devices,
                    bearing_cell: self.bearing_cell,
                    bearing_drag_nm_s_per_rad: motor_parameters.bearing_drag_nm_s_per_rad,
                    load_heat,
                }
            },
            &cancel,
        )?;
        Ok(CommutatedBatteryMotorEnthalpyCandidate {
            thermal,
            drive,
            report,
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn step_resolved_candidate<
        't,
        S,
        Out,
        E: From<ConverterThermalError>,
        const N: usize,
        C: Clone + PartialEq,
        Step: DriveReceipt,
    >(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        drive: &MaterialCommutatedBatteryMotorDrive<N, C>,
        intervals: &[CommutationInterval],
        commands: &[BatteryMotorSwitchCommand],
        substeps_per_interval: usize,
        load_heat_cell: Option<CellId>,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
        resolve: impl Fn(
            &DepositionSimulation<'t>,
            [f64; DRIVE_RESISTORS],
            Option<&MaterialCommutatedBatteryMotorDrive<N, C>>,
        ) -> Result<
            (
                CircuitThermalCoefficients<S>,
                [ConverterChannelMaterial; DRIVE_RESISTORS],
                [DeviceConductionMaterial; 4],
                C,
            ),
            CircuitEnthalpyError<E>,
        >,
        advance: impl Fn(
            BatteryMotorDrive,
            [SwitchDevice; 4],
            &C,
            &[f64; N],
        ) -> Result<([f64; N], Step), Error>,
        sample: impl Fn(
            S,
            Arc<[ConverterChannelMaterial; DRIVE_RESISTORS]>,
            Arc<[DeviceConductionMaterial; 4]>,
            C,
            Option<RegulatedLoadHeat>,
            DcMotorParameters,
        ) -> Out,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<
        (
            DepositionSimulation<'t>,
            MaterialCommutatedBatteryMotorDrive<N, C>,
            CircuitEnthalpyReport<Out, Step>,
        ),
        CircuitEnthalpyError<E>,
    > {
        check(cancel).map_err(map_drive)?;
        crate::circuit_enthalpy::admit_driver_unknowns(
            thermal,
            controls,
            initial_energy_increments_j,
        )?;
        let dt = finite(target_s - thermal.time_s())?;
        let integration_steps = intervals
            .len()
            .checked_mul(substeps_per_interval)
            .filter(|n| *n > 0)
            .ok_or_else(|| {
                map_drive(invalid(
                    "commutated drive requires bounded nonzero numerical work",
                ))
            })?;
        if drive.time_s.to_bits() != thermal.time_s().to_bits()
            || dt <= 0.0
            || self.bearing_cell != drive.bearing_cell
            || commands.len() != intervals.len()
        {
            return Err(map_drive(invalid(
                "commutated drive and heat need the same positive interval, bearing and gate commands",
            )));
        }
        self.admit_receiver(thermal, self.bearing_cell)
            .map_err(map_drive)?;
        if let Some(receiver) = load_heat_cell {
            self.admit_receiver(thermal, receiver).map_err(map_drive)?;
        }
        let fixed = resistances(
            drive.parameters,
            drive.battery_parameters,
            drive.lead_resistance_ohm,
            drive.motor_parameters,
        );
        let (_, original_resistors, original_devices, original_capacitors) =
            resolve(thermal, fixed, Some(drive))?;
        if original_resistors
            .iter()
            .zip(drive.resistors.iter())
            .any(|(a, b)| !same_channel(a, b))
            || original_devices != *drive.devices
            || original_capacitors != drive.capacitors
        {
            return Err(map_drive(invalid(
                "commutated drive must retain every accepted resistance, conduction and capacitance source",
            )));
        }
        // Every storage coefficient and the accepted state are fixed within this call.
        // Capacitor source/domain admission still runs on every thermal proposal.
        let cache = RefCell::new(None::<([u64; 33], [f64; N], Step)>);
        let (thermal, drive, report) = cell_driver_candidate(
            thermal,
            target_s,
            controls,
            initial_energy_increments_j,
            |trial, remaining| {
                check(cancel).map_err(map_drive)?;
                if remaining < integration_steps {
                    return Err(CircuitEnthalpyError::CircuitBudget);
                }
                let (coefficients, resistors, devices, capacitors) =
                    resolve(trial, fixed, Some(drive))?;
                let resistors = Arc::new(resistors);
                let devices = Arc::new(devices);
                let (parameters, battery_parameters, lead, motor_parameters) = parameters_at(
                    drive.parameters,
                    drive.battery_parameters,
                    drive.motor_parameters,
                    &resistors,
                    &devices,
                )
                .map_err(map_drive)?;
                let bits = coefficient_bits(parameters, battery_parameters, lead, motor_parameters);
                let cached = cache
                    .borrow()
                    .as_ref()
                    .filter(|entry| entry.0 == bits)
                    .cloned();
                let reused = cached.is_some();
                let (next, step) = if let Some((_, next, step)) = cached {
                    (next, step)
                } else {
                    let original = original_drive(
                        &drive.curve,
                        parameters,
                        battery_parameters,
                        lead,
                        motor_parameters,
                    )
                    .map_err(map_drive)?;
                    advance(original, parameters.devices, &capacitors, &drive.state)
                        .map_err(map_drive)?
                };
                let (heat, external) = physical_account(&step, controls.total_energy_tolerance_j)
                    .map_err(map_drive)?;
                let load_heat = load_heat_cell.map(|cell| RegulatedLoadHeat {
                    cell,
                    energy_j: step.account().load_work_j,
                });
                let routed = Dd::from_f64(load_heat.as_ref().map_or(0.0, |load| load.energy_j));
                let total_heat_j = finite((heat + routed).to_f64())?;
                let external_work_j = external + routed;
                finite(external_work_j.hi)?;
                finite(external_work_j.lo)?;
                check(cancel).map_err(map_drive)?;
                if !reused {
                    *cache.borrow_mut() = Some((bits, next, step.clone()));
                }
                Ok(DrivenCoefficients {
                    circuit: DriveDriver {
                        drive: MaterialCommutatedBatteryMotorDrive {
                            time_s: target_s,
                            curve: drive.curve.clone(),
                            battery_parameters,
                            parameters,
                            lead_resistance_ohm: lead,
                            motor_parameters,
                            state: next,
                            capacitors: capacitors.clone(),
                            resistors: Arc::clone(&resistors),
                            devices: Arc::clone(&devices),
                            bearing_cell: self.bearing_cell,
                        },
                        step,
                        load_heat: load_heat.clone(),
                        total_heat_j,
                        external_work_j,
                    },
                    sources_w: coefficients.sources_w,
                    conductances_w_k: coefficients.conductances_w_k,
                    source: sample(
                        coefficients.source,
                        resistors,
                        devices,
                        capacitors,
                        load_heat,
                        motor_parameters,
                    ),
                    circuit_trials: integration_steps,
                })
            },
            |driver, duration, checkpoint| {
                let heat = driver.step.account().heat;
                let mut by_cell: BTreeMap<CellId, Dd> = BTreeMap::new();
                let mut add = |cell, value| {
                    let sum = by_cell.entry(cell).or_insert(Dd::ZERO);
                    *sum = *sum + Dd::from_f64(value);
                };
                for (i, index) in [0, 5, 6, 7, 8].into_iter().enumerate() {
                    checkpoint()?;
                    add(driver.drive.resistors[i].cell(), heat[index]);
                }
                checkpoint()?;
                add(driver.drive.bearing_cell, heat[9]);
                for i in 0..4 {
                    checkpoint()?;
                    add(driver.drive.devices[i].cell(), heat[1 + i]);
                    add(driver.drive.devices[i].cell(), heat[10 + i]);
                }
                if let Some(load) = &driver.load_heat {
                    checkpoint()?;
                    add(load.cell, load.energy_j);
                }
                let mut powers = BTreeMap::new();
                for (cell, energy) in by_cell {
                    checkpoint()?;
                    let energy = finite(energy.to_f64())?;
                    let power = finite(energy / duration)?;
                    if energy > 0.0 && power == 0.0 {
                        return Err(map_drive(invalid(
                            "commutated drive heat is unrepresentable on this thermal interval",
                        )));
                    }
                    powers.insert(cell, power);
                }
                Ok(powers)
            },
            || (*cancel.borrow_mut())(),
        )?;
        check(cancel).map_err(map_drive)?;
        Ok((thermal, drive, report))
    }
}

/// Complete original heat assembly, including finite coolant/contact/radiation
/// transport, with the same five resistance, four device and bearing bindings.
/// Unowned electrical graph resistors/contacts refuse instead of losing work.
pub struct CommutatedBatteryMotorHeatNetwork<'a, A: HeatAssembly + ?Sized> {
    pub assembly: &'a A,
    pub resistors: &'a [ConverterThermalChannel<'a>; DRIVE_RESISTORS],
    pub devices: &'a [DeviceConductionBinding<'a>; 4],
    pub bearing_cell: CellId,
}
impl<A: HeatAssembly + ?Sized> CommutatedBatteryMotorHeatNetwork<'_, A> {
    fn network(
        &self,
        cx: &Cx<'_>,
    ) -> Result<CommutatedBatteryMotorThermalNetwork<'_>, AssemblyError<A::Error>> {
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        let prepared = self.assembly.prepared_heat_network();
        let bulk = prepared.thermal_material_network().map_err(heat_error)?;
        prepared
            .admit_field_independent_boundaries(cx, &[])
            .map_err(heat_error)?;
        Ok(CommutatedBatteryMotorThermalNetwork {
            cells: bulk.cells,
            links: bulk.links,
            resistors: self.resistors,
            devices: self.devices,
            bearing_cell: self.bearing_cell,
        })
    }
    fn resolve<const N: usize, C>(
        &self,
        cx: &Cx<'_>,
        network: &CommutatedBatteryMotorThermalNetwork<'_>,
        thermal: &DepositionSimulation<'_>,
        fixed: [f64; DRIVE_RESISTORS],
        retained: Option<&MaterialCommutatedBatteryMotorDrive<N, C>>,
        cancel: &RefCell<impl FnMut() -> bool>,
    ) -> Result<
        (
            CircuitThermalCoefficients<A::Sample>,
            [ConverterChannelMaterial; DRIVE_RESISTORS],
            [DeviceConductionMaterial; 4],
        ),
        AssemblyError<A::Error>,
    > {
        let coefficients = self
            .assembly
            .resolve_heat(cx, thermal)
            .map_err(|error| error.map_source(BatteryMotorHeatError::Heat))?;
        if !coefficients.contacts.is_empty() || !coefficients.bulk.resistors.is_empty() {
            return Err(heat_error(HeatNetworkError::Invalid(
                "commutated drive heat assembly cannot discard electrical contact or graph resistor work",
            )));
        }
        let (resistors, devices) = network
            .resolve_bindings(
                thermal,
                &coefficients.bulk.source.as_ref().material,
                fixed,
                retained,
                cancel,
            )
            .map_err(map_drive)?;
        Ok((coefficients.bulk, resistors, devices))
    }

    pub fn bind(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'_>,
        native: &BatteryMotorDrive,
        state: [f64; 9],
    ) -> Result<MaterialCommutatedBatteryMotorDrive, AssemblyError<A::Error>> {
        let network = self.network(cx)?;
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        network.bind_resolved(
            thermal,
            native,
            state,
            |fixed| {
                self.resolve::<9, ()>(cx, &network, thermal, fixed, None, &cancel)
                    .map(|(_, r, d)| (r, d, ()))
            },
            |state, _, _| admit_original_state(native, state),
            &cancel,
        )
    }

    /// The same source refresh, native trajectory cache, physical balances and
    /// deposition as the material-network entry point. Full heat-source samples
    /// and the original independent closed-radiation balance remain intact.
    /// Finite coolant transfer is internal; prescribed reservoir heat is external.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t>(
        &self,
        cx: &Cx<'_>,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        drive: &MaterialCommutatedBatteryMotorDrive,
        intervals: &[CommutationInterval],
        commands: &[BatteryMotorSwitchCommand],
        substeps_per_interval: usize,
        load_heat_cell: Option<CellId>,
        electrical_controls: PowerLoadControls,
        cycle_controls: SwitchCycleControls,
        controls: CircuitEnthalpySolve,
        initial_energy_increments_j: &[f64],
    ) -> Result<CommutatedBatteryMotorEnthalpyCandidate<'t, A::Sample>, AssemblyError<A::Error>>
    {
        cx.checkpoint()
            .map_err(|_| CircuitEnthalpyError::Cancelled)?;
        crate::circuit_enthalpy::admit_driver_unknowns(
            thermal,
            controls,
            initial_energy_increments_j,
        )?;
        let network = self.network(cx)?;
        let cancel = RefCell::new(|| cx.checkpoint().is_err());
        let duration_s = target_s - thermal.time_s();
        let (thermal, drive, report) = network.step_resolved_candidate(
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
                self.resolve(cx, &network, trial, fixed, retained, &cancel)
                    .map(|(coefficients, r, d)| (coefficients, r, d, ()))
            },
            |original, devices, _, state| {
                advance_commutated(
                    original,
                    devices,
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
            |network, resistors, devices, (), load_heat, motor_parameters| {
                CommutatedBatteryMotorThermalSample {
                    network,
                    resistors,
                    devices,
                    bearing_cell: self.bearing_cell,
                    bearing_drag_nm_s_per_rad: motor_parameters.bearing_drag_nm_s_per_rad,
                    load_heat,
                }
            },
            &cancel,
        )?;
        let candidate = CommutatedBatteryMotorEnthalpyCandidate {
            thermal,
            drive,
            report,
        };
        closed_balance(
            duration_s,
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
