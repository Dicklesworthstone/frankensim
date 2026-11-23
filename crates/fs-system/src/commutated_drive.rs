//! A battery-fed commutated motor drive whose resistances and switches heat
//! what their heat ports are joined to.
//!
//! [`CommutatedBatteryMotorDrive`] puts the drive of
//! `fs_phs::circuit::converter::commutation::capacitance::battery_motor`
//! behind `fs_ports::component::Component`. It joins the equivalent-circuit
//! battery and its lead to a two-leg converter and a DC motor whose rotor and
//! stator both turn. The converter has four switches, each a channel with
//! its body diode and with output charge storage. Nothing is integrated
//! here. A step is the owner's whole gate program
//! (`step_interval_commands_controlled`): every interval and substep in one
//! native workspace, on the owner built with the resistances and switches of
//! that step. A step whose length is not the program's period is refused.
//!
//! Ports, in this order:
//!
//! - `rotor`, `stator`: rotational, [`PortRole::EffortIn`]. Each receives the
//!   torque on its shaft from outside \[N m\], held through the program. It
//!   returns the shaft's mean angular velocity \[rad/s\], the owner's
//!   time-weighted mean of its discrete velocities, so torque times mean
//!   velocity times the step is the owner's shaft work.
//! - Ten heat ports, [`PortRole::EffortIn`], in [`HEAT_PORTS`] order:
//!   - the converter winding;
//!   - the four switches, a channel and its body diode together;
//!   - the battery's series and polarization resistances, the lead and the
//!     motor winding;
//!   - the bearing.
//!
//!   If the program has an auxiliary load, an eleventh,
//!   `heat_auxiliary_load`, carries the work the load took. Each heat port
//!   receives the temperature \[K\] where its heat leaves and returns minus
//!   that heat over the step divided by the step \[W\].
//!
//! State: the owner's eleven coordinates:
//! - the remaining chemical charge, polarization charge and battery bus
//!   charge \[C\];
//! - the converter's flux linkage \[Wb\], its two bus charges and its two
//!   switching-node charges \[C\];
//! - the motor's flux linkage \[Wb\] and its rotor and stator angular
//!   momenta \[N m s\].
//!
//! Its stored energy is the owner's Hamiltonian, the switches' output
//! storage included; no resistance or switch changes it.
//!
//! Switch storage ([`DeviceStorage`]) is either four constant output
//! capacitances or four voltage-charge curves whose energy is the integral
//! of voltage over charge. Either is fixed: storage that changed with
//! temperature would need its entropy, which this Component does not hold.
//!
//! Heat and work:
//! - Heat is the owner's fourteen accounts, routed to the heat ports. The
//!   auxiliary load's work leaves through its own port.
//! - Supplied work is the shaft work.
//! - The difference of heat and work from the stored change is booked as
//!   numerical loss. It is refused beyond the rounding of the step's
//!   arithmetic plus the declared energy tolerance, as is any difference
//!   between these accounts, the port works and the owner's native step
//!   record.
//!
//! Laws: a resistance or switch with no law keeps the value the drive was
//! built with. One following a [`ResistanceLaw`] or a [`SwitchLaw`] is given
//! the law's answer at the temperature its heat port receives in that step.
//! With a backward-Euler thermal timeline returning end-of-step
//! temperatures, that is the endpoint material resistance and conduction of
//! `fs_couple::circuit_enthalpy::converter::commutation::battery_motor`. A
//! law with no answer refuses the step ([`CommutatedDriveError::Law`]); it
//! is never clamped. Inductances, capacitances, the OCV, the motor's
//! coupling, its inertias and the bearing drag stay fixed.
//!
//! Class: [`Class::Smooth`], with no tangent, so the kernel differences the
//! port map. The map is smooth except where, as the inputs move, a body
//! diode starts or stops conducting, a switch reaches an operating limit, or
//! the auxiliary load changes branch. The module's own thermal Newton
//! differences the same map.
//!
//! Not claimed: torques that change within the program; gate charge and
//! gate-driver heat (the owner's fifteen-coordinate lane); reverse recovery
//! and avalanche; temperature-dependent storage, OCV or bearing drag.

use core::fmt;
use core::num::NonZeroUsize;
use std::sync::Arc;

use fs_exec::Cx;
use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_phs::circuit::battery::{BatteryOcvCurve, BatteryParameters, BatterySource};
use fs_phs::circuit::capacitance::CapacitorChargeCurve;
use fs_phs::circuit::converter::battery_fed::motor::BatteryMotorDrive;
use fs_phs::circuit::converter::commutation::capacitance::CapacitiveCommutationParameters;
use fs_phs::circuit::converter::commutation::capacitance::battery_motor::{
    CapacitiveBatteryMotorCycleStep, PreparedCapacitiveBatteryMotorCycle,
    PreparedChargeBatteryMotorCycle,
};
use fs_phs::circuit::converter::commutation::capacitance::charge::ChargeCommutationParameters;
use fs_phs::circuit::converter::commutation::{
    CommutationInterval, CommutationParameters, SwitchDevice,
};
use fs_phs::circuit::converter::power_load::{PowerLoadControls, RegulatedPowerLoad};
use fs_phs::circuit::converter::switching::SwitchCycleControls;
use fs_phs::circuit::converter::switching::battery_motor::BatteryMotorSwitchCommand;
use fs_phs::circuit::converter::{AveragedBuckBoost, BridgeModulation, BuckBoostParameters};
use fs_phs::circuit::motor::{DcMotor, DcMotorParameters};
use fs_phs::{PhsError, PreparedStepError, arithmetic_rounding_floor};
use fs_ports::component::{
    Class, Component, ContractViolation, PortLayout, PortRole, PortSpec, Proposal, Receipt,
};
use fs_ports::{CoupleError, PortKind, StableId};

use crate::circuit::ResistanceLaw;

/// The ten heat ports every drive has, in port order after `rotor` and
/// `stator`.
pub const HEAT_PORTS: [&str; 10] = [
    "heat_converter_winding",
    "heat_input_high",
    "heat_input_low",
    "heat_output_high",
    "heat_output_low",
    "heat_battery_series",
    "heat_battery_polarization",
    "heat_lead",
    "heat_motor_winding",
    "heat_bearing",
];

/// The port of the auxiliary load's work, after [`HEAT_PORTS`].
pub const LOAD_PORT: &str = "heat_auxiliary_load";

/// A switch as a function of the temperature of the package it heats.
pub trait SwitchLaw: fmt::Debug + Send + Sync {
    /// The switch's channel, body diode and operating limits at `kelvin`.
    ///
    /// # Errors
    /// Why the law has no answer there, for example a temperature outside
    /// the validity of the card behind it. The drive then refuses the step
    /// with [`CommutatedDriveError::Law`]; it never clamps the temperature.
    fn switch(&self, kelvin: f64) -> Result<SwitchDevice, String>;

    /// Words that identify the law; they join the drive's identity.
    fn identity_words(&self, words: &mut Vec<u64>);
}

/// One of the drive's five resistances.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DriveResistor {
    /// The converter's winding.
    ConverterWinding,
    /// The battery's series resistance.
    BatterySeries,
    /// The battery's polarization resistance.
    BatteryPolarization,
    /// The lead between the battery and the converter.
    Lead,
    /// The motor's winding.
    MotorWinding,
}

impl DriveResistor {
    const ALL: [Self; 5] = [
        Self::ConverterWinding,
        Self::BatterySeries,
        Self::BatteryPolarization,
        Self::Lead,
        Self::MotorWinding,
    ];

    const fn index(self) -> usize {
        self as usize
    }

    /// Its heat port's index in [`HEAT_PORTS`].
    const fn heat_port(self) -> usize {
        [0, 5, 6, 7, 8][self as usize]
    }
}

/// One of the converter's four switches, in the owner's order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SwitchPosition {
    /// The input leg's switch to the input bus.
    InputHigh,
    /// The input leg's switch to the return.
    InputLow,
    /// The output leg's switch to the output bus.
    OutputHigh,
    /// The output leg's switch to the return.
    OutputLow,
}

impl SwitchPosition {
    const fn index(self) -> usize {
        self as usize
    }
}

/// The switches' output storage, a fixed law in the owner's switch order.
#[derive(Debug, Clone)]
pub enum DeviceStorage {
    /// Four constant output capacitances \[F\].
    Capacitance([f64; 4]),
    /// Four charge-voltage curves. Each chord is traced through at most
    /// `maximum_cells_per_chord` cells of its curve.
    Charge {
        /// The curves, charge \[C\] against voltage \[V\].
        curves: [CapacitorChargeCurve; 4],
        /// The owner's bound on the cells one chord may cross.
        maximum_cells_per_chord: usize,
    },
}

/// The drive as authored: its battery, lead, converter and motor.
#[derive(Debug, Clone)]
pub struct DriveParameters {
    /// Usable capacity \[C\].
    pub capacity_c: f64,
    /// The OCV curve's knots: state of charge, volts.
    pub soc_voltage: Vec<(f64, f64)>,
    /// The battery's resistances and capacitances.
    pub battery: BatteryParameters,
    /// The lead's resistance \[ohm\].
    pub lead_resistance_ohm: f64,
    /// The converter's inductance, bus capacitances, winding and switches.
    pub converter: CommutationParameters,
    /// The motor.
    pub motor: DcMotorParameters,
}

/// The gate program one step runs, repeated every step.
#[derive(Debug, Clone)]
pub struct GateProgram {
    /// The gate states and their durations; their sum is the period.
    pub intervals: Vec<CommutationInterval>,
    /// The auxiliary load on the output bus in each interval, if any.
    pub auxiliary_loads: Vec<Option<RegulatedPowerLoad>>,
    /// The owner's substeps in each interval.
    pub substeps_per_interval: usize,
}

/// Refusal of the drive's owner or of its Component boundary.
#[derive(Debug)]
pub enum CommutatedDriveError {
    /// The owner refused the parameters, the storage, the program or a state.
    Model(PhsError),
    /// The owner refused the step, unchanged.
    Step(PreparedStepError),
    /// A law had no answer at the temperature its heat port received.
    Law {
        /// The heat port whose temperature the law was asked at.
        port: &'static str,
        /// The law's reason.
        reason: String,
    },
    /// Invalid or repeated port metadata.
    Port(ContractViolation),
    /// An invalid port identifier.
    Identifier(CoupleError),
    /// A shape, a value or a construction that is refused.
    Invalid {
        /// What was refused.
        what: &'static str,
    },
    /// Work, heat and storage do not close within the declared tolerance and
    /// rounding, or the owner's accounts of the step disagree.
    Balance {
        /// The largest defect \[J\].
        defect_j: f64,
        /// The allowance \[J\].
        tolerance_j: f64,
    },
    /// Cancellation was requested.
    Cancelled,
}

impl fmt::Display for CommutatedDriveError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Model(error) => write!(f, "drive owner refused: {error:?}"),
            Self::Step(error) => write!(f, "drive step refused: {error:?}"),
            Self::Law { port, reason } => write!(f, "the law at {port} has no answer: {reason}"),
            Self::Port(error) => write!(f, "drive port refused: {error:?}"),
            Self::Identifier(error) => write!(f, "drive port identifier refused: {error:?}"),
            Self::Invalid { what } => write!(f, "drive component refused: {what}"),
            Self::Balance {
                defect_j,
                tolerance_j,
            } => write!(
                f,
                "drive balance defect {defect_j} J exceeds {tolerance_j} J"
            ),
            Self::Cancelled => write!(f, "cancelled"),
        }
    }
}

impl std::error::Error for CommutatedDriveError {}

fn invalid(what: &'static str) -> CommutatedDriveError {
    CommutatedDriveError::Invalid { what }
}

fn step_error(error: PreparedStepError) -> CommutatedDriveError {
    match error {
        PreparedStepError::Cancelled => CommutatedDriveError::Cancelled,
        other @ PreparedStepError::Solver(_) => CommutatedDriveError::Step(other),
    }
}

/// The converter's storage law with the converter it belongs to.
#[derive(Debug, Clone)]
enum Storage {
    Capacitance(CapacitiveCommutationParameters),
    Charge(ChargeCommutationParameters),
}

impl Storage {
    fn energy(&self, converter: &[f64; 5]) -> Result<f64, PhsError> {
        match self {
            Self::Capacitance(p) => p.stored_energy_j(converter),
            Self::Charge(p) => p.stored_energy_j(converter),
        }
    }
}

/// A prepared native cycle of either storage law.
enum Prepared {
    Capacitance(PreparedCapacitiveBatteryMotorCycle),
    Charge(PreparedChargeBatteryMotorCycle),
}

impl Prepared {
    fn period_s(&self) -> f64 {
        match self {
            Self::Capacitance(native) => native.period_s(),
            Self::Charge(native) => native.period_s(),
        }
    }

    fn step(
        &mut self,
        old: &[f64; 11],
        commands: &[BatteryMotorSwitchCommand],
        next: &mut [f64; 11],
        cancelled: impl FnMut() -> bool,
    ) -> Result<CapacitiveBatteryMotorCycleStep, PreparedStepError> {
        match self {
            Self::Capacitance(native) => {
                native.step_interval_commands_controlled(old, commands, next, cancelled)
            }
            Self::Charge(native) => {
                native.step_interval_commands_controlled(old, commands, next, cancelled)
            }
        }
    }
}

/// The battery-fed commutated motor drive of fs-phs as one Component.
///
/// See the module documentation for the ports, the state, the heat and the
/// laws.
pub struct CommutatedBatteryMotorDrive {
    parameters: DriveParameters,
    curve: BatteryOcvCurve,
    storage: Storage,
    program: GateProgram,
    period_s: f64,
    electrical: PowerLoadControls,
    cycle: SwitchCycleControls,
    resistor_laws: [Option<Arc<dyn ResistanceLaw>>; 5],
    switch_laws: [Option<Arc<dyn SwitchLaw>>; 4],
    /// The battery and the motor at the authored resistances, for energies,
    /// which no resistance changes.
    battery: BatterySource,
    motor: DcMotor,
    ports: Vec<PortSpec>,
    energy_tolerance_j: f64,
}

impl fmt::Debug for CommutatedBatteryMotorDrive {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("CommutatedBatteryMotorDrive")
            .field("parameters", &self.parameters)
            .field("storage", &self.storage)
            .field("program", &self.program)
            .field("resistor_laws", &self.resistor_laws)
            .field("switch_laws", &self.switch_laws)
            .field("energy_tolerance_j", &self.energy_tolerance_j)
            .finish_non_exhaustive()
    }
}

impl CommutatedBatteryMotorDrive {
    /// The drive of `parameters` whose switches store charge by `storage`,
    /// running `program` once per step under the owner's `electrical` and
    /// `cycle` controls.
    ///
    /// `energy_tolerance_j` \[J\] is the largest balance defect a step may
    /// leave beyond the rounding of its own arithmetic: the owner's solve
    /// residuals times its efforts, which this boundary cannot bound by
    /// itself. A defect within it is booked as numerical loss.
    ///
    /// # Errors
    /// The owner's refusal of the curve, the parameters, the storage or the
    /// program ([`CommutatedDriveError::Model`],
    /// [`CommutatedDriveError::Step`]); one auxiliary load entry per
    /// interval is required; a tolerance that is not positive and finite.
    pub fn new(
        parameters: DriveParameters,
        storage: DeviceStorage,
        program: GateProgram,
        electrical: PowerLoadControls,
        cycle: SwitchCycleControls,
        energy_tolerance_j: f64,
    ) -> Result<Self, CommutatedDriveError> {
        if !(energy_tolerance_j.is_finite() && energy_tolerance_j > 0.0) {
            return Err(invalid("energy tolerance must be positive and finite"));
        }
        if program.auxiliary_loads.len() != program.intervals.len() {
            return Err(invalid("one auxiliary load entry per interval"));
        }
        let model = CommutatedDriveError::Model;
        let curve =
            BatteryOcvCurve::new(parameters.capacity_c, &parameters.soc_voltage).map_err(model)?;
        let battery = BatterySource::new(curve.clone(), parameters.battery).map_err(model)?;
        let motor = DcMotor::new(parameters.motor).map_err(model)?;
        let storage = match storage {
            DeviceStorage::Capacitance(output_capacitance_f) => {
                Storage::Capacitance(CapacitiveCommutationParameters {
                    converter: parameters.converter,
                    output_capacitance_f,
                })
            }
            DeviceStorage::Charge {
                curves,
                maximum_cells_per_chord,
            } => Storage::Charge(ChargeCommutationParameters {
                converter: parameters.converter,
                output_charge_curves: curves,
                maximum_cells_per_chord,
            }),
        };
        let mut ports = Vec::with_capacity(13);
        for name in ["rotor", "stator"] {
            ports.push(
                PortSpec::power(
                    StableId::new(name).map_err(CommutatedDriveError::Identifier)?,
                    PortKind::RotationalTorqueAngularVelocity,
                    PortRole::EffortIn,
                    NonZeroUsize::MIN,
                )
                .map_err(CommutatedDriveError::Port)?,
            );
        }
        let loaded = program.auxiliary_loads.iter().any(Option::is_some);
        for name in HEAT_PORTS.into_iter().chain(loaded.then_some(LOAD_PORT)) {
            let id = StableId::new(name).map_err(CommutatedDriveError::Identifier)?;
            ports.push(PortSpec::heat(id, PortRole::EffortIn).map_err(CommutatedDriveError::Port)?);
        }
        PortLayout::new(&ports).map_err(CommutatedDriveError::Port)?;
        let mut drive = Self {
            parameters,
            curve,
            storage,
            program,
            period_s: 0.0,
            electrical,
            cycle,
            resistor_laws: Default::default(),
            switch_laws: Default::default(),
            battery,
            motor,
            ports,
            energy_tolerance_j,
        };
        // The owner admits the program and the storage at the authored
        // resistances and switches; its period is every step's length.
        let p = &drive.parameters;
        let original = drive.original(&p.converter, p.battery, p.lead_resistance_ohm, p.motor)?;
        drive.period_s = drive
            .prepare(original, p.converter.devices, || false)?
            .period_s();
        Ok(drive)
    }

    /// The same drive with `resistor` following `law` at the temperature
    /// its heat port receives, in every proposal.
    ///
    /// # Errors
    /// That resistance already follows a law, or the law follows the
    /// operating point ([`ResistanceLaw::follows_operating_point`]), which
    /// the drive does not solve: it asks its laws at temperature alone.
    pub fn following(
        mut self,
        resistor: DriveResistor,
        law: Arc<dyn ResistanceLaw>,
    ) -> Result<Self, CommutatedDriveError> {
        if law.follows_operating_point() {
            return Err(invalid(
                "a drive resistance follows temperature, not the operating point",
            ));
        }
        if self.resistor_laws[resistor.index()].replace(law).is_some() {
            return Err(invalid("a resistance follows one law only"));
        }
        Ok(self)
    }

    /// The same drive with the switch at `position` following `law` at the
    /// temperature its heat port receives, in every proposal.
    ///
    /// # Errors
    /// That switch already follows a law.
    pub fn switch_following(
        mut self,
        position: SwitchPosition,
        law: Arc<dyn SwitchLaw>,
    ) -> Result<Self, CommutatedDriveError> {
        if self.switch_laws[position.index()].replace(law).is_some() {
            return Err(invalid("a switch follows one law only"));
        }
        Ok(self)
    }

    /// The parameters the drive was built with.
    #[must_use]
    pub const fn parameters(&self) -> &DriveParameters {
        &self.parameters
    }

    /// The program's period \[s\], the length of every step.
    #[must_use]
    pub const fn period_s(&self) -> f64 {
        self.period_s
    }

    fn energy(&self, state: &[f64; 11]) -> Result<f64, CommutatedDriveError> {
        let converter = core::array::from_fn(|i| state[3 + i]);
        let parts = [
            self.battery.system().hamiltonian(&state[..3]),
            self.storage
                .energy(&converter)
                .map_err(CommutatedDriveError::Model)?,
            self.motor.system().hamiltonian(&state[8..]),
        ];
        let energy = parts[0] + parts[1] + parts[2];
        if state.iter().all(|value| value.is_finite()) && energy.is_finite() {
            Ok(energy)
        } else {
            Err(invalid(
                "the drive's state is outside its laws or its energy is not finite",
            ))
        }
    }

    /// The owner's original drive at these resistances. Its averaged bridge
    /// is replaced by the switches when the cycle is prepared.
    fn original(
        &self,
        converter: &CommutationParameters,
        battery: BatteryParameters,
        lead_resistance_ohm: f64,
        motor: DcMotorParameters,
    ) -> Result<BatteryMotorDrive, CommutatedDriveError> {
        let model = CommutatedDriveError::Model;
        let on = converter.devices.map(|device| device.on_resistance_ohm);
        AveragedBuckBoost::new(
            BuckBoostParameters {
                inductance_h: converter.inductance_h,
                bus_capacitance_f: converter.bus_capacitance_f,
                winding_resistance_ohm: converter.winding_resistance_ohm,
                input_on_resistance_ohm: [on[0], on[1]],
                output_on_resistance_ohm: [on[2], on[3]],
            },
            BridgeModulation {
                input_high_fraction: 0.0,
                output_high_fraction: 0.0,
            },
        )
        .map_err(model)?
        .with_battery_source(
            BatterySource::new(self.curve.clone(), battery).map_err(model)?,
            lead_resistance_ohm,
        )
        .map_err(model)?
        .with_dc_motor(DcMotor::new(motor).map_err(model)?)
        .map_err(model)
    }

    fn prepare(
        &self,
        original: BatteryMotorDrive,
        switches: [SwitchDevice; 4],
        cancelled: impl FnMut() -> bool,
    ) -> Result<Prepared, CommutatedDriveError> {
        let program = &self.program;
        match &self.storage {
            Storage::Capacitance(p) => PreparedCapacitiveBatteryMotorCycle::new_controlled(
                original,
                switches,
                p.output_capacitance_f,
                &program.intervals,
                program.substeps_per_interval,
                self.electrical,
                self.cycle,
                cancelled,
            )
            .map(Prepared::Capacitance),
            Storage::Charge(p) => PreparedChargeBatteryMotorCycle::new_charge_controlled(
                original,
                switches,
                p.output_charge_curves.clone(),
                p.maximum_cells_per_chord,
                &program.intervals,
                program.substeps_per_interval,
                self.electrical,
                self.cycle,
                cancelled,
            )
            .map(Prepared::Charge),
        }
        .map_err(step_error)
    }

    /// The converter, battery, lead and motor of a step whose heat ports
    /// received `kelvin`.
    fn laws_at(
        &self,
        kelvin: &[f64],
    ) -> Result<
        (
            CommutationParameters,
            BatteryParameters,
            f64,
            DcMotorParameters,
        ),
        CommutatedDriveError,
    > {
        let p = &self.parameters;
        let mut ohms = [
            p.converter.winding_resistance_ohm,
            p.battery.series_resistance_ohm,
            p.battery.polarization_resistance_ohm,
            p.lead_resistance_ohm,
            p.motor.winding_resistance_ohm,
        ];
        for resistor in DriveResistor::ALL {
            let Some(law) = &self.resistor_laws[resistor.index()] else {
                continue;
            };
            let port = HEAT_PORTS[resistor.heat_port()];
            let temperature = kelvin[resistor.heat_port()];
            let answer = law
                .ohms(temperature)
                .map_err(|reason| CommutatedDriveError::Law { port, reason })?;
            if !(answer.is_finite() && answer > 0.0) {
                return Err(CommutatedDriveError::Law {
                    port,
                    reason: format!("{answer} ohm at {temperature} K is not positive and finite"),
                });
            }
            ohms[resistor.index()] = answer;
        }
        let mut converter = p.converter;
        converter.winding_resistance_ohm = ohms[0];
        for (i, law) in self.switch_laws.iter().enumerate() {
            let Some(law) = law else { continue };
            let port = HEAT_PORTS[1 + i];
            converter.devices[i] = law
                .switch(kelvin[1 + i])
                .map_err(|reason| CommutatedDriveError::Law { port, reason })?;
        }
        converter.validate().map_err(CommutatedDriveError::Model)?;
        let battery = BatteryParameters {
            series_resistance_ohm: ohms[1],
            polarization_resistance_ohm: ohms[2],
            ..p.battery
        };
        let motor = DcMotorParameters {
            winding_resistance_ohm: ohms[4],
            ..p.motor
        };
        Ok((converter, battery, ohms[3], motor))
    }

    /// The step's receipt, heat rates and mean shaft velocities, refused if
    /// the owner's accounts, its native record and the port works disagree.
    fn account(
        &self,
        step: &CapacitiveBatteryMotorCycleStep,
        old: &[f64; 11],
        next: &[f64; 11],
        torques: [f64; 2],
        kelvin: &[f64],
        dt_s: f64,
    ) -> Result<(Receipt, [f64; 11], [f64; 2]), CommutatedDriveError> {
        let heats = port_heats(step);
        if heats.iter().any(|heat| !(heat.is_finite() && *heat >= 0.0)) {
            return Err(invalid(
                "a heat or the load's work is not finite and nonnegative",
            ));
        }
        let rates = heats.map(|heat| heat / dt_s);
        if rates
            .iter()
            .zip(heats)
            .any(|(rate, heat)| !rate.is_finite() || (heat != 0.0 && *rate == 0.0))
        {
            return Err(invalid("a heat rate is not representable"));
        }
        let speeds = step.mean_shaft_speed_rad_s;
        let port_work_j = [0, 1].map(|i| torques[i] * speeds[i] * dt_s);
        let shaft_j = step.shaft_work_j[0] + step.shaft_work_j[1];
        let heat_j = step.heat_joules().iter().sum::<f64>();
        let stored_change_j = step.total_storage_change_j();
        let numerical_loss_j = -step.balance_residual_j();
        let energy_old = self.energy(old)?;
        let energy_new = self.energy(next)?;
        let magnitude = energy_old.abs()
            + energy_new.abs()
            + step.shaft_work_j[0].abs()
            + step.shaft_work_j[1].abs()
            + heat_j
            + step.load_work_j
            + step.native.delta_h.abs()
            + step.native.supplied.abs()
            + step.native.dissipated.abs();
        let rounding = arithmetic_rounding_floor(64 * step.integration_steps.max(1), magnitude);
        let allowance = rounding + self.energy_tolerance_j;
        // The owner's accounts against its native record and the port works.
        let defect_j = [
            numerical_loss_j,
            stored_change_j - step.native.delta_h,
            heat_j - step.native.dissipated,
            shaft_j - step.load_work_j - step.native.supplied,
            port_work_j[0] - step.shaft_work_j[0],
            port_work_j[1] - step.shaft_work_j[1],
        ]
        .into_iter()
        .map(f64::abs)
        .fold(0.0, f64::max);
        let entropy = heats
            .iter()
            .zip(kelvin)
            .map(|(heat, temperature)| heat / temperature)
            .sum::<f64>();
        if ![
            rounding,
            defect_j,
            stored_change_j,
            entropy,
            speeds[0],
            speeds[1],
        ]
        .iter()
        .all(|value| value.is_finite())
        {
            return Err(invalid("the drive's energy account is not finite"));
        }
        if defect_j > allowance {
            return Err(CommutatedDriveError::Balance {
                defect_j,
                tolerance_j: allowance,
            });
        }
        // The owner bounds every substep's Newton residual by its own budget
        // and carries no single certificate for the cycle: the receipt states
        // the largest substep residual as met.
        let receipt = Receipt {
            stored_change_j,
            dissipated_j: heat_j + step.load_work_j,
            numerical_loss_j,
            entropy_generated_j_per_k: entropy,
            energy_tolerance_j: rounding,
            residual: step.native.solver_residual,
            tolerance: step.native.solver_residual,
        };
        Ok((receipt, rates, speeds))
    }

    /// The identity words of the storage, the program and the controls.
    fn program_words(&self, words: &mut Vec<u64>) {
        match &self.storage {
            Storage::Capacitance(s) => {
                words.push(0);
                words.extend(s.output_capacitance_f.map(f64::to_bits));
            }
            Storage::Charge(s) => {
                words.push(1);
                words.push(s.maximum_cells_per_chord as u64);
                for curve in &s.output_charge_curves {
                    words.push(curve.knots().len() as u64);
                    for &(q, v) in curve.knots() {
                        words.extend([q.to_bits(), v.to_bits()]);
                    }
                }
            }
        }
        let program = &self.program;
        words.extend([
            program.intervals.len() as u64,
            program.substeps_per_interval as u64,
        ]);
        for (interval, load) in program.intervals.iter().zip(&program.auxiliary_loads) {
            let gates = interval
                .gates
                .on
                .iter()
                .enumerate()
                .fold(0_u64, |mask, (i, &on)| mask | (u64::from(on) << i));
            words.extend([gates, interval.duration_s.to_bits()]);
            match load {
                None => words.push(0),
                Some(load) => {
                    words.push(1);
                    words.extend(
                        [
                            load.requested_power_w(),
                            load.maximum_current_a(),
                            load.foldback_voltage_v(),
                            load.maximum_voltage_v(),
                        ]
                        .map(f64::to_bits),
                    );
                }
            }
        }
        let (e, s) = (self.electrical, self.cycle);
        words.extend([
            e.solve.relative.to_bits(),
            e.solve.absolute.to_bits(),
            e.solve.max_iterations as u64,
            u64::from(e.backtracking_halvings),
            e.absolute_energy_tolerance_j.to_bits(),
            e.relative_energy_tolerance.to_bits(),
            s.max_intervals as u64,
            s.absolute_energy_tolerance_j.to_bits(),
            s.relative_energy_tolerance.to_bits(),
        ]);
    }
}

/// A law's identity, its words counted so that one law's words cannot read
/// as another's; a single zero where there is no law.
fn counted(words: &mut Vec<u64>, identity: Option<impl FnOnce(&mut Vec<u64>)>) {
    let Some(identity) = identity else {
        words.push(0);
        return;
    };
    words.push(1);
    let count = words.len();
    words.push(0);
    identity(words);
    words[count] = (words.len() - count - 1) as u64;
}

/// The heat of each heat port over the step \[J\]: the owner's fourteen
/// accounts with each switch's channel and diode together, then the load's
/// work.
fn port_heats(step: &CapacitiveBatteryMotorCycleStep) -> [f64; 11] {
    let h = step.heat_joules();
    [
        h[0],
        h[1] + h[10],
        h[2] + h[11],
        h[3] + h[12],
        h[4] + h[13],
        h[5],
        h[6],
        h[7],
        h[8],
        h[9],
        step.load_work_j,
    ]
}

impl Component for CommutatedBatteryMotorDrive {
    type State = [f64; 11];
    type Error = CommutatedDriveError;

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-system::commutated_drive::CommutatedBatteryMotorDrive",
        state_family: "fs_phs commutated battery-motor drive coordinates",
        schema_grammar: "[f64; 11]: battery chemical, polarization and bus charges [C]; converter \
            flux linkage [Wb], two bus and two switching-node charges [C]; motor flux linkage \
            [Wb], rotor and stator angular momenta [N m s]; finite energy",
        codec_grammar: "one word per coordinate, the f64 bits, in schema order",
        codec_version: 1,
    });

    fn encode_state(&self, state: &[f64; 11], words: &mut Vec<u64>) -> bool {
        if self.energy(state).is_err() {
            return false;
        }
        words.extend(state.iter().map(|value| value.to_bits()));
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<[f64; 11]> {
        let words: &[u64; 11] = words.try_into().ok()?;
        let state = words.map(f64::from_bits);
        self.energy(&state).is_ok().then_some(state)
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        let p = &self.parameters;
        words.push(p.capacity_c.to_bits());
        words.push(p.soc_voltage.len() as u64);
        for &(soc, volts) in &p.soc_voltage {
            words.extend([soc.to_bits(), volts.to_bits()]);
        }
        let b = p.battery;
        let c = &p.converter;
        let m = p.motor;
        let values = [
            b.series_resistance_ohm,
            b.polarization_resistance_ohm,
            b.polarization_capacitance_f,
            b.bus_capacitance_f,
            p.lead_resistance_ohm,
            c.inductance_h,
            c.bus_capacitance_f[0],
            c.bus_capacitance_f[1],
            c.winding_resistance_ohm,
            m.inductance_h,
            m.winding_resistance_ohm,
            m.rotor_inertia_kg_m2,
            m.stator_inertia_kg_m2,
            m.coupling_nm_per_a,
            m.bearing_drag_nm_s_per_rad,
        ];
        words.extend(values.map(f64::to_bits));
        for d in &c.devices {
            words.extend(
                [
                    d.on_resistance_ohm,
                    d.off_resistance_ohm,
                    d.body_diode.forward_voltage_v,
                    d.body_diode.differential_resistance_ohm,
                    d.bounds.maximum_reverse_voltage_v,
                    d.bounds.maximum_channel_current_a,
                    d.bounds.maximum_diode_current_a,
                ]
                .map(f64::to_bits),
            );
        }
        self.program_words(words);
        // Each law's words are counted, so that one law's words cannot read
        // as another's.
        for law in &self.resistor_laws {
            counted(
                words,
                law.as_ref()
                    .map(|law| |w: &mut Vec<u64>| law.identity_words(w)),
            );
        }
        for law in &self.switch_laws {
            counted(
                words,
                law.as_ref()
                    .map(|law| |w: &mut Vec<u64>| law.identity_words(w)),
            );
        }
        words.push(self.energy_tolerance_j.to_bits());
    }

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, state: &[f64; 11]) -> Result<f64, CommutatedDriveError> {
        self.energy(state)
    }

    fn propose(
        &self,
        old: &[f64; 11],
        dt_s: f64,
        inputs: &[f64],
        new: &mut [f64; 11],
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, CommutatedDriveError> {
        cx.checkpoint()
            .map_err(|_| CommutatedDriveError::Cancelled)?;
        let width = self.ports.len();
        if !(dt_s.is_finite() && dt_s > 0.0)
            || inputs.len() != width
            || outputs.len() != width
            || inputs.iter().any(|value| !value.is_finite())
            || inputs[2..].iter().any(|&kelvin| kelvin <= 0.0)
        {
            return Err(invalid(
                "positive duration, complete finite ports and positive temperatures required",
            ));
        }
        let allowance_s = 64.0 * f64::EPSILON * self.period_s;
        if (dt_s - self.period_s).abs() > allowance_s {
            return Err(invalid("a step of the drive is one period of its program"));
        }
        let (torques, kelvin) = ([inputs[0], inputs[1]], &inputs[2..]);
        let (converter, battery, lead, motor) = self.laws_at(kelvin)?;
        let original = self.original(&converter, battery, lead, motor)?;
        let mut native = self.prepare(original, converter.devices, || cx.is_cancel_requested())?;
        let commands: Vec<_> = self
            .program
            .auxiliary_loads
            .iter()
            .map(|&auxiliary_load| BatteryMotorSwitchCommand {
                shaft_torque_nm: torques,
                auxiliary_load,
            })
            .collect();
        let mut next = [0.0; 11];
        let step = native
            .step(old, &commands, &mut next, || cx.is_cancel_requested())
            .map_err(step_error)?;

        let (receipt, rates, speeds) = self.account(&step, old, &next, torques, kelvin, dt_s)?;
        cx.checkpoint()
            .map_err(|_| CommutatedDriveError::Cancelled)?;
        outputs[0] = speeds[0];
        outputs[1] = speeds[1];
        for (output, rate) in outputs[2..].iter_mut().zip(rates) {
            *output = -rate;
        }
        *new = next;
        Ok(Proposal::Step(receipt))
    }
}
