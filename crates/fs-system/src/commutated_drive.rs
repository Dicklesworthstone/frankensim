//! A battery-fed commutated motor drive whose resistances and switches heat
//! what their heat ports are joined to.
//!
//! [`CommutatedBatteryMotorDrive`] puts the drive of
//! `fs_phs::circuit::converter::commutation::battery_motor` and of its
//! `capacitance` lanes behind `fs_ports::component::Component`. It joins the
//! equivalent-circuit battery and its lead to a two-leg converter and a DC
//! motor whose rotor and stator both turn. The converter has four switches,
//! each a channel with its body diode, with or without output charge storage,
//! or with driven gates ([`DriveLane`]). Nothing is integrated
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
//! - In the gate lane ([`GateStorage`]), eight rails, electrical,
//!   [`PortRole::EffortIn`], in [`RAIL_PORTS`] order: each switch's on rail,
//!   then its off rail. Each receives its voltage \[V\], held through the
//!   program, and returns the mean current its gate drew from it \[A\], so
//!   voltage times mean current times the step is the rail's work.
//! - Ten heat ports, [`PortRole::EffortIn`], in [`HEAT_PORTS`] order:
//!   - the converter winding;
//!   - the four switches, a channel and its body diode together;
//!   - the battery's series and polarization resistances, the lead and the
//!     motor winding;
//!   - the bearing.
//!
//!   In the gate lane, four driver heat ports follow, in
//!   [`DRIVER_HEAT_PORTS`] order: each gate's on and off driver resistances
//!   together. If the program has an auxiliary load, a last heat port,
//!   `heat_auxiliary_load`, carries the work the load took. Each heat port
//!   receives the temperature \[K\] where its heat leaves and returns minus
//!   that heat over the step divided by the step \[W\].
//!
//! State: the owner's coordinates, which the drive's lane fixes:
//! - the remaining chemical charge, polarization charge and battery bus
//!   charge \[C\];
//! - the converter's flux linkage \[Wb\] and its two bus charges, and with
//!   switch storage its two switching-node charges \[C\];
//! - the motor's flux linkage \[Wb\] and its rotor and stator angular
//!   momenta \[N m s\].
//!
//! That is nine coordinates without switch storage ([`NoDeviceStorage`]),
//! eleven with it ([`DeviceStorage`], the default lane), and fifteen with
//! driven gates ([`GateStorage`]): there the node charges are grouped with
//! their high-side gate electrodes and four gate charges follow. Its stored
//! energy is the owner's Hamiltonian, the switches' output and gate storage
//! included; no resistance, switch or channel changes it.
//!
//! Switch storage ([`DeviceStorage`]) is either four constant output
//! capacitances or four voltage-charge curves whose energy is the integral
//! of voltage over charge. Either is fixed: storage that changed with
//! temperature would need its entropy, which this Component does not hold.
//!
//! Heat and work:
//! - Heat is the owner's fourteen accounts, routed to the heat ports. The
//!   auxiliary load's work leaves through its own port.
//! - Supplied work is the shaft work and, in the gate lane, the rails' work.
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
//! is never clamped. In the gate lane a [`ChannelLaw`] gives each switch's
//! channel at that temperature too, and a driver resistance following a
//! [`ResistanceLaw`] ([`CommutatedBatteryMotorDrive::driver_following`]) is
//! given its answer at its gate's driver heat port.
//! Inductances, capacitances, the OCV, the motor's coupling, its inertias
//! and the bearing drag stay fixed.
//!
//! Programs: a drive runs the program it was built with every step, unless
//! it was given further programs with
//! [`CommutatedBatteryMotorDrive::scheduled`]. It then has a signal input
//! `program` after its heat ports. The signal is the index of the step's
//! program, an exact whole number held over the step, as `IdealSwitch`
//! takes its `closed` command; a `SignalSource` of `Waveform::held_table`
//! schedules it. Every program has the drive's period.
//!
//! Class: [`Class::Smooth`], with no tangent, so the kernel differences the
//! port map. The map is smooth except where, as the inputs move, a body
//! diode starts or stops conducting, a switch reaches an operating limit, or
//! the auxiliary load changes branch. The module's own thermal Newton
//! differences the same map.
//!
//! Not claimed: torques and rail voltages that change within the program;
//! reverse recovery and avalanche; temperature-dependent storage, OCV or
//! bearing drag.

use core::fmt;
use core::marker::PhantomData;
use core::num::NonZeroUsize;
use std::sync::Arc;

use fs_exec::Cx;
use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_phs::circuit::battery::{BatteryOcvCurve, BatteryParameters, BatterySource};
use fs_phs::circuit::capacitance::CapacitorChargeCurve;
use fs_phs::circuit::converter::battery_fed::motor::BatteryMotorDrive;
use fs_phs::circuit::converter::commutation::battery_motor::PreparedCommutatedBatteryMotorCycle;
use fs_phs::circuit::converter::commutation::capacitance::CapacitiveCommutationParameters;
use fs_phs::circuit::converter::commutation::capacitance::battery_motor::{
    PreparedCapacitiveBatteryMotorCycle, PreparedChargeBatteryMotorCycle,
    PreparedGateBatteryMotorCycle,
};
use fs_phs::circuit::converter::commutation::capacitance::charge::ChargeCommutationParameters;
use fs_phs::circuit::converter::commutation::capacitance::gate::GateCommutationParameters;
use fs_phs::circuit::converter::commutation::gate::{
    GateChannelLaw, GateChannelPoint, GateChannelResponse, GateDriver,
};
use fs_phs::circuit::converter::commutation::{
    CommutationInterval, CommutationParameters, SwitchDevice,
};
use fs_phs::circuit::converter::power_load::{PowerLoadControls, RegulatedPowerLoad};
use fs_phs::circuit::converter::switching::SwitchCycleControls;
use fs_phs::circuit::converter::switching::battery_motor::BatteryMotorSwitchCommand;
use fs_phs::circuit::converter::{AveragedBuckBoost, BridgeModulation, BuckBoostParameters};
use fs_phs::circuit::motor::{DcMotor, DcMotorParameters};
use fs_phs::{PhsError, PreparedStepError, PreparedStepRecord, arithmetic_rounding_floor};
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

/// The gate lane's rails, after `rotor` and `stator`: each switch's on rail,
/// then its off rail, in the owner's switch order.
pub const RAIL_PORTS: [&str; 8] = [
    "gate_input_high_on",
    "gate_input_high_off",
    "gate_input_low_on",
    "gate_input_low_off",
    "gate_output_high_on",
    "gate_output_high_off",
    "gate_output_low_on",
    "gate_output_low_off",
];

/// The gate lane's driver heat ports, after [`HEAT_PORTS`]: one per gate,
/// for its on and off resistances together.
pub const DRIVER_HEAT_PORTS: [&str; 4] = [
    "heat_gate_input_high",
    "heat_gate_input_low",
    "heat_gate_output_high",
    "heat_gate_output_low",
];

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

/// A driven switch's channel as a function of the temperature of the package
/// it heats, in the gate lane ([`GateStorage`]). The channel's current comes
/// from this law alone; a [`SwitchLaw`] there gives the body diode and the
/// limits.
pub trait ChannelLaw: fmt::Debug + Send + Sync {
    /// The channel's response to its drain-source and gate-source voltages
    /// at `kelvin`.
    ///
    /// # Errors
    /// Why the law has no answer there. The drive then refuses the step with
    /// [`CommutatedDriveError::Law`]; it never clamps the temperature.
    fn channel(&self, kelvin: f64) -> Result<GateChannelLaw, String>;

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

/// Switches without output storage: the owner's nine-coordinate lane,
/// `fs_phs::circuit::converter::commutation::battery_motor`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct NoDeviceStorage;

/// One of a gate driver's two rails: the one its command selects to turn the
/// switch on, or the one that turns it off.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GateRail {
    /// The rail that turns the switch on.
    On,
    /// The rail that turns the switch off.
    Off,
}

/// A gate's driver: its on and off resistances and the current it supports.
/// Its rail voltages are the efforts its two rail ports receive.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GateDriverResistances {
    /// Between the on rail and the gate \[ohm\].
    pub on_resistance_ohm: f64,
    /// Between the off rail and the gate \[ohm\].
    pub off_resistance_ohm: f64,
    /// The largest current either branch supports \[A\].
    pub maximum_current_a: f64,
}

/// Switches with driven gates: the owner's fifteen-coordinate lane,
/// `PreparedGateBatteryMotorCycle`. Each switch keeps its intrinsic
/// drain-source charge curve, a gate-source and a gate-drain (Miller)
/// capacitance, and a [`ChannelLaw`]; each gate a finite driver. The storage
/// is fixed.
#[derive(Debug, Clone)]
pub struct GateStorage {
    /// The intrinsic drain-source curves, charge \[C\] against voltage \[V\].
    pub intrinsic_drain_charge: [CapacitorChargeCurve; 4],
    /// Positive gate-source capacitances \[F\].
    pub gate_source_capacitance_f: [f64; 4],
    /// Nonnegative gate-drain capacitances \[F\]; zero is the exact limit
    /// without Miller coupling.
    pub gate_drain_capacitance_f: [f64; 4],
    /// The owner's bound on the cells one chord may cross.
    pub maximum_cells_per_chord: usize,
    /// Each switch's channel, asked at its heat port's temperature.
    pub channels: [Arc<dyn ChannelLaw>; 4],
    /// Each gate's driver.
    pub drivers: [GateDriverResistances; 4],
}

/// fs-phs builds the gate storage from a whole parameter set, though neither
/// the channel nor the drivers enter its energy. The drive's energy, and its
/// admission of each program when built, use this channel, which refuses
/// every evaluation, so no step can run on it.
#[derive(Debug, Clone, Copy)]
struct StandInChannel;

impl GateChannelResponse for StandInChannel {
    fn gate_voltage_range_v(&self) -> [f64; 2] {
        [-1.0, 1.0]
    }

    fn evaluate(&self, _: f64, _: f64) -> Result<GateChannelPoint, PhsError> {
        Err(PhsError::Unadmitted {
            what: "the drive's stand-in channel is never evaluated",
        })
    }
}

/// What the switches store, which fixes the drive's state:
/// [`NoDeviceStorage`], nine coordinates, [`DeviceStorage`], eleven, or
/// [`GateStorage`], fifteen. Sealed: the lanes are the owner's.
pub trait DriveLane: lane::Sealed + fmt::Debug + Clone + Send + Sync + 'static {}

impl DriveLane for NoDeviceStorage {}
impl DriveLane for DeviceStorage {}
impl DriveLane for GateStorage {}

mod lane {
    use super::{
        Arc, CapacitiveCommutationParameters, ChannelLaw, ChargeCommutationParameters,
        CommutatedDriveError, CommutationParameters, DeviceStorage, GateChannelLaw,
        GateCommutationParameters, GateDriver, GateStorage, NoDeviceStorage, ResistanceLaw,
        StandInChannel, StateIdentityCharterV2, averaged, fmt, invalid,
    };

    /// The lane's state and storage law, for the drive alone.
    pub trait Sealed {
        /// The owner's coordinates.
        type State: Copy
            + Default
            + AsRef<[f64]>
            + AsMut<[f64]>
            + fmt::Debug
            + PartialEq
            + Send
            + Sync
            + 'static;
        /// The state's identity charter.
        const CHARTER: StateIdentityCharterV2;
        /// Whether the lane drives its gates, with rail and driver heat ports.
        const GATED: bool = false;
        /// The switches' storage law on `converter`.
        fn law(self, converter: CommutationParameters) -> Result<Storage, CommutatedDriveError>;
    }

    /// The gate lane's storage, with the stand-in channel and the drivers at
    /// zero rail voltage, and its channel laws.
    #[derive(Debug, Clone)]
    pub struct Gated {
        /// The storage the drive's energy is taken from.
        pub storage: GateCommutationParameters<GateChannelLaw>,
        /// Each switch's channel.
        pub channels: [Arc<dyn ChannelLaw>; 4],
        /// Each gate's on and off driver resistance laws, if any.
        pub driver_laws: [[Option<Arc<dyn ResistanceLaw>>; 2]; 4],
    }

    /// The converter's storage law with the converter it belongs to.
    #[derive(Debug, Clone)]
    pub enum Storage {
        /// No switch storage.
        Plain(CommutationParameters),
        /// Constant output capacitances.
        Capacitance(CapacitiveCommutationParameters),
        /// Voltage-charge curves.
        Charge(ChargeCommutationParameters),
        /// Driven gates.
        Gate(Box<Gated>),
    }

    impl Storage {
        /// The converter's stored energy at its coordinates.
        pub fn energy(&self, converter: &[f64]) -> Result<f64, CommutatedDriveError> {
            let model = CommutatedDriveError::Model;
            let shape = || invalid("the converter's coordinates are the lane's");
            match self {
                Self::Plain(p) => {
                    let flux_and_buses: &[f64; 3] = converter.try_into().map_err(|_| shape())?;
                    Ok(averaged(p)
                        .map_err(model)?
                        .system()
                        .hamiltonian(flux_and_buses))
                }
                Self::Capacitance(p) => p
                    .stored_energy_j(converter.try_into().map_err(|_| shape())?)
                    .map_err(model),
                Self::Charge(p) => p
                    .stored_energy_j(converter.try_into().map_err(|_| shape())?)
                    .map_err(model),
                Self::Gate(g) => g
                    .storage
                    .stored_energy_j(converter.try_into().map_err(|_| shape())?)
                    .map_err(model),
            }
        }
    }

    impl Sealed for NoDeviceStorage {
        type State = [f64; 9];
        const CHARTER: StateIdentityCharterV2 = StateIdentityCharterV2 {
            owner: "fs-system::commutated_drive::CommutatedBatteryMotorDrive<NoDeviceStorage>",
            state_family: "fs_phs commutated battery-motor drive coordinates, no switch storage",
            schema_grammar: "[f64; 9]: battery chemical, polarization and bus charges [C]; \
                converter flux linkage [Wb] and two bus charges [C]; motor flux linkage [Wb], \
                rotor and stator angular momenta [N m s]; finite energy",
            codec_grammar: "one word per coordinate, the f64 bits, in schema order",
            codec_version: 1,
        };
        fn law(self, converter: CommutationParameters) -> Result<Storage, CommutatedDriveError> {
            Ok(Storage::Plain(converter))
        }
    }

    impl Sealed for DeviceStorage {
        type State = [f64; 11];
        const CHARTER: StateIdentityCharterV2 = StateIdentityCharterV2 {
            owner: "fs-system::commutated_drive::CommutatedBatteryMotorDrive",
            state_family: "fs_phs commutated battery-motor drive coordinates",
            schema_grammar: "[f64; 11]: battery chemical, polarization and bus charges [C]; \
                converter flux linkage [Wb], two bus and two switching-node charges [C]; motor \
                flux linkage [Wb], rotor and stator angular momenta [N m s]; finite energy",
            codec_grammar: "one word per coordinate, the f64 bits, in schema order",
            codec_version: 1,
        };
        fn law(self, converter: CommutationParameters) -> Result<Storage, CommutatedDriveError> {
            Ok(match self {
                Self::Capacitance(output_capacitance_f) => {
                    Storage::Capacitance(CapacitiveCommutationParameters {
                        converter,
                        output_capacitance_f,
                    })
                }
                Self::Charge {
                    curves,
                    maximum_cells_per_chord,
                } => Storage::Charge(ChargeCommutationParameters {
                    converter,
                    output_charge_curves: curves,
                    maximum_cells_per_chord,
                }),
            })
        }
    }

    impl Sealed for GateStorage {
        type State = [f64; 15];
        const CHARTER: StateIdentityCharterV2 = StateIdentityCharterV2 {
            owner: "fs-system::commutated_drive::CommutatedBatteryMotorDrive<GateStorage>",
            state_family: "fs_phs commutated battery-motor drive coordinates with driven gates",
            schema_grammar: "[f64; 15]: battery chemical, polarization and bus charges [C]; \
                converter flux linkage [Wb], two bus charges, two switching-node charges with \
                their high-side gate electrodes and four gate charges [C]; motor flux linkage \
                [Wb], rotor and stator angular momenta [N m s]; finite energy",
            codec_grammar: "one word per coordinate, the f64 bits, in schema order",
            codec_version: 1,
        };
        const GATED: bool = true;
        fn law(self, converter: CommutationParameters) -> Result<Storage, CommutatedDriveError> {
            let model = CommutatedDriveError::Model;
            let stand_in = GateChannelLaw::new(StandInChannel).map_err(model)?;
            let storage = GateCommutationParameters {
                converter,
                intrinsic_drain_charge_curves: self.intrinsic_drain_charge,
                gate_source_capacitance_f: self.gate_source_capacitance_f,
                gate_drain_capacitance_f: self.gate_drain_capacitance_f,
                channel_conductance: core::array::from_fn(|_| stand_in.clone()),
                gate_drivers: self.drivers.map(|d| GateDriver {
                    on_voltage_v: 0.0,
                    off_voltage_v: 0.0,
                    on_resistance_ohm: d.on_resistance_ohm,
                    off_resistance_ohm: d.off_resistance_ohm,
                    maximum_current_a: d.maximum_current_a,
                }),
                maximum_cells_per_chord: self.maximum_cells_per_chord,
            };
            storage.validate().map_err(model)?;
            Ok(Storage::Gate(Box::new(Gated {
                storage,
                channels: self.channels,
                driver_laws: Default::default(),
            })))
        }
    }
}

use lane::Storage;

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

/// The converter of `converter` with its averaged bridge, whose system holds
/// the converter's storage.
fn averaged(converter: &CommutationParameters) -> Result<AveragedBuckBoost, PhsError> {
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
}

/// What the drive reads of the owner's cycle receipt, the same in every lane.
struct Cycle {
    heat_j: [f64; 14],
    stored_change_j: f64,
    balance_residual_j: f64,
    shaft_work_j: [f64; 2],
    load_work_j: f64,
    mean_shaft_speed_rad_s: [f64; 2],
    integration_steps: usize,
    native: PreparedStepRecord,
    /// The gate lane's driver accounts.
    gate: Option<GateCycle>,
}

/// The gate lane's driver accounts over the program.
struct GateCycle {
    /// Each driver's resistor heat \[J\].
    heat_j: [f64; 4],
    /// Each gate's source work \[J\].
    source_work_j: [f64; 4],
    /// Each gate's charge from its on rail, then its off rail \[C\].
    rail_charge_c: [[f64; 2]; 4],
}

/// A [`Cycle`] from either owner's receipt, whose fields share their names.
macro_rules! cycle {
    ($step:expr) => {{
        let step = $step;
        Cycle {
            heat_j: step.heat_joules(),
            stored_change_j: step.total_storage_change_j(),
            balance_residual_j: step.balance_residual_j(),
            shaft_work_j: step.shaft_work_j,
            load_work_j: step.load_work_j,
            mean_shaft_speed_rad_s: step.mean_shaft_speed_rad_s,
            integration_steps: step.integration_steps,
            native: step.native,
            gate: None,
        }
    }};
}

/// The gate lane's channels and drivers for one step.
type Gates = ([GateChannelLaw; 4], [GateDriver; 4]);

/// The gate lane's rails over a step of `dt_s`: each rail's mean current,
/// each gate's two rail works less the owner's source work for it, and the
/// magnitudes of those works for the rounding floor. A rail's voltage is held
/// through the program, so voltage times mean current times the step is its
/// work. The other lanes have no rails.
fn rail_accounts(gate: Option<&GateCycle>, rails: &[f64], dt_s: f64) -> (Vec<f64>, Vec<f64>, f64) {
    let Some(gate) = gate else {
        return (Vec::new(), Vec::new(), 0.0);
    };
    let currents: Vec<f64> = gate
        .rail_charge_c
        .as_flattened()
        .iter()
        .map(|charge| charge / dt_s)
        .collect();
    let work_j: Vec<f64> = rails
        .iter()
        .zip(&currents)
        .map(|(volts, amperes)| volts * amperes * dt_s)
        .collect();
    let defects_j = work_j
        .as_chunks::<2>()
        .0
        .iter()
        .zip(gate.source_work_j)
        .map(|(&[on, off], source)| on + off - source)
        .collect();
    let magnitude_j = gate
        .source_work_j
        .iter()
        .chain(&work_j)
        .map(|w| w.abs())
        .sum();
    (currents, defects_j, magnitude_j)
}

/// A prepared native cycle of the drive's lane.
enum Prepared {
    Plain(PreparedCommutatedBatteryMotorCycle),
    Capacitance(PreparedCapacitiveBatteryMotorCycle),
    Charge(PreparedChargeBatteryMotorCycle),
    Gate(Box<PreparedGateBatteryMotorCycle<GateChannelLaw>>),
}

impl Prepared {
    fn period_s(&self) -> f64 {
        match self {
            Self::Plain(native) => native.period_s(),
            Self::Capacitance(native) => native.period_s(),
            Self::Charge(native) => native.period_s(),
            Self::Gate(native) => native.period_s(),
        }
    }

    fn step(
        &mut self,
        old: &[f64],
        commands: &[BatteryMotorSwitchCommand],
        next: &mut [f64],
        cancelled: impl FnMut() -> bool,
    ) -> Result<Cycle, CommutatedDriveError> {
        let shape = || invalid("the state has the lane's coordinates");
        match self {
            Self::Plain(native) => native
                .step_interval_commands_controlled(
                    old.try_into().map_err(|_| shape())?,
                    commands,
                    next.try_into().map_err(|_| shape())?,
                    cancelled,
                )
                .map(|step| cycle!(step)),
            Self::Capacitance(native) => native
                .step_interval_commands_controlled(
                    old.try_into().map_err(|_| shape())?,
                    commands,
                    next.try_into().map_err(|_| shape())?,
                    cancelled,
                )
                .map(|step| cycle!(step)),
            Self::Charge(native) => native
                .step_interval_commands_controlled(
                    old.try_into().map_err(|_| shape())?,
                    commands,
                    next.try_into().map_err(|_| shape())?,
                    cancelled,
                )
                .map(|step| cycle!(step)),
            Self::Gate(native) => native
                .step_interval_commands_controlled(
                    old.try_into().map_err(|_| shape())?,
                    commands,
                    next.try_into().map_err(|_| shape())?,
                    cancelled,
                )
                .map(|step| {
                    let gate = step.converter.gate.as_ref().map(|g| GateCycle {
                        heat_j: g.resistor_heat_j,
                        source_work_j: g.source_work_j,
                        rail_charge_c: g.rail_charge_c,
                    });
                    Cycle {
                        gate,
                        ..cycle!(step)
                    }
                }),
        }
        .map_err(step_error)
    }
}

/// The battery-fed commutated motor drive of fs-phs as one Component, in the
/// lane `L` ([`DriveLane`]).
///
/// See the module documentation for the ports, the state, the heat and the
/// laws.
pub struct CommutatedBatteryMotorDrive<L: DriveLane = DeviceStorage> {
    parameters: DriveParameters,
    curve: BatteryOcvCurve,
    storage: Storage,
    program: GateProgram,
    /// The programs the `program` signal names after the first, if any.
    scheduled: Vec<GateProgram>,
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
    lane: PhantomData<L>,
}

impl<L: DriveLane> fmt::Debug for CommutatedBatteryMotorDrive<L> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("CommutatedBatteryMotorDrive")
            .field("parameters", &self.parameters)
            .field("storage", &self.storage)
            .field("program", &self.program)
            .field("scheduled", &self.scheduled)
            .field("resistor_laws", &self.resistor_laws)
            .field("switch_laws", &self.switch_laws)
            .field("energy_tolerance_j", &self.energy_tolerance_j)
            .finish_non_exhaustive()
    }
}

impl CommutatedBatteryMotorDrive<GateStorage> {
    /// The same drive with the resistance between the gate of `switch` and
    /// its `rail` following `law` at the temperature that gate's driver heat
    /// port receives, in every proposal.
    ///
    /// # Errors
    /// That resistance already follows a law, or the law follows the
    /// operating point, as for [`Self::following`].
    pub fn driver_following(
        mut self,
        switch: SwitchPosition,
        rail: GateRail,
        law: Arc<dyn ResistanceLaw>,
    ) -> Result<Self, CommutatedDriveError> {
        if law.follows_operating_point() {
            return Err(invalid(
                "a gate driver resistance follows temperature, not the operating point",
            ));
        }
        let Storage::Gate(g) = &mut self.storage else {
            return Err(invalid("the gate lane's storage holds its drivers"));
        };
        if g.driver_laws[switch.index()][rail as usize]
            .replace(law)
            .is_some()
        {
            return Err(invalid("a driver resistance follows one law only"));
        }
        Ok(self)
    }
}

impl<L: DriveLane> CommutatedBatteryMotorDrive<L> {
    /// The drive of `parameters` whose switches store charge as `storage`
    /// says, or not at all ([`NoDeviceStorage`]),
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
        storage: L,
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
        let storage = storage.law(parameters.converter)?;
        let loaded = program.auxiliary_loads.iter().any(Option::is_some);
        let ports = port_specs(loaded, false, L::GATED)?;
        let mut drive = Self {
            parameters,
            curve,
            storage,
            program,
            scheduled: Vec::new(),
            period_s: 0.0,
            electrical,
            cycle,
            resistor_laws: Default::default(),
            switch_laws: Default::default(),
            battery,
            motor,
            ports,
            energy_tolerance_j,
            lane: PhantomData,
        };
        // The owner admits the program and the storage at the authored
        // resistances and switches; its period is every step's length.
        let p = &drive.parameters;
        let original = drive.original(&p.converter, p.battery, p.lead_resistance_ohm, p.motor)?;
        drive.period_s = drive
            .prepare(original, &p.converter, &drive.program, None, || false)?
            .period_s();
        Ok(drive)
    }

    /// The same drive choosing each step's program by its `program` signal:
    /// `0` runs the program it was built with, `k` runs `programs[k - 1]`.
    /// It gains the signal input after its heat ports, and the load port if
    /// one of the programs has a load.
    ///
    /// # Errors
    /// No program, or a drive already given its programs; a program without
    /// one auxiliary load entry per interval; the owner's refusal of a
    /// program; a program whose period is not the drive's.
    pub fn scheduled(mut self, programs: Vec<GateProgram>) -> Result<Self, CommutatedDriveError> {
        if programs.is_empty() || !self.scheduled.is_empty() {
            return Err(invalid(
                "a drive is given one nonempty list of further programs",
            ));
        }
        let p = &self.parameters;
        for program in &programs {
            if program.auxiliary_loads.len() != program.intervals.len() {
                return Err(invalid("one auxiliary load entry per interval"));
            }
            let original =
                self.original(&p.converter, p.battery, p.lead_resistance_ohm, p.motor)?;
            let period_s = self
                .prepare(original, &p.converter, program, None, || false)?
                .period_s();
            if (period_s - self.period_s).abs() > 64.0 * f64::EPSILON * self.period_s {
                return Err(invalid("every program has the drive's period"));
            }
        }
        let loaded = core::iter::once(&self.program)
            .chain(&programs)
            .flat_map(|program| &program.auxiliary_loads)
            .any(Option::is_some);
        self.ports = port_specs(loaded, true, L::GATED)?;
        self.scheduled = programs;
        Ok(self)
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

    fn energy(&self, state: &[f64]) -> Result<f64, CommutatedDriveError> {
        // The battery's three coordinates, the converter's, the motor's three.
        let motor = state.len() - 3;
        let parts = [
            self.battery.system().hamiltonian(&state[..3]),
            self.storage.energy(&state[3..motor])?,
            self.motor.system().hamiltonian(&state[motor..]),
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
        averaged(converter)
            .map_err(model)?
            .with_battery_source(
                BatterySource::new(self.curve.clone(), battery).map_err(model)?,
                lead_resistance_ohm,
            )
            .map_err(model)?
            .with_dc_motor(DcMotor::new(motor).map_err(model)?)
            .map_err(model)
    }

    /// The owner of `program` on `original` and `converter`; in the gate
    /// lane, with `gates`' channels and drivers, or when built, the stand-in
    /// channel and drivers at zero rail voltage.
    fn prepare(
        &self,
        original: BatteryMotorDrive,
        converter: &CommutationParameters,
        program: &GateProgram,
        gates: Option<Gates>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Prepared, CommutatedDriveError> {
        let switches = converter.devices;
        match &self.storage {
            Storage::Plain(_) => PreparedCommutatedBatteryMotorCycle::new_controlled(
                original,
                switches,
                &program.intervals,
                program.substeps_per_interval,
                self.electrical,
                self.cycle,
                cancelled,
            )
            .map(Prepared::Plain),
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
            Storage::Gate(g) => {
                let mut parameters = g.storage.clone();
                parameters.converter = *converter;
                if let Some((channels, drivers)) = gates {
                    parameters.channel_conductance = channels;
                    parameters.gate_drivers = drivers;
                }
                PreparedGateBatteryMotorCycle::new_gate_controlled(
                    original,
                    parameters,
                    &program.intervals,
                    program.substeps_per_interval,
                    self.electrical,
                    self.cycle,
                    cancelled,
                )
                .map(|native| Prepared::Gate(Box::new(native)))
            }
        }
        .map_err(step_error)
    }

    /// The gate lane's channels at the temperatures their switches' heat
    /// ports received, and its drivers at the rails' voltages and at what
    /// their laws answer at the driver heat ports' temperatures.
    fn gates_at(
        &self,
        kelvin: &[f64],
        rails: &[f64],
    ) -> Result<Option<Gates>, CommutatedDriveError> {
        let Storage::Gate(g) = &self.storage else {
            return Ok(None);
        };
        let mut channels = Vec::with_capacity(4);
        for (i, law) in g.channels.iter().enumerate() {
            let port = HEAT_PORTS[1 + i];
            let channel = law
                .channel(kelvin[1 + i])
                .map_err(|reason| CommutatedDriveError::Law { port, reason })?;
            channels.push(channel);
        }
        let channels = channels
            .try_into()
            .map_err(|_| invalid("one channel per switch"))?;
        let mut drivers = g.storage.gate_drivers;
        for (k, driver) in drivers.iter_mut().enumerate() {
            driver.on_voltage_v = rails[2 * k];
            driver.off_voltage_v = rails[2 * k + 1];
            let (port, temperature) = (DRIVER_HEAT_PORTS[k], kelvin[HEAT_PORTS.len() + k]);
            let ohms = [
                &mut driver.on_resistance_ohm,
                &mut driver.off_resistance_ohm,
            ];
            for (law, ohms) in g.driver_laws[k].iter().zip(ohms) {
                let Some(law) = law else { continue };
                let answer = law
                    .ohms(temperature)
                    .map_err(|reason| CommutatedDriveError::Law { port, reason })?;
                if !(answer.is_finite() && answer > 0.0) {
                    return Err(CommutatedDriveError::Law {
                        port,
                        reason: format!(
                            "{answer} ohm at {temperature} K is not positive and finite"
                        ),
                    });
                }
                *ohms = answer;
            }
        }
        Ok(Some((channels, drivers)))
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

    /// The step's receipt, its mean shaft velocities and the outputs after
    /// them: in the gate lane each rail's mean current, then minus each heat
    /// port's heat rate. `efforts` are the shafts' torques, then in the gate
    /// lane the rails' voltages. Refused if the owner's accounts, its native
    /// record and the port works disagree.
    fn account(
        &self,
        step: &Cycle,
        old: &[f64],
        next: &[f64],
        efforts: &[f64],
        kelvin: &[f64],
        dt_s: f64,
    ) -> Result<(Receipt, [f64; 2], Vec<f64>), CommutatedDriveError> {
        let (torques, rails) = ([efforts[0], efforts[1]], &efforts[2..]);
        let gate = match (&step.gate, L::GATED) {
            (Some(gate), true) => Some(gate),
            (None, false) => None,
            _ => return Err(invalid("the owner's gate accounts are the lane's")),
        };
        let mut heats = port_heats(step).to_vec();
        if let Some(gate) = gate {
            // The driver heat ports follow the ten.
            heats.splice(10..10, gate.heat_j);
        }
        if heats.iter().any(|heat| !(heat.is_finite() && *heat >= 0.0)) {
            return Err(invalid(
                "a heat or the load's work is not finite and nonnegative",
            ));
        }
        let rates: Vec<f64> = heats.iter().map(|heat| heat / dt_s).collect();
        if rates
            .iter()
            .zip(&heats)
            .any(|(rate, heat)| !rate.is_finite() || (*heat != 0.0 && *rate == 0.0))
        {
            return Err(invalid("a heat rate is not representable"));
        }
        let (currents, gate_defects_j, gate_magnitude_j) = rail_accounts(gate, rails, dt_s);
        let speeds = step.mean_shaft_speed_rad_s;
        let port_work_j = [0, 1].map(|i| torques[i] * speeds[i] * dt_s);
        let shaft_j = step.shaft_work_j[0] + step.shaft_work_j[1];
        let gate_j = gate.map_or(0.0, |g| g.source_work_j.iter().sum());
        let heat_j = step.heat_j.iter().sum::<f64>() + gate.map_or(0.0, |g| g.heat_j.iter().sum());
        let stored_change_j = step.stored_change_j;
        let numerical_loss_j = -step.balance_residual_j;
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
            + step.native.dissipated.abs()
            + gate_magnitude_j;
        let rounding = arithmetic_rounding_floor(64 * step.integration_steps.max(1), magnitude);
        let allowance = rounding + self.energy_tolerance_j;
        // The owner's accounts against its native record and the port works.
        let defect_j = [
            numerical_loss_j,
            stored_change_j - step.native.delta_h,
            heat_j - step.native.dissipated,
            shaft_j + gate_j - step.load_work_j - step.native.supplied,
            port_work_j[0] - step.shaft_work_j[0],
            port_work_j[1] - step.shaft_work_j[1],
        ]
        .into_iter()
        .chain(gate_defects_j)
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
        .chain(&currents)
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
        let outputs = currents
            .into_iter()
            .chain(rates.iter().map(|rate| -rate))
            .collect();
        Ok((receipt, speeds, outputs))
    }

    /// The identity words of the storage, the program and the controls.
    fn program_words(&self, words: &mut Vec<u64>) {
        match &self.storage {
            // The lane without switch storage: a tag of its own.
            Storage::Plain(_) => words.push(2),
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
            // Driven gates: the storage, the drivers, then each channel law.
            Storage::Gate(g) => {
                let s = &g.storage;
                words.push(3);
                words.push(s.maximum_cells_per_chord as u64);
                for curve in &s.intrinsic_drain_charge_curves {
                    words.push(curve.knots().len() as u64);
                    for &(q, v) in curve.knots() {
                        words.extend([q.to_bits(), v.to_bits()]);
                    }
                }
                words.extend(s.gate_source_capacitance_f.map(f64::to_bits));
                words.extend(s.gate_drain_capacitance_f.map(f64::to_bits));
                for d in &s.gate_drivers {
                    let values = [
                        d.on_resistance_ohm,
                        d.off_resistance_ohm,
                        d.maximum_current_a,
                    ];
                    words.extend(values.map(f64::to_bits));
                }
                for law in &g.channels {
                    counted(words, Some(|w: &mut Vec<u64>| law.identity_words(w)));
                }
                for law in g.driver_laws.as_flattened() {
                    counted(
                        words,
                        law.as_ref()
                            .map(|law| |w: &mut Vec<u64>| law.identity_words(w)),
                    );
                }
            }
        }
        gate_program_words(&self.program, words);
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

/// The words of one gate program: its intervals, gates, durations and loads.
fn gate_program_words(program: &GateProgram, words: &mut Vec<u64>) {
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
}

/// The drive's ports: the two shafts, the gate lane's eight rails, the heat
/// ports with the gate lane's four driver heat ports, the load's port if a
/// program has a load, and the `program` signal if it chooses programs.
fn port_specs(
    loaded: bool,
    signalled: bool,
    gated: bool,
) -> Result<Vec<PortSpec>, CommutatedDriveError> {
    let mut ports = Vec::with_capacity(26);
    let shafts = ["rotor", "stator"].map(|name| (name, PortKind::RotationalTorqueAngularVelocity));
    let rails = RAIL_PORTS
        .into_iter()
        .filter(|_| gated)
        .map(|name| (name, PortKind::ElectricalVoltageCurrent));
    for (name, kind) in shafts.into_iter().chain(rails) {
        ports.push(
            PortSpec::power(
                StableId::new(name).map_err(CommutatedDriveError::Identifier)?,
                kind,
                PortRole::EffortIn,
                NonZeroUsize::MIN,
            )
            .map_err(CommutatedDriveError::Port)?,
        );
    }
    let drivers = DRIVER_HEAT_PORTS.into_iter().filter(|_| gated);
    for name in HEAT_PORTS
        .into_iter()
        .chain(drivers)
        .chain(loaded.then_some(LOAD_PORT))
    {
        let id = StableId::new(name).map_err(CommutatedDriveError::Identifier)?;
        ports.push(PortSpec::heat(id, PortRole::EffortIn).map_err(CommutatedDriveError::Port)?);
    }
    if signalled {
        let id = StableId::new("program").map_err(CommutatedDriveError::Identifier)?;
        ports.push(PortSpec::signal_in(id, NonZeroUsize::MIN));
    }
    PortLayout::new(&ports).map_err(CommutatedDriveError::Port)?;
    Ok(ports)
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
fn port_heats(step: &Cycle) -> [f64; 11] {
    let h = step.heat_j;
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

impl<L: DriveLane> Component for CommutatedBatteryMotorDrive<L> {
    type State = L::State;
    type Error = CommutatedDriveError;

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(L::CHARTER);

    fn encode_state(&self, state: &L::State, words: &mut Vec<u64>) -> bool {
        if self.energy(state.as_ref()).is_err() {
            return false;
        }
        words.extend(state.as_ref().iter().map(|value| value.to_bits()));
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<L::State> {
        let mut state = L::State::default();
        if words.len() != state.as_ref().len() {
            return None;
        }
        for (value, word) in state.as_mut().iter_mut().zip(words) {
            *value = f64::from_bits(*word);
        }
        self.energy(state.as_ref()).is_ok().then_some(state)
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
        // A drive without further programs keeps the words above alone.
        if !self.scheduled.is_empty() {
            words.push(self.scheduled.len() as u64);
            for program in &self.scheduled {
                gate_program_words(program, words);
            }
        }
    }

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, state: &L::State) -> Result<f64, CommutatedDriveError> {
        self.energy(state.as_ref())
    }

    fn propose(
        &self,
        old: &L::State,
        dt_s: f64,
        inputs: &[f64],
        new: &mut L::State,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, CommutatedDriveError> {
        cx.checkpoint()
            .map_err(|_| CommutatedDriveError::Cancelled)?;
        let width = self.ports.len();
        // The heat ports begin after the shafts and any rails, and end where
        // the `program` signal, if any, begins.
        let heat_start = 2 + if L::GATED { RAIL_PORTS.len() } else { 0 };
        let heat_end = width - usize::from(!self.scheduled.is_empty());
        if !(dt_s.is_finite() && dt_s > 0.0)
            || inputs.len() != width
            || outputs.len() != heat_end
            || inputs.iter().any(|value| !value.is_finite())
            || inputs[heat_start..heat_end]
                .iter()
                .any(|&kelvin| kelvin <= 0.0)
        {
            return Err(invalid(
                "positive duration, complete finite ports and positive temperatures required",
            ));
        }
        let allowance_s = 64.0 * f64::EPSILON * self.period_s;
        if (dt_s - self.period_s).abs() > allowance_s {
            return Err(invalid("a step of the drive is one period of its program"));
        }
        let program = match inputs.get(heat_end) {
            None => &self.program,
            // An exact whole number naming a program.
            Some(&index) => match (0..=self.scheduled.len()).find(|&k| {
                u32::try_from(k).is_ok_and(|k| f64::from(k).to_bits() == index.to_bits())
            }) {
                Some(0) => &self.program,
                Some(k) => &self.scheduled[k - 1],
                None => {
                    return Err(invalid(
                        "the program signal names one of the drive's programs",
                    ));
                }
            },
        };
        let (efforts, kelvin) = (&inputs[..heat_start], &inputs[heat_start..heat_end]);
        let torques = [efforts[0], efforts[1]];
        let (converter, battery, lead, motor) = self.laws_at(kelvin)?;
        let gates = self.gates_at(kelvin, &efforts[2..])?;
        let original = self.original(&converter, battery, lead, motor)?;
        let cancelled = || cx.is_cancel_requested();
        let mut native = self.prepare(original, &converter, program, gates, cancelled)?;
        let commands: Vec<_> = program
            .auxiliary_loads
            .iter()
            .map(|&auxiliary_load| BatteryMotorSwitchCommand {
                shaft_torque_nm: torques,
                auxiliary_load,
            })
            .collect();
        let mut next = L::State::default();
        let step = native.step(old.as_ref(), &commands, next.as_mut(), || {
            cx.is_cancel_requested()
        })?;

        let (receipt, speeds, rest) =
            self.account(&step, old.as_ref(), next.as_ref(), efforts, kelvin, dt_s)?;
        cx.checkpoint()
            .map_err(|_| CommutatedDriveError::Cancelled)?;
        outputs[0] = speeds[0];
        outputs[1] = speeds[1];
        for (output, value) in outputs[2..].iter_mut().zip(rest) {
            *output = value;
        }
        *new = next;
        Ok(Proposal::Step(receipt))
    }
}
