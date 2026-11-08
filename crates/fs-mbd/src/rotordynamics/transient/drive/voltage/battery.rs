//! Finite electrical supply coupled to the ORIGINAL full-coordinate rotor.
//!
//! The winding/rotor midpoint equations are condensed to their affine terminal
//! current and inserted in the battery/converter's native nonlinear equation.
//! Accepted voltage then recovers the original mechanical motion. There is one
//! mechanical inertia model, no sampled bus voltage, no fixed-point interface
//! loop, and no postsolve correction of current, motion or energy.
//!
//! The converter owns conversion losses. The winding connection therefore
//! requires unit motoring/regenerative bridge efficiencies; a second empirical
//! bridge is not silently inserted. A limit requiring voltage curtailment
//! refuses the whole step. These are the existing averaged/discrete-gradient
//! electrical and constant-operating-point mechanical equations: no switching,
//! protection event, thermal drift, continuous peak or stability claim is added.
//! Motor power caps apply to its held-discrete-voltage/current reconstruction,
//! not to upstream battery chemical power or unresolved bus-voltage peaks.

use super::{
    Dd, RotorVoltageControl, RotorVoltageError, RotorVoltageMotor, RotorVoltageStep,
    RotorWindingState, poll,
};
use fs_phs::PreparedStepError;
use fs_phs::circuit::converter::battery_fed::{
    AffineOutputCurrent, BatteryFedBuckBoost, BatteryFedStep, PreparedBatteryFedBuckBoost,
};
use fs_phs::circuit::converter::power_load::{PowerLoadControls, RegulatedPowerLoad};

fn sum(values: impl IntoIterator<Item = f64>) -> f64 {
    values
        .into_iter()
        .fold(Dd::ZERO, |total, value| total + Dd::from_f64(value))
        .to_f64()
}

/// Continuing source charges/flux, winding current and original rotor state.
#[derive(Clone, Debug, PartialEq)]
pub struct BatteryRotorState {
    /// Original six-state `BatteryFedBuckBoost` ordering, without copied inertia.
    pub electrical: [f64; 6],
    pub motor: RotorWindingState,
}

/// One coupled interval retaining each original owner's physical accounts.
#[derive(Clone, Debug)]
pub struct BatteryRotorStep {
    pub electrical_after: [f64; 6],
    /// Original six-state solver with the condensed terminal in its equation.
    /// Its output work is INTERNAL to this complete electromechanical system.
    pub electrical: BatteryFedStep,
    /// Full native motion and its signed C/K/gyroscopic work, plus winding heat.
    /// Active damping/nonconservative stiffness are not relabeled as heat.
    pub motor: RotorVoltageStep,
    /// Motor midpoint current plus current entering the electrical output [A].
    pub terminal_current_residual_a: f64,
    /// Signed source output work plus winding terminal input work [J].
    pub terminal_work_residual_j: f64,
    /// Electrical/magnetic/kinetic changes, original mechanical operator work,
    /// electrical/winding heat and auxiliary export minus external load work.
    pub balance_residual_j: f64,
}

impl BatteryRotorStep {
    /// Original eight electrical channels followed by winding Joule heat.
    /// Mechanical C/K work remains in `motor.motion`, including its sign.
    pub fn heat_joules(&self) -> [f64; 9] {
        let electrical = self.electrical.heat_joules();
        core::array::from_fn(|i| {
            if i == 8 {
                self.motor.winding_heat_j
            } else {
                electrical[i]
            }
        })
    }

    pub fn into_state(self) -> BatteryRotorState {
        BatteryRotorState {
            electrical: self.electrical_after,
            motor: self.motor.into_state(),
        }
    }
}

/// Atomic refusal: no source debit or mechanical state is published.
#[derive(Debug)]
pub enum BatteryRotorError {
    Electrical(PreparedStepError),
    Motor(RotorVoltageError),
    Invalid(&'static str),
    /// The original motor's limits would change the solved bus voltage.
    VoltageCurtailment {
        bus_voltage_v: f64,
        applied_voltage_v: f64,
    },
    CurrentMismatch {
        residual_a: f64,
        tolerance_a: f64,
    },
    EnergyBalance {
        boundary: &'static str,
        residual_j: f64,
        tolerance_j: f64,
    },
}

impl From<PreparedStepError> for BatteryRotorError {
    fn from(error: PreparedStepError) -> Self {
        Self::Electrical(error)
    }
}
impl From<RotorVoltageError> for BatteryRotorError {
    fn from(error: RotorVoltageError) -> Self {
        Self::Motor(error)
    }
}
impl std::fmt::Display for BatteryRotorError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Electrical(error) => write!(f, "{error}"),
            Self::Motor(error) => write!(f, "{error}"),
            Self::Invalid(reason) => write!(f, "invalid battery rotor drive: {reason}"),
            Self::VoltageCurtailment {
                bus_voltage_v,
                applied_voltage_v,
            } => write!(
                f,
                "rotor limits change bus {bus_voltage_v} V to {applied_voltage_v} V"
            ),
            Self::CurrentMismatch {
                residual_a,
                tolerance_a,
            } => {
                write!(
                    f,
                    "rotor terminal current mismatch {residual_a} A exceeds {tolerance_a} A"
                )
            }
            Self::EnergyBalance {
                boundary,
                residual_j,
                tolerance_j,
            } => {
                write!(
                    f,
                    "{boundary} residual {residual_j} J exceeds {tolerance_j} J"
                )
            }
        }
    }
}
impl std::error::Error for BatteryRotorError {}

/// Cached original mechanical factors and one original electrical workspace.
pub struct BatteryRotorDrive<'step, 'model> {
    motor: RotorVoltageMotor<'step, 'model>,
    electrical: PreparedBatteryFedBuckBoost,
    electrical_controls: PowerLoadControls,
    motor_controls: RotorVoltageControl,
}

impl<'step, 'model> RotorVoltageMotor<'step, 'model> {
    /// Connect a finite battery/converter directly to this winding terminal.
    ///
    /// Reuses this prepared duration and mechanical factors. Source preparation
    /// allocates its existing six-state workspace; no mechanical states or
    /// inertia are appended there. Current agreement uses the supplied current
    /// allowance; terminal and full balances use the electrical energy controls.
    /// Winding-terminal limits and endpoint current capabilities retain the
    /// original motor rules. Auxiliary demand and conversion loss are additional
    /// battery demand, not included in that winding-terminal power capability.
    pub fn with_battery_converter(
        self,
        electrical: BatteryFedBuckBoost,
        electrical_controls: PowerLoadControls,
        motor_controls: RotorVoltageControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BatteryRotorDrive<'step, 'model>, BatteryRotorError> {
        poll(&mut cancelled).map_err(RotorVoltageError::from)?;
        if self.motor.drive.motoring_efficiency != 1.0
            || self.motor.drive.regenerative_efficiency != 1.0
        {
            return Err(BatteryRotorError::Invalid(
                "direct winding connection requires unit bridge efficiencies; converter owns losses",
            ));
        }
        self.motor
            .drive
            .admit(self.stepper.model.n, motor_controls.source, 0.0)
            .map_err(RotorVoltageError::from)?;
        if !motor_controls.absolute_current_tolerance_a.is_finite()
            || motor_controls.absolute_current_tolerance_a < 0.0
        {
            return Err(BatteryRotorError::Invalid(
                "finite nonnegative current allowance required",
            ));
        }
        let electrical = electrical.prepare_controlled(
            self.duration_s(),
            electrical_controls,
            &mut cancelled,
        )?;
        Ok(BatteryRotorDrive {
            motor: self,
            electrical,
            electrical_controls,
            motor_controls,
        })
    }
}

impl BatteryRotorDrive<'_, '_> {
    /// Immutable access to original source parameters and stored-energy model.
    pub fn electrical(&self) -> &PreparedBatteryFedBuckBoost {
        &self.electrical
    }

    fn admit_energy(
        &self,
        boundary: &'static str,
        residual_j: f64,
        scale_j: f64,
    ) -> Result<(), BatteryRotorError> {
        let tolerance_j = self.electrical_controls.absolute_energy_tolerance_j
            + self.electrical_controls.relative_energy_tolerance * scale_j;
        if !scale_j.is_finite()
            || !residual_j.is_finite()
            || !tolerance_j.is_finite()
            || residual_j.abs() > tolerance_j
        {
            return Err(BatteryRotorError::EnergyBalance {
                boundary,
                residual_j,
                tolerance_j,
            });
        }
        Ok(())
    }

    /// Advance both owners with feedback through their exact discrete terminal.
    /// External generalized loads remain uncurtailed. Optional auxiliary demand
    /// shares the same finite output bus. Failure/cancellation leaves the input
    /// state untouched; the result is published only after all original and
    /// combined balances and limits accept. Preparation remains usable on retry.
    pub fn step_controlled(
        &mut self,
        initial: &BatteryRotorState,
        external_force: &[f64],
        auxiliary_load: Option<RegulatedPowerLoad>,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<BatteryRotorStep, BatteryRotorError> {
        poll(&mut cancelled).map_err(RotorVoltageError::from)?;
        let (constant, slope) =
            self.motor
                .midpoint_current_response(&initial.motor, external_force, &mut cancelled)?;
        let terminal = AffineOutputCurrent {
            constant_a: -constant,
            slope_a_per_v: -slope,
        };
        let mut electrical_after = [0.0; 6];
        let electrical = self.electrical.step_affine_output_current_controlled(
            &initial.electrical,
            terminal,
            auxiliary_load,
            &mut electrical_after,
            &mut cancelled,
        )?;
        let voltage_v = electrical.discrete_bus_voltage_v[2];
        let motor = self.motor.step(
            &initial.motor,
            external_force,
            voltage_v,
            self.motor_controls,
            &mut cancelled,
        )?;
        if motor.applied_voltage_v != voltage_v {
            return Err(BatteryRotorError::VoltageCurtailment {
                bus_voltage_v: voltage_v,
                applied_voltage_v: motor.applied_voltage_v,
            });
        }
        let terminal_current_residual_a = sum([
            motor.midpoint_current_a,
            terminal
                .current_into_bus_a(voltage_v)
                .map_err(PreparedStepError::from)?,
        ]);
        if !terminal_current_residual_a.is_finite()
            || terminal_current_residual_a.abs() > self.motor_controls.absolute_current_tolerance_a
        {
            return Err(BatteryRotorError::CurrentMismatch {
                residual_a: terminal_current_residual_a,
                tolerance_a: self.motor_controls.absolute_current_tolerance_a,
            });
        }
        let terminal_work_residual_j = sum([
            electrical.external_work_j,
            motor.electrical_draw_j,
            -motor.electrical_return_j,
        ]);
        self.admit_energy(
            "electrical/winding terminal",
            terminal_work_residual_j,
            sum([
                electrical.external_work_j.abs(),
                motor.electrical_draw_j,
                motor.electrical_return_j,
            ]),
        )?;
        let terms = [
            electrical.interval().storage_change_j(),
            motor.magnetic_energy_change_j,
            motor.motion.kinetic_energy_change_j,
            motor.motion.velocity_force_work_j,
            motor.motion.stiffness_force_work_j,
            motor.motion.gyroscopic_work_j,
            sum(electrical.heat_joules()),
            motor.winding_heat_j,
            electrical.load_work_j,
            -motor.external_work_j,
        ];
        let balance_residual_j = sum(terms);
        self.admit_energy(
            "battery/converter/winding/rotor",
            balance_residual_j,
            sum(terms.into_iter().map(f64::abs)),
        )?;
        poll(&mut cancelled).map_err(RotorVoltageError::from)?;
        Ok(BatteryRotorStep {
            electrical_after,
            electrical,
            motor,
            terminal_current_residual_a,
            terminal_work_residual_j,
            balance_residual_j,
        })
    }
}
