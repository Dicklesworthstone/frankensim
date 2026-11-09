//! Finite electrical spin actuation of a balanced axisymmetric rigid rotor.
//!
//! The admitted rotor supplies its ORIGINAL polar inertia to the native
//! current/voltage/power-limited shaft owner. Its accepted midpoint speed then
//! drives the original four-coordinate lateral gyroscopic equation. Shaft-axis
//! torque cancels the transverse Euler term, as derived in the prescribed-ramp
//! owner. Neither a tilt coordinate nor a caller-invented inertia is used as
//! the shaft. Gyroscopic forces do no lateral work; axial motor work changes
//! spin storage and supplies separately authored shaft drag and shaft loads.
//!
//! This retains the existing small-tilt truncation and constant bearing
//! coefficients. Residual unbalance is deliberately absent: its reciprocal
//! spin reaction requires a coupled eccentric-inertia model. Current is held
//! during a step; finite inductance, motor saturation and thermal evolution
//! remain the corresponding electrical/material owners' responsibilities.

use super::{
    RigidRotor, RotorBearingModel, RotorError, RotorSpinDriveAxis, RotorSpinState, finite, poll,
};
use crate::rotordynamics::transient::drive::{
    RotorCurrentDrive, RotorCurrentDriveStep, RotorDriveError, RotorDrivePowerControl,
};
use crate::rotordynamics::transient::{RotorTransientState, RotorTransientStep};
use fs_math::dd::Dd;

/// Original balanced rotor and its reciprocal finite-current spin motor.
#[derive(Clone, Debug)]
pub struct RigidRotorCurrentMotor<'a> {
    rotor: &'a RigidRotor,
    shaft: RotorBearingModel,
    motor: RotorCurrentDrive,
}

/// Electrical command and externally supplied mechanical loads for one step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RotorSpinCurrentCommand {
    /// Held winding-current command [A], subject to the native drive limits.
    pub current_a: f64,
    /// Signed external torque along the shaft [N m], never curtailed with current.
    pub shaft_load_torque_nm: f64,
    /// Stationary `[Fx,Fy,tau_x,tau_y]` midpoint loads in `[N,N,Nm,Nm]`.
    /// These are independent external work ports, not motor-powered unbalance.
    pub lateral_force: [f64; 4],
}

/// Shared native duration, momentum admission and electrical capability.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RotorSpinDriveControl {
    pub duration_seconds: f64,
    pub relative_tolerance: f64,
    /// The power allowance also admits the combined energy residual per second.
    pub electrical: RotorDrivePowerControl,
}

/// Complete accepted spin/lateral motion and its original electrical accounts.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidRotorCurrentStep {
    /// Native shaft motion, torque, back-EMF, source draw/return and separate
    /// winding/bridge heat. Its kinetic change is the physical spin storage;
    /// its velocity-force work is the nonnegative authored shaft-drag loss.
    pub motor: RotorCurrentDriveStep,
    /// Original four-coordinate motion and signed bearing/gyroscopic work.
    pub lateral: RotorTransientStep,
    pub spin_after: RotorSpinState,
    /// Accepted midpoint spin used in the original lateral gyroscopic operator.
    pub midpoint_spin_rad_s: f64,
    /// Actual native endpoint speed increment divided by duration [rad/s²].
    pub angular_acceleration_rad_s2: f64,
    /// Spin/lateral kinetic changes plus operator work and electrical heat,
    /// minus net electrical input and independently supplied external work [J].
    pub balance_residual_j: f64,
}

/// Atomic refusal: no spin, lateral motion or electrical debit is published.
#[derive(Clone, Debug, PartialEq)]
pub enum RotorSpinDriveError {
    Native(RotorDriveError),
    Invalid(&'static str),
    EnergyBalance { residual_j: f64, tolerance_j: f64 },
}

impl From<RotorDriveError> for RotorSpinDriveError {
    fn from(error: RotorDriveError) -> Self {
        Self::Native(error)
    }
}
impl From<RotorError> for RotorSpinDriveError {
    fn from(error: RotorError) -> Self {
        Self::Native(error.into())
    }
}
impl core::fmt::Display for RotorSpinDriveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "rigid rotor spin drive refused: {self:?}")
    }
}
impl std::error::Error for RotorSpinDriveError {}

impl RigidRotor {
    /// Bind the physical polar inertia to a finite-current spin motor.
    ///
    /// `drive_axis` must explicitly be `Shaft`. The motor's `shaft_coordinate`
    /// must be zero in its own one-coordinate spin model; it never indexes the
    /// four lateral coordinates. Nonnegative `shaft_drag_nm_s` supplies a
    /// passive torque `-drag*Omega`. All remaining motor parameters and actual
    /// operating limits are admitted by the original drive at each step.
    pub fn current_spin_motor(
        &self,
        motor: RotorCurrentDrive,
        drive_axis: RotorSpinDriveAxis,
        shaft_drag_nm_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RigidRotorCurrentMotor<'_>, RotorSpinDriveError> {
        poll(&mut cancelled)?;
        if drive_axis != RotorSpinDriveAxis::Shaft || motor.shaft_coordinate != 0 {
            return Err(RotorSpinDriveError::Invalid(
                "shaft-axis drive and spin coordinate zero required",
            ));
        }
        if !shaft_drag_nm_s.is_finite() || shaft_drag_nm_s < 0.0 {
            return Err(RotorSpinDriveError::Invalid(
                "finite nonnegative shaft drag required",
            ));
        }
        // The rigid assembly already lowered and admitted this original Jp.
        // G[theta_x,theta_y]=Jp; lateral M contains Jd, not this spin storage.
        let shaft = RotorBearingModel::new(
            1,
            &[self.model.gyro[11]],
            &[shaft_drag_nm_s],
            &[0.0],
            &[0.0],
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        Ok(RigidRotorCurrentMotor {
            rotor: self,
            shaft,
            motor,
        })
    }
}

impl RigidRotorCurrentMotor<'_> {
    /// Advance real spin and lateral motion under finite electrical capability.
    ///
    /// One native scalar shaft preparation and the original four-coordinate
    /// preparation suffice; no eigenproblem, inverse or new integration kernel
    /// is introduced. The changing gyroscopic coefficient is sampled using
    /// the accepted shaft endpoints, including stall and reversal. For smooth
    /// resolved forcing this composition is second order. Sharp controller or
    /// limit transitions require the caller's event/time refinement.
    ///
    /// Spin phase comes directly from the native shaft displacement update;
    /// it is not reconstructed as current speed times absolute time. All
    /// current, voltage, speed, source-power and native physical-work gates
    /// apply before the complete energy balance and final cancellation poll.
    pub fn step(
        &self,
        lateral_state: &RotorTransientState,
        spin: RotorSpinState,
        command: RotorSpinCurrentCommand,
        control: RotorSpinDriveControl,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RigidRotorCurrentStep, RotorSpinDriveError> {
        poll(&mut cancelled)?;
        let h = control.duration_seconds;
        finite(spin.time_seconds)?;
        let time_after = finite(spin.time_seconds + h)?;
        if time_after <= spin.time_seconds {
            return Err(RotorSpinDriveError::Invalid(
                "spin-drive clock must advance",
            ));
        }
        let shaft_state = RotorTransientState {
            displacement: vec![spin.phase_rad],
            velocity: vec![spin.spin_rad_s],
        };
        let motor = self
            .shaft
            .at_speed(0.0, &mut cancelled)?
            .midpoint_stepper(h, control.relative_tolerance, &mut cancelled)?
            .step_current_drive(
                &shaft_state,
                &[command.shaft_load_torque_nm],
                command.current_a,
                self.motor,
                control.electrical,
                &mut cancelled,
            )?;
        let speed_after = motor.motion.state_after.velocity[0];
        let midpoint_spin_rad_s = finite(0.5 * spin.spin_rad_s + 0.5 * speed_after)?;
        let angular_acceleration_rad_s2 = finite((speed_after - spin.spin_rad_s) / h)?;
        let lateral = self
            .rotor
            .model
            .at_speed(midpoint_spin_rad_s, &mut cancelled)?
            .midpoint_stepper(h, control.relative_tolerance, &mut cancelled)?
            .step(lateral_state, &command.lateral_force, &mut cancelled)?;
        let balance_residual_j = finite(
            [
                motor.motion.kinetic_energy_change_j,
                motor.motion.velocity_force_work_j,
                motor.motion.stiffness_force_work_j,
                motor.motion.gyroscopic_work_j,
                lateral.kinetic_energy_change_j,
                lateral.velocity_force_work_j,
                lateral.stiffness_force_work_j,
                lateral.gyroscopic_work_j,
                motor.winding_heat_j,
                motor.bridge_heat_j,
                motor.electrical_return_j,
                -motor.electrical_draw_j,
                -motor.external_work_j,
                -lateral.input_work_j,
            ]
            .into_iter()
            .fold(Dd::ZERO, |total, term| total + Dd::from_f64(term))
            .to_f64(),
        )?;
        let tolerance_j = finite(h * control.electrical.absolute_power_tolerance_w)?;
        if balance_residual_j.abs() > tolerance_j {
            return Err(RotorSpinDriveError::EnergyBalance {
                residual_j: balance_residual_j,
                tolerance_j,
            });
        }
        let spin_after = RotorSpinState {
            time_seconds: time_after,
            phase_rad: motor.motion.state_after.displacement[0],
            spin_rad_s: speed_after,
        };
        poll(&mut cancelled)?;
        Ok(RigidRotorCurrentStep {
            motor,
            lateral,
            spin_after,
            midpoint_spin_rad_s,
            angular_acceleration_rad_s2,
            balance_residual_j,
        })
    }
}
