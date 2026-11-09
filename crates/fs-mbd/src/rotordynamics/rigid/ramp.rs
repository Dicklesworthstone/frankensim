//! Prescribed angular acceleration in the original stationary rotor coordinates.
//!
//! With `theta = [theta_x, theta_y]` and shaft direction
//! `n = [theta_y, -theta_x, 1]`, transverse angular momentum is
//! `H_perp = Jd theta_dot + Jp Omega [theta_y, -theta_x]`. Differentiating
//! adds BOTH `Omega G theta_dot` and `alpha G theta`. A drive torque along
//! the tilted shaft supplies the latter transverse moment and cancels it;
//! a drive along the nominal stationary +z axis does not. This physical
//! boundary choice must be explicit, including at zero spin and reversal.
//!
//! The same native full-matrix midpoint solve advances the lateral dynamics.
//! Speed is affine and phase quadratic within each prescribed interval.
//! The source is sampled at its actual temporal midpoint, including tangential
//! unbalance. This is second order under temporal refinement; it is not a
//! finite-torque motor, spin-energy reservoir or flexible-shaft model.

use super::{PreparedRotorUnbalance, RigidRotor, RotorError, finite, poll};
use crate::rotordynamics::transient::{RotorTransientState, RotorTransientStep};
use fs_math::dd::Dd;
use std::borrow::Cow;

/// Direction of the torque prescribing the nominal rotor's spin acceleration.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RotorSpinDriveAxis {
    /// Torque follows the actual tilted shaft. Its transverse component
    /// cancels `alpha G q`, leaving `M q'' + (C + Omega G) q' + K q = f`.
    Shaft,
    /// Torque remains parallel to nominal stationary +z. The full equation
    /// is `M q'' + (C + Omega G) q' + (K + alpha G) q = f`.
    Stationary,
}

/// Continuous physical clock and common unbalance phase for one rotor.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RotorSpinState {
    /// Physical time [s]. An accepted step must advance this represented clock.
    pub time_seconds: f64,
    /// Common shaft phase [rad], added to each prepared unbalance phase offset.
    /// Retained continuously, rather than recomputed as current speed times time.
    pub phase_rad: f64,
    /// Signed shaft speed about nominal +z [rad/s].
    pub spin_rad_s: f64,
}

/// One constant-acceleration interval and its explicit numerical admission.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RotorSpinRamp {
    /// Signed prescribed angular acceleration [rad/s²]. Braking and reversal
    /// need no special handling or division by shaft speed.
    pub angular_acceleration_rad_s2: f64,
    /// Physical interval duration [s], finite and strictly positive.
    pub duration_seconds: f64,
    /// Physical torque direction; no direction is silently assumed.
    pub drive_axis: RotorSpinDriveAxis,
    /// Relative original momentum residual allowed in `(0,1)`.
    pub relative_tolerance: f64,
}

/// Accepted lateral motion and exactly integrated prescribed spin kinematics.
#[derive(Clone, Debug, PartialEq)]
pub struct RigidRotorSpinRampStep {
    /// Native full-coordinate state, momentum residual and independent work.
    /// Its stiffness work includes the acceleration coupling when present;
    /// the two separately reported work terms below partition that work.
    pub motion: RotorTransientStep,
    /// Reusable physical time, phase and signed spin at the accepted endpoint.
    pub spin_after: RotorSpinState,
    /// Actual temporal midpoint used to sample the time-dependent operators
    /// and unbalance, including quadratic rather than linearly averaged phase.
    pub spin_midpoint: RotorSpinState,
    /// Generalized unbalance `[Fx,Fy,tau_x,tau_y]` in `[N,N,Nm,Nm]`, including
    /// the acceleration reaction. This is part of `motion.input_work_j`.
    pub unbalance_force: [f64; 4],
    /// Signed work against the original bearing stiffness `K q` [J].
    pub bearing_stiffness_work_j: f64,
    /// Signed work against `alpha G q` [J] for stationary drive torque,
    /// zero for shaft-following torque. This is energy exchange with the
    /// prescribed spin boundary, not heat or stored elastic energy.
    pub acceleration_coupling_work_j: f64,
}

impl RigidRotor {
    /// Advance lateral motion through one prescribed physical spin ramp.
    ///
    /// The ORIGINAL full `M,C,K,G` are retained, including active and
    /// cross-coupled bearings. A four-coordinate midpoint factorization is
    /// prepared at the actual interval midpoint. Unbalance preparation is
    /// reusable; this call neither traverses its source stations nor computes
    /// modes. Constant acceleration zero reduces to the steady native step.
    ///
    /// `midpoint_force` is an additional force/torque in the same stationary
    /// coordinates, sampled at `time_seconds + duration_seconds/2`. Both its
    /// work and unbalance work enter the native total input work. Bearing
    /// coefficients stay those authored in this rotor; reassemble when their
    /// physical operating conditions change. Nominal small-tilt and small-
    /// unbalance assumptions still apply throughout the resolved trajectory.
    ///
    /// Time, speed, phase and physical operators must remain representable.
    /// No state or partial spin advance is exposed after input, solve,
    /// original momentum-balance or cancellation refusal.
    pub fn spin_ramp_step(
        &self,
        state: &RotorTransientState,
        spin: RotorSpinState,
        unbalance: &PreparedRotorUnbalance,
        midpoint_force: [f64; 4],
        ramp: RotorSpinRamp,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<RigidRotorSpinRampStep, RotorError> {
        poll(&mut cancelled)?;
        let h = ramp.duration_seconds;
        let half = 0.5 * h;
        let alpha = finite(ramp.angular_acceleration_rad_s2)?;
        if !h.is_finite() || half <= 0.0 {
            return Err(RotorError::Invalid(
                "positive representable spin-ramp duration required",
            ));
        }
        finite(spin.time_seconds)?;
        finite(spin.phase_rad)?;
        finite(spin.spin_rad_s)?;
        let time_after = finite(spin.time_seconds + h)?;
        if time_after <= spin.time_seconds {
            return Err(RotorError::Invalid("spin-ramp clock must advance"));
        }
        let midpoint_speed = finite(alpha.mul_add(half, spin.spin_rad_s))?;
        let spin_midpoint = RotorSpinState {
            time_seconds: finite(spin.time_seconds + half)?,
            phase_rad: finite(half.mul_add(
                finite(alpha.mul_add(0.5 * half, spin.spin_rad_s))?,
                spin.phase_rad,
            ))?,
            spin_rad_s: midpoint_speed,
        };
        let spin_after = RotorSpinState {
            time_seconds: time_after,
            phase_rad: finite(h.mul_add(midpoint_speed, spin.phase_rad))?,
            spin_rad_s: finite(alpha.mul_add(h, spin.spin_rad_s))?,
        };
        let unbalance_force = unbalance.force_at(
            spin_midpoint.phase_rad,
            spin_midpoint.spin_rad_s,
            alpha,
            &mut cancelled,
        )?;
        let mut force = [0.0; 4];
        for i in 0..4 {
            force[i] = finite(finite(midpoint_force[i])? + unbalance_force[i])?;
        }
        let operating_point = self.model.at_speed(midpoint_speed, &mut cancelled)?;
        let has_euler_coupling = ramp.drive_axis == RotorSpinDriveAxis::Stationary && alpha != 0.0;
        let stepper = if has_euler_coupling {
            let stiffness = self
                .model
                .stiffness
                .iter()
                .zip(&self.model.gyro)
                .map(|(&k, &g)| finite(alpha.mul_add(g, k)))
                .collect::<Result<Vec<_>, _>>()?;
            operating_point.midpoint_stepper_with_stiffness(
                h,
                ramp.relative_tolerance,
                Cow::Owned(stiffness),
                &mut cancelled,
            )?
        } else {
            operating_point.midpoint_stepper(h, ramp.relative_tolerance, &mut cancelled)?
        };
        let motion = stepper.step(state, &force, &mut cancelled)?;
        let (bearing_stiffness_work_j, acceleration_coupling_work_j) = if has_euler_coupling {
            let mut bearing_work = Dd::ZERO;
            let mut acceleration_work = Dd::ZERO;
            for i in 0..4 {
                poll(&mut cancelled)?;
                let v = 0.5 * state.velocity[i] + 0.5 * motion.state_after.velocity[i];
                let mut bearing_force = 0.0;
                let mut acceleration_force = 0.0;
                for j in 0..4 {
                    let q = 0.5 * state.displacement[j] + 0.5 * motion.state_after.displacement[j];
                    bearing_force = self.model.stiffness[4 * i + j].mul_add(q, bearing_force);
                    acceleration_force = self.model.gyro[4 * i + j].mul_add(q, acceleration_force);
                }
                let bearing_impulse = finite(h * finite(bearing_force)?)?;
                let acceleration_impulse =
                    finite(h * finite(alpha * finite(acceleration_force)?)?)?;
                bearing_work = bearing_work + Dd::from_f64(bearing_impulse) * Dd::from_f64(v);
                acceleration_work =
                    acceleration_work + Dd::from_f64(acceleration_impulse) * Dd::from_f64(v);
            }
            (
                finite(bearing_work.to_f64())?,
                finite(acceleration_work.to_f64())?,
            )
        } else {
            (motion.stiffness_force_work_j, 0.0)
        };
        poll(&mut cancelled)?;
        Ok(RigidRotorSpinRampStep {
            motion,
            spin_after,
            spin_midpoint,
            unbalance_force,
            bearing_stiffness_work_j,
            acceleration_coupling_work_j,
        })
    }
}
