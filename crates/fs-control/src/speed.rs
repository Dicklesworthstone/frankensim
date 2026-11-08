//! Sampled relative-shaft speed control around a machine-owned torque loop.
//!
//! Electromagnetic torque acts equally and oppositely on rotor and stator;
//! its mechanical power is therefore `torque * (omega_rotor-omega_stator)`.
//! Neither absolute rotor speed nor an assumed fixed stator defines this
//! controller's speed error or power envelope. Motoring and regeneration limits
//! follow that signed power in either direction of rotation.
//!
//! The digital PI law predicts `z_p = z + Ts*Ki*error`, requests
//! `u = feedforward + Kp*error + z_p`, then applies the authored torque/power
//! envelope. Pass that request to the ORIGINAL machine's torque/current
//! allocator. `SpeedDemand::finish_controlled` tracks its achievable torque:
//! `z_next = z_p + Ts/(Taw+Ts)*(achievable-u)`. A zero tracking time requests
//! immediate tracking. The allocator, not nominal torque/current constants in
//! this controller, retains authority over saturation and magnetic support.
//!
//! Finish the outer proposal and inner current proposal before advancing the
//! plant. Adopt both proposed states only after that held plant interval is
//! accepted. No controller proposal injects work or changes physical storage.
//! This fixed-period law allocates nothing and owns no scheduling or plant integrator.
//!
//! Limits constrain REQUESTED torque and power at the sampled speed. At rest,
//! mechanical power is zero and finite starting torque remains available; no
//! arbitrary speed floor divides the power limits. Future speed, current-loop
//! dynamics, bus sag and stored magnetic energy can change actual interval
//! power. These limits are not electrical overcurrent, battery-charge protection
//! or a bound on unresolved transient power. Allocator anti-windup does not by
//! itself detect a voltage-infeasible current request. Gains and limits are
//! authored; no stability, thermal, field-weakening or hardware claim is made.

use crate::ControlError;

/// Authored scalar PI data in the torque/relative-speed coordinate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpeedControllerParameters {
    pub sample_period_s: f64,
    /// Proportional gain [N m s/rad], finite and nonnegative.
    pub proportional_gain_nm_s_per_rad: f64,
    /// Integral gain [N m/rad], finite and nonnegative.
    pub integral_gain_nm_per_rad: f64,
    /// Tracking time [s], nonnegative; zero means immediate anti-windup.
    pub anti_windup_time_s: f64,
}

/// Caller-owned accepted PI torque contribution. This is not physical energy.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SpeedControllerState {
    pub integral_torque_nm: f64,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpeedCommand {
    /// Requested rotor speed relative to the stator [rad/s], either sign.
    pub relative_speed_rad_s: f64,
    /// Authored signed electromagnetic-torque feedforward [N m].
    pub feedforward_torque_nm: f64,
}

/// Runtime command capability. All four magnitudes are finite/nonnegative;
/// zero removes that capability. These may reflect an externally owned thermal
/// or source policy, but this controller does not infer such a policy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ShaftPowerLimits {
    pub positive_torque_limit_nm: f64,
    pub negative_torque_limit_nm: f64,
    pub motoring_power_w: f64,
    pub regeneration_power_w: f64,
}

/// Torque request awaiting the original machine's allocation. Finishing this
/// demand creates a candidate state; it does not accept a controller sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpeedDemand {
    pub measured_relative_speed_rad_s: f64,
    pub speed_error_rad_s: f64,
    /// PI and feedforward before any torque or power limit [N m].
    pub unlimited_torque_nm: f64,
    /// Pass this signed request to the machine-owned torque allocator [N m].
    pub requested_torque_nm: f64,
    /// Feasible negative/positive request endpoints at this sampled speed.
    pub torque_bounds_nm: [f64; 2],
    pub requested_shaft_power_w: f64,
    pub torque_limited: bool,
    pub power_limited: bool,
    shaft_power_bounds_w: [f64; 2],
    predicted_integral_torque_nm: f64,
    tracking_fraction: f64,
}

/// Candidate outer state to accept only with the matching current-loop/plant
/// interval. Achievable torque is a constitutive command, not measured future
/// torque or a substitute for the native work and heat records.
#[derive(Clone, Copy, Debug, PartialEq)]
#[must_use = "accept next_state only with the corresponding successful plant interval"]
pub struct SpeedControlProposal {
    pub next_state: SpeedControllerState,
    pub achievable_torque_nm: f64,
    pub achievable_shaft_power_w: f64,
    /// Includes any representable allocation difference, including roundoff.
    pub allocation_limited: bool,
}

/// Prepared fixed-work PI coefficients, independent of any particular plant.
#[derive(Clone, Debug)]
pub struct SpeedController {
    parameters: SpeedControllerParameters,
    integral_increment_gain_nm_s_per_rad: f64,
    tracking_fraction: f64,
}

fn finite(values: &[f64], what: &'static str) -> Result<(), ControlError> {
    if values.iter().all(|v| v.is_finite()) {
        Ok(())
    } else {
        Err(ControlError::NumericalRange { what })
    }
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), ControlError> {
    if cancelled() {
        Err(ControlError::Cancelled)
    } else {
        Ok(())
    }
}

// A division overflow means the power cap is above every representable torque
// and the finite authored torque bound remains active. At exactly zero speed
// no division is made. Positive caps below torque resolution conservatively
// admit zero torque; no epsilon speed or invented positive power is introduced.
fn power_bound(torque: f64, power: f64, speed: f64) -> f64 {
    if speed == 0.0 {
        torque
    } else {
        torque.min(power / speed)
    }
}

impl SpeedController {
    pub fn new(parameters: SpeedControllerParameters) -> Result<Self, ControlError> {
        if !parameters.sample_period_s.is_finite() || parameters.sample_period_s <= 0.0 {
            return Err(ControlError::InvalidParameter {
                what: "speed sample period must be finite and positive",
            });
        }
        if [
            parameters.proportional_gain_nm_s_per_rad,
            parameters.integral_gain_nm_per_rad,
            parameters.anti_windup_time_s,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(ControlError::InvalidParameter {
                what: "speed PI gains and tracking time must be finite nonnegative",
            });
        }
        let ts = parameters.sample_period_s;
        let tau = parameters.anti_windup_time_s;
        let integral_increment_gain_nm_s_per_rad = ts * parameters.integral_gain_nm_per_rad;
        let tracking_fraction = if tau <= ts {
            1.0 / (1.0 + tau / ts)
        } else {
            let r = ts / tau;
            r / (1.0 + r)
        };
        finite(
            &[integral_increment_gain_nm_s_per_rad, tracking_fraction],
            "sampled speed gains",
        )?;
        if (parameters.integral_gain_nm_per_rad > 0.0
            && integral_increment_gain_nm_s_per_rad == 0.0)
            || tracking_fraction == 0.0
        {
            return Err(ControlError::NumericalRange {
                what: "sampled speed gains underflow",
            });
        }
        Ok(Self {
            parameters,
            integral_increment_gain_nm_s_per_rad,
            tracking_fraction,
        })
    }

    pub fn parameters(&self) -> SpeedControllerParameters {
        self.parameters
    }

    /// Compute one sampled request without changing any accepted state. Both
    /// speeds are absolute; equal common rotation changes neither error nor
    /// request. Hold the resulting torque/current command for exactly Ts.
    pub fn demand_controlled(
        &self,
        state: &SpeedControllerState,
        shaft_speed_rad_s: [f64; 2],
        command: SpeedCommand,
        limits: ShaftPowerLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SpeedDemand, ControlError> {
        poll(&mut cancelled)?;
        if [
            state.integral_torque_nm,
            shaft_speed_rad_s[0],
            shaft_speed_rad_s[1],
            command.relative_speed_rad_s,
            command.feedforward_torque_nm,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(ControlError::InvalidSample {
                what: "speed state, command and both shafts must be finite",
            });
        }
        if [
            limits.positive_torque_limit_nm,
            limits.negative_torque_limit_nm,
            limits.motoring_power_w,
            limits.regeneration_power_w,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(ControlError::InvalidSample {
                what: "torque and power capabilities must be finite nonnegative",
            });
        }
        let speed = shaft_speed_rad_s[0] - shaft_speed_rad_s[1];
        let error = command.relative_speed_rad_s - speed;
        let predicted =
            state.integral_torque_nm + self.integral_increment_gain_nm_s_per_rad * error;
        let unlimited = command.feedforward_torque_nm
            + self.parameters.proportional_gain_nm_s_per_rad * error
            + predicted;
        finite(
            &[speed, error, predicted, unlimited],
            "speed error and PI torque",
        )?;
        let (positive_power, negative_power) = if speed >= 0.0 {
            (limits.motoring_power_w, limits.regeneration_power_w)
        } else {
            (limits.regeneration_power_w, limits.motoring_power_w)
        };
        let bounds = [
            -power_bound(limits.negative_torque_limit_nm, negative_power, speed.abs()),
            power_bound(limits.positive_torque_limit_nm, positive_power, speed.abs()),
        ];
        let torque_clipped = unlimited.clamp(
            -limits.negative_torque_limit_nm,
            limits.positive_torque_limit_nm,
        );
        let requested = torque_clipped.clamp(bounds[0], bounds[1]);
        let power = requested * speed;
        finite(&[power], "requested sampled shaft power")?;
        let result = SpeedDemand {
            measured_relative_speed_rad_s: speed,
            speed_error_rad_s: error,
            unlimited_torque_nm: unlimited,
            requested_torque_nm: requested,
            torque_bounds_nm: bounds,
            requested_shaft_power_w: power,
            torque_limited: torque_clipped != unlimited,
            power_limited: requested != torque_clipped,
            shaft_power_bounds_w: [-limits.regeneration_power_w, limits.motoring_power_w],
            predicted_integral_torque_nm: predicted,
            tracking_fraction: self.tracking_fraction,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }
}

impl SpeedDemand {
    /// Track the ORIGINAL allocator's achievable torque, including a current-
    /// or magnetic-support-limited result. This must remain inside the sampled
    /// command envelope; numerical endpoint slack is 32 machine epsilons times
    /// that endpoint's own magnitude. Zero capabilities remain exact. Sampled
    /// power is checked against each corresponding power bound independently.
    /// The value is reported unchanged, never clipped or made into heat.
    ///
    /// Finish before the corresponding native plant step. A cancellation or
    /// refusal discards this proposal; no accepted PI state has been mutated.
    pub fn finish_controlled(
        &self,
        achievable_torque_nm: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SpeedControlProposal, ControlError> {
        poll(&mut cancelled)?;
        if !achievable_torque_nm.is_finite() {
            return Err(ControlError::InvalidSample {
                what: "allocator achievable torque must be finite",
            });
        }
        let allowance = self
            .torque_bounds_nm
            .map(|bound| (32.0 * f64::EPSILON) * bound.abs());
        if achievable_torque_nm < self.torque_bounds_nm[0]
            && self.torque_bounds_nm[0] - achievable_torque_nm > allowance[0]
            || achievable_torque_nm > self.torque_bounds_nm[1]
                && achievable_torque_nm - self.torque_bounds_nm[1] > allowance[1]
        {
            return Err(ControlError::InvalidSample {
                what: "allocator torque leaves the sampled command envelope",
            });
        }
        let integral = self.predicted_integral_torque_nm
            + self.tracking_fraction * (achievable_torque_nm - self.unlimited_torque_nm);
        let power = achievable_torque_nm * self.measured_relative_speed_rad_s;
        finite(
            &[integral, power],
            "speed anti-windup torque and achievable sampled power",
        )?;
        let power_allowance = self
            .shaft_power_bounds_w
            .map(|bound| (32.0 * f64::EPSILON) * bound.abs());
        if power < self.shaft_power_bounds_w[0]
            && self.shaft_power_bounds_w[0] - power > power_allowance[0]
            || power > self.shaft_power_bounds_w[1]
                && power - self.shaft_power_bounds_w[1] > power_allowance[1]
        {
            return Err(ControlError::InvalidSample {
                what: "allocator torque leaves its directional sampled power capability",
            });
        }
        let proposal = SpeedControlProposal {
            next_state: SpeedControllerState {
                integral_torque_nm: integral,
            },
            achievable_torque_nm,
            achievable_shaft_power_w: power,
            allocation_limited: achievable_torque_nm != self.requested_torque_nm,
        };
        poll(&mut cancelled)?;
        Ok(proposal)
    }
}
