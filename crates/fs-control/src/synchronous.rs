//! Sampled field-oriented current control for a balanced synchronous machine.
//!
//! All alpha/beta and d/q values use POWER-INVARIANT Clarke/Park normalization:
//! alpha = sqrt(2/3)*(a-(b+c)/2), beta = (b-c)/sqrt(2), with positive a,b,c
//! phase order. Positive q is a counterclockwise quarter turn from the rotor d
//! axis. Electrical angle is `pole_pairs * relative_mechanical_angle_rad` and
//! electrical speed is `pole_pairs * (omega_rotor - omega_stator)`. There is no
//! 3/2 multiplier and no fixed-stator assumption.
//!
//! The controller uses measured ACTUAL flux in the motional-voltage terms
//! `[-omega_e*lambda_q, +omega_e*lambda_d]`. Consequently its current-regulation
//! path supports both linear and saturating constitutive plants without
//! reconstructing flux from nominal inductance or the commanded current. A
//! machine-owned torque allocator can produce the requested current vector;
//! this module does not replace its nonlinear torque law.
//!
//! The exact balanced two-level inverter voltage hexagon is enforced by phase
//! span, followed by centered common-mode duty injection. This is an averaged
//! held-duty command: no PWM event schedule, dead time or device voltage drop
//! is inferred. The native plant owns actual finite-bus evolution and losses.
//!
//! References for the physical conventions and control structure:
//! * MIT 6.685, Notes 7, section 3.2 (amplitude-invariant in the reference;
//!   converted here to power invariance):
//!   <https://ocw.mit.edu/courses/6-685-electric-machines-fall-2013/c5edfb0f28f7c1b1b4d1588f179074c6_MIT6_685F13_chapter7.pdf>.
//! * Microchip MCAF, field-oriented current control, voltage saturation and
//!   anti-windup: <https://microchiptech.github.io/mcaf-doc/7.0.2/components/foc.html>.
//! * TI SPRA524, physical space-vector geometry (the reference's transform is
//!   amplitude normalized): <https://www.ti.com/lit/an/spra524/spra524.pdf>.
//!
//! The precise predictor/back-calculation recurrence below is this runtime's
//! declared digital law; those references do not certify its gains or a
//! particular closed-loop plant.

use crate::ControlError;
use fs_math::det;

/// Authored fixed-period current-loop data. There are no inferred motor gains,
/// bandwidths or voltage limits; voltage headroom comes from each actual sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CurrentControllerParameters {
    /// Duration of each held command [s], finite and strictly positive.
    pub sample_period_s: f64,
    /// Positive pole-pair count; the input angle is mechanical, not electrical.
    pub pole_pairs: u32,
    /// d/q proportional gains [V/A = ohm], finite and nonnegative.
    pub proportional_gain_ohm: [f64; 2],
    /// d/q integral gains [V/(A s)], finite and nonnegative.
    pub integral_gain_ohm_per_s: [f64; 2],
    /// d/q back-calculation tracking times [s], finite and nonnegative. Zero
    /// requests immediate tracking; a positive value uses Ts/(Taw+Ts).
    pub anti_windup_time_s: [f64; 2],
    /// Per-phase winding resistance [ohm], finite and nonnegative. With power
    /// invariance it is also the d/q resistance. Semiconductor losses are not
    /// folded into this feedforward term.
    pub winding_resistance_ohm: f64,
    /// Maximum magnitude sqrt(id^2+iq^2) of the POWER-INVARIANT current request
    /// [A], finite and strictly positive. It is not a phase-peak or RMS rating.
    /// This command limit does not certify measured overcurrent prevention.
    pub current_limit_a: f64,
}

/// Caller-owned accepted integrator voltage, expressed in the rotating d/q
/// axes. Retain its numerical components as the rotor advances: it is a d/q
/// controller state, not a stationary alpha/beta vector to rotate each sample.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CurrentControllerState {
    /// d/q integrator contribution [V]. Restored values must be finite.
    pub integral_voltage_dq_v: [f64; 2],
}

/// One simultaneous measurement at the START of a held controller interval.
/// Alpha/beta axes move with the stator and use power-invariant normalization.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SynchronousMeasurement {
    /// Actual power-invariant winding current [A].
    pub current_alpha_beta_a: [f64; 2],
    /// Actual power-invariant total flux linkage [Wb], including PM linkage.
    /// It must come from the plant or an explicitly owned estimator, not a
    /// nominal L*commanded_current substitution.
    pub flux_alpha_beta_wb: [f64; 2],
    /// Rotor angle relative to the stator [mechanical rad].
    pub relative_mechanical_angle_rad: f64,
    /// Absolute rotor and stator angular velocities [rad/s], in that order.
    pub shaft_speed_rad_s: [f64; 2],
    /// Sampled physical DC-link voltage [V], finite and nonnegative. Zero is
    /// admitted and emits zero differential voltage with active anti-windup.
    pub dc_bus_voltage_v: f64,
}

/// A complete deterministic sample proposal. Hold `duty_high_fraction` for one
/// authored sample period, then adopt `next_state` only if that plant interval
/// is accepted. Discard the entire proposal if the interval fails or is
/// cancelled; an identical retry from the previous state gives identical bits.
#[derive(Clone, Copy, Debug, PartialEq)]
#[must_use = "hold the proposed duties and accept next_state only with the plant interval"]
pub struct CurrentControlProposal {
    /// Candidate accepted PI state after this held command.
    pub next_state: CurrentControllerState,
    /// Physical high-side duty fractions [a,b,c], each in [0,1]. The low-side
    /// fractions are their complements in the averaged synchronous bridge.
    pub duty_high_fraction: [f64; 3],
    /// Original power-invariant d/q current command [A].
    pub requested_current_dq_a: [f64; 2],
    /// Admitted current command [A]. A request outside the configured circle
    /// is radially scaled, preserving its direction.
    pub current_reference_dq_a: [f64; 2],
    /// Measured current after the actual-angle Park transform [A].
    pub measured_current_dq_a: [f64; 2],
    /// Actual measured flux after the same Park transform [Wb].
    pub actual_flux_dq_wb: [f64; 2],
    /// Motional-voltage feedforward from actual flux and relative speed [V].
    pub speed_voltage_dq_v: [f64; 2],
    /// PI plus winding/speed feedforward before the physical hexagon limit [V].
    pub unlimited_voltage_dq_v: [f64; 2],
    /// Ideal bridge d/q voltage reconstructed from the actual proposed duties
    /// and sampled DC voltage [V]. This precedes device drops and future bus
    /// evolution, so it is not a prediction of held-interval terminal voltage.
    pub commanded_voltage_dq_v: [f64; 2],
    /// The same ideal bridge command in stationary stator alpha/beta axes [V].
    pub commanded_voltage_alpha_beta_v: [f64; 2],
    /// Whether the requested current vector exceeded its authored magnitude.
    pub current_limited: bool,
    /// Whether the unlimited voltage exceeded the sampled physical hexagon.
    pub voltage_limited: bool,
}

/// Immutable prepared current-loop coefficients. Each sample uses only fixed
/// arrays and bounded arithmetic, with no allocation, solver or hidden state.
#[derive(Clone, Debug)]
pub struct SynchronousCurrentController {
    parameters: CurrentControllerParameters,
    integral_increment_gain_ohm: [f64; 2],
    tracking_fraction: [f64; 2],
}

fn finite(values: &[f64], what: &'static str) -> Result<(), ControlError> {
    if values.iter().all(|x| x.is_finite()) {
        Ok(())
    } else {
        Err(ControlError::NumericalRange { what })
    }
}

fn park(vector: [f64; 2], sin_phi: f64, cos_phi: f64) -> [f64; 2] {
    [
        cos_phi * vector[0] + sin_phi * vector[1],
        -sin_phi * vector[0] + cos_phi * vector[1],
    ]
}

fn inverse_park(vector: [f64; 2], sin_phi: f64, cos_phi: f64) -> [f64; 2] {
    [
        cos_phi * vector[0] - sin_phi * vector[1],
        sin_phi * vector[0] + cos_phi * vector[1],
    ]
}

/// Map a balanced alpha/beta request to the exact physical voltage hexagon.
/// A common offset does not change winding voltage; choose it to center the
/// largest and smallest duties. If the request exceeds available phase span,
/// one common radial factor scales the entire requested differential vector.
fn modulate(
    requested: [f64; 2],
    dc_voltage: f64,
) -> Result<([f64; 3], [f64; 2], bool), ControlError> {
    let alpha_scale = det::sqrt(2.0 / 3.0);
    let beta_scale = det::sqrt(0.5);
    let phase = [
        alpha_scale * requested[0],
        -0.5 * alpha_scale * requested[0] + beta_scale * requested[1],
        -0.5 * alpha_scale * requested[0] - beta_scale * requested[1],
    ];
    finite(&phase, "inverse Clarke voltage")?;
    let low = phase[0].min(phase[1]).min(phase[2]);
    let high = phase[0].max(phase[1]).max(phase[2]);
    let span = high - low;
    finite(&[span], "phase voltage span")?;
    let limited = span > dc_voltage;
    if dc_voltage == 0.0 {
        return Ok(([0.5; 3], [0.0; 2], limited));
    }
    let divisor = if limited { span } else { dc_voltage };
    let padding = if limited {
        0.0
    } else {
        0.5 * (1.0 - span / dc_voltage)
    };
    // Using differences from the minimum avoids independent phase clipping:
    // every difference is in [0,span], and the same divisor owns all legs.
    let duty = phase.map(|v| padding + (v - low) / divisor);
    if duty
        .iter()
        .any(|v| !v.is_finite() || !(0.0..=1.0).contains(v))
    {
        return Err(ControlError::NumericalRange {
            what: "physical inverter duty fractions",
        });
    }
    let commanded = [
        dc_voltage * (alpha_scale * (duty[0] - f64::midpoint(duty[1], duty[2]))),
        dc_voltage * (beta_scale * (duty[1] - duty[2])),
    ];
    finite(&commanded, "commanded bridge voltage")?;
    Ok((duty, commanded, limited))
}

impl SynchronousCurrentController {
    /// Prepare the authored digital gains. This verifies physical domains and
    /// representability, not closed-loop stability or machine current ratings.
    pub fn new(parameters: CurrentControllerParameters) -> Result<Self, ControlError> {
        let ts = parameters.sample_period_s;
        if !ts.is_finite() || ts <= 0.0 {
            return Err(ControlError::InvalidParameter {
                what: "sample_period_s must be finite and strictly positive",
            });
        }
        if parameters.pole_pairs == 0 {
            return Err(ControlError::InvalidParameter {
                what: "pole_pairs must be positive",
            });
        }
        if !parameters.current_limit_a.is_finite() || parameters.current_limit_a <= 0.0 {
            return Err(ControlError::InvalidParameter {
                what: "power-invariant current_limit_a must be finite and positive",
            });
        }
        if parameters
            .proportional_gain_ohm
            .iter()
            .chain(&parameters.integral_gain_ohm_per_s)
            .chain(&parameters.anti_windup_time_s)
            .chain(std::iter::once(&parameters.winding_resistance_ohm))
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(ControlError::InvalidParameter {
                what: "PI gains, tracking times and winding resistance must be finite nonnegative",
            });
        }
        let integral_increment_gain_ohm = parameters.integral_gain_ohm_per_s.map(|v| ts * v);
        let tracking_fraction = parameters.anti_windup_time_s.map(|tau| {
            // Equivalent to Ts/(Taw+Ts), without overflowing the sum or ratio.
            if tau <= ts {
                1.0 / (1.0 + tau / ts)
            } else {
                let ratio = ts / tau;
                ratio / (1.0 + ratio)
            }
        });
        finite(&integral_increment_gain_ohm, "sampled integral gain")?;
        for axis in 0..2 {
            if (parameters.integral_gain_ohm_per_s[axis] > 0.0
                && integral_increment_gain_ohm[axis] == 0.0)
                || tracking_fraction[axis] <= 0.0
            {
                return Err(ControlError::NumericalRange {
                    what: "sampled integral or tracking gain underflow",
                });
            }
        }
        Ok(Self {
            parameters,
            integral_increment_gain_ohm,
            tracking_fraction,
        })
    }

    /// Authored fixed-period data, including the current-command limit.
    #[must_use]
    pub fn parameters(&self) -> CurrentControllerParameters {
        self.parameters
    }

    /// Propose one fixed-period command without changing accepted state.
    ///
    /// With e = i_ref-i, z_pred = z+Ts*Ki*e, u_raw = R*i_ref +
    /// [-omega_e*lambda_q,+omega_e*lambda_d] + Kp*e + z_pred. Physical
    /// modulation returns u_cmd. Back calculation uses
    /// z_next = z_pred + gamma*(u_cmd-u_raw), gamma=Ts/(Taw+Ts).
    /// With error, feedforward and saturated commanded voltage held fixed, the
    /// tracking homogeneous pole is 1-gamma in [0,1): a large sample/tracking-
    /// time ratio cannot destabilize this scalar tracking recurrence. The
    /// changing hexagon and coupled plant do not inherit a stability certificate.
    ///
    /// Cancellation is polled before arithmetic and immediately before return.
    /// There is no unbounded work, heap allocation or mutation between polls.
    /// Hold the returned duties for exactly `sample_period_s`; adapting plant
    /// substeps within that interval does not authorize another PI update.
    pub fn propose_controlled(
        &self,
        state: &CurrentControllerState,
        measurement: &SynchronousMeasurement,
        requested_current_dq_a: [f64; 2],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CurrentControlProposal, ControlError> {
        if cancelled() {
            return Err(ControlError::Cancelled);
        }
        if state
            .integral_voltage_dq_v
            .iter()
            .chain(&measurement.current_alpha_beta_a)
            .chain(&measurement.flux_alpha_beta_wb)
            .chain(&measurement.shaft_speed_rad_s)
            .chain(&requested_current_dq_a)
            .any(|v| !v.is_finite())
            || !measurement.relative_mechanical_angle_rad.is_finite()
            || !measurement.dc_bus_voltage_v.is_finite()
            || measurement.dc_bus_voltage_v < 0.0
        {
            return Err(ControlError::InvalidSample {
                what: "state, request and measurements must be finite; DC voltage must be nonnegative",
            });
        }
        let magnitude = det::hypot(requested_current_dq_a[0], requested_current_dq_a[1]);
        finite(&[magnitude], "requested current magnitude")?;
        let current_limited = magnitude > self.parameters.current_limit_a;
        let current_reference_dq_a = if current_limited {
            let scale = self.parameters.current_limit_a / magnitude;
            requested_current_dq_a.map(|v| v * scale)
        } else {
            requested_current_dq_a
        };
        let pole_pairs = f64::from(self.parameters.pole_pairs);
        let phi = pole_pairs * measurement.relative_mechanical_angle_rad;
        let omega_e =
            pole_pairs * (measurement.shaft_speed_rad_s[0] - measurement.shaft_speed_rad_s[1]);
        finite(&[phi, omega_e], "electrical angle and relative speed")?;
        let (sin_phi, cos_phi) = det::sin_cos(phi);
        let measured_current_dq_a = park(measurement.current_alpha_beta_a, sin_phi, cos_phi);
        let actual_flux_dq_wb = park(measurement.flux_alpha_beta_wb, sin_phi, cos_phi);
        let speed_voltage_dq_v = [
            -omega_e * actual_flux_dq_wb[1],
            omega_e * actual_flux_dq_wb[0],
        ];
        let mut predicted_integral = [0.0; 2];
        let mut unlimited_voltage_dq_v = [0.0; 2];
        for axis in 0..2 {
            let error = current_reference_dq_a[axis] - measured_current_dq_a[axis];
            predicted_integral[axis] =
                state.integral_voltage_dq_v[axis] + self.integral_increment_gain_ohm[axis] * error;
            unlimited_voltage_dq_v[axis] = self.parameters.winding_resistance_ohm
                * current_reference_dq_a[axis]
                + speed_voltage_dq_v[axis]
                + self.parameters.proportional_gain_ohm[axis] * error
                + predicted_integral[axis];
        }
        finite(&measured_current_dq_a, "measured current Park transform")?;
        finite(&actual_flux_dq_wb, "actual flux Park transform")?;
        finite(&predicted_integral, "predicted integrator voltage")?;
        finite(&unlimited_voltage_dq_v, "unlimited controller voltage")?;
        let requested_alpha_beta = inverse_park(unlimited_voltage_dq_v, sin_phi, cos_phi);
        let (duty_high_fraction, commanded_voltage_alpha_beta_v, voltage_limited) =
            modulate(requested_alpha_beta, measurement.dc_bus_voltage_v)?;
        let commanded_voltage_dq_v = park(commanded_voltage_alpha_beta_v, sin_phi, cos_phi);
        let mut integral_voltage_dq_v = [0.0; 2];
        for axis in 0..2 {
            integral_voltage_dq_v[axis] = predicted_integral[axis]
                + self.tracking_fraction[axis]
                    * (commanded_voltage_dq_v[axis] - unlimited_voltage_dq_v[axis]);
        }
        finite(&integral_voltage_dq_v, "anti-windup integrator voltage")?;
        let proposal = CurrentControlProposal {
            next_state: CurrentControllerState {
                integral_voltage_dq_v,
            },
            duty_high_fraction,
            requested_current_dq_a,
            current_reference_dq_a,
            measured_current_dq_a,
            actual_flux_dq_wb,
            speed_voltage_dq_v,
            unlimited_voltage_dq_v,
            commanded_voltage_dq_v,
            commanded_voltage_alpha_beta_v,
            current_limited,
            voltage_limited,
        };
        if cancelled() {
            return Err(ControlError::Cancelled);
        }
        Ok(proposal)
    }
}
