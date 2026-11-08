//! Exact linear response to the piecewise-linear acceleration actually supplied.
//!
//! On each interval acceleration satisfies the homogeneous oscillator equation.
//! Its analytic zeros locate relative-velocity extrema and partition velocity
//! into monotone pieces for displacement extrema. Absolute acceleration is
//! `-2*zeta*omega*v - omega²*u`; zeros of its homogeneous second derivative
//! partition its first derivative for the same bounded extremum search.
//! This resolves intersample peaks, not unmeasured ground-motion frequencies.

use fs_math::det;
use fs_modal::time::{LinearForceTransition, ModalAcousticState, ModalAcousticTimeError};

/// Explicit work limits and numerical extremum-localization tolerance.
#[derive(Clone, Copy, Debug)]
pub struct LinearSdofControls {
    /// Maximum supplied acceleration samples, checked before scanning the record.
    pub maximum_samples: usize,
    /// Maximum interior analytic zeros in either acceleration or its second derivative
    /// in any one input interval. Each partition has this cap independently.
    pub maximum_acceleration_zeros_per_interval: usize,
    /// Total full-interval and search transition evaluations.
    pub maximum_transition_evaluations: usize,
    /// Bisection ceiling for each displacement or absolute-acceleration extremum.
    pub maximum_bisections_per_extremum: usize,
    /// Maximum width [s] of an extremum bracket before evaluating its midpoint.
    pub extremum_time_tolerance_s: f64,
}

/// Physical relative motion at a record boundary.
///
/// Displacement and velocity are in SI units, with no modal mass normalization.
/// The zero/default state is rest. State admission occurs at the response call.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct LinearSdofState {
    /// Relative displacement [m].
    pub displacement_m: f64,
    /// Relative velocity [m/s].
    pub velocity_m_s: f64,
}

/// Response over exactly the supplied record duration, with segment-local times.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearSdofResponse {
    /// Maximum absolute relative displacement [m], including intersample extrema.
    pub peak_displacement_m: f64,
    /// Signed displacement [m] at the selected extremum or record endpoint.
    pub signed_peak_displacement_m: f64,
    /// Time [s] from this segment's start of the selected peak; ties select the first.
    pub peak_time_s: f64,
    /// Largest final velocity-zero bracket width [s].
    pub maximum_extremum_bracket_s: f64,
    /// Maximum localization allowance [m] from bracket width and endpoint speeds.
    /// This excludes floating-point transition/root errors and input/model errors;
    /// it is not a certified enclosure of the complete physical response.
    pub peak_localization_allowance_m: f64,
    /// Maximum absolute relative velocity [m/s], including analytic acceleration
    /// zeros. This is not pseudo-spectral velocity `omega * peak_displacement`.
    pub peak_velocity_m_s: f64,
    /// Signed relative velocity [m/s] at its selected peak.
    pub signed_peak_velocity_m_s: f64,
    /// Time [s] from this segment's start of the velocity peak; ties select the first.
    /// Its analytic localization still has floating-point transition/root error.
    pub peak_velocity_time_s: f64,
    /// Maximum absolute inertial acceleration [m/s²], including intersample
    /// extrema of `-2*zeta*omega*v - omega²*u`.
    pub peak_absolute_acceleration_m_s2: f64,
    /// Signed inertial acceleration [m/s²] at its selected peak.
    pub signed_peak_absolute_acceleration_m_s2: f64,
    /// Time [s] from this segment's start of the acceleration peak; ties select the first.
    pub peak_absolute_acceleration_time_s: f64,
    /// Largest final bracket width [s] for an absolute-acceleration extremum.
    pub maximum_absolute_acceleration_extremum_bracket_s: f64,
    /// Maximum localization allowance [m/s²] from bracket widths and endpoint
    /// acceleration rates. Excludes floating-point and input/model errors.
    pub peak_absolute_acceleration_localization_allowance_m_s2: f64,
    /// Endpoint relative displacement [m].
    pub final_displacement_m: f64,
    /// Endpoint relative velocity [m/s].
    pub final_velocity_m_s: f64,
    /// Exactly `(sample_count - 1) * dt` for a nonempty record [s].
    pub duration_s: f64,
    /// Number of original input intervals advanced.
    pub intervals: usize,
    /// Full steps and partial-step evaluations actually performed.
    pub transition_evaluations: usize,
}

/// Invalid inputs, exhausted budgets and cancellation return no partial response.
#[derive(Debug)]
pub enum LinearSdofError {
    /// Invalid or unrepresentable physical/numerical input.
    Invalid(&'static str),
    /// An explicit caller work limit was exhausted.
    Budget(&'static str),
    /// Requested cancellation was observed before publication.
    Cancelled,
    /// Original exact oscillator owner refused the transition.
    Transition(ModalAcousticTimeError),
}

impl core::fmt::Display for LinearSdofError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-UQ-LINEAR-SDOF: {self:?}")
    }
}
impl std::error::Error for LinearSdofError {}
impl From<ModalAcousticTimeError> for LinearSdofError {
    fn from(error: ModalAcousticTimeError) -> Self {
        Self::Transition(error)
    }
}

struct Work<C> {
    controls: LinearSdofControls,
    evaluations: usize,
    cancelled: C,
}
impl<C: FnMut() -> bool> Work<C> {
    fn poll(&mut self) -> Result<(), LinearSdofError> {
        if (self.cancelled)() {
            Err(LinearSdofError::Cancelled)
        } else {
            Ok(())
        }
    }
    fn evaluation(&mut self) -> Result<(), LinearSdofError> {
        self.poll()?;
        if self.evaluations == self.controls.maximum_transition_evaluations {
            return Err(LinearSdofError::Budget("transition evaluations"));
        }
        self.evaluations += 1;
        Ok(())
    }
}

struct Interval<'a, C> {
    start: ModalAcousticState,
    end: ModalAcousticState,
    force_start: f64,
    force_end: f64,
    dt: f64,
    omega: f64,
    zeta: f64,
    work: &'a mut Work<C>,
}

#[derive(Clone, Copy)]
enum Extremum {
    Displacement,
    AbsoluteAcceleration,
}

impl<C: FnMut() -> bool> Interval<'_, C> {
    fn at(&mut self, t: f64) -> Result<ModalAcousticState, LinearSdofError> {
        if t == 0.0 {
            return Ok(self.start);
        }
        if t == self.dt {
            return Ok(self.end);
        }
        self.work.evaluation()?;
        let fraction = t / self.dt;
        let force = (1.0 - fraction) * self.force_start + fraction * self.force_end;
        Ok(
            LinearForceTransition::try_new(self.omega, self.zeta, t)?.advance(
                self.start,
                self.force_start,
                force,
            )?,
        )
    }

    fn derivative(
        &self,
        quantity: Extremum,
        state: ModalAcousticState,
        time: f64,
    ) -> Result<f64, LinearSdofError> {
        match quantity {
            Extremum::Displacement => Ok(state.velocity_m_sqrt_kg_per_s),
            Extremum::AbsoluteAcceleration => {
                let damping = 2.0 * self.omega * self.zeta;
                let stiffness = self.omega * self.omega;
                let fraction = time / self.dt;
                let force = finite((1.0 - fraction) * self.force_start + fraction * self.force_end)?;
                let acceleration = finite(force + absolute_acceleration(state, damping, stiffness)?)?;
                finite(-damping * acceleration - stiffness * state.velocity_m_sqrt_kg_per_s)
            }
        }
    }

    fn extremum(
        &mut self,
        quantity: Extremum,
        mut low: f64,
        mut high: f64,
        mut low_derivative: f64,
        mut high_derivative: f64,
    ) -> Result<(f64, ModalAcousticState, f64, f64), LinearSdofError> {
        for _ in 0..self.work.controls.maximum_bisections_per_extremum {
            self.work.poll()?;
            let width = high - low;
            let middle = low + 0.5 * width;
            if width <= self.work.controls.extremum_time_tolerance_s {
                let allowance = 0.5 * width * low_derivative.abs().max(high_derivative.abs());
                return Ok((middle, self.at(middle)?, width, finite(allowance)?));
            }
            if middle == low || middle == high {
                return Err(LinearSdofError::Invalid(
                    "extremum tolerance is below representable time resolution",
                ));
            }
            let state = self.at(middle)?;
            let derivative = self.derivative(quantity, state, middle)?;
            if derivative == 0.0 {
                return Ok((middle, state, 0.0, 0.0));
            }
            if derivative.is_sign_positive() == low_derivative.is_sign_positive() {
                low = middle;
                low_derivative = derivative;
            } else {
                high = middle;
                high_derivative = derivative;
            }
        }
        Err(LinearSdofError::Budget("extremum bisections"))
    }
}

fn finite(value: f64) -> Result<f64, LinearSdofError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(LinearSdofError::Invalid(
            "unrepresentable oscillator/extremum arithmetic",
        ))
    }
}

fn absolute_acceleration(
    state: ModalAcousticState,
    damping: f64,
    stiffness: f64,
) -> Result<f64, LinearSdofError> {
    finite(-damping * state.velocity_m_sqrt_kg_per_s - stiffness * state.displacement_m_sqrt_kg)
}

// atanh(x)/x with its regular x=0 limit, avoiding subtractive logarithms near zero.
fn atanh_over_x(x: f64) -> f64 {
    if x < 0.125 {
        let square = x * x;
        let mut value = 1.0 / 19.0;
        for denominator in (1_u8..=17).step_by(2).rev() {
            value = 1.0 / f64::from(denominator) + square * value;
        }
        value
    } else {
        0.5 * (det::ln(1.0 + x) - det::ln(1.0 - x)) / x
    }
}

// Zeros of the derivative that partitions the selected extremum's first
// derivative into monotone pieces. For displacement this is u''; for inertial
// acceleration y = u'' - f, this is y'' = u'''' because the force is affine.
fn partition_zeros(
    quantity: Extremum,
    state: ModalAcousticState,
    f0: f64,
    f1: f64,
    dt: f64,
    omega: f64,
    zeta: f64,
) -> Result<(f64, f64), LinearSdofError> {
    let alpha = omega * zeta;
    let stiffness = omega * omega;
    let velocity = state.velocity_m_sqrt_kg_per_s;
    let a = finite(f0 + absolute_acceleration(state, 2.0 * alpha, stiffness)?)?;
    let slope = finite((f1 - f0) / dt)?;
    let (a, d) = match quantity {
        Extremum::Displacement => {
            // a'(0) + alpha*a(0), without cancelling 2*alpha*a.
            (a, finite(slope - alpha * a - stiffness * velocity)?)
        }
        Extremum::AbsoluteAcceleration => {
            let acceleration_rate = finite(slope - 2.0 * alpha * a - stiffness * velocity)?;
            let curvature = finite(-2.0 * alpha * acceleration_rate - stiffness * a)?;
            // u''''' + alpha*u'''' = -alpha*u'''' - omega²*u'''.
            (curvature, finite(-alpha * curvature - stiffness * acceleration_rate)?)
        }
    };
    homogeneous_zeros(a, d, omega, zeta)
}

// For x'' + 2*alpha*x' + omega²*x = 0, a = x(0), d = x'(0)+alpha*x(0).
// Returns the first strictly positive zero and its spacing. Critical and
// overdamped responses have at most one zero (spacing = infinity).
fn homogeneous_zeros(
    a: f64,
    d: f64,
    omega: f64,
    zeta: f64,
) -> Result<(f64, f64), LinearSdofError> {
    let alpha = omega * zeta;
    if a == 0.0 && d == 0.0 {
        return Ok((f64::INFINITY, f64::INFINITY));
    }
    if zeta < 1.0 {
        let frequency = finite(omega * ((1.0 - zeta) * (1.0 + zeta)).sqrt())?;
        let b = finite(d / frequency)?;
        let mut phase = det::atan2(-a, b).rem_euclid(core::f64::consts::PI);
        if phase == 0.0 {
            phase = core::f64::consts::PI;
        }
        Ok((phase / frequency, core::f64::consts::PI / frequency))
    } else {
        if d == 0.0 || a == 0.0 || a.is_sign_positive() == d.is_sign_positive() {
            return Ok((f64::INFINITY, f64::INFINITY));
        }
        let critical_time = -a / d;
        if zeta == 1.0 {
            return Ok((critical_time, f64::INFINITY));
        }
        let frequency = if zeta < 2.0 {
            omega * (zeta - 1.0).sqrt() * (zeta + 1.0).sqrt()
        } else {
            alpha * ((1.0 - 1.0 / zeta) * (1.0 + 1.0 / zeta)).sqrt()
        };
        let ratio = critical_time * frequency;
        if ratio >= 1.0 {
            return Ok((f64::INFINITY, f64::INFINITY));
        }
        Ok((finite(critical_time * atanh_over_x(ratio))?, f64::INFINITY))
    }
}

impl LinearSdofResponse {
    /// Physical endpoint state for continuing with another supplied segment.
    ///
    /// Continuation preserves motion; peak statistics and times belong to each
    /// segment separately. No forcing interval or decay tail is appended.
    #[must_use]
    pub const fn final_state(&self) -> LinearSdofState {
        LinearSdofState {
            displacement_m: self.final_displacement_m,
            velocity_m_s: self.final_velocity_m_s,
        }
    }

    fn observe(
        &mut self,
        state: ModalAcousticState,
        time: f64,
        damping: f64,
        stiffness: f64,
    ) -> Result<(), LinearSdofError> {
        let value = state.displacement_m_sqrt_kg;
        if value.abs() > self.peak_displacement_m
            || (value.abs() == self.peak_displacement_m && time < self.peak_time_s)
        {
            self.peak_displacement_m = value.abs();
            self.signed_peak_displacement_m = value;
            self.peak_time_s = time;
        }
        let velocity = state.velocity_m_sqrt_kg_per_s;
        if velocity.abs() > self.peak_velocity_m_s
            || (velocity.abs() == self.peak_velocity_m_s && time < self.peak_velocity_time_s)
        {
            self.peak_velocity_m_s = velocity.abs();
            self.signed_peak_velocity_m_s = velocity;
            self.peak_velocity_time_s = time;
        }
        let acceleration = absolute_acceleration(state, damping, stiffness)?;
        if acceleration.abs() > self.peak_absolute_acceleration_m_s2
            || (acceleration.abs() == self.peak_absolute_acceleration_m_s2
                && time < self.peak_absolute_acceleration_time_s)
        {
            self.peak_absolute_acceleration_m_s2 = acceleration.abs();
            self.signed_peak_absolute_acceleration_m_s2 = acceleration;
            self.peak_absolute_acceleration_time_s = time;
        }
        Ok(())
    }
}

/// Solve `u'' + 2*zeta*omega*u' + omega²*u = -ground_acceleration` from rest.
///
/// Uses the shared exact first-order-hold oscillator transition, cached for full
/// input intervals. The acceleration samples are endpoints of straight segments:
/// N samples cover N-1 intervals, with no implicit decay tail. Empty or single-
/// sample records return zero motion. Units are seconds, rad/s and m/s².
///
/// Displacement and absolute-acceleration extrema are sought on analytic
/// monotone-derivative pieces. Relative-velocity extrema use analytic acceleration
/// zeros directly. The tolerance controls root localization, not input bandwidth
/// or model error. Physical velocity and acceleration peaks are distinct from
/// the pseudo-spectral quantities obtained by scaling peak displacement.
/// Work is bounded and storage is constant, independent of the record length.
/// Deterministic repeatability is on the same ISA/toolchain; no interval-certified
/// extrema or cross-ISA claim is made. The legacy Newmark APIs remain separate.
///
/// # Errors
/// Refuses invalid inputs, unrepresentable extrema arithmetic, exhausted work
/// limits, cancellation, or an exact-transition refusal. Returns no partial result.
pub fn linear_sdof_response<C: FnMut() -> bool>(
    ground_acceleration_m_s2: &[f64],
    dt_s: f64,
    omega_rad_s: f64,
    damping_ratio: f64,
    controls: LinearSdofControls,
    cancelled: C,
) -> Result<LinearSdofResponse, LinearSdofError> {
    linear_sdof_response_from_state(
        ground_acceleration_m_s2,
        dt_s,
        omega_rad_s,
        damping_ratio,
        LinearSdofState::default(),
        controls,
        cancelled,
    )
}

/// Continue exact piecewise-linear response from supplied physical motion.
///
/// Uses the same transitions, continuous extrema and work limits as
/// [`linear_sdof_response`]. All peak searches include the initial state at
/// time zero. Peak times are relative to this segment, and empty/single-sample
/// segments preserve and report the supplied motion without advancing time.
///
/// To split a record, share its boundary sample between adjacent segments and
/// pass [`LinearSdofResponse::final_state`] into the next call. This retains
/// displacement and velocity without silently resetting the oscillator. Each
/// call has its own work budget and peak statistics. Explicit zero-acceleration
/// samples model free decay; no gap or tail is inferred between records.
///
/// # Errors
/// Refuses nonfinite initial motion in addition to the invalid inputs,
/// unrepresentable arithmetic, work limits and cancellation of the rest API.
#[allow(clippy::cast_precision_loss)]
pub fn linear_sdof_response_from_state<C: FnMut() -> bool>(
    ground_acceleration_m_s2: &[f64],
    dt_s: f64,
    omega_rad_s: f64,
    damping_ratio: f64,
    initial_state: LinearSdofState,
    controls: LinearSdofControls,
    cancelled: C,
) -> Result<LinearSdofResponse, LinearSdofError> {
    let mut work = Work {
        controls,
        evaluations: 0,
        cancelled,
    };
    work.poll()?;
    if ground_acceleration_m_s2.len() > controls.maximum_samples {
        return Err(LinearSdofError::Budget("acceleration samples"));
    }
    if !dt_s.is_finite()
        || dt_s <= 0.0
        || !omega_rad_s.is_finite()
        || omega_rad_s <= 0.0
        || !damping_ratio.is_finite()
        || damping_ratio < 0.0
        || !controls.extremum_time_tolerance_s.is_finite()
        || controls.extremum_time_tolerance_s <= 0.0
        || controls.maximum_transition_evaluations == 0
        || controls.maximum_bisections_per_extremum == 0
        || omega_rad_s * omega_rad_s == 0.0
    {
        return Err(LinearSdofError::Invalid(
            "positive finite clock/frequency/tolerance and valid damping/work limits required",
        ));
    }
    if !initial_state.displacement_m.is_finite() || !initial_state.velocity_m_s.is_finite() {
        return Err(LinearSdofError::Invalid(
            "initial displacement and velocity must be finite",
        ));
    }
    let stiffness = finite(omega_rad_s * omega_rad_s)?;
    let damping = finite(2.0 * omega_rad_s * damping_ratio)?;
    let intervals = ground_acceleration_m_s2.len().saturating_sub(1);
    let duration_s = finite(intervals as f64 * dt_s)?;
    for &sample in ground_acceleration_m_s2 {
        work.poll()?;
        if !sample.is_finite() {
            return Err(LinearSdofError::Invalid(
                "ground acceleration must be finite",
            ));
        }
    }
    let full = LinearForceTransition::try_new(omega_rad_s, damping_ratio, dt_s)?;
    let mut state = ModalAcousticState {
        displacement_m_sqrt_kg: initial_state.displacement_m,
        velocity_m_sqrt_kg_per_s: initial_state.velocity_m_s,
    };
    let mut report = LinearSdofResponse {
        peak_displacement_m: 0.0,
        signed_peak_displacement_m: 0.0,
        peak_time_s: 0.0,
        maximum_extremum_bracket_s: 0.0,
        peak_localization_allowance_m: 0.0,
        peak_velocity_m_s: 0.0,
        signed_peak_velocity_m_s: 0.0,
        peak_velocity_time_s: 0.0,
        peak_absolute_acceleration_m_s2: 0.0,
        signed_peak_absolute_acceleration_m_s2: 0.0,
        peak_absolute_acceleration_time_s: 0.0,
        maximum_absolute_acceleration_extremum_bracket_s: 0.0,
        peak_absolute_acceleration_localization_allowance_m_s2: 0.0,
        final_displacement_m: 0.0,
        final_velocity_m_s: 0.0,
        duration_s,
        intervals,
        transition_evaluations: 0,
    };
    report.observe(state, 0.0, damping, stiffness)?;
    for (index, samples) in ground_acceleration_m_s2.windows(2).enumerate() {
        work.evaluation()?;
        let (f0, f1) = (-samples[0], -samples[1]);
        let end = full.advance(state, f0, f1)?;
        let mut interval = Interval {
            start: state,
            end,
            force_start: f0,
            force_end: f1,
            dt: dt_s,
            omega: omega_rad_s,
            zeta: damping_ratio,
            work: &mut work,
        };
        let base_time = index as f64 * dt_s;
        for quantity in [Extremum::Displacement, Extremum::AbsoluteAcceleration] {
            // With no damping y=-omega²*u, so the displacement search already
            // visits every acceleration extremum; do not repeat its transitions.
            if matches!(quantity, Extremum::AbsoluteAcceleration) && damping == 0.0 {
                continue;
            }
            let (first_zero, spacing) =
                partition_zeros(quantity, state, f0, f1, dt_s, omega_rad_s, damping_ratio)?;
            let mut low = 0.0;
            let mut low_derivative = interval.derivative(quantity, state, low)?;
            let mut zero_index = 0usize;
            loop {
                interval.work.poll()?;
                let next_zero = if zero_index == 0 {
                    first_zero
                } else {
                    first_zero + zero_index as f64 * spacing
                };
                let high = next_zero.min(dt_s);
                if high <= low {
                    return Err(LinearSdofError::Invalid(
                        "derivative-zero spacing is below representable time resolution",
                    ));
                }
                if high < dt_s && zero_index >= controls.maximum_acceleration_zeros_per_interval {
                    return Err(LinearSdofError::Budget(match quantity {
                        Extremum::Displacement => "acceleration zeros per input interval",
                        Extremum::AbsoluteAcceleration => {
                            "absolute-acceleration curvature zeros per input interval"
                        }
                    }));
                }
                let high_state = interval.at(high)?;
                report.observe(high_state, base_time + high, damping, stiffness)?;
                let high_derivative = interval.derivative(quantity, high_state, high)?;
                if low_derivative != 0.0 && high_derivative != 0.0
                    && low_derivative.is_sign_positive() != high_derivative.is_sign_positive()
                {
                    let (time, peak, width, allowance) =
                        interval.extremum(quantity, low, high, low_derivative, high_derivative)?;
                    report.observe(peak, base_time + time, damping, stiffness)?;
                    match quantity {
                        Extremum::Displacement => {
                            report.maximum_extremum_bracket_s =
                                report.maximum_extremum_bracket_s.max(width);
                            report.peak_localization_allowance_m =
                                report.peak_localization_allowance_m.max(allowance);
                        }
                        Extremum::AbsoluteAcceleration => {
                            report.maximum_absolute_acceleration_extremum_bracket_s =
                                report.maximum_absolute_acceleration_extremum_bracket_s.max(width);
                            report.peak_absolute_acceleration_localization_allowance_m_s2 =
                                report.peak_absolute_acceleration_localization_allowance_m_s2.max(allowance);
                        }
                    }
                }
                if high == dt_s {
                    break;
                }
                low = high;
                low_derivative = high_derivative;
                zero_index += 1;
            }
        }
        state = end;
    }
    work.poll()?;
    report.final_displacement_m = state.displacement_m_sqrt_kg;
    report.final_velocity_m_s = state.velocity_m_sqrt_kg_per_s;
    if damping == 0.0 {
        report.maximum_absolute_acceleration_extremum_bracket_s = report.maximum_extremum_bracket_s;
        report.peak_absolute_acceleration_localization_allowance_m_s2 =
            finite(stiffness * report.peak_localization_allowance_m)?;
    }
    report.transition_evaluations = work.evaluations;
    Ok(report)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn controls() -> LinearSdofControls {
        LinearSdofControls {
            maximum_samples: 10_000,
            maximum_acceleration_zeros_per_interval: 256,
            maximum_transition_evaluations: 100_000,
            maximum_bisections_per_extremum: 128,
            extremum_time_tolerance_s: 1.0e-11,
        }
    }
    fn solve(samples: &[f64], dt: f64, omega: f64, zeta: f64) -> LinearSdofResponse {
        linear_sdof_response(samples, dt, omega, zeta, controls(), || false).unwrap()
    }
    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
        );
    }

    #[test]
    fn initial_velocity_free_decay_resolves_intersample_peaks() {
        let omega = 4.0_f64;
        let speed = 3.0_f64;
        let initial = LinearSdofState {
            displacement_m: 0.0,
            velocity_m_s: speed,
        };
        for zeta in [0.0_f64, 0.2] {
            let alpha = zeta * omega;
            let wd = omega * (1.0 - zeta * zeta).sqrt();
            let duration = core::f64::consts::PI / wd;
            let response = linear_sdof_response_from_state(
                &[0.0, 0.0], duration, omega, zeta, initial, controls(), || false,
            )
            .unwrap();
            // q(t) = V exp(-alpha*t) sin(wd*t)/wd: the sampled
            // displacements vanish while the displacement peak lies inside.
            let displacement_time = zeta.acos() / wd;
            let displacement = speed / omega * (-alpha * displacement_time).exp();
            close(response.final_displacement_m, 0.0, 2.0e-14);
            close(response.final_velocity_m_s, -speed * (-alpha * duration).exp(), 2.0e-13);
            close(response.peak_displacement_m, displacement, 2.0e-13);
            close(response.signed_peak_displacement_m, displacement, 2.0e-13);
            close(response.peak_time_s, displacement_time, 2.0e-11);
            close(response.peak_velocity_m_s, speed, 2.0e-13);
            // Differentiate the closed-form free response independently:
            // a/V = exp(-alpha*t)[-2*alpha*cos + (2*alpha^2-omega^2)*sin/wd].
            let angle = ((omega * omega - 4.0 * alpha * alpha) * wd)
                .atan2(alpha * (3.0 * omega * omega - 4.0 * alpha * alpha));
            let acceleration_time = angle / wd;
            let acceleration = speed * (-alpha * acceleration_time).exp()
                * (-2.0 * alpha * angle.cos()
                    + (2.0 * alpha * alpha - omega * omega) * angle.sin() / wd);
            close(response.peak_absolute_acceleration_m_s2, acceleration.abs(), 2.0e-12);
            close(response.signed_peak_absolute_acceleration_m_s2, acceleration, 2.0e-12);
            close(response.peak_absolute_acceleration_time_s, acceleration_time, 2.0e-11);
            assert_eq!(response.duration_s, duration);
            assert_eq!(response.intervals, 1);
            assert_eq!(
                response,
                linear_sdof_response_from_state(
                    &[0.0, 0.0], duration, omega, zeta, initial, controls(), || false,
                )
                .unwrap(),
            );
        }
    }

    #[test]
    fn critical_and_overdamped_initial_eigenstates_keep_initial_peaks() {
        let omega = 4.0_f64;
        let duration = 1.25;
        for zeta in [1.0_f64, 2.0] {
            // Select one exact decaying eigenmode, including the repeated root.
            let rate = -omega / (zeta + (zeta * zeta - 1.0).sqrt());
            let initial = LinearSdofState {
                displacement_m: 0.25,
                velocity_m_s: 0.25 * rate,
            };
            let response = linear_sdof_response_from_state(
                &[0.0, 0.0], duration, omega, zeta, initial, controls(), || false,
            )
            .unwrap();
            let displacement = initial.displacement_m * (rate * duration).exp();
            close(response.final_displacement_m, displacement, 2.0e-14);
            close(response.final_velocity_m_s, rate * displacement, 2.0e-13);
            close(response.peak_displacement_m, initial.displacement_m, 2.0e-14);
            close(response.signed_peak_velocity_m_s, initial.velocity_m_s, 2.0e-13);
            close(
                response.signed_peak_absolute_acceleration_m_s2,
                rate * rate * initial.displacement_m,
                2.0e-12,
            );
            assert_eq!(response.peak_time_s, 0.0);
            assert_eq!(response.peak_velocity_time_s, 0.0);
            assert_eq!(response.peak_absolute_acceleration_time_s, 0.0);
        }
    }

    #[test]
    fn equilibrium_and_zero_duration_segments_observe_initial_state() {
        let initial = LinearSdofState {
            displacement_m: 0.25,
            velocity_m_s: 0.0,
        };
        let mut limited = controls();
        limited.maximum_acceleration_zeros_per_interval = 0;
        limited.maximum_transition_evaluations = 2;
        let equilibrium = linear_sdof_response_from_state(
            &[-4.0, -4.0, -4.0], 0.25, 4.0, 0.3, initial, limited, || false,
        )
        .unwrap();
        assert_eq!(equilibrium.final_state(), initial);
        assert_eq!(equilibrium.peak_displacement_m, initial.displacement_m);
        assert_eq!(equilibrium.peak_velocity_m_s, 0.0);
        assert_eq!(equilibrium.signed_peak_absolute_acceleration_m_s2, -4.0);
        assert_eq!(equilibrium.peak_time_s, 0.0);
        assert_eq!(equilibrium.peak_absolute_acceleration_time_s, 0.0);
        assert_eq!(equilibrium.transition_evaluations, 2);
        assert_eq!(equilibrium.duration_s, 0.5);

        let moving = LinearSdofState {
            displacement_m: -0.4,
            velocity_m_s: 0.7,
        };
        for samples in [&[][..], &[9.0][..]] {
            let response = linear_sdof_response_from_state(
                samples, 0.2, 3.0, 0.2, moving, controls(), || false,
            )
            .unwrap();
            assert_eq!(response.final_state(), moving);
            assert_eq!(response.signed_peak_displacement_m, moving.displacement_m);
            assert_eq!(response.signed_peak_velocity_m_s, moving.velocity_m_s);
            close(response.signed_peak_absolute_acceleration_m_s2, 2.76, 2.0e-15);
            assert_eq!(response.peak_time_s, 0.0);
            assert_eq!(response.peak_velocity_time_s, 0.0);
            assert_eq!(response.peak_absolute_acceleration_time_s, 0.0);
            assert_eq!(response.duration_s, 0.0);
            assert_eq!(response.intervals, 0);
            assert_eq!(response.transition_evaluations, 0);
        }
    }

    #[test]
    fn initial_state_admission_and_continuation_limits_are_checked() {
        for initial in [
            LinearSdofState { displacement_m: f64::NAN, velocity_m_s: 0.0 },
            LinearSdofState { displacement_m: 0.0, velocity_m_s: f64::INFINITY },
            LinearSdofState { displacement_m: f64::MAX, velocity_m_s: 0.0 },
        ] {
            assert!(matches!(
                linear_sdof_response_from_state(
                    &[], 0.25, 4.0, 0.1, initial, controls(), || false,
                ),
                Err(LinearSdofError::Invalid(_)),
            ));
        }
        let initial = LinearSdofState { displacement_m: 0.1, velocity_m_s: 2.0 };
        let samples = [0.0; 5];
        let mut limited = controls();
        limited.maximum_transition_evaluations = 1;
        assert!(matches!(
            linear_sdof_response_from_state(
                &samples, 0.25, 4.0, 0.1, initial, limited, || false,
            ),
            Err(LinearSdofError::Budget("transition evaluations")),
        ));
        let mut polls = 0;
        assert!(matches!(
            linear_sdof_response_from_state(
                &samples, 0.25, 4.0, 0.1, initial, controls(), || {
                    polls += 1;
                    polls > 10
                },
            ),
            Err(LinearSdofError::Cancelled),
        ));
        assert_eq!(initial, LinearSdofState { displacement_m: 0.1, velocity_m_s: 2.0 });
        assert_eq!(samples, [0.0; 5]);
    }

    #[test]
    fn held_force_intersample_peak_survives_identical_endpoints() {
        let omega = 4.0;
        let dt = core::f64::consts::TAU / omega;
        let response = solve(&[3.0, 3.0], dt, omega, 0.0);
        close(response.final_displacement_m, 0.0, 2.0e-14);
        close(response.final_velocity_m_s, 0.0, 2.0e-14);
        close(response.peak_displacement_m, 6.0 / (omega * omega), 2.0e-14);
        close(
            response.signed_peak_displacement_m,
            -6.0 / (omega * omega),
            2.0e-14,
        );
        close(response.peak_time_s, core::f64::consts::PI / omega, 1.0e-11);
        close(response.peak_velocity_m_s, 3.0 / omega, 2.0e-14);
        // The full period has two tied velocity magnitudes. Use a half period
        // with one unique extremum to check its sign and physical time.
        let half_period = solve(&[3.0, 3.0], 0.5 * dt, omega, 0.0);
        close(half_period.signed_peak_velocity_m_s, -3.0 / omega, 2.0e-14);
        close(half_period.peak_velocity_time_s, core::f64::consts::FRAC_PI_2 / omega, 1.0e-11);
        close(response.peak_absolute_acceleration_m_s2, 6.0, 2.0e-13);
        close(response.signed_peak_absolute_acceleration_m_s2, 6.0, 2.0e-13);
        close(response.peak_absolute_acceleration_time_s, response.peak_time_s, 1.0e-11);
        assert_eq!(response.peak_absolute_acceleration_localization_allowance_m_s2,
            omega * omega * response.peak_localization_allowance_m);
        assert!(response.peak_localization_allowance_m < 1.0e-18);
        assert!(response.transition_evaluations > 1);
        assert_eq!(response.intervals, 1);
        assert_eq!(response.duration_s, dt);
    }

    #[test]
    fn damped_overshoot_and_monotone_critical_overdamped_limits() {
        let (omega, zeta, force) = (3.0_f64, 0.2_f64, 2.0_f64);
        let wd = omega * (1.0 - zeta * zeta).sqrt();
        let t_peak = core::f64::consts::PI / wd;
        let response = solve(&[force, force], 2.0 * t_peak, omega, zeta);
        let expected = force / (omega * omega) * (1.0 + (-zeta * omega * t_peak).exp());
        close(response.peak_displacement_m, expected, 3.0e-14);
        close(response.peak_time_s, t_peak, 1.0e-11);
        let velocity_time = zeta.acos() / wd;
        let velocity = force / omega * (-zeta * omega * velocity_time).exp();
        let acceleration_time = (core::f64::consts::PI - 2.0 * zeta.asin()) / wd;
        let acceleration = force * (1.0 + (-zeta * omega * acceleration_time).exp());
        close(response.peak_velocity_m_s, velocity, 3.0e-14);
        close(response.signed_peak_velocity_m_s, -velocity, 3.0e-14);
        close(response.peak_velocity_time_s, velocity_time, 1.0e-11);
        close(response.peak_absolute_acceleration_m_s2, acceleration, 3.0e-13);
        close(response.signed_peak_absolute_acceleration_m_s2, acceleration, 3.0e-13);
        close(response.peak_absolute_acceleration_time_s, acceleration_time, 1.0e-11);
        assert!(response.peak_absolute_acceleration_m_s2 > omega * omega * response.peak_displacement_m);
        assert!(response.peak_velocity_m_s < omega * response.peak_displacement_m - 0.1);

        let time = 1.7;
        let critical = solve(&[force, force], time, omega, 1.0);
        let q = -force / (omega * omega) * (1.0 - (1.0 + omega * time) * (-omega * time).exp());
        let v = -force * time * (-omega * time).exp();
        close(critical.final_displacement_m, q, 2.0e-14);
        close(critical.final_velocity_m_s, v, 2.0e-14);
        close(critical.peak_displacement_m, q.abs(), 2.0e-14);
        assert_eq!(critical.peak_time_s, time);
        close(critical.peak_velocity_m_s, force / (omega * core::f64::consts::E), 2.0e-14);
        close(critical.signed_peak_velocity_m_s, -critical.peak_velocity_m_s, 0.0);
        close(critical.peak_velocity_time_s, 1.0 / omega, 1.0e-11);
        close(critical.peak_absolute_acceleration_m_s2, force * (1.0 + (-2.0_f64).exp()), 2.0e-13);
        close(critical.peak_absolute_acceleration_time_s, 2.0 / omega, 1.0e-11);

        let zeta = 2.0_f64;
        let rs = -omega / (zeta + (zeta * zeta - 1.0).sqrt());
        let rf = -omega * (zeta + (zeta * zeta - 1.0).sqrt());
        let h = (rf * (rs * time).exp() - rs * (rf * time).exp()) / (rf - rs);
        let q = -force / (omega * omega) * (1.0 - h);
        let v = force * ((rs * time).exp() - (rf * time).exp()) / (rf - rs);
        let overdamped = solve(&[force, force], time, omega, zeta);
        close(overdamped.final_displacement_m, q, 3.0e-14);
        close(overdamped.final_velocity_m_s, v, 3.0e-14);
        close(overdamped.peak_displacement_m, q.abs(), 3.0e-14);
        assert_eq!(overdamped.peak_time_s, time);
        let velocity_time = (rf / rs).ln() / (rs - rf);
        let velocity = force * ((rs * velocity_time).exp() - (rf * velocity_time).exp()) / (rf - rs);
        let acceleration_time = 2.0 * velocity_time;
        let es = (rs * acceleration_time).exp();
        let ef = (rf * acceleration_time).exp();
        let acceleration = force * (1.0 - (rs * es - rf * ef) / (rs - rf));
        close(overdamped.peak_velocity_m_s, velocity.abs(), 3.0e-14);
        close(overdamped.peak_velocity_time_s, velocity_time, 1.0e-11);
        close(overdamped.peak_absolute_acceleration_m_s2, acceleration, 3.0e-13);
        close(overdamped.peak_absolute_acceleration_time_s, acceleration_time, 1.0e-11);
    }

    #[test]
    fn ramp_matches_analytic_response_and_has_no_appended_tail() {
        let (omega, slope, time) = (7.0_f64, 1.3_f64, 2.7_f64);
        let response = solve(&[0.0, slope * time], time, omega, 0.0);
        let q = -slope / (omega * omega) * (time - (omega * time).sin() / omega);
        let v = -slope / (omega * omega) * (1.0 - (omega * time).cos());
        close(response.final_displacement_m, q, 5.0e-14);
        close(response.final_velocity_m_s, v, 5.0e-14);
        close(response.peak_displacement_m, q.abs(), 5.0e-14);
        assert_eq!(response.peak_time_s, time);
        assert_eq!(response.duration_s, time);
        for samples in [&[][..], &[9.0][..]] {
            let empty = solve(samples, 0.2, omega, 0.4);
            assert_eq!(empty.peak_displacement_m, 0.0);
            assert_eq!(empty.peak_velocity_m_s, 0.0);
            assert_eq!(empty.peak_absolute_acceleration_m_s2, 0.0);
            assert_eq!(empty.duration_s, 0.0);
            assert_eq!(empty.intervals, 0);
            assert_eq!(empty.transition_evaluations, 0);
        }
    }

    #[test]
    fn exact_linear_record_subdivision_preserves_endpoints_and_extrema() {
        let samples = [0.0, 3.0, -1.0, 2.0, 0.0];
        let subdivisions = 7;
        let mut refined = Vec::new();
        for segment in samples.windows(2) {
            for j in 0..subdivisions {
                let fraction = f64::from(j) / f64::from(subdivisions);
                refined.push((1.0 - fraction) * segment[0] + fraction * segment[1]);
            }
        }
        refined.push(samples[4]);
        for zeta in [0.0, 0.08, 1.0 - 1.0e-10, 1.0, 1.0 + 1.0e-10, 3.0] {
            let coarse = solve(&samples, 0.9, 10.0, zeta);
            let fine = solve(&refined, 0.9 / f64::from(subdivisions), 10.0, zeta);
            close(
                coarse.final_displacement_m,
                fine.final_displacement_m,
                3.0e-13,
            );
            close(coarse.final_velocity_m_s, fine.final_velocity_m_s, 3.0e-12);
            close(
                coarse.peak_displacement_m,
                fine.peak_displacement_m,
                3.0e-13,
            );
            close(coarse.peak_time_s, fine.peak_time_s, 3.0e-10);
            close(coarse.peak_velocity_m_s, fine.peak_velocity_m_s, 3.0e-12);
            close(coarse.signed_peak_velocity_m_s, fine.signed_peak_velocity_m_s, 3.0e-12);
            close(coarse.peak_velocity_time_s, fine.peak_velocity_time_s, 3.0e-10);
            close(coarse.peak_absolute_acceleration_m_s2, fine.peak_absolute_acceleration_m_s2, 3.0e-11);
            close(coarse.signed_peak_absolute_acceleration_m_s2, fine.signed_peak_absolute_acceleration_m_s2, 3.0e-11);
            close(coarse.peak_absolute_acceleration_time_s, fine.peak_absolute_acceleration_time_s, 3.0e-10);
            close(coarse.duration_s, fine.duration_s, 1.0e-15);
            assert!(coarse.maximum_extremum_bracket_s <= controls().extremum_time_tolerance_s);
            assert!(coarse.maximum_absolute_acceleration_extremum_bracket_s
                <= controls().extremum_time_tolerance_s);
            assert_eq!(coarse, solve(&samples, 0.9, 10.0, zeta));
        }
    }

    // Independent physical-equation RK4 oracle; it does not use the exact kernel
    // or the acceleration-zero partition. Doubling its grid separately checks
    // the numerical reference, including critical and overdamped turning points.
    fn rk4(
        samples: &[f64],
        dt: f64,
        omega: f64,
        zeta: f64,
        subdivisions: usize,
    ) -> (f64, f64, f64, f64, f64) {
        let (mut q, mut v, mut peak) = (0.0_f64, 0.0_f64, 0.0_f64);
        let (mut velocity_peak, mut acceleration_peak) = (0.0_f64, 0.0_f64);
        let h = dt / subdivisions as f64;
        for segment in samples.windows(2) {
            let acceleration = |t: f64, x: f64, velocity: f64| {
                let ground = segment[0] + (segment[1] - segment[0]) * (t / dt);
                -ground - 2.0 * zeta * omega * velocity - omega * omega * x
            };
            for j in 0..subdivisions {
                let t = j as f64 * h;
                let (k1q, k1v) = (v, acceleration(t, q, v));
                let (k2q, k2v) = (
                    v + h * k1v / 2.0,
                    acceleration(t + h / 2.0, q + h * k1q / 2.0, v + h * k1v / 2.0),
                );
                let (k3q, k3v) = (
                    v + h * k2v / 2.0,
                    acceleration(t + h / 2.0, q + h * k2q / 2.0, v + h * k2v / 2.0),
                );
                let (k4q, k4v) = (v + h * k3v, acceleration(t + h, q + h * k3q, v + h * k3v));
                q += h / 6.0 * (k1q + 2.0 * k2q + 2.0 * k3q + k4q);
                v += h / 6.0 * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
                peak = peak.max(q.abs());
                velocity_peak = velocity_peak.max(v.abs());
                acceleration_peak = acceleration_peak.max((-2.0 * zeta * omega * v - omega * omega * q).abs());
            }
        }
        (q, v, peak, velocity_peak, acceleration_peak)
    }

    #[test]
    fn changing_forcing_matches_independent_equations_in_all_damping_regimes() {
        let samples = [2.0, -1.0, -3.0, 0.0, 2.0];
        for zeta in [0.0, 0.15, 1.0, 2.0, 20.0] {
            let response = solve(&samples, 0.7, 4.0, zeta);
            let coarse = rk4(&samples, 0.7, 4.0, zeta, 4096);
            let fine = rk4(&samples, 0.7, 4.0, zeta, 8192);
            close(coarse.0, fine.0, 2.0e-12);
            close(coarse.1, fine.1, 2.0e-12);
            close(response.final_displacement_m, fine.0, 2.0e-12);
            close(response.final_velocity_m_s, fine.1, 2.0e-12);
            // The independent grid underestimates a true peak by O(h²).
            close(response.peak_displacement_m, fine.2, 5.0e-8);
            assert!(response.peak_displacement_m + 2.0e-13 >= fine.2);
            close(coarse.3, fine.3, 2.0e-7);
            close(coarse.4, fine.4, 2.0e-6);
            close(response.peak_velocity_m_s, fine.3, 2.0e-7);
            close(response.peak_absolute_acceleration_m_s2, fine.4, 2.0e-6);
            assert!(response.peak_velocity_m_s + 2.0e-12 >= fine.3);
            assert!(response.peak_absolute_acceleration_m_s2 + 2.0e-11 >= fine.4);
        }
    }

    #[test]
    fn malformed_budget_and_cancelled_requests_never_publish_partial_peaks() {
        let samples = [1.0, 1.0];
        for (dt, omega, zeta) in [
            (0.0, 2.0, 0.1),
            (1.0, 0.0, 0.1),
            (1.0, 2.0, -0.1),
            (f64::NAN, 2.0, 0.1),
            (1.0, f64::INFINITY, 0.1),
        ] {
            assert!(linear_sdof_response(&samples, dt, omega, zeta, controls(), || false).is_err());
        }
        assert!(
            linear_sdof_response(&[0.0, f64::NAN], 1.0, 2.0, 0.1, controls(), || false).is_err()
        );
        let mut limited = controls();
        limited.maximum_samples = 1;
        assert!(matches!(
            linear_sdof_response(&samples, 2.0, 4.0, 0.0, limited, || false),
            Err(LinearSdofError::Budget("acceleration samples"))
        ));
        limited = controls();
        limited.maximum_acceleration_zeros_per_interval = 0;
        assert!(matches!(
            linear_sdof_response(&samples, 2.0, 4.0, 0.0, limited, || false),
            Err(LinearSdofError::Budget(
                "acceleration zeros per input interval"
            ))
        ));
        limited = controls();
        limited.maximum_transition_evaluations = 1;
        assert!(matches!(
            linear_sdof_response(&samples, 2.0, 4.0, 0.0, limited, || false),
            Err(LinearSdofError::Budget("transition evaluations"))
        ));
        limited = controls();
        limited.maximum_bisections_per_extremum = 1;
        assert!(matches!(
            linear_sdof_response(&samples, 2.0, 4.0, 0.2, limited, || false),
            Err(LinearSdofError::Budget("extremum bisections"))
        ));
        assert!(matches!(
            linear_sdof_response(&samples, 2.0, 4.0, 0.0, controls(), || true),
            Err(LinearSdofError::Cancelled)
        ));
        let mut polls = 0;
        assert!(matches!(
            linear_sdof_response(&samples, 2.0, 4.0, 0.2, controls(), || {
                polls += 1;
                polls > 12
            }),
            Err(LinearSdofError::Cancelled)
        ));
        assert_eq!(samples, [1.0, 1.0]);
    }

    #[test]
    fn absolute_acceleration_search_obeys_its_partition_and_shared_work_limits() {
        let samples = [1.0, 1.0];
        let mut limited = controls();
        // This interval contains one acceleration zero but two curvature zeros.
        limited.maximum_acceleration_zeros_per_interval = 1;
        assert!(matches!(
            linear_sdof_response(&samples, 1.1, 4.0, 0.2, limited, || false),
            Err(LinearSdofError::Budget("absolute-acceleration curvature zeros per input interval"))
        ));
        let full = solve(&samples, 1.1, 4.0, 0.2);
        limited = controls();
        limited.maximum_transition_evaluations = full.transition_evaluations - 1;
        assert!(matches!(
            linear_sdof_response(&samples, 1.1, 4.0, 0.2, limited, || false),
            Err(LinearSdofError::Budget("transition evaluations"))
        ));
        assert!(linear_sdof_response(&samples, 1e-110, 1e110, 1.0, controls(), || false).is_err());
    }
}
