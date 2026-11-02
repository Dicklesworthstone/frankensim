//! Exact linear response to the piecewise-linear acceleration actually supplied.
//!
//! On each interval acceleration satisfies the homogeneous oscillator equation.
//! Its analytic zeros partition velocity into monotone pieces; bracketing velocity
//! zeros on those pieces finds displacement extrema without a sampling grid.
//! This resolves intersample peaks, not unmeasured ground-motion frequencies.

use fs_math::det;
use fs_modal::time::{LinearForceTransition, ModalAcousticState, ModalAcousticTimeError};

/// Explicit work limits and numerical extremum-localization tolerance.
#[derive(Clone, Copy, Debug)]
pub struct LinearSdofControls {
    /// Maximum supplied acceleration samples, checked before scanning the record.
    pub maximum_samples: usize,
    /// Maximum acceleration zeros in any one input interval.
    pub maximum_acceleration_zeros_per_interval: usize,
    /// Total full-interval and search transition evaluations.
    pub maximum_transition_evaluations: usize,
    /// Bisection ceiling for each displacement extremum.
    pub maximum_bisections_per_extremum: usize,
    /// Maximum width [s] of a velocity-zero bracket before evaluating its midpoint.
    pub extremum_time_tolerance_s: f64,
}

/// Response from initial rest over exactly the supplied record duration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinearSdofResponse {
    /// Maximum absolute relative displacement [m], including intersample extrema.
    pub peak_displacement_m: f64,
    /// Signed displacement [m] at the selected extremum or record endpoint.
    pub signed_peak_displacement_m: f64,
    /// Physical time [s] of the selected evaluated peak; ties select the first.
    pub peak_time_s: f64,
    /// Largest final velocity-zero bracket width [s].
    pub maximum_extremum_bracket_s: f64,
    /// Maximum localization allowance [m] from bracket width and endpoint speeds.
    /// This excludes floating-point transition/root errors and input/model errors;
    /// it is not a certified enclosure of the complete physical response.
    pub peak_localization_allowance_m: f64,
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

    fn displacement_extremum(
        &mut self,
        mut low: f64,
        mut high: f64,
        mut low_state: ModalAcousticState,
        mut high_state: ModalAcousticState,
    ) -> Result<(f64, ModalAcousticState, f64, f64), LinearSdofError> {
        for _ in 0..self.work.controls.maximum_bisections_per_extremum {
            self.work.poll()?;
            let width = high - low;
            let middle = low + 0.5 * width;
            if width <= self.work.controls.extremum_time_tolerance_s {
                let allowance = 0.5
                    * width
                    * low_state
                        .velocity_m_sqrt_kg_per_s
                        .abs()
                        .max(high_state.velocity_m_sqrt_kg_per_s.abs());
                return Ok((middle, self.at(middle)?, width, finite(allowance)?));
            }
            if middle == low || middle == high {
                return Err(LinearSdofError::Invalid(
                    "extremum tolerance is below representable time resolution",
                ));
            }
            let state = self.at(middle)?;
            let velocity = state.velocity_m_sqrt_kg_per_s;
            if velocity == 0.0 {
                return Ok((middle, state, 0.0, 0.0));
            }
            if velocity.is_sign_positive() == low_state.velocity_m_sqrt_kg_per_s.is_sign_positive()
            {
                low = middle;
                low_state = state;
            } else {
                high = middle;
                high_state = state;
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

// Returns the first strictly positive acceleration zero and its spacing.
// Critical/overdamped acceleration has at most one zero (spacing = infinity).
fn acceleration_zeros(
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
    let a = finite(f0 - 2.0 * alpha * velocity - stiffness * state.displacement_m_sqrt_kg)?;
    let slope = finite((f1 - f0) / dt)?;
    // a'(0) + alpha*a(0); compute directly without cancelling 2*alpha*a.
    let d = finite(slope - alpha * a - stiffness * velocity)?;
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
    fn observe(&mut self, state: ModalAcousticState, time: f64) {
        let value = state.displacement_m_sqrt_kg;
        if value.abs() > self.peak_displacement_m
            || (value.abs() == self.peak_displacement_m && time < self.peak_time_s)
        {
            self.peak_displacement_m = value.abs();
            self.signed_peak_displacement_m = value;
            self.peak_time_s = time;
        }
    }
}

/// Solve `u'' + 2*zeta*omega*u' + omega²*u = -ground_acceleration` from rest.
///
/// Uses the shared exact first-order-hold oscillator transition, cached for full
/// input intervals. The acceleration samples are endpoints of straight segments:
/// N samples cover N-1 intervals, with no implicit decay tail. Empty or single-
/// sample records return zero motion. Units are seconds, rad/s and m/s².
///
/// Every displacement extremum is sought on analytic monotone-velocity pieces.
/// The tolerance controls root localization, not input bandwidth or model error.
/// Work is bounded and storage is constant, independent of the record length.
/// Deterministic repeatability is on the same ISA/toolchain; no interval-certified
/// extrema or cross-ISA claim is made. The legacy Newmark APIs remain separate.
///
/// # Errors
/// Refuses invalid inputs, unrepresentable extrema arithmetic, exhausted work
/// limits, cancellation, or an exact-transition refusal. Returns no partial result.
#[allow(clippy::cast_precision_loss)]
pub fn linear_sdof_response<C: FnMut() -> bool>(
    ground_acceleration_m_s2: &[f64],
    dt_s: f64,
    omega_rad_s: f64,
    damping_ratio: f64,
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
    finite(omega_rad_s * omega_rad_s)?;
    finite(2.0 * omega_rad_s * damping_ratio)?;
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
    let mut state = ModalAcousticState::default();
    let mut report = LinearSdofResponse {
        peak_displacement_m: 0.0,
        signed_peak_displacement_m: 0.0,
        peak_time_s: 0.0,
        maximum_extremum_bracket_s: 0.0,
        peak_localization_allowance_m: 0.0,
        final_displacement_m: 0.0,
        final_velocity_m_s: 0.0,
        duration_s,
        intervals,
        transition_evaluations: 0,
    };
    for (index, samples) in ground_acceleration_m_s2.windows(2).enumerate() {
        work.evaluation()?;
        let (f0, f1) = (-samples[0], -samples[1]);
        let end = full.advance(state, f0, f1)?;
        let (first_zero, spacing) =
            acceleration_zeros(state, f0, f1, dt_s, omega_rad_s, damping_ratio)?;
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
        let (mut low, mut low_state) = (0.0, state);
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
                    "acceleration-zero spacing is below representable time resolution",
                ));
            }
            if high < dt_s && zero_index >= controls.maximum_acceleration_zeros_per_interval {
                return Err(LinearSdofError::Budget(
                    "acceleration zeros per input interval",
                ));
            }
            let high_state = interval.at(high)?;
            report.observe(high_state, base_time + high);
            let (vl, vh) = (
                low_state.velocity_m_sqrt_kg_per_s,
                high_state.velocity_m_sqrt_kg_per_s,
            );
            if vl != 0.0 && vh != 0.0 && vl.is_sign_positive() != vh.is_sign_positive() {
                let (time, peak, width, allowance) =
                    interval.displacement_extremum(low, high, low_state, high_state)?;
                report.observe(peak, base_time + time);
                report.maximum_extremum_bracket_s = report.maximum_extremum_bracket_s.max(width);
                report.peak_localization_allowance_m =
                    report.peak_localization_allowance_m.max(allowance);
            }
            if high == dt_s {
                break;
            }
            low = high;
            low_state = high_state;
            zero_index += 1;
        }
        state = end;
    }
    work.poll()?;
    report.final_displacement_m = state.displacement_m_sqrt_kg;
    report.final_velocity_m_s = state.velocity_m_sqrt_kg_per_s;
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

        let time = 1.7;
        let critical = solve(&[force, force], time, omega, 1.0);
        let q = -force / (omega * omega) * (1.0 - (1.0 + omega * time) * (-omega * time).exp());
        let v = -force * time * (-omega * time).exp();
        close(critical.final_displacement_m, q, 2.0e-14);
        close(critical.final_velocity_m_s, v, 2.0e-14);
        close(critical.peak_displacement_m, q.abs(), 2.0e-14);
        assert_eq!(critical.peak_time_s, time);

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
            close(coarse.duration_s, fine.duration_s, 1.0e-15);
            assert!(coarse.maximum_extremum_bracket_s <= controls().extremum_time_tolerance_s);
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
    ) -> (f64, f64, f64) {
        let (mut q, mut v, mut peak) = (0.0_f64, 0.0_f64, 0.0_f64);
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
            }
        }
        (q, v, peak)
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
}
