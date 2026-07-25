//! Duty cycles: a declared power-versus-time schedule driving the transient
//! march, with a windowed energy audit.
//!
//! Bead `frankensim-extreal-program-f85xj.5.12`, the second item staged by
//! `f85xj.5.9`. A power map is a snapshot; hardware runs a cycle, and the
//! questions that make transient analysis worth doing — peak junction
//! temperature, time above a limit, excursion amplitude — are statements
//! about a HISTORY, not a steady value.
//!
//! # The schedule is a scale on an existing power map
//!
//! A [`DutyCycle`] is a piecewise profile of a dimensionless multiplier on the
//! volumetric source the caller already validated through
//! [`crate::power::PowerMap`]. Keeping it a scale rather than a second
//! description of the hardware means the cycle cannot disagree with the map
//! about which components exist or what they dissipate at full load — there
//! is one declaration of the power distribution and one of its time profile.
//!
//! # Segments must TILE the window
//!
//! Segments are contiguous by construction: each declares a duration, and the
//! window is their sum. Marching past the end of the declared window REFUSES
//! rather than holding the last value. Holding is the tempting default and it
//! is wrong — a schedule that ends is a statement that the model has nothing
//! to say afterwards, and silently extending it invents load history that
//! nobody declared.
//!
//! # The windowed energy audit, and when it is exact
//!
//! The theta method injects, over one step, the load weighted
//! `θ s^{n+1} + (1−θ) s^n`. Summed over the window that is a quadrature of
//! `∫ s dt`:
//!
//! - at `theta = 0.5` it is exactly the TRAPEZOID rule, which is EXACT for a
//!   piecewise-linear profile whenever step boundaries fall on segment
//!   boundaries. So Crank-Nicolson with aligned steps delivers exactly the
//!   declared cycle energy, to floating point.
//! - at `theta = 1` it is the right-endpoint rule, first order, and the
//!   residual is real rather than roundoff.
//!
//! [`WindowedEnergyAudit`] reports declared and delivered energy and the
//! residual between them; it does not silently normalize one to the other,
//! because the gap is exactly the integration error the caller chose by
//! picking a scheme and a step.

use crate::ConductionError;
use crate::assemble::{DofMap, assemble_operator, reduce_matrix_and_lift};
use crate::transient::{TransientConfig, TransientProblem, assemble_capacitance};
use fs_exec::Cx;
use fs_solver::CsrOp;
use fs_solver::krylov::CgState;
use fs_sparse::Csr;

/// Maximum segments admitted in one declared cycle.
pub const MAX_SEGMENTS: usize = 4_096;

/// How the scale varies across one segment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SegmentInterpolation {
    /// Held at the start value for the whole segment. A step load.
    Constant,
    /// Linear from the start value to the end value.
    Linear,
}

/// One segment of a declared duty cycle.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DutySegment {
    duration_s: f64,
    start_scale: f64,
    end_scale: f64,
    interpolation: SegmentInterpolation,
}

impl DutySegment {
    /// A segment held at one scale.
    ///
    /// # Errors
    /// A non-positive or non-finite duration, or a negative/non-finite scale.
    pub fn constant(duration_s: f64, scale: f64) -> Result<Self, ConductionError> {
        Self::admit(duration_s, scale, scale, SegmentInterpolation::Constant)
    }

    /// A segment ramping linearly between two scales.
    ///
    /// # Errors
    /// As [`DutySegment::constant`].
    pub fn ramp(
        duration_s: f64,
        start_scale: f64,
        end_scale: f64,
    ) -> Result<Self, ConductionError> {
        Self::admit(
            duration_s,
            start_scale,
            end_scale,
            SegmentInterpolation::Linear,
        )
    }

    fn admit(
        duration_s: f64,
        start_scale: f64,
        end_scale: f64,
        interpolation: SegmentInterpolation,
    ) -> Result<Self, ConductionError> {
        if !duration_s.is_finite() || duration_s <= 0.0 {
            return Err(duty_error(
                "segment duration",
                format!("{duration_s} s is not finite and positive"),
            ));
        }
        for (value, field) in [(start_scale, "start scale"), (end_scale, "end scale")] {
            if !value.is_finite() || value < 0.0 {
                return Err(duty_error(
                    field,
                    format!(
                        "{value} is not finite and non-negative; a negative scale would make the heat source a sink, which is a boundary condition rather than a duty state"
                    ),
                ));
            }
        }
        Ok(Self {
            duration_s,
            start_scale,
            end_scale,
            interpolation,
        })
    }

    /// Segment duration, s.
    #[must_use]
    pub const fn duration_s(self) -> f64 {
        self.duration_s
    }

    /// Scale at the segment start.
    #[must_use]
    pub const fn start_scale(self) -> f64 {
        self.start_scale
    }

    /// Scale at the segment end.
    #[must_use]
    pub const fn end_scale(self) -> f64 {
        self.end_scale
    }

    /// Declared interpolation.
    #[must_use]
    pub const fn interpolation(self) -> SegmentInterpolation {
        self.interpolation
    }

    /// Scale at a local offset into this segment.
    fn scale_at_offset(self, offset_s: f64) -> f64 {
        match self.interpolation {
            SegmentInterpolation::Constant => self.start_scale,
            SegmentInterpolation::Linear => {
                let fraction = (offset_s / self.duration_s).clamp(0.0, 1.0);
                self.end_scale
                    .mul_add(fraction, self.start_scale * (1.0 - fraction))
            }
        }
    }

    /// Exact `∫ scale dt` over the whole segment.
    fn energy_scale(self) -> f64 {
        match self.interpolation {
            SegmentInterpolation::Constant => self.start_scale * self.duration_s,
            SegmentInterpolation::Linear => {
                f64::midpoint(self.start_scale, self.end_scale) * self.duration_s
            }
        }
    }
}

/// A declared power-versus-time profile over a bounded window.
#[derive(Debug, Clone, PartialEq)]
pub struct DutyCycle {
    segments: Vec<DutySegment>,
    boundaries_s: Vec<f64>,
    window_s: f64,
}

impl DutyCycle {
    /// Admit a cycle from contiguous segments.
    ///
    /// # Errors
    /// An empty or oversized segment list, or a non-finite accumulated
    /// window.
    pub fn new(segments: Vec<DutySegment>) -> Result<Self, ConductionError> {
        if segments.is_empty() {
            return Err(duty_error(
                "duty cycle",
                "declares no segments; a cycle with no history cannot drive a march".to_string(),
            ));
        }
        if segments.len() > MAX_SEGMENTS {
            return Err(duty_error(
                "duty cycle",
                format!(
                    "declares {} segments, above the admitted maximum {MAX_SEGMENTS}",
                    segments.len()
                ),
            ));
        }
        let mut boundaries_s = Vec::with_capacity(segments.len() + 1);
        let mut accumulated = 0.0f64;
        boundaries_s.push(0.0);
        for segment in &segments {
            accumulated += segment.duration_s;
            if !accumulated.is_finite() {
                return Err(duty_error(
                    "duty cycle window",
                    "accumulated duration left the finite range".to_string(),
                ));
            }
            boundaries_s.push(accumulated);
        }
        Ok(Self {
            segments,
            boundaries_s,
            window_s: accumulated,
        })
    }

    /// Declared segments, in order.
    #[must_use]
    pub fn segments(&self) -> &[DutySegment] {
        &self.segments
    }

    /// Segment boundary times, `len() + 1` entries starting at zero.
    #[must_use]
    pub fn boundaries_s(&self) -> &[f64] {
        &self.boundaries_s
    }

    /// Total declared window, s.
    #[must_use]
    pub const fn window_s(&self) -> f64 {
        self.window_s
    }

    /// Scale at a time inside the declared window.
    ///
    /// # Errors
    /// A non-finite, negative, or past-the-window time. The cycle does NOT
    /// hold its last value: a schedule that ends is a statement that the
    /// model has nothing to say afterwards.
    pub fn scale_at(&self, time_s: f64) -> Result<f64, ConductionError> {
        if !time_s.is_finite() || time_s < 0.0 {
            return Err(duty_error(
                "duty cycle time",
                format!("{time_s} s is not finite and non-negative"),
            ));
        }
        if time_s > self.window_s * (1.0 + 1.0e-12) {
            return Err(duty_error(
                "duty cycle time",
                format!(
                    "{time_s} s is past the declared {} s window; the cycle does not hold its last value, because extending a finished schedule invents load history nobody declared",
                    self.window_s
                ),
            ));
        }
        let clamped = time_s.min(self.window_s);
        for (index, segment) in self.segments.iter().enumerate() {
            let end = self.boundaries_s[index + 1];
            if clamped <= end || index + 1 == self.segments.len() {
                return Ok(segment.scale_at_offset(clamped - self.boundaries_s[index]));
            }
        }
        unreachable!("boundaries cover the window")
    }

    /// Exact `∫ scale dt` over the whole declared window.
    ///
    /// Analytic per segment, so it is the reference the quadrature in
    /// [`WindowedEnergyAudit`] is judged against rather than another
    /// approximation of the same integral.
    #[must_use]
    pub fn energy_scale_seconds(&self) -> f64 {
        self.segments
            .iter()
            .fold(0.0f64, |acc, segment| acc + segment.energy_scale())
    }

    /// Whether every boundary falls on a multiple of `dt`, within tolerance.
    ///
    /// Alignment is what makes the Crank-Nicolson energy audit exact: the
    /// trapezoid rule reproduces a piecewise-linear profile only when no step
    /// straddles a breakpoint.
    #[must_use]
    pub fn steps_align(&self, dt_s: f64, tolerance: f64) -> bool {
        self.boundaries_s.iter().all(|boundary| {
            let ratio = boundary / dt_s;
            (ratio - ratio.round()).abs() <= tolerance
        })
    }
}

/// Declared-versus-delivered energy over one marched window.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WindowedEnergyAudit {
    declared_j: f64,
    delivered_j: f64,
    window_s: f64,
    aligned: bool,
}

impl WindowedEnergyAudit {
    /// `∫ P dt` from the declared schedule, J.
    #[must_use]
    pub const fn declared_j(self) -> f64 {
        self.declared_j
    }

    /// The energy the marched steps actually injected, J.
    #[must_use]
    pub const fn delivered_j(self) -> f64 {
        self.delivered_j
    }

    /// Signed `delivered − declared`, J.
    ///
    /// NOT normalized away: the gap is the integration error the caller chose
    /// by picking a scheme and a step, and reporting it is how that choice
    /// stays visible.
    #[must_use]
    pub fn residual_j(self) -> f64 {
        self.delivered_j - self.declared_j
    }

    /// Marched window, s.
    #[must_use]
    pub const fn window_s(self) -> f64 {
        self.window_s
    }

    /// Whether every segment boundary fell on a step boundary.
    #[must_use]
    pub const fn steps_aligned(self) -> bool {
        self.aligned
    }
}

/// Cycle-level quantities extracted from a marched window.
#[derive(Debug, Clone, PartialEq)]
pub struct CycleSummary {
    peak_temperature_k: f64,
    peak_time_s: f64,
    final_temperature_k: f64,
    excursion_k: f64,
    steps_above_limit: usize,
    time_above_limit_s: f64,
}

impl CycleSummary {
    /// Highest nodal temperature seen at any recorded step end, K.
    ///
    /// A step-sampled maximum: the true continuous peak can fall BETWEEN
    /// steps and is not resolved here. Refining the step tightens this;
    /// nothing in the record claims the sampled peak is the continuous one.
    #[must_use]
    pub const fn peak_temperature_k(&self) -> f64 {
        self.peak_temperature_k
    }

    /// Time at which the sampled peak occurred, s.
    #[must_use]
    pub const fn peak_time_s(&self) -> f64 {
        self.peak_time_s
    }

    /// Highest nodal temperature at the end of the window, K.
    #[must_use]
    pub const fn final_temperature_k(&self) -> f64 {
        self.final_temperature_k
    }

    /// Sampled peak minus the initial maximum, K.
    ///
    /// An amplitude, NOT a fatigue quantity. Cycle-counting, damage models,
    /// and acceleration factors are explicitly out of scope.
    #[must_use]
    pub const fn excursion_k(&self) -> f64 {
        self.excursion_k
    }

    /// Recorded steps whose maximum exceeded the declared limit.
    #[must_use]
    pub const fn steps_above_limit(&self) -> usize {
        self.steps_above_limit
    }

    /// Step-sampled time above the declared limit, s.
    ///
    /// Counted per whole step, so it is quantized to the step size and is an
    /// over- or under-estimate by up to one step at each crossing.
    #[must_use]
    pub const fn time_above_limit_s(&self) -> f64 {
        self.time_above_limit_s
    }
}

/// A completed duty-cycle march.
#[derive(Debug, Clone, PartialEq)]
pub struct DutyCycleSolution {
    /// Nodal temperature at the end of the window, K.
    pub temperature: Vec<f64>,
    /// The windowed energy audit.
    pub energy: WindowedEnergyAudit,
    /// Cycle-level summary quantities.
    pub summary: CycleSummary,
}

/// March a declared duty cycle over its whole window.
///
/// `problem.source` is the FULL-LOAD volumetric source — typically the field
/// [`crate::power::PowerMap::volumetric_source`] produced — and the cycle
/// scales it. `base_power_w` is that map's total delivered power, used only
/// to convert the schedule's dimensionless integral into joules for the
/// audit.
///
/// The march covers exactly `ceil(window / dt)` steps and the window must be
/// an exact multiple of the step, so no step straddles the end.
///
/// # Errors
/// [`ConductionError::ScenarioRow`] for a window that is not a whole number
/// of steps or a non-finite base power, [`ConductionError::Config`] for a
/// temperature-dependent conductivity, [`ConductionError::FieldLength`] for a
/// mismatched initial vector, [`ConductionError::LinearSolveFailed`] on a
/// non-converging step, and [`ConductionError::Cancelled`] at a step
/// boundary.
pub fn march_duty_cycle(
    cx: &Cx<'_>,
    problem: TransientProblem<'_>,
    config: &TransientConfig,
    cycle: &DutyCycle,
    base_power_w: f64,
    initial: &[f64],
    limit_k: f64,
) -> Result<DutyCycleSolution, ConductionError> {
    let TransientProblem {
        mesh,
        boundary,
        material,
        source,
        capacity,
    } = problem;
    if material.is_temperature_dependent() {
        return Err(ConductionError::Config {
            parameter: "conductivity",
            what: "duty-cycle marching received a temperature-dependent conductivity; the transient path admits constant conductivity only".to_string(),
        });
    }
    if !base_power_w.is_finite() || base_power_w < 0.0 {
        return Err(duty_error(
            "base power",
            format!("{base_power_w} W is not finite and non-negative"),
        ));
    }
    let n = mesh.vertex_count();
    if initial.len() != n {
        return Err(ConductionError::FieldLength {
            field: "initial temperature",
            expected: n,
            found: initial.len(),
        });
    }

    let dt = config.dt_s();
    let steps = admit_whole_steps(cycle, dt)?;

    let dofs = DofMap::new(boundary, n)?;
    let capacitance = assemble_capacitance(cx, mesh, capacity)?;
    let system = assemble_operator(cx, mesh, boundary, material, source, initial)?;

    let inverse_dt = 1.0 / dt;
    let theta = config.theta();
    let lhs = scaled_sum(&capacitance, inverse_dt, &system.operator, theta);
    let (a_ff, lift) = reduce_matrix_and_lift(&lhs, &dofs);
    let precond = crate::solve::spd_preconditioner(&a_ff);

    let mut temperature = initial.to_vec();
    for &vertex in dofs.fixed() {
        temperature[vertex] = dofs.prescribed()[vertex];
    }

    let initial_peak = temperature.iter().fold(f64::NEG_INFINITY, |a, b| a.max(*b));
    let mut peak = initial_peak;
    let mut peak_time = 0.0f64;
    let mut steps_above = 0usize;
    let mut delivered_scale_seconds = 0.0f64;

    let mut scratch_c = vec![0.0f64; n];
    let mut scratch_k = vec![0.0f64; n];

    for step in 0..steps {
        cx.checkpoint().map_err(|_| ConductionError::Cancelled {
            stage: "duty-cycle-step",
            at: step,
        })?;
        #[allow(clippy::cast_precision_loss)]
        let t_old = dt * (step as f64);
        #[allow(clippy::cast_precision_loss)]
        let t_new = dt * ((step + 1) as f64);
        let scale_old = cycle.scale_at(t_old)?;
        let scale_new = cycle.scale_at(t_new)?;
        // The theta method's own load weighting: at theta = 0.5 this is the
        // trapezoid rule, which is exact for a piecewise-linear profile when
        // no step straddles a segment boundary.
        let weight = theta.mul_add(scale_new, (1.0 - theta) * scale_old);
        delivered_scale_seconds += weight * dt;

        capacitance.spmv(&temperature, &mut scratch_c);
        system.operator.spmv(&temperature, &mut scratch_k);
        let mut rhs = Vec::with_capacity(dofs.n());
        for (slot, &vertex) in dofs.free().iter().enumerate() {
            let value = scratch_c[vertex].mul_add(
                inverse_dt,
                (1.0 - theta).mul_add(-scratch_k[vertex], weight * system.load[vertex]),
            ) + lift[slot];
            if !value.is_finite() {
                return Err(ConductionError::NonFinite {
                    field: "duty-cycle right-hand side",
                    bits: value.to_bits(),
                });
            }
            rhs.push(value);
        }

        let op = CsrOp::symmetric(a_ff.clone());
        let mut cg = CgState::new(&op, &precond, &rhs);
        let report = cg.run(
            &op,
            &precond,
            config.linear_tolerance(),
            config.linear_max_iterations(),
        );
        let truth = relative_residual(&a_ff, &rhs, &cg.x);
        if truth.is_nan() || truth >= config.linear_tolerance() {
            return Err(ConductionError::LinearSolveFailed {
                iteration: step,
                krylov_iterations: report.iters,
                true_relative_residual: truth,
                tolerance: config.linear_tolerance(),
            });
        }
        for (slot, &vertex) in dofs.free().iter().enumerate() {
            temperature[vertex] = cg.x[slot];
        }

        let step_peak = temperature.iter().fold(f64::NEG_INFINITY, |a, b| a.max(*b));
        if step_peak > peak {
            peak = step_peak;
            peak_time = t_new;
        }
        if step_peak > limit_k {
            steps_above += 1;
        }
    }

    cx.checkpoint().map_err(|_| ConductionError::Cancelled {
        stage: "duty-cycle-publish",
        at: steps,
    })?;

    let energy = WindowedEnergyAudit {
        declared_j: cycle.energy_scale_seconds() * base_power_w,
        delivered_j: delivered_scale_seconds * base_power_w,
        window_s: cycle.window_s(),
        aligned: cycle.steps_align(dt, 1.0e-9),
    };
    #[allow(clippy::cast_precision_loss)]
    let time_above_limit_s = dt * (steps_above as f64);
    let summary = CycleSummary {
        peak_temperature_k: peak,
        peak_time_s: peak_time,
        final_temperature_k: temperature.iter().fold(f64::NEG_INFINITY, |a, b| a.max(*b)),
        excursion_k: peak - initial_peak,
        steps_above_limit: steps_above,
        time_above_limit_s,
    };

    Ok(DutyCycleSolution {
        temperature,
        energy,
        summary,
    })
}

/// Refuse a window that is not a whole number of steps.
///
/// A step straddling the window end would inject load the schedule does not
/// declare, so the mismatch refuses rather than being truncated or extended.
fn admit_whole_steps(cycle: &DutyCycle, dt: f64) -> Result<usize, ConductionError> {
    let raw = cycle.window_s() / dt;
    let rounded = raw.round();
    if (raw - rounded).abs() > 1.0e-9 || rounded < 1.0 {
        return Err(duty_error(
            "duty cycle window",
            format!(
                "the {} s window is not a whole number of {dt} s steps ({raw} of them); a step straddling the window end would inject load the schedule does not declare",
                cycle.window_s()
            ),
        ));
    }
    #[allow(clippy::cast_sign_loss, clippy::cast_possible_truncation)]
    Ok(rounded as usize)
}

fn scaled_sum(a: &Csr, alpha: f64, b: &Csr, beta: f64) -> Csr {
    let n = a.nrows();
    let mut coo = fs_sparse::Coo::new(n, n);
    for row in 0..n {
        let (columns, values) = a.row(row);
        for (&column, &value) in columns.iter().zip(values.iter()) {
            coo.push(row, column, alpha * value);
        }
        let (columns, values) = b.row(row);
        for (&column, &value) in columns.iter().zip(values.iter()) {
            coo.push(row, column, beta * value);
        }
    }
    coo.assemble()
}

fn relative_residual(a: &Csr, b: &[f64], x: &[f64]) -> f64 {
    let mut ax = vec![0.0f64; b.len()];
    a.spmv(x, &mut ax);
    let residual: Vec<f64> = ax.iter().zip(b.iter()).map(|(v, t)| t - v).collect();
    let denominator = fs_solver::norm2(b);
    if denominator > 0.0 {
        fs_solver::norm2(&residual) / denominator
    } else {
        fs_solver::norm2(&residual)
    }
}

fn duty_error(field: &str, what: String) -> ConductionError {
    ConductionError::ScenarioRow {
        region: field.to_string(),
        what,
        fix: "correct the declared duty cycle; a schedule is a declaration, not a hint".to_string(),
    }
}
