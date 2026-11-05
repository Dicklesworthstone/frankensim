//! Embedded Bogacki--Shampine 3(2) integration of the original connected wake.
//! Four reusable vertex buffers; the accepted endpoint derivative is reused.

use super::super::{WakeAdaptiveReport, WakeTimeSpec};
use super::{
    ImagePlane, PreparedSources, Refusal, SelfInducedShed, SourceTree, StageSettings,
    WakeAdvectionReport, WakeTreeSpec, WorkMeter, admit_line_step, admit_point, allocation,
    checkpoint, displacement, refusal, translated, velocity_at,
};

#[allow(clippy::too_many_arguments)]
pub(crate) fn integrate_adaptive_nodes(
    points: &mut [[f64; 3]],
    config: SelfInducedShed,
    settings: StageSettings,
    tree_spec: Option<WakeTreeSpec>,
    time: WakeTimeSpec,
    line_displacement_m: f64,
    assemble: impl FnMut(&[[f64; 3]], f64, &mut PreparedSources) -> Result<(), Refusal>,
    variation: &mut impl FnMut([f64; 3], f64) -> Result<[f64; 3], Refusal>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<WakeAdaptiveReport, Refusal> {
    let dt = settings.dt_s;
    let max_step = dt / config.substeps as f64;
    if !time.absolute_position_tolerance_m.is_finite()
        || time.absolute_position_tolerance_m <= 0.0
        || !time.relative_position_tolerance.is_finite()
        || time.relative_position_tolerance < 0.0
        || !time.min_step_s.is_finite()
        || time.min_step_s <= 0.0
        || time.min_step_s > max_step
        || time.max_trials == 0
        || !max_step.is_finite()
        || max_step <= 0.0
    {
        return Err(refusal(
            "wake-advection-time-invalid",
            "finite positive absolute tolerance and minimum step, nonnegative relative tolerance, \
             and positive trial cap required; minimum step must not exceed the initial step",
        ));
    }
    if tree_spec.is_some_and(|s| {
        !s.absolute_velocity_error_mps.is_finite()
            || s.absolute_velocity_error_mps < 0.0
            || !(1..=256).contains(&s.leaf_capacity)
    }) {
        return Err(refusal(
            "wake-tree-invalid",
            "finite nonnegative accuracy and leaf size 1..=256 required",
        ));
    }
    if !line_displacement_m.is_finite()
        || line_displacement_m < 0.0
        || settings.image_plane.is_some_and(|p| !p.z0_m.is_finite())
    {
        return Err(refusal(
            "wake-advection-invalid",
            "invalid stage boundary or line displacement",
        ));
    }
    let bound = settings
        .source_bound
        .checked_mul(if settings.image_plane.is_some() { 2 } else { 1 })
        .ok_or_else(work_exhausted)?;
    let targets = u64::try_from(points.len()).map_err(|_| work_exhausted())?;
    let sources = u64::try_from(bound).map_err(|_| work_exhausted())?;
    if sources
        .checked_add(targets)
        .and_then(|n| n.checked_add(2))
        .is_none_or(|n| n > config.max_work)
    {
        return Err(work_exhausted());
    }
    let mut work = WorkMeter {
        limit: config.max_work,
        used: 0,
        metered: true,
    };
    work.charge(targets + 1)?;
    for (i, &point) in points.iter().enumerate() {
        if i % 256 == 0 {
            checkpoint(cancelled)?;
        }
        admit_point(point, settings.image_plane)?;
    }
    let mut prepared = PreparedSources {
        sources: Vec::new(),
        bound,
    };
    prepared
        .sources
        .try_reserve_exact(bound)
        .map_err(|_| allocation())?;
    let mut stages = Stages {
        config,
        settings,
        tree_spec,
        sources: prepared,
        resolved: Vec::new(),
        tree: SourceTree::default(),
        work,
        report: WakeAdvectionReport::default(),
        assemble,
        variation,
        cancelled,
    };
    let mut scratch = Scratch::new(points.len())?;
    let mut start_limit = stages
        .prepare(points, 0.0, true)
        .map_err(TrialError::refusal)?;
    for (p, k) in points.iter().zip(&mut scratch.first) {
        *k = stages.velocity(*p, 0.0)?;
    }

    let mut report = WakeAdaptiveReport::default();
    let mut elapsed = 0.0;
    let mut next_step = max_step;
    for _ in 0..time.max_trials {
        if elapsed == dt {
            break;
        }
        checkpoint(stages.cancelled)?;
        stages.work.charge(1)?;
        let remaining = dt - elapsed;
        let proposed_end = elapsed + next_step;
        // Integrate the actual final remainder instead of leaving a tiny
        // roundoff tail with no representable intermediate stage times.
        // This can exceed the nominal maximum step only by time roundoff.
        let end = if next_step >= remaining || dt - proposed_end <= 8.0 * f64::EPSILON * dt {
            dt
        } else {
            proposed_end
        };
        let h = end - elapsed;
        let middle = elapsed + 0.5 * h;
        let third = elapsed + 0.75 * h;
        if !(elapsed < middle && middle < third && third < end) {
            return Err(step_exhausted());
        }
        let trial = stages.trial(
            points,
            &mut scratch,
            elapsed,
            end,
            start_limit,
            line_displacement_m * (h / dt),
            time,
        );
        let retry_factor = match trial {
            Ok(estimate) if estimate.ratio <= 1.0 => {
                stages.work.charge(targets)?;
                points.copy_from_slice(&scratch.points);
                // k4 belongs to the accepted endpoint and is the next k1.
                // Rejected trials never overwrite the retained starting k1.
                std::mem::swap(&mut scratch.first, &mut scratch.second);
                start_limit = estimate.endpoint_limit;
                report.min_accepted_step_s = if report.accepted_substeps == 0 {
                    h
                } else {
                    report.min_accepted_step_s.min(h)
                };
                report.max_accepted_step_s = report.max_accepted_step_s.max(h);
                report.max_accepted_position_error_m =
                    report.max_accepted_position_error_m.max(estimate.error_m);
                report.max_accepted_error_ratio =
                    report.max_accepted_error_ratio.max(estimate.ratio);
                report.accepted_substeps += 1;
                elapsed = end;
                next_step = (h * step_factor(estimate.ratio)).clamp(time.min_step_s, max_step);
                continue;
            }
            Ok(estimate) => step_factor(estimate.ratio).min(0.5),
            Err(TrialError::Resolution) => 0.5,
            Err(TrialError::Fatal(error)) => return Err(error),
        };
        report.rejected_trials += 1;
        next_step = (h * retry_factor).max(time.min_step_s);
        if next_step >= h {
            return Err(step_exhausted());
        }
    }
    if elapsed != dt {
        return Err(refusal(
            "wake-advection-trials-exhausted",
            "adaptive wake evolution exhausted its trial cap before the interval endpoint",
        ));
    }
    // The final derivative was evaluated on the accepted endpoint sheet, so
    // its source geometry, evolved cores and boundary were already admitted.
    checkpoint(stages.cancelled)?;
    report.advection = stages.report;
    report.advection.work_used = stages.work.used;
    Ok(report)
}

struct Scratch {
    points: Vec<[f64; 3]>,
    first: Vec<[f64; 3]>,
    second: Vec<[f64; 3]>,
    third: Vec<[f64; 3]>,
}

impl Scratch {
    fn new(count: usize) -> Result<Self, Refusal> {
        fn bank(count: usize) -> Result<Vec<[f64; 3]>, Refusal> {
            let mut bank = Vec::new();
            bank.try_reserve_exact(count).map_err(|_| allocation())?;
            bank.resize(count, [0.0; 3]);
            Ok(bank)
        }
        Ok(Self {
            points: bank(count)?,
            first: bank(count)?,
            second: bank(count)?,
            third: bank(count)?,
        })
    }
}

struct Estimate {
    error_m: f64,
    ratio: f64,
    endpoint_limit: f64,
}

// Only this integrator's own resolution checks can request a retry. In
// particular, a caller's refusal retains its meaning even if its code happens
// to match one of the wake's public refusal codes.
enum TrialError {
    Fatal(Refusal),
    Resolution,
}

impl From<Refusal> for TrialError {
    fn from(error: Refusal) -> Self {
        Self::Fatal(error)
    }
}

impl TrialError {
    fn refusal(self) -> Refusal {
        match self {
            Self::Fatal(error) => error,
            Self::Resolution => refusal(
                "wake-advection-displacement",
                "a resolved edge collapsed during an adaptive trial",
            ),
        }
    }
}

struct Stages<'a, A, V, C> {
    config: SelfInducedShed,
    settings: StageSettings,
    tree_spec: Option<WakeTreeSpec>,
    sources: PreparedSources,
    resolved: Vec<bool>,
    tree: SourceTree,
    work: WorkMeter,
    report: WakeAdvectionReport,
    assemble: A,
    variation: &'a mut V,
    cancelled: &'a mut C,
}

impl<A, V, C> Stages<'_, A, V, C>
where
    A: FnMut(&[[f64; 3]], f64, &mut PreparedSources) -> Result<(), Refusal>,
    V: FnMut([f64; 3], f64) -> Result<[f64; 3], Refusal>,
    C: FnMut() -> bool,
{
    fn prepare(
        &mut self,
        points: &[[f64; 3]],
        fraction: f64,
        initial: bool,
    ) -> Result<f64, TrialError> {
        checkpoint(self.cancelled)?;
        self.work.charge(self.sources.bound as u64 + 1)?;
        self.sources.sources.clear();
        (self.assemble)(points, fraction, &mut self.sources)?;
        self.sources.images(self.settings.image_plane)?;
        if initial {
            self.resolved
                .try_reserve_exact(self.sources.sources.len())
                .map_err(|_| allocation())?;
            self.resolved.extend(
                self.sources
                    .sources
                    .iter()
                    .map(|s| s.geometry.start != s.geometry.end),
            );
        } else if self.resolved.len() == self.sources.sources.len()
            && self
                .sources
                .sources
                .iter()
                .zip(&self.resolved)
                .any(|(s, active)| *active && s.geometry.length == 0.0)
        {
            // A trial can numerically collapse a previously resolved edge.
            // Its initial geometry and positive core were already admitted;
            // retry the time resolution, not invalid source data or topology.
            return Err(TrialError::Resolution);
        }
        let limit =
            self.sources.scale(Some(&self.resolved))? * self.config.max_displacement_fraction;
        if let Some(spec) = self
            .tree_spec
            .filter(|s| s.absolute_velocity_error_mps > 0.0)
        {
            self.tree
                .rebuild(&self.sources.sources, spec, &mut self.work, self.cancelled)?;
            self.report.peak_tree_nodes = self.report.peak_tree_nodes.max(self.tree.node_count());
        }
        Ok(limit)
    }

    fn velocity(&mut self, point: [f64; 3], elapsed: f64) -> Result<[f64; 3], Refusal> {
        velocity_at(
            point,
            elapsed,
            &self.sources,
            self.config.ambient_velocity_mps,
            self.tree_spec,
            &mut self.tree,
            &mut self.work,
            &mut self.report,
            self.variation,
            self.cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn trial(
        &mut self,
        points: &[[f64; 3]],
        scratch: &mut Scratch,
        start: f64,
        end: f64,
        start_limit: f64,
        line_step: f64,
        time: WakeTimeSpec,
    ) -> Result<Estimate, TrialError> {
        let h = end - start;
        let dt = self.settings.dt_s;
        let plane = self.settings.image_plane;
        resolution(line_step, start_limit)?;
        self.work.charge(points.len() as u64)?;
        for (i, &point) in points.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint(self.cancelled)?;
            }
            let delta = resolved_displacement(scratch.first[i], h, start_limit)?;
            scratch.points[i] = resolved_translation(point, delta.map(|v| 0.5 * v), plane)?;
        }
        let middle_time = start + 0.5 * h;
        let mut limit = start_limit.min(self.prepare(&scratch.points, middle_time / dt, false)?);
        resolution(line_step, limit)?;
        for (p, k) in scratch.points.iter().zip(&mut scratch.second) {
            *k = self.velocity(*p, middle_time)?;
        }
        self.work.charge(points.len() as u64)?;
        for (i, &point) in points.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint(self.cancelled)?;
            }
            let delta = resolved_displacement(scratch.second[i], h, limit)?;
            scratch.points[i] = resolved_translation(point, delta.map(|v| 0.75 * v), plane)?;
        }
        let third_time = start + 0.75 * h;
        limit = limit.min(self.prepare(&scratch.points, third_time / dt, false)?);
        resolution(line_step, limit)?;
        for (p, k) in scratch.points.iter().zip(&mut scratch.third) {
            *k = self.velocity(*p, third_time)?;
        }
        self.work.charge(points.len() as u64)?;
        for (i, &point) in points.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint(self.cancelled)?;
            }
            resolved_displacement(scratch.third[i], h, limit)?;
            let delta = third_order_displacement(scratch, i, h);
            resolution(norm(delta), limit)?;
            scratch.points[i] = resolved_translation(point, delta, plane)?;
        }
        let endpoint_limit = self.prepare(&scratch.points, end / dt, false)?;
        limit = limit.min(endpoint_limit);
        resolution(line_step, limit)?;
        self.work.charge(points.len() as u64)?;
        let mut estimate = Estimate {
            error_m: 0.0,
            ratio: 0.0,
            endpoint_limit,
        };
        for i in 0..points.len() {
            let last = self.velocity(scratch.points[i], end)?;
            let d4 = resolved_displacement(last, h, limit)?;
            let high = third_order_displacement(scratch, i, h);
            resolution(norm(high), limit)?;
            // b3-b2 = (-5/72, 1/12, 1/9, -1/8). Difference form cancels
            // a constant field exactly and avoids subtracting world positions.
            // Scale before subtracting so large opposing velocities do not
            // overflow an otherwise representable local displacement estimate.
            let error = std::array::from_fn(|c| {
                (-5.0 / 72.0 * (h * scratch.first[i][c]) - (-5.0 / 72.0 * d4[c]))
                    + (1.0 / 12.0 * (h * scratch.second[i][c]) - 1.0 / 12.0 * d4[c])
                    + (1.0 / 9.0 * (h * scratch.third[i][c]) - 1.0 / 9.0 * d4[c])
            });
            let low = std::array::from_fn(|c| high[c] - error[c]);
            let error_m = norm(error);
            let tolerance = time.absolute_position_tolerance_m
                + time.relative_position_tolerance * norm(high).max(norm(low));
            if !error_m.is_finite() || !tolerance.is_finite() {
                return Err(refusal(
                    "wake-advection-nonfinite",
                    "adaptive position estimate exceeds the numerical domain",
                )
                .into());
            }
            estimate.error_m = estimate.error_m.max(error_m);
            estimate.ratio = estimate.ratio.max(error_m / tolerance);
            // Keep k1 until acceptance; a rejected trial can retry without
            // resampling its unchanged initial field or rebuilding its tree.
            scratch.second[i] = last;
        }
        Ok(estimate)
    }
}

fn third_order_displacement(scratch: &Scratch, i: usize, h: f64) -> [f64; 3] {
    std::array::from_fn(|c| {
        2.0 / 9.0 * (h * scratch.first[i][c])
            + 1.0 / 3.0 * (h * scratch.second[i][c])
            + 4.0 / 9.0 * (h * scratch.third[i][c])
    })
}

fn norm(v: [f64; 3]) -> f64 {
    v[0].hypot(v[1]).hypot(v[2])
}

fn resolved_displacement(velocity: [f64; 3], h: f64, limit: f64) -> Result<[f64; 3], TrialError> {
    displacement(velocity, h, limit).map_err(|error| {
        if error.code == "wake-advection-displacement" {
            TrialError::Resolution
        } else {
            TrialError::Fatal(error)
        }
    })
}

fn resolved_translation(
    point: [f64; 3],
    delta: [f64; 3],
    plane: Option<ImagePlane>,
) -> Result<[f64; 3], TrialError> {
    translated(point, delta, plane).map_err(|error| {
        if error.code == "image-plane-invalid" {
            TrialError::Resolution
        } else {
            TrialError::Fatal(error)
        }
    })
}

fn resolution(distance: f64, limit: f64) -> Result<(), TrialError> {
    if !distance.is_finite() {
        return Err(refusal(
            "wake-advection-nonfinite",
            "adaptive displacement exceeds the numerical domain",
        )
        .into());
    }
    admit_line_step(distance, limit).map_err(|_| TrialError::Resolution)
}

fn step_factor(ratio: f64) -> f64 {
    if ratio == 0.0 {
        5.0
    } else {
        (0.9 / ratio.cbrt()).clamp(0.1, 5.0)
    }
}

fn step_exhausted() -> Refusal {
    refusal(
        "wake-advection-time-step-exhausted",
        "adaptive accuracy or geometry resolution requires a smaller representable timestep",
    )
}

fn work_exhausted() -> Refusal {
    refusal(
        "wake-advection-work-exhausted",
        "adaptive source/target admission exceeds its work cap",
    )
}
