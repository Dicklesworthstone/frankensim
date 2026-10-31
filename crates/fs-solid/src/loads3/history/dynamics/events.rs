//! Event-aligned transient histories on the original Newmark/load-jump owner.
//! The schedule is borrowed immutably; finite knots are indexed once. Exact
//! left limits alter held VALUES, never the physical frame-evaluation time.
//! Periodic vector tables require an unrolled finite schedule or the explicit
//! load-jump API; this lane refuses them instead of silently missing wraps.
use super::*;
use crate::loads3::transient::events::LoadJump3;
use fs_scenario::frame::{Frame, FrameMotion};
use fs_scenario::payload::{OutsideDomainPolicy, SampleSource, TableInterpolation};
use fs_scenario::signal::Interp;
use fs_qty::QtyAny;

fn invalid(what: &'static str) -> TransientHistoryError3 {
    TransientError3::Invalid(what).into()
}
fn checkpoint(cx: &Cx<'_>) -> Result<(), TransientHistoryError3> {
    cx.checkpoint().map_err(|_| TransientError3::Cancelled.into())
}

/// Finite authored event index, including vector/scalar knots and frame tilts.
/// The borrowed history and its frame tree cannot change during an active run.
/// Indexing does not certify resolution of a smooth signal between its knots.
pub struct EventLoadHistory3<'a> {
    history: TetLoadHistory3<'a>,
    knots: Vec<f64>,
}
impl<'a> EventLoadHistory3<'a> {
    /// Bound setup by an aggregate scalar-entry allowance (hard ceiling 2^20).
    /// Scalar tables are validated before binary search; validated vector
    /// payloads keep their own constructor invariants. Duplicate knot times are
    /// coalesced exactly, not within an epsilon. At most 256 frame rows and 1024
    /// surface histories are admitted. Setup polls between bounded scans.
    pub fn new(history: TetLoadHistory3<'a>, maximum_schedule_entries: usize, cx: &Cx<'_>)
        -> Result<Self, TransientHistoryError3> {
        checkpoint(cx)?;
        if maximum_schedule_entries == 0 || maximum_schedule_entries > 1_048_576
            || history.surfaces.len() > 1024 || history.frames.frames.len() > 256 {
            return Err(invalid("event schedule exceeds its explicit setup envelope"));
        }
        let mut index = Index { knots: Vec::new(), remaining: maximum_schedule_entries, cx };
        for row in history.surfaces {
            checkpoint(cx)?;
            match row {
                SurfaceHistory3::Fixed(_) => {}
                SurfaceHistory3::Traction { value, .. } => index.vector(value)?,
                SurfaceHistory3::Pressure { value, .. } => index.signal(value)?,
            }
        }
        match history.volume {
            VolumeHistory3::Fixed(_) | VolumeHistory3::FixedAcceleration(_) => {}
            VolumeHistory3::ForceDensity(v) | VolumeHistory3::Acceleration(v) => index.vector(v)?,
        }
        for frame in &history.frames.frames {
            checkpoint(cx)?;
            if frame.name.len() > 4096 { return Err(invalid("event frame names exceed 4096 bytes")); }
            if let FrameMotion::Tilt { angle, .. } = &frame.motion { index.signal(angle)?; }
        }
        checkpoint(cx)?;
        index.knots.sort_by(f64::total_cmp);
        index.knots.dedup_by(|a, b| *a == *b);
        checkpoint(cx)?;
        Ok(Self { history, knots: index.knots })
    }
    pub fn knots_s(&self) -> &[f64] { &self.knots }

    /// Initialize with the actual RIGHT-continuous load. Operators and frozen
    /// thermal data remain owned by the original sparse transient implementation.
    #[allow(clippy::too_many_arguments)]
    pub fn begin_transient<'h, 'g, 'm, 'p, 's>(&'h self, geometry: &'g TetLoadGeometry3<'m>,
        problem: &'p TetLinearElasticProblem<'s>, thermal: Option<TetThermalStrainField<'_>>,
        time_s: f64, displacement_m: &[[f64; 3]], velocity_m_s: &[[f64; 3]],
        config: TransientConfig3, cx: &Cx<'_>)
        -> Result<EventTransientRun3<'h, 'a, 'g, 'm, 'p, 's>, TransientHistoryError3> {
        let run = self.history.begin_transient(geometry, problem, thermal, time_s,
            displacement_m, velocity_m_s, config, cx)?;
        Ok(EventTransientRun3 { schedule: self, run })
    }

    fn endpoint(&self, start: f64, target: f64, maximum_step: f64) -> Result<f64, TransientHistoryError3> {
        let span = target - start;
        if !target.is_finite() || !span.is_finite() || span <= 0.0
            || !maximum_step.is_finite() || maximum_step <= 0.0 {
            return Err(invalid("event stepping needs a finite later target and positive finite step cap"));
        }
        let mut end = if span <= maximum_step { target } else { start + maximum_step };
        if end - start > maximum_step { end = end.next_down(); }
        if let Some(&knot) = self.knots.get(self.knots.partition_point(|&t| t <= start)) {
            end = end.min(knot);
        }
        if end <= start || !end.is_finite() || start + (end - start) != end {
            return Err(invalid("event endpoint cannot be represented by a positive exact clock advance"));
        }
        Ok(end)
    }

    // Sample small constant replacements only where a held vector changes.
    // Never clone an entire payload table or approximate t-minus by next_down.
    fn left_at<'g, 'm, 'p, 's>(&self, geometry: &'g TetLoadGeometry3<'m>,
        problem: &'p TetLinearElasticProblem<'s>, time: f64, cx: &Cx<'_>)
        -> Result<TetLoadCase3<'g, 'm>, TransientHistoryError3> {
        let mut frames = FrameTree::new();
        for frame in &self.history.frames.frames {
            checkpoint(cx)?;
            let motion = match &frame.motion {
                FrameMotion::Tilt { axis, center, angle } => FrameMotion::Tilt {
                    axis: *axis, center: *center, angle: TimeSignal::Constant(before(angle, time)?),
                },
                other => other.clone(),
            };
            frames.frames.push(Frame { id: frame.id, parent: frame.parent, name: frame.name.clone(), motion });
        }
        let mut vectors = Vec::with_capacity(self.history.surfaces.len());
        let mut pressures = Vec::with_capacity(self.history.surfaces.len());
        for row in self.history.surfaces {
            checkpoint(cx)?;
            vectors.push(match row { SurfaceHistory3::Traction { value, .. } => vector_before(value, time)?, _ => None });
            pressures.push(match row { SurfaceHistory3::Pressure { value, .. } => Some(TimeSignal::Constant(before(value, time)?)), _ => None });
        }
        let surfaces: Vec<_> = self.history.surfaces.iter().enumerate().map(|(i, row)| match *row {
            SurfaceHistory3::Fixed(value) => SurfaceHistory3::Fixed(value),
            SurfaceHistory3::Traction { faces, value } => SurfaceHistory3::Traction { faces, value: vectors[i].as_ref().unwrap_or(value) },
            SurfaceHistory3::Pressure { faces, value } => SurfaceHistory3::Pressure { faces, value: pressures[i].as_ref().unwrap_or(value) },
        }).collect();
        let volume_value = match self.history.volume {
            VolumeHistory3::Fixed(_) | VolumeHistory3::FixedAcceleration(_) => None,
            VolumeHistory3::ForceDensity(v) | VolumeHistory3::Acceleration(v) => vector_before(v, time)?,
        };
        let volume = match self.history.volume {
            VolumeHistory3::Fixed(v) => VolumeHistory3::Fixed(v),
            VolumeHistory3::FixedAcceleration(v) => VolumeHistory3::FixedAcceleration(v),
            VolumeHistory3::ForceDensity(v) => VolumeHistory3::ForceDensity(volume_value.as_ref().unwrap_or(v)),
            VolumeHistory3::Acceleration(v) => VolumeHistory3::Acceleration(volume_value.as_ref().unwrap_or(v)),
        };
        let left = TetLoadHistory3 { surfaces: &surfaces, volume, frames: &frames,
            nodal_forces_n: self.history.nodal_forces_n };
        checkpoint(cx)?;
        Ok(left.at(geometry, problem, time, cx)?.case)
    }
}

struct Index<'a, 'cx> { knots: Vec<f64>, remaining: usize, cx: &'a Cx<'cx> }
impl Index<'_, '_> {
    fn charge(&mut self, count: usize) -> Result<(), TransientHistoryError3> {
        checkpoint(self.cx)?;
        self.remaining = self.remaining.checked_sub(count).ok_or_else(|| invalid("event schedule scalar-entry budget exhausted"))?;
        Ok(())
    }
    fn signal(&mut self, signal: &TimeSignal) -> Result<(), TransientHistoryError3> {
        let (count, probe) = match signal {
            TimeSignal::Constant(_) => (1, 0.0),
            TimeSignal::Ramp { t_start, .. } => (4, *t_start),
            TimeSignal::Table { times, values, .. } => (times.len().checked_add(values.len())
                .ok_or_else(|| invalid("event table size overflow"))?, times.first().copied().unwrap_or(0.0)),
            TimeSignal::Chebfun(p) => (p.cheb.coeffs().len().checked_add(2).ok_or_else(|| invalid("event profile size overflow"))?, p.cheb.domain().0),
        };
        self.charge(count)?;
        signal.eval(probe).map_err(|e| LoadHistoryError3::Scenario(e.to_string()))?;
        match signal {
            TimeSignal::Constant(_) => {}
            TimeSignal::Ramp { t_start, t_end, .. } => self.knots.extend([*t_start, *t_end]),
            TimeSignal::Table { times, .. } => {
                for (i, &time) in times.iter().enumerate() {
                    if i % 256 == 0 { checkpoint(self.cx)?; }
                    self.knots.push(time);
                }
            }
            TimeSignal::Chebfun(p) => { let (a, b) = p.cheb.domain(); self.knots.extend([a, b]); }
        }
        Ok(())
    }
    fn vector(&mut self, vector: &VectorPayload) -> Result<(), TransientHistoryError3> {
        if vector.components() != 3 { return Err(invalid("event loads need three-component vectors")); }
        match vector.source() {
            SampleSource::Fixed(values) => self.charge(values.len())?,
            SampleSource::Distribution { .. } => return Err(invalid("event loads require explicit histories, not distributions")),
            SampleSource::Table { times, outside, .. } => {
                if *outside == OutsideDomainPolicy::Periodic {
                    return Err(invalid("unroll periodic vector histories into a finite table or use explicit load jumps"));
                }
                self.charge(times.len().checked_mul(4).ok_or_else(|| invalid("event vector size overflow"))?)?;
                for (i, time) in times.iter().enumerate() {
                    if i % 256 == 0 { checkpoint(self.cx)?; }
                    self.knots.push(time.value);
                }
            }
        }
        Ok(())
    }
}
fn before(signal: &TimeSignal, time: f64) -> Result<QtyAny, TransientHistoryError3> {
    let value = signal.eval(time).map_err(|e| LoadHistoryError3::Scenario(e.to_string()))?;
    if let TimeSignal::Table { times, values, dims, interp: Interp::Hold } = signal {
        let i = times.partition_point(|&t| t < time);
        if i > 0 && times.get(i) == Some(&time) { return Ok(QtyAny::new(values[i-1], *dims)); }
    }
    Ok(value)
}
fn vector_before(vector: &VectorPayload, time: f64) -> Result<Option<VectorPayload>, TransientHistoryError3> {
    if let SampleSource::Table { times, values, interpolation: TableInterpolation::StepLeft, .. } = vector.source() {
        let i = times.partition_point(|t| t.value < time);
        if i > 0 && times.get(i).is_some_and(|t| t.value == time) {
            return VectorPayload::new(vector.meta().clone(), SampleSource::Fixed(values[i-1].clone()))
                .map(Some).map_err(|e| LoadHistoryError3::Scenario(e.to_string()).into());
        }
    }
    Ok(None)
}

/// Original accepted interval, plus an optional finite endpoint force jump.
#[derive(Debug, Clone, PartialEq)]
pub struct ScheduledStep3 {
    pub start_time_s: f64,
    pub end_time_s: f64,
    pub interval: TransientStep3,
    pub jump: Option<LoadJump3>,
}

/// A trajectory bound to one immutable event schedule. Clone forks accepted
/// motion while sharing the original operators. No unbounded report history.
#[derive(Clone)]
pub struct EventTransientRun3<'h, 'a, 'g, 'm, 'p, 's> {
    schedule: &'h EventLoadHistory3<'a>,
    run: TetTransientRun3<'g, 'm, 'p, 's>,
}
impl<'g, 'm, 'p, 's> EventTransientRun3<'_, '_, 'g, 'm, 'p, 's> {
    /// Read-only access to motion, response recovery and the existing diagnostics.
    pub fn run(&self) -> &TetTransientRun3<'g, 'm, 'p, 's> { &self.run }

    /// Advance ONE interval to the next authored knot, target or step cap.
    /// Call repeatedly to reach target. Left/right loads and frame tilts use the
    /// same exact endpoint. Finite jumps change acceleration without impulse;
    /// failure in either side/solve leaves the complete accepted run unchanged.
    /// Smooth forcing still requires a timestep study; no time-error certificate.
    pub fn advance_to(&mut self, target_s: f64, maximum_step_s: f64, cx: &Cx<'_>)
        -> Result<ScheduledStep3, TransientHistoryError3> {
        checkpoint(cx)?;
        let start_time_s = self.run.state().time_s;
        let end_time_s = self.schedule.endpoint(start_time_s, target_s, maximum_step_s)?;
        let geometry = self.run.geometry; let problem = self.run.problem;
        let left = self.schedule.left_at(geometry, problem, end_time_s, cx)?;
        let right = self.schedule.history.at(geometry, problem, end_time_s, cx)?;
        let mut candidate = self.run.clone();
        let interval = candidate.advance(&left, end_time_s-start_time_s, cx)?;
        let jump = if left.nodal_forces_n() != right.case.nodal_forces_n() {
            Some(candidate.apply_load_jump(&right.case, cx)?)
        } else { None };
        checkpoint(cx)?;
        self.run = candidate;
        Ok(ScheduledStep3 { start_time_s, end_time_s, interval, jump })
    }
}

/// Event-preserving local motion-error control over the original integrator.
pub mod adaptive;
