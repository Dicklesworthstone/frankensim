//! Prescribed support translation during a bowed/plucked material performance.
//! This samples a declared piecewise-linear geometry path, not a force envelope
//! or a new mechanical integrator. Each resulting endpoint joins the SAME step
//! as due bow and force controls. The contact owner computes actuator work.

use super::*;
use fs_exec::CancelGate;
use crate::render::schedule::force::coupled::ModalCouplingError;
use crate::render::schedule::force::coupled::contact::friction::{ModalFrictionError, multiple::MultiFrictionFrame};
pub use crate::render::schedule::force::coupled::contact::friction::multiple::obstacles::motion::ObstacleMotionLimits;

/// One geometry endpoint on the original physical clock.
#[derive(Clone, Debug, PartialEq)]
pub struct StopMotionKnot {
    /// Elapsed samples since this performance starts. 1 is the endpoint of the
    /// first rendered step, NOT its start. The implicit sample-zero knot is zero
    /// translation at every authored stop; a different initial gap belongs in
    /// its source geometry and contact-aware preload.
    pub sample: u64,
    /// One signed normal translation [m] per stop in source order. Positive
    /// opens the gap; negative presses in, regardless of the normal's direction.
    /// These are absolute offsets, not increments, pitches or modal weights.
    pub offsets_m: Vec<f64>,
}

/// A bounded continuous path in the fixed normal coordinates of existing stops.
#[derive(Clone, Debug, PartialEq)]
pub struct StopMotionProgram {
    /// Strictly increasing positive sample coordinates, within the performance.
    /// Linear interpolation connects knots; the final position is then held.
    /// A discontinuous position jump cannot be declared at one instant.
    pub knots: Vec<StopMotionKnot>,
    /// Original contact owner's translation and step-average speed limits.
    pub limits: ObstacleMotionLimits,
    /// Maximum authored knots, checked before allocating runtime buffers.
    pub max_knots: usize,
    /// Maximum scalar coordinates across all knots, not a per-knot allowance.
    pub max_coordinate_values: usize,
}

fn poll(gate: &CancelGate) -> Result<(), RenderError> {
    if gate.is_requested() { Err(RenderError::Coupled(ModalCouplingError::Cancelled)) }
    else { Ok(()) }
}
impl StopMotionProgram {
    /// Validate the complete path against an explicit physical clock and stop count.
    /// File/scene loaders can use this same admission before material I/O. The
    /// binding repeats it against the actual bank; every step also checks its
    /// rounded endpoint against the original physical limits.
    pub fn validate(&self, stops: usize, samples: u64, rate: u32, gate: &CancelGate) -> Result<(), RenderError> {
        poll(gate)?;
        self.limits.validate().map_err(RenderError::Coupled)?;
        if stops == 0 || samples == 0 || rate == 0 || samples > (1_u64 << 53) || self.knots.is_empty()
            || self.knots.len() > self.max_knots
            || stops.checked_mul(self.knots.len()).is_none_or(|n| n > self.max_coordinate_values) {
            return Err(render_invalid("stop motion needs existing stops, an exactly representable clock and bounded nonempty knots"));
        }
        let mut previous: Option<&StopMotionKnot> = None;
        for knot in &self.knots {
            poll(gate)?;
            let start = previous.map_or(0, |p| p.sample);
            if knot.sample <= start || knot.sample > samples || knot.offsets_m.len() != stops {
                return Err(render_invalid("stop-motion knots must be ordered complete endpoints inside the performance"));
            }
            let duration = (knot.sample-start) as f64 / f64::from(rate);
            for (i, &end) in knot.offsets_m.iter().enumerate() {
                let begin = previous.map_or(0.0, |p| p.offsets_m[i]);
                let velocity = (end-begin)/duration;
                if !end.is_finite() || end.abs() > self.limits.maximum_abs_offset_m
                    || !velocity.is_finite() || velocity.abs() > self.limits.maximum_speed_m_s {
                    return Err(render_invalid("stop trajectory exceeds its physical translation or speed limit"));
                }
            }
            previous = Some(knot);
        }
        Ok(())
    }
}

pub(super) struct BoundStopMotion {
    source: StopMotionProgram,
    pub(super) candidate: Vec<f64>,
    pub(super) next: usize,
}
impl BoundStopMotion {
    // Endpoint arithmetic is bounded by admission. Use exact authored values
    // at knots and hold the final one; only the last accepted cursor is reused.
    // No per-sample allocation or control-history expansion.
    pub(super) fn prepare(&mut self, endpoint: u64) -> usize {
        prepare_path(&self.source.knots, self.next, endpoint, None, &mut self.candidate)
    }
}

// One endpoint sampler for normal translation and longitudinal station paths.
// Their units/source admission remain distinct; their clock and interpolation
// arithmetic are shared. None preserves the original implicit zero offsets.
pub(super) trait GeometryKnot {
    fn sample(&self) -> u64;
    fn coordinates(&self) -> &[f64];
}
impl GeometryKnot for StopMotionKnot {
    fn sample(&self) -> u64 { self.sample }
    fn coordinates(&self) -> &[f64] { &self.offsets_m }
}
pub(super) fn prepare_path<K: GeometryKnot>(knots: &[K], mut cursor: usize, endpoint: u64,
    initial: Option<&[f64]>, candidate: &mut [f64]) -> usize
{
        while knots.get(cursor).is_some_and(|k| k.sample() < endpoint) { cursor += 1; }
        if let Some(right) = knots.get(cursor) {
            if endpoint == right.sample() {
                candidate.copy_from_slice(right.coordinates());
            } else {
                let left = cursor.checked_sub(1).map(|i| &knots[i]);
                let start = left.map_or(0, |k| k.sample());
                let fraction = (endpoint-start) as f64/(right.sample()-start) as f64;
                for (i, value) in candidate.iter_mut().enumerate() {
                    let begin = left.map_or_else(|| initial.map_or(0.0, |v| v[i]), |k| k.coordinates()[i]);
                    *value = begin+(right.coordinates()[i]-begin)*fraction;
                }
            }
        } else {
            candidate.copy_from_slice(knots[knots.len()-1].coordinates());
        }
        cursor
}

impl MultiBowStringBank {
    /// Enable prescribed normal translation of this bank's original stops.
    /// Materials, physical directions, receivers, cavity and preload are kept.
    /// The generic contact owner retains limits and the actual accepted offsets.
    pub fn with_stop_motion(mut self, limits: ObstacleMotionLimits, gate: &CancelGate) -> Result<Self, RenderError> {
        self.system = self.system.with_obstacle_motion(limits, gate).map_err(RenderError::Friction)?;
        Ok(self)
    }
    /// Direct physical stepping with final stop offsets in `source().stops()` order.
    /// Bow load/speed and external forces act in the same complete transaction.
    pub fn step_with_stop_offsets(&mut self, external: &[f64], drives: &[FrictionDrive], offsets: &[f64])
        -> Result<&MultiFrictionFrame, ModalFrictionError>
    {
        if self.bow_stations.is_some() {
            self.step_with_bow_geometry(external, drives, None, None, Some(offsets), &CancelGate::new_clock_free())
        } else { self.system.step_with_obstacle_offsets(external, drives, offsets) }
    }
    /// Same physical transaction with cancellation polling.
    pub fn step_with_stop_offsets_under_gate(&mut self, external: &[f64], drives: &[FrictionDrive],
        offsets: &[f64], gate: &CancelGate) -> Result<&MultiFrictionFrame, ModalFrictionError>
    {
        if self.bow_stations.is_some() {
            self.step_with_bow_geometry(external, drives, None, None, Some(offsets), gate)
        } else { self.system.step_with_obstacle_offsets_under_gate(external, drives, offsets, gate) }
    }
}

impl MultiBowPerformance {
    /// Bind a complete support-motion path before any physical sample.
    ///
    /// Static contact preload and initial material-body vibration are preserved.
    /// Force preload may also be attached afterwards; translation is still zero
    /// at initial time. All motion, bow and force controls affect ONE mechanical
    /// step; its real pressure goes through the original PressureRenderer/PCM
    /// path. Changing callback sizes changes neither this path nor its clock.
    ///
    /// This supports pressing and withdrawing compliant contacts along their
    /// fixed normals. `with_stop_station_motion` can additionally prescribe
    /// longitudinal positions on the same clock. Tissue motion, moving frames,
    /// patch-area evolution and pitch assignment are not inferred.
    /// A later physical/pressure refusal poisons the performance as before,
    /// leaving its bank at the last jointly accepted sample and stop positions.
    pub fn with_stop_motion(mut self, program: StopMotionProgram, gate: &CancelGate) -> Result<Self, RenderError> {
        poll(gate)?;
        self.validate_sample_count(0)?;
        if self.bank.system.samples_rendered() != 0 || self.stop_motion.is_some()
            || self.bank.system.obstacle_motion_limits().is_some() {
            return Err(render_invalid("stop motion binds once before the physical performance begins"));
        }
        let stops = self.bank.system.obstacle_count();
        program.validate(stops, self.config.samples, self.bank.source.bank().options().sample_rate_hz, gate)?;
        let bound = BoundStopMotion { source: program, candidate: vec![0.0; stops], next: 0 };
        self.bank.system = self.bank.system.with_obstacle_motion(bound.source.limits, gate)
            .map_err(RenderError::Friction)?;
        self.stop_motion = Some(bound);
        Ok(self)
    }
    /// Original normal geometry path and physical limits, without mutable bypass.
    #[must_use]
    pub fn stop_motion_program(&self) -> Option<&StopMotionProgram> {
        self.stop_motion.as_ref().map(|m| &m.source)
    }
}
