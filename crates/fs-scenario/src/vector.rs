//! Execute the existing typed vector payloads in Cartesian coordinates.
//!
//! This is a time sample, not a stochastic realization, event processor or
//! trajectory derivative. Metadata remains authoritative: only the explicit
//! `basis/world-cartesian` (world frame) and `basis/frame-cartesian` conventions
//! are executable here. Reference resets and distributions require their own
//! execution owners and refuse rather than becoming ordinary fixed values.
//! Table endpoint policies affect the VALUE time, not the frame's physical
//! evaluation time. Proper rigid rotations act on both orientation parities;
//! translation never acts on a free vector. No interval/phase-error claim.

use crate::frame::{FrameTree, WORLD};
use crate::payload::{
    OrientationParity, OutsideDomainPolicy, QuantityContract, ReferenceSemantics,
    SampleSource, TableInterpolation, VectorPayload,
};
use crate::ScenarioError;
use fs_ga::{Quat, Vec3};
use fs_qty::QtyAny;
use std::ops::ControlFlow;

/// Quantity and parity survive rotation; only component coordinates change.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WorldVector3 {
    pub values: [f64; 3],
    pub contract: QuantityContract,
    pub orientation: OrientationParity,
}

/// No value is returned after a policy, arithmetic, frame or caller refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum VectorEvaluationError {
    Unsupported(&'static str),
    InvalidTime,
    OutsideDomain { time_s: f64, first_s: f64, last_s: f64 },
    Arithmetic(&'static str),
    Frame(ScenarioError),
    Cancelled,
}
impl std::fmt::Display for VectorEvaluationError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Cartesian vector evaluation refused: {self:?}")
    }
}
impl std::error::Error for VectorEvaluationError {}

fn poll(checkpoint: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), VectorEvaluationError> {
    if checkpoint().is_break() { Err(VectorEvaluationError::Cancelled) } else { Ok(()) }
}
fn components(sample: &[QtyAny]) -> [f64; 3] {
    [sample[0].value, sample[1].value, sample[2].value]
}

impl VectorPayload {
    /// Evaluate a validated three-component payload in WORLD coordinates.
    ///
    /// Fixed, StepLeft and Linear sources use the existing immutable payload
    /// constructors' shape/unit guarantees. Search is logarithmic in table size;
    /// no table is cloned or revalidated per query. Exact knots retain their
    /// stored sample, including the final knot. Periodic reduction is applied
    /// only OUTSIDE the closed table interval and wraps to the first knot at a
    /// complete period. Unrepresentable time differences/periods refuse.
    ///
    /// Frame evaluation reuses `FrameTree::world_pose` (including moving
    /// ancestors), then extracts its rotor rather than subtracting translated
    /// points. Frames are borrowed, never cached across caller mutations.
    /// Checkpoints surround sampling/frame evaluation and precede publication;
    /// existing frame/signal evaluations are indivisible, not newly certified
    /// cancellation-latency bounds. Frame trees above 256 rows refuse this lane.
    pub fn world_vector_at(
        &self, frames: &FrameTree, time_s: f64,
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<WorldVector3, VectorEvaluationError> {
        use VectorEvaluationError as E;
        poll(checkpoint)?;
        if !time_s.is_finite() { return Err(E::InvalidTime); }
        if self.components() != 3 { return Err(E::Unsupported("exactly three Cartesian components required")); }
        if !matches!(self.meta().reference(), ReferenceSemantics::Continuous) {
            return Err(E::Unsupported("reference-reset payloads require an event-aware owner"));
        }
        match self.meta().basis_key().as_str() {
            "basis/world-cartesian" if self.meta().frame() == WORLD => {}
            "basis/frame-cartesian" => {}
            _ => return Err(E::Unsupported("declare world-cartesian in WORLD or frame-cartesian in its named frame")),
        }
        if frames.frames.len() > 256 { return Err(E::Unsupported("Cartesian evaluation admits at most 256 frame rows")); }
        let local = match self.source() {
            SampleSource::Fixed(value) => components(value),
            SampleSource::Distribution { .. } => return Err(E::Unsupported("a distribution is not an explicit time history")),
            SampleSource::Table { times, values, interpolation, outside } => {
                table_value(times, values, *interpolation, *outside, time_s)?
            }
        };
        poll(checkpoint)?;
        let pose = frames.world_pose(self.meta().frame(), time_s).map_err(E::Frame)?;
        poll(checkpoint)?;
        let q = Quat::from_rotor(&pose);
        let rotated = q.rotate(Vec3::new(local[0], local[1], local[2]));
        let values = [rotated.x, rotated.y, rotated.z];
        if !values.iter().all(|v| v.is_finite()) {
            return Err(E::Arithmetic("nonfinite rotated vector"));
        }
        poll(checkpoint)?;
        Ok(WorldVector3 { values, contract: self.meta().contract(), orientation: self.meta().orientation() })
    }
}

fn table_value(
    times: &[QtyAny], values: &[Vec<QtyAny>], interpolation: TableInterpolation,
    outside: OutsideDomainPolicy, time_s: f64,
) -> Result<[f64; 3], VectorEvaluationError> {
    use VectorEvaluationError as E;
    // Nonempty, increasing, finite time and equal shapes are constructor invariants.
    let first = times[0].value;
    let last = times[times.len() - 1].value;
    let mut time = time_s;
    if time < first || time > last {
        time = match outside {
            OutsideDomainPolicy::Refuse => return Err(E::OutsideDomain { time_s, first_s: first, last_s: last }),
            OutsideDomainPolicy::Clamp => time.clamp(first, last),
            OutsideDomainPolicy::Periodic => {
                let period = last - first;
                let delta = time - first;
                if !period.is_finite() || period <= 0.0 || !delta.is_finite() {
                    return Err(E::Arithmetic("unrepresentable periodic time reduction"));
                }
                let phase = delta.rem_euclid(period);
                let wrapped = first + phase;
                if !wrapped.is_finite() { return Err(E::Arithmetic("nonfinite wrapped time")); }
                wrapped.clamp(first, last)
            }
        };
    }
    let left = times.partition_point(|t| t.value <= time).saturating_sub(1);
    if left + 1 == times.len() || time == times[left].value
        || interpolation == TableInterpolation::StepLeft {
        return Ok(components(&values[left]));
    }
    let width = times[left + 1].value - times[left].value;
    if !width.is_finite() || width <= 0.0 { return Err(E::Arithmetic("unrepresentable interpolation span")); }
    let weight = (time - times[left].value) / width;
    if !weight.is_finite() || !(0.0..=1.0).contains(&weight) {
        return Err(E::Arithmetic("unrepresentable interpolation weight"));
    }
    // Convex form avoids overflow in (right-left) for opposite large samples.
    let out = std::array::from_fn(|i|
        (1.0 - weight) * values[left][i].value + weight * values[left + 1][i].value);
    if out.iter().all(|v| v.is_finite()) { Ok(out) }
    else { Err(E::Arithmetic("nonfinite interpolated component")) }
}
