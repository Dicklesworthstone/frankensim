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
    OrientationParity, OutsideDomainPolicy, Payload, QuantityContract, ReferenceSemantics,
    SampleSource, TableInterpolation, VectorPayload,
};
use crate::{BcKind, BcValue, BoundaryCondition, Physics, ScenarioError};
use fs_ga::{Quat, Vec3};
use fs_qty::{Dims, QtyAny};
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

impl BoundaryCondition {
    /// Evaluate prescribed fluid velocity in WORLD Cartesian components [m/s].
    ///
    /// This consumes the typed-vector form of IncompressibleFlow/Dirichlet.
    /// Components describe absolute velocity relative to the inertial WORLD,
    /// expressed in the payload's declared axes. A moving frame rotates those
    /// components at `time_s`; it does not add a frame-origin or wall velocity.
    /// Wall-relative velocity requires an explicit kinematic composition by the
    /// owning solver. No direction is guessed from a legacy scalar Dirichlet.
    ///
    /// The original vector owner supplies fixed/StepLeft/Linear sampling,
    /// outside-domain policies, physical-time frame rotation and cancellation.
    /// The payload is borrowed without copying or rescanning its history. This
    /// returns a spatially uniform velocity on the named patch, not an integrated
    /// mass flow or a boundary-geometry resolution. The solver still owns its
    /// dimensional scaling, boundary enforcement and stability admission.
    ///
    /// # Errors
    /// Refuses a different boundary role, missing patch, compatibility tag,
    /// directionless value, inconsistent frame, non-velocity quantity or odd
    /// parity. The vector owner also refuses unsupported bases, component
    /// counts, source/reference policies, invalid time/frames and cancellation.
    pub fn velocity_world_at(
        &self,
        frames: &FrameTree,
        time_s: f64,
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<[f64; 3], VectorEvaluationError> {
        use VectorEvaluationError as E;
        poll(checkpoint)?;
        if self.physics != Physics::IncompressibleFlow || self.kind != BcKind::Dirichlet {
            return Err(E::Unsupported(
                "prescribed velocity requires IncompressibleFlow/Dirichlet",
            ));
        }
        if self.region.is_empty() || self.compatibility.is_some() {
            return Err(E::Unsupported("prescribed velocity requires a named patch and no total-mass-flow compatibility tag"));
        }
        let Some(BcValue::Typed(Payload::Vector(value))) = &self.value else {
            return Err(E::Unsupported("prescribed velocity requires a typed vector; scalar values have no implied direction"));
        };
        if value.meta().frame().0 != self.frame {
            return Err(E::Unsupported(
                "velocity payload and boundary must declare the same component frame",
            ));
        }
        // The sealed QuantityKind catalog has angular velocity but no
        // translational-velocity kind. Do not reinterpret a future semantic
        // contract solely because it happens to share m/s dimensions.
        if value.meta().contract() != QuantityContract::Dimensions(Dims([1, 0, -1, 0, 0, 0])) {
            return Err(E::Unsupported(
                "prescribed velocity requires a dimensions-only m/s vector",
            ));
        }
        if value.meta().orientation() != OrientationParity::Even {
            return Err(E::Unsupported(
                "prescribed velocity requires ordinary even-parity vector components",
            ));
        }
        value
            .world_vector_at(frames, time_s, checkpoint)
            .map(|value| value.values)
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
