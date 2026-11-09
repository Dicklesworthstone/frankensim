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
    OrientationParity, Payload, QuantityContract, ReferenceSemantics, VectorPayload,
};
use crate::{BcKind, BcValue, BoundaryCondition, Physics};
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

pub use crate::sample::PayloadEvaluationError as VectorEvaluationError;

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
        let local = crate::sample::value(self.source(), time_s, |v| components(v))?;
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
    /// Sample physical heat-flux density in WORLD Cartesian components [W/m²].
    ///
    /// The native boundary owner must contract these components with each
    /// original outward unit normal. Translation does not act on this free
    /// vector; moving axes rotate it at physical time. The row prescribes no
    /// total power, normal direction, mesh motion or stochastic realization.
    pub fn heat_flux_world_at(
        &self,
        frames: &FrameTree,
        time_s: f64,
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<[f64; 3], VectorEvaluationError> {
        use VectorEvaluationError as E;
        poll(checkpoint)?;
        if self.physics != Physics::Thermal || self.kind != BcKind::Neumann {
            return Err(E::Unsupported("heat-flux vector requires Thermal/Neumann"));
        }
        if self.region.is_empty() || self.compatibility.is_some() {
            return Err(E::Unsupported("heat flux requires a named patch and no mass-flow compatibility tag"));
        }
        let Some(BcValue::Typed(Payload::Vector(value))) = &self.value else {
            return Err(E::Unsupported("physical heat flux requires a typed vector"));
        };
        if value.meta().frame().0 != self.frame {
            return Err(E::Unsupported("heat-flux payload and boundary must declare the same component frame"));
        }
        if value.meta().contract() != QuantityContract::Dimensions(Dims([0, 1, -3, 0, 0, 0])) {
            return Err(E::Unsupported("heat flux requires a dimensions-only W/m² vector"));
        }
        if value.meta().orientation() != OrientationParity::Even {
            return Err(E::Unsupported("heat flux requires ordinary even-parity vector components"));
        }
        value.world_vector_at(frames, time_s, checkpoint).map(|sample| sample.values)
    }

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
