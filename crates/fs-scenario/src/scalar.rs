//! Execute typed scalar histories without erasing their physical metadata.

use crate::frame::{FrameTree, WORLD};
use crate::payload::{OrientationParity, QuantityContract, ReferenceSemantics, ScalarPayload};
use std::ops::ControlFlow;

pub use crate::sample::PayloadEvaluationError as ScalarEvaluationError;

/// The original quantity contract accompanies the sampled SI value.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScalarSample {
    pub value: f64,
    pub contract: QuantityContract,
    /// Reflection parity remains explicit even though frames use proper rotations.
    pub orientation: OrientationParity,
}

impl ScalarPayload {
    /// Sample an ordinary scalar at explicit physical time, without allocation.
    ///
    /// Fixed, held and linear sources share the vector evaluator's logarithmic
    /// lookup, exact-knot behavior and outside-domain policy. Scalar values do
    /// not rotate or translate, but their declared frame is resolved at physical
    /// time, including when the value table clamps or wraps. This introduces
    /// no motional field, voltage-reference shift or spatial-profile inference.
    ///
    /// Admits `basis/scalar`, or the existing Cartesian basis conventions.
    /// Odd parity, reset references, distributions, unknown bases and frame
    /// trees above 256 rows refuse; the original metadata is never discarded.
    pub fn scalar_at(
        &self,
        frames: &FrameTree,
        time_s: f64,
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<ScalarSample, ScalarEvaluationError> {
        if checkpoint().is_break() {
            return Err(ScalarEvaluationError::Cancelled);
        }
        if self.meta().orientation() != OrientationParity::Even {
            return Err(ScalarEvaluationError::Unsupported(
                "ordinary scalar evaluation requires even parity",
            ));
        }
        self.oriented_scalar_at(frames, time_s, checkpoint)
    }

    /// Sample a scalar or pseudoscalar without erasing its reflection parity.
    ///
    /// Proper rigid frame transformations leave both values unchanged. The
    /// returned parity distinguishes ordinary scalars from oriented quantities
    /// such as axial magnetic flux density dotted with a polar outward normal.
    /// This does not admit reflected frames or infer a boundary-normal sign.
    /// Source, basis, frame, reference and cancellation policies are shared
    /// with [`Self::scalar_at`]; consumers must check the required parity.
    pub fn oriented_scalar_at(
        &self,
        frames: &FrameTree,
        time_s: f64,
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<ScalarSample, ScalarEvaluationError> {
        use ScalarEvaluationError as E;
        let poll = |check: &mut dyn FnMut() -> ControlFlow<()>| {
            if check().is_break() {
                Err(E::Cancelled)
            } else {
                Ok(())
            }
        };
        poll(checkpoint)?;
        if !time_s.is_finite() {
            return Err(E::InvalidTime);
        }
        if !matches!(self.meta().reference(), ReferenceSemantics::Continuous) {
            return Err(E::Unsupported(
                "reference-reset payloads require an event-aware owner",
            ));
        }
        match self.meta().basis_key().as_str() {
            "basis/scalar" | "basis/frame-cartesian" => {}
            "basis/world-cartesian" if self.meta().frame() == WORLD => {}
            _ => {
                return Err(E::Unsupported(
                    "declare a scalar basis or an explicit Cartesian frame convention",
                ));
            }
        }
        if frames.frames.len() > 256 {
            return Err(E::Unsupported(
                "scalar evaluation admits at most 256 frame rows",
            ));
        }
        let [value] = crate::sample::value(self.source(), time_s, |q| [q.value])?;
        poll(checkpoint)?;
        frames
            .world_pose(self.meta().frame(), time_s)
            .map_err(E::Frame)?;
        poll(checkpoint)?;
        Ok(ScalarSample {
            value,
            contract: self.meta().contract(),
            orientation: self.meta().orientation(),
        })
    }
}
