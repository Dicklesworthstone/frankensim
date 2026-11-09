//! Execute typed scalar histories without erasing their physical metadata.

use crate::frame::{FrameTree, WORLD};
use crate::payload::{
    OrientationParity, Payload, QuantityContract, ReferenceSemantics, ScalarPayload,
};
use crate::{BcKind, BcValue, BoundaryCondition, Physics};
use fs_qty::Pressure;
use fs_qty::semantic::{QuantityKind, ValueForm};
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

/// Fluid pressure is not solid stress or an acoustic amplitude, even though
/// all three have Pa dimensions. Gauge/absolute reference remains explicit
/// at the receiving physical owner.
pub(crate) fn pressure_contract(contract: QuantityContract) -> bool {
    match contract {
        QuantityContract::Dimensions(dims) => dims == Pressure::DIMS,
        QuantityContract::Semantic(semantic) => {
            semantic.kind() == QuantityKind::Pressure
                && matches!(
                    semantic.form(),
                    ValueForm::Static | ValueForm::Instantaneous
                )
        }
        QuantityContract::Heterogeneous => false,
    }
}

/// Thermal quantity contracts share the boundary vocabulary's dimensions.
/// Absolute temperature cannot be replaced by a same-dimension interval.
pub(crate) fn thermal_contract(kind: BcKind, contract: QuantityContract) -> bool {
    use crate::bc::{Expectation, expectation};
    let expected = match expectation(Physics::Thermal, kind) {
        Expectation::Value(dims) | Expectation::ValueOrScalar(dims)
        | Expectation::ValueOrScalarOrVector(dims) => dims,
        _ => return false,
    };
    match contract {
        QuantityContract::Dimensions(dims) => dims == expected,
        QuantityContract::Semantic(semantic) => {
            kind == BcKind::Dirichlet
                && semantic.kind() == QuantityKind::AbsoluteTemperature
                && semantic.form() == ValueForm::Static
        }
        QuantityContract::Heterogeneous => false,
    }
}

impl BoundaryCondition {
    /// Sample typed temperature, normal heat-flux density or heat-transfer
    /// coefficient at explicit physical time. A temperature difference cannot
    /// prescribe an absolute temperature. The receiver owns the normal-flux
    /// sign and supplies a separate Robin reference temperature; neither is
    /// inferred from the scalar or its frame. This does not move the patch.
    ///
    /// The ordinary scalar owner retains table/tail policies, frame validation,
    /// even parity and cancellation without copying the declared history.
    pub fn thermal_at(
        &self,
        frames: &FrameTree,
        time_s: f64,
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<ScalarSample, ScalarEvaluationError> {
        use ScalarEvaluationError as E;
        if checkpoint().is_break() {
            return Err(E::Cancelled);
        }
        if self.physics != Physics::Thermal
            || !matches!(self.kind, BcKind::Dirichlet | BcKind::Neumann | BcKind::Robin)
        {
            return Err(E::Unsupported(
                "thermal scalar requires Thermal/Dirichlet, Neumann or Robin",
            ));
        }
        if self.region.is_empty() || self.compatibility.is_some() {
            return Err(E::Unsupported(
                "thermal scalar requires a named patch and no mass-flow compatibility tag",
            ));
        }
        let Some(BcValue::Typed(Payload::Scalar(value))) = &self.value else {
            return Err(E::Unsupported("typed thermal boundary requires a scalar payload"));
        };
        if value.meta().frame().0 != self.frame {
            return Err(E::Unsupported(
                "thermal payload and boundary must declare the same frame",
            ));
        }
        if !thermal_contract(self.kind, value.meta().contract()) {
            return Err(E::Unsupported(
                "thermal scalar requires absolute K, normal W/m^2 or W/(m^2 K), matching its boundary kind",
            ));
        }
        value.scalar_at(frames, time_s, checkpoint)
    }

    /// Sample a typed fluid pressure outlet, retaining its original contract.
    ///
    /// The scalar owner supplies fixed/held/linear histories, declared table
    /// tails, physical-time frame validation and cancellation. Only ordinary
    /// Pa values or static/instantaneous fluid-pressure semantics are admitted.
    /// Frames do not add hydrostatic or dynamic pressure. The receiving solver
    /// must explicitly interpret absolute/gauge pressure and bind the patch.
    /// No history is copied or rescanned during sampling.
    pub fn pressure_at(
        &self,
        frames: &FrameTree,
        time_s: f64,
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<ScalarSample, ScalarEvaluationError> {
        use ScalarEvaluationError as E;
        if checkpoint().is_break() {
            return Err(E::Cancelled);
        }
        if self.physics != Physics::IncompressibleFlow || self.kind != BcKind::PressureOutlet {
            return Err(E::Unsupported(
                "fluid pressure requires IncompressibleFlow/PressureOutlet",
            ));
        }
        if self.region.is_empty() || self.compatibility.is_some() {
            return Err(E::Unsupported(
                "fluid pressure requires a named patch and no mass-flow compatibility tag",
            ));
        }
        let Some(BcValue::Typed(Payload::Scalar(value))) = &self.value else {
            return Err(E::Unsupported(
                "typed fluid pressure requires a scalar payload",
            ));
        };
        if value.meta().frame().0 != self.frame {
            return Err(E::Unsupported(
                "pressure payload and boundary must declare the same frame",
            ));
        }
        if !pressure_contract(value.meta().contract()) {
            return Err(E::Unsupported(
                "fluid pressure requires Pa or static/instantaneous Pressure semantics",
            ));
        }
        value.scalar_at(frames, time_s, checkpoint)
    }
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
