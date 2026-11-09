//! Execute the existing second-rank Cartesian tensor histories.
//!
//! Both indices use the payload's declared orthonormal basis: a proper frame
//! rotation acts as `R T R^T`. Translation does not act on a tensor. This owner
//! samples components; physical roles such as stress or conductivity, tensor
//! symmetry and constitutive admissibility remain the consuming solver's job.

use crate::frame::{FrameTree, WORLD};
use crate::payload::{OrientationParity, QuantityContract, ReferenceSemantics, TensorPayload};
use fs_ga::{Quat, Vec3};
use fs_qty::Pressure;
use fs_qty::semantic::{QuantityKind, ValueForm};
use std::ops::ControlFlow;

pub use crate::sample::PayloadEvaluationError as TensorEvaluationError;

/// Row-major WORLD components retaining the original quantity and parity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WorldTensor3 {
    pub values: [[f64; 3]; 3],
    pub contract: QuantityContract,
    pub orientation: OrientationParity,
}

/// Quantity admission shared by prescribed traction rows and native stress
/// loaders. Pa dimensions or static/instantaneous Stress are supported;
/// Pressure semantics and acoustic amplitudes describe different quantities.
#[must_use]
pub fn stress_contract(contract: QuantityContract) -> bool {
    match contract {
        QuantityContract::Dimensions(dims) => dims == Pressure::DIMS,
        QuantityContract::Semantic(kind) => kind.kind() == QuantityKind::Stress
            && matches!(kind.form(), ValueForm::Static | ValueForm::Instantaneous),
        QuantityContract::Heterogeneous => false,
    }
}

fn poll(checkpoint: &mut impl FnMut() -> ControlFlow<()>) -> Result<(), TensorEvaluationError> {
    if checkpoint().is_break() {
        Err(TensorEvaluationError::Cancelled)
    } else {
        Ok(())
    }
}

impl TensorPayload {
    /// Sample a 3-by-3 Cartesian tensor at physical time, without allocation.
    ///
    /// Fixed, held and linear sources use the shared scalar/vector sampler,
    /// including exact knots and explicit outside-domain policies. Frame time
    /// remains `time_s` when the component table clamps or wraps. The original
    /// metadata, full tensor shape and both transformation indices survive.
    /// Exactly symmetric sampled inputs retain exact symmetry after rotation;
    /// nonsymmetric inputs are never silently symmetrized.
    ///
    /// Only `basis/world-cartesian` in WORLD or `basis/frame-cartesian` in its
    /// named frame are admitted. Reference resets, stochastic distributions,
    /// other shapes/bases, nonfinite times and frame trees over 256 rows refuse.
    /// Checkpoints surround sampling, frame resolution and publication, like
    /// the existing vector owner; frame evaluation remains indivisible.
    pub fn world_tensor_at(
        &self,
        frames: &FrameTree,
        time_s: f64,
        checkpoint: &mut impl FnMut() -> ControlFlow<()>,
    ) -> Result<WorldTensor3, TensorEvaluationError> {
        use TensorEvaluationError as E;
        poll(checkpoint)?;
        if !time_s.is_finite() {
            return Err(E::InvalidTime);
        }
        if self.rows() != 3 || self.columns() != 3 {
            return Err(E::Unsupported("exactly 3-by-3 Cartesian tensor components required"));
        }
        if !matches!(self.meta().reference(), ReferenceSemantics::Continuous) {
            return Err(E::Unsupported("reference-reset payloads require an event-aware owner"));
        }
        match self.meta().basis_key().as_str() {
            "basis/world-cartesian" if self.meta().frame() == WORLD => {}
            "basis/frame-cartesian" => {}
            _ => return Err(E::Unsupported("declare world-cartesian in WORLD or frame-cartesian in its named frame")),
        }
        if frames.frames.len() > 256 {
            return Err(E::Unsupported("Cartesian evaluation admits at most 256 frame rows"));
        }
        let local: [f64; 9] = crate::sample::value(self.source(), time_s, |sample| {
            std::array::from_fn(|i| sample[i].value)
        })?;
        poll(checkpoint)?;
        let pose = frames.world_pose(self.meta().frame(), time_s).map_err(E::Frame)?;
        poll(checkpoint)?;
        let rotor = Quat::from_rotor(&pose);
        let axes = [Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0), Vec3::new(0.0, 0.0, 1.0)]
            .map(|axis| rotor.rotate(axis));
        let rotation = [axes.map(|v| v.x), axes.map(|v| v.y), axes.map(|v| v.z)];
        let left: [[f64; 3]; 3] = std::array::from_fn(|i| std::array::from_fn(|j| {
            (0..3).map(|k| rotation[i][k] * local[3 * k + j]).sum()
        }));
        let symmetric = local[1] == local[3] && local[2] == local[6] && local[5] == local[7];
        let mut values = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in (if symmetric { i } else { 0 })..3 {
                let component: f64 = (0..3).map(|k| left[i][k] * rotation[j][k]).sum();
                if !component.is_finite() {
                    return Err(E::Arithmetic("nonfinite rotated tensor"));
                }
                values[i][j] = component;
                if symmetric {
                    values[j][i] = component;
                }
            }
        }
        poll(checkpoint)?;
        Ok(WorldTensor3 {
            values,
            contract: self.meta().contract(),
            orientation: self.meta().orientation(),
        })
    }
}
