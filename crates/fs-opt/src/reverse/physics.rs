//! Explicit physics bindings for scalar PDE study residuals.
//!
//! The IR owns composition; the bound model owns the solve, its adjoint, and
//! its numerical/work admission. Binding does not turn provider assertions
//! into certified physics. No foreign runtime, finite differences, or hidden
//! study-name lookup is installed by this interface.

pub mod box_map;

use super::{ReverseError, Slot, capacity};
use crate::{Expr, NodeId, Problem, VarId};
use fs_exec::Cx;
use fs_qty::Dims;

/// Signature that must agree with the exact bound IR node and its variable.
#[derive(Debug, Clone, Copy)]
pub struct PhysicsSignature<'a> {
    /// Exact study reference, not a fuzzy model name or an implicit default.
    pub study: &'a str,
    /// Ambient point coordinates of the single `over` variable.
    pub parameter_count: usize,
    /// Physical dimensions of each input coordinate.
    pub parameter_dims: Dims,
    /// Physical dimensions of the scalar study residual.
    pub residual_dims: Dims,
}

/// One coherent solved residual and its total derivative at the supplied point.
#[derive(Debug, Clone, PartialEq)]
pub struct PhysicsSample {
    /// Scalar residual of the named study (before any IR composition).
    pub value: f64,
    /// Total ambient derivative with respect to the declared `over` variable.
    /// Includes the implicit solved-state dependence, not just explicit terms.
    pub gradient: Vec<f64>,
}

/// Common typed failure protocol for bound physics models.
/// Models retain precise solver residuals as bits rather than flattening them
/// into an unstructured string. A failed call returns no sample.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PhysicsError {
    /// Malformed model inputs or unsupported model configuration.
    InvalidInput(&'static str),
    /// This point is outside the model's declared numerical/physical domain.
    Domain(&'static str),
    /// A numerical quantity became unrepresentable.
    NonFinite { quantity: &'static str, bits: u64 },
    /// The requested work was not funded. Not a convergence diagnosis.
    Budget { resource: &'static str, used: usize, limit: usize },
    /// A solved-state or adjoint residual failed its acceptance gate.
    NotConverged {
        phase: &'static str,
        iterations: usize,
        residual_bits: u64,
        tolerance_bits: u64,
    },
    /// Cancellation, translated to the existing OptError::Cancelled boundary.
    Cancelled,
}

impl core::fmt::Display for PhysicsError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "physics model refused: {self:?}")
    }
}
impl std::error::Error for PhysicsError {}

/// A deterministic scalar study and its first-order adjoint contract.
///
/// For its entire binding lifetime, the provider must keep its study meaning,
/// fixed model data and signature unchanged. Every dependency on design data
/// must enter through `point`; hidden dependencies on other IR variables are
/// not permitted. The evaluator verifies signatures, shapes and finiteness,
/// not adjoint correctness or a continuum error bound. Providers must enforce
/// their own primal/adjoint solve tolerances and work/memory caps and propagate
/// cancellation. Arbitrary provider panics or side effects are not contained.
/// Second-order derivatives require the separate, explicit [`PhysicsHessian`]
/// capability. Existing first-order providers need not implement it.
pub trait PhysicsModel: core::fmt::Debug {
    /// Static scientific meaning of this binding.
    fn signature(&self) -> PhysicsSignature<'_>;

    /// Solve at the supplied point and return a matching first-order sample.
    /// `None` explicitly disables cooperative cancellation for this call.
    fn value_gradient(
        &self,
        point: &[f64],
        cx: Option<&Cx<'_>>,
    ) -> Result<PhysicsSample, PhysicsError>;

    /// Optional total second-order derivative of this SAME scalar model.
    /// Returning `None` preserves the first-order-only refusal. This capability
    /// and the model meaning must remain fixed for the binding lifetime.
    fn second_order(&self) -> Option<&dyn PhysicsHessian> { None }
}

/// Explicit ambient Hessian action for a bound, solved scalar study.
///
/// The result is `d(gradient(point))[direction]` in the binding's declared
/// coordinates. It must include implicit state/adjoint dependence, not just
/// explicit parameter terms or a Gauss--Newton approximation. The IR supplies
/// the outer composition terms itself; providers must not add them twice.
///
/// Reusing immutable primal/adjoint linearizations is encouraged. A cache must
/// be keyed to the supplied point, NOT the last evaluated trial: rejected
/// trials and cloned Newton studies can revisit older accepted points. Solves,
/// attempted-work/memory ceilings and cancellation remain provider-owned.
/// No finite-difference fallback, global Hessian definiteness, derivative
/// certificate, provider panic containment or side-effect rollback is implied.
pub trait PhysicsHessian: core::fmt::Debug {
    /// Apply the scalar study Hessian at `point` to `direction`. Both slices
    /// have the admitted parameter count and contain only finite values.
    /// Return exactly that many finite components or the original typed error.
    fn hessian_vector_product(
        &self,
        point: &[f64],
        direction: &[f64],
        cx: Option<&Cx<'_>>,
    ) -> Result<Vec<f64>, PhysicsError>;
}

/// Explicit association with one PDE node. Multiple nodes may share a model,
/// but each distinct reachable node is evaluated once per primal sweep.
#[derive(Debug, Clone, Copy)]
pub struct PhysicsBinding<'a> {
    /// Exact node in the immutable problem.
    pub node: NodeId,
    /// Provider retained for this compiled program.
    pub model: &'a dyn PhysicsModel,
}

#[derive(Debug)]
pub(super) struct CompiledPhysics<'a> {
    pub node: NodeId,
    pub model: &'a dyn PhysicsModel,
    pub variable: VarId,
    pub gradient: Slot,
}

pub(super) fn admitted_bindings<'a>(
    problem: &Problem,
    bindings: &[PhysicsBinding<'a>],
) -> Result<Vec<PhysicsBinding<'a>>, ReverseError> {
    let mut ordered = capacity(bindings.len())?;
    ordered.extend_from_slice(bindings);
    ordered.sort_unstable_by_key(|binding| binding.node);
    for pair in ordered.windows(2) {
        if pair[0].node == pair[1].node {
            return Err(ReverseError::PhysicsBinding {
                node: pair[0].node, what: "duplicate physics binding",
            });
        }
    }
    for binding in &ordered {
        let refuse = |what| ReverseError::PhysicsBinding { node: binding.node, what };
        let Expr::PdeResidual { study, over, adjoint_available, dims } = problem.expr(binding.node)? else {
            return Err(refuse("physics binding must name a PDE residual node"));
        };
        if !adjoint_available {
            return Err(refuse("PDE node does not declare an available adjoint"));
        }
        let variable = problem.variable(*over)?;
        let signature = binding.model.signature();
        if signature.study != study.as_str() {
            return Err(refuse("physics study reference mismatch"));
        }
        if Some(signature.parameter_count as u64) != variable.manifold.point_dim().map(u64::from) {
            return Err(refuse("physics parameter count mismatch"));
        }
        if signature.parameter_dims != variable.dims || signature.residual_dims != *dims {
            return Err(refuse("physics input or residual dimensions mismatch"));
        }
    }
    Ok(ordered)
}
