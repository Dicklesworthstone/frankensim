//! fs-opdsl — the typed mathematical operator IR (plan patch Rev C):
//! physics operators represented symbolically over FEEC building
//! blocks and LOWERED so that residual apply, Jacobian-vector product,
//! vector-Jacobian product / DISCRETE ADJOINT, DWR indicators,
//! preconditioner structure hints, and MMS studies all come from ONE
//! SOURCE OF TRUTH. Layer: L3 FLUX.
//!
//! Why this crate exists (the patch's words): FrankenSim should not
//! merely contain kernels — it should contain a kernel-generating
//! algebra. Deriving primal/adjoint/estimator from a single typed
//! expression makes "the primal changed but the adjoint didn't" rot
//! structurally impossible, and lets the Gauntlet test generated
//! adjoints against generated JVPs mechanically.
//!
//! Design lineage: generalizes fs-vskeleton's `EdgeLaw` seed (one
//! definition → both primal and adjoint) to a registry of linear
//! atoms (exterior derivatives, Whitney masses, external matrices)
//! and opaque-but-differentiable pointwise laws, with degree AND
//! Qty-dimension checking at expression construction.
//!
//! Escape hatch policy: hand-written operators enter as `External`
//! atoms — allowed, but they must supply their transpose (or declare
//! symmetry) and they pass the SAME consistency gates; the plan
//! report marks each atom `derived` or `hand` so the tradeoff stays
//! visible.

pub mod atoms;
pub mod expr;
pub mod fixtures;
pub mod kernels;
pub mod law;
pub mod mms;
pub mod plan;
pub mod system;

pub use atoms::{Atom, AtomId, Transpose};
pub use expr::{Expr, OperatorDef, Space, TypeError};
pub use law::{CubicReaction, LawId, PointwiseLaw};
pub use mms::{MmsReport, mms_poisson_study};
pub use plan::{LoweredOperator, PlanReport, dwr_indicators};
pub use system::{
    AdmittedSystem, AtomSignature, BlockEquation, ClockTransferSignature, ConventionRef,
    CoordinateConvention, FieldDecl, FieldId, FieldQuantity, ParameterRole, SYSTEM_IR_VERSION,
    ScalarConvention, SpatialSupport, StateOwnership, SystemDef, SystemExpr, SystemId,
    SystemTypeError, TransformSignature,
};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
