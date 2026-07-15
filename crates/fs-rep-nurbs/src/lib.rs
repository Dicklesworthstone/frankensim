//! fs-rep-nurbs (plan §7.2): rational B-spline charts — bounded EXACT spline
//! algebra (knot insertion and degree elevation in i128 rational arithmetic,
//! so "identical before/after refinement" is provable equality rather than a
//! tolerance while every intermediate remains representable), trimmed patches with CERTIFIED point
//! classification (convex-hull winding with exact subdivision), measured
//! closest-point bracket estimates (f64 branch-and-bound pending an
//! outward-rounded upgrade), and the
//! HONEST Boolean position: route through SDF by default, refuse direct
//! B-rep Booleans without a certificate.
//!
//! Layer: L2 (MORPH). Runtime deps: `std`, fs-evidence, fs-exec, fs-geom,
//! fs-ivl, fs-math.

pub mod basis;
pub mod boolean;
pub mod closest;
pub mod curve;
pub mod rat;
#[cfg(feature = "nurbs-refit")]
pub mod refit;
#[cfg(feature = "nurbs-sdf")]
pub mod sdf;
pub mod surface;
pub mod trim;

pub use basis::{AdmittedKnotVector, BasisRun, KnotAdmissionRun, KnotSpanRun, KnotVector, Scalar};
pub use boolean::{BooleanOp, BooleanPolicy, BooleanRefusal, boolean};
pub use closest::{
    ClosestPointRun, DistanceBracketEstimate, closest_point_curve, closest_point_surface,
};
pub use curve::{
    AdmittedNurbsCurve, CurveAdmissionRun, CurveBezierRun, CurveCloneRun, CurveDerivativesRun,
    CurveEvaluationRun, CurveInsertionRun, CurveRemovalRun, CurveSpanBoxesRun, NurbsCurve,
};
pub use rat::Rat;
pub use surface::{
    AdmittedNurbsSurface, NurbsSurface, SurfaceAdmissionRun, SurfaceCloneRun, SurfaceEvaluationRun,
    SurfaceInsertionRun, SurfacePartialsRun, SurfaceSpanBoxesRun,
};
pub use trim::{
    AdmittedTrimLoop, AdmittedTrimmedPatch, Classification, TrimClassificationRun, TrimLoop,
    TrimLoopAdmissionRun, TrimLoopReversalRun, TrimmedPatch, TrimmedPatchAdmissionRun,
};

use core::fmt;

/// Crate version (compile-time stamp).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Structured spline failures (Decalogue P10).
#[derive(Debug, Clone, PartialEq)]
pub enum NurbsError {
    /// A structurally invalid spline (knot ordering, count mismatches,
    /// non-positive weights).
    Structure {
        /// Diagnosis.
        what: String,
    },
    /// A parameter-domain violation or defensive work/memory envelope refusal.
    Domain {
        /// Diagnosis.
        what: String,
    },
    /// The exact-arithmetic domain was exceeded. Fallible exact helpers return
    /// this refusal; current `Rat` operator traits cannot transport it and fail
    /// with a named panic at that boundary rather than wrapping.
    Exactness {
        /// Diagnosis.
        what: String,
    },
}

impl fmt::Display for NurbsError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NurbsError::Structure { what } => write!(f, "invalid spline: {what}"),
            NurbsError::Domain { what } => write!(f, "parameter out of domain: {what}"),
            NurbsError::Exactness { what } => write!(f, "exact-arithmetic overflow: {what}"),
        }
    }
}

impl std::error::Error for NurbsError {}
