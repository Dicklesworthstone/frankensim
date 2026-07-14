//! fs-rep-mesh — mesh charts (plan §7.2). Layer: L2.
//!
//! - [`HalfEdgeMesh`]: manifold connectivity with boundary loops, a
//!   flip-edit core, and an invariant checker the property battery runs
//!   after every random edit;
//! - [`TetComplex`]: oriented volume elements whose signed incidence
//!   operators satisfy δδ = 0 EXACTLY (integer arithmetic) — the
//!   pre-FEEC sanity that makes fs-feec's exact sequences possible;
//! - [`Soup`] + generalized winding numbers ([`winding_exact`],
//!   [`WindingOctree`]): robust inside/outside on broken input, with the
//!   dipole octree's error MEASURED against exact;
//! - [`MeshChart`]: BVH-accelerated point-triangle magnitude with a
//!   generalized-winding sign and Woop watertight raycasts; raw soup exposes
//!   `NoClaim`, no Lipschitz bound, and Estimate/NoClaim numerical evidence;
//! - [`repair`]: dedupe/degenerate/orientation/hole pipeline with
//!   structured receipts;
//! - [`dual_contour`] + [`bracket_certificate`]: bounded uniform-grid
//!   extraction and an `ExactDistance`-only, enclosure-backed whole-triangle
//!   proximity certificate with structured refusal/cancellation progress;
//! - [`shapes`]: the public mesh fixture vocabulary (cube, icosphere,
//!   deterministic corruption).

mod chart;
mod complex;
mod contour;
mod convert;
mod halfedge;
mod repair;
pub mod shapes;
mod winding;

pub use chart::{Bvh, MeshChart, point_triangle_distance, ray_triangle_watertight};
pub use complex::{HexComplex, Incidence, TetComplex};
pub use contour::{
    BracketCertificateError, BracketEvidenceIssue, BracketFailure, BracketGeometryStage,
    BracketReport, ContourArithmeticStage, ContourError, ContourSampleStage, DC_MAX_CELLS_PER_AXIS,
    DcOptions, DcStats, NoLipschitz, Placement, bracket_certificate, dual_contour,
    dual_contour_clipped,
};
pub use convert::{IncrementalMeshSdf, MeshQuality, MeshSdfError, assess_quality, mesh_to_sdf};
pub use halfedge::{HalfEdge, HalfEdgeMesh, MeshBuildError, NO_FACE};
pub use repair::{RepairOutcome, RepairReceipt, repair};
pub use winding::{Soup, WindingOctree, triangle_winding, winding_exact};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    #[test]
    fn version_is_stamped() {
        assert!(!super::VERSION.is_empty());
    }
}
