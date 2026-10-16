//! fs-rep-mesh — mesh charts (plan §7.2). Layer: L2.
//!
//! - [`HalfEdgeMesh`]: manifold connectivity with boundary loops, a
//!   flip-edit core, and an invariant checker the property battery runs
//!   after every random edit;
//! - [`TetComplex`]: oriented volume elements whose signed incidence
//!   operators satisfy δδ = 0 EXACTLY (integer arithmetic) — the
//!   pre-FEEC sanity that makes fs-feec's exact sequences possible;
//! - [`TriComplex2`]: a genuine oriented 2-D complex with exact d1*d0,
//!   typed stable feature IDs, selected-side trace maps, and explicit planar
//!   or axisymmetric measure metadata;
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
mod polygon;
mod repair;
pub mod shapes;
/// Bounded diffuse enclosure geometry with opaque visibility and explicit
/// sampling/projection diagnostics; estimates, not topology certificates.
pub mod view_factors;
mod winding;

pub use chart::{Bvh, MeshChart, point_triangle_distance, ray_triangle_watertight};
pub use complex::{
    HexComplex, Incidence, Metric2, Metric2Error, TetComplex, TraceEdge2, TraceMap2, TriComplex2,
    TriComplex2Error, TriComplex2LineageId, TriComplex2LineageSchema, TriFeatureId,
    TriFeatureSchema, tri_complex2_lineage_id,
};
pub use contour::{
    BracketCertificateError, BracketEvidenceIssue, BracketFailure, BracketGeometryStage,
    BracketReport, ContourArithmeticStage, ContourError, ContourSampleStage, DC_MAX_CELLS_PER_AXIS,
    DcOptions, DcStats, NoLipschitz, Placement, bracket_certificate, dual_contour,
    dual_contour_clipped,
};
pub use convert::{IncrementalMeshSdf, MeshQuality, MeshSdfError, assess_quality, mesh_to_sdf};
pub use halfedge::{HalfEdge, HalfEdgeMesh, MeshBuildError, NO_FACE};
pub use polygon::{MAX_POLYGON_VERTICES, MAX_POLYGON_WORK, PolygonError, triangulate_polygon};
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
