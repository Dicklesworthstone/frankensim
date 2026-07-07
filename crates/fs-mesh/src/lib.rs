//! fs-mesh — body-fitted tet meshing (plan §7.5). Layer: L2.
//!
//! When a body-fitted mesh is WANTED (final verification, shells,
//! export): BRIO-ordered incremental Delaunay tetrahedralization on
//! EXACT predicates (fs-ivl `orient3d`/`insphere`, SoS tie-breaking) —
//! remembering CutFEM-on-SDF exists precisely so that meshing stays
//! optional inside optimization loops.
//!
//! What this crate certifies, it AUDITS with the same exact predicates
//! it builds with ([`Tetrahedralization::audit`]): the local Delaunay
//! property on every internal facet (the Delaunay lemma makes local ⇒
//! global), positive orientation of every tet, mutual adjacency, the
//! Euler-characteristic ball check, and exact convexity of the boundary
//! hull. Degenerate inputs (grids: massively cospherical/coplanar
//! configurations) complete correctly BECAUSE the predicates are exact —
//! cospherical ties resolve deterministically (`Zero` = not in
//! conflict), so identical input bytes give identical meshes (P2).
//!
//! v1 kernel scope: sequential Bowyer–Watson with ghost tets (hull at
//! infinity), jump-and-walk location with BRIO locality hints, cavity
//! GROWTH repair for degenerate visibility, radius-edge quality
//! refinement by circumcenter insertion. Constrained boundary recovery
//! (PLC conformity), full Ruppert with local feature size, sliver
//! exudation, and parallel domain coloring are the successor bead —
//! recorded as CONTRACT no-claims, not silently absent.

mod delaunay;
mod refine;

pub use delaunay::{
    AuditReport, DelaunayStats, MeshError, Tetrahedralization, delaunay, GHOST,
};
pub use refine::{RefineOptions, RefineStats, refine};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    #[test]
    fn version_is_stamped() {
        assert!(!super::VERSION.is_empty());
    }
}
