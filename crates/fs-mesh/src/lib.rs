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
//! refinement by circumcenter insertion, sliver exudation, and
//! deterministic parallel domain coloring (read-parallel rounds,
//! canonical application — bitwise thread-count-invariant).
//! Constrained boundary recovery (PLC conformity) and full-Ruppert
//! quality remain successor scope — recorded as CONTRACT no-claims,
//! not silently absent.

mod delaunay;
mod exude;
mod parallel;
mod recovery;
mod refine;
mod remesh;

pub use delaunay::{AuditReport, DelaunayStats, GHOST, MeshError, Tetrahedralization, delaunay};
pub use exude::{ExudeOptions, ExudeStats, exude};
pub use parallel::{ColoredStats, delaunay_colored, delaunay_colored_reversed};
pub use recovery::{Correspondence, RecoveryOptions, RecoveryStats, recover_segments};
pub use refine::{RefineOptions, RefineStats, refine};
pub use remesh::{MetricField, RemeshOptions, RemeshStats, UniformMetric, remesh};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    #[test]
    fn version_is_stamped() {
        assert!(!super::VERSION.is_empty());
    }
}
