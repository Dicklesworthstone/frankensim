//! fs-rep-sdf — signed-distance-field charts (plan §7.2). Layer: L2.
//!
//! Four representations, one honesty discipline:
//! - [`TiledSdf`] — dense Morton-tiled grids (fs-substrate fields), f32
//!   STORAGE / f64 EVALUATION, C¹ triquadratic B-spline reconstruction
//!   (shape optimization differentiates through samples), a CONSTRUCTED
//!   error enclosure, MEASURED eikonal statistics, and sphere-traced
//!   raycasts that respect the chart's own bounds;
//! - [`VdbGrid`] — FrankenVDB, the in-house sparse hierarchical tile tree
//!   (root map → 32³ internal → 8³ bitmasked leaves): the shared spatial
//!   substrate LBM lattices, CutFEM background grids, and voxel charts
//!   sit on, with deterministic iteration and measured footprint stats;
//! - [`NarrowBand`] — band-limited level sets on the VDB for level-set
//!   evolution (advect + reinitialize + rebuild, drift measured);
//! - [`AdaptiveSdf`] — octree with per-cell trilinear fits, refined where
//!   probe residuals exceed tolerance (Estimate-grade error, ledgered).
//!
//! Everything implements (or serves) fs-geom's [`fs_geom::Chart`], so the
//! agreement checker and conversion receipts apply unchanged.

mod adaptive;
mod band;
mod dense;
mod vdb;

pub use adaptive::{ADAPTIVE_MAX_NODES, AdaptiveSdf, AdaptiveStats};
pub use band::{
    BandStats, NARROW_BAND_MAX_SAMPLES_PER_AXIS, NARROW_BAND_MAX_SCAN_SAMPLES, NarrowBand,
};
pub use dense::{DENSE_MAX_SAMPLES_PER_AXIS, EikonalStats, SdfBuildError, TiledSdf};
pub use vdb::{VdbGrid, VdbStats};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    #[test]
    fn version_is_stamped() {
        assert!(!super::VERSION.is_empty());
    }
}
