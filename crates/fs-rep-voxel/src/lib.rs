//! fs-rep-voxel (plan §7.2): occupancy/multi-material voxel charts on the
//! sparse VDB substrate (shared with fs-rep-sdf), exact Euclidean
//! distance transforms, point clouds with estimated normals (the FITTING
//! TARGET role in scan-to-Region workflows), and explicit lattice/strut
//! graphs (FrankenNetworkx) with watertight solid realization.
//!
//! Layer: L2 (MORPH). Runtime deps: `std`, fs-rep-sdf (VdbGrid), fs-geom
//! (Chart), fs-exec (Cx), fs-evidence, fs-math, fnx-classes/fnx-runtime
//! (constellation).

pub mod chart;
pub mod cloud;
pub mod dt;
pub mod field;
pub mod lattice;

pub use chart::OccupancyChart;
pub use cloud::PointCloud;
pub use dt::{DistanceField, euclidean_dt};
pub use field::{DensityField, MaterialField, OccupancyField};
pub use lattice::{LatticeGraph, LatticeNode, Strut};

use core::fmt;

/// Crate version (compile-time stamp).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Structured voxel-representation failures (Decalogue P10).
#[derive(Debug, Clone, PartialEq)]
pub enum VoxelError {
    /// A field parameter is inadmissible (voxel size, fraction range…).
    Parameters {
        /// Diagnosis.
        what: String,
    },
    /// A lattice graph is structurally degenerate.
    Lattice {
        /// Diagnosis.
        what: String,
    },
    /// A point-cloud query cannot be answered as posed.
    Cloud {
        /// Diagnosis.
        what: String,
    },
    /// FrankenNetworkx round-trip failure.
    Graph {
        /// Diagnosis.
        what: String,
    },
}

impl fmt::Display for VoxelError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            VoxelError::Parameters { what } => write!(f, "bad voxel parameters: {what}"),
            VoxelError::Lattice { what } => write!(f, "degenerate lattice: {what}"),
            VoxelError::Cloud { what } => write!(f, "point-cloud query failed: {what}"),
            VoxelError::Graph { what } => write!(f, "graph round-trip failed: {what}"),
        }
    }
}

impl std::error::Error for VoxelError {}
