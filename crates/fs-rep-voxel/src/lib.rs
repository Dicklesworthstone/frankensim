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
    /// Two occupancy fields use different world-space lattice frames.
    FrameMismatch {
        /// Boolean operation that refused the mixed frames.
        operation: &'static str,
        /// Left-hand voxel edge length.
        left_voxel_size: f64,
        /// Right-hand voxel edge length.
        right_voxel_size: f64,
        /// Left-hand world origin.
        left_origin: [f64; 3],
        /// Right-hand world origin.
        right_origin: [f64; 3],
    },
    /// A requested coordinate box cannot be represented safely.
    CoordinateRange {
        /// Operation that required the box.
        operation: &'static str,
        /// Axis whose padded range was not representable.
        axis: usize,
        /// Unpadded lower coordinate.
        min: i32,
        /// Unpadded upper coordinate.
        max: i32,
        /// Halo requested on both sides.
        halo: u32,
    },
    /// A dense operation exceeded its explicit voxel budget.
    VoxelBudgetExceeded {
        /// Operation that requested dense work.
        operation: &'static str,
        /// Dense voxels required by the checked bounding box.
        required: u128,
        /// Maximum dense voxels authorized by the caller.
        maximum: usize,
    },
    /// A checked dense layout cannot be represented by this target.
    DenseVolumeOverflow {
        /// Operation that requested the dense layout.
        operation: &'static str,
        /// Checked per-axis dimensions before target-size conversion.
        dims: [u64; 3],
    },
    /// An operation requires at least one active occupancy voxel.
    EmptyOccupancy {
        /// Operation that refused the empty field.
        operation: &'static str,
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
            VoxelError::FrameMismatch {
                operation,
                left_voxel_size,
                right_voxel_size,
                left_origin,
                right_origin,
            } => write!(
                f,
                "{operation} refused mismatched voxel frames: left size/origin \
                 {left_voxel_size}/{left_origin:?}, right {right_voxel_size}/{right_origin:?}"
            ),
            VoxelError::CoordinateRange {
                operation,
                axis,
                min,
                max,
                halo,
            } => write!(
                f,
                "{operation} coordinate range on axis {axis} cannot include halo {halo}: \
                 [{min}, {max}]"
            ),
            VoxelError::VoxelBudgetExceeded {
                operation,
                required,
                maximum,
            } => write!(
                f,
                "{operation} requires {required} dense voxels, exceeding budget {maximum}"
            ),
            VoxelError::DenseVolumeOverflow { operation, dims } => write!(
                f,
                "{operation} dense dimensions {dims:?} cannot be represented on this target"
            ),
            VoxelError::EmptyOccupancy { operation } => {
                write!(f, "{operation} requires at least one active voxel")
            }
            VoxelError::Lattice { what } => write!(f, "degenerate lattice: {what}"),
            VoxelError::Cloud { what } => write!(f, "point-cloud query failed: {what}"),
            VoxelError::Graph { what } => write!(f, "graph round-trip failed: {what}"),
        }
    }
}

impl std::error::Error for VoxelError {}
