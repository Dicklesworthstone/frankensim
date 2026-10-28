//! Exact durable checkpoints for the single-load level-set optimizer.
//!
//! The fixed-run optimizer remains unchanged. This module reuses the same
//! corrected numerical evolution in a segment primitive that accepts an explicit
//! global iteration ordinal and current augmented-Lagrange multiplier, allowing
//! durable continuation without replaying prior geometry updates.

use crate::fim::{RedistanceAudit, redistance};
use crate::gridsdf::GridSdf;
use crate::optimize::{
    Cantilever, OptimizeReport, OptimizeSettings, SOLVER_MAX_ITERS, SOLVER_TOL,
    cantilever_support, fnv, invalid_input, mass_stiffness, material_volume,
    retain_cantilever_load_pad, validated_plane_strain_material,
};
use crate::topder::{NucleationEvent, nucleate, topological_derivative};
use crate::veloext::extend_velocity;
use crate::weno::{Velocity, advect, build_band};
use fs_cutfem::{
    BoundaryTraction, CutElasticity, ControlledElasticitySolution, CutFemError, EdgeBand, Quadtree,
};
use fs_material::IsotropicElastic;
use std::fmt::Write as _;

#[path = "optimize/engine_checkpoint.rs"]
mod engine;
#[path = "optimize/stateful.rs"]
mod stateful;

pub use stateful::OptimizeCheckpoint;

/// Cooperative interruption boundaries within one unpublished optimizer update.
///
/// CG counts are cumulative within the named solve, including correction passes.
/// Assembly and the other named stages are checked at their boundaries, not
/// preempted internally. This does not promise a wall-clock latency bound.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CheckpointStage {
    /// Before cloning the retained state or allocating discretization buffers.
    Prepare,
    /// Canonical solve of the retained design; zero also brackets assembly.
    InitialSolve(usize),
    /// Before sensitivity evaluation and Sobolev smoothing.
    Direction,
    /// Before advection, redistancing and optional hole nucleation.
    Evolution,
    /// Canonical solve of the proposed design; zero also brackets assembly.
    CandidateSolve(usize),
    /// Both solves succeeded, but the candidate has not replaced durable state.
    Publish,
}
