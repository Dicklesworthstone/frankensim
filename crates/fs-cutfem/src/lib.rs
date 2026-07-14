//! fs-cutfem — CutFEM on SDFs (plan §8.1 frontend 2, bead tfz.8): THE
//! MARQUEE BRIDGE. Level-set/SDF geometry simulated with FEM-grade
//! accuracy and ZERO MESHING — the property that lets topology
//! optimization loops rerun physics per design perturbation without a
//! mesher in the loop (Bet 3's biggest single payoff).
//!
//! Layer: L3. The 2D quadtree instantiation of the octree design — the
//! quadtree IS the octree restricted to two axes, sharing the
//! FrankenVDB dyadic-tile alignment (cells at level ℓ are 2⁻ℓ dyadic
//! boxes; a leaf at the tile depth is one FrankenVDB leaf face). The
//! 3D octree instantiation is a recorded no-claim (CONTRACT.md).
//!
//! The pipeline, module by module:
//! - [`sdf`]: the [`CutSdf`] trait — a level-set function plus a
//!   CERTIFIED per-cell enclosure (fs-ivl outward-rounded intervals).
//! - [`grid`]: [`Quadtree`] background grids — 2:1-balanced dyadic
//!   refinement, hanging-node constraints handled in the element space
//!   (non-uniform resolution comes FREE from the tree).
//! - [`quad`]: cut-cell bulk and interface quadrature — recursive
//!   subdivision with certified classification at every level, exact
//!   crossings by bisection, degree-2-exact polygon rules; error
//!   control by subdivision depth.
//! - [`fem`]: Q1 spaces on active cells, Nitsche weak embedded
//!   Dirichlet conditions, GHOST-PENALTY stabilization (small-cut
//!   conditioning), assembly to fs-sparse, fs-solver CG.
//! - [`elastic`]: vector Q1 small-strain elasticity over the same cut
//!   rules, with symmetric vector Nitsche terms, cut-independent
//!   penalty scaling, componentwise ghost stabilization, and an
//!   optional fs-adjoint VJP registration.
//! - [`agg`]: aggregated-element fallback — small-cut cells lend their
//!   ill-supported DOFs to a well-cut anchor by polynomial extension
//!   (belt + suspenders with ghost penalty; policy documented).
//! - [`cond`]: dense conditioning probe for the conformance batteries
//!   (eigenvalue-verified conditioning independence of cut fraction).
//!
//! Determinism: BTree-ordered traversal everywhere, straight-line IEEE
//! arithmetic, fs-solver's deterministic CG — bit-deterministic across
//! runs by construction. Ambition tag: the plan marks CutFEM-on-SDF
//! [F]; per the crate-granular form of the gating rule (the fs-feec
//! precedent) the frontier surface ships as this standalone crate and
//! consumers opt in by depending on it.

pub mod agg;
pub mod cond;
pub mod elastic;
pub mod fem;
pub mod grid;
pub mod quad;
pub mod sdf;

pub use agg::AggPolicy;
pub use cond::{CondReport, condition_estimate};
pub use elastic::MAX_PLANE_STRAIN_STIFFNESS_RATIO;
pub use elastic::{
    BoundaryTraction, CutElasticity, CutElasticityOperator, CutElasticitySolution, DesignBoxEdge,
    EdgeBand,
};
#[cfg(feature = "adjoint-vjp")]
pub use elastic::{
    ELASTICITY_APPLY_VJP_OP, elasticity_apply_vjp_key, register_elasticity_apply_vjp,
};
pub use fem::{BuildStats, CellClass, FemParams, ScalarSample, Solution, Space};
pub use grid::{CellKey, NodeKey, Quadtree};
pub use quad::{CutRules, cut_cell_rules};
pub use sdf::{Circle, CutSdf, HalfPlane};

/// Crate version, re-exported for provenance stamping (the Five
/// Explicits' "versions" pillar reaches down to individual crates).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Teaching errors (P10): every refusal names the violated assumption
/// and the repair.
#[derive(Debug, Clone, PartialEq)]
pub enum CutFemError {
    /// No active cell survived classification — the level set does not
    /// enter the background domain (or the domain is entirely outside).
    /// Repair: check the SDF sign convention (negative = inside) and
    /// the grid extent.
    EmptyDomain,
    /// A cut cell has an active face-neighbor at a DIFFERENT tree
    /// level while ghost stabilization is enabled. Ghost-penalty faces
    /// are assembled between equal-level cells only; the interface band
    /// must be uniformly refined. Ghost-free paths do not raise this error.
    /// Repair: call [`Quadtree::refine_toward_interface`] with the
    /// tree's max level before building the space.
    CutBandNotUniform {
        /// The offending cut cell (level, i, j).
        cell: CellKey,
        /// Its differently-leveled active neighbor.
        neighbor: CellKey,
    },
    /// A scalar CutFEM parameter cannot define the documented finite,
    /// nonnegative stabilization path.
    InvalidFemInput {
        /// Actionable description of the invalid field/value.
        what: String,
    },
    /// A vector-elasticity parameter or callback returned a value that
    /// cannot define a finite coercive discrete problem.
    InvalidElasticityInput {
        /// Actionable description of the invalid field/value.
        what: String,
    },
    /// Aggregation found no well-supported anchor cell within the
    /// search radius of a small-cut node. Repair: refine the interface
    /// band (isolated sliver islands shrink), or raise
    /// [`AggPolicy::good_fraction`] tolerance, or enable ghost penalty
    /// instead.
    AggregationNoAnchor {
        /// The orphaned node (lattice coordinates at max level).
        node: (u32, u32),
    },
    /// A constraint chain (hanging/aggregation) failed to terminate —
    /// this indicates a corrupted constraint graph, not user error.
    ConstraintCycle {
        /// The node where the cycle was detected.
        node: (u32, u32),
    },
    /// The CG solve did not reach the residual gate. Repair: raise
    /// `solver_max_iters`, enable ghost penalty or aggregation (the
    /// usual cause is a small-cut-degenerate spectrum), or loosen
    /// `solver_tol`.
    SolveNotConverged {
        /// Iterations performed.
        iters: usize,
        /// Relative residual reached.
        rel_residual: f64,
    },
}

impl core::fmt::Display for CutFemError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            CutFemError::EmptyDomain => write!(
                f,
                "no active cells: the level set never enters the grid \
                 (check the negative-inside sign convention and extent)"
            ),
            CutFemError::CutBandNotUniform { cell, neighbor } => write!(
                f,
                "cut cell {cell:?} has active neighbor {neighbor:?} at a \
                 different level; refine the interface band uniformly \
                 (Quadtree::refine_toward_interface) before building"
            ),
            CutFemError::InvalidFemInput { what } => {
                write!(f, "invalid scalar CutFEM input: {what}")
            }
            CutFemError::InvalidElasticityInput { what } => {
                write!(f, "invalid vector CutFEM elasticity input: {what}")
            }
            CutFemError::AggregationNoAnchor { node } => write!(
                f,
                "aggregation found no well-cut anchor near node {node:?}; \
                 refine the interface band or enable ghost penalty"
            ),
            CutFemError::ConstraintCycle { node } => write!(
                f,
                "constraint chain through node {node:?} does not terminate \
                 (corrupted constraint graph)"
            ),
            CutFemError::SolveNotConverged {
                iters,
                rel_residual,
            } => write!(
                f,
                "CG stalled at relative residual {rel_residual:.3e} after \
                 {iters} iterations; enable ghost penalty/aggregation or \
                 raise solver_max_iters"
            ),
        }
    }
}

impl std::error::Error for CutFemError {}

#[cfg(test)]
mod tests {
    #[test]
    fn version_is_stamped() {
        assert!(!super::VERSION.is_empty());
    }
}
