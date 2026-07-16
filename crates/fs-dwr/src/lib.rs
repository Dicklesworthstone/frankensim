//! fs-dwr — dual-weighted-residual goal-oriented adaptivity (plan
//! §8.6 [F], bead tfz.23): in a design-optimization system, "accurate
//! simulation" is not the goal — ACCURATE OBJECTIVE AND GRADIENT is,
//! and DWR is the estimator that knows the difference. The adjoint
//! solution weights local residuals by their influence on the quantity
//! being optimized; refinement follows the product.
//!
//! Layer: L3. One signal, four mechanisms:
//! - [`estimate`]: the DWR core on the fs-cutfem stack — primal and
//!   adjoint solves on the coarse active space, an ENRICHED adjoint
//!   (one-level-finer solve; the documented higher-resolution enrichment
//!   option, with patch recovery the recorded alternative), signed per-cell
//!   indicators from the full residual (interior + Nitsche interface terms),
//!   and effectivity against known-truth fixtures.
//! - [`mark`]: Dörfler fixed-energy marking with DETERMINISTIC
//!   tie-breaking (indicator desc, cell key asc) — two runs mark
//!   bitwise-identically.
//! - [`adapt`]: the octree h-refinement loop (mechanism 1) — solve →
//!   estimate → mark → constraint-aware split admission → rebalance.
//!   A scalar mark may deliberately raise the declared cut-band target; a
//!   mismatch-only repair at the existing target is deferred in deterministic
//!   indicator order, so local accuracy does not silently spend a global
//!   refinement budget. Accuracy-per-DOF curves are the ledgered output.
//! - [`aniso`]: anisotropic METRIC SYNTHESIS (mechanism 2) — recovered
//!   Hessians weighted by adjoint magnitude, normalized to a target
//!   complexity, exported as an fs-mesh `MetricField`-shaped tensor
//!   per cell (the conformance battery drives fs-mesh's remesher with
//!   it end-to-end).
//! - [`tiles`]: wavelet-style Haar coefficient THRESHOLDING
//!   (mechanism 4) — compression-as-adaptivity with DWR-weighted local
//!   budgets: spend accuracy where the adjoint says the goal cannot
//!   see it.
//! - [`hvsp`]: the h-vs-p DECISION signal (mechanism 3) — smoothness
//!   classification of where p-enrichment beats h-refinement;
//!   EXECUTING local p awaits the high-order FEEC families (recorded
//!   no-claim).
//!
//! Adjoint doctrine: the goal problems here are symmetric, so each active
//! space needs one transposed(=same-operator) adjoint solve — exactly
//! fs-adjoint's implicit-function-theorem discipline; nonsymmetric operators
//! wire through fs-adjoint::ift when they arrive. Scalar two-level DWR solves
//! coarse and enriched adjoints because the active CutFEM spaces are not
//! nodally nested.
//!
//! Vector elasticity compliance uses the same doctrine in [`elasticity`].
//! Its bulk/Nitsche/traction residual is explicit, while its ghost term is a
//! named coarse consistent-energy correction because the active CutFEM spaces
//! are non-nested.

pub mod adapt;
pub mod aniso;
pub mod elasticity;
pub mod estimate;
pub mod hvsp;
pub mod mark;
pub mod tiles;

pub use adapt::{AdaptStep, adapt_loop};
pub use aniso::synthesize_metric;
pub use elasticity::{
    ElasticityDwrEstimate, ElasticityGhostMethod, ElasticityResidualTerms,
    estimate_elasticity_compliance,
};
pub use estimate::{DwrEstimate, GoalContext, estimate, goal_value};
pub use hvsp::{Decision, h_vs_p};
pub use mark::dorfler;
pub use tiles::{ThresholdOutcome, haar_threshold};

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[cfg(test)]
mod tests {
    #[test]
    fn version_is_stamped() {
        assert!(!super::VERSION.is_empty());
    }
}
