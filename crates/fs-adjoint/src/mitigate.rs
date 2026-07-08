//! NON-DIFFERENTIABLE MESHING MITIGATIONS (addendum Proposal 1, bead
//! bk0o.2; [F] — behind the `diff-mitigations` feature until its
//! Gauntlet tier is green). Remeshing is genuinely non-differentiable
//! across topology events; rather than pretend otherwise, three
//! mitigations apply IN ORDER OF PREFERENCE:
//!
//! 1. DIFFERENTIABILITY AS A ROUTING REQUIREMENT — the Rep Router's
//!    fitness gains a differentiability term (implemented at the cost
//!    -oracle seam: non-differentiable edges carry a penalty when a
//!    query requests gradients), so SDF/spline paths are PREFERRED and
//!    multi-representation becomes the asset that provides
//!    differentiable paths where mesh-locked pipelines have none;
//! 2. HADAMARD boundary-form shape derivatives that avoid mesh
//!    sensitivities entirely where applicable (the base crate's
//!    `hadamard` module, wired as the mesh-free path);
//! 3. where a remesh event is UNAVOIDABLE inside the differentiation
//!    path: the gradient is emitted with ESTIMATED color (Proposal 3)
//!    plus a declared DISCONTINUITY FLAG — never a silently-wrong
//!    verified gradient.

use fs_evidence::Color;
use fs_geom::{CostOracle, RoutePlan, RouteRefusal, RouteRequest, Router};
use std::collections::BTreeSet;

/// The routing penalty applied to non-differentiable edges when a query
/// requests gradients (large enough to dominate any realistic cost).
pub const NON_DIFF_PENALTY_S: f64 = 1e6;

/// A cost-oracle wrapper that adds the DIFFERENTIABILITY TERM to the
/// router's fitness: when `gradients_requested`, edges named in the
/// non-differentiable set report a huge measured cost, so the Pareto
/// planner prefers smooth (SDF/spline) paths. Recording is a no-op —
/// this wrapper is a PLANNING view; actuals go to the real oracle.
pub struct DiffAwareOracle<'a> {
    inner: &'a dyn CostOracle,
    non_differentiable: &'a BTreeSet<String>,
    gradients_requested: bool,
}

impl<'a> DiffAwareOracle<'a> {
    /// Wrap an oracle with the differentiability term.
    #[must_use]
    pub fn new(
        inner: &'a dyn CostOracle,
        non_differentiable: &'a BTreeSet<String>,
        gradients_requested: bool,
    ) -> Self {
        DiffAwareOracle {
            inner,
            non_differentiable,
            gradients_requested,
        }
    }
}

impl CostOracle for DiffAwareOracle<'_> {
    fn measured_cost_s(&self, edge: &str) -> Option<f64> {
        let base = self.inner.measured_cost_s(edge);
        if self.gradients_requested && self.non_differentiable.contains(edge) {
            Some(base.unwrap_or(0.0) + NON_DIFF_PENALTY_S)
        } else {
            base
        }
    }

    fn measured_error_abs(&self, edge: &str) -> Option<f64> {
        self.inner.measured_error_abs(edge)
    }

    fn record(&mut self, _edge: &str, _cost_s: f64, _error_abs: f64) {
        // Planning view only: actuals belong to the real oracle.
    }
}

/// How a gradient answer must be graded (Proposal 3 colors + the
/// discontinuity flag).
#[derive(Debug, Clone, PartialEq)]
pub enum GradientGrade {
    /// The chosen path is smooth: the gradient may carry whatever color
    /// its numerics earn (verified/validated per the certificate).
    Smooth {
        /// The route taken (edge names).
        route: Vec<String>,
    },
    /// A remesh/topology event is UNAVOIDABLE on every viable path:
    /// the gradient is Estimated with a declared discontinuity.
    EstimatedWithDiscontinuity {
        /// The route taken.
        route: Vec<String>,
        /// The offending edge(s).
        crossing: Vec<String>,
        /// The Proposal-3 color the gradient must carry.
        color: Color,
    },
}

impl GradientGrade {
    /// The discontinuity flag, if any.
    #[must_use]
    pub fn discontinuity(&self) -> Option<&[String]> {
        match self {
            GradientGrade::Smooth { .. } => None,
            GradientGrade::EstimatedWithDiscontinuity { crossing, .. } => Some(crossing),
        }
    }
}

/// Plan a conversion route UNDER a gradient request and grade the
/// resulting gradient honestly:
///
/// - if a differentiable path exists, the penalty steers the router
///   onto it and the answer is [`GradientGrade::Smooth`];
/// - if every viable path crosses a non-differentiable edge, the
///   cheapest such path is taken and the answer is
///   estimated-with-discontinuity — NEVER a silently-verified gradient
///   across a topology event (the review-round-3 boundary case).
///
/// # Errors
/// Propagates the router's structured refusals (no route at all).
pub fn plan_gradient_route(
    router: &Router,
    req: &RouteRequest,
    oracle: &dyn CostOracle,
    non_differentiable: &BTreeSet<String>,
) -> Result<(RoutePlan, GradientGrade), RouteRefusal> {
    let wrapped = DiffAwareOracle::new(oracle, non_differentiable, true);
    let plan = router.plan(req, &wrapped)?;
    let crossing: Vec<String> = plan
        .edges
        .iter()
        .filter(|e| non_differentiable.contains(*e))
        .cloned()
        .collect();
    let grade = if crossing.is_empty() {
        GradientGrade::Smooth {
            route: plan.edges.clone(),
        }
    } else {
        GradientGrade::EstimatedWithDiscontinuity {
            route: plan.edges.clone(),
            crossing: crossing.clone(),
            color: Color::Estimated {
                estimator: format!(
                    "gradient across non-differentiable edge(s) {crossing:?}: remesh/topology \
                     event inside the differentiation path"
                ),
                dispersion: f64::INFINITY,
            },
        }
    };
    Ok((plan, grade))
}

/// Grade a DIRECT (non-routed) differentiation path by its op names —
/// the same honesty rule for tape-level chains: any op in the declared
/// non-differentiable set forces estimated + flag.
#[must_use]
pub fn grade_ops(ops: &[&str], non_differentiable: &BTreeSet<String>) -> GradientGrade {
    let crossing: Vec<String> = ops
        .iter()
        .filter(|o| non_differentiable.contains(**o))
        .map(|o| (*o).to_string())
        .collect();
    if crossing.is_empty() {
        GradientGrade::Smooth {
            route: ops.iter().map(|o| (*o).to_string()).collect(),
        }
    } else {
        GradientGrade::EstimatedWithDiscontinuity {
            route: ops.iter().map(|o| (*o).to_string()).collect(),
            crossing: crossing.clone(),
            color: Color::Estimated {
                estimator: format!("non-differentiable op(s) {crossing:?} in the path"),
                dispersion: f64::INFINITY,
            },
        }
    }
}
