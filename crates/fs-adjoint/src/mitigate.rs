//! NON-DIFFERENTIABLE MESHING MITIGATIONS (addendum Proposal 1, bead
//! bk0o.2; [F] — behind the `diff-mitigations` feature until its
//! Gauntlet tier is green). Remeshing is genuinely non-differentiable
//! across topology events; rather than pretend otherwise, three
//! mitigations apply IN ORDER OF PREFERENCE:
//!
//! 1. DIFFERENTIABILITY AS A ROUTING REQUIREMENT — a gradient query first
//!    plans over the subgraph containing only differentiable edges. The
//!    original graph is considered only when no smooth route satisfies the
//!    request's real cost/error budgets, so SDF/spline paths are PREFERRED
//!    without disguising differentiability policy as wall-clock cost;
//! 2. HADAMARD boundary-form shape derivatives that avoid mesh
//!    sensitivities entirely where applicable (the base crate's
//!    `hadamard` module, wired as the mesh-free path);
//! 3. where a remesh event is UNAVOIDABLE inside the differentiation
//!    path: the gradient is emitted with ESTIMATED color (Proposal 3)
//!    plus a declared DISCONTINUITY FLAG — never a silently-wrong
//!    verified gradient.

use fs_evidence::Color;
use fs_geom::{CostOracle, RoutePlan, RoutePlanError, RouteRequest, Router};
use std::collections::BTreeSet;

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
/// - if a differentiable path satisfies the request's original budgets, the
///   deterministic router selects within that smooth-only subgraph and the
///   answer is [`GradientGrade::Smooth`];
/// - otherwise the original graph is planned under the unchanged budgets; if
///   its winning path crosses a non-differentiable edge, the answer is
///   estimated-with-discontinuity — NEVER a silently-verified gradient
///   across a topology event (the review-round-3 boundary case).
///
/// # Errors
/// Propagates malformed requests, invalid oracle evidence, arithmetic failures,
/// and the router's structured refusal when no route is admissible.
pub fn plan_gradient_route(
    router: &Router,
    req: &RouteRequest,
    oracle: &dyn CostOracle,
    non_differentiable: &BTreeSet<String>,
) -> Result<(RoutePlan, GradientGrade), RoutePlanError> {
    let plan = router
        .plan_prefer_edge_filter(req, oracle, |spec| !non_differentiable.contains(&spec.name))?;
    let route = plan.edges().to_vec();
    let crossing: Vec<String> = route
        .iter()
        .filter(|e| non_differentiable.contains(*e))
        .cloned()
        .collect();
    let grade = if crossing.is_empty() {
        GradientGrade::Smooth { route }
    } else {
        let estimator = format!(
            "gradient across non-differentiable edge(s) {crossing:?}: remesh/topology \
             event inside the differentiation path"
        );
        GradientGrade::EstimatedWithDiscontinuity {
            route,
            crossing,
            color: Color::Estimated {
                estimator,
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
        let estimator = format!("non-differentiable op(s) {crossing:?} in the path");
        GradientGrade::EstimatedWithDiscontinuity {
            route: ops.iter().map(|o| (*o).to_string()).collect(),
            crossing,
            color: Color::Estimated {
                estimator,
                dispersion: f64::INFINITY,
            },
        }
    }
}
