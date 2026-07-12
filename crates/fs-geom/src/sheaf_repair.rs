//! SHEAF REPAIR (patch Rev L, bead wqd.14; [M] — behind the
//! `sheaf-repair` feature until certifier trials pass): upgrade the sheaf
//! machinery from diagnosis to CONSTRUCTIVE repair. The interface
//! mismatch cochain Hodge-decomposes as `m = exact ⊕ coexact ⊕ harmonic`
//! over the patch-adjacency complex, and each component has an
//! INTERPRETATION CONTRACT:
//!
//! - EXACT (`δ⁰c`): repairable by local patch/gauge adjustments — the
//!   auto-repair target, bounded by each chart's declared error budget
//!   (a repair must never silently distort geometry beyond budget);
//! - COEXACT (`δ¹ᵀw`): circulation-like inconsistency around triple
//!   junctions — usually orientation/trace errors, which points at
//!   CONVERTER bugs, not geometry;
//! - HARMONIC (the remainder): a true topological obstruction — no local
//!   fix exists, reported honestly with the minimal interface cut-set.
//!
//! Repairs are PROPOSALS, ranked with expected post-repair norms and
//! optional Rep-Router reroute costs; they apply only under an explicit
//! budget.

use crate::router::{CostOracle, RoutePlanError, RouteRequest, Router};
use crate::sheaf::SheafComplex;
use std::fmt::Write as _;

/// The complex skeleton the decomposition runs over (extractable from a
/// [`SheafComplex`] or built directly for controlled fixtures).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SheafSkeleton {
    /// Patch count.
    pub n_patches: usize,
    /// Interfaces as (u, v) with u < v (edge k orients u → v).
    pub edges: Vec<(usize, usize)>,
    /// Triple junctions (a, b, c) sorted; boundary = +e_ab + e_bc − e_ac.
    pub triangles: Vec<(usize, usize, usize)>,
}

impl SheafSkeleton {
    /// Extract the skeleton of a built complex.
    #[must_use]
    pub fn of(complex: &SheafComplex) -> SheafSkeleton {
        SheafSkeleton {
            n_patches: complex.n_patches,
            edges: complex.interfaces.iter().map(|i| i.patches).collect(),
            triangles: complex.triples.iter().map(|t| t.patches).collect(),
        }
    }

    fn edge_index(&self, a: usize, b: usize) -> Option<usize> {
        let key = (a.min(b), a.max(b));
        self.edges.iter().position(|&e| e == key)
    }

    /// Apply δ⁰ to a vertex cochain: `(δ⁰c)_e = c_v − c_u`.
    #[must_use]
    pub fn d0(&self, c: &[f64]) -> Vec<f64> {
        self.edges.iter().map(|&(u, v)| c[v] - c[u]).collect()
    }

    /// Apply δ⁰ᵀ to an edge cochain.
    #[must_use]
    pub fn d0t(&self, m: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0f64; self.n_patches];
        for (k, &(u, v)) in self.edges.iter().enumerate() {
            out[u] -= m[k];
            out[v] += m[k];
        }
        out
    }

    /// Apply δ¹ to an edge cochain: signed sum around each triangle.
    #[must_use]
    pub fn d1(&self, m: &[f64]) -> Vec<f64> {
        self.triangles
            .iter()
            .map(|&(a, b, c)| {
                let eab = self.edge_index(a, b).expect("triangle implies edge");
                let ebc = self.edge_index(b, c).expect("triangle implies edge");
                let eac = self.edge_index(a, c).expect("triangle implies edge");
                m[eab] + m[ebc] - m[eac]
            })
            .collect()
    }

    /// Apply δ¹ᵀ to a triangle cochain.
    #[must_use]
    pub fn d1t(&self, w: &[f64]) -> Vec<f64> {
        let mut out = vec![0.0f64; self.edges.len()];
        for (t, &(a, b, c)) in self.triangles.iter().enumerate() {
            let eab = self.edge_index(a, b).expect("triangle implies edge");
            let ebc = self.edge_index(b, c).expect("triangle implies edge");
            let eac = self.edge_index(a, c).expect("triangle implies edge");
            out[eab] += w[t];
            out[ebc] += w[t];
            out[eac] -= w[t];
        }
        out
    }
}

/// The Hodge decomposition of an edge mismatch cochain.
#[derive(Debug, Clone, PartialEq)]
pub struct HodgeSplit {
    /// The exact (coboundary) component `δ⁰c`.
    pub exact: Vec<f64>,
    /// The vertex potential `c` (gauge offsets; `c[0]` pinned to 0).
    pub potential: Vec<f64>,
    /// The coexact component `δ¹ᵀw`.
    pub coexact: Vec<f64>,
    /// The harmonic remainder.
    pub harmonic: Vec<f64>,
    /// Energy fractions (exact, coexact, harmonic) of ‖m‖².
    pub fractions: (f64, f64, f64),
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn norm2(a: &[f64]) -> f64 {
    dot(a, a)
}

/// Least squares `min ‖m − A x‖²` via Gauss–Seidel on the normal
/// equations, with `apply`/`apply_t` as the operator (small complexes;
/// deterministic sweep order; component 0 optionally pinned).
fn least_squares(
    m: &[f64],
    n_unknowns: usize,
    apply: impl Fn(&[f64]) -> Vec<f64>,
    apply_t: impl Fn(&[f64]) -> Vec<f64>,
    pin_first: bool,
) -> Vec<f64> {
    let mut x = vec![0.0f64; n_unknowns];
    let rhs = apply_t(m);
    // Diagonal of AᵀA via unit vectors (small n — fine and exact).
    let mut diag = vec![0.0f64; n_unknowns];
    for (i, d) in diag.iter_mut().enumerate() {
        let mut e = vec![0.0f64; n_unknowns];
        e[i] = 1.0;
        *d = norm2(&apply(&e));
    }
    for _ in 0..400 {
        for i in 0..n_unknowns {
            if pin_first && i == 0 {
                continue;
            }
            if diag[i] <= 0.0 {
                continue;
            }
            // Residual of the normal equations at coordinate i.
            let ax = apply(&x);
            let grad_i = {
                let atax = apply_t(&ax);
                atax[i] - rhs[i]
            };
            x[i] -= grad_i / diag[i];
        }
    }
    x
}

/// Hodge-decompose an edge cochain over a skeleton. Verified against a
/// dense oracle in conformance; orthogonality residuals are returned by
/// the components themselves (`d0t(harmonic) ≈ 0`, `d1(harmonic) ≈ 0`).
#[must_use]
pub fn hodge_decompose(skeleton: &SheafSkeleton, m: &[f64]) -> HodgeSplit {
    assert_eq!(m.len(), skeleton.edges.len(), "cochain size");
    // Exact: project onto im δ⁰.
    let c = least_squares(
        m,
        skeleton.n_patches,
        |x| skeleton.d0(x),
        |y| skeleton.d0t(y),
        true,
    );
    let exact = skeleton.d0(&c);
    let r1: Vec<f64> = m.iter().zip(&exact).map(|(a, b)| a - b).collect();
    // Coexact: project the remainder onto im δ¹ᵀ.
    let coexact = if skeleton.triangles.is_empty() {
        vec![0.0; m.len()]
    } else {
        let w = least_squares(
            &r1,
            skeleton.triangles.len(),
            |x| skeleton.d1t(x),
            |y| skeleton.d1(y),
            false,
        );
        skeleton.d1t(&w)
    };
    let harmonic: Vec<f64> = r1.iter().zip(&coexact).map(|(a, b)| a - b).collect();
    let total = norm2(m).max(f64::MIN_POSITIVE);
    HodgeSplit {
        fractions: (
            norm2(&exact) / total,
            norm2(&coexact) / total,
            norm2(&harmonic) / total,
        ),
        exact,
        potential: c,
        coexact,
        harmonic,
    }
}

/// One ranked repair proposal (the agent-facing format).
#[derive(Debug, Clone, PartialEq)]
pub struct RepairProposal {
    /// What to do, concretely.
    pub action: String,
    /// Expected post-repair worst interface mismatch.
    pub expected_post_norm: f64,
    /// Cost estimate in seconds (router-modeled where applicable).
    pub cost_s: f64,
}

/// The repair verdict for one model.
#[derive(Debug, Clone, PartialEq)]
pub struct RepairPlan {
    /// The decomposition driving the plan.
    pub split: HodgeSplit,
    /// Ranked proposals (best first).
    pub proposals: Vec<RepairProposal>,
    /// Gauge offsets the auto-repair would apply (per patch).
    pub gauge: Vec<f64>,
    /// True when the exact-component repair fits EVERY patch budget.
    pub auto_repairable: bool,
    /// Interfaces in the harmonic support (the minimal cut-set needing
    /// topology-level changes) with their magnitudes.
    pub obstruction_cutset: Vec<((usize, usize), f64)>,
    /// Structured reason an optional router alternative could not be planned.
    /// `None` means no reroute was requested or a proposal was produced.
    pub reroute_error: Option<RoutePlanError>,
}

/// Threshold below which a component is treated as absent (fractions).
pub const COMPONENT_FLOOR: f64 = 1e-6;

/// Build the repair plan: decompose, interpret, rank. `budgets` is each
/// patch's declared error budget — the exact-component gauge repair is
/// only auto-appliable when |offset| stays within it for EVERY patch
/// (a repair must never silently distort geometry beyond budget).
/// `reroute` optionally consults the Rep Router for a conversion-based
/// alternative for the worst-offending patch.
#[must_use]
pub fn plan_repair(
    skeleton: &SheafSkeleton,
    mismatch: &[f64],
    budgets: &[f64],
    reroute: Option<(&Router, &dyn CostOracle, &RouteRequest)>,
) -> RepairPlan {
    let split = hodge_decompose(skeleton, mismatch);
    let residual_after_exact: Vec<f64> = mismatch
        .iter()
        .zip(&split.exact)
        .map(|(m, e)| m - e)
        .collect();
    let expected_after_gauge = residual_after_exact
        .iter()
        .fold(0.0f64, |a, &b| a.max(b.abs()));
    let auto_repairable = split.fractions.0 > COMPONENT_FLOOR
        && split
            .potential
            .iter()
            .zip(budgets)
            .all(|(off, budget)| off.abs() <= *budget);
    let mut proposals = Vec::new();
    if split.fractions.0 > COMPONENT_FLOOR {
        proposals.push(gauge_proposal(
            &split,
            auto_repairable,
            expected_after_gauge,
        ));
    }
    if split.fractions.1 > COMPONENT_FLOOR {
        proposals.push(coexact_proposal(skeleton, mismatch, &split));
    }
    let mut obstruction_cutset: Vec<((usize, usize), f64)> = skeleton
        .edges
        .iter()
        .zip(&split.harmonic)
        .filter(|(_, h)| h.abs() > COMPONENT_FLOOR.sqrt())
        .map(|(&e, &h)| (e, h.abs()))
        .collect();
    obstruction_cutset.sort_by(|a, b| b.1.total_cmp(&a.1));
    if !obstruction_cutset.is_empty() {
        proposals.push(RepairProposal {
            action: format!(
                "harmonic obstruction: NO local fix exists; topology-level change needed \
                 on the cut-set {:?}",
                obstruction_cutset
                    .iter()
                    .map(|(e, _)| *e)
                    .collect::<Vec<_>>()
            ),
            expected_post_norm: obstruction_cutset[0].1,
            cost_s: f64::INFINITY,
        });
    }
    let mut reroute_error = None;
    if let Some((router, oracle, req)) = reroute {
        match router.plan(req, oracle) {
            Ok(route) => proposals.push(RepairProposal {
                action: format!(
                    "reroute worst patch {} -> {} via [{}] (router-planned alternative chart)",
                    req.from,
                    req.to,
                    route.edges().join(", ")
                ),
                expected_post_norm: route.composed_abs_error(),
                cost_s: route.predicted_cost_s(),
            }),
            Err(error) => reroute_error = Some(error),
        }
    }
    proposals.sort_by(|a, b| {
        a.expected_post_norm
            .total_cmp(&b.expected_post_norm)
            .then(a.cost_s.total_cmp(&b.cost_s))
    });
    RepairPlan {
        gauge: split.potential.clone(),
        split,
        proposals,
        auto_repairable,
        obstruction_cutset,
        reroute_error,
    }
}

/// The exact-component proposal: the concrete per-patch gauge projection.
fn gauge_proposal(split: &HodgeSplit, auto_repairable: bool, expected: f64) -> RepairProposal {
    let worst = split
        .potential
        .iter()
        .enumerate()
        .max_by(|a, b| a.1.abs().total_cmp(&b.1.abs()))
        .map_or(0, |(i, _)| i);
    let mut action = format!(
        "project patch P{worst} gauge by {:+.3e} (exact-component section \
         projection; offsets per patch: [",
        split.potential[worst]
    );
    for (i, off) in split.potential.iter().enumerate() {
        if i > 0 {
            action.push_str(", ");
        }
        let _ = write!(action, "{off:+.3e}");
    }
    action.push_str("])");
    if !auto_repairable {
        let _ = write!(action, " — EXCEEDS a patch budget; needs acceptance");
    }
    RepairProposal {
        action,
        expected_post_norm: expected,
        cost_s: 0.001, // local gauge arithmetic
    }
}

/// The coexact-component proposal: converter-side diagnosis, localized
/// to the worst triple junction.
fn coexact_proposal(
    skeleton: &SheafSkeleton,
    mismatch: &[f64],
    split: &HodgeSplit,
) -> RepairProposal {
    let d1m = skeleton.d1(mismatch);
    let worst_tri = skeleton
        .triangles
        .iter()
        .enumerate()
        .max_by(|a, b| d1m[a.0].abs().total_cmp(&d1m[b.0].abs()))
        .map(|(_, t)| *t);
    RepairProposal {
        action: format!(
            "coexact circulation detected around triple junction {worst_tri:?}: check \
             CONVERTER orientation/trace conventions (not a geometry edit)"
        ),
        expected_post_norm: (norm2(&split.harmonic) + norm2(&split.exact)).sqrt(),
        cost_s: 0.0,
    }
}

/// Apply the gauge repair to an edge cochain (the constructive step):
/// `m ← m − δ⁰c`. Idempotent on repaired models (the residual has no
/// exact component left).
#[must_use]
pub fn apply_gauge(skeleton: &SheafSkeleton, mismatch: &[f64], gauge: &[f64]) -> Vec<f64> {
    let correction = skeleton.d0(gauge);
    mismatch
        .iter()
        .zip(&correction)
        .map(|(m, c)| m - c)
        .collect()
}
