//! Infeasibility DIAGNOSIS: elastic-relaxation solves, minimal unsat
//! cores by deletion filtering, and RANKED repairs with feasibility
//! estimates — the machinery that turns "optimizer failed" into a
//! design conversation.
//!
//! The elastic solve minimizes total hinge violation `Σ max(gᵢ, 0)`
//! over a domain box with multi-start projected subgradient descent
//! (deterministic LCG starts). Feasible ⟺ the elastic optimum's total
//! violation is ~0. The unsat core starts from the elastic support
//! (violated constraints at the optimum) and is refined by the
//! DELETION FILTER. The support is verified jointly infeasible before
//! filtering; otherwise the seed expands to the full, already-proven
//! infeasible set. The result is MINIMAL: dropping ANY member restores
//! feasibility — a property the conformance battery checks against
//! brute-force enumeration.

use crate::{ConError, ConstraintSpec, DomainError, DomainRangeError, push_json_string, scalar_at};
use fs_exec::Cx;
use fs_opt::{Manifold, Problem};

/// Per-component design-domain box.
#[derive(Debug, Clone, PartialEq)]
pub struct DomainBox {
    /// `(lo, hi)` per component of the sole `Rn` design variable. Admission
    /// requires exact dimension, finite ordered endpoints, and finite spans;
    /// `lo == hi` denotes a valid fixed coordinate.
    pub ranges: Vec<(f64, f64)>,
}

/// The elastic-relaxation solve's outcome.
#[derive(Debug, Clone)]
pub struct ElasticReport {
    /// The minimizer of total violation.
    pub x: Vec<f64>,
    /// Total hinge violation at the optimum (~0 ⟺ feasible).
    pub total_violation: f64,
    /// Per-constraint violations at the optimum.
    pub violations: Vec<f64>,
    /// Objective evaluations spent.
    pub evals: u64,
}

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0
    }

    fn unit(&mut self) -> f64 {
        ((self.next() >> 11) as f64) / (1u64 << 53) as f64
    }
}

/// Feasibility tolerance for the elastic optimum.
pub(crate) const FEAS_TOL: f64 = 1e-6;

/// Defensive penalty for an evaluator backend that carries a NON-FINITE scalar
/// instead of returning fs-opt's current typed `EvalNonFinite` refusal. It is
/// finite so a raw `NaN.max(0.0)` can never be dropped into false feasibility.
const NONFINITE_PENALTY: f64 = 1e30;

const CANCELLATION_STRIDE: usize = 64;

fn checkpoint(cx: &Cx<'_>) -> Result<(), ConError> {
    cx.checkpoint()
        .map_err(|_| ConError::Eval(fs_opt::OptError::Cancelled))
}

fn violation_contribution(value: f64) -> f64 {
    if value.is_finite() {
        value.max(0.0)
    } else {
        NONFINITE_PENALTY
    }
}

fn checked_total_violation(total: f64, contribution: f64) -> Result<f64, ConError> {
    let next = total + contribution;
    if next.is_finite() {
        Ok(next)
    } else {
        Err(ConError::BadParam {
            what: "elastic total violation",
            value: next,
        })
    }
}

fn validate_domain(problem: &Problem, domain: &DomainBox, cx: &Cx<'_>) -> Result<usize, ConError> {
    if problem.vars().len() != 1 {
        return Err(ConError::InvalidDomain(DomainError::HostVariableCount {
            got: problem.vars().len(),
        }));
    }
    let variable = &problem.vars()[0];
    let Manifold::Rn { dim } = variable.manifold else {
        return Err(ConError::InvalidDomain(DomainError::HostVariableManifold {
            got: variable.manifold,
        }));
    };
    let expected = usize::try_from(dim).map_err(|_| {
        ConError::InvalidDomain(DomainError::PointDimensionUnrepresentable { declared: dim })
    })?;
    if domain.ranges.len() != expected {
        return Err(ConError::InvalidDomain(DomainError::DimensionMismatch {
            expected,
            got: domain.ranges.len(),
        }));
    }
    for (axis, &(lo, hi)) in domain.ranges.iter().enumerate() {
        if axis % CANCELLATION_STRIDE == 0 {
            checkpoint(cx)?;
        }
        if !lo.is_finite() || !hi.is_finite() {
            return Err(ConError::InvalidDomain(DomainError::InvalidRange {
                axis,
                lo,
                hi,
                reason: DomainRangeError::NonFiniteEndpoint,
            }));
        }
        if lo > hi {
            return Err(ConError::InvalidDomain(DomainError::InvalidRange {
                axis,
                lo,
                hi,
                reason: DomainRangeError::Reversed,
            }));
        }
        if !(hi - lo).is_finite() {
            return Err(ConError::InvalidDomain(DomainError::InvalidRange {
                axis,
                lo,
                hi,
                reason: DomainRangeError::UnrepresentableSpan,
            }));
        }
    }
    Ok(expected)
}

/// Minimize `Σ max(gᵢ(x), 0)` over the box: multi-start projected
/// subgradient descent (deterministic). Small-fixture machinery — the
/// production restoration solver is a later ASCENT bead.
///
/// # Errors
/// [`ConError::InvalidDomain`] before allocation/evaluation, evaluation
/// teaching errors carried through, or cancellation at a documented poll.
pub fn elastic_solve(
    problem: &Problem,
    specs: &[ConstraintSpec],
    domain: &DomainBox,
    skip: &[usize],
    cx: &Cx<'_>,
) -> Result<ElasticReport, ConError> {
    let n = validate_domain(problem, domain, cx)?;
    checkpoint(cx)?;
    let mut active = Vec::new();
    for i in 0..specs.len() {
        if i % CANCELLATION_STRIDE == 0 {
            checkpoint(cx)?;
        }
        if !skip.contains(&i) {
            active.push(i);
        }
    }
    let mut evals = 0u64;
    let total = |x: &[f64], evals: &mut u64| -> Result<f64, ConError> {
        let mut t = 0.0;
        for (ordinal, &i) in active.iter().enumerate() {
            if ordinal % CANCELLATION_STRIDE == 0 {
                checkpoint(cx)?;
            }
            let gi = scalar_at(problem, specs[i].node, x)?;
            t = checked_total_violation(t, violation_contribution(gi))?;
            *evals += 1;
        }
        Ok(t)
    };
    let clamp = |x: &mut [f64]| {
        for (v, &(lo, hi)) in x.iter_mut().zip(&domain.ranges) {
            *v = v.clamp(lo, hi);
        }
    };
    let mut rng = Lcg(0x1001_2026_0707_0001);
    let mut best_x: Vec<f64> = domain
        .ranges
        .iter()
        .map(|&(lo, hi)| f64::midpoint(lo, hi))
        .collect();
    let mut best_v = total(&best_x, &mut evals)?;
    for start in 0..8 {
        checkpoint(cx)?;
        let mut x: Vec<f64> = if start == 0 {
            best_x.clone()
        } else {
            domain
                .ranges
                .iter()
                .map(|&(lo, hi)| lo + (hi - lo) * rng.unit())
                .collect()
        };
        let mut v = total(&x, &mut evals)?;
        let diam: f64 = domain
            .ranges
            .iter()
            .map(|&(lo, hi)| hi - lo)
            .fold(0.0, f64::max);
        for step in 0..300 {
            checkpoint(cx)?;
            if v <= FEAS_TOL {
                break;
            }
            // FD subgradient of the hinge sum.
            let h = 1e-6 * diam.max(1.0);
            let mut g = vec![0.0; n];
            for (k, gk) in g.iter_mut().enumerate() {
                if k % CANCELLATION_STRIDE == 0 {
                    checkpoint(cx)?;
                }
                let mut xp = x.clone();
                xp[k] += h;
                clamp(&mut xp);
                let mut xm = x.clone();
                xm[k] -= h;
                clamp(&mut xm);
                *gk = (total(&xp, &mut evals)? - total(&xm, &mut evals)?)
                    / (xp[k] - xm[k]).max(1e-300);
            }
            let gn = g.iter().map(|v| v * v).sum::<f64>().sqrt();
            if gn < 1e-14 {
                break;
            }
            let lr = 0.3 * diam / (1.0 + f64::from(step) * 0.05) / gn;
            for (xv, gv) in x.iter_mut().zip(&g) {
                *xv -= lr * gv;
            }
            clamp(&mut x);
            v = total(&x, &mut evals)?;
        }
        if v < best_v {
            best_v = v;
            best_x = x;
        }
    }
    let mut violations = Vec::with_capacity(specs.len());
    for (i, spec) in specs.iter().enumerate() {
        if i % CANCELLATION_STRIDE == 0 {
            checkpoint(cx)?;
        }
        if skip.contains(&i) {
            violations.push(0.0);
        } else {
            violations.push(violation_contribution(scalar_at(
                problem, spec.node, &best_x,
            )?));
            evals += 1;
        }
    }
    // The published component vector is the authority source for the published
    // total. Recompute it in the same canonical order instead of trusting the
    // optimizer-carried `best_v`, which may be stale if final evidence
    // evaluation evolves independently of the search loop.
    let total_violation = violations.iter().try_fold(0.0, |total, &violation| {
        checked_total_violation(total, violation)
    })?;
    Ok(ElasticReport {
        x: best_x,
        total_violation,
        violations,
        evals,
    })
}

/// One suggested repair.
#[derive(Debug, Clone, PartialEq)]
pub struct RepairAction {
    /// What to do, in words (agent-facing).
    pub description: String,
    /// Structured form.
    pub kind: RepairKind,
    /// Estimated probability the repaired space is feasible
    /// (Monte-Carlo over the domain; calibrated in the battery).
    pub feasibility_estimate: f64,
}

/// Structured repair kinds.
#[derive(Debug, Clone, PartialEq)]
pub enum RepairKind {
    /// Relax constraint `index` by adding `slack` to its bound
    /// (`g ≤ 0` becomes `g ≤ slack`).
    RelaxBound {
        /// Which constraint.
        index: usize,
        /// Added slack.
        slack: f64,
    },
    /// Drop a soft constraint entirely.
    DropSoft {
        /// Which constraint.
        index: usize,
    },
}

/// The full diagnosis payload (the agent-facing artifact).
#[derive(Debug, Clone)]
pub struct Diagnosis {
    /// Whether the space is feasible as stated.
    pub feasible: bool,
    /// A feasible point when one exists.
    pub witness: Option<Vec<f64>>,
    /// MINIMAL unsat core (constraint indices), empty when feasible.
    pub core: Vec<usize>,
    /// Ranked repairs (best first), empty when feasible.
    pub repairs: Vec<RepairAction>,
    /// Elastic-solve evidence.
    pub elastic: ElasticReport,
}

impl Diagnosis {
    fn invalid_reason(&self, specs: &[ConstraintSpec]) -> Option<&'static str> {
        if !self.elastic.total_violation.is_finite() {
            return Some("nonfinite-total-violation");
        }
        if self.elastic.total_violation < 0.0 {
            return Some("negative-total-violation");
        }
        if self.elastic.x.iter().any(|value| !value.is_finite()) {
            return Some("nonfinite-elastic-point");
        }
        if self.elastic.violations.len() != specs.len() {
            return Some("component-violation-count-mismatch");
        }
        if self
            .elastic
            .violations
            .iter()
            .any(|value| !value.is_finite())
        {
            return Some("nonfinite-component-violation");
        }
        if self.elastic.violations.iter().any(|&value| value < 0.0) {
            return Some("negative-component-violation");
        }
        let component_total = self.elastic.violations.iter().sum::<f64>();
        if !component_total.is_finite() {
            return Some("nonfinite-component-violation-total");
        }
        if component_total != self.elastic.total_violation {
            return Some("total-component-violation-mismatch");
        }
        if self.core.iter().any(|&index| index >= specs.len()) {
            return Some("unknown-core-constraint");
        }
        if self.core.windows(2).any(|pair| pair[0] >= pair[1]) {
            return Some("noncanonical-core-order");
        }
        if self
            .witness
            .as_ref()
            .is_some_and(|point| point.iter().any(|value| !value.is_finite()))
        {
            return Some("nonfinite-witness");
        }
        for repair in &self.repairs {
            if !repair.feasibility_estimate.is_finite()
                || !(0.0..=1.0).contains(&repair.feasibility_estimate)
            {
                return Some("invalid-feasibility-estimate");
            }
            match &repair.kind {
                RepairKind::RelaxBound { index, slack } => {
                    if *index >= specs.len() {
                        return Some("unknown-repair-constraint");
                    }
                    if !slack.is_finite() || *slack < 0.0 {
                        return Some("invalid-repair-slack");
                    }
                }
                RepairKind::DropSoft { index } => {
                    if *index >= specs.len() {
                        return Some("unknown-repair-constraint");
                    }
                }
            }
        }
        if self.feasible {
            if self.elastic.total_violation > FEAS_TOL {
                return Some("feasible-claim-has-positive-violation");
            }
            if self.witness.is_none() {
                return Some("feasible-claim-missing-witness");
            }
            if self.witness.as_deref() != Some(self.elastic.x.as_slice()) {
                return Some("witness-does-not-match-elastic-point");
            }
            if self
                .elastic
                .violations
                .iter()
                .any(|&violation| violation > FEAS_TOL)
            {
                return Some("feasible-claim-has-component-violation");
            }
            if !self.core.is_empty() || !self.repairs.is_empty() {
                return Some("feasible-claim-has-conflict-evidence");
            }
        } else {
            if self.elastic.total_violation <= FEAS_TOL {
                return Some("infeasible-claim-lacks-positive-violation");
            }
            if self.witness.is_some() {
                return Some("infeasible-claim-has-witness");
            }
            if self.core.is_empty() {
                return Some("infeasible-claim-missing-core");
            }
        }
        None
    }

    /// Canonical JSON payload for the ledger/session surface. Dynamic text is
    /// escaped. A publicly forged inconsistent or non-finite diagnosis emits a
    /// deterministic invalid/no-claim object; it never retains `feasible:true`
    /// while silently replacing required evidence with `null`.
    #[must_use]
    pub fn to_json(&self, specs: &[ConstraintSpec]) -> String {
        use std::fmt::Write as _;

        if let Some(reason) = self.invalid_reason(specs) {
            let mut invalid = "{\"valid\":false,\"reason\":".to_string();
            push_json_string(&mut invalid, reason);
            invalid.push_str(
                ",\"feasible\":false,\"total_violation\":null,\"core\":[],\"repairs\":[]}",
            );
            return invalid;
        }

        let mut s = format!("{{\"feasible\":{},\"total_violation\":", self.feasible);
        let _ = write!(s, "{:.3e}", self.elastic.total_violation);
        s.push_str(",\"core\":[");
        for (k, &i) in self.core.iter().enumerate() {
            if k > 0 {
                s.push(',');
            }
            push_json_string(&mut s, &specs[i].name);
        }
        s.push_str("],\"repairs\":[");
        for (k, r) in self.repairs.iter().enumerate() {
            if k > 0 {
                s.push(',');
            }
            s.push_str("{\"action\":");
            push_json_string(&mut s, &r.description);
            s.push_str(",\"est_feasible\":");
            let _ = write!(s, "{:.2}", r.feasibility_estimate);
            s.push('}');
        }
        s.push_str("]}");
        s
    }
}

/// Monte-Carlo feasible-volume fraction with constraint `relax[i]`
/// slack applied (the repair feasibility estimator; deterministic).
fn feasible_fraction(
    problem: &Problem,
    specs: &[ConstraintSpec],
    domain: &DomainBox,
    relax: &[(usize, f64)],
    drop: Option<usize>,
    samples: u32,
    cx: &Cx<'_>,
) -> Result<f64, ConError> {
    let mut rng = Lcg(0x1001_2026_0707_0002);
    let mut hits = 0u32;
    for sample in 0..samples {
        if sample % u32::try_from(CANCELLATION_STRIDE).expect("small stride") == 0 {
            checkpoint(cx)?;
        }
        let x: Vec<f64> = domain
            .ranges
            .iter()
            .map(|&(lo, hi)| lo + (hi - lo) * rng.unit())
            .collect();
        let mut ok = true;
        for (i, spec) in specs.iter().enumerate() {
            if i % CANCELLATION_STRIDE == 0 {
                checkpoint(cx)?;
            }
            if Some(i) == drop {
                continue;
            }
            let slack = relax.iter().find(|(j, _)| *j == i).map_or(0.0, |(_, s)| *s);
            // A non-finite constraint value is undefined here, hence NOT feasible
            // — `NaN > slack` is false, which would otherwise count the sample as
            // feasible and inflate the feasibility estimate.
            let gi = scalar_at(problem, spec.node, &x)?;
            if !gi.is_finite() || gi > slack {
                ok = false;
                break;
            }
        }
        if ok {
            hits += 1;
        }
    }
    Ok(f64::from(hits) / f64::from(samples))
}

fn elastic_solve_subset(
    problem: &Problem,
    specs: &[ConstraintSpec],
    domain: &DomainBox,
    members: &[usize],
    cx: &Cx<'_>,
) -> Result<ElasticReport, ConError> {
    let skip: Vec<usize> = (0..specs.len())
        .filter(|index| !members.contains(index))
        .collect();
    elastic_solve(problem, specs, domain, &skip, cx)
}

/// Diagnose a constraint set over a domain: feasibility, MINIMAL unsat
/// core (deletion-filtered), and ranked repairs with feasibility
/// estimates.
///
/// # Errors
/// Evaluation teaching errors; cancellation polls inside the solves.
pub fn diagnose_infeasibility(
    problem: &Problem,
    specs: &[ConstraintSpec],
    domain: &DomainBox,
    cx: &Cx<'_>,
) -> Result<Diagnosis, ConError> {
    checkpoint(cx)?;
    let elastic = elastic_solve(problem, specs, domain, &[], cx)?;
    if elastic.total_violation <= FEAS_TOL {
        return Ok(Diagnosis {
            feasible: true,
            witness: Some(elastic.x.clone()),
            core: Vec::new(),
            repairs: Vec::new(),
            elastic,
        });
    }
    // Candidate core: the elastic support (violated at the optimum). A
    // support identifies the sum-optimum's active trade-off, but need not be
    // jointly infeasible by itself. Verify it before deletion filtering and
    // deterministically expand to the full, already-proven infeasible set when
    // the support is feasible.
    let mut core: Vec<usize> = elastic
        .violations
        .iter()
        .enumerate()
        .filter(|&(_, &v)| v > FEAS_TOL)
        .map(|(i, _)| i)
        .collect();
    let support = elastic_solve_subset(problem, specs, domain, &core, cx)?;
    if support.total_violation <= FEAS_TOL {
        core = (0..specs.len()).collect();
    }

    // Deletion filter for MINIMALITY. The current core is jointly infeasible
    // on entry. A removal is installed only when the resulting subset is also
    // jointly infeasible, so that invariant is preserved at every step.
    let mut k = 0;
    while k < core.len() {
        checkpoint(cx)?;
        let mut without_members = core.clone();
        without_members.remove(k);
        let without = elastic_solve_subset(problem, specs, domain, &without_members, cx)?;
        if without.total_violation <= FEAS_TOL {
            k += 1; // necessary: dropping it restores feasibility
        } else {
            core = without_members; // redundant: still infeasible without it
        }
    }
    let verified_core = elastic_solve_subset(problem, specs, domain, &core, cx)?;
    assert!(
        verified_core.total_violation > FEAS_TOL,
        "deletion filtering must not publish a jointly feasible unsat core"
    );
    // Repairs: relax each core member by graded slacks, or drop it if
    // soft; estimate feasibility by Monte-Carlo volume; rank.
    let mut repairs = Vec::new();
    for &i in &core {
        checkpoint(cx)?;
        let scale = elastic.violations[i].max(FEAS_TOL);
        for factor in [1.1, 1.5] {
            let slack = scale * factor;
            if !slack.is_finite() {
                return Err(ConError::BadParam {
                    what: "repair slack",
                    value: slack,
                });
            }
            let est = feasible_fraction(problem, specs, domain, &[(i, slack)], None, 400, cx)?;
            repairs.push(RepairAction {
                description: format!("relax `{}` by {slack:.3} (g <= {slack:.3})", specs[i].name),
                kind: RepairKind::RelaxBound { index: i, slack },
                feasibility_estimate: est,
            });
        }
        if matches!(specs[i].kind, crate::ConstraintKind::Soft(_)) {
            let est = feasible_fraction(problem, specs, domain, &[], Some(i), 400, cx)?;
            repairs.push(RepairAction {
                description: format!("drop soft constraint `{}`", specs[i].name),
                kind: RepairKind::DropSoft { index: i },
                feasibility_estimate: est,
            });
        }
    }
    repairs.sort_by(|a, b| {
        b.feasibility_estimate
            .partial_cmp(&a.feasibility_estimate)
            .expect("estimates are finite")
            .then_with(|| a.description.cmp(&b.description))
    });
    Ok(Diagnosis {
        feasible: false,
        witness: None,
        core,
        repairs,
        elastic,
    })
}
