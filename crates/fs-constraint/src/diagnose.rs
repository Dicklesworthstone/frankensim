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
//! DELETION FILTER, so the result is MINIMAL: dropping ANY member
//! restores feasibility — a property the conformance battery checks
//! against brute-force enumeration.

use crate::{ConError, ConstraintSpec, scalar_at};
use fs_exec::Cx;
use fs_opt::Problem;

/// Per-component design-domain box.
#[derive(Debug, Clone, PartialEq)]
pub struct DomainBox {
    /// (lo, hi) per component of the (single, Rn) design variable.
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

/// Minimize `Σ max(gᵢ(x), 0)` over the box: multi-start projected
/// subgradient descent (deterministic). Small-fixture machinery — the
/// production restoration solver is a later ASCENT bead.
///
/// # Errors
/// Evaluation teaching errors carried through; cancellation polls.
pub fn elastic_solve(
    problem: &Problem,
    specs: &[ConstraintSpec],
    domain: &DomainBox,
    skip: &[usize],
    cx: &Cx<'_>,
) -> Result<ElasticReport, ConError> {
    let n = domain.ranges.len();
    let active: Vec<usize> = (0..specs.len()).filter(|i| !skip.contains(i)).collect();
    let mut evals = 0u64;
    let mut total = |x: &[f64], evals: &mut u64| -> Result<f64, ConError> {
        let mut t = 0.0;
        for &i in &active {
            t += scalar_at(problem, specs[i].node, x)?.max(0.0);
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
        if cx.checkpoint().is_err() {
            return Err(ConError::Eval(fs_opt::OptError::Cancelled));
        }
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
            if v <= FEAS_TOL {
                break;
            }
            // FD subgradient of the hinge sum.
            let h = 1e-6 * diam.max(1.0);
            let mut g = vec![0.0; n];
            for (k, gk) in g.iter_mut().enumerate() {
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
        if skip.contains(&i) {
            violations.push(0.0);
        } else {
            violations.push(scalar_at(problem, spec.node, &best_x)?.max(0.0));
        }
    }
    Ok(ElasticReport {
        x: best_x,
        total_violation: best_v,
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
    /// Canonical JSON payload for the ledger/session surface.
    #[must_use]
    pub fn to_json(&self, specs: &[ConstraintSpec]) -> String {
        use std::fmt::Write as _;
        let mut s = format!(
            "{{\"feasible\":{},\"total_violation\":{:.3e},\"core\":[",
            self.feasible, self.elastic.total_violation
        );
        for (k, &i) in self.core.iter().enumerate() {
            if k > 0 {
                s.push(',');
            }
            let _ = write!(s, "\"{}\"", specs[i].name);
        }
        s.push_str("],\"repairs\":[");
        for (k, r) in self.repairs.iter().enumerate() {
            if k > 0 {
                s.push(',');
            }
            let _ = write!(
                s,
                "{{\"action\":\"{}\",\"est_feasible\":{:.2}}}",
                r.description, r.feasibility_estimate
            );
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
) -> Result<f64, ConError> {
    let mut rng = Lcg(0x1001_2026_0707_0002);
    let mut hits = 0u32;
    for _ in 0..samples {
        let x: Vec<f64> = domain
            .ranges
            .iter()
            .map(|&(lo, hi)| lo + (hi - lo) * rng.unit())
            .collect();
        let mut ok = true;
        for (i, spec) in specs.iter().enumerate() {
            if Some(i) == drop {
                continue;
            }
            let slack = relax
                .iter()
                .find(|(j, _)| *j == i)
                .map_or(0.0, |(_, s)| *s);
            if scalar_at(problem, spec.node, &x)? > slack {
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
    // Candidate core: the elastic support (violated at the optimum).
    let mut core: Vec<usize> = elastic
        .violations
        .iter()
        .enumerate()
        .filter(|&(_, &v)| v > FEAS_TOL)
        .map(|(i, _)| i)
        .collect();
    // Deletion filter for MINIMALITY: a member stays only if the rest
    // of the core remains infeasible WITH that member removed (against
    // the full non-core context, which is feasible by construction of
    // the support — we test against all constraints minus the member).
    let mut k = 0;
    while k < core.len() {
        let member = core[k];
        let mut skip: Vec<usize> = (0..specs.len())
            .filter(|i| !core.contains(i))
            .collect();
        skip.push(member);
        let without = elastic_solve(problem, specs, domain, &skip, cx)?;
        if without.total_violation <= FEAS_TOL {
            k += 1; // necessary: dropping it restores feasibility
        } else {
            core.remove(k); // redundant: still infeasible without it
        }
    }
    // Repairs: relax each core member by graded slacks, or drop it if
    // soft; estimate feasibility by Monte-Carlo volume; rank.
    let mut repairs = Vec::new();
    for &i in &core {
        let scale = elastic.violations[i].max(FEAS_TOL);
        for factor in [1.1, 1.5] {
            let slack = scale * factor;
            let est = feasible_fraction(problem, specs, domain, &[(i, slack)], None, 400)?;
            repairs.push(RepairAction {
                description: format!(
                    "relax `{}` by {slack:.3} (g <= {slack:.3})",
                    specs[i].name
                ),
                kind: RepairKind::RelaxBound { index: i, slack },
                feasibility_estimate: est,
            });
        }
        if matches!(specs[i].kind, crate::ConstraintKind::Soft(_)) {
            let est = feasible_fraction(problem, specs, domain, &[], Some(i), 400)?;
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
