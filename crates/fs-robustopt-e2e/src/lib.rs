//! fs-robustopt-e2e — ProofRobust: certified-global-optimal AND robust design.
//! Layer: L4 (ASCENT).
//!
//! # The campaign
//!
//! A local optimizer hands you a point and a shrug — no proof it is the global
//! best, no account of what a manufacturing perturbation does to it. This picks
//! a design with BOTH guarantees, from crates never designed to meet:
//!
//! - **Global optimality, PROVEN** ([`fs_sos`]): each design family's nominal
//!   cost is a convex quadratic `p(x) = a x² + b x + c`. `certify_quadratic`
//!   returns the exact global minimum together with an executable
//!   sum-of-squares certificate `p(x) − p* = (√a·x + b/2√a)²` — a machine-checked
//!   proof, verified by matching polynomial coefficients, that NO `x` does
//!   better. There is no local-optimum ambiguity.
//! - **Worst-case robustness** ([`fs_robust`]): the realized design deviates
//!   from `x*` by a manufacturing tolerance; the perturbed cost is
//!   `p(x*+δ) = p* + a·δ²`, so a STEEP family (large `a`) is punished. CVaR over
//!   a deterministic tolerance grid gives each family's worst-case cost.
//! - **Honest colors** ([`fs_evidence`]): the nominal optimum is `Verified` (SOS
//!   proof); the robust CVaR is `Estimated` (finite sample). The headline claim
//!   never outranks the weakest input.
//!
//! The punchline: the family with the LOWEST nominal cost is not the robust
//! winner — a flatter family with a higher nominal cost wins under CVaR, and the
//! campaign PROVES both nominal optima globally while ranking them robustly.
//! Deterministic; no dependencies beyond the composed crates.

use fs_evidence::{Color, ColorRank};
use fs_robust::{ColoredObjective, cvar};
use fs_sos::{Poly, certify_quadratic};

/// A design family with nominal cost `p(x) = a·x² + b·x + c` (`a > 0`, convex).
#[derive(Debug, Clone)]
pub struct Family {
    /// Name.
    pub name: String,
    /// Quadratic coefficient (curvature; the perturbation-sensitivity).
    pub a: f64,
    /// Linear coefficient.
    pub b: f64,
    /// Constant.
    pub c: f64,
}

impl Family {
    /// A family.
    #[must_use]
    pub fn new(name: impl Into<String>, a: f64, b: f64, c: f64) -> Family {
        Family {
            name: name.into(),
            a,
            b,
            c,
        }
    }

    fn poly(&self) -> Poly {
        Poly::new(vec![self.c, self.b, self.a])
    }

    /// The optimizing design `x* = −b/2a`.
    #[must_use]
    pub fn x_star(&self) -> f64 {
        -self.b / (2.0 * self.a)
    }

    /// Deterministic perturbed costs over a tolerance grid `±sigma` (`n` points):
    /// `p(x* + δ)`.
    #[must_use]
    pub fn perturbed_costs(&self, sigma: f64, n: usize) -> Vec<f64> {
        let xs = self.x_star();
        (0..n)
            .map(|k| {
                let t = if n <= 1 {
                    0.0
                } else {
                    2.0 * (k as f64 / (n - 1) as f64) - 1.0
                };
                let x = xs + t * sigma;
                self.a.mul_add(x * x, self.b.mul_add(x, self.c))
            })
            .collect()
    }
}

/// A per-family verdict.
#[derive(Debug, Clone)]
pub struct FamilyVerdict {
    /// Name.
    pub name: String,
    /// The proven global-minimum cost.
    pub nominal_cost: f64,
    /// The optimizing design.
    pub x_star: f64,
    /// The nominal cost's color (`Verified` iff the SOS proof checks).
    pub nominal_color: Color,
    /// The worst-case (CVaR) cost under manufacturing perturbation.
    pub robust_cost: f64,
}

/// The campaign report.
#[derive(Debug, Clone)]
pub struct RobustOptReport {
    /// Per-family verdicts (input order).
    pub families: Vec<FamilyVerdict>,
    /// The name of the family with the lowest NOMINAL cost.
    pub nominal_winner: String,
    /// The name of the family with the lowest ROBUST (CVaR) cost.
    pub robust_winner: String,
    /// True when robustness reorders the ranking (the interesting case).
    pub robustness_reorders: bool,
    /// Number of families whose nominal optimum is SOS-certified.
    pub certified_count: usize,
    /// The headline color rank (weakest across the robust winner's inputs).
    pub headline_rank: ColorRank,
}

/// Run the ProofRobust campaign: prove every family's global optimum, rank by
/// worst-case CVaR at level `alpha` under a `±sigma` (`n`-point) tolerance grid.
///
/// # Panics
/// If `families` is empty.
#[must_use]
pub fn run_campaign(families: &[Family], alpha: f64, sigma: f64, n: usize) -> RobustOptReport {
    assert!(!families.is_empty(), "need at least one family");
    let mut verdicts = Vec::new();
    let mut objectives = Vec::new();
    let mut certified_count = 0usize;

    for fam in families {
        let poly = fam.poly();
        // Prove the global optimum with a sum-of-squares certificate.
        let (nominal_cost, nominal_color) = match certify_quadratic(fam.a, fam.b, fam.c) {
            Some(cert) => match cert.certified_bound(&poly, 1e-9) {
                Some(bound) => {
                    certified_count += 1;
                    let slack = cert.residual(&poly).max(1e-12);
                    (
                        bound,
                        Color::Verified {
                            lo: bound - slack,
                            hi: bound + slack,
                        },
                    )
                }
                None => (
                    poly.eval(fam.x_star()),
                    Color::Estimated {
                        estimator: "sos-unverified".to_string(),
                        dispersion: 0.0,
                    },
                ),
            },
            None => (
                f64::NEG_INFINITY,
                Color::Estimated {
                    estimator: "unbounded-below".to_string(),
                    dispersion: f64::INFINITY,
                },
            ),
        };
        // Worst-case cost under manufacturing perturbation (CVaR of the grid).
        let costs = fam.perturbed_costs(sigma, n);
        let robust_cost = cvar(&costs, alpha).unwrap_or(f64::INFINITY);
        // A colored objective: the nominal proof is the input evidence, the CVaR
        // is the (Estimated) robust value.
        objectives.push(ColoredObjective::new(
            fam.name.clone(),
            costs,
            vec![nominal_color.clone()],
        ));
        verdicts.push(FamilyVerdict {
            name: fam.name.clone(),
            nominal_cost,
            x_star: fam.x_star(),
            nominal_color,
            robust_cost,
        });
    }

    let nominal_winner = verdicts
        .iter()
        .min_by(|a, b| a.nominal_cost.total_cmp(&b.nominal_cost))
        .map(|v| v.name.clone())
        .unwrap_or_default();
    let robust_idx = verdicts
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.robust_cost.total_cmp(&b.robust_cost))
        .map_or(0, |(i, _)| i);
    let robust_winner = verdicts[robust_idx].name.clone();
    // Headline rank = weakest input color of the robust winner (CVaR is a sample
    // statistic, so the robust claim is Estimated — no laundering).
    let headline_rank = objectives[robust_idx]
        .headline_color()
        .map_or(ColorRank::Estimated, |c| c.rank())
        .min(ColorRank::Estimated);

    RobustOptReport {
        robustness_reorders: nominal_winner != robust_winner,
        nominal_winner,
        robust_winner,
        certified_count,
        headline_rank,
        families: verdicts,
    }
}

/// The worked demonstration: three families with the SAME `x* = 2` but different
/// curvature — the lowest-nominal family is not the robust winner.
#[must_use]
pub fn demo_families() -> Vec<Family> {
    vec![
        Family::new("champion", 1.0, -4.0, 5.2), // nominal min 1.2 (lowest), a=1.0
        Family::new("flat", 0.5, -2.0, 4.0),     // nominal min 2.0, a=0.5 (flattest)
        Family::new("sharp", 2.0, -8.0, 10.0),   // nominal min 2.0, a=2.0 (steepest)
    ]
}
