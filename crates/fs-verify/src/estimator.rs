//! The VERIFIER: equilibrated-flux a-posteriori bounds (Prager–Synge,
//! 1D elliptic class), interval-evaluated to VERIFIED color.
//!
//! The rigor structure: ANY σ with `σ′ = −f` yields the guaranteed
//! bound `‖(u − u_h)′‖ ≤ ‖σ − u_h′‖` — so the free constant in
//! `σ = c − F` is optimized in plain f64 for TIGHTNESS while the bound
//! itself is evaluated with outward-rounded intervals over exact Gauss
//! quadrature (polynomial data ⇒ the quadrature identity is exact;
//! only rounding needs enclosing). An unbounded/NaN enclosure FAILS
//! CLOSED: reject, no color, ever.

use crate::fem1d::{MmsProblem, gauss5, true_energy_error};
use crate::interval::Iv;
use fs_evidence::Color;
use std::fmt::Write as _;

/// Estimator families (Proposal D's independence escalation needs at
/// least two registered per class).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EstimatorFamily {
    /// Equilibrated flux (guaranteed, constant-free — the verifier).
    EquilibratedFlux,
    /// Hierarchical (refined-mesh comparison — independent, NOT
    /// guaranteed; the falsifier's cross-check).
    Hierarchical,
}

impl EstimatorFamily {
    /// Stable id for ledger rows.
    #[must_use]
    pub fn id(self) -> &'static str {
        match self {
            EstimatorFamily::EquilibratedFlux => "equilibrated-flux-1d",
            EstimatorFamily::Hierarchical => "hierarchical-h2",
        }
    }
}

/// The verifier's verdict on one candidate.
#[derive(Debug, Clone)]
pub struct VerifierReport {
    /// The certified error-bound enclosure (energy norm).
    pub bound: Iv,
    /// Accept ⟺ `bound.hi ≤ tolerance` (fail closed on unbounded).
    pub accept: bool,
    /// The verified color carried by an ACCEPT (`None` on reject —
    /// never a badge without the bound).
    pub color: Option<Color>,
    /// The tolerance tested against (feeds the planner).
    pub tolerance: f64,
    /// Estimator family id.
    pub family: &'static str,
    /// FNV hash of the reconstructed flux (ledger identity).
    pub flux_hash: u64,
}

impl VerifierReport {
    /// The review-round-3 ledger row (structured, never stdout).
    #[must_use]
    pub fn to_row(&self, problem: &str, oracle_error: f64) -> String {
        let eff = if oracle_error > 0.0 {
            self.bound.hi / oracle_error
        } else {
            1.0
        };
        let mut s = String::new();
        let _ = write!(
            s,
            "{{\"problem\":\"{problem}\",\"estimator_family_id\":\"{}\",\
             \"flux_hash\":\"{:016X}\",\"bound_lo\":{:.6e},\"bound_hi\":{:.6e},\
             \"oracle_true_error\":{oracle_error:.6e},\"effectivity\":{eff:.4},\
             \"verdict\":\"{}\",\"tolerance\":{:.6e}}}",
            self.family,
            self.flux_hash,
            self.bound.lo,
            self.bound.hi,
            if self.accept { "accept" } else { "reject" },
            self.tolerance
        );
        s
    }
}

fn fnv(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

/// The equilibrated-flux VERIFIER: certify (or reject) a candidate's
/// nodal values against `tolerance`. The returned bound is a TRUE
/// upper bound on `‖(u − u_h)′‖` whenever the candidate satisfies the
/// boundary conditions; the enclosure is rigorous by outward rounding.
#[must_use]
pub fn verify(problem: &MmsProblem, candidate: &[f64], tolerance: f64) -> VerifierReport {
    let m = &problem.mesh;
    let n = m.len();
    // Optimal free constant (f64 — tightness only; ANY c is sound):
    // c* = mean over (0,1) of (F + u_h′).
    let mut mean = 0.0;
    for e in 0..n - 1 {
        let (x0, x1) = (m[e], m[e + 1]);
        let h = x1 - x0;
        let slope = (candidate[e + 1] - candidate[e]) / h;
        for (gx, gw) in gauss5(x0, x1) {
            mean += gw * (problem.big_f.eval(gx) + slope);
        }
    }
    let c_star = mean; // domain length is 1
    // η² = ∫ (c − F − u_h′)², interval-evaluated, exact quadrature.
    let mut eta_sq = Iv::zero();
    for e in 0..n - 1 {
        let (x0, x1) = (m[e], m[e + 1]);
        let h = x1 - x0;
        let slope = (candidate[e + 1] - candidate[e]) / h;
        let mid = f64::midpoint(x0, x1);
        let half = 0.5 * (x1 - x0);
        for (gn, gw) in GAUSS5_REF {
            // Interval node/weight enclosures (nudged constants).
            let node = Iv::point(mid).add(Iv::point(half).mul(iv_c(gn)));
            let weight = iv_c(gw).scale_pos(half.max(1e-300));
            let big_f = problem.big_f.eval_iv(node);
            let r = Iv::point(c_star).sub(big_f).sub(Iv::point(slope));
            eta_sq = eta_sq.add(weight.mul(r.sq()));
        }
    }
    let bound = eta_sq.sqrt();
    // FAIL CLOSED: unbounded/NaN enclosures never accept, never color.
    let sound = !bound.is_unbounded() && !bound.lo.is_nan() && !bound.hi.is_nan();
    let accept = sound && bound.hi <= tolerance;
    let color = if accept {
        Some(Color::Verified {
            lo: 0.0,
            hi: bound.hi,
        })
    } else {
        None
    };
    let mut flux_bytes = Vec::new();
    flux_bytes.extend_from_slice(&c_star.to_bits().to_le_bytes());
    for c in &problem.big_f.0 {
        flux_bytes.extend_from_slice(&c.to_bits().to_le_bytes());
    }
    VerifierReport {
        bound,
        accept,
        color,
        tolerance,
        family: EstimatorFamily::EquilibratedFlux.id(),
        flux_hash: fnv(&flux_bytes),
    }
}

const GAUSS5_REF: [(f64, f64); 5] = [
    (-0.906_179_845_938_664, 0.236_926_885_056_189),
    (-0.538_469_310_105_683, 0.478_628_670_499_366),
    (0.0, 0.568_888_888_888_889),
    (0.538_469_310_105_683, 0.478_628_670_499_366),
    (0.906_179_845_938_664, 0.236_926_885_056_189),
];

/// One-ulp-widened constant (the tabulated Gauss data carries ~1 ulp
/// of transcription error; widening keeps enclosures honest).
fn iv_c(v: f64) -> Iv {
    Iv {
        lo: crate::interval::down(v),
        hi: crate::interval::up(v),
    }
}

/// The INDEPENDENT second family: hierarchical estimate from a
/// uniformly refined solve (`h/2`). Not guaranteed — the falsifier's
/// cross-check, never a color source.
#[must_use]
pub fn hierarchical_estimate(problem: &MmsProblem, candidate: &[f64]) -> f64 {
    let mut fine_mesh = Vec::with_capacity(problem.mesh.len() * 2 - 1);
    for w in problem.mesh.windows(2) {
        fine_mesh.push(w[0]);
        fine_mesh.push(f64::midpoint(w[0], w[1]));
    }
    fine_mesh.push(*problem.mesh.last().expect("nonempty mesh"));
    let fine = MmsProblem::new(&problem.name, problem.u.clone(), fine_mesh);
    let fine_u = crate::fem1d::solve_p1(&fine);
    // ‖u_{h/2}′ − u_h′‖ over the fine mesh.
    let mut acc = 0.0;
    for e in 0..fine.mesh.len() - 1 {
        let (x0, x1) = (fine.mesh[e], fine.mesh[e + 1]);
        let h = x1 - x0;
        let fine_slope = (fine_u[e + 1] - fine_u[e]) / h;
        // The coarse element containing this fine element.
        let coarse_e = e / 2;
        let ch = problem.mesh[coarse_e + 1] - problem.mesh[coarse_e];
        let coarse_slope = (candidate[coarse_e + 1] - candidate[coarse_e]) / ch;
        let d = fine_slope - coarse_slope;
        acc += h * d * d;
    }
    acc.sqrt()
}

/// The nonlinear WARM-START fallback: the candidate is accepted only
/// as a starting point; the measured value is iteration savings and
/// the color is ESTIMATED, never verified (the honest R1 boundary).
#[derive(Debug, Clone)]
pub struct WarmStartReport {
    /// Newton iterations from a cold start (zero).
    pub cold_iterations: u32,
    /// Newton iterations from the candidate.
    pub warm_iterations: u32,
    /// The color of the claim (always `Estimated`).
    pub color: Color,
}

/// Measure warm-start savings on the toy nonlinear class.
#[must_use]
pub fn warm_start(problem: &MmsProblem, candidate: &[f64], max_iter: u32) -> WarmStartReport {
    let zero = vec![0.0; problem.mesh.len()];
    let (_, cold) = crate::fem1d::solve_nonlinear(problem, &zero, max_iter);
    let (_, warm) = crate::fem1d::solve_nonlinear(problem, candidate, max_iter);
    WarmStartReport {
        cold_iterations: cold,
        warm_iterations: warm,
        color: Color::Estimated {
            estimator: "warm-start-iteration-savings".to_string(),
            dispersion: f64::INFINITY,
        },
    }
}

/// Convenience for the batteries: effectivity of a report against the
/// oracle.
#[must_use]
pub fn effectivity(problem: &MmsProblem, candidate: &[f64], report: &VerifierReport) -> f64 {
    let truth = true_energy_error(problem, candidate);
    if truth > 0.0 {
        report.bound.hi / truth
    } else {
        1.0
    }
}
