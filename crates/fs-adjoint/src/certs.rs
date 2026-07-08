//! GRADIENT CERTIFICATES (addendum Proposal 1, bead bk0o.3; [S], behind
//! the `gradient-certs` feature): gradients are CLAIMS like any other
//! and get COLORS like any other — a gradient without a certificate is
//! folklore. Three pieces:
//!
//! - an INTERVAL BOUND on the adjoint consistency residual
//!   (`⟨Av, w⟩ − ⟨v, Aᵀw⟩` evaluated in outward-rounded fs-ivl
//!   arithmetic — a verified-color enclosure where the path is
//!   differentiable);
//! - MANDATORY finite-difference spot checks along seeded random
//!   directions via the falsifier registry's `adjoint-gradient` →
//!   `finite-difference-spot-check` pairing (Proposal 6): the
//!   independent cross-check a transpose bug or sign error trips;
//! - the MERGE GATE: no gradient merges without a passing check —
//!   the base plan's CI gradient-gate discipline extended to
//!   seam-crossing, ledger-transposed gradients.

use crate::mitigate::GradientGrade;
use crate::transpose::{FdVerdict, fd_falsifier};
use fs_evidence::{Color, ValidityDomain};
use fs_ivl::Interval;

/// A sparse linear operator (rows of `(col, weight)`) with point and
/// INTERVAL applies — the certifiable form of restriction/conversion
/// seams.
#[derive(Debug, Clone)]
pub struct SparseLinear {
    /// Rows of (column, weight).
    pub rows: Vec<Vec<(usize, f64)>>,
    /// Column count.
    pub ncols: usize,
}

impl SparseLinear {
    fn apply_iv(&self, x: &[Interval]) -> Vec<Interval> {
        self.rows
            .iter()
            .map(|row| {
                row.iter().fold(Interval::point(0.0), |acc, &(c, w)| {
                    acc + Interval::point(w) * x[c]
                })
            })
            .collect()
    }

    fn apply_t_iv(&self, y: &[Interval]) -> Vec<Interval> {
        let mut out = vec![Interval::point(0.0); self.ncols];
        for (r, row) in self.rows.iter().enumerate() {
            for &(c, w) in row {
                out[c] = out[c] + Interval::point(w) * y[r];
            }
        }
        out
    }
}

fn iv_dot(a: &[Interval], b: &[Interval]) -> Interval {
    a.iter()
        .zip(b)
        .fold(Interval::point(0.0), |acc, (x, y)| acc + *x * *y)
}

/// A VERIFIED enclosure of the worst adjoint consistency residual
/// `|⟨Av, w⟩ − ⟨v, Aᵀw⟩|` over seeded probes, computed entirely in
/// outward-rounded interval arithmetic — the true fp residual of the
/// registered transpose pair lies INSIDE the returned bound.
#[must_use]
pub fn adjoint_residual_bound(op: &SparseLinear, probes: usize) -> f64 {
    let n_out = op.rows.len();
    let mut state = 0x0dd5_eed5_u64;
    let mut lcg = move || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 11) as f64) / (1u64 << 53) as f64 - 0.5
    };
    let mut worst = 0.0f64;
    for _ in 0..probes {
        let v: Vec<Interval> = (0..op.ncols).map(|_| Interval::point(lcg())).collect();
        let w: Vec<Interval> = (0..n_out).map(|_| Interval::point(lcg())).collect();
        let lhs = iv_dot(&op.apply_iv(&v), &w);
        let rhs = iv_dot(&v, &op.apply_t_iv(&w));
        let residual = (lhs - rhs).abs();
        // The enclosure's upper end bounds the residual soundly.
        worst = worst.max(residual.midpoint() + residual.width() / 2.0);
    }
    worst
}

/// One gradient's certificate: color + evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct GradientCertificate {
    /// The Proposal-3 color.
    pub color: Color,
    /// The verified adjoint-consistency residual bound, where computed.
    pub residual_bound: Option<f64>,
    /// The FD spot-check verdicts (the falsifier's evidence).
    pub fd_checks: Vec<FdVerdict>,
    /// The discontinuity flag inherited from the routing grade.
    pub discontinuity: Option<Vec<String>>,
}

impl GradientCertificate {
    /// True when every FD spot check agreed.
    #[must_use]
    pub fn fd_all_consistent(&self) -> bool {
        !self.fd_checks.is_empty() && self.fd_checks.iter().all(|v| v.consistent)
    }
}

/// The mandatory FD spot checks: seeded random directions through the
/// falsifier pairing (`adjoint-gradient` → `finite-difference-spot-
/// check`). Deterministic directions; conditioning-aware tolerances.
#[must_use]
pub fn fd_spot_checks(
    f: &dyn Fn(&[f64]) -> f64,
    x: &[f64],
    grad: &[f64],
    directions: usize,
    seed: u64,
) -> Vec<FdVerdict> {
    let mut state = seed;
    let mut lcg = move || {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((state >> 11) as f64) / (1u64 << 53) as f64 - 0.5
    };
    (0..directions)
        .map(|_| {
            let dir: Vec<f64> = (0..x.len()).map(|_| lcg()).collect();
            let dd: f64 = grad.iter().zip(&dir).map(|(g, d)| g * d).sum();
            fd_falsifier(f, x, &dir, dd, 1e-5, 1e-7)
        })
        .collect()
}

/// Anchoring evidence for validated-color gradients (Proposal 11:
/// assimilated experimental data).
#[derive(Debug, Clone)]
pub struct Anchor {
    /// The anchoring dataset's identity.
    pub dataset: String,
    /// The regime in which the anchoring holds.
    pub regime: ValidityDomain,
}

/// Assign the gradient's color from its path grade, residual evidence,
/// and (optional) experimental anchoring:
///
/// - flagged remesh / surrogate path → the grade's ESTIMATED color
///   (never upgraded here — that would be laundering);
/// - smooth + anchored → VALIDATED in the anchor's regime;
/// - smooth + interval-bounded residual → VERIFIED with the bound;
/// - smooth with NO residual evidence → ESTIMATED (a gradient without
///   a certificate is folklore).
#[must_use]
pub fn certify(
    grade: &GradientGrade,
    residual_bound: Option<f64>,
    fd_checks: Vec<FdVerdict>,
    anchor: Option<&Anchor>,
) -> GradientCertificate {
    let (color, discontinuity) = match grade {
        GradientGrade::EstimatedWithDiscontinuity {
            color, crossing, ..
        } => (color.clone(), Some(crossing.clone())),
        GradientGrade::Smooth { .. } => {
            if let Some(a) = anchor {
                (
                    Color::Validated {
                        regime: a.regime.clone(),
                        dataset: a.dataset.clone(),
                    },
                    None,
                )
            } else if let Some(bound) = residual_bound {
                (Color::Verified { lo: 0.0, hi: bound }, None)
            } else {
                (
                    Color::Estimated {
                        estimator: "gradient without a residual certificate".to_string(),
                        dispersion: f64::INFINITY,
                    },
                    None,
                )
            }
        }
    };
    GradientCertificate {
        color,
        residual_bound,
        fd_checks,
        discontinuity,
    }
}

/// A merge-gate refusal: teaches what is missing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GateRefusal {
    /// What failed.
    pub what: String,
}

impl std::fmt::Display for GateRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "gradient merge gate: {} — a solver (or seam) without a passing gradient \
             check cannot merge",
            self.what
        )
    }
}

/// THE CI GRADIENT GATE, extended to seam-crossing gradients: no merge
/// without a passing check.
///
/// # Errors
/// [`GateRefusal`] when FD checks are missing or failing.
pub fn merge_gate(cert: &GradientCertificate) -> Result<(), GateRefusal> {
    if cert.fd_checks.is_empty() {
        return Err(GateRefusal {
            what: "no FD spot checks were run (the falsifier pairing is mandatory)".to_string(),
        });
    }
    if let Some(bad) = cert.fd_checks.iter().find(|v| !v.consistent) {
        return Err(GateRefusal {
            what: format!(
                "an FD spot check failed (adjoint dd {} vs FD {} beyond tolerance {}): a \
                 transpose or sign bug is the likely cause",
                bad.adjoint_dd, bad.fd_fine, bad.tolerance
            ),
        });
    }
    Ok(())
}
