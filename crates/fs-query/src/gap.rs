//! SDF-pair local gap oracle (bead rjnd, E1 query upgrades, part 4)
//! plus the strictly-scoped pointwise overlap-inradius witness
//! (part 7's `max(φ_A, φ_B)` quantity).
//!
//! [`ImplicitGapOracle`] pairs two exact-distance charts and answers
//! pointwise contact-adjacent queries with certificates where they
//! exist and honest estimate/no-claim labels where they do not:
//!
//! - `sum ∈ [lo, hi]` encloses `φ_A(p) + φ_B(p)`. When the point is
//!   certified outside BOTH bodies, the triangle inequality makes
//!   `hi` a rigorous UPPER bound on the distance between the two
//!   bodies (`dist(A,B) ≤ dist(p,A) + dist(p,B)`).
//! - `overlap_inradius`: when `max(φ_A, φ_B)` is certified negative,
//!   the ball of that radius around `p` lies inside BOTH bodies — a
//!   pointwise overlap witness with a certified inradius, and nothing
//!   more (it is NOT a penetration depth and never upgrades to one).
//! - `normal`: the Estimate-class contact-axis direction
//!   `∇φ_A - ∇φ_B`, normalized; absent whenever either chart honestly
//!   declines a gradient or the difference is degenerate. No
//!   certificate accompanies it.

use crate::QueryError;
use fs_evidence::{NumericalCertificate, NumericalKind};
use fs_exec::Cx;
use fs_geom::{Chart, Point3, TraceStepClaim, Vec3};

/// One pointwise answer from [`ImplicitGapOracle::gap_at`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GapSample {
    /// Outward-rounded enclosure of `φ_A(p) + φ_B(p)`.
    pub sum_lo: f64,
    /// Outward-rounded enclosure of `φ_A(p) + φ_B(p)`.
    pub sum_hi: f64,
    /// Rigorous upper bound on `dist(A, B)`, present exactly when the
    /// point is certified outside both bodies.
    pub separation_upper: Option<f64>,
    /// Certified radius of a ball around `p` contained in BOTH bodies,
    /// present exactly when `max(φ_A, φ_B)` is certified negative.
    pub overlap_inradius: Option<f64>,
    /// Estimate-class contact-axis direction (unit `∇φ_A - ∇φ_B`).
    /// Carries no certificate; `None` when either gradient is honestly
    /// absent or the difference is degenerate/non-finite.
    pub normal: Option<[f64; 3]>,
}

/// A pair of exact-distance charts answering local gap queries.
pub struct ImplicitGapOracle<'a> {
    a: &'a dyn Chart,
    b: &'a dyn Chart,
}

impl<'a> ImplicitGapOracle<'a> {
    /// Pair two charts. The certified quantities need the exact
    /// Euclidean signed-distance theorem from BOTH inputs, so weaker
    /// trace claims refuse here, at construction.
    ///
    /// # Errors
    /// [`QueryError::SeparationRequiresExactDistance`] naming the
    /// offending input (`"a"` or `"b"`) and its weaker claim.
    pub fn new(a: &'a dyn Chart, b: &'a dyn Chart) -> Result<ImplicitGapOracle<'a>, QueryError> {
        for (label, chart) in [("a", a), ("b", b)] {
            let claim = chart.trace_step_claim();
            if claim != TraceStepClaim::ExactDistance {
                return Err(QueryError::SeparationRequiresExactDistance {
                    input: if label == "a" { "a" } else { "b" },
                    claim,
                });
            }
        }
        Ok(ImplicitGapOracle { a, b })
    }

    /// Evaluate the local gap quantities at `p`.
    ///
    /// # Errors
    /// [`QueryError::InvalidPointSample`] for a non-finite query
    /// point; [`QueryError::InvalidTraceSample`] when either chart's
    /// enclosure is missing, non-finite, inverted, or weaker than
    /// Exact/Enclosure class; [`QueryError::Cancelled`] on
    /// cancellation (checked after each chart call).
    pub fn gap_at(&self, p: Point3, cx: &Cx<'_>) -> Result<GapSample, QueryError> {
        if !(p.x.is_finite() && p.y.is_finite() && p.z.is_finite()) {
            return Err(QueryError::InvalidPointSample {
                at: [p.x, p.y, p.z],
            });
        }
        let sample_a = self.a.eval(p, cx);
        if cx.checkpoint().is_err() {
            return Err(QueryError::Cancelled);
        }
        let sample_b = self.b.eval(p, cx);
        if cx.checkpoint().is_err() {
            return Err(QueryError::Cancelled);
        }
        let enc_a = self.a.trace_value_enclosure(p, &sample_a, cx);
        let enc_b = self.b.trace_value_enclosure(p, &sample_b, cx);
        validate_gap_enclosure(&enc_a, p)?;
        validate_gap_enclosure(&enc_b, p)?;
        let sum_lo = (enc_a.lo + enc_b.lo).next_down();
        let sum_hi = (enc_a.hi + enc_b.hi).next_up();
        let separation_upper = if enc_a.lo > 0.0 && enc_b.lo > 0.0 {
            Some(sum_hi)
        } else {
            None
        };
        let max_hi = enc_a.hi.max(enc_b.hi);
        let overlap_inradius = if max_hi < 0.0 {
            Some((-max_hi).next_down())
        } else {
            None
        };
        let normal = match (sample_a.gradient, sample_b.gradient) {
            (Some(ga), Some(gb)) => {
                let d = Vec3::new(ga.x - gb.x, ga.y - gb.y, ga.z - gb.z);
                let n = d.norm();
                if n.is_finite() && n > 1e-12 {
                    let u = d.scale(1.0 / n);
                    if u.x.is_finite() && u.y.is_finite() && u.z.is_finite() {
                        Some([u.x, u.y, u.z])
                    } else {
                        None
                    }
                } else {
                    None
                }
            }
            _ => None,
        };
        Ok(GapSample {
            sum_lo,
            sum_hi,
            separation_upper,
            overlap_inradius,
            normal,
        })
    }
}

fn validate_gap_enclosure(enclosure: &NumericalCertificate, p: Point3) -> Result<(), QueryError> {
    let sound = matches!(
        enclosure.kind,
        NumericalKind::Exact | NumericalKind::Enclosure
    ) && enclosure.lo.is_finite()
        && enclosure.hi.is_finite()
        && enclosure.lo <= enclosure.hi;
    if sound {
        Ok(())
    } else {
        Err(QueryError::InvalidTraceSample {
            at: [p.x, p.y, p.z],
        })
    }
}
