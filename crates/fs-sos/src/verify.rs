//! Rigorous verification of SOS solutions — the step that turns a numerical
//! SDP answer into a THEOREM, independent of how the answer was found.
//!
//! ## What is proved
//!
//! For every identity `constant + Σ_k g_k·d_k ≡ 0` of an [`SosProgram`] the
//! verifier designates one SOS term `c·σ_s` with a nonzero constant point
//! multiplier `c` as the SLACK. Every other decision is fixed at its exact
//! floating-point value: scalars and free polynomials as given, SOS Gram
//! matrices as the symmetric matrices defined by their upper triangles. The
//! residual `R = constant + Σ_{k≠s} g_k·d_k` is enclosed coefficient-wise with
//! `fs-ivl` outward-rounded intervals, so `σ_s = −R/c` holds for an exact real
//! polynomial whose coefficients lie in the computed intervals.
//!
//! The slack's Gram matrix is then rebuilt as an INTERVAL matrix `[G]`: for
//! each monomial `α` one designated basis pair `(i, j)` with `mᵢmⱼ = x^α`
//! absorbs the exact requirement (all other pairs keep their float values).
//! Hence some real symmetric `G* ∈ [G]` reproduces `σ_s` exactly.
//!
//! Positive definiteness is proved with the interval Cholesky method: if the
//! algorithm runs to completion on `[G]` with every pivot's lower bound
//! strictly positive, then EVERY symmetric matrix in `[G]` is positive
//! definite (Alefeld & Mayer, "The Cholesky method for interval data",
//! Linear Algebra Appl. 194 (1993); the point Cholesky of any member stays
//! inside the interval iterates by inclusion isotonicity). Non-slack Gram
//! matrices are checked the same way as point-interval matrices.
//!
//! Therefore every SOS decision is a genuine sum of squares and every
//! identity holds exactly — the certificate is sound up to the correctness of
//! IEEE-754 basic operations and `fs-ivl`'s outward rounding (basic ops only:
//! `+ − × ÷ √`).
//!
//! A failure is a refusal with a reason, never a weaker claim.

use std::collections::BTreeMap;

use fs_ivl::Interval;

use crate::mpoly::{IPoly, Monomial, isqr};
use crate::program::{DecisionId, DecisionKind, DecisionValue, SosProgram};

/// Why a solution could not be certified.
#[derive(Debug, Clone, PartialEq)]
pub enum VerifyError {
    /// The value vector does not match the program's decisions.
    ShapeMismatch {
        /// Description.
        what: String,
    },
    /// No SOS term with a constant nonzero point multiplier exists in an
    /// identity (or each candidate is already used as another identity's
    /// slack / appears in several identities).
    NoSlack {
        /// Identity index.
        identity: usize,
    },
    /// A residual monomial cannot be represented by the slack's basis.
    Unrepresentable {
        /// Identity index.
        identity: usize,
        /// The monomial's exponents.
        monomial: Vec<u32>,
    },
    /// The interval Cholesky test failed (the matrix may be indefinite or
    /// merely too close to the PSD boundary to prove definiteness).
    NotPositiveDefinite {
        /// The decision whose Gram matrix failed.
        decision: usize,
        /// The pivot where it failed.
        pivot: usize,
    },
    /// A value is NaN or infinite.
    NonFinite {
        /// The decision.
        decision: usize,
    },
}

impl core::fmt::Display for VerifyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            VerifyError::ShapeMismatch { what } => write!(f, "solution shape mismatch: {what}"),
            VerifyError::NoSlack { identity } => write!(
                f,
                "identity {identity}: no SOS term with a constant nonzero multiplier \
                 (used in no other identity) can absorb the residual"
            ),
            VerifyError::Unrepresentable { identity, monomial } => write!(
                f,
                "identity {identity}: residual monomial {monomial:?} is outside the slack basis"
            ),
            VerifyError::NotPositiveDefinite { decision, pivot } => write!(
                f,
                "decision {decision}: interval Cholesky failed at pivot {pivot} \
                 (indefinite, or too close to the PSD boundary — re-solve centred)"
            ),
            VerifyError::NonFinite { decision } => {
                write!(f, "decision {decision}: non-finite value")
            }
        }
    }
}

impl std::error::Error for VerifyError {}

/// Proof that every identity of a program holds with every SOS decision a
/// genuine sum of squares.
#[derive(Debug, Clone, PartialEq)]
pub struct Certificate {
    /// For each identity, the SOS decision used as slack.
    pub slacks: Vec<DecisionId>,
    /// For each SOS decision (by index), the smallest certified lower bound
    /// of its Cholesky pivots squared (a positive-definiteness margin proxy).
    pub pivot_margins: BTreeMap<usize, f64>,
    /// Largest radius among the absorbed slack Gram entries (how much
    /// rounding the proof had to carry).
    pub max_absorbed_radius: f64,
}

/// Interval Cholesky test on a symmetric interval matrix (row-major; only the
/// lower triangle `i ≥ j` is read). Returns the smallest pivot lower bound
/// (`> 0`) on success, or the failing pivot index.
///
/// # Errors
/// The index of the first pivot whose lower bound is not strictly positive.
pub fn interval_cholesky_pd(g: &[Interval], n: usize) -> Result<f64, usize> {
    let mut l: Vec<Interval> = vec![Interval::point(0.0); n * n];
    let mut min_pivot = f64::INFINITY;
    for j in 0..n {
        let mut s = g[j * n + j];
        for k in 0..j {
            s = s - isqr(l[j * n + k]);
        }
        if !(s.lo() > 0.0) || !s.hi().is_finite() {
            return Err(j);
        }
        min_pivot = min_pivot.min(s.lo());
        let d = s.sqrt();
        l[j * n + j] = d;
        for i in (j + 1)..n {
            let mut t = g[i * n + j];
            for k in 0..j {
                t = t - l[i * n + k] * l[j * n + k];
            }
            l[i * n + j] = t / d;
        }
    }
    Ok(min_pivot)
}

/// The SOS polynomial `Σ_{i≤j} w_ij Q_ij mᵢmⱼ` (w = 1 on the diagonal, 2 off
/// it) enclosed exactly from the upper triangle of `gram`.
fn gram_poly(nvars: usize, basis: &[Monomial], gram: &[f64]) -> IPoly {
    let n = basis.len();
    let mut out = IPoly::zero(nvars);
    for i in 0..n {
        for j in i..n {
            let q = gram[i * n + j];
            let w = if i == j { 1.0 } else { 2.0 };
            out.add_term(
                basis[i].mul(&basis[j]),
                Interval::point(q) * Interval::point(w),
            );
        }
    }
    out
}

/// A point-valued symmetric interval matrix from the upper triangle.
fn point_matrix(gram: &[f64], n: usize) -> Vec<Interval> {
    let mut out = vec![Interval::point(0.0); n * n];
    for i in 0..n {
        for j in i..n {
            let v = Interval::point(gram[i * n + j]);
            out[i * n + j] = v;
            out[j * n + i] = v;
        }
    }
    out
}

/// The exact value of a decision as an interval polynomial.
fn value_poly(nvars: usize, v: &DecisionValue) -> IPoly {
    match v {
        DecisionValue::Scalar(s) => {
            let mut p = IPoly::zero(nvars);
            p.add_term(Monomial::one(nvars), Interval::point(*s));
            p
        }
        DecisionValue::Free(p) => p.to_interval(),
        DecisionValue::Sos { basis, gram } => gram_poly(nvars, basis, gram),
    }
}

fn constant_point(p: &IPoly) -> Option<f64> {
    let mut it = p.terms();
    let (m, c) = it.next()?;
    if it.next().is_some()
        || m.degree() != 0
        || c.lo().to_bits() != c.hi().to_bits()
        || c.lo() == 0.0
    {
        return None;
    }
    Some(c.lo())
}

/// Verify a solution of `program` (see the module docs for the theorem).
///
/// # Errors
/// [`VerifyError`] naming why no proof was found. A refusal is never a claim
/// that the program is infeasible.
#[allow(clippy::too_many_lines)]
pub fn verify(program: &SosProgram, values: &[DecisionValue]) -> Result<Certificate, VerifyError> {
    let nvars = program.nvars();
    let decisions = program.decisions();
    if values.len() != decisions.len() {
        return Err(VerifyError::ShapeMismatch {
            what: format!("{} values for {} decisions", values.len(), decisions.len()),
        });
    }
    for (k, (d, v)) in decisions.iter().zip(values).enumerate() {
        let ok = match (&d.kind, v) {
            (DecisionKind::Scalar, DecisionValue::Scalar(s)) => {
                if !s.is_finite() {
                    return Err(VerifyError::NonFinite { decision: k });
                }
                true
            }
            (DecisionKind::Free(b), DecisionValue::Free(p)) => {
                if p.terms().any(|(_, c)| !c.is_finite()) {
                    return Err(VerifyError::NonFinite { decision: k });
                }
                p.terms().all(|(m, _)| b.contains(m))
            }
            (DecisionKind::Sos(b), DecisionValue::Sos { basis, gram }) => {
                if gram.iter().any(|c| !c.is_finite()) {
                    return Err(VerifyError::NonFinite { decision: k });
                }
                b == basis && gram.len() == b.len() * b.len()
            }
            _ => false,
        };
        if !ok {
            return Err(VerifyError::ShapeMismatch {
                what: format!("decision {k} value does not match its kind/basis"),
            });
        }
    }
    // How many identities each decision appears in.
    let mut appearances = vec![0usize; decisions.len()];
    for ident in program.identities() {
        let mut seen = std::collections::BTreeSet::new();
        for (_, d) in ident.terms() {
            if seen.insert(d.0) {
                appearances[d.0] += 1;
            }
        }
    }
    // Choose slacks: an SOS decision appearing exactly once in exactly this
    // identity, with a constant nonzero point multiplier; prefer the largest
    // basis (most absorbing pairs), then the lowest index (deterministic).
    let mut slacks = Vec::with_capacity(program.identities().len());
    let mut is_slack = vec![false; decisions.len()];
    for (ki, ident) in program.identities().iter().enumerate() {
        let mut best: Option<(usize, usize)> = None; // (term index, basis length)
        for (ti, (mult, d)) in ident.terms().iter().enumerate() {
            let DecisionKind::Sos(b) = &decisions[d.0].kind else {
                continue;
            };
            if appearances[d.0] != 1 || is_slack[d.0] {
                continue;
            }
            if ident.terms().iter().filter(|(_, e)| e == d).count() != 1 {
                continue;
            }
            if constant_point(mult).is_none() {
                continue;
            }
            if best.is_none_or(|(_, bl)| b.len() > bl) {
                best = Some((ti, b.len()));
            }
        }
        let Some((ti, _)) = best else {
            return Err(VerifyError::NoSlack { identity: ki });
        };
        let d = ident.terms()[ti].1;
        is_slack[d.0] = true;
        slacks.push(d);
    }
    let mut pivot_margins = BTreeMap::new();
    // Non-slack SOS Gram matrices: point interval Cholesky.
    for (k, (d, v)) in decisions.iter().zip(values).enumerate() {
        if is_slack[k] {
            continue;
        }
        if let (DecisionKind::Sos(b), DecisionValue::Sos { gram, .. }) = (&d.kind, v) {
            let n = b.len();
            let m = point_matrix(gram, n);
            match interval_cholesky_pd(&m, n) {
                Ok(p) => {
                    pivot_margins.insert(k, p);
                }
                Err(pivot) => return Err(VerifyError::NotPositiveDefinite { decision: k, pivot }),
            }
        }
    }
    // Each identity: enclose the residual, absorb it into the slack Gram.
    let mut max_absorbed_radius = 0.0f64;
    for (ki, ident) in program.identities().iter().enumerate() {
        let slack = slacks[ki];
        let mut r = ident.constant().clone();
        let mut c_slack = 0.0;
        for (mult, d) in ident.terms() {
            if *d == slack {
                c_slack = constant_point(mult).unwrap_or_else(|| unreachable!("slack chosen"));
                continue;
            }
            r = r.add(&mult.mul(&value_poly(nvars, &values[d.0])));
        }
        // target = −R / c, the exact coefficients σ_s must have.
        let inv = Interval::point(-1.0) / Interval::point(c_slack);
        let target = r.scale(inv);
        let DecisionValue::Sos { basis, gram } = &values[slack.0] else {
            unreachable!("slack is SOS")
        };
        let n = basis.len();
        // Pairs per monomial.
        let mut pairs: BTreeMap<Monomial, Vec<(usize, usize)>> = BTreeMap::new();
        for i in 0..n {
            for j in i..n {
                pairs
                    .entry(basis[i].mul(&basis[j]))
                    .or_default()
                    .push((i, j));
            }
        }
        for (m, c) in target.terms() {
            if !pairs.contains_key(m) && (c.lo() != 0.0 || c.hi() != 0.0) {
                return Err(VerifyError::Unrepresentable {
                    identity: ki,
                    monomial: m.exponents().to_vec(),
                });
            }
        }
        let mut g = point_matrix(gram, n);
        for (m, ps) in &pairs {
            // Designate a diagonal pair when one exists (absorption on the
            // diagonal perturbs definiteness least).
            let di = ps.iter().position(|(i, j)| i == j).unwrap_or(0);
            let (a, b) = ps[di];
            let mut s = Interval::point(0.0);
            for (k, &(i, j)) in ps.iter().enumerate() {
                if k == di {
                    continue;
                }
                let w = if i == j { 1.0 } else { 2.0 };
                s = s + Interval::point(gram[i * n + j]) * Interval::point(w);
            }
            let wd = if a == b { 1.0 } else { 2.0 };
            let entry = (target.coeff(m) - s) / Interval::point(wd);
            max_absorbed_radius = max_absorbed_radius.max(entry.rad());
            g[a * n + b] = entry;
            g[b * n + a] = entry;
        }
        match interval_cholesky_pd(&g, n) {
            Ok(p) => {
                pivot_margins.insert(slack.0, p);
            }
            Err(pivot) => {
                return Err(VerifyError::NotPositiveDefinite {
                    decision: slack.0,
                    pivot,
                });
            }
        }
    }
    Ok(Certificate {
        slacks,
        pivot_margins,
        max_absorbed_radius,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interval_cholesky_accepts_pd_and_rejects_indefinite() {
        let pd = point_matrix(&[2.0, 1.0, 1.0, 2.0], 2);
        assert!(interval_cholesky_pd(&pd, 2).is_ok());
        let indef = point_matrix(&[1.0, 2.0, 2.0, 1.0], 2);
        assert_eq!(interval_cholesky_pd(&indef, 2), Err(1));
        // A wide interval that contains an indefinite member must fail.
        let mut wide = point_matrix(&[1.0, 0.0, 0.0, 1.0], 2);
        wide[1] = Interval::new(-1.5, 1.5);
        wide[2] = Interval::new(-1.5, 1.5);
        assert!(interval_cholesky_pd(&wide, 2).is_err());
        // Singular PSD (on the boundary) is refused: definiteness unprovable.
        let sing = point_matrix(&[1.0, 1.0, 1.0, 1.0], 2);
        assert!(interval_cholesky_pd(&sing, 2).is_err());
    }
}
