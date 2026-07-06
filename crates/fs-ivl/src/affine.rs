//! Affine arithmetic (plan §6.4): values as `c₀ + Σ cᵢ·εᵢ` with noise
//! symbols εᵢ ∈ [−1, 1]. Because the SAME symbol appears on both sides of
//! a correlated expression, x − x collapses to (nearly) zero instead of
//! doubling like plain intervals do — the dependency problem that makes
//! deep F-rep DAG evaluation infeasible for interval arithmetic is exactly
//! what this form kills.
//!
//! RIGOR: every rounding error and every nonlinear-term approximation is
//! ABSORBED into the `err` slack (an anonymous ±err symbol), with the
//! absorption itself rounded up. The containment law is the same as
//! [`crate::Interval`]'s: `to_interval()` always encloses the true value
//! set. Tested against the interval path and the dd oracle.

use crate::Interval;

/// One rounding-error bound for a computed value: ½ ulp of the magnitude,
/// rounded up, plus a subnormal floor. Deliberately generous — err is
/// slack, not signal.
fn round_err(v: f64) -> f64 {
    (v.abs() * (f64::EPSILON * 0.5) + f64::MIN_POSITIVE).next_up()
}

/// Deterministic noise-symbol supply. Symbol identity IS correlation:
/// affine forms only cancel where they share symbols, so all correlated
/// inputs must be built through the same context (deterministic ids — no
/// global state, replay-stable).
#[derive(Debug, Clone, Default)]
pub struct AffineCtx {
    next: u32,
}

impl AffineCtx {
    /// Fresh context (symbols start at 0).
    #[must_use]
    pub fn new() -> AffineCtx {
        AffineCtx::default()
    }

    /// Lift an interval to an affine form with a FRESH noise symbol.
    /// The affine form rigorously encloses the interval: center and radius
    /// are computed with outward absorption.
    pub fn from_interval(&mut self, x: Interval) -> Affine {
        let mid = x.midpoint();
        // Radius covering both endpoints from the computed midpoint,
        // rounded up (mid is not exactly central).
        let rad = (x.hi() - mid).max(mid - x.lo()).next_up();
        let sym = self.next;
        self.next += 1;
        Affine { center: mid, terms: vec![(sym, rad)], err: 0.0 }
    }
}

/// An affine form `center + Σ coeff·ε_sym + [−err, +err]`.
/// Terms are sorted by symbol id (canonical; merges are deterministic).
#[derive(Debug, Clone, PartialEq)]
pub struct Affine {
    center: f64,
    /// (symbol id, coefficient), strictly ascending by id.
    terms: Vec<(u32, f64)>,
    /// Anonymous slack: rounding + nonlinear absorption, always ≥ 0.
    err: f64,
}

impl Affine {
    /// A constant (no noise symbols, no slack).
    #[must_use]
    pub fn constant(c: f64) -> Affine {
        Affine { center: c, terms: Vec::new(), err: 0.0 }
    }

    /// The center value.
    #[must_use]
    pub fn center(&self) -> f64 {
        self.center
    }

    /// Total radius: Σ|coeff| + err, rounded up per accumulation step.
    #[must_use]
    pub fn radius(&self) -> f64 {
        let mut r = self.err;
        for &(_, c) in &self.terms {
            r = (r + c.abs()).next_up();
        }
        r
    }

    /// Collapse to a rigorous interval enclosure.
    #[must_use]
    pub fn to_interval(&self) -> Interval {
        let r = self.radius();
        Interval::new((self.center - r).next_down(), (self.center + r).next_up())
    }

    /// Negation (exact).
    #[must_use]
    pub fn neg(&self) -> Affine {
        Affine {
            center: -self.center,
            terms: self.terms.iter().map(|&(s, c)| (s, -c)).collect(),
            err: self.err,
        }
    }

    /// Addition: symbol-wise coefficient merge; every computed coefficient
    /// absorbs its rounding error into `err`.
    #[must_use]
    pub fn add(&self, o: &Affine) -> Affine {
        let center = self.center + o.center;
        let mut err = (self.err + o.err).next_up();
        err = (err + round_err(center)).next_up();
        let mut terms = Vec::with_capacity(self.terms.len() + o.terms.len());
        let (mut i, mut j) = (0usize, 0usize);
        while i < self.terms.len() || j < o.terms.len() {
            let sa = self.terms.get(i).map_or(u32::MAX, |t| t.0);
            let sb = o.terms.get(j).map_or(u32::MAX, |t| t.0);
            let s = sa.min(sb);
            let mut c = 0.0f64;
            if sa == s {
                c = self.terms[i].1;
                i += 1;
            }
            if sb == s {
                let merged = c + o.terms[j].1;
                err = (err + round_err(merged)).next_up();
                c = merged;
                j += 1;
            }
            if c != 0.0 {
                terms.push((s, c));
            }
        }
        Affine { center, terms, err }
    }

    /// Subtraction. THE point of affine arithmetic: `x.sub(&x)` cancels
    /// symbol-wise to (nearly) exactly zero.
    #[must_use]
    pub fn sub(&self, o: &Affine) -> Affine {
        self.add(&o.neg())
    }

    /// Scale by a constant.
    #[must_use]
    pub fn scale(&self, k: f64) -> Affine {
        let center = self.center * k;
        let mut err = (self.err * k.abs()).next_up();
        err = (err + round_err(center)).next_up();
        let terms = self
            .terms
            .iter()
            .map(|&(s, c)| {
                let ck = c * k;
                err = (err + round_err(ck)).next_up();
                (s, ck)
            })
            .collect();
        Affine { center, terms, err }
    }

    /// Multiplication (standard first-order AA): the bilinear noise·noise
    /// residue is bounded by rad(x)·rad(y) and absorbed into `err`.
    #[must_use]
    pub fn mul(&self, o: &Affine) -> Affine {
        let center = self.center * o.center;
        let rx = self.radius();
        let ry = o.radius();
        // Nonlinear absorption + center rounding.
        let mut err = (rx * ry).next_up();
        err = (err + round_err(center)).next_up();
        // Cross terms: |x0|·err_y + |y0|·err_x are inside rx·ry only for the
        // noise parts; the err components must be added explicitly.
        err = (err + self.center.abs() * o.err).next_up();
        err = (err + o.center.abs() * self.err).next_up();
        // Linear terms: x0·yi + y0·xi, merged symbol-wise.
        let sx = self.scale_terms_only(o.center, &mut err);
        let sy = o.scale_terms_only(self.center, &mut err);
        let merged = merge_terms(&sx, &sy, &mut err);
        Affine { center, terms: merged, err }
    }

    /// Helper: terms scaled by k (center untouched), rounding absorbed.
    fn scale_terms_only(&self, k: f64, err: &mut f64) -> Vec<(u32, f64)> {
        self.terms
            .iter()
            .map(|&(s, c)| {
                let ck = c * k;
                *err = (*err + round_err(ck)).next_up();
                (s, ck)
            })
            .collect()
    }
}

/// Merge two sorted term lists, absorbing merge rounding into `err`.
fn merge_terms(a: &[(u32, f64)], b: &[(u32, f64)], err: &mut f64) -> Vec<(u32, f64)> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0usize, 0usize);
    while i < a.len() || j < b.len() {
        let sa = a.get(i).map_or(u32::MAX, |t| t.0);
        let sb = b.get(j).map_or(u32::MAX, |t| t.0);
        let s = sa.min(sb);
        let mut c = 0.0f64;
        if sa == s {
            c = a[i].1;
            i += 1;
        }
        if sb == s {
            let merged = c + b[j].1;
            *err = (*err + round_err(merged)).next_up();
            c = merged;
            j += 1;
        }
        if c != 0.0 {
            out.push((s, c));
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        ((*seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
    }

    #[test]
    fn x_minus_x_collapses() {
        let mut ctx = AffineCtx::new();
        let x = ctx.from_interval(Interval::new(1.0, 2.0));
        let diff = x.sub(&x).to_interval();
        // Plain intervals give width 1; affine must be ~machine-epsilon.
        assert!(diff.contains(0.0), "must contain the true value 0");
        assert!(diff.width() < 1e-13, "affine x - x too wide: {diff:?}");
        let plain = Interval::new(1.0, 2.0).sub(Interval::new(1.0, 2.0));
        let ratio = plain.width() / diff.width().max(f64::MIN_POSITIVE);
        assert!(ratio > 1e10, "tightness ratio only {ratio}");
        println!(
            "{{\"suite\":\"fs-ivl\",\"case\":\"affine-cancel\",\"verdict\":\"pass\",\"detail\":\"x-x tightness ratio {ratio:.2e} vs plain intervals\"}}"
        );
    }

    #[test]
    fn correlated_chain_is_tighter_than_intervals() {
        // f(x) = (1 + x)·(1 − x) on x ∈ [−0.5, 0.5]: true range [0.75, 1].
        // Plain IA: (1+x) ∈ [0.5, 1.5], (1−x) ∈ [0.5, 1.5] → product
        // [0.25, 2.25] (width 2). AA keeps the ±x correlation.
        let mut ctx = AffineCtx::new();
        let x = ctx.from_interval(Interval::new(-0.5, 0.5));
        let one = Affine::constant(1.0);
        let aa = one.add(&x).mul(&one.sub(&x)).to_interval();
        let xi = Interval::new(-0.5, 0.5);
        let ia = Interval::point(1.0).add(xi).mul(Interval::point(1.0).sub(xi));
        // Containment first (both must hold), tightness second.
        for t in 0..=10 {
            let p = -0.5 + f64::from(t) / 10.0;
            let truth = (1.0 + p) * (1.0 - p);
            assert!(aa.contains(truth), "AA lost truth at {p}");
            assert!(ia.contains(truth), "IA lost truth at {p}");
        }
        assert!(
            aa.width() < ia.width() * 0.7,
            "AA ({}) not meaningfully tighter than IA ({})",
            aa.width(),
            ia.width()
        );
        println!(
            "{{\"suite\":\"fs-ivl\",\"case\":\"affine-tightness\",\"verdict\":\"pass\",\"detail\":\"(1+x)(1-x): AA width {:.3} vs IA width {:.3}\"}}",
            aa.width(),
            ia.width()
        );
    }

    #[test]
    fn affine_containment_vs_point_oracle() {
        // Random affine expression DAGs must contain point evaluations.
        let mut seed = 0xAFF_u64;
        for _ in 0..5_000 {
            let (a, b) = (lcg(&mut seed) * 4.0, lcg(&mut seed) * 4.0);
            let (wa, wb) =
                ((lcg(&mut seed) + 0.5).abs() + 0.01, (lcg(&mut seed) + 0.5).abs() + 0.01);
            let ix = Interval::new(a - wa, a + wa);
            let iy = Interval::new(b - wb, b + wb);
            let mut ctx = AffineCtx::new();
            let x = ctx.from_interval(ix);
            let y = ctx.from_interval(iy);
            // f = (x+y)·(x−y) − x·x + y·y  ≡ 0 identically!
            let f = x.add(&y).mul(&x.sub(&y)).sub(&x.mul(&x)).add(&y.mul(&y));
            let enc = f.to_interval();
            assert!(enc.contains(0.0), "identity-zero not contained: {enc:?}");
            // AA is FIRST-order: linear (ε) terms cancel exactly, so the
            // residual width is O(radius²) — and, decisively, independent of
            // the CENTER magnitude (plain IA's error grows with |center|).
            let rad_budget = 20.0 * (wa + wb) * (wa + wb) + 1e-9;
            assert!(
                enc.width() < rad_budget,
                "identity-zero wider than the quadratic-in-radius bound: {enc:?} \
                 (radii {wa}, {wb}, centers {a}, {b})"
            );
        }
    }

    #[test]
    fn symbol_identity_is_correlation() {
        // Two DIFFERENT symbols over the same interval must NOT cancel.
        let mut ctx = AffineCtx::new();
        let x1 = ctx.from_interval(Interval::new(1.0, 2.0));
        let x2 = ctx.from_interval(Interval::new(1.0, 2.0));
        let diff = x1.sub(&x2).to_interval();
        assert!(diff.width() > 0.9, "independent symbols must not cancel: {diff:?}");
        assert!(diff.contains(0.9) && diff.contains(-0.9));
    }

    #[test]
    fn scale_and_constants() {
        let mut ctx = AffineCtx::new();
        let x = ctx.from_interval(Interval::new(-1.0, 3.0));
        let s = x.scale(-2.0).to_interval();
        assert!(s.contains(-6.0) && s.contains(2.0));
        assert!(s.lo() <= -6.0 && s.hi() >= 2.0 && s.width() < 8.1);
        let c = Affine::constant(2.5);
        assert_eq!(c.to_interval().midpoint(), 2.5);
        assert!(c.to_interval().width() < 1e-15);
    }
}
