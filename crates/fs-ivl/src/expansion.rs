//! Floating-point expansion arithmetic (Shewchuk 1997, plan §6.4): EXACT
//! operations on numbers represented as sums of nonoverlapping f64
//! components, built on the shared EFT primitives in `fs_math::eft`.
//!
//! Representation: an expansion is a slice of NONZERO components in
//! increasing magnitude whose exact real sum is the value; the empty slice
//! is exactly zero (all constructors here zero-eliminate). For expansions
//! produced by these operations the components are nonoverlapping, so the
//! SIGN of the value is the sign of the largest (last) component — the
//! property the exact predicates stand on.
//!
//! Every operation is error-free: the output's exact sum equals the exact
//! real-arithmetic result. That is bitwise-testable algebra (the G0 suite
//! reconstructs sums through a double-double ladder and an i128 lattice
//! oracle), not approximation.
//!
//! Determinism: pure +,−,×,`mul_add` — cross-ISA bit-deterministic by
//! construction.

use fs_math::eft::{quick_two_sum, two_prod, two_sum};

/// Exact difference: `(d, tail)` with `d = fl(a−b)` and `d + tail = a − b`
/// exactly (Shewchuk's Two-Diff, via [`two_sum`] on the negation).
#[must_use]
pub fn two_diff(a: f64, b: f64) -> (f64, f64) {
    two_sum(a, -b)
}

/// Push a component, eliminating zeros (empty output = exact zero).
fn push_nonzero(h: &mut Vec<f64>, x: f64) {
    if x != 0.0 {
        h.push(x);
    }
}

/// Exact sum of two expansions (Shewchuk's FAST-EXPANSION-SUM with zero
/// elimination). Inputs must be valid expansions (nonoverlapping,
/// increasing magnitude); the output is one too.
#[must_use]
pub fn fast_expansion_sum_zeroelim(e: &[f64], f: &[f64]) -> Vec<f64> {
    if e.is_empty() {
        return f.to_vec();
    }
    if f.is_empty() {
        return e.to_vec();
    }
    // Merge by increasing magnitude. `(fnow > enow) == (fnow > -enow)`
    // is Shewchuk's branchless |fnow| > |enow| test.
    let mut h = Vec::with_capacity(e.len() + f.len());
    let (mut ei, mut fi) = (0usize, 0usize);
    let smaller_is_e = |ei: usize, fi: usize| -> bool {
        // Shewchuk's branchless |f| > |e| test: true when e[ei] has the
        // smaller magnitude (take it first).
        (f[fi] > e[ei]) == (f[fi] > -e[ei])
    };
    let mut q = if smaller_is_e(ei, fi) {
        ei += 1;
        e[ei - 1]
    } else {
        fi += 1;
        f[fi - 1]
    };
    let mut first = true;
    while ei < e.len() && fi < f.len() {
        let now = if smaller_is_e(ei, fi) {
            ei += 1;
            e[ei - 1]
        } else {
            fi += 1;
            f[fi - 1]
        };
        let (qnew, hh) = if first {
            // The second-smallest component dominates the smallest.
            quick_two_sum(now, q)
        } else {
            two_sum(q, now)
        };
        first = false;
        q = qnew;
        push_nonzero(&mut h, hh);
    }
    for &now in &e[ei..] {
        let (qnew, hh) = two_sum(q, now);
        q = qnew;
        push_nonzero(&mut h, hh);
    }
    for &now in &f[fi..] {
        let (qnew, hh) = two_sum(q, now);
        q = qnew;
        push_nonzero(&mut h, hh);
    }
    push_nonzero(&mut h, q);
    h
}

/// Exact product of an expansion by one f64 (Shewchuk's SCALE-EXPANSION
/// with zero elimination). FMA-backed: no Dekker splitting.
#[must_use]
pub fn scale_expansion_zeroelim(e: &[f64], b: f64) -> Vec<f64> {
    let Some((&e0, rest)) = e.split_first() else {
        return Vec::new();
    };
    let mut h = Vec::with_capacity(2 * e.len());
    let (mut q, hh) = two_prod(e0, b);
    push_nonzero(&mut h, hh);
    for &enow in rest {
        let (product1, product0) = two_prod(enow, b);
        let (sum, hh) = two_sum(q, product0);
        push_nonzero(&mut h, hh);
        // |product1| >= |sum| holds by Shewchuk's proof (Theorem 19).
        let (qnew, hh) = quick_two_sum(product1, sum);
        q = qnew;
        push_nonzero(&mut h, hh);
    }
    push_nonzero(&mut h, q);
    h
}

/// Exact product of two expansions: distribute [`scale_expansion_zeroelim`]
/// over one operand's components and distill with pairwise exact sums.
/// Sizes grow as O(|a|·|b|) — this is the cold exact tail of the predicate
/// ladder, never the hot path.
#[must_use]
pub fn expansion_product(a: &[f64], b: &[f64]) -> Vec<f64> {
    match a.len() {
        0 => Vec::new(),
        1 => scale_expansion_zeroelim(b, a[0]),
        _ => {
            let mid = a.len() / 2;
            let lo = expansion_product(&a[..mid], b);
            let hi = expansion_product(&a[mid..], b);
            fast_expansion_sum_zeroelim(&lo, &hi)
        }
    }
}

/// Exact difference of two expansions (`a − b`).
#[must_use]
pub fn expansion_diff(a: &[f64], b: &[f64]) -> Vec<f64> {
    let neg: Vec<f64> = b.iter().map(|&x| -x).collect();
    fast_expansion_sum_zeroelim(a, &neg)
}

/// The sign of an expansion's exact value: the sign of its largest
/// component (valid because components are nonoverlapping), 0 when empty.
#[must_use]
pub fn expansion_sign(e: &[f64]) -> i32 {
    match e.last() {
        Some(&x) if x > 0.0 => 1,
        Some(&x) if x < 0.0 => -1,
        _ => 0,
    }
}

/// Approximate value (plain left-to-right sum; adequate for error-bound
/// comparisons in the adaptive stages).
#[must_use]
pub fn estimate(e: &[f64]) -> f64 {
    e.iter().sum()
}

/// The exact 2×2 determinant `a·b − c·d` as an expansion of ≤ 4 components
/// (Shewchuk's Two-Product difference, the workhorse minor).
#[must_use]
pub fn prod_diff(a: f64, b: f64, c: f64, d: f64) -> Vec<f64> {
    let (p1, p0) = two_prod(a, b);
    let (q1, q0) = two_prod(c, d);
    fast_expansion_sum_zeroelim(&nonzero2(p0, p1), &nonzero2(-q0, -q1))
}

/// Two-component expansion `(lo, hi)` with zeros eliminated.
fn nonzero2(lo: f64, hi: f64) -> Vec<f64> {
    let mut v = Vec::with_capacity(2);
    if lo != 0.0 {
        v.push(lo);
    }
    if hi != 0.0 {
        v.push(hi);
    }
    v
}

/// An exact coordinate difference as a ≤ 2-component expansion.
#[must_use]
pub fn diff_expansion(a: f64, b: f64) -> Vec<f64> {
    let (d, tail) = two_diff(a, b);
    nonzero2(tail, d)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
    }

    /// The exactness law in its native form: an EXACT residual. If op
    /// results are error-free, `result − inputs` distills to the empty
    /// expansion — no tolerance, no oracle precision ceiling.
    fn assert_exactly_zero(e: &[f64], law: &str) {
        assert_eq!(expansion_sign(e), 0, "{law}: nonzero residual {e:?}");
        assert!(
            e.is_empty(),
            "{law}: zero-sign residual with components {e:?}"
        );
    }

    #[test]
    fn sum_scale_product_satisfy_exact_residual_laws() {
        let mut seed = 0x5EED_E4BA_0000_0001u64;
        for round in 0..500 {
            // Wildly mixed magnitudes force nontrivial tails everywhere.
            let a = lcg(&mut seed);
            let b = lcg(&mut seed) * 1e-17;
            let c = lcg(&mut seed);
            let d = lcg(&mut seed) * 1e17;
            let s = lcg(&mut seed) * 1e5;
            let ea = diff_expansion(a, b);
            let eb = diff_expansion(c, d);
            // Sum law: (ea + eb) − ea − eb == 0 exactly.
            let sum = fast_expansion_sum_zeroelim(&ea, &eb);
            let residual = expansion_diff(&expansion_diff(&sum, &ea), &eb);
            assert_exactly_zero(&residual, &format!("sum law, round {round}"));
            // Scale ≡ product-by-singleton, exactly.
            let scaled = scale_expansion_zeroelim(&ea, s);
            let prod = expansion_product(&ea, &[s]);
            assert_exactly_zero(
                &expansion_diff(&scaled, &prod),
                &format!("scale/product law, round {round}"),
            );
            // Distributivity: (ea + eb)·ec == ea·ec + eb·ec, exactly.
            let ec = diff_expansion(s, a);
            let lhs = expansion_product(&sum, &ec);
            let rhs = fast_expansion_sum_zeroelim(
                &expansion_product(&ea, &ec),
                &expansion_product(&eb, &ec),
            );
            assert_exactly_zero(
                &expansion_diff(&lhs, &rhs),
                &format!("distributivity, round {round}"),
            );
        }
    }

    #[test]
    fn prod_diff_catches_catastrophic_cancellation() {
        let (a, b) = (1.0 + f64::EPSILON, 1.0 - f64::EPSILON);
        // ab − ba is exactly zero.
        assert_eq!(expansion_sign(&prod_diff(a, b, b, a)), 0);
        // Exact a·b = 1 − ε², but fl(a·b) = 1.0: the difference is the
        // 2⁻¹⁰⁴ the naive float determinant silently throws away.
        let fl_ab = a * b;
        assert_eq!(
            fl_ab.to_bits(),
            1.0f64.to_bits(),
            "premise: float product rounds to 1"
        );
        assert_eq!(expansion_sign(&prod_diff(a, b, 1.0, fl_ab)), -1);
        assert_eq!(expansion_sign(&prod_diff(1.0, fl_ab, a, b)), 1);
    }

    #[test]
    fn empty_expansion_is_exact_zero() {
        assert_eq!(expansion_sign(&[]), 0);
        assert!(fast_expansion_sum_zeroelim(&[], &[]).is_empty());
        assert!(scale_expansion_zeroelim(&[], 7.0).is_empty());
        assert!(expansion_product(&[], &[1.0]).is_empty());
        let x = diff_expansion(3.5, 3.5);
        assert!(x.is_empty(), "a − a is exactly zero: {x:?}");
    }
}
