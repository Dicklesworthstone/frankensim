//! Scalar twins: the portable correctness reference for every SIMD capsule
//! (Tier 0). Every vector primitive has EXACTLY one semantic definition —
//! this one — and capsules must match it bitwise (elementwise ops) or within
//! the documented reduction-shape bound (dot/sum).
//!
//! FMA policy (coordinates with fs-math): elementwise multiply-add uses
//! `f64::mul_add` (FUSED) so scalar and NEON/AVX-FMA tiers agree BITWISE.
//! Unfused fallback would silently diverge per-element — that divergence
//! class belongs to the G5 cross-ISA report, not inside one machine.

/// y[i] = a * x[i] + y[i] (fused).
pub fn axpy(a: f64, x: &[f64], y: &mut [f64]) {
    assert_eq!(x.len(), y.len(), "axpy length mismatch (programmer error)");
    for i in 0..x.len() {
        y[i] = a.mul_add(x[i], y[i]);
    }
}

/// x[i] *= a.
pub fn scale(a: f64, x: &mut [f64]) {
    for v in x {
        *v *= a;
    }
}

/// out[i] = a[i] * b[i].
pub fn mul_elem(a: &[f64], b: &[f64], out: &mut [f64]) {
    assert_eq!(a.len(), b.len(), "mul_elem length mismatch");
    assert_eq!(a.len(), out.len(), "mul_elem length mismatch");
    for i in 0..a.len() {
        out[i] = a[i] * b[i];
    }
}

/// out[i] = a[i] * b[i] + c[i] (fused).
pub fn fma3(a: &[f64], b: &[f64], c: &[f64], out: &mut [f64]) {
    assert_eq!(a.len(), b.len(), "fma3 length mismatch");
    assert_eq!(a.len(), c.len(), "fma3 length mismatch");
    assert_eq!(a.len(), out.len(), "fma3 length mismatch");
    for i in 0..a.len() {
        out[i] = a[i].mul_add(b[i], c[i]);
    }
}

/// Σ x[i]·y[i], SEQUENTIAL accumulation in index order — the scalar tier's
/// fixed reduction shape (each tier's shape is fixed; shapes differ ACROSS
/// tiers within a documented ULP envelope).
#[must_use]
pub fn dot(x: &[f64], y: &[f64]) -> f64 {
    assert_eq!(x.len(), y.len(), "dot length mismatch");
    let mut s = 0.0;
    for i in 0..x.len() {
        s = x[i].mul_add(y[i], s);
    }
    s
}

/// Σ x[i], sequential in index order (the scalar fixed shape).
#[must_use]
pub fn sum(x: &[f64]) -> f64 {
    let mut s = 0.0;
    for &v in x {
        s += v;
    }
    s
}
