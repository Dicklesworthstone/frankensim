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

/// The 8×4 f64 GEMM register microkernel (scalar twin): accumulate
/// `acc[r][s] += Σ_kk a_panel[kk·8 + r] · b_panel[kk·4 + s]` with k
/// ascending and fused `mul_add` per element — the bit-contract shape
/// fs-la's packed GEMM is built on. Panels are packed k-fastest
/// (fs-la `pack_a`/`pack_b` layout); `acc` is NOT zeroed here so KC
/// chunks can fold in caller-chosen order.
pub fn mk8x4_f64(a_panel: &[f64], b_panel: &[f64], kc: usize, acc: &mut [[f64; 4]; 8]) {
    assert!(
        a_panel.len() >= kc * 8 && b_panel.len() >= kc * 4,
        "mk8x4 panel length mismatch (programmer error)"
    );
    for kk in 0..kc {
        let av = &a_panel[kk * 8..kk * 8 + 8];
        let bv = &b_panel[kk * 4..kk * 4 + 4];
        for (accr, &ar) in acc.iter_mut().zip(av) {
            for (slot, &bs) in accr.iter_mut().zip(bv) {
                *slot = ar.mul_add(bs, *slot);
            }
        }
    }
}
