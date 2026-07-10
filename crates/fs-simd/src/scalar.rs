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

/// Batched-GEMM 4×4 entry-tile microkernel (scalar twin, bead 9ekv)
/// over plane-SoA batches: for tile rows i ∈ i0..i0+4 and columns
/// j ∈ j0..j0+4, `dst[(ti·4+tj)·mb + m] = Σ_l a[(i·k+l)·stride + m0 + m]
/// · b[(l·k+j)·stride + m0 + m]`, l ascending from a zero start, fused —
/// SIMD lanes run across the batch (independent matrices), so capsules
/// must match this BITWISE. α/β write-back stays with the caller.
///
/// # Panics
/// If the plane buffers or `dst` are too short for the tile.
#[allow(clippy::too_many_arguments)] // plane-SoA layout bundle (see fs-la::batched)
pub fn btile4x4_f64(
    a: &[f64],
    b: &[f64],
    i0: usize,
    j0: usize,
    stride: usize,
    k: usize,
    m0: usize,
    mb: usize,
    dst: &mut [f64],
) {
    let bounds = crate::checked_btile4x4_lengths(i0, j0, stride, k, m0, mb);
    assert!(
        matches!(
            bounds,
            Some((a_len, b_len, dst_len))
                if a_len <= a.len() && b_len <= b.len() && dst_len <= dst.len()
        ),
        "btile4x4 plane bounds (programmer error)"
    );
    for ti in 0..4 {
        for tj in 0..4 {
            let drow = &mut dst[(ti * 4 + tj) * mb..(ti * 4 + tj + 1) * mb];
            drow.fill(0.0);
            for l in 0..k {
                let ap = &a[((i0 + ti) * k + l) * stride + m0..][..mb];
                let bp = &b[(l * k + j0 + tj) * stride + m0..][..mb];
                for ((s, &am), &bm) in drow.iter_mut().zip(ap).zip(bp) {
                    *s = am.mul_add(bm, *s);
                }
            }
        }
    }
}

/// The 8×4 f64 GEMM register microkernel (scalar twin): accumulate
/// `acc[r][s] += Σ_kk a_panel[kk·8 + r] · b_panel[kk·4 + s]` with k
/// ascending and fused `mul_add` per element — the bit-contract shape
/// fs-la's packed GEMM is built on. Panels are packed k-fastest
/// (fs-la `pack_a`/`pack_b` layout); `acc` is NOT zeroed here so KC
/// chunks can fold in caller-chosen order.
pub fn mk8x4_f64(a_panel: &[f64], b_panel: &[f64], kc: usize, acc: &mut [[f64; 4]; 8]) {
    let lengths = crate::checked_mk8x4_lengths(kc);
    assert!(
        matches!(lengths, Some((a_len, b_len)) if a_len <= a_panel.len() && b_len <= b_panel.len()),
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

/// Radix-4 Stockham q-run butterfly (scalar twin, bead 27d3) over
/// INTERLEAVED complex rows (re, im pairs): `a..d` are the four source
/// runs of one (stage, p) butterfly group, `out` the contiguous block
/// of the four destination runs (X0|X1|X2|X3). Twiddles arrive as
/// [w1re, w1im, w2re, w2im, w3re, w3im] (already conjugated for the
/// inverse); `inverse` flips the ∓i rotation of (b − d). Per element
/// this is EXACTLY fs-fft's C64 add/sub/mul composition (fused re part,
/// same operand order), so capsules must match it bitwise.
///
/// # Panics
/// If run lengths mismatch or are odd.
pub fn r4qrun_f64(
    a: &[f64],
    b: &[f64],
    c: &[f64],
    d: &[f64],
    out: &mut [f64],
    w: &[f64; 6],
    inverse: bool,
) {
    let s2 = a.len();
    assert!(
        s2.is_multiple_of(2)
            && b.len() == s2
            && c.len() == s2
            && d.len() == s2
            && out.len() == 4 * s2,
        "r4qrun run-length mismatch (programmer error)"
    );
    let (o01, o23) = out.split_at_mut(2 * s2);
    let (o0, o1) = o01.split_at_mut(s2);
    let (o2, o3) = o23.split_at_mut(s2);
    for q in 0..s2 / 2 {
        let (i0, i1) = (2 * q, 2 * q + 1);
        let (t0re, t0im) = (a[i0] + c[i0], a[i1] + c[i1]);
        let (t1re, t1im) = (a[i0] - c[i0], a[i1] - c[i1]);
        let (t2re, t2im) = (b[i0] + d[i0], b[i1] + d[i1]);
        let (t3re, t3im) = (b[i0] - d[i0], b[i1] - d[i1]);
        // ∓i·t3: forward (t3im, −t3re), inverse (−t3im, t3re).
        let (t3ire, t3iim) = if inverse {
            (-t3im, t3re)
        } else {
            (t3im, -t3re)
        };
        o0[i0] = t0re + t2re;
        o0[i1] = t0im + t2im;
        let (u1re, u1im) = (t1re + t3ire, t1im + t3iim);
        o1[i0] = u1re.mul_add(w[0], -(u1im * w[1]));
        o1[i1] = u1re.mul_add(w[1], u1im * w[0]);
        let (u2re, u2im) = (t0re - t2re, t0im - t2im);
        o2[i0] = u2re.mul_add(w[2], -(u2im * w[3]));
        o2[i1] = u2re.mul_add(w[3], u2im * w[2]);
        let (u3re, u3im) = (t1re - t3ire, t1im - t3iim);
        o3[i0] = u3re.mul_add(w[4], -(u3im * w[5]));
        o3[i1] = u3re.mul_add(w[5], u3im * w[4]);
    }
}

/// PACKED batched-GEMM 4×4 tile microkernel (scalar twin, bead 9ekv
/// slice 2) over l-CONTIGUOUS packed operands: A packed i-major
/// (`a[(i·k + l)·mb + lane]`), B packed j-major (`b[(j·k + l)·mb +
/// lane]`), so BOTH walks advance by `mb` per l — contiguous streams,
/// no strided-plane TLB pressure. Per element identical to
/// [`btile4x4_f64`]: zero start, l-ascending fused accumulate —
/// capsules must match bitwise. `dst` is 16 rows of `mb` (row t =
/// tile (t/4, t%4)).
///
/// # Panics
/// If the packed buffers or `dst` are too short for the tile.
#[allow(clippy::too_many_arguments)] // packed-layout bundle (see fs-la::batched)
pub fn btile4x4p_f64(
    a: &[f64],
    b: &[f64],
    i0: usize,
    j0: usize,
    k: usize,
    mb: usize,
    dst: &mut [f64],
) {
    let bounds = crate::checked_btile4x4p_lengths(i0, j0, k, mb);
    assert!(
        matches!(
            bounds,
            Some((a_len, b_len, dst_len))
                if a_len <= a.len() && b_len <= b.len() && dst_len <= dst.len()
        ),
        "btile4x4p packed bounds (programmer error)"
    );
    for ti in 0..4 {
        for tj in 0..4 {
            let drow = &mut dst[(ti * 4 + tj) * mb..(ti * 4 + tj + 1) * mb];
            drow.fill(0.0);
            for l in 0..k {
                let ap = &a[((i0 + ti) * k + l) * mb..][..mb];
                let bp = &b[((j0 + tj) * k + l) * mb..][..mb];
                for ((s, &am), &bm) in drow.iter_mut().zip(ap).zip(bp) {
                    *s = am.mul_add(bm, *s);
                }
            }
        }
    }
}
