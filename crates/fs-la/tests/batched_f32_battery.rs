//! Batched f32 + mixed-precision battery (bead 9ekv, scope e): bitwise
//! scalar-oracle equality, batch-membership invariance, the
//! narrow-once mixed contract (measurably tighter than pure f32 on a
//! cancellation fixture), β = 0 NaN-overwrite semantics, and the frozen
//! golden (registered in golden-couplings.json against
//! `fs-la:batched-f32-bits`).

use fs_la::batched_f32::{BatchMatF32, batch_gemm_f32, batch_gemm_mixed};

fn assert_panics_with(expected: &str, f: impl FnOnce()) {
    let payload = std::panic::catch_unwind(std::panic::AssertUnwindSafe(f))
        .expect_err("overflowing shape unexpectedly succeeded");
    let message = payload
        .downcast_ref::<&str>()
        .copied()
        .or_else(|| payload.downcast_ref::<String>().map(String::as_str))
        .unwrap_or("non-string panic payload");
    assert!(
        message.contains(expected),
        "panic {message:?} did not contain {expected:?}"
    );
}

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
}

#[allow(clippy::cast_possible_truncation)]
fn fixture(k: usize, n: usize, seed: u64) -> (BatchMatF32, BatchMatF32, BatchMatF32) {
    let mut s = seed;
    let mut make = |_m: usize, _i: usize, _j: usize| (lcg(&mut s) * 4.0) as f32;
    let a = BatchMatF32::from_fn(k, n, &mut make);
    let b = BatchMatF32::from_fn(k, n, &mut make);
    let c = BatchMatF32::from_fn(k, n, &mut make);
    (a, b, c)
}

/// The mixed kernel's per-element oracle: the same chain, recomputed
/// scalar per matrix — MUST match bitwise.
#[allow(clippy::cast_possible_truncation)]
fn mixed_oracle(
    (alpha, beta): (f64, f64),
    a: &BatchMatF32,
    b: &BatchMatF32,
    c_old: &BatchMatF32,
    (m, i, j): (usize, usize, usize),
) -> f32 {
    let k = a.k();
    let mut s = 0.0f64;
    for l in 0..k {
        s = f64::from(a.get(m, i, l)).mul_add(f64::from(b.get(m, l, j)), s);
    }
    if beta == 0.0 {
        (alpha * s) as f32
    } else {
        alpha.mul_add(s, beta * f64::from(c_old.get(m, i, j))) as f32
    }
}

fn f32_oracle(
    (alpha, beta): (f32, f32),
    a: &BatchMatF32,
    b: &BatchMatF32,
    c_old: &BatchMatF32,
    (m, i, j): (usize, usize, usize),
) -> f32 {
    let mut s = 0.0f32;
    for l in 0..a.k() {
        s = a.get(m, i, l).mul_add(b.get(m, l, j), s);
    }
    if beta == 0.0 {
        alpha * s
    } else {
        alpha.mul_add(s, beta * c_old.get(m, i, j))
    }
}

#[test]
fn mixed_matches_scalar_oracle_bitwise() {
    for &(k, n) in &[(4usize, 33usize), (8, 64), (12, 7), (24, 40)] {
        let (a, b, c0) = fixture(k, n, 0xF32_u64 + k as u64);
        for &(alpha, beta) in &[(1.0f64, 0.0f64), (0.5, 1.25)] {
            let mut c = c0.clone();
            batch_gemm_mixed(alpha, &a, &b, beta, &mut c);
            for m in 0..n {
                for i in 0..k {
                    for j in 0..k {
                        let want = mixed_oracle((alpha, beta), &a, &b, &c0, (m, i, j));
                        assert_eq!(
                            c.get(m, i, j).to_bits(),
                            want.to_bits(),
                            "mixed k={k} m={m} ({i},{j}) α={alpha} β={beta}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn f32_matches_scalar_oracle_bitwise() {
    for &(k, n) in &[(4usize, 33usize), (8, 17), (12, 7)] {
        let (a, b, c0) = fixture(k, n, 0xF320_u64 + k as u64);
        for &(alpha, beta) in &[(1.0f32, 0.0f32), (0.5, 1.25), (-0.75, -0.25)] {
            let mut c = c0.clone();
            batch_gemm_f32(alpha, &a, &b, beta, &mut c);
            for m in 0..n {
                for i in 0..k {
                    for j in 0..k {
                        let want = f32_oracle((alpha, beta), &a, &b, &c0, (m, i, j));
                        assert_eq!(
                            c.get(m, i, j).to_bits(),
                            want.to_bits(),
                            "f32 k={k} m={m} ({i},{j}) alpha={alpha} beta={beta}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn f32_and_mixed_are_batch_membership_invariant() {
    let (k, n) = (8usize, 48usize);
    let (a, b, c0) = fixture(k, n, 0xACE);
    let mut c = c0.clone();
    batch_gemm_f32(1.0, &a, &b, 0.0, &mut c);
    // Membership invariance: matrix m alone in a batch of 1 must equal
    // matrix m inside the batch of n, bitwise — for BOTH kernels.
    for m in [0usize, 17, n - 1] {
        let a1 = BatchMatF32::from_fn(k, 1, |_, i, j| a.get(m, i, j));
        let b1 = BatchMatF32::from_fn(k, 1, |_, i, j| b.get(m, i, j));
        let mut c1 = BatchMatF32::zeros(k, 1);
        batch_gemm_f32(1.0, &a1, &b1, 0.0, &mut c1);
        let mut cm = BatchMatF32::from_fn(k, 1, |_, i, j| c0.get(m, i, j));
        batch_gemm_mixed(0.75, &a1, &b1, 0.5, &mut cm);
        let mut c_full = c0.clone();
        batch_gemm_mixed(0.75, &a, &b, 0.5, &mut c_full);
        for i in 0..k {
            for j in 0..k {
                assert_eq!(
                    c1.get(0, i, j).to_bits(),
                    c.get(m, i, j).to_bits(),
                    "f32 membership m={m} ({i},{j})"
                );
                assert_eq!(
                    cm.get(0, i, j).to_bits(),
                    c_full.get(m, i, j).to_bits(),
                    "mixed membership m={m} ({i},{j})"
                );
            }
        }
    }
}

#[test]
#[allow(clippy::cast_possible_truncation)]
fn mixed_is_tighter_than_pure_f32_on_cancellation() {
    // Fixture engineered for f32 accumulation error: entries alternate
    // large +/- values whose products nearly cancel; the f64 chain keeps
    // ~29 extra mantissa bits before the single narrow.
    let (k, n) = (24usize, 16usize);
    let a = BatchMatF32::from_fn(k, n, |m, i, l| {
        let sign = if l % 2 == 0 { 1.0f32 } else { -1.0 };
        sign * (1000.0 + (m + i + l) as f32)
    });
    let b = BatchMatF32::from_fn(k, n, |m, l, j| 1.0 + ((m + l + j) as f32) * 1e-3);
    let mut c32 = BatchMatF32::zeros(k, n);
    let mut cmx = BatchMatF32::zeros(k, n);
    batch_gemm_f32(1.0, &a, &b, 0.0, &mut c32);
    batch_gemm_mixed(1.0, &a, &b, 0.0, &mut cmx);
    let (mut err32, mut errmx) = (0.0f64, 0.0f64);
    for m in 0..n {
        for i in 0..k {
            for j in 0..k {
                // f64 reference with the same chain shape.
                let mut s = 0.0f64;
                for l in 0..k {
                    s = f64::from(a.get(m, i, l)).mul_add(f64::from(b.get(m, l, j)), s);
                }
                err32 += (f64::from(c32.get(m, i, j)) - s).abs();
                errmx += (f64::from(cmx.get(m, i, j)) - s).abs();
            }
        }
    }
    assert!(
        errmx < err32 * 0.51,
        "mixed accumulation must at least halve the f32 chain error on the \
         cancellation fixture: mixed {errmx:.3e} vs f32 {err32:.3e}"
    );
}

#[test]
fn beta_zero_overwrites_nan() {
    let (k, n) = (6usize, 5usize);
    let (a, b, _) = fixture(k, n, 0xBAD);
    let mut c = BatchMatF32::from_fn(k, n, |_, _, _| f32::NAN);
    batch_gemm_f32(1.0, &a, &b, 0.0, &mut c);
    let mut cm = BatchMatF32::from_fn(k, n, |_, _, _| f32::NAN);
    batch_gemm_mixed(1.0, &a, &b, 0.0, &mut cm);
    for m in 0..n {
        for i in 0..k {
            for j in 0..k {
                assert!(
                    c.get(m, i, j).is_finite() && cm.get(m, i, j).is_finite(),
                    "β = 0 must overwrite NaN in C (BLAS convention) at m={m} ({i},{j})"
                );
            }
        }
    }
}

#[test]
fn alpha_zero_does_not_read_f32_operands() {
    let (k, n) = (4usize, 3usize);
    let a = BatchMatF32::from_fn(k, n, |_, _, _| f32::NAN);
    let b = BatchMatF32::from_fn(k, n, |_, _, _| f32::INFINITY);
    let c0 = BatchMatF32::from_fn(k, n, |m, i, j| (m * 17 + i * 5 + j) as f32 + 0.25);
    for alpha in [0.0f32, -0.0] {
        for beta in [0.0f32, 1.0, -0.75] {
            let mut c = c0.clone();
            batch_gemm_f32(alpha, &a, &b, beta, &mut c);
            let mut mixed = c0.clone();
            batch_gemm_mixed(f64::from(alpha), &a, &b, f64::from(beta), &mut mixed);
            for m in 0..n {
                for i in 0..k {
                    for j in 0..k {
                        let old = c0.get(m, i, j);
                        let want_f32 = if beta == 0.0 { 0.0 } else { beta * old };
                        let want_mixed = if beta == 0.0 {
                            0.0
                        } else {
                            (f64::from(beta) * f64::from(old)) as f32
                        };
                        assert_eq!(c.get(m, i, j).to_bits(), want_f32.to_bits());
                        assert_eq!(mixed.get(m, i, j).to_bits(), want_mixed.to_bits());
                    }
                }
            }
        }
    }
}

#[test]
fn f32_batch_shape_overflow_is_refused() {
    assert_panics_with("batch stride overflow", || {
        let _ = BatchMatF32::zeros(1, usize::MAX);
    });
    assert_panics_with("batch matrix shape overflow", || {
        let _ = BatchMatF32::zeros(usize::MAX, 1);
    });
}

/// Recorded on aarch64-apple (M4 Pro); f32/f64 fused arithmetic only, so
/// identical across build modes and ISAs by the det doctrine.
const GOLDEN_HASH: u64 = 0x5600_7cfe_6a6d_1f9a;

#[test]
fn batched_f32_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed32 = |v: f32| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    for &(k, n) in &[(4usize, 19usize), (8, 33), (16, 21), (32, 9)] {
        let (a, b, c0) = fixture(k, n, 0x601D + k as u64);
        let mut c = c0.clone();
        batch_gemm_f32(1.25, &a, &b, 0.5, &mut c);
        let mut cm = c0.clone();
        batch_gemm_mixed(1.25, &a, &b, 0.5, &mut cm);
        for m in 0..n {
            for i in 0..k {
                feed32(c.get(m, i, i.min(k - 1)));
                feed32(cm.get(m, i, (i + 1) % k));
            }
        }
    }
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"batched-f32-golden\",\"verdict\":\"info\",\"detail\":\"{acc:#018x}\"}}"
    );
    assert_eq!(
        acc, GOLDEN_HASH,
        "batched f32/mixed bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with \
         semantic justification (golden-evidence policy)"
    );
}
