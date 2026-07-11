//! fs-la batched battery (6ys.4): G0 properties per op (L·Lᵀ
//! reconstruction, PA = LU residual, inv·A = I, det invariants) over
//! random well- and ill-conditioned batches, singular members FLAGGED
//! while the batch continues, batch-size/position bitwise invariance
//! (the structural determinism claim of SoA-across-batch), closed-form
//! eigh3 vs per-matrix Jacobi, plane alignment inherited from fs-soa,
//! and the cross-ISA golden hash.

use fs_la::batched::{
    BatchMat, BatchVec, batch_cholesky, batch_cholesky_solve, batch_det, batch_eigh,
    batch_eigh3_values, batch_gemm, batch_inv, batch_lu,
};
use fs_la::factor::FactorError;
use fs_rand::StreamKey;

fn log(case: &str, verdict: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"{case}\",\"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}"
    );
}

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

fn stream(tile: u32) -> fs_rand::Stream {
    StreamKey {
        seed: 64,
        kernel: 0xBA7C,
        tile,
    }
    .stream()
}

/// Random batch with entries in [-1, 1).
fn random_batch(k: usize, n: usize, tile: u32) -> BatchMat {
    let mut s = stream(tile);
    BatchMat::from_fn(k, n, |_, _, _| 2.0f64.mul_add(s.next_f64(), -1.0))
}

/// Random SPD batch: A = Mᵀ·M + k·I (well-conditioned members).
fn random_spd(k: usize, n: usize, tile: u32) -> BatchMat {
    let m = random_batch(k, n, tile);
    let mut a = BatchMat::zeros(k, n);
    for idx in 0..n {
        let g = m.gather(idx);
        let mut s = vec![0.0f64; k * k];
        for i in 0..k {
            for j in 0..k {
                let mut acc = 0.0f64;
                for l in 0..k {
                    acc = g[l * k + i].mul_add(g[l * k + j], acc);
                }
                s[i * k + j] = acc + if i == j { k as f64 } else { 0.0 };
            }
        }
        a.scatter(idx, &s);
    }
    a
}

#[test]
fn gemm_matches_scalar_reference_across_size_classes() {
    for &k in &[3usize, 4, 6, 8, 12, 16] {
        let n = 37; // odd batch length exercises the non-quantum tail
        let a = random_batch(k, n, 1000 + u32::try_from(k).expect("small"));
        let b = random_batch(k, n, 2000 + u32::try_from(k).expect("small"));
        let mut c = BatchMat::zeros(k, n);
        batch_gemm(1.5, &a, &b, 0.0, &mut c);
        let mut worst = 0.0f64;
        for m in [0usize, n / 2, n - 1] {
            let (ga, gb, gc) = (a.gather(m), b.gather(m), c.gather(m));
            for i in 0..k {
                for j in 0..k {
                    let mut acc = 0.0f64;
                    for l in 0..k {
                        acc = ga[i * k + l].mul_add(gb[l * k + j], acc);
                    }
                    worst = worst.max((1.5 * acc - gc[i * k + j]).abs());
                }
            }
        }
        assert!(
            worst < 1e-13,
            "k={k}: gemm deviates from scalar reference by {worst:.3e}"
        );
        // beta-accumulate path: C := 1·A·B + 2·C doubles then adds.
        let snapshot = c.gather(0);
        let mut c2 = c.clone();
        batch_gemm(1.5, &a, &b, 1.0, &mut c2);
        let after = c2.gather(0);
        let dev = (0..k * k)
            .map(|i| (after[i] - 2.0 * snapshot[i]).abs())
            .fold(0.0f64, f64::max);
        assert!(dev < 1e-12, "beta path wrong by {dev:.3e}");
    }
    log(
        "gemm-reference",
        "pass",
        "size classes 3,4,6,8,12,16 vs scalar oracle",
    );
}

#[test]
fn gemm_alpha_zero_does_not_read_operands() {
    let (k, n) = (4usize, 3usize);
    let a = BatchMat::from_fn(k, n, |_, _, _| f64::NAN);
    let b = BatchMat::from_fn(k, n, |_, _, _| f64::INFINITY);
    let c0 = BatchMat::from_fn(k, n, |m, i, j| (m * 17 + i * 5 + j) as f64 + 0.25);
    for alpha in [0.0, -0.0] {
        for beta in [0.0, 1.0, -0.75] {
            let mut c = c0.clone();
            batch_gemm(alpha, &a, &b, beta, &mut c);
            for m in 0..n {
                for i in 0..k {
                    for j in 0..k {
                        let want = if beta == 0.0 {
                            0.0
                        } else {
                            beta * c0.get(m, i, j)
                        };
                        assert_eq!(c.get(m, i, j).to_bits(), want.to_bits());
                    }
                }
            }
        }
    }
}

#[test]
fn batch_shape_overflow_is_refused() {
    assert_panics_with("batch stride overflow", || {
        let _ = BatchMat::zeros(1, usize::MAX);
    });
    assert_panics_with("batch matrix shape overflow", || {
        let _ = BatchMat::zeros(usize::MAX, 1);
    });
    assert_panics_with("batch stride overflow", || {
        let _ = BatchVec::zeros(1, usize::MAX);
    });
}

#[test]
fn batch_membership_is_bitwise_irrelevant() {
    // THE structural determinism claim: matrix m computed in a batch
    // of N is bitwise-equal to the same matrix alone in a batch of 1
    // — for every op family.
    let k = 6;
    let n = 17;
    let spd = random_spd(k, n, 3000);
    let general = random_batch(k, n, 3001);
    let (l_full, fl) = batch_cholesky(&spd);
    assert!(fl.is_empty());
    let lu_full = batch_lu(&general);
    for m in [0usize, 7, 16] {
        // Singleton batches carrying just matrix m.
        let spd1 = {
            let mut b = BatchMat::zeros(k, 1);
            b.scatter(0, &spd.gather(m));
            b
        };
        let gen1 = {
            let mut b = BatchMat::zeros(k, 1);
            b.scatter(0, &general.gather(m));
            b
        };
        let (l_one, _) = batch_cholesky(&spd1);
        for i in 0..k {
            for j in 0..=i {
                assert_eq!(
                    l_full.get(m, i, j).to_bits(),
                    l_one.get(0, i, j).to_bits(),
                    "cholesky bits differ at ({i},{j}) for member {m}"
                );
            }
        }
        let lu_one = batch_lu(&gen1);
        for i in 0..k {
            for j in 0..k {
                assert_eq!(
                    lu_full.lu.get(m, i, j).to_bits(),
                    lu_one.lu.get(0, i, j).to_bits(),
                    "LU bits differ at ({i},{j}) for member {m}"
                );
            }
            assert_eq!(
                lu_full.perm[i * n + m],
                lu_one.perm[i],
                "pivot differs at step {i} for member {m}"
            );
        }
    }
    log(
        "batch-invariance",
        "pass",
        "cholesky+LU bitwise across batch membership",
    );
}

#[test]
fn cholesky_reconstructs_and_solves() {
    for &k in &[4usize, 6, 12, 24, 48] {
        let n = 21;
        let a = random_spd(k, n, 4000 + u32::try_from(k).expect("small"));
        let (l, flags) = batch_cholesky(&a);
        assert!(flags.is_empty(), "k={k}: unexpected flags {flags:?}");
        // L·Lᵀ = A on sampled members.
        let mut worst = 0.0f64;
        for m in [0usize, n - 1] {
            let (gl, ga) = (l.gather(m), a.gather(m));
            for i in 0..k {
                for j in 0..k {
                    let mut acc = 0.0f64;
                    for p in 0..k {
                        acc = gl[i * k + p].mul_add(gl[j * k + p], acc);
                    }
                    worst = worst.max((acc - ga[i * k + j]).abs());
                }
            }
        }
        assert!(worst < 1e-10 * k as f64, "k={k}: LLt residual {worst:.3e}");
        // Solve against a known x.
        let x_true = BatchVec::from_fn(k, n, |m, i| {
            f64::from(u32::try_from(m + i).expect("small")) * 0.25 + 1.0
        });
        let mut b = BatchVec::zeros(k, n);
        for i in 0..k {
            let mut acc = vec![0.0f64; n];
            for j in 0..k {
                let ap = a.plane(i, j);
                let xp = x_true.plane(j);
                for m in 0..n {
                    acc[m] = ap[m].mul_add(xp[m], acc[m]);
                }
            }
            b.plane_mut(i).copy_from_slice(&acc);
        }
        batch_cholesky_solve(&l, &mut b);
        let mut err = 0.0f64;
        for i in 0..k {
            for m in 0..n {
                err = err.max((b.get(m, i) - x_true.get(m, i)).abs());
            }
        }
        assert!(err < 1e-8, "k={k}: solve error {err:.3e}");
    }
    log("cholesky", "pass", "LLt + solve for k in 4,6,12,24,48");
}

#[test]
fn cholesky_flags_non_spd_and_batch_continues() {
    let k = 6;
    let n = 9;
    let mut a = random_spd(k, n, 5000);
    // Poison member 4 with an indefinite matrix (negative diagonal).
    let mut bad = vec![0.0f64; k * k];
    for i in 0..k {
        bad[i * k + i] = -1.0;
    }
    a.scatter(4, &bad);
    let (l, flags) = batch_cholesky(&a);
    assert_eq!(flags.len(), 1);
    assert_eq!(flags[0].0, 4);
    assert!(matches!(flags[0].1, FactorError::NotSpd { index: 0 }));
    // Every healthy member still factors correctly.
    for m in (0..n).filter(|&m| m != 4) {
        let (gl, ga) = (l.gather(m), a.gather(m));
        let mut worst = 0.0f64;
        for i in 0..k {
            for j in 0..k {
                let mut acc = 0.0f64;
                for p in 0..k {
                    acc = gl[i * k + p].mul_add(gl[j * k + p], acc);
                }
                worst = worst.max((acc - ga[i * k + j]).abs());
            }
        }
        assert!(
            worst < 1e-10,
            "member {m} damaged by poisoned neighbor: {worst:.3e}"
        );
        // Flagged member's factor stays finite (no NaN storm).
        assert!(l.gather(4).iter().all(|v| v.is_finite()));
    }
    log(
        "cholesky-flags",
        "pass",
        "NotSpd flagged, 8/9 healthy members unaffected",
    );
}

#[test]
fn lu_residual_and_singular_flags() {
    let k = 8;
    let n = 13;
    let mut a = random_batch(k, n, 6000);
    // Poison member 2 with an exactly singular matrix (zero column 3).
    let mut g = a.gather(2);
    for i in 0..k {
        g[i * k + 3] = 0.0;
    }
    a.scatter(2, &g);
    let f = batch_lu(&a);
    assert_eq!(f.flags.len(), 1);
    assert_eq!(f.flags[0].0, 2);
    assert!(matches!(f.flags[0].1, FactorError::Singular { .. }));
    // PA = LU residual on healthy members.
    for m in (0..n).filter(|&m| m != 2) {
        let glu = f.lu.gather(m);
        let ga = a.gather(m);
        let mut worst = 0.0f64;
        for i in 0..k {
            for j in 0..k {
                // (L·U)[i][j] with unit-lower L: p runs to min(i, j),
                // hitting the unit diagonal only when i <= j.
                let mut acc = 0.0f64;
                for p in 0..=i.min(j) {
                    let lv = if p == i { 1.0 } else { glu[i * k + p] };
                    acc = lv.mul_add(glu[p * k + j], acc);
                }
                let pa = ga[(f.perm[i * n + m] as usize) * k + j];
                worst = worst.max((acc - pa).abs());
            }
        }
        assert!(worst < 1e-11, "member {m}: PA-LU residual {worst:.3e}");
    }
    // Solve on a healthy member.
    let x_true = BatchVec::from_fn(k, n, |m, i| {
        0.5f64.mul_add(f64::from(u32::try_from(i).expect("small")), 1.0) + m as f64 * 0.01
    });
    let mut b = BatchVec::zeros(k, n);
    for i in 0..k {
        let mut acc = vec![0.0f64; n];
        for j in 0..k {
            let ap = a.plane(i, j);
            let xp = x_true.plane(j);
            for m in 0..n {
                acc[m] = ap[m].mul_add(xp[m], acc[m]);
            }
        }
        b.plane_mut(i).copy_from_slice(&acc);
    }
    f.solve(&mut b);
    let mut err = 0.0f64;
    for i in 0..k {
        for m in (0..n).filter(|&m| m != 2) {
            err = err.max((b.get(m, i) - x_true.get(m, i)).abs());
        }
    }
    assert!(err < 1e-8, "LU solve error {err:.3e}");
    log(
        "lu",
        "pass",
        &format!("PA=LU + solve, singular member flagged, err={err:.2e}"),
    );
}

#[test]
fn det_and_inverse_small() {
    // det: alternating-sign invariant under row swap; inverse: A·A⁻¹ = I.
    for &k in &[1usize, 2, 3, 4] {
        let n = 25;
        let a = random_batch(k, n, 7000 + u32::try_from(k).expect("small"));
        let dets = batch_det(&a);
        // Cross-check determinant against batch_lu diag product for k where both apply.
        if k >= 2 {
            let f = batch_lu(&a);
            for m in [0usize, 12, 24] {
                let mut prod = 1.0f64;
                for i in 0..k {
                    prod *= f.lu.get(m, i, i);
                }
                // Permutation parity.
                let mut perm: Vec<usize> = (0..k).map(|s| f.perm[s * n + m] as usize).collect();
                let mut swaps = 0usize;
                for i in 0..k {
                    while perm[i] != i {
                        let t = perm[i];
                        perm.swap(i, t);
                        swaps += 1;
                    }
                }
                if swaps % 2 == 1 {
                    prod = -prod;
                }
                let rel = (dets[m] - prod).abs() / dets[m].abs().max(1e-30);
                assert!(
                    rel < 1e-9,
                    "k={k} m={m}: det {} vs LU {prod} (rel {rel:.3e})",
                    dets[m]
                );
            }
        }
        let mut inv = BatchMat::zeros(k, n);
        let flags = batch_inv(&a, &mut inv);
        for &(m, _) in &flags {
            assert!(
                dets[m].to_bits() == 0.0f64.to_bits(),
                "flagged member {m} must have zero det: {}",
                dets[m]
            );
        }
        for m in [0usize, 24] {
            if flags.iter().any(|&(fm, _)| fm == m) {
                continue;
            }
            let (ga, gi) = (a.gather(m), inv.gather(m));
            let mut worst = 0.0f64;
            for i in 0..k {
                for j in 0..k {
                    let mut acc = 0.0f64;
                    for l in 0..k {
                        acc = ga[i * k + l].mul_add(gi[l * k + j], acc);
                    }
                    let expect = if i == j { 1.0 } else { 0.0 };
                    worst = worst.max((acc - expect).abs());
                }
            }
            assert!(worst < 1e-9, "k={k} m={m}: A*inv deviates {worst:.3e}");
        }
    }
    // Exactly singular member flagged.
    let mut a = random_batch(3, 5, 7100);
    a.scatter(1, &[1.0, 2.0, 3.0, 2.0, 4.0, 6.0, 0.5, 1.0, 1.5]); // rank 1
    let mut inv = BatchMat::zeros(3, 5);
    let flags = batch_inv(&a, &mut inv);
    assert_eq!(flags.len(), 1);
    assert_eq!(flags[0].0, 1);
    log(
        "det-inv",
        "pass",
        "closed forms k<=4 vs LU cross-check, singular flagged",
    );
}

#[test]
fn eigh3_closed_form_matches_jacobi() {
    let n = 40;
    let a = random_spd(3, n, 8000);
    let vals = batch_eigh3_values(&a);
    let (jvals, jvecs) = batch_eigh(&a);
    for m in 0..n {
        // Compare as sorted sets (Jacobi ordering is its own).
        let mut cf: Vec<f64> = (0..3).map(|i| vals.get(m, i)).collect();
        let mut jj: Vec<f64> = (0..3).map(|i| jvals.get(m, i)).collect();
        cf.sort_by(f64::total_cmp);
        jj.sort_by(f64::total_cmp);
        for i in 0..3 {
            let rel = (cf[i] - jj[i]).abs() / jj[i].abs().max(1e-30);
            assert!(
                rel < 1e-9,
                "m={m} lambda[{i}]: closed {} vs jacobi {} ",
                cf[i],
                jj[i]
            );
        }
        assert!(
            cf[0] <= cf[1] && cf[1] <= cf[2],
            "closed form not ascending at m={m}"
        );
        // Trace and det invariants.
        let g = a.gather(m);
        let tr = g[0] + g[4] + g[8];
        let sum: f64 = cf.iter().sum();
        assert!(
            (sum - tr).abs() < 1e-9 * tr.abs().max(1.0),
            "trace identity violated"
        );
        // Jacobi vectors are genuine eigenvectors: ||A v - lambda v|| small.
        for j in 0..3 {
            let gv = jvecs.gather(m);
            let lambda = jvals.get(m, j);
            let mut resid = 0.0f64;
            for i in 0..3 {
                let mut av = 0.0f64;
                for l in 0..3 {
                    av = g[i * 3 + l].mul_add(gv[l * 3 + j], av);
                }
                resid = resid.max((av - lambda * gv[i * 3 + j]).abs());
            }
            assert!(resid < 1e-9, "m={m}: eigenvector residual {resid:.3e}");
        }
    }
    // Degenerate: exact multiple of identity.
    let mut iso = BatchMat::zeros(3, 2);
    iso.scatter(0, &[2.5, 0.0, 0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 2.5]);
    iso.scatter(1, &[1.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 3.0]);
    let v = batch_eigh3_values(&iso);
    assert_eq!(v.get(0, 0).to_bits(), 2.5f64.to_bits());
    assert_eq!(v.get(0, 2).to_bits(), 2.5f64.to_bits());
    assert!((v.get(1, 0) - 1.0).abs() < 1e-12 && (v.get(1, 2) - 3.0).abs() < 1e-12);
    log(
        "eigh3",
        "pass",
        "closed form vs Jacobi on 40 SPD members + degenerate cases",
    );
}

#[test]
fn planes_are_128_byte_aligned() {
    let a = BatchMat::zeros(6, 37);
    for i in 0..6 {
        for j in 0..6 {
            let addr = a.plane(i, j).as_ptr().addr();
            assert_eq!(addr % 128, 0, "plane ({i},{j}) misaligned");
        }
    }
    let v = BatchVec::zeros(6, 37);
    for i in 0..6 {
        assert_eq!(
            v.plane(i).as_ptr().addr() % 128,
            0,
            "vec plane {i} misaligned"
        );
    }
    log(
        "alignment",
        "pass",
        "all planes 128-byte aligned (fs-soa substrate)",
    );
}

const GOLDEN_HASH: u64 = 0x0377_a8c9_5992_aee9; // recorded at 6ys.4 landing, frozen

#[test]
fn batched_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let k = 6;
    let n = 11;
    let spd = random_spd(k, n, 9000);
    let general = random_batch(k, n, 9001);
    // GEMM.
    let mut c = BatchMat::zeros(k, n);
    batch_gemm(0.75, &spd, &general, 0.0, &mut c);
    for m in [0usize, 10] {
        for v in c.gather(m) {
            feed(v);
        }
    }
    // Cholesky + solve.
    let (l, _) = batch_cholesky(&spd);
    let mut b = BatchVec::from_fn(k, n, |m, i| {
        f64::from(u32::try_from(m + 2 * i).expect("small")) * 0.125 + 1.0
    });
    batch_cholesky_solve(&l, &mut b);
    for i in 0..k {
        feed(b.get(5, i));
    }
    // LU.
    let f = batch_lu(&general);
    for i in 0..k {
        feed(f.lu.get(3, i, i));
        feed(f64::from(f.perm[i * n + 3]));
    }
    // det/inv on 3x3.
    let a3 = random_batch(3, 7, 9002);
    for v in batch_det(&a3) {
        feed(v);
    }
    let mut inv3 = BatchMat::zeros(3, 7);
    let _ = batch_inv(&a3, &mut inv3);
    for v in inv3.gather(6) {
        feed(v);
    }
    // eigh3 closed form.
    let s3 = random_spd(3, 7, 9003);
    let vals = batch_eigh3_values(&s3);
    for m in 0..7 {
        for i in 0..3 {
            feed(vals.get(m, i));
        }
    }
    log("batched-golden", "info", &format!("{acc:#018x}"));
    assert_eq!(
        acc, GOLDEN_HASH,
        "batched bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}
