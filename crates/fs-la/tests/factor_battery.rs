//! Factorization battery: factor-reconstruct residuals on random +
//! adversarial (Kahan, Hilbert) matrices, permutation/orthogonality
//! invariants, deterministic pivot tie-breaks, degenerate shapes, and the
//! cross-ISA golden hash — fs-la's G0 suite for the factorization layer.

use fs_la::factor::{FactorError, cholesky, lu, qr, svd_jacobi, tsqr_r};
use fs_la::gemm_f64;

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
}

fn rand_mat(rows: usize, cols: usize, seed: u64) -> Vec<f64> {
    let mut s = seed;
    (0..rows * cols).map(|_| lcg(&mut s)).collect()
}

/// Random SPD: B·Bᵀ + n·I.
fn rand_spd(n: usize, seed: u64) -> Vec<f64> {
    let b = rand_mat(n, n, seed);
    let bt: Vec<f64> = (0..n * n).map(|i| b[(i % n) * n + i / n]).collect();
    let mut a = vec![0.0; n * n];
    gemm_f64(n, n, n, 1.0, &b, &bt, 0.0, &mut a);
    for i in 0..n {
        a[i * n + i] += n as f64;
    }
    a
}

fn frob(m: &[f64]) -> f64 {
    m.iter().map(|v| v * v).sum::<f64>().sqrt()
}

/// Hilbert matrix H[i][j] = 1/(i+j+1) — the classic ill-conditioned case.
fn hilbert(n: usize) -> Vec<f64> {
    (0..n * n).map(|i| 1.0 / ((i / n + i % n + 1) as f64)).collect()
}

/// Kahan matrix (upper triangular, pivot-adversarial): K[i][j] =
/// cⁱ·(1 if i==j else −s), s²+c²=1.
fn kahan(n: usize, theta: f64) -> Vec<f64> {
    let (s, c) = (theta.sin(), theta.cos());
    let mut k = vec![0.0; n * n];
    for i in 0..n {
        let ci = c.powi(i32::try_from(i).unwrap());
        for j in 0..n {
            if j == i {
                k[i * n + j] = ci;
            } else if j > i {
                k[i * n + j] = -s * ci;
            }
        }
    }
    k
}

#[test]
fn cholesky_residual_and_solve() {
    for (n, seed) in [(1usize, 1u64), (7, 2), (33, 3), (96, 4)] {
        let a = rand_spd(n, seed);
        let f = cholesky(&a, n).expect("SPD must factor");
        // ‖A − L·Lᵀ‖_F / ‖A‖_F ≤ c·n·eps.
        let mut llt = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                let mut acc = 0.0f64;
                for k in 0..=i.min(j) {
                    acc = f.l(i, k).mul_add(f.l(j, k), acc);
                }
                llt[i * n + j] = acc;
            }
        }
        let diff: Vec<f64> = a.iter().zip(&llt).map(|(x, y)| x - y).collect();
        let rel = frob(&diff) / frob(&a);
        assert!(rel <= 1e-14 * (n as f64), "chol residual {rel:.2e} at n={n}");
        // Solve check: A·x = b round trip.
        let x_true: Vec<f64> = (0..n).map(|i| (i as f64).sin() + 1.0).collect();
        let mut b = vec![0.0; n];
        for i in 0..n {
            let mut acc = 0.0;
            for j in 0..n {
                acc = a[i * n + j].mul_add(x_true[j], acc);
            }
            b[i] = acc;
        }
        f.solve(&mut b);
        for i in 0..n {
            assert!(
                (b[i] - x_true[i]).abs() < 1e-10,
                "chol solve at {i}: {} vs {}",
                b[i],
                x_true[i]
            );
        }
    }
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"cholesky\",\"verdict\":\"pass\",\"detail\":\"residual+solve n in 1..96\"}}"
    );
}

#[test]
fn cholesky_rejects_indefinite_with_index() {
    // Indefinite: negative eigenvalue reachable at pivot 1.
    let a = vec![4.0, 10.0, 10.0, 4.0]; // det < 0
    match cholesky(&a, 2) {
        Err(FactorError::NotSpd { index }) => assert_eq!(index, 1),
        other => panic!("expected NotSpd, got {other:?}"),
    }
}

#[test]
fn lu_residual_permutation_and_tiebreak() {
    for (n, seed) in [(1usize, 11u64), (17, 12), (64, 13)] {
        let a = rand_mat(n, n, seed);
        let f = lu(&a, n).expect("random square is a.s. nonsingular");
        // Permutation validity: perm is a bijection.
        let mut seen = vec![false; n];
        for &p in f.perm() {
            assert!(!seen[p], "duplicate perm entry");
            seen[p] = true;
        }
        // Solve residual: ‖A·x − b‖ small.
        let x_true: Vec<f64> = (0..n).map(|i| 1.0 + (i as f64) * 0.25).collect();
        let mut b = vec![0.0; n];
        for i in 0..n {
            let mut acc = 0.0;
            for j in 0..n {
                acc = a[i * n + j].mul_add(x_true[j], acc);
            }
            b[i] = acc;
        }
        f.solve(&mut b);
        for i in 0..n {
            assert!((b[i] - x_true[i]).abs() < 1e-9, "lu solve at {i} (n={n})");
        }
    }
    // Deterministic tie-break: column of equal magnitudes must pick the
    // LOWEST row index (P2). Matrix [[1,2],[−1,4]]: |a00| == |a10| → no swap.
    let f = lu(&[1.0, 2.0, -1.0, 4.0], 2).unwrap();
    assert_eq!(f.perm(), &[0, 1], "tie must keep the lowest index");
    // Singularity is typed, with the failing step.
    match lu(&[1.0, 2.0, 2.0, 4.0], 2) {
        Err(FactorError::Singular { index }) => assert_eq!(index, 1),
        other => panic!("expected Singular, got {other:?}"),
    }
    // Kahan matrix factors with bounded growth (ledgered statistic).
    let k = kahan(24, 1.2);
    let fk = lu(&k, 24).unwrap();
    assert!(fk.growth.is_finite() && fk.growth >= 1.0);
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"lu\",\"verdict\":\"pass\",\"detail\":\"residual+perm+tiebreak; kahan growth {:.2}\"}}",
        fk.growth
    );
}

#[test]
fn qr_orthogonality_reconstruction_least_squares() {
    for (m, n, seed) in [(1usize, 1usize, 21u64), (9, 4, 22), (60, 24, 23), (65, 33, 24)] {
        let a = rand_mat(m, n, seed);
        let f = qr(&a, m, n);
        // Orthogonality: ‖Qᵀ·Q − I‖ via applying Qᵀ then Q to unit vectors.
        for j in 0..n.min(8) {
            let mut e = vec![0.0; m];
            e[j] = 1.0;
            f.apply_qt(&mut e);
            f.apply_q(&mut e);
            for (i, &v) in e.iter().enumerate() {
                let want = if i == j { 1.0 } else { 0.0 };
                assert!((v - want).abs() < 1e-13, "Q not orthogonal at ({i},{j})");
            }
        }
        // Reconstruction: column j of A equals Q·[R eⱼ; 0].
        for j in 0..n {
            let mut col = vec![0.0; m];
            for i in 0..=j {
                col[i] = f.r(i, j);
            }
            f.apply_q(&mut col);
            for (i, &cv) in col.iter().enumerate() {
                let want = a[i * n + j];
                assert!(
                    (cv - want).abs() < 1e-12 * (1.0 + want.abs()),
                    "A != QR at ({i},{j}) for {m}x{n}"
                );
            }
        }
    }
    // Least squares vs normal equations on a well-conditioned system.
    let (m, n) = (40usize, 7usize);
    let a = rand_mat(m, n, 25);
    let b: Vec<f64> = (0..m).map(|i| (i as f64) * 0.1 - 2.0).collect();
    let x_qr = qr(&a, m, n).solve_ls(&b);
    // Normal equations: (AᵀA)x = Aᵀb via Cholesky.
    let at: Vec<f64> = (0..n * m).map(|i| a[(i % m) * n + i / m]).collect();
    let mut ata = vec![0.0; n * n];
    gemm_f64(n, n, m, 1.0, &at, &a, 0.0, &mut ata);
    let mut atb = vec![0.0; n];
    for i in 0..n {
        let mut acc = 0.0;
        for k in 0..m {
            acc = at[i * m + k].mul_add(b[k], acc);
        }
        atb[i] = acc;
    }
    let mut x_ne = atb;
    cholesky(&ata, n).unwrap().solve(&mut x_ne);
    for i in 0..n {
        assert!(
            (x_qr[i] - x_ne[i]).abs() < 1e-9,
            "LS mismatch at {i}: {} vs {}",
            x_qr[i],
            x_ne[i]
        );
    }
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"qr\",\"verdict\":\"pass\",\"detail\":\"orthogonality+reconstruct+LS, 4 shapes\"}}"
    );
}

#[test]
fn tsqr_matches_direct_qr_and_is_deterministic() {
    let (m, n) = (300usize, 6usize);
    let a = rand_mat(m, n, 31);
    let direct = {
        let f = qr(&a, m, n);
        let mut r = vec![0.0; n * n];
        for i in 0..n {
            let flip = if f.r(i, i) < 0.0 { -1.0 } else { 1.0 };
            for j in i..n {
                r[i * n + j] = flip * f.r(i, j);
            }
        }
        r
    };
    for block in [8usize, 50, 300] {
        let r = tsqr_r(&a, m, n, block.max(n));
        for i in 0..n * n {
            assert!(
                (r[i] - direct[i]).abs() < 1e-10 * (1.0 + direct[i].abs()),
                "TSQR(block={block}) vs direct at {i}: {} vs {}",
                r[i],
                direct[i]
            );
        }
    }
    // Bit-determinism for a fixed tree: rerun equality.
    let r1 = tsqr_r(&a, m, n, 50);
    let r2 = tsqr_r(&a, m, n, 50);
    assert!(r1.iter().zip(&r2).all(|(x, y)| x.to_bits() == y.to_bits()));
    // RᵀR == AᵀA (the Gram identity) within tolerance.
    let at: Vec<f64> = (0..n * m).map(|i| a[(i % m) * n + i / m]).collect();
    let mut ata = vec![0.0; n * n];
    gemm_f64(n, n, m, 1.0, &at, &a, 0.0, &mut ata);
    let r = tsqr_r(&a, m, n, 50);
    let rt: Vec<f64> = (0..n * n).map(|i| r[(i % n) * n + i / n]).collect();
    let mut rtr = vec![0.0; n * n];
    gemm_f64(n, n, n, 1.0, &rt, &r, 0.0, &mut rtr);
    for i in 0..n * n {
        assert!((rtr[i] - ata[i]).abs() < 1e-10 * (1.0 + ata[i].abs()), "Gram identity at {i}");
    }
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"tsqr\",\"verdict\":\"pass\",\"detail\":\"3 tree shapes agree with direct QR; fixed tree bitwise stable\"}}"
    );
}

#[test]
fn svd_reconstruction_orthogonality_and_hilbert() {
    // Known case: diagonal matrix → σ = |diag| sorted descending.
    let d = [3.0, -7.0, 0.5, 2.0];
    let n = 4;
    let mut a = vec![0.0; n * n];
    for (i, &v) in d.iter().enumerate() {
        a[i * n + i] = v;
    }
    let s = svd_jacobi(&a, n, n);
    let want = [7.0, 3.0, 2.0, 0.5];
    for (got, w) in s.sigma.iter().zip(&want) {
        assert!((got - w).abs() < 1e-12, "diag sigma {got} vs {w}");
    }
    // Random reconstruction + orthogonality.
    let (m, nn) = (30usize, 12usize);
    let a = rand_mat(m, nn, 41);
    let s = svd_jacobi(&a, m, nn);
    // ‖A − U·Σ·Vᵀ‖/‖A‖.
    let mut usv = vec![0.0; m * nn];
    for i in 0..m {
        for j in 0..nn {
            let mut acc = 0.0f64;
            for k in 0..nn {
                acc = (s.u[i * nn + k] * s.sigma[k]).mul_add(s.v[j * nn + k], acc);
            }
            usv[i * nn + j] = acc;
        }
    }
    let diff: Vec<f64> = a.iter().zip(&usv).map(|(x, y)| x - y).collect();
    assert!(frob(&diff) / frob(&a) < 1e-13, "SVD reconstruction");
    // UᵀU = I, VᵀV = I.
    for p in 0..nn {
        for q in 0..nn {
            let (mut du, mut dv) = (0.0f64, 0.0f64);
            for i in 0..m {
                du = s.u[i * nn + p].mul_add(s.u[i * nn + q], du);
            }
            for i in 0..nn {
                dv = s.v[i * nn + p].mul_add(s.v[i * nn + q], dv);
            }
            let want = if p == q { 1.0 } else { 0.0 };
            assert!((du - want).abs() < 1e-13, "U orthogonality ({p},{q})");
            assert!((dv - want).abs() < 1e-13, "V orthogonality ({p},{q})");
        }
    }
    // Hilbert 8: spectral condition ≈ 1.53e10 — the small-singular-value
    // relative-accuracy claim in action.
    let h = hilbert(8);
    let sh = svd_jacobi(&h, 8, 8);
    let cond = sh.sigma[0] / sh.sigma[7];
    assert!(
        (1.0e10..3.0e10).contains(&cond),
        "Hilbert-8 spectral condition {cond:.3e} outside the known ~1.5e10 band"
    );
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"svd\",\"verdict\":\"pass\",\"detail\":\"reconstruct+orthogonality; hilbert8 cond {cond:.3e}\"}}"
    );
}

#[test]
fn condition_estimate_hilbert() {
    let n = 8;
    let h = hilbert(n);
    let f = lu(&h, n).unwrap();
    let est = f.condition_1(&h);
    // κ₁(H₈) ≈ 3.39e10; the Hager estimator is a lower bound within a
    // small factor. Accept a broad but decisive band.
    assert!(
        (1.0e9..1.0e12).contains(&est),
        "Hilbert-8 kappa_1 estimate {est:.3e} outside plausibility band"
    );
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"cond\",\"verdict\":\"pass\",\"detail\":\"hilbert8 kappa1 est {est:.3e} (true ~3.4e10)\"}}"
    );
}

#[test]
fn degenerate_shapes() {
    // n = 0: all factorizations accept and do nothing.
    assert!(cholesky(&[], 0).is_ok());
    assert!(lu(&[], 0).is_ok());
    let f = qr(&[], 0, 0);
    assert_eq!(f.r(0, 0).to_bits(), 0.0f64.to_bits());
    // 1×1.
    let c = cholesky(&[9.0], 1).unwrap();
    assert!((c.l(0, 0) - 3.0).abs() < 1e-15);
    // Single column QR: R(0,0) = ±‖a‖.
    let a = [3.0, 4.0];
    let f = qr(&a, 2, 1);
    assert!((f.r(0, 0).abs() - 5.0).abs() < 1e-14);
    let x = f.solve_ls(&[3.0, 4.0]);
    assert!((x[0] - 1.0).abs() < 1e-14);
}

/// Recorded on aarch64-apple (M4 Pro); must match on x86-64 (trj).
const GOLDEN_HASH: u64 = 0x0; // placeholder: set from first run

#[test]
fn factorization_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for b in v.to_bits().to_le_bytes() {
            acc ^= u64::from(b);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let n = 48;
    let spd = rand_spd(n, 0x77);
    let ch = cholesky(&spd, n).unwrap();
    for i in 0..n {
        for j in 0..=i {
            feed(ch.l(i, j));
        }
    }
    let g = rand_mat(n, n, 0x78);
    let f = lu(&g, n).unwrap();
    let mut b: Vec<f64> = (0..n).map(|i| (i as f64) * 0.5 - 3.0).collect();
    f.solve(&mut b);
    for &v in &b {
        feed(v);
    }
    let a = rand_mat(120, 9, 0x79);
    for &v in &tsqr_r(&a, 120, 9, 40) {
        feed(v);
    }
    let s = svd_jacobi(&rand_mat(20, 8, 0x7A), 20, 8);
    for &v in &s.sigma {
        feed(v);
    }
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"factor-golden\",\"verdict\":\"info\",\"detail\":\"{acc:#018x}\"}}"
    );
    assert_eq!(
        acc, GOLDEN_HASH,
        "factorization bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with \
         semantic justification (golden-evidence policy)"
    );
}
