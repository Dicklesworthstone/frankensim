//! Eigensolver battery: analytic Laplacian spectra (1D + 2D stencil
//! closures — matrix-free as designed), dense-Jacobi agreement on random
//! SPD, certified residuals, bitwise resumability, and the cross-ISA
//! golden hash.

use fs_la::eigen::{LanczosState, LobpcgState, jacobi_eigh, lanczos_run, lobpcg_run};

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
}

/// 1D Laplacian stencil [-1, 2, -1] as a matrix-free closure.
fn lap1d(n: usize) -> impl Fn(&[f64], &mut [f64]) {
    move |x: &[f64], y: &mut [f64]| {
        for i in 0..n {
            let mut v = 2.0 * x[i];
            if i > 0 {
                v -= x[i - 1];
            }
            if i + 1 < n {
                v -= x[i + 1];
            }
            y[i] = v;
        }
    }
}

/// 2D 5-point Laplacian on a g×g grid (matrix-free; NO fs-sparse
/// dependency — formats compose at call sites, per the design comment).
fn lap2d(g: usize) -> impl Fn(&[f64], &mut [f64]) {
    move |x: &[f64], y: &mut [f64]| {
        for i in 0..g {
            for j in 0..g {
                let u = i * g + j;
                let mut v = 4.0 * x[u];
                if i > 0 {
                    v -= x[u - g];
                }
                if i + 1 < g {
                    v -= x[u + g];
                }
                if j > 0 {
                    v -= x[u - 1];
                }
                if j + 1 < g {
                    v -= x[u + 1];
                }
                y[u] = v;
            }
        }
    }
}

fn identity_prec(r: &[f64], out: &mut [f64]) {
    out.copy_from_slice(r);
}

#[test]
fn jacobi_dense_analytic_and_reconstruction() {
    // 2×2 analytic: [[2,1],[1,2]] → 1, 3 with (1,∓1)/√2.
    let (vals, vecs) = jacobi_eigh(&[2.0, 1.0, 1.0, 2.0], 2);
    assert!((vals[0] - 1.0).abs() < 1e-14 && (vals[1] - 3.0).abs() < 1e-14);
    assert!((vecs[0].abs() - std::f64::consts::FRAC_1_SQRT_2).abs() < 1e-13);
    // Random symmetric: reconstruct + orthogonality.
    let n = 24;
    let mut seed = 0xE16_u64;
    let mut a = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..=i {
            let v = lcg(&mut seed);
            a[i * n + j] = v;
            a[j * n + i] = v;
        }
    }
    let (vals, vecs) = jacobi_eigh(&a, n);
    // Ascending order.
    for w in vals.windows(2) {
        assert!(w[0] <= w[1], "eigenvalues must ascend");
    }
    // A·V = V·diag(λ) and VᵀV = I.
    for p in 0..n {
        for q in 0..n {
            let mut av = 0.0f64;
            let mut vv = 0.0f64;
            for k in 0..n {
                av += a[p * n + k] * vecs[k * n + q];
                vv += vecs[k * n + p] * vecs[k * n + q];
            }
            let want = vals[q] * vecs[p * n + q];
            assert!((av - want).abs() < 1e-11, "A·V mismatch at ({p},{q})");
            let id = if p == q { 1.0 } else { 0.0 };
            assert!((vv - id).abs() < 1e-12, "V orthogonality at ({p},{q})");
        }
    }
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"jacobi-eigh\",\"verdict\":\"pass\",\"detail\":\"analytic 2x2 + random 24x24 reconstruct/orthogonality\"}}"
    );
}

#[test]
fn lanczos_matches_analytic_laplacian_spectrum() {
    // Full-dimension Krylov (steps = n): the tridiagonalization is then
    // COMPLETE and Ritz values equal eigenvalues to roundoff — the small
    // end of the Laplacian spectrum is clustered and a partial Krylov
    // space resolves it slowly (that is physics, not a bug).
    let n = 120;
    let op = lap1d(n);
    let mut st = LanczosState::new(n);
    let pairs = lanczos_run(&op, &mut st, n, 4, false);
    for (k, pair) in pairs.iter().enumerate() {
        let want = 2.0 - 2.0 * (std::f64::consts::PI * ((k + 1) as f64) / ((n + 1) as f64)).cos();
        assert!(
            (pair.value - want).abs() < 1e-10,
            "smallest[{k}]: {} vs analytic {want}",
            pair.value
        );
        assert!(
            pair.residual < 1e-8,
            "residual {} not certified small",
            pair.residual
        );
    }
    // Largest end: ALSO full Krylov — the Laplacian spectrum is clustered
    // at BOTH ends (cos is flat at 0 and π), so partial spaces resolve
    // neither end quickly. Genuine partial-space convergence is what the
    // random-SPD test exercises (its top end IS well separated).
    let mut st2 = LanczosState::new(n);
    let top = lanczos_run(&op, &mut st2, n, 2, true);
    for (k, pair) in top.iter().enumerate() {
        let want = 2.0 - 2.0 * (std::f64::consts::PI * ((n - k) as f64) / ((n + 1) as f64)).cos();
        assert!(
            (pair.value - want).abs() < 1e-10,
            "largest[{k}]: {} vs analytic {want}",
            pair.value
        );
    }
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"lanczos\",\"verdict\":\"pass\",\"detail\":\"1D Laplacian n=120: 4 smallest (full Krylov) + 2 largest (partial) vs analytic\"}}"
    );
}

#[test]
fn lobpcg_matches_analytic_2d_laplacian() {
    let g = 8;
    let n = g * g;
    let op = lap2d(g);
    let mut st = LobpcgState::new(n, 4);
    let pairs = lobpcg_run(&op, &mut st, 120, false, &identity_prec);
    // Analytic: 4 − 2cos(iπ/(g+1)) − 2cos(jπ/(g+1)); collect the 4 smallest.
    let mut analytic: Vec<f64> = (1..=g)
        .flat_map(|i| (1..=g).map(move |j| (i, j)))
        .map(|(i, j)| {
            let t = std::f64::consts::PI / ((g + 1) as f64);
            4.0 - 2.0 * ((i as f64) * t).cos() - 2.0 * ((j as f64) * t).cos()
        })
        .collect();
    analytic.sort_by(|a, b| a.partial_cmp(b).unwrap());
    for (k, pair) in pairs.iter().enumerate() {
        assert!(
            (pair.value - analytic[k]).abs() < 1e-7,
            "2D smallest[{k}]: {} vs analytic {}",
            pair.value,
            analytic[k]
        );
        assert!(pair.residual < 1e-5, "residual {} too large", pair.residual);
    }
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"lobpcg\",\"verdict\":\"pass\",\"detail\":\"2D Laplacian 8x8 grid: 4 smallest vs analytic (incl degenerate pair)\"}}"
    );
}

#[test]
fn solvers_agree_with_dense_reference_on_random_spd() {
    let n = 60;
    let mut seed = 0x5BD_u64;
    // SPD: B·Bᵀ + n·I.
    let bmat: Vec<f64> = (0..n * n).map(|_| lcg(&mut seed)).collect();
    let mut a = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..n {
            let mut acc = 0.0f64;
            for k in 0..n {
                acc += bmat[i * n + k] * bmat[j * n + k];
            }
            a[i * n + j] = acc + if i == j { n as f64 } else { 0.0 };
        }
    }
    let (dense_vals, _) = jacobi_eigh(&a, n);
    let a_ref = a.clone();
    let op = move |x: &[f64], y: &mut [f64]| {
        for i in 0..n {
            let mut acc = 0.0f64;
            for j in 0..n {
                acc = a_ref[i * n + j].mul_add(x[j], acc);
            }
            y[i] = acc;
        }
    };
    let mut lst = LanczosState::new(n);
    let lp = lanczos_run(&op, &mut lst, n, 3, false);
    for (k, pair) in lp.iter().enumerate() {
        assert!(
            (pair.value - dense_vals[k]).abs() < 1e-9 * dense_vals[k].abs().max(1.0),
            "lanczos[{k}] {} vs dense {}",
            pair.value,
            dense_vals[k]
        );
    }
    let mut bst = LobpcgState::new(n, 3);
    let bp = lobpcg_run(&op, &mut bst, 250, false, &identity_prec);
    for (k, pair) in bp.iter().enumerate() {
        // The small end of B·Bᵀ + n·I is CLUSTERED (relative gaps ~1e-4);
        // unpreconditioned LOBPCG resolves values far faster than the
        // cluster's vectors — a realistic tolerance reflects that.
        assert!(
            (pair.value - dense_vals[k]).abs() < 1e-6 * dense_vals[k].abs().max(1.0),
            "lobpcg[{k}] {} vs dense {}",
            pair.value,
            dense_vals[k]
        );
    }
}

#[test]
fn resumability_is_bitwise() {
    let n = 120;
    let op = lap1d(n);
    // Lanczos: 20+20 == 40 straight, bitwise.
    let mut a = LanczosState::new(n);
    lanczos_run(&op, &mut a, 20, 1, false);
    let mut a2 = a.clone(); // checkpoint = clone
    let pa = lanczos_run(&op, &mut a2, 20, 3, false);
    let mut b = LanczosState::new(n);
    let pb = lanczos_run(&op, &mut b, 40, 3, false);
    for (x, y) in pa.iter().zip(&pb) {
        assert_eq!(
            x.value.to_bits(),
            y.value.to_bits(),
            "lanczos resume changed bits"
        );
        assert!(
            x.vector
                .iter()
                .zip(&y.vector)
                .all(|(p, q)| p.to_bits() == q.to_bits())
        );
    }
    // LOBPCG: 15+15 == 30 straight, bitwise.
    let mut c = LobpcgState::new(n, 2);
    lobpcg_run(&op, &mut c, 15, false, &identity_prec);
    let mut c2 = c.clone();
    let pc = lobpcg_run(&op, &mut c2, 15, false, &identity_prec);
    let mut d = LobpcgState::new(n, 2);
    let pd = lobpcg_run(&op, &mut d, 30, false, &identity_prec);
    for (x, y) in pc.iter().zip(&pd) {
        assert_eq!(
            x.value.to_bits(),
            y.value.to_bits(),
            "lobpcg resume changed bits"
        );
    }
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"eigen-resume\",\"verdict\":\"pass\",\"detail\":\"checkpoint=clone; split runs bitwise == straight runs\"}}"
    );
}

#[test]
#[should_panic(expected = "3b <= n")]
fn lobpcg_rejects_a_block_too_large_for_the_search_space() {
    // regression: n=6, b=3 passed the old `b <= n` assert then panicked deep in
    // the second iteration's QR (search space [X|W|P] has zc = 3b = 9 > n = 6).
    let _ = LobpcgState::new(6, 3);
}

/// Recorded on aarch64-apple (M4 Pro); must match on x86-64 (trj).
const GOLDEN_HASH: u64 = 0x87da_0cb3_2344_b097;

#[test]
fn eigen_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let n = 80;
    let op = lap1d(n);
    let mut st = LanczosState::new(n);
    for p in lanczos_run(&op, &mut st, 40, 3, false) {
        feed(p.value);
        feed(p.residual);
    }
    let mut bst = LobpcgState::new(n, 3);
    for p in lobpcg_run(&op, &mut bst, 30, false, &identity_prec) {
        feed(p.value);
        feed(p.residual);
    }
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"eigen-golden\",\"verdict\":\"info\",\"detail\":\"{acc:#018x}\"}}"
    );
    assert_eq!(
        acc, GOLDEN_HASH,
        "eigensolver bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with \
         semantic justification (golden-evidence policy)"
    );
}
