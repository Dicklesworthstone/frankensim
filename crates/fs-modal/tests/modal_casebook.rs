//! E2E casebook for the vibration eigenproblem facility (bead
//! frankensim-fsim-vibration-eig-jw6yq): slice windows of a ≥100k-DoF
//! finite-difference plate pencil with inertia-certified counts checked
//! against the ANALYTIC spectrum, emit JSON-line evidence per window, and
//! cross-check the dense strategy against the FrankenScipy sibling
//! eigensolver (dev-dependency oracle, per workspace convention).

use fs_modal::{SliceOptions, eigh_gen_dense, slice_window};
use fs_sparse::{Coo, Csr};

fn laplacian_2d(s: usize) -> Csr {
    let n = s * s;
    let mut coo = Coo::new(n, n);
    for i in 0..s {
        for j in 0..s {
            let u = i * s + j;
            coo.push(u, u, 4.0);
            if i > 0 {
                coo.push(u, u - s, -1.0);
            }
            if i + 1 < s {
                coo.push(u, u + s, -1.0);
            }
            if j > 0 {
                coo.push(u, u - 1, -1.0);
            }
            if j + 1 < s {
                coo.push(u, u + 1, -1.0);
            }
        }
    }
    coo.assemble()
}

fn identity(n: usize) -> Csr {
    let mut coo = Coo::new(n, n);
    for i in 0..n {
        coo.push(i, i, 1.0);
    }
    coo.assemble()
}

/// Smallest `count` analytic grid-Laplacian eigenvalues, ascending.
fn analytic_smallest(s: usize, count: usize) -> Vec<f64> {
    let mut vals: Vec<f64> = Vec::new();
    // The smallest eigenvalues come from small (i, j); scanning a corner
    // block is exact as long as the block's worst value exceeds the cutoff.
    let block = 40.min(s);
    for i in 1..=block {
        for j in 1..=block {
            let li = 2.0 - 2.0 * (i as f64 * std::f64::consts::PI / (s as f64 + 1.0)).cos();
            let lj = 2.0 - 2.0 * (j as f64 * std::f64::consts::PI / (s as f64 + 1.0)).cos();
            vals.push(li + lj);
        }
    }
    vals.sort_by(f64::total_cmp);
    vals.truncate(count);
    vals
}

#[test]
fn hundred_thousand_dof_windows_certified_against_analytic_counts() {
    let s = 317; // 100_489 unknowns ≥ the bead's 100k-DoF gate
    let n = s * s;
    let k = laplacian_2d(s);
    let m = identity(n);
    let truth = analytic_smallest(s, 40);

    // Window endpoints must split DISTINCT consecutive analytic values —
    // the spectrum is full of degenerate (i,j)/(j,i) pairs, and a midpoint
    // of a tied pair IS an eigenvalue (the endpoint factorization would
    // correctly refuse as singular). Walk forward to the next strict gap.
    let split_after = |idx: usize| -> f64 {
        let mut i = idx;
        while truth[i + 1] <= truth[i] {
            i += 1;
        }
        f64::midpoint(truth[i], truth[i + 1])
    };
    // Two windows over the low spectrum; both span DEGENERATE pairs,
    // exercising the deflation-restart path at scale.
    let windows = [
        (split_after(2), split_after(7)),
        (split_after(11), split_after(17)),
    ];
    for (widx, &(lo, hi)) in windows.iter().enumerate() {
        let analytic_count = truth.iter().filter(|&&v| v > lo && v <= hi).count();
        let t0 = std::time::Instant::now();
        let rep = slice_window(&k, &m, (lo, hi), &SliceOptions::default()).expect("slice");
        let wall_ms = t0.elapsed().as_secs_f64() * 1e3;
        assert_eq!(
            rep.expected, analytic_count,
            "inertia count must equal the analytic count in window {widx}"
        );
        assert_eq!(rep.modes.len(), rep.expected);
        let mut truth_in: Vec<f64> = truth
            .iter()
            .copied()
            .filter(|&v| v > lo && v <= hi)
            .collect();
        truth_in.sort_by(f64::total_cmp);
        for (mode, want) in rep.modes.iter().zip(&truth_in) {
            assert!(
                mode.interval.0 <= *want && *want <= mode.interval.1,
                "window {widx}: certified interval [{}, {}] must contain analytic {want}",
                mode.interval.0,
                mode.interval.1
            );
        }
        println!(
            "{{\"suite\":\"fs-modal-casebook\",\"stage\":\"window\",\"n\":{n},\"window\":[{lo},{hi}],\"expected\":{},\"below_low\":{},\"below_high\":{},\"shift\":{},\"factorizations\":{},\"lanczos_iters\":{},\"restarts\":{},\"factor_nnz_l\":{},\"factor_peak_bytes\":{},\"pivots_delayed\":{},\"max_residual\":{:.3e},\"wall_ms\":{wall_ms:.0}}}",
            rep.expected,
            rep.below_low,
            rep.below_high,
            rep.stats.shift,
            rep.stats.factorizations,
            rep.stats.lanczos_iters,
            rep.stats.restarts,
            rep.stats.factor_nnz_l,
            rep.stats.factor_peak_bytes,
            rep.stats.pivots_delayed,
            rep.modes.iter().fold(0.0f64, |mx, md| mx.max(md.residual)),
        );
    }
}

#[test]
fn dense_strategy_agrees_with_frankenscipy_oracle() {
    // Independent-oracle stage (workspace convention): the fs-modal dense
    // strategy on (C, I) must agree with the FrankenScipy sibling's
    // symmetric eigensolver on C. C is a deterministic seeded dense
    // symmetric matrix.
    let n = 24;
    let mut seed = 0x0DA1_5EEDu64;
    let mut lcg = || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
    };
    let mut c = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..=i {
            let v = lcg();
            c[i * n + j] = v;
            c[j * n + i] = v;
        }
        c[i * n + i] += 4.0;
    }
    let mut ident = vec![0.0f64; n * n];
    for i in 0..n {
        ident[i * n + i] = 1.0;
    }
    let mine = eigh_gen_dense(&c, &ident, n).expect("dense path");
    let rows: Vec<Vec<f64>> = (0..n).map(|i| c[i * n..(i + 1) * n].to_vec()).collect();
    let oracle =
        fsci_linalg::eigvalsh(&rows, fsci_linalg::DecompOptions::default()).expect("fsci eigvalsh");
    assert_eq!(mine.len(), oracle.len());
    let mut worst = 0.0f64;
    for (mode, want) in mine.iter().zip(&oracle) {
        let delta = (mode.lambda - want).abs();
        worst = worst.max(delta);
        assert!(
            delta <= 1e-9,
            "fs-modal {} vs fsci oracle {want} (delta {delta})",
            mode.lambda
        );
    }
    println!(
        "{{\"suite\":\"fs-modal-casebook\",\"stage\":\"fsci-oracle\",\"n\":{n},\"worst_delta\":{worst:.3e},\"verdict\":\"pass\"}}"
    );
}
