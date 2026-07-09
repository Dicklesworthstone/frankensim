//! Sparse-GP battery (tzeh lane c): EXACTNESS RECOVERY at Z = X
//! (predictions and ELBO match the exact GP — the tight-bound
//! identity); the Titsias ELBO lower-bounds the exact LML for every
//! m < n (the variational guarantee, checked instance-wise); the
//! accuracy-vs-m ladder against exact predictions; golden.

use fs_bo::{Gp, Kernel, Matern, SparseGp, farthest_point_inducing};
use fs_rand::StreamKey;

fn log(case: &str, verdict: &str, detail: &str) {
    println!("{{\"suite\":\"fs-bo-sparse\",\"case\":\"{case}\",\"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}");
}

fn rand_pts(n: usize, d: usize, tile: u32) -> Vec<Vec<f64>> {
    let mut s = StreamKey { seed: 141, kernel: 0x0561, tile }.stream();
    (0..n).map(|_| (0..d).map(|_| s.next_f64()).collect()).collect()
}

fn target(x: &[f64]) -> f64 {
    fs_math::det::sin(4.0 * x[0]) + 0.5 * fs_math::det::cos(6.0 * x[1]) + 0.3 * x[0] * x[1]
}

fn kernel() -> Kernel {
    Kernel {
        family: Matern::FiveHalves,
        signal: 1.0,
        lengthscales: vec![0.3, 0.3],
    }
}

#[test]
fn exactness_recovery_at_z_equals_x() {
    let x = rand_pts(40, 2, 1);
    let y: Vec<f64> = x.iter().map(|p| target(p)).collect();
    let noise = 1e-4;
    let exact = Gp::fit(&x, &y, kernel(), noise);
    let sparse = SparseGp::fit(&x, &y, kernel(), noise, x.clone());
    // Predictions match at held-out probes.
    let probes = rand_pts(20, 2, 2);
    let mut worst_mean = 0.0f64;
    let mut worst_var = 0.0f64;
    for p in &probes {
        let (me, ve) = exact.predict(p);
        let (ms, vs) = sparse.predict(p);
        worst_mean = worst_mean.max((me - ms).abs());
        worst_var = worst_var.max((ve - vs).abs());
    }
    assert!(worst_mean < 1e-5, "mean mismatch at Z=X: {worst_mean:.2e}");
    assert!(worst_var < 1e-5, "variance mismatch at Z=X: {worst_var:.2e}");
    // ELBO tight: equals the exact LML.
    let gap = (sparse.elbo - exact.lml).abs();
    assert!(gap < 1e-6, "ELBO not tight at Z=X: gap {gap:.2e}");
    log(
        "exactness",
        "pass",
        &format!("mean {worst_mean:.1e}, var {worst_var:.1e}, ELBO gap {gap:.1e}"),
    );
}

#[test]
fn elbo_lower_bounds_exact_lml() {
    let x = rand_pts(120, 2, 3);
    let y: Vec<f64> = x.iter().map(|p| target(p)).collect();
    let noise = 1e-3;
    let exact = Gp::fit(&x, &y, kernel(), noise);
    let mut prev_elbo = f64::NEG_INFINITY;
    let mut line = String::new();
    for &m in &[5usize, 15, 40, 120] {
        let z = if m == 120 {
            x.clone()
        } else {
            farthest_point_inducing(&x, m)
        };
        let sparse = SparseGp::fit(&x, &y, kernel(), noise, z);
        assert!(
            sparse.elbo <= exact.lml + 1e-6,
            "ELBO must lower-bound the exact LML: {} vs {} at m={m}",
            sparse.elbo,
            exact.lml
        );
        assert!(
            sparse.elbo >= prev_elbo - 1e-6,
            "ELBO should improve with more inducing points (farthest-point nesting)"
        );
        prev_elbo = sparse.elbo;
        line.push_str(&format!("m={m}: {:.2} ", sparse.elbo));
    }
    log(
        "elbo-bound",
        "pass",
        &format!("{line}<= LML {:.2}", exact.lml),
    );
}

#[test]
fn accuracy_ladder_vs_exact() {
    let n = 300usize;
    let x = rand_pts(n, 2, 4);
    let y: Vec<f64> = x.iter().map(|p| target(p)).collect();
    let noise = 1e-3;
    let exact = Gp::fit(&x, &y, kernel(), noise);
    let probes = rand_pts(60, 2, 5);
    let exact_means: Vec<f64> = probes.iter().map(|p| exact.predict(p).0).collect();
    let mut prev_rmse = f64::INFINITY;
    let mut line = String::new();
    for &m in &[10usize, 30, 80] {
        let z = farthest_point_inducing(&x, m);
        let sparse = SparseGp::fit(&x, &y, kernel(), noise, z);
        let rmse: f64 = fs_math::det::sqrt(
            probes
                .iter()
                .zip(&exact_means)
                .map(|(p, em)| {
                    let d = sparse.predict(p).0 - em;
                    d * d
                })
                .sum::<f64>()
                / probes.len() as f64,
        );
        assert!(
            rmse < prev_rmse + 1e-9,
            "RMSE-vs-exact should not increase with m: {line} then m={m}: {rmse:.4}"
        );
        prev_rmse = rmse;
        line.push_str(&format!("m={m}: {rmse:.4} "));
    }
    assert!(
        prev_rmse < 0.02,
        "m=80 should track the exact GP closely: {prev_rmse:.4}"
    );
    log("accuracy-ladder", "pass", line.trim());
}

const GOLDEN_HASH: u64 = 0; // recorded on first run, then frozen

#[test]
fn sparse_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let x = rand_pts(60, 2, 6);
    let y: Vec<f64> = x.iter().map(|p| target(p)).collect();
    let z = farthest_point_inducing(&x, 12);
    for zi in z.iter().step_by(3) {
        for v in zi {
            feed(*v);
        }
    }
    let sparse = SparseGp::fit(&x, &y, kernel(), 1e-3, z);
    feed(sparse.elbo);
    for p in rand_pts(8, 2, 7) {
        let (m, v) = sparse.predict(&p);
        feed(m);
        feed(v);
    }
    log("sparse-golden", "info", &format!("{acc:#018x}"));
    assert_eq!(
        acc, GOLDEN_HASH,
        "sparse bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}
