//! fs-fmm conformance battery (bead tfz.20).
//!
//! - fmm-001: accuracy vs interpolation order against the direct
//!   oracle — the error curve must fall (near-exponentially in p) and
//!   is ledgered.
//! - fmm-002 G3: translation invariance — a rigidly shifted cloud
//!   produces the same potentials to tight tolerance.
//! - fmm-003: scaling trend — measured time vs N fitted exponent well
//!   below the direct method's 2 (the 10⁷-point wall-clock target is
//!   the perf lanes' scope, ledgered here as a trend).

use fs_fmm::{Fmm, Laplace3d};
use std::fmt::Write as _;
use std::time::Instant;

fn verdict(name: &str, pass: bool, details: &str) {
    println!(
        "{{\"test\":\"{name}\",\"verdict\":\"{}\",{details}}}",
        if pass { "pass" } else { "fail" }
    );
    assert!(pass, "{name} failed: {details}");
}

struct Lcg(u64);

impl Lcg {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0
    }
    #[allow(clippy::cast_precision_loss)]
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
}

fn cloud(n: usize, seed: u64) -> (Vec<[f64; 3]>, Vec<f64>) {
    let mut lcg = Lcg(seed);
    let pts = (0..n)
        .map(|_| [lcg.unit(), lcg.unit(), lcg.unit()])
        .collect();
    let q = (0..n).map(|_| lcg.unit() - 0.5).collect();
    (pts, q)
}

// ------------------------------------------------------------------ fmm-001

#[test]
fn fmm_001_accuracy_vs_order() {
    let (pts, q) = cloud(1500, 0x1001_2026_0708_0001);
    let kernel = Laplace3d;
    let oracle = Fmm::new(&kernel, pts.clone(), 2, 32)
        .expect("admitted fixture")
        .direct(&q)
        .expect("finite fixture");
    let scale = oracle.iter().map(|v| v * v).sum::<f64>().sqrt();
    let mut errs = Vec::new();
    let mut rows = String::new();
    for p in [3usize, 5, 7] {
        let fmm = Fmm::new(&kernel, pts.clone(), p, 32).expect("admitted order sweep");
        let got = fmm.potentials(&q).expect("finite fixture");
        let err = got
            .iter()
            .zip(&oracle)
            .map(|(a, b)| (a - b) * (a - b))
            .sum::<f64>()
            .sqrt()
            / scale;
        let _ = write!(
            rows,
            "{{\"order\":{p},\"rel_l2\":{err:.3e},\"tree\":{}}},",
            fmm.stats()
        );
        errs.push(err);
    }
    let monotone = errs[1] < errs[0] && errs[2] < errs[1];
    let pass = monotone && errs[2] < 1e-5 && errs[0] < 1e-1;
    verdict(
        "fmm-001",
        pass,
        &format!(
            "\"detail\":\"Chebyshev order sweep vs direct oracle, 1500 pts\",\
             \"rows\":[{}]",
            rows.trim_end_matches(',')
        ),
    );
}

// ------------------------------------------------------------------ fmm-002

#[test]
fn fmm_002_translation_invariance() {
    let (pts, q) = cloud(1200, 0x1001_2026_0708_0002);
    let kernel = Laplace3d;
    let base = Fmm::new(&kernel, pts.clone(), 6, 32)
        .expect("admitted fixture")
        .potentials(&q)
        .expect("finite fixture");
    let shift = [17.25, -4.5, 9.75]; // dyadic-friendly rigid shift
    let moved: Vec<[f64; 3]> = pts
        .iter()
        .map(|p| [p[0] + shift[0], p[1] + shift[1], p[2] + shift[2]])
        .collect();
    let shifted = Fmm::new(&kernel, moved, 6, 32)
        .expect("admitted shifted fixture")
        .potentials(&q)
        .expect("finite fixture");
    let mut worst = 0.0f64;
    for (a, b) in base.iter().zip(&shifted) {
        worst = worst.max((a - b).abs() / a.abs().max(1e-12));
    }
    verdict(
        "fmm-002",
        worst < 1e-9,
        &format!(
            "\"detail\":\"G3: rigidly shifted cloud, same potentials\",\
             \"worst_rel\":{worst:.3e}"
        ),
    );
}

// ------------------------------------------------------------------ fmm-003

#[test]
fn fmm_003_scaling_trend() {
    let kernel = Laplace3d;
    let sizes = [4096usize, 8192, 16384, 32768];
    let mut times = Vec::new();
    let mut rows = String::new();
    for &n in &sizes {
        let (pts, q) = cloud(n, 0x1001_2026_0708_0003);
        let fmm = Fmm::new(&kernel, pts, 4, 48).expect("admitted scaling fixture");
        let t0 = Instant::now();
        let out = fmm.potentials(&q).expect("finite fixture");
        let dt = t0.elapsed().as_secs_f64();
        assert!(out.iter().all(|v| v.is_finite()), "finite potentials");
        let _ = write!(rows, "{{\"n\":{n},\"seconds\":{dt:.3}}},");
        times.push(dt);
    }
    // Fitted exponent over the doubling ladder.
    let mut exps = Vec::new();
    for w in times.windows(2) {
        exps.push((w[1] / w[0]).log2());
    }
    #[allow(clippy::cast_precision_loss)]
    let mean_exp = exps.iter().sum::<f64>() / exps.len() as f64;
    // O(N log N)-class: comfortably below the direct method's 2.
    let pass = mean_exp < 1.6;
    verdict(
        "fmm-003",
        pass,
        &format!(
            "\"detail\":\"time-vs-N trend (order 4); 1e7-point wall-clock is perf-lane scope\",\
             \"rows\":[{}],\"fitted_exponent\":{mean_exp:.2}",
            rows.trim_end_matches(',')
        ),
    );
}

// ------------------------------------------------------------------ totality

#[test]
fn fmm_rejects_invalid_inputs_and_work_before_evaluation() {
    let kernel = Laplace3d;
    assert!(Fmm::new(&kernel, Vec::new(), 4, 32).is_err());
    assert!(Fmm::new(&kernel, vec![[0.0; 3]], usize::MAX, 32).is_err());
    assert!(Fmm::new(&kernel, vec![[f64::NAN, 0.0, 0.0]], 4, 32).is_err());
    assert!(Fmm::new(&kernel, vec![[0.0; 3]], 4, 0).is_err());
    assert!(Fmm::new(&kernel, vec![[0.0; 3]; 1_500], 12, 32).is_err());
    assert!(Fmm::new(&kernel, vec![[0.0; 3]; 20_000], 2, 32).is_err());

    let fmm =
        Fmm::new(&kernel, vec![[0.0; 3], [1.0, 0.0, 0.0]], 4, 2).expect("small valid fixture");
    assert!(fmm.potentials(&[1.0]).is_err());
    assert!(fmm.potentials(&[1.0, f64::INFINITY]).is_err());
}
