//! Sinkhorn OT battery (vcia slice e): marginal feasibility; the 1D
//! quadratic-cost KNOWN ANSWER (monotone coupling) approached as the
//! regularization shrinks (measured ladder); symmetry; translation
//! covariance (W₂² of equal translates = |t|² exactly in the limit);
//! determinism; golden.

use fs_dfo::{cost_sq_1d, monotone_cost_1d, sinkhorn};
use fs_rand::StreamKey;

const SUITE: &str = "fs-dfo-ot";
const INPUT_SEED: u64 = 111;
const STREAM_KERNEL: u32 = 0x0007;
const MARGINAL_X_TILE: u32 = 1;
const MARGINAL_Y_TILE: u32 = 2;
const LADDER_X_TILE: u32 = 3;
const LADDER_Y_TILE: u32 = 4;
const TRANSLATION_TILE: u32 = 5;
const GOLDEN_X_TILE: u32 = 6;
const GOLDEN_Y_TILE: u32 = 7;

fn verdict(case: &str, pass: bool, detail: &str, seed: u64) {
    let mut emitter = fs_obs::Emitter::new(SUITE, case);
    let event = emitter.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: case.to_string(),
            pass,
            detail: detail.to_string(),
            seed,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("OT verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("OT verdict must use the fs-obs wire schema");
    println!("{line}");
    assert!(pass, "case {case}: {detail}");
}

fn measurement(case: &str, json: String) {
    let mut emitter = fs_obs::Emitter::new(SUITE, format!("{case}/measurement"));
    let event = emitter.emit(
        fs_obs::Severity::Info,
        fs_obs::EventKind::Custom {
            name: case.to_string(),
            json,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("OT measurement must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("OT measurement must use the fs-obs wire schema");
    println!("{line}");
}

fn rand_pts(n: usize, tile: u32) -> Vec<f64> {
    let mut s = StreamKey {
        seed: INPUT_SEED,
        kernel: STREAM_KERNEL,
        tile,
    }
    .stream();
    (0..n).map(|_| s.next_f64()).collect()
}

#[test]
fn marginals_and_symmetry() {
    let x = rand_pts(12, MARGINAL_X_TILE);
    let y = rand_pts(15, MARGINAL_Y_TILE);
    let a = vec![1.0 / 12.0; 12];
    let b = vec![1.0 / 15.0; 15];
    let c = cost_sq_1d(&x, &y);
    let rep = sinkhorn(&a, &b, &c, 0.01, 2000);
    assert!(
        rep.marginal_residual < 1e-8,
        "marginals violated: {:.3e}",
        rep.marginal_residual
    );
    // Symmetry: transpose the problem.
    let mut ct = vec![0.0f64; c.len()];
    for i in 0..12 {
        for j in 0..15 {
            ct[j * 12 + i] = c[i * 15 + j];
        }
    }
    let rep_t = sinkhorn(&b, &a, &ct, 0.01, 2000);
    let rel = (rep.cost - rep_t.cost).abs() / rep.cost.max(1e-30);
    assert!(
        rel < 1e-8,
        "W(a,b) != W(b,a): {} vs {}",
        rep.cost,
        rep_t.cost
    );
    verdict(
        "marginals-symmetry",
        true,
        &format!(
            "residual {:.1e}, sym rel {rel:.1e}, iters {}; input seed {INPUT_SEED}, \
             stream kernel {STREAM_KERNEL:#x}, tiles {MARGINAL_X_TILE}/{MARGINAL_Y_TILE}",
            rep.marginal_residual, rep.iters
        ),
        INPUT_SEED,
    );
}

#[test]
fn epsilon_ladder_approaches_monotone_coupling() {
    let n = 16usize;
    let x = rand_pts(n, LADDER_X_TILE);
    let y: Vec<f64> = rand_pts(n, LADDER_Y_TILE).iter().map(|v| v + 0.4).collect();
    let a = vec![1.0 / n as f64; n];
    let c = cost_sq_1d(&x, &y);
    let truth = monotone_cost_1d(&x, &y);
    let mut prev_gap = f64::INFINITY;
    let mut gaps = Vec::new();
    for &eps in &[0.05f64, 0.01, 0.002] {
        let rep = sinkhorn(&a, &a, &c, eps, 20_000);
        let gap = (rep.cost - truth).abs() / truth;
        gaps.push(format!("eps={eps}: {gap:.4}"));
        assert!(
            gap < prev_gap + 1e-9,
            "entropic cost must approach the monotone optimum as eps drops: {gaps:?}"
        );
        prev_gap = gap;
    }
    assert!(
        prev_gap < 0.02,
        "smallest-eps cost still far from the closed form: {gaps:?}"
    );
    verdict(
        "eps-ladder",
        true,
        &format!(
            "{}; input seed {INPUT_SEED}, stream kernel {STREAM_KERNEL:#x}, tiles \
             {LADDER_X_TILE}/{LADDER_Y_TILE}",
            gaps.join(", ")
        ),
        INPUT_SEED,
    );
}

#[test]
fn translation_covariance() {
    // W₂²(μ, μ + t) = t² for equal translates (every coupling moves
    // mass exactly t in the monotone limit).
    let n = 20usize;
    let x = rand_pts(n, TRANSLATION_TILE);
    let t = 0.7f64;
    let y: Vec<f64> = x.iter().map(|v| v + t).collect();
    let a = vec![1.0 / n as f64; n];
    let c = cost_sq_1d(&x, &y);
    let rep = sinkhorn(&a, &a, &c, 0.002, 20_000);
    let rel = (rep.cost - t * t).abs() / (t * t);
    assert!(rel < 0.02, "translate cost {} vs t^2 {}", rep.cost, t * t);
    verdict(
        "translation",
        true,
        &format!(
            "cost {:.5} vs {:.5}, rel {rel:.4}; input seed {INPUT_SEED}, stream kernel \
             {STREAM_KERNEL:#x}, tile {TRANSLATION_TILE}",
            rep.cost,
            t * t
        ),
        INPUT_SEED,
    );
}

const GOLDEN_HASH: u64 = 0x58eb_8443_224c_a689; // recorded at vcia slice e, frozen

#[test]
fn ot_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    let x = rand_pts(10, GOLDEN_X_TILE);
    let y = rand_pts(8, GOLDEN_Y_TILE);
    let a = vec![0.1f64; 10];
    let b = vec![0.125f64; 8];
    let c = cost_sq_1d(&x, &y);
    let rep = sinkhorn(&a, &b, &c, 0.02, 3000);
    feed(rep.cost);
    feed(rep.marginal_residual);
    for v in rep.plan.iter().step_by(7) {
        feed(*v);
    }
    measurement(
        "ot-golden",
        format!(
            "{{\"actual\":\"{acc:#018x}\",\"expected\":\"{GOLDEN_HASH:#018x}\",\
             \"input_seed\":{INPUT_SEED},\"stream_kernel\":{STREAM_KERNEL},\
             \"x_stream_tile\":{GOLDEN_X_TILE},\"y_stream_tile\":{GOLDEN_Y_TILE}}}"
        ),
    );
    assert_eq!(
        acc, GOLDEN_HASH,
        "ot bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}

#[test]
#[should_panic(expected = "epsilon must be positive")]
fn sinkhorn_rejects_a_nonpositive_epsilon() {
    // Regression: the Gibbs kernel exp(-c/epsilon) is 0/inf -> NaN potentials
    // for epsilon <= 0. Fail closed instead of returning a NaN plan.
    let a = [0.5, 0.5];
    let b = [0.5, 0.5];
    let c = [0.0, 1.0, 1.0, 0.0]; // 2x2 cost matrix
    let _ = sinkhorn(&a, &b, &c, 0.0, 100);
}
