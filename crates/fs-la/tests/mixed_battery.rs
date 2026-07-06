//! Mixed-precision battery: condition-controlled matrices (A = Q₁·D·Q₂ᵀ
//! with chosen spectrum), verifying the ladder policy, escalation
//! triggers, dd-residual forward-error gains, honest run-dry reporting,
//! and cross-ISA determinism of the whole decision process.

use fs_la::factor::qr;
use fs_la::mixed::{Ladder, ResidualTarget, solve_adaptive};

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
}

/// Random orthogonal n×n (columns of Q from a QR of a random matrix).
fn orthogonal(n: usize, seed: u64) -> Vec<f64> {
    let mut s = seed;
    let a: Vec<f64> = (0..n * n).map(|_| lcg(&mut s)).collect();
    let f = qr(&a, n, n);
    let mut q = vec![0.0; n * n];
    for j in 0..n {
        let mut e = vec![0.0; n];
        e[j] = 1.0;
        f.apply_q(&mut e);
        for i in 0..n {
            q[i * n + j] = e[i];
        }
    }
    q
}

/// A = Q₁·diag(σ)·Q₂ᵀ with σ log-spaced from 1 down to 1/κ — spectral
/// condition exactly κ.
fn conditioned(n: usize, kappa: f64, seed: u64) -> Vec<f64> {
    let q1 = orthogonal(n, seed);
    let q2 = orthogonal(n, seed ^ 0xFFFF);
    let mut a = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            let mut acc = 0.0f64;
            for k in 0..n {
                let sigma = kappa.powf(-(k as f64) / ((n - 1) as f64));
                acc = (q1[i * n + k] * sigma).mul_add(q2[j * n + k], acc);
            }
            a[i * n + j] = acc;
        }
    }
    a
}

fn rhs_for(a: &[f64], n: usize, x_true: &[f64]) -> Vec<f64> {
    let mut b = vec![0.0; n];
    for i in 0..n {
        let mut acc = 0.0f64;
        for (j, &xj) in x_true.iter().enumerate() {
            acc = a[i * n + j].mul_add(xj, acc);
        }
        b[i] = acc;
    }
    b
}

fn forward_error(x: &[f64], x_true: &[f64]) -> f64 {
    let num = x
        .iter()
        .zip(x_true)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, f64::max);
    let den = x_true.iter().fold(0.0f64, |m, &v| m.max(v.abs()));
    num / den
}

const N: usize = 24;

#[test]
fn f32_rung_converges_on_well_conditioned() {
    for kappa in [1e2, 1e4] {
        let a = conditioned(N, kappa, 0xAB);
        let x_true: Vec<f64> = (0..N).map(|i| (i as f64).sin() + 1.5).collect();
        let b = rhs_for(&a, N, &x_true);
        let (x, rep) = solve_adaptive(
            &a,
            N,
            &b,
            ResidualTarget {
                backward: 1e-14,
                forward: None,
            },
        )
        .expect("nonsingular");
        assert_eq!(rep.ladder, Ladder::F32Refine, "kappa={kappa}: {rep:?}");
        assert!(rep.converged && !rep.escalated, "kappa={kappa}: {rep:?}");
        assert!(rep.achieved <= 1e-14);
        let fe = forward_error(&x, &x_true);
        assert!(
            fe <= kappa * 1e-12,
            "forward error {fe:.2e} at kappa={kappa}"
        );
        // Trajectory is decreasing (refinement contracts).
        for w in rep.trajectory.windows(2) {
            assert!(
                w[1] <= w[0] * 1.5,
                "trajectory not contracting: {:?}",
                rep.trajectory
            );
        }
        println!(
            "{{\"suite\":\"fs-la\",\"case\":\"mixed-f32\",\"verdict\":\"pass\",\"detail\":\"kappa={kappa:.0e} steps={} achieved={:.2e}\"}}",
            rep.steps, rep.achieved
        );
    }
}

#[test]
fn policy_vetoes_f32_on_ill_conditioned() {
    let a = conditioned(N, 1e10, 0xCD);
    let x_true: Vec<f64> = (0..N).map(|i| 1.0 + (i as f64) * 0.1).collect();
    let b = rhs_for(&a, N, &x_true);
    let (_, rep) = solve_adaptive(
        &a,
        N,
        &b,
        ResidualTarget {
            backward: 1e-13,
            forward: None,
        },
    )
    .expect("nonsingular");
    assert!(
        rep.escalated,
        "condition evidence must veto the f32 rung: {rep:?}"
    );
    assert_ne!(rep.ladder, Ladder::F32Refine);
    assert!(
        rep.converged,
        "f64 direct meets a 1e-13 backward target: {rep:?}"
    );
    assert!(
        rep.condition_estimate > 1e8,
        "estimate {:.2e} must expose the ill-conditioning",
        rep.condition_estimate
    );
}

#[test]
fn dd_rung_improves_forward_error() {
    let a = conditioned(N, 1e10, 0xEF);
    let x_seed: Vec<f64> = (0..N).map(|i| ((i as f64) * 0.7).cos() + 2.0).collect();
    let b = rhs_for(&a, N, &x_seed);
    // GROUND TRUTH is the exact solution of the STORED (A, b): rounding b
    // to f64 already moved the true solution by ~kappa*eps away from
    // x_seed, and no refinement can undo that. Compute the reference by
    // running dd refinement far past convergence.
    let f = fs_la::factor::lu(&a, N).unwrap();
    let mut x_star = b.clone();
    f.solve(&mut x_star);
    for _ in 0..12 {
        let mut r = vec![0.0f64; N];
        for i in 0..N {
            let mut acc = fs_math::dd::Dd::from_f64(b[i]);
            for (j, &xj) in x_star.iter().enumerate() {
                acc = acc - fs_math::dd::Dd::from_f64(a[i * N + j]) * fs_math::dd::Dd::from_f64(xj);
            }
            r[i] = acc.to_f64();
        }
        f.solve(&mut r);
        for (xi, di) in x_star.iter_mut().zip(&r) {
            *xi += di;
        }
    }
    let x_true = x_star;
    let mut x_direct = b.clone();
    f.solve(&mut x_direct);
    let fe_direct = forward_error(&x_direct, &x_true);
    assert!(
        fe_direct > 1e-12,
        "direct must visibly miss at kappa=1e10: {fe_direct:.2e}"
    );
    // A forward target below κ·eps engages the dd rung (backward alone is
    // already met by any stable solve — that's the point of the test).
    let target = ResidualTarget {
        backward: 1e-14,
        forward: Some(1e-12),
    };
    let (x_dd, rep) = solve_adaptive(&a, N, &b, target).expect("nonsingular");
    assert_eq!(rep.ladder, Ladder::F64DdRefine, "{rep:?}");
    assert!(
        rep.converged,
        "dd refinement reaches the forward target: {rep:?}"
    );
    assert!(
        rep.forward_estimate.unwrap() <= 1e-12,
        "forward estimate must certify the target: {rep:?}"
    );
    let fe_dd = forward_error(&x_dd, &x_true);
    assert!(
        fe_dd < fe_direct * 1e-2,
        "dd refinement must beat direct forward error by >=100x: {fe_dd:.2e} vs {fe_direct:.2e}"
    );
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"mixed-dd\",\"verdict\":\"pass\",\"detail\":\"kappa=1e10 forward {fe_direct:.2e} -> {fe_dd:.2e} in {} steps\"}}",
        rep.steps
    );
}

#[test]
fn singular_in_f32_escalates_to_f64() {
    // 1 + 1e-12 collapses to 1.0 in f32 → exactly singular there; f64 is
    // fine (kappa ~ 4e12, and a 1e-10 backward target is easy for direct).
    let a = vec![1.0, 1.0, 1.0, 1.0 + 1e-12];
    let b = vec![2.0, 2.0 + 2e-12];
    let (x, rep) = solve_adaptive(
        &a,
        2,
        &b,
        ResidualTarget {
            backward: 1e-10,
            forward: None,
        },
    )
    .expect("f64 nonsingular");
    assert!(rep.escalated, "f32 singularity must escalate: {rep:?}");
    assert!(rep.converged);
    // True solution of this system is x = [0, 2] (the 1e-12 perturbation
    // sits in row 2 only).
    assert!(
        (x[0] - 0.0).abs() < 1e-3 && (x[1] - 2.0).abs() < 1e-3,
        "x = {x:?}"
    );
}

#[test]
fn impossible_target_reports_honestly() {
    let a = conditioned(N, 1e6, 0x11);
    let x_true = vec![1.0; N];
    let b = rhs_for(&a, N, &x_true);
    let (_, rep) = solve_adaptive(
        &a,
        N,
        &b,
        ResidualTarget {
            backward: 1e-40,
            forward: None,
        },
    )
    .expect("nonsingular");
    assert!(!rep.converged, "1e-40 backward is unreachable: {rep:?}");
    assert_eq!(
        rep.ladder,
        Ladder::F64DdRefine,
        "must have climbed the whole ladder"
    );
    assert!(rep.escalated);
    assert!(
        rep.achieved > 0.0 && rep.achieved < 1e-13,
        "best-achieved recorded"
    );
    assert!(!rep.trajectory.is_empty());
}

#[test]
fn policy_and_solution_are_deterministic() {
    let a = conditioned(N, 1e5, 0x22);
    let x_true: Vec<f64> = (0..N).map(|i| 1.0 + (i as f64) * 0.01).collect();
    let b = rhs_for(&a, N, &x_true);
    let (x1, r1) = solve_adaptive(
        &a,
        N,
        &b,
        ResidualTarget {
            backward: 1e-14,
            forward: None,
        },
    )
    .unwrap();
    let (x2, r2) = solve_adaptive(
        &a,
        N,
        &b,
        ResidualTarget {
            backward: 1e-14,
            forward: None,
        },
    )
    .unwrap();
    assert_eq!(r1.ladder, r2.ladder);
    assert_eq!(r1.steps, r2.steps);
    assert!(x1.iter().zip(&x2).all(|(p, q)| p.to_bits() == q.to_bits()));
    assert!(
        r1.trajectory
            .iter()
            .zip(&r2.trajectory)
            .all(|(p, q)| p.to_bits() == q.to_bits()),
        "trajectory must be bit-stable"
    );
}

/// Recorded on aarch64-apple (M4 Pro); must match on x86-64 (trj).
const GOLDEN_HASH: u64 = 0xfeaa_02d5_a8b3_5aa9;

#[test]
fn mixed_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    for (kappa, seed) in [(1e3, 0x31u64), (1e7, 0x32), (1e11, 0x33)] {
        let a = conditioned(N, kappa, seed);
        let x_true: Vec<f64> = (0..N).map(|i| 1.0 + (i as f64) * 0.05).collect();
        let b = rhs_for(&a, N, &x_true);
        let (x, rep) = solve_adaptive(
            &a,
            N,
            &b,
            ResidualTarget {
                backward: 1e-15,
                forward: Some(1e-11),
            },
        )
        .unwrap();
        for &v in &x {
            feed(v);
        }
        feed(rep.achieved);
        feed(rep.steps as f64);
        feed(match rep.ladder {
            Ladder::F32Refine => 1.0,
            Ladder::F64Direct => 2.0,
            Ladder::F64DdRefine => 3.0,
        });
    }
    println!(
        "{{\"suite\":\"fs-la\",\"case\":\"mixed-golden\",\"verdict\":\"info\",\"detail\":\"{acc:#018x}\"}}"
    );
    assert_eq!(
        acc, GOLDEN_HASH,
        "mixed-precision decision bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only \
         with semantic justification (golden-evidence policy)"
    );
}
