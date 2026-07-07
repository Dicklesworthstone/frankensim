//! fs-time battery (tfz.12): symplectic energy boundedness vs RK4's
//! secular drift, the variational-integrator equivalence, Lie-group
//! norm preservation + gyroscope physics, generalized-α spectral
//! behavior against Chung–Hulbert theory, IMEX stiff stability +
//! second-order convergence, exponential-Euler exactness on linears,
//! embedded-pair adaptivity with bitwise resumability, the discrete
//! adjoint gradcheck, and the cross-ISA golden hash.
//!
//! Test-side ORACLES may use std/platform libm (disjoint-path rule);
//! everything feeding the golden hash flows through solver code only.

use fs_time::{
    AdaptiveState, ExpEuler, GeneralizedAlpha, Imex2, PiController, galpha_step, imex2_step,
    quat_exp_step, quat_rotate, rigid_body_step, rk45_adaptive, verlet_adjoint, verlet_step,
};

fn log(case: &str, verdict: &str, detail: &str) {
    println!("{{\"suite\":\"fs-time\",\"case\":\"{case}\",\"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}");
}

// ---------------------------------------------------------------- symplectic

/// Harmonic-oscillator energy for unit mass/stiffness.
fn ho_energy(q: &[f64], p: &[f64]) -> f64 {
    0.5 * (q[0] * q[0] + p[0] * p[0])
}

/// Classic RK4 on the first-order form of q̈ = −q (the drift comparator).
fn rk4_ho_step(q: &mut f64, p: &mut f64, h: f64) {
    let f = |q: f64, p: f64| (p, -q);
    let (k1q, k1p) = f(*q, *p);
    let (k2q, k2p) = f(*q + 0.5 * h * k1q, *p + 0.5 * h * k1p);
    let (k3q, k3p) = f(*q + 0.5 * h * k2q, *p + 0.5 * h * k2p);
    let (k4q, k4p) = f(*q + h * k3q, *p + h * k3p);
    *q += h / 6.0 * (k1q + 2.0 * k2q + 2.0 * k3q + k4q);
    *p += h / 6.0 * (k1p + 2.0 * k2p + 2.0 * k3p + k4p);
}

#[test]
fn verlet_energy_bounded_1e6_steps_vs_rk4_drift() {
    let h = 0.1;
    let steps = 1_000_000usize;
    let force = |q: &[f64], out: &mut [f64]| out[0] = -q[0];
    let (mut q, mut p) = (vec![1.0f64], vec![0.0f64]);
    let mut scratch = vec![0.0f64; 1];
    let e0 = ho_energy(&q, &p);
    let (mut max_dev_first, mut max_dev_second) = (0.0f64, 0.0f64);
    for k in 0..steps {
        verlet_step(&mut q, &mut p, h, &force, &mut scratch);
        let dev = (ho_energy(&q, &p) - e0).abs();
        if k < steps / 2 {
            max_dev_first = max_dev_first.max(dev);
        } else {
            max_dev_second = max_dev_second.max(dev);
        }
    }
    // Symplectic ⇒ BOUNDED oscillating energy error: the second half of
    // the run must not exceed the first half beyond roundoff slack.
    assert!(
        max_dev_second <= 1.05 * max_dev_first + 1e-12,
        "secular drift in Verlet: first-half max {max_dev_first:.3e}, second {max_dev_second:.3e}"
    );
    assert!(max_dev_first < 2e-3, "Verlet energy deviation too large: {max_dev_first:.3e}");
    // RK4 comparator: same h, same span — visible SECULAR energy decay.
    let (mut qr, mut pr) = (1.0f64, 0.0f64);
    for _ in 0..steps {
        rk4_ho_step(&mut qr, &mut pr, h);
    }
    let rk4_drift = (0.5 * (qr * qr + pr * pr) - e0).abs();
    assert!(
        rk4_drift > 5.0 * max_dev_second,
        "RK4 should drift visibly: {rk4_drift:.3e} vs Verlet bound {max_dev_second:.3e}"
    );
    log(
        "verlet-vs-rk4",
        "pass",
        &format!("verlet_bound={max_dev_second:.3e} rk4_drift={rk4_drift:.3e}"),
    );
}

#[test]
fn verlet_energy_bounded_kepler_eccentric() {
    // e = 0.6 Kepler orbit, ~16 revolutions: energy bounded, no drift.
    let h = 1e-3;
    let steps = 100_000usize;
    let force = |q: &[f64], out: &mut [f64]| {
        let r2 = q[0] * q[0] + q[1] * q[1];
        let r3 = r2 * r2.sqrt();
        out[0] = -q[0] / r3;
        out[1] = -q[1] / r3;
    };
    let energy = |q: &[f64], p: &[f64]| {
        0.5 * (p[0] * p[0] + p[1] * p[1]) - 1.0 / (q[0] * q[0] + q[1] * q[1]).sqrt()
    };
    let (mut q, mut p) = (vec![0.4f64, 0.0], vec![0.0f64, 2.0]);
    let mut scratch = vec![0.0f64; 2];
    let e0 = energy(&q, &p);
    let mut max_dev = 0.0f64;
    for _ in 0..steps {
        verlet_step(&mut q, &mut p, h, &force, &mut scratch);
        max_dev = max_dev.max((energy(&q, &p) - e0).abs());
    }
    assert!(max_dev < 5e-4, "Kepler energy deviation {max_dev:.3e}");
    log("verlet-kepler", "pass", &format!("e=0.6 max_energy_dev={max_dev:.3e}"));
}

#[test]
fn verlet_is_the_variational_integrator() {
    // Marsden–West: extremizing Σ h·[½‖(q_{k+1}−q_k)/h‖² − ½(V_k+V_{k+1})]
    // gives the position two-step q_{k+1} = 2q_k − q_{k−1} + h²F(q_k).
    // The kick–drift–kick positions must satisfy it (same method, two
    // derivations) to accumulated-roundoff accuracy on a NONLINEAR V.
    let h = 0.05;
    let force = |q: &[f64], out: &mut [f64]| out[0] = -q[0].sin(); // pendulum
    let (mut q, mut p) = (vec![1.2f64], vec![0.3f64]);
    let mut scratch = vec![0.0f64; 1];
    let mut positions = vec![q[0]];
    for _ in 0..200 {
        verlet_step(&mut q, &mut p, h, &force, &mut scratch);
        positions.push(q[0]);
    }
    let mut worst = 0.0f64;
    for k in 1..200 {
        let predicted = 2.0 * positions[k] - positions[k - 1] - h * h * positions[k].sin();
        worst = worst.max((positions[k + 1] - predicted).abs());
    }
    assert!(worst < 1e-12, "discrete Euler–Lagrange residual {worst:.3e}");
    log("verlet-variational", "pass", &format!("max_dEL_residual={worst:.3e}"));
}

// ----------------------------------------------------------------- Lie/SO(3)

#[test]
fn quat_norm_preserved_1e5_steps() {
    // Time-varying ω from a deterministic recurrence; exp-map updates
    // must keep ‖q‖ = 1 to a roundoff random walk (~√N·ε ≈ 4e-14).
    let mut q = [1.0f64, 0.0, 0.0, 0.0];
    let mut w = [0.3f64, -0.2, 0.9];
    for k in 0..100_000usize {
        q = quat_exp_step(q, w, 0.01);
        let s = if k % 2 == 0 { 1.0 } else { -1.0 };
        w = [w[1], w[2], w[0] + s * 1e-4];
    }
    let norm = (q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]).sqrt();
    assert!((norm - 1.0).abs() < 1e-12, "quaternion norm drift: {:.3e}", (norm - 1.0).abs());
    log("quat-norm", "pass", &format!("drift={:.3e}", (norm - 1.0).abs()));
}

#[test]
fn rigid_body_gyroscope_physics() {
    // Symmetric top I = (2, 2, 1), ω = (0.1, 0, 5): Euler's equations
    // give ω₃ EXACTLY constant and (ω₁, ω₂) precessing in the body
    // frame at Ω = (I₃ − I₁)/I₁·ω₃ = −2.5 rad/s. Second-order method ⇒
    // both hold to O(h²) over T = 10.
    let inertia = [2.0f64, 2.0, 1.0];
    let h = 1e-3;
    let steps = 10_000usize;
    let mut q = [1.0f64, 0.0, 0.0, 0.0];
    let mut w = [0.1f64, 0.0, 5.0];
    let l0 = quat_rotate(q, [inertia[0] * w[0], inertia[1] * w[1], inertia[2] * w[2]]);
    let e0 = 0.5 * (inertia[0] * w[0] * w[0] + inertia[1] * w[1] * w[1] + inertia[2] * w[2] * w[2]);
    for _ in 0..steps {
        let (qn, wn) = rigid_body_step(q, w, inertia, h);
        q = qn;
        w = wn;
    }
    // (a) ω₃ constant.
    assert!((w[2] - 5.0).abs() < 1e-9, "omega3 drift {:.3e}", (w[2] - 5.0).abs());
    // (b) body-frame precession phase after T = 10: Ω·T = −25 rad.
    let phase = w[1].atan2(w[0]);
    let expected = -25.0f64.rem_euclid(2.0 * std::f64::consts::PI);
    let mut diff = (phase - (2.0 * std::f64::consts::PI - expected)).abs();
    // Fold to [0, π].
    while diff > std::f64::consts::PI {
        diff = (diff - 2.0 * std::f64::consts::PI).abs();
    }
    assert!(diff < 1e-3, "precession phase error {diff:.3e}");
    // (c) energy and SPATIAL angular momentum conserved to O(h²).
    let e1 = 0.5 * (inertia[0] * w[0] * w[0] + inertia[1] * w[1] * w[1] + inertia[2] * w[2] * w[2]);
    assert!((e1 - e0).abs() / e0 < 1e-6, "energy drift {:.3e}", (e1 - e0).abs() / e0);
    let l1 = quat_rotate(q, [inertia[0] * w[0], inertia[1] * w[1], inertia[2] * w[2]]);
    let ldev = (0..3).map(|i| (l1[i] - l0[i]).abs()).fold(0.0f64, f64::max);
    assert!(ldev < 1e-3, "spatial L drift {ldev:.3e}");
    log(
        "gyroscope",
        "pass",
        &format!("phase_err={diff:.3e} L_dev={ldev:.3e} E_rel={:.3e}", (e1 - e0).abs() / e0),
    );
}

// ------------------------------------------------------------- generalized-α

#[test]
fn galpha_high_frequency_dissipation_matches_rho_inf() {
    // 1-DOF, ωh = 10³ (deep in the high-frequency limit): the numerical
    // amplification's asymptotic per-step contraction must approach ρ∞.
    for &rho in &[0.0f64, 0.5, 0.9] {
        let ga = GeneralizedAlpha::new(&[1.0], &[0.0], &[1.0e6], 1, 1.0, rho);
        let (mut q, mut v, mut a) = (vec![1.0f64], vec![0.0f64], vec![-1.0e6f64]);
        let f = vec![0.0f64];
        let mut prev_norm = 0.0f64;
        let mut ratio = 0.0f64;
        for k in 0..60 {
            galpha_step(&ga, &mut q, &mut v, &mut a, &f);
            let norm = (q[0] * q[0] + (v[0] / 1.0e3) * (v[0] / 1.0e3)).sqrt();
            if k >= 40 && prev_norm > 0.0 {
                ratio = norm / prev_norm;
            }
            prev_norm = norm;
        }
        if rho == 0.0 {
            // Asymptotic annihilation: state collapses towards zero.
            assert!(prev_norm < 1e-30, "rho=0 should annihilate: {prev_norm:.3e}");
        } else {
            assert!(
                (ratio - rho).abs() < 0.02,
                "spectral radius {ratio:.4} vs rho_inf {rho} at omega*h=1e3"
            );
        }
        log("galpha-spectral", "pass", &format!("rho_inf={rho} measured={ratio:.4}"));
    }
}

#[test]
fn galpha_no_dissipation_at_rho_one_and_order_two() {
    // ρ∞ = 1 (trapezoidal limit): free vibration amplitude preserved.
    let ga = GeneralizedAlpha::new(&[1.0], &[0.0], &[1.0], 1, 0.1, 1.0);
    let (mut q, mut v, mut a) = (vec![1.0f64], vec![0.0f64], vec![-1.0f64]);
    let f = vec![0.0f64];
    for _ in 0..10_000 {
        galpha_step(&ga, &mut q, &mut v, &mut a, &f);
    }
    let energy = 0.5 * (q[0] * q[0] + v[0] * v[0]);
    assert!((energy - 0.5).abs() < 1e-3, "rho=1 energy drift {:.3e}", (energy - 0.5).abs());
    // O(h²): error in q(2π) vs cos, halving h → error / 4 (±30%).
    let err_at = |h: f64| -> f64 {
        let steps = (2.0 * std::f64::consts::PI / h).round() as usize;
        let hh = 2.0 * std::f64::consts::PI / steps as f64;
        let ga = GeneralizedAlpha::new(&[1.0], &[0.0], &[1.0], 1, hh, 0.8);
        let (mut q, mut v, mut a) = (vec![1.0f64], vec![0.0f64], vec![-1.0f64]);
        let f = vec![0.0f64];
        for _ in 0..steps {
            galpha_step(&ga, &mut q, &mut v, &mut a, &f);
        }
        (q[0] - 1.0).abs()
    };
    let (e1, e2) = (err_at(0.02), err_at(0.01));
    let order = (e1 / e2).log2();
    assert!((order - 2.0).abs() < 0.4, "generalized-alpha order {order:.2} (errors {e1:.3e}/{e2:.3e})");
    log("galpha-order", "pass", &format!("order={order:.2}"));
}

// --------------------------------------------------------------------- stiff

#[test]
fn imex_stiff_stability() {
    // hλ = −100 on the implicit part: explicit methods explode, the
    // ARS(2,2,2) step contracts monotonically. Mild monotone nonlinearity.
    let im = Imex2::new(&[-1.0e4], 1, 0.01);
    let nonlin = |u: &[f64], out: &mut [f64]| out[0] = -u[0] * u[0] * u[0];
    let mut u = vec![1.0f64];
    let mut prev = u[0].abs();
    for _ in 0..200 {
        imex2_step(&im, &mut u, &nonlin);
        assert!(u[0].abs() <= prev + 1e-15, "IMEX not contracting: {} -> {}", prev, u[0]);
        prev = u[0].abs();
    }
    assert!(u[0].abs() < 1e-8, "stiff mode not damped: {:.3e}", u[0]);
    log("imex-stability", "pass", &format!("u_final={:.3e}", u[0]));
}

#[test]
fn imex_second_order_on_logistic() {
    // u′ = −u + u², u₀ = ½: exact u(t) = u₀ / (u₀ + (1−u₀)eᵗ).
    // L = −1 implicit, N(u) = u² explicit — order 2 requires the ARS
    // (δ, 1−δ) weights (trapezoidal ½,½ degrades to order 1 HERE).
    let exact = |t: f64| 0.5 / (0.5 + 0.5 * t.exp());
    let err_at = |h: f64| -> f64 {
        let steps = (1.0 / h).round() as usize;
        let im = Imex2::new(&[-1.0], 1, h);
        let nonlin = |u: &[f64], out: &mut [f64]| out[0] = u[0] * u[0];
        let mut u = vec![0.5f64];
        for _ in 0..steps {
            imex2_step(&im, &mut u, &nonlin);
        }
        (u[0] - exact(1.0)).abs()
    };
    let (e1, e2, e3) = (err_at(0.1), err_at(0.05), err_at(0.025));
    let (o1, o2) = ((e1 / e2).log2(), (e2 / e3).log2());
    assert!(
        (o1 - 2.0).abs() < 0.35 && (o2 - 2.0).abs() < 0.35,
        "IMEX order ratios {o1:.2}, {o2:.2} (errors {e1:.3e}/{e2:.3e}/{e3:.3e})"
    );
    log("imex-order", "pass", &format!("orders={o1:.2},{o2:.2}"));
}

#[test]
fn exp_euler_exact_on_linear() {
    // A = QᵀDQ (rotation θ = 0.3, D = diag(−2, −3)); with N ≡ 0 the
    // exponential integrator must reproduce u(t) = Qᵀe^{tD}Qu₀ to
    // roundoff over 100 steps — the oracle path (std exp + hand
    // rotation) is DISJOINT from the solver path (Jacobi eigenbasis).
    let (c, s) = (0.3f64.cos(), 0.3f64.sin());
    // A = Qᵀ D Q with Q = [[c, s], [−s, c]].
    let d = [-2.0f64, -3.0];
    let a = [
        c * d[0] * c + s * d[1] * s,
        c * d[0] * s - s * d[1] * c,
        s * d[0] * c - c * d[1] * s,
        s * d[0] * s + c * d[1] * c,
    ];
    let h = 0.05;
    let ee = ExpEuler::new(&a, 2, h);
    let zero = |_: &[f64], out: &mut [f64]| out.fill(0.0);
    let mut u = vec![0.7f64, -0.4];
    for _ in 0..100 {
        ee.step(&mut u, &zero);
    }
    let t = h * 100.0;
    let u0 = [0.7f64, -0.4];
    let w = [c * u0[0] + s * u0[1], -s * u0[0] + c * u0[1]];
    let w = [w[0] * (d[0] * t).exp(), w[1] * (d[1] * t).exp()];
    let expect = [c * w[0] - s * w[1], s * w[0] + c * w[1]];
    let err = (u[0] - expect[0]).abs().max((u[1] - expect[1]).abs());
    assert!(err < 1e-13, "exp-Euler not exact on linear: {err:.3e}");
    log("expeuler-exact", "pass", &format!("err={err:.3e}"));
}

#[test]
fn exp_euler_first_order_with_nonlinearity() {
    // ETD1 is order 1 in N; verify the convergence rate on the logistic
    // problem (same exact solution as the IMEX test).
    let exact = |t: f64| 0.5 / (0.5 + 0.5 * t.exp());
    let err_at = |h: f64| -> f64 {
        let steps = (1.0 / h).round() as usize;
        let ee = ExpEuler::new(&[-1.0], 1, h);
        let nonlin = |u: &[f64], out: &mut [f64]| out[0] = u[0] * u[0];
        let mut u = vec![0.5f64];
        for _ in 0..steps {
            ee.step(&mut u, &nonlin);
        }
        (u[0] - exact(1.0)).abs()
    };
    let (e1, e2) = (err_at(0.01), err_at(0.005));
    let order = (e1 / e2).log2();
    assert!((order - 1.0).abs() < 0.25, "ETD1 order {order:.2}");
    log("expeuler-order", "pass", &format!("order={order:.2}"));
}

// ------------------------------------------------------------------ adaptive

fn ho_rhs(_t: f64, u: &[f64], out: &mut [f64]) {
    out[0] = u[1];
    out[1] = -u[0];
}

#[test]
fn rk45_accuracy_tracks_tolerance() {
    let pi = PiController::default();
    let mut st = AdaptiveState::new(0.0, &[1.0, 0.0], 0.1);
    rk45_adaptive(&mut st, &ho_rhs, 10.0, 1e-8, 1e-10, &pi, 100_000);
    assert!((st.t - 10.0).abs() < 1e-12, "did not reach t_end: {}", st.t);
    let err = (st.u[0] - 10.0f64.cos()).abs().max((st.u[1] + 10.0f64.sin()).abs());
    assert!(err < 1e-6, "RK45 error {err:.3e} at rtol=1e-8");
    assert!(st.accepted > 20 && st.accepted < 2_000, "step count off: {}", st.accepted);
    log(
        "rk45-accuracy",
        "pass",
        &format!("err={err:.3e} accepted={} rejected={}", st.accepted, st.rejected),
    );
}

#[test]
fn rk45_rejection_recovers_from_huge_h0() {
    let pi = PiController::default();
    let mut st = AdaptiveState::new(0.0, &[1.0, 0.0], 50.0); // absurd h₀
    rk45_adaptive(&mut st, &ho_rhs, 10.0, 1e-8, 1e-10, &pi, 100_000);
    assert!(st.rejected >= 1, "expected at least one rejection from h0=50");
    let err = (st.u[0] - 10.0f64.cos()).abs().max((st.u[1] + 10.0f64.sin()).abs());
    assert!(err < 1e-6, "post-rejection accuracy {err:.3e}");
    log("rk45-reject", "pass", &format!("rejected={} err={err:.3e}", st.rejected));
}

#[test]
fn rk45_resumable_split_run_bitwise() {
    // The P7 contract: state is COMPLETE. Interrupting after any number
    // of attempts (max_steps), checkpointing by clone, and resuming
    // must be bitwise-identical to the straight run — controller memory
    // (err_prev) and counters included.
    let pi = PiController::default();
    let mut straight = AdaptiveState::new(0.0, &[1.0, 0.0], 0.1);
    rk45_adaptive(&mut straight, &ho_rhs, 10.0, 1e-9, 1e-12, &pi, 100_000);
    for &cut in &[1usize, 7, 23, 61] {
        let mut first = AdaptiveState::new(0.0, &[1.0, 0.0], 0.1);
        rk45_adaptive(&mut first, &ho_rhs, 10.0, 1e-9, 1e-12, &pi, cut);
        let mut resumed = first.clone(); // checkpoint = clone
        rk45_adaptive(&mut resumed, &ho_rhs, 10.0, 1e-9, 1e-12, &pi, 100_000);
        assert_eq!(resumed.t.to_bits(), straight.t.to_bits(), "t differs at cut {cut}");
        assert_eq!(resumed.h.to_bits(), straight.h.to_bits(), "h differs at cut {cut}");
        assert_eq!(
            resumed.err_prev.to_bits(),
            straight.err_prev.to_bits(),
            "controller memory differs at cut {cut}"
        );
        for i in 0..2 {
            assert_eq!(resumed.u[i].to_bits(), straight.u[i].to_bits(), "u[{i}] at cut {cut}");
        }
        assert_eq!(
            (resumed.accepted, resumed.rejected),
            (straight.accepted, straight.rejected),
            "counters at cut {cut}"
        );
    }
    log("rk45-resume", "pass", "4 split points bitwise == straight run");
}

// ------------------------------------------------------------------- adjoint

#[test]
fn verlet_adjoint_gradcheck_vs_central_fd() {
    // F(q) = −q − 0.3‖q‖²q (symmetric Jacobian, as the adjoint's
    // Jᵀv = Jv shortcut requires); J = ½‖q_N‖² + ½‖p_N‖².
    let n = 3usize;
    let (h, steps) = (0.01f64, 64usize);
    let force = |q: &[f64], out: &mut [f64]| {
        let r2: f64 = q.iter().map(|x| x * x).sum();
        for i in 0..q.len() {
            out[i] = -q[i] - 0.3 * r2 * q[i];
        }
    };
    let force_jvp = |q: &[f64], v: &[f64], out: &mut [f64]| {
        let r2: f64 = q.iter().map(|x| x * x).sum();
        let qv: f64 = q.iter().zip(v).map(|(a, b)| a * b).sum();
        for i in 0..q.len() {
            out[i] = -v[i] - 0.3 * (r2 * v[i] + 2.0 * qv * q[i]);
        }
    };
    let q0 = vec![0.4f64, -0.7, 0.2];
    let p0 = vec![0.1f64, 0.3, -0.5];
    let cost = |q0: &[f64], p0: &[f64]| -> f64 {
        let (mut q, mut p) = (q0.to_vec(), p0.to_vec());
        let mut scratch = vec![0.0f64; n];
        for _ in 0..steps {
            verlet_step(&mut q, &mut p, h, &force, &mut scratch);
        }
        0.5 * (q.iter().map(|x| x * x).sum::<f64>() + p.iter().map(|x| x * x).sum::<f64>())
    };
    // Terminal cotangent: (q_N, p_N).
    let (mut qn, mut pn) = (q0.clone(), p0.clone());
    let mut scratch = vec![0.0f64; n];
    for _ in 0..steps {
        verlet_step(&mut qn, &mut pn, h, &force, &mut scratch);
    }
    let (bar_q0, bar_p0) = verlet_adjoint(&q0, &p0, h, steps, &force, &force_jvp, &qn, &pn);
    let eps = 1e-6;
    let mut worst = 0.0f64;
    for i in 0..n {
        let mut qp = q0.clone();
        qp[i] += eps;
        let mut qm = q0.clone();
        qm[i] -= eps;
        let fd = (cost(&qp, &p0) - cost(&qm, &p0)) / (2.0 * eps);
        worst = worst.max((fd - bar_q0[i]).abs() / fd.abs().max(1.0));
        let mut pp = p0.clone();
        pp[i] += eps;
        let mut pm = p0.clone();
        pm[i] -= eps;
        let fd = (cost(&q0, &pp) - cost(&q0, &pm)) / (2.0 * eps);
        worst = worst.max((fd - bar_p0[i]).abs() / fd.abs().max(1.0));
    }
    assert!(worst < 1e-7, "adjoint gradcheck worst rel err {worst:.3e}");
    log("verlet-adjoint", "pass", &format!("worst_rel={worst:.3e}"));
}

// --------------------------------------------------------------- golden hash

const GOLDEN_HASH: u64 = 0; // recorded on first run, then frozen

#[test]
fn time_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    // Verlet: 10k pendulum steps.
    let force = |q: &[f64], out: &mut [f64]| out[0] = -fs_math::det::sin(q[0]);
    let (mut q, mut p) = (vec![1.2f64], vec![0.3f64]);
    let mut scratch = vec![0.0f64; 1];
    for _ in 0..10_000 {
        verlet_step(&mut q, &mut p, 0.01, &force, &mut scratch);
    }
    feed(q[0]);
    feed(p[0]);
    // Lie: asymmetric rigid body, 10k steps.
    let mut quat = [1.0f64, 0.0, 0.0, 0.0];
    let mut w = [0.3f64, 1.5, -0.4];
    for _ in 0..10_000 {
        let (qn, wn) = rigid_body_step(quat, w, [1.0, 2.0, 3.0], 1e-3);
        quat = qn;
        w = wn;
    }
    for v in quat {
        feed(v);
    }
    for v in w {
        feed(v);
    }
    // Generalized-α: damped 2-DOF, 1k steps.
    let m = [2.0f64, 0.0, 0.0, 1.0];
    let c = [0.1f64, 0.0, 0.0, 0.1];
    let k = [4.0f64, -1.0, -1.0, 3.0];
    let ga = GeneralizedAlpha::new(&m, &c, &k, 2, 0.05, 0.8);
    let (mut gq, mut gv, mut ga_acc) = (vec![1.0f64, -0.5], vec![0.0f64; 2], vec![0.0f64; 2]);
    // Consistent initial acceleration a = M⁻¹(f − Cv − Kq), f = 0.
    ga_acc[0] = -(k[0] * gq[0] + k[1] * gq[1]) / m[0];
    ga_acc[1] = -(k[2] * gq[0] + k[3] * gq[1]) / m[3];
    let f = vec![0.0f64; 2];
    for _ in 0..1_000 {
        galpha_step(&ga, &mut gq, &mut gv, &mut ga_acc, &f);
    }
    for v in gq.iter().chain(gv.iter()) {
        feed(*v);
    }
    // IMEX: 2-D stiff/nonstiff mix, 500 steps.
    let im = Imex2::new(&[-50.0, 1.0, 0.0, -0.5], 2, 0.02);
    let nl = |u: &[f64], out: &mut [f64]| {
        out[0] = -u[0] * u[1];
        out[1] = u[0] * u[0];
    };
    let mut iu = vec![1.0f64, 0.2];
    for _ in 0..500 {
        imex2_step(&im, &mut iu, &nl);
    }
    feed(iu[0]);
    feed(iu[1]);
    // Exponential Euler: symmetric A, 500 steps.
    let ee = ExpEuler::new(&[-2.0, 0.7, 0.7, -3.0], 2, 0.02);
    let mut eu = vec![0.9f64, -0.6];
    for _ in 0..500 {
        ee.step(&mut eu, &nl);
    }
    feed(eu[0]);
    feed(eu[1]);
    // RK45 + PI controller: full trajectory state including controller
    // memory and counters (the step SEQUENCE is part of the contract).
    let pi = PiController::default();
    let mut st = AdaptiveState::new(0.0, &[1.0, 0.0], 0.1);
    rk45_adaptive(&mut st, &ho_rhs, 10.0, 1e-9, 1e-12, &pi, 100_000);
    feed(st.u[0]);
    feed(st.u[1]);
    feed(st.h);
    feed(st.err_prev);
    #[allow(clippy::cast_precision_loss)]
    feed(st.accepted as f64);
    #[allow(clippy::cast_precision_loss)]
    feed(st.rejected as f64);
    log("time-golden", "info", &format!("{acc:#018x}"));
    assert_eq!(
        acc, GOLDEN_HASH,
        "fs-time bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}
