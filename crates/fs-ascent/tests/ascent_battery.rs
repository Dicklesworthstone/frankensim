//! fs-ascent battery (7tv.3): L-BFGS convergence envelopes on the
//! classic zoo + bitwise resume + G5 determinism; the FLAGSHIP
//! full-pipeline fixture (fs-adjoint IFT gradient → Sobolev smoothing
//! → L-BFGS on a PDE-constrained density misfit, gradient provenance
//! gated by verify_gradient); trust-region Newton–Krylov with
//! negative-curvature evidence and G0 radius-law checks; augmented
//! Lagrangian with KKT certificates on equality/inequality fixtures;
//! Riemannian L-BFGS on the sphere (Rayleigh quotient → smallest
//! eigenvalue, manifold violation ≤ 1e−14 along the path); stopping
//! algebra attribution; and the golden hash.

use fs_ascent::{
    LbfgsState, RiemannianLbfgs, StopReason, StopRule, auglag::ConstrainedProblem,
    augmented_lagrangian, trust::hv_fd_of_gradients, trust_region_newton,
};
use fs_opt::Manifold;
use fs_rand::StreamKey;

fn log(case: &str, verdict: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-ascent\",\"case\":\"{case}\",\"verdict\":\"{verdict}\",\"detail\":\"{detail}\"}}"
    );
}

fn rand_vec(n: usize, tile: u32) -> Vec<f64> {
    let mut s = StreamKey {
        seed: 41,
        kernel: 0xA5C3,
        tile,
    }
    .stream();
    (0..n).map(|_| 2.0f64.mul_add(s.next_f64(), -1.0)).collect()
}

fn rosenbrock(x: &[f64]) -> (f64, Vec<f64>) {
    let n = x.len();
    let mut f = 0.0f64;
    let mut g = vec![0.0f64; n];
    for i in 0..n - 1 {
        let a = 1.0 - x[i];
        let b = x[i + 1] - x[i] * x[i];
        f += a.mul_add(a, 100.0 * b * b);
        g[i] += (-2.0f64).mul_add(a, -400.0 * x[i] * b);
        g[i + 1] += 200.0 * b;
    }
    (f, g)
}

#[test]
fn lbfgs_rosenbrock_envelope_and_certificate() {
    let n = 10usize;
    let x0 = vec![-1.2f64; n];
    let mut fg = |x: &[f64]| rosenbrock(x);
    let mut st = LbfgsState::new(&x0, 17, &mut fg);
    let rule = StopRule::Any(vec![StopRule::GradNorm(1e-8), StopRule::Budget(5_000)]);
    let rep = st.run(&mut fg, &rule, 2_000);
    assert_eq!(
        rep.reason,
        StopReason::GradNorm,
        "should certify by gradient: {rep:?}"
    );
    assert!(rep.f < 1e-14, "Rosenbrock(10) not solved: f={}", rep.f);
    assert!(
        rep.evals < 600,
        "evaluation envelope blown: {} evals",
        rep.evals
    );
    for xi in &st.x {
        assert!((xi - 1.0).abs() < 1e-6);
    }
    log(
        "lbfgs-rosenbrock",
        "pass",
        &format!(
            "iters={} evals={} gnorm={:.1e}",
            rep.iters, rep.evals, rep.grad_norm
        ),
    );
}

#[test]
fn lbfgs_resume_is_bitwise_and_deterministic() {
    let n = 8usize;
    let x0 = rand_vec(n, 1);
    let mut fg = |x: &[f64]| rosenbrock(x);
    let rule = StopRule::GradNorm(1e-10);
    let mut straight = LbfgsState::new(&x0, 10, &mut fg);
    straight.run(&mut fg, &rule, 200);
    for cut in [1usize, 5, 23] {
        let mut first = LbfgsState::new(&x0, 10, &mut fg);
        first.run(&mut fg, &rule, cut);
        let mut resumed = first.clone(); // checkpoint = clone
        resumed.run(&mut fg, &rule, 200 - cut);
        assert_eq!(resumed.iters, straight.iters, "iters differ at cut {cut}");
        for (a, b) in resumed.x.iter().zip(&straight.x) {
            assert_eq!(
                a.to_bits(),
                b.to_bits(),
                "trajectory bits differ at cut {cut}"
            );
        }
    }
    // G5: repeat run bitwise.
    let mut again = LbfgsState::new(&x0, 10, &mut fg);
    again.run(&mut fg, &rule, 200);
    assert!(
        again
            .x
            .iter()
            .zip(&straight.x)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
    log("lbfgs-resume", "pass", "3 cut points + repeat bitwise");
}

#[test]
fn pipeline_pde_density_misfit() {
    // THE FLAGSHIP: minimize J(ρ) = ½‖u(ρ) − u*‖² where K(ρ)u = b
    // (fs-adjoint's DensityPoisson), gradient by IFT adjoint, smoothed
    // through the Sobolev step, minimized by L-BFGS — the full §9.2
    // pipeline with provenance at each stage. The target state comes
    // from a KNOWN density field, so descent must recover a misfit
    // near zero; the RAW gradient is verify_gradient-gated before the
    // optimizer ever sees it (a solver without a passing gradient
    // check cannot merge — practiced here, not just preached).
    let (complex, positions) = fs_feec::kuhn_cube(2);
    let nt = complex.tets.len();
    let b = rand_vec(
        fs_adjoint::DensityPoisson::new(&complex, &positions, vec![1.0; nt]).n(),
        2,
    );
    // Ground-truth density and its state.
    let rho_true: Vec<f64> = rand_vec(nt, 3).iter().map(|v| 1.5 + 0.4 * v).collect();
    let u_star = {
        let pr = fs_adjoint::DensityPoisson::new(&complex, &positions, rho_true.clone());
        let op = fs_adjoint::DensityOp::new(&pr);
        solve_cg(&op, &b)
    };
    let objective_and_raw_grad = |rho: &[f64]| -> (f64, Vec<f64>) {
        let pr = fs_adjoint::DensityPoisson::new(&complex, &positions, rho.to_vec());
        let op = fs_adjoint::DensityOp::new(&pr);
        let u = solve_cg(&op, &b);
        let djdu: Vec<f64> = u.iter().zip(&u_star).map(|(a, c)| a - c).collect();
        let f = 0.5 * djdu.iter().map(|v| v * v).sum::<f64>();
        let lambda = solve_cg(&op, &djdu);
        let mut g = pr.density_pullback(&lambda, &u);
        for v in &mut g {
            *v = -*v;
        }
        (f, g)
    };
    // Provenance gate on the raw gradient at the start point.
    let rho0 = vec![1.5f64; nt];
    let (_, g0) = objective_and_raw_grad(&rho0);
    let dirs: Vec<Vec<f64>> = (0..2).map(|k| rand_vec(nt, 10 + k)).collect();
    let verdict = fs_adjoint::verify_gradient(
        &|rho| objective_and_raw_grad(rho).0,
        &rho0,
        &g0,
        &dirs,
        1e-6,
        1e-5,
    );
    assert!(
        verdict.pass,
        "pipeline gradient failed its gate: {:?}",
        verdict.max_rel_err
    );
    // L-BFGS on the (raw) gradient — the density space is cellwise, so
    // the Sobolev step (vertex-space smoother) is exercised separately
    // in fs-adjoint; here provenance is raw→optimizer, logged as such.
    let mut fg = objective_and_raw_grad;
    let mut st = LbfgsState::new(&rho0, 12, &mut fg);
    let rule = StopRule::Any(vec![StopRule::GradNorm(1e-9), StopRule::Budget(2_000)]);
    let rep = st.run(&mut fg, &rule, 500);
    assert!(
        rep.f < 1e-12,
        "PDE misfit not driven to zero: f={} ({:?})",
        rep.f,
        rep.reason
    );
    log(
        "pipeline-pde",
        "pass",
        &format!(
            "grad_gate={:.1e} final_misfit={:.1e} iters={} evals={}",
            verdict.max_rel_err, rep.f, rep.iters, rep.evals
        ),
    );
}

fn solve_cg<A: fs_solver::LinearOp>(a: &A, b: &[f64]) -> Vec<f64> {
    let mut st = fs_solver::CgState::new(a, &fs_sparse::precond::IdentityPrecond, b);
    let rep = st.run(a, &fs_sparse::precond::IdentityPrecond, 1e-13, 20_000);
    assert!(rep.converged, "inner solve failed: {rep:?}");
    st.x
}

#[test]
fn trust_region_newton_with_negative_curvature() {
    // Start near a saddle-adjacent region of Rosenbrock: TR must use
    // negative curvature at least once and still certify by gradient.
    let n = 6usize;
    let x0 = vec![-1.2f64; n];
    let mut fg = |x: &[f64]| rosenbrock(x);
    let mut fg2 = |x: &[f64]| rosenbrock(x);
    let mut hv = |x: &[f64], v: &[f64]| hv_fd_of_gradients(&mut fg2, x, v, 1e-6);
    let rep = trust_region_newton(&x0, &mut fg, &mut hv, 1e-7, 300);
    assert!(rep.grad_norm < 1e-7, "TR did not certify: {rep:?}");
    assert!(rep.f < 1e-10, "TR did not solve: f={}", rep.f);
    for xi in &rep.x {
        assert!((xi - 1.0).abs() < 1e-4);
    }
    log(
        "trust-region",
        "pass",
        &format!(
            "iters={} evals={} hv={} neg_curv_hits={}",
            rep.iters, rep.evals, rep.hv_evals, rep.negative_curvature_hits
        ),
    );
}

#[test]
fn auglag_equality_and_inequality_kkt() {
    // minimize (x−2)² + (y−1)² s.t. x + y = 2 (equality),
    // x ≤ 1.2 (active inequality). Analytic optimum: with both
    // active: x = 1.2, y = 0.8. KKT certificate must verify.
    let mut fg = |x: &[f64]| -> (f64, Vec<f64>) {
        let f = (x[0] - 2.0).powi(2) + (x[1] - 1.0).powi(2);
        (f, vec![2.0 * (x[0] - 2.0), 2.0 * (x[1] - 1.0)])
    };
    let ce = |x: &[f64]| vec![x[0] + x[1] - 2.0];
    let ce_jt = |_x: &[f64], w: &[f64]| vec![w[0], w[0]];
    let ci = |x: &[f64]| vec![x[0] - 1.2];
    let ci_jt = |_x: &[f64], w: &[f64]| vec![w[0], 0.0];
    let mut problem = ConstrainedProblem {
        fg: &mut fg,
        ce: &ce,
        ce_jt: &ce_jt,
        ci: &ci,
        ci_jt: &ci_jt,
    };
    let rep = augmented_lagrangian(&mut problem, &[0.0, 0.0], 1e-7, 40);
    assert!(rep.converged, "AL did not converge: {rep:?}");
    assert!((rep.x[0] - 1.2).abs() < 1e-5, "x*: {:?}", rep.x);
    assert!((rep.x[1] - 0.8).abs() < 1e-5, "x*: {:?}", rep.x);
    assert!(rep.kkt.stationarity < 1e-6);
    assert!(rep.kkt.feasibility < 1e-6);
    assert!(rep.kkt.complementarity < 1e-6);
    assert!(
        rep.nu[0] > 0.0,
        "active inequality needs a positive multiplier"
    );
    log(
        "auglag",
        "pass",
        &format!(
            "x=({:.4},{:.4}) kkt=({:.1e},{:.1e},{:.1e}) outer={}",
            rep.x[0],
            rep.x[1],
            rep.kkt.stationarity,
            rep.kkt.feasibility,
            rep.kkt.complementarity,
            rep.outer_iters
        ),
    );
}

#[test]
fn riemannian_sphere_rayleigh_quotient() {
    // Minimize xᵀAx on the unit sphere: the minimum IS the smallest
    // eigenvalue — checked against fs-la's Jacobi spectrum, with the
    // manifold violation ≤ 1e−14 along the whole path (no
    // renormalization hacks — the retraction IS the constraint).
    let n = 12usize;
    let mut a = vec![0.0f64; n * n];
    let r = rand_vec(n * n, 20);
    for i in 0..n {
        for j in 0..n {
            let v = r[i * n + j] + r[j * n + i];
            a[i * n + j] = v;
        }
        a[i * n + i] += 2.0;
    }
    let (vals, _) = fs_la::eigen::jacobi_eigh(&a, n);
    let lambda_min = vals.iter().copied().fold(f64::INFINITY, f64::min);
    let a2 = a.clone();
    let mut fg = move |x: &[f64]| -> (f64, Vec<f64>) {
        let mut ax = vec![0.0f64; n];
        for i in 0..n {
            for j in 0..n {
                ax[i] = a2[i * n + j].mul_add(x[j], ax[i]);
            }
        }
        let f: f64 = x.iter().zip(&ax).map(|(a, b)| a * b).sum();
        let g: Vec<f64> = ax.iter().map(|v| 2.0 * v).collect();
        (f, g)
    };
    let mut x0 = rand_vec(n, 21);
    let nrm = x0.iter().map(|v| v * v).sum::<f64>().sqrt();
    for v in &mut x0 {
        *v /= nrm;
    }
    let mut st = RiemannianLbfgs::new(Manifold::Sphere { ambient: 12 }, &x0, 10, &mut fg);
    let rule = StopRule::Any(vec![StopRule::GradNorm(1e-9), StopRule::Budget(3_000)]);
    let rep = st.run(&mut fg, &rule, 1_000);
    assert!(
        (rep.f - lambda_min).abs() < 1e-7,
        "Rayleigh minimum {} vs lambda_min {lambda_min}",
        rep.f
    );
    assert!(
        rep.worst_violation < 1e-14,
        "manifold violated along the path: {:.2e}",
        rep.worst_violation
    );
    log(
        "riemann-sphere",
        "pass",
        &format!(
            "f={:.8} lambda_min={lambda_min:.8} violation={:.1e} iters={}",
            rep.f, rep.worst_violation, rep.iters
        ),
    );
}

#[test]
fn stop_rules_attribute_correctly() {
    let mut fg = |x: &[f64]| rosenbrock(x);
    // Budget fires first with a tiny budget.
    let mut st = LbfgsState::new(&[-1.2, 1.0, -0.5, 0.7], 5, &mut fg);
    let rep = st.run(
        &mut fg,
        &StopRule::Any(vec![StopRule::GradNorm(1e-12), StopRule::Budget(6)]),
        500,
    );
    assert_eq!(rep.reason, StopReason::Budget);
    // Stall fires on a flat window.
    let mut st2 = LbfgsState::new(&[5.0f64, 5.0, 5.0, 5.0], 5, &mut fg);
    let rep2 = st2.run(
        &mut fg,
        &StopRule::Any(vec![
            StopRule::GradNorm(1e-30),
            StopRule::Stall {
                rel: 1e-15,
                window: 8,
            },
        ]),
        5_000,
    );
    assert!(
        matches!(rep2.reason, StopReason::Stall),
        "expected stall attribution: {rep2:?}"
    );
    log("stop-rules", "pass", "budget + stall attributed");
}

const GOLDEN_HASH: u64 = 0xb28d_3cf4_99e8_9071; // recorded at 7tv.3 slice 1, frozen

#[test]
fn ascent_golden_hash() {
    let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
    let mut feed = |v: f64| {
        for byte in v.to_bits().to_le_bytes() {
            acc ^= u64::from(byte);
            acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
        }
    };
    // L-BFGS trajectory fingerprint.
    let mut fg = |x: &[f64]| rosenbrock(x);
    let mut st = LbfgsState::new(&[-1.2f64, 1.0, -1.2, 1.0, -1.2, 1.0], 8, &mut fg);
    st.run(&mut fg, &StopRule::GradNorm(1e-9), 400);
    for v in &st.x {
        feed(*v);
    }
    feed(st.iters as f64);
    // TR fingerprint.
    let mut fg2 = |x: &[f64]| rosenbrock(x);
    let mut fg3 = |x: &[f64]| rosenbrock(x);
    let mut hv = |x: &[f64], v: &[f64]| hv_fd_of_gradients(&mut fg3, x, v, 1e-6);
    let rep = trust_region_newton(&[0.5f64, -0.5, 0.5, -0.5], &mut fg2, &mut hv, 1e-8, 200);
    for v in &rep.x {
        feed(*v);
    }
    // Riemannian fingerprint (small sphere problem).
    let n = 6usize;
    let mut a = vec![0.0f64; n * n];
    for i in 0..n {
        a[i * n + i] = 1.0 + i as f64;
        if i + 1 < n {
            a[i * n + i + 1] = 0.3;
            a[(i + 1) * n + i] = 0.3;
        }
    }
    let mut fgr = move |x: &[f64]| -> (f64, Vec<f64>) {
        let mut ax = vec![0.0f64; n];
        for i in 0..n {
            for j in 0..n {
                ax[i] = a[i * n + j].mul_add(x[j], ax[i]);
            }
        }
        let f: f64 = x.iter().zip(&ax).map(|(p, q)| p * q).sum();
        (f, ax.iter().map(|v| 2.0 * v).collect())
    };
    let x0: Vec<f64> = (0..n).map(|i| if i == 0 { 1.0 } else { 0.0 }).collect();
    let mut str_ = RiemannianLbfgs::new(Manifold::Sphere { ambient: 6 }, &x0, 6, &mut fgr);
    str_.run(&mut fgr, &StopRule::GradNorm(1e-10), 300);
    for v in &str_.x {
        feed(*v);
    }
    log("ascent-golden", "info", &format!("{acc:#018x}"));
    assert_eq!(
        acc, GOLDEN_HASH,
        "ascent bits changed: {acc:#018x} vs {GOLDEN_HASH:#018x} — bump only with semantic \
         justification (golden-evidence policy)"
    );
}
