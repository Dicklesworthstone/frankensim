//! fs-bem conformance battery (bead tfz.20).
//!
//! - bem-001 G0: the uniform-sphere Gauss identity — the assembled
//!   Neumann operator applied to ones gives −1 at every centroid
//!   (sign conventions and self-terms cannot silently drift).
//! - bem-002 G2: sphere exterior flow — surface speed vs the analytic
//!   1.5·U·sinθ at two refinements (order ledgered), Cp curve.
//! - bem-003: the FMM path — matvec and transpose match dense to
//!   interpolation tolerance; GMRES(FMM) reproduces the dense-LU solution;
//!   iterations ledgered.
//! - bem-004: 2D Hess–Smith — thin-airfoil lift slope band, Cp sanity
//!   (stagnation near LE, smooth TE), and the ADJOINT dCl/dα against
//!   central FD.
//! - bem-005: impulsive-start free wake — Wagner-like transient toward
//!   the pressure-derived screening circulation asymptote, Kelvin bookkeeping,
//!   bounded stable roll-up,
//!   bitwise determinism.

use fs_bem::panel3d::{SpherePanels, solve_exterior, surface_velocity};
use fs_bem::{WakeSim, naca4_symmetric, panel2d};
use fs_la::factor::lu;
use std::fmt::Write as _;

fn verdict(name: &str, pass: bool, details: &str) {
    println!(
        "{{\"test\":\"{name}\",\"model\":\"inviscid-screening\",\"verdict\":\"{}\",{details}}}",
        if pass { "pass" } else { "fail" }
    );
    assert!(pass, "{name} failed: {details}");
}

// ------------------------------------------------------------------ bem-001

#[test]
fn bem_001_gauss_identity() {
    let panels = SpherePanels::icosphere(1.0, 2).expect("valid fixture");
    let n = panels.centroids().len();
    let a = panels.dense_matrix().expect("admitted dense fixture");
    let mut worst = 0.0f64;
    for i in 0..n {
        let row: f64 = (0..n).map(|j| a[i * n + j]).sum();
        worst = worst.max((row - (-1.0)).abs());
    }
    verdict(
        "bem-001",
        worst < 0.05,
        &format!(
            "\"detail\":\"uniform source sheet on the sphere: row action = -1 (Gauss)\",\
             \"panels\":{n},\"worst_dev\":{worst:.3e}"
        ),
    );
}

// ------------------------------------------------------------------ bem-002

fn sphere_speed_error(subdivisions: u32) -> (f64, usize) {
    let panels = SpherePanels::icosphere(1.0, subdivisions).expect("valid fixture");
    let n = panels.centroids().len();
    let a = panels.dense_matrix().expect("admitted dense fixture");
    let u_inf = [1.0, 0.0, 0.0];
    let mut rhs: Vec<f64> = (0..n)
        .map(|i| {
            -(u_inf[0] * panels.normals()[i][0]
                + u_inf[1] * panels.normals()[i][1]
                + u_inf[2] * panels.normals()[i][2])
        })
        .collect();
    let f = lu(&a, n).expect("dense solve");
    f.solve(&mut rhs);
    let vel = surface_velocity(&panels, &rhs, u_inf, 6).expect("valid velocity request");
    let mut err_sum = 0.0;
    for (i, got_vel) in vel.iter().enumerate() {
        let c = panels.centroids()[i];
        let r = (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt();
        let sin_theta = (1.0 - (c[0] / r) * (c[0] / r)).max(0.0).sqrt();
        let want = 1.5 * sin_theta;
        let got =
            (got_vel[0] * got_vel[0] + got_vel[1] * got_vel[1] + got_vel[2] * got_vel[2]).sqrt();
        err_sum += (got - want).abs();
    }
    #[allow(clippy::cast_precision_loss)]
    (err_sum / n as f64, n)
}

#[test]
fn bem_002_sphere_analytic() {
    let (e2, n2) = sphere_speed_error(2);
    let (e3, n3) = sphere_speed_error(3);
    let order = (e2 / e3).log2() / ((n3 as f64 / n2 as f64).log2() / 2.0);
    let pass = e3 < 0.03 && e3 < e2;
    verdict(
        "bem-002",
        pass,
        &format!(
            "\"detail\":\"sphere surface speed vs 1.5 U sin(theta)\",\
             \"mean_abs_err\":[{{\"panels\":{n2},\"err\":{e2:.4}}},{{\"panels\":{n3},\"err\":{e3:.4}}}],\
             \"order_per_h\":{order:.2}"
        ),
    );
}

// ------------------------------------------------------------------ bem-003

#[test]
fn bem_003_fmm_path_matches_dense() {
    let panels = SpherePanels::icosphere(1.0, 3).expect("valid fixture");
    let n = panels.centroids().len();
    // Matvec consistency on a deterministic test vector.
    #[allow(clippy::cast_precision_loss)]
    let x: Vec<f64> = (0..n).map(|i| ((i as f64) * 0.37).sin()).collect();
    let a = panels.dense_matrix().expect("admitted dense fixture");
    let mut dense = vec![0.0f64; n];
    for i in 0..n {
        for j in 0..n {
            dense[i] += a[i * n + j] * x[j];
        }
    }
    let fast = panels.fmm_matvec(&x, 6).expect("admitted FMM request");
    let scale = dense.iter().map(|v| v * v).sum::<f64>().sqrt();
    let dev = dense
        .iter()
        .zip(&fast)
        .map(|(d, f)| (d - f) * (d - f))
        .sum::<f64>()
        .sqrt()
        / scale;
    let mut dense_t = vec![0.0f64; n];
    for i in 0..n {
        for j in 0..n {
            dense_t[j] += a[i * n + j] * x[i];
        }
    }
    let fast_t = panels
        .fmm_transpose_matvec(&x, 6)
        .expect("admitted transpose request");
    let scale_t = dense_t.iter().map(|v| v * v).sum::<f64>().sqrt();
    let dev_t = dense_t
        .iter()
        .zip(&fast_t)
        .map(|(d, f)| (d - f) * (d - f))
        .sum::<f64>()
        .sqrt()
        / scale_t;
    // GMRES(FMM) vs dense LU.
    let u_inf = [1.0, 0.0, 0.0];
    let exterior = solve_exterior(&panels, u_inf, 6, 1e-8).expect("fixture must converge");
    let sigma_fmm = exterior.sigma;
    let iters = exterior.report.iters;
    let rr = exterior.report.rel_residual;
    let mut sigma_dense: Vec<f64> = (0..n)
        .map(|i| -(u_inf[0] * panels.normals()[i][0]))
        .collect();
    let f = lu(&a, n).expect("dense solve");
    f.solve(&mut sigma_dense);
    let sscale = sigma_dense.iter().map(|v| v * v).sum::<f64>().sqrt();
    let sdev = sigma_dense
        .iter()
        .zip(&sigma_fmm)
        .map(|(d, f)| (d - f) * (d - f))
        .sum::<f64>()
        .sqrt()
        / sscale;
    let pass = dev < 1e-4 && dev_t < 1e-4 && sdev < 1e-3 && rr < 1e-6;
    verdict(
        "bem-003",
        pass,
        &format!(
            "\"detail\":\"FMM matvec/transpose + GMRES(FMM) vs dense oracle, 1280 panels\",\
             \"matvec_rel\":{dev:.3e},\"transpose_rel\":{dev_t:.3e},\"solution_rel\":{sdev:.3e},\
             \"gmres_iters\":{iters},\"gmres_residual\":{rr:.3e}"
        ),
    );
}

// ------------------------------------------------------------------ bem-004

#[test]
fn bem_004_hess_smith_kutta_and_adjoint() {
    let foil = naca4_symmetric(0.10, 120).expect("valid NACA fixture");
    // Lift slope from ±2°.
    let a2 = 2.0f64.to_radians();
    let cl_p = panel2d::solve(&foil, a2).expect("valid solve").cl;
    let cl_m = panel2d::solve(&foil, -a2).expect("valid solve").cl;
    let slope = (cl_p - cl_m) / (2.0 * a2);
    let thin = 2.0 * std::f64::consts::PI;
    // Thickness raises the inviscid slope: 2π(1 + 0.77 t) is the
    // classic correction — gate against BOTH bands.
    let thick_corrected = thin * (1.0 + 0.77 * 0.10);
    let slope_rel = (slope - thick_corrected).abs() / thick_corrected;
    // Cp sanity at 4°: stagnation near the LE (Cp → 1), TE speeds
    // matched (Kutta).
    let sol = panel2d::solve(&foil, 4.0f64.to_radians()).expect("valid solve");
    let cp_max = sol
        .vt
        .iter()
        .map(|v| 1.0 - v * v)
        .fold(f64::NEG_INFINITY, f64::max);
    let n = sol.vt.len();
    let kutta_dev = (sol.vt[0] + sol.vt[n - 1]).abs();
    // Adjoint gate.
    let alpha0 = 3.0f64.to_radians();
    let adj = panel2d::dcl_dalpha_adjoint(&foil, alpha0).expect("valid adjoint");
    let h = 1e-5;
    let fd = (panel2d::solve(&foil, alpha0 + h).expect("valid solve").cl
        - panel2d::solve(&foil, alpha0 - h).expect("valid solve").cl)
        / (2.0 * h);
    let adj_rel = (adj - fd).abs() / fd.abs().max(1e-30);
    let pass = slope_rel < 0.05
        && slope > thin
        && cp_max > 0.95
        && cp_max < 1.05
        && kutta_dev < 0.05
        && adj_rel < 1e-6;
    verdict(
        "bem-004",
        pass,
        &format!(
            "\"detail\":\"NACA0010 Hess-Smith: thin-airfoil slope band, Cp sanity, adjoint gate \
             (inviscid screening honesty label)\",\
             \"dcl_dalpha\":{slope:.4},\"thin_airfoil\":{thin:.4},\"thickness_corrected\":{thick_corrected:.4},\"slope_rel\":{slope_rel:.3},\
             \"cp_stagnation\":{cp_max:.3},\"kutta_dev\":{kutta_dev:.3e},\
             \"adjoint\":{adj:.5},\"fd\":{fd:.5},\"adjoint_rel\":{adj_rel:.3e}"
        ),
    );
}

// ------------------------------------------------------------------ bem-005

#[test]
fn bem_005_impulsive_start_free_wake() {
    let run = || {
        let foil = naca4_symmetric(0.08, 80).expect("valid NACA fixture");
        let mut sim = WakeSim::new(&foil, 5.0f64.to_radians(), 0.05, 0.05).expect("valid wake");
        let mut kelvin_worst = 0.0f64;
        for _ in 0..200 {
            sim.step().expect("bounded finite step");
            let total = sim.history().last().expect("step").bound
                + sim.wake().iter().map(|vortex| vortex.gamma).sum::<f64>();
            kelvin_worst = kelvin_worst.max(total.abs());
        }
        (sim, kelvin_worst)
    };
    let (sim, kelvin_worst) = run();
    let (sim_b, kelvin_worst_b) = run();
    let history_deterministic = sim.history().iter().zip(sim_b.history()).all(|(a, b)| {
        a.t.to_bits() == b.t.to_bits()
            && a.bound.to_bits() == b.bound.to_bits()
            && a.vortices == b.vortices
            && a.peak_speed.to_bits() == b.peak_speed.to_bits()
    });
    let wake_deterministic = sim.wake().iter().zip(sim_b.wake()).all(|(a, b)| {
        a.pos[0].to_bits() == b.pos[0].to_bits()
            && a.pos[1].to_bits() == b.pos[1].to_bits()
            && a.gamma.to_bits() == b.gamma.to_bits()
    });
    let trace = sim.trace_json(40).expect("valid trace");
    let trace_b = sim_b.trace_json(40).expect("valid trace");
    let deterministic = sim.history().len() == sim_b.history().len()
        && sim.wake().len() == sim_b.wake().len()
        && history_deterministic
        && wake_deterministic
        && kelvin_worst.to_bits() == kelvin_worst_b.to_bits()
        && trace == trace_b;
    let steady = sim.steady_circulation();
    let first = sim.history()[0].bound / steady;
    let last = sim.history().last().expect("steps").bound / steady;
    // Stability: all wake positions finite and bounded; peak speeds
    // bounded; growth monotone-ish (no oscillatory blowup).
    let bounded = sim.wake().iter().all(|w| {
        w.pos[0].is_finite()
            && w.pos[1].is_finite()
            && w.pos[0].abs() < 50.0
            && w.pos[1].abs() < 50.0
    });
    let peak = sim
        .history()
        .iter()
        .map(|s| s.peak_speed)
        .fold(0.0f64, f64::max);
    // The lumped starting vortex passing the control point causes
    // real early-transient dips (ledgered); the QUALITATIVE gate is
    // the coarse-grained trend: stride-40 samples nondecreasing.
    let mut monotone_violations = 0usize;
    let mut backslide = 0.0f64;
    for w in sim.history().windows(2) {
        if w[1].bound < w[0].bound - 1e-3 * steady {
            monotone_violations += 1;
            backslide += w[0].bound - w[1].bound;
        }
    }
    let coarse: Vec<f64> = sim.history().iter().step_by(40).map(|s| s.bound).collect();
    let coarse_monotone = coarse.windows(2).all(|w| w[1] >= w[0]);
    let wagner = (0.3..=0.7).contains(&first);
    let asymptote = last > 0.9 && last < 1.05;
    let kelvin = kelvin_worst < 1e-12;
    let pass =
        wagner && asymptote && bounded && peak < 5.0 && coarse_monotone && kelvin && deterministic;
    let mut tail = String::new();
    let _ = write!(tail, "{trace}");
    verdict(
        "bem-005",
        pass,
        &format!(
            "\"detail\":\"impulsive start: Wagner-like transient, stable roll-up, determinism\",\
             \"first_over_steady\":{first:.3},\"last_over_steady\":{last:.3},\
             \"peak_speed\":{peak:.3},\"backslide\":{backslide:.3e},\"coarse_monotone\":{coarse_monotone},\"wagner\":{wagner},\"asymptote\":{asymptote},\"bounded\":{bounded},\"kelvin_worst\":{kelvin_worst:.3e},\"kelvin\":{kelvin},\"monotone_violations\":{monotone_violations},\"vortices\":{},\"deterministic\":{deterministic},\
             \"trace\":{tail}",
            sim.wake().len()
        ),
    );
}

// ------------------------------------------------------------------ totality

#[test]
fn bem_rejects_invalid_geometry_work_and_trace_requests() {
    assert!(naca4_symmetric(f64::NAN, 80).is_err());
    assert!(naca4_symmetric(0.1, usize::MAX).is_err());
    assert!(panel2d::Airfoil2d::new(vec![[0.0, 0.0]; 4]).is_err());

    let foil = naca4_symmetric(0.1, 80).expect("valid fixture");
    assert!(panel2d::solve(&foil, f64::NAN).is_err());
    assert!(WakeSim::new(&foil, 0.1, 0.0, 0.05).is_err());
    let sim = WakeSim::new(&foil, 0.1, 0.05, 0.05).expect("valid wake");
    assert!(sim.trace_json(0).is_err());

    assert!(SpherePanels::icosphere(1.0, u32::MAX).is_err());
    let too_dense = SpherePanels::icosphere(1.0, 4).expect("bounded panelization");
    assert!(too_dense.dense_matrix().is_err());
    let small = SpherePanels::icosphere(1.0, 1).expect("valid fixture");
    assert!(
        small
            .fmm_matvec(&vec![1.0; small.centroids().len()], usize::MAX)
            .is_err()
    );
    assert!(solve_exterior(&small, [1.0, 0.0, 0.0], 4, 0.0).is_err());

    let zero = solve_exterior(&small, [0.0; 3], 4, 1e-8).expect("zero flow is exact");
    assert!(zero.report.converged);
    assert_eq!(zero.report.rel_residual, 0.0);
    assert!(zero.sigma.iter().all(|value| *value == 0.0));
}

#[test]
fn unconverged_exterior_iterate_is_not_published_as_success() {
    let panels = SpherePanels::icosphere(1.0, 1).expect("valid fixture");
    let error = solve_exterior(&panels, [1.0, 0.0, 0.0], 4, f64::MIN_POSITIVE)
        .expect_err("subnormal tolerance must not be reported as converged");
    match error {
        fs_bem::panel3d::ExteriorSolveError::NotConverged { sigma, report } => {
            assert!(!report.converged);
            assert_eq!(sigma.len(), panels.centroids().len());
            assert!(report.rel_residual.is_finite());
        }
        other => panic!("expected an unconverged report, got {other}"),
    }
}
