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
//! - bem-006 G2: NACA 0012 pre-stall lift slope against the corrected
//!   NASA TM-4074 force table, with an explicit ten-degree admission boundary.

use fs_bem::panel3d::{SpherePanels, solve_exterior, surface_velocity};
use fs_bem::{
    BemError, NACA0012_PRESTALL_MAX_ALPHA_RAD, WakeSim, naca4_symmetric, panel2d,
    solve_naca0012_prestall,
};
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

// ------------------------------------------------------------------ bem-006

#[allow(clippy::cast_precision_loss)]
fn least_squares_slope<const N: usize>(x: &[f64; N], y: &[f64; N]) -> f64 {
    let count = N as f64;
    let mean_x = x.iter().sum::<f64>() / count;
    let mean_y = y.iter().sum::<f64>() / count;
    let covariance = x
        .iter()
        .zip(y)
        .map(|(x, y)| (x - mean_x) * (y - mean_y))
        .sum::<f64>();
    let variance = x.iter().map(|x| (x - mean_x).powi(2)).sum::<f64>();
    covariance / variance
}

#[test]
fn bem_006_naca0012_ladson_prestall_envelope() {
    // Charles L. Ladson, NASA TM-4074 (1988), table I: M=0.15,
    // Re=5.97e6, free transition. The report states that all tabulated force
    // coefficients include the standard low-speed tunnel-wall correction
    // (about two percent). The seven linear-range rows below are copied
    // verbatim; their least-squares slope is the independent experimental
    // reference rather than a thin-airfoil formula.
    // https://ntrs.nasa.gov/api/citations/19880019495/downloads/19880019495.pdf
    // SHA-256: 8e466706cbdf54b3c778ea2b089c4f52d87686bf7f6b1bd10d224b42c2d06902
    const LADSON_ALPHA_DEG: [f64; 7] = [-4.05, -2.00, 0.05, 1.98, 4.18, 6.20, 8.22];
    const LADSON_CL: [f64; 7] = [-0.4280, -0.2150, 0.0040, 0.2080, 0.4520, 0.6630, 0.8800];
    let ladson_alpha_rad = LADSON_ALPHA_DEG.map(f64::to_radians);
    let measured_slope = least_squares_slope(&ladson_alpha_rad, &LADSON_CL);

    let model_alpha_rad = [-8.0f64, -4.0, 0.0, 4.0, 8.0].map(f64::to_radians);
    let model_cl = model_alpha_rad.map(|alpha| {
        solve_naca0012_prestall(120, alpha)
            .expect("linear-range NACA 0012 request must be admitted")
            .cl
    });
    let model_slope = least_squares_slope(&model_alpha_rad, &model_cl);
    let relative_bias = model_slope / measured_slope - 1.0;
    let odd_symmetry_worst = (model_cl[0] + model_cl[4])
        .abs()
        .max((model_cl[1] + model_cl[3]).abs())
        .max(model_cl[2].abs());

    // TM-4074 reports that inviscid theory overpredicts its measured slope.
    // This gate therefore admits a one-sided, explicitly screening-grade
    // discrepancy of at most 20%; it does not reinterpret agreement as a
    // viscous or stall prediction.
    let pass = measured_slope.is_finite()
        && model_slope.is_finite()
        && (0.0..=0.20).contains(&relative_bias)
        && odd_symmetry_worst < 1e-10;
    verdict(
        "bem-006",
        pass,
        &format!(
            "\"detail\":\"NACA 0012 pre-stall lift slope vs NASA TM-4074 table I; inviscid screening only\",\
             \"source\":\"NASA-TM-4074-table-I-M0.15-Re5.97e6-free-transition\",\
             \"measured_dcl_dalpha_per_rad\":{measured_slope:.6},\
             \"panel_dcl_dalpha_per_rad\":{model_slope:.6},\
             \"panel_relative_bias\":{relative_bias:.6},\
             \"admitted_bias\":[0.0,0.20],\"odd_symmetry_worst\":{odd_symmetry_worst:.3e},\
             \"max_abs_alpha_deg\":10.0,\"stall_claim\":false"
        ),
    );

    solve_naca0012_prestall(120, NACA0012_PRESTALL_MAX_ALPHA_RAD)
        .expect("the documented ten-degree boundary is inclusive");
    for alpha in [
        10.000_001f64.to_radians(),
        -10.000_001f64.to_radians(),
        f64::NAN,
    ] {
        let error = solve_naca0012_prestall(120, alpha)
            .expect_err("outside-envelope NACA 0012 request must be refused");
        assert!(matches!(
            error,
            BemError::InvalidScalar {
                name: "NACA 0012 validation angle of attack",
                ..
            }
        ));
    }
}

// ------------------------------------------------------------------ totality

#[test]
fn bem_rejects_invalid_geometry_work_and_trace_requests() {
    assert!(naca4_symmetric(f64::NAN, 80).is_err());
    assert!(naca4_symmetric(0.1, usize::MAX).is_err());
    assert!(panel2d::Airfoil2d::new(vec![[0.0, 0.0]; 4]).is_err());
    let offset = 1.0e12 - 2.0;
    let translated_square = panel2d::Airfoil2d::new(vec![
        [offset, offset],
        [offset, offset + 1.0],
        [offset + 1.0, offset + 1.0],
        [offset + 1.0, offset],
    ]);
    assert!(
        translated_square.is_ok(),
        "orientation and area validation must be translation-stable"
    );

    let foil = naca4_symmetric(0.1, 80).expect("valid fixture");
    assert!(panel2d::solve(&foil, f64::NAN).is_err());
    assert!(WakeSim::new(&foil, 0.1, 0.0, 0.05).is_err());
    assert!(WakeSim::new(&foil, 0.1, 0.05, 5.0e-7).is_err());
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
    assert_eq!(zero.report.rel_residual.to_bits(), 0.0f64.to_bits());
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
