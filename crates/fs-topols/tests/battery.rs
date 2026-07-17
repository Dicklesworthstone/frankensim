//! fs-topols conformance battery (bead 7tv.12).
//!
//! - tls-001 G1: WENO advection order on a smooth rotating field.
//! - tls-002 G0: FIM redistancing — zero-level-set Hausdorff drift,
//!   |∇φ|−1 statistics, idempotency.
//! - tls-003 G0: velocity extension is constant along interface
//!   normals (radial fixture).
//! - tls-004: volume conservation under rigid rotation + the
//!   redistancing-frequency drift policy (audited accumulation).
//! - tls-005 G0: NUMERICAL GATES on the sensitivities — the
//!   topological derivative predicts a real punched hole's compliance
//!   change, and the shape velocity predicts a real uniform boundary
//!   motion's compliance change (adjoint-vs-FD in the bead's sense).
//! - tls-006: the short-cantilever descent — volume convergence,
//!   trajectory stabilization, hole nucleation firing with positive
//!   predicted gain and a genuine topology change, drift audits under
//!   budget, bitwise-deterministic snapshots, and the optimized
//!   topology beating the trivial uniform band at EQUAL volume.

#![allow(clippy::cast_possible_wrap)] // lattice indices are tiny; i64 stencil arithmetic is exact

use fs_cutfem::{
    BoundaryTraction, CutElasticity, CutElasticitySolution, CutFemError, DesignBoxEdge, EdgeBand,
    Quadtree,
};
use fs_material::IsotropicElastic;
use fs_topols::optimize::{Cantilever, material_volume};
use fs_topols::{
    GridSdf, OptimizeSettings, Velocity, advect, extend_velocity, hausdorff, nucleate,
    optimize_compliance, redistance, topological_derivative, zero_crossings,
};
use std::fmt::Write as _;

fn verdict(name: &str, pass: bool, details: &str) {
    println!(
        "{{\"test\":\"{name}\",\"verdict\":\"{}\",{details}}}",
        if pass { "pass" } else { "fail" }
    );
    assert!(pass, "{name} failed: {details}");
}

fn circle_sdf(cx: f64, cy: f64, r: f64) -> impl Fn(f64, f64) -> f64 {
    move |x: f64, y: f64| ((x - cx) * (x - cx) + (y - cy) * (y - cy)).sqrt() - r
}

// ------------------------------------------------------------------ tls-001

#[test]
fn tls_001_weno_advection_order() {
    // Rigid rotation about the center, half a revolution: the exact
    // solution is the initial circle rotated by π.
    let omega = std::f64::consts::TAU;
    let u = move |x: f64, y: f64| [-omega * (y - 0.5), omega * (x - 0.5)];
    let mut errs = Vec::new();
    for n in [24usize, 48, 96] {
        let mut phi = GridSdf::from_fn(n, &circle_sdf(0.68, 0.5, 0.14));
        let band = vec![true; (n + 1) * (n + 1)];
        advect(&mut phi, &band, &Velocity::Linear(&u), 0.5, 0.4);
        // After t = 0.5 (half turn) the circle sits at (0.32, 0.5).
        let exact = circle_sdf(0.32, 0.5, 0.14);
        let mut worst = 0.0f64;
        for j in 0..=n {
            for i in 0..=n {
                let p = phi.pos(i, j);
                let e = exact(p[0], p[1]);
                // Measure near the interface where the SDF is smooth.
                if e.abs() < 3.0 * phi.h() {
                    worst = worst.max((phi.node(i, j) - e).abs());
                }
            }
        }
        errs.push(worst);
    }
    let o1 = (errs[0] / errs[1]).log2();
    let o2 = (errs[1] / errs[2]).log2();
    let pass = o2 > 1.6 && errs[2] < 5e-4;
    verdict(
        "tls-001",
        pass,
        &format!(
            "\"detail\":\"WENO5+RK3 rigid rotation, half revolution\",\
             \"errs\":[{:.3e},{:.3e},{:.3e}],\"orders\":[{o1:.2},{o2:.2}]",
            errs[0], errs[1], errs[2]
        ),
    );
}

// ------------------------------------------------------------------ tls-002

#[test]
fn tls_002_fim_redistancing() {
    let n = 64usize;
    // A level set with badly non-unit gradient but a known zero set.
    let mut phi = GridSdf::from_fn(n, &|x, y| {
        let d = circle_sdf(0.5, 0.5, 0.3)(x, y);
        d * (2.0 + (5.0 * x).sin().abs() + 3.0 * d.abs())
    });
    let truth = GridSdf::from_fn(n, &circle_sdf(0.5, 0.5, 0.3));
    let before = zero_crossings(&phi);
    let want = zero_crossings(&truth);
    let pre_match = hausdorff(&before, &want) / phi.h();
    let audit = redistance(&mut phi, 6.0);
    // Idempotency: a second pass barely moves anything.
    let audit2 = redistance(&mut phi, 6.0);
    let pass = audit.interface_drift_h < 0.2
        && audit.grad_dev_mean < 0.05
        && audit.grad_dev_max < 0.35
        && audit2.interface_drift_h < 0.05
        && pre_match < 0.5;
    verdict(
        "tls-002",
        pass,
        &format!(
            "\"detail\":\"FIM restores |grad phi|=1 without moving the zero set\",\
             \"audit\":{},\"audit_repeat\":{}",
            audit.to_json(),
            audit2.to_json()
        ),
    );
}

// ------------------------------------------------------------------ tls-003

#[test]
fn tls_003_normal_extension_radial() {
    let n = 48usize;
    let phi = GridSdf::from_fn(n, &circle_sdf(0.5, 0.5, 0.25));
    let stride = n + 1;
    let mut vals = vec![0.0f64; stride * stride];
    let mut seeds = vec![false; stride * stride];
    let f = |theta: f64| (2.0 * theta).sin();
    for j in 0..=n {
        for i in 0..=n {
            let k = i + j * stride;
            if phi.node(i, j).abs() <= 1.5 * phi.h() {
                let p = phi.pos(i, j);
                seeds[k] = true;
                vals[k] = f((p[1] - 0.5).atan2(p[0] - 0.5));
            }
        }
    }
    extend_velocity(&phi, &mut vals, &seeds);
    // Along each ray, the extension equals the interface value.
    let mut worst = 0.0f64;
    for j in 0..=n {
        for i in 0..=n {
            let p = phi.pos(i, j);
            let r = (p[0] - 0.5).hypot(p[1] - 0.5);
            if (0.12..=0.38).contains(&r) {
                let want = f((p[1] - 0.5).atan2(p[0] - 0.5));
                worst = worst.max((vals[i + j * stride] - want).abs());
            }
        }
    }
    let pass = worst < 0.15;
    verdict(
        "tls-003",
        pass,
        &format!(
            "\"detail\":\"extension constant along circle normals in a 0.13-wide annulus\",\
             \"worst_dev\":{worst:.3e}"
        ),
    );
}

// ------------------------------------------------------------------ tls-004

#[test]
fn tls_004_volume_conservation_and_drift_policy() {
    let n = 64usize;
    let level = 6u32;
    let grid = Quadtree::uniform(level);
    let omega = std::f64::consts::TAU;
    let u = move |x: f64, y: f64| [-omega * (y - 0.5), omega * (x - 0.5)];
    let mut phi = GridSdf::from_fn(n, &circle_sdf(0.65, 0.5, 0.18));
    let v0 = material_volume(&grid, &phi);
    // Five advect+redistance cycles (a tenth of a revolution each):
    // the drift POLICY ledger — per-cycle audits stay under budget and
    // volume is conserved through the full half revolution.
    let mut drifts = Vec::new();
    for _ in 0..5 {
        let band = vec![true; (n + 1) * (n + 1)];
        advect(&mut phi, &band, &Velocity::Linear(&u), 0.1, 0.4);
        let audit = redistance(&mut phi, 6.0);
        drifts.push(audit.interface_drift_h);
    }
    let v1 = material_volume(&grid, &phi);
    let rel = ((v1 - v0) / v0).abs();
    let drift_max = drifts.iter().copied().fold(0.0f64, f64::max);
    let pass = rel < 0.01 && drift_max < 0.25;
    let mut rows = String::new();
    for d in &drifts {
        let _ = write!(rows, "{d:.3e},");
    }
    verdict(
        "tls-004",
        pass,
        &format!(
            "\"detail\":\"half-revolution volume conservation under 5 redistance cycles\",\
             \"v0\":{v0:.5},\"v1\":{v1:.5},\"rel_drift\":{rel:.3e},\
             \"per_cycle_interface_drift_h\":[{}]",
            rows.trim_end_matches(',')
        ),
    );
}

// ------------------------------------------------------------------ tls-005

fn try_cantilever_solution(
    grid: &Quadtree,
    phi: &GridSdf,
    load: f64,
    band: f64,
) -> Result<CutElasticitySolution, CutFemError> {
    if !(load.is_finite() && load > 0.0) {
        return Err(CutFemError::InvalidElasticityInput {
            what: format!("cantilever load magnitude {load} must be finite and strictly positive"),
        });
    }
    let support = EdgeBand::new(DesignBoxEdge::Right, 0.5 - band, 0.5 + band).map_err(|error| {
        CutFemError::InvalidElasticityInput {
            what: format!(
                "cantilever load-band half-width {band} must be finite and lie in [0, 0.5]: {error}"
            ),
        }
    })?;
    let material = IsotropicElastic::new(1.0, 0.3, 1.0).map_err(|error| {
        CutFemError::InvalidElasticityInput {
            what: format!("cantilever material card was refused: {error}"),
        }
    })?;
    let (lambda, mu) = material.lame();
    let legacy_stiffness_scale = (lambda + 2.0 * mu) / mu;
    let clamp = |x: f64, _y: f64| x < 1e-9;
    let traction = move |_: f64, _: f64| [0.0, -load];
    let solver = CutElasticity {
        grid,
        sdf: phi,
        material: &material,
        nitsche_beta: 20.0 * legacy_stiffness_scale,
        ghost_gamma: 0.5 * legacy_stiffness_scale,
        quad_depth: 2,
        clamp: Some(&clamp),
        boundary_traction: None,
        traction_free_interface: true,
        solver_tol: 1e-12,
        solver_max_iters: 60_000,
    };
    solver.solve_with_boundary_traction(
        &|_, _| [0.0, 0.0],
        &|_, _| [0.0, 0.0],
        BoundaryTraction::EdgeBand {
            support,
            value: &traction,
        },
    )
}

/// Solve the cantilever on a given level set; returns the exact assembled-load
/// compliance `b^T u` for this zero-body-force, zero-interface-data fixture.
fn cantilever_compliance(grid: &Quadtree, phi: &GridSdf, load: f64, band: f64) -> f64 {
    try_cantilever_solution(grid, phi, load, band)
        .expect("cantilever solves")
        .compliance()
}

#[test]
fn typed_cantilever_support_succeeds_only_when_the_loaded_band_is_uncut() {
    let level = 3;
    let n = 1usize << level;
    let grid = Quadtree::uniform(level);
    let wide_beam = GridSdf::from_fn(n, &|_, y| (y - 0.5).abs() - 0.4);
    let first = try_cantilever_solution(&grid, &wide_beam, 1.0, 0.125)
        .expect("grid-aligned support disjoint from both SDF crossings must solve");
    let replay = try_cantilever_solution(&grid, &wide_beam, 1.0, 0.125)
        .expect("typed support solve must replay");
    assert!(first.compliance().is_finite() && first.compliance() > 0.0);
    assert_eq!(first.compliance().to_bits(), replay.compliance().to_bits());

    let crossing = GridSdf::from_fn(n, &|_, y| y - 0.5);
    assert!(matches!(
        try_cantilever_solution(&grid, &crossing, 1.0, 0.125),
        Err(CutFemError::InvalidElasticityInput { .. })
    ));
}

#[test]
fn invalid_cantilever_and_material_settings_fail_before_mutation() {
    let level = 3;
    let n = 1usize << level;
    let original = GridSdf::from_fn(n, &|_, y| (y - 0.5).abs() - 0.4);
    let original_bits = original
        .nodes()
        .iter()
        .map(|value| value.to_bits())
        .collect::<Vec<_>>();
    let settings = OptimizeSettings {
        level,
        iterations: 0,
        ..OptimizeSettings::default()
    };
    for fixture in [
        Cantilever {
            load: f64::NAN,
            band: 0.125,
        },
        Cantilever {
            load: -1.0,
            band: 0.125,
        },
        Cantilever {
            load: 0.0,
            band: 0.125,
        },
        Cantilever {
            load: 1.0,
            band: f64::NAN,
        },
        Cantilever {
            load: 1.0,
            band: -0.125,
        },
        Cantilever {
            load: 1.0,
            band: 0.625,
        },
    ] {
        let mut phi = original.clone();
        assert!(matches!(
            optimize_compliance(&mut phi, fixture, settings),
            Err(CutFemError::InvalidElasticityInput { .. })
        ));
        assert_eq!(
            phi.nodes()
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>(),
            original_bits
        );
    }

    let fixture = Cantilever {
        load: 1.0,
        band: 0.125,
    };
    for invalid_settings in [
        OptimizeSettings {
            youngs: f64::NAN,
            ..settings
        },
        OptimizeSettings {
            youngs: 0.0,
            ..settings
        },
        OptimizeSettings {
            poisson: f64::NAN,
            ..settings
        },
        OptimizeSettings {
            poisson: -1.0,
            ..settings
        },
        OptimizeSettings {
            poisson: 0.49,
            ..settings
        },
    ] {
        let mut phi = original.clone();
        assert!(matches!(
            optimize_compliance(&mut phi, fixture, invalid_settings),
            Err(CutFemError::InvalidElasticityInput { .. })
        ));
        assert_eq!(
            phi.nodes()
                .iter()
                .map(|value| value.to_bits())
                .collect::<Vec<_>>(),
            original_bits
        );
    }
}

#[test]
fn tls_005_sensitivity_numerical_gates() {
    let level = 5u32;
    let n = 1usize << level;
    let grid = Quadtree::uniform(level);
    // A wide beam (interface away from clamp/load).
    let half = 0.42;
    let mut phi = GridSdf::from_fn(n, &move |_, y: f64| (y - 0.5).abs() - half);
    let _ = redistance(&mut phi, 6.0);
    let j0 = cantilever_compliance(&grid, &phi, 1.0, 0.12);
    // Gate A: topological derivative vs a real punched hole.
    let probe = [0.5, 0.5];
    let (lambda, mu) = IsotropicElastic::new(1.0, 0.3, 1.0)
        .expect("valid sensitivity fixture material")
        .lame();
    // Strain at the probe from a fresh solve (reuse the optimizer's
    // energy plumbing indirectly: finite differences of displacement).
    let load = 1.0f64;
    let sol = try_cantilever_solution(&grid, &phi, load, 0.12).expect("typed cantilever solves");
    let eps = strain_probe(&grid, &sol, probe);
    let sxx = (lambda + 2.0 * mu) * eps[0] + lambda * eps[1];
    let syy = lambda * eps[0] + (lambda + 2.0 * mu) * eps[1];
    let sxy = 2.0 * mu * eps[2];
    let dt = topological_derivative(lambda, mu, [sxx, syy, sxy], eps);
    let rho = 3.0 * phi.h();
    let mut holed = phi.clone();
    let holes = nucleate(
        &mut holed,
        &vec![dt; (n + 1) * (n + 1)],
        f64::INFINITY,
        rho,
        0.1,
        1,
    );
    // nucleate picks the best candidate everywhere; force OUR probe by
    // punching manually instead.
    let _ = holes;
    let mut holed = phi.clone();
    for j in 0..=n {
        for i in 0..=n {
            let p = holed.pos(i, j);
            let hole = rho - (p[0] - probe[0]).hypot(p[1] - probe[1]);
            let v = holed.node(i, j);
            *holed.node_mut(i, j) = v.max(hole);
        }
    }
    let _ = redistance(&mut holed, 6.0);
    let j_holed = cantilever_compliance(&grid, &holed, 1.0, 0.12);
    let dj_measured = j_holed - j0;
    let dj_predicted = dt * std::f64::consts::PI * rho * rho;
    let ratio_t = dj_measured / dj_predicted;
    // Gate B: shape velocity — uniform inward motion δ shrinks the
    // beam; predicted dJ = +∫Γ w δ (removing material raises J).
    let delta = 0.75 * phi.h();
    let mut shrunk = GridSdf::from_fn(n, &move |_, y: f64| (y - 0.5).abs() - (half - delta));
    let _ = redistance(&mut shrunk, 6.0);
    let j_shrunk = cantilever_compliance(&grid, &shrunk, 1.0, 0.12);
    // Interface energy density sampled just inside the top/bottom
    // faces, integrated over both (length ≈ 2).
    let mut w_sum = 0.0;
    let mut w_count = 0usize;
    for i in 0..=n {
        #[allow(clippy::cast_precision_loss)]
        let x = i as f64 / n as f64;
        for yy in [0.5 - half + 1.5 * phi.h(), 0.5 + half - 1.5 * phi.h()] {
            let eps = strain_probe(&grid, &sol, [x, yy]);
            let sxx = (lambda + 2.0 * mu) * eps[0] + lambda * eps[1];
            let syy = lambda * eps[0] + (lambda + 2.0 * mu) * eps[1];
            let sxy = 2.0 * mu * eps[2];
            w_sum += 0.5 * (sxx * eps[0] + syy * eps[1] + 2.0 * sxy * eps[2]);
            w_count += 1;
        }
    }
    #[allow(clippy::cast_precision_loss)]
    let w_mean = w_sum / w_count as f64;
    let dj_shape_pred = w_mean * 2.0 * delta; // two faces × length 1
    let ratio_s = (j_shrunk - j0) / dj_shape_pred;
    let pass =
        dj_measured > 0.0 && (0.25..=4.0).contains(&ratio_t) && (0.25..=4.0).contains(&ratio_s);
    verdict(
        "tls-005",
        pass,
        &format!(
            "\"detail\":\"numerical gates: DT vs punched hole; shape velocity vs FD\",\
             \"j0\":{j0:.6e},\"dj_hole_measured\":{dj_measured:.3e},\
             \"dj_hole_predicted\":{dj_predicted:.3e},\"ratio_topological\":{ratio_t:.2},\
             \"dj_shape_measured\":{:.3e},\"dj_shape_predicted\":{dj_shape_pred:.3e},\
             \"ratio_shape\":{ratio_s:.2}",
            j_shrunk - j0
        ),
    );
}

fn strain_probe(grid: &Quadtree, sol: &CutElasticitySolution, p: [f64; 2]) -> [f64; 3] {
    let level = grid.max_level();
    let nf = f64::from(1u32 << level);
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let ci = ((p[0] * nf).floor().clamp(0.0, nf - 1.0)) as u32;
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let cj = ((p[1] * nf).floor().clamp(0.0, nf - 1.0)) as u32;
    let cell = (level, ci, cj);
    let (lo, hi) = grid.rect(cell);
    let corners = grid.corner_nodes(cell);
    let nodal = sol.nodal();
    let mut vals = [[0.0f64; 2]; 4];
    for (a, c) in corners.iter().enumerate() {
        vals[a] = nodal.get(c).copied().unwrap_or([0.0, 0.0]);
    }
    let hx = hi[0] - lo[0];
    let hy = hi[1] - lo[1];
    let xi = ((p[0] - lo[0]) / hx).clamp(0.0, 1.0);
    let et = ((p[1] - lo[1]) / hy).clamp(0.0, 1.0);
    let g = [
        [-(1.0 - et) / hx, -(1.0 - xi) / hy],
        [(1.0 - et) / hx, -xi / hy],
        [et / hx, xi / hy],
        [-et / hx, (1.0 - xi) / hy],
    ];
    let mut gu = [[0.0f64; 2]; 2];
    for a in 0..4 {
        for c in 0..2 {
            gu[c][0] += g[a][0] * vals[a][c];
            gu[c][1] += g[a][1] * vals[a][c];
        }
    }
    [gu[0][0], gu[1][1], f64::midpoint(gu[0][1], gu[1][0])]
}

// ------------------------------------------------------------------ tls-006

#[test]
#[allow(clippy::too_many_lines)]
fn tls_006_cantilever_descent() {
    let settings = OptimizeSettings {
        level: 5,
        volfrac: 0.45,
        iterations: 30,
        nucleation_period: 9,
        mu_al: 1.2,
        ..OptimizeSettings::default()
    };
    let n = 1usize << settings.level;
    let fixture = Cantilever {
        load: 1.0,
        band: 0.12,
    };
    let init = |_: f64, y: f64| (y - 0.5).abs() - 0.42;
    let run = || {
        let mut phi = GridSdf::from_fn(n, &init);
        let _ = redistance(&mut phi, 6.0);
        let report = optimize_compliance(&mut phi, fixture, settings).expect("descent runs");
        (phi, report)
    };
    let (phi, report) = run();
    let (_, report_b) = run();
    let deterministic = report.snapshots == report_b.snapshots;
    let load_pad_deterministic = report.load_pad_nodes == report_b.load_pad_nodes;
    let load_pad_tracked = report.load_pad_nodes.len() == settings.iterations;
    let load_pad_intervened = report.load_pad_nodes.iter().any(|&count| count > 0);
    // Volume convergence.
    let v_final = *report.volume.last().expect("volume history");
    // Trajectory stabilization: compliance variation over the last 3.
    let tail = &report.compliance[report.compliance.len() - 3..];
    let tmax = tail.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let tmin = tail.iter().copied().fold(f64::INFINITY, f64::min);
    let stabilized = (tmax - tmin) / tmax.abs().max(1e-30) < 0.2;
    // Nucleation fired with positive predicted gain; topology changed:
    // an interior void component not touching the box boundary.
    let fired = report.events.iter().any(|e| e.predicted_gain > 0.0);
    let interior_hole = has_interior_void(&phi);
    // Drift audits under budget.
    let drift_ok = report.audits.iter().all(|a| a.interface_drift_h < 0.5);
    // Quality: beats the TRIVIAL design (uniform band) atEQUAL volume.
    let grid = Quadtree::uniform(settings.level);
    let final_support_valid =
        try_cantilever_solution(&grid, &phi, fixture.load, fixture.band).is_ok();
    let half_trivial = 0.5 * v_final;
    let mut trivial = GridSdf::from_fn(n, &move |_, y: f64| (y - 0.5).abs() - half_trivial);
    let _ = redistance(&mut trivial, 6.0);
    let j_trivial = cantilever_compliance(&grid, &trivial, fixture.load, fixture.band);
    let j_final = *report.compliance.last().expect("compliance history");
    let beats_trivial = j_final < j_trivial;
    let pass = (v_final - settings.volfrac).abs() < 0.05
        && stabilized
        && fired
        && interior_hole
        && drift_ok
        && deterministic
        && load_pad_deterministic
        && load_pad_tracked
        && load_pad_intervened
        && final_support_valid
        && beats_trivial;
    let mut rows = String::new();
    for r in report.rows.iter().step_by(4) {
        let _ = write!(rows, "{r},");
    }
    let mut ev = String::new();
    for e in &report.events {
        let _ = write!(ev, "{},", e.to_json());
    }
    verdict(
        "tls-006",
        pass,
        &format!(
            "\"detail\":\"cantilever descent: volume, stabilization, nucleation, \
             determinism, beats-trivial\",\"v_final\":{v_final:.3},\
             \"j_final\":{j_final:.5e},\"j_trivial_same_volume\":{j_trivial:.5e},\
             \"deterministic\":{deterministic},\"interior_hole\":{interior_hole},\
             \"load_pad_deterministic\":{load_pad_deterministic},\
             \"load_pad_tracked\":{load_pad_tracked},\
             \"load_pad_intervened\":{load_pad_intervened},\
             \"final_support_valid\":{final_support_valid},\
             \"events\":[{}],\"trajectory\":[{}]",
            ev.trim_end_matches(','),
            rows.trim_end_matches(',')
        ),
    );
}

/// Flood fill: is there a void (φ > 0) component not touching the box
/// boundary?
fn has_interior_void(phi: &GridSdf) -> bool {
    let n = phi.n();
    let stride = n + 1;
    let void: Vec<bool> = phi.nodes().iter().map(|&v| v > 0.0).collect();
    let mut comp = vec![usize::MAX; stride * stride];
    let mut ncomp = 0usize;
    let mut touches = Vec::new();
    for start in 0..stride * stride {
        if !void[start] || comp[start] != usize::MAX {
            continue;
        }
        let mut stack = vec![start];
        let mut touch = false;
        comp[start] = ncomp;
        while let Some(k) = stack.pop() {
            let (i, j) = (k % stride, k / stride);
            if i == 0 || j == 0 || i == n || j == n {
                touch = true;
            }
            for (di, dj) in [(1i64, 0i64), (-1, 0), (0, 1), (0, -1)] {
                let (ni, nj) = (i as i64 + di, j as i64 + dj);
                if ni < 0 || nj < 0 || ni > n as i64 || nj > n as i64 {
                    continue;
                }
                #[allow(clippy::cast_sign_loss)]
                let nk = ni as usize + nj as usize * stride;
                if void[nk] && comp[nk] == usize::MAX {
                    comp[nk] = ncomp;
                    stack.push(nk);
                }
            }
        }
        touches.push(touch);
        ncomp += 1;
    }
    touches.iter().any(|t| !t)
}
