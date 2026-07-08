//! fs-solid stability battery (bead tfz.15).
//!
//! - stab-001 G0+G2: pencil symmetry, mode K-orthogonality, Euler
//!   column critical load vs the analytic reference, Richardson
//!   discretization indicator on λ.
//! - stab-002: pseudo-arclength continuation VALIDATED AGAINST CLOSED
//!   FORM — the shallow von Mises truss: limit points at the analytic
//!   locations, snap-through traversed, load control demonstrably
//!   failing where arclength does not; checkpoint/resume determinism.
//! - stab-003: continuum snap-through — a shallow hyperelastic arch
//!   traced through both limit points on one path (the fixture load
//!   control cannot finish).
//! - stab-004: branch-point detection on the compressed column at the
//!   pencil's predicted load; branch switching lands on the bent
//!   post-buckling branch.
//! - stab-005: eigenvalue derivative gradient gate (direct pencil
//!   derivative vs FD at frozen prebuckling stress) and the documented
//!   clustered-eigenvalue trap: min() kinks across a crossing while
//!   the KS aggregate stays smooth and its derivative matches FD.

use fs_solid::continuation::{ArcSettings, PathEvent, PathResidual, PathState, advance};
use fs_solid::linear::{Formulation, LinearProblem, PlaneKind};
use fs_solid::{
    HyperProblem, Mesh2, NewtonSettings, Patch, buckling_loads, eigenvalue_derivative,
    expand_mode, group_stiffness, ks_aggregate, ks_aggregate_derivative, lambda_indicator,
    reduced_pencil, switch_branch,
};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_solid::SolidError;
use fs_sparse::{Coo, Csr};
use std::fmt::Write as _;

fn verdict(name: &str, pass: bool, details: &str) {
    println!(
        "{{\"test\":\"{name}\",\"verdict\":\"{}\",{details}}}",
        if pass { "pass" } else { "fail" }
    );
    assert!(pass, "{name} failed: {details}");
}

// --------------------------------------------------------------- fixtures

const COL_L: f64 = 8.0;
const COL_T: f64 = 0.4;
const COL_P: f64 = 1e-4; // reference edge traction (per length)

fn column_problem(mesh: &Mesh2) -> LinearProblem<'_> {
    LinearProblem {
        mesh,
        youngs: 1.0,
        poisson: 0.0, // beam-like: no Poisson stiffening of the Euler value
        plane: PlaneKind::Stress,
        formulation: Formulation::Standard,
        body_force: None,
        dirichlet: vec![(Patch::Left, &|_, _| [0.0, 0.0])],
        traction: vec![(Patch::Right, &|_, _| [-COL_P, 0.0])],
    }
}

fn euler_load() -> f64 {
    // Clamped-free strut: Pcr = π²EI/(4L²), I = t³/12, E = 1.
    let i = COL_T * COL_T * COL_T / 12.0;
    std::f64::consts::PI.powi(2) * i / (4.0 * COL_L * COL_L)
}

// ------------------------------------------------------------------ stab-001

#[test]
fn stab_001_euler_column_pencil() {
    let coarse_mesh = Mesh2::quads(COL_L, COL_T, 20, 2);
    let fine_mesh = Mesh2::quads(COL_L, COL_T, 40, 4);
    let mut lambdas = Vec::new();
    let mut sym_err = 0.0f64;
    let mut ortho = 0.0f64;
    for mesh in [&coarse_mesh, &fine_mesh] {
        let problem = column_problem(mesh);
        let (k, kg, dof_map, _u0) = reduced_pencil(&problem).expect("pencil builds");
        // G0: K_G symmetry.
        for r in 0..kg.nrows() {
            let (cols, vals) = kg.row(r);
            for (c, v) in cols.iter().zip(vals) {
                sym_err = sym_err.max((v - kg.get(*c, r)).abs());
            }
        }
        let res = buckling_loads(&k, &kg, &dof_map, 2, 400).expect("pencil solves");
        // G0: K-orthogonality of distinct modes.
        if res.modes.len() >= 2 {
            let n = res.modes[0].len();
            let mut t = vec![0.0f64; n];
            k.spmv(&res.modes[0], &mut t);
            let x12: f64 = res.modes[1].iter().zip(&t).map(|(a, b)| a * b).sum();
            let x11: f64 = res.modes[0].iter().zip(&t).map(|(a, b)| a * b).sum();
            ortho = ortho.max((x12 / x11).abs());
        }
        lambdas.push(res.loads[0]);
    }
    // Critical LOAD (traction × edge height) vs Euler.
    let p_ref = COL_P * COL_T;
    let (extrap, indicator) = lambda_indicator(lambdas[0], lambdas[1]);
    let fem_pcr = lambdas[1] * p_ref;
    let rel_fine = (fem_pcr - euler_load()).abs() / euler_load();
    // Q1 parasitic shear inflates the raw fine value; the Richardson
    // EXTRAPOLATION lands on Euler and the indicator honestly covers
    // the fine-value gap — which is precisely the bead's "λ with
    // discretization-error indicator" deliverable.
    let rel = (extrap * p_ref - euler_load()).abs() / euler_load();
    let covered = indicator >= 0.9 * (extrap - lambdas[1]).abs();
    let pass = sym_err < 1e-10 && ortho < 1e-6 && rel < 0.03 && covered;
    verdict(
        "stab-001",
        pass,
        &format!(
            "\"detail\":\"Euler clamped-free strut, LOBPCG pencil\",\
             \"fem_pcr\":{fem_pcr:.6e},\"euler\":{:.6e},\"rel_extrap\":{rel:.4},\"rel_fine\":{rel_fine:.4},\
             \"kg_sym\":{sym_err:.2e},\"mode_k_ortho\":{ortho:.2e},\
             \"lambda_extrap\":{extrap:.4},\"lambda_indicator\":{indicator:.3e}",
            euler_load()
        ),
    );
}

// ------------------------------------------------------------------ stab-002

/// The shallow von Mises truss: one DOF (crown drop w), closed-form
/// equilibrium — the continuation algorithm's exact oracle.
struct VonMises {
    b: f64,
    h: f64,
    ea: f64,
}

impl VonMises {
    fn internal(&self, w: f64) -> f64 {
        let l0 = (self.b * self.b + self.h * self.h).sqrt();
        let z = self.h - w;
        let l = (self.b * self.b + z * z).sqrt();
        let strain = (l - l0) / l0;
        // dΠ_int/dw with Π_int = 2·(EA/2L0)(L−L0)².
        -2.0 * self.ea * strain * z / l
    }

    /// λ(w) on the equilibrium manifold (load = internal force).
    fn lambda_of(&self, w: f64) -> f64 {
        self.internal(w)
    }

    /// The analytic limit points: extrema of λ(w) on (0, 2h).
    fn limit_ws(&self) -> (f64, f64) {
        // Dense scan + refinement (test-side oracle; exactness by
        // bisection on dλ/dw sign changes).
        let dw = 1e-6;
        let dldw = |w: f64| (self.lambda_of(w + dw) - self.lambda_of(w - dw)) / (2.0 * dw);
        let mut roots = Vec::new();
        let mut prev = dldw(0.0);
        let mut wprev = 0.0;
        let steps = 4000;
        for i in 1..=steps {
            #[allow(clippy::cast_precision_loss)]
            let w = 2.0 * self.h * (i as f64) / (steps as f64);
            let d = dldw(w);
            if prev * d < 0.0 {
                let (mut a, mut c) = (wprev, w);
                for _ in 0..60 {
                    let m = 0.5 * (a + c);
                    if dldw(a) * dldw(m) <= 0.0 {
                        c = m;
                    } else {
                        a = m;
                    }
                }
                roots.push(0.5 * (a + c));
            }
            prev = d;
            wprev = w;
        }
        (roots[0], roots[1])
    }
}

impl PathResidual for VonMises {
    fn ndof(&self) -> usize {
        1
    }
    fn residual(&self, u: &[f64], lambda: f64) -> Result<Vec<f64>, SolidError> {
        Ok(vec![self.internal(u[0]) - lambda])
    }
    fn tangent(&self, u: &[f64], _lambda: f64) -> Result<Csr, SolidError> {
        let dw = 1e-7;
        let d = (self.internal(u[0] + dw) - self.internal(u[0] - dw)) / (2.0 * dw);
        let mut coo = Coo::new(1, 1);
        coo.push(0, 0, d);
        Ok(coo.assemble())
    }
    fn load_vector(&self) -> Vec<f64> {
        vec![1.0]
    }
}

#[test]
fn stab_002_von_mises_truss_closed_form() {
    let truss = VonMises {
        b: 1.0,
        h: 0.2,
        ea: 1.0,
    };
    let settings = ArcSettings {
        ds: 0.005,
        ds_max: 0.02,
        ..ArcSettings::default()
    };
    let mut path = PathState::start(1, &settings);
    advance(&truss, &mut path, &settings, 400).expect("truss path advances");
    let (w1, w2) = truss.limit_ws();
    let (l1, l2) = (truss.lambda_of(w1), truss.lambda_of(w2));
    // Limit-point events near the analytic locations.
    let limit_lambdas: Vec<f64> = path
        .events
        .iter()
        .filter_map(|e| match e {
            PathEvent::LimitPoint { lambda, .. } => Some(*lambda),
            PathEvent::BranchPoint { .. } => None,
        })
        .collect();
    let hit1 = limit_lambdas
        .iter()
        .any(|l| (l - l1).abs() < 0.05 * l1.abs().max(1e-12));
    let hit2 = limit_lambdas
        .iter()
        .any(|l| (l - l2).abs() < 0.10 * l2.abs().max(1e-12));
    // Snap-through traversed: the path reaches beyond w = 2h.
    let w_end = path.u[0];
    // Every accepted point sits on the analytic manifold.
    let manifold_dev = path
        .trace
        .iter()
        .zip(0..)
        .map(|((l, w), _)| (truss.lambda_of(*w) - l).abs())
        .fold(0.0f64, f64::max);
    // Load control just above the first limit load CANNOT trace the
    // near branch — it either diverges or JUMPS discontinuously to the
    // far branch (w beyond the second limit point): the snap that
    // arclength traces continuously.
    let target = 1.05 * l1;
    let mut w = 0.0f64;
    let mut newton_failed = true;
    for _ in 0..60 {
        let r = truss.internal(w) - target;
        if r.abs() < 1e-12 {
            newton_failed = w > w2; // converged = jumped past the snap
            break;
        }
        let dw = 1e-7;
        let d = (truss.internal(w + dw) - truss.internal(w - dw)) / (2.0 * dw);
        let step = r / d;
        if !step.is_finite() || step.abs() > 10.0 {
            break;
        }
        w -= step;
    }
    // Checkpoint/resume determinism: 110+110 equals 220 bitwise.
    let mut a = PathState::start(1, &settings);
    advance(&truss, &mut a, &settings, 200).expect("first half");
    let mut b = a.clone();
    advance(&truss, &mut b, &settings, 200).expect("second half");
    let resumed_identical = b.u[0].to_bits() == path.u[0].to_bits()
        && b.lambda.to_bits() == path.lambda.to_bits();
    let pass = hit1 && hit2 && w_end > 2.0 * truss.h && manifold_dev < 1e-7 && newton_failed
        && resumed_identical;
    verdict(
        "stab-002",
        pass,
        &format!(
            "\"detail\":\"closed-form von Mises truss oracle\",\
             \"analytic_limits\":[{l1:.6e},{l2:.6e}],\
             \"detected_limit_lambdas\":{limit_lambdas:?},\
             \"w_end\":{w_end:.4},\"manifold_dev\":{manifold_dev:.2e},\
             \"load_control_fails\":{newton_failed},\
             \"resume_bitwise\":{resumed_identical}"
        ),
    );
}

// ------------------------------------------------------------------ stab-003

#[test]
fn stab_003_continuum_arch_snap_through() {
    // A shallow clamped-clamped hyperelastic arch under a crown load.
    // A continuum von Mises TENT: two inclined strips meeting at a
    // crown. The snap mechanism is AXIAL shortening (locking-proof) —
    // a curved bending-dominated arch under full-integration Q1 goes
    // membrane-stiff and refuses to snap (the sol-005 pathology).
    let half_span = 1.0;
    let rise = 0.3;
    let thick = 0.08;
    let mesh = Mesh2::mapped_quads(32, 2, &|s, t| {
        [
            2.0 * half_span * (s - 0.5),
            rise * (1.0 - (2.0 * s - 1.0).abs()) + thick * (t - 0.5),
        ]
    });
    let card = Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 1.0,
            lambda: 3.0,
        },
        4.0,
    )
    .expect("card");
    let problem = HyperProblem {
        mesh: &mesh,
        material: &card,
        dirichlet: vec![
            (Patch::Left, &|_, _| [0.0, 0.0]),
            (Patch::Right, &|_, _| [0.0, 0.0]),
        ],
        traction: vec![(Patch::Top, &|x: f64, _y: f64| {
            // Crown-concentrated downward load band.
            if x.abs() < 0.15 { [0.0, -0.002] } else { [0.0, 0.0] }
        })],
        settings: NewtonSettings::default(),
    };
    // SELF-CALIBRATING load scale: a fixed-load Newton probe marches
    // up until it fails — locating the snap region AND demonstrating
    // load-control failure on the continuum fixture (the thing
    // arclength exists to fix).
    let fixed_newton = |u0: &[f64], lam: f64| -> Option<Vec<f64>> {
        let mut u = u0.to_vec();
        for _ in 0..25 {
            let r = problem.residual(&u, lam).ok()?;
            let rn: f64 = r.iter().map(|x| x * x).sum::<f64>().sqrt();
            if rn < 1e-9 {
                return Some(u);
            }
            let k = problem.tangent(&u, lam).ok()?;
            let f = fs_la::factor::lu(&k.to_dense(), k.nrows()).ok()?;
            let mut d: Vec<f64> = r.iter().map(|x| -x).collect();
            f.solve(&mut d);
            if !d.iter().all(|x| x.is_finite()) {
                return None;
            }
            for (ui, di) in u.iter_mut().zip(&d) {
                *ui += di;
            }
        }
        None
    };
    let mut u_warm = vec![0.0f64; PathResidual::ndof(&problem)];
    let mut lambda_fail = None;
    let mut lam = 0.0f64;
    for _ in 0..60 {
        lam += 1.0;
        match fixed_newton(&u_warm, lam) {
            Some(u) => u_warm = u,
            None => {
                lambda_fail = Some(lam);
                break;
            }
        }
    }
    let Some(lam_fail) = lambda_fail else {
        verdict(
            "stab-003",
            false,
            "\"detail\":\"fixture never lost load control - no snap present\"",
        );
        return;
    };
    let settings = ArcSettings {
        ds: lam_fail / 60.0,
        ds_max: lam_fail / 15.0,
        ..ArcSettings::default()
    };
    let mut path = PathState::start(PathResidual::ndof(&problem), &settings);
    eprintln!("stab-003: load control fails at lambda = {lam_fail}");
    advance(&problem, &mut path, &settings, 100).expect("arch path advances (1)");
    eprintln!("stab-003: 100 steps, lambda = {}", path.lambda);
    advance(&problem, &mut path, &settings, 100).expect("arch path advances (2)");
    eprintln!("stab-003: 200 steps, lambda = {}", path.lambda);
    let limits: Vec<f64> = path
        .events
        .iter()
        .filter_map(|e| match e {
            PathEvent::LimitPoint { lambda, .. } => Some(*lambda),
            PathEvent::BranchPoint { .. } => None,
        })
        .collect();
    // Trace-level snap evidence: the load direction reverses at least
    // twice along one continuous path.
    let mut flips = 0;
    let mut prev_dir = 0.0f64;
    for pair in path.trace.windows(2) {
        let d = pair[1].0 - pair[0].0;
        if d * prev_dir < 0.0 {
            flips += 1;
        }
        if d != 0.0 {
            prev_dir = d;
        }
    }
    let max_defl = path.trace.iter().map(|&(_, d)| d).fold(0.0f64, f64::max);
        let lam_end = path.lambda;
    let pass = !limits.is_empty() && flips >= 2 && max_defl > 0.8 * rise && lam_end > lam_fail;
    verdict(
        "stab-003",
        pass,
        &format!(
            "\"detail\":\"shallow arch: load control fails, arclength traces through\",\
             \"load_control_fails_at\":{lam_fail:.3},\"limit_lambdas\":{limits:?},\
             \"trace_direction_flips\":{flips},\"max_defl\":{max_defl:.4},\
             \"rise\":{rise:.4},\"lambda_end\":{lam_end:.3},\"steps\":{}",
            path.step
        ),
    );
}

// ------------------------------------------------------------------ stab-004

#[test]
fn stab_004_branch_point_and_switching() {
    // The compressed column, now hyperelastic, followed along the
    // fundamental branch past the pencil's predicted bifurcation.
    let mesh = Mesh2::quads(COL_L, COL_T, 30, 3);
    let linear = column_problem(&mesh);
    let (k, kg, dof_map, _) = reduced_pencil(&linear).expect("pencil builds");
    let pencil = buckling_loads(&k, &kg, &dof_map, 1, 400).expect("pencil solves");
    let lambda_cr = pencil.loads[0];
    let mode = expand_mode(&pencil.modes[0], &pencil.dof_map);
    let card = Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 0.5,
            lambda: 0.0,
        },
        4.0,
    )
    .expect("card"); // E = 3μ·(2·0.5/…) — matches E=1, ν=0 small-strain
    let problem = HyperProblem {
        mesh: &mesh,
        material: &card,
        dirichlet: vec![(Patch::Left, &|_, _| [0.0, 0.0])],
        traction: vec![(Patch::Right, &|_, _| [-COL_P, 0.0])],
        settings: NewtonSettings::default(),
    };
    let settings = ArcSettings {
        ds: 0.4,
        ds_max: 1.2,
        ..ArcSettings::default()
    };
    let mut path = PathState::start(PathResidual::ndof(&problem), &settings);
    advance(&problem, &mut path, &settings, 60).expect("fundamental branch advances");
    let branch_at = path.events.iter().find_map(|e| match e {
        PathEvent::BranchPoint { lambda, .. } => Some(*lambda),
        PathEvent::LimitPoint { .. } => None,
    });
    let pre_transverse = path
        .u
        .iter()
        .skip(1)
        .step_by(2)
        .fold(0.0f64, |m, &v| m.max(v.abs()));
    // Switch onto the bent branch and keep following.
    let detected = branch_at.unwrap_or(f64::NAN);
    let mut pass = false;
    let mut post_transverse = 0.0f64;
    if let Some(lam_bp) = branch_at {
        switch_branch(&mut path, &mode, 0.08);
        if advance(&problem, &mut path, &settings, 25).is_ok() {
            post_transverse = path
                .u
                .iter()
                .skip(1)
                .step_by(2)
                .fold(0.0f64, |m, &v| m.max(v.abs()));
            let rel = (lam_bp - lambda_cr).abs() / lambda_cr;
            pass = rel < 0.15 && pre_transverse < 1e-4 && post_transverse > 20.0 * pre_transverse.max(1e-9);
        }
    }
    verdict(
        "stab-004",
        pass,
        &format!(
            "\"detail\":\"branch detection at pencil load + switch onto bent branch\",\
             \"pencil_lambda_cr\":{lambda_cr:.4},\"detected_at\":{detected:.4},\
             \"pre_transverse\":{pre_transverse:.3e},\
             \"post_transverse\":{post_transverse:.3e}"
        ),
    );
}

// ------------------------------------------------------------------ stab-005

#[test]
fn stab_005_eigenvalue_derivatives_and_cluster_trap() {
    let mesh = Mesh2::quads(COL_L, COL_T, 30, 3);
    let problem = column_problem(&mesh);
    let (k, kg, dof_map, _) = reduced_pencil(&problem).expect("pencil builds");
    let res = buckling_loads(&k, &kg, &dof_map, 1, 400).expect("pencil solves");
    // Parameter: Young's scale on the clamped half of the strut.
    let group = |e: usize| (e % 30) < 15; // left half columns of each row
    let dk = group_stiffness(&problem, &dof_map, &group);
    let dlam = eigenvalue_derivative(&res.modes[0], &dk, &kg);
    // FD at frozen K_G: λ(s) from the pencil (K + (s−1)·K_group, K_G).
    let fd = {
        let mut lams = Vec::new();
        for s in [0.999f64, 1.001] {
            let mut coo = Coo::new(k.nrows(), k.nrows());
            for r in 0..k.nrows() {
                let (cols, vals) = k.row(r);
                for (c, v) in cols.iter().zip(vals) {
                    coo.push(r, *c, *v);
                }
                let (cols, vals) = dk.row(r);
                for (c, v) in cols.iter().zip(vals) {
                    coo.push(r, *c, (s - 1.0) * v);
                }
            }
            let ks = coo.assemble();
            let r = buckling_loads(&ks, &kg, &dof_map, 1, 400).expect("perturbed pencil");
            lams.push(r.loads[0]);
        }
        (lams[1] - lams[0]) / 0.002
    };
    let grad_rel = (dlam - fd).abs() / fd.abs().max(1e-30);
    // The clustered trap, demonstrated where it lives: two branches
    // crossing under a parameter. min() has a kink; the KS aggregate
    // is smooth and its derivative matches FD through the crossing.
    let branches = |s: f64| [1.0 + s, 2.0 - s];
    let dbranches = [1.0, -1.0];
    let rho = 200.0;
    let s0 = 0.5; // the crossing
    let eps = 1e-5;
    let ks_fd = (ks_aggregate(&branches(s0 + eps), rho)
        - ks_aggregate(&branches(s0 - eps), rho))
        / (2.0 * eps);
    let ks_an = ks_aggregate_derivative(&branches(s0), &dbranches, rho);
    let ks_rel = (ks_fd - ks_an).abs() / ks_an.abs().max(1e-30);
    let min_fd_left = (branches(s0).iter().copied().fold(f64::INFINITY, f64::min)
        - branches(s0 - eps).iter().copied().fold(f64::INFINITY, f64::min))
        / eps;
    let min_fd_right = (branches(s0 + eps).iter().copied().fold(f64::INFINITY, f64::min)
        - branches(s0).iter().copied().fold(f64::INFINITY, f64::min))
        / eps;
    let kink = (min_fd_left - min_fd_right).abs(); // ≈ 2: the trap
    let conservative = ks_aggregate(&branches(s0), rho)
        <= branches(s0).iter().copied().fold(f64::INFINITY, f64::min);
    let pass = grad_rel < 1e-3 && ks_rel < 1e-6 && kink > 1.5 && conservative;
    let mut details = String::new();
    let _ = write!(
        details,
        "\"detail\":\"direct pencil derivative vs FD; KS aggregate through a crossing\",\
         \"dlambda_direct\":{dlam:.6e},\"dlambda_fd\":{fd:.6e},\"grad_rel\":{grad_rel:.3e},\
         \"ks_deriv_rel\":{ks_rel:.3e},\"min_kink\":{kink:.3},\"ks_conservative\":{conservative}"
    );
    verdict("stab-005", pass, &details);
}
