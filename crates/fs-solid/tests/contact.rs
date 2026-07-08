//! fs-solid contact battery (bead tfz.16): barrier calculus gates,
//! the intersection-free torture audit, resting force balance,
//! tight-clearance squeeze, the friction cone, and the differentiable
//! contact gradient — all against exact SDF charts.

use fs_cutfem::sdf::{Circle, CutSdf, HalfPlane};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_solid::contact::{Barrier, ContactProblem, Friction};
use fs_solid::hyper2d::{HyperProblem, NewtonSettings};
use fs_solid::mesh2::{Mesh2, Patch};

fn verdict(name: &str, pass: bool, details: &str) {
    println!("{{\"test\":\"{name}\",\"pass\":{pass},\"details\":\"{details}\"}}");
    assert!(pass, "{name}: {details}");
}

fn card() -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 1.0,
            lambda: 3.0,
        },
        3.0,
    )
    .expect("card")
}

/// Bottom-middle x-pins so a frictionless floor cannot let the block
/// float sideways (the elastic probe is unconstrained).
fn x_pins(mesh: &Mesh2) -> Vec<(usize, f64)> {
    let mut pins = Vec::new();
    for (i, p) in mesh.nodes.iter().enumerate() {
        if p[1] <= 1e-12 && (p[0] - 0.5).abs() < 0.26 {
            pins.push((2 * i, 0.0));
        }
    }
    pins
}

/// cnt-001: barrier calculus — zero at and beyond d̂, C¹ at the
/// activation boundary, monotone decreasing on (0, d̂), divergent as
/// d → 0⁺, and derivatives consistent with finite differences.
#[test]
fn cnt_001_barrier_calculus() {
    let b = Barrier {
        kappa: 2.0,
        dhat: 0.1,
    };
    let smooth_at_dhat = b.value(0.1).abs() < 1e-15
        && b.d1(0.0999).abs() < 1e-3
        && b.value(0.15).abs() < 1e-30;
    let mut monotone = true;
    let mut prev = f64::INFINITY;
    for k in 1..100 {
        let d = 0.1 * f64::from(k) / 100.0;
        let v = b.value(d);
        if v > prev {
            monotone = false;
        }
        prev = v;
    }
    let diverges = b.value(1e-9) > 1e6;
    let mut worst_fd = 0.0f64;
    for &d in &[0.02f64, 0.05, 0.08] {
        let h = 1e-7;
        let fd1 = (b.value(d + h) - b.value(d - h)) / (2.0 * h);
        let fd2 = (b.d1(d + h) - b.d1(d - h)) / (2.0 * h);
        worst_fd = worst_fd
            .max((fd1 - b.d1(d)).abs() / b.d1(d).abs().max(1.0))
            .max((fd2 - b.d2(d)).abs() / b.d2(d).abs().max(1.0));
    }
    verdict(
        "cnt-001-barrier-calculus",
        smooth_at_dhat && monotone && diverges && worst_fd < 1e-5,
        &format!(
            "C1 at dhat, monotone, b(1e-9)={:.1e}, worst FD dev {worst_fd:.1e}",
            b.value(1e-9)
        ),
    );
}

/// cnt-002/003: a block pressed HARD onto a frictionless floor —
/// every accepted Newton iterate keeps positive gap (the audit trail
/// IS the IPC guarantee), and at rest the summed contact reactions
/// balance the applied load exactly.
#[test]
fn cnt_002_003_torture_and_balance() {
    let mesh = Mesh2::quads(1.0, 1.0, 6, 6);
    let material = card();
    let q = 0.25; // strong downward traction per unit length
    let press = move |_: f64, _: f64| [0.0, -q];
    let hyper = HyperProblem {
        mesh: &mesh,
        material: &material,
        dirichlet: vec![],
        traction: vec![(Patch::Top, &press)],
        settings: NewtonSettings::default(),
    };
    let floor = HalfPlane {
        normal: [0.0, 1.0],
        offset: -0.02, // floor at y = −0.02: initial gap 0.02
    };
    let problem = ContactProblem {
        hyper: &hyper,
        sdf: &floor,
        barrier: Barrier {
            kappa: 10.0,
            dhat: 0.05,
        },
        friction: None,
        pins: x_pins(&mesh),
        max_newton: 200,
        tol: 1e-9,
    };
    let sol = problem.solve(1.0).expect("contact equilibrium");
    verdict(
        "cnt-002-intersection-free",
        sol.min_gap_ever > 0.0 && sol.min_gap > 0.0 && sol.residual < 1e-9,
        &format!(
            "min gap ever {:.3e}, at solution {:.3e}, residual {:.1e}, {} Newton iters",
            sol.min_gap_ever, sol.min_gap, sol.residual, sol.iterations
        ),
    );
    let total_reaction: f64 = sol.reactions.iter().map(|r| r[1]).sum();
    let applied = q * 1.0; // traction × top-edge length
    verdict(
        "cnt-003-force-balance",
        (total_reaction - applied).abs() < 1e-6 * applied,
        &format!(
            "sum of contact reactions {total_reaction:.8} vs applied load {applied:.8} (barrier energy {:.3e})",
            sol.barrier_energy
        ),
    );
}

/// cnt-004: tight clearance — the block squeezed between floor and
/// ceiling keeps positive gap on BOTH sides, and a curved obstacle
/// (circle) exercises the non-planar gradient path with the same
/// audit.
#[test]
fn cnt_004_tight_clearance_and_circle() {
    let mesh = Mesh2::quads(1.0, 1.0, 6, 6);
    let material = card();
    let hyper = HyperProblem {
        mesh: &mesh,
        material: &material,
        dirichlet: vec![],
        traction: vec![(Patch::Top, &|_, _| [0.0, -0.15])],
        settings: NewtonSettings::default(),
    };
    /// min() of two charts: exact where the closer obstacle wins; the
    /// kink sits mid-channel, far outside barrier support.
    struct Channel {
        floor: HalfPlane,
        ceiling: HalfPlane,
    }
    impl CutSdf for Channel {
        fn value(&self, p: [f64; 2]) -> f64 {
            self.floor.value(p).min(self.ceiling.value(p))
        }
        fn gradient(&self, p: [f64; 2]) -> [f64; 2] {
            if self.floor.value(p) <= self.ceiling.value(p) {
                self.floor.gradient(p)
            } else {
                self.ceiling.gradient(p)
            }
        }
        fn enclose(&self, lo: [f64; 2], hi: [f64; 2]) -> fs_ivl::Interval {
            let a = self.floor.enclose(lo, hi);
            let b = self.ceiling.enclose(lo, hi);
            fs_ivl::Interval::new(a.lo().min(b.lo()), a.hi().min(b.hi()))
        }
    }
    let channel = Channel {
        floor: HalfPlane {
            normal: [0.0, 1.0],
            offset: -0.02,
        },
        // Ceiling at y = 1.02: φ = 1.02 − y.
        ceiling: HalfPlane {
            normal: [0.0, -1.0],
            offset: -1.02,
        },
    };
    let problem = ContactProblem {
        hyper: &hyper,
        sdf: &channel,
        barrier: Barrier {
            kappa: 10.0,
            dhat: 0.015,
        },
        friction: None,
        pins: x_pins(&mesh),
        max_newton: 200,
        tol: 1e-9,
    };
    let sol = problem.solve(1.0).expect("squeeze equilibrium");
    verdict(
        "cnt-004-tight-clearance",
        sol.min_gap_ever > 0.0 && sol.min_gap > 0.0,
        &format!(
            "channel squeeze: min gap ever {:.3e} (audited every iterate)",
            sol.min_gap_ever
        ),
    );
    // Curved obstacle: bump under the block's midline.
    let bump = Circle {
        center: [0.5, -0.45],
        radius: 0.44, // apex at y = −0.01: gap 0.01 under the center
    };
    let problem2 = ContactProblem {
        hyper: &hyper,
        sdf: &bump,
        barrier: Barrier {
            kappa: 10.0,
            dhat: 0.02,
        },
        friction: None,
        pins: x_pins(&mesh),
        max_newton: 200,
        tol: 1e-9,
    };
    let sol2 = problem2.solve(1.0).expect("circle equilibrium");
    verdict(
        "cnt-004-curved-obstacle",
        sol2.min_gap_ever > 0.0 && sol2.residual < 1e-9,
        &format!(
            "circle bump: min gap ever {:.3e}, residual {:.1e}",
            sol2.min_gap_ever, sol2.residual
        ),
    );
}

/// cnt-005: the friction cone — the same tangentially-loaded block
/// slips an order of magnitude more when the tangential load exceeds
/// μ times the normal load than when it sits inside the cone.
#[test]
fn cnt_005_friction_cone() {
    let mesh = Mesh2::quads(1.0, 1.0, 6, 6);
    let material = card();
    let mu = 0.5;
    let q = 0.2;
    let run = |tx: f64| -> f64 {
        let pull = move |_: f64, _: f64| [tx, -q];
        let hyper = HyperProblem {
            mesh: &mesh,
            material: &material,
            dirichlet: vec![],
            traction: vec![(Patch::Top, &pull)],
            settings: NewtonSettings::default(),
        };
        let floor = HalfPlane {
            normal: [0.0, 1.0],
            offset: -0.02,
        };
        let problem = ContactProblem {
            hyper: &hyper,
            sdf: &floor,
            barrier: Barrier {
                kappa: 10.0,
                dhat: 0.05,
            },
            friction: Some(Friction {
                mu,
                eps_v: 1e-4,
                rounds: 4,
            }),
            pins: Vec::new(), // friction grounds the tangential mode
            max_newton: 300,
            tol: 1e-8,
        };
        let sol = problem.solve(1.0).expect("friction equilibrium");
        assert!(sol.min_gap_ever > 0.0, "intersection-free under friction");
        // Mean tangential displacement of the bottom row.
        let mut acc = 0.0;
        let mut cnt = 0.0;
        for (i, p) in mesh.nodes.iter().enumerate() {
            if p[1] <= 1e-12 {
                acc += sol.u[2 * i];
                cnt += 1.0;
            }
        }
        (acc / cnt).abs()
    };
    let stick = run(0.5 * mu * q); // half the cone
    let slip = run(2.0 * mu * q); // double the cone
    verdict(
        "cnt-005-friction-cone",
        slip > 10.0 * stick.max(1e-12) && stick < 5e-3,
        &format!("bottom-row slip: inside cone {stick:.3e}, outside cone {slip:.3e} (ratio {:.1})", slip / stick.max(1e-12)),
    );
}

/// cnt-006: differentiable contact — dJ/dh for a rigid floor
/// translation matches central finite differences (the plane makes
/// the dropped-curvature term exact).
#[test]
fn cnt_006_adjoint_gradient() {
    let mesh = Mesh2::quads(1.0, 1.0, 4, 4);
    let material = card();
    let hyper = HyperProblem {
        mesh: &mesh,
        material: &material,
        dirichlet: vec![],
        traction: vec![(Patch::Top, &|_, _| [0.0, -0.2])],
        settings: NewtonSettings::default(),
    };
    let solve_at = |h: f64| -> (fs_solid::contact::ContactSolution, ContactProblem<'_>) {
        let floor = HalfPlane {
            normal: [0.0, 1.0],
            offset: -0.02 + h,
        };
        let floor: &'static HalfPlane = Box::leak(Box::new(floor));
        let problem = ContactProblem {
            hyper: &hyper,
            sdf: floor,
            barrier: Barrier {
                kappa: 10.0,
                dhat: 0.05,
            },
            friction: None,
            pins: x_pins(&mesh),
            max_newton: 200,
            tol: 1e-11,
        };
        let sol = problem.solve(1.0).expect("equilibrium");
        (sol, problem)
    };
    // J = mean vertical displacement of the top edge.
    let n = mesh.nodes.len();
    let mut j = vec![0.0f64; 2 * n];
    let mut cnt = 0.0;
    for (i, p) in mesh.nodes.iter().enumerate() {
        if (p[1] - 1.0).abs() < 1e-12 {
            j[2 * i + 1] = 1.0;
            cnt += 1.0;
        }
    }
    for v in &mut j {
        *v /= cnt;
    }
    let (sol0, problem0) = solve_at(0.0);
    // Floor translated UP by h means the surface y = −0.02 + h moves
    // with +e_y... the offset enters as φ = y + 0.02 − h, so the
    // translation direction is e = (0, 1).
    let dj_adj = problem0
        .translation_gradient(&sol0, 1.0, &j, [0.0, 1.0])
        .expect("adjoint");
    let delta = 1e-6;
    let jval = |sol: &fs_solid::contact::ContactSolution| -> f64 {
        j.iter().zip(&sol.u).map(|(a, b)| a * b).sum()
    };
    let (sp, _) = solve_at(delta);
    let (sm, _) = solve_at(-delta);
    let dj_fd = (jval(&sp) - jval(&sm)) / (2.0 * delta);
    let rel = (dj_adj - dj_fd).abs() / dj_fd.abs().max(1e-30);
    verdict(
        "cnt-006-adjoint-gradient",
        rel < 1e-4,
        &format!("dJ/dh adjoint {dj_adj:.8e} vs FD {dj_fd:.8e} rel {rel:.2e}"),
    );
}
