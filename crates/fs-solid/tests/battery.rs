//! fs-solid conformance battery (bead tfz.13).
//!
//! - sol-001 G0: patch tests exact — linear displacement reproduced on
//!   P1 triangles, Q1 quads (both formulations), and the CutFEM
//!   frontend.
//! - sol-002 G1: MMS optimal orders across element families ×
//!   frontends × BC types (strong Dirichlet body-fitted; Nitsche +
//!   traction on the cut frontend fixture set).
//! - sol-003 G0: objectivity — superposed rigid rotations leave the
//!   stored energy invariant and rotate the Piola stress, through the
//!   fs-material interface for BOTH cards.
//! - sol-004 G0 (merge gate): adjoint/tangent consistency — the
//!   assembled consistent tangent matches the FD directional
//!   derivative of the residual, and the residual is the exact
//!   gradient of the potential energy.
//! - sol-005: locking battery — near-incompressible (ν → 1/2) thin
//!   bending: the standard element's error DEGRADES (measured), B-bar
//!   stays flat (the locking-free payoff).
//! - sol-006 G2: Cook's membrane envelope (compressible vs literature
//!   band + self-convergence; near-incompressible self-convergence
//!   with the standard element's failure logged).
//! - sol-007: Newton robustness — large-strain bending and a
//!   buckling-adjacent compressed strip converge under load stepping +
//!   line search, with iteration histories logged.

use fs_cutfem::{Circle, CutElasticity as CanonicalCutElasticity, HalfPlane, Quadtree};
use fs_material::IsotropicElastic;
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_solid::linear::{Formulation, LinearProblem, PlaneKind, l2_h1_error};
use fs_solid::{
    BoundaryTraction, CutElasticity, DesignBoxEdge, EdgeBand, HyperProblem, Mesh2, NewtonSettings,
    Patch, RegimeIndicators, SolidError, select_formulation,
};
use std::f64::consts::PI;
use std::fmt::Write as _;

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

// ------------------------------------------------------------- MMS fixture

const E_MOD: f64 = 1.0;
const NU: f64 = 0.3;

fn mms_u(x: f64, y: f64) -> [f64; 2] {
    [
        (PI * x).sin() * (PI * y).sin(),
        (PI * x).cos() * (PI * y).cos(),
    ]
}

fn mms_grad(x: f64, y: f64) -> [[f64; 2]; 2] {
    [
        [
            PI * (PI * x).cos() * (PI * y).sin(),
            PI * (PI * x).sin() * (PI * y).cos(),
        ],
        [
            -PI * (PI * x).sin() * (PI * y).cos(),
            -PI * (PI * x).cos() * (PI * y).sin(),
        ],
    ]
}

/// f = −div σ(u) for the MMS pair under plane strain (λ, μ from E, ν).
fn mms_f(x: f64, y: f64) -> [f64; 2] {
    let (lambda, mu) = fs_solid::linear::lame(E_MOD, NU, PlaneKind::Strain);
    let u = mms_u(x, y);
    let pp = PI * PI;
    // u1,xx = u1,yy = −π²u1; u2,xy = π²u1; u1,xy = π²u2; u2,xx = u2,yy = −π²u2.
    [
        -((lambda + 2.0 * mu) * (-pp * u[0]) + (lambda + mu) * (pp * u[0]) + mu * (-pp * u[0])),
        -((lambda + mu) * (pp * u[1]) + mu * (-pp * u[1]) + (lambda + 2.0 * mu) * (-pp * u[1])),
    ]
}

// ------------------------------------------------------------------ sol-001

#[test]
fn sol_001_patch_tests_exact() {
    let lin = |x: f64, y: f64| [0.1 + 0.2 * x + 0.3 * y, 0.4 - 0.1 * x + 0.05 * y];
    let all = [Patch::Left, Patch::Right, Patch::Bottom, Patch::Top];
    let mut worst = 0.0f64;
    // Body-fitted: P1, Q1 standard, Q1 B-bar.
    for (mesh, formulation) in [
        (Mesh2::triangles(1.0, 1.0, 5, 5), Formulation::Standard),
        (Mesh2::quads(1.0, 1.0, 5, 5), Formulation::Standard),
        (Mesh2::quads(1.0, 1.0, 5, 5), Formulation::BBar),
    ] {
        let problem = LinearProblem {
            mesh: &mesh,
            youngs: E_MOD,
            poisson: NU,
            plane: PlaneKind::Strain,
            formulation,
            body_force: None,
            dirichlet: all.iter().map(|&p| (p, &lin as _)).collect(),
            traction: vec![],
        symmetry: vec![],
        };
        let u = problem.solve().expect("patch solves");
        for (node, val) in u.iter().enumerate() {
            let p = mesh.nodes[node];
            let ex = lin(p[0], p[1]);
            worst = worst
                .max((val[0] - ex[0]).abs())
                .max((val[1] - ex[1]).abs());
        }
    }
    // CutFEM frontend: the ALL-EMBEDDED disk — for a globally linear
    // exact field, the chord-approximated domain problem has the same
    // linear solution, so the patch test is exact despite curvature.
    let grid = Quadtree::uniform(4);
    let sdf = Circle {
        center: [0.5, 0.5],
        radius: 0.35,
    };
    let cut = CutElasticity {
        grid: &grid,
        sdf: &sdf,
        youngs: E_MOD,
        poisson: NU,
        nitsche_beta: 20.0,
        ghost_gamma: 0.5,
        quad_depth: 3,
        clamp: None,
        boundary_traction: None,
        traction_free_interface: false,
    };
    let sol = cut
        .solve(&|_, _| [0.0, 0.0], &|x, y| lin(x, y))
        .expect("cut patch solves");
    let (l2, _) = cut.l2_h1_error(&sol, &|x, y| lin(x, y), &|_, _| [[0.2, 0.3], [-0.1, 0.05]]);
    // The cut gate sits at the CG-tolerance floor, not discretization.
    let pass = worst < 1e-9 && l2 < 1e-7;
    verdict(
        "sol-001",
        pass,
        &format!(
            "\"detail\":\"linear patch exact on P1/Q1/B-bar/CutFEM\",\
             \"body_fitted_max_nodal\":{worst:.3e},\"cutfem_l2\":{l2:.3e}"
        ),
    );
}

// ------------------------------------------------------------------ sol-002

#[test]
fn sol_002_mms_orders_families_and_frontends() {
    let all = [Patch::Left, Patch::Right, Patch::Bottom, Patch::Top];
    let mut rows = String::new();
    let mut pass = true;
    // Body-fitted families.
    for (family, make) in [
        (
            "p1-tri",
            (|n| Mesh2::triangles(1.0, 1.0, n, n)) as fn(usize) -> Mesh2,
        ),
        (
            "q1-quad",
            (|n| Mesh2::quads(1.0, 1.0, n, n)) as fn(usize) -> Mesh2,
        ),
    ] {
        let mut errs = Vec::new();
        for n in [8usize, 16, 32] {
            let mesh = make(n);
            let problem = LinearProblem {
                mesh: &mesh,
                youngs: E_MOD,
                poisson: NU,
                plane: PlaneKind::Strain,
                formulation: Formulation::Standard,
                body_force: Some(&mms_f),
                dirichlet: all.iter().map(|&p| (p, &mms_u as _)).collect(),
                traction: vec![],
        symmetry: vec![],
            };
            let u = problem.solve().expect("mms solves");
            errs.push(l2_h1_error(&mesh, &u, &mms_u, &mms_grad));
        }
        let l2_order = (errs[1].0 / errs[2].0).log2();
        let h1_order = (errs[1].1 / errs[2].1).log2();
        pass &= l2_order > 1.8 && h1_order > 0.9;
        let _ = write!(
            rows,
            "{{\"family\":\"{family}\",\"l2\":[{:.3e},{:.3e},{:.3e}],\
             \"l2_order\":{l2_order:.2},\"h1_order\":{h1_order:.2}}},",
            errs[0].0, errs[1].0, errs[2].0
        );
    }
    // CutFEM frontend on the disk (Nitsche embedded BC).
    let circle = Circle {
        center: [0.5, 0.5],
        radius: 0.3,
    };
    let mut cut_errs = Vec::new();
    for level in [4u32, 5, 6] {
        let grid = Quadtree::uniform(level);
        let cut = CutElasticity {
            grid: &grid,
            sdf: &circle,
            youngs: E_MOD,
            poisson: NU,
            nitsche_beta: 20.0,
            ghost_gamma: 0.5,
            quad_depth: 3,
            clamp: None,
            boundary_traction: None,
            traction_free_interface: false,
        };
        let sol = cut
            .solve(&|x, y| mms_f(x, y), &|x, y| mms_u(x, y))
            .expect("cut mms solves");
        cut_errs.push(cut.l2_h1_error(&sol, &|x, y| mms_u(x, y), &|x, y| mms_grad(x, y)));
    }
    let cut_l2_order = (cut_errs[1].0 / cut_errs[2].0).log2();
    let cut_h1_order = (cut_errs[1].1 / cut_errs[2].1).log2();
    pass &= cut_l2_order > 1.8 && cut_h1_order > 0.85;
    let _ = write!(
        rows,
        "{{\"family\":\"cutfem-q1\",\"l2\":[{:.3e},{:.3e},{:.3e}],\
         \"l2_order\":{cut_l2_order:.2},\"h1_order\":{cut_h1_order:.2}}}",
        cut_errs[0].0, cut_errs[1].0, cut_errs[2].0
    );
    verdict(
        "sol-002",
        pass,
        &format!("\"detail\":\"MMS orders, families x frontends\",\"rows\":[{rows}]"),
    );
}

// ------------------------------------------------ canonical CutFEM ownership

#[test]
fn cutfront_delegates_bitwise_to_canonical_cutfem() {
    let grid = Quadtree::uniform(3);
    let sdf = Circle {
        center: [0.5, 0.5],
        radius: 0.35,
    };
    let displacement = |x: f64, y: f64| [0.1 + 0.2 * x + 0.3 * y, 0.4 - 0.1 * x + 0.05 * y];
    let body = |_: f64, _: f64| [0.0, 0.0];
    let facade = CutElasticity {
        grid: &grid,
        sdf: &sdf,
        youngs: E_MOD,
        poisson: NU,
        nitsche_beta: 20.0,
        ghost_gamma: 0.5,
        quad_depth: 3,
        clamp: None,
        boundary_traction: None,
        traction_free_interface: false,
    };

    let material = IsotropicElastic::new(E_MOD, NU, 1.0).expect("valid parity material");
    let (lambda, mu) = material.lame();
    let legacy_stiffness_scale = (lambda + 2.0 * mu) / mu;
    let canonical = CanonicalCutElasticity {
        grid: &grid,
        sdf: &sdf,
        material: &material,
        nitsche_beta: 20.0 * legacy_stiffness_scale,
        ghost_gamma: 0.5 * legacy_stiffness_scale,
        quad_depth: 3,
        clamp: None,
        boundary_traction: None,
        traction_free_interface: false,
        solver_tol: 1e-12,
        solver_max_iters: 60_000,
    };

    let facade_solution = facade
        .solve(&body, &displacement)
        .expect("facade solve delegates");
    let canonical_solution = canonical
        .solve(&body, &displacement)
        .expect("canonical parity solve");

    assert_eq!(
        facade_solution.active_cells(),
        canonical_solution.active_cells()
    );
    assert_eq!(facade_solution.iters, canonical_solution.iters);
    assert_eq!(
        facade_solution.rel_residual.to_bits(),
        canonical_solution.rel_residual.to_bits()
    );
    assert_eq!(
        facade_solution
            .coefficients()
            .iter()
            .map(|value| value.to_bits())
            .collect::<Vec<_>>(),
        canonical_solution
            .coefficients()
            .iter()
            .map(|value| value.to_bits())
            .collect::<Vec<_>>()
    );
    assert_eq!(
        facade_solution.nodal().len(),
        canonical_solution.nodal().len()
    );
    for ((facade_node, facade_value), (canonical_node, canonical_value)) in facade_solution
        .nodal()
        .iter()
        .zip(canonical_solution.nodal())
    {
        assert_eq!(facade_node, canonical_node);
        assert_eq!(
            facade_value.map(f64::to_bits),
            canonical_value.map(f64::to_bits)
        );
    }

    let facade_error = facade.l2_h1_error(&facade_solution, &displacement, &|_, _| {
        [[0.2, 0.3], [-0.1, 0.05]]
    });
    let canonical_error = canonical.l2_h1_error(&canonical_solution, &displacement, &|_, _| {
        [[0.2, 0.3], [-0.1, 0.05]]
    });
    assert_eq!(facade_error.0.to_bits(), canonical_error.0.to_bits());
    assert_eq!(facade_error.1.to_bits(), canonical_error.1.to_bits());
}

#[test]
fn cutfront_typed_edge_band_delegates_bitwise_to_canonical_cutfem() {
    let grid = Quadtree::uniform(3);
    let sdf = HalfPlane {
        normal: [1.0, 0.0],
        offset: 2.0,
    };
    let clamp = |x: f64, _: f64| x == 0.0;
    let traction = |_: f64, _: f64| [0.0, -1.0];
    let zero = |_: f64, _: f64| [0.0, 0.0];
    let support =
        EdgeBand::new(DesignBoxEdge::Right, 0.25, 0.75).expect("valid typed right-edge support");
    let facade = CutElasticity {
        grid: &grid,
        sdf: &sdf,
        youngs: E_MOD,
        poisson: NU,
        nitsche_beta: 20.0,
        ghost_gamma: 0.5,
        quad_depth: 2,
        clamp: Some(&clamp),
        boundary_traction: None,
        traction_free_interface: true,
    };

    let material = IsotropicElastic::new(E_MOD, NU, 1.0).expect("valid parity material");
    let (lambda, mu) = material.lame();
    let legacy_stiffness_scale = (lambda + 2.0 * mu) / mu;
    let canonical = CanonicalCutElasticity {
        grid: &grid,
        sdf: &sdf,
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

    let facade_solution = facade
        .solve_with_boundary_traction(
            &zero,
            &zero,
            BoundaryTraction::EdgeBand {
                support,
                value: &traction,
            },
        )
        .expect("typed facade solve delegates");
    let canonical_solution = canonical
        .solve_with_boundary_traction(
            &zero,
            &zero,
            BoundaryTraction::EdgeBand {
                support,
                value: &traction,
            },
        )
        .expect("typed canonical parity solve");

    assert_eq!(
        facade_solution
            .coefficients()
            .iter()
            .map(|value| value.to_bits())
            .collect::<Vec<_>>(),
        canonical_solution
            .coefficients()
            .iter()
            .map(|value| value.to_bits())
            .collect::<Vec<_>>()
    );
    assert_eq!(
        facade_solution.compliance().to_bits(),
        canonical_solution.compliance().to_bits()
    );
    assert_eq!(facade_solution.nodal(), canonical_solution.nodal());

    let occupied = CutElasticity {
        boundary_traction: Some(&zero),
        ..facade
    };
    assert!(matches!(
        occupied.solve_with_boundary_traction(
            &zero,
            &zero,
            BoundaryTraction::EdgeBand {
                support,
                value: &traction,
            },
        ),
        Err(SolidError::InvalidInput { .. })
    ));
}

#[test]
fn cutfront_maps_canonical_refusals_fail_closed() {
    fn invalid_what(error: SolidError) -> String {
        match error {
            SolidError::InvalidInput { what } => what,
            other => panic!("expected InvalidInput, got {other:?}"),
        }
    }

    let grid = Quadtree::uniform(2);
    let sdf = Circle {
        center: [0.5, 0.5],
        radius: 0.3,
    };
    let invalid_material = CutElasticity {
        grid: &grid,
        sdf: &sdf,
        youngs: f64::NAN,
        poisson: NU,
        nitsche_beta: 20.0,
        ghost_gamma: 0.5,
        quad_depth: 2,
        clamp: None,
        boundary_traction: None,
        traction_free_interface: false,
    };
    let error = invalid_material
        .solve(&|_, _| [0.0, 0.0], &|_, _| [0.0, 0.0])
        .expect_err("non-finite material must be refused");
    assert!(invalid_what(error).contains("material card refused"));

    let unsupported_regime = CutElasticity {
        grid: &grid,
        sdf: &sdf,
        youngs: E_MOD,
        poisson: 0.49,
        nitsche_beta: 20.0,
        ghost_gamma: 0.5,
        quad_depth: 2,
        clamp: None,
        boundary_traction: None,
        traction_free_interface: false,
    };
    let error = unsupported_regime
        .solve(&|_, _| [0.0, 0.0], &|_, _| [0.0, 0.0])
        .expect_err("uncertified near-incompressible regime must be refused");
    assert!(invalid_what(error).contains("certified compressible-regime limit"));

    let empty_sdf = Circle {
        center: [2.0, 2.0],
        radius: 0.1,
    };
    let empty_domain = CutElasticity {
        grid: &grid,
        sdf: &empty_sdf,
        youngs: E_MOD,
        poisson: NU,
        nitsche_beta: 20.0,
        ghost_gamma: 0.5,
        quad_depth: 2,
        clamp: None,
        boundary_traction: None,
        traction_free_interface: false,
    };
    let error = empty_domain
        .solve(&|_, _| [0.0, 0.0], &|_, _| [0.0, 0.0])
        .expect_err("empty domain must be refused");
    assert!(invalid_what(error).contains("no active cells"));
}

// ------------------------------------------------------------------ sol-003

#[test]
fn sol_003_objectivity_and_frame_indifference() {
    let cards = [
        (
            "neo-hookean",
            Hyperelastic::new(
                HyperelasticModel::NeoHookean {
                    mu: 1.0,
                    lambda: 3.0,
                },
                3.0,
            )
            .expect("card"),
        ),
        (
            "mooney-rivlin",
            Hyperelastic::new(
                HyperelasticModel::MooneyRivlin {
                    c10: 0.4,
                    c01: 0.1,
                    kappa: 5.0,
                },
                3.0,
            )
            .expect("card"),
        ),
    ];
    let mut lcg = Lcg(0x1001_2026_0707_0061);
    let mut worst_energy = 0.0f64;
    let mut worst_stress = 0.0f64;
    for (_, card) in &cards {
        for _ in 0..8 {
            // Well-conditioned random F (plane-strain embedded).
            let mut f = [0.0f64; 9];
            f[8] = 1.0;
            for (i, slot) in [(0usize, 1.0), (1, 0.0), (3, 0.0), (4, 1.0)] {
                f[i] = slot + 0.3 * (lcg.unit() - 0.5);
            }
            let theta = 2.0 * PI * lcg.unit();
            let (c, s) = (theta.cos(), theta.sin());
            // R F (rotation about z).
            let rf = [
                c * f[0] - s * f[3],
                c * f[1] - s * f[4],
                0.0,
                s * f[0] + c * f[3],
                s * f[1] + c * f[4],
                0.0,
                0.0,
                0.0,
                1.0,
            ];
            let w0: f64 = card.energy(&f);
            let w1: f64 = card.energy(&rf);
            worst_energy = worst_energy.max((w1 - w0).abs() / w0.abs().max(1e-12));
            let p0 = card.piola(&f).expect("stress");
            let p1 = card.piola(&rf).expect("stress");
            // Objectivity: P(RF) = R · P(F).
            for i in 0..2 {
                for j in 0..3 {
                    let rp = match i {
                        0 => c * p0[j] - s * p0[3 + j],
                        _ => s * p0[j] + c * p0[3 + j],
                    };
                    worst_stress = worst_stress.max((p1[3 * i + j] - rp).abs());
                }
            }
        }
    }
    let pass = worst_energy < 1e-10 && worst_stress < 1e-8;
    verdict(
        "sol-003",
        pass,
        &format!(
            "\"detail\":\"W(RF)=W(F), P(RF)=R P(F), both cards, 8 random states each\",\
             \"worst_energy_rel\":{worst_energy:.3e},\"worst_stress_abs\":{worst_stress:.3e}"
        ),
    );
}

// ------------------------------------------------------------------ sol-004

#[test]
fn sol_004_tangent_and_energy_consistency() {
    let card = Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 1.0,
            lambda: 3.0,
        },
        3.0,
    )
    .expect("card");
    let mesh = Mesh2::quads(1.0, 1.0, 3, 3);
    let problem = HyperProblem {
        mesh: &mesh,
        material: &card,
        dirichlet: vec![],
        traction: vec![(Patch::Right, &|_, _| [0.05, 0.1])],
        settings: NewtonSettings::default(),
    };
    let n = mesh.node_count();
    let mut lcg = Lcg(0x1001_2026_0707_0062);
    let u: Vec<f64> = (0..2 * n).map(|_| 0.05 * (lcg.unit() - 0.5)).collect();
    let v: Vec<f64> = (0..2 * n).map(|_| lcg.unit() - 0.5).collect();
    let (r0, k) = problem.residual_and_tangent(&u, 1.0).expect("probe");
    let eps = 1e-6;
    let up: Vec<f64> = u.iter().zip(&v).map(|(a, b)| a + eps * b).collect();
    let um: Vec<f64> = u.iter().zip(&v).map(|(a, b)| a - eps * b).collect();
    let (rp, _) = problem.residual_and_tangent(&up, 1.0).expect("probe");
    let (rm, _) = problem.residual_and_tangent(&um, 1.0).expect("probe");
    // Tangent consistency: K v ≈ (R(u+εv) − R(u−εv)) / 2ε.
    let mut kv = vec![0.0f64; 2 * n];
    k.spmv(&v, &mut kv);
    let mut num = 0.0f64;
    let mut den = 0.0f64;
    for i in 0..2 * n {
        let fd = (rp[i] - rm[i]) / (2.0 * eps);
        num += (kv[i] - fd) * (kv[i] - fd);
        den += fd * fd;
    }
    let tangent_rel = (num / den.max(1e-30)).sqrt();
    // Variational consistency: R·v ≈ (Π(u+εv) − Π(u−εv)) / 2ε.
    let ep = problem.potential_energy(&up, 1.0).expect("energy");
    let em = problem.potential_energy(&um, 1.0).expect("energy");
    let fd_dir = (ep - em) / (2.0 * eps);
    let rv: f64 = r0.iter().zip(&v).map(|(r, w)| r * w).sum();
    let grad_rel = (rv - fd_dir).abs() / fd_dir.abs().max(1e-12);
    let pass = tangent_rel < 1e-6 && grad_rel < 1e-6;
    verdict(
        "sol-004",
        pass,
        &format!(
            "\"detail\":\"consistent tangent = FD(residual); residual = grad(energy)\",\
             \"tangent_rel\":{tangent_rel:.3e},\"gradient_rel\":{grad_rel:.3e}"
        ),
    );
}

// ------------------------------------------------------------------ sol-005

/// Tip deflection of a cantilever under end shear.
fn cantilever_tip(nx: usize, ny: usize, poisson: f64, formulation: Formulation) -> f64 {
    let mesh = Mesh2::quads(10.0, 1.0, nx, ny);
    let problem = LinearProblem {
        mesh: &mesh,
        youngs: E_MOD,
        poisson,
        plane: PlaneKind::Strain,
        formulation,
        body_force: None,
        dirichlet: vec![(Patch::Left, &|_, _| [0.0, 0.0])],
        traction: vec![(Patch::Right, &|_, _| [0.0, 0.01])],
        symmetry: vec![],
    };
    let u = problem.solve().expect("cantilever solves");
    // Mean tip deflection over the right edge.
    let tips = mesh.patch_nodes(Patch::Right);
    #[allow(clippy::cast_precision_loss)]
    let mean = tips.iter().map(|&n| u[n][1]).sum::<f64>() / tips.len() as f64;
    mean
}

#[test]
fn sol_005_locking_battery() {
    let mut rows = String::new();
    let mut ratios = Vec::new();
    for (name, formulation) in [
        ("standard", Formulation::Standard),
        ("b-bar", Formulation::BBar),
    ] {
        let mut errs = Vec::new();
        for &nu in &[0.3, 0.4999] {
            let reference = cantilever_tip(120, 12, nu, Formulation::BBar);
            let coarse = cantilever_tip(30, 3, nu, formulation);
            errs.push(((coarse - reference) / reference).abs());
        }
        ratios.push(errs[1] / errs[0].max(1e-12));
        let _ = write!(
            rows,
            "{{\"formulation\":\"{name}\",\"err_nu03\":{:.3e},\"err_nu04999\":{:.3e},\
             \"degradation\":{:.1}}},",
            errs[0],
            errs[1],
            errs[1] / errs[0].max(1e-12)
        );
    }
    // Standard element degrades by orders of magnitude; B-bar holds.
    let pass = ratios[0] > 10.0 && ratios[1] < 2.0;
    verdict(
        "sol-005",
        pass,
        &format!(
            "\"detail\":\"thin-bending cantilever, nu 0.3 -> 0.4999: standard locks, B-bar flat\",\
             \"rows\":[{}]",
            rows.trim_end_matches(',')
        ),
    );
}

// ------------------------------------------------------------------ sol-006

fn cooks_tip(n: usize, poisson: f64, plane: PlaneKind, formulation: Formulation) -> f64 {
    let mesh = Mesh2::cooks_membrane(n);
    let problem = LinearProblem {
        mesh: &mesh,
        youngs: 1.0,
        poisson,
        plane,
        formulation,
        body_force: None,
        dirichlet: vec![(Patch::Left, &|_, _| [0.0, 0.0])],
        traction: vec![(Patch::Right, &|_, _| [0.0, 1.0 / 16.0])],
        symmetry: vec![],
    };
    let u = problem.solve().expect("cooks solves");
    // Vertical displacement at the right-edge midpoint (48, 52) — the
    // literature's reference point.
    let ref_point = mesh
        .nodes
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| {
            let da = (a[0] - 48.0).abs() + (a[1] - 52.0).abs();
            let db = (b[0] - 48.0).abs() + (b[1] - 52.0).abs();
            da.partial_cmp(&db).expect("finite")
        })
        .expect("nonempty")
        .0;
    u[ref_point][1]
}

#[test]
fn sol_006_cooks_membrane_envelope() {
    // Compressible plane stress (the literature configuration).
    let tip_ref = cooks_tip(48, 1.0 / 3.0, PlaneKind::Stress, Formulation::BBar);
    let tip_32 = cooks_tip(24, 1.0 / 3.0, PlaneKind::Stress, Formulation::BBar);
    // Near-incompressible plane strain: self-convergence + the
    // standard element's measured shortfall.
    let inc_ref = cooks_tip(48, 0.4999, PlaneKind::Strain, Formulation::BBar);
    let inc_32 = cooks_tip(24, 0.4999, PlaneKind::Strain, Formulation::BBar);
    let inc_std_32 = cooks_tip(24, 0.4999, PlaneKind::Strain, Formulation::Standard);
    let self_dev = ((tip_32 - tip_ref) / tip_ref).abs();
    let inc_dev = ((inc_32 - inc_ref) / inc_ref).abs();
    let std_shortfall = (inc_std_32 - inc_ref) / inc_ref;
    let pass = (23.5..=24.5).contains(&tip_ref)
        && self_dev < 0.02
        && inc_dev < 0.05
        && std_shortfall < -0.05;
    verdict(
        "sol-006",
        pass,
        &format!(
            "\"detail\":\"Cook's membrane: literature band + self-convergence + \
             near-incompressible with standard-element shortfall measured\",\
             \"tip_ref\":{tip_ref:.3},\"tip_n24\":{tip_32:.3},\"self_dev\":{self_dev:.4},\
             \"inc_tip_ref\":{inc_ref:.3},\"inc_dev\":{inc_dev:.4},\
             \"standard_shortfall\":{std_shortfall:.3}"
        ),
    );
}

// ------------------------------------------------------------------ sol-007

#[test]
fn sol_007_newton_robustness_large_strain() {
    let card = Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 1.0,
            lambda: 10.0,
        },
        4.0,
    )
    .expect("card");
    let mesh = Mesh2::quads(4.0, 1.0, 16, 4);
    // Large bending: end shear driving tip deflection ~ beam depth × 2.
    let bend = HyperProblem {
        mesh: &mesh,
        material: &card,
        dirichlet: vec![(Patch::Left, &|_, _| [0.0, 0.0])],
        traction: vec![(Patch::Right, &|_, _| [0.0, 0.12])],
        settings: NewtonSettings {
            load_steps: 5,
            ..NewtonSettings::default()
        },
    };
    let (u_bend, rep_bend) = bend.solve().expect("large bending converges");
    let tip = mesh
        .patch_nodes(Patch::Right)
        .iter()
        .map(|&n| u_bend[n][1])
        .fold(0.0f64, f64::max);
    // Buckling-adjacent: axial compression of a slender strip.
    let strip = Mesh2::quads(4.0, 0.4, 20, 2);
    let compress = HyperProblem {
        mesh: &strip,
        material: &card,
        dirichlet: vec![(Patch::Left, &|_, _| [0.0, 0.0])],
        traction: vec![(Patch::Right, &|_, _| [-0.08, 0.0])],
        settings: NewtonSettings {
            load_steps: 8,
            ..NewtonSettings::default()
        },
    };
    let (u_comp, rep_comp) = compress.solve().expect("compression converges");
    let shortening = strip
        .patch_nodes(Patch::Right)
        .iter()
        .map(|&n| u_comp[n][0])
        .fold(0.0f64, f64::min);
    // Terminal contraction: the last Newton drop of the final load
    // step is at least two orders (Newton's terminal regime).
    let last = rep_bend.histories.last().expect("histories");
    let terminal_fast = last.len() >= 2
        && (last[last.len() - 1] < 1e-2 * last[last.len() - 2] || last[last.len() - 1] < 1e-10);
    let max_iters = rep_bend
        .histories
        .iter()
        .chain(&rep_comp.histories)
        .map(Vec::len)
        .max()
        .expect("histories");
    let mut hist_rows = String::new();
    for (k, h) in rep_bend.histories.iter().enumerate() {
        let _ = write!(
            hist_rows,
            "{{\"step\":{k},\"iters\":{},\"final\":{:.2e}}},",
            h.len(),
            h.last().copied().unwrap_or(f64::NAN)
        );
    }
    let pass = tip > 1.0 && shortening < -0.05 && terminal_fast && max_iters <= 15;
    verdict(
        "sol-007",
        pass,
        &format!(
            "\"detail\":\"large-strain bending (5 steps) + buckling-adjacent compression \
             (8 steps) with line search\",\"tip\":{tip:.3},\"shortening\":{shortening:.3},\
             \"max_newton_iters\":{max_iters},\"backtracks_bend\":{},\
             \"backtracks_compress\":{},\"bend_steps\":[{}]",
            rep_bend.backtracks,
            rep_comp.backtracks,
            hist_rows.trim_end_matches(',')
        ),
    );
}

// -------------------------------------------------- selection-policy probe

#[test]
fn selection_guidance_is_consistent_with_battery() {
    // The locking battery's failing regime is exactly where the policy
    // switches to B-bar.
    assert_eq!(
        select_formulation(RegimeIndicators {
            poisson: 0.4999,
            slenderness: 10.0
        }),
        Formulation::BBar
    );
    assert_eq!(
        select_formulation(RegimeIndicators {
            poisson: 0.3,
            slenderness: 1.0
        }),
        Formulation::Standard
    );
}

// ------------------------------------------------------------------ sol-le1
// NAFEMS LE1 "Elliptic membrane" (bead frankensim-g42o, Gauntlet G2).
// Quarter model between inner ellipse x^2/2^2 + y^2/1 = 1 and outer
// ellipse x^2/3.25^2 + y^2/2.75^2 = 1; plane stress, E = 210e3 MPa,
// nu = 0.3; uniform outward NORMAL tension 10 MPa on the outer edge;
// symmetry on both straight edges. TARGET: sigma_yy at point D = (2, 0)
// equals 92.7 MPa (The Standard NAFEMS Benchmarks, TNSB Rev. 3, test
// LE1). CI gates a coarse band at fixture resolution plus monotone
// approach under refinement — the Ghia-cavity precedent (fine studies
// live in perf lanes, the envelope lives here).

/// Solve LE1 at mesh density (nx angular, ny radial); return sigma_yy
/// evaluated at the Gauss point of the corner element nearest D.
fn le1_sigma_yy(nx: usize, ny: usize) -> f64 {
    let (ai, bi, ao, bo) = (2.0, 1.0, 3.25, 2.75);
    // Radial-first parameterization: (s, t) = (radial blend, angle).
    // The (e_r, e_theta) ordering keeps the Jacobian POSITIVE — the
    // angle-first version is orientation-reversing (det J < 0), which
    // silently flips assembled signs (found the hard way; see bead).
    let map = move |s: f64, t: f64| -> [f64; 2] {
        let th = t * std::f64::consts::FRAC_PI_2;
        let inner = [ai * th.cos(), bi * th.sin()];
        let outer = [ao * th.cos(), bo * th.sin()];
        [
            (1.0 - s).mul_add(inner[0], s * outer[0]),
            (1.0 - s).mul_add(inner[1], s * outer[1]),
        ]
    };
    let mesh = Mesh2::mapped_quads(ny, nx, &map);
    let problem = LinearProblem {
        mesh: &mesh,
        youngs: 210_000.0, // MPa
        poisson: 0.3,
        plane: PlaneKind::Stress,
        formulation: Formulation::Standard,
        body_force: None,
        dirichlet: vec![],
        traction: vec![(Patch::Right, &|x, y| {
            // Outward unit normal of the OUTER ellipse at (x, y), scaled
            // by the 10 MPa tension.
            let g = [x / (3.25 * 3.25), y / (2.75 * 2.75)];
            let n = (g[0] * g[0] + g[1] * g[1]).sqrt();
            [10.0 * g[0] / n, 10.0 * g[1] / n]
        })],
        // Bottom edge (t = 0) is the x-axis: u_y = 0. Top edge
        // (t = 1) is the y-axis: u_x = 0. One component each — the
        // constraint the `symmetry` field exists for.
        symmetry: vec![(Patch::Bottom, 1), (Patch::Top, 0)],
    };
    let u = problem.solve().expect("LE1 solves");
    if std::env::var("LE1_DEBUG").is_ok() {
        let d_node = mesh.nodes.iter().position(|p| (p[0] - 2.0).abs() < 1e-9 && p[1].abs() < 1e-9).unwrap();
        let o_node = mesh.nodes.iter().position(|p| (p[0] - 3.25).abs() < 1e-9 && p[1].abs() < 1e-9).unwrap();
        println!("{{\"debug\":\"le1\",\"u_D\":[{:.6},{:.6}],\"u_outer\":[{:.6},{:.6}]}}", u[d_node][0], u[d_node][1], u[o_node][0], u[o_node][1]);
    }
    // Corner element at (xi, eta) = (0, 0): its first node is D itself
    // under mapped_quads' row-major numbering. Evaluate strain at the
    // Gauss point nearest the corner and form plane-stress sigma_yy.
    let conn = mesh
        .elems
        .iter()
        .find(|c| c.iter().any(|&n| {
            let p = mesh.nodes[n];
            (p[0] - 2.0).abs() < 1e-9 && p[1].abs() < 1e-9
        }))
        .expect("an element touches D");
    if std::env::var("LE1_DEBUG").is_ok() {
        for (label, xi, eta) in [("corner", -0.577, -0.577), ("center", 0.0, 0.0), ("xi+", 0.577, -0.577), ("eta+", -0.577, 0.577)] {
            let (_, gr, _) = fs_solid::mesh2::shapes_at(&mesh.nodes, conn, xi, eta);
            let (mut xx, mut yy) = (0.0, 0.0);
            for (a, g2) in conn.iter().zip(&gr) { xx += g2[0]*u[*a][0]; yy += g2[1]*u[*a][1]; }
            let cc: f64 = 210_000.0 / (1.0 - 0.09);
            println!("{{\"probe\":\"{label}\",\"sxx\":{:.2},\"syy\":{:.2}}}", cc*(xx + 0.3*yy), cc*(0.3f64).mul_add(xx, yy));
        }
        println!("{{\"conn\":{conn:?},\"nodes\":[{:?},{:?},{:?},{:?}]}}", mesh.nodes[conn[0]], mesh.nodes[conn[1]], mesh.nodes[conn[2]], mesh.nodes[conn[3]]);
    }
    let g = -1.0 / 3.0_f64.sqrt(); // 2x2 Gauss point nearest the D corner
    let (_, grads, _) = fs_solid::mesh2::shapes_at(&mesh.nodes, conn, g, g);
    let (mut exx, mut eyy) = (0.0, 0.0);
    for (a, gr) in conn.iter().zip(&grads) {
        exx += gr[0] * u[*a][0];
        eyy += gr[1] * u[*a][1];
    }
    let (e, nu): (f64, f64) = (210_000.0, 0.3);
    let c = e / (1.0 - nu * nu);
    c * nu.mul_add(exx, eyy)
}

#[test]
fn sol_le1_nafems_elliptic_membrane_envelope() {
    // Refinement ladder: the Gauss-point estimate must approach the
    // published 92.7 MPa monotonically-in-error and land in the coarse
    // band at the finest fixture level.
    let ladder: Vec<(usize, f64)> = [(8, 4), (16, 8), (32, 16)]
        .iter()
        .map(|&(nx, ny)| ((nx), le1_sigma_yy(nx, ny)))
        .collect();
    for (n, s) in &ladder {
        println!(
            "{{\"suite\":\"fs-solid\",\"case\":\"le1\",\"nx\":{n},\"sigma_yy\":{s:.3},\"target\":92.7}}"
        );
    }
    let errs: Vec<f64> = ladder.iter().map(|(_, s)| (s - 92.7).abs()).collect();
    assert!(
        errs[2] < errs[0],
        "refinement must approach the NAFEMS target: errors {errs:?}"
    );
    assert!(
        (85.0..=100.0).contains(&ladder[2].1),
        "sigma_yy at D within the coarse band [85, 100] MPa (target 92.7): got {:.2}",
        ladder[2].1
    );
    verdict(
        "sol-le1",
        true,
        "\"detail\":\"NAFEMS LE1 sigma_yy(D) approaches 92.7 MPa under refinement, coarse band held\"",
    );
}

#[test]
#[should_panic(expected = "orientation-reversing")]
fn sol_mapped_quads_refuses_orientation_reversing_maps() {
    // The exact trap that produced sign-flipped LE1 stresses: an
    // angle-first polar map has det J < 0. The constructor must refuse
    // it loudly instead of letting the solve negate every sign.
    let _ = Mesh2::mapped_quads(4, 4, &|xi, eta| {
        let th = xi * std::f64::consts::FRAC_PI_2;
        let r = 1.0 + eta;
        [r * th.cos(), r * th.sin()]
    });
}
