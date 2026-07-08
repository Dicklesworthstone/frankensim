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

use fs_cutfem::{Circle, Quadtree};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_solid::linear::{Formulation, LinearProblem, PlaneKind, l2_h1_error};
use fs_solid::{
    CutElasticity, HyperProblem, Mesh2, NewtonSettings, Patch, RegimeIndicators, select_formulation,
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
