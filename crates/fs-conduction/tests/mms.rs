//! G1 — manufactured solutions and convergence-order verification for
//! the steady conduction solve, gated by `fs_mms::OrderGate`.
//!
//! Every case runs the same ladder shape: Kuhn/Freudenthal subdivisions
//! of the unit cube at `n = 4, 6, 8, 10` (so `h = 1/n` is strictly
//! decreasing over four rungs — a two-point slope would hide a
//! pre-asymptotic bend, which is why `fs-mms` refuses ladders shorter
//! than three).
//!
//! **What a passing gate means and does not mean.** It means: on THESE
//! meshes, for THESE coefficients and data, the fitted log-log slope of
//! the error sits within `fs_mms::ORDER_GATE_TOLERANCE` of the P₁
//! theoretical rate. It is an OBSERVED order. It is not a proof that the
//! scheme attains that rate on other meshes, other coefficients, or
//! non-smooth data, and it is not an error bound on any particular
//! solve.
//!
//! Error norms come from an EXACT quadrature (see `support`), so the
//! fitted slopes carry no quadrature artefact.
//!
//! The two Dirichlet ladders use a QUARTIC manufactured solution on
//! purpose. P₁ on these meshes reproduces quadratic and cubic solutions
//! at the nodes exactly (pinned in `conformance.rs`), so a ladder run on
//! one of those measures the interpolation error `‖T − I_h T‖` — a
//! number independent of `K`. The quartic ladders fit the scheme's own
//! discretization error; the Neumann and Robin ladders stay on lower-order
//! solutions because those keep the BOUNDARY DATA exactly representable
//! in the P₁ trace space, which is what those two cases are there to
//! test. That trade is stated here rather than left to be discovered.

mod support;

use fs_conduction::assemble::DofMap;
use fs_conduction::bc::{ThermalBc, ThermalBoundaryBuilder};
use fs_conduction::field::ScalarField;
use fs_conduction::fixtures::{on_box_face, unit_cube};
use fs_conduction::material::{ConductivityModel, ConductivityTable};
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::solve::{
    ConductionProblem, InitialGuess, LinearConfig, Nonlinearity, SolveConfig, StopRule, solve,
};
use fs_mms::{Coverage, LadderSide, MmsMatrix, MmsMatrixRow, OrderGate, RefinementLadder};
use support::{FaceLinearQuadratic, FullQuadratic, Quartic, h1_error, l2_error, with_cx};

const GRIDS: [usize; 4] = [4, 6, 8, 10];
const ISOTROPIC_K: f64 = 2.5;
const ANISOTROPIC_K: [[f64; 3]; 3] = [[3.0, 0.5, 0.25], [0.5, 2.0, 0.75], [0.25, 0.75, 1.5]];

fn linear_config() -> SolveConfig {
    SolveConfig {
        // A linear material makes Picard with ω = 1 the exact one-shot
        // solve, and its operator is SPD, so the ladder runs on CG.
        nonlinearity: Nonlinearity::FixedPoint {
            relaxation: 1.0,
            max_backtracks: 8,
        },
        stop: StopRule {
            residual_rtol: 1e-11,
            residual_atol: 1e-24,
            step_atol: 0.0,
            max_iterations: 12,
        },
        linear: LinearConfig {
            tolerance: 1e-13,
            max_iterations: 40_000,
            restart: 60,
        },
        initial: InitialGuess::Uniform(300.0),
    }
}

fn nodal(mesh: &ConductionMesh, f: &dyn Fn([f64; 3]) -> f64) -> ScalarField {
    ScalarField::Nodal(mesh.positions().iter().map(|&p| f(p)).collect())
}

fn report(case: &str, side: LadderSide, theoretical: f64, hs: &[f64], errors: &[f64]) {
    let ladder = RefinementLadder::new(hs.to_vec(), errors.to_vec())
        .unwrap_or_else(|e| panic!("{case}: inadmissible ladder: {e}"));
    let gate = OrderGate { theoretical };
    match gate.check(case, side, &ladder) {
        Ok(verdict) => {
            println!("{}", verdict.json_line(true));
            println!(
                "{{\"mms\":\"ladder\",\"case\":\"{}\",\"h\":{hs:?},\"errors\":{errors:?}}}",
                support::json_escape(case)
            );
        }
        Err(e) => panic!("{case}: order gate refused: {e}\n  h = {hs:?}\n  err = {errors:?}"),
    }
}

// ---------------------------------------------------------------- case 1

fn run_isotropic_dirichlet(n: usize) -> (f64, f64) {
    let (complex, positions) = unit_cube(n);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let material = ConductivityModel::isotropic_declared(ISOTROPIC_K).expect("material");
    let k = material.tensor_at(0.0).expect("tensor");
    let source = nodal(&mesh, &|p| Quartic::source(k, p));
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "all",
            |_| true,
            ThermalBc::Dirichlet {
                temperature: nodal(&mesh, &Quartic::value),
            },
        )
        .expect("dirichlet region")
        .finish()
        .expect("boundary");
    let solution = with_cx(|cx| {
        solve(
            cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                source: &source,
            },
            linear_config(),
        )
        .expect("solve")
    });
    (
        l2_error(&mesh, &solution.temperature, &Quartic::value),
        h1_error(&mesh, &solution.temperature, &Quartic::gradient),
    )
}

#[test]
fn mms_isotropic_dirichlet_orders() {
    let mut hs = Vec::new();
    let mut l2 = Vec::new();
    let mut h1 = Vec::new();
    for &n in &GRIDS {
        let (e2, e1) = run_isotropic_dirichlet(n);
        hs.push(1.0 / n as f64);
        l2.push(e2);
        h1.push(e1);
    }
    report(
        "conduction/mms/isotropic-dirichlet/l2",
        LadderSide::Primal,
        2.0,
        &hs,
        &l2,
    );
    report(
        "conduction/mms/isotropic-dirichlet/h1",
        LadderSide::Primal,
        1.0,
        &hs,
        &h1,
    );
}

// ---------------------------------------------------------------- case 2

fn run_anisotropic_dirichlet(n: usize) -> f64 {
    let (complex, positions) = unit_cube(n);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let material = ConductivityModel::constant_tensor(ANISOTROPIC_K).expect("material");
    let source = nodal(&mesh, &|p| Quartic::source(ANISOTROPIC_K, p));
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "all",
            |_| true,
            ThermalBc::Dirichlet {
                temperature: nodal(&mesh, &Quartic::value),
            },
        )
        .expect("dirichlet region")
        .finish()
        .expect("boundary");
    let solution = with_cx(|cx| {
        solve(
            cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                source: &source,
            },
            linear_config(),
        )
        .expect("solve")
    });
    l2_error(&mesh, &solution.temperature, &Quartic::value)
}

#[test]
fn mms_anisotropic_dirichlet_order() {
    let mut hs = Vec::new();
    let mut l2 = Vec::new();
    for &n in &GRIDS {
        hs.push(1.0 / n as f64);
        l2.push(run_anisotropic_dirichlet(n));
    }
    report(
        "conduction/mms/anisotropic-dirichlet/l2",
        LadderSide::Primal,
        2.0,
        &hs,
        &l2,
    );
}

// ---------------------------------------------------------------- case 3

fn run_mixed_neumann(n: usize) -> f64 {
    let (complex, positions) = unit_cube(n);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let material = ConductivityModel::isotropic_declared(ISOTROPIC_K).expect("material");
    let k = material.tensor_at(0.0).expect("tensor");
    let source = ScalarField::Uniform(FullQuadratic::source(k));
    // Outward flux on x = 0 (n = −e_x) is +k ∂T/∂x; on x = 1 it is
    // −k ∂T/∂x. Both are LINEAR in the face coordinates, so the nodal
    // A(1+δ)/12 rule integrates them exactly.
    let flux_lo = nodal(&mesh, &|p| ISOTROPIC_K * FullQuadratic::gradient(p)[0]);
    let flux_hi = nodal(&mesh, &|p| -ISOTROPIC_K * FullQuadratic::gradient(p)[0]);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "x-lo",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::Neumann {
                outward_flux: flux_lo,
            },
        )
        .expect("x-lo")
        .region(
            "x-hi",
            |f| on_box_face(f.centroid[0], 1.0),
            ThermalBc::Neumann {
                outward_flux: flux_hi,
            },
        )
        .expect("x-hi")
        .remainder(
            "dirichlet",
            ThermalBc::Dirichlet {
                temperature: nodal(&mesh, &FullQuadratic::value),
            },
        )
        .expect("remainder")
        .finish()
        .expect("boundary");
    let solution = with_cx(|cx| {
        solve(
            cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                source: &source,
            },
            linear_config(),
        )
        .expect("solve")
    });
    l2_error(&mesh, &solution.temperature, &FullQuadratic::value)
}

#[test]
fn mms_mixed_neumann_order() {
    let mut hs = Vec::new();
    let mut l2 = Vec::new();
    for &n in &GRIDS {
        hs.push(1.0 / n as f64);
        l2.push(run_mixed_neumann(n));
    }
    report(
        "conduction/mms/mixed-neumann/l2",
        LadderSide::Primal,
        2.0,
        &hs,
        &l2,
    );
}

// ---------------------------------------------------------------- case 4

const ROBIN_HTC: f64 = 8.0;

fn run_robin(n: usize) -> f64 {
    let (complex, positions) = unit_cube(n);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let material = ConductivityModel::isotropic_declared(ISOTROPIC_K).expect("material");
    let source = ScalarField::Uniform(-ISOTROPIC_K * FaceLinearQuadratic::LAPLACIAN);
    // On x = 0 the outward flux is +k ∂T/∂x, and this manufactured
    // solution is LINEAR on that face, so T_ref = T − q_n/h is linear
    // too: the Robin data lives exactly in the P₁ trace space and the
    // ladder measures the discretization, not its own boundary data.
    let t_ref = nodal(&mesh, &|p| {
        FaceLinearQuadratic::value(p)
            - ISOTROPIC_K * FaceLinearQuadratic::gradient(p)[0] / ROBIN_HTC
    });
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "convective",
            |f| on_box_face(f.centroid[0], 0.0),
            ThermalBc::Robin {
                htc: ScalarField::Uniform(ROBIN_HTC),
                t_ref,
            },
        )
        .expect("robin region")
        .remainder(
            "dirichlet",
            ThermalBc::Dirichlet {
                temperature: nodal(&mesh, &FaceLinearQuadratic::value),
            },
        )
        .expect("remainder")
        .finish()
        .expect("boundary");
    let solution = with_cx(|cx| {
        solve(
            cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                source: &source,
            },
            linear_config(),
        )
        .expect("solve")
    });
    l2_error(&mesh, &solution.temperature, &FaceLinearQuadratic::value)
}

#[test]
fn mms_robin_order() {
    let mut hs = Vec::new();
    let mut l2 = Vec::new();
    for &n in &GRIDS {
        hs.push(1.0 / n as f64);
        l2.push(run_robin(n));
    }
    report("conduction/mms/robin/l2", LadderSide::Primal, 2.0, &hs, &l2);
}

// ---------------------------------------------------------------- case 5

/// `k(T) = k₀ (1 + β (T − T₀))` — LINEAR in T, so the two-knot
/// piecewise-linear table reproduces it exactly and the ladder measures
/// the nonlinear discretization rather than a table-interpolation error.
const K0: f64 = 2.5;
const BETA: f64 = 0.01;
const T0: f64 = 300.0;

fn conductivity_of(t: f64) -> f64 {
    K0 * BETA.mul_add(t - T0, 1.0)
}

fn nonlinear_material() -> ConductivityModel {
    ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![
            (280.0, conductivity_of(280.0)),
            (340.0, conductivity_of(340.0)),
        ])
        .expect("curve"),
    )
}

fn nonlinear_source(p: [f64; 3]) -> f64 {
    let g = FaceLinearQuadratic::gradient(p);
    let grad2 = g[0].mul_add(g[0], g[1].mul_add(g[1], g[2] * g[2]));
    let t = FaceLinearQuadratic::value(p);
    -(K0 * BETA).mul_add(grad2, conductivity_of(t) * FaceLinearQuadratic::LAPLACIAN)
}

fn run_nonlinear(n: usize) -> f64 {
    let (complex, positions) = unit_cube(n);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let material = nonlinear_material();
    let source = nodal(&mesh, &nonlinear_source);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "all",
            |_| true,
            ThermalBc::Dirichlet {
                temperature: nodal(&mesh, &FaceLinearQuadratic::value),
            },
        )
        .expect("dirichlet region")
        .finish()
        .expect("boundary");
    let config = SolveConfig {
        nonlinearity: Nonlinearity::default(),
        stop: StopRule {
            residual_rtol: 1e-11,
            residual_atol: 1e-24,
            step_atol: 0.0,
            max_iterations: 25,
        },
        linear: LinearConfig {
            tolerance: 1e-13,
            max_iterations: 40_000,
            restart: 60,
        },
        initial: InitialGuess::DirichletMean,
    };
    let solution = with_cx(|cx| {
        solve(
            cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                source: &source,
            },
            config,
        )
        .expect("nonlinear solve")
    });
    assert!(
        solution.report.iterations >= 2,
        "a k(T) problem that converges in one iteration is not exercising the nonlinearity"
    );
    l2_error(&mesh, &solution.temperature, &FaceLinearQuadratic::value)
}

#[test]
fn mms_nonlinear_conductivity_order() {
    let mut hs = Vec::new();
    let mut l2 = Vec::new();
    for &n in &GRIDS {
        hs.push(1.0 / n as f64);
        l2.push(run_nonlinear(n));
    }
    report(
        "conduction/mms/nonlinear-kt/l2",
        LadderSide::Primal,
        2.0,
        &hs,
        &l2,
    );
}

// ----------------------------------------------------------- the matrix

/// The declared battery matrix: coverage in data, gaps visible with a
/// reason. `fs-mms` exists so a coverage hole is lintable rather than
/// silently absent.
#[test]
fn mms_battery_matrix_is_declared() {
    let matrix = MmsMatrix {
        rows: vec![
            row(
                "p1-simplicial",
                "dirichlet-isotropic-quartic",
                Coverage::Covered {
                    test: "tests/mms.rs::mms_isotropic_dirichlet_orders".to_string(),
                },
            ),
            row(
                "p1-simplicial",
                "dirichlet-anisotropic-quartic",
                Coverage::Covered {
                    test: "tests/mms.rs::mms_anisotropic_dirichlet_order".to_string(),
                },
            ),
            row(
                "p1-simplicial",
                "neumann-mixed",
                Coverage::Covered {
                    test: "tests/mms.rs::mms_mixed_neumann_order".to_string(),
                },
            ),
            row(
                "p1-simplicial",
                "robin-convective",
                Coverage::Covered {
                    test: "tests/mms.rs::mms_robin_order".to_string(),
                },
            ),
            row(
                "p1-simplicial",
                "dirichlet-nonlinear-kt",
                Coverage::Covered {
                    test: "tests/mms.rs::mms_nonlinear_conductivity_order".to_string(),
                },
            ),
            row(
                "p1-simplicial",
                "adjoint-order",
                Coverage::Gap {
                    reason: "the adjoint ladder needs a dual manufactured solution for the \
                         QoI functional; the IFT gradient is verified against central \
                         finite differences in tests/adjoint.rs instead, which checks \
                         consistency but fits no dual convergence order"
                        .to_string(),
                },
            ),
            row(
                "p2-simplicial",
                "any",
                Coverage::Gap {
                    reason: "this crate discretizes the P1 (FEEC 0-form) space only; \
                         higher-order thermal elements are not built"
                        .to_string(),
                },
            ),
            row(
                "cut-p1",
                "any",
                Coverage::Gap {
                    reason: "the CutFEM thermal frontend is a separate bead; this crate is \
                         body-fitted only"
                        .to_string(),
                },
            ),
        ],
    };
    for line in matrix.json_lines() {
        println!("{line}");
    }
    let gaps = matrix.gaps();
    assert_eq!(gaps.len(), 3, "every declared gap must carry a reason");
    for gap in gaps {
        match &gap.coverage {
            Coverage::Gap { reason } => assert!(!reason.is_empty()),
            Coverage::Covered { .. } => unreachable!("filtered to gaps"),
        }
    }
}

fn row(family: &str, bc: &str, coverage: Coverage) -> MmsMatrixRow {
    MmsMatrixRow {
        frontend: "feec-body-fitted-conduction".to_string(),
        family: family.to_string(),
        bc: bc.to_string(),
        coverage,
    }
}

/// The reduced system is what the Krylov method actually sees, so the
/// ladder's claim rests on the elimination being consistent: the free
/// rows of `A·T_full − b` must be the same vector the reduced system
/// reports as its residual.
#[test]
fn dof_map_partitions_the_unit_cube() {
    let (complex, positions) = unit_cube(4);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("all", |_| true, ThermalBc::dirichlet(300.0).expect("bc"))
        .expect("region")
        .finish()
        .expect("boundary");
    let dofs = DofMap::new(&boundary, mesh.vertex_count()).expect("dofs");
    // 5³ vertices, 3³ interior.
    assert_eq!(mesh.vertex_count(), 125);
    assert_eq!(dofs.n(), 27);
    assert_eq!(dofs.fixed().len(), 98);
    assert!(dofs.fixed().windows(2).all(|w| w[0] < w[1]));
}
