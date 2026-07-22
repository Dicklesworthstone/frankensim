//! fs-conduction conformance suite.
//!
//! Tiers exercised here:
//!
//! - **G0** — algebraic laws: operator symmetry and positive
//!   definiteness, the Dirichlet elimination identity, the energy-balance
//!   identity, the Newton Jacobian against central differences, and the
//!   typed refusal surface (every refusal in the error model that a
//!   caller can actually trigger).
//! - **G4** — cancellation drills: a cancelled assembly and a cancelled
//!   mid-iteration solve both return a STRUCTURED refusal naming the
//!   stage, never a partial operator.
//! - **G5** — determinism: a replayed solve is bitwise identical, and a
//!   run paused, sealed into `fs-exec`'s versioned envelope, restored,
//!   and resumed reproduces the uninterrupted trajectory bitwise.
//!
//! Every verdict prints one JSON line so a failure is reproducible from
//! its log line alone.

mod support;

use fs_blake3::hash_bytes;
use fs_conduction::adjoint::ConductivityDesign;
use fs_conduction::assemble::{
    DofMap, assemble_jacobian, assemble_operator, element_stiffness, reduce, residual,
};
use fs_conduction::bc::{ThermalBc, ThermalBoundaryBuilder};
use fs_conduction::field::ScalarField;
use fs_conduction::fixtures::{on_box_face, unit_cube};
use fs_conduction::material::{
    CONDUCTIVITY_DIMS, ConductivityModel, ConductivityTable, ProvenanceClass,
};
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::solve::{
    ConductionProblem, ConductionSolver, ConductionState, InitialGuess, LinearConfig, Nonlinearity,
    SolveConfig, StopRule, element_heat_flux, solve,
};
use fs_conduction::{ConductionError, TEMPERATURE_DIMS};
use fs_evidence::ValidityDomain;
use fs_exec::solver::SolverState;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, ObservationDataset, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{Dims, QtyAny};
use fs_rep_mesh::TetComplex;
use fs_scenario::bc::{BcKind, BcValue, BoundaryCondition, Physics};
use fs_scenario::scenario::Environment;
use support::{
    Cubic, FullQuadratic, Quartic, max_nodal_error, with_cancelled_cx, with_cx, with_gate,
};

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-conduction/conformance\",\"case\":\"{case}\",\
         \"verdict\":\"pass\",\"detail\":\"{}\"}}",
        support::json_escape(detail)
    );
}

const ISO_K: f64 = 2.5;
const ANISO_K: [[f64; 3]; 3] = [[3.0, 0.5, 0.25], [0.5, 2.0, 0.75], [0.25, 0.75, 1.5]];

fn nodal(mesh: &ConductionMesh, f: &dyn Fn([f64; 3]) -> f64) -> ScalarField {
    ScalarField::Nodal(mesh.positions().iter().map(|&p| f(p)).collect())
}

fn linear_config() -> SolveConfig {
    SolveConfig {
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

/// The reference fixture every determinism drill replays: an anisotropic
/// block with a Dirichlet face, a Neumann face, a convective face, and a
/// volumetric source.
struct Fixture {
    mesh: ConductionMesh,
    boundary: fs_conduction::bc::ThermalBoundary,
    material: ConductivityModel,
    source: ScalarField,
}

fn enclosure_block(n: usize) -> Fixture {
    let (complex, positions) = unit_cube(n);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let material = ConductivityModel::constant_tensor(ANISO_K).expect("material");
    let source = ScalarField::Uniform(5.0e3);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "cold-plate",
            |f| on_box_face(f.centroid[2], 0.0),
            ThermalBc::dirichlet(300.0).expect("dirichlet"),
        )
        .expect("cold plate")
        .region(
            "heated-face",
            |f| on_box_face(f.centroid[2], 1.0),
            ThermalBc::neumann(-1.2e3).expect("neumann"),
        )
        .expect("heated face")
        .region(
            "convective-side",
            |f| on_box_face(f.centroid[0], 1.0),
            ThermalBc::robin(25.0, 295.0).expect("robin"),
        )
        .expect("convective side")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    Fixture {
        mesh,
        boundary,
        material,
        source,
    }
}

impl Fixture {
    fn problem(&self) -> ConductionProblem<'_> {
        ConductionProblem {
            mesh: &self.mesh,
            boundary: &self.boundary,
            material: &self.material,
            source: &self.source,
        }
    }
}

// ------------------------------------------------------------------- G0

#[test]
fn assembled_operator_is_symmetric_and_positive_definite() {
    let fixture = enclosure_block(3);
    let dofs = DofMap::new(&fixture.boundary, fixture.mesh.vertex_count()).expect("dofs");
    let zero = vec![300.0f64; fixture.mesh.vertex_count()];
    let (a, _) = with_cx(|cx| {
        let system = assemble_operator(
            cx,
            &fixture.mesh,
            &fixture.boundary,
            &fixture.material,
            &fixture.source,
            &zero,
        )
        .expect("assemble");
        reduce(&system, &dofs)
    });

    let n = a.nrows();
    let mut worst_asymmetry = 0.0f64;
    let mut scale = 0.0f64;
    for i in 0..n {
        let (cols, vals) = a.row(i);
        for (&j, &v) in cols.iter().zip(vals) {
            scale = scale.max(v.abs());
            worst_asymmetry = worst_asymmetry.max((v - a.get(j, i)).abs());
        }
    }
    assert!(
        worst_asymmetry <= 1e-12 * scale,
        "operator asymmetry {worst_asymmetry:e} exceeds 1e-12 x {scale:e}"
    );

    // Positive definiteness on a deterministic spanning set of probes:
    // every unit vector, plus alternating-sign and ramp vectors.
    let mut smallest = f64::INFINITY;
    let mut probes: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let mut v = vec![0.0f64; n];
            v[i] = 1.0;
            v
        })
        .collect();
    probes.push(
        (0..n)
            .map(|i| if i % 2 == 0 { 1.0 } else { -1.0 })
            .collect(),
    );
    probes.push((0..n).map(|i| (i as f64 + 1.0) / n as f64).collect());
    probes.push(vec![1.0; n]);
    for v in &probes {
        let mut av = vec![0.0f64; n];
        a.spmv(v, &mut av);
        let q: f64 = v.iter().zip(&av).map(|(x, y)| x * y).sum();
        let norm: f64 = v.iter().map(|x| x * x).sum();
        smallest = smallest.min(q / norm);
    }
    assert!(
        smallest > 0.0,
        "Rayleigh quotient {smallest:e} is not positive: the conduction + Robin \
         operator must be SPD after Dirichlet elimination"
    );
    verdict(
        "operator-spd",
        &format!("n={n} asymmetry={worst_asymmetry:e} min_rayleigh={smallest:e}"),
    );
}

#[test]
fn dirichlet_elimination_matches_the_full_residual() {
    let fixture = enclosure_block(3);
    let solution = with_cx(|cx| solve(cx, fixture.problem(), linear_config()).expect("solve"));
    let dofs = DofMap::new(&fixture.boundary, fixture.mesh.vertex_count()).expect("dofs");
    // Prescribed values survive the round trip untouched.
    for &(v, value) in fixture.boundary.dirichlet() {
        assert_eq!(
            solution.temperature[v].to_bits(),
            value.to_bits(),
            "Dirichlet vertex {v} was not held exactly"
        );
    }
    let system = with_cx(|cx| {
        assemble_operator(
            cx,
            &fixture.mesh,
            &fixture.boundary,
            &fixture.material,
            &fixture.source,
            &solution.temperature,
        )
        .expect("assemble")
    });
    let r = residual(&system, &dofs, &solution.temperature);
    let worst = r.iter().fold(0.0f64, |a, b| a.max(b.abs()));
    let scale = system.load.iter().fold(0.0f64, |a, b| a.max(b.abs()));
    assert!(
        worst <= 1e-8 * scale.max(1.0),
        "free-row residual {worst:e} is not small against load scale {scale:e}"
    );
    verdict(
        "dirichlet-elimination",
        &format!("free_residual_inf={worst:e} load_scale={scale:e}"),
    );
}

/// P₁ Galerkin on these Kuhn meshes reproduces exact solutions up to
/// CUBIC at the nodes; a QUARTIC one it does not. That is a real
/// property of the scheme and it is also the reason a quadratic or cubic
/// MMS ladder measures the INTERPOLATION error, not the scheme's own
/// approximation error — a number independent of `K`. Pinning it here
/// keeps the G1 evidence honestly scoped instead of letting a degenerate
/// ladder read as extra strength, and it is why the Dirichlet ladders in
/// `tests/mms.rs` use the quartic solution.
#[test]
fn polynomial_reproduction_at_the_nodes() {
    let n = 4;
    let (complex, positions) = unit_cube(n);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let material = ConductivityModel::isotropic_declared(ISO_K).expect("material");
    let k = material.tensor_at(0.0).expect("tensor");
    let source = ScalarField::Uniform(FullQuadratic::source(k));
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "all",
            |_| true,
            ThermalBc::Dirichlet {
                temperature: nodal(&mesh, &FullQuadratic::value),
            },
        )
        .expect("region")
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
    let quadratic = max_nodal_error(&mesh, &solution.temperature, &FullQuadratic::value);
    assert!(
        quadratic < 1e-9,
        "expected nodal exactness for a quadratic solution, got {quadratic:e}"
    );

    // A CUBIC solution is reproduced at the nodes as well: the fourth
    // derivative that drives the truncation error vanishes.
    let source3 = nodal(&mesh, &|p| Cubic::source(k, p));
    let boundary3 = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "all",
            |_| true,
            ThermalBc::Dirichlet {
                temperature: nodal(&mesh, &Cubic::value),
            },
        )
        .expect("region")
        .finish()
        .expect("boundary");
    let solution3 = with_cx(|cx| {
        solve(
            cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary3,
                material: &material,
                source: &source3,
            },
            linear_config(),
        )
        .expect("solve")
    });
    let cubic = max_nodal_error(&mesh, &solution3.temperature, &Cubic::value);
    assert!(
        cubic < 1e-9,
        "expected nodal reproduction for a cubic solution, got {cubic:e}"
    );

    // A QUARTIC solution is NOT reproduced: that is what makes the
    // quartic ladders in tests/mms.rs non-degenerate.
    let source4 = nodal(&mesh, &|p| Quartic::source(k, p));
    let boundary4 = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "all",
            |_| true,
            ThermalBc::Dirichlet {
                temperature: nodal(&mesh, &Quartic::value),
            },
        )
        .expect("region")
        .finish()
        .expect("boundary");
    let solution4 = with_cx(|cx| {
        solve(
            cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary4,
                material: &material,
                source: &source4,
            },
            linear_config(),
        )
        .expect("solve")
    });
    let quartic = max_nodal_error(&mesh, &solution4.temperature, &Quartic::value);
    assert!(
        quartic > 1e-5,
        "a quartic solution must NOT be nodally reproduced, else every ladder is \
         degenerate; got {quartic:e}"
    );
    verdict(
        "nodal-reproduction",
        &format!("quadratic={quadratic:e} cubic={cubic:e} quartic={quartic:e}"),
    );
}

#[test]
fn newton_jacobian_matches_central_differences() {
    let n = 2;
    let (complex, positions) = unit_cube(n);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let material = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(250.0, 1.5), (450.0, 4.5)]).expect("curve"),
    );
    let source = ScalarField::Uniform(2.0e3);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "hot",
            |f| on_box_face(f.centroid[2], 0.0),
            ThermalBc::dirichlet(320.0).expect("bc"),
        )
        .expect("hot")
        .region(
            "cool",
            |f| on_box_face(f.centroid[2], 1.0),
            ThermalBc::dirichlet(300.0).expect("bc"),
        )
        .expect("cool")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let dofs = DofMap::new(&boundary, mesh.vertex_count()).expect("dofs");

    // A deterministic, non-uniform iterate so every element sees a real
    // temperature gradient (a constant iterate would zero the K'·∇T term
    // and the check would pass vacuously).
    let free: Vec<f64> = (0..dofs.n())
        .map(|i| 305.0 + 7.0 * ((i % 5) as f64) - 3.0 * ((i % 3) as f64))
        .collect();

    let residual_at = |free: &[f64]| -> Vec<f64> {
        let full = dofs.scatter(free);
        with_cx(|cx| {
            let system =
                assemble_operator(cx, &mesh, &boundary, &material, &source, &full).expect("asm");
            residual(&system, &dofs, &full)
        })
    };
    let base = residual_at(&free);
    assert!(
        base.iter().any(|v| v.abs() > 1e-6),
        "the probe iterate must produce a non-trivial residual"
    );

    let jacobian = with_cx(|cx| {
        let full = dofs.scatter(&free);
        assemble_jacobian(cx, &mesh, &boundary, &material, &full).expect("jacobian")
    });
    let (j_ff, _) = fs_conduction::assemble::reduce_matrix_and_lift(&jacobian, &dofs);

    let eps = 1e-4;
    let mut worst = 0.0f64;
    let mut scale = 0.0f64;
    for c in 0..dofs.n() {
        let mut plus = free.clone();
        let mut minus = free.clone();
        plus[c] += eps;
        minus[c] -= eps;
        let rp = residual_at(&plus);
        let rm = residual_at(&minus);
        for row in 0..dofs.n() {
            let fd = (rp[row] - rm[row]) / (2.0 * eps);
            let analytic = j_ff.get(row, c);
            scale = scale.max(fd.abs()).max(analytic.abs());
            worst = worst.max((fd - analytic).abs());
        }
    }
    assert!(
        worst <= 1e-5 * scale,
        "Jacobian deviates from central differences by {worst:e} against scale {scale:e}"
    );
    verdict(
        "jacobian-vs-central-differences",
        &format!("n={} worst={worst:e} scale={scale:e}", dofs.n()),
    );
}

#[test]
fn energy_balance_closes() {
    let fixture = enclosure_block(4);
    let solution = with_cx(|cx| solve(cx, fixture.problem(), linear_config()).expect("solve"));
    let e = solution.report.energy;
    assert!(
        e.relative_closure() < 1e-9,
        "energy balance did not close: {e:?} (relative {:e})",
        e.relative_closure()
    );
    assert!(e.source_w > 0.0, "the fixture has a positive source");
    verdict(
        "energy-balance",
        &format!(
            "source={:e} neumann_out={:e} robin_out={:e} dirichlet_in={:e} closure_rel={:e}",
            e.source_w,
            e.neumann_out_w,
            e.robin_out_w,
            e.dirichlet_in_w,
            e.relative_closure()
        ),
    );
}

#[test]
fn linear_solve_evidence_recomputes_the_euclidean_residual() {
    let fixture = enclosure_block(3);
    let solution = with_cx(|cx| solve(cx, fixture.problem(), linear_config()).expect("solve"));
    assert!(!solution.report.linear.is_empty());
    for record in &solution.report.linear {
        // CG's own claim is a RECURSIVE ESTIMATE, so `euclidean()` must
        // refuse to hand it back under the Euclidean name.
        assert_eq!(
            record.reported.euclidean(),
            None,
            "CG must not present a recurrence estimate as a Euclidean residual"
        );
        assert!(
            record.converged_true,
            "this crate gates on its own recomputed residual, which must have converged"
        );
        assert!(record.true_relative_residual.is_finite());
        println!(
            "{{\"suite\":\"fs-conduction/conformance\",\"case\":\"linear-evidence\",\
             \"method\":\"{}\",\"iters\":{},\"reported\":\"{}\",\"true_rel\":{:e}}}",
            record.method,
            record.iterations,
            record.reported.provenance(),
            record.true_relative_residual
        );
    }
}

#[test]
fn recovered_flux_matches_the_analytic_gradient() {
    // A pure linear temperature field: P₁ reproduces it exactly, so the
    // recovered element flux must equal −K∇T to round-off.
    let (complex, positions) = unit_cube(3);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let material = ConductivityModel::constant_tensor(ANISO_K).expect("material");
    let source = ScalarField::Uniform(0.0);
    let exact = |p: [f64; 3]| 300.0 + 20.0 * p[0] - 5.0 * p[1] + 2.0 * p[2];
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "all",
            |_| true,
            ThermalBc::Dirichlet {
                temperature: nodal(&mesh, &exact),
            },
        )
        .expect("region")
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
    let grad = [20.0f64, -5.0, 2.0];
    let mut want = [0.0f64; 3];
    for (i, w) in want.iter_mut().enumerate() {
        *w = -(ANISO_K[i][0] * grad[0] + ANISO_K[i][1] * grad[1] + ANISO_K[i][2] * grad[2]);
    }
    let flux = element_heat_flux(&mesh, &material, &solution.temperature).expect("flux");
    let mut worst = 0.0f64;
    for q in &flux {
        for i in 0..3 {
            worst = worst.max((q[i] - want[i]).abs());
        }
    }
    assert!(
        worst < 1e-8,
        "recovered flux deviates by {worst:e} from the analytic −K∇T = {want:?}"
    );
    verdict(
        "flux-recovery-linear-field",
        &format!("worst={worst:e} q={want:?}"),
    );
}

// ------------------------------------------------------------------- G4

#[test]
fn cancelled_assembly_is_a_structured_refusal() {
    let fixture = enclosure_block(3);
    let zero = vec![300.0f64; fixture.mesh.vertex_count()];
    let err = with_cancelled_cx(|cx| {
        assemble_operator(
            cx,
            &fixture.mesh,
            &fixture.boundary,
            &fixture.material,
            &fixture.source,
            &zero,
        )
        .expect_err("a cancelled assembly must refuse")
    });
    match err {
        ConductionError::Cancelled { stage, at } => {
            assert_eq!(stage, "assemble-elements");
            assert_eq!(at, 0);
        }
        other => panic!("expected a cancellation refusal, got {other:?}"),
    }
    assert_eq!(err.rule(), "conduction-cancelled");
    verdict("cancel-assembly", &format!("{err}"));
}

#[test]
fn cancellation_mid_iteration_drains_with_a_structured_refusal() {
    let n = 3;
    let (complex, positions) = unit_cube(n);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let material = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(250.0, 1.5), (450.0, 4.5)]).expect("curve"),
    );
    let source = ScalarField::Uniform(2.0e3);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("all", |_| true, ThermalBc::dirichlet(310.0).expect("bc"))
        .expect("region")
        .finish()
        .expect("boundary");
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        source: &source,
    };
    let config = SolveConfig {
        nonlinearity: Nonlinearity::default(),
        stop: StopRule {
            residual_rtol: 1e-12,
            residual_atol: 1e-24,
            step_atol: 0.0,
            max_iterations: 30,
        },
        linear: LinearConfig::default(),
        initial: InitialGuess::Uniform(305.0),
    };
    let (iterations_before, err) = with_gate(|gate, cx| {
        let mut solver = ConductionSolver::new(problem, config).expect("solver");
        let first = solver.step(cx).expect("first step");
        assert!(!first.converged, "the fixture must need more than one step");
        let before = solver.state().iteration;
        gate.request();
        let err = solver
            .step(cx)
            .expect_err("a cancelled iteration must refuse");
        // The state is UNCHANGED: cancellation drains, it does not
        // half-apply an update.
        assert_eq!(solver.state().iteration, before);
        (before, err)
    });
    match err {
        ConductionError::Cancelled { stage, at } => {
            assert_eq!(stage, "nonlinear-iteration");
            assert_eq!(at, iterations_before);
        }
        other => panic!("expected a cancellation refusal, got {other:?}"),
    }
    verdict("cancel-mid-iteration", &format!("{err}"));
}

// ------------------------------------------------------------------- G5

#[test]
fn replay_is_bitwise_identical() {
    let fixture = enclosure_block(4);
    let a = with_cx(|cx| solve(cx, fixture.problem(), linear_config()).expect("solve"));
    let b = with_cx(|cx| solve(cx, fixture.problem(), linear_config()).expect("solve"));
    assert_eq!(a.temperature.len(), b.temperature.len());
    for (i, (x, y)) in a.temperature.iter().zip(&b.temperature).enumerate() {
        assert_eq!(
            x.to_bits(),
            y.to_bits(),
            "vertex {i} differs between replays: {x} vs {y}"
        );
    }
    assert_eq!(a.report.iterations, b.report.iterations);
    assert_eq!(
        a.report.energy.closure_w.to_bits(),
        b.report.energy.closure_w.to_bits()
    );
    verdict(
        "replay-bitwise",
        &format!("dofs={} iters={}", a.temperature.len(), a.report.iterations),
    );
}

#[test]
fn snapshot_resume_reproduces_the_uninterrupted_trajectory() {
    let n = 3;
    let (complex, positions) = unit_cube(n);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let material = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(250.0, 1.5), (450.0, 4.5)]).expect("curve"),
    );
    let source = ScalarField::Uniform(2.0e3);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("all", |_| true, ThermalBc::dirichlet(310.0).expect("bc"))
        .expect("region")
        .finish()
        .expect("boundary");
    let problem = ConductionProblem {
        mesh: &mesh,
        boundary: &boundary,
        material: &material,
        source: &source,
    };
    let config = SolveConfig {
        nonlinearity: Nonlinearity::default(),
        stop: StopRule {
            residual_rtol: 1e-12,
            residual_atol: 1e-24,
            step_atol: 0.0,
            max_iterations: 30,
        },
        linear: LinearConfig::default(),
        initial: InitialGuess::Uniform(305.0),
    };

    let straight = with_cx(|cx| {
        ConductionSolver::new(problem, config.clone())
            .expect("solver")
            .run(cx)
    })
    .expect("straight run");

    let (sealed, split) = with_cx(|cx| {
        let mut solver = ConductionSolver::new(problem, config.clone()).expect("solver");
        solver.step(cx).expect("step 1");
        let sealed = solver.snapshot(0xABCD_EF01);
        let mut resumed = ConductionSolver::new(problem, config.clone()).expect("solver");
        let provenance = resumed.restore(&sealed).expect("restore");
        assert_eq!(provenance, 0xABCD_EF01);
        (sealed, resumed.run(cx).expect("resumed run"))
    });

    for (i, (x, y)) in straight
        .temperature
        .iter()
        .zip(&split.temperature)
        .enumerate()
    {
        assert_eq!(
            x.to_bits(),
            y.to_bits(),
            "vertex {i}: split run {y} != straight run {x}"
        );
    }
    assert_eq!(straight.report.iterations, split.report.iterations);
    verdict(
        "snapshot-resume-bitwise",
        &format!(
            "envelope_bytes={} iters={}",
            sealed.len(),
            straight.report.iterations
        ),
    );
}

#[test]
fn snapshot_envelope_refuses_tampered_bytes() {
    let state = ConductionState {
        iteration: 3,
        free_temperature: vec![301.5, 302.25, 299.125],
        residual_history: vec![1.0e3, 4.0e1],
        last_step_norm: 0.5,
    };
    let sealed = state.seal(0x1234);
    let (round, provenance) = ConductionState::unseal(&sealed).expect("round trip");
    assert_eq!(round, state);
    assert_eq!(provenance, 0x1234);

    // A flipped payload byte must be refused by the envelope checksum,
    // never decoded into plausible-but-wrong state.
    let mut tampered = sealed.clone();
    let last = tampered.len() - 1;
    tampered[last] ^= 0x01;
    assert!(ConductionState::unseal(&tampered).is_err());

    // Truncation and appended garbage are refused too.
    assert!(ConductionState::unseal(&sealed[..sealed.len() - 4]).is_err());
    let mut appended = sealed.clone();
    appended.push(0x00);
    assert!(ConductionState::unseal(&appended).is_err());
    verdict(
        "snapshot-envelope",
        &format!(
            "type_id=0x{:016x} bytes={}",
            ConductionState::TYPE_ID,
            sealed.len()
        ),
    );
}

// -------------------------------------------------- material provenance

fn conductivity_claims() -> ClaimSet {
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "AA6061-T6, longitudinal, as-extruded".to_string(),
            method: "ASTM E1461 laser flash / frame A".to_string(),
            artifact: hash_bytes(b"fs-conduction-conformance-thermal-conductivity-table"),
            caveats: "conductivity measured jointly with diffusivity; no censoring".to_string(),
            provenance: Provenance {
                source: "fs-conduction conformance fixture".to_string(),
                license: "CC-BY-4.0".to_string(),
                artifact: None,
            },
        })
        .expect("observation");
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("thermal_conductivity", CONDUCTIVITY_DIMS),
            value: PropertyValue::Curve {
                abscissa: "T".to_string(),
                abscissa_dims: TEMPERATURE_DIMS,
                knots: vec![
                    (250.0, 160.0),
                    (300.0, 167.0),
                    (400.0, 177.0),
                    (500.0, 186.0),
                ],
                dims: CONDUCTIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 250.0, 500.0),
            uncertainty: UncertaintyModel::RelativeHalfWidth {
                fraction: 0.03,
                confidence: 0.95,
            },
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![observation],
            provenance: Provenance {
                source: "fs-conduction conformance fixture".to_string(),
                license: "CC-BY-4.0".to_string(),
                artifact: None,
            },
        })
        .expect("claim");
    claims
}

#[test]
fn matdb_receipts_travel_with_the_solve() {
    let claims = conductivity_claims();
    let grid = [280.0f64, 300.0, 320.0, 340.0, 360.0];
    let table = ConductivityTable::from_claims(
        &claims,
        "thermal_conductivity",
        &grid,
        SelectionPolicy::PreferObservationBacked,
    )
    .expect("table");
    assert_eq!(table.receipts().len(), grid.len());
    assert_eq!(table.provenance(), ProvenanceClass::MatdbReceipts);
    // Every retained receipt must still verify against the claim set it
    // came from — a receipt that cannot be replayed is not provenance.
    for receipt in table.receipts() {
        claims.verify_receipt(receipt).expect("receipt replays");
        assert!(receipt.observation_backed);
    }

    let material = ConductivityModel::isotropic(table);
    let (complex, positions) = unit_cube(3);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let source = ScalarField::Uniform(1.0e5);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("all", |_| true, ThermalBc::dirichlet(300.0).expect("bc"))
        .expect("region")
        .finish()
        .expect("boundary");
    let config = SolveConfig {
        nonlinearity: Nonlinearity::default(),
        stop: StopRule {
            residual_rtol: 1e-11,
            residual_atol: 1e-24,
            step_atol: 0.0,
            max_iterations: 30,
        },
        linear: LinearConfig::default(),
        initial: InitialGuess::Uniform(300.0),
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
        .expect("solve")
    });
    assert_eq!(solution.report.material_receipts, grid.len());
    assert_eq!(
        solution.report.material_provenance,
        ProvenanceClass::MatdbReceipts
    );
    verdict(
        "matdb-receipts",
        &format!(
            "receipts={} provenance={} iterations={}",
            solution.report.material_receipts,
            solution.report.material_provenance.tag(),
            solution.report.iterations
        ),
    );
}

#[test]
fn declared_conductivity_never_claims_matdb_provenance() {
    let material = ConductivityModel::isotropic_declared(200.0).expect("material");
    assert_eq!(material.provenance(), ProvenanceClass::Declared);
    assert!(material.receipts().is_empty());
    let tensor = ConductivityModel::constant_tensor(ANISO_K).expect("tensor");
    assert_eq!(tensor.provenance(), ProvenanceClass::Declared);
    verdict("declared-provenance", "no receipts, no matdb claim");
}

#[test]
fn material_queries_outside_validity_refuse() {
    let claims = conductivity_claims();
    // 600 K is outside the claim's validity box.
    let err = ConductivityTable::from_claims(
        &claims,
        "thermal_conductivity",
        &[500.0, 600.0],
        SelectionPolicy::PreferObservationBacked,
    )
    .expect_err("out-of-validity sampling must refuse");
    assert_eq!(err.rule(), "conduction-material-query");

    // A table refuses to extrapolate beyond the grid it was sampled on.
    let table = ConductivityTable::from_claims(
        &claims,
        "thermal_conductivity",
        &[300.0, 400.0],
        SelectionPolicy::PreferObservationBacked,
    )
    .expect("table");
    assert!(table.eval(299.0).is_err());
    assert!(table.eval(401.0).is_err());
    let mid = table.eval(350.0).expect("inside the span");
    assert!(mid > 167.0 && mid < 177.0);
    verdict("material-refusals", &format!("{err}"));
}

#[test]
fn non_spd_conductivity_tensors_refuse() {
    // Asymmetric.
    assert!(
        ConductivityModel::constant_tensor([[3.0, 0.5, 0.0], [0.9, 2.0, 0.0], [0.0, 0.0, 1.0]])
            .is_err()
    );
    // Symmetric but indefinite.
    assert!(
        ConductivityModel::constant_tensor([[1.0, 2.0, 0.0], [2.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
            .is_err()
    );
    // Negative conductivity.
    assert!(ConductivityModel::isotropic_declared(-1.0).is_err());
    verdict(
        "conductivity-admission",
        "asymmetric/indefinite/negative refused",
    );
}

// ------------------------------------------------------- scenario seams

#[test]
fn scenario_thermal_rows_lower_into_conditions() {
    let env = Environment::earth_lab();

    let dirichlet = BoundaryCondition {
        region: "base".to_string(),
        physics: Physics::Thermal,
        kind: BcKind::Dirichlet,
        value: Some(BcValue::Uniform(QtyAny::new(350.0, TEMPERATURE_DIMS))),
        compatibility: None,
        frame: 0,
    };
    assert!(matches!(
        ThermalBc::from_scenario_row(&dirichlet, &env, None).expect("dirichlet"),
        ThermalBc::Dirichlet { .. }
    ));

    let neumann = BoundaryCondition {
        region: "die".to_string(),
        physics: Physics::Thermal,
        kind: BcKind::Neumann,
        value: Some(BcValue::Uniform(QtyAny::new(
            -2.5e4,
            fs_conduction::HEAT_FLUX_DIMS,
        ))),
        compatibility: None,
        frame: 0,
    };
    assert!(matches!(
        ThermalBc::from_scenario_row(&neumann, &env, None).expect("neumann"),
        ThermalBc::Neumann { .. }
    ));

    // THE ROBIN SEAM: fs-scenario's row carries h ONLY, so T_ref must be
    // named at the lowering call. Absent an explicit value the scenario
    // ENVIRONMENT's ambient temperature is used — a declaration, not a
    // hidden default.
    let robin = BoundaryCondition {
        region: "fins".to_string(),
        physics: Physics::Thermal,
        kind: BcKind::Robin,
        value: Some(BcValue::Uniform(QtyAny::new(45.0, fs_conduction::HTC_DIMS))),
        compatibility: None,
        frame: 0,
    };
    let from_env = ThermalBc::from_scenario_row(&robin, &env, None).expect("robin");
    match &from_env {
        ThermalBc::Robin { htc, t_ref } => {
            assert_eq!(htc.at(0).to_bits(), 45.0f64.to_bits());
            assert_eq!(
                t_ref.at(0).to_bits(),
                env.ambient_temperature.value.to_bits()
            );
        }
        other => panic!("expected a Robin row, got {other:?}"),
    }
    let explicit = ThermalBc::from_scenario_row(&robin, &env, Some(313.15)).expect("robin");
    match &explicit {
        ThermalBc::Robin { t_ref, .. } => {
            assert_eq!(t_ref.at(0).to_bits(), 313.15f64.to_bits());
        }
        other => panic!("expected a Robin row, got {other:?}"),
    }

    // Refusals: wrong physics, wrong dimensions, and a time signal
    // (this crate is steady-only).
    let wrong_physics = BoundaryCondition {
        physics: Physics::Elasticity,
        ..dirichlet.clone()
    };
    assert_eq!(
        ThermalBc::from_scenario_row(&wrong_physics, &env, None)
            .expect_err("non-thermal")
            .rule(),
        "conduction-scenario-row"
    );
    let wrong_dims = BoundaryCondition {
        value: Some(BcValue::Uniform(QtyAny::new(
            350.0,
            Dims([1, 0, 0, 0, 0, 0]),
        ))),
        ..dirichlet.clone()
    };
    assert_eq!(
        ThermalBc::from_scenario_row(&wrong_dims, &env, None)
            .expect_err("wrong dims")
            .rule(),
        "conduction-dimensions"
    );
    let no_value = BoundaryCondition {
        value: None,
        ..dirichlet.clone()
    };
    assert_eq!(
        ThermalBc::from_scenario_row(&no_value, &env, None)
            .expect_err("no value")
            .rule(),
        "conduction-scenario-row"
    );
    verdict(
        "scenario-lowering",
        "dirichlet/neumann/robin lowered; T_ref named at the seam; four refusals typed",
    );
}

// -------------------------------------------------------- refusal suite

#[test]
fn boundary_partition_refusals_are_typed() {
    let (complex, positions) = unit_cube(2);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");

    // Gap: leftover faces without an explicit adiabatic remainder.
    let gap = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "z-lo",
            |f| on_box_face(f.centroid[2], 0.0),
            ThermalBc::dirichlet(300.0).expect("bc"),
        )
        .expect("region")
        .finish()
        .expect_err("untagged faces must refuse");
    assert_eq!(gap.rule(), "conduction-untagged-boundary");

    // Overlap: a second region claiming a face the first owns.
    let overlap = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "z-lo",
            |f| on_box_face(f.centroid[2], 0.0),
            ThermalBc::dirichlet(300.0).expect("bc"),
        )
        .expect("region")
        .region("everything", |_| true, ThermalBc::adiabatic())
        .err()
        .expect("overlap must refuse");
    assert_eq!(overlap.rule(), "conduction-overlapping-region");

    // Duplicate name.
    let duplicate = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "z-lo",
            |f| on_box_face(f.centroid[2], 0.0),
            ThermalBc::dirichlet(300.0).expect("bc"),
        )
        .expect("region")
        .region(
            "z-lo",
            |f| on_box_face(f.centroid[2], 1.0),
            ThermalBc::dirichlet(310.0).expect("bc"),
        )
        .err()
        .expect("duplicate name must refuse");
    assert_eq!(duplicate.rule(), "conduction-duplicate-region");

    // A non-positive transfer coefficient is not a Robin row.
    assert!(ThermalBc::robin(0.0, 300.0).is_err());
    verdict(
        "boundary-refusals",
        &format!("{gap} | {overlap} | {duplicate}"),
    );
}

#[test]
fn pure_neumann_problem_refuses_instead_of_returning_a_plausible_field() {
    let (complex, positions) = unit_cube(2);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let material = ConductivityModel::isotropic_declared(ISO_K).expect("material");
    let source = ScalarField::Uniform(0.0);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("all", |_| true, ThermalBc::adiabatic())
        .expect("region")
        .finish()
        .expect("boundary");
    let err = with_cx(|cx| {
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
        .expect_err("a pure-Neumann steady problem is singular")
    });
    assert_eq!(err.rule(), "conduction-singular-pure-neumann");
    verdict("pure-neumann", &format!("{err}"));
}

#[test]
fn degenerate_and_mismatched_meshes_refuse() {
    // A flat tet: zero volume.
    let complex = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let flat = vec![
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.5, 0.5, 0.0],
    ];
    let err = ConductionMesh::new(complex, flat).expect_err("degenerate tet must refuse");
    assert_eq!(err.rule(), "conduction-degenerate-element");

    // Wrong number of positions.
    let complex = TetComplex::from_tets(4, vec![[0, 1, 2, 3]]);
    let short = vec![[0.0, 0.0, 0.0]];
    let err2 = ConductionMesh::new(complex, short).expect_err("size mismatch must refuse");
    assert_eq!(err2.rule(), "conduction-mesh");
    verdict("mesh-refusals", &format!("{err} | {err2}"));
}

#[test]
fn adjoint_hook_refuses_a_temperature_dependent_material() {
    let (complex, positions) = unit_cube(2);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let material = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(250.0, 1.5), (450.0, 4.5)]).expect("curve"),
    );
    let source = ScalarField::Uniform(1.0e4);
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("all", |_| true, ThermalBc::dirichlet(300.0).expect("bc"))
        .expect("region")
        .finish()
        .expect("boundary");
    let err = ConductivityDesign::new(
        ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            source: &source,
        },
        LinearConfig::default(),
    )
    .err()
    .expect("the IFT hook covers the linear case only");
    assert_eq!(err.rule(), "conduction-conductivity");
    verdict("adjoint-scope", &format!("{err}"));
}

#[test]
fn element_stiffness_rows_sum_to_zero() {
    // A constant temperature field carries no flux, so every element
    // matrix must annihilate the constant vector — the discrete form of
    // `∇·(K∇c) = 0`. This is what makes the conduction block singular
    // before boundary conditions, and why a pure-Neumann problem refuses.
    let (complex, positions) = unit_cube(2);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let mut worst = 0.0f64;
    let mut scale = 0.0f64;
    for e in 0..mesh.element_count() {
        let ke = element_stiffness(&mesh, e, &ANISO_K);
        for row in &ke {
            let sum: f64 = row.iter().sum();
            worst = worst.max(sum.abs());
            scale = scale.max(row.iter().fold(0.0f64, |a, b| a.max(b.abs())));
        }
    }
    assert!(
        worst <= 1e-12 * scale,
        "element row sums {worst:e} exceed 1e-12 x {scale:e}"
    );
    verdict(
        "element-nullspace",
        &format!("worst_row_sum={worst:e} scale={scale:e}"),
    );
}
