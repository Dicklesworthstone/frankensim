//! Real thermal assembly, a separate dense elimination oracle, and refusals
//! for the stored-system goal analyzer. No continuum bound is inferred.

mod support;
#[path = "algebraic_goal/controlled.rs"]
mod controlled;

use fs_conduction::adjoint::{
    LinearGoalAnalysisConfig, LinearGoalAnalyzer, analyze_linear_goal,
};
use fs_conduction::assemble::{DofMap, assemble_operator, reduce};
use fs_conduction::bc::{ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
use fs_conduction::field::ScalarField;
use fs_conduction::fixtures::{on_box_face, unit_cube};
use fs_conduction::material::{ConductivityModel, ConductivityTable};
use fs_conduction::mesh::ConductionMesh;
use fs_conduction::solve::{ConductionProblem, LinearConfig};
use fs_conduction::ConductionError;
use fs_solver::goal::{GoalBoundStatus, GoalResidualLimits};
use fs_sparse::Csr;
use support::{with_cancelled_cx, with_cx};

struct Fixture {
    mesh: ConductionMesh,
    boundary: ThermalBoundary,
    material: ConductivityModel,
    source: ScalarField,
}

impl Fixture {
    fn new(all_fixed_faces: bool, robin: bool) -> Self {
        let (complex, positions) = unit_cube(2);
        let mesh = ConductionMesh::new(complex, positions).unwrap();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("cold", |f| all_fixed_faces || on_box_face(f.centroid[0], 0.0),
                ThermalBc::dirichlet(300.0).unwrap()).unwrap();
        let boundary = if robin {
            boundary.region("exchange", |f| on_box_face(f.centroid[0], 1.0),
                ThermalBc::robin(2.0, 295.0).unwrap()).unwrap()
        } else { boundary };
        let boundary = boundary.adiabatic_remainder().finish().unwrap();
        Self { mesh, boundary, material: ConductivityModel::isotropic_declared(2.0).unwrap(),
            source: ScalarField::Uniform(6.0) }
    }

    fn problem(&self) -> ConductionProblem<'_> {
        ConductionProblem { mesh: &self.mesh, boundary: &self.boundary,
            material: &self.material, source: &self.source, element_materials: None }
    }

    fn initial(&self) -> Vec<f64> { vec![300.0; self.mesh.vertex_count()] }

    fn weights(&self) -> Vec<f64> {
        let dofs = DofMap::new(&self.boundary, self.mesh.vertex_count()).unwrap();
        let mut weights = vec![0.0; self.mesh.vertex_count()];
        for &v in dofs.free() { weights[v] = 1.0 / dofs.n() as f64; }
        weights
    }
}

fn linear() -> LinearConfig {
    LinearConfig { tolerance: 1e-11, max_iterations: 1024, restart: 40 }
}
fn config() -> LinearGoalAnalysisConfig {
    LinearGoalAnalysisConfig {
        residual_limits: GoalResidualLimits { max_rows: 1024, max_nonzeros: 65536 },
        max_stability_iterations: 1024,
    }
}

// Partial-pivot Gaussian elimination, independent of every production Krylov
// recurrence, preconditioner, inverse-dominance test and interval evaluator.
fn dense_solve(matrix: &Csr, rhs: &[f64]) -> Vec<f64> {
    let n = rhs.len();
    let mut rows: Vec<Vec<f64>> = (0..n).map(|i| {
        let mut row: Vec<_> = (0..n).map(|j| matrix.get(i, j)).collect();
        row.push(rhs[i]); row
    }).collect();
    for k in 0..n {
        let pivot = (k..n).max_by(|&i, &j| rows[i][k].abs().total_cmp(&rows[j][k].abs())).unwrap();
        rows.swap(k, pivot);
        assert!(rows[k][k].abs() > 1e-14);
        let pivot_row = rows[k].clone();
        for row in rows.iter_mut().skip(k+1) {
            let ratio = row[k] / pivot_row[k];
            for j in k..=n { row[j] -= ratio * pivot_row[j]; }
        }
    }
    let mut x = vec![0.0; n];
    for i in (0..n).rev() {
        let tail: f64 = (i+1..n).map(|j| rows[i][j] * x[j]).sum();
        x[i] = (rows[i][n] - tail) / rows[i][i];
    }
    x
}

fn oracle(f: &Fixture, cx: &fs_exec::Cx<'_>) -> Vec<f64> {
    let initial = f.initial();
    let system = assemble_operator(cx, &f.mesh, &f.boundary, &f.material, &f.source, &initial).unwrap();
    let dofs = DofMap::new(&f.boundary, f.mesh.vertex_count()).unwrap();
    let (matrix, rhs) = reduce(&system, &dofs);
    dofs.scatter(&dense_solve(&matrix, &rhs))
}

#[test]
fn scalar_thermal_goal_encloses_true_error_and_detects_residual_sign() {
    let f = Fixture::new(true, false);
    let initial = f.initial();
    let weights = f.weights();
    with_cx(|cx| {
        let analyzer = LinearGoalAnalyzer::new(cx, f.problem(), None, linear(), &initial, &weights, config()).unwrap();
        assert_eq!(analyzer.dofs().n(), 1);
        let report = analyzer.analyze(cx, &initial).unwrap();
        let exact = oracle(&f, cx);
        let v = analyzer.dofs().free()[0];
        let error = exact[v] - initial[v];
        assert!(error > 0.01);
        let bound = report.enclosure.goal_error().unwrap();
        assert!(bound.lower() <= error && error <= bound.upper(), "{error} {bound:?}");
        assert!(bound.lower() > 0.0, "residual sign reversal must be detected");
        assert!(bound.upper() - bound.lower() < 1e-8);
        assert!(report.enclosure.evaluation_roundoff_upper() > 0.0);
        assert!(report.enclosure.dual_error_upper().unwrap() < 1e-9);
        assert_eq!(report.stability_iterations, 0, "one scalar is already strictly dominant");
        assert!(!report.meets_absolute_tolerance(0.001));
        assert!(report.meets_absolute_tolerance(1.0));
        for bad in [0.0, -1.0, f64::INFINITY, f64::NAN] {
            assert!(!report.meets_absolute_tolerance(bad));
        }
    });
}

#[test]
fn cached_dual_and_checked_scaling_assess_multiple_real_thermal_iterates() {
    let f = Fixture::new(false, false);
    let initial = f.initial();
    let weights = f.weights();
    with_cx(|cx| {
        let mut no_proposal = config(); no_proposal.max_stability_iterations = 0;
        let missing = analyze_linear_goal(cx, f.problem(), None, linear(), &initial, &weights, no_proposal).unwrap();
        assert_eq!(missing.enclosure.status(), GoalBoundStatus::InverseBoundUnavailable);
        assert!(!missing.meets_absolute_tolerance(1e10), "NO-DATA is not zero error");
        let analyzer = LinearGoalAnalyzer::new(cx, f.problem(), None, linear(), &initial, &weights, config()).unwrap();
        let first = analyzer.analyze(cx, &initial).unwrap();
        assert_eq!(first.enclosure.status(), GoalBoundStatus::Enclosed);
        assert!(first.stability_iterations > 0);
        assert!(analyzer.stability_scaling().unwrap().iter().all(|w| *w > 0.0));
        let exact = oracle(&f, cx);
        let dual = analyzer.free_dual().to_vec();
        for fraction in [0.0, 0.25, 0.75] {
            let field: Vec<_> = initial.iter().zip(&exact).map(|(a,b)| a + fraction*(b-a)).collect();
            let report = analyzer.analyze(cx, &field).unwrap();
            let error: f64 = weights.iter().zip(exact.iter().zip(&field)).map(|(w,(a,b))| w*(a-b)).sum();
            let bound = report.enclosure.goal_error().unwrap();
            assert!(bound.lower() <= error && error <= bound.upper(), "{error} {bound:?}");
            assert_eq!(report.dual_iterations, first.dual_iterations);
            assert_eq!(report.stability_iterations, first.stability_iterations);
            assert_eq!(analyzer.free_dual(), dual);
        }
        assert_eq!(first, analyzer.analyze(cx, &initial).unwrap());
    });
}

#[test]
fn fixed_robin_operator_is_retained_in_the_goal_correction() {
    let f = Fixture::new(false, true);
    with_cx(|cx| {
        let initial = f.initial();
        let weights = f.weights();
        let report = analyze_linear_goal(cx, f.problem(), None, linear(), &initial, &weights, config()).unwrap();
        let exact = oracle(&f, cx);
        let expected: f64 = weights.iter().zip(exact.iter().zip(&initial)).map(|(w,(a,b))| w*(a-b)).sum();
        assert!(expected.abs() > 0.1);
        assert!((report.enclosure.nominal_correction() - expected).abs() < 1e-8);
        assert!(report.dual_relative_residual < linear().tolerance);
        // Consistent Robin mass can destroy M-matrix dominance. Its absence
        // is NOT repaired by pretending a diagonal lumped matrix was solved.
        if let Some(bound) = report.enclosure.goal_error() {
            assert!(bound.lower() <= expected && expected <= bound.upper());
        } else {
            assert!(!report.meets_absolute_tolerance(1e10));
        }
    });
}

#[test]
fn prescribed_goal_weights_cancel_before_adjoint_scaling() {
    let f = Fixture::new(true, false);
    let initial = f.initial();
    let weights = f.weights();
    with_cx(|cx| {
        let baseline = analyze_linear_goal(cx, f.problem(), None, linear(), &initial, &weights, config()).unwrap();
        let dofs = DofMap::new(&f.boundary, f.mesh.vertex_count()).unwrap();
        let mut changed = weights.clone();
        for &v in dofs.fixed() { changed[v] = 1e200; }
        let with_constant = analyze_linear_goal(cx, f.problem(), None, linear(), &initial, &changed, config()).unwrap();
        assert_eq!(baseline, with_constant);
        for factor in [1e-100, -2.0, 1e100] {
            let scaled: Vec<_> = weights.iter().map(|w| w*factor).collect();
            let report = analyze_linear_goal(cx, f.problem(), None, linear(), &initial, &scaled, config()).unwrap();
            assert!((report.enclosure.nominal_correction()/factor - baseline.enclosure.nominal_correction()).abs() < 1e-10);
            assert!(report.enclosure.goal_error().is_some());
        }
    });
}

#[test]
fn nonlinear_conductivity_refuses_and_sampled_support_is_rechecked() {
    let mut f = Fixture::new(true, false);
    f.material = ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![(250.0, 1.0), (400.0, 2.0)]).unwrap());
    with_cx(|cx| {
        assert!(matches!(analyze_linear_goal(cx, f.problem(), None, linear(), &f.initial(), &f.weights(), config()), Err(ConductionError::Config { .. })));
    });
    f.material = ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![(290.0, 2.0), (310.0, 2.0)]).unwrap());
    with_cx(|cx| {
        let mut initial = f.initial();
        let analyzer = LinearGoalAnalyzer::new(cx, f.problem(), None, linear(), &initial, &f.weights(), config()).unwrap();
        let free = analyzer.dofs().free()[0];
        initial[free] = 500.0; // adjacent element mean is 350 K, outside support
        assert!(matches!(analyzer.analyze(cx, &initial), Err(ConductionError::OutsideTemperatureSpan { .. })));
    });
}

#[test]
fn shape_prescribed_temperature_and_structural_budgets_refuse() {
    let f = Fixture::new(true, false);
    with_cx(|cx| {
        let initial = f.initial();
        let weights = f.weights();
        let analyzer = LinearGoalAnalyzer::new(cx, f.problem(), None, linear(), &initial, &weights, config()).unwrap();
        assert!(matches!(analyzer.analyze(cx, &[]), Err(ConductionError::FieldLength { .. })));
        let mut changed = initial.clone(); changed[analyzer.dofs().fixed()[0]] += 1.0;
        assert!(matches!(analyzer.analyze(cx, &changed), Err(ConductionError::Config { .. })));
        let mut nonfinite = initial.clone(); nonfinite[analyzer.dofs().free()[0]] = f64::NAN;
        assert!(matches!(analyzer.analyze(cx, &nonfinite), Err(ConductionError::NonFinite { .. })));
        for limits in [GoalResidualLimits { max_rows: 0, max_nonzeros: 100 }, GoalResidualLimits { max_rows: 100, max_nonzeros: 0 }] {
            let mut bounded = config(); bounded.residual_limits = limits;
            assert!(matches!(analyze_linear_goal(cx, f.problem(), None, linear(), &initial, &weights, bounded), Err(ConductionError::Config { .. })));
        }
    });
}

#[test]
fn cancellation_never_publishes_partial_authority_and_retry_is_identical() {
    let f = Fixture::new(false, false);
    let initial = f.initial();
    let weights = f.weights();
    with_cancelled_cx(|cx| {
        assert!(matches!(analyze_linear_goal(cx, f.problem(), None, linear(), &initial, &weights, config()), Err(ConductionError::Cancelled { .. })));
    });
    let analyzer = with_cx(|cx| LinearGoalAnalyzer::new(cx, f.problem(), None, linear(), &initial, &weights, config()).unwrap());
    let first = with_cx(|cx| analyzer.analyze(cx, &initial).unwrap());
    with_cancelled_cx(|cx| assert!(matches!(analyzer.analyze(cx, &initial), Err(ConductionError::Cancelled { .. }))));
    let retried = with_cx(|cx| analyzer.analyze(cx, &initial).unwrap());
    assert_eq!(first, retried);
    assert!(initial.iter().all(|&v| v == 300.0));
}

fn contact_model<R>(resistance_value: f64,
    f: impl FnOnce(ConductionProblem<'_>, &fs_conduction::ThermalInterfaces) -> R) -> R {
    use fs_conduction::{InterfaceFacePair, InterfaceResistance, InterfaceSurface, ThermalInterfaces,
        AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP};
    use fs_conduction::fixtures::box_grid;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId,
        PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy,
        SurfaceSpec, SystemContext, UncertaintyModel};
    use fs_rep_mesh::TetComplex;
    let mut positions = Vec::new();
    let mut tets = Vec::new();
    for side in 0..2 {
        let (complex, points) = box_grid([2, 1, 1], [1.0, 1.0, 1.0]);
        let offset = positions.len() as u32;
        tets.extend(complex.tets.into_iter().map(|tet| tet.map(|v| v + offset)));
        positions.extend(points.into_iter().map(|[x, y, z]| [x + f64::from(side), y, z]));
    }
    let mesh = ConductionMesh::new(TetComplex::from_tets(positions.len(), tets), positions).unwrap();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("hot", |face| on_box_face(face.centroid[0], 0.0), ThermalBc::robin(20.0, 330.0).unwrap()).unwrap()
        .region("cold", |face| on_box_face(face.centroid[0], 2.0), ThermalBc::robin(40.0, 300.0).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap();
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim { key: PropertyKey::new(RP, RD),
        value: PropertyValue::Scalar { value: resistance_value, dims: RD }, validity: ValidityDomain::unconstrained(),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: Vec::new(), provenance: Provenance { source: "declared algebraic contact fixture".into(),
            license: "internal-test-use".into(), artifact: None } }).unwrap();
    let side = |name: &str| SurfaceSpec { material: MaterialStateId { chemistry: name.into(),
        phase: "solid".into(), process: "as-fixtured".into(), revision: 0 }, texture_frame: "declared".into() };
    let card = InterfaceSystemCard::assemble(side("a"), side("b"), SystemContext {
        medium: "declared".into(), third_body: None, environment: "declared".into(), history: "declared".into(),
    }, claims, Vec::new()).unwrap();
    let resistance = InterfaceResistance::from_card("bond", &card, &QueryPoint::new(), SelectionPolicy::SingleClaimOnly).unwrap();
    let pairs = ThermalInterfaces::coincident_face_pairs(&mesh).unwrap().into_iter().map(|pair| {
        if mesh.boundary()[pair.side_a].outward_normal[0] > 0.0 { pair }
        else { InterfaceFacePair { side_a: pair.side_b, side_b: pair.side_a } }
    }).collect();
    let interfaces = ThermalInterfaces::new(&mesh, &boundary,
        vec![InterfaceSurface::new("bond", pairs, resistance).unwrap()]).unwrap();
    let material = ConductivityModel::isotropic_declared(10.0).unwrap();
    let source = ScalarField::Uniform(0.0);
    f(ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
        element_materials: None, source: &source }, &interfaces)
}

#[test]
fn contact_resistance_reaches_the_goal_and_missing_pairs_refuse() {
    with_cx(|cx| {
        for resistance in [0.1, 0.2] {
            contact_model(resistance, |problem, interfaces| {
                let initial = vec![315.0; problem.mesh.vertex_count()];
                let mut weights = vec![0.0; initial.len()];
                let cold = problem.mesh.positions().iter().position(|p| on_box_face(p[0], 2.0)).unwrap();
                weights[cold] = 1.0;
                let report = analyze_linear_goal(cx, problem, Some(interfaces), linear(), &initial, &weights, config()).unwrap();
                // Independent one-dimensional series resistance oracle: two
                // 1 m / 10 W/(m K) slabs, contact, and both Robin films.
                let heat = 30.0 / (1.0/20.0 + 0.1 + resistance + 0.1 + 1.0/40.0);
                let expected = 300.0 + heat/40.0 - initial[cold];
                assert!((report.enclosure.nominal_correction() - expected).abs() < 1e-8);
                assert!(report.enclosure.weighted_residual().upper() < -10.0);
                assert!(analyze_linear_goal(cx, problem, None, linear(), &initial, &weights, config()).is_err());
            });
        }
    });
}
