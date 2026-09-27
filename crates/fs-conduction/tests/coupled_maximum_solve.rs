//! Coupled corrections through real FEM assembly, with an independent dense
//! oracle. The frozen-solid system is deliberately not the target equation.
mod support;

use fs_conduction::adjoint::{LinearGoalAnalysisConfig, LinearGoalAnalyzer,
    LinearGoalSolveConfig, LinearGoalStop, RobinFeedbackAnalysisConfig};
use fs_conduction::assemble::DofMap;
use fs_conduction::fixtures::{on_box_face, unit_cube};
use fs_conduction::{ConductionMesh, ConductionProblem, ConductivityModel, LinearConfig,
    ScalarField, ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
use fs_solver::goal::GoalResidualLimits;
use fs_solver::goal::feedback::{FeedbackInverseMethod, FeedbackResidualLimits};
use fs_sparse::Csr;
use support::{with_cx, with_gate};

struct Fixture {
    mesh: ConductionMesh, boundary: ThermalBoundary,
    material: ConductivityModel, source: ScalarField,
}
impl Fixture {
    fn new(fixed: bool, split: bool) -> Self {
        let (complex, positions) = unit_cube(1);
        let mesh = ConductionMesh::new(complex, positions).unwrap();
        let mut builder = ThermalBoundaryBuilder::new(&mesh);
        if fixed {
            builder = builder.region("cold", |f| on_box_face(f.centroid[0], 0.0),
                ThermalBc::dirichlet(300.0).unwrap()).unwrap();
        }
        let boundary = if split {
            builder.region("first", |f| f.centroid[0] < 0.5,
                ThermalBc::robin(2.0, 300.0).unwrap()).unwrap()
                .region("second", |f| f.centroid[0] >= 0.5,
                    ThermalBc::robin(4.0, 300.0).unwrap()).unwrap()
                .finish().unwrap()
        } else {
            builder.region("air", |f| !fixed || !on_box_face(f.centroid[0], 0.0),
                ThermalBc::robin(2.0, 300.0).unwrap()).unwrap().finish().unwrap()
        };
        Self { mesh, boundary, material: ConductivityModel::isotropic_declared(10.0).unwrap(),
            source: ScalarField::Uniform(0.0) }
    }
    fn base(&self, cx: &fs_exec::Cx<'_>) -> LinearGoalAnalyzer<'_> {
        LinearGoalAnalyzer::new_for_maximum(cx, ConductionProblem {
            mesh: &self.mesh, boundary: &self.boundary, material: &self.material,
            element_materials: None, source: &self.source,
        }, None, LinearConfig { tolerance: 1e-12, max_iterations: 1000, restart: 20 },
            &vec![300.0; self.mesh.vertex_count()], LinearGoalAnalysisConfig {
                residual_limits: feedback_config().residual.solid,
                max_stability_iterations: 1000,
            }).unwrap()
    }
    fn vertices(&self) -> Vec<usize> { (0..self.mesh.vertex_count()).collect() }
}
fn feedback_config() -> RobinFeedbackAnalysisConfig {
    RobinFeedbackAnalysisConfig {
        residual: FeedbackResidualLimits {
            solid: GoalResidualLimits { max_rows: 100, max_nonzeros: 100_000 },
            max_ports: 4, max_transfer_nonzeros: 10_000,
            max_response_entries: 400, max_verification_entries: 10_000_000,
        }, max_response_iterations: 1000, max_lowering_entries: 100_000,
    }
}
fn controls() -> LinearGoalSolveConfig {
    LinearGoalSolveConfig { absolute_tolerance: 1e-6, max_primal_iterations: 160,
        check_every: 8, max_defect_corrections: 2 }
}
fn dense(a: &Csr, rhs: &[f64], b: &Csr, c: &Csr, d: &[f64]) -> Vec<f64> {
    let n = rhs.len();
    let mut rows = vec![vec![0.0; n+1]; n];
    for i in 0..n {
        rows[i][n] = rhs[i];
        for k in 0..d.len() { rows[i][n] += b.get(i,k)*d[k]; }
        for j in 0..n {
            rows[i][j] = a.get(i,j);
            for k in 0..d.len() { rows[i][j] -= b.get(i,k)*c.get(k,j); }
        }
    }
    for k in 0..n {
        let pivot = (k..n).max_by(|&i,&j| rows[i][k].abs().total_cmp(&rows[j][k].abs())).unwrap();
        rows.swap(k,pivot);
        let diagonal = rows[k][k]; assert!(diagonal.abs() > 1e-12);
        for j in k..=n { rows[k][j] /= diagonal; }
        for i in 0..n { if i != k {
            let factor = rows[i][k];
            for j in k..=n { let value = rows[k][j]; rows[i][j] -= factor*value; }
        } }
    }
    rows.iter().map(|r| r[n]).collect()
}

#[test]
fn coupled_corrections_solve_even_a_noncontractive_reference_law() {
    let fixture = Fixture::new(false, false);
    with_cx(|cx| {
        for slope in [0.4, -2.0, 2.0] {
            let analyzer = fixture.base(cx).with_robin_feedback(cx, &["air"],
                &[330.0*(1.0-slope)], &[slope], feedback_config()).unwrap();
            let initial = vec![300.0; fixture.mesh.vertex_count()];
            let prepared = analyzer.response_iterations();
            let result = analyzer.solve_maximum_to_goal(cx, &initial, &fixture.vertices(), controls()).unwrap();
            assert_eq!(result.stop, LinearGoalStop::GoalTolerance, "{result:?}");
            assert!(result.primal_iterations > 0 && result.primal_iterations <= controls().max_primal_iterations);
            assert!(result.analysis.meets_absolute_tolerance(controls().absolute_tolerance));
            assert_eq!(result.analysis.response_iterations(), prepared);
            for &temperature in &result.temperature { assert!((temperature-330.0).abs() < 1e-6); }
            if slope.abs() > 1.0 {
                assert_eq!(result.analysis.coupled().inverse_method(), Some(FeedbackInverseMethod::PortSchurDominance));
            }
            assert_eq!(result.analysis, analyzer.analyze_maximum(cx, &result.temperature, &fixture.vertices()).unwrap());
            assert_eq!(initial, vec![300.0; initial.len()]);
        }
    });
}

#[test]
fn nonsymmetric_upstream_feedback_and_relocating_maximum_match_dense_solution() {
    let mut fixture = Fixture::new(false, true);
    fixture.source = ScalarField::Nodal(fixture.mesh.positions().iter().map(|p| 4.0+2.0*p[0]).collect());
    with_cx(|cx| {
        let analyzer = fixture.base(cx).with_robin_feedback(cx, &["first", "second"],
            &[198.0,132.0], &[0.4,0.0,0.2,0.4], feedback_config()).unwrap();
        let (a,rhs,b,c,d) = analyzer.stored_system();
        let exact = dense(a,rhs,b,c,d);
        let exact_max = exact.iter().copied().fold(f64::NEG_INFINITY,f64::max);
        let cold = (0..exact.len()).min_by(|&i,&j| exact[i].total_cmp(&exact[j])).unwrap();
        assert!(exact_max-exact[cold] > 1e-3);
        let mut initial = vec![300.0; exact.len()]; initial[cold] = 350.0;
        let mut asymmetric = false;
        for i in 0..exact.len() { for j in 0..exact.len() {
            let ij = a.get(i,j)-(0..d.len()).map(|k| b.get(i,k)*c.get(k,j)).sum::<f64>();
            let ji = a.get(j,i)-(0..d.len()).map(|k| b.get(j,k)*c.get(k,i)).sum::<f64>();
            asymmetric |= (ij-ji).abs()>1e-4;
        } }
        assert!(asymmetric, "this test must require a nonsymmetric operator");
        let result = analyzer.solve_maximum_to_goal(cx,&initial,&fixture.vertices(),controls()).unwrap();
        assert_eq!(result.stop,LinearGoalStop::GoalTolerance,"{result:?}");
        let [lo,hi] = result.analysis.interval_k().unwrap();
        assert!(lo <= exact_max && exact_max <= hi);
        for (&found,&expected) in result.temperature.iter().zip(&exact) { assert!((found-expected).abs()<1e-6); }
        assert!(result.analysis.nominal_k()-result.temperature[cold]>1e-3);
        let mut reversed = fixture.vertices(); reversed.reverse();
        assert_eq!(result,analyzer.solve_maximum_to_goal(cx,&initial,&reversed,controls()).unwrap());
    });
}

#[test]
fn initial_goal_and_prescribed_selection_can_finish_without_krylov_work() {
    for fixed in [false,true] {
        let fixture = Fixture::new(fixed,false);
        with_cx(|cx| {
            let dofs = DofMap::new(&fixture.boundary,fixture.mesh.vertex_count()).unwrap();
            let analyzer = fixture.base(cx).with_robin_feedback(cx,&["air"],&[198.0],&[0.4],feedback_config()).unwrap();
            let initial = vec![300.0;fixture.mesh.vertex_count()];
            let selection = if fixed { dofs.fixed().to_vec() } else { fixture.vertices() };
            let bound = analyzer.analyze_maximum(cx,&initial,&selection).unwrap().algebraic_half_width_k().unwrap();
            let mut config = controls(); config.absolute_tolerance = (2.0*bound).max(1e-6); config.max_primal_iterations=0;
            let result = analyzer.solve_maximum_to_goal(cx,&initial,&selection,config).unwrap();
            assert_eq!(result.stop,LinearGoalStop::GoalTolerance);
            assert_eq!(result.primal_iterations,0); assert_eq!(result.goal_checks,1);
            assert_eq!(result.temperature,initial);
            if fixed { assert_eq!(result.analysis.algebraic_half_width_k(),Some(0.0)); }
        });
    }
}

#[test]
fn exhausted_work_keeps_the_checked_field_and_does_not_repeat_preparation() {
    let fixture=Fixture::new(true,false);
    with_cx(|cx| {
        let dofs=DofMap::new(&fixture.boundary,fixture.mesh.vertex_count()).unwrap();
        let analyzer=fixture.base(cx).with_robin_feedback(cx,&["air"],&[198.0],&[0.4],feedback_config()).unwrap();
        let initial=vec![300.0;fixture.mesh.vertex_count()];
        let before=analyzer.analyze_maximum(cx,&initial,&fixture.vertices()).unwrap();
        for cap in [0,1,3,17] {
            let mut config=controls();config.max_primal_iterations=cap;config.absolute_tolerance=1e-28;config.max_defect_corrections=1;
            let result=analyzer.solve_maximum_to_goal(cx,&initial,&fixture.vertices(),config).unwrap();
            assert_ne!(result.stop,LinearGoalStop::GoalTolerance);
            assert!(result.primal_iterations<=cap);assert!(result.defect_corrections<=1);
            assert!(result.analysis.algebraic_half_width_k().unwrap()<=before.algebraic_half_width_k().unwrap());
            assert_eq!(result.analysis.response_iterations(),before.response_iterations());
            assert_eq!(result.analysis,analyzer.analyze_maximum(cx,&result.temperature,&fixture.vertices()).unwrap());
            for &v in dofs.fixed() { assert_eq!(result.temperature[v],300.0); }
            if cap==0 { assert_eq!(result.temperature,initial); }
        }
    });
}

#[test]
fn missing_coupled_inverse_does_not_fall_back_to_frozen_solid_success() {
    let fixture=Fixture::new(false,false);
    with_cx(|cx| {
        // Unit reference feedback leaves the constant-temperature mode free.
        let analyzer=fixture.base(cx).with_robin_feedback(cx,&["air"],&[0.0],&[1.0],feedback_config()).unwrap();
        let initial=vec![300.0;fixture.mesh.vertex_count()];
        let mut config=controls();config.absolute_tolerance=1e100;
        let result=analyzer.solve_maximum_to_goal(cx,&initial,&fixture.vertices(),config).unwrap();
        assert_eq!(result.stop,LinearGoalStop::BoundUnavailable);
        assert_eq!(result.primal_iterations,0);assert_eq!(result.temperature,initial);
        assert!(result.analysis.interval_k().is_none());
    });
}

#[test]
fn observer_cancellation_refuses_even_a_passing_candidate_and_retry_replays() {
    let fixture=Fixture::new(false,false);
    with_cx(|cx| {
        let analyzer=fixture.base(cx).with_robin_feedback(cx,&["air"],&[198.0],&[0.4],feedback_config()).unwrap();
        let initial=vec![300.0;fixture.mesh.vertex_count()];
        let expected=analyzer.solve_maximum_to_goal(cx,&initial,&fixture.vertices(),controls()).unwrap();
        for cancel_initial in [true,false] {
            let mut observed=false;
            with_gate(|gate,cancelled| {
                let mut config=controls();if cancel_initial { config.absolute_tolerance=1e9; }
                let stopped=analyzer.solve_maximum_to_goal_observed(cancelled,&initial,&fixture.vertices(),config,|iterations,_| {
                    if cancel_initial || iterations>0 { observed=true;gate.request(); }
                });
                assert!(matches!(stopped,Err(fs_conduction::ConductionError::Cancelled{..})));
            });
            assert!(observed);
            assert_eq!(expected,analyzer.solve_maximum_to_goal(cx,&initial,&fixture.vertices(),controls()).unwrap());
        }
    });
}

#[test]
fn malformed_controls_fields_and_selections_refuse_before_observation() {
    let fixture=Fixture::new(true,false);
    with_cx(|cx| {
        let analyzer=fixture.base(cx).with_robin_feedback(cx,&["air"],&[198.0],&[0.4],feedback_config()).unwrap();
        let initial=vec![300.0;fixture.mesh.vertex_count()];
        for (tol,every) in [(f64::NAN,8),(0.0,8),(1e-6,0),(1e-6,33)] {
            let mut config=controls();config.absolute_tolerance=tol;config.check_every=every;
            assert!(analyzer.solve_maximum_to_goal_observed(cx,&initial,&fixture.vertices(),config,|_,_|panic!("invalid policy observed")).is_err());
        }
        for selected in [vec![],vec![0,0],vec![usize::MAX]] {
            assert!(analyzer.solve_maximum_to_goal(cx,&initial,&selected,controls()).is_err());
        }
        let dofs=DofMap::new(&fixture.boundary,fixture.mesh.vertex_count()).unwrap();
        let mut changed=initial;changed[dofs.fixed()[0]]=301.0;
        assert!(analyzer.solve_maximum_to_goal(cx,&changed,&fixture.vertices(),controls()).is_err());
    });
}
