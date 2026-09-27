//! Real coupled FEM correction must not stop at a loose goal before a
//! consumer's physical residual gate. No replacement numerical solver.
use fs_conduction::adjoint::{LinearGoalAnalysisConfig, LinearGoalAnalyzer,
    LinearGoalSolveConfig, LinearGoalStop, LinearRobinFeedbackAnalyzer,
    RobinFeedbackAnalysisConfig};
use fs_conduction::{ConductionError, ConductionMesh, ConductionProblem, ConductivityModel,
    LinearConfig, ScalarField, ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solver::goal::{GoalResidualLimits, feedback::FeedbackResidualLimits};

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&gate, &Cx::new(&gate, arena,
            StreamKey { seed: 7306, kernel_id: 73, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic));
    });
}
struct Fixture {
    mesh: ConductionMesh, boundary: ThermalBoundary,
    material: ConductivityModel, source: ScalarField,
}
impl Fixture {
    fn new() -> Self {
        let (complex, points) = fs_conduction::fixtures::unit_cube(1);
        let mesh = ConductionMesh::new(complex, points).unwrap();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("air", |_| true, ThermalBc::robin(2.0, 300.0).unwrap())
            .unwrap().finish().unwrap();
        Self { mesh, boundary, material: ConductivityModel::isotropic_declared(10.0).unwrap(),
            source: ScalarField::Uniform(6.0) }
    }
    fn analyzer(&self, cx: &Cx<'_>, slope: f64) -> LinearRobinFeedbackAnalyzer<'_> {
        let solid = GoalResidualLimits { max_rows: 100, max_nonzeros: 100_000 };
        LinearGoalAnalyzer::new_for_maximum(cx, ConductionProblem {
            mesh: &self.mesh, boundary: &self.boundary, material: &self.material,
            source: &self.source, element_materials: None,
        }, None, LinearConfig { tolerance: 1e-12, max_iterations: 1000, restart: 20 },
            &[300.0; 8], LinearGoalAnalysisConfig {
                residual_limits: solid, max_stability_iterations: 1000,
            }).unwrap().with_robin_feedback(cx, &["air"], &[330.0*(1.0-slope)], &[slope],
                RobinFeedbackAnalysisConfig { residual: FeedbackResidualLimits {
                    solid, max_ports: 4, max_transfer_nonzeros: 10_000,
                    max_response_entries: 400, max_verification_entries: 10_000_000,
                }, max_response_iterations: 1000, max_lowering_entries: 100_000 }).unwrap()
    }
}
fn control(cap: usize) -> LinearGoalSolveConfig {
    LinearGoalSolveConfig { absolute_tolerance: 1e9, max_primal_iterations: cap,
        check_every: 8, max_defect_corrections: 2 }
}
fn residual(analyzer: &LinearRobinFeedbackAnalyzer<'_>, x: &[f64]) -> f64 {
    let (a, rhs, b, c, d) = analyzer.stored_system();
    let ports: Vec<_> = (0..d.len()).map(|k| d[k]
        + (0..x.len()).map(|j| c.get(k,j)*x[j]).sum::<f64>()).collect();
    (0..x.len()).map(|i| {
        let ax = (0..x.len()).map(|j| a.get(i,j)*x[j]).sum::<f64>();
        (rhs[i] + (0..ports.len()).map(|k| b.get(i,k)*ports[k]).sum::<f64>() - ax).abs()
    }).fold(0.0, f64::max)
}

#[test]
fn physical_gate_continues_past_a_goal_accepted_initial_field() {
    let fixture = Fixture::new();
    with_cx(|_, cx| {
        for slope in [0.4, -2.0, 2.0] {
            let analyzer = fixture.analyzer(cx, slope);
            let initial = [300.0; 8];
            let region: Vec<_> = (0..8).collect();
            let ordinary = analyzer.solve_maximum_to_goal(cx, &initial, &region, control(160)).unwrap();
            assert_eq!(ordinary.stop, LinearGoalStop::GoalTolerance);
            assert_eq!(ordinary.primal_iterations, 0);
            assert!(residual(&analyzer, &ordinary.temperature) > 1.0);
            let prepared = analyzer.response_iterations();
            let mut calls = 0;
            let solved = analyzer.solve_maximum_to_goal_admitted(cx, &initial, &region, control(160), |t, a| {
                calls += 1;
                assert_eq!(*a, analyzer.analyze_maximum(cx, t, &region).unwrap());
                Ok::<_, ConductionError>(residual(&analyzer, t) <= 1e-8)
            }).unwrap();
            assert_eq!(solved.stop, LinearGoalStop::GoalTolerance, "{solved:?}");
            assert!(solved.primal_iterations > 0 && solved.primal_iterations <= 160);
            assert!(calls > 1 && calls <= solved.goal_checks);
            assert!(residual(&analyzer, &solved.temperature) <= 1e-8);
            assert_eq!(prepared, analyzer.response_iterations());
            assert_eq!(initial, [300.0; 8]);
        }
    });
}

#[test]
fn failed_gate_never_claims_success_or_resets_the_iteration_budget() {
    let fixture = Fixture::new();
    with_cx(|_, cx| {
        let analyzer = fixture.analyzer(cx, 0.4);
        for cap in [0, 1, 3, 17] {
            let mut calls = 0;
            let result = analyzer.solve_maximum_to_goal_admitted(cx, &[300.0; 8],
                &(0..8).collect::<Vec<_>>(), control(cap), |_, _| {
                    calls += 1; Ok::<_, ConductionError>(false)
                }).unwrap();
            assert_ne!(result.stop, LinearGoalStop::GoalTolerance);
            assert!(result.primal_iterations <= cap);
            assert!(result.defect_corrections <= control(cap).max_defect_corrections);
            assert!(calls > 0 && calls <= result.goal_checks);
            if cap == 0 {
                assert_eq!(calls, 1);
                assert_eq!(result.temperature, vec![300.0; 8]);
                assert_eq!(result.stop, LinearGoalStop::IterationBudget);
            }
        }
    });
}

#[test]
fn always_admit_preserves_the_existing_numerical_path_bit_for_bit() {
    let fixture = Fixture::new();
    with_cx(|_, cx| {
        let analyzer = fixture.analyzer(cx, 0.4);
        let region: Vec<_> = (0..8).collect();
        for tolerance in [1e9, 1e-6, 1e-28] {
            let cfg = LinearGoalSolveConfig { absolute_tolerance: tolerance, ..control(160) };
            let old = analyzer.solve_maximum_to_goal(cx, &[300.0; 8], &region, cfg).unwrap();
            let new = analyzer.solve_maximum_to_goal_admitted(cx, &[300.0; 8], &region, cfg,
                |_, _| Ok::<_, ConductionError>(true)).unwrap();
            assert_eq!(old, new);
        }
    });
}

#[test]
fn cancellation_after_a_gate_overrides_acceptance_rejection_and_callback_error() {
    let fixture = Fixture::new();
    for kind in 0..3 {
        with_cx(|gate, cx| {
            let analyzer = fixture.analyzer(cx, 0.4);
            let result = analyzer.solve_maximum_to_goal_admitted(cx, &[300.0; 8],
                &(0..8).collect::<Vec<_>>(), control(160), |_, _| {
                    gate.request();
                    if kind == 2 { Err(ConductionError::Config {
                        parameter: "test gate", what: "failed callback".into() }) }
                    else { Ok(kind == 0) }
                });
            assert!(matches!(result, Err(ConductionError::Cancelled { .. })));
        });
    }
}

#[test]
fn typed_gate_errors_propagate_and_retry_replays_without_repreparing() {
    #[derive(Debug, PartialEq)]
    enum Error { Numerical(ConductionError), SensorUnavailable }
    impl From<ConductionError> for Error { fn from(e: ConductionError) -> Self { Self::Numerical(e) } }
    let fixture = Fixture::new();
    with_cx(|_, cx| {
        let analyzer = fixture.analyzer(cx, 0.4);
        let region: Vec<_> = (0..8).collect();
        let first = analyzer.solve_maximum_to_goal(cx, &[300.0; 8], &region, control(160)).unwrap();
        let failed = analyzer.solve_maximum_to_goal_admitted(cx, &[300.0; 8], &region, control(160),
            |_, _| Err::<bool, _>(Error::SensorUnavailable));
        assert_eq!(failed.unwrap_err(), Error::SensorUnavailable);
        assert_eq!(first, analyzer.solve_maximum_to_goal(cx, &[300.0; 8], &region, control(160)).unwrap());
    });
}

#[test]
fn a_missing_inverse_never_reaches_the_consumer_success_gate() {
    let fixture = Fixture::new();
    with_cx(|_, cx| {
        let analyzer = fixture.analyzer(cx, 1.0);
        let result = analyzer.solve_maximum_to_goal_admitted(cx, &[300.0; 8],
            &(0..8).collect::<Vec<_>>(), control(160), |_, _| -> Result<bool, ConductionError> {
                panic!("a missing coupled bound cannot be admitted")
            }).unwrap();
        assert_eq!(result.stop, LinearGoalStop::BoundUnavailable);
        assert_eq!(result.primal_iterations, 0);
    });
}
