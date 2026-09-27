//! Complete production solid/air correction, not a frozen-solid stand-in.
use fs_airflow::conjugate::{AirPath, AirSegment};
use fs_airflow::conjugate::goal::maximum::{
    SpectralMaximumControl, solve_linear_maximum, solve_linear_maximum_with_spectral,
};
use fs_conduction::adjoint::{
    LinearGoalAnalysisConfig, LinearGoalSolveConfig, LinearGoalStop,
    RobinFeedbackAnalysisConfig, SpectralInverseLimits, SpectralStop,
};
use fs_conduction::{
    ConductionMesh, ConductionProblem, ConductivityModel, LinearConfig, ScalarField,
    ThermalBc, ThermalBoundary, ThermalBoundaryBuilder,
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solver::goal::{GoalResidualLimits, feedback::FeedbackResidualLimits};

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena,
            StreamKey { seed: 7307, kernel_id: 73, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        f(&gate, &cx);
    });
}
struct Fixture {
    mesh: ConductionMesh, boundary: ThermalBoundary, material: ConductivityModel,
    source: ScalarField, paths: Vec<AirPath>,
}
impl Fixture {
    fn new(cells: usize, shear: f64) -> Self {
        let (complex, mut points) = fs_conduction::fixtures::unit_cube(cells);
        for point in &mut points { point[0] += shear * point[1]; }
        let mesh = ConductionMesh::new(complex, points).unwrap();
        let area: f64 = mesh.boundary().iter().map(|face| face.area).sum();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("air", |_| true, ThermalBc::robin(2.0, 300.0).unwrap())
            .unwrap().finish().unwrap();
        Self { mesh, boundary, material: ConductivityModel::isotropic_declared(10.0).unwrap(),
            source: ScalarField::Uniform(0.0),
            paths: vec![AirPath::new(330.0, 1.0, 20.0 * area,
                vec![AirSegment::new("air", area, 2.0).unwrap()]).unwrap()] }
    }
    fn problem(&self) -> ConductionProblem<'_> {
        ConductionProblem { mesh: &self.mesh, boundary: &self.boundary,
            material: &self.material, source: &self.source, element_materials: None }
    }
}
fn configs() -> (LinearConfig, LinearGoalAnalysisConfig, RobinFeedbackAnalysisConfig,
    LinearGoalSolveConfig, SpectralMaximumControl) {
    let system = GoalResidualLimits { max_rows: 1000, max_nonzeros: 100_000 };
    (LinearConfig { tolerance: 1e-12, max_iterations: 2000, restart: 32 },
     LinearGoalAnalysisConfig { residual_limits: system, max_stability_iterations: 0 },
     RobinFeedbackAnalysisConfig {
        residual: FeedbackResidualLimits { solid: system, max_ports: 4,
            max_transfer_nonzeros: 10_000, max_response_entries: 10_000,
            max_verification_entries: 20_000_000 },
        max_response_iterations: 2000, max_lowering_entries: 100_000 },
     LinearGoalSolveConfig { absolute_tolerance: 1e-4, max_primal_iterations: 1000,
        check_every: 16, max_defect_corrections: 4 },
     SpectralMaximumControl { initial_shift: 0.001,
        limits: SpectralInverseLimits { system, max_storage_entries: 1_000_000,
            max_work_entries: 20_000_000, max_shift_attempts: 16 } })
}

#[test]
fn spectral_cooling_correction_reaches_the_coupled_equilibrium_and_rebuilds_air() {
    let fixture = Fixture::new(6, 1.5);
    with_cx(|_, cx| {
        let (linear, solid, feedback, control, spectral) = configs();
        let initial = vec![300.0; fixture.mesh.vertex_count()];
        let vertices: Vec<_> = (0..initial.len()).collect();
        let old = solve_linear_maximum(cx, fixture.problem(), None, &fixture.paths,
            linear, &initial, &vertices, solid, feedback, control).unwrap();
        assert_eq!(old.solid.stop, LinearGoalStop::BoundUnavailable);
        let out = solve_linear_maximum_with_spectral(cx, fixture.problem(), None, &fixture.paths,
            linear, &initial, &vertices, solid, feedback, control, spectral).unwrap();
        let preparation = out.preparation.as_ref().expect("must prepare the missing inverse");
        assert_eq!(preparation.stop, SpectralStop::Certified);
        assert!(preparation.work_entries > 0 && preparation.work_entries <= spectral.limits.max_work_entries);
        assert!(out.initial_analysis.algebraic_half_width_k().is_none());
        assert_eq!(out.solution.solid.stop, LinearGoalStop::GoalTolerance);
        assert!(out.solution.solid.primal_iterations > 0);
        assert!(out.solution.solid.primal_iterations <= control.max_primal_iterations);
        assert!(out.solution.solid.analysis.meets_absolute_tolerance(control.absolute_tolerance));
        assert!(out.solution.solid.goal_checks > 1);
        let reused = out.solution.solid.analysis.coupled().solid_spectral().unwrap();
        assert_eq!(reused.work_entries(), 0, "correction assessments must not refactor");
        assert_eq!(reused.shift_attempts(), 0);
        assert!(out.solution.solid.temperature.iter().all(|t| (*t - 330.0).abs() < 1e-4));
        assert!(initial.iter().all(|t| *t == 300.0));
        assert!((out.solution.wall_temperatures_k[0] - 330.0).abs() < 1e-4);
        assert_eq!(out.solution.air[0], fixture.paths[0].march(&out.solution.wall_temperatures_k).unwrap());
        assert_ne!(out.solution.air, old.air, "air states must belong to the corrected field");
    });
}

#[test]
fn failed_preparation_and_zero_primal_allowance_keep_honest_partial_results() {
    let fixture = Fixture::new(1, 1.5);
    with_cx(|_, cx| {
        let (linear, solid, feedback, mut control, mut spectral) = configs();
        let initial = vec![300.0; fixture.mesh.vertex_count()];
        let vertices: Vec<_> = (0..initial.len()).collect();
        spectral.limits.max_work_entries = 0;
        let stopped = solve_linear_maximum_with_spectral(cx, fixture.problem(), None, &fixture.paths,
            linear, &initial, &vertices, solid, feedback, control, spectral).unwrap();
        assert_eq!(stopped.preparation.unwrap().stop, SpectralStop::WorkLimit);
        assert_eq!(stopped.solution.solid.stop, LinearGoalStop::BoundUnavailable);
        assert_eq!(stopped.solution.solid.primal_iterations, 0);
        assert_eq!(stopped.solution.solid.temperature, initial);
        spectral = configs().4;
        control.max_primal_iterations = 0;
        let stopped = solve_linear_maximum_with_spectral(cx, fixture.problem(), None, &fixture.paths,
            linear, &initial, &vertices, solid, feedback, control, spectral).unwrap();
        assert_eq!(stopped.preparation.unwrap().stop, SpectralStop::Certified);
        assert_eq!(stopped.solution.solid.stop, LinearGoalStop::IterationBudget);
        assert_eq!(stopped.solution.solid.temperature, initial);
        assert_eq!(stopped.solution.solid.primal_iterations, 0);
        assert!(stopped.solution.solid.analysis.algebraic_half_width_k().unwrap() >= 30.0);
    });
}

#[test]
fn existing_inverse_route_matches_the_original_workflow_without_preparation() {
    let fixture = Fixture::new(1, 0.0);
    with_cx(|_, cx| {
        let (linear, solid, feedback, control, spectral) = configs();
        let initial = vec![300.0; fixture.mesh.vertex_count()];
        let vertices: Vec<_> = (0..initial.len()).collect();
        let old = solve_linear_maximum(cx, fixture.problem(), None, &fixture.paths,
            linear, &initial, &vertices, solid, feedback, control).unwrap();
        let new = solve_linear_maximum_with_spectral(cx, fixture.problem(), None, &fixture.paths,
            linear, &initial, &vertices, solid, feedback, control, spectral).unwrap();
        assert!(new.preparation.is_none());
        assert_eq!(new.solution, old);
    });
}

#[test]
fn invalid_proof_controls_and_cancellation_never_publish_a_temperature_bundle() {
    let fixture = Fixture::new(1, 1.5);
    with_cx(|gate, cx| {
        let (linear, solid, feedback, control, mut spectral) = configs();
        let initial = vec![300.0; fixture.mesh.vertex_count()];
        let vertices: Vec<_> = (0..initial.len()).collect();
        for shift in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            spectral.initial_shift = shift;
            assert!(solve_linear_maximum_with_spectral(cx, fixture.problem(), None, &fixture.paths,
                linear, &initial, &vertices, solid, feedback, control, spectral).is_err());
        }
        spectral = configs().4;
        spectral.limits.max_shift_attempts = 0;
        assert!(solve_linear_maximum_with_spectral(cx, fixture.problem(), None, &fixture.paths,
            linear, &initial, &vertices, solid, feedback, control, spectral).is_err());
        gate.request();
        assert!(solve_linear_maximum_with_spectral(cx, fixture.problem(), None, &fixture.paths,
            linear, &initial, &vertices, solid, feedback, control, configs().4).is_err());
        assert!(initial.iter().all(|t| *t == 300.0));
    });
}
