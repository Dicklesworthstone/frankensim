//! Physical cooling gates participate in the real correction loop, rather
//! than discarding the first numerically goal-accepted field afterwards.
use fs_airflow::conjugate::{AirPath, AirSegment};
use fs_airflow::conjugate::goal::maximum::{SpectralMaximumControl,
    physical::{PhysicalAirMaximumPolish, PhysicalCoolingGates, PhysicalCoolingRefusal,
        polish_linear_maximum_with_spectral}};
use fs_conduction::{ConductionMesh, ConductionProblem, ConductionSolution, ConductivityModel,
    InitialGuess, LinearConfig, ScalarField, SolveConfig, ThermalBc, ThermalBoundary,
    ThermalBoundaryBuilder};
use fs_conduction::adjoint::{LinearGoalAnalysisConfig, LinearGoalSolveConfig, LinearGoalStop,
    RobinFeedbackAnalysisConfig, SpectralInverseLimits};
use fs_solver::goal::{GoalResidualLimits, feedback::FeedbackResidualLimits};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&gate, &Cx::new(&gate, arena,
            StreamKey { seed: 7308, kernel_id: 73, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic));
    });
}
struct Fixture {
    mesh: ConductionMesh, boundary: ThermalBoundary,
    material: ConductivityModel, source: ScalarField, paths: Vec<AirPath>,
}
impl Fixture {
    fn new() -> Self {
        let (complex, points) = fs_conduction::fixtures::unit_cube(1);
        let mesh = ConductionMesh::new(complex, points).unwrap();
        let area = mesh.boundary().iter().map(|face| face.area).sum();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("air", |_| true, ThermalBc::robin(2.0, 300.0).unwrap())
            .unwrap().finish().unwrap();
        let paths = vec![AirPath::new(330.0, 1.0, 24.0,
            vec![AirSegment::new("air", area, 2.0).unwrap()]).unwrap()];
        Self { mesh, boundary, paths, material: ConductivityModel::isotropic_declared(10.0).unwrap(),
            source: ScalarField::Uniform(6.0) }
    }
    fn problem(&self) -> ConductionProblem<'_> {
        ConductionProblem { mesh: &self.mesh, boundary: &self.boundary,
            material: &self.material, source: &self.source, element_materials: None }
    }
    fn original(&self, cx: &Cx<'_>) -> ConductionSolution {
        let mut config = SolveConfig { initial: InitialGuess::Uniform(300.0), ..SolveConfig::default() };
        config.stop.residual_rtol = 1e-8;
        fs_conduction::solve(cx, self.problem(), config).unwrap()
    }
    fn polish(&self, cx: &Cx<'_>, original: &ConductionSolution, goal: f64,
        iterations: usize, gates: PhysicalCoolingGates) -> PhysicalAirMaximumPolish
    {
        let limits = GoalResidualLimits { max_rows: 100, max_nonzeros: 100_000 };
        polish_linear_maximum_with_spectral(cx, self.problem(), None, &self.paths,
            LinearConfig { tolerance: 1e-11, max_iterations: 500, restart: 20 }, original,
            &(0..self.mesh.vertex_count()).collect::<Vec<_>>(),
            LinearGoalAnalysisConfig { residual_limits: limits, max_stability_iterations: 500 },
            RobinFeedbackAnalysisConfig { residual: FeedbackResidualLimits {
                solid: limits, max_ports: 8, max_transfer_nonzeros: 10_000,
                max_response_entries: 1000, max_verification_entries: 1_000_000 },
                max_response_iterations: 500, max_lowering_entries: 100_000 },
            LinearGoalSolveConfig { absolute_tolerance: goal, max_primal_iterations: iterations,
                check_every: 8, max_defect_corrections: 3 },
            SpectralMaximumControl { initial_shift: 0.001, limits: SpectralInverseLimits {
                system: limits, max_storage_entries: 100_000, max_work_entries: 1_000_000,
                max_shift_attempts: 12 } }, gates).unwrap()
    }
}
fn gates() -> PhysicalCoolingGates {
    PhysicalCoolingGates { energy_relative_tolerance: 1e-6, reference_tolerance_k: 1e-8,
        balance_tolerance_w: 1e-8, balance_relative_tolerance: 1e-7 }
}

#[test]
fn loose_temperature_goal_cannot_stop_before_physical_cooling_is_solved() {
    let f = Fixture::new();
    with_cx(|_, cx| {
        let original = f.original(cx);
        let before = original.clone();
        let result = f.polish(cx, &original, 1e9, 160, gates());
        assert!(result.correction.initial_analysis.meets_absolute_tolerance(1e9));
        let numerical = &result.correction.solution.solid;
        let accepted = result.accepted.as_ref().expect("continue after initial physical rejection");
        assert_eq!(numerical.stop, LinearGoalStop::GoalTolerance);
        assert!(numerical.primal_iterations > 0 && numerical.primal_iterations <= 160);
        assert!(result.physical_checks >= 2);
        assert!(result.physical_rejections >= 1);
        assert!(result.physical_checks <= numerical.goal_checks + 1);
        assert!(result.physical_refusal.is_none());
        assert_eq!(accepted.solid.temperature, numerical.temperature);
        assert!(accepted.solid.report.final_residual <= original.report.residual_threshold);
        assert!(accepted.solid.report.energy.relative_closure() <= gates().energy_relative_tolerance);
        assert!((accepted.air[0].total_heat_rate_w - 6.0).abs() < 1e-5);
        assert!((accepted.air[0].outlet_temperature_k - 330.25).abs() < 1e-6);
        assert!(accepted.stored_goal_met && accepted.temperature_changed);
        assert_eq!(accepted.solid.report.linear, original.report.linear);
        assert_eq!(original, before);
        assert_eq!(result, f.polish(cx, &original, 1e9, 160, gates()));
    });
}

#[test]
fn energy_gate_drives_correction_even_when_the_watt_residual_gate_is_loose() {
    let f = Fixture::new();
    with_cx(|_, cx| {
        let mut original = f.original(cx);
        original.report.residual_threshold = 1e6;
        let zero = f.polish(cx, &original, 1e9, 0, gates());
        assert!(matches!(zero.physical_refusal, Some(PhysicalCoolingRefusal::Solid(
            fs_conduction::ConductionError::Config { parameter: "corrected field energy", .. }))));
        assert_ne!(zero.correction.solution.solid.stop, LinearGoalStop::GoalTolerance);
        let result = f.polish(cx, &original, 1e9, 160, gates());
        let accepted = result.accepted.expect("energy rejection must continue correction");
        assert!(result.physical_rejections > 0);
        assert!(result.correction.solution.solid.primal_iterations > 0);
        assert!(accepted.solid.report.energy.relative_closure() <= gates().energy_relative_tolerance);
        assert!((accepted.air[0].total_heat_rate_w - 6.0).abs() < 1e-5);
    });
}

#[test]
fn no_primal_allowance_retains_the_real_rejection_not_a_goal_success() {
    let f = Fixture::new();
    with_cx(|_, cx| {
        let original = f.original(cx);
        let before = original.clone();
        let result = f.polish(cx, &original, 1e9, 0, gates());
        assert!(result.correction.solution.solid.analysis.meets_absolute_tolerance(1e9));
        assert_eq!(result.correction.solution.solid.stop, LinearGoalStop::IterationBudget);
        assert_eq!(result.correction.solution.solid.primal_iterations, 0);
        assert!(result.accepted.is_none() && result.physical_refusal.is_some());
        assert_eq!(result.physical_checks, 2); // Initial gate and final best candidate.
        assert_eq!(result.physical_rejections, 2);
        assert_eq!(result.correction.solution.solid.temperature, original.temperature);
        assert_eq!(original, before);
    });
}

#[test]
fn budget_candidate_can_be_physically_accepted_without_claiming_requested_goal_accuracy() {
    let f = Fixture::new();
    with_cx(|_, cx| {
        let original = f.original(cx);
        let result = f.polish(cx, &original, 1e-28, 16, gates());
        let accepted = result.accepted.expect("a full small-system solve passes physical gates");
        assert_ne!(result.correction.solution.solid.stop, LinearGoalStop::GoalTolerance);
        assert!(result.correction.solution.solid.primal_iterations <= 16);
        assert!(!accepted.stored_goal_met);
        assert_eq!(result.physical_checks, 1);
        assert_eq!(result.physical_rejections, 0);
        assert_eq!(accepted.solid.temperature, result.correction.solution.solid.temperature);
        assert!(accepted.solid.report.final_residual <= original.report.residual_threshold);
    });
}
