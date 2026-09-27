//! End-to-end consumer of full-feedback correction AND physical publication.
use fs_airflow::conjugate::{AirPath, AirSegment};
use fs_airflow::conjugate::goal::maximum::{SpectralMaximumControl,
    physical::{PhysicalCoolingGates, PhysicalCoolingRefusal, polish_linear_maximum_with_spectral}};
use fs_conduction::{ConductionMesh, ConductionProblem, ConductionSolution, ConductivityModel,
    InitialGuess, LinearConfig, ScalarField, SolveConfig, ThermalBc, ThermalBoundary, ThermalBoundaryBuilder};
use fs_conduction::adjoint::{LinearGoalAnalysisConfig, LinearGoalSolveConfig,
    RobinFeedbackAnalysisConfig, SpectralInverseLimits};
use fs_solver::goal::{GoalResidualLimits, feedback::FeedbackResidualLimits};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};

fn with_cx(f: impl FnOnce(&CancelGate, &Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&gate, &Cx::new(
        &gate, arena, StreamKey { seed: 7307, kernel_id: 73, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic,
    )));
}
struct Fixture {
    mesh: ConductionMesh, boundary: ThermalBoundary, material: ConductivityModel,
    source: ScalarField, paths: Vec<AirPath>,
}
impl Fixture {
    fn new(multiple: bool) -> Self {
        let (mesh, points) = fs_conduction::fixtures::unit_cube(1);
        let mesh = ConductionMesh::new(mesh, points).unwrap();
        let boundary = if multiple {
            ThermalBoundaryBuilder::new(&mesh)
                .region("upstream", |f| f.centroid[0] < 1e-9, ThermalBc::robin(2.0, 300.0).unwrap()).unwrap()
                .region("downstream", |f| f.centroid[0] > 1.0 - 1e-9, ThermalBc::robin(2.0, 300.0).unwrap()).unwrap()
                .region("other", |f| f.centroid[1] < 1e-9, ThermalBc::robin(3.0, 300.0).unwrap()).unwrap()
                .region("static", |f| f.centroid[1] > 1.0 - 1e-9, ThermalBc::robin(1.0, 310.0).unwrap()).unwrap()
                .adiabatic_remainder().finish().unwrap()
        } else {
            ThermalBoundaryBuilder::new(&mesh)
                .region("air", |_| true, ThermalBc::robin(2.0, 300.0).unwrap()).unwrap().finish().unwrap()
        };
        let area = |name: &str| mesh.boundary().iter().enumerate().filter_map(|(slot, face)| {
            boundary.region_for(slot).filter(|&r| boundary.region_names()[r] == name).map(|_| face.area)
        }).sum();
        let segment = |name: &str, h| AirSegment::new(name, area(name), h).unwrap();
        let paths = if multiple { vec![
            AirPath::new(330.0, 1.0, 12.0, vec![segment("upstream", 2.0), segment("downstream", 2.0)]).unwrap(),
            AirPath::new(290.0, 1.0, 8.0, vec![segment("other", 3.0)]).unwrap(),
        ] } else { vec![AirPath::new(330.0, 1.0, 24.0, vec![segment("air", 2.0)]).unwrap()] };
        Self { mesh, boundary, paths, material: ConductivityModel::isotropic_declared(10.0).unwrap(),
            source: ScalarField::Uniform(6.0) }
    }
    fn problem(&self) -> ConductionProblem<'_> {
        ConductionProblem { mesh: &self.mesh, boundary: &self.boundary, material: &self.material,
            source: &self.source, element_materials: None }
    }
    fn original(&self, cx: &Cx<'_>) -> ConductionSolution {
        let mut config = SolveConfig { initial: InitialGuess::Uniform(300.0), ..SolveConfig::default() };
        config.stop.residual_rtol = 1e-8;
        fs_conduction::solve(cx, self.problem(), config).unwrap()
    }
    fn polish(&self, cx: &Cx<'_>, original: &ConductionSolution, iterations: usize, gates: PhysicalCoolingGates)
        -> Result<fs_airflow::conjugate::goal::maximum::physical::PhysicalAirMaximumPolish,
            fs_airflow::conjugate::goal::CoupledGoalError>
    {
        let limits = GoalResidualLimits { max_rows: 100, max_nonzeros: 100_000 };
        polish_linear_maximum_with_spectral(cx, self.problem(), None, &self.paths,
            LinearConfig { tolerance: 1e-11, max_iterations: 500, restart: 20 }, original,
            &(0..self.mesh.vertex_count()).collect::<Vec<_>>(),
            LinearGoalAnalysisConfig { residual_limits: limits, max_stability_iterations: 500 },
            RobinFeedbackAnalysisConfig { residual: FeedbackResidualLimits { solid: limits,
                max_ports: 8, max_transfer_nonzeros: 10_000, max_response_entries: 1000,
                max_verification_entries: 1_000_000 }, max_response_iterations: 500, max_lowering_entries: 100_000 },
            LinearGoalSolveConfig { absolute_tolerance: 1e-7, max_primal_iterations: iterations,
                check_every: 8, max_defect_corrections: 3 },
            SpectralMaximumControl { initial_shift: 0.001, limits: SpectralInverseLimits {
                system: limits, max_storage_entries: 100_000, max_work_entries: 1_000_000, max_shift_attempts: 12 } },
            gates)
    }
}
fn gates() -> PhysicalCoolingGates {
    PhysicalCoolingGates { energy_relative_tolerance: 1e-6, reference_tolerance_k: 1e-8,
        balance_tolerance_w: 1e-8, balance_relative_tolerance: 1e-7 }
}

#[test]
fn physical_solution_uses_corrected_boundary_fluxes_and_air_not_the_stale_frozen_report() {
    let f = Fixture::new(false);
    with_cx(|_, cx| {
        let original = f.original(cx);
        let saved = original.clone();
        let result = f.polish(cx, &original, 500, gates()).unwrap();
        let adopted = result.accepted.as_ref().expect("real corrected field passes physical gates");
        assert!(adopted.temperature_changed && adopted.stored_goal_met);
        assert!(result.physical_refusal.is_none());
        assert!(result.correction.solution.solid.primal_iterations > 0);
        assert_eq!(adopted.solid.temperature, result.correction.solution.solid.temperature);
        assert!(adopted.solid.report.final_residual <= original.report.residual_threshold);
        assert!(adopted.solid.report.energy.relative_closure() <= gates().energy_relative_tolerance);
        assert!((adopted.air[0].total_heat_rate_w - 6.0).abs() < 1e-6);
        assert!((adopted.air[0].outlet_temperature_k - 330.25).abs() < 1e-6);
        assert_eq!(adopted.solid.report.linear, original.report.linear);
        let fresh = fs_conduction::solve(cx, ConductionProblem { boundary: &adopted.boundary, ..f.problem() },
            SolveConfig { initial: InitialGuess::Uniform(330.0), ..SolveConfig::default() }).unwrap();
        for (a, b) in fresh.temperature.iter().zip(&adopted.solid.temperature) { assert!((a-b).abs() < 1e-6); }
        assert_eq!(original, saved);
        assert_eq!(result, f.polish(cx, &original, 500, gates()).unwrap());
    });
}

#[test]
fn independent_branches_and_static_robin_rows_share_one_fresh_physical_solution() {
    let f = Fixture::new(true);
    with_cx(|_, cx| {
        let original = f.original(cx);
        let result = f.polish(cx, &original, 500, gates()).unwrap();
        let adopted = result.accepted.expect("multi-branch physical result");
        assert!(adopted.stored_goal_met);
        assert_eq!(adopted.branches.len(), 2);
        assert_eq!(adopted.air[0].segments[0].inlet_temperature_k, 330.0);
        assert_eq!(adopted.air[1].segments[0].inlet_temperature_k, 290.0);
        let mut offset = 0;
        for (path, actual) in f.paths.iter().zip(&adopted.air) {
            let end = offset + path.segments().len();
            assert_eq!(*actual, path.march(&adopted.wall_temperatures_k[offset..end]).unwrap());
            offset = end;
        }
        let static_row = adopted.solid.report.robin_fluxes.iter().find(|r| r.region == "static").unwrap();
        assert!((static_row.mean_reference_temperature_k - 310.0).abs() < 1e-10);
        let total: f64 = adopted.solid.report.robin_fluxes.iter().map(|r| r.heat_rate_w).sum();
        assert!((total - adopted.solid.report.energy.robin_out_w).abs() <= 1e-7);
        for branch in &adopted.branches {
            assert!(branch.interface_imbalance_w <= branch.watt_limit);
            assert!(branch.enthalpy_imbalance_w <= branch.watt_limit);
        }
    });
}

#[test]
fn uncorrected_stale_field_fails_physical_gates_without_losing_paid_work() {
    let f = Fixture::new(false);
    with_cx(|_, cx| {
        let original = f.original(cx);
        let saved = original.clone();
        let result = f.polish(cx, &original, 0, gates()).unwrap();
        assert!(result.accepted.is_none());
        assert!(matches!(result.physical_refusal, Some(PhysicalCoolingRefusal::Solid(
            fs_conduction::ConductionError::NotConverged { .. }))));
        assert_eq!(result.correction.solution.solid.primal_iterations, 0);
        assert!(result.correction.solution.solid.analysis.response_iterations() > 0);
        assert_eq!(original, saved);
        // Passing a loose residual alone cannot authorize a stale energy report.
        let mut loose = original.clone(); loose.report.residual_threshold = 1e6;
        let result = f.polish(cx, &loose, 0, gates()).unwrap();
        assert!(result.accepted.is_none());
        assert!(matches!(result.physical_refusal, Some(PhysicalCoolingRefusal::Solid(
            fs_conduction::ConductionError::Config { parameter: "corrected field energy", .. }))));
    });
}

#[test]
fn invalid_acceptance_and_cancellation_return_no_partial_bundle() {
    let f = Fixture::new(false);
    with_cx(|gate, cx| {
        let original = f.original(cx);
        for changed in [PhysicalCoolingGates { reference_tolerance_k: 0.0, ..gates() },
            PhysicalCoolingGates { energy_relative_tolerance: f64::NAN, ..gates() },
            PhysicalCoolingGates { balance_relative_tolerance: 1.0, ..gates() }]
        { assert!(f.polish(cx, &original, 500, changed).is_err()); }
        gate.request();
        assert!(f.polish(cx, &original, 500, gates()).is_err());
    });
}
