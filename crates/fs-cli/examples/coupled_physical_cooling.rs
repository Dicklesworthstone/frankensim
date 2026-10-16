//! Synthetic 6 W cube: correct the full cooling model and publish only after
//! production residual, energy and air checks. No experimental-validation claim.
//!
//! cargo run -p fs-cli --example coupled_physical_cooling

use fs_airflow::conjugate::{AirPath, AirSegment};
use fs_airflow::conjugate::goal::maximum::{SpectralMaximumControl,
    physical::{PhysicalCoolingGates, polish_linear_maximum_with_spectral}};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_conduction::{ConductionMesh, ConductionProblem, ConductivityModel, InitialGuess,
    LinearConfig, ScalarField, SolveConfig, ThermalBc, ThermalBoundaryBuilder};
use fs_conduction::adjoint::{LinearGoalAnalysisConfig, LinearGoalSolveConfig,
    RobinFeedbackAnalysisConfig, SpectralInverseLimits};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solver::goal::{GoalResidualLimits, feedback::FeedbackResidualLimits};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena,
            StreamKey { seed: 7307, kernel_id: 73, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        run(&cx)
    })
}

fn run(cx: &Cx<'_>) -> Result<(), Box<dyn std::error::Error>> {
    let (complex, positions) = fs_conduction::fixtures::unit_cube(1);
    let mesh = ConductionMesh::new(complex, positions)?;
    let area: f64 = mesh.boundary().iter().map(|face| face.area).sum();
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("air", |_| true, ThermalBc::robin(2.0, 300.0)?)?.finish()?;
    let material = ConductivityModel::isotropic_declared(10.0)?;
    let source = ScalarField::Uniform(6.0); // W/m^3 on a 1 m^3 synthetic cube.
    let problem = ConductionProblem { mesh: &mesh, boundary: &boundary, material: &material,
        element_materials: None, source: &source };
    let mut config = SolveConfig { initial: InitialGuess::Uniform(300.0), ..SolveConfig::default() };
    config.stop.residual_rtol = 1e-8;
    let original = fs_conduction::solve(cx, problem, config)?;
    // Independent inlet 330 K; capacity rate 24 W/K; h=2 W/(m^2 K).
    let paths = [AirPath::new(330.0, 1.0, 24.0, vec![AirSegment::new("air", area, 2.0)?])?];
    let limits = GoalResidualLimits { max_rows: 100, max_nonzeros: 100_000 };
    let result = polish_linear_maximum_with_spectral(cx, problem, None, &paths,
        LinearConfig { tolerance: 1e-11, max_iterations: 500, restart: 20 }, &original,
        &(0..mesh.vertex_count()).collect::<Vec<_>>(),
        LinearGoalAnalysisConfig { residual_limits: limits, max_stability_iterations: 500 },
        RobinFeedbackAnalysisConfig { residual: FeedbackResidualLimits { solid: limits,
            max_ports: 8, max_transfer_nonzeros: 10_000, max_response_entries: 1000,
            max_verification_entries: 1_000_000 }, max_response_iterations: 500, max_lowering_entries: 100_000 },
        LinearGoalSolveConfig { absolute_tolerance: 1e-7, max_primal_iterations: 500,
            check_every: 8, max_defect_corrections: 3 },
        SpectralMaximumControl { initial_shift: 0.001, limits: SpectralInverseLimits {
            system: limits, max_storage_entries: 100_000, max_work_entries: 1_000_000, max_shift_attempts: 12 } },
        PhysicalCoolingGates { energy_relative_tolerance: 1e-6, reference_tolerance_k: 1e-8,
            balance_tolerance_w: 1e-8, balance_relative_tolerance: 1e-7 })?;
    let accepted = result.accepted.as_ref().ok_or_else(||
        format!("physical cooling refused: {:?}", result.physical_refusal))?;
    if !accepted.stored_goal_met {
        return Err(format!("physical field accepted but goal unresolved: {:?}", result.correction.solution.solid.stop).into());
    }
    println!(concat!("{{\"schema\":\"frankensim.example.physical-cooling.v1\",",
        "\"model\":\"synthetic-numerical\",\"goal_met\":true,\"temperature_changed\":{},",
        "\"maximum_k\":{},\"physical_residual_w\":{},\"relative_energy_closure\":{},",
        "\"air_heat_w\":{},\"air_outlet_k\":{},\"correction_iterations\":{}}}"),
        accepted.temperature_changed, result.correction.solution.solid.analysis.nominal_k(),
        accepted.solid.report.final_residual, accepted.solid.report.energy.relative_closure(),
        accepted.air[0].total_heat_rate_w, accepted.air[0].outlet_temperature_k,
        result.correction.solution.solid.primal_iterations);
    Ok(())
}
