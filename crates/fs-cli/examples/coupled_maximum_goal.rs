//! Run with: cargo run -p fs-cli --example coupled_maximum_goal
//! Synthetic linear cube/air demonstration, not a validated cooling design.
use fs_airflow::conjugate::{AirPath, AirSegment};
use fs_airflow::conjugate::goal::maximum::solve_linear_maximum;
use fs_conduction::adjoint::{LinearGoalAnalysisConfig, LinearGoalSolveConfig,
    LinearGoalStop, RobinFeedbackAnalysisConfig};
use fs_conduction::{ConductionMesh, ConductionProblem, ConductivityModel,
    LinearConfig, ScalarField, ThermalBc, ThermalBoundaryBuilder};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_solver::goal::GoalResidualLimits;
use fs_solver::goal::feedback::FeedbackResidualLimits;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(&gate, arena,
            StreamKey { seed: 73, kernel_id: 740, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        run(&cx)
    })
}
fn run(cx: &Cx<'_>) -> Result<(), Box<dyn std::error::Error>> {
    let (complex, positions) = fs_conduction::fixtures::unit_cube(1);
    let mesh = ConductionMesh::new(complex, positions)?;
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region("air", |_| true, ThermalBc::robin(2.0, 300.0)?)?.finish()?;
    let material = ConductivityModel::isotropic_declared(10.0)?;
    let source = ScalarField::Uniform(0.0);
    let area = mesh.boundary().iter().map(|f| f.area).sum();
    let paths = [AirPath::new(330.0, 1.0, 24.0, vec![AirSegment::new("air", area, 2.0)?])?];
    let initial = vec![300.0; mesh.vertex_count()];
    let vertices: Vec<_> = (0..mesh.vertex_count()).collect();
    let limits = GoalResidualLimits { max_rows: 100, max_nonzeros: 100_000 };
    let result = solve_linear_maximum(cx, ConductionProblem {
        mesh: &mesh, boundary: &boundary, material: &material, element_materials: None, source: &source,
    }, None, &paths, LinearConfig { tolerance: 1e-12, max_iterations: 1000, restart: 20 },
        &initial, &vertices,
        LinearGoalAnalysisConfig { residual_limits: limits, max_stability_iterations: 1000 },
        RobinFeedbackAnalysisConfig {
            residual: FeedbackResidualLimits { solid: limits, max_ports: 4,
                max_transfer_nonzeros: 10_000, max_response_entries: 400,
                max_verification_entries: 1_000_000 },
            max_response_iterations: 1000, max_lowering_entries: 100_000,
        }, LinearGoalSolveConfig { absolute_tolerance: 1e-6, max_primal_iterations: 160,
            check_every: 8, max_defect_corrections: 2 })?;
    println!("stop={:?}; iterations={}; maximum={} K; stored-system error={:?} K",
        result.solid.stop, result.solid.primal_iterations,
        result.solid.analysis.nominal_k(), result.solid.analysis.algebraic_half_width_k());
    println!("fresh outlet={} K; air heat={} W; scope=stored linear coupled model",
        result.air[0].outlet_temperature_k, result.air[0].total_heat_rate_w);
    if result.solid.stop != LinearGoalStop::GoalTolerance {
        return Err("the declared maximum-goal budget was not met".into());
    }
    Ok(())
}
