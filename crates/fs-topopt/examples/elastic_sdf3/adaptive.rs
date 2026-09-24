//! The public adaptive continuation API on the same implicit cantilever.
use super::*;
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::elastic3::{ElasticityError3, adaptive::AdaptiveElasticity3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::QuadratureWork3;
use fs_solver::op::two_level::TwoLevelBudget;
use fs_topopt::sdf3::adaptive_continuation::{
    AdaptiveContinuationOptions3, controlled_adaptive_sdf3_continuation,
};
use fs_topopt::sdf3_goal::{
    GoalPreconditioner3, GoalReferenceLoad3, GoalRefinementError3, GoalRefinementOptions3,
};

pub fn run(updates: usize, iterations: usize) -> Result<(), Box<dyn std::error::Error>> {
    let mut tree = Octree3::uniform(1, 5, 2048)?;
    let domain = HexCell::try_new([0.0; 3], [1.0; 3])?;
    let material = IsotropicElastic::new(1.0, 0.3, 1.0)?;
    let limits = QuadratureOptions3::default();
    let mut geometry = QuadratureWork3::default();
    let mut build = |tree: &Octree3,
                     checkpoint: &mut dyn FnMut() -> ControlFlow<()>|
     -> Result<AdaptiveSolveSpace3, GoalRefinementError3> {
        // Keep one cumulative allowance across initial, probe and candidate
        // builds, including failed geometry. Each callback uses SolveControl.
        let options = QuadratureOptions3 {
            max_boxes: limits.max_boxes - geometry.boxes,
            max_points: limits.max_points - geometry.points,
            ..limits
        };
        let mut poll = |_| checkpoint();
        let mut quadrature =
            QuadratureControl3::new(options, &mut poll).map_err(ElasticityError3::from)?;
        let operator = AdaptiveElasticity3::build(
            domain,
            tree,
            &CurvedCantilever,
            &material,
            &|p| p[0] == 0.0,
            ElasticityOptions3::default(),
            &mut quadrature,
        );
        let spent = quadrature.work();
        geometry.boxes += spent.boxes;
        geometry.points += spent.points;
        geometry.field_evaluations += spent.field_evaluations;
        Ok(AdaptiveSolveSpace3::jacobi(operator?, 100_000_000))
    };
    let operator = build(&tree, &mut || ControlFlow::Continue(()))?;
    let mut study = CutDensityStudy3::new(operator, 0.15, SimpParams::default());
    let raw = vec![0.5; study.cells()];
    let y = |_: [f64; 3]| [0.0, -1.0, 0.0];
    let z = |_: [f64; 3]| [0.0, 0.0, -1.0];
    let loads = [
        GoalReferenceLoad3 {
            load: ReferenceLoad3::body(&y),
            weight: 0.3,
        },
        GoalReferenceLoad3 {
            load: ReferenceLoad3::body(&z),
            weight: 0.7,
        },
    ];
    let schedule = [(1.0, 1.0), (2.0, 2.0), (3.0, 8.0)].map(|(penal, beta)| SimpParams {
        penal,
        beta,
        ..Default::default()
    });
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(
        SolveBudget {
            total_iterations: iterations,
            ..Default::default()
        },
        &mut poll,
    );
    let report = controlled_adaptive_sdf3_continuation(
        &mut study,
        &mut tree,
        &loads,
        &raw,
        &schedule,
        AdaptiveContinuationOptions3 {
            optimization: MultiLoadOcOptions {
                max_iterations: updates,
                ..Default::default()
            },
            enrichment: GoalRefinementOptions3 {
                preconditioner: GoalPreconditioner3::TwoLevel {
                    budget: TwoLevelBudget::default(),
                    max_diagonal_contributions: 100_000_000,
                },
                ..Default::default()
            },
            ..Default::default()
        },
        &mut control,
        &mut build,
    );
    println!("stage,penal,beta,iteration,compliance,volume_fraction,max_change");
    for stage in &report.continuation.stages {
        for row in &stage.history {
            println!(
                "{},{},{},{},{:.17e},{:.17e},{:.17e}",
                stage.stage,
                stage.params.penal,
                stage.params.beta,
                row.iteration,
                row.compliance,
                row.volume_fraction,
                row.max_change
            );
        }
        let check = stage.gradient_check.as_ref().expect("gradient-gated stage");
        for probe in &check.probes {
            eprintln!(
                "stage={} direction={:?} gradient_relative_error={:.9e} volume_gradient_relative_error={:.9e} passed={}",
                stage.stage,
                probe.direction,
                probe.compliance_relative_error,
                probe.volume_relative_error,
                check.passed()
            );
        }
        eprintln!(
            "stage={} incoming_volume={:.9e} restoration_scale={:.9e}",
            stage.stage, stage.incoming_volume_fraction, stage.restoration_scale
        );
    }
    for step in &report.refinements {
        eprintln!(
            "destination_stage={} source_active_cells={} target_active_cells={} marked={:?} marking_fraction={:.6} marking_target_met={} two_grid_compliance_change={:.9e} installed={}",
            step.destination_stage,
            step.source_active_cells,
            step.target_active_cells,
            step.marking.marked,
            step.marking.achieved_fraction,
            step.marking.target_met,
            step.estimate.correction,
            step.installed
        );
    }
    let summary = &report.continuation;
    eprintln!(
        "continuation_stop={:?}; stopped_stage={:?}; evaluation_stop={:?}; refinement_error={:?}; cumulative_linear_iterations={}; geometry_points={}; background_refinements={}; implicit_boundary_rebuilt=false; continuum_certified=false; cross_model_descent_claimed=false",
        summary.termination,
        summary.stopped_stage,
        summary.evaluation_stop,
        report.refinement_error,
        summary.work.linear_iterations,
        geometry.points,
        report.refinements.iter().filter(|r| r.installed).count()
    );
    if let Some(last) = &summary.last {
        println!("level,cell_x,cell_y,cell_z,raw_density,projected_density");
        for ((cell, raw), projected) in study
            .operator()
            .elasticity()
            .leaves()
            .iter()
            .zip(&last.rho)
            .zip(&last.projected_rho)
        {
            let [x, y, z] = cell.index();
            println!("{},{x},{y},{z},{raw:.17e},{projected:.17e}", cell.level());
        }
    }
    if summary.termination != ContinuationTermination::ScheduleComplete {
        return Err(
            "adaptive continuation stopped; output contains only the retained solved prefix".into(),
        );
    }
    Ok(())
}
