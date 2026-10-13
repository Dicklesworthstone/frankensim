//! Two selective robust design refinements after the material schedule.
use super::*;
use fs_cutfem::elastic3::{ElasticityError3, adaptive::AdaptiveElasticity3, surface::ReferenceLoad3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::QuadratureWork3;
use fs_topopt::sdf3::robust::adaptive::{RobustAdaptiveOptions3, RobustAdaptiveStop3, controlled_adaptive_robust_design};
use fs_topopt::sdf3_goal::{GoalPreconditioner3, GoalReferenceLoad3, GoalRefinementError3};

pub fn run(updates: usize, iterations: usize) -> Result<(), Box<dyn std::error::Error>> {
    let mut tree = Octree3::uniform(1, 5, 2048)?;
    let domain = HexCell::try_new([0.0; 3], [1.0; 3])?;
    let material = IsotropicElastic::new(1.0, 0.3, 1.0)?;
    let limits = QuadratureOptions3::default();
    let mut geometry = QuadratureWork3::default();
    let mut build = |tree: &Octree3, cp: &mut dyn FnMut() -> ControlFlow<()>|
        -> Result<AdaptiveSolveSpace3, GoalRefinementError3> {
        let remaining = QuadratureOptions3 {
            max_boxes: limits.max_boxes.saturating_sub(geometry.boxes),
            max_points: limits.max_points.saturating_sub(geometry.points), ..limits
        };
        let mut poll = |_| cp();
        let mut q = QuadratureControl3::new(remaining, &mut poll).map_err(ElasticityError3::from)?;
        let result = AdaptiveElasticity3::build(domain, tree, &CurvedCantilever, &material,
            &|p| p[0] == 0.0, ElasticityOptions3::default(), &mut q);
        let spent = q.work();
        geometry.boxes += spent.boxes; geometry.points += spent.points;
        geometry.field_evaluations += spent.field_evaluations;
        Ok(AdaptiveSolveSpace3::jacobi(result?, 100_000_000))
    };
    let schedule = [(1.0, 1.0), (2.0, 2.0), (3.0, 8.0)]
        .map(|(penal, beta)| SimpParams { penal, beta, ..Default::default() });
    let mut study = CutDensityStudy3::new(build(&tree, &mut || ControlFlow::Continue(()))?, 0.15, schedule[0]);
    let raw = vec![0.3; study.cells()];
    let y = |_: [f64; 3]| [0.0, -1.0, 0.0];
    let z = |_: [f64; 3]| [0.0, 0.0, -1.0];
    let loads = [GoalReferenceLoad3 { load: ReferenceLoad3::body(&y), weight: 0.3 },
        GoalReferenceLoad3 { load: ReferenceLoad3::body(&z), weight: 0.7 }];
    let mut options = RobustAdaptiveOptions3::default();
    options.refit.continuation.steps_per_stage = updates;
    // A declared fixed numerical conditioning scale, not a physical law or a
    // value recomputed mid-optimization. All grids use this SAME scale.
    options.refit.continuation.design.compliance_scale = 1000.0;
    options.enrichment.preconditioner = GoalPreconditioner3::Jacobi { max_contributions: 100_000_000 };
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: iterations, ..Default::default() }, &mut poll);
    let report = controlled_adaptive_robust_design(&mut study, &mut tree, &loads, &raw, &schedule,
        options, &mut control, &mut build)?;
    drop(build);
    eprintln!("units=dimensionless; seeds=none; filter_radius=0.15; initial_schedule={schedule:?}; adaptive_policy={options:?}; total_krylov_budget={iterations}; geometry_limits={limits:?}");
    println!("grid,background_cells,active_cells,stage,penal,beta,iteration,eroded_compliance,dilated_volume,feasible,selected");
    for (grid_index, grid) in report.grids.iter().enumerate() {
        for stage in &grid.design.stages {
            for (i, row) in stage.history.iter().enumerate() {
                println!("{grid_index},{},{},{},{},{},{},{:.17e},{:.17e},{},{}",
                    grid.background_cells, grid.active_cells, stage.stage, stage.params.penal, stage.params.beta,
                    row.iteration, row.eroded_compliance, row.dilated_volume_fraction,
                    row.feasible, i == stage.selected_iteration);
            }
            eprintln!("grid={grid_index}; stage={}; incoming_dilated_volume={:.17e}; restoration_scale={:.17e}; gradient_passed={}; optimizer_outcome={:?}",
                stage.stage, stage.incoming_dilated_volume_fraction, stage.restoration_scale,
                stage.gradient_check.numerical.passed(), stage.outcome);
        }
    }
    println!("source_grid,eroded_coarse_compliance,inherited_eroded_fine_compliance,max_case_relative_gap,installed");
    for probe in &report.probes {
        println!("{},{:.17e},{:.17e},{:.17e},{}", probe.source_grid, probe.estimate.coarse_value,
            probe.estimate.fine_value, probe.maximum_relative_gap, probe.installed);
        eprintln!("probe_source_grid={}; marking={:?}; target_background_cells={:?}",
            probe.source_grid, probe.marking, probe.target_background_cells);
    }
    eprintln!("adaptive_stop={:?}; error={:?}; rejected_candidate_stop={:?}; solve_work={:?}; geometry_work={geometry:?}; cross_grid_descent_claimed=false; refinement_limit_is_accuracy=false; continuum_bound=false; minimum_feature_size_certified=false",
        report.stop, report.error, report.rejected_design.as_ref().map(|r| r.stop), report.work);
    let Some(last) = report.last() else {
        return Err("no robust baseline was admitted; no design is exported".into());
    };
    let grid = report.grids.last().expect("retained point belongs to a grid");
    eprintln!("exported_model={:?}; exported_grid={}; exported_design=best-feasible; rejected_candidate_exported=false",
        grid.design.params, report.grids.len() - 1);
    println!("level,cell_x,cell_y,cell_z,raw_density,eroded_density,nominal_density,dilated_density");
    let op = study.operator().elasticity();
    for (i, key) in op.leaves().iter().enumerate() {
        let [x,y,z] = key.index();
        println!("{},{x},{y},{z},{:.17e},{:.17e},{:.17e},{:.17e}", key.level(), last.rho[i],
            last.evaluation.eroded.projected_rho[i], last.evaluation.nominal_rho[i], last.evaluation.dilated_rho[i]);
    }
    println!("node_x,node_y,node_z,eroded_uy_x,eroded_uy_y,eroded_uy_z,eroded_uz_x,eroded_uz_y,eroded_uz_z");
    for (i, p) in op.nodes().iter().enumerate() {
        let uy = &last.evaluation.eroded.objective.displacements[0];
        let uz = &last.evaluation.eroded.objective.displacements[1];
        println!("{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            p[0], p[1], p[2], uy[3*i], uy[3*i+1], uy[3*i+2], uz[3*i], uz[3*i+1], uz[3*i+2]);
    }
    if !matches!(report.stop, RobustAdaptiveStop3::RefinementLimit | RobustAdaptiveStop3::ErodedComplianceTolerance) {
        return Err("robust adaptation stopped early; exported fields describe the installed feasible prefix".into());
    }
    Ok(())
}
