//! Feasible robust continuation, then a same-eroded-material two-grid probe.
use super::*;
use fs_cutfem::elastic3::adaptive::{AdaptiveElasticity3,
    enrichment::precondition::AdaptiveSolveSpace3};
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::octree3::Octree3;
use fs_topopt::sdf3::robust::{RobustDesignOptions3, RobustOptions3};
use fs_topopt::sdf3::robust::continuation::{RobustContinuationOptions3,
    RobustContinuationStop3, controlled_sdf3_robust_continuation};
use fs_topopt::sdf3_goal::{GoalPreconditioner3, GoalReferenceLoad3, GoalRefinementOptions3};

pub fn run(updates: usize, iterations: usize) -> Result<(), Box<dyn std::error::Error>> {
    let tree = Octree3::uniform(1, 5, 2048)?;
    let domain = HexCell::try_new([0.0; 3], [1.0; 3])?;
    let material = IsotropicElastic::new(1.0, 0.3, 1.0)?;
    let limits = QuadratureOptions3::default();
    let mut geometry_poll = |_| ControlFlow::Continue(());
    // The initial and diagnostic builds share ONE geometry allowance.
    let mut quadrature = QuadratureControl3::new(limits, &mut geometry_poll)?;
    let operator = AdaptiveElasticity3::build(domain, &tree, &CurvedCantilever,
        &material, &|p| p[0] == 0.0, ElasticityOptions3::default(), &mut quadrature)?;
    let y = |_: [f64; 3]| [0.0, -1.0, 0.0];
    let z = |_: [f64; 3]| [0.0, 0.0, -1.0];
    let fy = operator.body_load(&y, || ControlFlow::Continue(()))?;
    let fz = operator.body_load(&z, || ControlFlow::Continue(()))?;
    let keys = operator.leaves().to_vec(); let nodes = operator.nodes().to_vec();
    let schedule = [(1.0, 1.0), (2.0, 2.0), (3.0, 8.0)]
        .map(|(penal, beta)| SimpParams { penal, beta, ..Default::default() });
    let mut study = CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(operator, 100_000_000),
        0.15, schedule[0]);
    let raw = vec![0.3; study.cells()];
    let loads = [LoadCase { force: &fy, weight: 0.3 }, LoadCase { force: &fz, weight: 0.7 }];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: iterations,
        ..Default::default() }, &mut poll);
    let robustness = RobustOptions3::default();
    // Explicit conditioning scale fixed ONCE from a charged physical solve.
    let baseline = study.evaluate_robust(&raw, &loads, robustness, &mut control)?;
    let design = RobustDesignOptions3 { robustness, volume_limit: 0.5,
        compliance_scale: baseline.eroded.objective.compliance, ..Default::default() };
    let options = RobustContinuationOptions3 { design, steps_per_stage: updates, ..Default::default() };
    let report = controlled_sdf3_robust_continuation(&mut study, &loads, &raw, &schedule, options, &mut control)?;
    eprintln!("units=dimensionless; seeds=none; filter_radius=0.15; threshold_offset={}; density_floor={}; dilated_volume_cap={}; compliance_scale={:.17e}; steps_per_stage={updates}; optimizer_evaluations_per_stage={}; tolerance={}; total_krylov_budget={iterations}; geometry_box_budget={}; geometry_point_budget={}",
        robustness.threshold_offset, design.density_floor, design.volume_limit, design.compliance_scale,
        design.optimizer.max_evaluations, design.optimizer.tolerance, limits.max_boxes, limits.max_points);
    println!("stage,penal,beta,iteration,eroded_compliance,eroded_volume,nominal_volume,dilated_volume,feasible,selected");
    for stage in &report.stages {
        for (index, row) in stage.history.iter().enumerate() {
            println!("{},{},{},{},{:.17e},{:.17e},{:.17e},{:.17e},{},{}", stage.stage,
                stage.params.penal, stage.params.beta, row.iteration, row.eroded_compliance,
                row.eroded_volume_fraction, row.nominal_volume_fraction, row.dilated_volume_fraction,
                row.feasible, index == stage.selected_iteration);
        }
        eprintln!("stage={}; eta={}; e_min={}; incoming_dilated_volume={:.17e}; restoration_scale={:.17e}; gradient_passed={}; final_al_outcome={:?}",
            stage.stage, stage.params.eta, stage.params.e_min, stage.incoming_dilated_volume_fraction,
            stage.restoration_scale, stage.gradient_check.numerical.passed(), stage.outcome);
    }
    eprintln!("continuation_stop={:?}; stopped_stage={:?}; error={:?}; rejected_gradient_check={:?}; solve_work={:?}; geometry_work={:?}; cross_model_descent_claimed=false; schedule_completion_is_convergence=false",
        report.stop, report.stopped_stage, report.error, report.rejected_gradient_check, report.work, quadrature.work());
    let Some(last) = &report.last else {
        return Err("no feasible gradient-checked robust baseline was admitted; no design is exported".into());
    };
    let evaluated = &last.evaluation;
    eprintln!("exported_stage={}; exported_penal={}; exported_beta={}; exported_eta={}; exported_e_min={}; exported_design=latest-stage-best-feasible; minimum_feature_size_certified=false; continuum_safety_certified=false; global_optimum_certified=false",
        report.stages.last().expect("admitted point has a stage").stage, report.params.penal,
        report.params.beta, report.params.eta, report.params.e_min);
    println!("level,cell_x,cell_y,cell_z,raw_density,eroded_density,nominal_density,dilated_density");
    for (i, key) in keys.iter().enumerate() {
        let [x, y, z] = key.index();
        println!("{},{x},{y},{z},{:.17e},{:.17e},{:.17e},{:.17e}", key.level(), last.rho[i],
            evaluated.eroded.projected_rho[i], evaluated.nominal_rho[i], evaluated.dilated_rho[i]);
    }
    println!("node_x,node_y,node_z,eroded_uy_x,eroded_uy_y,eroded_uy_z,eroded_uz_x,eroded_uz_y,eroded_uz_z");
    for (i, p) in nodes.iter().enumerate() {
        let uy = &evaluated.eroded.objective.displacements[0];
        let uz = &evaluated.eroded.objective.displacements[1];
        println!("{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            p[0], p[1], p[2], uy[3*i], uy[3*i+1], uy[3*i+2], uz[3*i], uz[3*i+1], uz[3*i+2]);
    }
    // Publish the aligned accepted prefix before any interruption or later
    // diagnostic build/solve can fail. Never assess a failed transition model.
    if report.stop != RobustContinuationStop3::ScheduleComplete {
        return Err("robust continuation stopped early; exported fields describe the retained feasible prefix".into());
    }
    let fine_tree = tree.refined(&tree.leaves().iter().copied().collect::<Vec<_>>(),
        || ControlFlow::Continue(()))?;
    let fine = AdaptiveElasticity3::build(domain, &fine_tree, &CurvedCantilever,
        &material, &|p| p[0] == 0.0, ElasticityOptions3::default(), &mut quadrature);
    eprintln!("diagnostic_geometry_work={:?}", quadrature.work());
    let laws = [GoalReferenceLoad3 { load: ReferenceLoad3::body(&y), weight: 0.3 },
        GoalReferenceLoad3 { load: ReferenceLoad3::body(&z), weight: 0.7 }];
    let assessment = study.estimate_robust_reference_compliance_enrichment(fine?, &laws,
        &last.rho, report.params, robustness, GoalRefinementOptions3 {
            preconditioner: GoalPreconditioner3::Jacobi { max_contributions: 100_000_000 },
            ..Default::default() }, &mut control);
    eprintln!("total_solve_work={:?}", control.work());
    let assessment = assessment?;
    let estimate = &assessment.estimate;
    println!("eroded_coarse_compliance,eroded_fine_compliance,two_grid_correction");
    println!("{:.17e},{:.17e},{:.17e}", estimate.coarse_value, estimate.fine_value, estimate.correction);
    let marks = estimate.mark(0.5, 2, || ControlFlow::Continue(()))?;
    eprintln!("probe_material=inherited-eroded-stiffness; fine_raw_refiltered=false; marking_target=0.5; marking_cap=2; marking_fraction={:.17e}; marking_target_met={}; fine_reoptimization_performed=false; continuum_error_bound=false; dilated_volume_error_estimated=false",
        marks.achieved_fraction, marks.target_met);
    println!("marked_level,marked_x,marked_y,marked_z,eroded_goal_mass");
    for key in &marks.marked {
        let [x, y, z] = key.index();
        println!("{},{x},{y},{z},{:.17e}", key.level(), estimate.marking_mass[key]);
    }
    Ok(())
}
