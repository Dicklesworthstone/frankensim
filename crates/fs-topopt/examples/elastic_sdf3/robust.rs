//! Robust mode of the same raw-implicit curved cantilever example.
use super::*;
use fs_topopt::sdf3::robust::{RobustDesignOptions3, RobustDesignStudy3, RobustOptions3};

pub fn run(updates: usize, iterations: usize) -> Result<(), Box<dyn std::error::Error>> {
    let mut poll = |_| ControlFlow::Continue(());
    let mut quadrature = QuadratureControl3::new(QuadratureOptions3::default(), &mut poll)?;
    let op = CutElasticity3::build(HexCell::try_new([0.0; 3], [1.0; 3])?, [4, 2, 2], &CurvedCantilever,
        &IsotropicElastic::new(1.0, 0.3, 1.0)?, &|p| p[0] == 0.0, ElasticityOptions3::default(), &mut quadrature)?;
    let y = op.body_load(&|_| [0.0, -1.0, 0.0], || ControlFlow::Continue(()))?;
    let z = op.body_load(&|_| [0.0, 0.0, -1.0], || ControlFlow::Continue(()))?;
    let keys = op.cell_keys(); let nodes = op.nodes().to_vec();
    let params = SimpParams { beta: 3.0, ..Default::default() };
    let mut study = CutDensityStudy3::new(op, 0.15, params);
    let rho = vec![0.3; study.cells()];
    let loads = [LoadCase { force: &y, weight: 0.3 }, LoadCase { force: &z, weight: 0.7 }];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: iterations, ..Default::default() }, &mut poll);
    let robustness = RobustOptions3::default();
    // The named example explicitly fixes its conditioning scale from a real
    // baseline; that extra solve is charged to the same cumulative work limit.
    let baseline = study.evaluate_robust(&rho, &loads, robustness, &mut control)?;
    let options = RobustDesignOptions3 { robustness, compliance_scale: baseline.eroded.objective.compliance,
        volume_limit: 0.5, ..Default::default() };
    let mut design = RobustDesignStudy3::new(&mut study, &loads, &rho, options, &mut control)?;
    let result = design.run(updates);
    println!("iteration,eroded_compliance,eroded_volume,nominal_volume,dilated_volume,constraint_violation,feasible");
    for row in design.history() {
        println!("{},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{}", row.iteration, row.eroded_compliance,
            row.eroded_volume_fraction, row.nominal_volume_fraction, row.dilated_volume_fraction,
            row.constraint_violation, row.feasible);
    }
    eprintln!("units=dimensionless; seeds=none; p={}; beta={}; eta={}; threshold_offset={}; filter_radius=0.15; density_floor={}; volume_limit={}; compliance_scale={:.17e}; solve_budget={iterations}; max_evaluations={}; tolerance={}",
        params.penal, params.beta, params.eta, robustness.threshold_offset, options.density_floor,
        options.volume_limit, options.compliance_scale, options.optimizer.max_evaluations, options.optimizer.tolerance);
    eprintln!("optimizer_result={result:?}; accepted_feasible={}; optimizer_work={:?}; solve_work={:?}",
        design.feasible(), design.optimizer_work(), design.work());
    let feasible = design.best_feasible().is_some();
    let retained = design.best_feasible().unwrap_or_else(|| design.accepted());
    let e = &retained.evaluation;
    eprintln!("exported_design={}; exported_feasible={feasible}; geometry_rebuilt=false; geometric_offset_certified=false; minimum_feature_size_certified=false; continuum_safety_certified=false; global_optimum_certified=false",
        if feasible { "best-feasible" } else { "last-accepted" });
    println!("cell_x,cell_y,cell_z,raw_density,eroded_density,nominal_density,dilated_density");
    for (i, key) in keys.iter().enumerate() {
        println!("{},{},{},{:.17e},{:.17e},{:.17e},{:.17e}", key[0], key[1], key[2], retained.rho[i],
            e.eroded.projected_rho[i], e.nominal_rho[i], e.dilated_rho[i]);
    }
    println!("node_x,node_y,node_z,eroded_uy_x,eroded_uy_y,eroded_uy_z,eroded_uz_x,eroded_uz_y,eroded_uz_z");
    for (i, p) in nodes.iter().enumerate() {
        let uy = &e.eroded.objective.displacements[0]; let uz = &e.eroded.objective.displacements[1];
        println!("{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            p[0], p[1], p[2], uy[3*i], uy[3*i+1], uy[3*i+2], uz[3*i], uz[3*i+1], uz[3*i+2]);
    }
    result?;
    if !feasible { return Err("no accepted design meets the dilated-volume tolerance; exported point is infeasible".into()); }
    Ok(())
}
