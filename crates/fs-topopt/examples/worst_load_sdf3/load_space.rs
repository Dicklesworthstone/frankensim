//! Optimize simultaneous combinations, not a finite list of exclusive cases.
//! --binary-material CAP [UPDATES [ITERATIONS [BINARY_DILATED_VOLUME_LIMIT]]]
//! adds full-load-set solid/void selection after minimum-material optimization.
#[path = "load_space/material.rs"]
mod material;
#[path = "load_space/refinement.rs"]
mod refinement;
use super::*;
use fs_topopt::sdf3::robust::envelope::load_space::{LoadSpace3, design::LoadSpaceDesignStudy3};

pub fn run(args: &[String]) -> Result<(), Box<dyn std::error::Error>> {
    if args.first().is_some_and(|s| s == "--minimum-material") { return material::run(&args[1..]); }
    if args.first().is_some_and(|s| s == "--binary-material") { return material::run_binary(&args[1..]); }
    if args.first().is_some_and(|s| s == "--refine-material") { return refinement::run(&args[1..]); }
    if args.len() > 2 { return Err("usage: worst_load_sdf3 --load-space [UPDATES [TOTAL_KRYLOV_ITERATIONS]]".into()); }
    let updates = args.first().map_or(Ok(200), |s| s.parse::<usize>())?;
    let iterations = args.get(1).map_or(Ok(250_000), |s| s.parse::<usize>())?;
    let mut poll = |_| ControlFlow::Continue(());
    let mut geometry = QuadratureControl3::new(QuadratureOptions3::default(), &mut poll)?;
    let op = CutElasticity3::build(HexCell::try_new([0.0;3], [1.0;3])?, [3,2,2], &Slab,
        &IsotropicElastic::new(1.0,0.3,1.0)?, &|p| p[0] == 0.0, ElasticityOptions3::default(), &mut geometry)?;
    let y = op.body_load(&|p| [0.0,-1.0-p[1],0.0], || ControlFlow::Continue(()))?;
    // Correlation and relative amplitudes are part of the authored LOAD SET.
    let second = op.body_load(&|p| [0.0,-0.6*(1.0+p[1]),-1.0-p[0]], || ControlFlow::Continue(()))?;
    let keys = op.cell_keys(); let nodes = op.nodes().to_vec();
    let params = SimpParams { penal: 2.0, beta: 2.0, e_min: 1e-3, ..Default::default() };
    let mut study = CutDensityStudy3::new(op, 0.15, params);
    let raw = vec![0.3; study.cells()]; let basis = [&y[..], &second[..]];
    let options = EnvelopeDesignOptions3::default();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: iterations, ..Default::default() }, &mut poll);
    let baseline = study.evaluate_robust_load_space(&raw, &basis, options.robustness,
        &LoadSpace3::new(2,1.0,0.05)?, &mut control)?;
    // A single charged baseline fixes the reference for the WHOLE set, once.
    let policy = LoadSpace3::new(2, baseline.worst_compliance, 0.05)?;
    let mut design = LoadSpaceDesignStudy3::new(&mut study, &basis, &raw, policy.clone(), options, &mut control)?;
    let outcome = design.run(updates);
    println!("iteration,smooth_spectral_maximum,worst_normalized_compliance,dilated_volume,feasible");
    for row in design.history() {
        println!("{},{:.17e},{:.17e},{:.17e},{}", row.iteration, row.smooth_maximum,
            row.maximum, row.dilated_volume, row.feasible);
    }
    let recovery = if outcome.is_ok() { Some(design.recover_feasible()) } else { None };
    let incumbent = design.best_feasible();
    let recovered = recovery.as_ref().and_then(|r| r.as_ref().ok()).filter(|r|
        incumbent.is_none_or(|b| r.evaluation.spectrum.smooth_maximum < b.evaluation.spectrum.smooth_maximum));
    let feasible = recovered.is_some() || incumbent.is_some();
    let point = recovered.or(incumbent).unwrap_or_else(|| design.accepted());
    let e = &point.evaluation;
    let source = if recovered.is_some() { "recovered-terminal" } else if incumbent.is_some() { "accepted-feasible" } else { "infeasible-accepted" };
    eprintln!("load_set=f0*q0+f1*q1; coefficient_norm_limit=1; nominal_load_offset=none; reference={}; temperature={}; algebraic_bias_bound={}; params={params:?}; design_policy={options:?}; exported_source={source}; feasible={feasible}; outcome={outcome:?}",
        policy.reference(), policy.temperature(), policy.bias_bound());
    eprintln!("worst_physical_compliance={:.17e}; normalized_spectrum={:?}; reciprocity_defect={:.3e}; eigenpair_residual={:.3e}; orthogonality_defect={:.3e}; gradient_passed={}; physical_work={:?}; optimizer_work={:?}; geometry_work={:?}; probability_model=false; continuum_certificate=false",
        e.worst_compliance, e.spectrum.eigenvalues, e.spectrum.reciprocity_defect,
        e.spectrum.eigenpair_residual, e.spectrum.orthogonality_defect,
        design.gradient_check().numerical.passed(), design.work(), design.optimizer_work(), geometry.work());
    println!("basis_index,worst_combination_coefficient,basis_compliance");
    for (i, (&q, &c)) in e.worst_coefficients.iter().zip(&e.physical.eroded.objective.case_compliances).enumerate() {
        println!("{i},{q:.17e},{c:.17e}");
    }
    println!("cell_x,cell_y,cell_z,raw_density,eroded_density,nominal_density,dilated_density");
    for (i,k) in keys.iter().enumerate() {
        println!("{},{},{},{:.17e},{:.17e},{:.17e},{:.17e}", k[0], k[1], k[2], point.rho[i],
            e.physical.eroded.projected_rho[i], e.physical.nominal_rho[i], e.physical.dilated_rho[i]);
    }
    println!("node_x,node_y,node_z,worst_fx,worst_fy,worst_fz,worst_ux,worst_uy,worst_uz");
    for (i,p) in nodes.iter().enumerate() {
        let f = &e.worst_force; let u = &e.worst_displacement;
        println!("{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            p[0],p[1],p[2],f[3*i],f[3*i+1],f[3*i+2],u[3*i],u[3*i+1],u[3*i+2]);
    }
    let terminal = outcome?;
    if let Some(recovery) = recovery { recovery?; }
    if !feasible { return Err("no material-feasible load-space design; exported accepted point is infeasible".into()); }
    if !matches!(terminal.stop, ProjectedAlStop::Converged | ProjectedAlStop::IterationLimit) {
        return Err("load-space optimization stopped before its allowance; exported fields retain the solved prefix".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn command_runs_real_load_space_physics_and_refuses_bad_inputs() {
        let args = ["0".to_string(), "250000".to_string()];
        super::run(&args).unwrap();
        assert!(super::run(&["not-a-number".to_string()]).is_err());
        assert!(super::run(&["0".into(),"250000".into(),"extra".into()]).is_err());
    }
}
