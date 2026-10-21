//! Minimize dilated material under an authored simultaneous-load requirement.
#[path = "material/binary.rs"]
mod binary;
use super::*;
use fs_topopt::sdf3::robust::envelope::load_space::design::{LoadSpaceDesignPoint3,
    limits::{LoadSpaceCap3, LoadSpaceLimitedDesign3, LoadSpaceLimitedOptions3,
        start::{LoadSpaceStartOptions3, LoadSpaceStartStop3, controlled_load_space_feasible_start}}};

pub fn run(args: &[String]) -> Result<(), Box<dyn std::error::Error>> { run_mode(args, false) }
pub fn run_binary(args: &[String]) -> Result<(), Box<dyn std::error::Error>> { run_mode(args, true) }
fn run_mode(args: &[String], binary_mode: bool) -> Result<(), Box<dyn std::error::Error>> {
    let maximum = if binary_mode { 4 } else { 3 };
    if !(1..=maximum).contains(&args.len()) {
        let mode = if binary_mode { "--binary-material" } else { "--minimum-material" };
        return Err(format!("usage: worst_load_sdf3 --load-space {mode} CAP [UPDATES [TOTAL_KRYLOV_ITERATIONS{}]]",
            if binary_mode { " [BINARY_DILATED_VOLUME_LIMIT]" } else { "" }).into());
    }
    // Authored physical requirement, never a baseline-derived normalization.
    let cap = LoadSpaceCap3::new(2, args[0].parse()?, 0.03)?;
    let updates = args.get(1).map_or(Ok(200), |s| s.parse::<usize>())?;
    let iterations = args.get(2).map_or(Ok(250_000), |s| s.parse::<usize>())?;
    let binary_volume_limit = args.get(3).map_or(Ok(1.0), |s| s.parse::<f64>())?;
    if !binary_volume_limit.is_finite() || !(0.0..=1.0).contains(&binary_volume_limit) {
        return Err("binary dilated-volume limit must be finite and in [0,1]".into());
    }
    let mut poll = |_| ControlFlow::Continue(());
    let mut geometry = QuadratureControl3::new(QuadratureOptions3::default(), &mut poll)?;
    let op = CutElasticity3::build(HexCell::try_new([0.0;3], [1.0;3])?, [3,2,2], &Slab,
        &IsotropicElastic::new(1.0,0.3,1.0)?, &|p| p[0] == 0.0, ElasticityOptions3::default(), &mut geometry)?;
    let y = op.body_load(&|p| [0.0,-1.0-p[1],0.0], || ControlFlow::Continue(()))?;
    let second = op.body_load(&|p| [0.0,-0.6*(1.0+p[1]),-1.0-p[0]], || ControlFlow::Continue(()))?;
    let keys = op.cell_keys(); let nodes = op.nodes().to_vec();
    let params = SimpParams { penal: 2.0, beta: 2.0, e_min: 1e-3, ..Default::default() };
    let mut study = CutDensityStudy3::new(op, 0.15, params);
    let raw = vec![0.3; study.cells()]; let basis = [&y[..], &second[..]];
    let options = LoadSpaceLimitedOptions3::default(); let start_options = LoadSpaceStartOptions3::default();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: iterations, ..Default::default() }, &mut poll);
    eprintln!("units=dimensionless; load_set=Fq; coefficient_norm_limit=1; load_offset=none; authored_cap={}; spectral_temperature={}; smoothing_bias_bound={}; params={params:?}; options={options:?}; start_options={start_options:?}; total_krylov_budget={iterations}; geometry_work={:?}",
        cap.cap(), cap.policy().temperature(), cap.policy().bias_bound(), geometry.work());
    let start = controlled_load_space_feasible_start(&mut study, &basis, &raw, &cap, options, start_options, &mut control)?;
    eprintln!("restoration_stop={:?}; restoration_error={:?}; blend={:?}; evaluations={}; cumulative_work={:?}",
        start.stop, start.error, start.blend, start.evaluations, start.work);
    let Some(seed) = &start.point else {
        return Err("no checked combined-load feasible seed; full-density refusal is not a global infeasibility proof".into());
    };
    if !matches!(start.stop, LoadSpaceStartStop3::IncomingSatisfied | LoadSpaceStartStop3::Restored) {
        export(seed, &cap, &keys, &nodes, "interrupted-restoration-seed");
        return Err("restoration stopped; only its checked feasible prefix was exported".into());
    }
    let mut session = match LoadSpaceLimitedDesign3::new(&mut study, &basis, &seed.rho, cap.clone(), options, &mut control) {
        Ok(session) => session,
        Err(error) => { export(seed, &cap, &keys, &nodes, "restored-seed"); return Err(error.into()); }
    };
    let outcome = session.run(updates);
    println!("iteration,dilated_volume,smooth_normalized_maximum,actual_normalized_maximum,combined_load_feasible");
    for row in session.history() {
        println!("{},{:.17e},{:.17e},{:.17e},{}", row.iteration, row.dilated_volume,
            row.smooth_maximum, row.maximum, row.feasible);
    }
    eprintln!("optimizer_outcome={outcome:?}; gradient_passed={}; physical_work={:?}; optimizer_work={:?}; terminal_kkt_is_incumbent_certificate=false",
        session.gradient_check().numerical.passed(), session.work(), session.optimizer_work());
    let Some(best) = session.best_feasible() else {
        export(seed, &cap, &keys, &nodes, "restored-seed");
        return Err("optimizer has no feasible incumbent; exported restoration seed is not optimized material".into());
    };
    // Publish solved fields before a final solve or interrupted run can fail.
    export(best, &cap, &keys, &nodes, "least-volume-accepted-feasible");
    let terminal = outcome?;
    let checked = session.recheck_best();
    eprintln!("fresh_full_load_set_check_passed={}; final_physical_work={:?}", matches!(&checked, Ok(Some(_))), session.work());
    let checked = checked?.ok_or("no retained design available for combined-load recheck")?;
    if !matches!(terminal.stop, ProjectedAlStop::Converged | ProjectedAlStop::IterationLimit) {
        return Err("optimizer stopped before the requested allowance; exported result is the retained feasible prefix".into());
    }
    if binary_mode {
        drop(session);
        return binary::finish(&mut study, &basis, &checked.rho, &cap, options,
            binary_volume_limit, &keys, &nodes, &mut control);
    }
    Ok(())
}

fn export(point: &LoadSpaceDesignPoint3, cap: &LoadSpaceCap3, keys: &[[usize;3]], nodes: &[[f64;3]], source: &str) {
    let e = &point.evaluation;
    eprintln!("exported_source={source}; worst_physical_compliance={:.17e}; physical_cap={:.17e}; physical_cap_satisfied={}; normalized_spectrum={:?}; eigenpair_residual={:.3e}; reciprocity_defect={:.3e}; fixed_grid=true; continuum_certificate=false; minimum_feature_size=false; global_optimum=false",
        e.worst_compliance, cap.cap(), e.worst_compliance <= cap.cap(), e.spectrum.eigenvalues,
        e.spectrum.eigenpair_residual, e.spectrum.reciprocity_defect);
    println!("basis_index,maximizing_coefficient,basis_compliance_diagnostic");
    for (i, (&q, &c)) in e.worst_coefficients.iter().zip(&e.physical.eroded.objective.case_compliances).enumerate() {
        println!("{i},{q:.17e},{c:.17e}");
    }
    println!("cell_x,cell_y,cell_z,raw_density,eroded_density,nominal_density,dilated_density");
    for (i,k) in keys.iter().enumerate() {
        println!("{},{},{},{:.17e},{:.17e},{:.17e},{:.17e}", k[0],k[1],k[2],point.rho[i],
            e.physical.eroded.projected_rho[i],e.physical.nominal_rho[i],e.physical.dilated_rho[i]);
    }
    println!("node_x,node_y,node_z,worst_fx,worst_fy,worst_fz,worst_ux,worst_uy,worst_uz");
    for (i,p) in nodes.iter().enumerate() {
        let f = &e.worst_force; let u = &e.worst_displacement;
        println!("{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            p[0],p[1],p[2],f[3*i],f[3*i+1],f[3*i+2],u[3*i],u[3*i+1],u[3*i+2]);
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn authored_cap_and_iteration_policy_are_checked_before_geometry() {
        for args in [vec![], vec!["NaN"], vec!["0"], vec!["1000", "bad-count"], vec!["1000", "0", "1", "extra"]] {
            let args: Vec<String> = args.into_iter().map(String::from).collect();
            assert!(super::run(&args).is_err());
        }
    }
    #[test]
    fn minimum_material_load_space_command_uses_complete_physics() {
        // Loose authored cap isolates integration, not optimization quality.
        let args = ["--minimum-material", "100000000", "0", "250000"].map(|s| s.to_string());
        super::super::run(&args).unwrap();
    }
}
