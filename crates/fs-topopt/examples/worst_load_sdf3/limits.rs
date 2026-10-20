//! Minimum-material mode: physical caps are authored, never inferred from a baseline.
use super::*;
use fs_topopt::sdf3::robust::envelope::limits::{ComplianceCaps3, ComplianceLimitedDesign3,
    ComplianceLimitedOptions3, start::{ComplianceStartOptions3, ComplianceStartStop3,
    controlled_compliance_feasible_start}};
use fs_topopt::sdf3::robust::envelope::optimize::EnvelopeDesignPoint3;

pub fn run(args: &[String]) -> Result<(), Box<dyn std::error::Error>> { run_mode(args, false) }
pub fn run_binary(args: &[String]) -> Result<(), Box<dyn std::error::Error>> { run_mode(args, true) }
fn run_mode(args: &[String], binary: bool) -> Result<(), Box<dyn std::error::Error>> {
    if !(2..=4).contains(&args.len()) {
        let mode = if binary { "--binary-material" } else { "--minimum-material" };
        return Err(format!("usage: worst_load_sdf3 {mode} CAP_Y CAP_Z [ACCEPTED_UPDATES [TOTAL_KRYLOV_ITERATIONS]]").into());
    }
    // Admit authored requirements before constructing geometry or load fields.
    let caps = ComplianceCaps3::new(vec![args[0].parse()?, args[1].parse()?], 0.03, 2)?;
    let updates = args.get(2).map_or(Ok(200), |s| s.parse::<usize>())?;
    let iterations = args.get(3).map_or(Ok(250_000), |s| s.parse::<usize>())?;
    let options = ComplianceLimitedOptions3::default();
    let start_options = ComplianceStartOptions3::default();
    let mut poll = |_| ControlFlow::Continue(());
    let mut geometry = QuadratureControl3::new(QuadratureOptions3::default(), &mut poll)?;
    let op = CutElasticity3::build(HexCell::try_new([0.0; 3], [1.0; 3])?, [3, 2, 2], &Slab,
        &IsotropicElastic::new(1.0, 0.3, 1.0)?, &|p| p[0] == 0.0, ElasticityOptions3::default(), &mut geometry)?;
    let y = op.body_load(&|p| [0.0, -1.0-p[1], 0.0], || ControlFlow::Continue(()))?;
    let z = op.body_load(&|p| [0.0, 0.0, -1.0-p[0]], || ControlFlow::Continue(()))?;
    let keys = op.cell_keys(); let nodes = op.nodes().to_vec();
    let params = SimpParams { penal: 2.0, beta: 2.0, e_min: 1e-3, ..Default::default() };
    let mut study = CutDensityStudy3::new(op, 0.15, params);
    let incoming = vec![0.3; study.cells()]; let forces = [&y[..], &z[..]];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: iterations, ..Default::default() }, &mut poll);
    eprintln!("units=dimensionless; seeds=none; mode=minimum-dilated-material; params={params:?}; options={options:?}; start_policy={start_options:?}; caps={:?}; temperature={}; smoothing_bias_bound={:.17e}; solve_budget={iterations}; geometry_work={:?}",
        caps.caps(), caps.envelope().temperature(), caps.envelope().bias_bound(), geometry.work());
    let start = controlled_compliance_feasible_start(&mut study, &forces, &incoming, &caps,
        options, start_options, &mut control)?;
    eprintln!("start_stop={:?}; start_error={:?}; blend={:?}; start_evaluations={}; full_density_envelope={:?}; full_density_maximum_ratio={:?}; solve_work={:?}",
        start.stop, start.error, start.blend, start.evaluations,
        start.full_density_envelope, start.full_density_maximum_ratio, start.work);
    let Some(seed) = &start.point else {
        return Err("no cap-feasible seed was found; a refused full endpoint is not a global infeasibility proof".into());
    };
    if !matches!(start.stop, ComplianceStartStop3::IncomingSatisfied | ComplianceStartStop3::Restored) {
        export(seed, &caps, &keys, &nodes, "interrupted-restoration-seed");
        return Err("restoration stopped; exported point is the retained solved feasible endpoint".into());
    }
    let mut design = match ComplianceLimitedDesign3::new(&mut study, &forces, &seed.rho,
        caps.clone(), options, &mut control) {
        Ok(design) => design,
        Err(error) => {
            export(seed, &caps, &keys, &nodes, "restored-seed-only");
            return Err(error.into());
        }
    };
    let outcome = design.run(updates);
    println!("iteration,dilated_volume,envelope,maximum_ratio,limiting_case,all_caps_feasible");
    for row in design.history() {
        println!("{},{:.17e},{:.17e},{:.17e},{},{}", row.iteration, row.dilated_volume,
            row.envelope, row.maximum_ratio, row.limiting_case, row.feasible);
    }
    eprintln!("gradient_passed={}; terminal_optimizer={outcome:?}; optimizer_work={:?}; solve_work={:?}",
        design.gradient_check().numerical.passed(), design.optimizer_work(), design.work());
    let Some(best) = design.best_feasible() else {
        export(seed, &caps, &keys, &nodes, "restored-seed-only");
        return Err("optimizer retained no cap-feasible point; exporting the original restoration seed".into());
    };
    // Emit the complete previously solved prefix BEFORE any late error/recheck.
    export(best, &caps, &keys, &nodes, "minimum-material-feasible-incumbent");
    let terminal = outcome?;
    let recheck = design.recheck_best();
    eprintln!("fresh_resolve_passed={}; final_solve_work={:?}; terminal_kkt_is_incumbent_certificate=false",
        matches!(&recheck, Ok(Some(_))), design.work());
    let checked = recheck?.ok_or("no incumbent available for fresh cap verification")?;
    if !matches!(terminal.stop, ProjectedAlStop::Converged | ProjectedAlStop::IterationLimit) {
        return Err("optimizer stopped before its requested allowance; the exported feasible prefix remains available".into());
    }
    if binary {
        drop(design);
        return super::binary::finish(&mut study, &forces, &checked.rho, &caps,
            options, &keys, &nodes, &mut control);
    }
    Ok(())
}

fn export(point: &EnvelopeDesignPoint3, caps: &ComplianceCaps3,
    keys: &[[usize; 3]], nodes: &[[f64; 3]], label: &str) {
    let e = &point.evaluation;
    eprintln!("exported_design={label}; each_physical_cap_satisfied={}; continuum_load_guarantee=false; global_optimum=false; geometry_rebuilt=false; durable_resume=false",
        caps.satisfied_by(&e.physical.eroded.objective.case_compliances));
    println!("load_case,eroded_compliance,physical_cap,ratio,cap_satisfied");
    for (i, (&c, &limit)) in e.physical.eroded.objective.case_compliances.iter().zip(caps.caps()).enumerate() {
        println!("{i},{c:.17e},{limit:.17e},{:.17e},{}", c/limit, c <= limit);
    }
    println!("cell_x,cell_y,cell_z,raw_density,eroded_density,nominal_density,dilated_density");
    for (i, key) in keys.iter().enumerate() {
        println!("{},{},{},{:.17e},{:.17e},{:.17e},{:.17e}", key[0], key[1], key[2], point.rho[i],
            e.physical.eroded.projected_rho[i], e.physical.nominal_rho[i], e.physical.dilated_rho[i]);
    }
    println!("node_x,node_y,node_z,eroded_uy_x,eroded_uy_y,eroded_uy_z,eroded_uz_x,eroded_uz_y,eroded_uz_z");
    for (i, p) in nodes.iter().enumerate() {
        let uy = &e.physical.eroded.objective.displacements[0]; let uz = &e.physical.eroded.objective.displacements[1];
        println!("{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            p[0], p[1], p[2], uy[3*i], uy[3*i+1], uy[3*i+2], uz[3*i], uz[3*i+1], uz[3*i+2]);
    }
}
