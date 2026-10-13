//! Optimize, assess every cap, selectively refit, then probe the new design.
use super::*;
use fs_cutfem::elastic3::{adaptive::AdaptiveElasticity3, surface::ReferenceLoad3};
use fs_cutfem::octree3::{Octant3, Octree3};
use fs_topopt::sdf3::robust::envelope::limits::{ComplianceCaps3, ComplianceLimitedDesign3,
    refit::ComplianceCapRefitOptions3, refinement::ComplianceCapRefinement3,
    start::{ComplianceStartStop3, controlled_compliance_feasible_start}};
use fs_topopt::sdf3::robust::envelope::optimize::EnvelopeDesignPoint3;
use fs_topopt::sdf3_goal::{GoalPreconditioner3, GoalRefinementOptions3};

type Error = Box<dyn std::error::Error>;
const CHANGE_TOLERANCE: f64 = 0.02;
fn parse(args: &[String]) -> Result<(ComplianceCaps3, usize, usize), Error> {
    if !(2..=4).contains(&args.len()) {
        return Err("usage: worst_load_sdf3 --refine-material CAP_Y CAP_Z [UPDATES_PER_GRID [TOTAL_KRYLOV_ITERATIONS]]".into());
    }
    let caps = ComplianceCaps3::new(vec![args[0].parse()?, args[1].parse()?], 0.03, 2)?;
    Ok((caps, args.get(2).map_or(Ok(200), |s| s.parse::<usize>())?,
        args.get(3).map_or(Ok(250_000), |s| s.parse::<usize>())?))
}
fn y(p: [f64; 3]) -> [f64; 3] { [0.0, -1.0-p[1], 0.0] }
fn z(p: [f64; 3]) -> [f64; 3] { [0.0, 0.0, -1.0-p[0]] }
fn laws() -> [ReferenceLoad3<'static>; 2] { [ReferenceLoad3::body(&y), ReferenceLoad3::body(&z)] }
fn build(tree: &Octree3, geometry: &mut QuadratureControl3<'_>) -> Result<AdaptiveElasticity3, Error> {
    Ok(AdaptiveElasticity3::build(HexCell::try_new([0.0; 3], [1.0; 3])?, tree, &Slab,
        &IsotropicElastic::new(1.0, 0.3, 1.0)?, &|p| p[0] == 0.0,
        ElasticityOptions3::default(), geometry)?)
}

#[allow(clippy::too_many_arguments)]
fn assess(study: &mut CutDensityStudy3<AdaptiveElasticity3>, tree: &Octree3,
    point: &EnvelopeDesignPoint3, caps: &ComplianceCaps3, options: ComplianceCapRefitOptions3,
    control: &mut SolveControl<'_>, geometry: &mut QuadratureControl3<'_>, grid: usize,
) -> Result<ComplianceCapRefinement3, Error> {
    let probe_tree = tree.refined(study.operator().leaves(), || ControlFlow::Continue(()))?;
    let operator = build(&probe_tree, geometry);
    eprintln!("assessed_grid={grid}; cumulative_geometry_work={:?}", geometry.work());
    let params = study.params();
    let result = study.estimate_compliance_cap_enrichment(operator?, &laws(), &point.rho,
        params, options.design.robustness, caps, GoalRefinementOptions3 {
            preconditioner: GoalPreconditioner3::Jacobi { max_contributions: 100_000_000 },
            max_transfer_terms: options.max_transfer_terms, max_load_cases: options.max_load_cases,
            ..Default::default()
        }, control);
    eprintln!("assessed_grid={grid}; cumulative_solve_work={:?}", control.work());
    let checked = result?;
    println!("grid,case,cap,coarse_compliance,fine_compliance,cap_normalized_change,both_caps_satisfied");
    for (i, check) in checked.checks.iter().enumerate() {
        println!("{grid},{i},{:.17e},{:.17e},{:.17e},{:.17e},{}", check.cap,
            check.coarse, check.fine, check.normalized_change, check.satisfied());
    }
    eprintln!("assessed_grid={grid}; coarse_envelope={:.17e}; fine_envelope={:.17e}; two_grid_pass={}; continuum_bound=false; dilated_volume_error_estimated=false",
        checked.coarse_envelope, checked.fine_envelope,
        checked.stable_and_feasible(options.design.optimizer.tolerance, CHANGE_TOLERANCE));
    Ok(checked)
}

pub fn run(args: &[String]) -> Result<(), Error> {
    // Admission precedes geometry; caps are never inferred from solved fields.
    let (caps, updates, iterations) = parse(args)?;
    let options = ComplianceCapRefitOptions3 { accepted_steps: updates, ..Default::default() };
    let params = SimpParams { penal: 2.0, beta: 2.0, e_min: 1e-3, ..Default::default() };
    let tree = Octree3::uniform(1, 4, 4096)?;
    let geometry_limits = QuadratureOptions3::default();
    let mut poll = |_| ControlFlow::Continue(());
    let mut geometry = QuadratureControl3::new(geometry_limits, &mut poll)?;
    let operator = build(&tree, &mut geometry)?;
    let keys = operator.leaves().to_vec(); let nodes = operator.nodes().to_vec();
    let mut study = CutDensityStudy3::new(operator, 0.15, params);
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: iterations, ..Default::default() }, &mut poll);
    eprintln!("units=dimensionless; seeds=none; params={params:?}; filter_radius=0.15; cap_policy={caps:?}; refit_policy={options:?}; total_krylov_budget={iterations}; geometry_limits={geometry_limits:?}; max_refits=1; marking_fraction=0.5; max_marks=2; cap_change_tolerance={CHANGE_TOLERANCE}; cross_grid_descent_claimed=false");
    let forces: Vec<_> = laws().iter().map(|&law|
        study.operator().reference_load(law, || ControlFlow::Continue(()))).collect::<Result<_, _>>()?;
    let borrowed: Vec<_> = forces.iter().map(Vec::as_slice).collect();
    let start = controlled_compliance_feasible_start(&mut study, &borrowed, &vec![0.3; keys.len()],
        &caps, options.design, options.start, &mut control)?;
    eprintln!("grid=0; restoration_stop={:?}; restoration_error={:?}; blend={:?}; work={:?}",
        start.stop, start.error, start.blend, start.work);
    let seed = start.point.as_ref().ok_or("no cap-feasible coarse seed; full-endpoint rejection is not a global infeasibility proof")?;
    export(seed, &caps, &keys, &nodes, 0, "restoration-seed");
    if !matches!(start.stop, ComplianceStartStop3::IncomingSatisfied | ComplianceStartStop3::Restored) {
        return Err("coarse restoration interrupted; only the solved seed is exported".into());
    }
    let mut design = ComplianceLimitedDesign3::new(&mut study, &borrowed, &seed.rho,
        caps.clone(), options.design, &mut control)?;
    let outcome = design.run(updates);
    println!("grid,iteration,dilated_volume,envelope,maximum_ratio,each_cap_feasible");
    for row in design.history() {
        println!("0,{},{:.17e},{:.17e},{:.17e},{}", row.iteration, row.dilated_volume,
            row.envelope, row.maximum_ratio, row.feasible);
    }
    let best = design.best_feasible().cloned();
    eprintln!("grid=0; optimizer_outcome={outcome:?}; gradient_passed={}; work={:?}",
        design.gradient_check().numerical.passed(), design.work());
    if let Some(best) = &best { export(best, &caps, &keys, &nodes, 0, "optimized-incumbent"); }
    drop(design);
    let terminal = outcome?;
    if !matches!(terminal.stop, ProjectedAlStop::Converged | ProjectedAlStop::IterationLimit) {
        return Err("coarse optimizer stopped early; the solved prefix remains exported".into());
    }
    let best = best.ok_or("coarse optimizer retained no cap-feasible design")?;
    let check = assess(&mut study, &tree, &best, &caps, options, &mut control, &mut geometry, 0)?;
    if check.stable_and_feasible(options.design.optimizer.tolerance, CHANGE_TOLERANCE) {
        eprintln!("refit_performed=false; current_two_grid_cap_check_passed=true");
        return Ok(());
    }
    let marks = check.mark(0.5, 2, || ControlFlow::Continue(()))?;
    eprintln!("source_grid=0; marks={marks:?}");
    if marks.marked.is_empty() { return Err("no refinement signal; cap accuracy is not established".into()); }
    // Refine the SOURCE background, not the already enriched diagnostic grid.
    let target_tree = tree.refined(&marks.marked, || ControlFlow::Continue(()))?;
    let target = build(&target_tree, &mut geometry);
    eprintln!("candidate_geometry_work={:?}", geometry.work());
    let refitted = study.refit_compliance_caps(target?, &best.rho, params, &laws(), &caps, options, &mut control);
    eprintln!("candidate_solve_work={:?}", control.work());
    let mut refitted = refitted?;
    let op = refitted.study.operator();
    eprintln!("grid=1; restoration_stop={:?}; restoration_error={:?}; optimizer_outcome={:?}; fresh_recheck={:?}",
        refitted.restoration.stop, refitted.restoration.error, refitted.outcome,
        refitted.recheck.as_ref().map(|r| r.as_ref().map(|_| ())));
    if let Some(seed) = &refitted.restoration.point {
        export(seed, &caps, op.leaves(), op.nodes(), 1, "restoration-seed");
    }
    println!("grid,iteration,dilated_volume,envelope,maximum_ratio,each_cap_feasible");
    for row in &refitted.history {
        println!("1,{},{:.17e},{:.17e},{:.17e},{}", row.iteration, row.dilated_volume,
            row.envelope, row.maximum_ratio, row.feasible);
    }
    if let Some(best) = &refitted.best { export(best, &caps, op.leaves(), op.nodes(), 1, "optimized-incumbent"); }
    if !refitted.completed_and_rechecked() {
        return Err("candidate refit did not complete and recheck every cap; source and solved candidate prefixes remain exported".into());
    }
    let point = refitted.recheck.as_ref().and_then(|r| r.as_ref().ok())
        .ok_or("completed refit lacks a checked point")?.clone();
    // Always assess the NEW design; the preceding probe described old material.
    let check = assess(&mut refitted.study, &target_tree, &point, &caps, options,
        &mut control, &mut geometry, 1)?;
    if !check.stable_and_feasible(options.design.optimizer.tolerance, CHANGE_TOLERANCE) {
        return Err("one-refit limit reached without a current two-grid cap pass; no continuum or accuracy guarantee is claimed".into());
    }
    Ok(())
}

fn export(point: &EnvelopeDesignPoint3, caps: &ComplianceCaps3, keys: &[Octant3],
    nodes: &[[f64; 3]], grid: usize, label: &str) {
    let p = &point.evaluation.physical;
    eprintln!("exported_grid={grid}; exported_kind={label}; physical_caps_satisfied={}; dilated_volume={:.17e}; continuum_load_guarantee=false; global_optimum=false",
        caps.satisfied_by(&p.eroded.objective.case_compliances), p.dilated_volume_fraction);
    println!("grid,kind,case,eroded_compliance,cap,cap_satisfied");
    for (i, (&c, &cap)) in p.eroded.objective.case_compliances.iter().zip(caps.caps()).enumerate() {
        println!("{grid},{label},{i},{c:.17e},{cap:.17e},{}", c <= cap);
    }
    println!("grid,kind,level,x,y,z,raw,eroded,nominal,dilated");
    for (i, key) in keys.iter().enumerate() {
        let [x, y, z] = key.index();
        println!("{grid},{label},{},{x},{y},{z},{:.17e},{:.17e},{:.17e},{:.17e}", key.level(),
            point.rho[i], p.eroded.projected_rho[i], p.nominal_rho[i], p.dilated_rho[i]);
    }
    println!("grid,kind,case,node_x,node_y,node_z,ux,uy,uz");
    for (case, displacement) in p.eroded.objective.displacements.iter().enumerate() {
        for (i, xyz) in nodes.iter().enumerate() {
            println!("{grid},{label},{case},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                xyz[0], xyz[1], xyz[2], displacement[3*i], displacement[3*i+1], displacement[3*i+2]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn args(xs: &[&str]) -> Vec<String> { xs.iter().map(|x| (*x).to_owned()).collect() }
    #[test]
    fn authored_caps_and_zero_updates_are_preserved() {
        let (caps, steps, budget) = parse(&args(&["1000", "1500", "0", "12345"])).unwrap();
        assert_eq!(caps.caps(), &[1000.0, 1500.0]);
        assert_eq!((steps, budget), (0, 12345));
    }
    #[test]
    fn malformed_policies_refuse_before_geometry() {
        for xs in [vec![], vec!["1000"], vec!["NaN", "1"], vec!["1", "inf"],
            vec!["0", "1"], vec!["-1", "1"], vec!["1", "1", "bad"],
            vec!["1", "1", "1", "bad"], vec!["1", "1", "1", "1", "extra"]] {
            assert!(parse(&args(&xs)).is_err());
        }
    }
    #[test]
    fn broad_caps_execute_the_actual_coarse_and_fine_path_without_updates() {
        run(&args(&["1e6", "2e6", "0", "250000"])).unwrap();
    }
}
