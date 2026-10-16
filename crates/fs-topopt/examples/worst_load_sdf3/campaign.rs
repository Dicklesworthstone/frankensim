//! Repeated minimum-material refits with unchanged authored load requirements.
use super::*;
use fs_cutfem::elastic3::{ElasticityError3, adaptive::AdaptiveElasticity3, surface::ReferenceLoad3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::QuadratureWork3;
use fs_topopt::sdf3::AdaptiveSdf3Elasticity;
use fs_topopt::sdf3::robust::envelope::limits::{ComplianceCaps3,
    refit::{ComplianceCapRefit3, ComplianceCapRefitOptions3},
    campaign::{ComplianceCapCampaign3, ComplianceCapCampaignOptions3, ComplianceCapCampaignStop3 as Stop}};
use fs_topopt::sdf3_goal::{GoalPreconditioner3, GoalRefinementError3};

type Error = Box<dyn std::error::Error>;
fn parse(args: &[String]) -> Result<(ComplianceCaps3, ComplianceCapCampaignOptions3, usize, usize), Error> {
    if !(3..=6).contains(&args.len()) {
        return Err("usage: worst_load_sdf3 --adapt-material CAP_Y CAP_Z MAX_REFITS [UPDATES_PER_GRID [TOTAL_KRYLOV_ITERATIONS [CAP_CHANGE_TOLERANCE]]]".into());
    }
    let caps = ComplianceCaps3::new(vec![args[0].parse()?, args[1].parse()?], 0.03, 2)?;
    let refits = args[2].parse::<usize>()?;
    let tolerance = args.get(5).map_or(Ok(0.02), |s| s.parse::<f64>())?;
    if refits > 4 || !tolerance.is_finite() || tolerance < 0.0 {
        return Err("example requires 0..=4 refits and a finite nonnegative cap-change tolerance".into());
    }
    let mut options = ComplianceCapCampaignOptions3 { max_refits: refits,
        max_candidate_attempts: refits, change_tolerance: tolerance, ..Default::default() };
    options.enrichment.preconditioner = GoalPreconditioner3::Jacobi { max_contributions: 100_000_000 };
    Ok((caps, options, args.get(3).map_or(Ok(200), |s| s.parse::<usize>())?,
        args.get(4).map_or(Ok(250_000), |s| s.parse::<usize>())?))
}

pub fn run(args: &[String]) -> Result<(), Error> {
    let (caps, options, steps, budget) = parse(args)?;
    let refit = ComplianceCapRefitOptions3 { accepted_steps: steps, ..Default::default() };
    let params = SimpParams { penal: 2.0, beta: 2.0, e_min: 1e-3, ..Default::default() };
    // Reserve one level for the probe of the last allowed design grid.
    let tree = Octree3::uniform(1, options.max_refits as u8 + 2, 4096)?;
    let limits = QuadratureOptions3::default();
    let mut spent = QuadratureWork3::default();
    let mut build = |tree: &Octree3, cp: &mut dyn FnMut() -> ControlFlow<()>|
        -> Result<AdaptiveSolveSpace3, GoalRefinementError3> {
        let remaining = QuadratureOptions3 { max_boxes: limits.max_boxes.saturating_sub(spent.boxes),
            max_points: limits.max_points.saturating_sub(spent.points), ..limits };
        let mut poll = |_| cp();
        let mut q = QuadratureControl3::new(remaining, &mut poll).map_err(ElasticityError3::from)?;
        let result = AdaptiveElasticity3::build(HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(), tree, &Slab,
            &IsotropicElastic::new(1.0, 0.3, 1.0).unwrap(), &|p| p[0] == 0.0,
            ElasticityOptions3::default(), &mut q);
        let work = q.work();
        spent.boxes += work.boxes; spent.points += work.points; spent.field_evaluations += work.field_evaluations;
        Ok(AdaptiveSolveSpace3::jacobi(result?, 100_000_000))
    };
    let op = build(&tree, &mut || ControlFlow::Continue(()))?;
    let study = CutDensityStudy3::new(op, 0.15, params);
    let raw = vec![0.3; study.cells()];
    let y = |p: [f64; 3]| [0.0, -1.0-p[1], 0.0];
    let z = |p: [f64; 3]| [0.0, 0.0, -1.0-p[0]];
    let laws = [ReferenceLoad3::body(&y), ReferenceLoad3::body(&z)];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: budget, ..Default::default() }, &mut poll);
    eprintln!("units=dimensionless; seeds=none; params={params:?}; caps={caps:?}; refit_policy={refit:?}; campaign_policy={options:?}; geometry_limits={limits:?}; total_krylov_budget={budget}; filter_radius=0.15");
    let initial = study.solve_reference_compliance_caps(&raw, &laws, &caps, refit, &mut control)?;
    export(&initial, 0);
    if !initial.completed_and_rechecked() {
        return Err("initial solve did not finish and recheck its caps; any solved prefix was exported".into());
    }
    let mut campaign = ComplianceCapCampaign3::new(initial, tree, &laws, options, &mut control, &mut build)?;
    let (mut printed_grids, mut printed_probes) = (0, 0);
    // Exercise the segmented API: each call promotes at most one design grid,
    // retaining its completed probe for the next segment rather than replaying it.
    let outcome = loop {
        let result = campaign.run(1);
        for (i, grid) in campaign.grids().iter().enumerate().skip(printed_grids) {
            eprintln!("grid={i}; active_cells={}; background_cells={}; dilated_volume={:.17e}; rechecked_compliances={:?}",
                grid.active_cells, grid.background_cells, grid.dilated_volume, grid.compliances);
        }
        for probe in campaign.probes().iter().skip(printed_probes) {
            println!("grid,case,cap,coarse_compliance,fine_compliance,cap_normalized_change,both_caps_satisfied");
            for (i, c) in probe.checks.iter().enumerate() {
                println!("{},{i},{:.17e},{:.17e},{:.17e},{:.17e},{}", probe.source_grid,
                    c.cap, c.coarse, c.fine, c.normalized_change, c.satisfied());
            }
            eprintln!("assessed_grid={}; coarse_envelope={:.17e}; fine_envelope={:.17e}; current_two_grid_pass={}",
                probe.source_grid, probe.coarse_envelope, probe.fine_envelope, probe.passed);
        }
        printed_grids = campaign.grids().len(); printed_probes = campaign.probes().len();
        if matches!(result, Ok(Stop::Paused)) { continue; }
        break result;
    };
    // Preserve installed fields before propagating any build, solve or final
    // probe error. Never export pending/rejected fields under the installed grid.
    export(campaign.current(), campaign.grids().len()-1);
    eprintln!("campaign_outcome={outcome:?}; candidate_attempts={}; pending_complete_candidate={}; rejected_candidate={}; cumulative_solve_work={:?}; continuum_load_guarantee=false; volume_error_estimated=false; cross_grid_descent_claimed=false; durable_resume=false",
        campaign.candidate_attempts(), campaign.pending().is_some(), campaign.rejected().is_some(), campaign.work());
    drop(campaign); drop(build);
    eprintln!("cumulative_geometry_work={spent:?}");
    if outcome? != Stop::StableAndFeasible {
        return Err("campaign stopped without a current two-grid cap pass; the installed feasible design remains exported".into());
    }
    Ok(())
}

fn export<O: AdaptiveSdf3Elasticity>(solved: &ComplianceCapRefit3<O>, grid: usize) {
    let point = solved.recheck.as_ref().and_then(|r| r.as_ref().ok())
        .or(solved.best.as_ref()).or(solved.restoration.point.as_ref());
    let Some(point) = point else { return; };
    let label = if solved.best.is_some() { "optimized-incumbent" } else { "restoration-seed" };
    let physical = &point.evaluation.physical; let op = solved.study.operator().adaptive();
    eprintln!("exported_grid={grid}; kind={label}; completed_and_rechecked={}; restoration={:?}; optimizer_outcome={:?}; physical_caps_satisfied={}",
        solved.completed_and_rechecked(), solved.restoration.stop, solved.outcome,
        solved.caps.satisfied_by(&physical.eroded.objective.case_compliances));
    println!("grid,kind,level,x,y,z,raw,eroded,nominal,dilated");
    for (i, key) in op.leaves().iter().enumerate() {
        let [x,y,z] = key.index();
        println!("{grid},{label},{},{x},{y},{z},{:.17e},{:.17e},{:.17e},{:.17e}", key.level(), point.rho[i],
            physical.eroded.projected_rho[i], physical.nominal_rho[i], physical.dilated_rho[i]);
    }
    println!("grid,kind,case,node_x,node_y,node_z,ux,uy,uz");
    for (case, u) in physical.eroded.objective.displacements.iter().enumerate() {
        for (i, p) in op.nodes().iter().enumerate() {
            println!("{grid},{label},{case},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                p[0], p[1], p[2], u[3*i], u[3*i+1], u[3*i+2]);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn args(xs: &[&str]) -> Vec<String> { xs.iter().map(|s| (*s).to_owned()).collect() }
    #[test]
    fn authored_caps_round_count_and_tolerance_are_explicit() {
        let (caps, options, steps, budget) = parse(&args(&["1000", "1500", "3", "0", "1234", "0.01"])).unwrap();
        assert_eq!(caps.caps(), &[1000.0, 1500.0]);
        assert_eq!((options.max_refits, options.max_candidate_attempts), (3, 3));
        assert_eq!(options.change_tolerance, 0.01); assert_eq!((steps, budget), (0, 1234));
    }
    #[test]
    fn malformed_or_unbounded_campaigns_refuse_before_geometry() {
        for xs in [vec![], vec!["1", "2"], vec!["NaN", "2", "1"], vec!["1", "2", "5"],
            vec!["1", "2", "1", "0", "10", "NaN"], vec!["1", "2", "1", "0", "10", "-1"],
            vec!["1", "2", "1", "0", "10", "0.1", "extra"]] { assert!(parse(&args(&xs)).is_err()); }
    }
    #[test]
    fn broad_caps_run_real_initial_and_current_grid_checks() {
        run(&args(&["1e6", "2e6", "0", "0", "250000"])).unwrap();
    }
}
