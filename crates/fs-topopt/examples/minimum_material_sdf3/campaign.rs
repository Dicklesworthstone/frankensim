//! Repeated selective design, with a new per-limit probe of the final endpoint.
use super::*;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_topopt::sdf3::response::limits::refinement::*;
use fs_topopt::sdf3::response::limits::refinement::adaptive::*;

type Error = Box<dyn std::error::Error>;
#[derive(Clone, Copy)]
pub(super) struct Config { grids: usize, response_tolerance: f64, volume_tolerance: f64 }

/// Leave legacy options alone unless a campaign is explicitly selected.
pub(super) fn parse(args: &mut Vec<String>) -> Result<Option<Config>, Error> {
    let flags = ["--adapt-rounds", "--response-tol", "--volume-tol"];
    if !args.iter().any(|a| a == flags[0]) {
        if args.iter().any(|a| flags[1..].contains(&a.as_str())) {
            return Err("adaptive tolerances require --adapt-rounds".into());
        }
        return Ok(None);
    }
    let defaults = AdaptiveMinimumMaterialOptions3::default();
    let mut config = Config { grids: 2, response_tolerance: defaults.response_tolerance,
        volume_tolerance: defaults.volume_tolerance };
    let mut seen = [false; 3]; let mut positional = Vec::new(); let mut i = 0;
    while i < args.len() {
        let arg = &args[i];
        if let Some(flag) = flags.iter().position(|f| *f == arg.as_str()) {
            if seen[flag] { return Err(format!("duplicate option: {arg}").into()); }
            seen[flag] = true;
            let value = args.get(i+1).ok_or_else(|| format!("missing value for {arg}"))?;
            if flag == 0 { config.grids = value.parse()?; }
            else {
                let tolerance: f64 = value.parse()?;
                if !tolerance.is_finite() || tolerance < 0.0 { return Err(format!("invalid {arg}").into()); }
                if flag == 1 { config.response_tolerance = tolerance; } else { config.volume_tolerance = tolerance; }
            }
            i += 2;
        } else if arg.starts_with("--") {
            return Err(format!("unsupported or conflicting campaign option: {arg}").into());
        } else { positional.push(arg.clone()); i += 1; }
    }
    if !(2..=4).contains(&config.grids) { return Err("--adapt-rounds must be 2..=4 total retained grids".into()); }
    *args = positional; Ok(Some(config))
}
fn build(tree: &Octree3, geometry: &mut QuadratureControl3<'_>) -> Result<AdaptiveSolveSpace3, Error> {
    let op = AdaptiveElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0; 3], [1.0; 3])?,
        tree, &Slab, &IsotropicElastic::new(1.0, 0.0, 1.0)?, &|_| false, &|_, _| true,
        ElasticityOptions3::default(), Default::default(), Default::default(), geometry)?;
    Ok(AdaptiveSolveSpace3::jacobi(op, 100_000_000))
}
pub(super) fn run(steps: usize, budget: usize, lo: f64, hi: f64, config: Config) -> Result<(), Error> {
    let options = AdaptiveMinimumMaterialOptions3 { steps_per_grid: steps, max_refinements: config.grids-1,
        max_marks: 2, response_tolerance: config.response_tolerance,
        volume_tolerance: config.volume_tolerance, ..Default::default() };
    let mut gp = |_| ControlFlow::Continue(());
    let mut geometry = QuadratureControl3::new(QuadratureOptions3::default(), &mut gp)?;
    // Reserve a level for assessment AFTER the final permitted selective refit.
    let tree = Octree3::uniform(1, 5, 4096)?;
    let mut study = CutDensityStudy3::new(build(&tree, &mut geometry)?, 0.15, SimpParams::default());
    let observe = |_:[f64; 3]| [1.0, 0.0, 0.0];
    let tension = [
        ReferenceResponseLimit3 { response: ReferenceBoundedResponse3::Reaction(&right),
            sense: ResponseLimitSense3::AtLeast, bound: lo, scale: 0.002 },
        ReferenceResponseLimit3 { response: ReferenceBoundedResponse3::Reaction(&right),
            sense: ResponseLimitSense3::AtMost, bound: hi, scale: 0.002 },
        ReferenceResponseLimit3 { response: ReferenceBoundedResponse3::Displacement(ReferenceLoad3::body(&observe)),
            sense: ResponseLimitSense3::AtMost, bound: 0.006, scale: 0.005 },
    ];
    let compressed = [ReferenceResponseLimit3 { sense: ResponseLimitSense3::AtMost, bound: -lo, ..tension[0] },
        ReferenceResponseLimit3 { sense: ResponseLimitSense3::AtLeast, bound: -hi, ..tension[1] },
        ReferenceResponseLimit3 { sense: ResponseLimitSense3::AtLeast, bound: -0.006, ..tension[2] }];
    let cases = [ReferenceResponseLimitCase3 { load: ReferenceLoad3::default(), prescribed: Some(&motion), limits: &tension },
        ReferenceResponseLimitCase3 { load: ReferenceLoad3::default(), prescribed: Some(&compression), limits: &compressed }];
    let rho = vec![0.5; study.cells()];
    let mut sp = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: budget, ..Default::default() }, &mut sp);
    let mut session = AdaptiveMinimumMaterialStudy3::new(&mut study, tree, &cases, &rho, options, &mut control)?;
    let outcome = session.run(options.max_refinements, |tree, _| build(tree, &mut geometry), &mut control);
    // The retained source stays useful even after a late candidate/probe failure.
    println!("round,source_cells,retained_cells,source_volume,new_baseline_volume,retained_volume,marking_fraction,mark_target_met");
    for (i, row) in session.rounds().iter().enumerate() {
        println!("{},{},{},{:.17e},{:.17e},{:.17e},{:.17e},{}", i+1, row.cells[0], row.cells[1],
            row.volumes[0], row.volumes[1], row.volumes[2], row.marking.achieved_fraction, row.marking.target_met);
        eprintln!("round={}; source_assessment={:?}; latest_al_stop={:?}; latest_al_kkt={:?}; selected_latest={}; work={:?}",
            i+1, row.source_assessment, row.optimizer.stop, row.optimizer.kkt, row.selected_latest, row.work);
    }
    let accepted = session.accepted();
    println!("retained_case,limit,sense,value,bound,scale,residual");
    for (i, case) in cases.iter().enumerate() { for (j, limit) in case.limits.iter().enumerate() {
        println!("{i},{j},{:?},{:.17e},{:.17e},{:.17e},{:.17e}", limit.sense,
            accepted.measurements[i][j], limit.bound, limit.scale, accepted.residuals[i][j]);
    } }
    if let Some(probe) = session.probe() {
        println!("probe_case,limit,source_value,probe_value,source_residual,probe_residual");
        for (i, case) in probe.cases.iter().enumerate() { for (j, l) in case.limits.iter().enumerate() {
            println!("{i},{j},{:.17e},{:.17e},{:.17e},{:.17e}", l.goal.coarse_value, l.goal.fine_value,
                l.coarse_residual, l.fine_residual);
        } }
    }
    println!("cell_level,cell_i,cell_j,cell_k,raw_density,projected_density");
    for ((leaf, raw), projected) in session.study().operator().elasticity().leaves().iter()
        .zip(&accepted.rho).zip(&accepted.projected_rho) {
        let [i,j,k] = leaf.index(); println!("{},{i},{j},{k},{raw:.17e},{projected:.17e}", leaf.level());
    }
    eprintln!("retained_grids={}; source_feasible={}; final_assessment={:?}; response_tolerance={}; volume_tolerance={}; constraint_tolerance={}; latest_al_outcome={:?}; unpromoted_candidate={}; work={:?}; geometry_points={}; continuum_certified=false",
        session.rounds().len()+1, accepted.feasible(options.fit.optimizer.tolerance), session.assessment(),
        options.response_tolerance, options.volume_tolerance, options.fit.optimizer.tolerance,
        session.fit().outcome, session.unpromoted().is_some(), control.work(), geometry.work().points);
    let stop = outcome?; eprintln!("adaptive_stop={stop:?}");
    if stop != AdaptiveMinimumMaterialStop3::TwoGridLimitsSatisfied {
        return Err(format!("adaptive material design stopped without passing final two-grid limits: {stop:?}").into());
    }
    Ok(())
}
