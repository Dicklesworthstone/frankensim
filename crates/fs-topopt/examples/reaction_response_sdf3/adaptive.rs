//! Bounded mixed-response campaigns using the library's actual adaptive owner.
use super::*;
use fs_topopt::sdf3::response::adaptive::{AdaptiveResponseOptions3, AdaptiveResponseStudy3};

type Error = Box<dyn std::error::Error>;

#[derive(Debug, Clone, Copy)]
pub(super) struct Config {
    rounds: usize,
    objective_tolerance: f64,
    response_tolerance: f64,
    volume_tolerance: f64,
}
impl Default for Config {
    fn default() -> Self {
        let defaults = AdaptiveResponseOptions3::default();
        Self { rounds: 3, objective_tolerance: defaults.objective_tolerance,
            response_tolerance: defaults.response_tolerance, volume_tolerance: defaults.volume_tolerance }
    }
}

// Keep the legacy positional/one-shot parser unchanged. Consume new flags only
// for an explicitly requested campaign, and never silently ignore a tolerance.
pub(super) fn parse(args: &mut Vec<String>) -> Result<Option<Config>, Error> {
    let flags = ["--adapt-rounds", "--objective-tol", "--response-tol", "--volume-tol"];
    if !args.iter().any(|s| s == flags[0]) {
        if args.iter().any(|s| flags[1..].contains(&s.as_str())) {
            return Err("adaptive tolerances require --adapt-rounds".into());
        }
        return Ok(None);
    }
    let mut config = Config::default();
    let mut seen = [false; 4];
    let mut positional = Vec::new();
    let mut i = 0;
    while i < args.len() {
        let arg = &args[i];
        if let Some(flag) = flags.iter().position(|f| *f == arg.as_str()) {
            if seen[flag] { return Err(format!("duplicate adaptive flag: {arg}").into()); }
            seen[flag] = true;
            let value = args.get(i+1).ok_or_else(|| format!("missing value after {arg}"))?;
            if flag == 0 {
                config.rounds = value.parse()?;
            } else {
                let tolerance: f64 = value.parse()?;
                if !tolerance.is_finite() || tolerance < 0.0 {
                    return Err(format!("{arg} requires a finite nonnegative tolerance").into());
                }
                match flag { 1 => config.objective_tolerance = tolerance,
                    2 => config.response_tolerance = tolerance, _ => config.volume_tolerance = tolerance }
            }
            i += 2;
        } else if arg.starts_with("--") {
            return Err(format!("unsupported or conflicting campaign flag: {arg}; do not combine --adapt-rounds with --refit or --estimate").into());
        } else {
            positional.push(arg.clone());
            i += 1;
        }
    }
    if !(2..=4).contains(&config.rounds) {
        return Err("--adapt-rounds requires 2..=4 total grids including the initial grid".into());
    }
    *args = positional;
    Ok(Some(config))
}

fn build(tree: &Octree3, geometry: &mut QuadratureControl3<'_>) -> Result<AdaptiveSolveSpace3, Error> {
    let op = AdaptiveElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0; 3], [1.0; 3])?,
        tree, &Slab, &IsotropicElastic::new(1.0, 0.0, 1.0)?, &|_| false, &|_, _| true,
        ElasticityOptions3::default(), Default::default(), Default::default(), geometry)?;
    Ok(AdaptiveSolveSpace3::jacobi(op, 100_000_000))
}

pub(super) fn run(steps: usize, iterations: usize, reaction: f64, config: Config) -> Result<(), Error> {
    let options = AdaptiveResponseOptions3 {
        steps_per_grid: steps, max_refinements: config.rounds-1, max_marks: 2,
        objective_tolerance: config.objective_tolerance, response_tolerance: config.response_tolerance,
        volume_tolerance: config.volume_tolerance, ..Default::default()
    };
    eprintln!("mode=adaptive-reaction; max_grids={}; updates_per_grid={steps}; krylov_budget={iterations}; objective_scale={}; objective_tolerance={}; response_tolerance={}; volume_tolerance={}; max_marks={}; theta={}; max_level=5; max_leaves=4096",
        config.rounds, options.fit.objective_scale, options.objective_tolerance,
        options.response_tolerance, options.volume_tolerance, options.max_marks, options.theta);
    let mut gp = |_| ControlFlow::Continue(());
    let mut geometry = QuadratureControl3::new(QuadratureOptions3::default(), &mut gp)?;
    // Reserve an extra level for a probe after the final allowed promotion.
    // The existing leaf, quadrature and Krylov caps still refuse excess work.
    let tree = Octree3::uniform(1, 5, 4096)?;
    let mut study = CutDensityStudy3::new(build(&tree, &mut geometry)?, 0.15, SimpParams::default());
    let observe = |_:[f64; 3]| [1.0, 0.0, 0.0];
    let targets = [ReferenceResponseTarget3 { observation: ReferenceLoad3::body(&observe),
        target: 0.02*0.66*0.66/2.0, scale: 0.005, weight: 0.1 }];
    let reactions = [ReactionTarget3 { mode: &right, target: reaction, scale: 0.002, weight: 1.0 }];
    let rows: [&[ReactionTarget3<'_>]; 1] = [&reactions];
    let cases = [ReferenceResponseCase3 { load: ReferenceLoad3::default(), prescribed: Some(&motion), targets: &targets }];
    let initial = vec![0.5; study.cells()];
    let mut sp = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: iterations, ..Default::default() }, &mut sp);
    let mut session = AdaptiveResponseStudy3::new(&mut study, tree, &cases, Some(&rows), &initial, options, &mut control)?;
    let outcome = session.run(options.max_refinements, |tree, _| build(tree, &mut geometry), &mut control);
    // Emit the LAST promoted feasible design before returning any late error.
    // Candidate failures never relabel that source's fields or baseline.
    println!("adaptive_round,source_cells,retained_cells,source_objective,new_baseline_objective,retained_objective,marking_fraction,marking_target_met,selected_latest");
    for (i, row) in session.rounds().iter().enumerate() {
        println!("{},{},{},{:.17e},{:.17e},{:.17e},{:.17e},{},{}", i+1, row.cells[0], row.cells[1],
            row.objectives[0], row.objectives[1], row.objectives[2], row.marking.achieved_fraction,
            row.marking.target_met, row.selected_latest);
        eprintln!("source_round={i}; two_grid_changes={:?}; marked={:?}; fine_al_stop={:?}; fine_al_kkt={:?}; selected_latest={}; cumulative_work={:?}; cross_grid_descent_claimed=false",
            row.source_change, row.marking.marked, row.optimizer.stop, row.optimizer.kkt, row.selected_latest, row.work);
    }
    let round = session.rounds().len();
    let accepted = session.accepted();
    println!("retained_round,cell_level,cell_i,cell_j,cell_k,raw_density,projected_density");
    for ((leaf, raw), physical) in session.study().operator().elasticity().leaves().iter()
        .zip(&accepted.rho).zip(&accepted.projected_rho) {
        let [i, j, k] = leaf.index();
        println!("{round},{},{i},{j},{k},{raw:.17e},{physical:.17e}", leaf.level());
    }
    println!("retained_round,reaction,target_reaction,displacement_integral,target_displacement_integral,objective,volume_fraction");
    println!("{round},{:.17e},{reaction:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
        accepted.reaction_responses[0][0], accepted.responses[0][0], targets[0].target,
        accepted.objective, accepted.volume_fraction);
    let violation = (accepted.volume_fraction-options.fit.volume_cap).max(0.0);
    eprintln!("retained_round={round}; volume_violation={violation:.9e}; selected_latest={}; latest_al_outcome={:?}; current_probe_changes={:?}; unpromoted_candidate_retained={}; work={:?}; geometry_points={}; targets_source=declared_design_objective; continuum_certified=false; unique_recovery_claimed=false; actuator_energy_claimed=false",
        accepted.rho == session.fit().accepted.rho, session.fit().outcome, session.measured_change(),
        session.unpromoted().is_some(), control.work(), geometry.work().points);
    let stop = outcome?;
    eprintln!("adaptive_stop={stop:?}; numerical_tolerance_pass={}; optimizer_convergence_is_separate=true",
        stop == fs_topopt::sdf3::response::adaptive::AdaptiveResponseStop3::TwoGridToleranceMet);
    if violation > options.fit.optimizer.tolerance { return Err("retained adaptive design violates volume cap".into()); }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    fn args(s: &str) -> Vec<String> { s.split_whitespace().map(str::to_owned).collect() }
    #[test]
    fn g0_adaptive_flags_keep_positional_arguments_and_independent_tolerances() {
        let mut a = args("40 --adapt-rounds 3 250000 0.0012 --response-tol 0.002 --objective-tol 0 --volume-tol 0.0001");
        let config = parse(&mut a).unwrap().unwrap();
        assert_eq!(a, args("40 250000 0.0012")); assert_eq!(config.rounds, 3);
        assert_eq!(config.response_tolerance, 0.002); assert_eq!(config.objective_tolerance, 0.0);
        assert_eq!(config.volume_tolerance, 0.0001);
        let defaults = parse(&mut args("--adapt-rounds 2")).unwrap().unwrap();
        assert_eq!(defaults.response_tolerance, AdaptiveResponseOptions3::default().response_tolerance);
    }
    #[test]
    fn g0_invalid_campaign_controls_are_refused_without_modifying_the_arguments() {
        for s in ["--response-tol 0.1", "--adapt-rounds", "--adapt-rounds 1", "--adapt-rounds 5",
            "--adapt-rounds 3 --adapt-rounds 2", "--adapt-rounds 3 --response-tol NaN",
            "--adapt-rounds 3 --volume-tol -1", "--adapt-rounds 3 --objective-tol inf",
            "--adapt-rounds 3 --refit", "--estimate --adapt-rounds 3", "--adapt-rounds 2 --unknown",
            "--adapt-rounds 3 --response-tol 1 --response-tol 2"] {
            let mut a = args(s); let before = a.clone();
            assert!(parse(&mut a).is_err(), "{s}"); assert_eq!(a, before);
        }
    }
    #[test]
    fn g5_legacy_fit_estimate_and_refit_arguments_are_untouched() {
        for s in ["", "40 250000 0.0012", "40 --estimate", "40 250000 -0.001 --refit"] {
            let mut a = args(s); let before = a.clone();
            assert!(parse(&mut a).unwrap().is_none()); assert_eq!(a, before);
        }
    }
}
