//! Repeat the existing full-spectrum assessment/refit, keeping one work budget.
use super::*;
use fs_cutfem::elastic3::ElasticityError3;
use fs_cutfem::quad3::QuadratureWork3;
use fs_topopt::sdf3::robust::envelope::load_space::design::limits::campaign::{
    LoadSpaceCapCampaign3, LoadSpaceCapCampaignOptions3, LoadSpaceCapCampaignStop3,
};

pub(super) fn run(args: &[String]) -> CommandResult {
    if !(1..=5).contains(&args.len()) {
        return Err("usage: worst_load_sdf3 --load-space --adapt-material CAP [UPDATES [TOTAL_KRYLOV_ITERATIONS [GRAM_CHANGE_TOLERANCE [MAX_REFITS]]]]".into());
    }
    let cap = LoadSpaceCap3::new(2,args[0].parse()?,0.03)?;
    let updates = args.get(1).map_or(Ok(200),|s|s.parse::<usize>())?;
    let iterations = args.get(2).map_or(Ok(250_000),|s|s.parse::<usize>())?;
    let change = args.get(3).map_or(Ok(0.02),|s|s.parse::<f64>())?;
    let max_refits = args.get(4).map_or(Ok(2),|s|s.parse::<usize>())?;
    if !change.is_finite() || change <= 0.0 {
        return Err("Gram-change tolerance must be finite and positive".into());
    }
    let tree = Octree3::uniform(1,4,4096)?;
    let params = SimpParams { penal:2.0,beta:2.0,e_min:1e-3,..Default::default() };
    let refit = LoadSpaceCapRefitOptions3 {accepted_steps:updates,..Default::default()};
    let policy = LoadSpaceCapCampaignOptions3 {
        max_refits, max_candidate_attempts:max_refits, change_tolerance:change,
        enrichment:GoalRefinementOptions3 { preconditioner:GoalPreconditioner3::Jacobi {
            max_contributions:DIAGONAL_LIMIT }, ..Default::default() },
    };
    let y = |p: [f64;3]| [0.0,-1.0-p[1],0.0];
    let yz = |p: [f64;3]| [0.0,-0.6*(1.0+p[1]),-1.0-p[0]];
    let laws = [ReferenceLoad3::body(&y),ReferenceLoad3::body(&yz)];
    let domain = HexCell::try_new([0.0;3],[1.0;3])?;
    let material = IsotropicElastic::new(1.0,0.3,1.0)?;
    let limits = QuadratureOptions3::default();
    let mut geometry = QuadratureWork3::default();
    let mut build = |tree: &Octree3, cp: &mut dyn FnMut() -> ControlFlow<()>|
        -> Result<AdaptiveSolveSpace3,GoalRefinementError3> {
        let remaining = QuadratureOptions3 {
            max_boxes:limits.max_boxes.saturating_sub(geometry.boxes),
            max_points:limits.max_points.saturating_sub(geometry.points),..limits
        };
        let mut poll = |_| cp();
        let mut q = QuadratureControl3::new(remaining,&mut poll).map_err(ElasticityError3::from)?;
        let op = AdaptiveElasticity3::build(domain,tree,&Slab,&material,&|p|p[0]==0.0,
            ElasticityOptions3::default(),&mut q);
        // Failed builds spend geometry work too. No fresh allowance per grid.
        let spent = q.work();
        geometry.boxes += spent.boxes; geometry.points += spent.points;
        geometry.field_evaluations += spent.field_evaluations;
        Ok(AdaptiveSolveSpace3::jacobi(op?,DIAGONAL_LIMIT))
    };
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget {total_iterations:iterations,..Default::default()},&mut poll);
    eprintln!("units=dimensionless; load_set=Fq; coefficient_norm_limit=1; authored_cap={}; params={params:?}; refit_policy={refit:?}; campaign_policy={policy:?}; refine_all_active_cells=true; filter_radius=0.15; total_krylov_budget={iterations}; geometry_limits={limits:?}; max_tree_level=4; max_background_leaves=4096; continuum_certificate=false",
        cap.cap());
    let outcome = (|| -> CommandResult {
        control.checkpoint("sdf3-load-space-adapt-command-start")?;
        let op = build(&tree,&mut || {
            if control.checkpoint("sdf3-load-space-adapt-command-build").is_ok() {ControlFlow::Continue(())}
            else {ControlFlow::Break(())}
        })?;
        let study = CutDensityStudy3::new(op,0.15,params);
        let raw = vec![0.3;study.cells()];
        let initial = study.solve_reference_load_space_cap(&raw,&laws,&cap,refit,&mut control)?;
        export(&initial,0);
        if !initial.completed_and_rechecked() {
            return Err("initial combined-load design stopped; exported fields are only its solved prefix".into());
        }
        let mut campaign = LoadSpaceCapCampaign3::new(initial,tree,&laws,policy,&mut control,&mut build)?;
        let (mut printed_grid,mut printed_probes) = (0,0);
        loop {
            // One segment per grid makes complete installed fields available
            // for export before the next adaptation can fail.
            let result = campaign.run(1);
            for p in &campaign.probes()[printed_probes..] {
                println!("assessed_grid,coarse_worst_compliance,fine_worst_compliance,cap,normalized_full_gram_change,passed");
                println!("{},{:.17e},{:.17e},{:.17e},{:.17e},{}",p.source_grid,
                    p.coarse_worst_compliance,p.fine_worst_compliance,cap.cap(),p.normalized_gram_change,p.passed);
            }
            printed_probes = campaign.probes().len();
            let grid = campaign.grids().len()-1;
            if grid != printed_grid { export(campaign.current(),grid); printed_grid = grid; }
            eprintln!("adaptation_outcome={result:?}; current_grid={grid}; candidate_attempts={}; pending_candidate={}; rejected_candidate={}; physical_work={:?}; refit_limit_is_accuracy=false; optimizer_stop_is_accuracy=false",
                campaign.candidate_attempts(),campaign.pending().is_some(),campaign.rejected().is_some(),campaign.work());
            match result? {
                LoadSpaceCapCampaignStop3::Paused => continue,
                LoadSpaceCapCampaignStop3::StableAndFeasible => return Ok(()),
                stop => return Err(format!("combined-load adaptation stopped at {stop:?}; installed fields were exported but the two-grid requirement did not pass").into()),
            }
        }
    })();
    drop(build);
    eprintln!("final_physical_work={:?}; total_geometry_work={geometry:?}; binary_material_checked=false; cross_grid_descent_claimed=false",
        control.work());
    outcome
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn input_admission_precedes_geometry_and_zero_refits_still_assess_the_current_grid() {
        for args in [vec![],vec!["0"],vec!["NaN"],vec!["1000","bad"],
            vec!["1000","0","250000","NaN"],vec!["1000","0","250000","0"],
            vec!["1000","0","250000","0.1","bad"],vec!["1000","0","250000","0.1","1","extra"]] {
            assert!(run(&args.into_iter().map(String::from).collect::<Vec<_>>()).is_err());
        }
        // The full outer route reaches actual physics, even with no refit budget.
        super::super::super::run(&["--adapt-material","1000000000000","0","250000","1","0"]
            .map(String::from)).unwrap();
    }
    #[test]
    fn a_failed_check_at_the_refit_limit_is_not_command_success() {
        let error = run(&["1000000000000","0","250000","1e-30","0"].map(String::from)).unwrap_err();
        assert!(error.to_string().contains("RefinementLimit"),"{error}");
        assert!(run(&["1000000000000","0","1","1","1"].map(String::from)).is_err());
    }
}
