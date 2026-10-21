//! Minimum-material gray design -> exhaustive binary two-grid qualification.
use super::*;
use fs_cutfem::elastic3::{ElasticityError3,
    adaptive::{AdaptiveElasticity3, enrichment::precondition::AdaptiveSolveSpace3},
    surface::ReferenceLoad3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3, QuadratureWork3};
use fs_topopt::sdf3::robust::envelope::load_space::design::limits::{
    LoadSpaceCap3, refit::LoadSpaceCapRefitOptions3,
    binary::BinaryLoadSpaceScanStop3,
    refinement::binary::scan::BinaryRefinedScanOptions3,
};
use fs_topopt::sdf3_goal::{GoalPreconditioner3, GoalRefinementError3, GoalRefinementOptions3};

pub fn run(args: &[String]) -> CommandResult {
    if !(1..=5).contains(&args.len()) {
        return Err("usage: worst_load_sdf3 --load-space --binary-refined CAP [UPDATES [TOTAL_KRYLOV_ITERATIONS [GRAM_CHANGE_TOLERANCE [BINARY_DILATED_VOLUME_LIMIT]]]]".into());
    }
    let cap = LoadSpaceCap3::new(2,args[0].parse()?,0.03)?;
    let updates = args.get(1).map_or(Ok(200),|s|s.parse::<usize>())?;
    let iterations = args.get(2).map_or(Ok(250_000),|s|s.parse::<usize>())?;
    let change = args.get(3).map_or(Ok(0.02),|s|s.parse::<f64>())?;
    let volume_limit = args.get(4).map_or(Ok(1.0),|s|s.parse::<f64>())?;
    if !change.is_finite() || change < 0.0 {
        return Err("Gram-change tolerance must be finite and nonnegative".into());
    }
    if !volume_limit.is_finite() || !(0.0..=1.0).contains(&volume_limit) {
        return Err("binary dilated-volume limit must be finite and in [0,1]".into());
    }

    let params = SimpParams { penal:2.0,beta:2.0,e_min:1e-3,..Default::default() };
    let gray = LoadSpaceCapRefitOptions3 { accepted_steps:updates,..Default::default() };
    let tree = Octree3::uniform(1,4,4096)?;
    let y = |p:[f64;3]| [0.0,-1.0-p[1],0.0];
    let yz = |p:[f64;3]| [0.0,-0.6*(1.0+p[1]),-1.0-p[0]];
    let laws = [ReferenceLoad3::body(&y),ReferenceLoad3::body(&yz)];
    let domain = HexCell::try_new([0.0;3],[1.0;3])?;
    let material = IsotropicElastic::new(1.0,0.3,1.0)?;
    let limits = QuadratureOptions3::default();
    let mut geometry_work = QuadratureWork3::default();

    let mut build = |tree:&Octree3, checkpoint:&mut dyn FnMut()->ControlFlow<()>|
        -> Result<AdaptiveElasticity3,GoalRefinementError3> {
        let remaining = QuadratureOptions3 {
            max_boxes:limits.max_boxes.saturating_sub(geometry_work.boxes),
            max_points:limits.max_points.saturating_sub(geometry_work.points),
            ..limits
        };
        let mut poll = |_| checkpoint();
        let mut quadrature = QuadratureControl3::new(remaining,&mut poll)
            .map_err(ElasticityError3::from)?;
        let result = AdaptiveElasticity3::build(domain,tree,&Slab,&material,&|p|p[0]==0.0,
            ElasticityOptions3::default(),&mut quadrature);
        let spent = quadrature.work();
        geometry_work.boxes = geometry_work.boxes.saturating_add(spent.boxes);
        geometry_work.points = geometry_work.points.saturating_add(spent.points);
        geometry_work.field_evaluations =
            geometry_work.field_evaluations.saturating_add(spent.field_evaluations);
        Ok(result?)
    };

    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(
        SolveBudget { total_iterations:iterations,..Default::default() },&mut poll);

    let mut stopped = None;
    let initial = build(&tree,&mut || match control.checkpoint("sdf3-binary-refined-command-build") {
        Ok(()) => ControlFlow::Continue(()),
        Err(error) => { stopped = Some(error); ControlFlow::Break(()) },
    });
    if let Some(error) = stopped { return Err(error.into()); }
    let study = CutDensityStudy3::new(
        AdaptiveSolveSpace3::jacobi(initial?,DIAGONAL_LIMIT),0.15,params);
    let raw = vec![0.3;study.cells()];
    let mut gray_design = study.solve_reference_load_space_cap(&raw,&laws,&cap,gray,&mut control)?;
    super::export(&gray_design,0);
    if !gray_design.completed_and_rechecked() {
        return Err("initial combined-load design stopped; exported fields are only its solved prefix".into());
    }
    let rho = gray_design.recheck.as_ref().unwrap().as_ref().unwrap().rho.clone();

    let mut tree_stop = None;
    let refined = tree.refined(gray_design.study.operator().elasticity().leaves(),|| {
        match control.checkpoint("sdf3-binary-refined-command-tree") {
            Ok(()) => ControlFlow::Continue(()),
            Err(error) => { tree_stop = Some(error); ControlFlow::Break(()) },
        }
    });
    if let Some(error) = tree_stop { return Err(error.into()); }
    let fine_tree = refined?;

    let mut scan = BinaryRefinedScanOptions3::default();
    scan.selection.robustness = gray.design.robustness;
    scan.selection.void_stiffness = 1e-3;
    scan.selection.max_load_cases = 2;
    scan.selection.max_candidates = 64;
    scan.selection.dilated_volume_limit = volume_limit;
    scan.selection.smooth_tolerance = gray.design.optimizer.tolerance;
    scan.change_tolerance = change;
    scan.enrichment = GoalRefinementOptions3 {
        preconditioner:GoalPreconditioner3::Jacobi { max_contributions:DIAGONAL_LIMIT },
        ..Default::default()
    };

    eprintln!("mode=exhaustive-binary-two-grid-load-space; authored_cap={}; params={params:?}; gray_policy={gray:?}; binary_policy={scan:?}; total_krylov_budget={iterations}; geometry_limits={limits:?}; continuum_certificate=false; zero_void_stiffness=false",
        cap.cap());

    let report = gray_design.study.scan_all_binary_load_space_cap_enrichment(
        &rho,&laws,params,&cap,scan,&mut control,
        |checkpoint| build(&fine_tree,checkpoint))?;
    drop(build);

    println!("candidate,threshold,dilated_volume,coarse_worst_compliance,fine_worst_compliance,cap,normalized_full_gram_change,two_grid_passed,feasible");
    for row in &report.candidates {
        println!("{},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{},{}",
            row.index,row.threshold,row.dilated_volume,row.coarse_worst_compliance,
            row.fine_worst_compliance,cap.cap(),row.normalized_gram_change,
            row.two_grid_passed,row.feasible);
    }
    eprintln!("binary_refined_stop={:?}; error={:?}; attempted_candidates={}; selected_index={:?}; physical_work={:?}; geometry_work={geometry_work:?}",
        report.stop,report.error,report.attempted_candidates,report.selected_index,report.work);

    let Some(best) = &report.best else {
        return Err("no distinct binary material state passed both physical grids, full-Gram change, and material-volume gates".into());
    };
    let threshold = best.physical.options.threshold;
    eprintln!("exported_source=least-dilated-two-grid-binary-endpoint; threshold={threshold}; coarse_worst={:.17e}; fine_worst={:.17e}; normalized_full_gram_change={:.17e}; cap={:.17e}",
        best.coarse.worst_compliance,best.fine.worst_compliance,best.normalized_gram_change,cap.cap());
    println!("basis_index,coarse_worst_coefficient,fine_worst_coefficient");
    for i in 0..best.coarse.worst_coefficients.len() {
        println!("{i},{:.17e},{:.17e}",
            best.coarse.worst_coefficients[i],best.fine.worst_coefficients[i]);
    }
    let op = gray_design.study.operator().elasticity();
    let scales = best.physical.eroded_stiffness_scales();
    println!("level,cell_x,cell_y,cell_z,eroded_solid,nominal_solid,dilated_solid,binary_eroded_stiffness");
    for (i,key) in op.leaves().iter().enumerate() {
        let [x,y,z] = key.index();
        println!("{},{x},{y},{z},{},{},{},{:.17e}",key.level(),
            u8::from(best.physical.eroded.solid[i]),
            u8::from(best.physical.nominal.solid[i]),
            u8::from(best.physical.dilated.solid[i]),scales[i]);
    }
    if report.stop != BinaryLoadSpaceScanStop3::ScheduleComplete {
        return Err("two-grid binary scan stopped; exported endpoint belongs only to the completed prefix".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn binary_refined_command_admits_inputs_before_geometry_and_runs_real_physics() {
        for args in [vec![],vec!["0"],vec!["NaN"],vec!["1000","bad"],
            vec!["1000","0","250000","NaN"],vec!["1000","0","250000","0.02","1.1"],
            vec!["1000","0","250000","0.02","1","extra"]] {
            assert!(super::run(&args.into_iter().map(String::from).collect::<Vec<_>>()).is_err());
        }
        // Loose requirement and permissive grid-change gate exercise the full
        // gray solve plus exhaustive binary two-grid path without asserting
        // design quality or a manufactured threshold.
        super::run(&["1000000000000","0","250000","1","1"].map(String::from)).unwrap();
    }
}
