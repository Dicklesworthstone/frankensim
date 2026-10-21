//! Assess -> refit the complete load-space cap -> assess the NEW material.
#[path = "refinement/adaptive.rs"]
mod adaptive;
pub fn run_adaptive(args: &[String]) -> CommandResult { adaptive::run(args) }
use super::*;
use fs_cutfem::elastic3::{adaptive::{AdaptiveElasticity3,
    enrichment::precondition::AdaptiveSolveSpace3}, surface::ReferenceLoad3};
use fs_cutfem::octree3::Octree3;
use fs_topopt::sdf3::robust::envelope::load_space::design::limits::{LoadSpaceCap3,
    refit::{LoadSpaceCapRefit3, LoadSpaceCapRefitOptions3}, refinement::LoadSpaceCapRefinement3};
use fs_topopt::sdf3_goal::{GoalPreconditioner3, GoalRefinementError3, GoalRefinementOptions3};

type Design = LoadSpaceCapRefit3<AdaptiveSolveSpace3>;
type CommandResult = Result<(), Box<dyn std::error::Error>>;
const DIAGONAL_LIMIT: usize = 100_000_000;

pub fn run(args: &[String]) -> CommandResult {
    if !(1..=4).contains(&args.len()) {
        return Err("usage: worst_load_sdf3 --load-space --refine-material CAP [UPDATES [TOTAL_KRYLOV_ITERATIONS [GRAM_CHANGE_TOLERANCE]]]".into());
    }
    let cap = LoadSpaceCap3::new(2, args[0].parse()?, 0.03)?;
    let updates = args.get(1).map_or(Ok(200), |s| s.parse::<usize>())?;
    let iterations = args.get(2).map_or(Ok(250_000), |s| s.parse::<usize>())?;
    let change = args.get(3).map_or(Ok(0.02), |s| s.parse::<f64>())?;
    if !change.is_finite() || change <= 0.0 {
        return Err("Gram-change tolerance must be finite and positive".into());
    }
    let params = SimpParams { penal:2.0, beta:2.0, e_min:1e-3, ..Default::default() };
    let options = LoadSpaceCapRefitOptions3 { accepted_steps:updates, ..Default::default() };
    let tree = Octree3::uniform(1,4,4096)?;
    let y = |p: [f64;3]| [0.0,-1.0-p[1],0.0];
    let yz = |p: [f64;3]| [0.0,-0.6*(1.0+p[1]),-1.0-p[0]];
    let laws = [ReferenceLoad3::body(&y),ReferenceLoad3::body(&yz)];
    let limits = QuadratureOptions3::default();
    let mut geometry_poll = |_| ControlFlow::Continue(());
    // One allowance shared by initial, probe and refit builds, including failures.
    let mut geometry = QuadratureControl3::new(limits,&mut geometry_poll)?;
    let domain = HexCell::try_new([0.0;3],[1.0;3])?;
    let material = IsotropicElastic::new(1.0,0.3,1.0)?;
    let mut build = |tree: &Octree3| -> Result<AdaptiveElasticity3,GoalRefinementError3> {
        Ok(AdaptiveElasticity3::build(domain,tree,&Slab,&material,&|p| p[0]==0.0,
            ElasticityOptions3::default(),&mut geometry)?)
    };
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations:iterations,..Default::default() },&mut poll);
    eprintln!("units=dimensionless; load_set=Fq; coefficient_norm_limit=1; cap={}; params={params:?}; refit_policy={options:?}; gram_change_tolerance={change}; max_refits=1; refine_all_active_cells=true; filter_radius=0.15; total_krylov_budget={iterations}; geometry_limits={limits:?}; continuum_certificate=false",
        cap.cap());
    let result = execute(&tree,&laws,&cap,params,options,change,&mut control,&mut build);
    drop(build);
    eprintln!("total_physical_work={:?}; total_geometry_work={:?}",control.work(),geometry.work());
    result
}

#[allow(clippy::too_many_arguments)]
fn execute(tree: &Octree3, laws: &[ReferenceLoad3<'_>], cap: &LoadSpaceCap3,
    params: SimpParams, options: LoadSpaceCapRefitOptions3, change: f64,
    control: &mut SolveControl<'_>,
    build: &mut impl FnMut(&Octree3) -> Result<AdaptiveElasticity3,GoalRefinementError3>) -> CommandResult {
    control.checkpoint("sdf3-load-space-command-build")?;
    let study = CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(build(tree)?,DIAGONAL_LIMIT),0.15,params);
    let raw = vec![0.3;study.cells()];
    let mut coarse = study.solve_reference_load_space_cap(&raw,laws,cap,options,control)?;
    export(&coarse,0);
    if !coarse.completed_and_rechecked() {
        return Err("initial combined-load design stopped; exported fields are only its solved prefix".into());
    }
    let coarse_rho = coarse.recheck.as_ref().unwrap().as_ref().unwrap().rho.clone();
    let next_tree = tree.refined(coarse.study.operator().elasticity().leaves(),||ControlFlow::Continue(()))?;
    let goal = GoalRefinementOptions3 {
        preconditioner:GoalPreconditioner3::Jacobi { max_contributions:DIAGONAL_LIMIT },..Default::default()
    };
    control.checkpoint("sdf3-load-space-command-build")?;
    let first = coarse.study.estimate_load_space_cap_enrichment(build(&next_tree)?,laws,&coarse_rho,
        params,options.design.robustness,cap,goal,control)?;
    let passed = assessment(&first,0,options.design.optimizer.tolerance,change);
    if passed { return Ok(()); }
    // The full Gram check is global. Basis-local DWR masses are not falsely
    // presented as error localization for the complete ellipsoid.
    control.checkpoint("sdf3-load-space-command-build")?;
    let mut refined = coarse.study.refit_load_space_cap(
        AdaptiveSolveSpace3::jacobi(build(&next_tree)?,DIAGONAL_LIMIT),&coarse_rho,
        params,laws,cap,options,control)?;
    export(&refined,1);
    if !refined.completed_and_rechecked() {
        return Err("refined combined-load design stopped; the previous design and candidate prefix were exported".into());
    }
    let refined_rho = refined.recheck.as_ref().unwrap().as_ref().unwrap().rho.clone();
    let probe_tree = next_tree.refined(refined.study.operator().elasticity().leaves(),||ControlFlow::Continue(()))?;
    control.checkpoint("sdf3-load-space-command-build")?;
    let second = refined.study.estimate_load_space_cap_enrichment(build(&probe_tree)?,laws,&refined_rho,
        params,options.design.robustness,cap,goal,control)?;
    if !assessment(&second,1,options.design.optimizer.tolerance,change) {
        return Err("refitted design did not pass the full two-grid cap/change check; no tolerance or cap was relaxed".into());
    }
    Ok(())
}

fn assessment(a: &LoadSpaceCapRefinement3, grid: usize, smooth: f64, change: f64) -> bool {
    let passed = a.stable_and_feasible(smooth,change);
    println!("assessed_grid,coarse_worst_compliance,fine_worst_compliance,cap,normalized_full_gram_change,passed");
    println!("{grid},{:.17e},{:.17e},{:.17e},{:.17e},{passed}",
        a.coarse.worst_compliance,a.fine.worst_compliance,a.cap.cap(),a.normalized_gram_change);
    eprintln!("assessed_grid={grid}; coarse_spectrum={:?}; fine_spectrum={:?}; coarse_worst_coefficients={:?}; fine_worst_coefficients={:?}; probe_inherits_physical_stiffness=true; fine_raw_refiltered=false; change_is_continuum_bound=false",
        a.coarse.spectrum.eigenvalues,a.fine.spectrum.eigenvalues,a.coarse.worst_coefficients,a.fine.worst_coefficients);
    passed
}
fn export(result: &Design, grid: usize) {
    eprintln!("grid={grid}; restoration_stop={:?}; restoration_error={:?}; blend={:?}; optimizer_outcome={:?}; recheck_error={:?}; completed_and_rechecked={}; work={:?}; cross_grid_descent_claimed=false",
        result.restoration.stop,result.restoration.error,result.restoration.blend,result.outcome,
        result.recheck.as_ref().and_then(|r|r.as_ref().err()),result.completed_and_rechecked(),result.work);
    println!("grid,iteration,dilated_volume,smooth_normalized_maximum,actual_normalized_maximum,feasible");
    for row in &result.history {
        println!("{grid},{},{:.17e},{:.17e},{:.17e},{}",row.iteration,row.dilated_volume,row.smooth_maximum,row.maximum,row.feasible);
    }
    let checked = result.recheck.as_ref().and_then(|r|r.as_ref().ok());
    let Some(point) = checked.or(result.best.as_ref()).or(result.restoration.point.as_ref()) else {
        eprintln!("grid={grid}; exported_design=none"); return;
    };
    let source = if checked.is_some() { "rechecked-incumbent" }
        else if result.best.is_some() { "accepted-incumbent" } else { "restoration-seed-not-optimized" };
    eprintln!("grid={grid}; exported_source={source}; physical_cap_satisfied={}",
        result.cap.admits(&point.evaluation,result.options.design.optimizer.tolerance));
    let op = result.study.operator().elasticity(); let e = &point.evaluation;
    println!("grid,level,cell_x,cell_y,cell_z,raw_density,eroded_density,nominal_density,dilated_density");
    for (i,key) in op.leaves().iter().enumerate() {
        let [x,y,z] = key.index();
        println!("{grid},{},{x},{y},{z},{:.17e},{:.17e},{:.17e},{:.17e}",key.level(),point.rho[i],
            e.physical.eroded.projected_rho[i],e.physical.nominal_rho[i],e.physical.dilated_rho[i]);
    }
    println!("grid,basis_index,worst_coefficient");
    for (i,q) in e.worst_coefficients.iter().enumerate() { println!("{grid},{i},{q:.17e}"); }
    println!("grid,node_x,node_y,node_z,worst_fx,worst_fy,worst_fz,worst_ux,worst_uy,worst_uz");
    for (i,p) in op.nodes().iter().enumerate() {
        let f = &e.worst_force; let u = &e.worst_displacement;
        println!("{grid},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            p[0],p[1],p[2],f[3*i],f[3*i+1],f[3*i+2],u[3*i],u[3*i+1],u[3*i+2]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn command_admits_authored_inputs_before_building_and_can_pass_a_real_probe() {
        for args in [vec![],vec!["0"],vec!["NaN"],vec!["1000","bad"],
            vec!["1000","0","250000","-1"],vec!["1000","0","250000","NaN"],
            vec!["1000","0","250000","0.02","extra"]] {
            assert!(run(&args.into_iter().map(String::from).collect::<Vec<_>>()).is_err());
        }
        // Loose authored requirement isolates command plumbing and real physics.
        super::super::run(&["--refine-material","1000000000000","0","250000","1"]
            .map(String::from)).unwrap();
    }
    #[test]
    fn failed_first_probe_refits_then_checks_the_new_material_without_relaxation() {
        let tree = Octree3::uniform(1,4,4096).unwrap();
        let cap = LoadSpaceCap3::new(2,1e12,0.03).unwrap();
        let params = SimpParams { penal:2.0,beta:2.0,e_min:1e-3,..Default::default() };
        let y = |_: [f64;3]| [0.0,-1.0,0.0]; let z = |_: [f64;3]| [0.0,-0.4,-1.0];
        let laws = [ReferenceLoad3::body(&y),ReferenceLoad3::body(&z)];
        let mut builds = 0;
        let mut build = |tree: &Octree3| -> Result<AdaptiveElasticity3,GoalRefinementError3> {
            builds += 1;
            let mut poll = |_| ControlFlow::Continue(());
            let mut q = QuadratureControl3::new(QuadratureOptions3 {depth:1,..Default::default()},&mut poll).unwrap();
            Ok(AdaptiveElasticity3::build(HexCell::try_new([0.0;3],[1.0;3]).unwrap(),tree,&Slab,
                &IsotropicElastic::new(1.0,0.3,1.0).unwrap(),&|p|p[0]==0.0,
                ElasticityOptions3::default(),&mut q)?)
        };
        let mut poll = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(),&mut poll);
        // Zero internal change allowance forces both nonidentical-grid probes
        // to be assessed. The public command requires a positive allowance.
        let result = execute(&tree,&laws,&cap,params,
            LoadSpaceCapRefitOptions3 {accepted_steps:0,..Default::default()},0.0,&mut control,&mut build);
        assert_eq!(result.unwrap_err().to_string(),
            "refitted design did not pass the full two-grid cap/change check; no tolerance or cap was relaxed");
        assert_eq!(builds,4); // initial, diagnostic, actual refit, new diagnostic
        assert!(control.work().linear_solves > 0);
    }
}
