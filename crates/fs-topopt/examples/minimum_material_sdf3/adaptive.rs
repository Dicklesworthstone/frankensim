//! One limit-driven estimate/mark/refit using the existing minimum-material API.
use super::*;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_topopt::sdf3::response::limits::refinement::*;

type Study = CutDensityStudy3<AdaptiveSolveSpace3>;
type Error = Box<dyn std::error::Error>;
fn observe(_: [f64;3])->[f64;3] {[1.0,0.0,0.0]}
fn build(tree:&Octree3, geometry:&mut QuadratureControl3<'_>)->Result<AdaptiveSolveSpace3,Error> {
    let op=AdaptiveElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0;3],[1.0;3])?,tree,&Slab,
        &IsotropicElastic::new(1.0,0.0,1.0)?,&|_|false,&|_,_|true,
        ElasticityOptions3::default(),Default::default(),Default::default(),geometry)?;
    Ok(AdaptiveSolveSpace3::jacobi(op,100_000_000))
}
#[allow(clippy::too_many_arguments)]
pub(super) fn run(study:&mut Study, tree:&Octree3, accepted:&ResponseLimitEvaluation3,
    cases:&[ResponseLimitCase3<'_>], options:MinimumMaterialOptions3, steps:usize, refit:bool,
    geometry:&mut QuadratureControl3<'_>, control:&mut SolveControl<'_>)->Result<(),Error> {
    // These reference laws are exactly this fixture's zero external force and
    // volume-integral observation, not a generic reconstruction from nodal q.
    let rows:Vec<Vec<_>>=cases.iter().map(|case|case.limits.iter().map(|limit|ReferenceResponseLimit3 {
        sense:limit.sense,bound:limit.bound,scale:limit.scale,response:match limit.response {
            BoundedResponse3::Displacement(_)=>ReferenceBoundedResponse3::Displacement(ReferenceLoad3::body(&observe)),
            BoundedResponse3::Reaction(mode)=>ReferenceBoundedResponse3::Reaction(mode),
        },
    }).collect()).collect();
    let laws:Vec<_>=cases.iter().zip(&rows).map(|(case,limits)|ReferenceResponseLimitCase3 {
        load:ReferenceLoad3::default(),prescribed:case.prescribed,limits,
    }).collect();
    let probe=tree.refined(study.operator().elasticity().leaves(),||ControlFlow::Continue(()))?;
    let evidence=study.estimate_response_limit_enrichment(build(&probe,geometry)?,accepted,&laws,
        ResponseLimitRefinementOptions3 {limits:options.limits,objective:options.objective,
            volume_scale:options.volume_scale,..Default::default()},control)?;
    println!("estimate_case,limit,sense,coarse_value,enriched_value,bound,scale,coarse_residual,enriched_residual,reconstructed_change,reaction_offset_transfer");
    for (i,case) in evidence.cases.iter().enumerate() {for (j,l) in case.limits.iter().enumerate() {
        println!("{i},{j},{:?},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
            l.sense,l.goal.coarse_value,l.goal.fine_value,l.bound,l.scale,l.coarse_residual,l.fine_residual,l.correction(),l.goal.offset_transfer);
    }}
    let marks=evidence.mark(0.5,2,||ControlFlow::Continue(()))?;
    let feasible=evidence.fine_feasible(options.optimizer.tolerance);
    eprintln!("phase=limit-estimate; inherited_physical_stiffness=true; fine_feasible={feasible}; smooth_constraint={:.17e}; worst_original={:.17e}; maximum_response_change={:.17e}; marked={:?}; achieved_fraction={:.6}; marking_target_met={}; work={:?}; continuum_certified=false",
        evidence.fine_constraint,evidence.fine_worst_residual,evidence.max_normalized_change(),marks.marked,marks.achieved_fraction,marks.target_met,evidence.work);
    if !refit {
        if !feasible {return Err("enriched physical design violates a declared numerical response limit".into());}
        return Ok(());
    }
    if marks.marked.is_empty() {return Err("no localized refinement signal; this is not a continuum feasibility certificate".into());}
    let marked_tree=tree.refined(&marks.marked,||ControlFlow::Continue(()))?;
    let candidate=CutDensityStudy3::new(build(&marked_tree,geometry)?,0.15,study.params());
    let next=study.refit_reference_minimum_material(candidate,accepted,&laws,options,steps,2_000_000,control)?;
    println!("refit_iteration,volume_fraction,smooth_constraint,worst_original_residual");
    for row in &next.fit.history {println!("{},{:.17e},{:.17e},{:.17e}",row.iteration,row.volume_fraction,row.constraint,row.worst_residual);}
    // Export the retained candidate before propagating a late AL/cancellation
    // error. The source stays installed until a feasible endpoint is selected.
    let selected=next.fit.best_feasible.as_ref().unwrap_or(&next.fit.accepted);
    println!("refit_case,limit,value,bound,scale,original_residual");
    for (i,case) in laws.iter().enumerate() {for (j,l) in case.limits.iter().enumerate() {
        println!("{i},{j},{:.17e},{:.17e},{:.17e},{:.17e}",selected.measurements[i][j],l.bound,l.scale,selected.residuals[i][j]);
    }}
    println!("refit_cell,inherited_raw_density,selected_raw_density,projected_density");
    for (i,raw) in selected.rho.iter().enumerate() {println!("{i},{:.17e},{raw:.17e},{:.17e}",next.inherited_rho[i],selected.projected_rho[i]);}
    eprintln!("phase=marked-refit; baseline_feasible={}; selected_feasible={}; selected_latest={}; source_volume={:.17e}; new_baseline_volume={:.17e}; selected_volume={:.17e}; work={:?}; cross_grid_descent_claimed=false",
        next.fit.initial.feasible(options.optimizer.tolerance),selected.feasible(options.optimizer.tolerance),
        selected.rho==next.fit.accepted.rho,accepted.volume_fraction,next.fit.initial.volume_fraction,selected.volume_fraction,next.fit.work);
    let report=next.fit.outcome.clone()?;
    eprintln!("latest_al_stop={:?}; latest_al_kkt={:?}; final_design_needs_its_own_error_probe=true",report.stop,report.kkt);
    let (candidate,_selected)=next.into_best_feasible(control)?;
    *study=candidate;
    Ok(())
}
