//! Selective response refitting; no second optimizer or boundary/goal integrator.
use super::*;
use fs_topopt::sdf3::response::ResponseEvaluation3;
use fs_topopt::sdf3::response::refinement::ResponseRefinement3;

type Study = CutDensityStudy3<AdaptiveSolveSpace3>;
type Error = Box<dyn std::error::Error>;

pub(super) fn validate(rounds: usize, projected: bool) -> Result<(), Error> {
    if !(2..=4).contains(&rounds) || !projected {
        return Err("--adapt-rounds requires 2..=4 total grids and explicit --projected".into());
    }
    Ok(())
}
fn build(tree: &Octree3, geometry: &mut QuadratureControl3<'_>) -> Result<AdaptiveSolveSpace3, Error> {
    let op = AdaptiveElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0;3],[1.0;3])?,
        tree,&Slab,&IsotropicElastic::new(1.0,0.3,1.0)?,&|_| false,&|_,_| true,
        ElasticityOptions3::default(),Default::default(),Default::default(),geometry)?;
    Ok(AdaptiveSolveSpace3::jacobi(op,100_000_000))
}
fn estimate(study: &mut Study, tree: &Octree3, accepted: &ResponseEvaluation3,
    cases: &[ReferenceResponseCase3<'_>], options: ProjectedResponseOptions3,
    geometry: &mut QuadratureControl3<'_>, control: &mut SolveControl<'_>) -> Result<ResponseRefinement3, Error> {
    // The global refinement is ONLY an estimation probe. The retained grid
    // below is refined from the original partition and the marked coarse keys.
    let probe = tree.refined(&tree.leaves().iter().copied().collect::<Vec<_>>(),||ControlFlow::Continue(()))?;
    let evidence = study.estimate_response_enrichment(build(&probe,geometry)?,accepted,cases,
        ResponseRefinementOptions3 { response: options.response, ..Default::default() },control)?;
    Ok(evidence)
}
fn retained(round: usize, study: &Study, accepted: &ResponseEvaluation3, values: &[f64]) {
    println!("retained_round,cell_level,cell_i,cell_j,cell_k,raw_density,projected_density");
    for ((leaf,rho),physical) in study.operator().elasticity().leaves().iter().zip(&accepted.rho).zip(&accepted.projected_rho) {
        let [i,j,k] = leaf.index();
        println!("{round},{},{i},{j},{k},{rho:.17e},{physical:.17e}",leaf.level());
    }
    println!("retained_round,case,observation,response,target");
    for (i,case) in accepted.responses.iter().enumerate() { for (j,value) in case.iter().enumerate() {
        println!("{round},{i},{j},{value:.17e},{:.17e}",values[2*i+j]);
    } }
}

// All physical laws and target scalars stay fixed across the entire campaign.
// The same SolveControl and QuadratureControl are never reset between rounds.
// Error paths print the LAST promoted design, not a rejected fine candidate.
#[allow(clippy::too_many_arguments)]
pub(super) fn run(study: &mut Study, mut accepted: ResponseEvaluation3, level: u8, steps: usize,
    rounds: usize, values: &[f64], final_estimate: bool, options: ProjectedResponseOptions3,
    geometry: &mut QuadratureControl3<'_>, control: &mut SolveControl<'_>) -> Result<(), Error> {
    validate(rounds,true)?;
    if values.len()!=4 || values.iter().any(|v| !v.is_finite()) { return Err("four finite frozen response targets required".into()); }
    let f = |_:[f64;3]| [0.002,0.0,-0.003];
    let f_other = |p| f(p).map(|v| -0.5*v);
    let observe_q = |p:[f64;3]| [0.0,0.0,1.0+p[0]];
    let observe_r = |p:[f64;3]| [1.0+p[1],0.0,0.0];
    let t0 = [
        ReferenceResponseTarget3 { observation: ReferenceLoad3::body(&observe_q), target: values[0], scale: 0.02, weight: 0.7 },
        ReferenceResponseTarget3 { observation: ReferenceLoad3::body(&observe_r), target: values[1], scale: 0.02, weight: 0.3 },
    ];
    let t1 = [ReferenceResponseTarget3 { target: values[2], ..t0[0] },ReferenceResponseTarget3 { target: values[3], ..t0[1] }];
    let cases = [
        ReferenceResponseCase3 { load: ReferenceLoad3::body(&f), prescribed: Some(&motion), targets: &t0 },
        ReferenceResponseCase3 { load: ReferenceLoad3::body(&f_other), prescribed: Some(&other), targets: &t1 },
    ];
    let mut tree = Octree3::uniform(level,4,4096)?;
    let mut promoted_round = 0;
    let result = (|| -> Result<(), Error> {
        for round in 1..rounds {
            let evidence = estimate(study,&tree,&accepted,&cases,options,geometry,control)?;
            let marks = evidence.mark(0.5,2,||ControlFlow::Continue(()))?;
            eprintln!("source_round={}; coarse_objective={:.17e}; probe_objective={:.17e}; two_grid_change={:.17e}; marking_fraction={:.6}; target_met={}; marked={:?}",
                promoted_round,evidence.coarse_objective,evidence.fine_objective,evidence.correction(),marks.achieved_fraction,marks.target_met,marks.marked);
            if marks.marked.is_empty() {
                eprintln!("stop=no_refinement_signal; accuracy_certified=false");
                break;
            }
            let refined_tree = tree.refined(&marks.marked,||ControlFlow::Continue(()))?;
            let candidate = CutDensityStudy3::new(build(&refined_tree,geometry)?,0.15,study.params());
            let next = study.refit_reference_responses(candidate,&accepted,&cases,options,steps,2_000_000,control)?;
            let violation = (next.fit.accepted.volume_fraction-options.volume_cap).max(0.0);
            match &next.fit.outcome {
                Err(error) => return Err(format!("candidate_round={round} rejected; source retained; {error}").into()),
                Ok(report) => eprintln!("candidate_round={round}; stop={:?}; numerical_projected_kkt={:?}; evaluations={}; volume_violation={violation:.9e}",
                    report.stop,report.kkt,report.work.evaluations),
            }
            if violation>options.optimizer.tolerance {
                return Err(format!("candidate_round={round} remains infeasible; source retained").into());
            }
            // A new-grid baseline is not an optimization step and its changed
            // objective is not counted as descent from the old discretization.
            eprintln!("promoting_round={round}; source_objective={:.17e}; new_grid_baseline={:.17e}; refitted_objective={:.17e}; background_cells={}; active_cells={}; cross_grid_descent_claimed=false",
                accepted.objective,next.fit.initial.objective,next.fit.accepted.objective,refined_tree.leaves().len(),next.study.cells());
            println!("round,iteration,objective,volume_fraction,constraint_violation");
            for row in &next.fit.history {
                println!("{round},{},{:.17e},{:.17e},{:.17e}",row.iteration,row.objective,row.volume_fraction,row.constraint_violation);
            }
            *study = next.study;
            accepted = next.fit.accepted;
            tree = refined_tree;
            promoted_round = round;
        }
        if final_estimate && promoted_round+1==rounds {
            let final_check = estimate(study,&tree,&accepted,&cases,options,geometry,control)?;
            eprintln!("final_round={promoted_round}; retained_objective={:.17e}; probe_objective={:.17e}; final_two_grid_change={:.17e}; continuum_certified=false",
                final_check.coarse_objective,final_check.fine_objective,final_check.correction());
        }
        Ok(())
    })();
    retained(promoted_round,study,&accepted,values);
    let work = control.work();
    eprintln!("retained_round={promoted_round}; schedule_complete={}; linear_iterations={}; setup_applications={}; galerkin_products={}; geometry_points={}; continuum_certified=false",
        result.is_ok() && promoted_round+1==rounds,work.linear_iterations,work.preconditioner_operator_applications,
        work.preconditioner_galerkin_products,geometry.work().points);
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g0_adaptive_mode_requires_explicit_optimizer_and_bounded_round_count() {
        assert!(validate(2,true).is_ok()); assert!(validate(4,true).is_ok());
        for n in [0,1,5,usize::MAX] { assert!(validate(n,true).is_err()); }
        assert!(validate(2,false).is_err());
    }

    #[test]
    fn g4_estimation_exhaustion_leaves_the_initial_physical_model_untouched() {
        let mut gp = |_| ControlFlow::Continue(());
        let mut geometry = QuadratureControl3::new(QuadratureOptions3::default(),&mut gp).unwrap();
        let mut study = CutDensityStudy3::new(build(&Octree3::uniform(1,4,4096).unwrap(),&mut geometry).unwrap(),0.15,SimpParams::default());
        let op = study.operator().elasticity();
        let force = op.body_load(&|_| [0.002,0.0,-0.003],||ControlFlow::Continue(())).unwrap();
        let other_force: Vec<_> = force.iter().map(|v| -0.5*v).collect();
        let q = op.body_load(&|p| [0.0,0.0,1.0+p[0]],||ControlFlow::Continue(())).unwrap();
        let r = op.body_load(&|p| [1.0+p[1],0.0,0.0],||ControlFlow::Continue(())).unwrap();
        let targets = [ResponseTarget3 { q: &q, target: 0.0, scale: 0.02, weight: 0.7 },ResponseTarget3 { q: &r, target: 0.0, scale: 0.02, weight: 0.3 }];
        let cases = [ResponseCase3 { force: &force, prescribed: Some(&motion), targets: &targets },ResponseCase3 { force: &other_force, prescribed: Some(&other), targets: &targets }];
        let mut cp = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(),&mut cp);
        let rho = vec![0.45;study.cells()];
        let session = ProjectedResponseStudy3::new(&mut study,&cases,&rho,Default::default(),&mut control).unwrap();
        let accepted = session.accepted().clone(); drop(session);
        let scales = study.operator().elasticity().scales().to_vec();
        let leaves = study.operator().elasticity().leaves().to_vec();
        let mut cp = |_| ControlFlow::Continue(());
        let mut limited = SolveControl::new(SolveBudget { total_iterations: 1, ..Default::default() },&mut cp);
        assert!(run(&mut study,accepted,1,2,2,&[0.0;4],false,Default::default(),&mut geometry,&mut limited).is_err());
        assert_eq!(study.operator().elasticity().scales(),scales);
        assert_eq!(study.operator().elasticity().leaves(),leaves);
        assert_eq!(limited.work().linear_iterations,1);
    }
}
