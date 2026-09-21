//! G1/G3/G4/G5: real reference experiments and transactional selective refits.
use std::{cell::Cell, ops::ControlFlow};
use fs_ascent::projected_al::{ProjectedAlError, ProjectedAlStop};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{adaptive::AdaptiveElasticity3, ElasticityOptions3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{SimpParams, SolveBudget, SolveControl, SolveProgress};
use fs_topopt::sdf3::CutDensityStudy3;
use fs_topopt::sdf3::response::{ProjectedResponseOptions3, ProjectedResponseStudy3, ResponseCase3, ResponseTarget3};
use fs_topopt::sdf3::response::refinement::{ReferenceResponseCase3, ReferenceResponseTarget3, ResponseRefinementOptions3};

struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64;3]) -> f64 { (p[0]-0.17)*(p[0]-0.83) }
    fn enclose(&self, lo: [f64;3], hi: [f64;3]) -> Interval {
        let x = Interval::new(lo[0],hi[0]);
        (x-Interval::new(0.17,0.17))*(x-Interval::new(0.83,0.83))
    }
    fn derivative_enclose(&self, lo: [f64;3], hi: [f64;3], axis: HeightAxis) -> Interval {
        if axis == HeightAxis::X { Interval::new(2.0,2.0)*Interval::new(lo[0],hi[0])-Interval::new(1.0,1.0) }
        else { Interval::new(0.0,0.0) }
    }
}
fn tree(level: u8) -> Octree3 { Octree3::uniform(level,4,4096).unwrap() }
fn backend(tree: &Octree3) -> AdaptiveSolveSpace3 {
    let mut poll = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3::default(),&mut poll).unwrap();
    let op = AdaptiveElasticity3::build_with_embedded_dirichlet(
        HexCell::try_new([0.0;3],[1.0;3]).unwrap(),tree,&Slab,&IsotropicElastic::new(1.0,0.3,1.0).unwrap(),
        &|_| false,&|_,_| true,ElasticityOptions3::default(),Default::default(),Default::default(),&mut q).unwrap();
    AdaptiveSolveSpace3::jacobi(op,100_000_000)
}
fn study(tree: &Octree3) -> CutDensityStudy3<AdaptiveSolveSpace3> {
    CutDensityStudy3::new(backend(tree),0.15,SimpParams::default())
}
fn force(_: [f64;3]) -> [f64;3] { [0.002,0.0,-0.003] }
fn observe(p: [f64;3]) -> [f64;3] { [0.0,0.0,1.0+p[0]] }
fn motion(p: [f64;3], n: [f64;3]) -> [f64;3] {
    if n[0]>0.0 { [0.02,0.0,0.0] } else { [0.0,0.0,0.007*p[1]] }
}
fn targets() -> [ReferenceResponseTarget3<'static>;1] {
    [ReferenceResponseTarget3 { observation: ReferenceLoad3::body(&observe), target: 0.003, scale: 0.02, weight: 0.001 }]
}
fn options() -> ProjectedResponseOptions3 {
    let mut o = ProjectedResponseOptions3::default();
    o.density_floor = 0.2; o.response.volume_weight = 1.0; o
}

#[test]
fn g5_reference_binding_matches_the_existing_nodal_fit_exactly() {
    let mut a = study(&tree(1)); let mut b = study(&tree(1)); let ts = targets();
    let cases = [ReferenceResponseCase3 { load: ReferenceLoad3::body(&force), prescribed: Some(&motion), targets: &ts }];
    let f = b.operator().elasticity().reference_load(cases[0].load,||ControlFlow::Continue(())).unwrap();
    let q = b.operator().elasticity().reference_load(ts[0].observation,||ControlFlow::Continue(())).unwrap();
    let ns = [ResponseTarget3 { q: &q, target: ts[0].target, scale: ts[0].scale, weight: ts[0].weight }];
    let nodal = [ResponseCase3 { force: &f, prescribed: Some(&motion), targets: &ns }];
    let rho = vec![0.45;a.cells()];
    let mut pa = |_| ControlFlow::Continue(()); let mut ca = SolveControl::new(SolveBudget::default(),&mut pa);
    let mut pb = |_| ControlFlow::Continue(()); let mut cb = SolveControl::new(SolveBudget::default(),&mut pb);
    let fit = a.fit_reference_responses(&cases,&rho,options(),2,&mut ca).unwrap();
    let mut original = ProjectedResponseStudy3::new(&mut b,&nodal,&rho,options(),&mut cb).unwrap();
    let report = original.run(2).unwrap();
    assert_eq!(fit.outcome.unwrap().stop,report.stop);
    assert_eq!(fit.accepted.rho,original.point());
    assert_eq!(fit.accepted.displacements,original.accepted().displacements);
    assert_eq!(fit.accepted.gradient,original.accepted().gradient);
    assert_eq!(fit.work,original.work());
    assert!(fit.accepted.objective < fit.initial.objective);
}

#[test]
fn g1_goal_marked_refit_reassembles_experiments_and_solves_a_fresh_feasible_baseline() {
    let grid = tree(1); let mut coarse = study(&grid); let ts = targets();
    let cases = [ReferenceResponseCase3 { load: ReferenceLoad3::body(&force), prescribed: Some(&motion), targets: &ts }];
    let mut poll = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(),&mut poll);
    let fit = coarse.fit_reference_responses(&cases,&vec![0.8;coarse.cells()],options(),0,&mut c).unwrap();
    let evidence = coarse.estimate_response_enrichment(backend(&tree(2)),&fit.accepted,&cases,
        ResponseRefinementOptions3 { response: options().response, ..Default::default() },&mut c).unwrap();
    let marked = evidence.mark(0.5,1,||ControlFlow::Continue(())).unwrap();
    assert_eq!(marked.marked.len(),1);
    let fine_grid = grid.refined(&marked.marked,||ControlFlow::Continue(())).unwrap();
    let previous = coarse.operator().elasticity().scales().to_vec(); let work = c.work();
    let next = coarse.refit_reference_responses(study(&fine_grid),&fit.accepted,&cases,options(),
        2,2_000_000,&mut c).unwrap();
    assert!(next.fit.outcome.is_ok());
    assert!(next.fit.history.len()>1);
    assert!(next.inherited_rho.iter().all(|r| *r==0.8));
    assert!(next.starting_rho.iter().all(|r| *r>=0.2 && *r<0.8));
    assert!(next.fit.initial.volume_fraction<=0.5);
    assert!(next.fit.accepted.objective<next.fit.initial.objective);
    assert!(next.study.cells()>coarse.cells());
    assert!(next.fit.initial.displacements[0].len()>fit.accepted.displacements[0].len());
    let op = next.study.operator().elasticity();
    let rhs = op.reference_load_with_motion(cases[0].load,cases[0].prescribed,||ControlFlow::Continue(())).unwrap();
    assert!(op.field_residual(&next.fit.accepted.displacements[0],&rhs,||ControlFlow::Continue(())).unwrap()<1e-8);
    assert_eq!(coarse.operator().elasticity().scales(),previous);
    assert!(c.work().linear_iterations>work.linear_iterations);
}

#[test]
fn g4_interrupted_candidate_retains_its_baseline_and_does_not_replace_source() {
    let mut coarse = study(&tree(1)); let ts = targets();
    let cases = [ReferenceResponseCase3 { load: ReferenceLoad3::body(&force), prescribed: Some(&motion), targets: &ts }];
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(),&mut p);
    let fit = coarse.fit_reference_responses(&cases,&vec![0.45;coarse.cells()],options(),0,&mut c).unwrap();
    let scales = coarse.operator().elasticity().scales().to_vec(); let mut evaluations = 0;
    let mut p = |s: SolveProgress| {
        if s.stage=="response-projected-evaluate" { evaluations+=1; }
        if evaluations==2 { ControlFlow::Break(()) } else { ControlFlow::Continue(()) }
    };
    let mut c = SolveControl::new(SolveBudget::default(),&mut p);
    let next = coarse.refit_reference_responses(study(&tree(2)),&fit.accepted,&cases,options(),3,2_000_000,&mut c).unwrap();
    assert!(next.fit.outcome.is_err());
    assert_eq!(next.fit.accepted.rho,next.fit.initial.rho);
    assert_eq!(next.fit.accepted.displacements,next.fit.initial.displacements);
    assert!(next.fit.work.linear_iterations>0);
    assert_eq!(coarse.operator().elasticity().scales(),scales);
}

#[test]
fn g0_policy_and_grid_refusals_do_not_invoke_unadmitted_load_laws() {
    let mut source = study(&tree(1)); let calls = Cell::new(0);
    let law = |p| { calls.set(calls.get()+1); force(p) };
    let mut ts = targets(); ts[0].scale = f64::NAN;
    let cases = [ReferenceResponseCase3 { load: ReferenceLoad3::body(&law), prescribed: Some(&motion), targets: &ts }];
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(),&mut p);
    assert!(source.fit_reference_responses(&cases,&vec![0.45;source.cells()],options(),1,&mut c).is_err());
    assert_eq!(calls.get(),0); assert_eq!(c.work().linear_solves,0);
    let ts = targets(); let cases = [ReferenceResponseCase3 { load: ReferenceLoad3::body(&law), prescribed: Some(&motion), targets: &ts }];
    let fit = source.fit_reference_responses(&cases,&vec![0.45;source.cells()],options(),0,&mut c).unwrap();
    let count = calls.get(); let scales = source.operator().elasticity().scales().to_vec();
    assert!(source.refit_reference_responses(study(&tree(1)),&fit.accepted,&cases,options(),1,2_000_000,&mut c).is_err());
    assert_eq!(calls.get(),count);
    let mut stale = fit.accepted.clone(); stale.rho[0] = 0.7;
    assert!(source.refit_reference_responses(study(&tree(2)),&stale,&cases,options(),1,2_000_000,&mut c).is_err());
    assert_eq!(source.operator().elasticity().scales(),scales);
}

#[test]
fn g4_unfunded_new_grid_cannot_replace_source_or_claim_optimizer_convergence() {
    let mut source = study(&tree(1)); let ts = targets();
    let cases = [ReferenceResponseCase3 { load: ReferenceLoad3::body(&force), prescribed: Some(&motion), targets: &ts }];
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(),&mut p);
    let fit = source.fit_reference_responses(&cases,&vec![0.45;source.cells()],options(),0,&mut c).unwrap();
    let scales = source.operator().elasticity().scales().to_vec();
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget { total_iterations: 1, ..Default::default() },&mut p);
    assert!(matches!(source.refit_reference_responses(study(&tree(2)),&fit.accepted,&cases,options(),3,2_000_000,&mut c),Err(ProjectedAlError::Evaluation(_))));
    assert_eq!(c.work().linear_iterations,1); assert_eq!(source.operator().elasticity().scales(),scales);
    let mut o = options(); o.optimizer.max_evaluations = 1;
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(),&mut p);
    let limited = source.fit_reference_responses(&cases,&vec![0.8;source.cells()],o,100,&mut c).unwrap();
    let report = limited.outcome.unwrap();
    assert_eq!(report.stop,ProjectedAlStop::EvaluationLimit);
    assert!(!report.kkt.within_tolerance(o.optimizer.tolerance));
    assert!(limited.accepted.volume_fraction>o.volume_cap);
}
