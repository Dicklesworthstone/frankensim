//! Real adaptive elasticity consumers, not substituted multigrid operators.
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3,HeightAxis,HexCell};
use fs_cutfem::elastic3::{ElasticityOptions3,adaptive::AdaptiveElasticity3};
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::elastic3::adaptive::enrichment::precondition::{
    AdaptiveMultilevelOptions3,AdaptiveSolveSpace3,AdaptiveSetupWork3,AdaptivePreconditionError3,
};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3,QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_solver::{CgState,LinearOp};
use fs_solver::op::multilevel::MultilevelError;
use fs_sparse::precond::Precond;
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self,p:[f64;3])->f64 {p[2]-0.73}
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval {Interval::new(lo[2],hi[2])-Interval::new(0.73,0.73)}
    fn derivative_enclose(&self,_:[f64;3],_:[f64;3],a:HeightAxis)->Interval {
        let d=if a==HeightAxis::Z{1.0}else{0.0};Interval::new(d,d)
    }
}
fn build(t:&Octree3)->AdaptiveElasticity3 {
    let mut poll=|_|ControlFlow::Continue(());
    let mut q=QuadratureControl3::new(QuadratureOptions3{depth:1,..Default::default()},&mut poll).unwrap();
    AdaptiveElasticity3::build(HexCell::try_new([0.0;3],[1.0;3]).unwrap(),t,&Slab,
        &IsotropicElastic::new(1.0,0.3,1.0).unwrap(),&|p|p[0]==0.0,ElasticityOptions3::default(),&mut q).unwrap()
}
fn tree(level:u8)->Octree3 {Octree3::uniform(level,5,8192).unwrap()}
fn action(p:&impl Precond,n:usize)->Vec<f64> {
    let r:Vec<f64>=(0..n).map(|i|((i*17)%29) as f64/29.0-0.5).collect();let mut z=vec![0.0;n];p.apply(&r,&mut z);z
}
#[test]
fn g0_single_sparse_coarse_level_matches_existing_dense_two_level() {
    let coarse=build(&tree(1));let fine_tree=tree(1).refined(&[*tree(1).leaves().iter().next().unwrap()],||ControlFlow::Continue(())).unwrap();
    let mut space=AdaptiveSolveSpace3::multilevel(build(&fine_tree),&[&coarse],AdaptiveMultilevelOptions3::default(),||ControlFlow::Continue(())).unwrap();
    let scales:Vec<f64>=(0..space.elasticity().cells()).map(|i|0.3+0.1*(i%5) as f64).collect();space.set_scales(&scales).unwrap();
    let prepared=space.prepare_with_work(|_|ControlFlow::Continue(())).unwrap();
    let transfer=AdaptiveTransfer3::new(&coarse,space.elasticity(),4_000_000,||ControlFlow::Continue(())).unwrap();
    let dense=transfer.prepare_two_level(Default::default(),100_000_000,|_|ControlFlow::Continue(())).unwrap();
    let a=action(&prepared,space.n());let b=action(&dense,space.n());let scale=b.iter().map(|v|v.abs()).fold(1e-30_f64,f64::max);
    assert!(a.iter().zip(&b).all(|(a,b)|(a-b).abs()<1e-9*scale));
    assert!(prepared.setup_work().galerkin_products>0);assert_eq!(prepared.work().operator_applications,0);
}
#[test]
fn g1_mixed_cut_grid_solves_with_more_than_512_first_coarse_coordinates() {
    let c3=build(&tree(3));let c2=build(&tree(2));let c1=build(&tree(1));
    let t=tree(3);let fine=t.refined(&[*t.leaves().iter().next().unwrap()],||ControlFlow::Continue(())).unwrap();
    let space=AdaptiveSolveSpace3::multilevel(build(&fine),&[&c3,&c2,&c1],AdaptiveMultilevelOptions3::default(),||ControlFlow::Continue(())).unwrap();
    let sizes=space.level_sizes();assert!(sizes[1]>512);assert!(sizes.last().unwrap()<&192);assert_eq!(sizes.len(),4);
    let p=space.prepare_with_work(|_|ControlFlow::Continue(())).unwrap();
    let rhs=space.elasticity().body_load(&|_|[0.0,0.0,-1.0],||ControlFlow::Continue(())).unwrap();
    let mut cg=CgState::new(&space,&p,&rhs);assert!(cg.run(&space,&p,1e-11,2000).converged);
    assert!(space.elasticity().field_residual(&cg.x,&rhs,||ControlFlow::Continue(())).unwrap()<1e-8);
    assert!(p.setup_work().galerkin_products>0);
}
#[test]
fn g3_density_changes_rebuild_sparse_hierarchy_without_changing_geometry() {
    let c1=build(&tree(1));let c0=build(&tree(0));
    let mut space=AdaptiveSolveSpace3::multilevel(build(&tree(2)),&[&c1,&c0],AdaptiveMultilevelOptions3::default(),||ControlFlow::Continue(())).unwrap();
    let sizes=space.level_sizes();let nodes=space.elasticity().physical_nodes().to_vec();let mut previous=None;
    for stage in 0..3 {
        let scales:Vec<f64>=(0..space.elasticity().cells()).map(|i|0.2+0.1*((i+stage)%7) as f64).collect();space.set_scales(&scales).unwrap();
        let p=space.prepare_with_work(|_|ControlFlow::Continue(())).unwrap();let a=action(&p,space.n());
        let q=space.prepare_with_work(|_|ControlFlow::Continue(())).unwrap();assert_eq!(a,action(&q,space.n()));
        if let Some(old)=previous {assert_ne!(old,a);}previous=Some(a);
        assert_eq!(p.setup_work(),q.setup_work());assert_eq!(space.level_sizes(),sizes);assert_eq!(space.elasticity().physical_nodes(),nodes);
    }
}
#[test]
fn g4_cancelled_and_exhausted_setup_reports_spent_products_without_partial_action() {
    let coarse=build(&tree(0));
    let mut options=AdaptiveMultilevelOptions3::default();options.hierarchy.max_galerkin_products=3;
    let space=AdaptiveSolveSpace3::multilevel(build(&tree(1)),&[&coarse],options,||ControlFlow::Continue(())).unwrap();
    let scales=space.elasticity().scales().to_vec();let mut spent=AdaptiveSetupWork3::default();
    assert!(matches!(space.prepare_with_work(|w|{spent=w;ControlFlow::Continue(())}),Err(AdaptivePreconditionError3::Hierarchy(MultilevelError::Budget(_)))));
    assert_eq!(spent.galerkin_products,3);assert_eq!(spent.operator_applications,0);assert_eq!(space.elasticity().scales(),scales);
    let space=AdaptiveSolveSpace3::multilevel(build(&tree(1)),&[&coarse],AdaptiveMultilevelOptions3::default(),||ControlFlow::Continue(())).unwrap();
    assert!(matches!(space.prepare_with_work(|w|if w.galerkin_products>0{ControlFlow::Break(())}else{ControlFlow::Continue(())}),Err(AdaptivePreconditionError3::Hierarchy(MultilevelError::Cancelled))));
    assert!(space.prepare_with_work(|_|ControlFlow::Continue(())).is_ok());
}
#[test]
fn g0_non_nested_or_truncated_hierarchies_refuse_instead_of_rediscretizing() {
    let c1=build(&tree(1));let c2=build(&tree(2));
    assert!(AdaptiveSolveSpace3::multilevel(build(&tree(1)),&[&c2],AdaptiveMultilevelOptions3::default(),||ControlFlow::Continue(())).is_err());
    let mut options=AdaptiveMultilevelOptions3::default();options.hierarchy.max_coarsest_dofs=1;
    assert!(matches!(AdaptiveSolveSpace3::multilevel(build(&tree(2)),&[&c1],options,||ControlFlow::Continue(())),Err(AdaptivePreconditionError3::Hierarchy(MultilevelError::Budget(_)))));
    assert!(AdaptiveSolveSpace3::multilevel(build(&tree(1)),&[],AdaptiveMultilevelOptions3::default(),||ControlFlow::Continue(())).is_err());
}
