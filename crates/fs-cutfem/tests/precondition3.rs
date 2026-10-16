use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{ElasticityOptions3, adaptive::{AdaptiveElasticity3, enrichment::AdaptiveTransfer3}};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptivePreconditionError3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_solver::{CgState, LinearOp};
use fs_solver::op::two_level::{TwoLevelBudget, TwoLevelError};
use fs_sparse::precond::{IdentityPrecond, Precond};
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self,p:[f64;3])->f64 {p[2]-0.73}
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval {Interval::new(lo[2],hi[2])-Interval::new(0.73,0.73)}
    fn derivative_enclose(&self,_:[f64;3],_:[f64;3],a:HeightAxis)->Interval {
        let d=if a==HeightAxis::Z {1.0}else{0.0};Interval::new(d,d)
    }
}
fn build(tree:&Octree3)->AdaptiveElasticity3 {
    let mut p=|_|ControlFlow::Continue(());
    let mut q=QuadratureControl3::new(QuadratureOptions3{depth:1,..Default::default()},&mut p).unwrap();
    AdaptiveElasticity3::build(HexCell::try_new([0.0;3],[1.0;3]).unwrap(),tree,&Slab,
        &IsotropicElastic::new(1.0,0.3,1.0).unwrap(),&|p|p[0]==0.0,ElasticityOptions3::default(),&mut q).unwrap()
}
fn mixed_tree()->Octree3 {
    let tree=Octree3::uniform(1,4,4096).unwrap();
    let mark=*tree.leaves().iter().find(|c|c.index()==[1,0,0]).unwrap();
    tree.refined(&[mark],||ControlFlow::Continue(())).unwrap()
}
#[test]
fn constrained_diagonal_matches_full_operator_probing_and_tracks_density() {
    let mut op=build(&mixed_tree());
    for pass in 0..2 {
        if pass==1 {op.set_scales(&vec![0.2;op.cells()]).unwrap();}
        let jacobi=op.prepare_jacobi(20_000_000,||ControlFlow::Continue(())).unwrap();
        let mut basis=vec![0.0;op.n()];let mut column=basis.clone();
        for i in 0..op.n() {
            basis.fill(0.0);basis[i]=1.0;op.apply(&basis,&mut column);
            assert!((1.0/jacobi.inverse_diagonal()[i]-column[i]).abs()<1e-12*column[i].abs().max(1.0));
        }
        assert!(jacobi.contributions()>op.n());
    }
}
#[test]
fn geometric_two_level_reduces_work_on_a_real_mixed_level_cut_domain() {
    let tree=mixed_tree();let coarse=build(&tree);
    let fine_tree=tree.refined(&tree.leaves().iter().copied().collect::<Vec<_>>(),||ControlFlow::Continue(())).unwrap();
    let fine=build(&fine_tree);
    let transfer=AdaptiveTransfer3::new(&coarse,&fine,100_000,||ControlFlow::Continue(())).unwrap();
    let pc=transfer.prepare_two_level(TwoLevelBudget::default(),20_000_000,|_|ControlFlow::Continue(())).unwrap();
    let rhs=fine.body_load(&|_|[0.0,0.0,-1.0],||ControlFlow::Continue(())).unwrap();
    let mut plain=CgState::new(&fine,&IdentityPrecond,&rhs);
    assert!(plain.run(&fine,&IdentityPrecond,1e-12,2000).converged);
    let mut accelerated=CgState::new(&fine,&pc,&rhs);
    assert!(accelerated.run(&fine,&pc,1e-12,2000).converged);
    assert!(accelerated.iters<plain.iters,"{} vs {}",accelerated.iters,plain.iters);
    assert!(fine.field_residual(&accelerated.x,&rhs,||ControlFlow::Continue(())).unwrap()<1e-9);
    let c0:f64=plain.x.iter().zip(&rhs).map(|(u,f)|u*f).sum();
    let c1:f64=accelerated.x.iter().zip(&rhs).map(|(u,f)|u*f).sum();
    assert!((c0-c1).abs()<1e-8*c0.abs());
    assert_eq!(pc.work().operator_applications,pc.coarse_dofs());
}
#[test]
fn coarse_geometry_supplies_a_space_not_a_substituted_stiffness() {
    let tree=Octree3::uniform(1,3,512).unwrap();let mut coarse=build(&tree);
    let fine_tree=tree.refined(&tree.leaves().iter().copied().collect::<Vec<_>>(),||ControlFlow::Continue(())).unwrap();
    let fine=build(&fine_tree);let r=vec![1.0;fine.n()];let mut outputs=Vec::new();
    for scale in [1.0,0.1] {
        coarse.set_scales(&vec![scale;coarse.cells()]).unwrap();
        let t=AdaptiveTransfer3::new(&coarse,&fine,100_000,||ControlFlow::Continue(())).unwrap();
        let p=t.prepare_two_level(TwoLevelBudget::default(),20_000_000,|_|ControlFlow::Continue(())).unwrap();
        let mut z=vec![0.0;fine.n()];p.apply(&r,&mut z);outputs.push(z);
    }
    assert_eq!(outputs[0],outputs[1]);
}
#[test]
fn prepared_actions_are_positive_symmetric_and_replayable() {
    let tree=Octree3::uniform(1,3,512).unwrap();let coarse=build(&tree);
    let ft=tree.refined(&tree.leaves().iter().copied().collect::<Vec<_>>(),||ControlFlow::Continue(())).unwrap();
    let fine=build(&ft);let t=AdaptiveTransfer3::new(&coarse,&fine,100_000,||ControlFlow::Continue(())).unwrap();
    let p=t.prepare_two_level(TwoLevelBudget::default(),20_000_000,|_|ControlFlow::Continue(())).unwrap();
    let x:Vec<f64>=(0..fine.n()).map(|i|(i%13)as f64-6.0).collect();
    let y:Vec<f64>=(0..fine.n()).map(|i|(i%17)as f64-8.0).collect();
    let mut bx=vec![0.0;fine.n()];let mut by=bx.clone();let mut replay=bx.clone();
    p.apply(&x,&mut bx);p.apply(&y,&mut by);p.apply(&x,&mut replay);assert_eq!(bx,replay);
    let xx:f64=x.iter().zip(&bx).map(|(a,b)|a*b).sum();assert!(xx>0.0);
    let xy:f64=x.iter().zip(&by).map(|(a,b)|a*b).sum();let yx:f64=y.iter().zip(&bx).map(|(a,b)|a*b).sum();
    assert!((xy-yx).abs()<1e-10*xy.abs().max(yx.abs()).max(1.0));
}
#[test]
fn setup_refusals_and_cancellation_do_not_mutate_material() {
    let tree=mixed_tree();let coarse=build(&tree);
    let ft=tree.refined(&tree.leaves().iter().copied().collect::<Vec<_>>(),||ControlFlow::Continue(())).unwrap();
    let fine=build(&ft);let scales=fine.scales().to_vec();
    assert!(fine.prepare_jacobi(0,||ControlFlow::Continue(())).is_err());
    assert!(fine.prepare_jacobi(20_000_000,||ControlFlow::Break(())).is_err());
    let t=AdaptiveTransfer3::new(&coarse,&fine,100_000,||ControlFlow::Continue(())).unwrap();
    assert!(matches!(t.prepare_two_level(TwoLevelBudget{max_coarse_dofs:1,..Default::default()},20_000_000,
        |_|ControlFlow::Continue(())),Err(AdaptivePreconditionError3::Coarse(TwoLevelError::Budget(_)))));
    assert!(matches!(t.prepare_two_level(TwoLevelBudget::default(),20_000_000,
        |w|if w.operator_applications>1{ControlFlow::Break(())}else{ControlFlow::Continue(())}),
        Err(AdaptivePreconditionError3::Coarse(TwoLevelError::Cancelled))));
    assert_eq!(fine.scales(),scales);
}
