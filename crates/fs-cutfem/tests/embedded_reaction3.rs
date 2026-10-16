//! Real embedded supports: signs, balance, state/direct/adjoint derivatives.
use std::ops::ControlFlow;
use fs_cutfem::{CutSdf3,HeightAxis,HexCell};
use fs_cutfem::elastic3::{CutElasticity3,ElasticityError3,ElasticityOptions3,adaptive::AdaptiveElasticity3};
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3,QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_solver::LinearOp;
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self,p:[f64;3])->f64 {(p[0]-0.17)*(p[0]-0.83)}
    fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval {
        let x=Interval::new(lo[0],hi[0]);(x-Interval::new(0.17,0.17))*(x-Interval::new(0.83,0.83))
    }
    fn derivative_enclose(&self,lo:[f64;3],hi:[f64;3],axis:HeightAxis)->Interval {
        if axis==HeightAxis::X {Interval::new(2.0,2.0)*Interval::new(lo[0],hi[0])-Interval::new(1.0,1.0)}
        else {Interval::new(0.0,0.0)}
    }
}
fn build(mixed:bool,nu:f64)->AdaptiveElasticity3 {
    let tree=Octree3::uniform(1,4,4096).unwrap();
    let tree=if mixed {tree.refined(&[*tree.leaves().iter().next().unwrap()],||ControlFlow::Continue(())).unwrap()}else{tree};
    let mut p=|_|ControlFlow::Continue(());let mut q=QuadratureControl3::new(QuadratureOptions3::default(),&mut p).unwrap();
    AdaptiveElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0;3],[1.0;3]).unwrap(),&tree,&Slab,
        &IsotropicElastic::new(1.0,nu,1.0).unwrap(),&|_|false,&|_,_|true,
        ElasticityOptions3::default(),Default::default(),Default::default(),&mut q).unwrap()
}
fn right(_:[f64;3],n:[f64;3])->[f64;3] {if n[0]>0.0 {[1.0,0.0,0.0]}else{[0.0;3]}}
fn motion(p:[f64;3],n:[f64;3])->[f64;3] {if n[0]>0.0 {[0.02,0.003*p[1],0.0]}else{[0.0,0.0,0.002*p[2]]}}
fn solve(op:&AdaptiveElasticity3,b:&[f64])->Vec<f64> {
    op.solve_controlled(b,1e-11,20000,32,|_|ControlFlow::Continue(())).unwrap().coefficients().to_vec()
}
fn dot(a:&[f64],b:&[f64])->f64 {a.iter().zip(b).map(|(a,b)|a*b).sum()}
#[test]
fn g1_affine_extension_has_the_correct_cartesian_and_mixed_grid_reaction_sign() {
    let g=|p:[f64;3],_:[f64;3]|[0.02*(p[0]-0.17),0.0,0.0];
    for mixed in [false,true] {
        let mut op=build(mixed,0.0);op.set_scales(&vec![0.6;op.cells()]).unwrap();
        let b=op.prescribed_displacement_load(&g,||ControlFlow::Continue(())).unwrap();let u=solve(&op,&b);
        let r=op.embedded_reaction(&u,Some(&g),&right,||ControlFlow::Continue(())).unwrap();
        assert!((r.value-0.012).abs()<1e-8,"{}",r.value);
        let left=|_:[f64;3],n:[f64;3]|if n[0]<0.0{[1.0,0.0,0.0]}else{[0.0;3]};
        assert!((op.embedded_reaction(&u,Some(&g),&left,||ControlFlow::Continue(())).unwrap().value+0.012).abs()<1e-8);
    }
    let mut poll=|_|ControlFlow::Continue(());let mut q=QuadratureControl3::new(QuadratureOptions3::default(),&mut poll).unwrap();
    let op=CutElasticity3::build_with_embedded_dirichlet(HexCell::try_new([0.0;3],[1.0;3]).unwrap(),[2;3],&Slab,
        &IsotropicElastic::new(1.0,0.0,1.0).unwrap(),&|_|false,&|_,_|true,ElasticityOptions3::default(),Default::default(),Default::default(),&mut q).unwrap();
    let b=op.prescribed_displacement_load(&g,||ControlFlow::Continue(())).unwrap();
    let u=op.solve_controlled(&b,1e-11,20000,32,|_|ControlFlow::Continue(())).unwrap();
    assert!((op.embedded_reaction(u.coefficients(),Some(&g),&right,||ControlFlow::Continue(())).unwrap().value-0.02).abs()<1e-8);
}
#[test]
fn g0_force_and_moment_balance_hold_for_rigid_virtual_modes() {
    let op=build(true,0.3);let force=[0.01,0.02,-0.03];
    let b=op.body_load(&|_|force,||ControlFlow::Continue(())).unwrap();let u=solve(&op,&b);
    let volume: f64=op.volumes().iter().sum();let resultant=force.map(|v|v*volume);
    let moment=[0.5*(resultant[2]-resultant[1]),0.5*(resultant[0]-resultant[2]),0.5*(resultant[1]-resultant[0])];
    for axis in 0..3 {
        let h=move |_: [f64;3],_: [f64;3]|std::array::from_fn(|i|if i==axis{1.0}else{0.0});
        assert!((op.embedded_reaction(&u,None,&h,||ControlFlow::Continue(())).unwrap().value+resultant[axis]).abs()<1e-8);
        let rotation=move |p:[f64;3],_:[f64;3]|match axis {0=>[0.0,-p[2],p[1]],1=>[p[2],0.0,-p[0]],_=>[-p[1],p[0],0.0]};
        assert!((op.embedded_reaction(&u,None,&rotation,||ControlFlow::Continue(())).unwrap().value+moment[axis]).abs()<1e-8);
    }
}
#[test]
fn g3_state_and_direct_scale_derivatives_match_fixed_field_differences() {
    let mut op=build(true,0.3);let scales=vec![0.55;op.cells()];op.set_scales(&scales).unwrap();
    let u:Vec<f64>=(0..op.n()).map(|i|0.002*(i%11)as f64).collect();
    let reaction=op.embedded_reaction(&u,Some(&motion),&right,||ControlFlow::Continue(())).unwrap();
    let direction:Vec<f64>=(0..op.n()).map(|i|((i*7)%13)as f64/13.0-0.5).collect();let h=1e-5;
    let a:Vec<f64>=u.iter().zip(&direction).map(|(u,d)|u+h*d).collect();let b:Vec<f64>=u.iter().zip(&direction).map(|(u,d)|u-h*d).collect();
    let fd=(op.embedded_reaction(&a,Some(&motion),&right,||ControlFlow::Continue(())).unwrap().value
        -op.embedded_reaction(&b,Some(&motion),&right,||ControlFlow::Continue(())).unwrap().value)/(2.0*h);
    assert!((fd-dot(&reaction.displacement_gradient,&direction)).abs()<1e-8);
    for i in 0..op.cells() {
        let mut s=scales.clone();s[i]+=h;op.set_scales(&s).unwrap();let a=op.embedded_reaction(&u,Some(&motion),&right,||ControlFlow::Continue(())).unwrap().value;
        s[i]-=2.0*h;op.set_scales(&s).unwrap();let b=op.embedded_reaction(&u,Some(&motion),&right,||ControlFlow::Continue(())).unwrap().value;
        assert!(((a-b)/(2.0*h)-reaction.scale_gradient[i]).abs()<1e-8);
    }
}
#[test]
fn g3_total_adjoint_derivative_includes_direct_reaction_and_motion_lifting() {
    let mut op=build(false,0.3);let scales:Vec<f64>=(0..op.cells()).map(|i|0.3+0.04*i as f64).collect();op.set_scales(&scales).unwrap();
    let f=op.body_load(&|_|[0.001,0.0,-0.002],||ControlFlow::Continue(())).unwrap();
    let forward=|op:&AdaptiveElasticity3| {
        let mut b=op.prescribed_displacement_load(&motion,||ControlFlow::Continue(())).unwrap();for (b,f) in b.iter_mut().zip(&f){*b+=f;}
        solve(op,&b)
    };
    let u=forward(&op);let r=op.embedded_reaction(&u,Some(&motion),&right,||ControlFlow::Continue(())).unwrap();
    let z=solve(&op,&r.displacement_gradient);let k=op.scale_bilinear_forms(&z,&u,||ControlFlow::Continue(())).unwrap();
    let b=op.prescribed_displacement_scale_work(&motion,&z,||ControlFlow::Continue(())).unwrap();
    for i in 0..op.cells() {
        let h=1e-4;let mut s=scales.clone();s[i]+=h;op.set_scales(&s).unwrap();let plus=op.embedded_reaction(&forward(&op),Some(&motion),&right,||ControlFlow::Continue(())).unwrap().value;
        s[i]-=2.0*h;op.set_scales(&s).unwrap();let minus=op.embedded_reaction(&forward(&op),Some(&motion),&right,||ControlFlow::Continue(())).unwrap().value;
        let fd=(plus-minus)/(2.0*h);assert!((fd-(r.scale_gradient[i]+b[i]-k[i])).abs()<1e-7);
    }
}
#[test]
fn g4_invalid_or_cancelled_reactions_never_publish_partial_values_or_mutate_physics() {
    let op=build(false,0.3);let u=vec![0.0;op.n()];let scales=op.scales().to_vec();
    assert!(op.embedded_reaction(&[0.0],None,&right,||ControlFlow::Continue(())).is_err());
    assert!(op.embedded_reaction(&u,None,&|_,_|[f64::NAN;3],||ControlFlow::Continue(())).is_err());
    let mut calls=0;let result=op.embedded_reaction(&u,Some(&motion),&right,||{calls+=1;if calls>12{ControlFlow::Break(())}else{ControlFlow::Continue(())}});
    assert!(matches!(result,Err(ElasticityError3::Cancelled)));assert_eq!(op.scales(),scales);
    let a=op.embedded_reaction(&u,Some(&motion),&right,||ControlFlow::Continue(())).unwrap();let b=op.embedded_reaction(&u,Some(&motion),&right,||ControlFlow::Continue(())).unwrap();
    assert_eq!(a.value,b.value);assert_eq!(a.displacement_gradient,b.displacement_gradient);assert_eq!(a.scale_gradient,b.scale_gradient);
}
