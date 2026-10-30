//! Source-matrix derivatives, not finite-difference production sensitivities.
use fs_la::pole_response::PoleResponseOptions;
use fs_math::c64::C64;
use fs_phs::circuit::{Branch,CircuitGraph,assemble_circuit};
use fs_phs::harmonic::{HarmonicBudget,HarmonicModel,HarmonicPole,PoleOptions};
use fs_phs::harmonic::response::{HarmonicParameterDirection as Direction,HarmonicPoleResponseOptions};
use fs_phs::{PortHamiltonian,QuadraticStorage,RelaxationBranch};
use std::cell::Cell;
fn budget()->HarmonicBudget { HarmonicBudget {max_dimension:32,max_dense_bytes:2_000_000,
    max_dense_work:10_000_000,relative_residual_tolerance:1e-10} }
fn options(scale:f64)->HarmonicPoleResponseOptions { HarmonicPoleResponseOptions {
    numerical:PoleResponseOptions {spectral_scale:scale,pole_residual_tolerance:1e-12,
        pole_step_tolerance:1e-12,linear_residual_tolerance:1e-11,max_dimension:32,
        max_vectors:8,max_dense_bytes:2_000_000,max_work:10_000_000},
    max_physical_bytes:2_000_000,max_physical_work:10_000_000} }
fn pole(model:&HarmonicModel,s:C64,v:&[C64],scale:f64)->HarmonicPole {
    model.resonance_pole(s,v,&PoleOptions {spectral_scale:scale,step_tolerance:1e-13,
        residual_tolerance:1e-13,max_iterations:40,max_dimension:32,max_dense_bytes:2_000_000,
        max_work:10_000_000},||false).unwrap()
}
fn close(a:C64,b:C64,t:f64) {assert!((a-b).abs()<=t*(1.0+b.abs()),"{a:?} vs {b:?}");}
fn mechanical()->(HarmonicModel,HarmonicPole) {
    let source=fs_phs::mass_spring_damper(2.5,1600.,0.6).unwrap();
    let model=HarmonicModel::from_ode(&source,budget(),||false).unwrap();
    let s=C64::new(-0.12,-(640.0_f64-0.12*0.12).sqrt());
    let p=pole(&model,s,&[C64::ONE,s.scale(2.5)],26.);(model,p)
}
#[test]
fn g1_mass_stiffness_loss_and_port_directions_match_exact_mechanics() {
    let (model,p)=mechanical();
    let prepared=model.prepare_pole_response(&p,options(26.),||false).unwrap();
    let dm=[0.,0.,0.,-1.0/(2.5*2.5)];let dk=[1.,0.,0.,0.];let dd=[0.,0.,0.,1.];let dg=[0.,1.];
    let directions=[Direction {storage:Some(&dm),..Default::default()},
        Direction {storage:Some(&dk),..Default::default()},
        Direction {dissipation:Some(&dd),..Default::default()},
        Direction {ports:Some(&dg),..Default::default()}];
    let values=prepared.parameter_derivatives(&directions,||false).unwrap();
    let s=p.laplace_pole_per_s;let denominator=s.scale(5.)+C64::from_re(0.6);
    let expected=[-(s*s)/denominator,-denominator.recip(),-s/denominator,C64::ZERO];
    for (d,want) in values.iter().zip(expected) {
        close(d.laplace_pole_per_s,want,1e-10);
        close(d.angular_frequency_rad_s,C64::new(-want.im,want.re),1e-10);
        assert!(d.relative_residual<1e-10);
        assert!(d.decay_power_residual_w.abs()<1e-9);
        assert!(d.state[prepared.normalization_pivot()].abs()<1e-12);
    }
    assert_eq!(prepared.normalization_pivot(),1);
    close(values[0].output[0],C64::from_re(-1.0/(2.5*2.5)),1e-11);
    close(values[3].output[0],prepared.mode_output()[0],1e-12);
    assert!(values[3].state.iter().all(|v|v.abs()==0.));
}
#[test]
fn g1_original_maxwell_strength_and_time_chain_rules_keep_relaxation_poles() {
    let kb=50_000.0_f64;let tau=1.0/300.0;
    let source=fs_phs::mass_spring_damper(1.,60_000.,0.).unwrap()
        .with_relaxation_branches(vec![RelaxationBranch {projection:vec![1.,0.],
            stiffness:kb,relaxation_time_s:tau}]).unwrap();
    let model=HarmonicModel::from_ode(&source,budget(),||false).unwrap();
    let dq=[1.,0.,-0.5/kb.sqrt(), 0.,0.,0., -0.5/kb.sqrt(),0.,0.];
    let mut dr=[0.;9];dr[8]=-1.0/(tau*tau);
    for s in [C64::new(-50.,-87_500.0_f64.sqrt()),C64::from_re(-200.)] {
        let z=C64::from_re(kb.sqrt())/(C64::ONE+s.scale(tau));
        let p=pole(&model,s,&[C64::ONE,s,z],400.);
        let prepared=model.prepare_pole_response(&p,options(400.),||false).unwrap();
        let result=prepared.parameter_derivatives(&[
            Direction {storage:Some(&dq),..Default::default()},
            Direction {dissipation:Some(&dr),..Default::default()},
        ],||false).unwrap();
        let d=C64::ONE+s.scale(tau);
        let slope=s.scale(2.)+C64::from_re(kb*tau)/(d*d);
        close(result[0].laplace_pole_per_s,-s.scale(tau)/d/slope,1e-10);
        close(result[1].laplace_pole_per_s,-s.scale(kb)/(d*d)/slope,1e-9);
        for value in result {
            assert_eq!(value.state.len(),3);
            assert!(value.relative_residual<1e-10);
            assert!(value.decay_power_residual_w.abs()<1e-6*(1.0+value.dissipated_power_w.abs()));
        }
    }
}
#[test]
fn g1_inductive_storage_derivative_preserves_kirchhoff_algebraic_rows() {
    let (r,l,c)=(0.2_f64,0.1_f64,2.2e-7_f64);
    let source=assemble_circuit(&CircuitGraph {node_count:4,branches:vec![
        (1,0,Branch::VoltageSource{port:0}),(1,2,Branch::Inductor{henries:l}),
        (2,3,Branch::Capacitor{farads:c}),(3,0,Branch::Resistor{ohms:r}),
    ],transformers:vec![]}).unwrap();
    let model=HarmonicModel::from_descriptor(&source.system,budget(),||false).unwrap();
    let omega=1.0/(l*c).sqrt();let s=C64::new(-1.,-(omega*omega-1.).sqrt());
    let seed=model.factor(omega,||false).unwrap().response(&[C64::ONE],||false).unwrap().state;
    let p=pole(&model,s,&seed,omega);
    let prepared=model.prepare_pole_response(&p,options(omega),||false).unwrap();
    let n=seed.len();let mut dq=vec![0.;n*n];let flux=source.flux_index[0];
    dq[flux*n+flux]=-1.0/(l*l);
    let d=prepared.parameter_derivatives(&[Direction{storage:Some(&dq),..Default::default()}],||false).unwrap().remove(0);
    close(d.laplace_pole_per_s,-(s*s)/(s.scale(2.*l)+C64::from_re(r)),1e-9);
    assert!(d.state[source.node_potential_index[0]].abs()<1e-10);
    let bad=source.node_potential_index[2];dq[bad*n+bad]=1.;
    assert!(prepared.parameter_derivatives(&[Direction{storage:Some(&dq),..Default::default()}],||false).is_err());
}
fn coupled(j:Vec<f64>,r:Vec<f64>,q:Vec<f64>,g:Vec<f64>)->HarmonicModel {
    let source=PortHamiltonian::new(4,1,j,r,g,Box::new(QuadraticStorage::new(q,4).unwrap())).unwrap();
    HarmonicModel::from_ode(&source,budget(),||false).unwrap()
}
#[test]
fn g3_full_cross_storage_loss_interconnection_and_observation_match_rebuilt_models() {
    let j=vec![0.,0.,1.,0., 0.,0.,0.,1., -1.,0.,0.,0.2, 0.,-1.,-0.2,0.];
    let mut r=vec![0.;16];r[10]=0.4;r[15]=0.6;r[11]=0.05;r[14]=0.05;
    let mut q=vec![0.;16];q[0]=100.;q[5]=225.;q[10]=1.;q[15]=1.;q[1]=-15.;q[4]=-15.;
    let g=vec![0.,0.,1.,0.4];
    let model=coupled(j.clone(),r.clone(),q.clone(),g.clone());
    let seed=model.factor(9.8,||false).unwrap().response(&[C64::ONE],||false).unwrap().state;
    let p=pole(&model,C64::new(-0.2,-9.8),&seed,15.);
    let prepared=model.prepare_pole_response(&p,options(15.),||false).unwrap();
    let mut dj=vec![0.;16];dj[11]=0.3;dj[14]=-0.3;
    let mut dr=vec![0.;16];dr[11]=0.1;dr[14]=0.1;
    let mut dq=vec![0.;16];dq[1]=7.;dq[4]=7.;dq[11]=-0.04;dq[14]=-0.04;
    let dg=vec![0.2,0.,0.3,-0.1];
    let directions=[Direction{interconnection:Some(&dj),..Default::default()},
        Direction{dissipation:Some(&dr),..Default::default()},
        Direction{storage:Some(&dq),..Default::default()},
        Direction{ports:Some(&dg),..Default::default()}];
    let result=prepared.parameter_derivatives(&directions,||false).unwrap();
    for (dir,actual) in directions.iter().zip(result) {
        let h=2.5e-4;
        let perturb=|a:&[f64],d:Option<&[f64]>,t:f64|->Vec<f64> {
            a.iter().enumerate().map(|(i,&x)|x+t*d.map_or(0.,|v|v[i])).collect()
        };
        let resolved=|t:f64| {
            let model=coupled(perturb(&j,dir.interconnection,t),perturb(&r,dir.dissipation,t),
                perturb(&q,dir.storage,t),perturb(&g,dir.ports,t));
            let p=pole(&model,p.laplace_pole_per_s,prepared.right_state(),15.);
            let gauge=p.state[prepared.normalization_pivot()];
            (p.laplace_pole_per_s,p.state.iter().map(|&v|v/gauge).collect::<Vec<_>>(),
                p.effort.iter().map(|&v|v/gauge).collect::<Vec<_>>(),p.port_output[0]/gauge)
        };
        let a=resolved(h);let b=resolved(-h);
        close(actual.laplace_pole_per_s,(a.0-b.0).scale(0.5/h),1e-6);
        for i in 0..4 {
            close(actual.state[i],(a.1[i]-b.1[i]).scale(0.5/h),1e-6);
            close(actual.effort[i],(a.2[i]-b.2[i]).scale(0.5/h),1e-6);
        }
        close(actual.output[0],(a.3-b.3).scale(0.5/h),1e-6);
        assert!(actual.relative_residual<1e-10);
    }
}
#[test]
fn g0_g4_signed_tangents_admit_but_wrong_shapes_symmetries_and_cancelled_batches_refuse() {
    let (model,p)=mechanical();let prepared=model.prepare_pole_response(&p,options(26.),||false).unwrap();
    let signed=[0.,0.,0.,-1.];
    let direction=Direction{dissipation:Some(&signed),..Default::default()};
    let polls=Cell::new(0);
    let before=prepared.parameter_derivatives(&[direction,direction],||{polls.set(polls.get()+1);false}).unwrap();
    let total=polls.get();
    for cut in [1,total/2,total] {
        let count=Cell::new(0);
        assert!(prepared.parameter_derivatives(&[direction,direction],||{count.set(count.get()+1);count.get()>=cut}).is_err());
    }
    let after=prepared.parameter_derivatives(&[direction],||false).unwrap();
    assert_eq!(before[0].state,after[0].state);
    assert_eq!(before[0].output,after[0].output);
    let skew_bad=[1.,0.,0.,0.];let asymmetric=[0.,1.,0.,0.];
    for bad in [Direction{interconnection:Some(&skew_bad),..Default::default()},
        Direction{storage:Some(&asymmetric),..Default::default()},
        Direction{dissipation:Some(&[0.]),..Default::default()},
        Direction{ports:Some(&[f64::NAN,0.]),..Default::default()}] {
        assert!(prepared.parameter_derivatives(&[bad],||false).is_err());
    }
    assert!(prepared.parameter_derivatives(&[],||false).is_err());
}

#[test]
fn g1_a_zero_stiffness_mode_has_an_exact_tangent_without_inverting_storage() {
    let source=PortHamiltonian::new(2,1,vec![0.,1.,-1.,0.],vec![0.,0.,0.,0.6],
        vec![0.,1.],Box::new(QuadraticStorage::new(vec![0.,0.,0.,0.4],2).unwrap())).unwrap();
    let model=HarmonicModel::from_ode(&source,budget(),||false).unwrap();
    let p=pole(&model,C64::ZERO,&[C64::ONE,C64::ZERO],1.);
    let prepared=model.prepare_pole_response(&p,options(1.),||false).unwrap();
    let dq=[1.,0.,0.,0.];
    let result=prepared.parameter_derivatives(&[Direction{storage:Some(&dq),..Default::default()}],||false).unwrap().remove(0);
    // m*s^2+d*s+k=0 at k=s=0 gives ds/dk=-1/d.
    close(result.laplace_pole_per_s,C64::from_re(-1.0/0.6),1e-12);
    close(result.state[1],C64::from_re(-2.5/0.6),1e-12);
    assert!((result.stored_energy_j-0.5).abs()<1e-12);
}
