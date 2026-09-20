use fs_couple::modal_acoustic_time::{ModalAcousticMode,ModalAcousticState,ModalAcousticTimeBudget,ModalAcousticTimeModel};
use fs_couple::render::schedule::force::coupled::{CoupledModalSystem,ModalAttachment,ModalConnection,ModalCouplingConfig,ModalCouplingError};
use fs_couple::render::schedule::force::coupled::contact::{ModalContact,ModalContactConfig};
use fs_couple::render::schedule::force::coupled::contact::multiple::MultiContactConfig;
use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::{EquilibriumLinearization,SensitivityBudget};
use fs_dcontact::{Obstacle,OpeningContactStep};
use fs_exec::CancelGate;
use fs_math::c64::C64;

fn config()->ModalCouplingConfig {
    ModalCouplingConfig {max_modes:16,max_connections:4,max_setup_terms:8192,nyquist_guard_fraction:0.9,
        maximum_total_energy_j:1000.0,maximum_abs_pressure_pa:1e6,maximum_abs_connection_force_n:1e6,
        solve_relative_tolerance:1e-10,energy_absolute_tolerance_j:1e-11,energy_relative_tolerance:1e-9}
}
fn cconfig()->ModalContactConfig {
    ModalContactConfig {max_iterations:128,maximum_force_n:1e6,maximum_penetration_m:1.0,
        force_absolute_tolerance_n:1e-12,force_relative_tolerance:1e-13}
}
fn joint()->MultiContactConfig {MultiContactConfig {max_contacts:4,max_sweeps:128,max_setup_terms:8192}}
fn budget()->SensitivityBudget {SensitivityBudget {max_contacts:4,max_setup_terms:8192,max_query_terms:8192,minimum_contact_margin_m:1e-8}}
fn model(rate:u32,w:f64)->ModalAcousticTimeModel {
    ModalAcousticTimeModel::try_new(rate,vec![ModalAcousticMode {angular_frequency_rad_s:w,damping_ratio:0.05,
        pressure_per_modal_velocity:C64::ZERO}],ModalAcousticTimeBudget::audible_reference()).unwrap()
}
// omega0, omega1, spring k, spring rest, spring left/right shape,
// contact k, gap, weight, contact left/right shape, force0, force1.
fn parameters()->[f64;13] {[80.0,110.0,300.0,0.0002,1.0,0.8,3e6,0.0001,0.7,1.1,0.9,20.0,-5.0]}
fn descriptions(p:&[f64;13])->(Vec<ModalConnection>,Vec<(ModalContact,ModalContactConfig)>) {
    let spring=ModalConnection {left:ModalAttachment {component:0,shapes:vec![p[4]]},right:ModalAttachment {component:1,shapes:vec![p[5]]},
        stiffness_n_m:p[2],damping_n_s_m:5.0,rest_extension_m:p[3]};
    let contact=ModalContact {left:ModalAttachment {component:0,shapes:vec![p[9]]},right:ModalAttachment {component:1,shapes:vec![p[10]]},
        law:Obstacle::new(vec![-1.0],1,1,vec![p[7]],vec![p[8]],p[6],1.5,"sensitivity fixture".into()).unwrap()};
    (vec![spring],vec![(contact,cconfig())])
}
fn equilibrium(p:&[f64;13],rate:u32)->(CoupledModalSystem,Vec<(ModalContact,ModalContactConfig)>) {
    let (springs,contacts)=descriptions(p);
    let mut n=CoupledModalSystem::new(vec![model(rate,p[0]),model(rate,p[1])],springs,config(),&CancelGate::new()).unwrap();
    n.initialize_contact_equilibrium(&p[11..],&contacts,joint(),&CancelGate::new()).unwrap();
    (n,contacts)
}
fn q(n:&CoupledModalSystem)->Vec<f64> {n.components().iter().flat_map(|m|m.states()).map(|s|s.displacement_m_sqrt_kg).collect()}
fn goal(p:&[f64;13])->f64 {let (n,_)=equilibrium(p,48000);let q=q(&n);0.7*q[0]-0.25*q[1]}
fn dot(a:&[f64],b:&[f64])->f64 {a.iter().zip(b).map(|(a,b)|a*b).sum()}
fn near(a:f64,b:f64) {assert!((a-b).abs()<=3e-5*a.abs().max(b.abs())+1e-11,"{a:.16e} != {b:.16e}");}

#[test]
fn tangent_action_and_inverse_match_independent_two_by_two_stiffness() {
    let p=parameters();let (n,c)=equilibrium(&p,48000);let gate=CancelGate::new();
    let linear=EquilibriumLinearization::new(&n,&p[11..],&c,budget(),&gate).unwrap();
    assert_eq!(linear.report().active_contacts,1);
    let q=linear.displacement();let x=p[9]*q[0]-p[10]*q[1];
    let kt=OpeningContactStep::new(&c[0].0.law,-x).unwrap().static_differential(0.0).unwrap().closure_stiffness_n_m;
    let a=p[0]*p[0]+p[2]*p[4]*p[4]+kt*p[9]*p[9];
    let b=-p[2]*p[4]*p[5]-kt*p[9]*p[10];
    let d=p[1]*p[1]+p[2]*p[5]*p[5]+kt*p[10]*p[10];
    let rhs=[0.7,-0.25];let solved=linear.solve(&rhs,&gate).unwrap();let det=a*d-b*b;
    near(solved.values[0],(d*rhs[0]-b*rhs[1])/det);
    near(solved.values[1],(a*rhs[1]-b*rhs[0])/det);
    let x=[0.1,-0.3];let applied=linear.apply(&x,&gate).unwrap();near(applied[0],a*x[0]+b*x[1]);near(applied[1],b*x[0]+d*x[1]);
    let y=[-0.7,0.9];let ay=linear.apply(&y,&gate).unwrap();near(dot(&x,&ay),dot(&y,&applied));
    assert!(solved.relative_residual<config().solve_relative_tolerance);
}

#[test]
fn one_adjoint_matches_re_solved_differences_for_all_parameter_families_and_both_shape_signs() {
    let p=parameters();let (n,c)=equilibrium(&p,48000);let gate=CancelGate::new();
    let linear=EquilibriumLinearization::new(&n,&p[11..],&c,budget(),&gate).unwrap();
    let adj=linear.solve(&[0.7,-0.25],&gate).unwrap();let r=linear.parameter_pullback(&adj.values,&gate).unwrap();
    // TOTAL gradient = minus the residual pullback; right-side attachment
    // parameters additionally negate their contribution to the signed column.
    let gradient=[-r.angular_frequencies[0].unwrap(),-r.angular_frequencies[1].unwrap(),
        -r.springs[0].stiffness,-r.springs[0].rest_extension,-r.springs[0].column[0],r.springs[0].column[1],
        -r.contacts[0].stiffness,-r.contacts[0].gap,-r.contacts[0].weight,-r.contacts[0].column[0],r.contacts[0].column[1],
        -r.external_forces[0],-r.external_forces[1]];
    for j in 0..p.len() {
        let h=(p[j].abs()*1e-4).max(1e-8);let mut lo=p;let mut hi=p;lo[j]-=h;hi[j]+=h;
        near(gradient[j],(goal(&hi)-goal(&lo))/(2.0*h));
    }
    assert_eq!(q(&n),linear.displacement());assert_eq!(n.samples_rendered(),0);
}

#[test]
fn separated_and_zero_stiffness_contacts_keep_their_correct_distinct_derivatives() {
    let gate=CancelGate::new();
    for (stiffness,gap,active) in [(3e6,0.5,false),(0.0,0.0001,true)] {
        let mut p=parameters();p[6]=stiffness;p[7]=gap;
        let (n,c)=equilibrium(&p,48000);let l=EquilibriumLinearization::new(&n,&p[11..],&c,budget(),&gate).unwrap();
        let z=l.solve(&[1.0,0.0],&gate).unwrap();let r=l.parameter_pullback(&z.values,&gate).unwrap();
        assert_eq!(l.report().active_contacts,usize::from(active));
        assert_eq!(r.contacts[0].gap,0.0);
        if active {assert!(r.contacts[0].stiffness.abs()>0.0);} else {assert_eq!(r.contacts[0].stiffness,0.0);}
    }
}

#[test]
fn physically_supported_free_mass_has_a_load_sensitivity_without_an_invented_frequency() {
    let gate=CancelGate::new();
    let make=|force:f64| {
        let mass=ModalAcousticTimeModel::try_free_mass(48000,0.04,0.0,0.0,ModalAcousticTimeBudget::audible_reference()).unwrap();
        let support=ModalConnection {left:ModalAttachment {component:0,shapes:vec![5.0]},right:ModalAttachment {component:0,shapes:vec![0.0]},
            stiffness_n_m:400.0,damping_n_s_m:0.0,rest_extension_m:0.0};
        let contact=ModalContact {left:ModalAttachment {component:0,shapes:vec![5.0]},right:ModalAttachment {component:1,shapes:vec![1.0]},
            law:Obstacle::new(vec![-1.0],1,1,vec![0.0001],vec![1.0],3e6,1.5,"supported fixture".into()).unwrap()};
        let contacts=vec![(contact,cconfig())];
        let mut n=CoupledModalSystem::new(vec![mass,model(48000,100.0)],vec![support],config(),&gate).unwrap();
        n.initialize_contact_equilibrium(&[force*5.0,0.0],&contacts,joint(),&gate).unwrap();(n,contacts)
    };
    let (n,c)=make(2.0);let l=EquilibriumLinearization::new(&n,&[10.0,0.0],&c,budget(),&gate).unwrap();
    let z=l.solve(&[5.0,0.0],&gate).unwrap();let r=l.parameter_pullback(&z.values,&gate).unwrap();
    assert_eq!(r.angular_frequencies[0],None);
    let h=1e-4;let (lo,_)=make(2.0-h);let (hi,_)=make(2.0+h);
    near(z.values[0]*5.0,(q(&hi)[0]-q(&lo)[0])*5.0/(2.0*h));
}

#[test]
fn stale_loads_moving_states_and_activity_switches_refuse_without_mutation() {
    let p=parameters();let (n,c)=equilibrium(&p,48000);let gate=CancelGate::new();let before=q(&n);
    assert!(EquilibriumLinearization::new(&n,&[21.0,-5.0],&c,budget(),&gate).is_err());
    assert!(EquilibriumLinearization::new(&n,&p[11..],&[],budget(),&gate).is_err());
    assert_eq!(q(&n),before);
    let mut p=parameters();p[11]=0.0;p[12]=0.0;p[3]=0.0;p[7]=0.0;
    let (n,c)=equilibrium(&p,48000);
    assert!(EquilibriumLinearization::new(&n,&[0.0,0.0],&c,budget(),&gate).is_err());
    let mut moving=model(48000,100.0);moving.restore_states(&[ModalAcousticState {displacement_m_sqrt_kg:0.0,velocity_m_sqrt_kg_per_s:1.0}]).unwrap();
    let n=CoupledModalSystem::new(vec![moving],vec![],config(),&gate).unwrap();
    assert!(EquilibriumLinearization::new(&n,&[0.0],&[],budget(),&gate).is_err());
}

#[test]
fn derivative_budgets_cancellation_and_bad_queries_preserve_replay() {
    let p=parameters();let (n,c)=equilibrium(&p,48000);let gate=CancelGate::new();
    let exact=SensitivityBudget {max_contacts:1,max_setup_terms:45,max_query_terms:15,minimum_contact_margin_m:1e-8};
    assert!(EquilibriumLinearization::new(&n,&p[11..],&c,SensitivityBudget {max_setup_terms:44,..exact},&gate).is_err());
    let l=EquilibriumLinearization::new(&n,&p[11..],&c,exact,&gate).unwrap();let expected=l.solve(&[1.0,0.0],&gate).unwrap();
    let cancelled=CancelGate::new();cancelled.request();
    assert!(matches!(l.solve(&[1.0,0.0],&cancelled),Err(ModalCouplingError::Cancelled)));
    assert!(l.solve(&[],&gate).is_err());assert!(l.parameter_pullback(&[f64::NAN,0.0],&gate).is_err());
    assert_eq!(l.solve(&[1.0,0.0],&gate).unwrap(),expected);
    let capped=EquilibriumLinearization::new(&n,&p[11..],&c,SensitivityBudget {max_query_terms:14,..exact},&gate).unwrap();
    assert!(capped.solve(&[1.0,0.0],&gate).is_err());
}

#[test]
fn static_tangent_does_not_depend_on_the_audio_clock() {
    let p=parameters();let gate=CancelGate::new();let mut previous=None;
    for rate in [24000,48000,192000] {
        let (n,c)=equilibrium(&p,rate);let l=EquilibriumLinearization::new(&n,&p[11..],&c,budget(),&gate).unwrap();
        let z=l.solve(&[0.7,-0.25],&gate).unwrap();let r=l.parameter_pullback(&z.values,&gate).unwrap();
        if let Some(old)=previous {assert_eq!(old,r);}previous=Some(r);
    }
}

use fs_couple::render::schedule::force::coupled::equilibrium::sensitivity::objective::DisplacementTarget;
fn physical_target(p:&[f64;13])->DisplacementTarget {
    DisplacementTarget {attachment:ModalAttachment {component:1,shapes:vec![p[10]]},target_m:0.0004,scale_m:0.001,weight:1.3}
}
fn physical_objective(p:&[f64;13])->f64 {
    let (n,c)=equilibrium(p,48000);let gate=CancelGate::new();
    let l=EquilibriumLinearization::new(&n,&p[11..],&c,budget(),&gate).unwrap();
    l.displacement_objective(&[physical_target(p)],&[],1,&gate).unwrap().value
}

#[test]
fn physical_objective_includes_actuator_units_and_shared_observation_shape_terms() {
    let p=parameters();let (n,c)=equilibrium(&p,48000);let gate=CancelGate::new();
    let l=EquilibriumLinearization::new(&n,&p[11..],&c,budget(),&gate).unwrap();
    let actuator=ModalAttachment {component:0,shapes:vec![2.0]};
    let result=l.displacement_objective(&[physical_target(&p)],&[actuator],2,&gate).unwrap();
    let h=1e-4;let mut lo=p;let mut hi=p;lo[11]-=2.0*h;hi[11]+=2.0*h;
    near(result.physical_force_gradient[0],(physical_objective(&hi)-physical_objective(&lo))/(2.0*h));
    // This parameter changes the contact's right force/closure map AND the
    // sensor's physical observation. Neither derivative can be silently lost.
    lo=p;hi=p;lo[10]-=h;hi[10]+=h;
    let combined=result.residual_pullback.contacts[0].column[1]+result.observation_shape_gradient[0][0];
    near(combined,(physical_objective(&hi)-physical_objective(&lo))/(2.0*h));
    assert_eq!(n.samples_rendered(),0);
}

#[test]
fn weighted_target_scale_and_target_partials_match_explicit_objective_differences() {
    let p=parameters();let (n,c)=equilibrium(&p,48000);let gate=CancelGate::new();
    let l=EquilibriumLinearization::new(&n,&p[11..],&c,budget(),&gate).unwrap();
    let target=physical_target(&p);let r=l.displacement_objective(&[target.clone()],&[],1,&gate).unwrap();
    let h=1e-8;
    let eval=|t:DisplacementTarget|l.displacement_objective(&[t],&[],1,&gate).unwrap().value;
    let mut lo=target.clone();let mut hi=target.clone();lo.target_m-=h;hi.target_m+=h;
    near(r.target_gradient[0],(eval(hi)-eval(lo))/(2.0*h));
    let mut lo=target.clone();let mut hi=target.clone();lo.scale_m-=h;hi.scale_m+=h;
    near(r.scale_gradient[0],(eval(hi)-eval(lo))/(2.0*h));
    let mut lo=target.clone();let mut hi=target.clone();lo.weight-=1e-4;hi.weight+=1e-4;
    near(r.weight_gradient[0],(eval(hi)-eval(lo))/2e-4);
    let repeated=l.displacement_objective(&[target.clone(),target],&[],2,&gate).unwrap();
    near(repeated.value,2.0*r.value);
    for (a,b) in repeated.adjoint.values.iter().zip(r.adjoint.values) {near(*a,2.0*b);}
}

#[test]
fn objective_admission_and_cancellation_never_return_partial_gradient_families() {
    let p=parameters();let (n,c)=equilibrium(&p,48000);let gate=CancelGate::new();
    let l=EquilibriumLinearization::new(&n,&p[11..],&c,budget(),&gate).unwrap();let t=physical_target(&p);
    let expected=l.displacement_objective(&[t.clone()],&[],1,&gate).unwrap();
    let mut bad=t.clone();bad.scale_m=0.0;
    assert!(l.displacement_objective(&[t.clone(),bad],&[],2,&gate).is_err());
    let mut bad=t.clone();bad.attachment.component=99;
    assert!(l.displacement_objective(&[bad],&[],1,&gate).is_err());
    assert!(l.displacement_objective(&[t.clone()],&[],0,&gate).is_err());
    assert!(l.displacement_objective(&[],&[],1,&gate).is_err());
    let cancelled=CancelGate::new();cancelled.request();
    assert!(matches!(l.displacement_objective(&[t.clone()],&[],1,&cancelled),Err(ModalCouplingError::Cancelled)));
    assert_eq!(l.displacement_objective(&[t],&[],1,&gate).unwrap(),expected);
}
