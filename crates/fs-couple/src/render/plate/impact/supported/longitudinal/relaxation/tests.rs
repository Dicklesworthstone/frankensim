use super::*;
use crate::render::plate::impact::{ImpactConfig,ImpactSystem,VolumeSpring};
use crate::modal_acoustic_time::ModalAcousticState;
use fs_exec::CancelGate;

fn support(ends:StringEndMotion)->SupportedStrings {
    let wire=WireSpan {endpoints_m:[[0.0,0.0],[0.2,0.0]],linear_density_kg_m:0.012,
        tension_n:8.0,bending_stiffness_n_m2:1e-7,damping_per_s:vec![0.0]};
    let geometry=StringSupportGeometry::try_new([[0.0;3],[0.2,0.0,0.0]],
        [0.0,1.0,0.0],[0.6,0.8,0.0]).unwrap();
    SupportedStrings::with_longitudinal_motion(&[wire],
        &[Some(StringStretching {axial_rigidity_n:100.0,maximum_slope:0.2})],&[ends],&[geometry],
        &[LongitudinalSupportModes {solid_linear_density_kg_m:0.01,modes:vec![1,2],maximum_strain:0.02}],
        TranslatingSupport {mass_kg:0.01,stiffness_n_m:12.0,damping_n_s_m:0.0,
            initial_position_m:2e-5,initial_velocity_m_s:0.0,maximum_travel_m:0.01,maximum_slope:0.2}).unwrap()
}
fn arms()->Vec<SupportedRelaxationArm> {
    vec![SupportedRelaxationArm {string:0,axial_rigidity_n:50.0,relaxation_time_s:0.001},
        SupportedRelaxationArm {string:0,axial_rigidity_n:20.0,relaxation_time_s:0.003}]
}
fn config()->ImpactConfig {ImpactConfig {dt_s:1e-6,max_steps:128,maximum_energy_j:1.0,
    energy_absolute_tolerance_j:1e-12,energy_relative_tolerance:1e-8,maximum_generalized_force:1e3}}
fn initial()->Vec<ModalAcousticState> {(0..4).map(|i|ModalAcousticState {
    displacement_m_sqrt_kg:if i<2 {1e-6*(i as f64+1.0)} else {1e-8*(i as f64-1.2)},
    velocity_m_sqrt_kg_per_s:0.0,
}).collect()}
fn near(a:f64,b:f64,relative:f64) {
    assert!((a-b).abs()<1e-11+relative*a.abs().max(b.abs()),"{a} != {b}");
}

// Independent spatial quadrature: reconstruct absolute fields and their
// canonical derivatives directly from sine functions and the kinetic cross
// integral. No production strain-polynomial or support pullback is called.
fn spatial(s:&SupportedStrings,state:&[f64])->(f64,Vec<f64>,f64) {
    let ends=s.0.end_motion[0];let sign=match ends {StringEndMotion::Both=>0.0,
        StringEndMotion::First=>-1.0,StringEndMotion::Second=>1.0};
    let length=0.2;let weight=s.force_weight();let x=weight*state[8];
    let mass=[0.012,0.012,0.01,0.01];let directions=[0.8,0.0,0.6,0.6];
    let indices=[1usize,1,1,2];let mut c=[0.0;4];let mut relative=[0.0;4];
    for i in 0..4 {
        let k=indices[i] as f64*core::f64::consts::PI;
        let integral=match ends {StringEndMotion::First=>1.0/k,
            StringEndMotion::Second=>-k.cos()/k,StringEndMotion::Both=>(1.0-k.cos())/k};
        c[i]=mass[i]*length*(2.0_f64/(mass[i]*length)).sqrt()*directions[i]*integral;
        relative[i]=state[2*i]-c[i]*x;
    }
    let mut energy=0.0;let mut power=0.0;let mut gradient=vec![0.0;state.len()];
    let cells=4096;let dx=length/cells as f64;
    for cell in 0..cells {
        let station=(cell as f64+0.5)/cells as f64;let mut shapes=[0.0;4];
        for i in 0..4 {
            let phase=indices[i] as f64*core::f64::consts::PI;
            shapes[i]=(2.0/(mass[i]*length)).sqrt()*phase/length*(phase*station).cos();
        }
        let y=0.8*sign*x/length+shapes[0]*relative[0];let z=shapes[1]*relative[1];
        let epsilon=0.6*sign*x/length+shapes[2]*relative[2]+shapes[3]*relative[3]+0.5*(y*y+z*z);
        let b=[y*shapes[0],z*shapes[1],shapes[2],shapes[3]];
        let endpoint=weight*(sign/length*(0.6+0.8*y)-c.iter().zip(b).map(|(c,b)|c*b).sum::<f64>());
        for (arm_index,arm) in arms().iter().enumerate() {
            let mut viscous=0.0;let mut basis=[1.0;3];
            for channel in 0..3 {
                if channel>0 {basis[channel]=(channel as f64*core::f64::consts::PI*station).cos();}
                let coefficient=arm.axial_rigidity_n*length*if channel==0 {1.0} else {0.5};
                viscous+=state[10+3*arm_index+channel]/coefficient.sqrt()*basis[channel];
            }
            let residual=epsilon-viscous;let stress=arm.axial_rigidity_n*residual;
            energy+=0.5*stress*residual*dx;power+=stress*residual/arm.relaxation_time_s*dx;
            for i in 0..4 {gradient[2*i]+=stress*b[i]*dx;}
            gradient[8]+=stress*endpoint*dx;
            for channel in 0..3 {
                let coefficient=arm.axial_rigidity_n*length*if channel==0 {1.0} else {0.5};
                gradient[10+3*arm_index+channel]-=stress*basis[channel]/coefficient.sqrt()*dx;
            }
        }
    }
    (energy,gradient,power)
}

#[test]
fn g1_hereditary_energy_force_tangent_and_oblique_endpoint_match_spatial_work() {
    for ends in [StringEndMotion::First,StringEndMotion::Second,StringEndMotion::Both] {
        let base=support(ends);let loaded=base.clone().with_longitudinal_relaxation(arms(),InitialMemory::Unrelaxed,6).unwrap();
        let plain=ImpactSystem::new(vec![base.body_with_relative_states(&initial()).unwrap()],vec![],vec![],vec![],config()).unwrap();
        let mut run=ImpactSystem::new(vec![loaded.body_with_relative_states(&initial()).unwrap()],vec![],vec![],vec![],config()).unwrap();
        assert_eq!(run.relaxation_branch_count(),6);assert_eq!(run.state().len(),16);
        for i in 0..6 {run.x[10+i]=1e-5*(i as f64-2.1);}
        run.relaxation.as_ref().unwrap().admit_state(run.state()).unwrap();
        let (energy,gradient,power)=spatial(&loaded,run.state());
        let receipt=run.relaxation_observation().unwrap();near(receipt.stored_energy_j,energy,1e-11);
        near(receipt.dissipated_power_w,power,1e-11);
        near(run.stored_energy_j()-plain.stored_energy_j(),energy,1e-10);
        let actual=run.system.effort(run.state());let original=plain.system.effort(plain.state());
        for i in 0..16 {near(actual[i]-original.get(i).copied().unwrap_or(0.0),gradient[i],2e-10);}
        let direction:Vec<_>=(0..16).map(|i|if i<10 {1e-4*(i as f64-3.2)} else {0.03*(i as f64-12.1)}).collect();
        let mut tangent=vec![0.0;16];let mut old=vec![0.0;10];
        assert!(run.hessian_vector(run.state(),&direction,&mut tangent));
        assert!(plain.hessian_vector(plain.state(),&direction[..10],&mut old));
        let h=1e-7;let plus:Vec<_>=run.state().iter().zip(&direction).map(|(x,d)|x+h*d).collect();
        let minus:Vec<_>=run.state().iter().zip(&direction).map(|(x,d)|x-h*d).collect();
        let gp=spatial(&loaded,&plus).1;let gm=spatial(&loaded,&minus).1;
        for i in 0..16 {near(tangent[i]-old.get(i).copied().unwrap_or(0.0),(gp[i]-gm[i])/(2.0*h),3e-6);}
    }
}

#[test]
fn g3_material_history_changes_shared_motion_with_work_and_atomic_refusal() {
    let gate=CancelGate::new_clock_free();let base=support(StringEndMotion::Second);
    let loaded=base.clone().with_longitudinal_relaxation(arms(),InitialMemory::Relaxed,6).unwrap();
    let build=|s:&SupportedStrings|ImpactSystem::new(vec![s.body_with_relative_states(&initial()).unwrap()],vec![],vec![],vec![],config()).unwrap();
    let mut plain=build(&base).prepare_analytic().unwrap();
    let mut reference=build(&loaded);let mut run=build(&loaded).prepare_analytic().unwrap();
    let before=run.stored_energy_j();let(mut supplied,mut loss)=(0.0,0.0);let mut force=vec![0.0;5];
    for step in 0..64 {
        force[4]=if step<12 {0.02*loaded.force_weight()} else {0.0};
        let previous=run.state().to_vec();let frame=run.step(&force,&gate).unwrap();
        reference.step(&force,&gate).unwrap();plain.step(&force,&gate).unwrap();
        supplied+=frame.supplied_work_j;loss+=frame.dissipated_energy_j;
        near(frame.supplied_work_j,force[4]*(run.state()[8]-previous[8]),1e-9);
        assert!((run.stored_energy_j()+loss-before-supplied).abs()<1e-10);
        for (a,b) in run.state().iter().zip(reference.state()) {near(*a,*b,2e-6);}
    }
    assert!(loss>0.0);assert!((run.state()[8]-plain.state()[8]).abs()>1e-13);
    let before=run.state().to_vec();let sample=run.samples();force[0]=2e3;
    assert!(run.step(&force,&gate).is_err());assert_eq!(run.state(),before);assert_eq!(run.samples(),sample);
    let mut invalid=reference.state().to_vec();invalid[10]=(50.0_f64*0.2).sqrt()*0.021;
    assert!(reference.relaxation.as_ref().unwrap().admit_state(&invalid).is_err());
    assert_eq!(reference.samples(),64);
}

#[test]
fn g0_body_composition_retains_source_domain_clock_budget_and_elastic_limit() {
    let base=support(StringEndMotion::Second);
    assert!(base.clone().with_longitudinal_relaxation(arms(),InitialMemory::Relaxed,5).is_err());
    let loaded=base.clone().with_longitudinal_relaxation(arms(),InitialMemory::Relaxed,6).unwrap();
    let upper=loaded.clock_bounds().iter().copied().fold(0.0_f64,f64::max)*1.01;
    let loaded=loaded.bind_longitudinal_relaxation_domain(upper,config().dt_s).unwrap();
    let body=loaded.body_with_relative_states(&initial()).unwrap();
    let mut wrong=config();wrong.dt_s*=0.5;
    assert!(ImpactSystem::new(vec![body.clone()],vec![],vec![],vec![],wrong).is_err());
    // A later stiff volume must not bypass the original source spectrum.
    let volume=VolumeSpring {bulk_modulus_pa:4.0*upper*upper,volume_m3:1.0,areas:vec![0.0,0.0,0.0,0.0,1.0]};
    assert!(ImpactSystem::new(vec![body],vec![],vec![],vec![volume],config()).is_err());
    let empty=base.clone().with_longitudinal_relaxation(vec![],InitialMemory::Relaxed,1).unwrap();
    assert_eq!(base.clock_bounds(),empty.clock_bounds());
    let a=ImpactSystem::new(vec![base.body_with_relative_states(&initial()).unwrap()],vec![],vec![],vec![],config()).unwrap();
    let b=ImpactSystem::new(vec![empty.body_with_relative_states(&initial()).unwrap()],vec![],vec![],vec![],config()).unwrap();
    assert_eq!(a.state(),b.state());assert_eq!(a.stored_energy_j().to_bits(),b.stored_energy_j().to_bits());
    assert_eq!(b.relaxation_branch_count(),0);
    // An explicit elastic source selection still carries frequency authority;
    // zero dynamic histories cannot erase it before a later cold attachment.
    let upper=empty.clock_bounds().iter().copied().fold(0.0_f64,f64::max)*1.01;
    let bounded=empty.bind_longitudinal_relaxation_domain(upper,config().dt_s).unwrap();
    let run=ImpactSystem::new(vec![bounded.body_with_relative_states(&initial()).unwrap()],vec![],vec![],vec![],config()).unwrap();
    assert_eq!(run.state(),a.state());assert_eq!(run.relaxation_branch_count(),0);
    let outside=crate::render::plate::impact::radiation::Model {ports:1,poles:vec![
        crate::render::plate::impact::radiation::Pole {omega:2.0*upper,zeta:0.2,coupling:vec![1.0]},
    ]};
    assert!(run.with_radiation_load(&outside,&[4],1).is_err());
}

#[test]
fn g3_original_radiation_stays_reciprocal_and_rechecks_material_domain() {
    use crate::render::plate::impact::radiation::{Model,Pole};
    use crate::render::plate::impact::squeeze::{ResistiveFilm,FilmCell,FilmChannel,FilmLimits,GapPort};
    use crate::render::plate::impact::gas_film::AmbientGas;
    let s=support(StringEndMotion::Second).with_longitudinal_relaxation(arms(),InitialMemory::Relaxed,6).unwrap();
    let bound=s.clock_bounds().iter().copied().fold(0.0_f64,f64::max);
    let s=s.bind_longitudinal_relaxation_domain(bound*1.01,config().dt_s).unwrap();
    let build=||ImpactSystem::new(vec![s.body_with_relative_states(&initial()).unwrap()],vec![],vec![],vec![],config()).unwrap();
    let model=Model {ports:1,poles:vec![Pole {omega:400.0,zeta:0.2,coupling:vec![1.0]}]};
    let before=build().state().to_vec();
    let mut run=build().with_radiation_load(&model,&[4],1).unwrap().prepare_analytic().unwrap();
    assert_eq!(&run.state()[..before.len()],before.as_slice());assert_eq!(run.state().len(),18);
    let energy=run.stored_energy_j();let mut loss=0.0;
    for _ in 0..32 {loss+=run.step(&[0.0;5],&CancelGate::new_clock_free()).unwrap().dissipated_energy_j;}
    assert!(run.state()[16..].iter().any(|v|*v!=0.0));
    assert!((run.stored_energy_j()+loss-energy).abs()<1e-10);
    let outside=Model {ports:1,poles:vec![Pole {omega:2.0*bound,zeta:0.2,coupling:vec![1.0]}]};
    assert!(build().with_radiation_load(&outside,&[4],1).is_err());
    let zero=Model {ports:1,poles:vec![Pole {omega:2.0*bound,zeta:0.2,coupling:vec![0.0]}]};
    assert_eq!(build().with_radiation_load(&zero,&[4],1).unwrap().state(),before);
    let gap=GapPort {reference_m:0.001,closure:vec![0.0,0.0,0.0,0.0,s.force_weight()]};
    let film=ResistiveFilm::new(vec![FilmCell {area_m2:0.002,gap:gap.clone()}],
        vec![FilmChannel {from:0,to:None,width_m:0.05,length_m:0.01,gap}],5,1.8e-5,
        FilmLimits {minimum_cell_gap_m:1e-6,maximum_gap_m:0.01,maximum_pressure_pa:1000.0}).unwrap();
    assert!(build().with_squeeze_film(film.clone()).is_err());
    assert!(build().with_compressible_squeeze_film(film,AmbientGas {pressure_pa:101325.0,
        temperature_k:293.15,specific_gas_constant_j_kg_k:287.0}).is_err());
}
