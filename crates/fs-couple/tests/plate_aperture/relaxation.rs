use super::*;
use fs_couple::bernoulli_aperture::dynamic::{DynamicAperture,ApertureDrive,ApertureState};
use fs_couple::bernoulli_aperture::dynamic::relaxation::{PlateRelaxationSpec,PlateRelaxationRegion,InitialApertureMemory};
use fs_couple::bernoulli_aperture::plate::closure::PlateClosureSpec;
use fs_couple::bernoulli_aperture::tube::{ApertureTube,UniformTubeSpec,TubeDrive,TubeFrame};
use fs_material::visco::GeneralizedMaxwell;
use fs_dcontact::Obstacle;
const DT:f64=1e-5;
fn plate(tension:f64)->PlateApertureReduction {
    let (c,mut o)=fixture(4e9,900.0);o.damping_ratio=0.0;o.assembly.pretension=tension;
    PlateApertureReduction::from_chart(c,o,&CancelGate::new()).unwrap()
}
fn spec(p:&PlateApertureReduction)->PlateRelaxationSpec {
    PlateRelaxationSpec {regions:vec![PlateRelaxationRegion {
        triangles:(0..p.chart().mesh.tris.len()).collect(),
        material:GeneralizedMaxwell::new(4e9,vec![(2e9,0.001),(1e9,0.01)]).unwrap(),
        poisson_ratio:0.3,band_hz:(0.0,1000.0),provenance:"synthetic supplied Maxwell material".into(),
    }],max_branches:64,max_dt_over_tau:0.1,max_angular_step:0.1}
}
fn base(p:PlateApertureReduction,z:f64,initial:ApertureState)->DynamicAperture {
    let lay=Obstacle::new(vec![-1.0],1,1,vec![0.0],vec![1.0],0.0,2.0,"explicit no-contact control".into()).unwrap();
    DynamicAperture::from_plate(p,1.2,z,DT,2048,initial,lay).unwrap()
}
fn rest(p:&PlateApertureReduction)->ApertureState {ApertureState {opening_m:p.options().rest_opening_m,opening_velocity_m_s:0.0}}
fn tube_spec()->UniformTubeSpec {
    UniformTubeSpec {length_m:32.0*(343.0*DT),radius_m:0.007,sound_speed_m_s:343.0,
        terminal_reflection:-0.8,max_length_error_m:1e-12,max_wave_memory_bytes:1<<20}
}
fn profiled(p:PlateApertureReduction)->DynamicAperture {
    let s=spec(&p);let initial=rest(&p);
    let closure=PlateClosureSpec {
        nodal_rest_gap_m:p.chart().mesh.nodes.iter().map(|q|0.0002*(0.05+1.9*q.1/0.01)).collect(),
        lay_triangles:(0..p.chart().mesh.tris.len()).collect(),stiffness_pa_per_m_alpha:1e12,
        alpha:2.0,internal_loss_s_per_m:0.5,provenance:"synthetic nonuniform lay".into(),max_penetration_m:0.0002,
    };
    DynamicAperture::from_plate_with_closure(p,closure,1.2,tube_spec().characteristic_impedance(1.2).unwrap(),
        DT,2048,initial).unwrap().with_plate_relaxation(s,InitialApertureMemory::Unrelaxed).unwrap()
}
fn near_scaled(a:f64,b:f64,scale:f64) {assert!((a-b).abs()<=2e-8*scale.max(1e-20),"{a:e} != {b:e}");}

#[test]
fn actual_regional_bending_projection_excludes_membrane_prestress() {
    let p=plate(10.0);let mut s=spec(&p);let half=p.chart().mesh.tris.len()/2;
    let mut right=s.regions[0].clone();s.regions[0].triangles=(0..half).collect();
    right.triangles=(half..p.chart().mesh.tris.len()).collect();right.material.terms=vec![(1e9,0.005)];
    s.regions.push(right);
    let m=base(p.clone(),1e6,rest(&p)).with_plate_relaxation(s,InitialApertureMemory::Relaxed).unwrap();
    let mut options=p.options().assembly.clone();options.pretension=0.0;
    let bending=p.chart().assemble(&[],&options).unwrap();
    let mut q=vec![0.0;bending.free];
    for (n,v) in p.shape_per_opening().iter().enumerate() {for c in 0..3 {
        if let Some(i)=bending.dof_map[3*n+c] {q[i]=v[c];}
    }}
    let mut kq=vec![0.0;q.len()];bending.k.spmv(&q,&mut kq);
    let expected=q.iter().zip(kq).map(|(q,k)|q*k).sum::<f64>();
    let memory=m.relaxation().unwrap();
    near_scaled(memory.region_bending_stiffness_n_m().iter().sum(),expected,expected);
    assert!(expected<p.stiffness_n_m());assert_eq!(memory.branches().len(),3);
    near_scaled(memory.branches()[0].stiffness,0.5*memory.region_bending_stiffness_n_m()[0],expected);
    near_scaled(memory.branches()[2].stiffness,0.25*memory.region_bending_stiffness_n_m()[1],expected);
    assert_eq!(m.spec().stiffness_n_m,p.stiffness_n_m());
    assert_eq!(m.plate_reduction().unwrap().pressure_area_m2(),p.pressure_area_m2());
}

#[test]
fn coupled_junction_matches_full_phs_material_mechanics_without_a_lagged_force() {
    let p=plate(0.0);let s=spec(&p);let h=p.options().rest_opening_m;
    let initial=ApertureState {opening_m:h-3e-6,opening_velocity_m_s:0.0};
    let mut actual=base(p.clone(),1e6,initial).with_plate_relaxation(s,InitialApertureMemory::Unrelaxed).unwrap();
    let arms=actual.relaxation().unwrap().branches().iter().map(|a|fs_phs::RelaxationBranch {
        projection:vec![1.0,0.0],stiffness:a.stiffness,relaxation_time_s:a.relaxation_time_s,
    }).collect();
    let full=fs_phs::mass_spring_damper(p.mass_kg(),p.stiffness_n_m(),0.0).unwrap()
        .with_relaxation_branches(arms).unwrap();
    let mut x=vec![initial.opening_m-h,0.0,0.0,0.0];
    for n in 0..128 {
        let force=if n<64 {-1e-4}else{1e-4};
        let result=fs_phs::step(&full,&x,&[force],DT).unwrap();
        let next=h+result.x[0];let vm=f64::midpoint(x[1],result.x[1])/p.mass_kg();
        let area=p.pressure_area_m2();let dp=-force/area;
        let opening=f64::midpoint(h+x[0],next);
        let jet=p.width_m()*opening*dp.signum()*(2.0*dp.abs()/1.2).sqrt();
        let outgoing=1e6*(jet-area*vm);
        let frame=actual.step(ApertureDrive {upstream_pressure_pa:outgoing+dp,incoming_pressure_pa:0.0,body_flow_m3_s:0.0}).unwrap();
        near_scaled(frame.state.opening_m,next,h);
        near_scaled(frame.state.opening_velocity_m_s,result.x[1]/p.mass_kg(),0.1);
        for (&z,&expected) in actual.relaxation().unwrap().memory_sqrt_j().iter().zip(&result.x[2..]) {
            near_scaled(z,expected,1e-3);
        }
        let scale=frame.stored_energy_j+frame.storage_change_j.abs()+frame.dissipated_energy_j+frame.pressure_work_j.abs();
        assert!(frame.balance_residual_j().abs()<2e-8*scale.max(1e-20));
        assert!(frame.relaxation_loss_j>=0.0);
        x=result.x;
    }
}

#[test]
fn prior_material_history_changes_motion_even_at_identical_initial_opening_and_velocity() {
    let p=plate(0.0);let s=spec(&p);let initial=ApertureState {opening_m:p.options().rest_opening_m-1e-5,opening_velocity_m_s:0.0};
    let mut relaxed=base(p.clone(),1e6,initial).with_plate_relaxation(s.clone(),InitialApertureMemory::Relaxed).unwrap();
    let mut unrelaxed=base(p.clone(),1e6,initial).with_plate_relaxation(s,InitialApertureMemory::Unrelaxed).unwrap();
    assert_eq!(relaxed.state(),unrelaxed.state());assert_eq!(relaxed.relaxation().unwrap().stored_energy_j(),0.0);
    assert!(unrelaxed.relaxation().unwrap().stored_energy_j()>0.0);
    let mut difference=0.0_f64;
    for _ in 0..512 {
        let a=relaxed.step(ApertureDrive::default()).unwrap();let b=unrelaxed.step(ApertureDrive::default()).unwrap();
        difference=difference.max((a.state.opening_m-b.state.opening_m).abs());
    }
    assert!(difference>1e-8);
    assert!(relaxed.relaxation().unwrap().dissipated_energy_j()>0.0);
    assert!(unrelaxed.relaxation().unwrap().dissipated_energy_j()>0.0);
}

#[test]
fn profiled_contact_tube_and_material_memory_resume_with_identical_complete_history() {
    let p=plate(0.0);let mut reference=ApertureTube::new(profiled(p.clone()),tube_spec()).unwrap();
    let mut resumed=ApertureTube::new(profiled(p),tube_spec()).unwrap();
    let drives:Vec<_>=(0..513).map(|n|TubeDrive {upstream_pressure_pa:if n<256 {5.0}else{0.0},body_flow_m3_s:0.0}).collect();
    let expected:Vec<_>=drives.iter().map(|&d|reference.step(d).unwrap()).collect();
    let mut frames=vec![TubeFrame::default();513];let gate=CancelGate::new();
    resumed.advance_block(&drives[..37],&mut frames[..37],&gate).unwrap();
    let before=resumed.aperture().relaxation().unwrap().memory_sqrt_j().to_vec();
    let cancel=CancelGate::new();cancel.request();
    assert_eq!(resumed.advance_block(&drives[37..],&mut frames[37..],&cancel).unwrap().completed,0);
    assert_eq!(resumed.aperture().relaxation().unwrap().memory_sqrt_j(),before);
    assert!(resumed.step(TubeDrive {upstream_pressure_pa:f64::NAN,..TubeDrive::default()}).is_err());
    assert_eq!(resumed.aperture().relaxation().unwrap().memory_sqrt_j(),before);
    for (input,out) in drives[37..].chunks(17).zip(frames[37..].chunks_mut(17)) {resumed.advance_block(input,out,&gate).unwrap();}
    assert_eq!(frames,expected);
    assert_eq!(resumed.aperture().relaxation().unwrap().memory_sqrt_j(),reference.aperture().relaxation().unwrap().memory_sqrt_j());
    for f in frames {let scale=f.stored_energy_j+f.dissipated_energy_j+f.upstream_work_j.abs()+f.storage_change_j.abs();
        assert!(f.balance_residual_j().abs()<2e-8*scale.max(1e-20));}
    assert!(resumed.aperture().relaxation().unwrap().dissipated_energy_j()>0.0);
    assert!(resumed.aperture().plate_closure().is_some());
}

#[test]
fn empty_spectrum_keeps_the_original_physical_trajectory_exactly() {
    let p=plate(0.0);let mut s=spec(&p);s.regions[0].material.terms.clear();
    let mut a=base(p.clone(),1e6,rest(&p));let mut b=base(p.clone(),1e6,rest(&p))
        .with_plate_relaxation(s,InitialApertureMemory::Relaxed).unwrap();
    for n in 0..256 {let drive=ApertureDrive {upstream_pressure_pa:if n<64 {5.0}else{0.0},..ApertureDrive::default()};
        assert_eq!(a.step(drive).unwrap(),b.step(drive).unwrap());}
}

#[test]
fn unmatched_material_missing_regions_unresolved_clocks_and_late_history_are_rejected() {
    let p=plate(0.0);let good=spec(&p);let mut cases=Vec::new();
    let mut s=good.clone();s.regions[0].material.e_inf=5e9;cases.push(s);
    let mut s=good.clone();s.regions[0].triangles.pop();cases.push(s);
    let mut s=good.clone();s.regions[0].triangles.push(0);cases.push(s);
    let mut s=good.clone();s.regions[0].material.terms[0].1=1e-10;cases.push(s);
    let mut s=good.clone();s.regions[0].band_hz=(0.0,1.0);cases.push(s);
    let mut s=good.clone();s.max_angular_step=1e-10;cases.push(s);
    let mut s=good.clone();s.max_branches=1;cases.push(s);
    let mut s=good.clone();s.regions[0].material.terms[0].0=f64::NAN;cases.push(s);
    for s in cases {assert!(base(p.clone(),1e6,rest(&p)).with_plate_relaxation(s,InitialApertureMemory::Relaxed).is_err());}
    assert!(base(p.clone(),1e6,rest(&p)).with_plate_relaxation(good.clone(),InitialApertureMemory::ViscousDisplacementM(vec![])).is_err());
    let mut started=base(p.clone(),1e6,rest(&p));started.step(ApertureDrive::default()).unwrap();
    assert!(started.with_plate_relaxation(good.clone(),InitialApertureMemory::Relaxed).is_err());
    let damped=reduction(4e9,900.0);
    assert!(base(damped.clone(),1e6,rest(&damped)).with_plate_relaxation(good,InitialApertureMemory::Relaxed).is_err());
}

#[test]
fn a_late_slope_refusal_does_not_consume_viscous_or_reflected_wave_history() {
    let (chart,mut options)=fixture(4e9,900.0);options.damping_ratio=0.0;options.max_slope=1e-6;
    let p=PlateApertureReduction::from_chart(chart,options,&CancelGate::new()).unwrap();
    let pressure=p.closing_pressure_pa();let s=spec(&p);
    let build=||ApertureTube::new(base(p.clone(),tube_spec().characteristic_impedance(1.2).unwrap(),rest(&p))
        .with_plate_relaxation(s.clone(),InitialApertureMemory::Unrelaxed).unwrap(),tube_spec()).unwrap();
    let mut a=build();let mut b=build();
    let gentle=TubeDrive {upstream_pressure_pa:pressure*1e-7,body_flow_m3_s:0.0};
    for _ in 0..70 {assert_eq!(a.step(gentle).unwrap(),b.step(gentle).unwrap());}
    let state=b.aperture().state();let memory=b.aperture().relaxation().unwrap().memory_sqrt_j().to_vec();
    let energy=b.stored_energy_j();let loss=b.aperture().relaxation().unwrap().dissipated_energy_j();
    let error=b.step(TubeDrive {upstream_pressure_pa:pressure*100.0,body_flow_m3_s:0.0}).unwrap_err();
    assert!(matches!(error,fs_couple::acoustic_realize::AcousticRealizeError::InvalidDescription {
        what:"plate aperture exceeds its declared linear-slope domain"
    }),"{error:?}");
    assert_eq!(b.aperture().state(),state);assert_eq!(b.stored_energy_j(),energy);
    assert_eq!(b.aperture().relaxation().unwrap().memory_sqrt_j(),memory);
    assert_eq!(b.aperture().relaxation().unwrap().dissipated_energy_j(),loss);
    b.extend_step_budget(4096).unwrap();
    for _ in 0..32 {assert_eq!(a.step(gentle).unwrap(),b.step(gentle).unwrap());}
}
