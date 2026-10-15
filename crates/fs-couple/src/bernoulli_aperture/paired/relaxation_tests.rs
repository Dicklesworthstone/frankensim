//! Constitutive and coupling regressions; no measured-instrument timbre claim.
use super::*;
use crate::bernoulli_aperture::dynamic::relaxation::PlateRelaxationRegion;
use crate::bernoulli_aperture::plate::PlateApertureOptions;
use fs_exec::CancelGate;
use fs_material::{gas::GasSpec, visco::GeneralizedMaxwell};
use fs_plate::{AssemblyOptions, EdgeSupport, PlateChart, PlateMesh, PlateSection, SliceOptions};

const DT: f64 = 1e-5;
fn plate(e: f64, density: f64) -> PlateApertureReduction {
    let mesh=PlateMesh::rectangle(0.025,0.01,6,2);
    let root=mesh.nodes.iter().enumerate().filter(|(_,p)|p.0==0.0).map(|(i,_)|i).collect();
    let edges=mesh.boundary_edges().into_iter().filter(|&(a,b)|
        (mesh.nodes[a].0-0.025).abs()<1e-12 && (mesh.nodes[b].0-0.025).abs()<1e-12)
        .map(|(a,b)|[a,b]).collect();
    let section=PlateSection::isotropic(e,0.3,0.0003,density).unwrap();
    let chart=PlateChart::with_boundary_and_regions(mesh,section,root,vec![]).unwrap();
    let options=PlateApertureOptions {
        assembly:AssemblyOptions {pretension:0.0,support:EdgeSupport::Clamped},
        eigenvalue_window:(1.0,(core::f64::consts::TAU*450.0).powi(2)),mode_index:0,
        eigensolver:SliceOptions::default(),max_nodes:100,max_triangles:100,
        slit_edges:edges,rest_opening_m:0.0002,damping_ratio:0.0,
        max_slit_mode_variation:0.25,max_slope:0.1,
    };
    PlateApertureReduction::from_chart(chart,options,&CancelGate::new_clock_free()).unwrap()
}
fn law(p: &PlateApertureReduction, e: f64, terms: Vec<(f64,f64)>) -> PlateRelaxationSpec {
    PlateRelaxationSpec {regions:vec![PlateRelaxationRegion {
        triangles:(0..p.chart().mesh.tris.len()).collect(),
        material:GeneralizedMaxwell::new(e,terms).unwrap(),poisson_ratio:0.3,
        band_hz:(0.0,1000.0),provenance:"synthetic independent-blade material law".into(),
    }],max_branches:64,max_dt_over_tau:0.1,max_angular_step:0.1}
}
fn base(plates: [PlateApertureReduction;2], states: [BladeState;2]) -> PairedAperture {
    let gas=GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,101325.0).unwrap();
    let contact=Obstacle::new(vec![-1.0],1,1,vec![0.0],vec![1.0],0.0,2.0,
        "explicit no-contact material comparison".into()).unwrap();
    PairedAperture::from_plates(plates,gas,1e6,DT,4096,
        PairedApertureLimits {maximum_abs_displacement_m:0.001,maximum_abs_velocity_m_s:10.0,
            maximum_energy_j:1.0,max_evaluations:256,velocity_tolerance_m_s:1e-15,
            relative_tolerance:1e-13},states,contact).unwrap()
}
fn material_pair() -> PairedAperture {
    let plates=[plate(4e9,900.0),plate(8e9,1300.0)];
    let laws=[law(&plates[0],4e9,vec![(2e9,0.001),(1e9,0.01)]),
        law(&plates[1],8e9,vec![(3e9,0.003)])];
    base(plates,[BladeState {displacement_m:-3e-6,velocity_m_s:0.0},
        BladeState {displacement_m:2e-6,velocity_m_s:0.0}])
        .with_blade_relaxation(0,laws[0].clone(),InitialApertureMemory::Unrelaxed).unwrap()
        .with_blade_relaxation(1,laws[1].clone(),InitialApertureMemory::ViscousDisplacementM(vec![1e-6])).unwrap()
}
fn near(a: f64, b: f64, scale: f64) {
    assert!((a-b).abs()<=3e-8*scale.max(1e-20),"{a:e} != {b:e}");
}

#[test]
fn condensed_material_response_is_the_existing_owner_not_a_lagged_damping_fit() {
    let model=material_pair();
    for i in 0..2 {
        let m=model.blade_relaxation(i).unwrap();let old=model.blades[i].displacement_m;
        let (bias,resistance)=m.midpoint_response(old,DT).unwrap();
        assert!(resistance>0.0);
        for velocity in [-0.4,-0.001,0.0,0.02,0.4] {
            let direct=m.restoring_force(old,old+DT*velocity,DT).unwrap();
            near(bias+resistance*velocity,direct,bias.abs()+resistance*velocity.abs());
        }
    }
    assert_eq!(model.blade_relaxation(0).unwrap().branches().len(),2);
    assert_eq!(model.blade_relaxation(1).unwrap().branches().len(),1);
}

#[test]
fn both_independent_histories_match_full_phs_mechanics_and_the_same_fluid_equation() {
    let mut model=material_pair();let s=*model.spec();
    let reference:Vec<_>=(0..2).map(|i| {
        let b=s.blades[i];
        let arms=model.blade_relaxation(i).unwrap().branches().iter().map(|a|fs_phs::RelaxationBranch {
            projection:vec![1.0,0.0],stiffness:a.stiffness,relaxation_time_s:a.relaxation_time_s,
        }).collect();
        fs_phs::mass_spring_damper(b.mass_kg,b.stiffness_n_m,0.0).unwrap()
            .with_relaxation_branches(arms).unwrap()
    }).collect();
    let mut states:Vec<Vec<f64>>=(0..2).map(|i| {
        let mut x=vec![model.blades[i].displacement_m,s.blades[i].mass_kg*model.blades[i].velocity_m_s];
        x.extend_from_slice(model.blade_relaxation(i).unwrap().memory_sqrt_j());x
    }).collect();
    for n in 0..128 {
        let dp=if n<64 {2.0_f64} else {-1.0};
        let forces=if n<80 {[1e-4,-2e-4]} else {[-1e-4,3e-4]};
        let results:Vec<_>=(0..2).map(|i|fs_phs::step(&reference[i],&states[i],
            &[forces[i]-s.blades[i].pressure_area_m2*dp],DT).unwrap()).collect();
        let velocities:[f64;2]=std::array::from_fn(|i|
            f64::midpoint(states[i][1],results[i].x[1])/s.blades[i].mass_kg);
        let old_gap=s.rest_opening_m+states[0][0]+states[1][0];
        let next_gap=s.rest_opening_m+results[0].x[0]+results[1].x[0];
        let jet=s.width_m*f64::midpoint(old_gap,next_gap).max(0.0)*dp.signum()
            *(2.0*dp.abs()/s.density_kg_m3).sqrt();
        let swept=-(0..2).map(|i|s.blades[i].pressure_area_m2*velocities[i]).sum::<f64>();
        let incoming=0.3;let body=1e-8;
        let frame=model.step_with_forces(ApertureDrive {
            upstream_pressure_pa:dp+2.0*incoming+s.impedance_pa_s_m3*(jet+swept+body),
            incoming_pressure_pa:incoming,body_flow_m3_s:body,
        },forces).unwrap();
        near(frame.pressure_drop_pa,dp,dp.abs());
        let mut loss=0.0;let mut energy=0.0;
        for i in 0..2 {
            near(frame.blades[i].displacement_m,results[i].x[0],s.rest_opening_m);
            near(frame.blades[i].velocity_m_s,results[i].x[1]/s.blades[i].mass_kg,0.1);
            assert!(frame.momentum_residual_n[i].abs()<1e-9);
            let m=model.blade_relaxation(i).unwrap();
            for (&z,&expected) in m.memory_sqrt_j().iter().zip(&results[i].x[2..]) {
                near(z,expected,1e-3);
            }
            energy+=m.stored_energy_j();loss+=results[i].dissipated;
            states[i]=results[i].x.clone();
        }
        assert_eq!(frame.aperture.relaxation_energy_j,energy);
        near(frame.aperture.relaxation_loss_j,loss,frame.aperture.stored_energy_j+loss);
        let scale=frame.aperture.stored_energy_j+frame.aperture.storage_change_j.abs()
            +frame.aperture.dissipated_energy_j+frame.aperture.pressure_work_j.abs()
            +frame.aperture.mechanical_work_j.abs();
        near(frame.aperture.balance_residual_j(),0.0,scale);
        assert!(frame.aperture.relaxation_loss_j>=0.0);
    }
}

#[test]
fn a_one_blade_material_change_alters_pressure_without_retuning_equilibrium_mechanics() {
    let p=plate(4e9,900.0);
    let short=law(&p,4e9,vec![(2e9,0.001)]);
    let long=law(&p,4e9,vec![(2e9,0.02)]);
    let build=|s|base([p.clone(),p.clone()],[BladeState::default();2])
        .with_blade_relaxation(1,s,InitialApertureMemory::Relaxed).unwrap();
    let mut a=build(short);let mut b=build(long);
    assert_eq!(a.spec(),b.spec()); // No E_inf, mass, area, density or clock retune.
    assert!(a.blade_relaxation(0).is_none());
    let mut changed=false;
    for n in 0..256 {
        let d=ApertureDrive {upstream_pressure_pa:if n<128 {5.0}else{0.0},..ApertureDrive::default()};
        let fa=a.step(d).unwrap();let fb=b.step(d).unwrap();
        changed|=(fa.aperture.bore_pressure_pa-fb.aperture.bore_pressure_pa).abs()>1e-10;
    }
    assert!(changed);
    assert!(a.blade_relaxation(1).unwrap().dissipated_energy_j()>0.0);
    assert!(b.blade_relaxation(1).unwrap().dissipated_energy_j()>0.0);
}

#[test]
fn zero_arm_laws_preserve_the_original_elastic_trajectory_exactly() {
    let p=plate(4e9,900.0);let empty=law(&p,4e9,vec![]);
    let mut a=base([p.clone(),p.clone()],[BladeState::default();2]);
    let mut b=base([p.clone(),p],[BladeState::default();2])
        .with_blade_relaxation(0,empty.clone(),InitialApertureMemory::Relaxed).unwrap()
        .with_blade_relaxation(1,empty,InitialApertureMemory::Unrelaxed).unwrap();
    for n in 0..128 {
        let d=ApertureDrive {upstream_pressure_pa:if n<64 {3.0}else{-1.0},..ApertureDrive::default()};
        assert_eq!(a.step(d).unwrap(),b.step(d).unwrap());
    }
}

#[test]
fn immutable_trials_and_refusals_preserve_both_histories_and_energy() {
    let mut model=material_pair();
    let before=*model.blade_states();let energy=model.stored_energy_j();
    let history:Vec<_>=(0..2).map(|i|model.blade_relaxation(i).unwrap().memory_sqrt_j().to_vec()).collect();
    let drive=ApertureDrive {upstream_pressure_pa:3.0,..ApertureDrive::default()};
    let trial=model.preview(drive,[1e-4,-1e-4]).unwrap();
    assert_eq!(model.blade_states(),&before);assert_eq!(model.stored_energy_j(),energy);
    assert!(model.step_with_forces(drive,[f64::NAN,0.0]).is_err());
    for i in 0..2 {
        assert_eq!(model.blade_relaxation(i).unwrap().memory_sqrt_j(),history[i]);
        assert_eq!(model.blade_relaxation(i).unwrap().dissipated_energy_j(),0.0);
    }
    let expected=trial.frame;model.accept(trial);
    assert_eq!(model.accepted_steps(),1);
    assert_eq!(model.stored_energy_j(),expected.aperture.stored_energy_j);
    let history:Vec<_>=(0..2).map(|i|model.blade_relaxation(i).unwrap().memory_sqrt_j().to_vec()).collect();
    model.spec.max_steps=1;
    assert!(model.step(drive).is_err());
    for i in 0..2 {assert_eq!(model.blade_relaxation(i).unwrap().memory_sqrt_j(),history[i]);}
    model.extend_step_budget(2).unwrap();model.step(drive).unwrap();
}

#[test]
fn source_clock_damping_history_and_combined_energy_admission_are_not_bypassed() {
    let p=plate(4e9,900.0);let good=law(&p,4e9,vec![(2e9,0.001)]);
    let build=||base([p.clone(),p.clone()],[BladeState::default();2]);
    assert!(build().with_blade_relaxation(2,good.clone(),InitialApertureMemory::Relaxed).is_err());
    let mut bad=good.clone();bad.regions[0].triangles.pop();
    assert!(build().with_blade_relaxation(0,bad,InitialApertureMemory::Relaxed).is_err());
    let mut bad=good.clone();bad.regions[0].material.e_inf=8e9;
    assert!(build().with_blade_relaxation(0,bad,InitialApertureMemory::Relaxed).is_err());
    let mut bad=good.clone();bad.regions[0].material.terms[0].1=1e-12;
    assert!(build().with_blade_relaxation(1,bad,InitialApertureMemory::Relaxed).is_err());
    assert!(build().with_blade_relaxation(0,good.clone(),InitialApertureMemory::ViscousDisplacementM(vec![])).is_err());
    let mut damped=build();damped.spec.blades[1].damping_n_s_m=0.01;
    assert!(damped.with_blade_relaxation(1,good.clone(),InitialApertureMemory::Relaxed).is_err());
    let once=build().with_blade_relaxation(0,good.clone(),InitialApertureMemory::Relaxed).unwrap();
    assert!(once.with_blade_relaxation(0,good.clone(),InitialApertureMemory::Relaxed).is_err());
    let mut late=build();late.step(ApertureDrive::default()).unwrap();
    assert!(late.with_blade_relaxation(1,good.clone(),InitialApertureMemory::Relaxed).is_err());
    let mut limited=base([p.clone(),p],[BladeState {displacement_m:1e-5,velocity_m_s:0.0};2]);
    limited.spec.limits.maximum_energy_j=limited.stored_energy_j()*1.01;
    assert!(limited.with_blade_relaxation(0,good,InitialApertureMemory::Unrelaxed).is_err());
}
