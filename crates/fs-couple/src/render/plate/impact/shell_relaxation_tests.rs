use super::*;
use crate::render::plate::impact::{ImpactBody,ImpactConfig};
use fs_plate::{ModePair,PlateSection};
use fs_plate::shell::{ShellMesh,ShellSupport,assemble_shell_sections};
use fs_plate::shell::reduction::{ShellReduction,ReductionBudget};
use fs_exec::CancelGate;
fn shell()->ShellReduction {
    let mesh=ShellMesh::new(vec![[0.,0.,0.],[0.1,0.,0.01],[0.,0.1,0.02],[0.08,0.09,0.035]],
        vec![[0,1,3],[0,3,2]]).unwrap();
    let sections=vec![PlateSection::isotropic(112.6e9,0.342,0.001,8607.).unwrap(),
        PlateSection::isotropic(112.6e9,0.342,0.0007,8607.).unwrap()];
    let model=assemble_shell_sections(&mesh,&sections,&[0,1,2],ShellSupport::Clamped).unwrap();
    let n=model.free;let mut k=vec![0.;n*n];let mut m=k.clone();
    for i in 0..n {for j in 0..n {k[i*n+j]=model.k.get(i,j);m[i*n+j]=model.m.get(i,j);}}
    let modes:Vec<ModePair>=fs_modal::eigh_gen_dense(&k,&m,n).unwrap().into_iter().take(2).collect();
    ShellReduction::new(&mesh,&sections,&model,&modes,ReductionBudget {
        max_modes:12,max_facet_modes:100,relative_tolerance:1e-6}).unwrap()
}
fn make()->ImpactSystem {
    let s=shell();let n=s.mode_count();let mut initial=vec![crate::modal_acoustic_time::ModalAcousticState::default();n];
    initial[0].displacement_m_sqrt_kg=2e-8;initial[1].velocity_m_sqrt_kg_per_s=0.001;
    ImpactSystem::new(vec![ImpactBody {potential:BodyPotential::Shell(s),initial,damping_per_s:vec![0.;n]}],
        vec![],vec![],vec![],ImpactConfig{dt_s:1e-8,max_steps:1000,maximum_energy_j:20.,
            energy_absolute_tolerance_j:1e-12,energy_relative_tolerance:1e-7,maximum_generalized_force:1000.}).unwrap()
}
fn spectrum()->ShellBendingSpectrum {ShellBendingSpectrum {body:0,branches:vec![(0.15,1e-6),(0.05,1e-5)],band_hz:[0.,1e8]}}
#[test]
fn geometric_bending_memory_changes_motion_and_closes_total_work_loss() {
    let base=make();let s=shell();let kb=s.bending_stiffness(&[1.,1.],1000).unwrap();
    let x=base.state();let extra=0.1*(kb[0]*x[0]*x[0]+2.*kb[1]*x[0]*x[2]+kb[3]*x[2]*x[2]);
    let initial=base.stored_energy_j();
    let unrelaxed=make().with_shell_bending_relaxation(&[spectrum()],InitialMemory::Unrelaxed,256).unwrap();
    assert!((unrelaxed.stored_energy_j()-initial-extra).abs()<1e-15);
    let mut loaded=base.with_shell_bending_relaxation(&[spectrum()],InitialMemory::Relaxed,256).unwrap().prepare_analytic().unwrap();
    let mut bare=make().prepare_analytic().unwrap();
    assert_eq!(&loaded.state()[..4],bare.state());assert_eq!(loaded.relaxation_observation().unwrap().states,4);
    assert_eq!(loaded.relaxation_observation().unwrap().stored_energy_j,0.);
    let gate=CancelGate::new_clock_free();let mut balance=0.;let mut memory_loss=0.0_f64;
    for _ in 0..256 {
        let f=loaded.step(&[0.,0.],&gate).unwrap();bare.step(&[0.,0.],&gate).unwrap();balance+=f.dissipated_energy_j;
        memory_loss=memory_loss.max(loaded.relaxation_observation().unwrap().dissipated_power_w);
    }
    assert!((loaded.stored_energy_j()+balance-initial).abs()<1e-9);
    assert!(memory_loss>0.);assert!(balance>0.);
    assert!(loaded.state()[..4].iter().zip(bare.state()).any(|(a,b)|(a-b).abs()>1e-12));
}
#[test]
fn cancellation_and_acoustic_attachment_preserve_the_complete_material_state() {
    use crate::render::plate::impact::{radiation::{Model,Pole},ImpactSubstepConfig};
    let load=Model {ports:2,poles:vec![Pole{omega:10000.,zeta:0.1,coupling:vec![20.,-10.]}]};
    for air_first in [false,true] {
        let build=|| {
            let mut s=make();if air_first {s=s.with_radiation_load(&load,&[0,1],1).unwrap();}
            s=s.with_shell_bending_relaxation(&[spectrum()],InitialMemory::Relaxed,256).unwrap();
            if !air_first {s=s.with_radiation_load(&load,&[0,1],1).unwrap();}
            s.prepare_analytic().unwrap().with_substeps(ImpactSubstepConfig{max_depth:2,max_attempts:7}).unwrap()
        };
        let mut trial=build();let mut clean=build();let gate=CancelGate::new_clock_free();
        for _ in 0..16 {trial.step(&[0.,0.],&gate).unwrap();clean.step(&[0.,0.],&gate).unwrap();}
        let state=trial.state().to_vec();let memory=trial.relaxation_observation().unwrap();
        let cancelled=CancelGate::new_clock_free();cancelled.request();assert!(trial.step(&[0.,0.],&cancelled).is_err());
        assert_eq!(trial.state(),state);assert_eq!(trial.relaxation_observation().unwrap(),memory);
        trial.step(&[0.,0.],&gate).unwrap();clean.step(&[0.,0.],&gate).unwrap();assert_eq!(trial.state(),clean.state());
    }
}
#[test]
fn invalid_spectra_and_work_limits_refuse_without_substituting_a_material() {
    let mut bad=spectrum();bad.body=1;assert!(make().with_shell_bending_relaxation(&[bad],InitialMemory::Relaxed,256).is_err());
    let mut bad=spectrum();bad.branches[0].0=-1.;assert!(make().with_shell_bending_relaxation(&[bad],InitialMemory::Relaxed,256).is_err());
    let mut bad=spectrum();bad.branches[0].1=1e-12;assert!(make().with_shell_bending_relaxation(&[bad],InitialMemory::Relaxed,256).is_err());
    let mut bad=spectrum();bad.band_hz=[0.,1.];assert!(make().with_shell_bending_relaxation(&[bad],InitialMemory::Relaxed,256).is_err());
    assert!(make().with_shell_bending_relaxation(&[spectrum(),spectrum()],InitialMemory::Relaxed,256).is_err());
    assert!(make().with_shell_bending_relaxation(&[spectrum()],InitialMemory::Relaxed,3).is_err());
    let mut zero=spectrum();zero.branches=vec![(0.,1e-6)];
    let mut a=make().prepare_analytic().unwrap();let mut b=make().with_shell_bending_relaxation(&[zero],InitialMemory::Relaxed,256).unwrap().prepare_analytic().unwrap();
    let gate=CancelGate::new_clock_free();for _ in 0..16 {a.step(&[0.,0.],&gate).unwrap();b.step(&[0.,0.],&gate).unwrap();assert_eq!(a.state(),b.state());}
}
