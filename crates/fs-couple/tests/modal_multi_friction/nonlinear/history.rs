//! Constitutive history shares the original single/multi-contact transaction.
use super::*;
use fs_couple::render::schedule::force::coupled::contact::friction::nonlinear::{
    MappedModalPotentialBinding, ModalCoordinate,
    history::{MappedModalHistoryBinding, ModalMaterialHistory},
};
use std::sync::{Arc, atomic::{AtomicBool, AtomicUsize, Ordering}};

// Manufactured Maxwell material, independent of the string implementation.
// Its deliberately corruptible heat receipt tests the HOST's work check after
// candidate history has already been written. These are not instrument data.
struct Maxwell {
    corrupt_heat:Arc<AtomicBool>,
    stages:Arc<AtomicUsize>,
    stiffness_bound:f64,
}
impl ModalMaterialHistory for Maxwell {
    fn coordinate_count(&self)->usize {2}
    fn memory_count(&self)->usize {2}
    fn stiffness_upper_bound_per_s2(&self,_:&[f64])->Result<f64,ModalCouplingError> {Ok(self.stiffness_bound)}
    fn relaxed_memory(&self,q:&[f64])->Result<Vec<f64>,ModalCouplingError> {Ok(q.to_vec())}
    fn energy(&self,q:&[f64],z:&[f64])->Result<f64,ModalCouplingError> {
        Ok(q.iter().zip(z).map(|(q,z)|0.5e6*(q-z).powi(2)).sum())
    }
    fn reaction(&self,q:&[f64],z:&[f64],out:&mut[f64])->Result<(),ModalCouplingError> {
        for ((out,q),z) in out.iter_mut().zip(q).zip(z) {*out=-1e6*(q-z);}
        Ok(())
    }
    fn stage(&self,a:&[f64],b:&[f64],z:&[f64],dt:f64,next:&mut[f64],out:&mut[f64])
        ->Result<f64,ModalCouplingError> {
        self.stages.fetch_add(1,Ordering::SeqCst);
        let rate=dt/0.001;let mut heat=0.0;
        for i in 0..2 {
            let strain=(0.5*(a[i]+b[i])-z[i])/(1.0+0.5*rate);
            next[i]=z[i]+rate*strain;out[i]=-1e6*strain;
            heat+=1e6*rate*strain*strain;
        }
        if self.corrupt_heat.load(Ordering::SeqCst) {heat+=1.0;}
        Ok(heat)
    }
}
fn core(corrupt_heat:Arc<AtomicBool>,stages:Arc<AtomicUsize>)->FrictionModalSystem {
    let gate=CancelGate::new_clock_free();
    let coordinates=vec![ModalCoordinate {component:0,mode:0},ModalCoordinate {component:1,mode:0}];
    system(true).with_mapped_nonlinear_potentials(vec![MappedModalPotentialBinding {
        coordinates:coordinates.clone(),replaces_linear:vec![],potential:Box::new(Quartic(1e8)),
    }],NonlinearModalConfig {max_iterations:64,maximum_generalized_force:1e5,
        force_absolute_tolerance:1e-11,force_relative_tolerance:1e-10},&gate).unwrap()
        .with_mapped_material_histories(vec![MappedModalHistoryBinding {coordinates,
            law:Box::new(Maxwell {corrupt_heat,stages,stiffness_bound:2e6}),initial_memory:Some(vec![1e-4,-2e-4])}],2,&gate).unwrap()
}
fn memory(s:&FrictionModalSystem)->Vec<Vec<f64>> {
    s.nonlinear_material_memory().into_iter().map(|m|m.to_vec()).collect()
}
fn joint_memory(s:&MultiFrictionModalSystem)->Vec<Vec<f64>> {
    s.nonlinear_material_memory().into_iter().map(|m|m.to_vec()).collect()
}

fn additional(bound:f64)->MappedModalHistoryBinding {
    MappedModalHistoryBinding {coordinates:vec![ModalCoordinate {component:0,mode:0},
        ModalCoordinate {component:1,mode:0}],initial_memory:None,
        law:Box::new(Maxwell {corrupt_heat:Arc::new(AtomicBool::new(false)),
            stages:Arc::new(AtomicUsize::new(0)),stiffness_bound:bound})}
}

#[test]
fn shared_band_sums_histories_and_reads_the_current_linear_operator() {
    let gate=CancelGate::new_clock_free();
    let raw=system(true).linear_frequency_upper_bound_rad_s();
    let actual=core(Arc::new(AtomicBool::new(false)),Arc::new(AtomicUsize::new(0)));
    let before=actual.linear_frequency_upper_bound_rad_s();
    assert!((before*before-raw*raw-2e6).abs()<1e-6);
    let original_memory=memory(&actual);
    let actual=actual.with_mapped_material_histories(vec![additional(2e6)],4,&gate).unwrap();
    assert!((actual.linear_frequency_upper_bound_rad_s().powi(2)-raw*raw-4e6).abs()<1e-6);
    assert_eq!(memory(&actual)[0],original_memory[0]);
    let branches=vec![fs_phs::RelaxationBranch {projection:vec![1.0,0.0,0.0,0.0],
        stiffness:3e6,relaxation_time_s:0.001}];
    let baseline=system(true).with_material_relaxation(branches.clone(),None,&gate).unwrap();
    let actual=actual.with_material_relaxation(branches,None,&gate).unwrap();
    let after=actual.linear_frequency_upper_bound_rad_s();
    assert!((after*after-baseline.linear_frequency_upper_bound_rad_s().powi(2)-4e6).abs()<1e-6);
    assert!(after>before);assert_eq!(memory(&actual)[0],original_memory[0]);
    let multiple=actual.with_additional_ports(vec![(ModalFrictionPort {
        left:attachment(1),right:None,law:law()},config())],joint_config(),&gate).unwrap();
    assert_eq!(multiple.linear_frequency_upper_bound_rad_s(),after);
    assert_eq!(system(true).linear_frequency_upper_bound_rad_s().to_bits(),raw.to_bits());
}

#[test]
fn absent_or_unresolved_history_stiffness_cannot_enter_the_owner() {
    for bound in [f64::NAN,f64::INFINITY,-1.0,1e20] {
        let actual=core(Arc::new(AtomicBool::new(false)),Arc::new(AtomicUsize::new(0)));
        assert!(actual.with_mapped_material_histories(vec![additional(bound)],4,
            &CancelGate::new_clock_free()).is_err());
    }
}

#[test]
fn single_contact_commits_history_once_and_rolls_back_a_late_bad_heat_receipt() {
    let bad=Arc::new(AtomicBool::new(false));let stages=Arc::new(AtomicUsize::new(0));
    let mut actual=core(bad.clone(),stages.clone());
    let mut control=core(Arc::new(AtomicBool::new(false)),Arc::new(AtomicUsize::new(0)));
    let initial=actual.total_energy_j().unwrap();let mut work=0.0;let mut loss=0.0;
    let forces=[0.02,-0.03];let bow=drive(0.2,0.01);
    for _ in 0..8 {
        let previous=memory(&actual);
        let frame=actual.step(&forces,bow).unwrap().clone();control.step(&forces,bow).unwrap();
        assert_ne!(memory(&actual),previous);
        assert_eq!(actual.last_frame(),control.last_frame());
        assert_eq!(memory(&actual),memory(&control));
        work+=frame.external_work_j+frame.drive_work_j;
        loss+=frame.network_dissipation_j+frame.friction_dissipation_j;
        assert!((actual.total_energy_j().unwrap()-initial+loss-work).abs()<1e-9);
    }
    let q=states(&actual);let z=memory(&actual);let last=actual.last_frame().cloned();
    let clock=actual.samples_rendered();let calls=stages.load(Ordering::SeqCst);
    bad.store(true,Ordering::SeqCst);
    assert!(actual.step(&forces,bow).is_err());
    assert!(stages.load(Ordering::SeqCst)>calls,"failure occurs after writing a trial history");
    assert_eq!(states(&actual),q);assert_eq!(memory(&actual),z);
    assert_eq!(actual.samples_rendered(),clock);assert_eq!(actual.last_frame(),last.as_ref());
    bad.store(false,Ordering::SeqCst);
    actual.step(&forces,bow).unwrap();control.step(&forces,bow).unwrap();
    assert_eq!(states(&actual),states(&control));assert_eq!(memory(&actual),memory(&control));
    assert_eq!(actual.last_frame(),control.last_frame());
}

#[test]
fn simultaneous_contacts_preserve_one_history_through_repeated_trials_and_failure() {
    let bad=Arc::new(AtomicBool::new(false));let stages=Arc::new(AtomicUsize::new(0));
    let joint=|core:FrictionModalSystem|core.with_additional_ports(vec![(ModalFrictionPort {
        left:attachment(1),right:None,law:law(),
    },config())],joint_config(),&CancelGate::new_clock_free()).unwrap();
    let mut actual=joint(core(bad.clone(),stages.clone()));
    let mut control=joint(core(Arc::new(AtomicBool::new(false)),Arc::new(AtomicUsize::new(0))));
    let initial=actual.total_energy_j().unwrap();let mut work=0.0;let mut loss=0.0;
    let forces=[0.02,-0.03];let bows=[drive(0.2,0.01),drive(-0.1,0.02)];
    for _ in 0..8 {
        let previous=joint_memory(&actual);
        let frame=actual.step(&forces,&bows).unwrap().clone();control.step(&forces,&bows).unwrap();
        assert_ne!(joint_memory(&actual),previous);
        assert_eq!(actual.last_frame(),control.last_frame());
        assert_eq!(joint_memory(&actual),joint_memory(&control));
        work+=frame.external_work_j+frame.drive_work_j;
        loss+=frame.network_dissipation_j+frame.friction_dissipation_j;
        assert!((actual.total_energy_j().unwrap()-initial+loss-work).abs()<1e-9);
    }
    let q=all_states(&actual);let z=joint_memory(&actual);let last=actual.last_frame().cloned();
    let clock=actual.samples_rendered();let calls=stages.load(Ordering::SeqCst);
    bad.store(true,Ordering::SeqCst);
    assert!(actual.step(&forces,&bows).is_err());
    assert!(stages.load(Ordering::SeqCst)>calls);
    assert_eq!(all_states(&actual),q);assert_eq!(joint_memory(&actual),z);
    assert_eq!(actual.samples_rendered(),clock);assert_eq!(actual.last_frame(),last.as_ref());
    bad.store(false,Ordering::SeqCst);
    actual.step(&forces,&bows).unwrap();control.step(&forces,&bows).unwrap();
    assert_eq!(all_states(&actual),all_states(&control));
    assert_eq!(joint_memory(&actual),joint_memory(&control));
    assert_eq!(actual.last_frame(),control.last_frame());
}
