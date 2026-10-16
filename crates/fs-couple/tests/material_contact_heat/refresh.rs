//! In-flight source changes must reach energy, force, tangent and retained memory.
use super::*;
use fs_couple::render::plate::impact::material::refresh::ContactSourceUpdate;
use fs_couple::render::plate::impact::relaxation::InitialMemory;
use fs_phs::RelaxationBranch;

fn update(young:f64,chi:Option<f64>)->ContactSourceUpdate {
    ContactSourceUpdate {contact:"touch".into(),source:crate::common::source(young,chi)}
}

#[test]
fn unloaded_refresh_changes_the_next_impact_without_resetting_motion_or_clock() {
    let gate=CancelGate::new_clock_free();
    for analytic in [false,true] {
        let mut live=collision(Some(0.4),0.0,false).prepare().unwrap();
        let mut expected=collision(Some(0.4),0.0,false);
        expected.refresh_unloaded_contacts(vec![update(3e7,Some(0.4))],&gate).unwrap();
        let mut expected=expected.prepare().unwrap();
        live.set_analytic_newton(analytic);expected.set_analytic_newton(analytic);
        for _ in 0..100 {live.step(&[0.0;2],&gate).unwrap();expected.step(&[0.0;2],&gate).unwrap();}
        assert_eq!(live.state(),expected.state());
        let x=live.state().to_vec();let h=live.stored_energy_j();let count=live.samples();
        let map=live.material_contacts()[0].weights().to_vec();
        assert_eq!(live.refresh_unloaded_contacts(vec![update(3e7,Some(0.4))],&gate).unwrap(),1);
        assert_eq!(live.state(),x);assert_eq!(live.samples(),count);assert_eq!(live.stored_energy_j(),h);
        assert_eq!(live.material_contacts()[0].weights(),map);
        let mut peak=0.0_f64;
        for _ in 0..1500 {
            let a=live.step(&[0.0;2],&gate).unwrap();let b=expected.step(&[0.0;2],&gate).unwrap();
            assert_eq!(live.state(),expected.state());assert_eq!(a.stored_energy_j,b.stored_energy_j);
            assert_eq!(a.dissipated_energy_j,b.dissipated_energy_j);
            peak=peak.max(live.material_contact_observation(0).unwrap().normal_force_n);
        }
        assert!(peak>0.05);assert_eq!(live.material_contacts(),expected.material_contacts());
    }
}

#[test]
fn source_refresh_reaches_wrapped_storage_and_keeps_relaxation_history() {
    let gate=CancelGate::new_clock_free();
    let attach=|s:ImpactSystem|s.with_relaxation_branches(vec![RelaxationBranch {
        projection:vec![1.0,0.0,-1.0,0.0],stiffness:10.0,relaxation_time_s:0.01,
    }],InitialMemory::Unrelaxed,1).unwrap();
    let mut original=collision(Some(0.4),0.0,false);
    original.refresh_unloaded_contacts(vec![update(3e7,Some(0.4))],&gate).unwrap();
    let mut expected=attach(original).prepare_analytic().unwrap();
    let mut live=attach(collision(Some(0.4),0.0,false)).prepare_analytic().unwrap();
    for _ in 0..100 {live.step(&[0.0;2],&gate).unwrap();expected.step(&[0.0;2],&gate).unwrap();}
    let before=live.state().to_vec();let energy=live.stored_energy_j();let memory=live.relaxation_observation();
    live.refresh_unloaded_contacts(vec![update(3e7,Some(0.4))],&gate).unwrap();
    assert_eq!(live.state(),before);assert_eq!(live.stored_energy_j(),energy);
    assert_eq!(live.relaxation_observation(),memory);
    for _ in 0..1400 {
        live.step(&[0.0;2],&gate).unwrap();expected.step(&[0.0;2],&gate).unwrap();
        assert_eq!(live.state(),expected.state());assert_eq!(live.stored_energy_j(),expected.stored_energy_j());
    }
}

#[test]
fn unloaded_refresh_can_activate_contact_loss_in_reference_and_prepared_equations() {
    let gate=CancelGate::new_clock_free();
    let mut live=collision(None,0.0,false);let mut expected=collision(Some(0.4),0.0,false);
    for _ in 0..100 {live.step(&[0.0;2],&gate).unwrap();expected.step(&[0.0;2],&gate).unwrap();}
    live.refresh_unloaded_contacts(vec![update(1e7,Some(0.4))],&gate).unwrap();
    let mut loss=0.0;
    for _ in 0..1300 {
        loss+=live.step(&[0.0;2],&gate).unwrap().dissipated_energy_j;
        expected.step(&[0.0;2],&gate).unwrap();assert_eq!(live.state(),expected.state());
    }
    assert!(loss>1e-9);
    // Reference/prepared conversion retains the current source and all state.
    let mut live=live.prepare_analytic().unwrap();let mut expected=expected.prepare_analytic().unwrap();
    for _ in 0..100 {live.step(&[0.0;2],&gate).unwrap();expected.step(&[0.0;2],&gate).unwrap();}
    assert_eq!(live.state(),expected.state());assert_eq!(live.material_contacts(),expected.material_contacts());
}

#[test]
fn loaded_contact_cannot_be_reparameterized_as_free_energy() {
    let gate=CancelGate::new_clock_free();let mut s=collision(Some(0.4),0.0,false).prepare_analytic().unwrap();
    while s.material_contact_observation(0).unwrap().approach_m==0.0 {s.step(&[0.0;2],&gate).unwrap();}
    let x=s.state().to_vec();let h=s.stored_energy_j();let n=s.samples();let source=s.material_contacts().to_vec();
    assert!(s.refresh_unloaded_contacts(vec![update(3e7,Some(0.4))],&gate).is_err());
    assert_eq!(s.state(),x);assert_eq!(s.stored_energy_j(),h);assert_eq!(s.samples(),n);
    assert_eq!(s.material_contacts(),source);
}

#[test]
fn failed_batch_geometry_rate_or_cancellation_preserves_every_source() {
    let gate=CancelGate::new_clock_free();let mut s=collision(Some(0.4),0.0,false).prepare_analytic().unwrap();
    let x=s.state().to_vec();let sources=s.material_contacts().to_vec();
    let good=update(3e7,Some(0.4));let mut unknown=good.clone();unknown.contact="missing".into();
    assert!(s.refresh_unloaded_contacts(vec![good.clone(),unknown],&gate).is_err());
    assert!(s.refresh_unloaded_contacts(vec![good.clone(),good.clone()],&gate).is_err());
    let mut geometry=good.clone();geometry.source=SourcePointContact::new(good.source.model().clone(),
        PointContactGeometry::SpherePlane {effective_radius_m:0.2},good.source.envelope()).unwrap();
    assert!(s.refresh_unloaded_contacts(vec![geometry],&gate).is_err());
    let mut bounded=good.clone();let mut limit=good.source.envelope();limit.maximum_absolute_rate_m_per_s=0.01;
    bounded.source=SourcePointContact::new(good.source.model().clone(),good.source.geometry(),limit).unwrap();
    assert!(s.refresh_unloaded_contacts(vec![bounded],&gate).is_err());
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(s.refresh_unloaded_contacts(vec![good.clone()],&cancelled).is_err());
    assert_eq!(s.state(),x);assert_eq!(s.material_contacts(),sources);assert_eq!(s.samples(),0);
    s.refresh_unloaded_contacts(vec![good],&gate).unwrap();s.step(&[0.0;2],&gate).unwrap();
    assert_eq!(s.samples(),1);
}

#[path="feedback.rs"]
mod feedback;
