//! Actual contact work becomes material heat without a second collision solver.
#[path="../../fs-contact/tests/source_point/common.rs"]
mod common;
use common::*;
use fs_blake3::ContentHash;
use fs_contact::normal_patch::source::{PointContactGeometry, SourcePointContact};
use fs_couple::render::plate::impact::{ImpactBody,ImpactConfig,ImpactSystem};
use fs_couple::render::plate::impact::material::{MaterialContact,heating::{
    CaloricReservoir,ContactHeatRoute,HeatedImpactSystem}};
use fs_exec::CancelGate;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve,EnthalpyPhaseKnot,SolidLiquidPhase};

fn config()->ImpactConfig {ImpactConfig {dt_s:2e-6,max_steps:3000,maximum_energy_j:1.0,
    energy_absolute_tolerance_j:1e-12,energy_relative_tolerance:1e-7,maximum_generalized_force:1.0}}
fn collision(chi:Option<f64>, damping:f64, second_contact:bool)->ImpactSystem {
    let a=solid("A",1000.0,1e7);let b=solid("B",2000.0,2e7);
    let (mut left,wa)=ImpactBody::free_mass(0.001,-2e-5,0.04).unwrap();
    left.damping_per_s[0]=damping;
    let (right,wb)=ImpactBody::free_mass(0.002,0.0,0.0).unwrap();
    let source=SourcePointContact::new(bind(&a,&b,chi,0.0),
        PointContactGeometry::SpherePlane {effective_radius_m:0.1},envelope()).unwrap();
    let mut contacts=vec![MaterialContact::new("touch".into(),source.clone(),vec![wa,-wb],0.0).unwrap()];
    if second_contact {contacts.push(MaterialContact::new("other".into(),source,vec![wa,-wb],1e-6).unwrap());}
    ImpactSystem::new_with_material_contacts(vec![left,right],contacts,vec![],vec![],config()).unwrap()
}
fn curve(cp:f64,upper:f64)->EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(ContentHash([7;32]),SolidLiquidPhase::Solid,vec![
        EnthalpyPhaseKnot {specific_enthalpy_j_kg:0.0,temperature_k:300.0,liquid_mass_fraction:0.0,bulk_density_kg_m3:1000.0},
        EnthalpyPhaseKnot {specific_enthalpy_j_kg:upper,temperature_k:300.0+upper/cp,liquid_mass_fraction:0.0,bulk_density_kg_m3:1000.0},
    ]).unwrap()
}
fn stores(c:&EquilibriumEnthalpyPhaseCurve)->Vec<CaloricReservoir<'_>> {
    ["left","right"].into_iter().map(|key|CaloricReservoir::new(key.into(),c,
        c.state_at_specific_enthalpy(0.0).unwrap(),0.001,1e-6,0.0006).unwrap()).collect()
}
fn route(contact:&str)->ContactHeatRoute {ContactHeatRoute {contact:contact.into(),
    to_a:"left".into(),to_b:"right".into(),fraction_to_a:0.3}}
fn near(a:f64,b:f64,tol:f64) {assert!((a-b).abs()<=tol,"{a} vs {b}");}

#[test]
fn actual_step_loss_heats_both_material_stores_without_changing_collision() {
    let c=curve(0.02,1.0);let gate=CancelGate::new_clock_free();
    for analytic in [false,true] {
        let mut plain=collision(Some(0.4),0.0,false).prepare().unwrap();
        plain.set_analytic_newton(analytic);
        let mut heated=HeatedImpactSystem::new(collision(Some(0.4),0.0,false).prepare().unwrap(),
            stores(&c),vec![route("touch")],1e-12).unwrap();
        heated.set_analytic_newton(analytic);
        let initial=plain.stored_energy_j();let mut deposited=0.0;
        for _ in 0..2000 {
            let expected=plain.step(&[0.0;2],&gate).unwrap();let f=heated.step(&[0.0;2],&gate).unwrap();
            assert_eq!(plain.state(),heated.mechanics().state());
            assert_eq!(expected.dissipated_energy_j,f.mechanical.dissipated_energy_j);
            near(f.deposited_heat_j,f.mechanical.dissipated_energy_j,1e-18);
            assert!(f.combined_energy_residual_j.abs()<1e-10);deposited+=f.deposited_heat_j;
        }
        assert!(deposited>1e-9);
        let r=heated.reservoirs();
        near(r[0].received_heat_j(),0.3*deposited,1e-18);
        near(r[1].received_heat_j(),0.7*deposited,1e-18);
        near(r[0].state().temperature_k()-300.0,0.3*deposited/(0.001*0.02),1e-10);
        near(heated.mechanics().stored_energy_j()+r.iter().map(|s|s.storage_change_j()).sum::<f64>(),initial,1e-9);
        assert!(std::ptr::eq(r[0].curve(),&c));
    }
}

#[test]
fn unrelated_viscous_loss_is_not_misassigned_to_the_contact_or_counted_twice() {
    let c=curve(1.0,1.0);let gate=CancelGate::new_clock_free();
    let mut heated=HeatedImpactSystem::new(collision(Some(0.4),20.0,false).prepare_analytic().unwrap(),
        stores(&c),vec![route("touch")],1e-12).unwrap();
    let initial=heated.mechanics().stored_energy_j();let mut remaining=0.0;let mut contact=0.0;
    for _ in 0..2000 {let f=heated.step(&[0.0;2],&gate).unwrap();
        remaining+=f.remaining_mechanical_loss_j;contact+=f.deposited_heat_j;}
    assert!(remaining>1e-9 && contact>1e-10);
    near(heated.mechanics().stored_energy_j()+heated.reservoirs().iter().map(|s|s.storage_change_j()).sum::<f64>()+remaining,initial,1e-9);
}

#[test]
fn multiple_contacts_share_stores_in_canonical_order_and_retained_parts_resume() {
    let c=curve(1.0,1.0);let gate=CancelGate::new_clock_free();
    let make=|routes|HeatedImpactSystem::new(collision(Some(0.4),0.0,true).prepare_analytic().unwrap(),stores(&c),routes,1e-12).unwrap();
    let mut a=make(vec![route("touch"),route("other")]);
    let mut b=make(vec![route("other"),route("touch")]);
    for _ in 0..800 {a.step(&[0.0;2],&gate).unwrap();b.step(&[0.0;2],&gate).unwrap();}
    assert_eq!(a.reservoirs(),b.reservoirs());assert_eq!(a.mechanics().state(),b.mechanics().state());
    let (mechanics,stores)=a.into_parts();
    let mut a=HeatedImpactSystem::new(mechanics,stores,vec![route("touch"),route("other")],1e-12).unwrap();
    for _ in 0..400 {a.step(&[0.0;2],&gate).unwrap();b.step(&[0.0;2],&gate).unwrap();}
    assert_eq!(a.reservoirs(),b.reservoirs());assert_eq!(a.mechanics().samples(),b.mechanics().samples());
}

#[test]
fn caloric_latent_storage_changes_fraction_not_temperature() {
    // Explicit thermal absorber, NOT a claim that the rigid bodies melt.
    let c=EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([8;32]),vec![
        EnthalpyPhaseKnot {specific_enthalpy_j_kg:0.0,temperature_k:300.0,liquid_mass_fraction:0.0,bulk_density_kg_m3:1000.0},
        EnthalpyPhaseKnot {specific_enthalpy_j_kg:1.0,temperature_k:300.0,liquid_mass_fraction:1.0,bulk_density_kg_m3:900.0},
    ]).unwrap();
    let mut heated=HeatedImpactSystem::new(collision(Some(0.4),0.0,false).prepare_analytic().unwrap(),stores(&c),vec![route("touch")],1e-12).unwrap();
    for _ in 0..2000 {heated.step(&[0.0;2],&CancelGate::new_clock_free()).unwrap();}
    for s in heated.reservoirs() {assert_eq!(s.state().temperature_k(),300.0);
        assert!(s.state().liquid_mass_fraction()>0.0);
        near(s.state().liquid_mass_fraction(),s.received_heat_j()/s.mass_kg(),1e-15);}
}

#[test]
fn exhausted_thermal_domain_refuses_before_motion_and_other_store_publication() {
    let big=curve(1.0,1.0);let small=curve(1e-10,1e-12);let gate=CancelGate::new_clock_free();
    let mut rs=stores(&big);rs[1]=CaloricReservoir::new("right".into(),&small,
        small.state_at_specific_enthalpy(0.0).unwrap(),0.001,1e-6,0.0006).unwrap();
    let mut heated=HeatedImpactSystem::new(collision(Some(0.4),0.0,false).prepare_analytic().unwrap(),rs,vec![route("touch")],1e-12).unwrap();
    let mut plain=collision(Some(0.4),0.0,false).prepare_analytic().unwrap();let mut refused=false;
    for _ in 0..2000 {
        let before=heated.reservoirs().to_vec();let motion=heated.mechanics().state().to_vec();let sample=heated.mechanics().samples();
        if heated.step(&[0.0;2],&gate).is_err() {
            assert_eq!(heated.reservoirs(),before);assert_eq!(heated.mechanics().state(),motion);
            assert_eq!(heated.mechanics().samples(),sample);refused=true;break;
        }
        plain.step(&[0.0;2],&gate).unwrap();
    }
    assert!(refused);assert_eq!(heated.mechanics().state(),plain.state());
}

#[test]
fn elastic_cancellation_solver_refusal_and_route_validation_do_not_create_heat() {
    let c=curve(1.0,1.0);let gate=CancelGate::new_clock_free();
    let make=|routes|HeatedImpactSystem::new(collision(None,0.0,false).prepare_analytic().unwrap(),stores(&c),routes,1e-12);
    assert!(make(vec![route("touch"),route("touch")]).is_err());
    assert!(make(vec![route("unknown")]).is_err());
    let mut bad=route("touch");bad.fraction_to_a=f64::NAN;assert!(make(vec![bad]).is_err());
    let mut bad=route("touch");bad.to_a="missing".into();assert!(make(vec![bad]).is_err());
    let mut heated=make(vec![route("touch")]).unwrap();let saved=heated.reservoirs().to_vec();
    for _ in 0..600 {heated.step(&[0.0;2],&gate).unwrap();}
    assert_eq!(heated.reservoirs(),saved);
    let motion=heated.mechanics().state().to_vec();let sample=heated.mechanics().samples();
    heated.set_iteration_limit(0).unwrap();assert!(heated.step(&[0.0;2],&gate).is_err());
    let cancelled=CancelGate::new_clock_free();cancelled.request();assert!(heated.step(&[0.0;2],&cancelled).is_err());
    assert!(heated.step(&[f64::NAN,0.0],&gate).is_err());
    assert_eq!(heated.reservoirs(),saved);assert_eq!(heated.mechanics().state(),motion);assert_eq!(heated.mechanics().samples(),sample);
}
