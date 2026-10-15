//! Caloric data changes a later impact through the same sourced contact ladder.
use super::*;
use fs_contact::normal_patch::source::state::PointContactMaterialCards;
use fs_contact::interface_binding::NormalContactModelSelection;
use fs_couple::render::plate::impact::material::heating::feedback::resolve_isothermal_contact;
use fs_couple::render::plate::impact::material::heating::transport::{CoupledThermalStep,ImpactThermalTrial};
use fs_evidence::ValidityDomain;
use fs_matdb::{MaterialCard,MaterialStateId,InterfaceSystemCard,ClaimSet,PropertyClaim,PropertyKey,
    PropertyValue,Provenance,UncertaintyModel,InterpolationPolicy,QueryPoint};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Density,Pressure,Temperature,Dims};
use fs_tribo::InputAuthority;

fn material(cp:f64)->MaterialCard {
    let mut claims=ClaimSet::new();
    for (name,value,dims) in [("density",1000.0,Density::DIMS),("young_modulus",1e7,Pressure::DIMS),
        ("poisson_ratio",0.25,Dims::NONE),("yield_stress",1e10,Pressure::DIMS),
        ("specific_heat_capacity",cp,Dims([2,0,-2,-1,0,0]))] {
        let young=name=="young_modulus";
        claims.insert_claim(PropertyClaim {key:PropertyKey::new(name,dims),
            value:if young {PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
                knots:vec![(280.0,1e7),(320.0,4e7)],dims}}else{PropertyValue::Scalar {value,dims}},
            validity:ValidityDomain::unconstrained().with("T",280.0,320.0),
            interpolation:if young{InterpolationPolicy::LinearInside}else{InterpolationPolicy::ConstantWithinValidity},
            uncertainty:UncertaintyModel::Unstated,observations:vec![],
            provenance:Provenance {source:"synthetic temperature-contact fixture, not measured material".into(),
                license:"CC0-1.0".into(),artifact:None}}).unwrap();
    }
    MaterialCard::assemble(MaterialStateId {chemistry:"thermal fixture".into(),phase:"solid".into(),
        process:"synthetic".into(),revision:0},claims,vec![]).unwrap()
}
fn interface(card:&MaterialCard)->InterfaceSystemCard {
    crate::common::bind(card,card,Some(0.4),0.0).normal_state().elastic().interface().card().clone()
}
fn cards<'a>(card:&'a MaterialCard,interface:&'a InterfaceSystemCard)->PointContactMaterialCards<'a> {
    PointContactMaterialCards {surface_a:card,surface_b:card,interface,
        bulk_selection:[MaterialPropertySelection::SingleClaimOnly,MaterialPropertySelection::SingleClaimOnly],
        interface_selection:MaterialPropertySelection::SingleClaimOnly,
        model_selection:NormalContactModelSelection::SingleSupported,authority:InputAuthority::SyntheticFixture}
}
fn point()->QueryPoint {QueryPoint::new().with("T",300.0).unwrap()}
fn caloric(card:&MaterialCard,cp:f64)->EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(card.content_hash(),SolidLiquidPhase::Solid,vec![
        EnthalpyPhaseKnot {specific_enthalpy_j_kg:0.0,temperature_k:300.0,liquid_mass_fraction:0.0,bulk_density_kg_m3:1000.0},
        EnthalpyPhaseKnot {specific_enthalpy_j_kg:20.0*cp,temperature_k:320.0,liquid_mass_fraction:0.0,bulk_density_kg_m3:1000.0},
    ]).unwrap()
}
fn heated<'a>(cards:&PointContactMaterialCards<'_>,curve:&'a EquilibriumEnthalpyPhaseCurve,h:f64)->HeatedImpactSystem<'a> {
    let source=cards.resolve(&point(),PointContactGeometry::SpherePlane {effective_radius_m:0.1},envelope()).unwrap();
    let (a,wa)=ImpactBody::free_mass(0.001,-2e-5,0.04).unwrap();
    let (b,wb)=ImpactBody::free_mass(0.001,0.0,0.0).unwrap();
    let contact=MaterialContact::new("touch".into(),source,vec![wa,-wb],0.0).unwrap();
    let mut config=config();config.max_steps=8000;
    let mechanics=ImpactSystem::new_with_material_contacts(vec![a,b],vec![contact],vec![],vec![],config).unwrap();
    let stores=["left","right"].into_iter().map(|key|CaloricReservoir::new(key.into(),curve,
        curve.state_at_specific_enthalpy(h).unwrap(),0.001,1e-6,0.0006).unwrap()).collect();
    let mut route=route("touch");route.fraction_to_a=0.5;
    HeatedImpactSystem::new(mechanics.prepare_analytic().unwrap(),stores,vec![route],1e-11).unwrap()
}
fn resolve(c:&MaterialContact,s:&[CaloricReservoir<'_>],cards:&PointContactMaterialCards<'_>)
    ->Result<Option<SourcePointContact>,ImpactError> {
    if c.source().model().normal_state().elastic().state_coordinate("T")==Some(s[0].state().temperature_k()) {
        return Ok(None);
    }
    resolve_isothermal_contact(c,cards,&point(),[&s[0],&s[1]]).map(Some)
}
fn insulated(trial:ImpactThermalTrial<'_>)->Result<CoupledThermalStep<()>,ImpactError> {
    Ok(CoupledThermalStep {transfers:trial.reservoirs().iter().map(|r|adiabatic(r.thermal)).collect(),boundary_heat_j:0.0})
}

#[test]
fn heat_capacity_only_changes_the_next_strike_not_the_first_contact_episode() {
    let a=material(1e-6);let b=material(2e-6);let ia=interface(&a);let ib=interface(&b);
    let ca=caloric(&a,1e-6);let cb=caloric(&b,2e-6);let ma=cards(&a,&ia);let mb=cards(&b,&ib);
    let mut fast=heated(&ma,&ca,0.0);let mut slow=heated(&mb,&cb,0.0);let gate=CancelGate::new_clock_free();
    let initial=fast.mechanics().stored_energy_j();let mut work=0.0;
    let mut peaks=[[0.0_f64;2];2];let mut refreshed=0;let mut held=0;
    for sample in 0..6000 {
        let forces=if (2000..2500).contains(&sample){[1.0,-1.0]}else{[0.0;2]};
        let fa=fast.step_with_material_feedback(&forces,&gate,|c,s|resolve(c,s,&ma),insulated).unwrap();
        slow.step_with_material_feedback(&forces,&gate,|c,s|resolve(c,s,&mb),insulated).unwrap();
        refreshed+=fa.refreshed_contacts;held+=fa.held_loaded_contacts;
        work+=fa.thermal.impact.mechanical.supplied_work_j;
        let episode=usize::from(sample>=2000);
        peaks[0][episode]=peaks[0][episode].max(fast.mechanics().material_contact_observation(0).unwrap().normal_force_n);
        peaks[1][episode]=peaks[1][episode].max(slow.mechanics().material_contact_observation(0).unwrap().normal_force_n);
        if sample<2000 {assert_eq!(fast.mechanics().state(),slow.mechanics().state());}
        assert!(fa.thermal.impact.combined_energy_residual_j.abs()<1e-9);
    }
    assert_eq!(peaks[0][0],peaks[1][0]);assert!(peaks[0][1]>1.005*peaks[1][1] && peaks[1][1]>0.0);
    assert!(refreshed>=2 && held>0);
    assert!(fast.reservoirs()[0].state().temperature_k()>slow.reservoirs()[0].state().temperature_k());
    assert!(fast.mechanics().material_contacts()[0].source().stiffness_n_per_m_3_2()
        >slow.mechanics().material_contacts()[0].source().stiffness_n_per_m_3_2());
    near(fast.mechanics().stored_energy_j()+fast.reservoirs().iter().map(|s|s.storage_change_j()).sum::<f64>()-work,initial,1e-9);
}

#[test]
fn thermal_and_solver_refusals_restore_the_previous_source_as_well_as_motion_and_heat() {
    let card=material(0.02);let interface=interface(&card);let cards=cards(&card,&interface);let curve=caloric(&card,0.02);
    let mut s=heated(&cards,&curve,0.2);let gate=CancelGate::new_clock_free();
    let sources=s.mechanics().material_contacts().to_vec();let motion=s.mechanics().state().to_vec();let stores=s.reservoirs().to_vec();
    let result=s.step_with_material_feedback(&[0.0;2],&gate,|c,r|resolve(c,r,&cards),
        |_|Err::<CoupledThermalStep<()>,_>(ImpactError::Invalid("injected thermal refusal")));
    assert!(result.is_err());assert_eq!(s.mechanics().material_contacts(),sources);
    s.set_iteration_limit(0).unwrap();
    assert!(s.step_with_material_feedback(&[0.0;2],&gate,|c,r|resolve(c,r,&cards),insulated).is_err());
    s.set_iteration_limit(50).unwrap();let late=CancelGate::new_clock_free();
    assert!(s.step_with_material_feedback(&[0.0;2],&late,|c,r|resolve(c,r,&cards),|trial| {
        let result=insulated(trial);late.request();result
    }).is_err());
    assert_eq!(s.mechanics().material_contacts(),sources);assert_eq!(s.mechanics().state(),motion);
    assert_eq!(s.reservoirs(),stores);assert_eq!(s.mechanics().samples(),0);
    let frame=s.step_with_material_feedback(&[0.0;2],&gate,|c,r|resolve(c,r,&cards),insulated).unwrap();
    assert_eq!(frame.refreshed_contacts,1);assert_eq!(s.mechanics().samples(),1);
    assert_ne!(s.mechanics().material_contacts(),sources);
}

#[test]
fn thermal_source_binding_rejects_foreign_cards_unequal_temperatures_phase_and_context_changes() {
    let card=material(0.02);let interface=interface(&card);let cards=cards(&card,&interface);let curve=caloric(&card,0.02);
    let s=heated(&cards,&curve,0.2);let c=&s.mechanics().material_contacts()[0];let stores=s.reservoirs();
    let resolved=resolve_isothermal_contact(c,&cards,&point(),[&stores[0],&stores[1]]).unwrap();
    assert_eq!(resolved.model().normal_state().elastic().state_coordinate("T"),Some(310.0));
    let other=material(0.04);let other_interface=interface_from(&other);
    let other_cards=PointContactMaterialCards {surface_a:&other,surface_b:&other,interface:&other_interface,
        bulk_selection:[MaterialPropertySelection::SingleClaimOnly,MaterialPropertySelection::SingleClaimOnly],
        interface_selection:MaterialPropertySelection::SingleClaimOnly,
        model_selection:NormalContactModelSelection::SingleSupported,authority:InputAuthority::SyntheticFixture};
    assert!(resolve_isothermal_contact(c,&other_cards,&point(),[&stores[0],&stores[1]]).is_err());
    let cold=CaloricReservoir::new("cold".into(),&curve,curve.state_at_specific_enthalpy(0.0).unwrap(),0.001,1e-6,0.0006).unwrap();
    assert!(resolve_isothermal_contact(c,&cards,&point(),[&stores[0],&cold]).is_err());
    let extra=point().with("p",1e5).unwrap();
    assert!(resolve_isothermal_contact(c,&cards,&extra,[&stores[0],&stores[1]]).is_err());
    let latent=EquilibriumEnthalpyPhaseCurve::try_new(card.content_hash(),vec![
        EnthalpyPhaseKnot {specific_enthalpy_j_kg:0.0,temperature_k:310.0,liquid_mass_fraction:0.0,bulk_density_kg_m3:1000.0},
        EnthalpyPhaseKnot {specific_enthalpy_j_kg:1.0,temperature_k:310.0,liquid_mass_fraction:1.0,bulk_density_kg_m3:1000.0},
    ]).unwrap();
    let mixed=CaloricReservoir::new("mixed".into(),&latent,latent.state_at_specific_enthalpy(0.5).unwrap(),0.001,1e-6,0.0006).unwrap();
    assert!(resolve_isothermal_contact(c,&cards,&point(),[&mixed,&mixed]).is_err());
    assert!(cards.resolve(&point().with("T",330.0).unwrap(),c.source().geometry(),c.source().envelope()).is_err());
}
fn interface_from(card:&MaterialCard)->InterfaceSystemCard {interface(card)}

#[path="nonisothermal.rs"]
mod nonisothermal;
