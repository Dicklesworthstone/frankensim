//! Manufactured references; no empirical adsorbent validation is asserted.
use super::*;
use crate::phase::{HeatCapacitySourceQuery, SolidLiquidPhase};
use crate::thermochemical_gas::{energy::GasEnergySolve,tests::{model,point}};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet,InterpolationPolicy,MaterialStateId,PropertyClaim,PropertyKey,
    PropertyValue,Provenance,SelectionPolicy,UncertaintyModel};

fn claim(name:&str,dims:Dims,value:f64)->PropertyClaim {
    PropertyClaim { key:PropertyKey::new(name,dims),value:PropertyValue::Scalar{value,dims},
        validity:ValidityDomain::unconstrained().with_quantity("temperature",temp_quantity(),200.0,800.0)
            .with_quantity("pressure",QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure,ValueForm::Static)),1000.0,1e7),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
        observations:vec![],provenance:Provenance{source:"manufactured competitive adsorption".into(),license:"CC0-1.0".into(),artifact:None} }
}
fn card(name:&str,claims:Vec<PropertyClaim>)->MaterialCard {
    let mut set=ClaimSet::new(); for c in claims { set.insert_claim(c).unwrap(); }
    MaterialCard::assemble(MaterialStateId{chemistry:name.into(),phase:"solid".into(),process:"analytic competitive fixture".into(),revision:0},set,vec![]).unwrap()
}
fn chart(name:&str,cp:f64)->ResolvedHeatCapacityCurve {
    let c=Dims([2,0,-2,-1,0,0]); let r=Dims([-3,1,0,0,0,0]);
    let source=card(name,vec![claim("specific_heat_capacity",c,cp),claim("density",r,1000.0)]);
    ResolvedHeatCapacityCurve::try_new(&source,HeatCapacitySourceQuery{
        heat_capacity:&PropertyKey::new("specific_heat_capacity",c),density:&PropertyKey::new("density",r),
        temperature_axis:"temperature",lower:&point(200.0,1e5),upper:&point(800.0,1e5),selection:SelectionPolicy::SingleClaimOnly,
    },SolidLiquidPhase::Solid,0.0,1e-6).unwrap()
}
pub(super) fn fixture(loadings:[f64;2])->(CompetitiveSolid,ClosedGasParcel,ThermochemicalGasModel) {
    let catalog=model(vec![0.5,0.5],0.0); let dry=chart("shared dry skeleton",500.0);
    let mut laws=Vec::new();
    for (name,mass,cp,affinity,heat) in [("A",0.028,3.5,2e-5,6000.0),("B",0.004,2.5,8e-5,9000.0)] {
        let species=SpeciesId::new(name).unwrap(); let bound=chart(name,cp*R/mass);
        let hg=catalog.species_specific_enthalpy_j_kg(&species,300.0).unwrap();
        let hb=bound.state_at_temperature(300.0).unwrap().specific_enthalpy_j_kg();
        let source=AdsorptionSource { material:card(name,vec![
            claim(LANGMUIR_CAPACITY_PROPERTY,LANGMUIR_CAPACITY_DIMS,0.5),
            claim(LANGMUIR_REFERENCE_AFFINITY_PROPERTY,LANGMUIR_AFFINITY_DIMS,affinity),
            claim(LANGMUIR_DESORPTION_RATE_PROPERTY,LANGMUIR_RATE_DIMS,0.4),
            claim("adsorbed_surface_diffusivity",Dims([2,0,-1,0,0,0]),if name=="A"{1e-8}else{2e-8}),
        ]),point:point(300.0,1e5),selection:MaterialPropertySelection::SingleClaimOnly,reference_temperature_k:300.0 };
        laws.push(LangmuirAdsorbent::try_new(&catalog,species,source,AdsorptionCaloricBinding {
            dry_skeleton:dry.clone(),bound_adsorbate:bound,bound_enthalpy_offset_j_kg:hg-hb-heat/mass,
        }).unwrap());
    }
    let law=Arc::new(CompetitiveAdsorbent::try_new(laws).unwrap());
    let solid=CompetitiveSolid::try_new(law,0.1,&[(SpeciesId::new("A").unwrap(),loadings[0]),
        (SpeciesId::new("B").unwrap(),loadings[1])],300.0).unwrap();
    let gas=ClosedGasParcel::new(catalog.clone(),point(300.0,1e5),0.001).unwrap();
    (solid,gas,catalog)
}
pub(super) fn controls()->CompetitiveSolve {
    CompetitiveSolve { gas:GasSpeciesBoundarySolve { gas:GasEnergySolve {
        minimum_temperature_k:200.0,maximum_temperature_k:800.0,energy_tolerance_j_kg:1e-7,maximum_iterations:128,
    },mass_tolerance_kg:1e-13,energy_tolerance_j:1e-8 },maximum_gas_gross_fraction:0.1,maximum_site_turnover:0.25,
        maximum_gas_surface_temperature_difference_k:5.0,maximum_temperature_change_k:2.0,
        entropy_tolerance_j_k:1e-10,maximum_iterations:128 }
}
fn near(a:f64,b:f64,tol:f64) { assert!((a-b).abs()<=tol,"{a:.17e} != {b:.17e}"); }
#[test]
fn both_species_heat_one_skeleton_and_share_one_gas_reconstruction() {
    let (solid,gas,catalog)=fixture([0.0,0.0]); let saved=(solid.clone(),gas.clone());
    let result=advance_candidate(&solid,&gas,&catalog,0.001,controls(),||false).unwrap();
    assert_eq!((solid.clone(),gas.clone()),saved);
    assert!(result.solid().temperature_k()>300.0 && result.gas().state().state().temperature<300.0);
    // Equal desorption rates make the common vacancy C/(1+dt*sum(a_i)).
    // BE retains d_i as well: n_i=C*dt*a_i/(1+dt*(d+sum(a_i))).
    let c=0.05; let dt=0.001; let a=[0.4,1.6];
    for i in 0..2 {
        let expected=c*dt*a[i]/(1.0+dt*(0.4+a[0]+a[1]));
        near(result.solid().moles[i],expected,2e-16);
        let law=&solid.model.channels[i];
        near(gas.species_mass_kg(&law.species).unwrap()-result.gas().species_mass_kg(&law.species).unwrap(),
            expected*law.molar_mass(),1e-13);
    }
    let n=&result.solid().moles;
    let expected_heat=6000.0*n[0]+9000.0*n[1];
    let expected_capacity=0.1*500.0+R*(3.5*n[0]+2.5*n[1]);
    near(result.solid().temperature_k(),300.0+expected_heat/expected_capacity,2e-10);
    near(result.step().solid_energy_change_j+result.step().gas_energy_change_j,0.0,2e-8);
    assert!(result.step().entropy_change_j_k>=-1e-10);
    assert!(result.step().thermal_iterations<=3);
}
#[test]
fn common_equilibrium_retains_both_gross_streams_and_no_temperature_change() {
    let (solid,gas,catalog)=fixture([0.05/6.0,0.05*4.0/6.0]);
    let result=advance_candidate(&solid,&gas,&catalog,0.001,controls(),||false).unwrap();
    for (a,b) in solid.moles.iter().zip(&result.solid().moles) { near(*a,*b,2e-16); }
    near(result.solid().temperature_k(),300.0,1e-10);
    for row in &result.step().species { assert!(row.equilibrium_coverage[0]>0.0); }
    assert!(result.step().sites.gross_adsorption_mol.iter().all(|&n|n>0.0));
    near(result.step().entropy_change_j_k,0.0,1e-10);
}
#[test]
fn weak_species_is_displaced_and_all_source_orderings_are_canonical() {
    let (solid,gas,catalog)=fixture([0.04,0.0]);
    let a=advance_candidate(&solid,&gas,&catalog,0.001,controls(),||false).unwrap();
    assert!(a.solid().moles[0]<solid.moles[0] && a.solid().moles[1]>0.0);
    let mut laws=solid.model.channels.clone(); laws.reverse();
    let reversed=CompetitiveSolid::try_new(Arc::new(CompetitiveAdsorbent::try_new(laws).unwrap()),0.1,
        &[(SpeciesId::new("B").unwrap(),0.0),(SpeciesId::new("A").unwrap(),0.04)],300.0).unwrap();
    let b=advance_candidate(&reversed,&gas,&catalog,0.001,controls(),||false).unwrap();
    assert_eq!(a,b);
}
#[test]
fn common_capacity_source_identity_and_atomic_refusals() {
    let (solid,gas,catalog)=fixture([0.0,0.0]); let saved=(solid.clone(),gas.clone());
    let mut laws=solid.model.channels.clone(); laws[1].capacity_mol_kg_dry*=2.0;
    assert!(CompetitiveAdsorbent::try_new(laws).is_err());
    let duplicate=vec![solid.model.channels[0].clone();2]; assert!(CompetitiveAdsorbent::try_new(duplicate).is_err());
    assert!(CompetitiveSolid::try_new(Arc::clone(&solid.model),0.1,&[(SpeciesId::new("A").unwrap(),0.04),
        (SpeciesId::new("B").unwrap(),0.04)],300.0).is_err());
    let mut c=controls(); c.maximum_site_turnover=1e-10;
    assert!(matches!(advance_candidate(&solid,&gas,&catalog,0.1,c,||false),Err(AdsorptionError::InventoryLimit{..})));
    let changed=model(vec![0.5,0.5],0.001);
    assert!(advance_candidate(&solid,&gas,&changed,0.001,controls(),||false).is_err());
    let mut polls=0; assert!(matches!(advance_candidate(&solid,&gas,&catalog,0.001,controls(),||{
        polls+=1; polls>10 }),Err(AdsorptionError::Cancelled)));
    assert_eq!((solid,gas),saved);
}
#[test]
fn zero_duration_is_exact_and_advancing_temperature_is_checked() {
    let (solid,gas,catalog)=fixture([0.0,0.0]);
    let result=advance_candidate(&solid,&gas,&catalog,0.0,controls(),||false).unwrap();
    assert_eq!(result.solid(),&solid); assert_eq!(result.gas(),&gas);
    assert_eq!(result.step().entropy_change_j_k,0.0);
    let mut c=controls(); c.maximum_temperature_change_k=0.0;
    assert!(advance_candidate(&solid,&gas,&catalog,0.001,c,||false).is_err());
}
