//! Synthetic independent-permeant fixtures exercise the original gas owners.
use super::*;
use super::super::mixture::*;
use super::super::reservoir::{IsothermalGasMembraneSolve,IsothermalGasMembraneError as MixtureError};
use crate::thermochemical_gas::{ThermochemicalGasModel,tests::{species_with_energy_reference,
    transport_in_domain,point as gas_point}};
use crate::thermochemical_gas::energy::{GasEnergySolve,parcel::ClosedGasParcel};
use crate::thermochemical_gas::energy::parcel::transfer::species::boundary::GasSpeciesBoundarySolve;
use fs_thermochem::{Composition,CompositionBasis,mixture::FrozenIdealGasMixtureModelV1,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use fs_qty::semantic::{CompositionBasis as SemanticCompositionBasis,MoistureBasis};

fn id(name:&str)->SpeciesId {SpeciesId::new(name).unwrap()}
fn solve()->IsothermalGasMembraneSolve {
    IsothermalGasMembraneSolve {membrane:controls(),gas:[GasSpeciesBoundarySolve {
        gas:GasEnergySolve {minimum_temperature_k:200.0,maximum_temperature_k:1000.0,
            energy_tolerance_j_kg:1e-7,maximum_iterations:128},
        mass_tolerance_kg:1e-16,energy_tolerance_j:1e-8};2],mass_tolerance_kg:1e-16}
}
fn gas_model(fractions:[f64;2],reference:f64)->ThermochemicalGasModel {
    ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(vec![
        species_with_energy_reference("A",0.028,3.5,0.001,reference),
        species_with_energy_reference("B",0.044,4.0,0.002,reference)],
        Composition::new(CompositionBasis::MoleFraction,fractions.to_vec()).unwrap()).unwrap(),
        vec![transport_in_domain("A",1.7e-5,0.025,[200.0,1000.0],[100.0,1e7]),
            transport_in_domain("B",2.2e-5,0.04,[200.0,1000.0],[100.0,1e7])],None).unwrap()
}
fn gas(model:ThermochemicalGasModel,p:f64,volume:f64)->ClosedGasParcel {
    ClosedGasParcel::new(model,gas_point(300.0,p).with("specimen-lot",17.0).unwrap(),volume).unwrap()
}
fn card(name:&str,d:f64,s:f64,bounds:[f64;2])->MaterialCard {
    let mut claims=ClaimSet::new();
    claims.insert_claim(claim(MEMBRANE_DIFFUSIVITY_PROPERTY,d,MEMBRANE_DIFFUSIVITY_DIMS,bounds)).unwrap();
    claims.insert_claim(claim(MOLAR_HENRY_SOLUBILITY_PROPERTY,s,MOLAR_HENRY_SOLUBILITY_DIMS,bounds)).unwrap();
    MaterialCard::assemble(MaterialStateId {chemistry:name.into(),phase:"dense solid".into(),
        process:"independent mixture analytic fixture".into(),revision:0},claims,vec![]).unwrap()
}
fn channel<'a>(name:&str,card:&'a MaterialCard,q:&'a QueryPoint,loading:f64)->MembraneMixtureChannelSource<'a> {
    let mut layer=layer(card,q,8);layer.material.species=id(name);
    MembraneMixtureChannelSource {layers:vec![layer],interfaces:vec![],initial_partial_pressures_pa:vec![loading;8]}
}
fn near(a:f64,b:f64,tol:f64) {assert!((a-b).abs()<=tol,"{a} != {b}");}
fn audit(state:&IsothermalGasMembraneMixture<'_>,step:&IsothermalGasMembraneMixtureStep) {
    assert!(step.mass_residual_kg.abs()<=solve().mass_tolerance_kg);
    assert!(step.total_mass_residual_kg.abs()<=solve().mass_tolerance_kg);
    for row in &step.species {
        assert!(row.mass_residual_kg.abs()<=solve().mass_tolerance_kg);
        assert!(row.total_mass_residual_kg.abs()<=solve().mass_tolerance_kg);
        for side in 0..2 {
            assert_eq!(row.actual_partial_pressures_pa[side],
                state.gases()[side].species_partial_pressure_pa(&row.species).unwrap());
            assert_eq!(state.finite_state(&row.species).unwrap().reservoir_pressures_pa()[side],
                row.actual_partial_pressures_pa[side]);
        }
    }
    for gas in state.gases() {
        assert_eq!(gas.state().state().temperature,300.0);
        assert_eq!(&gas.model().evaluate(gas.point()).unwrap(),gas.state());
        near(gas.state().state().density*gas.volume_m3(),gas.mass_kg(),1e-16);
        assert_eq!(gas.point().axes().get("specimen-lot"),Some(&17.0));
    }
    assert_eq!(step.gas_energy_residual_j,0.0);
}

#[test]
fn mixtures_separate_against_total_pressure_gradient_with_original_gas_reconstruction() {
    let (a,b)=(card("fast A",1e-9,2e-5,[90000.0,210000.0]),
        card("slow B",1e-10,1e-5,[90000.0,210000.0]));
    let q=point(150000.0);
    let sources=[channel("A",&a,&q,60000.0),channel("B",&b,&q,90000.0)];
    let gases=[gas(gas_model([0.8,0.2],0.0),1e5,1e-6),
        gas(gas_model([0.2,0.8],0.0),2e5,2e-6)];
    let state=IsothermalGasMembraneMixture::new(300.0,0.01,&sources,gases.clone(),solve(),||false).unwrap();
    // A partial pressures (80/40 kPa) are below this claim's 90 kPa total-
    // pressure floor. Correct mixture admission uses totals (100/200 kPa).
    let source=state.source_support(&id("A")).unwrap();
    assert_eq!(source.layers[0][0].state_coordinate("pressure"),Some(100000.0));
    assert_eq!(source.layers[0][1].state_coordinate("pressure"),Some(200000.0));
    let mut work=state.workspace().unwrap();let next=state.step_candidate(0.1,&mut work,||false).unwrap();
    let report=next.step();audit(next.state(),report);
    assert!(report.species[0].actual_gas_mass_changes_kg[0]<0.0);
    assert!(report.species[0].actual_gas_mass_changes_kg[1]>0.0);
    assert!(report.species[1].actual_gas_mass_changes_kg[0]>0.0);
    assert!(report.species[1].actual_gas_mass_changes_kg[1]<0.0);
    let before=gases[1].species_partial_pressure_pa(&id("A")).unwrap()/gases[1].state().state().pressure;
    let after=next.state().gases()[1].species_partial_pressure_pa(&id("A")).unwrap()
        /next.state().gases()[1].state().state().pressure;
    assert!(after>before);
    for side in 0..2 {
        assert_eq!(report.gas[side].species.len(),2);
        assert_eq!(next.state().gases()[side].model().thermodynamics().components(),
            gases[side].model().thermodynamics().components());
        assert_eq!(next.state().gases()[side].model().transport_materials(),gases[side].model().transport_materials());
    }
    assert_eq!(state.gases(),&gases);
}

#[test]
fn blocked_species_keeps_gas_and_trapped_storage_and_reference_and_order_do_not_change_solution() {
    let (a,b)=(card("A",1e-10,2e-5,[0.0,3e5]),card("blocked B",0.0,1e-5,[0.0,3e5]));
    let q=point(100000.0);let mut expected=None;
    for (reference,reversed) in [(0.0,false),(1e24,false),(1e24,true)] {
        let mut sources=vec![channel("A",&a,&q,50000.0),channel("B",&b,&q,12345.0)];
        if reversed {sources.reverse();}
        let gases=[gas(gas_model([0.8,0.2],reference),1e5,1e-7),
            gas(gas_model([0.2,0.8],reference),1e5,3e-7)];
        let saved=[gases[0].species_mass_kg(&id("B")).unwrap(),gases[1].species_mass_kg(&id("B")).unwrap()];
        let state=IsothermalGasMembraneMixture::new(300.0,0.01,&sources,gases,solve(),||false).unwrap();
        let trapped=state.finite_state(&id("B")).unwrap().membrane().stored_amount_mol().unwrap();
        let mut work=state.workspace().unwrap();let next=state.step_candidate(1.0,&mut work,||false).unwrap();
        audit(next.state(),next.step());
        let blocked=&next.step().species[1];
        assert_eq!(blocked.species,id("B"));assert_eq!(blocked.sorbed_mass_change_kg,0.0);
        assert_eq!(blocked.diffusion.membrane.face_molar_fluxes_mol_m2_s,vec![0.0;9]);
        assert_eq!(next.state().finite_state(&id("B")).unwrap().membrane().stored_amount_mol().unwrap(),trapped);
        assert_eq!(next.state().finite_state(&id("B")).unwrap().membrane().partial_pressures_pa(),&[12345.0;8]);
        for side in 0..2 {near(next.state().gases()[side].species_mass_kg(&id("B")).unwrap(),saved[side],1e-20);}
        let values=(next.state().gases().each_ref().map(|g|g.mass_kg()),
            next.state().gases().each_ref().map(|g|g.state().state().pressure),
            next.state().finite_state(&id("A")).unwrap().membrane().partial_pressures_pa().to_vec(),
            next.step().thermostat_heat_into_gas_j);
        if let Some(expected)=&expected {assert_eq!(&values,expected);}else{expected=Some(values);}
        if reference>0.0 {assert!(next.step().stream_enthalpy_into_gas_j[0].abs()>1e12);}
    }
}

#[test]
fn each_species_equilibrates_using_gas_plus_its_own_sorbed_capacity() {
    let (a,b)=(card("A",1e-9,2e-5,[0.0,4e5]),card("B",2e-10,5e-5,[0.0,4e5]));
    let q=point(100000.0);let sources=[channel("A",&a,&q,0.0),channel("B",&b,&q,25000.0)];
    let gases=[gas(gas_model([0.8,0.2],0.0),1e5,1e-7),gas(gas_model([0.1,0.9],0.0),2e5,3e-7)];
    let mut state=IsothermalGasMembraneMixture::new(300.0,0.01,&sources,gases,solve(),||false).unwrap();
    let initial=state.total_mass_kg().unwrap();let mut expected=[0.0;2];
    for (i,name) in ["A","B"].iter().enumerate() {
        let total=state.finite_state(&id(name)).unwrap().total_amount_mol().unwrap();
        expected[i]=total/((1e-7+3e-7)/(R*300.0)+0.01*1e-4*[2e-5,5e-5][i]);
    }
    let mut work=state.workspace().unwrap();let mut converged=false;
    for _ in 0..80 {
        let next=state.step_candidate(100.0,&mut work,||false).unwrap();audit(next.state(),next.step());
        state=next.into_parts().0;near(state.total_mass_kg().unwrap(),initial,2e-15);
        converged=["A","B"].iter().enumerate().all(|(i,name)| {
            let finite=state.finite_state(&id(name)).unwrap();
            finite.reservoir_pressures_pa().iter().copied()
                .chain(finite.membrane().partial_pressures_pa().iter().copied())
                .all(|p|(p-expected[i]).abs()<1e-5*expected[i])
        });
        if converged {break;}
    }
    assert!(converged);
}

#[test]
fn absent_endpoint_species_enters_from_the_retained_catalog() {
    let (a,b)=(card("A",1e-9,2e-5,[0.0,3e5]),card("B",1e-10,1e-5,[0.0,3e5]));
    let q=point(100000.0);let sources=[channel("A",&a,&q,0.0),channel("B",&b,&q,50000.0)];
    let bank=gas_model([0.5,0.5],0.0);
    let pure_b=bank.with_mole_amounts(&[fs_qty::Amount::new(0.0),fs_qty::Amount::new(1.0)]).unwrap();
    let gases=[gas(bank,1e5,1e-6),gas(pure_b,1e5,2e-6)];
    assert_eq!(gases[1].species_mass_kg(&id("A")).unwrap(),0.0);
    let state=IsothermalGasMembraneMixture::new(300.0,0.01,&sources,gases,solve(),||false).unwrap();
    let mut work=state.workspace().unwrap();let next=state.step_candidate(10.0,&mut work,||false).unwrap();
    audit(next.state(),next.step());
    assert!(next.state().gases()[1].species_mass_kg(&id("A")).unwrap()>0.0);
    assert_eq!(next.state().catalog(),state.catalog());
    assert_eq!(next.state().gases()[1].model().thermodynamics().components(),state.catalog().thermodynamics().components());
}

#[test]
fn total_pressure_dynamic_conditions_channel_coverage_and_cancel_are_real_refusals() {
    let (a,b)=(card("A",1e-10,2e-5,[0.0,3e5]),card("B",1e-10,1e-5,[0.0,3e5]));
    let q=point(100000.0);let sources=[channel("A",&a,&q,50000.0),channel("B",&b,&q,50000.0)];
    let gases=[gas(gas_model([0.5,0.5],0.0),1e5,1e-7),gas(gas_model([0.5,0.5],0.0),2e5,3e-7)];
    assert!(IsothermalGasMembraneMixture::new(300.0,0.01,&sources[..1],gases.clone(),solve(),||false).is_err());
    let mut wrong=sources.clone();wrong[1].layers[0].material.species=id("A");
    assert!(IsothermalGasMembraneMixture::new(300.0,0.01,&wrong,gases.clone(),solve(),||false).is_err());
    let mut wrong=sources.clone();wrong[1].layers[0].thickness_m*=2.0;
    assert!(IsothermalGasMembraneMixture::new(300.0,0.01,&wrong,gases.clone(),solve(),||false).is_err());
    let narrow=card("narrow total",1e-10,2e-5,[0.0,150000.0]);
    let wrong=[channel("A",&narrow,&q,50000.0),sources[1].clone()];
    // Each partial pressure is <=100 kPa, but the actual total is 200 kPa.
    assert!(IsothermalGasMembraneMixture::new(300.0,0.01,&wrong,gases.clone(),solve(),||false).is_err());
    let equal_gases=[gas(gas_model([0.5,0.5],0.0),1e5,1e-7),gas(gas_model([0.5,0.5],0.0),1e5,3e-7)];
    let interior=[channel("A",&narrow,&q,100000.0),channel("B",&b,&q,100000.0)];
    assert!(IsothermalGasMembraneMixture::new(300.0,0.01,&interior,equal_gases.clone(),solve(),||false).is_err());
    // Both cells and gas endpoints total 100 kPa. Opposing species profiles
    // with complementary layer resistances make the shared contact ~180 kPa.
    let fast=card("contact fast",1e-8,1e-5,[0.0,150000.0]);
    let slow=card("contact slow",1e-12,1e-5,[0.0,150000.0]);
    let mut contact=Vec::new();
    for (name,left,right,loading) in [("A",&fast,&slow,vec![90000.0,10000.0]),
        ("B",&slow,&fast,vec![10000.0,90000.0])] {
        let mut l=layer(left,&q,1);l.material.species=id(name);
        let mut r=layer(right,&q,1);r.material.species=id(name);
        contact.push(MembraneMixtureChannelSource {layers:vec![l,r],interfaces:vec![None],initial_partial_pressures_pa:loading});
    }
    assert!(IsothermalGasMembraneMixture::new(300.0,0.01,&contact,equal_gases,solve(),||false).is_err());
    for kind in [QuantityKind::Composition(SemanticCompositionBasis::MoleFraction),
        QuantityKind::Moisture(MoistureBasis::RelativeHumidity)] {
        let mut claims=ClaimSet::new();
        for key in [MEMBRANE_DIFFUSIVITY_PROPERTY,MOLAR_HENRY_SOLUBILITY_PROPERTY] {
            let mut property=a.claims_for(key)[0].1.clone();
            property.validity=property.validity.with_quantity("mixture-condition",typed(kind),0.0,1.0);
            claims.insert_claim(property).unwrap();
        }
        let restricted=MaterialCard::assemble(a.id().clone(),claims,vec![]).unwrap();
        let context=q.clone().with_quantity("mixture-condition",typed(kind),0.5).unwrap();
        let wrong=[channel("A",&restricted,&context,50000.0),sources[1].clone()];
        assert!(IsothermalGasMembraneMixture::new(300.0,0.01,&wrong,gases.clone(),solve(),||false).is_err());
    }
    let state=IsothermalGasMembraneMixture::new(300.0,0.01,&sources,gases.clone(),solve(),||false).unwrap();
    let saved=state.finite_state(&id("A")).unwrap().clone();let mut work=state.workspace().unwrap();let mut polls=0;
    assert!(matches!(state.step_candidate(1.0,&mut work,||{polls+=1;polls==30}),Err(MixtureError::Cancelled)));
    assert_eq!(state.gases(),&gases);assert_eq!(state.finite_state(&id("A")).unwrap(),&saved);
}
