//! G1/G3/G4: actual gas transport and finite sites, no material presets.
use super::*;
use super::super::tests::{fixture, controls as pair_controls};
use crate::thermochemical_gas::tests::point;
use crate::thermochemical_gas::energy::parcel::transfer::nozzle::{
    GasNozzle, network::snapshot::SnapshotNozzleLink,
};
use gas_network::SpeciesNetworkSolve;

fn controls() -> GasTransportSolve {
    GasTransportSolve { exchange: SpeciesNetworkSolve { maximum_species_fraction:0.1,
        mass_tolerance_kg:1e-13, energy_tolerance_j:2e-9 }, entropy_tolerance_j_k:1e-10 }
}
fn endpoints(gases:&BTreeMap<u64,ClosedGasParcel>)->BTreeMap<u64,GasEnergySolve> {
    gases.keys().map(|&id|(id,pair_controls().gas.gas)).collect()
}
fn link()->SnapshotNozzleLink {
    SnapshotNozzleLink { id:1, parcels:[0,1], nozzle:GasNozzle { area_m2:1e-8,discharge_coefficient:1.0 },
        nozzle_solve:pair_controls().gas.gas }
}
fn two_gases(catalog:&ThermochemicalGasModel)->BTreeMap<u64,ClosedGasParcel> {
    [(0,ClosedGasParcel::new(catalog.clone(),point(300.0,110_000.0),0.001).unwrap()),
        (1,ClosedGasParcel::new(catalog.clone(),point(300.0,100_000.0),0.001).unwrap())].into()
}
fn near(a:f64,b:f64,tol:f64) { assert!((a-b).abs()<=tol,"{a:.17e} != {b:.17e}"); }

#[test]
fn one_bound_surface_reuses_original_kinetics_and_finite_caloric_endpoint() {
    let (solid,gas,catalog)=fixture([0.0,0.0]);
    let pair=super::super::advance_candidate(&solid,&gas,&catalog,0.001,pair_controls(),||false).unwrap();
    let gases=BTreeMap::from([(0,gas.clone())]);
    let binding=AdsorptionNetworkBinding { id:7,parcel:0,solid:&solid,solve:pair_controls() };
    let candidate=step_candidate(&gases,GasTransportPaths::default(),&[binding],&catalog,0.001,
        &endpoints(&gases),controls(),||false).unwrap();
    assert_eq!(&candidate.solids()[&7],pair.solid());
    near(candidate.parcels()[&0].mass_kg(),pair.gas().mass_kg(),1e-13);
    near(candidate.parcels()[&0].state().state().temperature,pair.gas().state().state().temperature,2e-7);
    near(candidate.step().entropy_change_j_k,pair.step().entropy_change_j_k,1e-9);
    assert!(candidate.step().entropy_change_j_k>0.0);
    assert_eq!(candidate.step().exchange.boundaries.len(),2);
    near(candidate.step().total_energy_residual_j,0.0,controls().exchange.energy_tolerance_j);
    assert_eq!(gases[&0],gas);
}

#[test]
fn multiple_sites_use_one_old_bath_and_are_canonical_under_binding_order() {
    let (first,gas,catalog)=fixture([0.0,0.0]); let (second,_,_)=fixture([0.004,0.001]);
    let gases=BTreeMap::from([(0,gas)]);
    let bindings=[AdsorptionNetworkBinding { id:9,parcel:0,solid:&first,solve:pair_controls() },
        AdsorptionNetworkBinding { id:3,parcel:0,solid:&second,solve:pair_controls() }];
    let run=|bindings:&[AdsorptionNetworkBinding<'_>]|step_candidate(&gases,GasTransportPaths::default(),
        bindings,&catalog,0.001,&endpoints(&gases),controls(),||false).unwrap();
    let a=run(&bindings); let b=run(&[bindings[1],bindings[0]]); assert_eq!(a,b);
    for binding in &bindings {
        let isolated=super::super::advance_candidate(binding.solid,&gases[&0],&catalog,0.001,pair_controls(),||false).unwrap();
        assert_eq!(&a.solids()[&binding.id],isolated.solid());
    }
    assert_eq!(a.step().exchange.parcels.len(),1);
    assert_eq!(a.step().solids.iter().map(|r|r.id).collect::<Vec<_>>(),vec![3,9]);
    assert!(a.step().maximum_species_mass_residual_kg<=1e-13);
    assert!(a.step().entropy_change_j_k>=0.0);
}

#[test]
fn two_individually_admitted_surfaces_compete_for_one_original_gas_allowance() {
    let (first,gas,catalog)=fixture([0.0,0.0]); let (second,_,_)=fixture([0.0,0.0]);
    let gases=BTreeMap::from([(0,gas)]);
    let bindings=[AdsorptionNetworkBinding { id:0,parcel:0,solid:&first,solve:pair_controls() },
        AdsorptionNetworkBinding { id:1,parcel:0,solid:&second,solve:pair_controls() }];
    let mut c=controls(); c.exchange.maximum_species_fraction=0.006;
    for binding in bindings { assert!(step_candidate(&gases,GasTransportPaths::default(),&[binding],
        &catalog,0.001,&endpoints(&gases),c,||false).is_ok()); }
    assert!(matches!(step_candidate(&gases,GasTransportPaths::default(),&bindings,&catalog,0.001,
        &endpoints(&gases),c,||false),Err(AdsorptionNetworkError::Gas(SpeciesNetworkError::Withdrawal { .. }))));
}

#[test]
fn nozzle_and_adsorption_sample_once_and_all_final_species_energy_close() {
    let (solid,_,catalog)=fixture([0.04,0.0]); let gases=two_gases(&catalog); let saved=(gases.clone(),solid.clone());
    let nozzles=[link()]; let paths=GasTransportPaths { nozzles:&nozzles,..GasTransportPaths::default() };
    let binding=AdsorptionNetworkBinding { id:42,parcel:0,solid:&solid,solve:pair_controls() };
    let coupled=step_candidate(&gases,paths,&[binding],&catalog,0.001,&endpoints(&gases),controls(),||false).unwrap();
    let transport_only=transport::step_candidate_with_paths(&gases,paths,0.001,&endpoints(&gases),controls(),||false).unwrap();
    // Desorbed A does not alter this interval's originally sampled nozzle flow.
    assert_eq!(coupled.step().transport.nozzles,transport_only.step().nozzles);
    assert_eq!(coupled.step().exchange.fluxes,transport_only.step().exchange.fluxes);
    assert!(coupled.solids()[&42].moles[0]<solid.moles[0]);
    for (i,law) in solid.model.channels.iter().enumerate() {
        let gas_change:f64=gases.iter().map(|(&id,old)|coupled.parcels()[&id].species_mass_kg(&law.species).unwrap()
            -old.species_mass_kg(&law.species).unwrap()).sum();
        near(gas_change+(coupled.solids()[&42].moles[i]-solid.moles[i])*law.molar_mass(),0.0,1e-13);
    }
    let gas_energy:f64=coupled.step().exchange.parcels.iter().map(|r|r.energy_change_j).sum();
    let solid_energy:f64=coupled.step().solids.iter().map(|r|r.solid_energy_change_j).sum();
    near(gas_energy+solid_energy,0.0,2e-8);
    assert!(coupled.step().entropy_change_j_k>=0.0); assert_eq!((gases,solid),saved);
}

#[test]
fn empty_site_list_keeps_original_closed_gas_transport_endpoints() {
    let (_,_,catalog)=fixture([0.0,0.0]); let gases=two_gases(&catalog); let nozzles=[link()];
    let paths=GasTransportPaths { nozzles:&nozzles,..GasTransportPaths::default() };
    let old=transport::step_candidate_with_paths(&gases,paths,0.001,&endpoints(&gases),controls(),||false).unwrap();
    let new=step_candidate(&gases,paths,&[],&catalog,0.001,&endpoints(&gases),controls(),||false).unwrap();
    assert_eq!(old.parcels(),new.parcels()); assert_eq!(old.step().exchange,new.step().exchange);
    near(old.step().entropy_change_j_k,new.step().entropy_change_j_k,1e-10);
    assert!(new.solids().is_empty());
}

#[test]
fn late_joint_thermal_failure_aliasing_and_cancellation_publish_nothing() {
    let (first,gas,catalog)=fixture([0.0,0.0]); let (second,_,_)=fixture([0.0,0.0]);
    let gases=BTreeMap::from([(0,gas)]); let saved=(gases.clone(),first.clone(),second.clone());
    let mut local=pair_controls(); local.maximum_temperature_change_k=0.5;
    let a=AdsorptionNetworkBinding { id:0,parcel:0,solid:&first,solve:local };
    let b=AdsorptionNetworkBinding { id:1,parcel:0,solid:&second,solve:local };
    let result=step_candidate(&gases,GasTransportPaths::default(),&[a,b],&catalog,0.001,
        &endpoints(&gases),controls(),||false);
    assert!(matches!(result,Err(AdsorptionNetworkError::Adsorption {
        error:AdsorptionError::Invalid("joint gas endpoint exceeds the adsorption thermal regime"),.. })));
    let alias=AdsorptionNetworkBinding { id:2,..a };
    assert!(step_candidate(&gases,GasTransportPaths::default(),&[a,alias],&catalog,0.001,
        &endpoints(&gases),controls(),||false).is_err());
    let mut calls=0;
    assert!(matches!(step_candidate(&gases,GasTransportPaths::default(),&[a],&catalog,0.001,
        &endpoints(&gases),controls(),||{ calls+=1; calls>15 }),Err(AdsorptionNetworkError::Cancelled)));
    assert_eq!((gases,first,second),saved);
}

#[test]
fn continuing_network_uses_the_new_bath_without_resetting_site_history() {
    let (solid,gas,catalog)=fixture([0.0,0.0]); let mut gases=BTreeMap::from([(0,gas)]); let mut current=solid.clone();
    for _ in 0..4 {
        let binding=AdsorptionNetworkBinding { id:8,parcel:0,solid:&current,solve:pair_controls() };
        let candidate=step_candidate(&gases,GasTransportPaths::default(),&[binding],&catalog,0.00025,
            &endpoints(&gases),controls(),||false).unwrap();
        assert!(candidate.solids()[&8].moles[0]>current.moles[0]);
        let (new_gas,new_solids,_)=candidate.into_parts(); gases=new_gas; current=new_solids[&8].clone();
    }
    assert!(current.temperature_k>solid.temperature_k && gases[&0].state().state().temperature<300.0);
}

#[test]
fn source_bound_dusty_gas_permeation_and_uptake_use_the_same_original_pore_state() {
    use crate::porous::permeability::{DarcyGasBarrier,IntrinsicPermeability,PorousGasLayer,PorousGasSupport};
    use crate::thermochemical_gas::energy::parcel::transfer::species::diffusion::dusty_gas::{
        DustyGasFaceSolve,DustyGasMedium,DustyGasPores,
    };
    use transport::dusty_gas::DustyGasNetworkLink;
    use fs_thermochem::{Composition,CompositionBasis};
    let (fixture,_,catalog)=fixture([0.0,0.0]);
    let solid=CompetitiveSolid::ordered(Arc::new(CompetitiveAdsorbent::try_new(
        vec![fixture.model.channels[0].clone()]).unwrap()),0.1,vec![0.0],300.0).unwrap();
    let pure=catalog.with_composition(Composition::new(CompositionBasis::MoleFraction,vec![1.0,0.0]).unwrap()).unwrap();
    let gases=two_gases(&pure); let context=point(300.0,100_000.0);
    let layer=PorousGasLayer::try_new(IntrinsicPermeability::isotropic(2e-14).unwrap(),0.01,0.0,
        "synthetic adsorbing pore permeability").unwrap();
    let barrier=DarcyGasBarrier::try_new(vec![layer],[1.0,0.0,0.0],PorousGasSupport {
        temperature_k:(200.0,800.0),pressure_pa:(1000.0,1e7),maximum_superficial_speed_m_s:100.0 }).unwrap();
    let medium=DustyGasMedium::try_new(barrier,DustyGasPores::try_new(1e-7,0.4,2.0,
        "synthetic dry pore geometry").unwrap()).unwrap();
    let links=[DustyGasNetworkLink { id:5,parcels:[0,1],medium:&medium,area_m2:0.001,pairs:&[],
        gas_context:&context,solve:DustyGasFaceSolve { maximum_temperature_difference_k:1.0,
            maximum_relative_pressure_difference:0.25,gradient_tolerance_per_m:1e-8 } }];
    let paths=GasTransportPaths { dusty_gas:&links,..GasTransportPaths::default() };
    let binding=AdsorptionNetworkBinding { id:6,parcel:0,solid:&solid,solve:pair_controls() };
    let result=step_candidate(&gases,paths,&[binding],&catalog,1e-4,&endpoints(&gases),controls(),||false).unwrap();
    let only=transport::step_candidate_with_paths(&gases,paths,1e-4,&endpoints(&gases),controls(),||false).unwrap();
    assert_eq!(result.step().transport.dusty_gas,only.step().dusty_gas);
    assert_eq!(result.step().exchange.fluxes,only.step().exchange.fluxes);
    assert!(result.parcels()[&1].mass_kg()>gases[&1].mass_kg());
    assert!(result.solids()[&6].moles[0]>0.0);
    near(result.step().maximum_species_mass_residual_kg,0.0,1e-13);
    near(result.step().total_energy_residual_j,0.0,2e-9);
}
