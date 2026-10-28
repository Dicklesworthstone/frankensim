//! G1/G3/G4: the actual gas owner, manufactured open-phase streams.
use super::*;
use crate::thermochemical_gas::tests::{model, point, species_with_energy_reference, transport};
use fs_thermochem::{Composition, CompositionBasis, mixture::FrozenIdealGasMixtureModelV1};

const M: f64 = 0.028;
const CV: f64 = 2.5 * R / M;
const CP: f64 = 3.5 * R / M;
fn gases(reference: f64) -> BTreeMap<u64, ClosedGasParcel> {
    let model = ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(
        vec![species_with_energy_reference("A", M, 3.5, 0.0, reference)],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap()).unwrap(),
        vec![transport("A", 2e-5, 0.03)], None).unwrap();
    [300.0, 340.0, 320.0].into_iter().enumerate().map(|(i,t)|
        (i as u64, ClosedGasParcel::new(model.clone(), point(t,1e5),0.001).unwrap())).collect()
}
fn controls() -> SpeciesNetworkSolve {
    SpeciesNetworkSolve { maximum_species_fraction:0.2, mass_tolerance_kg:1e-13, energy_tolerance_j:1e-7 }
}
fn endpoints(gases: &BTreeMap<u64,ClosedGasParcel>) -> BTreeMap<u64,GasEnergySolve> {
    gases.keys().map(|&id|(id,GasEnergySolve { minimum_temperature_k:200.0,
        maximum_temperature_k:1800.0, energy_tolerance_j_kg:1e-7, maximum_iterations:128 })).collect()
}
fn row(incoming:f64,outgoing:f64,relative:f64) -> BoundarySpeciesTransfer {
    BoundarySpeciesTransfer { species:SpeciesId::new("A").unwrap(),incoming_mass_kg:incoming,
        outgoing_mass_kg:outgoing,net_mass_into_gas_kg:incoming-outgoing,relative_energy_into_gas_j:relative }
}
fn edge(mass:f64) -> SpeciesNetworkTransfer {
    SpeciesNetworkTransfer { donor:0,receiver:1,species:SpeciesId::new("A").unwrap(),mass_kg:mass }
}
fn near(a:f64,b:f64,t:f64) { assert!((a-b).abs()<=t,"{a:.17e} != {b:.17e}"); }

#[test]
fn simultaneous_boundary_internal_flow_and_heat_reconstruct_each_actual_endpoint() {
    let gases=gases(0.0); let saved=gases.clone(); let dm=2e-6; let q=0.025;
    let uptake=[row(0.0,1e-6,0.0)]; let desorb=[row(3e-6,0.0,3e-6*CP*20.0)];
    let bindings=[SpeciesNetworkBoundary { id:10,parcel:0,catalog:gases[&0].model(),species:&uptake },
        SpeciesNetworkBoundary { id:11,parcel:1,catalog:gases[&1].model(),species:&desorb }];
    let result=exchange_candidate_with_boundaries(&gases,&[edge(dm)],
        &[SpeciesNetworkHeatTransfer { parcels:[0,1],heat_left_to_right_j:q }],
        &bindings,&endpoints(&gases),controls(),||false).unwrap();
    let m0=gases[&0].mass_kg(); let m1=gases[&1].mass_kg();
    near(result.parcels()[&0].state().state().temperature,
        (m0*CV*300.0-(dm+1e-6)*CP*300.0-q)/((m0-dm-1e-6)*CV),2e-8);
    near(result.parcels()[&1].state().state().temperature,
        (m1*CV*340.0+dm*CP*300.0+3e-6*CP*360.0+q)/((m1+dm+3e-6)*CV),2e-8);
    assert_eq!(result.parcels()[&2],gases[&2]); assert_eq!(gases,saved);
    assert_eq!(result.step().boundaries.len(),2);
    near(result.step().total_mass_residual_kg,0.0,1e-13);
    near(result.step().total_energy_residual_j,0.0,1e-7);
}

#[test]
fn counterdesorption_cannot_hide_shared_gross_withdrawal() {
    let gases=gases(0.0); let m=gases[&0].mass_kg();
    let stream=[row(0.1*m,0.1*m,0.0)];
    let boundary=[SpeciesNetworkBoundary { id:7,parcel:0,catalog:gases[&0].model(),species:&stream }];
    let result=exchange_candidate_with_boundaries(&gases,&[edge(0.15*m)],&[],&boundary,
        &endpoints(&gases),controls(),||false);
    assert!(matches!(result,Err(SpeciesNetworkError::Withdrawal { id:0,.. })));
    assert!(exchange_candidate_with_boundaries(&gases,&[],&[],&boundary,&endpoints(&gases),controls(),||false).is_ok());
    assert!(exchange_candidate_with_heat(&gases,&[edge(0.15*m)],&[],&endpoints(&gases),controls(),||false).is_ok());
}

#[test]
fn thermal_target_is_invariant_under_large_species_formation_offset() {
    let run=|reference| {
        let gases=gases(reference); let streams=[row(4e-6,2e-6,0.0125)];
        let boundaries=[SpeciesNetworkBoundary { id:4,parcel:1,catalog:gases[&1].model(),species:&streams }];
        exchange_candidate_with_boundaries(&gases,&[edge(1e-6)],&[],&boundaries,
            &endpoints(&gases),controls(),||false).unwrap()
    };
    let a=run(0.0); let b=run(1e12);
    for id in 0..3 {
        near(a.parcels()[&id].state().state().temperature,b.parcels()[&id].state().state().temperature,2e-8);
        near(a.parcels()[&id].mass_kg(),b.parcels()[&id].mass_kg(),1e-15);
    }
    near(a.step().total_energy_residual_j,b.step().total_energy_residual_j,1e-7);
}

#[test]
fn previously_absent_desorbate_enters_through_its_original_source_catalog() {
    let catalog=model(vec![0.5,0.5],0.0);
    let pure=catalog.with_composition(Composition::new(CompositionBasis::MoleFraction,vec![1.0,0.0]).unwrap()).unwrap();
    let gas=ClosedGasParcel::new(pure,point(300.0,1e5),0.001).unwrap();
    let gases=BTreeMap::from([(0,gas)]); let b=SpeciesId::new("B").unwrap();
    let streams=[BoundarySpeciesTransfer { species:b.clone(),incoming_mass_kg:1e-8,
        outgoing_mass_kg:0.0,net_mass_into_gas_kg:1e-8,relative_energy_into_gas_j:0.0 }];
    let boundaries=[SpeciesNetworkBoundary { id:1,parcel:0,catalog:&catalog,species:&streams }];
    let result=exchange_candidate_with_boundaries(&gases,&[],&[],&boundaries,&endpoints(&gases),controls(),||false).unwrap();
    near(result.parcels()[&0].species_mass_kg(&b).unwrap(),1e-8,1e-13);
    assert_eq!(gases[&0].species_mass_kg(&b).unwrap(),0.0);
    assert_eq!(result.parcels()[&0].volume_m3(),gases[&0].volume_m3());
}

#[test]
fn external_owner_order_is_canonical_and_empty_boundaries_preserve_old_path() {
    let gases=gases(0.0); let first=[row(1e-6,0.0,0.01)]; let second=[row(0.0,2e-6,-0.02)];
    let a=SpeciesNetworkBoundary { id:12,parcel:0,catalog:gases[&0].model(),species:&first };
    let b=SpeciesNetworkBoundary { id:3,parcel:0,catalog:gases[&0].model(),species:&second };
    let run=|bindings:&[SpeciesNetworkBoundary<'_>]| exchange_candidate_with_boundaries(
        &gases,&[],&[],bindings,&endpoints(&gases),controls(),||false).unwrap();
    assert_eq!(run(&[a,b]),run(&[b,a]));
    assert_eq!(run(&[]),exchange_candidate(&gases,&[],&endpoints(&gases),controls(),||false).unwrap());
}

#[test]
fn invalid_sources_streams_and_cancellation_leave_every_gas_unchanged() {
    let changed=gases(1e6);
    let gases=gases(0.0); let saved=gases.clone(); let streams=[row(1e-6,0.0,0.0)];
    let bad=[SpeciesNetworkBoundary { id:1,parcel:0,catalog:changed[&0].model(),species:&streams }];
    assert!(exchange_candidate_with_boundaries(&gases,&[],&[],&bad,&endpoints(&gases),controls(),||false).is_err());
    let good=SpeciesNetworkBoundary { id:1,parcel:0,catalog:gases[&0].model(),species:&streams };
    assert!(exchange_candidate_with_boundaries(&gases,&[],&[],&[good,good],&endpoints(&gases),controls(),||false).is_err());
    let mut polls=0;
    assert!(matches!(exchange_candidate_with_boundaries(&gases,&[],&[],&[good],&endpoints(&gases),controls(),||{
        polls+=1; polls>8 }),Err(SpeciesNetworkError::Cancelled)));
    assert_eq!(gases,saved);
}
