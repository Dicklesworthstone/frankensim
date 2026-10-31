//! Multiple material diffusion faces on one unchanged physical snapshot.

mod adaptive;

#[test]
fn graph_and_combined_transport_cancel_arbitrary_species_entropy_datums() {
    use crate::thermochemical_gas::energy::parcel::transfer::species::network::transport::{self, GasTransportSolve};
    let d = card(2e-5); let p = point(400.0, 100000.0);
    let pairs = [pair("A", "B", &d, &p)];
    let links = [link(1, [0, 1], &pairs), link(2, [1, 2], &pairs)];
    let gases = |offsets| BTreeMap::from([(0, entropy_gas(0.75, offsets)),
        (1, entropy_gas(0.5, offsets)), (2, entropy_gas(0.25, offsets))]);
    let ordinary = gases([0.0; 2]);
    let shifted = gases([1e24, -3e24]);
    let expected = network::step_candidate(&ordinary, &links, 0.1,
        &endpoints(&ordinary), solve(), || false).unwrap();
    let actual = network::step_candidate(&shifted, &links, 0.1,
        &endpoints(&shifted), solve(), || false).unwrap();
    let combined = transport::step_candidate(&shifted, &[], &links, 0.1,
        &endpoints(&shifted), GasTransportSolve { exchange: solve().exchange,
            entropy_tolerance_j_k: solve().entropy_tolerance_j_k }, || false).unwrap();
    assert!(expected.step().entropy_change_j_k > 0.0);
    assert_eq!(actual.step().entropy_change_j_k.to_bits(), expected.step().entropy_change_j_k.to_bits());
    assert_eq!(combined.step().entropy_change_j_k.to_bits(), expected.step().entropy_change_j_k.to_bits());
    assert_eq!(combined.parcels(), actual.parcels());
    for id in ordinary.keys() {
        assert_eq!(expected.parcels()[id].state().state(), actual.parcels()[id].state().state());
        assert_eq!(expected.parcels()[id].mass_kg().to_bits(), actual.parcels()[id].mass_kg().to_bits());
    }
}

#[test]
fn disconnected_gas_systems_keep_distinct_source_models_for_the_same_species_name() {
    let d = card(2e-5); let p = point(400.0, 100000.0);
    let pairs = [pair("A", "B", &d, &p)];
    let first = BTreeMap::from([(0, entropy_gas(0.75, [0.0; 2])),
        (1, entropy_gas(0.25, [0.0; 2]))]);
    let expected = network::step_candidate(&first, &[link(1, [0, 1], &pairs)], 0.1,
        &endpoints(&first), solve(), || false).unwrap();
    let mut both = first;
    both.insert(2, entropy_gas(0.75, [1e24, -3e24]));
    both.insert(3, entropy_gas(0.25, [1e24, -3e24]));
    let actual = network::step_candidate(&both,
        &[link(1, [0, 1], &pairs), link(2, [2, 3], &pairs)], 0.1,
        &endpoints(&both), solve(), || false).unwrap();
    near(actual.step().entropy_change_j_k, 2.0 * expected.step().entropy_change_j_k, 1e-13);
    assert!(actual.step().exchange.maximum_species_mass_residual_kg <= solve().exchange.mass_tolerance_kg);
}

use super::*;
use super::super::network::{self,DiffusionFaceSolve,DiffusionNetworkLink,DiffusionNetworkSolve,DiffusionNetworkError};
use crate::thermochemical_gas::energy::parcel::transfer::species::network::{SpeciesNetworkSolve,SpeciesNetworkError};
fn face()->DiffusionFaceSolve {
    let c=controls(); DiffusionFaceSolve {maximum_temperature_difference_k:c.maximum_temperature_difference_k,
        maximum_pressure_difference_pa:c.maximum_pressure_difference_pa,gradient_tolerance_per_m:c.gradient_tolerance_per_m,
        molar_flux_tolerance_mol_m2_s:c.molar_flux_tolerance_mol_m2_s}
}
fn solve()->DiffusionNetworkSolve {
    let c=controls(); DiffusionNetworkSolve {exchange:SpeciesNetworkSolve {
        maximum_species_fraction:c.exchange.maximum_species_fraction,mass_tolerance_kg:c.exchange.mass_tolerance_kg,
        energy_tolerance_j:c.exchange.energy_tolerance_j},entropy_tolerance_j_k:c.entropy_tolerance_j_k}
}
fn endpoints(gases:&BTreeMap<u64,ClosedGasParcel>)->BTreeMap<u64,GasEnergySolve> {
    gases.keys().map(|&id|(id,controls().exchange.gas[0])).collect()
}
fn link<'a>(id:u64,parcels:[u64;2],pairs:&'a[BinaryDiffusionMaterial<'a>])->DiffusionNetworkLink<'a> {
    DiffusionNetworkLink {id,parcels,geometry:geometry(),pairs,solve:face()}
}
#[test]
fn single_face_matches_pair_and_symmetric_outlets_ignore_reporting_order() {
    let d=card(2e-5);let p=point(400.0,100000.0);let pairs=[pair("A","B",&d,&p)];
    let two=BTreeMap::from([(0,gas(0.75)),(1,gas(0.25))]);
    let expected=two[&0].diffuse_with(&two[&1],geometry(),0.1,&pairs,controls(),||false).unwrap();
    let actual=network::step_candidate(&two,&[link(1,[0,1],&pairs)],0.1,&endpoints(&two),solve(),||false).unwrap();
    for (id,pair_state) in [(0,expected.exchange().left()),(1,expected.exchange().right())] {
        near(actual.parcels()[&id].state().state().temperature,pair_state.state().state().temperature,1e-7);
        near(actual.parcels()[&id].state().state().sound_speed,pair_state.state().state().sound_speed,1e-7);
    }
    let mut three=two;three.insert(2,three[&1].clone());
    let links=[link(9,[0,1],&pairs),link(2,[0,2],&pairs)];
    let trial=network::step_candidate(&three,&links,0.1,&endpoints(&three),solve(),||false).unwrap();
    assert_eq!(trial.parcels()[&1],trial.parcels()[&2]);
    let renamed=[link(300,[0,2],&pairs),link(1,[0,1],&pairs)];
    let other=network::step_candidate(&three,&renamed,0.1,&endpoints(&three),solve(),||false).unwrap();
    assert_eq!(trial.parcels(),other.parcels());
    assert_eq!(trial.step().exchange,other.step().exchange);
}
#[test]
fn chain_relaxation_refines_to_the_original_diffusion_operator_eigenmode() {
    let d=card(2e-5);let p=point(400.0,100000.0);let pairs=[pair("A","B",&d,&p)];
    let links=[link(1,[0,1],&pairs),link(2,[1,2],&pairs)];
    let mut errors=Vec::new();
    for count in [16,32,64] {
        let mut gases=BTreeMap::from([(0,gas(0.75)),(1,gas(0.5)),(2,gas(0.25))]);
        let bounds=endpoints(&gases);let initial=gases.clone();
        for _ in 0..count {
            let trial=network::step_candidate(&gases,&links,10.0/f64::from(count),&bounds,solve(),||false).unwrap();
            assert!(trial.step().entropy_change_j_k>=-solve().entropy_tolerance_j_k);
            assert!(trial.step().exchange.total_energy_residual_j.abs()<solve().exchange.energy_tolerance_j);
            gases=trial.into_parts().0;
        }
        let x=|id|gases[&id].model().thermodynamics().mole_composition().fractions()[0];
        errors.push((x(0)-x(2)-0.5*(-0.2_f64).exp()).abs());
        near(x(1),0.5,1e-10);
        assert_ne!(gases[&0].state().state().sound_speed,initial[&0].state().state().sound_speed);
    }
    assert!(errors.windows(2).all(|e|e[0]/e[1]>1.9 && e[0]/e[1]<2.1),"{errors:?}");
}
#[test]
fn jointly_excessive_faces_and_late_cancellation_return_no_network() {
    let d=card(2e-5);let p=point(400.0,100000.0);let pairs=[pair("A","B",&d,&p)];
    let gases=BTreeMap::from([(0,gas(0.75)),(1,gas(0.25)),(2,gas(0.25))]);
    let bounds=endpoints(&gases);let saved=gases.clone();let links=[link(1,[0,1],&pairs),link(2,[0,2],&pairs)];
    for one in &links {
        network::step_candidate(&gases,std::slice::from_ref(one),11.25,&bounds,solve(),||false).unwrap();
    }
    assert!(matches!(network::step_candidate(&gases,&links,11.25,&bounds,solve(),||false),
        Err(DiffusionNetworkError::Exchange(SpeciesNetworkError::Withdrawal {id:0,..}))));
    let mut count=0;
    network::step_candidate(&gases,&links,0.1,&bounds,solve(),||{count+=1;false}).unwrap();
    let mut calls=0;
    assert_eq!(network::step_candidate(&gases,&links,0.1,&bounds,solve(),||{calls+=1;calls==count}),Err(DiffusionNetworkError::Cancelled));
    assert_eq!(gases,saved);
}
#[test]
fn closed_faces_leave_states_bitwise_unchanged_without_missing_pair_defaults() {
    let gases=BTreeMap::from([(0,gas(0.75)),(1,gas(0.25))]);
    let closed=DiffusionNetworkLink {id:1,parcels:[0,1],geometry:GasDiffusionGeometry {area_m2:0.0,..geometry()},pairs:&[],solve:face()};
    let trial=network::step_candidate(&gases,&[closed],0.1,&endpoints(&gases),solve(),||false).unwrap();
    assert_eq!(trial.parcels(),&gases);
    assert_eq!(trial.step().entropy_change_j_k,0.0);
    assert_eq!(trial.step().exchange.fluxes.len(),0);
}
