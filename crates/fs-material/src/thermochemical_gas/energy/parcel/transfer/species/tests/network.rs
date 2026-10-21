//! Same source fixtures as the pair tests; exercise the shared multiport update.
use super::*;
use super::super::network::{self,SpeciesNetworkTransfer as Transfer,SpeciesNetworkSolve,SpeciesNetworkError};
fn solve() -> SpeciesNetworkSolve {
    let c=controls(); SpeciesNetworkSolve {maximum_species_fraction:c.maximum_species_fraction,
        mass_tolerance_kg:c.mass_tolerance_kg,energy_tolerance_j:c.energy_tolerance_j}
}
fn flux(donor:u64,receiver:u64,name:&str,mass_kg:f64)->Transfer {
    Transfer {donor,receiver,species:SpeciesId::new(name).unwrap(),mass_kg}
}
fn endpoints(gas:&BTreeMap<u64,ClosedGasParcel>)->BTreeMap<u64,GasEnergySolve> {
    gas.keys().map(|&id|(id,energy())).collect()
}
#[test]
fn selective_network_matches_pair_and_is_independent_of_transfer_order() {
    let gas=BTreeMap::from([(1,parcel(0.75,500.0,0.001)),(2,parcel(0.25,350.0,0.001))]);
    let pair=gas[&1].exchange_species_with(&gas[&2],&[row("A",1e-5),row("B",-2e-6)],controls(),||false).unwrap();
    let rows=[flux(1,2,"A",1e-5),flux(2,1,"B",2e-6)];
    let trial=network::exchange_candidate(&gas,&rows,&endpoints(&gas),solve(),||false).unwrap();
    let reversed=[rows[1].clone(),rows[0].clone()];
    assert_eq!(trial,network::exchange_candidate(&gas,&reversed,&endpoints(&gas),solve(),||false).unwrap());
    for (id,reference) in [(1,pair.left()),(2,pair.right())] {
        let actual=&trial.parcels()[&id];
        near(actual.mass_kg(),reference.mass_kg(),1e-14);
        near(actual.state().state().temperature,reference.state().state().temperature,1e-7);
        near(actual.state().state().pressure,reference.state().state().pressure,1e-5);
        assert_eq!(actual.point().axes().get("local-lot"),gas[&id].point().axes().get("local-lot"));
    }
    assert!(trial.step().total_energy_residual_j.abs()<solve().energy_tolerance_j);
}
#[test]
fn incoming_species_cannot_depart_a_second_face_in_the_same_snapshot() {
    let gas=BTreeMap::from([(0,pure("A",0.028,3.5)),(1,pure("B",0.004,2.5)),(2,pure("B",0.004,2.5))]);
    let rows=[flux(0,1,"A",0.02*gas[&0].mass_kg()),flux(1,2,"B",0.03*gas[&1].mass_kg())];
    let trial=network::exchange_candidate(&gas,&rows,&endpoints(&gas),solve(),||false).unwrap();
    assert_eq!(trial.parcels()[&1].model().thermodynamics().components().len(),2);
    assert_eq!(trial.parcels()[&2].model().thermodynamics().components().len(),1);
    let forbidden=[rows[0].clone(),flux(1,2,"A",0.001*gas[&0].mass_kg())];
    assert!(network::exchange_candidate(&gas,&forbidden,&endpoints(&gas),solve(),||false).is_err());
}
#[test]
fn multiple_outlets_cannot_hide_gross_species_depletion_with_inflow() {
    let gas=BTreeMap::from([(0,parcel(0.75,400.0,0.001)),(1,parcel(0.75,400.0,0.001)),(2,parcel(0.75,400.0,0.001))]);
    let a=inventory(gas[&0].model(),gas[&0].mass_kg()).unwrap()[&SpeciesId::new("A").unwrap()];
    let rows=[flux(0,1,"A",0.15*a),flux(0,2,"A",0.15*a),flux(1,0,"A",0.1*a)];
    for one in &rows[..2] {
        network::exchange_candidate(&gas,std::slice::from_ref(one),&endpoints(&gas),solve(),||false).unwrap();
    }
    assert!(matches!(network::exchange_candidate(&gas,&rows,&endpoints(&gas),solve(),||false),
        Err(SpeciesNetworkError::Withdrawal {id:0,..})));
}
#[test]
fn net_zero_exchange_preserves_state_and_late_failures_preserve_every_owner() {
    let gas=BTreeMap::from([(0,parcel(0.75,400.0,0.001)),(1,parcel(0.75,400.0,0.001)),(2,parcel(0.25,400.0,0.001))]);
    let saved=gas.clone(); let rows=[flux(0,1,"A",1e-6),flux(1,0,"A",1e-6)];
    let trial=network::exchange_candidate(&gas,&rows,&endpoints(&gas),solve(),||false).unwrap();
    assert_eq!(trial.parcels(),&gas);
    assert!(trial.step().parcels.iter().all(|p|p.iterations==0));
    let rows=[flux(0,1,"A",1e-6),flux(1,2,"B",1e-7)];
    let mut count=0;
    network::exchange_candidate(&gas,&rows,&endpoints(&gas),solve(),||{count+=1;false}).unwrap();
    let mut calls=0;
    assert_eq!(network::exchange_candidate(&gas,&rows,&endpoints(&gas),solve(),||{calls+=1;calls==count}),Err(SpeciesNetworkError::Cancelled));
    let mut bounds=endpoints(&gas); bounds.get_mut(&2).unwrap().maximum_temperature_k=250.0;
    assert!(network::exchange_candidate(&gas,&rows,&bounds,solve(),||false).is_err());
    assert_eq!(gas,saved);
}
