//! Actual material closures and original transport owners; synthetic data only.
use super::*;
use crate::thermochemical_gas::tests::{model, point};
use crate::thermochemical_gas::energy::parcel::transfer::snapshot::GasSnapshotSolve;
use crate::thermochemical_gas::energy::parcel::transfer::nozzle::GasNozzle;
use crate::thermochemical_gas::energy::parcel::transfer::species::{GasSpeciesExchangeSolve,
    diffusion::{BinaryDiffusionMaterial, GasDiffusionGeometry, GasDiffusionSolve,
        BINARY_DIFFUSIVITY_DIMS, BINARY_DIFFUSIVITY_PROPERTY,
        network::DiffusionFaceSolve}};
use crate::state_point::MaterialPropertySelection;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_qty::{Temperature, QuantitySpec, semantic::{QuantityKind, SemanticType, ValueForm}};
use fs_thermochem::SpeciesId;

fn energy() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 }
}
fn solve() -> GasTransportSolve {
    GasTransportSolve { exchange: SpeciesNetworkSolve { maximum_species_fraction: 0.1,
        mass_tolerance_kg: 1e-13, energy_tolerance_j: 1e-8 }, entropy_tolerance_j_k: 1e-10 }
}
fn gas(x: f64, p: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(model(vec![x, 1.0-x], 0.0), point(400.0,p), 0.001).unwrap()
}
fn gases() -> BTreeMap<u64, ClosedGasParcel> {
    BTreeMap::from([(0,gas(0.75,100000.0)),(1,gas(0.25,100000.0)),
        (2,gas(0.25,50000.0)),(3,gas(0.75,50000.0))])
}
fn endpoints(p: &BTreeMap<u64, ClosedGasParcel>) -> BTreeMap<u64, GasEnergySolve> {
    p.keys().map(|&id|(id,energy())).collect()
}
fn nozzle(id: u64, parcels: [u64;2]) -> SnapshotNozzleLink {
    SnapshotNozzleLink { id, parcels, nozzle: GasNozzle { area_m2: 2e-8, discharge_coefficient: 0.8 },
        nozzle_solve: energy() }
}
fn card() -> MaterialCard {
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(BINARY_DIFFUSIVITY_PROPERTY,BINARY_DIFFUSIVITY_DIMS),
        value: PropertyValue::Curve { abscissa:"temperature".into(),abscissa_dims:Temperature::DIMS,
            knots:vec![(200.0,2e-5),(1800.0,2e-5)],dims:BINARY_DIFFUSIVITY_DIMS },
        // Typed like the query's temperature axis, as diffusion/tests.rs types the same claim.
        validity:ValidityDomain::unconstrained()
            .with_quantity("temperature",QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature,ValueForm::Static)),200.0,1800.0),
        uncertainty:UncertaintyModel::Unstated,
        interpolation:InterpolationPolicy::LinearInside, observations:vec![],
        provenance:Provenance {source:"synthetic mixed-path diffusivity".into(),license:"CC0-1.0".into(),artifact:None},
    }).unwrap();
    MaterialCard::assemble(MaterialStateId {chemistry:"synthetic A/B pair".into(),phase:"gas".into(),
        process:"test".into(),revision:0},claims,vec![]).unwrap()
}
fn pair<'a>(d: &'a MaterialCard, p: &'a QueryPoint) -> [BinaryDiffusionMaterial<'a>;1] {
    [BinaryDiffusionMaterial {species:[SpeciesId::new("A").unwrap(),SpeciesId::new("B").unwrap()],
        material:d,context:p,selection:MaterialPropertySelection::SingleClaimOnly}]
}
fn face<'a>(pairs: &'a [BinaryDiffusionMaterial<'a>]) -> DiffusionNetworkLink<'a> {
    DiffusionNetworkLink {id:7,parcels:[0,1],pairs,
        geometry:GasDiffusionGeometry {area_m2:0.01,distance_m:0.01},
        solve:DiffusionFaceSolve {maximum_temperature_difference_k:1e-5,maximum_pressure_difference_pa:0.01,
            gradient_tolerance_per_m:1e-9,molar_flux_tolerance_mol_m2_s:1e-12}}
}
fn molecular_solve() -> DiffusionNetworkSolve {
    DiffusionNetworkSolve {exchange:solve().exchange,entropy_tolerance_j_k:solve().entropy_tolerance_j_k}
}
fn near(a:f64,b:f64,tol:f64) {assert!((a-b).abs()<=tol,"{a} != {b}");}

#[test]
fn empty_and_single_mechanism_paths_match_original_physical_owners() {
    let p=gases();let bounds=endpoints(&p);let d=card();let q=point(400.0,100000.0);let pairs=pair(&d,&q);
    let idle=step_candidate(&p,&[],&[],0.01,&bounds,solve(),||false).unwrap();
    assert_eq!(idle.parcels(),&p);
    let diff=molecular::step_candidate(&p,&[face(&pairs)],0.01,&bounds,molecular_solve(),||false).unwrap();
    let only=step_candidate(&p,&[],&[face(&pairs)],0.01,&bounds,solve(),||false).unwrap();
    assert_eq!(only.parcels(),diff.parcels());assert_eq!(only.step().exchange,diff.step().exchange);
    let link=nozzle(9,[0,2]);
    let native=pressure::step_candidate(&p,&[link],0.01,&bounds,GasSnapshotSolve {
        maximum_withdrawn_fraction:0.1,mass_tolerance_kg:1e-13,energy_tolerance_j:1e-8},||false).unwrap();
    let only=step_candidate(&p,&[link],&[],0.01,&bounds,solve(),||false).unwrap();
    assert_eq!(only.step().nozzles,native.step().links);
    for id in p.keys() {
        near(only.parcels()[id].mass_kg(),native.parcels()[id].mass_kg(),1e-13);
        near(only.parcels()[id].state().state().temperature,native.parcels()[id].state().state().temperature,1e-7);
        near(only.parcels()[id].state().state().pressure,native.parcels()[id].state().state().pressure,1e-5);
    }
}

#[test]
fn both_mechanisms_read_original_gases_and_close_one_species_energy_update() {
    let p=gases();let bounds=endpoints(&p);let d=card();let q=point(400.0,100000.0);let pairs=pair(&d,&q);
    let n=nozzle(9,[0,2]);let f=face(&pairs);let dt=0.01;
    let candidate=step_candidate(&p,&[n],&[f.clone()],dt,&bounds,solve(),||false).unwrap();
    let original_nozzle=p[&0].nozzle_flow_to(50000.0,n.nozzle,energy(),||false).unwrap();
    assert_eq!(candidate.step().nozzles[0].flow,original_nozzle);
    let original_diff=p[&0].diffusion_flux_to(&p[&1],f.geometry,&pairs,GasDiffusionSolve {
        maximum_temperature_difference_k:f.solve.maximum_temperature_difference_k,
        maximum_pressure_difference_pa:f.solve.maximum_pressure_difference_pa,
        gradient_tolerance_per_m:f.solve.gradient_tolerance_per_m,
        molar_flux_tolerance_mol_m2_s:f.solve.molar_flux_tolerance_mol_m2_s,
        entropy_tolerance_j_k:solve().entropy_tolerance_j_k,
        exchange:GasSpeciesExchangeSolve {gas:[energy();2],maximum_species_fraction:0.1,
            mass_tolerance_kg:1e-13,energy_tolerance_j:1e-8}},||false).unwrap();
    assert_eq!(candidate.step().diffusion[0].flux,original_diff);
    let mut rows=inventory(p[&0].model(),original_nozzle.mass_flow_kg_s*dt).unwrap().into_iter()
        .map(|(species,mass_kg)|SpeciesNetworkTransfer {donor:0,receiver:2,species,mass_kg}).collect::<Vec<_>>();
    for flux in &original_diff.species {
        let dm=flux.mass_rate_left_to_right_kg_s*dt;
        if dm!=0.0 {rows.push(SpeciesNetworkTransfer {donor:if dm>0.0 {0}else{1},
            receiver:if dm>0.0 {1}else{0},species:flux.species.clone(),mass_kg:dm.abs()});}
    }
    let expected=super::super::exchange_candidate(&p,&rows,&bounds,solve().exchange,||false).unwrap();
    assert_eq!(candidate.parcels(),expected.parcels());assert_eq!(candidate.step().exchange,*expected.step());
    for id in p.keys() {assert_eq!(candidate.parcels()[id].volume_m3(),p[id].volume_m3());}
    assert!(candidate.step().exchange.total_energy_residual_j.abs()<=solve().exchange.energy_tolerance_j);
    assert!(candidate.step().entropy_change_j_k>=-solve().entropy_tolerance_j_k);
    assert_eq!(candidate.parcels()[&3],p[&3]);
}

#[test]
fn separate_safe_mechanisms_cannot_jointly_spend_the_same_species_inventory() {
    let p=gases();let bounds=endpoints(&p);let d=card();let q=point(400.0,100000.0);let pairs=pair(&d,&q);
    let dt=0.1;let mut n=nozzle(9,[0,2]);let mut f=face(&pairs);
    let sample=step_candidate(&p,&[n],&[f.clone()],dt,&bounds,solve(),||false).unwrap();
    n.nozzle.area_m2*=0.06*p[&0].mass_kg()/sample.step().nozzles[0].transferred_mass_kg;
    let a=SpeciesId::new("A").unwrap();
    let initial_a=inventory(p[&0].model(),p[&0].mass_kg()).unwrap()[&a];
    let rate=sample.step().diffusion[0].flux.species.iter().find(|s|s.species==a).unwrap().mass_rate_left_to_right_kg_s;
    f.geometry.area_m2*=0.06*initial_a/(rate*dt);
    step_candidate(&p,&[n],&[],dt,&bounds,solve(),||false).unwrap();
    step_candidate(&p,&[],&[f.clone()],dt,&bounds,solve(),||false).unwrap();
    assert!(matches!(step_candidate(&p,&[n],&[f],dt,&bounds,solve(),||false),
        Err(GasTransportError::Exchange(SpeciesNetworkError::Withdrawal {id:0,..}))));
}

#[test]
fn reporting_ids_and_slice_order_do_not_select_a_transport_operator_order() {
    let p=gases();let bounds=endpoints(&p);let d=card();let q=point(400.0,100000.0);let pairs=pair(&d,&q);
    let mut second=face(&pairs);second.id=8;
    let a=step_candidate(&p,&[nozzle(9,[0,2]),nozzle(4,[1,3])],&[face(&pairs),second.clone()],
        0.01,&bounds,solve(),||false).unwrap();
    let mut first=face(&pairs);first.id=900;second.id=2;
    let b=step_candidate(&p,&[nozzle(1,[1,3]),nozzle(500,[0,2])],&[second,first],
        0.01,&bounds,solve(),||false).unwrap();
    assert_eq!(a.parcels(),b.parcels());assert_eq!(a.step().exchange,b.step().exchange);
}

#[test]
fn invalid_source_late_endpoint_and_cancellation_preserve_all_inputs() {
    let p=gases();let saved=p.clone();let bounds=endpoints(&p);let d=card();let q=point(400.0,100000.0);let pairs=pair(&d,&q);
    let n=nozzle(9,[0,2]);let f=face(&pairs);
    let mut missing=f.clone();missing.pairs=&[];
    assert!(matches!(step_candidate(&p,&[n],&[missing],0.01,&bounds,solve(),||false),
        Err(GasTransportError::Diffusion(_))));
    let mut narrow=bounds.clone();narrow.get_mut(&2).unwrap().maximum_temperature_k=250.0;
    assert!(step_candidate(&p,&[n],&[f.clone()],0.01,&narrow,solve(),||false).is_err());
    let mut polls=0;
    step_candidate(&p,&[n],&[f.clone()],0.01,&bounds,solve(),||{polls+=1;false}).unwrap();
    for stop in [1,polls/2,polls] {
        let mut calls=0;
        assert_eq!(step_candidate(&p,&[n],&[f.clone()],0.01,&bounds,solve(),||{calls+=1;calls==stop}),
            Err(GasTransportError::Cancelled));
    }
    assert_eq!(p,saved);
}

#[test]
fn combined_material_states_feed_the_original_heat_work_and_load_owner() {
    let p=gases();let bounds=endpoints(&p);let d=card();let q=point(400.0,100000.0);let pairs=pair(&d,&q);
    let result=step_candidate(&p,&[nozzle(9,[0,2])],&[face(&pairs)],0.01,&bounds,solve(),||false).unwrap();
    let gas=&result.parcels()[&0];
    let moved=gas.midpoint_heat_volume(0.00105,1.0,energy(),||false).unwrap();
    let gradient=[0.01];
    let load=moved.discrete_mechanical_load(0.00105,&gradient,100000.0).unwrap();
    assert!(load.generalized_forces().unwrap()[0].is_finite());
    near(moved.step().stored_energy_change_j+moved.step().work_by_gas_j,1.0,1e-8);
    assert_eq!(moved.endpoint().model(),gas.model());
}
