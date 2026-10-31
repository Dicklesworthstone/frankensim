//! G0/G1/G3 analytic synthetic data, not validation for a named polymer.
use super::*;
use crate::thermochemical_gas::{ThermochemicalGasModel, tests::{model,point,species,transport}};
use crate::thermochemical_gas::energy::GasEnergySolve;
use crate::thermochemical_gas::energy::parcel::transfer::nozzle::{
    GasNozzle, network::snapshot::SnapshotNozzleLink};
use crate::thermochemical_gas::energy::parcel::transfer::species::network::{
    SpeciesNetworkError,SpeciesNetworkSolve};
use super::super::{GasTransportPaths,GasTransportSolve,GasTransportError,step_candidate_with_paths};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet,InterpolationPolicy,MaterialStateId,PropertyClaim,PropertyKey,
    PropertyValue,Provenance,UncertaintyModel};
use fs_qty::{Dims,Pressure,QuantitySpec,semantic::{QuantityKind,SemanticType,ValueForm}};
use fs_thermochem::{Composition,CompositionBasis,mixture::FrozenIdealGasMixtureModelV1};

fn typed(k:QuantityKind)->QuantitySpec {QuantitySpec::semantic(SemanticType::new(k,ValueForm::Static))}
fn source(value:PropertyValue,dims:Dims,policy:InterpolationPolicy,pressures:[f64;2])->MaterialCard {
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {key:PropertyKey::new(MOLAR_PERMEABILITY_PROPERTY,dims),value,
        validity:ValidityDomain::unconstrained()
            .with_quantity("temperature",typed(QuantityKind::AbsoluteTemperature),200.0,1800.0)
            .with_quantity("pressure",typed(QuantityKind::Pressure),pressures[0],pressures[1]),
        uncertainty:UncertaintyModel::Unstated,interpolation:policy,observations:vec![],
        provenance:Provenance {source:"synthetic Henry membrane".into(),license:"CC0-1.0".into(),artifact:None}
    }).unwrap();
    MaterialCard::assemble(MaterialStateId {chemistry:"synthetic membrane/permeant".into(),
        phase:"dense solid".into(),process:"analytic test".into(),revision:0},claims,vec![]).unwrap()
}
fn card(value:f64)->MaterialCard {
    source(PropertyValue::Scalar {value,dims:MOLAR_PERMEABILITY_DIMS},
        MOLAR_PERMEABILITY_DIMS,InterpolationPolicy::ConstantWithinValidity,[1000.0,1e7])
}
fn bindings<'a>(cards:[&'a MaterialCard;2],q:&'a QueryPoint)->[MembranePermeabilityMaterial<'a>;2] {
    ["A","B"].into_iter().zip(cards).map(|(name,material)|MembranePermeabilityMaterial {
        species:SpeciesId::new(name).unwrap(),material,context:q,
        selection:MaterialPropertySelection::SingleClaimOnly
    }).collect::<Vec<_>>().try_into().unwrap()
}
fn gas(x:f64,p:f64)->ClosedGasParcel {
    ClosedGasParcel::new(model(vec![x,1.0-x],0.0),point(400.0,p),0.001).unwrap()
}
fn pure(name:&str,mass:f64,cp:f64)->ClosedGasParcel {
    let model=ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(
        vec![species(name,mass,cp,0.0)],
        Composition::new(CompositionBasis::MoleFraction,vec![1.0]).unwrap()).unwrap(),
        vec![transport(name,2e-5,0.03)],None).unwrap();
    ClosedGasParcel::new(model,point(400.0,100000.0),0.001).unwrap()
}
fn geometry()->MembraneGeometry {MembraneGeometry {area_m2:0.02,thickness_m:1e-4}}
fn energy()->GasEnergySolve {
    GasEnergySolve {minimum_temperature_k:200.0,maximum_temperature_k:1800.0,
        energy_tolerance_j_kg:1e-7,maximum_iterations:128}
}
fn solve()->GasTransportSolve {
    GasTransportSolve {exchange:SpeciesNetworkSolve {maximum_species_fraction:0.1,
        mass_tolerance_kg:1e-13,energy_tolerance_j:1e-8},entropy_tolerance_j_k:1e-10}
}
fn link<'a>(materials:&'a [MembranePermeabilityMaterial<'a>])->MembraneNetworkLink<'a> {
    MembraneNetworkLink {id:9,parcels:[0,1],path:MembranePath::Film {geometry:geometry(),materials},
        maximum_temperature_difference_k:1e-8}
}
fn near(a:f64,b:f64) {assert!((a-b).abs()<=2e-12*b.abs().max(1e-20),"{a} != {b}");}

#[test]
fn molar_basis_requires_thickness_and_refuses_nonphysical_conversions() {
    let p=MolarPermeability::try_new(2e-14).unwrap();
    near(p.permeance(1e-4).unwrap().value_si(),2e-10);
    assert_ne!(MOLAR_PERMEABILITY_DIMS,MOLAR_PERMEANCE_DIMS);
    assert_ne!(MOLAR_PERMEABILITY_DIMS,Dims([2,0,0,0,0,0]));
    for v in [-1.0,f64::NAN,f64::INFINITY] {assert!(MolarPermeability::try_new(v).is_err());}
    for v in [0.0,-1.0,f64::NAN,f64::INFINITY] {assert!(p.permeance(v).is_err());}
    assert_eq!(MolarPermeability::try_new(0.0).unwrap().permeance(1e-4).unwrap().value_si(),0.0);
}

#[test]
fn analytic_selectivity_counterpermeation_reversal_and_thickness_scaling() {
    let (left,right)=(gas(0.8,100000.0),gas(0.2,100000.0));
    let (a,b,q)=(card(2e-14),card(5e-14),point(400.0,100000.0));
    let mats=bindings([&a,&b],&q);
    let flux=left.membrane_flux_to(&right,geometry(),&mats,0.0,||false).unwrap();
    near(flux.species[0].molar_flux_mol_m2_s,2e-14/1e-4*60000.0);
    near(flux.species[1].molar_flux_mol_m2_s,-5e-14/1e-4*60000.0);
    near(flux.species[0].mass_rate_left_to_right_kg_s,1.2e-5*0.028*0.02);
    near(flux.species[1].mass_rate_left_to_right_kg_s,-3e-5*0.004*0.02);
    let reverse=right.membrane_flux_to(&left,geometry(),&mats,0.0,||false).unwrap();
    let thick=left.membrane_flux_to(&right,MembraneGeometry {thickness_m:2e-4,..geometry()},
        &mats,0.0,||false).unwrap();
    for i in 0..2 {
        near(reverse.species[i].mass_rate_left_to_right_kg_s,-flux.species[i].mass_rate_left_to_right_kg_s);
        near(thick.species[i].mass_rate_left_to_right_kg_s,0.5*flux.species[i].mass_rate_left_to_right_kg_s);
    }
    assert_ne!(flux.species.iter().map(|s|s.molar_flux_mol_m2_s).sum::<f64>(),0.0);
}

#[test]
fn partial_pressure_drives_flow_against_the_total_pressure_difference() {
    let (left,right)=(gas(0.1,200000.0),gas(0.9,100000.0));
    let (a,b,q)=(card(2e-14),card(5e-14),point(400.0,100000.0));
    let mats=bindings([&a,&b],&q);
    let flux=left.membrane_flux_to(&right,geometry(),&mats,0.0,||false).unwrap();
    assert_eq!(flux.species[0].partial_pressures_pa,[20000.0,90000.0]);
    assert!(flux.species[0].mass_rate_left_to_right_kg_s<0.0);
    assert!(flux.species[1].mass_rate_left_to_right_kg_s>0.0);
}

#[test]
fn absent_species_closed_paths_and_explicit_impermeability_are_exact() {
    let (left,right)=(pure("A",0.028,3.5),pure("B",0.004,2.5));
    let (a,b,q)=(card(2e-14),card(0.0),point(400.0,100000.0));
    let mats=bindings([&a,&b],&q);
    let flux=left.membrane_flux_to(&right,geometry(),&mats,0.0,||false).unwrap();
    assert_eq!(flux.species[0].partial_pressures_pa,[100000.0,0.0]);
    assert_eq!(flux.species[1].partial_pressures_pa,[0.0,100000.0]);
    assert_eq!(flux.species[1].mass_rate_left_to_right_kg_s,0.0);
    assert!(left.membrane_flux_to(&right,geometry(),&mats[..1],0.0,||false).is_err());
    assert!(left.membrane_flux_to(&right,MembraneGeometry {area_m2:0.0,..geometry()},
        &[],0.0,||false).unwrap().species.is_empty());
}

#[test]
fn source_support_checks_both_boundaries_and_rejects_wrong_units_and_pressure_laws() {
    let (left,right)=(gas(0.75,200000.0),gas(0.25,10000.0));
    let (b,q)=(card(5e-14),point(400.0,105000.0));
    let narrow=source(PropertyValue::Scalar {value:2e-14,dims:MOLAR_PERMEABILITY_DIMS},
        MOLAR_PERMEABILITY_DIMS,InterpolationPolicy::ConstantWithinValidity,[90000.0,110000.0]);
    assert!(left.membrane_flux_to(&right,geometry(),&bindings([&narrow,&b],&q),0.0,||false).is_err());
    let wrong=source(PropertyValue::Scalar {value:2e-14,dims:Dims([2,0,0,0,0,0])},
        Dims([2,0,0,0,0,0]),InterpolationPolicy::ConstantWithinValidity,[1000.0,1e7]);
    assert!(left.membrane_flux_to(&right,geometry(),&bindings([&wrong,&b],&q),0.0,||false).is_err());
    let pressure=source(PropertyValue::Curve {abscissa:"pressure".into(),abscissa_dims:Pressure::DIMS,
        knots:vec![(1000.0,2e-14),(1e7,2e-14)],dims:MOLAR_PERMEABILITY_DIMS},
        MOLAR_PERMEABILITY_DIMS,InterpolationPolicy::LinearInside,[1000.0,1e7]);
    assert!(left.membrane_flux_to(&right,geometry(),&bindings([&pressure,&b],&q),0.0,||false).is_err());
    // A stale legacy T=300 must not qualify actual gas at canonical T=400.
    let base=card(2e-14);
    let mut claim=base.claims_for(MOLAR_PERMEABILITY_PROPERTY)[0].1.clone();
    claim.validity=claim.validity.with("T",290.0,310.0);
    let mut claims=ClaimSet::new();claims.insert_claim(claim).unwrap();
    let legacy=MaterialCard::assemble(base.id().clone(),claims,vec![]).unwrap();
    let stale=q.clone().with("T",300.0).unwrap();
    assert!(left.membrane_flux_to(&right,geometry(),&bindings([&legacy,&b],&stale),0.0,||false).is_err());
}

#[test]
fn combined_owner_conserves_selective_mass_species_and_partial_enthalpy() {
    let parcels=BTreeMap::from([(0,pure("A",0.028,3.5)),(1,pure("B",0.004,2.5))]);
    let saved=parcels.clone();
    let endpoints=parcels.keys().map(|&id|(id,energy())).collect();
    let (a,b,q)=(card(2e-14),card(5e-14),point(400.0,100000.0));
    let mats=bindings([&a,&b],&q);let links=[link(&mats)];let dt=0.01;
    let candidate=step_candidate_with_paths(&parcels,GasTransportPaths {nozzles:&[],
        diffusion:&[],porous:&[],membranes:&links,..Default::default()},dt,&endpoints,solve(),||false).unwrap();
    let step=candidate.step();
    assert_eq!(step.exchange.fluxes.len(),2);
    assert!(step.exchange.total_mass_residual_kg.abs()<1e-13);
    assert!(step.exchange.maximum_species_mass_residual_kg<1e-13);
    assert!(step.exchange.total_energy_residual_j.abs()<1e-8);
    assert!(step.entropy_change_j_k>=-1e-10);
    for row in &step.exchange.fluxes {
        let upstream=&parcels[&row.transfer.donor];
        let species=&upstream.model().thermodynamics().components()[0];
        let h=species.evaluate(fs_qty::Temperature::new(400.0)).unwrap().properties().h().value()
            /species.molar_mass().value();
        near(row.source_enthalpy_j_kg,h);
        near(row.carried_enthalpy_j,row.transfer.mass_kg*h);
    }
    for p in candidate.parcels().values() {assert_eq!(p.model().thermodynamics().components().len(),2);}
    assert_eq!(parcels,saved);
}

#[test]
fn membrane_and_nozzle_cannot_jointly_spend_an_original_species_allowance() {
    let parcels=BTreeMap::from([(0,gas(0.75,100000.0)),(1,gas(0.25,100000.0)),
        (2,gas(0.75,50000.0))]);
    let endpoints=parcels.keys().map(|&id|(id,energy())).collect();
    let (a,b,q)=(card(2e-14),card(0.0),point(400.0,100000.0));
    let mats=bindings([&a,&b],&q);let dt=0.01;let mut film=link(&mats);
    let flux=parcels[&0].membrane_flux_to(&parcels[&1],geometry(),&mats,0.0,||false).unwrap();
    let initial_a=parcels[&0].mass_kg()*0.75*0.028/(0.75*0.028+0.25*0.004);
    let MembranePath::Film { geometry, .. } = &mut film.path else { unreachable!() };
    geometry.area_m2*=0.06*initial_a/(flux.species[0].mass_rate_left_to_right_kg_s*dt);
    let mut nozzle=SnapshotNozzleLink {id:3,parcels:[0,2],
        nozzle:GasNozzle {area_m2:1e-8,discharge_coefficient:0.8},nozzle_solve:energy()};
    let rate=parcels[&0].nozzle_flow_to(50000.0,nozzle.nozzle,energy(),||false).unwrap().mass_flow_kg_s;
    nozzle.nozzle.area_m2*=0.06*parcels[&0].mass_kg()/(rate*dt);
    step_candidate_with_paths(&parcels,GasTransportPaths {nozzles:&[nozzle],
        diffusion:&[],porous:&[],membranes:&[],..Default::default()},dt,&endpoints,solve(),||false).unwrap();
    step_candidate_with_paths(&parcels,GasTransportPaths {nozzles:&[],
        diffusion:&[],porous:&[],membranes:&[film.clone()],..Default::default()},dt,&endpoints,solve(),||false).unwrap();
    assert!(matches!(step_candidate_with_paths(&parcels,GasTransportPaths {nozzles:&[nozzle],
        diffusion:&[],porous:&[],membranes:&[film],..Default::default()},dt,&endpoints,solve(),||false),
        Err(GasTransportError::Exchange(SpeciesNetworkError::Withdrawal {id:0,..}))));
}

#[test]
fn final_partial_pressure_crossing_and_cancellation_preserve_inputs() {
    let parcels=BTreeMap::from([(0,gas(0.8,100000.0)),(1,gas(0.2,100000.0))]);
    let saved=parcels.clone();let (a,b,q)=(card(2e-14),card(5e-14),point(400.0,100000.0));
    let mats=bindings([&a,&b],&q);let links=[link(&mats)];
    let (_,reports)=sample_transfers(&parcels,&links,0.01,||false).unwrap();
    admit_partial_pressures(&parcels,&reports,||false).unwrap();
    let crossed=BTreeMap::from([(0,gas(0.2,100000.0)),(1,gas(0.8,100000.0))]);
    assert!(matches!(admit_partial_pressures(&crossed,&reports,||false),
        Err(MembraneNetworkError::EquilibriumCrossed {id:9,..})));
    assert_eq!(sample_transfers(&parcels,&links,0.01,||true),Err(MembraneNetworkError::Cancelled));
    let mut polls=0;
    assert_eq!(sample_transfers(&parcels,&links,0.01,||{polls+=1;polls==8}),
        Err(MembraneNetworkError::Cancelled));
    assert_eq!(parcels,saved);
}

mod series;

