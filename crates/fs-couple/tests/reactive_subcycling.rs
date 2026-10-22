//! Source-backed chemistry and the original piston; synthetic isomers only.
#[path = "material_exchange/gas.rs"]
#[allow(dead_code)]
mod source;
use fs_couple::cavity_phs::piston::{GasPiston, reactive::{GasReaction, GasReactionError, GasReactionSolve}};
use fs_couple::cavity_phs::piston::reactive::subcycling::ReactionSubcycling;
use fs_couple::cavity_phs::piston::reactive::material_rate::{MaterialRateScale, ACTIVITY_RATE_SCALE_PROPERTY, ACTIVITY_RATE_SCALE_DIMS};
use fs_material::thermochemical_gas::{ThermochemicalGasModel, energy::parcel::ClosedGasParcel};
use fs_material::state_point::MaterialPropertySelection;
use fs_evidence::ValidityDomain;
use fs_matdb::*;
use fs_qty::{Pressure, Temperature, QuantitySpec, semantic::{QuantityKind, SemanticType, ValueForm}};
use fs_thermochem::*;
use fs_ports::PortedModule;
use fs_thermochem::equilibrium::*;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;

fn q(k: QuantityKind) -> QuantitySpec { QuantitySpec::semantic(SemanticType::new(k, ValueForm::Static)) }
fn channels(scale: f64, formation: [f64; 3]) -> Vec<GasReaction> {
    let components = ["A", "B", "C"].into_iter().zip(formation)
        .map(|(id, f)| source::component(id, 0.028, f)).collect::<Vec<_>>();
    let species = components.iter().map(|(s, _)| s.species().clone()).collect::<Vec<_>>();
    let ids = ["ab", "ac"].map(|id| ReactionId::new(id).unwrap());
    let matrix = StoichiometricMatrix::new(species.clone(), ids.to_vec(),
        vec![vec![-1,-1],vec![1,0],vec![0,1]]).unwrap();
    let elements = ElementalMatrix::new(vec![ElementId::new("X").unwrap()],species.clone(),vec![vec![1;3]]).unwrap();
    let charge = ChargeVector::new(species,vec![0;3]).unwrap();
    let certificate = verify_conservation(&elements,&matrix,&charge).unwrap();
    let catalog = ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(
        components.iter().map(|(s,_)|s.clone()).collect(),
        Composition::new(CompositionBasis::MoleFraction,vec![0.5,0.25,0.25]).unwrap()).unwrap(),
        components.iter().map(|(_,t)|t.clone()).collect(),None).unwrap();
    ids.into_iter().enumerate().map(|(i,id)| {
        let equilibrium = IdealGasReactionEquilibriumV1::new(matrix.clone(),certificate,id,
            vec![components[0].0.clone(),components[i+1].0.clone()]).unwrap();
        GasReaction::new(catalog.clone(),DeclaredStoichiometricMassActionV1::new(
            ThermodynamicReverseRateClosureV1::new(equilibrium,
                ReactionProgressRateScaleV1::new(ReactionProgressRateQuantityV1::new(scale)).unwrap()),
            MassActionKineticOrderConventionV1::CallerDeclaredStoichiometricCoefficients)).unwrap()
    }).collect()
}
fn gas(law: &GasReaction, fractions: [f64;3], temperature: f64) -> ClosedGasParcel {
    let model = law.catalog().with_composition(Composition::new(CompositionBasis::MoleFraction,fractions.to_vec()).unwrap()).unwrap();
    let point = QueryPoint::new().with_quantity("temperature",q(QuantityKind::AbsoluteTemperature),temperature).unwrap()
        .with_quantity("pressure",q(QuantityKind::Pressure),100000.0).unwrap().with("catalyst-lot",7.0).unwrap();
    ClosedGasParcel::new(model,point,0.001).unwrap()
}
fn controls() -> GasReactionSolve {
    GasReactionSolve { gas: fs_material::thermochemical_gas::energy::GasEnergySolve {
        minimum_temperature_k:200.0,maximum_temperature_k:1800.0,energy_tolerance_j_kg:1e-7,maximum_iterations:128 },
        maximum_consumed_fraction:0.2,species_tolerance_mol:1e-12,entropy_tolerance_j_k:1e-10 }
}
fn policy() -> ReactionSubcycling {
    ReactionSubcycling {absolute_amount_mol:1e-8,absolute_temperature_k:0.01,absolute_pressure_pa:1.0,relative_tolerance:1e-3,
        maximum_substeps:256,maximum_evaluations:511}
}
fn rate_source(scale: f64) -> MaterialRateScale {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key:PropertyKey::new(ACTIVITY_RATE_SCALE_PROPERTY,ACTIVITY_RATE_SCALE_DIMS),
        value:PropertyValue::Curve {abscissa:"temperature".into(),abscissa_dims:Temperature::DIMS,
            knots:vec![(200.0,0.5*scale),(800.0,2.0*scale)],dims:ACTIVITY_RATE_SCALE_DIMS},
        validity:ValidityDomain::unconstrained().with_quantity("temperature",q(QuantityKind::AbsoluteTemperature),200.0,800.0)
            .with_quantity("pressure",q(QuantityKind::Pressure),1000.0,1e7).with("catalyst-lot",7.0,7.0),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,observations:vec![],
        provenance:Provenance {source:"synthetic subcycling rate".into(),license:"CC0-1.0".into(),artifact:None},
    }).unwrap();
    MaterialRateScale::new(MaterialCard::assemble(MaterialStateId {chemistry:"synthetic".into(),phase:"gas".into(),
        process:"test".into(),revision:0},claims,vec![]).unwrap(),MaterialPropertySelection::SingleClaimOnly,
        Pressure::new(100000.0)).unwrap()
}
fn near(a:f64,b:f64,tolerance:f64) {assert!((a-b).abs()<=tolerance,"{a} != {b}");}

#[test]
fn g1_fast_reversible_material_reaction_refines_instead_of_clipping_inventory() {
    let raw = channels(1000.0,[0.0;3]).remove(0);
    let old = gas(&raw,[1.0,0.0,0.0],400.0);
    assert!(matches!(raw.parcel_candidate(&old,0.01,controls(),||false),Err(GasReactionError::ConsumptionLimit {..})));
    let refined = raw.clone().with_subcycling(policy()).unwrap();
    let result = refined.parcel_candidate(&old,0.01,controls(),||false).unwrap();
    let report = result.step().subcycling.as_ref().unwrap();
    assert!(report.steps.len()>2 && report.steps.len()<=policy().maximum_substeps);
    assert!(report.estimated_error_ratio<=1.0 && report.evaluations<=policy().maximum_evaluations);
    let amount = old.mass_kg()/0.028;
    let lambda = 1000.0*UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K*400.0/100000.0;
    let exact = amount*(0.5+0.5*(-2.0*lambda*0.01).exp());
    near(result.step().species[0].final_mol,exact,8e-6);
    near(result.gas().state().state().temperature,400.0,1e-7);
    assert_eq!(result.gas().mass_kg(),old.mass_kg());
    assert_eq!(result.gas().volume_m3(),old.volume_m3());
    assert_eq!(report.steps.first().unwrap().start_s,0.0);
    assert_eq!(report.steps.last().unwrap().end_s,0.01);
    assert!(report.steps.windows(2).all(|p|p[0].end_s==p[1].start_s));
    assert!(report.steps.iter().all(|s|s.reaction.maximum_consumed_fraction<=controls().maximum_consumed_fraction));
    assert!(result.step().maximum_species_residual_mol<=controls().species_tolerance_mol);
    assert!(result.step().stored_energy_change_j.abs()<=old.mass_kg()*controls().gas.energy_tolerance_j_kg);
    assert!(raw.subcycling().is_none());
    assert!(raw.parcel_candidate(&old,0.0001,controls(),||false).unwrap().step().subcycling.is_none());
}

#[test]
fn g3_material_rates_refresh_on_each_actual_substep_and_replay_exactly() {
    let raw = GasReaction::combine(channels(1.0,[0.0,-3000.0,-6000.0]).into_iter().enumerate()
        .map(|(i,c)|c.with_material_rate(rate_source(500.0*(i+1) as f64)).unwrap()).collect()).unwrap();
    let old = gas(&raw,[1.0,0.0,0.0],400.0);
    let refined = raw.clone().with_subcycling(ReactionSubcycling {relative_tolerance:0.01,..policy()}).unwrap();
    let result = refined.parcel_candidate(&old,0.001,controls(),||false).unwrap();
    let history = &result.step().subcycling.as_ref().unwrap().steps;
    let mut replay = old.clone();
    let mut inner = controls();
    inner.gas.energy_tolerance_j_kg/=history.len() as f64;
    inner.species_tolerance_mol/=history.len() as f64;
    inner.entropy_tolerance_j_k/=history.len() as f64;
    let mut extents=[0.0;2];
    for interval in history {
        let expected=raw.parcel_candidate(&replay,interval.end_s-interval.start_s,inner,||false).unwrap();
        assert_eq!(&interval.reaction,expected.step());
        let scale=interval.reaction.rate_material.as_ref().unwrap().property(ACTIVITY_RATE_SCALE_PROPERTY).unwrap().value_si();
        near(scale,500.0*replay.state().state().temperature/400.0,1e-9);
        for (i,(_,extent)) in interval.reaction.reaction_extents().enumerate() {extents[i]+=extent;}
        replay=expected.into_parts().0;
    }
    assert_eq!(&replay,result.gas());
    for ((_,actual),expected) in result.step().reaction_extents().zip(extents) {near(actual,expected,1e-14);}
    assert!(replay.state().state().temperature>400.0);
    assert!(history.last().unwrap().reaction.rate_material.as_ref().unwrap().property(ACTIVITY_RATE_SCALE_PROPERTY).unwrap().value_si()>500.0);
    assert!(result.step().entropy_change_j_k>0.0);
    assert_eq!(replay.point().axes().get("catalyst-lot"),Some(&7.0));
}

#[test]
fn g3_tighter_error_policy_refines_the_same_conserved_trajectory() {
    let raw=channels(1000.0,[0.0;3]).remove(0);
    let old=gas(&raw,[1.0,0.0,0.0],400.0);
    let solve=|relative_tolerance|raw.clone().with_subcycling(ReactionSubcycling {relative_tolerance,..policy()}).unwrap()
        .parcel_candidate(&old,0.01,controls(),||false).unwrap();
    let loose=solve(0.004);let tight=solve(0.001);
    let lambda=1000.0*UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K*400.0/100000.0;
    let exact=(old.mass_kg()/0.028)*(0.5+0.5*(-2.0*lambda*0.01).exp());
    assert!(tight.step().subcycling.as_ref().unwrap().steps.len()>loose.step().subcycling.as_ref().unwrap().steps.len());
    assert!((tight.step().species[0].final_mol-exact).abs()<(loose.step().species[0].final_mol-exact).abs());
}

#[test]
fn g4_budget_source_and_late_cancellation_return_no_partial_trajectory() {
    let raw=channels(1000.0,[0.0;3]).remove(0);
    let old=gas(&raw,[1.0,0.0,0.0],400.0);let saved=old.clone();
    let limited=raw.clone().with_subcycling(ReactionSubcycling {maximum_evaluations:3,..policy()}).unwrap();
    assert!(matches!(limited.parcel_candidate(&old,0.01,controls(),||false),
        Err(GasReactionError::SubcyclingLimit {evaluations:3,..})));
    let refined=raw.clone().with_subcycling(ReactionSubcycling {relative_tolerance:0.01,..policy()}).unwrap();
    let mut calls=0;
    refined.parcel_candidate(&old,0.0001,controls(),||{calls+=1;false}).unwrap();
    for stop in [1,calls/2,calls] {
        let mut seen=0;
        assert_eq!(refined.parcel_candidate(&old,0.0001,controls(),||{seen+=1;seen==stop}),Err(GasReactionError::Cancelled));
        assert_eq!(old,saved);
    }
    let sourced=raw.with_material_rate(rate_source(1000.0)).unwrap().with_subcycling(policy()).unwrap();
    let hot=gas(&sourced,[1.0,0.0,0.0],900.0);
    assert!(matches!(sourced.parcel_candidate(&hot,0.0001,controls(),||false),Err(GasReactionError::RateMaterial(_))));
    assert_eq!(old,saved);
}

#[test]
fn g3_existing_piston_consumes_one_complete_refined_chemical_interval() {
    let raw=GasReaction::combine(channels(500.0,[0.0,-3000.0,-6000.0])).unwrap();
    let refined=raw.with_subcycling(ReactionSubcycling {relative_tolerance:0.01,..policy()}).unwrap();
    let old=gas(&refined,[1.0,0.0,0.0],400.0);
    let piston=GasPiston::new(old.clone(),1.0,0.01,0.0).unwrap();
    let ch=source::chamber(10,&piston,0.001);
    let chemistry=refined.parcel_candidate(&old,0.001,controls(),||false).unwrap();
    let mut expected=GasPiston::new(chemistry.gas().clone(),1.0,0.01,0.0).unwrap();
    expected.step(ch.input,ch.piston_solve,||false).unwrap();
    let actual=piston.reactive_candidate(&refined,ch.input,controls(),ch.piston_solve,||false).unwrap();
    assert_eq!(&actual.step().reaction,chemistry.step());
    assert_eq!(actual.piston().gas(),expected.gas());
    assert_eq!(actual.piston().velocity_m_s(),expected.velocity_m_s());
    assert_eq!(actual.piston().last_window_duration_s(),Some(ch.input.duration_s));
    assert!(actual.step().total_energy_residual_j.abs()<=ch.piston_solve.energy_tolerance_j);
    assert_eq!(piston.gas(),&old);
}

#[test]
fn g0_subcycling_is_one_explicit_mechanism_policy_and_equilibrium_stays_exact() {
    let channels=channels(1000.0,[0.0;3]);
    let bound=channels[0].clone().with_subcycling(policy()).unwrap();
    assert!(bound.clone().with_subcycling(policy()).is_err());
    assert_eq!(GasReaction::combine(vec![bound.clone()]).unwrap(),bound);
    assert!(GasReaction::combine(vec![bound.clone(),channels[1].clone()]).is_err());
    for maximum_substeps in [0,1,3,512] {
        assert!(channels[0].clone().with_subcycling(ReactionSubcycling {maximum_substeps,..policy()}).is_err());
    }
    let old=gas(&bound,[0.5,0.5,0.0],400.0);
    let result=bound.parcel_candidate(&old,0.01,controls(),||false).unwrap();
    assert_eq!(result.gas(),&old);
    assert!(!result.step().has_progress());
    assert_eq!(result.step().subcycling.as_ref().unwrap().steps.len(),2);
}
