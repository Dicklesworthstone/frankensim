//! Actual source-backed kinetics, caloric reconstruction and existing mechanics.
//! Manufactured isomers, NOT a measured combustion mechanism.
#[path = "material_exchange/gas.rs"]
#[allow(dead_code)]
mod source;
use fs_couple::cavity_phs::piston::{GasPiston, reactive::{GasReaction, GasReactionError, GasReactionSolve}};
use fs_couple::cavity_phs::piston::reactive::material_rate::{MaterialRateScale, ACTIVITY_RATE_SCALE_PROPERTY, ACTIVITY_RATE_SCALE_DIMS};
use fs_material::thermochemical_gas::{ThermochemicalGasModel, energy::parcel::ClosedGasParcel};
use fs_material::state_point::MaterialPropertySelection;
use fs_matdb::*;
use fs_evidence::ValidityDomain;
use fs_qty::{Pressure, Temperature, QuantitySpec, semantic::{QuantityKind, SemanticType, ValueForm}};
use fs_thermochem::*;
use fs_thermochem::equilibrium::*;
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use std::collections::BTreeMap;

const DT: f64 = 0.00005;

#[path = "reactive_mechanism/implicit.rs"]
mod implicit;
fn quantity(k: QuantityKind) -> QuantitySpec { QuantitySpec::semantic(SemanticType::new(k, ValueForm::Static)) }
fn laws(edges: &[(&str, usize, usize, f64)], formation: [f64; 3]) -> Vec<GasReaction> {
    let components = ["A", "B", "C"].into_iter().zip(formation)
        .map(|(name, f)| source::component(name, 0.028, f)).collect::<Vec<_>>();
    let species = components.iter().map(|(s, _)| s.species().clone()).collect::<Vec<_>>();
    let reactions = edges.iter().map(|e| ReactionId::new(e.0).unwrap()).collect();
    let rows = (0..3).map(|i| edges.iter().map(|e|
        if i == e.1 { -1 } else if i == e.2 { 1 } else { 0 }).collect()).collect();
    let matrix = StoichiometricMatrix::new(species.clone(), reactions, rows).unwrap();
    let elements = ElementalMatrix::new(vec![ElementId::new("X").unwrap()], species.clone(), vec![vec![1; 3]]).unwrap();
    let charge = ChargeVector::new(species, vec![0; 3]).unwrap();
    let conserved = verify_conservation(&elements, &matrix, &charge).unwrap();
    let catalog = ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(
        components.iter().map(|(s, _)| s.clone()).collect(),
        Composition::new(CompositionBasis::MoleFraction, vec![0.25, 0.25, 0.5]).unwrap()).unwrap(),
        components.iter().map(|(_, t)| t.clone()).collect(), None).unwrap();
    edges.iter().map(|&(name, a, b, scale)| {
        let equilibrium = IdealGasReactionEquilibriumV1::new(matrix.clone(), conserved,
            ReactionId::new(name).unwrap(), vec![components[a].0.clone(), components[b].0.clone()]).unwrap();
        GasReaction::new(catalog.clone(), DeclaredStoichiometricMassActionV1::new(
            ThermodynamicReverseRateClosureV1::new(equilibrium,
                ReactionProgressRateScaleV1::new(ReactionProgressRateQuantityV1::new(scale)).unwrap()),
            MassActionKineticOrderConventionV1::CallerDeclaredStoichiometricCoefficients)).unwrap()
    }).collect()
}
fn gas(law: &GasReaction, x: [f64; 3], temperature: f64) -> ClosedGasParcel {
    let model = law.catalog().with_composition(Composition::new(CompositionBasis::MoleFraction, x.to_vec()).unwrap()).unwrap();
    let point = QueryPoint::new().with_quantity("temperature", quantity(QuantityKind::AbsoluteTemperature), temperature).unwrap()
        .with_quantity("pressure", quantity(QuantityKind::Pressure), 100000.0).unwrap()
        .with("catalyst-lot", 7.0).unwrap();
    ClosedGasParcel::new(model, point, 0.001).unwrap()
}
fn controls() -> GasReactionSolve {
    GasReactionSolve { gas: fs_material::thermochemical_gas::energy::GasEnergySolve {
        minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 },
        maximum_consumed_fraction: 0.2, species_tolerance_mol: 1e-12, entropy_tolerance_j_k: 1e-10 }
}
fn near(a: f64, b: f64, abs: f64) { assert!((a-b).abs() <= abs, "{a} != {b}"); }

#[test]
fn shared_inventory_and_formation_energy_match_the_independent_isomer_balance() {
    let laws = laws(&[("ab", 0, 1, 100.0), ("ac", 0, 2, 200.0)], [0.0, -3000.0, -6000.0]);
    let old = gas(&laws[0], [1.0, 0.0, 0.0], 400.0);
    let one = GasReaction::combine(vec![laws[0].clone()]).unwrap();
    assert_eq!(one, laws[0]);
    assert_eq!(one.parcel_candidate(&old, DT, controls(), || false).unwrap(),
        laws[0].parcel_candidate(&old, DT, controls(), || false).unwrap());
    let mechanism = GasReaction::combine(laws.clone()).unwrap();
    let reversed = GasReaction::combine(laws.into_iter().rev().collect()).unwrap();
    assert_eq!(mechanism, reversed);
    let trial = mechanism.parcel_candidate(&old, DT, controls(), || false).unwrap();
    let report = trial.step();
    let amount = old.mass_kg()/0.028;
    let a = 100.0*0.001*DT;
    let b = 200.0*0.001*DT;
    near(report.extent_mol, a, 1e-16);
    near(report.additional_reactions[0].extent_mol, b, 1e-16);
    for (row, expected) in report.species.iter().zip([amount-a-b, a, b]) { near(row.final_mol, expected, 1e-12); }
    let expected_t = 400.0+(3000.0*a+6000.0*b)/(amount*2.5*UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K);
    near(trial.gas().state().state().temperature, expected_t, 1e-7);
    near(trial.gas().state().state().pressure, 100000.0*expected_t/400.0, 1e-5);
    assert_eq!(trial.gas().mass_kg(), old.mass_kg());
    assert_eq!(trial.gas().volume_m3(), old.volume_m3());
    assert_eq!(trial.gas().point().axes().get("catalyst-lot"), Some(&7.0));
    assert!(trial.gas().state().state().sound_speed > old.state().state().sound_speed);
    assert!(report.has_progress());
    assert!(report.stored_energy_change_j.abs() <= old.mass_kg()*controls().gas.energy_tolerance_j_kg);
}

#[test]
fn newly_created_intermediate_is_not_consumed_until_the_next_physical_interval() {
    let mechanism = GasReaction::combine(laws(&[("ab", 0, 1, 200.0), ("bc", 1, 2, 300.0)], [0.0; 3])).unwrap();
    let old = gas(&mechanism, [1.0, 0.0, 0.0], 400.0);
    let first = mechanism.parcel_candidate(&old, DT, controls(), || false).unwrap();
    assert!(first.step().species[1].final_mol > 0.0);
    assert_eq!(first.step().additional_reactions[0].extent_mol, 0.0);
    assert_eq!(first.step().species[2].final_mol, 0.0);
    let second = mechanism.parcel_candidate(first.gas(), DT, controls(), || false).unwrap();
    assert!(second.step().additional_reactions[0].extent_mol > 0.0);
    assert!(second.step().species[2].final_mol > 0.0);
    near(second.gas().state().state().temperature, 400.0, 1e-7);
}

#[test]
fn individually_admissible_channels_cannot_spend_the_same_reactant_twice() {
    let initial = laws(&[("ab", 0, 1, 1.0), ("ac", 0, 2, 1.0)], [0.0; 3]);
    let old = gas(&initial[0], [1.0, 0.0, 0.0], 400.0);
    let k = 0.1*(old.mass_kg()/0.028)/(DT*old.volume_m3());
    let channels = laws(&[("ab", 0, 1, k), ("ac", 0, 2, k)], [0.0; 3]);
    let limit = GasReactionSolve { maximum_consumed_fraction: 0.15, ..controls() };
    for channel in &channels { channel.parcel_candidate(&old, DT, limit, || false).unwrap(); }
    let mechanism = GasReaction::combine(channels).unwrap();
    let saved = old.clone();
    match mechanism.parcel_candidate(&old, DT, limit, || false) {
        Err(GasReactionError::ConsumptionLimit { species, fraction, .. }) => {
            assert_eq!(species, SpeciesId::new("A").unwrap()); near(fraction, 0.2, 1e-12);
        },
        other => panic!("aggregate consumption must refuse: {other:?}"),
    }
    assert_eq!(old, saved);
    // A third channel produces A, but cannot conceal the two original withdrawals.
    let amount = old.mass_kg()/0.028;
    let base = 0.1*amount/(DT*old.volume_m3());
    let channels = laws(&[("ab1", 0, 1, base), ("ab2", 0, 1, base), ("ca", 2, 0, 2.0*base)], [0.0; 3]);
    let mixed = gas(&channels[0], [0.4, 0.0, 0.6], 400.0);
    for channel in &channels { channel.parcel_candidate(&mixed, DT, limit, || false).unwrap(); }
    assert!(matches!(GasReaction::combine(channels).unwrap().parcel_candidate(&mixed, DT, limit, || false),
        Err(GasReactionError::ConsumptionLimit { .. })));
}

#[test]
fn invalid_mechanisms_and_cancelled_or_unsupported_endpoints_publish_nothing() {
    let channels = laws(&[("ab", 0, 1, 100.0), ("ac", 0, 2, 200.0)], [0.0, -3000.0, -6000.0]);
    assert!(GasReaction::combine(vec![]).is_err());
    assert!(GasReaction::combine(vec![channels[0].clone(); 33]).is_err());
    assert!(GasReaction::combine(vec![channels[0].clone(), channels[0].clone()]).is_err());
    let conflict = laws(&[("ab", 0, 1, 100.0), ("ac", 0, 2, 200.0)], [0.0, -3000.0, -7000.0]);
    assert!(GasReaction::combine(vec![channels[0].clone(), conflict[1].clone()]).is_err());
    let unrelated_matrix = laws(&[("different", 0, 2, 200.0)], [0.0, -3000.0, -6000.0]);
    assert!(GasReaction::combine(vec![channels[0].clone(), unrelated_matrix[0].clone()]).is_err());
    let mechanism = GasReaction::combine(channels).unwrap();
    let old = gas(&mechanism, [1.0, 0.0, 0.0], 400.0);
    let saved = old.clone();
    let mut total = 0;
    mechanism.parcel_candidate(&old, DT, controls(), || { total += 1; false }).unwrap();
    for stop in [1, total/2, total] {
        let mut polls = 0;
        assert_eq!(mechanism.parcel_candidate(&old, DT, controls(), || { polls += 1; polls == stop }),
            Err(GasReactionError::Cancelled));
        assert_eq!(old, saved);
    }
    let narrow = GasReactionSolve { gas: fs_material::thermochemical_gas::energy::GasEnergySolve {
        maximum_temperature_k: 400.001, ..controls().gas }, ..controls() };
    assert!(mechanism.parcel_candidate(&old, DT, narrow, || false).is_err());
    assert_eq!(old, saved);
}

fn rate_source(scale_at_400: f64) -> MaterialRateScale {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(ACTIVITY_RATE_SCALE_PROPERTY, ACTIVITY_RATE_SCALE_DIMS),
        value: PropertyValue::Curve { abscissa: "temperature".into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0, 0.5*scale_at_400), (800.0, 2.0*scale_at_400)], dims: ACTIVITY_RATE_SCALE_DIMS },
        validity: ValidityDomain::unconstrained().with_quantity("temperature", quantity(QuantityKind::AbsoluteTemperature), 200.0, 800.0)
            .with_quantity("pressure", quantity(QuantityKind::Pressure), 1000.0, 1e7).with("catalyst-lot", 7.0, 7.0),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: Provenance { source: "synthetic per-channel kinetics".into(), license: "CC0-1.0".into(), artifact: None },
    }).unwrap();
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "synthetic rate".into(), phase: "gas".into(), process: "test".into(), revision: 0 }, claims, vec![]).unwrap();
    MaterialRateScale::new(card, MaterialPropertySelection::SingleClaimOnly, Pressure::new(100000.0)).unwrap()
}

#[test]
fn each_channel_uses_its_own_live_material_rate_without_changing_the_reverse_closure() {
    let channels = laws(&[("ab", 0, 1, 1.0), ("ac", 0, 2, 1.0)], [0.0, -3000.0, -6000.0]);
    let mechanism = GasReaction::combine(channels.into_iter().zip([100.0, 300.0])
        .map(|(c, k)| c.with_material_rate(rate_source(k)).unwrap()).collect()).unwrap();
    assert!(mechanism.clone().with_material_rate(rate_source(1.0)).is_err());
    for (t, multiplier) in [(400.0, 1.0), (600.0, 1.5)] {
        let old = gas(&mechanism, [0.5, 0.25, 0.25], t);
        let report = mechanism.parcel_candidate(&old, DT, controls(), || false).unwrap();
        let reference = GasReaction::combine(laws(&[("ab", 0, 1, 100.0*multiplier), ("ac", 0, 2, 300.0*multiplier)],
            [0.0, -3000.0, -6000.0])).unwrap().parcel_candidate(&old, DT, controls(), || false).unwrap();
        // Exact rate equality exercises the original reverse closure too.
        assert_eq!(report.step().rate, reference.step().rate);
        assert_eq!(report.step().additional_reactions[0].rate, reference.step().additional_reactions[0].rate);
        for (sample, k) in [(report.step().rate_material.as_ref().unwrap(), 100.0),
            (report.step().additional_reactions[0].rate_material.as_ref().unwrap(), 300.0)] {
            near(sample.property(ACTIVITY_RATE_SCALE_PROPERTY).unwrap().value_si(), multiplier*k, 1e-10);
            assert!(sample.query_point().contains(&("catalyst-lot".into(), 7.0)));
        }
        assert_eq!(report.gas(), reference.gas());
    }
    let unsupported = gas(&mechanism, [1.0, 0.0, 0.0], 900.0);
    assert!(matches!(mechanism.parcel_candidate(&unsupported, DT, controls(), || false), Err(GasReactionError::RateMaterial(_))));
}

#[test]
fn secondary_reaction_reaches_existing_piston_and_closed_network_energy_accounts() {
    use fs_couple::cavity_phs::piston::exchange::graph::{reactive, MovingGasDrive, MovingGasSolve};
    use fs_material::thermochemical_gas::energy::parcel::transfer::snapshot::GasSnapshotSolve;
    use fs_ports::{ModuleWindowAudit, PortedModule};
    // Canonical first reaction is idle: B and C are both absent initially.
    let mechanism = GasReaction::combine(laws(&[("idle", 1, 2, 100.0), ("working", 0, 1, 200.0)],
        [0.0, -3000.0, -6000.0])).unwrap();
    let piston = GasPiston::new(gas(&mechanism, [1.0, 0.0, 0.0], 400.0), 1.0, 0.01, 0.0).unwrap();
    let ch = source::chamber(10, &piston, DT);
    let standalone = piston.reactive_candidate(&mechanism, ch.input, controls(), ch.piston_solve, || false).unwrap();
    assert_eq!(standalone.step().reaction.extent_mol, 0.0);
    assert!(standalone.step().reaction.has_progress());
    let pistons = BTreeMap::from([(10, piston.clone())]);
    let trial = reactive::step_candidate(&pistons, &BTreeMap::new(), &[], DT,
        &BTreeMap::from([(10, controls().gas)]), &BTreeMap::from([(10, MovingGasDrive { input: ch.input, solve: ch.piston_solve })]),
        &[reactive::ChamberReaction { chamber: 10, reaction: &mechanism, controls: controls() }],
        MovingGasSolve { exchange: GasSnapshotSolve { maximum_withdrawn_fraction: 0.1, mass_tolerance_kg: 1e-13,
            energy_tolerance_j: 1e-7 }, total_energy_tolerance_j: 1e-6 }, || false).unwrap();
    assert_eq!(trial.pistons()[&10].gas(), standalone.piston().gas());
    near(trial.pistons()[&10].velocity_m_s(), standalone.piston().velocity_m_s(), 1e-12);
    assert_eq!(&trial.reactions()[&10], &standalone.step().reaction);
    assert!(trial.step().total_energy_residual_j.abs() < 1e-6);
    assert_eq!(trial.pistons()[&10].last_window_duration_s(), Some(DT));
    let ModuleWindowAudit::Available(window) = trial.pistons()[&10].audit_last_window().unwrap() else {
        panic!("the closed reacting graph must retain a complete account");
    };
    assert!(window.violations().is_empty());
    assert_eq!(pistons[&10], piston);
}

#[test]
fn coupled_nonisomeric_stoichiometry_changes_molecule_count_without_changing_mass() {
    let components = [source::component("A", 0.014, 0.0), source::component("B", 0.028, -60000.0),
        source::component("C", 0.014, 0.0)];
    let species = components.iter().map(|(s, _)| s.species().clone()).collect::<Vec<_>>();
    let ids = [ReactionId::new("dimerize").unwrap(), ReactionId::new("split").unwrap()];
    let matrix = StoichiometricMatrix::new(species.clone(), ids.to_vec(),
        vec![vec![-2, 0], vec![1, -1], vec![0, 2]]).unwrap();
    let conserved = verify_conservation(&ElementalMatrix::new(vec![ElementId::new("X").unwrap()],
        species.clone(), vec![vec![1, 2, 1]]).unwrap(), &matrix, &ChargeVector::new(species, vec![0; 3]).unwrap()).unwrap();
    let catalog = ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(
        components.iter().map(|(s, _)| s.clone()).collect(),
        Composition::new(CompositionBasis::MoleFraction, vec![0.25, 0.25, 0.5]).unwrap()).unwrap(),
        components.iter().map(|(_, t)| t.clone()).collect(), None).unwrap();
    let channels = [(0, [0, 1], 30.0), (1, [1, 2], 50.0)].map(|(j, terms, k)| {
        let equilibrium = IdealGasReactionEquilibriumV1::new(matrix.clone(), conserved, ids[j].clone(),
            terms.map(|i| components[i].0.clone()).to_vec()).unwrap();
        GasReaction::new(catalog.clone(), DeclaredStoichiometricMassActionV1::new(
            ThermodynamicReverseRateClosureV1::new(equilibrium,
                ReactionProgressRateScaleV1::new(ReactionProgressRateQuantityV1::new(k)).unwrap()),
            MassActionKineticOrderConventionV1::CallerDeclaredStoichiometricCoefficients)).unwrap()
    });
    let mechanism = GasReaction::combine(channels.to_vec()).unwrap();
    let old = gas(&mechanism, [1.0, 0.0, 0.0], 400.0);
    let trial = mechanism.parcel_candidate(&old, DT, controls(), || false).unwrap();
    let n = old.mass_kg()/0.014;
    let xi = trial.step().extent_mol;
    near(xi, 30.0*0.001*DT, 1e-16);
    near(trial.step().species[0].final_mol, n-2.0*xi, 1e-12);
    near(trial.step().species[1].final_mol, xi, 1e-12);
    assert_eq!(trial.step().species[2].final_mol, 0.0);
    let t = (n*2.5*UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K*400.0+60000.0*xi)
        / ((n-xi)*2.5*UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K);
    near(trial.gas().state().state().temperature, t, 1e-7);
    assert_eq!(trial.gas().mass_kg(), old.mass_kg());
}

