//! Actual material-card diffusivity and original thermochemical endpoint paths.
use super::*;
use crate::thermochemical_gas::tests::{model, point, species, transport};
use crate::thermochemical_gas::energy::GasEnergySolve;
use crate::thermochemical_gas::ThermochemicalGasModel;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, UncertaintyModel};
use fs_qty::{Temperature, QuantitySpec, semantic::{QuantityKind, SemanticType, ValueForm}};
use fs_thermochem::{Composition, CompositionBasis};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;

fn controls() -> GasDiffusionSolve {
    let energy = GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 };
    GasDiffusionSolve { maximum_temperature_difference_k: 1e-5, maximum_pressure_difference_pa: 0.01,
        gradient_tolerance_per_m: 1e-9, molar_flux_tolerance_mol_m2_s: 1e-12, entropy_tolerance_j_k: 1e-10,
        exchange: GasSpeciesExchangeSolve { gas: [energy; 2], maximum_species_fraction: 0.2,
            mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8 } }
}
fn geometry() -> GasDiffusionGeometry { GasDiffusionGeometry { area_m2: 0.01, distance_m: 0.01 } }
fn gas(x: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(model(vec![x, 1.0-x], 0.001), point(400.0, 100000.0), 0.001).unwrap()
}
fn card(value: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(BINARY_DIFFUSIVITY_PROPERTY, BINARY_DIFFUSIVITY_DIMS),
        value: PropertyValue::Curve { abscissa: "temperature".into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0, value), (1800.0, value)], dims: BINARY_DIFFUSIVITY_DIMS },
        validity: ValidityDomain::unconstrained()
            .with_quantity("temperature", QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature, ValueForm::Static)), 200.0, 1800.0)
            .with_quantity("pressure", QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)), 1000.0, 1e7),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: Provenance { source: "synthetic binary diffusion".into(),
            license: "CC0-1.0".into(), artifact: None },
    }).unwrap();
    MaterialCard::assemble(MaterialStateId { chemistry: "synthetic pair".into(), phase: "gas".into(),
        process: "diffusion fixture".into(), revision: 0 }, claims, vec![]).unwrap()
}
fn pair<'a>(a: &str, b: &str, card: &'a MaterialCard, context: &'a QueryPoint) -> BinaryDiffusionMaterial<'a> {
    BinaryDiffusionMaterial { species: [SpeciesId::new(a).unwrap(), SpeciesId::new(b).unwrap()],
        material: card, context, selection: MaterialPropertySelection::SingleClaimOnly }
}
fn near(a: f64, b: f64, tolerance: f64) { assert!((a-b).abs() <= tolerance, "{a} != {b}"); }

#[test]
fn binary_limit_is_fick_and_is_molar_not_mass_average_counterflow() {
    let (a, b) = (gas(0.75), gas(0.25)); let d = card(2e-5); let p = point(300.0, 90000.0);
    let flux = a.diffusion_flux_to(&b, geometry(), &[pair("A", "B", &d, &p)], controls(), || false).unwrap();
    let expected = 100000.0/(UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K*400.0)*2e-5*0.5/0.01;
    near(flux.species[0].molar_flux_mol_m2_s, expected, 1e-13);
    near(flux.species[1].molar_flux_mol_m2_s, -expected, 1e-13);
    assert!((flux.species[0].mass_rate_left_to_right_kg_s
        + flux.species[1].mass_rate_left_to_right_kg_s).abs() > 1e-6);
    assert_eq!(flux.pairs[0].diffusivity_m2_s, 2e-5);
    assert!(flux.pairs[0].material.query_point().contains(&("temperature".into(), 400.0)));
    let changed = card(4e-5);
    let faster = a.diffusion_flux_to(&b, geometry(), &[pair("B", "A", &changed, &p)], controls(), || false).unwrap();
    near(faster.species[0].molar_flux_mol_m2_s, 2.0*expected, 1e-13);
    let reverse = b.diffusion_flux_to(&a, geometry(), &[pair("B", "A", &d, &p)], controls(), || false).unwrap();
    near(reverse.species[0].molar_flux_mol_m2_s, -expected, 1e-13);
}

fn ternary(x: Vec<f64>) -> ClosedGasParcel {
    let model = ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(
        vec![species("A", 0.028, 3.5, 0.001), species("B", 0.004, 2.5, 0.0), species("C", 0.044, 4.0, 0.0)],
        Composition::new(CompositionBasis::MoleFraction, x).unwrap()).unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2e-5, 0.15), transport("C", 1.5e-5, 0.02)], None).unwrap();
    ClosedGasParcel::new(model, point(400.0, 100000.0), 0.001).unwrap()
}
#[test]
fn ternary_cross_diffusion_moves_a_species_with_zero_own_gradient_and_ignores_pair_order() {
    let (a, b) = (ternary(vec![0.6, 0.2, 0.2]), ternary(vec![0.2, 0.6, 0.2]));
    let cards = [card(1e-5), card(2e-5), card(5e-5)]; let context = point(400.0, 100000.0);
    let pairs = [pair("A", "B", &cards[0], &context), pair("A", "C", &cards[1], &context),
        pair("B", "C", &cards[2], &context)];
    let flux = a.diffusion_flux_to(&b, geometry(), &pairs, controls(), || false).unwrap();
    let reordered = [pairs[2].clone(), pairs[0].clone(), pairs[1].clone()];
    assert_eq!(flux, a.diffusion_flux_to(&b, geometry(), &reordered, controls(), || false).unwrap());
    // Independent two-unknown elimination of J_C = -J_A-J_B.
    let c = 100000.0/(UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K*400.0);
    let aa = 0.4/1e-5+0.6/2e-5; let ab = -0.4/1e-5+0.4/2e-5;
    let ba = -0.4/1e-5+0.4/5e-5; let bb = 0.4/1e-5+0.6/5e-5;
    let r = c*40.0; let det = aa*bb-ab*ba;
    let ja = (r*bb+ab*r)/det; let jb = (-aa*r-ba*r)/det;
    near(flux.species[0].molar_flux_mol_m2_s, ja, 1e-13);
    near(flux.species[1].molar_flux_mol_m2_s, jb, 1e-13);
    near(flux.species[2].molar_flux_mol_m2_s, -ja-jb, 1e-13);
    assert!(flux.species[2].molar_flux_mol_m2_s.abs() > 1e-4);
}

#[test]
fn equal_pressure_mixing_changes_shared_sound_and_transport_without_creating_heat() {
    let (a, b) = (gas(0.75), gas(0.25)); let d = card(2e-5); let p = point(400.0, 100000.0);
    let trial = a.diffuse_with(&b, geometry(), 0.1, &[pair("A", "B", &d, &p)], controls(), || false).unwrap();
    assert!(trial.entropy_change_j_k() > 0.0);
    for (old, next) in [(&a, trial.exchange().left()), (&b, trial.exchange().right())] {
        near(next.state().state().temperature, 400.0, 1e-7);
        near(next.state().state().pressure, 100000.0, 1e-5);
        assert_ne!(next.state().state().sound_speed, old.state().state().sound_speed);
        assert_ne!(next.state().state().thermal_conductivity, old.state().state().thermal_conductivity);
    }
    assert!(trial.exchange().step().energy_residual_j.abs() < 1e-8);
    assert!(trial.exchange().step().maximum_species_mass_residual_kg < 1e-14);
}

#[test]
fn binary_pair_relaxation_refines_to_the_independent_exponential_solution() {
    let d = card(2e-5); let p = point(400.0, 100000.0);
    let pairs = [pair("A", "B", &d, &p)]; let mut errors = Vec::new();
    for count in [16, 32, 64] {
        let (mut a, mut b) = (gas(0.75), gas(0.25));
        for _ in 0..count {
            let proposal = a.diffuse_with(&b, geometry(), 10.0/f64::from(count), &pairs, controls(), || false).unwrap();
            assert!(proposal.entropy_change_j_k() >= -controls().entropy_tolerance_j_k);
            let (exchange, _, _) = proposal.into_parts(); let (next_a, next_b, _) = exchange.into_parts();
            a = next_a; b = next_b;
        }
        let gap = a.model().thermodynamics().mole_composition().fractions()[0]
            - b.model().thermodynamics().mole_composition().fractions()[0];
        let exact = 0.5*(-0.4_f64).exp();
        errors.push((gap-exact).abs());
    }
    assert!(errors.windows(2).all(|e| e[0]/e[1] > 1.9 && e[0]/e[1] < 2.1));
}

#[test]
fn source_support_pair_coverage_and_unresolved_intervals_refuse_atomically() {
    let (a, b) = (gas(0.75), gas(0.25)); let d = card(2e-5); let p = point(400.0, 100000.0);
    let saved = (a.clone(), b.clone()); let pairs = [pair("A", "B", &d, &p)];
    assert!(a.diffusion_flux_to(&b, geometry(), &[], controls(), || false).is_err());
    assert!(a.diffusion_flux_to(&b, geometry(), &[pair("A", "A", &d, &p)], controls(), || false).is_err());
    let mut hot = b.clone(); hot.add_heat(1.0, controls().exchange.gas[1], || false).unwrap();
    assert!(a.diffusion_flux_to(&hot, geometry(), &pairs, controls(), || false).is_err());
    assert!(a.diffuse_with(&b, geometry(), 1e8, &pairs, controls(), || false).is_err());
    let bad = card(-2e-5);
    assert!(a.diffusion_flux_to(&b, geometry(), &[pair("A", "B", &bad, &p)], controls(), || false).is_err());
    let mut calls = 0;
    a.diffuse_with(&b, geometry(), 0.1, &pairs, controls(), || { calls += 1; false }).unwrap();
    for stop in [1, calls/2, calls] {
        let mut count = 0;
        assert_eq!(a.diffuse_with(&b, geometry(), 0.1, &pairs, controls(), || { count += 1; count == stop }),
            Err(GasEnergyError::Cancelled));
        assert_eq!((&a, &b), (&saved.0, &saved.1));
    }
}

#[test]
fn closed_path_and_identical_compositions_are_bitwise_noops() {
    let (a, b) = (gas(0.75), gas(0.25));
    let trial = a.diffuse_with(&b, GasDiffusionGeometry { area_m2: 0.0, ..geometry() },
        0.1, &[], controls(), || false).unwrap();
    assert_eq!(trial.exchange().left(), &a); assert_eq!(trial.exchange().right(), &b);
    let d = card(2e-5); let p = point(400.0, 100000.0);
    let trial = a.diffuse_with(&a, geometry(), 0.1, &[pair("A", "B", &d, &p)], controls(), || false).unwrap();
    assert_eq!(trial.exchange().left(), &a); assert_eq!(trial.exchange().right(), &a);
    assert_eq!(trial.entropy_change_j_k(), 0.0);
}

mod network;
