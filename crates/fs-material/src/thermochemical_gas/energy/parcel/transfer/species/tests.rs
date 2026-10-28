//! Real source-backed endpoints; synthetic coefficients are not measured gases.
use super::*;
use crate::thermochemical_gas::tests::{model, point, species, transport};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;

fn energy() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 }
}
fn controls() -> GasSpeciesExchangeSolve {
    GasSpeciesExchangeSolve { gas: [energy(); 2], maximum_species_fraction: 0.2,
        mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8 }
}
fn parcel(x: f64, t: f64, slope: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(model(vec![x, 1.0-x], slope),
        point(t, 100000.0).with("local-lot", x).unwrap(), 0.001).unwrap()
}
fn pure(name: &str, mass: f64, cp: f64) -> ClosedGasParcel {
    let model = ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(
        vec![species(name, mass, cp, 0.0)],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap()).unwrap(),
        vec![transport(name, 2e-5, 0.03)], None).unwrap();
    ClosedGasParcel::new(model, point(400.0, 100000.0), 0.001).unwrap()
}
fn row(name: &str, mass: f64) -> SpeciesMassTransfer {
    SpeciesMassTransfer { species: SpeciesId::new(name).unwrap(), mass_left_to_right_kg: mass }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance, "{a} != {b}");
}

#[test]
fn g1_selective_counterflow_uses_partial_enthalpy_with_reference_invariant_temperatures() {
    use super::super::tests::reference_parcel;
    use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
    let pair = |reference| [
        reference_parcel(0.75, 0.0, 500.0, 100000.0, reference),
        reference_parcel(0.25, 0.0, 350.0, 100000.0, reference),
    ];
    let original = pair([0.0; 2]);
    let shifted = pair([1e24, -2e24]);
    let rows = [row("A", 1e-5), row("B", -2e-6)];
    let baseline = original[0].exchange_species_with(&original[1], &rows, controls(), || false).unwrap();
    let trial = shifted[0].exchange_species_with(&shifted[1], &rows, controls(), || false).unwrap();
    let cv = [2.5 * R / 0.028, 1.5 * R / 0.004];
    let cp = [3.5 * R / 0.028, 2.5 * R / 0.004];
    let flux = rows[0].mass_left_to_right_kg * cp[0] * 500.0
        + rows[1].mass_left_to_right_kg * cp[1] * 350.0;
    for (side, sign, actual, expected) in [
        (0, -1.0, trial.left(), baseline.left()), (1, 1.0, trial.right(), baseline.right()),
    ] {
        let original_inventory = inventory(original[side].model(), original[side].mass_kg()).unwrap();
        let a = original_inventory[&SpeciesId::new("A").unwrap()];
        let b = original_inventory[&SpeciesId::new("B").unwrap()];
        let initial_u = (a * cv[0] + b * cv[1]) * original[side].state().state().temperature;
        let final_capacity = (a + sign * rows[0].mass_left_to_right_kg) * cv[0]
            + (b + sign * rows[1].mass_left_to_right_kg) * cv[1];
        let analytic_t = (initial_u + sign * flux) / final_capacity;
        near(actual.state().state().temperature, analytic_t, 2e-8);
        assert_eq!(actual.state().state(), expected.state().state());
        assert_eq!(actual.mass_kg(), expected.mass_kg());
    }
    assert_eq!(trial.step().energy_residual_j, baseline.step().energy_residual_j);
    assert!(trial.step().energy_residual_j.abs() <= controls().energy_tolerance_j);
    assert!(trial.step().maximum_species_mass_residual_kg <= controls().mass_tolerance_kg);
}

#[test]
fn counterflow_carries_each_original_partial_enthalpy_not_mixture_enthalpy() {
    let (a, b) = (parcel(0.75, 500.0, 0.001), parcel(0.25, 350.0, 0.001));
    let rows = [row("A", 1e-5), row("B", -2e-6)];
    let trial = a.exchange_species_with(&b, &rows, controls(), || false).unwrap();
    let mut expected = 0.0;
    for (i, upstream) in [(0, &a), (1, &b)] {
        let pure = upstream.model().thermodynamics().components()[i].evaluate(
            Temperature::new(upstream.state().state().temperature)).unwrap();
        let h = pure.properties().h().value()
            / upstream.model().thermodynamics().components()[i].molar_mass().value();
        near(trial.step().species[i].source_enthalpy_j_kg, h, 1e-7);
        assert!((h - upstream.state().specific_enthalpy_j_kg()).abs() > 1.0);
        expected += rows[i].mass_left_to_right_kg * h;
    }
    near(trial.step().enthalpy_left_to_right_j, expected, 1e-12);
    near(trial.step().energy_change_j[0], -expected, 1e-8);
    near(trial.step().energy_change_j[1], expected, 1e-8);
    assert!(trial.step().maximum_species_mass_residual_kg < 1e-14);
    assert_eq!(trial.left().point().axes().get("local-lot"), Some(&0.75));
    assert_eq!(trial.right().point().axes().get("local-lot"), Some(&0.25));
}

#[test]
fn equimolar_exchange_introduces_species_without_changing_equal_temperature_pressure() {
    let (a, b) = (pure("A", 0.028, 3.5), pure("B", 0.004, 2.5));
    let n = 0.002;
    let trial = a.exchange_species_with(&b, &[row("A", n*0.028), row("B", -n*0.004)],
        controls(), || false).unwrap();
    for (old, next) in [(&a, trial.left()), (&b, trial.right())] {
        near(next.state().state().temperature, 400.0, 1e-7);
        near(next.state().state().pressure, 100000.0, 1e-5);
        assert_eq!(next.model().thermodynamics().components().len(), 2);
        assert_eq!(next.volume_m3(), old.volume_m3());
        assert_ne!(next.state().state().density, old.state().state().density);
        assert_ne!(next.state().state().sound_speed, old.state().state().sound_speed);
        assert_ne!(next.state().state().specific_heat_cp, old.state().state().specific_heat_cp);
    }
    near(trial.left().mass_kg()-a.mass_kg(), -n*(0.028-0.004), 1e-14);
    near(trial.step().energy_residual_j, 0.0, 1e-8);
}

#[test]
fn whole_mixture_rows_reproduce_the_existing_upwind_pair_transfer() {
    let (a, b) = (parcel(0.75, 500.0, 0.001), parcel(0.25, 400.0, 0.001));
    let amount = 0.01*a.mass_kg();
    let reference = a.transfer_mass_to(&b, amount, GasMassTransferSolve {
        donor: energy(), receiver: energy(), mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8 }, || false).unwrap();
    let rows = inventory(a.model(), amount).unwrap().into_iter().map(|(species, mass)|
        SpeciesMassTransfer { species, mass_left_to_right_kg: mass }).collect::<Vec<_>>();
    let trial = a.exchange_species_with(&b, &rows, controls(), || false).unwrap();
    for (actual, expected) in [(trial.left(), reference.donor()), (trial.right(), reference.receiver())] {
        near(actual.mass_kg(), expected.mass_kg(), 1e-14);
        near(actual.state().state().temperature, expected.state().state().temperature, 1e-7);
        near(actual.state().state().pressure, expected.state().state().pressure, 1e-5);
    }
}

#[test]
fn canonical_rows_and_zero_exchange_preserve_physical_owners() {
    let (a, b) = (parcel(0.75, 500.0, 0.001), parcel(0.25, 400.0, 0.001));
    let rows = [row("A", 1e-5), row("B", -1e-6)];
    let first = a.exchange_species_with(&b, &rows, controls(), || false).unwrap();
    let reverse_rows = [rows[1].clone(), rows[0].clone()];
    assert_eq!(first, a.exchange_species_with(&b, &reverse_rows, controls(), || false).unwrap());
    let zero = a.exchange_species_with(&b, &[], controls(), || false).unwrap();
    assert_eq!(zero.left(), &a); assert_eq!(zero.right(), &b);
    assert_eq!(zero.step().iterations, 0);
}

#[test]
fn species_limits_conflicting_sources_and_late_failure_publish_neither_endpoint() {
    let (a, b) = (parcel(0.75, 500.0, 0.001), parcel(0.25, 400.0, 0.001));
    let saved = (a.clone(), b.clone());
    for rows in [vec![row("A", f64::NAN)], vec![row("foreign", 0.0)],
        vec![row("A", 1.0)], vec![row("A", 1e-5), row("A", 1e-5)]] {
        assert!(a.exchange_species_with(&b, &rows, controls(), || false).is_err());
    }
    assert!(a.exchange_species_with(&parcel(0.25, 400.0, 0.002), &[row("A", 1e-5)],
        controls(), || false).is_err());
    let narrow = GasSpeciesExchangeSolve { gas: [energy(), GasEnergySolve {
        maximum_temperature_k: 250.0, ..energy() }], ..controls() };
    assert!(a.exchange_species_with(&b, &[row("A", 1e-5)], narrow, || false).is_err());
    let mut calls = 0;
    a.exchange_species_with(&b, &[row("A", 1e-5)], controls(), || { calls += 1; false }).unwrap();
    for stop in [1, calls/2, calls] {
        let mut count = 0;
        assert_eq!(a.exchange_species_with(&b, &[row("A", 1e-5)], controls(), || {
            count += 1; count == stop }), Err(GasEnergyError::Cancelled));
        assert_eq!((&a, &b), (&saved.0, &saved.1));
    }
}

#[test]
fn selective_endpoints_continue_through_the_original_heat_work_and_load_paths() {
    let (a, b) = (parcel(0.75, 500.0, 0.001), parcel(0.25, 400.0, 0.001));
    let trial = a.exchange_species_with(&b, &[row("A", 1e-5), row("B", -1e-6)], controls(), || false).unwrap();
    let motion = trial.right().midpoint_heat_volume(0.0011, 3.0, energy(), || false).unwrap();
    let gradient = [0.01];
    let load = motion.discrete_mechanical_load(0.0011, &gradient, 100000.0).unwrap();
    assert!(load.generalized_forces().unwrap()[0].is_finite());
    near(stored_change(&a, trial.left()) + stored_change(&b, motion.endpoint())
        + motion.step().work_by_gas_j, 3.0, 1e-8);
}

mod network;
