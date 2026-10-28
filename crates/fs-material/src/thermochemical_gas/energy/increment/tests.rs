use super::*;
use crate::thermochemical_gas::energy::IsochoricGasParcel;
use crate::thermochemical_gas::tests::{close, model, point, species_in_domain,
    species_with_energy_reference, transport};
use fs_thermochem::{Composition, CompositionBasis, MolarEnergyQuantityV1,
    MolarThermalQuantityV1, Nasa9StandardStateModelV1,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;

fn controls() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1500.0,
        energy_tolerance_j_kg: 1e-8, maximum_iterations: 128 }
}

fn shifted(reference: f64) -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(
        vec![species_with_energy_reference("A", 0.028, 3.5, 0.001, reference),
            species_with_energy_reference("B", 0.004, 2.5, 0.0, -reference)],
        Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
    ).unwrap(), vec![transport("A", 1.7e-5, 0.025), transport("B", 2e-5, 0.15)], None).unwrap()
}

fn delta(t0: f64, t1: f64, internal: bool) -> f64 {
    let constant = if internal { 2.25 } else { 3.25 };
    (R / (0.75 * 0.028 + 0.25 * 0.004)) * (t1 - t0)
        * (constant + 0.000375 * (t1 + t0))
}

#[test]
fn g1_isochoric_heat_is_reference_invariant_in_the_actual_parcel_owner() {
    for target in [300.0 + 1e-6, 700.0] {
        let limits = GasEnergySolve { energy_tolerance_j_kg: 1e-12, ..controls() };
        let mut baseline = IsochoricGasParcel::new(shifted(0.0),
            point(300.0, 100000.0).with("material-lot", 42.0).unwrap(), 0.01).unwrap();
        let heat = baseline.mass_kg() * delta(300.0, target, true);
        // Large steps need an ordinary absolute residual target; small ones
        // deliberately use a tolerance far below the absolute-reference ULP.
        let limits = if target > 400.0 { controls() } else { limits };
        let expected = baseline.add_heat(heat, limits, || false).unwrap();
        close(baseline.state().state().temperature, target);
        for reference in [-1e24, 1e24] {
            let mut parcel = IsochoricGasParcel::new(shifted(reference),
                point(300.0, 100000.0).with("material-lot", 42.0).unwrap(), 0.01).unwrap();
            let initial = parcel.clone();
            let step = parcel.add_heat(heat, limits, || false).unwrap();
            assert_eq!(step, expected);
            assert_eq!(parcel.state().state(), baseline.state().state());
            assert_eq!(parcel.point().axes().get("material-lot"), Some(&42.0));
            assert_eq!(parcel.mass_kg(), initial.mass_kg());
            assert_eq!(parcel.volume_m3(), initial.volume_m3());
            assert!(parcel.state().specific_internal_energy_j_kg().abs() > 1e20);
            assert!(step.residual_j.abs() <= parcel.mass_kg() * limits.energy_tolerance_j_kg);
            close(step.stored_energy_change_j, heat);
        }
    }
}

#[test]
fn g1_isobaric_increment_uses_cp_and_preserves_pressure_and_sources() {
    let context = point(300.0, 100000.0);
    let expected = shifted(0.0).reconstruct_pressure_enthalpy_increment(
        &context, delta(300.0, 800.0, false), controls(), || false,
    ).unwrap();
    for reference in [-1e24, 0.0, 1e24] {
        let gas = shifted(reference);
        let solved = gas.reconstruct_pressure_enthalpy_increment(
            &context, delta(300.0, 800.0, false), controls(), || false,
        ).unwrap();
        close(solved.state().state().temperature, 800.0);
        assert_eq!(solved.state().state().pressure, 100000.0);
        assert_eq!(solved.state().state(), expected.state().state());
        close(solved.enthalpy_change_j_kg(), delta(300.0, 800.0, false));
        assert_eq!(solved.state().components().len(), 2);
        let cooled = gas.reconstruct_pressure_enthalpy_increment(
            &point(solved.state().state().temperature, 100000.0),
            -solved.enthalpy_change_j_kg(), controls(), || false,
        ).unwrap();
        close(cooled.state().state().temperature, 300.0);
    }
}

fn two_regions() -> ThermochemicalGasModel {
    let low = species_in_domain("A", 0.028, 3.5, 0.0, 0.0, [200.0, 600.0]);
    let high = species_in_domain("A", 0.028, 4.0, 0.0, 0.0, [600.0, 2000.0]);
    let species = Nasa9StandardStateModelV1::new(low.species().clone(), low.molar_mass(),
        low.convention().clone(), vec![low.regions()[0].clone(), high.regions()[0].clone()]).unwrap();
    ThermochemicalGasModel::new(FrozenIdealGasMixtureModelV1::new(vec![species],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap()).unwrap(),
        vec![transport("A", 1.7e-5, 0.025)], None).unwrap()
}

#[test]
fn g3_increment_requires_continuous_sources_across_a_caloric_seam() {
    let raw = two_regions();
    let context = point(300.0, 100000.0);
    let heat = (R / 0.028) * (3.5 * 300.0 + 4.0 * 200.0);
    assert!(matches!(raw.reconstruct_pressure_enthalpy_increment(
        &context, heat, controls(), || false,
    ), Err(GasEnergyError::Invalid(_))));
    let continuous = raw.with_caloric_continuity(Temperature::new(300.0),
        MolarEnergyQuantityV1::new(1e5), MolarThermalQuantityV1::new(1e3)).unwrap();
    let result = continuous.reconstruct_pressure_enthalpy_increment(
        &context, heat, controls(), || false,
    ).unwrap();
    close(result.state().state().temperature, 800.0);
    let density = continuous.evaluate(&context).unwrap().state().density;
    let result = continuous.reconstruct_density_internal_energy_increment(
        &context, density, (R / 0.028) * (2.5 * 300.0 + 3.0 * 200.0), controls(), || false,
    ).unwrap();
    close(result.state().state().temperature, 800.0);
}

#[test]
fn g4_increment_refusals_do_not_publish_partial_parcels() {
    let mut parcel = IsochoricGasParcel::new(shifted(1e24), point(300.0, 100000.0), 0.01).unwrap();
    let saved = parcel.clone();
    for (heat, limits) in [
        (f64::NAN, controls()),
        (1e9, controls()),
        (100.0, GasEnergySolve { maximum_iterations: 1, ..controls() }),
        (100.0, GasEnergySolve { minimum_temperature_k: 400.0, ..controls() }),
        (100.0, GasEnergySolve { maximum_temperature_k: 2100.0, ..controls() }),
    ] {
        assert!(parcel.add_heat(heat, limits, || false).is_err());
        assert_eq!(parcel, saved);
    }
    let mut completed_polls = 0;
    let mut probe = parcel.clone();
    probe.add_heat(100.0, controls(), || { completed_polls += 1; false }).unwrap();
    for stop in [1, completed_polls / 2, completed_polls] {
        let mut polls = 0;
        assert!(matches!(parcel.add_heat(100.0, controls(), || {
            polls += 1; polls == stop
        }), Err(GasEnergyError::Cancelled)));
        assert_eq!(parcel, saved);
    }
}

#[test]
fn g0_zero_increments_preserve_reference_temperature_and_validate_inputs() {
    let gas = model(vec![0.75, 0.25], 0.001);
    let context = point(525.0, 100000.0);
    let solved = gas.reconstruct_pressure_enthalpy_increment(&context, 0.0, controls(), || false).unwrap();
    assert_eq!(solved.state().state().temperature, 525.0);
    assert_eq!(solved.enthalpy_change_j_kg(), 0.0);
    assert_eq!(solved.iterations(), 0);
    assert!(gas.reconstruct_density_internal_energy_increment(&context, 0.0, 0.0, controls(), || false).is_err());
    assert!(gas.reconstruct_pressure_enthalpy_increment(&context, 0.0,
        GasEnergySolve { maximum_iterations: 0, ..controls() }, || false).is_err());
}
