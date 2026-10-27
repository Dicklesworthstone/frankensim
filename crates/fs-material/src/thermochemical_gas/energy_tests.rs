//! Numerical closure tests against independent synthetic caloric formulas.
use super::*;
use super::super::tests::{close, model, point};
use crate::fluid_transport::FluidTransport;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn controls() -> GasEnergySolve {
    GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1500.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 }
}
fn specific_u(t: f64, a_fraction: f64, slope: f64) -> f64 {
    let mass = a_fraction * 0.028 + (1.0 - a_fraction) * 0.004;
    let a = a_fraction * 3.5 + (1.0 - a_fraction) * 2.5;
    let b = a_fraction * slope;
    (R / mass) * ((a - 1.0) * t + 0.5 * b * t * t)
}

#[test]
fn source_isochoric_slope_bounds_actual_cv_without_replacing_nasa_energy() {
    for fractions in [vec![0.75, 0.25], vec![0.25, 0.75]] {
        let gas = model(fractions, 0.001);
        let bound = gas.maximum_isochoric_temperature_slope(controls()).unwrap();
        for temperature in [200.0, 500.0, 1000.0, 1500.0] {
            let state = gas.evaluate(&point(temperature, 100000.0)).unwrap();
            let cv = state.thermodynamics().properties().mass_specific_cv().value();
            assert!(bound >= 1.0 / cv);
        }
        let outside = GasEnergySolve { maximum_temperature_k: 2100.0, ..controls() };
        assert!(gas.maximum_isochoric_temperature_slope(outside).is_err());
    }
}

#[test]
fn fixed_volume_heat_changes_pressure_and_the_complete_shared_medium() {
    let mut parcel = IsochoricGasParcel::new(model(vec![0.75, 0.25], 0.0), point(300.0, 100000.0), 0.01).unwrap();
    let before = parcel.state().clone();
    let mass = parcel.mass_kg();
    let q = mass * (specific_u(600.0, 0.75, 0.0) - specific_u(300.0, 0.75, 0.0));
    let report = parcel.add_heat(q, controls(), || false).unwrap();
    close(parcel.state().state().temperature, 600.0);
    close(parcel.state().state().pressure, 200000.0);
    close(parcel.state().state().density, before.state().density);
    close(parcel.state().state().sound_speed / before.state().sound_speed, 2.0f64.sqrt());
    assert!(parcel.state().state().dynamic_viscosity > before.state().dynamic_viscosity);
    assert!(parcel.state().transport_state().unwrap().thermal_diffusivity_m2_s() >
        before.transport_state().unwrap().thermal_diffusivity_m2_s());
    assert_eq!(parcel.mass_kg(), mass);
    assert_eq!(parcel.volume_m3(), 0.01);
    assert!(report.residual_j.abs() <= mass * controls().energy_tolerance_j_kg);
    close(report.stored_energy_change_j, q);
    close(*parcel.point().axes().get("temperature").unwrap(), 600.0);
}

#[test]
fn nonlinear_caloric_source_changes_heat_response_without_a_new_solver() {
    let mut parcel = IsochoricGasParcel::new(model(vec![0.75, 0.25], 0.001), point(300.0, 100000.0), 0.01).unwrap();
    let initial = parcel.state().clone();
    let du = specific_u(700.0, 0.75, 0.001) - specific_u(300.0, 0.75, 0.001);
    let q = parcel.mass_kg() * du;
    let cv0 = initial.thermodynamics().properties().mass_specific_cv().value();
    assert!((300.0 + du / cv0 - 700.0).abs() > 10.0);
    let result = parcel.add_heat(q, controls(), || false).unwrap();
    close(parcel.state().state().temperature, 700.0);
    close(parcel.state().state().pressure, 100000.0 * 700.0 / 300.0);
    assert!(parcel.state().state().gamma < initial.state().gamma);
    assert!(result.iterations > 1);
    let inverse = parcel.add_heat(-q, controls(), || false).unwrap();
    close(parcel.state().state().temperature, 300.0);
    assert!(inverse.residual_j.abs() <= parcel.mass_kg() * controls().energy_tolerance_j_kg);
}

#[test]
fn conserved_density_energy_reconstruction_roundtrips_and_preserves_other_conditions() {
    let model = model(vec![0.25, 0.75], 0.001);
    let context = point(300.0, 100000.0).with("material-lot", 12.0).unwrap();
    let density = model.evaluate(&context).unwrap().state().density;
    let result = model.reconstruct_density_energy(&context, density,
        specific_u(725.0, 0.25, 0.001), controls(), || false).unwrap();
    close(result.state().state().temperature, 725.0);
    close(result.state().state().density, density);
    let mut parcel = IsochoricGasParcel::new(model, context, 0.01).unwrap();
    parcel.add_heat(100.0, controls(), || false).unwrap();
    assert_eq!(parcel.point().axes().get("material-lot"), Some(&12.0));
}

#[test]
fn failure_cancellation_and_work_exhaustion_are_atomic() {
    let mut parcel = IsochoricGasParcel::new(model(vec![0.75, 0.25], 0.001), point(300.0, 100000.0), 0.01).unwrap();
    let saved = parcel.clone();
    assert!(matches!(parcel.add_heat(10.0, controls(), || true), Err(GasEnergyError::Cancelled)));
    assert_eq!(parcel, saved);
    assert!(parcel.add_heat(f64::NAN, controls(), || false).is_err());
    assert_eq!(parcel, saved);
    assert!(matches!(parcel.add_heat(1e9, controls(), || false), Err(GasEnergyError::OutsideBracket { .. })));
    assert_eq!(parcel, saved);
    let limited = GasEnergySolve { maximum_iterations: 1, ..controls() };
    assert!(matches!(parcel.add_heat(100.0, limited, || false), Err(GasEnergyError::NotConverged { .. })));
    assert_eq!(parcel, saved);
    let mut polls = 0;
    assert!(matches!(parcel.add_heat(100.0, controls(), || { polls += 1; polls == 4 }), Err(GasEnergyError::Cancelled)));
    assert_eq!(parcel, saved);
    // Target is inside the thermodynamic bracket but pressure crosses the
    // pure-transport source pressure limit. No tentative state may become observable.
    let mut dense = IsochoricGasParcel::new(model(vec![0.75, 0.25], 0.001), point(300.0, 8e6), 0.01).unwrap();
    let saved_dense = dense.clone();
    let q = dense.mass_kg() * (specific_u(600.0, 0.75, 0.001) - specific_u(300.0, 0.75, 0.001));
    assert!(dense.add_heat(q, controls(), || false).is_err());
    assert_eq!(dense, saved_dense);
}

#[test]
fn zero_heat_is_bitwise_unchanged_and_bad_controls_are_not_ignored() {
    let mut parcel = IsochoricGasParcel::new(model(vec![0.75, 0.25], 0.001), point(300.0, 100000.0), 0.01).unwrap();
    let saved = parcel.clone();
    let step = parcel.add_heat(0.0, controls(), || false).unwrap();
    assert_eq!(parcel, saved);
    assert_eq!(step.stored_energy_change_j, 0.0);
    for invalid in [GasEnergySolve { maximum_iterations: 0, ..controls() },
        GasEnergySolve { energy_tolerance_j_kg: 0.0, ..controls() },
        GasEnergySolve { minimum_temperature_k: 1600.0, ..controls() }] {
        assert!(parcel.add_heat(0.0, invalid, || false).is_err());
        assert_eq!(parcel, saved);
    }
}

#[test]
fn source_supported_heating_and_compression_cross_legacy_air_temperature_limits() {
    use super::super::tests::model_in_domains;
    use super::entropy::GasEntropySolve;
    use super::parcel::ClosedGasParcel;

    let model = model_in_domains(
        vec![0.75, 0.25],
        0.001,
        [200.0, 6000.0],
        [200.0, 5000.0],
        [1000.0, 1.0e7],
    );
    let mut parcel = ClosedGasParcel::new(model, point(1800.0, 100000.0), 0.01).unwrap();
    let mass = parcel.mass_kg();
    let energy = GasEnergySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 4900.0,
        energy_tolerance_j_kg: 1e-6,
        maximum_iterations: 128,
    };
    let entropy = GasEntropySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 4900.0,
        entropy_tolerance_j_kg_k: 1e-8,
        maximum_iterations: 128,
    };
    let heat = mass * (specific_u(3200.0, 0.75, 0.001) - specific_u(1800.0, 0.75, 0.001));
    let report = parcel.add_heat(heat, energy, || false).unwrap();
    close(parcel.state().state().temperature, 3200.0);
    close(parcel.state().state().pressure, 100000.0 * 3200.0 / 1800.0);
    close(report.stored_energy_change_j, heat);
    parcel.state().state().validate_coherence().unwrap();
    assert_eq!(
        parcel.state().state().transport_state().unwrap(),
        parcel.state().transport_state().unwrap()
    );

    let heated = parcel.state().clone();
    let work = parcel
        .set_volume_isentropic(0.005, entropy, || false)
        .unwrap();
    let compressed = parcel.state().state();
    assert!(compressed.temperature > 3200.0);
    assert!(compressed.sound_speed > heated.state().sound_speed);
    close(compressed.density, 2.0 * heated.state().density);
    assert!(work.work_by_gas_j < 0.0);
    assert!(work.entropy_change_j_k.abs() <= mass * entropy.entropy_tolerance_j_kg_k);
    // Independent integral of cv(T)/T - R_specific d(ln rho) for cp/R = a+bT.
    let specific_r = R / (0.75 * 0.028 + 0.25 * 0.004);
    let a = 0.75 * 3.5 + 0.25 * 2.5;
    let b = 0.75 * 0.001;
    let ds = specific_r
        * ((a - 1.0) * (compressed.temperature / 3200.0).ln()
            + b * (compressed.temperature - 3200.0)
            - 2.0f64.ln());
    assert!(ds.abs() < 1e-8);
    compressed.validate_coherence().unwrap();

    // Returning to the original volume restores the entire thermal/acoustic state.
    parcel
        .set_volume_isentropic(0.01, entropy, || false)
        .unwrap();
    close(parcel.state().state().temperature, 3200.0);
    close(
        parcel.state().state().sound_speed,
        heated.state().sound_speed,
    );
    close(parcel.mass_kg(), mass);
    // Thermodynamic support alone cannot authorize unsupported transport.
    let saved = parcel.clone();
    let invalid_heat =
        mass * (specific_u(5500.0, 0.75, 0.001) - parcel.state().specific_internal_energy_j_kg());
    let wider = GasEnergySolve {
        maximum_temperature_k: 5900.0,
        ..energy
    };
    assert!(matches!(
        parcel.add_heat(invalid_heat, wider, || false),
        Err(GasEnergyError::Material(_))
    ));
    assert_eq!(parcel, saved);
}

#[test]
fn sensible_energy_reconstruction_is_invariant_to_a_species_reference_offset() {
    use super::super::tests::{species_with_energy_reference, transport};
    use fs_thermochem::{Composition, CompositionBasis};
    use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;

    let context = point(300.0, 100000.0).with("material-lot", 12.0).unwrap();
    let target_temperature = 725.0;
    let baseline = model(vec![0.75, 0.25], 0.001);
    let initial = baseline.evaluate(&context).unwrap();
    let density = initial.state().density;
    let cv_initial = initial.thermodynamics().properties().mass_specific_cv().value();
    let sensible_energy = specific_u(target_temperature, 0.75, 0.001)
        - specific_u(300.0, 0.75, 0.001)
        + cv_initial * 300.0;
    let expected = baseline
        .evaluate(&point(target_temperature, 100000.0 * target_temperature / 300.0))
        .unwrap();

    for offset in [-2.0e6, 0.0, 8.0e6] {
        let source = ThermochemicalGasModel::new(
            FrozenIdealGasMixtureModelV1::new(
                vec![
                    species_with_energy_reference("A", 0.028, 3.5, 0.001, offset),
                    species_with_energy_reference("B", 0.004, 2.5, 0.0, offset),
                ],
                Composition::new(CompositionBasis::MoleFraction, vec![0.75, 0.25]).unwrap(),
            ).unwrap(),
            vec![transport("B", 2.0e-5, 0.15), transport("A", 1.7e-5, 0.025)],
            None,
        ).unwrap();
        let reference = source.evaluate(&context).unwrap().specific_internal_energy_j_kg()
            - cv_initial * 300.0;
        let result = source.reconstruct_density_energy_change(
            &context, density, reference, sensible_energy, controls(), || false,
        ).unwrap();
        let accepted = result.state();
        close(accepted.state().temperature, target_temperature);
        close(accepted.state().pressure, expected.state().pressure);
        close(accepted.state().density, density);
        close(accepted.state().gamma, expected.state().gamma);
        close(accepted.state().sound_speed, expected.state().sound_speed);
        close(accepted.state().thermal_conductivity, expected.state().thermal_conductivity);
        assert!(result.residual_j_kg().abs() <= controls().energy_tolerance_j_kg);
        assert_eq!(accepted.components()[0].query().axes().get("material-lot"), Some(&12.0));
        assert!((accepted.specific_internal_energy_j_kg()
            - expected.specific_internal_energy_j_kg() - offset).abs() < 1e-6);
    }
}

#[test]
fn relative_energy_closure_keeps_signed_changes_and_refuses_invalid_inputs() {
    let source = model(vec![0.75, 0.25], 0.001);
    let context = point(600.0, 200000.0);
    let initial = source.evaluate(&context).unwrap();
    let density = initial.state().density;
    let reference = initial.specific_internal_energy_j_kg();
    let cooling = specific_u(300.0, 0.75, 0.001) - reference;
    let result = source.reconstruct_density_energy_change(
        &context, density, reference, cooling, controls(), || false,
    ).unwrap();
    close(result.state().state().temperature, 300.0);
    close(result.state().state().pressure, 100000.0);
    assert!(matches!(
        source.reconstruct_density_energy_change(
            &context, density, reference, cooling, controls(), || true,
        ),
        Err(GasEnergyError::Cancelled)
    ));
    for (invalid_reference, invalid_change) in [
        (f64::NAN, cooling),
        (reference, f64::INFINITY),
    ] {
        assert!(matches!(
            source.reconstruct_density_energy_change(
                &context, density, invalid_reference, invalid_change, controls(), || false,
            ),
            Err(GasEnergyError::Invalid(_))
        ));
    }
    assert_eq!(source.evaluate(&context).unwrap(), initial);
}
