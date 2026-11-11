//! G0/G1 original-gas conservation, lifecycle boundaries and analytical limits.
use super::*;
use crate::sublimation::gas::{GasSurfaceBinding, tests as fixture};
use fs_thermal::sublimation::{SublimationError, SublimationLaw, SublimationProperties};
use fs_thermochem::SpeciesId;

const SUPPORT_MASS: f64 = 0.02;
const SUPPORT_CP: f64 = 500.0;
const AREA: f64 = 1e-4;

fn controls() -> GasSupportedLifecycleSolve {
    let original = fixture::controls();
    GasSupportedLifecycleSolve {
        gas: original.gas,
        storage: SolverConfig {
            relative_tolerance: 1e-12,
            absolute_energy_tolerance_j: 1e-10,
            max_iterations: 128,
        },
        latent_compatibility_tolerance_j_kg: original.latent_compatibility_tolerance_j_kg,
        maximum_deposit_volume_fraction: 0.01,
    }
}

fn support(h: f64) -> Result<f64, SublimationError> {
    Ok(400.0 + h / SUPPORT_CP)
}

fn surface<'a>(
    temperature_k: f64,
    mass_kg: f64,
    law: &'a dyn SublimationLaw,
    binding: &'a GasSurfaceBinding,
    deposit_enthalpy: &'a dyn Fn(f64) -> Result<f64, SublimationError>,
) -> GasSupportedSurface<'a> {
    GasSupportedSurface {
        id: 17,
        law,
        binding,
        state: SupportedSurfaceState {
            substrate_mass_kg: SUPPORT_MASS,
            substrate_specific_enthalpy_j_kg: SUPPORT_CP * (temperature_k - 400.0),
            deposit_mass_kg: mass_kg,
            area_m2: AREA,
        },
        substrate_temperature: &support,
        deposit_enthalpy,
    }
}

fn balances(
    initial: &GasSupportedSurface<'_>,
    gas: &ClosedGasParcel,
    result: &GasSupportedLifecycleCandidate<'_>,
) {
    let row = result.endpoint().surfaces()[0];
    let t0 =
        (initial.substrate_temperature)(initial.state.substrate_specific_enthalpy_j_kg).unwrap();
    let h0 = (initial.deposit_enthalpy)(t0).unwrap();
    let h1 = (initial.deposit_enthalpy)(row.exchange.surface_temperature_k[1]).unwrap();
    let solid0 = initial.state.substrate_mass_kg * initial.state.substrate_specific_enthalpy_j_kg
        + initial.state.deposit_mass_kg * (h0 + initial.binding.solid_enthalpy_offset_j_kg);
    let solid1 = row.state.substrate_mass_kg * row.state.substrate_specific_enthalpy_j_kg
        + row.state.deposit_mass_kg * (h1 + initial.binding.solid_enthalpy_offset_j_kg);
    let gas_change = result.gas().mass_kg() * result.gas().state().specific_internal_energy_j_kg()
        - gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
    fixture::near(solid1 - solid0 + gas_change, 0.0, 2e-7);
    fixture::near(
        row.substrate_energy_change_j + row.deposit_energy_change_j,
        solid1 - solid0,
        2e-7,
    );
    fixture::near(result.endpoint().step().energy_residual_j, 0.0, 1e-7);
    fixture::near(
        row.state.deposit_mass_kg - initial.state.deposit_mass_kg + result.gas().mass_kg()
            - gas.mass_kg(),
        0.0,
        1e-14,
    );
    let carrier = SpeciesId::new("A").unwrap();
    fixture::near(
        result.gas().species_mass_kg(&carrier).unwrap(),
        gas.species_mass_kg(&carrier).unwrap(),
        1e-14,
    );
    assert_eq!(row.state.substrate_mass_kg, initial.state.substrate_mass_kg);
    assert_eq!(row.state.area_m2, initial.state.area_m2);
    assert_eq!(result.gas().volume_m3(), gas.volume_m3());
    assert_eq!(result.gas().point().axes().get("local-lot"), Some(&9.0));
    assert_eq!(result.catalog(), &fixture::catalog(0.0));
}

#[test]
fn bare_supersaturation_grows_with_concurrent_emission_and_original_storage() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let law = fixture::phase();
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let initial = surface(300.0, 0.0, &law, &binding, &deposit);
    let dt = 0.1;
    let result =
        advance_supported_lifecycle(&initial, None, &gas, &catalog, dt, controls(), || false)
            .unwrap();
    let streams = result.turnover();
    assert!(streams.outgoing_mass_kg > 0.0);
    assert!(streams.incoming_mass_kg > streams.outgoing_mass_kg);
    assert!(result.state().deposit_mass_kg > 0.0);
    assert_eq!(streams.extinction_time_s, None);

    // Independent closed-form finite-reservoir integral, at the original T.
    let r = fs_thermal::sublimation::MOLAR_GAS_CONSTANT / 0.018015;
    let outgoing = binding.accommodation_coefficient
        * law.properties(300.0).unwrap().saturation_pressure_pa
        / (2.0 * core::f64::consts::PI * r * 300.0).sqrt()
        * AREA;
    let incoming_coefficient =
        binding.accommodation_coefficient / (2.0 * core::f64::consts::PI * r * 400.0).sqrt() * AREA;
    let decay = incoming_coefficient * r * 400.0 / gas.volume_m3();
    let vapor0 = gas.species_mass_kg(&fixture::species()).unwrap();
    let dm = (outgoing / decay - vapor0) * -(-decay * dt).exp_m1();
    fixture::near(streams.transferred_mass_kg, dm, 1e-17);
    fixture::near(streams.outgoing_mass_kg, outgoing * dt, 1e-17);
    let heat =
        -dm * (200000.0 - 0.5 * r * 300.0) - streams.incoming_mass_kg * 3.5 * r * (300.0 - 400.0);
    let expected_t = 300.0 + heat / (SUPPORT_MASS * SUPPORT_CP - dm * chart.1);
    fixture::near(
        result.endpoint().surfaces()[0]
            .exchange
            .surface_temperature_k[1],
        expected_t,
        1e-8,
    );
    balances(&initial, &gas, &result);
}

#[test]
fn full_extinction_keeps_the_support_and_continues_bare_sensible_exchange() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let law = fixture::phase();
    let gas = fixture::gas(&catalog, 400.0, 100.0);
    let initial = surface(450.0, 1e-8, &law, &binding, &deposit);
    let dt = 0.02;
    let result =
        advance_supported_lifecycle(&initial, None, &gas, &catalog, dt, controls(), || false)
            .unwrap();
    assert_eq!(result.state().deposit_mass_kg, 0.0);
    assert_eq!(
        result.turnover().transferred_mass_kg,
        initial.state.deposit_mass_kg
    );
    let event = result.turnover().extinction_time_s.unwrap();
    assert!(event > 0.0 && event < dt);
    assert!(result.turnover().outgoing_mass_kg > initial.state.deposit_mass_kg);
    balances(&initial, &gas, &result);

    let bare = GasSupportedSurface {
        state: result.state(),
        ..initial
    };
    let continued =
        advance_supported_lifecycle(&bare, None, result.gas(), &catalog, dt, controls(), || {
            false
        })
        .unwrap();
    assert_eq!(continued.state().deposit_mass_kg, 0.0);
    assert_eq!(continued.turnover().transferred_mass_kg, 0.0);
    assert!(continued.turnover().incoming_mass_kg > 0.0);
    assert_eq!(
        continued.turnover().outgoing_mass_kg,
        continued.turnover().incoming_mass_kg
    );
    assert!(continued.endpoint().step().energy_into_gas_j > 0.0);
    assert!(
        continued.state().substrate_specific_enthalpy_j_kg
            < bare.state.substrate_specific_enthalpy_j_kg
    );
    balances(&bare, result.gas(), &continued);
}

#[test]
fn absent_vapor_is_emitted_mixed_and_reabsorbed_in_the_same_interval() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let law = fixture::phase();
    let gas = fixture::gas(&catalog, 400.0, 0.0);
    let initial = surface(450.0, 1e-6, &law, &binding, &deposit);
    assert_eq!(gas.species_mass_kg(&fixture::species()).unwrap(), 0.0);
    let result =
        advance_supported_lifecycle(&initial, None, &gas, &catalog, 0.02, controls(), || false)
            .unwrap();
    assert!(result.turnover().outgoing_mass_kg > result.turnover().incoming_mass_kg);
    assert!(result.turnover().incoming_mass_kg > 0.0);
    assert!(
        result.endpoint().surfaces()[0]
            .exchange
            .incoming_specific_energy_j_kg
            .is_some()
    );
    assert!(result.gas().species_mass_kg(&fixture::species()).unwrap() > 0.0);
    assert!(result.state().deposit_mass_kg > 0.0);
    balances(&initial, &gas, &result);
}

#[test]
fn positive_film_small_steps_match_the_original_explicit_chamber_limit() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let law = fixture::phase();
    let gas = fixture::gas(&catalog, 420.0, 2000.0);
    let initial = surface(380.0, 1e-5, &law, &binding, &deposit);
    let dt = 1e-7;
    let result =
        advance_supported_lifecycle(&initial, None, &gas, &catalog, dt, controls(), || false)
            .unwrap();
    let explicit = super::super::advance_supported_chamber(
        &[initial],
        None,
        &gas,
        &catalog,
        dt,
        controls().chamber_controls(),
        || false,
    )
    .unwrap();
    fixture::near(
        result.state().deposit_mass_kg,
        explicit.surfaces()[0].state.deposit_mass_kg,
        1e-18,
    );
    fixture::near(
        result.gas().state().state().temperature,
        explicit.gas().state().state().temperature,
        1e-8,
    );
    fixture::near(
        result.endpoint().surfaces()[0]
            .exchange
            .surface_temperature_k[1],
        explicit.surfaces()[0].exchange.surface_temperature_k[1],
        1e-8,
    );
    balances(&initial, &gas, &result);
}

#[test]
fn kinetic_and_caloric_results_do_not_depend_on_energy_datums() {
    let run = |species_datum: f64, deposit_offset: f64, support_datum: f64| {
        let catalog = fixture::catalog(species_datum);
        let binding = fixture::binding(deposit_offset);
        let chart = fixture::chart(&catalog, &binding);
        let deposit = |t| Ok(chart.0 + chart.1 * (t - 400.0));
        let substrate = |h| Ok(400.0 + (h - support_datum) / SUPPORT_CP);
        let law = fixture::phase();
        let gas = fixture::gas(&catalog, 400.0, 3000.0);
        let mut initial = surface(350.0, 0.0, &law, &binding, &deposit);
        initial.substrate_temperature = &substrate;
        initial.state.substrate_specific_enthalpy_j_kg += support_datum;
        // The original enthalpy encoding at a 1e9 J/kg datum has a finite
        // extensive-energy quantum; admit that quantum without changing laws.
        let controls = GasSupportedLifecycleSolve {
            storage: SolverConfig {
                absolute_energy_tolerance_j: 1e-8,
                ..controls().storage
            },
            ..controls()
        };
        let result =
            advance_supported_lifecycle(&initial, None, &gas, &catalog, 0.1, controls, || false)
                .unwrap();
        [
            result.state().deposit_mass_kg,
            result.endpoint().surfaces()[0]
                .exchange
                .surface_temperature_k[1],
            result.gas().state().state().temperature,
        ]
    };
    let original = run(0.0, 0.0, 0.0);
    let shifted = run(1e6, 1e5, 1e9);
    fixture::near(original[0], shifted[0], 1e-15);
    fixture::near(original[1], shifted[1], 1e-7);
    fixture::near(original[2], shifted[2], 1e-7);
}

#[test]
fn exact_stationary_sealed_and_zero_time_states_retain_original_owners() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let law = fixture::phase();
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    for mass in [0.0, 1e-5] {
        let initial = surface(400.0, mass, &law, &binding, &deposit);
        let result =
            advance_supported_lifecycle(&initial, None, &gas, &catalog, 20.0, controls(), || false)
                .unwrap();
        assert_eq!(result.state(), initial.state);
        assert_eq!(result.gas(), &gas);
        assert_eq!(
            result.turnover().outgoing_mass_kg,
            result.turnover().incoming_mass_kg
        );
        assert_eq!(result.turnover().transferred_mass_kg, 0.0);
    }
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let initial = surface(350.0, 0.0, &law, &binding, &deposit);
    for (film, dt) in [(None, 0.0), (Some(0.0), 0.1)] {
        let result =
            advance_supported_lifecycle(&initial, film, &gas, &catalog, dt, controls(), || false)
                .unwrap();
        assert_eq!(result.state(), initial.state);
        assert_eq!(result.gas(), &gas);
        assert_eq!(result.turnover().outgoing_mass_kg, 0.0);
        assert_eq!(result.turnover().incoming_mass_kg, 0.0);
        assert_eq!(result.turnover().extinction_time_s, None);
    }
}

#[test]
fn sources_volume_controls_and_late_cancellation_refuse_transactionally() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let law = fixture::phase();
    let gas = fixture::gas(&catalog, 400.0, 3000.0);
    let original_gas = gas.clone();
    let initial = surface(350.0, 1e-5, &law, &binding, &deposit);
    let foreign = fixture::catalog(1.0);
    assert!(
        advance_supported_lifecycle(&initial, None, &gas, &foreign, 0.0, controls(), || false,)
            .is_err()
    );
    let tiny_volume = GasSupportedLifecycleSolve {
        maximum_deposit_volume_fraction: 1e-9,
        ..controls()
    };
    assert!(
        advance_supported_lifecycle(&initial, None, &gas, &catalog, 0.0, tiny_volume, || false,)
            .is_err()
    );
    for film in [Some(-1.0), Some(f64::NAN)] {
        assert!(
            advance_supported_lifecycle(&initial, film, &gas, &catalog, 0.0, controls(), || false,)
                .is_err()
        );
    }
    let mut polls = 0;
    let cancelled =
        advance_supported_lifecycle(&initial, None, &gas, &catalog, 0.1, controls(), || {
            polls += 1;
            polls >= 6
        });
    assert!(matches!(
        cancelled,
        Err(GasSurfaceError::Surface(SublimationError::Cancelled))
    ));
    assert!(polls >= 6);
    assert_eq!(gas, original_gas);
    assert_eq!(initial.state.deposit_mass_kg, 1e-5);
}

#[test]
fn positive_asymptotic_vapor_is_not_deleted_by_net_transfer_rounding() {
    struct ZeroSaturation;
    impl SublimationLaw for ZeroSaturation {
        fn temperature_range_k(&self) -> [f64; 2] {
            [200.0, 900.0]
        }
        fn properties(&self, t: f64) -> Result<SublimationProperties, SublimationError> {
            Ok(SublimationProperties {
                saturation_pressure_pa: 0.0,
                saturation_pressure_derivative_pa_k: 0.0,
                ..fixture::phase().properties(t)?
            })
        }
    }
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = fixture::chart(&catalog, &binding);
    let deposit = |t| Ok(chart.0 + chart.1 * (t - 400.0));
    let gas = fixture::gas(&catalog, 400.0, 1000.0);
    let initial = surface(400.0, 0.0, &ZeroSaturation, &binding, &deposit);
    let result =
        advance_supported_lifecycle(&initial, None, &gas, &catalog, 2000.0, controls(), || false);
    assert!(matches!(
        result,
        Err(GasSurfaceError::Invalid(
            "integrated molecular net transfer cannot represent a positive inventory tail"
        ))
    ));
}

#[test]
fn a_net_transfer_encodes_a_tail_only_in_the_tails_own_inventory_class() {
    // (old inventory plus net transfer, integrated tail, encoded). The first
    // row is the printed case of the test above: -6.776263578034403e-21 kg
    // against a 7.04e-154 kg tail. A zero tail keeps its old rule: only an
    // exactly zero sum encodes it, so a negative sum still refuses there.
    let tail = 7.040127103708638e-154;
    for (sum, tail, encoded) in [
        (-6.776263578034403e-21, tail, false),
        (0.0, tail, false),
        (-0.0, tail, false),
        (1e-30, tail, true),
        (-1e-30, 0.0, false),
        (1e-30, 0.0, false),
        (0.0, 0.0, true),
        (-0.0, 0.0, true),
        (f64::NAN, 0.0, false),
        (f64::INFINITY, 1.0, false),
        (1.0, f64::NAN, false),
    ] {
        assert_eq!(
            net_transfer_encodes_tail(sum, tail),
            encoded,
            "sum {sum:e}, tail {tail:e}"
        );
    }
}
