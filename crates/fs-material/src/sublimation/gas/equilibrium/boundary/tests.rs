//! Independent sharp-fusion, molecular-energy and explicit-handoff controls.
use super::*;
use crate::sublimation::gas::{
    equilibrium::{advance_equilibrium_surface, tests as phase_fixture},
    tests as fixture,
};
use crate::thermochemical_gas::energy::GasEnergyError;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;

const FUSION: f64 = 100000.0;
const MOLAR_MASS: f64 = 0.018015;

#[test]
fn reflected_heat_completes_melting_at_its_analytic_time_including_extreme_exposures() {
    for (area, duration, liquid_fraction) in [
        (1e-4, 1.0, 0.9001600000000001),
        (1e100, 1e308, 0.999),
        (1e-300, 1e300, 0.999),
    ] {
        let catalog = fixture::catalog(0.0);
        let mut binding = fixture::binding(0.0);
        binding.accommodation_coefficient = 0.0;
        let chart = phase_fixture::curve(&catalog, &binding);
        let inventory = EquilibriumSublimation::new(&chart).unwrap();
        let initial = SolidSurfaceState {
            area_m2: area,
            ..phase_fixture::bulk(&chart, liquid_fraction)
        };
        let gas = fixture::gas(&catalog, 500.0, 0.0);
        let reflection = phase_fixture::accommodation(&catalog);
        let initial_solid = inventory.state(initial).unwrap().solid_mass_kg;
        // Only carrier A is present; its original constant Cp is 3.5 R.
        let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.02897;
        let heat_rate = 0.5 * 100000.0 / (2.0 * core::f64::consts::PI * r * 500.0).sqrt()
            * area
            * (3.0 * r)
            * (400.0 - 500.0);
        let expected_time = initial_solid * FUSION / -heat_rate;
        let candidate = advance_to_melting(
            inventory,
            &fixture::phase(),
            &binding,
            initial,
            &gas,
            &catalog,
            Some(&reflection),
            duration,
            fixture::controls(),
            || false,
        )
        .unwrap();
        assert!(candidate.boundary().melting_completed);
        fixture::near(
            candidate.boundary().elapsed_duration_s,
            expected_time,
            expected_time * 1e-10,
        );
        assert_eq!(
            candidate.boundary().remaining_duration_s,
            duration - candidate.boundary().elapsed_duration_s
        );
        assert_eq!(candidate.exchange().solid().mass_kg, initial.mass_kg);
        assert_eq!(candidate.exchange().gas().mass_kg(), gas.mass_kg());
        assert_eq!(candidate.exchange().step().transferred_mass_kg, 0.0);
        assert_eq!(
            candidate.exchange().solid().specific_enthalpy_j_kg,
            chart.knots()[3].specific_enthalpy_j_kg
        );
        assert_eq!(
            candidate.phase().final_state.equilibrium.temperature_k(),
            400.0
        );
        assert_eq!(candidate.phase().final_state.solid_mass_kg, 0.0);
        assert_eq!(
            candidate.phase().final_state.liquid_mass_kg,
            initial.mass_kg
        );
        fixture::near(
            candidate.phase().fusion_energy_j,
            -initial_solid * FUSION,
            1e-10,
        );
        fixture::near(
            candidate.reflection().unwrap().heat_into_gas_j,
            -initial_solid * FUSION,
            1e-10,
        );
        fixture::near(
            candidate.exchange().step().solid_energy_change_j
                + candidate.exchange().step().gas_energy_change_j,
            0.0,
            1e-7,
        );
        assert!(inventory.state(candidate.exchange().solid()).is_err());

        // The event root can round beyond this replay's request even though
        // the original represented update lands EXACTLY on the liquid knot.
        let replay = advance_to_melting(
            inventory,
            &fixture::phase(),
            &binding,
            initial,
            &gas,
            &catalog,
            Some(&reflection),
            candidate.boundary().elapsed_duration_s,
            fixture::controls(),
            || false,
        )
        .unwrap();
        assert!(replay.boundary().melting_completed);
        assert_eq!(
            replay.exchange().solid().specific_enthalpy_j_kg,
            candidate.exchange().solid().specific_enthalpy_j_kg
        );
        assert!(
            advance_equilibrium_surface(
                inventory,
                &fixture::phase(),
                &binding,
                initial,
                &gas,
                &catalog,
                Some(&reflection),
                duration,
                fixture::controls(),
                || false
            )
            .is_err()
        );
    }
}

#[test]
fn deposition_melts_initial_and_new_solid_before_the_full_request_exhausts_its_donor() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = phase_fixture::curve(&catalog, &binding);
    let inventory = EquilibriumSublimation::new(&chart).unwrap();
    let initial = phase_fixture::bulk(&chart, 0.999);
    let gas = fixture::gas(&catalog, 400.0, 100000.0);
    let initial_solid = inventory.state(initial).unwrap().solid_mass_kg;
    let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / MOLAR_MASS;
    let carried_above_solid = 200000.0 - 0.5 * r * 400.0;
    let expected_deposit = initial_solid * FUSION / (carried_above_solid - FUSION);
    let rate = 0.1 * (1000.0 - 100000.0) * initial.area_m2
        / (2.0 * core::f64::consts::PI * r * 400.0).sqrt();
    let candidate = advance_to_melting(
        inventory,
        &fixture::phase(),
        &binding,
        initial,
        &gas,
        &catalog,
        None,
        1.0,
        fixture::controls(),
        || false,
    )
    .unwrap();
    assert!(candidate.boundary().melting_completed);
    fixture::near(
        candidate.boundary().elapsed_duration_s,
        expected_deposit / -rate,
        1e-10,
    );
    fixture::near(
        candidate.exchange().step().transferred_mass_kg,
        -expected_deposit,
        1e-13,
    );
    fixture::near(
        candidate.exchange().solid().mass_kg,
        initial.mass_kg + expected_deposit,
        1e-13,
    );
    fixture::near(
        candidate.phase().freezing_mass_kg,
        -(initial_solid + expected_deposit),
        1e-13,
    );
    assert_eq!(candidate.phase().final_state.solid_mass_kg, 0.0);
    fixture::near(candidate.phase().solid_partition_residual_kg, 0.0, 1e-14);
    fixture::near(candidate.exchange().step().mass_residual_kg, 0.0, 1e-14);
    fixture::near(
        candidate.exchange().step().solid_energy_change_j
            + candidate.exchange().step().gas_energy_change_j,
        0.0,
        1e-7,
    );
    assert!(
        candidate.exchange().step().outgoing_mass_kg
            <= fixture::controls().maximum_gross_inventory_fraction * initial_solid
    );
    assert!(matches!(
        advance_equilibrium_surface(
            inventory,
            &fixture::phase(),
            &binding,
            initial,
            &gas,
            &catalog,
            None,
            1.0,
            fixture::controls(),
            || false
        ),
        Err(GasSurfaceError::InventoryLimit { .. })
    ));
}

#[test]
fn intervals_before_completion_reproduce_the_existing_scalar_endpoint() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = phase_fixture::curve(&catalog, &binding);
    let inventory = EquilibriumSublimation::new(&chart).unwrap();
    let initial = phase_fixture::bulk(&chart, 0.5);
    let reflection = phase_fixture::accommodation(&catalog);
    for (temperature, pressure, duration) in [
        (350.0, 500.0, 0.0),
        (350.0, 500.0, 0.001),
        (500.0, 500.0, 0.0001),
        (400.0, 1000.0, 10000.0),
    ] {
        let gas = fixture::gas(&catalog, temperature, pressure);
        let ordinary = advance_equilibrium_surface(
            inventory,
            &fixture::phase(),
            &binding,
            initial,
            &gas,
            &catalog,
            Some(&reflection),
            duration,
            fixture::controls(),
            || false,
        )
        .unwrap();
        let bounded = advance_to_melting(
            inventory,
            &fixture::phase(),
            &binding,
            initial,
            &gas,
            &catalog,
            Some(&reflection),
            duration,
            fixture::controls(),
            || false,
        )
        .unwrap();
        assert!(!bounded.boundary().melting_completed);
        assert_eq!(bounded.boundary().elapsed_duration_s, duration);
        assert_eq!(bounded.boundary().remaining_duration_s, 0.0);
        assert_eq!(bounded.endpoint(), &ordinary);
    }
}

#[test]
fn event_keeps_original_donor_source_and_cancellation_refusals_atomic() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = phase_fixture::curve(&catalog, &binding);
    let inventory = EquilibriumSublimation::new(&chart).unwrap();
    let initial = phase_fixture::bulk(&chart, 0.999);
    let gas = fixture::gas(&catalog, 400.0, 100000.0);
    let saved = gas.clone();
    let mut tight = fixture::controls();
    tight.maximum_gross_inventory_fraction = 0.01;
    assert!(matches!(
        advance_to_melting(
            inventory,
            &fixture::phase(),
            &binding,
            initial,
            &gas,
            &catalog,
            None,
            1.0,
            tight,
            || false
        ),
        Err(GasSurfaceError::InventoryLimit { .. })
    ));
    let wrong = GasSurfaceBinding {
        solid_enthalpy_offset_j_kg: 1000.0,
        ..binding.clone()
    };
    assert!(matches!(
        advance_to_melting(
            inventory,
            &fixture::phase(),
            &wrong,
            initial,
            &gas,
            &catalog,
            None,
            1.0,
            fixture::controls(),
            || false
        ),
        Err(GasSurfaceError::CaloricMismatch { .. })
    ));
    assert!(
        advance_to_melting(
            inventory,
            &fixture::phase(),
            &binding,
            phase_fixture::bulk(&chart, 1.0),
            &gas,
            &catalog,
            None,
            1.0,
            fixture::controls(),
            || false
        )
        .is_err()
    );
    let mut calls = 0;
    advance_to_melting(
        inventory,
        &fixture::phase(),
        &binding,
        initial,
        &gas,
        &catalog,
        None,
        1.0,
        fixture::controls(),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    // Each stop asserts the one cancellation that its first firing poll gives. With no
    // accommodation law the polls run in this order. First boundary.rs:117, the probe at
    // boundary.rs:138, the event re-evaluation at boundary.rs:138 when the event shortens the
    // interval, and equilibrium.rs:262: b <= 4 polls through gas.rs:204-209, each refusing
    // with `Surface(SublimationError::Cancelled)`. Then the gas owner,
    // `exchange_species_boundary_relative` in transfer/species/boundary.rs: :447 and :456,
    // then :478 and :484 on its idle path, or :511, :557, :569, its caloric reconstruction,
    // :609 per species and :627, so g >= 4 polls, each refusing with
    // `GasEnergyError::Cancelled`. Last, gas.rs:920 in `finish_surface` and boundary.rs:260.
    // Stop 1 fires at boundary.rs:117 and stop `calls` at boundary.rs:260. Since g >= b, stop
    // calls / 2 = (b + g + 2) / 2 lies in (b, b + g], inside the gas owner, whose refusal
    // reaches this owner through `From<GasEnergyError> for GasSurfaceError` (gas.rs:117-121).
    // The sibling sublimation::gas::tests, over the same three stops, accepts either owner
    // (sublimation/gas/tests.rs:648-652):
    //     Err(GasSurfaceError::Surface(SublimationError::Cancelled))
    //         | Err(GasSurfaceError::Gas(GasEnergyError::Cancelled))
    for (stop, gas_owner_fires) in [(1, false), (calls / 2, true), (calls, false)] {
        let mut count = 0;
        let result = advance_to_melting(
            inventory,
            &fixture::phase(),
            &binding,
            initial,
            &gas,
            &catalog,
            None,
            1.0,
            fixture::controls(),
            || {
                count += 1;
                count == stop
            },
        );
        if gas_owner_fires {
            assert!(
                matches!(result, Err(GasSurfaceError::Gas(GasEnergyError::Cancelled))),
                "stop {stop} of {calls}: {result:?}"
            );
        } else {
            assert!(
                matches!(
                    result,
                    Err(GasSurfaceError::Surface(SublimationError::Cancelled))
                ),
                "stop {stop} of {calls}: {result:?}"
            );
        }
        assert_eq!(gas, saved);
        assert_eq!(initial, phase_fixture::bulk(&chart, 0.999));
    }
}

#[test]
fn melting_time_retains_admitted_alignment_and_large_formation_reference_shifts() {
    let base = fixture::catalog(0.0);
    let chart = phase_fixture::curve(&base, &fixture::binding(0.0));
    let inventory = EquilibriumSublimation::new(&chart).unwrap();
    let initial = phase_fixture::bulk(&chart, 0.999);
    let initial_solid = inventory.state(initial).unwrap().solid_mass_kg;
    let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / MOLAR_MASS;
    let rate = 0.1 * (1000.0 - 100000.0) * initial.area_m2
        / (2.0 * core::f64::consts::PI * r * 400.0).sqrt();
    for (reference, misalignment) in [(0.0, 0.0), (0.0, 50.0), (-1e14, 0.0), (1e14, 0.0)] {
        let catalog = fixture::catalog(reference);
        let binding = fixture::binding(reference + misalignment);
        let gas = fixture::gas(&catalog, 400.0, 100000.0);
        let mut controls = fixture::controls();
        controls.latent_compatibility_tolerance_j_kg = 64.0;
        let alignment = (catalog
            .species_specific_enthalpy_j_kg(&fixture::species(), 400.0)
            .unwrap()
            - (chart.knots()[1].specific_enthalpy_j_kg + binding.solid_enthalpy_offset_j_kg))
            - 200000.0;
        let effective_latent = 200000.0 + alignment - 0.5 * r * 400.0;
        let expected_mass = initial_solid * FUSION / (effective_latent - FUSION);
        let candidate = advance_to_melting(
            inventory,
            &fixture::phase(),
            &binding,
            initial,
            &gas,
            &catalog,
            None,
            1.0,
            controls,
            || false,
        )
        .unwrap();
        assert!(candidate.boundary().melting_completed);
        fixture::near(
            candidate.boundary().elapsed_duration_s,
            expected_mass / -rate,
            1e-10,
        );
        fixture::near(
            candidate.exchange().step().transferred_mass_kg,
            -expected_mass,
            1e-13,
        );
        assert_eq!(candidate.inventory().curve().identity(), chart.identity());
        assert_eq!(
            candidate.phase().final_state.equilibrium.temperature_k(),
            400.0
        );
        assert!(
            candidate.boundary().condensed_energy_residual_j.abs()
                + candidate
                    .exchange()
                    .gas_exchange()
                    .step()
                    .energy_residual_j
                    .abs()
                <= controls.gas.energy_tolerance_j
        );
    }
}

#[test]
fn explicit_liquid_handoff_keeps_the_same_chart_mass_gas_and_energy_datum() {
    use crate::evaporation::gas::{GasEvaporationSolve, GasLiquidBinding, equilibrium as liquid};
    use fs_thermal::evaporation::{
        EvaporationError, EvaporationLaw, EvaporationProperties, LiquidSurfaceState,
    };

    // Separate, explicitly chosen continuum law. Its latent heat is h_v-h_l,
    // whereas the molecular source uses h_v-h_s on this SAME original chart.
    struct LiquidPhase;
    impl EvaporationLaw for LiquidPhase {
        fn temperature_range_k(&self) -> [f64; 2] {
            [200.0, 900.0]
        }
        fn properties(&self, t: f64) -> Result<EvaporationProperties, EvaporationError> {
            let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / MOLAR_MASS;
            let pressure = 1000.0 * fs_math::det::exp(FUSION / r * (1.0 / 400.0 - 1.0 / t));
            Ok(EvaporationProperties {
                saturation_pressure_pa: pressure,
                saturation_pressure_derivative_pa_k: pressure * FUSION / (r * t * t),
                latent_heat_j_kg: FUSION,
                liquid_density_kg_m3: 917.0,
                molar_mass_kg_mol: MOLAR_MASS,
            })
        }
    }

    let catalog = fixture::catalog(0.0);
    let mut binding = fixture::binding(0.0);
    binding.accommodation_coefficient = 0.0;
    let chart = phase_fixture::curve(&catalog, &binding);
    let inventory = EquilibriumSublimation::new(&chart).unwrap();
    let initial = phase_fixture::bulk(&chart, 0.999);
    let gas = fixture::gas(&catalog, 500.0, 0.0);
    let reflection = phase_fixture::accommodation(&catalog);
    let melted = advance_to_melting(
        inventory,
        &fixture::phase(),
        &binding,
        initial,
        &gas,
        &catalog,
        Some(&reflection),
        1.0,
        fixture::controls(),
        || false,
    )
    .unwrap();
    let endpoint = melted.exchange().solid();
    let liquid_state = LiquidSurfaceState {
        mass_kg: endpoint.mass_kg,
        specific_enthalpy_j_kg: endpoint.specific_enthalpy_j_kg,
        area_m2: initial.area_m2,
    };
    let liquid_binding = GasLiquidBinding {
        species: binding.species.clone(),
        liquid_enthalpy_offset_j_kg: binding.solid_enthalpy_offset_j_kg,
    };
    let controls = GasEvaporationSolve {
        gas: fixture::controls().gas,
        maximum_inventory_fraction: 0.25,
        maximum_temperature_increment_k: 5.0,
        latent_compatibility_tolerance_j_kg: 1e-5,
    };
    let continued = liquid::advance_equilibrium_surface(
        liquid::EquilibriumEvaporation::new(&chart).unwrap(),
        &LiquidPhase,
        &liquid_binding,
        liquid_state,
        melted.exchange().gas(),
        &catalog,
        0.001,
        0.01,
        controls,
        || false,
    )
    .unwrap();
    assert!(continued.exchange().step().transferred_mass_kg > 0.0);
    assert!(continued.phase().freezing_mass_kg > 0.0);
    assert_eq!(continued.inventory().curve().identity(), chart.identity());
    fixture::near(
        (continued.exchange().liquid().mass_kg - initial.mass_kg)
            + (continued.exchange().gas().mass_kg() - gas.mass_kg()),
        0.0,
        1e-14,
    );
    let first = melted.exchange().step();
    let second = continued.exchange().step();
    fixture::near(
        first.solid_energy_change_j
            + first.gas_energy_change_j
            + second.liquid_enthalpy_change_j
            + second.gas_internal_energy_change_j,
        0.0,
        1e-7,
    );
    assert_eq!(
        continued.phase().final_state.equilibrium.temperature_k(),
        400.0
    );
}
