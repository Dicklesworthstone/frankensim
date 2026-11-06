//! G1 analytic phase/energy balances and G3 original-source/datum comparisons.
use super::super::{
    advance_surface_with_gas,
    tests::{Phase, binding, controls, gas},
};
use super::*;
use crate::phase::EnthalpyPhaseKnot;
use crate::thermochemical_gas::tests::{species_with_energy_reference, transport};
use fs_blake3::ContentHash;
use fs_thermal::evaporation::{EvaporationError, EvaporationProperties};
use fs_thermochem::{Composition, CompositionBasis, mixture::FrozenIdealGasMixtureModelV1};

const MOLAR_MASS: f64 = 0.018015;
const R: f64 = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / MOLAR_MASS;
pub(in crate::evaporation::gas) const CP: f64 = 4.0 * R;
pub(in crate::evaporation::gas) const MELTING: f64 = 300.0;
pub(in crate::evaporation::gas) const FUSION: f64 = 100_000.0;
pub(in crate::evaporation::gas) const VAPORIZATION: f64 = 200_000.0;
pub(in crate::evaporation::gas) const TRANSFER: f64 = 0.02;

pub(in crate::evaporation::gas) fn catalog(reference: f64) -> ThermochemicalGasModel {
    ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![
                species_with_energy_reference("A", 0.02897, 3.5, 0.0, 0.0),
                species_with_energy_reference("B", MOLAR_MASS, 4.0, 0.0, reference),
            ],
            Composition::new(CompositionBasis::MoleFraction, vec![0.99, 0.01]).unwrap(),
        )
        .unwrap(),
        vec![transport("A", 1.7e-5, 0.025), transport("B", 2e-5, 0.03)],
        Some(&binding(reference).species),
    )
    .unwrap()
}

pub(in crate::evaporation::gas) fn curve() -> EquilibriumEnthalpyPhaseCurve {
    let liquid_h = CP * MELTING - VAPORIZATION;
    let solid_h = liquid_h - FUSION;
    let knot = |h, t, fraction: f64| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h,
        temperature_k: t,
        liquid_mass_fraction: fraction,
        bulk_density_kg_m3: 1.0 / ((1.0 - fraction) / 900.0 + fraction / 1000.0),
    };
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0x72; 32]),
        vec![
            knot(solid_h - 2.0 * R * 100.0, 200.0, 0.0),
            knot(solid_h, MELTING, 0.0),
            knot(solid_h + 0.5 * FUSION, MELTING, 0.5),
            knot(liquid_h, MELTING, 1.0),
            knot(liquid_h + CP * 600.0, 900.0, 1.0),
        ],
    )
    .unwrap()
}

pub(in crate::evaporation::gas) fn bulk(
    curve: &EquilibriumEnthalpyPhaseCurve,
    fraction: f64,
) -> LiquidSurfaceState {
    LiquidSurfaceState {
        mass_kg: 0.01,
        area_m2: 1e-3,
        specific_enthalpy_j_kg: curve.knots()[1].specific_enthalpy_j_kg + fraction * FUSION,
    }
}

pub(in crate::evaporation::gas) fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a} != {b}; tolerance {tolerance}"
    );
}

fn advance<'curve, 'catalog>(
    curve: &'curve EquilibriumEnthalpyPhaseCurve,
    bulk: LiquidSurfaceState,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    duration: f64,
) -> EquilibriumSurfaceCandidate<'curve, 'catalog> {
    advance_equilibrium_surface(
        EquilibriumEvaporation::new(curve).unwrap(),
        &Phase,
        &binding(0.0),
        bulk,
        gas,
        catalog,
        TRANSFER,
        duration,
        controls(),
        || false,
    )
    .unwrap()
}

#[test]
fn evaporation_freezes_liquid_on_the_original_plateau_and_conserves_total_energy() {
    let catalog = catalog(0.0);
    let chart = curve();
    let initial = bulk(&chart, 0.5);
    let gas = gas(&catalog, MELTING, 0.0, 0.01);
    let result = advance(&chart, initial, &gas, &catalog, 1.0);
    let phase = result.phase();
    let step = result.exchange().step();
    // Independently integrated Stefan film, then the isothermal constituent
    // balance: evaporating dm requires another dm*L_v/L_f to freeze.
    let expected_dm = initial.area_m2 * TRANSFER * 100_000.0 / (R * MELTING)
        * (100_000.0_f64 / (100_000.0 - 1500.0)).ln();
    let frozen = expected_dm * VAPORIZATION / FUSION;
    near(step.transferred_mass_kg, expected_dm, 1e-18);
    near(phase.freezing_mass_kg, frozen, 2e-17);
    near(phase.final_state.solid_mass_kg, 0.005 + frozen, 2e-17);
    near(
        phase.final_state.liquid_mass_kg,
        0.005 - expected_dm - frozen,
        2e-17,
    );
    assert_eq!(step.surface_temperature_k, [MELTING, MELTING]);
    near(phase.fusion_energy_j, expected_dm * VAPORIZATION, 2e-12);
    near(step.latent_energy_j, expected_dm * VAPORIZATION, 1e-12);
    near(
        phase.liquid_partition_residual_kg,
        0.0,
        controls().gas.mass_tolerance_kg,
    );
    let final_bulk = result.exchange().liquid();
    let original_total = initial.mass_kg * initial.specific_enthalpy_j_kg
        + gas.mass_kg() * gas.state().specific_internal_energy_j_kg();
    let final_total = final_bulk.mass_kg * final_bulk.specific_enthalpy_j_kg
        + result.exchange().gas().mass_kg()
            * result
                .exchange()
                .gas()
                .state()
                .specific_internal_energy_j_kg();
    near(final_total, original_total, 1e-7);
    near(
        step.energy_residual_j,
        0.0,
        controls().gas.energy_tolerance_j,
    );
    assert_eq!(
        phase.final_state.equilibrium.phase_curve_identity(),
        chart.identity()
    );
    assert!(std::ptr::eq(result.exchange().catalog(), &catalog));
    assert_ne!(phase.initial.equilibrium.bulk_density_kg_m3(), 1000.0);
}

#[test]
fn hot_vapor_condensation_melts_solid_with_both_latent_and_upstream_sensible_heat() {
    let catalog = catalog(0.0);
    let chart = curve();
    let initial = bulk(&chart, 0.5);
    let gas = gas(&catalog, 320.0, 5000.0, 0.01);
    let result = advance(&chart, initial, &gas, &catalog, 1.0);
    let step = result.exchange().step();
    let phase = result.phase();
    assert!(step.transferred_mass_kg < 0.0);
    let frozen = step.transferred_mass_kg * (VAPORIZATION + CP * 20.0) / FUSION;
    assert!(frozen < 0.0);
    near(phase.freezing_mass_kg, frozen, 2e-17);
    near(phase.final_state.solid_mass_kg, 0.005 + frozen, 2e-17);
    near(
        phase.final_state.liquid_mass_kg,
        0.005 - step.transferred_mass_kg - frozen,
        2e-17,
    );
    assert_eq!(step.surface_temperature_k, [MELTING, MELTING]);
    near(
        step.energy_into_gas_j,
        step.transferred_mass_kg * CP * 320.0,
        1e-9,
    );
    near(
        step.liquid_enthalpy_change_j + step.gas_internal_energy_change_j,
        0.0,
        1e-7,
    );
    near(step.mass_residual_kg, 0.0, controls().gas.mass_tolerance_kg);
}

#[test]
fn pure_liquid_cooling_can_enter_fusion_without_temperature_clipping() {
    let catalog = catalog(0.0);
    let chart = curve();
    let mut initial = bulk(&chart, 1.0);
    initial.specific_enthalpy_j_kg += CP * 0.001;
    let gas = gas(&catalog, MELTING, 0.0, 0.01);
    let result = advance(&chart, initial, &gas, &catalog, 1.0);
    assert_eq!(result.phase().initial.solid_mass_kg, 0.0);
    assert!(result.phase().final_state.solid_mass_kg > 0.0);
    assert_eq!(
        result.phase().final_state.equilibrium.temperature_k(),
        MELTING
    );
    let dm = result.exchange().step().transferred_mass_kg;
    let sensible_release = (initial.mass_kg - dm) * CP * 0.001;
    near(
        result.phase().fusion_energy_j,
        dm * VAPORIZATION - sensible_release,
        2e-12,
    );
    near(
        result.exchange().step().energy_residual_j,
        0.0,
        controls().gas.energy_tolerance_j,
    );
}

#[test]
fn complete_melting_continues_on_the_original_liquid_branch_and_matches_legacy() {
    let catalog = catalog(0.0);
    let chart = curve();
    let initial = bulk(&chart, 1.0 - 1e-6);
    let gas = gas(&catalog, MELTING, 5000.0, 0.01);
    let melted = advance(&chart, initial, &gas, &catalog, 1.0);
    assert_eq!(melted.phase().final_state.solid_mass_kg, 0.0);
    assert!(melted.phase().final_state.equilibrium.temperature_k() > MELTING);
    near(
        melted.phase().freezing_mass_kg,
        -melted.phase().initial.solid_mass_kg,
        1e-18,
    );
    let mixed = advance(
        &chart,
        melted.exchange().liquid(),
        melted.exchange().gas(),
        &catalog,
        0.1,
    );
    let pure = advance_surface_with_gas(
        &Phase,
        &binding(0.0),
        melted.exchange().liquid(),
        melted.exchange().gas(),
        &catalog,
        |h| {
            chart
                .temperature_at_specific_enthalpy(h)
                .map_err(|_| EvaporationError::Invalid("original liquid caloric source"))
        },
        TRANSFER,
        0.1,
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(mixed.exchange(), &pure);
    assert_eq!(mixed.phase().freezing_mass_kg, 0.0);
}

#[test]
fn old_solid_inventory_cannot_finance_liquid_evaporation() {
    let catalog = catalog(0.0);
    let chart = curve();
    let initial = bulk(&chart, 1e-6);
    let gas = gas(&catalog, MELTING, 0.0, 0.01);
    let refused = advance_equilibrium_surface(
        EquilibriumEvaporation::new(&chart).unwrap(),
        &Phase,
        &binding(0.0),
        initial,
        &gas,
        &catalog,
        TRANSFER,
        1.0,
        controls(),
        || false,
    );
    let Err(GasEvaporationError::StepSizeLimit {
        maximum_duration_s, ..
    }) = refused
    else {
        panic!("evaporation must bound the original liquid inventory before its caloric update");
    };
    assert!(maximum_duration_s > 0.0 && maximum_duration_s < 0.1);
    let result = advance(&chart, initial, &gas, &catalog, 0.99 * maximum_duration_s);
    assert!(
        result.exchange().step().transferred_mass_kg
            <= controls().maximum_inventory_fraction * result.phase().initial.liquid_mass_kg
    );
    assert!(result.phase().final_state.liquid_mass_kg > 0.0);
}

#[test]
fn large_gas_datum_keeps_fusion_and_source_alignment_in_the_original_balance() {
    let chart = curve();
    let initial = bulk(&chart, 0.5);
    let inventory = EquilibriumEvaporation::new(&chart).unwrap();
    let mut baseline = None;
    for reference in [0.0, -1e14, 1e14] {
        let catalog = catalog(reference);
        let binding = binding(reference);
        let gas = gas(&catalog, MELTING, 0.0, 0.01);
        let mut budget = controls();
        budget.latent_compatibility_tolerance_j_kg = 0.125;
        let result = advance_equilibrium_surface(
            inventory,
            &Phase,
            &binding,
            initial,
            &gas,
            &catalog,
            TRANSFER,
            1.0,
            budget,
            || false,
        )
        .unwrap();
        let step = result.exchange().step();
        let alignment = (catalog
            .species_specific_enthalpy_j_kg(&binding.species, MELTING)
            .unwrap()
            - (chart.knots()[3].specific_enthalpy_j_kg + reference))
            - VAPORIZATION;
        let frozen = step.transferred_mass_kg * (VAPORIZATION + alignment) / FUSION;
        near(result.phase().freezing_mass_kg, frozen, 2e-17);
        near(step.energy_residual_j, 0.0, budget.gas.energy_tolerance_j);
        let state = (step.transferred_mass_kg, step.gas_temperature_k[1]);
        if let Some((mass, temperature)) = baseline {
            assert_eq!(state.0, mass);
            near(state.1, temperature, 1e-7);
        } else {
            baseline = Some(state);
        }
    }
}

#[test]
fn absent_liquid_mushy_charts_and_incompatible_constituents_refuse() {
    let chart = curve();
    let inventory = EquilibriumEvaporation::new(&chart).unwrap();
    assert!(inventory.state(bulk(&chart, 0.0)).is_err());
    let mut knots = chart.knots().to_vec();
    knots[2].temperature_k += 1.0;
    knots[3].temperature_k += 2.0;
    let mushy = EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x73; 32]), knots).unwrap();
    assert!(EquilibriumEvaporation::new(&mushy).is_err());

    struct WrongDensity;
    impl EvaporationLaw for WrongDensity {
        fn temperature_range_k(&self) -> [f64; 2] {
            Phase.temperature_range_k()
        }
        fn properties(&self, t: f64) -> Result<EvaporationProperties, EvaporationError> {
            Ok(EvaporationProperties {
                liquid_density_kg_m3: 900.0,
                ..Phase.properties(t)?
            })
        }
    }
    let catalog = catalog(0.0);
    let gas = gas(&catalog, MELTING, 0.0, 0.01);
    let initial = bulk(&chart, 0.5);
    assert!(matches!(
        advance_equilibrium_surface(
            inventory,
            &WrongDensity,
            &binding(0.0),
            initial,
            &gas,
            &catalog,
            TRANSFER,
            1.0,
            controls(),
            || false,
        ),
        Err(GasEvaporationError::Invalid(_))
    ));
    assert!(matches!(
        advance_equilibrium_surface(
            inventory,
            &Phase,
            &binding(1000.0),
            initial,
            &gas,
            &catalog,
            TRANSFER,
            1.0,
            controls(),
            || false,
        ),
        Err(GasEvaporationError::CaloricMismatch { .. })
    ));
}

#[test]
fn sealed_surface_and_zero_duration_keep_the_same_phase_and_gas() {
    let catalog = catalog(0.0);
    let chart = curve();
    let initial = bulk(&chart, 0.5);
    let gas = gas(&catalog, MELTING, 0.0, 0.01);
    for (coefficient, duration) in [(0.0, 1.0), (TRANSFER, 0.0)] {
        let result = advance_equilibrium_surface(
            EquilibriumEvaporation::new(&chart).unwrap(),
            &Phase,
            &binding(0.0),
            initial,
            &gas,
            &catalog,
            coefficient,
            duration,
            controls(),
            || false,
        )
        .unwrap();
        assert_eq!(result.exchange().liquid(), initial);
        assert_eq!(result.exchange().gas(), &gas);
        assert_eq!(result.phase().initial, result.phase().final_state);
        assert_eq!(result.phase().fusion_energy_j, 0.0);
    }
}

#[test]
fn late_cancellation_and_complete_freezing_publish_no_partial_endpoint() {
    let catalog = catalog(0.0);
    let chart = curve();
    let inventory = EquilibriumEvaporation::new(&chart).unwrap();
    let initial = bulk(&chart, 0.5);
    let gas = gas(&catalog, MELTING, 0.0, 0.01);
    let saved = gas.clone();
    let mut calls = 0;
    advance_equilibrium_surface(
        inventory,
        &Phase,
        &binding(0.0),
        initial,
        &gas,
        &catalog,
        TRANSFER,
        1.0,
        controls(),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    assert!(calls > 1);
    let mut attempt = 0;
    assert!(matches!(
        advance_equilibrium_surface(
            inventory,
            &Phase,
            &binding(0.0),
            initial,
            &gas,
            &catalog,
            TRANSFER,
            1.0,
            controls(),
            || {
                attempt += 1;
                attempt == calls
            },
        ),
        Err(GasEvaporationError::Cancelled)
    ));

    // Looser inventory admission cannot turn the remaining liquid into a solid
    // interface implicitly. An event owner must locate freezing completion.
    let tiny_liquid = bulk(&chart, 1e-6);
    let mut budget = controls();
    budget.maximum_inventory_fraction = 0.9;
    assert!(
        advance_equilibrium_surface(
            inventory,
            &Phase,
            &binding(0.0),
            tiny_liquid,
            &gas,
            &catalog,
            TRANSFER,
            0.02,
            budget,
            || false,
        )
        .is_err()
    );
    assert_eq!(gas, saved);
    assert_eq!(initial, bulk(&chart, 0.5));
}
