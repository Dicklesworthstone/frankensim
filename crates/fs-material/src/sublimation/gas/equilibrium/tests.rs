//! Manufactured sharp-phase charts with independent mass/energy references.
use super::*;
use crate::sublimation::gas::{advance_surface_with_gas, tests as fixture, SublimationError};
use fs_blake3::ContentHash;
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K};

const FUSION: f64 = 100000.0;

pub(super) fn curve(catalog: &ThermochemicalGasModel, binding: &GasSurfaceBinding) -> EquilibriumEnthalpyPhaseCurve {
    let (h, cp) = fixture::chart(catalog, binding);
    let knot = |h, t, fraction| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h, temperature_k: t,
        liquid_mass_fraction: fraction, bulk_density_kg_m3: 917.0,
    };
    EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x61; 32]), vec![
        knot(h - 200.0 * cp, 200.0, 0.0), knot(h, 400.0, 0.0),
        knot(h + 0.5 * FUSION, 400.0, 0.5), knot(h + FUSION, 400.0, 1.0),
        knot(h + FUSION + 500.0 * cp, 900.0, 1.0),
    ]).unwrap()
}
pub(super) fn bulk(curve: &EquilibriumEnthalpyPhaseCurve, fraction: f64) -> SolidSurfaceState {
    SolidSurfaceState {
        mass_kg: 0.01, area_m2: 1e-4,
        specific_enthalpy_j_kg: curve.knots()[1].specific_enthalpy_j_kg + fraction * FUSION,
    }
}
pub(super) fn accommodation(catalog: &ThermochemicalGasModel) -> ThermalAccommodation {
    ThermalAccommodation::new(catalog, &[
        (SpeciesId::new("A").unwrap(), 0.5), (fixture::species(), 0.5),
    ], [200.0, 900.0]).unwrap()
}

#[test]
fn sublimation_freezes_liquid_at_constant_temperature_and_accounts_fusion_once() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = curve(&catalog, &binding);
    let inventory = EquilibriumSublimation::new(&chart).unwrap();
    let initial = bulk(&chart, 0.5);
    let gas = fixture::gas(&catalog, 400.0, 0.0);
    let dt = 0.01;
    let result = advance_equilibrium_surface(inventory, &fixture::phase(), &binding,
        initial, &gas, &catalog, None, dt, fixture::controls(), || false).unwrap();
    let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.018015;
    let dm = 0.1 * 1000.0 / (2.0 * core::f64::consts::PI * r * 400.0).sqrt() * initial.area_m2 * dt;
    let solid_h = chart.knots()[1].specific_enthalpy_j_kg;
    let molecular_h = solid_h + 200000.0 - 0.5 * r * 400.0;
    let frozen = dm * (molecular_h - solid_h) / FUSION;
    let phase = result.phase();
    assert_eq!(phase.initial.equilibrium.temperature_k(), 400.0);
    assert_eq!(phase.final_state.equilibrium.temperature_k(), 400.0);
    fixture::near(result.exchange().step().transferred_mass_kg, dm, 1e-18);
    fixture::near(phase.freezing_mass_kg, frozen, 1e-15);
    fixture::near(phase.final_state.liquid_mass_kg, 0.005 - frozen, 1e-15);
    fixture::near(phase.final_state.solid_mass_kg, 0.005 - dm + frozen, 1e-15);
    fixture::near(result.exchange().step().energy_into_gas_j, dm * molecular_h, 1e-10);
    fixture::near(phase.fusion_energy_j, frozen * FUSION, 1e-10);
    fixture::near(result.exchange().step().energy_residual_j, 0.0, 1e-7);
    fixture::near(phase.solid_partition_residual_kg, 0.0, 1e-14);
    assert_eq!(result.inventory().curve().identity(), chart.identity());
    assert!(result.exchange().gas().species_mass_kg(&fixture::species()).unwrap() > 0.0);
}

#[test]
fn deposition_melts_liquid_inventory_without_becoming_a_second_latent_source() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = curve(&catalog, &binding);
    let inventory = EquilibriumSublimation::new(&chart).unwrap();
    let initial = bulk(&chart, 0.5);
    let gas = fixture::gas(&catalog, 400.0, 5000.0);
    let result = advance_equilibrium_surface(inventory, &fixture::phase(), &binding,
        initial, &gas, &catalog, None, 0.01, fixture::controls(), || false).unwrap();
    let step = result.exchange().step();
    let solid_h = chart.knots()[1].specific_enthalpy_j_kg;
    let expected = step.transferred_mass_kg * (step.outgoing_specific_energy_j_kg - solid_h) / FUSION;
    assert!(step.transferred_mass_kg < 0.0);
    assert!(result.phase().freezing_mass_kg < 0.0);
    assert_eq!(result.phase().final_state.equilibrium.temperature_k(), 400.0);
    fixture::near(result.phase().freezing_mass_kg, expected, 1e-15);
    fixture::near(step.solid_energy_change_j + step.gas_energy_change_j, 0.0, 1e-7);
    fixture::near(step.mass_residual_kg, 0.0, 1e-14);
}

#[test]
fn reflection_changes_phase_fraction_even_with_zero_mass_and_zero_temperature_change() {
    let catalog = fixture::catalog(0.0);
    let mut binding = fixture::binding(0.0);
    binding.accommodation_coefficient = 0.0;
    let chart = curve(&catalog, &binding);
    let inventory = EquilibriumSublimation::new(&chart).unwrap();
    let initial = bulk(&chart, 0.5);
    let gas = fixture::gas(&catalog, 350.0, 0.0);
    let coefficients = accommodation(&catalog);
    let result = advance_equilibrium_surface(inventory, &fixture::phase(), &binding,
        initial, &gas, &catalog, Some(&coefficients), 0.001, fixture::controls(), || false).unwrap();
    let q = result.reflection().unwrap().heat_into_gas_j;
    assert!(q > 0.0);
    assert_eq!(result.exchange().solid().mass_kg, initial.mass_kg);
    assert_eq!(result.exchange().step().transferred_mass_kg, 0.0);
    assert_eq!(result.phase().final_state.equilibrium.temperature_k(), 400.0);
    fixture::near(result.phase().freezing_mass_kg, q / FUSION, 1e-15);
    fixture::near(result.phase().fusion_energy_j, q, 1e-10);
    fixture::near(result.exchange().step().energy_residual_j, 0.0, 1e-7);
    assert_eq!(result.exchange().gas().mass_kg(), gas.mass_kg());
}

#[test]
fn initially_liquid_mass_cannot_fund_gross_solid_emission() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = curve(&catalog, &binding);
    let inventory = EquilibriumSublimation::new(&chart).unwrap();
    let initial = bulk(&chart, 0.99999);
    let gas = fixture::gas(&catalog, 400.0, 0.0);
    let refused = advance_equilibrium_surface(inventory, &fixture::phase(), &binding,
        initial, &gas, &catalog, None, 1.0, fixture::controls(), || false);
    let Err(GasSurfaceError::InventoryLimit { maximum_duration_s }) = refused else {
        panic!("gross withdrawal must use old solid mass, not the complete condensed mass");
    };
    assert!(maximum_duration_s > 0.0 && maximum_duration_s < 0.01);
    let accepted = advance_equilibrium_surface(inventory, &fixture::phase(), &binding,
        initial, &gas, &catalog, None, 0.99 * maximum_duration_s, fixture::controls(), || false).unwrap();
    assert!(accepted.exchange().step().outgoing_mass_kg
        <= fixture::controls().maximum_gross_inventory_fraction * accepted.phase().initial.solid_mass_kg);
}

#[test]
fn latent_plateau_can_finish_freezing_and_continue_on_the_original_solid_branch() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = curve(&catalog, &binding);
    let inventory = EquilibriumSublimation::new(&chart).unwrap();
    let initial = bulk(&chart, 1e-6);
    let gas = fixture::gas(&catalog, 400.0, 0.0);
    let first = advance_equilibrium_surface(inventory, &fixture::phase(), &binding,
        initial, &gas, &catalog, None, 0.1, fixture::controls(), || false).unwrap();
    assert!(first.phase().final_state.equilibrium.is_fully_solid());
    assert!(first.phase().final_state.equilibrium.temperature_k() < 400.0);
    fixture::near(first.phase().freezing_mass_kg, first.phase().initial.liquid_mass_kg, 1e-15);
    let mixed = advance_equilibrium_surface(inventory, &fixture::phase(), &binding,
        first.exchange().solid(), first.exchange().gas(), &catalog, None,
        0.01, fixture::controls(), || false).unwrap();
    let pure = advance_surface_with_gas(&fixture::phase(), &binding, first.exchange().solid(),
        first.exchange().gas(), &catalog,
        |h| chart.temperature_at_specific_enthalpy(h).map_err(|_| SublimationError::Invalid("source chart")),
        0.01, fixture::controls(), || false).unwrap();
    assert_eq!(mixed.exchange(), &pure);
}

#[test]
fn consistent_energy_reference_changes_preserve_phase_mass_and_temperature() {
    let mut outputs = Vec::new();
    for offset in [0.0, 1e7, -1e7] {
        let catalog = fixture::catalog(offset);
        let binding = fixture::binding(offset);
        let chart = curve(&catalog, &binding);
        let gas = fixture::gas(&catalog, 350.0, 500.0);
        let result = advance_equilibrium_surface(EquilibriumSublimation::new(&chart).unwrap(),
            &fixture::phase(), &binding, bulk(&chart, 0.5), &gas, &catalog, None,
            0.01, fixture::controls(), || false).unwrap();
        outputs.push((result.phase().final_state.solid_mass_kg,
            result.phase().final_state.liquid_mass_kg, result.exchange().gas().state().state().temperature));
    }
    for value in &outputs[1..] {
        fixture::near(value.0, outputs[0].0, 1e-14);
        fixture::near(value.1, outputs[0].1, 1e-14);
        fixture::near(value.2, outputs[0].2, 1e-7);
    }
}

#[test]
fn unsupported_mushy_charts_and_complete_melting_do_not_change_accepted_owners() {
    let catalog = fixture::catalog(0.0);
    let mut binding = fixture::binding(0.0);
    binding.accommodation_coefficient = 0.0;
    let chart = curve(&catalog, &binding);
    let mut knots = chart.knots().to_vec();
    knots[2].temperature_k = 410.0;
    knots[3].temperature_k = 420.0;
    let mushy = EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x62; 32]), knots).unwrap();
    assert!(EquilibriumSublimation::new(&mushy).is_err());
    let mut knots = chart.knots().to_vec();
    knots[2].liquid_mass_fraction = 0.4;
    let wrong = EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x63; 32]), knots).unwrap();
    assert!(EquilibriumSublimation::new(&wrong).is_err());
    let inventory = EquilibriumSublimation::new(&chart).unwrap();
    assert!(inventory.state(bulk(&chart, 1.0)).is_err());
    let initial = bulk(&chart, 0.999999);
    let gas = fixture::gas(&catalog, 500.0, 0.0);
    let saved = gas.clone();
    let coefficients = accommodation(&catalog);
    assert!(advance_equilibrium_surface(inventory, &fixture::phase(), &binding,
        initial, &gas, &catalog, Some(&coefficients), 0.001, fixture::controls(), || false).is_err());
    assert_eq!(gas, saved);
    assert_eq!(initial, bulk(&chart, 0.999999));
}

#[test]
fn late_cancellation_and_inconsistent_latent_references_publish_no_partition() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = curve(&catalog, &binding);
    let inventory = EquilibriumSublimation::new(&chart).unwrap();
    let initial = bulk(&chart, 0.5);
    let gas = fixture::gas(&catalog, 350.0, 500.0);
    let saved = gas.clone();
    let wrong = GasSurfaceBinding { solid_enthalpy_offset_j_kg: 1000.0, ..binding.clone() };
    assert!(matches!(advance_equilibrium_surface(inventory, &fixture::phase(), &wrong,
        initial, &gas, &catalog, None, 0.01, fixture::controls(), || false),
        Err(GasSurfaceError::CaloricMismatch { .. })));
    let mut calls = 0;
    advance_equilibrium_surface(inventory, &fixture::phase(), &binding, initial, &gas,
        &catalog, None, 0.01, fixture::controls(), || { calls += 1; false }).unwrap();
    for stop in [1, calls / 2, calls] {
        let mut count = 0;
        assert!(advance_equilibrium_surface(inventory, &fixture::phase(), &binding,
            initial, &gas, &catalog, None, 0.01, fixture::controls(),
            || { count += 1; count == stop }).is_err());
        assert_eq!(gas, saved);
        assert_eq!(initial, bulk(&chart, 0.5));
    }
}
