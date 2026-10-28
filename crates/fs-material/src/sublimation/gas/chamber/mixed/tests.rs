//! Manufactured caloric/phase fixtures, not empirical coexistence validation.
use super::*;
use crate::phase::EnthalpyPhaseKnot;
use crate::sublimation::gas::{GasSurfaceBinding, tests as fixture};
use fs_blake3::ContentHash;
use fs_thermal::sublimation::SublimationError;
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

pub(super) const FUSION: f64 = 80000.0;
pub(super) fn curve(catalog: &ThermochemicalGasModel, binding: &GasSurfaceBinding) -> EquilibriumEnthalpyPhaseCurve {
    let (hs, cp) = fixture::chart(catalog, binding);
    EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([57; 32]), vec![
        EnthalpyPhaseKnot { specific_enthalpy_j_kg: hs - 200.0 * cp, temperature_k: 200.0,
            liquid_mass_fraction: 0.0, bulk_density_kg_m3: 917.0 },
        EnthalpyPhaseKnot { specific_enthalpy_j_kg: hs, temperature_k: 400.0,
            liquid_mass_fraction: 0.0, bulk_density_kg_m3: 917.0 },
        EnthalpyPhaseKnot { specific_enthalpy_j_kg: hs + FUSION, temperature_k: 400.0,
            liquid_mass_fraction: 1.0, bulk_density_kg_m3: 1000.0 },
        EnthalpyPhaseKnot { specific_enthalpy_j_kg: hs + FUSION + 500.0 * cp, temperature_k: 900.0,
            liquid_mass_fraction: 1.0, bulk_density_kg_m3: 1000.0 },
    ]).unwrap()
}
pub(super) fn invert(curve: &EquilibriumEnthalpyPhaseCurve, h: f64) -> Result<f64, SublimationError> {
    curve.temperature_at_specific_enthalpy(h).map_err(|_| SublimationError::Invalid("fixture caloric domain"))
}
pub(super) fn controls() -> GasChamberSolve {
    GasChamberSolve { surface: fixture::controls(), maximum_surfaces: 8 }
}
pub(super) fn state(curve: &EquilibriumEnthalpyPhaseCurve, fraction: f64) -> SolidSurfaceState {
    SolidSurfaceState { mass_kg: 0.01, specific_enthalpy_j_kg: curve.knots()[1].specific_enthalpy_j_kg + fraction * FUSION,
        area_m2: 1e-4 }
}
pub(super) fn coefficients(catalog: &ThermochemicalGasModel, alpha: f64) -> ThermalAccommodation {
    ThermalAccommodation::new(catalog, &[(SpeciesId::new("A").unwrap(), alpha), (fixture::species(), alpha)],
        [200.0, 900.0]).unwrap()
}

#[test]
fn sublimation_freezes_liquid_without_cooling_or_double_counting_fusion() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = curve(&catalog, &binding);
    let caloric = SharpMeltingCaloric::new(&chart).unwrap();
    let phase = fixture::phase();
    let temperature = |h| invert(&chart, h);
    let bulk = state(&chart, 0.5);
    let surface = GasChamberSurface { id: 1, law: &phase, binding: &binding, solid: bulk,
        temperature_from_enthalpy: &temperature };
    let gas = fixture::gas(&catalog, 400.0, 0.0);
    let sources = BTreeMap::from([(1, caloric)]);
    let trial = advance_partitioned_chamber(&[surface], &sources, None, &gas, &catalog,
        0.1, controls(), || false).unwrap();
    let row = &trial.exchange().surfaces()[0];
    let report = &trial.phases()[&1];
    let r = R / 0.018015;
    let dm = 0.1 * 1000.0 / (2.0 * core::f64::consts::PI * r * 400.0).sqrt() * 1e-4 * 0.1;
    let effective_latent = 200000.0 - 0.5 * r * 400.0;
    let liquid = 0.005 - dm * effective_latent / FUSION;
    fixture::near(row.transferred_mass_kg, dm, 1e-18);
    fixture::near(report.final_inventory.liquid_mass_kg, liquid, 1e-16);
    fixture::near(report.final_inventory.solid_mass_kg, 0.01 - dm - liquid, 1e-16);
    assert_eq!(row.surface_temperature_k, [400.0, 400.0]);
    assert!(report.melted_mass_kg < 0.0);
    assert!(report.final_inventory.solid_mass_kg > report.initial.solid_mass_kg);
    assert_eq!(report.final_inventory.caloric.phase_curve_identity(), chart.identity());
    fixture::near(report.solid_mass_residual_kg, 0.0, 1e-14);
    fixture::near(trial.exchange().step().energy_residual_j, 0.0, 1e-7);
    // The legacy pure-solid API must NOT silently reinterpret bulk enthalpy.
    assert!(super::super::advance_chamber_with_gas(&[surface], &gas, &catalog,
        0.1, controls(), || false).is_err());
}

#[test]
fn deposition_melts_the_original_solid_partition_and_carriers_remain_unchanged() {
    let catalog = fixture::catalog(-1e7);
    let binding = fixture::binding(-1e7);
    let chart = curve(&catalog, &binding);
    let phase = fixture::phase();
    let temperature = |h| invert(&chart, h);
    let surface = GasChamberSurface { id: 1, law: &phase, binding: &binding,
        solid: state(&chart, 0.5), temperature_from_enthalpy: &temperature };
    let gas = fixture::gas(&catalog, 400.0, 5000.0);
    let sources = BTreeMap::from([(1, SharpMeltingCaloric::new(&chart).unwrap())]);
    let candidate = advance_partitioned_chamber(&[surface], &sources, None, &gas, &catalog,
        0.01, controls(), || false).unwrap();
    let row = &candidate.exchange().surfaces()[0];
    let inventory = &candidate.phases()[&1];
    assert!(row.transferred_mass_kg < 0.0);
    assert!(row.outgoing_mass_kg > 0.0);
    assert!(inventory.melted_mass_kg > 0.0);
    assert_eq!(row.surface_temperature_k[1], 400.0);
    let latent = 200000.0 - 0.5 * (R / 0.018015) * 400.0;
    fixture::near(inventory.melted_mass_kg, -row.transferred_mass_kg * latent / FUSION, 1e-14);
    let carrier = SpeciesId::new("A").unwrap();
    fixture::near(candidate.exchange().gas().species_mass_kg(&carrier).unwrap(), gas.species_mass_kg(&carrier).unwrap(), 1e-14);
    fixture::near(candidate.exchange().step().energy_residual_j, 0.0, 1e-7);
}

#[test]
fn liquid_inventory_cannot_fund_gross_solid_emission() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = curve(&catalog, &binding);
    let phase = fixture::phase();
    let temperature = |h| invert(&chart, h);
    let surface = GasChamberSurface { id: 1, law: &phase, binding: &binding,
        solid: state(&chart, 0.999), temperature_from_enthalpy: &temperature };
    let gas = fixture::gas(&catalog, 400.0, 0.0);
    let sources = BTreeMap::from([(1, SharpMeltingCaloric::new(&chart).unwrap())]);
    let result = advance_partitioned_chamber(&[surface], &sources, None, &gas, &catalog,
        1.0, controls(), || false);
    let Err(GasSurfaceError::InventoryLimit { maximum_duration_s }) = result else {
        panic!("the small old SOLID inventory must limit the interval");
    };
    assert!(maximum_duration_s > 0.0 && maximum_duration_s < 1.0);
    let accepted = advance_partitioned_chamber(&[surface], &sources, None, &gas, &catalog,
        0.9 * maximum_duration_s, controls(), || false).unwrap();
    assert!(accepted.exchange().surfaces()[0].outgoing_mass_kg
        <= controls().surface.maximum_gross_inventory_fraction * accepted.phases()[&1].initial.solid_mass_kg);
}

#[test]
fn reflected_heat_changes_fusion_at_exactly_zero_phase_mass_flux() {
    let catalog = fixture::catalog(0.0);
    let mut binding = fixture::binding(0.0);
    binding.accommodation_coefficient = 0.0;
    let chart = curve(&catalog, &binding);
    let phase = fixture::phase();
    let temperature = |h| invert(&chart, h);
    let surface = GasChamberSurface { id: 1, law: &phase, binding: &binding,
        solid: state(&chart, 0.5), temperature_from_enthalpy: &temperature };
    let gas = fixture::gas(&catalog, 500.0, 0.0);
    let sources = BTreeMap::from([(1, SharpMeltingCaloric::new(&chart).unwrap())]);
    let reflection = BTreeMap::from([(1, coefficients(&catalog, 0.5))]);
    let candidate = advance_partitioned_chamber(&[surface], &sources, Some(&reflection), &gas,
        &catalog, 0.001, controls(), || false).unwrap();
    let heat = candidate.reflection()[&1].heat_into_gas_j;
    assert!(heat < 0.0);
    assert_eq!(candidate.exchange().step().transferred_mass_kg, 0.0);
    assert_eq!(candidate.exchange().surfaces()[0].surface_temperature_k, [400.0, 400.0]);
    fixture::near(candidate.phases()[&1].melted_mass_kg, -heat / FUSION, 1e-16);
    fixture::near(candidate.exchange().step().energy_residual_j, 0.0, 1e-7);
}

#[test]
fn pure_solid_limit_preserves_original_chamber_bits_and_input_order() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = curve(&catalog, &binding);
    let phase = fixture::phase();
    let temperature = |h| invert(&chart, h);
    let rows = [350.0, 300.0].map(|t| GasChamberSurface { id: t as u64,
        law: &phase, binding: &binding, solid: fixture::solid(t, fixture::chart(&catalog, &binding)),
        temperature_from_enthalpy: &temperature });
    let sources = BTreeMap::from([(350, SharpMeltingCaloric::new(&chart).unwrap()),
        (300, SharpMeltingCaloric::new(&chart).unwrap())]);
    let gas = fixture::gas(&catalog, 350.0, 500.0);
    let original = super::super::advance_chamber_with_gas(&rows, &gas, &catalog, 0.001, controls(), || false).unwrap();
    let mixed = advance_partitioned_chamber(&rows, &sources, None, &gas, &catalog, 0.001, controls(), || false).unwrap();
    assert_eq!(mixed.exchange(), &original);
    let reversed = advance_partitioned_chamber(&[rows[1], rows[0]], &sources, None,
        &gas, &catalog, 0.001, controls(), || false).unwrap();
    assert_eq!(mixed, reversed);
}

#[test]
fn source_shape_no_solid_and_late_cancellation_are_transactional() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let chart = curve(&catalog, &binding);
    let mut mushy = chart.knots().to_vec();
    mushy[2].temperature_k = 410.0;
    let mushy = EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([59; 32]), mushy).unwrap();
    assert!(SharpMeltingCaloric::new(&mushy).is_err());
    let source = SharpMeltingCaloric::new(&chart).unwrap();
    assert!(source.resolve(chart.knots()[2].specific_enthalpy_j_kg + 1.0).is_err());
    let phase = fixture::phase();
    let temperature = |h| invert(&chart, h);
    let surface = GasChamberSurface { id: 1, law: &phase, binding: &binding,
        solid: state(&chart, 0.5), temperature_from_enthalpy: &temperature };
    let gas = fixture::gas(&catalog, 400.0, 0.0);
    let saved = gas.clone();
    let sources = BTreeMap::from([(1, source)]);
    let all_liquid = GasChamberSurface { solid: state(&chart, 1.0), ..surface };
    assert!(advance_partitioned_chamber(&[all_liquid], &sources, None, &gas, &catalog,
        0.001, controls(), || false).is_err());
    assert!(advance_partitioned_chamber(&[surface], &BTreeMap::new(), None, &gas, &catalog,
        0.001, controls(), || false).is_err());
    let mut calls = 0;
    advance_partitioned_chamber(&[surface], &sources, None, &gas, &catalog,
        0.001, controls(), || { calls += 1; false }).unwrap();
    for stop in [1, calls / 2, calls] {
        let mut count = 0;
        assert!(advance_partitioned_chamber(&[surface], &sources, None, &gas, &catalog,
            0.001, controls(), || { count += 1; count == stop }).is_err());
        assert_eq!(gas, saved);
        assert_eq!(surface.solid, state(&chart, 0.5));
    }
}
