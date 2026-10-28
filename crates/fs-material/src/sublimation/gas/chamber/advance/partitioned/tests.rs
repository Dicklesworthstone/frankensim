//! Real phase/gas owners with manufactured sharp-phase calorics and kinetic oracles.
use super::*;
use crate::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use crate::sublimation::gas::{
    GasSurfaceBinding, SublimationError, SublimationLaw,
    chamber::{GasChamberSolve, advance::advance_chamber_adaptive, mixed::advance_partitioned_chamber},
    equilibrium::{EquilibriumSublimation, advance_equilibrium_surface},
    tests as fixture,
};
use fs_blake3::ContentHash;
use fs_thermal::sublimation::SolidSurfaceState;
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K};

const FUSION: f64 = 100000.0;
fn chart(catalog: &ThermochemicalGasModel, binding: &GasSurfaceBinding) -> EquilibriumEnthalpyPhaseCurve {
    let (h, cp) = fixture::chart(catalog, binding);
    let knot = |h, t, f| EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h, temperature_k: t,
        liquid_mass_fraction: f, bulk_density_kg_m3: 917.0,
    };
    EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([0x67; 32]), vec![
        knot(h - 200.0 * cp, 200.0, 0.0), knot(h, 400.0, 0.0),
        knot(h + FUSION, 400.0, 1.0), knot(h + FUSION + 500.0 * cp, 900.0, 1.0),
    ]).unwrap()
}
fn bulk(chart: &EquilibriumEnthalpyPhaseCurve, fraction: f64) -> SolidSurfaceState {
    SolidSurfaceState {
        mass_kg: 0.01, area_m2: 1e-4,
        specific_enthalpy_j_kg: chart.knots()[1].specific_enthalpy_j_kg + fraction * FUSION,
    }
}
fn temperature(chart: &EquilibriumEnthalpyPhaseCurve, h: f64) -> Result<f64, SublimationError> {
    chart.temperature_at_specific_enthalpy(h)
        .map_err(|_| SublimationError::Invalid("original sharp fusion chart refused"))
}
fn surface<'a>(id: u64, law: &'a dyn SublimationLaw, binding: &'a GasSurfaceBinding,
    solid: SolidSurfaceState, caloric: &'a dyn Fn(f64) -> Result<f64, SublimationError>) -> GasChamberSurface<'a> {
    GasChamberSurface { id, law, binding, solid, temperature_from_enthalpy: caloric }
}
fn controls(mass_tolerance: f64) -> GasChamberAdvanceSolve {
    GasChamberAdvanceSolve {
        step: GasChamberSolve { surface: fixture::controls(), maximum_surfaces: 8 },
        absolute_mass_tolerance_kg: mass_tolerance,
        absolute_temperature_tolerance_k: 1e-4,
        relative_tolerance: 0.0, minimum_dt_s: 1e-10,
        maximum_dt_s: 0.2, maximum_trials: 100000,
    }
}
fn reflection(catalog: &ThermochemicalGasModel, alpha: f64) -> ThermalAccommodation {
    ThermalAccommodation::new(catalog, &[
        (SpeciesId::new("A").unwrap(), alpha), (fixture::species(), alpha),
    ], [200.0, 900.0]).unwrap()
}

#[test]
fn phase_mass_error_alone_resolves_reflected_heat_on_a_constant_temperature_plateau() {
    let catalog = fixture::catalog(0.0);
    let mut binding = fixture::binding(0.0);
    binding.accommodation_coefficient = 0.0;
    let curve = chart(&catalog, &binding);
    let caloric = |h| temperature(&curve, h);
    let phase = fixture::phase();
    let initial = bulk(&curve, 0.5);
    let rows = [surface(1, &phase, &binding, initial, &caloric)];
    let sources = BTreeMap::from([(1, SharpMeltingCaloric::new(&curve).unwrap())]);
    let accommodation = BTreeMap::from([(1, reflection(&catalog, 0.5))]);
    let gas = fixture::gas(&catalog, 350.0, 0.0);
    let duration = 0.2;
    // Ts is fixed by the latent plateau, not an external thermostat. Finite gas
    // warming consumes actual fusion energy in the retained condensed store.
    // For a diatomic constant-cp carrier: dTg/dt=K*sqrt(Tg)*(Ts-Tg)/Cg.
    let r = UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.02897;
    let cg = gas.mass_kg() * 2.5 * r;
    let k = 0.5 * initial.area_m2 * gas.mass_kg() / gas.volume_m3()
        * (r / (2.0 * core::f64::consts::PI)).sqrt() * 3.0 * r;
    let x = (350.0_f64 / 400.0).sqrt();
    let z = 0.5 * fs_math::det::ln((1.0 + x) / (1.0 - x))
        + k * 400.0_f64.sqrt() * duration / (2.0 * cg);
    let v = fs_math::det::tanh(z);
    let exact_gas_temperature = 400.0 * v * v;
    let exact_liquid = 0.005 - cg * (exact_gas_temperature - 350.0) / FUSION;
    let mut results = Vec::new();
    for tolerance in [1e-5, 1e-10] {
        let mut solve = controls(tolerance);
        // Temperature cannot force refinement; total and gas masses are constant.
        solve.absolute_temperature_tolerance_k = 1e6;
        let result = advance_partitioned_chamber_adaptive(&rows, &sources, Some(&accommodation),
            &gas, &catalog, duration, solve, || false).unwrap();
        let report = result.history().audit();
        let end = &result.phases()[&1];
        assert_eq!(end.final_inventory.caloric.temperature_k(), 400.0);
        assert_eq!(result.history().endpoint().surfaces()[0].solid.mass_kg, initial.mass_kg);
        assert_eq!(result.gas().mass_kg(), gas.mass_kg());
        assert_eq!(report.surfaces[0].transferred_mass_kg, 0.0);
        assert!(report.maximum_accepted_error_ratio <= 1.0);
        fixture::near(report.reflected_heat_into_gas_j, -FUSION * end.melted_mass_kg, 1e-7);
        fixture::near(report.energy_residual_j, 0.0, 1e-7);
        fixture::near(report.gas_stream_energy_residual_j, 0.0, 1e-7);
        fixture::near(end.solid_mass_residual_kg, 0.0, 1e-14);
        results.push(((end.final_inventory.liquid_mass_kg - exact_liquid).abs(), report.accepted_intervals));
    }
    assert!(results[1].0 < results[0].0 / 100.0, "{results:?}");
    assert!(results[1].0 < 5e-8, "{results:?}");
    assert!(results[1].1 > results[0].1);
    assert_eq!(gas.state().state().temperature, 350.0);
    assert_eq!(rows[0].solid, initial);
}

#[test]
fn entirely_solid_history_keeps_the_original_adaptive_endpoint_and_audit_bits() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let curve = chart(&catalog, &binding);
    let caloric = |h| temperature(&curve, h);
    let phase = fixture::phase();
    let initial = fixture::solid(300.0, fixture::chart(&catalog, &binding));
    let rows = [surface(7, &phase, &binding, initial, &caloric)];
    let sources = BTreeMap::from([(7, SharpMeltingCaloric::new(&curve).unwrap())]);
    let gas = fixture::gas(&catalog, 300.0, 0.0);
    let mut solve = controls(1e-10);
    solve.maximum_dt_s = 0.005;
    let original = advance_chamber_adaptive(&rows, &gas, &catalog, 0.01, solve, || false).unwrap();
    let mixed = advance_partitioned_chamber_adaptive(&rows, &sources, None,
        &gas, &catalog, 0.01, solve, || false).unwrap();
    assert_eq!(mixed.history(), &original);
    assert_eq!(mixed.phases()[&7].melted_mass_kg, 0.0);
    assert_eq!(mixed.phases()[&7].final_inventory.liquid_mass_kg, 0.0);
}

#[test]
fn two_condensed_bodies_share_one_gas_with_canonical_conservative_phase_histories() {
    let catalog = fixture::catalog(-1e7);
    let binding = fixture::binding(-1e7);
    let curve = chart(&catalog, &binding);
    let caloric = |h| temperature(&curve, h);
    let phase = fixture::phase();
    let rows = [surface(20, &phase, &binding, bulk(&curve, 0.5), &caloric),
        surface(10, &phase, &binding, fixture::solid(300.0, fixture::chart(&catalog, &binding)), &caloric)];
    let source = SharpMeltingCaloric::new(&curve).unwrap();
    let sources = BTreeMap::from([(10, source), (20, source)]);
    let accommodation = BTreeMap::from([(20, reflection(&catalog, 0.4)), (10, reflection(&catalog, 0.2))]);
    let gas = fixture::gas(&catalog, 350.0, 800.0);
    let solve = controls(1e-9);
    let forward = advance_partitioned_chamber_adaptive(&rows, &sources, Some(&accommodation),
        &gas, &catalog, 0.04, solve, || false).unwrap();
    let reverse = advance_partitioned_chamber_adaptive(&[rows[1], rows[0]], &sources, Some(&accommodation),
        &gas, &catalog, 0.04, solve, || false).unwrap();
    assert_eq!(forward, reverse);
    let report = forward.history().audit();
    assert!(report.surfaces[1].transferred_mass_kg > 0.0);
    assert!(report.surfaces[0].transferred_mass_kg < 0.0);
    assert!(forward.phases()[&20].melted_mass_kg < 0.0);
    for account in forward.phases().values() {
        fixture::near(account.solid_mass_residual_kg, 0.0, 1e-14);
        assert_eq!(account.final_inventory.caloric.phase_curve_identity(), curve.identity());
    }
    fixture::near(report.energy_residual_j, 0.0, 1e-7);
    fixture::near(report.mass_residual_kg, 0.0, 1e-14);
    fixture::near(forward.gas().species_mass_kg(&SpeciesId::new("A").unwrap()).unwrap(),
        gas.species_mass_kg(&SpeciesId::new("A").unwrap()).unwrap(), 1e-14);
    fixture::near(report.energy_into_gas_j, report.gas_energy_change_j, 1e-7);
    assert_eq!(forward.gas().point().axes().get("local-lot"), Some(&9.0));
}

#[test]
fn tiny_solid_inventory_drives_shared_retry_without_spending_liquid_mass() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let curve = chart(&catalog, &binding);
    let caloric = |h| temperature(&curve, h);
    let phase = fixture::phase();
    let initial = bulk(&curve, 0.99999);
    let rows = [surface(1, &phase, &binding, initial, &caloric)];
    let sources = BTreeMap::from([(1, SharpMeltingCaloric::new(&curve).unwrap())]);
    let gas = fixture::gas(&catalog, 400.0, 0.0);
    let result = advance_partitioned_chamber_adaptive(&rows, &sources, None,
        &gas, &catalog, 0.01, controls(1e-9), || false).unwrap();
    assert!(result.history().audit().rejected_trials > 0);
    assert_eq!(result.history().audit().elapsed_s, 0.01);
    assert!(result.phases()[&1].initial.solid_mass_kg < 1.1e-7);
    assert!(result.phases()[&1].final_inventory.solid_mass_kg > 0.0);
    assert!(result.phases()[&1].melted_mass_kg < 0.0);
    fixture::near(result.history().audit().mass_residual_kg, 0.0, 1e-14);
    fixture::near(result.history().audit().energy_residual_j, 0.0, 1e-7);
    assert_eq!(initial, rows[0].solid);
}

#[test]
fn scalar_facade_and_partitioned_chamber_use_identical_source_and_storage_owners() {
    let catalog = fixture::catalog(0.0);
    let binding = fixture::binding(0.0);
    let curve = chart(&catalog, &binding);
    let caloric = |h| temperature(&curve, h);
    let phase = fixture::phase();
    let initial = bulk(&curve, 0.5);
    let gas = fixture::gas(&catalog, 350.0, 500.0);
    let inventory = EquilibriumSublimation::new(&curve).unwrap();
    let thermal = reflection(&catalog, 0.5);
    let scalar = advance_equilibrium_surface(inventory, &phase, &binding, initial,
        &gas, &catalog, Some(&thermal), 0.001, fixture::controls(), || false).unwrap();
    let rows = [surface(9, &phase, &binding, initial, &caloric)];
    let sources = BTreeMap::from([(9, inventory.source())]);
    let accommodation = BTreeMap::from([(9, thermal)]);
    let joint = advance_partitioned_chamber(&rows, &sources, Some(&accommodation),
        &gas, &catalog, 0.001, controls(1e-9).step, || false).unwrap();
    assert_eq!(scalar.exchange().gas_exchange(), joint.exchange().gas_exchange());
    assert_eq!(scalar.exchange().solid(), joint.exchange().surfaces()[0].solid);
    fixture::near(scalar.phase().freezing_mass_kg, -joint.phases()[&9].melted_mass_kg, 1e-18);
    assert_eq!(scalar.reflection(), joint.reflection().get(&9));
}

#[test]
fn zero_time_admission_and_failed_or_cancelled_histories_preserve_every_input() {
    let catalog = fixture::catalog(0.0);
    let mut binding = fixture::binding(0.0);
    binding.accommodation_coefficient = 0.0;
    let curve = chart(&catalog, &binding);
    let caloric = |h| temperature(&curve, h);
    let phase = fixture::phase();
    let initial = bulk(&curve, 0.5);
    let rows = [surface(1, &phase, &binding, initial, &caloric)];
    let source = SharpMeltingCaloric::new(&curve).unwrap();
    let sources = BTreeMap::from([(1, source)]);
    let accommodation = BTreeMap::from([(1, reflection(&catalog, 0.5))]);
    let gas = fixture::gas(&catalog, 350.0, 0.0);
    let saved = gas.clone();
    let solve = controls(1e-10);
    let zero = advance_partitioned_chamber_adaptive(&rows, &sources, Some(&accommodation),
        &gas, &catalog, 0.0, solve, || false).unwrap();
    assert_eq!(zero.gas(), &gas);
    assert_eq!(zero.history().audit().accepted_intervals, 0);
    assert_eq!(zero.history().audit().reflected_heat_into_gas_j, 0.0);
    assert_eq!(zero.phases()[&1].initial, zero.phases()[&1].final_inventory);
    for missing in [BTreeMap::new(), BTreeMap::from([(2, source)])] {
        assert!(advance_partitioned_chamber_adaptive(&rows, &missing, None,
            &gas, &catalog, 0.0, solve, || false).is_err());
    }
    let limited = GasChamberAdvanceSolve { maximum_trials: 1, ..solve };
    assert!(matches!(advance_partitioned_chamber_adaptive(&rows, &sources, Some(&accommodation),
        &gas, &catalog, 0.2, limited, || false), Err(GasChamberAdvanceError::TrialBudget { .. })));
    let liquid = [surface(1, &phase, &binding, bulk(&curve, 1.0), &caloric)];
    assert!(advance_partitioned_chamber_adaptive(&liquid, &sources, None,
        &gas, &catalog, 0.0, solve, || false).is_err());
    let mut calls = 0;
    advance_partitioned_chamber_adaptive(&rows, &sources, Some(&accommodation),
        &gas, &catalog, 0.01, solve, || { calls += 1; false }).unwrap();
    for stop in [1, calls / 2, calls] {
        let mut count = 0;
        assert!(advance_partitioned_chamber_adaptive(&rows, &sources, Some(&accommodation),
            &gas, &catalog, 0.01, solve, || { count += 1; count == stop }).is_err());
        assert_eq!(gas, saved);
        assert_eq!(rows[0].solid, initial);
    }
}
