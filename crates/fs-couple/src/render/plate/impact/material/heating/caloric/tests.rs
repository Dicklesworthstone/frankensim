use super::*;
use fs_blake3::ContentHash;
use fs_material::phase::{EnthalpyPhaseKnot, SolidLiquidPhase};

fn knot(h: f64, t: f64, f: f64) -> EnthalpyPhaseKnot {
    EnthalpyPhaseKnot {
        specific_enthalpy_j_kg: h,
        temperature_k: t,
        liquid_mass_fraction: f,
        bulk_density_kg_m3: 1000.0,
    }
}
fn store<'a>(curve: &'a EquilibriumEnthalpyPhaseCurve, initial_h: f64) -> CaloricReservoir<'a> {
    CaloricReservoir::new(
        "specimen".into(),
        curve,
        curve.state_at_specific_enthalpy(initial_h).unwrap(),
        0.1,
        1e-4,
        0.01,
    )
    .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e} within {tolerance:e}"
    );
}

#[test]
fn original_chart_stages_sensible_and_latent_heat_before_publication() {
    let curve = EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([71; 32]),
        vec![
            knot(0.0, 300.0, 0.0),
            knot(100.0, 310.0, 0.0),
            knot(300.0, 310.0, 1.0),
            knot(500.0, 330.0, 1.0),
        ],
    )
    .unwrap();
    let mut reservoir = store(&curve, 50.0);
    let mut deposited = 0.0;
    for (heat, h, t, f) in [
        (5.0, 100.0, 310.0, 0.0),
        (10.0, 200.0, 310.0, 0.5),
        (10.0, 300.0, 310.0, 1.0),
        (15.0, 450.0, 325.0, 1.0),
    ] {
        let before = reservoir.clone();
        let candidate = reservoir.stage_deposit(heat, 1e-12).unwrap();
        assert_eq!(reservoir, before);
        near(candidate.energy_change_j, heat, 1e-12);
        reservoir.accept_deposit(candidate);
        deposited += heat;
        near(reservoir.state().specific_enthalpy_j_kg(), h, 1e-12);
        near(reservoir.state().temperature_k(), t, 1e-12);
        near(reservoir.state().liquid_mass_fraction(), f, 1e-12);
        near(reservoir.received_heat_j(), deposited, 1e-12);
        near(reservoir.storage_change_j(), deposited, 1e-12);
    }
}

#[test]
fn compensated_small_deposits_retain_the_original_large_enthalpy_datum() {
    let datum = 1e12;
    let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([72; 32]),
        SolidLiquidPhase::Solid,
        vec![knot(datum, 300.0, 0.0), knot(datum + 1.0, 301.0, 0.0)],
    )
    .unwrap();
    let mut incremental = store(&curve, datum);
    let mut one = incremental.clone();
    let count = 8192;
    for _ in 0..count {
        let candidate = incremental.stage_deposit(1e-6, 1e-5).unwrap();
        incremental.accept_deposit(candidate);
    }
    let total = count as f64 * 1e-6;
    let candidate = one.stage_deposit(total, 1e-5).unwrap();
    one.accept_deposit(candidate);
    assert_eq!(incremental.state(), one.state());
    assert_eq!(incremental.received_heat_j(), total);
    near(incremental.storage_change_j(), total, 1e-5);
    assert!(incremental.state().temperature_k() > 300.0);
}

#[test]
fn finite_chart_refusal_never_changes_inventory_and_retry_matches_fresh_owner() {
    let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
        ContentHash([73; 32]),
        SolidLiquidPhase::Solid,
        vec![knot(0.0, 300.0, 0.0), knot(1.0, 301.0, 0.0)],
    )
    .unwrap();
    let mut reservoir = store(&curve, 0.0);
    let original = reservoir.clone();
    for heat in [1.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(reservoir.stage_deposit(heat, 1e-12).is_err());
        assert_eq!(reservoir, original);
    }
    assert!(reservoir.stage_deposit(0.01, 0.0).is_err());
    assert_eq!(reservoir, original);
    let candidate = reservoir.stage_deposit(0.01, 1e-12).unwrap();
    reservoir.accept_deposit(candidate);
    let mut fresh = original.clone();
    fresh.accept_deposit(fresh.stage_deposit(0.01, 1e-12).unwrap());
    assert_eq!(reservoir, fresh);
    reservoir.accept_deposit(reservoir.stage_deposit(0.0, 1e-12).unwrap());
    assert_eq!(reservoir, fresh);
}
