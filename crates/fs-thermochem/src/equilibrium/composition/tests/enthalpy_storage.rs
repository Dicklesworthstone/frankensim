//! G1/G3/G4 repeatable reacting calorics for finite isobaric heat storage.
use super::super::enthalpy::isobaric::EquilibriumIsobaricCaloricV1;
use super::*;

fn controls() -> EquilibriumEnthalpySolveV1 {
    let mut c = hp_controls();
    c.equilibrium.relative_balance_tolerance = 1e-13;
    c.equilibrium.dimensionless_potential_tolerance = 1e-12;
    c.enthalpy_tolerance_j = 1e-8;
    c
}
fn feed() -> [(SpeciesId, Amount); 1] { [(id("A"), Amount::new(2.0))] }
fn fraction(t: f64) -> f64 { 1.0 / (1.0 + (-20000.0 / (R * t)).exp()) }
fn enthalpy(t: f64) -> f64 { 2.0 * (4.0 * R * t - 20000.0 * fraction(t)) }
fn entropy(t: f64) -> f64 {
    let y = fraction(t);
    2.0 * R * (4.0 * t.ln() - y * y.ln() - (1.0 - y) * (1.0 - y).ln())
}

#[test]
fn g1_reacting_storage_forward_inverse_and_entropy_share_original_sources() {
    let chemistry = isomers(20000.0);
    let mut law = EquilibriumIsobaricCaloricV1::new(
        &chemistry, Pressure::new(1e5), &feed(), controls(), || false,
    ).unwrap();
    let cold = controls().minimum_temperature_k;
    for t in [400.0, 900.0, 1300.0] {
        let h = (enthalpy(t) - enthalpy(cold)) / 0.02;
        let direct = law.sample_temperature(Temperature::new(t), || false).unwrap();
        close(direct.enthalpy_change_j, h * 0.02, 1e-8);
        close(direct.entropy_change_j_k, entropy(t) - entropy(cold), 2e-9);
        let inverse = law.sample_specific_enthalpy(h, || false).unwrap();
        close(inverse.temperature.value(), t, 2e-9);
        let y = fraction(t);
        let cp = 2.0 * (4.0 * R + 20000.0_f64.powi(2) * y * (1.0 - y) / (R * t * t));
        close(inverse.response.heat_capacity_p_j_k, cp, 2e-9);
        let state = law.publish_current(|| false).unwrap();
        assert!(state.equilibrium().is_bound_to(&chemistry));
        assert_eq!(state.equilibrium().target_enthalpy(), None);
    }
    let a = (enthalpy(900.0) - enthalpy(cold)) / 0.02;
    let b = (enthalpy(1300.0) - enthalpy(cold)) / 0.02;
    let change = law.change_between_specific_enthalpies(a, b, || false).unwrap();
    close(change.enthalpy_change_j, enthalpy(1300.0) - enthalpy(900.0), 2e-8);
    close(change.entropy_change_j_k, entropy(1300.0) - entropy(900.0), 2e-9);
}

#[test]
fn g1_dissociation_storage_density_and_cp_include_changing_molecular_count() {
    let chemistry = dissociation();
    let pressure = 30000.0;
    let alpha = |t: f64| {
        let k = (t / REFERENCE_T).powf(1.5);
        (k / (k + 4.0 * pressure / 1e5)).sqrt()
    };
    let h = |t: f64| R * t * (3.5 + 1.5 * alpha(t));
    let mut law = EquilibriumIsobaricCaloricV1::new(
        &chemistry, Pressure::new(pressure), &[(id("A2"), Amount::new(1.0))],
        controls(), || false,
    ).unwrap();
    for t in [400.0, 1000.0, 1800.0] {
        let state = law.sample_temperature(Temperature::new(t), || false).unwrap();
        close(state.enthalpy_change_j, h(t) - h(controls().minimum_temperature_k), 1e-8);
        close(state.response.density_kg_m3, 0.002 * pressure / ((1.0 + alpha(t)) * R * t), 2e-12);
        let cp = (h(t + 0.001) - h(t - 0.001)) / 0.002;
        close(state.response.heat_capacity_p_j_k, cp, 2e-7);
        let recovered = law.sample_specific_enthalpy(state.enthalpy_change_j / 0.002, || false).unwrap();
        close(recovered.temperature.value(), t, 2e-9);
    }
}

#[test]
fn g3_original_numerical_anchor_removes_query_order_and_failure_history() {
    let chemistry = isomers(20000.0);
    let mut law = EquilibriumIsobaricCaloricV1::new(
        &chemistry, Pressure::new(1e5), &feed(), controls(), || false,
    ).unwrap();
    let h = (enthalpy(900.0) - enthalpy(controls().minimum_temperature_k)) / 0.02;
    let initial = law.sample_specific_enthalpy(h, || false).unwrap();
    let native = law.publish_current(|| false).unwrap();
    for other in [0.0, h * 1.3, h * 0.25, law.specific_enthalpy_domain_j_kg()[1]] {
        law.sample_specific_enthalpy(other, || false).unwrap();
        assert_eq!(law.sample_specific_enthalpy(h, || false).unwrap(), initial);
        assert_eq!(law.publish_current(|| false).unwrap(), native);
    }
    assert!(law.sample_temperature(Temperature::new(500.0), || true).is_err());
    assert!(law.publish_current(|| false).is_err());
    assert_eq!(law.sample_specific_enthalpy(h, || false).unwrap(), initial);
    // Rewinding must also restore h, g and Cp; a cached cold sample may not
    // publish the previous hot source evaluations under a cold temperature.
    law.sample_specific_enthalpy(0.0, || false).unwrap();
    let cold = law.publish_current(|| false).unwrap();
    assert_eq!(cold.entropy_j_k(), law.reference_state().entropy_j_k());
    assert_eq!(cold.equilibrium().enthalpy(), law.reference_state().equilibrium().enthalpy());
    assert_eq!(cold.response(), law.reference_state().response());
}

#[test]
fn g3_local_caloric_differences_retain_small_heat_far_from_cold_reference() {
    for offset in [-1e11, 0.0, 1e11] {
        let chemistry = owner(
            vec![model("A", 0.01, 4.0, offset, offset / (R * REFERENCE_T))],
            &["H"], vec![vec![1]],
        );
        let mut law = EquilibriumIsobaricCaloricV1::new(
            &chemistry, Pressure::new(1e5), &feed(), controls(), || false,
        ).unwrap();
        let initial_h = 8.0 * R * (900.0 - controls().minimum_temperature_k) / 0.02;
        let final_h = initial_h + 0.001 / 0.02;
        let change = law.change_between_specific_enthalpies(initial_h, final_h, || false).unwrap();
        close(change.enthalpy_change_j, 0.001, 2e-8);
        let t1 = 900.0 + 0.001 / (8.0 * R);
        close(change.entropy_change_j_k, 8.0 * R * (t1 / 900.0).ln(), 2e-11);
        let reverse = law.change_between_specific_enthalpies(final_h, initial_h, || false).unwrap();
        close(reverse.enthalpy_change_j, -change.enthalpy_change_j, 1e-13);
        close(reverse.entropy_change_j_k, -change.entropy_change_j_k, 1e-13);
    }
}

#[test]
fn g4_caloric_domain_work_and_late_interval_cancellation_are_atomic() {
    let chemistry = isomers(20000.0);
    let make = |c| EquilibriumIsobaricCaloricV1::new(
        &chemistry, Pressure::new(1e5), &feed(), c, || false,
    ).unwrap();
    let mut law = make(controls());
    let max = law.specific_enthalpy_domain_j_kg()[1];
    for h in [-1.0, max.next_up(), f64::NAN, f64::INFINITY] {
        assert!(law.sample_specific_enthalpy(h, || false).is_err());
        assert!(law.publish_current(|| false).is_err());
    }
    let (a, b) = (0.2 * max, 0.6 * max);
    let mut count = 0;
    law.change_between_specific_enthalpies(a, b, || { count += 1; false }).unwrap();
    let mut cancelled = make(controls());
    let mut polls = 0;
    assert!(matches!(
        cancelled.change_between_specific_enthalpies(a, b, || { polls += 1; polls == count }),
        Err(EquilibriumExpansionErrorV1::Equilibrium(CompositionEquilibriumErrorV1::Cancelled))
    ));
    assert_eq!(polls, count);
    assert!(cancelled.publish_current(|| false).is_err());
    assert!(cancelled.sample_specific_enthalpy(a, || false).is_ok());
    let mut c = controls();
    c.maximum_temperature_iterations = 1;
    let mut limited = make(c);
    assert!(matches!(
        limited.sample_specific_enthalpy(b, || false),
        Err(EquilibriumExpansionErrorV1::Equilibrium(
            CompositionEquilibriumErrorV1::EnthalpyNoConvergence { iterations: 1, .. }
        ))
    ));
    assert!(limited.publish_current(|| false).is_err());
}
