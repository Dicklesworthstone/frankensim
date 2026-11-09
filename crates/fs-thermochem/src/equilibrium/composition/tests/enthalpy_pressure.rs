//! G1/G3 pressure/enthalpy endpoints against independent mass-action equations.
use super::super::enthalpy::EquilibriumEnthalpyPressurePathV1;
use super::*;

fn controls() -> EquilibriumEnthalpySolveV1 {
    let mut c = hp_controls();
    c.equilibrium.relative_balance_tolerance = 1e-13;
    c.equilibrium.dimensionless_potential_tolerance = 1e-12;
    c.enthalpy_tolerance_j = 1e-7;
    c
}

fn bisect(target: f64, energy: impl Fn(f64) -> f64) -> f64 {
    let (mut low, mut high) = (200.0, 2000.0);
    for _ in 0..80 {
        let t = 0.5 * (low + high);
        if energy(t) < target {
            low = t;
        } else {
            high = t;
        }
    }
    0.5 * (low + high)
}

#[test]
fn g1_signed_heating_uses_reacting_cp_and_original_chemical_energy() {
    let chemistry = isomers(20000.0);
    let y = |t: f64| 1.0 / (1.0 + (-20000.0 / (R * t)).exp());
    let h = |t: f64| 4.0 * R * t - 20000.0 * y(t);
    let mut path = EquilibriumEnthalpyPressurePathV1::new(
        &chemistry,
        Temperature::new(900.0),
        Pressure::new(100000.0),
        &[(id("A"), Amount::new(2.0))],
        controls(),
        || false,
    )
    .unwrap();
    for (heat, pressure) in [(4000.0, 100000.0), (8000.0, 200000.0), (-6000.0, 100000.0)] {
        let sample = path
            .sample(Pressure::new(pressure), Energy::new(heat), || false)
            .unwrap();
        let expected = bisect(h(900.0) + heat / 2.0, h);
        close(sample.temperature.value(), expected, 2e-10);
        close(sample.enthalpy_change_j, heat, 1e-8);
        assert!(sample.enthalpy_residual_j.abs() <= controls().enthalpy_tolerance_j);
        let fraction = y(expected);
        let cp = 2.0
            * (4.0 * R
                + 20000.0_f64.powi(2) * fraction * (1.0 - fraction) / (R * expected.powi(2)));
        close(sample.response.heat_capacity_p_j_k, cp, 2e-10);
        let mixing = |x: f64| -R * (x * x.ln() + (1.0 - x) * (1.0 - x).ln());
        let entropy = 2.0
            * (4.0 * R * (expected / 900.0).ln() - R * (pressure / 100000.0).ln()
                + mixing(fraction)
                - mixing(y(900.0)));
        close(sample.entropy_change_j_k, entropy, 2e-9);
        let state = path.publish_current(|| false).unwrap();
        assert!(state.equilibrium().is_bound_to(&chemistry));
        close(
            species(state.equilibrium(), "B").amount.value(),
            2.0 * fraction,
            2e-11,
        );
        assert_eq!(state.equilibrium().target_enthalpy(), None);
        assert_eq!(state.response(), sample.response);
        assert!(
            state
                .equilibrium()
                .species()
                .iter()
                .all(|s| s.source.is_some())
        );
    }
}

#[test]
fn g1_isenthalpic_pressure_drop_cools_dissociating_gas_and_increases_entropy() {
    let chemistry = dissociation();
    let alpha = |t: f64, p: f64| {
        let kp = (t / REFERENCE_T).powf(1.5);
        (kp / (kp + 4.0 * p / 100000.0)).sqrt()
    };
    let h = |t: f64, p: f64| R * t * (3.5 + 1.5 * alpha(t, p));
    let mut path = EquilibriumEnthalpyPressurePathV1::new(
        &chemistry,
        Temperature::new(1000.0),
        Pressure::new(100000.0),
        &[(id("A2"), Amount::new(1.0))],
        controls(),
        || false,
    )
    .unwrap();
    let sample = path
        .sample(Pressure::new(10000.0), Energy::new(0.0), || false)
        .unwrap();
    let expected = bisect(h(1000.0, 100000.0), |t| h(t, 10000.0));
    close(sample.temperature.value(), expected, 2e-10);
    assert!(sample.temperature.value() < 1000.0);
    assert!(sample.entropy_change_j_k > 0.0);
    assert!(sample.enthalpy_change_j.abs() <= controls().enthalpy_tolerance_j);
    let state = path.publish_current(|| false).unwrap();
    close(
        species(state.equilibrium(), "A").amount.value(),
        2.0 * alpha(expected, 10000.0),
        2e-10,
    );
}

#[test]
fn g1_pressure_response_matches_equilibrium_joule_thomson_derivative() {
    let chemistry = dissociation();
    let mut path = EquilibriumEnthalpyPressurePathV1::new(
        &chemistry,
        Temperature::new(1000.0),
        Pressure::new(100000.0),
        &[(id("A2"), Amount::new(1.0))],
        controls(),
        || false,
    )
    .unwrap();
    let state = path.reference_state();
    let response = state.response();
    let volume = state.equilibrium().total_amount().value() * R * 1000.0 / 100000.0;
    // dH=Cp dT+V(1-T alpha)dP at fixed conserved elemental inventory.
    let expected =
        volume * (1000.0 * response.thermal_expansion_per_k - 1.0) / response.heat_capacity_p_j_k;
    let low = path
        .sample(Pressure::new(99900.0), Energy::new(0.0), || false)
        .unwrap();
    let high = path
        .sample(Pressure::new(100100.0), Energy::new(0.0), || false)
        .unwrap();
    let observed = (high.temperature.value() - low.temperature.value()) / 200.0;
    assert!(expected > 0.0);
    assert!((observed / expected - 1.0).abs() < 2e-6);
}

#[test]
fn g3_small_heat_survives_large_absolute_reference_and_extensive_scaling() {
    let heat = 0.001;
    for offset in [-1e11, 0.0, 1e11] {
        let chemistry = owner(
            vec![model("A", 0.01, 4.0, offset, offset / (R * REFERENCE_T))],
            &["H"],
            vec![vec![1]],
        );
        for scale in [0.25, 1.0, 4.0] {
            let mut c = controls();
            c.enthalpy_tolerance_j = 1e-8 * scale;
            c.maximum_temperature_iterations = 1;
            let mut path = EquilibriumEnthalpyPressurePathV1::new(
                &chemistry,
                Temperature::new(900.0),
                Pressure::new(100000.0),
                &[(id("A"), Amount::new(2.0 * scale))],
                c,
                || false,
            )
            .unwrap();
            let step = path
                .sample(Pressure::new(30000.0), Energy::new(heat * scale), || false)
                .unwrap();
            assert_eq!(step.temperature_iterations, 1);
            close(step.temperature.value(), 900.0 + heat / (8.0 * R), 1e-12);
            close(step.enthalpy_change_j, heat * scale, 1e-8 * scale);
            close(path.feed_mass_kg(), 0.02 * scale, 1e-15);
            let same = path
                .sample(Pressure::new(30000.0), Energy::new(heat * scale), || false)
                .unwrap();
            assert_eq!(same.composition_iterations, 0);
            assert_eq!(same.temperature_iterations, 0);
            assert_eq!(same.enthalpy_change_j, step.enthalpy_change_j);
        }
    }
}

#[test]
fn g0_g4_failed_queries_invalidate_publication_and_preserve_the_reference() {
    let chemistry = isomers(20000.0);
    let mut path = EquilibriumEnthalpyPressurePathV1::new(
        &chemistry,
        Temperature::new(900.0),
        Pressure::new(100000.0),
        &[(id("A"), Amount::new(2.0))],
        controls(),
        || false,
    )
    .unwrap();
    let reference = path.reference_state().clone();
    for (pressure, heat, cancel) in [
        (0.0, 0.0, false),
        (100000.0, 1e20, false),
        (100000.0, 4000.0, true),
    ] {
        assert!(
            path.sample(Pressure::new(pressure), Energy::new(heat), || cancel)
                .is_err()
        );
        assert!(path.publish_current(|| false).is_err());
        assert_eq!(path.reference_state(), &reference);
        let retry = path
            .sample(Pressure::new(100000.0), Energy::new(4000.0), || false)
            .unwrap();
        assert!(retry.enthalpy_residual_j.abs() <= controls().enthalpy_tolerance_j);
        assert!(path.publish_current(|| false).is_ok());
    }
}

#[test]
fn g4_final_poll_cached_queries_and_work_limit_cannot_publish_partial_states() {
    let chemistry = isomers(20000.0);
    let feed = [(id("A"), Amount::new(2.0))];
    let make = |c| {
        EquilibriumEnthalpyPressurePathV1::new(
            &chemistry,
            Temperature::new(900.0),
            Pressure::new(100000.0),
            &feed,
            c,
            || false,
        )
        .unwrap()
    };
    let mut complete = make(controls());
    let mut count = 0;
    complete
        .sample(Pressure::new(80000.0), Energy::new(4000.0), || {
            count += 1;
            false
        })
        .unwrap();
    let mut path = make(controls());
    let original = path.reference_state().clone();
    let mut polls = 0;
    assert!(matches!(
        path.sample(Pressure::new(80000.0), Energy::new(4000.0), || {
            polls += 1;
            polls == count
        }),
        Err(EquilibriumExpansionErrorV1::Equilibrium(
            CompositionEquilibriumErrorV1::Cancelled
        ))
    ));
    assert_eq!(polls, count);
    assert!(path.publish_current(|| false).is_err());
    assert_eq!(path.reference_state(), &original);
    path.sample(Pressure::new(80000.0), Energy::new(4000.0), || false)
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        path.sample(Pressure::new(80000.0), Energy::new(4000.0), || {
            polls += 1;
            polls == 2
        }),
        Err(EquilibriumExpansionErrorV1::Equilibrium(
            CompositionEquilibriumErrorV1::Cancelled
        ))
    ));
    assert!(path.publish_current(|| false).is_err());
    let mut capped = controls();
    capped.maximum_temperature_iterations = 1;
    let mut path = make(capped);
    assert!(matches!(
        path.sample(Pressure::new(100000.0), Energy::new(4000.0), || false),
        Err(EquilibriumExpansionErrorV1::Equilibrium(
            CompositionEquilibriumErrorV1::EnthalpyNoConvergence { iterations: 1, .. }
        ))
    ));
    assert!(path.publish_current(|| false).is_err());
}

#[test]
fn g1_original_isentrope_and_enthalpy_path_share_reacting_source_integrals() {
    use super::super::isentrope::EquilibriumIsentropeV1;
    let chemistry = isomers(20000.0);
    let feed = [(id("A"), Amount::new(2.0))];
    let y = |t: f64| 1.0 / (1.0 + (-20000.0 / (R * t)).exp());
    let h = |t: f64| 4.0 * R * t - 20000.0 * y(t);
    let mixing = |x: f64| -R * (x * x.ln() + (1.0 - x) * (1.0 - x).ln());
    let pressure = 50000.0;
    let entropy_change = |t: f64| {
        2.0 * (4.0 * R * (t / 900.0).ln() - R * (pressure / 100000.0_f64).ln() + mixing(y(t))
            - mixing(y(900.0)))
    };
    let expected = bisect(0.0, entropy_change);
    let dh = 2.0 * (h(expected) - h(900.0));
    let c = controls();
    let mut isentrope = EquilibriumIsentropeV1::new(
        &chemistry,
        Temperature::new(900.0),
        Pressure::new(100000.0),
        &feed,
        EquilibriumEntropySolveV1 {
            equilibrium: c.equilibrium,
            minimum_temperature_k: c.minimum_temperature_k,
            maximum_temperature_k: c.maximum_temperature_k,
            entropy_tolerance_j_k: 1e-10,
            maximum_temperature_iterations: c.maximum_temperature_iterations,
        },
        || false,
    )
    .unwrap();
    let s = isentrope.sample(Pressure::new(pressure), || false).unwrap();
    close(s.temperature.value(), expected, 2e-10);
    close(s.specific_enthalpy_change_j_kg, dh / 0.02, 1e-8);
    let mut path = EquilibriumEnthalpyPressurePathV1::new(
        &chemistry,
        Temperature::new(900.0),
        Pressure::new(100000.0),
        &feed,
        c,
        || false,
    )
    .unwrap();
    let p = path
        .sample(Pressure::new(pressure), Energy::new(dh), || false)
        .unwrap();
    close(p.temperature.value(), expected, 2e-10);
    close(p.entropy_change_j_k, 0.0, 1e-8);
}
