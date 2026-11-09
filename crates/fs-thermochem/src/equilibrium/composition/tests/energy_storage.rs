//! Original-source reacting internal-energy storage at fixed conserved density.
use super::super::volume::energy_path::isochoric::EquilibriumIsochoricCaloricV1;
use super::*;

fn controls() -> EquilibriumInternalEnergySolveV1 {
    let mut c = uv_controls();
    c.maximum_temperature_k = 2000.0;
    c.volume.equilibrium.relative_balance_tolerance = 1e-13;
    c.volume.equilibrium.dimensionless_potential_tolerance = 1e-12;
    c.volume.relative_volume_tolerance = 1e-12;
    c.energy_tolerance_j = 1e-8;
    c
}
fn feed() -> [(SpeciesId, Amount); 1] {
    [(id("A"), Amount::new(2.0))]
}
fn fraction(t: f64) -> f64 {
    1.0 / (1.0 + (-20000.0 / (R * t)).exp())
}
fn energy(t: f64) -> f64 {
    2.0 * (3.0 * R * t - 20000.0 * fraction(t))
}
fn entropy(t: f64) -> f64 {
    let y = fraction(t);
    // At fixed V and fixed molecular count, p is proportional to T. Include
    // that original pressure entropy once, leaving 3R ln(T), not 4R ln(T).
    2.0 * R * (3.0 * t.ln() - y * y.ln() - (1.0 - y) * (1.0 - y).ln())
}

#[test]
fn g1_isochoric_storage_uses_reacting_cv_pressure_and_source_entropy() {
    let chemistry = isomers(20000.0);
    let mut law = EquilibriumIsochoricCaloricV1::new(
        &chemistry,
        Volume::new(0.1),
        &feed(),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(law.feed_mass_kg(), 0.02);
    assert_eq!(law.density_kg_m3(), 0.02 / 0.1);
    let cold = controls().minimum_temperature_k;
    for t in [400.0, 900.0, 1300.0] {
        let u = (energy(t) - energy(cold)) / 0.02;
        let direct = law
            .sample_temperature(Temperature::new(t), || false)
            .unwrap();
        close(
            direct.internal_energy_change_j,
            energy(t) - energy(cold),
            1e-10,
        );
        close(direct.pressure.value(), 2.0 * R * t / 0.1, 2e-12);
        close(direct.entropy_change_j_k, entropy(t) - entropy(cold), 2e-11);
        let inverse = law.sample_specific_internal_energy(u, || false).unwrap();
        close(inverse.temperature.value(), t, 2e-11);
        let y = fraction(t);
        let cv = 2.0 * (3.0 * R + 20000.0_f64.powi(2) * y * (1.0 - y) / (R * t * t));
        close(inverse.response.heat_capacity_v_j_k, cv, 2e-11);
        let state = law.publish_current(|| false).unwrap();
        assert!(state.equilibrium().is_bound_to(&chemistry));
        assert_eq!(state.target_internal_energy(), None);
        assert_eq!(state.volume(), Volume::new(0.1));
    }
    let a = (energy(900.0) - energy(cold)) / 0.02;
    let b = (energy(1300.0) - energy(cold)) / 0.02;
    let change = law
        .change_between_specific_internal_energies(a, b, || false)
        .unwrap();
    close(
        change.internal_energy_change_j,
        energy(1300.0) - energy(900.0),
        2e-10,
    );
    close(
        change.entropy_change_j_k,
        entropy(1300.0) - entropy(900.0),
        2e-11,
    );
}

#[test]
fn g1_dissociation_storage_closes_fixed_volume_and_extensive_scaling() {
    let chemistry = dissociation();
    let volume = 0.1;
    let alpha = |t: f64| {
        let k = (t / REFERENCE_T).powf(1.5);
        let q = R * t / (volume * 1e5);
        2.0 / (1.0 + (1.0 + 16.0 * q / k).sqrt())
    };
    let u = |t: f64| R * t * (2.5 + 0.5 * alpha(t));
    for scale in [0.001, 1.0, 1000.0] {
        let mut c = controls();
        c.energy_tolerance_j *= scale;
        let mut law = EquilibriumIsochoricCaloricV1::new(
            &chemistry,
            Volume::new(volume * scale),
            &[(id("A2"), Amount::new(scale))],
            c,
            || false,
        )
        .unwrap();
        for t in [400.0, 1000.0, 1800.0] {
            let state = law
                .sample_temperature(Temperature::new(t), || false)
                .unwrap();
            let a = alpha(t);
            let da = 0.5 * a * (1.0 - a) / (t * (2.0 - a));
            let cv = R * (2.5 + 0.5 * a + 0.5 * t * da);
            close(
                state.internal_energy_change_j / scale,
                u(t) - u(c.minimum_temperature_k),
                2e-10,
            );
            close(state.pressure.value(), (1.0 + a) * R * t / volume, 2e-11);
            close(state.response.heat_capacity_v_j_k / scale, cv, 2e-10);
            close(state.response.density_kg_m3, 0.002 / volume, 2e-12);
            // Reacting Cp-Cv is not the frozen molecular-count R shortcut.
            assert!((state.response.heat_capacity_p_j_k / scale - cv - (1.0 + a) * R).abs() > 1e-3);
            let inverse = law
                .sample_specific_internal_energy(
                    (u(t) - u(c.minimum_temperature_k)) / 0.002,
                    || false,
                )
                .unwrap();
            close(inverse.temperature.value(), t, 2e-10);
        }
    }
}

#[test]
fn g3_isochoric_queries_restore_the_same_source_and_pressure_anchor() {
    let chemistry = isomers(20000.0);
    let mut law = EquilibriumIsochoricCaloricV1::new(
        &chemistry,
        Volume::new(0.1),
        &feed(),
        controls(),
        || false,
    )
    .unwrap();
    let u = (energy(900.0) - energy(controls().minimum_temperature_k)) / 0.02;
    let sample = law.sample_specific_internal_energy(u, || false).unwrap();
    let state = law.publish_current(|| false).unwrap();
    for other in [
        0.0,
        u * 1.3,
        u * 0.25,
        law.specific_internal_energy_domain_j_kg()[1],
    ] {
        law.sample_specific_internal_energy(other, || false)
            .unwrap();
        assert_eq!(
            law.sample_specific_internal_energy(u, || false).unwrap(),
            sample
        );
        assert_eq!(law.publish_current(|| false).unwrap(), state);
    }
    assert!(
        law.sample_temperature(Temperature::new(500.0), || true)
            .is_err()
    );
    assert!(law.publish_current(|| false).is_err());
    assert_eq!(
        law.sample_specific_internal_energy(u, || false).unwrap(),
        sample
    );
    law.sample_specific_internal_energy(0.0, || false).unwrap();
    let cold = law.publish_current(|| false).unwrap();
    assert_eq!(
        cold.internal_energy(),
        law.reference_state().internal_energy()
    );
    assert_eq!(
        cold.thermodynamic().entropy_j_k(),
        law.reference_state().thermodynamic().entropy_j_k()
    );
    assert_eq!(
        cold.thermodynamic().response(),
        law.reference_state().thermodynamic().response()
    );
}

#[test]
fn g3_local_isochoric_heat_retains_small_changes_beside_formation_offsets() {
    for offset in [-1e11, 0.0, 1e11] {
        let chemistry = owner(
            vec![model("A", 0.01, 4.0, offset, offset / (R * REFERENCE_T))],
            &["H"],
            vec![vec![1]],
        );
        let mut law = EquilibriumIsochoricCaloricV1::new(
            &chemistry,
            Volume::new(0.1),
            &feed(),
            controls(),
            || false,
        )
        .unwrap();
        let a = 6.0 * R * (900.0 - controls().minimum_temperature_k) / 0.02;
        let b = a + 0.001 / 0.02;
        let change = law
            .change_between_specific_internal_energies(a, b, || false)
            .unwrap();
        close(change.internal_energy_change_j, 0.001, 2e-8);
        let t1 = 900.0 + 0.001 / (6.0 * R);
        close(
            change.entropy_change_j_k,
            6.0 * R * (t1 / 900.0).ln(),
            2e-11,
        );
        let reverse = law
            .change_between_specific_internal_energies(b, a, || false)
            .unwrap();
        close(
            reverse.internal_energy_change_j,
            -change.internal_energy_change_j,
            1e-13,
        );
        close(
            reverse.entropy_change_j_k,
            -change.entropy_change_j_k,
            1e-13,
        );
    }
}

#[test]
fn g4_isochoric_caloric_bounds_work_and_late_cancellation_refuse_atomically() {
    let chemistry = isomers(20000.0);
    let make = |c| {
        EquilibriumIsochoricCaloricV1::new(&chemistry, Volume::new(0.1), &feed(), c, || false)
            .unwrap()
    };
    let mut law = make(controls());
    let upper = law.specific_internal_energy_domain_j_kg()[1];
    for u in [-1.0, upper.next_up(), f64::NAN, f64::INFINITY] {
        assert!(law.sample_specific_internal_energy(u, || false).is_err());
        assert!(law.publish_current(|| false).is_err());
    }
    let (a, b) = (0.2 * upper, 0.6 * upper);
    let mut count = 0;
    law.change_between_specific_internal_energies(a, b, || {
        count += 1;
        false
    })
    .unwrap();
    let mut cancelled = make(controls());
    let mut polls = 0;
    assert!(matches!(
        cancelled.change_between_specific_internal_energies(a, b, || {
            polls += 1;
            polls == count
        }),
        Err(EquilibriumVolumeErrorV1::Equilibrium(
            CompositionEquilibriumErrorV1::Cancelled
        ))
    ));
    assert_eq!(polls, count);
    assert!(cancelled.publish_current(|| false).is_err());
    assert!(
        cancelled
            .sample_specific_internal_energy(a, || false)
            .is_ok()
    );
    let mut c = controls();
    c.maximum_temperature_iterations = 1;
    let mut limited = make(c);
    assert!(matches!(
        limited.sample_specific_internal_energy(b, || false),
        Err(EquilibriumVolumeErrorV1::InternalEnergyNoConvergence { iterations: 1, .. })
    ));
    assert!(limited.publish_current(|| false).is_err());
}
