//! Prepared UV regressions using independent mass action and energy equations.
use super::super::volume::energy_path::EquilibriumEnergyVolumePathV1;
use super::*;

fn controls() -> EquilibriumInternalEnergySolveV1 {
    let mut c = uv_controls();
    c.volume.equilibrium.relative_balance_tolerance = 1e-13;
    c.volume.equilibrium.dimensionless_potential_tolerance = 1e-12;
    c.volume.relative_volume_tolerance = 1e-11;
    c.energy_tolerance_j = 1e-7;
    c
}
fn fraction(t: f64) -> f64 {
    1.0 / (1.0 + (-20000.0 / (R * t)).exp())
}
fn iso_u(t: f64) -> f64 {
    3.0 * R * t - 20000.0 * fraction(t)
}
fn bisect(target: f64, energy: impl Fn(f64) -> f64) -> f64 {
    let (mut low, mut high) = (200.0, 2500.0);
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
fn g1_prepared_uv_matches_isomer_energy_and_variable_composition_cv() {
    let chemistry = isomers(20000.0);
    let feed = [(id("A"), Amount::new(2.0))];
    let mut path = EquilibriumEnergyVolumePathV1::new(
        &chemistry,
        Temperature::new(900.0),
        Volume::new(0.02),
        &feed,
        controls(),
        || false,
    )
    .unwrap();
    for (energy, volume) in [(4000.0, 0.02), (8000.0, 0.04), (-6000.0, 0.01), (0.0, 0.02)] {
        let sample = path
            .sample(Volume::new(volume), Energy::new(energy), || false)
            .unwrap();
        let expected = bisect(iso_u(900.0) + energy / 2.0, iso_u);
        close(sample.temperature.value(), expected, 2e-10);
        close(sample.pressure.value(), 2.0 * R * expected / volume, 2e-10);
        assert!(sample.energy_residual_j.abs() <= controls().energy_tolerance_j);
        close(sample.internal_energy_change_j, energy, 1e-7);
        let y = fraction(expected);
        let expected_cv =
            2.0 * (3.0 * R + 20000.0_f64.powi(2) * y * (1.0 - y) / (R * expected * expected));
        close(sample.response.heat_capacity_v_j_k, expected_cv, 2e-10);
        let state = path.publish_current(|| false).unwrap();
        close(
            species(state.equilibrium(), "B").amount.value(),
            2.0 * y,
            2e-11,
        );
        assert_eq!(state.volume(), Volume::new(volume));
        assert_eq!(state.target_internal_energy(), None);
        assert!(
            state
                .equilibrium()
                .species()
                .iter()
                .all(|s| s.source.is_some())
        );
        close(path.feed_mass_kg(), 0.02, 1e-15);
    }
}

#[test]
fn g1_free_expansion_changes_dissociation_temperature_at_conserved_energy() {
    let chemistry = dissociation();
    let alpha = |t: f64, v: f64| {
        let kp = (t / REFERENCE_T).powf(1.5);
        2.0 / (1.0 + (1.0 + 16.0 * R * t / (v * 100000.0 * kp)).sqrt())
    };
    let u = |t: f64, v: f64| R * t * (2.5 + 0.5 * alpha(t, v));
    let feed = [(id("A2"), Amount::new(1.0))];
    let mut path = EquilibriumEnergyVolumePathV1::new(
        &chemistry,
        Temperature::new(1000.0),
        Volume::new(0.02),
        &feed,
        controls(),
        || false,
    )
    .unwrap();
    let sample = path
        .sample(Volume::new(0.2), Energy::new(0.0), || false)
        .unwrap();
    let expected = bisect(u(1000.0, 0.02), |t| u(t, 0.2));
    close(sample.temperature.value(), expected, 2e-10);
    assert!(sample.temperature.value() < 1000.0);
    assert!(sample.entropy_change_j_k > 0.0);
    assert!(sample.internal_energy_change_j.abs() <= controls().energy_tolerance_j);
    let state = path.publish_current(|| false).unwrap();
    close(
        species(state.equilibrium(), "A").amount.value(),
        2.0 * alpha(expected, 0.2),
        2e-10,
    );
    close(
        sample.pressure.value(),
        (1.0 + alpha(expected, 0.2)) * R * expected / 0.2,
        2e-10,
    );
}

#[test]
fn g3_relative_energy_retains_small_heat_and_extensive_scaling() {
    let q = 0.001;
    for offset in [0.0, 1e11] {
        let chemistry = owner(
            vec![model("A", 0.01, 4.0, offset, offset / (R * REFERENCE_T))],
            &["H"],
            vec![vec![1]],
        );
        for scale in [0.25, 1.0, 4.0] {
            let mut c = controls();
            c.energy_tolerance_j = 1e-8 * scale;
            c.maximum_temperature_iterations = 1;
            let mut path = EquilibriumEnergyVolumePathV1::new(
                &chemistry,
                Temperature::new(900.0),
                Volume::new(0.02 * scale),
                &[(id("A"), Amount::new(2.0 * scale))],
                c,
                || false,
            )
            .unwrap();
            let step = path
                .sample(Volume::new(0.03 * scale), Energy::new(q * scale), || false)
                .unwrap();
            assert_eq!(step.temperature_iterations, 1);
            close(step.temperature.value(), 900.0 + q / (6.0 * R), 1e-12);
            close(step.internal_energy_change_j, q * scale, 1e-8 * scale);
            close(
                step.pressure.value(),
                2.0 * R * step.temperature.value() / 0.03,
                1e-12,
            );
        }
    }
}

#[test]
fn g0_g4_failed_queries_invalidate_publication_and_retry_preserves_reference() {
    let chemistry = isomers(20000.0);
    let feed = [(id("A"), Amount::new(2.0))];
    let mut path = EquilibriumEnergyVolumePathV1::new(
        &chemistry,
        Temperature::new(900.0),
        Volume::new(0.02),
        &feed,
        controls(),
        || false,
    )
    .unwrap();
    let reference = path.reference_state().clone();
    let same = path
        .sample(Volume::new(0.02), Energy::new(0.0), || false)
        .unwrap();
    assert_eq!(
        (same.pressure_iterations, same.composition_iterations, same.temperature_iterations),
        (0, 0, 0)
    );
    assert!(
        path.sample(Volume::new(0.03), Energy::new(4000.0), || true)
            .is_err()
    );
    assert!(path.publish_current(|| false).is_err());
    for (v, q) in [(0.0, 0.0), (0.02, f64::NAN), (0.02, 1e12)] {
        assert!(
            path.sample(Volume::new(v), Energy::new(q), || false)
                .is_err()
        );
        assert!(path.publish_current(|| false).is_err());
        assert_eq!(path.reference_state(), &reference);
    }
    let accepted = path
        .sample(Volume::new(0.03), Energy::new(4000.0), || false)
        .unwrap();
    close(
        accepted.temperature.value(),
        bisect(iso_u(900.0) + 2000.0, iso_u),
        2e-10,
    );
    assert!(
        accepted.pressure_iterations
            <= (controls().maximum_temperature_iterations + 2)
                * controls().volume.maximum_pressure_iterations
    );
    let mut capped = controls();
    capped.maximum_temperature_iterations = 1;
    let mut path = EquilibriumEnergyVolumePathV1::new(
        &chemistry,
        Temperature::new(900.0),
        Volume::new(0.02),
        &feed,
        capped,
        || false,
    )
    .unwrap();
    assert!(matches!(
        path.sample(Volume::new(0.03), Energy::new(10000.0), || false),
        Err(EquilibriumVolumeErrorV1::InternalEnergyNoConvergence { iterations: 1, .. })
    ));
    assert!(path.publish_current(|| false).is_err());
}
