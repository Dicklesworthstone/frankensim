//! G1/G3: exact scalar work derivatives retain source phases and signed paths.
#[path = "normal_work_derivative/fixtures.rs"]
mod fixtures;
use fixtures::{close, curve, dmt, point};
use fs_contact::normal_patch::{NormalLawSource, NormalLawSourceError};

#[test]
fn hertz_derivative_retains_inactive_time_and_hunt_crossley_reversible_storage() {
    let d = 1e-6;
    for chi in [None, Some(0.4)] {
        let point = point(chi);
        let k = point.stiffness_n_per_m_3_2();
        let source = NormalLawSource::from(point.clone());
        let force = |x: f64| k * x.max(0.0).powf(1.5);
        let energy = |x: f64| 0.4 * k * x.max(0.0).powf(2.5);
        for (a, b) in [
            (-2.0, 1.0),
            (1.0, -2.0),
            (0.2, 2.3),
            (2.3, 0.2),
            (0.0, 1.7),
            (1.7, 0.0),
            (-2.0, -1.0),
        ] {
            let (a, b) = (a * d, b * d);
            let mean = (energy(b) - energy(a)) / (b - a);
            let expected = (force(b) - mean) / (b - a);
            close(
                source.mean_force_end_derivative(a, b).unwrap().unwrap(),
                expected,
                k * d.sqrt(),
                3e-14,
            );
        }
        for x in [0.2 * d, d, 3.0 * d] {
            let expected = 0.75 * k * x.sqrt();
            close(
                source.mean_force_end_derivative(x, x).unwrap().unwrap(),
                expected,
                expected,
                3e-15,
            );
            close(
                source
                    .mean_force_end_derivative(x, x.next_up())
                    .unwrap()
                    .unwrap(),
                expected,
                expected,
                3e-15,
            );
            let (u, f, t) = source.potential(x).unwrap();
            close(u, energy(x), energy(x), 3e-15);
            close(f, force(x), force(x), 3e-15);
            close(t, 1.5 * k * x.sqrt(), k * x.sqrt(), 3e-15);
        }
        assert_eq!(
            source.mean_force_end_derivative(0.0, 0.0).unwrap(),
            Some(0.0)
        );
        assert!(
            source
                .mean_force_end_derivative(d, point.envelope().maximum_indentation_m.next_up())
                .is_err()
        );
        assert!(source.mean_force_end_derivative(f64::NAN, d).is_err());
    }
}

#[test]
fn dmt_derivative_integrates_negative_tail_stiffness_and_both_crossings() {
    let law = dmt();
    let z = law.interaction_range_m();
    let source = NormalLawSource::from(law.clone());
    let scale = law.force_controlled_pull_off_n().abs() / z;
    for (a, b) in [
        (-2.0, -1.3),
        (-1.8, -0.4),
        (-0.4, -1.8),
        (-0.9, -0.1),
        (-0.1, -0.9),
        (-0.5, 0.4),
        (0.4, -0.5),
        (-1.8, 0.6),
        (0.6, -1.8),
        (0.1, 0.8),
        (0.8, 0.1),
        (-1.0, 0.0),
        (0.0, -1.0),
    ] {
        let (a, b) = (a * z, b * z);
        let ra = law.evaluate(a).unwrap();
        let rb = law.evaluate(b).unwrap();
        // Independent endpoint potential identity on finite, resolved paths.
        let mean = (rb.reversible_energy_j - ra.reversible_energy_j) / (b - a);
        let expected = (rb.normal_force_n - mean) / (b - a);
        close(
            source.mean_force_end_derivative(a, b).unwrap().unwrap(),
            expected,
            scale,
            3e-13,
        );
    }
    for at in [-0.6 * z, 0.4 * z] {
        let half = 0.5 * law.evaluate(at).unwrap().tangent_n_per_m.unwrap();
        close(
            source
                .mean_force_end_derivative(at, at.next_up())
                .unwrap()
                .unwrap(),
            half,
            scale,
            3e-14,
        );
    }
    // At a coincident corner, use the explicitly documented right Newton
    // representative; the original smooth potential view still refuses it.
    for at in [-z, 0.0] {
        assert_eq!(
            source.mean_force_end_derivative(at, at).unwrap(),
            Some(0.5 * law.evaluate(at).unwrap().one_sided_tangents_n_per_m[1])
        );
        assert!(matches!(
            source.potential(at),
            Err(NormalLawSourceError::Dmt(_))
        ));
    }
    close(
        source
            .mean_force_end_derivative(-z, (-z).next_up())
            .unwrap()
            .unwrap(),
        -0.5 * scale,
        scale,
        3e-15,
    );
    assert!(
        source
            .mean_force_end_derivative(0.0, law.envelope().maximum_indentation_m.next_up())
            .is_err()
    );
}

#[test]
fn finite_gap_derivative_integrates_original_segments_and_unilateral_fraction() {
    let d = 1e-6;
    let table = curve(d, 0.01);
    let source = NormalLawSource::from(table.clone());
    let scale = 0.01 / d;
    for (a, b) in [
        (-0.5, 0.7),
        (0.7, -0.5),
        (0.1, 0.8),
        (0.8, 0.1),
        (0.2, 0.5),
        (0.5, 0.2),
        (0.0, 0.9),
        (0.9, 0.0),
    ] {
        let (a, b) = (a * d, b * d);
        let ra = table.evaluate(a.max(0.0)).unwrap();
        let rb = table.evaluate(b.max(0.0)).unwrap();
        let mean = (rb.reversible_energy_j - ra.reversible_energy_j) / (b - a);
        let end_force = if b < 0.0 { 0.0 } else { rb.normal_force_n };
        let expected = (end_force - mean) / (b - a);
        close(
            source.mean_force_end_derivative(a, b).unwrap().unwrap(),
            expected,
            scale,
            3e-14,
        );
    }
    for at in [0.0, 0.1 * d, 0.2 * d, 0.5 * d, d] {
        assert_eq!(
            source.mean_force_end_derivative(at, at).unwrap(),
            Some(0.5 * table.evaluate(at).unwrap().normal_tangent_n_per_m)
        );
    }
    let at = 0.3 * d;
    close(
        source
            .mean_force_end_derivative(at, at.next_up())
            .unwrap()
            .unwrap(),
        0.5 * table.evaluate(at).unwrap().normal_tangent_n_per_m,
        scale,
        3e-15,
    );
    assert_eq!(
        source.mean_force_end_derivative(-d, 0.0).unwrap(),
        Some(0.0)
    );
    assert_eq!(
        source.mean_force_end_derivative(0.0, -d).unwrap(),
        Some(0.0)
    );
    assert!(source.mean_force_end_derivative(0.0, d.next_up()).is_err());
    let mut invalid = table;
    invalid.nodes[2].curve_reversible_energy_j += 1e-3;
    assert!(
        NormalLawSource::from(invalid)
            .mean_force_end_derivative(0.1 * d, 0.9 * d)
            .is_err()
    );
}
