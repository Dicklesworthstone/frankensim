//! Original-potential work, coupling, tangent and complete-path admission.
use super::{close, source};
use fs_contact::normal_patch::adhesive::mixed_mode::{
    SourceJkrMixedMode,
    combined::{JkrCombinedResponse, path::JkrCombinedPathError},
};
use fs_exec::CancelGate;

fn fixture(source: &SourceJkrMixedMode) -> [f64; 4] {
    let normal = source.normal_source();
    let a = normal
        .evaluate_force_controlled(0.0)
        .unwrap()
        .contact_radius_m;
    let e = normal.normal_state().elastic().reduced_modulus_pa();
    let w = normal.normal_state().adhesion_energy_j_per_m2();
    let shear =
        (0.12 * core::f64::consts::PI * w * a / (2.0 * source.reduced_shear_modulus_pa())).sqrt();
    let twist = (0.18 * core::f64::consts::PI * w
        / (4.0 * source.torsional_reduced_shear_modulus_pa() * a))
        .sqrt();
    [
        a * a / normal.effective_radius_m()
            - (0.70 * 2.0 * core::f64::consts::PI * w * a / e).sqrt(),
        0.6 * shear,
        0.8 * shear,
        twist,
    ]
}
fn storage(point: JkrCombinedResponse) -> f64 {
    point.normal_elastic_energy_j
        + point.tangential_elastic_energy_j
        + point.torsional_elastic_energy_j
        - point.surface_energy_j
}
fn work_close(actual: f64, expected: f64, scale: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * scale,
        "{actual:e} != {expected:e}; scale={scale:e}"
    );
}
fn dot(a: [f64; 4], b: [f64; 4]) -> f64 {
    (0..4).map(|i| a[i] * b[i]).sum()
}

#[test]
fn one_area_discrete_work_agrees_with_original_forces_and_reverses() {
    let source = source();
    let gate = CancelGate::new();
    let q = fixture(&source);
    for (from, to) in [
        (
            [0.85 * q[0], 0.5 * q[1], 0.4 * q[2], 0.3 * q[3]],
            [1.10 * q[0], 0.9 * q[1], 0.8 * q[2], 0.9 * q[3]],
        ),
        (
            [q[0], -0.7 * q[1], 0.8 * q[2], -0.6 * q[3]],
            [1.03 * q[0], 0.8 * q[1], -0.6 * q[2], 0.7 * q[3]],
        ),
    ] {
        let path = source.combined_discrete_path(from, to, &gate).unwrap();
        let reverse = source.combined_discrete_path(to, from, &gate).unwrap();
        let scale = storage(path.before) + storage(path.contact);
        work_close(
            path.reversible_work_j,
            path.contact.reversible_energy_j - path.before.reversible_energy_j,
            scale,
            5e-14,
        );
        work_close(path.work_residual_j, 0.0, scale, 5e-14);
        work_close(
            reverse.reversible_work_j,
            -path.reversible_work_j,
            scale,
            2e-15,
        );
        assert_eq!(path.mean_effort, reverse.mean_effort);
        assert!(path.evaluations <= 5);
        assert_eq!(path.contact.model_identity(), source.identity());
        assert_eq!(
            path.contact.normal_state_identity(),
            source.normal_source().normal_state().identity()
        );

        // Independent composite Simpson work of POINTWISE original efforts.
        let increment = core::array::from_fn(|i| to[i] - from[i]);
        let mut integral = 0.0;
        for i in 0..=256 {
            let t = i as f64 / 256.0;
            let point = source
                .evaluate_combined_displacement_controlled(core::array::from_fn(|j| {
                    from[j] + t * increment[j]
                }))
                .unwrap();
            assert!(point.contact_radius_m >= path.minimum_contact_radius_m);
            assert!(point.contact_radius_m <= path.maximum_contact_radius_m);
            let weight = if i == 0 || i == 256 {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            integral += weight * dot(point.force, increment);
        }
        work_close(
            path.reversible_work_j,
            integral / (3.0 * 256.0),
            scale,
            3e-10,
        );

        let mut previous = from;
        let mut subdivided = 0.0;
        for t in [0.07, 0.24, 0.63, 1.0] {
            let next = core::array::from_fn(|i| from[i] + t * increment[i]);
            subdivided += source
                .combined_discrete_path(previous, next, &gate)
                .unwrap()
                .reversible_work_j;
            previous = next;
        }
        work_close(subdivided, path.reversible_work_j, scale, 5e-14);
    }
}

#[test]
fn all_sixteen_endpoint_derivatives_and_zero_increment_limit_are_consistent() {
    let source = source();
    let gate = CancelGate::new();
    let q = fixture(&source);
    let from = [0.85 * q[0], 0.4 * q[1], 0.5 * q[2], 0.3 * q[3]];
    let to = [1.1 * q[0], 0.8 * q[1], 0.7 * q[2], 0.8 * q[3]];
    let path = source.combined_discrete_path(from, to, &gate).unwrap();
    for j in 0..4 {
        let step = 1e-5 * q[j];
        let mut plus = to;
        let mut minus = to;
        plus[j] += step;
        minus[j] -= step;
        let p = source.combined_discrete_path(from, plus, &gate).unwrap();
        let m = source.combined_discrete_path(from, minus, &gate).unwrap();
        for i in 0..4 {
            close(
                (p.mean_effort[i] - m.mean_effort[i]) / (2.0 * step),
                path.end_tangent[i][j],
                4e-7,
            );
        }
    }
    let same = source.combined_discrete_path(to, to, &gate).unwrap();
    let original = source
        .evaluate_combined_displacement_controlled(to)
        .unwrap();
    assert_eq!(same.mean_effort, original.force);
    assert_eq!(
        same.end_tangent,
        original
            .coupled_stiffness
            .unwrap()
            .map(|row| row.map(|x| 0.5 * x))
    );
    assert_eq!(same.reversible_work_j, 0.0);
    assert_eq!(same.evaluations, 1);
}

#[test]
fn sealed_cache_retains_single_ulp_work_and_refused_trials_do_not_mutate_history() {
    let source = source();
    let gate = CancelGate::new();
    let q = fixture(&source);
    let prepared = source.prepare_combined_discrete_path(q, &gate).unwrap();
    let original = source.evaluate_combined_displacement_controlled(q).unwrap();
    assert!(core::ptr::eq(prepared.source(), &source));
    for i in 0..4 {
        for upper in [false, true] {
            let mut next = q;
            next[i] = if upper {
                next[i].next_up()
            } else {
                next[i].next_down()
            };
            let trial = prepared.trial(next, &gate).unwrap();
            let expected = original.force[i] * (next[i] - q[i]);
            close(trial.path().reversible_work_j, expected, 5e-12);
            assert_ne!(trial.path().reversible_work_j, 0.0);
            assert!(trial.path().evaluations <= 4);
            assert_eq!(trial.endpoint().coordinates(), next);
            assert_eq!(prepared.coordinates(), q);
        }
    }
    let same = prepared.trial(q, &gate).unwrap();
    assert_eq!(same.path().evaluations, 0);
    assert_eq!(same.endpoint(), prepared);
    let mut invalid = q;
    invalid[2] = f64::NAN;
    assert!(prepared.trial(invalid, &gate).is_err());
    assert_eq!(prepared.coordinates(), q);
    gate.request();
    assert_eq!(
        prepared.trial(q, &gate),
        Err(JkrCombinedPathError::Cancelled)
    );
    assert_eq!(prepared.coordinates(), q);
}

#[test]
fn simultaneous_shear_and_twist_reversal_encloses_hidden_area_expansion() {
    let source = source();
    let gate = CancelGate::new();
    let q = fixture(&source);
    let end = [q[0], -q[1], -q[2], -q[3]];
    let path = source.combined_discrete_path(q, end, &gate).unwrap();
    let zero = source
        .evaluate_combined_displacement_controlled([q[0], 0.0, 0.0, 0.0])
        .unwrap();
    assert!(zero.contact_area_m2 > path.before.contact_area_m2);
    assert!(path.maximum_contact_radius_m >= zero.contact_radius_m);
    assert_eq!(&path.mean_effort[1..], &[0.0; 3]);
    assert_eq!(path.reversible_work_j, 0.0);
    assert_eq!(
        path.before.reversible_energy_j,
        path.contact.reversible_energy_j
    );

    // A change of tangent-frame axes preserves work and rotates both shears.
    let rotate = |p: [f64; 4]| [p[0], -p[2], p[1], p[3]];
    let rotated = source
        .combined_discrete_path(rotate(q), rotate(end), &gate)
        .unwrap();
    close(rotated.mean_effort[0], path.mean_effort[0], 3e-14);
    assert_eq!(&rotated.mean_effort[1..], &[0.0; 3]);
    assert_eq!(rotated.reversible_work_j, path.reversible_work_j);
    close(
        rotated.contact.contact_radius_m,
        path.contact.contact_radius_m,
        3e-14,
    );
}

#[test]
fn original_pure_limits_keep_the_same_potential_work() {
    let source = source();
    let gate = CancelGate::new();
    let q = fixture(&source);
    let from = [0.8 * q[0], 0.0, 0.0, 0.3 * q[3]];
    let to = [1.1 * q[0], 0.0, 0.0, 0.8 * q[3]];
    let combined = source.combined_discrete_path(from, to, &gate).unwrap();
    let torsion = source
        .torsional_discrete_path([from[0], from[3]], [to[0], to[3]], &gate)
        .unwrap();
    let scale = storage(combined.before) + storage(combined.contact);
    work_close(
        combined.reversible_work_j,
        torsion.reversible_work_j,
        scale,
        3e-14,
    );
    close(combined.mean_effort[0], torsion.mean_effort[0], 3e-13);
    close(combined.mean_effort[3], torsion.mean_effort[1], 3e-13);
    assert_eq!(combined.mean_effort[1], 0.0);
    assert_eq!(combined.mean_effort[2], 0.0);

    let from = [0.8 * q[0], 0.3 * q[1], -0.2 * q[2], 0.0];
    let to = [1.1 * q[0], 0.8 * q[1], 0.7 * q[2], 0.0];
    let combined = source.combined_discrete_path(from, to, &gate).unwrap();
    let shear = source
        .discrete_path([from[0], from[1], from[2]], [to[0], to[1], to[2]], &gate)
        .unwrap();
    let scale = storage(combined.before) + storage(combined.contact);
    work_close(
        combined.reversible_work_j,
        shear.reversible_work_j,
        scale,
        3e-14,
    );
    for i in 0..3 {
        close(combined.mean_effort[i], shear.mean_force_n[i], 3e-13);
    }
    assert_eq!(combined.mean_effort[3], 0.0);
}

#[test]
fn endpoint_attachment_does_not_admit_a_combined_path_crossing_the_opening_edge() {
    let source = source();
    let gate = CancelGate::new();
    let q = fixture(&source);
    let low = 0.25 * q[0];
    let phi_low = 0.99
        * source
            .torsional_displacement_limit(low)
            .unwrap()
            .unwrap()
            .twist_angle_rad;
    let phi_high = 0.99
        * source
            .torsional_displacement_limit(q[0])
            .unwrap()
            .unwrap()
            .twist_angle_rad;
    let from = [low, 0.02 * q[1], -0.01 * q[2], phi_low];
    let to = [q[0], -0.03 * q[1], 0.015 * q[2], phi_high];
    source
        .evaluate_combined_displacement_controlled(from)
        .unwrap();
    source
        .evaluate_combined_displacement_controlled(to)
        .unwrap();
    assert!(
        source
            .evaluate_combined_displacement_controlled(core::array::from_fn(|i| {
                from[i].midpoint(to[i])
            }))
            .is_err()
    );
    assert_eq!(
        source.combined_discrete_path(from, to, &gate),
        Err(JkrCombinedPathError::UnresolvedEnvelope)
    );
    let fold = [
        source.normal_source().minimum_attached_indentation_m(),
        0.0,
        0.0,
        0.0,
    ];
    assert_eq!(
        source.combined_discrete_path(fold, fold, &gate),
        Err(JkrCombinedPathError::UnresolvedSecant)
    );
}
