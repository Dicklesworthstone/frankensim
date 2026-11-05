//! Independent work, derivative, covariance and whole-trajectory regressions.
use super::*;
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_contact::normal_patch::adhesive::mixed_mode::path::JkrMixedModePathError;
use fs_exec::CancelGate;

#[path = "../adhesive_dmt_coupled/common.rs"]
mod sourced;

fn q_force(law: &SourceJkrMixedMode, q: [f64; 3]) -> [f64; 3] {
    let point = law
        .evaluate_displacement_controlled(q[0], [q[1], q[2]])
        .unwrap();
    [
        point.normal_force_n,
        point.tangential_force_n[0],
        point.tangential_force_n[1],
    ]
}
fn interpolate(a: [f64; 3], b: [f64; 3], s: f64) -> [f64; 3] {
    core::array::from_fn(|i| a[i] + s * (b[i] - a[i]))
}
fn quadratic(points: [[f64; 3]; 3], s: f64) -> [f64; 3] {
    interpolate(
        interpolate(points[0], points[1], s),
        interpolate(points[1], points[2], s),
        s,
    )
}
fn energy(law: &SourceJkrMixedMode, q: [f64; 3]) -> f64 {
    law.evaluate_displacement_controlled(q[0], [q[1], q[2]])
        .unwrap()
        .reversible_energy_j
}
pub(super) fn sourced_law(strength: f64) -> (SourceJkrMixedMode, SourceAdhesiveShear) {
    let (normal, shear) = sourced::sources(0.05, strength, 0.5 * strength);
    let e = normal.normal_state().elastic().reduced_modulus_pa();
    let spinodal = (PI * 0.05 * 0.01 * 0.01 / (8.0 * e)).cbrt();
    let mut envelope = source(1.0).envelope();
    envelope.maximum_indentation_m = 100.0 * spinodal * spinodal / 0.01;
    let normal =
        SourceJkrSphereContact::new(normal.normal_state().clone(), 0.01, 1e-12, envelope).unwrap();
    (
        SourceJkrMixedMode::new(normal, shear_envelope()).unwrap(),
        shear,
    )
}

#[test]
fn g1_discrete_gradient_closes_original_vector_work_and_reversal() {
    let law = model(1.0);
    let gate = CancelGate::new();
    for (from, to) in [
        ([1e-6, 0.15e-6, -0.2e-6], [5e-6, 0.8e-6, 0.35e-6]),
        ([2e-6, 0.9e-6, 0.0], [2e-6, 0.0, -0.9e-6]),
        ([4e-6, -0.3e-6, 0.6e-6], [1e-6, 0.1e-6, -0.2e-6]),
    ] {
        let path = law.discrete_path(from, to, &gate).unwrap();
        let reverse = law.discrete_path(to, from, &gate).unwrap();
        let scale = path.before.normal_elastic_energy_j
            + path.contact.normal_elastic_energy_j
            + path.before.tangential_elastic_energy_j
            + path.contact.tangential_elastic_energy_j
            - path.before.surface_energy_j
            - path.contact.surface_energy_j;
        close(
            path.reversible_work_j,
            energy(&law, to) - energy(&law, from),
            scale,
            2e-14,
        );
        close(path.work_residual_j, 0.0, scale, 2e-14);
        close(
            reverse.reversible_work_j,
            -path.reversible_work_j,
            scale,
            2e-15,
        );
        assert_eq!(path.mean_force_n, reverse.mean_force_n);
        assert_eq!(path.evaluations, 2);
        assert_eq!(path.before.model_identity(), law.identity());
        assert_eq!(
            path.contact.normal_state_identity(),
            law.normal_source().normal_state().identity()
        );

        // Independent force-line quadrature checks all three conjugate terms;
        // production uses neither this quadrature nor endpoint energy subtraction.
        let mut integrated = 0.0;
        let n = 256;
        for i in 0..=n {
            let point = interpolate(from, to, i as f64 / n as f64);
            let f = q_force(&law, point);
            let weight = if i == 0 || i == n {
                1.0
            } else if i % 2 == 0 {
                2.0
            } else {
                4.0
            };
            integrated += weight * (0..3).map(|j| f[j] * (to[j] - from[j])).sum::<f64>();
        }
        close(
            path.reversible_work_j,
            integrated / (3.0 * n as f64),
            scale,
            2e-10,
        );
    }
}

#[test]
fn prepared_paths_reuse_only_sealed_original_endpoints_and_preserve_failed_trials() {
    let law = model(1.0);
    let gate = CancelGate::new();
    let initial = [2e-6, 0.2e-6, -0.1e-6];
    let mut prepared = law.prepare_discrete_path(initial, &gate).unwrap();
    assert!(core::ptr::eq(prepared.source(), &law));
    for next in [initial, [3e-6, -0.4e-6, 0.3e-6], [1e-6, 0.1e-6, 0.2e-6]] {
        let original = prepared;
        let uncached = law
            .discrete_path(prepared.position_m(), next, &gate)
            .unwrap();
        let trial = prepared.trial(next, &gate).unwrap();
        let mut comparable = *trial.path();
        assert_eq!(
            comparable.evaluations,
            usize::from(next != prepared.position_m())
        );
        comparable.evaluations += 1;
        assert_eq!(comparable, uncached);
        let (endpoint, mut diagnostic) = trial.into_parts();
        diagnostic.contact.contact_radius_m *= 10.0;
        diagnostic.contact.tangential_displacement_m = [99.0; 2];
        assert_eq!(endpoint.contact(), &uncached.contact);
        assert_eq!(endpoint.position_m(), next);
        assert_eq!(prepared, original);
        assert_ne!(diagnostic.contact, *endpoint.contact());
        prepared = endpoint;
    }
    let retained = prepared;
    assert!(prepared.trial([f64::NAN, 0.0, 0.0], &gate).is_err());
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(matches!(
        prepared.trial(initial, &cancelled),
        Err(JkrMixedModePathError::Cancelled)
    ));
    assert!(matches!(
        law.prepare_discrete_path(initial, &cancelled),
        Err(JkrMixedModePathError::Cancelled)
    ));
    assert_eq!(prepared, retained);
    assert_eq!(
        prepared.trial(initial, &gate).unwrap(),
        retained.trial(initial, &gate).unwrap()
    );
}

#[test]
fn g1_analytic_discrete_tangent_matches_independent_endpoint_differences() {
    let law = model(1.0);
    let gate = CancelGate::new();
    let from = [1.5e-6, -0.25e-6, 0.3e-6];
    for to in [[4e-6, 0.7e-6, -0.2e-6], [0.5e-6, -0.35e-6, -0.2e-6], from] {
        let path = law.discrete_path(from, to, &gate).unwrap();
        for j in 0..3 {
            let h = 2e-11;
            let mut plus = to;
            let mut minus = to;
            plus[j] += h;
            minus[j] -= h;
            let fp = law.discrete_path(from, plus, &gate).unwrap().mean_force_n;
            let fm = law.discrete_path(from, minus, &gate).unwrap().mean_force_n;
            for i in 0..3 {
                let fd = (fp[i] - fm[i]) / (2.0 * h);
                close(path.end_tangent_n_per_m[i][j], fd, 3000.0, 3e-8);
            }
        }
    }
}

#[test]
fn g1_zero_and_adjacent_motion_keep_original_scalar_work_and_finite_tangents() {
    let law = model(1.0);
    let gate = CancelGate::new();
    let q = [2e-6, 0.3e-6, -0.2e-6];
    let still = law.discrete_path(q, q, &gate).unwrap();
    assert_eq!(still.mean_force_n, q_force(&law, q));
    assert_eq!(still.reversible_work_j, 0.0);
    assert_eq!(still.work_residual_j, 0.0);
    assert_eq!(still.evaluations, 1);
    assert_eq!(
        still.end_tangent_n_per_m,
        still
            .contact
            .coupled_stiffness_n_per_m
            .unwrap()
            .map(|r| r.map(|x| 0.5 * x))
    );
    for to in [q.map(f64::next_up), [q[0].next_down(), q[1], q[2]]] {
        let tiny = law.discrete_path(q, to, &gate).unwrap();
        for i in 0..3 {
            close(tiny.mean_force_n[i], still.mean_force_n[i], 0.05, 2e-13);
        }
        assert!(
            tiny.end_tangent_n_per_m
                .iter()
                .flatten()
                .all(|x| x.is_finite())
        );
    }
    let minimum = law.normal_source().minimum_attached_indentation_m();
    for (from, to) in [(0.0, 4e-6), (4e-6, 0.0), (0.99 * minimum, 0.95 * minimum)] {
        let path = law
            .discrete_path([from, 0.0, 0.0], [to, 0.0, 0.0], &gate)
            .unwrap();
        let normal = law.normal_source();
        let mean = normal.attached_path_average(from, to).unwrap();
        let tangent = normal.attached_path_end_derivative(from, to).unwrap();
        close(path.mean_force_n[0], mean, mean.abs().max(0.01), 2e-12);
        close(
            path.end_tangent_n_per_m[0][0],
            tangent,
            tangent.abs().max(3000.0),
            2e-10,
        );
        assert_eq!(path.mean_force_n[1..], [0.0; 2]);
    }
    assert_eq!(
        law.discrete_path([minimum, 0.0, 0.0], [minimum, 0.0, 0.0], &gate),
        Err(JkrMixedModePathError::UnresolvedSecant)
    );
}

#[test]
fn g3_discrete_force_tangent_and_work_rotate_in_the_actual_shear_plane() {
    let law = model(1.0);
    let gate = CancelGate::new();
    let from = [1e-6, 0.2e-6, -0.3e-6];
    let to = [3e-6, -0.4e-6, -0.1e-6];
    let angle: f64 = 0.713;
    let (s, c) = angle.sin_cos();
    let rotation = [[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]];
    let rotate =
        |q: [f64; 3]| core::array::from_fn(|i| (0..3).map(|j| rotation[i][j] * q[j]).sum());
    let original = law.discrete_path(from, to, &gate).unwrap();
    let moved = law.discrete_path(rotate(from), rotate(to), &gate).unwrap();
    let expected = rotate(original.mean_force_n);
    for i in 0..3 {
        close(moved.mean_force_n[i], expected[i], 0.05, 2e-14);
        for j in 0..3 {
            let mut value = 0.0;
            for k in 0..3 {
                for l in 0..3 {
                    value += rotation[i][k] * original.end_tangent_n_per_m[k][l] * rotation[j][l];
                }
            }
            close(moved.end_tangent_n_per_m[i][j], value, 3000.0, 2e-13);
        }
    }
    close(
        moved.reversible_work_j,
        original.reversible_work_j,
        1e-8,
        1e-13,
    );
}

#[test]
fn g1_subdivided_quadratic_hulls_enclose_all_actual_radius_shear_and_force() {
    let (law, shear) = sourced_law(1e9);
    let h = law.normal_source().indentation_scale_m();
    let from = [0.0, 0.0, 0.0];
    let to = [
        20.0 * h,
        0.8 * law.critical_tangential_displacement_m(20.0 * h).unwrap(),
        0.0,
    ];
    let points = [from, interpolate(from, to, 0.5), to];
    let gate = CancelGate::new();
    let bounds = law
        .admit_quadratic_path(&shear, points, 512, &gate)
        .unwrap();
    assert!(
        bounds.evaluations > 4,
        "opposite radius corner must require subdivision"
    );
    assert!(bounds.evaluations <= 512);
    for i in 0..=256 {
        let q = quadratic(points, i as f64 / 256.0);
        let response = law
            .evaluate_displacement_controlled(q[0], [q[1], q[2]])
            .unwrap();
        assert!(response.contact_radius_m >= bounds.minimum_radius_m);
        assert!(response.contact_radius_m <= bounds.maximum_radius_m);
        assert!(
            q[1].hypot(q[2]) / response.contact_radius_m
                <= bounds.maximum_shear_displacement_to_radius
        );
        assert!(
            response
                .normal_force_n
                .hypot(response.tangential_force_n[0])
                .hypot(response.tangential_force_n[1])
                <= bounds.maximum_force_n
        );
    }
    assert_eq!(
        law.admit_quadratic_path(&shear, points, bounds.evaluations - 1, &gate),
        Err(JkrMixedModePathError::Budget)
    );

    // Noncollinear motion reaches a larger contact radius between its ends;
    // the bound covers that actual area rather than holding endpoint area.
    let u = 0.6 * law.critical_tangential_displacement_m(2.0 * h).unwrap();
    let rotating = [
        [2.0 * h, u, 0.0],
        [2.0 * h, 0.0, 0.2 * u],
        [2.0 * h, -u, 0.0],
    ];
    let bounds = law
        .admit_quadratic_path(&shear, rotating, 512, &gate)
        .unwrap();
    let middle = quadratic(rotating, 0.5);
    let center = law
        .evaluate_displacement_controlled(middle[0], [middle[1], middle[2]])
        .unwrap();
    let end = law
        .evaluate_displacement_controlled(rotating[0][0], [u, 0.0])
        .unwrap();
    assert!(center.contact_area_m2 > end.contact_area_m2);
    assert!(bounds.maximum_radius_m >= center.contact_radius_m);
}

#[test]
fn g4_hidden_fold_capacity_source_excursions_and_cancellation_never_publish() {
    let (law, shear) = sourced_law(1e9);
    let saved = law.clone();
    let gate = CancelGate::new();
    let critical = law.critical_tangential_displacement_m(0.0).unwrap();
    let zero = [0.0; 3];
    let hidden_fold = [zero, [0.0, 3.0 * critical, 0.0], zero];
    assert!(matches!(
        law.admit_quadratic_path(&shear, hidden_fold, 512, &gate),
        Err(JkrMixedModePathError::Source(
            JkrMixedModeError::NoFullStickEquilibrium { .. }
        ))
    ));
    let maximum = law.normal_source().envelope().maximum_indentation_m;
    assert!(matches!(
        law.admit_quadratic_path(&shear, [zero, [3.0 * maximum, 0.0, 0.0], zero], 512, &gate),
        Err(JkrMixedModePathError::Source(JkrMixedModeError::Normal(_)))
    ));

    let u = 0.4 * critical;
    let response = law.evaluate_displacement_controlled(0.0, [u, 0.0]).unwrap();
    let capacity = 0.5 * response.tangential_force_n[0] / response.contact_area_m2;
    let (weak, weak_shear) = sourced_law(capacity);
    assert!(matches!(
        weak.admit_quadratic_path(&weak_shear, [zero, [0.0, 2.0 * u, 0.0], zero], 512, &gate),
        Err(JkrMixedModePathError::StaticCapacity { .. })
    ));
    assert_eq!(
        law.admit_quadratic_path(&weak_shear, [zero; 3], 512, &gate),
        Err(JkrMixedModePathError::StateMismatch)
    );
    assert_eq!(
        law.admit_quadratic_path(&shear, [zero; 3], 1, &gate),
        Err(JkrMixedModePathError::Budget)
    );
    let cancelled = CancelGate::new();
    cancelled.request();
    assert_eq!(
        law.discrete_path(zero, zero, &cancelled),
        Err(JkrMixedModePathError::Cancelled)
    );
    assert_eq!(
        law.admit_quadratic_path(&shear, [zero; 3], 512, &cancelled),
        Err(JkrMixedModePathError::Cancelled)
    );
    assert_eq!(law, saved);
}
