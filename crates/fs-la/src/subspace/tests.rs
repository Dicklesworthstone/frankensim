use super::*;

fn controls() -> SubspaceControls {
    SubspaceControls {
        maximum_dimension: 16,
        maximum_constraint_residual: 1e-12,
        maximum_orthogonality_error: 1e-12,
    }
}
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-13 * a.abs().max(b.abs()).max(1.0),
        "{a:e} != {b:e}"
    );
}
fn gear() -> LinearSubspace {
    let q = 1.0 / 2.0_f64.sqrt();
    LinearSubspace::try_new(
        &[0.1, 0.2],
        1,
        vec![0.1 * q, 0.2 * q],
        &[-2.0, 1.0],
        controls(),
        || false,
    )
    .unwrap()
}

#[test]
fn geared_mass_stiffness_and_damping_match_independent_analytic_energy() {
    let map = gear();
    let mass = map
        .symmetric_congruence(&[2.0, 0.0, 0.0, 3.0], || false)
        .unwrap();
    let stiffness = map.congruence(&[100.0, 0.0, 0.0, 200.0], || false).unwrap();
    let damping = map.congruence(&[7.0, 0.0, 0.0, 11.0], || false).unwrap();
    near(mass[0], 0.07);
    near(stiffness[0], 4.5);
    near(damping[0], 0.255);
    near(
        stiffness[0] / mass[0],
        (100.0 + 4.0 * 200.0) / (2.0 + 4.0 * 3.0),
    );
    let velocity = map.lift(&[0.7], || false).unwrap();
    near(
        0.5 * mass[0] * 0.7_f64.powi(2),
        0.5 * (2.0 * velocity[0].powi(2) + 3.0 * velocity[1].powi(2)),
    );
}

#[test]
fn physical_force_projection_is_work_conjugate_and_reactions_disappear() {
    let map = gear();
    let velocity = map.lift(&[-0.4], || false).unwrap();
    let force = [5.0, -7.0];
    let reduced = map.project(&force, || false).unwrap();
    near(
        reduced[0] * (-0.4),
        velocity[0] * force[0] + velocity[1] * force[1],
    );
    near(map.project(&[-10.0, 5.0], || false).unwrap()[0], 0.0);
    near(-2.0 * velocity[0] + velocity[1], 0.0);
}

#[test]
fn signed_row_rescaling_preserves_admission_and_original_map() {
    let q = 1.0 / 2.0_f64.sqrt();
    let expected = gear();
    for scale in [1e-150, -1.0, 1e150] {
        let actual = LinearSubspace::try_new(
            &[0.1, 0.2],
            1,
            vec![0.1 * q, 0.2 * q],
            &[-2.0 * scale, scale],
            controls(),
            || false,
        )
        .unwrap();
        assert_eq!(actual.basis(), expected.basis());
        assert_eq!(
            actual.project(&[3.0, 4.0], || false),
            expected.project(&[3.0, 4.0], || false)
        );
    }
}

#[test]
fn gyro_and_unstable_stiffness_are_not_symmetrized_or_clamped() {
    let map = LinearSubspace::try_new(
        &[1.0; 3],
        2,
        vec![1., 0., 0., 1., 0., 0.],
        &[0., 0., 1.],
        controls(),
        || false,
    )
    .unwrap();
    let c = map
        .congruence(&[4., -12., 99., 12., 5., 17., 8., 9., 6.], || false)
        .unwrap();
    assert_eq!(c, [4., -12., 12., 5.]);
    let k = map
        .congruence(&[-19.62, 0., 15., 0., 72., 11., 3., 4., 100.], || false)
        .unwrap();
    assert_eq!(k, [-19.62, 0., 0., 72.]);
    let v = [0.7, -0.4];
    near(
        v[0] * (c[0] * v[0] + c[1] * v[1]) + v[1] * (c[2] * v[0] + c[3] * v[1]),
        4.0 * v[0] * v[0] + 5.0 * v[1] * v[1],
    );
    assert!(
        map.symmetric_congruence(&[1., 2., 0., 0., 1., 0., 0., 0., 1.], || false)
            .is_err()
    );
}

#[test]
fn full_space_locked_space_and_empty_body_have_their_exact_meanings() {
    let full = LinearSubspace::try_new(&[2., 3.], 2, vec![2., 0., 0., 3.], &[], controls(), || {
        false
    })
    .unwrap();
    assert_eq!(full.lift(&[4., 5.], || false).unwrap(), [8., 15.]);
    assert_eq!(
        full.congruence(&[1., 2., 3., 4.], || false).unwrap(),
        [4., 12., 18., 36.]
    );
    let locked =
        LinearSubspace::try_new(&[2., 3.], 0, vec![], &[1., 0., 0., 1.], controls(), || {
            false
        })
        .unwrap();
    assert_eq!(locked.lift(&[], || false).unwrap(), [0., 0.]);
    assert!(locked.project(&[8., 9.], || false).unwrap().is_empty());
    assert!(
        locked
            .congruence(&[1., 2., 3., 4.], || false)
            .unwrap()
            .is_empty()
    );
    let empty = LinearSubspace::try_new(&[], 0, vec![], &[], controls(), || false).unwrap();
    assert!(empty.lift(&[], || false).unwrap().is_empty());
    assert!(
        empty
            .symmetric_congruence(&[], || false)
            .unwrap()
            .is_empty()
    );
}

#[test]
fn original_constraint_and_column_independence_are_independently_gated() {
    assert!(matches!(
        LinearSubspace::try_new(&[1., 1.], 1, vec![1., 0.], &[1., 0.], controls(), || false),
        Err(SubspaceError::Constraint { .. })
    ));
    assert!(matches!(
        LinearSubspace::try_new(&[1., 1.], 2, vec![1., 1., 0., 0.], &[], controls(), || {
            false
        }),
        Err(SubspaceError::Orthogonality { .. })
    ));
    assert!(
        LinearSubspace::try_new(&[1., 1.], 1, vec![1., 0.], &[0., 0.], controls(), || false)
            .is_err()
    );
}

#[test]
fn source_range_loss_is_refused_instead_of_erasing_work() {
    let map =
        LinearSubspace::try_new(&[1e-300], 1, vec![1e-300], &[], controls(), || false).unwrap();
    assert!(map.lift(&[1e-300], || false).is_err());
    assert!(map.project(&[1e-300], || false).is_err());
    assert!(map.congruence(&[1e-300], || false).is_err());
    assert!(
        LinearSubspace::try_new(&[1e300], 1, vec![1e300], &[1e300], controls(), || false).is_err()
    );
}

#[test]
fn explicit_dimension_and_cancellation_bounds_publish_no_partial_map() {
    let mut cap = controls();
    cap.maximum_dimension = 1;
    assert!(
        LinearSubspace::try_new(&[1., 1.], 2, vec![1., 0., 0., 1.], &[], cap, || false).is_err()
    );
    assert_eq!(
        LinearSubspace::try_new(&[1.], 1, vec![1.], &[], controls(), || true),
        Err(SubspaceError::Cancelled)
    );
    let map = gear();
    let before = map.clone();
    assert_eq!(
        map.project(&[3., 4.], || true),
        Err(SubspaceError::Cancelled)
    );
    assert_eq!(map, before);
    assert!(map.project(&[3., 4.], || false).is_ok());
}
