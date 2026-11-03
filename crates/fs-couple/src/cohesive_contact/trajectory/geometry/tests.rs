//! Geometry and time-retry identities for the actual affine contact adapter.
use super::*;

fn row<'a>(id: u64, gap_m: f64, trace: [&'a [TraceTerm]; 3]) -> FrictionalContact<'a> {
    FrictionalContact {
        id,
        gap_m,
        activation_tolerance_m: 1e-6,
        friction_coefficient: 0.0,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.5 },
        prescribed_velocity_m_s: [0.0; 3],
        trace,
    }
}
fn one_body() -> [[TraceTerm; 1]; 3] {
    core::array::from_fn(|dof| {
        [TraceTerm {
            body: 0,
            dof,
            coefficient: 1.0,
        }]
    })
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-13 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}
fn admitted(contacts: &[FrictionalContact<'_>], reference: &[&[f64]]) {
    admit(contacts, reference, 128, 4096, &mut || false).unwrap();
}

#[test]
fn displacement_differences_precede_compensated_trace_sums() {
    let normal = core::array::from_fn::<_, 3, _>(|dof| TraceTerm {
        body: 0,
        dof,
        coefficient: 1.0,
    });
    let contact = row(7, 0.25, [&normal, &[], &[]]);
    let reference = [1e16, 0.0, -1e16];
    let position = [1e16, 1.0, -1e16];
    admitted(&[contact], &[&reference]);
    let sampled = sample(
        &[contact],
        &[&reference],
        &[&position],
        &[&[0.0; 3]],
        0.0,
        &mut || false,
    )
    .unwrap();
    close(sampled.gaps_m[0], 1.25);
    // Forming the complete current trace first would lose the unit change.
    assert_eq!((position[0] + position[1]) + position[2], 0.0);

    let reference = [0.0; 3];
    admitted(&[contact], &[&reference]);
    let sampled = sample(
        &[contact],
        &[&reference],
        &[&position],
        &[&[1e16, 2.0, -1e16]],
        0.0,
        &mut || false,
    )
    .unwrap();
    close(sampled.gaps_m[0], 1.25);
    close(sampled.normal_velocities_m_s[0], 2.0);
}

#[test]
fn cross_body_traces_include_affine_support_motion_in_gap_and_speed() {
    let rows: [[TraceTerm; 2]; 3] = core::array::from_fn(|dof| {
        [
            TraceTerm {
                body: 0,
                dof,
                coefficient: -1.0,
            },
            TraceTerm {
                body: 1,
                dof,
                coefficient: 1.0,
            },
        ]
    });
    let mut contact = row(11, 0.4, [&rows[0], &rows[1], &rows[2]]);
    contact.prescribed_velocity_m_s = [-2.0, 0.25, -0.5];
    let reference = [[3.0, 4.0, 5.0], [7.0, 8.0, 9.0]];
    let positions = [[3.1, 4.2, 5.0], [7.5, 7.9, 9.0]];
    let velocities = [[2.0, 1.0, 0.0], [-1.0, 0.5, 0.0]];
    admitted(&[contact], &[&reference[0], &reference[1]]);
    let sampled = sample(
        &[contact],
        &[&reference[0], &reference[1]],
        &[&positions[0], &positions[1]],
        &[&velocities[0], &velocities[1]],
        0.1,
        &mut || false,
    )
    .unwrap();
    close(sampled.gaps_m[0], 0.6);
    close(sampled.normal_velocities_m_s[0], -5.0);
    close(sampled.minimum_gap_m, 0.6);
    assert_eq!(sampled.closing_contact(&[contact], 1e-9), None);
}

#[test]
fn nonnegative_endpoints_admit_the_entire_declared_affine_path_without_gap_projection() {
    let rows = one_body();
    let mut contact = row(3, 1.0, [&rows[0], &rows[1], &rows[2]]);
    contact.prescribed_velocity_m_s[0] = -0.1;
    let reference = [0.0; 3];
    admitted(&[contact], &[&reference]);
    for fraction in [0.0, 0.25, 0.5, 0.75, 1.0] {
        let position = [-0.8 * fraction, 0.0, 0.0];
        let sampled = sample(
            &[contact],
            &[&reference],
            &[&position],
            &[&[-0.8, 0.0, 0.0]],
            fraction,
            &mut || false,
        )
        .unwrap();
        close(sampled.gaps_m[0], 1.0 - 0.9 * fraction);
        assert!(sampled.minimum_gap_m >= 0.0);
    }
    let penetrating = sample(
        &[contact],
        &[&reference],
        &[&[-1.0, 0.0, 0.0]],
        &[&[-1.0, 0.0, 0.0]],
        1.0,
        &mut || false,
    )
    .unwrap();
    close(penetrating.gaps_m[0], -0.1);
    assert!(penetrating.minimum_gap_m < 0.0);
    assert_eq!(penetrating.closing_contact(&[contact], 0.0), None);

    let mut initially_penetrating = contact;
    initially_penetrating.gap_m = -1e-14;
    assert!(matches!(
        admit(
            &[initially_penetrating],
            &[&reference],
            128,
            4096,
            &mut || false
        ),
        Err(GeometryError::Penetration { contact_id: 3, .. })
    ));
}

#[test]
fn active_closing_selection_uses_actual_nonnegative_gaps_and_stable_ids() {
    let rows = one_body();
    let contacts = [
        row(20, 0.5e-6, [&rows[0], &rows[1], &rows[2]]),
        row(5, 0.0, [&rows[0], &rows[1], &rows[2]]),
        row(1, 0.1, [&rows[0], &rows[1], &rows[2]]),
    ];
    let reference = [0.0; 3];
    admitted(&contacts, &[&reference]);
    let sampled = sample(
        &contacts,
        &[&reference],
        &[&reference],
        &[&[-1.0, 0.0, 0.0]],
        0.0,
        &mut || false,
    )
    .unwrap();
    assert_eq!(sampled.closing_contact(&contacts, 0.9), Some(1));
    assert_eq!(sampled.closing_contact(&contacts, 1.0), None);
    let reversed = [contacts[2], contacts[1], contacts[0]];
    let replay = sample(
        &reversed,
        &[&reference],
        &[&reference],
        &[&[-1.0, 0.0, 0.0]],
        0.0,
        &mut || false,
    )
    .unwrap();
    assert_eq!(
        reversed[replay.closing_contact(&reversed, 0.0).unwrap()].id,
        5
    );
    let departing = sample(
        &contacts,
        &[&reference],
        &[&reference],
        &[&[1.0, 0.0, 0.0]],
        0.0,
        &mut || false,
    )
    .unwrap();
    assert_eq!(departing.closing_contact(&contacts, 0.0), None);
}

#[test]
fn all_three_trace_rows_and_body_inputs_are_admitted_before_motion_work() {
    let rows = one_body();
    let contact = row(1, 0.2, [&rows[0], &rows[1], &rows[2]]);
    let reference = [0.0; 3];
    admitted(&[contact], &[&reference]);
    for component in 0..3 {
        for bad in [
            TraceTerm {
                body: 1,
                dof: 0,
                coefficient: 1.0,
            },
            TraceTerm {
                body: 0,
                dof: 3,
                coefficient: 1.0,
            },
            TraceTerm {
                body: 0,
                dof: 0,
                coefficient: 0.0,
            },
            TraceTerm {
                body: 0,
                dof: 0,
                coefficient: f64::NAN,
            },
        ] {
            let bad_row = [bad];
            let mut malformed = contact;
            malformed.trace[component] = &bad_row;
            assert!(matches!(
                admit(&[malformed], &[&reference], 128, 4096, &mut || false),
                Err(GeometryError::Invalid(_))
            ));
        }
        for bad_row in [[rows[0][0], rows[0][0]], [rows[1][0], rows[0][0]]] {
            let mut malformed = contact;
            malformed.trace[component] = &bad_row;
            assert!(matches!(
                admit(&[malformed], &[&reference], 128, 4096, &mut || false),
                Err(GeometryError::Invalid(_))
            ));
        }
        let mut malformed = contact;
        malformed.prescribed_velocity_m_s[component] = f64::INFINITY;
        assert!(matches!(
            admit(&[malformed], &[&reference], 128, 4096, &mut || false),
            Err(GeometryError::Invalid(_))
        ));
    }
    assert!(matches!(
        admit(&[contact, contact], &[&reference], 128, 4096, &mut || false),
        Err(GeometryError::Invalid(_))
    ));
    assert!(matches!(
        admit(&[contact], &[&[0.0, f64::NAN, 0.0]], 128, 4096, &mut || {
            false
        }),
        Err(GeometryError::NonFinite(_))
    ));
    assert!(matches!(
        sample(
            &[contact],
            &[&reference],
            &[&[0.0; 2]],
            &[&reference],
            0.0,
            &mut || false
        ),
        Err(GeometryError::Invalid(_))
    ));
    assert!(matches!(
        sample(
            &[contact],
            &[&reference],
            &[&reference],
            &[&[0.0, 0.0, f64::NAN]],
            0.0,
            &mut || false
        ),
        Err(GeometryError::NonFinite(_))
    ));
}

#[test]
fn retry_duration_requires_recomputed_native_geometry_not_interpolated_state() {
    let rows = one_body();
    let contact = row(9, 0.4, [&rows[0], &rows[1], &rows[2]]);
    let reference = [0.0; 3];
    admitted(&[contact], &[&reference]);
    // An independent constant-acceleration trajectory has
    // g(h)=.4-h-5h^2, and hits at h=.2. The geometry helper owns no integrator.
    let h = 0.3;
    let old_position = [-h - 5.0 * h * h, 0.0, 0.0];
    let old_velocity = [-1.0 - 10.0 * h, 0.0, 0.0];
    let rejected = sample(
        &[contact],
        &[&reference],
        &[&old_position],
        &[&old_velocity],
        h,
        &mut || false,
    )
    .unwrap();
    close(rejected.gaps_m[0], -0.35);
    let shorter = shorter_step(&[contact], &[0.4], &rejected.gaps_m, h).unwrap();
    assert!(shorter > 0.0 && shorter <= 0.5 * h);
    let position = [-shorter - 5.0 * shorter * shorter, 0.0, 0.0];
    let velocity = [-1.0 - 10.0 * shorter, 0.0, 0.0];
    let recomputed = sample(
        &[contact],
        &[&reference],
        &[&position],
        &[&velocity],
        shorter,
        &mut || false,
    )
    .unwrap();
    close(
        recomputed.gaps_m[0],
        0.4 - shorter - 5.0 * shorter * shorter,
    );
    assert!(recomputed.minimum_gap_m >= 0.0);
    let interpolated_gap = 0.4 + (rejected.gaps_m[0] - 0.4) * shorter / h;
    assert!((recomputed.gaps_m[0] - interpolated_gap).abs() > 1e-3);
    close(
        shorter_step(&[contact], &[f64::MAX], &[-f64::MAX], 1.0).unwrap(),
        0.5,
    );
    // Near-contact retries still shorten; they never create a zero or enlarged step.
    close(shorter_step(&[contact], &[0.0], &[-1.0], 1.0).unwrap(), 0.5);
    assert!(matches!(
        shorter_step(&[contact], &[0.0], &[-1.0], f64::from_bits(1)),
        Err(GeometryError::NoProgress { .. })
    ));
    assert!(shorter_step(&[contact], &[0.4], &[0.1], h).is_err());
}

#[test]
fn work_bounds_bad_clocks_and_cancellation_return_no_geometry_candidate() {
    let rows = one_body();
    let contact = row(1, 0.2, [&rows[0], &rows[1], &rows[2]]);
    let reference = [0.0; 3];
    assert!(matches!(
        admit(&[contact], &[&reference], 2, 4096, &mut || false),
        Err(GeometryError::Budget(_))
    ));
    assert!(matches!(
        admit(&[contact], &[&reference], 128, 0, &mut || false),
        Err(GeometryError::Budget(_))
    ));
    assert!(matches!(
        admit(&[contact], &[&reference], 128, 4096, &mut || true),
        Err(GeometryError::Cancelled)
    ));
    for elapsed in [-1.0, f64::INFINITY, f64::NAN] {
        assert!(matches!(
            sample(
                &[contact],
                &[&reference],
                &[&reference],
                &[&reference],
                elapsed,
                &mut || false
            ),
            Err(GeometryError::Invalid(_))
        ));
    }
    let mut polls = 0usize;
    sample(
        &[contact],
        &[&reference],
        &[&reference],
        &[&reference],
        0.0,
        &mut || {
            polls += 1;
            false
        },
    )
    .unwrap();
    for stop in [1, polls / 2, polls] {
        let mut observed = 0;
        assert!(matches!(
            sample(
                &[contact],
                &[&reference],
                &[&reference],
                &[&reference],
                0.0,
                &mut || {
                    observed += 1;
                    observed == stop
                }
            ),
            Err(GeometryError::Cancelled)
        ));
    }
    let empty = sample(
        &[],
        &[&reference],
        &[&reference],
        &[&reference],
        0.0,
        &mut || false,
    )
    .unwrap();
    assert_eq!(empty.minimum_gap_m, 0.0);
    assert!(empty.gaps_m.is_empty());
}
