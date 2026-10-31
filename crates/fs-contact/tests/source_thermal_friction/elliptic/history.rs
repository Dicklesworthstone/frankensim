//! Independently admitted motion updates original slip/source history once.
use super::*;
use fs_contact::interface_binding::thermal_friction::elliptic::history::*;
use fs_contact::interface_binding::thermal_friction::elliptic::impulse::{
    SourceFrictionConditions, SourceImpulseCoefficientRole,
};
use fs_tribo::tangential::transport::{ObjectiveTangentialLaw, TangentialHistoryRotation};

fn history_controls() -> SourcePostcontactMotionControls {
    SourcePostcontactMotionControls {
        absolute_slip_tolerance_m: 1e-13,
        absolute_path_tolerance_m: 1e-13,
    }
}

fn boundary() -> SourceFrictionConditions {
    SourceFrictionConditions {
        query: endpoint(),
        contact_temperature_k: Some(320.0),
        nominal_pressure_pa: Some(2e5),
    }
}

#[test]
fn finite_motion_after_impact_refreshes_pinned_sources_and_unsigned_path_through_reversal() {
    let card = motion_card([0.3, 0.5], [0.2, 0.4]);
    let law = motion_law(&card, [1.0, 0.0]);
    let initial = law.initial_state();
    let sampled = law
        .sample_preimpact_coefficients(
            &initial,
            [2.0, 0.0],
            &boundary(),
            SourceImpulseCoefficientRole::Kinetic,
            &mut || false,
        )
        .unwrap();
    // The instantaneous stage still owns no travel, temperature evolution or
    // clock advancement. Only the admitted finite phase below records slip.
    assert_eq!(
        sampled.state.constitutive_state(),
        initial.constitutive_state()
    );
    close(
        sampled.state.friction_sources()[0].query().axes()["path"],
        0.25,
    );
    let mut state = sampled.state;
    for (velocity, expected_slip, expected_path) in [
        ([0.6, 0.8], [0.3, 0.4], 0.75),
        ([-0.6, -0.8], [0.0, 0.0], 1.25),
    ] {
        let saved = state.clone();
        let report = law
            .accept_postcontact_motion(
                &state,
                velocity,
                0.5,
                &boundary(),
                history_controls(),
                &mut || false,
            )
            .unwrap();
        close(report.duration_s, 0.5);
        close(report.path_increment_m, 0.5);
        close(report.queried_speed_m_s, 1.0);
        for axis in 0..2 {
            close(report.slip_increment_m[axis], 0.5 * velocity[axis]);
            close(
                report.state.constitutive_state().slip_m()[axis],
                expected_slip[axis],
            );
            close(report.slip_residual_m[axis], 0.0);
        }
        close(report.maximum_path_coordinate_residual_m, 0.0);
        assert_eq!(
            report.state.constitutive_state().principal_axis(),
            [1.0, 0.0]
        );
        assert_eq!(
            report.state.constitutive_state().contact_temperature_k(),
            Some(320.0)
        );
        close(
            report.state.constitutive_law().kinetic_coefficients()[0],
            0.4,
        );
        close(
            report.state.constitutive_law().kinetic_coefficients()[1],
            0.2 + 0.1 * expected_path,
        );
        for (axis, source) in report.state.friction_sources().iter().enumerate() {
            let query = source.query();
            close(query.axes()["speed"], 1.0);
            close(query.axes()["path"], expected_path);
            for name in ["T", "Ta", "Tb", "p", "ambient"] {
                assert_eq!(query.axes()[name], endpoint().axes()[name]);
            }
            for property in source.properties() {
                let original = initial.friction_sources()[axis]
                    .property(property.requirement().name())
                    .unwrap();
                assert_eq!(
                    property.answer().receipt.selected,
                    original.answer().receipt.selected
                );
                card.claims()
                    .verify_receipt(&property.answer().receipt)
                    .unwrap();
            }
        }
        assert_eq!(state, saved);
        state = report.state;
    }
    let stopped = law
        .accept_postcontact_motion(
            &state,
            [0.0; 2],
            0.25,
            &boundary(),
            history_controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(stopped.state.constitutive_state().slip_m(), [0.0; 2]);
    assert_eq!(stopped.path_increment_m, 0.0);
    assert_eq!(stopped.queried_speed_m_s, 0.0);
    close(
        stopped.state.friction_sources()[0].query().axes()["path"],
        1.25,
    );
    assert_eq!(initial, law.initial_state());
}

#[test]
fn held_motion_subdivision_and_objective_gauge_change_preserve_physical_history() {
    let law = motion_law(&motion_card([0.3, 0.5], [0.2, 0.4]), [1.0, 0.0]);
    let initial = law.initial_state();
    let whole = law
        .accept_postcontact_motion(
            &initial,
            [0.5, 0.0],
            0.5,
            &boundary(),
            history_controls(),
            &mut || false,
        )
        .unwrap();
    let half = law
        .accept_postcontact_motion(
            &initial,
            [0.5, 0.0],
            0.25,
            &boundary(),
            history_controls(),
            &mut || false,
        )
        .unwrap();
    let second = law
        .accept_postcontact_motion(
            &half.state,
            [0.5, 0.0],
            0.25,
            &boundary(),
            history_controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(whole.state, second.state);
    close(
        whole.path_increment_m,
        half.path_increment_m + second.path_increment_m,
    );
    let old_frame = TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap();
    let new_frame = TangentFrame::new([0.0, 0.0, 1.0], [0.0, 1.0, 0.0]).unwrap();
    let rotation = TangentialHistoryRotation::between(old_frame, new_frame, 0.0).unwrap();
    let changed = law.transport_state(&initial, rotation).unwrap();
    let in_new_frame = law
        .accept_postcontact_motion(
            &changed,
            rotation.apply([0.5, 0.0]).unwrap(),
            0.5,
            &boundary(),
            history_controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(
        law.transport_state(&whole.state, rotation).unwrap(),
        in_new_frame.state
    );
    assert_eq!(
        in_new_frame.state.constitutive_state().slip_m(),
        [0.0, -0.25]
    );
    assert_eq!(
        in_new_frame.state.constitutive_state().principal_axis(),
        [0.0, -1.0]
    );
    close(in_new_frame.path_increment_m, whole.path_increment_m);
    assert_eq!(
        in_new_frame.state.friction_sources(),
        whole.state.friction_sources()
    );
}

#[test]
fn source_boundary_late_refusal_and_cancellation_preserve_the_accepted_checkpoint() {
    let card = motion_card([0.3, 0.5], [0.2, 0.4]);
    let law = motion_law(&card, [1.0, 0.0]);
    let state = law.initial_state();
    let saved = state.clone();
    let call = |velocity, dt, boundary: &SourceFrictionConditions, controls| {
        law.accept_postcontact_motion(&state, velocity, dt, boundary, controls, &mut || false)
    };
    for dt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(call([0.5, 0.0], dt, &boundary(), history_controls()).is_err());
    }
    for velocity in [[f64::NAN, 0.0], [f64::INFINITY, 0.0], [f64::MAX, f64::MAX]] {
        assert!(call(velocity, 1.0, &boundary(), history_controls()).is_err());
    }
    assert!(call([1e-300, 0.0], 1e-300, &boundary(), history_controls()).is_err());
    for bad in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(
            call(
                [0.5, 0.0],
                0.5,
                &boundary(),
                SourcePostcontactMotionControls {
                    absolute_path_tolerance_m: bad,
                    ..history_controls()
                }
            )
            .is_err()
        );
        assert!(
            call(
                [0.5, 0.0],
                0.5,
                &boundary(),
                SourcePostcontactMotionControls {
                    absolute_slip_tolerance_m: bad,
                    ..history_controls()
                }
            )
            .is_err()
        );
    }
    let mut wrong = boundary();
    wrong.nominal_pressure_pa = Some(3e5);
    assert!(call([0.5, 0.0], 0.5, &wrong, history_controls()).is_err());
    wrong = boundary();
    wrong.query = wrong.query.with("ambient", 296.0).unwrap();
    assert!(call([0.5, 0.0], 0.5, &wrong, history_controls()).is_err());
    // First principal's speed=1 is valid; the second principal's requested
    // path=2.25 lies beyond its ORIGINAL claim. No first-source update leaks.
    assert!(matches!(
        call([1.0, 0.0], 2.0, &boundary(), history_controls()),
        Err(SourcePostcontactMotionError::Source(
            SourceThermalCoulombError::Material(_)
        ))
    ));
    let foreign = motion_law(&card, [0.0, 1.0]);
    assert!(
        foreign
            .accept_postcontact_motion(
                &state,
                [0.5, 0.0],
                0.5,
                &boundary(),
                history_controls(),
                &mut || false,
            )
            .is_err()
    );
    assert!(matches!(
        law.accept_postcontact_motion(
            &state,
            [0.5, 0.0],
            0.5,
            &boundary(),
            history_controls(),
            &mut || true,
        ),
        Err(SourcePostcontactMotionError::Cancelled)
    ));
    let mut polls = 0;
    assert!(matches!(
        law.accept_postcontact_motion(
            &state,
            [0.5, 0.0],
            0.5,
            &boundary(),
            history_controls(),
            &mut || {
                polls += 1;
                polls == 8
            },
        ),
        Err(SourcePostcontactMotionError::Cancelled)
    ));
    assert_eq!(state, saved);
    assert_eq!(state, law.initial_state());
}

#[test]
fn source_path_cannot_silently_lose_or_misrepresent_a_resolved_increment() {
    let card = directional_card();
    let principal = [[STATIC, KINETIC], [STATIC_Y, KINETIC_Y]].map(|names| {
        SourceThermalCoulomb::resolve_with_query_axes(
            bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap(),
            point().with("path", 1e16).unwrap(),
            requirement(names[0]),
            requirement(names[1]),
            MaterialPropertySelection::SingleClaimOnly,
            ThermalFrictionAxes {
                temperature: "T".into(),
                nominal_pressure: Some("p".into()),
            },
            [
                ("Ta", FrictionQueryAxisKind::Temperature),
                ("Tb", FrictionQueryAxisKind::Temperature),
                ("speed", FrictionQueryAxisKind::SlidingSpeed),
                ("path", FrictionQueryAxisKind::AccumulatedDistance),
            ]
            .into_iter()
            .map(|(name, kind)| FrictionQueryAxis {
                name: name.into(),
                kind,
            })
            .collect(),
        )
        .unwrap()
    });
    let law = SourceThermalEllipticCoulomb::from_principal_sources(principal, [1.0, 0.0]).unwrap();
    let state = law.initial_state();
    assert!(
        law.accept_postcontact_motion(
            &state,
            [0.5, 0.0],
            1.0,
            &boundary(),
            SourcePostcontactMotionControls {
                absolute_path_tolerance_m: 10.0,
                ..history_controls()
            },
            &mut || false,
        )
        .is_err()
    );
    assert!(matches!(
        law.accept_postcontact_motion(
            &state,
            [1.5, 0.0],
            1.0,
            &boundary(),
            SourcePostcontactMotionControls {
                absolute_path_tolerance_m: 0.25,
                ..history_controls()
            },
            &mut || false,
        ),
        Err(SourcePostcontactMotionError::PathResolution {
            principal_axis: 0,
            requested_increment_m: 1.5,
            represented_increment_m: 2.0,
            ..
        })
    ));
    let explicit_roundoff = law
        .accept_postcontact_motion(
            &state,
            [1.5, 0.0],
            1.0,
            &boundary(),
            SourcePostcontactMotionControls {
                absolute_path_tolerance_m: 0.5,
                ..history_controls()
            },
            &mut || false,
        )
        .unwrap();
    assert_eq!(explicit_roundoff.maximum_path_coordinate_residual_m, 0.5);
    assert_eq!(
        explicit_roundoff.state.constitutive_state().slip_m(),
        [1.5, 0.0]
    );
    assert_eq!(state, law.initial_state());
}
