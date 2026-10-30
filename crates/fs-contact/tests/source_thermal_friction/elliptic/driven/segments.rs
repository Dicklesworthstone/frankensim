//! G1/G4: source trials stay private; accepted physical segments retain work.
use super::*;
use fs_contact::interface_binding::thermal_friction::elliptic::driven::SourceDrivenSegmentError;
use std::cell::Cell;
use std::convert::Infallible;

#[test]
fn speed_dependent_breakaway_visits_only_the_accepted_original_source_trajectory() {
    let card = motion_card([0.2, 0.6], [0.2, 0.4]);
    let law = motion_law(&card, [1.0, 0.0]);
    let initial = law.initial_state();
    let input = drive([0.0; 2], [13.0, 0.0], EllipticDrivenMode::Resting);
    // F=2+20r amplifies distance-root error into the force assertion below.
    let mut motion = solve([0.0, 0.2]);
    motion.absolute_distance_tolerance_m = 1e-14;
    motion.relative_distance_tolerance = 0.0;
    let mut native_controls = driven_controls();
    native_controls.solve.force_absolute_tolerance_n = 1e-13;
    native_controls.solve.force_relative_tolerance = 1e-15;
    let mut segments = Vec::new();
    let result = law
        .advance_driven_motion_consistent_with_segments(
            &initial,
            &input,
            &endpoint(),
            motion,
            native_controls,
            &mut || false,
            &mut |segment| {
                segments.push(segment.clone());
                Ok::<_, Infallible>(())
            },
        )
        .unwrap();
    assert!(result.evaluations > 1);
    assert_eq!(
        segments.len(),
        1,
        "source root trials must never visit poses"
    );
    let step = &segments[0];
    close(step.elapsed_before_s, 0.0);
    close(step.dt_s, 0.1);
    close(step.traction_n[0], 3.0);
    close(step.velocity_after_m_s[0], 1.0);
    close(step.slip_increment_m[0], 0.05);
    close(step.distance_m, 0.05);
    close(step.dissipated_j, 0.15);
    assert_eq!(step.mode_before, EllipticDrivenMode::Resting);
    assert_eq!(step.mode_after, EllipticDrivenMode::Sliding);
    close(result.driven.drive_work_j, 0.65);
    assert_eq!(result.state.constitutive_state(), &result.driven.state);
    assert_eq!(initial, law.initial_state());
    for (axis, source) in result.state.friction_sources().iter().enumerate() {
        close(source.query().axes()["speed"], 0.5);
        close(source.query().axes()["path"], 0.3);
        for property in source.properties() {
            assert_eq!(
                property.answer().receipt.selected,
                initial.friction_sources()[axis]
                    .property(property.requirement().name())
                    .unwrap()
                    .answer()
                    .receipt
                    .selected,
            );
            card.claims()
                .verify_receipt(&property.answer().receipt)
                .unwrap();
        }
    }
}

#[test]
fn zero_drive_curved_source_queries_close_the_same_segmented_path_that_moves_bodies() {
    let law = motion_law(&motion_card([0.3; 2], [0.6; 2]), [0.6, 0.8]);
    let initial = law.initial_state();
    let native = EllipticCoulomb::new([1.2; 2], [0.3, 0.6], [0.6, 0.8]).unwrap();
    let mut input = drive([0.6, 0.8], [0.0; 2], EllipticDrivenMode::Sliding);
    input.mobility_per_kg = [[1.0, 0.35], [0.35, 2.0]];
    input.dt_s = 0.05;
    let expected = native
        .advance_driven_segments(
            &native.initial_state(),
            &input,
            driven_controls(),
            &mut || false,
            &mut |_| Ok::<_, Infallible>(()),
        )
        .unwrap();
    let mut time = 0.0;
    let mut velocity = input.velocity_m_s;
    let mut impulse = [0.0; 2];
    let mut displacement = [0.0; 2];
    let mut path = 0.0;
    let mut heat = 0.0;
    let mut visits = 0;
    let result = law
        .advance_driven_motion_consistent_with_segments(
            &initial,
            &input,
            &endpoint(),
            solve([0.0, 0.1]),
            driven_controls(),
            &mut || false,
            &mut |segment| {
                close(segment.elapsed_before_s, time);
                for i in 0..2 {
                    close(segment.velocity_before_m_s[i], velocity[i]);
                    let recoil = input.mobility_per_kg[i][0] * segment.traction_n[0]
                        + input.mobility_per_kg[i][1] * segment.traction_n[1];
                    let after = velocity[i] - segment.dt_s * recoil;
                    let slip = 0.5 * segment.dt_s * (velocity[i] + after);
                    close(segment.velocity_after_m_s[i], after);
                    close(segment.slip_increment_m[i], slip);
                    displacement[i] += slip;
                    impulse[i] += segment.dt_s * segment.traction_n[i];
                    heat += segment.traction_n[i] * slip;
                    velocity[i] = after;
                }
                time += segment.dt_s;
                path += segment.distance_m;
                visits += 1;
                Ok::<_, Infallible>(())
            },
        )
        .unwrap();
    assert!(visits > 1);
    close(time, input.dt_s);
    close(result.actual_distance_m, expected.distance_m);
    close(path, result.actual_distance_m);
    close(heat, result.driven.dissipated_j);
    close(result.driven.dissipated_j, expected.dissipated_j);
    for i in 0..2 {
        close(velocity[i], result.driven.velocity_after_m_s[i]);
        close(impulse[i], result.driven.traction_impulse_ns[i]);
        close(displacement[i], result.driven.slip_increment_m[i]);
    }
    assert!(path > norm(displacement));
    close(
        result.state.friction_sources()[0].query().axes()["speed"],
        path / input.dt_s,
    );
    close(
        result.state.friction_sources()[0].query().axes()["path"],
        0.25 + path,
    );
    close(result.driven.energy_residual_j, 0.0);
    assert_eq!(initial, law.initial_state());
}

#[test]
fn source_refusals_never_visit_and_failed_replays_leave_original_history_retryable() {
    let law = motion_law(&motion_card([0.3; 2], [0.6; 2]), [0.6, 0.8]);
    let initial = law.initial_state();
    let mut input = drive([0.6, 0.8], [0.4, -0.2], EllipticDrivenMode::Sliding);
    input.mobility_per_kg = [[1.0, 0.35], [0.35, 2.0]];
    input.dt_s = 0.05;
    for failure in 0..4 {
        let mut motion = solve([0.0, 0.1]);
        let mut query = endpoint();
        match failure {
            0 => motion.maximum_evaluations = 1,
            1 => query = query.with("ambient", 296.0).unwrap(),
            _ => {}
        }
        let visits = Cell::new(0);
        let cancel = Cell::new(false);
        let result = law.advance_driven_motion_consistent_with_segments(
            &initial,
            &input,
            &query,
            motion,
            driven_controls(),
            &mut || cancel.get(),
            &mut |_| {
                visits.set(visits.get() + 1);
                if failure == 2 {
                    return Err("pose owner refused");
                }
                if failure == 3 {
                    cancel.set(true);
                }
                Ok(())
            },
        );
        assert!(result.is_err());
        match failure {
            0 | 1 => assert_eq!(visits.get(), 0),
            2 => {
                assert_eq!(visits.get(), 1);
                assert!(matches!(
                    result,
                    Err(SourceDrivenSegmentError::Segments(
                        EllipticDrivenSegmentError::Visitor("pose owner refused")
                    ))
                ));
            }
            3 => {
                assert_eq!(visits.get(), 1);
                assert!(matches!(
                    result,
                    Err(SourceDrivenSegmentError::Segments(
                        EllipticDrivenSegmentError::Driven(EllipticDrivenError::Cancelled)
                    ))
                ));
            }
            _ => unreachable!(),
        }
        assert_eq!(initial, law.initial_state());
    }
    let result = law
        .advance_driven_motion_consistent_with_segments(
            &initial,
            &input,
            &endpoint(),
            solve([0.0, 0.1]),
            driven_controls(),
            &mut || false,
            &mut |_| Ok::<_, Infallible>(()),
        )
        .unwrap();
    assert!(result.actual_distance_m > 0.0);
    assert!(result.driven.dissipated_j > 0.0);
}
