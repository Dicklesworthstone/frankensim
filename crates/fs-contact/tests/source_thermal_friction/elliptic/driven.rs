//! G1/G3/G4: original source claims drive static hold, release and physical work.
use super::*;
use fs_tribo::tangential::elliptic::EllipticCoulomb;
use fs_tribo::tangential::elliptic::driven::*;
use fs_tribo::tangential::elliptic::relaxation::EllipticRelaxationControls;

fn driven_controls() -> EllipticDrivenControls {
    EllipticDrivenControls {
        maximum_step_s: 0.005,
        minimum_step_s: 1e-10,
        maximum_substeps: 20_000,
        rest_velocity_tolerance_m_s: 1e-12,
        solve: TangentialSolveConfig {
            max_iterations: 128,
            maximum_force_n: 100.0,
            force_absolute_tolerance_n: 1e-11,
            force_relative_tolerance: 1e-13,
        },
        relaxation: EllipticRelaxationControls {
            time_tolerance_s: 1e-11,
            travel_tolerance_m: 1e-11,
            maximum_evaluations: 200_000,
            maximum_root_iterations: 96,
        },
    }
}
fn drive(
    velocity: [f64; 2],
    acceleration: [f64; 2],
    mode: EllipticDrivenMode,
) -> EllipticDrivenInput {
    EllipticDrivenInput {
        velocity_m_s: velocity,
        driving_acceleration_m_s2: acceleration,
        mode,
        mobility_per_kg: [[1.0, 0.0], [0.0, 2.0]],
        normal_n: 10.0,
        dt_s: 0.1,
        maximum_force_n: 100.0,
        contact_temperature_k: Some(320.0),
        nominal_pressure_pa: Some(2e5),
    }
}

#[test]
fn accepted_static_mode_holds_between_limits_and_sourced_breakaway_closes_speed_work() {
    let card = motion_card([0.2, 0.6], [0.2, 0.4]);
    let law = motion_law(&card, [1.0, 0.0]);
    let initial = law.initial_state();
    let held = law
        .advance_driven_motion_consistent(
            &initial,
            &drive([0.0; 2], [10.0, 0.0], EllipticDrivenMode::Resting),
            &endpoint(),
            solve([0.0, 0.2]),
            driven_controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(held.driven.mode, EllipticDrivenMode::Resting);
    assert_eq!(held.actual_distance_m, 0.0);
    assert_eq!(held.driven.dissipated_j, 0.0);
    assert_eq!(held.driven.drive_work_j, 0.0);
    assert_eq!(held.evaluations, 1);
    close(held.driven.traction_impulse_ns[0], 1.0);
    close(
        held.state.friction_sources()[0].query().axes()["path"],
        0.25,
    );
    close(
        held.state.friction_sources()[0].query().axes()["speed"],
        0.0,
    );

    // mu_k=.2+.2*(r/dt). At rest with a=13 and W=1:
    // r=.5*.1²*(13-2-20r), hence r=.05 and v_end=1.
    let moving = law
        .advance_driven_motion_consistent(
            &held.state,
            &drive([0.0; 2], [13.0, 0.0], held.driven.mode),
            &endpoint(),
            solve([0.0, 0.2]),
            driven_controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(moving.driven.mode, EllipticDrivenMode::Sliding);
    close(moving.actual_distance_m, 0.05);
    close(moving.driven.velocity_after_m_s[0], 1.0);
    close(moving.driven.dissipated_j, 0.15);
    close(moving.driven.drive_work_j, 0.65);
    close(moving.driven.kinetic_energy_change_j, 0.5);
    close(moving.driven.energy_residual_j, 0.0);
    assert!(moving.distance_residual_m <= 0.0);
    for (i, source) in moving.state.friction_sources().iter().enumerate() {
        close(source.query().axes()["speed"], 0.5);
        close(source.query().axes()["path"], 0.3);
        for name in ["T", "Ta", "Tb", "p", "ambient"] {
            assert_eq!(source.query().axes()[name], endpoint().axes()[name]);
        }
        for property in source.properties() {
            assert_eq!(
                property.answer().receipt.selected,
                initial.friction_sources()[i]
                    .property(property.requirement().name())
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
            card.claims()
                .verify_receipt(&property.answer().receipt)
                .unwrap();
        }
    }
    assert_eq!(moving.state.constitutive_state(), &moving.driven.state);
    // The same force between static/kinetic limits cannot capture motion that
    // the mechanical owner has already admitted as sliding.
    let sliding = law
        .advance_driven_motion_consistent(
            &initial,
            &drive([0.01, 0.0], [10.0, 0.0], EllipticDrivenMode::Sliding),
            &endpoint(),
            solve([0.0, 0.2]),
            driven_controls(),
            &mut || false,
        )
        .unwrap();
    assert!(sliding.actual_distance_m > 0.03);
    assert_eq!(sliding.driven.mode, EllipticDrivenMode::Sliding);
}

#[test]
fn original_temperature_dependent_static_limit_releases_without_switching_claims() {
    let card = directional_card();
    let law = super::super::law(&card, [1.0, 0.0]);
    let initial = law.initial_state();
    let mut input = drive([0.0; 2], [5.0, 0.0], EllipticDrivenMode::Resting);
    input.nominal_pressure_pa = None;
    input.contact_temperature_k = Some(280.0);
    let cold = law
        .advance_driven_motion_consistent(
            &initial,
            &input,
            &query(280.0),
            solve([0.0, 1.0]),
            driven_controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(cold.driven.mode, EllipticDrivenMode::Resting);
    input.contact_temperature_k = Some(320.0);
    let hot = law
        .advance_driven_motion_consistent(
            &cold.state,
            &input,
            &query(320.0),
            solve([0.0, 1.0]),
            driven_controls(),
            &mut || false,
        )
        .unwrap();
    assert_eq!(hot.driven.mode, EllipticDrivenMode::Sliding);
    close(hot.actual_distance_m, 0.015);
    close(hot.driven.velocity_after_m_s[0], 0.3);
    close(hot.driven.dissipated_j, 0.03);
    close(hot.driven.drive_work_j, 0.075);
    assert_eq!(hot.evaluations, 1);
    for i in 0..2 {
        for property in hot.state.friction_sources()[i].properties() {
            assert_eq!(
                property.answer().receipt.selected,
                initial.friction_sources()[i]
                    .property(property.requirement().name())
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
        }
    }
}

#[test]
fn full_matrix_driven_curve_matches_original_native_path_and_energy() {
    let law = motion_law(&motion_card([0.3; 2], [0.6; 2]), [0.6, 0.8]);
    let native = EllipticCoulomb::new([1.2; 2], [0.3, 0.6], [0.6, 0.8]).unwrap();
    let mut input = drive([0.6, 0.8], [0.4, -0.2], EllipticDrivenMode::Sliding);
    input.mobility_per_kg = [[1.0, 0.35], [0.35, 2.0]];
    input.dt_s = 0.05;
    let expected = native
        .advance_driven(
            &native.initial_state(),
            &input,
            driven_controls(),
            &mut || false,
        )
        .unwrap();
    let actual = law
        .advance_driven_motion_consistent(
            &law.initial_state(),
            &input,
            &endpoint(),
            solve([0.0, 0.1]),
            driven_controls(),
            &mut || false,
        )
        .unwrap();
    close(actual.actual_distance_m, expected.distance_m);
    close(actual.driven.drive_work_j, expected.drive_work_j);
    close(actual.driven.dissipated_j, expected.dissipated_j);
    for i in 0..2 {
        close(
            actual.driven.slip_increment_m[i],
            expected.slip_increment_m[i],
        );
        close(
            actual.driven.traction_impulse_ns[i],
            expected.traction_impulse_ns[i],
        );
    }
    assert!(actual.actual_distance_m > norm(actual.driven.slip_increment_m));
    close(
        actual.state.friction_sources()[0].query().axes()["speed"],
        actual.actual_distance_m / input.dt_s,
    );
    close(
        actual.state.friction_sources()[0].query().axes()["path"],
        0.25 + actual.actual_distance_m,
    );
}

#[test]
fn source_driven_failures_preserve_original_history_and_retry() {
    let law = motion_law(&motion_card([0.2, 0.6], [0.2, 0.4]), [1.0, 0.0]);
    let initial = law.initial_state();
    let before = initial.clone();
    for case in 0..5 {
        let mut input = drive([0.0; 2], [13.0, 0.0], EllipticDrivenMode::Resting);
        let mut point = endpoint();
        let mut closure = solve([0.0, 0.2]);
        let mut control = driven_controls();
        match case {
            0 => closure.maximum_evaluations = 1,
            1 => closure.distance_bracket_m = [0.0, 1e-6],
            2 => point = point.with("ambient", 296.0).unwrap(),
            3 => control.maximum_substeps = 0,
            _ => input.velocity_m_s = [0.5, 0.0],
        }
        assert!(
            law.advance_driven_motion_consistent(
                &initial,
                &input,
                &point,
                closure,
                control,
                &mut || false
            )
            .is_err()
        );
        assert_eq!(initial, before);
    }
    let mut polls = 0;
    assert!(
        law.advance_driven_motion_consistent(
            &initial,
            &drive([0.0; 2], [13.0, 0.0], EllipticDrivenMode::Resting),
            &endpoint(),
            solve([0.0, 0.2]),
            driven_controls(),
            &mut || {
                polls += 1;
                polls >= 3
            }
        )
        .is_err()
    );
    assert_eq!(initial, before);
    let retry = law
        .advance_driven_motion_consistent(
            &initial,
            &drive([0.0; 2], [13.0, 0.0], EllipticDrivenMode::Resting),
            &endpoint(),
            solve([0.0, 0.2]),
            driven_controls(),
            &mut || false,
        )
        .unwrap();
    close(retry.actual_distance_m, 0.05);
}
