//! G0/G1/G3: configured source work, real geared flight and finite supply.
use fs_airscrew::drive::ShaftDrive;
use fs_flyer::aircraft::{OpenLoopDesign, wright_openloop_v1};
use fs_flyer::propcoupling::CANDIDATE_B;
use fs_flyer::spine::SixDofState;
use fs_flyer::spine::rotors::geared::{GearedAssembly, GearedRotor};
use fs_flyer::spine::rotors::propulsion::drive::{
    DriveInputs, DriveSolveControls, DrivenPropulsion, DrivenState,
};
use fs_flyer::spine::rotors::propulsion::{GearedPropulsion, PropellerHandedness, PropellerMount};
use fs_ga::Mat3;

const ZERO: [f64; 3] = [0.0; 3];
const DT: f64 = 0.002;

fn design(bearing_nms: f64) -> OpenLoopDesign {
    let mut design = wright_openloop_v1()
        .with_shaft_drive(
            ShaftDrive::new(10000.0, 200.0, 4.0, 0.8, bearing_nms).unwrap(),
            0.25,
        )
        .unwrap();
    design.coupling = CANDIDATE_B;
    design
}

fn assembly(mass_kg: f64, output_per_gear: f64) -> GearedAssembly<2> {
    // Explicit synthetic inertias. Changing the generalized coordinate does
    // not change either physical rotor, its absolute spin or the aircraft.
    GearedAssembly::new(
        mass_kg,
        Mat3 {
            m: [1000.0, 0.0, 0.0, 0.0, 400.0, 0.0, 0.0, 0.0, 1200.0],
        },
        [
            GearedRotor {
                axis_body: [1.0, 0.0, 0.0],
                axial_inertia_kgm2: 1.5,
                ratio: output_per_gear,
            },
            GearedRotor {
                axis_body: [1.0, 0.0, 0.0],
                axial_inertia_kgm2: 2.5,
                ratio: -output_per_gear,
            },
        ],
    )
    .unwrap()
}

fn mounts() -> [PropellerMount; 2] {
    [
        PropellerMount {
            rotor_index: 0,
            handedness: PropellerHandedness::Positive,
            thrust_direction_body: [1.0, 0.0, 0.0],
        },
        PropellerMount {
            rotor_index: 1,
            handedness: PropellerHandedness::Negative,
            thrust_direction_body: [1.0, 0.0, 0.0],
        },
    ]
}

fn state(
    assembly: &GearedAssembly<2>,
    output_per_gear: f64,
    remaining_engine_energy_j: Option<f64>,
) -> DrivenState {
    DrivenState {
        motion: assembly
            .state_from_speed(
                SixDofState {
                    pos_m: [0.0, 0.0, -10.0],
                    vel_mps: [12.0 * 0.08_f64.cos(), 0.0, 12.0 * 0.08_f64.sin()],
                    quat: [1.0, 0.0, 0.0, 0.0],
                    omega_body: ZERO,
                },
                0.0,
                35.0 / output_per_gear,
            )
            .unwrap(),
        remaining_engine_energy_j,
    }
}

fn input(command: f64) -> DriveInputs {
    DriveInputs {
        density_kg_m3: 1.225,
        canard_rad: 0.03,
        command,
        gravity_world_mps2: ZERO,
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a} != {b}");
}

#[test]
fn g1_held_drive_work_uses_real_travel_including_startup_and_net_braking() {
    let drive = ShaftDrive::new(1024.0, 64.0, 4.0, 0.75, 2.0).unwrap();
    for (speed, command, angle, engine, output, transmission, bearing) in [
        (0.0, 0.5, 0.125, 4.0, 3.0, 1.0, 0.0),
        (8.0, 0.5, 0.25, 8.0, 2.0, 2.0, 4.0),
        (32.0, 1.0, 0.5, 16.0, -20.0, 4.0, 32.0),
        (8.0, 0.0, 0.25, 0.0, -4.0, 0.0, 4.0),
    ] {
        let work = drive.work_over_angle(speed, command, angle).unwrap();
        assert_eq!(work.engine_work_j, engine);
        assert_eq!(work.prop_work_j, output);
        assert_eq!(work.transmission_heat_j, transmission);
        assert_eq!(work.bearing_heat_j, bearing);
        assert_eq!(work.work_residual_j, 0.0);
    }
    assert_eq!(drive.evaluate(0.0, 0.5).unwrap().engine_power_w, 0.0);
    for angle in [-0.01, f64::NAN, f64::INFINITY] {
        assert_eq!(
            drive.work_over_angle(8.0, 0.5, angle).unwrap_err().code,
            "shaft-drive-invalid"
        );
    }
}

#[test]
fn g3_runtime_commands_reach_both_real_stages_and_debit_only_accepted_work() {
    let design = design(0.5);
    let assembly = assembly(design.gross_mass_kg, 2.0);
    let aero = GearedPropulsion::new(&design, &assembly, mounts()).unwrap();
    let model = DrivenPropulsion::new(&aero, 2.0).unwrap();
    let initial = state(&assembly, 2.0, Some(1000.0));
    let before = initial;
    let mut stages = Vec::new();
    let result = model
        .step(
            &initial,
            2.0,
            DT,
            ZERO,
            None,
            |time, motion| {
                let command = 0.6 + 50.0 * (time - 2.0);
                let point = design
                    .shaft_drive
                    .evaluate(2.0 * assembly.speed(motion)?, command)
                    .unwrap();
                stages.push((time, point));
                Ok(input(command))
            },
            |_, _| Ok(ZERO),
        )
        .unwrap();
    assert_eq!(initial, before);
    assert_eq!(stages.len(), 2);
    assert_eq!([stages[0].0, stages[1].0], [2.0, 2.0 + DT / 2.0]);
    assert_eq!(result.initial_drive, stages[0].1);
    assert_eq!(result.midpoint_drive, stages[1].1);
    assert_ne!(result.initial_drive, result.midpoint_drive);
    assert_eq!(
        result.propulsion.midpoint.loads.motor_nm,
        2.0 * result.midpoint_drive.prop_torque_nm
    );
    let expected_angle = 2.0
        * DT
        * (0.5 * assembly.speed(&initial.motion).unwrap()
            + 0.5 * result.propulsion.motion.speed_rad_s);
    close(result.drive_output_travel_rad, expected_angle, 1e-15);
    close(
        result.work.engine_work_j,
        result.midpoint_drive.engine_torque_nm * 4.0 * expected_angle,
        1e-12,
    );
    close(
        1000.0 - result.state.remaining_engine_energy_j.unwrap(),
        result.work.engine_work_j,
        1e-12,
    );
    close(
        result.work.prop_work_j,
        result.propulsion.motion.motor_work_j,
        1e-12,
    );
    assert!(result.work.transmission_heat_j > 0.0);
    assert!(result.work.bearing_heat_j > 0.0);
    assert!(result.work.work_residual_j.abs() < 1e-12);
    assert!(result.supply_residual_j.abs() < 1e-12);
    assert!(result.drive_work_residual_j.abs() < 1e-12);
    for residual in result.propulsion.motion.angular_momentum_residual {
        assert!(residual.abs() < 1e-9);
    }
    assert_eq!(result.state.motion, result.propulsion.motion.state);
}

#[test]
fn g3_coordinate_rescaling_preserves_physical_flight_source_work_and_losses() {
    let design = design(0.5);
    let mut reference: Option<(SixDofState, f64, f64, f64)> = None;
    for ratio in [1.0, 2.0, 0.5] {
        let assembly = assembly(design.gross_mass_kg, ratio);
        let aero = GearedPropulsion::new(&design, &assembly, mounts()).unwrap();
        let model = DrivenPropulsion::new(&aero, ratio).unwrap();
        let result = model
            .step(
                &state(&assembly, ratio, Some(1000.0)),
                0.0,
                DT,
                ZERO,
                None,
                |_, _| Ok(input(0.8)),
                |_, _| Ok(ZERO),
            )
            .unwrap();
        let values = (
            result.state.motion.rigid,
            ratio * result.propulsion.motion.speed_rad_s,
            result.work.engine_work_j,
            result.work.bearing_heat_j,
        );
        if let Some(reference) = reference {
            for (actual, expected) in values.0.vel_mps.into_iter().zip(reference.0.vel_mps) {
                close(actual, expected, 1e-12);
            }
            for (actual, expected) in values.0.omega_body.into_iter().zip(reference.0.omega_body) {
                close(actual, expected, 1e-12);
            }
            close(values.1, reference.1, 1e-12);
            close(values.2, reference.2, 1e-12);
            close(values.3, reference.3, 1e-12);
        } else {
            reference = Some(values);
        }
    }
}

#[test]
fn g0_finite_supply_refuses_overdraw_then_exact_exhaustion_coasts_with_bearing_heat() {
    let design = design(0.5);
    let assembly = assembly(design.gross_mass_kg, 1.0);
    let aero = GearedPropulsion::new(&design, &assembly, mounts()).unwrap();
    let model = DrivenPropulsion::new(&aero, 1.0).unwrap();
    let external = state(&assembly, 1.0, None);
    let run = |state: &DrivenState| {
        model.step(
            state,
            0.0,
            DT,
            ZERO,
            None,
            |_, _| Ok(input(0.8)),
            |_, _| Ok(ZERO),
        )
    };
    let supplied = run(&external).unwrap();
    let required = supplied.work.engine_work_j;
    assert!(required > 0.0);
    let insufficient = DrivenState {
        remaining_engine_energy_j: Some(required * 0.5),
        ..external
    };
    let before = insufficient;
    assert_eq!(
        run(&insufficient).unwrap_err().code,
        "shaft-energy-exhausted"
    );
    assert_eq!(insufficient, before);
    let exact = DrivenState {
        remaining_engine_energy_j: Some(required),
        ..external
    };
    let exhausted = run(&exact).unwrap();
    assert_eq!(exhausted.state.remaining_engine_energy_j, Some(0.0));
    assert_eq!(exhausted.state.motion, supplied.state.motion);
    let coast = run(&exhausted.state).unwrap();
    assert_eq!(coast.initial_drive.engine_torque_nm, 0.0);
    assert_eq!(coast.midpoint_drive.engine_torque_nm, 0.0);
    assert_eq!(coast.work.engine_work_j, 0.0);
    assert_eq!(coast.state.remaining_engine_energy_j, Some(0.0));
    assert!(coast.work.bearing_heat_j > 0.0);
    assert!(coast.propulsion.motion.motor_work_j < 0.0);
    assert!(
        coast
            .propulsion
            .midpoint
            .coupled
            .torque_nm
            .iter()
            .all(|q| q.abs() > 0.0)
    );
}

#[test]
fn g3_net_braking_still_debits_positive_engine_work_and_command_zero_preserves_supply() {
    let design = design(10.0);
    let assembly = assembly(design.gross_mass_kg, 1.0);
    let aero = GearedPropulsion::new(&design, &assembly, mounts()).unwrap();
    let model = DrivenPropulsion::new(&aero, 1.0).unwrap();
    let initial = state(&assembly, 1.0, Some(1000.0));
    for command in [1.0, 0.0] {
        let result = model
            .step(
                &initial,
                0.0,
                DT,
                ZERO,
                None,
                |_, _| Ok(input(command)),
                |_, _| Ok(ZERO),
            )
            .unwrap();
        assert!(result.propulsion.motion.motor_work_j < 0.0);
        assert!(result.work.bearing_heat_j > 0.0);
        if command == 0.0 {
            assert_eq!(result.work.engine_work_j, 0.0);
            assert_eq!(result.state.remaining_engine_energy_j, Some(1000.0));
        } else {
            assert!(result.work.engine_work_j > 0.0);
            assert!(result.state.remaining_engine_energy_j.unwrap() < 1000.0);
        }
        close(
            result.work.engine_work_j,
            result.propulsion.motion.motor_work_j
                + result.work.transmission_heat_j
                + result.work.bearing_heat_j,
            1e-11,
        );
    }
}

#[test]
fn g0_invalid_supply_or_late_command_refuses_without_publishing_partial_state() {
    let design = design(0.5);
    let assembly = assembly(design.gross_mass_kg, 1.0);
    let aero = GearedPropulsion::new(&design, &assembly, mounts()).unwrap();
    let model = DrivenPropulsion::new(&aero, 1.0).unwrap();
    for ratio in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(DrivenPropulsion::new(&aero, ratio).is_err());
    }
    for remaining in [-1.0, f64::NAN, f64::INFINITY] {
        let mut calls = 0;
        let result = model.step(
            &state(&assembly, 1.0, Some(remaining)),
            0.0,
            DT,
            ZERO,
            None,
            |_, _| {
                calls += 1;
                Ok(input(1.0))
            },
            |_, _| Ok(ZERO),
        );
        assert_eq!(result.unwrap_err().code, "driven-propulsion-invalid");
        assert_eq!(calls, 0);
    }
    for remaining in [Some(1000.0), Some(0.0)] {
        let initial = state(&assembly, 1.0, remaining);
        let before = initial;
        let mut calls = 0;
        let refused = model.step(
            &initial,
            0.0,
            DT,
            ZERO,
            None,
            |time, _| {
                calls += 1;
                Ok(input(if time == 0.0 { 0.8 } else { f64::NAN }))
            },
            |_, _| Ok(ZERO),
        );
        assert_eq!(refused.unwrap_err().code, "shaft-drive-invalid");
        assert_eq!(calls, 2);
        assert_eq!(initial, before);
    }
}

fn source_controls() -> DriveSolveControls {
    DriveSolveControls {
        maximum_iterations: 64,
        speed_tolerance_rad_s: 1e-11,
        energy_tolerance_j: 1e-9,
    }
}

fn braking_design(knee_rad_s: f64) -> OpenLoopDesign {
    let mut design = design(20.0)
        .with_shaft_drive(
            ShaftDrive::new(100.0, 4.0 * knee_rad_s, 4.0, 0.8, 20.0).unwrap(),
            1.0,
        )
        .unwrap();
    // Keep the aerodynamic residual below the independent shaft-source check.
    design.coupling.tol = 1e-10;
    design.coupling.cap = 64;
    design
}

#[test]
fn g1_accepted_source_work_respects_power_while_braking_and_changes_real_motion() {
    let design = braking_design(5.0);
    let assembly = assembly(design.gross_mass_kg, 1.0);
    let aero = GearedPropulsion::new(&design, &assembly, mounts()).unwrap();
    let model = DrivenPropulsion::new(&aero, 1.0).unwrap();
    let initial = state(&assembly, 1.0, Some(3.0));
    let before = initial;
    let dt = 0.02;
    let command = |time| 1.0 - 10.0 * time;
    let empirical = model
        .step(
            &initial,
            0.0,
            dt,
            ZERO,
            None,
            |time, _| Ok(input(command(time))),
            |_, _| Ok(ZERO),
        )
        .unwrap();
    let rated_work = command(0.5 * dt) * 100.0 * dt;
    // Explicit midpoint samples a slower predictor, then performs more angle:
    // the old source's power cap therefore does not bound its accepted work.
    assert!(empirical.work.engine_work_j > rated_work + 1e-7);
    let mut calls = Vec::new();
    let solved = model
        .step_power_limited(
            &initial,
            0.0,
            dt,
            ZERO,
            None,
            source_controls(),
            |time, _| {
                calls.push(time);
                Ok(input(command(time)))
            },
            |_, _| Ok(ZERO),
        )
        .unwrap();
    assert_eq!(calls, [0.0, dt / 2.0]);
    assert!(solved.source_iterations > 1 && solved.source_iterations <= 64);
    assert!(solved.mean_speed_residual_rad_s.abs() <= source_controls().speed_tolerance_rad_s);
    let result = solved.driven;
    assert_eq!(initial, before);
    assert_eq!(
        result.propulsion.initial.coupled,
        empirical.propulsion.initial.coupled
    );
    assert_eq!(
        result.propulsion.midpoint.coupled,
        empirical.propulsion.midpoint.coupled
    );
    assert_ne!(
        result.state.motion.momentum,
        empirical.state.motion.momentum
    );
    close(result.work.engine_work_j, rated_work, 1e-9);
    close(
        result.state.remaining_engine_energy_j.unwrap(),
        3.0 - rated_work,
        1e-9,
    );
    let actual_mean =
        0.5 * assembly.speed(&initial.motion).unwrap() + 0.5 * result.propulsion.motion.speed_rad_s;
    close(
        result.midpoint_drive.engine_speed_rad_s / 4.0,
        actual_mean,
        1e-11,
    );
    close(result.work.transmission_heat_j, 0.2 * rated_work, 1e-9);
    close(
        result.work.bearing_heat_j,
        20.0 * actual_mean * actual_mean * dt,
        1e-9,
    );
    close(
        result.work.prop_work_j,
        result.propulsion.motion.motor_work_j,
        1e-10,
    );
    assert!(result.work.prop_work_j < 0.0);
    let loads = result.propulsion.midpoint.loads;
    close(
        result.state.motion.momentum,
        initial.motion.momentum
            + dt * (loads.motor_nm + loads.aerodynamic_nm[0] - loads.aerodynamic_nm[1]),
        1e-12,
    );
}

#[test]
fn g3_power_solve_preserves_gear_coordinates_and_crosses_the_rating_knee() {
    for knee in [5.0, 34.0] {
        let design = braking_design(knee);
        let mut reference = None;
        for ratio in [1.0, 2.0] {
            let assembly = assembly(design.gross_mass_kg, ratio);
            let aero = GearedPropulsion::new(&design, &assembly, mounts()).unwrap();
            let model = DrivenPropulsion::new(&aero, ratio).unwrap();
            let initial = state(&assembly, ratio, Some(10.0));
            let dt = 0.02;
            let result = model
                .step_power_limited(
                    &initial,
                    0.0,
                    dt,
                    ZERO,
                    None,
                    source_controls(),
                    |_, _| Ok(input(1.0)),
                    |_, _| Ok(ZERO),
                )
                .unwrap()
                .driven;
            let mean = ratio
                * (0.5 * assembly.speed(&initial.motion).unwrap()
                    + 0.5 * result.propulsion.motion.speed_rad_s);
            if knee == 34.0 {
                assert!(35.0 > knee && mean < knee);
                close(result.work.engine_work_j, (100.0 / knee) * mean * dt, 1e-9);
                assert!(result.work.engine_work_j < 100.0 * dt);
            } else {
                assert!(mean > knee);
                close(result.work.engine_work_j, 100.0 * dt, 1e-9);
            }
            let values = (
                result.state.motion.rigid,
                mean,
                result.work.engine_work_j,
                result.work.bearing_heat_j,
            );
            if let Some((rigid, rate, engine, heat)) = reference {
                assert_eq!(values.0, rigid);
                close(values.1, rate, 1e-11);
                close(values.2, engine, 1e-9);
                close(values.3, heat, 1e-9);
            } else {
                reference = Some(values);
            }
        }
    }
}

#[test]
fn g0_power_solve_keeps_empty_coasting_and_refuses_trial_budget_or_energy_overdraw() {
    let design = braking_design(5.0);
    let assembly = assembly(design.gross_mass_kg, 1.0);
    let aero = GearedPropulsion::new(&design, &assembly, mounts()).unwrap();
    let model = DrivenPropulsion::new(&aero, 1.0).unwrap();
    let dt = 0.02;
    let empty = state(&assembly, 1.0, Some(0.0));
    let coast = model
        .step_power_limited(
            &empty,
            0.0,
            dt,
            ZERO,
            None,
            source_controls(),
            |_, _| Ok(input(1.0)),
            |_, _| Ok(ZERO),
        )
        .unwrap()
        .driven;
    assert_eq!(coast.work.engine_work_j, 0.0);
    assert_eq!(coast.state.remaining_engine_energy_j, Some(0.0));
    assert!(coast.work.bearing_heat_j > 0.0 && coast.work.prop_work_j < 0.0);
    close(
        coast.work.bearing_heat_j + coast.propulsion.motion.motor_work_j,
        0.0,
        1e-10,
    );

    for (remaining, iterations, code) in [
        (3.0, 1, "shaft-drive-not-converged"),
        (0.5, 64, "shaft-energy-exhausted"),
    ] {
        let initial = state(&assembly, 1.0, Some(remaining));
        let before = initial;
        let mut calls = 0;
        let result = model.step_power_limited(
            &initial,
            0.0,
            dt,
            ZERO,
            None,
            DriveSolveControls {
                maximum_iterations: iterations,
                ..source_controls()
            },
            |_, _| {
                calls += 1;
                Ok(input(1.0))
            },
            |_, _| Ok(ZERO),
        );
        assert_eq!(result.unwrap_err().code, code);
        assert_eq!(calls, 2);
        assert_eq!(initial, before);
    }
    let mut calls = 0;
    let invalid = model.step_power_limited(
        &empty,
        0.0,
        dt,
        ZERO,
        None,
        DriveSolveControls {
            maximum_iterations: 65,
            ..source_controls()
        },
        |_, _| {
            calls += 1;
            Ok(input(1.0))
        },
        |_, _| Ok(ZERO),
    );
    assert_eq!(invalid.unwrap_err().code, "driven-propulsion-invalid");
    assert_eq!(calls, 0);
}
