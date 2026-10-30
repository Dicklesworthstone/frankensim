//! G1/G3/G4: external body work, static hold/breakaway and original state owners.
use super::*;
use fs_mbd::Wrench;
use fs_mbd::tangential_friction::directional::driven::*;
use fs_mbd::tangential_friction::directional::trajectory::apply_rigid_directional_trajectory;
use fs_tribo::tangential::elliptic::driven::*;
use fs_tribo::tangential::elliptic::relaxation::EllipticRelaxationControls;

fn driven_controls() -> EllipticDrivenControls {
    EllipticDrivenControls {
        maximum_step_s: 0.002,
        minimum_step_s: 1e-10,
        maximum_substeps: 100_000,
        rest_velocity_tolerance_m_s: 1e-10,
        solve: controls().contact,
        relaxation: EllipticRelaxationControls {
            time_tolerance_s: 1e-11,
            travel_tolerance_m: 1e-11,
            maximum_evaluations: 200_000,
            maximum_root_iterations: 96,
        },
    }
}

fn simple_pair(
    speed: f64,
) -> (
    [RigidBodyState; 3],
    [MassProperties; 3],
    RigidDirectionalContact,
) {
    let a = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(2.0, 2.0, 2.0)).unwrap();
    let b = MassProperties::new(3.0, Vec3::ZERO, Vec3::new(3.0, 3.0, 3.0)).unwrap();
    (
        [
            RigidBodyState::new(
                Pose::identity(),
                Vec3::new(2.0 * speed, 0.0, 0.0),
                Vec3::ZERO,
            )
            .unwrap(),
            RigidBodyState::new(
                Pose::identity(),
                Vec3::new(3.0 * speed, 0.0, 0.0),
                Vec3::ZERO,
            )
            .unwrap(),
            RigidBodyState::new(Pose::identity(), Vec3::new(7.0, 8.0, 9.0), Vec3::ZERO).unwrap(),
        ],
        [a, b, a],
        contact([Vec3::ZERO; 2]),
    )
}

#[test]
fn subthreshold_drive_preserves_stiction_and_common_body_work_then_breaks_away() {
    let law = EllipticCoulomb::new([0.9, 0.8], [0.2, 0.4], [1.0, 0.0]).unwrap();
    let (mut bodies, props, c) = simple_pair(1.0);
    let initial = bodies;
    let mut history = law.initial_state();
    let mut mode = EllipticDrivenMode::Resting;
    let report = apply_rigid_directional_driven(
        &mut bodies,
        &props,
        c,
        [
            Wrench {
                force_world: Vec3::new(2.0, 0.0, 0.0),
                torque_body: Vec3::ZERO,
            },
            Wrench::ZERO,
        ],
        &law,
        &mut history,
        &mut mode,
        0.2,
        controls(),
        driven_controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(mode, EllipticDrivenMode::Resting);
    close(report.driven.distance_m, 0.0);
    close(report.driven.dissipated_j, 0.0);
    close(report.driven.drive_work_j, 0.0);
    close(report.driven.traction_impulse_ns[0], 0.24);
    vector(relative(&bodies, &props, c), Vec3::ZERO);
    close(
        bodies[0].center_of_mass_velocity_world(props[0]).unwrap().x,
        1.08,
    );
    close(
        bodies[1].center_of_mass_velocity_world(props[1]).unwrap().x,
        1.08,
    );
    // The external force acts on a moving body even with zero contact slip.
    // Its work is 2N * (1m/s*.2s + .5*.4m/s²*.2²s²) = .416J.
    close(report.external_work_j, 0.416);
    close(report.kinetic_energy_change_j, 0.416);
    close(report.audit.total_energy_residual_j, 0.0);
    assert_eq!(bodies[2], initial[2]);
    assert_eq!(history.slip_m(), [0.0; 2]);

    let (mut bodies, props, c) = simple_pair(0.0);
    let mut history = law.initial_state();
    let mut mode = EllipticDrivenMode::Resting;
    let report = apply_rigid_directional_driven(
        &mut bodies,
        &props,
        c,
        [
            Wrench {
                force_world: Vec3::new(20.0, 0.0, 0.0),
                torque_body: Vec3::ZERO,
            },
            Wrench::ZERO,
        ],
        &law,
        &mut history,
        &mut mode,
        0.1,
        controls(),
        driven_controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(mode, EllipticDrivenMode::Sliding);
    close(bodies[0].linear_momentum_world().x, 1.8);
    close(bodies[1].linear_momentum_world().x, 0.2);
    close(report.driven.slip_increment_m[0], 1.0 / 24.0);
    close(report.driven.dissipated_j, 1.0 / 12.0);
    close(report.external_work_j, 0.9);
    close(report.kinetic_energy_change_j, 0.9 - 1.0 / 12.0);
    close(report.audit.total_energy_residual_j, 0.0);
}

#[test]
fn full_tensor_driven_response_is_covariant_and_keeps_actual_external_momenta() {
    let law = EllipticCoulomb::new([0.9, 0.8], [0.2, 0.4], [0.6, 0.8]).unwrap();
    let force = [Vec3::new(0.3, -0.1, 0.4), Vec3::new(-0.2, 0.4, -0.1)];
    let torque_world = [Vec3::new(0.1, -0.04, 0.03), Vec3::new(-0.05, 0.1, 0.04)];
    let mut reference: Option<([RigidBodyState; 3], Vec3, f64, f64, f64)> = None;
    for gauged in [false, true] {
        let (mut bodies, props, c) = full_pair(gauged);
        let initial = bodies;
        let external = [0, 1].map(|i| Wrench {
            force_world: force[i],
            torque_body: bodies[i]
                .pose()
                .orientation()
                .rotate_world_to_body(torque_world[i]),
        });
        let mut history = law.initial_state();
        let mut mode = EllipticDrivenMode::Sliding;
        let dt = 0.05;
        let report = apply_rigid_directional_driven(
            &mut bodies,
            &props,
            c,
            external,
            &law,
            &mut history,
            &mut mode,
            dt,
            controls(),
            driven_controls(),
            &CancelGate::new(),
        )
        .unwrap();
        assert!(report.mobility_per_kg[0][1].abs() > 1e-3);
        assert!(report.driven.dissipated_j > 0.0);
        let mut dp_sum = Vec3::ZERO;
        let mut dl_sum = Vec3::ZERO;
        for i in 0..2 {
            let dp = bodies[i]
                .linear_momentum_world()
                .sub(initial[i].linear_momentum_world());
            let dl = bodies[i]
                .angular_momentum_body()
                .sub(initial[i].angular_momentum_body());
            dp_sum = dp_sum.add(dp);
            dl_sum = dl_sum
                .add(initial[i].pose().position_world().cross(dp))
                .add(initial[i].pose().orientation().rotate_body_to_world(dl));
            assert_eq!(bodies[i].pose(), initial[i].pose());
        }
        vector(dp_sum, force[0].add(force[1]).scale(dt));
        vector(
            dl_sum,
            [0, 1].into_iter().fold(Vec3::ZERO, |sum, i| {
                sum.add(
                    initial[i]
                        .pose()
                        .position_world()
                        .cross(force[i])
                        .add(torque_world[i])
                        .scale(dt),
                )
            }),
        );
        close(
            report.kinetic_energy_change_j + report.driven.dissipated_j,
            report.external_work_j,
        );
        close(report.audit.drive_work_residual_j, 0.0);
        close(report.audit.relative_energy_residual_j, 0.0);
        assert_eq!(bodies[2], initial[2]);
        let values = (
            bodies,
            report.force_impulse_on_a_world_ns,
            report.driven.distance_m,
            report.driven.dissipated_j,
            report.external_work_j,
        );
        if let Some((old_bodies, old_impulse, path, heat, work)) = reference {
            for i in 0..2 {
                vector(
                    bodies[i].linear_momentum_world(),
                    old_bodies[i].linear_momentum_world(),
                );
                vector(
                    bodies[i]
                        .pose()
                        .orientation()
                        .rotate_body_to_world(bodies[i].angular_momentum_body()),
                    old_bodies[i]
                        .pose()
                        .orientation()
                        .rotate_body_to_world(old_bodies[i].angular_momentum_body()),
                );
            }
            vector(report.force_impulse_on_a_world_ns, old_impulse);
            close(report.driven.distance_m, path);
            close(report.driven.dissipated_j, heat);
            close(report.external_work_j, work);
        } else {
            reference = Some(values);
        }
    }
}

#[test]
fn zero_drive_uses_original_unforced_curve_and_finite_arrest() {
    let (mut bodies, props, c) = full_pair(false);
    let mut reference = bodies;
    let law = EllipticCoulomb::new([0.9, 0.8], [0.2, 0.4], [0.6, 0.8]).unwrap();
    let mut history = law.initial_state();
    let mut reference_history = history.clone();
    let mut mode = EllipticDrivenMode::Sliding;
    let original = apply_rigid_directional_trajectory(
        &mut reference,
        &props,
        c,
        &law,
        &mut reference_history,
        5.0,
        controls(),
        driven_controls().relaxation,
        &CancelGate::new(),
    )
    .unwrap();
    let driven = apply_rigid_directional_driven(
        &mut bodies,
        &props,
        c,
        [Wrench::ZERO; 2],
        &law,
        &mut history,
        &mut mode,
        5.0,
        controls(),
        driven_controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(bodies, reference);
    assert_eq!(history, reference_history);
    assert_eq!(mode, EllipticDrivenMode::Resting);
    close(driven.driven.distance_m, original.relaxation.distance_m);
    close(driven.driven.dissipated_j, original.relaxation.dissipated_j);
    close(driven.external_work_j, 0.0);
    close(driven.audit.total_energy_residual_j, 0.0);
}

#[test]
fn callback_drive_and_work_audits_preserve_body_history_and_mode_on_failure() {
    let (original, props, c) = full_pair(false);
    let law = EllipticCoulomb::new([0.9, 0.8], [0.2, 0.4], [0.6, 0.8]).unwrap();
    let history = (law.initial_state(), 7_u64);
    let external = [
        Wrench {
            force_world: Vec3::new(0.4, -0.2, 0.1),
            torque_body: Vec3::new(0.1, 0.2, -0.1),
        },
        Wrench::ZERO,
    ];
    for refusal in 0..6 {
        let mut bodies = original;
        let mut state = history.clone();
        let mut mode = EllipticDrivenMode::Sliding;
        let gate = CancelGate::new();
        let result = apply_rigid_directional_driven_with(
            &mut bodies,
            &props,
            c,
            external,
            &mut state,
            &mut mode,
            0.02,
            controls(),
            |old, input| {
                assert_eq!(old, &history);
                assert_eq!(input.mode, EllipticDrivenMode::Sliding);
                if refusal == 1 {
                    return Err("original drive-dependent source unavailable");
                }
                let mut native = law
                    .advance_driven(&old.0, input, driven_controls(), &mut || false)
                    .unwrap();
                let next = (native.state.clone(), old.1 + 1);
                if refusal == 2 {
                    native.drive_work_j += 1.0;
                }
                if refusal == 3 {
                    native.driving_force_n[0] += 1.0;
                }
                if refusal == 4 {
                    gate.request();
                }
                if refusal == 5 {
                    native.mode = EllipticDrivenMode::Resting;
                }
                Ok((native, next))
            },
            &gate,
        );
        if refusal == 0 {
            let accepted = result.unwrap();
            assert_eq!(state.1, 8);
            assert_eq!(state.0, accepted.driven.state);
            assert_eq!(bodies[2], original[2]);
            close(accepted.audit.total_energy_residual_j, 0.0);
        } else {
            assert!(result.is_err());
            assert_eq!(bodies, original);
            assert_eq!(state, history);
            assert_eq!(mode, EllipticDrivenMode::Sliding);
        }
    }
}
