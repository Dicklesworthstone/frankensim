//! Independent two-sphere collision and complete-advance conservation oracles.

use fs_mbd::sphere_pair::{
    SpherePairAdvance, SpherePairBody, SpherePairError, SpherePairScene, SpherePairTorqueConfig,
    advance_sphere_pair, advance_sphere_pair_with_torques,
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

const MASSES: [f64; 2] = [2.0, 3.0];
const RADII: [f64; 2] = [0.5, 0.25];
const INERTIAS: [f64; 2] = [0.2, 0.075];

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 5e-10,
        "{actual:e} != {expected:e}"
    );
}

fn close_vector(actual: Vec3, expected: Vec3) {
    assert!(
        actual.sub(expected).norm_squared().sqrt() < 5e-10,
        "{actual:?} != {expected:?}"
    );
}

fn scene() -> SpherePairScene {
    SpherePairScene {
        bodies: std::array::from_fn(|i| SpherePairBody {
            properties: MassProperties::new(
                MASSES[i],
                Vec3::ZERO,
                Vec3::new(INERTIAS[i], INERTIAS[i], INERTIAS[i]),
            )
            .unwrap(),
            radius: RADII[i],
        }),
        acceleration_world: Vec3::ZERO,
        restitution: 0.5,
        friction_coefficient: 0.0,
    }
}

fn body(
    index: usize,
    position: Vec3,
    velocity: Vec3,
    spin_world: Vec3,
    orientation: UnitQuaternion,
) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(position, orientation).unwrap(),
        velocity.scale(MASSES[index]),
        orientation.rotate_world_to_body(spin_world.scale(INERTIAS[index])),
    )
    .unwrap()
}

fn velocity(body: RigidBodyState, index: usize) -> Vec3 {
    body.linear_momentum_world().scale(1.0 / MASSES[index])
}

fn spin(body: RigidBodyState, index: usize) -> Vec3 {
    body.pose()
        .orientation()
        .rotate_body_to_world(body.angular_momentum_body())
        .scale(1.0 / INERTIAS[index])
}

fn kinetic(states: [RigidBodyState; 2]) -> f64 {
    (0..2)
        .map(|i| {
            states[i].linear_momentum_world().norm_squared() / (2.0 * MASSES[i])
                + states[i].angular_momentum_body().norm_squared() / (2.0 * INERTIAS[i])
        })
        .sum()
}

fn linear_momentum(states: [RigidBodyState; 2]) -> Vec3 {
    states[0]
        .linear_momentum_world()
        .add(states[1].linear_momentum_world())
}

fn angular_momentum(states: [RigidBodyState; 2]) -> Vec3 {
    states.iter().fold(Vec3::ZERO, |total, body| {
        total
            .add(
                body.pose()
                    .orientation()
                    .rotate_body_to_world(body.angular_momentum_body()),
            )
            .add(
                body.pose()
                    .position_world()
                    .cross(body.linear_momentum_world()),
            )
    })
}

fn balances(
    scene: &SpherePairScene,
    initial: [RigidBodyState; 2],
    duration: f64,
    result: &SpherePairAdvance,
) {
    let after = result.states_after;
    let a = scene.acceleration_world;
    let total_mass = MASSES[0] + MASSES[1];
    close_vector(
        linear_momentum(after).sub(linear_momentum(initial)),
        a.scale(total_mass * duration),
    );
    let force_work: f64 = (0..2)
        .map(|i| {
            a.scale(MASSES[i]).dot(
                after[i]
                    .pose()
                    .position_world()
                    .sub(initial[i].pose().position_world()),
            )
        })
        .sum();
    close(result.external_force_work, force_work);
    close(
        kinetic(after) - kinetic(initial),
        force_work - result.dissipated_energy,
    );
    // Internal contact does not alter the mass-weighted centre trajectory.
    // Integrating its cross product with common acceleration gives the exact
    // applied angular impulse about the fixed world origin.
    let weighted_position = initial[0]
        .pose()
        .position_world()
        .scale(MASSES[0])
        .add(initial[1].pose().position_world().scale(MASSES[1]));
    let applied_angular_impulse = weighted_position
        .scale(duration)
        .add(linear_momentum(initial).scale(0.5 * duration * duration))
        .cross(a);
    close_vector(
        angular_momentum(after).sub(angular_momentum(initial)),
        applied_angular_impulse,
    );
    close(result.energy_balance_residual, 0.0);
    close_vector(result.linear_momentum_residual, Vec3::ZERO);
    close_vector(result.angular_momentum_residual, Vec3::ZERO);
    let separation = after[1]
        .pose()
        .position_world()
        .sub(after[0].pose().position_world());
    close(
        result.gap_after,
        separation.norm_squared().sqrt() - RADII[0] - RADII[1],
    );
    assert!(result.gap_after >= -5e-10);
    assert!(result.dissipated_energy >= 0.0);
}

fn head_on() -> [RigidBodyState; 2] {
    [
        body(
            0,
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
        ),
        body(
            1,
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
        ),
    ]
}

fn spinning_contact() -> [RigidBodyState; 2] {
    [
        body(
            0,
            Vec3::new(-RADII[0], 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(3.0, 1.2, 1.6),
            UnitQuaternion::IDENTITY,
        ),
        body(
            1,
            Vec3::new(RADII[1], 0.0, 0.0),
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(-2.0, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
        ),
    ]
}

#[test]
fn unequal_spheres_exchange_exact_recoil_at_ballistic_contact() {
    let initial = head_on();
    for restitution in [0.0, 0.5, 1.0] {
        let scene = SpherePairScene {
            restitution,
            ..scene()
        };
        let result = advance_sphere_pair(&scene, initial, 0.75, || false).unwrap();
        let event = result.impact.as_ref().unwrap();
        close(event.time_s, 5.0 / 12.0);
        close_vector(event.normal_world, Vec3::new(1.0, 0.0, 0.0));
        close_vector(event.points_world[0], Vec3::new(1.0 / 3.0, 0.0, 0.0));
        close_vector(event.points_world[1], event.points_world[0]);
        let impulse = (18.0 / 5.0) * (1.0 + restitution);
        close(event.normal_impulse, impulse);
        close_vector(event.tangent_impulse_world, Vec3::ZERO);
        close(event.relative_normal_speed_before, -3.0);
        close(event.relative_normal_speed_after, 3.0 * restitution);
        let speeds = [2.0 - impulse / 2.0, -1.0 + impulse / 3.0];
        for i in 0..2 {
            close_vector(
                velocity(result.states_after[i], i),
                Vec3::new(speeds[i], 0.0, 0.0),
            );
            close_vector(spin(result.states_after[i], i), Vec3::ZERO);
        }
        close(
            result.states_after[0].pose().position_world().x,
            -1.0 / 6.0 + speeds[0] / 3.0,
        );
        close(
            result.states_after[1].pose().position_world().x,
            7.0 / 12.0 + speeds[1] / 3.0,
        );
        close(
            result.dissipated_energy,
            (27.0 / 5.0) * (1.0 - restitution * restitution),
        );
        balances(&scene, initial, 0.75, &result);
    }
}

#[test]
fn both_spheres_spin_and_recoil_under_one_two_dimensional_coulomb_impulse() {
    let initial = spinning_contact();
    let slip_direction = Vec3::new(0.0, -0.8, 0.6);
    let tangent_direction = slip_direction.scale(-1.0);
    let normal = Vec3::new(1.0, 0.0, 0.0);
    for (mu, tangent_magnitude) in [(0.2, 12.0 / 35.0), (0.05, 0.18)] {
        let scene = SpherePairScene {
            friction_coefficient: mu,
            ..scene()
        };
        let result = advance_sphere_pair(&scene, initial, 0.1, || false).unwrap();
        let event = result.impact.as_ref().unwrap();
        close(event.time_s, 0.0);
        close(event.normal_impulse, 18.0 / 5.0);
        close_vector(event.points_world[0], Vec3::ZERO);
        close_vector(event.points_world[1], Vec3::ZERO);
        let tangent_impulse = tangent_direction.scale(tangent_magnitude);
        close_vector(event.tangent_impulse_world, tangent_impulse);
        let impulse_on_one = normal.scale(18.0 / 5.0).add(tangent_impulse);
        for i in 0..2 {
            let signed_impulse = impulse_on_one.scale(if i == 0 { -1.0 } else { 1.0 });
            close_vector(
                result.states_after[i]
                    .linear_momentum_world()
                    .sub(initial[i].linear_momentum_world()),
                signed_impulse,
            );
            close_vector(
                spin(result.states_after[i], i),
                spin(initial[i], i)
                    .sub(normal.cross(tangent_impulse).scale(RADII[i] / INERTIAS[i])),
            );
        }
        // The contact begins at the world origin, so spin and orbital changes
        // cancel without an external angular impulse, even for unequal radii.
        close_vector(
            angular_momentum(result.states_after),
            angular_momentum(initial),
        );
        let contact_velocity_zero = velocity(result.states_after[0], 0)
            .add(spin(result.states_after[0], 0).cross(normal.scale(RADII[0])));
        let contact_velocity_one = velocity(result.states_after[1], 1)
            .add(spin(result.states_after[1], 1).cross(normal.scale(-RADII[1])));
        let relative = contact_velocity_one.sub(contact_velocity_zero);
        close(relative.x, 1.0);
        close_vector(
            relative.sub(normal.scale(relative.x)),
            slip_direction.scale(1.0 - (35.0 / 12.0) * tangent_magnitude),
        );
        close(
            result.dissipated_energy,
            1.8 + tangent_magnitude - (35.0 / 24.0) * tangent_magnitude * tangent_magnitude,
        );
        balances(&scene, initial, 0.1, &result);
    }
}

#[test]
fn continuous_detection_catches_crossing_and_preserves_noncolliding_flight() {
    let scene = SpherePairScene {
        restitution: 1.0,
        ..scene()
    };
    let initial = [
        body(
            0,
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(4.0, 0.0, 0.0),
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
        ),
        body(
            1,
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(-4.0, 0.0, 0.0),
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
        ),
    ];
    // Without CCD both step endpoints would be separated after passing through.
    let result = advance_sphere_pair(&scene, initial, 0.5, || false).unwrap();
    close(result.impact.as_ref().unwrap().time_s, 5.0 / 32.0);
    close(result.impact.as_ref().unwrap().normal_impulse, 96.0 / 5.0);
    close(result.states_after[0].pose().position_world().x, -2.3);
    close(result.states_after[1].pose().position_world().x, 1.2);
    close(result.dissipated_energy, 0.0);
    balances(&scene, initial, 0.5, &result);

    for separating in [false, true] {
        let miss = [
            body(
                0,
                Vec3::new(-1.0, 0.0, 0.0),
                Vec3::new(if separating { -4.0 } else { 4.0 }, 0.0, 0.0),
                Vec3::ZERO,
                UnitQuaternion::IDENTITY,
            ),
            body(
                1,
                Vec3::new(1.0, if separating { 0.0 } else { 0.8 }, 0.0),
                Vec3::new(if separating { 4.0 } else { -4.0 }, 0.0, 0.0),
                Vec3::ZERO,
                UnitQuaternion::IDENTITY,
            ),
        ];
        let result = advance_sphere_pair(&scene, miss, 0.5, || false).unwrap();
        assert!(result.impact.is_none());
        for i in 0..2 {
            close_vector(
                result.states_after[i].pose().position_world(),
                miss[i]
                    .pose()
                    .position_world()
                    .add(velocity(miss[i], i).scale(0.5)),
            );
            close_vector(velocity(result.states_after[i], i), velocity(miss[i], i));
        }
        balances(&scene, miss, 0.5, &result);
    }
}

#[test]
fn common_acceleration_and_world_frame_changes_preserve_internal_collision() {
    let baseline_scene = SpherePairScene {
        friction_coefficient: 0.2,
        ..scene()
    };
    let initial = spinning_contact();
    let duration = 0.1;
    let baseline = advance_sphere_pair(&baseline_scene, initial, duration, || false).unwrap();
    let rotation = UnitQuaternion::new(0.5, 0.5, 0.5, 0.5).unwrap();
    let translation = Vec3::new(3.0, -4.0, 2.0);
    let boost = Vec3::new(0.3, -0.2, 0.1);
    let acceleration = Vec3::new(0.3, -2.0, 0.4);
    let transformed_scene = SpherePairScene {
        acceleration_world: acceleration,
        ..baseline_scene
    };
    let transformed_initial = std::array::from_fn(|i| {
        body(
            i,
            rotation
                .rotate_body_to_world(initial[i].pose().position_world())
                .add(translation),
            rotation
                .rotate_body_to_world(velocity(initial[i], i))
                .add(boost),
            rotation.rotate_body_to_world(spin(initial[i], i)),
            rotation,
        )
    });
    let transformed =
        advance_sphere_pair(&transformed_scene, transformed_initial, duration, || false).unwrap();
    for i in 0..2 {
        close_vector(
            transformed.states_after[i].pose().position_world(),
            rotation
                .rotate_body_to_world(baseline.states_after[i].pose().position_world())
                .add(translation)
                .add(boost.scale(duration))
                .add(acceleration.scale(0.5 * duration * duration)),
        );
        close_vector(
            velocity(transformed.states_after[i], i),
            rotation
                .rotate_body_to_world(velocity(baseline.states_after[i], i))
                .add(boost)
                .add(acceleration.scale(duration)),
        );
        close_vector(
            spin(transformed.states_after[i], i),
            rotation.rotate_body_to_world(spin(baseline.states_after[i], i)),
        );
        for basis in [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ] {
            close_vector(
                transformed.states_after[i]
                    .pose()
                    .orientation()
                    .rotate_body_to_world(basis),
                rotation.rotate_body_to_world(
                    baseline.states_after[i]
                        .pose()
                        .orientation()
                        .rotate_body_to_world(basis),
                ),
            );
        }
    }
    close(transformed.dissipated_energy, baseline.dissipated_energy);
    close(
        transformed.impact.as_ref().unwrap().normal_impulse,
        baseline.impact.as_ref().unwrap().normal_impulse,
    );
    close_vector(
        transformed.impact.as_ref().unwrap().tangent_impulse_world,
        rotation.rotate_body_to_world(baseline.impact.as_ref().unwrap().tangent_impulse_world),
    );
    balances(
        &transformed_scene,
        transformed_initial,
        duration,
        &transformed,
    );
}

#[test]
fn invalid_overlap_and_cancelled_advances_do_not_publish_partial_states() {
    let scene = scene();
    let initial = head_on();
    let before = initial;
    assert!(matches!(
        advance_sphere_pair(&scene, initial, 0.75, || true),
        Err(SpherePairError::Cancelled)
    ));
    assert!(matches!(
        advance_sphere_pair(&scene, initial, 0.0, || false),
        Err(SpherePairError::Invalid(_))
    ));
    let overlap = [
        body(
            0,
            Vec3::new(-0.5, 0.0, 0.0),
            Vec3::ZERO,
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
        ),
        body(
            1,
            Vec3::new(0.1, 0.0, 0.0),
            Vec3::ZERO,
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
        ),
    ];
    assert!(matches!(
        advance_sphere_pair(&scene, overlap, 0.1, || false),
        Err(SpherePairError::Penetration { .. })
    ));
    let mut anisotropic = scene.clone();
    anisotropic.bodies[0].properties =
        MassProperties::new(MASSES[0], Vec3::ZERO, Vec3::new(0.2, 0.21, 0.22)).unwrap();
    assert!(matches!(
        advance_sphere_pair(&anisotropic, initial, 0.75, || false),
        Err(SpherePairError::Invalid(_))
    ));
    assert_eq!(initial, before);
    let accepted = advance_sphere_pair(&scene, initial, 0.75, || false).unwrap();
    close(accepted.states_after[0].pose().position_world().x, -0.4);
    close(accepted.states_after[1].pose().position_world().x, 0.85);
}


// G1/G3: integrate each prescribed world couple independently on the two
// analytic flight intervals. The original contact impulse changes the initial
// spin for the second interval; it is not work supplied by the torque source.
fn torque_balances(
    scene: &SpherePairScene,
    initial: [RigidBodyState; 2],
    torques: [Vec3; 2],
    duration: f64,
    result: &SpherePairAdvance,
) {
    let before_hit = result.impact.map_or(duration, |event| event.time_s);
    let mut expected_work = 0.0;
    for i in 0..2 {
        let omega = spin(initial[i], i);
        let alpha = torques[i].scale(1.0 / INERTIAS[i]);
        expected_work += torques[i].dot(omega.add(alpha.scale(0.5 * before_hit))) * before_hit;
        if let Some(event) = result.impact {
            let arm = event
                .normal_world
                .scale(if i == 0 { RADII[i] } else { -RADII[i] });
            let impulse = event
                .normal_world
                .scale(event.normal_impulse)
                .add(event.tangent_impulse_world)
                .scale(if i == 0 { -1.0 } else { 1.0 });
            let omega_after = omega
                .add(alpha.scale(before_hit))
                .add(arm.cross(impulse).scale(1.0 / INERTIAS[i]));
            let rest = duration - before_hit;
            expected_work += torques[i].dot(omega_after.add(alpha.scale(0.5 * rest))) * rest;
        }
    }
    close(result.external_torque_work, expected_work);
    let expected_torque_impulse = torques[0].add(torques[1]).scale(duration);
    close_vector(result.external_torque_impulse_world, expected_torque_impulse);
    let force_work: f64 = (0..2)
        .map(|i| {
            scene.acceleration_world.scale(MASSES[i]).dot(
                result.states_after[i]
                    .pose()
                    .position_world()
                    .sub(initial[i].pose().position_world()),
            )
        })
        .sum();
    close(result.external_force_work, force_work);
    close(
        kinetic(result.states_after) - kinetic(initial),
        force_work + expected_work - result.dissipated_energy,
    );
    let weighted_position = initial[0]
        .pose()
        .position_world()
        .scale(MASSES[0])
        .add(initial[1].pose().position_world().scale(MASSES[1]));
    let orbital_impulse = weighted_position
        .scale(duration)
        .add(linear_momentum(initial).scale(0.5 * duration * duration))
        .cross(scene.acceleration_world);
    close_vector(
        angular_momentum(result.states_after).sub(angular_momentum(initial)),
        orbital_impulse.add(expected_torque_impulse),
    );
    close_vector(
        linear_momentum(result.states_after).sub(linear_momentum(initial)),
        scene.acceleration_world.scale((MASSES[0] + MASSES[1]) * duration),
    );
    close(result.energy_balance_residual, 0.0);
    close_vector(result.linear_momentum_residual, Vec3::ZERO);
    close_vector(result.angular_momentum_residual, Vec3::ZERO);
}

#[test]
fn g1_held_world_torque_spinup_braking_reversal_and_zero_load_parity() {
    let scene = scene();
    let orientation = UnitQuaternion::new(0.7, 0.2, -0.3, 0.4).unwrap();
    for (omega, torque, duration) in [(0.0, 0.4, 0.5), (3.0, -0.4, 0.5), (1.0, -0.4, 1.0)] {
        let initial = [
            body(
                0,
                Vec3::new(-2.0, 0.0, 0.0),
                Vec3::ZERO,
                Vec3::new(0.0, 0.0, omega),
                orientation,
            ),
            body(
                1,
                Vec3::new(2.0, 0.0, 0.0),
                Vec3::ZERO,
                Vec3::ZERO,
                UnitQuaternion::IDENTITY,
            ),
        ];
        let torques = [Vec3::new(0.0, 0.0, torque), Vec3::ZERO];
        let result = advance_sphere_pair_with_torques(
            &scene,
            initial,
            torques,
            duration,
            SpherePairTorqueConfig::default(),
            || false,
        )
        .unwrap();
        assert!(result.impact.is_none());
        close_vector(
            spin(result.states_after[0], 0),
            Vec3::new(0.0, 0.0, omega + torque * duration / INERTIAS[0]),
        );
        let angle = omega * duration + 0.5 * torque * duration * duration / INERTIAS[0];
        let expected = orientation
            .right_exp(orientation.rotate_world_to_body(Vec3::new(0.0, 0.0, angle)))
            .unwrap();
        for basis in [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ] {
            close_vector(
                result.states_after[0]
                    .pose()
                    .orientation()
                    .rotate_body_to_world(basis),
                expected.rotate_body_to_world(basis),
            );
        }
        if omega == 3.0 {
            assert!(result.external_torque_work < 0.0);
        }
        torque_balances(&scene, initial, torques, duration, &result);
        let plain = advance_sphere_pair(&scene, initial, duration, || false).unwrap();
        let zero = advance_sphere_pair_with_torques(
            &scene,
            initial,
            [Vec3::ZERO; 2],
            duration,
            SpherePairTorqueConfig::default(),
            || false,
        )
        .unwrap();
        assert_eq!(plain, zero);
        assert_eq!(zero.attitude_substeps, 0);
    }
}

#[test]
fn g1_preimpact_torque_reaches_friction_and_retains_work_and_frame_covariance() {
    let scene = SpherePairScene {
        acceleration_world: Vec3::new(0.3, -2.0, 0.4),
        friction_coefficient: 0.2,
        ..scene()
    };
    let initial = head_on();
    let torques = [
        Vec3::new(0.0, 0.0, 0.24),
        Vec3::new(0.0, 0.0, -0.06),
    ];
    let duration = 0.75;
    let result = advance_sphere_pair_with_torques(
        &scene,
        initial,
        torques,
        duration,
        SpherePairTorqueConfig::default(),
        || false,
    )
    .unwrap();
    let event = result.impact.unwrap();
    close(event.time_s, 5.0 / 12.0);
    close(event.normal_impulse, 27.0 / 5.0);
    // Independent scalar slip: -(r0*alpha0 + r1*alpha1)*t = -1/6.
    // Native two-body tangent mobility is 35/12; sticking needs Jy = 2/35.
    close_vector(event.tangent_impulse_world, Vec3::new(0.0, 2.0 / 35.0, 0.0));
    close(result.dissipated_energy, 81.0 / 20.0 + 1.0 / 210.0);
    for i in 0..2 {
        close_vector(
            spin(result.states_after[i], i),
            Vec3::new(
                0.0,
                0.0,
                (torques[i].z * duration - RADII[i] * 2.0 / 35.0) / INERTIAS[i],
            ),
        );
    }
    torque_balances(&scene, initial, torques, duration, &result);
    let unforced = advance_sphere_pair(&scene, initial, duration, || false).unwrap();
    close_vector(unforced.impact.unwrap().tangent_impulse_world, Vec3::ZERO);

    let rotation = UnitQuaternion::new(0.5, 0.5, 0.5, 0.5).unwrap();
    let shift = Vec3::new(3.0, -4.0, 2.0);
    let moved_scene = SpherePairScene {
        acceleration_world: rotation.rotate_body_to_world(scene.acceleration_world),
        ..scene
    };
    let moved_initial = std::array::from_fn(|i| {
        body(
            i,
            rotation
                .rotate_body_to_world(initial[i].pose().position_world())
                .add(shift),
            rotation.rotate_body_to_world(velocity(initial[i], i)),
            rotation.rotate_body_to_world(spin(initial[i], i)),
            rotation,
        )
    });
    let moved_torques = torques.map(|torque| rotation.rotate_body_to_world(torque));
    let moved = advance_sphere_pair_with_torques(
        &moved_scene,
        moved_initial,
        moved_torques,
        duration,
        SpherePairTorqueConfig::default(),
        || false,
    )
    .unwrap();
    for i in 0..2 {
        close_vector(
            moved.states_after[i].pose().position_world(),
            rotation
                .rotate_body_to_world(result.states_after[i].pose().position_world())
                .add(shift),
        );
        close_vector(
            velocity(moved.states_after[i], i),
            rotation.rotate_body_to_world(velocity(result.states_after[i], i)),
        );
        close_vector(
            spin(moved.states_after[i], i),
            rotation.rotate_body_to_world(spin(result.states_after[i], i)),
        );
    }
    close(moved.external_torque_work, result.external_torque_work);
    close(moved.dissipated_energy, result.dissipated_energy);
    torque_balances(&moved_scene, moved_initial, moved_torques, duration, &moved);
}

// Independent fourth-order integration of a transported marker, x' = omega(t)
// cross x. This is not the production quaternion-exponential discretization.
fn marker_reference(omega: Vec3, alpha: Vec3) -> Vec3 {
    let mut marker = Vec3::new(1.0, 0.0, 0.0);
    let h = 1.0 / 2048.0;
    for i in 0..2048 {
        let t = i as f64 * h;
        let f = |time: f64, x: Vec3| omega.add(alpha.scale(time)).cross(x);
        let k1 = f(t, marker);
        let k2 = f(t + 0.5 * h, marker.add(k1.scale(0.5 * h)));
        let k3 = f(t + 0.5 * h, marker.add(k2.scale(0.5 * h)));
        let k4 = f(t + h, marker.add(k3.scale(h)));
        marker = marker.add(
            k1.add(k2.scale(2.0))
                .add(k3.scale(2.0))
                .add(k4)
                .scale(h / 6.0),
        );
    }
    marker
}

#[test]
fn g3_changing_spin_axis_has_second_order_attitude_and_exact_angular_impulse() {
    let scene = scene();
    let omega = Vec3::new(1.0, -0.6, 0.3);
    let alpha = Vec3::new(0.8, 0.5, -0.6);
    let initial = [
        body(
            0,
            Vec3::new(-2.0, 0.0, 0.0),
            Vec3::ZERO,
            omega,
            UnitQuaternion::IDENTITY,
        ),
        body(
            1,
            Vec3::new(2.0, 0.0, 0.0),
            Vec3::ZERO,
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
        ),
    ];
    let torques = [alpha.scale(INERTIAS[0]), Vec3::ZERO];
    let reference = marker_reference(omega, alpha);
    let mut errors = Vec::new();
    for steps in [16usize, 32, 64] {
        let result = advance_sphere_pair_with_torques(
            &scene,
            initial,
            torques,
            1.0,
            SpherePairTorqueConfig {
                max_attitude_step_s: 1.0 / steps as f64,
                max_attitude_substeps: steps,
            },
            || false,
        )
        .unwrap();
        assert_eq!(result.attitude_substeps, steps);
        close_vector(spin(result.states_after[0], 0), omega.add(alpha));
        let actual = result.states_after[0]
            .pose()
            .orientation()
            .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
        errors.push(actual.sub(reference).norm_squared().sqrt());
        torque_balances(&scene, initial, torques, 1.0, &result);
    }
    for pair in errors.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!(
            (3.8..4.2).contains(&ratio),
            "attitude refinement ratio {ratio:?}: {errors:?}"
        );
    }
}

#[test]
fn g4_torque_attitude_budget_and_late_cancellation_publish_no_partial_pair() {
    let scene = scene();
    let omega = Vec3::new(1.0, -0.6, 0.3);
    let alpha = Vec3::new(0.8, 0.5, -0.6);
    let initial = [
        body(
            0,
            Vec3::new(-2.0, 0.0, 0.0),
            Vec3::ZERO,
            omega,
            UnitQuaternion::IDENTITY,
        ),
        body(
            1,
            Vec3::new(2.0, 0.0, 0.0),
            Vec3::ZERO,
            omega,
            UnitQuaternion::IDENTITY,
        ),
    ];
    let torques = INERTIAS.map(|inertia| alpha.scale(inertia));
    let config = SpherePairTorqueConfig {
        max_attitude_step_s: 1.0 / 64.0,
        max_attitude_substeps: 128,
    };
    let before = initial;
    assert!(matches!(
        advance_sphere_pair_with_torques(
            &scene,
            initial,
            torques,
            1.0,
            SpherePairTorqueConfig {
                max_attitude_substeps: 127,
                ..config
            },
            || false,
        ),
        Err(SpherePairError::AttitudeBudget)
    ));
    let mut polls = 0;
    assert!(matches!(
        advance_sphere_pair_with_torques(&scene, initial, torques, 1.0, config, || {
            polls += 1;
            polls == 80
        }),
        Err(SpherePairError::Cancelled)
    ));
    assert!(matches!(
        advance_sphere_pair_with_torques(
            &scene,
            initial,
            [Vec3::new(f64::NAN, 0.0, 0.0), Vec3::ZERO],
            1.0,
            config,
            || false,
        ),
        Err(SpherePairError::Invalid(_))
    ));
    assert!(matches!(
        advance_sphere_pair_with_torques(
            &scene,
            initial,
            torques,
            1.0,
            SpherePairTorqueConfig {
                max_attitude_step_s: 0.0,
                ..config
            },
            || false,
        ),
        Err(SpherePairError::Invalid(_))
    ));
    assert_eq!(initial, before);
    let retry =
        advance_sphere_pair_with_torques(&scene, initial, torques, 1.0, config, || false).unwrap();
    assert_eq!(retry.attitude_substeps, 128);
    assert_eq!(
        retry,
        advance_sphere_pair_with_torques(&scene, initial, torques, 1.0, config, || false).unwrap()
    );
    torque_balances(&scene, initial, torques, 1.0, &retry);
}
