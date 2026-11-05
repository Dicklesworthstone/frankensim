//! Directional friction on automatically discovered native sphere contacts.
use super::*;

fn directional(mu: [f64; 2], direction: Vec3) -> SphereSceneEllipticLaw {
    SphereSceneEllipticLaw {
        coefficients_mu: mu,
        texture_direction: SphereSceneTextureDirection::World(direction),
        restitution: 0.0,
    }
}

#[test]
fn equal_principal_coefficients_recover_scalar_scene_and_native_work() {
    for mu in [0.0, 0.25, 1.0] {
        let initial = [body(
            Vec3::new(0.0, 0.0, RADIUS + CLEARANCE),
            Vec3::new(3.0, 1.0, 0.0),
            1.0,
        )];
        let mut scalar = initial;
        let mut elliptic = initial;
        let gravity = Gravity::new(Vec3::new(0.0, 0.0, -8.0)).unwrap();
        let a = advance_sphere_scene(
            &mut scalar,
            &[property(1.0)],
            &[sphere(10)],
            &[plane(Vec3::ZERO)],
            &[Wrench::ZERO],
            gravity,
            0.0,
            0.5,
            controls(1.0 / 32.0),
            |_| Some(law(mu, 0.0)),
            &CancelGate::new(),
        )
        .unwrap();
        let b = advance_elliptic_sphere_scene(
            &mut elliptic,
            &[property(1.0)],
            &[sphere(10)],
            &[plane(Vec3::ZERO)],
            &[Wrench::ZERO],
            gravity,
            0.0,
            0.5,
            controls(1.0 / 32.0),
            // This is the scalar scene's original normalized tangent frame.
            |_| Some(directional([mu; 2], Vec3::new(0.0, 1.0, 0.0))),
            &CancelGate::new(),
        )
        .unwrap();
        balances(&b, 0.5);
        assert_eq!(scalar, elliptic);
        assert_eq!(a.steps.len(), b.steps.len());
        assert_eq!(a.rejected_steps, b.rejected_steps);
        assert_eq!(a.force_impulse_work_j, b.force_impulse_work_j);
        assert_eq!(a.force_displacement_work_j, b.force_displacement_work_j);
        assert_eq!(a.support_work_j, b.support_work_j);
        assert_eq!(a.unallocated_contact_loss_j, b.unallocated_contact_loss_j);
        for (a, b) in a.steps.iter().zip(&b.steps) {
            assert_eq!(a.contact_keys, b.contact_keys);
            assert_eq!(
                a.native.impact.contact.impulses,
                b.native.impact.contact.impulses
            );
            assert_eq!(
                a.native.impact.contact.relative_velocities_after,
                b.native.impact.contact.relative_velocities_after,
            );
        }
    }
}

#[test]
fn one_free_texture_axis_keeps_sliding_while_the_other_reaches_rolling() {
    let mut bodies = [body(
        Vec3::new(0.0, 0.0, RADIUS + CLEARANCE),
        Vec3::new(3.0, 4.0, 0.0),
        1.0,
    )];
    let report = advance_elliptic_sphere_scene(
        &mut bodies,
        &[property(1.0)],
        &[sphere(10)],
        &[plane(Vec3::ZERO)],
        &[Wrench::ZERO],
        Gravity::new(Vec3::new(0.0, 0.0, -8.0)).unwrap(),
        0.0,
        0.5,
        controls(1.0 / 32.0),
        |_| Some(directional([0.0, 0.5], Vec3::new(1.0, 0.0, 0.0))),
        &CancelGate::new(),
    )
    .unwrap();
    balances(&report, 0.5);
    // For a uniform sphere, tangent mobility is 1/m + R^2/I = 7/2.
    // Stopping the y slip requires total impulse -4/(7/2) = -8/7.
    vector(
        bodies[0].linear_momentum_world(),
        Vec3::new(3.0, 20.0 / 7.0, 0.0),
        2e-8,
    );
    vector(spin_world(bodies[0]), Vec3::new(-4.0 / 7.0, 0.0, 0.0), 2e-8);
    close(bodies[0].pose().position_world().x, 1.5, 2e-9);
    for step in &report.steps {
        let contact = &step.native.impact.contact;
        close(contact.impulses[0][1], 0.0, 1e-12);
        close(contact.relative_velocities_after[0][1], 3.0, 2e-8);
    }
    close(
        report
            .steps
            .last()
            .unwrap()
            .native
            .impact
            .contact
            .relative_velocities_after[0][2],
        0.0,
        2e-8,
    );
}

#[test]
fn rotated_texture_scene_preserves_physical_impulse_and_moving_support_work() {
    let turn = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.7).unwrap();
    let texture = Vec3::new(0.8, 0.6, 0.0);
    let boost = Vec3::new(0.3, -0.2, 0.1);
    let shift = Vec3::new(0.4, -0.7, 0.2);
    let dt = 0.125;
    for rotation in [UnitQuaternion::IDENTITY, turn] {
        let normal = rotation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
        let tangent = rotation.rotate_body_to_world(texture);
        let initial_position = shift.add(normal.scale(RADIUS + CLEARANCE));
        let mut bodies = [RigidBodyState::new(
            Pose::new(initial_position, rotation).unwrap(),
            boost.add(tangent.scale(4.0)),
            Vec3::ZERO,
        )
        .unwrap()];
        let support = SphereScenePlane {
            id: 99,
            point_world: shift,
            normal_world: normal,
            velocity_world: boost,
        };
        let report = advance_elliptic_sphere_scene(
            &mut bodies,
            &[property(1.0)],
            &[sphere(10)],
            &[support],
            &[Wrench::ZERO],
            Gravity::new(normal.scale(-8.0)).unwrap(),
            0.0,
            dt,
            controls(dt),
            // A normal component in an authored direction must not rotate
            // the physical projected ellipse or change its coefficients.
            |_| Some(directional([0.25, 0.5], tangent.add(normal.scale(0.3)))),
            &CancelGate::new(),
        )
        .unwrap();
        balances(&report, dt);
        assert_eq!(report.steps.len(), 1);
        let impulse = normal.sub(tangent.scale(0.25));
        let endpoint_velocity = boost.add(tangent.scale(3.75));
        vector(report.support_impulse_ns, impulse, 2e-8);
        close(report.support_work_j, impulse.dot(boost), 2e-8);
        vector(bodies[0].linear_momentum_world(), endpoint_velocity, 2e-8);
        vector(
            spin_world(bodies[0]),
            normal.cross(tangent).scale(0.125),
            2e-8,
        );
        vector(
            bodies[0].pose().position_world(),
            initial_position.add(endpoint_velocity.scale(dt)),
            2e-8,
        );
        // Native normal-kick loss 1/2 plus tangent loss
        // (1/4) * (4 + 4 - (7/2)*(1/4)) / 2.
        close(report.unallocated_contact_loss_j, 1.390625, 2e-8);
    }
}

#[test]
fn body_texture_rotates_between_intervals_instead_of_staying_world_frozen() {
    let dt = 0.125;
    let spin = core::f64::consts::FRAC_PI_2 / dt;
    let inertia = property(1.0).principal_inertia_body().z;
    let initial = RigidBodyState::new(
        Pose::new(
            Vec3::new(0.0, 0.0, RADIUS + CLEARANCE),
            UnitQuaternion::IDENTITY,
        )
        .unwrap(),
        Vec3::new(4.0, 0.0, 0.0),
        Vec3::new(0.0, 0.0, inertia * spin),
    )
    .unwrap();
    for attached in [false, true] {
        let mut bodies = [initial];
        let report = advance_elliptic_sphere_scene(
            &mut bodies,
            &[property(1.0)],
            &[sphere(10)],
            &[plane(Vec3::ZERO)],
            &[Wrench::ZERO],
            Gravity::new(Vec3::new(0.0, 0.0, -8.0)).unwrap(),
            0.0,
            2.0 * dt,
            controls(dt),
            |_| {
                let mut law = directional([0.0, 0.5], Vec3::new(1.0, 0.0, 0.0));
                if attached {
                    law.texture_direction = SphereSceneTextureDirection::Sphere {
                        id: 10,
                        direction_body: Vec3::new(1.0, 0.0, 0.0),
                    };
                }
                Some(law)
            },
            &CancelGate::new(),
        )
        .unwrap();
        balances(&report, 2.0 * dt);
        assert_eq!(report.steps.len(), 2);
        // The first interval slips freely along body x. After a quarter turn
        // the attached strong axis opposes world x, supplying 0.5 N s.
        close(
            bodies[0].linear_momentum_world().x,
            if attached { 3.5 } else { 4.0 },
            2e-8,
        );
        close(
            spin_world(bodies[0]).y,
            if attached { 0.25 } else { 0.0 },
            2e-8,
        );
        close(spin_world(bodies[0]).z, inertia * spin, 2e-8);
    }
}

#[test]
fn detected_sphere_pair_keeps_finite_mass_recoil_and_attached_texture_identity() {
    let initial = [
        body(
            Vec3::new(-RADIUS - CLEARANCE, 0.0, 0.0),
            Vec3::new(1.0, 4.0, 0.0),
            1.0,
        ),
        body(Vec3::new(RADIUS, 0.0, 0.0), Vec3::new(-1.0, 0.0, 0.0), 1.0),
    ];
    // IDs are deliberately in the opposite order to the supplied body array.
    // Bind texture to sphere 20, which is the detected contact's counterpart.
    let spheres = [sphere(20), sphere(10)];
    let mut bodies = initial;
    let report = advance_elliptic_sphere_scene(
        &mut bodies,
        &[property(1.0); 2],
        &spheres,
        &[],
        &[Wrench::ZERO; 2],
        Gravity::ZERO,
        0.0,
        0.03125,
        controls(0.03125),
        |key| {
            assert_eq!(
                key,
                SphereSceneContactKey::Spheres {
                    first: 10,
                    second: 20
                }
            );
            Some(SphereSceneEllipticLaw {
                coefficients_mu: [0.25, 0.5],
                texture_direction: SphereSceneTextureDirection::Sphere {
                    id: 20,
                    direction_body: Vec3::new(0.0, 1.0, 0.0),
                },
                restitution: 0.0,
            })
        },
        &CancelGate::new(),
    )
    .unwrap();
    balances(&report, 0.03125);
    assert_eq!(report.steps[0].native.impact.contact.islands.len(), 1);
    vector(
        bodies[0].linear_momentum_world(),
        Vec3::new(0.0, 3.75, 0.0),
        2e-8,
    );
    vector(
        bodies[1].linear_momentum_world(),
        Vec3::new(0.0, 0.25, 0.0),
        2e-8,
    );
    for body in bodies {
        vector(spin_world(body), Vec3::new(0.0, 0.0, -0.125), 2e-8);
    }
    vector(report.support_impulse_ns, Vec3::ZERO, 1e-12);
    close(report.support_work_j, 0.0, 1e-12);
}

#[test]
fn invalid_texture_and_late_scene_refusal_leave_original_native_states_unchanged() {
    let initial = [body(
        Vec3::new(0.0, 0.0, RADIUS + CLEARANCE),
        Vec3::new(4.0, 0.0, 0.0),
        1.0,
    )];
    let dt = 0.125;
    for texture in [
        SphereSceneTextureDirection::World(Vec3::ZERO),
        SphereSceneTextureDirection::World(Vec3::new(0.0, 0.0, 1.0)),
        SphereSceneTextureDirection::World(Vec3::new(f64::NAN, 0.0, 0.0)),
        SphereSceneTextureDirection::Sphere {
            id: 999,
            direction_body: Vec3::new(1.0, 0.0, 0.0),
        },
        // This valid source reaches one admitted interval before the step cap.
        SphereSceneTextureDirection::World(Vec3::new(1.0, 0.0, 0.0)),
    ] {
        let mut bodies = initial;
        let mut limited = controls(dt);
        limited.maximum_steps = 1;
        let failure = advance_elliptic_sphere_scene(
            &mut bodies,
            &[property(1.0)],
            &[sphere(10)],
            &[plane(Vec3::ZERO)],
            &[Wrench::ZERO],
            Gravity::new(Vec3::new(0.0, 0.0, -8.0)).unwrap(),
            0.0,
            2.0 * dt,
            limited,
            |_| {
                Some(SphereSceneEllipticLaw {
                    coefficients_mu: [0.25, 0.5],
                    texture_direction: texture,
                    restitution: 0.0,
                })
            },
            &CancelGate::new(),
        );
        assert!(failure.is_err());
        if texture == SphereSceneTextureDirection::World(Vec3::new(1.0, 0.0, 0.0)) {
            assert!(matches!(failure, Err(SphereSceneError::StepBudget)));
        }
        assert_eq!(bodies, initial);
    }
}
