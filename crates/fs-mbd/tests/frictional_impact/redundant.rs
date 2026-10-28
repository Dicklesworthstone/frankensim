//! G1/G3: native redundant support has determined motion, not invented local forces.
use super::*;
use fs_contact::frictional_multicontact::FrictionalContactDeterminacy;

fn tripod() -> [RigidFrictionalContact; 3] {
    [(-1.0, -1.0), (1.0, -1.0), (0.0, 1.0)].map(|(x, y)| {
        let mut c = plane(
            RigidImpactPoint {
                body: 0,
                arm_body: Vec3::new(x, y, -0.25),
            },
            0.1,
            0.0,
        );
        c.id = if y > 0.0 {
            30
        } else if x > 0.0 {
            20
        } else {
            10
        };
        c
    })
}

#[test]
fn g1_tripod_stick_and_slide_determine_native_motion_with_redundant_shear() {
    // A finite-height body has three noncollinear support points. The nine
    // contact rows act on six body velocities; independent shear impulses
    // cannot be inferred from rigid-body mechanics alone.
    let properties = [property(3.0, Vec3::new(1.0, 1.0, 1.0))];
    let contacts = tripod();
    for speed in [0.02, 3.0] {
        let initial = [body(
            Vec3::new(0.0, 0.0, 0.25),
            UnitQuaternion::IDENTITY,
            Vec3::new(3.0 * speed, 0.0, -3.0),
        )];
        let mut bodies = initial;
        let report = apply_rigid_frictional_impacts(
            &mut bodies,
            &properties,
            &contacts,
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
        let FrictionalContactDeterminacy::UniqueVelocity(witness) = &report.contact.determinacy
        else {
            panic!("redundant tripod must report velocity-only determinacy");
        };
        assert_eq!(witness.rank, 6);
        assert_eq!(witness.dimension, 9);
        assert!(witness.friction_feedback_bound < 1.0);

        // Independently: zero normal speed at all three supports fixes vz,
        // wx and wy. Normal impulse totals 3 N s. Its centre of pressure
        // counters the shear moment at height 1/4, leaving zero angular speed.
        // Sticking removes all horizontal momentum; sliding removes mu*3.
        let impulse_x = -(3.0 * speed).min(0.3);
        close_vector(
            bodies[0].linear_momentum_world(),
            Vec3::new(3.0 * speed + impulse_x, 0.0, 0.0),
        );
        close_vector(bodies[0].angular_momentum_body(), Vec3::ZERO);
        let normal = &report.contact.impulses;
        close(normal[2][0], 1.5);
        close(normal[1][0] - normal[0][0], -0.25 * impulse_x);
        close(normal.iter().map(|p| p[0]).sum(), 3.0);
        close(normal.iter().map(|p| p[1]).sum(), impulse_x);
        let expected_loss = 1.5 * (1.0 + speed * speed) - (3.0 * speed + impulse_x).powi(2) / 6.0;
        close(report.contact.dissipated_energy, expected_loss);
        close(report.support_work_j, 0.0);
        close_vector(report.linear_momentum_residual_ns, Vec3::ZERO);
        close_vector(report.angular_momentum_residual_nms, Vec3::ZERO);

        let mut reordered = initial;
        let permuted = [contacts[2], contacts[0], contacts[1]];
        let repeated = apply_rigid_frictional_impacts(
            &mut reordered,
            &properties,
            &permuted,
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(bodies, reordered);
        assert_eq!(report.contact, repeated.contact);
    }
}

#[test]
fn g1_redundant_frictionless_normal_reactions_preserve_restitution() {
    let properties = [property(2.0, Vec3::new(1.0, 1.0, 1.0))];
    let mut bodies = [body(
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::new(0.0, 0.0, -2.0),
    )];
    let first = plane(
        RigidImpactPoint {
            body: 0,
            arm_body: Vec3::ZERO,
        },
        0.0,
        0.5,
    );
    let mut duplicate = first;
    duplicate.id = 2;
    let report = apply_rigid_frictional_impacts(
        &mut bodies,
        &properties,
        &[first, duplicate],
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert!(matches!(
        report.contact.determinacy,
        FrictionalContactDeterminacy::UniqueVelocity(_)
    ));
    close_vector(bodies[0].linear_momentum_world(), Vec3::new(0.0, 0.0, 1.0));
    close(report.contact.impulses.iter().map(|p| p[0]).sum(), 3.0);
    close(report.contact.dissipated_energy, 0.75);
}

#[test]
fn g3_common_support_work_is_determined_but_independent_support_work_refuses() {
    let properties = [property(3.0, Vec3::new(1.0, 1.0, 1.0))];
    let initial = [body(
        Vec3::new(0.0, 0.0, 0.25),
        UnitQuaternion::IDENTITY,
        Vec3::new(0.06, 0.0, -3.0),
    )];
    let mut contacts = tripod();
    for c in &mut contacts {
        c.other = RigidImpactCounterpart::PrescribedPlane {
            point_world: Vec3::ZERO,
            velocity_world: Vec3::new(0.01, 0.0, 0.0),
        };
    }
    let mut bodies = initial;
    let report = apply_rigid_frictional_impacts(
        &mut bodies,
        &properties,
        &contacts,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    close_vector(bodies[0].linear_momentum_world(), Vec3::new(0.03, 0.0, 0.0));
    close_vector(bodies[0].angular_momentum_body(), Vec3::ZERO);
    close(report.support_work_j, -0.0003);
    close(report.contact.dissipated_energy, 1.50015);

    // Frictionless redundant normals still have unique body motion. Different
    // prescribed translations require a separate support-work witness before
    // this native owner can publish its combined energy report.
    let properties = [property(2.0, Vec3::new(1.0, 1.0, 1.0))];
    let initial = [body(
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::new(0.0, 0.0, -2.0),
    )];
    let first = plane(
        RigidImpactPoint {
            body: 0,
            arm_body: Vec3::ZERO,
        },
        0.0,
        0.0,
    );
    let mut second = first;
    second.id = 2;
    second.other = RigidImpactCounterpart::PrescribedPlane {
        point_world: Vec3::ZERO,
        velocity_world: Vec3::new(0.01, 0.0, 0.0),
    };
    let mut bodies = initial;
    assert!(matches!(
        apply_rigid_frictional_impacts(
            &mut bodies,
            &properties,
            &[first, second],
            controls(),
            &CancelGate::new()
        ),
        Err(RigidFrictionalImpactError::UnresolvedSupportWork)
    ));
    assert_eq!(bodies, initial);
}

#[test]
fn g1_sustained_tripod_support_resolves_load_transfer_without_a_supplied_pressure_centre() {
    use fs_mbd::frictional_impact::step::{
        RigidContactStepControls, step_rigid_frictional_contacts,
    };
    use fs_mbd::{Gravity, Wrench};

    let properties = [property(3.0, Vec3::new(1.0, 1.0, 1.0))];
    let mut contacts = tripod();
    for c in &mut contacts {
        c.activation_tolerance_m = 1e-7;
    }
    let config = RigidContactStepControls {
        impact: controls(),
        drift_energy_tolerance_j: 1e-10,
        drift_angular_momentum_tolerance_nms: 1e-10,
        work_energy_tolerance_j: 3e-8,
        linear_momentum_tolerance_ns: 1e-9,
        angular_momentum_tolerance_nms: 1e-9,
    };
    let mut sliding_errors = [0.0; 3];
    for (rung, count) in [8_usize, 16, 32].into_iter().enumerate() {
        let dt = 1.0 / count as f64;
        for horizontal_acceleration in [0.05_f64, 0.1, 0.2] {
            let mut bodies = [body(
                Vec3::new(0.0, 0.0, 0.25 + 1e-8),
                UnitQuaternion::IDENTITY,
                Vec3::ZERO,
            )];
            for step_index in 0..count {
                let report = step_rigid_frictional_contacts(
                    &mut bodies,
                    &properties,
                    &[Wrench::ZERO],
                    Gravity::new(Vec3::new(horizontal_acceleration, 0.0, -1.0)).unwrap(),
                    &contacts,
                    dt,
                    config,
                    |path| -> Result<(), &'static str> {
                        // Complete fixture geometry lies in this asymmetric
                        // triangular prism. Its height admits the declared mass
                        // distribution with the original COM at zero and I = 1.
                        // For every time, native COM height is affine and rotation
                        // moves each vertex by no more than R*|omega|*dt. This bound
                        // admits the whole path, rather than sampled endpoints.
                        let q = path.post_contact[0].pose().orientation();
                        let mut lowest = f64::INFINITY;
                        for (x, y) in [(-1.0, -1.0), (1.0, -1.0), (0.0, 1.0)] {
                            for z in [-0.25, 0.75] {
                                lowest = lowest.min(q.rotate_body_to_world(Vec3::new(x, y, z)).z);
                            }
                        }
                        let start_z = path.post_contact[0].pose().position_world().z;
                        let end_z = path.endpoint[0].pose().position_world().z;
                        // Original principal inertia is the identity, so body
                        // angular momentum equals angular speed throughout drift.
                        let omega = path.post_contact[0]
                            .angular_momentum_body()
                            .norm_squared()
                            .sqrt();
                        let lower = start_z.min(end_z) + lowest
                            - 2.5625_f64.sqrt() * omega * path.duration_s;
                        if lower < 0.0 {
                            return Err("whole prism path crosses its plane");
                        }
                        Ok(())
                    },
                    &CancelGate::new(),
                )
                .unwrap_or_else(|error| {
                    panic!(
                        "count={count} gx={horizontal_acceleration} step={step_index}: {error:?}"
                    )
                });
                assert!(matches!(
                    report.impact.contact.determinacy,
                    FrictionalContactDeterminacy::UniqueVelocity(_)
                ));
                close(report.total_energy_residual_j, 0.0);
                close(
                    report.impact.contact.impulses.iter().map(|p| p[0]).sum(),
                    3.0 * dt,
                );
            }
            let acceleration = (horizontal_acceleration - 0.1).max(0.0);
            close_vector(
                bodies[0].linear_momentum_world(),
                Vec3::new(3.0 * acceleration, 0.0, 0.0),
            );
            close_vector(bodies[0].angular_momentum_body(), Vec3::ZERO);
            close(
                bodies[0].pose().position_world().x,
                0.5 * acceleration * (1.0 + dt),
            );
            close(bodies[0].pose().position_world().z, 0.25 + 1e-8);
            if acceleration > 0.0 {
                sliding_errors[rung] = bodies[0].pose().position_world().x - 0.5 * acceleration;
                close(sliding_errors[rung], 0.5 * acceleration * dt);
            }
        }
    }
    for errors in sliding_errors.windows(2) {
        close(errors[0], 2.0 * errors[1]);
    }
}
