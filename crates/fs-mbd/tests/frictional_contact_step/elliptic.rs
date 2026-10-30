//! G1/G3/G4: directional laws act on the original loaded native trajectories.
use super::*;
use fs_contact::frictional_multicontact::elliptic::EllipticFriction;
use fs_mbd::frictional_impact::step::elliptic::*;

fn law() -> EllipticFriction {
    EllipticFriction {
        coefficients_mu: [0.25, 0.5],
        principal_axis: [1.0, 0.0],
    }
}
fn with_friction<F>(c: RigidFrictionalContact, friction: F) -> RigidFrictionalContact<F> {
    RigidFrictionalContact {
        id: c.id,
        point: c.point,
        other: c.other,
        normal_world: c.normal_world,
        tangent_world: c.tangent_world,
        activation_tolerance_m: c.activation_tolerance_m,
        friction_coefficient: friction,
        normal_law: c.normal_law,
    }
}
fn directed_plane(direction: Vec3, drive: f64, mu: f64) -> RigidFrictionalContact {
    let mut contact = plane(drive);
    contact.point.arm_body = direction
        .scale(HALF * (drive / 8.0).min(mu))
        .add(Vec3::new(0.0, 0.0, -HALF));
    contact
}

#[test]
fn g1_directional_loaded_block_sticks_slides_and_refines_actual_force_work() {
    for (direction, mu) in [
        (Vec3::new(1.0, 0.0, 0.0), 0.25),
        (Vec3::new(0.0, 1.0, 0.0), 0.5),
    ] {
        for drive in [4.0 * mu, 8.0 * mu, 6.0] {
            let mut bodies = [block()];
            let contact = with_friction(directed_plane(direction, drive, mu), law());
            let load = Gravity::new(direction.scale(drive).add(Vec3::new(0.0, 0.0, -8.0))).unwrap();
            for _ in 0..32 {
                let r = step_rigid_elliptic_contacts(
                    &mut bodies,
                    &[property()],
                    &[Wrench::ZERO],
                    load,
                    &[contact],
                    1.0 / 32.0,
                    controls(),
                    admit_box,
                    &CancelGate::new(),
                )
                .unwrap();
                close(r.impact.contact.impulses[0][0], 8.0 / 32.0, 2e-8);
                close(r.total_energy_residual_j, 0.0, 3e-8);
                assert!(r.path_admission >= 0.0 && r.endpoint_gaps[0].signed_gap_m >= 0.0);
            }
            let acceleration = (drive - 8.0 * mu).max(0.0);
            close_vec(
                bodies[0].linear_momentum_world(),
                direction.scale(acceleration),
                3e-7,
            );
            close_vec(bodies[0].angular_momentum_body(), Vec3::ZERO, 3e-9);
            close_vec(
                bodies[0].pose().position_world(),
                block()
                    .pose()
                    .position_world()
                    .add(direction.scale(0.5 * acceleration * (1.0 + 1.0 / 32.0))),
                3e-7,
            );
        }
    }
    let mut errors = Vec::new();
    for count in [8, 16, 32] {
        let dt = 1.0 / f64::from(count);
        let mut bodies = [block()];
        let mut kick = 0.0;
        let mut displacement_work = 0.0;
        let mut loss = 0.0;
        let contact = with_friction(directed_plane(Vec3::new(1.0, 0.0, 0.0), 6.0, 0.25), law());
        for _ in 0..count {
            let r = step_rigid_elliptic_contacts(
                &mut bodies,
                &[property()],
                &[Wrench::ZERO],
                gravity(6.0),
                &[contact],
                dt,
                controls(),
                admit_box,
                &CancelGate::new(),
            )
            .unwrap();
            kick += r.force_impulse_work_j;
            displacement_work += r.force_displacement_work_j;
            loss += r.impact.contact.dissipated_energy;
        }
        // Independent sums for a=6-0.25*8=4, x(1)=2 and v(1)=4:
        // kick=12+38dt; event loss=4+38dt includes the normal kick loss.
        close(bodies[0].pose().position_world().x, 2.0 * (1.0 + dt), 3e-7);
        close(bodies[0].linear_momentum_world().x, 4.0, 3e-7);
        close(kick, 12.0 + 38.0 * dt, 3e-6);
        close(loss, 4.0 + 38.0 * dt, 3e-6);
        close(displacement_work, 12.0 * (1.0 + dt), 3e-6);
        close(kick - displacement_work, 26.0 * dt, 3e-6);
        errors.push((bodies[0].pose().position_world().x - 2.0).abs());
    }
    for pair in errors.windows(2) {
        close(pair[0] / pair[1], 2.0, 2e-5);
    }
}

#[test]
fn g3_rotated_texture_and_moving_support_keep_original_pose_and_work() {
    let direction = Vec3::new(0.8, 0.6, 0.0);
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.7).unwrap();
    let boost = Vec3::new(0.3, -0.2, 0.1);
    let shift = Vec3::new(0.4, -0.7, 0.2);
    for q in [UnitQuaternion::IDENTITY, rotation] {
        let start = shift.add(q.rotate_body_to_world(Vec3::new(0.0, 0.0, HALF + CLEARANCE)));
        let mut bodies =
            [RigidBodyState::new(Pose::new(start, q).unwrap(), boost, Vec3::ZERO).unwrap()];
        let mut c = with_friction(
            directed_plane(direction, 6.0, 0.25),
            EllipticFriction {
                coefficients_mu: [0.25, 0.5],
                principal_axis: [0.8, 0.6],
            },
        );
        c.normal_world = q.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
        c.tangent_world = q.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
        let acceleration =
            q.rotate_body_to_world(direction.scale(6.0).add(Vec3::new(0.0, 0.0, -8.0)));
        let load = Gravity::new(acceleration).unwrap();
        for step in 0..8 {
            c.other = RigidImpactCounterpart::PrescribedPlane {
                point_world: shift.add(boost.scale(f64::from(step) / 8.0)),
                velocity_world: boost,
            };
            let r = step_rigid_elliptic_contacts(
                &mut bodies,
                &[property()],
                &[Wrench::ZERO],
                load,
                &[c],
                1.0 / 8.0,
                controls(),
                |path| {
                    let lower = admit_box(path)?;
                    let half = path.sample(0.5).map_err(|_| "original half drift")?;
                    close_vec(
                        half[0].pose().position_world(),
                        path.post_contact[0].pose().position_world().add(
                            path.post_contact[0]
                                .linear_momentum_world()
                                .scale(path.duration_s * 0.5),
                        ),
                        3e-10,
                    );
                    Ok::<_, &'static str>(lower)
                },
                &CancelGate::new(),
            )
            .unwrap();
            let expected_impulse =
                q.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0).sub(direction.scale(0.25)));
            close_vec(r.impact.support_impulse_ns, expected_impulse, 3e-8);
            close(r.impact.support_work_j, expected_impulse.dot(boost), 3e-8);
            close_vec(r.angular_momentum_residual_nms, Vec3::ZERO, 1e-9);
            close(r.total_energy_residual_j, 0.0, 3e-8);
        }
        close_vec(
            bodies[0].linear_momentum_world(),
            boost.add(q.rotate_body_to_world(direction.scale(4.0))),
            3e-7,
        );
        close_vec(
            bodies[0].pose().position_world(),
            start
                .add(boost)
                .add(q.rotate_body_to_world(direction.scale(2.25))),
            3e-7,
        );
        close_vec(bodies[0].angular_momentum_body(), Vec3::ZERO, 3e-9);
    }
}

#[test]
fn g3_equal_coefficients_preserve_scalar_step_states_and_every_work_account() {
    for mu in [0.0, 0.5] {
        let mut circular = [block()];
        let mut directional = circular;
        let mut c = directed_plane(Vec3::new(1.0, 0.0, 0.0), 6.0, mu);
        c.friction_coefficient = mu;
        let a = step_rigid_frictional_contacts(
            &mut circular,
            &[property()],
            &[Wrench::ZERO],
            gravity(6.0),
            &[c],
            0.125,
            controls(),
            admit_box,
            &CancelGate::new(),
        )
        .unwrap();
        let b = step_rigid_elliptic_contacts(
            &mut directional,
            &[property()],
            &[Wrench::ZERO],
            gravity(6.0),
            &[with_friction(
                c,
                EllipticFriction {
                    coefficients_mu: [mu; 2],
                    principal_axis: [0.6, 0.8],
                },
            )],
            0.125,
            controls(),
            admit_box,
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(circular, directional);
        assert_eq!(a.impact.contact.impulses, b.impact.contact.impulses);
        assert_eq!(
            a.impact.contact.relative_velocities_after,
            b.impact.contact.relative_velocities_after
        );
        assert_eq!(a.initial_gaps, b.initial_gaps);
        assert_eq!(a.endpoint_gaps, b.endpoint_gaps);
        assert_eq!(a.path_admission, b.path_admission);
        assert_eq!(
            [
                a.force_impulse_work_j,
                a.force_displacement_work_j,
                a.force_splitting_work_j,
                a.torque_impulse_work_j,
                a.kick_work_energy_residual_j,
                a.drift_energy_defect_j,
                a.total_energy_residual_j,
                a.impact.support_work_j,
                a.impact.contact.dissipated_energy
            ],
            [
                b.force_impulse_work_j,
                b.force_displacement_work_j,
                b.force_splitting_work_j,
                b.torque_impulse_work_j,
                b.kick_work_energy_residual_j,
                b.drift_energy_defect_j,
                b.total_energy_residual_j,
                b.impact.support_work_j,
                b.impact.contact.dissipated_energy
            ],
        );
        assert_eq!(a.linear_momentum_residual_ns, b.linear_momentum_residual_ns);
        assert_eq!(
            a.angular_momentum_residual_nms,
            b.angular_momentum_residual_nms
        );
    }
}

#[test]
fn g4_source_candidates_wait_for_complete_native_step_acceptance() {
    #[derive(Debug, PartialEq)]
    struct History(u64); // Deliberately no Clone requirement.
    let active = with_friction(directed_plane(Vec3::new(1.0, 0.0, 0.0), 6.0, 0.25), ());
    let mut inactive = active;
    inactive.id = 0;
    inactive.other = RigidImpactCounterpart::PrescribedPlane {
        point_world: Vec3::new(0.0, 0.0, -10.0),
        velocity_world: Vec3::ZERO,
    };
    let initial = [block()];
    let mut accepted = None;
    for outcome in [0, 1, 2, 3, 3] {
        let mut bodies = initial;
        let mut histories = [History(7), History(40)];
        let gate = CancelGate::new();
        let queried = Cell::new(0);
        let admitted = Cell::new(0);
        let result = step_rigid_elliptic_contacts_with(
            &mut bodies,
            &[property()],
            &[Wrench::ZERO],
            gravity(6.0),
            &[active, inactive],
            &mut histories,
            0.125,
            controls(),
            |state, query| {
                queried.set(queried.get() + 1);
                assert_eq!((query.contact_id, query.original_index, state.0), (1, 0, 7));
                close(query.relative_velocity_m_s[0], -1.0, 2e-14);
                close(query.relative_velocity_m_s[1], 0.75, 2e-14);
                close(query.relative_velocity_m_s[2], 0.0, 2e-14);
                close_vec(
                    query.point_world,
                    initial[0]
                        .pose()
                        .position_world()
                        .add(active.point.arm_body),
                    2e-14,
                );
                let friction = if outcome == 2 {
                    EllipticFriction {
                        coefficients_mu: [0.0, 0.5],
                        ..law()
                    }
                } else {
                    law()
                };
                Ok::<_, &'static str>((friction, History(state.0 + 1)))
            },
            |path| {
                admitted.set(admitted.get() + 1);
                if outcome == 0 {
                    return Err("original geometry owner declined path");
                }
                // The source's second candidate is inactive; one active box-plane
                // contact defines the complete collision boundary in this fixture.
                assert_eq!(path.contacts.len(), 2);
                assert!(admit_box(path)? >= 0.0);
                if outcome == 1 {
                    gate.request();
                }
                Ok::<_, &'static str>(17)
            },
            &gate,
        );
        assert_eq!(queried.get(), 1);
        match outcome {
            0 => assert!(matches!(result, Err(RigidContactStepError::Path(_)))),
            1 => assert!(matches!(result, Err(RigidContactStepError::Cancelled))),
            2 => assert!(matches!(result, Err(RigidContactStepError::Impact(_)))),
            _ => {
                let report = result.unwrap();
                assert_eq!(report.path_admission, 17);
                assert_eq!(histories, [History(8), History(40)]);
                assert!(bodies[0].pose().position_world().x > initial[0].pose().position_world().x);
                if let Some(before) = accepted {
                    assert_eq!(bodies, before);
                } else {
                    accepted = Some(bodies);
                }
                continue;
            }
        }
        assert_eq!(bodies, initial);
        assert_eq!(histories, [History(7), History(40)]);
        assert_eq!(admitted.get(), u32::from(outcome != 2));
    }
}
