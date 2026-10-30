//! G1/G3/G4: directional impulses retain original body inertia and transactions.
use super::*;
use fs_contact::frictional_multicontact::FrictionalContactDeterminacy;
use fs_contact::frictional_multicontact::elliptic::EllipticFriction;
use fs_mbd::frictional_impact::elliptic::*;

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
fn ellipse(coefficients_mu: [f64; 2], principal_axis: [f64; 2]) -> EllipticFriction {
    EllipticFriction {
        coefficients_mu,
        principal_axis,
    }
}

#[test]
fn g1_rotated_ellipse_full_inertia_and_offcentre_normal_coupling_match_physical_oracle() {
    let principal = UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), 0.45).unwrap();
    let properties = [MassProperties::from_principal_frame(
        2.0,
        Vec3::ZERO,
        Vec3::new(1.0, 1.4, 1.8),
        principal,
    )
    .unwrap()];
    let arm = Vec3::new(0.3, -0.2, -0.5);
    // Manufacture the endpoint law independently: p_t = -p_n Q v_t /
    // sqrt(v_t^T Q v_t), where Q = R diag(mu_1²,mu_2²) R^T.
    let mu = [0.3_f64, 0.6];
    let axis = [0.8, 0.6];
    let end = [2.0, 1.0];
    let principal_velocity = [
        axis[0] * end[0] + axis[1] * end[1],
        -axis[1] * end[0] + axis[0] * end[1],
    ];
    let support = (mu[0] * principal_velocity[0]).hypot(mu[1] * principal_velocity[1]);
    let scaled = [
        mu[0].powi(2) * principal_velocity[0],
        mu[1].powi(2) * principal_velocity[1],
    ];
    let pn = 1.2;
    let impulse_body = Vec3::new(
        -pn * (axis[0] * scaled[0] - axis[1] * scaled[1]) / support,
        -pn * (axis[1] * scaled[0] + axis[0] * scaled[1]) / support,
        pn,
    );
    let torque = arm.cross(impulse_body);
    let principal_torque = principal.rotate_world_to_body(torque);
    let omega = principal.rotate_body_to_world(Vec3::new(
        principal_torque.x,
        principal_torque.y / 1.4,
        principal_torque.z / 1.8,
    ));
    let point_increment = impulse_body.scale(0.5).add(omega.cross(arm));
    let before_relative = Vec3::new(end[0], end[1], 0.0).sub(point_increment);
    assert!(
        before_relative.z < -0.1,
        "fixture must be an incoming impact"
    );
    // Neither a circular cone nor normal-first clipping gives this impulse.
    assert!(impulse_body.y.abs() < 0.01);
    assert!((pn + 2.0 * before_relative.z).abs() > 0.05);
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.7).unwrap();
    for q in [UnitQuaternion::IDENTITY, rotation] {
        let boost = Vec3::new(0.2, -0.1, 0.3);
        let shift = Vec3::new(1.0, -2.0, 0.5);
        let initial = [body(
            shift.sub(q.rotate_body_to_world(arm)),
            q,
            q.rotate_body_to_world(before_relative)
                .add(boost)
                .scale(2.0),
        )];
        let mut c = plane(
            RigidImpactPoint {
                body: 0,
                arm_body: arm,
            },
            0.0,
            0.0,
        );
        c.normal_world = q.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
        c.tangent_world = q.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
        let actual_point = initial[0]
            .point_kinematics(properties[0], arm)
            .unwrap()
            .point_world;
        c.other = RigidImpactCounterpart::PrescribedPlane {
            // An infinite plane's reference anchor need not be its contact
            // point. The source must still receive the actual local position.
            point_world: if q == UnitQuaternion::IDENTITY {
                actual_point.sub(Vec3::new(1.0, 0.0, 0.0))
            } else {
                actual_point
            },
            velocity_world: boost,
        };
        let mut bodies = initial;
        let mut history = [7];
        let report = apply_rigid_elliptic_impacts_with(
            &mut bodies,
            &properties,
            &[with_friction(c, ())],
            &mut history,
            controls(),
            |state, query| {
                close_vector(query.point_world, actual_point);
                close_vector(query.counterpart_point_world, actual_point);
                Ok::<_, &str>((ellipse(mu, axis), state + 1))
            },
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(history, [8]);
        let actual = report.contact.impulses[0];
        close(actual[0], pn);
        close(actual[1], impulse_body.x);
        close(actual[2], impulse_body.y);
        for (a, expected) in report.contact.relative_velocities_after[0]
            .into_iter()
            .zip([0.0, end[0], end[1]])
        {
            close(a, expected);
        }
        close_vector(
            bodies[0].linear_momentum_world(),
            initial[0]
                .linear_momentum_world()
                .add(q.rotate_body_to_world(impulse_body)),
        );
        close_vector(bodies[0].angular_momentum_body(), torque);
        assert_eq!(bodies[0].pose(), initial[0].pose());
        let physical_work =
            impulse_body.dot(before_relative.add(Vec3::new(end[0], end[1], 0.0))) * 0.5;
        close(report.relative_impulse_work_j, physical_work);
        close(report.contact.dissipated_energy, -physical_work);
        close(
            report.support_work_j,
            q.rotate_body_to_world(impulse_body).dot(boost),
        );
        close(report.total_energy_residual_j, 0.0);
        close_vector(report.linear_momentum_residual_ns, Vec3::ZERO);
        close_vector(report.angular_momentum_residual_nms, Vec3::ZERO);
    }
}

fn chain() -> (
    [RigidBodyState; 3],
    [MassProperties; 3],
    [RigidFrictionalContact; 2],
) {
    let properties = [property(1.0, Vec3::new(1.0, 1.0, 1.0)); 3];
    let initial = [
        body(
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            Vec3::new(2.0, 1.0, 0.0),
        ),
        body(
            Vec3::new(2.0, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
        ),
        body(
            Vec3::new(4.0, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
        ),
    ];
    let contacts = std::array::from_fn(|i| RigidFrictionalContact {
        id: 10 + 10 * i as u64,
        point: RigidImpactPoint {
            body: i,
            arm_body: Vec3::new(1.0, 0.0, 0.0),
        },
        other: RigidImpactCounterpart::Body(RigidImpactPoint {
            body: i + 1,
            arm_body: Vec3::new(-1.0, 0.0, 0.0),
        }),
        normal_world: Vec3::new(-1.0, 0.0, 0.0),
        tangent_world: Vec3::new(0.0, 1.0, 0.0),
        activation_tolerance_m: 1e-10,
        friction_coefficient: 0.15,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
    });
    (initial, properties, contacts)
}

#[test]
fn g1_three_body_ellipse_chain_is_simultaneous_and_circular_limit_is_exact() {
    let (initial, properties, contacts) = chain();
    let directional = contacts.map(|c| with_friction(c, ellipse([0.15, 0.4], [1.0, 0.0])));
    let mut a = initial;
    let mut b = initial;
    let first = apply_rigid_elliptic_impacts(
        &mut a,
        &properties,
        &directional,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let second = apply_rigid_elliptic_impacts(
        &mut b,
        &properties,
        &[directional[1], directional[0]],
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(a, b);
    assert_eq!(first.contact_ids, second.contact_ids);
    assert_eq!(first.contact, second.contact);
    // Normal block [[2,-1],[-1,2]]; tangent block diag(4,4).
    // First shear slides at .15*(4/3)=.2; the second contact sticks.
    close(first.contact.impulses[0][0], 4.0 / 3.0);
    close(first.contact.impulses[1][0], 2.0 / 3.0);
    close(first.contact.impulses[0][1], -0.2);
    close(first.contact.impulses[1][1], 0.0);
    for state in a {
        close(state.linear_momentum_world().x, 2.0 / 3.0);
    }
    close(a[0].linear_momentum_world().y, 0.8);
    close(a[1].linear_momentum_world().y, 0.2);
    close(a[2].linear_momentum_world().y, 0.0);
    close(first.contact.dissipated_energy, 109.0 / 75.0);
    close_vector(first.linear_momentum_residual_ns, Vec3::ZERO);
    close_vector(first.angular_momentum_residual_nms, Vec3::ZERO);
    let mut circular = initial;
    let original = apply_rigid_frictional_impacts(
        &mut circular,
        &properties,
        &contacts,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let mut identical = initial;
    let isotropic = contacts.map(|c| with_friction(c, ellipse([0.15; 2], [0.8, 0.6])));
    let limit = apply_rigid_elliptic_impacts(
        &mut identical,
        &properties,
        &isotropic,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(circular, identical);
    assert_eq!(original.contact, limit.contact.transformed);
}

#[test]
fn g3_redundant_directional_support_retains_velocity_only_work_authority() {
    let properties = [property(3.0, Vec3::new(1.0, 1.0, 1.0))];
    let boost = Vec3::new(0.2, 0.0, 0.0);
    let mut bodies = [body(
        Vec3::new(0.0, 0.0, 0.25),
        UnitQuaternion::IDENTITY,
        Vec3::new(9.0, 0.0, -3.0).add(boost.scale(3.0)),
    )];
    let contacts: Vec<_> = [(-1.0, -1.0), (1.0, -1.0), (0.0, 1.0)]
        .into_iter()
        .enumerate()
        .map(|(i, (x, y))| {
            let mut c = plane(
                RigidImpactPoint {
                    body: 0,
                    arm_body: Vec3::new(x, y, -0.25),
                },
                0.0,
                0.0,
            );
            c.id = i as u64;
            c.other = RigidImpactCounterpart::PrescribedPlane {
                point_world: Vec3::ZERO,
                velocity_world: boost,
            };
            with_friction(c, ellipse([0.1, 0.2], [1.0, 0.0]))
        })
        .collect();
    let report = apply_rigid_elliptic_impacts(
        &mut bodies,
        &properties,
        &contacts,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert!(matches!(
        report.contact.transformed.determinacy,
        FrictionalContactDeterminacy::UniqueVelocity(_)
    ));
    close_vector(
        bodies[0].linear_momentum_world(),
        Vec3::new(8.7, 0.0, 0.0).add(boost.scale(3.0)),
    );
    close_vector(bodies[0].angular_momentum_body(), Vec3::ZERO);
    close(report.contact.impulses.iter().map(|p| p[0]).sum(), 3.0);
    close(report.contact.impulses.iter().map(|p| p[1]).sum(), -0.3);
    close(report.support_work_j, -0.3 * boost.x);
    close(
        report.contact.dissipated_energy,
        15.0 - 8.7_f64.powi(2) / 6.0,
    );
    close(report.total_energy_residual_j, 0.0);

    // Different support translations still require individual load authority.
    let properties = [property(2.0, Vec3::new(1.0, 1.0, 1.0))];
    let initial = [body(
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::new(0.0, 0.0, -2.0),
    )];
    let c = with_friction(
        plane(
            RigidImpactPoint {
                body: 0,
                arm_body: Vec3::ZERO,
            },
            0.0,
            0.0,
        ),
        ellipse([0.0; 2], [1.0, 0.0]),
    );
    let mut other = c;
    other.id = 2;
    other.other = RigidImpactCounterpart::PrescribedPlane {
        point_world: Vec3::ZERO,
        velocity_world: boost,
    };
    let mut bodies = initial;
    assert!(matches!(
        apply_rigid_elliptic_impacts(
            &mut bodies,
            &properties,
            &[c, other],
            controls(),
            &CancelGate::new()
        ),
        Err(RigidFrictionalImpactError::UnresolvedSupportWork)
    ));
    assert_eq!(bodies, initial);
}

#[test]
fn g4_source_queries_use_original_kinematics_and_all_failures_preserve_history() {
    let (initial, properties, contacts) = chain();
    let mut inactive = contacts[0];
    inactive.id = 5;
    inactive.point.arm_body = Vec3::ZERO;
    inactive.other = RigidImpactCounterpart::Body(RigidImpactPoint {
        body: 1,
        arm_body: Vec3::ZERO,
    });
    let geometry = [contacts[1], inactive, contacts[0]].map(|c| with_friction(c, ()));
    let accepted_history = [100_usize, 200, 300];
    let law = ellipse([0.15, 0.4], [1.0, 0.0]);
    for failure in 0..3 {
        let mut bodies = initial;
        let mut histories = accepted_history;
        let gate = CancelGate::new();
        let mut calls = Vec::new();
        let mut config = controls();
        if failure == 2 {
            config.contact.maximum_evaluations = 1;
        }
        let result = apply_rigid_elliptic_impacts_with(
            &mut bodies,
            &properties,
            &geometry,
            &mut histories,
            config,
            |history, query| {
                calls.push((query.contact_id, query.original_index));
                assert_eq!(
                    query.axes_world,
                    [
                        Vec3::new(-1.0, 0.0, 0.0),
                        Vec3::new(0.0, 1.0, 0.0),
                        Vec3::new(0.0, 0.0, -1.0)
                    ]
                );
                assert_eq!(
                    query.relative_velocity_m_s,
                    if query.contact_id == 10 {
                        [-2.0, 1.0, 0.0]
                    } else {
                        [0.0; 3]
                    }
                );
                if calls.len() == 2 && failure == 0 {
                    return Err("second source refused");
                }
                if calls.len() == 2 && failure == 1 {
                    gate.request();
                }
                Ok((law, history + 1))
            },
            &gate,
        );
        assert!(result.is_err(), "failure mode {failure} must refuse");
        assert_eq!(calls, [(10, 2), (20, 0)]);
        assert_eq!(bodies, initial);
        assert_eq!(histories, accepted_history);
    }
    let mut bodies = initial;
    let mut histories = accepted_history;
    let report = apply_rigid_elliptic_impacts_with(
        &mut bodies,
        &properties,
        &geometry,
        &mut histories,
        controls(),
        |history, _| Ok::<_, &str>((law, history + 1)),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(histories, [101, 200, 301]);
    assert_eq!(report.contact_ids, [5, 10, 20]);
    assert_eq!(report.contact.impulses[0], [0.0; 3]);
    close(report.contact.relative_velocities_after[0][0], 0.0);
    close(report.contact.relative_velocities_after[0][1], 0.6);

    // The contact law can pass while the original finite witness separation
    // creates a real angular couple. The late body audit must discard source state.
    let properties = [property(1.0, Vec3::new(1.0, 1.0, 1.0)); 2];
    let initial = [
        body(
            Vec3::new(0.0, 0.0, 1e-5),
            UnitQuaternion::IDENTITY,
            Vec3::new(3.0, 0.0, -1.0),
        ),
        body(Vec3::ZERO, UnitQuaternion::IDENTITY, Vec3::ZERO),
    ];
    let mut c = plane(
        RigidImpactPoint {
            body: 0,
            arm_body: Vec3::ZERO,
        },
        0.0,
        0.0,
    );
    c.other = RigidImpactCounterpart::Body(RigidImpactPoint {
        body: 1,
        arm_body: Vec3::ZERO,
    });
    c.activation_tolerance_m = 1e-4;
    let mut config = controls();
    config.geometry_tolerance_m = 1e-4;
    config.angular_momentum_tolerance_nms = 1e-12;
    let mut bodies = initial;
    let mut histories = [7];
    let result = apply_rigid_elliptic_impacts_with(
        &mut bodies,
        &properties,
        &[with_friction(c, ())],
        &mut histories,
        config,
        |history, _| Ok::<_, &str>((law, history + 1)),
        &CancelGate::new(),
    );
    assert!(matches!(result, Err(RigidFrictionalImpactError::Audit(_))));
    assert_eq!(bodies, initial);
    assert_eq!(histories, [7]);
}
