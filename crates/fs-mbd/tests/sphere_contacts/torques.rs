//! G1/G3/G4: native world couples through simultaneous sphere collisions.
use super::*;

fn chain() -> ([MassProperties; 3], [ContactSphere; 3], [RigidBodyState; 3]) {
    (
        [property(1.0, 0.5); 3],
        [10, 20, 30].map(|id| ContactSphere { id, radius_m: 0.5 }),
        [
            body(
                Vec3::new(-2.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                1.0,
                Vec3::ZERO,
            ),
            body(Vec3::ZERO, Vec3::ZERO, 1.0, Vec3::ZERO),
            body(
                Vec3::new(2.0, 0.0, 0.0),
                Vec3::new(-1.0, 0.0, 0.0),
                1.0,
                Vec3::ZERO,
            ),
        ],
    )
}

fn rotation_controls(step: f64) -> SphereTorqueControls {
    SphereTorqueControls {
        maximum_step_s: step,
        maximum_substeps: 16_384,
    }
}

fn chain_advance(
    bodies: &mut [RigidBodyState; 3],
    properties: &[MassProperties; 3],
    spheres: &[ContactSphere; 3],
    torques: &[Vec3; 3],
    duration: f64,
    rotation: SphereTorqueControls,
    gate: &CancelGate,
) -> Result<SphereContactReport, SphereContactError> {
    advance_wrenched_sphere_contacts(
        bodies,
        properties,
        spheres,
        &laws(spheres, 0.2, 1.0),
        Gravity::new(Vec3::new(0.0, -0.3, 0.0)).unwrap(),
        &[Vec3::new(0.0, 0.0, 0.2); 3],
        torques,
        duration,
        controls(),
        rotation,
        gate,
    )
}

#[test]
fn g1_torque_driven_simultaneous_chain_has_independent_friction_and_work() {
    let (properties, spheres, mut bodies) = chain();
    let torques = [Vec3::ZERO, Vec3::new(0.0, 0.0, 0.1), Vec3::ZERO];
    let duration = 1.25;
    let tail = duration - 1.0;
    let report = chain_advance(
        &mut bodies,
        &properties,
        &spheres,
        &torques,
        duration,
        rotation_controls(0.125),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(report.events.len(), 1);
    assert_eq!(report.events[0].pairs, vec![[10, 20], [20, 30]]);
    close(report.events[0].elapsed_s, 1.0);
    // Two simultaneous normal impulses are 2 N s each. The centre has
    // omega_z=1 just before impact. Both tangent slips are 1/2 m/s.
    // Their coupled mobility is [[7,3/2],[3/2,7]], giving j=1/17 at
    // each interface and total friction work -1/34 J. Pairwise independent
    // friction would instead use j=1/14 and fail this complete trajectory.
    let j = 1.0 / 17.0;
    let y = -0.15 * duration * duration;
    let z = 0.1 * duration * duration;
    for (i, expected) in [
        Vec3::new(-1.0 - tail, y - j * tail, z),
        Vec3::new(0.0, y, z),
        Vec3::new(1.0 + tail, y + j * tail, z),
    ]
    .into_iter()
    .enumerate()
    {
        vector(bodies[i].pose().position_world(), expected);
    }
    for (i, velocity) in [
        Vec3::new(-1.0, -j - 0.3 * duration, 0.2 * duration),
        Vec3::new(0.0, -0.3 * duration, 0.2 * duration),
        Vec3::new(1.0, j - 0.3 * duration, 0.2 * duration),
    ]
    .into_iter()
    .enumerate()
    {
        vector(bodies[i].linear_momentum_world(), velocity);
    }
    vector(spin_world(bodies[0]), Vec3::new(0.0, 0.0, -0.5 * j));
    vector(
        spin_world(bodies[1]),
        Vec3::new(0.0, 0.0, 0.1 * duration - j),
    );
    vector(spin_world(bodies[2]), Vec3::new(0.0, 0.0, -0.5 * j));
    let torque_work = 0.05 + (0.1 - j) * tail + 0.05 * tail * tail;
    close(report.couple_work_j, torque_work);
    vector(
        report.couple_angular_impulse_nms,
        Vec3::new(0.0, 0.0, 0.1 * duration),
    );
    close(report.unallocated_impact_loss_j, 1.0 / 34.0);
    close(report.gravity_work_j, 0.135 * duration * duration);
    close(report.centre_force_work_j, 0.06 * duration * duration);
    vector(
        report.gravity_angular_impulse_nms,
        Vec3::new(0.03 * duration.powi(3), 0.0, 0.0),
    );
    vector(
        report.centre_force_angular_impulse_nms,
        Vec3::new(-0.03 * duration.powi(3), 0.0, 0.0),
    );
    let angle = torque_work / 0.1;
    vector(
        bodies[1]
            .pose()
            .orientation()
            .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
        Vec3::new(angle.cos(), angle.sin(), 0.0),
    );
    close(report.total_energy_residual_j, 0.0);
    close(report.drift_energy_defect_j, 0.0);
    vector(report.linear_momentum_residual_ns, Vec3::ZERO);
    vector(report.angular_momentum_residual_nms, Vec3::ZERO);
    assert!(report.torque_substeps >= 10);
    assert!(report.minimum_path_gap_m >= 0.0);
}

#[test]
fn g3_world_frame_permutation_split_continuation_and_zero_couple_limits() {
    let (properties, spheres, initial) = chain();
    let torques = [Vec3::ZERO, Vec3::new(0.0, 0.0, 0.1), Vec3::ZERO];
    let mut reference = initial;
    let report = chain_advance(
        &mut reference,
        &properties,
        &spheres,
        &torques,
        1.25,
        rotation_controls(0.0625),
        &CancelGate::new(),
    )
    .unwrap();
    let q = UnitQuaternion::from_axis_angle(Vec3::new(0.3, -0.5, 0.7), 0.6).unwrap();
    let origin = Vec3::new(2.0, -3.0, 4.0);
    let mut transformed = initial.map(|b| {
        RigidBodyState::new(
            Pose::new(
                origin.add(q.rotate_body_to_world(b.pose().position_world())),
                q,
            )
            .unwrap(),
            q.rotate_body_to_world(b.linear_momentum_world()),
            b.angular_momentum_body(),
        )
        .unwrap()
    });
    let transformed_report = advance_wrenched_sphere_contacts(
        &mut transformed,
        &properties,
        &spheres,
        &laws(&spheres, 0.2, 1.0),
        Gravity::new(q.rotate_body_to_world(Vec3::new(0.0, -0.3, 0.0))).unwrap(),
        &[q.rotate_body_to_world(Vec3::new(0.0, 0.0, 0.2)); 3],
        &torques.map(|v| q.rotate_body_to_world(v)),
        1.25,
        controls(),
        rotation_controls(0.0625),
        &CancelGate::new(),
    )
    .unwrap();
    for i in 0..3 {
        vector(
            transformed[i].pose().position_world(),
            origin.add(q.rotate_body_to_world(reference[i].pose().position_world())),
        );
        vector(
            transformed[i].linear_momentum_world(),
            q.rotate_body_to_world(reference[i].linear_momentum_world()),
        );
        vector(
            spin_world(transformed[i]),
            q.rotate_body_to_world(spin_world(reference[i])),
        );
        for axis in [Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0)] {
            vector(
                transformed[i]
                    .pose()
                    .orientation()
                    .rotate_body_to_world(axis),
                q.rotate_body_to_world(
                    reference[i].pose().orientation().rotate_body_to_world(axis),
                ),
            );
        }
    }
    close(transformed_report.couple_work_j, report.couple_work_j);
    close(
        transformed_report.unallocated_impact_loss_j,
        report.unallocated_impact_loss_j,
    );
    vector(
        transformed_report.couple_angular_impulse_nms,
        q.rotate_body_to_world(report.couple_angular_impulse_nms),
    );

    let order = [2, 0, 1];
    let mut permuted = order.map(|i| initial[i]);
    let shuffled = chain_advance(
        &mut permuted,
        &order.map(|i| properties[i]),
        &order.map(|i| spheres[i]),
        &order.map(|i| torques[i]),
        1.25,
        rotation_controls(0.0625),
        &CancelGate::new(),
    )
    .unwrap();
    for (i, original) in order.into_iter().enumerate() {
        vector(
            permuted[i].pose().position_world(),
            reference[original].pose().position_world(),
        );
        vector(spin_world(permuted[i]), spin_world(reference[original]));
    }
    assert_eq!(shuffled.events[0].pairs, report.events[0].pairs);
    close(shuffled.couple_work_j, report.couple_work_j);

    let mut split = initial;
    let first = chain_advance(
        &mut split,
        &properties,
        &spheres,
        &torques,
        0.5,
        rotation_controls(0.0625),
        &CancelGate::new(),
    )
    .unwrap();
    let second = chain_advance(
        &mut split,
        &properties,
        &spheres,
        &torques,
        0.75,
        rotation_controls(0.0625),
        &CancelGate::new(),
    )
    .unwrap();
    for i in 0..3 {
        vector(
            split[i].pose().position_world(),
            reference[i].pose().position_world(),
        );
        vector(spin_world(split[i]), spin_world(reference[i]));
    }
    close(
        first.couple_work_j + second.couple_work_j,
        report.couple_work_j,
    );
    close(second.events[0].elapsed_s + 0.5, report.events[0].elapsed_s);

    let mut old = initial;
    let ordinary = advance_forced_sphere_contacts(
        &mut old,
        &properties,
        &spheres,
        &laws(&spheres, 0.2, 1.0),
        Gravity::new(Vec3::new(0.0, -0.3, 0.0)).unwrap(),
        &[Vec3::new(0.0, 0.0, 0.2); 3],
        1.25,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let mut zero = initial;
    let zero_report = chain_advance(
        &mut zero,
        &properties,
        &spheres,
        &[Vec3::ZERO; 3],
        1.25,
        SphereTorqueControls {
            maximum_step_s: 1e-20,
            maximum_substeps: 1,
        },
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(zero, old);
    assert_eq!(
        zero_report.total_energy_residual_j,
        ordinary.total_energy_residual_j
    );
    assert_eq!(
        zero_report.events[0].elapsed_s,
        ordinary.events[0].elapsed_s
    );
    assert_eq!(zero_report.torque_substeps, 0);
    assert_eq!(zero_report.couple_work_j, 0.0);
    assert_eq!(zero_report.couple_angular_impulse_nms, Vec3::ZERO);
}

fn rotate(q: [f64; 4], v: Vec3) -> Vec3 {
    let a = Vec3::new(q[1], q[2], q[3]);
    v.add(a.cross(a.cross(v).add(v.scale(q[0]))).scale(2.0))
}

fn euler_rhs(y: [f64; 8], inertia: Vec3, torque: Vec3) -> [f64; 8] {
    let q = [y[0], y[1], y[2], y[3]];
    let qv = Vec3::new(q[1], q[2], q[3]);
    let momentum = Vec3::new(y[4], y[5], y[6]);
    let omega = Vec3::new(
        momentum.x / inertia.x,
        momentum.y / inertia.y,
        momentum.z / inertia.z,
    );
    let body_torque = rotate([q[0], -q[1], -q[2], -q[3]], torque);
    let qdot = omega.scale(q[0]).add(qv.cross(omega)).scale(0.5);
    let ldot = body_torque.sub(omega.cross(momentum));
    [
        -0.5 * qv.dot(omega),
        qdot.x,
        qdot.y,
        qdot.z,
        ldot.x,
        ldot.y,
        ldot.z,
        body_torque.dot(omega),
    ]
}

fn euler_reference(mut state: [f64; 8], inertia: Vec3, torque: Vec3, duration: f64) -> [f64; 8] {
    let h = duration / 8192.0;
    let shifted =
        |a: [f64; 8], b: [f64; 8], factor: f64| std::array::from_fn(|i| a[i] + factor * b[i]);
    for _ in 0..8192 {
        let a = euler_rhs(state, inertia, torque);
        let b = euler_rhs(shifted(state, a, 0.5 * h), inertia, torque);
        let c = euler_rhs(shifted(state, b, 0.5 * h), inertia, torque);
        let d = euler_rhs(shifted(state, c, h), inertia, torque);
        state =
            std::array::from_fn(|i| state[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0);
    }
    state
}

#[test]
fn g1_full_tensor_driven_precession_refines_against_independent_euler_motion() {
    let axes = UnitQuaternion::from_axis_angle(Vec3::new(0.2, -0.3, 0.7), 0.6).unwrap();
    let inertia = Vec3::new(0.8, 1.2, 1.6);
    let properties =
        [MassProperties::from_principal_frame(2.0, Vec3::ZERO, inertia, axes).unwrap()];
    let spheres = [ContactSphere {
        id: 10,
        radius_m: 1.2,
    }];
    let lp = Vec3::new(0.5, -0.7, 0.9);
    let initial = RigidBodyState::new(
        Pose::new(Vec3::new(0.4, -0.2, 0.8), UnitQuaternion::IDENTITY).unwrap(),
        Vec3::new(0.2, -0.3, 0.4),
        axes.rotate_body_to_world(lp),
    )
    .unwrap();
    let torque = Vec3::new(0.4, 0.2, -0.1);
    let force = Vec3::new(0.1, -0.2, 0.05);
    let gravity = Gravity::new(Vec3::new(0.0, 0.0, -0.3)).unwrap();
    let duration = 1.25;
    let q = axes.components();
    let reference = euler_reference(
        [q[0], q[1], q[2], q[3], lp.x, lp.y, lp.z, 0.0],
        inertia,
        torque,
        duration,
    );
    let mut gates = controls();
    // Native asymmetric free rotation has a reported second-order energy
    // defect. This explicit bound admits the coarse reference refinement;
    // neither it nor phase error is reclassified as friction or heat.
    gates.drift_energy_tolerance_j = 1e-3;
    gates.total_energy_tolerance_j = 1e-3;
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let mut actual = [initial];
        let report = advance_wrenched_sphere_contacts(
            &mut actual,
            &properties,
            &spheres,
            &[],
            gravity,
            &[force],
            &[torque],
            duration,
            gates,
            rotation_controls(duration / n as f64),
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(report.events.len(), 0);
        assert_eq!(report.torque_substeps, n);
        assert_eq!(report.unallocated_impact_loss_j, 0.0);
        close(report.total_energy_residual_j, report.drift_energy_defect_j);
        vector(
            spin_world(actual[0]),
            spin_world(initial).add(torque.scale(duration)),
        );
        let acceleration = gravity.acceleration_world().add(force.scale(0.5));
        vector(
            actual[0].pose().position_world(),
            initial
                .pose()
                .position_world()
                .add(initial.linear_momentum_world().scale(0.5 * duration))
                .add(acceleration.scale(0.5 * duration * duration)),
        );
        let mut error = (report.couple_work_j - reference[7]).abs();
        for axis in [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ] {
            let got = actual[0]
                .pose()
                .orientation()
                .rotate_body_to_world(axes.rotate_body_to_world(axis));
            let expected = rotate(
                [reference[0], reference[1], reference[2], reference[3]],
                axis,
            );
            error += got.sub(expected).norm_squared().sqrt();
        }
        errors.push(error);
        if n == 32 {
            let diagonal = [MassProperties::new(2.0, Vec3::ZERO, inertia).unwrap()];
            let mut principal = [RigidBodyState::new(
                Pose::new(initial.pose().position_world(), axes).unwrap(),
                initial.linear_momentum_world(),
                lp,
            )
            .unwrap()];
            let other = advance_wrenched_sphere_contacts(
                &mut principal,
                &diagonal,
                &spheres,
                &[],
                gravity,
                &[force],
                &[torque],
                duration,
                gates,
                rotation_controls(duration / n as f64),
                &CancelGate::new(),
            )
            .unwrap();
            vector(spin_world(principal[0]), spin_world(actual[0]));
            vector(
                principal[0].pose().position_world(),
                actual[0].pose().position_world(),
            );
            close(other.couple_work_j, report.couple_work_j);
        }
    }
    for ratio in [errors[0] / errors[1], errors[1] / errors[2]] {
        assert!(
            (3.8..4.2).contains(&ratio),
            "second-order ratio {ratio}; errors {errors:?}"
        );
    }
    assert!(errors[2] < 2e-4, "refined error {errors:?}");
}

#[test]
fn g4_torque_budget_invalid_input_cancellation_and_retry_are_atomic() {
    let (properties, spheres, initial) = chain();
    let torques = [Vec3::ZERO, Vec3::new(0.0, 0.0, 0.1), Vec3::ZERO];
    let mut accepted = initial;
    let report = chain_advance(
        &mut accepted,
        &properties,
        &spheres,
        &torques,
        1.25,
        rotation_controls(0.125),
        &CancelGate::new(),
    )
    .unwrap();
    let mut failed = initial;
    assert!(matches!(
        chain_advance(
            &mut failed,
            &properties,
            &spheres,
            &torques,
            1.25,
            SphereTorqueControls {
                maximum_step_s: 0.125,
                maximum_substeps: report.torque_substeps - 1,
            },
            &CancelGate::new(),
        ),
        Err(SphereContactError::TorqueBudget)
    ));
    assert_eq!(failed, initial);
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        chain_advance(
            &mut failed,
            &properties,
            &spheres,
            &torques,
            1.25,
            rotation_controls(0.125),
            &gate,
        ),
        Err(SphereContactError::Cancelled)
    ));
    assert_eq!(failed, initial);
    for invalid in [
        SphereTorqueControls {
            maximum_step_s: 0.0,
            maximum_substeps: 1,
        },
        SphereTorqueControls {
            maximum_step_s: f64::NAN,
            maximum_substeps: 1,
        },
        SphereTorqueControls {
            maximum_step_s: 0.1,
            maximum_substeps: 0,
        },
        SphereTorqueControls {
            maximum_step_s: 0.1,
            maximum_substeps: 1_000_001,
        },
    ] {
        assert!(matches!(
            chain_advance(
                &mut failed,
                &properties,
                &spheres,
                &torques,
                1.25,
                invalid,
                &CancelGate::new(),
            ),
            Err(SphereContactError::Invalid(_))
        ));
        assert_eq!(failed, initial);
    }
    let bad_torques = [Vec3::ZERO, Vec3::new(0.0, f64::NAN, 0.1), Vec3::ZERO];
    assert!(matches!(
        chain_advance(
            &mut failed,
            &properties,
            &spheres,
            &bad_torques,
            1.25,
            rotation_controls(0.125),
            &CancelGate::new(),
        ),
        Err(SphereContactError::Invalid(_))
    ));
    assert_eq!(failed, initial);

    // This budget spans bodies as well as events: a completed first body's
    // torque rotation must not leak if the next body's required work refuses.
    let far_spheres = [
        ContactSphere {
            id: 1,
            radius_m: 0.5,
        },
        ContactSphere {
            id: 2,
            radius_m: 0.5,
        },
    ];
    let far_initial = [
        body(Vec3::ZERO, Vec3::ZERO, 1.0, Vec3::ZERO),
        body(Vec3::new(10.0, 0.0, 0.0), Vec3::ZERO, 1.0, Vec3::ZERO),
    ];
    let mut far = far_initial;
    assert!(matches!(
        advance_wrenched_sphere_contacts(
            &mut far,
            &[property(1.0, 0.5); 2],
            &far_spheres,
            &laws(&far_spheres, 0.0, 1.0),
            Gravity::ZERO,
            &[Vec3::ZERO; 2],
            &[Vec3::new(0.0, 0.0, 0.1); 2],
            0.5,
            controls(),
            SphereTorqueControls {
                maximum_step_s: 0.125,
                maximum_substeps: 7
            },
            &CancelGate::new(),
        ),
        Err(SphereContactError::TorqueBudget)
    ));
    assert_eq!(far, far_initial);
    assert!(matches!(
        advance_wrenched_sphere_contacts(
            &mut far,
            &[property(1.0, 0.5); 2],
            &far_spheres,
            &laws(&far_spheres, 0.0, 1.0),
            Gravity::ZERO,
            &[Vec3::ZERO; 2],
            &[Vec3::ZERO],
            0.5,
            controls(),
            rotation_controls(0.125),
            &CancelGate::new(),
        ),
        Err(SphereContactError::Invalid(_))
    ));
    assert_eq!(far, far_initial);
    // Equal opposite external couples have zero total angular impulse, but
    // both brake their original spins. Work must remain signed and body-local.
    let mut braking = [
        body(Vec3::ZERO, Vec3::ZERO, 1.0, Vec3::new(0.0, 0.0, 0.1)),
        body(
            Vec3::new(10.0, 0.0, 0.0),
            Vec3::ZERO,
            1.0,
            Vec3::new(0.0, 0.0, -0.1),
        ),
    ];
    let brake = advance_wrenched_sphere_contacts(
        &mut braking,
        &[property(1.0, 0.5); 2],
        &far_spheres,
        &laws(&far_spheres, 0.0, 1.0),
        Gravity::ZERO,
        &[Vec3::ZERO; 2],
        &[Vec3::new(0.0, 0.0, -0.1), Vec3::new(0.0, 0.0, 0.1)],
        0.5,
        controls(),
        SphereTorqueControls {
            maximum_step_s: 0.125,
            maximum_substeps: 8,
        },
        &CancelGate::new(),
    )
    .unwrap();
    vector(spin_world(braking[0]), Vec3::new(0.0, 0.0, 0.05));
    vector(spin_world(braking[1]), Vec3::new(0.0, 0.0, -0.05));
    close(brake.couple_work_j, -0.075);
    close(brake.kinetic_after_j - brake.kinetic_before_j, -0.075);
    assert_eq!(brake.couple_angular_impulse_nms, Vec3::ZERO);
    assert_eq!(brake.unallocated_impact_loss_j, 0.0);
    assert_eq!(brake.torque_substeps, 8);
    close(brake.total_energy_residual_j, 0.0);
    chain_advance(
        &mut failed,
        &properties,
        &spheres,
        &torques,
        1.25,
        rotation_controls(0.125),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(failed, accepted);
}
