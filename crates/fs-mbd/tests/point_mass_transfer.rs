use fs_mbd::assembly::capture::{
    PointMassCapture, PointMassMomentumCapture, capture_point_masses,
    capture_point_masses_with_momenta,
};
use fs_mbd::assembly::release::PointMassRelease;
use fs_mbd::assembly::transfer::{
    PointMassTransfer, PointMassTransferControls, PointMassTransferResult, transfer_point_masses,
};
use fs_mbd::assembly::{AssemblyControls, AssemblyError};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

fn controls(n: usize, impulses: usize) -> PointMassTransferControls {
    PointMassTransferControls {
        assembly: AssemblyControls {
            maximum_parts: n + 2,
            maximum_impulses: impulses,
            linear_momentum_tolerance_ns: 1e-10,
            angular_momentum_tolerance_nms: 1e-9,
            energy_tolerance_j: 1e-9,
        },
        attachment_tolerance_m: 1e-12,
        mass_tolerance_kg: 1e-12,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:.17e} != {b:.17e}");
}
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    close(a.x, b.x, tolerance);
    close(a.y, b.y, tolerance);
    close(a.z, b.z, tolerance);
}
fn rotation(axis: Vec3, angle: f64) -> UnitQuaternion {
    let scale = (0.5 * angle).sin() / axis.dot(axis).sqrt();
    UnitQuaternion::new(
        (0.5 * angle).cos(),
        axis.x * scale,
        axis.y * scale,
        axis.z * scale,
    )
    .unwrap()
}
fn compose(a: UnitQuaternion, b: UnitQuaternion) -> UnitQuaternion {
    let q = a.multiply(b);
    UnitQuaternion::new(q.w, q.x, q.y, q.z).unwrap()
}
fn state(p: MassProperties, x: Vec3, q: UnitQuaternion, v: Vec3, w: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(x, q).unwrap(),
        v.scale(p.mass()),
        p.angular_momentum_body_checked(w).unwrap(),
    )
    .unwrap()
}
fn simple(m: f64) -> MassProperties {
    MassProperties::new(m, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap()
}
fn transfer(m: f64, arm: Vec3, receiver_arm: Vec3, impulse: Vec3) -> PointMassTransfer {
    PointMassTransfer {
        release: PointMassRelease {
            mass_kg: m,
            attachment_arm_body_m: arm,
            impulse_on_parcel_body_ns: impulse,
        },
        receiver_attachment_arm_body_m: receiver_arm,
    }
}
fn inertia(a: MassProperties, b: MassProperties, tolerance: f64) {
    close(a.mass(), b.mass(), tolerance);
    for (a, b) in a
        .inertia_tensor_body()
        .m
        .into_iter()
        .zip(b.inertia_tensor_body().m)
    {
        close(a, b, tolerance);
    }
}
fn conservation(
    before: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    out: &PointMassTransferResult,
) {
    let after = [out.donor_state, out.receiver.state];
    let after_properties = [out.donor_properties, out.receiver.properties];
    let observer = Vec3::new(-0.7, 1.3, -0.9);
    let totals = |bodies: [RigidBodyState; 2], masses: [MassProperties; 2]| {
        let mut mass = 0.0;
        let mut p = Vec3::ZERO;
        let mut h = Vec3::ZERO;
        let mut energy = 0.0;
        for i in 0..2 {
            mass += masses[i].mass();
            p = p.add(bodies[i].linear_momentum_world());
            h = h
                .add(
                    bodies[i]
                        .pose()
                        .orientation()
                        .rotate_body_to_world(bodies[i].angular_momentum_body()),
                )
                .add(
                    bodies[i]
                        .pose()
                        .position_world()
                        .sub(observer)
                        .cross(bodies[i].linear_momentum_world()),
                );
            energy += bodies[i].kinetic_energy(masses[i]).unwrap();
        }
        (mass, p, h, energy)
    };
    let a = totals(before, properties);
    let b = totals(after, after_properties);
    close(b.0, a.0, 2e-13);
    vector(b.1, a.1, 2e-11);
    vector(b.2, a.2, 2e-9);
    close(
        b.3 - a.3,
        out.supplied_work_j - out.receiver.removed_relative_energy_j,
        2e-9,
    );
}

#[test]
fn canonical_handoff_preserves_momentum_that_velocity_roundtrip_changes() {
    let donor_properties = simple(98.0);
    let receiver_properties = simple(49.0);
    let donor = state(
        donor_properties,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Vec3::ZERO,
    );
    let receiver = state(
        receiver_properties,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Vec3::ZERO,
    );
    let input = transfer(49.0, Vec3::ZERO, Vec3::ZERO, Vec3::new(1.0, 0.0, 0.0));
    assert_ne!((1.0 / 49.0) * 49.0, 1.0);
    let mut c = controls(1, 2);
    c.assembly.linear_momentum_tolerance_ns = 0.0;
    let out = transfer_point_masses(
        &donor,
        donor_properties,
        &receiver,
        receiver_properties,
        &[input],
        c,
        || false,
    )
    .unwrap();
    assert_eq!(
        out.donor_state.linear_momentum_world(),
        Vec3::new(-1.0, 0.0, 0.0)
    );
    assert_eq!(
        out.receiver.state.linear_momentum_world(),
        Vec3::new(1.0, 0.0, 0.0)
    );
    close(out.supplied_work_j, 1.0 / 49.0, 2e-17);
    close(out.receiver.removed_relative_energy_j, 1.0 / 196.0, 2e-17);
    assert_eq!(out.balance.linear_momentum_residual_world_ns, Vec3::ZERO);
    conservation(
        [donor, receiver],
        [donor_properties, receiver_properties],
        &out,
    );
}

#[test]
fn finite_central_bodies_exchange_mass_recoil_and_separately_owned_energy() {
    let donor_properties = simple(4.0);
    let receiver_properties = simple(2.0);
    let donor = state(
        donor_properties,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::ZERO,
    );
    let receiver = state(
        receiver_properties,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::new(-2.0, 0.0, 0.0),
        Vec3::ZERO,
    );
    let input = transfer(1.0, Vec3::ZERO, Vec3::ZERO, Vec3::new(2.0, 0.0, 0.0));
    let out = transfer_point_masses(
        &donor,
        donor_properties,
        &receiver,
        receiver_properties,
        &[input],
        controls(1, 2),
        || false,
    )
    .unwrap();
    assert_eq!(out.donor_properties.mass(), 3.0);
    assert_eq!(out.receiver.properties.mass(), 3.0);
    assert_eq!(
        out.donor_state.linear_momentum_world(),
        Vec3::new(1.0, 0.0, 0.0)
    );
    assert_eq!(
        out.receiver.state.linear_momentum_world(),
        Vec3::new(-1.0, 0.0, 0.0)
    );
    close(out.supplied_work_j, 8.0 / 3.0, 2e-15);
    close(out.receiver.removed_relative_energy_j, 25.0 / 3.0, 8e-15);
    assert_eq!(out.mass_residual_kg, 0.0);
    assert_eq!(out.maximum_attachment_mismatch_m, 0.0);
    conservation(
        [donor, receiver],
        [donor_properties, receiver_properties],
        &out,
    );

    let boost = Vec3::new(1_099_511_627_776.0, 0.0, 0.0);
    let moved = |b: RigidBodyState, p: MassProperties| {
        RigidBodyState::new(
            b.pose(),
            b.linear_momentum_world().add(boost.scale(p.mass())),
            b.angular_momentum_body(),
        )
        .unwrap()
    };
    let boosted = transfer_point_masses(
        &moved(donor, donor_properties),
        donor_properties,
        &moved(receiver, receiver_properties),
        receiver_properties,
        &[input],
        controls(1, 2),
        || false,
    )
    .unwrap();
    assert_eq!(boosted.supplied_work_j, out.supplied_work_j);
    assert_eq!(
        boosted.receiver.removed_relative_energy_j,
        out.receiver.removed_relative_energy_j
    );
    close(boosted.balance.energy_residual_j, 0.0, 2e-14);
    assert_eq!(
        boosted.balance.linear_momentum_residual_world_ns,
        Vec3::ZERO
    );
    assert_eq!(
        boosted.balance.angular_momentum_residual_world_nms,
        Vec3::ZERO
    );
}

fn off_center() -> (
    [RigidBodyState; 2],
    [MassProperties; 2],
    [PointMassTransfer; 2],
) {
    let donor_properties = MassProperties::from_inertia_tensor(
        6.0,
        Vec3::ZERO,
        fs_ga::Mat3 {
            m: [
                22.0 / 3.0,
                -2.0 / 3.0,
                0.0,
                -2.0 / 3.0,
                19.0 / 3.0,
                0.0,
                0.0,
                0.0,
                32.0 / 3.0,
            ],
        },
    )
    .unwrap();
    let receiver_properties =
        MassProperties::new(5.0, Vec3::ZERO, Vec3::new(6.0, 7.0, 8.0)).unwrap();
    let donor = state(
        donor_properties,
        Vec3::new(2.0, -3.0, 1.0),
        rotation(Vec3::new(0.3, -0.4, 0.8), 0.6),
        Vec3::new(-0.5, 0.7, 0.8),
        Vec3::new(0.2, -0.3, 0.4),
    );
    let receiver = state(
        receiver_properties,
        Vec3::new(0.5, -1.0, 2.0),
        rotation(Vec3::new(-0.5, 0.9, 0.3), -0.7),
        Vec3::new(0.6, -0.2, 0.1),
        Vec3::new(-0.1, 0.2, -0.3),
    );
    let arm = |r: Vec3| {
        receiver.pose().orientation().rotate_world_to_body(
            donor
                .pose()
                .point_world_from_body(r)
                .unwrap()
                .sub(receiver.pose().position_world()),
        )
    };
    let a = Vec3::new(1.0, 0.0, 0.0);
    let b = Vec3::new(0.0, 1.0, 0.0);
    (
        [donor, receiver],
        [donor_properties, receiver_properties],
        [
            transfer(1.0, a, arm(a), Vec3::new(0.2, -0.1, 0.5)),
            transfer(2.0, b, arm(b), Vec3::new(-0.4, 0.7, -0.6)),
        ],
    )
}

#[test]
fn off_center_transfer_handles_independent_rotating_frames_and_full_inertias() {
    let (b, p, inputs) = off_center();
    let out =
        transfer_point_masses(&b[0], p[0], &b[1], p[1], &inputs, controls(2, 4), || false).unwrap();
    assert_eq!(out.donor_properties.mass(), 3.0);
    assert_eq!(out.receiver.properties.mass(), 8.0);
    vector(
        out.donor_center_shift_body_m,
        Vec3::new(-1.0 / 3.0, -2.0 / 3.0, 0.0),
        2e-15,
    );
    assert_eq!(
        out.donor_state.pose().orientation(),
        b[0].pose().orientation()
    );
    assert_eq!(
        out.receiver.state.pose().orientation(),
        b[1].pose().orientation()
    );
    assert!(out.receiver.removed_relative_energy_j > 0.0);
    conservation(b, p, &out);

    let turn = rotation(Vec3::new(0.2, 0.6, -0.3), 0.8);
    let shift = Vec3::new(-4.0, 2.0, 3.0);
    let boost = Vec3::new(16.0, -8.0, 4.0);
    let transformed = core::array::from_fn(|i| {
        RigidBodyState::new(
            Pose::new(
                turn.rotate_body_to_world(b[i].pose().position_world())
                    .add(shift),
                compose(turn, b[i].pose().orientation()),
            )
            .unwrap(),
            turn.rotate_body_to_world(b[i].linear_momentum_world())
                .add(boost.scale(p[i].mass())),
            b[i].angular_momentum_body(),
        )
        .unwrap()
    });
    let moved = transfer_point_masses(
        &transformed[0],
        p[0],
        &transformed[1],
        p[1],
        &inputs,
        controls(2, 4),
        || false,
    )
    .unwrap();
    inertia(moved.donor_properties, out.donor_properties, 0.0);
    inertia(moved.receiver.properties, out.receiver.properties, 0.0);
    vector(
        moved.donor_state.angular_momentum_body(),
        out.donor_state.angular_momentum_body(),
        2e-12,
    );
    vector(
        moved.receiver.state.angular_momentum_body(),
        out.receiver.state.angular_momentum_body(),
        3e-11,
    );
    close(moved.supplied_work_j, out.supplied_work_j, 0.0);
    close(
        moved.receiver.removed_relative_energy_j,
        out.receiver.removed_relative_energy_j,
        4e-11,
    );
    vector(
        moved.receiver.state.linear_momentum_world(),
        turn.rotate_body_to_world(out.receiver.state.linear_momentum_world())
            .add(boost.scale(out.receiver.properties.mass())),
        6e-13,
    );
    conservation(transformed, p, &moved);
}

#[test]
fn batch_permutation_and_parcel_subdivision_preserve_one_transfer_event() {
    let (b, p, inputs) = off_center();
    let run = |items: &[PointMassTransfer]| {
        transfer_point_masses(
            &b[0],
            p[0],
            &b[1],
            p[1],
            items,
            controls(items.len(), 2 * items.len()),
            || false,
        )
        .unwrap()
    };
    let out = run(&inputs);
    let split = |fraction: f64| PointMassTransfer {
        release: PointMassRelease {
            mass_kg: inputs[0].release.mass_kg * fraction,
            impulse_on_parcel_body_ns: inputs[0].release.impulse_on_parcel_body_ns.scale(fraction),
            ..inputs[0].release
        },
        ..inputs[0]
    };
    for changed in [
        run(&[inputs[1], inputs[0]]),
        run(&[split(0.25), inputs[1], split(0.75)]),
    ] {
        inertia(changed.donor_properties, out.donor_properties, 4e-12);
        inertia(changed.receiver.properties, out.receiver.properties, 5e-11);
        vector(
            changed.donor_state.linear_momentum_world(),
            out.donor_state.linear_momentum_world(),
            1e-13,
        );
        vector(
            changed.receiver.state.linear_momentum_world(),
            out.receiver.state.linear_momentum_world(),
            1e-13,
        );
        vector(
            changed.receiver.state.angular_momentum_body(),
            out.receiver.state.angular_momentum_body(),
            3e-11,
        );
        close(changed.supplied_work_j, out.supplied_work_j, 2e-12);
        close(
            changed.receiver.removed_relative_energy_j,
            out.receiver.removed_relative_energy_j,
            4e-11,
        );
        conservation(b, p, &changed);
    }
}

#[test]
fn attachment_admission_and_combined_angular_gate_reject_teleportation() {
    let p = simple(4.0);
    let donor = state(
        p,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Vec3::ZERO,
    );
    let displaced = state(
        p,
        Vec3::new(1e-6, 0.0, 0.0),
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Vec3::ZERO,
    );
    let input = transfer(1.0, Vec3::ZERO, Vec3::ZERO, Vec3::new(0.0, 1.0, 0.0));
    assert!(
        transfer_point_masses(&donor, p, &displaced, p, &[input], controls(1, 2), || false)
            .is_err()
    );
    let mut loose_geometry = controls(1, 2);
    loose_geometry.attachment_tolerance_m = 1e-5;
    loose_geometry.assembly.angular_momentum_tolerance_nms = 0.0;
    let error = transfer_point_masses(&donor, p, &displaced, p, &[input], loose_geometry, || false)
        .unwrap_err();
    let AssemblyError::Balance(balance) = error else {
        panic!("expected actual angular gate, got {error:?}");
    };
    vector(
        balance.angular_momentum_residual_world_nms,
        Vec3::new(0.0, 0.0, 1e-6),
        2e-21,
    );

    let far = state(
        p,
        Vec3::new(36_028_797_018_963_968.0, 0.0, 0.0),
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Vec3::ZERO,
    );
    let hidden_gap = transfer(1.0, Vec3::ZERO, Vec3::new(1.0, 0.0, 0.0), Vec3::ZERO);
    assert_eq!(
        far.pose().point_world_from_body(Vec3::ZERO).unwrap(),
        far.pose()
            .point_world_from_body(hidden_gap.receiver_attachment_arm_body_m)
            .unwrap()
    );
    assert!(
        transfer_point_masses(&far, p, &far, p, &[hidden_gap], controls(1, 1), || false).is_err()
    );
}

#[test]
fn identities_invalid_inventory_global_budgets_and_late_cancellation_are_atomic() {
    let (b, p, inputs) = off_center();
    let empty =
        transfer_point_masses(&b[0], p[0], &b[1], p[1], &[], controls(0, 0), || false).unwrap();
    assert_eq!(empty.donor_state, b[0]);
    assert_eq!(empty.donor_properties, p[0]);
    assert_eq!(empty.receiver.state, b[1]);
    assert_eq!(empty.receiver.properties, p[1]);
    assert_eq!(empty.supplied_work_j, 0.0);
    assert_eq!(empty.receiver.removed_relative_energy_j, 0.0);
    let zero = transfer(0.0, Vec3::ZERO, Vec3::ZERO, Vec3::ZERO);
    let no_op =
        transfer_point_masses(&b[0], p[0], &b[1], p[1], &[zero], controls(1, 0), || false).unwrap();
    assert_eq!(no_op, empty);
    assert!(matches!(
        transfer_point_masses(&b[0], p[0], &b[1], p[1], &[zero], controls(0, 0), || false),
        Err(AssemblyError::Budget {
            requested: 3,
            maximum: 2
        })
    ));
    assert!(matches!(
        transfer_point_masses(&b[0], p[0], &b[1], p[1], &inputs, controls(2, 3), || false),
        Err(AssemblyError::Budget {
            requested: 4,
            maximum: 3
        })
    ));
    for bad in [
        transfer(-1.0, Vec3::ZERO, Vec3::ZERO, Vec3::ZERO),
        transfer(0.0, Vec3::ZERO, Vec3::ZERO, Vec3::new(1.0, 0.0, 0.0)),
        transfer(0.0, Vec3::ZERO, Vec3::new(f64::NAN, 0.0, 0.0), Vec3::ZERO),
        transfer(f64::NAN, Vec3::ZERO, Vec3::ZERO, Vec3::ZERO),
    ] {
        assert!(
            transfer_point_masses(&b[0], p[0], &b[1], p[1], &[bad], controls(1, 2), || false)
                .is_err()
        );
    }
    let mut bad = controls(2, 4);
    bad.mass_tolerance_kg = -1.0;
    assert!(transfer_point_masses(&b[0], p[0], &b[1], p[1], &inputs, bad, || false).is_err());
    let mut all_polls = 0;
    let result = transfer_point_masses(&b[0], p[0], &b[1], p[1], &inputs, controls(2, 4), || {
        all_polls += 1;
        false
    })
    .unwrap();
    assert!(all_polls > 8 * inputs.len());
    let mut polls = 0;
    assert_eq!(
        transfer_point_masses(&b[0], p[0], &b[1], p[1], &inputs, controls(2, 4), || {
            polls += 1;
            polls == all_polls
        }),
        Err(AssemblyError::Cancelled)
    );
    assert_eq!(polls, all_polls);
    assert_eq!(
        transfer_point_masses(&b[0], p[0], &b[1], p[1], &inputs, controls(2, 4), || false).unwrap(),
        result
    );
    assert_eq!(b, off_center().0);
    assert_eq!(p, off_center().1);
}

#[test]
fn canonical_capture_shares_physics_and_rejects_momentum_without_mass() {
    let p = simple(2.0);
    let b = state(
        p,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::new(0.5, -0.5, 0.0),
        Vec3::new(0.0, 0.0, 0.25),
    );
    let velocity = [
        PointMassCapture {
            mass_kg: 1.0,
            attachment_arm_body_m: Vec3::new(1.0, 0.0, 0.0),
            velocity_world_m_s: Vec3::new(2.0, 1.0, 0.0),
        },
        PointMassCapture {
            mass_kg: 2.0,
            attachment_arm_body_m: Vec3::new(0.0, 1.0, 0.0),
            velocity_world_m_s: Vec3::new(-1.0, 0.5, 0.0),
        },
    ];
    let momenta = velocity.map(|v| PointMassMomentumCapture {
        mass_kg: v.mass_kg,
        attachment_arm_body_m: v.attachment_arm_body_m,
        linear_momentum_world_ns: v.velocity_world_m_s.scale(v.mass_kg),
    });
    let a = capture_point_masses(&b, p, &velocity, controls(2, 2).assembly, || false).unwrap();
    let c = capture_point_masses_with_momenta(&b, p, &momenta, controls(2, 2).assembly, || false)
        .unwrap();
    inertia(c.properties, a.properties, 0.0);
    vector(
        c.state.linear_momentum_world(),
        a.state.linear_momentum_world(),
        0.0,
    );
    vector(
        c.state.angular_momentum_body(),
        a.state.angular_momentum_body(),
        2e-14,
    );
    close(
        c.removed_relative_energy_j,
        a.removed_relative_energy_j,
        2e-14,
    );
    let zero = PointMassMomentumCapture {
        mass_kg: 0.0,
        attachment_arm_body_m: Vec3::ZERO,
        linear_momentum_world_ns: Vec3::ZERO,
    };
    let identity =
        capture_point_masses_with_momenta(&b, p, &[zero], controls(1, 0).assembly, || false)
            .unwrap();
    assert_eq!(identity.state, b);
    assert_eq!(identity.properties, p);
    for invalid in [
        PointMassMomentumCapture {
            linear_momentum_world_ns: Vec3::new(1.0, 0.0, 0.0),
            ..zero
        },
        PointMassMomentumCapture {
            mass_kg: -1.0,
            ..zero
        },
        PointMassMomentumCapture {
            mass_kg: 1.0,
            linear_momentum_world_ns: Vec3::new(f64::NAN, 0.0, 0.0),
            ..zero
        },
    ] {
        assert!(
            capture_point_masses_with_momenta(&b, p, &[invalid], controls(1, 1).assembly, || false)
                .is_err()
        );
    }
}
