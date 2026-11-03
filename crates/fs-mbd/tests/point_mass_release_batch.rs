use fs_mbd::assembly::capture::{PointMassCapture, capture_point_mass, capture_point_masses};
use fs_mbd::assembly::release::{
    PointMassRelease, PointMassReleaseBatchResult, release_point_mass, release_point_masses,
};
use fs_mbd::assembly::{AssemblyControls, AssemblyError};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

fn controls(inputs: usize) -> AssemblyControls {
    AssemblyControls {
        maximum_parts: inputs + 1,
        maximum_impulses: inputs,
        linear_momentum_tolerance_ns: 1e-9,
        angular_momentum_tolerance_nms: 1e-8,
        energy_tolerance_j: 1e-9,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}
fn vector(actual: Vec3, expected: Vec3, tolerance: f64) {
    close(actual.x, expected.x, tolerance);
    close(actual.y, expected.y, tolerance);
    close(actual.z, expected.z, tolerance);
}
fn inertia(actual: MassProperties, expected: MassProperties, tolerance: f64) {
    close(actual.mass(), expected.mass(), tolerance);
    for (a, b) in actual
        .inertia_tensor_body()
        .m
        .into_iter()
        .zip(expected.inertia_tensor_body().m)
    {
        close(a, b, tolerance);
    }
}
fn state(p: MassProperties, x: Vec3, q: UnitQuaternion, v: Vec3, w: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(x, q).unwrap(),
        v.scale(p.mass()),
        p.angular_momentum_body_checked(w).unwrap(),
    )
    .unwrap()
}
fn rotation(axis: Vec3, angle: f64) -> UnitQuaternion {
    UnitQuaternion::from_axis_angle(axis, angle).unwrap()
}
fn compose(q: UnitQuaternion, r: UnitQuaternion) -> UnitQuaternion {
    let [a, b, c, d] = q.components();
    let [e, f, g, h] = r.components();
    UnitQuaternion::new(
        a * e - b * f - c * g - d * h,
        a * f + b * e + c * h - d * g,
        a * g - b * h + c * e + d * f,
        a * h + b * g - c * f + d * e,
    )
    .unwrap()
}
fn balanced(spin: f64) -> (MassProperties, RigidBodyState, [PointMassRelease; 2]) {
    let p = MassProperties::new(4.0, Vec3::ZERO, Vec3::new(3.0, 6.0, 7.0)).unwrap();
    let b = state(
        p,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::new(3.0, -2.0, 1.0),
        Vec3::new(0.0, 0.0, spin),
    );
    let inputs = [
        PointMassRelease {
            mass_kg: 1.0,
            attachment_arm_body_m: Vec3::new(1.0, 0.0, 0.0),
            impulse_on_parcel_body_ns: Vec3::new(0.0, 1.0, 0.0),
        },
        PointMassRelease {
            mass_kg: 1.0,
            attachment_arm_body_m: Vec3::new(-1.0, 0.0, 0.0),
            impulse_on_parcel_body_ns: Vec3::new(0.0, -1.0, 0.0),
        },
    ];
    (p, b, inputs)
}
fn asymmetric_properties() -> MassProperties {
    // Retained M=3, COM=(-1/3,-2/3,0), I=diag(4,5,6), plus
    // m1=1 at ex and m2=2 at ey. Total COM is zero and total mass is six.
    let mut tensor = MassProperties::new(
        6.0,
        Vec3::ZERO,
        Vec3::new(22.0 / 3.0, 19.0 / 3.0, 32.0 / 3.0),
    )
    .unwrap()
    .inertia_tensor_body();
    tensor.m[1] = -2.0 / 3.0;
    tensor.m[3] = -2.0 / 3.0;
    MassProperties::from_inertia_tensor(6.0, Vec3::ZERO, tensor).unwrap()
}
fn asymmetric_inputs() -> [PointMassRelease; 2] {
    [
        PointMassRelease {
            mass_kg: 1.0,
            attachment_arm_body_m: Vec3::new(1.0, 0.0, 0.0),
            impulse_on_parcel_body_ns: Vec3::new(0.0, 0.0, 3.0),
        },
        PointMassRelease {
            mass_kg: 2.0,
            attachment_arm_body_m: Vec3::new(0.0, 1.0, 0.0),
            impulse_on_parcel_body_ns: Vec3::ZERO,
        },
    ]
}
fn release(
    b: &RigidBodyState,
    p: MassProperties,
    inputs: &[PointMassRelease],
) -> PointMassReleaseBatchResult {
    release_point_masses(b, p, inputs, controls(inputs.len()), || false).unwrap()
}
// Direct fixed-observer sums from all actual returned states. These do not
// use the event's work gate or its internal center/inertia formulas.
fn conservation(before: RigidBodyState, p: MassProperties, out: &PointMassReleaseBatchResult) {
    let observer = Vec3::new(-0.7, 1.1, 0.4);
    let before_angular = before
        .pose()
        .orientation()
        .rotate_body_to_world(before.angular_momentum_body())
        .add(
            before
                .pose()
                .position_world()
                .sub(observer)
                .cross(before.linear_momentum_world()),
        );
    let mut linear = out.state.linear_momentum_world();
    let mut angular = out
        .state
        .pose()
        .orientation()
        .rotate_body_to_world(out.state.angular_momentum_body())
        .add(
            out.state
                .pose()
                .position_world()
                .sub(observer)
                .cross(out.state.linear_momentum_world()),
        );
    let mut energy = out.state.kinetic_energy(out.properties).unwrap();
    let mut mass = out.properties.mass();
    for parcel in &out.parcels {
        mass += parcel.mass_kg;
        linear = linear.add(parcel.linear_momentum_world_ns);
        angular = angular.add(
            parcel
                .position_world_m
                .sub(observer)
                .cross(parcel.linear_momentum_world_ns),
        );
        energy += 0.5
            * parcel
                .linear_momentum_world_ns
                .dot(parcel.velocity_world_m_s().unwrap());
    }
    close(mass, p.mass(), 2e-14);
    vector(linear, before.linear_momentum_world(), 5e-12);
    vector(angular, before_angular, 2e-9);
    close(
        energy,
        before.kinetic_energy(p).unwrap() + out.supplied_work_j,
        2e-9,
    );
}

#[test]
fn simultaneous_jets_keep_shared_translation_and_rotation_work_cross_terms() {
    for spin in [0.0, 1.0] {
        let (p, b, inputs) = balanced(spin);
        let out = release(&b, p, &inputs);
        let retained = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(3.0, 4.0, 5.0)).unwrap();
        inertia(out.properties, retained, 2e-14);
        vector(out.center_shift_body_m, Vec3::ZERO, 0.0);
        vector(
            out.state.linear_momentum_world(),
            Vec3::new(6.0, -4.0, 2.0),
            2e-14,
        );
        vector(
            out.state.angular_momentum_body(),
            Vec3::new(0.0, 0.0, 5.0 * spin - 2.0),
            2e-14,
        );
        vector(
            out.parcels[0].linear_momentum_world_ns,
            Vec3::new(3.0, -1.0 + spin, 1.0),
            2e-14,
        );
        vector(
            out.parcels[1].linear_momentum_world_ns,
            Vec3::new(3.0, -3.0 - spin, 1.0),
            2e-14,
        );
        // Opposed impulses: parcel energy 1, retained translation zero,
        // retained rotation 4/(2*5). Independent impulse work would be 1.7.
        close(out.supplied_work_j, 7.0 / 5.0, 2e-14);
        assert!((out.supplied_work_j - 1.7).abs() > 0.29);
        conservation(b, p, &out);

        let co_directed = [
            inputs[0],
            PointMassRelease {
                impulse_on_parcel_body_ns: inputs[0].impulse_on_parcel_body_ns,
                ..inputs[1]
            },
        ];
        let co = release(&b, p, &co_directed);
        vector(
            co.state.linear_momentum_world(),
            Vec3::new(6.0, -6.0, 2.0),
            2e-14,
        );
        vector(
            co.state.angular_momentum_body(),
            Vec3::new(0.0, 0.0, 5.0 * spin),
            2e-14,
        );
        // Co-directed impulses: parcel energy 1, retained translation 1,
        // retained rotation zero. The cross terms now increase total work.
        close(co.supplied_work_j, 2.0, 2e-14);
        conservation(b, p, &co);
    }
    // Sequential releases are different physical events: the second parcel
    // inherits an already recoiling body, not the original common snapshot.
    let (p, b, inputs) = balanced(0.0);
    let first = release_point_mass(&b, p, inputs[0], controls(1), || false).unwrap();
    let second = release_point_mass(
        &first.state,
        first.properties,
        PointMassRelease {
            attachment_arm_body_m: inputs[1]
                .attachment_arm_body_m
                .sub(first.center_shift_body_m),
            ..inputs[1]
        },
        controls(1),
        || false,
    )
    .unwrap();
    close(
        first.supplied_work_j + second.supplied_work_j,
        569.0 / 340.0,
        3e-14,
    );
    assert!((first.supplied_work_j + second.supplied_work_j - 7.0 / 5.0).abs() > 0.27);
}

#[test]
fn asymmetric_detachment_shifts_com_and_applies_recoil_about_retained_center() {
    let p = asymmetric_properties();
    let b = state(
        p,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Vec3::ZERO,
    );
    let inputs = asymmetric_inputs();
    let mut c = controls(2);
    c.maximum_impulses = 1; // The second parcel detaches without an impulse.
    let out = release_point_masses(&b, p, &inputs, c, || false).unwrap();
    inertia(
        out.properties,
        MassProperties::new(3.0, Vec3::ZERO, Vec3::new(4.0, 5.0, 6.0)).unwrap(),
        3e-12,
    );
    vector(
        out.center_shift_body_m,
        Vec3::new(-1.0 / 3.0, -2.0 / 3.0, 0.0),
        2e-15,
    );
    vector(
        out.state.linear_momentum_world(),
        Vec3::new(0.0, 0.0, -3.0),
        0.0,
    );
    vector(
        out.state.angular_momentum_body(),
        Vec3::new(-2.0, 4.0, 0.0),
        3e-14,
    );
    vector(
        out.properties
            .angular_velocity_body_checked(out.state.angular_momentum_body())
            .unwrap(),
        Vec3::new(-0.5, 0.8, 0.0),
        4e-13,
    );
    // Parcel recoil 4.5 J + retained translation 1.5 J + rotation 2.1 J.
    close(out.supplied_work_j, 81.0 / 10.0, 4e-12);
    for (input, parcel) in inputs.iter().zip(&out.parcels) {
        vector(
            out.state
                .pose()
                .point_world_from_body(input.attachment_arm_body_m.sub(out.center_shift_body_m))
                .unwrap(),
            parcel.position_world_m,
            2e-15,
        );
    }
    conservation(b, p, &out);
}

#[test]
fn zero_impulse_batch_inherits_one_rigid_velocity_field_and_singleton_agrees() {
    let (p, b, inputs) = balanced(1.0);
    let unforced = inputs.map(|input| PointMassRelease {
        impulse_on_parcel_body_ns: Vec3::ZERO,
        ..input
    });
    let mut c = controls(2);
    c.maximum_impulses = 0;
    let out = release_point_masses(&b, p, &unforced, c, || false).unwrap();
    assert_eq!(out.supplied_work_j, 0.0);
    for (input, parcel) in unforced.iter().zip(&out.parcels) {
        let old = b.point_kinematics(p, input.attachment_arm_body_m).unwrap();
        vector(parcel.position_world_m, old.point_world, 0.0);
        vector(
            parcel.velocity_world_m_s().unwrap(),
            old.point_velocity_world,
            2e-14,
        );
    }
    vector(
        out.properties
            .angular_velocity_body_checked(out.state.angular_momentum_body())
            .unwrap(),
        p.angular_velocity_body_checked(b.angular_momentum_body())
            .unwrap(),
        2e-14,
    );
    conservation(b, p, &out);
    for input in [
        inputs[0],
        PointMassRelease {
            attachment_arm_body_m: Vec3::ZERO,
            ..inputs[0]
        },
    ] {
        let batch = release(&b, p, &[input]);
        let single = release_point_mass(&b, p, input, controls(1), || false).unwrap();
        inertia(batch.properties, single.properties, 3e-13);
        vector(
            batch.state.pose().position_world(),
            single.state.pose().position_world(),
            2e-15,
        );
        vector(
            batch.state.linear_momentum_world(),
            single.state.linear_momentum_world(),
            3e-14,
        );
        vector(
            batch.state.angular_momentum_body(),
            single.state.angular_momentum_body(),
            3e-13,
        );
        vector(
            batch.parcels[0].linear_momentum_world_ns,
            single.parcel.unwrap().linear_momentum_world_ns,
            2e-14,
        );
        close(batch.supplied_work_j, single.supplied_work_j, 3e-13);
    }
}

#[test]
fn full_tensor_batch_is_covariant_and_recapture_restores_original_body() {
    let p = asymmetric_properties();
    let b = state(
        p,
        Vec3::new(2.0, -3.0, 1.0),
        rotation(Vec3::new(0.3, -0.4, 0.8), 0.6),
        Vec3::new(-0.5, 0.7, 0.8),
        Vec3::new(0.2, -0.3, 0.4),
    );
    let mut inputs = asymmetric_inputs();
    inputs[0].impulse_on_parcel_body_ns = Vec3::new(0.2, -0.1, 0.5);
    inputs[1].impulse_on_parcel_body_ns = Vec3::new(-0.4, 0.7, -0.6);
    let out = release(&b, p, &inputs);
    conservation(b, p, &out);

    // Simultaneous recapture uses ONE retained-body COM for every parcel.
    let incoming: Vec<_> = inputs
        .iter()
        .zip(&out.parcels)
        .map(|(input, parcel)| PointMassCapture {
            mass_kg: parcel.mass_kg,
            attachment_arm_body_m: input.attachment_arm_body_m.sub(out.center_shift_body_m),
            velocity_world_m_s: parcel.velocity_world_m_s().unwrap(),
        })
        .collect();
    let captured = capture_point_masses(
        &out.state,
        out.properties,
        &incoming,
        controls(incoming.len()),
        || false,
    )
    .unwrap();
    inertia(captured.properties, p, 8e-12);
    vector(
        captured.state.pose().position_world(),
        b.pose().position_world(),
        3e-14,
    );
    vector(
        captured.state.linear_momentum_world(),
        b.linear_momentum_world(),
        3e-14,
    );
    vector(
        captured.state.angular_momentum_body(),
        b.angular_momentum_body(),
        8e-12,
    );
    close(
        captured.removed_relative_energy_j,
        out.supplied_work_j,
        3e-12,
    );

    let mut recovered = out.state;
    let mut recovered_properties = out.properties;
    let mut shift = out.center_shift_body_m;
    let mut removed_work = 0.0;
    for (input, parcel) in inputs.iter().zip(&out.parcels) {
        let captured = capture_point_mass(
            &recovered,
            recovered_properties,
            PointMassCapture {
                mass_kg: parcel.mass_kg,
                attachment_arm_body_m: input.attachment_arm_body_m.sub(shift),
                velocity_world_m_s: parcel.velocity_world_m_s().unwrap(),
            },
            controls(1),
            || false,
        )
        .unwrap();
        recovered = captured.state;
        recovered_properties = captured.properties;
        shift = shift.add(captured.center_shift_body_m);
        removed_work += captured.removed_relative_energy_j;
    }
    inertia(recovered_properties, p, 8e-12);
    vector(
        recovered.pose().position_world(),
        b.pose().position_world(),
        3e-14,
    );
    vector(
        recovered.linear_momentum_world(),
        b.linear_momentum_world(),
        3e-14,
    );
    vector(
        recovered.angular_momentum_body(),
        b.angular_momentum_body(),
        8e-12,
    );
    close(removed_work, out.supplied_work_j, 3e-12);

    let turn = rotation(Vec3::new(-0.5, 0.9, 0.3), -0.7);
    let translation = Vec3::new(4.0, -2.0, 3.0);
    let boost = Vec3::new(64.0, -32.0, 16.0);
    let transformed = RigidBodyState::new(
        Pose::new(
            turn.rotate_body_to_world(b.pose().position_world())
                .add(translation),
            compose(turn, b.pose().orientation()),
        )
        .unwrap(),
        turn.rotate_body_to_world(b.linear_momentum_world())
            .add(boost.scale(p.mass())),
        b.angular_momentum_body(),
    )
    .unwrap();
    let moved = release(&transformed, p, &inputs);
    assert_eq!(moved.properties, out.properties);
    assert_eq!(
        moved.state.angular_momentum_body(),
        out.state.angular_momentum_body()
    );
    assert_eq!(moved.supplied_work_j, out.supplied_work_j);
    vector(
        moved.state.pose().position_world(),
        turn.rotate_body_to_world(out.state.pose().position_world())
            .add(translation),
        4e-14,
    );
    vector(
        moved.state.linear_momentum_world(),
        turn.rotate_body_to_world(out.state.linear_momentum_world())
            .add(boost.scale(out.properties.mass())),
        2e-13,
    );
    for (a, b) in moved.parcels.iter().zip(&out.parcels) {
        vector(
            a.position_world_m,
            turn.rotate_body_to_world(b.position_world_m)
                .add(translation),
            4e-14,
        );
        vector(
            a.linear_momentum_world_ns,
            turn.rotate_body_to_world(b.linear_momentum_world_ns)
                .add(boost.scale(b.mass_kg)),
            1e-13,
        );
    }
    conservation(transformed, p, &moved);
}

#[test]
fn subdivision_and_permutation_preserve_simultaneous_release() {
    let p = asymmetric_properties();
    let b = state(
        p,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::new(0.5, 0.25, -0.5),
        Vec3::new(0.2, 0.3, -0.1),
    );
    let inputs = asymmetric_inputs();
    let original = release(&b, p, &inputs);
    let reversed = release(&b, p, &[inputs[1], inputs[0]]);
    inertia(reversed.properties, original.properties, 1e-13);
    vector(
        reversed.state.linear_momentum_world(),
        original.state.linear_momentum_world(),
        3e-14,
    );
    vector(
        reversed.state.angular_momentum_body(),
        original.state.angular_momentum_body(),
        3e-13,
    );
    close(reversed.supplied_work_j, original.supplied_work_j, 3e-13);
    assert_eq!(reversed.parcels[0], original.parcels[1]);
    assert_eq!(reversed.parcels[1], original.parcels[0]);

    let split = [
        PointMassRelease {
            mass_kg: 0.25,
            impulse_on_parcel_body_ns: inputs[0].impulse_on_parcel_body_ns.scale(0.25),
            ..inputs[0]
        },
        inputs[1],
        PointMassRelease {
            mass_kg: 0.75,
            impulse_on_parcel_body_ns: inputs[0].impulse_on_parcel_body_ns.scale(0.75),
            ..inputs[0]
        },
    ];
    let divided = release(&b, p, &split);
    inertia(divided.properties, original.properties, 2e-12);
    vector(
        divided.center_shift_body_m,
        original.center_shift_body_m,
        2e-15,
    );
    vector(
        divided.state.linear_momentum_world(),
        original.state.linear_momentum_world(),
        3e-14,
    );
    vector(
        divided.state.angular_momentum_body(),
        original.state.angular_momentum_body(),
        2e-12,
    );
    vector(
        divided.parcels[0]
            .linear_momentum_world_ns
            .add(divided.parcels[2].linear_momentum_world_ns),
        original.parcels[0].linear_momentum_world_ns,
        2e-14,
    );
    close(divided.supplied_work_j, original.supplied_work_j, 2e-12);
    conservation(b, p, &divided);
}

#[test]
fn common_boost_keeps_small_work_and_actual_first_moment_detects_lost_geometry() {
    let (p, b, inputs) = balanced(1.0);
    let plain = release(&b, p, &inputs);
    let boost = Vec3::new(1_099_511_627_776.0, 0.0, 0.0);
    let boosted = RigidBodyState::new(
        b.pose(),
        b.linear_momentum_world().add(boost.scale(p.mass())),
        b.angular_momentum_body(),
    )
    .unwrap();
    let moved = release(&boosted, p, &inputs);
    assert_eq!(moved.supplied_work_j, plain.supplied_work_j);
    assert_eq!(moved.balance, plain.balance);
    assert_eq!(
        boosted.kinetic_energy(p).unwrap(),
        boosted.kinetic_energy(p).unwrap() + moved.supplied_work_j
    );

    // Zero spin and zero impulses: ONLY the actual returned first-moment
    // defect times the common boost exposes this unrepresentable geometry.
    let p = MassProperties::new(4.0, Vec3::ZERO, Vec3::new(3.0, 40.0, 41.0)).unwrap();
    let inputs = [
        PointMassRelease {
            mass_kg: 1.0,
            attachment_arm_body_m: Vec3::new(5.0, 0.0, 0.0),
            impulse_on_parcel_body_ns: Vec3::ZERO,
        },
        PointMassRelease {
            mass_kg: 1.0,
            attachment_arm_body_m: Vec3::new(-3.0, 0.0, 0.0),
            impulse_on_parcel_body_ns: Vec3::ZERO,
        },
    ];
    let b = state(
        p,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::ZERO,
    );
    let good = release(&b, p, &inputs);
    vector(good.center_shift_body_m, Vec3::new(-1.0, 0.0, 0.0), 0.0);
    assert_eq!(good.supplied_work_j, 0.0);
    let distant = state(
        p,
        Vec3::new(36_028_797_018_963_968.0, 0.0, 0.0),
        UnitQuaternion::IDENTITY,
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::ZERO,
    );
    let error = release_point_masses(&distant, p, &inputs, controls(2), || false).unwrap_err();
    let AssemblyError::Balance(balance) = error else {
        panic!("expected actual geometry balance, got {error:?}")
    };
    vector(balance.linear_momentum_residual_world_ns, Vec3::ZERO, 0.0);
    close(balance.energy_residual_j, 0.0, 0.0);
    vector(
        balance.angular_momentum_residual_world_nms,
        Vec3::new(0.0, 0.0, 4.0),
        0.0,
    );
}

#[test]
fn identities_budgets_invalid_inertia_and_late_cancellation_are_atomic() {
    let (p, b, inputs) = balanced(1.0);
    let empty = release_point_masses(&b, p, &[], controls(0), || false).unwrap();
    assert_eq!(empty.state, b);
    assert_eq!(empty.properties, p);
    assert!(empty.parcels.is_empty());
    assert_eq!(empty.supplied_work_j, 0.0);
    let zero = PointMassRelease {
        mass_kg: 0.0,
        impulse_on_parcel_body_ns: Vec3::ZERO,
        ..inputs[0]
    };
    let zeros = release_point_masses(&b, p, &[zero, zero], controls(2), || false).unwrap();
    assert_eq!(zeros, empty);
    let mixed = release(&b, p, &[zero, inputs[0], zero, inputs[1]]);
    let expected = release(&b, p, &inputs);
    assert_eq!(mixed, expected);
    assert!(matches!(
        release_point_masses(&b, p, &[zero, zero], controls(1), || false),
        Err(AssemblyError::Budget {
            requested: 3,
            maximum: 2
        }),
    ));
    let mut c = controls(2);
    c.maximum_impulses = 1;
    assert!(matches!(
        release_point_masses(&b, p, &inputs, c, || false),
        Err(AssemblyError::Budget { .. })
    ));
    for invalid in [
        PointMassRelease {
            mass_kg: -1.0,
            ..inputs[0]
        },
        PointMassRelease {
            mass_kg: 0.0,
            ..inputs[0]
        },
        PointMassRelease {
            mass_kg: p.mass(),
            ..inputs[0]
        },
        PointMassRelease {
            mass_kg: f64::NAN,
            ..inputs[0]
        },
        PointMassRelease {
            mass_kg: f64::EPSILON / 8.0,
            ..inputs[0]
        },
        PointMassRelease {
            attachment_arm_body_m: Vec3::new(f64::NAN, 0.0, 0.0),
            ..zero
        },
        PointMassRelease {
            attachment_arm_body_m: Vec3::new(100.0, 0.0, 0.0),
            ..inputs[0]
        },
        PointMassRelease {
            impulse_on_parcel_body_ns: Vec3::new(f64::INFINITY, 0.0, 0.0),
            ..inputs[0]
        },
    ] {
        assert!(release_point_masses(&b, p, &[invalid], controls(1), || false).is_err());
    }
    assert!(
        release_point_masses(
            &b,
            p,
            &[
                PointMassRelease {
                    mass_kg: 2.0,
                    ..inputs[0]
                },
                PointMassRelease {
                    mass_kg: 2.0,
                    ..inputs[1]
                },
            ],
            controls(2),
            || false
        )
        .is_err()
    );
    let triangle = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(1.9, 1.0, 1.0)).unwrap();
    let resting = state(
        triangle,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Vec3::ZERO,
    );
    let parcel = PointMassRelease {
        mass_kg: 0.5,
        attachment_arm_body_m: Vec3::new(0.5, 0.0, 0.0),
        impulse_on_parcel_body_ns: Vec3::ZERO,
    };
    assert!(
        release_point_masses(&resting, triangle, &[parcel, parcel], controls(2), || false).is_err()
    );

    let mut all_polls = 0;
    let completed = release_point_masses(&b, p, &inputs, controls(2), || {
        all_polls += 1;
        false
    })
    .unwrap();
    assert!(all_polls >= 4 * inputs.len() + 2);
    let mut polls = 0;
    assert_eq!(
        release_point_masses(&b, p, &inputs, controls(2), || {
            polls += 1;
            polls == all_polls
        }),
        Err(AssemblyError::Cancelled)
    );
    assert_eq!(polls, all_polls);
    assert_eq!(release(&b, p, &inputs), completed);
    assert_eq!(b, balanced(1.0).1);
    assert_eq!(p, balanced(1.0).0);
}
