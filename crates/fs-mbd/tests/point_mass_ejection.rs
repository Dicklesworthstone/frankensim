use fs_mbd::assembly::release::{
    PointMassEjection, PointMassEjectionBatchResult, PointMassEjectionControls,
    PointMassEjectionError, PointMassRelease, eject_point_masses, release_point_masses,
};
use fs_mbd::assembly::{AssemblyControls, AssemblyError};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

fn controls(count: usize) -> PointMassEjectionControls {
    PointMassEjectionControls {
        assembly: AssemblyControls {
            maximum_parts: count + 1,
            maximum_impulses: count,
            linear_momentum_tolerance_ns: 1e-9,
            angular_momentum_tolerance_nms: 1e-8,
            energy_tolerance_j: 1e-8,
        },
        relative_velocity_tolerance_m_s: 1e-10,
        maximum_supplied_work_j: 1e6,
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
fn body(p: MassProperties, x: Vec3, q: UnitQuaternion, v: Vec3, w: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(x, q).unwrap(),
        v.scale(p.mass()),
        p.angular_momentum_body_checked(w).unwrap(),
    )
    .unwrap()
}
fn resting(p: MassProperties) -> RigidBodyState {
    body(
        p,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Vec3::ZERO,
    )
}
fn ejection(mass: f64, arm: Vec3, velocity: Vec3) -> PointMassEjection {
    PointMassEjection {
        mass_kg: mass,
        attachment_arm_body_m: arm,
        relative_velocity_body_m_s: velocity,
    }
}
fn eject(
    b: &RigidBodyState,
    p: MassProperties,
    inputs: &[PointMassEjection],
) -> PointMassEjectionBatchResult {
    eject_point_masses(b, p, inputs, controls(inputs.len()), || false).unwrap()
}
// Independent native point-kinematics evaluation at the returned physical point.
fn velocities(inputs: &[PointMassEjection], out: &PointMassEjectionBatchResult, tolerance: f64) {
    let release = &out.release;
    let q = release.state.pose().orientation();
    for (input, parcel) in inputs
        .iter()
        .filter(|input| input.mass_kg > 0.0)
        .zip(&release.parcels)
    {
        let arm = q.rotate_world_to_body(
            parcel
                .position_world_m
                .sub(release.state.pose().position_world()),
        );
        let rigid = release
            .state
            .point_kinematics(release.properties, arm)
            .unwrap();
        vector(
            parcel
                .velocity_world_m_s()
                .unwrap()
                .sub(rigid.point_velocity_world),
            q.rotate_body_to_world(input.relative_velocity_body_m_s),
            tolerance,
        );
    }
}
// Direct energy and spin-plus-orbit sums about an unrelated fixed observer.
fn conservation(b: RigidBodyState, p: MassProperties, out: &PointMassEjectionBatchResult) {
    let observer = Vec3::new(-0.7, 1.1, 0.4);
    let angular = |state: RigidBodyState| {
        state
            .pose()
            .orientation()
            .rotate_body_to_world(state.angular_momentum_body())
            .add(
                state
                    .pose()
                    .position_world()
                    .sub(observer)
                    .cross(state.linear_momentum_world()),
            )
    };
    let r = &out.release;
    let mut mass = r.properties.mass();
    let mut linear = r.state.linear_momentum_world();
    let mut spin_orbit = angular(r.state);
    let mut energy = r.state.kinetic_energy(r.properties).unwrap();
    for parcel in &r.parcels {
        mass += parcel.mass_kg;
        linear = linear.add(parcel.linear_momentum_world_ns);
        spin_orbit = spin_orbit.add(
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
    vector(linear, b.linear_momentum_world(), 3e-12);
    vector(spin_orbit, angular(b), 2e-10);
    close(
        energy,
        b.kinetic_energy(p).unwrap() + r.supplied_work_j,
        3e-10,
    );
}
fn asymmetric() -> MassProperties {
    // Retained M=3, centre=(-1/3,-2/3,0), I=diag(4,5,6), and
    // m1=1 at ex plus m2=2 at ey. Original COM is zero.
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

#[test]
fn central_jet_uses_reduced_mass_and_retains_nearly_depleted_body_recoil() {
    let p = MassProperties::new(4.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap();
    let b = body(
        p,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::new(3.0, -1.0, 2.0),
        Vec3::new(0.2, -0.3, 0.4),
    );
    let inputs = [ejection(1.0, Vec3::ZERO, Vec3::new(4.0, -2.0, 1.0))];
    let out = eject(&b, p, &inputs);
    // mu=m*M/(m+M)=3/4, J=mu*c, W=mu*|c|^2/2.
    vector(
        out.impulses[0].impulse_on_parcel_body_ns,
        Vec3::new(3.0, -1.5, 0.75),
        2e-15,
    );
    close(out.release.supplied_work_j, 63.0 / 8.0, 8e-15);
    assert_eq!(
        out.release.state.angular_momentum_body(),
        b.angular_momentum_body()
    );
    velocities(&inputs, &out, 2e-14);
    conservation(b, p, &out);

    let total = 10_000_000_000_000_000.0;
    let p = MassProperties::new(total, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap();
    let inputs = [ejection(total - 2.0, Vec3::ZERO, Vec3::new(1.0, 0.0, 0.0))];
    let out = eject(&resting(p), p, &inputs);
    let reduced = (total - 2.0) * (2.0 / total);
    close(out.release.properties.mass(), 2.0, 0.0);
    close(out.impulses[0].impulse_on_parcel_body_ns.x, reduced, 1e-15);
    assert!(((total - 2.0) * (1.0 - (total - 2.0) / total) - reduced).abs() > 0.2);
    close(out.release.supplied_work_j, 0.5 * reduced, 2e-15);
    velocities(&inputs, &out, 5e-15);
}

#[test]
fn simultaneous_jets_share_retained_translation_and_rotation() {
    let p = MassProperties::new(4.0, Vec3::ZERO, Vec3::new(3.0, 6.0, 7.0)).unwrap();
    for original_spin in [0.0, 1.0] {
        let b = body(
            p,
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            Vec3::new(3.0, -2.0, 1.0),
            Vec3::new(0.0, 0.0, original_spin),
        );
        let opposed = [
            ejection(
                1.0,
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 7.0 / 5.0, 0.0),
            ),
            ejection(
                1.0,
                Vec3::new(-1.0, 0.0, 0.0),
                Vec3::new(0.0, -7.0 / 5.0, 0.0),
            ),
        ];
        let out = eject(&b, p, &opposed);
        vector(
            out.impulses[0].impulse_on_parcel_body_ns,
            Vec3::new(0.0, 1.0, 0.0),
            8e-15,
        );
        vector(
            out.impulses[1].impulse_on_parcel_body_ns,
            Vec3::new(0.0, -1.0, 0.0),
            8e-15,
        );
        close(out.release.supplied_work_j, 7.0 / 5.0, 2e-14);
        vector(
            out.release
                .properties
                .angular_velocity_body_checked(out.release.state.angular_momentum_body())
                .unwrap(),
            Vec3::new(0.0, 0.0, original_spin - 2.0 / 5.0),
            8e-15,
        );
        velocities(&opposed, &out, 2e-14);
        conservation(b, p, &out);

        let parallel = opposed.map(|input| PointMassEjection {
            relative_velocity_body_m_s: Vec3::new(0.0, 2.0, 0.0),
            ..input
        });
        let out = eject(&b, p, &parallel);
        for impulse in &out.impulses {
            vector(
                impulse.impulse_on_parcel_body_ns,
                Vec3::new(0.0, 1.0, 0.0),
                1e-14,
            );
        }
        close(out.release.supplied_work_j, 2.0, 2e-14);
        velocities(&parallel, &out, 2e-14);
        conservation(b, p, &out);
    }
}

#[test]
fn full_tensor_zero_relative_jet_still_receives_shared_recoil_impulse() {
    let p = asymmetric();
    let b = resting(p);
    let inputs = [
        ejection(1.0, Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 3.0)),
        ejection(2.0, Vec3::new(0.0, 1.0, 0.0), Vec3::ZERO),
    ];
    let out = eject(&b, p, &inputs);
    // The original full xy inertia has determinant 46. Conservation gives
    // dw=(1/23,11/23,0), J1=93/46 ez, J2=-21/23 ez.
    vector(
        out.impulses[0].impulse_on_parcel_body_ns,
        Vec3::new(0.0, 0.0, 93.0 / 46.0),
        3e-13,
    );
    vector(
        out.impulses[1].impulse_on_parcel_body_ns,
        Vec3::new(0.0, 0.0, -21.0 / 23.0),
        3e-13,
    );
    vector(
        out.release
            .properties
            .angular_velocity_body_checked(out.release.state.angular_momentum_body())
            .unwrap(),
        Vec3::new(1.0 / 23.0, 11.0 / 23.0, 0.0),
        3e-13,
    );
    vector(
        out.release.center_shift_body_m,
        Vec3::new(-1.0 / 3.0, -2.0 / 3.0, 0.0),
        2e-15,
    );
    // Independent impulse/relative-speed work: sum(J_i.c_i)/2.
    close(out.release.supplied_work_j, 279.0 / 92.0, 8e-13);
    velocities(&inputs, &out, 5e-13);
    conservation(b, p, &out);
    let mut c = controls(2);
    c.assembly.maximum_impulses = 1;
    assert!(matches!(
        eject_point_masses(&b, p, &inputs, c, || false),
        Err(PointMassEjectionError::Assembly(AssemblyError::Budget {
            requested: 2,
            maximum: 1
        }))
    ));
}

#[test]
fn covariance_subdivision_and_permutation_keep_the_coupled_event() {
    let p = asymmetric();
    let b = body(
        p,
        Vec3::new(2.0, -3.0, 1.0),
        UnitQuaternion::IDENTITY,
        Vec3::new(-0.5, 0.7, 0.8),
        Vec3::new(0.2, -0.3, 0.4),
    );
    let inputs = [
        ejection(1.0, Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.2, -0.1, 0.5)),
        ejection(2.0, Vec3::new(0.0, 1.0, 0.0), Vec3::new(-0.4, 0.7, -0.6)),
    ];
    let plain = eject(&b, p, &inputs);
    velocities(&inputs, &plain, 3e-13);
    conservation(b, p, &plain);
    let reversed = eject(&b, p, &[inputs[1], inputs[0]]);
    vector(
        reversed.release.state.linear_momentum_world(),
        plain.release.state.linear_momentum_world(),
        3e-14,
    );
    vector(
        reversed.release.state.angular_momentum_body(),
        plain.release.state.angular_momentum_body(),
        3e-13,
    );
    close(
        reversed.release.supplied_work_j,
        plain.release.supplied_work_j,
        3e-13,
    );
    vector(
        reversed.impulses[0].impulse_on_parcel_body_ns,
        plain.impulses[1].impulse_on_parcel_body_ns,
        3e-14,
    );

    let split = [
        PointMassEjection {
            mass_kg: 0.25,
            ..inputs[0]
        },
        inputs[1],
        PointMassEjection {
            mass_kg: 0.75,
            ..inputs[0]
        },
    ];
    let divided = eject(&b, p, &split);
    vector(
        divided.release.state.linear_momentum_world(),
        plain.release.state.linear_momentum_world(),
        3e-14,
    );
    vector(
        divided.release.state.angular_momentum_body(),
        plain.release.state.angular_momentum_body(),
        5e-13,
    );
    vector(
        divided.impulses[0]
            .impulse_on_parcel_body_ns
            .add(divided.impulses[2].impulse_on_parcel_body_ns),
        plain.impulses[0].impulse_on_parcel_body_ns,
        3e-14,
    );
    close(
        divided.release.supplied_work_j,
        plain.release.supplied_work_j,
        5e-13,
    );
    velocities(&split, &divided, 4e-13);

    let q = UnitQuaternion::from_axis_angle(Vec3::new(0.3, -0.4, 0.8), 0.6).unwrap();
    let shift = Vec3::new(4.0, -2.0, 3.0);
    let boost = Vec3::new(64.0, -32.0, 16.0);
    let transformed = RigidBodyState::new(
        Pose::new(
            q.rotate_body_to_world(b.pose().position_world()).add(shift),
            q,
        )
        .unwrap(),
        q.rotate_body_to_world(b.linear_momentum_world())
            .add(boost.scale(p.mass())),
        b.angular_momentum_body(),
    )
    .unwrap();
    let moved = eject(&transformed, p, &inputs);
    assert_eq!(moved.impulses, plain.impulses);
    assert_eq!(moved.release.properties, plain.release.properties);
    assert_eq!(
        moved.release.state.angular_momentum_body(),
        plain.release.state.angular_momentum_body()
    );
    assert_eq!(moved.release.supplied_work_j, plain.release.supplied_work_j);
    vector(
        moved.release.state.pose().position_world(),
        q.rotate_body_to_world(plain.release.state.pose().position_world())
            .add(shift),
        4e-14,
    );
    vector(
        moved.release.state.linear_momentum_world(),
        q.rotate_body_to_world(plain.release.state.linear_momentum_world())
            .add(boost.scale(plain.release.properties.mass())),
        3e-13,
    );
    velocities(&inputs, &moved, 5e-13);
    conservation(transformed, p, &moved);
}

#[test]
fn common_boost_keeps_recoil_and_erased_actual_velocity_fails_its_own_gate() {
    let p = MassProperties::new(4.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap();
    let inputs = [ejection(1.0, Vec3::ZERO, Vec3::new(4.0, 0.0, 0.0))];
    let plain = eject(&resting(p), p, &inputs);
    let moved = body(
        p,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::new(1_099_511_627_776.0, 0.0, 0.0),
        Vec3::ZERO,
    );
    let boosted = eject(&moved, p, &inputs);
    assert_eq!(boosted.impulses, plain.impulses);
    assert_eq!(
        boosted.release.supplied_work_j,
        plain.release.supplied_work_j
    );
    assert_eq!(boosted.maximum_relative_velocity_residual_m_s, 0.0);
    velocities(&inputs, &boosted, 0.0);

    // Permit the native energy defect here ONLY to isolate the independent
    // velocity boundary-condition gate. Every canonical recoil increment is
    // erased by rounding at this much larger common velocity.
    let too_fast = body(
        p,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::new(1_152_921_504_606_846_976.0, 0.0, 0.0),
        Vec3::ZERO,
    );
    let mut c = controls(1);
    c.assembly.energy_tolerance_j = 10.0;
    let error = eject_point_masses(&too_fast, p, &inputs, c, || false).unwrap_err();
    let PointMassEjectionError::RelativeVelocity {
        input_index,
        residual_body_m_s,
        ..
    } = error
    else {
        panic!("expected actual velocity refusal, got {error:?}");
    };
    assert_eq!(input_index, 0);
    vector(residual_body_m_s, Vec3::new(-4.0, 0.0, 0.0), 0.0);
}

#[test]
fn passive_parity_source_budget_descriptor_limits_and_late_cancel_are_atomic() {
    let p = asymmetric();
    let b = body(
        p,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::new(1.0, -2.0, 3.0),
        Vec3::new(0.2, -0.3, 0.4),
    );
    let inputs = [
        ejection(1.0, Vec3::new(1.0, 0.0, 0.0), Vec3::ZERO),
        ejection(2.0, Vec3::new(0.0, 1.0, 0.0), Vec3::ZERO),
    ];
    let native: Vec<_> = inputs
        .iter()
        .map(|input| PointMassRelease {
            mass_kg: input.mass_kg,
            attachment_arm_body_m: input.attachment_arm_body_m,
            impulse_on_parcel_body_ns: Vec3::ZERO,
        })
        .collect();
    let mut c = controls(2);
    c.maximum_supplied_work_j = 0.0;
    c.assembly.maximum_impulses = 0;
    let passive = eject_point_masses(&b, p, &inputs, c, || false).unwrap();
    assert_eq!(
        passive.release,
        release_point_masses(&b, p, &native, c.assembly, || false).unwrap()
    );
    velocities(&inputs, &passive, 4e-13);
    assert_eq!(eject(&b, p, &[]).release.state, b);
    let zero = ejection(0.0, Vec3::ZERO, Vec3::new(1e300, 0.0, 0.0));
    let identity = eject_point_masses(&b, p, &[zero, zero], c, || false).unwrap();
    assert_eq!(identity.release.state, b);
    assert!(identity.release.parcels.is_empty());
    assert!(matches!(
        eject_point_masses(&b, p, &[zero, zero], controls(1), || false),
        Err(PointMassEjectionError::Assembly(
            AssemblyError::Budget { .. }
        ))
    ));

    let p = MassProperties::new(4.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0)).unwrap();
    let b = resting(p);
    let active = [ejection(1.0, Vec3::ZERO, Vec3::new(4.0, 0.0, 0.0))];
    let mut c = controls(1);
    c.maximum_supplied_work_j = 5.0;
    let error = eject_point_masses(&b, p, &active, c, || false).unwrap_err();
    let PointMassEjectionError::SourceWork {
        required_j,
        maximum_j,
    } = error
    else {
        panic!("expected source-work refusal, got {error:?}");
    };
    close(required_j, 6.0, 3e-14);
    assert_eq!(maximum_j, 5.0);
    for invalid in [
        PointMassEjection {
            mass_kg: -1.0,
            ..active[0]
        },
        PointMassEjection {
            mass_kg: p.mass(),
            ..active[0]
        },
        PointMassEjection {
            attachment_arm_body_m: Vec3::new(100.0, 0.0, 0.0),
            ..active[0]
        },
        PointMassEjection {
            relative_velocity_body_m_s: Vec3::new(f64::NAN, 0.0, 0.0),
            ..zero
        },
    ] {
        assert!(eject_point_masses(&b, p, &[invalid], controls(1), || false).is_err());
    }
    c = controls(1);
    c.relative_velocity_tolerance_m_s = f64::INFINITY;
    assert!(eject_point_masses(&b, p, &active, c, || false).is_err());
    c = controls(1);
    c.maximum_supplied_work_j = -1.0;
    assert!(eject_point_masses(&b, p, &active, c, || false).is_err());

    let mut total_polls = 0;
    let expected = eject_point_masses(&b, p, &active, controls(1), || {
        total_polls += 1;
        false
    })
    .unwrap();
    let mut polls = 0;
    assert_eq!(
        eject_point_masses(&b, p, &active, controls(1), || {
            polls += 1;
            polls == total_polls
        }),
        Err(PointMassEjectionError::Assembly(AssemblyError::Cancelled))
    );
    assert_eq!(polls, total_polls);
    assert_eq!(eject(&b, p, &active), expected);
    assert_eq!(b, resting(p));
}
