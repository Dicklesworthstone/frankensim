use fs_ga::Mat3;
use fs_mbd::rocket::RocketState;
use fs_mbd::rocket::rigid::staging::RocketStage;
use fs_mbd::rocket::rigid::tensor::{
    SpecificInertia, TensorRocketForcing, TensorRocketInertia, TensorRocketState,
};
use fs_mbd::rocket::rigid::{RigidRocketError, RigidRocketSolve};
use fs_mbd::{UnitQuaternion, Vec3};

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}
fn near_vec(a: Vec3, b: Vec3, tolerance: f64) {
    for (a, b) in [a.x, a.y, a.z].into_iter().zip([b.x, b.y, b.z]) {
        near(a, b, tolerance);
    }
}
fn add(a: Mat3, b: Mat3) -> Mat3 {
    Mat3 {
        m: core::array::from_fn(|i| a.m[i] + b.m[i]),
    }
}
fn scale(a: Mat3, factor: f64) -> Mat3 {
    Mat3 {
        m: a.m.map(|x| factor * x),
    }
}
fn mul(a: Mat3, v: Vec3) -> Vec3 {
    Vec3::new(
        a.m[0] * v.x + a.m[1] * v.y + a.m[2] * v.z,
        a.m[3] * v.x + a.m[4] * v.y + a.m[5] * v.z,
        a.m[6] * v.x + a.m[7] * v.y + a.m[8] * v.z,
    )
}
// Independent reciprocal-basis solve, without native principal-axis admission.
fn inverse_mul(a: Mat3, v: Vec3) -> Vec3 {
    let x = Vec3::new(a.m[0], a.m[3], a.m[6]);
    let y = Vec3::new(a.m[1], a.m[4], a.m[7]);
    let z = Vec3::new(a.m[2], a.m[5], a.m[8]);
    Vec3::new(v.dot(y.cross(z)), v.dot(z.cross(x)), v.dot(x.cross(y)))
        .scale(x.dot(y.cross(z)).recip())
}
fn shifted(inertia: Mat3, mass: f64, r: Vec3) -> Mat3 {
    let v = [r.x, r.y, r.z];
    Mat3 {
        m: core::array::from_fn(|i| {
            let row = i / 3;
            let col = i % 3;
            inertia.m[i] + mass * (if row == col { r.dot(r) } else { 0.0 } - v[row] * v[col])
        }),
    }
}

fn fixture() -> (
    TensorRocketState,
    [RocketStage<TensorRocketInertia>; 2],
    Vec3,
) {
    let r = Vec3::new(-0.4, 0.3, -0.2);
    let dry = [
        Mat3 {
            m: [1.1, 0.08, -0.05, 0.08, 1.4, 0.06, -0.05, 0.06, 1.7],
        },
        Mat3 {
            m: [0.8, -0.06, 0.04, -0.06, 1.0, 0.03, 0.04, 0.03, 1.2],
        },
    ];
    let propellant = [
        SpecificInertia::from_principal_frame(
            Vec3::new(0.12, 0.18, 0.24),
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -0.5), 0.4).unwrap(),
        )
        .unwrap(),
        SpecificInertia::from_principal_frame(
            Vec3::new(0.09, 0.14, 0.19),
            UnitQuaternion::from_axis_angle(Vec3::new(-0.4, 1.0, 0.7), -0.55).unwrap(),
        )
        .unwrap(),
    ];
    let stages = [
        RocketStage {
            dry_mass_kg: 3.0,
            propellant_mass_kg: 3.0,
            center_of_mass_body_m: r,
            inertia: TensorRocketInertia::new(dry[0], propellant[0]).unwrap(),
        },
        RocketStage {
            dry_mass_kg: 2.0,
            propellant_mass_kg: 2.0,
            center_of_mass_body_m: r.scale(-1.5),
            inertia: TensorRocketInertia::new(dry[1], propellant[1]).unwrap(),
        },
    ];
    // Dry and propellant COMs both coincide with the parent COM. Each tensor
    // assembly is independently formed in the common material coordinates.
    let dry_parent = add(shifted(dry[0], 3.0, r), shifted(dry[1], 2.0, r.scale(-1.5)));
    let prop_parent = scale(
        add(
            shifted(scale(propellant[0].tensor_body_m2(), 3.0), 3.0, r),
            shifted(
                scale(propellant[1].tensor_body_m2(), 2.0),
                2.0,
                r.scale(-1.5),
            ),
        ),
        0.2,
    );
    let model = TensorRocketInertia::new(
        dry_parent,
        SpecificInertia::from_inertia_tensor(prop_parent).unwrap(),
    )
    .unwrap();
    let omega = Vec3::new(0.2, -0.35, 0.45);
    let parent = TensorRocketState::new(
        RocketState {
            position_m: Vec3::new(3.0, -4.0, 2.0),
            velocity_m_s: Vec3::new(0.7, -0.2, 0.4),
            dry_mass_kg: 5.0,
            propellant_mass_kg: 5.0,
        },
        UnitQuaternion::from_axis_angle(Vec3::new(0.3, -0.2, 1.0), 0.6).unwrap(),
        mul(add(dry_parent, scale(prop_parent, 5.0)), omega),
        model,
    )
    .unwrap();
    (parent, stages, omega)
}

// G1: complete spin, orbital motion and recoil follow from the two original
// child tensors, including products of inertia and nonprincipal angular speed.
#[test]
fn full_tensor_staging_preserves_spatial_momentum_and_physical_recoil_work() {
    let (parent, stages, omega) = fixture();
    let impulse = Vec3::new(0.6, -0.4, 0.3);
    let point = Vec3::new(0.1, 0.05, -0.2);
    let result = parent
        .separate(stages[0], stages[1], impulse, point, 1e-12, || false)
        .unwrap();
    let rotation = parent.orientation();
    let mut linear = Vec3::ZERO;
    let mut angular = Vec3::ZERO;
    let mut recoil_work = 0.0;
    for (i, child) in [result.retained, result.released].into_iter().enumerate() {
        let properties = child.mass_properties().unwrap();
        let inertia = add(
            stages[i].inertia.dry_inertia_body_kg_m2(),
            scale(
                stages[i]
                    .inertia
                    .propellant_specific_inertia()
                    .tensor_body_m2(),
                stages[i].propellant_mass_kg,
            ),
        );
        let r = stages[i].center_of_mass_body_m;
        let j = impulse.scale(if i == 0 { 1.0 } else { -1.0 });
        let spin_kick = point.sub(r).cross(j);
        let inherited_velocity = parent
            .translation()
            .velocity_m_s
            .add(rotation.rotate_body_to_world(omega.cross(r)));
        let expected_p = inherited_velocity
            .scale(properties.mass())
            .add(rotation.rotate_body_to_world(j));
        let expected_h = mul(inertia, omega).add(spin_kick);
        let body = child.rigid_body().unwrap();
        near_vec(
            child.translation().position_m,
            parent
                .translation()
                .position_m
                .add(rotation.rotate_body_to_world(r)),
            3e-14,
        );
        near_vec(body.linear_momentum_world(), expected_p, 3e-12);
        near_vec(child.angular_momentum_body_kg_m2_s(), expected_h, 3e-12);
        linear = linear.add(body.linear_momentum_world());
        angular = angular
            .add(rotation.rotate_body_to_world(child.angular_momentum_body_kg_m2_s()))
            .add(
                rotation
                    .rotate_body_to_world(r)
                    .cross(body.linear_momentum_world()),
            );
        recoil_work +=
            0.5 * (j.dot(j) / properties.mass() + spin_kick.dot(inverse_mul(inertia, spin_kick)));
        assert_eq!(child.inertia(), stages[i].inertia);
        assert_eq!(
            child.translation().propellant_mass_kg,
            stages[i].propellant_mass_kg
        );
    }
    near_vec(
        linear,
        parent.rigid_body().unwrap().linear_momentum_world(),
        4e-12,
    );
    near_vec(
        angular,
        rotation.rotate_body_to_world(parent.angular_momentum_body_kg_m2_s()),
        4e-12,
    );
    near(result.separation_work_j, recoil_work, 3e-12);
    near(result.energy_balance_residual_j, 0.0, 3e-12);
    near_vec(
        result.linear_momentum_residual_world_kg_m_s,
        Vec3::ZERO,
        4e-12,
    );
    near_vec(
        result.angular_momentum_residual_world_kg_m2_s,
        Vec3::ZERO,
        4e-12,
    );
    assert_eq!(result.mass_residual_kg, 0.0);
}

// G1: no-recoil release inherits the whole rigid velocity field. Both child
// inventories then continue through the existing full-tensor propulsion owner.
#[test]
fn released_tensor_stages_keep_their_own_inertia_and_continue_powered_flight() {
    let (parent, stages, omega) = fixture();
    let result = parent
        .separate(stages[0], stages[1], Vec3::ZERO, Vec3::ZERO, 1e-12, || {
            false
        })
        .unwrap();
    near(result.separation_work_j, 0.0, 1e-14);
    for (i, child) in [result.retained, result.released].into_iter().enumerate() {
        near_vec(
            child
                .mass_properties()
                .unwrap()
                .angular_velocity_body_checked(child.angular_momentum_body_kg_m2_s())
                .unwrap(),
            omega,
            3e-12,
        );
        near_vec(
            child.translation().velocity_m_s,
            parent.translation().velocity_m_s.add(
                parent
                    .orientation()
                    .rotate_body_to_world(omega.cross(stages[i].center_of_mass_body_m)),
            ),
            3e-12,
        );
        let forcing = TensorRocketForcing {
            mass_flow_kg_s: 0.4,
            effective_exhaust_velocity_body_m_s: Vec3::new(2.0, 0.1, -0.2),
            thrust_application_arm_body_m: Vec3::new(0.0, 0.1, 0.05),
            co_rotating_exhaust_specific_inertia_m2: child.inertia().propellant_specific_inertia(),
            exhaust_swirl_specific_angular_momentum_body_m2_s: Vec3::ZERO,
            gravity_world_m_s2: Vec3::ZERO,
            external_force_world_n: Vec3::ZERO,
            external_torque_body_n_m: Vec3::ZERO,
        };
        let flight = child
            .advance(
                forcing,
                0.5,
                RigidRocketSolve {
                    maximum_substep_s: 0.025,
                    maximum_substeps: 100,
                },
                || false,
            )
            .unwrap();
        near(flight.burned_propellant_mass_kg, 0.2, 2e-14);
        near(
            flight.state_after.translation().propellant_mass_kg,
            stages[i].propellant_mass_kg - 0.2,
            2e-14,
        );
        near_vec(
            flight.angular_momentum_residual_world_kg_m2_s,
            Vec3::ZERO,
            2e-10,
        );
        assert_eq!(flight.state_after.inertia(), child.inertia());
        assert!(flight.thrust_impulse_world_n_s.norm_squared() > 0.01);
    }
}

// G3/G4: diagonal-only admission would miss this valid off-diagonal change.
// Neither failed admission nor a late cancellation changes the source owner.
#[test]
fn off_diagonal_mismatch_and_late_cancellation_refuse_atomically() {
    let (parent, stages, _) = fixture();
    let mut wrong = stages[0];
    let mut dry = wrong.inertia.dry_inertia_body_kg_m2();
    dry.m[1] += 0.01;
    dry.m[3] += 0.01;
    wrong.inertia =
        TensorRocketInertia::new(dry, wrong.inertia.propellant_specific_inertia()).unwrap();
    assert!(matches!(
        parent.separate(wrong, stages[1], Vec3::ZERO, Vec3::ZERO, 1e-12, || false),
        Err(RigidRocketError::Invalid(
            "stage inertia must reproduce all parent parallel-axis entries"
        ))
    ));
    let mut checkpoints = 0;
    let accepted = parent
        .separate(stages[0], stages[1], Vec3::ZERO, Vec3::ZERO, 1e-12, || {
            checkpoints += 1;
            false
        })
        .unwrap();
    let original = parent;
    let mut calls = 0;
    assert_eq!(
        parent.separate(stages[0], stages[1], Vec3::ZERO, Vec3::ZERO, 1e-12, || {
            calls += 1;
            calls == checkpoints
        }),
        Err(RigidRocketError::Cancelled)
    );
    assert_eq!(parent, original);
    assert_eq!(
        parent
            .separate(stages[0], stages[1], Vec3::ZERO, Vec3::ZERO, 1e-12, || {
                false
            })
            .unwrap(),
        accepted
    );
}
