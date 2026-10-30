use fs_ga::Mat3;
use fs_mbd::rocket::rigid::tensor::{
    SpecificInertia, TensorRocketForcing, TensorRocketInertia, TensorRocketState,
};
use fs_mbd::rocket::rigid::{RigidRocketError, RigidRocketSolve, RocketInertiaModel};
use fs_mbd::rocket::{RocketForcing, RocketState};
use fs_mbd::{MassProperties, UnitQuaternion, Vec3};

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}; error {}",
        (a - b).abs()
    );
}
fn vector_close(a: Vec3, b: Vec3, tolerance: f64) {
    for (a, b) in [a.x, a.y, a.z].into_iter().zip([b.x, b.y, b.z]) {
        close(a, b, tolerance);
    }
}
fn max_norm(v: Vec3) -> f64 {
    v.x.abs().max(v.y.abs()).max(v.z.abs())
}
fn mul(a: Mat3, v: Vec3) -> Vec3 {
    Vec3::new(
        a.m[0] * v.x + a.m[1] * v.y + a.m[2] * v.z,
        a.m[3] * v.x + a.m[4] * v.y + a.m[5] * v.z,
        a.m[6] * v.x + a.m[7] * v.y + a.m[8] * v.z,
    )
}
fn add(a: Mat3, b: Mat3, scale: f64) -> Mat3 {
    Mat3 {
        m: core::array::from_fn(|i| a.m[i] + scale * b.m[i]),
    }
}
fn diagonal(x: f64, y: f64, z: f64) -> Mat3 {
    Mat3 {
        m: [x, 0.0, 0.0, 0.0, y, 0.0, 0.0, 0.0, z],
    }
}
fn point() -> RocketState {
    RocketState {
        position_m: Vec3::new(0.4, -0.2, 0.6),
        velocity_m_s: Vec3::new(0.8, -0.4, 0.2),
        dry_mass_kg: 6.0,
        propellant_mass_kg: 4.0,
    }
}
fn controls(h: f64) -> RigidRocketSolve {
    RigidRocketSolve {
        maximum_substep_s: h,
        maximum_substeps: 10000,
    }
}
fn input() -> TensorRocketForcing {
    TensorRocketForcing {
        mass_flow_kg_s: 0.75,
        effective_exhaust_velocity_body_m_s: Vec3::new(0.4, -0.2, 1.3),
        thrust_application_arm_body_m: Vec3::new(0.1, -0.2, 0.05),
        co_rotating_exhaust_specific_inertia_m2: SpecificInertia::ZERO,
        exhaust_swirl_specific_angular_momentum_body_m2_s: Vec3::new(0.01, -0.02, 0.03),
        gravity_world_m_s2: Vec3::new(0.0, 0.0, -1.2),
        external_force_world_n: Vec3::new(0.2, -0.1, 0.3),
        external_torque_body_n_m: Vec3::new(0.05, 0.03, -0.02),
    }
}
fn fixture() -> (TensorRocketState, TensorRocketForcing) {
    let dry = Mat3 {
        m: [3.2, 0.3, -0.15, 0.3, 4.1, 0.2, -0.15, 0.2, 5.0],
    };
    let prop = SpecificInertia::from_principal_frame(
        Vec3::new(0.25, 0.4, 0.55),
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 0.5), 0.7).unwrap(),
    )
    .unwrap();
    let model = TensorRocketInertia::new(dry, prop).unwrap();
    let state = TensorRocketState::new(
        point(),
        UnitQuaternion::from_axis_angle(Vec3::new(-0.2, 0.3, 1.0), 0.4).unwrap(),
        Vec3::new(1.1, -0.7, 1.8),
        model,
    )
    .unwrap();
    let forcing = TensorRocketForcing {
        co_rotating_exhaust_specific_inertia_m2: SpecificInertia::from_principal_frame(
            Vec3::new(0.1, 0.2, 0.25),
            UnitQuaternion::from_axis_angle(Vec3::new(0.2, 1.0, -0.7), -0.8).unwrap(),
        )
        .unwrap(),
        ..input()
    };
    (state, forcing)
}

#[test]
fn g1_oblique_corotating_depletion_retains_rate_and_exact_burnout_translation() {
    let axes = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -0.7), 0.65).unwrap();
    let dry = MassProperties::from_principal_frame(6.0, Vec3::ZERO, Vec3::new(2.0, 3.0, 4.0), axes)
        .unwrap()
        .inertia_tensor_body();
    let prop = SpecificInertia::from_principal_frame(Vec3::new(0.3, 0.4, 0.5), axes).unwrap();
    let model = TensorRocketInertia::new(dry, prop).unwrap();
    let omega = axes.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.2));
    let properties = model.mass_properties(point()).unwrap();
    let body = TensorRocketState::new(
        point(),
        UnitQuaternion::IDENTITY,
        properties.angular_momentum_body_checked(omega).unwrap(),
        model,
    )
    .unwrap();
    let forcing = TensorRocketForcing {
        mass_flow_kg_s: 2.0,
        effective_exhaust_velocity_body_m_s: omega.scale(10.0),
        thrust_application_arm_body_m: Vec3::ZERO,
        co_rotating_exhaust_specific_inertia_m2: prop,
        exhaust_swirl_specific_angular_momentum_body_m2_s: Vec3::ZERO,
        external_torque_body_n_m: Vec3::ZERO,
        ..input()
    };
    let result = body
        .advance(forcing, 3.0, controls(0.125), || false)
        .unwrap();
    let after = result.state_after;
    vector_close(
        after
            .mass_properties()
            .unwrap()
            .angular_velocity_body_checked(after.angular_momentum_body_kg_m2_s())
            .unwrap(),
        omega,
        3e-13,
    );
    vector_close(
        result.exported_angular_momentum_world_kg_m2_s,
        mul(prop.tensor_body_m2(), omega).scale(4.0),
        5e-13,
    );
    vector_close(
        result.angular_momentum_residual_world_kg_m2_s,
        Vec3::ZERO,
        5e-13,
    );
    assert_eq!(result.burned_propellant_mass_kg, 4.0);
    assert_eq!(result.powered_duration_s, 2.0);
    assert_eq!(result.coast_duration_s, 1.0);
    assert_eq!(after.translation().propellant_mass_kg, 0.0);
    let translation = point()
        .advance(
            RocketForcing {
                mass_flow_kg_s: 2.0,
                effective_exhaust_velocity_m_s: omega.scale(10.0),
                gravity_m_s2: forcing.gravity_world_m_s2,
                external_force_n: forcing.external_force_world_n,
            },
            3.0,
            || false,
        )
        .unwrap();
    vector_close(
        after.translation().position_m,
        translation.state_after.position_m,
        2e-13,
    );
    vector_close(
        after.translation().velocity_m_s,
        translation.state_after.velocity_m_s,
        2e-13,
    );
}

#[test]
fn g1_noncommuting_outflow_is_simultaneous_and_tracks_all_tensor_components() {
    let (body, forcing) = fixture();
    let model = body.inertia();
    let before = body.mass_properties().unwrap();
    let after = model
        .mass_properties(RocketState {
            propellant_mass_kg: 2.5,
            ..point()
        })
        .unwrap();
    let i0 = before.inertia_tensor_body();
    let i1 = after.inertia_tensor_body();
    let e = forcing
        .co_rotating_exhaust_specific_inertia_m2
        .tensor_body_m2();
    let commutator = add(i0.compose(e), e.compose(i0), -1.0);
    assert!(
        commutator.max_abs() > 0.01,
        "fixture must distinguish a scalar/diagonal treatment"
    );
    let midpoint_omega = Vec3::new(0.4, -0.9, 1.2);
    let mid_h = mul(add(i0, i1, 1.0), midpoint_omega).scale(0.5);
    let export = mul(e, midpoint_omega).scale(1.5);
    let h0 = mid_h.add(export.scale(0.5));
    let (h1, actual_export) = model
        .exhaust_spin_step(
            h0,
            before,
            after,
            forcing.co_rotating_exhaust_specific_inertia_m2,
            1.5,
        )
        .unwrap();
    vector_close(h1, mid_h.sub(export.scale(0.5)), 3e-14);
    vector_close(actual_export, export, 2e-15);

    // If removed mass carries precisely its own co-rotation inertia, the
    // frozen-attitude submap retains arbitrary omega even for moving axes.
    let (h1, _) = model
        .exhaust_spin_step(
            mul(i0, midpoint_omega),
            before,
            after,
            model.propellant_specific_inertia(),
            1.5,
        )
        .unwrap();
    vector_close(h1, mul(i1, midpoint_omega), 2e-14);
}

// Independent continuous ODE reference: cofactor inverse and raw quaternion
// differential equation, without native stepping, principal axes or kick maps.
fn inverse_apply(a: Mat3, v: Vec3) -> Vec3 {
    let m = a.m;
    let c00 = m[4] * m[8] - m[5] * m[7];
    let c01 = m[2] * m[7] - m[1] * m[8];
    let c02 = m[1] * m[5] - m[2] * m[4];
    let c10 = m[5] * m[6] - m[3] * m[8];
    let c11 = m[0] * m[8] - m[2] * m[6];
    let c12 = m[2] * m[3] - m[0] * m[5];
    let c20 = m[3] * m[7] - m[4] * m[6];
    let c21 = m[1] * m[6] - m[0] * m[7];
    let c22 = m[0] * m[4] - m[1] * m[3];
    let d = m[0] * c00 + m[1] * c10 + m[2] * c20;
    Vec3::new(
        c00 * v.x + c01 * v.y + c02 * v.z,
        c10 * v.x + c11 * v.y + c12 * v.z,
        c20 * v.x + c21 * v.y + c22 * v.z,
    )
    .scale(1.0 / d)
}
fn derivative(
    t: f64,
    y: [f64; 13],
    body: TensorRocketState,
    input: TensorRocketForcing,
) -> [f64; 13] {
    let model = body.inertia();
    let propellant = body.translation().propellant_mass_kg - input.mass_flow_kg_s * t;
    let inertia = add(
        model.dry_inertia_body_kg_m2(),
        model.propellant_specific_inertia().tensor_body_m2(),
        propellant,
    );
    let h = Vec3::new(y[4], y[5], y[6]);
    let omega = inverse_apply(inertia, h);
    let thrust = input
        .effective_exhaust_velocity_body_m_s
        .scale(input.mass_flow_kg_s);
    let torque = input
        .external_torque_body_n_m
        .add(input.thrust_application_arm_body_m.cross(thrust));
    let flux = mul(
        input
            .co_rotating_exhaust_specific_inertia_m2
            .tensor_body_m2(),
        omega,
    )
    .add(input.exhaust_swirl_specific_angular_momentum_body_m2_s)
    .scale(input.mass_flow_kg_s);
    let dh = torque.sub(omega.cross(h)).sub(flux);
    let qv = Vec3::new(y[1], y[2], y[3]);
    let dqv = omega.scale(y[0]).add(qv.cross(omega)).scale(0.5);
    // Quaternion vector rotation without the native orientation algebra.
    let world_thrust = thrust
        .add(qv.cross(thrust).scale(2.0 * y[0]))
        .add(qv.cross(qv.cross(thrust)).scale(2.0));
    let accel = world_thrust
        .add(input.external_force_world_n)
        .scale(1.0 / (body.translation().dry_mass_kg + propellant))
        .add(input.gravity_world_m_s2);
    [
        -0.5 * qv.dot(omega),
        dqv.x,
        dqv.y,
        dqv.z,
        dh.x,
        dh.y,
        dh.z,
        y[10],
        y[11],
        y[12],
        accel.x,
        accel.y,
        accel.z,
    ]
}
fn continuous_reference(
    body: TensorRocketState,
    input: TensorRocketForcing,
    duration: f64,
    steps: usize,
) -> [f64; 13] {
    let q = body.orientation().components();
    let h = body.angular_momentum_body_kg_m2_s();
    let p = body.translation().position_m;
    let v = body.translation().velocity_m_s;
    let mut y = [
        q[0], q[1], q[2], q[3], h.x, h.y, h.z, p.x, p.y, p.z, v.x, v.y, v.z,
    ];
    let dt = duration / steps as f64;
    for j in 0..steps {
        let t = j as f64 * dt;
        let k1 = derivative(t, y, body, input);
        let k2 = derivative(
            t + 0.5 * dt,
            core::array::from_fn(|i| y[i] + 0.5 * dt * k1[i]),
            body,
            input,
        );
        let k3 = derivative(
            t + 0.5 * dt,
            core::array::from_fn(|i| y[i] + 0.5 * dt * k2[i]),
            body,
            input,
        );
        let k4 = derivative(
            t + dt,
            core::array::from_fn(|i| y[i] + dt * k3[i]),
            body,
            input,
        );
        y = core::array::from_fn(|i| y[i] + dt / 6.0 * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]));
    }
    y
}

#[test]
fn g1_full_tensor_flight_converges_to_independent_variable_mass_ode() {
    let (body, input) = fixture();
    let reference = continuous_reference(body, input, 1.0, 4096);
    let check = continuous_reference(body, input, 1.0, 2048);
    for (a, b) in reference.into_iter().zip(check) {
        close(a, b, 2e-13);
    }
    let mut errors = Vec::new();
    for steps in [8, 16, 32] {
        let result = body
            .advance(input, 1.0, controls(1.0 / steps as f64), || false)
            .unwrap();
        let after = result.state_after;
        let q = after.orientation().components();
        let h = after.angular_momentum_body_kg_m2_s();
        let p = after.translation().position_m;
        let v = after.translation().velocity_m_s;
        let actual = [
            q[0], q[1], q[2], q[3], h.x, h.y, h.z, p.x, p.y, p.z, v.x, v.y, v.z,
        ];
        errors.push(
            actual
                .into_iter()
                .zip(reference)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max),
        );
        vector_close(
            result.angular_momentum_residual_world_kg_m2_s,
            Vec3::ZERO,
            8e-13,
        );
    }
    println!("full tensor rocket maximum state errors: {errors:?}");
    assert!(errors[2] < 2e-5);
    for pair in errors.windows(2) {
        assert!(pair[0] / pair[1] > 3.7 && pair[0] / pair[1] < 4.3);
    }
    let wet = body.mass_properties().unwrap().principal_axes_to_body();
    let dry = body
        .inertia()
        .mass_properties(RocketState {
            propellant_mass_kg: 0.0,
            ..point()
        })
        .unwrap()
        .principal_axes_to_body();
    assert!(
        wet.components()
            .into_iter()
            .zip(dry.components())
            .map(|(a, b)| (a - b).abs())
            .sum::<f64>()
            > 0.01
    );
}

fn rotation(q: UnitQuaternion) -> Mat3 {
    let x = q.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
    let y = q.rotate_body_to_world(Vec3::new(0.0, 1.0, 0.0));
    let z = q.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
    Mat3 {
        m: [x.x, y.x, z.x, x.y, y.y, z.y, x.z, y.z, z.z],
    }
}
fn quaternion_product(a: UnitQuaternion, b: UnitQuaternion) -> UnitQuaternion {
    let [w, x, y, z] = a.components();
    let [s, i, j, k] = b.components();
    UnitQuaternion::new(
        w * s - x * i - y * j - z * k,
        w * i + x * s + y * k - z * j,
        w * j - x * k + y * s + z * i,
        w * k + x * j - y * i + z * s,
    )
    .unwrap()
}

#[test]
fn g3_body_coordinate_rotation_world_rotation_and_common_boost_preserve_motion() {
    let (body, input) = fixture();
    let duration = 0.75;
    let baseline = body
        .advance(input, duration, controls(0.025), || false)
        .unwrap();
    let s = UnitQuaternion::from_axis_angle(Vec3::new(0.4, -0.3, 1.0), 0.6).unwrap();
    let r = rotation(s);
    let model = TensorRocketInertia::new(
        r.compose(body.inertia().dry_inertia_body_kg_m2())
            .compose(r.transpose()),
        SpecificInertia::from_inertia_tensor(
            r.compose(
                body.inertia()
                    .propellant_specific_inertia()
                    .tensor_body_m2(),
            )
            .compose(r.transpose()),
        )
        .unwrap(),
    )
    .unwrap();
    let [w, x, y, z] = s.components();
    let inverse_s = UnitQuaternion::new(w, -x, -y, -z).unwrap();
    let rebased = TensorRocketState::new(
        *body.translation(),
        quaternion_product(body.orientation(), inverse_s),
        s.rotate_body_to_world(body.angular_momentum_body_kg_m2_s()),
        model,
    )
    .unwrap();
    let rebased_input = TensorRocketForcing {
        effective_exhaust_velocity_body_m_s: s
            .rotate_body_to_world(input.effective_exhaust_velocity_body_m_s),
        thrust_application_arm_body_m: s.rotate_body_to_world(input.thrust_application_arm_body_m),
        co_rotating_exhaust_specific_inertia_m2: SpecificInertia::from_inertia_tensor(
            r.compose(
                input
                    .co_rotating_exhaust_specific_inertia_m2
                    .tensor_body_m2(),
            )
            .compose(r.transpose()),
        )
        .unwrap(),
        exhaust_swirl_specific_angular_momentum_body_m2_s: s
            .rotate_body_to_world(input.exhaust_swirl_specific_angular_momentum_body_m2_s),
        external_torque_body_n_m: s.rotate_body_to_world(input.external_torque_body_n_m),
        ..input
    };
    let result = rebased
        .advance(rebased_input, duration, controls(0.025), || false)
        .unwrap();
    vector_close(
        result.state_after.translation().position_m,
        baseline.state_after.translation().position_m,
        3e-13,
    );
    vector_close(
        result.state_after.translation().velocity_m_s,
        baseline.state_after.translation().velocity_m_s,
        3e-13,
    );
    vector_close(
        result.state_after.angular_momentum_body_kg_m2_s(),
        s.rotate_body_to_world(baseline.state_after.angular_momentum_body_kg_m2_s()),
        2e-12,
    );
    vector_close(
        result.exported_angular_momentum_world_kg_m2_s,
        baseline.exported_angular_momentum_world_kg_m2_s,
        2e-13,
    );
    vector_close(
        result
            .state_after
            .orientation()
            .rotate_body_to_world(s.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0))),
        baseline
            .state_after
            .orientation()
            .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
        3e-13,
    );

    let boost = Vec3::new(8.0, -2.0, 3.0);
    let offset = Vec3::new(10.0, -4.0, 1.0);
    let world = TensorRocketState::new(
        RocketState {
            position_m: s
                .rotate_body_to_world(body.translation().position_m)
                .add(offset),
            velocity_m_s: s
                .rotate_body_to_world(body.translation().velocity_m_s)
                .add(boost),
            ..*body.translation()
        },
        quaternion_product(s, body.orientation()),
        body.angular_momentum_body_kg_m2_s(),
        body.inertia(),
    )
    .unwrap();
    let world_input = TensorRocketForcing {
        gravity_world_m_s2: s.rotate_body_to_world(input.gravity_world_m_s2),
        external_force_world_n: s.rotate_body_to_world(input.external_force_world_n),
        ..input
    };
    let result = world
        .advance(world_input, duration, controls(0.025), || false)
        .unwrap();
    vector_close(
        result.state_after.translation().position_m,
        s.rotate_body_to_world(baseline.state_after.translation().position_m)
            .add(offset)
            .add(boost.scale(duration)),
        5e-13,
    );
    vector_close(
        result.state_after.translation().velocity_m_s,
        s.rotate_body_to_world(baseline.state_after.translation().velocity_m_s)
            .add(boost),
        5e-13,
    );
    vector_close(
        result.state_after.angular_momentum_body_kg_m2_s(),
        baseline.state_after.angular_momentum_body_kg_m2_s(),
        3e-13,
    );
    vector_close(
        result.exported_angular_momentum_world_kg_m2_s,
        s.rotate_body_to_world(baseline.exported_angular_momentum_world_kg_m2_s),
        3e-13,
    );
}

#[test]
fn g0_zero_spin_export_preserves_world_momentum_while_inertia_axes_move() {
    let (body, _) = fixture();
    let forcing = TensorRocketForcing {
        mass_flow_kg_s: 1.0,
        effective_exhaust_velocity_body_m_s: Vec3::ZERO,
        exhaust_swirl_specific_angular_momentum_body_m2_s: Vec3::ZERO,
        external_torque_body_n_m: Vec3::ZERO,
        co_rotating_exhaust_specific_inertia_m2: SpecificInertia::ZERO,
        ..input()
    };
    let result = body
        .advance(forcing, 4.5, controls(0.025), || false)
        .unwrap();
    assert_eq!(result.exported_angular_momentum_world_kg_m2_s, Vec3::ZERO);
    vector_close(
        result
            .state_after
            .orientation()
            .rotate_body_to_world(result.state_after.angular_momentum_body_kg_m2_s()),
        body.orientation()
            .rotate_body_to_world(body.angular_momentum_body_kg_m2_s()),
        2e-12,
    );
    assert_eq!(result.state_after.translation().propellant_mass_kg, 0.0);
    let first = body
        .advance(forcing, 1.0, controls(0.025), || false)
        .unwrap();
    let second = first
        .state_after
        .advance(forcing, 3.5, controls(0.025), || false)
        .unwrap();
    vector_close(
        second.state_after.angular_momentum_body_kg_m2_s(),
        result.state_after.angular_momentum_body_kg_m2_s(),
        1e-12,
    );
}

#[test]
fn g4_invalid_models_semidefinite_limits_and_late_failure_remain_atomic() {
    let identity = UnitQuaternion::IDENTITY;
    assert!(SpecificInertia::from_principal_frame(Vec3::new(-0.1, 1.0, 1.0), identity).is_err());
    assert!(SpecificInertia::from_principal_frame(Vec3::new(3.0, 1.0, 1.0), identity).is_err());
    assert!(
        SpecificInertia::from_principal_frame(Vec3::new(f64::NAN, 1.0, 1.0), identity).is_err()
    );
    assert!(
        SpecificInertia::from_inertia_tensor(Mat3 {
            m: [1.0, 0.2, 0.0, 0.1, 1.0, 0.0, 0.0, 0.0, 1.0]
        })
        .is_err()
    );
    assert!(TensorRocketInertia::new(diagonal(4.0, 1.0, 1.0), SpecificInertia::ZERO).is_err());
    assert!(TensorRocketInertia::new(diagonal(0.0, 1.0, 1.0), SpecificInertia::ZERO).is_err());
    let line = SpecificInertia::from_principal_frame(
        Vec3::new(0.0, 0.4, 0.4),
        UnitQuaternion::from_axis_angle(Vec3::new(0.3, -1.0, 0.1), 0.7).unwrap(),
    )
    .unwrap();
    let model = TensorRocketInertia::new(diagonal(2.0, 2.0, 2.0), line).unwrap();
    assert!(model.mass_properties(point()).is_ok());
    let (body, forcing) = fixture();
    let original = body;
    let baseline = body
        .advance(forcing, 0.5, controls(0.025), || false)
        .unwrap();
    let mut polls = 0;
    assert!(matches!(
        body.advance(forcing, 0.5, controls(0.025), || {
            polls += 1;
            polls == 45
        }),
        Err(RigidRocketError::Cancelled)
            | Err(RigidRocketError::Translation(
                fs_mbd::rocket::RocketError::Cancelled
            ))
    ));
    assert!(polls >= 45);
    assert!(matches!(
        body.advance(
            forcing,
            0.5,
            RigidRocketSolve {
                maximum_substep_s: 0.025,
                maximum_substeps: 19
            },
            || false
        ),
        Err(RigidRocketError::SubstepBudget {
            completed: 19,
            maximum: 19
        })
    ));
    assert_eq!(body, original);
    assert_eq!(
        body.advance(forcing, 0.5, controls(0.025), || false)
            .unwrap(),
        baseline
    );
    assert_eq!(
        body.advance(forcing, 0.0, controls(1.0), || false)
            .unwrap()
            .state_after,
        body
    );
    assert!(max_norm(baseline.angular_momentum_residual_world_kg_m2_s) < 1e-12);
}
