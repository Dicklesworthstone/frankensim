//! G1/G3/G4: internal wrench transport on original finite full-inertia bodies.

use fs_mbd::spatial_loads::internal_impulse::{
    InternalWrenchImpulse, InternalWrenchImpulseControls, InternalWrenchImpulseError,
    apply_internal_wrench_impulse,
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-10, "{a:e} != {b:e}");
}
fn vector(a: Vec3, b: Vec3) {
    close(a.x, b.x);
    close(a.y, b.y);
    close(a.z, b.z);
}
fn controls(origin: Vec3) -> InternalWrenchImpulseControls {
    InternalWrenchImpulseControls {
        angular_origin_world_m: origin,
        linear_tolerance_ns: 1e-10,
        angular_tolerance_nms: 1e-10,
        work_tolerance_j: 1e-10,
    }
}
fn body(position: Vec3, orientation: UnitQuaternion, p: Vec3, h: Vec3) -> RigidBodyState {
    RigidBodyState::new(Pose::new(position, orientation).unwrap(), p, h).unwrap()
}
fn properties() -> [MassProperties; 2] {
    [2.0, 3.0].map(|mass| {
        // I_body = [[2.5,-0.5,0],[-0.5,2.5,0],[0,0,4]], not diagonal.
        MassProperties::from_principal_frame(
            mass,
            Vec3::ZERO,
            Vec3::new(2.0, 3.0, 4.0),
            UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), core::f64::consts::FRAC_PI_4)
                .unwrap(),
        )
        .unwrap()
    })
}
fn states() -> [RigidBodyState; 2] {
    [
        body(
            Vec3::new(-1.0, 2.0, 0.5),
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.6).unwrap(),
            Vec3::new(0.5, -1.0, 2.0),
            Vec3::new(0.8, -0.4, 0.3),
        ),
        body(
            Vec3::new(2.0, -0.5, 1.5),
            UnitQuaternion::from_axis_angle(Vec3::new(-2.0, 1.0, 1.0), -0.4).unwrap(),
            Vec3::new(-0.7, 0.2, -0.5),
            Vec3::new(0.1, 0.5, -0.6),
        ),
    ]
}
fn wrench() -> InternalWrenchImpulse {
    InternalWrenchImpulse {
        arms_body_m: [Vec3::new(0.3, -0.2, 0.1), Vec3::new(-0.1, 0.4, 0.2)],
        impulse_on_second_world_ns: Vec3::new(0.5, -0.8, 0.3),
        couple_on_second_world_nms: Vec3::new(0.2, 0.1, -0.4),
    }
}
fn angular(states: [RigidBodyState; 2], origin: Vec3) -> Vec3 {
    states.iter().fold(Vec3::ZERO, |sum, state| {
        sum.add(
            state
                .pose()
                .position_world()
                .sub(origin)
                .cross(state.linear_momentum_world()),
        )
        .add(
            state
                .pose()
                .orientation()
                .rotate_body_to_world(state.angular_momentum_body()),
        )
    })
}
fn linear(states: [RigidBodyState; 2]) -> Vec3 {
    states[0]
        .linear_momentum_world()
        .add(states[1].linear_momentum_world())
}
fn kinetic(states: [RigidBodyState; 2]) -> f64 {
    // Independent inverse of the explicit full tensor above, plus original masses.
    states
        .iter()
        .zip([2.0, 3.0])
        .map(|(state, mass)| {
            let h = state.angular_momentum_body();
            state.linear_momentum_world().norm_squared() / (2.0 * mass)
                + 0.5 * ((2.5 * h.x * h.x + h.x * h.y + 2.5 * h.y * h.y) / 6.0 + h.z * h.z / 4.0)
        })
        .sum()
}

#[test]
fn g1_noncoincident_wrench_retains_linear_spin_orbital_momentum_and_actual_work() {
    let before = states();
    let props = properties();
    let load = wrench();
    let arms = core::array::from_fn::<_, 2, _>(|i| {
        before[i]
            .pose()
            .orientation()
            .rotate_body_to_world(load.arms_body_m[i])
    });
    let separation = before[1]
        .pose()
        .position_world()
        .add(arms[1])
        .sub(before[0].pose().position_world().add(arms[0]));
    let j = load.impulse_on_second_world_ns;
    let c = load.couple_on_second_world_nms;
    let opposite = c.add(separation.cross(j)).scale(-1.0);
    let angular_impulses = [
        arms[0].cross(j.scale(-1.0)).add(opposite),
        arms[1].cross(j).add(c),
    ];
    // Equal opposite point forces ALONE deliberately retain their orbital couple.
    let point_only = core::array::from_fn(|i| {
        before[i]
            .apply_impulse_at_body_point(
                props[i],
                load.arms_body_m[i],
                j.scale(if i == 0 { -1.0 } else { 1.0 }),
            )
            .unwrap()
            .state_after
    });
    vector(
        angular(point_only, Vec3::ZERO).sub(angular(before, Vec3::ZERO)),
        separation.cross(j),
    );
    assert!(separation.cross(j).norm_squared() > 0.1);
    for origin in [
        Vec3::ZERO,
        Vec3::new(8.0, -5.0, 3.0),
        Vec3::new(-40.0, 20.0, 10.0),
    ] {
        let result =
            apply_internal_wrench_impulse(before, props, load, controls(origin), || false).unwrap();
        let after = result.bodies.map(|r| r.state_after);
        vector(result.attachment_separation_world_m, separation);
        vector(result.attachment_couples_world_nms[0], opposite);
        vector(result.attachment_couples_world_nms[1], c);
        for i in 0..2 {
            assert_eq!(after[i].pose(), before[i].pose());
            vector(
                after[i]
                    .linear_momentum_world()
                    .sub(before[i].linear_momentum_world()),
                j.scale(if i == 0 { -1.0 } else { 1.0 }),
            );
            vector(
                after[i]
                    .angular_momentum_body()
                    .sub(before[i].angular_momentum_body()),
                before[i]
                    .pose()
                    .orientation()
                    .rotate_world_to_body(angular_impulses[i]),
            );
        }
        vector(linear(after), linear(before));
        vector(angular(after, origin), angular(before, origin));
        vector(result.linear_momentum_change_world_ns, Vec3::ZERO);
        vector(result.angular_momentum_change_world_nms, Vec3::ZERO);
        close(
            result.kinetic_energy_change_j,
            kinetic(after) - kinetic(before),
        );
        close(result.connector_work_j, kinetic(after) - kinetic(before));
        close(result.work_energy_residual_j, 0.0);
    }
}

#[test]
fn g1_coincident_point_limit_and_reversible_motor_energy_exchange() {
    let props = [MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap(); 2];
    let before = [-1.0, 1.0].map(|x| {
        body(
            Vec3::new(x, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
            Vec3::ZERO,
        )
    });
    let load = InternalWrenchImpulse {
        arms_body_m: [Vec3::new(1.0, 0.0, 0.0), Vec3::new(-1.0, 0.0, 0.0)],
        impulse_on_second_world_ns: Vec3::new(0.0, 0.5, 0.0),
        couple_on_second_world_nms: Vec3::ZERO,
    };
    let result =
        apply_internal_wrench_impulse(before, props, load, controls(Vec3::ZERO), || false).unwrap();
    assert_eq!(result.attachment_couples_world_nms, [Vec3::ZERO; 2]);
    for i in 0..2 {
        assert_eq!(
            result.bodies[i].state_after,
            before[i]
                .apply_impulse_at_body_point(
                    props[i],
                    load.arms_body_m[i],
                    load.impulse_on_second_world_ns
                        .scale(if i == 0 { -1.0 } else { 1.0 })
                )
                .unwrap()
                .state_after
        );
    }
    let motor = InternalWrenchImpulse {
        arms_body_m: [Vec3::ZERO; 2],
        impulse_on_second_world_ns: Vec3::ZERO,
        couple_on_second_world_nms: Vec3::new(0.0, 0.0, 2.0),
    };
    let driven =
        apply_internal_wrench_impulse(before, props, motor, controls(Vec3::ZERO), || false)
            .unwrap();
    close(driven.connector_work_j, 4.0);
    vector(
        driven.bodies[0].state_after.angular_momentum_body(),
        Vec3::new(0.0, 0.0, -2.0),
    );
    vector(
        driven.bodies[1].state_after.angular_momentum_body(),
        Vec3::new(0.0, 0.0, 2.0),
    );
    let reversed = apply_internal_wrench_impulse(
        driven.bodies.map(|r| r.state_after),
        props,
        InternalWrenchImpulse {
            couple_on_second_world_nms: motor.couple_on_second_world_nms.scale(-1.0),
            ..motor
        },
        controls(Vec3::ZERO),
        || false,
    )
    .unwrap();
    close(reversed.connector_work_j, -4.0);
    assert_eq!(reversed.bodies.map(|r| r.state_after), before);
}

fn multiply(a: UnitQuaternion, b: UnitQuaternion) -> UnitQuaternion {
    let [w, x, y, z] = a.components();
    let [v, p, q, r] = b.components();
    UnitQuaternion::new(
        w * v - x * p - y * q - z * r,
        w * p + x * v + y * r - z * q,
        w * q - x * r + y * v + z * p,
        w * r + x * q - y * p + z * v,
    )
    .unwrap()
}

#[test]
fn g3_world_frame_origin_and_common_velocity_covariance() {
    let before = states();
    let props = properties();
    let load = wrench();
    let original =
        apply_internal_wrench_impulse(before, props, load, controls(Vec3::ZERO), || false).unwrap();
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(2.0, -1.0, 3.0), 0.7).unwrap();
    let shift = Vec3::new(32.0, -16.0, 8.0);
    let boost = Vec3::new(16.0, -8.0, 4.0);
    let moved = core::array::from_fn(|i| {
        body(
            rotation
                .rotate_body_to_world(before[i].pose().position_world())
                .add(shift),
            multiply(rotation, before[i].pose().orientation()),
            rotation
                .rotate_body_to_world(before[i].linear_momentum_world())
                .add(boost.scale(props[i].mass())),
            before[i].angular_momentum_body(),
        )
    });
    let transformed = apply_internal_wrench_impulse(
        moved,
        props,
        InternalWrenchImpulse {
            impulse_on_second_world_ns: rotation
                .rotate_body_to_world(load.impulse_on_second_world_ns),
            couple_on_second_world_nms: rotation
                .rotate_body_to_world(load.couple_on_second_world_nms),
            ..load
        },
        controls(shift),
        || false,
    )
    .unwrap();
    for i in 0..2 {
        vector(
            transformed.bodies[i].state_after.linear_momentum_world(),
            rotation
                .rotate_body_to_world(original.bodies[i].state_after.linear_momentum_world())
                .add(boost.scale(props[i].mass())),
        );
        vector(
            transformed.bodies[i].state_after.angular_momentum_body(),
            original.bodies[i].state_after.angular_momentum_body(),
        );
    }
    vector(
        transformed.attachment_separation_world_m,
        rotation.rotate_body_to_world(original.attachment_separation_world_m),
    );
    close(transformed.connector_work_j, original.connector_work_j);
    close(
        transformed.kinetic_energy_change_j,
        original.kinetic_energy_change_j,
    );
}

#[test]
fn g3_small_resolved_transfer_survives_large_background_energy() {
    let props = [MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap(); 2];
    let boost = 2.0_f64.powi(40);
    let before = [0.0, 2.0].map(|x| {
        body(
            Vec3::new(x, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::new(boost, 0.0, 0.0),
            Vec3::ZERO,
        )
    });
    let load = InternalWrenchImpulse {
        arms_body_m: [Vec3::ZERO; 2],
        impulse_on_second_world_ns: Vec3::new(1.0, 0.0, 0.0),
        couple_on_second_world_nms: Vec3::ZERO,
    };
    let result = apply_internal_wrench_impulse(
        before,
        props,
        load,
        controls(Vec3::new(0.0, 8.0, 0.0)),
        || false,
    )
    .unwrap();
    assert_eq!(
        result.bodies[0].linear_increment_world,
        Vec3::new(-1.0, 0.0, 0.0)
    );
    assert_eq!(
        result.bodies[1].linear_increment_world,
        Vec3::new(1.0, 0.0, 0.0)
    );
    assert_eq!(result.linear_momentum_change_world_ns, Vec3::ZERO);
    assert_eq!(result.angular_momentum_change_world_nms, Vec3::ZERO);
    assert_eq!(result.kinetic_energy_change_j, 1.0);
    assert_eq!(result.connector_work_j, 1.0);
    // Rounded large total energies cannot resolve this physically real exchange.
    let rounded_change: f64 = result
        .bodies
        .iter()
        .map(|r| r.work.kinetic_energy_after_j - r.work.kinetic_energy_before_j)
        .sum();
    assert_eq!(rounded_change, 0.0);
}

#[test]
fn g4_swallowed_transfer_invalid_input_and_every_cancellation_phase_refuse() {
    let props = [MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap(); 2];
    let before = [0.0, 2.0].map(|x| {
        body(
            Vec3::new(x, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::new(2.0_f64.powi(54), 0.0, 0.0),
            Vec3::ZERO,
        )
    });
    let load = InternalWrenchImpulse {
        arms_body_m: [Vec3::ZERO; 2],
        impulse_on_second_world_ns: Vec3::new(1.0, 0.0, 0.0),
        couple_on_second_world_nms: Vec3::ZERO,
    };
    assert!(matches!(
        apply_internal_wrench_impulse(before, props, load, controls(Vec3::ZERO), || false),
        Err(InternalWrenchImpulseError::Unresolved {
            linear_residual_ns: 1.0,
            ..
        })
    ));
    // The reverse scale separation matters too: a large kick can erase old
    // momentum, and a rounded f64 endpoint subtraction conceals that loss.
    let old_linear = [1.0, -1.0].map(|p| {
        body(
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            Vec3::new(p, 0.0, 0.0),
            Vec3::ZERO,
        )
    });
    assert!(matches!(
        apply_internal_wrench_impulse(
            old_linear,
            props,
            InternalWrenchImpulse {
                impulse_on_second_world_ns: Vec3::new(2.0_f64.powi(54), 0.0, 0.0),
                ..load
            },
            controls(Vec3::ZERO),
            || false
        ),
        Err(InternalWrenchImpulseError::Unresolved {
            linear_residual_ns: 1.0,
            ..
        })
    ));
    let old_spin = [1.0, -1.0].map(|h| {
        body(
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
            Vec3::new(0.0, 0.0, h),
        )
    });
    assert!(matches!(
        apply_internal_wrench_impulse(
            old_spin,
            props,
            InternalWrenchImpulse {
                impulse_on_second_world_ns: Vec3::ZERO,
                couple_on_second_world_nms: Vec3::new(0.0, 0.0, 2.0_f64.powi(54)),
                ..load
            },
            controls(Vec3::ZERO),
            || false
        ),
        Err(InternalWrenchImpulseError::Unresolved {
            angular_residual_nms: 1.0,
            ..
        })
    ));
    // A unit free couple must not disappear behind a large attachment-force
    // moment before the per-body gate, even though the pair's net change is 0.
    for orientation in [
        UnitQuaternion::IDENTITY,
        UnitQuaternion::new(0.5, 0.5, 0.5, 0.5).unwrap(),
    ] {
        let large_arm = InternalWrenchImpulse {
            arms_body_m: [Vec3::new(2.0_f64.powi(54), 0.0, 0.0); 2],
            impulse_on_second_world_ns: orientation.rotate_body_to_world(Vec3::new(0.0, 1.0, 0.0)),
            couple_on_second_world_nms: orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0)),
        };
        let rest = [body(Vec3::ZERO, orientation, Vec3::ZERO, Vec3::ZERO); 2];
        assert!(matches!(
            apply_internal_wrench_impulse(rest, props, large_arm, controls(Vec3::ZERO), || false),
            Err(InternalWrenchImpulseError::Unresolved { angular_residual_nms, .. }) if angular_residual_nms >= 1.0
        ));
    }
    // Opposite frame-conversion losses can also cancel globally. Admission
    // compares EACH actual world increment to the original supplied couple.
    let orientation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.7).unwrap();
    let rest = [body(Vec3::ZERO, orientation, Vec3::ZERO, Vec3::ZERO); 2];
    let large_couple = InternalWrenchImpulse {
        arms_body_m: [Vec3::ZERO; 2],
        impulse_on_second_world_ns: Vec3::ZERO,
        couple_on_second_world_nms: Vec3::new(
            2.0_f64.powi(54),
            -2.0_f64.powi(53),
            2.0_f64.powi(52),
        ),
    };
    assert!(matches!(
        apply_internal_wrench_impulse(rest, props, large_couple,
            InternalWrenchImpulseControls { work_tolerance_j: 1e30, ..controls(Vec3::ZERO) }, || false),
        Err(InternalWrenchImpulseError::Unresolved { angular_residual_nms, .. }) if angular_residual_nms > 1e-3
    ));
    let ordinary = states();
    for cancel_at in 1..=6 {
        let mut calls = 0;
        let result = apply_internal_wrench_impulse(
            ordinary,
            properties(),
            wrench(),
            controls(Vec3::ZERO),
            || {
                calls += 1;
                calls == cancel_at
            },
        );
        assert_eq!(result, Err(InternalWrenchImpulseError::Cancelled));
        assert_eq!(ordinary, states());
    }
    let invalid = InternalWrenchImpulseControls {
        linear_tolerance_ns: f64::NAN,
        ..controls(Vec3::ZERO)
    };
    assert_eq!(
        apply_internal_wrench_impulse(ordinary, properties(), wrench(), invalid, || false),
        Err(InternalWrenchImpulseError::InvalidControls)
    );
    let invalid = InternalWrenchImpulse {
        arms_body_m: [Vec3::new(f64::NAN, 0.0, 0.0), Vec3::ZERO],
        ..wrench()
    };
    assert!(matches!(
        apply_internal_wrench_impulse(
            ordinary,
            properties(),
            invalid,
            controls(Vec3::ZERO),
            || false
        ),
        Err(InternalWrenchImpulseError::Dynamics(_))
    ));
    let retry = apply_internal_wrench_impulse(
        ordinary,
        properties(),
        wrench(),
        controls(Vec3::ZERO),
        || false,
    )
    .unwrap();
    assert_eq!(
        retry,
        apply_internal_wrench_impulse(
            ordinary,
            properties(),
            wrench(),
            controls(Vec3::ZERO),
            || false
        )
        .unwrap()
    );
}
