//! G1/G3/G4: finite internally driven motion with independent trajectory oracles.

use fs_mbd::spatial_loads::internal_impulse::{
    InternalWrenchImpulseControls, InternalWrenchImpulseError,
};
use fs_mbd::spatial_loads::internal_motion::{
    InternalWrenchStepError, PrescribedInternalWrench, step_internal_wrench,
};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() < tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    assert!(a.sub(b).norm_squared().sqrt() < tolerance, "{a:?} != {b:?}");
}
fn controls(origin: Vec3) -> InternalWrenchImpulseControls {
    InternalWrenchImpulseControls {
        angular_origin_world_m: origin,
        linear_tolerance_ns: 2e-10,
        angular_tolerance_nms: 2e-10,
        work_tolerance_j: 2e-10,
    }
}
fn state(x: Vec3, q: UnitQuaternion, p: Vec3, h: Vec3) -> RigidBodyState {
    RigidBodyState::new(Pose::new(x, q).unwrap(), p, h).unwrap()
}
fn properties() -> [MassProperties; 2] {
    core::array::from_fn(|i| {
        MassProperties::from_principal_frame(
            [2.0, 3.0][i],
            Vec3::ZERO,
            [Vec3::new(2.0, 3.0, 4.0), Vec3::new(3.0, 4.0, 5.0)][i],
            UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), core::f64::consts::FRAC_PI_4)
                .unwrap(),
        )
        .unwrap()
    })
}
fn initial() -> [RigidBodyState; 2] {
    [
        state(
            Vec3::new(-1.0, 0.2, -0.3),
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.4).unwrap(),
            Vec3::new(0.4, -0.1, 0.2),
            Vec3::new(0.4, 0.2, -0.3),
        ),
        state(
            Vec3::new(1.4, -0.5, 0.7),
            UnitQuaternion::from_axis_angle(Vec3::new(-2.0, 1.0, 0.5), -0.6).unwrap(),
            Vec3::new(-0.3, 0.7, -0.4),
            Vec3::new(-0.2, 0.5, 0.1),
        ),
    ]
}
fn drive() -> PrescribedInternalWrench {
    PrescribedInternalWrench {
        arms_body_m: [Vec3::new(0.3, -0.1, 0.2), Vec3::new(-0.2, 0.4, 0.1)],
        force_on_second_world_n: Vec3::new(0.9, -0.7, 0.5),
        couple_on_second_world_nm: Vec3::new(0.1, 0.2, -0.3),
    }
}
fn no_drive() -> PrescribedInternalWrench {
    PrescribedInternalWrench {
        arms_body_m: [Vec3::ZERO; 2],
        force_on_second_world_n: Vec3::ZERO,
        couple_on_second_world_nm: Vec3::ZERO,
    }
}
fn linear(bodies: [RigidBodyState; 2]) -> Vec3 {
    bodies[0]
        .linear_momentum_world()
        .add(bodies[1].linear_momentum_world())
}
fn angular(bodies: [RigidBodyState; 2], origin: Vec3) -> Vec3 {
    bodies.iter().fold(Vec3::ZERO, |sum, b| {
        sum.add(
            b.pose()
                .position_world()
                .sub(origin)
                .cross(b.linear_momentum_world()),
        )
        .add(
            b.pose()
                .orientation()
                .rotate_body_to_world(b.angular_momentum_body()),
        )
    })
}
// Independent inverses of the explicit nonprincipal body tensors:
// A=[[2.5,-.5,0],[-.5,2.5,0],[0,0,4]], B=[[3.5,-.5,0],[-.5,3.5,0],[0,0,5]].
fn omega(h: Vec3, i: usize) -> Vec3 {
    let a = [2.5, 3.5][i];
    let determinant = [6.0, 12.0][i];
    Vec3::new(
        (a * h.x + 0.5 * h.y) / determinant,
        (0.5 * h.x + a * h.y) / determinant,
        h.z / [4.0, 5.0][i],
    )
}
fn kinetic(bodies: [RigidBodyState; 2]) -> f64 {
    (0..2)
        .map(|i| {
            let b = bodies[i];
            b.linear_momentum_world().norm_squared() / (2.0 * [2.0, 3.0][i])
                + 0.5
                    * b.angular_momentum_body()
                        .dot(omega(b.angular_momentum_body(), i))
        })
        .sum()
}

#[test]
fn g1_internal_force_advances_both_coms_with_exact_unequal_mass_recoil() {
    let props = properties();
    let before = [
        state(
            Vec3::new(-1.0, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::ZERO,
        ),
        state(
            Vec3::new(1.0, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::new(-0.3, 0.0, 0.0),
            Vec3::ZERO,
        ),
    ];
    let load = PrescribedInternalWrench {
        force_on_second_world_n: Vec3::new(3.0, 0.0, 0.0),
        ..no_drive()
    };
    let h = 0.4;
    let result = step_internal_wrench(
        before,
        props,
        load,
        h,
        controls(Vec3::new(2.0, -3.0, 4.0)),
        || false,
    )
    .unwrap();
    for i in 0..2 {
        let force = load
            .force_on_second_world_n
            .scale(if i == 0 { -1.0 } else { 1.0 });
        let displacement = before[i]
            .linear_momentum_world()
            .scale(h / props[i].mass())
            .add(force.scale(0.5 * h * h / props[i].mass()));
        vector(
            result.states_after[i].pose().position_world(),
            before[i].pose().position_world().add(displacement),
            1e-12,
        );
        vector(
            result.states_after[i].linear_momentum_world(),
            before[i].linear_momentum_world().add(force.scale(h)),
            1e-12,
        );
        vector(
            result.states_after[i].angular_momentum_body(),
            Vec3::ZERO,
            1e-12,
        );
    }
    let relative_displacement = result.states_after[1]
        .pose()
        .position_world()
        .sub(before[1].pose().position_world())
        .sub(
            result.states_after[0]
                .pose()
                .position_world()
                .sub(before[0].pose().position_world()),
        );
    close(
        result.actuator_work_j,
        load.force_on_second_world_n.dot(relative_displacement),
        1e-12,
    );
    close(
        result.kinetic_energy_change_j,
        kinetic(result.states_after) - kinetic(before),
        1e-12,
    );
    close(result.work_energy_residual_j, 0.0, 1e-12);
    assert_ne!(result.states_after[0].pose(), before[0].pose());
    assert_ne!(
        result.half_impulses[0].attachment_separation_world_m,
        result.half_impulses[1].attachment_separation_world_m
    );
}

#[test]
fn g1_finite_motor_couples_exchange_spin_angle_and_recover_energy() {
    let props = properties();
    let before = [
        state(
            Vec3::new(-1.0, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
            Vec3::ZERO,
        ),
        state(
            Vec3::new(1.0, 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
            Vec3::ZERO,
        ),
    ];
    let load = PrescribedInternalWrench {
        couple_on_second_world_nm: Vec3::new(0.0, 0.0, 0.8),
        ..no_drive()
    };
    let h = 0.6;
    let driven =
        step_internal_wrench(before, props, load, h, controls(Vec3::ZERO), || false).unwrap();
    let mut relative_angle = 0.0;
    for i in 0..2 {
        let sign = if i == 0 { -1.0 } else { 1.0 };
        let angle = sign * 0.5 * 0.8 * h * h / [4.0, 5.0][i];
        vector(
            driven.states_after[i]
                .pose()
                .orientation()
                .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
            Vec3::new(angle.cos(), angle.sin(), 0.0),
            2e-12,
        );
        vector(
            driven.states_after[i].angular_momentum_body(),
            Vec3::new(0.0, 0.0, sign * 0.8 * h),
            1e-12,
        );
        relative_angle += sign * angle;
    }
    close(driven.actuator_work_j, 0.8 * relative_angle, 1e-12);
    close(driven.actuator_work_j, kinetic(driven.states_after), 1e-12);
    vector(angular(driven.states_after, Vec3::ZERO), Vec3::ZERO, 1e-12);
    let braking = PrescribedInternalWrench {
        couple_on_second_world_nm: load.couple_on_second_world_nm.scale(-1.0),
        ..load
    };
    let stopped = step_internal_wrench(
        driven.states_after,
        props,
        braking,
        h,
        controls(Vec3::ZERO),
        || false,
    )
    .unwrap();
    assert!(stopped.actuator_work_j < 0.0);
    close(stopped.actuator_work_j, -driven.actuator_work_j, 1e-12);
    for i in 0..2 {
        vector(
            stopped.states_after[i].angular_momentum_body(),
            Vec3::ZERO,
            1e-12,
        );
        assert_ne!(
            stopped.states_after[i].pose().orientation(),
            UnitQuaternion::IDENTITY
        );
    }
}

// Independent continuous ODE reference. No native time step or impulse update
// is called: translations are analytic, and RK4 integrates quaternion kinematics,
// Euler's full-tensor equations and instantaneous mechanical actuator power.
#[derive(Clone, Copy)]
struct Reference {
    q: [[f64; 4]; 2],
    h: [Vec3; 2],
    work: f64,
}
fn qmul(a: [f64; 4], b: [f64; 4]) -> [f64; 4] {
    let [w, x, y, z] = a;
    let [v, p, q, r] = b;
    [
        w * v - x * p - y * q - z * r,
        w * p + x * v + y * r - z * q,
        w * q - x * r + y * v + z * p,
        w * r + x * q - y * p + z * v,
    ]
}
fn rotate(q: [f64; 4], v: Vec3) -> Vec3 {
    let [w, x, y, z] = q;
    Vec3::new(
        (1.0 - 2.0 * (y * y + z * z)) * v.x
            + 2.0 * (x * y - w * z) * v.y
            + 2.0 * (x * z + w * y) * v.z,
        2.0 * (x * y + w * z) * v.x
            + (1.0 - 2.0 * (x * x + z * z)) * v.y
            + 2.0 * (y * z - w * x) * v.z,
        2.0 * (x * z - w * y) * v.x
            + 2.0 * (y * z + w * x) * v.y
            + (1.0 - 2.0 * (x * x + y * y)) * v.z,
    )
}
fn conjugate([w, x, y, z]: [f64; 4]) -> [f64; 4] {
    [w, -x, -y, -z]
}
fn shifted(s: Reference, k: Reference, scale: f64) -> Reference {
    Reference {
        q: core::array::from_fn(|i| core::array::from_fn(|j| s.q[i][j] + scale * k.q[i][j])),
        h: core::array::from_fn(|i| s.h[i].add(k.h[i].scale(scale))),
        work: s.work + scale * k.work,
    }
}
fn rhs(
    s: Reference,
    t: f64,
    initial: [RigidBodyState; 2],
    load: PrescribedInternalWrench,
) -> Reference {
    let velocity: [Vec3; 2] = core::array::from_fn(|i| {
        let f = load
            .force_on_second_world_n
            .scale(if i == 0 { -1.0 } else { 1.0 });
        initial[i]
            .linear_momentum_world()
            .add(f.scale(t))
            .scale(1.0 / [2.0, 3.0][i])
    });
    let x: [Vec3; 2] = core::array::from_fn(|i| {
        let f = load
            .force_on_second_world_n
            .scale(if i == 0 { -1.0 } else { 1.0 });
        initial[i]
            .pose()
            .position_world()
            .add(initial[i].linear_momentum_world().scale(t / [2.0, 3.0][i]))
            .add(f.scale(0.5 * t * t / [2.0, 3.0][i]))
    });
    let rb = rotate(s.q[1], load.arms_body_m[1]);
    let torque = [
        x[1].sub(x[0])
            .add(rb)
            .cross(load.force_on_second_world_n)
            .add(load.couple_on_second_world_nm)
            .scale(-1.0),
        rb.cross(load.force_on_second_world_n)
            .add(load.couple_on_second_world_nm),
    ];
    let w: [Vec3; 2] = core::array::from_fn(|i| omega(s.h[i], i));
    Reference {
        q: core::array::from_fn(|i| qmul(s.q[i], [0.0, w[i].x, w[i].y, w[i].z]).map(|a| 0.5 * a)),
        h: core::array::from_fn(|i| s.h[i].cross(w[i]).add(rotate(conjugate(s.q[i]), torque[i]))),
        work: load
            .force_on_second_world_n
            .dot(velocity[1].sub(velocity[0]))
            + torque[0].dot(rotate(s.q[0], w[0]))
            + torque[1].dot(rotate(s.q[1], w[1])),
    }
}
fn reference(
    initial: [RigidBodyState; 2],
    load: PrescribedInternalWrench,
    time: f64,
    steps: usize,
) -> Reference {
    let mut s = Reference {
        q: initial.map(|b| b.pose().orientation().components()),
        h: initial.map(|b| b.angular_momentum_body()),
        work: 0.0,
    };
    let h = time / steps as f64;
    for j in 0..steps {
        let t = j as f64 * h;
        let a = rhs(s, t, initial, load);
        let b = rhs(shifted(s, a, 0.5 * h), t + 0.5 * h, initial, load);
        let c = rhs(shifted(s, b, 0.5 * h), t + 0.5 * h, initial, load);
        let d = rhs(shifted(s, c, h), t + h, initial, load);
        s = shifted(
            shifted(shifted(shifted(s, a, h / 6.0), b, h / 3.0), c, h / 3.0),
            d,
            h / 6.0,
        );
        for q in &mut s.q {
            let norm = q.iter().map(|x| x * x).sum::<f64>().sqrt();
            for x in q {
                *x /= norm;
            }
        }
    }
    s
}
fn qerror(a: [f64; 4], b: [f64; 4]) -> f64 {
    let minus = (0..4).map(|i| (a[i] - b[i]).powi(2)).sum::<f64>().sqrt();
    let plus = (0..4).map(|i| (a[i] + b[i]).powi(2)).sum::<f64>().sqrt();
    minus.min(plus)
}
fn integrate(
    initial: [RigidBodyState; 2],
    load: PrescribedInternalWrench,
    duration: f64,
    steps: usize,
) -> ([RigidBodyState; 2], f64) {
    let mut state = initial;
    let mut work = 0.0;
    for _ in 0..steps {
        let next = step_internal_wrench(
            state,
            properties(),
            load,
            duration / steps as f64,
            controls(Vec3::new(3.0, -2.0, 1.0)),
            || false,
        )
        .unwrap();
        vector(next.linear_momentum_residual_ns, Vec3::ZERO, 2e-10);
        vector(next.angular_momentum_residual_nms, Vec3::ZERO, 2e-10);
        close(next.work_energy_residual_j, 0.0, 2e-10);
        close(next.native_drift_energy_defect_j, 0.0, 2e-10);
        work += next.actuator_work_j;
        state = next.states_after;
    }
    (state, work)
}

#[test]
fn g1_full_inertia_offcentre_driven_trajectory_refines_quadratically_to_independent_ode() {
    let start = initial();
    let load = drive();
    let duration = 1.0;
    let truth = reference(start, load, duration, 8192);
    let check = reference(start, load, duration, 4096);
    for i in 0..2 {
        close(qerror(truth.q[i], check.q[i]), 0.0, 1e-11);
        vector(truth.h[i], check.h[i], 1e-11);
    }
    close(truth.work, check.work, 1e-11);
    let mut errors = Vec::new();
    for steps in [8, 16, 32] {
        let (end, work) = integrate(start, load, duration, steps);
        vector(linear(end), linear(start), 2e-11);
        for origin in [Vec3::ZERO, Vec3::new(8.0, -5.0, 3.0)] {
            vector(angular(end, origin), angular(start, origin), 2e-10);
        }
        close(work, kinetic(end) - kinetic(start), 2e-10);
        let angle = (0..2)
            .map(|i| qerror(end[i].pose().orientation().components(), truth.q[i]))
            .fold(0.0_f64, f64::max);
        let spin = (0..2)
            .map(|i| {
                end[i]
                    .angular_momentum_body()
                    .sub(truth.h[i])
                    .norm_squared()
                    .sqrt()
            })
            .fold(0.0_f64, f64::max);
        errors.push((angle, spin, (work - truth.work).abs()));
    }
    eprintln!("internal actuator N=8/16/32 orientation, spin, work errors: {errors:?}");
    for pair in errors.windows(2) {
        assert!(pair[0].0 / pair[1].0 > 3.7 && pair[0].0 / pair[1].0 < 4.3);
        assert!(pair[0].1 / pair[1].1 > 3.7 && pair[0].1 / pair[1].1 < 4.3);
        assert!(pair[0].2 / pair[1].2 > 3.7 && pair[0].2 / pair[1].2 < 4.3);
    }
    assert!(errors[2].0 < 1e-4);
    assert!(errors[2].1 < 1e-4);
    assert!(errors[2].2 < 1e-4);
}

fn multiply(a: UnitQuaternion, b: UnitQuaternion) -> UnitQuaternion {
    let [w, x, y, z] = qmul(a.components(), b.components());
    UnitQuaternion::new(w, x, y, z).unwrap()
}
#[test]
fn g3_motion_is_frame_boost_covariant_and_reversible_for_the_held_command() {
    let before = initial();
    let props = properties();
    let load = drive();
    let h = 0.12;
    let original =
        step_internal_wrench(before, props, load, h, controls(Vec3::ZERO), || false).unwrap();
    let q = UnitQuaternion::from_axis_angle(Vec3::new(2.0, -1.0, 3.0), 0.7).unwrap();
    let shift = Vec3::new(32.0, -16.0, 8.0);
    let boost = Vec3::new(4.0, -2.0, 1.0);
    let moved = core::array::from_fn(|i| {
        state(
            q.rotate_body_to_world(before[i].pose().position_world())
                .add(shift),
            multiply(q, before[i].pose().orientation()),
            q.rotate_body_to_world(before[i].linear_momentum_world())
                .add(boost.scale(props[i].mass())),
            before[i].angular_momentum_body(),
        )
    });
    let transformed = step_internal_wrench(
        moved,
        props,
        PrescribedInternalWrench {
            force_on_second_world_n: q.rotate_body_to_world(load.force_on_second_world_n),
            couple_on_second_world_nm: q.rotate_body_to_world(load.couple_on_second_world_nm),
            ..load
        },
        h,
        controls(shift),
        || false,
    )
    .unwrap();
    for i in 0..2 {
        let original_after = original.states_after[i];
        let actual = transformed.states_after[i];
        vector(
            actual.pose().position_world(),
            q.rotate_body_to_world(original_after.pose().position_world())
                .add(shift)
                .add(boost.scale(h)),
            2e-11,
        );
        vector(
            actual.linear_momentum_world(),
            q.rotate_body_to_world(original_after.linear_momentum_world())
                .add(boost.scale(props[i].mass())),
            2e-11,
        );
        vector(
            actual.angular_momentum_body(),
            original_after.angular_momentum_body(),
            2e-11,
        );
        close(
            qerror(
                actual.pose().orientation().components(),
                multiply(q, original_after.pose().orientation()).components(),
            ),
            0.0,
            2e-11,
        );
    }
    close(transformed.actuator_work_j, original.actuator_work_j, 2e-11);
    let reversed = original.states_after.map(|b| {
        state(
            b.pose().position_world(),
            b.pose().orientation(),
            b.linear_momentum_world().scale(-1.0),
            b.angular_momentum_body().scale(-1.0),
        )
    });
    let back =
        step_internal_wrench(reversed, props, load, h, controls(Vec3::ZERO), || false).unwrap();
    for i in 0..2 {
        vector(
            back.states_after[i].pose().position_world(),
            before[i].pose().position_world(),
            2e-11,
        );
        close(
            qerror(
                back.states_after[i].pose().orientation().components(),
                before[i].pose().orientation().components(),
            ),
            0.0,
            2e-11,
        );
        vector(
            back.states_after[i].linear_momentum_world(),
            before[i].linear_momentum_world().scale(-1.0),
            2e-11,
        );
        vector(
            back.states_after[i].angular_momentum_body(),
            before[i].angular_momentum_body().scale(-1.0),
            2e-11,
        );
    }
    close(back.actuator_work_j, -original.actuator_work_j, 2e-11);
}

#[test]
fn g3_zero_drive_reuses_native_free_drift_and_rechecks_each_body() {
    let before = initial();
    let props = properties();
    let h = 0.08;
    let result =
        step_internal_wrench(before, props, no_drive(), h, controls(Vec3::ZERO), || false).unwrap();
    let native = core::array::from_fn(|i| {
        RigidBodyIntegrator::new(Gravity::ZERO)
            .step(before[i], props[i], Wrench::ZERO, h)
            .unwrap()
            .state_after
    });
    assert_eq!(result.states_after, native);
    assert_eq!(result.actuator_work_j, 0.0);
    close(
        result.kinetic_energy_change_j,
        result.native_drift_energy_defect_j,
        1e-14,
    );
    // At a distant origin the two representability-induced orbital defects
    // would cancel in the pair's total. Each free body must still pass its gate.
    let unit = [MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap(); 2];
    let distant = [1.0, -1.0].map(|y| {
        state(
            Vec3::new(2.0_f64.powi(54), 0.0, 0.0),
            UnitQuaternion::IDENTITY,
            Vec3::new(1.0, y, 0.0),
            Vec3::ZERO,
        )
    });
    assert!(
        matches!(step_internal_wrench(distant,unit,no_drive(),1.0,controls(Vec3::ZERO),||false),
        Err(InternalWrenchStepError::Balance {body:Some(0),angular_residual_nms,..}) if angular_residual_nms>=1.0)
    );
}

#[test]
fn g4_late_cancellation_and_unrepresentable_drive_publish_no_partial_motion() {
    let before = initial();
    let props = properties();
    let load = drive();
    let mut count = 0;
    let accepted = step_internal_wrench(before, props, load, 0.1, controls(Vec3::ZERO), || {
        count += 1;
        false
    })
    .unwrap();
    assert!(count >= 18);
    for cancel_at in 1..=count {
        let mut calls = 0;
        assert_eq!(
            step_internal_wrench(before, props, load, 0.1, controls(Vec3::ZERO), || {
                calls += 1;
                calls == cancel_at
            }),
            Err(InternalWrenchStepError::Cancelled)
        );
        assert_eq!(before, initial());
    }
    for h in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert_eq!(
            step_internal_wrench(before, props, load, h, controls(Vec3::ZERO), || false),
            Err(InternalWrenchStepError::Dynamics(
                fs_mbd::DynamicsError::InvalidStepDuration
            ))
        );
    }
    for (force, h) in [
        (f64::MAX, 4.0),
        (f64::from_bits(1), 1.0),
        (1.0, f64::from_bits(1)),
    ] {
        let invalid = PrescribedInternalWrench {
            force_on_second_world_n: Vec3::new(force, 0.0, 0.0),
            ..load
        };
        assert_eq!(
            step_internal_wrench(before, props, invalid, h, controls(Vec3::ZERO), || false),
            Err(InternalWrenchStepError::UnrepresentableImpulse)
        );
    }
    let invalid = PrescribedInternalWrench {
        arms_body_m: [Vec3::new(f64::NAN, 0.0, 0.0), Vec3::ZERO],
        ..load
    };
    assert!(matches!(
        step_internal_wrench(before, props, invalid, 0.1, controls(Vec3::ZERO), || false),
        Err(InternalWrenchStepError::Impulse(
            InternalWrenchImpulseError::Dynamics(_)
        ))
    ));
    assert_eq!(
        accepted,
        step_internal_wrench(before, props, load, 0.1, controls(Vec3::ZERO), || false).unwrap()
    );
}
