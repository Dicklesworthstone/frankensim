//! Full passive coupling uses the original supported dynamics and work accounts.

use fs_ga::{Mat3, Se3, So3, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::FreeFloatingMotionLoads;
use fs_mbd::articulated::motion::prescribed::{PrescribedBaseMotionState, PrescribedBasePath};
use fs_mbd::articulated::motion::spatial_loads::SpatialMotionLoads;
use fs_mbd::articulated::motion::spatial_loads::bushing::coupled::{
    self as full, BushingCoordinateScales, CoupledSpatialBushing, PassiveBushingMatrix,
};
use fs_mbd::articulated::motion::spatial_loads::bushing::{MaterialFrame, SpatialBushing};
use fs_mbd::articulated::motion::spatial_loads::prescribed::{
    PrescribedSpatialMotionStep, step_prescribed_base_bushings,
};
use fs_mbd::articulated::{ArticulatedModel, BaseState, JointModel, Link, SpatialInertia};

const K: [[f64; 2]; 2] = [[12.0, 4.0], [4.0, 8.0]];
const C: [[f64; 2]; 2] = [[0.6, 0.1], [0.1, 0.4]];

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 3e-12, "{a:.16e} != {b:.16e}");
}

fn model() -> ArticulatedModel {
    let inertia = |mass| {
        SpatialInertia::new(
            mass,
            Vec3::default(),
            Mat3 {
                m: [2.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0.0, 0.0, 2.0],
            },
        )
        .unwrap()
    };
    // Both moving centres lie on the slider axis. The independent effective
    // masses are exactly 2 kg and 2 kg m^2; only the constitutive law couples them.
    ArticulatedModel::new(vec![
        Link::new(
            "support",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(3.0),
        ),
        Link::new(
            "slider",
            Some(0),
            Se3::identity(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia(1.0),
        ),
        Link::new(
            "rotor",
            Some(1),
            Se3::identity(),
            JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
            inertia(1.0),
        ),
    ])
    .unwrap()
}

fn base(time: f64, speed: f64) -> BaseState {
    BaseState::prescribed(
        Se3::from_parts(So3::identity(), Vec3::new(time * speed, 0.0, 0.0)).unwrap(),
        Twist::new(Vec3::default(), Vec3::new(speed, 0.0, 0.0)),
        Twist::zero(),
    )
}

fn initial(speed: f64) -> PrescribedBaseMotionState {
    PrescribedBaseMotionState {
        time_s: 2.0,
        base: base(2.0, speed),
        joint_position: vec![0.3, -0.2],
        joint_velocity: vec![-0.1, 0.4],
    }
}

fn path(state: &PrescribedBaseMotionState, h: f64, speed: f64) -> PrescribedBasePath {
    PrescribedBasePath {
        midpoint: base(state.time_s + h / 2.0, speed),
        endpoint: base(state.time_s + h, speed),
    }
}

fn diagonal() -> SpatialBushing {
    SpatialBushing {
        a: MaterialFrame {
            link: 0,
            link_from_attachment: Se3::identity(),
        },
        b: MaterialFrame {
            link: 2,
            link_from_attachment: Se3::identity(),
        },
        rest_translation_a: Vec3::default(),
        translation_stiffness: Vec3::new(K[0][0], 0.0, 0.0),
        translation_damping: Vec3::new(C[0][0], 0.0, 0.0),
        rotation_stiffness: Vec3::new(0.0, 0.0, K[1][1]),
        rotation_damping: Vec3::new(0.0, 0.0, C[1][1]),
    }
}

fn scales() -> BushingCoordinateScales {
    BushingCoordinateScales {
        translation_m: 0.25,
        rotation_rad: 1.0,
    }
}

fn coupled() -> CoupledSpatialBushing {
    let embed = |a: [[f64; 2]; 2]| {
        let mut full = [[0.0; 6]; 6];
        for i in 0..2 {
            for j in 0..2 {
                full[5 * i][5 * j] = a[i][j];
            }
        }
        PassiveBushingMatrix::try_new(full, scales()).unwrap()
    };
    let b = diagonal();
    CoupledSpatialBushing::try_new(b.a, b.b, b.rest_translation_a, embed(K), embed(C)).unwrap()
}

const ZERO_WRENCHES: [Wrench; 3] = [Wrench::new(
    Vec3 {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    },
    Vec3 {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    },
); 3];

fn loads(effort: &[f64]) -> SpatialMotionLoads<'_> {
    SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: effort,
            external_body_wrench: &ZERO_WRENCHES,
            gravity_world: Vec3::default(),
        },
        world_points: &[],
        connections: &[],
    }
}

fn force(q: [f64; 2], v: [f64; 2]) -> [f64; 2] {
    core::array::from_fn(|i| (0..2).map(|j| K[i][j] * q[j] + C[i][j] * v[j]).sum())
}

fn check(
    step: &PrescribedSpatialMotionStep,
    q: [f64; 2],
    v: [f64; 2],
    qh: [f64; 2],
    vh: [f64; 2],
    uh: [f64; 2],
    speed: f64,
    h: f64,
) {
    let f = force(qh, vh);
    for i in 0..2 {
        close(step.motion.state_after.joint_position[i], q[i] + h * vh[i]);
        close(
            step.motion.state_after.joint_velocity[i],
            v[i] + h * (uh[i] - f[i]) / 2.0,
        );
    }
    close(step.motion.support_wrench_body.force.x, uh[0] - f[0]);
    close(step.motion.support_wrench_body.torque.z, uh[1] - f[1]);
    close(step.motion.support_work, h * speed * (uh[0] - f[0]));
    close(
        step.motion.actuator_work,
        h * (uh[0] * vh[0] + uh[1] * vh[1]),
    );
    close(
        step.motion.internal_body_work,
        -h * (f[0] * vh[0] + f[1] * vh[1]),
    );
    close(
        step.motion.linear_momentum_change_world.x,
        h * (uh[0] - f[0]),
    );
    close(
        step.motion.angular_momentum_change_world.z,
        h * (uh[1] - f[1]),
    );
    close(step.motion.linear_momentum_residual_world.x, 0.0);
    close(step.motion.angular_momentum_residual_world.z, 0.0);
    let loss: f64 = (0..2)
        .map(|i| (0..2).map(|j| vh[i] * C[i][j] * vh[j]).sum::<f64>())
        .sum();
    close(step.dissipated_energy, h * loss);
    close(step.motion.external_work, 0.0);
}

#[test]
fn cross_translation_rotation_coupling_drives_joint_and_support_momentum() {
    let model = model();
    for speed in [0.0, 0.7] {
        let state = initial(speed);
        let before = state.clone();
        let q = [0.3, -0.2];
        let v = [-0.1, 0.4];
        let u = [0.3, -0.2];
        let h = 0.04;
        let f = force(q, v);
        let qh = core::array::from_fn(|i| q[i] + h * v[i] / 2.0);
        let vh = core::array::from_fn(|i| v[i] + h * (u[i] - f[i]) / 4.0);
        let step = full::step_prescribed_base_bushings(
            &model,
            &state,
            path(&state, h, speed),
            loads(&u),
            &[coupled()],
            h,
            || false,
        )
        .unwrap();
        check(&step, q, v, qh, vh, u, speed, h);
        assert_eq!(state, before);
        assert_eq!(step.motion.state_after.base, base(state.time_s + h, speed));
    }
}

#[test]
fn feedback_and_full_tensor_use_the_same_actual_supported_stages() {
    let model = model();
    let speed = 0.7;
    let state = initial(speed);
    let h = 0.03;
    let q = [0.3, -0.2];
    let v = [-0.1, 0.4];
    let held = [0.3, -0.2];
    let feedback = |t: f64, q: [f64; 2], v: [f64; 2]| {
        [
            held[0] - 3.0 * q[0] - 0.2 * v[0] + t,
            held[1] - 2.0 * q[1] + 0.5 * t,
        ]
    };
    let u0 = feedback(0.0, q, v);
    let f0 = force(q, v);
    let qh = core::array::from_fn(|i| q[i] + h * v[i] / 2.0);
    let vh = core::array::from_fn(|i| v[i] + h * (u0[i] - f0[i]) / 4.0);
    let uh = feedback(h / 2.0, qh, vh);
    let mut calls = 0;
    let step = full::step_prescribed_base_bushings_with_effort(
        &model,
        &state,
        path(&state, h, speed),
        loads(&held),
        &[coupled()],
        h,
        |elapsed, stage, out| {
            assert_eq!(out, &held);
            let (expected_q, expected_v) = if calls == 0 { (q, v) } else { (qh, vh) };
            close(elapsed, if calls == 0 { 0.0 } else { h / 2.0 });
            close(stage.time_s, state.time_s + elapsed);
            assert_eq!(stage.base, base(state.time_s + elapsed, speed));
            for i in 0..2 {
                close(stage.joint_position[i], expected_q[i]);
                close(stage.joint_velocity[i], expected_v[i]);
            }
            out.copy_from_slice(&feedback(elapsed, expected_q, expected_v));
            calls += 1;
            Ok(())
        },
        || false,
    )
    .unwrap();
    assert_eq!(calls, 2);
    check(&step, q, v, qh, vh, uh, speed, h);
}

#[test]
fn diagonal_full_tensor_retains_the_original_supported_trajectory_exactly() {
    let model = model();
    let state = initial(0.7);
    let h = 0.03;
    let b = diagonal();
    let full = CoupledSpatialBushing::try_from_diagonal(b, scales()).unwrap();
    let u = [0.3, -0.2];
    assert_eq!(
        step_prescribed_base_bushings(
            &model,
            &state,
            path(&state, h, 0.7),
            loads(&u),
            &[b],
            h,
            || false,
        )
        .unwrap(),
        full::step_prescribed_base_bushings(
            &model,
            &state,
            path(&state, h, 0.7),
            loads(&u),
            &[full],
            h,
            || false,
        )
        .unwrap(),
    );
}

#[test]
fn stiff_cross_coupling_and_feedback_share_one_implicit_support_solve() {
    use fs_mbd::articulated::motion::spatial_loads::implicit::ImplicitSpatialConfig;

    let model = model();
    let speed = 0.7;
    let state = initial(speed);
    let before = state.clone();
    let q = [0.3, -0.2];
    let v = [-0.1, 0.4];
    let held = [0.3, -0.2];
    // This step is far beyond the explicit oscillator's stability interval.
    // Solve the independent two-coordinate linear midpoint equations directly.
    let h = 3.0;
    let effective_k = [[14.0, 4.0], [4.0, 10.0]];
    let effective_c = [[0.9, 0.1], [0.1, 0.7]];
    let a: [[f64; 2]; 2] = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            (if i == j { 1.0 } else { 0.0 })
                + h * effective_c[i][j] / 4.0
                + h * h * effective_k[i][j] / 8.0
        })
    });
    let b: [f64; 2] = core::array::from_fn(|i| {
        v[i] + h / 4.0 * (held[i] - (0..2).map(|j| effective_k[i][j] * q[j]).sum::<f64>())
    });
    let determinant = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    let vh = [
        (a[1][1] * b[0] - a[0][1] * b[1]) / determinant,
        (a[0][0] * b[1] - a[1][0] * b[0]) / determinant,
    ];
    let qh = core::array::from_fn(|i| q[i] + h * vh[i] / 2.0);
    let uh = core::array::from_fn(|i| held[i] - 2.0 * qh[i] - 0.3 * vh[i]);
    let mut calls = 0;
    let step = full::step_prescribed_base_bushings_implicit_with_effort(
        &model,
        &state,
        path(&state, h, speed),
        loads(&held),
        &[coupled()],
        h,
        ImplicitSpatialConfig {
            residual_tolerance: 1e-12,
            ..Default::default()
        },
        |elapsed, stage, out| {
            close(elapsed, h / 2.0);
            close(stage.time_s, state.time_s + h / 2.0);
            assert_eq!(stage.base, base(state.time_s + h / 2.0, speed));
            assert_eq!(out, &held);
            for i in 0..2 {
                out[i] -= 2.0 * stage.joint_position[i] + 0.3 * stage.joint_velocity[i];
            }
            calls += 1;
            Ok(())
        },
        || false,
    )
    .unwrap();
    check(&step.spatial, q, v, qh, vh, uh, speed, h);
    for i in 0..2 {
        close(step.midpoint_state.joint_position[i], qh[i]);
        close(step.midpoint_state.joint_velocity[i], vh[i]);
    }
    assert!(step.nonlinear.converged);
    assert!(step.scaled_residual_norm <= 1e-12);
    assert_eq!(calls, step.residual_evaluations);
    assert!(calls > 2, "feedback must be inside the nonlinear solve");
    close(step.spatial.mechanical_work_residual, 0.0);
    assert!(step.spatial.dissipated_energy > 0.0);
    assert_eq!(state, before);
}

#[test]
fn implicit_diagonal_full_tensor_retains_native_supported_parity() {
    use fs_mbd::articulated::motion::spatial_loads::prescribed::implicit::{
        ImplicitSpatialConfig, step_prescribed_base_bushings_implicit,
    };
    let model = model();
    let state = initial(0.7);
    let h = 3.0;
    let b = diagonal();
    let full = CoupledSpatialBushing::try_from_diagonal(b, scales()).unwrap();
    let u = [0.3, -0.2];
    let config = ImplicitSpatialConfig::default();
    assert_eq!(
        step_prescribed_base_bushings_implicit(
            &model,
            &state,
            path(&state, h, 0.7),
            loads(&u),
            &[b],
            h,
            config,
            || false,
        )
        .unwrap(),
        full::step_prescribed_base_bushings_implicit(
            &model,
            &state,
            path(&state, h, 0.7),
            loads(&u),
            &[full],
            h,
            config,
            || false,
        )
        .unwrap(),
    );
}
