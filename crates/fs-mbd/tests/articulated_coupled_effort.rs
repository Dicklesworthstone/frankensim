//! Full passive constitutive coupling and internal control on a recoiling tree.

use fs_ga::{Mat3, Se3, Twist, Vec3, Wrench};
use fs_mbd::articulated::motion::spatial_loads::bushing::coupled::{
    self as full, BushingCoordinateScales, CoupledSpatialBushing, PassiveBushingMatrix,
};
use fs_mbd::articulated::motion::spatial_loads::bushing::{self, MaterialFrame, SpatialBushing};
use fs_mbd::articulated::motion::spatial_loads::implicit::ImplicitSpatialConfig;
use fs_mbd::articulated::motion::spatial_loads::{SpatialMotionLoads, SpatialMotionStep};
use fs_mbd::articulated::motion::{
    FreeFloatingMotionError, FreeFloatingMotionLoads, FreeFloatingMotionState,
};
use fs_mbd::articulated::{
    ArticulatedModel, FreeFloatingBaseState, JointModel, Link, SpatialInertia,
};

const ZERO: Vec3 = Vec3 {
    x: 0.0,
    y: 0.0,
    z: 0.0,
};
const WRENCHES: [Wrench; 3] = [Wrench::new(ZERO, ZERO); 3];
const HELD: [f64; 2] = [0.3, -0.2];
const K: [[f64; 2]; 2] = [[12.0, 4.0], [4.0, 8.0]];
const C: [[f64; 2]; 2] = [[0.6, 0.1], [0.1, 0.4]];
const MASS: [f64; 2] = [6.0 / 5.0, 10.0 / 7.0];

fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-10, "{a:.16e} != {b:.16e}");
}

fn model() -> ArticulatedModel {
    let inertia = |mass, spin| {
        SpatialInertia::new(
            mass,
            ZERO,
            Mat3 {
                m: [spin, 0.0, 0.0, 0.0, spin, 0.0, 0.0, 0.0, spin],
            },
        )
        .unwrap()
    };
    // Translation and spin share the x axis, so their physical kinematics
    // commute. The carrier recoils in both modes. Independent reduced inertia:
    // mu=3*2/(3+2), J=(3+2)*2/(3+2+2). Only full K/C couple the modes.
    ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(3.0, 3.0),
        ),
        Link::new(
            "slider",
            Some(0),
            Se3::identity(),
            JointModel::prismatic(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia(1.0, 2.0),
        ),
        Link::new(
            "rotor",
            Some(1),
            Se3::identity(),
            JointModel::revolute(Vec3::new(1.0, 0.0, 0.0), None).unwrap(),
            inertia(1.0, 2.0),
        ),
    ])
    .unwrap()
}

fn initial() -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::new(
            Se3::identity(),
            Twist::new(
                Vec3::new(-2.0 * 0.4 / 7.0, 0.0, 0.0),
                Vec3::new(2.0 * 0.1 / 5.0, 0.0, 0.0),
            ),
        ),
        joint_position: vec![0.3, -0.2],
        joint_velocity: vec![-0.1, 0.4],
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
        rest_translation_a: ZERO,
        translation_stiffness: Vec3::new(K[0][0], 0.0, 0.0),
        translation_damping: Vec3::new(C[0][0], 0.0, 0.0),
        rotation_stiffness: Vec3::new(K[1][1], 0.0, 0.0),
        rotation_damping: Vec3::new(C[1][1], 0.0, 0.0),
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
        let mut matrix = [[0.0; 6]; 6];
        for i in 0..2 {
            for j in 0..2 {
                matrix[3 * i][3 * j] = a[i][j];
            }
        }
        PassiveBushingMatrix::try_new(matrix, scales()).unwrap()
    };
    let b = diagonal();
    CoupledSpatialBushing::try_new(b.a, b.b, b.rest_translation_a, embed(K), embed(C)).unwrap()
}

fn loads() -> SpatialMotionLoads<'static> {
    SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &HELD,
            external_body_wrench: &WRENCHES,
            gravity_world: ZERO,
        },
        world_points: &[],
        connections: &[],
    }
}

fn force(q: [f64; 2], v: [f64; 2]) -> [f64; 2] {
    core::array::from_fn(|i| (0..2).map(|j| K[i][j] * q[j] + C[i][j] * v[j]).sum())
}

fn effort(q: [f64; 2], v: [f64; 2]) -> [f64; 2] {
    core::array::from_fn(|i| HELD[i] - 2.0 * q[i] - 0.3 * v[i])
}

fn control(
    _: f64,
    stage: &FreeFloatingMotionState,
    output: &mut [f64],
) -> Result<(), FreeFloatingMotionError> {
    assert_eq!(output, &HELD);
    output.copy_from_slice(&effort(
        [stage.joint_position[0], stage.joint_position[1]],
        [stage.joint_velocity[0], stage.joint_velocity[1]],
    ));
    Ok(())
}

fn check(step: &SpatialMotionStep, qm: [f64; 2], vm: [f64; 2], h: f64) {
    let q0 = [0.3, -0.2];
    let v0 = [-0.1, 0.4];
    let u = effort(qm, vm);
    let f = force(qm, vm);
    let q1: [f64; 2] = core::array::from_fn(|i| q0[i] + h * vm[i]);
    let v1: [f64; 2] = core::array::from_fn(|i| v0[i] + h * (u[i] - f[i]) / MASS[i]);
    for i in 0..2 {
        close(step.motion.state_after.joint_position[i], q1[i]);
        close(step.motion.state_after.joint_velocity[i], v1[i]);
    }
    close(
        step.motion.state_after.base.twist_body.linear.x,
        -2.0 * v1[0] / 5.0,
    );
    close(
        step.motion.state_after.base.twist_body.angular.x,
        -2.0 * v1[1] / 7.0,
    );
    close(step.motion.momentum_after.linear_momentum_world.norm(), 0.0);
    close(
        step.motion.momentum_after.angular_momentum_world.norm(),
        0.0,
    );
    close(step.motion.center_of_mass_world_after.x, 2.0 * q0[0] / 5.0);
    close(step.motion.actuator_work, h * (u[0] * vm[0] + u[1] * vm[1]));
    let quadratic = |a: [[f64; 2]; 2], x: [f64; 2]| -> f64 {
        (0..2)
            .map(|i| (0..2).map(|j| x[i] * a[i][j] * x[j]).sum::<f64>())
            .sum()
    };
    let kinetic = |v: [f64; 2]| 0.5 * (MASS[0] * v[0] * v[0] + MASS[1] * v[1] * v[1]);
    close(step.motion.momentum_after.kinetic_energy, kinetic(v1));
    close(step.stored_energy_before, 0.5 * quadratic(K, q0));
    close(step.stored_energy_after, 0.5 * quadratic(K, q1));
    close(step.dissipated_energy, h * quadratic(C, vm));
    close(
        step.mechanical_work_residual,
        kinetic(v1) - kinetic(v0) + 0.5 * quadratic(K, q1) - 0.5 * quadratic(K, q0)
            + h * quadratic(C, vm)
            - h * (u[0] * vm[0] + u[1] * vm[1]),
    );
}

#[test]
fn explicit_feedback_uses_actual_coupled_stages_and_both_carrier_recoils() {
    let model = model();
    let state = initial();
    let saved = state.clone();
    let q = [0.3, -0.2];
    let v = [-0.1, 0.4];
    let u = effort(q, v);
    let f = force(q, v);
    let h = 0.04;
    let qm = core::array::from_fn(|i| q[i] + h * v[i] / 2.0);
    let vm = core::array::from_fn(|i| v[i] + h * (u[i] - f[i]) / (2.0 * MASS[i]));
    let mut calls = 0;
    let step = full::step_free_floating_bushings_with_effort(
        &model,
        &state,
        loads(),
        &[coupled()],
        h,
        |elapsed, stage, out| {
            let (expected_q, expected_v) = if calls == 0 { (q, v) } else { (qm, vm) };
            close(elapsed, if calls == 0 { 0.0 } else { h / 2.0 });
            for i in 0..2 {
                close(stage.joint_position[i], expected_q[i]);
                close(stage.joint_velocity[i], expected_v[i]);
            }
            close(stage.base.twist_body.linear.x, -2.0 * expected_v[0] / 5.0);
            close(stage.base.twist_body.angular.x, -2.0 * expected_v[1] / 7.0);
            calls += 1;
            control(elapsed, stage, out)
        },
        || false,
    )
    .unwrap();
    assert_eq!(calls, 2);
    check(&step, qm, vm, h);
    assert_eq!(state, saved);
}

#[test]
fn stiff_full_tensor_and_feedback_match_independent_free_inertia_midpoint() {
    let model = model();
    let state = initial();
    let q = [0.3, -0.2];
    let v = [-0.1, 0.4];
    let h = 3.0;
    // Independent exact 2x2 midpoint map for the two reduced inertias.
    let a: [[f64; 2]; 2] = core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            let diagonal = if i == j { 1.0 } else { 0.0 };
            diagonal
                + h * (C[i][j] + 0.3 * diagonal) / (2.0 * MASS[i])
                + h * h * (K[i][j] + 2.0 * diagonal) / (4.0 * MASS[i])
        })
    });
    let rhs: [f64; 2] = core::array::from_fn(|i| {
        let kq: f64 = (0..2).map(|j| K[i][j] * q[j]).sum();
        v[i] + h * (HELD[i] - kq - 2.0 * q[i]) / (2.0 * MASS[i])
    });
    let determinant = a[0][0] * a[1][1] - a[0][1] * a[1][0];
    let vm = [
        (a[1][1] * rhs[0] - a[0][1] * rhs[1]) / determinant,
        (a[0][0] * rhs[1] - a[1][0] * rhs[0]) / determinant,
    ];
    let qm = core::array::from_fn(|i| q[i] + h * vm[i] / 2.0);
    let mut calls = 0;
    let step = full::step_free_floating_bushings_implicit_with_effort(
        &model,
        &state,
        loads(),
        &[coupled()],
        h,
        ImplicitSpatialConfig {
            residual_tolerance: 1e-12,
            ..Default::default()
        },
        |elapsed, stage, out| {
            close(elapsed, h / 2.0);
            calls += 1;
            control(elapsed, stage, out)
        },
        || false,
    )
    .unwrap();
    check(&step.spatial, qm, vm, h);
    for i in 0..2 {
        close(step.midpoint_state.joint_position[i], qm[i]);
        close(step.midpoint_state.joint_velocity[i], vm[i]);
    }
    assert!(calls > 2);
    assert_eq!(calls, step.residual_evaluations);
    assert!(step.scaled_residual_norm <= 1e-12);
    close(step.spatial.mechanical_work_residual, 0.0);
}

#[test]
fn full_diagonal_feedback_retains_original_explicit_and_implicit_results() {
    let model = model();
    let state = initial();
    let b = diagonal();
    let full = CoupledSpatialBushing::try_from_diagonal(b, scales()).unwrap();
    let h = 0.02;
    assert_eq!(
        bushing::step_free_floating_bushings_with_effort(
            &model,
            &state,
            loads(),
            &[b],
            h,
            control,
            || false,
        )
        .unwrap(),
        full::step_free_floating_bushings_with_effort(
            &model,
            &state,
            loads(),
            std::slice::from_ref(&full),
            h,
            control,
            || false,
        )
        .unwrap(),
    );
    let config = ImplicitSpatialConfig::default();
    assert_eq!(
        bushing::step_free_floating_bushings_implicit_with_effort(
            &model,
            &state,
            loads(),
            &[b],
            h,
            config,
            control,
            || false,
        )
        .unwrap(),
        full::step_free_floating_bushings_implicit_with_effort(
            &model,
            &state,
            loads(),
            &[full],
            h,
            config,
            control,
            || false,
        )
        .unwrap(),
    );
}
