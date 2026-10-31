//! Integrated dry/fluid inertia, recoil and midpoint physical oracles.
use fs_flyer::Refusal;
use fs_flyer::addedmass::articulated::{LinkedStrip, assemble_articulated_strips};
use fs_flyer::addedmass::{Strip, solve_articulated_joint};
use fs_flyer::partitioned::MemoryState;
use fs_flyer::partitioned::articulated::{
    ArticulatedLoads, ArticulatedState, ArticulatedStepConfig, JointImpedance, step,
};
use fs_ga::{Mat3, Se3, So3, So3Tangent, Twist, Vec3, Wrench};
use fs_mbd::articulated::generalized::{GeneralizedKinematics, free_floating_terms};
use fs_mbd::articulated::motion::FreeFloatingMotionState;
use fs_mbd::articulated::{
    ArticulatedModel, BaseState, FreeFloatingBaseState, JointLimits, JointModel, Link,
    SpatialInertia, forward_kinematics, free_floating_forward_dynamics, free_floating_momentum,
};
use std::cell::Cell;

const ZERO_IMPEDANCE: JointImpedance = JointImpedance {
    damping: 0.0,
    stiffness: 0.0,
    rest_position: 0.0,
};

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e} (tol {tolerance:e})"
    );
}
fn vector(a: [f64; 3]) -> Vec3 {
    Vec3::new(a[0], a[1], a[2])
}
fn pose(rotation: [f64; 3], position: [f64; 3]) -> Se3 {
    Se3::from_parts(
        So3::exp(So3Tangent::new(vector(rotation))).unwrap(),
        vector(position),
    )
    .unwrap()
}
fn inertia(mass: f64, center: [f64; 3], diagonal: [f64; 3]) -> SpatialInertia {
    let [x, y, z] = diagonal;
    SpatialInertia::new(
        mass,
        vector(center),
        Mat3 {
            m: [x, 0.0, 0.0, 0.0, y, 0.0, 0.0, 0.0, z],
        },
    )
    .unwrap()
}
fn config() -> ArticulatedStepConfig<'static> {
    ArticulatedStepConfig {
        max_iterations: 64,
        linear_tolerance_mps: 1e-12,
        angular_tolerance_rad_s: 1e-12,
        joint_velocity_tolerances: &[1e-12; 2],
    }
}
fn plate(chord: f64, span: f64, position: [f64; 3], normal: [f64; 3], axis: [f64; 3]) -> Strip {
    Strip {
        name: "attached plate",
        chord_m: chord,
        span_m: span,
        position_m: position,
        normal,
        span_axis: axis,
        control_coord: None,
        control_gain: 0.0,
        control_angular_gain: [0.0; 3],
    }
}
fn loads(
    map: &GeneralizedKinematics,
    strips: &[LinkedStrip],
    rho: f64,
) -> Result<ArticulatedLoads, Refusal> {
    Ok(ArticulatedLoads {
        fluid: assemble_articulated_strips(rho, strips, map)?,
        world_wrenches: vec![Wrench::default(); map.link_count()],
        gravity_world: Vec3::default(),
        joint_effort: vec![0.0; map.dof_count()],
    })
}

fn oblique_tree() -> ArticulatedModel {
    let first = SpatialInertia::new(
        1.8,
        vector([0.25, 0.05, -0.1]),
        Mat3 {
            m: [0.15, 0.012, -0.008, 0.012, 0.2, 0.01, -0.008, 0.01, 0.25],
        },
    )
    .unwrap();
    ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(4.0, [0.1, -0.05, 0.02], [2.0, 2.5, 3.0]),
        ),
        Link::new(
            "oblique hinge",
            Some(0),
            pose([0.1, -0.2, 0.15], [0.5, 0.2, -0.1]),
            JointModel::revolute(vector([0.3, 0.8, 0.4]), None).unwrap(),
            first,
        ),
        Link::new(
            "moving screw",
            Some(1),
            pose([-0.15, 0.1, 0.2], [0.4, -0.2, 0.3]),
            JointModel::helical(vector([-0.4, 0.3, 0.8]), 0.12, None).unwrap(),
            inertia(1.1, [0.2, -0.1, 0.15], [0.1, 0.12, 0.16]),
        ),
    ])
    .unwrap()
}
fn oblique_state() -> ArticulatedState {
    ArticulatedState {
        motion: FreeFloatingMotionState {
            base: FreeFloatingBaseState::new(
                pose([-0.2, 0.3, 0.1], [0.4, -0.2, 0.8]),
                Twist::new(vector([0.3, -0.25, 0.2]), vector([0.8, -0.3, 0.2])),
            ),
            joint_position: vec![0.35, -0.25],
            joint_velocity: vec![0.7, -0.45],
        },
        memory: vec![],
    }
}

#[test]
fn full_mass_adapter_matches_native_aba_with_bias_gravity_and_link_loads() {
    let model = oblique_tree();
    let state = oblique_state().motion;
    let dry = free_floating_terms(
        &model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
        || false,
    )
    .unwrap();
    let external = [
        Wrench::new(vector([0.2, -0.1, 0.3]), vector([0.4, 0.1, -0.2])),
        Wrench::new(vector([-0.1, 0.3, 0.2]), vector([0.2, -0.5, 0.4])),
        Wrench::new(vector([0.15, -0.2, 0.1]), vector([-0.3, 0.4, 0.2])),
    ];
    let gravity = vector([0.3, -0.2, -9.81]);
    let effort = [0.7, -0.4];
    let native = free_floating_forward_dynamics(
        &model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
        &effort,
        gravity,
        &external,
    )
    .unwrap();
    let mut applied = dry.kinematics.project_body_wrenches(&external).unwrap();
    for (force, weight) in applied
        .iter_mut()
        .zip(dry.kinematics.project_uniform_gravity(gravity).unwrap())
    {
        *force += weight;
    }
    let mut fluid = assemble_articulated_strips(1.2, &[], &dry.kinematics).unwrap();
    fluid.q_rigid_nonaccel[..3].copy_from_slice(&applied[3..6]);
    fluid.q_rigid_nonaccel[3..].copy_from_slice(&applied[..3]);
    for k in 0..2 {
        fluid.q_control_nonaccel[k] = applied[6 + k] + effort[k];
    }
    let (body, joint) = solve_articulated_joint(&dry, &fluid).unwrap();
    let expected = native.base_spatial_acceleration_body.to_array();
    for i in 0..3 {
        close(body[i], expected[i + 3], 2e-11);
        close(body[i + 3], expected[i], 2e-11);
    }
    for (got, expected) in joint.into_iter().zip(native.generalized_acceleration) {
        close(got, expected, 2e-11);
    }
}

fn rotor_tree(limits: Option<JointLimits>) -> ArticulatedModel {
    let joint = JointModel::revolute(vector([0.0, 0.0, 1.0]), limits).unwrap();
    ArticulatedModel::new(vec![
        Link::new(
            "carrier",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(4.0, [0.0; 3], [3.0, 3.0, 4.0]),
        ),
        Link::new(
            "first rotor",
            Some(0),
            Se3::identity(),
            joint,
            inertia(2.0, [0.0; 3], [0.8, 0.8, 1.0]),
        ),
        Link::new(
            "second rotor",
            Some(0),
            Se3::identity(),
            joint,
            inertia(1.0, [0.0; 3], [0.4, 0.4, 0.5]),
        ),
    ])
    .unwrap()
}
fn rotor_state() -> ArticulatedState {
    ArticulatedState {
        motion: FreeFloatingMotionState {
            base: FreeFloatingBaseState::stationary(Se3::identity()),
            joint_position: vec![0.0; 2],
            joint_velocity: vec![0.0; 2],
        },
        memory: vec![],
    }
}

#[test]
fn two_stiff_rotors_match_independent_schur_cayley_recoil_and_damping_work() {
    let model = rotor_tree(None);
    let rho = 1.2;
    let dimensions = [(1.3_f64, 1.8_f64), (0.9_f64, 1.2_f64)];
    let strips: Vec<_> = dimensions
        .iter()
        .enumerate()
        .map(|(i, (chord, span))| {
            LinkedStrip::new(
                i + 1,
                plate(*chord, *span, [0.0; 3], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]),
            )
            .unwrap()
        })
        .collect();
    let j: [f64; 2] = core::array::from_fn(|i| {
        [1.0, 0.5][i]
            + rho * std::f64::consts::PI * dimensions[i].0.powi(4) * dimensions[i].1 / 128.0
    });
    // Independent angular mass for [carrier yaw, relative spin1, spin2].
    // Eliminate the unforced carrier analytically; no production mass query.
    let a = 4.0 + j[0] + j[1];
    let reduced = [
        [j[0] - j[0] * j[0] / a, -j[0] * j[1] / a],
        [-j[0] * j[1] / a, j[1] - j[1] * j[1] / a],
    ];
    let impedance = [
        JointImpedance {
            damping: 4.0,
            stiffness: 200_000.0,
            rest_position: 0.0,
        },
        JointImpedance {
            damping: 3.0,
            stiffness: 130_000.0,
            rest_position: 0.0,
        },
    ];
    let mut state = rotor_state();
    state.motion.joint_position = vec![0.4, -0.25];
    state.motion.joint_velocity = vec![0.3, -0.2];
    state.motion.base.twist_body.angular.z = -(j[0] * 0.3 - j[1] * 0.2) / a;
    let energy = |q: &[f64], v: &[f64]| {
        0.5 * (reduced[0][0] * v[0] * v[0]
            + 2.0 * reduced[0][1] * v[0] * v[1]
            + reduced[1][1] * v[1] * v[1]
            + impedance[0].stiffness * q[0] * q[0]
            + impedance[1].stiffness * q[1] * q[1])
    };
    let dt = 0.03;
    let h = dt / 2.0;
    for tick in 0..24 {
        let q = &state.motion.joint_position;
        let v = &state.motion.joint_velocity;
        let aa = reduced[0][0] + h * impedance[0].damping + h * h * impedance[0].stiffness;
        let bb = reduced[0][1];
        let dd = reduced[1][1] + h * impedance[1].damping + h * h * impedance[1].stiffness;
        let rhs: [f64; 2] = core::array::from_fn(|i| {
            -dt * (impedance[i].damping * v[i] + impedance[i].stiffness * (q[i] + h * v[i]))
        });
        let det = aa * dd - bb * bb;
        let dv = [
            (dd * rhs[0] - bb * rhs[1]) / det,
            (aa * rhs[1] - bb * rhs[0]) / det,
        ];
        let want_v = [v[0] + dv[0], v[1] + dv[1]];
        let want_q = [q[0] + h * (v[0] + want_v[0]), q[1] + h * (v[1] + want_v[1])];
        let out = step(
            &model,
            &impedance,
            &state,
            &[],
            tick as f64 * dt,
            dt,
            config(),
            |_, _, map| loads(map, &strips, rho),
            || false,
        )
        .unwrap();
        assert!(out.iterations <= 4, "{} iterations", out.iterations);
        for i in 0..2 {
            close(out.state.motion.joint_position[i], want_q[i], 4e-12);
            close(out.state.motion.joint_velocity[i], want_v[i], 4e-10);
        }
        close(
            out.state.motion.base.twist_body.angular.z,
            -(j[0] * want_v[0] + j[1] * want_v[1]) / a,
            4e-10,
        );
        let loss: f64 = (0..2)
            .map(|i| {
                dt * impedance[i].damping
                    * (0.5 * (v[i] + out.state.motion.joint_velocity[i])).powi(2)
            })
            .sum();
        close(
            energy(
                &out.state.motion.joint_position,
                &out.state.motion.joint_velocity,
            ) + loss,
            energy(q, v),
            2e-7,
        );
        state = out.state;
    }
}

fn complete_invariants(
    model: &ArticulatedModel,
    state: &ArticulatedState,
    strips: &[(usize, Strip)],
    rho: f64,
) -> (f64, Vec3, Vec3) {
    let state = &state.motion;
    let dry = free_floating_momentum(
        model,
        state.base,
        &state.joint_position,
        &state.joint_velocity,
    )
    .unwrap();
    let kin = forward_kinematics(
        model,
        BaseState::prescribed(
            state.base.world_from_base,
            state.base.twist_body,
            Twist::zero(),
        ),
        &state.joint_position,
        &state.joint_velocity,
    )
    .unwrap();
    let (mut energy, mut linear, mut angular) = (
        dry.kinetic_energy,
        dry.linear_momentum_world,
        dry.angular_momentum_world,
    );
    // Independent local plate law; no assembled generalized matrix/bias or
    // coupled-solver receipt is used to measure the physical balances.
    for (link, strip) in strips {
        let v = kin.body_twist[*link];
        let r = vector(strip.position_m);
        let normal = vector(strip.normal);
        let span = vector(strip.span_axis);
        let chord = span.cross(normal);
        let m = rho * std::f64::consts::PI * strip.chord_m.powi(2) * strip.span_m / 4.0;
        let pitch = m * strip.chord_m.powi(2) / 32.0;
        let spanwise = m * strip.span_m.powi(2) / 12.0;
        let speed = normal.dot(v.linear + v.angular.cross(r));
        let omega_pitch = span.dot(v.angular);
        let omega_span = chord.dot(v.angular);
        let force = normal.scale(m * speed);
        let torque =
            r.cross(force) + span.scale(pitch * omega_pitch) + chord.scale(spanwise * omega_span);
        let world = Wrench::new(torque, force).transform_by(&kin.world_from_link[*link]);
        energy += 0.5
            * (m * speed * speed
                + pitch * omega_pitch * omega_pitch
                + spanwise * omega_span * omega_span);
        linear = linear + world.force;
        angular = angular + world.torque;
    }
    (energy, linear, angular)
}

#[test]
fn serial_moving_com_and_fluid_have_second_order_energy_and_world_momentum_balance() {
    let model = oblique_tree();
    let initial = oblique_state();
    let rho = 0.7;
    let geometry = [
        (
            1,
            plate(
                1.2,
                1.4,
                [0.15, 0.05, -0.1],
                [0.0, 0.0, 1.0],
                [0.0, 1.0, 0.0],
            ),
        ),
        (
            2,
            plate(
                0.8,
                1.2,
                [0.05, -0.15, 0.1],
                [1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
            ),
        ),
    ];
    let strips: Vec<_> = geometry
        .iter()
        .map(|(link, strip)| LinkedStrip::new(*link, strip.clone()).unwrap())
        .collect();
    let before = complete_invariants(&model, &initial, &geometry, rho);
    let errors = [0.02_f64, 0.01, 0.005].map(|dt| {
        let mut state = initial.clone();
        for tick in 0..(0.4 / dt).round() as usize {
            state = step(
                &model,
                &[ZERO_IMPEDANCE; 2],
                &state,
                &[],
                tick as f64 * dt,
                dt,
                config(),
                |_, _, map| loads(map, &strips, rho),
                || false,
            )
            .unwrap()
            .state;
        }
        let after = complete_invariants(&model, &state, &geometry, rho);
        [
            (after.0 - before.0).abs(),
            (after.1 - before.1).norm(),
            (after.2 - before.2).norm(),
        ]
    });
    println!("articulated energy/linear/angular errors at .02/.01/.005: {errors:?}");
    for channel in 0..3 {
        for pair in errors.windows(2) {
            assert!(
                (3.4..4.6).contains(&(pair[0][channel] / pair[1][channel])),
                "{errors:?}, channel {channel}"
            );
        }
    }
    assert!(errors[2][0] < 2e-5 * (1.0 + before.0.abs()));
    assert!(errors[2][1] < 2e-5 * (1.0 + before.1.norm()));
    assert!(errors[2][2] < 2e-5 * (1.0 + before.2.norm()));
}

#[test]
fn memory_sampling_and_atomic_refusals_preserve_native_effort_ownership() {
    let model = rotor_tree(Some(JointLimits::new(-2.0, 2.0, 100.0, 0.2).unwrap()));
    let mut state = rotor_state();
    state.memory = vec![MemoryState {
        tau_s: 0.02,
        value: 0.4,
    }];
    let saved = state.clone();
    let dt = 0.01;
    let sampled = Cell::new(0);
    let out = step(
        &model,
        &[ZERO_IMPEDANCE; 2],
        &state,
        &[1.0],
        2.0,
        dt,
        config(),
        |time, mid, map| {
            close(time, 2.0 + 0.5 * dt, 0.0);
            close(
                mid.memory[0].value,
                1.0 - 0.6 * (-0.5 * dt / 0.02).exp(),
                2e-14,
            );
            sampled.set(sampled.get() + 1);
            let mut result = loads(map, &[], 1.2)?;
            result.joint_effort[0] = 0.2;
            result.world_wrenches[0].force = vector([100.0, 0.0, 0.0]);
            Ok(result)
        },
        || false,
    )
    .unwrap();
    assert!(sampled.get() >= 2);
    close(
        out.state.memory[0].value,
        1.0 - 0.6 * (-dt / 0.02).exp(),
        2e-14,
    );
    assert_eq!(state, saved);

    let cancel = Cell::new(false);
    let error = step(
        &model,
        &[ZERO_IMPEDANCE; 2],
        &state,
        &[1.0],
        0.0,
        dt,
        config(),
        |_, _, map| {
            cancel.set(true);
            loads(map, &[], 1.2)
        },
        || cancel.get(),
    )
    .unwrap_err();
    assert_eq!(error.code, "articulated-step-cancelled");
    assert_eq!(state, saved);
    for (duplicate, excessive) in [(true, false), (false, true)] {
        let error = step(
            &model,
            &[ZERO_IMPEDANCE; 2],
            &state,
            &[1.0],
            0.0,
            dt,
            config(),
            |_, _, map| {
                let mut result = loads(map, &[], 1.2)?;
                if duplicate {
                    result.fluid.q_control_nonaccel[0] = 1.0;
                }
                if excessive {
                    result.joint_effort[0] = 0.2_f64.next_up();
                }
                Ok(result)
            },
            || false,
        )
        .unwrap_err();
        assert_eq!(
            error.code,
            if duplicate {
                "coupled-load-ownership"
            } else {
                "articulated-step-invalid"
            }
        );
        assert_eq!(state, saved);
    }
    let mut spring_state = state.clone();
    spring_state.motion.joint_position[0] = 0.1;
    let impedance = [
        JointImpedance {
            stiffness: 3.0,
            ..ZERO_IMPEDANCE
        },
        ZERO_IMPEDANCE,
    ];
    assert!(
        step(
            &model,
            &impedance,
            &spring_state,
            &[1.0],
            0.0,
            dt,
            config(),
            |_, _, map| loads(map, &[], 1.2),
            || false
        )
        .is_err()
    );
    let mut one_iteration = config();
    one_iteration.max_iterations = 1;
    let error = step(
        &model,
        &[ZERO_IMPEDANCE; 2],
        &state,
        &[1.0],
        0.0,
        dt,
        one_iteration,
        |_, _, map| {
            let mut result = loads(map, &[], 1.2)?;
            result.joint_effort[0] = 0.1;
            Ok(result)
        },
        || false,
    )
    .unwrap_err();
    assert_eq!(error.code, "articulated-step-not-converged");
    assert_eq!(state, saved);
}
