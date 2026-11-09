//! G1/G3/G4: reciprocal eccentric bearing forces on the original full body.
use fs_ga::{Mat3, Se3, So3, So3Tangent, Vec3 as GaVec3};
use fs_mbd::articulated::SpatialInertia;
use fs_mbd::rotordynamics::rigid::{
    ElasticEccentricRotor, ElasticEccentricRotorError, RadialBearing,
};
use fs_mbd::{Gravity, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3};
use std::cell::Cell;

const MASS: f64 = 3.0;
const ECCENTRICITY: f64 = 0.25;
const POLAR: f64 = 0.8;

fn source(offset: GaVec3, tensor: [f64; 9]) -> SpatialInertia {
    SpatialInertia::new(MASS, offset, Mat3 { m: tensor }).unwrap()
}
fn bearing() -> RadialBearing {
    RadialBearing {
        station_m: 0.0,
        stiffness_n_per_m: [[30.0, 5.0], [5.0, 18.0]],
        damping_ns_per_m: [[0.0; 2]; 2],
    }
}
fn rotor(bearings: &[RadialBearing]) -> ElasticEccentricRotor {
    ElasticEccentricRotor::new(
        source(
            GaVec3::new(ECCENTRICITY, 0.0, 0.0),
            [1.2, 0.0, 0.0, 0.0, 1.2, 0.0, 0.0, 0.0, POLAR],
        ),
        Se3::identity(),
        bearings,
        || false,
    )
    .unwrap()
}
fn integrator() -> RigidBodyIntegrator {
    RigidBodyIntegrator::new(Gravity::ZERO)
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + b.abs()),
        "{a:e} != {b:e}"
    );
}
fn components(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn state(position: Vec3, phase: f64, velocity: Vec3, speed: f64) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(
            position,
            UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), phase).unwrap(),
        )
        .unwrap(),
        velocity.scale(MASS),
        Vec3::new(0.0, 0.0, POLAR * speed),
    )
    .unwrap()
}

#[test]
fn g1_shaft_acceleration_solves_the_reciprocal_planar_mass_equations() {
    let rotor = rotor(&[bearing()]);
    let (phase, speed, torque) = (0.7_f64, 2.4, -0.35);
    let radial = Vec3::new(phase.cos(), phase.sin(), 0.0);
    let tangent = Vec3::new(-phase.sin(), phase.cos(), 0.0);
    let shaft = Vec3::new(0.2, -0.1, 0.0);
    let initial = state(
        shaft.add(radial.scale(ECCENTRICITY)),
        phase,
        Vec3::ZERO,
        speed,
    );
    let motion = rotor
        .prepare_motion(integrator(), initial, torque, || false)
        .unwrap()
        .point_motion(rotor.shaft_arm_body(0.0).unwrap())
        .unwrap();
    let force = Vec3::new(-5.5, 0.8, 0.0); // Full original K, including off-diagonal terms.
    let u = MASS * ECCENTRICITY;
    let mut equations = [
        [
            MASS,
            0.0,
            u * tangent.x,
            force.x + u * speed * speed * radial.x,
        ],
        [
            0.0,
            MASS,
            u * tangent.y,
            force.y + u * speed * speed * radial.y,
        ],
        [
            u * tangent.x,
            u * tangent.y,
            POLAR + MASS * ECCENTRICITY.powi(2),
            torque,
        ],
    ];
    // Independent dense solve in shaft-centre coordinates, including the
    // parallel-axis term and both off-diagonal inertial couplings.
    for pivot in 0..3 {
        let d = equations[pivot][pivot];
        for j in pivot..4 {
            equations[pivot][j] /= d;
        }
        let pivot_row = equations[pivot];
        for row in 0..3 {
            if row != pivot {
                let factor = equations[row][pivot];
                for j in pivot..4 {
                    equations[row][j] -= factor * pivot_row[j];
                }
            }
        }
    }
    near(motion.point_acceleration_world.x, equations[0][3], 2e-14);
    near(motion.point_acceleration_world.y, equations[1][3], 2e-14);
    near(motion.angular_acceleration_body.z, equations[2][3], 2e-14);
    near(
        motion.center_of_mass_acceleration_world.x,
        force.x / MASS,
        2e-14,
    );
    near(
        motion.center_of_mass_acceleration_world.y,
        force.y / MASS,
        2e-14,
    );
    assert!((motion.angular_acceleration_body.z - torque / POLAR).abs() > 1.0);
    near(
        rotor.mass_properties().inertia_tensor_body().m[8],
        POLAR,
        2e-15,
    );
    assert_eq!(rotor.source_inertia().center_of_mass().x, ECCENTRICITY);
}

#[test]
fn g1_finite_torque_reverses_free_spin_and_retains_tangential_shaft_acceleration() {
    let rotor = rotor(&[]);
    let initial = state(Vec3::new(ECCENTRICITY, 0.0, 0.0), 0.0, Vec3::ZERO, -1.0);
    let receipt = rotor
        .step(integrator(), initial, 2.0, 1.0, || false)
        .unwrap();
    let after = receipt.native.state_after;
    let motion = rotor
        .prepare_motion(integrator(), after, 2.0, || false)
        .unwrap()
        .point_motion(rotor.shaft_arm_body(0.0).unwrap())
        .unwrap();
    near(motion.kinematics.angular_velocity_body.z, 1.5, 3e-15);
    near(motion.angular_acceleration_body.z, 2.5, 3e-15);
    assert_eq!(
        after.pose().position_world(),
        initial.pose().position_world()
    );
    let (s, c) = 0.25_f64.sin_cos();
    near(
        motion.kinematics.point_world.x,
        ECCENTRICITY * (1.0 - c),
        3e-15,
    );
    near(motion.kinematics.point_world.y, -ECCENTRICITY * s, 3e-15);
    near(
        motion.point_acceleration_world.x,
        ECCENTRICITY * (2.5 * s + 2.25 * c),
        3e-15,
    );
    near(
        motion.point_acceleration_world.y,
        ECCENTRICITY * (-2.5 * c + 2.25 * s),
        3e-15,
    );
    near(receipt.shaft_work_j, 0.5, 3e-15);
    near(receipt.balance_residual_j, 0.0, 3e-15);
    assert_eq!(receipt.bearing_potential_change_j, 0.0);
}

// Independent COM-coordinate ODE for the planar eccentric body. Only this
// oracle uses RK4; production uses the existing native rigid-body Verlet map.
fn derivative(y: [f64; 6], torque: f64) -> [f64; 6] {
    let (s, c) = y[2].sin_cos();
    let q = [y[0] - ECCENTRICITY * c, y[1] - ECCENTRICITY * s];
    let force = [-30.0 * q[0] - 5.0 * q[1], -5.0 * q[0] - 18.0 * q[1]];
    [
        y[3],
        y[4],
        y[5],
        force[0] / MASS,
        force[1] / MASS,
        (torque - ECCENTRICITY * (-s * force[0] + c * force[1])) / POLAR,
    ]
}
fn reference(torque: f64) -> [f64; 6] {
    let mut y = [ECCENTRICITY, 0.0, 0.0, 0.0, 0.0, 3.0];
    let h = 2.0 / 8000.0;
    for _ in 0..8000 {
        let k1 = derivative(y, torque);
        let k2 = derivative(core::array::from_fn(|i| y[i] + h * k1[i] / 2.0), torque);
        let k3 = derivative(core::array::from_fn(|i| y[i] + h * k2[i] / 2.0), torque);
        let k4 = derivative(core::array::from_fn(|i| y[i] + h * k3[i]), torque);
        y = core::array::from_fn(|i| y[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0);
    }
    y
}
fn trajectory(steps: usize, torque: f64, reference: [f64; 6]) -> (f64, f64, f64) {
    let rotor = rotor(&[bearing()]);
    let mut current = state(Vec3::new(ECCENTRICITY, 0.0, 0.0), 0.0, Vec3::ZERO, 3.0);
    let mut defect = 0.0;
    for _ in 0..steps {
        let step = rotor
            .step(integrator(), current, torque, 2.0 / steps as f64, || false)
            .unwrap();
        defect += step.balance_residual_j;
        current = step.native.state_after;
    }
    let p = current.pose().position_world();
    let v = current
        .center_of_mass_velocity_world(rotor.mass_properties())
        .unwrap();
    let w = rotor
        .mass_properties()
        .angular_velocity_body_checked(current.angular_momentum_body())
        .unwrap()
        .z;
    let axis = current
        .pose()
        .orientation()
        .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
    let error = (p.x - reference[0]).abs()
        + (p.y - reference[1]).abs()
        + (v.x - reference[3]).abs()
        + (v.y - reference[4]).abs()
        + (w - reference[5]).abs()
        + (axis.x - reference[2].cos()).abs()
        + (axis.y - reference[2].sin()).abs();
    (error, defect.abs(), w)
}

#[test]
fn g3_bearing_orbit_changes_spin_and_the_native_trajectory_converges_at_second_order() {
    for torque in [0.0, 0.12] {
        let oracle = reference(torque);
        let coarse = trajectory(50, torque, oracle);
        let medium = trajectory(100, torque, oracle);
        let fine = trajectory(200, torque, oracle);
        assert!((fine.2 - (3.0 + torque * 2.0 / POLAR)).abs() > 0.01);
        assert!(
            coarse.0 / medium.0 > 3.8 && medium.0 / fine.0 > 3.8,
            "{coarse:?} {medium:?} {fine:?}"
        );
        assert!(
            coarse.1 / medium.1 > 3.7 && medium.1 / fine.1 > 3.7,
            "{coarse:?} {medium:?} {fine:?}"
        );
        assert!(fine.0 < 1.5e-4 && fine.1 < 3e-4, "{fine:?}");
    }
}

#[test]
fn g3_full_inertia_and_finite_rotation_bearing_jacobian_remain_original() {
    let tensor = [1.3, 0.1, 0.05, 0.1, 1.1, 0.08, 0.05, 0.08, 0.9];
    let rotation = So3::exp(So3Tangent::new(GaVec3::new(0.2, -0.4, 0.6))).unwrap();
    let frame = Se3::from_parts(rotation, GaVec3::new(1.0, -2.0, 0.5)).unwrap();
    let mut second = bearing();
    second.station_m = 0.7;
    let rotor = ElasticEccentricRotor::new(
        source(GaVec3::new(0.2, -0.1, 0.15), tensor),
        frame,
        &[bearing(), second],
        || false,
    )
    .unwrap();
    let props = rotor.mass_properties();
    let reference = RigidBodyState::new(rotor.reference_pose(), Vec3::ZERO, Vec3::ZERO).unwrap();
    near(
        rotor
            .loads(reference, 0.0, || false)
            .unwrap()
            .bearing_potential_j,
        0.0,
        1e-27,
    );
    let pose = Pose::new(
        rotor
            .reference_pose()
            .position_world()
            .add(Vec3::new(0.1, 0.03, -0.02)),
        rotor
            .reference_pose()
            .orientation()
            .right_exp(Vec3::new(0.3, -0.2, 0.4))
            .unwrap(),
    )
    .unwrap();
    let omega = Vec3::new(0.4, -0.7, 2.0);
    let momentum = props.angular_momentum_body_checked(omega).unwrap();
    let current = RigidBodyState::new(pose, Vec3::ZERO, momentum).unwrap();
    let load = rotor.loads(current, 0.0, || false).unwrap();
    let motion = rotor
        .prepare_motion(integrator(), current, 0.0, || false)
        .unwrap()
        .point_motion(Vec3::ZERO)
        .unwrap();
    let alpha = components(motion.angular_acceleration_body);
    let gyroscopic = components(omega.cross(momentum));
    let moment = components(load.wrench.torque_body);
    for i in 0..3 {
        near(
            (0..3).map(|j| tensor[3 * i + j] * alpha[j]).sum::<f64>() + gyroscopic[i],
            moment[i],
            3e-14,
        );
    }
    let axes = [
        Vec3::new(1.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        Vec3::new(0.0, 0.0, 1.0),
    ];
    for (i, axis) in axes.into_iter().enumerate() {
        for angular in [false, true] {
            let potential = |epsilon: f64| {
                let p = if angular {
                    Pose::new(
                        pose.position_world(),
                        pose.orientation().right_exp(axis.scale(epsilon)).unwrap(),
                    )
                    .unwrap()
                } else {
                    Pose::new(
                        pose.position_world().add(axis.scale(epsilon)),
                        pose.orientation(),
                    )
                    .unwrap()
                };
                rotor
                    .loads(
                        RigidBodyState::new(p, Vec3::ZERO, momentum).unwrap(),
                        0.0,
                        || false,
                    )
                    .unwrap()
                    .bearing_potential_j
            };
            let gradient = (potential(1e-6) - potential(-1e-6)) / 2e-6;
            let expected = if angular {
                moment[i]
            } else {
                components(load.wrench.force_world)[i]
            };
            near(-gradient, expected, 2e-8);
        }
    }
}

#[test]
fn g4_unsupported_bearings_and_late_cancellation_publish_no_candidate() {
    for (stiffness, damping) in [
        ([[1.0, 0.1], [0.0, 1.0]], [[0.0; 2]; 2]),
        ([[1.0, 2.0], [2.0, 1.0]], [[0.0; 2]; 2]),
        ([[1.0, 0.0], [0.0, -1.0]], [[0.0; 2]; 2]),
        ([[1.0, 0.0], [0.0, 1.0]], [[0.1, 0.0], [0.0, 0.0]]),
    ] {
        let good = rotor(&[]);
        let bad = RadialBearing {
            station_m: 0.0,
            stiffness_n_per_m: stiffness,
            damping_ns_per_m: damping,
        };
        assert!(
            ElasticEccentricRotor::new(good.source_inertia(), Se3::identity(), &[bad], || false)
                .is_err()
        );
    }
    let rotor = rotor(&[bearing()]);
    let initial = state(Vec3::new(ECCENTRICITY, 0.0, 0.0), 0.0, Vec3::ZERO, 3.0);
    let polls = Cell::new(0);
    let expected = rotor
        .step(integrator(), initial, 0.1, 0.01, || {
            polls.set(polls.get() + 1);
            false
        })
        .unwrap();
    let refusal_at = polls.get();
    polls.set(0);
    assert_eq!(
        rotor.step(integrator(), initial, 0.1, 0.01, || {
            polls.set(polls.get() + 1);
            polls.get() == refusal_at
        }),
        Err(ElasticEccentricRotorError::Cancelled)
    );
    assert_eq!(
        rotor
            .step(integrator(), initial, 0.1, 0.01, || false)
            .unwrap(),
        expected
    );
    assert!(
        rotor
            .step(integrator(), initial, f64::NAN, 0.01, || false)
            .is_err()
    );
    assert!(
        rotor
            .step(integrator(), initial, 0.0, 0.0, || false)
            .is_err()
    );
}
