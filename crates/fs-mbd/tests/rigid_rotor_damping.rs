//! G1/G3/G4: physical bearing heat and reciprocal damping on the original body.
use fs_ga::{Mat3, Se3, So3, So3Tangent, Vec3 as GaVec3};
use fs_mbd::articulated::SpatialInertia;
use fs_mbd::rotordynamics::rigid::{
    DampedEccentricRotor, EccentricRotorDampingControl, EccentricRotorDampingError, RadialBearing,
};
use fs_mbd::{Gravity, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3};
use std::cell::Cell;

const MASS: f64 = 3.0;
const E: f64 = 0.25;
const POLAR: f64 = 0.8;

fn bearing() -> RadialBearing {
    RadialBearing {
        station_m: 0.0,
        stiffness_n_per_m: [[30.0, 5.0], [5.0, 18.0]],
        damping_ns_per_m: [[0.8, 0.1], [0.1, 0.5]],
    }
}
fn rotor(bearings: &[RadialBearing]) -> DampedEccentricRotor {
    DampedEccentricRotor::new(
        SpatialInertia::new(
            MASS,
            GaVec3::new(E, 0.0, 0.0),
            Mat3 {
                m: [1.2, 0.0, 0.0, 0.0, 1.2, 0.0, 0.0, 0.0, POLAR],
            },
        )
        .unwrap(),
        Se3::identity(),
        bearings,
        || false,
    )
    .unwrap()
}
fn integrator() -> RigidBodyIntegrator {
    RigidBodyIntegrator::new(Gravity::ZERO)
}
fn control() -> EccentricRotorDampingControl {
    EccentricRotorDampingControl {
        relative_tolerance: 1e-11,
        maximum_modal_step: 1.0,
        linear_impulse_tolerance_ns: 1e-12,
        angular_impulse_tolerance_nms: 1e-12,
        linear_speed_tolerance_m_s: 1e-11,
        angular_speed_tolerance_rad_s: 1e-11,
        source_energy_tolerance_j: 1e-10,
        total_energy_tolerance_j: 0.01,
    }
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
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance * (1.0 + b.abs()), "{a:e} != {b:e}");
}
fn rates(rotor: &DampedEccentricRotor, state: RigidBodyState) -> [f64; 6] {
    let props = rotor.elastic_component().mass_properties();
    let v = state.center_of_mass_velocity_world(props).unwrap();
    let w = props.angular_velocity_body_checked(state.angular_momentum_body()).unwrap();
    [v.x, v.y, v.z, w.x, w.y, w.z]
}

#[test]
fn g1_eccentric_bearing_damping_brakes_spin_and_returns_exact_midpoint_heat() {
    let coefficient = 4.0;
    let rotor = rotor(&[RadialBearing {
        stiffness_n_per_m: [[0.0; 2]; 2],
        damping_ns_per_m: [[0.0, 0.0], [0.0, coefficient]],
        ..bearing()
    }]);
    let initial = state(Vec3::new(E, 0.0, 0.0), 0.0, Vec3::ZERO, 3.0);
    let h = 0.1;
    // The material shaft velocity is u=vy-e*omega, so its original mobility
    // includes both COM translation and reciprocal central angular acceleration.
    let u0 = -E * 3.0;
    let rate = coefficient * (1.0 / MASS + E * E / POLAR);
    let u1 = u0 * (1.0 - 0.5 * h * rate) / (1.0 + 0.5 * h * rate);
    let midpoint = 0.5 * (u0 + u1);
    let impulse = -h * coefficient * midpoint;
    let heat = h * coefficient * midpoint * midpoint;
    let step = rotor.relax_damping(initial, h, control(), || false).unwrap();
    let after = step.impulse.state_after;
    assert_eq!(after.pose(), initial.pose());
    near(after.linear_momentum_world().y, impulse, 2e-14);
    near(after.angular_momentum_body().z, POLAR * 3.0 - E * impulse, 2e-14);
    near(rates(&rotor, after)[1] - E * rates(&rotor, after)[5], u1, 2e-14);
    near(step.decay_rate_bound_per_s, rate, 1e-13);
    near(step.dissipated_energy_j, heat, 2e-14);
    near(step.bearing_heat_j[0], heat, 2e-14);
    near(step.impulse.work.kinetic_energy_change_j, -heat, 2e-14);
    near(step.kinetic_heat_residual_j, 0.0, 2e-14);
    near(step.impulse_heat_residual_j, 0.0, 2e-14);
    assert!(step.native.is_some() && step.dissipated_energy_j > 0.0);
    assert!(after.angular_momentum_body().z < initial.angular_momentum_body().z);

    let load = rotor.loads(initial, 0.0, || false).unwrap();
    near(load.wrench.force_world.y, -coefficient * u0, 2e-14);
    near(load.wrench.torque_body.z, E * coefficient * u0, 2e-14);
    near(load.dissipation_power_w, coefficient * u0 * u0, 2e-14);
    let motion = rotor.prepare_motion(integrator(), initial, 0.0, || false)
        .unwrap().point_motion(Vec3::new(-E, 0.0, 0.0)).unwrap();
    near(motion.point_acceleration_world.x, E * 9.0, 2e-14);
    near(motion.point_acceleration_world.y, -rate * u0, 2e-14);
    near(motion.angular_acceleration_body.z, E * coefficient * u0 / POLAR, 2e-14);

    // The undamped null direction remains free: no artificial damping floor.
    let null = state(Vec3::new(E, 0.0, 0.0), 0.0, Vec3::new(2.0, E * 3.0, 0.0), 3.0);
    let zero = rotor.relax_damping(null, h, control(), || false).unwrap();
    near(zero.dissipated_energy_j, 0.0, 1e-27);
    for (actual, expected) in rates(&rotor, zero.impulse.state_after).into_iter().zip(rates(&rotor, null)) {
        near(actual, expected, 2e-14);
    }
}

// Independent pivoted dense solve for the physical test oracle only.
fn solve(matrix: [[f64; 6]; 6], rhs: [f64; 6]) -> [f64; 6] {
    let mut a: [[f64; 7]; 6] = core::array::from_fn(|i| {
        core::array::from_fn(|j| if j == 6 { rhs[i] } else { matrix[i][j] })
    });
    for k in 0..6 {
        let pivot = (k..6).max_by(|&i, &j| a[i][k].abs().total_cmp(&a[j][k].abs())).unwrap();
        a.swap(k, pivot);
        let diagonal = a[k][k];
        assert!(diagonal.abs() > 1e-8);
        for j in k..7 { a[k][j] /= diagonal; }
        let row = a[k];
        for i in 0..6 {
            if i != k {
                let scale = a[i][k];
                for j in k..7 { a[i][j] -= scale * row[j]; }
            }
        }
    }
    core::array::from_fn(|i| a[i][6])
}

#[test]
fn g1_full_tensor_and_cross_bearings_solve_one_conjugate_source_at_the_actual_pose() {
    let tensor = [1.3, 0.1, 0.05, 0.1, 1.1, 0.08, 0.05, 0.08, 0.9];
    let offset = GaVec3::new(0.2, -0.1, 0.15);
    let frame = Se3::from_parts(
        So3::exp(So3Tangent::new(GaVec3::new(0.2, -0.4, 0.6))).unwrap(),
        GaVec3::new(1.0, -2.0, 0.5),
    ).unwrap();
    let bearings = [
        RadialBearing { station_m: -0.4, damping_ns_per_m: [[4.0, 1.0], [1.0, 2.0]], ..bearing() },
        RadialBearing { station_m: 0.7, damping_ns_per_m: [[3.0, -0.5], [-0.5, 1.0]], ..bearing() },
    ];
    let rotor = DampedEccentricRotor::new(
        SpatialInertia::new(MASS, offset, Mat3 { m: tensor }).unwrap(),
        frame, &bearings, || false,
    ).unwrap();
    let original = rotor.elastic_component();
    let pose = Pose::new(
        original.reference_pose().position_world().add(Vec3::new(0.1, 0.03, -0.02)),
        original.reference_pose().orientation().right_exp(Vec3::new(0.3, -0.2, 0.4)).unwrap(),
    ).unwrap();
    let before = [0.6, -0.3, 0.2, 0.4, -0.7, 2.0];
    let initial = RigidBodyState::new(
        pose, Vec3::new(before[0], before[1], before[2]).scale(MASS),
        original.mass_properties().angular_momentum_body_checked(Vec3::new(before[3], before[4], before[5])).unwrap(),
    ).unwrap();
    let basis = [Vec3::new(1.0, 0.0, 0.0), Vec3::new(0.0, 1.0, 0.0), Vec3::new(0.0, 0.0, 1.0)];
    let axes = [basis[0], basis[1]].map(|e| original.reference_pose().orientation().rotate_body_to_world(e));
    // Obtain each row from forward world-point velocities, independently of
    // production's reciprocal arm x inverse-rotated-axis construction.
    let jacobians: [[[f64; 6]; 2]; 2] = core::array::from_fn(|b| {
        let arm = Vec3::new(-offset.x, -offset.y, bearings[b].station_m - offset.z);
        core::array::from_fn(|a| core::array::from_fn(|k| {
            let world = if k < 3 { basis[k] } else { pose.orientation().rotate_body_to_world(basis[k - 3].cross(arm)) };
            axes[a].dot(world)
        }))
    });
    let mass: [[f64; 6]; 6] = core::array::from_fn(|i| core::array::from_fn(|j| {
        if i < 3 && j < 3 { if i == j { MASS } else { 0.0 } }
        else if i >= 3 && j >= 3 { tensor[3 * (i - 3) + j - 3] } else { 0.0 }
    }));
    let mut damping = [[0.0; 6]; 6];
    for (b, j) in jacobians.iter().enumerate() {
        assert_eq!(rotor.damping_coefficients(b), Some(bearings[b].damping_ns_per_m));
        for row in 0..6 { for col in 0..6 { for a in 0..2 { for d in 0..2 {
            damping[row][col] += j[a][row] * bearings[b].damping_ns_per_m[a][d] * j[d][col];
        }}}}
    }
    let h = 0.02;
    let equation = core::array::from_fn(|i| core::array::from_fn(|j| mass[i][j] + 0.5 * h * damping[i][j]));
    let rhs = core::array::from_fn(|i| -h * (0..6).map(|j| damping[i][j] * before[j]).sum::<f64>());
    let delta = solve(equation, rhs);
    let expected: [f64; 6] = core::array::from_fn(|i| before[i] + delta[i]);
    let step = rotor.relax_damping(initial, h, control(), || false).unwrap();
    let actual = rates(&rotor, step.impulse.state_after);
    assert_eq!(step.impulse.state_after.pose(), initial.pose());
    for i in 0..6 { near(actual[i], expected[i], 3e-13); }
    let midpoint: [f64; 6] = core::array::from_fn(|i| 0.5 * (before[i] + expected[i]));
    for b in 0..2 {
        let u: [f64; 2] = core::array::from_fn(|a| (0..6).map(|k| jacobians[b][a][k] * midpoint[k]).sum());
        let c = bearings[b].damping_ns_per_m;
        let heat = h * (u[0] * (c[0][0] * u[0] + c[0][1] * u[1]) + u[1] * (c[1][0] * u[0] + c[1][1] * u[1]));
        assert!(heat > 0.0);
        near(step.bearing_heat_j[b], heat, 3e-13);
    }
    let trace: f64 = (0..6).map(|column| solve(mass, core::array::from_fn(|row| damping[row][column]))[column]).sum();
    near(step.decay_rate_bound_per_s, trace * (1.0 + 256.0 * f64::EPSILON), 3e-13);
    near(step.kinetic_heat_residual_j, 0.0, 3e-13);
    let load = rotor.loads(initial, 0.0, || false).unwrap();
    let elastic = original.loads(initial, 0.0, || false).unwrap();
    let f = load.wrench.force_world.sub(elastic.wrench.force_world);
    let t = load.wrench.torque_body.sub(elastic.wrench.torque_body);
    let power = f.dot(Vec3::new(before[0], before[1], before[2])) + t.dot(Vec3::new(before[3], before[4], before[5]));
    near(power + load.dissipation_power_w, 0.0, 3e-13);
}

// Independent physical COM equations; the seventh coordinate integrates the
// continuous constitutive heat. Production uses no RK4 or duplicate body state.
fn derivative(y: [f64; 7], torque: f64) -> [f64; 7] {
    let (s, c) = y[2].sin_cos();
    let q = [y[0] - E * c, y[1] - E * s];
    let u = [y[3] + E * y[5] * s, y[4] - E * y[5] * c];
    let force = [-30.0 * q[0] - 5.0 * q[1] - 0.8 * u[0] - 0.1 * u[1], -5.0 * q[0] - 18.0 * q[1] - 0.1 * u[0] - 0.5 * u[1]];
    [y[3], y[4], y[5], force[0] / MASS, force[1] / MASS,
     (torque - E * (-s * force[0] + c * force[1])) / POLAR,
     u[0] * (0.8 * u[0] + 0.1 * u[1]) + u[1] * (0.1 * u[0] + 0.5 * u[1])]
}
fn reference(torque: f64) -> [f64; 7] {
    let mut y = [E, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0];
    let h = 2.0 / 8000.0;
    for _ in 0..8000 {
        let k1 = derivative(y, torque);
        let k2 = derivative(core::array::from_fn(|i| y[i] + 0.5 * h * k1[i]), torque);
        let k3 = derivative(core::array::from_fn(|i| y[i] + 0.5 * h * k2[i]), torque);
        let k4 = derivative(core::array::from_fn(|i| y[i] + h * k3[i]), torque);
        y = core::array::from_fn(|i| y[i] + h * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]) / 6.0);
    }
    y
}
fn trajectory(steps: usize, torque: f64, oracle: [f64; 7]) -> (f64, f64, f64) {
    let rotor = rotor(&[bearing()]);
    let mut current = state(Vec3::new(E, 0.0, 0.0), 0.0, Vec3::ZERO, 3.0);
    let mut heat = 0.0;
    let mut defect = 0.0;
    for _ in 0..steps {
        let step = rotor.step(integrator(), current, torque, 2.0 / steps as f64, control(), || false).unwrap();
        for source in &step.sources { assert!(source.dissipated_energy_j >= 0.0); }
        heat += step.dissipated_energy_j;
        defect += step.total_energy_residual_j;
        current = step.state_after;
    }
    let p = current.pose().position_world();
    let v = rates(&rotor, current);
    let axis = current.pose().orientation().rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
    let error = (p.x - oracle[0]).abs() + (p.y - oracle[1]).abs()
        + (v[0] - oracle[3]).abs() + (v[1] - oracle[4]).abs() + (v[5] - oracle[5]).abs()
        + (axis.x - oracle[2].cos()).abs() + (axis.y - oracle[2].sin()).abs();
    (error, (heat - oracle[6]).abs(), defect.abs())
}

#[test]
fn g3_reciprocal_damped_orbit_and_independent_heat_converge_at_second_order() {
    for torque in [0.0, 0.12] {
        let oracle = reference(torque);
        assert!(oracle[6] > 0.1);
        let coarse = trajectory(50, torque, oracle);
        let medium = trajectory(100, torque, oracle);
        let fine = trajectory(200, torque, oracle);
        assert!(coarse.0 / medium.0 > 3.6 && medium.0 / fine.0 > 3.6, "{coarse:?} {medium:?} {fine:?}");
        assert!(coarse.1 / medium.1 > 3.6 && medium.1 / fine.1 > 3.6, "{coarse:?} {medium:?} {fine:?}");
        assert!(fine.0 < 2e-4 && fine.1 < 4e-5 && fine.2 < 3e-4, "{fine:?}");
    }
}

#[test]
fn g4_original_passivity_and_modal_resolution_refuse_unphysical_sources() {
    let valid = rotor(&[]);
    for c in [ [[1.0, 0.1], [0.0, 1.0]], [[1.0, 2.0], [2.0, 1.0]], [[1.0, 0.0], [0.0, -1.0]], [[f64::NAN, 0.0], [0.0, 1.0]] ] {
        assert!(matches!(DampedEccentricRotor::new(valid.elastic_component().source_inertia(), Se3::identity(), &[RadialBearing { damping_ns_per_m: c, ..bearing() }], || false), Err(EccentricRotorDampingError::Material(_))));
    }
    let rotor = rotor(&[RadialBearing { stiffness_n_per_m: [[0.0; 2]; 2], damping_ns_per_m: [[0.0, 0.0], [0.0, 4.0]], ..bearing() }]);
    let initial = state(Vec3::new(E, 0.0, 0.0), 0.0, Vec3::ZERO, 3.0);
    let resolved = EccentricRotorDampingControl { maximum_modal_step: 2.0, ..control() };
    assert!(matches!(rotor.relax_damping(initial, 2.0, resolved, || false), Err(EccentricRotorDampingError::ModalStep { .. })));
    let accepted = rotor.relax_damping(initial, 1.0, resolved, || false).unwrap();
    let v = rates(&rotor, accepted.impulse.state_after);
    assert!(v[1] - E * v[5] < 0.0); // Resolved passive decay cannot reverse this scalar port.
    let invalid = EccentricRotorDampingControl { maximum_modal_step: 2.01, ..control() };
    assert!(rotor.relax_damping(initial, 0.01, invalid, || false).is_err());
    let identity = rotor.relax_damping(initial, 0.0, control(), || false).unwrap();
    assert_eq!(identity.impulse.state_after, initial);
    assert_eq!(identity.dissipated_energy_j, 0.0);
    assert!(identity.native.is_none());
}

#[test]
fn g4_zero_damping_reduces_to_original_elastic_motion_and_late_cancel_is_atomic() {
    let zero = rotor(&[RadialBearing { damping_ns_per_m: [[0.0; 2]; 2], ..bearing() }]);
    let initial = state(Vec3::new(E, 0.0, 0.0), 0.0, Vec3::ZERO, 3.0);
    let elastic = zero.elastic_component().step(integrator(), initial, 0.1, 0.01, || false).unwrap();
    let complete = zero.step(integrator(), initial, 0.1, 0.01, control(), || false).unwrap();
    assert_eq!(complete.state_after, elastic.native.state_after);
    assert_eq!(complete.elastic, elastic);
    assert_eq!(complete.dissipated_energy_j, 0.0);
    assert!(complete.sources.iter().all(|s| s.native.is_none()));

    let rotor = rotor(&[bearing()]);
    let polls = Cell::new(0);
    let expected = rotor.step(integrator(), initial, 0.1, 0.01, control(), || { polls.set(polls.get() + 1); false }).unwrap();
    let refusal_at = polls.get();
    polls.set(0);
    assert_eq!(rotor.step(integrator(), initial, 0.1, 0.01, control(), || { polls.set(polls.get() + 1); polls.get() == refusal_at }), Err(EccentricRotorDampingError::Cancelled));
    assert_eq!(rotor.step(integrator(), initial, 0.1, 0.01, control(), || false).unwrap(), expected);
    assert_eq!(rotor.damping_coefficients(0), Some(bearing().damping_ns_per_m));
}
