//! Finite revolute actuation against independent rotor momentum references.

use fs_mbd::rigid_rattle::distance::graph::{
    AngularDrive, DistanceDrive, DistanceGraphControls, DistanceGraphError, DistanceLink,
    DrivenJoint, JointReactionImpulse, PrismaticJoint, RevoluteMotorJoint,
    advance_driven_joint_graph, advance_spatial_driven_joint_graph,
};
use fs_mbd::rigid_rattle::distance::{DistanceStepControls, DistanceStepError};
use fs_mbd::spatial_loads::{PreparedSpatialLoads, WorldPointForce};
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

const X: Vec3 = Vec3::new(1.0, 0.0, 0.0);
const Y: Vec3 = Vec3::new(0.0, 1.0, 0.0);
const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tol={tolerance:e}"
    );
}
fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    near(a.x, b.x, tolerance);
    near(a.y, b.y, tolerance);
    near(a.z, b.z, tolerance);
}
fn free() -> RigidBodyIntegrator {
    RigidBodyIntegrator::new(Gravity::ZERO)
}
fn controls() -> DistanceGraphControls {
    DistanceGraphControls {
        step: DistanceStepControls {
            maximum_impulse_ns: 20.0,
            maximum_rotation_increment_rad: 0.25,
            maximum_position_iterations: 50,
            maximum_backtracks: 20,
        },
        position_tolerance_m: 5e-13,
        velocity_tolerance_mps: 5e-13,
        minimum_scaled_pivot: 1e-10,
        maximum_bodies: 64,
        maximum_constraints: 32,
        maximum_trial_evaluations: 200,
        maximum_live_bytes: 2_000_000,
    }
}
fn spherical(mass: f64, moment: f64) -> MassProperties {
    MassProperties::new(mass, Vec3::ZERO, Vec3::new(moment, moment, moment)).unwrap()
}
fn state(p: MassProperties, x: Vec3, q: UnitQuaternion, v: Vec3, w_body: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(x, q).unwrap(),
        v.scale(p.mass()),
        p.angular_momentum_body_checked(w_body).unwrap(),
    )
    .unwrap()
}
fn motor(bodies: [usize; 2], rate: f64, acceleration: f64) -> DrivenJoint {
    DrivenJoint::RevoluteMotor {
        joint: RevoluteMotorJoint::new(bodies, [Vec3::ZERO; 2], [UnitQuaternion::IDENTITY; 2], 0.8)
            .unwrap(),
        drive: AngularDrive {
            angular_rate_rad_s: rate,
            angular_acceleration_rad_s2: acceleration,
        },
    }
}
fn reverse(s: RigidBodyState) -> RigidBodyState {
    RigidBodyState::new(
        s.pose(),
        s.linear_momentum_world().scale(-1.0),
        s.angular_momentum_body().scale(-1.0),
    )
    .unwrap()
}
fn reverse_joint(j: DrivenJoint) -> DrivenJoint {
    match j {
        DrivenJoint::RevoluteMotor { joint, drive } => DrivenJoint::RevoluteMotor {
            joint,
            drive: AngularDrive {
                angular_rate_rad_s: -drive.angular_rate_rad_s,
                ..drive
            },
        },
        _ => j,
    }
}
fn compose(a: UnitQuaternion, b: UnitQuaternion) -> UnitQuaternion {
    let [w, x, y, z] = a.components();
    let [s, u, v, t] = b.components();
    UnitQuaternion::new(
        w * s - x * u - y * v - z * t,
        w * u + x * s + y * t - z * v,
        w * v - x * t + y * s + z * u,
        w * t + x * v - y * u + z * s,
    )
    .unwrap()
}
fn state_error(a: RigidBodyState, b: RigidBodyState) -> f64 {
    let mut error = norm(a.pose().position_world().sub(b.pose().position_world()))
        .max(norm(
            a.linear_momentum_world().sub(b.linear_momentum_world()),
        ))
        .max(norm(
            a.angular_momentum_body().sub(b.angular_momentum_body()),
        ));
    for axis in [X, Y, Z] {
        error = error.max(norm(
            a.pose()
                .orientation()
                .rotate_body_to_world(axis)
                .sub(b.pose().orientation().rotate_body_to_world(axis)),
        ));
    }
    error
}
fn totals(s: &[RigidBodyState], p: &[MassProperties]) -> (Vec3, Vec3, f64) {
    s.iter().zip(p).fold(
        (Vec3::ZERO, Vec3::ZERO, 0.0),
        |(linear, angular, energy), (&s, &p)| {
            let d = free().diagnostics(s, p).unwrap();
            (
                linear.add(d.linear_momentum_world),
                angular.add(d.angular_momentum_world),
                energy + d.mechanical_energy,
            )
        },
    )
}

#[test]
fn g1_loaded_three_rotor_chain_has_analytic_recoil_signed_work_and_multi_turn_continuation() {
    // H = sum(I_i omega_i), omega_1 - omega_0 = nu_1 and
    // omega_2 - omega_1 = nu_2 determine every angular acceleration directly.
    // No constrained solver is used in this reference, including motor torques.
    let mass = [1.0, 2.0, 4.0];
    let moment = [2.0, 3.0, 5.0];
    let p = core::array::from_fn::<_, 3, _>(|i| spherical(mass[i], moment[i]));
    let initial_rate = [-0.69, 0.51, 0.11];
    let acceleration = [-0.07, 0.03, -0.03];
    let torque = [0.2, -0.1, -0.3];
    let force = [1.2, -0.7, 0.4];
    let boost = Vec3::new(0.2, -0.1, 0.05);
    let initial = core::array::from_fn::<_, 3, _>(|i| {
        state(
            p[i],
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            boost,
            Z.scale(initial_rate[i]),
        )
    });
    let loads = core::array::from_fn::<_, 3, _>(|i| Wrench {
        force_world: X.scale(force[i]),
        torque_body: Z.scale(torque[i]),
    });
    let mut joints = vec![motor([0, 1], 1.2, 0.1), motor([1, 2], -0.4, -0.06)];
    let mut s = initial.to_vec();
    let before = totals(&s, &p);
    let mut work = [0.0; 2];
    let mut defect = 0.0;
    let h = 0.02;
    for k in 1..=300 {
        let out =
            advance_driven_joint_graph(&s, &p, &joints, free(), &loads, h, controls(), || false)
                .unwrap();
        for (i, receipt) in out.joints.iter().enumerate() {
            assert_eq!(receipt.rows.row_count, 6);
            assert_eq!(receipt.rows.first_row, 6 * i);
            let command = receipt.angular.unwrap();
            assert!(command.maximum_excursion_rad < 0.04);
            assert!(command.displacement_representation_residual_rad.abs() < 1e-15);
            assert!(command.transverse_frame_residual_rad < 1e-15);
            work[i] += receipt.actuator_work_j;
        }
        defect += out.constraint_work_defect_j;
        near(
            out.constraint_work_j,
            out.actuator_work_j + out.constraint_work_defect_j,
            2e-16,
        );
        joints = out.joints.iter().map(|j| j.joint_after).collect();
        s = out.step.states_after;
        let t = f64::from(k) * h;
        for i in 0..3 {
            let expected = state(
                p[i],
                boost.scale(t).add(X.scale(0.5 * (0.9 / 7.0) * t * t)),
                UnitQuaternion::from_axis_angle(
                    Z,
                    initial_rate[i] * t + 0.5 * acceleration[i] * t * t,
                )
                .unwrap(),
                boost.add(X.scale((0.9 / 7.0) * t)),
                Z.scale(initial_rate[i] + acceleration[i] * t),
            );
            assert!(
                state_error(s[i], expected) < 2e-10,
                "step {k}, body {i}: {}",
                state_error(s[i], expected)
            );
        }
    }
    // The first actuator passes a full revolution, including canonical frame
    // sign changes; the second recovers energy while the first supplies it.
    near(work[0], 0.34 * 9.0, 4e-10);
    near(work[1], 0.15 * -3.48, 4e-10);
    assert!(work[0] > 0.0 && work[1] < 0.0);
    let after = totals(&s, &p);
    vector(after.0.sub(before.0), X.scale(0.9 * 6.0), 4e-11);
    // The common translating center contributes its known orbital momentum.
    let x = s[0].pose().position_world();
    let intrinsic = after.1.sub(x.cross(after.0));
    vector(intrinsic, Z.scale(0.7 - 0.2 * 6.0), 4e-11);
    let source = (0..3)
        .map(|i| torque[i] * (initial_rate[i] * 6.0 + 18.0 * acceleration[i]))
        .sum::<f64>()
        + 0.9 * x.x;
    near(
        after.2 - before.2,
        work.iter().sum::<f64>() + source + defect,
        4e-10,
    );
    assert!(defect.abs() < 4e-10);

    s = s.into_iter().map(reverse).collect();
    joints = joints.into_iter().map(reverse_joint).collect();
    for _ in 0..300 {
        let out =
            advance_driven_joint_graph(&s, &p, &joints, free(), &loads, h, controls(), || false)
                .unwrap();
        joints = out.joints.iter().map(|j| j.joint_after).collect();
        s = out.step.states_after;
    }
    for i in 0..3 {
        assert!(state_error(reverse(s[i]), initial[i]) < 5e-10);
    }
}

// Independent 3x3 tensor arithmetic and raw-quaternion ODE. The production
// nonlinear solver, midpoint/free step and angular Jacobian are not used here.
type Matrix = [[f64; 3]; 3];
fn components(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn from_components(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn mat_vec(a: Matrix, x: [f64; 3]) -> [f64; 3] {
    a.map(|r| (0..3).map(|k| r[k] * x[k]).sum())
}
fn transpose(a: Matrix) -> Matrix {
    core::array::from_fn(|i| core::array::from_fn(|j| a[j][i]))
}
fn product(a: Matrix, b: Matrix) -> Matrix {
    core::array::from_fn(|i| core::array::from_fn(|j| (0..3).map(|k| a[i][k] * b[k][j]).sum()))
}
fn solve(mut a: Matrix, mut b: [f64; 3]) -> [f64; 3] {
    for k in 0..3 {
        let pivot = (k..3)
            .max_by(|&i, &j| a[i][k].abs().total_cmp(&a[j][k].abs()))
            .unwrap();
        a.swap(k, pivot);
        b.swap(k, pivot);
        for i in k + 1..3 {
            let scale = a[i][k] / a[k][k];
            for j in k..3 {
                a[i][j] -= scale * a[k][j];
            }
            b[i] -= scale * b[k];
        }
    }
    let mut x = [0.0; 3];
    for i in (0..3).rev() {
        x[i] = (b[i] - (i + 1..3).map(|j| a[i][j] * x[j]).sum::<f64>()) / a[i][i];
    }
    x
}
fn rotation(mut q: [f64; 4]) -> Matrix {
    let norm = q.iter().map(|x| x * x).sum::<f64>().sqrt();
    for x in &mut q {
        *x /= norm;
    }
    let [w, x, y, z] = q;
    [
        [
            1.0 - 2.0 * (y * y + z * z),
            2.0 * (x * y - w * z),
            2.0 * (x * z + w * y),
        ],
        [
            2.0 * (x * y + w * z),
            1.0 - 2.0 * (x * x + z * z),
            2.0 * (y * z - w * x),
        ],
        [
            2.0 * (x * z - w * y),
            2.0 * (y * z + w * x),
            1.0 - 2.0 * (x * x + y * y),
        ],
    ]
}
fn tensor(p: MassProperties) -> Matrix {
    transpose([X, Y, Z].map(|axis| components(p.angular_momentum_body_checked(axis).unwrap())))
}
struct TensorRotors {
    inertia: [Matrix; 2],
    momentum: [f64; 3],
    rate: f64,
    acceleration: f64,
}
impl TensorRotors {
    fn velocity(&self, t: f64, q: [f64; 4]) -> [[f64; 3]; 2] {
        let angle = self.rate * t + 0.5 * self.acceleration * t * t;
        let rate = self.rate + self.acceleration * t;
        let (sin, cos) = angle.sin_cos();
        let relative = [[cos, -sin, 0.0], [sin, cos, 0.0], [0.0, 0.0, 1.0]];
        let turned = product(product(relative, self.inertia[1]), transpose(relative));
        let combined = core::array::from_fn(|i| {
            core::array::from_fn(|j| self.inertia[0][i][j] + turned[i][j])
        });
        let bias = mat_vec(relative, mat_vec(self.inertia[1], [0.0, 0.0, rate]));
        let local_h = mat_vec(transpose(rotation(q)), self.momentum);
        let first = solve(combined, core::array::from_fn(|i| local_h[i] - bias[i]));
        let mut second = mat_vec(transpose(relative), first);
        second[2] += rate;
        [first, second]
    }
    fn rhs(&self, t: f64, q: [f64; 4]) -> [f64; 4] {
        let [w, x, y, z] = q;
        let [a, b, c] = self.velocity(t, q)[0];
        [
            -0.5 * (x * a + y * b + z * c),
            0.5 * (w * a + y * c - z * b),
            0.5 * (w * b + z * a - x * c),
            0.5 * (w * c + x * b - y * a),
        ]
    }
    fn reference(&self, q_initial: UnitQuaternion, duration: f64, count: u32) -> UnitQuaternion {
        let mut q = q_initial.components();
        let h = duration / f64::from(count);
        for k in 0..count {
            let t = f64::from(k) * h;
            let a = self.rhs(t, q);
            let b = self.rhs(t + 0.5 * h, core::array::from_fn(|i| q[i] + 0.5 * h * a[i]));
            let c = self.rhs(t + 0.5 * h, core::array::from_fn(|i| q[i] + 0.5 * h * b[i]));
            let d = self.rhs(t + h, core::array::from_fn(|i| q[i] + h * c[i]));
            q = core::array::from_fn(|i| {
                q[i] + (h / 6.0) * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i])
            });
            let length = q.iter().map(|x| x * x).sum::<f64>().sqrt();
            for x in &mut q {
                *x /= length;
            }
        }
        UnitQuaternion::new(q[0], q[1], q[2], q[3]).unwrap()
    }
    fn states(
        &self,
        t: f64,
        q: UnitQuaternion,
        p: [MassProperties; 2],
        boost: Vec3,
    ) -> [RigidBodyState; 2] {
        let rate = self.velocity(t, q.components());
        let angle = self.rate * t + 0.5 * self.acceleration * t * t;
        let [w, x, y, z] = q.components();
        let (s, c) = (0.5 * angle).sin_cos();
        let second =
            UnitQuaternion::new(w * c - z * s, x * c + y * s, y * c - x * s, z * c + w * s)
                .unwrap();
        let orientation = [q, second];
        core::array::from_fn(|i| {
            state(
                p[i],
                boost.scale(t),
                orientation[i],
                boost,
                from_components(rate[i]),
            )
        })
    }
}

#[test]
fn g1_full_tensor_motor_converges_to_independent_conserved_angular_momentum_ode() {
    let p = [
        MassProperties::from_principal_frame(
            2.0,
            Vec3::ZERO,
            Vec3::new(0.7, 1.1, 1.5),
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, -0.4, 0.3), 0.43).unwrap(),
        )
        .unwrap(),
        MassProperties::from_principal_frame(
            3.0,
            Vec3::ZERO,
            Vec3::new(0.9, 1.4, 1.9),
            UnitQuaternion::from_axis_angle(Vec3::new(-0.5, 1.0, 0.6), -0.37).unwrap(),
        )
        .unwrap(),
    ];
    let model = TensorRotors {
        inertia: p.map(tensor),
        momentum: [0.3, 1.1, -0.4],
        rate: 0.4,
        acceleration: 0.2,
    };
    let q0 = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 0.4), 0.3).unwrap();
    let boost = Vec3::new(0.1, -0.2, 0.3);
    let initial = model.states(0.0, q0, p, boost);
    let before = totals(&initial, &p);
    vector(before.1, from_components(model.momentum), 2e-15);
    let duration = 0.8;
    let expected = model.states(duration, model.reference(q0, duration, 8192), p, boost);
    let reference_work = totals(&expected, &p).2 - before.2;
    let mut errors = Vec::new();
    let mut work_errors = Vec::new();
    for count in [16, 32, 64] {
        let h = duration / f64::from(count);
        let mut s = initial.to_vec();
        let mut joints = vec![motor([0, 1], model.rate, model.acceleration)];
        let (mut work, mut defect, mut free_energy) = (0.0, 0.0, 0.0);
        for _ in 0..count {
            let out = advance_driven_joint_graph(
                &s,
                &p,
                &joints,
                free(),
                &[Wrench::ZERO; 2],
                h,
                controls(),
                || false,
            )
            .unwrap();
            work += out.actuator_work_j;
            defect += out.constraint_work_defect_j;
            for receipt in &out.step.free_steps {
                free_energy += receipt.diagnostics_after.mechanical_energy
                    - receipt.diagnostics_before.mechanical_energy;
            }
            vector(out.constraint_angular_impulse_world_nms, Vec3::ZERO, 2e-16);
            joints = out.joints.iter().map(|j| j.joint_after).collect();
            s = out.step.states_after;
        }
        let after = totals(&s, &p);
        vector(after.0, before.0, 2e-13);
        vector(after.1, before.1, 2e-11);
        near(after.2 - before.2, work + defect + free_energy, 2e-12);
        errors.push(
            (0..2)
                .map(|i| state_error(s[i], expected[i]))
                .fold(0.0_f64, f64::max),
        );
        work_errors.push((work - reference_work).abs());
        if count == 64 {
            s = s.into_iter().map(reverse).collect();
            joints = joints.into_iter().map(reverse_joint).collect();
            for _ in 0..count {
                let out = advance_driven_joint_graph(
                    &s,
                    &p,
                    &joints,
                    free(),
                    &[Wrench::ZERO; 2],
                    h,
                    controls(),
                    || false,
                )
                .unwrap();
                joints = out.joints.iter().map(|j| j.joint_after).collect();
                s = out.step.states_after;
            }
            for i in 0..2 {
                assert!(state_error(reverse(s[i]), initial[i]) < 5e-10);
            }
        }
    }
    for pair in errors.windows(2) {
        assert!(
            (3.6..4.4).contains(&(pair[0] / pair[1])),
            "state refinement: {errors:?}"
        );
    }
    for pair in work_errors.windows(2) {
        assert!(
            (3.6..4.4).contains(&(pair[0] / pair[1])),
            "work refinement: {work_errors:?}"
        );
    }
    assert!(errors[2] < 2e-5, "{errors:?}");
    assert!(work_errors[2] < 2e-6, "{work_errors:?}");
}

fn spatial_mechanism() -> ([MassProperties; 3], [RigidBodyState; 3], [DrivenJoint; 3]) {
    let p = core::array::from_fn(|i| {
        MassProperties::from_principal_frame(
            1.3 + i as f64,
            Vec3::ZERO,
            Vec3::new(0.7, 1.1, 1.5).scale(1.0 + 0.2 * i as f64),
            UnitQuaternion::from_axis_angle(
                Vec3::new(0.4 + i as f64, -0.7, 1.0),
                0.35 + 0.1 * i as f64,
            )
            .unwrap(),
        )
        .unwrap()
    });
    let q = UnitQuaternion::from_axis_angle(Vec3::new(0.2, -0.5, 1.0), 0.4).unwrap();
    let x = [
        Vec3::new(-0.2, 0.1, -0.3),
        Vec3::new(0.4, -0.2, 0.1),
        Vec3::new(1.3, -0.2, 0.1),
    ]
    .map(|v| q.rotate_body_to_world(v));
    let w = Vec3::new(0.2, -0.1, 0.3);
    let world_w = q.rotate_body_to_world(w);
    let boost = Vec3::new(0.1, 0.2, -0.15);
    let s = core::array::from_fn(|i| state(p[i], x[i], q, boost.add(world_w.cross(x[i])), w));
    let motor_frame = UnitQuaternion::from_axis_angle(Vec3::new(0.4, -0.3, 1.0), 0.3).unwrap();
    let rail_frame = UnitQuaternion::from_axis_angle(Y, core::f64::consts::FRAC_PI_2).unwrap();
    let joints = [
        DrivenJoint::RevoluteMotor {
            joint: RevoluteMotorJoint::new(
                [0, 1],
                [Vec3::new(0.3, 0.1, 0.3), Vec3::new(-0.3, 0.4, -0.1)],
                [motor_frame; 2],
                0.8,
            )
            .unwrap(),
            drive: AngularDrive {
                angular_rate_rad_s: 0.0,
                angular_acceleration_rad_s2: 0.2,
            },
        },
        DrivenJoint::Prismatic(
            PrismaticJoint::new([1, 2], [Vec3::ZERO; 2], [rail_frame; 2], 0.9).unwrap(),
        ),
        DrivenJoint::Distance {
            link: DistanceLink {
                bodies: [0, 2],
                points_body_m: [Vec3::ZERO; 2],
                length_m: norm(x[2].sub(x[0])),
            },
            drive: DistanceDrive {
                range_rate_mps: 0.0,
                range_acceleration_mps2: 0.1,
            },
        },
    ];
    (p, s, joints)
}
fn spatial_loads(rotation: UnitQuaternion) -> [PreparedSpatialLoads; 3] {
    core::array::from_fn(|i| {
        let scale = 1.0 + 0.3 * i as f64;
        PreparedSpatialLoads::new(
            Wrench {
                force_world: rotation.rotate_body_to_world(Vec3::new(0.1, -0.2, 0.05).scale(scale)),
                torque_body: Vec3::new(0.03, -0.02, 0.04).scale(scale),
            },
            rotation.rotate_body_to_world(Vec3::new(-0.04, 0.08, 0.025).scale(scale)),
            &[WorldPointForce {
                arm_body: Vec3::new(0.2, 0.1, -0.12),
                force_world: rotation
                    .rotate_body_to_world(Vec3::new(0.18, -0.07, 0.09).scale(scale)),
            }],
            || false,
        )
        .unwrap()
    })
}
fn transformed(
    s: RigidBodyState,
    p: MassProperties,
    rotation: UnitQuaternion,
    boost: Vec3,
    time: f64,
) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(
            rotation
                .rotate_body_to_world(s.pose().position_world())
                .add(boost.scale(time)),
            compose(rotation, s.pose().orientation()),
        )
        .unwrap(),
        rotation
            .rotate_body_to_world(s.linear_momentum_world())
            .add(boost.scale(p.mass())),
        s.angular_momentum_body(),
    )
    .unwrap()
}

#[test]
fn g2_off_center_tensor_motor_rail_and_radial_drive_share_spatial_reactions_and_work() {
    let (p, initial, original) = spatial_mechanism();
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(-0.4, 0.2, 0.7), 1.1).unwrap();
    let boost = Vec3::new(0.7, -0.4, 0.5);
    let mut finals = Vec::new();
    let mut accounts = Vec::new();
    for (rotation, boost) in [(UnitQuaternion::IDENTITY, Vec3::ZERO), (rotation, boost)] {
        let loads = spatial_loads(rotation);
        let mut s = core::array::from_fn::<_, 3, _>(|i| {
            transformed(initial[i], p[i], rotation, boost, 0.0)
        })
        .to_vec();
        let mut joints = original.to_vec();
        let before = totals(&s, &p);
        let (mut linear, mut angular) = (Vec3::ZERO, Vec3::ZERO);
        let (mut source, mut motor_work, mut radial_work, mut constraint, mut defect) =
            (0.0, 0.0, 0.0, 0.0, 0.0);
        for _ in 0..60 {
            let out = advance_spatial_driven_joint_graph(
                &s,
                &p,
                &joints,
                free(),
                &loads,
                0.003,
                controls(),
                || false,
            )
            .unwrap();
            assert_eq!(
                out.joints
                    .iter()
                    .map(|j| j.rows.row_count)
                    .collect::<Vec<_>>(),
                [6, 5, 1]
            );
            assert_eq!(out.joints[1].actuator_work_j, 0.0);
            assert!(out.joints[0].angular.is_some());
            assert!(out.joints[1..].iter().all(|j| j.angular.is_none()));
            let DrivenJoint::RevoluteMotor { joint, drive } = joints[0] else {
                unreachable!()
            };
            let axis = joint.frames_body()[0].rotate_body_to_world(Z);
            let mut actual_motor_work = 0.0;
            for (events, rate) in [
                (&out.step.position_impulses, drive.angular_rate_rad_s),
                (
                    &out.step.velocity_impulses,
                    drive.angular_rate_rad_s + 0.003 * drive.angular_acceleration_rad_s2,
                ),
            ] {
                for event in &events[..3] {
                    let JointReactionImpulse::Point(event) = event else {
                        panic!("motor point closure")
                    };
                    vector(
                        event.first.point_kinematics_before.point_world,
                        event.second.point_kinematics_before.point_world,
                        1e-12,
                    );
                }
                for event in &events[3..6] {
                    let JointReactionImpulse::Couple(event) = event else {
                        panic!("motor pure couple")
                    };
                    vector(event.first.linear_increment_world, Vec3::ZERO, 0.0);
                    vector(event.second.linear_increment_world, Vec3::ZERO, 0.0);
                    actual_motor_work += event.couple_on_second_world_nms.dot(
                        event
                            .first
                            .state_before
                            .pose()
                            .orientation()
                            .rotate_body_to_world(axis),
                    ) * rate;
                }
            }
            near(actual_motor_work, out.joints[0].actuator_work_j, 3e-17);
            near(
                out.constraint_work_j,
                out.actuator_work_j + out.constraint_work_defect_j,
                2e-16,
            );
            angular = angular.add(out.constraint_angular_impulse_world_nms);
            motor_work += out.joints[0].actuator_work_j;
            radial_work += out.joints[2].actuator_work_j;
            constraint += out.constraint_work_defect_j;
            for (i, receipt) in out.step.free_steps.iter().enumerate() {
                assert_eq!(
                    *receipt,
                    free()
                        .step_spatial(receipt.step.state_before, p[i], &loads[i], 0.003)
                        .unwrap()
                );
                linear = linear.add(receipt.linear_impulse_world);
                angular = angular.add(receipt.angular_impulse_world);
                source += receipt.force_work_j + receipt.couple_work_j;
                defect += receipt.work_energy_residual_j;
            }
            joints = out.joints.iter().map(|j| j.joint_after).collect();
            s = out.step.states_after;
        }
        let after = totals(&s, &p);
        vector(after.0.sub(before.0), linear, 2e-11);
        vector(after.1.sub(before.1), angular, 3e-11);
        near(
            after.2 - before.2,
            source + motor_work + radial_work + constraint + defect,
            3e-11,
        );
        assert!(motor_work.abs() > 1e-5 && radial_work.abs() > 1e-5);
        finals.push(s);
        accounts.push((linear, [source, motor_work, radial_work, constraint]));
    }
    for i in 0..3 {
        assert!(
            state_error(
                finals[1][i],
                transformed(finals[0][i], p[i], rotation, boost, 0.18)
            ) < 3e-10
        );
    }
    near(
        accounts[1].1[0] - accounts[0].1[0],
        boost.dot(accounts[1].0),
        3e-11,
    );
    for i in 1..4 {
        near(accounts[0].1[i], accounts[1].1[i], 3e-11);
    }
}

#[test]
fn g4_motor_counts_six_rows_and_admits_exact_payload_and_complete_couple_bounds() {
    let p = [spherical(1.0, 1.0); 8];
    let s = core::array::from_fn::<_, 8, _>(|i| {
        state(
            p[i],
            X.scale(i.saturating_sub(5) as f64),
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
            Vec3::ZERO,
        )
    });
    let mut joints = (0..5)
        .map(|i| motor([i, i + 1], 0.0, 0.05 + 0.02 * i as f64))
        .collect::<Vec<_>>();
    for i in [5, 6] {
        joints.push(DrivenJoint::Distance {
            link: DistanceLink {
                bodies: [i, i + 1],
                points_body_m: [Vec3::ZERO; 2],
                length_m: 1.0,
            },
            drive: DistanceDrive {
                range_rate_mps: 0.0,
                range_acceleration_mps2: 0.02 - 0.03 * (i - 5) as f64,
            },
        });
    }
    let run = |joints: &[DrivenJoint], cap| {
        advance_driven_joint_graph(
            &s,
            &p,
            joints,
            free(),
            &[Wrench::ZERO; 8],
            0.02,
            cap,
            || false,
        )
    };
    let out = run(&joints, controls()).unwrap();
    assert_eq!(out.step.residuals.len(), 32);
    assert_eq!(out.joints.len(), 7);
    let mut cap = controls();
    cap.maximum_live_bytes = out.planned_live_bytes;
    assert_eq!(run(&joints, cap).unwrap(), out);
    cap.maximum_live_bytes -= 1;
    assert!(matches!(
        run(&joints, cap),
        Err(DistanceGraphError::Limit {
            resource: "live bytes",
            ..
        })
    ));
    joints.push(joints[5]);
    assert!(matches!(
        run(&joints, controls()),
        Err(DistanceGraphError::Limit {
            resource: "constraints",
            required: 33,
            ..
        })
    ));

    let pair = [s[0]; 2];
    let joint = [motor([0, 1], 0.0, 0.0)];
    let loads = [
        Wrench {
            force_world: Vec3::ZERO,
            torque_body: Vec3::new(1.0, 1.0, 1.0),
        },
        Wrench::ZERO,
    ];
    let evaluate = |cap| {
        advance_driven_joint_graph(&pair, &p[..2], &joint, free(), &loads, 0.01, cap, || false)
    };
    let baseline = evaluate(controls()).unwrap();
    assert_eq!(baseline.actuator_work_j, 0.0);
    for events in [
        &baseline.step.position_impulses,
        &baseline.step.velocity_impulses,
    ] {
        let mut complete = Vec3::ZERO;
        for event in &events[3..6] {
            let JointReactionImpulse::Couple(event) = event else {
                unreachable!()
            };
            assert!(norm(event.couple_on_second_world_nms) < 0.0035);
            complete = complete.add(event.couple_on_second_world_nms);
        }
        assert!(norm(complete) > 0.0035);
    }
    cap = controls();
    cap.step.maximum_impulse_ns = 0.0035 / 0.8;
    assert!(matches!(
        evaluate(cap),
        Err(DistanceGraphError::Step(
            DistanceStepError::LocalBranch | DistanceStepError::NotConverged
        ))
    ));
}

#[test]
fn g4_command_admits_whole_excursion_and_resolved_hidden_rate_before_state_changes() {
    let p = [spherical(2.0, 1.0); 2];
    let pair = [state(
        p[0],
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Vec3::ZERO,
    ); 2];
    let evaluate = |s: &[RigidBodyState], joints: &[DrivenJoint]| {
        advance_driven_joint_graph(
            s,
            &p,
            joints,
            free(),
            &[Wrench::ZERO; 2],
            0.01,
            controls(),
            || false,
        )
    };
    for (rate, acceleration) in [
        (core::f64::consts::TAU / 0.01, 0.0), // a hidden complete turn
        (400.0, -80_000.0), // zero endpoint displacement, one-radian interior peak
        (f64::MAX, f64::MAX),
        (f64::NAN, 0.0),
    ] {
        assert!(matches!(
            evaluate(&pair, &[motor([0, 1], rate, acceleration)]),
            Err(DistanceGraphError::PrescribedMotion { constraint: 3, .. })
        ));
    }
    let mut turn = pair;
    turn[0] = state(
        p[0],
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Z.scale(-0.15),
    );
    turn[1] = state(
        p[1],
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Z.scale(0.15),
    );
    let out = evaluate(&turn, &[motor([0, 1], 0.3, -60.0)]).unwrap();
    let report = out.joints[0].angular.unwrap();
    near(report.maximum_excursion_rad, 0.00075, 2e-18);
    near(report.displacement_rad, 0.0, 2e-18);
    near(out.actuator_work_j, 0.0, 2e-15);
    vector(
        out.step.states_after[0].angular_momentum_body(),
        Z.scale(0.15),
        2e-13,
    );
    vector(
        out.step.states_after[1].angular_momentum_body(),
        Z.scale(-0.15),
        2e-13,
    );
    let boost = X.scale(2_f64.powi(54));
    let mut hidden = core::array::from_fn::<_, 2, _>(|i| {
        state(
            p[i],
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            boost,
            Z.scale(if i == 0 { 0.0 } else { 0.301 }),
        )
    });
    assert!(matches!(
        evaluate(&hidden, &[motor([0, 1], 0.3, 0.0)]),
        Err(DistanceGraphError::InitialConstraint { constraint: 5 })
    ));
    hidden[1] = state(
        p[1],
        Vec3::ZERO,
        UnitQuaternion::from_axis_angle(Z, core::f64::consts::PI).unwrap(),
        boost,
        Z.scale(0.3),
    );
    assert!(matches!(
        evaluate(&hidden, &[motor([0, 1], 0.3, 0.0)]),
        Err(DistanceGraphError::Step(DistanceStepError::LocalBranch))
    ));
    assert_eq!(
        pair,
        [state(
            p[0],
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
            Vec3::ZERO
        ); 2]
    );
}

#[test]
fn g4_every_spatial_motor_cancellation_is_atomic_and_retry_is_identical() {
    let (p, s, joints) = spatial_mechanism();
    let loads = spatial_loads(UnitQuaternion::IDENTITY);
    let mut calls = 0;
    let reference = advance_spatial_driven_joint_graph(
        &s,
        &p,
        &joints,
        free(),
        &loads,
        0.02,
        controls(),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    assert!(calls > 20);
    for stop in 1..=calls {
        let mut current = 0;
        let result = advance_spatial_driven_joint_graph(
            &s,
            &p,
            &joints,
            free(),
            &loads,
            0.02,
            controls(),
            || {
                current += 1;
                current == stop
            },
        );
        assert!(
            matches!(
                result,
                Err(DistanceGraphError::Step(DistanceStepError::Cancelled))
            ),
            "checkpoint {stop}/{calls}: {result:?}"
        );
        assert_eq!((p, s, joints), spatial_mechanism());
        assert_eq!(loads, spatial_loads(UnitQuaternion::IDENTITY));
    }
    assert_eq!(
        advance_spatial_driven_joint_graph(
            &s,
            &p,
            &joints,
            free(),
            &loads,
            0.02,
            controls(),
            || false
        )
        .unwrap(),
        reference
    );
}
