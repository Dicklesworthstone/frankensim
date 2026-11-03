//! Finite sliding mechanisms, with independent continuous recoil references.

use fs_mbd::rigid_rattle::distance::graph::{
    ArticulatedJoint, DistanceDrive, DistanceGraphControls, DistanceGraphError, DistanceLink,
    DrivenJoint, DrivenJointGraphStep, JointReactionImpulse, PointPlaneJoint, PrismaticJoint,
    RevoluteJoint, advance_articulated_graph, advance_driven_joint_graph,
    advance_spatial_articulated_graph, advance_spatial_driven_joint_graph,
};
use fs_mbd::rigid_rattle::distance::{DistanceStepControls, DistanceStepError};
use fs_mbd::spatial_loads::{PreparedSpatialLoads, SpatialStepReceipt, WorldPointForce};
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
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    near(a.x, b.x, tolerance);
    near(a.y, b.y, tolerance);
    near(a.z, b.z, tolerance);
}
fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
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
fn free() -> RigidBodyIntegrator {
    RigidBodyIntegrator::new(Gravity::ZERO)
}
fn spherical(mass: f64, moment: f64) -> MassProperties {
    MassProperties::new(mass, Vec3::ZERO, Vec3::new(moment, moment, moment)).unwrap()
}
fn state(
    p: MassProperties,
    x: Vec3,
    q: UnitQuaternion,
    v: Vec3,
    world_omega: Vec3,
) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(x, q).unwrap(),
        v.scale(p.mass()),
        p.angular_momentum_body_checked(q.rotate_world_to_body(world_omega))
            .unwrap(),
    )
    .unwrap()
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
fn inverse(q: UnitQuaternion) -> UnitQuaternion {
    let [w, x, y, z] = q.components();
    UnitQuaternion::new(w, -x, -y, -z).unwrap()
}
fn reverse(s: RigidBodyState) -> RigidBodyState {
    RigidBodyState::new(
        s.pose(),
        s.linear_momentum_world().scale(-1.0),
        s.angular_momentum_body().scale(-1.0),
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
fn rotor(angle: f64) -> UnitQuaternion {
    UnitQuaternion::from_axis_angle(Z, angle).unwrap()
}
fn x_rail() -> UnitQuaternion {
    UnitQuaternion::from_axis_angle(Y, core::f64::consts::FRAC_PI_2).unwrap()
}
fn radial(angle: f64) -> Vec3 {
    Vec3::new(angle.cos(), angle.sin(), 0.0)
}

#[test]
fn g1_finite_slider_shares_transverse_load_and_axial_torque_but_keeps_free_travel() {
    let p = [spherical(2.0, 1.3), spherical(5.0, 2.1)];
    let fractions = [5.0 / 7.0, -2.0 / 7.0];
    let boost = Vec3::new(0.2, -0.1, 0.3);
    let joint =
        PrismaticJoint::new([0, 1], [Vec3::ZERO; 2], [UnitQuaternion::IDENTITY; 2], 0.8).unwrap();
    assert_eq!(joint.bodies(), [0, 1]);
    assert_eq!(joint.points_body_m(), [Vec3::ZERO; 2]);
    assert_eq!(joint.frames_body(), [UnitQuaternion::IDENTITY; 2]);
    assert_eq!(joint.angular_length_m(), 0.8);
    for (force, slide_rate) in [(1.1, 0.0), (0.0, 0.4)] {
        let initial = core::array::from_fn::<_, 2, _>(|i| {
            state(
                p[i],
                Z.scale(fractions[i] * 1.2),
                rotor(0.2),
                boost.add(Z.scale(fractions[i] * slide_rate)),
                Z.scale(0.3),
            )
        });
        let loads = [
            PreparedSpatialLoads::new(
                Wrench {
                    force_world: X.scale(force),
                    torque_body: Vec3::ZERO,
                },
                Z.scale(0.7),
                &[],
                || false,
            )
            .unwrap(),
            // This known source couple cancels the transverse reaction moment
            // on the carrier, giving an exact translating, axially spinning pair.
            PreparedSpatialLoads::new(
                Wrench::ZERO,
                Y.scale(-force * fractions[0] * 1.2),
                &[],
                || false,
            )
            .unwrap(),
        ];
        let mut s = initial.to_vec();
        let before = totals(&s, &p);
        let (mut linear, mut angular, mut work, mut defect) = (Vec3::ZERO, Vec3::ZERO, 0.0, 0.0);
        for step_index in 0..25 {
            let out = advance_spatial_articulated_graph(
                &s,
                &p,
                &[ArticulatedJoint::Prismatic(joint)],
                free(),
                &loads,
                0.02,
                controls(),
                || false,
            )
            .unwrap();
            assert_eq!(out.joint_rows[0].row_count, 5);
            angular = angular.add(out.constraint_angular_impulse_world_nms);
            work += out.constraint_work_j;
            for (i, step) in out.step.free_steps.iter().enumerate() {
                assert_eq!(
                    *step,
                    free()
                        .step_spatial(step.step.state_before, p[i], &loads[i], 0.02)
                        .unwrap()
                );
                linear = linear.add(step.linear_impulse_world);
                angular = angular.add(step.angular_impulse_world);
                work += step.force_work_j + step.couple_work_j;
                defect += step.work_energy_residual_j;
            }
            s = out.step.states_after;
            let t = f64::from(step_index + 1) * 0.02;
            for i in 0..2 {
                let expected = state(
                    p[i],
                    boost
                        .scale(t)
                        .add(Z.scale(fractions[i] * (1.2 + slide_rate * t)))
                        .add(X.scale(0.5 * force * t * t / 7.0)),
                    rotor(0.2 + 0.3 * t + 0.5 * 0.7 * t * t / 3.4),
                    boost
                        .add(Z.scale(fractions[i] * slide_rate))
                        .add(X.scale(force * t / 7.0)),
                    Z.scale(0.3 + 0.7 * t / 3.4),
                );
                assert!(
                    state_error(s[i], expected) < 2e-10,
                    "body={i}, t={t}, error={}",
                    state_error(s[i], expected)
                );
            }
        }
        let after = totals(&s, &p);
        vector(after.0.sub(before.0), linear, 3e-12);
        vector(after.1.sub(before.1), angular, 2e-11);
        near(after.2 - before.2, work + defect, 8e-12);
    }
}

#[test]
fn g1_rotating_finite_rail_matches_continuous_centrifugal_motion_and_reverses() {
    let p = [spherical(2.0, 0.7), spherical(3.0, 1.1)];
    let mu = 6.0 / 5.0;
    let joint = ArticulatedJoint::Prismatic(
        PrismaticJoint::new([0, 1], [Vec3::ZERO; 2], [x_rail(); 2], 0.9).unwrap(),
    );
    let from_reference = |y: [f64; 4]| {
        let u = radial(y[2]);
        let d = u.scale(y[1]).add(Z.cross(u).scale(y[0] * y[3]));
        [0.6, -0.4]
            .map(|fraction| (u.scale(fraction * y[0]), d.scale(fraction)))
            .into_iter()
            .enumerate()
            .map(|(i, (x, v))| state(p[i], x, rotor(y[2]), v, Z.scale(y[3])))
            .collect::<Vec<_>>()
    };
    let original = [1.2, 0.15, 0.2, 0.6];
    let initial = from_reference(original);
    let time = 0.8;
    let rhs = |y: [f64; 4]| {
        [
            y[1],
            y[0] * y[3] * y[3],
            y[3],
            -2.0 * mu * y[0] * y[1] * y[3] / (1.8 + mu * y[0] * y[0]),
        ]
    };
    let mut reference = original;
    let h = time / 8192.0;
    for _ in 0..8192 {
        let a = rhs(reference);
        let b = rhs(core::array::from_fn(|i| reference[i] + 0.5 * h * a[i]));
        let c = rhs(core::array::from_fn(|i| reference[i] + 0.5 * h * b[i]));
        let d = rhs(core::array::from_fn(|i| reference[i] + h * c[i]));
        reference = core::array::from_fn(|i| {
            reference[i] + h * (a[i] + 2.0 * b[i] + 2.0 * c[i] + d[i]) / 6.0
        });
    }
    let expected = from_reference(reference);
    let before = totals(&initial, &p);
    let mut errors = Vec::new();
    for count in [16, 32, 64] {
        let h = time / f64::from(count);
        let mut s = initial.clone();
        let mut work = 0.0;
        for _ in 0..count {
            let out = advance_articulated_graph(
                &s,
                &p,
                &[joint],
                free(),
                &[Wrench::ZERO; 2],
                h,
                controls(),
                || false,
            )
            .unwrap();
            work += out.constraint_work_j;
            s = out.step.states_after;
        }
        errors.push(
            s.iter()
                .zip(&expected)
                .map(|(&a, &b)| state_error(a, b))
                .fold(0.0_f64, f64::max),
        );
        let after = totals(&s, &p);
        vector(after.0, before.0, 8e-12);
        vector(after.1, before.1, 2e-11);
        near(after.2 - before.2, work, 8e-12);
        s = s.into_iter().map(reverse).collect();
        for _ in 0..count {
            s = advance_articulated_graph(
                &s,
                &p,
                &[joint],
                free(),
                &[Wrench::ZERO; 2],
                h,
                controls(),
                || false,
            )
            .unwrap()
            .step
            .states_after;
        }
        for i in 0..2 {
            assert!(state_error(reverse(s[i]), initial[i]) < 2e-9);
        }
    }
    for pair in errors.windows(2) {
        assert!((3.7..4.3).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
}

// Independent free-floating slider-crank reduction. Carrier, finite crank and
// finite slider are joined by a massless connecting rod; all three bodies move.
// Eliminating the center of mass gives T=(A theta'^2+2 B theta' phi'+C phi'^2)/2.
struct Crank {
    mass: [f64; 3],
    moment: [f64; 3],
    crank: f64,
    rod: f64,
    anchor: f64,
    length: f64,
    rate: f64,
    acceleration: f64,
    angular_momentum: f64,
}
impl Crank {
    fn shape(&self, phi: f64) -> ([Vec3; 3], [Vec3; 3], f64, f64, f64) {
        let (s, c) = phi.sin_cos();
        let root = (self.rod * self.rod - self.crank * self.crank * s * s).sqrt();
        let x = self.crank * c + root;
        let dx = -self.crank * s - self.crank * self.crank * s * c / root;
        let raw = [Vec3::ZERO, radial(phi).scale(0.5 * self.crank), X.scale(x)];
        let tangent = [
            Vec3::ZERO,
            Z.cross(radial(phi)).scale(0.5 * self.crank),
            X.scale(dx),
        ];
        let mass = self.mass.iter().sum::<f64>();
        let center = (0..3).fold(Vec3::ZERO, |sum, i| {
            sum.add(raw[i].scale(self.mass[i] / mass))
        });
        let center_rate = (0..3).fold(Vec3::ZERO, |sum, i| {
            sum.add(tangent[i].scale(self.mass[i] / mass))
        });
        let r = raw.map(|p| p.sub(center));
        let d = tangent.map(|p| p.sub(center_rate));
        let a = self.moment.iter().sum::<f64>()
            + (0..3).map(|i| self.mass[i] * r[i].dot(r[i])).sum::<f64>();
        let b = self.moment[1]
            + (0..3)
                .map(|i| self.mass[i] * Z.cross(r[i]).dot(d[i]))
                .sum::<f64>();
        let c = self.moment[1] + (0..3).map(|i| self.mass[i] * d[i].dot(d[i])).sum::<f64>();
        (r, d, a, b, c)
    }
    fn command(&self, t: f64) -> (f64, f64) {
        let length = self.length + self.rate * t + 0.5 * self.acceleration * t * t;
        let rate = self.rate + self.acceleration * t;
        let phi = ((self.crank * self.crank + self.anchor * self.anchor - length * length)
            / (2.0 * self.crank * self.anchor))
            .asin();
        let phi_rate = -length * rate / (self.crank * self.anchor * phi.cos());
        (phi, phi_rate)
    }
    fn carrier_rate(&self, t: f64) -> f64 {
        let (phi, rate) = self.command(t);
        let (_, _, a, b, _) = self.shape(phi);
        (self.angular_momentum - b * rate) / a
    }
    fn state(&self, t: f64, theta: f64) -> [RigidBodyState; 3] {
        let (phi, rate) = self.command(t);
        let carrier_rate = self.carrier_rate(t);
        let (r, d, _, _, _) = self.shape(phi);
        let q = rotor(theta);
        core::array::from_fn(|i| {
            state(
                spherical(self.mass[i], self.moment[i]),
                q.rotate_body_to_world(r[i]),
                rotor(theta + if i == 1 { phi } else { 0.0 }),
                q.rotate_body_to_world(Z.cross(r[i]).scale(carrier_rate).add(d[i].scale(rate))),
                Z.scale(carrier_rate + if i == 1 { rate } else { 0.0 }),
            )
        })
    }
}

#[test]
fn g1_actuated_closed_slider_crank_matches_independent_finite_recoil_quadrature() {
    let mut model = Crank {
        mass: [4.0, 1.5, 2.0],
        moment: [1.1, 0.3, 0.5],
        crank: 0.7,
        rod: 1.5,
        anchor: 0.5,
        length: (0.74 - 0.7 * 0.3_f64.sin()).sqrt(),
        rate: -0.08,
        acceleration: 0.04,
        angular_momentum: 0.0,
    };
    let (phi, rate) = model.command(0.0);
    let (_, _, a, b, _) = model.shape(phi);
    model.angular_momentum = a * 0.25 + b * rate;
    let initial = model.state(0.0, 0.2);
    let p = core::array::from_fn::<_, 3, _>(|i| spherical(model.mass[i], model.moment[i]));
    let original = [
        DrivenJoint::Revolute(
            RevoluteJoint::new(
                [0, 1],
                [Vec3::ZERO, X.scale(-0.5 * model.crank)],
                [Z; 2],
                0.8,
            )
            .unwrap(),
        ),
        DrivenJoint::Prismatic(
            PrismaticJoint::new([2, 0], [Vec3::ZERO; 2], [x_rail(); 2], 0.9).unwrap(),
        ),
        DrivenJoint::Distance {
            link: DistanceLink {
                bodies: [1, 2],
                points_body_m: [X.scale(0.5 * model.crank), Vec3::ZERO],
                length_m: model.rod,
            },
            drive: DistanceDrive::default(),
        },
        DrivenJoint::Distance {
            link: DistanceLink {
                bodies: [0, 1],
                points_body_m: [Y.scale(model.anchor), X.scale(0.5 * model.crank)],
                length_m: model.length,
            },
            drive: DistanceDrive {
                range_rate_mps: model.rate,
                range_acceleration_mps2: model.acceleration,
            },
        },
    ];
    let time = 0.6;
    let h = time / 8192.0;
    let mut theta = 0.2;
    for i in 0..8192 {
        let t = f64::from(i) * h;
        theta += h
            * (model.carrier_rate(t)
                + 4.0 * model.carrier_rate(t + 0.5 * h)
                + model.carrier_rate(t + h))
            / 6.0;
    }
    let expected = model.state(time, theta);
    let before = totals(&initial, &p);
    let reference_work = totals(&expected, &p).2 - before.2;
    let mut errors = Vec::new();
    let mut work_errors = Vec::new();
    for count in [15, 30, 60] {
        let h = time / f64::from(count);
        let mut s = initial.to_vec();
        let mut joints = original.to_vec();
        let (mut work, mut defect) = (0.0, 0.0);
        for _ in 0..count {
            let out = advance_driven_joint_graph(
                &s,
                &p,
                &joints,
                free(),
                &[Wrench::ZERO; 3],
                h,
                controls(),
                || false,
            )
            .unwrap();
            assert_eq!(out.step.residuals.len(), 12);
            for passive in &out.joints[..3] {
                assert_eq!(passive.actuator_work_j, 0.0);
            }
            work += out.actuator_work_j;
            defect += out.constraint_work_defect_j;
            near(
                out.constraint_work_j,
                out.actuator_work_j + out.constraint_work_defect_j,
                3e-13,
            );
            joints = out.joints.iter().map(|j| j.joint_after).collect();
            s = out.step.states_after;
        }
        errors.push(
            s.iter()
                .zip(expected)
                .map(|(&a, b)| state_error(a, b))
                .fold(0.0_f64, f64::max),
        );
        work_errors.push((work - reference_work).abs());
        let after = totals(&s, &p);
        vector(after.0, before.0, 8e-12);
        vector(after.1, before.1, 3e-11);
        near(after.2 - before.2, work + defect, 1e-11);
    }
    for errors in [errors, work_errors] {
        for pair in errors.windows(2) {
            assert!((3.6..4.4).contains(&(pair[0] / pair[1])), "{errors:?}");
        }
    }
}

fn spatial_mechanism() -> ([MassProperties; 3], [RigidBodyState; 3], [DrivenJoint; 3]) {
    let p = core::array::from_fn(|i| {
        MassProperties::from_principal_frame(
            1.5 + 0.7 * i as f64,
            Vec3::ZERO,
            Vec3::new(0.7, 1.0, 1.3).scale(1.0 + 0.2 * i as f64),
            UnitQuaternion::from_axis_angle(Vec3::new(0.4, -1.0, 0.7), 0.23 + 0.14 * i as f64)
                .unwrap(),
        )
        .unwrap()
    });
    let q = [
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, -0.2, 0.5), 0.4).unwrap(),
        UnitQuaternion::from_axis_angle(Vec3::new(-0.3, 0.8, 0.4), -0.3).unwrap(),
        UnitQuaternion::from_axis_angle(Vec3::new(0.2, 0.3, -1.0), 0.5).unwrap(),
    ];
    let frame = UnitQuaternion::from_axis_angle(Vec3::new(0.7, 0.3, -0.2), 0.6).unwrap();
    let points = [Vec3::new(-0.15, 0.12, -0.1), Vec3::new(0.1, -0.2, 0.05)];
    let hinge_points = [Vec3::new(-0.3, 0.25, 0.15), Vec3::new(0.2, -0.1, 0.25)];
    let carrier = Vec3::new(0.2, -0.1, 0.3);
    let anchor = carrier.add(q[1].rotate_body_to_world(points[1]));
    let hinge = carrier.add(q[1].rotate_body_to_world(hinge_points[0]));
    let x = [
        anchor
            .add(frame.rotate_body_to_world(Z).scale(0.8))
            .sub(q[0].rotate_body_to_world(points[0])),
        carrier,
        hinge.sub(q[2].rotate_body_to_world(hinge_points[1])),
    ];
    let omega = Vec3::new(0.3, -0.2, 0.4);
    let velocity = Vec3::new(0.1, 0.2, -0.15);
    let s =
        core::array::from_fn(|i| state(p[i], x[i], q[i], velocity.add(omega.cross(x[i])), omega));
    let drive_points = [Vec3::new(0.2, -0.1, 0.05), Vec3::new(-0.1, 0.15, -0.2)];
    let drive_separation = x[2]
        .add(q[2].rotate_body_to_world(drive_points[1]))
        .sub(x[0].add(q[0].rotate_body_to_world(drive_points[0])));
    let joints = [
        DrivenJoint::Prismatic(
            PrismaticJoint::new(
                [0, 1],
                points,
                [compose(inverse(q[0]), frame), compose(inverse(q[1]), frame)],
                0.8,
            )
            .unwrap(),
        ),
        DrivenJoint::Revolute(
            RevoluteJoint::new(
                [1, 2],
                hinge_points,
                [1, 2].map(|i| q[i].rotate_world_to_body(frame.rotate_body_to_world(X))),
                0.7,
            )
            .unwrap(),
        ),
        DrivenJoint::Distance {
            link: DistanceLink {
                bodies: [0, 2],
                points_body_m: drive_points,
                length_m: norm(drive_separation),
            },
            drive: DistanceDrive {
                range_rate_mps: 0.0,
                range_acceleration_mps2: 0.12,
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
fn g2_full_tensor_spatial_rail_hinge_and_actuator_close_balances_and_transform_together() {
    let (p, initial, original) = spatial_mechanism();
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(-0.4, 0.2, 0.7), 1.1).unwrap();
    let boost = Vec3::new(0.7, -0.4, 0.5);
    let configurations = [(UnitQuaternion::IDENTITY, Vec3::ZERO), (rotation, boost)];
    let mut final_states = Vec::new();
    let mut impulse_accounts = Vec::new();
    let mut work_accounts = Vec::new();
    for (rotation, boost) in configurations {
        let loads = spatial_loads(rotation);
        let mut s = core::array::from_fn::<_, 3, _>(|i| {
            transformed(initial[i], p[i], rotation, boost, 0.0)
        })
        .to_vec();
        let mut joints = original.to_vec();
        let before = totals(&s, &p);
        let (mut linear, mut angular) = (Vec3::ZERO, Vec3::ZERO);
        let (mut source, mut actuator, mut constraint, mut defect) = (0.0, 0.0, 0.0, 0.0);
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
                [5, 5, 1]
            );
            for event in [&out.step.position_impulses, &out.step.velocity_impulses] {
                for point in &event[..2] {
                    let JointReactionImpulse::Point(point) = point else {
                        panic!("transverse native point impulse required")
                    };
                    vector(
                        point.first.point_kinematics_before.point_world,
                        point.second.point_kinematics_before.point_world,
                        1e-12,
                    );
                    let DrivenJoint::Prismatic(joint) = joints[0] else {
                        unreachable!()
                    };
                    // The reaction acts where the carriage is NOW, a substantial
                    // axial distance from the authored carrier rail anchor.
                    assert!(norm(point.second.arm_body.sub(joint.points_body_m()[1])) > 0.7);
                }
                for couple in &event[2..5] {
                    let JointReactionImpulse::Couple(couple) = couple else {
                        panic!("three native frame-lock couples required")
                    };
                    vector(couple.first.linear_increment_world, Vec3::ZERO, 0.0);
                    vector(couple.second.linear_increment_world, Vec3::ZERO, 0.0);
                }
            }
            angular = angular.add(out.constraint_angular_impulse_world_nms);
            actuator += out.actuator_work_j;
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
        vector(after.0.sub(before.0), linear, 8e-12);
        vector(after.1.sub(before.1), angular, 2e-11);
        near(
            after.2 - before.2,
            source + actuator + constraint + defect,
            2e-11,
        );
        assert!(actuator.abs() > 1e-5);
        final_states.push(s);
        impulse_accounts.push(linear);
        work_accounts.push([source, actuator, constraint]);
    }
    for i in 0..3 {
        assert!(
            state_error(
                final_states[1][i],
                transformed(final_states[0][i], p[i], rotation, boost, 0.18)
            ) < 2e-10
        );
    }
    near(
        work_accounts[1][0] - work_accounts[0][0],
        boost.dot(impulse_accounts[1]),
        2e-11,
    );
    for i in 1..3 {
        near(work_accounts[0][i], work_accounts[1][i], 2e-11);
    }
}

#[test]
fn g4_prismatic_32_row_graph_admits_exact_payload_and_refuses_branches_and_vector_caps() {
    let p = [spherical(1.0, 1.0); 7];
    let s = core::array::from_fn::<_, 7, _>(|i| {
        state(
            p[i],
            Z.scale(i as f64),
            UnitQuaternion::IDENTITY,
            Vec3::ZERO,
            Z.scale(0.2),
        )
    });
    let mut joints = (0..6)
        .map(|i| {
            DrivenJoint::Prismatic(
                PrismaticJoint::new(
                    [i, i + 1],
                    [Vec3::ZERO; 2],
                    [UnitQuaternion::IDENTITY; 2],
                    1.0,
                )
                .unwrap(),
            )
        })
        .collect::<Vec<_>>();
    for (i, acceleration) in [(0, 0.1), (5, -0.06)] {
        joints.push(DrivenJoint::Distance {
            link: DistanceLink {
                bodies: [i, i + 1],
                points_body_m: [Vec3::ZERO; 2],
                length_m: 1.0,
            },
            drive: DistanceDrive {
                range_rate_mps: 0.0,
                range_acceleration_mps2: acceleration,
            },
        });
    }
    let evaluate = |joints: &[DrivenJoint], cap| {
        advance_driven_joint_graph(&s, &p, joints, free(), &[Wrench::ZERO; 7], 0.1, cap, || {
            false
        })
    };
    let out = evaluate(&joints, controls()).unwrap();
    assert_eq!(out.step.residuals.len(), 32);
    assert_eq!(out.joints.len(), 8);
    assert!(out.actuator_work_j > 0.0);
    let mut cap = controls();
    cap.maximum_live_bytes = out.planned_live_bytes;
    assert_eq!(evaluate(&joints, cap).unwrap(), out);
    cap.maximum_live_bytes -= 1;
    assert!(matches!(
        evaluate(&joints, cap),
        Err(DistanceGraphError::Limit {
            resource: "live bytes",
            ..
        })
    ));
    assert!(matches!(
        evaluate(&[joints[0], joints[0]], controls()),
        Err(DistanceGraphError::DependentConstraints { .. })
    ));
    joints.push(DrivenJoint::PointPlane(
        PointPlaneJoint::new([0, 6], Vec3::ZERO, Vec3::ZERO, X).unwrap(),
    ));
    assert!(matches!(
        evaluate(&joints, controls()),
        Err(DistanceGraphError::Limit {
            resource: "constraints",
            required: 33,
            ..
        })
    ));

    let make_joint = |length| {
        PrismaticJoint::new(
            [0, 1],
            [Vec3::ZERO; 2],
            [UnitQuaternion::IDENTITY; 2],
            length,
        )
    };
    for length in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(make_joint(length).is_err());
    }
    assert!(
        PrismaticJoint::new([0, 0], [Vec3::ZERO; 2], [UnitQuaternion::IDENTITY; 2], 1.0).is_err()
    );
    assert!(
        PrismaticJoint::new(
            [0, 1],
            [Vec3::ZERO, Vec3::new(f64::NAN, 0.0, 0.0)],
            [UnitQuaternion::IDENTITY; 2],
            1.0
        )
        .is_err()
    );
    let pair = [state(
        p[0],
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        Vec3::ZERO,
        Vec3::ZERO,
    ); 2];
    let run = |states: &[RigidBodyState], length, cap, loads: &[Wrench]| {
        advance_articulated_graph(
            states,
            &p[..2],
            &[ArticulatedJoint::Prismatic(make_joint(length).unwrap())],
            free(),
            loads,
            0.01,
            cap,
            || false,
        )
    };
    assert!(matches!(
        run(
            &pair,
            controls().position_tolerance_m,
            controls(),
            &[Wrench::ZERO; 2]
        ),
        Err(DistanceGraphError::Step(DistanceStepError::Invalid(_)))
    ));
    let mut bad = pair;
    bad[1] = state(
        p[1],
        Vec3::ZERO,
        rotor(core::f64::consts::PI),
        Vec3::ZERO,
        Vec3::ZERO,
    );
    assert!(matches!(
        run(&bad, 1.0, controls(), &[Wrench::ZERO; 2]),
        Err(DistanceGraphError::Step(DistanceStepError::LocalBranch))
    ));
    for i in 0..2 {
        bad[i] = state(
            p[i],
            Vec3::ZERO,
            UnitQuaternion::IDENTITY,
            X.scale(2_f64.powi(54)),
            if i == 0 { Vec3::ZERO } else { Z.scale(0.01) },
        );
    }
    assert!(matches!(
        run(&bad, 1.0, controls(), &[Wrench::ZERO; 2]),
        Err(DistanceGraphError::InitialConstraint { constraint: 4 })
    ));
    // Each Cartesian component fits its scalar tolerance. The complete
    // three-row frame error/rate must still satisfy the declared vector norm.
    for (orientation, omega) in [
        (
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, 1.0, 1.0), 3_f64.sqrt() * 4e-13)
                .unwrap(),
            Vec3::ZERO,
        ),
        (UnitQuaternion::IDENTITY, Vec3::new(4e-13, 4e-13, 4e-13)),
    ] {
        bad = pair;
        bad[1] = state(p[1], Vec3::ZERO, orientation, Vec3::ZERO, omega);
        assert!(matches!(
            run(&bad, 1.0, controls(), &[Wrench::ZERO; 2]),
            Err(DistanceGraphError::InitialConstraint { constraint: 2 })
        ));
    }

    let loads = [
        Wrench {
            force_world: Vec3::ZERO,
            torque_body: Vec3::new(1.0, 1.0, 1.0),
        },
        Wrench::ZERO,
    ];
    let baseline = run(&pair, 1.0, controls(), &loads).unwrap();
    for events in [
        &baseline.step.position_impulses,
        &baseline.step.velocity_impulses,
    ] {
        let mut total = Vec3::ZERO;
        for event in &events[2..5] {
            let JointReactionImpulse::Couple(event) = event else {
                unreachable!()
            };
            assert!(norm(event.couple_on_second_world_nms) < 0.0035);
            total = total.add(event.couple_on_second_world_nms);
        }
        assert!(norm(total) > 0.0035);
    }
    cap = controls();
    cap.step.maximum_impulse_ns = 0.0035;
    assert!(matches!(
        run(&pair, 1.0, cap, &loads),
        Err(DistanceGraphError::Step(
            DistanceStepError::LocalBranch | DistanceStepError::NotConverged
        ))
    ));
}

#[test]
fn g4_every_spatial_prismatic_cancellation_is_atomic_and_retry_is_identical() {
    let (p, s, joints) = spatial_mechanism();
    let loads = spatial_loads(UnitQuaternion::IDENTITY);
    let mut calls = 0;
    let reference: DrivenJointGraphStep<SpatialStepReceipt> = advance_spatial_driven_joint_graph(
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
            || false,
        )
        .unwrap(),
        reference
    );
}
