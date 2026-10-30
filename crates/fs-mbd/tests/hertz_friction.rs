//! Two finite original bodies: coupled Hertz force/spin and actual endpoint work.

use core::f64::consts::PI;
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::RigidContactBinding;
use fs_mbd::tangential_friction::hertz::*;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_tribo::partial_slip::TangentFrame;
use fs_tribo::sliding_patch::contact::*;

const RADIUS: f64 = 0.03;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; tolerance {tolerance:e}"
    );
}
fn close_vec(a: Vec3, b: Vec3, tolerance: f64) {
    close(a.x, b.x, tolerance);
    close(a.y, b.y, tolerance);
    close(a.z, b.z, tolerance);
}
fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn norm(v: Vec3) -> f64 {
    v.x.hypot(v.y).hypot(v.z)
}
fn mul(c: [[f64; 3]; 3], x: [f64; 3]) -> [f64; 3] {
    c.map(|row| row.into_iter().zip(x).map(|(a, b)| a * b).sum())
}
fn controls() -> RigidHertzControls {
    RigidHertzControls {
        maximum_bodies: 8,
        geometry_tolerance_m: 2e-12,
        velocity_tolerance_m_s: 2e-11,
        angular_velocity_tolerance_rad_s: 2e-10,
        linear_momentum_tolerance_ns: 2e-11,
        angular_momentum_tolerance_nms: 2e-11,
        energy_tolerance_j: 2e-10,
        contact: HertzWrenchConfig {
            maximum_static_iterations: 96,
            maximum_iterations: 256,
            maximum_evaluations: 512,
            force_tolerance_n: 2e-9,
            moment_tolerance_nm: 2e-10,
            displacement_tolerance_m: 2e-13,
            energy_tolerance_j: 2e-11,
            relative_matrix_tolerance: 1e-13,
            maximum_force_n: 20.0,
            maximum_spin_moment_nm: 2.0,
        },
    }
}
fn law() -> HertzWrenchContact {
    HertzWrenchContact::new(HertzWrenchParameters {
        normal_force_n: 10.0,
        contact_radius_m: RADIUS,
        static_mu: 0.6,
        kinetic_mu: 0.4,
    })
    .unwrap()
}
fn state(
    pose: Pose,
    properties: MassProperties,
    velocity: Vec3,
    omega_world: Vec3,
) -> RigidBodyState {
    let omega = pose.orientation().rotate_world_to_body(omega_world);
    let inertia = properties.principal_inertia_body();
    RigidBodyState::new(
        pose,
        velocity.scale(properties.mass()),
        Vec3::new(
            inertia.x * omega.x,
            inertia.y * omega.y,
            inertia.z * omega.z,
        ),
    )
    .unwrap()
}

struct Fixture {
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    binding: RigidContactBinding,
    frame: TangentFrame,
    law: HertzWrenchContact,
}
impl Fixture {
    fn contact(&self) -> RigidHertzContact<'_> {
        RigidHertzContact {
            binding: self.binding,
            frame: self.frame,
            law: &self.law,
        }
    }
    fn advance(
        &self,
        bodies: &mut [RigidBodyState],
        dt: f64,
        cfg: RigidHertzControls,
    ) -> Result<RigidHertzReport, RigidHertzError> {
        apply_rigid_hertz(
            bodies,
            &self.properties,
            self.contact(),
            dt,
            cfg,
            &CancelGate::new(),
        )
    }
}

fn axial(slip: f64, spin: f64) -> Fixture {
    let properties = [
        MassProperties::new(2.0, Vec3::ZERO, Vec3::new(0.08, 0.10, 0.12)).unwrap(),
        MassProperties::new(3.0, Vec3::ZERO, Vec3::new(0.15, 0.20, 0.25)).unwrap(),
    ];
    let poses = [
        Pose::new(Vec3::new(0.0, 0.0, 0.2), UnitQuaternion::IDENTITY).unwrap(),
        Pose::new(Vec3::new(0.0, 0.0, -0.15), UnitQuaternion::IDENTITY).unwrap(),
    ];
    Fixture {
        bodies: [
            state(
                poses[0],
                properties[0],
                Vec3::new(slip, 0.0, 0.0),
                Vec3::new(0.0, 0.0, spin),
            ),
            state(poses[1], properties[1], Vec3::ZERO, Vec3::ZERO),
        ],
        properties,
        binding: RigidContactBinding {
            body_a: 0,
            body_b: 1,
            arm_a_body: Vec3::new(0.0, 0.0, -0.2),
            arm_b_body: Vec3::new(0.0, 0.0, 0.15),
        },
        frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
        law: law(),
    }
}

#[test]
fn moving_high_static_capacity_uses_kinetic_traction_through_true_arrest() {
    for (slip, spin) in [(1.0, 0.0), (0.0, 1.0)] {
        let mut f = axial(slip, spin);
        f.law = HertzWrenchContact::new(HertzWrenchParameters {
            static_mu: 1.5,
            ..f.law.parameters()
        })
        .unwrap();
        let (capacity, mobility) = if slip != 0.0 {
            (
                4.0,
                1.0 / 2.0 + 1.0 / 3.0 + 0.2 * 0.2 / 0.10 + 0.15 * 0.15 / 0.20,
            )
        } else {
            (4.0 * 3.0 * PI * RADIUS / 16.0, 1.0 / 0.12 + 1.0 / 0.25)
        };
        let stop = 1.0 / (capacity * mobility);
        for fraction in [0.75, 1.0] {
            let mut bodies = f.bodies;
            let report = f.advance(&mut bodies, fraction * stop, controls()).unwrap();
            assert_eq!(report.contact.mode, HertzWrenchMode::Sliding);
            let after = if slip != 0.0 {
                report.slip_velocity_after_m_s[0]
            } else {
                report.spin_after_rad_s
            };
            close(after, 1.0 - fraction, 3e-14);
            close(
                report.contact.dissipated_energy_j,
                (fraction - 0.5 * fraction * fraction) / mobility,
                3e-14,
            );
            close(report.total_energy_residual_j, 0.0, 3e-14);
        }
        assert_eq!(f.law.parameters().static_mu, 1.5);
    }
}

#[test]
fn positive_static_but_zero_kinetic_capacity_does_not_capture_existing_motion() {
    let mut f = axial(1e-20, 0.0);
    f.law = HertzWrenchContact::new(HertzWrenchParameters {
        static_mu: 1.5,
        kinetic_mu: 0.0,
        ..f.law.parameters()
    })
    .unwrap();
    let mut bodies = f.bodies;
    let report = f.advance(&mut bodies, 1.0, controls()).unwrap();
    assert_eq!(bodies, f.bodies);
    assert_eq!(report.contact.mode, HertzWrenchMode::Frictionless);
    assert_eq!(report.contact.dissipated_energy_j, 0.0);
    close(report.slip_velocity_after_m_s[0], 1e-20, 1e-35);
}

#[test]
fn exact_kinetic_translation_stops_relative_point_slip_with_finite_support_recoil() {
    let f = axial(1.0, 0.0);
    // Independent scalar mechanics from both original masses and lever inertias.
    let inverse_mass = 1.0 / 2.0 + 1.0 / 3.0 + 0.2 * 0.2 / 0.10 + 0.15 * 0.15 / 0.20;
    let stop = 1.0 / (4.0 * inverse_mass);
    let mut bodies = f.bodies;
    let report = f.advance(&mut bodies, stop, controls()).unwrap();
    assert_eq!(report.contact.mode, HertzWrenchMode::Sliding);
    close(report.contact.force_n[0], -4.0, 2e-15);
    close(report.contact.spin_moment_nm, 0.0, 1e-16);
    close(report.scaled_mobility_per_kg[0][0], inverse_mass, 1e-15);
    close(report.slip_velocity_after_m_s[0], 0.0, 5e-16);
    close(report.dissipated_energy_j, 0.5 / inverse_mass, 2e-15);
    close(
        report.kinetic_after_j - report.kinetic_before_j,
        -0.5 / inverse_mass,
        2e-15,
    );
    close(bodies[0].linear_momentum_world().x, 2.0 - 4.0 * stop, 2e-15);
    close(bodies[1].linear_momentum_world().x, 4.0 * stop, 2e-15);
    close(bodies[0].angular_momentum_body().y, 0.2 * 4.0 * stop, 2e-15);
    close(
        bodies[1].angular_momentum_body().y,
        0.15 * 4.0 * stop,
        2e-15,
    );
    assert!(bodies[1].linear_momentum_world().x > 0.0);
    assert!(norm(report.linear_momentum_change_ns) < 1e-15);
    assert!(norm(report.angular_momentum_change_nms) < 1e-15);
    for i in 0..2 {
        assert_eq!(bodies[i].pose(), f.bodies[i].pose());
    }
    let held = f.advance(&mut bodies, 0.01, controls()).unwrap();
    assert_eq!(held.contact.mode, HertzWrenchMode::Sticking);
    assert_eq!(held.dissipated_energy_j, 0.0);
    close(held.slip_velocity_after_m_s[0], 0.0, 1e-14);
}

#[test]
fn pure_spin_stopping_preserves_total_spin_and_heats_only_relative_rotation() {
    let f = axial(0.0, 0.3);
    let inverse_inertia = 1.0 / 0.12 + 1.0 / 0.25;
    let capacity = (3.0 * PI / 16.0) * 4.0 * RADIUS;
    let stop = 0.3 / (inverse_inertia * capacity);
    let mut bodies = f.bodies;
    let report = f.advance(&mut bodies, stop, controls()).unwrap();
    close(report.contact.spin_moment_nm, -capacity, 3e-17);
    assert_eq!(report.contact.force_n, [0.0; 2]);
    close(report.spin_after_rad_s, 0.0, 1e-16);
    close(
        report.dissipated_energy_j,
        0.5 * 0.3 * 0.3 / inverse_inertia,
        5e-18,
    );
    let common = 0.12 * 0.3 / (0.12 + 0.25);
    close(bodies[0].angular_momentum_body().z / 0.12, common, 1e-16);
    close(bodies[1].angular_momentum_body().z / 0.25, common, 1e-16);
    close(report.total_energy_residual_j, 0.0, 6e-18);
    assert!(norm(report.angular_momentum_change_nms) < 1e-17);
}

// Independent direct inverse-inertia bilinear form; no native Gram factor,
// effective mass API or returned mechanical matrix enters this reference.
fn mobility(f: &Fixture) -> [[f64; 3]; 3] {
    let axes = [
        vector(f.frame.longitudinal()),
        vector(f.frame.lateral()),
        vector(f.frame.normal()),
    ];
    let mut result = [[0.0; 3]; 3];
    for (body, arm) in [(0, f.binding.arm_a_body), (1, f.binding.arm_b_body)] {
        let q = f.bodies[body].pose().orientation();
        let torques = [
            arm.cross(q.rotate_world_to_body(axes[0])),
            arm.cross(q.rotate_world_to_body(axes[1])),
            q.rotate_world_to_body(axes[2]).scale(RADIUS),
        ];
        let inertia = f.properties[body].principal_inertia_body();
        for i in 0..3 {
            for j in 0..3 {
                let b = torques[j];
                result[i][j] +=
                    torques[i].dot(Vec3::new(b.x / inertia.x, b.y / inertia.y, b.z / inertia.z));
                if i < 2 && j < 2 {
                    result[i][j] += axes[i].dot(axes[j]) / f.properties[body].mass();
                }
            }
        }
    }
    result
}

// Manufacture midpoint slip/spin ratio epsilon=1, whose independent closed
// Hertz result is F/(mu*N)=9*pi/32 and M/(mu*N*a)=9*pi/128.
fn coupled(dt: f64) -> (Fixture, Vec3, f64) {
    let properties = [
        MassProperties::new(2.0, Vec3::ZERO, Vec3::new(0.09, 0.13, 0.17)).unwrap(),
        MassProperties::new(5.0, Vec3::ZERO, Vec3::new(0.20, 0.27, 0.31)).unwrap(),
    ];
    let orientations = [
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -1.0), 0.63).unwrap(),
        UnitQuaternion::from_axis_angle(Vec3::new(-2.0, 1.0, 1.0), -0.43).unwrap(),
    ];
    let arms = [Vec3::new(0.13, -0.08, -0.20), Vec3::new(-0.06, 0.11, 0.18)];
    let common = Vec3::new(0.4, -0.3, 0.2);
    let poses: [Pose; 2] = core::array::from_fn(|i| {
        Pose::new(
            common.sub(orientations[i].rotate_body_to_world(arms[i])),
            orientations[i],
        )
        .unwrap()
    });
    let mut f = Fixture {
        bodies: [
            state(poses[0], properties[0], Vec3::ZERO, Vec3::ZERO),
            state(poses[1], properties[1], Vec3::ZERO, Vec3::ZERO),
        ],
        properties,
        binding: RigidContactBinding {
            body_a: 0,
            body_b: 1,
            arm_a_body: arms[0],
            arm_b_body: arms[1],
        },
        frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
        law: law(),
    };
    let force = Vec3::new(-0.6, 0.8, 0.0).scale(4.0 * 9.0 * PI / 32.0);
    let moment = -4.0 * RADIUS * 9.0 * PI / 128.0;
    let response = mul(mobility(&f), [force.x, force.y, moment / RADIUS]);
    let initial = [
        0.3 - 0.5 * dt * response[0],
        -0.4 - 0.5 * dt * response[1],
        0.5 - 0.5 * dt * response[2],
    ];
    let omega = Vec3::new(0.0, 0.0, initial[2] / RADIUS);
    let velocity = Vec3::new(initial[0], initial[1], 0.0)
        .sub(omega.cross(orientations[0].rotate_body_to_world(arms[0])));
    f.bodies[0] = state(poses[0], properties[0], velocity, omega);
    (f, force, moment)
}

fn independent_endpoint(f: &Fixture, force: Vec3, moment: f64, dt: f64) -> [RigidBodyState; 2] {
    let mut after = f.bodies;
    for (i, arm, sign) in [
        (0, f.binding.arm_a_body, 1.0),
        (1, f.binding.arm_b_body, -1.0),
    ] {
        let body = f.bodies[i];
        let q = body.pose().orientation();
        let impulse = force.scale(sign * dt);
        let spin_couple = Vec3::new(0.0, 0.0, sign * moment * dt);
        let angular = arm
            .cross(q.rotate_world_to_body(impulse))
            .add(q.rotate_world_to_body(spin_couple));
        after[i] = RigidBodyState::new(
            body.pose(),
            body.linear_momentum_world().add(impulse),
            body.angular_momentum_body().add(angular),
        )
        .unwrap();
    }
    after
}
fn kinetic(bodies: &[RigidBodyState; 2], properties: &[MassProperties; 2]) -> f64 {
    bodies
        .iter()
        .zip(properties)
        .map(|(b, p)| {
            let linear = b.linear_momentum_world();
            let angular = b.angular_momentum_body();
            let inertia = p.principal_inertia_body();
            0.5 * (linear.dot(linear) / p.mass()
                + angular.x * angular.x / inertia.x
                + angular.y * angular.y / inertia.y
                + angular.z * angular.z / inertia.z)
        })
        .sum()
}

#[test]
fn full_anisotropic_three_axis_response_matches_independent_force_spin_and_endpoint_work() {
    let dt = 0.005;
    let (f, force, moment) = coupled(dt);
    let expected = independent_endpoint(&f, force, moment, dt);
    let expected_w = mobility(&f);
    assert!(expected_w[0][1].abs() > 1e-3);
    assert!(expected_w[0][2].abs() > 1e-3);
    assert!(expected_w[1][2].abs() > 1e-3);
    let mut bodies = f.bodies;
    let report = f.advance(&mut bodies, dt, controls()).unwrap();
    close(report.contact.force_n[0], force.x, 3e-9);
    close(report.contact.force_n[1], force.y, 3e-9);
    close(report.contact.spin_moment_nm, moment, 2e-10);
    for i in 0..3 {
        for j in 0..3 {
            close(report.scaled_mobility_per_kg[i][j], expected_w[i][j], 1e-15);
        }
    }
    for i in 0..2 {
        close_vec(
            bodies[i].linear_momentum_world(),
            expected[i].linear_momentum_world(),
            2e-11,
        );
        close_vec(
            bodies[i].angular_momentum_body(),
            expected[i].angular_momentum_body(),
            5e-12,
        );
        assert_eq!(bodies[i].pose(), f.bodies[i].pose());
    }
    close(report.actual_midpoint_displacement_m[0], 0.3 * dt, 2e-13);
    close(report.actual_midpoint_displacement_m[1], -0.4 * dt, 2e-13);
    close(
        report.actual_midpoint_spin_angle_rad,
        (0.5 / RADIUS) * dt,
        2e-12,
    );
    let work = dt * (force.dot(Vec3::new(0.3, -0.4, 0.0)) + moment * 0.5 / RADIUS);
    close(report.midpoint_impulse_work_j, work, 2e-11);
    close(report.dissipated_energy_j, -work, 2e-11);
    close(
        kinetic(&bodies, &f.properties) - kinetic(&f.bodies, &f.properties),
        work,
        2e-11,
    );
    assert!(report.native_force_residual_n <= controls().contact.force_tolerance_n);
    assert!(report.native_spin_moment_residual_nm <= controls().contact.moment_tolerance_nm);
    assert!(norm(report.linear_momentum_change_ns) < 1e-14);
    assert!(norm(report.angular_momentum_change_nms) < 1e-14);
}

fn compose(q: UnitQuaternion, r: UnitQuaternion) -> UnitQuaternion {
    let [a, b, c, d] = q.components();
    let [e, f, g, h] = r.components();
    UnitQuaternion::new(
        a * e - b * f - c * g - d * h,
        a * f + b * e + c * h - d * g,
        a * g - b * h + c * e + d * f,
        a * h + b * g - c * f + d * e,
    )
    .unwrap()
}

#[test]
fn nonprincipal_body_coordinates_preserve_hertz_force_spin_recoil_and_heat() {
    let dt = 0.005;
    let (f, force, moment) = coupled(dt);
    let expected = independent_endpoint(&f, force, moment, dt);
    let expected_mobility = mobility(&f);
    let labels = [
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -0.7), 0.81).unwrap(),
        UnitQuaternion::from_axis_angle(Vec3::new(-0.4, 0.7, 1.3), -0.63).unwrap(),
    ];
    let (mut dense, _, _) = coupled(dt);
    // Preserve the independently manufactured physical state while changing
    // each body's geometry frame. The original diagonal oracle remains
    // independent of the new full-tensor mobility implementation.
    for k in 0..2 {
        let [w, x, y, z] = labels[k].components();
        let inverse = UnitQuaternion::new(w, -x, -y, -z).unwrap();
        dense.properties[k] = MassProperties::from_principal_frame(
            f.properties[k].mass(),
            Vec3::ZERO,
            f.properties[k].principal_inertia_body(),
            labels[k],
        )
        .unwrap();
        dense.bodies[k] = RigidBodyState::new(
            Pose::new(
                f.bodies[k].pose().position_world(),
                compose(f.bodies[k].pose().orientation(), inverse),
            )
            .unwrap(),
            f.bodies[k].linear_momentum_world(),
            labels[k].rotate_body_to_world(f.bodies[k].angular_momentum_body()),
        )
        .unwrap();
        let m = dense.properties[k].inertia_tensor_body().m;
        assert!(m[1].abs() + m[2].abs() + m[5].abs() > 0.01);
    }
    dense.binding.arm_a_body = labels[0].rotate_body_to_world(f.binding.arm_a_body);
    dense.binding.arm_b_body = labels[1].rotate_body_to_world(f.binding.arm_b_body);
    let mut bodies = dense.bodies;
    let report = dense.advance(&mut bodies, dt, controls()).unwrap();
    close(report.contact.force_n[0], force.x, 3e-9);
    close(report.contact.force_n[1], force.y, 3e-9);
    close(report.contact.spin_moment_nm, moment, 2e-10);
    for i in 0..3 {
        for j in 0..3 {
            close(
                report.scaled_mobility_per_kg[i][j],
                expected_mobility[i][j],
                4e-15,
            );
        }
    }
    for k in 0..2 {
        close_vec(
            bodies[k].linear_momentum_world(),
            expected[k].linear_momentum_world(),
            2e-11,
        );
        close_vec(
            bodies[k]
                .pose()
                .orientation()
                .rotate_body_to_world(bodies[k].angular_momentum_body()),
            expected[k]
                .pose()
                .orientation()
                .rotate_body_to_world(expected[k].angular_momentum_body()),
            5e-12,
        );
        assert_eq!(bodies[k].pose(), dense.bodies[k].pose());
    }
    let work = dt * (force.dot(Vec3::new(0.3, -0.4, 0.0)) + moment * 0.5 / RADIUS);
    close(report.midpoint_impulse_work_j, work, 2e-11);
    close(report.dissipated_energy_j, -work, 2e-11);
    close(
        report.kinetic_after_j,
        kinetic(&expected, &f.properties),
        2e-11,
    );
    close(report.actual_midpoint_displacement_m[0], 0.3 * dt, 2e-13);
    close(report.actual_midpoint_displacement_m[1], -0.4 * dt, 2e-13);
    close(
        report.actual_midpoint_spin_angle_rad,
        (0.5 / RADIUS) * dt,
        2e-12,
    );
    assert!(norm(report.linear_momentum_change_ns) < 1e-14);
    assert!(norm(report.angular_momentum_change_nms) < 1e-14);
    assert!(report.total_energy_residual_j.abs() < 2e-11);
}

#[test]
fn world_covariance_common_boost_and_body_index_permutation_preserve_the_physical_solve() {
    let dt = 0.005;
    let (f, _, _) = coupled(dt);
    let mut original = f.bodies;
    let reference = f.advance(&mut original, dt, controls()).unwrap();
    let q = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -1.0, 2.0), 0.43).unwrap();
    let shift = Vec3::new(1.2, -0.4, 0.3);
    let boost = Vec3::new(0.3, -0.2, 0.1);
    let mut rotated = f.bodies;
    for i in 0..2 {
        rotated[i] = RigidBodyState::new(
            Pose::new(
                q.rotate_body_to_world(f.bodies[i].pose().position_world())
                    .add(shift),
                compose(q, f.bodies[i].pose().orientation()),
            )
            .unwrap(),
            q.rotate_body_to_world(f.bodies[i].linear_momentum_world())
                .add(boost.scale(f.properties[i].mass())),
            f.bodies[i].angular_momentum_body(),
        )
        .unwrap();
    }
    let frame = TangentFrame::new(
        array(q.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0))),
        array(q.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0))),
    )
    .unwrap();
    let report = apply_rigid_hertz(
        &mut rotated,
        &f.properties,
        RigidHertzContact {
            frame,
            ..f.contact()
        },
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    for i in 0..2 {
        close_vec(
            rotated[i].linear_momentum_world(),
            q.rotate_body_to_world(original[i].linear_momentum_world())
                .add(boost.scale(f.properties[i].mass())),
            3e-11,
        );
        close_vec(
            rotated[i].angular_momentum_body(),
            original[i].angular_momentum_body(),
            7e-12,
        );
    }
    close_vec(
        report.force_impulse_world_ns,
        q.rotate_body_to_world(reference.force_impulse_world_ns),
        3e-11,
    );
    close_vec(
        report.spin_angular_impulse_world_nms,
        q.rotate_body_to_world(reference.spin_angular_impulse_world_nms),
        2e-12,
    );
    close(
        report.dissipated_energy_j,
        reference.dissipated_energy_j,
        3e-11,
    );
    let mut swapped = [f.bodies[1], f.bodies[0]];
    let mut binding = f.binding;
    binding.body_a = 1;
    binding.body_b = 0;
    let reversed = apply_rigid_hertz(
        &mut swapped,
        &[f.properties[1], f.properties[0]],
        RigidHertzContact {
            binding,
            ..f.contact()
        },
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(swapped, [original[1], original[0]]);
    assert_eq!(reversed, reference);
}

#[test]
fn coarse_midpoint_reversal_late_angular_audit_and_cancellation_publish_nothing() {
    let f = axial(1.0, 0.0);
    let inverse_mass = 1.0 / 2.0 + 1.0 / 3.0 + 0.2 * 0.2 / 0.10 + 0.15 * 0.15 / 0.20;
    let stop = 1.0 / (4.0 * inverse_mass);
    let mut coarse = f.bodies;
    assert!(matches!(
        f.advance(&mut coarse, 1.5 * stop, controls()),
        Err(RigidHertzError::UnresolvedSlipEvent { .. })
    ));
    assert_eq!(coarse, f.bodies);
    let report = f.advance(&mut coarse, stop, controls()).unwrap();
    close(report.slip_velocity_after_m_s[0], 0.0, 5e-16);

    let dt = 0.005;
    let (mut mixed, _, _) = coupled(dt);
    let saved = mixed.bodies;
    let b = mixed.bodies[1];
    mixed.bodies[1] = RigidBodyState::new(
        Pose::new(
            b.pose().position_world().add(Vec3::new(0.0, 0.0, 1e-6)),
            b.pose().orientation(),
        )
        .unwrap(),
        b.linear_momentum_world(),
        b.angular_momentum_body(),
    )
    .unwrap();
    let mut candidate = mixed.bodies;
    let mut allowed_geometry = controls();
    allowed_geometry.geometry_tolerance_m = 2e-6;
    assert!(matches!(
        mixed.advance(&mut candidate, dt, allowed_geometry),
        Err(RigidHertzError::Audit(_))
    ));
    assert_eq!(candidate, mixed.bodies);
    mixed.bodies = saved;
    candidate = saved;
    let mut limited = controls();
    limited.contact.maximum_force_n = 0.01;
    assert!(matches!(
        mixed.advance(&mut candidate, dt, limited),
        Err(RigidHertzError::Native(
            HertzWrenchError::EffortBudget { .. }
        ))
    ));
    assert_eq!(candidate, saved);
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        apply_rigid_hertz(
            &mut candidate,
            &mixed.properties,
            mixed.contact(),
            dt,
            controls(),
            &gate
        ),
        Err(RigidHertzError::Cancelled)
    ));
    assert_eq!(candidate, saved);
    let clean = mixed.advance(&mut candidate, dt, controls()).unwrap();
    let mut retry = saved;
    assert_eq!(mixed.advance(&mut retry, dt, controls()).unwrap(), clean);
    assert_eq!(retry, candidate);
}

#[test]
fn exact_rest_and_zero_load_do_not_create_heat_or_discard_any_body_momentum() {
    let rest = axial(0.0, 0.0);
    let mut bodies = rest.bodies;
    let report = rest.advance(&mut bodies, 0.1, controls()).unwrap();
    assert_eq!(bodies, rest.bodies);
    assert_eq!(report.contact.mode, HertzWrenchMode::Sticking);
    assert_eq!(report.dissipated_energy_j, 0.0);
    let (mut f, _, _) = coupled(0.005);
    let mut parameters = f.law.parameters();
    parameters.normal_force_n = 0.0;
    f.law = HertzWrenchContact::new(parameters).unwrap();
    bodies = f.bodies;
    let report = f.advance(&mut bodies, 0.005, controls()).unwrap();
    assert_eq!(bodies, f.bodies);
    assert_eq!(report.contact.mode, HertzWrenchMode::Frictionless);
    assert_eq!(report.dissipated_energy_j, 0.0);
    assert_eq!(report.force_impulse_world_ns, Vec3::ZERO);
    assert_eq!(report.spin_angular_impulse_world_nms, Vec3::ZERO);
}

#[test]
fn prepared_anisotropic_pair_reuses_the_original_snapshot_and_matches_native_publication() {
    let dt = 0.005;
    let (f, force, moment) = coupled(dt);
    let expected = independent_endpoint(&f, force, moment, dt);
    let expected_mobility = mobility(&f);
    let mut original = f.bodies;
    let prepared = PreparedRigidHertzPair::new(
        &original,
        &f.properties,
        f.contact(),
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let first = prepared.advance(&f.law, &CancelGate::new()).unwrap();
    assert_eq!(original, f.bodies);
    for i in 0..2 {
        close_vec(
            first.bodies_after[i].linear_momentum_world(),
            expected[i].linear_momentum_world(),
            2e-11,
        );
        close_vec(
            first.bodies_after[i].angular_momentum_body(),
            expected[i].angular_momentum_body(),
            5e-12,
        );
        close_vec(
            prepared.contact_points_world()[i],
            Vec3::new(0.4, -0.3, 0.2),
            2e-16,
        );
    }
    for i in 0..3 {
        for j in 0..3 {
            close(
                prepared.scaled_mobility_per_kg()[i][j],
                expected_mobility[i][j],
                1e-15,
            );
            close(
                prepared.scaled_compliance_m_per_n()[i][j],
                0.5 * dt * dt * expected_mobility[i][j],
                1e-20,
            );
        }
    }
    assert_eq!(prepared.binding().body_a, f.binding.body_a);
    assert_eq!(prepared.binding().body_b, f.binding.body_b);
    assert_eq!(prepared.binding().arm_a_body, f.binding.arm_a_body);
    assert_eq!(prepared.binding().arm_b_body, f.binding.arm_b_body);
    assert_eq!(prepared.frame(), f.frame);
    assert_eq!(prepared.duration_s(), dt);
    assert_eq!(prepared.normal_force_n(), 10.0);
    assert_eq!(prepared.contact_radius_m(), RADIUS);
    assert_eq!(prepared.controls().contact, controls().contact);

    // Caller containers are not retained: later changes cannot mutate the
    // admitted native snapshot or silently alter its cached mechanical response.
    original[0] = RigidBodyState::new(original[0].pose(), Vec3::ZERO, Vec3::ZERO).unwrap();
    assert_eq!(prepared.advance(&f.law, &CancelGate::new()).unwrap(), first);
    let mut published = f.bodies;
    let report = f.advance(&mut published, dt, controls()).unwrap();
    assert_eq!(published, first.bodies_after);
    assert_eq!(report, first.report);
}

#[test]
fn prepared_coefficient_trials_keep_one_interval_and_match_independent_finite_recoil_heat() {
    let f = axial(1.0, 0.0);
    let dt = 0.02;
    let prepared = PreparedRigidHertzPair::new(
        &f.bodies,
        &f.properties,
        f.contact(),
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let w = 1.0 / 2.0 + 1.0 / 3.0 + 0.2 * 0.2 / 0.10 + 0.15 * 0.15 / 0.20;
    let saved = f.bodies;
    let mut first = None;
    for (i, mu) in [0.4, 0.2, 0.5, 0.4].into_iter().enumerate() {
        let law = HertzWrenchContact::new(HertzWrenchParameters {
            static_mu: mu + 0.1,
            kinetic_mu: mu,
            ..f.law.parameters()
        })
        .unwrap();
        let trial = prepared.advance(&law, &CancelGate::new()).unwrap();
        let force = 10.0 * mu;
        close(trial.report.contact.force_n[0], -force, 1e-15);
        close(
            trial.report.slip_velocity_after_m_s[0],
            1.0 - dt * w * force,
            2e-15,
        );
        close(
            trial.report.dissipated_energy_j,
            force * dt * (1.0 - 0.5 * dt * w * force),
            2e-15,
        );
        close(
            trial.bodies_after[0].linear_momentum_world().x,
            2.0 - dt * force,
            1e-15,
        );
        close(
            trial.bodies_after[1].linear_momentum_world().x,
            dt * force,
            1e-15,
        );
        assert_eq!(
            trial.report.scaled_mobility_per_kg,
            *prepared.scaled_mobility_per_kg()
        );
        assert_eq!(
            trial.report.scaled_compliance_m_per_n,
            *prepared.scaled_compliance_m_per_n()
        );
        if i == 0 {
            first = Some(trial);
        } else if i == 3 {
            assert_eq!(Some(trial), first);
        }
        assert_eq!(f.bodies, saved);
    }
}

#[test]
fn prepared_patch_effort_event_and_cancellation_refusals_leave_exact_retry() {
    let f = axial(1.0, 0.0);
    let prepared = PreparedRigidHertzPair::new(
        &f.bodies,
        &f.properties,
        f.contact(),
        0.1,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let gentle = HertzWrenchContact::new(HertzWrenchParameters {
        static_mu: 0.3,
        kinetic_mu: 0.2,
        ..f.law.parameters()
    })
    .unwrap();
    let expected = prepared.advance(&gentle, &CancelGate::new()).unwrap();
    for parameters in [
        HertzWrenchParameters {
            normal_force_n: 11.0,
            ..gentle.parameters()
        },
        HertzWrenchParameters {
            contact_radius_m: 1.01 * RADIUS,
            ..gentle.parameters()
        },
    ] {
        assert!(matches!(
            prepared.advance(
                &HertzWrenchContact::new(parameters).unwrap(),
                &CancelGate::new()
            ),
            Err(RigidHertzError::Invalid(
                "prepared Hertz normal force and radius must remain unchanged"
            ))
        ));
    }
    // The same physical interval is too coarse at the stronger coefficients.
    // An admitted midpoint-static wrench must not bypass endpoint reversal.
    let strong = HertzWrenchContact::new(HertzWrenchParameters {
        static_mu: 1.5,
        kinetic_mu: 1.0,
        ..f.law.parameters()
    })
    .unwrap();
    assert!(matches!(
        prepared.advance(&strong, &CancelGate::new()),
        Err(RigidHertzError::UnresolvedSlipEvent { .. })
    ));
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        prepared.advance(&gentle, &gate),
        Err(RigidHertzError::Cancelled)
    ));
    assert_eq!(
        prepared.advance(&gentle, &CancelGate::new()).unwrap(),
        expected
    );

    let mut bounded = controls();
    bounded.contact.maximum_force_n = 3.0;
    let bounded = PreparedRigidHertzPair::new(
        &f.bodies,
        &f.properties,
        f.contact(),
        0.1,
        bounded,
        &CancelGate::new(),
    )
    .unwrap();
    let before = bounded.advance(&gentle, &CancelGate::new()).unwrap();
    assert!(matches!(
        bounded.advance(&f.law, &CancelGate::new()),
        Err(RigidHertzError::Native(
            HertzWrenchError::EffortBudget { .. }
        ))
    ));
    assert_eq!(
        bounded.advance(&gentle, &CancelGate::new()).unwrap(),
        before
    );
    assert_eq!(f.bodies, axial(1.0, 0.0).bodies);
}
