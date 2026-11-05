//! G1/G3/G4: finite Hertz impulses through original full-inertia native bodies.

use core::f64::consts::PI;
use fs_ga::Mat3;
use fs_mbd::frictional_impact::RigidImpactPoint;
use fs_mbd::frictional_impact::hertz_patch::{
    PreparedRigidHertzImpulse, RigidHertzImpulseContact, RigidHertzImpulseControls,
    RigidHertzImpulseCounterpart, RigidHertzImpulseError, RigidHertzImpulseReport,
    apply_rigid_hertz_impulse,
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_tribo::sliding_patch::contact::HertzWrenchMode;
use fs_tribo::sliding_patch::contact::impulse::{
    HertzImpulseControls, HertzImpulseError, HertzImpulsePatch,
};

const X: Vec3 = Vec3::new(1.0, 0.0, 0.0);
const Y: Vec3 = Vec3::new(0.0, 1.0, 0.0);
const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance:e}"
    );
}
fn vector(actual: Vec3, expected: Vec3, tolerance: f64) {
    near(actual.x, expected.x, tolerance);
    near(actual.y, expected.y, tolerance);
    near(actual.z, expected.z, tolerance);
}
fn properties(mass: f64, inertia: f64) -> MassProperties {
    MassProperties::new(mass, Vec3::ZERO, Vec3::new(inertia, inertia, inertia)).unwrap()
}
fn state(property: MassProperties, position: Vec3, velocity: Vec3, omega: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(position, UnitQuaternion::IDENTITY).unwrap(),
        velocity.scale(property.mass()),
        property.angular_momentum_body_checked(omega).unwrap(),
    )
    .unwrap()
}
fn controls() -> RigidHertzImpulseControls {
    RigidHertzImpulseControls {
        maximum_bodies: 8,
        geometry_tolerance_m: 1e-10,
        velocity_tolerance_m_s: 3e-9,
        spin_tolerance_rad_s: 3e-9,
        linear_momentum_tolerance_ns: 2e-9,
        angular_momentum_tolerance_nms: 2e-9,
        energy_tolerance_j: 2e-8,
        contact: HertzImpulseControls {
            maximum_static_iterations: 96,
            maximum_iterations: 256,
            maximum_evaluations: 2048,
            impulse_tolerance_ns: 2e-9,
            spin_impulse_tolerance_nms: 5e-10,
            velocity_tolerance_m_s: 2e-13,
            energy_tolerance_j: 2e-10,
            relative_matrix_tolerance: 1e-13,
            maximum_impulse_ns: 20.0,
            maximum_spin_impulse_nms: 10.0,
        },
    }
}
fn pair_contact(law: &HertzImpulsePatch, arms: [Vec3; 2]) -> RigidHertzImpulseContact<'_> {
    RigidHertzImpulseContact {
        point: RigidImpactPoint {
            body: 0,
            arm_body: arms[0],
        },
        other: RigidHertzImpulseCounterpart::Body(RigidImpactPoint {
            body: 1,
            arm_body: arms[1],
        }),
        normal_world: Y,
        tangent_world: X,
        law,
    }
}
fn support_contact(
    law: &HertzImpulsePatch,
    arm: Vec3,
    velocity: Vec3,
    omega: Vec3,
) -> RigidHertzImpulseContact<'_> {
    RigidHertzImpulseContact {
        point: RigidImpactPoint {
            body: 0,
            arm_body: arm,
        },
        other: RigidHertzImpulseCounterpart::PrescribedPatch {
            point_world: Vec3::ZERO,
            velocity_world: velocity,
            angular_velocity_world: omega,
        },
        normal_world: Y,
        tangent_world: X,
        law,
    }
}
fn balance(report: &RigidHertzImpulseReport) {
    assert!(report.dissipated_energy_j >= 0.0);
    near(report.work_energy_residual_j, 0.0, 2e-8);
    near(report.total_energy_residual_j, 0.0, 2e-8);
    near(report.contact_work_residual_j, 0.0, 2e-8);
    vector(report.linear_momentum_residual_ns, Vec3::ZERO, 2e-9);
    vector(report.angular_momentum_residual_nms, Vec3::ZERO, 2e-9);
    near(report.endpoint_tangent_residual_m_s, 0.0, 3e-9);
    near(report.endpoint_spin_residual_rad_s, 0.0, 3e-9);
    near(report.normal_response_residual_m_s, 0.0, 3e-9);
}

#[test]
fn g1_finite_translation_recoils_both_bodies_and_arrests_without_rebound() {
    let p = [properties(2.0, 0.8), properties(3.0, 1.2)];
    let initial = [
        state(p[0], Y, X, Vec3::ZERO),
        state(p[1], Y.scale(-1.0), X.scale(-1.0), Vec3::ZERO),
    ];
    let patch = HertzImpulsePatch::new(0.25, 0.5).unwrap();
    let contact = pair_contact(&patch, [Y.scale(-1.0), Y]);
    // Exact scalar response: 1/2 + 1/3 + 1/.8 + 1/1.2 = 35/12.
    for (normal, expected_impulse, expected_mode) in [
        (4.0, -24.0 / 35.0, HertzWrenchMode::Sticking),
        (0.5, -0.25, HertzWrenchMode::Sliding),
    ] {
        let mut bodies = initial;
        let report =
            apply_rigid_hertz_impulse(&mut bodies, &p, contact, normal, controls(), || false)
                .unwrap();
        let after_slip = 2.0 + (35.0 / 12.0) * expected_impulse;
        vector(
            report.force_impulse_world_ns,
            X.scale(expected_impulse),
            2e-13,
        );
        vector(report.spin_impulse_world_nms, Vec3::ZERO, 1e-14);
        assert_eq!(report.contact.mode, expected_mode);
        near(report.slip_velocity_after_m_s[0], after_slip, 3e-13);
        near(
            report.dissipated_energy_j,
            -expected_impulse * (2.0 + after_slip) / 2.0,
            2e-13,
        );
        vector(
            bodies[0].linear_momentum_world(),
            X.scale(2.0 + expected_impulse),
            2e-13,
        );
        vector(
            bodies[1].linear_momentum_world(),
            X.scale(-3.0 - expected_impulse),
            2e-13,
        );
        vector(
            bodies[0].angular_momentum_body(),
            Z.scale(expected_impulse),
            2e-13,
        );
        vector(
            bodies[1].angular_momentum_body(),
            Z.scale(expected_impulse),
            2e-13,
        );
        assert_eq!(bodies[0].pose(), initial[0].pose());
        assert_eq!(bodies[1].pose(), initial[1].pose());
        assert!(report.counterpart_impulse.is_some());
        balance(&report);
    }
}

#[test]
fn g1_spin_impulse_uses_finite_rotor_recoil_and_the_hertz_radial_moment() {
    let p = [properties(2.0, 0.8), properties(3.0, 1.2)];
    let initial = [
        state(p[0], Y, Vec3::ZERO, Y.scale(4.0)),
        state(p[1], Y.scale(-1.0), Vec3::ZERO, Y.scale(-1.0)),
    ];
    let patch = HertzImpulsePatch::new(0.25, 0.5).unwrap();
    let contact = pair_contact(&patch, [Y.scale(-1.0), Y]);
    for normal in [64.0, 16.0] {
        let mut bodies = initial;
        let report =
            apply_rigid_hertz_impulse(&mut bodies, &p, contact, normal, controls(), || false)
                .unwrap();
        let capacity = (3.0 * PI / 16.0) * 0.5 * normal * 0.25;
        let expected = -capacity.min(2.4);
        let relative_after = 5.0 + (1.0 / 0.8 + 1.0 / 1.2) * expected;
        vector(report.force_impulse_world_ns, Vec3::ZERO, 1e-14);
        vector(report.spin_impulse_world_nms, Y.scale(expected), 3e-13);
        near(report.spin_after_rad_s, relative_after, 3e-13);
        near(
            report.dissipated_energy_j,
            -expected * (5.0 + relative_after) / 2.0,
            3e-13,
        );
        if normal == 64.0 {
            near(report.dissipated_energy_j, 6.0, 3e-13);
            vector(
                p[0].angular_velocity_body(bodies[0].angular_momentum_body()),
                Y,
                3e-13,
            );
            vector(
                p[1].angular_velocity_body(bodies[1].angular_momentum_body()),
                Y,
                3e-13,
            );
        }
        balance(&report);
    }
}

#[test]
fn g1_moving_and_spinning_supports_supply_or_recover_signed_work() {
    let property = properties(2.0, 0.8);
    let patch = HertzImpulsePatch::new(0.25, 1.0).unwrap();
    for (velocity, impulse, kinetic, support) in [
        (0.0, 4.0 / 7.0, 2.0 / 7.0, 4.0 / 7.0),
        (2.0, -4.0 / 7.0, -6.0 / 7.0, -4.0 / 7.0),
    ] {
        let mut bodies = [state(property, Y, X.scale(velocity), Vec3::ZERO)];
        let contact = support_contact(&patch, Y.scale(-1.0), X, Vec3::ZERO);
        let report =
            apply_rigid_hertz_impulse(&mut bodies, &[property], contact, 2.0, controls(), || false)
                .unwrap();
        vector(report.force_impulse_world_ns, X.scale(impulse), 2e-13);
        near(report.kinetic_change_j, kinetic, 2e-13);
        near(report.support_work_j, support, 2e-13);
        near(report.dissipated_energy_j, 2.0 / 7.0, 2e-13);
        vector(report.support_impulse_ns, X.scale(impulse), 2e-13);
        vector(report.support_angular_impulse_nms, Vec3::ZERO, 1e-14);
        assert!(report.counterpart_impulse.is_none());
        balance(&report);
    }
    for (spin, expected_impulse, expected_kinetic) in [(0.0, 1.6, 1.6), (4.0, -1.6, -4.8)] {
        let mut bodies = [state(property, Y, Vec3::ZERO, Y.scale(spin))];
        let contact = support_contact(&patch, Y.scale(-1.0), Vec3::ZERO, Y.scale(2.0));
        let report =
            apply_rigid_hertz_impulse(&mut bodies, &[property], contact, 32.0, controls(), || {
                false
            })
            .unwrap();
        vector(
            report.spin_impulse_world_nms,
            Y.scale(expected_impulse),
            3e-13,
        );
        near(report.kinetic_change_j, expected_kinetic, 3e-13);
        near(report.support_work_j, 2.0 * expected_impulse, 3e-13);
        near(report.dissipated_energy_j, 1.6, 3e-13);
        near(report.spin_after_rad_s, 0.0, 3e-13);
        balance(&report);
    }
}

fn inverse(a: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut cofactors = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            cofactors[i][j] = a[(i + 1) % 3][(j + 1) % 3] * a[(i + 2) % 3][(j + 2) % 3]
                - a[(i + 1) % 3][(j + 2) % 3] * a[(i + 2) % 3][(j + 1) % 3];
        }
    }
    let determinant = (0..3).map(|j| a[0][j] * cofactors[0][j]).sum::<f64>();
    core::array::from_fn(|i| core::array::from_fn(|j| cofactors[j][i] / determinant))
}
fn multiply(matrix: [[f64; 3]; 3], vector: Vec3) -> Vec3 {
    let v = [vector.x, vector.y, vector.z];
    let result = matrix.map(|row| (0..3).map(|j| row[j] * v[j]).sum::<f64>());
    Vec3::new(result[0], result[1], result[2])
}
fn tensor_fixture() -> (RigidBodyState, MassProperties, [[f64; 3]; 3], Vec3) {
    let inertia = [[2.0, 0.4, 0.0], [0.4, 3.0, 0.3], [0.0, 0.3, 4.0]];
    let property = MassProperties::from_inertia_tensor(
        2.0,
        Vec3::ZERO,
        Mat3 {
            m: [2.0, 0.4, 0.0, 0.4, 3.0, 0.3, 0.0, 0.3, 4.0],
        },
    )
    .unwrap();
    let arm = Vec3::new(0.25, -0.5, 0.125);
    // An independently authored point impulse plus free couple generates the
    // entire incoming state; its negative must arrest all six body components.
    let body = RigidBodyState::new(
        Pose::new(arm.scale(-1.0), UnitQuaternion::IDENTITY).unwrap(),
        Vec3::new(0.5, 0.0, 0.25),
        Vec3::new(-0.125, 0.03125, 0.25),
    )
    .unwrap();
    (body, property, inertia, arm)
}

/// Independent local pressure integration at the actual outgoing velocity.
/// r=a sin(phi) makes the Hertz square-root edge smooth for Simpson quadrature.
fn pressure_impulses(v: [f64; 2], spin: f64, radius: f64, capacity: f64) -> [f64; 3] {
    let radial = 128;
    let angular = 256;
    let dphi = (PI / 2.0) / radial as f64;
    let dtheta = 2.0 * PI / angular as f64;
    let mut result = [0.0; 3];
    for i in 0..=radial {
        let phi = i as f64 * dphi;
        let r = radius * phi.sin();
        let simpson = if i == 0 || i == radial {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        let weight = capacity * 3.0 / (2.0 * PI) * phi.sin() * phi.cos().powi(2) * simpson * dphi
            / 3.0
            * dtheta;
        for j in 0..angular {
            let theta = (j as f64 + 0.5) * dtheta;
            let x = r * theta.cos();
            let y = r * theta.sin();
            let vx = v[0] - spin * y;
            let vy = v[1] + spin * x;
            let speed = vx.hypot(vy);
            let px = -weight * vx / speed;
            let py = -weight * vy / speed;
            result[0] += px;
            result[1] += py;
            result[2] += x * py - y * px;
        }
    }
    result
}

#[test]
fn g1_full_tensor_mixed_impulse_matches_independent_mechanics_and_pressure() {
    let (initial, property, inertia, arm) = tensor_fixture();
    let inverse_inertia = inverse(inertia);
    let patch = HertzImpulsePatch::new(0.25, 1.0).unwrap();
    let contact = support_contact(&patch, arm, Vec3::ZERO, Vec3::ZERO);
    let prepared =
        PreparedRigidHertzImpulse::new(&[initial], &[property], contact, controls(), || false)
            .unwrap();
    let force_axes = [X, Z.scale(-1.0), Vec3::ZERO];
    let torque_axes = [
        arm.cross(force_axes[0]),
        arm.cross(force_axes[1]),
        Y.scale(0.25),
    ];
    for i in 0..3 {
        for j in 0..3 {
            let expected = force_axes[i].dot(force_axes[j]) / 2.0
                + torque_axes[i].dot(multiply(inverse_inertia, torque_axes[j]));
            near(prepared.scaled_mobility_per_kg()[i][j], expected, 2e-12);
        }
        let expected = arm.cross(Y).dot(multiply(inverse_inertia, torque_axes[i]));
        near(prepared.normal_cross_mobility_per_kg()[i], expected, 2e-12);
    }
    assert!(prepared.scaled_mobility_per_kg()[0][2].abs() > 0.001);
    for normal in [2.0, 0.2] {
        let mut bodies = [initial];
        let report = apply_rigid_hertz_impulse(
            &mut bodies,
            &[property],
            contact,
            normal,
            controls(),
            || false,
        )
        .unwrap();
        let impulse = report.force_impulse_world_ns;
        let couple = report.spin_impulse_world_nms;
        let expected_p = initial.linear_momentum_world().add(impulse);
        let expected_l = initial
            .angular_momentum_body()
            .add(arm.cross(impulse))
            .add(couple);
        vector(bodies[0].linear_momentum_world(), expected_p, 2e-14);
        vector(bodies[0].angular_momentum_body(), expected_l, 2e-14);
        let omega_before = multiply(inverse_inertia, initial.angular_momentum_body());
        let omega_after = multiply(inverse_inertia, expected_l);
        let velocity_after = expected_p.scale(0.5).add(omega_after.cross(arm));
        near(report.slip_velocity_after_m_s[0], velocity_after.x, 2e-12);
        near(report.slip_velocity_after_m_s[1], -velocity_after.z, 2e-12);
        near(report.spin_after_rad_s, omega_after.y, 2e-12);
        near(report.normal_velocity_after_m_s, velocity_after.y, 2e-12);
        let kinetic_before = initial.linear_momentum_world().norm_squared() / 4.0
            + 0.5 * initial.angular_momentum_body().dot(omega_before);
        let kinetic_after = expected_p.norm_squared() / 4.0 + 0.5 * expected_l.dot(omega_after);
        near(
            report.kinetic_change_j,
            kinetic_after - kinetic_before,
            2e-12,
        );
        near(
            report.dissipated_energy_j,
            kinetic_before - kinetic_after,
            2e-12,
        );
        assert!(report.normal_velocity_before_m_s.abs() > 0.01);
        if normal == 2.0 {
            assert_eq!(report.contact.mode, HertzWrenchMode::Sticking);
            vector(bodies[0].linear_momentum_world(), Vec3::ZERO, 2e-12);
            vector(bodies[0].angular_momentum_body(), Vec3::ZERO, 2e-12);
            near(report.normal_velocity_after_m_s, 0.0, 2e-12);
        } else {
            assert_eq!(report.contact.mode, HertzWrenchMode::Sliding);
            let reference = pressure_impulses(
                report.slip_velocity_after_m_s,
                report.spin_after_rad_s,
                0.25,
                normal,
            );
            near(report.contact.impulse_ns[0], reference[0], 3e-8);
            near(report.contact.impulse_ns[1], reference[1], 3e-8);
            near(report.contact.spin_impulse_nms, reference[2], 3e-8);
        }
        balance(&report);
    }
}

#[test]
fn g3_native_result_is_covariant_under_frame_boost_and_endpoint_reversal() {
    let (_, first_property, _, arm_a) = tensor_fixture();
    let second_property = MassProperties::from_principal_frame(
        3.0,
        Vec3::ZERO,
        Vec3::new(1.5, 2.0, 2.5),
        UnitQuaternion::from_axis_angle(Vec3::new(2.0, -1.0, 1.0), 0.4).unwrap(),
    )
    .unwrap();
    let p = [first_property, second_property];
    let arm_b = Vec3::new(-0.375, 0.25, 0.125);
    let initial = [
        RigidBodyState::new(
            Pose::new(arm_a.scale(-1.0), UnitQuaternion::IDENTITY).unwrap(),
            Vec3::new(1.0, 0.5, -0.5),
            Vec3::new(0.2, -0.1, 0.3),
        )
        .unwrap(),
        RigidBodyState::new(
            Pose::new(arm_b.scale(-1.0), UnitQuaternion::IDENTITY).unwrap(),
            Vec3::new(-0.2, 0.4, 0.8),
            Vec3::new(-0.05, 0.15, -0.2),
        )
        .unwrap(),
    ];
    let patch = HertzImpulsePatch::new(0.125, 0.4).unwrap();
    let contact = pair_contact(&patch, [arm_a, arm_b]);
    let mut baseline = initial;
    let reference =
        apply_rigid_hertz_impulse(&mut baseline, &p, contact, 2.0, controls(), || false).unwrap();
    let mut reversed = initial;
    let mut swapped = contact;
    swapped.point = RigidImpactPoint {
        body: 1,
        arm_body: arm_b,
    };
    swapped.other = RigidHertzImpulseCounterpart::Body(contact.point);
    swapped.normal_world = Y.scale(-1.0);
    let reverse =
        apply_rigid_hertz_impulse(&mut reversed, &p, swapped, 2.0, controls(), || false).unwrap();
    for i in 0..2 {
        vector(
            reversed[i].linear_momentum_world(),
            baseline[i].linear_momentum_world(),
            5e-10,
        );
        vector(
            reversed[i].angular_momentum_body(),
            baseline[i].angular_momentum_body(),
            5e-10,
        );
    }
    vector(
        reverse.force_impulse_world_ns,
        reference.force_impulse_world_ns.scale(-1.0),
        5e-10,
    );
    vector(
        reverse.spin_impulse_world_nms,
        reference.spin_impulse_world_nms.scale(-1.0),
        5e-10,
    );
    near(
        reverse.dissipated_energy_j,
        reference.dissipated_energy_j,
        5e-10,
    );
    balance(&reverse);

    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(2.0, -1.0, 3.0), 0.7).unwrap();
    let shift = Vec3::new(3.0, -7.0, 11.0);
    let boost = Vec3::new(8.0, -16.0, 32.0);
    let mut moved = core::array::from_fn::<_, 2, _>(|i| {
        RigidBodyState::new(
            Pose::new(
                rotation
                    .rotate_body_to_world(initial[i].pose().position_world())
                    .add(shift),
                rotation,
            )
            .unwrap(),
            rotation
                .rotate_body_to_world(initial[i].linear_momentum_world())
                .add(boost.scale(p[i].mass())),
            initial[i].angular_momentum_body(),
        )
        .unwrap()
    });
    let mut rotated_contact = contact;
    rotated_contact.normal_world = rotation.rotate_body_to_world(Y);
    rotated_contact.tangent_world = rotation.rotate_body_to_world(X);
    let rotated =
        apply_rigid_hertz_impulse(&mut moved, &p, rotated_contact, 2.0, controls(), || false)
            .unwrap();
    for i in 0..2 {
        vector(
            moved[i].linear_momentum_world(),
            rotation
                .rotate_body_to_world(baseline[i].linear_momentum_world())
                .add(boost.scale(p[i].mass())),
            2e-9,
        );
        vector(
            moved[i].angular_momentum_body(),
            baseline[i].angular_momentum_body(),
            2e-9,
        );
    }
    vector(
        rotated.force_impulse_world_ns,
        rotation.rotate_body_to_world(reference.force_impulse_world_ns),
        2e-9,
    );
    vector(
        rotated.spin_impulse_world_nms,
        rotation.rotate_body_to_world(reference.spin_impulse_world_nms),
        2e-9,
    );
    near(rotated.kinetic_change_j, reference.kinetic_change_j, 3e-9);
    near(
        rotated.dissipated_energy_j,
        reference.dissipated_energy_j,
        3e-9,
    );
    balance(&rotated);
}

#[test]
fn g3_zero_limits_and_reusable_normal_budget_trials_keep_native_identity() {
    let (first, p0, _, arm) = tensor_fixture();
    let second_property = properties(3.0, 1.2);
    let second = state(second_property, Y.scale(-1.0), X.scale(-0.3), Z.scale(0.2));
    let untouched = state(properties(2.0, 1.0), Vec3::new(5.0, 3.0, 2.0), Y, X);
    let initial = [first, second, untouched];
    let p = [p0, second_property, properties(2.0, 1.0)];
    let patch = HertzImpulsePatch::new(0.25, 0.5).unwrap();
    let contact = pair_contact(&patch, [arm, Y]);
    let prepared =
        PreparedRigidHertzImpulse::new(&initial, &p, contact, controls(), || false).unwrap();
    for (normal, coefficient) in [(0.0, 0.5), (2.0, 0.0)] {
        let zero_law = HertzImpulsePatch::new(0.25, coefficient).unwrap();
        let mut c = contact;
        c.law = &zero_law;
        let mut bodies = initial;
        let report =
            apply_rigid_hertz_impulse(&mut bodies, &p, c, normal, controls(), || false).unwrap();
        assert_eq!(bodies, initial);
        assert_eq!(report.kinetic_change_j, 0.0);
        assert_eq!(report.dissipated_energy_j, 0.0);
        balance(&report);
    }
    let changed_coefficient = HertzImpulsePatch::new(0.25, 0.3).unwrap();
    for normal in [0.25, 2.0] {
        let a = prepared
            .solve(&changed_coefficient, normal, || false)
            .unwrap();
        let b = prepared
            .solve(&changed_coefficient, normal, || false)
            .unwrap();
        assert_eq!(a, b);
        let mut c = contact;
        c.law = &changed_coefficient;
        let mut bodies = initial;
        let fresh =
            apply_rigid_hertz_impulse(&mut bodies, &p, c, normal, controls(), || false).unwrap();
        assert_eq!(a, fresh);
        assert_eq!(bodies[2], untouched);
        assert_eq!(initial[0], first);
    }
    let wrong_radius = HertzImpulsePatch::new(0.5, 0.3).unwrap();
    assert!(matches!(
        prepared.solve(&wrong_radius, 2.0, || false),
        Err(RigidHertzImpulseError::Invalid(_))
    ));
}

#[test]
fn g4_every_cancellation_and_late_native_refusal_is_atomic() {
    let property = properties(2.0, 0.8);
    let patch = HertzImpulsePatch::new(0.25, 1.0).unwrap();
    let contact = support_contact(&patch, Y.scale(-1.0), X, Vec3::ZERO);
    let initial = [state(property, Y, Vec3::ZERO, Vec3::ZERO)];
    let mut calls = 0;
    apply_rigid_hertz_impulse(
        &mut initial.clone(),
        &[property],
        contact,
        2.0,
        controls(),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    assert!(calls > 10);
    for cancel_at in 1..=calls {
        let mut bodies = initial;
        let mut count = 0;
        let result =
            apply_rigid_hertz_impulse(&mut bodies, &[property], contact, 2.0, controls(), || {
                count += 1;
                count == cancel_at
            });
        assert!(
            matches!(result, Err(RigidHertzImpulseError::Cancelled)),
            "{cancel_at}: {result:?}"
        );
        assert_eq!(bodies, initial);
    }
    let mut limited = controls();
    limited.contact.maximum_impulse_ns = 1e-6;
    let mut bodies = initial;
    assert!(matches!(
        apply_rigid_hertz_impulse(&mut bodies, &[property], contact, 2.0, limited, || false),
        Err(RigidHertzImpulseError::Contact(
            HertzImpulseError::ImpulseBudget { .. }
        ))
    ));
    assert_eq!(bodies, initial);

    // A small requested impulse cannot be hidden inside a much larger stored
    // COM momentum, even if its opposite reaction or total energies cancel.
    let huge = 2.0_f64.powi(54);
    let p = properties(1.0, 1.0);
    let impossible = [state(p, Y, X.scale(huge), Vec3::ZERO)];
    let law = HertzImpulsePatch::new(0.25, 0.25).unwrap();
    let c = support_contact(&law, Y.scale(-1.0), X.scale(huge - 4.0), Vec3::ZERO);
    let mut bodies = impossible;
    assert!(matches!(
        apply_rigid_hertz_impulse(&mut bodies, &[p], c, 1.0, controls(), || false),
        Err(RigidHertzImpulseError::Audit(_))
    ));
    assert_eq!(bodies, impossible);

    // Rounded absolute points coincide, but their original material offset is
    // genuinely one metre. Preparation must retain that mismatch.
    let far = [state(p, X.scale(huge), Vec3::ZERO, Vec3::ZERO)];
    let mut c = support_contact(&law, X, Vec3::ZERO, Vec3::ZERO);
    c.other = RigidHertzImpulseCounterpart::PrescribedPatch {
        point_world: X.scale(huge),
        velocity_world: Vec3::ZERO,
        angular_velocity_world: Vec3::ZERO,
    };
    let mut bodies = far;
    assert!(matches!(
        apply_rigid_hertz_impulse(&mut bodies, &[p], c, 1.0, controls(), || false),
        Err(RigidHertzImpulseError::Invalid(_))
    ));
    assert_eq!(bodies, far);
}
