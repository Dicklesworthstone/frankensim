//! G1/G3/G4: restitution, finite-patch friction and original-body recoil.

use fs_contact::tangential::impact::hertz::HertzImpactControls;
use fs_ga::Mat3;
use fs_mbd::frictional_impact::RigidImpactPoint;
use fs_mbd::frictional_impact::hertz_patch::impact::{
    RigidHertzImpactControls, RigidHertzImpactError, RigidHertzImpactReport,
    apply_rigid_hertz_impact,
};
use fs_mbd::frictional_impact::hertz_patch::{
    PreparedRigidHertzImpulse, RigidHertzImpulseContact, RigidHertzImpulseControls,
    RigidHertzImpulseCounterpart,
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_solver::orthant::active_set::SpdOrthantConfig;
use fs_tribo::sliding_patch::contact::impulse::{HertzImpulseControls, HertzImpulsePatch};

const X: Vec3 = Vec3::new(1.0, 0.0, 0.0);
const Y: Vec3 = Vec3::new(0.0, 1.0, 0.0);
const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);

fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; tolerance {tolerance:e}"
    );
}
fn vector(a: Vec3, b: Vec3, tolerance: f64) {
    near(a.x, b.x, tolerance);
    near(a.y, b.y, tolerance);
    near(a.z, b.z, tolerance);
}
fn property(mass: f64, inertia: f64) -> MassProperties {
    MassProperties::new(mass, Vec3::ZERO, Vec3::new(inertia, inertia, inertia)).unwrap()
}
fn state(p: MassProperties, position: Vec3, velocity: Vec3, omega: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(position, UnitQuaternion::IDENTITY).unwrap(),
        velocity.scale(p.mass()),
        p.angular_momentum_body_checked(omega).unwrap(),
    )
    .unwrap()
}
fn controls() -> RigidHertzImpactControls {
    let patch = HertzImpulseControls {
        maximum_static_iterations: 128,
        maximum_iterations: 256,
        maximum_evaluations: 4096,
        impulse_tolerance_ns: 2e-9,
        spin_impulse_tolerance_nms: 5e-10,
        velocity_tolerance_m_s: 2e-13,
        energy_tolerance_j: 2e-10,
        relative_matrix_tolerance: 1e-13,
        maximum_impulse_ns: 30.0,
        maximum_spin_impulse_nms: 10.0,
    };
    RigidHertzImpactControls {
        native: RigidHertzImpulseControls {
            maximum_bodies: 8,
            geometry_tolerance_m: 1e-10,
            velocity_tolerance_m_s: 5e-9,
            spin_tolerance_rad_s: 5e-9,
            linear_momentum_tolerance_ns: 3e-9,
            angular_momentum_tolerance_nms: 3e-9,
            energy_tolerance_j: 3e-8,
            contact: patch,
        },
        impact: HertzImpactControls {
            patch,
            matrix: SpdOrthantConfig {
                max_variables: 4,
                ..SpdOrthantConfig::default()
            },
            maximum_iterations: 96,
            maximum_evaluations: 200_000,
            normal_impulse_tolerance_ns: 1e-10,
            normal_velocity_tolerance_m_s: 2e-9,
            energy_tolerance_j: 2e-8,
            maximum_normal_impulse_ns: 40.0,
            minimum_determinacy_margin: 1e-8,
        },
    }
}
fn support<'a>(
    law: &'a HertzImpulsePatch,
    arm: Vec3,
    velocity: Vec3,
    spin: Vec3,
) -> RigidHertzImpulseContact<'a> {
    RigidHertzImpulseContact {
        point: RigidImpactPoint {
            body: 0,
            arm_body: arm,
        },
        other: RigidHertzImpulseCounterpart::PrescribedPatch {
            point_world: Vec3::ZERO,
            velocity_world: velocity,
            angular_velocity_world: spin,
        },
        normal_world: Z,
        tangent_world: X,
        law,
    }
}
fn audit(report: &RigidHertzImpactReport) {
    let m = &report.mechanics;
    near(m.work_energy_residual_j, 0.0, 3e-8);
    near(m.total_energy_residual_j, 0.0, 3e-8);
    near(m.contact_work_residual_j, 0.0, 3e-8);
    near(report.normal_complementarity_j, 0.0, 2e-8);
    near(m.normal_response_residual_m_s, 0.0, 5e-9);
    near(m.endpoint_tangent_residual_m_s, 0.0, 5e-9);
    near(m.endpoint_spin_residual_rad_s, 0.0, 5e-9);
    vector(m.linear_momentum_residual_ns, Vec3::ZERO, 3e-9);
    vector(m.angular_momentum_residual_nms, Vec3::ZERO, 3e-9);
    assert!(m.dissipated_energy_j >= -controls().impact.energy_tolerance_j);
    near(m.dissipated_energy_j, report.impact.dissipated_j, 3e-8);
    near(
        m.normal_impulse_applied_ns,
        report.impact.scaled_impulse_ns[0],
        0.0,
    );
}

#[test]
fn g1_two_finite_masses_share_one_normal_impulse_and_newton_restitution() {
    let props = [property(2.0, 1.0), property(3.0, 2.0)];
    let initial = [
        state(props[0], Vec3::ZERO, Z.scale(-3.0), Vec3::ZERO),
        state(props[1], Vec3::ZERO, Z.scale(2.0), Vec3::ZERO),
    ];
    let law = HertzImpulsePatch::new(0.25, 0.0).unwrap();
    let mut contact = support(&law, Vec3::ZERO, Vec3::ZERO, Vec3::ZERO);
    contact.other = RigidHertzImpulseCounterpart::Body(RigidImpactPoint {
        body: 1,
        arm_body: Vec3::ZERO,
    });
    for e in [0.0, 0.5, 1.0] {
        let mut bodies = initial;
        let report =
            apply_rigid_hertz_impact(&mut bodies, &props, contact, e, controls(), || false)
                .unwrap();
        let normal = 6.0 * (1.0 + e);
        near(report.impact.scaled_impulse_ns[0], normal, 2e-9);
        near(report.mechanics.normal_velocity_after_m_s, 5.0 * e, 2e-9);
        vector(
            bodies[0].linear_momentum_world(),
            Z.scale(-6.0 + normal),
            2e-9,
        );
        vector(
            bodies[1].linear_momentum_world(),
            Z.scale(6.0 - normal),
            2e-9,
        );
        near(report.impact.dissipated_j, 15.0 * (1.0 - e * e), 2e-8);
        vector(
            report.mechanics.force_impulse_world_ns,
            Z.scale(normal),
            2e-9,
        );
        vector(report.mechanics.spin_impulse_world_nms, Vec3::ZERO, 0.0);
        assert_eq!(bodies[0].pose(), initial[0].pose());
        assert_eq!(bodies[1].pose(), initial[1].pose());
        audit(&report);
    }
}

#[test]
fn g1_off_center_sliding_changes_the_normal_impulse_through_the_same_body() {
    let props = [property(1.0, 1.0)];
    let arm = Vec3::new(1.0, 0.0, 1.0);
    let mut bodies = [state(
        props[0],
        arm.scale(-1.0),
        Vec3::new(4.0, 0.0, -2.0),
        Vec3::ZERO,
    )];
    let law = HertzImpulsePatch::new(0.2, 0.3).unwrap();
    let contact = support(&law, arm, Vec3::ZERO, Vec3::ZERO);
    let report =
        apply_rigid_hertz_impact(&mut bodies, &props, contact, 0.25, controls(), || false).unwrap();
    // Independent point mechanics: Wnn=Wxx=2, Wnx=-1. Sliding gives
    // px=-.3*Jn, so .5=-2+(2+.3)*Jn. A normal-first solve would give 1.25.
    let normal = 25.0 / 23.0;
    let tangent = -0.3 * normal;
    near(report.impact.scaled_impulse_ns[0], normal, 2e-9);
    near(report.impact.scaled_impulse_ns[1], tangent, 2e-9);
    near(report.mechanics.normal_velocity_after_m_s, 0.5, 2e-9);
    near(
        report.mechanics.slip_velocity_after_m_s[0],
        4.0 - 1.6 * normal,
        3e-9,
    );
    vector(
        bodies[0].linear_momentum_world(),
        Vec3::new(4.0 + tangent, 0.0, -2.0 + normal),
        3e-9,
    );
    vector(
        bodies[0].angular_momentum_body(),
        Y.scale(tangent - normal),
        3e-9,
    );
    near(report.mechanics.spin_after_rad_s, 0.0, 1e-12);
    audit(&report);
}

fn inverse(a: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let cofactors = core::array::from_fn::<_, 3, _>(|i| {
        core::array::from_fn::<_, 3, _>(|j| {
            a[(i + 1) % 3][(j + 1) % 3] * a[(i + 2) % 3][(j + 2) % 3]
                - a[(i + 1) % 3][(j + 2) % 3] * a[(i + 2) % 3][(j + 1) % 3]
        })
    });
    let determinant = (0..3).map(|j| a[0][j] * cofactors[0][j]).sum::<f64>();
    core::array::from_fn(|i| core::array::from_fn(|j| cofactors[j][i] / determinant))
}
fn multiply(a: [[f64; 3]; 3], v: Vec3) -> Vec3 {
    let v = [v.x, v.y, v.z];
    let r = a.map(|row| (0..3).map(|j| row[j] * v[j]).sum::<f64>());
    Vec3::new(r[0], r[1], r[2])
}

#[test]
fn g3_full_inertia_mixed_impact_retains_support_work_and_world_covariance() {
    let inertia = [[2.4, 0.2, 0.15], [0.2, 1.8, -0.1], [0.15, -0.1, 2.0]];
    let props = [MassProperties::from_inertia_tensor(
        2.0,
        Vec3::ZERO,
        Mat3 {
            m: [2.4, 0.2, 0.15, 0.2, 1.8, -0.1, 0.15, -0.1, 2.0],
        },
    )
    .unwrap()];
    let arm = Vec3::new(0.3, -0.2, 0.6);
    let initial = state(
        props[0],
        arm.scale(-1.0),
        Vec3::new(1.6, -0.4, -1.5),
        Vec3::new(0.8, -0.3, 2.1),
    );
    let support_v = Vec3::new(-0.2, 0.15, 0.25);
    let support_w = Vec3::new(-0.1, 0.2, -0.4);
    let law = HertzImpulsePatch::new(0.35, 0.12).unwrap();
    let contact = support(&law, arm, support_v, support_w);
    let prepared =
        PreparedRigidHertzImpulse::new(&[initial], &props, contact, controls().native, || false)
            .unwrap();
    let response = prepared.impact_mobility_per_kg();
    let forces = [Z, X, Y, Vec3::ZERO];
    let moments = [
        arm.cross(Z),
        arm.cross(X),
        arm.cross(Y),
        Z.scale(law.radius_m()),
    ];
    let inv = inverse(inertia);
    for i in 0..4 {
        for j in 0..4 {
            let expected =
                forces[i].dot(forces[j]) / 2.0 + moments[i].dot(multiply(inv, moments[j]));
            near(response[i][j], expected, 4e-15);
        }
    }
    assert!(response[0][3].abs() > 1e-4);
    let mut bodies = [initial];
    let report =
        apply_rigid_hertz_impact(&mut bodies, &props, contact, 0.2, controls(), || false).unwrap();
    let m = &report.mechanics;
    assert!(m.force_impulse_world_ns.x.abs() > 1e-4);
    assert!(m.spin_impulse_world_nms.z.abs() > 1e-4);
    let expected_p = initial
        .linear_momentum_world()
        .add(m.force_impulse_world_ns);
    let expected_l = initial
        .angular_momentum_body()
        .add(arm.cross(m.force_impulse_world_ns))
        .add(m.spin_impulse_world_nms);
    vector(bodies[0].linear_momentum_world(), expected_p, 2e-14);
    vector(bodies[0].angular_momentum_body(), expected_l, 2e-14);
    near(
        m.support_work_j,
        m.force_impulse_world_ns.dot(support_v) + m.spin_impulse_world_nms.dot(support_w),
        2e-13,
    );
    let kinetic = |p: Vec3, l: Vec3| 0.25 * p.dot(p) + 0.5 * l.dot(multiply(inv, l));
    near(
        m.kinetic_change_j,
        kinetic(expected_p, expected_l)
            - kinetic(
                initial.linear_momentum_world(),
                initial.angular_momentum_body(),
            ),
        3e-13,
    );
    audit(&report);

    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(2.0, -1.0, 3.0), 0.7).unwrap();
    let origin = Vec3::new(7.0, -3.0, 5.0);
    let boost = Vec3::new(1.0, -2.0, 3.0);
    let mut transformed = [RigidBodyState::new(
        Pose::new(
            rotation
                .rotate_body_to_world(initial.pose().position_world())
                .add(origin),
            rotation,
        )
        .unwrap(),
        rotation
            .rotate_body_to_world(initial.linear_momentum_world())
            .add(boost.scale(2.0)),
        initial.angular_momentum_body(),
    )
    .unwrap()];
    let mut moved_contact = contact;
    moved_contact.normal_world = rotation.rotate_body_to_world(Z);
    moved_contact.tangent_world = rotation.rotate_body_to_world(X);
    moved_contact.other = RigidHertzImpulseCounterpart::PrescribedPatch {
        point_world: origin,
        velocity_world: rotation.rotate_body_to_world(support_v).add(boost),
        angular_velocity_world: rotation.rotate_body_to_world(support_w),
    };
    let moved = apply_rigid_hertz_impact(
        &mut transformed,
        &props,
        moved_contact,
        0.2,
        controls(),
        || false,
    )
    .unwrap();
    vector(
        transformed[0].linear_momentum_world(),
        rotation
            .rotate_body_to_world(bodies[0].linear_momentum_world())
            .add(boost.scale(2.0)),
        2e-8,
    );
    vector(
        transformed[0].angular_momentum_body(),
        bodies[0].angular_momentum_body(),
        2e-8,
    );
    near(moved.impact.dissipated_j, report.impact.dissipated_j, 2e-8);
    vector(
        moved.mechanics.force_impulse_world_ns,
        rotation.rotate_body_to_world(m.force_impulse_world_ns),
        2e-8,
    );
    audit(&moved);
}

#[test]
fn g3_separating_contact_and_prepared_trials_preserve_the_original_state() {
    let props = [property(1.0, 1.0)];
    let initial = [state(
        props[0],
        Vec3::ZERO,
        Vec3::new(2.0, 1.0, 3.0),
        Z.scale(4.0),
    )];
    let law = HertzImpulsePatch::new(0.2, 0.4).unwrap();
    let contact = support(&law, Vec3::ZERO, Vec3::ZERO, Vec3::ZERO);
    let prepared =
        PreparedRigidHertzImpulse::new(&initial, &props, contact, controls().native, || false)
            .unwrap();
    let first = prepared
        .solve_impact(&law, 0.5, controls().impact, || false)
        .unwrap();
    assert_eq!(first.impact.scaled_impulse_ns, [0.0; 4]);
    assert_eq!(first.mechanics.point_impulse.state_after, initial[0]);
    assert_eq!(
        first,
        prepared
            .solve_impact(&law, 0.5, controls().impact, || false)
            .unwrap()
    );
    let changed_radius = HertzImpulsePatch::new(0.3, 0.4).unwrap();
    assert!(
        prepared
            .solve_impact(&changed_radius, 0.5, controls().impact, || false)
            .is_err()
    );
    let mut bodies = initial;
    apply_rigid_hertz_impact(&mut bodies, &props, contact, 0.5, controls(), || false).unwrap();
    assert_eq!(bodies, initial);
    audit(&first);
}

#[test]
fn g4_every_cancellation_and_exhausted_normal_budget_publish_no_partial_momentum() {
    let props = [property(2.0, 1.0)];
    let initial = [state(props[0], Vec3::ZERO, Z.scale(-2.0), Vec3::ZERO)];
    let law = HertzImpulsePatch::new(0.25, 0.0).unwrap();
    let contact = support(&law, Vec3::ZERO, Vec3::ZERO, Vec3::ZERO);
    let mut total_polls = 0;
    let mut bodies = initial;
    apply_rigid_hertz_impact(&mut bodies, &props, contact, 0.5, controls(), || {
        total_polls += 1;
        false
    })
    .unwrap();
    assert!(total_polls > 4);
    for stop in 1..=total_polls {
        let mut bodies = initial;
        let mut polls = 0;
        let result =
            apply_rigid_hertz_impact(&mut bodies, &props, contact, 0.5, controls(), || {
                polls += 1;
                polls >= stop
            });
        assert!(
            matches!(result, Err(RigidHertzImpactError::Cancelled)),
            "poll {stop}: {result:?}"
        );
        assert_eq!(bodies, initial);
    }
    let mut too_small = controls();
    too_small.impact.maximum_normal_impulse_ns = 5.0;
    let mut bodies = initial;
    assert!(
        apply_rigid_hertz_impact(&mut bodies, &props, contact, 0.5, too_small, || false).is_err()
    );
    assert_eq!(bodies, initial);
}
