//! Curved sourced adhesive traction recoils both original finite bodies.
#[path = "../../fs-contact/tests/adhesive_dmt_coupled/common.rs"]
mod common;
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::{
    RigidContactBinding,
    adhesive::{transient::*, *},
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_tribo::adhesive::spinning::contact::AdhesiveWrenchConfig;
use fs_tribo::partial_slip::TangentFrame;

fn controls() -> RigidAdhesiveControls {
    RigidAdhesiveControls {
        maximum_bodies: 4,
        geometry_tolerance_m: 1e-18,
        velocity_tolerance_m_s: 1e-12,
        angular_velocity_tolerance_rad_s: 1e-4,
        linear_momentum_tolerance_ns: 1e-30,
        angular_momentum_tolerance_nms: 1e-36,
        energy_tolerance_j: 2e-31,
        contact: AdhesiveWrenchConfig {
            maximum_static_iterations: 64,
            maximum_iterations: 128,
            maximum_evaluations: 256,
            force_tolerance_n: 1e-19,
            moment_tolerance_nm: 1e-27,
            displacement_tolerance_m: 1e-27,
            energy_tolerance_j: 1e-32,
            relative_matrix_tolerance: 1e-13,
            maximum_force_n: 1e-6,
            maximum_spin_moment_nm: 1e-12,
        },
    }
}
fn transient() -> AdhesiveTransientControls {
    AdhesiveTransientControls {
        maximum_step_s: 1e-12,
        minimum_step_s: 1e-25,
        maximum_attempts: 100_000,
        maximum_evaluations: 400_000,
        scaled_rate_tolerance_m_s: 1e-12,
        scaled_impulse_tolerance_ns: 1e-28,
        scaled_travel_tolerance_m: 1e-23,
        heat_tolerance_j: 1e-33,
        stopping_time_tolerance_s: 1e-17,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    common::close(a, b, tolerance, 2e-9);
}
fn close_vec(a: Vec3, b: Vec3, tolerance: f64) {
    close(a.x, b.x, tolerance);
    close(a.y, b.y, tolerance);
    close(a.z, b.z, tolerance);
}
struct Fixture {
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    binding: RigidContactBinding,
    frame: TangentFrame,
}
fn fixture(radius: f64) -> Fixture {
    let arms = [
        Vec3::new(0.2, -0.3, 0.0).scale(radius),
        Vec3::new(-0.15, 0.1, 0.0).scale(radius),
    ];
    let properties = [2e-15, 3e-15].map(|m| {
        MassProperties::new(
            m,
            Vec3::ZERO,
            Vec3::new(0.7, 0.9, 1.1).scale(m * radius * radius),
        )
        .unwrap()
    });
    let poses = arms.map(|a| Pose::new(a.scale(-1.0), UnitQuaternion::IDENTITY).unwrap());
    let omega = Vec3::new(0.0, 0.0, 1000.0);
    let com = Vec3::new(0.8, -0.3, 0.0)
        .scale(radius * 1000.0)
        .sub(omega.cross(arms[0]));
    let bodies = [
        RigidBodyState::new(
            poses[0],
            com.scale(properties[0].mass()),
            Vec3::new(0.0, 0.0, 1000.0 * properties[0].principal_inertia_body().z),
        )
        .unwrap(),
        RigidBodyState::new(poses[1], Vec3::ZERO, Vec3::ZERO).unwrap(),
    ];
    Fixture {
        bodies,
        properties,
        binding: RigidContactBinding {
            body_a: 0,
            body_b: 1,
            arm_a_body: arms[0],
            arm_b_body: arms[1],
        },
        frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
    }
}
fn kinetic(bodies: [RigidBodyState; 2], properties: [MassProperties; 2]) -> f64 {
    (0..2)
        .map(|i| {
            let p = bodies[i].linear_momentum_world();
            let l = bodies[i].angular_momentum_body();
            let inertia = properties[i].principal_inertia_body();
            0.5 * (p.dot(p) / properties[i].mass()
                + l.x * l.x / inertia.x
                + l.y * l.y / inertia.y
                + l.z * l.z / inertia.z)
        })
        .sum()
}

#[test]
fn curved_motion_updates_both_bodies_with_actual_force_spin_work_and_gauge_covariance() {
    let (normal, shear) = common::sources(0.05, 3e6, 2e6);
    let response = normal.evaluate_contact_load(0.0).unwrap();
    let source = shear
        .sliding_patch_dmt(normal.normal_state(), &response)
        .unwrap();
    let f = fixture(response.contact_radius_m);
    let mut original: Option<[RigidBodyState; 2]> = None;
    for frame in [
        f.frame,
        TangentFrame::new([0.0, 0.0, 1.0], [0.0, 1.0, 0.0]).unwrap(),
    ] {
        let contact = RigidAdhesiveContact {
            binding: f.binding,
            frame,
            source,
        };
        let mut bodies = f.bodies;
        let report = apply_rigid_adhesive_transient(
            &mut bodies,
            &f.properties,
            contact,
            1e-12,
            controls(),
            transient(),
            &CancelGate::new(),
        )
        .unwrap();
        let path = report.transient.as_ref().unwrap();
        assert!(!path.exact_ray && path.stopping_time_bracket_s.is_none());
        assert!(path.accepted_steps > 1 && path.dissipated_energy_j > 0.0);
        let w = report.scaled_mobility_per_kg.unwrap();
        assert!(w[0][1] != 0.0 && w[0][2] != 0.0 && w[1][2] != 0.0);
        close(
            kinetic(bodies, f.properties) - kinetic(f.bodies, f.properties),
            -report.dissipated_energy_j,
            2e-31,
        );
        close(
            report.kinetic_energy_change_j,
            -report.dissipated_energy_j,
            2e-31,
        );
        close(report.impulse_work_j, -report.dissipated_energy_j, 2e-31);
        let impulse = report.force_impulse_on_a_world_ns;
        let couple = report.spin_angular_impulse_on_a_world_nms;
        for (i, sign) in [(0, 1.0), (1, -1.0)] {
            let arm = [f.binding.arm_a_body, f.binding.arm_b_body][i];
            close_vec(
                bodies[i]
                    .linear_momentum_world()
                    .sub(f.bodies[i].linear_momentum_world()),
                impulse.scale(sign),
                1e-32,
            );
            close_vec(
                bodies[i]
                    .angular_momentum_body()
                    .sub(f.bodies[i].angular_momentum_body()),
                arm.cross(impulse).add(couple).scale(sign),
                1e-40,
            );
            assert_eq!(bodies[i].pose(), f.bodies[i].pose());
        }
        assert_ne!(bodies[1].linear_momentum_world(), Vec3::ZERO);
        assert_eq!(report.audit.unwrap().normal_velocity_bound_m_s, 0.0);
        if let Some(reference) = original {
            for i in 0..2 {
                close_vec(
                    bodies[i].linear_momentum_world(),
                    reference[i].linear_momentum_world(),
                    2e-28,
                );
                close_vec(
                    bodies[i].angular_momentum_body(),
                    reference[i].angular_momentum_body(),
                    2e-36,
                );
            }
        }
        original = Some(bodies);
    }
}

#[test]
fn full_curved_arrest_retains_native_common_motion_and_only_relative_heat() {
    let (normal, shear) = common::sources(0.05, 3e6, 2e6);
    let response = normal.evaluate_contact_load(0.0).unwrap();
    let source = shear
        .sliding_patch_dmt(normal.normal_state(), &response)
        .unwrap();
    let f = fixture(response.contact_radius_m);
    let contact = RigidAdhesiveContact {
        binding: f.binding,
        frame: f.frame,
        source,
    };
    let mut bodies = f.bodies;
    let report = apply_rigid_adhesive_transient(
        &mut bodies,
        &f.properties,
        contact,
        1e-8,
        controls(),
        transient(),
        &CancelGate::new(),
    )
    .unwrap();
    let path = report.transient.as_ref().unwrap();
    let [lo, hi] = path.stopping_time_bracket_s.unwrap();
    assert!(!path.exact_ray && lo < hi && hi < 1e-8);
    assert!(path.stopping_tail_heat_j > 0.0);
    close(
        report.slip_velocity_m_s[1][0],
        0.0,
        controls().velocity_tolerance_m_s,
    );
    close(
        report.slip_velocity_m_s[1][1],
        0.0,
        controls().velocity_tolerance_m_s,
    );
    close(
        report.spin_rad_s[1],
        0.0,
        controls().angular_velocity_tolerance_rad_s,
    );
    let initial = kinetic(f.bodies, f.properties);
    let final_energy = kinetic(bodies, f.properties);
    assert!(final_energy > 0.0 && final_energy < initial);
    close(
        initial - final_energy,
        report.dissipated_energy_j,
        controls().energy_tolerance_j,
    );
    close_vec(
        bodies[0]
            .linear_momentum_world()
            .add(bodies[1].linear_momentum_world()),
        f.bodies[0]
            .linear_momentum_world()
            .add(f.bodies[1].linear_momentum_world()),
        1e-30,
    );
    assert!(path.prefix_work_residual_j.abs() <= controls().contact.energy_tolerance_j);
}

#[test]
fn normal_path_cancellation_budget_and_detached_cases_never_publish_partial_bodies() {
    let (normal, shear) = common::sources(0.05, 3e6, 2e6);
    let response = normal.evaluate_contact_load(0.0).unwrap();
    let source = shear
        .sliding_patch_dmt(normal.normal_state(), &response)
        .unwrap();
    let f = fixture(response.contact_radius_m);
    let contact = RigidAdhesiveContact {
        binding: f.binding,
        frame: f.frame,
        source,
    };
    let mut bodies = f.bodies;
    let gate = CancelGate::new();
    gate.request();
    assert!(
        apply_rigid_adhesive_transient(
            &mut bodies,
            &f.properties,
            contact,
            1e-12,
            controls(),
            transient(),
            &gate
        )
        .is_err()
    );
    assert_eq!(bodies, f.bodies);
    assert!(
        apply_rigid_adhesive_transient(
            &mut bodies,
            &f.properties,
            contact,
            1e-12,
            controls(),
            AdhesiveTransientControls {
                maximum_evaluations: 2,
                ..transient()
            },
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(bodies, f.bodies);
    let mut coupled = contact;
    coupled.binding.arm_a_body =
        coupled
            .binding
            .arm_a_body
            .add(Vec3::new(0.0, 0.0, 0.15 * response.contact_radius_m));
    let pose = Pose::new(
        coupled.binding.arm_a_body.scale(-1.0),
        UnitQuaternion::IDENTITY,
    )
    .unwrap();
    bodies[0] = RigidBodyState::new(
        pose,
        bodies[0].linear_momentum_world(),
        bodies[0].angular_momentum_body(),
    )
    .unwrap();
    let original = bodies;
    assert!(matches!(
        apply_rigid_adhesive_transient(
            &mut bodies,
            &f.properties,
            coupled,
            1e-12,
            controls(),
            transient(),
            &CancelGate::new()
        ),
        Err(RigidAdhesiveTransientError::Invalid(_))
    ));
    assert_eq!(bodies, original);
    let first = prepare_rigid_adhesive_transient(
        &f.bodies,
        &f.properties,
        contact,
        1e-12,
        controls(),
        transient(),
        &CancelGate::new(),
    )
    .unwrap();
    let retry = prepare_rigid_adhesive_transient(
        &f.bodies,
        &f.properties,
        contact,
        1e-12,
        controls(),
        transient(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(first, retry);
    let detached = normal.evaluate(-0.5 * common::Z).unwrap();
    let source = shear
        .sliding_patch_dmt(normal.normal_state(), &detached)
        .unwrap();
    let report = apply_rigid_adhesive_transient(
        &mut bodies,
        &f.properties,
        RigidAdhesiveContact { source, ..coupled },
        1.0,
        controls(),
        transient(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(bodies, original);
    assert!(report.transient.is_none() && report.audit.is_none());
    assert_eq!(report.dissipated_energy_j, 0.0);
}
