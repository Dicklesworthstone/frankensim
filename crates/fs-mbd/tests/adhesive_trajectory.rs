//! Actual finite bodies receive only the impulse up to a located adhesive stop.
#[path = "../../fs-contact/tests/adhesive_dmt_coupled/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::{
    RigidContactBinding,
    adhesive::{trajectory::*, *},
};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_tribo::adhesive::spinning::contact::AdhesiveWrenchConfig;
use fs_tribo::adhesive::spinning::relaxation::AdhesiveRelaxationError;
use fs_tribo::partial_slip::TangentFrame;

const STRENGTH: f64 = 2e6;

fn close(actual: f64, expected: f64, absolute: f64) {
    common::close(actual, expected, absolute, 3e-10);
}
fn close_vec(actual: Vec3, expected: Vec3, absolute: f64) {
    close(actual.x, expected.x, absolute);
    close(actual.y, expected.y, absolute);
    close(actual.z, expected.z, absolute);
}
fn controls() -> RigidAdhesiveControls {
    RigidAdhesiveControls {
        maximum_bodies: 4,
        geometry_tolerance_m: 1e-18,
        velocity_tolerance_m_s: 1e-13,
        angular_velocity_tolerance_rad_s: 1e-7,
        linear_momentum_tolerance_ns: 1e-30,
        angular_momentum_tolerance_nms: 1e-34,
        energy_tolerance_j: 1e-31,
        contact: AdhesiveWrenchConfig {
            maximum_static_iterations: 96,
            maximum_iterations: 256,
            maximum_evaluations: 512,
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
fn state(pose: Pose, property: MassProperties, velocity: Vec3, spin: f64) -> RigidBodyState {
    RigidBodyState::new(
        pose,
        velocity.scale(property.mass()),
        Vec3::new(0.0, 0.0, spin * property.principal_inertia_body().z),
    )
    .unwrap()
}
struct Fixture {
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    binding: RigidContactBinding,
    frame: TangentFrame,
}
fn fixture(radius: f64, coupled: bool) -> Fixture {
    let arm_fraction = if coupled { 0.2 } else { 0.0 };
    let arm = Vec3::new(0.0, -arm_fraction * radius, 0.0);
    let pose = Pose::new(arm.scale(-1.0), UnitQuaternion::IDENTITY).unwrap();
    let properties = [2e-15, 3e-15].map(|mass| {
        // At the rim F=3 M/a. With arm=b*a and I=m*a²*(1+2b-3b²)/3,
        // the original complete mobility maps [F,0,M/a] onto [1,0,1].
        let inertia = if coupled {
            mass * radius * radius * (1.0 + 2.0 * arm_fraction - 3.0 * arm_fraction * arm_fraction)
                / 3.0
        } else {
            mass * 1e-12
        };
        MassProperties::new(mass, Vec3::ZERO, Vec3::new(inertia, inertia, inertia)).unwrap()
    });
    let spin = 1000.0;
    let velocity = if coupled {
        Vec3::new((1.0 - arm_fraction) * radius * spin, 0.0, 0.0)
    } else {
        Vec3::ZERO
    };
    Fixture {
        bodies: [
            state(pose, properties[0], velocity, spin),
            state(pose, properties[1], Vec3::ZERO, 0.0),
        ],
        properties,
        binding: RigidContactBinding {
            body_a: 0,
            body_b: 1,
            arm_a_body: arm,
            arm_b_body: arm,
        },
        frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
    }
}

#[test]
fn zero_and_tensile_load_spin_stops_then_dwells_without_extra_impulse_or_heat() {
    let (normal, shear) = common::sources(0.05, 3e6, STRENGTH);
    let pull_off =
        2.0 * PI * normal.effective_radius_m() * normal.normal_state().adhesion_energy_j_per_m2();
    for load in [0.0, -0.5 * pull_off] {
        let response = normal.evaluate_contact_load(load).unwrap();
        let source = shear
            .sliding_patch_dmt(normal.normal_state(), &response)
            .unwrap();
        let f = fixture(response.contact_radius_m, false);
        let ia = f.properties[0].principal_inertia_body().z;
        let ib = f.properties[1].principal_inertia_body().z;
        let inverse = 1.0 / ia + 1.0 / ib;
        let moment = (2.0 / 3.0) * STRENGTH * response.contact_area_m2 * response.contact_radius_m;
        let stop = 1000.0 / (moment * inverse);
        let contact = RigidAdhesiveContact {
            binding: f.binding,
            frame: f.frame,
            source,
        };
        let mut old_attempt = f.bodies;
        assert!(matches!(
            apply_rigid_adhesive(
                &mut old_attempt,
                &f.properties,
                contact,
                1.5 * stop,
                controls(),
                &CancelGate::new()
            ),
            Err(RigidAdhesiveError::UnresolvedSlipEvent { .. })
        ));
        assert_eq!(old_attempt, f.bodies);
        let mut accepted: Option<[RigidBodyState; 2]> = None;
        for duration in [1.5 * stop, 5.0 * stop] {
            let mut bodies = f.bodies;
            let report = apply_rigid_adhesive_trajectory(
                &mut bodies,
                &f.properties,
                contact,
                duration,
                controls(),
                &CancelGate::new(),
            )
            .unwrap();
            let path = report.relaxation.unwrap();
            close(path.kinetic_phase_end_s.unwrap(), stop, 1e-24);
            close(path.stationary_duration_s, duration - stop, 1e-24);
            close(path.spin_increment_rad, 0.5 * 1000.0 * stop, 1e-21);
            close(report.spin_after_rad_s, 0.0, 1e-9);
            close(
                bodies[0].angular_momentum_body().z,
                ia * 1000.0 - moment * stop,
                1e-35,
            );
            close(bodies[1].angular_momentum_body().z, moment * stop, 1e-35);
            let heat = 0.5 * 1000.0_f64.powi(2) / inverse;
            close(report.dissipated_energy_j, heat, 1e-33);
            close(report.kinetic_energy_change_j, -heat, 1e-33);
            close(report.impulse_work_j, -heat, 1e-33);
            assert!(report.audit.total_energy_residual_j.abs() < 1e-33);
            assert_eq!(report.force_impulse_on_a_world_ns, Vec3::ZERO);
            assert!(bodies[1].angular_momentum_body().z > 0.0);
            for i in 0..2 {
                assert_eq!(bodies[i].pose(), f.bodies[i].pose());
            }
            if let Some(previous) = accepted {
                assert_eq!(bodies, previous);
            }
            accepted = Some(bodies);
        }
    }
}

#[test]
fn coupled_rim_ray_uses_full_native_recoil_and_is_tangent_gauge_covariant() {
    let (normal, shear) = common::sources(0.05, 3e6, STRENGTH);
    let response = normal.evaluate_contact_load(0.0).unwrap();
    let radius = response.contact_radius_m;
    let source = shear
        .sliding_patch_dmt(normal.normal_state(), &response)
        .unwrap();
    let f = fixture(radius, true);
    let capacity = STRENGTH * response.contact_area_m2;
    let force = capacity * 8.0 / (3.0 * PI);
    let moment = capacity * radius * 8.0 / (9.0 * PI);
    let inverse_mass = 1.0 / f.properties[0].mass() + 1.0 / f.properties[1].mass();
    let inverse_inertia = 1.0 / f.properties[0].principal_inertia_body().z
        + 1.0 / f.properties[1].principal_inertia_body().z;
    let b = 0.2 * radius;
    let deceleration =
        (inverse_mass + b * b * inverse_inertia) * force + b * inverse_inertia * moment;
    let stop = radius * 1000.0 / deceleration;
    let expected_impulse = Vec3::new(-force * stop, 0.0, 0.0);
    let expected_couple = Vec3::new(0.0, 0.0, -moment * stop);
    let heat = 0.5 * (force * radius * 1000.0 + moment * 1000.0) * stop;
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
        let candidate = prepare_rigid_adhesive_trajectory(
            &f.bodies,
            &f.properties,
            contact,
            3.0 * stop,
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(f.bodies[1].linear_momentum_world(), Vec3::ZERO);
        let report = &candidate.report;
        let path = report.relaxation.unwrap();
        close(path.kinetic_phase_end_s.unwrap(), stop, 1e-24);
        close_vec(report.force_impulse_on_a_world_ns, expected_impulse, 1e-33);
        close_vec(
            report.spin_angular_impulse_on_a_world_nms,
            expected_couple,
            1e-40,
        );
        close(report.dissipated_energy_j, heat, 1e-34);
        close(report.kinetic_energy_change_j, -heat, 1e-34);
        close(report.spin_after_rad_s, 0.0, 1e-8);
        for (i, sign) in [(0, 1.0), (1, -1.0)] {
            let linear = f.bodies[i]
                .linear_momentum_world()
                .add(expected_impulse.scale(sign));
            let angular = f.bodies[i].angular_momentum_body().add(
                f.binding
                    .arm_a_body
                    .cross(expected_impulse)
                    .add(expected_couple)
                    .scale(sign),
            );
            close_vec(
                candidate.bodies_after[i].linear_momentum_world(),
                linear,
                1e-33,
            );
            close_vec(
                candidate.bodies_after[i].angular_momentum_body(),
                angular,
                1e-40,
            );
        }
        let w = report.scaled_mobility_per_kg.unwrap();
        assert!(w[0][2] != 0.0 || w[1][2] != 0.0);
        if let Some(previous) = original {
            for i in 0..2 {
                close_vec(
                    candidate.bodies_after[i].linear_momentum_world(),
                    previous[i].linear_momentum_world(),
                    1e-33,
                );
                close_vec(
                    candidate.bodies_after[i].angular_momentum_body(),
                    previous[i].angular_momentum_body(),
                    1e-40,
                );
            }
        }
        original = Some(candidate.bodies_after);
    }
}

#[test]
fn detached_attractive_gap_is_an_exact_noop_for_the_complete_requested_time() {
    let (normal, shear) = common::sources(0.05, 3e6, STRENGTH);
    let response = normal.evaluate(-0.5 * common::Z).unwrap();
    assert!(response.normal_force_n < 0.0 && response.contact_area_m2 == 0.0);
    let source = shear
        .sliding_patch_dmt(normal.normal_state(), &response)
        .unwrap();
    let f = fixture(1e-8, false);
    let mut separated = f.binding;
    separated.arm_a_body = Vec3::new(1e-6, 0.0, 0.0);
    let mut bodies = f.bodies;
    let report = apply_rigid_adhesive_trajectory(
        &mut bodies,
        &f.properties,
        RigidAdhesiveContact {
            binding: separated,
            frame: f.frame,
            source,
        },
        1.0,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(bodies, f.bodies);
    assert!(report.relaxation.is_none());
    assert!(report.scaled_mobility_per_kg.is_none());
    assert_eq!(report.force_impulse_on_a_world_ns, Vec3::ZERO);
    assert_eq!(report.spin_angular_impulse_on_a_world_nms, Vec3::ZERO);
    assert_eq!(report.dissipated_energy_j, 0.0);
    assert_eq!(report.kinetic_energy_change_j, 0.0);
}

#[test]
fn curved_motion_cancellation_and_late_native_refusal_are_atomic() {
    let (normal, shear) = common::sources(0.05, 3e6, STRENGTH);
    let response = normal.evaluate_contact_load(0.0).unwrap();
    let source = shear
        .sliding_patch_dmt(normal.normal_state(), &response)
        .unwrap();
    let mut f = fixture(response.contact_radius_m, false);
    let mut contact = RigidAdhesiveContact {
        binding: f.binding,
        frame: f.frame,
        source,
    };
    // Translation plus spin with this unrelated inertia is a curved path. A
    // short interval cannot authorize substituting a straight stop-and-hold.
    f.bodies[0] = state(
        f.bodies[0].pose(),
        f.properties[0],
        Vec3::new(1e-3, 0.0, 0.0),
        1000.0,
    );
    let mut bodies = f.bodies;
    assert!(matches!(
        apply_rigid_adhesive_trajectory(
            &mut bodies,
            &f.properties,
            contact,
            1e-15,
            controls(),
            &CancelGate::new()
        ),
        Err(RigidAdhesiveTrajectoryError::Relaxation(
            AdhesiveRelaxationError::CurvedTrajectory { .. }
        ))
    ));
    assert_eq!(bodies, f.bodies);
    f.bodies[0] = state(
        f.bodies[0].pose(),
        f.properties[0],
        Vec3::new(1e-3, 0.0, 0.0),
        0.0,
    );
    bodies = f.bodies;
    let gate = CancelGate::new();
    gate.request();
    assert!(
        apply_rigid_adhesive_trajectory(
            &mut bodies,
            &f.properties,
            contact,
            1e-6,
            controls(),
            &gate
        )
        .is_err()
    );
    assert_eq!(bodies, f.bodies);
    contact.binding.arm_a_body = Vec3::new(0.0, 0.0, 1e-9);
    let permissive_geometry = RigidAdhesiveControls {
        geometry_tolerance_m: 2e-9,
        ..controls()
    };
    assert!(matches!(
        apply_rigid_adhesive_trajectory(
            &mut bodies,
            &f.properties,
            contact,
            1e-6,
            permissive_geometry,
            &CancelGate::new()
        ),
        Err(RigidAdhesiveTrajectoryError::Audit(_))
    ));
    assert_eq!(bodies, f.bodies);
    contact.binding = f.binding;
    let first = prepare_rigid_adhesive_trajectory(
        &bodies,
        &f.properties,
        contact,
        1e-6,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let retry = prepare_rigid_adhesive_trajectory(
        &bodies,
        &f.properties,
        contact,
        1e-6,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(first, retry);
    assert_eq!(bodies, f.bodies);
}
