//! Sourced adhesive disk traction changes both original finite body momenta.
//!
//! G1: independent pure-spin and combined rim-slip solutions with native work.
//! G3: tangent-gauge covariance. G4: cancellation and late-refusal rollback.
#[path = "../../fs-contact/tests/adhesive_dmt_coupled/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::{RigidContactBinding, adhesive::*};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_tribo::adhesive::spinning::contact::{AdhesiveWrenchConfig, AdhesiveWrenchMode};
use fs_tribo::partial_slip::TangentFrame;

const LENGTH: f64 = 1e-6;
const MASS: f64 = 1e-15;
const INERTIA: f64 = MASS * LENGTH * LENGTH;
const STRENGTH: f64 = 2e6;

fn close(actual: f64, expected: f64, absolute: f64) {
    common::close(actual, expected, absolute, 2e-9);
}
fn close_vec(actual: Vec3, expected: Vec3, tolerance: f64) {
    close(actual.x, expected.x, tolerance);
    close(actual.y, expected.y, tolerance);
    close(actual.z, expected.z, tolerance);
}
fn norm(value: Vec3) -> f64 {
    value.x.hypot(value.y).hypot(value.z)
}
fn vector(value: [f64; 3]) -> Vec3 {
    Vec3::new(value[0], value[1], value[2])
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
fn state(pose: Pose, mass: MassProperties, velocity: Vec3, omega: Vec3) -> RigidBodyState {
    let omega = pose.orientation().rotate_world_to_body(omega);
    let inertia = mass.principal_inertia_body();
    RigidBodyState::new(
        pose,
        velocity.scale(mass.mass()),
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
}
fn axial(spin: f64) -> Fixture {
    let properties = [
        MassProperties::new(
            2.0 * MASS,
            Vec3::ZERO,
            Vec3::new(0.8, 1.0, 1.2).scale(INERTIA),
        )
        .unwrap(),
        MassProperties::new(
            3.0 * MASS,
            Vec3::ZERO,
            Vec3::new(1.5, 2.0, 2.5).scale(INERTIA),
        )
        .unwrap(),
    ];
    let arms = [
        Vec3::new(0.0, 0.0, -0.2 * LENGTH),
        Vec3::new(0.0, 0.0, 0.15 * LENGTH),
    ];
    let poses = arms.map(|arm| Pose::new(arm.scale(-1.0), UnitQuaternion::IDENTITY).unwrap());
    Fixture {
        bodies: [
            state(
                poses[0],
                properties[0],
                Vec3::ZERO,
                Vec3::new(0.0, 0.0, spin),
            ),
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
    }
}

#[test]
fn common_normal_translation_does_not_hide_source_owned_junction_heat() {
    let (normal, shear) = common::sources(0.05, 3e6, STRENGTH);
    let response = normal.evaluate_contact_load(0.0).unwrap();
    let source = shear
        .sliding_patch_dmt(normal.normal_state(), &response)
        .unwrap();
    let f = axial(1000.0);
    let ia = f.properties[0].principal_inertia_body().z;
    let ib = f.properties[1].principal_inertia_body().z;
    let inverse = 1.0 / ia + 1.0 / ib;
    let moment = (2.0 / 3.0) * STRENGTH * response.contact_area_m2 * response.contact_radius_m;
    let dt = 0.25 * 1000.0 / (moment * inverse);
    let run = |boost: f64| {
        let mut bodies = core::array::from_fn::<_, 2, _>(|k| {
            RigidBodyState::new(
                f.bodies[k].pose(),
                Vec3::new(0.0, 0.0, boost * f.properties[k].mass()),
                f.bodies[k].angular_momentum_body(),
            )
            .unwrap()
        });
        let report = apply_rigid_adhesive(
            &mut bodies,
            &f.properties,
            RigidAdhesiveContact {
                binding: f.binding,
                frame: f.frame,
                source,
            },
            dt,
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
        (bodies, report)
    };
    let reference = run(0.0);
    let boosted = run(2.0_f64.powi(20));
    assert_eq!(reference.1.contact, boosted.1.contact);
    assert_eq!(reference.1.kinetic_change_j, boosted.1.kinetic_change_j);
    assert_eq!(boosted.1.kinetic_before_j, boosted.1.kinetic_after_j);
    assert!(boosted.1.dissipated_energy_j > 0.0);
    close(
        boosted.1.kinetic_change_j,
        -boosted.1.dissipated_energy_j,
        1e-33,
    );
    for k in 0..2 {
        assert_eq!(
            reference.0[k].angular_momentum_body(),
            boosted.0[k].angular_momentum_body()
        );
    }
}

#[test]
fn zero_and_tensile_normal_load_spin_both_finite_bodies_and_dissipate_only_relative_energy() {
    let (normal, shear) = common::sources(0.05, 3e6, STRENGTH);
    let pull_off =
        2.0 * PI * normal.effective_radius_m() * normal.normal_state().adhesion_energy_j_per_m2();
    for load in [0.0, -0.5 * pull_off] {
        let response = normal.evaluate_contact_load(load).unwrap();
        let source = shear
            .sliding_patch_dmt(normal.normal_state(), &response)
            .unwrap();
        let f = axial(1000.0);
        let ia = f.properties[0].principal_inertia_body().z;
        let ib = f.properties[1].principal_inertia_body().z;
        let inverse = 1.0 / ia + 1.0 / ib;
        let moment = (2.0 / 3.0) * STRENGTH * response.contact_area_m2 * response.contact_radius_m;
        let stop_time = 1000.0 / (moment * inverse);
        for fraction in [0.25, 1.0] {
            let dt = fraction * stop_time;
            let mut bodies = f.bodies;
            let report = apply_rigid_adhesive(
                &mut bodies,
                &f.properties,
                RigidAdhesiveContact {
                    binding: f.binding,
                    frame: f.frame,
                    source,
                },
                dt,
                controls(),
                &CancelGate::new(),
            )
            .unwrap();
            let contact = report.contact.unwrap();
            assert_eq!(contact.mode, AdhesiveWrenchMode::Sliding);
            assert_eq!(contact.force_n, [0.0; 2]);
            close(contact.spin_moment_nm, -moment, 1e-29);
            close(report.spin_after_rad_s, (1.0 - fraction) * 1000.0, 2e-7);
            close(
                bodies[0].angular_momentum_body().z,
                ia * 1000.0 - moment * dt,
                1e-35,
            );
            close(bodies[1].angular_momentum_body().z, moment * dt, 1e-35);
            assert!(bodies[1].angular_momentum_body().z > 0.0);
            let heat = (fraction - 0.5 * fraction * fraction) * 1000.0_f64.powi(2) / inverse;
            close(report.dissipated_energy_j, heat, 1e-33);
            close(
                report.kinetic_after_j - report.kinetic_before_j,
                -heat,
                1e-33,
            );
            close(report.midpoint_impulse_work_j, -heat, 1e-33);
            assert!(report.energy_residual_j.abs() < 1e-33);
            assert!(norm(report.linear_momentum_change_ns) < 1e-32);
            assert!(norm(report.angular_momentum_change_nms) < 1e-35);
            for i in 0..2 {
                assert_eq!(bodies[i].pose(), f.bodies[i].pose());
            }
        }
        // A midpoint solve can still dissipate while its frozen moment drives
        // the endpoint backwards. That interval needs the real stopping event.
        let mut bodies = f.bodies;
        assert!(matches!(
            apply_rigid_adhesive(
                &mut bodies,
                &f.properties,
                RigidAdhesiveContact {
                    binding: f.binding,
                    frame: f.frame,
                    source,
                },
                1.5 * stop_time,
                controls(),
                &CancelGate::new(),
            ),
            Err(RigidAdhesiveError::UnresolvedSlipEvent { .. })
        ));
        assert_eq!(bodies, f.bodies);
    }
}

// Independent inverse-inertia bilinear form; the production Gram matrix and
// finite-body impulse helpers are not used by this reference.
fn mobility(f: &Fixture, radius: f64) -> [[f64; 3]; 3] {
    let axes = [
        vector(f.frame.longitudinal()),
        vector(f.frame.lateral()),
        vector(f.frame.normal()),
    ];
    let mut result = [[0.0; 3]; 3];
    for (i, arm) in [(0, f.binding.arm_a_body), (1, f.binding.arm_b_body)] {
        let q = f.bodies[i].pose().orientation();
        let torques = [
            arm.cross(q.rotate_world_to_body(axes[0])),
            arm.cross(q.rotate_world_to_body(axes[1])),
            q.rotate_world_to_body(axes[2]).scale(radius),
        ];
        let inertia = f.properties[i].principal_inertia_body();
        for k in 0..3 {
            for l in 0..3 {
                let t = torques[l];
                result[k][l] +=
                    torques[k].dot(Vec3::new(t.x / inertia.x, t.y / inertia.y, t.z / inertia.z));
                if k < 2 && l < 2 {
                    result[k][l] += axes[k].dot(axes[l]) / f.properties[i].mass();
                }
            }
        }
    }
    result
}

// Prescribe midpoint |u|=a*omega. The independent uniform-disk rim solution is
// F/(tau*A)=8/(3*pi), M/(tau*A*a)=8/(9*pi).
fn coupled(radius: f64, capacity: f64, dt: f64) -> (Fixture, Vec3, f64) {
    let mut f = axial(0.0);
    let orientations = [
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -1.0), 0.63).unwrap(),
        UnitQuaternion::from_axis_angle(Vec3::new(-2.0, 1.0, 1.0), -0.43).unwrap(),
    ];
    let arms = [
        Vec3::new(0.13, -0.08, -0.20).scale(LENGTH),
        Vec3::new(-0.06, 0.11, 0.18).scale(LENGTH),
    ];
    let poses: [Pose; 2] = core::array::from_fn(|i| {
        Pose::new(
            orientations[i].rotate_body_to_world(arms[i]).scale(-1.0),
            orientations[i],
        )
        .unwrap()
    });
    f.binding.arm_a_body = arms[0];
    f.binding.arm_b_body = arms[1];
    f.bodies = core::array::from_fn(|i| state(poses[i], f.properties[i], Vec3::ZERO, Vec3::ZERO));
    let force = Vec3::new(-0.6, 0.8, 0.0).scale(capacity * 8.0 / (3.0 * PI));
    let moment = -capacity * radius * 8.0 / (9.0 * PI);
    let response = common::mul(mobility(&f, radius), [force.x, force.y, moment / radius]);
    let edge_speed = radius * 1000.0;
    let initial = [
        0.6 * edge_speed - 0.5 * dt * response[0],
        -0.8 * edge_speed - 0.5 * dt * response[1],
        edge_speed - 0.5 * dt * response[2],
    ];
    let omega = Vec3::new(0.0, 0.0, initial[2] / radius);
    let velocity = Vec3::new(initial[0], initial[1], 0.0)
        .sub(omega.cross(orientations[0].rotate_body_to_world(arms[0])));
    f.bodies[0] = state(poses[0], f.properties[0], velocity, omega);
    (f, force, moment)
}

fn independent_endpoint(f: &Fixture, force: Vec3, moment: f64, dt: f64) -> [RigidBodyState; 2] {
    core::array::from_fn(|i| {
        let before = f.bodies[i];
        let q = before.pose().orientation();
        let arm = [f.binding.arm_a_body, f.binding.arm_b_body][i];
        let sign = [1.0, -1.0][i];
        let impulse = force.scale(sign * dt);
        let angular = arm
            .cross(q.rotate_world_to_body(impulse))
            .add(q.rotate_world_to_body(Vec3::new(0.0, 0.0, sign * moment * dt)));
        RigidBodyState::new(
            before.pose(),
            before.linear_momentum_world().add(impulse),
            before.angular_momentum_body().add(angular),
        )
        .unwrap()
    })
}

#[test]
fn combined_junction_shear_spin_uses_all_original_inertia_couplings_and_actual_work() {
    let (normal, shear) = common::sources(0.05, 3e6, STRENGTH);
    let response = normal.evaluate_contact_load(0.0).unwrap();
    let source = shear
        .sliding_patch_dmt(normal.normal_state(), &response)
        .unwrap();
    let dt = 1e-12;
    let (f, force, moment) = coupled(
        response.contact_radius_m,
        STRENGTH * response.contact_area_m2,
        dt,
    );
    let expected = independent_endpoint(&f, force, moment, dt);
    let matrix = mobility(&f, response.contact_radius_m);
    assert!(matrix[0][1].abs() > 1e10 && matrix[0][2].abs() > 1e10 && matrix[1][2].abs() > 1e10);
    let mut bodies = f.bodies;
    let report = apply_rigid_adhesive(
        &mut bodies,
        &f.properties,
        RigidAdhesiveContact {
            binding: f.binding,
            frame: f.frame,
            source,
        },
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let actual_matrix = report.scaled_mobility_per_kg.unwrap();
    for i in 0..3 {
        for j in 0..3 {
            close(actual_matrix[i][j], matrix[i][j], 1.0);
        }
    }
    let contact = report.contact.unwrap();
    close(contact.force_n[0], force.x, 2e-19);
    close(contact.force_n[1], force.y, 2e-19);
    close(contact.spin_moment_nm, moment, 2e-27);
    for i in 0..2 {
        close_vec(
            bodies[i].linear_momentum_world(),
            expected[i].linear_momentum_world(),
            3e-31,
        );
        close_vec(
            bodies[i].angular_momentum_body(),
            expected[i].angular_momentum_body(),
            3e-37,
        );
        assert_eq!(bodies[i].pose(), f.bodies[i].pose());
    }
    let speed = response.contact_radius_m * 1000.0;
    let work = dt * (force.dot(Vec3::new(0.6 * speed, -0.8 * speed, 0.0)) + moment * 1000.0);
    close(report.midpoint_impulse_work_j, work, 3e-35);
    close(report.dissipated_energy_j, -work, 3e-35);
    assert!(report.force_residual_n <= controls().contact.force_tolerance_n);
    assert!(report.spin_moment_residual_nm <= controls().contact.moment_tolerance_nm);
    assert!(report.energy_residual_j.abs() < 1e-33);
    assert!(norm(report.angular_momentum_change_nms) < 1e-35);
}

#[test]
fn changing_tangent_gauge_preserves_the_physical_two_body_wrench() {
    let (normal, shear) = common::sources(0.05, 3e6, STRENGTH);
    let response = normal.evaluate_contact_load(0.0).unwrap();
    let source = shear
        .sliding_patch_dmt(normal.normal_state(), &response)
        .unwrap();
    let (f, _, _) = coupled(
        response.contact_radius_m,
        STRENGTH * response.contact_area_m2,
        1e-12,
    );
    let mut a = f.bodies;
    let mut b = f.bodies;
    let run = |bodies: &mut [RigidBodyState], frame| {
        apply_rigid_adhesive(
            bodies,
            &f.properties,
            RigidAdhesiveContact {
                binding: f.binding,
                frame,
                source,
            },
            1e-12,
            controls(),
            &CancelGate::new(),
        )
        .unwrap()
    };
    let first = run(&mut a, f.frame);
    let rotated = run(
        &mut b,
        TangentFrame::new([0.0, 0.0, 1.0], [0.0, 1.0, 0.0]).unwrap(),
    );
    for i in 0..2 {
        close_vec(
            a[i].linear_momentum_world(),
            b[i].linear_momentum_world(),
            3e-31,
        );
        close_vec(
            a[i].angular_momentum_body(),
            b[i].angular_momentum_body(),
            3e-37,
        );
    }
    close_vec(
        first.force_impulse_world_ns,
        rotated.force_impulse_world_ns,
        3e-31,
    );
    close_vec(
        first.spin_angular_impulse_world_nms,
        rotated.spin_angular_impulse_world_nms,
        3e-37,
    );
    close(
        first.dissipated_energy_j,
        rotated.dissipated_energy_j,
        3e-35,
    );
}

#[test]
fn attractive_noncontact_gap_is_a_noop_even_for_separated_moving_material_points() {
    let (normal, shear) = common::sources(0.05, 3e6, STRENGTH);
    let response = normal.evaluate(-0.5 * common::Z).unwrap();
    assert!(response.normal_force_n < 0.0 && response.contact_area_m2 == 0.0);
    let source = shear
        .sliding_patch_dmt(normal.normal_state(), &response)
        .unwrap();
    let f = axial(1000.0);
    let mut separated = f.binding;
    separated.arm_a_body = separated.arm_a_body.add(Vec3::new(0.2 * LENGTH, 0.0, 0.0));
    let mut bodies = f.bodies;
    let report = apply_rigid_adhesive(
        &mut bodies,
        &f.properties,
        RigidAdhesiveContact {
            binding: separated,
            frame: f.frame,
            source,
        },
        1e-12,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(bodies, f.bodies);
    assert!(report.contact.is_none() && report.scaled_mobility_per_kg.is_none());
    assert_eq!(report.force_impulse_world_ns, Vec3::ZERO);
    assert_eq!(report.spin_angular_impulse_world_nms, Vec3::ZERO);
    assert_eq!(report.dissipated_energy_j, 0.0);
    assert_eq!(report.kinetic_before_j, report.kinetic_after_j);
}

#[test]
fn cancellation_and_late_angular_momentum_refusal_preserve_both_original_bodies() {
    let (normal, shear) = common::sources(0.05, 3e6, STRENGTH);
    let response = normal.evaluate_contact_load(0.0).unwrap();
    let source = shear
        .sliding_patch_dmt(normal.normal_state(), &response)
        .unwrap();
    let (f, _, _) = coupled(
        response.contact_radius_m,
        STRENGTH * response.contact_area_m2,
        1e-12,
    );
    let contact = RigidAdhesiveContact {
        binding: f.binding,
        frame: f.frame,
        source,
    };
    let gate = CancelGate::new();
    gate.request();
    let mut bodies = f.bodies;
    assert!(
        apply_rigid_adhesive(
            &mut bodies,
            &f.properties,
            contact,
            1e-12,
            controls(),
            &gate
        )
        .is_err()
    );
    assert_eq!(bodies, f.bodies);
    // Proximity admission does not excuse a real net angular impulse from
    // equal/opposite shear applied at separated material points.
    let mut separated = contact;
    separated.binding.arm_a_body = separated.binding.arm_a_body.add(Vec3::new(0.0, 0.0, 1e-9));
    let mut permissive_geometry = controls();
    permissive_geometry.geometry_tolerance_m = 2e-9;
    permissive_geometry.angular_momentum_tolerance_nms = 1e-33;
    assert!(matches!(
        apply_rigid_adhesive(
            &mut bodies,
            &f.properties,
            separated,
            1e-12,
            permissive_geometry,
            &CancelGate::new()
        ),
        Err(RigidAdhesiveError::Audit(_))
    ));
    assert_eq!(bodies, f.bodies);
    let first = prepare_rigid_adhesive(
        &bodies,
        &f.properties,
        contact,
        1e-12,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let retry = prepare_rigid_adhesive(
        &bodies,
        &f.properties,
        contact,
        1e-12,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(first, retry);
    assert_eq!(bodies, f.bodies);
}
