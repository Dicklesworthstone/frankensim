//! G1/G3/G4: bonded torsion changes original finite-body momenta and closes work.
#[path = "../../fs-contact/tests/source_point/common.rs"]
mod common;

#[path = "adhesive_torsion/normal.rs"]
mod normal;

use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::torsion::{JkrTorsionControls, SourceJkrTorsionContact};
use fs_contact::normal_patch::adhesive::{JkrEnvelope, SourceJkrSphereContact};
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::{RigidContactBinding, adhesive::torsion::*};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3, Wrench};
use fs_tribo::partial_slip::TangentFrame;

fn source(tensile: bool) -> SourceJkrTorsionContact {
    let a = common::solid("torsion-A", 1000.0, 1e7);
    let b = common::solid("torsion-B", 2000.0, 2e9);
    let normal = SourceJkrSphereContact::new(
        common::bind(&a, &b, None, 0.05).normal_state().clone(),
        0.01,
        1e-9,
        JkrEnvelope {
            maximum_indentation_m: 1e-5,
            half_space_depth_m: 0.1,
            layer_thickness_m: 0.1,
            maximum_patch_to_radius: 0.1,
            maximum_patch_to_depth: 0.1,
            maximum_patch_to_layer: 0.1,
            maximum_nominal_strain: 0.1,
            maximum_repulsive_pressure_to_yield: 0.1,
            maximum_process_zone_to_patch: 0.01,
            minimum_tabor_parameter: 5.0,
        },
    )
    .unwrap();
    let load = if tensile {
        -0.75 * PI * 0.05 * 0.01
    } else {
        0.0
    };
    let approach = normal
        .evaluate_force_controlled(load)
        .unwrap()
        .indentation_m;
    SourceJkrTorsionContact::new(
        normal,
        approach,
        JkrTorsionControls {
            critical_twist_angle_rad: 1e-5,
            maximum_critical_twist_angle_rad: 0.01,
            maximum_mode_iii_energy_release_to_adhesion: 0.01,
        },
    )
    .unwrap()
}
fn controls() -> RigidAdhesiveTorsionControls {
    RigidAdhesiveTorsionControls {
        maximum_bodies: 4,
        maximum_step_s: 0.01,
        maximum_angle_increment_rad: 1e-3,
        maximum_iterations: 96,
        angle_tolerance_rad: 1e-18,
        geometry_tolerance_m: 1e-15,
        angular_velocity_tolerance_rad_s: 1e-10,
        linear_momentum_tolerance_ns: 1e-25,
        angular_momentum_tolerance_nms: 1e-25,
        energy_tolerance_j: 1e-25,
    }
}
fn close(a: f64, b: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= absolute + 2e-9 * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}
fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}
fn vec_close(a: Vec3, b: Vec3, tolerance: f64) {
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
impl Fixture {
    fn contact<'a>(&self, source: &'a SourceJkrTorsionContact) -> RigidAdhesiveTorsionContact<'a> {
        RigidAdhesiveTorsionContact {
            binding: self.binding,
            frame: self.frame,
            source,
        }
    }
    fn mobility(&self) -> f64 {
        let axis = vector(self.frame.normal());
        (0..2)
            .map(|i| {
                let body_axis = self.bodies[i]
                    .pose()
                    .orientation()
                    .rotate_world_to_body(axis);
                body_axis.dot(
                    self.properties[i]
                        .angular_velocity_body_checked(body_axis)
                        .unwrap(),
                )
            })
            .sum()
    }
}
fn fixture(rate: f64, boost: f64, rotated: bool) -> Fixture {
    let frame = if rotated {
        TangentFrame::new([-2.0, 1.0, 3.0], [1.0, 2.0, 0.0]).unwrap()
    } else {
        TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap()
    };
    let axis = vector(frame.normal());
    let orientation = if rotated {
        [
            UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), 0.63).unwrap(),
            UnitQuaternion::from_axis_angle(Vec3::new(0.0, 1.0, 0.0), -0.41).unwrap(),
        ]
    } else {
        [UnitQuaternion::IDENTITY; 2]
    };
    let principal = if rotated {
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 0.0, 0.0), 0.47).unwrap()
    } else {
        UnitQuaternion::IDENTITY
    };
    let properties = [
        MassProperties::from_principal_frame(
            2e-5,
            Vec3::ZERO,
            Vec3::new(0.8e-9, 1e-9, 1.2e-9),
            principal,
        )
        .unwrap(),
        MassProperties::from_principal_frame(
            3e-5,
            Vec3::ZERO,
            Vec3::new(1.5e-9, 2e-9, 2.5e-9),
            principal,
        )
        .unwrap(),
    ];
    let arms = [Vec3::new(0.0, 0.0, 0.01), Vec3::new(0.0, 0.0, -0.01)];
    let bodies = core::array::from_fn(|i| {
        let pose = Pose::new(
            orientation[i].rotate_body_to_world(arms[i]).scale(-1.0),
            orientation[i],
        )
        .unwrap();
        let omega =
            orientation[i].rotate_world_to_body(axis.scale(if i == 0 { rate } else { 0.0 }));
        RigidBodyState::new(
            pose,
            Vec3::new(0.0, 0.0, boost * properties[i].mass()),
            properties[i].angular_momentum_body_checked(omega).unwrap(),
        )
        .unwrap()
    });
    Fixture {
        bodies,
        properties,
        binding: RigidContactBinding {
            body_a: 0,
            body_b: 1,
            arm_a_body: arms[0],
            arm_b_body: arms[1],
        },
        frame,
    }
}

#[test]
fn g1_elastic_pair_uses_normal_axis_and_full_native_inertia() {
    let source = source(false);
    let c = source.law().parameters().critical_twist_angle_rad;
    let k = source.law().stiffness_nm_per_rad();
    let dt = 1e-4;
    let rate = 0.25 * c / dt;
    for rotated in [false, true] {
        let f = fixture(rate, 0.0, rotated);
        let contact = f.contact(&source);
        let state = contact
            .initial_state(&f.bodies, &f.properties, controls(), &CancelGate::new())
            .unwrap();
        let candidate = prepare_rigid_adhesive_torsion(
            &f.bodies,
            &f.properties,
            &state,
            contact,
            dt,
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
        let angle = dt * rate / (1.0 + 0.25 * dt * dt * f.mobility() * k);
        let mean = -0.5 * k * angle;
        let r = &candidate.report;
        close(r.angle_increment_rad, angle, 1e-18);
        close(r.rotational_mobility_per_kg_m2, f.mobility(), 1e-5);
        close(r.contact.torsion().mean_moment_n_m, mean, 1e-24);
        close(
            r.twisting_rate_after_rad_s,
            rate + dt * f.mobility() * mean,
            1e-12,
        );
        close(r.kinetic_energy_change_j, -0.5 * k * angle * angle, 1e-25);
        assert_eq!(r.contact.torsion().irreversible_loss_j, 0.0);
        assert_eq!(r.evaluations, 1);
        vec_close(
            r.angular_impulse_on_a_world_nms,
            vector(f.frame.normal()).scale(dt * mean),
            1e-25,
        );
        close(
            r.angular_impulse_on_a_world_nms
                .dot(vector(f.frame.longitudinal())),
            0.0,
            1e-25,
        );
        for i in 0..2 {
            let expected = f.bodies[i].angular_momentum_body().add(
                f.bodies[i].pose().orientation().rotate_world_to_body(
                    r.angular_impulse_on_a_world_nms
                        .scale(if i == 0 { 1.0 } else { -1.0 }),
                ),
            );
            vec_close(
                candidate.bodies_after[i].angular_momentum_body(),
                expected,
                1e-25,
            );
            assert_eq!(candidate.bodies_after[i].pose(), f.bodies[i].pose());
            assert_eq!(
                candidate.bodies_after[i].linear_momentum_world(),
                f.bodies[i].linear_momentum_world()
            );
        }
        vec_close(r.angular_momentum_change_nms, Vec3::ZERO, 1e-25);
    }
}

#[test]
fn g1_tensile_junction_yields_and_common_translation_cannot_hide_its_loss() {
    let source = source(true);
    assert!(source.normal_response().normal_force_n < 0.0);
    let c = source.law().parameters().critical_twist_angle_rad;
    let dt = 1e-4;
    let run = |boost| {
        let f = fixture(4.0 * c / dt, boost, false);
        let contact = f.contact(&source);
        let mut state = contact
            .initial_state(&f.bodies, &f.properties, controls(), &CancelGate::new())
            .unwrap();
        let mut bodies = f.bodies;
        let report = apply_rigid_adhesive_torsion(
            &mut bodies,
            &f.properties,
            &mut state,
            contact,
            dt,
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
        close(state.elastic_twist_rad(), c, 1e-20);
        assert_ne!(
            bodies[0].angular_momentum_body(),
            f.bodies[0].angular_momentum_body()
        );
        assert_ne!(
            bodies[1].angular_momentum_body(),
            f.bodies[1].angular_momentum_body()
        );
        report
    };
    let rest = run(0.0);
    let moving = run(1e12);
    let law = rest.contact.torsion();
    assert!(law.irreversible_loss_j > 0.0 && law.storage_after_j > 0.0);
    close(
        law.endpoint_moment_n_m,
        -source.law().critical_moment_n_m(),
        1e-23,
    );
    close(
        -rest.kinetic_energy_change_j,
        law.storage_change_j + law.irreversible_loss_j,
        1e-25,
    );
    assert_eq!(rest.contact, moving.contact);
    assert_eq!(rest.kinetic_energy_change_j, moving.kinetic_energy_change_j);
}

#[test]
fn g3_elastic_turns_are_admitted_but_hidden_and_ambiguous_yield_refuse() {
    let source = source(false);
    let c = source.law().parameters().critical_twist_angle_rad;
    let k = source.law().stiffness_nm_per_rad();
    for rotated in [false, true] {
        let frequency = (k * fixture(0.0, 0.0, rotated).mobility()).sqrt();
        let dt = 4.0 / frequency;
        // Large authored audit tolerances cannot waive unresolved physical yield.
        let controls = RigidAdhesiveTorsionControls {
            maximum_step_s: 1.1 * dt,
            energy_tolerance_j: 1.0,
            ..controls()
        };
        for speed in [0.8, 1.0, 1.1] {
            let f = fixture(speed * c * frequency, 1e12, rotated);
            let contact = f.contact(&source);
            let initial = contact
                .initial_state(&f.bodies, &f.properties, controls, &CancelGate::new())
                .unwrap();
            let mut state = initial;
            let mut bodies = f.bodies;
            let result = apply_rigid_adhesive_torsion(
                &mut bodies,
                &f.properties,
                &mut state,
                contact,
                dt,
                controls,
                &CancelGate::new(),
            );
            if speed == 0.8 {
                let report = result.unwrap();
                assert!(report.twisting_rate_after_rad_s < 0.0);
                close(report.angle_increment_rad, 0.64 * c, 1e-18);
                close(
                    report.twisting_rate_after_rad_s,
                    -0.48 * c * frequency,
                    1e-12,
                );
                close(report.total_energy_residual_j, 0.0, 1e-25);
                assert_eq!(report.contact.torsion().irreversible_loss_j, 0.0);
                assert_eq!(report.evaluations, 1);
            } else {
                assert!(
                    matches!(
                        result,
                        Err(RigidAdhesiveTorsionError::UnresolvedTorsionReversal)
                    ),
                    "{result:?}"
                );
                assert_eq!(bodies, f.bodies);
                assert_eq!(state, initial);
            }
        }
    }
}

#[test]
fn g4_cancel_foreign_source_geometry_and_budget_refuse_atomically() {
    let foreign = source(true);
    let source = source(false);
    let f = fixture(0.4, 0.0, false);
    let contact = f.contact(&source);
    let initial = contact
        .initial_state(&f.bodies, &f.properties, controls(), &CancelGate::new())
        .unwrap();
    let mut state = initial;
    let mut bodies = f.bodies;
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        apply_rigid_adhesive_torsion(
            &mut bodies,
            &f.properties,
            &mut state,
            contact,
            1e-4,
            controls(),
            &cancelled
        )
        .is_err()
    );
    assert_eq!(state, initial);
    assert_eq!(bodies, f.bodies);
    let changed = RigidAdhesiveTorsionContact {
        source: &foreign,
        ..contact
    };
    assert!(matches!(
        apply_rigid_adhesive_torsion(
            &mut bodies,
            &f.properties,
            &mut state,
            changed,
            1e-4,
            controls(),
            &CancelGate::new()
        ),
        Err(RigidAdhesiveTorsionError::Source(_))
    ));
    let changed = RigidAdhesiveTorsionContact {
        frame: TangentFrame::new([0.0, 1.0, 0.0], [1.0, 0.0, 0.0]).unwrap(),
        ..contact
    };
    assert!(matches!(
        apply_rigid_adhesive_torsion(
            &mut bodies,
            &f.properties,
            &mut state,
            changed,
            1e-4,
            controls(),
            &CancelGate::new()
        ),
        Err(RigidAdhesiveTorsionError::GeometryMismatch)
    ));
    let restricted = RigidAdhesiveTorsionControls {
        maximum_angle_increment_rad: 1e-12,
        ..controls()
    };
    assert!(matches!(
        apply_rigid_adhesive_torsion(
            &mut bodies,
            &f.properties,
            &mut state,
            contact,
            1e-4,
            restricted,
            &CancelGate::new()
        ),
        Err(RigidAdhesiveTorsionError::StepLimit)
    ));
    let short = RigidAdhesiveTorsionControls {
        maximum_iterations: 1,
        angle_tolerance_rad: 1e-30,
        ..controls()
    };
    assert!(matches!(
        apply_rigid_adhesive_torsion(
            &mut bodies,
            &f.properties,
            &mut state,
            contact,
            1e-4,
            short,
            &CancelGate::new()
        ),
        Err(RigidAdhesiveTorsionError::Unresolved { .. })
    ));
    assert_eq!(state, initial);
    assert_eq!(bodies, f.bodies);
}

#[test]
fn g1_loaded_torsion_uses_actual_inertia_and_accounts_for_applied_work() {
    let source = source(true);
    let dt = 1e-4;
    let critical = source.law().parameters().critical_twist_angle_rad;
    let stiffness = source.law().stiffness_nm_per_rad();
    let f = fixture(0.1 * critical / dt, 0.0, true);
    let axis = vector(f.frame.normal());
    let contact = f.contact(&source);
    let state = contact
        .initial_state(&f.bodies, &f.properties, controls(), &CancelGate::new())
        .unwrap();
    let wrenches = [
        Wrench {
            force_world: Vec3::new(1e-8, -2e-8, 3e-8),
            torque_body: Vec3::new(2e-8, -1e-8, 3e-8),
        },
        Wrench {
            force_world: Vec3::new(-2e-8, 3e-8, 1e-8),
            torque_body: Vec3::new(-1e-8, 2e-8, -2e-8),
        },
    ];
    let acceleration: [Vec3; 2] = core::array::from_fn(|i| {
        f.bodies[i].pose().orientation().rotate_body_to_world(
            f.properties[i]
                .angular_velocity_body_checked(wrenches[i].torque_body)
                .unwrap(),
        )
    });
    let alpha = acceleration[0].sub(acceleration[1]).dot(axis);
    let expected = (0.1 * critical + 0.5 * dt * dt * alpha)
        / (1.0 + 0.25 * dt * dt * f.mobility() * stiffness);
    let candidate = prepare_rigid_adhesive_torsion_with_wrenches(
        &f.bodies,
        &f.properties,
        &state,
        contact,
        wrenches,
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    close(candidate.report.angle_increment_rad, expected, 1e-18);
    close(
        candidate.report.contact.torsion().mean_moment_n_m,
        -0.5 * stiffness * expected,
        1e-23,
    );
    assert_eq!(candidate.report.contact.torsion().irreversible_loss_j, 0.0);
    let mut work = 0.0;
    let mut angular_impulse = Vec3::ZERO;
    for i in 0..2 {
        let before = f.bodies[i];
        let after = candidate.bodies_after[i];
        let mean_velocity = before
            .center_of_mass_velocity_world(f.properties[i])
            .unwrap()
            .add(
                after
                    .center_of_mass_velocity_world(f.properties[i])
                    .unwrap(),
            )
            .scale(0.5);
        let mean_omega = f.properties[i]
            .angular_velocity_body_checked(before.angular_momentum_body())
            .unwrap()
            .add(
                f.properties[i]
                    .angular_velocity_body_checked(after.angular_momentum_body())
                    .unwrap(),
            )
            .scale(0.5);
        work += dt
            * (wrenches[i].force_world.dot(mean_velocity)
                + wrenches[i].torque_body.dot(mean_omega));
        angular_impulse = angular_impulse
            .add(
                before
                    .pose()
                    .position_world()
                    .cross(wrenches[i].force_world.scale(dt)),
            )
            .add(
                before
                    .pose()
                    .orientation()
                    .rotate_body_to_world(wrenches[i].torque_body.scale(dt)),
            );
    }
    close(candidate.report.external_work_j, work, 1e-25);
    close(
        candidate.report.kinetic_energy_change_j
            + candidate.report.contact.torsion().storage_change_j,
        work,
        1e-25,
    );
    vec_close(
        candidate.report.angular_momentum_change_nms,
        angular_impulse,
        1e-25,
    );
    vec_close(
        candidate.report.angular_momentum_residual_nms,
        Vec3::ZERO,
        1e-25,
    );
}

#[test]
fn g1_loaded_twist_reaches_yield_and_invalid_loads_cannot_publish() {
    let source = source(false);
    let f = fixture(0.0, 0.0, false);
    let contact = f.contact(&source);
    let initial = contact
        .initial_state(&f.bodies, &f.properties, controls(), &CancelGate::new())
        .unwrap();
    let mut bodies = f.bodies;
    let mut state = initial;
    let dt = 1e-3;
    let alpha = 6.0 * source.law().parameters().critical_twist_angle_rad / (dt * dt)
        + 2.0 * f.mobility() * source.law().critical_moment_n_m();
    let wrenches = [
        Wrench {
            force_world: Vec3::ZERO,
            torque_body: f.properties[0]
                .angular_momentum_body_checked(Vec3::new(0.0, 0.0, alpha))
                .unwrap(),
        },
        Wrench::ZERO,
    ];
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        apply_rigid_adhesive_torsion_with_wrenches(
            &mut bodies,
            &f.properties,
            &mut state,
            contact,
            wrenches,
            dt,
            controls(),
            &cancelled,
        )
        .is_err()
    );
    assert_eq!(bodies, f.bodies);
    assert_eq!(state, initial);
    let invalid = [
        Wrench {
            force_world: Vec3::new(f64::INFINITY, 0.0, 0.0),
            ..wrenches[0]
        },
        Wrench::ZERO,
    ];
    assert!(
        apply_rigid_adhesive_torsion_with_wrenches(
            &mut bodies,
            &f.properties,
            &mut state,
            contact,
            invalid,
            dt,
            controls(),
            &CancelGate::new(),
        )
        .is_err()
    );
    assert_eq!(bodies, f.bodies);
    assert_eq!(state, initial);
    let expected = prepare_rigid_adhesive_torsion_with_wrenches(
        &bodies,
        &f.properties,
        &state,
        contact,
        wrenches,
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let report = apply_rigid_adhesive_torsion_with_wrenches(
        &mut bodies,
        &f.properties,
        &mut state,
        contact,
        wrenches,
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(report, expected.report);
    assert_eq!(bodies, expected.bodies_after);
    assert_eq!(state, expected.state_after);
    assert!(report.contact.torsion().irreversible_loss_j > 0.0);
    close(
        report.kinetic_energy_change_j
            + report.contact.torsion().storage_change_j
            + report.contact.torsion().irreversible_loss_j,
        report.external_work_j,
        1e-25,
    );
}
