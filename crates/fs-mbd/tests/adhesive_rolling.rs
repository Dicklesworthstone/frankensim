//! G1/G3/G4: source-owned adhesive rolling changes actual finite-body momenta.
#[path = "../../fs-contact/tests/source_point/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::rolling::{JkrRollingControls, SourceJkrRollingContact};
use fs_contact::normal_patch::adhesive::{JkrEnvelope, SourceJkrSphereContact};
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::{RigidContactBinding, adhesive::rolling::*};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3, Wrench};
use fs_tribo::partial_slip::TangentFrame;

#[path = "adhesive_rolling/events.rs"]
mod events;
#[path = "adhesive_rolling/normal.rs"]
mod normal;

fn source(tensile: bool) -> SourceJkrRollingContact {
    let a = common::solid("rolling-A", 1000.0, 1e7);
    let b = common::solid("rolling-B", 2000.0, 2e7);
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
    SourceJkrRollingContact::new(
        normal,
        approach,
        JkrRollingControls {
            critical_rolling_displacement_m: 1e-8,
            maximum_critical_displacement_to_contact_radius: 0.01,
        },
    )
    .unwrap()
}
fn controls() -> RigidAdhesiveRollingControls {
    RigidAdhesiveRollingControls {
        maximum_bodies: 4,
        maximum_step_s: 0.01,
        maximum_angle_increment_rad: 1e-4,
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
fn vec_close(a: Vec3, b: Vec3, absolute: f64) {
    close(a.x, b.x, absolute);
    close(a.y, b.y, absolute);
    close(a.z, b.z, absolute);
}

struct Fixture {
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    binding: RigidContactBinding,
    frame: TangentFrame,
}
impl Fixture {
    fn contact<'a>(&self, source: &'a SourceJkrRollingContact) -> RigidAdhesiveRollingContact<'a> {
        RigidAdhesiveRollingContact {
            binding: self.binding,
            frame: self.frame,
            source,
        }
    }
    fn mobility(&self) -> f64 {
        let axis = self.frame.longitudinal();
        let axis = Vec3::new(axis[0], axis[1], axis[2]);
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
        TangentFrame::new([-2.0, 1.0, 0.0], [1.0, 2.0, 3.0]).unwrap()
    } else {
        TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap()
    };
    let axes = frame.longitudinal();
    let axis = Vec3::new(axes[0], axes[1], axes[2]);
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
fn g1_elastic_rotor_pair_matches_independent_midpoint_solution() {
    let source = source(false);
    let dt = 1e-4;
    let rate = 0.25 * source.law().critical_rolling_angle_rad() / dt;
    let f = fixture(rate, 0.0, false);
    let contact = f.contact(&source);
    let state = contact
        .initial_state(&f.bodies, &f.properties, controls(), &CancelGate::new())
        .unwrap();
    let candidate = prepare_rigid_adhesive_rolling(
        &f.bodies,
        &f.properties,
        &state,
        contact,
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let r = source.normal_source().effective_radius_m();
    let stiffness = source.law().stiffness_n_per_m() * r * r;
    let expected_angle = dt * rate / (1.0 + 0.25 * dt * dt * f.mobility() * stiffness);
    let expected_moment = -0.5 * stiffness * expected_angle;
    assert_eq!(candidate.report.evaluations, 1);
    close(candidate.report.angle_increment_rad, expected_angle, 1e-18);
    close(
        candidate.report.contact.rolling().mean_moment_n_m,
        expected_moment,
        1e-23,
    );
    close(
        candidate.report.rolling_rate_after_rad_s,
        rate + dt * f.mobility() * expected_moment,
        1e-12,
    );
    assert_eq!(candidate.report.contact.rolling().irreversible_loss_j, 0.0);
    assert!(candidate.report.contact.rolling().storage_after_j > 0.0);
    close(
        candidate.report.kinetic_energy_change_j,
        -0.5 * stiffness * expected_angle.powi(2),
        1e-26,
    );
    for i in 0..2 {
        assert_eq!(candidate.bodies_after[i].pose(), f.bodies[i].pose());
        assert_eq!(
            candidate.bodies_after[i].linear_momentum_world(),
            f.bodies[i].linear_momentum_world()
        );
    }
    vec_close(
        candidate.report.angular_momentum_change_nms,
        Vec3::ZERO,
        1e-25,
    );
}

#[test]
fn g1_applied_torques_load_the_original_elastic_bond_with_full_inertia() {
    let source = source(true);
    let dt = 1e-4;
    let critical = source.law().critical_rolling_angle_rad();
    let stiffness = source.law().critical_moment_n_m() / critical;
    for rotated in [false, true] {
        let f = fixture(0.0, 0.0, rotated);
        let axis = f.frame.longitudinal();
        let axis = Vec3::new(axis[0], axis[1], axis[2]);
        let contact = f.contact(&source);
        let state = contact
            .initial_state(&f.bodies, &f.properties, controls(), &CancelGate::new())
            .unwrap();
        let alpha = 0.5 * critical / (dt * dt);
        let torques: [Vec3; 2] = core::array::from_fn(|i| {
            let desired = axis.scale(if i == 0 { alpha } else { -0.25 * alpha });
            f.properties[i]
                .angular_momentum_body_checked(
                    f.bodies[i]
                        .pose()
                        .orientation()
                        .rotate_world_to_body(desired),
                )
                .unwrap()
        });
        let wrenches = [
            Wrench {
                force_world: Vec3::new(1e-8, -2e-8, 3e-8),
                torque_body: torques[0],
            },
            Wrench {
                force_world: Vec3::new(-2e-8, 1e-8, -1e-8),
                torque_body: torques[1],
            },
        ];
        let candidate = prepare_rigid_adhesive_rolling_with_wrenches(
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
        let expected_angle =
            0.5 * dt * dt * (1.25 * alpha) / (1.0 + 0.25 * dt * dt * f.mobility() * stiffness);
        close(candidate.report.angle_increment_rad, expected_angle, 1e-18);
        close(
            candidate.report.rolling_rate_after_rad_s,
            2.0 * expected_angle / dt,
            1e-12,
        );
        assert_eq!(candidate.report.evaluations, 1);
        assert_eq!(candidate.report.contact.rolling().irreversible_loss_j, 0.0);
        let mut external_work = 0.0;
        for i in 0..2 {
            let after = candidate.bodies_after[i];
            let expected_angular = f.bodies[i]
                .pose()
                .orientation()
                .rotate_world_to_body(
                    candidate
                        .report
                        .angular_impulse_on_a_world_nms
                        .scale(if i == 0 { 1.0 } else { -1.0 }),
                )
                .add(torques[i].scale(dt));
            vec_close(after.angular_momentum_body(), expected_angular, 1e-25);
            vec_close(
                after.linear_momentum_world(),
                wrenches[i].force_world.scale(dt),
                1e-25,
            );
            let travel = after
                .center_of_mass_velocity_world(f.properties[i])
                .unwrap()
                .scale(0.5 * dt);
            let angular = f.properties[i]
                .angular_velocity_body_checked(after.angular_momentum_body())
                .unwrap()
                .scale(0.5 * dt);
            external_work += wrenches[i].force_world.dot(travel) + torques[i].dot(angular);
            vec_close(candidate.report.body_displacement_world_m[i], travel, 1e-20);
            vec_close(
                candidate.report.angular_displacement_body_rad[i],
                angular,
                1e-18,
            );
        }
        close(candidate.report.external_work_j, external_work, 1e-25);
        close(
            candidate.report.kinetic_energy_change_j
                + candidate.report.contact.rolling().storage_change_j,
            external_work,
            1e-25,
        );
        vec_close(
            candidate.report.linear_momentum_residual_ns,
            Vec3::ZERO,
            1e-25,
        );
        vec_close(
            candidate.report.angular_momentum_residual_nms,
            Vec3::ZERO,
            1e-25,
        );
    }
}

#[test]
fn g1_sustained_applied_torque_drives_rolling_yield_from_rest() {
    let source = source(true);
    let f = fixture(0.0, 0.0, false);
    let contact = f.contact(&source);
    let mut state = contact
        .initial_state(&f.bodies, &f.properties, controls(), &CancelGate::new())
        .unwrap();
    let mut bodies = f.bodies;
    let dt = 1e-3;
    let torque = 10.0 * source.law().critical_moment_n_m();
    let wrenches = [
        Wrench {
            force_world: Vec3::ZERO,
            torque_body: Vec3::new(torque, 0.0, 0.0),
        },
        Wrench::ZERO,
    ];
    let mut input = 0.0;
    let mut loss = 0.0;
    let mut storage = 0.0;
    for _ in 0..8 {
        let report = apply_rigid_adhesive_rolling_with_wrenches(
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
        assert!(report.rolling_rate_after_rad_s > report.rolling_rate_before_rad_s);
        input += report.external_work_j;
        loss += report.contact.rolling().irreversible_loss_j;
        storage = report.contact.rolling().storage_after_j;
    }
    assert!(loss > 0.0);
    close(
        state.elastic_shift_m(),
        source.law().parameters().critical_rolling_displacement_m,
        1e-22,
    );
    let kinetic: f64 = (0..2)
        .map(|i| bodies[i].kinetic_energy(f.properties[i]).unwrap())
        .sum();
    close(kinetic + storage + loss, input, 1e-24);
}

#[test]
fn g3_zero_load_equivalence_and_common_translation_do_not_change_adhesion() {
    let source = source(false);
    let f = fixture(1e-3, 0.0, false);
    let contact = f.contact(&source);
    let state = contact
        .initial_state(&f.bodies, &f.properties, controls(), &CancelGate::new())
        .unwrap();
    let dt = 1e-4;
    let original = prepare_rigid_adhesive_rolling(
        &f.bodies,
        &f.properties,
        &state,
        contact,
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let zero = prepare_rigid_adhesive_rolling_with_wrenches(
        &f.bodies,
        &f.properties,
        &state,
        contact,
        [Wrench::ZERO; 2],
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(original, zero);
    let acceleration = Vec3::new(0.0, 0.0, 0.02);
    let wrenches = f.properties.map(|p| Wrench {
        force_world: acceleration.scale(p.mass()),
        torque_body: Vec3::ZERO,
    });
    let driven = prepare_rigid_adhesive_rolling_with_wrenches(
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
    assert_eq!(driven.report.contact, original.report.contact);
    assert_eq!(driven.state_after, original.state_after);
    close(
        driven.report.external_work_j,
        0.5 * (f.properties[0].mass() + f.properties[1].mass())
            * acceleration.dot(acceleration)
            * dt
            * dt,
        1e-26,
    );
}

#[test]
fn g4_applied_load_refusals_are_atomic_and_retryable() {
    let source = source(false);
    let f = fixture(1e-3, 0.0, false);
    let contact = f.contact(&source);
    let original = contact
        .initial_state(&f.bodies, &f.properties, controls(), &CancelGate::new())
        .unwrap();
    let mut bodies = f.bodies;
    let mut state = original;
    let dt = 1e-4;
    let invalid = [
        Wrench {
            force_world: Vec3::ZERO,
            torque_body: Vec3::new(f64::NAN, 0.0, 0.0),
        },
        Wrench::ZERO,
    ];
    assert!(
        apply_rigid_adhesive_rolling_with_wrenches(
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
    assert_eq!(state, original);
    // A large negative applied torque turns this trajectory. Its net angle
    // cannot authorize an unrepresented out-and-back plastic loading path.
    let reversing = [
        Wrench {
            force_world: Vec3::ZERO,
            torque_body: Vec3::new(-1e-7, 0.0, 0.0),
        },
        Wrench::ZERO,
    ];
    assert!(matches!(
        apply_rigid_adhesive_rolling_with_wrenches(
            &mut bodies,
            &f.properties,
            &mut state,
            contact,
            reversing,
            dt,
            controls(),
            &CancelGate::new(),
        ),
        Err(RigidAdhesiveRollingError::UnresolvedRollingReversal)
    ));
    assert_eq!(bodies, f.bodies);
    assert_eq!(state, original);
    let valid = [
        Wrench {
            force_world: Vec3::ZERO,
            torque_body: Vec3::new(1e-10, 0.0, 0.0),
        },
        Wrench::ZERO,
    ];
    let expected = prepare_rigid_adhesive_rolling_with_wrenches(
        &bodies,
        &f.properties,
        &state,
        contact,
        valid,
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let report = apply_rigid_adhesive_rolling_with_wrenches(
        &mut bodies,
        &f.properties,
        &mut state,
        contact,
        valid,
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(report, expected.report);
    assert_eq!(bodies, expected.bodies_after);
    assert_eq!(state, expected.state_after);
}

#[test]
fn g1_yielded_tensile_contact_retains_storage_and_resolves_physical_loss() {
    let source = source(true);
    assert!(source.normal_response().normal_force_n < 0.0);
    let dt = 1e-4;
    let f = fixture(
        4.0 * source.law().critical_rolling_angle_rad() / dt,
        0.0,
        false,
    );
    let contact = f.contact(&source);
    let mut state = contact
        .initial_state(&f.bodies, &f.properties, controls(), &CancelGate::new())
        .unwrap();
    let mut bodies = f.bodies;
    let report = apply_rigid_adhesive_rolling(
        &mut bodies,
        &f.properties,
        &mut state,
        contact,
        dt,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let rolling = report.contact.rolling();
    assert!(rolling.irreversible_loss_j > 0.0 && rolling.storage_after_j > 0.0);
    close(
        rolling.endpoint_moment_n_m,
        -source.law().critical_moment_n_m(),
        1e-24,
    );
    close(
        -report.kinetic_energy_change_j,
        rolling.storage_change_j + rolling.irreversible_loss_j,
        1e-25,
    );
    close(
        state.elastic_shift_m(),
        source.law().parameters().critical_rolling_displacement_m,
        1e-22,
    );
    for i in 0..2 {
        let actual = bodies[i]
            .angular_momentum_body()
            .sub(f.bodies[i].angular_momentum_body());
        let expected = report
            .angular_impulse_on_a_world_nms
            .scale(if i == 0 { 1.0 } else { -1.0 });
        vec_close(actual, expected, 1e-25);
    }
}

#[test]
fn g3_arbitrary_body_axes_retain_full_inertia_and_conjugate_work() {
    let source = source(false);
    let f = fixture(0.003, 0.0, true);
    let contact = f.contact(&source);
    let state = contact
        .initial_state(&f.bodies, &f.properties, controls(), &CancelGate::new())
        .unwrap();
    let candidate = prepare_rigid_adhesive_rolling(
        &f.bodies,
        &f.properties,
        &state,
        contact,
        1e-4,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    close(
        candidate.report.rotational_mobility_per_kg_m2,
        f.mobility(),
        1e-5,
    );
    let axis = f.frame.longitudinal();
    let axis = Vec3::new(axis[0], axis[1], axis[2]);
    for i in 0..2 {
        let before = f.bodies[i];
        let expected = before.angular_momentum_body().add(
            before.pose().orientation().rotate_world_to_body(
                candidate
                    .report
                    .angular_impulse_on_a_world_nms
                    .scale(if i == 0 { 1.0 } else { -1.0 }),
            ),
        );
        vec_close(
            candidate.bodies_after[i].angular_momentum_body(),
            expected,
            1e-25,
        );
    }
    close(
        candidate.report.angular_impulse_on_a_world_nms.dot(axis),
        1e-4 * candidate.report.contact.rolling().mean_moment_n_m,
        1e-25,
    );
    close(candidate.report.total_energy_residual_j, 0.0, 1e-25);
}

#[test]
fn g3_common_translation_cannot_hide_small_rolling_loss() {
    let source = source(false);
    let run = |boost| {
        let f = fixture(0.04, boost, false);
        let contact = f.contact(&source);
        let state = contact
            .initial_state(&f.bodies, &f.properties, controls(), &CancelGate::new())
            .unwrap();
        prepare_rigid_adhesive_rolling(
            &f.bodies,
            &f.properties,
            &state,
            contact,
            1e-4,
            controls(),
            &CancelGate::new(),
        )
        .unwrap()
    };
    let rest = run(0.0);
    let moving = run(1e12);
    assert!(rest.report.contact.rolling().irreversible_loss_j > 0.0);
    assert_eq!(
        rest.report.kinetic_energy_change_j,
        moving.report.kinetic_energy_change_j
    );
    assert_eq!(rest.report.contact, moving.report.contact);
    assert!(moving.report.kinetic_energy_change_j < 0.0);
}

#[test]
fn g4_cancellation_bounds_and_changed_geometry_leave_body_and_history_unchanged() {
    let source = source(false);
    let f = fixture(0.04, 0.0, false);
    let contact = f.contact(&source);
    let initial = contact
        .initial_state(&f.bodies, &f.properties, controls(), &CancelGate::new())
        .unwrap();
    let mut state = initial;
    let mut bodies = f.bodies;
    let gate = CancelGate::new();
    gate.request();
    assert!(
        apply_rigid_adhesive_rolling(
            &mut bodies,
            &f.properties,
            &mut state,
            contact,
            1e-4,
            controls(),
            &gate
        )
        .is_err()
    );
    let restricted = RigidAdhesiveRollingControls {
        maximum_angle_increment_rad: 1e-8,
        ..controls()
    };
    assert!(matches!(
        apply_rigid_adhesive_rolling(
            &mut bodies,
            &f.properties,
            &mut state,
            contact,
            1e-4,
            restricted,
            &CancelGate::new()
        ),
        Err(RigidAdhesiveRollingError::StepLimit)
    ));
    let changed = RigidAdhesiveRollingContact {
        frame: TangentFrame::new([0.0, 0.0, 1.0], [0.0, 1.0, 0.0]).unwrap(),
        ..contact
    };
    assert!(matches!(
        apply_rigid_adhesive_rolling(
            &mut bodies,
            &f.properties,
            &mut state,
            changed,
            1e-4,
            controls(),
            &CancelGate::new()
        ),
        Err(RigidAdhesiveRollingError::GeometryMismatch)
    ));
    let short = RigidAdhesiveRollingControls {
        maximum_iterations: 1,
        angle_tolerance_rad: 1e-30,
        ..controls()
    };
    let unresolved = apply_rigid_adhesive_rolling(
        &mut bodies,
        &f.properties,
        &mut state,
        contact,
        1e-4,
        short,
        &CancelGate::new(),
    );
    assert!(
        matches!(
            unresolved,
            Err(RigidAdhesiveRollingError::Unresolved { evaluations: 3, .. })
        ),
        "{unresolved:?}"
    );
    assert_eq!(bodies, f.bodies);
    assert_eq!(state, initial);
}

#[test]
fn g3_elastic_looking_net_step_cannot_hide_yield_before_a_turning_point() {
    let source = source(false);
    let r = source.normal_source().effective_radius_m();
    let stiffness = source.law().stiffness_n_per_m() * r * r;
    let critical = source.law().critical_rolling_angle_rad();
    let reference = fixture(0.0, 0.0, false);
    let frequency = (stiffness * reference.mobility()).sqrt();
    // In units K = mobility = theta_c = 1, omega0=1.1 and dt=4
    // give an elastic midpoint net angle 0.88 and endpoint rate -0.66.
    // But the initial kinetic energy exceeds yield storage: the intervening
    // physical path reaches yield before turning and must not get zero loss.
    let dt = 4.0 / frequency;
    let f = fixture(1.1 * critical * frequency, 0.0, false);
    let contact = f.contact(&source);
    let initial = contact
        .initial_state(&f.bodies, &f.properties, controls(), &CancelGate::new())
        .unwrap();
    let trial = source
        .advance(&initial.source_state(), 0.88 * critical)
        .unwrap();
    assert_eq!(trial.rolling().irreversible_loss_j, 0.0);
    let mut bodies = f.bodies;
    let mut state = initial;
    let result = apply_rigid_adhesive_rolling(
        &mut bodies,
        &f.properties,
        &mut state,
        contact,
        dt,
        RigidAdhesiveRollingControls {
            energy_tolerance_j: 1.0,
            ..controls()
        },
        &CancelGate::new(),
    );
    assert!(
        matches!(
            result,
            Err(RigidAdhesiveRollingError::UnresolvedRollingReversal)
        ),
        "{result:?}"
    );
    assert_eq!(bodies, f.bodies);
    assert_eq!(state, initial);
}

#[test]
fn g1_guarded_elastic_orbit_admits_rocking_but_refuses_ambiguous_yield() {
    let source = source(false);
    let radius = source.normal_source().effective_radius_m();
    let stiffness = source.law().stiffness_n_per_m() * radius * radius;
    let critical = source.law().critical_rolling_angle_rad();
    for rotated in [false, true] {
        let frequency = (stiffness * fixture(0.0, 0.0, rotated).mobility()).sqrt();
        let dt = 4.0 / frequency;
        let policy = RigidAdhesiveRollingControls {
            maximum_step_s: 1.1 * dt,
            energy_tolerance_j: 1.0,
            ..controls()
        };
        for speed in [0.8, 1.0] {
            let f = fixture(speed * critical * frequency, 1e12, rotated);
            let contact = f.contact(&source);
            let state = contact
                .initial_state(&f.bodies, &f.properties, policy, &CancelGate::new())
                .unwrap();
            let result = prepare_rigid_adhesive_rolling(
                &f.bodies,
                &f.properties,
                &state,
                contact,
                dt,
                policy,
                &CancelGate::new(),
            );
            if speed == 1.0 {
                assert!(
                    matches!(
                        result,
                        Err(RigidAdhesiveRollingError::UnresolvedRollingReversal)
                    ),
                    "{result:?}"
                );
                continue;
            }
            let candidate = result.unwrap();
            let report = &candidate.report;
            assert!(
                report.rolling_rate_before_rad_s > 0.0 && report.rolling_rate_after_rad_s < 0.0
            );
            close(report.angle_increment_rad, 0.64 * critical, 1e-18);
            close(
                report.rolling_rate_after_rad_s,
                -0.48 * critical * frequency,
                1e-13,
            );
            assert_eq!(report.contact.rolling().irreversible_loss_j, 0.0);
            assert_eq!(report.evaluations, 1);
            close(report.total_energy_residual_j, 0.0, 1e-25);
            close(
                report.contact.rolling().storage_after_j,
                0.5 * stiffness * (0.64 * critical).powi(2),
                1e-25,
            );
        }
    }
}
