//! G1/G3/G4: source-owned adhesive rolling changes actual finite-body momenta.
#[path = "../../fs-contact/tests/source_point/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::rolling::{JkrRollingControls, SourceJkrRollingContact};
use fs_contact::normal_patch::adhesive::{JkrEnvelope, SourceJkrSphereContact};
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::{RigidContactBinding, adhesive::rolling::*};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_tribo::partial_slip::TangentFrame;

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
