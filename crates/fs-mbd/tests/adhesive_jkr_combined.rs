//! One real disk, four source coordinates, and independent actual-body balances.
#[path = "../../fs-contact/tests/adhesive_dmt_coupled/common.rs"]
mod common;

#[path = "adhesive_jkr_combined/temperature.rs"]
mod temperature;

use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::{
    JkrEnvelope, SourceJkrSphereContact,
    mixed_mode::{
        JkrMixedModeEnvelope, SourceJkrMixedMode, combined::motion::JkrCombinedMotionBudget,
    },
};
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::{RigidContactBinding, adhesive::jkr_combined::*};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3, Wrench};
use fs_tribo::partial_slip::TangentFrame;

const MASS: f64 = 1e-5;
const ARM: f64 = 0.01;
const INERTIA: f64 = MASS * ARM * ARM;
const FORCE: f64 = PI * 0.05 * 0.01 / 6.0;

fn source() -> SourceJkrMixedMode {
    let (original, _) = common::sources(0.05, 1e9, 5e8);
    let normal = SourceJkrSphereContact::new(
        original.normal_state().clone(),
        0.01,
        1e-12,
        JkrEnvelope {
            maximum_indentation_m: 1e-7,
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
    SourceJkrMixedMode::new(
        normal,
        JkrMixedModeEnvelope {
            maximum_tangential_displacement_to_radius: 0.1,
            maximum_mean_shear_stress_to_yield: 0.1,
        },
    )
    .unwrap()
}

fn scales(source: &SourceJkrMixedMode) -> (f64, f64, f64, f64, f64) {
    let length = source.normal_source().indentation_scale_m();
    let speed = (FORCE * length / (1.2 * MASS)).sqrt();
    let stiffness = source
        .evaluate_combined_displacement_controlled([0.0; 4])
        .unwrap()
        .coupled_stiffness
        .unwrap();
    (
        length,
        speed,
        length / speed,
        (FORCE * length / stiffness[1][1]).sqrt(),
        (FORCE * length / stiffness[3][3]).sqrt(),
    )
}

fn controls(source: &SourceJkrMixedMode) -> RigidJkrCombinedControls {
    let (length, speed, time, shear, twist) = scales(source);
    let torque = FORCE * length / twist;
    RigidJkrCombinedControls {
        motion: JkrCombinedMotionBudget {
            maximum_step_s: 2.0 * time,
            maximum_coordinate_increment: [20.0 * length, 20.0 * shear, 20.0 * shear, 20.0 * twist],
            maximum_rate: [
                100.0 * speed,
                100.0 * shear / time,
                100.0 * shear / time,
                100.0 * twist / time,
            ],
            maximum_effort: [1e4 * FORCE, 1e4 * FORCE, 1e4 * FORCE, 1e4 * torque],
            maximum_iterations: 64,
            maximum_evaluations: 8192,
            maximum_backtracks: 24,
            coordinate_absolute_tolerance: [
                1e-13 * length,
                1e-13 * shear,
                1e-13 * shear,
                1e-13 * twist,
            ],
            coordinate_relative_tolerance: 2e-12,
            rate_absolute_tolerance: [
                1e-12 * speed,
                1e-12 * shear / time,
                1e-12 * shear / time,
                1e-12 * twist / time,
            ],
            rate_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-10 * FORCE * length,
            energy_relative_tolerance: 1e-10,
        },
        maximum_bodies: 4,
        geometry_tolerance_m: 2e-17,
        displacement_tolerance_m: 1e-9 * length.max(shear),
        twist_tolerance_rad: 1e-9 * twist,
        velocity_tolerance_m_s: 1e-9 * speed.max(shear / time),
        angular_velocity_tolerance_rad_s: 1e-9 * twist / time,
        linear_momentum_tolerance_ns: 1e-11 * MASS * speed,
        angular_momentum_tolerance_nms: 1e-11 * MASS * speed * ARM,
        energy_tolerance_j: 1e-8 * FORCE * length,
    }
}

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "actual={a:e}, expected={b:e}, tolerance={tolerance:e}"
    );
}
fn close_vec(a: Vec3, b: Vec3, tolerance: f64) {
    close(a.x, b.x, tolerance);
    close(a.y, b.y, tolerance);
    close(a.z, b.z, tolerance);
}
fn vector(v: [f64; 3]) -> Vec3 {
    Vec3::new(v[0], v[1], v[2])
}

struct Fixture {
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    arms: [Vec3; 2],
    frame: TangentFrame,
}
impl Fixture {
    fn contact<'a>(&self, source: &'a SourceJkrMixedMode) -> RigidJkrCombinedContact<'a> {
        RigidJkrCombinedContact {
            binding: RigidContactBinding {
                body_a: 0,
                body_b: 1,
                arm_a_body: self.arms[0],
                arm_b_body: self.arms[1],
            },
            frame: self.frame,
            reference_approach_m: 0.0,
            source,
        }
    }
    fn prepared<'a>(&self, source: &'a SourceJkrMixedMode) -> PreparedRigidJkrCombinedPair<'a> {
        PreparedRigidJkrCombinedPair::new(
            &self.bodies,
            &self.properties,
            self.contact(source),
            controls(source),
            &CancelGate::new(),
        )
        .unwrap()
    }
    fn set_rates(&mut self, i: usize, velocity: Vec3, omega_body: Vec3) {
        self.bodies[i] = RigidBodyState::new(
            self.bodies[i].pose(),
            velocity.scale(self.properties[i].mass()),
            self.properties[i]
                .angular_momentum_body_checked(omega_body)
                .unwrap(),
        )
        .unwrap();
    }
}

fn fixture(arms: [Vec3; 2]) -> Fixture {
    Fixture {
        bodies: arms.map(|arm| {
            RigidBodyState::new(
                Pose::new(arm.scale(-1.0), UnitQuaternion::IDENTITY).unwrap(),
                Vec3::ZERO,
                Vec3::ZERO,
            )
            .unwrap()
        }),
        properties: [2.0, 3.0].map(|m| {
            MassProperties::new(
                m * MASS,
                Vec3::ZERO,
                Vec3::new(0.8, 1.1, 1.4).scale(INERTIA),
            )
            .unwrap()
        }),
        arms,
        frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
    }
}

fn facing() -> Fixture {
    // Equal inertias and opposite facing arms turn shear reactions into common
    // transverse spin. Relative rolling remains zero without suppressing torque.
    fixture([Vec3::new(0.0, 0.0, ARM), Vec3::new(0.0, 0.0, -ARM)])
}

fn check_native(f: &Fixture, source: &SourceJkrMixedMode, c: &RigidJkrCombinedCandidate) {
    let limits = controls(source);
    let dt = c.state_after.time_s() - c.report.motion.initial_state().time_s();
    let reaction = f.prepared(source).reaction_arms_body();
    let mut kinetic = 0.0;
    let mut external_work = 0.0;
    let mut linear = Vec3::ZERO;
    let mut angular = Vec3::ZERO;
    let mut expected_linear = Vec3::ZERO;
    let mut expected_angular = Vec3::ZERO;
    for (i, sign) in [(0, 1.0), (1, -1.0)] {
        let before = f.bodies[i];
        let after = c.bodies_after[i];
        let pose = before.pose();
        let impulse = c.report.force_impulse_on_a_world_ns.scale(sign);
        let couple = c.report.angular_impulse_on_a_world_nms.scale(sign);
        let external = c.report.external_wrenches[i];
        let dp = after
            .linear_momentum_world()
            .sub(before.linear_momentum_world());
        let dl = after
            .angular_momentum_body()
            .sub(before.angular_momentum_body());
        close_vec(
            dp,
            impulse.add(external.force_world.scale(dt)),
            limits.linear_momentum_tolerance_ns,
        );
        close_vec(
            dl,
            reaction[i]
                .cross(pose.orientation().rotate_world_to_body(impulse))
                .add(pose.orientation().rotate_world_to_body(couple))
                .add(external.torque_body.scale(dt)),
            limits.angular_momentum_tolerance_nms,
        );
        assert_eq!(before.pose(), after.pose());
        kinetic += before
            .kinetic_energy_change_to(f.properties[i], after)
            .unwrap();
        let k0 = before
            .point_kinematics(f.properties[i], Vec3::ZERO)
            .unwrap();
        let k1 = after.point_kinematics(f.properties[i], Vec3::ZERO).unwrap();
        let travel = k0
            .center_of_mass_velocity_world
            .add(k1.center_of_mass_velocity_world)
            .scale(0.5 * dt);
        let angle = k0
            .angular_velocity_body
            .add(k1.angular_velocity_body)
            .scale(0.5 * dt);
        close_vec(
            c.report.body_displacement_world_m[i],
            travel,
            limits.displacement_tolerance_m,
        );
        close_vec(
            c.report.angular_displacement_body_rad[i],
            angle,
            limits.twist_tolerance_rad,
        );
        external_work += external.force_world.dot(travel) + external.torque_body.dot(angle);
        linear = linear.add(dp);
        angular = angular
            .add(pose.position_world().cross(dp))
            .add(pose.orientation().rotate_body_to_world(dl));
        expected_linear = expected_linear.add(external.force_world.scale(dt));
        expected_angular = expected_angular
            .add(pose.position_world().cross(external.force_world.scale(dt)))
            .add(
                pose.orientation()
                    .rotate_body_to_world(external.torque_body.scale(dt)),
            );
    }
    close(
        kinetic,
        c.report.kinetic_energy_change_j,
        limits.energy_tolerance_j,
    );
    close(
        external_work,
        c.report.external_work_j,
        limits.energy_tolerance_j,
    );
    close(
        kinetic + c.report.contact_storage_change_j - external_work,
        0.0,
        limits.energy_tolerance_j,
    );
    close(
        c.report.contact_impulse_work_j + c.report.contact_storage_change_j,
        0.0,
        limits.energy_tolerance_j,
    );
    close_vec(linear, expected_linear, limits.linear_momentum_tolerance_ns);
    close_vec(
        angular,
        expected_angular,
        limits.angular_momentum_tolerance_nms,
    );
    let end = c.state_after.coordinates();
    let actual = source
        .evaluate_combined_displacement_controlled(end)
        .unwrap();
    assert_eq!(actual.model_identity(), source.identity());
    close(
        actual.reversible_energy_j - c.report.motion.initial_contact.reversible_energy_j,
        c.report.contact_storage_change_j,
        limits.energy_tolerance_j,
    );
    for i in 0..4 {
        close(
            actual.force[i],
            c.report.motion.contact.force[i],
            1e-11 * actual.force[i].abs().max(1e-25),
        );
    }
}

#[test]
fn simultaneous_shear_twist_and_normal_recoil_has_one_area_and_actual_external_work() {
    let source = source();
    let (_, speed, time, shear, twist) = scales(&source);
    let torque = FORCE * scales(&source).0 / twist;
    let mut f = facing();
    f.set_rates(
        0,
        Vec3::new(0.02 * shear / time, -0.01 * shear / time, 0.03 * speed),
        Vec3::new(0.0, 0.0, 0.02 * twist / time),
    );
    let pair = f.prepared(&source);
    let state = pair
        .initialize(
            [0.18 * shear, -0.12 * shear, 0.2 * twist],
            &CancelGate::new(),
        )
        .unwrap();
    let wrenches = [
        Wrench {
            force_world: Vec3::new(0.06 * FORCE, -0.03 * FORCE, 0.25 * FORCE),
            torque_body: Vec3::new(0.0, 0.0, 0.3 * torque),
        },
        Wrench {
            force_world: Vec3::new(-0.02 * FORCE, 0.01 * FORCE, -0.1 * FORCE),
            torque_body: Vec3::new(0.0, 0.0, -0.1 * torque),
        },
    ];
    let c = pair
        .advance_with_wrenches(state, wrenches, 0.03 * time, &CancelGate::new())
        .unwrap();
    check_native(&f, &source, &c);
    assert_ne!(c.report.external_work_j, 0.0);
    assert!(c.report.coordinate_increment.iter().all(|v| v.abs() > 0.0));
    let end = c.state_after.coordinates();
    let normal = source.normal_source().evaluate_attached(end[0]).unwrap();
    let torsion = source
        .evaluate_torsion_displacement_controlled(end[0], end[3])
        .unwrap();
    assert!(c.report.motion.contact.contact_area_m2 < torsion.contact_area_m2);
    assert!(torsion.contact_area_m2 < normal.contact_area_m2);
    assert_ne!(
        c.report.motion.contact.contact_area_m2,
        c.report.motion.initial_contact.contact_area_m2
    );
    let k = c.bodies_after.map(|body| body.angular_momentum_body());
    assert!(k[0].x.abs() + k[0].y.abs() > 0.0);
    close(
        k[0].x,
        k[1].x,
        controls(&source).angular_momentum_tolerance_nms,
    );
    close(
        k[0].y,
        k[1].y,
        controls(&source).angular_momentum_tolerance_nms,
    );
    assert_eq!(
        c.state_after.pending_com_travel_m(),
        c.report.body_displacement_world_m
    );
    assert_eq!(
        c.state_after.pending_angular_displacement_body_rad(),
        c.report.angular_displacement_body_rad
    );
}

#[test]
fn full_inertia_and_all_six_cross_mobilities_match_actual_impulses_and_held_preload() {
    let source = source();
    let mut f = fixture([
        Vec3::new(ARM, 0.3 * ARM, 0.2 * ARM),
        Vec3::new(0.4 * ARM, -0.2 * ARM, -0.5 * ARM),
    ]);
    for i in 0..2 {
        f.properties[i] = MassProperties::from_principal_frame(
            f.properties[i].mass(),
            Vec3::ZERO,
            f.properties[i].principal_inertia_body(),
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -0.7), 0.45 + 0.1 * i as f64)
                .unwrap(),
        )
        .unwrap();
    }
    let pair = f.prepared(&source);
    let w = pair.mobility();
    let axes = [f.frame.normal(), f.frame.longitudinal(), f.frame.lateral()].map(vector);
    let reaction = pair.reaction_arms_body();
    for row in 0..4 {
        for column in 0..row {
            assert_eq!(w[row][column], w[column][row]);
            assert!(w[row][column].abs() > 1e-5 * (w[row][row] * w[column][column]).sqrt());
        }
    }
    // Independent unit impulses use the existing body owner, not this module's
    // Gram construction, so this checks all point/couple signs and cross terms.
    for column in 0..4 {
        let mut endpoint = f.bodies;
        for (i, sign) in [(0, 1.0), (1, -1.0)] {
            let force = if column < 3 {
                axes[column].scale(sign)
            } else {
                Vec3::ZERO
            };
            let couple = if column == 3 {
                axes[0].scale(sign)
            } else {
                Vec3::ZERO
            };
            let orientation = endpoint[i].pose().orientation();
            let angular = reaction[i]
                .cross(orientation.rotate_world_to_body(force))
                .add(orientation.rotate_world_to_body(couple));
            endpoint[i] = endpoint[i]
                .apply_momentum_impulse(f.properties[i], force, angular)
                .unwrap()
                .state_after;
        }
        let a = endpoint[0]
            .point_kinematics(f.properties[0], reaction[0])
            .unwrap();
        let b = endpoint[1]
            .point_kinematics(f.properties[1], reaction[1])
            .unwrap();
        let velocity = a.relative_velocity_world(b).unwrap();
        let spin = a.angular_velocity_world.sub(b.angular_velocity_world);
        let projected = [
            velocity.dot(axes[0]),
            velocity.dot(axes[1]),
            velocity.dot(axes[2]),
            spin.dot(axes[0]),
        ];
        for row in 0..4 {
            close(
                projected[row],
                w[row][column],
                2e-12 * (w[row][row] * w[column][column]).sqrt(),
            );
        }
    }
    let (_, _, time, shear, twist) = scales(&source);
    let old = pair
        .initialize(
            [0.16 * shear, -0.1 * shear, 0.15 * twist],
            &CancelGate::new(),
        )
        .unwrap();
    let response = old.motion_state().contact();
    let force = axes[0]
        .scale(response.force[0])
        .add(axes[1].scale(response.force[1]))
        .add(axes[2].scale(response.force[2]));
    let wrenches = core::array::from_fn(|i| {
        let sign = if i == 0 { 1.0 } else { -1.0 };
        let orientation = f.bodies[i].pose().orientation();
        Wrench {
            force_world: force.scale(sign),
            torque_body: reaction[i]
                .cross(orientation.rotate_world_to_body(force.scale(sign)))
                .add(orientation.rotate_world_to_body(axes[0].scale(sign * response.force[3]))),
        }
    });
    let c = pair
        .advance_with_wrenches(old, wrenches, 0.02 * time, &CancelGate::new())
        .unwrap();
    check_native(&f, &source, &c);
    for i in 0..4 {
        let tolerance = if i == 3 {
            controls(&source).twist_tolerance_rad
        } else {
            controls(&source).displacement_tolerance_m
        };
        close(
            c.state_after.coordinates()[i],
            old.coordinates()[i],
            tolerance,
        );
        close(
            c.report.driving_acceleration[i],
            w[i].iter().zip(response.force).map(|(a, b)| a * b).sum(),
            1e-11 * c.report.driving_acceleration[i].abs().max(1e-20),
        );
    }
    close(
        c.report.kinetic_energy_change_j,
        0.0,
        controls(&source).energy_tolerance_j,
    );
}

#[test]
fn elastic_assisting_work_and_world_rotation_preserve_actual_reactions() {
    let source = source();
    let (_, _, time, shear, twist) = scales(&source);
    let f = facing();
    let pair = f.prepared(&source);
    let old = pair
        .initialize(
            [0.25 * shear, -0.18 * shear, 0.25 * twist],
            &CancelGate::new(),
        )
        .unwrap();
    let c = pair.advance(old, 0.025 * time, &CancelGate::new()).unwrap();
    check_native(&f, &source, &c);
    assert!(c.report.kinetic_energy_change_j > 0.0);
    assert!(c.report.contact_storage_change_j < 0.0);
    assert!(c.report.contact_impulse_work_j > 0.0);
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -1.0), 0.73).unwrap();
    let mut rotated = facing();
    for i in 0..2 {
        rotated.bodies[i] = RigidBodyState::new(
            Pose::new(
                rotation.rotate_body_to_world(f.bodies[i].pose().position_world()),
                rotation,
            )
            .unwrap(),
            Vec3::ZERO,
            Vec3::ZERO,
        )
        .unwrap();
    }
    let n = rotation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0));
    let t = rotation.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
    rotated.frame = TangentFrame::new([n.x, n.y, n.z], [t.x, t.y, t.z]).unwrap();
    let pair_rotated = rotated.prepared(&source);
    let old_rotated = pair_rotated
        .initialize(
            [0.25 * shear, -0.18 * shear, 0.25 * twist],
            &CancelGate::new(),
        )
        .unwrap();
    let rotated_c = pair_rotated
        .advance(old_rotated, 0.025 * time, &CancelGate::new())
        .unwrap();
    check_native(&rotated, &source, &rotated_c);
    close_vec(
        rotated_c.report.force_impulse_on_a_world_ns,
        rotation.rotate_body_to_world(c.report.force_impulse_on_a_world_ns),
        controls(&source).linear_momentum_tolerance_ns,
    );
    close_vec(
        rotated_c.report.angular_impulse_on_a_world_nms,
        rotation.rotate_body_to_world(c.report.angular_impulse_on_a_world_nms),
        controls(&source).angular_momentum_tolerance_nms,
    );
    close(
        rotated_c.report.contact_storage_change_j,
        c.report.contact_storage_change_j,
        controls(&source).energy_tolerance_j,
    );
    for i in 0..4 {
        close(
            rotated_c.state_after.coordinates()[i],
            c.state_after.coordinates()[i],
            if i == 3 {
                controls(&source).twist_tolerance_rad
            } else {
                controls(&source).displacement_tolerance_m
            },
        );
    }
}

#[test]
fn zero_shear_reduces_to_original_torsion_and_zero_twist_to_pure_normal() {
    let source = source();
    let (_, _, time, _, twist) = scales(&source);
    let f = facing();
    let pair = f.prepared(&source);
    for preload in [0.0, 0.2 * twist] {
        let old = pair
            .initialize([0.0, 0.0, preload], &CancelGate::new())
            .unwrap();
        let c = pair.advance(old, 0.02 * time, &CancelGate::new()).unwrap();
        check_native(&f, &source, &c);
        let q = c.state_after.coordinates();
        let torsion = source
            .evaluate_torsion_displacement_controlled(q[0], q[3])
            .unwrap();
        close(
            c.report.motion.contact.contact_radius_m,
            torsion.contact_radius_m,
            1e-11 * torsion.contact_radius_m,
        );
        close(
            c.report.motion.contact.force[0],
            torsion.normal_force_n,
            1e-11 * FORCE,
        );
        close(
            c.report.motion.contact.force[3],
            torsion.torsional_moment_n_m,
            1e-11 * FORCE * scales(&source).0 / twist,
        );
        assert_eq!([q[1], q[2]], [0.0; 2]);
        assert_eq!(
            [c.state_after.rates()[1], c.state_after.rates()[2]],
            [0.0; 2]
        );
        assert_eq!(
            [
                c.report.motion.mean_effort[1],
                c.report.motion.mean_effort[2]
            ],
            [0.0; 2]
        );
        if preload == 0.0 {
            assert_eq!(q[3], 0.0);
            assert_eq!(c.report.angular_impulse_on_a_world_nms, Vec3::ZERO);
            let normal = source.normal_source().evaluate_attached(q[0]).unwrap();
            close(
                c.report.motion.contact.contact_area_m2,
                normal.contact_area_m2,
                1e-11 * normal.contact_area_m2,
            );
        }
    }
}

#[test]
fn changed_history_stale_geometry_rolling_and_cancellation_never_publish_partial_bodies() {
    let source = source();
    let (_, speed, time, shear, twist) = scales(&source);
    let f = facing();
    let pair = f.prepared(&source);
    let old = pair
        .initialize(
            [0.2 * shear, -0.15 * shear, 0.15 * twist],
            &CancelGate::new(),
        )
        .unwrap();
    let expected = pair.advance(old, 0.02 * time, &CancelGate::new()).unwrap();
    let mut bodies = f.bodies;
    let mut state = old;
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        apply_rigid_jkr_combined(
            &mut bodies,
            &f.properties,
            f.contact(&source),
            &mut state,
            [Wrench::ZERO; 2],
            0.02 * time,
            controls(&source),
            &cancelled
        )
        .is_err()
    );
    assert_eq!(bodies, f.bodies);
    assert_eq!(state, old);
    let mut envelope = source.envelope();
    envelope.maximum_tangential_displacement_to_radius *= 0.9;
    let other = SourceJkrMixedMode::new(source.normal_source().clone(), envelope).unwrap();
    assert!(
        apply_rigid_jkr_combined(
            &mut bodies,
            &f.properties,
            f.contact(&other),
            &mut state,
            [Wrench::ZERO; 2],
            0.02 * time,
            controls(&source),
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(bodies, f.bodies);
    assert_eq!(state, old);
    let report = apply_rigid_jkr_combined(
        &mut bodies,
        &f.properties,
        f.contact(&source),
        &mut state,
        [Wrench::ZERO; 2],
        0.02 * time,
        controls(&source),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(report, expected.report);
    assert_eq!(bodies, expected.bodies_after);
    assert_eq!(state, expected.state_after);
    // A normal source step does not silently move the held material markers.
    let stale_bodies = bodies;
    let stale_state = state;
    assert!(matches!(
        apply_rigid_jkr_combined(
            &mut bodies,
            &f.properties,
            f.contact(&source),
            &mut state,
            [Wrench::ZERO; 2],
            0.02 * time,
            controls(&source),
            &CancelGate::new()
        ),
        Err(RigidJkrCombinedError::GeometryMismatch { .. })
    ));
    assert_eq!(bodies, stale_bodies);
    assert_eq!(state, stale_state);

    let mut rolling = facing();
    rolling.set_rates(0, Vec3::ZERO, Vec3::new(speed / ARM, 0.0, 0.0));
    assert!(matches!(
        rolling
            .prepared(&source)
            .initialize([0.0; 3], &CancelGate::new()),
        Err(RigidJkrCombinedError::UnsupportedRolling { .. })
    ));
    let mut asymmetric = facing();
    asymmetric.properties[1] = MassProperties::new(
        3.0 * MASS,
        Vec3::ZERO,
        Vec3::new(1.3, 1.7, 2.0).scale(INERTIA),
    )
    .unwrap();
    let old_asymmetric = asymmetric
        .prepared(&source)
        .initialize([0.2 * shear, -0.15 * shear, 0.0], &CancelGate::new())
        .unwrap();
    let mut asymmetric_state = old_asymmetric;
    let mut asymmetric_bodies = asymmetric.bodies;
    assert!(matches!(
        apply_rigid_jkr_combined(
            &mut asymmetric_bodies,
            &asymmetric.properties,
            asymmetric.contact(&source),
            &mut asymmetric_state,
            [Wrench::ZERO; 2],
            0.02 * time,
            controls(&source),
            &CancelGate::new()
        ),
        Err(RigidJkrCombinedError::UnsupportedRolling { .. })
    ));
    assert_eq!(asymmetric_bodies, asymmetric.bodies);
    assert_eq!(asymmetric_state, old_asymmetric);
    let mut changed = f.bodies;
    changed[0] = RigidBodyState::new(
        changed[0].pose(),
        Vec3::new(MASS * speed, 0.0, 0.0),
        Vec3::ZERO,
    )
    .unwrap();
    assert!(matches!(
        PreparedRigidJkrCombinedPair::new(
            &changed,
            &f.properties,
            f.contact(&source),
            controls(&source),
            &CancelGate::new()
        )
        .unwrap()
        .advance(old, 0.02 * time, &CancelGate::new()),
        Err(RigidJkrCombinedError::StateMismatch)
    ));
}
