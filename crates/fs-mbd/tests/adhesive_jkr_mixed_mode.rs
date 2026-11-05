//! Original mixed-mode area/work and actual full-inertia paired-body recoil.
#[path = "../../fs-contact/tests/adhesive_dmt_coupled/common.rs"]
mod common;
#[path = "adhesive_jkr_mixed_mode/geometry.rs"]
mod geometry;
#[path = "adhesive_jkr_mixed_mode/temperature.rs"]
mod temperature;

use core::f64::consts::PI;
use fs_contact::interface_binding::adhesive_shear::SourceAdhesiveShear;
use fs_contact::normal_patch::adhesive::{
    JkrEnvelope, SourceJkrSphereContact,
    mixed_mode::{JkrMixedModeEnvelope, SourceJkrMixedMode, motion::JkrMixedModeMotionBudget},
    motion::{JkrNormalMotion, JkrNormalMotionBudget, MAX_JKR_NORMAL_EVENTS},
};
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::{RigidContactBinding, adhesive::jkr_mixed_mode::*};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3, Wrench};
use fs_tribo::partial_slip::TangentFrame;

const MASS: f64 = 1e-5;
const ARM: f64 = 0.01;
const INERTIA: f64 = MASS * ARM * ARM;
const FORCE: f64 = PI * 0.05 * 0.01 / 6.0;

fn sources(strength_pa: f64) -> (SourceJkrMixedMode, SourceAdhesiveShear) {
    let (original, shear) = common::sources(0.05, strength_pa, 0.5 * strength_pa);
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
    (
        SourceJkrMixedMode::new(
            normal,
            JkrMixedModeEnvelope {
                maximum_tangential_displacement_to_radius: 0.1,
                maximum_mean_shear_stress_to_yield: 0.1,
            },
        )
        .unwrap(),
        shear,
    )
}
fn scales(source: &SourceJkrMixedMode) -> (f64, f64, f64) {
    let length = source.normal_source().indentation_scale_m();
    let speed = (FORCE * length / (1.2 * MASS)).sqrt();
    (length, speed, length / speed)
}
fn controls(source: &SourceJkrMixedMode) -> RigidJkrMixedModeControls {
    let (length, speed, time) = scales(source);
    RigidJkrMixedModeControls {
        motion: JkrMixedModeMotionBudget {
            maximum_step_s: 2.0 * time,
            maximum_coordinate_increment_m: 20.0 * length,
            maximum_speed_m_s: 100.0 * speed,
            maximum_force_n: 1e4 * FORCE,
            maximum_iterations: 64,
            maximum_evaluations: 2048,
            maximum_backtracks: 24,
            coordinate_absolute_tolerance_m: 1e-13 * length,
            coordinate_relative_tolerance: 2e-12,
            velocity_absolute_tolerance_m_s: 1e-12 * speed,
            velocity_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-10 * FORCE * length,
            energy_relative_tolerance: 1e-10,
        },
        maximum_bodies: 4,
        geometry_tolerance_m: 2e-18,
        displacement_tolerance_m: 1e-9 * length,
        velocity_tolerance_m_s: 1e-9 * speed,
        angular_velocity_tolerance_rad_s: 1e-9 * speed / ARM,
        linear_momentum_tolerance_ns: 1e-12 * MASS * speed,
        angular_momentum_tolerance_nms: 1e-12 * MASS * speed * ARM,
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
    fn contact<'a>(
        &self,
        source: &'a SourceJkrMixedMode,
        shear: &'a SourceAdhesiveShear,
    ) -> RigidJkrMixedModeContact<'a> {
        RigidJkrMixedModeContact {
            binding: RigidContactBinding {
                body_a: 0,
                body_b: 1,
                arm_a_body: self.arms[0],
                arm_b_body: self.arms[1],
            },
            frame: self.frame,
            reference_approach_m: 0.0,
            source,
            shear,
        }
    }
    fn prepared<'a>(
        &self,
        source: &'a SourceJkrMixedMode,
        shear: &'a SourceAdhesiveShear,
    ) -> PreparedRigidJkrMixedModePair<'a> {
        PreparedRigidJkrMixedModePair::new(
            &self.bodies,
            &self.properties,
            self.contact(source, shear),
            controls(source),
            &CancelGate::new(),
        )
        .unwrap()
    }
    fn set_velocity(&mut self, i: usize, v: Vec3) {
        self.bodies[i] = RigidBodyState::new(
            self.bodies[i].pose(),
            v.scale(self.properties[i].mass()),
            Vec3::ZERO,
        )
        .unwrap();
    }
}
fn fixture(arms: [Vec3; 2]) -> Fixture {
    let properties = [
        MassProperties::new(
            2.0 * MASS,
            Vec3::ZERO,
            Vec3::new(0.8, 1.1, 1.4).scale(INERTIA),
        )
        .unwrap(),
        MassProperties::new(
            3.0 * MASS,
            Vec3::ZERO,
            Vec3::new(1.5, 1.7, 2.0).scale(INERTIA),
        )
        .unwrap(),
    ];
    Fixture {
        bodies: arms.map(|arm| {
            RigidBodyState::new(
                Pose::new(arm.scale(-1.0), UnitQuaternion::IDENTITY).unwrap(),
                Vec3::ZERO,
                Vec3::ZERO,
            )
            .unwrap()
        }),
        properties,
        arms,
        frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
    }
}
fn axial() -> Fixture {
    fixture([Vec3::new(0.0, 0.0, ARM), Vec3::new(0.0, 0.0, -ARM)])
}
fn check_native(f: &Fixture, source: &SourceJkrMixedMode, c: &RigidJkrMixedModeCandidate) {
    let limits = controls(source);
    let dt = c.state_after.time_s() - c.report.motion.initial_state().time_s();
    let mut change = 0.0;
    let mut external_work = 0.0;
    let mut linear = Vec3::ZERO;
    let mut angular = Vec3::ZERO;
    let mut expected_linear = Vec3::ZERO;
    let mut expected_angular = Vec3::ZERO;
    for (i, sign) in [(0, 1.0), (1, -1.0)] {
        let before = f.bodies[i];
        let after = c.bodies_after[i];
        let pose = before.pose();
        let contact = c.report.force_impulse_on_a_world_ns.scale(sign);
        let external = c.report.external_wrenches[i];
        let dp = after
            .linear_momentum_world()
            .sub(before.linear_momentum_world());
        let dl = after
            .angular_momentum_body()
            .sub(before.angular_momentum_body());
        close_vec(
            dp,
            contact.add(external.force_world.scale(dt)),
            limits.linear_momentum_tolerance_ns,
        );
        close_vec(
            dl,
            f.arms[i]
                .cross(pose.orientation().rotate_world_to_body(contact))
                .add(external.torque_body.scale(dt)),
            limits.angular_momentum_tolerance_nms,
        );
        assert_eq!(before.pose(), after.pose());
        change += after.kinetic_energy(f.properties[i]).unwrap()
            - before.kinetic_energy(f.properties[i]).unwrap();
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
            limits.displacement_tolerance_m / ARM,
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
        change,
        c.report.kinetic_energy_change_j,
        limits.energy_tolerance_j,
    );
    close(
        external_work,
        c.report.external_work_j,
        limits.energy_tolerance_j,
    );
    let q0 = c.report.motion.initial_state().position_m();
    let q1 = c.state_after.position_m();
    let u0 = source
        .evaluate_displacement_controlled(q0[0], [q0[1], q0[2]])
        .unwrap()
        .reversible_energy_j;
    let u1 = source
        .evaluate_displacement_controlled(q1[0], [q1[1], q1[2]])
        .unwrap()
        .reversible_energy_j;
    close(change + u1 - u0, external_work, limits.energy_tolerance_j);
    close(
        c.report.contact_storage_change_j,
        u1 - u0,
        limits.energy_tolerance_j,
    );
    close_vec(linear, expected_linear, limits.linear_momentum_tolerance_ns);
    close_vec(
        angular,
        expected_angular,
        limits.angular_momentum_tolerance_nms,
    );
    for i in 0..3 {
        close(
            q1[i] - q0[i],
            c.report.relative_displacement_m[i],
            limits.displacement_tolerance_m,
        );
        close(
            c.report.motion.force_time_moment_n_s2[i],
            0.5 * dt * c.report.motion.force_impulse_n_s[i],
            1e-12 * MASS * scales(source).0,
        );
    }
}

#[test]
fn zero_shear_recovers_original_normal_discrete_motion_and_zero_load_api() {
    let (source, shear) = sources(1e9);
    let (length, speed, time) = scales(&source);
    let f = axial();
    let mut pair = f.prepared(&source, &shear);
    let initial = pair.initialize([0.0; 2], &CancelGate::new()).unwrap();
    let dt = 0.1 * time;
    let c = pair.advance(initial, dt, &CancelGate::new()).unwrap();
    assert_eq!(
        c,
        pair.advance_with_wrenches(initial, [Wrench::ZERO; 2], dt, &CancelGate::new())
            .unwrap()
    );
    let original = JkrNormalMotion::new(
        source.normal_source(),
        pair.mobility_per_kg()[0][0].recip(),
        JkrNormalMotionBudget {
            maximum_step_s: 2.0 * time,
            maximum_substep_s: dt,
            maximum_substeps: 16,
            maximum_root_iterations: 64,
            maximum_events: MAX_JKR_NORMAL_EVENTS,
            minimum_indentation_m: -20.0 * length,
            maximum_accumulated_travel_m: 20.0 * length,
            maximum_normal_speed_m_s: 100.0 * speed,
            maximum_force_n: 1e4 * FORCE,
            coordinate_absolute_tolerance_m: 1e-13 * length,
            coordinate_relative_tolerance: 2e-12,
            energy_absolute_tolerance_j: 1e-10 * FORCE * length,
            energy_relative_tolerance: 1e-10,
            impulse_absolute_tolerance_n_s: 1e-12 * MASS * speed,
            impulse_relative_tolerance: 1e-11,
        },
    )
    .unwrap();
    let pure = original
        .initialize(
            source
                .normal_source()
                .evaluate_attached(pair.approach_m())
                .unwrap()
                .state(),
            0.0,
        )
        .unwrap();
    let pure = original.advance(pure, 0.0, dt, &CancelGate::new()).unwrap();
    close(
        c.state_after.position_m()[0],
        pure.state().contact_state().indentation_m(),
        1e-9 * length,
    );
    close(
        c.state_after.velocity_m_s()[0],
        pure.state().velocity_m_s(),
        1e-9 * speed,
    );
    close(
        c.report.motion.contact.contact_radius_m,
        pure.contact.contact_radius_m,
        1e-10 * pure.contact.contact_radius_m,
    );
    assert_eq!(&c.state_after.position_m()[1..], &[0.0, 0.0]);
    assert_eq!(&c.report.motion.mean_force_n[1..], &[0.0, 0.0]);
    check_native(&f, &source, &c);
}

#[test]
fn elastic_shear_unloading_returns_native_energy_and_restores_real_contact_area() {
    let (source, shear) = sources(1e9);
    let (length, speed, time) = scales(&source);
    let mut f = axial();
    f.set_velocity(0, Vec3::new(-0.8 * speed, 0.0, 0.0));
    let mut pair = f.prepared(&source, &shear);
    let initial = pair
        .initialize([0.8 * length, 0.0], &CancelGate::new())
        .unwrap();
    let before = source
        .evaluate_displacement_controlled(pair.approach_m(), [0.8 * length, 0.0])
        .unwrap();
    let preload = [
        Wrench {
            force_world: Vec3::new(0.0, 0.0, before.normal_force_n),
            torque_body: Vec3::ZERO,
        },
        Wrench {
            force_world: Vec3::new(0.0, 0.0, -before.normal_force_n),
            torque_body: Vec3::ZERO,
        },
    ];
    let c = pair
        .advance_with_wrenches(initial, preload, 0.1 * time, &CancelGate::new())
        .unwrap();
    assert!(c.state_after.position_m()[1] < initial.position_m()[1]);
    assert!(c.report.motion.contact.contact_area_m2 > before.contact_area_m2);
    assert!(c.report.motion.mean_force_n[1] * c.report.relative_displacement_m[1] < 0.0);
    assert!(
        c.report.force_impulse_on_a_world_ns.x
            * (c.report.relative_velocity_m_s[0][1] + c.report.relative_velocity_m_s[1][1])
            > 0.0
    );
    assert!(c.report.kinetic_energy_change_j > 0.0);
    assert!(c.report.contact_storage_change_j < 0.0);
    check_native(&f, &source, &c);
}

#[test]
fn all_three_point_axes_keep_full_tensor_cross_recoil_and_world_covariance() {
    let (source, shear) = sources(1e9);
    let (length, speed, time) = scales(&source);
    // Opposite arms and proportionate full tensors give identical angular
    // recoil on the two bodies, so no unrepresented RELATIVE twist is driven.
    let kappa = 1.5;
    let arm = Vec3::new(0.4, 0.25, 1.0).scale(ARM);
    let mut f = fixture([arm, arm.scale(-kappa)]);
    let principal = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -1.0, 2.0), 0.6).unwrap();
    for (i, scale) in [(0, 1.0), (1, kappa)] {
        f.properties[i] = MassProperties::from_principal_frame(
            (2.0 + i as f64) * MASS,
            Vec3::ZERO,
            Vec3::new(0.8, 1.1, 1.4).scale(scale * INERTIA),
            principal,
        )
        .unwrap();
    }
    f.set_velocity(0, Vec3::new(0.8, -0.4, 0.2).scale(speed));
    f.set_velocity(1, Vec3::new(-0.1, 0.15, -0.05).scale(speed));
    let torque = Vec3::new(0.3, -0.4, 0.2).scale(FORCE * ARM);
    let wrenches = [
        Wrench {
            force_world: Vec3::new(0.6, -0.3, 1.0).scale(FORCE),
            torque_body: torque,
        },
        Wrench {
            force_world: Vec3::new(-0.2, 0.1, -0.5).scale(FORCE),
            torque_body: torque.scale(kappa),
        },
    ];
    let mut pair = f.prepared(&source, &shear);
    let axes = [f.frame.normal(), f.frame.longitudinal(), f.frame.lateral()].map(vector);
    let mut expected = [[0.0; 3]; 3];
    for body in 0..2 {
        let moments = f.properties[body].principal_inertia_body();
        let inverse = |value: Vec3| {
            let value = principal.rotate_world_to_body(value);
            principal.rotate_body_to_world(Vec3::new(
                value.x / moments.x,
                value.y / moments.y,
                value.z / moments.z,
            ))
        };
        for i in 0..3 {
            for j in 0..3 {
                expected[i][j] += axes[i].dot(axes[j]) / f.properties[body].mass()
                    + f.arms[body]
                        .cross(axes[i])
                        .dot(inverse(f.arms[body].cross(axes[j])));
            }
        }
    }
    for i in 0..3 {
        for j in 0..3 {
            close(pair.mobility_per_kg()[i][j], expected[i][j], 1e-9);
        }
    }
    assert!(expected[0][1].abs() > 0.05 * expected[0][0]);
    assert!(expected[0][2].abs() > 0.05 * expected[0][0]);
    let initial = pair
        .initialize([0.6 * length, -0.3 * length], &CancelGate::new())
        .unwrap();
    let dt = 0.1 * time;
    let c = pair
        .advance_with_wrenches(initial, wrenches, dt, &CancelGate::new())
        .unwrap();
    let force = c.report.motion.mean_force_n;
    let shear_recoil = -dt * (expected[0][1] * force[1] + expected[0][2] * force[2]);
    assert!(shear_recoil.abs() > 0.01 * speed);
    close(
        c.report.relative_velocity_m_s[1][0]
            - c.report.relative_velocity_m_s[0][0]
            - dt * c.report.driving_acceleration_m_s2[0]
            + dt * expected[0][0] * force[0],
        shear_recoil,
        1e-10 * speed,
    );
    let pure = source
        .normal_source()
        .evaluate_attached(c.state_after.position_m()[0])
        .unwrap();
    assert!(c.report.motion.contact.contact_area_m2 < pure.contact_area_m2);
    check_native(&f, &source, &c);
    let mut swapped = f.contact(&source, &shear);
    swapped.binding.body_a = 1;
    swapped.binding.body_b = 0;
    // The physical result is independent of container order, but a retained
    // native token cannot silently acquire a different ordered body identity.
    assert!(
        prepare_rigid_jkr_mixed_mode(
            &[f.bodies[1], f.bodies[0]],
            &[f.properties[1], f.properties[0]],
            swapped,
            initial,
            wrenches,
            dt,
            controls(&source),
            &CancelGate::new()
        )
        .is_err()
    );
    let mut swapped_pair = PreparedRigidJkrMixedModePair::new(
        &[f.bodies[1], f.bodies[0]],
        &[f.properties[1], f.properties[0]],
        swapped,
        controls(&source),
        &CancelGate::new(),
    )
    .unwrap();
    let swapped_initial = swapped_pair
        .initialize([0.6 * length, -0.3 * length], &CancelGate::new())
        .unwrap();
    let reordered = swapped_pair
        .advance_with_wrenches(swapped_initial, wrenches, dt, &CancelGate::new())
        .unwrap();
    assert_eq!(reordered.bodies_after, c.bodies_after);
    assert_eq!(
        reordered.state_after.motion_state(),
        c.state_after.motion_state()
    );
    assert_eq!(reordered.report, c.report);

    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -1.0), 0.7).unwrap();
    let mut rotated = fixture(f.arms);
    rotated.properties = f.properties;
    rotated.frame = TangentFrame::new(
        {
            let n = rotation.rotate_body_to_world(axes[0]);
            [n.x, n.y, n.z]
        },
        {
            let t = rotation.rotate_body_to_world(axes[1]);
            [t.x, t.y, t.z]
        },
    )
    .unwrap();
    for i in 0..2 {
        rotated.bodies[i] = RigidBodyState::new(
            Pose::new(
                rotation.rotate_body_to_world(f.arms[i]).scale(-1.0),
                rotation,
            )
            .unwrap(),
            rotation.rotate_body_to_world(f.bodies[i].linear_momentum_world()),
            f.bodies[i].angular_momentum_body(),
        )
        .unwrap();
    }
    let rotated_wrenches = wrenches.map(|w| Wrench {
        force_world: rotation.rotate_body_to_world(w.force_world),
        torque_body: w.torque_body,
    });
    let mut rotated_pair = rotated.prepared(&source, &shear);
    assert_eq!(rotated_pair.approach_m(), pair.approach_m());
    let rotated_initial = rotated_pair
        .initialize([0.6 * length, -0.3 * length], &CancelGate::new())
        .unwrap();
    let rotated_c = rotated_pair
        .advance_with_wrenches(rotated_initial, rotated_wrenches, dt, &CancelGate::new())
        .unwrap();
    for i in 0..2 {
        close_vec(
            rotated_c.bodies_after[i].linear_momentum_world(),
            rotation.rotate_body_to_world(c.bodies_after[i].linear_momentum_world()),
            1e-10 * MASS * speed,
        );
        close_vec(
            rotated_c.bodies_after[i].angular_momentum_body(),
            c.bodies_after[i].angular_momentum_body(),
            1e-10 * MASS * speed * ARM,
        );
        close_vec(
            rotated_c.report.body_displacement_world_m[i],
            rotation.rotate_body_to_world(c.report.body_displacement_world_m[i]),
            1e-10 * length,
        );
    }
    for i in 0..3 {
        close(
            rotated_c.state_after.position_m()[i],
            c.state_after.position_m()[i],
            1e-10 * length,
        );
    }
    check_native(&rotated, &source, &rotated_c);
}

#[test]
fn held_mixed_preload_preserves_shear_history_and_absolute_native_chronology() {
    let (source, shear) = sources(1e9);
    let (length, speed, time) = scales(&source);
    let f = axial();
    let mut pair = f.prepared(&source, &shear);
    let initial = pair
        .initialize([0.5 * length, -0.3 * length], &CancelGate::new())
        .unwrap();
    let response = source
        .evaluate_displacement_controlled(
            initial.position_m()[0],
            [initial.position_m()[1], initial.position_m()[2]],
        )
        .unwrap();
    let force = Vec3::new(
        response.tangential_force_n[0],
        response.tangential_force_n[1],
        response.normal_force_n,
    );
    let wrenches = core::array::from_fn(|i| {
        let force = force.scale(if i == 0 { 1.0 } else { -1.0 });
        Wrench {
            force_world: force,
            torque_body: f.arms[i].cross(force),
        }
    });
    let first = pair
        .advance_with_wrenches(initial, wrenches, 0.1 * time, &CancelGate::new())
        .unwrap();
    check_native(&f, &source, &first);
    let next = Fixture {
        bodies: first.bodies_after,
        properties: f.properties,
        arms: f.arms,
        frame: f.frame,
    };
    let second = next
        .prepared(&source, &shear)
        .advance_with_wrenches(first.state_after, wrenches, 0.1 * time, &CancelGate::new())
        .unwrap();
    assert_eq!(
        second.report.motion.initial_state(),
        first.state_after.motion_state()
    );
    assert_eq!(second.state_after.steps(), 2);
    assert_eq!(
        second.state_after.time_s(),
        first.state_after.time_s() + 0.1 * time
    );
    for i in 0..3 {
        close(
            second.state_after.position_m()[i],
            initial.position_m()[i],
            1e-10 * length,
        );
        close(second.state_after.velocity_m_s()[i], 0.0, 1e-10 * speed);
    }
    check_native(&next, &source, &second);
}

#[test]
fn frame_source_mobility_geometry_twist_and_load_refusals_are_atomic_and_retryable() {
    let (source, shear) = sources(1e9);
    let (length, _, time) = scales(&source);
    let f = axial();
    let initial = f
        .prepared(&source, &shear)
        .initialize([0.2 * length, 0.0], &CancelGate::new())
        .unwrap();
    let mut native = f.bodies;
    let mut history = initial;
    let mut contact = f.contact(&source, &shear);
    contact.frame = TangentFrame::new([0.0, 0.0, 1.0], [-1.0, 0.0, 0.0]).unwrap();
    assert!(matches!(
        apply_rigid_jkr_mixed_mode(
            &mut native,
            &f.properties,
            contact,
            &mut history,
            [Wrench::ZERO; 2],
            0.1 * time,
            controls(&source),
            &CancelGate::new()
        ),
        Err(RigidJkrMixedModeError::StateMismatch)
    ));
    assert_eq!(native, f.bodies);
    assert_eq!(history, initial);
    let mut changed = f.properties;
    changed[0] =
        MassProperties::new(2.1 * MASS, Vec3::ZERO, changed[0].principal_inertia_body()).unwrap();
    assert!(
        apply_rigid_jkr_mixed_mode(
            &mut native,
            &changed,
            f.contact(&source, &shear),
            &mut history,
            [Wrench::ZERO; 2],
            0.1 * time,
            controls(&source),
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(native, f.bodies);
    assert_eq!(history, initial);
    let (other_source, other_shear) = sources(1.1e9);
    assert!(
        apply_rigid_jkr_mixed_mode(
            &mut native,
            &f.properties,
            f.contact(&other_source, &other_shear),
            &mut history,
            [Wrench::ZERO; 2],
            0.1 * time,
            controls(&source),
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(native, f.bodies);
    assert_eq!(history, initial);
    let mut bad = [Wrench::ZERO; 2];
    bad[0].force_world.x = f64::NAN;
    assert!(
        apply_rigid_jkr_mixed_mode(
            &mut native,
            &f.properties,
            f.contact(&source, &shear),
            &mut history,
            bad,
            0.1 * time,
            controls(&source),
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(native, f.bodies);
    assert_eq!(history, initial);
    let mut twisting = [Wrench::ZERO; 2];
    twisting[0].torque_body.z = FORCE * ARM;
    assert!(matches!(
        apply_rigid_jkr_mixed_mode(
            &mut native,
            &f.properties,
            f.contact(&source, &shear),
            &mut history,
            twisting,
            0.1 * time,
            controls(&source),
            &CancelGate::new()
        ),
        Err(RigidJkrMixedModeError::UnsupportedTorsion { .. })
    ));
    assert_eq!(native, f.bodies);
    assert_eq!(history, initial);
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        apply_rigid_jkr_mixed_mode(
            &mut native,
            &f.properties,
            f.contact(&source, &shear),
            &mut history,
            [Wrench::ZERO; 2],
            0.1 * time,
            controls(&source),
            &cancelled
        )
        .is_err()
    );
    assert_eq!(native, f.bodies);
    assert_eq!(history, initial);
    let mut pair = f.prepared(&source, &shear);
    let expected = pair
        .advance(initial, 0.1 * time, &CancelGate::new())
        .unwrap();
    let report = apply_rigid_jkr_mixed_mode(
        &mut native,
        &f.properties,
        f.contact(&source, &shear),
        &mut history,
        [Wrench::ZERO; 2],
        0.1 * time,
        controls(&source),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(native, expected.bodies_after);
    assert_eq!(history, expected.state_after);
    assert_eq!(report, expected.report);
    assert!(matches!(
        prepare_rigid_jkr_mixed_mode(
            &native,
            &f.properties,
            f.contact(&source, &shear),
            history,
            [Wrench::ZERO; 2],
            0.1 * time,
            controls(&source),
            &CancelGate::new()
        ),
        Err(RigidJkrMixedModeError::GeometryMismatch { .. })
    ));
    assert!(
        pair.initialize(
            [
                2.0 * source
                    .critical_tangential_displacement_m(pair.approach_m())
                    .unwrap(),
                0.0
            ],
            &CancelGate::new()
        )
        .is_err()
    );
    let (weak_source, weak_shear) = sources(1e4);
    let mut weak_pair = f.prepared(&weak_source, &weak_shear);
    let weak_initial = weak_pair.initialize([0.0; 2], &CancelGate::new()).unwrap();
    let mut forced = [Wrench::ZERO; 2];
    forced[0].force_world.x = 10.0 * FORCE;
    let mut weak_native = f.bodies;
    let mut weak_history = weak_initial;
    assert!(
        apply_rigid_jkr_mixed_mode(
            &mut weak_native,
            &f.properties,
            f.contact(&weak_source, &weak_shear),
            &mut weak_history,
            forced,
            0.1 * time,
            controls(&weak_source),
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(weak_native, f.bodies);
    assert_eq!(weak_history, weak_initial);
    assert_eq!(
        weak_pair
            .advance(weak_initial, 0.1 * time, &CancelGate::new())
            .unwrap()
            .state_after
            .position_m()[1],
        0.0
    );
}
