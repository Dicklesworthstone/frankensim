//! Original changing-area I/III storage and actual paired-body mechanics.
#[path = "../../fs-contact/tests/adhesive_dmt_coupled/common.rs"]
mod common;
#[path = "adhesive_jkr_torsional/geometry.rs"]
mod geometry;
#[path = "adhesive_jkr_torsional/temperature.rs"]
mod temperature;

use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::{
    JkrEnvelope, SourceJkrSphereContact,
    mixed_mode::{
        JkrMixedModeEnvelope, SourceJkrMixedMode, torsion::motion::JkrTorsionalMotionBudget,
    },
};
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::{RigidContactBinding, adhesive::jkr_torsional::*};
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

fn scales(source: &SourceJkrMixedMode) -> (f64, f64, f64, f64) {
    let length = source.normal_source().indentation_scale_m();
    let speed = (FORCE * length / (1.2 * MASS)).sqrt();
    let zero = source
        .evaluate_torsion_displacement_controlled(0.0, 0.0)
        .unwrap();
    let stiffness = zero.coupled_stiffness.unwrap()[1][1];
    (
        length,
        speed,
        length / speed,
        (FORCE * length / stiffness).sqrt(),
    )
}

fn controls(source: &SourceJkrMixedMode) -> RigidJkrTorsionalControls {
    let (length, speed, time, twist) = scales(source);
    let torque = FORCE * length / twist;
    RigidJkrTorsionalControls {
        motion: JkrTorsionalMotionBudget {
            maximum_step_s: 2.0 * time,
            maximum_coordinate_increment: [20.0 * length, 20.0 * twist],
            maximum_rate: [100.0 * speed, 100.0 * twist / time],
            maximum_effort: [1e4 * FORCE, 1e4 * torque],
            maximum_iterations: 64,
            maximum_evaluations: 2048,
            maximum_backtracks: 24,
            coordinate_absolute_tolerance: [1e-13 * length, 1e-13 * twist],
            coordinate_relative_tolerance: 2e-12,
            rate_absolute_tolerance: [1e-12 * speed, 1e-12 * twist / time],
            rate_relative_tolerance: 1e-11,
            energy_absolute_tolerance_j: 1e-10 * FORCE * length,
            energy_relative_tolerance: 1e-10,
        },
        maximum_bodies: 4,
        geometry_tolerance_m: 4e-18,
        displacement_tolerance_m: 1e-9 * length,
        twist_tolerance_rad: 1e-9 * twist,
        velocity_tolerance_m_s: 1e-9 * speed,
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
    fn contact<'a>(&self, source: &'a SourceJkrMixedMode) -> RigidJkrTorsionalContact<'a> {
        RigidJkrTorsionalContact {
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
    fn prepared<'a>(&self, source: &'a SourceJkrMixedMode) -> PreparedRigidJkrTorsionalPair<'a> {
        PreparedRigidJkrTorsionalPair::new(
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

fn check_native(f: &Fixture, source: &SourceJkrMixedMode, c: &RigidJkrTorsionalCandidate) {
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
}

#[test]
fn native_normal_and_twist_share_actual_area_and_external_work() {
    let source = source();
    let (length, speed, time, twist) = scales(&source);
    let torque = FORCE * length / twist;
    let mut f = axial();
    f.set_rates(
        0,
        Vec3::new(0.0, 0.0, 0.03 * speed),
        Vec3::new(0.0, 0.0, 0.02 * twist / time),
    );
    let pair = f.prepared(&source);
    let state = pair.initialize(0.2 * twist, &CancelGate::new()).unwrap();
    let wrenches = [
        Wrench {
            force_world: Vec3::new(0.0, 0.0, 0.25 * FORCE),
            torque_body: Vec3::new(0.0, 0.0, 0.3 * torque),
        },
        Wrench {
            force_world: Vec3::new(0.0, 0.0, -0.1 * FORCE),
            torque_body: Vec3::new(0.0, 0.0, -0.1 * torque),
        },
    ];
    let c = pair
        .advance_with_wrenches(state, wrenches, 0.04 * time, &CancelGate::new())
        .unwrap();
    check_native(&f, &source, &c);
    let end = c.state_after.coordinates();
    let actual = source
        .evaluate_torsion_displacement_controlled(end[0], end[1])
        .unwrap();
    let normal = source.normal_source().evaluate_attached(end[0]).unwrap();
    assert!(actual.contact_area_m2 < normal.contact_area_m2);
    assert!(actual.normal_force_n > normal.normal_force_n);
    assert!(actual.torsional_moment_n_m > 0.0);
    assert_ne!(
        c.report.motion.contact.contact_area_m2,
        c.report.motion.initial_contact.contact_area_m2
    );
    assert_eq!(actual.model_identity(), source.identity());
    close(
        actual.contact_area_m2,
        c.report.motion.contact.contact_area_m2,
        1e-12 * actual.contact_area_m2,
    );
    close(
        actual.reversible_energy_j - c.report.motion.initial_contact.reversible_energy_j,
        c.report.contact_storage_change_j,
        controls(&source).energy_tolerance_j,
    );
    close(c.report.mobility[0][0], 1.0 / (1.2 * MASS), 1e-10 / MASS);
    close(
        c.report.mobility[1][1],
        1.0 / (1.4 * INERTIA) + 1.0 / (2.0 * INERTIA),
        1e-10 / INERTIA,
    );
    assert_eq!(c.report.mobility[0][1], 0.0);
}

#[test]
fn zero_twist_reduces_to_original_normal_storage_without_spurious_torque() {
    let source = source();
    let f = axial();
    let pair = f.prepared(&source);
    let old = pair.initialize(0.0, &CancelGate::new()).unwrap();
    let c = pair
        .advance(old, 0.02 * scales(&source).2, &CancelGate::new())
        .unwrap();
    check_native(&f, &source, &c);
    let normal = source
        .normal_source()
        .evaluate_attached(c.state_after.coordinates()[0])
        .unwrap();
    close(
        c.report.motion.contact.normal_force_n,
        normal.normal_force_n,
        1e-11 * FORCE,
    );
    close(
        c.report.motion.contact.contact_radius_m,
        normal.contact_radius_m,
        1e-11 * normal.contact_radius_m,
    );
    assert_eq!(c.state_after.coordinates()[1], 0.0);
    assert_eq!(c.state_after.rates()[1], 0.0);
    assert_eq!(c.report.angular_impulse_on_a_world_nms, Vec3::ZERO);
}

#[test]
fn full_inertia_cross_mobility_matches_independent_impulses_and_held_preload() {
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
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, 0.0, 0.0), 0.45 + 0.1 * i as f64)
                .unwrap(),
        )
        .unwrap();
    }
    let pair = f.prepared(&source);
    let w = pair.mobility();
    assert_eq!(w[0][1], w[1][0]);
    assert!(w[0][1].abs() > 1e-3 * (w[0][0] * w[1][1]).sqrt());
    let normal = vector(f.frame.normal());
    let reaction = pair.reaction_arms_body();
    for column in 0..2 {
        let mut endpoint = f.bodies;
        for (i, sign) in [(0, 1.0), (1, -1.0)] {
            let force = if column == 0 {
                normal.scale(sign)
            } else {
                Vec3::ZERO
            };
            let couple = if column == 1 {
                normal.scale(sign)
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
        let projected = [
            a.relative_velocity_world(b).unwrap().dot(normal),
            a.angular_velocity_world
                .sub(b.angular_velocity_world)
                .dot(normal),
        ];
        for row in 0..2 {
            close(projected[row], w[row][column], 3e-13 * w[row][column].abs());
        }
    }
    // Exact held elastic preload on an anisotropic off-axis body. The free
    // acceleration contains both cross terms; diagonal substitution moves it.
    let twist = 0.15 * scales(&source).3;
    let response = source
        .evaluate_torsion_displacement_controlled(pair.approach_m(), twist)
        .unwrap();
    let state = pair.initialize(twist, &CancelGate::new()).unwrap();
    let wrenches = core::array::from_fn(|i| {
        let sign = if i == 0 { 1.0 } else { -1.0 };
        let force = normal.scale(sign * response.normal_force_n);
        Wrench {
            force_world: force,
            torque_body: reaction[i]
                .cross(force)
                .add(normal.scale(sign * response.torsional_moment_n_m)),
        }
    });
    let c = pair
        .advance_with_wrenches(
            state,
            wrenches,
            0.02 * scales(&source).2,
            &CancelGate::new(),
        )
        .unwrap();
    check_native(&f, &source, &c);
    for (i, tolerance) in [
        controls(&source).displacement_tolerance_m,
        controls(&source).twist_tolerance_rad,
    ]
    .into_iter()
    .enumerate()
    {
        close(
            c.state_after.coordinates()[i],
            state.coordinates()[i],
            tolerance,
        );
    }
    close(
        c.report.kinetic_energy_change_j,
        0.0,
        controls(&source).energy_tolerance_j,
    );
}

#[test]
fn elastic_release_is_assisting_work_and_rotates_with_the_world_frame() {
    let source = source();
    let (length, _, time, twist) = scales(&source);
    let f = axial();
    let pair = f.prepared(&source);
    let state = pair.initialize(0.35 * twist, &CancelGate::new()).unwrap();
    let c = pair
        .advance(state, 0.03 * time, &CancelGate::new())
        .unwrap();
    check_native(&f, &source, &c);
    assert!(c.report.kinetic_energy_change_j > 0.0);
    assert!(c.report.contact_storage_change_j < 0.0);
    assert!(c.report.contact_impulse_work_j > 0.0);
    let rotation = UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, -1.0), 0.73).unwrap();
    let mut rotated = axial();
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
    let rotated_pair = rotated.prepared(&source);
    let rs = rotated_pair
        .initialize(0.35 * twist, &CancelGate::new())
        .unwrap();
    let rc = rotated_pair
        .advance(rs, 0.03 * time, &CancelGate::new())
        .unwrap();
    check_native(&rotated, &source, &rc);
    close(
        rc.state_after.coordinates()[0],
        c.state_after.coordinates()[0],
        1e-9 * length,
    );
    close(
        rc.state_after.coordinates()[1],
        c.state_after.coordinates()[1],
        1e-9 * twist,
    );
    close_vec(
        rc.report.force_impulse_on_a_world_ns,
        rotation.rotate_body_to_world(c.report.force_impulse_on_a_world_ns),
        controls(&source).linear_momentum_tolerance_ns,
    );
    close_vec(
        rc.report.angular_impulse_on_a_world_nms,
        rotation.rotate_body_to_world(c.report.angular_impulse_on_a_world_nms),
        controls(&source).angular_momentum_tolerance_nms,
    );
}

#[test]
fn stale_geometry_unrepresented_modes_and_late_refusals_publish_nothing() {
    let source = source();
    let (_, speed, time, twist) = scales(&source);
    let f = axial();
    let pair = f.prepared(&source);
    let original = pair.initialize(0.1 * twist, &CancelGate::new()).unwrap();
    let dt = 0.02 * time;
    let expected = pair.advance(original, dt, &CancelGate::new()).unwrap();
    let mut bodies = f.bodies;
    let mut state = original;
    let bad = [
        Wrench {
            force_world: Vec3::new(FORCE, 0.0, 0.0),
            torque_body: Vec3::ZERO,
        },
        Wrench::ZERO,
    ];
    assert!(matches!(
        apply_rigid_jkr_torsional(
            &mut bodies,
            &f.properties,
            f.contact(&source),
            &mut state,
            bad,
            dt,
            controls(&source),
            &CancelGate::new()
        ),
        Err(RigidJkrTorsionalError::UnsupportedShear { .. })
    ));
    assert_eq!(bodies, f.bodies);
    assert_eq!(state, original);
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        apply_rigid_jkr_torsional(
            &mut bodies,
            &f.properties,
            f.contact(&source),
            &mut state,
            [Wrench::ZERO; 2],
            dt,
            controls(&source),
            &cancelled
        )
        .is_err()
    );
    assert_eq!(bodies, f.bodies);
    assert_eq!(state, original);
    apply_rigid_jkr_torsional(
        &mut bodies,
        &f.properties,
        f.contact(&source),
        &mut state,
        [Wrench::ZERO; 2],
        dt,
        controls(&source),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(bodies, expected.bodies_after);
    assert_eq!(state, expected.state_after);
    let next = PreparedRigidJkrTorsionalPair::new(
        &bodies,
        &f.properties,
        f.contact(&source),
        controls(&source),
        &CancelGate::new(),
    )
    .unwrap();
    assert!(matches!(
        next.advance(state, dt, &CancelGate::new()),
        Err(RigidJkrTorsionalError::GeometryMismatch { .. })
    ));
    let mut shear = axial();
    shear.set_rates(0, Vec3::new(speed, 0.0, 0.0), Vec3::ZERO);
    assert!(matches!(
        shear.prepared(&source).initialize(0.0, &CancelGate::new()),
        Err(RigidJkrTorsionalError::UnsupportedShear { .. })
    ));
    let mut rolling = fixture([Vec3::ZERO; 2]);
    rolling.set_rates(0, Vec3::ZERO, Vec3::new(twist / time, 0.0, 0.0));
    assert!(matches!(
        rolling
            .prepared(&source)
            .initialize(0.0, &CancelGate::new()),
        Err(RigidJkrTorsionalError::UnsupportedRolling { .. })
    ));
}
