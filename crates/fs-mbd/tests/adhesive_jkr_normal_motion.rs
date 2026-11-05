//! G1/G3/G4: original JKR events drive actual finite-body momentum and travel.
#[path = "../../fs-contact/tests/source_point/common.rs"]
mod common;

use core::f64::consts::PI;
use fs_contact::normal_patch::adhesive::{
    JkrBranch, JkrEnvelope, SourceJkrSphereContact,
    motion::{JkrNormalMotionBudget, MAX_JKR_NORMAL_EVENTS},
};
use fs_exec::CancelGate;
use fs_mbd::tangential_friction::{RigidContactBinding, adhesive::jkr_normal_motion::*};
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};
use fs_tribo::partial_slip::TangentFrame;

const MASS: f64 = 1e-5;
const REDUCED_MASS: f64 = 1.2 * MASS;
const ARM: f64 = 0.01;
const INERTIA: f64 = MASS * ARM * ARM;
const FORCE: f64 = PI * 0.05 * 0.01 / 6.0;

fn source(young: f64) -> SourceJkrSphereContact {
    let a = common::solid("A", 1000.0, young);
    let b = common::solid("B", 2000.0, 2e7);
    SourceJkrSphereContact::new(
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
    .unwrap()
}
fn scales(normal: &SourceJkrSphereContact) -> (f64, f64, f64) {
    let length = normal.indentation_scale_m();
    let speed = (FORCE * length / REDUCED_MASS).sqrt();
    (length, speed, length / speed)
}
fn controls(normal: &SourceJkrSphereContact) -> RigidJkrControls {
    let (length, speed, time) = scales(normal);
    RigidJkrControls {
        motion: JkrNormalMotionBudget {
            maximum_step_s: 10.0 * time,
            maximum_substep_s: 0.01 * time,
            maximum_substeps: 4096,
            maximum_root_iterations: 48,
            maximum_events: MAX_JKR_NORMAL_EVENTS,
            minimum_indentation_m: -100.0 * length,
            maximum_accumulated_travel_m: 200.0 * length,
            maximum_normal_speed_m_s: 100.0 * speed,
            maximum_force_n: 1e4 * FORCE,
            coordinate_absolute_tolerance_m: 1e-14 * length,
            coordinate_relative_tolerance: 1e-13,
            energy_absolute_tolerance_j: 1e-11 * FORCE * length,
            energy_relative_tolerance: 1e-10,
            impulse_absolute_tolerance_n_s: 1e-12 * REDUCED_MASS * speed,
            impulse_relative_tolerance: 1e-11,
        },
        maximum_bodies: 4,
        geometry_tolerance_m: 2e-16,
        displacement_tolerance_m: 1e-9 * length,
        velocity_tolerance_m_s: 1e-9 * speed,
        linear_momentum_tolerance_ns: 1e-12 * MASS * speed,
        angular_momentum_tolerance_nms: 1e-12 * MASS * speed * ARM,
        energy_tolerance_j: 1e-8 * FORCE * length,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}
fn close_vec(actual: Vec3, expected: Vec3, tolerance: f64) {
    close(actual.x, expected.x, tolerance);
    close(actual.y, expected.y, tolerance);
    close(actual.z, expected.z, tolerance);
}

struct Fixture {
    bodies: [RigidBodyState; 2],
    properties: [MassProperties; 2],
    binding: RigidContactBinding,
    frame: TangentFrame,
}
impl Fixture {
    fn contact<'a>(&self, normal: &'a SourceJkrSphereContact) -> RigidJkrContact<'a> {
        RigidJkrContact {
            binding: self.binding,
            frame: self.frame,
            reference_approach_m: 0.0,
            normal,
        }
    }
    fn prepared<'a>(&self, normal: &'a SourceJkrSphereContact) -> PreparedRigidJkrPair<'a> {
        PreparedRigidJkrPair::new(
            &self.bodies,
            &self.properties,
            self.contact(normal),
            controls(normal),
            &CancelGate::new(),
        )
        .unwrap()
    }
    fn set_motion(&mut self, body: usize, velocity: Vec3, omega_body: Vec3) {
        self.bodies[body] = RigidBodyState::new(
            self.bodies[body].pose(),
            velocity.scale(self.properties[body].mass()),
            self.properties[body]
                .angular_momentum_body_checked(omega_body)
                .unwrap(),
        )
        .unwrap();
    }
    fn advance_poses(&mut self, candidate: &RigidJkrCandidate) {
        // This fixture's axial translation has no rotation. The production
        // phase deliberately leaves the real pose integration to its owner.
        for i in 0..2 {
            assert_eq!(
                candidate.report.angular_displacement_body_rad[i],
                Vec3::ZERO
            );
            let body = candidate.bodies_after[i];
            let pose = Pose::new(
                body.pose()
                    .position_world()
                    .add(candidate.report.body_displacement_world_m[i]),
                body.pose().orientation(),
            )
            .unwrap();
            self.bodies[i] = RigidBodyState::new(
                pose,
                body.linear_momentum_world(),
                body.angular_momentum_body(),
            )
            .unwrap();
        }
    }
}
fn fixture(approach: f64, arms: [Vec3; 2], orientation: [UnitQuaternion; 2]) -> Fixture {
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
    let bodies = core::array::from_fn(|i| {
        RigidBodyState::new(
            Pose::new(
                orientation[i].rotate_body_to_world(arms[i]).scale(-1.0),
                orientation[i],
            )
            .unwrap(),
            Vec3::ZERO,
            Vec3::ZERO,
        )
        .unwrap()
    });
    Fixture {
        bodies,
        properties,
        binding: RigidContactBinding {
            body_a: 0,
            body_b: 1,
            arm_a_body: arms[0].add(orientation[0].rotate_world_to_body(Vec3::new(
                0.0,
                0.0,
                0.5 * approach,
            ))),
            arm_b_body: arms[1].sub(orientation[1].rotate_world_to_body(Vec3::new(
                0.0,
                0.0,
                0.5 * approach,
            ))),
        },
        frame: TangentFrame::new([0.0, 0.0, 1.0], [1.0, 0.0, 0.0]).unwrap(),
    }
}
fn axial(approach: f64) -> Fixture {
    fixture(
        approach,
        [Vec3::new(0.0, 0.0, ARM), Vec3::new(0.0, 0.0, -ARM)],
        [UnitQuaternion::IDENTITY; 2],
    )
}
fn check_native(f: &Fixture, candidate: &RigidJkrCandidate, normal: &SourceJkrSphereContact) {
    let (length, speed, _) = scales(normal);
    let impulse = candidate.report.force_impulse_on_a_world_ns;
    let mut change = 0.0;
    for (i, sign) in [(0, 1.0), (1, -1.0)] {
        let before = f.bodies[i];
        let after = candidate.bodies_after[i];
        let expected = impulse.scale(sign);
        close_vec(
            after
                .linear_momentum_world()
                .sub(before.linear_momentum_world()),
            expected,
            1e-12 * MASS * speed,
        );
        let arm = if i == 0 {
            f.binding.arm_a_body
        } else {
            f.binding.arm_b_body
        };
        let torque = arm.cross(before.pose().orientation().rotate_world_to_body(expected));
        close_vec(
            after
                .angular_momentum_body()
                .sub(before.angular_momentum_body()),
            torque,
            1e-12 * MASS * speed * ARM,
        );
        assert_eq!(before.pose(), after.pose());
        change += after.kinetic_energy(f.properties[i]).unwrap()
            - before.kinetic_energy(f.properties[i]).unwrap();
    }
    close(
        change,
        candidate.report.kinetic_energy_change_j,
        1e-11 * FORCE * length,
    );
    close(
        change
            + candidate.report.normal_storage_change_j
            + candidate.report.motion.released_energy_j,
        0.0,
        1e-8 * FORCE * length,
    );
    let dt = candidate.report.motion.state().time_s();
    let center_travel = candidate.report.body_displacement_world_m[0]
        .scale(f.properties[0].mass())
        .add(candidate.report.body_displacement_world_m[1].scale(f.properties[1].mass()));
    close_vec(
        center_travel,
        f.bodies[0]
            .linear_momentum_world()
            .add(f.bodies[1].linear_momentum_world())
            .scale(dt),
        1e-11 * MASS * length,
    );
}

#[test]
fn detached_flight_preserves_free_tangent_and_spin_and_actual_travel() {
    let normal = source(1e7);
    let (length, speed, time) = scales(&normal);
    let mut f = axial(-8.0 * length);
    f.set_motion(
        0,
        Vec3::new(2.0 * speed, -speed, -speed),
        Vec3::new(0.0, 0.0, speed / ARM),
    );
    f.set_motion(1, Vec3::new(speed, 0.0, 0.0), Vec3::ZERO);
    let pair = f.prepared(&normal);
    let initial = normal.detached_state(pair.approach_m()).unwrap();
    let dt = 0.5 * time;
    let c = pair.advance(initial, dt, &CancelGate::new()).unwrap();
    assert_eq!(c.bodies_after, f.bodies);
    assert_eq!(c.state_after.branch(), JkrBranch::Detached);
    assert_eq!(c.report.motion.event_count(), 0);
    assert_eq!(c.report.motion.released_energy_j, 0.0);
    assert_eq!(c.report.force_impulse_on_a_world_ns, Vec3::ZERO);
    close(
        c.report.relative_normal_displacement_m,
        -speed * dt,
        1e-13 * length,
    );
    for i in 0..2 {
        close_vec(
            c.report.body_displacement_world_m[i],
            f.bodies[i]
                .center_of_mass_velocity_world(f.properties[i])
                .unwrap()
                .scale(dt),
            1e-14 * length,
        );
    }
    close(
        c.report.angular_displacement_body_rad[0].z,
        speed / ARM * dt,
        1e-15,
    );
    check_native(&f, &c, &normal);
}

#[test]
fn attached_attraction_uses_actual_reduced_mass_and_preserves_native_energy() {
    let normal = source(1e7);
    let (length, speed, time) = scales(&normal);
    let mut f = axial(0.0);
    f.set_motion(
        0,
        Vec3::new(speed, -speed, 0.0),
        Vec3::new(0.0, 0.0, speed / ARM),
    );
    let pair = f.prepared(&normal);
    close(pair.mobility_per_kg(), 1.0 / REDUCED_MASS, 1e-10);
    let initial = normal.evaluate_attached(pair.approach_m()).unwrap().state();
    let dt = 1e-4 * time;
    let c = pair.advance(initial, dt, &CancelGate::new()).unwrap();
    // Independent JKR identity at delta=0: F = -4*pi*W*R/3 = -8*scale.
    close(
        c.report.force_impulse_on_a_world_ns.z / dt,
        8.0 * FORCE,
        1e-7 * FORCE,
    );
    assert!(c.bodies_after[0].linear_momentum_world().z > 0.0);
    assert!(c.bodies_after[1].linear_momentum_world().z < 0.0);
    for i in 0..2 {
        assert_eq!(
            c.bodies_after[i].linear_momentum_world().x,
            f.bodies[i].linear_momentum_world().x
        );
        assert_eq!(
            c.bodies_after[i].linear_momentum_world().y,
            f.bodies[i].linear_momentum_world().y
        );
        assert_eq!(
            c.bodies_after[i].angular_momentum_body(),
            f.bodies[i].angular_momentum_body()
        );
    }
    assert!(c.report.kinetic_energy_change_j > 0.0);
    assert!(c.report.normal_storage_change_j < 0.0);
    assert_eq!(c.report.motion.event_count(), 0);
    close(
        c.report.relative_normal_velocity_m_s[1],
        c.report.force_impulse_on_a_world_ns.z / REDUCED_MASS,
        1e-12 * speed,
    );
    close(
        c.state_after.indentation_m(),
        c.report.relative_normal_displacement_m,
        1e-13 * length,
    );
    check_native(&f, &c, &normal);
}

#[test]
fn capture_and_rebound_transfer_original_snap_energy_without_artificial_impulses() {
    let normal = source(1e7);
    let (length, speed, time) = scales(&normal);
    for (speed_factor, events, branch) in [
        (5.0, 1, JkrBranch::Attached),
        (10.0, 2, JkrBranch::Detached),
    ] {
        let mut f = axial(-0.25 * length);
        f.set_motion(0, Vec3::new(0.0, 0.0, speed_factor * speed), Vec3::ZERO);
        let pair = f.prepared(&normal);
        let old = normal.detached_state(pair.approach_m()).unwrap();
        let c = pair.advance(old, 8.0 * time, &CancelGate::new()).unwrap();
        assert_eq!(c.state_after.branch(), branch);
        assert_eq!(c.report.motion.event_count(), events);
        let expected_release =
            FORCE * length * (3.6 * 4.0_f64.powf(4.0 / 3.0) + if events == 2 { 2.4 } else { 0.0 });
        close(
            c.report.motion.released_energy_j,
            expected_release,
            1e-11 * FORCE * length,
        );
        close(
            c.report.motion.events().next().unwrap().velocity_m_s,
            speed_factor * speed,
            1e-12 * speed,
        );
        if events == 2 {
            let expected_velocity =
                -((speed_factor * speed).powi(2) - 2.0 * expected_release / REDUCED_MASS).sqrt();
            close(
                c.report.relative_normal_velocity_m_s[1],
                expected_velocity,
                5e-9 * speed,
            );
        } else {
            assert!(c.report.relative_normal_velocity_m_s[1].abs() > 0.01 * speed);
        }
        assert!(
            c.report.motion.accumulated_travel_m > c.report.relative_normal_displacement_m.abs()
        );
        // h*I/2 would replace the whole curved/event trajectory by an affine
        // impulse phase and give a visibly different physical displacement.
        let midpoint_travel = 4.0
            * time
            * (c.report.relative_normal_velocity_m_s[0] + c.report.relative_normal_velocity_m_s[1]);
        assert!((midpoint_travel - c.report.relative_normal_displacement_m).abs() > 0.1 * length);
        check_native(&f, &c, &normal);
        assert_eq!(
            c,
            pair.advance(old, 8.0 * time, &CancelGate::new()).unwrap()
        );
    }
}

#[test]
fn off_axis_full_inertia_restores_all_native_travel_and_is_world_frame_covariant() {
    let normal = source(1e7);
    let (length, _, time) = scales(&normal);
    let orientations = [
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, 2.0, 3.0), 0.4).unwrap(),
        UnitQuaternion::from_axis_angle(Vec3::new(-2.0, 1.0, 1.0), -0.3).unwrap(),
    ];
    let arms = [
        Vec3::new(0.4, 0.2, 1.0).scale(ARM),
        Vec3::new(-0.3, 0.1, -1.0).scale(ARM),
    ];
    let mut f = fixture(0.0, arms, orientations);
    let principal = UnitQuaternion::from_axis_angle(Vec3::new(1.0, -1.0, 2.0), 0.6).unwrap();
    f.properties[0] = MassProperties::from_principal_frame(
        2.0 * MASS,
        Vec3::ZERO,
        Vec3::new(0.8, 1.1, 1.4).scale(INERTIA),
        principal,
    )
    .unwrap();
    let inverse = |i: usize, x: Vec3| {
        let p = f.properties[i];
        let rotation = if i == 0 {
            principal
        } else {
            UnitQuaternion::IDENTITY
        };
        let v = rotation.rotate_world_to_body(x);
        let moments = p.principal_inertia_body();
        rotation.rotate_body_to_world(Vec3::new(v.x / moments.x, v.y / moments.y, v.z / moments.z))
    };
    let n = Vec3::new(0.0, 0.0, 1.0);
    let independent_mobility: f64 = (0..2)
        .map(|i| {
            let torque = arms[i].cross(orientations[i].rotate_world_to_body(n));
            1.0 / f.properties[i].mass() + torque.dot(inverse(i, torque))
        })
        .sum();
    let pair = f.prepared(&normal);
    close(pair.mobility_per_kg(), independent_mobility, 2e-10);
    let initial = normal.evaluate_attached(pair.approach_m()).unwrap().state();
    let dt = 0.2 * time;
    let c = pair.advance(initial, dt, &CancelGate::new()).unwrap();
    let q = c
        .report
        .force_impulse_on_a_world_ns
        .scale(dt)
        .sub(c.report.force_time_moment_on_a_world_ns2);
    for (i, sign) in [(0, 1.0), (1, -1.0)] {
        let torque = arms[i].cross(orientations[i].rotate_world_to_body(q.scale(sign)));
        close_vec(
            c.report.angular_displacement_body_rad[i],
            inverse(i, torque),
            1e-12 * length / ARM,
        );
        assert!(c.report.angular_displacement_body_rad[i].norm_squared() > 0.0);
    }
    check_native(&f, &c, &normal);
    let world_rotation = Vec3::new(1.0, 1.0, -1.0).scale(0.7 / 3.0_f64.sqrt());
    let rotation = UnitQuaternion::IDENTITY.right_exp(world_rotation).unwrap();
    let mut rotated = fixture(0.0, arms, orientations);
    rotated.properties = f.properties;
    for i in 0..2 {
        let pose = f.bodies[i].pose();
        let orientation = pose
            .orientation()
            .right_exp(pose.orientation().rotate_world_to_body(world_rotation))
            .unwrap();
        rotated.bodies[i] = RigidBodyState::new(
            Pose::new(
                rotation.rotate_body_to_world(pose.position_world()),
                orientation,
            )
            .unwrap(),
            Vec3::ZERO,
            Vec3::ZERO,
        )
        .unwrap();
    }
    let normal_rotated = rotation.rotate_body_to_world(n);
    let tangent_rotated = rotation.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
    rotated.frame = TangentFrame::new(
        [normal_rotated.x, normal_rotated.y, normal_rotated.z],
        [tangent_rotated.x, tangent_rotated.y, tangent_rotated.z],
    )
    .unwrap();
    let rp = rotated.prepared(&normal);
    let rc = rp
        .advance(
            normal.evaluate_attached(rp.approach_m()).unwrap().state(),
            dt,
            &CancelGate::new(),
        )
        .unwrap();
    close(
        rc.report.relative_normal_displacement_m,
        c.report.relative_normal_displacement_m,
        1e-10 * length,
    );
    close_vec(
        rc.report.force_impulse_on_a_world_ns,
        rotation.rotate_body_to_world(c.report.force_impulse_on_a_world_ns),
        1e-12 * FORCE * time,
    );
    for i in 0..2 {
        close_vec(
            rc.report.body_displacement_world_m[i],
            rotation.rotate_body_to_world(c.report.body_displacement_world_m[i]),
            1e-10 * length,
        );
    }
}

#[test]
fn geometry_owned_continuation_and_body_indices_preserve_history_without_replaying_snap() {
    let normal = source(1e7);
    let (length, speed, time) = scales(&normal);
    let mut f = axial(-0.2 * length);
    f.set_motion(0, Vec3::new(0.0, 0.0, speed), Vec3::ZERO);
    let pair = f.prepared(&normal);
    let initial = normal.detached_state(pair.approach_m()).unwrap();
    let touch_time = -pair.approach_m() / speed;
    let first = pair
        .advance(initial, touch_time, &CancelGate::new())
        .unwrap();
    assert_eq!(first.report.motion.event_count(), 1);
    close_vec(
        first.report.force_impulse_on_a_world_ns,
        Vec3::ZERO,
        1e-12 * MASS * speed,
    );
    close(
        first.report.kinetic_energy_change_j,
        0.0,
        1e-12 * FORCE * length,
    );
    let mut stale = f.bodies;
    let mut history = initial;
    let report = apply_rigid_jkr_contact(
        &mut stale,
        &f.properties,
        f.contact(&normal),
        &mut history,
        touch_time,
        controls(&normal),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(report, first.report);
    assert_eq!(stale, first.bodies_after);
    assert_eq!(history, first.state_after);
    assert!(matches!(
        prepare_rigid_jkr_contact(
            &stale,
            &f.properties,
            f.contact(&normal),
            history,
            0.05 * time,
            controls(&normal),
            &CancelGate::new()
        ),
        Err(RigidJkrError::GeometryMismatch { .. })
    ));
    f.advance_poses(&first);
    let second = f
        .prepared(&normal)
        .advance(history, 0.05 * time, &CancelGate::new())
        .unwrap();
    assert_eq!(second.report.motion.event_count(), 0);
    assert_eq!(second.report.motion.released_energy_j, 0.0);
    let mut contact = f.contact(&normal);
    contact.binding.body_a = 1;
    contact.binding.body_b = 0;
    let swapped = prepare_rigid_jkr_contact(
        &[f.bodies[1], f.bodies[0]],
        &[f.properties[1], f.properties[0]],
        contact,
        history,
        0.05 * time,
        controls(&normal),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(second, swapped);
}

#[test]
fn source_geometry_budget_and_cancellation_refusals_are_atomic_and_retryable() {
    let normal = source(1e7);
    let other = source(1.1e7);
    let (_, _, time) = scales(&normal);
    let f = axial(0.0);
    let old = normal.evaluate_attached(0.0).unwrap().state();
    let pair = f.prepared(&normal);
    assert!(matches!(
        pair.advance(
            other.evaluate_attached(0.0).unwrap().state(),
            0.1 * time,
            &CancelGate::new()
        ),
        Err(RigidJkrError::Motion(_))
    ));
    let mut malformed = f.contact(&normal);
    malformed.binding.arm_a_body.x = 1e-3;
    assert!(
        PreparedRigidJkrPair::new(
            &f.bodies,
            &f.properties,
            malformed,
            controls(&normal),
            &CancelGate::new()
        )
        .is_err()
    );
    malformed = f.contact(&normal);
    malformed.binding.body_b = 0;
    assert!(
        PreparedRigidJkrPair::new(
            &f.bodies,
            &f.properties,
            malformed,
            controls(&normal),
            &CancelGate::new()
        )
        .is_err()
    );
    let mut native = f.bodies;
    let mut history = old;
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        apply_rigid_jkr_contact(
            &mut native,
            &f.properties,
            f.contact(&normal),
            &mut history,
            0.1 * time,
            controls(&normal),
            &cancelled
        )
        .is_err()
    );
    assert_eq!(native, f.bodies);
    assert_eq!(history, old);
    let mut limited = controls(&normal);
    limited.motion.maximum_substeps = 1;
    assert!(
        apply_rigid_jkr_contact(
            &mut native,
            &f.properties,
            f.contact(&normal),
            &mut history,
            time,
            limited,
            &CancelGate::new()
        )
        .is_err()
    );
    assert_eq!(native, f.bodies);
    assert_eq!(history, old);
    for invalid_dt in [0.0, -time, f64::NAN, 11.0 * time] {
        assert!(pair.advance(old, invalid_dt, &CancelGate::new()).is_err());
    }
    let first = pair.advance(old, 0.1 * time, &CancelGate::new()).unwrap();
    assert_eq!(
        first,
        pair.advance(old, 0.1 * time, &CancelGate::new()).unwrap()
    );
}
