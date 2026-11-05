//! Automatic sphere/plane discovery, coupled support and complete-path admission.
//! The oracles use the original rigid-body mass/inertia and the declared
//! first-order kick/contact/drift scheme, including its normal splitting loss.

#[path = "sphere_scene/connected.rs"]
mod connected;

use fs_contact::broad_phase::SweepConfig;
use fs_contact::frictional_multicontact::FrictionalContactConfig;
use fs_exec::CancelGate;
use fs_mbd::frictional_impact::RigidImpactControls;
use fs_mbd::frictional_impact::sphere_scene::*;
use fs_mbd::frictional_impact::spheres::ContactSphere;
use fs_mbd::frictional_impact::step::RigidContactStepControls;
use fs_mbd::{Gravity, MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3, Wrench};
use std::cell::Cell;

const RADIUS: f64 = 0.5;
const CLEARANCE: f64 = 1e-8;
const ACTIVATION: f64 = 1e-7;

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual {actual:.16e}, expected {expected:.16e}, tolerance {tolerance:.3e}"
    );
}

fn vector(actual: Vec3, expected: Vec3, tolerance: f64) {
    close(actual.x, expected.x, tolerance);
    close(actual.y, expected.y, tolerance);
    close(actual.z, expected.z, tolerance);
}

fn property(mass: f64) -> MassProperties {
    let inertia = 0.4 * mass * RADIUS * RADIUS;
    MassProperties::new(mass, Vec3::ZERO, Vec3::new(inertia, inertia, inertia)).unwrap()
}

fn body(position: Vec3, velocity: Vec3, mass: f64) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(position, UnitQuaternion::IDENTITY).unwrap(),
        velocity.scale(mass),
        Vec3::ZERO,
    )
    .unwrap()
}

fn sphere(id: u64) -> ContactSphere {
    ContactSphere {
        id,
        radius_m: RADIUS,
    }
}

fn plane(velocity: Vec3) -> SphereScenePlane {
    SphereScenePlane {
        id: 99,
        point_world: Vec3::ZERO,
        normal_world: Vec3::new(0.0, 0.0, 1.0),
        velocity_world: velocity,
    }
}

fn law(mu: f64, restitution: f64) -> SphereSceneLaw {
    SphereSceneLaw {
        friction_coefficient: mu,
        restitution,
    }
}

fn controls(maximum_step_s: f64) -> SphereSceneControls {
    SphereSceneControls {
        step: RigidContactStepControls {
            impact: RigidImpactControls {
                maximum_bodies: 8,
                maximum_contacts: 16,
                // Discovery's activation band must fit the native geometry gate.
                geometry_tolerance_m: ACTIVATION,
                velocity_tolerance_mps: 2e-8,
                linear_momentum_tolerance_ns: 1e-8,
                angular_momentum_tolerance_nms: 1e-8,
                energy_tolerance_j: 1e-7,
                contact: FrictionalContactConfig::default(),
            },
            drift_energy_tolerance_j: 1e-9,
            drift_angular_momentum_tolerance_nms: 1e-9,
            work_energy_tolerance_j: 2e-7,
            linear_momentum_tolerance_ns: 1e-8,
            angular_momentum_tolerance_nms: 1e-8,
        },
        broad_phase: SweepConfig {
            max_candidates: 4096,
            max_comparisons: 65_536,
        },
        maximum_step_s,
        maximum_steps: 2048,
        maximum_rejections: 2048,
        maximum_plane_queries: 100_000,
        activation_gap_m: ACTIVATION,
    }
}

fn spin_world(body: RigidBodyState) -> Vec3 {
    body.pose()
        .orientation()
        .rotate_body_to_world(body.angular_momentum_body())
}

fn partition_square_sum(report: &SphereSceneReport) -> f64 {
    report
        .steps
        .iter()
        .map(|step| step.native.duration_s * step.native.duration_s)
        .sum()
}

fn balances(report: &SphereSceneReport, duration: f64) {
    close(report.duration_s, duration, 1e-12);
    close(
        report.steps.iter().map(|step| step.native.duration_s).sum(),
        duration,
        1e-12,
    );
    assert!(report.minimum_path_gap_m >= 0.0);
    close(report.total_energy_residual_j, 0.0, 2e-7);
    vector(report.linear_momentum_residual_ns, Vec3::ZERO, 1e-8);
    vector(report.angular_momentum_residual_nms, Vec3::ZERO, 1e-8);
    for step in &report.steps {
        assert!(step.native.duration_s > 0.0);
        assert_eq!(
            step.contact_keys.len(),
            step.native.impact.contact.impulses.len()
        );
        close(step.native.total_energy_residual_j, 0.0, 2e-7);
        vector(step.native.linear_momentum_residual_ns, Vec3::ZERO, 1e-8);
        vector(step.native.angular_momentum_residual_nms, Vec3::ZERO, 1e-8);
    }
}

#[test]
fn g3_whole_scene_retains_small_impact_loss_under_a_large_tangent_boost() {
    let properties = [property(1.0)];
    let duration = 0.125;
    for boost in [0.0, 1e6] {
        let mut bodies = [body(
            Vec3::new(0.0, 0.0, RADIUS + CLEARANCE),
            Vec3::new(boost, 0.0, -0.3),
            1.0,
        )];
        let report = advance_sphere_scene(
            &mut bodies,
            &properties,
            &[sphere(10)],
            &[plane(Vec3::new(boost, 0.0, 0.0))],
            &[Wrench::ZERO],
            Gravity::ZERO,
            0.0,
            duration,
            controls(1.0 / 32.0),
            |_| Some(law(0.0, 0.5)),
            &CancelGate::new(),
        )
        .unwrap_or_else(|error| panic!("boost={boost}: {error:?}"));
        balances(&report, duration);
        vector(
            bodies[0].linear_momentum_world(),
            Vec3::new(boost, 0.0, 0.15),
            1e-12,
        );
        // The normal kinetic-energy loss is independent of tangent motion.
        close(report.unallocated_contact_loss_j, 0.03375, 1e-12);
        close(report.support_work_j, 0.0, 1e-12);
        close(report.force_impulse_work_j, 0.0, 1e-12);
        vector(report.support_impulse_ns, Vec3::new(0.0, 0.0, 0.45), 1e-12);
        assert_eq!(report.steps.len(), 4);
        if boost != 0.0 {
            // The two diagnostic energies cannot resolve the original loss
            // to the existing acceptance tolerance. Native increments can.
            assert!(
                (report.kinetic_after_j - report.kinetic_before_j
                    + report.unallocated_contact_loss_j)
                    .abs()
                    > controls(1.0 / 32.0).step.work_energy_tolerance_j
            );
        }
    }
}

#[test]
fn g3_whole_scene_retains_support_torque_beside_large_orbital_momentum() {
    let properties = [property(1.0)];
    let origin = Vec3::new(0.7, -0.4, 1e9);
    let boost = 1e6;
    let initial = body(
        origin.add(Vec3::new(0.0, 0.0, RADIUS)),
        Vec3::new(boost, 0.0, -0.3),
        1.0,
    );
    let mut bodies = [initial];
    let mut support = plane(Vec3::new(boost, 0.0, 0.0));
    support.point_world = origin;
    let duration = 0.125;
    let report = advance_sphere_scene(
        &mut bodies,
        &properties,
        &[sphere(10)],
        &[support],
        &[Wrench::ZERO],
        Gravity::ZERO,
        0.0,
        duration,
        controls(1.0 / 32.0),
        |_| Some(law(0.0, 0.0)),
        &CancelGate::new(),
    )
    .unwrap();
    balances(&report, duration);
    vector(
        bodies[0].linear_momentum_world(),
        Vec3::new(boost, 0.0, 0.0),
        1e-12,
    );
    close(report.unallocated_contact_loss_j, 0.045, 1e-12);
    vector(
        report.support_angular_impulse_nms,
        Vec3::new(-0.12, -0.21, 0.0),
        1e-12,
    );
    // The actual support couple is smaller than the absolute orbital totals'
    // resolution. It must still be retained in the whole-call balance.
    let naive_change = bodies[0]
        .pose()
        .position_world()
        .cross(bodies[0].linear_momentum_world())
        .sub(
            initial
                .pose()
                .position_world()
                .cross(initial.linear_momentum_world()),
        );
    assert!((naive_change.y - report.support_angular_impulse_nms.y).abs() > 0.01);
}

#[test]
fn whole_scene_accumulates_native_impulse_errors_before_atomic_publication() {
    // Each representable kick misses its command by about 6e-9 N s, inside
    // the unchanged per-step limit. Four same-sign misses exceed the global
    // limit even though every individual native step remains admissible.
    let mass = 2_f64.powi(26);
    let properties = [property(mass)];
    let initial = [body(Vec3::ZERO, Vec3::new(1.0, 0.0, 0.0), mass)];
    let wrench = Wrench {
        force_world: Vec3::new(0.1, 0.0, 0.0),
        torque_body: Vec3::ZERO,
    };
    let mut separate = initial;
    let mut residual = 0.0;
    for index in 0..4 {
        let report = advance_sphere_scene(
            &mut separate,
            &properties,
            &[sphere(10)],
            &[],
            &[wrench],
            Gravity::ZERO,
            f64::from(index),
            1.0,
            controls(1.0),
            |_| None,
            &CancelGate::new(),
        )
        .unwrap();
        residual += report.linear_momentum_residual_ns.x;
    }
    assert!(residual.abs() > controls(1.0).step.linear_momentum_tolerance_ns);

    let mut refused = initial;
    let failure = advance_sphere_scene(
        &mut refused,
        &properties,
        &[sphere(10)],
        &[],
        &[wrench],
        Gravity::ZERO,
        0.0,
        4.0,
        controls(1.0),
        |_| None,
        &CancelGate::new(),
    );
    assert!(matches!(
        failure,
        Err(SphereSceneError::Audit(
            "whole-call native momentum or energy"
        ))
    ));
    assert_eq!(refused, initial);
}

#[test]
fn thirty_two_spheres_share_a_moving_floor_without_a_global_dense_contact_limit() {
    const COUNT: usize = 32;
    let duration = 0.125;
    let velocity = Vec3::new(1.0, 0.0, -0.3);
    let support = plane(Vec3::new(0.25, 0.0, 0.0));
    let mut reference = [body(Vec3::new(0.0, 0.0, RADIUS + CLEARANCE), velocity, 1.0)];
    let reference_report = advance_sphere_scene(
        &mut reference,
        &[property(1.0)],
        &[sphere(10)],
        &[support],
        &[Wrench::ZERO],
        Gravity::ZERO,
        0.0,
        duration,
        controls(1.0 / 32.0),
        |_| Some(law(0.8, 0.5)),
        &CancelGate::new(),
    )
    .unwrap();
    // This is a real tangential stopping impulse with support work, rather
    // than a zero-friction case that can bypass the coupled Coulomb owner.
    assert!(reference_report.support_impulse_ns.x < -0.2);
    assert!(reference_report.support_work_j < -0.05);
    assert!(reference[0].angular_momentum_body().y > 0.1);

    let properties = vec![property(1.0); COUNT];
    let spheres: Vec<_> = (0..COUNT).map(|index| sphere(10 + index as u64)).collect();
    let initial: Vec<_> = (0..COUNT)
        .map(|index| {
            body(
                Vec3::new(2.0 * index as f64, 0.0, RADIUS + CLEARANCE),
                velocity,
                1.0,
            )
        })
        .collect();
    let mut bodies = initial.clone();
    let mut scene_controls = controls(1.0 / 32.0);
    scene_controls.step.impact.maximum_bodies = COUNT;
    scene_controls.step.impact.maximum_contacts = COUNT;
    // Raise only the scene cap. The original dense constitutive cap remains
    // unchanged for every connected physical island.
    assert_eq!(scene_controls.step.impact.contact.maximum_contacts, 21);
    assert!(COUNT > scene_controls.step.impact.contact.maximum_contacts);
    let report = advance_sphere_scene(
        &mut bodies,
        &properties,
        &spheres,
        &[support],
        &vec![Wrench::ZERO; COUNT],
        Gravity::ZERO,
        0.0,
        duration,
        scene_controls,
        |_| Some(law(0.8, 0.5)),
        &CancelGate::new(),
    )
    .unwrap();
    balances(&report, duration);
    assert_eq!(report.steps.len(), 4);
    let contact = &report.steps[0].native.impact.contact;
    assert_eq!(report.steps[0].contact_keys.len(), COUNT);
    assert_eq!(contact.impulses.len(), COUNT);
    assert_eq!(contact.islands.len(), COUNT);
    for (index, island) in contact.islands.iter().enumerate() {
        assert_eq!(island.body_indices, [index]);
        assert_eq!(island.contact_indices, [index]);
        for axis in 0..3 {
            close(
                island.impact.contact.impulses[0][axis],
                reference_report.steps[0].native.impact.contact.impulses[0][axis],
                1e-12,
            );
        }
        vector(
            bodies[index].linear_momentum_world(),
            reference[0].linear_momentum_world(),
            1e-12,
        );
        vector(spin_world(bodies[index]), spin_world(reference[0]), 1e-12);
        vector(
            bodies[index]
                .pose()
                .position_world()
                .sub(initial[index].pose().position_world()),
            reference[0]
                .pose()
                .position_world()
                .sub(Vec3::new(0.0, 0.0, RADIUS + CLEARANCE)),
            1e-12,
        );
    }
    close(
        report.unallocated_contact_loss_j,
        COUNT as f64 * reference_report.unallocated_contact_loss_j,
        1e-11,
    );
    close(
        report.support_work_j,
        COUNT as f64 * reference_report.support_work_j,
        1e-11,
    );
    vector(
        report.support_impulse_ns,
        reference_report.support_impulse_ns.scale(COUNT as f64),
        1e-11,
    );
}

#[test]
fn automatic_floor_settling_and_resting_support_do_not_create_micro_bounces() {
    let properties = [property(1.0)];
    let spheres = [sphere(10)];
    let planes = [plane(Vec3::ZERO)];
    let gravity = Gravity::new(Vec3::new(0.0, 0.0, -8.0)).unwrap();
    let initial = body(Vec3::new(0.0, 0.0, RADIUS + CLEARANCE), Vec3::ZERO, 1.0);
    let mut bodies = [initial];
    let report = advance_sphere_scene(
        &mut bodies,
        &properties,
        &spheres,
        &planes,
        &[Wrench::ZERO],
        gravity,
        0.0,
        1.0,
        controls(1.0 / 32.0),
        |_| Some(law(0.0, 0.8)),
        &CancelGate::new(),
    )
    .unwrap();
    balances(&report, 1.0);
    vector(
        bodies[0].pose().position_world(),
        initial.pose().position_world(),
        2e-9,
    );
    vector(bodies[0].linear_momentum_world(), Vec3::ZERO, 2e-8);
    vector(spin_world(bodies[0]), Vec3::ZERO, 2e-9);
    vector(report.support_impulse_ns, Vec3::new(0.0, 0.0, 8.0), 2e-8);
    // Each step kicks a resting unit mass to -8h, then support arrests it.
    // Both kick work and unallocated normal loss are 32h², not physical heat.
    let normal_splitting_loss = 32.0 * partition_square_sum(&report);
    close(report.force_impulse_work_j, normal_splitting_loss, 2e-8);
    close(
        report.unallocated_contact_loss_j,
        normal_splitting_loss,
        2e-8,
    );
    close(report.force_splitting_work_j, normal_splitting_loss, 2e-8);
    close(report.force_displacement_work_j, 0.0, 2e-8);
    close(report.kinetic_after_j, 0.0, 2e-8);

    // A genuine fall must discover the floor before an omitted drift penetrates.
    let mut falling = [body(
        Vec3::new(0.0, 0.0, RADIUS + 0.5 + CLEARANCE),
        Vec3::ZERO,
        1.0,
    )];
    let settled = advance_sphere_scene(
        &mut falling,
        &properties,
        &spheres,
        &planes,
        &[Wrench::ZERO],
        gravity,
        0.0,
        1.0,
        controls(1.0 / 64.0),
        |_| Some(law(0.0, 0.0)),
        &CancelGate::new(),
    )
    .unwrap();
    balances(&settled, 1.0);
    assert!(settled.rejected_steps > 0);
    let gap = falling[0].pose().position_world().z - RADIUS;
    assert!(
        (0.0..=2.0 * ACTIVATION).contains(&gap),
        "terminal gap {gap:e}"
    );
    vector(falling[0].linear_momentum_world(), Vec3::ZERO, 2e-8);
    vector(settled.support_impulse_ns, Vec3::new(0.0, 0.0, 8.0), 2e-8);
    close(settled.force_displacement_work_j, 4.0, 2e-6);
}

#[test]
fn automatic_stack_transfers_upper_weight_through_the_lower_sphere() {
    let properties = [property(2.0), property(1.0)];
    let spheres = [sphere(10), sphere(20)];
    let initial = [
        body(Vec3::new(0.0, 0.0, RADIUS + CLEARANCE), Vec3::ZERO, 2.0),
        body(
            Vec3::new(0.0, 0.0, 3.0 * RADIUS + 2.0 * CLEARANCE),
            Vec3::ZERO,
            1.0,
        ),
    ];
    let mut bodies = initial;
    let duration = 0.5;
    let report = advance_sphere_scene(
        &mut bodies,
        &properties,
        &spheres,
        &[plane(Vec3::ZERO)],
        &[Wrench::ZERO; 2],
        Gravity::new(Vec3::new(0.0, 0.0, -8.0)).unwrap(),
        0.0,
        duration,
        controls(1.0 / 32.0),
        |_| Some(law(0.0, 0.8)),
        &CancelGate::new(),
    )
    .unwrap();
    balances(&report, duration);
    let mut floor_load = 0.0;
    let mut internal_load = 0.0;
    for step in &report.steps {
        for (key, impulse) in step
            .contact_keys
            .iter()
            .zip(&step.native.impact.contact.impulses)
        {
            match key {
                SphereSceneContactKey::SpherePlane {
                    sphere: 10,
                    plane: 99,
                } => floor_load += impulse[0],
                SphereSceneContactKey::Spheres {
                    first: 10,
                    second: 20,
                } => internal_load += impulse[0],
                _ => close(impulse[0], 0.0, 2e-8),
            }
        }
    }
    close(floor_load, 24.0 * duration, 3e-8);
    close(internal_load, 8.0 * duration, 3e-8);
    vector(
        report.support_impulse_ns,
        Vec3::new(0.0, 0.0, 24.0 * duration),
        3e-8,
    );
    for (after, before) in bodies.iter().zip(initial) {
        vector(
            after.pose().position_world(),
            before.pose().position_world(),
            2e-9,
        );
        vector(after.linear_momentum_world(), Vec3::ZERO, 2e-8);
        vector(spin_world(*after), Vec3::ZERO, 2e-9);
    }
}

#[test]
fn discovered_incline_contacts_resolve_rolling_and_coulomb_slip() {
    // In plane-aligned coordinates g=(6,0,-8). A homogeneous sphere has
    // I=2mR²/5: rolling gives a=5gx/7 and requires mu>=3/14.
    for (mu, acceleration, angular_acceleration, slip_acceleration) in
        [(0.3, 30.0 / 7.0, 60.0 / 7.0, 0.0), (0.1, 5.2, 4.0, 3.2)]
    {
        let properties = [property(1.0)];
        let mut bodies = [body(
            Vec3::new(0.0, 0.0, RADIUS + CLEARANCE),
            Vec3::ZERO,
            1.0,
        )];
        let duration = 0.5;
        let report = advance_sphere_scene(
            &mut bodies,
            &properties,
            &[sphere(10)],
            &[plane(Vec3::ZERO)],
            &[Wrench::ZERO],
            Gravity::new(Vec3::new(6.0, 0.0, -8.0)).unwrap(),
            0.0,
            duration,
            controls(1.0 / 32.0),
            |_| Some(law(mu, 0.0)),
            &CancelGate::new(),
        )
        .unwrap();
        balances(&report, duration);
        let velocity = bodies[0]
            .center_of_mass_velocity_world(properties[0])
            .unwrap();
        let omega = spin_world(bodies[0]).scale(10.0);
        vector(velocity, Vec3::new(acceleration * duration, 0.0, 0.0), 2e-7);
        vector(
            omega,
            Vec3::new(0.0, angular_acceleration * duration, 0.0),
            2e-7,
        );
        close(
            velocity.x - RADIUS * omega.y,
            slip_acceleration * duration,
            2e-7,
        );
        // Exact endpoint for this declared first-order partition, not a
        // second-order continuum trajectory incorrectly claimed as exact.
        close(
            bodies[0].pose().position_world().x,
            0.5 * acceleration * (duration * duration + partition_square_sum(&report)),
            2e-7,
        );
        close(
            bodies[0].pose().position_world().z,
            RADIUS + CLEARANCE,
            2e-9,
        );
    }
}

#[test]
fn discovered_moving_conveyor_transfers_work_and_spin_to_a_finite_sphere() {
    let properties = [property(1.0)];
    let mut bodies = [body(
        Vec3::new(0.0, 0.0, RADIUS + CLEARANCE),
        Vec3::ZERO,
        1.0,
    )];
    let report = advance_sphere_scene(
        &mut bodies,
        &properties,
        &[sphere(10)],
        &[plane(Vec3::new(1.0, 0.0, 0.0))],
        &[Wrench::ZERO],
        Gravity::new(Vec3::new(0.0, 0.0, -8.0)).unwrap(),
        0.0,
        0.25,
        controls(1.0 / 32.0),
        |_| Some(law(0.25, 0.0)),
        &CancelGate::new(),
    )
    .unwrap();
    balances(&report, 0.25);
    // W_t=1/m+R²/I=7/2: the impulse arresting one m/s relative
    // surface slip is 2/7 Ns. The conveyor performs J*U=2/7 J.
    vector(
        bodies[0].linear_momentum_world(),
        Vec3::new(2.0 / 7.0, 0.0, 0.0),
        2e-8,
    );
    vector(spin_world(bodies[0]), Vec3::new(0.0, -1.0 / 7.0, 0.0), 2e-8);
    vector(
        report.support_impulse_ns,
        Vec3::new(2.0 / 7.0, 0.0, 2.0),
        2e-8,
    );
    close(report.support_work_j, 2.0 / 7.0, 2e-8);
    close(report.kinetic_after_j, 1.0 / 7.0, 2e-8);
    close(
        bodies[0].kinetic_energy(properties[0]).unwrap(),
        1.0 / 7.0,
        2e-8,
    );
    close(
        bodies[0].pose().position_world().z,
        RADIUS + CLEARANCE,
        2e-9,
    );
    close(
        report.unallocated_contact_loss_j,
        1.0 / 7.0 + 32.0 * partition_square_sum(&report),
        3e-8,
    );
}

#[test]
fn automatic_scene_retries_hidden_crossings_and_refuses_exhausted_retry_work() {
    let properties = [property(1.0); 2];
    let spheres = [sphere(10), sphere(20)];
    let initial = [
        body(Vec3::new(-1.5, 0.0, 0.0), Vec3::new(4.0, 0.0, 0.0), 1.0),
        body(Vec3::new(1.5, 0.0, 0.0), Vec3::new(-4.0, 0.0, 0.0), 1.0),
    ];
    let mut bodies = initial;
    let report = advance_sphere_scene(
        &mut bodies,
        &properties,
        &spheres,
        &[],
        &[Wrench::ZERO; 2],
        Gravity::ZERO,
        0.0,
        0.75,
        controls(0.75),
        |_| Some(law(0.0, 1.0)),
        &CancelGate::new(),
    )
    .unwrap();
    balances(&report, 0.75);
    assert!(report.rejected_steps > 0);
    // Uncollided endpoints are separated again after passing through. The
    // shared CCD hint should reach the activation band at t=1/4, rather than
    // first accepting the 3/16 step produced by repeated blind halving.
    close(report.steps[0].native.duration_s, 0.25, ACTIVATION);
    vector(
        bodies[0].pose().position_world(),
        Vec3::new(-2.5, 0.0, 0.0),
        2.0 * ACTIVATION,
    );
    vector(
        bodies[1].pose().position_world(),
        Vec3::new(2.5, 0.0, 0.0),
        2.0 * ACTIVATION,
    );
    vector(
        bodies[0].linear_momentum_world(),
        Vec3::new(-4.0, 0.0, 0.0),
        2e-8,
    );
    vector(
        bodies[1].linear_momentum_world(),
        Vec3::new(4.0, 0.0, 0.0),
        2e-8,
    );
    close(report.kinetic_before_j, 16.0, 2e-8);
    close(report.kinetic_after_j, 16.0, 2e-8);
    vector(report.support_impulse_ns, Vec3::ZERO, 2e-8);

    // The first 3/16 s drift is clear and changes both staged positions.
    // The next drift crosses. A zero retry budget must roll back the entire
    // call, including that already accepted first step, not just the trial.
    let mut limited = controls(3.0 / 16.0);
    limited.maximum_rejections = 0;
    let mut refused = initial;
    let failure = advance_sphere_scene(
        &mut refused,
        &properties,
        &spheres,
        &[],
        &[Wrench::ZERO; 2],
        Gravity::ZERO,
        0.0,
        0.75,
        limited,
        |_| Some(law(0.0, 1.0)),
        &CancelGate::new(),
    );
    assert!(matches!(failure, Err(SphereSceneError::RejectionBudget)));
    assert_eq!(refused, initial);
}

#[test]
fn plane_query_budget_is_whole_call_and_callback_cancellation_is_atomic() {
    let properties = [property(1.0)];
    let spheres = [sphere(10)];
    let planes = [plane(Vec3::ZERO)];
    let initial = [body(
        Vec3::new(0.0, 0.0, RADIUS + CLEARANCE),
        Vec3::new(1.0, 0.0, 0.0),
        1.0,
    )];
    let mut bodies = initial;
    let report = advance_sphere_scene(
        &mut bodies,
        &properties,
        &spheres,
        &planes,
        &[Wrench::ZERO],
        Gravity::ZERO,
        0.0,
        0.5,
        controls(0.125),
        |_| Some(law(0.0, 0.0)),
        &CancelGate::new(),
    )
    .unwrap();
    balances(&report, 0.5);
    close(bodies[0].pose().position_world().x, 0.5, 2e-9);
    assert!(report.plane_queries > 1);

    // Use the measured whole-call count, not an assumed internal query count.
    // Resetting the budget at each step would incorrectly admit this call.
    let mut limited = controls(0.125);
    limited.maximum_plane_queries = report.plane_queries - 1;
    let mut refused = initial;
    let failure = advance_sphere_scene(
        &mut refused,
        &properties,
        &spheres,
        &planes,
        &[Wrench::ZERO],
        Gravity::ZERO,
        0.0,
        0.5,
        limited,
        |_| Some(law(0.0, 0.0)),
        &CancelGate::new(),
    );
    assert!(failure.is_err());
    assert_eq!(refused, initial);

    let gate = CancelGate::new();
    let callback_ran = Cell::new(false);
    let mut cancelled = initial;
    let failure = advance_sphere_scene(
        &mut cancelled,
        &properties,
        &spheres,
        &planes,
        &[Wrench::ZERO],
        Gravity::ZERO,
        0.0,
        0.5,
        controls(0.125),
        |_| {
            callback_ran.set(true);
            gate.request();
            Some(law(0.0, 0.0))
        },
        &gate,
    );
    assert!(callback_ran.get());
    assert!(gate.is_requested());
    assert!(failure.is_err());
    assert_eq!(cancelled, initial);
}
