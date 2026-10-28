//! Native rigid contact trajectories, physical force work and geometric admission.
//! Synthetic cuboid/support witnesses exercise the actual rotational mobility;
//! their center of pressure is supplied analytically, not discovered by the solver.
use fs_contact::frictional_multicontact::{FrictionalContactConfig, NormalVelocityLaw};
use fs_exec::CancelGate;
use fs_mbd::frictional_impact::step::*;
use fs_mbd::frictional_impact::*;
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};
use std::cell::Cell;

const HALF: f64 = 0.25;
const INERTIA: f64 = 0.0625;
// A declared positive clearance inside the original contact activation band
// avoids pretending floating-point zero-gap certification or projection.
const CLEARANCE: f64 = 1e-8;

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}; tolerance {tolerance:e}"
    );
}
fn close_vec(a: Vec3, b: Vec3, tolerance: f64) {
    close(a.x, b.x, tolerance);
    close(a.y, b.y, tolerance);
    close(a.z, b.z, tolerance);
}
fn controls() -> RigidContactStepControls {
    RigidContactStepControls {
        impact: RigidImpactControls {
            maximum_bodies: 4,
            maximum_contacts: 4,
            geometry_tolerance_m: 1e-10,
            velocity_tolerance_mps: 2e-8,
            linear_momentum_tolerance_ns: 1e-9,
            angular_momentum_tolerance_nms: 1e-9,
            energy_tolerance_j: 2e-8,
            contact: FrictionalContactConfig::default(),
        },
        drift_energy_tolerance_j: 1e-10,
        drift_angular_momentum_tolerance_nms: 1e-10,
        work_energy_tolerance_j: 3e-8,
        linear_momentum_tolerance_ns: 1e-9,
        angular_momentum_tolerance_nms: 1e-9,
    }
}
fn property() -> MassProperties {
    MassProperties::new(1.0, Vec3::ZERO, Vec3::new(INERTIA, INERTIA, INERTIA)).unwrap()
}
fn block() -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(
            Vec3::new(0.0, 0.0, HALF + CLEARANCE),
            UnitQuaternion::IDENTITY,
        )
        .unwrap(),
        Vec3::ZERO,
        Vec3::ZERO,
    )
    .unwrap()
}
fn plane(gx: f64) -> RigidFrictionalContact {
    // For sticking the reaction balances gravity through the COM; for sliding
    // the limiting Coulomb reaction does. The finite support point lies within
    // the original cuboid footprint, and its actual lever arm enters mobility.
    let pressure_x = HALF * (gx / 8.0).min(0.5);
    RigidFrictionalContact {
        id: 1,
        point: RigidImpactPoint {
            body: 0,
            arm_body: Vec3::new(pressure_x, 0.0, -HALF),
        },
        other: RigidImpactCounterpart::PrescribedPlane {
            point_world: Vec3::ZERO,
            velocity_world: Vec3::ZERO,
        },
        normal_world: Vec3::new(0.0, 0.0, 1.0),
        tangent_world: Vec3::new(1.0, 0.0, 0.0),
        activation_tolerance_m: 1e-7,
        friction_coefficient: 0.5,
        normal_law: NormalVelocityLaw::NewtonRestitution { coefficient: 0.0 },
    }
}
fn gravity(gx: f64) -> Gravity {
    Gravity::new(Vec3::new(gx, 0.0, -8.0)).unwrap()
}

// This fixture's complete collision scene is one known cuboid and one infinite,
// translating plane. The native spherical top follows an exact constant-spin
// rotation. For EVERY time, each vertex's rotational displacement is bounded by
// R*|omega|*dt and its centre-plane height is affine. This is an interval-wide
// geometric bound, not acceptance inferred from a few nonnegative samples.
fn admit_box(path: &RigidContactPath<'_>) -> Result<f64, &'static str> {
    if path.initial.len() != 1 || path.contacts.len() != 1 {
        return Err("unexpected collision scene");
    }
    let c = path.contacts[0];
    if c.normal_world != Vec3::new(0.0, 0.0, 1.0) {
        return Err("unsupported fixture plane normal");
    }
    let (plane, velocity) = match c.other {
        RigidImpactCounterpart::PrescribedPlane {
            point_world,
            velocity_world,
        } => (point_world, velocity_world),
        _ => return Err("fixture needs its complete plane"),
    };
    let props = path.properties[0];
    if props.principal_inertia_body() != Vec3::new(INERTIA, INERTIA, INERTIA) {
        return Err("rotation bound needs original spherical inertia");
    }
    let q = path.post_contact[0].pose().orientation();
    let mut lowest = f64::INFINITY;
    for x in [-HALF, HALF] {
        for y in [-HALF, HALF] {
            for z in [-HALF, HALF] {
                lowest = lowest.min(q.rotate_body_to_world(Vec3::new(x, y, z)).z);
            }
        }
    }
    let initial_z = path.post_contact[0].pose().position_world().z - plane.z;
    let final_z =
        path.endpoint[0].pose().position_world().z - plane.z - velocity.z * path.duration_s;
    let omega = props
        .angular_velocity_body_checked(path.post_contact[0].angular_momentum_body())
        .map_err(|_| "angular speed")?;
    let rotation_bound = 3.0_f64.sqrt() * HALF * omega.norm_squared().sqrt() * path.duration_s;
    let lower = initial_z.min(final_z) + lowest - rotation_bound;
    if !lower.is_finite() || lower < 0.0 {
        return Err("whole cuboid path intersects plane");
    }
    Ok(lower)
}

#[test]
fn finite_height_block_obeys_incline_threshold_and_sustained_support() {
    // Rotating coordinates to align the plane leaves tan(theta)=gx/|gz|.
    // The original native orientation, nonzero COM arm and inertia are retained.
    for gx in [2.0, 4.0, 6.0] {
        let mut bodies = [block()];
        let properties = [property()];
        let contact = plane(gx);
        let mut accumulated_loss = 0.0;
        let mut accumulated_kick = 0.0;
        let initial = bodies[0];
        for _ in 0..32 {
            let result = step_rigid_frictional_contacts(
                &mut bodies,
                &properties,
                &[Wrench::ZERO],
                gravity(gx),
                &[contact],
                1.0 / 32.0,
                controls(),
                admit_box,
                &CancelGate::new(),
            )
            .unwrap();
            assert!(result.endpoint_gaps[0].signed_gap_m >= 0.0);
            assert!(result.path_admission >= 0.0);
            close(result.impact.contact.impulses[0][0], 8.0 / 32.0, 2e-8);
            close(result.total_energy_residual_j, 0.0, 3e-8);
            accumulated_loss += result.impact.contact.dissipated_energy;
            accumulated_kick += result.force_impulse_work_j;
        }
        let acceleration = (gx - 0.5 * 8.0).max(0.0);
        close(bodies[0].linear_momentum_world().x, acceleration, 2e-7);
        close(
            bodies[0].pose().position_world().x,
            0.5 * acceleration * (1.0 + 1.0 / 32.0),
            2e-7,
        );
        close(
            bodies[0].pose().position_world().z,
            initial.pose().position_world().z,
            2e-10,
        );
        close_vec(bodies[0].angular_momentum_body(), Vec3::ZERO, 2e-9);
        close(
            bodies[0].kinetic_energy(properties[0]).unwrap(),
            accumulated_kick - accumulated_loss,
            1e-7,
        );
        if gx <= 4.0 {
            close_vec(bodies[0].linear_momentum_world(), Vec3::ZERO, 2e-8);
            close_vec(
                bodies[0].pose().position_world(),
                initial.pose().position_world(),
                2e-9,
            );
        }
    }
}

#[test]
fn sliding_displacement_force_work_and_event_loss_converge_at_first_order() {
    let mut displacement_errors = Vec::new();
    let mut splitting_errors = Vec::new();
    let mut loss_errors = Vec::new();
    for count in [8usize, 16, 32] {
        let dt = 1.0 / count as f64;
        let mut bodies = [block()];
        let properties = [property()];
        let mut kick = 0.0;
        let mut force = 0.0;
        let mut defect = 0.0;
        let mut loss = 0.0;
        for _ in 0..count {
            let result = step_rigid_frictional_contacts(
                &mut bodies,
                &properties,
                &[Wrench::ZERO],
                gravity(6.0),
                &[plane(6.0)],
                dt,
                controls(),
                admit_box,
                &CancelGate::new(),
            )
            .unwrap();
            kick += result.force_impulse_work_j;
            force += result.force_displacement_work_j;
            defect += result.force_splitting_work_j;
            loss += result.impact.contact.dissipated_energy;
        }
        // Exact continuous sliding: a=6-mu*8=2, x(1)=1, v(1)=2,
        // gravity force work 6 J and Coulomb work 4 J. Kick/contact/drift
        // is first order: x=1+dt, kick work=6+44dt, event loss=4+44dt.
        // The event loss includes normal kick loss and MUST NOT all heat a body.
        let x = bodies[0].pose().position_world().x;
        close(x, 1.0 + dt, 2e-7);
        close(bodies[0].linear_momentum_world().x, 2.0, 2e-7);
        close(kick, 6.0 + 44.0 * dt, 3e-6);
        close(loss, 4.0 + 44.0 * dt, 3e-6);
        close(force, 6.0 * (1.0 + dt), 2e-6);
        close(defect, kick - force, 2e-7);
        displacement_errors.push((x - 1.0).abs());
        splitting_errors.push(defect.abs());
        loss_errors.push((loss - 4.0).abs());
    }
    for errors in [displacement_errors, splitting_errors, loss_errors] {
        for pair in errors.windows(2) {
            assert!(
                pair[1] < 0.51 * pair[0] && pair[1] > 0.49 * pair[0],
                "{errors:?}"
            );
        }
    }
}

#[test]
fn actual_midpoint_free_rotation_and_prescribed_plane_motion_are_retained() {
    let props = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.5, 2.0)).unwrap();
    let initial = RigidBodyState::new(
        Pose::new(
            Vec3::new(1.0, 2.0, 3.0),
            UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 3.0), 0.7).unwrap(),
        )
        .unwrap(),
        Vec3::new(0.3, -0.2, 0.4),
        Vec3::new(0.4, 0.7, 1.1),
    )
    .unwrap();
    let expected = RigidBodyIntegrator::new(Gravity::ZERO)
        .step(initial, props, Wrench::ZERO, 0.125)
        .unwrap()
        .state_after;
    let mut bodies = [initial];
    let report = step_rigid_frictional_contacts(
        &mut bodies,
        &[props],
        &[Wrench::ZERO],
        Gravity::ZERO,
        &[],
        0.125,
        controls(),
        |path| {
            assert!(path.contacts.is_empty()); // The fixture scene contains no collision geometry.
            assert_eq!(path.sample(0.0).unwrap(), path.post_contact);
            assert_eq!(path.sample(1.0).unwrap(), path.endpoint);
            let half = RigidBodyIntegrator::new(Gravity::ZERO)
                .step(initial, props, Wrench::ZERO, 0.0625)
                .unwrap()
                .state_after;
            assert_eq!(path.sample(0.5).unwrap(), vec![half]);
            Ok::<_, &'static str>(())
        },
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(bodies[0], expected);
    close(report.drift_energy_defect_j, 0.0, 1e-11);
    close_vec(report.drift_angular_momentum_defect_nms, Vec3::ZERO, 1e-11);

    let mut bodies =
        [RigidBodyState::new(block().pose(), Vec3::new(0.0, 0.0, -1.0), Vec3::ZERO).unwrap()];
    let mut c = plane(0.0);
    c.friction_coefficient = 0.0;
    c.other = RigidImpactCounterpart::PrescribedPlane {
        point_world: Vec3::ZERO,
        velocity_world: Vec3::new(0.0, 0.0, 0.5),
    };
    let result = step_rigid_frictional_contacts(
        &mut bodies,
        &[property()],
        &[Wrench::ZERO],
        Gravity::ZERO,
        &[c],
        0.125,
        controls(),
        admit_box,
        &CancelGate::new(),
    )
    .unwrap();
    close(bodies[0].linear_momentum_world().z, 0.5, 1e-10);
    close(result.endpoint_gaps[0].signed_gap_m, CLEARANCE, 1e-10);
    close(result.impact.support_work_j, 0.75, 1e-9);
    close(result.impact.contact.dissipated_energy, 1.125, 1e-9);
    close(result.total_energy_residual_j, 0.0, 1e-9);
}

#[test]
fn endpoint_penetration_unresolved_contact_and_whole_path_refusals_are_atomic() {
    let properties = [property()];
    let mut bodies = [block()];
    let original = bodies;
    let contacted = Cell::new(false);
    let failure = step_rigid_frictional_contacts(
        &mut bodies,
        &properties,
        &[Wrench::ZERO],
        gravity(6.0),
        &[plane(6.0)],
        0.125,
        controls(),
        |path| {
            contacted.set(true);
            assert!(path.endpoint[0].pose().position_world().x > 0.0);
            Err::<(), _>("unmodeled feature intersects the complete path")
        },
        &CancelGate::new(),
    );
    assert!(matches!(failure, Err(RigidContactStepError::Path(_))));
    assert!(contacted.get());
    assert_eq!(bodies, original);
    let gate = CancelGate::new();
    let failure = step_rigid_frictional_contacts(
        &mut bodies,
        &properties,
        &[Wrench::ZERO],
        gravity(6.0),
        &[plane(6.0)],
        0.125,
        controls(),
        |path| {
            admit_box(path)?;
            gate.request();
            Ok::<_, &'static str>(())
        },
        &gate,
    );
    assert!(matches!(failure, Err(RigidContactStepError::Cancelled)));
    assert_eq!(bodies, original);

    let mut distant = [RigidBodyState::new(
        Pose::new(Vec3::new(0.0, 0.0, 1.0), UnitQuaternion::IDENTITY).unwrap(),
        Vec3::new(0.0, 0.0, -2.0),
        Vec3::ZERO,
    )
    .unwrap()];
    let old = distant;
    let called = Cell::new(false);
    let failure = step_rigid_frictional_contacts(
        &mut distant,
        &properties,
        &[Wrench::ZERO],
        Gravity::ZERO,
        &[plane(0.0)],
        1.0,
        controls(),
        |_| {
            called.set(true);
            Ok::<_, &'static str>(())
        },
        &CancelGate::new(),
    );
    assert!(matches!(
        failure,
        Err(RigidContactStepError::Penetration { endpoint: true, .. })
    ));
    assert!(!called.get());
    assert_eq!(distant, old);
    let c = plane(6.0);
    let mut redundant = c;
    redundant.id = 2;
    // Identical redundant laws can now have determined motion. Different
    // friction coefficients at an unresolved shared normal reaction can give
    // different native motions, so this remains a genuine determinacy refusal.
    redundant.friction_coefficient = 0.1;
    let failure = step_rigid_frictional_contacts(
        &mut bodies,
        &properties,
        &[Wrench::ZERO],
        gravity(6.0),
        &[c, redundant],
        0.125,
        controls(),
        |_| Ok::<_, &'static str>(()),
        &CancelGate::new(),
    );
    assert!(matches!(failure, Err(RigidContactStepError::Impact(_))));
    assert_eq!(bodies, original);
    let mut penetrating = [RigidBodyState::new(
        Pose::new(Vec3::new(0.0, 0.0, HALF - 1e-4), UnitQuaternion::IDENTITY).unwrap(),
        Vec3::ZERO,
        Vec3::ZERO,
    )
    .unwrap()];
    let old = penetrating;
    let failure = step_rigid_frictional_contacts(
        &mut penetrating,
        &properties,
        &[Wrench::ZERO],
        gravity(6.0),
        &[plane(6.0)],
        0.125,
        controls(),
        |_| Ok::<_, &'static str>(()),
        &CancelGate::new(),
    );
    assert!(matches!(
        failure,
        Err(RigidContactStepError::Penetration {
            endpoint: false,
            ..
        })
    ));
    assert_eq!(penetrating, old);
    let retry = step_rigid_frictional_contacts(
        &mut bodies,
        &properties,
        &[Wrench::ZERO],
        gravity(6.0),
        &[plane(6.0)],
        0.125,
        controls(),
        admit_box,
        &CancelGate::new(),
    )
    .unwrap();
    let mut clean = original;
    let fresh = step_rigid_frictional_contacts(
        &mut clean,
        &properties,
        &[Wrench::ZERO],
        gravity(6.0),
        &[plane(6.0)],
        0.125,
        controls(),
        admit_box,
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(bodies, clean);
    assert_eq!(retry.impact.contact, fresh.impact.contact);
}
