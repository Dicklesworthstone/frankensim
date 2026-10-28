//! G1/G3/G4: automatic curved collision paths with original independent loads.
use super::*;

fn spheres(radius: f64) -> [ContactSphere; 2] {
    [
        ContactSphere {
            id: 10,
            radius_m: radius,
        },
        ContactSphere {
            id: 20,
            radius_m: radius,
        },
    ]
}

#[test]
fn g1_independent_force_collision_closes_gravity_and_orbital_force_work() {
    let spheres = spheres(0.5);
    let properties = [property(1.0, 0.5); 2];
    let mut bodies = [
        body(
            Vec3::new(3.0, 1.0, 0.0),
            Vec3::new(-4.0, 0.0, 0.0),
            1.0,
            Vec3::ZERO,
        ),
        body(Vec3::new(0.0, 1.0, 0.0), Vec3::ZERO, 1.0, Vec3::ZERO),
    ];
    let report = advance_forced_sphere_contacts(
        &mut bodies,
        &properties,
        &spheres,
        &laws(&spheres, 0.0, 1.0),
        Gravity::new(Vec3::new(0.0, -1.0, 0.0)).unwrap(),
        &[Vec3::new(2.0, 0.0, 0.0), Vec3::ZERO],
        1.0,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let root = 2.0_f64.sqrt();
    assert_eq!(report.events.len(), 1);
    close(report.events[0].elapsed_s, 2.0 - root);
    vector(
        bodies[0].pose().position_world(),
        Vec3::new(4.0 - 2.0 * root, 0.5, 0.0),
    );
    vector(
        bodies[1].pose().position_world(),
        Vec3::new(-4.0 + 2.0 * root, 0.5, 0.0),
    );
    vector(
        bodies[0].linear_momentum_world(),
        Vec3::new(2.0 * root - 2.0, -1.0, 0.0),
    );
    vector(
        bodies[1].linear_momentum_world(),
        Vec3::new(-2.0 * root, -1.0, 0.0),
    );
    close(report.gravity_work_j, 1.0);
    close(report.centre_force_work_j, 2.0 - 4.0 * root);
    vector(report.gravity_impulse_ns, Vec3::new(0.0, -2.0, 0.0));
    vector(report.centre_force_impulse_ns, Vec3::new(2.0, 0.0, 0.0));
    // The gravity torque includes the extra centre-force acceleration of the
    // total mass moment. The centre force sees y(t)=1-t²/2 across the collision.
    vector(
        report.gravity_angular_impulse_nms,
        Vec3::new(0.0, 0.0, -4.0 / 3.0),
    );
    vector(
        report.centre_force_angular_impulse_nms,
        Vec3::new(0.0, 0.0, -5.0 / 3.0),
    );
    close(report.unallocated_impact_loss_j, 0.0);
    close(report.total_energy_residual_j, 0.0);
    vector(report.angular_momentum_residual_nms, Vec3::ZERO);
    assert!(report.minimum_path_gap_m >= 0.0);
}

#[test]
fn g1_curved_oblique_entry_uses_original_three_dimensional_coulomb_impulse() {
    // Independently manufactured relative path r=(2-2t,t-t²,3/5).
    // At t=3/5, r=(4/5,6/25,3/5), defining the exact radius sum.
    let radius_sum = (1.0_f64 + 36.0 / 625.0).sqrt();
    let radius = 0.5 * radius_sum;
    let spheres = spheres(radius);
    let masses = [2.0, 3.0];
    let properties = [property(masses[0], radius), property(masses[1], radius)];
    let initial = [
        body(
            Vec3::new(2.0, 0.0, 0.6),
            Vec3::new(-2.0, 1.0, 0.0),
            2.0,
            Vec3::ZERO,
        ),
        body(Vec3::ZERO, Vec3::ZERO, 3.0, Vec3::ZERO),
    ];
    let mut bodies = initial;
    let report = advance_forced_sphere_contacts(
        &mut bodies,
        &properties,
        &spheres,
        &laws(&spheres, 0.2, 0.5),
        Gravity::ZERO,
        &[Vec3::new(0.0, -4.0, 0.0), Vec3::ZERO],
        0.75,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(report.events.len(), 1);
    close(report.events[0].elapsed_s, 0.6);
    let n = Vec3::new(0.8, 0.24, 0.6).scale(radius_sum.recip());
    let incoming = Vec3::new(-2.0, -0.2, 0.0);
    let normal_speed = incoming.dot(n);
    let tangent = incoming.sub(n.scale(normal_speed));
    let inv_mass = 1.0 / 2.0 + 1.0 / 3.0;
    let normal_impulse = -1.5 * normal_speed / inv_mass;
    // For each original solid sphere R²/I=5/(2m), so its actual tangential
    // contact mobility is 7/(2m), independent of the sphere's radius.
    let tangent_mobility = 3.5 * inv_mass;
    let tangent_speed = tangent.norm_squared().sqrt();
    let magnitude = (tangent_speed / tangent_mobility).min(0.2 * normal_impulse);
    let jt = tangent.scale(-magnitude / tangent_speed);
    let impulse = n.scale(normal_impulse).add(jt);
    let tail = 0.15;
    let force = Vec3::new(0.0, -4.0, 0.0);
    let velocity_a = incoming.add(impulse.scale(0.5));
    let velocity_b = impulse.scale(-1.0 / 3.0);
    let at_contact = Vec3::new(0.8, 0.24, 0.6);
    let end_a = at_contact
        .add(velocity_a.scale(tail))
        .add(force.scale(0.25 * tail * tail));
    vector(
        bodies[0].linear_momentum_world(),
        velocity_a.scale(2.0).add(force.scale(tail)),
    );
    vector(bodies[1].linear_momentum_world(), impulse.scale(-1.0));
    vector(spin_world(bodies[0]), n.scale(-radius).cross(jt));
    vector(spin_world(bodies[1]), n.scale(-radius).cross(jt));
    vector(bodies[0].pose().position_world(), end_a);
    vector(bodies[1].pose().position_world(), velocity_b.scale(tail));
    let loss = -(normal_impulse * normal_speed
        + 0.5 * inv_mass * normal_impulse * normal_impulse
        + jt.dot(tangent)
        + 0.5 * tangent_mobility * jt.norm_squared());
    close(report.unallocated_impact_loss_j, loss);
    close(report.centre_force_work_j, -4.0 * end_a.y);
    vector(
        report.centre_force_angular_impulse_nms,
        Vec3::new(1.44, 0.0, -3.36)
            .add(at_contact.cross(force).scale(tail))
            .add(velocity_a.cross(force).scale(0.5 * tail * tail)),
    );
    close(report.total_energy_residual_j, 0.0);
    vector(report.angular_momentum_residual_nms, Vec3::ZERO);
}

#[test]
fn g3_accelerated_return_permutation_and_split_continuation_preserve_motion() {
    let spheres = spheres(0.5);
    let properties = [property(1.0, 0.5); 2];
    let initial = [
        body(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
            1.0,
            Vec3::ZERO,
        ),
        body(Vec3::ZERO, Vec3::ZERO, 1.0, Vec3::ZERO),
    ];
    let forces = [Vec3::new(-2.0, 0.0, 0.0), Vec3::ZERO];
    let law = laws(&spheres, 0.2, 1.0);
    let mut full = initial;
    let report = advance_forced_sphere_contacts(
        &mut full,
        &properties,
        &spheres,
        &law,
        Gravity::ZERO,
        &forces,
        2.5,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(report.events.len(), 1);
    close(report.events[0].elapsed_s, 2.0);
    vector(full[0].pose().position_world(), Vec3::new(0.75, 0.0, 0.0));
    vector(full[1].pose().position_world(), Vec3::new(-1.0, 0.0, 0.0));
    let mut permuted = [initial[1], initial[0]];
    let reverse = advance_forced_sphere_contacts(
        &mut permuted,
        &properties,
        &[spheres[1], spheres[0]],
        &law,
        Gravity::ZERO,
        &[forces[1], forces[0]],
        2.5,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(permuted, [full[1], full[0]]);
    assert_eq!(reverse.events[0].pairs, report.events[0].pairs);
    let mut split = initial;
    let a = advance_forced_sphere_contacts(
        &mut split,
        &properties,
        &spheres,
        &law,
        Gravity::ZERO,
        &forces,
        1.0,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let b = advance_forced_sphere_contacts(
        &mut split,
        &properties,
        &spheres,
        &law,
        Gravity::ZERO,
        &forces,
        1.5,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    for i in 0..2 {
        vector(
            split[i].pose().position_world(),
            full[i].pose().position_world(),
        );
        vector(
            split[i].linear_momentum_world(),
            full[i].linear_momentum_world(),
        );
    }
    close(
        a.centre_force_work_j + b.centre_force_work_j,
        report.centre_force_work_j,
    );
    close(b.events[0].elapsed_s + 1.0, report.events[0].elapsed_s);
}

#[test]
fn g4_loaded_support_and_late_event_exhaustion_are_atomic() {
    let spheres = spheres(0.5);
    let properties = [property(1.0, 0.5); 2];
    let initial = [
        body(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
            1.0,
            Vec3::ZERO,
        ),
        body(Vec3::ZERO, Vec3::ZERO, 1.0, Vec3::ZERO),
    ];
    let forces = [Vec3::new(-2.0, 0.0, 0.0), Vec3::ZERO];
    let law = laws(&spheres, 0.2, 1.0);
    let mut bodies = initial;
    let mut limited = controls();
    limited.maximum_events = 1;
    assert!(matches!(
        advance_forced_sphere_contacts(
            &mut bodies,
            &properties,
            &spheres,
            &law,
            Gravity::ZERO,
            &forces,
            4.25,
            limited,
            &CancelGate::new()
        ),
        Err(SphereContactError::EventBudget)
    ));
    assert_eq!(bodies, initial);
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        advance_forced_sphere_contacts(
            &mut bodies,
            &properties,
            &spheres,
            &law,
            Gravity::ZERO,
            &forces,
            4.25,
            controls(),
            &gate
        ),
        Err(SphereContactError::Cancelled)
    ));
    assert_eq!(bodies, initial);
    let mut retry = initial;
    advance_forced_sphere_contacts(
        &mut bodies,
        &properties,
        &spheres,
        &law,
        Gravity::ZERO,
        &forces,
        4.25,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    advance_forced_sphere_contacts(
        &mut retry,
        &properties,
        &spheres,
        &law,
        Gravity::ZERO,
        &forces,
        4.25,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(bodies, retry);
    // Resting under a closing differential acceleration needs a sustained
    // normal reaction. This impact-only owner must not integrate through it.
    let resting = [
        body(Vec3::new(1.0, 0.0, 0.0), Vec3::ZERO, 1.0, Vec3::ZERO),
        initial[1],
    ];
    let mut bodies = resting;
    assert!(matches!(
        advance_forced_sphere_contacts(
            &mut bodies,
            &properties,
            &spheres,
            &law,
            Gravity::ZERO,
            &forces,
            0.5,
            controls(),
            &CancelGate::new()
        ),
        Err(SphereContactError::Unresolved(_))
    ));
    assert_eq!(bodies, resting);
    assert!(matches!(
        advance_forced_sphere_contacts(
            &mut bodies,
            &properties,
            &spheres,
            &law,
            Gravity::ZERO,
            &forces[..1],
            0.5,
            controls(),
            &CancelGate::new()
        ),
        Err(SphereContactError::Invalid(_))
    ));
    assert_eq!(bodies, resting);
}
