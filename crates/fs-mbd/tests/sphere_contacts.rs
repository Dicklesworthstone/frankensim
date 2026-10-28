//! G1/G3/G4: automatic native sphere trajectories and independent event oracles.
use fs_contact::frictional_multicontact::FrictionalContactConfig;
use fs_exec::CancelGate;
use fs_mbd::frictional_impact::RigidImpactControls;
use fs_mbd::frictional_impact::spheres::*;
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 3e-8 * a.abs().max(b.abs()).max(1.0),
        "{a:.16e} != {b:.16e}"
    );
}
fn vector(a: Vec3, b: Vec3) {
    close(a.x, b.x);
    close(a.y, b.y);
    close(a.z, b.z);
}
fn property(m: f64, r: f64) -> MassProperties {
    MassProperties::new(
        m,
        Vec3::ZERO,
        Vec3::new(0.4 * m * r * r, 0.4 * m * r * r, 0.4 * m * r * r),
    )
    .unwrap()
}
fn body(p: Vec3, velocity: Vec3, mass: f64, spin: Vec3) -> RigidBodyState {
    RigidBodyState::new(
        Pose::new(p, UnitQuaternion::IDENTITY).unwrap(),
        velocity.scale(mass),
        spin,
    )
    .unwrap()
}
fn controls() -> SphereContactControls {
    SphereContactControls {
        impact: RigidImpactControls {
            maximum_bodies: 16,
            maximum_contacts: 21,
            geometry_tolerance_m: 1e-9,
            velocity_tolerance_mps: 1e-8,
            linear_momentum_tolerance_ns: 1e-9,
            angular_momentum_tolerance_nms: 1e-9,
            energy_tolerance_j: 1e-8,
            contact: FrictionalContactConfig::default(),
        },
        maximum_pairs: 120,
        maximum_events: 64,
        event_gap_tolerance_m: 1e-10,
        drift_energy_tolerance_j: 1e-8,
        total_energy_tolerance_j: 1e-7,
        linear_momentum_tolerance_ns: 1e-9,
        angular_momentum_tolerance_nms: 1e-8,
    }
}
fn laws(spheres: &[ContactSphere], mu: f64, e: f64) -> Vec<SpherePairLaw> {
    let mut ids: Vec<_> = spheres.iter().map(|s| s.id).collect();
    ids.sort();
    let mut out = Vec::new();
    for (i, &first) in ids.iter().enumerate() {
        for &second in &ids[i + 1..] {
            out.push(SpherePairLaw {
                first,
                second,
                friction_coefficient: mu,
                restitution: e,
            });
        }
    }
    out
}
fn spin_world(b: RigidBodyState) -> Vec3 {
    b.pose()
        .orientation()
        .rotate_body_to_world(b.angular_momentum_body())
}

#[test]
fn g1_hidden_crossing_unequal_mass_collision_has_analytic_time_and_endpoint() {
    let spheres = [
        ContactSphere {
            id: 10,
            radius_m: 0.5,
        },
        ContactSphere {
            id: 20,
            radius_m: 0.5,
        },
    ];
    let properties = [property(1.0, 0.5), property(3.0, 0.5)];
    let mut bodies = [
        body(Vec3::ZERO, Vec3::new(4.0, 0.0, 0.0), 1.0, Vec3::ZERO),
        body(Vec3::new(4.0, 0.0, 0.0), Vec3::ZERO, 3.0, Vec3::ZERO),
    ];
    // At t=2 the uncollided endpoints are separated again; endpoint-only
    // collision detection would incorrectly let the bodies pass through.
    let report = advance_sphere_contacts(
        &mut bodies,
        &properties,
        &spheres,
        &laws(&spheres, 0.0, 0.5),
        2.0,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(report.events.len(), 1);
    close(report.events[0].elapsed_s, 0.75);
    vector(
        bodies[0]
            .center_of_mass_velocity_world(properties[0])
            .unwrap(),
        Vec3::new(-0.5, 0.0, 0.0),
    );
    vector(
        bodies[1]
            .center_of_mass_velocity_world(properties[1])
            .unwrap(),
        Vec3::new(1.5, 0.0, 0.0),
    );
    vector(
        bodies[0].pose().position_world(),
        Vec3::new(2.375, 0.0, 0.0),
    );
    vector(
        bodies[1].pose().position_world(),
        Vec3::new(5.875, 0.0, 0.0),
    );
    close(report.unallocated_impact_loss_j, 4.5);
    close(report.kinetic_before_j, 8.0);
    close(report.kinetic_after_j, 3.5);
    close(report.total_energy_residual_j, 0.0);
    assert!(report.minimum_path_gap_m >= 0.0);
    assert!(report.events[0].maximum_gap_m <= controls().event_gap_tolerance_m);
}

#[test]
fn g1_spin_driven_stick_and_slide_use_actual_sphere_inertia() {
    let spheres = [
        ContactSphere {
            id: 1,
            radius_m: 0.5,
        },
        ContactSphere {
            id: 2,
            radius_m: 0.5,
        },
    ];
    let properties = [property(1.0, 0.5); 2];
    // The geometric normal is central, but the surface slip comes entirely
    // from sphere A's spin. W_t = 2/m + 2*r^2/I = 7; u_t = 1.
    let initial = [
        body(
            Vec3::new(-0.5, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            1.0,
            Vec3::new(0.0, 0.0, 0.2),
        ),
        body(Vec3::new(0.5, 0.0, 0.0), Vec3::ZERO, 1.0, Vec3::ZERO),
    ];
    for (mu, jt, loss) in [(0.5, -1.0 / 7.0, 1.0 / 14.0), (0.05, -0.05, 0.04125)] {
        let mut bodies = initial;
        let report = advance_sphere_contacts(
            &mut bodies,
            &properties,
            &spheres,
            &laws(&spheres, mu, 1.0),
            0.5,
            controls(),
            &CancelGate::new(),
        )
        .unwrap();
        assert_eq!(report.events.len(), 1);
        close(report.events[0].elapsed_s, 0.0);
        vector(bodies[0].linear_momentum_world(), Vec3::new(0.0, jt, 0.0));
        vector(bodies[1].linear_momentum_world(), Vec3::new(1.0, -jt, 0.0));
        vector(spin_world(bodies[0]), Vec3::new(0.0, 0.0, 0.2 + 0.5 * jt));
        vector(spin_world(bodies[1]), Vec3::new(0.0, 0.0, 0.5 * jt));
        close(report.unallocated_impact_loss_j, loss);
        close(report.total_energy_residual_j, 0.0);
        vector(report.angular_momentum_residual_nms, Vec3::ZERO);
        let slip = 1.0 + 7.0 * jt;
        if mu > 0.1 {
            close(slip, 0.0);
        } else {
            close(slip, 0.65);
        }
    }
}

#[test]
fn g1_simultaneous_symmetric_chain_is_not_sequential_pair_resolution() {
    let spheres = [
        ContactSphere {
            id: 30,
            radius_m: 0.5,
        },
        ContactSphere {
            id: 10,
            radius_m: 0.5,
        },
        ContactSphere {
            id: 20,
            radius_m: 0.5,
        },
    ];
    let properties = [property(1.0, 0.5); 3];
    let initial = [
        body(
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            1.0,
            Vec3::ZERO,
        ),
        body(Vec3::ZERO, Vec3::ZERO, 1.0, Vec3::ZERO),
        body(
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(-1.0, 0.0, 0.0),
            1.0,
            Vec3::ZERO,
        ),
    ];
    let mut bodies = initial;
    let report = advance_sphere_contacts(
        &mut bodies,
        &properties,
        &spheres,
        &laws(&spheres, 0.0, 1.0),
        1.0,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(report.events.len(), 1);
    assert_eq!(report.events[0].pairs, vec![[10, 20], [10, 30]]);
    for (b, v) in bodies.iter().zip([-1.0, 0.0, 1.0]) {
        close(b.linear_momentum_world().x, v);
    }
    close(report.unallocated_impact_loss_j, 0.0);
    // The complete canonical geometry/law order survives a body-array permutation.
    let perm = [2, 0, 1];
    let mut reordered = perm.map(|i| initial[i]);
    let spheres_reordered = perm.map(|i| spheres[i]);
    let mut laws_reordered = laws(&spheres, 0.0, 1.0);
    laws_reordered.reverse();
    let other = advance_sphere_contacts(
        &mut reordered,
        &properties,
        &spheres_reordered,
        &laws_reordered,
        1.0,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    for (slot, &i) in perm.iter().enumerate() {
        assert_eq!(reordered[slot], bodies[i]);
    }
    assert_eq!(other.events[0].pairs, report.events[0].pairs);
}

#[test]
fn g1_sequential_collisions_and_split_continuation_recover_newton_cradle() {
    let spheres = [
        ContactSphere {
            id: 1,
            radius_m: 0.5,
        },
        ContactSphere {
            id: 2,
            radius_m: 0.5,
        },
        ContactSphere {
            id: 3,
            radius_m: 0.5,
        },
    ];
    let properties = [property(1.0, 0.5); 3];
    let initial = [
        body(Vec3::ZERO, Vec3::new(2.0, 0.0, 0.0), 1.0, Vec3::ZERO),
        body(Vec3::new(3.0, 0.0, 0.0), Vec3::ZERO, 1.0, Vec3::ZERO),
        body(Vec3::new(6.0, 0.0, 0.0), Vec3::ZERO, 1.0, Vec3::ZERO),
    ];
    let law = laws(&spheres, 0.4, 1.0);
    let mut complete = initial;
    let report = advance_sphere_contacts(
        &mut complete,
        &properties,
        &spheres,
        &law,
        3.0,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(report.events.len(), 2);
    close(report.events[0].elapsed_s, 1.0);
    close(report.events[1].elapsed_s, 2.0);
    for (b, x) in complete.iter().zip([2.0, 5.0, 8.0]) {
        close(b.pose().position_world().x, x);
    }
    for (b, v) in complete.iter().zip([0.0, 0.0, 2.0]) {
        close(b.linear_momentum_world().x, v);
    }
    let mut split = initial;
    advance_sphere_contacts(
        &mut split,
        &properties,
        &spheres,
        &law,
        1.5,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    advance_sphere_contacts(
        &mut split,
        &properties,
        &spheres,
        &law,
        1.5,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(split, complete);
}

#[test]
fn g1_tangent_touch_needs_no_impulse_and_free_asymmetric_drift_is_original() {
    let spheres = [
        ContactSphere {
            id: 1,
            radius_m: 0.5,
        },
        ContactSphere {
            id: 2,
            radius_m: 0.5,
        },
    ];
    let properties = [property(1.0, 0.5); 2];
    let mut tangent = [
        body(
            Vec3::new(-2.0, 1.0, 0.0),
            Vec3::new(4.0, 0.0, 0.0),
            1.0,
            Vec3::ZERO,
        ),
        body(Vec3::ZERO, Vec3::ZERO, 1.0, Vec3::ZERO),
    ];
    let report = advance_sphere_contacts(
        &mut tangent,
        &properties,
        &spheres,
        &laws(&spheres, 0.8, 0.0),
        1.0,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert!(report.events.is_empty());
    close(report.minimum_path_gap_m, 0.0);
    vector(tangent[0].pose().position_world(), Vec3::new(2.0, 1.0, 0.0));
    close(report.unallocated_impact_loss_j, 0.0);
    let props = [MassProperties::new(2.0, Vec3::ZERO, Vec3::new(0.3, 0.4, 0.5)).unwrap()];
    let sphere = [ContactSphere {
        id: 9,
        radius_m: 1.0,
    }];
    let initial = body(
        Vec3::new(1.0, 2.0, 3.0),
        Vec3::new(0.2, -0.1, 0.4),
        2.0,
        Vec3::new(0.13, 0.2, -0.07),
    );
    let expected = RigidBodyIntegrator::new(Gravity::ZERO)
        .step(initial, props[0], Wrench::ZERO, 0.2)
        .unwrap()
        .state_after;
    let mut actual = [initial];
    let report = advance_sphere_contacts(
        &mut actual,
        &props,
        &sphere,
        &[],
        0.2,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(actual[0], expected);
    assert!(report.events.is_empty());
}

#[test]
fn g4_initial_geometry_late_event_budget_and_cancellation_refuse_atomically() {
    let spheres = [
        ContactSphere {
            id: 1,
            radius_m: 0.5,
        },
        ContactSphere {
            id: 2,
            radius_m: 0.5,
        },
        ContactSphere {
            id: 3,
            radius_m: 0.5,
        },
    ];
    let properties = [property(1.0, 0.5); 3];
    let initial = [
        body(Vec3::ZERO, Vec3::new(2.0, 0.0, 0.0), 1.0, Vec3::ZERO),
        body(Vec3::new(3.0, 0.0, 0.0), Vec3::ZERO, 1.0, Vec3::ZERO),
        body(Vec3::new(6.0, 0.0, 0.0), Vec3::ZERO, 1.0, Vec3::ZERO),
    ];
    let law = laws(&spheres, 0.0, 1.0);
    let mut bodies = initial;
    let mut cap = controls();
    cap.maximum_events = 1;
    assert!(matches!(
        advance_sphere_contacts(
            &mut bodies,
            &properties,
            &spheres,
            &law,
            3.0,
            cap,
            &CancelGate::new()
        ),
        Err(SphereContactError::EventBudget)
    ));
    assert_eq!(bodies, initial);
    let gate = CancelGate::new();
    gate.request();
    assert!(matches!(
        advance_sphere_contacts(
            &mut bodies,
            &properties,
            &spheres,
            &law,
            3.0,
            controls(),
            &gate
        ),
        Err(SphereContactError::Cancelled)
    ));
    assert_eq!(bodies, initial);
    assert!(matches!(
        advance_sphere_contacts(
            &mut bodies,
            &properties,
            &spheres,
            &law[..2],
            3.0,
            controls(),
            &CancelGate::new()
        ),
        Err(SphereContactError::Invalid(_))
    ));
    assert_eq!(bodies, initial);
    let penetrating = [
        initial[0],
        body(
            Vec3::new(1.0 - 1e-12, 0.0, 0.0),
            Vec3::ZERO,
            1.0,
            Vec3::ZERO,
        ),
        initial[2],
    ];
    let mut bodies = penetrating;
    assert!(matches!(
        advance_sphere_contacts(
            &mut bodies,
            &properties,
            &spheres,
            &law,
            3.0,
            controls(),
            &CancelGate::new()
        ),
        Err(SphereContactError::Penetration { .. })
    ));
    assert_eq!(bodies, penetrating);
    let mut retry = initial;
    advance_sphere_contacts(
        &mut retry,
        &properties,
        &spheres,
        &law,
        3.0,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    close(retry[2].pose().position_world().x, 8.0);
}
