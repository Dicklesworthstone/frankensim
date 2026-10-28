//! G1/G3/G4: real accelerated native bodies, internal friction and gravity work.
use super::*;

#[test]
fn g1_falling_unequal_spheres_collide_at_the_relative_analytic_time() {
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
        body(
            Vec3::new(0.0, 0.0, 10.0),
            Vec3::new(4.0, 0.0, 0.0),
            1.0,
            Vec3::ZERO,
        ),
        body(Vec3::new(4.0, 0.0, 10.0), Vec3::ZERO, 3.0, Vec3::ZERO),
    ];
    let report = advance_sphere_contacts(
        &mut bodies,
        &properties,
        &spheres,
        &laws(&spheres, 0.0, 0.5),
        Gravity::new(Vec3::new(0.0, 0.0, -10.0)).unwrap(),
        2.0,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(report.events.len(), 1);
    close(report.events[0].elapsed_s, 0.75);
    vector(
        bodies[0].pose().position_world(),
        Vec3::new(2.375, 0.0, -10.0),
    );
    vector(
        bodies[1].pose().position_world(),
        Vec3::new(5.875, 0.0, -10.0),
    );
    vector(
        bodies[0].linear_momentum_world(),
        Vec3::new(-0.5, 0.0, -20.0),
    );
    vector(
        bodies[1].linear_momentum_world(),
        Vec3::new(4.5, 0.0, -60.0),
    );
    // m_total=4, each centre falls20m. Initially sum(m*x)=12 along x
    // and total p_x=4, so the exact integrated gravity torque is320 about y.
    close(report.gravity_work_j, 800.0);
    vector(report.gravity_impulse_ns, Vec3::new(0.0, 0.0, -80.0));
    vector(
        report.gravity_angular_impulse_nms,
        Vec3::new(0.0, 320.0, 0.0),
    );
    close(report.unallocated_impact_loss_j, 4.5);
    close(report.kinetic_after_j, 803.5);
    close(report.total_energy_residual_j, 0.0);
    vector(report.linear_momentum_residual_ns, Vec3::ZERO);
    vector(report.angular_momentum_residual_nms, Vec3::ZERO);
    assert!(report.minimum_path_gap_m >= 0.0);
}

#[test]
fn g3_free_fall_equivalence_preserves_spin_friction_and_multiple_impacts() {
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
        body(
            Vec3::ZERO,
            Vec3::new(2.0, 0.0, 0.0),
            1.0,
            Vec3::new(0.0, 0.0, 0.1),
        ),
        body(Vec3::new(1.0, 0.0, 0.0), Vec3::ZERO, 1.0, Vec3::ZERO),
        body(Vec3::new(4.0, 0.0, 0.0), Vec3::ZERO, 1.0, Vec3::ZERO),
    ];
    let law = laws(&spheres, 0.4, 1.0);
    let mut free = initial;
    let free_report = advance_sphere_contacts(
        &mut free,
        &properties,
        &spheres,
        &law,
        Gravity::ZERO,
        2.0,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert!(free_report.events.len() >= 2);
    assert!(free_report.unallocated_impact_loss_j > 0.0);
    let g = Vec3::new(0.2, -0.3, -1.7);
    let mut falling = initial;
    let report = advance_sphere_contacts(
        &mut falling,
        &properties,
        &spheres,
        &law,
        Gravity::new(g).unwrap(),
        2.0,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(report.events.len(), free_report.events.len());
    for (a, b) in report.events.iter().zip(&free_report.events) {
        close(a.elapsed_s, b.elapsed_s);
        assert_eq!(a.pairs, b.pairs);
    }
    for i in 0..3 {
        vector(
            falling[i]
                .pose()
                .position_world()
                .sub(free[i].pose().position_world()),
            g.scale(2.0),
        );
        vector(
            falling[i]
                .linear_momentum_world()
                .sub(free[i].linear_momentum_world()),
            g.scale(2.0),
        );
        vector(spin_world(falling[i]), spin_world(free[i]));
        for (a, b) in falling[i]
            .pose()
            .orientation()
            .components()
            .into_iter()
            .zip(free[i].pose().orientation().components())
        {
            close(a, b);
        }
    }
    close(
        report.unallocated_impact_loss_j,
        free_report.unallocated_impact_loss_j,
    );
    // Work = g dot P_initial * T + M|g|²T²/2; all impact exchanges are internal.
    close(report.gravity_work_j, 18.92);
    vector(
        report.gravity_angular_impulse_nms,
        Vec3::new(0.0, 23.8, -4.2),
    );
    close(report.total_energy_residual_j, 0.0);
}

#[test]
fn g1_asymmetric_gravity_drift_is_the_original_native_update() {
    let properties = [MassProperties::new(2.0, Vec3::ZERO, Vec3::new(0.3, 0.4, 0.5)).unwrap()];
    let spheres = [ContactSphere {
        id: 9,
        radius_m: 1.0,
    }];
    let initial = body(
        Vec3::new(1.0, 2.0, 3.0),
        Vec3::new(0.2, -0.1, 0.4),
        2.0,
        Vec3::new(0.13, 0.2, -0.07),
    );
    let gravity = Gravity::new(Vec3::new(0.4, -0.2, -9.8)).unwrap();
    let expected = RigidBodyIntegrator::new(gravity)
        .step(initial, properties[0], Wrench::ZERO, 0.2)
        .unwrap()
        .state_after;
    let mut actual = [initial];
    let report = advance_sphere_contacts(
        &mut actual,
        &properties,
        &spheres,
        &[],
        gravity,
        0.2,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    assert_eq!(actual[0], expected);
    assert!(report.events.is_empty());
    let displacement = Vec3::new(0.048, -0.024, -0.116);
    vector(
        actual[0]
            .pose()
            .position_world()
            .sub(initial.pose().position_world()),
        displacement,
    );
    close(
        report.gravity_work_j,
        2.0 * gravity.acceleration_world().dot(displacement),
    );
    close(report.drift_energy_defect_j, 0.0);
    vector(report.angular_momentum_residual_nms, Vec3::ZERO);
}

#[test]
fn g4_late_gravity_event_budget_is_atomic_and_split_continuation_agrees() {
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
    let gravity = Gravity::new(Vec3::new(0.0, 0.0, -9.8)).unwrap();
    let law = laws(&spheres, 0.4, 1.0);
    let mut bodies = initial;
    let mut capped = controls();
    capped.maximum_events = 1;
    assert!(matches!(
        advance_sphere_contacts(
            &mut bodies,
            &properties,
            &spheres,
            &law,
            gravity,
            3.0,
            capped,
            &CancelGate::new()
        ),
        Err(SphereContactError::EventBudget)
    ));
    assert_eq!(bodies, initial);
    let report = advance_sphere_contacts(
        &mut bodies,
        &properties,
        &spheres,
        &law,
        gravity,
        3.0,
        controls(),
        &CancelGate::new(),
    )
    .unwrap();
    let mut split = initial;
    let mut work = 0.0;
    for _ in 0..2 {
        work += advance_sphere_contacts(
            &mut split,
            &properties,
            &spheres,
            &law,
            gravity,
            1.5,
            controls(),
            &CancelGate::new(),
        )
        .unwrap()
        .gravity_work_j;
    }
    for i in 0..3 {
        vector(
            split[i].pose().position_world(),
            bodies[i].pose().position_world(),
        );
        vector(
            split[i].linear_momentum_world(),
            bodies[i].linear_momentum_world(),
        );
    }
    close(work, report.gravity_work_j);
    close(report.unallocated_impact_loss_j, 0.0);
    close(bodies[2].pose().position_world().x, 8.0);
    close(bodies[2].pose().position_world().z, -44.1);
}
