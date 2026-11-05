//! Automatic geometry and sustained load transfer across one large stack.
use super::*;

#[test]
fn thirty_two_connected_spheres_stop_and_transfer_the_entire_stack_weight() {
    const COUNT: usize = 32;
    let duration = 1.0 / 32.0;
    let dt = duration / 2.0;
    let properties = vec![property(1.0); COUNT];
    let spheres: Vec<_> = (0..COUNT).map(|i| sphere(10 + i as u64)).collect();
    // The original explicit activation clearance is retained. It does not
    // loosen momentum, energy or complete-path acceptance for this larger set.
    let initial: Vec<_> = (0..COUNT)
        .map(|i| {
            body(
                Vec3::new(
                    0.0,
                    0.0,
                    RADIUS + CLEARANCE + i as f64 * (2.0 * RADIUS + CLEARANCE),
                ),
                Vec3::new(0.0, 0.0, -0.5),
                1.0,
            )
        })
        .collect();
    let mut bodies = initial.clone();
    let mut scene_controls = controls(dt);
    scene_controls.step.impact.maximum_bodies = COUNT;
    scene_controls.step.impact.maximum_contacts = COUNT;
    scene_controls.step.impact.contact.maximum_contacts = COUNT;
    scene_controls.step.impact.contact.matrix.max_variables = 3 * COUNT;
    // This opt-in larger component also needs a larger linear work budget.
    scene_controls.step.impact.contact.newton.linear_restart = 64;
    scene_controls.step.impact.contact.newton.max_linear_cycles = 16;
    let report = advance_sphere_scene(
        &mut bodies,
        &properties,
        &spheres,
        &[plane(Vec3::ZERO)],
        &vec![Wrench::ZERO; COUNT],
        Gravity::new(Vec3::new(0.0, 0.0, -8.0)).unwrap(),
        0.0,
        duration,
        scene_controls,
        |_| Some(law(0.0, 0.0)),
        &CancelGate::new(),
    )
    .unwrap();
    balances(&report, duration);
    assert_eq!(report.steps.len(), 2);
    assert_eq!(report.rejected_steps, 0);
    for (step_index, step) in report.steps.iter().enumerate() {
        assert_eq!(step.contact_keys.len(), COUNT);
        assert_eq!(step.native.impact.contact.islands.len(), 1);
        let island = &step.native.impact.contact.islands[0];
        assert_eq!(island.body_indices, (0..COUNT).collect::<Vec<_>>());
        assert_eq!(island.contact_indices.len(), COUNT);
        // The first simultaneous solve stops initial approach plus the force
        // kick. The next solve transfers each upper body's actual weight.
        let stopping_speed = if step_index == 0 { 0.5 } else { 0.0 } + 8.0 * dt;
        for (key, impulse) in step
            .contact_keys
            .iter()
            .zip(&step.native.impact.contact.impulses)
        {
            let supported_bodies = match *key {
                SphereSceneContactKey::SpherePlane {
                    sphere: 10,
                    plane: 99,
                } => COUNT,
                SphereSceneContactKey::Spheres { first, second } => {
                    assert_eq!(second, first + 1);
                    COUNT - (second - 10) as usize
                }
                _ => panic!("unexpected stack contact {key:?}"),
            };
            close(impulse[0], supported_bodies as f64 * stopping_speed, 2e-8);
            close(impulse[1], 0.0, 2e-8);
            close(impulse[2], 0.0, 2e-8);
        }
    }
    for (after, before) in bodies.iter().zip(&initial) {
        vector(after.linear_momentum_world(), Vec3::ZERO, 2e-8);
        vector(spin_world(*after), Vec3::ZERO, 2e-9);
        vector(
            after.pose().position_world(),
            before.pose().position_world(),
            2e-9,
        );
    }
    // The declared kick/contact split contributes 2.5 J of force work,
    // balanced with the initial 4 J by 6.5 J of unallocated contact loss.
    // Stationary supports supply impulse but no energy or displacement work.
    vector(report.support_impulse_ns, Vec3::new(0.0, 0.0, 24.0), 3e-8);
    close(report.force_impulse_work_j, 2.5, 2e-8);
    close(report.force_displacement_work_j, 0.0, 2e-8);
    close(report.support_work_j, 0.0, 2e-8);
    close(report.unallocated_contact_loss_j, 6.5, 2e-8);
    close(report.kinetic_after_j, 0.0, 2e-8);
}
