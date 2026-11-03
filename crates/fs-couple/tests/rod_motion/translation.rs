//! Sub-ULP motion remains part of the original accepted trajectory.
use super::*;
use fs_couple::rod_motion::impulse::{RodImpulseControls, RodPointImpulse};

#[test]
fn tiny_uniform_translation_accumulates_through_clone_impulse_and_refusal() {
    let gate = CancelGate::new();
    let model = Arc::new(
        SpatialMaterialRod::try_straight(0.8, vec![section(1.0)])
            .unwrap()
            .into_motion_model(RodRootBoundary::Free, &gate)
            .unwrap(),
    );
    let mut pose = model.initial_pose();
    for point in &mut pose.positions_m {
        point[1] = 1.0;
    }
    // An exact uniform rigid translation has zero strain and no elastic force.
    // Each one-second displacement is only 1/8 ULP at the supplied y origin.
    let speed = f64::EPSILON / 8.0;
    let original = RodMotion::new(
        model,
        pose,
        &[[0.0, speed, 0.0]; 2],
        &[[0.0; 3]; 2],
        0.0,
        &gate,
    )
    .unwrap();
    let checkpoint = evolve(original.clone(), 1.0, 1);
    assert_eq!(checkpoint.pose(), original.pose());
    assert_eq!(checkpoint.accepted_steps(), 1);

    let direct = evolve(original.clone(), 1.0, 32);
    let resumed = evolve(checkpoint.clone(), 1.0, 31);
    assert_eq!(direct.pose(), resumed.pose());
    for point in &direct.pose().positions_m {
        assert_eq!(point[1], 1.0 + 32.0 * speed);
    }
    assert_eq!(
        direct.linear_momenta_kg_m_s(),
        original.linear_momenta_kg_m_s()
    );
    assert_eq!(direct.body_spins_kg_m2_s(), original.body_spins_kg_m2_s());

    // A real instantaneous impulse doubles the common velocity. Two more
    // intervals add 1/2 ULP. The carried 1/8 ULP must survive the impulse or the
    // resulting y rounds incorrectly back to 1.0 at this exact halfway case.
    let impulses: Vec<_> = checkpoint
        .model()
        .inertia()
        .iter()
        .enumerate()
        .map(|(node, inertia)| RodPointImpulse {
            node,
            point_body_m: [0.0; 3],
            linear_world_n_s: [0.0, inertia.mass_kg() * speed, 0.0],
            couple_world_nm_s: [0.0; 3],
        })
        .collect();
    let (accelerated, _) = checkpoint
        .apply_point_impulses(
            &impulses,
            RodImpulseControls {
                balance_origin_m: [0.0; 3],
                linear_tolerance_n_s: 1e-30,
                angular_tolerance_nm_s: 1e-30,
                work_tolerance_j: 1e-45,
            },
            &gate,
        )
        .unwrap();
    let accelerated = evolve(accelerated, 1.0, 2);
    for point in &accelerated.pose().positions_m {
        assert_eq!(point[1], 1.0 + f64::EPSILON);
    }

    // Refuse after advancing a private candidate through its drift and forces;
    // the accepted checkpoint must retain its small unrepresented translation.
    let mut strict = controls(1e-3);
    strict.energy_defect_limit_j = 0.0;
    let loads = [
        RodNodalWrench::default(),
        RodNodalWrench {
            force_n: [0.0, 1.0, 0.0],
            torque_nm: [0.0; 3],
        },
    ];
    assert!(matches!(
        checkpoint.step_candidate(1e-3, &loads, strict, &gate),
        Err(RodMotionError::Energy { .. })
    ));
    let retry = evolve(checkpoint.clone(), 1.0, 31);
    assert_eq!(retry.pose(), direct.pose());
    assert_eq!(
        retry.linear_momenta_kg_m_s(),
        direct.linear_momenta_kg_m_s()
    );
    assert_eq!(checkpoint.pose(), original.pose());
    assert_eq!(checkpoint.accepted_steps(), 1);
}
