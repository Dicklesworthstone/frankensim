//! Sub-ULP motion remains part of the original accepted trajectory.
use super::*;
use fs_couple::rod_motion::impulse::{RodImpulseControls, RodPointImpulse};

fn free_uniform(velocity: f64, spin: f64, y: f64) -> RodMotion {
    let gate = CancelGate::new();
    let model = Arc::new(
        SpatialMaterialRod::try_straight(0.8, vec![section(1.0)])
            .unwrap()
            .into_motion_model(RodRootBoundary::Free, &gate)
            .unwrap(),
    );
    let mut pose = model.initial_pose();
    for point in &mut pose.positions_m { point[1] = y; }
    RodMotion::new(
        model, pose, &[[0.0, velocity, 0.0]; 2], &[[spin, 0.0, 0.0]; 2], 0.0, &gate,
    ).unwrap()
}

#[test]
fn uniform_force_retains_the_full_impulse_across_both_native_kicks() {
    let gate = CancelGate::new();
    for sign in [-1.0, 1.0] {
        let original = free_uniform(sign * 0.25, 0.0, 0.0);
        let p0 = original.linear_momenta_kg_m_s()[0][1];
        let spacing = p0.abs().next_up() - p0.abs();
        let force = sign * 0.75 * spacing;
        let loads = [RodNodalWrench { force_n: [0.0, force, 0.0], torque_nm: [0.0; 3] }; 2];
        // Each separately rounded half-kick loses 3/8 ULP; the complete
        // physical impulse rounds to the adjacent representable momentum.
        assert_eq!(0.5_f64.mul_add(force, p0), p0);
        let (next, report) = original.step_candidate(1.0, &loads, controls(1.0), &gate).unwrap();
        for p in next.linear_momenta_kg_m_s() {
            assert_eq!(p[1], p0 + force);
            assert_eq!([p[0], p[2]], [0.0; 2]);
        }
        let expected_work = 2.0 * force * sign * 0.25;
        near(report.source_work_quadrature_j, expected_work,
            64.0 * f64::EPSILON * expected_work.abs());
        assert_eq!(next.body_spins_kg_m2_s(), original.body_spins_kg_m2_s());
    }
}

#[test]
fn retained_half_impulse_reaches_the_native_translation_after_cancellation() {
    let gate = CancelGate::new();
    let speed = 9_007_199_254_740_992.0;
    for sign in [-1.0, 1.0] {
        // Uniform translation leaves the actual source rod undeformed.
        // The large initial position cancels the large free drift exactly,
        // exposing the small displacement caused by the midpoint force.
        let original = free_uniform(sign * speed, 0.0, -sign * speed);
        let p0 = original.linear_momenta_kg_m_s()[0][1];
        let mass = original.model().inertia()[0].mass_kg();
        assert_eq!(p0 / mass, sign * speed);
        let force = sign * 0.75 * (p0.abs().next_up() - p0.abs());
        let loads = [RodNodalWrench { force_n: [0.0, force, 0.0], torque_nm: [0.0; 3] }; 2];
        let mut settings = controls(1.0);
        // The existing energy diagnostic subtracts two large f64 kinetic
        // totals. Bound that subtraction's scale explicitly in this fixture.
        settings.energy_defect_limit_j = 128.0 * f64::EPSILON
            * original.observe([0.0; 3]).unwrap().kinetic_energy_j;
        let (next, _) = original.step_candidate(1.0, &loads, settings, &gate).unwrap();
        for node in 0..2 {
            assert_eq!(next.linear_momenta_kg_m_s()[node][1], p0 + force);
            near(next.pose().positions_m[node][1], 0.5 * force / mass, 2e-14);
            assert_eq!(next.pose().positions_m[node][0], original.pose().positions_m[node][0]);
        }
    }
}

#[test]
fn longitudinal_torque_retains_full_spin_impulse_through_original_free_rotation() {
    let gate = CancelGate::new();
    let dt = 0.125;
    for sign in [-1.0, 1.0] {
        let original = free_uniform(0.0, sign, 0.0);
        let spin = original.body_spins_kg_m2_s()[0][0];
        let impulse = sign * 0.75 * (spin.abs().next_up() - spin.abs());
        let torque = impulse / dt;
        assert_eq!((0.5 * dt).mul_add(torque, spin), spin);
        let loads = [RodNodalWrench { force_n: [0.0; 3], torque_nm: [torque, 0.0, 0.0] }; 2];
        let mut settings = controls(dt);
        settings.maximum_rotation_rad = 0.5;
        let (next, report) = original.step_candidate(dt, &loads, settings, &gate).unwrap();
        let expected = So3::exp(So3Tangent::new(Vec3::new(sign * dt, 0.0, 0.0))).unwrap();
        for node in 0..2 {
            assert_eq!(next.body_spins_kg_m2_s()[node][0], spin + impulse);
            assert_eq!(next.pose().positions_m[node], original.pose().positions_m[node]);
            assert!(next.pose().frames[node].body_minus(expected).unwrap().angular.norm() < 2e-14);
        }
        let expected_work = 2.0 * impulse * sign;
        near(report.source_work_quadrature_j, expected_work,
            64.0 * f64::EPSILON * expected_work.abs());
        assert_eq!(next.linear_momenta_kg_m_s(), original.linear_momenta_kg_m_s());
    }
}

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
