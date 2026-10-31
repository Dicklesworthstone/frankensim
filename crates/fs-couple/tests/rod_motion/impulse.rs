//! G0/G1/G3/G4: physical attachment strikes on the original rod trajectory.
use super::*;
use fs_couple::rod_motion::SpatialRodMotion;
use fs_couple::rod_motion::impulse::{RodImpulseControls, RodPointImpulse};
use fs_ga::Mat3;

const OFFSET: [f64; 3] = [0.05, 0.08, -0.03];

fn impulse_controls(origin: [f64; 3]) -> RodImpulseControls {
    RodImpulseControls {
        balance_origin_m: origin,
        linear_tolerance_n_s: 1e-12,
        angular_tolerance_nm_s: 1e-12,
        work_tolerance_j: 1e-12,
    }
}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|i| a[i] + b[i])
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn rotate(r: So3, x: [f64; 3]) -> [f64; 3] {
    let y = r.rotate(Vec3::new(x[0], x[1], x[2])).unwrap();
    [y.x, y.y, y.z]
}
fn vector_near(a: [f64; 3], b: [f64; 3], tolerance: f64) {
    for c in 0..3 {
        near(a[c], b[c], tolerance);
    }
}
fn eccentric_source() -> SpatialMaterialRod {
    let material = section(1.0);
    let link = MaterialTensionLink::try_new(
        material.elastic().clone(),
        [
            TensionLinkAttachment {
                node: 0,
                offset_body_m: OFFSET,
            },
            TensionLinkAttachment {
                node: 1,
                offset_body_m: OFFSET,
            },
        ],
        0.9,
        1,
        0.2,
    )
    .unwrap();
    SpatialMaterialRod::try_straight(0.8, vec![material])
        .unwrap()
        .with_tension_links(vec![link], &CancelGate::new())
        .unwrap()
}

#[test]
fn eccentric_point_strike_matches_full_spatial_mass_and_spin_orbital_impulse() {
    let gate = CancelGate::new();
    let model = Arc::new(
        eccentric_source()
            .into_spatial_motion_model(RodRootBoundary::Free, &gate)
            .unwrap(),
    );
    let rotation = So3::exp(So3Tangent::new(Vec3::new(0.4, -0.2, 0.3))).unwrap();
    let mut pose = model.initial_pose();
    for node in 0..2 {
        pose.positions_m[node] = add(rotate(rotation, pose.positions_m[node]), [0.7, -0.2, 0.1]);
        pose.frames[node] = rotation;
    }
    let velocity = [[0.1, 0.2, -0.1]; 2];
    let omega = [[0.2, -0.1, 0.3]; 2];
    let state =
        SpatialRodMotion::new_spatial(model, pose.clone(), &velocity, &omega, 7.5, 1e-12, &gate)
            .unwrap();
    let force = RodPointImpulse {
        node: 1,
        point_body_m: [0.03, -0.04, 0.05],
        linear_world_n_s: [0.03, -0.01, 0.02],
        couple_world_nm_s: [0.0; 3],
    };
    let couple = RodPointImpulse {
        node: 1,
        point_body_m: [0.0; 3],
        linear_world_n_s: [0.0; 3],
        couple_world_nm_s: [0.002, 0.001, -0.003],
    };
    let origin = [-0.3, 0.5, 0.2];
    let (next, report) = state
        .apply_point_impulses(&[force, couple], impulse_controls(origin), &gate)
        .unwrap();

    // Independent six-by-six material-node mass solve retains translation/
    // rotation coupling; production uses its prelowered central tensor.
    let inertia = state.model().inertia()[1];
    let center = inertia.center_of_mass_body_m();
    let tensor = inertia.inertia_com_body_kg_m2();
    let spatial = fs_mbd::articulated::SpatialInertia::new(
        inertia.mass_kg(),
        Vec3::new(center[0], center[1], center[2]),
        Mat3 {
            m: std::array::from_fn(|i| tensor[i / 3][i % 3]),
        },
    )
    .unwrap()
    .matrix()
    .unwrap();
    let body_j = rotate(rotation.inverse(), force.linear_world_n_s);
    let body_k = add(
        cross(force.point_body_m, body_j),
        rotate(rotation.inverse(), couple.couple_world_nm_s),
    );
    let mut change = [
        body_k[0], body_k[1], body_k[2], body_j[0], body_j[1], body_j[2],
    ];
    fs_la::factor::lu(&spatial.m, 6).unwrap().solve(&mut change);
    let expected_omega = std::array::from_fn(|c| omega[1][c] + change[c]);
    let expected_velocity = add(
        velocity[1],
        rotate(rotation, [change[3], change[4], change[5]]),
    );
    vector_near(
        next.angular_velocities_body_rad_s().unwrap()[1],
        expected_omega,
        2e-11,
    );
    vector_near(next.velocities_m_s().unwrap()[1], expected_velocity, 2e-12);
    vector_near(
        next.linear_momenta_kg_m_s()[0],
        state.linear_momenta_kg_m_s()[0],
        0.0,
    );
    vector_near(
        next.body_spins_kg_m2_s()[0],
        state.body_spins_kg_m2_s()[0],
        0.0,
    );

    let point = add(pose.positions_m[1], rotate(rotation, force.point_body_m));
    let arm = std::array::from_fn(|c| point[c] - origin[c]);
    let moment = add(cross(arm, force.linear_world_n_s), couple.couple_world_nm_s);
    vector_near(
        report.applied_linear_impulse_n_s,
        force.linear_world_n_s,
        2e-13,
    );
    vector_near(report.applied_angular_impulse_nm_s, moment, 2e-13);
    vector_near(report.linear_balance_residual_n_s, [0.0; 3], 2e-13);
    vector_near(report.angular_balance_residual_nm_s, [0.0; 3], 2e-13);
    let before_point_velocity = add(
        velocity[1],
        rotate(rotation, cross(omega[1], force.point_body_m)),
    );
    let after_point_velocity = add(
        expected_velocity,
        rotate(rotation, cross(expected_omega, force.point_body_m)),
    );
    let world_omega_before = rotate(rotation, omega[1]);
    let world_omega_after = rotate(rotation, expected_omega);
    let exact_work: f64 = (0..3)
        .map(|c| {
            0.5 * force.linear_world_n_s[c] * (before_point_velocity[c] + after_point_velocity[c])
                + 0.5 * couple.couple_world_nm_s[c] * (world_omega_before[c] + world_omega_after[c])
        })
        .sum();
    near(report.applied_work_j, exact_work, 2e-13);
    near(report.kinetic_energy_change_j, exact_work, 2e-13);
    near(report.energy_residual_j, 0.0, 2e-13);
    assert_eq!(next.pose(), &pose);
    assert_eq!(next.time_s(), state.time_s());
    assert_eq!(next.accepted_steps(), state.accepted_steps());
    assert_eq!(
        report.after.potential_energy_j,
        report.before.potential_energy_j
    );
}

#[test]
fn eccentric_supports_take_the_correct_node_couple_and_do_no_work() {
    let gate = CancelGate::new();
    let source = eccentric_source()
        .with_supports(
            &[RodSupport {
                node: 1,
                position_m: [None, Some(0.0), None],
                frame: Some(So3::identity()),
            }],
            &gate,
        )
        .unwrap();
    let model = Arc::new(
        source
            .into_spatial_motion_model(RodRootBoundary::Clamped, &gate)
            .unwrap(),
    );
    let pose = model.initial_pose();
    let state = SpatialRodMotion::new_spatial(
        model,
        pose.clone(),
        &[[0.0; 3]; 2],
        &[[0.0; 3]; 2],
        0.0,
        1e-12,
        &gate,
    )
    .unwrap();
    let impulse = RodPointImpulse {
        node: 0,
        point_body_m: [0.03, -0.04, 0.05],
        linear_world_n_s: [0.003, -0.002, 0.004],
        couple_world_nm_s: [0.002, 0.001, -0.003],
    };
    let tip = RodPointImpulse { node: 1, ..impulse };
    let (next, report) = state
        .apply_point_impulses(&[impulse, tip], impulse_controls([0.0; 3]), &gate)
        .unwrap();
    let node_moment = add(
        cross(impulse.point_body_m, impulse.linear_world_n_s),
        impulse.couple_world_nm_s,
    );
    vector_near(
        report.support_impulses[0].linear_n_s,
        impulse.linear_world_n_s.map(|j| -j),
        1e-15,
    );
    vector_near(
        report.support_impulses[0].angular_nm_s,
        node_moment.map(|k| -k),
        1e-15,
    );
    vector_near(next.linear_momenta_kg_m_s()[0], [0.0; 3], 0.0);
    let transmitted = [tip.linear_world_n_s[0], 0.0, tip.linear_world_n_s[2]];
    let center = state.model().inertia()[1].center_of_mass_body_m();
    let expected_couple = add(node_moment.map(|k| -k), cross(center, transmitted));
    vector_near(
        report.support_impulses[1].linear_n_s,
        [0.0, -tip.linear_world_n_s[1], 0.0],
        1e-15,
    );
    vector_near(
        report.support_impulses[1].angular_nm_s,
        expected_couple,
        1e-15,
    );
    vector_near(next.linear_momenta_kg_m_s()[1], transmitted, 1e-15);
    assert_eq!(next.body_spins_kg_m2_s(), &[[0.0; 3]; 2]);
    let energy = transmitted
        .into_iter()
        .map(|j| 0.5 * j * j / state.model().inertia()[1].mass_kg())
        .sum::<f64>();
    near(report.applied_work_j, energy, 1e-15);
    near(report.after.kinetic_energy_j, energy, 1e-15);
    near(report.energy_residual_j, 0.0, 1e-15);
    vector_near(report.linear_balance_residual_n_s, [0.0; 3], 1e-15);
    vector_near(report.angular_balance_residual_nm_s, [0.0; 3], 1e-15);
    assert_eq!(next.pose(), &pose);
}

#[test]
fn reciprocal_point_strike_drives_bending_and_torsion_without_net_momentum() {
    let gate = CancelGate::new();
    let model = Arc::new(
        SpatialMaterialRod::try_straight(0.8, vec![section(1.0)])
            .unwrap()
            .into_motion_model(RodRootBoundary::Free, &gate)
            .unwrap(),
    );
    let state = RodMotion::new(
        model.clone(),
        model.initial_pose(),
        &[[0.0; 3]; 2],
        &[[0.0; 3]; 2],
        0.0,
        &gate,
    )
    .unwrap();
    // Equal/opposite impulses act at the SAME world point via two rigid arms.
    let impulses = [
        RodPointImpulse {
            node: 0,
            point_body_m: [0.4, 0.05, 0.0],
            linear_world_n_s: [0.0, 0.0, 0.001],
            couple_world_nm_s: [0.0; 3],
        },
        RodPointImpulse {
            node: 1,
            point_body_m: [-0.4, 0.05, 0.0],
            linear_world_n_s: [0.0, 0.0, -0.001],
            couple_world_nm_s: [0.0; 3],
        },
    ];
    let (struck, report) = state
        .apply_point_impulses(&impulses, impulse_controls([0.0; 3]), &gate)
        .unwrap();
    vector_near(report.after.linear_momentum_kg_m_s, [0.0; 3], 1e-15);
    vector_near(report.after.angular_momentum_kg_m2_s, [0.0; 3], 1e-15);
    assert!(struck.body_spins_kg_m2_s()[0][0] > 0.0);
    assert!(struck.body_spins_kg_m2_s()[0][1] < 0.0);
    assert!(report.applied_work_j > 0.0);
    let next = evolve(struck, 1e-4, 50);
    let after = next.observe([0.0; 3]).unwrap();
    assert!(after.potential_energy_j > report.after.potential_energy_j + 1e-14);
    vector_near(after.linear_momentum_kg_m_s, [0.0; 3], 2e-11);
    vector_near(after.angular_momentum_kg_m2_s, [0.0; 3], 2e-11);
    assert_eq!(next.accepted_steps(), 50);
}

#[test]
fn instantaneous_couple_does_not_advance_gravity_or_elastic_force() {
    let gate = CancelGate::new();
    let source = SpatialMaterialRod::try_straight(0.8, vec![section(1.0)])
        .unwrap()
        .with_gravity([0.0, 0.0, -9.0], &gate)
        .unwrap();
    let model = Arc::new(
        source
            .into_motion_model(RodRootBoundary::Free, &gate)
            .unwrap(),
    );
    let mut pose = model.initial_pose();
    pose.positions_m[1][0] += 0.008;
    let state = RodMotion::new(
        model,
        pose.clone(),
        &[[0.0; 3]; 2],
        &[[0.0; 3]; 2],
        3.0,
        &gate,
    )
    .unwrap();
    let impulse = RodPointImpulse {
        node: 1,
        point_body_m: [0.0; 3],
        linear_world_n_s: [0.0; 3],
        couple_world_nm_s: [1e-4, 0.0, 0.0],
    };
    let (next, report) = state
        .apply_point_impulses(&[impulse], impulse_controls([0.0; 3]), &gate)
        .unwrap();
    assert_eq!(next.pose(), &pose);
    assert_eq!(next.linear_momenta_kg_m_s(), &[[0.0; 3]; 2]);
    vector_near(
        next.body_spins_kg_m2_s()[1],
        impulse.couple_world_nm_s,
        1e-15,
    );
    assert_eq!(
        report.before.potential_energy_j,
        report.after.potential_energy_j
    );
    assert_eq!(next.time_s(), 3.0);
    assert_eq!(next.accepted_steps(), 0);
}

#[test]
fn invalid_late_impulse_overflow_and_cancellation_leave_the_trajectory_unchanged() {
    let gate = CancelGate::new();
    let model = axial_model(false);
    let state = RodMotion::new(
        model.clone(),
        model.initial_pose(),
        &[[0.0; 3]; 2],
        &[[0.0; 3]; 2],
        0.0,
        &gate,
    )
    .unwrap();
    let impulse = RodPointImpulse {
        node: 1,
        point_body_m: [0.0; 3],
        linear_world_n_s: [0.001, 0.0, 0.0],
        couple_world_nm_s: [0.0; 3],
    };
    let before = state.observe([0.0; 3]).unwrap();
    let (accepted, _) = state
        .apply_point_impulses(&[impulse], impulse_controls([0.0; 3]), &gate)
        .unwrap();
    for bad in [
        RodPointImpulse { node: 2, ..impulse },
        RodPointImpulse {
            point_body_m: [f64::NAN, 0.0, 0.0],
            ..impulse
        },
        RodPointImpulse {
            point_body_m: [f64::MAX, f64::MAX, 0.0],
            linear_world_n_s: [1.0, -1.0, 0.0],
            ..impulse
        },
        RodPointImpulse {
            linear_world_n_s: [f64::MAX, 0.0, 0.0],
            ..impulse
        },
    ] {
        assert!(
            state
                .apply_point_impulses(&[impulse, bad], impulse_controls([0.0; 3]), &gate)
                .is_err()
        );
        assert_eq!(state.observe([0.0; 3]).unwrap(), before);
    }
    let (replay, _) = state
        .apply_point_impulses(&[impulse], impulse_controls([0.0; 3]), &gate)
        .unwrap();
    assert_eq!(
        accepted.linear_momenta_kg_m_s(),
        replay.linear_momenta_kg_m_s()
    );
    assert_eq!(accepted.body_spins_kg_m2_s(), replay.body_spins_kg_m2_s());
    for tolerance in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let controls = impulse_controls([0.0; 3]);
        for bad in [
            RodImpulseControls {
                linear_tolerance_n_s: tolerance,
                ..controls
            },
            RodImpulseControls {
                angular_tolerance_nm_s: tolerance,
                ..controls
            },
            RodImpulseControls {
                work_tolerance_j: tolerance,
                ..controls
            },
        ] {
            assert!(state.apply_point_impulses(&[impulse], bad, &gate).is_err());
        }
    }
    gate.request();
    assert!(
        state
            .apply_point_impulses(&[impulse], impulse_controls([0.0; 3]), &gate)
            .is_err()
    );
    assert_eq!(state.observe([0.0; 3]).unwrap(), before);
}

#[test]
fn lost_parallel_impulses_refuse_while_small_resolved_transverse_work_is_retained() {
    let gate = CancelGate::new();
    let model = Arc::new(
        SpatialMaterialRod::try_straight(0.8, vec![section(1.0)])
            .unwrap()
            .into_motion_model(RodRootBoundary::Free, &gate)
            .unwrap(),
    );
    let state = RodMotion::new(
        model.clone(),
        model.initial_pose(),
        &[[1e20, 0.0, 0.0]; 2],
        &[[0.0; 3]; 2],
        0.0,
        &gate,
    )
    .unwrap();
    let before = state.observe([0.0; 3]).unwrap();
    let transverse = RodPointImpulse {
        node: 1,
        point_body_m: [0.0; 3],
        linear_world_n_s: [0.0, 0.001, 0.0],
        couple_world_nm_s: [0.0; 3],
    };
    let (resolved, report) = state
        .apply_point_impulses(&[transverse], impulse_controls([0.0; 3]), &gate)
        .unwrap();
    let exact = 0.5e-6 / model.inertia()[1].mass_kg();
    near(report.kinetic_energy_change_j, exact, 1e-20);
    near(report.applied_work_j, exact, 1e-20);
    near(report.energy_residual_j, 0.0, 1e-20);
    // Total kinetic storage cannot represent this increment, but the actual
    // transverse momentum and its stable quadratic increment are resolved.
    assert_eq!(
        report.after.kinetic_energy_j,
        report.before.kinetic_energy_j
    );
    near(resolved.linear_momenta_kg_m_s()[1][1], 0.001, 0.0);
    let parallel = RodPointImpulse {
        linear_world_n_s: [0.001, 0.0, 0.0],
        ..transverse
    };
    let opposite = RodPointImpulse {
        node: 0,
        linear_world_n_s: [-0.001, 0.0, 0.0],
        ..parallel
    };
    for batch in [vec![parallel], vec![parallel, opposite]] {
        assert!(matches!(
            state.apply_point_impulses(&batch, impulse_controls([0.0; 3]), &gate),
            Err(RodMotionError::ImpulseBalance { linear_defect_n_s, .. }) if linear_defect_n_s >= 0.001
        ));
    }
    assert_eq!(state.observe([0.0; 3]).unwrap(), before);
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(state.accepted_steps(), 0);
}
