//! G0/G1/G3/G4: shared finite inertia, joint contact law, and native continuation.
use super::*;
use fs_contact::multicontact::{MultiImpactConfig, MultiImpactDeterminacy, MultiImpactError};
use fs_couple::rod_motion::rigid_impact::multiple::RodRigidMultiImpactError;

fn pair_witnesses(restitution: f64) -> [RodRigidNormalContact; 2] {
    std::array::from_fn(|node| RodRigidNormalContact {
        node,
        point_body_m: [0.0; 3],
        rigid_arm_body_m: [if node == 0 { -0.4 } else { 0.4 }, 0.0, 0.0],
        normal_world: [0.0, 1.0, 0.0],
        activation_tolerance_m: 1e-8,
        restitution,
    })
}

#[test]
fn two_attachments_share_body_rotation_and_match_closed_form_simultaneous_collision() {
    let gate = CancelGate::new();
    let rod = free_rod([0.0; 3]);
    let m = rod.model().inertia()[0].mass_kg();
    let mass = inertia(2.3);
    let rigid = body(mass, [0.4, 0.0, 0.0], [0.0, 1.7, 0.0]);
    for restitution in [0.0, 0.4, 1.0] {
        let contacts = pair_witnesses(restitution);
        let (next, finite, report) = rod
            .apply_rigid_normal_impacts(
                rigid,
                mass,
                &contacts,
                impact_controls(),
                MultiImpactConfig::default(),
                &gate,
            )
            .unwrap();
        let impulse = (1.0 + restitution) * 1.7 / (1.0 / m + 2.0 / mass.mass());
        for &p in &report.contact.normal_impulses {
            near(p, impulse, 3e-12);
        }
        for v in next.velocities_m_s().unwrap() {
            near(v[1], impulse / m, 3e-12);
        }
        near(
            finite.center_of_mass_velocity_world(mass).unwrap().y,
            1.7 - 2.0 * impulse / mass.mass(),
            3e-12,
        );
        vector_near(array(finite.angular_momentum_body()), [0.0; 3], 3e-12);
        let reduced = 2.0 * m * mass.mass() / (2.0 * m + mass.mass());
        let loss = 0.5 * reduced * (1.0 - restitution * restitution) * 1.7 * 1.7;
        near(report.contact.dissipated_energy, loss, 3e-12);
        let balance = report.balance.unwrap();
        near(balance.kinetic_energy_change_j, -loss, 3e-12);
        vector_near(balance.linear_momentum_residual_n_s, [0.0; 3], 3e-12);
        vector_near(balance.angular_momentum_residual_nm_s, [0.0; 3], 3e-12);
        assert_eq!(next.pose(), rod.pose());
        assert_eq!(finite.pose(), rigid.pose());
        assert_eq!(next.time_s(), rod.time_s());
        assert_eq!(next.accepted_steps(), rod.accepted_steps());

        let reversed = rod
            .apply_rigid_normal_impacts(
                rigid,
                mass,
                &[contacts[1], contacts[0]],
                impact_controls(),
                MultiImpactConfig::default(),
                &gate,
            )
            .unwrap();
        for (a, b) in next
            .linear_momenta_kg_m_s()
            .iter()
            .zip(reversed.0.linear_momenta_kg_m_s())
        {
            vector_near(*a, *b, 3e-12);
        }
        vector_near(
            array(finite.linear_momentum_world()),
            array(reversed.1.linear_momentum_world()),
            3e-12,
        );

        // Independent sequential single-pair events do not solve this physical
        // simultaneous event: the first impulse rotates the body seen by the second.
        let first = rod
            .apply_rigid_normal_impact(rigid, mass, contacts[0], impact_controls(), &gate)
            .unwrap();
        let sequential = first
            .0
            .apply_rigid_normal_impact(first.1, mass, contacts[1], impact_controls(), &gate)
            .unwrap();
        assert!(
            (sequential.1.angular_momentum_body().z - finite.angular_momentum_body().z).abs()
                > 1e-3
        );
        let continued = next
            .step_candidate(1e-5, &[RodNodalWrench::default(); 2], controls(1e-5), &gate)
            .unwrap()
            .0;
        assert_eq!(continued.accepted_steps(), rod.accepted_steps() + 1);
    }
}

fn mixed_speeds() -> (RodMotion, RigidBodyState, MassProperties) {
    let base = free_rod([0.0; 3]);
    let model = Arc::new(base.model().clone());
    let rod = RodMotion::new(
        model.clone(),
        model.initial_pose(),
        &[[0.0; 3], [0.0, 1.1, 0.0]],
        &[[0.0; 3]; 2],
        0.0,
        &CancelGate::new(),
    )
    .unwrap();
    let mass = inertia(2.3);
    (rod, body(mass, [0.4, 0.0, 0.0], [0.0, 1.0, 0.0]), mass)
}

#[test]
fn initially_separating_active_contact_still_participates_in_the_joint_response() {
    let (rod, rigid, mass) = mixed_speeds();
    let m = rod.model().inertia()[0].mass_kg();
    let d = 1.0 / m + 1.0 / mass.mass() + 0.16 / mass.principal_inertia_body().z;
    let c = 1.0 / mass.mass() - 0.16 / mass.principal_inertia_body().z;
    let mut independent = [1.0, -0.1];
    fs_la::factor::lu(&[d, c, c, d], 2)
        .unwrap()
        .solve(&mut independent);
    assert!(independent.iter().all(|&p| p > 0.0));
    let (_, _, report) = rod
        .apply_rigid_normal_impacts(
            rigid,
            mass,
            &pair_witnesses(0.0),
            impact_controls(),
            MultiImpactConfig::default(),
            &CancelGate::new(),
        )
        .unwrap();
    for i in 0..2 {
        near(report.contact.normal_impulses[i], independent[i], 3e-12);
        near(report.relative_normal_speeds_after_m_s[i], 0.0, 3e-12);
    }
}

#[test]
fn shared_eccentric_node_and_full_tensors_match_independent_spatial_mass_response() {
    let gate = CancelGate::new();
    let model = Arc::new(
        eccentric_source()
            .into_spatial_motion_model(RodRootBoundary::Free, &gate)
            .unwrap(),
    );
    let nodal = model.inertia()[1];
    let tensor = nodal.inertia_com_body_kg_m2();
    let center = nodal.center_of_mass_body_m();
    let rod_mass = fs_mbd::articulated::SpatialInertia::new(
        nodal.mass_kg(),
        Vec3::new(center[0], center[1], center[2]),
        Mat3 {
            m: std::array::from_fn(|i| tensor[i / 3][i % 3]),
        },
    )
    .unwrap()
    .matrix()
    .unwrap();
    let body_tensor = Mat3 {
        m: [0.07, 0.01, -0.004, 0.01, 0.09, 0.006, -0.004, 0.006, 0.08],
    };
    let mass = MassProperties::from_inertia_tensor(2.4, BodyVec::ZERO, body_tensor).unwrap();
    let body_center = [0.65, 0.1, -0.1];
    let points = [[0.04, 0.06, -0.03], [-0.02, -0.05, 0.04]];
    let normals = [[1.0, 0.0, 0.0], [0.6, 0.8, 0.0]];
    let arms: [[f64; 3]; 2] = std::array::from_fn(|i| {
        sub(
            sub(add([0.8, 0.0, 0.0], points[i]), body_center),
            normals[i].map(|x| x * 1e-9),
        )
    });
    let mut rows = [[0.0; 6]; 2];
    let mut response = [[0.0; 6]; 2];
    let mut body_rows = [[0.0; 3]; 2];
    let mut body_response = [[0.0; 3]; 2];
    for i in 0..2 {
        let torque = cross(points[i], normals[i]);
        rows[i] = [
            torque[0],
            torque[1],
            torque[2],
            normals[i][0],
            normals[i][1],
            normals[i][2],
        ];
        response[i] = rows[i];
        fs_la::factor::lu(&rod_mass.m, 6)
            .unwrap()
            .solve(&mut response[i]);
        body_rows[i] = cross(arms[i], normals[i]);
        body_response[i] = body_rows[i];
        fs_la::factor::lu(&body_tensor.m, 3)
            .unwrap()
            .solve(&mut body_response[i]);
    }
    let mut mobility = [0.0; 4];
    for i in 0..2 {
        for j in 0..2 {
            mobility[2 * i + j] = (0..6).map(|k| rows[i][k] * response[j][k]).sum::<f64>()
                + (0..3)
                    .map(|k| {
                        normals[i][k] * normals[j][k] / mass.mass()
                            + body_rows[i][k] * body_response[j][k]
                    })
                    .sum::<f64>();
        }
    }
    let exact = [0.17, 0.23];
    let closing: [f64; 2] =
        std::array::from_fn(|i| mobility[2 * i] * exact[0] + mobility[2 * i + 1] * exact[1]);
    let initial_velocity = [closing[0], (closing[1] - 0.6 * closing[0]) / 0.8, 0.0];
    for (q, translation, boost) in [
        (So3::identity(), [0.0; 3], [0.0; 3]),
        (
            So3::exp(So3Tangent::new(Vec3::new(0.3, -0.2, 0.4))).unwrap(),
            [0.7, -0.5, 0.2],
            [0.3, -0.2, 0.1],
        ),
    ] {
        let mut pose = model.initial_pose();
        for i in 0..2 {
            pose.positions_m[i] = add(rotate(q, pose.positions_m[i]), translation);
            pose.frames[i] = q;
        }
        let rod = SpatialRodMotion::new_spatial(
            model.clone(),
            pose,
            &[boost; 2],
            &[[0.0; 3]; 2],
            0.0,
            1e-12,
            &gate,
        )
        .unwrap();
        let rigid = RigidBodyState::new(
            Pose::new(
                vector(add(rotate(q, body_center), translation)),
                orientation(q),
            )
            .unwrap(),
            vector(add(rotate(q, initial_velocity), boost)).scale(mass.mass()),
            BodyVec::ZERO,
        )
        .unwrap();
        let witnesses: [RodRigidNormalContact; 2] =
            std::array::from_fn(|i| RodRigidNormalContact {
                node: 1,
                point_body_m: points[i],
                rigid_arm_body_m: arms[i],
                normal_world: rotate(q, normals[i]),
                activation_tolerance_m: 1e-8,
                restitution: 0.0,
            });
        let (_, _, report) = rod
            .apply_rigid_normal_impacts(
                rigid,
                mass,
                &witnesses,
                impact_controls(),
                MultiImpactConfig::default(),
                &gate,
            )
            .unwrap();
        for i in 0..2 {
            near(report.contact.normal_impulses[i], exact[i], 3e-11);
        }
        for i in 0..4 {
            near(report.delassus_kg_inverse[i], mobility[i], 3e-11);
        }
        let balance = report.balance.unwrap();
        vector_near(balance.linear_momentum_residual_n_s, [0.0; 3], 3e-11);
        vector_near(balance.angular_momentum_residual_nm_s, [0.0; 3], 3e-11);
        near(balance.energy_residual_j, 0.0, 3e-11);
    }
}

#[test]
fn redundant_contacts_retain_determinacy_and_actual_stationary_support_reaction() {
    let gate = CancelGate::new();
    for clamped in [false, true] {
        let source = SpatialMaterialRod::try_straight(0.8, vec![section(1.0)]).unwrap();
        let source = if clamped {
            source
                .with_supports(
                    &[RodSupport {
                        node: 1,
                        position_m: [Some(0.8), Some(0.0), Some(0.0)],
                        frame: Some(So3::identity()),
                    }],
                    &gate,
                )
                .unwrap()
        } else {
            source
        };
        let model = Arc::new(
            source
                .into_motion_model(RodRootBoundary::Free, &gate)
                .unwrap(),
        );
        let rod = RodMotion::new(
            model.clone(),
            model.initial_pose(),
            &[[0.0; 3]; 2],
            &[[0.0; 3]; 2],
            0.0,
            &gate,
        )
        .unwrap();
        let mass = inertia(2.0);
        let rigid = body(mass, [0.8, 0.07, 0.0], [1.0, 0.0, 0.0]);
        let point = RodRigidNormalContact {
            point_body_m: [0.0, 0.07, 0.0],
            ..witness(1, 0.3)
        };
        let single = rod
            .apply_rigid_normal_impact(rigid, mass, point, impact_controls(), &gate)
            .unwrap();
        let (next, finite, report) = rod
            .apply_rigid_normal_impacts(
                rigid,
                mass,
                &[point, point],
                impact_controls(),
                MultiImpactConfig::default(),
                &gate,
            )
            .unwrap();
        assert!(matches!(
            report.contact.determinacy,
            MultiImpactDeterminacy::UniqueVelocity(_)
        ));
        near(
            report.contact.normal_impulses.iter().sum(),
            single.2.contact.normal_impulse,
            3e-11,
        );
        for (a, b) in next
            .linear_momenta_kg_m_s()
            .iter()
            .zip(single.0.linear_momenta_kg_m_s())
        {
            vector_near(*a, *b, 3e-11);
        }
        for (a, b) in next
            .body_spins_kg_m2_s()
            .iter()
            .zip(single.0.body_spins_kg_m2_s())
        {
            vector_near(*a, *b, 3e-11);
        }
        vector_near(
            array(finite.linear_momentum_world()),
            array(single.1.linear_momentum_world()),
            3e-11,
        );
        if clamped {
            let impulse = 1.3 * mass.mass();
            let support = &report.rod_impulse.as_ref().unwrap().support_impulses[1];
            vector_near(support.linear_n_s, [-impulse, 0.0, 0.0], 3e-11);
            near(support.angular_nm_s[2], 0.07 * impulse, 3e-11);
            assert_eq!(next.pose(), rod.pose());
            assert_eq!(next.linear_momenta_kg_m_s(), rod.linear_momenta_kg_m_s());
        }
    }
}

#[test]
fn energy_injecting_restitution_budget_and_cancellation_refuse_without_partial_motion() {
    let gate = CancelGate::new();
    let (rod, rigid, mass) = mixed_speeds();
    let mut bad = pair_witnesses(0.0);
    bad[0].restitution = 1.0;
    assert!(matches!(
        rod.apply_rigid_normal_impacts(
            rigid,
            mass,
            &bad,
            impact_controls(),
            MultiImpactConfig::default(),
            &gate
        ),
        Err(RodRigidMultiImpactError::Contact(
            MultiImpactError::EnergyIncrease { .. }
        ))
    ));
    let mut capped = MultiImpactConfig::default();
    capped.solver.max_variables = 1;
    assert!(matches!(
        rod.apply_rigid_normal_impacts(rigid, mass, &bad, impact_controls(), capped, &gate),
        Err(RodRigidMultiImpactError::Invalid(_))
    ));
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        rod.apply_rigid_normal_impacts(
            rigid,
            mass,
            &pair_witnesses(0.0),
            impact_controls(),
            MultiImpactConfig::default(),
            &cancelled
        )
        .is_err()
    );
    let first = rod
        .apply_rigid_normal_impacts(
            rigid,
            mass,
            &pair_witnesses(0.0),
            impact_controls(),
            MultiImpactConfig::default(),
            &gate,
        )
        .unwrap();
    let retry = rod
        .apply_rigid_normal_impacts(
            rigid,
            mass,
            &pair_witnesses(0.0),
            impact_controls(),
            MultiImpactConfig::default(),
            &gate,
        )
        .unwrap();
    assert_eq!(
        first.0.linear_momenta_kg_m_s(),
        retry.0.linear_momenta_kg_m_s()
    );
    assert_eq!(first.1, retry.1);
    for witnesses in [&[][..], &pair_witnesses(0.0)[..]] {
        let still = free_rod([0.0; 3]);
        let departing = body(mass, [0.4, 0.0, 0.0], [0.0, -1.0, 0.0]);
        let (next, finite, report) = still
            .apply_rigid_normal_impacts(
                departing,
                mass,
                witnesses,
                impact_controls(),
                MultiImpactConfig::default(),
                &gate,
            )
            .unwrap();
        assert!(
            report.balance.is_none()
                && report.rod_impulse.is_none()
                && report.rigid_impulse.is_none()
        );
        assert_eq!(next.linear_momenta_kg_m_s(), still.linear_momenta_kg_m_s());
        assert_eq!(next.pose(), still.pose());
        assert_eq!(finite, departing);
    }
}
