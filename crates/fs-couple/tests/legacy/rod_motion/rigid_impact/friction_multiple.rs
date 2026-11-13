//! G0/G1/G3/G4: one coupled Coulomb law and original native momentum response.

use super::*;
use fs_contact::frictional_multicontact::{
    FrictionalContact, FrictionalContactConfig, FrictionalContactDeterminacy,
    FrictionalContactError, NormalVelocityLaw, solve_frictional_contacts,
};
use fs_couple::rod_motion::rigid_impact::friction::{
    RodRigidFrictionalContact, RodRigidFrictionalControls, RodRigidFrictionalError,
};
use fs_couple::rod_motion::rigid_impact::friction_multiple::RodRigidFrictionalMultiImpactReport;

fn bounds() -> RodRigidFrictionalControls {
    let mut impact = impact_controls();
    impact.velocity_tolerance_m_s = 2e-9;
    impact.impulse.work_tolerance_j = 2e-9;
    RodRigidFrictionalControls {
        impact,
        contact: FrictionalContactConfig {
            velocity_tolerance: 2e-11,
            impulse_tolerance: 2e-11,
            complementarity_tolerance: 2e-11,
            energy_tolerance: 2e-11,
            ..Default::default()
        },
    }
}
fn contact(mu: f64) -> RodRigidFrictionalContact {
    RodRigidFrictionalContact {
        geometry: witness(1, 0.35),
        tangent_world: [0.0, 1.0, 0.0],
        friction_coefficient: mu,
    }
}
fn audit(report: &RodRigidFrictionalMultiImpactReport, witnesses: &[RodRigidFrictionalContact]) {
    for (i, witness) in witnesses.iter().enumerate() {
        let p = report.contact.impulses[i];
        let v = report.relative_velocities_after_m_s[i];
        let mu = witness.friction_coefficient;
        assert!(p[0] >= -2e-10);
        if report.gaps_m[i] <= witness.geometry.activation_tolerance_m {
            assert!(v[0] >= report.contact.normal_targets[i] - 2e-9);
            near(p[0] * (v[0] - report.contact.normal_targets[i]), 0.0, 2e-9);
        }
        assert!(p[1].hypot(p[2]) <= mu * p[0] + 2e-9);
        let slip = v[1].hypot(v[2]);
        if slip > 2e-9 {
            near(p[1], -mu * p[0] * v[1] / slip, 2e-9);
            near(p[2], -mu * p[0] * v[2] / slip, 2e-9);
        }
        vector_near(report.endpoint_velocity_residuals_m_s[i], [0.0; 3], 2e-9);
    }
    assert!(report.contact.dissipated_energy >= 0.0);
    if let Some(balance) = report.balance {
        vector_near(balance.linear_momentum_residual_n_s, [0.0; 3], 2e-11);
        vector_near(balance.angular_momentum_residual_nm_s, [0.0; 3], 2e-11);
        near(balance.energy_residual_j, 0.0, 2e-9);
        near(balance.work_energy_residual_j, 0.0, 2e-9);
        near(
            balance.relative_impulse_work_j,
            report.contact.relative_work,
            2e-9,
        );
    }
}

#[test]
fn singleton_matches_published_friction_and_both_states_continue() {
    let gate = CancelGate::new();
    let rod = free_rod([-0.6, 0.8, -0.6]);
    let mass = inertia(2.3);
    let rigid = body(mass, [0.8, 0.0, 0.0], [1.4, 0.0, 0.0]);
    for mu in [0.0, 0.1, 0.7] {
        let witness = contact(mu);
        let single = rod
            .apply_rigid_frictional_impact(rigid, mass, witness, bounds(), &gate)
            .unwrap();
        let (next, finite, report) = rod
            .apply_rigid_frictional_impacts(rigid, mass, &[witness], bounds(), &gate)
            .unwrap();
        vector_near(
            report.impulses_world_n_s[0],
            single.2.impulse_world_n_s,
            2e-10,
        );
        near(
            report.contact.dissipated_energy,
            single.2.contact.dissipated_energy,
            2e-10,
        );
        for node in 0..2 {
            vector_near(
                next.linear_momenta_kg_m_s()[node],
                single.0.linear_momenta_kg_m_s()[node],
                2e-10,
            );
            vector_near(
                next.body_spins_kg_m2_s()[node],
                single.0.body_spins_kg_m2_s()[node],
                2e-10,
            );
        }
        vector_near(
            array(finite.linear_momentum_world()),
            array(single.1.linear_momentum_world()),
            2e-10,
        );
        vector_near(
            array(finite.angular_momentum_body()),
            array(single.1.angular_momentum_body()),
            2e-10,
        );
        assert_eq!(next.pose(), rod.pose());
        assert_eq!(finite.pose(), rigid.pose());
        assert_eq!(next.time_s(), rod.time_s());
        assert_eq!(next.accepted_steps(), rod.accepted_steps());
        assert!(std::ptr::eq(next.model(), rod.model()));
        audit(&report, &[witness]);
        let continued = next
            .step_candidate(1e-6, &[RodNodalWrench::default(); 2], controls(1e-6), &gate)
            .unwrap()
            .0;
        assert_eq!(continued.accepted_steps(), rod.accepted_steps() + 1);
        let continued_body = fs_mbd::RigidBodyIntegrator::new(fs_mbd::Gravity::ZERO)
            .step(finite, mass, fs_mbd::Wrench::ZERO, 1e-6)
            .unwrap();
        assert_eq!(
            continued_body.state_after.linear_momentum_world(),
            finite.linear_momentum_world()
        );
        assert_ne!(
            continued_body.state_after.pose().position_world(),
            finite.pose().position_world()
        );
    }
}

fn pair_contacts(mu: f64, restitution: f64) -> [RodRigidFrictionalContact; 2] {
    std::array::from_fn(|node| RodRigidFrictionalContact {
        geometry: RodRigidNormalContact {
            node,
            point_body_m: [0.0; 3],
            rigid_arm_body_m: [if node == 0 { -0.4 } else { 0.4 }, 0.0, 0.0],
            normal_world: [0.0, 1.0, 0.0],
            activation_tolerance_m: 1e-8,
            restitution,
        },
        tangent_world: [1.0, 0.0, 0.0],
        friction_coefficient: mu,
    })
}

#[test]
fn two_nodes_share_finite_body_mobility_in_stick_and_slip_without_impact_order() {
    let gate = CancelGate::new();
    let rod = free_rod([0.0; 3]);
    let mass = inertia(2.3);
    let rigid = body(mass, [0.4, 0.0, 0.0], [0.8, 1.7, -0.6]);
    let m = rod.model().inertia()[0].mass_kg();
    let inverse_mass = 1.0 / m + 2.0 / mass.mass();
    for mu in [0.0, 0.1, 0.7] {
        let contacts = pair_contacts(mu, 0.4);
        let (next, finite, report) = rod
            .apply_rigid_frictional_impacts(rigid, mass, &contacts, bounds(), &gate)
            .unwrap();
        let normal = 1.4 * 1.7 / inverse_mass;
        let tangent = (1.0 / inverse_mass).min(mu * normal);
        let expected = [0.8 * tangent, normal, -0.6 * tangent];
        for i in 0..2 {
            vector_near(report.impulses_world_n_s[i], expected, 2e-9);
            vector_near(
                next.velocities_m_s().unwrap()[i],
                expected.map(|p| p / m),
                2e-9,
            );
            near(report.relative_velocities_after_m_s[i][0], 0.4 * 1.7, 2e-9);
        }
        vector_near(array(finite.angular_momentum_body()), [0.0; 3], 2e-9);
        let loss = (1.0 - 0.4 * 0.4) * 1.7 * 1.7 / inverse_mass + 2.0 * tangent
            - inverse_mass * tangent * tangent;
        near(report.contact.dissipated_energy, loss, 2e-9);
        let reversed = rod
            .apply_rigid_frictional_impacts(
                rigid,
                mass,
                &[contacts[1], contacts[0]],
                bounds(),
                &gate,
            )
            .unwrap();
        for i in 0..2 {
            vector_near(
                next.linear_momenta_kg_m_s()[i],
                reversed.0.linear_momenta_kg_m_s()[i],
                2e-9,
            );
        }
        vector_near(
            array(finite.linear_momentum_world()),
            array(reversed.1.linear_momentum_world()),
            2e-9,
        );
        audit(&report, &contacts);
        // A sequential pair update changes the shared body's spin before the
        // second contact and fails this simultaneous physical response.
        let first = rod
            .apply_rigid_frictional_impact(rigid, mass, contacts[0], bounds(), &gate)
            .unwrap();
        let sequential = first
            .0
            .apply_rigid_frictional_impact(first.1, mass, contacts[1], bounds(), &gate)
            .unwrap();
        assert!(
            (sequential.1.angular_momentum_body().z - finite.angular_momentum_body().z).abs()
                > 1e-3
        );
    }
}

#[test]
fn full_six_dimensional_origin_mass_oracle_covers_shared_nodes_rotations_and_tangent_frames() {
    let gate = CancelGate::new();
    let model = Arc::new(
        eccentric_source()
            .into_spatial_motion_model(RodRootBoundary::Free, &gate)
            .unwrap(),
    );
    let body_tensor = Mat3 {
        m: [0.07, 0.01, -0.004, 0.01, 0.09, 0.006, -0.004, 0.006, 0.08],
    };
    let mass = MassProperties::from_inertia_tensor(2.4, BodyVec::ZERO, body_tensor).unwrap();
    let body_matrix = fs_mbd::articulated::SpatialInertia::new(
        mass.mass(),
        Vec3::new(0.0, 0.0, 0.0),
        body_tensor,
    )
    .unwrap()
    .matrix()
    .unwrap();
    let body_factor = fs_la::factor::lu(&body_matrix.m, 6).unwrap();
    // Equal x offsets keep the two-point redundant force gauge in the tangent
    // plane. Normal loads remain observable through the rotational response.
    let points = [[0.03, 0.05, -0.04], [0.03, -0.05, 0.04]];
    let center = [0.4, 0.07, -0.03];
    for (nodes, sliding) in [([0, 1], false), ([0, 1], true), ([1, 1], false)] {
        let arms: [[f64; 3]; 2] = std::array::from_fn(|i| {
            sub(
                add([0.8 * nodes[i] as f64, 0.0, 0.0], points[i]),
                add(center, [1e-13, 0.0, 0.0]),
            )
        });
        let mut rod_rows = [[0.0; 6]; 6];
        let mut rod_response = [[0.0; 6]; 6];
        let mut body_rows = [[0.0; 6]; 6];
        let mut body_response = [[0.0; 6]; 6];
        for row in 0..6 {
            let contact = row / 3;
            let nodal = model.inertia()[nodes[contact]];
            let tensor = nodal.inertia_com_body_kg_m2();
            let centroid = nodal.center_of_mass_body_m();
            let matrix = fs_mbd::articulated::SpatialInertia::new(
                nodal.mass_kg(),
                Vec3::new(centroid[0], centroid[1], centroid[2]),
                Mat3 {
                    m: std::array::from_fn(|c| tensor[c / 3][c % 3]),
                },
            )
            .unwrap()
            .matrix()
            .unwrap();
            let axis: [f64; 3] = std::array::from_fn(|c| if c == row % 3 { 1.0 } else { 0.0 });
            let torque = cross(points[contact], axis);
            rod_rows[row] = [torque[0], torque[1], torque[2], axis[0], axis[1], axis[2]];
            rod_response[row] = rod_rows[row];
            fs_la::factor::lu(&matrix.m, 6)
                .unwrap()
                .solve(&mut rod_response[row]);
            let torque = cross(arms[contact], axis);
            body_rows[row] = [torque[0], torque[1], torque[2], axis[0], axis[1], axis[2]];
            body_response[row] = body_rows[row];
            body_factor.solve(&mut body_response[row]);
        }
        let mut mobility = [0.0; 36];
        for i in 0..6 {
            for j in 0..6 {
                mobility[6 * i + j] = (0..6)
                    .map(|k| {
                        body_rows[i][k] * body_response[j][k]
                            + if nodes[i / 3] == nodes[j / 3] {
                                rod_rows[i][k] * rod_response[j][k]
                            } else {
                                0.0
                            }
                    })
                    .sum();
            }
        }
        assert!(mobility[3..6].iter().any(|x| x.abs() > 1e-3));
        assert!((mobility[1] - mobility[2]).abs() > 1e-3);
        let mu = 0.2;
        let impulses: [[f64; 3]; 2] = [0.18, 0.24].map(|pn| {
            if sliding {
                [pn, -mu * pn * 0.6, -mu * pn * 0.8]
            } else {
                [pn, -0.1 * pn, 0.05 * pn]
            }
        });
        let final_velocity = if sliding { [0.0, 1.2, 1.6] } else { [0.0; 3] };
        let mut rod_before_v = [final_velocity; 2];
        let mut rod_before_w = [[0.0; 3]; 2];
        let mut body_before_v = [0.0; 3];
        let mut body_before_w = [0.0; 3];
        for row in 0..6 {
            let p = impulses[row / 3][row % 3];
            let node = nodes[row / 3];
            for c in 0..3 {
                rod_before_w[node][c] -= rod_response[row][c] * p;
                rod_before_v[node][c] -= rod_response[row][c + 3] * p;
                body_before_w[c] += body_response[row][c] * p;
                body_before_v[c] += body_response[row][c + 3] * p;
            }
        }
        let q = So3::exp(So3Tangent::new(Vec3::new(0.3, -0.2, 0.4))).unwrap();
        for (q, shift, boost, hints, order) in [
            (
                So3::identity(),
                [0.0; 3],
                [0.0; 3],
                [[0.0, 1.0, 0.0]; 2],
                [0, 1],
            ),
            (
                q,
                [0.7, -0.5, 0.2],
                [0.3, -0.2, 0.1],
                [[0.0, 0.6, 0.8], [0.0, 0.8, -0.6]],
                [1, 0],
            ),
        ] {
            let mut pose = model.initial_pose();
            for node in 0..2 {
                pose.positions_m[node] = add(rotate(q, pose.positions_m[node]), shift);
                pose.frames[node] = q;
            }
            let rod = SpatialRodMotion::new_spatial(
                model.clone(),
                pose,
                &rod_before_v.map(|v| add(rotate(q, v), boost)),
                &rod_before_w,
                0.0,
                1e-12,
                &gate,
            )
            .unwrap();
            let rigid = RigidBodyState::new(
                Pose::new(vector(add(rotate(q, center), shift)), orientation(q)).unwrap(),
                vector(add(rotate(q, body_before_v), boost)).scale(mass.mass()),
                mass.angular_momentum_body_checked(vector(body_before_w))
                    .unwrap(),
            )
            .unwrap();
            let contacts = order.map(|i| RodRigidFrictionalContact {
                geometry: RodRigidNormalContact {
                    node: nodes[i],
                    point_body_m: points[i],
                    rigid_arm_body_m: arms[i],
                    normal_world: rotate(q, [1.0, 0.0, 0.0]),
                    activation_tolerance_m: 1e-8,
                    restitution: 0.0,
                },
                tangent_world: rotate(q, hints[i]),
                friction_coefficient: mu,
            });
            let (next, finite, report) = rod
                .apply_rigid_frictional_impacts(rigid, mass, &contacts, bounds(), &gate)
                .unwrap();
            // Independently transform the full origin-mass oracle into each
            // caller tangent frame, including both off-diagonal contact blocks.
            for i in 0..6 {
                let a = rotate(q.inverse(), report.frames_world[i / 3][i % 3]);
                for j in 0..6 {
                    let b = rotate(q.inverse(), report.frames_world[j / 3][j % 3]);
                    let mut expected = 0.0;
                    for r in 0..3 {
                        for s in 0..3 {
                            expected += a[r]
                                * mobility[(3 * order[i / 3] + r) * 6 + 3 * order[j / 3] + s]
                                * b[s];
                        }
                    }
                    near(report.delassus_kg_inverse[6 * i + j], expected, 2e-11);
                }
            }
            for node in 0..2 {
                vector_near(
                    next.velocities_m_s().unwrap()[node],
                    add(rotate(q, final_velocity), boost),
                    2e-9,
                );
                vector_near(
                    next.angular_velocities_body_rad_s().unwrap()[node],
                    [0.0; 3],
                    2e-9,
                );
            }
            vector_near(
                array(finite.center_of_mass_velocity_world(mass).unwrap()),
                boost,
                2e-9,
            );
            vector_near(array(finite.angular_momentum_body()), [0.0; 3], 2e-9);
            if nodes[0] != nodes[1] {
                assert!(matches!(
                    report.contact.determinacy,
                    FrictionalContactDeterminacy::UniqueImpulses
                ));
                for i in 0..2 {
                    vector_near(
                        rotate(q.inverse(), report.impulses_world_n_s[i]),
                        impulses[order[i]],
                        2e-9,
                    );
                }
            } else {
                let FrictionalContactDeterminacy::UniqueVelocity(ref witness) =
                    report.contact.determinacy
                else {
                    panic!(
                        "two original points on one relative rigid motion have a tangent load gauge"
                    );
                };
                assert_eq!((witness.dimension, witness.rank), (6, 5));
                assert!(report.contact.impulse_work_witness.is_some());
            }
            audit(&report, &contacts);
        }
    }
}

#[test]
fn supported_and_redundant_sets_preserve_unique_motion_without_assigning_local_loads() {
    let gate = CancelGate::new();
    let source = SpatialMaterialRod::try_straight(0.8, vec![section(1.0)])
        .unwrap()
        .with_supports(
            &[RodSupport {
                node: 1,
                position_m: [Some(0.8), Some(0.0), Some(0.0)],
                frame: Some(So3::identity()),
            }],
            &gate,
        )
        .unwrap();
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
    let rigid = body(mass, [0.8, 0.0, 0.0], [1.0, -1.0, 0.0]);
    let contacts = [-0.07, 0.07].map(|arm| RodRigidFrictionalContact {
        geometry: RodRigidNormalContact {
            point_body_m: [0.0, arm, 0.0],
            rigid_arm_body_m: [0.0, arm, 0.0],
            ..witness(1, 0.0)
        },
        ..contact(0.2)
    });
    let (next, finite, report) = rod
        .apply_rigid_frictional_impacts(rigid, mass, &contacts, bounds(), &gate)
        .unwrap();
    assert_eq!(next.linear_momenta_kg_m_s(), rod.linear_momenta_kg_m_s());
    assert_eq!(next.body_spins_kg_m2_s(), rod.body_spins_kg_m2_s());
    vector_near(
        array(finite.center_of_mass_velocity_world(mass).unwrap()),
        [0.0, -0.8, 0.0],
        2e-9,
    );
    vector_near(array(finite.angular_momentum_body()), [0.0; 3], 2e-9);
    let support = &report.rod_impulse.as_ref().unwrap().support_impulses[1];
    vector_near(support.linear_n_s, [-2.0, 0.4, 0.0], 2e-9);
    vector_near(support.angular_nm_s, [0.0; 3], 2e-9);
    assert!(matches!(
        report.contact.determinacy,
        FrictionalContactDeterminacy::UniqueVelocity(_)
    ));
    audit(&report, &contacts);

    let rod = free_rod([-0.6, 0.8, -0.6]);
    let mass = inertia(2.3);
    let rigid = body(mass, [0.8, 0.0, 0.0], [1.4, 0.0, 0.0]);
    let duplicate = [contact(0.1); 2];
    let one = rod
        .apply_rigid_frictional_impacts(rigid, mass, &duplicate[..1], bounds(), &gate)
        .unwrap();
    let repeated = rod
        .apply_rigid_frictional_impacts(rigid, mass, &duplicate, bounds(), &gate)
        .unwrap();
    let FrictionalContactDeterminacy::UniqueVelocity(ref witness) = repeated.2.contact.determinacy
    else {
        panic!("duplicate physical contacts cannot determine local loads uniquely");
    };
    assert_eq!((witness.dimension, witness.rank), (6, 3));
    assert!(repeated.2.contact.impulse_work_witness.is_some());
    vector_near(
        repeated.0.linear_momenta_kg_m_s()[1],
        one.0.linear_momenta_kg_m_s()[1],
        2e-9,
    );
    vector_near(
        array(repeated.1.linear_momentum_world()),
        array(one.1.linear_momentum_world()),
        2e-9,
    );
    near(
        repeated.2.contact.dissipated_energy,
        one.2.contact.dissipated_energy,
        2e-9,
    );
    audit(&repeated.2, &duplicate);

    // Different friction laws on a redundant normal load may determine
    // different native motions. Preserve the original law's explicit refusal.
    let inconsistent = [contact(0.1), contact(0.2)];
    let original_contacts: Vec<_> = inconsistent
        .iter()
        .enumerate()
        .map(|(i, c)| FrictionalContact {
            id: i as u64,
            gap: 0.0,
            activation_tolerance: c.geometry.activation_tolerance_m,
            relative_velocity: repeated.2.relative_velocities_before_m_s[i],
            friction_coefficient: c.friction_coefficient,
            normal_law: NormalVelocityLaw::NewtonRestitution {
                coefficient: c.geometry.restitution,
            },
        })
        .collect();
    let original_error = solve_frictional_contacts(
        &original_contacts,
        &repeated.2.delassus_kg_inverse,
        bounds().contact,
        || false,
    )
    .unwrap_err();
    assert!(matches!(
        original_error,
        FrictionalContactError::Nonunique { .. }
            | FrictionalContactError::UnresolvedDeterminacy { .. }
    ));
    let adapter_error = rod
        .apply_rigid_frictional_impacts(rigid, mass, &inconsistent, bounds(), &gate)
        .unwrap_err();
    let RodRigidFrictionalError::Contact(adapter_error) = adapter_error else {
        panic!("preserve original law refusal");
    };
    assert_eq!(
        std::mem::discriminant(&adapter_error),
        std::mem::discriminant(&original_error)
    );
}

#[test]
fn initially_separating_and_inactive_points_keep_their_shared_body_response() {
    let gate = CancelGate::new();
    let base = free_rod([0.0; 3]);
    let model = Arc::new(base.model().clone());
    let rod = RodMotion::new(
        model.clone(),
        model.initial_pose(),
        &[[0.0; 3], [0.0, 1.1, 0.0]],
        &[[0.0; 3]; 2],
        0.0,
        &gate,
    )
    .unwrap();
    let mass = inertia(2.3);
    let rigid = body(mass, [0.4, 0.0, 0.0], [0.1, 1.0, 0.0]);
    let contacts = pair_contacts(0.05, 0.0);
    let (_, _, report) = rod
        .apply_rigid_frictional_impacts(rigid, mass, &contacts, bounds(), &gate)
        .unwrap();
    assert!(report.relative_velocities_before_m_s[1][0] > 0.0);
    assert!(report.contact.impulses[1][0] > 0.0);
    for v in &report.relative_velocities_after_m_s {
        near(v[0], 0.0, 2e-9);
    }
    audit(&report, &contacts);

    let mut inactive = contacts;
    inactive[1].geometry.rigid_arm_body_m[1] = -1e-4;
    let (_, _, report) = rod
        .apply_rigid_frictional_impacts(rigid, mass, &inactive, bounds(), &gate)
        .unwrap();
    assert_eq!(report.contact.active_indices, vec![0]);
    assert_eq!(report.contact.impulses[1], [0.0; 3]);
    assert!(
        (report.relative_velocities_after_m_s[1][0] - report.relative_velocities_before_m_s[1][0])
            .abs()
            > 1e-3
    );
    audit(&report, &inactive);
}

#[test]
fn bounded_cancellation_noops_and_late_geometry_refusal_are_atomic_and_retryable() {
    let gate = CancelGate::new();
    let rod = free_rod([-0.6, 0.8, -0.6]);
    let mass = inertia(2.3);
    let rigid = body(mass, [0.8, 0.0, 0.0], [1.4, 0.0, 0.0]);
    let witnesses = [contact(0.1); 2];
    let original_p = rod.linear_momenta_kg_m_s().to_vec();
    let original_spin = rod.body_spins_kg_m2_s().to_vec();
    for (other, contacts) in [
        (rigid, &witnesses[..0]),
        (
            body(mass, [0.79, 0.0, 0.0], [1.4, 0.0, 0.0]),
            &witnesses[..],
        ),
        (
            body(mass, [0.8, 0.0, 0.0], [-1.4, 0.0, 0.0]),
            &witnesses[..],
        ),
    ] {
        let (next, finite, report) = rod
            .apply_rigid_frictional_impacts(other, mass, contacts, bounds(), &gate)
            .unwrap();
        assert_eq!(next.linear_momenta_kg_m_s(), original_p);
        assert_eq!(next.body_spins_kg_m2_s(), original_spin);
        assert_eq!(next.pose(), rod.pose());
        assert_eq!(finite, other);
        assert!(
            report.rod_impulse.is_none()
                && report.rigid_impulse.is_none()
                && report.balance.is_none()
        );
    }
    for limit in [0, 1] {
        let mut limited = bounds();
        limited.contact.maximum_contacts = limit;
        assert!(
            rod.apply_rigid_frictional_impacts(rigid, mass, &witnesses, limited, &gate)
                .is_err()
        );
    }
    let mut scalar_budget = bounds();
    scalar_budget.contact.matrix.max_variables = 5;
    assert!(
        rod.apply_rigid_frictional_impacts(rigid, mass, &witnesses, scalar_budget, &gate)
            .is_err()
    );
    assert!(
        rod.apply_rigid_frictional_impacts(rigid, mass, &[contact(0.1); 22], bounds(), &gate)
            .is_err()
    );
    let mut evaluations = bounds();
    evaluations.contact.maximum_evaluations = 1;
    assert!(matches!(
        rod.apply_rigid_frictional_impacts(rigid, mass, &witnesses, evaluations, &gate),
        Err(RodRigidFrictionalError::Contact(
            FrictionalContactError::EvaluationBudget
        ))
    ));
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        rod.apply_rigid_frictional_impacts(rigid, mass, &witnesses, bounds(), &cancelled)
            .is_err()
    );
    let separated = body(mass, [0.8 - 1e-9, 0.0, 0.0], [1.4, 0.0, 0.0]);
    assert!(matches!(
        rod.apply_rigid_frictional_impacts(separated, mass, &witnesses, bounds(), &gate),
        Err(RodRigidFrictionalError::Impact(
            RodRigidImpactError::Invalid(_)
        ))
    ));
    let mut loose = bounds();
    loose.impact.geometry_tolerance_m = 1e-3;
    let crooked = body(mass, [0.8, 1e-4, 0.0], [1.4, 0.0, 0.0]);
    assert!(matches!(
        rod.apply_rigid_frictional_impacts(crooked, mass, &witnesses, loose, &gate),
        Err(RodRigidFrictionalError::Impact(
            RodRigidImpactError::Balance { .. }
        ))
    ));
    let penetrating = body(mass, [0.8 + 1e-4, 0.0, 0.0], [1.4, 0.0, 0.0]);
    assert!(matches!(
        rod.apply_rigid_frictional_impacts(penetrating, mass, &witnesses, bounds(), &gate),
        Err(RodRigidFrictionalError::Contact(
            FrictionalContactError::PenetratingContact { .. }
        ))
    ));
    assert_eq!(rod.linear_momenta_kg_m_s(), original_p);
    assert_eq!(rod.body_spins_kg_m2_s(), original_spin);
    let first = rod
        .apply_rigid_frictional_impacts(rigid, mass, &witnesses, bounds(), &gate)
        .unwrap();
    let retry = rod
        .apply_rigid_frictional_impacts(rigid, mass, &witnesses, bounds(), &gate)
        .unwrap();
    assert_eq!(
        first.0.linear_momenta_kg_m_s(),
        retry.0.linear_momenta_kg_m_s()
    );
    assert_eq!(first.0.body_spins_kg_m2_s(), retry.0.body_spins_kg_m2_s());
    assert_eq!(first.1, retry.1);
}
