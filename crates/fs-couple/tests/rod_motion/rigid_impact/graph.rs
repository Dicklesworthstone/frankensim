//! G0/G1/G3/G4: original momentum exchange across a mixed contact graph.

#[path = "graph/friction.rs"]
mod friction;

#[path = "graph/network.rs"]
mod network;

use super::*;
use fs_contact::multicontact::{MultiImpactConfig, MultiImpactDeterminacy};
use fs_couple::rod_motion::rigid_impact::graph::{
    RodRigidEndpoint, RodRigidGraphBudget, RodRigidGraphContact, RodRigidGraphImpactReport,
    RodRigidParticipant,
};
use fs_couple::rod_motion::rigid_impact::multiple::RodRigidMultiImpactError;

fn graph_config() -> MultiImpactConfig {
    let mut controls = MultiImpactConfig::default();
    controls.solver.gradient_tolerance = 1e-12;
    controls
}
fn rod_point(node: usize, point_body_m: [f64; 3]) -> RodRigidEndpoint {
    RodRigidEndpoint::Rod { node, point_body_m }
}
fn rigid_point(body: usize, arm_body_m: [f64; 3]) -> RodRigidEndpoint {
    RodRigidEndpoint::Rigid { body, arm_body_m }
}
fn edge(
    positive: RodRigidEndpoint,
    negative: RodRigidEndpoint,
    normal_world: [f64; 3],
    restitution: f64,
) -> RodRigidGraphContact {
    RodRigidGraphContact {
        positive,
        negative,
        normal_world,
        activation_tolerance_m: 1e-8,
        restitution,
    }
}
fn participant(mass: f64, position: [f64; 3], velocity: [f64; 3]) -> RodRigidParticipant {
    let inertia = inertia(mass);
    RodRigidParticipant {
        state: body(inertia, position, velocity),
        inertia,
    }
}
fn audit_graph(report: &RodRigidGraphImpactReport) {
    for row in &report.transfer.endpoint_velocity_residuals_m_s {
        near(row[0], 0.0, 2e-11);
    }
    if let Some(balance) = report.transfer.balance {
        vector_near(balance.linear_momentum_residual_n_s, [0.0; 3], 2e-11);
        vector_near(balance.angular_momentum_residual_nm_s, [0.0; 3], 2e-11);
        near(balance.energy_residual_j, 0.0, 2e-11);
        near(balance.work_energy_residual_j, 0.0, 2e-11);
        near(
            balance.relative_impulse_work_j,
            report.contact.relative_work,
            2e-11,
        );
    }
}

#[test]
fn singleton_and_reversed_orientation_match_original_pair_and_preserve_unused_bodies() {
    let gate = CancelGate::new();
    let rod = free_rod([-0.6, 0.0, 0.0]);
    let bodies = [
        participant(2.3, [0.8, 0.0, 0.0], [1.7, 0.0, 0.0]),
        participant(1.1, [7.0, -2.0, 3.0], [0.1, -0.2, 0.3]),
    ];
    for restitution in [0.0, 0.35, 1.0] {
        let contact = edge(
            rod_point(1, [0.0; 3]),
            rigid_point(0, [0.0; 3]),
            [1.0, 0.0, 0.0],
            restitution,
        );
        let scalar = rod
            .apply_rigid_normal_impact(
                bodies[0].state,
                bodies[0].inertia,
                witness(1, restitution),
                impact_controls(),
                &gate,
            )
            .unwrap();
        let (next, finite, report) = rod
            .apply_rigid_normal_contact_graph(
                &bodies,
                &[contact],
                impact_controls(),
                graph_config(),
                RodRigidGraphBudget::default(),
                &gate,
            )
            .unwrap();
        vector_near(
            next.linear_momenta_kg_m_s()[1],
            scalar.0.linear_momenta_kg_m_s()[1],
            2e-12,
        );
        vector_near(
            array(finite[0].linear_momentum_world()),
            array(scalar.1.linear_momentum_world()),
            2e-12,
        );
        near(
            report.contact.normal_impulses[0],
            scalar.2.contact.normal_impulse,
            2e-12,
        );
        assert_eq!(finite[1], bodies[1].state);
        assert!(report.transfer.rigid_impulses[1].is_none());
        assert_eq!(next.pose(), rod.pose());
        assert_eq!(next.time_s(), rod.time_s());
        assert_eq!(next.accepted_steps(), rod.accepted_steps());
        assert!(std::ptr::eq(next.model(), rod.model()));
        let reversed = edge(
            contact.negative,
            contact.positive,
            [-1.0, 0.0, 0.0],
            restitution,
        );
        let reverse = rod
            .apply_rigid_normal_contact_graph(
                &bodies,
                &[reversed],
                impact_controls(),
                graph_config(),
                RodRigidGraphBudget::default(),
                &gate,
            )
            .unwrap();
        vector_near(
            reverse.0.linear_momenta_kg_m_s()[1],
            next.linear_momenta_kg_m_s()[1],
            2e-12,
        );
        vector_near(
            array(reverse.1[0].linear_momentum_world()),
            array(finite[0].linear_momentum_world()),
            2e-12,
        );
        vector_near(
            reverse.2.transfer.impulses_world_n_s[0],
            report.transfer.impulses_world_n_s[0].map(|v| -v),
            2e-12,
        );
        audit_graph(&report);
        let continued = next
            .step_candidate(1e-5, &[RodNodalWrench::default(); 2], controls(1e-5), &gate)
            .unwrap()
            .0;
        assert_eq!(continued.accepted_steps(), rod.accepted_steps() + 1);
        let body_step = fs_mbd::RigidBodyIntegrator::new(fs_mbd::Gravity::ZERO)
            .step(finite[0], bodies[0].inertia, fs_mbd::Wrench::ZERO, 1e-5)
            .unwrap();
        assert_eq!(
            body_step.state_after.linear_momentum_world(),
            finite[0].linear_momentum_world()
        );
    }
}

#[test]
fn mixed_chain_couples_rod_self_contact_rod_body_and_body_body_without_pair_order() {
    let gate = CancelGate::new();
    let base = free_rod([0.0; 3]);
    let model = Arc::new(base.model().clone());
    let rod = RodMotion::new(
        model.clone(),
        model.initial_pose(),
        &[[3.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
        &[[0.0; 3]; 2],
        0.0,
        &gate,
    )
    .unwrap();
    let bodies = [
        participant(2.3, [0.4, 0.0, 0.0], [1.0, 0.0, 0.0]),
        participant(1.4, [0.4, 0.0, 0.0], [0.0; 3]),
    ];
    let masses = [
        model.inertia()[0].mass_kg(),
        model.inertia()[1].mass_kg(),
        2.3,
        1.4,
    ];
    let velocities = [3.0, 2.0, 1.0, 0.0];
    let mean = masses
        .iter()
        .zip(velocities)
        .map(|(m, v)| m * v)
        .sum::<f64>()
        / masses.iter().sum::<f64>();
    for restitution in [0.0, 0.4, 1.0] {
        let contacts = [
            edge(
                rod_point(1, [-0.4, 0.0, 0.0]),
                rod_point(0, [0.4, 0.0, 0.0]),
                [1.0, 0.0, 0.0],
                restitution,
            ),
            edge(
                rigid_point(0, [0.0; 3]),
                rod_point(1, [-0.4, 0.0, 0.0]),
                [1.0, 0.0, 0.0],
                restitution,
            ),
            edge(
                rigid_point(1, [0.0; 3]),
                rigid_point(0, [0.0; 3]),
                [1.0, 0.0, 0.0],
                restitution,
            ),
        ];
        let (next, finite, report) = rod
            .apply_rigid_normal_contact_graph(
                &bodies,
                &contacts,
                impact_controls(),
                graph_config(),
                RodRigidGraphBudget::default(),
                &gate,
            )
            .unwrap();
        let expected: [f64; 4] =
            std::array::from_fn(|i| (1.0 + restitution) * mean - restitution * velocities[i]);
        near(next.velocities_m_s().unwrap()[0][0], expected[0], 2e-11);
        near(next.velocities_m_s().unwrap()[1][0], expected[1], 2e-11);
        for i in 0..2 {
            near(
                finite[i]
                    .center_of_mass_velocity_world(bodies[i].inertia)
                    .unwrap()
                    .x,
                expected[i + 2],
                2e-11,
            );
            vector_near(array(finite[i].angular_momentum_body()), [0.0; 3], 2e-12);
        }
        let mut prefix = 0.0;
        for i in 0..3 {
            prefix += masses[i] * (velocities[i] - expected[i]);
            near(report.contact.normal_impulses[i], prefix, 2e-11);
            for j in 0..3 {
                let value = if i == j {
                    1.0 / masses[i] + 1.0 / masses[i + 1]
                } else if i + 1 == j {
                    -1.0 / masses[j]
                } else if j + 1 == i {
                    -1.0 / masses[i]
                } else {
                    0.0
                };
                near(report.delassus_kg_inverse[3 * i + j], value, 2e-12);
            }
        }
        let variance: f64 = masses
            .iter()
            .zip(velocities)
            .map(|(m, v)| 0.5 * m * (v - mean) * (v - mean))
            .sum();
        near(
            report.contact.dissipated_energy,
            (1.0 - restitution * restitution) * variance,
            2e-11,
        );
        let reindex = |p| match p {
            RodRigidEndpoint::Rigid { body, arm_body_m } => rigid_point(1 - body, arm_body_m),
            other => other,
        };
        let permuted: Vec<_> = [2, 0, 1]
            .into_iter()
            .map(|i| RodRigidGraphContact {
                positive: reindex(contacts[i].positive),
                negative: reindex(contacts[i].negative),
                ..contacts[i]
            })
            .collect();
        let reordered = rod
            .apply_rigid_normal_contact_graph(
                &[bodies[1], bodies[0]],
                &permuted,
                impact_controls(),
                graph_config(),
                RodRigidGraphBudget::default(),
                &gate,
            )
            .unwrap();
        for i in 0..2 {
            vector_near(
                reordered.0.linear_momenta_kg_m_s()[i],
                next.linear_momenta_kg_m_s()[i],
                2e-11,
            );
            vector_near(
                array(reordered.1[1 - i].linear_momentum_world()),
                array(finite[i].linear_momentum_world()),
                2e-11,
            );
        }
        audit_graph(&report);
        let first = rod
            .apply_rigid_normal_contact_graph(
                &bodies,
                &contacts[..1],
                impact_controls(),
                graph_config(),
                RodRigidGraphBudget::default(),
                &gate,
            )
            .unwrap();
        let body_states: Vec<_> = bodies
            .iter()
            .zip(first.1)
            .map(|(body, state)| RodRigidParticipant { state, ..*body })
            .collect();
        let sequential = first
            .0
            .apply_rigid_normal_contact_graph(
                &body_states,
                &contacts[1..],
                impact_controls(),
                graph_config(),
                RodRigidGraphBudget::default(),
                &gate,
            )
            .unwrap();
        assert!((sequential.0.velocities_m_s().unwrap()[0][0] - expected[0]).abs() > 1e-3);
    }
}

#[test]
fn two_finite_bodies_on_one_eccentric_node_match_independent_spatial_mass_and_full_tensors() {
    let gate = CancelGate::new();
    let model = Arc::new(
        eccentric_source()
            .into_spatial_motion_model(RodRootBoundary::Free, &gate)
            .unwrap(),
    );
    let k = model.inertia()[1];
    let center = k.center_of_mass_body_m();
    let tensor = k.inertia_com_body_kg_m2();
    let spatial = fs_mbd::articulated::SpatialInertia::new(
        k.mass_kg(),
        Vec3::new(center[0], center[1], center[2]),
        Mat3 {
            m: std::array::from_fn(|i| tensor[i / 3][i % 3]),
        },
    )
    .unwrap()
    .matrix()
    .unwrap();
    let factor = fs_la::factor::lu(&spatial.m, 6).unwrap();
    let rod_arms = [[0.06, -0.03, 0.04], [-0.02, 0.05, 0.03]];
    let body_arms = [[0.03, 0.02, -0.01], [-0.02, 0.04, 0.05]];
    let normals = [[1.0, 0.0, 0.0], [0.6, 0.8, 0.0]];
    let body_tensor = Mat3 {
        m: [0.07, 0.01, -0.004, 0.01, 0.09, 0.006, -0.004, 0.006, 0.08],
    };
    let properties = [
        MassProperties::from_inertia_tensor(1.7, BodyVec::ZERO, body_tensor).unwrap(),
        MassProperties::from_inertia_tensor(2.6, BodyVec::ZERO, body_tensor).unwrap(),
    ];
    let body_factor = fs_la::factor::lu(&body_tensor.m, 3).unwrap();
    let mut rows = [[0.0; 6]; 2];
    let mut responses = [[0.0; 6]; 2];
    let mut w = [[0.0; 2]; 2];
    for i in 0..2 {
        let torque = cross(rod_arms[i], normals[i]);
        rows[i] = [
            torque[0],
            torque[1],
            torque[2],
            normals[i][0],
            normals[i][1],
            normals[i][2],
        ];
        responses[i] = rows[i];
        factor.solve(&mut responses[i]);
    }
    for i in 0..2 {
        for j in 0..2 {
            w[i][j] = (0..6).map(|c| rows[i][c] * responses[j][c]).sum();
        }
        let torque = cross(body_arms[i], normals[i]);
        let mut response = torque;
        body_factor.solve(&mut response);
        w[i][i] +=
            1.0 / properties[i].mass() + (0..3).map(|c| torque[c] * response[c]).sum::<f64>();
    }
    let expected_impulses = [0.17, 0.23];
    let approach: [f64; 2] =
        std::array::from_fn(|i| (0..2).map(|j| w[i][j] * expected_impulses[j]).sum());
    assert!(approach.iter().all(|v| *v > 0.0));
    let mut response: [f64; 6] =
        std::array::from_fn(|c| (0..2).map(|i| rows[i][c] * expected_impulses[i]).sum());
    factor.solve(&mut response);
    let q = So3::exp(So3Tangent::new(Vec3::new(0.2, -0.1, 0.3))).unwrap();
    for (rotation, translation, boost) in [
        (So3::identity(), [0.0; 3], [0.0; 3]),
        (q, [0.75, -0.25, 0.125], [0.4, -0.3, 0.2]),
    ] {
        let mut pose = model.initial_pose();
        for node in 0..2 {
            pose.positions_m[node] = add(rotate(rotation, pose.positions_m[node]), translation);
            pose.frames[node] = rotation;
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
        let mut bodies = Vec::new();
        let mut contacts = Vec::new();
        for i in 0..2 {
            let normal = rotate(rotation, normals[i]);
            let rod_point_world = add(rod.pose().positions_m[1], rotate(rotation, rod_arms[i]));
            let position = sub(
                sub(rod_point_world, rotate(rotation, body_arms[i])),
                normal.map(|x| x * 1e-9),
            );
            let velocity = add(boost, normal.map(|x| x * approach[i]));
            let state = RigidBodyState::new(
                Pose::new(vector(position), orientation(rotation)).unwrap(),
                vector(velocity).scale(properties[i].mass()),
                BodyVec::ZERO,
            )
            .unwrap();
            bodies.push(RodRigidParticipant {
                state,
                inertia: properties[i],
            });
            contacts.push(edge(
                rod_point(1, rod_arms[i]),
                rigid_point(i, body_arms[i]),
                normal,
                0.0,
            ));
        }
        let (next, finite, report) = rod
            .apply_rigid_normal_contact_graph(
                &bodies,
                &contacts,
                impact_controls(),
                graph_config(),
                RodRigidGraphBudget::default(),
                &gate,
            )
            .unwrap();
        for i in 0..2 {
            near(
                report.contact.normal_impulses[i],
                expected_impulses[i],
                2e-11,
            );
            for j in 0..2 {
                near(report.delassus_kg_inverse[2 * i + j], w[i][j], 2e-11);
            }
            let impulse_body = normals[i].map(|v| v * expected_impulses[i]);
            vector_near(
                array(finite[i].linear_momentum_world()),
                sub(
                    array(bodies[i].state.linear_momentum_world()),
                    rotate(rotation, impulse_body),
                ),
                2e-11,
            );
            vector_near(
                array(finite[i].angular_momentum_body()),
                cross(body_arms[i], impulse_body).map(|v| -v),
                2e-11,
            );
        }
        vector_near(
            next.angular_velocities_body_rad_s().unwrap()[1],
            [response[0], response[1], response[2]],
            2e-11,
        );
        vector_near(
            next.velocities_m_s().unwrap()[1],
            add(
                boost,
                rotate(rotation, [response[3], response[4], response[5]]),
            ),
            2e-11,
        );
        audit_graph(&report);
    }
}

#[test]
fn rod_only_and_rigid_only_redundancy_preserve_native_motion_and_unaffected_participants() {
    let gate = CancelGate::new();
    let base = free_rod([0.0; 3]);
    let model = Arc::new(base.model().clone());
    let rod = RodMotion::new(
        model.clone(),
        model.initial_pose(),
        &[[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]],
        &[[0.0; 3]; 2],
        0.0,
        &gate,
    )
    .unwrap();
    let contact = edge(
        rod_point(1, [-0.4, 0.0, 0.0]),
        rod_point(0, [0.4, 0.0, 0.0]),
        [1.0, 0.0, 0.0],
        0.4,
    );
    let (next, bodies, report) = rod
        .apply_rigid_normal_contact_graph(
            &[],
            &[contact, contact],
            impact_controls(),
            graph_config(),
            RodRigidGraphBudget { maximum_bodies: 0 },
            &gate,
        )
        .unwrap();
    near(next.velocities_m_s().unwrap()[0][0], -0.4, 2e-11);
    near(next.velocities_m_s().unwrap()[1][0], 0.4, 2e-11);
    assert!(bodies.is_empty());
    assert!(matches!(
        report.contact.determinacy,
        MultiImpactDeterminacy::UniqueVelocity(_)
    ));
    audit_graph(&report);

    let bodies = [
        participant(1.0, [0.4, 0.0, 0.0], [1.0, 0.0, 0.0]),
        participant(2.0, [0.4, 0.0, 0.0], [-1.0, 0.0, 0.0]),
        participant(3.0, [0.4, 0.0, 0.0], [0.0; 3]),
    ];
    let contacts = [
        edge(
            rigid_point(0, [0.0; 3]),
            rigid_point(1, [0.0; 3]),
            [1.0, 0.0, 0.0],
            0.0,
        ),
        edge(
            rigid_point(1, [0.0; 3]),
            rigid_point(2, [0.0; 3]),
            [1.0, 0.0, 0.0],
            0.0,
        ),
        edge(
            rigid_point(2, [0.0; 3]),
            rigid_point(0, [0.0; 3]),
            [1.0, 0.0, 0.0],
            0.0,
        ),
    ];
    let (untouched, finite, report) = rod
        .apply_rigid_normal_contact_graph(
            &bodies,
            &contacts,
            impact_controls(),
            graph_config(),
            RodRigidGraphBudget::default(),
            &gate,
        )
        .unwrap();
    assert!(report.relative_normal_speeds_before_m_s[0] > 0.0);
    assert!(matches!(
        report.contact.determinacy,
        MultiImpactDeterminacy::UniqueVelocity(_)
    ));
    assert!(report.transfer.rod_impulse.is_none());
    assert_eq!(
        untouched.linear_momenta_kg_m_s(),
        rod.linear_momenta_kg_m_s()
    );
    assert_eq!(untouched.body_spins_kg_m2_s(), rod.body_spins_kg_m2_s());
    for i in 0..3 {
        near(
            finite[i]
                .center_of_mass_velocity_world(bodies[i].inertia)
                .unwrap()
                .x,
            -1.0 / 6.0,
            2e-11,
        );
    }
    near(report.contact.dissipated_energy, 1.5 - 1.0 / 12.0, 2e-11);
    audit_graph(&report);
}

#[test]
fn original_stationary_supports_receive_the_complete_multi_body_force_and_couple() {
    let gate = CancelGate::new();
    let source = eccentric_source()
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
            .into_spatial_motion_model(RodRootBoundary::Free, &gate)
            .unwrap(),
    );
    let rod = SpatialRodMotion::new_spatial(
        model.clone(),
        model.initial_pose(),
        &[[0.0; 3]; 2],
        &[[0.0; 3]; 2],
        0.0,
        1e-12,
        &gate,
    )
    .unwrap();
    let arms = [[0.0, 0.07, 0.0], [0.05, 0.0, 0.03]];
    let normals = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let bodies = [
        participant(2.3, add([0.8, 0.0, 0.0], arms[0]), [1.0, 0.0, 0.0]),
        participant(1.4, add([0.8, 0.0, 0.0], arms[1]), [0.0, 2.0, 0.0]),
    ];
    let contacts: Vec<_> = (0..2)
        .map(|i| {
            edge(
                rod_point(1, arms[i]),
                rigid_point(i, [0.0; 3]),
                normals[i],
                0.0,
            )
        })
        .collect();
    let (next, finite, report) = rod
        .apply_rigid_normal_contact_graph(
            &bodies,
            &contacts,
            impact_controls(),
            graph_config(),
            RodRigidGraphBudget::default(),
            &gate,
        )
        .unwrap();
    let expected = [[2.3, 0.0, 0.0], [0.0, 2.8, 0.0]];
    let support = &report
        .transfer
        .rod_impulse
        .as_ref()
        .unwrap()
        .support_impulses[1];
    vector_near(
        support.linear_n_s,
        add(expected[0], expected[1]).map(|v| -v),
        2e-11,
    );
    vector_near(
        support.angular_nm_s,
        add(cross(arms[0], expected[0]), cross(arms[1], expected[1])).map(|v| -v),
        2e-11,
    );
    assert_eq!(next.linear_momenta_kg_m_s(), rod.linear_momenta_kg_m_s());
    assert_eq!(next.body_spins_kg_m2_s(), rod.body_spins_kg_m2_s());
    for body in finite {
        vector_near(array(body.linear_momentum_world()), [0.0; 3], 2e-11);
    }
    audit_graph(&report);
}

#[test]
fn graph_limits_geometry_noops_and_per_body_resolution_refuse_atomically() {
    let gate = CancelGate::new();
    let rod = free_rod([0.0; 3]);
    let bodies = [participant(2.3, [0.8, 0.0, 0.0], [1.0, 0.0, 0.0])];
    let contact = edge(
        rod_point(1, [0.0; 3]),
        rigid_point(0, [0.0; 3]),
        [1.0, 0.0, 0.0],
        0.0,
    );
    let empty = rod
        .apply_rigid_normal_contact_graph(
            &bodies,
            &[],
            impact_controls(),
            graph_config(),
            RodRigidGraphBudget::default(),
            &gate,
        )
        .unwrap();
    assert_eq!(empty.1[0], bodies[0].state);
    assert_eq!(empty.0.linear_momenta_kg_m_s(), rod.linear_momenta_kg_m_s());
    assert!(empty.2.transfer.rod_impulse.is_none() && empty.2.transfer.balance.is_none());
    assert!(
        rod.apply_rigid_normal_contact_graph(
            &bodies,
            &[contact],
            impact_controls(),
            graph_config(),
            RodRigidGraphBudget { maximum_bodies: 0 },
            &gate
        )
        .is_err()
    );
    assert!(
        rod.apply_rigid_normal_contact_graph(
            &[bodies[0]; 129],
            &[contact],
            impact_controls(),
            graph_config(),
            RodRigidGraphBudget::default(),
            &gate
        )
        .is_err()
    );
    assert!(
        rod.apply_rigid_normal_contact_graph(
            &bodies,
            &[contact; 65],
            impact_controls(),
            graph_config(),
            RodRigidGraphBudget::default(),
            &gate
        )
        .is_err()
    );
    for malformed in [
        RodRigidGraphContact {
            positive: rigid_point(7, [0.0; 3]),
            ..contact
        },
        RodRigidGraphContact {
            negative: rod_point(1, [0.1, 0.0, 0.0]),
            ..contact
        },
        RodRigidGraphContact {
            positive: rigid_point(0, [0.1, 0.0, 0.0]),
            ..contact
        },
        RodRigidGraphContact {
            normal_world: [0.0; 3],
            ..contact
        },
    ] {
        assert!(
            rod.apply_rigid_normal_contact_graph(
                &bodies,
                &[malformed],
                impact_controls(),
                graph_config(),
                RodRigidGraphBudget::default(),
                &gate
            )
            .is_err()
        );
    }
    let penetration = [participant(2.3, [0.81, 0.0, 0.0], [1.0, 0.0, 0.0])];
    assert!(
        rod.apply_rigid_normal_contact_graph(
            &penetration,
            &[contact],
            impact_controls(),
            graph_config(),
            RodRigidGraphBudget::default(),
            &gate
        )
        .is_err()
    );
    let mut broad = impact_controls();
    broad.geometry_tolerance_m = 1e-3;
    let crooked = [participant(2.3, [0.8, 1e-4, 0.0], [1.0, 0.0, 0.0])];
    assert!(matches!(
        rod.apply_rigid_normal_contact_graph(
            &crooked,
            &[contact],
            broad,
            graph_config(),
            RodRigidGraphBudget::default(),
            &gate
        ),
        Err(RodRigidMultiImpactError::Native(
            RodRigidImpactError::Balance { .. }
        ))
    ));
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        rod.apply_rigid_normal_contact_graph(
            &bodies,
            &[contact],
            impact_controls(),
            graph_config(),
            RodRigidGraphBudget::default(),
            &cancelled
        )
        .is_err()
    );
    let first = rod
        .apply_rigid_normal_contact_graph(
            &bodies,
            &[contact],
            impact_controls(),
            graph_config(),
            RodRigidGraphBudget::default(),
            &gate,
        )
        .unwrap();
    let retry = rod
        .apply_rigid_normal_contact_graph(
            &bodies,
            &[contact],
            impact_controls(),
            graph_config(),
            RodRigidGraphBudget::default(),
            &gate,
        )
        .unwrap();
    assert_eq!(
        first.0.linear_momenta_kg_m_s(),
        retry.0.linear_momenta_kg_m_s()
    );
    assert_eq!(first.1, retry.1);

    // The native body owner records commanded versus represented impulses.
    // Opposite unresolved transfers cancel globally, but each original body
    // must still meet the declared absolute linear-transfer bound.
    let mass = inertia(2.3);
    let rotating: Vec<_> = [1.0, -1.0]
        .into_iter()
        .map(|spin| RodRigidParticipant {
            state: RigidBodyState::new(
                Pose::new(vector([0.8, 0.0, 0.0]), UnitQuaternion::IDENTITY).unwrap(),
                BodyVec::new(2.3e20, 0.0, 0.0),
                mass.angular_momentum_body_checked(BodyVec::new(0.0, 0.0, spin))
                    .unwrap(),
            )
            .unwrap(),
            inertia: mass,
        })
        .collect();
    let rotating_contact = edge(
        rigid_point(0, [0.0, 0.1, 0.0]),
        rigid_point(1, [0.0, 0.1, 0.0]),
        [1.0, 0.0, 0.0],
        0.0,
    );
    let mut lenient_work = impact_controls();
    lenient_work.velocity_tolerance_m_s = 1.0;
    lenient_work.impulse.work_tolerance_j = 1.0;
    match rod.apply_rigid_normal_contact_graph(
        &rotating,
        &[rotating_contact],
        lenient_work,
        graph_config(),
        RodRigidGraphBudget::default(),
        &gate,
    ) {
        Err(RodRigidMultiImpactError::Native(RodRigidImpactError::Balance {
            linear_defect_n_s,
            ..
        })) => assert!(linear_defect_n_s > 0.01),
        outcome => panic!(
            "per-body unresolved transfer was not refused: {:?}",
            outcome.map(|x| x.2)
        ),
    }
}
