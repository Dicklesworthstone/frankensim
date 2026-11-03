//! G0/G1/G3/G4: a coupled event across independently continuing native rods.

#[path = "network/friction.rs"]
mod friction;

use super::*;
use fs_couple::rod_motion::rigid_impact::graph::network::{
    RodRigidNetworkBudget, RodRigidNetworkContact, RodRigidNetworkEndpoint,
    RodRigidNetworkImpactReport, apply_rigid_normal_contact_network,
};

fn tip_rod(length: f64, tip: [f64; 3], velocity: [f64; 3], time: f64) -> RodMotion {
    let gate = CancelGate::new();
    let model = Arc::new(
        SpatialMaterialRod::try_straight(length, vec![section(1.0)])
            .unwrap()
            .into_motion_model(RodRootBoundary::Free, &gate)
            .unwrap(),
    );
    let mut pose = model.initial_pose();
    let shift = sub(tip, [length, 0.0, 0.0]);
    for point in &mut pose.positions_m {
        *point = add(*point, shift);
    }
    RodMotion::new(model, pose, &[velocity; 2], &[[0.0; 3]; 2], time, &gate).unwrap()
}

fn rod_tip(rod: usize, point_body_m: [f64; 3]) -> RodRigidNetworkEndpoint {
    RodRigidNetworkEndpoint::Rod {
        rod,
        node: 1,
        point_body_m,
    }
}
fn finite_point(body: usize, arm_body_m: [f64; 3]) -> RodRigidNetworkEndpoint {
    RodRigidNetworkEndpoint::Rigid { body, arm_body_m }
}
fn network_edge(
    positive: RodRigidNetworkEndpoint,
    negative: RodRigidNetworkEndpoint,
    normal_world: [f64; 3],
    restitution: f64,
) -> RodRigidNetworkContact {
    RodRigidNetworkContact {
        positive,
        negative,
        normal_world,
        activation_tolerance_m: 1e-8,
        restitution,
    }
}
fn audit_network(report: &RodRigidNetworkImpactReport) {
    for row in &report.transfer.endpoint_velocity_residuals_m_s {
        near(row[0], 0.0, 2e-11);
    }
    if let Some(balance) = report.transfer.balance {
        vector_near(balance.linear_momentum_residual_n_s, [0.0; 3], 2e-11);
        vector_near(balance.angular_momentum_residual_nm_s, [0.0; 3], 2e-11);
        near(balance.work_energy_residual_j, 0.0, 2e-11);
        near(balance.energy_residual_j, 0.0, 2e-11);
    }
}

#[test]
fn independent_rods_and_body_match_simultaneous_chain_and_preserve_native_owners() {
    let gate = CancelGate::new();
    let point = [0.8, 0.0, 0.0];
    let rods = [
        tip_rod(0.8, point, [3.0, 0.0, 0.0], 2.0),
        tip_rod(1.6, point, [2.0, 0.0, 0.0], 2.0),
        tip_rod(2.4, point, [1.0, 0.0, 0.0], 2.0),
        tip_rod(0.5, [7.0, 8.0, 9.0], [0.1, 0.2, 0.0], 19.0),
    ];
    let bodies = [
        participant(3.7, point, [0.0; 3]),
        participant(1.2, [5.0, 6.0, 7.0], [0.5, 0.2, -0.3]),
    ];
    let masses = [
        rods[0].model().inertia()[1].mass_kg(),
        rods[1].model().inertia()[1].mass_kg(),
        rods[2].model().inertia()[1].mass_kg(),
        bodies[0].inertia.mass(),
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
            network_edge(
                rod_tip(1, [0.0; 3]),
                rod_tip(0, [0.0; 3]),
                [1.0, 0.0, 0.0],
                restitution,
            ),
            network_edge(
                rod_tip(2, [0.0; 3]),
                rod_tip(1, [0.0; 3]),
                [1.0, 0.0, 0.0],
                restitution,
            ),
            network_edge(
                finite_point(0, [0.0; 3]),
                rod_tip(2, [0.0; 3]),
                [1.0, 0.0, 0.0],
                restitution,
            ),
        ];
        let (next, finite, report) = apply_rigid_normal_contact_network(
            &rods,
            &bodies,
            &contacts,
            impact_controls(),
            graph_config(),
            RodRigidNetworkBudget::default(),
            &gate,
        )
        .unwrap();
        let expected_velocity = velocities.map(|v| mean - restitution * (v - mean));
        for i in 0..3 {
            let impulse = (1.0 + restitution)
                * (0..=i)
                    .map(|j| masses[j] * (velocities[j] - mean))
                    .sum::<f64>();
            near(report.contact.normal_impulses[i], impulse, 2e-11);
            near(
                next[i].velocities_m_s().unwrap()[1][0],
                expected_velocity[i],
                2e-11,
            );
            assert!(report.transfer.rod_impulses[i].is_some());
            assert!(std::ptr::eq(next[i].model(), rods[i].model()));
            assert_eq!(next[i].pose(), rods[i].pose());
            assert_eq!(next[i].time_s(), rods[i].time_s());
            assert_eq!(next[i].accepted_steps(), rods[i].accepted_steps());
            assert_eq!(
                next[i].linear_momenta_kg_m_s()[0],
                rods[i].linear_momenta_kg_m_s()[0]
            );
            let continued = next[i]
                .step_candidate(1e-6, &[RodNodalWrench::default(); 2], controls(1e-6), &gate)
                .unwrap()
                .0;
            assert_eq!(continued.accepted_steps(), rods[i].accepted_steps() + 1);
        }
        near(
            finite[0]
                .center_of_mass_velocity_world(bodies[0].inertia)
                .unwrap()
                .x,
            expected_velocity[3],
            2e-11,
        );
        assert_eq!(finite[1], bodies[1].state);
        assert!(report.transfer.rigid_impulses[1].is_none());
        assert!(report.transfer.rod_impulses[3].is_none());
        assert_eq!(next[3].pose(), rods[3].pose());
        assert_eq!(next[3].time_s(), 19.0);
        assert_eq!(
            next[3].linear_momenta_kg_m_s(),
            rods[3].linear_momenta_kg_m_s()
        );
        let variance = masses
            .iter()
            .zip(velocities)
            .map(|(m, v)| m * (v - mean) * (v - mean))
            .sum::<f64>();
        near(
            report.contact.dissipated_energy,
            0.5 * (1.0 - restitution * restitution) * variance,
            2e-11,
        );
        audit_network(&report);

        // Neither the rod index nor the contact order defines a pair sequence.
        let order = [2, 0, 3, 1];
        let inverse = [1, 3, 0, 2];
        let remap = |endpoint| match endpoint {
            RodRigidNetworkEndpoint::Rod {
                rod,
                node,
                point_body_m,
            } => RodRigidNetworkEndpoint::Rod {
                rod: inverse[rod],
                node,
                point_body_m,
            },
            RodRigidNetworkEndpoint::Rigid { body, arm_body_m } => RodRigidNetworkEndpoint::Rigid {
                body: 1 - body,
                arm_body_m,
            },
        };
        let reordered: Vec<_> = order.into_iter().map(|i| rods[i].clone()).collect();
        let reversed: Vec<_> = contacts
            .iter()
            .rev()
            .map(|c| RodRigidNetworkContact {
                positive: remap(c.positive),
                negative: remap(c.negative),
                ..*c
            })
            .collect();
        let permuted = apply_rigid_normal_contact_network(
            &reordered,
            &[bodies[1], bodies[0]],
            &reversed,
            impact_controls(),
            graph_config(),
            RodRigidNetworkBudget::default(),
            &gate,
        )
        .unwrap();
        for i in 0..3 {
            vector_near(
                permuted.0[inverse[i]].linear_momenta_kg_m_s()[1],
                next[i].linear_momenta_kg_m_s()[1],
                2e-11,
            );
            near(
                permuted.2.contact.normal_impulses[2 - i],
                report.contact.normal_impulses[i],
                2e-11,
            );
        }
        vector_near(
            array(permuted.1[1].linear_momentum_world()),
            array(finite[0].linear_momentum_world()),
            2e-11,
        );

        if restitution == 0.0 {
            let mut sequential_rods = rods.to_vec();
            let mut sequential_bodies = bodies.to_vec();
            for contact in contacts {
                let (next, states, _) = apply_rigid_normal_contact_network(
                    &sequential_rods,
                    &sequential_bodies,
                    &[contact],
                    impact_controls(),
                    graph_config(),
                    RodRigidNetworkBudget::default(),
                    &gate,
                )
                .unwrap();
                sequential_rods = next;
                for (body, state) in sequential_bodies.iter_mut().zip(states) {
                    body.state = state;
                }
            }
            assert!((sequential_rods[0].velocities_m_s().unwrap()[1][0] - mean).abs() > 0.05);
        }
    }

    // Equal native masses carry exactly representable, almost equal momenta.
    // Their collision loses kick^2/m even though subtracting whole-trajectory
    // kinetic energies loses the event beneath the large common translation.
    let mass = masses[0];
    let common_momentum = 65_536.0;
    let kick = 1.0 / 4096.0;
    let moving = [
        tip_rod(0.8, point, [(common_momentum + kick) / mass, 0.0, 0.0], 2.0),
        tip_rod(0.8, point, [(common_momentum - kick) / mass, 0.0, 0.0], 2.0),
    ];
    assert_eq!(
        moving[0].linear_momenta_kg_m_s()[1][0],
        common_momentum + kick
    );
    assert_eq!(
        moving[1].linear_momenta_kg_m_s()[1][0],
        common_momentum - kick
    );
    let collision = network_edge(
        rod_tip(1, [0.0; 3]),
        rod_tip(0, [0.0; 3]),
        [1.0, 0.0, 0.0],
        0.0,
    );
    let (stopped, _, report) = apply_rigid_normal_contact_network(
        &moving,
        &[],
        &[collision],
        impact_controls(),
        graph_config(),
        RodRigidNetworkBudget::default(),
        &gate,
    )
    .unwrap();
    let loss = kick * kick / mass;
    near(report.contact.normal_impulses[0], kick, 2e-14);
    near(report.contact.dissipated_energy, loss, 2e-14);
    near(
        report.transfer.balance.unwrap().kinetic_energy_change_j,
        -loss,
        2e-14,
    );
    let naive = stopped
        .iter()
        .map(|rod| rod.observe([0.0; 3]).unwrap().kinetic_energy_j)
        .sum::<f64>()
        - moving
            .iter()
            .map(|rod| rod.observe([0.0; 3]).unwrap().kinetic_energy_j)
            .sum::<f64>();
    assert!((naive + loss).abs() > 0.5 * loss);
    audit_network(&report);
}

#[test]
fn two_eccentric_rods_and_tensor_body_match_independent_spatial_mass_blocks() {
    let gate = CancelGate::new();
    let model = Arc::new(
        eccentric_source()
            .into_spatial_motion_model(RodRootBoundary::Free, &gate)
            .unwrap(),
    );
    let k = model.inertia()[1];
    let tensor = k.inertia_com_body_kg_m2();
    let center = k.center_of_mass_body_m();
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
    let rod_factor = fs_la::factor::lu(&spatial.m, 6).unwrap();
    let body_tensor = Mat3 {
        m: [0.07, 0.01, -0.004, 0.01, 0.09, 0.006, -0.004, 0.006, 0.08],
    };
    let properties = MassProperties::from_inertia_tensor(1.7, BodyVec::ZERO, body_tensor).unwrap();
    let body_factor = fs_la::factor::lu(&body_tensor.m, 3).unwrap();
    let arms = [
        [0.06, -0.03, 0.04],
        [-0.02, 0.05, 0.03],
        [0.01, -0.03, 0.06],
        [0.03, 0.02, -0.01],
    ];
    let normals = [[1.0, 0.0, 0.0], [0.6, 0.8, 0.0]];
    let row = |arm, axis| {
        let torque = cross(arm, axis);
        [torque[0], torque[1], torque[2], axis[0], axis[1], axis[2]]
    };
    // Participant blocks are [rod0, rod1, finite body]. Rod1 is the shared
    // negative endpoint, with two DIFFERENT material-point lever arms.
    let rows = [
        [row(arms[0], normals[0]), [0.0; 6]],
        [
            row(arms[1], normals[0].map(|v| -v)),
            row(arms[2], normals[1].map(|v| -v)),
        ],
        [[0.0; 6], row(arms[3], normals[1])],
    ];
    let mut responses = rows;
    for p in 0..3 {
        for i in 0..2 {
            if p < 2 {
                rod_factor.solve(&mut responses[p][i]);
            } else {
                let mut angular = [responses[p][i][0], responses[p][i][1], responses[p][i][2]];
                body_factor.solve(&mut angular);
                responses[p][i][..3].copy_from_slice(&angular);
                for c in 3..6 {
                    responses[p][i][c] /= properties.mass();
                }
            }
        }
    }
    let w: [[f64; 2]; 2] = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            (0..3)
                .map(|p| {
                    (0..6)
                        .map(|c| rows[p][i][c] * responses[p][j][c])
                        .sum::<f64>()
                })
                .sum()
        })
    });
    let impulses = [0.17, 0.23];
    let approach: [f64; 2] = std::array::from_fn(|i| (0..2).map(|j| w[i][j] * impulses[j]).sum());
    assert!(approach.iter().all(|v| *v > 0.0));
    let increments: [[f64; 6]; 3] = std::array::from_fn(|p| {
        std::array::from_fn(|c| (0..2).map(|i| responses[p][i][c] * impulses[i]).sum())
    });
    let q = So3::exp(So3Tangent::new(Vec3::new(0.2, -0.1, 0.3))).unwrap();
    for (rotation, translation, boost) in [
        (So3::identity(), [0.0; 3], [0.0; 3]),
        (q, [0.75, -0.25, 0.125], [0.4, -0.3, 0.2]),
    ] {
        let shared_tip = [0.8, 0.0, 0.0];
        let rod0_shift = add(sub(arms[1], arms[0]), normals[0].map(|v| 1e-9 * v));
        let base_velocities = [normals[0].map(|v| -approach[0] * v), [0.0; 3]];
        let rods: Vec<_> = (0..2)
            .map(|index| {
                let mut pose = model.initial_pose();
                for node in 0..2 {
                    let shift = if index == 0 { rod0_shift } else { [0.0; 3] };
                    pose.positions_m[node] = add(
                        translation,
                        rotate(rotation, add(pose.positions_m[node], shift)),
                    );
                    pose.frames[node] = rotation;
                }
                SpatialRodMotion::new_spatial(
                    model.clone(),
                    pose,
                    &[add(boost, rotate(rotation, base_velocities[index])); 2],
                    &[[0.0; 3]; 2],
                    0.0,
                    1e-12,
                    &gate,
                )
                .unwrap()
            })
            .collect();
        let body_velocity = normals[1].map(|v| -approach[1] * v);
        let state = RigidBodyState::new(
            Pose::new(
                vector(add(
                    translation,
                    rotate(
                        rotation,
                        add(
                            sub(add(shared_tip, arms[2]), arms[3]),
                            normals[1].map(|v| 1e-9 * v),
                        ),
                    ),
                )),
                orientation(rotation),
            )
            .unwrap(),
            vector(add(boost, rotate(rotation, body_velocity))).scale(properties.mass()),
            BodyVec::ZERO,
        )
        .unwrap();
        let contacts = [
            network_edge(
                rod_tip(0, arms[0]),
                rod_tip(1, arms[1]),
                rotate(rotation, normals[0]),
                0.0,
            ),
            network_edge(
                finite_point(0, arms[3]),
                rod_tip(1, arms[2]),
                rotate(rotation, normals[1]),
                0.0,
            ),
        ];
        let (next, finite, report) = apply_rigid_normal_contact_network(
            &rods,
            &[RodRigidParticipant {
                state,
                inertia: properties,
            }],
            &contacts,
            impact_controls(),
            graph_config(),
            RodRigidNetworkBudget::default(),
            &gate,
        )
        .unwrap();
        for i in 0..2 {
            near(report.contact.normal_impulses[i], impulses[i], 2e-11);
            for j in 0..2 {
                near(report.delassus_kg_inverse[2 * i + j], w[i][j], 2e-11);
            }
            vector_near(
                next[i].angular_velocities_body_rad_s().unwrap()[1],
                [increments[i][0], increments[i][1], increments[i][2]],
                2e-11,
            );
            vector_near(
                next[i].velocities_m_s().unwrap()[1],
                add(
                    boost,
                    rotate(
                        rotation,
                        add(
                            base_velocities[i],
                            [increments[i][3], increments[i][4], increments[i][5]],
                        ),
                    ),
                ),
                2e-11,
            );
        }
        vector_near(
            array(
                properties
                    .angular_velocity_body_checked(finite[0].angular_momentum_body())
                    .unwrap(),
            ),
            [increments[2][0], increments[2][1], increments[2][2]],
            2e-11,
        );
        vector_near(
            array(finite[0].center_of_mass_velocity_world(properties).unwrap()),
            add(
                boost,
                rotate(
                    rotation,
                    add(
                        body_velocity,
                        [increments[2][3], increments[2][4], increments[2][5]],
                    ),
                ),
            ),
            2e-11,
        );
        audit_network(&report);
    }
}

#[test]
fn multiple_supported_rods_close_each_original_reaction_and_shared_body_couple() {
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
    let rods = [rod.clone(), rod];
    let arms = [[0.0, 0.07, 0.0], [0.05, 0.0, 0.03]];
    let normals = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let finite = participant(2.3, [0.8, 0.0, 0.0], [1.0, 2.0, 0.0]);
    let tensor = finite.inertia.inertia_tensor_body();
    let factor = fs_la::factor::lu(&tensor.m, 3).unwrap();
    let torques = [cross(arms[0], normals[0]), cross(arms[1], normals[1])];
    let mut responses = torques;
    for response in &mut responses {
        factor.solve(response);
    }
    let mut matrix = [0.0; 4];
    for i in 0..2 {
        for j in 0..2 {
            matrix[2 * i + j] = if i == j {
                1.0 / finite.inertia.mass()
            } else {
                0.0
            } + (0..3).map(|c| torques[i][c] * responses[j][c]).sum::<f64>();
        }
    }
    let mut impulses = [1.0, 2.0];
    fs_la::factor::lu(&matrix, 2).unwrap().solve(&mut impulses);
    let contacts: Vec<_> = (0..2)
        .map(|i| {
            network_edge(
                rod_tip(i, arms[i]),
                finite_point(0, arms[i]),
                normals[i],
                0.0,
            )
        })
        .collect();
    let (next, states, report) = apply_rigid_normal_contact_network(
        &rods,
        &[finite],
        &contacts,
        impact_controls(),
        graph_config(),
        RodRigidNetworkBudget::default(),
        &gate,
    )
    .unwrap();
    let mut force = [0.0; 3];
    let mut couple = [0.0; 3];
    for i in 0..2 {
        near(report.contact.normal_impulses[i], impulses[i], 2e-11);
        let applied = normals[i].map(|v| v * impulses[i]);
        let moment = cross(arms[i], applied);
        let support = &report.transfer.rod_impulses[i]
            .as_ref()
            .unwrap()
            .support_impulses[1];
        vector_near(support.linear_n_s, applied.map(|v| -v), 2e-11);
        vector_near(support.angular_nm_s, moment.map(|v| -v), 2e-11);
        assert_eq!(
            next[i].linear_momenta_kg_m_s(),
            rods[i].linear_momenta_kg_m_s()
        );
        assert_eq!(next[i].body_spins_kg_m2_s(), rods[i].body_spins_kg_m2_s());
        force = add(force, applied);
        couple = add(couple, moment);
    }
    vector_near(
        array(states[0].linear_momentum_world()),
        sub(array(finite.state.linear_momentum_world()), force),
        2e-11,
    );
    vector_near(
        array(states[0].angular_momentum_body()),
        couple.map(|v| -v),
        2e-11,
    );
    audit_network(&report);
}

#[test]
fn network_limits_clocks_native_compatibility_and_redundancy_are_atomic() {
    let gate = CancelGate::new();
    let point = [0.8, 0.0, 0.0];
    let rods = [
        tip_rod(0.8, point, [2.0, 0.0, 0.0], 2.0),
        tip_rod(1.6, point, [0.0; 3], 2.0),
    ];
    let contact = network_edge(
        rod_tip(1, [0.0; 3]),
        rod_tip(0, [0.0; 3]),
        [1.0, 0.0, 0.0],
        0.0,
    );
    let call =
        |rods: &[RodMotion], contacts: &[RodRigidNetworkContact], budget, gate: &CancelGate| {
            apply_rigid_normal_contact_network(
                rods,
                &[],
                contacts,
                impact_controls(),
                graph_config(),
                budget,
                gate,
            )
        };
    let empty = call(&rods, &[], RodRigidNetworkBudget::default(), &gate).unwrap();
    assert!(empty.2.transfer.balance.is_none());
    assert!(empty.2.transfer.rod_impulses.iter().all(Option::is_none));
    for i in 0..2 {
        assert_eq!(
            empty.0[i].linear_momenta_kg_m_s(),
            rods[i].linear_momenta_kg_m_s()
        );
    }
    for budget in [
        RodRigidNetworkBudget {
            maximum_rods: 1,
            ..RodRigidNetworkBudget::default()
        },
        RodRigidNetworkBudget {
            maximum_rods: 129,
            ..RodRigidNetworkBudget::default()
        },
        RodRigidNetworkBudget {
            maximum_rod_nodes: 3,
            ..RodRigidNetworkBudget::default()
        },
    ] {
        assert!(call(&rods, &[contact], budget, &gate).is_err());
    }
    assert!(
        call(
            &vec![rods[0].clone(); 129],
            &[],
            RodRigidNetworkBudget::default(),
            &gate
        )
        .is_err()
    );
    assert!(
        call(
            &rods,
            &[contact; 65],
            RodRigidNetworkBudget::default(),
            &gate
        )
        .is_err()
    );
    for bad in [
        RodRigidNetworkContact {
            positive: rod_tip(9, [0.0; 3]),
            ..contact
        },
        RodRigidNetworkContact {
            positive: rod_tip(0, [0.0; 3]),
            ..contact
        },
        RodRigidNetworkContact {
            positive: RodRigidNetworkEndpoint::Rod {
                rod: 1,
                node: 9,
                point_body_m: [0.0; 3],
            },
            ..contact
        },
    ] {
        assert!(call(&rods, &[bad], RodRigidNetworkBudget::default(), &gate).is_err());
    }
    let different_time = [rods[0].clone(), tip_rod(1.6, point, [0.0; 3], 2.001)];
    assert!(matches!(
        call(
            &different_time,
            &[contact],
            RodRigidNetworkBudget::default(),
            &gate
        ),
        Err(RodRigidMultiImpactError::Native(
            RodRigidImpactError::Invalid("contact network rod clocks differ")
        ))
    ));
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        call(
            &rods,
            &[contact],
            RodRigidNetworkBudget::default(),
            &cancelled
        )
        .is_err()
    );
    let first = call(&rods, &[contact], RodRigidNetworkBudget::default(), &gate).unwrap();
    let retry = call(&rods, &[contact], RodRigidNetworkBudget::default(), &gate).unwrap();
    for i in 0..2 {
        assert_eq!(
            first.0[i].linear_momenta_kg_m_s(),
            retry.0[i].linear_momenta_kg_m_s()
        );
    }
    audit_network(&first.2);

    // A rod-only redundant cycle determines native velocities, not unique
    // individual edge loads. Distinct rods may all use native node index one.
    let cycle_rods = [
        rods[0].clone(),
        rods[1].clone(),
        tip_rod(2.4, point, [-1.0, 0.0, 0.0], 2.0),
    ];
    let cycle = [
        contact,
        network_edge(
            rod_tip(2, [0.0; 3]),
            rod_tip(1, [0.0; 3]),
            [1.0, 0.0, 0.0],
            0.0,
        ),
        network_edge(
            rod_tip(2, [0.0; 3]),
            rod_tip(0, [0.0; 3]),
            [1.0, 0.0, 0.0],
            0.0,
        ),
    ];
    let cycle_result = call(&cycle_rods, &cycle, RodRigidNetworkBudget::default(), &gate).unwrap();
    assert!(matches!(
        cycle_result.2.contact.determinacy,
        MultiImpactDeterminacy::UniqueVelocity(_)
    ));
    let mean = cycle_rods
        .iter()
        .map(|r| r.linear_momenta_kg_m_s()[1][0])
        .sum::<f64>()
        / cycle_rods
            .iter()
            .map(|r| r.model().inertia()[1].mass_kg())
            .sum::<f64>();
    for rod in &cycle_result.0 {
        near(rod.velocities_m_s().unwrap()[1][0], mean, 2e-11);
    }
    audit_network(&cycle_result.2);

    let finite = [participant(2.3, point, [3.0, 0.0, 0.0])];
    let legacy = edge(
        rod_point(1, [0.0; 3]),
        rigid_point(0, [0.0; 3]),
        [1.0, 0.0, 0.0],
        0.4,
    );
    let old = rods[0]
        .apply_rigid_normal_contact_graph(
            &finite,
            &[legacy],
            impact_controls(),
            graph_config(),
            RodRigidGraphBudget::default(),
            &gate,
        )
        .unwrap();
    let single = apply_rigid_normal_contact_network(
        &rods[..1],
        &finite,
        &[legacy.into()],
        impact_controls(),
        graph_config(),
        RodRigidNetworkBudget::default(),
        &gate,
    )
    .unwrap();
    assert_eq!(
        old.0.linear_momenta_kg_m_s(),
        single.0[0].linear_momenta_kg_m_s()
    );
    assert_eq!(old.1, single.1);
    assert_eq!(
        old.2.contact.normal_impulses,
        single.2.contact.normal_impulses
    );
    let no_rods: &[RodMotion] = &[];
    let bodies = [
        participant(1.0, point, [1.0, 0.0, 0.0]),
        participant(2.0, point, [0.0; 3]),
    ];
    let only_bodies = apply_rigid_normal_contact_network(
        no_rods,
        &bodies,
        &[network_edge(
            finite_point(1, [0.0; 3]),
            finite_point(0, [0.0; 3]),
            [1.0, 0.0, 0.0],
            1.0,
        )],
        impact_controls(),
        graph_config(),
        RodRigidNetworkBudget::default(),
        &gate,
    )
    .unwrap();
    assert!(only_bodies.0.is_empty() && only_bodies.2.transfer.rod_impulses.is_empty());
    near(only_bodies.2.contact.normal_impulses[0], 4.0 / 3.0, 2e-11);
    audit_network(&only_bodies.2);
}
