//! G0/G1/G3/G4: coupled Coulomb impulses across independent native rod owners.

use super::super::super::*;
use fs_contact::frictional_multicontact::{
    FrictionalContact, FrictionalContactConfig, FrictionalContactDeterminacy,
    FrictionalContactError, NormalVelocityLaw, solve_frictional_contacts,
};
use fs_couple::rod_motion::rigid_impact::friction::{
    RodRigidFrictionalControls, RodRigidFrictionalError,
};
use fs_couple::rod_motion::rigid_impact::graph::friction::RodRigidFrictionalGraphContact;
use fs_couple::rod_motion::rigid_impact::graph::network::friction::{
    RodRigidFrictionalNetworkContact, RodRigidFrictionalNetworkReport,
    apply_rigid_frictional_contact_network as apply,
};
use fs_couple::rod_motion::rigid_impact::graph::network::{
    RodRigidNetworkBudget, RodRigidNetworkContact, RodRigidNetworkEndpoint,
};
use fs_couple::rod_motion::rigid_impact::graph::{
    RodRigidEndpoint, RodRigidGraphBudget, RodRigidGraphContact, RodRigidParticipant,
};

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

fn rod_point(rod: usize, node: usize, point_body_m: [f64; 3]) -> RodRigidNetworkEndpoint {
    RodRigidNetworkEndpoint::Rod {
        rod,
        node,
        point_body_m,
    }
}

fn rigid_point(body: usize, arm_body_m: [f64; 3]) -> RodRigidNetworkEndpoint {
    RodRigidNetworkEndpoint::Rigid { body, arm_body_m }
}

fn edge(
    positive: RodRigidNetworkEndpoint,
    negative: RodRigidNetworkEndpoint,
    restitution: f64,
    friction_coefficient: f64,
) -> RodRigidFrictionalNetworkContact {
    RodRigidFrictionalNetworkContact {
        geometry: RodRigidNetworkContact {
            positive,
            negative,
            normal_world: [1.0, 0.0, 0.0],
            activation_tolerance_m: 1e-8,
            restitution,
        },
        tangent_world: [0.0, 1.0, 0.0],
        friction_coefficient,
    }
}

fn centered(length: f64, tip: [f64; 3], velocity: [f64; 3], time: f64) -> RodMotion {
    let gate = CancelGate::new();
    let model = Arc::new(
        SpatialMaterialRod::try_straight(length, vec![section(1.0)])
            .unwrap()
            .into_motion_model(RodRootBoundary::Free, &gate)
            .unwrap(),
    );
    let mut pose = model.initial_pose();
    let shift = sub(tip, pose.positions_m[1]);
    for position in &mut pose.positions_m {
        *position = add(*position, shift);
    }
    RodMotion::new(model, pose, &[velocity; 2], &[[0.0; 3]; 2], time, &gate).unwrap()
}

fn participant(mass: f64, point: [f64; 3], velocity: [f64; 3]) -> RodRigidParticipant {
    let inertia = inertia(mass);
    RodRigidParticipant {
        state: body(inertia, point, velocity),
        inertia,
    }
}

fn bits<I: fs_couple::rod_motion::RodMotionInertia>(rod: &RodMotion<I>) -> Vec<u64> {
    rod.linear_momenta_kg_m_s()
        .iter()
        .chain(rod.body_spins_kg_m2_s())
        .flatten()
        .map(|x| x.to_bits())
        .collect()
}

fn unchanged<I: fs_couple::rod_motion::RodMotionInertia>(
    before: &RodMotion<I>,
    after: &RodMotion<I>,
) {
    assert_eq!(bits(before), bits(after));
    assert_eq!(before.pose(), after.pose());
    assert_eq!(before.time_s().to_bits(), after.time_s().to_bits());
    assert_eq!(before.accepted_steps(), after.accepted_steps());
    assert!(std::ptr::eq(before.model(), after.model()));
}

fn audit(report: &RodRigidFrictionalNetworkReport, contacts: &[RodRigidFrictionalNetworkContact]) {
    for (i, contact) in contacts.iter().enumerate() {
        let impulse = report.contact.impulses[i];
        let velocity = report.transfer.relative_velocities_after_m_s[i];
        let mu = contact.friction_coefficient;
        assert!(impulse[0] >= -2e-10);
        if report.gaps_m[i] <= contact.geometry.activation_tolerance_m {
            assert!(velocity[0] >= report.contact.normal_targets[i] - 2e-9);
            near(
                impulse[0] * (velocity[0] - report.contact.normal_targets[i]),
                0.0,
                2e-9,
            );
        }
        assert!(impulse[1].hypot(impulse[2]) <= mu * impulse[0] + 2e-9);
        let slip = velocity[1].hypot(velocity[2]);
        if slip > 2e-9 {
            near(impulse[1], -mu * impulse[0] * velocity[1] / slip, 2e-9);
            near(impulse[2], -mu * impulse[0] * velocity[2] / slip, 2e-9);
        }
        vector_near(
            report.transfer.endpoint_velocity_residuals_m_s[i],
            [0.0; 3],
            2e-9,
        );
    }
    assert!(report.contact.dissipated_energy >= 0.0);
    if let Some(balance) = report.transfer.balance {
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
fn two_independent_rods_and_a_body_resolve_global_stick_slip_and_continue() {
    let gate = CancelGate::new();
    let point = [0.8, 0.0, 0.0];
    let rods = [
        centered(0.8, point, [-2.0, 1.0, 0.0], 2.0),
        centered(1.6, point, [-0.5, -1.0, 0.0], 2.0),
        // This unused owner's independent clock cannot constrain this event.
        centered(0.4, [2.0, -1.0, 0.5], [-0.3, 0.2, 0.1], 7.0),
    ];
    let mass = rods[0].model().inertia()[1].mass_kg();
    near(rods[1].model().inertia()[1].mass_kg(), 2.0 * mass, 1e-12);
    let bodies = [participant(3.0 * mass, point, [1.0, 1.0 / 3.0, 0.0])];
    for restitution in [0.0, 0.4, 1.0] {
        for mu in [0.0, 0.1, 1.0] {
            let contacts = [
                edge(
                    rod_point(0, 1, [0.0; 3]),
                    rod_point(1, 1, [0.0; 3]),
                    restitution,
                    mu,
                ),
                edge(
                    rod_point(1, 1, [0.0; 3]),
                    rigid_point(0, [0.0; 3]),
                    restitution,
                    mu,
                ),
            ];
            let (next, states, report) = apply(
                &rods,
                &bodies,
                &contacts,
                bounds(),
                RodRigidNetworkBudget::default(),
                &gate,
            )
            .unwrap();
            // Independent three-mass collision equations. The tangential
            // solution either saturates both cones or reaches common speed.
            let normal = [
                2.0 * mass * (1.0 + restitution),
                3.0 * mass * (1.0 + restitution),
            ];
            let tangent = if mu < 1.0 {
                [-mu * normal[0], mu * normal[1]]
            } else {
                [-mass, mass]
            };
            let impulses = [[normal[0], tangent[0], 0.0], [normal[1], tangent[1], 0.0]];
            let velocity = [
                add([-2.0, 1.0, 0.0], impulses[0].map(|p| p / mass)),
                add(
                    [-0.5, -1.0, 0.0],
                    sub(impulses[1], impulses[0]).map(|p| p / (2.0 * mass)),
                ),
                sub([1.0, 1.0 / 3.0, 0.0], impulses[1].map(|p| p / (3.0 * mass))),
            ];
            for i in 0..2 {
                vector_near(report.transfer.impulses_world_n_s[i], impulses[i], 2e-9);
                vector_near(next[i].velocities_m_s().unwrap()[1], velocity[i], 2e-9);
                vector_near(next[i].body_spins_kg_m2_s()[1], [0.0; 3], 2e-11);
                near(
                    report.transfer.relative_velocities_after_m_s[i][0],
                    1.5 * restitution,
                    2e-9,
                );
                assert!(report.transfer.rod_impulses[i].is_some());
                assert_eq!(next[i].pose(), rods[i].pose());
                assert_eq!(next[i].time_s(), rods[i].time_s());
                assert_eq!(next[i].accepted_steps(), rods[i].accepted_steps());
                assert!(std::ptr::eq(next[i].model(), rods[i].model()));
            }
            vector_near(
                array(
                    states[0]
                        .center_of_mass_velocity_world(bodies[0].inertia)
                        .unwrap(),
                ),
                velocity[2],
                2e-9,
            );
            unchanged(&rods[2], &next[2]);
            assert!(report.transfer.rod_impulses[2].is_none());
            let energy = |v: [[f64; 3]; 3]| {
                (0..3)
                    .map(|i| {
                        0.5 * (i as f64 + 1.0) * mass * v[i].iter().map(|x| x * x).sum::<f64>()
                    })
                    .sum::<f64>()
            };
            near(
                report.contact.dissipated_energy,
                energy([[-2.0, 1.0, 0.0], [-0.5, -1.0, 0.0], [1.0, 1.0 / 3.0, 0.0]])
                    - energy(velocity),
                2e-9,
            );
            audit(&report, &contacts);

            if restitution == 0.4 && mu == 1.0 {
                let reordered = [rods[1].clone(), rods[0].clone(), rods[2].clone()];
                let reversed = [
                    edge(
                        rod_point(0, 1, [0.0; 3]),
                        rigid_point(0, [0.0; 3]),
                        restitution,
                        mu,
                    ),
                    edge(
                        rod_point(1, 1, [0.0; 3]),
                        rod_point(0, 1, [0.0; 3]),
                        restitution,
                        mu,
                    ),
                ];
                let permuted = apply(
                    &reordered,
                    &bodies,
                    &reversed,
                    bounds(),
                    RodRigidNetworkBudget::default(),
                    &gate,
                )
                .unwrap();
                for i in 0..2 {
                    vector_near(
                        permuted.0[1 - i].linear_momenta_kg_m_s()[1],
                        next[i].linear_momenta_kg_m_s()[1],
                        2e-9,
                    );
                }
                vector_near(
                    array(permuted.1[0].linear_momentum_world()),
                    array(states[0].linear_momentum_world()),
                    2e-9,
                );
                let first = apply(
                    &rods,
                    &bodies,
                    &contacts[..1],
                    bounds(),
                    RodRigidNetworkBudget::default(),
                    &gate,
                )
                .unwrap();
                let intermediate = [RodRigidParticipant {
                    state: first.1[0],
                    inertia: bodies[0].inertia,
                }];
                let sequential = apply(
                    &first.0,
                    &intermediate,
                    &contacts[1..],
                    bounds(),
                    RodRigidNetworkBudget::default(),
                    &gate,
                )
                .unwrap();
                assert!(
                    (sequential.0[0].linear_momenta_kg_m_s()[1][0]
                        - next[0].linear_momenta_kg_m_s()[1][0])
                        .abs()
                        > 0.1 * mass
                );
                for i in 0..2 {
                    assert_eq!(
                        next[i]
                            .step_candidate(
                                1e-6,
                                &[RodNodalWrench::default(); 2],
                                controls(1e-6),
                                &gate
                            )
                            .unwrap()
                            .0
                            .accepted_steps(),
                        rods[i].accepted_steps() + 1
                    );
                }
            }
        }
    }
}

#[test]
fn eccentric_cross_rod_network_matches_independent_origin_mass_and_frame_changes() {
    let gate = CancelGate::new();
    let models = [
        Arc::new(
            eccentric_source()
                .into_spatial_motion_model(RodRootBoundary::Free, &gate)
                .unwrap(),
        ),
        Arc::new(
            eccentric_source()
                .into_spatial_motion_model(RodRootBoundary::Free, &gate)
                .unwrap(),
        ),
    ];
    let body_tensor = Mat3 {
        m: [0.12, -0.008, 0.006, -0.008, 0.11, 0.005, 0.006, 0.005, 0.14],
    };
    let body_mass = MassProperties::from_inertia_tensor(3.1, BodyVec::ZERO, body_tensor).unwrap();
    let centers = [[0.8, 0.0, 0.0], [0.65, 0.1, -0.1], [0.9, -0.08, 0.05]];
    let axes = [
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [[0.6, 0.8, 0.0], [0.0, 0.0, 1.0], [0.8, -0.6, 0.0]],
    ];
    let arms = [
        [0.04, 0.06, -0.03],
        sub(
            sub(add(centers[0], [0.04, 0.06, -0.03]), centers[1]),
            [1e-13, 0.0, 0.0],
        ),
        [-0.04, 0.06, 0.05],
        sub(
            sub(add(centers[1], [-0.04, 0.06, 0.05]), centers[2]),
            axes[1][0].map(|x| 1e-13 * x),
        ),
    ];
    let mut matrices = [[0.0; 36]; 3];
    for owner in 0..2 {
        let nodal = models[owner].inertia()[1];
        let center = nodal.center_of_mass_body_m();
        let tensor = nodal.inertia_com_body_kg_m2();
        matrices[owner] = fs_mbd::articulated::SpatialInertia::new(
            nodal.mass_kg(),
            Vec3::new(center[0], center[1], center[2]),
            Mat3 {
                m: std::array::from_fn(|i| tensor[i / 3][i % 3]),
            },
        )
        .unwrap()
        .matrix()
        .unwrap()
        .m;
    }
    matrices[2] = fs_mbd::articulated::SpatialInertia::new(
        body_mass.mass(),
        Vec3::new(0.0, 0.0, 0.0),
        body_tensor,
    )
    .unwrap()
    .matrix()
    .unwrap()
    .m;
    // Independent full six-dimensional kinetic blocks about each original
    // material node. The production path acts on the native central tensors.
    let mut rows = [[[0.0; 6]; 6]; 3];
    let mut response = [[[0.0; 6]; 6]; 3];
    for row in 0..6 {
        let contact = row / 3;
        let axis = axes[contact][row % 3];
        let endpoints = if contact == 0 {
            [(0, arms[0], 1.0), (1, arms[1], -1.0)]
        } else {
            [(1, arms[2], 1.0), (2, arms[3], -1.0)]
        };
        for (owner, arm, sign) in endpoints {
            let force = axis.map(|x| sign * x);
            let torque = cross(arm, force);
            rows[owner][row] = [
                torque[0], torque[1], torque[2], force[0], force[1], force[2],
            ];
        }
    }
    for owner in 0..3 {
        let factor = fs_la::factor::lu(&matrices[owner], 6).unwrap();
        for row in 0..6 {
            response[owner][row] = rows[owner][row];
            factor.solve(&mut response[owner][row]);
        }
    }
    let mut mobility = [0.0; 36];
    for i in 0..6 {
        for j in 0..6 {
            mobility[6 * i + j] = (0..3)
                .flat_map(|owner| (0..6).map(move |c| rows[owner][i][c] * response[owner][j][c]))
                .sum();
        }
    }
    assert!(mobility[3..6].iter().any(|x| x.abs() > 1e-3));
    assert!(mobility[1].abs() + mobility[2].abs() > 1e-3);
    let contact_factor = fs_la::factor::lu(&mobility, 6).unwrap();
    let energy = |twists: [[f64; 6]; 3]| {
        (0..3)
            .map(|owner| {
                (0..6)
                    .flat_map(|i| {
                        (0..6).map(move |j| {
                            0.5 * twists[owner][i] * matrices[owner][6 * i + j] * twists[owner][j]
                        })
                    })
                    .sum::<f64>()
            })
            .sum::<f64>()
    };
    for sliding in [false, true] {
        let impulses = if sliding {
            [[0.3, -0.03, 0.0], [0.25, 0.0, 0.0375]]
        } else {
            [[0.3, -0.03, 0.015], [0.25, 0.02, -0.0125]]
        };
        let mu = if sliding { [0.1, 0.15] } else { [0.4; 2] };
        let mut post_coefficients = if sliding {
            [0.0, 0.7, 0.0, 0.0, 0.0, -0.5]
        } else {
            [0.0; 6]
        };
        contact_factor.solve(&mut post_coefficients);
        let mut before = [[0.0; 6]; 3];
        let mut after = [[0.0; 6]; 3];
        for owner in 0..3 {
            for c in 0..6 {
                for row in 0..6 {
                    after[owner][c] += response[owner][row][c] * post_coefficients[row];
                    before[owner][c] += response[owner][row][c]
                        * (post_coefficients[row] - impulses[row / 3][row % 3]);
                }
            }
        }
        for row in [0, 3] {
            let speed: f64 = (0..3)
                .flat_map(|owner| (0..6).map(move |c| rows[owner][row][c] * before[owner][c]))
                .sum();
            assert!(speed < -1e-3);
        }
        let q = So3::exp(So3Tangent::new(Vec3::new(0.4, -0.2, 0.3))).unwrap();
        for (rotation, shift, boost, angles, order, rod_order) in [
            (
                So3::identity(),
                [0.0; 3],
                [0.0; 3],
                [0.0_f64; 2],
                [0, 1],
                [0, 1],
            ),
            (
                q,
                [1.2, -0.4, 0.3],
                [0.3, -0.2, 0.5],
                [0.4_f64, -0.7],
                [1, 0],
                [1, 0],
            ),
        ] {
            let source_rods: [SpatialRodMotion; 2] = std::array::from_fn(|owner| {
                let mut pose = models[owner].initial_pose();
                let translation = sub(centers[owner], pose.positions_m[1]);
                for node in 0..2 {
                    pose.positions_m[node] = add(
                        rotate(rotation, add(pose.positions_m[node], translation)),
                        shift,
                    );
                    pose.frames[node] = rotation;
                }
                SpatialRodMotion::new_spatial(
                    models[owner].clone(),
                    pose,
                    &[
                        boost,
                        add(
                            rotate(
                                rotation,
                                [before[owner][3], before[owner][4], before[owner][5]],
                            ),
                            boost,
                        ),
                    ],
                    &[
                        [0.0; 3],
                        [before[owner][0], before[owner][1], before[owner][2]],
                    ],
                    0.0,
                    1e-12,
                    &gate,
                )
                .unwrap()
            });
            let rods = rod_order.map(|owner| source_rods[owner].clone());
            let rod_index: [usize; 2] =
                std::array::from_fn(|owner| rod_order.iter().position(|&i| i == owner).unwrap());
            let bodies = [RodRigidParticipant {
                inertia: body_mass,
                state: RigidBodyState::new(
                    Pose::new(
                        vector(add(rotate(rotation, centers[2]), shift)),
                        orientation(rotation),
                    )
                    .unwrap(),
                    vector(add(
                        rotate(rotation, [before[2][3], before[2][4], before[2][5]]),
                        boost,
                    ))
                    .scale(body_mass.mass()),
                    body_mass
                        .angular_momentum_body_checked(vector([
                            before[2][0],
                            before[2][1],
                            before[2][2],
                        ]))
                        .unwrap(),
                )
                .unwrap(),
            }];
            let endpoints = [
                (
                    rod_point(rod_index[0], 1, arms[0]),
                    rod_point(rod_index[1], 1, arms[1]),
                ),
                (rod_point(rod_index[1], 1, arms[2]), rigid_point(0, arms[3])),
            ];
            let contacts = order.map(|i| RodRigidFrictionalNetworkContact {
                geometry: RodRigidNetworkContact {
                    positive: endpoints[i].0,
                    negative: endpoints[i].1,
                    normal_world: rotate(rotation, axes[i][0]),
                    activation_tolerance_m: 1e-8,
                    restitution: 0.0,
                },
                tangent_world: rotate(
                    rotation,
                    add(
                        axes[i][1].map(|x| x * angles[i].cos()),
                        axes[i][2].map(|x| x * angles[i].sin()),
                    ),
                ),
                friction_coefficient: mu[i],
            });
            let (next, states, report) = apply(
                &rods,
                &bodies,
                &contacts,
                bounds(),
                RodRigidNetworkBudget::default(),
                &gate,
            )
            .unwrap();
            for row in 0..6 {
                let owner = order[row / 3];
                let a = rotate(rotation.inverse(), report.frames_world[row / 3][row % 3]);
                let left: [f64; 3] =
                    std::array::from_fn(|k| (0..3).map(|c| a[c] * axes[owner][k][c]).sum());
                for col in 0..6 {
                    let other = order[col / 3];
                    let b = rotate(rotation.inverse(), report.frames_world[col / 3][col % 3]);
                    let right: [f64; 3] =
                        std::array::from_fn(|k| (0..3).map(|c| b[c] * axes[other][k][c]).sum());
                    let expected: f64 = (0..3)
                        .flat_map(|i| {
                            (0..3).map(move |j| {
                                left[i] * mobility[(3 * owner + i) * 6 + 3 * other + j] * right[j]
                            })
                        })
                        .sum();
                    near(report.delassus_kg_inverse[6 * row + col], expected, 3e-10);
                }
            }
            for owner in 0..2 {
                vector_near(
                    next[rod_index[owner]].velocities_m_s().unwrap()[1],
                    add(
                        rotate(
                            rotation,
                            [after[owner][3], after[owner][4], after[owner][5]],
                        ),
                        boost,
                    ),
                    2e-9,
                );
                vector_near(
                    next[rod_index[owner]]
                        .angular_velocities_body_rad_s()
                        .unwrap()[1],
                    [after[owner][0], after[owner][1], after[owner][2]],
                    2e-9,
                );
                assert!(report.transfer.rod_impulses[rod_index[owner]].is_some());
                let contact = order[owner];
                let expected = std::array::from_fn(|c| {
                    (0..3)
                        .map(|a| axes[contact][a][c] * impulses[contact][a])
                        .sum()
                });
                vector_near(
                    report.transfer.impulses_world_n_s[owner],
                    rotate(rotation, expected),
                    2e-9,
                );
            }
            vector_near(
                array(states[0].center_of_mass_velocity_world(body_mass).unwrap()),
                add(
                    rotate(rotation, [after[2][3], after[2][4], after[2][5]]),
                    boost,
                ),
                2e-9,
            );
            vector_near(
                array(
                    body_mass
                        .angular_velocity_body_checked(states[0].angular_momentum_body())
                        .unwrap(),
                ),
                [after[2][0], after[2][1], after[2][2]],
                2e-9,
            );
            assert!(matches!(
                report.contact.determinacy,
                FrictionalContactDeterminacy::UniqueImpulses
            ));
            near(
                report.contact.dissipated_energy,
                energy(before) - energy(after),
                2e-9,
            );
            audit(&report, &contacts);
        }
    }
}

#[test]
fn cross_rod_supports_and_redundant_edges_retain_original_work_and_determinacy() {
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
    let rods = [
        RodMotion::new(
            model.clone(),
            model.initial_pose(),
            &[[0.0; 3]; 2],
            &[[0.0; 3]; 2],
            2.0,
            &gate,
        )
        .unwrap(),
        centered(0.8, [0.8, 0.07, 0.0], [1.0, 0.3, -0.2], 2.0),
    ];
    let bodies = [participant(2.1, [0.8, 0.07, 0.0], [2.0, 0.4, 0.1])];
    let contacts = [
        edge(
            rod_point(0, 1, [0.0, 0.07, 0.0]),
            rod_point(1, 1, [0.0; 3]),
            0.0,
            1.0,
        ),
        edge(
            rod_point(1, 1, [0.0; 3]),
            rigid_point(0, [0.0; 3]),
            0.0,
            1.0,
        ),
    ];
    let (next, states, report) = apply(
        &rods,
        &bodies,
        &contacts,
        bounds(),
        RodRigidNetworkBudget::default(),
        &gate,
    )
    .unwrap();
    let delivered = add(
        rods[1].linear_momenta_kg_m_s()[1],
        array(bodies[0].state.linear_momentum_world()),
    );
    let support = &report.transfer.rod_impulses[0]
        .as_ref()
        .unwrap()
        .support_impulses[1];
    vector_near(support.linear_n_s, delivered.map(|p| -p), 2e-9);
    vector_near(
        support.angular_nm_s,
        cross([0.0, 0.07, 0.0], delivered).map(|p| -p),
        2e-9,
    );
    unchanged(&rods[0], &next[0]);
    vector_near(next[1].linear_momenta_kg_m_s()[1], [0.0; 3], 2e-9);
    vector_near(next[1].body_spins_kg_m2_s()[1], [0.0; 3], 2e-9);
    vector_near(array(states[0].linear_momentum_world()), [0.0; 3], 2e-9);
    assert_eq!(
        next[1].linear_momenta_kg_m_s()[0],
        rods[1].linear_momenta_kg_m_s()[0]
    );
    let loss = 0.5 * rods[1].model().inertia()[1].mass_kg() * (1.0 + 0.3 * 0.3 + 0.2 * 0.2)
        + bodies[0].state.kinetic_energy(bodies[0].inertia).unwrap();
    near(report.contact.dissipated_energy, loss, 2e-9);
    audit(&report, &contacts);

    let rods = [
        centered(0.8, [0.8, 0.0, 0.0], [-0.6, 0.8, -0.6], 2.0),
        centered(1.6, [0.8, 0.0, 0.0], [1.4, 0.0, 0.0], 2.0),
    ];
    let contact = edge(
        rod_point(0, 1, [0.0; 3]),
        rod_point(1, 1, [0.0; 3]),
        0.35,
        0.1,
    );
    let one = apply(
        &rods,
        &[],
        &[contact],
        bounds(),
        RodRigidNetworkBudget::default(),
        &gate,
    )
    .unwrap();
    let repeated = apply(
        &rods,
        &[],
        &[contact; 2],
        bounds(),
        RodRigidNetworkBudget::default(),
        &gate,
    )
    .unwrap();
    let FrictionalContactDeterminacy::UniqueVelocity(ref range) = repeated.2.contact.determinacy
    else {
        panic!("duplicate cross-rod edges cannot determine individual loads");
    };
    assert_eq!((range.dimension, range.rank), (6, 3));
    assert!(repeated.2.contact.impulse_work_witness.is_some());
    for i in 0..2 {
        vector_near(
            repeated.0[i].linear_momenta_kg_m_s()[1],
            one.0[i].linear_momenta_kg_m_s()[1],
            2e-9,
        );
    }
    near(
        repeated.2.contact.dissipated_energy,
        one.2.contact.dissipated_energy,
        2e-9,
    );
    audit(&repeated.2, &[contact; 2]);
    let mut incompatible = [contact; 2];
    incompatible[1].friction_coefficient = 0.2;
    let original: Vec<_> = incompatible
        .iter()
        .enumerate()
        .map(|(i, edge)| FrictionalContact {
            id: i as u64,
            gap: repeated.2.gaps_m[i],
            activation_tolerance: edge.geometry.activation_tolerance_m,
            relative_velocity: repeated.2.relative_velocities_before_m_s[i],
            friction_coefficient: edge.friction_coefficient,
            normal_law: NormalVelocityLaw::NewtonRestitution {
                coefficient: edge.geometry.restitution,
            },
        })
        .collect();
    let error = solve_frictional_contacts(
        &original,
        &repeated.2.delassus_kg_inverse,
        bounds().contact,
        || false,
    )
    .unwrap_err();
    assert!(matches!(
        error,
        FrictionalContactError::Nonunique { .. }
            | FrictionalContactError::UnresolvedDeterminacy { .. }
    ));
    let RodRigidFrictionalError::Contact(refused) = apply(
        &rods,
        &[],
        &incompatible,
        bounds(),
        RodRigidNetworkBudget::default(),
        &gate,
    )
    .unwrap_err() else {
        panic!("the original law's ambiguity must remain a contact refusal");
    };
    assert_eq!(
        std::mem::discriminant(&error),
        std::mem::discriminant(&refused)
    );
}

#[test]
fn independent_clocks_budgets_refusals_and_empty_owner_limits_preserve_native_states() {
    let gate = CancelGate::new();
    let point = [0.8, 0.0, 0.0];
    let rods = [
        centered(0.8, point, [-0.6, 0.8, -0.6], 2.0),
        centered(1.6, point, [1.4, 0.0, 0.0], 2.0),
    ];
    let original = rods.each_ref().map(bits);
    let contact = edge(
        rod_point(0, 1, [0.0; 3]),
        rod_point(1, 1, [0.0; 3]),
        0.35,
        0.1,
    );
    for budget in [
        RodRigidNetworkBudget {
            maximum_rods: 1,
            ..Default::default()
        },
        RodRigidNetworkBudget {
            maximum_rod_nodes: 3,
            ..Default::default()
        },
        RodRigidNetworkBudget {
            maximum_rods: 129,
            ..Default::default()
        },
    ] {
        assert!(apply(&rods, &[], &[contact], bounds(), budget, &gate).is_err());
    }
    assert!(
        apply(
            &rods,
            &[],
            &[contact; 22],
            bounds(),
            RodRigidNetworkBudget::default(),
            &gate
        )
        .is_err()
    );
    let mut short = bounds();
    short.contact.matrix.max_variables = 2;
    assert!(
        apply(
            &rods,
            &[],
            &[contact],
            short,
            RodRigidNetworkBudget::default(),
            &gate
        )
        .is_err()
    );
    short = bounds();
    short.contact.maximum_evaluations = 1;
    assert!(
        apply(
            &rods,
            &[],
            &[contact],
            short,
            RodRigidNetworkBudget::default(),
            &gate
        )
        .is_err()
    );
    for endpoint in [
        rod_point(2, 1, [0.0; 3]),
        rod_point(1, 2, [0.0; 3]),
        rod_point(0, 1, [0.0; 3]),
    ] {
        let mut invalid = contact;
        invalid.geometry.negative = endpoint;
        assert!(
            apply(
                &rods,
                &[],
                &[invalid],
                bounds(),
                RodRigidNetworkBudget::default(),
                &gate
            )
            .is_err()
        );
    }
    let mut invalid = contact;
    invalid.tangent_world = invalid.geometry.normal_world;
    assert!(
        apply(
            &rods,
            &[],
            &[invalid],
            bounds(),
            RodRigidNetworkBudget::default(),
            &gate
        )
        .is_err()
    );
    invalid = contact;
    invalid.friction_coefficient = -0.1;
    assert!(matches!(
        apply(
            &rods,
            &[],
            &[invalid],
            bounds(),
            RodRigidNetworkBudget::default(),
            &gate
        ),
        Err(RodRigidFrictionalError::Contact(_))
    ));
    let mismatch = [rods[0].clone(), centered(1.6, point, [1.4, 0.0, 0.0], 3.0)];
    assert!(matches!(
        apply(
            &mismatch,
            &[],
            &[contact],
            bounds(),
            RodRigidNetworkBudget::default(),
            &gate
        ),
        Err(RodRigidFrictionalError::Impact(
            RodRigidImpactError::Invalid(_)
        ))
    ));
    // No event names these clocks, so each unused native trajectory is retained.
    let (next, states, empty) = apply(
        &mismatch,
        &[],
        &[],
        bounds(),
        RodRigidNetworkBudget::default(),
        &gate,
    )
    .unwrap();
    for i in 0..2 {
        unchanged(&mismatch[i], &next[i]);
    }
    assert!(states.is_empty() && empty.transfer.balance.is_none());
    assert!(empty.transfer.rod_impulses.iter().all(Option::is_none));
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        apply(
            &rods,
            &[],
            &[contact],
            bounds(),
            RodRigidNetworkBudget::default(),
            &cancelled
        )
        .is_err()
    );
    for tip in [[0.8 + 1e-4, 0.0, 0.0], [0.8 - 1e-9, 0.0, 0.0]] {
        let inadmissible = [rods[0].clone(), centered(1.6, tip, [1.4, 0.0, 0.0], 2.0)];
        assert!(
            apply(
                &inadmissible,
                &[],
                &[contact],
                bounds(),
                RodRigidNetworkBudget::default(),
                &gate
            )
            .is_err()
        );
    }
    let crooked = [
        rods[0].clone(),
        centered(1.6, [0.8, 1e-4, 0.0], [1.4, 0.0, 0.0], 2.0),
    ];
    let mut broad_geometry = bounds();
    broad_geometry.impact.geometry_tolerance_m = 1e-3;
    assert!(matches!(
        apply(
            &crooked,
            &[],
            &[contact],
            broad_geometry,
            RodRigidNetworkBudget::default(),
            &gate
        ),
        Err(RodRigidFrictionalError::Impact(
            RodRigidImpactError::Balance { .. }
        ))
    ));
    assert_eq!(original, rods.each_ref().map(bits));
    let first = apply(
        &rods,
        &[],
        &[contact],
        bounds(),
        RodRigidNetworkBudget::default(),
        &gate,
    )
    .unwrap();
    let retry = apply(
        &rods,
        &[],
        &[contact],
        bounds(),
        RodRigidNetworkBudget::default(),
        &gate,
    )
    .unwrap();
    for i in 0..2 {
        assert_eq!(bits(&first.0[i]), bits(&retry.0[i]));
    }
    for quiet in [
        centered(1.6, [0.79, 0.0, 0.0], [1.4, 0.0, 0.0], 2.0),
        centered(1.6, point, [-1.4, 0.0, 0.0], 2.0),
    ] {
        let input = [rods[0].clone(), quiet];
        let out = apply(
            &input,
            &[],
            &[contact],
            bounds(),
            RodRigidNetworkBudget::default(),
            &gate,
        )
        .unwrap();
        for i in 0..2 {
            unchanged(&input[i], &out.0[i]);
        }
        assert!(out.2.transfer.rod_impulses.iter().all(Option::is_none));
        assert!(out.2.transfer.balance.is_none());
    }

    // The same owner now supports a true rigid-only call with no dummy rod.
    let no_rods: [RodMotion; 0] = [];
    let bodies = [
        participant(0.8, point, [-1.0, 0.8, -0.6]),
        participant(1.9, point, [1.0, 0.0, 0.0]),
    ];
    let contact = edge(
        rigid_point(0, [0.0; 3]),
        rigid_point(1, [0.0; 3]),
        0.35,
        0.2,
    );
    let (next, states, report) = apply(
        &no_rods,
        &bodies,
        &[contact],
        bounds(),
        RodRigidNetworkBudget {
            maximum_rods: 0,
            maximum_rod_nodes: 0,
            ..Default::default()
        },
        &gate,
    )
    .unwrap();
    assert!(next.is_empty() && report.transfer.rod_impulses.is_empty());
    let reduced_mass: f64 = 1.0 / (1.0 / 0.8 + 1.0 / 1.9);
    let normal = 1.35 * 2.0 * reduced_mass;
    let tangent = (0.2 * normal).min(reduced_mass);
    let impulse = [normal, -0.8 * tangent, 0.6 * tangent];
    vector_near(report.transfer.impulses_world_n_s[0], impulse, 2e-9);
    for i in 0..2 {
        vector_near(
            array(states[i].linear_momentum_world()),
            add(
                array(bodies[i].state.linear_momentum_world()),
                impulse.map(|p| if i == 0 { p } else { -p }),
            ),
            2e-9,
        );
    }
    audit(&report, &[contact]);

    // Published one-rod API and indexed entry point use identical mechanics.
    let bodies = [participant(2.3, point, [1.4, 0.0, 0.0])];
    let single = RodRigidFrictionalGraphContact {
        geometry: RodRigidGraphContact {
            positive: RodRigidEndpoint::Rod {
                node: 1,
                point_body_m: [0.0; 3],
            },
            negative: RodRigidEndpoint::Rigid {
                body: 0,
                arm_body_m: [0.0; 3],
            },
            normal_world: [1.0, 0.0, 0.0],
            activation_tolerance_m: 1e-8,
            restitution: 0.35,
        },
        tangent_world: [0.0, 1.0, 0.0],
        friction_coefficient: 0.1,
    };
    let old = rods[0]
        .apply_rigid_frictional_contact_graph(
            &bodies,
            &[single],
            bounds(),
            RodRigidGraphBudget::default(),
            &gate,
        )
        .unwrap();
    let indexed = RodRigidFrictionalNetworkContact {
        geometry: single.geometry.into(),
        tangent_world: single.tangent_world,
        friction_coefficient: single.friction_coefficient,
    };
    let new = apply(
        core::slice::from_ref(&rods[0]),
        &bodies,
        &[indexed],
        bounds(),
        RodRigidNetworkBudget::default(),
        &gate,
    )
    .unwrap();
    assert_eq!(bits(&old.0), bits(&new.0[0]));
    assert_eq!(old.1, new.1);
    assert_eq!(old.2.delassus_kg_inverse, new.2.delassus_kg_inverse);
    assert_eq!(old.2.contact.impulses, new.2.contact.impulses);
}
