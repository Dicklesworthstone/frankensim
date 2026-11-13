//! G0/G1/G3/G4: original momentum and Coulomb work on mixed native contact graphs.

use super::super::*;
use fs_contact::frictional_multicontact::{
    FrictionalContact, FrictionalContactConfig, FrictionalContactDeterminacy,
    FrictionalContactError, NormalVelocityLaw, solve_frictional_contacts,
};
use fs_couple::rod_motion::rigid_impact::friction::{
    RodRigidFrictionalControls, RodRigidFrictionalError,
};
use fs_couple::rod_motion::rigid_impact::graph::friction::{
    RodRigidFrictionalGraphContact, RodRigidFrictionalGraphReport,
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

fn rod_point(node: usize, point_body_m: [f64; 3]) -> RodRigidEndpoint {
    RodRigidEndpoint::Rod { node, point_body_m }
}

fn rigid_point(body: usize, arm_body_m: [f64; 3]) -> RodRigidEndpoint {
    RodRigidEndpoint::Rigid { body, arm_body_m }
}

fn contact(
    positive: RodRigidEndpoint,
    negative: RodRigidEndpoint,
    restitution: f64,
    friction_coefficient: f64,
) -> RodRigidFrictionalGraphContact {
    RodRigidFrictionalGraphContact {
        geometry: RodRigidGraphContact {
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

fn participant(mass: f64, position: [f64; 3], velocity: [f64; 3]) -> RodRigidParticipant {
    let inertia = inertia(mass);
    RodRigidParticipant {
        state: body(inertia, position, velocity),
        inertia,
    }
}

fn audit(report: &RodRigidFrictionalGraphReport, contacts: &[RodRigidFrictionalGraphContact]) {
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

fn momentum_bits<I: fs_couple::rod_motion::RodMotionInertia>(
    rod: &fs_couple::rod_motion::RodMotion<I>,
) -> Vec<u64> {
    rod.linear_momenta_kg_m_s()
        .iter()
        .chain(rod.body_spins_kg_m2_s())
        .flatten()
        .map(|value| value.to_bits())
        .collect()
}

#[test]
fn mixed_three_mass_chain_resolves_global_stick_slip_and_keeps_participant_identity() {
    let gate = CancelGate::new();
    let rod = free_rod([-2.0, 1.0, 0.0]);
    let m = rod.model().inertia()[1].mass_kg();
    let bodies = [
        participant(2.0 * m, [0.8, 0.0, 0.0], [-0.5, -1.0, 0.0]),
        participant(3.0 * m, [0.8, 0.0, 0.0], [1.0, 1.0 / 3.0, 0.0]),
        participant(1.7, [2.0, -1.0, 0.5], [-0.3, 0.2, 0.1]),
    ];
    for restitution in [0.0, 0.4, 1.0] {
        for mu in [0.0, 0.1, 1.0] {
            let witnesses = [
                contact(
                    rod_point(1, [0.0; 3]),
                    rigid_point(0, [0.0; 3]),
                    restitution,
                    mu,
                ),
                contact(
                    rigid_point(0, [0.0; 3]),
                    rigid_point(1, [0.0; 3]),
                    restitution,
                    mu,
                ),
            ];
            let (next, states, report) = rod
                .apply_rigid_frictional_contact_graph(
                    &bodies,
                    &witnesses,
                    bounds(),
                    RodRigidGraphBudget::default(),
                    &gate,
                )
                .unwrap();
            let normal = [2.0 * m * (1.0 + restitution), 3.0 * m * (1.0 + restitution)];
            let tangent = if mu < 1.0 {
                [-mu * normal[0], mu * normal[1]]
            } else {
                [-m, m]
            };
            let impulses = [[normal[0], tangent[0], 0.0], [normal[1], tangent[1], 0.0]];
            let velocities = [
                add([-2.0, 1.0, 0.0], impulses[0].map(|p| p / m)),
                add(
                    [-0.5, -1.0, 0.0],
                    sub(impulses[1], impulses[0]).map(|p| p / (2.0 * m)),
                ),
                sub([1.0, 1.0 / 3.0, 0.0], impulses[1].map(|p| p / (3.0 * m))),
            ];
            for i in 0..2 {
                vector_near(report.transfer.impulses_world_n_s[i], impulses[i], 2e-9);
                near(
                    report.transfer.relative_velocities_after_m_s[i][0],
                    1.5 * restitution,
                    2e-9,
                );
                vector_near(
                    array(
                        states[i]
                            .center_of_mass_velocity_world(bodies[i].inertia)
                            .unwrap(),
                    ),
                    velocities[i + 1],
                    2e-9,
                );
                vector_near(array(states[i].angular_momentum_body()), [0.0; 3], 2e-11);
            }
            vector_near(next.velocities_m_s().unwrap()[1], velocities[0], 2e-9);
            assert_eq!(states[2], bodies[2].state);
            assert!(report.transfer.rigid_impulses[2].is_none());
            assert_eq!(next.pose(), rod.pose());
            assert_eq!(next.time_s(), rod.time_s());
            let before = [[-2.0, 1.0, 0.0], [-0.5, -1.0, 0.0], [1.0, 1.0 / 3.0, 0.0]];
            let energy = |v: [[f64; 3]; 3]| {
                (0..3)
                    .map(|i| 0.5 * (i as f64 + 1.0) * m * v[i].iter().map(|x| x * x).sum::<f64>())
                    .sum::<f64>()
            };
            near(
                report.contact.dissipated_energy,
                energy(before) - energy(velocities),
                2e-9,
            );
            audit(&report, &witnesses);

            if restitution == 0.4 && mu == 1.0 {
                let reordered = [bodies[1], bodies[0], bodies[2]];
                let edges = [
                    contact(
                        rigid_point(1, [0.0; 3]),
                        rigid_point(0, [0.0; 3]),
                        restitution,
                        mu,
                    ),
                    contact(
                        rod_point(1, [0.0; 3]),
                        rigid_point(1, [0.0; 3]),
                        restitution,
                        mu,
                    ),
                ];
                let permuted = rod
                    .apply_rigid_frictional_contact_graph(
                        &reordered,
                        &edges,
                        bounds(),
                        RodRigidGraphBudget::default(),
                        &gate,
                    )
                    .unwrap();
                vector_near(
                    permuted.0.linear_momenta_kg_m_s()[1],
                    next.linear_momenta_kg_m_s()[1],
                    2e-9,
                );
                vector_near(
                    array(permuted.1[1].linear_momentum_world()),
                    array(states[0].linear_momentum_world()),
                    2e-9,
                );
                vector_near(
                    array(permuted.1[0].linear_momentum_world()),
                    array(states[1].linear_momentum_world()),
                    2e-9,
                );
                let first = rod
                    .apply_rigid_frictional_contact_graph(
                        &bodies,
                        &witnesses[..1],
                        bounds(),
                        RodRigidGraphBudget::default(),
                        &gate,
                    )
                    .unwrap();
                let intermediate: Vec<_> = bodies
                    .iter()
                    .zip(first.1)
                    .map(|(old, state)| RodRigidParticipant {
                        state,
                        inertia: old.inertia,
                    })
                    .collect();
                let sequential = first
                    .0
                    .apply_rigid_frictional_contact_graph(
                        &intermediate,
                        &witnesses[1..],
                        bounds(),
                        RodRigidGraphBudget::default(),
                        &gate,
                    )
                    .unwrap();
                assert!(
                    (sequential.0.linear_momenta_kg_m_s()[1][0]
                        - next.linear_momenta_kg_m_s()[1][0])
                        .abs()
                        > 0.1 * m
                );
                assert_eq!(
                    next.step_candidate(
                        1e-6,
                        &[RodNodalWrench::default(); 2],
                        controls(1e-6),
                        &gate
                    )
                    .unwrap()
                    .0
                    .accepted_steps(),
                    rod.accepted_steps() + 1,
                );
                let continued = fs_mbd::RigidBodyIntegrator::new(fs_mbd::Gravity::ZERO)
                    .step(states[0], bodies[0].inertia, fs_mbd::Wrench::ZERO, 1e-6)
                    .unwrap();
                assert_ne!(continued.state_after.pose(), states[0].pose());
            }
        }
    }
}

#[test]
fn full_spatial_mass_network_oracle_preserves_eccentric_spin_and_frame_covariance() {
    let gate = CancelGate::new();
    let model = Arc::new(
        eccentric_source()
            .into_spatial_motion_model(RodRootBoundary::Free, &gate)
            .unwrap(),
    );
    let tensors = [
        Mat3 {
            m: [0.07, 0.01, -0.004, 0.01, 0.09, 0.006, -0.004, 0.006, 0.08],
        },
        Mat3 {
            m: [0.12, -0.008, 0.006, -0.008, 0.11, 0.005, 0.006, 0.005, 0.14],
        },
    ];
    let masses = [
        MassProperties::from_inertia_tensor(2.4, BodyVec::ZERO, tensors[0]).unwrap(),
        MassProperties::from_inertia_tensor(3.1, BodyVec::ZERO, tensors[1]).unwrap(),
    ];
    let centers = [[0.65, 0.1, -0.1], [0.9, -0.08, 0.05]];
    let axes = [
        [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        [[0.6, 0.8, 0.0], [0.0, 0.0, 1.0], [0.8, -0.6, 0.0]],
    ];
    let rod_arm = [0.04, 0.06, -0.03];
    let first_arm = sub(
        sub(add([0.8, 0.0, 0.0], rod_arm), centers[0]),
        [1e-13, 0.0, 0.0],
    );
    let second_arm = [-0.04, 0.06, 0.05];
    let last_arm = sub(
        sub(add(centers[0], second_arm), centers[1]),
        axes[1][0].map(|x| 1e-13 * x),
    );
    let nodal = model.inertia()[1];
    let centroid = nodal.center_of_mass_body_m();
    let tensor = nodal.inertia_com_body_kg_m2();
    let mut matrices = [[0.0; 36]; 3];
    matrices[0] = fs_mbd::articulated::SpatialInertia::new(
        nodal.mass_kg(),
        Vec3::new(centroid[0], centroid[1], centroid[2]),
        Mat3 {
            m: std::array::from_fn(|i| tensor[i / 3][i % 3]),
        },
    )
    .unwrap()
    .matrix()
    .unwrap()
    .m;
    for body in 0..2 {
        matrices[body + 1] = fs_mbd::articulated::SpatialInertia::new(
            masses[body].mass(),
            Vec3::new(0.0, 0.0, 0.0),
            tensors[body],
        )
        .unwrap()
        .matrix()
        .unwrap()
        .m;
    }
    // Independent origin-coordinate kinetic blocks and signed point-force
    // transposes. Production instead uses each native central tensor action.
    let mut rows = [[[0.0; 6]; 6]; 3];
    let mut response = [[[0.0; 6]; 6]; 3];
    for row in 0..6 {
        let edge = row / 3;
        let axis = axes[edge][row % 3];
        let endpoints = if edge == 0 {
            [(0, rod_arm, 1.0), (1, first_arm, -1.0)]
        } else {
            [(1, second_arm, 1.0), (2, last_arm, -1.0)]
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
        let post_contact = if sliding {
            [0.0, 0.7, 0.0, 0.0, 0.0, -0.5]
        } else {
            [0.0; 6]
        };
        let mut post_coefficients = post_contact;
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
        for (rotation, shift, boost, angles, order) in [
            (So3::identity(), [0.0; 3], [0.0; 3], [0.0_f64; 2], [0, 1]),
            (
                q,
                [1.2, -0.4, 0.3],
                [0.3, -0.2, 0.5],
                [0.4_f64, -0.7],
                [1, 0],
            ),
        ] {
            let mut pose = model.initial_pose();
            for node in 0..2 {
                pose.positions_m[node] = add(rotate(rotation, pose.positions_m[node]), shift);
                pose.frames[node] = rotation;
            }
            let rod = SpatialRodMotion::new_spatial(
                model.clone(),
                pose,
                &[
                    boost,
                    add(
                        rotate(rotation, [before[0][3], before[0][4], before[0][5]]),
                        boost,
                    ),
                ],
                &[[0.0; 3], [before[0][0], before[0][1], before[0][2]]],
                0.0,
                1e-12,
                &gate,
            )
            .unwrap();
            let bodies: [RodRigidParticipant; 2] = std::array::from_fn(|i| RodRigidParticipant {
                inertia: masses[i],
                state: RigidBodyState::new(
                    Pose::new(
                        vector(add(rotate(rotation, centers[i]), shift)),
                        orientation(rotation),
                    )
                    .unwrap(),
                    vector(add(
                        rotate(
                            rotation,
                            [before[i + 1][3], before[i + 1][4], before[i + 1][5]],
                        ),
                        boost,
                    ))
                    .scale(masses[i].mass()),
                    masses[i]
                        .angular_momentum_body_checked(vector([
                            before[i + 1][0],
                            before[i + 1][1],
                            before[i + 1][2],
                        ]))
                        .unwrap(),
                )
                .unwrap(),
            });
            let endpoints = [
                (rod_point(1, rod_arm), rigid_point(0, first_arm)),
                (rigid_point(0, second_arm), rigid_point(1, last_arm)),
            ];
            let witnesses = order.map(|i| RodRigidFrictionalGraphContact {
                geometry: RodRigidGraphContact {
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
            let (next, states, report) = rod
                .apply_rigid_frictional_contact_graph(
                    &bodies,
                    &witnesses,
                    bounds(),
                    RodRigidGraphBudget::default(),
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
            vector_near(
                next.velocities_m_s().unwrap()[1],
                add(
                    rotate(rotation, [after[0][3], after[0][4], after[0][5]]),
                    boost,
                ),
                2e-9,
            );
            vector_near(
                next.angular_velocities_body_rad_s().unwrap()[1],
                [after[0][0], after[0][1], after[0][2]],
                2e-9,
            );
            for i in 0..2 {
                vector_near(
                    array(states[i].center_of_mass_velocity_world(masses[i]).unwrap()),
                    add(
                        rotate(
                            rotation,
                            [after[i + 1][3], after[i + 1][4], after[i + 1][5]],
                        ),
                        boost,
                    ),
                    2e-9,
                );
                vector_near(
                    array(
                        masses[i]
                            .angular_velocity_body_checked(states[i].angular_momentum_body())
                            .unwrap(),
                    ),
                    [after[i + 1][0], after[i + 1][1], after[i + 1][2]],
                    2e-9,
                );
                let edge = order[i];
                let expected = std::array::from_fn(|c| {
                    (0..3).map(|a| axes[edge][a][c] * impulses[edge][a]).sum()
                });
                vector_near(
                    report.transfer.impulses_world_n_s[i],
                    rotate(rotation, expected),
                    2e-9,
                );
            }
            assert!(matches!(
                report.contact.determinacy,
                FrictionalContactDeterminacy::UniqueImpulses
            ));
            near(
                report.contact.dissipated_energy,
                energy(before) - energy(after),
                2e-9,
            );
            audit(&report, &witnesses);
        }
    }
}

#[test]
fn rigid_only_and_rod_self_contact_preserve_the_original_untouched_states() {
    let gate = CancelGate::new();
    let rod = free_rod([-0.0, 0.0, -0.0]);
    let before_bits = momentum_bits(&rod);
    let bodies = [
        participant(0.8, [1.2, 0.4, -0.1], [-0.6, 0.8, -0.6]),
        participant(1.9, [1.2, 0.4, -0.1], [1.4, 0.0, 0.0]),
        participant(2.7, [3.0, 0.0, 0.0], [-0.2, 0.3, 0.4]),
    ];
    let witness = contact(
        rigid_point(0, [0.0; 3]),
        rigid_point(1, [0.0; 3]),
        0.35,
        0.2,
    );
    let (next, states, report) = rod
        .apply_rigid_frictional_contact_graph(
            &bodies,
            &[witness],
            bounds(),
            RodRigidGraphBudget::default(),
            &gate,
        )
        .unwrap();
    let reduced_mass: f64 = 1.0 / (1.0 / 0.8 + 1.0 / 1.9);
    let normal = 1.35 * 2.0 * reduced_mass;
    let tangent = (0.2 * normal).min(reduced_mass);
    let impulse = [normal, -0.8 * tangent, 0.6 * tangent];
    vector_near(report.transfer.impulses_world_n_s[0], impulse, 2e-9);
    for i in 0..2 {
        let sign = if i == 0 { 1.0 } else { -1.0 };
        vector_near(
            array(states[i].linear_momentum_world()),
            add(
                array(bodies[i].state.linear_momentum_world()),
                impulse.map(|j| sign * j),
            ),
            2e-9,
        );
    }
    assert_eq!(before_bits, momentum_bits(&next));
    assert_eq!(next.pose(), rod.pose());
    assert_eq!(next.time_s(), rod.time_s());
    assert!(report.transfer.rod_impulse.is_none());
    assert_eq!(states[2], bodies[2].state);
    assert!(report.transfer.rigid_impulses[2].is_none());
    audit(&report, &[witness]);

    let model = Arc::new(rod.model().clone());
    let state = RodMotion::new(
        model.clone(),
        model.initial_pose(),
        &[[-1.0, 0.5, 0.0], [1.0, -0.5, 0.0]],
        &[[0.0; 3]; 2],
        0.0,
        &gate,
    )
    .unwrap();
    let m0 = model.inertia()[0].mass_kg();
    let m1 = model.inertia()[1].mass_kg();
    let inverse_mass = 1.0 / m0 + 1.0 / m1;
    let tangent_mobility = inverse_mass
        + 0.16 / model.inertia()[0].principal_kg_m2()[2]
        + 0.16 / model.inertia()[1].principal_kg_m2()[2];
    let normal = 1.4 * 2.0 / inverse_mass;
    let sticking = 1.0 / tangent_mobility;
    for mu in [0.0, 0.25 * sticking / normal, 2.0 * sticking / normal] {
        let witness = contact(
            rod_point(0, [0.4, 0.0, 0.0]),
            rod_point(1, [-0.4, 0.0, 0.0]),
            0.4,
            mu,
        );
        let (next, bodies, report) = state
            .apply_rigid_frictional_contact_graph(
                &[],
                &[witness],
                bounds(),
                RodRigidGraphBudget { maximum_bodies: 0 },
                &gate,
            )
            .unwrap();
        let tangent = -(mu * normal).min(sticking);
        let impulse = [normal, tangent, 0.0];
        vector_near(report.transfer.impulses_world_n_s[0], impulse, 2e-9);
        vector_near(
            next.linear_momenta_kg_m_s()[0],
            add(state.linear_momenta_kg_m_s()[0], impulse),
            2e-9,
        );
        vector_near(
            next.linear_momenta_kg_m_s()[1],
            sub(state.linear_momenta_kg_m_s()[1], impulse),
            2e-9,
        );
        for node in 0..2 {
            vector_near(
                next.body_spins_kg_m2_s()[node],
                [0.0, 0.0, 0.4 * tangent],
                2e-9,
            );
        }
        near(
            report.contact.dissipated_energy,
            2.0 * (1.0 - 0.4 * 0.4) / inverse_mass
                - tangent
                - 0.5 * tangent_mobility * tangent * tangent,
            2e-9,
        );
        assert!(bodies.is_empty() && report.transfer.rigid_impulses.is_empty());
        assert_eq!(next.pose(), state.pose());
        audit(&report, &[witness]);
    }
}

#[test]
fn original_supports_and_redundant_graph_edges_retain_their_physical_authority() {
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
    let bodies = [
        participant(1.3, [0.8, 0.07, 0.0], [1.0, 0.3, -0.2]),
        participant(2.1, [0.8, 0.07, 0.0], [2.0, 0.4, 0.1]),
    ];
    let witnesses = [
        contact(
            rod_point(1, [0.0, 0.07, 0.0]),
            rigid_point(0, [0.0; 3]),
            0.0,
            1.0,
        ),
        contact(rigid_point(0, [0.0; 3]), rigid_point(1, [0.0; 3]), 0.0, 1.0),
    ];
    let (next, states, report) = rod
        .apply_rigid_frictional_contact_graph(
            &bodies,
            &witnesses,
            bounds(),
            RodRigidGraphBudget::default(),
            &gate,
        )
        .unwrap();
    let delivered = add(
        array(bodies[0].state.linear_momentum_world()),
        array(bodies[1].state.linear_momentum_world()),
    );
    let support = &report
        .transfer
        .rod_impulse
        .as_ref()
        .unwrap()
        .support_impulses[1];
    vector_near(support.linear_n_s, delivered.map(|p| -p), 2e-9);
    vector_near(
        support.angular_nm_s,
        cross([0.0, 0.07, 0.0], delivered).map(|p| -p),
        2e-9,
    );
    assert_eq!(momentum_bits(&next), momentum_bits(&rod));
    for state in states {
        vector_near(array(state.linear_momentum_world()), [0.0; 3], 2e-9);
        vector_near(array(state.angular_momentum_body()), [0.0; 3], 2e-9);
    }
    near(
        report.contact.dissipated_energy,
        bodies
            .iter()
            .map(|p| p.state.kinetic_energy(p.inertia).unwrap())
            .sum(),
        2e-9,
    );
    audit(&report, &witnesses);

    let rod = free_rod([-0.6, 0.8, -0.6]);
    let bodies = [participant(2.3, [0.8, 0.0, 0.0], [1.4, 0.0, 0.0])];
    let edge = contact(rod_point(1, [0.0; 3]), rigid_point(0, [0.0; 3]), 0.35, 0.1);
    let one = rod
        .apply_rigid_frictional_contact_graph(
            &bodies,
            &[edge],
            bounds(),
            RodRigidGraphBudget::default(),
            &gate,
        )
        .unwrap();
    let isolated = rod
        .apply_rigid_frictional_impact(
            bodies[0].state,
            bodies[0].inertia,
            fs_couple::rod_motion::rigid_impact::friction::RodRigidFrictionalContact {
                geometry: witness(1, 0.35),
                tangent_world: edge.tangent_world,
                friction_coefficient: edge.friction_coefficient,
            },
            bounds(),
            &gate,
        )
        .unwrap();
    vector_near(
        one.2.transfer.impulses_world_n_s[0],
        isolated.2.impulse_world_n_s,
        2e-9,
    );
    let duplicate = [edge; 2];
    let repeated = rod
        .apply_rigid_frictional_contact_graph(
            &bodies,
            &duplicate,
            bounds(),
            RodRigidGraphBudget::default(),
            &gate,
        )
        .unwrap();
    let FrictionalContactDeterminacy::UniqueVelocity(ref range) = repeated.2.contact.determinacy
    else {
        panic!("duplicate contacts do not determine local impulses");
    };
    assert_eq!((range.dimension, range.rank), (6, 3));
    assert!(repeated.2.contact.impulse_work_witness.is_some());
    vector_near(
        repeated.0.linear_momenta_kg_m_s()[1],
        one.0.linear_momenta_kg_m_s()[1],
        2e-9,
    );
    vector_near(
        array(repeated.1[0].linear_momentum_world()),
        array(one.1[0].linear_momentum_world()),
        2e-9,
    );
    near(
        repeated.2.contact.dissipated_energy,
        one.2.contact.dissipated_energy,
        2e-9,
    );
    audit(&repeated.2, &duplicate);

    let mut incompatible = duplicate;
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
    let refused = rod
        .apply_rigid_frictional_contact_graph(
            &bodies,
            &incompatible,
            bounds(),
            RodRigidGraphBudget::default(),
            &gate,
        )
        .unwrap_err();
    let RodRigidFrictionalError::Contact(refused) = refused else {
        panic!("preserve original law refusal");
    };
    assert_eq!(
        std::mem::discriminant(&error),
        std::mem::discriminant(&refused)
    );
}

#[test]
fn budgets_geometry_and_cancellation_refuse_atomically_and_zero_graphs_preserve_bits() {
    let gate = CancelGate::new();
    let rod = free_rod([-0.6, 0.8, -0.6]);
    let bodies = [participant(2.3, [0.8, 0.0, 0.0], [1.4, 0.0, 0.0])];
    let edge = contact(rod_point(1, [0.0; 3]), rigid_point(0, [0.0; 3]), 0.35, 0.1);
    let original_bits = momentum_bits(&rod);
    let too_many = [edge; 22];
    assert!(
        rod.apply_rigid_frictional_contact_graph(
            &bodies,
            &too_many,
            bounds(),
            RodRigidGraphBudget::default(),
            &gate,
        )
        .is_err()
    );
    assert!(
        rod.apply_rigid_frictional_contact_graph(
            &bodies,
            &[edge],
            bounds(),
            RodRigidGraphBudget { maximum_bodies: 0 },
            &gate,
        )
        .is_err()
    );
    assert!(
        rod.apply_rigid_frictional_contact_graph(
            &[],
            &[],
            bounds(),
            RodRigidGraphBudget {
                maximum_bodies: 129
            },
            &gate,
        )
        .is_err()
    );
    let mut short = bounds();
    short.contact.matrix.max_variables = 2;
    assert!(
        rod.apply_rigid_frictional_contact_graph(
            &bodies,
            &[edge],
            short,
            RodRigidGraphBudget::default(),
            &gate,
        )
        .is_err()
    );
    short = bounds();
    short.contact.maximum_contacts = 1;
    assert!(
        rod.apply_rigid_frictional_contact_graph(
            &bodies,
            &[edge; 2],
            short,
            RodRigidGraphBudget::default(),
            &gate,
        )
        .is_err()
    );
    short = bounds();
    short.contact.maximum_evaluations = 1;
    assert!(
        rod.apply_rigid_frictional_contact_graph(
            &bodies,
            &[edge],
            short,
            RodRigidGraphBudget::default(),
            &gate,
        )
        .is_err()
    );
    let mut absent = edge;
    absent.geometry.negative = rigid_point(1, [0.0; 3]);
    assert!(
        rod.apply_rigid_frictional_contact_graph(
            &bodies,
            &[edge, absent],
            bounds(),
            RodRigidGraphBudget::default(),
            &gate,
        )
        .is_err()
    );
    for same in [rigid_point(0, [0.0; 3]), rod_point(1, [0.0; 3])] {
        let invalid = contact(same, same, 0.0, 0.1);
        assert!(
            rod.apply_rigid_frictional_contact_graph(
                &bodies,
                &[invalid],
                bounds(),
                RodRigidGraphBudget::default(),
                &gate,
            )
            .is_err()
        );
    }
    let mut bad_frame = edge;
    bad_frame.tangent_world = edge.geometry.normal_world;
    assert!(
        rod.apply_rigid_frictional_contact_graph(
            &bodies,
            &[bad_frame],
            bounds(),
            RodRigidGraphBudget::default(),
            &gate,
        )
        .is_err()
    );
    let mut bad_law = edge;
    bad_law.friction_coefficient = -0.1;
    assert!(matches!(
        rod.apply_rigid_frictional_contact_graph(
            &bodies,
            &[bad_law],
            bounds(),
            RodRigidGraphBudget::default(),
            &gate,
        ),
        Err(RodRigidFrictionalError::Contact(_))
    ));
    let cancelled = CancelGate::new();
    cancelled.request();
    assert!(
        rod.apply_rigid_frictional_contact_graph(
            &bodies,
            &[edge],
            bounds(),
            RodRigidGraphBudget::default(),
            &cancelled,
        )
        .is_err()
    );
    for position in [[0.8 + 1e-4, 0.0, 0.0], [0.8 - 1e-9, 0.0, 0.0]] {
        let separated = [participant(2.3, position, [1.4, 0.0, 0.0])];
        assert!(
            rod.apply_rigid_frictional_contact_graph(
                &separated,
                &[edge],
                bounds(),
                RodRigidGraphBudget::default(),
                &gate,
            )
            .is_err()
        );
    }
    let crooked = [participant(2.3, [0.8, 1e-4, 0.0], [1.4, 0.0, 0.0])];
    let mut loose = bounds();
    loose.impact.geometry_tolerance_m = 1e-3;
    assert!(matches!(
        rod.apply_rigid_frictional_contact_graph(
            &crooked,
            &[edge],
            loose,
            RodRigidGraphBudget::default(),
            &gate,
        ),
        Err(RodRigidFrictionalError::Impact(
            RodRigidImpactError::Balance { .. }
        ))
    ));
    assert_eq!(original_bits, momentum_bits(&rod));
    let first = rod
        .apply_rigid_frictional_contact_graph(
            &bodies,
            &[edge],
            bounds(),
            RodRigidGraphBudget::default(),
            &gate,
        )
        .unwrap();
    let retry = rod
        .apply_rigid_frictional_contact_graph(
            &bodies,
            &[edge],
            bounds(),
            RodRigidGraphBudget::default(),
            &gate,
        )
        .unwrap();
    assert_eq!(momentum_bits(&first.0), momentum_bits(&retry.0));
    assert_eq!(first.1, retry.1);
    for quiet in [
        participant(2.3, [0.79, 0.0, 0.0], [1.4, 0.0, 0.0]),
        participant(2.3, [0.8, 0.0, 0.0], [-1.4, 0.0, 0.0]),
    ] {
        let (next, states, report) = rod
            .apply_rigid_frictional_contact_graph(
                &[quiet],
                &[edge],
                bounds(),
                RodRigidGraphBudget::default(),
                &gate,
            )
            .unwrap();
        assert_eq!(momentum_bits(&next), original_bits);
        assert_eq!(states, vec![quiet.state]);
        assert!(report.transfer.rod_impulse.is_none());
        assert!(report.transfer.rigid_impulses.iter().all(Option::is_none));
        assert!(report.transfer.balance.is_none());
    }
    for participants in [&bodies[..], &[][..]] {
        let (next, states, report) = rod
            .apply_rigid_frictional_contact_graph(
                participants,
                &[],
                bounds(),
                RodRigidGraphBudget {
                    maximum_bodies: participants.len(),
                },
                &gate,
            )
            .unwrap();
        assert_eq!(momentum_bits(&next), original_bits);
        assert_eq!(
            states,
            participants
                .iter()
                .map(|body| body.state)
                .collect::<Vec<_>>()
        );
        assert!(report.contact.impulses.is_empty());
        assert!(report.transfer.balance.is_none());
    }
}
