//! G0/G1/G3/G4: distinct finite bodies coupled through the original plate mass.

use super::*;
use fs_contact::multicontact::{MultiImpactConfig, MultiImpactDeterminacy};
use fs_couple::plate_impact::bodies::{
    PlateImpactBody, PlateImpactBodyWitness, PlateRigidImpactNetwork, PlateRigidImpactNetworkBudget,
};
use fs_couple::plate_impact::multiple::{
    PlateImpactSite, PlateRigidImpactSet, PlateRigidImpactSetBudget, PlateRigidImpactSetControls,
};

fn network_controls() -> PlateRigidImpactSetControls {
    let single = controls(0.0);
    PlateRigidImpactSetControls {
        impulse: single.impulse,
        contact: MultiImpactConfig::default(),
        point_tolerance_m: single.point_tolerance_m,
        speed_tolerance_m_s: single.speed_tolerance_m_s,
        linear_impulse_tolerance_ns: single.linear_impulse_tolerance_ns,
        angular_impulse_tolerance_nms: single.angular_impulse_tolerance_nms,
        energy_tolerance_j: single.energy_tolerance_j,
    }
}
fn network_budget() -> PlateRigidImpactNetworkBudget {
    PlateRigidImpactNetworkBudget {
        contacts: PlateRigidImpactSetBudget {
            plate: budget(),
            max_contacts: 8,
            max_mobility_products: 1000,
        },
        max_bodies: 4,
    }
}
fn sites(points: &[[f64; 2]], restitution: f64) -> Vec<PlateImpactSite> {
    points
        .iter()
        .map(|&point_xy| PlateImpactSite {
            point_xy,
            activation_tolerance_m: 1e-9,
            restitution,
        })
        .collect()
}
// Independent exact P1 transpose on the original unit right triangle.
fn row(point: [f64; 2]) -> [f64; 9] {
    let [x, y] = [point[0] - 1.2, point[1] + 0.7];
    [1.0 - x - y, 0.0, 0.0, x, 0.0, 0.0, y, 0.0, 0.0]
}
fn at_point(point: [f64; 2], height: f64, mass: f64, vz: f64) -> PlateImpactBody {
    let properties = MassProperties::new(mass, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap();
    PlateImpactBody {
        state: RigidBodyState::new(
            Pose::new(
                Vec3::new(point[0], point[1], height),
                UnitQuaternion::IDENTITY,
            )
            .unwrap(),
            Vec3::new(0.0, 0.0, mass * vz),
            Vec3::ZERO,
        )
        .unwrap(),
        properties,
    }
}
fn identity_witnesses(count: usize) -> Vec<PlateImpactBodyWitness> {
    (0..count)
        .map(|body| PlateImpactBodyWitness {
            body,
            arm_body: Vec3::ZERO,
        })
        .collect()
}

#[test]
fn mixed_body_sites_match_independent_dense_mass_and_tensor_response() {
    let points = [[1.25, -0.65], [1.95, -0.6], [1.3, 0.05]];
    let owners = [0, 1, 0];
    for restitution in [0.0, 0.4, 1.0] {
        let mut reference: Option<(Vec<f64>, Vec<RigidBodyState>)> = None;
        for scene in 0..2 {
            let contact_order = if scene == 0 { [0, 1, 2] } else { [2, 0, 1] };
            let body_order = if scene == 0 { [0, 1, 2] } else { [2, 0, 1] };
            let shift = if scene == 0 { [0.0, 0.0] } else { [0.7, -0.4] };
            let (mut mesh, model, full) = fixture(false);
            for (x, y) in &mut mesh.nodes {
                *x += shift[0];
                *y += shift[1];
            }
            let dynamics = dynamics(&model, 1.3, 0.23);
            let ordered_points =
                contact_order.map(|i| [points[i][0] + shift[0], points[i][1] + shift[1]]);
            let rows = contact_order.map(|i| row(points[i]));
            let prepared = PlateRigidImpactNetwork::prepare(
                dynamics,
                &mesh,
                &full,
                &sites(&ordered_points, restitution),
                network_controls(),
                network_budget(),
                &mut || false,
            )
            .unwrap();
            let velocity = [0.02, 0.0, 0.0, -0.04, 0.0, 0.0, 0.01, 0.0, 0.0];
            let mut plate =
                SecondOrderState::new(0.7, &[0.0; 9], &velocity, &velocity.map(|v| -0.23 * v));
            let method =
                OperatorGeneralizedAlpha::new(9, 1e-3, 1.0, ImplicitSolveConfig::default());
            method.step(&mut plate, &dynamics, &[0.0; 9]).unwrap();
            let tensors = [TENSOR, TENSOR.map(|v| 0.65 * v)];
            let mut original = Vec::new();
            for index in 0..2 {
                let properties = MassProperties::from_inertia_tensor(
                    if index == 0 { 4.3 } else { 2.2 },
                    Vec3::ZERO,
                    Mat3 { m: tensors[index] },
                )
                .unwrap();
                let orientation = if index == 0 {
                    UnitQuaternion::new(0.94, 0.2, -0.1, 0.25).unwrap()
                } else {
                    UnitQuaternion::new(0.91, -0.12, 0.23, 0.17).unwrap()
                };
                let position = if index == 0 {
                    [1.3, -0.25, 0.35]
                } else {
                    [1.95, -0.6, 0.2]
                };
                let omega = if index == 0 {
                    Vec3::new(0.04, -0.06, 0.02)
                } else {
                    Vec3::new(-0.05, 0.03, 0.01)
                };
                original.push(PlateImpactBody {
                    state: RigidBodyState::new(
                        Pose::new(
                            Vec3::new(position[0] + shift[0], position[1] + shift[1], position[2]),
                            orientation,
                        )
                        .unwrap(),
                        Vec3::new(0.08, -0.03, if index == 0 { -0.9 } else { -1.2 })
                            .scale(properties.mass()),
                        properties
                            .angular_momentum_body_checked(orientation.rotate_world_to_body(omega))
                            .unwrap(),
                    )
                    .unwrap(),
                    properties,
                });
            }
            original.push(at_point([8.0, -4.0], 7.0, 2.0, 3.0));
            let mut bodies: Vec<_> = body_order.iter().map(|&i| original[i]).collect();
            let witnesses: Vec<_> = contact_order
                .iter()
                .enumerate()
                .map(|(i, &source)| {
                    let owner = owners[source];
                    let body = body_order.iter().position(|&b| b == owner).unwrap();
                    let participant = bodies[body];
                    let height = rows[i]
                        .iter()
                        .zip(&plate.q)
                        .map(|(b, q)| b * q)
                        .sum::<f64>()
                        + 1e-13;
                    PlateImpactBodyWitness {
                        body,
                        arm_body: participant.state.pose().orientation().rotate_world_to_body(
                            Vec3::new(ordered_points[i][0], ordered_points[i][1], height)
                                .sub(participant.state.pose().position_world()),
                        ),
                    }
                })
                .collect();
            let before = plate.clone();
            let before_bodies = bodies.clone();
            let dense: Vec<_> = full.to_dense().iter().map(|v| 1.3 * v).collect();
            let factor = lu(&dense, 9).unwrap();
            let mut response = rows;
            for r in &mut response {
                factor.solve(r);
            }
            let mut torque = [Vec3::ZERO; 3];
            let mut spin_response = [[0.0; 3]; 3];
            let mut speed = [0.0; 3];
            for i in 0..3 {
                let participant = bodies[witnesses[i].body];
                let normal = participant
                    .state
                    .pose()
                    .orientation()
                    .rotate_world_to_body(Vec3::new(0.0, 0.0, 1.0));
                torque[i] = witnesses[i].arm_body.cross(normal);
                spin_response[i] = [torque[i].x, torque[i].y, torque[i].z];
                lu(&tensors[owners[contact_order[i]]], 3)
                    .unwrap()
                    .solve(&mut spin_response[i]);
                speed[i] = participant
                    .state
                    .point_kinematics(participant.properties, witnesses[i].arm_body)
                    .unwrap()
                    .point_velocity_world
                    .z
                    - rows[i]
                        .iter()
                        .zip(&plate.v)
                        .map(|(b, v)| b * v)
                        .sum::<f64>();
            }
            let mut matrix = [0.0; 9];
            for i in 0..3 {
                for j in 0..3 {
                    let plate_response: f64 =
                        rows[i].iter().zip(response[j]).map(|(b, r)| b * r).sum();
                    close(prepared.plate_cross_mobility()[3 * i + j], plate_response);
                    matrix[3 * i + j] = plate_response;
                    if witnesses[i].body == witnesses[j].body {
                        matrix[3 * i + j] += 1.0 / bodies[witnesses[i].body].properties.mass()
                            + torque[i].dot(Vec3::new(
                                spin_response[j][0],
                                spin_response[j][1],
                                spin_response[j][2],
                            ));
                    }
                }
            }
            let mut expected = speed.map(|v| -(1.0 + restitution) * v);
            lu(&matrix, 3).unwrap().solve(&mut expected);
            assert!(expected.iter().all(|j| *j > 0.0));
            let report = prepared
                .apply(&mut plate, &mut bodies, &witnesses, &mut || false)
                .unwrap();
            assert!(matches!(
                report.contact.determinacy,
                MultiImpactDeterminacy::UniqueImpulses
            ));
            for i in 0..3 {
                close(report.contact.normal_impulses[i], expected[i]);
                close(
                    report.relative_normal_speeds_after[i],
                    -restitution * speed[i],
                );
                for j in 0..3 {
                    close(report.delassus_kg_inverse[3 * i + j], matrix[3 * i + j]);
                }
            }
            for r in 0..9 {
                close(
                    plate.v[r],
                    before.v[r] - (0..3).map(|i| expected[i] * response[i][r]).sum::<f64>(),
                );
                close(plate.a[r], -0.23 * plate.v[r]);
            }
            for b in 0..bodies.len() {
                let mut j = 0.0;
                let mut l = Vec3::ZERO;
                for i in 0..3 {
                    if witnesses[i].body == b {
                        j += expected[i];
                        l = l.add(torque[i].scale(expected[i]));
                    }
                }
                close_vec(
                    bodies[b].state.linear_momentum_world(),
                    before_bodies[b]
                        .state
                        .linear_momentum_world()
                        .add(Vec3::new(0.0, 0.0, j)),
                );
                close_vec(
                    bodies[b].state.angular_momentum_body(),
                    before_bodies[b].state.angular_momentum_body().add(l),
                );
                assert_eq!(bodies[b].state.pose(), before_bodies[b].state.pose());
                if body_order[b] == 2 {
                    assert_eq!(bodies[b], before_bodies[b]);
                    assert!(report.bodies[b].is_none());
                }
            }
            let balance = report.balance.unwrap();
            close(
                balance.kinetic_energy_change_j,
                0.5 * (1.0 - restitution) * (0..3).map(|i| expected[i] * speed[i]).sum::<f64>(),
            );
            close(balance.energy_residual_j, 0.0);
            close(balance.work_energy_residual_j, 0.0);
            close_vec(balance.linear_balance_residual, Vec3::ZERO);
            close_vec(balance.angular_balance_residual, Vec3::ZERO);
            assert_eq!(
                (&plate.q, plate.t, plate.steps, &plate.history),
                (&before.q, before.t, before.steps, &before.history)
            );
            let states: Vec<_> = (0..3)
                .map(|original| {
                    bodies[body_order.iter().position(|&b| b == original).unwrap()].state
                })
                .collect();
            if let Some((v, prior)) = &reference {
                for (a, b) in plate.v.iter().zip(v) {
                    close(*a, *b);
                }
                for i in 0..3 {
                    close_vec(
                        states[i].linear_momentum_world(),
                        prior[i].linear_momentum_world(),
                    );
                    close_vec(
                        states[i].angular_momentum_body(),
                        prior[i].angular_momentum_body(),
                    );
                }
            } else {
                reference = Some((plate.v.clone(), states));
            }
            if scene == 0 && restitution == 0.0 {
                let mut sequential_plate = before.clone();
                let mut sequential_bodies = before_bodies.clone();
                for i in 0..3 {
                    let single = PlateRigidImpact::prepare(
                        dynamics,
                        &mesh,
                        &full,
                        ordered_points[i],
                        controls(0.0),
                        budget(),
                        &mut || false,
                    )
                    .unwrap();
                    let participant = &mut sequential_bodies[witnesses[i].body];
                    single
                        .apply(
                            &mut sequential_plate,
                            &mut participant.state,
                            participant.properties,
                            witnesses[i].arm_body,
                            &mut || false,
                        )
                        .unwrap();
                }
                assert!(
                    plate
                        .v
                        .iter()
                        .zip(&sequential_plate.v)
                        .any(|(a, b)| (a - b).abs() > 1e-3)
                );
            }
            let post = plate.clone();
            method.step(&mut plate, &dynamics, &[0.0; 9]).unwrap();
            for i in 0..9 {
                close(
                    plate.v[i],
                    post.v[i] * (1.0 - 0.23e-3 / 2.0) / (1.0 + 0.23e-3 / 2.0),
                );
            }
            for participant in &bodies {
                let continued = RigidBodyIntegrator::new(Gravity::ZERO)
                    .step(
                        participant.state,
                        participant.properties,
                        Wrench::ZERO,
                        1e-3,
                    )
                    .unwrap();
                close_vec(
                    continued.state_after.linear_momentum_world(),
                    participant.state.linear_momentum_world(),
                );
            }
        }
    }
}

#[test]
fn distinct_bodies_share_plate_response_including_initially_separating_sites() {
    let (mesh, model, full) = fixture(false);
    let points = [[1.2, -0.7], [2.2, -0.7]];
    let prepared = PlateRigidImpactNetwork::prepare(
        dynamics(&model, 1.0, 0.0),
        &mesh,
        &full,
        &sites(&points, 0.0),
        network_controls(),
        network_budget(),
        &mut || false,
    )
    .unwrap();
    for gap in [1e-13, 0.01] {
        let mut plate = SecondOrderState::new(0.0, &[0.0; 9], &[0.0; 9], &[0.0; 9]);
        let mut bodies = [
            at_point(points[0], 1e-13, 1.0, -1.0),
            at_point(points[1], gap, 1.0, 0.1),
            at_point([8.0, 3.0], 4.0, 2.0, -7.0),
        ];
        let before = bodies;
        let report = prepared
            .apply(&mut plate, &mut bodies, &identity_witnesses(2), &mut || {
                false
            })
            .unwrap();
        for (a, b) in report
            .delassus_kg_inverse
            .iter()
            .zip([4.0, -1.0, -1.0, 4.0])
        {
            close(*a, b);
        }
        if gap < 1e-9 {
            close(report.contact.normal_impulses[0], 0.26);
            close(report.contact.normal_impulses[1], 0.04);
            for v in report.relative_normal_speeds_after {
                close(v, 0.0);
            }
        } else {
            close(report.contact.normal_impulses[0], 0.25);
            assert_eq!(report.contact.normal_impulses[1], 0.0);
            close(report.relative_normal_speeds_after[1], -0.15);
            assert_eq!(bodies[1], before[1]);
            assert!(report.bodies[1].is_none());
        }
        assert_eq!(bodies[2], before[2]);
        assert!(report.bodies[2].is_none());
    }
}

#[test]
fn one_body_network_retains_the_published_simultaneous_api_response() {
    let (mesh, model, full) = fixture(false);
    let points = [[1.25, -0.65], [2.1, -0.65]];
    let dynamics = dynamics(&model, 1.3, 0.23);
    let sites = sites(&points, 0.4);
    let old = PlateRigidImpactSet::prepare(
        dynamics,
        &mesh,
        &full,
        &sites,
        network_controls(),
        network_budget().contacts,
        &mut || false,
    )
    .unwrap();
    let new = PlateRigidImpactNetwork::prepare(
        dynamics,
        &mesh,
        &full,
        &sites,
        network_controls(),
        network_budget(),
        &mut || false,
    )
    .unwrap();
    let mass = MassProperties::from_inertia_tensor(4.3, Vec3::ZERO, Mat3 { m: TENSOR }).unwrap();
    let mut rigid = RigidBodyState::new(
        Pose::new(Vec3::new(1.675, -0.65, 0.3), UnitQuaternion::IDENTITY).unwrap(),
        Vec3::new(0.0, 0.0, -mass.mass()),
        Vec3::ZERO,
    )
    .unwrap();
    let arms: Vec<_> = points
        .iter()
        .map(|p| Vec3::new(p[0], p[1], 1e-13).sub(rigid.pose().position_world()))
        .collect();
    let mut plate = SecondOrderState::new(0.0, &[0.0; 9], &[0.0; 9], &[0.0; 9]);
    let mut network_plate = plate.clone();
    let mut bodies = [PlateImpactBody {
        state: rigid,
        properties: mass,
    }];
    let expected = old
        .apply(&mut plate, &mut rigid, mass, &arms, &mut || false)
        .unwrap();
    let witnesses: Vec<_> = arms
        .iter()
        .map(|&arm_body| PlateImpactBodyWitness { body: 0, arm_body })
        .collect();
    let actual = new
        .apply(&mut network_plate, &mut bodies, &witnesses, &mut || false)
        .unwrap();
    for (a, b) in actual
        .contact
        .normal_impulses
        .iter()
        .zip(&expected.contact.normal_impulses)
    {
        close(*a, *b);
    }
    for (a, b) in network_plate.v.iter().zip(&plate.v) {
        close(*a, *b);
    }
    close_vec(
        bodies[0].state.linear_momentum_world(),
        rigid.linear_momentum_world(),
    );
    close_vec(
        bodies[0].state.angular_momentum_body(),
        rigid.angular_momentum_body(),
    );
    close(
        actual.balance.unwrap().kinetic_energy_change_j,
        expected.balance.unwrap().kinetic_energy_change_j,
    );
}

#[test]
fn supported_multi_body_impacts_preserve_full_reactions_and_redundant_load_determinacy() {
    let (mesh, model, full) = fixture(true);
    let points = [[1.2, -0.7], [2.1, -0.65]];
    let dynamics = dynamics(&model, 1.7, 0.0);
    let initial = [
        at_point(points[0], 1e-13, 2.0, -1.0),
        at_point(points[1], 1e-13, 3.0, -1.0),
    ];
    let mut reference: Option<(Vec<f64>, [PlateImpactBody; 2], f64)> = None;
    for duplicated in [false, true] {
        let points = if duplicated {
            vec![points[0], points[0], points[1]]
        } else {
            points.to_vec()
        };
        let owners = if duplicated {
            vec![0, 0, 1]
        } else {
            vec![0, 1]
        };
        let prepared = PlateRigidImpactNetwork::prepare(
            dynamics,
            &mesh,
            &full,
            &sites(&points, 0.3),
            network_controls(),
            network_budget(),
            &mut || false,
        )
        .unwrap();
        let mut plate = SecondOrderState::new(0.0, &[0.0; 6], &[0.0; 6], &[0.0; 6]);
        let mut bodies = initial;
        let witnesses: Vec<_> = owners
            .iter()
            .map(|&body| PlateImpactBodyWitness {
                body,
                arm_body: Vec3::ZERO,
            })
            .collect();
        let report = prepared
            .apply(&mut plate, &mut bodies, &witnesses, &mut || false)
            .unwrap();
        if duplicated {
            let MultiImpactDeterminacy::UniqueVelocity(witness) = &report.contact.determinacy
            else {
                panic!("duplicate contacts need a velocity witness")
            };
            assert_eq!(witness.rank, 2);
        }
        let balance = report.balance.unwrap();
        let mut full_change = [0.0; 9];
        for (i, r) in model.dof_map.iter().enumerate() {
            if let Some(r) = r {
                full_change[i] = plate.v[*r];
            }
        }
        let dense = full.to_dense();
        let mut inertia_seen: f64 = 0.0;
        for &(i, support) in &balance.support_impulses {
            let inertia: f64 = (0..9)
                .map(|j| 1.7 * dense[9 * i + j] * full_change[j])
                .sum();
            let applied: f64 = points
                .iter()
                .zip(&report.contact.normal_impulses)
                .map(|(&p, &j)| -row(p)[i] * j)
                .sum();
            close(support, inertia - applied);
            inertia_seen = inertia_seen.max((support + applied).abs());
        }
        assert!(inertia_seen > 1e-3);
        close_vec(balance.linear_balance_residual, Vec3::ZERO);
        close_vec(balance.angular_balance_residual, Vec3::ZERO);
        close(balance.energy_residual_j, 0.0);
        if let Some((v, states, loss)) = &reference {
            for (a, b) in plate.v.iter().zip(v) {
                close(*a, *b);
            }
            for i in 0..2 {
                close_vec(
                    bodies[i].state.linear_momentum_world(),
                    states[i].state.linear_momentum_world(),
                );
            }
            close(report.contact.dissipated_energy, *loss);
        } else {
            reference = Some((plate.v.clone(), bodies, report.contact.dissipated_energy));
        }
    }
}

#[test]
fn collective_native_work_retains_small_loss_during_large_common_translation() {
    let (mut mesh, model, full) = fixture(false);
    mesh.nodes = vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)];
    let points = [[0.0, 0.0], [1.0, 0.0]];
    let prepared = PlateRigidImpactNetwork::prepare(
        dynamics(&model, 1.0, 0.0),
        &mesh,
        &full,
        &sites(&points, 0.0),
        network_controls(),
        network_budget(),
        &mut || false,
    )
    .unwrap();
    // The exact shared-plate matrix is [[4,-1],[-1,4]]. A dyadic unit
    // impulse is natively representable on the moving participants while the
    // resulting 3 * unit^2 loss is lost by subtracting their total energies.
    let common = 67_108_864.0;
    let unit = 1.0 / 1024.0;
    let loss = 3.0 * unit * unit;
    let velocity = [common, 0.0, 0.0, common, 0.0, 0.0, common, 0.0, 0.0];
    let mut plate = SecondOrderState::new(0.0, &[0.0; 9], &velocity, &[0.0; 9]);
    let mut bodies = [
        at_point(points[0], 1e-13, 1.0, common - 3.0 * unit),
        at_point(points[1], 1e-13, 1.0, common - 3.0 * unit),
    ];
    let before_bodies = bodies;
    let report = prepared
        .apply(&mut plate, &mut bodies, &identity_witnesses(2), &mut || {
            false
        })
        .unwrap();
    for impulse in report.contact.normal_impulses {
        close(impulse, unit);
    }
    let balance = report.balance.unwrap();
    close(balance.kinetic_energy_change_j, -loss);
    close(report.contact.dissipated_energy, loss);
    close(balance.work_energy_residual_j, 0.0);
    close(balance.energy_residual_j, 0.0);
    close_vec(balance.linear_balance_residual, Vec3::ZERO);
    let energy = |v: &[f64]| {
        let mut mv = vec![0.0; 9];
        full.spmv(v, &mut mv);
        0.5 * v.iter().zip(mv).map(|(v, p)| v * p).sum::<f64>()
    };
    let before = energy(&velocity)
        + before_bodies
            .iter()
            .map(|b| b.state.kinetic_energy(b.properties).unwrap())
            .sum::<f64>();
    let after = energy(&plate.v)
        + bodies
            .iter()
            .map(|b| b.state.kinetic_energy(b.properties).unwrap())
            .sum::<f64>();
    assert!(((after - before) - balance.kinetic_energy_change_j).abs() > 0.5 * loss);

    // The original native work gate remains authoritative. Larger background
    // work magnifies a few ulps in the solved request beyond the SAME 1e-10 J
    // numerical allowance, so the entire network event must be refused.
    let common = 1_073_741_824.0;
    let velocity = [common, 0.0, 0.0, common, 0.0, 0.0, common, 0.0, 0.0];
    let mut plate = SecondOrderState::new(0.0, &[0.0; 9], &velocity, &[0.0; 9]);
    let mut bodies = [
        at_point(points[0], 1e-13, 1.0, common - 3.0),
        at_point(points[1], 1e-13, 1.0, common - 3.0),
    ];
    let before = (plate.clone(), bodies);
    assert!(matches!(
        prepared.apply(&mut plate, &mut bodies, &identity_witnesses(2), &mut || {
            false
        }),
        Err(PlateRigidImpactError::Plate(
            fs_time::galpha::impulse::ImpulseError::Work { .. }
        ))
    ));
    assert_eq!((plate, bodies), before);
}

#[test]
fn participant_limits_noops_late_refusals_and_cancellation_publish_nothing_partial() {
    let (mesh, model, full) = fixture(false);
    let points = [[1.2, -0.7], [2.2, -0.7]];
    let dynamics = dynamics(&model, 1.0, 0.0);
    let sites = sites(&points, 0.4);
    let prepared = PlateRigidImpactNetwork::prepare(
        dynamics,
        &mesh,
        &full,
        &sites,
        network_controls(),
        network_budget(),
        &mut || false,
    )
    .unwrap();
    let initial = SecondOrderState::new(3.0, &[0.0; 9], &[0.0; 9], &[0.0; 9]);
    let initial_bodies = [
        at_point(points[0], 1e-13, 1.0, -1.0),
        at_point(points[1], 1e-13, 1.0, -1.0),
    ];
    let witnesses = identity_witnesses(2);
    let mut plate = initial.clone();
    let mut bodies = initial_bodies;
    let mut polls = 0;
    prepared
        .apply(&mut plate, &mut bodies, &witnesses, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    let expected = (plate.clone(), bodies);
    let mut plate = initial.clone();
    let mut bodies = initial_bodies;
    let mut seen = 0;
    assert!(matches!(
        prepared.apply(&mut plate, &mut bodies, &witnesses, &mut || {
            seen += 1;
            seen == polls
        }),
        Err(PlateRigidImpactError::Cancelled)
    ));
    assert_eq!(plate, initial);
    assert_eq!(bodies, initial_bodies);
    prepared
        .apply(&mut plate, &mut bodies, &witnesses, &mut || false)
        .unwrap();
    assert_eq!((plate, bodies), expected);
    for (gap, velocity) in [(0.1, -1.0), (1e-13, 1.0)] {
        let mut plate = initial.clone();
        let mut bodies = [
            at_point(points[0], gap, 1.0, velocity),
            at_point(points[1], gap, 1.0, velocity),
        ];
        let before = bodies;
        let report = prepared
            .apply(&mut plate, &mut bodies, &witnesses, &mut || false)
            .unwrap();
        assert!(report.plate.is_none() && report.balance.is_none());
        assert!(report.bodies.iter().all(Option::is_none));
        assert_eq!(plate, initial);
        assert_eq!(bodies, before);
    }
    let mut cap = network_budget();
    cap.max_bodies = 1;
    let capped = PlateRigidImpactNetwork::prepare(
        dynamics,
        &mesh,
        &full,
        &sites,
        network_controls(),
        cap,
        &mut || false,
    )
    .unwrap();
    let mut plate = initial.clone();
    let mut bodies = initial_bodies;
    assert!(matches!(
        capped.apply(&mut plate, &mut bodies, &witnesses, &mut || false),
        Err(PlateRigidImpactError::Budget(_))
    ));
    assert_eq!(plate, initial);
    assert_eq!(bodies, initial_bodies);
    cap.max_bodies = usize::MAX;
    assert!(matches!(
        PlateRigidImpactNetwork::prepare(
            dynamics,
            &mesh,
            &full,
            &sites,
            network_controls(),
            cap,
            &mut || false
        ),
        Err(PlateRigidImpactError::Budget(_))
    ));
    let mut invalid = witnesses.clone();
    invalid[1].body = 2;
    assert!(matches!(
        prepared.apply(&mut plate, &mut bodies, &invalid, &mut || false),
        Err(PlateRigidImpactError::Invalid(_))
    ));
    invalid = witnesses.clone();
    invalid[1].arm_body.z = -0.01;
    assert!(matches!(
        prepared.apply(&mut plate, &mut bodies, &invalid, &mut || false),
        Err(PlateRigidImpactError::MultipleContact(_))
    ));
    assert_eq!(plate, initial);
    assert_eq!(bodies, initial_bodies);
    let mut loose = network_controls();
    loose.point_tolerance_m = 1e-3;
    let loose = PlateRigidImpactNetwork::prepare(
        dynamics,
        &mesh,
        &full,
        &sites,
        loose,
        network_budget(),
        &mut || false,
    )
    .unwrap();
    let mut crooked = initial_bodies;
    crooked[1] = at_point([points[1][0] + 1e-4, points[1][1]], 1e-13, 1.0, -1.0);
    let before = crooked;
    assert!(matches!(
        loose.apply(&mut plate, &mut crooked, &witnesses, &mut || false),
        Err(PlateRigidImpactError::Balance { .. })
    ));
    assert_eq!(plate, initial);
    assert_eq!(crooked, before);
}
