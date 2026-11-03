//! G0/G1/G3/G4: original shared mass, independent dense contact equations.
use super::*;
use fs_contact::multicontact::{MultiImpactConfig, MultiImpactDeterminacy};
use fs_couple::plate_impact::multiple::{
    PlateImpactSite, PlateRigidImpactSet, PlateRigidImpactSetBudget, PlateRigidImpactSetControls,
};

const POINTS: [[f64; 2]; 2] = [[1.25, -0.65], [2.1, -0.65]];

fn set_controls() -> PlateRigidImpactSetControls {
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
fn set_budget() -> PlateRigidImpactSetBudget {
    PlateRigidImpactSetBudget {
        plate: budget(),
        max_contacts: 8,
        max_mobility_products: 1000,
    }
}
fn sites(points: &[[f64; 2]], e: f64) -> Vec<PlateImpactSite> {
    points
        .iter()
        .map(|&point_xy| PlateImpactSite {
            point_xy,
            activation_tolerance_m: 1e-9,
            restitution: e,
        })
        .collect()
}
// Independent P1 transpose on the manufactured triangle, including full rows.
fn row(point: [f64; 2]) -> [f64; 9] {
    let [x, y] = [point[0] - 1.2, point[1] + 0.7];
    [1.0 - x - y, 0.0, 0.0, x, 0.0, 0.0, y, 0.0, 0.0]
}
fn body_at(
    points: &[[f64; 2]],
    z: &[f64],
    mass: MassProperties,
    orientation: UnitQuaternion,
    velocity: Vec3,
    omega_world: Vec3,
) -> (RigidBodyState, Vec<Vec3>) {
    let position = Vec3::new(1.675, -0.65, 0.35);
    let body = RigidBodyState::new(
        Pose::new(position, orientation).unwrap(),
        velocity.scale(mass.mass()),
        mass.angular_momentum_body_checked(orientation.rotate_world_to_body(omega_world))
            .unwrap(),
    )
    .unwrap();
    let arms = points
        .iter()
        .zip(z)
        .map(|(p, &z)| orientation.rotate_world_to_body(Vec3::new(p[0], p[1], z).sub(position)))
        .collect();
    (body, arms)
}
fn tensor_mass() -> MassProperties {
    MassProperties::from_inertia_tensor(4.3, Vec3::ZERO, Mat3 { m: TENSOR }).unwrap()
}

#[test]
fn coupled_two_point_impact_matches_dense_equations_and_permutation() {
    let (mesh, model, full) = fixture(false);
    let dynamics = dynamics(&model, 1.3, 0.23);
    let mass = tensor_mass();
    let orientation = UnitQuaternion::new(0.94, 0.2, -0.1, 0.25).unwrap();
    for e in [0.0, 0.4, 1.0] {
        let mut reference: Option<(Vec<f64>, Vec3, Vec3)> = None;
        for order in [[0, 1], [1, 0]] {
            let points = order.map(|i| POINTS[i]);
            let prepared = PlateRigidImpactSet::prepare(
                dynamics,
                &mesh,
                &full,
                &sites(&points, e),
                set_controls(),
                set_budget(),
                &mut || false,
            )
            .unwrap();
            let velocity: Vec<_> = (0..9).map(|i| 0.02 * (i as f64 - 4.0)).collect();
            let mut plate = SecondOrderState::new(
                0.7,
                &[0.0; 9],
                &velocity,
                &velocity.iter().map(|v| -0.23 * v).collect::<Vec<_>>(),
            );
            let method =
                OperatorGeneralizedAlpha::new(9, 1e-3, 1.0, ImplicitSolveConfig::default());
            method.step(&mut plate, &dynamics, &[0.0; 9]).unwrap();
            let rows = points.map(row);
            let heights: Vec<_> = rows
                .iter()
                .map(|b| b.iter().zip(&plate.q).map(|(b, q)| b * q).sum::<f64>() + 1e-13)
                .collect();
            let (mut rigid, arms) = body_at(
                &points,
                &heights,
                mass,
                orientation,
                Vec3::new(-0.1, 0.05, -1.0),
                Vec3::new(0.08, -0.12, 0.03),
            );
            let before = plate.clone();
            let body_before = rigid;
            let dense: Vec<_> = full.to_dense().iter().map(|v| 1.3 * v).collect();
            let factor = lu(&dense, 9).unwrap();
            let mut responses = rows;
            for response in &mut responses {
                factor.solve(response);
            }
            let normal_body = orientation.rotate_world_to_body(Vec3::new(0.0, 0.0, 1.0));
            let torques: Vec<_> = arms.iter().map(|arm| arm.cross(normal_body)).collect();
            let mut spin_responses: Vec<_> = torques.iter().map(|t| [t.x, t.y, t.z]).collect();
            let tensor = lu(&TENSOR, 3).unwrap();
            for response in &mut spin_responses {
                tensor.solve(response);
            }
            let mut w = [0.0; 4];
            let mut speed = [0.0; 2];
            for i in 0..2 {
                speed[i] = rigid
                    .point_kinematics(mass, arms[i])
                    .unwrap()
                    .point_velocity_world
                    .z
                    - rows[i]
                        .iter()
                        .zip(&plate.v)
                        .map(|(b, v)| b * v)
                        .sum::<f64>();
                for j in 0..2 {
                    let plate_mobility: f64 =
                        rows[i].iter().zip(responses[j]).map(|(b, r)| b * r).sum();
                    close(prepared.plate_cross_mobility()[2 * i + j], plate_mobility);
                    let r = spin_responses[j];
                    w[2 * i + j] = plate_mobility
                        + 1.0 / mass.mass()
                        + torques[i].dot(Vec3::new(r[0], r[1], r[2]));
                }
            }
            let mut expected = speed.map(|v| -(1.0 + e) * v);
            lu(&w, 2).unwrap().solve(&mut expected);
            assert!(expected.iter().all(|j| *j > 0.0));
            let independent = [-(1.0 + e) * speed[0] / w[0], -(1.0 + e) * speed[1] / w[3]];
            assert!(
                (independent[0] - expected[0]).abs() > 0.01,
                "this fixture must detect independent scalar impacts"
            );
            let report = prepared
                .apply(&mut plate, &mut rigid, mass, &arms, &mut || false)
                .unwrap();
            assert!(matches!(
                report.contact.determinacy,
                MultiImpactDeterminacy::UniqueImpulses
            ));
            for i in 0..2 {
                close(report.contact.normal_impulses[i], expected[i]);
                close(report.relative_normal_speeds_after[i], -e * speed[i]);
            }
            for r in 0..9 {
                close(
                    plate.v[r],
                    before.v[r] - expected[0] * responses[0][r] - expected[1] * responses[1][r],
                );
                close(plate.a[r], -0.23 * plate.v[r]);
            }
            close_vec(
                rigid.linear_momentum_world(),
                body_before.linear_momentum_world().add(Vec3::new(
                    0.0,
                    0.0,
                    expected[0] + expected[1],
                )),
            );
            close_vec(
                rigid.angular_momentum_body(),
                body_before
                    .angular_momentum_body()
                    .add(torques[0].scale(expected[0]))
                    .add(torques[1].scale(expected[1])),
            );
            let balance = report.balance.unwrap();
            let loss = -0.5 * (1.0 - e) * (expected[0] * speed[0] + expected[1] * speed[1]);
            close(balance.kinetic_energy_change_j, -loss);
            close(balance.energy_residual_j, 0.0);
            close_vec(balance.linear_balance_residual, Vec3::ZERO);
            close_vec(balance.angular_balance_residual, Vec3::ZERO);
            assert_eq!(
                (&plate.q, plate.t, plate.steps, &plate.history),
                (&before.q, before.t, before.steps, &before.history)
            );
            assert_eq!(rigid.pose(), body_before.pose());
            if let Some((v, p, l)) = &reference {
                for (a, b) in plate.v.iter().zip(v) {
                    close(*a, *b);
                }
                close_vec(rigid.linear_momentum_world(), *p);
                close_vec(rigid.angular_momentum_body(), *l);
            } else {
                reference = Some((
                    plate.v.clone(),
                    rigid.linear_momentum_world(),
                    rigid.angular_momentum_body(),
                ));
            }
            let post = plate.clone();
            method.step(&mut plate, &dynamics, &[0.0; 9]).unwrap();
            for i in 0..9 {
                close(
                    plate.v[i],
                    post.v[i] * (1.0 - 0.23e-3 / 2.0) / (1.0 + 0.23e-3 / 2.0),
                );
            }
            let continued = RigidBodyIntegrator::new(Gravity::ZERO)
                .step(rigid, mass, Wrench::ZERO, 1e-3)
                .unwrap();
            close_vec(
                continued.state_after.linear_momentum_world(),
                rigid.linear_momentum_world(),
            );
        }
    }
}

#[test]
fn initially_separating_active_site_receives_the_shared_impulse() {
    let (mesh, model, full) = fixture(false);
    let prepared = PlateRigidImpactSet::prepare(
        dynamics(&model, 1.0, 0.0),
        &mesh,
        &full,
        &sites(&POINTS, 0.0),
        set_controls(),
        set_budget(),
        &mut || false,
    )
    .unwrap();
    let mass = tensor_mass();
    for second_gap in [1e-13, 0.01] {
        let mut plate = SecondOrderState::new(0.0, &[0.0; 9], &[0.0; 9], &[0.0; 9]);
        let (mut rigid, arms) = body_at(
            &POINTS,
            &[1e-13, second_gap],
            mass,
            UnitQuaternion::IDENTITY,
            Vec3::new(0.0, 0.0, -0.45),
            Vec3::new(0.0, -1.1 / 0.85, 0.0),
        );
        close(
            rigid
                .point_kinematics(mass, arms[0])
                .unwrap()
                .point_velocity_world
                .z,
            -1.0,
        );
        close(
            rigid
                .point_kinematics(mass, arms[1])
                .unwrap()
                .point_velocity_world
                .z,
            0.1,
        );
        let report = prepared
            .apply(&mut plate, &mut rigid, mass, &arms, &mut || false)
            .unwrap();
        assert!(report.contact.normal_impulses[0] > 0.0);
        if second_gap < 1e-9 {
            assert!(report.contact.normal_impulses[1] > 0.0);
            for speed in report.relative_normal_speeds_after {
                close(speed, 0.0);
            }
        } else {
            assert_eq!(report.contact.normal_impulses[1], 0.0);
            assert!(
                report.relative_normal_speeds_after[1] < 0.0,
                "inactive row still measures the true coupled cross response"
            );
        }
    }
}

#[test]
fn singleton_and_redundant_rows_preserve_the_same_native_motion_and_total_loss() {
    let (mesh, model, full) = fixture(false);
    let dynamics = dynamics(&model, 1.0, 0.0);
    let point = [1.4, -0.4];
    let mass = tensor_mass();
    let initial = SecondOrderState::new(0.0, &[0.0; 9], &[0.0; 9], &[0.0; 9]);
    let (body, arms) = body_at(
        &[point],
        &[1e-13],
        mass,
        UnitQuaternion::IDENTITY,
        Vec3::new(0.0, 0.0, -1.0),
        Vec3::ZERO,
    );
    let scalar = PlateRigidImpact::prepare(
        dynamics,
        &mesh,
        &full,
        point,
        controls(0.5),
        budget(),
        &mut || false,
    )
    .unwrap();
    let mut reference = initial.clone();
    let mut reference_body = body;
    let one = scalar
        .apply(
            &mut reference,
            &mut reference_body,
            mass,
            arms[0],
            &mut || false,
        )
        .unwrap();
    for count in [1, 3] {
        let prepared = PlateRigidImpactSet::prepare(
            dynamics,
            &mesh,
            &full,
            &sites(&vec![point; count], 0.5),
            set_controls(),
            set_budget(),
            &mut || false,
        )
        .unwrap();
        let mut plate = initial.clone();
        let mut rigid = body;
        let report = prepared
            .apply(
                &mut plate,
                &mut rigid,
                mass,
                &vec![arms[0]; count],
                &mut || false,
            )
            .unwrap();
        if count == 3 {
            let MultiImpactDeterminacy::UniqueVelocity(witness) = report.contact.determinacy else {
                panic!("redundant loads need a velocity-only witness")
            };
            assert_eq!(witness.rank, 1);
        }
        close(
            report.contact.normal_impulses.iter().sum(),
            one.contact.normal_impulse,
        );
        close(
            report.contact.dissipated_energy,
            one.contact.dissipated_energy,
        );
        for (a, b) in plate.v.iter().zip(&reference.v) {
            close(*a, *b);
        }
        close_vec(
            rigid.linear_momentum_world(),
            reference_body.linear_momentum_world(),
        );
        close_vec(
            rigid.angular_momentum_body(),
            reference_body.angular_momentum_body(),
        );
    }
}

#[test]
fn simultaneous_supported_points_recover_complete_cross_inertia_reactions() {
    let (mesh, model, full) = fixture(true);
    let points = [[1.2, -0.7], POINTS[1]];
    let prepared = PlateRigidImpactSet::prepare(
        dynamics(&model, 1.7, 0.0),
        &mesh,
        &full,
        &sites(&points, 0.3),
        set_controls(),
        set_budget(),
        &mut || false,
    )
    .unwrap();
    let mass = tensor_mass();
    let mut plate = SecondOrderState::new(0.0, &[0.0; 6], &[0.0; 6], &[0.0; 6]);
    let (mut rigid, arms) = body_at(
        &points,
        &[1e-13; 2],
        mass,
        UnitQuaternion::IDENTITY,
        Vec3::new(0.0, 0.0, -1.0),
        Vec3::ZERO,
    );
    let report = prepared
        .apply(&mut plate, &mut rigid, mass, &arms, &mut || false)
        .unwrap();
    let mut full_change = vec![0.0; 9];
    for (i, r) in model.dof_map.iter().enumerate() {
        if let Some(r) = r {
            full_change[i] = plate.v[*r];
        }
    }
    let dense = full.to_dense();
    let rows = points.map(row);
    let balance = report.balance.unwrap();
    let mut missing_inertia = 0.0_f64;
    for &(i, reaction) in &balance.support_impulses {
        let inertial: f64 = (0..9)
            .map(|j| 1.7 * dense[9 * i + j] * full_change[j])
            .sum();
        let applied = -rows[0][i] * report.contact.normal_impulses[0]
            - rows[1][i] * report.contact.normal_impulses[1];
        close(reaction, inertial - applied);
        missing_inertia = missing_inertia.max((reaction + applied).abs());
    }
    assert!(
        missing_inertia > 1e-3,
        "fixed-row applied load alone must fail this fixture"
    );
    close_vec(balance.linear_balance_residual, Vec3::ZERO);
    close_vec(balance.angular_balance_residual, Vec3::ZERO);
    close(balance.energy_residual_j, 0.0);
}

#[test]
fn simultaneous_budget_noop_and_late_cancellation_preserve_both_states() {
    let (mesh, model, full) = fixture(false);
    let dynamics = dynamics(&model, 1.0, 0.0);
    let sites = sites(&POINTS, 0.4);
    let prepared = PlateRigidImpactSet::prepare(
        dynamics,
        &mesh,
        &full,
        &sites,
        set_controls(),
        set_budget(),
        &mut || false,
    )
    .unwrap();
    let initial = SecondOrderState::new(3.0, &[0.0; 9], &[0.0; 9], &[0.0; 9]);
    let mass = tensor_mass();
    let (body, arms) = body_at(
        &POINTS,
        &[1e-13; 2],
        mass,
        UnitQuaternion::IDENTITY,
        Vec3::new(0.0, 0.0, -1.0),
        Vec3::ZERO,
    );
    let mut plate = initial.clone();
    let mut rigid = body;
    let mut polls = 0;
    prepared
        .apply(&mut plate, &mut rigid, mass, &arms, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    let accepted = plate.clone();
    let accepted_body = rigid;
    let mut plate = initial.clone();
    let mut rigid = body;
    let mut seen = 0;
    assert!(matches!(
        prepared.apply(&mut plate, &mut rigid, mass, &arms, &mut || {
            seen += 1;
            seen == polls
        }),
        Err(PlateRigidImpactError::Cancelled)
    ));
    assert_eq!(plate, initial);
    assert_eq!(rigid, body);
    prepared
        .apply(&mut plate, &mut rigid, mass, &arms, &mut || false)
        .unwrap();
    assert_eq!(plate, accepted);
    assert_eq!(rigid, accepted_body);
    for (gap, speed) in [(0.1, -1.0), (1e-13, 1.0)] {
        let mut plate = initial.clone();
        let (mut rigid, arms) = body_at(
            &POINTS,
            &[gap; 2],
            mass,
            UnitQuaternion::IDENTITY,
            Vec3::new(0.0, 0.0, speed),
            Vec3::ZERO,
        );
        let before = rigid;
        let report = prepared
            .apply(&mut plate, &mut rigid, mass, &arms, &mut || false)
            .unwrap();
        assert!(report.contact.normal_impulses.iter().all(|&j| j == 0.0));
        assert!(report.plate.is_none() && report.body.is_none() && report.balance.is_none());
        assert_eq!(plate, initial);
        assert_eq!(rigid, before);
    }
    let mut cap = set_budget();
    cap.max_mobility_products = 35;
    assert!(matches!(
        PlateRigidImpactSet::prepare(
            dynamics,
            &mesh,
            &full,
            &sites,
            set_controls(),
            cap,
            &mut || false
        ),
        Err(PlateRigidImpactError::Budget(_))
    ));
    let mut plate = initial.clone();
    let mut rigid = body;
    let mut penetrating = arms.clone();
    penetrating[1].z -= 0.01;
    assert!(matches!(
        prepared.apply(&mut plate, &mut rigid, mass, &penetrating, &mut || false),
        Err(PlateRigidImpactError::MultipleContact(_))
    ));
    assert_eq!(plate, initial);
    assert_eq!(rigid, body);
}
