use fs_couple::plate_impact::{
    PlateRigidImpact, PlateRigidImpactBudget, PlateRigidImpactControls, PlateRigidImpactError,
};
use fs_ga::Mat3;
use fs_la::factor::lu;
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, UnitQuaternion, Vec3,
    Wrench,
};
use fs_plate::loading::{PlateLoadBudget, PlatePointStencil};
use fs_plate::{PlateMesh, PlateModel};
use fs_plate_transient::{
    NoPlateLoad, PlateDynamics, PlateDynamicsBudget, PlateDynamicsParameters,
};
use fs_sparse::{Coo, Csr};
use fs_time::galpha::impulse::ImpulseControls;
use fs_time::galpha::initialization::InitialSolveConfig;
use fs_time::galpha::{ImplicitSolveConfig, OperatorGeneralizedAlpha, SecondOrderState};

#[path = "plate_impact/multiple.rs"]
mod multiple;

#[path = "plate_impact/bodies.rs"]
mod bodies;

#[path = "plate_impact/plates.rs"]
mod plates;

const TENSOR: [f64; 9] = [2.0, 0.2, -0.1, 0.2, 3.0, 0.3, -0.1, 0.3, 4.0];

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-10 * (1.0 + a.abs().max(b.abs())),
        "{a} != {b}"
    );
}
fn close_vec(a: Vec3, b: Vec3) {
    close(a.x, b.x);
    close(a.y, b.y);
    close(a.z, b.z);
}

// Independent exact P1 area integration on one unit right triangle: total
// transverse mass is 3 kg, M_ii=1/2 and M_ij=1/4. The slope blocks carry the
// corresponding centered rotary inertia. No modal or diagonal mass substitute.
fn fixture(supported: bool) -> (PlateMesh, PlateModel, Csr) {
    let mesh = PlateMesh {
        nodes: vec![(1.2, -0.7), (2.2, -0.7), (1.2, 0.3)],
        tris: vec![[0, 1, 2]],
    };
    let mut mass = Coo::new(9, 9);
    for i in 0..3 {
        for j in 0..3 {
            for component in 0..3 {
                let rotary = if component == 0 { 1.0 } else { 1.0 / 12.0 };
                mass.push(
                    3 * i + component,
                    3 * j + component,
                    rotary * if i == j { 0.5 } else { 0.25 },
                );
            }
        }
    }
    let full = mass.assemble();
    let mut free = 0;
    let map: Vec<_> = (0..9)
        .map(|i| {
            if supported && i < 3 {
                None
            } else {
                let r = free;
                free += 1;
                Some(r)
            }
        })
        .collect();
    let mut reduced = Coo::new(free, free);
    for i in 0..9 {
        if let Some(r) = map[i] {
            let (cols, values) = full.row(i);
            for (&j, &value) in cols.iter().zip(values) {
                if let Some(c) = map[j] {
                    reduced.push(r, c, value);
                }
            }
        }
    }
    let model = PlateModel {
        m: reduced.assemble(),
        k: Coo::new(free, free).assemble(),
        dof_map: map,
        free,
    };
    (mesh, model, full)
}
fn dynamics(model: &PlateModel, scale: f64, damping: f64) -> PlateDynamics<'_, NoPlateLoad> {
    PlateDynamics::new(
        model,
        PlateDynamicsParameters {
            stiffness_scale: 1.0,
            mass_scale: scale,
            mass_damping_per_s: damping,
            stiffness_damping_s: 0.0,
        },
        &NoPlateLoad,
        PlateDynamicsBudget {
            max_dofs: 9,
            max_nonzeros: 100,
            max_full_dofs: 9,
            max_load_parameters: 0,
        },
        &mut || false,
    )
    .unwrap()
}
fn controls(e: f64) -> PlateRigidImpactControls {
    PlateRigidImpactControls {
        impulse: ImpulseControls {
            mass_solve: InitialSolveConfig {
                restart: 9,
                max_cycles: 4,
                tolerance: 1e-13,
            },
            absolute_impulse_tolerance: 1e-11,
            relative_impulse_tolerance: 1e-11,
            work_tolerance: 1e-10,
        },
        restitution: e,
        activation_tolerance_m: 1e-9,
        point_tolerance_m: 1e-12,
        speed_tolerance_m_s: 1e-10,
        linear_impulse_tolerance_ns: 1e-10,
        angular_impulse_tolerance_nms: 1e-10,
        energy_tolerance_j: 1e-9,
    }
}
fn budget() -> PlateRigidImpactBudget {
    PlateRigidImpactBudget {
        location: PlateLoadBudget {
            max_nodes: 3,
            max_triangles: 1,
        },
        max_full_mass_nonzeros: 27,
        max_workspace_components: 10000,
    }
}
fn body(point: [f64; 3], properties: MassProperties) -> (RigidBodyState, Vec3) {
    let orientation = UnitQuaternion::new(0.94, 0.2, -0.1, 0.25).unwrap();
    let arm = Vec3::new(0.4, -0.2, -0.3);
    let position =
        Vec3::new(point[0], point[1], point[2]).sub(orientation.rotate_body_to_world(arm));
    (
        RigidBodyState::new(
            Pose::new(position, orientation).unwrap(),
            Vec3::new(0.1, -0.2, -1.3).scale(properties.mass()),
            Vec3::new(0.2, -0.1, 0.15),
        )
        .unwrap(),
        arm,
    )
}

#[test]
fn eccentric_full_tensor_impact_matches_dense_mobility_and_continues_both_states() {
    for restitution in [0.0, 0.65, 1.0] {
        let (mesh, model, full) = fixture(false);
        let dynamics = dynamics(&model, 1.3, 0.23);
        let velocity: Vec<_> = (0..9).map(|i| 0.03 * (i as f64 - 3.0)).collect();
        let acceleration: Vec<_> = velocity.iter().map(|&v| -0.23 * v).collect();
        let mut plate = SecondOrderState::new(0.7, &[0.0; 9], &velocity, &acceleration);
        let method = OperatorGeneralizedAlpha::new(9, 1e-3, 1.0, ImplicitSolveConfig::default());
        method.step(&mut plate, &dynamics, &[0.0; 9]).unwrap();
        let point = [1.4, -0.4];
        let stencil = PlatePointStencil::locate(&mesh, &model, point, budget().location).unwrap();
        let sample = stencil.sample(&plate.q, &plate.v).unwrap();
        let properties =
            MassProperties::from_inertia_tensor(4.3, Vec3::ZERO, Mat3 { m: TENSOR }).unwrap();
        let (mut rigid, arm) = body(
            [point[0], point[1], sample.displacement + 1e-13],
            properties,
        );
        let before = plate.clone();
        let rigid_before = rigid;
        let prepared = PlateRigidImpact::prepare(
            dynamics,
            &mesh,
            &full,
            point,
            controls(restitution),
            budget(),
            &mut || false,
        )
        .unwrap();
        // Independent dense inverse using the complete matrix and direct tensor
        // solve, not either prepared scalar effective-mass query.
        let mut unit = [0.0; 9];
        stencil.add_load([1.0, 0.0, 0.0], &mut unit).unwrap();
        let dense: Vec<_> = full.to_dense().iter().map(|&m| 1.3 * m).collect();
        let mut response = unit;
        lu(&dense, 9).unwrap().solve(&mut response);
        let wp: f64 = unit.iter().zip(response).map(|(&u, r)| u * r).sum();
        close(prepared.inverse_plate_effective_mass(), wp);
        let nb = rigid
            .pose()
            .orientation()
            .rotate_world_to_body(Vec3::new(0.0, 0.0, 1.0));
        let torque = arm.cross(nb);
        let mut spin = [torque.x, torque.y, torque.z];
        lu(&TENSOR, 3).unwrap().solve(&mut spin);
        let wb = 1.0 / properties.mass() + torque.dot(Vec3::new(spin[0], spin[1], spin[2]));
        let u = rigid
            .point_kinematics(properties, arm)
            .unwrap()
            .point_velocity_world
            .z
            - sample.velocity;
        let j = -(1.0 + restitution) * u / (wp + wb);
        let report = prepared
            .apply(&mut plate, &mut rigid, properties, arm, &mut || false)
            .unwrap();
        close(report.contact.normal_impulse, j);
        close(report.relative_normal_speed_after, -restitution * u);
        for i in 0..9 {
            close(plate.v[i], before.v[i] - j * response[i]);
            close(plate.a[i], -0.23 * plate.v[i]);
        }
        close_vec(
            rigid.linear_momentum_world(),
            rigid_before
                .linear_momentum_world()
                .add(Vec3::new(0.0, 0.0, j)),
        );
        close_vec(
            rigid.angular_momentum_body(),
            rigid_before.angular_momentum_body().add(torque.scale(j)),
        );
        close(
            report.kinetic_energy_change_j,
            -0.5 * (1.0 - restitution * restitution) * u * u / (wp + wb),
        );
        close_vec(report.linear_balance_residual, Vec3::ZERO);
        close_vec(report.angular_balance_residual, Vec3::ZERO);
        assert!(report.support_impulses.is_empty());
        assert_eq!(plate.q, before.q);
        assert_eq!(plate.t, before.t);
        assert_eq!(plate.steps, before.steps);
        assert_eq!(plate.history, before.history);
        assert_eq!(rigid.pose(), rigid_before.pose());
        let post = plate.clone();
        method.step(&mut plate, &dynamics, &[0.0; 9]).unwrap();
        for i in 0..9 {
            close(
                plate.v[i],
                post.v[i] * (1.0 - 0.23e-3 / 2.0) / (1.0 + 0.23e-3 / 2.0),
            );
        }
        assert_eq!(plate.steps, before.steps + 1);
        assert_eq!(plate.history.len(), before.history.len() + 1);
        let continued = RigidBodyIntegrator::new(Gravity::ZERO)
            .step(rigid, properties, Wrench::ZERO, 1e-3)
            .unwrap();
        close_vec(
            continued.state_after.linear_momentum_world(),
            rigid.linear_momentum_world(),
        );
        close_vec(
            continued.state_after.pose().position_world(),
            rigid.pose().position_world().add(
                rigid
                    .center_of_mass_velocity_world(properties)
                    .unwrap()
                    .scale(1e-3),
            ),
        );
    }
}

#[test]
fn supported_consistent_mass_recovers_inertial_reactions_and_fixed_point_limit() {
    for point in [[1.4, -0.4], [1.2, -0.7]] {
        let (mesh, model, full) = fixture(true);
        let dynamics = dynamics(&model, 1.7, 0.0);
        let prepared = PlateRigidImpact::prepare(
            dynamics,
            &mesh,
            &full,
            point,
            controls(0.4),
            budget(),
            &mut || false,
        )
        .unwrap();
        let mut plate = SecondOrderState::new(2.0, &[0.0; 6], &[0.0; 6], &[0.0; 6]);
        let properties =
            MassProperties::from_inertia_tensor(4.3, Vec3::ZERO, Mat3 { m: TENSOR }).unwrap();
        let (mut rigid, arm) = body([point[0], point[1], 1e-13], properties);
        let report = prepared
            .apply(&mut plate, &mut rigid, properties, arm, &mut || false)
            .unwrap();
        let j = report.contact.normal_impulse;
        let weight_fixed = 1.0 - (point[0] - 1.2) - (point[1] + 0.7);
        let expected = 1.7 * 0.25 * (plate.v[0] + plate.v[3]) + weight_fixed * j;
        close(report.support_impulses[0].1, expected);
        close_vec(
            report.support_linear_impulse_world,
            Vec3::new(0.0, 0.0, expected),
        );
        close_vec(
            report.support_angular_impulse_world,
            Vec3::new(-0.7 * expected, -1.2 * expected, 0.0),
        );
        close_vec(report.linear_balance_residual, Vec3::ZERO);
        close_vec(report.angular_balance_residual, Vec3::ZERO);
        if point[0] > 1.2 {
            assert!((expected - weight_fixed * j).abs() > 0.02 * j);
        } else {
            close(prepared.inverse_plate_effective_mass(), 0.0);
            assert_eq!(plate.v, vec![0.0; 6]);
            close(expected, j);
        }
    }
}

#[test]
fn compensated_pair_work_retains_loss_on_large_common_translation() {
    let (mut mesh, model, full) = fixture(false);
    mesh.nodes = vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)];
    let mut control = controls(0.0);
    // A few ulps in the iterative point mobility produce <1e-8 J per-body
    // commanded/actual work error here. The combined physical loss gate stays
    // at 1e-9 J, far below the independently known 0.125 J impact loss.
    control.impulse.work_tolerance = 1e-8;
    let prepared = PlateRigidImpact::prepare(
        dynamics(&model, 1.0, 0.0),
        &mesh,
        &full,
        [0.0, 0.0],
        control,
        budget(),
        &mut || false,
    )
    .unwrap();
    let common = 67_108_864.0;
    let mut velocity = [0.0; 9];
    for i in [0, 3, 6] {
        velocity[i] = common;
    }
    let mut plate = SecondOrderState::new(0.0, &[0.0; 9], &velocity, &[0.0; 9]);
    let properties = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap();
    let mut rigid = RigidBodyState::new(
        Pose::identity(),
        Vec3::new(0.0, 0.0, common - 1.0),
        Vec3::ZERO,
    )
    .unwrap();
    let report = prepared
        .apply(&mut plate, &mut rigid, properties, Vec3::ZERO, &mut || {
            false
        })
        .unwrap();
    close(report.contact.normal_impulse, 0.25);
    close(report.contact.dissipated_energy, 0.125);
    close(report.kinetic_energy_change_j, -0.125);
    close(report.energy_residual_j, 0.0);
}

#[test]
fn rotation_hidden_by_common_translation_cannot_be_misclassified_as_no_impact() {
    let (mut mesh, model, full) = fixture(false);
    mesh.nodes = vec![(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)];
    let prepared = PlateRigidImpact::prepare(
        dynamics(&model, 1.0, 0.0),
        &mesh,
        &full,
        [0.0, 0.0],
        controls(0.0),
        budget(),
        &mut || false,
    )
    .unwrap();
    let common = 1_152_921_504_606_846_976.0;
    let mut velocity = [0.0; 9];
    for i in [0, 3, 6] {
        velocity[i] = common;
    }
    let mut plate = SecondOrderState::new(0.0, &[0.0; 9], &velocity, &[0.0; 9]);
    let properties = MassProperties::new(1.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap();
    let arm = Vec3::new(0.0, 1.0, 0.0);
    let mut rigid = RigidBodyState::new(
        Pose::new(Vec3::new(0.0, -1.0, 0.0), UnitQuaternion::IDENTITY).unwrap(),
        Vec3::new(0.0, 0.0, common),
        Vec3::new(-1.0, 0.0, 0.0),
    )
    .unwrap();
    // Absolute point velocity rounds away the real -1 m/s rotational flow.
    assert_eq!(
        rigid
            .point_kinematics(properties, arm)
            .unwrap()
            .point_velocity_world
            .z,
        common
    );
    let before = plate.clone();
    let body_before = rigid;
    let error = prepared
        .apply(&mut plate, &mut rigid, properties, arm, &mut || false)
        .unwrap_err();
    // The contact does approach, but its finite selected impulse cannot be
    // represented in this stored velocity. Refusal is required, never a no-op.
    assert!(matches!(
        error,
        PlateRigidImpactError::Plate(fs_time::galpha::impulse::ImpulseError::Momentum { .. })
    ));
    assert_eq!(plate, before);
    assert_eq!(rigid, body_before);
}

#[test]
fn cancellation_noop_and_incompatible_mass_preserve_complete_original_states() {
    let (mesh, model, full) = fixture(false);
    let dynamics = dynamics(&model, 1.0, 0.0);
    let point = [1.4, -0.4];
    let prepared = PlateRigidImpact::prepare(
        dynamics,
        &mesh,
        &full,
        point,
        controls(0.4),
        budget(),
        &mut || false,
    )
    .unwrap();
    let original = SecondOrderState::new(3.0, &[0.0; 9], &[0.0; 9], &[0.0; 9]);
    let properties = MassProperties::new(2.0, Vec3::ZERO, Vec3::new(1.0, 1.0, 1.0)).unwrap();
    let (original_body, arm) = body([point[0], point[1], 1e-13], properties);
    let mut plate = original.clone();
    let mut rigid = original_body;
    let mut polls = 0;
    prepared
        .apply(&mut plate, &mut rigid, properties, arm, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    let mut plate = original.clone();
    let mut rigid = original_body;
    let mut count = 0;
    let err = prepared
        .apply(&mut plate, &mut rigid, properties, arm, &mut || {
            count += 1;
            count == polls
        })
        .unwrap_err();
    assert!(matches!(err, PlateRigidImpactError::Cancelled));
    assert_eq!(plate, original);
    assert_eq!(rigid, original_body);
    for (gap, speed) in [(0.1, -1.0), (0.0, 1.0)] {
        let mut rigid = RigidBodyState::new(
            Pose::new(Vec3::new(point[0], point[1], gap), UnitQuaternion::IDENTITY).unwrap(),
            Vec3::new(0.0, 0.0, 2.0 * speed),
            Vec3::ZERO,
        )
        .unwrap();
        let before = rigid;
        let report = prepared
            .apply(&mut plate, &mut rigid, properties, Vec3::ZERO, &mut || {
                false
            })
            .unwrap();
        assert!(!report.contact.impacted);
        assert_eq!(plate, original);
        assert_eq!(rigid, before);
    }
    let mut bad = Coo::new(9, 9);
    for i in 0..9 {
        let (cols, values) = full.row(i);
        for (&j, &v) in cols.iter().zip(values) {
            bad.push(i, j, v * 1.01);
        }
    }
    assert!(matches!(
        PlateRigidImpact::prepare(
            dynamics,
            &mesh,
            &bad.assemble(),
            point,
            controls(0.4),
            budget(),
            &mut || false
        ),
        Err(PlateRigidImpactError::Invalid(_))
    ));
    let mut cap = budget();
    cap.max_full_mass_nonzeros = 26;
    assert!(matches!(
        PlateRigidImpact::prepare(
            dynamics,
            &mesh,
            &full,
            point,
            controls(0.4),
            cap,
            &mut || false
        ),
        Err(PlateRigidImpactError::Budget(_))
    ));
}
