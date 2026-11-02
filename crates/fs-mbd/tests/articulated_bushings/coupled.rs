//! Full constitutive coupling through the original geometry and native motion.
use super::*;
use fs_mbd::articulated::motion::spatial_loads::bushing::coupled::{
    self as full, BushingCoordinateScales, CoupledSpatialBushing, PassiveBushingMatrix,
};

fn matrix(diagonal: f64, direction: [f64; 6]) -> [[f64; 6]; 6] {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| direction[i] * direction[j] + if i == j { diagonal } else { 0. })
    })
}
fn law(k: [[f64; 6]; 6], c: [[f64; 6]; 6], b: SpatialBushing) -> CoupledSpatialBushing {
    let scales = BushingCoordinateScales {
        translation_m: 0.25,
        rotation_rad: 1.,
    };
    CoupledSpatialBushing::try_new(
        b.a,
        b.b,
        b.rest_translation_a,
        PassiveBushingMatrix::try_new(k, scales).unwrap(),
        PassiveBushingMatrix::try_new(c, scales).unwrap(),
    )
    .unwrap()
}
fn strain(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    b: &CoupledSpatialBushing,
) -> [f64; 6] {
    let fk = forward_kinematics(
        model,
        BaseState::stationary(Se3::identity()),
        &state.joint_position,
        &state.joint_velocity,
    )
    .unwrap();
    let [a, c] = b.frames().map(|f| {
        fk.world_from_link[f.link]
            .compose(f.link_from_attachment)
            .unwrap()
    });
    let relative = a.inverse().unwrap().compose(c).unwrap();
    let q = relative.translation() - b.rest_translation_a();
    let r = relative.rotation().log().angular;
    [q.x, q.y, q.z, r.x, r.y, r.z]
}
fn quadratic(matrix: &[[f64; 6]; 6], x: [f64; 6]) -> f64 {
    (0..6)
        .map(|i| (0..6).map(|j| x[i] * matrix[i][j] * x[j]).sum::<f64>())
        .sum()
}
fn energy(
    model: &ArticulatedModel,
    state: &FreeFloatingMotionState,
    bs: &[CoupledSpatialBushing],
) -> f64 {
    bs.iter()
        .map(|b| 0.5 * quadratic(b.stiffness().native_coefficients(), strain(model, state, b)))
        .sum()
}

#[test]
fn cross_axis_and_translation_rotation_forces_preserve_native_work_and_objectivity() {
    let (model, state, diagonal) = mixed();
    let bs = diagonal.map(|b| {
        law(
            matrix(3., [1., 2., -1., 3., -2., 1.]),
            matrix(0.2, [0.4, -0.3, 0.2, 0.5, 0.2, -0.6]),
            b,
        )
    });
    let response = full::evaluate_bushing_forces(&model, &state, &bs, || false).unwrap();
    close(response.stored_energy, energy(&model, &state, &bs), 2e-12);
    let epsilon = 2e-6;
    let mut plus = state.clone();
    let mut minus = state.clone();
    for i in 0..model.dof_count() {
        plus.joint_position[i] += epsilon * state.joint_velocity[i];
        minus.joint_position[i] -= epsilon * state.joint_velocity[i];
    }
    let derivative = (energy(&model, &plus, &bs) - energy(&model, &minus, &bs)) / (2. * epsilon);
    let independent_loss: f64 = bs
        .iter()
        .map(|b| {
            let a = strain(&model, &plus, b);
            let z = strain(&model, &minus, b);
            quadratic(
                b.damping().native_coefficients(),
                core::array::from_fn(|i| (a[i] - z[i]) / (2. * epsilon)),
            )
        })
        .sum();
    close(response.dissipation_power, independent_loss, 3e-9);
    close(
        response.internal_power + derivative + response.dissipation_power,
        0.,
        3e-8,
    );
    assert!(response.dissipation_power > 0.);
    let fk = forward_kinematics(
        &model,
        BaseState::prescribed(
            state.base.world_from_base,
            state.base.twist_body,
            Twist::zero(),
        ),
        &state.joint_position,
        &state.joint_velocity,
    )
    .unwrap();
    let mut force = Z;
    let mut torque = Z;
    let mut power = 0.;
    for (i, w) in response.body_wrench.iter().enumerate() {
        let pose = fk.world_from_link[i];
        let f = pose.rotation().rotate(w.force).unwrap();
        force = force + f;
        torque = torque + pose.rotation().rotate(w.torque).unwrap() + pose.translation().cross(f);
        power += w.pairing(fk.body_twist[i]);
    }
    vector(force, Z, 3e-12);
    vector(torque, Z, 3e-12);
    close(power, response.internal_power, 3e-12);
    let mut moving = state.clone();
    moving.base = FreeFloatingBaseState::new(
        pose(Z, Vec3::new(1e12, -2e12, 3e12)),
        Twist::new(Vec3::new(1e8, 2e8, -3e8), Vec3::new(1e10, -2e10, 3e10)),
    );
    assert_eq!(
        response,
        full::evaluate_bushing_forces(&model, &moving, &bs, || false).unwrap()
    );
}

#[test]
fn prepared_diagonal_laws_reproduce_both_original_native_integrators() {
    let (model, state, bs) = mixed();
    let full: Vec<_> = bs
        .iter()
        .map(|b| {
            CoupledSpatialBushing::try_from_diagonal(
                *b,
                BushingCoordinateScales {
                    translation_m: 0.25,
                    rotation_rad: 1.,
                },
            )
            .unwrap()
        })
        .collect();
    assert_eq!(
        evaluate_bushing_forces(&model, &state, &bs, || false).unwrap(),
        full::evaluate_bushing_forces(&model, &state, &full, || false).unwrap()
    );
    let efforts = [0.; 3];
    let wrenches = [ZERO; 4];
    let loads = held(&efforts, &wrenches);
    let h = 0.002;
    assert_eq!(
        step_free_floating_bushings(&model, &state, loads, &bs, h, || false).unwrap(),
        full::step_free_floating_bushings(&model, &state, loads, &full, h, || false).unwrap()
    );
    assert_eq!(
        step_free_floating_bushings_implicit(
            &model,
            &state,
            loads,
            &bs,
            h,
            ImplicitSpatialConfig::default(),
            || false
        )
        .unwrap(),
        full::step_free_floating_bushings_implicit(
            &model,
            &state,
            loads,
            &full,
            h,
            ImplicitSpatialConfig::default(),
            || false
        )
        .unwrap()
    );
}

#[test]
fn full_matrices_drive_actual_coupled_motion_and_keep_physical_heat_separate() {
    let (model, state, diagonal) = mixed();
    let bs = diagonal.map(|b| {
        law(
            matrix(5., [1., 2., -1., 3., -2., 1.]),
            matrix(0.5, [0.4, -0.3, 0.2, 0.5, 0.2, -0.6]),
            b,
        )
    });
    let efforts = [0.; 3];
    let wrenches = [ZERO; 4];
    let loads = held(&efforts, &wrenches);
    let duration = 0.0002;
    let step = full::step_free_floating_bushings_implicit(
        &model,
        &state,
        loads,
        &bs,
        duration,
        ImplicitSpatialConfig {
            residual_tolerance: 1e-12,
            ..Default::default()
        },
        || false,
    )
    .unwrap();
    let midpoint =
        full::evaluate_bushing_forces(&model, &step.midpoint_state, &bs, || false).unwrap();
    close(
        step.spatial.dissipated_energy,
        duration * midpoint.dissipation_power,
        2e-15,
    );
    close(
        step.spatial.stored_energy_after,
        energy(&model, &step.spatial.motion.state_after, &bs),
        3e-12,
    );
    balance(&step.spatial.motion, 3e-12);
    assert!(step.spatial.dissipated_energy > 0.);
    assert!(step.spatial.mechanical_work_residual.abs() < 2e-8);
    assert_ne!(
        step.spatial.motion.state_after.joint_velocity,
        state.joint_velocity
    );
    // A physical cross tensor changes the trajectory; dropping cross entries
    // is not an equivalent orientation or principal-axis approximation.
    let uncoupled: Vec<_> = bs
        .iter()
        .zip(diagonal)
        .map(|(b, frame)| {
            let diag = |a: &[[f64; 6]; 6]| {
                core::array::from_fn(|i| {
                    core::array::from_fn(|j| if i == j { a[i][i] } else { 0. })
                })
            };
            law(
                diag(b.stiffness().native_coefficients()),
                diag(b.damping().native_coefficients()),
                frame,
            )
        })
        .collect();
    let other = full::step_free_floating_bushings_implicit(
        &model,
        &state,
        loads,
        &uncoupled,
        duration,
        ImplicitSpatialConfig::default(),
        || false,
    )
    .unwrap();
    assert!(
        step.spatial
            .motion
            .state_after
            .joint_velocity
            .iter()
            .zip(&other.spatial.motion.state_after.joint_velocity)
            .any(|(a, b)| (a - b).abs() > 1e-6)
    );
    let mut calls = 0;
    full::step_free_floating_bushings_implicit(
        &model,
        &state,
        loads,
        &bs,
        duration,
        ImplicitSpatialConfig::default(),
        || {
            calls += 1;
            false
        },
    )
    .unwrap();
    let mut visited = 0;
    assert!(
        full::step_free_floating_bushings_implicit(
            &model,
            &state,
            loads,
            &bs,
            duration,
            ImplicitSpatialConfig::default(),
            || {
                visited += 1;
                visited == calls
            }
        )
        .is_err()
    );
    assert_eq!(visited, calls);
}
