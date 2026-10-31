//! Actual rod motion under rotating material-point forces and follower couples.
use super::{controls, near, section};
use fs_couple::rod_motion::{
    PreparedRodLoads, RodMotion, RodNodalWrench, RodPointForce, SpatialRodMotion,
};
use fs_exec::CancelGate;
use fs_ga::{So3, So3Tangent, Vec3};
use fs_solid::rod::material_field::SpatialMaterialRod;
use fs_solid::rod::material_field::motion::{PreparedRodMotion, RodRootBoundary};
use fs_solid::rod::material_field::supports::RodSupport;
use fs_solid::rod::material_field::tension::{MaterialTensionLink, TensionLinkAttachment};
use std::sync::Arc;

fn rotation(v: [f64; 3]) -> So3 {
    So3::exp(So3Tangent::new(Vec3::new(v[0], v[1], v[2]))).unwrap()
}
fn rotate(r: So3, v: [f64; 3]) -> [f64; 3] {
    let v = r.rotate(Vec3::new(v[0], v[1], v[2])).unwrap();
    [v.x, v.y, v.z]
}
fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|c| a[c] + b[c])
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    (0..3).map(|c| a[c] * b[c]).sum()
}
fn vector_near(a: [f64; 3], b: [f64; 3], tol: f64) {
    for c in 0..3 { near(a[c], b[c], tol); }
}

fn torsion_model() -> Arc<PreparedRodMotion> {
    let gate = CancelGate::new();
    Arc::new(SpatialMaterialRod::try_straight(0.8, vec![section(1.0)]).unwrap()
        .with_supports(&[RodSupport { node: 1, position_m: [Some(0.8), Some(0.0), Some(0.0)], frame: None }], &gate).unwrap()
        .into_motion_model(RodRootBoundary::Clamped, &gate).unwrap())
}

#[test]
fn g0_g3_resultants_power_and_covariance_use_original_material_points() {
    let gate = CancelGate::new();
    let frame = rotation([0.4, -0.3, 0.2]);
    let global = rotation([-0.2, 0.6, 0.35]);
    let held = [RodNodalWrench { force_n: [0.2, -0.3, 0.1], torque_nm: [-0.03, 0.01, 0.02] }; 2];
    let body = [[0.04, -0.02, 0.01]; 2];
    let points = [
        RodPointForce { node: 1, offset_body_m: [0.03, 0.07, -0.02], force_world_n: [0.4, -0.2, 0.3] },
        RodPointForce { node: 1, offset_body_m: [-0.04, 0.01, 0.05], force_world_n: [-0.1, 0.3, 0.2] },
    ];
    let loads = PreparedRodLoads::new(&held, &body, &points, &gate).unwrap();
    let actual = loads.wrench_at(1, frame).unwrap();
    let mut force = held[1].force_n;
    let mut moment = add(held[1].torque_nm, rotate(frame, body[1]));
    let velocity = [0.3, -0.1, 0.2];
    let omega = [-0.4, 0.2, 0.3];
    let omega_world = rotate(frame, omega);
    let mut power = dot(held[1].force_n, velocity)
        + dot(held[1].torque_nm, omega_world) + dot(body[1], omega);
    for point in &points {
        force = add(force, point.force_world_n);
        moment = add(moment, cross(rotate(frame, point.offset_body_m), point.force_world_n));
        power += dot(point.force_world_n,
            add(velocity, rotate(frame, cross(omega, point.offset_body_m))));
    }
    vector_near(actual.force_n, force, 2e-14);
    vector_near(actual.torque_nm, moment, 2e-14);
    near(dot(actual.force_n, velocity) + dot(actual.torque_nm, omega_world), power, 2e-14);
    let rotated_held = held.map(|w| RodNodalWrench { force_n: rotate(global, w.force_n), torque_nm: rotate(global, w.torque_nm) });
    let rotated_points = points.map(|p| RodPointForce { force_world_n: rotate(global, p.force_world_n), ..p });
    let transformed = PreparedRodLoads::new(&rotated_held, &body, &rotated_points, &gate).unwrap()
        .wrench_at(1, global.compose(frame).unwrap()).unwrap();
    vector_near(transformed.force_n, rotate(global, force), 3e-14);
    vector_near(transformed.torque_nm, rotate(global, moment), 3e-14);
}

// Independent scalar ODE: I theta'' = 2 r F cos(theta) - (GJ/L) theta.
fn scalar_reference(k: f64, inertia: f64, amplitude: f64) -> [f64; 2] {
    let mut y = [0.45, 0.3];
    let h = 0.4 / 20_000.0;
    let rhs = |y: [f64; 2]| [y[1], (amplitude * y[0].cos() - k * y[0]) / inertia];
    for _ in 0..20_000 {
        let a = rhs(y);
        let b = rhs(std::array::from_fn(|c| y[c] + 0.5 * h * a[c]));
        let c = rhs(std::array::from_fn(|c| y[c] + 0.5 * h * b[c]));
        let d = rhs(std::array::from_fn(|j| y[j] + h * c[j]));
        for j in 0..2 { y[j] += h * (a[j] + 2.0 * b[j] + 2.0 * c[j] + d[j]) / 6.0; }
    }
    y
}

#[test]
fn g1_rotating_force_pair_drives_nonlinear_torsion_with_second_order_refinement() {
    let gate = CancelGate::new();
    let model = torsion_model();
    let k = model.source().sections()[0].torsional_stiffness_n_m2() / 0.8;
    let inertia = model.inertia()[1].principal_kg_m2()[0];
    let amplitude = 2.0 * 0.04 * 0.3;
    let exact = scalar_reference(k, inertia, amplitude);
    let points = [
        RodPointForce { node: 1, offset_body_m: [0.0, 0.04, 0.0], force_world_n: [0.0, 0.0, 0.3] },
        RodPointForce { node: 1, offset_body_m: [0.0, -0.04, 0.0], force_world_n: [0.0, 0.0, -0.3] },
    ];
    let loads = PreparedRodLoads::new(&[RodNodalWrench::default(); 2], &[[0.0; 3]; 2], &points, &gate).unwrap();
    let mut errors = Vec::new();
    let mut work_errors = Vec::new();
    for steps in [100, 200, 400] {
        let h = 0.4 / steps as f64;
        let mut pose = model.initial_pose();
        pose.frames[1] = rotation([0.45, 0.0, 0.0]);
        let mut state = RodMotion::new(model.clone(), pose, &[[0.0; 3]; 2], &[[0.0; 3], [0.3, 0.0, 0.0]], 0.0, &gate).unwrap();
        let mut work = 0.0;
        for _ in 0..steps {
            let (next, report) = state.step_spatial_candidate(h, &loads, controls(h), &gate).unwrap();
            work += report.source_work_quadrature_j;
            vector_near(report.angular_balance_residual_nm_s, [0.0; 3], 2e-12);
            state = next;
        }
        let angle = state.pose().frames[1].body_minus(So3::identity()).unwrap().angular.x;
        let omega = state.angular_velocities_body_rad_s().unwrap()[1][0];
        errors.push((angle - exact[0]).abs() + (omega - exact[1]).abs());
        let force_work = amplitude * (angle.sin() - 0.45_f64.sin());
        work_errors.push((work - force_work).abs());
        vector_near(state.pose().positions_m[1], [0.8, 0.0, 0.0], 0.0);
    }
    for i in 0..2 {
        assert!((3.7..4.3).contains(&(errors[i] / errors[i + 1])), "trajectory {errors:?}");
        assert!((3.7..4.3).contains(&(work_errors[i] / work_errors[i + 1])), "work {work_errors:?}");
    }
}

#[test]
fn g0_follower_couple_is_sampled_in_both_actual_rotating_frames() {
    let gate = CancelGate::new();
    let model = Arc::new(SpatialMaterialRod::try_straight(0.8, vec![section(1.0)]).unwrap()
        .into_motion_model(RodRootBoundary::Free, &gate).unwrap());
    let state = RodMotion::new(model.clone(), model.initial_pose(), &[[0.0; 3]; 2], &[[4.0, 0.0, 0.0]; 2], 0.0, &gate).unwrap();
    let body = [[0.0, 0.03, 0.0]; 2];
    let loads = PreparedRodLoads::new(&[RodNodalWrench::default(); 2], &body, &[], &gate).unwrap();
    let h = 0.002;
    let (next, report) = state.step_spatial_candidate(h, &loads, controls(h), &gate).unwrap();
    let mut expected = [0.0; 3];
    for (node, couple) in body.into_iter().enumerate() {
        expected = add(expected, add(rotate(state.pose().frames[node], couple),
            rotate(next.pose().frames[node], couple)).map(|v| 0.5 * h * v));
    }
    vector_near(report.external_angular_impulse_nm_s, expected, 1e-14);
    vector_near(report.angular_balance_residual_nm_s, [0.0; 3], 2e-12);
    assert!(expected[2].abs() > 1e-8, "the original fixed-world moment would miss this impulse");
    // A world-held couple with the same initial value retains its world axis.
    let held = [RodNodalWrench { force_n: [0.0; 3], torque_nm: body[0] }; 2];
    let world = PreparedRodLoads::new(&held, &[[0.0; 3]; 2], &[], &gate).unwrap();
    let (_, fixed) = state.step_spatial_candidate(h, &world, controls(h), &gate).unwrap();
    vector_near(fixed.external_angular_impulse_nm_s, [0.0, 2.0 * h * 0.03, 0.0], 1e-14);
}

#[test]
fn g1_point_force_at_eccentric_com_accelerates_without_false_spin() {
    let gate = CancelGate::new();
    let material = section(1.0);
    let offset = [0.02, 0.04, -0.03];
    let link = MaterialTensionLink::try_new(material.elastic().clone(), [
        TensionLinkAttachment { node: 0, offset_body_m: offset },
        TensionLinkAttachment { node: 1, offset_body_m: offset },
    ], 0.9, 1, 0.2).unwrap();
    let model = Arc::new(SpatialMaterialRod::try_straight(0.8, vec![material]).unwrap()
        .with_tension_links(vec![link], &gate).unwrap()
        .into_spatial_motion_model(RodRootBoundary::Free, &gate).unwrap());
    assert!(model.inertia()[0].inertia_com_body_kg_m2()[0][1] != 0.0);
    let frame = rotation([0.3, -0.2, 0.4]);
    let mut pose = model.initial_pose();
    for (p, r) in pose.positions_m.iter_mut().zip(&mut pose.frames) {
        *p = add(rotate(frame, *p), [0.3, -0.4, 0.1]);
        *r = frame;
    }
    let initial = pose.clone();
    let mut state = SpatialRodMotion::new_spatial(model.clone(), pose, &[[0.0; 3]; 2], &[[0.0; 3]; 2], 0.0, 1e-12, &gate).unwrap();
    let acceleration = [0.2, -0.1, 0.3];
    let points: Vec<_> = model.inertia().iter().enumerate().map(|(node, inertia)| RodPointForce {
        node, offset_body_m: inertia.center_of_mass_body_m(),
        force_world_n: acceleration.map(|a| a * inertia.mass_kg()),
    }).collect();
    let loads = PreparedRodLoads::new(&[RodNodalWrench::default(); 2], &[[0.0; 3]; 2], &points, &gate).unwrap();
    for _ in 0..10 {
        let (next, report) = state.step_spatial_candidate(0.001, &loads, controls(0.001), &gate).unwrap();
        vector_near(report.angular_balance_residual_nm_s, [0.0; 3], 2e-12);
        near(report.integration_defect_j, 0.0, 2e-12);
        state = next;
    }
    for node in 0..2 {
        vector_near(state.pose().positions_m[node], add(initial.positions_m[node], acceleration.map(|a| 0.5 * 0.01_f64.powi(2) * a)), 2e-12);
        vector_near(state.velocities_m_s().unwrap()[node], acceleration.map(|a| 0.01 * a), 2e-12);
        vector_near(state.angular_velocities_body_rad_s().unwrap()[node], [0.0; 3], 2e-11);
        near(state.pose().frames[node].body_minus(frame).unwrap().angular.norm(), 0.0, 2e-11);
    }
}

#[test]
fn g0_stationary_clamps_balance_actual_point_moments_without_work() {
    let gate = CancelGate::new();
    let model = Arc::new(SpatialMaterialRod::try_straight(0.8, vec![section(1.0)]).unwrap()
        .with_supports(&[RodSupport { node: 1, position_m: [Some(0.8), Some(0.0), Some(0.0)], frame: Some(So3::identity()) }], &gate).unwrap()
        .into_motion_model(RodRootBoundary::Clamped, &gate).unwrap());
    let state = RodMotion::new(model.clone(), model.initial_pose(), &[[0.0; 3]; 2], &[[0.0; 3]; 2], 0.0, &gate).unwrap();
    let point = RodPointForce { node: 1, offset_body_m: [0.01, 0.05, 0.02], force_world_n: [0.3, 0.2, -0.1] };
    let body = [[0.0; 3], [0.02, -0.03, 0.01]];
    let loads = PreparedRodLoads::new(&[RodNodalWrench::default(); 2], &body, &[point], &gate).unwrap();
    let h = 0.001;
    let (next, report) = state.step_spatial_candidate(h, &loads, controls(h), &gate).unwrap();
    assert_eq!(next.pose(), state.pose());
    assert_eq!(report.source_work_quadrature_j, 0.0);
    vector_near(report.support_impulses[1].linear_n_s, point.force_world_n.map(|f| -h * f), 1e-14);
    vector_near(report.support_impulses[1].angular_nm_s, add(body[1], cross(point.offset_body_m, point.force_world_n)).map(|t| -h * t), 1e-14);
    vector_near(report.angular_balance_residual_nm_s, [0.0; 3], 1e-14);
}

#[test]
fn g4_invalid_preparation_and_failed_candidates_preserve_accepted_motion() {
    let gate = CancelGate::new();
    let world = [RodNodalWrench::default(); 2];
    let body = [[0.0; 3]; 2];
    let point = RodPointForce { node: 1, offset_body_m: [0.0, 0.1, 0.0], force_world_n: [0.0, 0.0, 0.2] };
    assert!(PreparedRodLoads::new(&world, &body[..1], &[], &gate).is_err());
    assert!(PreparedRodLoads::new(&world, &body, &[RodPointForce { node: 2, ..point }], &gate).is_err());
    assert!(PreparedRodLoads::new(&world, &body, &[RodPointForce { force_world_n: [f64::NAN, 0.0, 0.0], ..point }], &gate).is_err());
    let model = torsion_model();
    let state = RodMotion::new(model.clone(), model.initial_pose(), &[[0.0; 3]; 2], &body, 0.0, &gate).unwrap();
    let loads = PreparedRodLoads::new(&world, &body, &[point], &gate).unwrap();
    let mut tight = controls(0.001);
    tight.energy_defect_limit_j = 0.0;
    assert!(state.step_spatial_candidate(0.001, &loads, tight, &gate).is_err());
    let (a, ra) = state.step_spatial_candidate(0.001, &loads, controls(0.001), &gate).unwrap();
    let (b, rb) = state.step_spatial_candidate(0.001, &loads, controls(0.001), &gate).unwrap();
    assert_eq!(a.pose(), b.pose());
    assert_eq!(a.body_spins_kg_m2_s(), b.body_spins_kg_m2_s());
    assert_eq!(ra.source_work_quadrature_j.to_bits(), rb.source_work_quadrature_j.to_bits());
    gate.request();
    assert!(PreparedRodLoads::new(&world, &body, &[point], &gate).is_err());
    assert!(state.step_spatial_candidate(0.001, &loads, controls(0.001), &gate).is_err());
    assert_eq!(state.pose(), &model.initial_pose());
    assert_eq!(state.accepted_steps(), 0);
    assert_eq!(state.time_s(), 0.0);
}
