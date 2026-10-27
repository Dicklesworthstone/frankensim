//! G1: native mass actions on independently known coupled joint inertia.
use fs_ga::{Mat3, Se3, Vec3, Wrench};
use fs_mbd::articulated::{
    ArticulatedError, ArticulatedModel, BaseState, JointLimits, JointModel, Link, SpatialInertia,
    forward_dynamics, stationary_base_joint_momentum, stationary_base_joint_velocity,
};

fn model() -> ArticulatedModel {
    let limits = JointLimits::new(-10., 10., 20., 1e-6).unwrap();
    let joint = JointModel::revolute(Vec3::new(0., 0., 1.), Some(limits)).unwrap();
    let inertia = |i| {
        SpatialInertia::new(
            2.,
            Vec3::default(),
            Mat3 {
                m: [i, 0., 0., 0., i, 0., 0., 0., i],
            },
        )
        .unwrap()
    };
    ArticulatedModel::new(vec![
        Link::new(
            "base",
            None,
            Se3::identity(),
            JointModel::FIXED,
            inertia(1.),
        ),
        Link::new("first", Some(0), Se3::identity(), joint, inertia(2.)),
        Link::new("second", Some(1), Se3::identity(), joint, inertia(3.)),
    ])
    .unwrap()
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-12, "{a} != {b}");
}

#[test]
fn native_recursions_apply_the_coupled_mass_law_without_effort_limit_confusion() {
    let model = model();
    let q = [0.4, -0.2];
    let v = [3., -1.];
    // T = 1/2 I1 v1^2 + 1/2 I2 (v1+v2)^2.
    // M = [[5,3],[3,3]], p = [12,6].
    let p = stationary_base_joint_momentum(&model, &q, &v).unwrap();
    close(p[0], 12.);
    close(p[1], 6.);
    let recovered = stationary_base_joint_velocity(&model, &q, &p).unwrap();
    close(recovered[0], 3.);
    close(recovered[1], -1.);
    // The same numbers as physical force MUST still violate actuator limits.
    assert!(matches!(
        forward_dynamics(
            &model,
            BaseState::stationary(Se3::identity()),
            &q,
            &[0.; 2],
            &p,
            Vec3::default(),
            &[Wrench::default(); 3]
        ),
        Err(ArticulatedError::JointEffortLimitViolation { .. })
    ));
}

#[test]
fn mass_actions_retain_position_rate_dimension_and_finite_admission() {
    let model = model();
    assert!(stationary_base_joint_momentum(&model, &[0.; 2], &[21., 0.]).is_err());
    assert!(stationary_base_joint_velocity(&model, &[11., 0.], &[0.; 2]).is_err());
    assert!(stationary_base_joint_velocity(&model, &[0.; 2], &[1.]).is_err());
    assert!(stationary_base_joint_velocity(&model, &[0.; 2], &[f64::NAN, 0.]).is_err());
    assert!(stationary_base_joint_velocity(&model, &[0.; 2], &[1e6, 0.]).is_err());
    assert_eq!(
        stationary_base_joint_velocity(&model, &[0.; 2], &[0.; 2]).unwrap(),
        vec![0.; 2]
    );
}
