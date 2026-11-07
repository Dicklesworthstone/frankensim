//! Multiple prescribed sources retain native dynamics and conservative work.
use fs_mbd::gravity::{
    CentralGravity, GravityError, GravityLoad, PointGravityField, PointGravitySample,
    RigidGravityField, SuperposedGravity, UniformGravity, step_in_field_with_loads,
};
use fs_mbd::{Gravity, MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3, Wrench};

const FRAME: u64 = 37;

fn body() -> MassProperties {
    MassProperties::new(3.0, Vec3::ZERO, Vec3::new(0.2, 0.3, 0.4)).unwrap()
}

fn pose(position: Vec3) -> Pose {
    Pose::new(position, UnitQuaternion::IDENTITY).unwrap()
}

fn values(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}

fn axis(i: usize, value: f64) -> Vec3 {
    let mut v = [0.0; 3];
    v[i] = value;
    Vec3::new(v[0], v[1], v[2])
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}, tolerance {tolerance:e}"
    );
}

#[test]
fn multiple_centers_and_uniform_field_differentiate_one_total_potential() {
    let a = CentralGravity::new(FRAME, Vec3::new(-2.0, 0.3, 0.5), 7.0).unwrap();
    let b = CentralGravity::new(FRAME, Vec3::new(2.5, -0.5, 1.0), 3.0).unwrap();
    let g = UniformGravity::new(
        FRAME,
        Gravity::new(Vec3::new(0.1, -0.2, -0.7)).unwrap(),
        Vec3::new(0.2, 0.5, -0.3),
    )
    .unwrap();
    let sources: [&dyn PointGravityField; 3] = [&a, &b, &g];
    let field = SuperposedGravity::new(FRAME, &sources, 3).unwrap();
    let point = Vec3::new(0.8, 1.1, -0.4);
    let sample = field.sample_point(FRAME, point).unwrap();
    let h = 1e-5;
    for j in 0..3 {
        let plus = field.sample_point(FRAME, point.add(axis(j, h))).unwrap();
        let minus = field.sample_point(FRAME, point.sub(axis(j, h))).unwrap();
        near(
            -(plus.specific_potential_j_kg - minus.specific_potential_j_kg) / (2.0 * h),
            values(sample.acceleration_world_m_s2)[j],
            1e-9,
        );
        let derivative = values(
            plus.acceleration_world_m_s2
                .sub(minus.acceleration_world_m_s2)
                .scale(0.5 / h),
        );
        for (i, value) in derivative.into_iter().enumerate() {
            near(value, sample.acceleration_gradient_s2.m[3 * i + j], 1e-9);
        }
    }
    let gradient = sample.acceleration_gradient_s2.m;
    near(gradient[0] + gradient[4] + gradient[8], 0.0, 2e-15);
    assert_eq!(field.source_count(), 3);
}

#[test]
fn compensation_retains_small_residual_between_large_balancing_fields() {
    let fields = [1e16, 1.0, -1e16].map(|g| {
        UniformGravity::new(
            FRAME,
            Gravity::new(Vec3::new(g, 0.0, 0.0)).unwrap(),
            Vec3::ZERO,
        )
        .unwrap()
    });
    let point = Vec3::new(2.0, 0.0, 0.0);
    for order in [[0, 1, 2], [2, 1, 0], [0, 2, 1]] {
        let sources = order.map(|i| &fields[i]);
        let field = SuperposedGravity::new(FRAME, &sources, 3).unwrap();
        let sample = field.sample_point(FRAME, point).unwrap();
        assert_eq!(sample.acceleration_world_m_s2, Vec3::new(1.0, 0.0, 0.0));
        assert_eq!(sample.specific_potential_j_kg, -2.0);
        let load = field.load(FRAME, pose(point), body()).unwrap();
        assert_eq!(load.wrench.force_world, Vec3::new(3.0, 0.0, 0.0));
        assert_eq!(load.potential_energy_j, -6.0);
        assert_eq!(
            field.specific_potential_increment(FRAME, point, Vec3::new(3.0, 0.0, 0.0)),
            Ok(Some(-1.0))
        );
    }
}

// An existing user field can supply endpoint potentials without a direct
// increment. Its fallback must not discard other sources' accurate work.
struct EndpointOnly(UniformGravity);
impl PointGravityField for EndpointOnly {
    fn sample_point(&self, frame: u64, position: Vec3) -> Result<PointGravitySample, GravityError> {
        self.0.sample_point(frame, position)
    }
}
impl RigidGravityField for EndpointOnly {
    fn load(
        &self,
        frame: u64,
        pose: Pose,
        mass: MassProperties,
    ) -> Result<GravityLoad, GravityError> {
        self.0.load(frame, pose, mass)
    }
}

#[test]
fn mixed_increment_support_preserves_work_hidden_by_a_large_potential_datum() {
    let direct = UniformGravity::new(
        FRAME,
        Gravity::new(Vec3::new(0.0, 0.0, -8.0)).unwrap(),
        Vec3::new(0.0, 0.0, -1e30),
    )
    .unwrap();
    let fallback = EndpointOnly(
        UniformGravity::new(
            FRAME,
            Gravity::new(Vec3::new(0.0, 0.0, -2.0)).unwrap(),
            Vec3::ZERO,
        )
        .unwrap(),
    );
    let before = Vec3::new(0.0, 0.0, 0.25);
    let after = Vec3::new(0.0, 0.0, 0.125);
    let sources: [&dyn PointGravityField; 2] = [&direct, &fallback];
    let field = SuperposedGravity::new(FRAME, &sources, 2).unwrap();
    assert_eq!(
        field
            .sample_point(FRAME, before)
            .unwrap()
            .specific_potential_j_kg,
        field
            .sample_point(FRAME, after)
            .unwrap()
            .specific_potential_j_kg
    );
    assert_eq!(
        field.specific_potential_increment(FRAME, before, after),
        Ok(Some(-1.25))
    );
    assert_eq!(
        field.specific_potential_increment(FRAME, after, before),
        Ok(Some(1.25))
    );
    let sources: [&dyn RigidGravityField; 2] = [&direct, &fallback];
    let field = SuperposedGravity::new(FRAME, &sources, 2).unwrap();
    assert_eq!(
        field.potential_energy_increment(FRAME, pose(before), pose(after), body()),
        Ok(Some(-3.75))
    );
}

#[test]
fn native_rigid_step_uses_combined_gravity_once_and_retains_physical_work() {
    let a = UniformGravity::new(
        FRAME,
        Gravity::new(Vec3::new(0.25, -0.5, -2.0)).unwrap(),
        Vec3::new(0.0, 0.0, -1e30),
    )
    .unwrap();
    let b = UniformGravity::new(
        FRAME,
        Gravity::new(Vec3::new(0.75, 0.5, -6.0)).unwrap(),
        Vec3::ZERO,
    )
    .unwrap();
    let sources = [&a, &b];
    let field = SuperposedGravity::new(FRAME, &sources, 2).unwrap();
    let initial = RigidBodyState::new(pose(Vec3::ZERO), Vec3::ZERO, Vec3::ZERO).unwrap();
    let step =
        step_in_field_with_loads(&field, FRAME, initial, body(), Wrench::ZERO, 0.25, || false)
            .unwrap();
    assert_eq!(
        step.step.state_after.pose().position_world(),
        Vec3::new(0.03125, 0.0, -0.25)
    );
    assert_eq!(
        step.step.state_after.linear_momentum_world(),
        Vec3::new(0.75, 0.0, -6.0)
    );
    assert_eq!(step.potential_energy_change_j, -6.09375);
    near(
        step.step.diagnostics_after.translational_kinetic_energy + step.potential_energy_change_j,
        0.0,
        2e-15,
    );
    let mut polls = 0;
    let refused =
        step_in_field_with_loads(&field, FRAME, initial, body(), Wrench::ZERO, 0.25, || {
            polls += 1;
            polls == 4
        });
    assert!(matches!(refused, Err(GravityError::Cancelled)));
}

#[test]
fn extended_body_sources_preserve_matching_translational_force_and_spin_torque() {
    let a = CentralGravity::new(FRAME, Vec3::ZERO, 60.0)
        .unwrap()
        .quadrupole(1.0, 0.2)
        .unwrap();
    let b = CentralGravity::new(FRAME, Vec3::new(-4.0, 8.0, 1.0), 20.0)
        .unwrap()
        .quadrupole(1.0, 0.2)
        .unwrap();
    let sources: [&dyn RigidGravityField; 2] = [&a, &b];
    let field = SuperposedGravity::new(FRAME, &sources, 2).unwrap();
    let point = Vec3::new(12.0, 3.0, -2.0);
    let total = field.load(FRAME, pose(point), body()).unwrap();
    let first = a.load(FRAME, pose(point), body()).unwrap();
    let second = b.load(FRAME, pose(point), body()).unwrap();
    let expected_torque = first.wrench.torque_body.add(second.wrench.torque_body);
    for (actual, expected) in values(total.wrench.torque_body)
        .into_iter()
        .zip(values(expected_torque))
    {
        near(actual, expected, 1e-17);
    }
    assert_ne!(total.wrench.torque_body, Vec3::ZERO);
    for i in 0..3 {
        let plus = field
            .load(FRAME, pose(point.add(axis(i, 1e-4))), body())
            .unwrap();
        let minus = field
            .load(FRAME, pose(point.sub(axis(i, 1e-4))), body())
            .unwrap();
        near(
            -(plus.potential_energy_j - minus.potential_energy_j) / 2e-4,
            values(total.wrench.force_world)[i],
            1e-9,
        );
    }
}

#[test]
fn empty_budget_frame_and_late_source_domains_are_explicit() {
    let sources: [&UniformGravity; 0] = [];
    let empty = SuperposedGravity::new(FRAME, &sources, 0).unwrap();
    let zero = empty.sample_point(FRAME, Vec3::ZERO).unwrap();
    assert_eq!(zero.acceleration_world_m_s2, Vec3::ZERO);
    assert_eq!(zero.acceleration_gradient_s2.m, [0.0; 9]);
    assert_eq!(zero.specific_potential_j_kg, 0.0);
    assert_eq!(
        empty.load(FRAME, pose(Vec3::ZERO), body()).unwrap().wrench,
        Wrench::ZERO
    );
    assert_eq!(
        empty.sample_point(FRAME + 1, Vec3::ZERO),
        Err(GravityError::FrameMismatch)
    );
    assert!(
        empty
            .sample_point(FRAME, Vec3::new(f64::NAN, 0.0, 0.0))
            .is_err()
    );

    let uniform = UniformGravity::new(FRAME, Gravity::ZERO, Vec3::ZERO).unwrap();
    let central = CentralGravity::new(FRAME, Vec3::ZERO, 1.0).unwrap();
    let sources: [&dyn PointGravityField; 2] = [&uniform, &central];
    assert!(SuperposedGravity::new(FRAME, &sources, 1).is_err());
    let field = SuperposedGravity::new(FRAME, &sources, 2).unwrap();
    let point = Vec3::new(1.0, 0.0, 0.0);
    assert!(field.sample_point(FRAME, point).is_ok());
    assert!(field.sample_point(FRAME, Vec3::ZERO).is_err());
    assert!(
        field
            .specific_potential_increment(FRAME, point, Vec3::ZERO)
            .is_err()
    );
    let foreign = CentralGravity::new(FRAME + 1, Vec3::ZERO, 1.0).unwrap();
    let sources: [&dyn PointGravityField; 2] = [&uniform, &foreign];
    let field = SuperposedGravity::new(FRAME, &sources, 2).unwrap();
    assert_eq!(
        field.sample_point(FRAME, point),
        Err(GravityError::FrameMismatch)
    );
}
