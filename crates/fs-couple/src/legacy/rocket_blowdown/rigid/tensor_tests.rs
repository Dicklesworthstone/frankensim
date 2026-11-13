use super::*;
use fs_ga::Mat3;
use fs_mbd::rocket::rigid::tensor::{SpecificInertia, TensorRocketInertia};

fn tensor_state() -> GasRigidRocketState<TensorRocketInertia> {
    let diagonal = initial(true);
    let propellant = SpecificInertia::from_principal_frame(
        Vec3::new(0.00025, 0.0004, 0.00055),
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 0.5), 0.7).unwrap(),
    )
    .unwrap();
    let inertia = TensorRocketInertia::new(
        Mat3 {
            m: [3.2e-5, 3e-6, -1.5e-6, 3e-6, 4.1e-5, 2e-6, -1.5e-6, 2e-6, 5e-5],
        },
        propellant,
    )
    .unwrap();
    GasRigidRocketState::new(
        diagonal.reservoir().clone(),
        0.01,
        Vec3::ZERO,
        Vec3::ZERO,
        UnitQuaternion::from_axis_angle(Vec3::new(-0.2, 0.3, 1.0), 0.4).unwrap(),
        Vec3::new(1.1e-4, -0.7e-4, 1.8e-4),
        inertia,
        7,
    )
    .unwrap()
}

fn tensor_environment() -> GasRigidRocketEnvironment<SpecificInertia> {
    let diagonal = environment();
    GasRigidRocketEnvironment {
        inertial_frame: diagonal.inertial_frame,
        ambient_pressure_pa: diagonal.ambient_pressure_pa,
        thrust_direction_body: diagonal.thrust_direction_body,
        thrust_application_arm_body_m: diagonal.thrust_application_arm_body_m,
        co_rotating_exhaust_specific_inertia_m2: SpecificInertia::from_principal_frame(
            Vec3::new(0.0001, 0.0002, 0.00025),
            UnitQuaternion::from_axis_angle(Vec3::new(0.2, 1.0, -0.7), -0.8).unwrap(),
        )
        .unwrap(),
        exhaust_swirl_specific_angular_momentum_body_m2_s: Vec3::new(1e-4, -2e-4, 3e-4),
        gravity_world_m_s2: diagonal.gravity_world_m_s2,
        external_force_world_n: diagonal.external_force_world_n,
        external_torque_body_n_m: diagonal.external_torque_body_n_m,
    }
}

#[test]
fn source_gas_carries_noncommuting_inertia_and_actual_angular_outflow() {
    let start = tensor_state();
    let env = tensor_environment();
    let dt = 0.03;
    let inertia = start.vehicle().mass_properties().unwrap().inertia_tensor_body();
    let exhaust = env.co_rotating_exhaust_specific_inertia_m2.tensor_body_m2();
    let ie = inertia.compose(exhaust);
    let ei = exhaust.compose(inertia);
    assert!((0..9).any(|i| (ie.m[i] - ei.m[i]).abs() > 1e-11));

    let gas = start.reservoir().rocket_blowdown(
        env.ambient_pressure_pa, nozzle(), dt, controls().gas, || false,
    ).unwrap();
    let dm = gas.reservoir().step().expelled_mass_kg;
    let actual = start.advance(env, nozzle(), dt, controls(), || false).unwrap();
    let expected = start.vehicle().advance(
        RigidRocketForcing {
            mass_flow_kg_s: dm / dt,
            effective_exhaust_velocity_body_m_s: env.thrust_direction_body
                .scale(gas.axial_impulse_n_s / dm),
            thrust_application_arm_body_m: env.thrust_application_arm_body_m,
            co_rotating_exhaust_specific_inertia_m2: env.co_rotating_exhaust_specific_inertia_m2,
            exhaust_swirl_specific_angular_momentum_body_m2_s:
                env.exhaust_swirl_specific_angular_momentum_body_m2_s,
            gravity_world_m_s2: env.gravity_world_m_s2,
            external_force_world_n: env.external_force_world_n,
            external_torque_body_n_m: env.external_torque_body_n_m,
        }, dt, controls().mechanics, || false,
    ).unwrap();
    assert_eq!(actual.endpoint().reservoir(), gas.reservoir().endpoint());
    assert_eq!(&actual.step().gas, gas.reservoir().step());
    assert_eq!(actual.step().mechanics, expected);
    assert_eq!(actual.endpoint().vehicle().inertia(), start.vehicle().inertia());
    near(actual.step().mechanics.burned_propellant_mass_kg, dm, 1e-17);
    near(actual.step().axial_impulse_residual_n_s, 0.0, 1e-15);
    near(actual.endpoint().reservoir().mass_kg(),
        start.reservoir().mass_kg() - dm, controls().gas.reservoir.mass_tolerance_kg);

    let before = start.vehicle().orientation().rotate_body_to_world(
        start.vehicle().angular_momentum_body_kg_m2_s());
    let after = actual.endpoint().vehicle();
    let retained = after.orientation().rotate_body_to_world(
        after.angular_momentum_body_kg_m2_s());
    let exported = actual.step().mechanics.exported_angular_momentum_world_kg_m2_s;
    assert!(exported.norm_squared() > 0.0);
    vector_near(retained.add(exported),
        before.add(actual.step().mechanics.applied_angular_impulse_world_kg_m2_s), 2e-14);
}

#[test]
fn tensor_coast_keeps_the_gas_and_total_spin_and_refusals_publish_neither_owner() {
    let start = tensor_state();
    let original = start.clone();
    let env = GasRigidRocketCoastEnvironment {
        inertial_frame: 7,
        gravity_world_m_s2: Vec3::new(0.0, 0.0, -9.0),
        external_force_world_n: Vec3::ZERO,
        external_torque_body_n_m: Vec3::ZERO,
    };
    let candidate = start.coast(env, 0.03, controls().mechanics, || false).unwrap();
    assert_eq!(candidate.endpoint().reservoir(), start.reservoir());
    assert_eq!(candidate.endpoint().vehicle().inertia(), start.vehicle().inertia());
    assert_eq!(candidate.mechanics().burned_propellant_mass_kg, 0.0);
    assert_eq!(candidate.mechanics().axial_thrust_impulse_n_s, 0.0);
    assert_eq!(candidate.mechanics().exported_angular_momentum_world_kg_m2_s, Vec3::ZERO);
    let before = start.vehicle().orientation().rotate_body_to_world(
        start.vehicle().angular_momentum_body_kg_m2_s());
    let after = candidate.endpoint().vehicle();
    vector_near(after.orientation().rotate_body_to_world(
        after.angular_momentum_body_kg_m2_s()), before, 2e-14);
    vector_near(after.translation().velocity_m_s, Vec3::new(0.0, 0.0, -0.27), 1e-13);

    let mut limited = controls();
    limited.mechanics.maximum_substeps = 1;
    assert!(matches!(
        start.advance(tensor_environment(), nozzle(), 0.03, limited, || false),
        Err(GasRigidRocketError::Mechanics(RigidRocketError::SubstepBudget { .. }))
    ));
    let mut polls = 0;
    start.advance(tensor_environment(), nozzle(), 0.03, controls(), || {
        polls += 1;
        false
    }).unwrap();
    let mut observed = 0;
    assert!(start.advance(tensor_environment(), nozzle(), 0.03, controls(), || {
        observed += 1;
        observed == polls
    }).is_err());
    assert_eq!(start, original);
}
