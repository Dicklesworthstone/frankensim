use super::*;
use fs_ga::Mat3;
use fs_mbd::rocket::rigid::tensor::{SpecificInertia, TensorRocketInertia};

fn tensor_state() -> GasRigidRocketState<TensorRocketInertia> {
    let original = state(1000.0, Vec3::new(3.0, 12.0, -2.0), 0.0);
    let propellant = SpecificInertia::from_principal_frame(
        Vec3::new(0.25, 0.4, 0.55),
        UnitQuaternion::from_axis_angle(Vec3::new(1.0, -2.0, 0.5), 0.7).unwrap(),
    ).unwrap();
    let inertia = TensorRocketInertia::new(Mat3 {
        m: [0.0032, 0.0003, -0.00015, 0.0003, 0.0041, 0.0002,
            -0.00015, 0.0002, 0.005],
    }, propellant).unwrap();
    GasRigidRocketState::new(
        original.reservoir().clone(), 0.1, original.vehicle().translation().position_m,
        original.vehicle().translation().velocity_m_s,
        UnitQuaternion::from_axis_angle(Vec3::new(-0.2, 0.3, 1.0), 0.4).unwrap(),
        Vec3::new(0.011, -0.007, 0.018), inertia, FRAME,
    ).unwrap()
}

fn tensor_forcing() -> RocketFlightForcing<SpecificInertia> {
    RocketFlightForcing {
        propulsion: RocketFlightPropulsion::Powered {
            nozzle: RocketNozzle { throat_area_m2: 1e-6, exit_area_m2: 2e-6 },
            thrust_direction_body: Vec3::new(1.0, 0.0, 0.0),
            thrust_application_arm_body_m: Vec3::new(-0.03, 0.0001, 0.0),
            co_rotating_exhaust_specific_inertia_m2: SpecificInertia::from_principal_frame(
                Vec3::new(0.1, 0.2, 0.25),
                UnitQuaternion::from_axis_angle(Vec3::new(0.2, 1.0, -0.7), -0.8).unwrap(),
            ).unwrap(),
            exhaust_swirl_specific_angular_momentum_body_m2_s: Vec3::new(0.001, -0.002, 0.003),
        },
        additional_force_world_n: Vec3::ZERO,
        additional_torque_body_n_m: Vec3::ZERO,
    }
}

// Independent cofactor inverse, without the native principal-axis preparation.
fn inverse_apply(a: Mat3, v: Vec3) -> Vec3 {
    let m = a.m;
    let c00 = m[4] * m[8] - m[5] * m[7];
    let c01 = m[2] * m[7] - m[1] * m[8];
    let c02 = m[1] * m[5] - m[2] * m[4];
    let c10 = m[5] * m[6] - m[3] * m[8];
    let c11 = m[0] * m[8] - m[2] * m[6];
    let c12 = m[2] * m[3] - m[0] * m[5];
    let c20 = m[3] * m[7] - m[4] * m[6];
    let c21 = m[1] * m[6] - m[0] * m[7];
    let c22 = m[0] * m[4] - m[1] * m[3];
    let determinant = m[0] * c00 + m[1] * c10 + m[2] * c20;
    Vec3::new(
        c00 * v.x + c01 * v.y + c02 * v.z,
        c10 * v.x + c11 * v.y + c12 * v.z,
        c20 * v.x + c21 * v.y + c22 * v.z,
    ).scale(determinant.recip())
}

#[test]
fn tensor_atmospheric_flight_retains_source_impulse_and_full_inertia_point_motion() {
    let source = atmosphere(true);
    let start = tensor_state();
    let command = tensor_forcing();
    let shape = RocketFlightGeometry {
        aerodynamic_application_arm_body_m: Vec3::new(-0.02, 0.01, 0.03),
        ..geometry()
    };
    let actual = advance_atmospheric_rocket(
        &start, &source, shape, command, 0.02, controls(), |_| Ok(0.45), || false,
    ).unwrap();
    let midpoint = actual.midpoint_environment();
    let RocketFlightPropulsion::Powered {
        nozzle, thrust_direction_body, thrust_application_arm_body_m,
        co_rotating_exhaust_specific_inertia_m2,
        exhaust_swirl_specific_angular_momentum_body_m2_s,
    } = command.propulsion else { panic!("powered fixture required") };
    let expected = start.advance(GasRigidRocketEnvironment {
        inertial_frame: FRAME,
        ambient_pressure_pa: midpoint.ambient_pressure_pa,
        thrust_direction_body, thrust_application_arm_body_m,
        co_rotating_exhaust_specific_inertia_m2,
        exhaust_swirl_specific_angular_momentum_body_m2_s,
        gravity_world_m_s2: midpoint.gravity_world_m_s2,
        external_force_world_n: midpoint.aerodynamic_load.force_world_n,
        external_torque_body_n_m: midpoint.aerodynamic_load.torque_body_n_m,
    }, nozzle, 0.02, controls(), || false).unwrap();
    assert_eq!(actual.endpoint(), expected.endpoint());
    let RocketFlightStep::Powered(report) = actual.step() else {
        panic!("powered tensor flight became coast")
    };
    assert_eq!(report, expected.step());
    near(report.axial_impulse_residual_n_s, 0.0, 1e-14);
    near(report.remaining_mass_residual_kg, 0.0, controls().gas.reservoir.mass_tolerance_kg);
    assert_eq!(actual.endpoint().vehicle().inertia(), start.vehicle().inertia());

    let vehicle = actual.endpoint().vehicle();
    let omega = inverse_apply(vehicle.mass_properties().unwrap().inertia_tensor_body(),
        vehicle.angular_momentum_body_kg_m2_s());
    let arm = shape.aerodynamic_application_arm_body_m;
    let point = vehicle.translation().position_m
        .add(vehicle.orientation().rotate_body_to_world(arm));
    let air = source.sample_at(FRAME, point).unwrap().air;
    let expected_flow = vehicle.translation().velocity_m_s.sub(air.wind_world_m_s())
        .add(vehicle.orientation().rotate_body_to_world(omega.cross(arm)));
    vector_near(actual.endpoint_environment().airflow.relative_velocity_world_m_s(),
        expected_flow, 2e-12);
}

#[test]
fn tensor_atmospheric_coast_and_late_domain_refusal_keep_original_owners() {
    let source = atmosphere(false);
    let start = tensor_state();
    let command = RocketFlightForcing::<SpecificInertia> {
        propulsion: RocketFlightPropulsion::ClosedValve,
        additional_force_world_n: Vec3::ZERO,
        additional_torque_body_n_m: Vec3::ZERO,
    };
    let candidate = advance_atmospheric_rocket(
        &start, &source, geometry(), command, 0.03, controls(), |_| Ok(0.0), || false,
    ).unwrap();
    assert_eq!(candidate.endpoint().reservoir(), start.reservoir());
    assert_eq!(candidate.endpoint().vehicle().inertia(), start.vehicle().inertia());
    let RocketFlightStep::Coast(report) = candidate.step() else {
        panic!("closed valve became powered")
    };
    assert_eq!(report.burned_propellant_mass_kg, 0.0);
    assert_eq!(report.exported_angular_momentum_world_kg_m2_s, Vec3::ZERO);
    let before = start.vehicle().orientation().rotate_body_to_world(
        start.vehicle().angular_momentum_body_kg_m2_s());
    let after = candidate.endpoint().vehicle();
    vector_near(after.orientation().rotate_body_to_world(
        after.angular_momentum_body_kg_m2_s()), before, 2e-13);

    let falling = GasRigidRocketState::new(
        start.reservoir().clone(), 0.1, Vec3::new(RADIUS + 10.0, 0.0, 0.0),
        Vec3::new(-60.0, 0.0, 0.0), start.vehicle().orientation(),
        start.vehicle().angular_momentum_body_kg_m2_s(), start.vehicle().inertia(), FRAME,
    ).unwrap();
    let original = falling.clone();
    assert!(matches!(advance_atmospheric_rocket(
        &falling, &source, geometry(), command, 0.2, controls(), |_| Ok(0.0), || false,
    ), Err(AtmosphericFlightError::Environment(_))));
    assert_eq!(falling, original);
}
