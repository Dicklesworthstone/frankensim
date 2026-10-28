use super::*;
use fs_material::gas::atmosphere::StandardAtmosphere1976;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::{
    blowdown::{GasBlowdownSolve, RocketBlowdownSolve},
    rocket::RocketNozzleSolve,
};
use fs_material::thermochemical_gas::energy::{
    GasEnergySolve, entropy::GasEntropySolve, parcel::ClosedGasParcel,
};
use fs_mbd::{
    UnitQuaternion,
    rocket::rigid::{RigidRocketInertia, RigidRocketSolve},
};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
#[path = "../../../../fs-material/examples/material_gas_spring/source.rs"]
mod source;
const MIXTURE: &str =
    include_str!("../../../../fs-material/examples/material_gas_spring/mixture-a.tsv");
const FRAME: u64 = 19;
const RADIUS: f64 = 6_371_000.0;

fn atmosphere(rotating: bool) -> SphericalStandardAtmosphere {
    SphericalStandardAtmosphere::new(
        FRAME,
        Vec3::ZERO,
        RADIUS,
        9.80665 * RADIUS * RADIUS,
        if rotating {
            Vec3::new(0.0, 0.0, 1e-5)
        } else {
            Vec3::ZERO
        },
        Vec3::ZERO,
    )
    .unwrap()
}
fn state(altitude: f64, velocity: Vec3, spin: f64) -> GasRigidRocketState {
    let model = source::load(
        MIXTURE,
        "manufactured gas for atmospheric rocket integration",
    )
    .unwrap();
    let gas = ClosedGasParcel::new(model, source::point(1000.0, 1e6).unwrap(), 0.001).unwrap();
    let inertia = RigidRocketInertia {
        dry_principal_kg_m2: Vec3::new(0.001, 0.002, 0.002),
        propellant_specific_principal_m2: Vec3::new(0.0001, 0.001, 0.001),
    };
    let angular = Vec3::new(0.0, 0.0, (0.002 + gas.mass_kg() * 0.001) * spin);
    GasRigidRocketState::new(
        gas,
        0.1,
        Vec3::new(RADIUS + altitude, 0.0, 0.0),
        velocity,
        UnitQuaternion::IDENTITY,
        angular,
        inertia,
        FRAME,
    )
    .unwrap()
}
fn geometry() -> RocketFlightGeometry {
    RocketFlightGeometry {
        aerodynamic_application_arm_body_m: Vec3::ZERO,
        aerodynamic_reference_length_m: 0.1,
        aerodynamic_reference_area_m2: 1e-4,
        nozzle_ambient_arm_body_m: Vec3::new(-0.03, 0.0, 0.0),
    }
}
fn controls() -> GasRigidRocketSolve {
    GasRigidRocketSolve {
        gas: RocketBlowdownSolve {
            nozzle: RocketNozzleSolve {
                nozzle: GasEnergySolve {
                    minimum_temperature_k: 200.0,
                    maximum_temperature_k: 1800.0,
                    energy_tolerance_j_kg: 1e-7,
                    maximum_iterations: 128,
                },
                relative_mass_flow_tolerance: 1e-10,
            },
            reservoir: GasBlowdownSolve {
                entropy: GasEntropySolve {
                    minimum_temperature_k: 200.0,
                    maximum_temperature_k: 1800.0,
                    entropy_tolerance_j_kg_k: 1e-11,
                    maximum_iterations: 128,
                },
                mass_tolerance_kg: 1e-14,
            },
            maximum_donor_fraction: 0.15,
            maximum_iterations: 32,
        },
        mechanics: RigidRocketSolve {
            maximum_substep_s: 0.01,
            maximum_substeps: 2000,
        },
    }
}
fn coast() -> RocketFlightForcing {
    RocketFlightForcing {
        propulsion: RocketFlightPropulsion::ClosedValve,
        additional_force_world_n: Vec3::ZERO,
        additional_torque_body_n_m: Vec3::ZERO,
    }
}
fn powered() -> RocketFlightForcing {
    RocketFlightForcing {
        propulsion: RocketFlightPropulsion::Powered {
            nozzle: RocketNozzle {
                throat_area_m2: 1e-6,
                exit_area_m2: 2e-6,
            },
            thrust_direction_body: Vec3::new(1.0, 0.0, 0.0),
            thrust_application_arm_body_m: Vec3::new(-0.03, 0.0001, 0.0),
            co_rotating_exhaust_specific_inertia_m2: Vec3::new(0.0001, 0.001, 0.001),
            exhaust_swirl_specific_angular_momentum_body_m2_s: Vec3::ZERO,
        },
        ..coast()
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}: {}",
        (a - b).abs()
    );
}
fn vector_near(a: Vec3, b: Vec3, tolerance: f64) {
    near(a.x, b.x, tolerance);
    near(a.y, b.y, tolerance);
    near(a.z, b.z, tolerance);
}

#[test]
fn midpoint_atmosphere_drag_and_central_gravity_refine_coast_at_second_order() {
    let source = atmosphere(true);
    let original = state(1000.0, Vec3::new(0.0, 120.0, 0.0), 0.0);
    let mut endpoints = Vec::new();
    for count in [16, 32, 64] {
        let mut current = original.clone();
        for _ in 0..count {
            let next = advance_atmospheric_rocket(
                &current,
                &source,
                geometry(),
                coast(),
                2.0 / f64::from(count),
                controls(),
                |flow| {
                    assert!(flow.mach_number() < 1.0);
                    Ok(0.45)
                },
                || false,
            )
            .unwrap();
            assert!(matches!(next.step(), RocketFlightStep::Coast(_)));
            assert!(
                next.midpoint_environment()
                    .aerodynamic_load
                    .relative_power_w
                    < 0.0
            );
            assert_eq!(next.endpoint().reservoir(), original.reservoir());
            current = next.into_parts().0;
        }
        endpoints.push(*current.vehicle().translation());
    }
    let errors = |a: usize, b: usize| {
        [
            endpoints[a]
                .position_m
                .sub(endpoints[b].position_m)
                .norm_squared()
                .sqrt(),
            endpoints[a]
                .velocity_m_s
                .sub(endpoints[b].velocity_m_s)
                .norm_squared()
                .sqrt(),
        ]
    };
    let coarse = errors(0, 1);
    let fine = errors(1, 2);
    for i in 0..2 {
        assert!(
            coarse[i] > 3.6 * fine[i],
            "component {i}: {coarse:?} versus {fine:?}"
        );
    }
    assert!(endpoints[2].velocity_m_s.x < -19.0);
    assert!(endpoints[2].velocity_m_s.y < 120.0);
}

#[test]
fn powered_journey_reuses_original_gas_and_actual_body_point_environment() {
    let source = atmosphere(true);
    let position = Vec3::new(RADIUS + 1000.0, 0.0, 0.0);
    let wind = source
        .sample_at(FRAME, position)
        .unwrap()
        .air
        .wind_world_m_s();
    let mut current = state(1000.0, wind, 0.2);
    let initial_mass = current.reservoir().mass_kg();
    let model_geometry = RocketFlightGeometry {
        aerodynamic_application_arm_body_m: Vec3::new(-0.02, 0.0, 0.0),
        ..geometry()
    };
    let mut expelled = 0.0;
    let mut impulse = 0.0;
    for _ in 0..20 {
        let candidate = advance_atmospheric_rocket(
            &current,
            &source,
            model_geometry,
            powered(),
            0.01,
            controls(),
            |flow| {
                assert!(flow.reynolds_number() >= 0.0);
                Ok(0.5)
            },
            || false,
        )
        .unwrap();
        let midpoint = candidate.midpoint_environment();
        near(
            midpoint.ambient_pressure_pa,
            StandardAtmosphere1976::new()
                .at_geometric_altitude(midpoint.nozzle_altitude_m)
                .unwrap()
                .gas
                .pressure,
            0.0,
        );
        assert!(midpoint.nozzle_altitude_m < midpoint.center_altitude_m);
        let RocketFlightStep::Powered(report) = candidate.step() else {
            panic!("powered command became coast")
        };
        expelled += report.gas.expelled_mass_kg;
        impulse += report.mechanics.axial_thrust_impulse_n_s;
        near(report.remaining_mass_residual_kg, 0.0, 1e-14);
        near(report.axial_impulse_residual_n_s, 0.0, 1e-15);
        let end = candidate.endpoint();
        let vehicle = end.vehicle();
        let arm_world = vehicle
            .orientation()
            .rotate_body_to_world(model_geometry.aerodynamic_application_arm_body_m);
        let actual_point = vehicle.translation().position_m.add(arm_world);
        let air = source.sample_at(FRAME, actual_point).unwrap().air;
        let omega = vehicle
            .mass_properties()
            .unwrap()
            .angular_velocity_body(vehicle.angular_momentum_body_kg_m2_s());
        let expected = vehicle
            .translation()
            .velocity_m_s
            .sub(air.wind_world_m_s())
            .add(vehicle.orientation().rotate_body_to_world(
                omega.cross(model_geometry.aerodynamic_application_arm_body_m),
            ));
        vector_near(
            candidate
                .endpoint_environment()
                .airflow
                .relative_velocity_world_m_s(),
            expected,
            1e-13,
        );
        assert_eq!(candidate.step().mechanics().state_after, *vehicle);
        current = candidate.into_parts().0;
    }
    near(
        current.reservoir().mass_kg() + expelled,
        initial_mass,
        2e-16,
    );
    assert!(impulse > 0.0);
    assert!(current.vehicle().translation().velocity_m_s.x > 0.0);
    assert!(current.reservoir().state().state().temperature < 1000.0);
    assert!(current.vehicle().orientation() != UnitQuaternion::IDENTITY);
}

#[test]
fn final_source_and_coefficient_admission_are_atomic_after_valid_midpoint_trials() {
    let source = atmosphere(false);
    let falling = state(10.0, Vec3::new(-60.0, 0.0, 0.0), 0.0);
    let original = falling.clone();
    // At the half step all points are above the surface; the full endpoint is below it.
    let refused = advance_atmospheric_rocket(
        &falling,
        &source,
        geometry(),
        coast(),
        0.2,
        controls(),
        |_| Ok(0.0),
        || false,
    );
    assert!(matches!(
        refused,
        Err(AtmosphericFlightError::Environment(_))
    ));
    assert_eq!(falling, original);
    let free = state(1000.0, Vec3::ZERO, 0.0);
    let old = free.clone();
    let domain = |flow: &RelativeAirflow| {
        if flow.speed_m_s() > 2.0 {
            Err(FlightEnvironmentError::Invalid(
                "declared coefficient domain ends at 2 m/s",
            ))
        } else {
            Ok(0.0)
        }
    };
    // Initial and midpoint speeds are in the coefficient domain; final speed is not.
    assert!(matches!(
        advance_atmospheric_rocket(
            &free,
            &source,
            geometry(),
            coast(),
            0.3,
            controls(),
            domain,
            || false
        ),
        Err(AtmosphericFlightError::Environment(
            FlightEnvironmentError::Invalid(_)
        ))
    ));
    assert_eq!(free, old);
}

#[test]
fn late_cancellation_and_frame_mismatch_cannot_publish_a_trial_or_cutoff() {
    let source = atmosphere(false);
    let start = state(1000.0, Vec3::ZERO, 0.1);
    let original = start.clone();
    let mut total = 0;
    advance_atmospheric_rocket(
        &start,
        &source,
        geometry(),
        powered(),
        0.02,
        controls(),
        |_| Ok(0.4),
        || {
            total += 1;
            false
        },
    )
    .unwrap();
    for stop in [total / 2, total - 1, total] {
        let mut polls = 0;
        assert!(
            advance_atmospheric_rocket(
                &start,
                &source,
                geometry(),
                powered(),
                0.02,
                controls(),
                |_| Ok(0.4),
                || {
                    polls += 1;
                    polls == stop
                }
            )
            .is_err()
        );
        assert_eq!(start, original);
    }
    let wrong = SphericalStandardAtmosphere::new(
        FRAME + 1,
        Vec3::ZERO,
        RADIUS,
        9.80665 * RADIUS * RADIUS,
        Vec3::ZERO,
        Vec3::ZERO,
    )
    .unwrap();
    assert!(
        advance_atmospheric_rocket(
            &start,
            &wrong,
            geometry(),
            powered(),
            0.02,
            controls(),
            |_| Ok(0.4),
            || false
        )
        .is_err()
    );
    let invalid = RocketFlightForcing {
        propulsion: RocketFlightPropulsion::Powered {
            nozzle: RocketNozzle {
                throat_area_m2: -1.0,
                exit_area_m2: 2e-6,
            },
            thrust_direction_body: Vec3::new(1.0, 0.0, 0.0),
            thrust_application_arm_body_m: Vec3::ZERO,
            co_rotating_exhaust_specific_inertia_m2: Vec3::ZERO,
            exhaust_swirl_specific_angular_momentum_body_m2_s: Vec3::ZERO,
        },
        ..coast()
    };
    assert!(matches!(
        advance_atmospheric_rocket(
            &start,
            &source,
            geometry(),
            invalid,
            0.02,
            controls(),
            |_| Ok(0.4),
            || false
        ),
        Err(AtmosphericFlightError::Coupling(_))
    ));
    assert_eq!(start, original);
}
