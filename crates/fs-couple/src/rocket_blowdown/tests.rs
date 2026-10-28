//! G0 cross-owner inventory/impulse and G1 smooth propulsion refinement.

use super::*;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::blowdown::GasBlowdownSolve;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::rocket::RocketNozzleSolve;
use fs_material::thermochemical_gas::energy::{GasEnergySolve, entropy::GasEntropySolve};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
// Existing authored synthetic source, original NASA-9/material-card loader;
// these tests make no measured-gas or production propellant calibration claim.
#[path = "../../../fs-material/examples/material_gas_spring/source.rs"]
mod source;
const MIXTURE: &str =
    include_str!("../../../fs-material/examples/material_gas_spring/mixture-a.tsv");

fn controls() -> RocketBlowdownSolve {
    RocketBlowdownSolve {
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
    }
}
fn nozzle() -> RocketNozzle {
    RocketNozzle {
        throat_area_m2: 1e-6,
        exit_area_m2: 2e-6,
    }
}
fn initial() -> GasRocketState {
    let model = source::load(MIXTURE, "rocket coupled finite reservoir test").unwrap();
    let gas = ClosedGasParcel::new(model, source::point(1000.0, 1e6).unwrap(), 0.001).unwrap();
    GasRocketState::new(gas, 0.01, Vec3::ZERO, Vec3::ZERO, 7).unwrap()
}
fn environment() -> GasRocketEnvironment {
    GasRocketEnvironment {
        inertial_frame: 7,
        thrust_direction_world: Vec3::new(0.6, 0.8, 0.0),
        ambient_pressure_pa: 0.0,
        gravity_m_s2: Vec3::ZERO,
        external_force_n: Vec3::ZERO,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1e-10),
        "{a} != {b}"
    );
}
fn vector_near(a: Vec3, b: Vec3, tolerance: f64) {
    for (a, b) in [(a.x, b.x), (a.y, b.y), (a.z, b.z)] {
        near(a, b, tolerance);
    }
}

#[test]
fn original_accepted_gas_mass_and_impulse_drive_one_analytic_variable_mass_step() {
    let initial = initial();
    assert_eq!(
        initial.vehicle().propellant_mass_kg,
        initial.reservoir().mass_kg()
    );
    let dt = 0.03;
    let standalone = initial
        .reservoir()
        .rocket_blowdown(0.0, nozzle(), dt, controls(), || false)
        .unwrap();
    let candidate = initial
        .advance(environment(), nozzle(), dt, controls(), || false)
        .unwrap();
    assert_eq!(
        candidate.endpoint().reservoir(),
        standalone.reservoir().endpoint()
    );
    assert_eq!(&candidate.step().gas, standalone.reservoir().step());
    assert_eq!(
        &candidate.step().midpoint_nozzle,
        standalone.midpoint_flow()
    );
    let dm = standalone.reservoir().step().expelled_mass_kg;
    let mass0 = initial.vehicle().dry_mass_kg + initial.reservoir().mass_kg();
    let mass1 = initial.vehicle().dry_mass_kg + candidate.endpoint().reservoir().mass_kg();
    let dv = (standalone.axial_impulse_n_s / dm) * (mass0 / mass1).ln();
    vector_near(
        candidate.endpoint().vehicle().velocity_m_s,
        environment().thrust_direction_world.scale(dv),
        1e-11,
    );
    vector_near(
        candidate.step().mechanics.thrust_impulse_n_s,
        environment()
            .thrust_direction_world
            .scale(standalone.axial_impulse_n_s),
        1e-14,
    );
    near(
        candidate.step().mechanics.burned_propellant_mass_kg,
        dm,
        1e-14,
    );
    assert_eq!(
        candidate.endpoint().vehicle().dry_mass_kg,
        initial.vehicle().dry_mass_kg
    );
    assert_eq!(
        candidate.endpoint().inertial_frame(),
        initial.inertial_frame()
    );
    assert!(
        candidate.step().remaining_mass_residual_kg.abs() <= controls().reservoir.mass_tolerance_kg
    );
    assert!(
        candidate.step().consumed_mass_residual_kg.abs() <= controls().reservoir.mass_tolerance_kg
    );
    let endpoint_rate = candidate
        .endpoint()
        .reservoir()
        .rocket_nozzle_flow(0.0, nozzle(), controls().nozzle, || false)
        .unwrap()
        .mass_flow_kg_s;
    assert!(
        (dm - endpoint_rate * dt).abs() > 1e-9,
        "endpoint resampling would consume a different mass"
    );
}

#[test]
fn ambient_pressure_thrust_is_included_exactly_once_in_vehicle_motion() {
    let initial = initial();
    let dt = 0.03;
    let vacuum = initial
        .advance(environment(), nozzle(), dt, controls(), || false)
        .unwrap();
    let pa = 0.5 * vacuum.step().midpoint_nozzle.exit().state().pressure;
    let ambient = initial
        .advance(
            GasRocketEnvironment {
                ambient_pressure_pa: pa,
                ..environment()
            },
            nozzle(),
            dt,
            controls(),
            || false,
        )
        .unwrap();
    assert_eq!(
        vacuum.endpoint().reservoir(),
        ambient.endpoint().reservoir()
    );
    let expected_loss = pa * nozzle().exit_area_m2 * dt;
    let actual_loss = vacuum
        .step()
        .mechanics
        .thrust_impulse_n_s
        .sub(ambient.step().mechanics.thrust_impulse_n_s);
    vector_near(
        actual_loss,
        environment().thrust_direction_world.scale(expected_loss),
        1e-11,
    );
    assert!(
        vacuum.endpoint().vehicle().velocity_m_s.x > ambient.endpoint().vehicle().velocity_m_s.x
    );
}

#[test]
fn variable_caloric_tank_thrust_and_motion_refine_second_order_without_inventory_drift() {
    let mut final_states = Vec::new();
    let environment = GasRocketEnvironment {
        gravity_m_s2: Vec3::new(0.0, 0.0, -9.80665),
        external_force_n: Vec3::new(-0.0001, 0.0, 0.0),
        ..environment()
    };
    for count in [8, 16, 32] {
        let mut state = initial();
        let initial_mass = state.reservoir().mass_kg();
        let mut expelled = 0.0;
        for _ in 0..count {
            let next = state
                .advance(
                    environment,
                    nozzle(),
                    0.3 / f64::from(count),
                    controls(),
                    || false,
                )
                .unwrap();
            assert!(
                next.step().remaining_mass_residual_kg.abs()
                    <= controls().reservoir.mass_tolerance_kg
            );
            let (endpoint, step) = next.into_parts();
            expelled += step.gas.expelled_mass_kg;
            state = endpoint;
        }
        near(state.reservoir().mass_kg() + expelled, initial_mass, 1e-12);
        near(
            state.vehicle().propellant_mass_kg,
            state.reservoir().mass_kg(),
            1e-12,
        );
        near(state.vehicle().velocity_m_s.z, -9.80665 * 0.3, 1e-13);
        final_states.push(state);
    }
    for component in [0, 1] {
        let value = |state: &GasRocketState| {
            if component == 0 {
                state.vehicle().velocity_m_s.x
            } else {
                state.vehicle().position_m.x
            }
        };
        let coarse = (value(&final_states[0]) - value(&final_states[1])).abs();
        let fine = (value(&final_states[1]) - value(&final_states[2])).abs();
        assert!(
            coarse > 3.7 * fine,
            "component {component}: {coarse} vs {fine}"
        );
    }
}

#[test]
fn source_failure_and_late_cancellation_preserve_both_original_owners() {
    let initial = initial();
    let original = initial.clone();
    let mut limited = controls();
    limited.reservoir.entropy.minimum_temperature_k = 999.999;
    assert!(
        initial
            .advance(environment(), nozzle(), 0.03, limited, || false)
            .is_err()
    );
    assert!(
        initial
            .advance(environment(), nozzle(), 1000.0, controls(), || false)
            .is_err()
    );
    let mut total = 0;
    initial
        .advance(environment(), nozzle(), 0.03, controls(), || {
            total += 1;
            false
        })
        .unwrap();
    for stop in [1, total / 2, total - 1, total] {
        let mut calls = 0;
        assert!(
            initial
                .advance(environment(), nozzle(), 0.03, controls(), || {
                    calls += 1;
                    calls == stop
                })
                .is_err()
        );
        assert_eq!(initial, original);
    }
}

#[test]
fn invalid_vehicle_direction_and_frame_are_rejected_without_rescaling_impulse() {
    let initial = initial();
    for dry_mass in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            GasRocketState::new(
                initial.reservoir().clone(),
                dry_mass,
                Vec3::ZERO,
                Vec3::ZERO,
                7
            )
            .is_err()
        );
    }
    for direction in [
        Vec3::ZERO,
        Vec3::new(2.0, 0.0, 0.0),
        Vec3::new(f64::NAN, 0.0, 0.0),
        Vec3::new(1e300, 0.0, 0.0),
    ] {
        assert!(matches!(
            initial.advance(
                GasRocketEnvironment {
                    thrust_direction_world: direction,
                    ..environment()
                },
                nozzle(),
                0.03,
                controls(),
                || false
            ),
            Err(GasRocketError::Invalid(_))
        ));
    }
    assert!(matches!(
        initial.advance(
            GasRocketEnvironment {
                inertial_frame: 8,
                ..environment()
            },
            nozzle(),
            0.03,
            controls(),
            || false
        ),
        Err(GasRocketError::Invalid(_))
    ));
    assert!(
        initial
            .advance(
                GasRocketEnvironment {
                    gravity_m_s2: Vec3::new(0.0, 0.0, f64::INFINITY),
                    ..environment()
                },
                nozzle(),
                0.03,
                controls(),
                || false
            )
            .is_err()
    );
}
