use super::*;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::blowdown::GasBlowdownSolve;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::rocket::RocketNozzleSolve;
use fs_material::thermochemical_gas::energy::{GasEnergySolve, entropy::GasEntropySolve};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
#[path = "../../../../fs-material/examples/material_gas_spring/source.rs"]
mod source;
const MIXTURE: &str =
    include_str!("../../../../fs-material/examples/material_gas_spring/mixture-a.tsv");

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
            maximum_substep_s: 0.001,
            maximum_substeps: 1000,
        },
    }
}
fn nozzle() -> RocketNozzle {
    RocketNozzle {
        throat_area_m2: 1e-6,
        exit_area_m2: 2e-6,
    }
}
fn environment() -> GasRigidRocketEnvironment {
    GasRigidRocketEnvironment {
        inertial_frame: 7,
        ambient_pressure_pa: 0.0,
        thrust_direction_body: Vec3::new(1.0, 0.0, 0.0),
        thrust_application_arm_body_m: Vec3::new(-0.03, 0.0002, 0.0),
        co_rotating_exhaust_specific_inertia_m2: Vec3::new(0.0001, 0.001, 0.001),
        exhaust_swirl_specific_angular_momentum_body_m2_s: Vec3::ZERO,
        gravity_world_m_s2: Vec3::ZERO,
        external_force_world_n: Vec3::ZERO,
        external_torque_body_n_m: Vec3::ZERO,
    }
}
fn initial(spin: bool) -> GasRigidRocketState {
    let model = source::load(MIXTURE, "manufactured source for coupled rigid flight").unwrap();
    let gas = ClosedGasParcel::new(model, source::point(1000.0, 1e6).unwrap(), 0.001).unwrap();
    let inertia = RigidRocketInertia {
        dry_principal_kg_m2: Vec3::new(1e-5, 2e-5, 2e-5),
        propellant_specific_principal_m2: Vec3::new(0.0001, 0.001, 0.001),
    };
    let angular = if spin {
        Vec3::new(0.0, 0.0, 20.0 * (2e-5 + gas.mass_kg() * 0.001))
    } else {
        Vec3::ZERO
    };
    GasRigidRocketState::new(
        gas,
        0.01,
        Vec3::ZERO,
        Vec3::ZERO,
        UnitQuaternion::IDENTITY,
        angular,
        inertia,
        7,
    )
    .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}, difference {}",
        (a - b).abs()
    );
}
fn vector_near(a: Vec3, b: Vec3, tolerance: f64) {
    near(a.x, b.x, tolerance);
    near(a.y, b.y, tolerance);
    near(a.z, b.z, tolerance);
}

#[test]
fn accepted_source_mass_and_scalar_impulse_drive_the_original_rotating_vehicle() {
    let start = initial(true);
    let env = environment();
    let dt = 0.03;
    let gas = start
        .reservoir()
        .rocket_blowdown(
            env.ambient_pressure_pa,
            nozzle(),
            dt,
            controls().gas,
            || false,
        )
        .unwrap();
    let actual = start
        .advance(env, nozzle(), dt, controls(), || false)
        .unwrap();
    assert_eq!(actual.endpoint().reservoir(), gas.reservoir().endpoint());
    assert_eq!(&actual.step().gas, gas.reservoir().step());
    assert_eq!(&actual.step().midpoint_nozzle, gas.midpoint_flow());
    let dm = gas.reservoir().step().expelled_mass_kg;
    let expected = start
        .vehicle()
        .advance(
            RigidRocketForcing {
                mass_flow_kg_s: dm / dt,
                effective_exhaust_velocity_body_m_s: env
                    .thrust_direction_body
                    .scale(gas.axial_impulse_n_s / dm),
                thrust_application_arm_body_m: env.thrust_application_arm_body_m,
                co_rotating_exhaust_specific_inertia_m2: env
                    .co_rotating_exhaust_specific_inertia_m2,
                exhaust_swirl_specific_angular_momentum_body_m2_s: env
                    .exhaust_swirl_specific_angular_momentum_body_m2_s,
                gravity_world_m_s2: env.gravity_world_m_s2,
                external_force_world_n: env.external_force_world_n,
                external_torque_body_n_m: env.external_torque_body_n_m,
            },
            dt,
            controls().mechanics,
            || false,
        )
        .unwrap();
    assert_eq!(actual.endpoint().vehicle(), &expected.state_after);
    assert_eq!(actual.step().mechanics, expected);
    let report = actual.step();
    near(report.mechanics.burned_propellant_mass_kg, dm, 1e-17);
    near(
        report.mechanics.axial_thrust_impulse_n_s,
        gas.axial_impulse_n_s,
        1e-15,
    );
    assert!(
        report
            .mechanics
            .thrust_impulse_world_n_s
            .norm_squared()
            .sqrt()
            < 0.999 * report.mechanics.axial_thrust_impulse_n_s,
        "a finite sweep must not be rejected by comparing a vector norm with axial impulse"
    );
    assert!(actual.endpoint().vehicle().translation().velocity_m_s.y > 0.0);
    vector_near(
        report.mechanics.angular_momentum_residual_world_kg_m2_s,
        Vec3::ZERO,
        1e-15,
    );
    assert!(
        report
            .mechanics
            .exported_angular_momentum_world_kg_m2_s
            .norm_squared()
            > 0.0
    );
    assert!(report.remaining_mass_residual_kg.abs() <= controls().gas.reservoir.mass_tolerance_kg);
    assert!(report.consumed_mass_residual_kg.abs() <= controls().gas.reservoir.mass_tolerance_kg);
}

#[test]
fn pressure_thrust_contributes_once_to_axial_impulse_and_native_lever_torque() {
    let start = initial(false);
    let env = GasRigidRocketEnvironment {
        co_rotating_exhaust_specific_inertia_m2: Vec3::ZERO,
        ..environment()
    };
    let dt = 0.03;
    let vacuum = start
        .advance(env, nozzle(), dt, controls(), || false)
        .unwrap();
    let ambient_pressure = 0.4 * vacuum.step().midpoint_nozzle.exit().state().pressure;
    let ambient = start
        .advance(
            GasRigidRocketEnvironment {
                ambient_pressure_pa: ambient_pressure,
                ..env
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
    let pressure_impulse = ambient_pressure * nozzle().exit_area_m2 * dt;
    near(
        vacuum.step().mechanics.axial_thrust_impulse_n_s
            - ambient.step().mechanics.axial_thrust_impulse_n_s,
        pressure_impulse,
        1e-14,
    );
    vector_near(
        vacuum
            .step()
            .mechanics
            .applied_angular_impulse_world_kg_m2_s
            .sub(
                ambient
                    .step()
                    .mechanics
                    .applied_angular_impulse_world_kg_m2_s,
            ),
        env.thrust_application_arm_body_m
            .cross(env.thrust_direction_body)
            .scale(pressure_impulse),
        1e-17,
    );
    assert!(
        vacuum
            .endpoint()
            .vehicle()
            .angular_momentum_body_kg_m2_s()
            .z
            < ambient
                .endpoint()
                .vehicle()
                .angular_momentum_body_kg_m2_s()
                .z
    );
}

#[test]
fn evolving_source_rigid_flight_refines_position_velocity_and_attitude_without_inventory_drift() {
    let mut endpoints = Vec::new();
    for count in [8, 16, 32] {
        let mut current = initial(true);
        let initial_mass = current.reservoir().mass_kg();
        let mut consumed = 0.0;
        let dt = 0.12 / f64::from(count);
        let controls = GasRigidRocketSolve {
            mechanics: RigidRocketSolve {
                maximum_substep_s: dt / 4.0,
                maximum_substeps: 8,
            },
            ..controls()
        };
        for _ in 0..count {
            let candidate = current
                .advance(environment(), nozzle(), dt, controls, || false)
                .unwrap();
            assert!(
                candidate.step().remaining_mass_residual_kg.abs()
                    <= controls.gas.reservoir.mass_tolerance_kg
            );
            near(candidate.step().axial_impulse_residual_n_s, 0.0, 1e-15);
            let (next, report) = candidate.into_parts();
            consumed += report.gas.expelled_mass_kg;
            current = next;
        }
        near(
            current.reservoir().mass_kg() + consumed,
            initial_mass,
            1e-16,
        );
        near(
            current.reservoir().mass_kg(),
            current.vehicle().translation().propellant_mass_kg,
            1e-15,
        );
        endpoints.push(*current.vehicle());
    }
    let errors = |a: &RigidRocketState, b: &RigidRocketState| {
        [
            a.translation()
                .position_m
                .sub(b.translation().position_m)
                .norm_squared()
                .sqrt(),
            a.translation()
                .velocity_m_s
                .sub(b.translation().velocity_m_s)
                .norm_squared()
                .sqrt(),
            a.orientation()
                .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0))
                .sub(
                    b.orientation()
                        .rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0)),
                )
                .norm_squared()
                .sqrt(),
        ]
    };
    let coarse = errors(&endpoints[0], &endpoints[1]);
    let fine = errors(&endpoints[1], &endpoints[2]);
    for i in 0..3 {
        assert!(
            coarse[i] > 3.6 * fine[i],
            "component {i}: {coarse:?} versus {fine:?}"
        );
    }
}

#[test]
fn source_failure_native_budget_and_late_cancellation_preserve_both_owners() {
    let start = initial(true);
    let original = start.clone();
    let mut limited = controls();
    limited.gas.reservoir.entropy.minimum_temperature_k = 999.999;
    assert!(
        start
            .advance(environment(), nozzle(), 0.03, limited, || false)
            .is_err()
    );
    limited = controls();
    limited.mechanics.maximum_substeps = 1;
    assert!(matches!(
        start.advance(environment(), nozzle(), 0.03, limited, || false),
        Err(GasRigidRocketError::Mechanics(
            RigidRocketError::SubstepBudget { .. }
        ))
    ));
    let mut total = 0;
    start
        .advance(environment(), nozzle(), 0.03, controls(), || {
            total += 1;
            false
        })
        .unwrap();
    for stop in [1, total / 2, total - 1, total] {
        let mut polls = 0;
        assert!(
            start
                .advance(environment(), nozzle(), 0.03, controls(), || {
                    polls += 1;
                    polls == stop
                })
                .is_err()
        );
        assert_eq!(start, original);
    }
}

#[test]
fn frame_or_body_direction_mismatch_cannot_rescale_or_redirect_the_accepted_impulse() {
    let start = initial(true);
    assert!(
        start
            .advance(
                GasRigidRocketEnvironment {
                    inertial_frame: 8,
                    ..environment()
                },
                nozzle(),
                0.03,
                controls(),
                || false
            )
            .is_err()
    );
    for direction in [
        Vec3::ZERO,
        Vec3::new(2.0, 0.0, 0.0),
        Vec3::new(f64::NAN, 0.0, 0.0),
    ] {
        assert!(
            start
                .advance(
                    GasRigidRocketEnvironment {
                        thrust_direction_body: direction,
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
}

#[test]
fn explicit_closed_valve_coasts_with_the_identical_gas_and_no_engine_impulse() {
    let start = initial(true);
    let env = GasRigidRocketCoastEnvironment {
        inertial_frame: 7,
        gravity_world_m_s2: Vec3::new(0.0, 0.0, -9.0),
        external_force_world_n: Vec3::new(0.0002, 0.0, 0.0),
        external_torque_body_n_m: Vec3::new(0.0, 0.0, 1e-6),
    };
    let dt = 0.3;
    let candidate = start
        .coast(env, dt, controls().mechanics, || false)
        .unwrap();
    assert_eq!(candidate.endpoint().reservoir(), start.reservoir());
    let mechanics = candidate.mechanics();
    assert_eq!(mechanics.burned_propellant_mass_kg, 0.0);
    assert_eq!(mechanics.axial_thrust_impulse_n_s, 0.0);
    assert_eq!(mechanics.thrust_impulse_world_n_s, Vec3::ZERO);
    assert_eq!(
        mechanics.exported_angular_momentum_world_kg_m2_s,
        Vec3::ZERO
    );
    let acceleration = env.gravity_world_m_s2.add(
        env.external_force_world_n
            .scale(start.vehicle().mass_properties().unwrap().mass().recip()),
    );
    vector_near(
        candidate.endpoint().vehicle().translation().velocity_m_s,
        acceleration.scale(dt),
        1e-12,
    );
    vector_near(
        candidate.endpoint().vehicle().translation().position_m,
        acceleration.scale(0.5 * dt * dt),
        1e-12,
    );
    vector_near(
        mechanics.angular_momentum_residual_world_kg_m2_s,
        Vec3::ZERO,
        1e-14,
    );
    let mut total = 0;
    start
        .coast(env, dt, controls().mechanics, || {
            total += 1;
            false
        })
        .unwrap();
    let mut polls = 0;
    assert!(
        start
            .coast(env, dt, controls().mechanics, || {
                polls += 1;
                polls == total
            })
            .is_err()
    );
    assert_eq!(start, initial(true));
    assert!(
        start
            .coast(
                GasRigidRocketCoastEnvironment {
                    inertial_frame: 8,
                    ..env
                },
                dt,
                controls().mechanics,
                || false
            )
            .is_err()
    );
}
