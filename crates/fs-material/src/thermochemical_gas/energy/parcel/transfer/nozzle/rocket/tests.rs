//! G1 analytic nozzle relations and G3 source/area/ambient metamorphics.
use super::*;
use crate::thermochemical_gas::tests::{model_in_domains, point};

fn controls() -> RocketNozzleSolve {
    RocketNozzleSolve {
        nozzle: GasEnergySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 4000.0,
            energy_tolerance_j_kg: 1e-6,
            maximum_iterations: 128,
        },
        relative_mass_flow_tolerance: 1e-11,
    }
}
fn gas(slope: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(
        model_in_domains(
            vec![0.75, 0.25],
            slope,
            [200.0, 4000.0],
            [200.0, 4000.0],
            [1.0, 1e8],
        ),
        point(1800.0, 3e6),
        0.01,
    )
    .unwrap()
}
fn nozzle() -> RocketNozzle {
    RocketNozzle {
        throat_area_m2: 1e-4,
        exit_area_m2: 4e-4,
    }
}
fn close(a: f64, b: f64, relative: f64) {
    assert!((a - b).abs() <= relative * b.abs().max(1e-15), "{a} != {b}");
}
fn area_ratio(mach: f64, gamma: f64) -> f64 {
    ((2.0 / (gamma + 1.0)) * (1.0 + 0.5 * (gamma - 1.0) * mach * mach))
        .powf((gamma + 1.0) / (2.0 * (gamma - 1.0)))
        / mach
}

#[test]
fn g1_supersonic_nozzle_matches_analytic_area_energy_pressure_and_thrust() {
    let reservoir = gas(0.0);
    let state = reservoir.state().state();
    let (gamma, r, t0, p0) = (
        state.gamma,
        state.specific_gas_constant,
        state.temperature,
        state.pressure,
    );
    for mach in [1.0, 1.2, 2.0, 3.0] {
        let geometry = RocketNozzle {
            exit_area_m2: 1e-4 * area_ratio(mach, gamma),
            ..nozzle()
        };
        let flow = reservoir
            .rocket_nozzle_flow(0.0, geometry, controls(), || false)
            .unwrap();
        let factor = 1.0 + 0.5 * (gamma - 1.0) * mach * mach;
        let temperature = t0 / factor;
        let pressure = p0 / factor.powf(gamma / (gamma - 1.0));
        let velocity = mach * (gamma * r * temperature).sqrt();
        let rate = geometry.throat_area_m2
            * p0
            * (gamma / (r * t0)).sqrt()
            * (2.0 / (gamma + 1.0)).powf((gamma + 1.0) / (2.0 * (gamma - 1.0)));
        let thrust = rate * velocity + pressure * geometry.exit_area_m2;
        close(flow.mass_flow_kg_s, rate, 2e-10);
        close(flow.exit_mach_number, mach, 2e-10);
        close(flow.exit().state().temperature, temperature, 2e-10);
        close(flow.exit().state().pressure, pressure, 2e-10);
        close(flow.exit_velocity_m_s, velocity, 2e-10);
        close(flow.thrust_n, thrust, 2e-10);
        close(
            flow.characteristic_velocity_m_s,
            p0 * geometry.throat_area_m2 / rate,
            2e-10,
        );
        close(
            flow.specific_impulse_s,
            thrust / (rate * STANDARD_GRAVITY_M_S2),
            2e-10,
        );
        assert!(flow.relative_mass_flow_residual.abs() <= controls().relative_mass_flow_tolerance);
        assert!(flow.energy_residual_j_kg.abs() < 2e-6);
        assert!(flow.entropy_energy_residual_j_kg.abs() < 2e-6);
    }
}

#[test]
fn g1_variable_cp_uses_original_enthalpy_and_entropy_integrals() {
    let reservoir = gas(0.002);
    let flow = reservoir
        .rocket_nozzle_flow(0.0, nozzle(), controls(), || false)
        .unwrap();
    let state = reservoir.state().state();
    let (t0, te, ts, r) = (
        state.temperature,
        flow.exit().state().temperature,
        flow.throat().state().temperature,
        state.specific_gas_constant,
    );
    // Independent exact integrals of cp/R = 3.25 + 0.0015 T.
    let (a, b) = (3.25, 0.0015);
    let delta_h = r * (a * (te - t0) + 0.5 * b * (te - t0) * (te + t0));
    close(flow.exit_velocity_m_s, (-2.0 * delta_h).sqrt(), 2e-11);
    close(
        flow.exit().state().pressure,
        state.pressure * (a * (te / t0).ln() + b * (te - t0)).exp(),
        2e-11,
    );
    let throat_cp_r = a + b * ts;
    let throat_v2 = throat_cp_r / (throat_cp_r - 1.0) * r * ts;
    close(
        -2.0 * r * (a * (ts - t0) + 0.5 * b * (ts - t0) * (ts + t0)),
        throat_v2,
        2e-11,
    );
    close(
        flow.exit().state().density * flow.exit_velocity_m_s * nozzle().exit_area_m2,
        flow.mass_flow_kg_s,
        2e-11,
    );
    let frozen_area = area_ratio(flow.exit_mach_number, state.gamma);
    assert!(
        (frozen_area - 4.0).abs() > 0.05,
        "upstream gamma substitution must be detectable"
    );
}

#[test]
fn g3_ambient_pressure_changes_only_pressure_thrust_and_areas_scale_extensive_outputs() {
    let reservoir = gas(0.002);
    let vacuum = reservoir
        .rocket_nozzle_flow(0.0, nozzle(), controls(), || false)
        .unwrap();
    let pa = 0.5 * vacuum.exit().state().pressure;
    let ambient = reservoir
        .rocket_nozzle_flow(pa, nozzle(), controls(), || false)
        .unwrap();
    assert_eq!(vacuum.exit(), ambient.exit());
    assert_eq!(vacuum.mass_flow_kg_s, ambient.mass_flow_kg_s);
    close(
        vacuum.thrust_n - ambient.thrust_n,
        pa * nozzle().exit_area_m2,
        1e-12,
    );
    let matched = reservoir
        .rocket_nozzle_flow(vacuum.exit().state().pressure, nozzle(), controls(), || {
            false
        })
        .unwrap();
    assert_eq!(matched.pressure_thrust_n, 0.0);
    close(
        matched.effective_exhaust_velocity_m_s,
        matched.exit_velocity_m_s,
        1e-14,
    );
    for scale in [0.125, 8.0] {
        let geometry = RocketNozzle {
            throat_area_m2: scale * nozzle().throat_area_m2,
            exit_area_m2: scale * nozzle().exit_area_m2,
        };
        let scaled = reservoir
            .rocket_nozzle_flow(0.0, geometry, controls(), || false)
            .unwrap();
        assert_eq!(scaled.exit(), vacuum.exit());
        close(scaled.mass_flow_kg_s, scale * vacuum.mass_flow_kg_s, 1e-13);
        close(scaled.thrust_n, scale * vacuum.thrust_n, 1e-13);
        close(scaled.specific_impulse_s, vacuum.specific_impulse_s, 1e-13);
    }
}

#[test]
fn unsupported_shocks_source_bounds_and_work_exhaustion_return_no_flow() {
    let reservoir = gas(0.002);
    let before = reservoir.clone();
    let flow = reservoir
        .rocket_nozzle_flow(0.0, nozzle(), controls(), || false)
        .unwrap();
    for pa in [
        -1.0,
        f64::NAN,
        f64::INFINITY,
        1.01 * flow.exit().state().pressure,
    ] {
        assert!(
            reservoir
                .rocket_nozzle_flow(pa, nozzle(), controls(), || false)
                .is_err()
        );
    }
    for geometry in [
        RocketNozzle {
            throat_area_m2: 0.0,
            ..nozzle()
        },
        RocketNozzle {
            exit_area_m2: 0.5 * nozzle().throat_area_m2,
            ..nozzle()
        },
        RocketNozzle {
            throat_area_m2: f64::NAN,
            ..nozzle()
        },
        RocketNozzle {
            exit_area_m2: f64::INFINITY,
            ..nozzle()
        },
    ] {
        assert!(
            reservoir
                .rocket_nozzle_flow(0.0, geometry, controls(), || false)
                .is_err()
        );
    }
    let narrow = RocketNozzleSolve {
        nozzle: GasEnergySolve {
            minimum_temperature_k: 1400.0,
            ..controls().nozzle
        },
        ..controls()
    };
    assert!(
        reservoir
            .rocket_nozzle_flow(0.0, nozzle(), narrow, || false)
            .is_err()
    );
    let limited = RocketNozzleSolve {
        nozzle: GasEnergySolve {
            maximum_iterations: 1,
            ..controls().nozzle
        },
        ..controls()
    };
    assert!(
        reservoir
            .rocket_nozzle_flow(0.0, nozzle(), limited, || false)
            .is_err()
    );
    assert_eq!(reservoir, before);
}

#[test]
fn g4_cancellation_is_polled_through_both_roots_and_before_publication() {
    let reservoir = gas(0.002);
    let mut calls = 0;
    reservoir
        .rocket_nozzle_flow(0.0, nozzle(), controls(), || {
            calls += 1;
            false
        })
        .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut seen = 0;
        assert_eq!(
            reservoir.rocket_nozzle_flow(0.0, nozzle(), controls(), || {
                seen += 1;
                seen == stop
            }),
            Err(GasEnergyError::Cancelled)
        );
    }
}
