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
fn invalid_inputs_source_bounds_and_work_exhaustion_return_no_flow() {
    let reservoir = gas(0.002);
    let before = reservoir.clone();
    for pa in [
        -1.0,
        f64::NAN,
        f64::INFINITY,
        reservoir.state().state().pressure,
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

// Independent calorically-perfect area/Mach and Rankine-Hugoniot oracles.
fn analytic_mach(area: f64, gamma: f64, supersonic: bool) -> f64 {
    let (mut lo, mut hi) = if supersonic { (1.0, 20.0) } else { (1e-8, 1.0) };
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        let residual = area_ratio(mid, gamma) - area;
        if (residual < 0.0) == supersonic {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}
fn perfect_shock(mach: f64, gamma: f64) -> (f64, f64, f64) {
    let m2 = mach * mach;
    (
        (gamma + 1.0) * m2 / ((gamma - 1.0) * m2 + 2.0),
        (2.0 * gamma * m2 - (gamma - 1.0)) / (gamma + 1.0),
        ((1.0 + 0.5 * (gamma - 1.0) * m2) / (gamma * m2 - 0.5 * (gamma - 1.0))).sqrt(),
    )
}

#[test]
fn g1_normal_shock_matches_independent_perfect_gas_jump() {
    let incoming =
        ClosedGasParcel::new(gas(0.0).model().clone(), point(600.0, 200000.0), 0.01).unwrap();
    let state = incoming.state().state();
    let gamma = state.gamma;
    for mach in [1.0001, 1.2, 2.0, 3.0, 5.0] {
        let shock = incoming
            .normal_shock_flow(mach * state.sound_speed, controls().nozzle, || false)
            .unwrap();
        let (rho_ratio, p_ratio, m2) = perfect_shock(mach, gamma);
        close(
            shock.downstream().state().density / state.density,
            rho_ratio,
            2e-8,
        );
        close(
            shock.downstream().state().pressure / state.pressure,
            p_ratio,
            2e-8,
        );
        close(
            shock.downstream().state().temperature / state.temperature,
            p_ratio / rho_ratio,
            2e-8,
        );
        close(shock.downstream_mach_number, m2, 2e-8);
        close(
            shock.downstream_velocity_m_s,
            shock.upstream_velocity_m_s / rho_ratio,
            2e-8,
        );
        let total_pressure_ratio =
            rho_ratio.powf(gamma / (gamma - 1.0)) / p_ratio.powf(1.0 / (gamma - 1.0));
        close(shock.stagnation_pressure_ratio, total_pressure_ratio, 2e-8);
        assert!(shock.relative_mass_flux_residual.abs() < 1e-13);
        assert!(shock.energy_residual_j_kg.abs() < 4e-6);
        assert!(shock.entropy_increase_j_kg_k > -1e-10);
    }
    assert_eq!(incoming.state().state().temperature, 600.0);
}

#[test]
fn g1_variable_cp_shock_conserves_actual_energy_and_increases_entropy() {
    let incoming =
        ClosedGasParcel::new(gas(0.002).model().clone(), point(700.0, 200000.0), 0.01).unwrap();
    let old = incoming.state().state();
    let u1 = 2.5 * old.sound_speed;
    let shock = incoming
        .normal_shock_flow(u1, controls().nozzle, || false)
        .unwrap();
    let end = shock.downstream().state();
    let r = old.specific_gas_constant;
    let dh = r
        * (3.25 * (end.temperature - old.temperature)
            + 0.00075 * (end.temperature - old.temperature) * (end.temperature + old.temperature));
    let ds = r
        * (3.25 * (end.temperature / old.temperature).ln()
            + 0.0015 * (end.temperature - old.temperature)
            - (end.pressure / old.pressure).ln());
    assert!(
        (dh + 0.5 * (shock.downstream_velocity_m_s - u1) * (shock.downstream_velocity_m_s + u1))
            .abs()
            < 5e-6
    );
    close(shock.entropy_increase_j_kg_k, ds, 1e-10);
    assert!(ds > 0.0);
    close(
        old.density * u1,
        end.density * shock.downstream_velocity_m_s,
        1e-13,
    );
    close(
        old.pressure + old.density * u1 * u1,
        end.pressure + end.density * shock.downstream_velocity_m_s.powi(2),
        1e-13,
    );
    let (_, constant_gamma_pressure, _) = perfect_shock(2.5, old.gamma);
    assert!((end.pressure / old.pressure - constant_gamma_pressure).abs() > 0.01);
}

#[test]
fn g1_internal_nozzle_shock_recovers_independently_specified_area() {
    let chamber = gas(0.0);
    let state = chamber.state().state();
    let gamma = state.gamma;
    let alpha = 0.5 * (gamma - 1.0);
    let exponent = gamma / (gamma - 1.0);
    for shock_area_ratio in [1.3, 2.0, 3.7] {
        let m1 = analytic_mach(shock_area_ratio, gamma, true);
        let (_, pressure_jump, m2) = perfect_shock(m1, gamma);
        let p1 = state.pressure / (1.0 + alpha * m1 * m1).powf(exponent);
        let p2 = p1 * pressure_jump;
        let me = analytic_mach(area_ratio(m2, gamma) * 4.0 / shock_area_ratio, gamma, false);
        let ambient = p2 * ((1.0 + alpha * m2 * m2) / (1.0 + alpha * me * me)).powf(exponent);
        let flow = chamber
            .rocket_nozzle_flow(ambient, nozzle(), controls(), || false)
            .unwrap();
        assert_eq!(flow.regime, RocketNozzleRegime::ChokedNormalShock);
        let shock = flow.normal_shock.as_ref().unwrap();
        close(
            shock.area_m2 / nozzle().throat_area_m2,
            shock_area_ratio,
            2e-8,
        );
        close(shock.flow.upstream_mach_number, m1, 2e-8);
        close(shock.flow.downstream_mach_number, m2, 2e-8);
        close(flow.exit_mach_number, me, 2e-8);
        close(flow.exit().state().pressure, ambient, 2e-11);
        close(
            flow.thrust_n,
            flow.mass_flow_kg_s * flow.exit_velocity_m_s,
            2e-9,
        );
        assert!(flow.entropy_increase_j_kg_k > 0.0);
        assert!(flow.energy_residual_j_kg.abs() < 2e-6);
        assert!(flow.entropy_energy_residual_j_kg.abs() < 4e-6);
    }
}

#[test]
fn g1_unchoked_flow_and_external_compression_have_distinct_physical_closures() {
    let chamber = gas(0.0);
    let state = chamber.state().state();
    let gamma = state.gamma;
    let mt = 0.6;
    let me = analytic_mach(4.0 * area_ratio(mt, gamma), gamma, false);
    let ambient =
        state.pressure / (1.0 + 0.5 * (gamma - 1.0) * me * me).powf(gamma / (gamma - 1.0));
    let unchoked = chamber
        .rocket_nozzle_flow(ambient, nozzle(), controls(), || false)
        .unwrap();
    assert_eq!(unchoked.regime, RocketNozzleRegime::UnchokedSubsonic);
    close(unchoked.throat_mach_number, mt, 2e-8);
    close(unchoked.exit_mach_number, me, 2e-8);
    assert!(unchoked.normal_shock.is_none());
    assert!(unchoked.back_pressure_residual_pa.unwrap().abs() < 3e-5);
    let vacuum = chamber
        .rocket_nozzle_flow(0.0, nozzle(), controls(), || false)
        .unwrap();
    assert!(unchoked.mass_flow_kg_s < vacuum.mass_flow_kg_s);
    let (_, jump, _) = perfect_shock(vacuum.exit_mach_number, gamma);
    let pa = 0.5 * vacuum.exit().state().pressure * (1.0 + jump);
    let external = chamber
        .rocket_nozzle_flow(pa, nozzle(), controls(), || false)
        .unwrap();
    assert_eq!(
        external.regime,
        RocketNozzleRegime::ChokedExternalCompression
    );
    assert!(external.normal_shock.is_none());
    assert!(external.back_pressure_residual_pa.is_none());
    assert!(external.pressure_thrust_n < 0.0);
    assert_eq!(external.exit(), vacuum.exit());
    close(
        vacuum.thrust_n - external.thrust_n,
        pa * nozzle().exit_area_m2,
        1e-12,
    );
}

#[test]
fn g3_thrust_is_continuous_at_shock_exit_and_unchoking_thresholds() {
    let chamber = gas(0.0);
    let gamma = chamber.state().state().gamma;
    let vacuum = chamber
        .rocket_nozzle_flow(0.0, nozzle(), controls(), || false)
        .unwrap();
    let (_, pressure_jump, _) = perfect_shock(vacuum.exit_mach_number, gamma);
    let exit_threshold = vacuum.exit().state().pressure * pressure_jump;
    let me = analytic_mach(4.0, gamma, false);
    let choke_threshold = chamber.state().state().pressure
        / (1.0 + 0.5 * (gamma - 1.0) * me * me).powf(gamma / (gamma - 1.0));
    for threshold in [exit_threshold, choke_threshold] {
        let below = chamber
            .rocket_nozzle_flow(threshold * (1.0 - 1e-7), nozzle(), controls(), || false)
            .unwrap();
        let above = chamber
            .rocket_nozzle_flow(threshold * (1.0 + 1e-7), nozzle(), controls(), || false)
            .unwrap();
        close(below.thrust_n, above.thrust_n, 1e-4);
        close(below.mass_flow_kg_s, above.mass_flow_kg_s, 1e-5);
    }
    // At a shock exactly on the exit plane, momentum flux conservation makes
    // upstream and downstream thrust identical, including pressure thrust once.
    let raw_exit = ClosedGasParcel::new(
        chamber.model().clone(),
        point(
            vacuum.exit().state().temperature,
            vacuum.exit().state().pressure,
        ),
        0.01,
    )
    .unwrap();
    let shock = raw_exit
        .normal_shock_flow(vacuum.exit_velocity_m_s, controls().nozzle, || false)
        .unwrap();
    let pa = shock.downstream().state().pressure;
    let upstream_thrust = vacuum.mass_flow_kg_s * vacuum.exit_velocity_m_s
        + (vacuum.exit().state().pressure - pa) * nozzle().exit_area_m2;
    let downstream_thrust = vacuum.mass_flow_kg_s * shock.downstream_velocity_m_s;
    close(upstream_thrust, downstream_thrust, 1e-10);
}

#[test]
fn g3_subsonic_source_admission_does_not_demand_hypothetical_sonic_transport() {
    let chamber = ClosedGasParcel::new(
        model_in_domains(
            vec![0.75, 0.25],
            0.0,
            [200.0, 4000.0],
            [200.0, 4000.0],
            [2.7e6, 3.1e6],
        ),
        point(1800.0, 3e6),
        0.01,
    )
    .unwrap();
    let flow = chamber
        .rocket_nozzle_flow(2.999e6, nozzle(), controls(), || false)
        .unwrap();
    assert_eq!(flow.regime, RocketNozzleRegime::UnchokedSubsonic);
    assert!(flow.throat().state().pressure > 2.7e6);
    assert!(
        chamber
            .rocket_nozzle_flow(0.0, nozzle(), controls(), || false)
            .is_err()
    );
}

#[test]
fn g3_subsonic_flow_does_not_require_a_source_valid_hypothetical_sonic_temperature() {
    let chamber = ClosedGasParcel::new(
        model_in_domains(
            vec![0.75, 0.25],
            0.0,
            [1700.0, 4000.0],
            [1700.0, 4000.0],
            [2.7e6, 3.1e6],
        ),
        point(1800.0, 3e6),
        0.01,
    )
    .unwrap();
    let limited = RocketNozzleSolve {
        nozzle: GasEnergySolve {
            minimum_temperature_k: 1700.0,
            ..controls().nozzle
        },
        ..controls()
    };
    let flow = chamber
        .rocket_nozzle_flow(2.999e6, nozzle(), limited, || false)
        .unwrap();
    assert_eq!(flow.regime, RocketNozzleRegime::UnchokedSubsonic);
    assert!(flow.throat().state().temperature > 1700.0);
    assert!(flow.exit().state().temperature > 1700.0);
    let wide = gas(0.0)
        .rocket_nozzle_flow(2.999e6, nozzle(), controls(), || false)
        .unwrap();
    close(flow.mass_flow_kg_s, wide.mass_flow_kg_s, 2e-8);
    close(flow.thrust_n, wide.thrust_n, 2e-8);
    assert!(
        chamber
            .rocket_nozzle_flow(0.0, nozzle(), limited, || false)
            .is_err()
    );
}

#[test]
fn g4_internal_shock_and_reusable_jump_cancel_without_publishing() {
    let chamber = gas(0.002);
    let initial = chamber.clone();
    let mut polls = 0;
    let flow = chamber
        .rocket_nozzle_flow(
            0.8 * chamber.state().state().pressure,
            nozzle(),
            controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert_eq!(flow.regime, RocketNozzleRegime::ChokedNormalShock);
    for at in [1, polls / 2, polls] {
        let mut calls = 0;
        assert_eq!(
            chamber.rocket_nozzle_flow(
                0.8 * chamber.state().state().pressure,
                nozzle(),
                controls(),
                || {
                    calls += 1;
                    calls == at
                }
            ),
            Err(GasEnergyError::Cancelled)
        );
    }
    let shock = flow.normal_shock.as_ref().unwrap();
    let incoming = ClosedGasParcel::new(
        chamber.model().clone(),
        point(
            shock.flow.upstream().state().temperature,
            shock.flow.upstream().state().pressure,
        ),
        0.01,
    )
    .unwrap();
    assert_eq!(
        incoming.normal_shock_flow(shock.flow.upstream_velocity_m_s, controls().nozzle, || true),
        Err(GasEnergyError::Cancelled)
    );
    assert!(
        incoming
            .normal_shock_flow(
                0.9 * incoming.state().state().sound_speed,
                controls().nozzle,
                || false
            )
            .is_err()
    );
    assert_eq!(chamber, initial);
}

#[test]
fn finite_blowdown_uses_the_same_internal_shock_regime() {
    use super::super::blowdown::{GasBlowdownSolve, RocketBlowdownSolve};
    use crate::thermochemical_gas::energy::entropy::GasEntropySolve;
    let chamber = gas(0.002);
    let step = chamber
        .rocket_blowdown(
            0.8 * chamber.state().state().pressure,
            nozzle(),
            1e-4,
            RocketBlowdownSolve {
                nozzle: controls(),
                reservoir: GasBlowdownSolve {
                    entropy: GasEntropySolve {
                        minimum_temperature_k: 200.0,
                        maximum_temperature_k: 4000.0,
                        entropy_tolerance_j_kg_k: 1e-9,
                        maximum_iterations: 128,
                    },
                    mass_tolerance_kg: 1e-11,
                },
                maximum_donor_fraction: 0.05,
                maximum_iterations: 128,
            },
            || false,
        )
        .unwrap();
    assert_eq!(
        step.midpoint_flow().regime,
        RocketNozzleRegime::ChokedNormalShock
    );
    assert!(step.reservoir().endpoint().mass_kg() < chamber.mass_kg());
    assert!(step.axial_impulse_n_s > 0.0);
    close(
        step.axial_impulse_n_s,
        1e-4 * step.midpoint_flow().thrust_n,
        1e-13,
    );
    assert!(step.mass_quadrature_residual_kg.abs() <= 1e-11);
}
