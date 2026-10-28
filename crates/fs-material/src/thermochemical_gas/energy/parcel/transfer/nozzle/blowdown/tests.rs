//! G0 inventory/source invariants and G1 analytic adiabatic blowdown/refinement.

use super::*;
use crate::thermochemical_gas::energy::GasEnergySolve;
use crate::thermochemical_gas::tests::{model, point};

fn controls() -> GasBlowdownSolve {
    GasBlowdownSolve {
        entropy: GasEntropySolve {
            minimum_temperature_k: 200.0,
            maximum_temperature_k: 1800.0,
            entropy_tolerance_j_kg_k: 1e-11,
            maximum_iterations: 128,
        },
        mass_tolerance_kg: 1e-14,
    }
}
fn rocket_controls() -> RocketBlowdownSolve {
    RocketBlowdownSolve {
        reservoir: controls(),
        maximum_donor_fraction: 0.15,
        maximum_iterations: 32,
        nozzle: RocketNozzleSolve {
            nozzle: GasEnergySolve {
                minimum_temperature_k: 200.0,
                maximum_temperature_k: 1800.0,
                energy_tolerance_j_kg: 1e-7,
                maximum_iterations: 128,
            },
            relative_mass_flow_tolerance: 1e-10,
        },
    }
}
fn nozzle() -> RocketNozzle {
    RocketNozzle {
        throat_area_m2: 1e-6,
        exit_area_m2: 2e-6,
    }
}
fn gas(slope: f64) -> ClosedGasParcel {
    ClosedGasParcel::new(
        model(vec![0.75, 0.25], slope),
        point(1000.0, 1000000.0),
        0.001,
    )
    .unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * b.abs().max(1e-12),
        "{a} != {b}"
    );
}

#[test]
fn finite_withdrawal_matches_exact_perfect_gas_depletion_and_inventory_energy() {
    let initial = gas(0.0);
    let before = initial.clone();
    let gamma = initial.state().state().gamma;
    for fraction in [0.01, 0.2, 0.6] {
        let mass = fraction * initial.mass_kg();
        let candidate = initial
            .exhaust_mass_adiabatic(mass, controls(), || false)
            .unwrap();
        let state = candidate.endpoint().state().state();
        near(
            state.temperature,
            1000.0 * (1.0 - fraction).powf(gamma - 1.0),
            1e-11,
        );
        near(
            state.pressure,
            1000000.0 * (1.0 - fraction).powf(gamma),
            1e-11,
        );
        near(
            candidate.endpoint().mass_kg(),
            initial.mass_kg() - mass,
            1e-14,
        );
        assert_eq!(candidate.endpoint().model(), initial.model());
        assert_eq!(candidate.endpoint().volume_m3(), initial.volume_m3());
        let lost = initial.mass_kg() * initial.state().specific_internal_energy_j_kg()
            - candidate.endpoint().mass_kg()
                * candidate.endpoint().state().specific_internal_energy_j_kg();
        near(candidate.step().exhausted_enthalpy_j, lost, 1e-11);
        near(
            candidate.step().species.iter().map(|s| s.mass_kg).sum(),
            mass,
            1e-14,
        );
        assert!(candidate.step().mass_residual_kg.abs() <= controls().mass_tolerance_kg);
        assert!(candidate.step().maximum_species_mass_residual_kg <= controls().mass_tolerance_kg);
        assert!(
            candidate.step().specific_entropy_residual_j_kg_k.abs()
                <= controls().entropy.entropy_tolerance_j_kg_k
        );
        let mean = candidate.step().mean_stagnation_enthalpy_j_kg.unwrap();
        assert!(mean > candidate.endpoint().state().specific_enthalpy_j_kg());
        assert!(mean < initial.state().specific_enthalpy_j_kg());
    }
    assert_eq!(initial, before);
}

#[test]
fn variable_heat_capacity_follows_source_entropy_instead_of_initial_gamma() {
    let initial = gas(0.002);
    let fraction = 0.4;
    let candidate = initial
        .exhaust_mass_adiabatic(fraction * initial.mass_kg(), controls(), || false)
        .unwrap();
    let actual = candidate.endpoint().state().state().temperature;
    // cp/R = 3.25 + 0.0015*T for the fixture's mixture, so constant-volume
    // reservoir entropy obeys 2.25 ln(T/T0)+0.0015(T-T0)=ln(m/m0).
    let (mut lo, mut hi) = (200.0, 1000.0);
    for _ in 0..80 {
        let t: f64 = 0.5 * (lo + hi);
        let residual = 2.25 * (t / 1000.0).ln() + 0.0015 * (t - 1000.0) - (1.0_f64 - fraction).ln();
        if residual < 0.0 {
            lo = t;
        } else {
            hi = t;
        }
    }
    near(actual, 0.5 * (lo + hi), 1e-11);
    let frozen_gamma = 1000.0 * (1.0 - fraction).powf(initial.state().state().gamma - 1.0);
    assert!((actual - frozen_gamma).abs() > 1.0);
    let expected_du = initial.state().state().specific_gas_constant
        * (2.25 * (actual - 1000.0) + 0.00075 * (actual * actual - 1000.0 * 1000.0));
    let actual_loss = candidate.step().exhausted_enthalpy_j;
    let expected_loss =
        fraction * initial.mass_kg() * initial.state().specific_internal_energy_j_kg()
            - candidate.endpoint().mass_kg() * expected_du;
    near(actual_loss, expected_loss, 1e-11);
}

#[test]
fn vacuum_transient_converges_second_order_to_analytic_mass_pressure_and_impulse() {
    let initial = gas(0.0);
    let flow = initial
        .rocket_nozzle_flow(0.0, nozzle(), rocket_controls().nozzle, || false)
        .unwrap();
    let n = 0.5 * (initial.state().state().gamma - 1.0);
    let duration = 0.3;
    let ratio = (1.0 + n * (flow.mass_flow_kg_s / initial.mass_kg()) * duration).powf(-1.0 / n);
    let exact_mass = initial.mass_kg() * ratio;
    let exact_pressure =
        initial.state().state().pressure * ratio.powf(initial.state().state().gamma);
    let exact_impulse = flow.effective_exhaust_velocity_m_s * initial.mass_kg() / (1.0 + n)
        * (1.0 - ratio.powf(1.0 + n));
    let mut mass_errors = Vec::new();
    let mut impulse_errors = Vec::new();
    for steps in [4, 8, 16] {
        let mut reservoir = initial.clone();
        let mut impulse = 0.0;
        let mut exhaust_energy = 0.0;
        for _ in 0..steps {
            let candidate = reservoir
                .rocket_blowdown(
                    0.0,
                    nozzle(),
                    duration / f64::from(steps),
                    rocket_controls(),
                    || false,
                )
                .unwrap();
            assert!(candidate.mass_quadrature_residual_kg.abs() <= controls().mass_tolerance_kg);
            impulse += candidate.axial_impulse_n_s;
            let (endpoint, step, _) = candidate.into_parts();
            assert!(endpoint.mass_kg() > 0.0);
            exhaust_energy += step.exhausted_enthalpy_j;
            reservoir = endpoint;
        }
        mass_errors.push((reservoir.mass_kg() - exact_mass).abs());
        impulse_errors.push((impulse - exact_impulse).abs());
        let actual_energy_loss = initial.mass_kg()
            * initial.state().specific_internal_energy_j_kg()
            - reservoir.mass_kg() * reservoir.state().specific_internal_energy_j_kg();
        near(exhaust_energy, actual_energy_loss, 1e-10);
        near(reservoir.state().state().pressure, exact_pressure, 3e-4);
    }
    for errors in [&mass_errors, &impulse_errors] {
        assert!(errors[0] > 3.8 * errors[1], "{errors:?}");
        assert!(errors[1] > 3.8 * errors[2], "{errors:?}");
    }
}

#[test]
fn zero_withdrawal_preserves_the_full_source_bearing_state() {
    let initial = gas(0.001);
    let candidate = initial
        .exhaust_mass_adiabatic(0.0, controls(), || false)
        .unwrap();
    assert_eq!(candidate.endpoint(), &initial);
    assert!(candidate.step().species.is_empty());
    assert_eq!(candidate.step().exhausted_enthalpy_j, 0.0);
    assert_eq!(candidate.step().mean_stagnation_enthalpy_j_kg, None);
}

#[test]
fn invalid_inventory_source_bounds_iteration_limits_and_cancellation_never_publish() {
    let initial = gas(0.001);
    let before = initial.clone();
    for amount in [
        -1.0,
        f64::NAN,
        f64::INFINITY,
        initial.mass_kg(),
        2.0 * initial.mass_kg(),
    ] {
        assert!(
            initial
                .exhaust_mass_adiabatic(amount, controls(), || false)
                .is_err()
        );
    }
    let narrow = GasBlowdownSolve {
        entropy: GasEntropySolve {
            minimum_temperature_k: 990.0,
            ..controls().entropy
        },
        ..controls()
    };
    assert!(
        initial
            .exhaust_mass_adiabatic(0.5 * initial.mass_kg(), narrow, || false)
            .is_err()
    );
    assert!(
        initial
            .rocket_blowdown(0.0, nozzle(), 1000.0, rocket_controls(), || false)
            .is_err()
    );
    assert!(matches!(
        initial.rocket_blowdown(
            0.0,
            nozzle(),
            0.03,
            RocketBlowdownSolve {
                maximum_iterations: 1,
                ..rocket_controls()
            },
            || false
        ),
        Err(GasBlowdownError::NotConverged { .. })
    ));
    for dt in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(
            initial
                .rocket_blowdown(0.0, nozzle(), dt, rocket_controls(), || false)
                .is_err()
        );
    }
    let mut total = 0;
    initial
        .rocket_blowdown(0.0, nozzle(), 0.01, rocket_controls(), || {
            total += 1;
            false
        })
        .unwrap();
    for stop in [1, total / 2, total] {
        let mut count = 0;
        assert!(
            initial
                .rocket_blowdown(0.0, nozzle(), 0.01, rocket_controls(), || {
                    count += 1;
                    count == stop
                })
                .is_err()
        );
        assert_eq!(initial, before);
    }
}

#[test]
fn endpoint_can_cross_into_supported_external_compression_without_losing_impulse() {
    use super::super::rocket::RocketNozzleRegime;
    let initial = gas(0.0);
    let initial_flow = initial
        .rocket_nozzle_flow(0.0, nozzle(), rocket_controls().nozzle, || false)
        .unwrap();
    let ambient = 0.98 * initial_flow.exit().state().pressure;
    assert!(
        initial
            .rocket_nozzle_flow(ambient, nozzle(), rocket_controls().nozzle, || false)
            .is_ok()
    );
    let candidate = initial
        .rocket_blowdown(ambient, nozzle(), 0.3, rocket_controls(), || false)
        .unwrap();
    let endpoint = candidate.reservoir().endpoint();
    let flow = endpoint
        .rocket_nozzle_flow(ambient, nozzle(), rocket_controls().nozzle, || false)
        .unwrap();
    assert_eq!(flow.regime, RocketNozzleRegime::ChokedExternalCompression);
    assert!(flow.pressure_thrust_n < 0.0);
    assert!(endpoint.mass_kg() < initial.mass_kg());
    assert!(candidate.reservoir().step().mass_residual_kg.abs() <= controls().mass_tolerance_kg);
    near(
        candidate.axial_impulse_n_s,
        0.3 * candidate.midpoint_flow().thrust_n,
        1e-13,
    );
}
