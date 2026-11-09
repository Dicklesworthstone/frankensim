use super::*;
use crate::multicomponent::Conserved;
use crate::multicomponent::tests::{gas, source};
use crate::multicomponent::wall_exchange::WallThermal;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn model() -> MixtureGas {
    // Original NASA calorics: R_specific=1, Cv=2.5, admitted T in .001..=50.
    gas(vec![source("A", R, 3.5, 0.0, 0.0)])
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn controls() -> SourceCoupledDragControl {
    SourceCoupledDragControl {
        absolute_slip_tolerance_m_s: 1e-11,
        relative_slip_tolerance: 1e-8,
        absolute_heating_tolerance_j: 1e-12,
        relative_heating_tolerance: 1e-8,
        ..SourceCoupledDragControl::default()
    }
}

fn contact(wall: f64) -> WallContact {
    WallContact {
        drag: WallDrag::MovingLinear {
            rate_per_s: 1.0,
            wall_velocity_m_s: wall,
        },
        thermal: WallThermal::Adiabatic,
    }
}

fn run(
    g: &MixtureGas,
    state: &Conserved,
    wall: f64,
    dt: f64,
    source_control: SourceCoupledDragControl,
    mut rate: impl FnMut(f64, f64) -> f64,
) -> SourceCoupledDragUpdate {
    exchange_cell_with_source_drag(
        g,
        state.as_cell(),
        0.2,
        Time::new(dt),
        contact(wall),
        WallControl::default(),
        source_control,
        |trial, cancelled| {
            let primitive = g.primitive(trial, cancelled)?;
            Ok(WallDrag::MovingLinear {
                rate_per_s: rate(primitive.temperature, primitive.velocity - wall),
                wall_velocity_m_s: wall,
            })
        },
        || false,
    )
    .unwrap()
}

#[test]
fn unchanged_original_source_retains_stiff_exact_decay_and_work() {
    let g = model();
    for wall in [0.0, 0.7, -0.3] {
        for signed_slip in [-0.4, 0.4] {
            let state = g.conserved(&[1.0], wall + signed_slip, 2.0).unwrap();
            for rate in [0.0, 4.0, 1e6] {
                let update = run(&g, &state, wall, 0.3, controls(), |_, _| rate);
                let expected = signed_slip * (-rate * 0.3).exp();
                close(update.gas.momentum, wall + expected, 2e-15);
                close(
                    update.gas.report.drag_heating_j,
                    0.1 * (signed_slip * signed_slip - expected * expected),
                    2e-16,
                );
                let work = wall * (update.gas.momentum - state.momentum) * 0.2;
                close(update.gas.report.mechanical_work_to_gas_j, work, 1e-16);
                close((update.gas.energy - state.energy) * 0.2, work, 2e-16);
                assert_eq!(update.source.accepted_steps, 1);
                assert_eq!(update.source.attempted_steps, 1);
            }
        }
    }
}

#[test]
fn viscosity_proportional_to_actual_temperature_matches_exact_thermal_feedback() {
    let g = model();
    let (t0, cv, slip0, a0, dt): (f64, f64, f64, f64, f64) = (2.0, 2.5, 0.4, 4.0, 0.3);
    // The reciprocal of squared slip obeys a linear ODE. This is an independent
    // continuum oracle, not a repeated evaluation of the integration scheme.
    let delta = slip0 * slip0 / (2.0 * cv * t0);
    let decay = (-a0 * (1.0 + delta) * dt).exp();
    let expected = slip0 * decay * ((1.0 + delta) / (1.0 + delta * decay * decay)).sqrt();
    for wall in [0.0, 0.7, -2.0] {
        for sign in [-1.0, 1.0] {
            let state = g.conserved(&[1.0], wall + sign * slip0, t0).unwrap();
            let update = run(&g, &state, wall, dt, controls(), |t, _| a0 * t / t0);
            close(update.gas.momentum - wall, sign * expected, 5e-9);
            close(
                update.gas.report.drag_heating_j,
                0.1 * (slip0 * slip0 - expected * expected),
                2e-10,
            );
            let endpoint = Cell {
                partial_densities: &state.partial_densities,
                momentum: update.gas.momentum,
                energy: update.gas.energy,
            };
            close(
                g.primitive(endpoint, || false).unwrap().temperature,
                t0 + (slip0 * slip0 - expected * expected) / (2.0 * cv),
                2e-9,
            );
            assert!(update.source.accepted_steps > 1);
            assert!(
                update.source.estimated_slip_error_m_s
                    <= controls().absolute_slip_tolerance_m_s
                        + controls().relative_slip_tolerance * slip0
            );
            assert!(update.gas.report.energy_residual_j.abs() < 2e-15);
        }
    }
    // The original entering-temperature freeze has a resolved physical error.
    assert!((slip0 * (-a0 * dt).exp() - expected).abs() > 1e-3);
}

#[test]
fn midpoint_ordering_prevents_stiff_decreasing_viscosity_false_agreement() {
    let g = model();
    let (slip0, low, high, q, dt): (f64, f64, f64, i32, f64) = (0.4, 0.001, 1000.0, 10, 1.0);
    let state = g.conserved(&[1.0], slip0, 2.0).unwrap();
    let control = SourceCoupledDragControl::default();
    let update = run(&g, &state, 0.0, dt, control, |_, slip| {
        low + (high - low) * (slip / slip0).powi(q)
    });
    // A smooth positive Bernoulli law. Coarse and two-half exponential
    // midpoints without the ordering gate both return s0*exp(-low*dt), hiding
    // the entire fast initial transient with an exactly zero embedded defect.
    let expected = slip0 * (-low * dt).exp()
        / (1.0 + ((high - low) / low) * -(-(q as f64) * low * dt).exp_m1()).powf(1.0 / q as f64);
    close(update.gas.momentum, expected, 1e-8);
    assert!((slip0 * (-low * dt).exp() - expected).abs() > 0.2);
    assert!(update.source.attempted_steps > update.source.accepted_steps);
    assert!(
        update.source.estimated_heating_error_j
            <= control.absolute_heating_tolerance_j
                + control.relative_heating_tolerance * update.gas.report.drag_heating_j
    );
}

#[test]
fn source_refusal_cancellation_and_whole_interval_work_are_transactional() {
    let g = model();
    let state = g.conserved(&[1.0], 0.4, 2.0).unwrap();
    let original = state.clone();
    let refused = WallExchangeError::TransportDomain {
        quantity: "manufactured viscosity temperature",
        value: 2.01,
        minimum: 0.001,
        maximum: 2.005,
    };
    let result = exchange_cell_with_source_drag(
        &g,
        state.as_cell(),
        0.2,
        Time::new(0.3),
        contact(0.0),
        WallControl::default(),
        controls(),
        |trial, cancelled| {
            let t = g.primitive(trial, cancelled)?.temperature;
            if t > 2.005 {
                return Err(refused.clone());
            }
            Ok(WallDrag::Linear { rate_per_s: 4.0 })
        },
        || false,
    );
    assert_eq!(result, Err(refused));
    let mut calls = 0;
    let exhausted = exchange_cell_with_source_drag(
        &g,
        state.as_cell(),
        0.2,
        Time::new(0.3),
        contact(0.0),
        WallControl::default(),
        SourceCoupledDragControl {
            maximum_source_evaluations: 1,
            ..controls()
        },
        |_, _| {
            calls += 1;
            Ok(WallDrag::Linear { rate_per_s: 4.0 })
        },
        || false,
    );
    assert!(matches!(
        exhausted,
        Err(WallExchangeError::SourceDragWork {
            source_evaluations: 1,
            ..
        })
    ));
    assert_eq!(calls, 1);
    let cancelled = exchange_cell_with_source_drag(
        &g,
        state.as_cell(),
        0.2,
        Time::new(0.3),
        contact(0.0),
        WallControl::default(),
        controls(),
        |_, _| panic!("cancelled source must not be evaluated"),
        || true,
    );
    assert_eq!(cancelled, Err(WallExchangeError::Cancelled));
    assert_eq!(state.partial_densities, original.partial_densities);
    assert_eq!(state.momentum, original.momentum);
    assert_eq!(state.energy, original.energy);
}

#[test]
fn identity_and_weak_drive_keep_resolved_impulse_without_sensible_double_counting() {
    let g = model();
    let at_rest = g.conserved(&[1.0], 0.7, 2.0).unwrap();
    let identity = run(&g, &at_rest, 0.7, 1.0, controls(), |_, _| 1e6);
    assert_eq!(identity.gas.momentum, at_rest.momentum);
    assert_eq!(identity.gas.energy, at_rest.energy);
    assert_eq!(identity.source.attempted_steps, 0);
    let state = g.conserved(&[1.0], 0.0, 2.0).unwrap();
    let weak = run(&g, &state, 1.0, 1.0, controls(), |_, _| 1e-18);
    close(weak.gas.momentum, 1e-18, 2e-33);
    close(weak.gas.report.drag_heating_j, 2e-19, 2e-34);
    assert!(weak.gas.momentum > 0.0);
    assert_eq!(weak.gas.report.heat_to_gas_j, 0.0);
}
