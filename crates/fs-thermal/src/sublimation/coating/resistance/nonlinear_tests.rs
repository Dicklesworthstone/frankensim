//! Independent analytic heat laws on the original two-temperature endpoint.
use super::*;
use crate::sublimation::{MOLAR_GAS_CONSTANT, SublimationProperties};

const CS: f64 = 500.0;
const CD: f64 = 2000.0;
const MS: f64 = 0.02;
const MD: f64 = 0.001;
const LATENT: f64 = 2.8e6;
const MOLAR: f64 = 0.018;

struct Phase;
impl SublimationLaw for Phase {
    fn temperature_range_k(&self) -> [f64; 2] {
        [100.0, 500.0]
    }
    fn properties(&self, t: f64) -> Result<SublimationProperties, SublimationError> {
        // A manufactured saturation curve makes the gross outgoing rate
        // temperature-independent; this supplies an independent vapor ODE.
        let pressure = 10.0 * (t / 230.0).sqrt();
        Ok(SublimationProperties {
            saturation_pressure_pa: pressure,
            saturation_pressure_derivative_pa_k: pressure / (2.0 * t),
            latent_heat_j_kg: LATENT,
            latent_heat_derivative_j_kg_k: 0.0,
            solid_density_kg_m3: 900.0,
            molar_mass_kg_mol: MOLAR,
        })
    }
}
fn state(ts: f64, td: f64) -> ResistiveSurfaceState {
    ResistiveSurfaceState {
        substrate_mass_kg: MS,
        substrate_specific_enthalpy_j_kg: CS * ts,
        deposit: SolidSurfaceState {
            mass_kg: MD,
            specific_enthalpy_j_kg: CD * td,
            area_m2: 0.001,
        },
    }
}
fn controls() -> SolverConfig {
    SolverConfig {
        relative_tolerance: 1e-12,
        absolute_energy_tolerance_j: 1e-9,
        max_iterations: 256,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}
fn heat(
    initial: ResistiveSurfaceState,
    dt: f64,
    law: &dyn ResistiveInterfaceLaw,
    qs: f64,
    qd: f64,
) -> ResistiveHeatStep {
    heat_resistive_surface_with_interface(
        &Phase,
        initial,
        dt,
        law,
        qs,
        qd,
        controls(),
        |h| Ok(h / CS),
        |h| Ok(h / CD),
        || Ok(()),
    )
    .unwrap()
}
fn balances(step: ResistiveHeatStep, qs: f64, qd: f64) {
    close(
        step.substrate_energy_change_j + step.deposit_energy_change_j,
        qs + qd,
        5e-9,
    );
    close(step.substrate_energy_residual_j, 0.0, 5e-9);
    close(step.deposit_energy_residual_j, 0.0, 5e-9);
    close(step.energy_residual_j, 0.0, 5e-9);
}

// For k(T) proportional to T, Kirchhoff heat flow is beta*(Ts^2-Td^2).
// Eliminate endpoint temperatures using the TWO independent energy balances,
// then solve the resulting quadratic analytically. No production root is used.
fn quadratic_heat(ts: f64, td: f64, dt: f64, beta: f64, qs: f64, qd: f64) -> f64 {
    let s = ts + qs / (MS * CS);
    let d = td + qd / (MD * CD);
    let a = 1.0 / (MS * CS) + 1.0 / (MD * CD);
    let b = 1.0 / (MD * CD) - 1.0 / (MS * CS);
    let aa = dt * beta * a * b;
    let bb = 1.0 - dt * beta * ((s - d) * b - a * (s + d));
    let cc = -dt * beta * (s - d) * (s + d);
    -2.0 * cc / (bb + (bb * bb - 4.0 * aa * cc).sqrt())
}

#[test]
fn kirchhoff_heat_matches_independent_quadratic_with_separate_caloric_datums() {
    let beta = 0.002;
    let interface = |ts: f64, td: f64| Ok(beta * (ts - td) * (ts + td));
    let (dt, qs, qd) = (0.5, 7.0, -3.0);
    let q = quadratic_heat(300.0, 200.0, dt, beta, qs, qd);
    for (support_datum, deposit_datum) in [(0.0, 0.0), (1e6, -8e5)] {
        let mut initial = state(300.0, 200.0);
        initial.substrate_specific_enthalpy_j_kg += support_datum;
        initial.deposit.specific_enthalpy_j_kg += deposit_datum;
        let step = heat_resistive_surface_with_interface(
            &Phase,
            initial,
            dt,
            &interface,
            qs,
            qd,
            controls(),
            |h| Ok((h - support_datum) / CS),
            |h| Ok((h - deposit_datum) / CD),
            || Ok(()),
        )
        .unwrap();
        close(step.interface_heat_j, q, 5e-9);
        close(
            step.substrate_temperature_k,
            300.0 + (qs - q) / (MS * CS),
            1e-8,
        );
        close(
            step.deposit_temperature_k,
            200.0 + (qd + q) / (MD * CD),
            1e-8,
        );
        assert!(step.interface_heat_j > 0.0);
        balances(step, qs, qd);
    }
}

#[test]
fn fourth_power_exchange_retains_original_nonlinear_caloric_storage() {
    // Independent manufactured endpoints on genuinely nonlinear caloric charts.
    let hs = |t: f64| 300.0 * t + 0.5 * t * t;
    let hd = |t: f64| 1000.0 * t + 2.0 * t * t;
    let ts: [f64; 2] = [360.0, 320.0];
    let td: [f64; 2] = [180.0, 240.0];
    let (dt, radiation) = (0.2, 2e-8);
    let q = dt * radiation * (ts[1].powi(4) - td[1].powi(4));
    let qs = MS * (hs(ts[1]) - hs(ts[0])) + q;
    let qd = MD * (hd(td[1]) - hd(td[0])) - q;
    let interface = |ts: f64, td: f64| Ok(radiation * (ts - td) * (ts + td) * (ts * ts + td * td));
    let initial = ResistiveSurfaceState {
        substrate_specific_enthalpy_j_kg: hs(ts[0]),
        deposit: SolidSurfaceState {
            specific_enthalpy_j_kg: hd(td[0]),
            ..state(360.0, 180.0).deposit
        },
        ..state(360.0, 180.0)
    };
    let step = heat_resistive_surface_with_interface(
        &Phase,
        initial,
        dt,
        &interface,
        qs,
        qd,
        controls(),
        |h| Ok(2.0 * h / (300.0 + (90000.0 + 2.0 * h).sqrt())),
        |h| Ok(2.0 * h / (1000.0 + (1e6 + 8.0 * h).sqrt())),
        || Ok(()),
    )
    .unwrap();
    close(step.substrate_temperature_k, ts[1], 1e-8);
    close(step.deposit_temperature_k, td[1], 1e-8);
    close(step.interface_heat_j, q, 1e-8);
    balances(step, qs, qd);
}

#[test]
fn finite_vapor_phase_exchange_keeps_the_same_latent_and_nonlinear_heat_accounts() {
    let beta = 0.001;
    let interface = |ts: f64, td: f64| Ok(beta * (ts - td) * (ts + td));
    let initial = state(230.0, 220.0);
    let r = MOLAR_GAS_CONSTANT / MOLAR;
    let (bath_t, volume, dt, alpha) = (230.0, 0.01, 0.2, 0.2);
    let b = alpha / (2.0 * core::f64::consts::PI * r * bath_t).sqrt();
    let decay = initial.deposit.area_m2 * b * r * bath_t / volume;
    let equilibrium_vapor = 10.0 * volume / (r * bath_t);
    for pressure in [0.0, 20.0] {
        let vapor = pressure * volume / (r * bath_t);
        let dm = (equilibrium_vapor - vapor) * (1.0 - (-decay * dt).exp());
        let q = quadratic_heat(230.0, 220.0, dt, beta, 0.0, -dm * LATENT);
        let step = advance_resistive_surface_with_interface(
            &Phase,
            SurfaceKinetics {
                accommodation_coefficient: alpha,
                mass_transfer_coefficient_m_s: None,
            },
            initial,
            VaporReservoir::Finite {
                vapor_mass_kg: vapor,
                volume_m3: volume,
                temperature_k: bath_t,
            },
            dt,
            &interface,
            0.0,
            0.0,
            controls(),
            |h| Ok(h / CS),
            |h| Ok(h / CD),
            || Ok(()),
        )
        .unwrap();
        close(step.phase.transferred_mass_kg, dm, 5e-20);
        close(step.interface_heat_j, q, 1e-8);
        close(step.substrate_temperature_k, 230.0 - q / (MS * CS), 1e-8);
        close(
            step.phase.surface_temperature_k,
            220.0 + (q - dm * LATENT) / (MD * CD),
            1e-8,
        );
        close(
            step.substrate_energy_change_j
                + step.deposit_energy_change_j
                + step.phase.vapor_energy_j,
            0.0,
            5e-9,
        );
        close(step.energy_residual_j, 0.0, 5e-9);
        close(step.mass_residual_kg, 0.0, 1e-18);
    }
}

#[test]
fn constant_interface_limit_matches_the_original_fixed_conductance_owner() {
    let initial = state(300.0, 200.0);
    // Dyadic dt and G retain exactly the same products under either grouping.
    let (dt, g, qs, qd) = (0.125, 4.0, 7.0, -3.0);
    let interface = |ts: f64, td: f64| Ok(g * (ts - td));
    let nonlinear = heat(initial, dt, &interface, qs, qd);
    let fixed = heat_resistive_surface(
        &Phase,
        initial,
        dt,
        g,
        qs,
        qd,
        controls(),
        |h| Ok(h / CS),
        |h| Ok(h / CD),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(nonlinear, fixed);
    let reservoir = VaporReservoir::Prescribed {
        partial_pressure_pa: 20.0,
        temperature_k: 230.0,
    };
    let kinetics = SurfaceKinetics {
        accommodation_coefficient: 0.2,
        mass_transfer_coefficient_m_s: None,
    };
    let nonlinear = advance_resistive_surface_with_interface(
        &Phase,
        kinetics,
        initial,
        reservoir,
        dt,
        &interface,
        qs,
        qd,
        controls(),
        |h| Ok(h / CS),
        |h| Ok(h / CD),
        || Ok(()),
    )
    .unwrap();
    let fixed = advance_resistive_surface_with_heat(
        &Phase,
        kinetics,
        initial,
        reservoir,
        dt,
        g,
        qs,
        qd,
        controls(),
        |h| Ok(h / CS),
        |h| Ok(h / CD),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(nonlinear, fixed);
}

#[test]
fn nonpassive_nonfinite_and_equal_temperature_heat_are_physical_refusals() {
    for (initial, kind) in [
        (state(300.0, 200.0), 0),
        (state(300.0, 200.0), 1),
        (state(300.0, 200.0), 2),
        (state(300.0, 300.0), 3),
    ] {
        let calls = Cell::new(0);
        let broken = |ts: f64, td: f64| {
            calls.set(calls.get() + 1);
            Ok(match kind {
                0 => td - ts,
                1 => f64::NAN,
                2 => f64::INFINITY,
                _ => 1.0,
            })
        };
        let result = heat_resistive_surface_with_interface(
            &Phase,
            initial,
            0.5,
            &broken,
            0.0,
            0.0,
            controls(),
            |h| Ok(h / CS),
            |h| Ok(h / CD),
            || Ok(()),
        );
        assert!(
            matches!(result, Err(SublimationError::Invalid(_))),
            "{result:?}"
        );
        assert_eq!(
            calls.get(),
            1,
            "a physical violation must not become a bracket search"
        );
    }
}

#[test]
fn zero_time_skips_interface_and_nonlinear_trials_share_one_budget_and_cancellation() {
    let initial = state(300.0, 200.0);
    let calls = Cell::new(0);
    let interface = |ts: f64, td: f64| {
        calls.set(calls.get() + 1);
        Ok(0.002 * (ts - td) * (ts + td))
    };
    let zero = heat(initial, 0.0, &interface, 7.0, -3.0);
    assert_eq!(calls.get(), 0);
    assert_eq!(zero.interface_heat_j, 0.0);
    balances(zero, 7.0, -3.0);
    let mut limited = controls();
    limited.max_iterations = 1;
    let result = heat_resistive_surface_with_interface(
        &Phase,
        initial,
        0.5,
        &interface,
        0.0,
        0.0,
        limited,
        |h| Ok(h / CS),
        |h| Ok(h / CD),
        || Ok(()),
    );
    assert_eq!(result, Err(SublimationError::NonConverged));
    assert_eq!(calls.get(), 1);
    calls.set(0);
    let mut checkpoints = 0;
    let result = heat_resistive_surface_with_interface(
        &Phase,
        initial,
        0.5,
        &interface,
        0.0,
        0.0,
        controls(),
        |h| Ok(h / CS),
        |h| Ok(h / CD),
        || {
            checkpoints += 1;
            if checkpoints == 2 {
                Err(SublimationError::Cancelled)
            } else {
                Ok(())
            }
        },
    );
    assert_eq!(result, Err(SublimationError::Cancelled));
    assert_eq!(calls.get(), 0);
    let cancelled = |_: f64, _: f64| Err(SublimationError::Cancelled);
    assert_eq!(
        heat_resistive_surface_with_interface(
            &Phase,
            initial,
            0.5,
            &cancelled,
            0.0,
            0.0,
            controls(),
            |h| Ok(h / CS),
            |h| Ok(h / CD),
            || Ok(())
        ),
        Err(SublimationError::Cancelled)
    );
}
