//! Independent two-body, molecular-relaxation and conservation oracles.
use super::*;
use crate::sublimation::{MOLAR_GAS_CONSTANT, SublimationProperties};

const TS: f64 = 230.0;
const CP_S: f64 = 500.0;
const CP_D: f64 = 2000.0;
const LATENT: f64 = 2.8e6;
const MOLAR: f64 = 0.018;

struct Law;
impl SublimationLaw for Law {
    fn temperature_range_k(&self) -> [f64; 2] {
        [100.0, 400.0]
    }
    fn properties(&self, t: f64) -> Result<SublimationProperties, SublimationError> {
        // The outgoing molecular mass rate is independent of temperature.
        let p = 10.0 * (t / TS).sqrt();
        Ok(SublimationProperties {
            saturation_pressure_pa: p,
            saturation_pressure_derivative_pa_k: p / (2.0 * t),
            latent_heat_j_kg: LATENT,
            latent_heat_derivative_j_kg_k: 0.0,
            solid_density_kg_m3: 900.0,
            molar_mass_kg_mol: MOLAR,
        })
    }
}
fn state(td: f64) -> ResistiveSurfaceState {
    ResistiveSurfaceState {
        substrate_mass_kg: 0.02,
        substrate_specific_enthalpy_j_kg: CP_S * TS,
        deposit: SolidSurfaceState {
            mass_kg: 0.001,
            specific_enthalpy_j_kg: CP_D * td,
            area_m2: 0.001,
        },
    }
}
fn config() -> SolverConfig {
    SolverConfig {
        relative_tolerance: 1e-12,
        absolute_energy_tolerance_j: 1e-9,
        max_iterations: 256,
    }
}
fn kinetics(film: Option<f64>) -> SurfaceKinetics {
    SurfaceKinetics {
        accommodation_coefficient: 1.0,
        mass_transfer_coefficient_m_s: film,
    }
}
fn bath(p: f64) -> VaporReservoir {
    VaporReservoir::Prescribed {
        partial_pressure_pa: p,
        temperature_k: TS,
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}; tolerance {tolerance:e}"
    );
}
fn heat(state: ResistiveSurfaceState, dt: f64, g: f64, qs: f64, qd: f64) -> ResistiveHeatStep {
    heat_resistive_surface(
        &Law,
        state,
        dt,
        g,
        qs,
        qd,
        config(),
        |h| Ok(h / CP_S),
        |h| Ok(h / CP_D),
        || Ok(()),
    )
    .unwrap()
}
fn phase(
    state: ResistiveSurfaceState,
    reservoir: VaporReservoir,
    film: Option<f64>,
    dt: f64,
    g: f64,
    qs: f64,
    qd: f64,
) -> ResistiveSurfaceStep {
    advance_resistive_surface_with_heat(
        &Law,
        kinetics(film),
        state,
        reservoir,
        dt,
        g,
        qs,
        qd,
        config(),
        |h| Ok(h / CP_S),
        |h| Ok(h / CP_D),
        || Ok(()),
    )
    .unwrap()
}
// Independent solution of the 2x2 backward-Euler caloric/contact equations.
fn exact(initial: ResistiveSurfaceState, dt: f64, g: f64, qs: f64, qd: f64) -> [f64; 3] {
    let cs = initial.substrate_mass_kg * CP_S;
    let cd = initial.deposit.mass_kg * CP_D;
    let ts = initial.substrate_specific_enthalpy_j_kg / CP_S;
    let td = initial.deposit.specific_enthalpy_j_kg / CP_D;
    let contact = dt * g;
    let q = contact * (ts - td + qs / cs - qd / cd) / (1.0 + contact * (1.0 / cs + 1.0 / cd));
    [ts + (qs - q) / cs, td + (qd + q) / cd, q]
}
fn heat_balance(step: &ResistiveHeatStep, qs: f64, qd: f64) {
    close(
        step.substrate_energy_change_j + step.deposit_energy_change_j,
        qs + qd,
        3e-9,
    );
    close(step.substrate_energy_residual_j, 0.0, 3e-9);
    close(step.deposit_energy_residual_j, 0.0, 3e-9);
    close(step.energy_residual_j, 0.0, 3e-9);
}
fn phase_balance(step: &ResistiveSurfaceStep, qs: f64, qd: f64) {
    close(
        step.substrate_energy_change_j + step.deposit_energy_change_j + step.phase.vapor_energy_j,
        qs + qd,
        3e-9,
    );
    close(step.substrate_energy_residual_j, 0.0, 3e-9);
    close(step.phase.energy_residual_j, 0.0, 3e-9);
    close(step.energy_residual_j, 0.0, 3e-9);
    close(step.mass_residual_kg, 0.0, 1e-18);
}

#[test]
fn two_body_heat_matches_independent_backward_euler_and_flows_down_temperature() {
    let initial = state(260.0);
    let step = heat(initial, 3.0, 4.0, 7.0, -3.0);
    let [ts, td, q] = exact(initial, 3.0, 4.0, 7.0, -3.0);
    close(step.substrate_temperature_k, ts, 1e-8);
    close(step.deposit_temperature_k, td, 1e-8);
    close(step.interface_heat_j, q, 3e-9);
    assert!(
        step.interface_heat_j
            * (1.0 / step.deposit_temperature_k - 1.0 / step.substrate_temperature_k)
            >= 0.0
    );
    assert_eq!(step.state.substrate_mass_kg, initial.substrate_mass_kg);
    assert_eq!(step.state.deposit.mass_kg, initial.deposit.mass_kg);
    heat_balance(&step, 7.0, -3.0);
}

#[test]
fn refinement_converges_to_two_body_exponential_without_overshoot() {
    let initial = state(260.0);
    let rate: f64 = 4.0 * (1.0 / 10.0 + 1.0 / 2.0);
    let exact_gap = -30.0 * (-rate).exp();
    let mut errors = Vec::new();
    for n in [8, 16, 32] {
        let mut current = initial;
        for _ in 0..n {
            let step = heat(current, 1.0 / n as f64, 4.0, 0.0, 0.0);
            assert!(step.substrate_temperature_k <= step.deposit_temperature_k);
            assert!(step.substrate_temperature_k >= TS);
            assert!(step.deposit_temperature_k <= 260.0);
            heat_balance(&step, 0.0, 0.0);
            current = step.state;
        }
        let gap = current.substrate_specific_enthalpy_j_kg / CP_S
            - current.deposit.specific_enthalpy_j_kg / CP_D;
        errors.push((gap - exact_gap).abs());
    }
    assert!(
        errors[1] < 0.56 * errors[0] && errors[2] < 0.56 * errors[1],
        "{errors:?}"
    );
}

#[test]
fn insulated_and_zero_time_heat_keep_the_unheated_support_exactly_unchanged() {
    let initial = state(260.0);
    for (dt, g) in [(1.0, 0.0), (0.0, 1e6)] {
        let step = heat(initial, dt, g, 0.0, 9.0);
        assert_eq!(
            step.state.substrate_specific_enthalpy_j_kg,
            initial.substrate_specific_enthalpy_j_kg
        );
        assert_eq!(step.substrate_energy_change_j, 0.0);
        assert_eq!(step.interface_heat_j, 0.0);
        close(step.deposit_temperature_k, 264.5, 1e-12);
        assert_eq!(step.evaluations, 1);
        heat_balance(&step, 0.0, 9.0);
    }
    let step = phase(initial, bath(20.0), None, 0.1, 0.0, 0.0, 0.0);
    assert_eq!(
        step.state.substrate_specific_enthalpy_j_kg,
        initial.substrate_specific_enthalpy_j_kg
    );
    assert_eq!(step.substrate_energy_change_j, 0.0);
    assert_eq!(step.interface_heat_j, 0.0);
    phase_balance(&step, 0.0, 0.0);
}

#[test]
fn finite_conductance_recovers_common_temperature_limit_with_local_energy_admission() {
    let initial = state(TS);
    let common = TS + 10.0 / 12.0;
    let mut errors = Vec::new();
    for g in [10.0, 100.0, 1000.0, 10000.0] {
        let step = heat(initial, 0.1, g, 0.0, 10.0);
        let [ts, td, q] = exact(initial, 0.1, g, 0.0, 10.0);
        close(step.substrate_temperature_k, ts, 1e-8);
        close(step.deposit_temperature_k, td, 1e-8);
        close(step.interface_heat_j, q, 3e-9);
        heat_balance(&step, 0.0, 10.0);
        errors.push((step.deposit_temperature_k - common).abs());
    }
    assert!(errors.windows(2).all(|pair| pair[1] < 0.25 * pair[0]));
    assert!(errors[3] < 0.007);
    // An unresolved temperature jump must not silently become perfect contact
    // or inherit a large tolerance from the stiff predictor's fictitious heat.
    assert!(
        heat_resistive_surface(
            &Law,
            initial,
            0.1,
            1e18,
            0.0,
            10.0,
            config(),
            |h| Ok(h / CP_S),
            |h| Ok(h / CP_D),
            || Ok(())
        )
        .is_err()
    );
}

#[test]
fn thermal_resistance_heats_a_deposit_and_reduces_its_growth_through_vapor_pressure() {
    struct Clausius;
    impl SublimationLaw for Clausius {
        fn temperature_range_k(&self) -> [f64; 2] {
            [100.0, 400.0]
        }
        fn properties(&self, t: f64) -> Result<SublimationProperties, SublimationError> {
            let exponent = LATENT * MOLAR / MOLAR_GAS_CONSTANT;
            let pressure = 10.0 * fs_math::det::exp(exponent * (1.0 / TS - 1.0 / t));
            Ok(SublimationProperties {
                saturation_pressure_pa: pressure,
                saturation_pressure_derivative_pa_k: pressure * exponent / (t * t),
                ..Law.properties(t)?
            })
        }
    }
    let mut initial = state(TS);
    initial.deposit.mass_kg = 1e-4;
    let run = |conductance| {
        advance_resistive_surface_with_heat(
            &Clausius,
            kinetics(None),
            initial,
            bath(20.0),
            0.5,
            conductance,
            0.0,
            0.0,
            config(),
            |h| Ok(h / CP_S),
            |h| Ok(h / CP_D),
            || Ok(()),
        )
        .unwrap()
    };
    let insulated = run(0.02);
    let cooled = run(20.0);
    assert!(insulated.phase.surface_temperature_k > cooled.phase.surface_temperature_k);
    assert!(insulated.substrate_temperature_k < cooled.substrate_temperature_k);
    assert!(
        cooled.state.deposit.mass_kg - initial.deposit.mass_kg
            > 2.0 * (insulated.state.deposit.mass_kg - initial.deposit.mass_kg)
    );
    for step in [insulated, cooled] {
        assert!(step.phase.transferred_mass_kg < 0.0 && step.interface_heat_j < 0.0);
        let td = step.phase.surface_temperature_k;
        let pressure = 10.0 * (LATENT * MOLAR / MOLAR_GAS_CONSTANT * (1.0 / TS - 1.0 / td)).exp();
        let dm = 0.5
            * initial.deposit.area_m2
            * (MOLAR / (2.0 * core::f64::consts::PI * MOLAR_GAS_CONSTANT)).sqrt()
            * (pressure / td.sqrt() - 20.0 / TS.sqrt());
        close(step.phase.transferred_mass_kg, dm, 1e-17);
        phase_balance(&step, 0.0, 0.0);
    }
}

#[test]
fn deposition_and_sublimation_latent_heat_enter_the_deposit_before_the_support() {
    let initial = state(250.0);
    let coefficient = (MOLAR / (2.0 * core::f64::consts::PI * MOLAR_GAS_CONSTANT * TS)).sqrt();
    for pressure in [0.0, 20.0] {
        let dm = initial.deposit.area_m2 * 0.1 * coefficient * (10.0 - pressure);
        let step = phase(initial, bath(pressure), None, 0.1, 5.0, 3.0, -1.0);
        let [ts, td, q] = exact(initial, 0.1, 5.0, 3.0, -1.0 - dm * LATENT);
        close(step.phase.transferred_mass_kg, dm, 1e-18);
        close(
            step.state.deposit.mass_kg,
            initial.deposit.mass_kg - dm,
            1e-18,
        );
        close(step.substrate_temperature_k, ts, 1e-8);
        close(step.phase.surface_temperature_k, td, 1e-8);
        close(step.interface_heat_j, q, 3e-9);
        assert_eq!(step.state.substrate_mass_kg, initial.substrate_mass_kg);
        assert_ne!(
            step.substrate_temperature_k,
            step.phase.surface_temperature_k
        );
        phase_balance(&step, 3.0, -1.0);
    }
}

#[test]
fn finite_bath_and_gas_film_reuse_the_exact_molecular_relaxation() {
    let initial = state(250.0);
    let beta = MOLAR_GAS_CONSTANT * TS / (MOLAR * 0.001);
    let reservoir = VaporReservoir::Finite {
        vapor_mass_kg: 20.0 / beta,
        volume_m3: 0.001,
        temperature_k: TS,
    };
    let kinetic = (MOLAR / (2.0 * core::f64::consts::PI * MOLAR_GAS_CONSTANT * TS)).sqrt();
    let film = 0.4 * MOLAR / (MOLAR_GAS_CONSTANT * TS);
    let coefficient = kinetic / (1.0 + kinetic / film);
    let dm = (10.0 - 20.0) / beta * (1.0 - (-initial.deposit.area_m2 * coefficient * beta).exp());
    let step = phase(initial, reservoir, Some(0.4), 1.0, 2.0, 0.0, 0.0);
    let [ts, td, q] = exact(initial, 1.0, 2.0, 0.0, -dm * LATENT);
    close(step.phase.transferred_mass_kg, dm, 1e-20);
    close(step.substrate_temperature_k, ts, 1e-8);
    close(step.phase.surface_temperature_k, td, 1e-8);
    close(step.interface_heat_j, q, 3e-9);
    let VaporReservoir::Finite { vapor_mass_kg, .. } = step.phase.reservoir else {
        panic!()
    };
    close(vapor_mass_kg, 20.0 / beta + dm, 1e-20);
    assert!(vapor_mass_kg >= 10.0 / beta);
    phase_balance(&step, 0.0, 0.0);
}

#[test]
fn nonlinear_original_calorics_keep_independent_datums_and_exact_total_storage() {
    let mut temperatures = Vec::new();
    for (substrate_datum, deposit_datum) in [(0.0, 0.0), (1e6, -2e6)] {
        let initial = ResistiveSurfaceState {
            substrate_mass_kg: 0.02,
            substrate_specific_enthalpy_j_kg: substrate_datum + 5.0 * 50.0f64.powi(2),
            deposit: SolidSurfaceState {
                mass_kg: 0.001,
                specific_enthalpy_j_kg: deposit_datum + 5.0 * 200.0f64.powi(2),
                area_m2: 0.001,
            },
        };
        let step = heat_resistive_surface(
            &Law,
            initial,
            0.5,
            2.0,
            8.0,
            5.0,
            config(),
            |h| Ok(200.0 + ((h - substrate_datum) / 5.0).sqrt()),
            |h| Ok(100.0 + ((h - deposit_datum) / 5.0).sqrt()),
            || Ok(()),
        )
        .unwrap();
        heat_balance(&step, 8.0, 5.0);
        close(
            step.interface_heat_j,
            step.substrate_temperature_k - step.deposit_temperature_k,
            1e-12,
        );
        temperatures.push([step.substrate_temperature_k, step.deposit_temperature_k]);
    }
    for i in 0..2 {
        close(temperatures[0][i], temperatures[1][i], 1e-8);
    }
}

#[test]
fn equilibrium_is_exact_and_bare_or_extinct_interfaces_require_a_transition() {
    let initial = state(TS);
    let step = phase(initial, bath(10.0), None, 1.0, 10.0, 0.0, 0.0);
    assert_eq!(step.state, initial);
    assert_eq!(step.interface_heat_j, 0.0);
    assert_eq!(step.phase.transferred_mass_kg, 0.0);
    assert_eq!(step.energy_residual_j, 0.0);
    let mut bare = initial;
    bare.deposit.mass_kg = 0.0;
    assert!(matches!(
        heat_resistive_surface(
            &Law,
            bare,
            1.0,
            1.0,
            0.0,
            0.0,
            config(),
            |h| Ok(h / CP_S),
            |h| Ok(h / CP_D),
            || Ok(())
        ),
        Err(SublimationError::Invalid(_))
    ));
    let extinct = advance_resistive_surface_with_heat(
        &Law,
        kinetics(None),
        initial,
        bath(0.0),
        1000.0,
        0.0,
        0.0,
        initial.deposit.mass_kg * LATENT,
        config(),
        |h| Ok(h / CP_S),
        |h| Ok(h / CP_D),
        || Ok(()),
    );
    assert_eq!(
        extinct,
        Err(SublimationError::Invalid(
            "resistive coating extinction needs an interface transition"
        ))
    );
}

#[test]
fn original_domains_invalid_conductance_work_and_late_cancellation_refuse_atomically() {
    let initial = state(260.0);
    for g in [-1.0, f64::INFINITY, f64::NAN] {
        assert!(
            heat_resistive_surface(
                &Law,
                initial,
                1.0,
                g,
                0.0,
                0.0,
                config(),
                |h| Ok(h / CP_S),
                |h| Ok(h / CP_D),
                || Ok(())
            )
            .is_err()
        );
    }
    let limited = SolverConfig {
        max_iterations: 1,
        ..config()
    };
    assert_eq!(
        heat_resistive_surface(
            &Law,
            initial,
            1.0,
            2.0,
            0.0,
            0.0,
            limited,
            |h| Ok(h / CP_S),
            |h| Ok(h / CP_D),
            || Ok(())
        ),
        Err(SublimationError::NonConverged)
    );
    let mut polls = 0;
    heat_resistive_surface(
        &Law,
        initial,
        1.0,
        2.0,
        0.0,
        0.0,
        config(),
        |h| Ok(h / CP_S),
        |h| Ok(h / CP_D),
        || {
            polls += 1;
            Ok(())
        },
    )
    .unwrap();
    for stop in [1, polls] {
        let mut calls = 0;
        assert_eq!(
            heat_resistive_surface(
                &Law,
                initial,
                1.0,
                2.0,
                0.0,
                0.0,
                config(),
                |h| Ok(h / CP_S),
                |h| Ok(h / CP_D),
                || {
                    calls += 1;
                    if calls == stop {
                        Err(SublimationError::Cancelled)
                    } else {
                        Ok(())
                    }
                }
            ),
            Err(SublimationError::Cancelled)
        );
    }
    assert!(
        heat_resistive_surface(
            &Law,
            initial,
            1.0,
            2.0,
            0.0,
            0.0,
            config(),
            |_| Err(SublimationError::Invalid("unsupported substrate chart")),
            |h| Ok(h / CP_D),
            || Ok(())
        )
        .is_err()
    );
    assert_eq!(initial, state(260.0));
}
