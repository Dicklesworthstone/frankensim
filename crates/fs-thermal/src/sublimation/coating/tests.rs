//! G0/G1/G3: bare growth, separate inventories, exact storage and datum shifts.
#[path = "heat_tests.rs"]
mod heat_tests;

use super::*;
use crate::sublimation::{MOLAR_GAS_CONSTANT, SublimationProperties};

// Manufactured saturation chart makes p_sat/sqrt(T) constant, so the accepted
// exchange is independent of endpoint T. The heat balance then has a separate
// closed-form solution without reproducing the nonlinear solver.
struct ManufacturedLaw;
impl SublimationLaw for ManufacturedLaw {
    fn temperature_range_k(&self) -> [f64; 2] {
        [100.0, 400.0]
    }
    fn properties(&self, t: f64) -> Result<SublimationProperties, SublimationError> {
        let p = 10.0 * (t / 230.0).sqrt();
        Ok(SublimationProperties {
            saturation_pressure_pa: p,
            saturation_pressure_derivative_pa_k: p / (2.0 * t),
            latent_heat_j_kg: 2.8e6,
            latent_heat_derivative_j_kg_k: 0.0,
            solid_density_kg_m3: 900.0,
            molar_mass_kg_mol: 0.018,
        })
    }
}
fn kinetics() -> SurfaceKinetics {
    SurfaceKinetics {
        accommodation_coefficient: 1.0,
        mass_transfer_coefficient_m_s: None,
    }
}
fn state() -> SupportedSurfaceState {
    SupportedSurfaceState {
        substrate_mass_kg: 0.02,
        substrate_specific_enthalpy_j_kg: 500.0 * 230.0,
        deposit_mass_kg: 0.0,
        area_m2: 0.001,
    }
}
fn bath(pressure: f64) -> VaporReservoir {
    VaporReservoir::Prescribed {
        partial_pressure_pa: pressure,
        temperature_k: 230.0,
    }
}
fn run(state: SupportedSurfaceState, reservoir: VaporReservoir, dt: f64) -> SupportedSurfaceStep {
    advance_supported_surface(
        &ManufacturedLaw,
        kinetics(),
        state,
        reservoir,
        dt,
        SolverConfig::default(),
        |h| Ok(h / 500.0),
        |t| Ok(2000.0 * t),
        || Ok(()),
    )
    .unwrap()
}
fn near(a: f64, b: f64, atol: f64, rtol: f64) {
    assert!(
        (a - b).abs() <= atol + rtol * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}
fn closed_balance(step: SupportedSurfaceStep) {
    near(
        step.substrate_energy_change_j + step.deposit_energy_change_j + step.phase.vapor_energy_j,
        0.0,
        2.0e-8,
        0.0,
    );
    near(step.phase.energy_residual_j, 0.0, 2.0e-8, 0.0);
}

#[test]
fn bare_substrate_grows_without_a_numerical_seed_and_releases_latent_heat() {
    let initial = state();
    let step = run(initial, bath(20.0), 1.0);
    let coefficient = (0.018 / (2.0 * core::f64::consts::PI * MOLAR_GAS_CONSTANT * 230.0)).sqrt();
    let dm = -coefficient * 10.0 * initial.area_m2;
    near(step.phase.transferred_mass_kg, dm, 1e-20, 1e-13);
    near(step.state.deposit_mass_kg, -dm, 1e-20, 1e-13);
    near(
        step.phase.surface_temperature_k,
        230.0 - dm * 2.8e6 / 10.0,
        1e-9,
        0.0,
    );
    assert_eq!(step.state.substrate_mass_kg, initial.substrate_mass_kg);
    assert!(step.substrate_energy_change_j > 0.0);
    assert!(step.phase.recession_m < 0.0);
    closed_balance(step);
}

#[test]
fn finite_vapor_reservoir_limits_the_same_growth_that_heats_the_substrate() {
    let volume = 0.001;
    let equilibrium = 10.0 * volume * 0.018 / (MOLAR_GAS_CONSTANT * 230.0);
    let reservoir = VaporReservoir::Finite {
        vapor_mass_kg: 2.0 * equilibrium,
        volume_m3: volume,
        temperature_k: 230.0,
    };
    let step = run(state(), reservoir, 1e4);
    near(step.state.deposit_mass_kg, equilibrium, 1e-20, 1e-13);
    let VaporReservoir::Finite { vapor_mass_kg, .. } = step.phase.reservoir else {
        panic!()
    };
    near(
        vapor_mass_kg + step.state.deposit_mass_kg,
        2.0 * equilibrium,
        1e-20,
        1e-13,
    );
    near(
        step.phase.surface_temperature_k,
        230.0 + equilibrium * 2.8e6 / 10.0,
        1e-9,
        0.0,
    );
    closed_balance(step);
}

#[test]
fn sublimation_can_exhaust_the_deposit_but_never_consumes_substrate() {
    let initial = SupportedSurfaceState {
        deposit_mass_kg: 1e-5,
        ..state()
    };
    let step = run(initial, bath(0.0), 1e4);
    assert_eq!(step.state.deposit_mass_kg, 0.0);
    assert_eq!(step.phase.transferred_mass_kg, initial.deposit_mass_kg);
    assert_eq!(step.state.substrate_mass_kg, initial.substrate_mass_kg);
    let capacity = 10.0 + initial.deposit_mass_kg * 2000.0;
    near(
        step.phase.surface_temperature_k,
        230.0 - 1e-5 * 2.8e6 / capacity,
        1e-9,
        0.0,
    );
    closed_balance(step);
    let empty = run(step.state, bath(0.0), 1e4);
    assert_eq!(empty.state, step.state);
    assert_eq!(empty.phase.transferred_mass_kg, 0.0);
    let regrown = run(empty.state, bath(20.0), 0.1);
    assert!(regrown.state.deposit_mass_kg > 0.0);
    closed_balance(regrown);
}

#[test]
fn sealed_zero_time_zero_area_and_undersaturated_bare_surfaces_are_noops() {
    for (initial, reservoir, dt, kinetics) in [
        (state(), bath(0.0), 1.0, kinetics()),
        (state(), bath(20.0), 0.0, kinetics()),
        (
            SupportedSurfaceState {
                area_m2: 0.0,
                ..state()
            },
            bath(20.0),
            1.0,
            kinetics(),
        ),
        (
            state(),
            bath(20.0),
            1.0,
            SurfaceKinetics {
                accommodation_coefficient: 0.0,
                ..kinetics()
            },
        ),
    ] {
        let step = advance_supported_surface(
            &ManufacturedLaw,
            kinetics,
            initial,
            reservoir,
            dt,
            SolverConfig::default(),
            |h| Ok(h / 500.0),
            |t| Ok(2000.0 * t),
            || Ok(()),
        )
        .unwrap();
        assert_eq!(step.state, initial);
        assert_eq!(step.phase.transferred_mass_kg, 0.0);
        closed_balance(step);
    }
}

#[test]
fn independent_substrate_and_phase_datums_preserve_the_physical_solution() {
    let original = run(state(), bath(20.0), 1.0);
    let substrate_shift = 1e6;
    let phase_shift = -3e6;
    let shifted = advance_supported_surface(
        &ManufacturedLaw,
        kinetics(),
        SupportedSurfaceState {
            substrate_specific_enthalpy_j_kg: state().substrate_specific_enthalpy_j_kg
                + substrate_shift,
            ..state()
        },
        bath(20.0),
        1.0,
        SolverConfig::default(),
        |h| Ok((h - substrate_shift) / 500.0),
        |t| Ok(2000.0 * t + phase_shift),
        || Ok(()),
    )
    .unwrap();
    near(
        shifted.phase.surface_temperature_k,
        original.phase.surface_temperature_k,
        1e-9,
        0.0,
    );
    near(
        shifted.phase.transferred_mass_kg,
        original.phase.transferred_mass_kg,
        1e-20,
        1e-13,
    );
    near(
        shifted.phase.vapor_energy_j - original.phase.vapor_energy_j,
        original.phase.transferred_mass_kg * phase_shift,
        1e-10,
        1e-13,
    );
    closed_balance(shifted);
}

#[test]
fn invalid_inventory_chart_cancellation_and_work_budget_refuse() {
    for initial in [
        SupportedSurfaceState {
            substrate_mass_kg: 0.0,
            ..state()
        },
        SupportedSurfaceState {
            deposit_mass_kg: -1.0,
            ..state()
        },
    ] {
        assert!(
            advance_supported_surface(
                &ManufacturedLaw,
                kinetics(),
                initial,
                bath(20.0),
                1.0,
                SolverConfig::default(),
                |h| Ok(h / 500.0),
                |t| Ok(2000.0 * t),
                || Ok(()),
            )
            .is_err()
        );
    }
    assert_eq!(
        advance_supported_surface(
            &ManufacturedLaw,
            kinetics(),
            state(),
            bath(20.0),
            1.0,
            SolverConfig::default(),
            |h| Ok(h / 500.0),
            |_| Err(SublimationError::Invalid("deposit chart not admitted")),
            || Ok(()),
        ),
        Err(SublimationError::Invalid("deposit chart not admitted"))
    );
    assert_eq!(
        advance_supported_surface(
            &ManufacturedLaw,
            kinetics(),
            state(),
            bath(20.0),
            1.0,
            SolverConfig::default(),
            |h| Ok(h / 500.0),
            |t| Ok(2000.0 * t),
            || Err(SublimationError::Cancelled),
        ),
        Err(SublimationError::Cancelled)
    );
    assert_eq!(
        advance_supported_surface(
            &ManufacturedLaw,
            kinetics(),
            state(),
            bath(20.0),
            1.0,
            SolverConfig {
                max_iterations: 1,
                ..SolverConfig::default()
            },
            |h| Ok(h / 500.0),
            |t| Ok(2000.0 * t),
            || Ok(()),
        ),
        Err(SublimationError::NonConverged)
    );
}
