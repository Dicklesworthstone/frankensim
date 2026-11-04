//! G0/G3: exact shared storage, domain backtracking and transactional refusal.
use super::*;
use core::cell::Cell;

fn heat(initial: SupportedSurfaceState, heat_j: f64) -> SupportedHeatStep {
    heat_supported_surface(
        &ManufacturedLaw,
        initial,
        heat_j,
        SolverConfig::default(),
        |h| Ok(h / 500.0),
        |t| Ok(2000.0 * t),
        || Ok(()),
    )
    .unwrap()
}

#[test]
fn heating_and_cooling_share_energy_between_both_original_caloric_charts() {
    let initial = SupportedSurfaceState {
        deposit_mass_kg: 0.005,
        ..state()
    };
    for heat_j in [-200.0, 200.0] {
        let step = heat(initial, heat_j);
        let capacity = initial.substrate_mass_kg * 500.0 + initial.deposit_mass_kg * 2000.0;
        near(
            step.surface_temperature_k,
            230.0 + heat_j / capacity,
            1e-9,
            0.0,
        );
        near(step.substrate_energy_change_j, heat_j / 2.0, 1e-8, 0.0);
        near(step.deposit_energy_change_j, heat_j / 2.0, 1e-8, 0.0);
        near(step.energy_residual_j, 0.0, 3e-8, 0.0);
        assert_eq!(step.state.substrate_mass_kg, initial.substrate_mass_kg);
        assert_eq!(step.state.deposit_mass_kg, initial.deposit_mass_kg);
        assert_eq!(step.state.area_m2, initial.area_m2);
    }
}

#[test]
fn nonlinear_deposit_storage_is_integrated_in_enthalpy_without_a_frozen_heat_capacity() {
    let initial = SupportedSurfaceState {
        deposit_mass_kg: 0.01,
        ..state()
    };
    for target in [150.0, 310.0] {
        let substrate_heat = initial.substrate_mass_kg * 500.0 * (target - 230.0);
        let deposit_heat = initial.deposit_mass_kg * 2.0 * (target * target - 230.0 * 230.0);
        let step = heat_supported_surface(
            &ManufacturedLaw,
            initial,
            substrate_heat + deposit_heat,
            SolverConfig::default(),
            |h| Ok(h / 500.0),
            |t| Ok(2.0 * t * t),
            || Ok(()),
        )
        .unwrap();
        near(step.surface_temperature_k, target, 1e-8, 0.0);
        near(step.substrate_energy_change_j, substrate_heat, 1e-7, 0.0);
        near(step.deposit_energy_change_j, deposit_heat, 1e-7, 0.0);
        near(step.energy_residual_j, 0.0, 2e-7, 0.0);
    }
}

#[test]
fn dominant_deposit_capacity_backtracks_unsupported_outer_temperature_trials() {
    let initial = SupportedSurfaceState {
        deposit_mass_kg: 0.1,
        ..state()
    };
    for heat_j in [-21000.0, 21000.0] {
        let unsupported = Cell::new(0);
        let step = heat_supported_surface(
            &ManufacturedLaw,
            initial,
            heat_j,
            SolverConfig::default(),
            |h| {
                let t = h / 500.0;
                if !(100.0..=400.0).contains(&t) {
                    unsupported.set(unsupported.get() + 1);
                }
                Ok(t)
            },
            |t| Ok(2000.0 * t),
            || Ok(()),
        )
        .unwrap();
        assert!(unsupported.get() > 0);
        near(
            step.surface_temperature_k,
            230.0 + heat_j / 210.0,
            2e-8,
            0.0,
        );
        near(
            step.substrate_energy_change_j + step.deposit_energy_change_j,
            heat_j,
            3e-6,
            0.0,
        );
    }
}

#[test]
fn bare_or_zero_area_support_stores_heat_and_zero_heat_preserves_state_exactly() {
    let bare = SupportedSurfaceState {
        area_m2: 0.0,
        ..state()
    };
    let step = heat(bare, 25.0);
    near(step.surface_temperature_k, 232.5, 1e-10, 0.0);
    near(step.substrate_energy_change_j, 25.0, 1e-10, 0.0);
    assert_eq!(step.deposit_energy_change_j, 0.0);
    assert_eq!(step.state.deposit_mass_kg, 0.0);
    assert_eq!(step.state.area_m2, 0.0);
    for initial in [
        bare,
        SupportedSurfaceState {
            deposit_mass_kg: 0.01,
            ..state()
        },
    ] {
        let zero = heat(initial, 0.0);
        assert_eq!(zero.state, initial);
        assert_eq!(zero.substrate_energy_change_j, 0.0);
        assert_eq!(zero.deposit_energy_change_j, 0.0);
        assert_eq!(zero.energy_residual_j, 0.0);
    }
}

#[test]
fn independent_large_caloric_datums_do_not_enter_heat_storage() {
    let initial = SupportedSurfaceState {
        deposit_mass_kg: 0.005,
        ..state()
    };
    let reference = heat(initial, 200.0);
    let substrate_shift = 1e12;
    let deposit_shift = -3e12;
    let shifted = heat_supported_surface(
        &ManufacturedLaw,
        SupportedSurfaceState {
            substrate_specific_enthalpy_j_kg: initial.substrate_specific_enthalpy_j_kg
                + substrate_shift,
            ..initial
        },
        200.0,
        SolverConfig::default(),
        |h| Ok((h - substrate_shift) / 500.0),
        |t| Ok(2000.0 * t + deposit_shift),
        || Ok(()),
    )
    .unwrap();
    near(
        shifted.surface_temperature_k,
        reference.surface_temperature_k,
        1e-10,
        0.0,
    );
    near(
        shifted.substrate_energy_change_j,
        reference.substrate_energy_change_j,
        1e-10,
        0.0,
    );
    near(
        shifted.deposit_energy_change_j,
        reference.deposit_energy_change_j,
        1e-10,
        0.0,
    );
    near(shifted.energy_residual_j, 0.0, 1e-10, 0.0);
}

#[test]
fn unrepresentable_heat_is_refused_or_explicitly_bounded_by_energy_tolerance() {
    let datum = 1e30;
    let initial = SupportedSurfaceState {
        substrate_specific_enthalpy_j_kg: datum,
        ..state()
    };
    let run = |heat_j| {
        heat_supported_surface(
            &ManufacturedLaw,
            initial,
            heat_j,
            SolverConfig::default(),
            |h| Ok(230.0 + (h - datum) / 500.0),
            |t| Ok(2000.0 * t),
            || Ok(()),
        )
    };
    assert_eq!(run(1e-3), Err(SublimationError::NonConverged));
    assert_eq!(run(-1e-3), Err(SublimationError::NonConverged));
    let below_tolerance = run(1e-12).unwrap();
    assert_eq!(below_tolerance.state, initial);
    assert_eq!(below_tolerance.energy_residual_j, -1e-12);
}

#[test]
fn invalid_states_failed_solves_and_cancellation_publish_no_candidate() {
    let initial = SupportedSurfaceState {
        deposit_mass_kg: 0.005,
        ..state()
    };
    let saved = initial;
    for invalid in [
        SupportedSurfaceState {
            substrate_mass_kg: 0.0,
            ..initial
        },
        SupportedSurfaceState {
            deposit_mass_kg: -1.0,
            ..initial
        },
        SupportedSurfaceState {
            area_m2: -1.0,
            ..initial
        },
    ] {
        assert!(
            heat_supported_surface(
                &ManufacturedLaw,
                invalid,
                0.0,
                SolverConfig::default(),
                |h| Ok(h / 500.0),
                |t| Ok(2000.0 * t),
                || Ok(()),
            )
            .is_err()
        );
    }
    for heat_j in [f64::NAN, f64::INFINITY, -f64::INFINITY] {
        assert!(
            heat_supported_surface(
                &ManufacturedLaw,
                initial,
                heat_j,
                SolverConfig::default(),
                |h| Ok(h / 500.0),
                |t| Ok(2000.0 * t),
                || Ok(()),
            )
            .is_err()
        );
    }
    assert_eq!(
        heat_supported_surface(
            &ManufacturedLaw,
            initial,
            200.0,
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
    assert!(
        heat_supported_surface(
            &ManufacturedLaw,
            initial,
            1e8,
            SolverConfig::default(),
            |h| Ok(h / 500.0),
            |t| Ok(2000.0 * t),
            || Ok(()),
        )
        .is_err()
    );
    let mut calls = 0;
    heat_supported_surface(
        &ManufacturedLaw,
        initial,
        200.0,
        SolverConfig::default(),
        |h| Ok(h / 500.0),
        |t| Ok(2000.0 * t),
        || {
            calls += 1;
            Ok(())
        },
    )
    .unwrap();
    for stop in [1, calls / 2, calls] {
        let mut count = 0;
        assert_eq!(
            heat_supported_surface(
                &ManufacturedLaw,
                initial,
                200.0,
                SolverConfig::default(),
                |h| Ok(h / 500.0),
                |t| Ok(2000.0 * t),
                || {
                    count += 1;
                    if count == stop {
                        Err(SublimationError::Cancelled)
                    } else {
                        Ok(())
                    }
                },
            ),
            Err(SublimationError::Cancelled)
        );
        assert_eq!(initial, saved);
    }
}
