use super::*;

fn parcel(mass: f64, initial: f64, bounds: [f64; 2]) -> CaloricInventory {
    CaloricInventory {
        mass_kg: mass,
        initial_specific_energy_j_kg: initial,
        specific_energy_bounds_j_kg: bounds,
    }
}
fn state(temperature: f64, slope: f64) -> CaloricState {
    CaloricState {
        temperature_k: temperature,
        temperature_slope_kg_k_j: slope,
    }
}
fn controls() -> CaloricEquilibriumControls {
    CaloricEquilibriumControls {
        energy_tolerance_j: 1e-8,
        ..CaloricEquilibriumControls::default()
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}
fn latent(energy: f64) -> CaloricState {
    state(
        if energy < 300.0 {
            energy
        } else if energy <= 500.0 {
            300.0
        } else {
            energy - 200.0
        },
        // Deliberately select the sensible one-sided derivative at BOTH corners.
        if energy > 300.0 && energy < 500.0 {
            0.0
        } else {
            1.0
        },
    )
}

#[test]
fn unequal_masses_and_original_capacities_conserve_energy_not_mean_temperature() {
    let inventory = [
        parcel(3.0, 1000.0, [200.0, 2000.0]),
        parcel(1.0, 1000.0, [500.0, 5000.0]),
    ];
    let result = equilibrate_caloric_inventory(
        &inventory,
        controls(),
        |i, energy| Ok::<_, &str>(state(energy / [2.0, 5.0][i], 1.0 / [2.0, 5.0][i])),
        || Ok(()),
    )
    .unwrap();
    let exact = 4000.0 / 11.0;
    for temperature in result.temperatures_k {
        close(temperature, exact, 2e-10);
    }
    close(result.energy_residual_j, 0.0, controls().energy_tolerance_j);
    assert!(!result.retained_initial_state);
}

#[test]
fn nonlinear_original_calorics_recover_the_energy_defined_temperature() {
    let inventory = [
        parcel(1.0, 40_000.0, [10_000.0, 1_000_000.0]),
        parcel(1.0, 160_000.0, [10_000.0, 1_000_000.0]),
    ];
    let result = equilibrate_caloric_inventory(
        &inventory,
        controls(),
        |_, energy| Ok::<_, &str>(state(energy.sqrt(), 0.5 / energy.sqrt())),
        || Ok(()),
    )
    .unwrap();
    for temperature in result.temperatures_k {
        close(temperature, 100_000.0_f64.sqrt(), 1e-10);
    }
    close(result.energy_residual_j, 0.0, controls().energy_tolerance_j);
    assert!(result.temperature_iterations > 0);
}

#[test]
fn large_independent_energy_datums_do_not_erase_small_conserved_exchange() {
    let datum = [1e16, -1e16];
    let inventory = [
        parcel(1.0, datum[0] + 600.0, [datum[0] + 100.0, datum[0] + 1000.0]),
        parcel(3.0, datum[1] + 200.0, [datum[1] + 100.0, datum[1] + 1000.0]),
    ];
    let result = equilibrate_caloric_inventory(
        &inventory,
        controls(),
        |i, energy| Ok::<_, &str>(state(energy - datum[i], 1.0)),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(result.temperatures_k, [300.0, 300.0]);
    assert_eq!(result.energy_residual_j, 0.0);
    assert_eq!(
        result.specific_energies_j_kg,
        [datum[0] + 300.0, datum[1] + 300.0]
    );
}

#[test]
fn one_latent_parcel_gets_only_the_energy_left_by_the_original_sensible_inverse() {
    let inventory = [
        parcel(2.0, 525.0, [100.0, 1000.0]),
        parcel(1.0, 200.0, [100.0, 1000.0]),
    ];
    let result = equilibrate_caloric_inventory(
        &inventory,
        controls(),
        |i, energy| {
            Ok::<_, &str>(if i == 0 {
                latent(energy)
            } else {
                state(energy, 1.0)
            })
        },
        || Ok(()),
    )
    .unwrap();
    assert_eq!(result.temperatures_k, [300.0, 300.0]);
    close(result.specific_energies_j_kg[0], 475.0, 1e-12);
    close(result.specific_energies_j_kg[1], 300.0, 1e-12);
    assert_eq!(result.energy_residual_j, 0.0);
}

#[test]
fn several_latent_parcels_refuse_an_unselected_phase_allocation() {
    let inventory = [
        parcel(1.0, 100.0, [100.0, 1000.0]),
        parcel(1.0, 700.0, [100.0, 1000.0]),
    ];
    assert!(matches!(
        equilibrate_caloric_inventory(
            &inventory,
            controls(),
            |_, energy| Ok::<_, &str>(latent(energy)),
            || Ok(())
        ),
        Err(CaloricEquilibriumError::UnderdeterminedPhaseAllocation { cells: 2 })
    ));
}

#[test]
fn loose_energy_tolerance_cannot_choose_an_underdetermined_latent_endpoint() {
    let inventory = [
        parcel(1.0, 100.0, [100.0, 1000.0]),
        parcel(1.0, 500.125, [100.0, 1000.0]),
    ];
    let mut loose = controls();
    loose.energy_tolerance_j = 1.0;
    assert!(matches!(
        equilibrate_caloric_inventory(
            &inventory,
            loose,
            |_, energy| Ok::<_, &str>(latent(energy)),
            || Ok(())
        ),
        Err(CaloricEquilibriumError::UnderdeterminedPhaseAllocation { cells: 2 })
    ));
}

#[test]
fn already_equilibrated_latent_allocations_are_retained_exactly() {
    let inventory = [
        parcel(1.0, 350.125, [100.0, 1000.0]),
        parcel(3.0, 499.75, [100.0, 1000.0]),
    ];
    let result = equilibrate_caloric_inventory(
        &inventory,
        controls(),
        |_, energy| Ok::<_, &str>(latent(energy)),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(result.specific_energies_j_kg, [350.125, 499.75]);
    assert_eq!(result.temperature_iterations, 0);
    assert_eq!(result.caloric_evaluations, 4);
    assert_eq!(result.energy_residual_j, 0.0);
    assert!(result.retained_initial_state);
}

#[test]
fn all_solid_or_all_liquid_latent_endpoints_have_unique_energy_allocation() {
    for (initial, expected) in [([100.0, 500.0], 300.0), ([600.0, 400.0], 500.0)] {
        let inventory = initial.map(|energy| parcel(1.0, energy, [100.0, 1000.0]));
        let result = equilibrate_caloric_inventory(
            &inventory,
            controls(),
            |_, energy| Ok::<_, &str>(latent(energy)),
            || Ok(()),
        )
        .unwrap();
        close(result.specific_energies_j_kg[0], expected, 1e-12);
        close(result.specific_energies_j_kg[1], expected, 1e-12);
        assert_eq!(result.energy_residual_j, 0.0);
    }
}

#[test]
fn incompatible_source_domains_and_energy_ranges_refuse_without_extrapolation() {
    let inventory = [
        parcel(1.0, 150.0, [100.0, 200.0]),
        parcel(1.0, 350.0, [300.0, 400.0]),
    ];
    assert!(matches!(
        equilibrate_caloric_inventory(
            &inventory,
            controls(),
            |_, energy| Ok::<_, &str>(state(energy, 1.0)),
            || Ok(())
        ),
        Err(CaloricEquilibriumError::NoCommonTemperature)
    ));
    let inventory = [
        parcel(100.0, 100.0, [100.0, 200.0]),
        parcel(1.0, 200.0, [150.0, 250.0]),
    ];
    assert!(matches!(
        equilibrate_caloric_inventory(
            &inventory,
            controls(),
            |_, energy| Ok::<_, &str>(state(energy, 1.0)),
            || Ok(())
        ),
        Err(CaloricEquilibriumError::NoConservedEnergyRoot)
    ));
}

#[test]
fn cancellation_source_refusal_and_exact_work_budgets_do_not_change_inventory() {
    let inventory = [
        parcel(1.0, 200.0, [100.0, 1000.0]),
        parcel(1.0, 600.0, [100.0, 1000.0]),
    ];
    let before = inventory;
    let mut limited = controls();
    limited.maximum_caloric_evaluations = 1;
    assert!(matches!(
        equilibrate_caloric_inventory(
            &inventory,
            limited,
            |_, energy| Ok::<_, &str>(state(energy, 1.0)),
            || Ok(())
        ),
        Err(CaloricEquilibriumError::EvaluationBudget {
            used: 1,
            maximum: 1
        })
    ));
    assert!(matches!(
        equilibrate_caloric_inventory(
            &inventory,
            controls(),
            |_, energy| Ok::<_, &str>(state(energy, 1.0)),
            || Err("cancelled")
        ),
        Err(CaloricEquilibriumError::Source("cancelled"))
    ));
    assert!(matches!(
        equilibrate_caloric_inventory(
            &inventory,
            controls(),
            |_, _| Err::<CaloricState, _>("source-domain"),
            || Ok(())
        ),
        Err(CaloricEquilibriumError::Source("source-domain"))
    ));
    assert_eq!(inventory, before);
}

#[test]
fn unrepresentable_conserved_exchange_is_not_rounded_into_equilibrium() {
    let inventory = [
        parcel(f64::MIN_POSITIVE, 1e-100, [0.0, 1e-99]),
        parcel(f64::MIN_POSITIVE, 3e-100, [0.0, 1e-99]),
    ];
    assert!(matches!(
        equilibrate_caloric_inventory(
            &inventory,
            controls(),
            |_, energy| Ok::<_, &str>(state(100.0 + energy * 1e100, 1e100)),
            || Ok(())
        ),
        Err(CaloricEquilibriumError::Unrepresentable(_))
    ));
}
