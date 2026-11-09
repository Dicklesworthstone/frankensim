//! Finite-window accuracy and source-owned withdrawal/candidate admission.

use super::*;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::{
    SpeciesNetworkAdvanceControls, transport::GasTransportSolve,
};

fn temporal(tolerance: f64) -> SpeciesNetworkAdvanceControls {
    SpeciesNetworkAdvanceControls {
        maximum_step_s: 2.0,
        minimum_step_s: 1e-6,
        maximum_evaluations: 10_000,
        absolute_species_tolerance_kg: tolerance,
        relative_species_tolerance: 0.0,
        temperature_tolerance_k: 1.0,
    }
}

fn physical() -> GasTransportSolve {
    GasTransportSolve {
        exchange: controls(),
        entropy_tolerance_j_k: 1e-10,
    }
}

#[test]
fn refined_material_windows_approach_the_independent_species_relaxation() {
    let domain = VoxelDomain::graded([vec![0.05, 0.15], vec![0.1], vec![0.1]]).unwrap();
    let original = pair(&domain, 0.0);
    let accepted = carrier(&domain, &original);
    let field = flow(&domain, &[0.0001, 0.0002]);
    let duration = 2.0;
    let solve = |tolerance| {
        field
            .advance_turbulent_gas_species(
                &domain,
                &accepted,
                &original,
                0.7,
                duration,
                &endpoints(&original),
                physical(),
                temporal(tolerance),
                &CancelGate::new_clock_free(),
            )
            .unwrap()
    };
    let coarse = solve(1e-8);
    let fine = solve(1e-10);
    let mass = [original[&0].mass_kg(), original[&1].mass_kg()];
    let y = [fractions(&original[&0]), fractions(&original[&1])];
    let rho = [mass[0] / domain.volume(0), mass[1] / domain.volume(1)];
    let conductance = 0.01 / (0.025 / (rho[0] * 0.0001 / 0.7) + 0.075 / (rho[1] * 0.0002 / 0.7));
    let decay = (-conductance * (1.0 / mass[0] + 1.0 / mass[1]) * duration).exp();
    // Source temperature and partial enthalpy evolve, but these frozen common
    // mass-frame coefficients still have an exact independent linear species
    // solution. No numerical source state supplies this expected composition.
    let expected: [[f64; 2]; 2] = core::array::from_fn(|cell| {
        core::array::from_fn(|k| {
            let equilibrium = (mass[0] * y[0][k] + mass[1] * y[1][k]) / (mass[0] + mass[1]);
            mass[cell] * (equilibrium + (y[cell][k] - equilibrium) * decay)
        })
    });
    let error = |result: &fs_lbm::conjugate::TurbulentGasSpeciesAdvance| {
        (0..2)
            .flat_map(|cell| {
                let endpoint = &result.exchange.parcels()[&(cell as u64)];
                let actual = fractions(endpoint).map(|y| y * endpoint.mass_kg());
                (0..2).map(move |k| (actual[k] - expected[cell][k]).abs())
            })
            .fold(0.0_f64, f64::max)
    };
    assert!(error(&fine) < 0.3 * error(&coarse));
    assert!(error(&fine) < 1e-8);
    assert!(
        fine.exchange.report().accepted_intervals > coarse.exchange.report().accepted_intervals
    );
    for result in [&coarse, &fine] {
        let report = result.exchange.report();
        assert_eq!(report.duration_s, duration);
        assert_eq!(report.accepted_substeps, 2 * report.accepted_intervals);
        assert!(report.maximum_species_error_ratio <= 1.0);
        assert!(report.maximum_temperature_error_ratio <= 1.0);
        assert!(report.heat.is_empty());
        assert!(report.total_energy_residual_j.abs() <= controls().energy_tolerance_j);
        assert!(report.endpoint_energy_change_j.abs() <= controls().energy_tolerance_j);
        assert!(report.maximum_species_mass_residual_kg <= controls().mass_tolerance_kg);
        assert!(result.absolute_carrier_mass_defect_kg <= controls().mass_tolerance_kg);
    }
}

#[test]
fn typed_withdrawal_refusal_shrinks_the_window_and_keeps_refreshed_sources() {
    let domain = VoxelDomain::graded([vec![0.05, 0.15], vec![0.1], vec![0.1]]).unwrap();
    let original = pair(&domain, 1e18);
    let accepted = carrier(&domain, &original);
    let field = flow(&domain, &[0.0001, 0.0002]);
    let mut budget = physical();
    budget.exchange.maximum_species_fraction = 0.005;
    let gate = CancelGate::new_clock_free();
    assert!(
        field
            .diffuse_turbulent_gas_species(
                &domain,
                &accepted,
                &original,
                0.7,
                2.0,
                &endpoints(&original),
                budget.exchange,
                &gate,
            )
            .is_err()
    );
    let mut time = temporal(1.0);
    time.temperature_tolerance_k = 1000.0;
    let result = field
        .advance_turbulent_gas_species(
            &domain,
            &accepted,
            &original,
            0.7,
            2.0,
            &endpoints(&original),
            budget,
            time,
            &gate,
        )
        .unwrap();
    let report = result.exchange.report();
    assert_eq!(report.duration_s, 2.0);
    assert!(report.rejected_intervals > 0);
    assert!(report.accepted_intervals > 1);
    assert!(report.endpoint_energy_change_j.abs() <= budget.exchange.energy_tolerance_j);
    assert!(!report.transfers.is_empty());
    for (&id, after) in result.exchange.parcels() {
        assert_eq!(
            after.model().source_packs(),
            original[&id].model().source_packs()
        );
        assert_eq!(after.volume_m3(), original[&id].volume_m3());
        assert_ne!(
            after.state().state().temperature,
            original[&id].state().state().temperature
        );
    }
}

#[test]
fn work_exhaustion_after_an_accepted_prefix_publishes_no_partial_state() {
    let domain = VoxelDomain::graded([vec![0.05, 0.15], vec![0.1], vec![0.1]]).unwrap();
    let original = pair(&domain, 0.0);
    let accepted = carrier(&domain, &original);
    let field = flow(&domain, &[0.0001, 0.0002]);
    let saved = (original.clone(), accepted.clone(), field.clone());
    let mut time = temporal(1.0);
    time.maximum_step_s = 0.1;
    time.maximum_evaluations = 4;
    time.temperature_tolerance_k = 1000.0;
    let error = field
        .advance_turbulent_gas_species(
            &domain,
            &accepted,
            &original,
            0.7,
            1.0,
            &endpoints(&original),
            physical(),
            time,
            &CancelGate::new_clock_free(),
        )
        .unwrap_err();
    assert!(error.to_string().contains("WorkLimit"));
    assert_eq!(original, saved.0);
    assert_eq!(accepted, saved.1);
    assert_eq!(field, saved.2);
}
