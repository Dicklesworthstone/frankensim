//! Independent kinetic trajectory, finite heat and transactional window tests.
use super::*;
use fs_couple::cavitation::dissolved_gas::trajectory::{
    ThermalDissolutionTimeControls, ThermalDissolutionTrajectoryError as Error, advance_thermal_to,
};

fn time_controls(moles: f64) -> ThermalDissolutionTimeControls {
    ThermalDissolutionTimeControls {
        initial_step_s: 1000.0,
        minimum_step_s: 1e-5,
        maximum_step_s: 1000.0,
        maximum_evaluations: 10_000,
        relative_tolerance: 0.0,
        temperature_tolerance_k: 1e-4,
        pressure_tolerance_pa: 100.0,
        species_tolerance_mol: moles,
        energy_tolerance_j: 5e-5,
        mass_tolerance_kg: 1e-11,
    }
}

#[test]
fn g1_tighter_thermal_trajectories_resolve_warming_induced_exsolution() {
    // Independent one-coordinate ODE: total U eliminates T algebraically.
    // Gas diffusivity is proportional to T, so Dg/(R T delta_g) is constant;
    // the two-film resistance still changes through the original H(T).
    let total_gas = P0 * VG / (R * T0);
    let total_solute = 0.01 * total_gas;
    let (end, power) = (8000.0, 20.0);
    let temperature = |n: f64, t: f64| {
        T0 + (power * t + n * (SOLUTION_HEAT - R * T0))
            / (2000.0 + (total_gas - n) * 2.5 * R + n * M * CP_DISSOLVED)
    };
    let rate = |n: f64, t: f64| {
        let tk = temperature(n, t);
        let h = H * (1.0 - 0.004 * (tk - T0));
        let conductance = 0.01 / (1e-5 / D + h / (H * D / 1e-5));
        conductance * (h * (total_solute - n) * R * tk / VG - n / VL)
    };
    let mut oracle = 0.0;
    let dt = end / 4096.0;
    for i in 0..4096 {
        let t = i as f64 * dt;
        let k1 = rate(oracle, t);
        let k2 = rate(oracle + 0.5 * dt * k1, t + 0.5 * dt);
        let k3 = rate(oracle + 0.5 * dt * k2, t + 0.5 * dt);
        let k4 = rate(oracle + dt * k3, t + dt);
        oracle += dt / 6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }
    assert!(rate(0.0, 0.0) > 0.0 && rate(oracle, end) < 0.0);
    let mut errors = Vec::new();
    let mut work = Vec::new();
    for tolerance in [1e-9, 1e-10] {
        let mut gas = gas_with_fractions([0.01, 0.99]);
        let carrier = gas.species_mass_kg(&id("B")).unwrap();
        let (law, _) = sourced_gas_film(law(&gas, "A", SOLUTION_HEAT, true, 500.0), 500.0);
        let mut laws = [law];
        let mut liquids = [DissolvedGasInventory::new(VL, 0.0, 0.0).unwrap()];
        let film = DissolutionFilm {
            gas_side_transfer: None,
            ..film()
        };
        let report = advance_thermal_to(
            &mut laws,
            &mut gas,
            &mut liquids,
            &[film],
            &solvent(),
            end,
            power,
            controls(),
            time_controls(tolerance),
            || false,
        )
        .unwrap();
        let n = liquids[0].dissolved_moles();
        errors.push((n - oracle).abs());
        work.push(report.evaluations);
        assert!(report.rejected_intervals > 0 && report.maximum_accepted_error_ratio <= 1.0);
        assert!(report.species[0].absorbed_moles > 0.0 && report.species[0].exsolved_moles > 0.0);
        assert!(report.last_substep.endpoint_sources[0].molar_flux_mol_m2_s < 0.0);
        near(
            report.species[0].absorbed_moles - report.species[0].exsolved_moles,
            n,
            2e-15,
        );
        near(gas.state().state().temperature, temperature(n, end), 2e-6);
        near(gas.species_mass_kg(&id("B")).unwrap(), carrier, 1e-14);
        near(report.integrated_external_heat_j, power * end, 1e-8);
        near(
            report.energy_residual_j,
            0.0,
            time_controls(tolerance).energy_tolerance_j,
        );
        near(
            laws[0].gas_film_diffusivity_m2_s().unwrap(),
            2e-10 * gas.state().state().temperature / T0,
            1e-24,
        );
        assert_eq!(report.interval_s, [0.0, end]);
        assert_eq!(liquids[0].time_s(), end);
    }
    assert!(
        errors[1] < 0.5 * errors[0] && errors[1] < 7e-9,
        "{errors:?}"
    );
    assert!(work[1] > work[0], "{work:?}");
}

#[test]
fn g1_sealed_heat_uses_increment_retries_and_one_complete_heat_account() {
    let mut gas = gas();
    let original = gas.clone();
    let mut laws = [law(&gas, "A", SOLUTION_HEAT, true, 500.0)];
    let mut liquids = [DissolvedGasInventory::new(VL, 0.001, 0.0).unwrap()];
    let film = DissolutionFilm {
        area_m2: 0.0,
        ..film()
    };
    let mut step = controls();
    step.maximum_temperature_increment_k = 0.4;
    let report = advance_thermal_to(
        &mut laws,
        &mut gas,
        &mut liquids,
        &[film],
        &solvent(),
        100.0,
        100.0,
        step,
        time_controls(1e-9),
        || false,
    )
    .unwrap();
    let capacity = 2000.0 + P0 * VG / (R * T0) * 2.5 * R + 0.001 * M * CP_DISSOLVED;
    near(
        gas.state().state().temperature,
        T0 + 10_000.0 / capacity,
        1e-9,
    );
    assert!(report.rejected_intervals > 0);
    assert_eq!(gas.mass_kg(), original.mass_kg());
    assert_eq!(gas.model(), original.model());
    assert_eq!(liquids[0].dissolved_moles(), 0.001);
    assert_eq!(report.species[0].absorbed_moles, 0.0);
    assert_eq!(report.species[0].exsolved_moles, 0.0);
    assert_eq!(report.external_heat_j, 10_000.0);
    near(report.integrated_external_heat_j, 10_000.0, 1e-10);
    near(
        report.energy_residual_j,
        0.0,
        time_controls(1e-9).energy_tolerance_j,
    );
}

#[test]
fn g3_opposite_species_share_canonical_adaptive_work_and_gross_transfers() {
    let original = gas();
    let laws = [
        law(&original, "A", SOLUTION_HEAT, true, 500.0),
        law(&original, "B", 2.0 * SOLUTION_HEAT, true, 500.0),
    ];
    let liquids = [
        DissolvedGasInventory::new(VL, 0.0, 0.0).unwrap(),
        DissolvedGasInventory::new(VL, 0.003, 0.0).unwrap(),
    ];
    let mut a = (laws.clone(), original.clone(), liquids);
    let mut b = (
        [laws[1].clone(), laws[0].clone()],
        original,
        [liquids[1], liquids[0]],
    );
    let first = advance_thermal_to(
        &mut a.0,
        &mut a.1,
        &mut a.2,
        &[film(); 2],
        &solvent(),
        300.0,
        20.0,
        controls(),
        time_controls(1e-8),
        || false,
    )
    .unwrap();
    let second = advance_thermal_to(
        &mut b.0,
        &mut b.1,
        &mut b.2,
        &[film(); 2],
        &solvent(),
        300.0,
        20.0,
        controls(),
        time_controls(1e-8),
        || false,
    )
    .unwrap();
    assert_eq!(first, second);
    assert_eq!(a.1, b.1);
    assert_eq!(a.2, [b.2[1], b.2[0]]);
    assert!(first.species[0].absorbed_moles > 0.0 && first.species[1].exsolved_moles > 0.0);
    assert!(first.maximum_species_mass_residual_kg < time_controls(1e-8).mass_tolerance_kg);
}

#[test]
fn g3_actual_unequal_halves_at_large_clocks_control_the_fine_error() {
    let start = 2.0_f64.powi(52);
    let end = start + 3.0;
    let middle = start + 2.0; // Rounded midpoint: actual intervals are 2 s and 1 s.
    let original = gas();
    let laws = [law(&original, "A", SOLUTION_HEAT, true, 500.0)];
    let liquids = [DissolvedGasInventory::new(VL, 0.0, start).unwrap()];
    let mut coarse = (laws.clone(), original.clone(), liquids);
    let mut fine = coarse.clone();
    advance_thermal_mixture(
        &mut coarse.0,
        &mut coarse.1,
        &mut coarse.2,
        &[film()],
        &solvent(),
        end,
        6000.0,
        controls(),
        || false,
    )
    .unwrap();
    advance_thermal_mixture(
        &mut fine.0,
        &mut fine.1,
        &mut fine.2,
        &[film()],
        &solvent(),
        middle,
        4000.0,
        controls(),
        || false,
    )
    .unwrap();
    advance_thermal_mixture(
        &mut fine.0,
        &mut fine.1,
        &mut fine.2,
        &[film()],
        &solvent(),
        end,
        2000.0,
        controls(),
        || false,
    )
    .unwrap();
    let mut difference = (fine.2[0].dissolved_moles() - coarse.2[0].dissolved_moles()).abs();
    for name in ["A", "B"] {
        difference = difference.max(
            (fine.1.species_mass_kg(&id(name)).unwrap() / M
                - coarse.1.species_mass_kg(&id(name)).unwrap() / M)
                .abs(),
        );
    }
    assert!(difference > 0.0);
    let mut accuracy = time_controls(2.0 * difference);
    accuracy.initial_step_s = 3.0;
    accuracy.maximum_step_s = 3.0;
    accuracy.minimum_step_s = 1.0;
    accuracy.temperature_tolerance_k = 1e6;
    accuracy.pressure_tolerance_pa = 1e10;
    let mut result = (laws.clone(), original.clone(), liquids);
    let report = advance_thermal_to(
        &mut result.0,
        &mut result.1,
        &mut result.2,
        &[film()],
        &solvent(),
        end,
        2000.0,
        controls(),
        accuracy,
        || false,
    )
    .unwrap();
    assert_eq!(result, fine);
    near(report.maximum_accepted_error_ratio, 0.625, 1e-12); // (2²+1²)/(2*2*1) / 2
    assert_eq!(report.minimum_accepted_step_s, 1.0);
    assert_eq!(report.maximum_accepted_step_s, 2.0);
    assert_eq!(report.evaluations, 3);
    let mut refused = (laws.clone(), original.clone(), liquids);
    accuracy.initial_step_s = 1.0;
    accuracy.maximum_step_s = 1.0;
    assert!(matches!(
        advance_thermal_to(
            &mut refused.0,
            &mut refused.1,
            &mut refused.2,
            &[film()],
            &solvent(),
            start + 1.0,
            2000.0,
            controls(),
            accuracy,
            || false
        ),
        Err(Error::TimeResolution { .. })
    ));
    assert_eq!(refused, (laws, original, liquids));
}

#[test]
fn g4_window_limits_late_sources_and_final_cancellation_preserve_all_original_state() {
    let gas = gas();
    let laws = [law(&gas, "A", SOLUTION_HEAT, true, 500.0)];
    let liquids = [DissolvedGasInventory::new(VL, 0.0, 0.0).unwrap()];
    let before = (laws, gas, liquids);
    let mut accuracy = time_controls(1e-7);
    accuracy.initial_step_s = 60.0;
    accuracy.maximum_step_s = 60.0;
    let mut successful = before.clone();
    let mut checkpoints = 0;
    let expected = advance_thermal_to(
        &mut successful.0,
        &mut successful.1,
        &mut successful.2,
        &[film()],
        &solvent(),
        120.0,
        20.0,
        controls(),
        accuracy,
        || {
            checkpoints += 1;
            false
        },
    )
    .unwrap();
    let mut trial = before.clone();
    let mut calls = 0;
    assert!(matches!(
        advance_thermal_to(
            &mut trial.0,
            &mut trial.1,
            &mut trial.2,
            &[film()],
            &solvent(),
            120.0,
            20.0,
            controls(),
            accuracy,
            || {
                calls += 1;
                calls == checkpoints
            }
        ),
        Err(Error::Cancelled)
    ));
    assert_eq!(trial, before);
    assert_eq!(
        advance_thermal_to(
            &mut trial.0,
            &mut trial.1,
            &mut trial.2,
            &[film()],
            &solvent(),
            120.0,
            20.0,
            controls(),
            accuracy,
            || false
        )
        .unwrap(),
        expected
    );
    for minimum in [false, true] {
        let mut trial = before.clone();
        let mut limited = accuracy;
        if minimum {
            limited.minimum_step_s = 60.0;
            limited.species_tolerance_mol = 1e-18;
        } else {
            limited.maximum_evaluations = 3;
        }
        let error = advance_thermal_to(
            &mut trial.0,
            &mut trial.1,
            &mut trial.2,
            &[film()],
            &solvent(),
            120.0,
            20.0,
            controls(),
            limited,
            || false,
        )
        .unwrap_err();
        if minimum {
            assert!(matches!(error, Error::Accuracy { .. }));
        } else {
            assert!(matches!(
                error,
                Error::WorkLimit {
                    accepted_intervals: 1,
                    ..
                }
            ));
        }
        assert_eq!(trial, before);
    }
    let mut trial = before.clone();
    trial.0 = [law(&trial.1, "A", SOLUTION_HEAT, true, 310.0)];
    let bounded = trial.clone();
    accuracy.initial_step_s = 10.0;
    accuracy.maximum_step_s = 10.0;
    accuracy.species_tolerance_mol = 1e-5;
    accuracy.temperature_tolerance_k = 10.0;
    let error = advance_thermal_to(
        &mut trial.0,
        &mut trial.1,
        &mut trial.2,
        &[film()],
        &solvent(),
        100.0,
        1000.0,
        controls(),
        accuracy,
        || false,
    )
    .unwrap_err();
    assert!(
        matches!(error, Error::Step { time_s, error: ThermalDissolutionError::Dissolution(
        DissolvedGasError::Material(_)), .. } if time_s >= 10.0)
    );
    assert_eq!(trial, bounded);
}
