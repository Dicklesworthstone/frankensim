//! Exact applied power on the same source-qualified Cp/Cv cell.
//! All references below integrate the declared affine capacity independently.

use super::*;
use fs_conduction::lumped::exchange::caloric::solve_lumped_caloric_driven_exchange as solve_driven;

fn affine_impulse_mean(initial_k: f64, final_k: f64) -> f64 {
    // Integral T*C(T)dT / integral C(T)dT, for C(T)=200+0.5T.
    // Multiplying the entire capacity by two for Cp cancels from this ratio.
    let average = 0.5 * (initial_k + final_k);
    (200.0 * average
        + 0.5 * (initial_k * initial_k + initial_k * final_k + final_k * final_k) / 3.0)
        / (200.0 + 0.5 * average)
}

#[test]
fn g1_constant_power_crosses_cp_cv_knots_and_allows_negative_virtual_equilibrium() {
    let (cv_source, cp_source) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (cv_source.cv(), cp_source.cp());
    let cv_conduction =
        IsochoricConductionLaw::try_new(&cv_source.card, &cv, &cv_source.conductivity, SELECTION)
            .unwrap();
    let cp_conduction =
        IsobaricConductionLaw::try_new(&cp_source.card, &cp, &cp_source.conductivity, SELECTION)
            .unwrap();
    let cv_layer = CaloricConductionLayer::from_isochoric(&cv_conduction).unwrap();
    let cp_layer = CaloricConductionLayer::from_isobaric(&cp_conduction).unwrap();
    with_cx(|cx| {
        for (layer, isochoric) in [(&cv_layer, true), (&cp_layer, false)] {
            for (initial_k, final_k, bath_k, power_w) in [
                (360.0, 600.0, 320.0, 800.0),
                (650.0, 350.0, 400.0, -1200.0),
            ] {
                let mut body = cell(energy(isochoric, initial_k), bath_k);
                body.boundary.source_w = power_w;
                // This absolute virtual equilibrium appears ONLY in the
                // independent analytic oracle. The physical path is positive
                // and supported even when the virtual value is -200 K.
                let virtual_k = bath_k + power_w / 2.0;
                let duration =
                    body.mass_kg * exposure(isochoric, initial_k, final_k, virtual_k) / 2.0;
                let expected_heat =
                    body.mass_kg * (energy(isochoric, final_k) - energy(isochoric, initial_k));
                let expected_mean = virtual_k - expected_heat / (2.0 * duration);
                let result = solve_driven(cx, &body, layer, GATE, config(duration)).unwrap();
                assert_eq!(
                    result.caloric.storage_basis,
                    if isochoric {
                        ThermalStorageBasis::InternalEnergy
                    } else {
                        ThermalStorageBasis::Enthalpy
                    }
                );
                assert_eq!(
                    result.caloric_identity,
                    layer.caloric_law().caloric_identity().unwrap()
                );
                near(result.caloric.temperature_k, final_k, 2e-9);
                near(result.caloric.mean_temperature_k, expected_mean, 2e-9);
                near(result.heat_into_body_j, expected_heat, 2e-7);
                near(result.stored_energy_change_j, expected_heat, 2e-7);
                near(result.source_heat_j, power_w * duration, 2e-7);
                near(
                    result.ambient_heat_j,
                    2.0 * duration * (bath_k - expected_mean),
                    2e-7,
                );
                assert_eq!(result.reservoir_heat_j, 0.0);
                near(
                    result.source_heat_j + result.ambient_heat_j,
                    result.heat_into_body_j,
                    config(duration).energy_tolerance_j,
                );
                assert!(result.heat_residual_j.abs() <= config(duration).energy_tolerance_j);
                assert!(result.heat_partition_residual_j.abs() <= config(duration).energy_tolerance_j);
                assert!(result.energy_balance_residual_j.abs() <= config(duration).energy_tolerance_j);

                let first = solve_driven(cx, &body, layer, GATE, config(0.5 * duration)).unwrap();
                let mut midpoint = body.clone();
                midpoint.specific_h_j_kg = first.caloric.final_specific_energy_j_kg;
                let second =
                    solve_driven(cx, &midpoint, layer, GATE, config(0.5 * duration)).unwrap();
                near(second.caloric.temperature_k, result.caloric.temperature_k, 2e-9);
                near(
                    first.heat_into_body_j + second.heat_into_body_j,
                    result.heat_into_body_j,
                    2e-7,
                );
            }
        }
    });
}

#[test]
fn g1_isolated_source_adds_signed_native_energy_and_exact_capacity_weighted_mean() {
    let (cv_source, cp_source) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (cv_source.cv(), cp_source.cp());
    let cv_conduction =
        IsochoricConductionLaw::try_new(&cv_source.card, &cv, &cv_source.conductivity, SELECTION)
            .unwrap();
    let cp_conduction =
        IsobaricConductionLaw::try_new(&cp_source.card, &cp, &cp_source.conductivity, SELECTION)
            .unwrap();
    let cv_layer = CaloricConductionLayer::from_isochoric(&cv_conduction).unwrap();
    let cp_layer = CaloricConductionLayer::from_isobaric(&cp_conduction).unwrap();
    with_cx(|cx| {
        for (layer, isochoric) in [(&cv_layer, true), (&cp_layer, false)] {
            for (initial_k, final_k, power_w) in [(350.0, 650.0, 600.0), (650.0, 350.0, -600.0)] {
                let mut body = cell(energy(isochoric, initial_k), 500.0);
                body.boundary.htc_w_m2_k = 0.0;
                body.boundary.source_w = power_w;
                let expected_heat =
                    body.mass_kg * (energy(isochoric, final_k) - energy(isochoric, initial_k));
                let duration = expected_heat / power_w;
                let result = solve_driven(cx, &body, layer, GATE, config(duration)).unwrap();
                assert_eq!(result.total_conductance_w_k, 0.0);
                assert_eq!(result.maximum_biot, 0.0);
                near(result.caloric.temperature_k, final_k, 2e-9);
                near(
                    result.caloric.mean_temperature_k,
                    affine_impulse_mean(initial_k, final_k),
                    2e-9,
                );
                near(
                    result.caloric.final_specific_energy_j_kg,
                    body.specific_h_j_kg + power_w * duration / body.mass_kg,
                    3e-7,
                );
                near(result.heat_into_body_j, expected_heat, 2e-7);
                near(result.source_heat_j, expected_heat, 2e-7);
                assert_eq!(result.ambient_heat_j, 0.0);
                assert_eq!(result.reservoir_heat_j, 0.0);
                assert_eq!(result.heat_initial_energy_derivative_kg, 0.0);
                assert_eq!(result.heat_ambient_temperature_derivative_j_k, 0.0);
                assert_eq!(result.heat_reservoir_temperature_derivative_j_k, 0.0);
                near(result.heat_source_power_derivative_s, duration, 1e-12);
                near(result.heat_duration_derivative_w, power_w, 1e-12);
                assert!(result.energy_balance_residual_j.abs() <= config(duration).energy_tolerance_j);
            }
        }
    });
}

#[test]
fn g1_applied_source_and_two_baths_keep_separate_heat_accounts_at_rest_and_in_motion() {
    let source = Fixture::new(true);
    let cv = source.cv();
    let conduction =
        IsochoricConductionLaw::try_new(&source.card, &cv, &source.conductivity, SELECTION).unwrap();
    let layer = CaloricConductionLayer::from_isochoric(&conduction).unwrap();
    let mut body = cell(energy(true, 360.0), 320.0);
    body.boundary.reservoir = Some(Reservoir {
        area_m2: 0.0625,
        conductance_w_k: 0.5,
        temperature_k: 680.0,
    });
    body.boundary.source_w = 300.0;
    // The actual bath reduction is 410 K; source power moves the virtual
    // equilibrium to 560 K. No source joules belong to either reservoir.
    let duration = body.mass_kg * exposure(true, 360.0, 500.0, 560.0) / 2.0;
    let heat = body.mass_kg * (energy(true, 500.0) - energy(true, 360.0));
    let mean = 560.0 - heat / (2.0 * duration);
    with_cx(|cx| {
        let result = solve_driven(cx, &body, &layer, GATE, config(duration)).unwrap();
        near(result.caloric.temperature_k, 500.0, 2e-9);
        near(result.heat_into_body_j, heat, 1e-7);
        near(result.source_heat_j, 300.0 * duration, 1e-7);
        near(result.ambient_heat_j, 1.5 * duration * (320.0 - mean), 1e-7);
        near(result.reservoir_heat_j, 0.5 * duration * (680.0 - mean), 1e-7);
        assert!(result.source_heat_j > 0.0);
        assert!(result.ambient_heat_j < 0.0);
        assert!(result.reservoir_heat_j > 0.0);
        near(
            result.source_heat_j + result.ambient_heat_j + result.reservoir_heat_j,
            result.stored_energy_change_j,
            config(duration).energy_tolerance_j,
        );
        assert!(result.heat_partition_residual_j.abs() <= config(duration).energy_tolerance_j);
        assert!(result.energy_balance_residual_j.abs() <= config(duration).energy_tolerance_j);

        body.specific_h_j_kg = energy(true, 560.0);
        let steady = solve_driven(cx, &body, &layer, GATE, config(50.0)).unwrap();
        assert_eq!(steady.caloric.temperature_k, 560.0);
        assert_eq!(steady.caloric.final_specific_energy_j_kg, body.specific_h_j_kg);
        assert_eq!(steady.heat_into_body_j, 0.0);
        near(steady.source_heat_j, 15000.0, 1e-9);
        near(steady.ambient_heat_j, -18000.0, 1e-9);
        near(steady.reservoir_heat_j, 3000.0, 1e-9);
        near(steady.heat_duration_derivative_w, 0.0, 1e-12);
        near(steady.energy_balance_residual_j, 0.0, 1e-9);
    });
}

#[test]
fn g1_source_power_and_boundary_jacobian_follow_the_same_physical_forcing() {
    let source = Fixture::new(true);
    let cv = source.cv();
    let conduction =
        IsochoricConductionLaw::try_new(&source.card, &cv, &source.conductivity, SELECTION).unwrap();
    let layer = CaloricConductionLayer::from_isochoric(&conduction).unwrap();
    let mut body = cell(energy(true, 360.0), 320.0);
    body.boundary.reservoir = Some(Reservoir {
        area_m2: 0.0625,
        conductance_w_k: 0.5,
        temperature_k: 680.0,
    });
    body.boundary.source_w = 300.0;
    with_cx(|cx| {
        let initial = body.specific_h_j_kg;
        let duration = 100.0;
        let result = solve_driven(cx, &body, &layer, GATE, config(duration)).unwrap();
        let heat = |e: f64, ambient: f64, reservoir: f64, power: f64, time: f64| {
            let mut perturbed = body.clone();
            perturbed.specific_h_j_kg = e;
            perturbed.boundary.ambient_k = ambient;
            perturbed.boundary.reservoir.as_mut().unwrap().temperature_k = reservoir;
            perturbed.boundary.source_w = power;
            solve_driven(cx, &perturbed, &layer, GATE, config(time))
                .unwrap()
                .heat_into_body_j
        };
        let differences = [
            heat(initial + 0.5, 320.0, 680.0, 300.0, duration)
                - heat(initial - 0.5, 320.0, 680.0, 300.0, duration),
            (heat(initial, 320.01, 680.0, 300.0, duration)
                - heat(initial, 319.99, 680.0, 300.0, duration))
                / 0.02,
            (heat(initial, 320.0, 680.01, 300.0, duration)
                - heat(initial, 320.0, 679.99, 300.0, duration))
                / 0.02,
            (heat(initial, 320.0, 680.0, 300.01, duration)
                - heat(initial, 320.0, 680.0, 299.99, duration))
                / 0.02,
            (heat(initial, 320.0, 680.0, 300.0, duration + 0.01)
                - heat(initial, 320.0, 680.0, 300.0, duration - 0.01))
                / 0.02,
        ];
        let derivatives = [
            result.heat_initial_energy_derivative_kg,
            result.heat_ambient_temperature_derivative_j_k,
            result.heat_reservoir_temperature_derivative_j_k,
            result.heat_source_power_derivative_s,
            result.heat_duration_derivative_w,
        ];
        for (analytic, measured) in derivatives.into_iter().zip(differences) {
            near(analytic, measured, 5e-7 * measured.abs().max(1.0));
        }
        assert!(result.heat_source_power_derivative_s > 0.0);
        assert!(result.heat_source_power_derivative_s < duration);
        // A bath perturbation supplies G*dT watts of the same original forcing.
        near(
            result.heat_ambient_temperature_derivative_j_k,
            1.5 * result.heat_source_power_derivative_s,
            1e-9,
        );
        near(
            result.heat_reservoir_temperature_derivative_j_k,
            0.5 * result.heat_source_power_derivative_s,
            1e-9,
        );
        let final_k = result.caloric.temperature_k;
        near(
            result.heat_duration_derivative_w,
            300.0 + 1.5 * (320.0 - final_k) + 0.5 * (680.0 - final_k),
            1e-8,
        );
    });
}

#[test]
fn g0_driven_paths_must_stay_in_source_support_and_meet_the_native_energy_budget() {
    let source = Fixture::new(true);
    let cv = source.cv();
    let conduction =
        IsochoricConductionLaw::try_new(&source.card, &cv, &source.conductivity, SELECTION).unwrap();
    let layer = CaloricConductionLayer::from_isochoric(&conduction).unwrap();
    with_cx(|cx| {
        for (initial_k, ambient_k, htc, power, duration) in [
            (350.0, 500.0, 0.0, -1000.0, 10.0),
            (650.0, 500.0, 0.0, 1000.0, 20.0),
            (650.0, 320.0, 8.0, 800.0, 1000.0),
            (350.0, 400.0, 8.0, -1200.0, 1000.0),
        ] {
            let mut body = cell(energy(true, initial_k), ambient_k);
            body.boundary.htc_w_m2_k = htc;
            body.boundary.source_w = power;
            assert!(solve_driven(cx, &body, &layer, GATE, config(duration)).is_err());
            let stopped = solve_driven(cx, &body, &layer, GATE, config(0.0)).unwrap();
            assert_eq!(stopped.caloric.final_specific_energy_j_kg, body.specific_h_j_kg);
            assert_eq!(stopped.heat_into_body_j, 0.0);
            assert_eq!(stopped.source_heat_j, 0.0);
        }

        // At this datum a 1e-12 J impulse is below the cell's native storage
        // resolution. A tighter budget must reject the lost stored increment.
        let mut weak = cell(energy(true, 500.0), 500.0);
        weak.boundary.htc_w_m2_k = 0.0;
        weak.boundary.source_w = 1e-12;
        let mut strict = config(1.0);
        strict.energy_tolerance_j = 1e-15;
        assert!(solve_driven(cx, &weak, &layer, GATE, strict).is_err());
    });
}

#[test]
fn g1_weak_applied_power_keeps_heat_and_bath_loss_when_absolute_temperature_cannot_move() {
    let source = Fixture::new(true);
    let cv = IsochoricHeatCapacityLaw::try_new(
        &source.card,
        IsochoricHeatCapacityQuery {
            specific_heat_cv: &source.capacity,
            temperature_axis: "T",
            density_axis: "rho",
            lower: &source.lower,
            upper: &source.upper,
            selection: SELECTION,
        },
        500.0,
    )
    .unwrap();
    let conduction =
        IsochoricConductionLaw::try_new(&source.card, &cv, &source.conductivity, SELECTION).unwrap();
    let layer = CaloricConductionLayer::from_isochoric(&conduction).unwrap();
    let duration = 1.0_f64;
    let attenuation = 2.0 * duration / (0.4 * 450.0);
    let retained_fraction = -(-attenuation).exp_m1() / attenuation;
    // Independent expansion of 1-(1-exp(-x))/x. The first omitted term
    // is below one part in 1e16 of the bath heat at x=1/90.
    let lost_fraction = attenuation
        * (0.5 + attenuation
            * (-1.0 / 6.0 + attenuation
                * (1.0 / 24.0 + attenuation
                    * (-1.0 / 120.0 + attenuation
                        * (1.0 / 720.0 - attenuation / 5040.0)))));
    with_cx(|cx| {
        for power in [-1e-12, 1e-12] {
            let mut body = cell(0.0, 500.0);
            body.boundary.source_w = power;
            let result = solve_driven(cx, &body, &layer, GATE, config(duration)).unwrap();
            let expected_heat = power * duration * retained_fraction;
            let expected_bath = -power * duration * lost_fraction;
            let expected_power = power * (-attenuation).exp();
            assert_eq!(result.caloric.temperature_k, 500.0);
            assert_ne!(result.caloric.final_specific_energy_j_kg, 0.0);
            assert_eq!(result.heat_into_body_j.is_sign_positive(), power.is_sign_positive());
            assert_eq!(result.ambient_heat_j.is_sign_positive(), !power.is_sign_positive());
            near(result.heat_into_body_j, expected_heat, 256.0 * f64::EPSILON * expected_heat.abs());
            near(result.ambient_heat_j, expected_bath, 256.0 * f64::EPSILON * expected_bath.abs());
            near(result.source_heat_j, power * duration, 64.0 * f64::EPSILON * power.abs());
            near(
                result.heat_duration_derivative_w,
                expected_power,
                256.0 * f64::EPSILON * expected_power.abs(),
            );
            near(result.heat_source_power_derivative_s, duration * retained_fraction, 1e-13);
            near(
                result.caloric.mean_temperature_change_k,
                -expected_bath / (2.0 * duration),
                256.0 * f64::EPSILON * (expected_bath / (2.0 * duration)).abs(),
            );
            near(
                result.source_heat_j + result.ambient_heat_j,
                result.heat_into_body_j,
                128.0 * f64::EPSILON * expected_heat.abs(),
            );
            near(
                result.energy_balance_residual_j,
                0.0,
                128.0 * f64::EPSILON * expected_heat.abs(),
            );
        }
    });
}

#[test]
fn g1_passive_saturation_keeps_finite_heat_when_integrated_initial_drive_overflows() {
    let source = Fixture::new(true);
    let cv = source.cv();
    let conduction =
        IsochoricConductionLaw::try_new(&source.card, &cv, &source.conductivity, SELECTION).unwrap();
    let layer = CaloricConductionLayer::from_isochoric(&conduction).unwrap();
    with_cx(|cx| {
        for (duration, finite_initial_drive) in [(1e21, true), (1e306, false)] {
            for (initial_k, bath_k) in [(360.0, 680.0), (650.0, 320.0)] {
                let body = cell(energy(true, initial_k), bath_k);
                let exposure = (2.0 / body.mass_kg) * duration;
                assert!(exposure.is_finite());
                assert_eq!(
                    (exposure * (bath_k - initial_k)).is_finite(),
                    finite_initial_drive,
                );
                // Heat is bounded by the finite energy needed to reach the
                // supported bath, whether the initial integrated drive is
                // enormous-but-finite or requires the passive range fallback.
                let expected_heat = body.mass_kg * (energy(true, bath_k) - energy(true, initial_k));
                let result = solve(cx, &body, &layer, GATE, config(duration)).unwrap();
                assert_eq!(result.caloric.temperature_k, bath_k);
                near(result.caloric.final_specific_energy_j_kg, energy(true, bath_k), 1e-7);
                near(result.heat_into_body_j, expected_heat, 1e-7);
                near(result.stored_energy_change_j, expected_heat, 1e-7);
                near(result.ambient_heat_j, expected_heat, 1e-7);
                assert_eq!(result.reservoir_heat_j, 0.0);
                near(result.heat_initial_energy_derivative_kg, -body.mass_kg, 1e-15);
                near(
                    result.heat_ambient_temperature_derivative_j_k,
                    body.mass_kg * (200.0 + 0.5 * bath_k),
                    1e-10,
                );
                assert_eq!(result.heat_duration_derivative_w, 0.0);
                assert!(result.heat_residual_j.abs() <= config(duration).energy_tolerance_j);
                assert!(result.heat_partition_residual_j.abs() <= config(duration).energy_tolerance_j);

                if finite_initial_drive {
                    let single = solve_driven(cx, &body, &layer, GATE, config(duration)).unwrap();
                    let mut shared = body.clone();
                    shared.boundary.reservoir = Some(Reservoir {
                        area_m2: 0.0625,
                        conductance_w_k: 0.5,
                        temperature_k: bath_k,
                    });
                    let shared = solve_driven(cx, &shared, &layer, GATE, config(duration)).unwrap();
                    // Equal baths split the same G=2 into ambient G=1.5 and
                    // support G=.5. The rounded mean is at the bath, but each
                    // original bath must retain its finite share of the heat.
                    for (driven, ambient_fraction) in [(single, 1.0), (shared, 0.75)] {
                        assert_eq!(driven.caloric.temperature_k, bath_k);
                        near(driven.caloric.mean_temperature_k, bath_k, 1e-12);
                        near(driven.heat_into_body_j, expected_heat, 1e-7);
                        near(driven.stored_energy_change_j, expected_heat, 1e-7);
                        near(driven.ambient_heat_j, ambient_fraction * expected_heat, 1e-7);
                        near(
                            driven.reservoir_heat_j,
                            (1.0 - ambient_fraction) * expected_heat,
                            1e-7,
                        );
                        assert_eq!(driven.source_heat_j, 0.0);
                        near(
                            driven.ambient_heat_j + driven.reservoir_heat_j,
                            driven.heat_into_body_j,
                            config(duration).energy_tolerance_j,
                        );
                        assert!(driven.heat_partition_residual_j.abs() <= config(duration).energy_tolerance_j);
                        assert!(driven.energy_balance_residual_j.abs() <= config(duration).energy_tolerance_j);
                    }
                }
            }
            // Neither numerical range excuses original source or inventory errors.
            let unsupported = cell(energy(true, 500.0), 900.0);
            assert!(solve(cx, &unsupported, &layer, GATE, config(duration)).is_err());
            let mut wrong_inventory = cell(energy(true, 360.0), 680.0);
            wrong_inventory.mass_kg *= 1.1;
            assert!(solve(cx, &wrong_inventory, &layer, GATE, config(duration)).is_err());
        }
    });
}
