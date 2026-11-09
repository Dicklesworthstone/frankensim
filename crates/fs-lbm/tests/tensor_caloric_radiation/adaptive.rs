//! G1/G3/G4: adaptive radiation keeps original phase energy and accepted heat.
use super::*;
use fs_lbm::conjugate::{
    TensorCaloricAdaptiveConfig, TensorCaloricAdaptiveError, TensorCaloricForcing,
    advance_caloric_tensor_conduction_radiating_adaptive,
};

fn adaptive(absolute_j_kg: f64) -> TensorCaloricAdaptiveConfig {
    TensorCaloricAdaptiveConfig {
        step: controls(),
        initial_dt_s: 0.5,
        minimum_dt_s: 1e-8,
        maximum_dt_s: 0.5,
        maximum_attempts: 4096,
        energy_absolute_tolerance_j_kg: absolute_j_kg,
        relative_tolerance: 0.0,
        energy_tolerance_j: 1e-5,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

/// Independent explicit RK4 reference for the single finite radiating body.
/// This uses the scalar physical law, not either native SDIR stage or estimator.
fn cooling_reference(duration: f64, count: usize) -> f64 {
    let rate = |t: f64| -((t - 300.0) + radiant_power(t, 300.0)) / (MASS * CP);
    let dt = duration / count as f64;
    let mut temperature = 400.0;
    for _ in 0..count {
        let a = rate(temperature);
        let b = rate(temperature + 0.5 * dt * a);
        let c = rate(temperature + 0.5 * dt * b);
        let d = rate(temperature + dt * c);
        temperature += dt * (a + 2.0 * b + 2.0 * c + d) / 6.0;
    }
    temperature
}

fn shifted_latent(datum: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        ContentHash([0xD3; 32]),
        [
            (-200.0, 200.0, 0.0),
            (0.0, 300.0, 0.0),
            (100.0, 300.0, 1.0),
            (500.0, 500.0, 1.0),
        ]
        .map(|(h, t, fraction)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: datum + h,
            temperature_k: t,
            liquid_mass_fraction: fraction,
            bulk_density_kg_m3: MASS,
        })
        .to_vec(),
    )
    .unwrap()
}

#[test]
fn adaptive_radiative_cooling_refines_against_an_independent_scalar_flow() {
    let (domain, solids, setup, radiation) = exposed_pair(300.0);
    let law = constant_cp(0.0);
    let laws: [&dyn TensorCaloricLaw; 2] = [&law, &law];
    let k = conductivity();
    let previous = [200.0, 0.0];
    let reference = cooling_reference(0.5, 20_000);
    near(reference, cooling_reference(0.5, 10_000), 1e-9);
    let mut errors = Vec::new();
    let mut counts = Vec::new();
    for tolerance in [0.2, 0.002] {
        let mut calls = 0;
        let result = advance_caloric_tensor_conduction_radiating_adaptive(
            &domain,
            &solids,
            &[&k, &k],
            TensorCaloricStorage {
                laws: &laws,
                reference_mass_kg: &[MASS; 2],
                previous_specific_energy_j_kg: &previous,
            },
            0.0,
            0.5,
            &[],
            &adaptive(tolerance),
            &radiation,
            |sample| {
                calls += 1;
                assert!(!sample.left_limit);
                Ok(TensorCaloricForcing::borrowed(&setup, &[(1, 0.0)]))
            },
            &CancelGate::new(),
        )
        .unwrap();
        let thermal = &result.thermal;
        assert_eq!(calls, 2 * thermal.report.attempts);
        assert!(thermal.report.accuracy_rejections > 0);
        assert_eq!(thermal.report.solver_rejections, 0);
        assert_eq!(thermal.records.len(), thermal.report.accepted_steps);
        assert_eq!(result.radiation.len(), thermal.records.len());
        assert!(thermal.report.maximum_accepted_normalized_error <= 1.0);
        assert!(
            thermal
                .records
                .windows(2)
                .all(|r| r[0].time_s < r[1].time_s)
        );
        let mut radiated_energy = 0.0;
        let mut prior_time = 0.0;
        let mut accepted_solves = 0;
        for (record, reports) in thermal.records.iter().zip(&result.radiation) {
            let dt = record.time_s - prior_time;
            radiated_energy += dt
                * ((1.0 - TENSOR_SDIRK2_GAMMA) * reports[0].radiated_w
                    + TENSOR_SDIRK2_GAMMA * reports[1].radiated_w);
            accepted_solves += reports.iter().map(|r| r.iterations).sum::<usize>();
            prior_time = record.time_s;
        }
        near(result.radiated_energy_j, radiated_energy, 1e-9);
        assert!(result.completed_trial_radiation_solves > accepted_solves);
        let temperature = thermal.endpoint.solution.solution.temperature[0];
        near(
            thermal.record.stored_energy_change_j,
            MASS * CP * (temperature - 400.0),
            1e-8,
        );
        near(thermal.record.source_j, 0.0, 1e-12);
        near(
            thermal.record.boundary_outflow_j,
            -thermal.record.stored_energy_change_j,
            1e-6,
        );
        near(thermal.record.closure_j, 0.0, 1e-6);
        near(
            thermal.record.radiated_w,
            radiant_power(temperature, 300.0),
            1e-7,
        );
        assert_eq!(thermal.record.time_s, 0.5);
        errors.push((temperature - reference).abs());
        counts.push(thermal.report.accepted_steps);
    }
    assert!(errors[1] < errors[0] / 5.0, "{errors:?}");
    assert!(errors[1] < 0.01, "{errors:?}");
    assert!(counts[1] > counts[0]);
    assert_eq!(previous, [200.0, 0.0]);
}

#[test]
fn latent_energy_controls_adaptation_while_radiation_and_temperature_stay_constant() {
    let (domain, solids, _, radiation) = exposed_pair(250.0);
    let k = conductivity();
    let power = radiant_power(300.0, 250.0);
    for datum in [0.0, 1e6] {
        let law = shifted_latent(datum);
        let laws: [&dyn TensorCaloricLaw; 2] = [&law, &law];
        let previous = [datum + 20.0, datum];
        let mut config = adaptive(0.01);
        config.initial_dt_s = 1.0;
        config.maximum_dt_s = 1.0;
        config.step.linear.tolerance = 1e-8;
        let result = advance_caloric_tensor_conduction_radiating_adaptive(
            &domain,
            &solids,
            &[&k, &k],
            TensorCaloricStorage {
                laws: &laws,
                reference_mass_kg: &[MASS; 2],
                previous_specific_energy_j_kg: &previous,
            },
            0.0,
            1.0,
            &[],
            &config,
            &radiation,
            |sample| {
                let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
                setup.power_w = vec![power + 20.0 * sample.time_s, 0.0];
                setup.fixed_temperature = vec![(1, 300.0)];
                Ok(TensorCaloricForcing::owned(setup, vec![(1, datum)]))
            },
            &CancelGate::new(),
        )
        .unwrap();
        let thermal = &result.thermal;
        assert!(thermal.report.accuracy_rejections > 0);
        assert!(thermal.report.accepted_steps > 1);
        assert!(thermal.records.iter().all(|r| r.max_temperature_k == 300.0));
        assert_eq!(
            thermal.endpoint.solution.solution.temperature,
            vec![300.0; 2]
        );
        near(thermal.endpoint.specific_energy_j_kg[0] - datum, 25.0, 2e-5);
        near(thermal.record.source_j, power + 10.0, 1e-7);
        near(thermal.record.boundary_outflow_j, power, 1e-7);
        near(result.radiated_energy_j, power, 1e-7);
        near(thermal.record.stored_energy_change_j, 10.0, 2e-5);
        near(thermal.record.closure_j, 0.0, config.energy_tolerance_j);
        near(
            law.values_at_specific_enthalpy(thermal.endpoint.specific_energy_j_kg[0])
                .unwrap()
                .liquid_mass_fraction,
            0.25,
            2e-7,
        );
        assert_eq!(previous, [datum + 20.0, datum]);
    }
}

#[test]
fn transition_left_limits_keep_prescribed_phase_heat_and_radiation_separate() {
    let (domain, solids, _, radiation) = exposed_pair(250.0);
    let law = melting();
    let laws: [&dyn TensorCaloricLaw; 2] = [&law, &law];
    let k = conductivity();
    let previous = [80.0, 0.0];
    let power = radiant_power(300.0, 250.0);
    let mut config = adaptive(0.001);
    config.initial_dt_s = 1.0;
    config.maximum_dt_s = 1.0;
    let mut left_limits = 0;
    let result = advance_caloric_tensor_conduction_radiating_adaptive(
        &domain,
        &solids,
        &[&k, &k],
        TensorCaloricStorage {
            laws: &laws,
            reference_mass_kg: &[MASS; 2],
            previous_specific_energy_j_kg: &previous,
        },
        0.0,
        1.0,
        &[0.5],
        &config,
        &radiation,
        |sample| {
            assert!(sample.interval_s[1] <= 0.5 || sample.interval_s[0] >= 0.5);
            let before = sample.time_s < 0.5 || (sample.time_s == 0.5 && sample.left_limit);
            if sample.left_limit {
                left_limits += 1;
                assert_eq!(sample.time_s, 0.5);
            }
            let mut setup = ThermalSetup::new([ThermalFace::Adiabatic; 6]);
            setup.power_w = vec![power + 4.0, 0.0];
            setup.fixed_temperature = vec![(1, 300.0)];
            let mut prescribed = vec![(1, 0.0)];
            if !before {
                setup.fixed_temperature.push((0, 300.0));
                prescribed.push((0, 90.0));
            }
            Ok(TensorCaloricForcing::owned(setup, prescribed))
        },
        &CancelGate::new(),
    )
    .unwrap();
    let thermal = &result.thermal;
    assert_eq!(left_limits, 1);
    assert_eq!(thermal.records.len(), 2);
    assert_eq!(thermal.records[0].time_s, 0.5);
    assert_eq!(thermal.records[1].time_s, 1.0);
    assert_eq!(thermal.endpoint.specific_energy_j_kg, vec![90.0, 0.0]);
    assert_eq!(
        thermal.endpoint.solution.solution.temperature,
        vec![300.0; 2]
    );
    near(thermal.records[0].stored_energy_change_j, 2.0, 1e-8);
    near(thermal.records[1].stored_energy_change_j, 18.0, 1e-8);
    near(thermal.record.source_j, 0.5 * power + 2.0, 1e-8);
    near(thermal.record.boundary_outflow_j, 0.5 * power - 18.0, 1e-8);
    near(thermal.record.stored_energy_change_j, 20.0, 1e-8);
    near(thermal.record.closure_j, 0.0, 1e-8);
    near(result.radiated_energy_j, power, 1e-8);
    assert_eq!(previous, [80.0, 0.0]);
}

#[test]
fn exhausted_work_and_late_refusals_publish_no_adaptive_radiation_prefix() {
    let (domain, solids, mut setup, radiation) = exposed_pair(250.0);
    let law = melting();
    let laws: [&dyn TensorCaloricLaw; 2] = [&law, &law];
    let k = conductivity();
    let previous = [80.0, 0.0];
    let power = radiant_power(300.0, 250.0);
    setup.power_w = vec![power + 1.0, 0.0];
    let saved_setup = setup.clone();
    let storage = TensorCaloricStorage {
        laws: &laws,
        reference_mass_kg: &[MASS; 2],
        previous_specific_energy_j_kg: &previous,
    };
    let mut config = adaptive(0.01);
    config.initial_dt_s = 0.25;
    config.maximum_dt_s = 0.25;
    let mut one_attempt = config;
    one_attempt.maximum_attempts = 1;
    let limited = advance_caloric_tensor_conduction_radiating_adaptive(
        &domain,
        &solids,
        &[&k, &k],
        storage,
        0.0,
        1.0,
        &[],
        &one_attempt,
        &radiation,
        |_| Ok(TensorCaloricForcing::borrowed(&setup, &[(1, 0.0)])),
        &CancelGate::new(),
    );
    assert!(matches!(
        limited,
        Err(TensorCaloricAdaptiveError::AttemptBudget {
            time_s: 0.25,
            attempts: 1,
        })
    ));
    let mut one_heat_solve = radiation.clone();
    one_heat_solve.max_iterations = 1;
    let mut floor = config;
    floor.minimum_dt_s = 0.25;
    assert!(matches!(
        advance_caloric_tensor_conduction_radiating_adaptive(
            &domain,
            &solids,
            &[&k, &k],
            storage,
            0.0,
            1.0,
            &[],
            &floor,
            &one_heat_solve,
            |_| Ok(TensorCaloricForcing::borrowed(&setup, &[(1, 0.0)])),
            &CancelGate::new(),
        ),
        Err(TensorCaloricAdaptiveError::Solve(
            ChtError::SolverNotConverged { .. }
        ))
    ));

    for cancel in [false, true] {
        let gate = CancelGate::new();
        let mut calls = 0;
        let result = advance_caloric_tensor_conduction_radiating_adaptive(
            &domain,
            &solids,
            &[&k, &k],
            storage,
            0.0,
            1.0,
            &[],
            &config,
            &radiation,
            |sample| {
                calls += 1;
                if sample.interval_s[0] >= 0.25 {
                    if cancel {
                        gate.request();
                    } else {
                        return Err(ChtError::InvalidInput {
                            field: "test.radiation_source",
                            reason: "physical source support ended".into(),
                        });
                    }
                }
                Ok(TensorCaloricForcing::borrowed(&setup, &[(1, 0.0)]))
            },
            &gate,
        );
        assert_eq!(calls, 3);
        if cancel {
            assert!(matches!(
                result,
                Err(TensorCaloricAdaptiveError::Solve(ChtError::Cancelled))
            ));
        } else {
            assert!(matches!(
                result,
                Err(TensorCaloricAdaptiveError::Solve(
                    ChtError::InvalidInput { .. }
                ))
            ));
        }
    }
    assert_eq!(setup, saved_setup);
    assert_eq!(previous, [80.0, 0.0]);
}

#[test]
fn fatal_caloric_refusal_in_a_newton_trial_preserves_the_original_diagnostic() {
    use fs_conduction::enthalpy::{DepositionError, ThermalLaw};
    use std::cell::Cell;

    struct TrialRefusal<'a> {
        original: &'a dyn TensorCaloricLaw,
        refusals: Cell<usize>,
        refused_energy: Cell<f64>,
    }

    impl ThermalLaw for TrialRefusal<'_> {
        fn temperature_k(&self, energy: f64) -> Result<f64, DepositionError> {
            let temperature = self.original.temperature_k(energy)?;
            // The original opening h=0 and both source bounds remain admitted.
            // Only a positive interior Newton trial triggers this source budget.
            if energy > 0.0 && energy < 100.0 {
                self.refusals.set(self.refusals.get() + 1);
                self.refused_energy.set(energy);
                return Err(DepositionError::Budget("trial caloric work exhausted"));
            }
            Ok(temperature)
        }

        fn maximum_temperature_slope(&self) -> f64 {
            self.original.maximum_temperature_slope()
        }
    }

    impl TensorCaloricLaw for TrialRefusal<'_> {
        fn specific_energy_bounds_j_kg(&self) -> [f64; 2] {
            self.original.specific_energy_bounds_j_kg()
        }

        fn temperature_slope_kg_k_j(&self, energy: f64) -> Result<f64, DepositionError> {
            self.original.temperature_slope_kg_k_j(energy)
        }
    }

    let (domain, solids, mut setup, radiation) = exposed_pair(300.0);
    setup.power_w = vec![80.0, 0.0];
    let original = constant_cp(0.0);
    let guarded = TrialRefusal {
        original: &original,
        refusals: Cell::new(0),
        refused_energy: Cell::new(f64::NAN),
    };
    let laws: [&dyn TensorCaloricLaw; 2] = [&guarded, &original];
    let k = conductivity();
    let previous = [0.0; 2];
    let mut config = adaptive(0.01);
    config.initial_dt_s = 0.25;
    config.minimum_dt_s = 0.25;
    config.maximum_dt_s = 0.25;
    config.maximum_attempts = 4;
    config.step.max_backtracks = 4;
    let mut samples = 0;
    let refusal = advance_caloric_tensor_conduction_radiating_adaptive(
        &domain,
        &solids,
        &[&k, &k],
        TensorCaloricStorage {
            laws: &laws,
            reference_mass_kg: &[MASS; 2],
            previous_specific_energy_j_kg: &previous,
        },
        0.0,
        0.25,
        &[],
        &config,
        &radiation,
        |_| {
            samples += 1;
            Ok(TensorCaloricForcing::borrowed(&setup, &[(1, 0.0)]))
        },
        &CancelGate::new(),
    )
    .unwrap_err();
    match refusal {
        TensorCaloricAdaptiveError::Solve(ChtError::InvalidInput { field, reason }) => {
            assert_eq!(field, "thermal.tensor_caloric.original_source_refusal");
            assert_eq!(
                reason,
                DepositionError::Budget("trial caloric work exhausted").to_string(),
            );
        }
        other => panic!("original source refusal was swallowed: {other:?}"),
    }
    assert_eq!(
        samples, 2,
        "one attempt samples exactly two physical stages"
    );
    assert_eq!(
        guarded.refusals.get(),
        1,
        "fatal source work is never retried"
    );
    // At T=300 K the first radiative tangent has zero offset power. Its
    // positive-energy Newton correction has this independent scalar value.
    let a = TENSOR_SDIRK2_GAMMA * config.initial_dt_s;
    let radiation_slope = 4.0 * EMISSIVITY * STEFAN_BOLTZMANN * 300.0_f64.powi(3);
    let first_trial = 80.0 / (MASS / a + (1.0 + radiation_slope) / CP);
    near(guarded.refused_energy.get(), first_trial, 1e-8);
    assert_eq!(previous, [0.0; 2]);
    assert_eq!(setup.power_w, [80.0, 0.0]);
}
