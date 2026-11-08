//! G1/G3: original ramp histories identify a passive material on an authored ladder.
use fs_material::visco::{
    GeneralizedMaxwell,
    history_fit::{
        PronyHistoryFitBudget, PronyHistoryFitError, PronyHistorySample, fit_prony_history,
        fit_prony_history_with_cancel,
    },
};
use fs_math::det;

fn law() -> GeneralizedMaxwell {
    GeneralizedMaxwell::new(1.2e9, vec![(3.5e8, 0.008), (7.0e8, 0.11), (2.3e8, 1.7)]).unwrap()
}

fn protocol(phase: f64) -> Vec<(f64, f64)> {
    let mut time = 0.0;
    let mut strain = 0.0004;
    let mut rows = vec![(time, strain)];
    let intervals = [0.003, 0.01, 0.03, 0.08, 0.2, 0.5, 0.9];
    for i in 0..100 {
        time += intervals[i % intervals.len()];
        if i % 5 != 0 {
            strain =
                0.0015 * det::sin(0.49 * i as f64 + phase) + 0.0004 * det::cos(0.17 * i as f64);
        }
        rows.push((time, strain));
    }
    rows
}

fn observations(model: &GeneralizedMaxwell, protocol: &[(f64, f64)]) -> Vec<PronyHistorySample> {
    let mut state = model.state();
    state.strain = protocol[0].1;
    state.stress = model.e_inf * state.strain;
    protocol
        .iter()
        .enumerate()
        .map(|(i, &(time_s, strain))| {
            let stress_pa = if i == 0 {
                state.stress
            } else {
                model.step(&mut state, strain, time_s - protocol[i - 1].0)
            };
            PronyHistorySample {
                time_s,
                strain,
                stress_pa,
                stress_error_pa: 0.05,
            }
        })
        .collect()
}

fn budget() -> PronyHistoryFitBudget {
    PronyHistoryFitBudget {
        max_samples: 512,
        max_branches: 8,
    }
}

fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual / expected - 1.0).abs() <= tolerance,
        "actual {actual:e}, expected {expected:e}"
    );
}

#[test]
fn g1_nonuniform_loading_reversal_recovers_passive_moduli_and_held_out_heat() {
    let source = law();
    let samples = observations(&source, &protocol(0.0));
    let fit = fit_prony_history(&samples, &[0.008, 0.11, 0.45, 1.7], budget()).unwrap();
    relative(fit.model.e_inf, source.e_inf, 2e-8);
    for (i, j) in [(0, 0), (1, 1), (3, 2)] {
        relative(fit.model.terms[i].0, source.terms[j].0, 2e-8);
    }
    assert!(
        fit.model.terms[2].0 < 1.0,
        "an absent arm must not acquire material loss"
    );
    assert!(fit.max_error_ratio <= 1.0 && fit.rms_error_ratio <= fit.max_error_ratio);
    for (prediction, observation) in fit.predicted_stress_pa.iter().zip(&samples) {
        assert!((prediction - observation.stress_pa).abs() <= observation.stress_error_pa);
    }
    // Independent experiment: different ramps, same constitutive step and true
    // heat/work ledger. No training design-matrix values enter this replay.
    let held_out = protocol(1.13);
    let mut truth = source.state();
    let mut recovered = fit.model.state();
    truth.strain = held_out[0].1;
    recovered.strain = held_out[0].1;
    for pair in held_out.windows(2) {
        let dt = pair[1].0 - pair[0].0;
        let expected = source.step(&mut truth, pair[1].1, dt);
        let actual = fit.model.step(&mut recovered, pair[1].1, dt);
        assert!((actual - expected).abs() < 0.05);
        assert!(recovered.dissipated >= 0.0);
    }
    relative(recovered.dissipated, truth.dissipated, 2e-8);
    relative(recovered.work, truth.work, 2e-8);
}

#[test]
fn g1_authored_observation_errors_weight_noisy_data_without_hiding_bad_rows() {
    let source = law();
    let mut samples = observations(&source, &protocol(0.2));
    for (i, sample) in samples.iter_mut().enumerate() {
        sample.stress_error_pa = if i == 17 { 2e5 } else { 0.1 };
        sample.stress_pa += sample.stress_error_pa
            * if i == 17 {
                0.75
            } else {
                0.2 * ((i % 3) as f64 - 1.0)
            };
    }
    let fit = fit_prony_history(&samples, &[0.008, 0.11, 1.7], budget()).unwrap();
    relative(fit.model.e_inf, source.e_inf, 2e-6);
    for (actual, expected) in fit.model.terms.iter().zip(&source.terms) {
        relative(actual.0, expected.0, 2e-6);
    }
    assert!(fit.max_error_ratio > 0.7 && fit.max_error_ratio < 0.8);
    samples[17].stress_error_pa = 0.1;
    assert!(matches!(
        fit_prony_history(&samples, &[0.008, 0.11, 1.7], budget()),
        Err(PronyHistoryFitError::StressResidual { .. })
    ));
}

#[test]
fn g3_stress_units_clock_units_and_branch_order_preserve_the_fitted_law() {
    let source = law();
    let samples = observations(&source, &protocol(0.4));
    for (stress_scale, clock_scale) in [(1e-150, 1e-6), (1.0, 1.0), (1e150, 1e6)] {
        let scaled: Vec<_> = samples
            .iter()
            .map(|s| PronyHistorySample {
                time_s: s.time_s * clock_scale,
                strain: s.strain,
                stress_pa: s.stress_pa * stress_scale,
                stress_error_pa: s.stress_error_pa * stress_scale,
            })
            .collect();
        let taus = [1.7 * clock_scale, 0.008 * clock_scale, 0.11 * clock_scale];
        let fit = fit_prony_history(&scaled, &taus, budget()).unwrap();
        relative(fit.model.e_inf / stress_scale, source.e_inf, 2e-8);
        for (i, j) in [(0, 2), (1, 0), (2, 1)] {
            relative(fit.model.terms[i].0 / stress_scale, source.terms[j].0, 2e-8);
            assert_eq!(fit.model.terms[i].1.to_bits(), taus[i].to_bits());
        }
    }
}

#[test]
fn g0_relaxed_initial_strain_and_pure_elastic_limit_are_explicit() {
    let source = GeneralizedMaxwell::new(3.7e9, vec![]).unwrap();
    let samples = observations(&source, &protocol(0.0));
    let fit = fit_prony_history(
        &samples,
        &[],
        PronyHistoryFitBudget {
            max_samples: samples.len(),
            max_branches: 0,
        },
    )
    .unwrap();
    relative(fit.model.e_inf, source.e_inf, 2e-14);
    assert!(fit.model.terms.is_empty());
    relative(
        fit.predicted_stress_pa[0],
        source.e_inf * samples[0].strain,
        2e-14,
    );

    // The same first strain with an instantaneous, UNRELAXED initial jump
    // has different stress. That prehistory is not silently inferred.
    let source = law();
    let mut samples = observations(&source, &protocol(0.0));
    samples[0].stress_pa += source.terms.iter().map(|t| t.0).sum::<f64>() * samples[0].strain;
    assert!(fit_prony_history(&samples, &[0.008, 0.11, 1.7], budget()).is_err());
}

#[test]
fn g0_unexcited_or_active_experiments_do_not_create_passive_parameters() {
    let mut samples = observations(&law(), &protocol(0.0));
    for s in &mut samples {
        s.strain = 0.0;
        s.stress_pa = 0.0;
    }
    assert_eq!(
        fit_prony_history(&samples, &[0.1], budget()).unwrap_err(),
        PronyHistoryFitError::UnexcitedColumn(0)
    );
    let samples = [
        PronyHistorySample {
            time_s: 0.0,
            strain: 0.0,
            stress_pa: 0.0,
            stress_error_pa: 0.1,
        },
        PronyHistorySample {
            time_s: 1.0,
            strain: 0.001,
            stress_pa: -100.0,
            stress_error_pa: 0.1,
        },
        PronyHistorySample {
            time_s: 2.0,
            strain: 0.002,
            stress_pa: -200.0,
            stress_error_pa: 0.1,
        },
    ];
    assert!(fit_prony_history(&samples, &[0.1], budget()).is_err());
    // A fully relaxed endpoint column cannot identify a too-fast arm.
    assert_eq!(
        fit_prony_history(&samples, &[1e-320], budget()).unwrap_err(),
        PronyHistoryFitError::UnexcitedColumn(1)
    );
}

#[test]
fn g0_admission_capacity_and_cancellation_return_no_partial_fit() {
    let samples = observations(&law(), &protocol(0.0));
    for taus in [&[0.1, 0.1][..], &[0.0], &[f64::NAN], &[f64::INFINITY]] {
        assert!(matches!(
            fit_prony_history(&samples, taus, budget()),
            Err(PronyHistoryFitError::Invalid(_))
        ));
    }
    let mut bad = samples.clone();
    bad[1].time_s = bad[0].time_s;
    assert!(fit_prony_history(&bad, &[0.1], budget()).is_err());
    bad = samples.clone();
    bad[3].stress_error_pa = 0.0;
    assert!(fit_prony_history(&bad, &[0.1], budget()).is_err());
    assert!(matches!(
        fit_prony_history(
            &samples,
            &[0.1],
            PronyHistoryFitBudget {
                max_samples: samples.len() - 1,
                max_branches: 1,
            }
        ),
        Err(PronyHistoryFitError::Capacity {
            quantity: "observations",
            ..
        })
    ));
    assert!(matches!(
        fit_prony_history(
            &samples,
            &[0.1, 0.2],
            PronyHistoryFitBudget {
                max_samples: samples.len(),
                max_branches: 1,
            }
        ),
        Err(PronyHistoryFitError::Capacity {
            quantity: "Maxwell branches",
            ..
        })
    ));
    let mut checks = 0;
    let result = fit_prony_history_with_cancel(&samples, &[0.1], budget(), || {
        checks += 1;
        checks == 7
    });
    assert_eq!(result.unwrap_err(), PronyHistoryFitError::Cancelled);
}
