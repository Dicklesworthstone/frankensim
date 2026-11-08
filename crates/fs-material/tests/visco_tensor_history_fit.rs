//! G1/G3/G4: joint full-tensor calibration, original material replay and refusals.
use fs_material::visco_tensor::history_fit::{
    TensorPronyHistoryFitBudget as Budget, TensorPronyHistoryFitError as Error,
    TensorPronyHistorySample as Sample, TensorPronyHistoryTerm as Term,
    fit_tensor_prony_histories as fit, fit_tensor_prony_histories_with_cancel as fit_cancel,
};
use fs_material::visco_tensor::{TensorProny, TensorPronyBranch};
use fs_material::{Tangent6, Voigt};

const W: Voigt = [
    1.0,
    1.0,
    1.0,
    std::f64::consts::SQRT_2,
    std::f64::consts::SQRT_2,
    std::f64::consts::SQRT_2,
];

fn budget() -> Budget {
    Budget {
        max_experiments: 8,
        max_samples: 512,
        max_terms: 16,
    }
}
fn rank_one(row: Voigt) -> Tangent6 {
    let mut result = [[0.0; 6]; 6];
    result[0] = row;
    result
}
fn factor(term: &Term) -> Tangent6 {
    match term {
        Term::Equilibrium { gram_rows }
        | Term::Maxwell { gram_rows, .. }
        | Term::Kelvin { gram_rows } => *gram_rows,
    }
}
fn gram(factor: &Tangent6) -> Tangent6 {
    core::array::from_fn(|i| {
        core::array::from_fn(|j| (0..6).map(|k| factor[k][i] * factor[k][j]).sum())
    })
}
fn apply(a: &Tangent6, x: Voigt) -> Voigt {
    core::array::from_fn(|i| (0..6).map(|j| a[i][j] * x[j]).sum())
}
fn mandel(x: Voigt) -> Voigt {
    core::array::from_fn(|i| x[i] * W[i])
}
fn tensor(x: Voigt) -> Voigt {
    core::array::from_fn(|i| x[i] / W[i])
}
fn dot(a: Voigt, b: Voigt) -> f64 {
    (0..6).map(|i| a[i] * b[i]).sum()
}
fn relative(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual / expected - 1.0).abs() < tolerance,
        "actual={actual:e}, expected={expected:e}"
    );
}

fn terms() -> Vec<Term> {
    vec![
        Term::Equilibrium {
            gram_rows: core::array::from_fn(|i| {
                core::array::from_fn(|j| if i == j { 0.8 + 0.1 * i as f64 } else { 0.0 })
            }),
        },
        Term::Equilibrium {
            gram_rows: rank_one([0.8, -0.6, 0.2, 0.7, 0.0, 0.1]),
        },
        Term::Maxwell {
            gram_rows: rank_one([1.0, 0.2, 0.0, -0.5, 0.0, 0.0]),
            relaxation_time_s: 0.4,
        },
        Term::Maxwell {
            gram_rows: rank_one([0.0, 0.7, -0.3, 0.2, 0.6, 0.1]),
            relaxation_time_s: 1.7,
        },
        Term::Kelvin {
            gram_rows: rank_one([0.3, -0.2, 0.4, 0.1, 0.7, -0.5]),
        },
        Term::Maxwell {
            gram_rows: rank_one([0.4, 0.1, 0.5, 0.7, -0.2, 0.0]),
            relaxation_time_s: 0.9,
        },
    ]
}
fn strengths() -> Vec<f64> {
    vec![1200.0, 370.0, 650.0, 230.0, 17.0, 0.0]
}

fn protocol(phase: f64, initial: bool) -> Vec<(f64, Voigt)> {
    let mut time = 0.0;
    let mut rows = Vec::new();
    for i in 0..38 {
        if i > 0 {
            time += [0.015, 0.06, 0.11, 0.4, 1.2][(i - 1) % 5];
        }
        let strain = core::array::from_fn(|k| {
            if i == 0 && !initial {
                0.0
            } else {
                0.001 * (0.37 * i as f64 + 0.71 * k as f64 + phase).sin()
                    + 0.0004 * (0.19 * i as f64 - 0.41 * k as f64).cos()
            }
        });
        rows.push((time, strain));
    }
    rows
}

// Independent analytic strain-memory recurrence and assembled ORIGINAL Gram
// tensors. No fitting design matrix or PronyRamp implementation enters this oracle.
fn observations(terms: &[Term], strengths: &[f64], protocol: &[(f64, Voigt)]) -> Vec<Sample> {
    let mut memories = vec![[0.0; 6]; terms.len()];
    let mut output = Vec::new();
    for (i, &(time_s, strain)) in protocol.iter().enumerate() {
        let e = mandel(strain);
        let (dt, delta) = if i == 0 {
            (0.0, [0.0; 6])
        } else {
            (
                time_s - protocol[i - 1].0,
                core::array::from_fn(|k| e[k] - mandel(protocol[i - 1].1)[k]),
            )
        };
        let mut stress = [0.0; 6];
        for (j, (term, &strength)) in terms.iter().zip(strengths).enumerate() {
            let argument = match term {
                Term::Equilibrium { .. } => e,
                Term::Maxwell {
                    relaxation_time_s: tau,
                    ..
                } => {
                    if i > 0 {
                        let r = dt / tau;
                        let decay = (-r).exp();
                        let ramp = -(-r).exp_m1() / r;
                        for (q, d) in memories[j].iter_mut().zip(delta) {
                            *q = decay * *q + ramp * d;
                        }
                    }
                    memories[j]
                }
                Term::Kelvin { .. } => {
                    if i == 0 {
                        [0.0; 6]
                    } else {
                        delta.map(|v| v / dt)
                    }
                }
            };
            let sigma = apply(&gram(&factor(term)), argument);
            for k in 0..6 {
                stress[k] += strength * sigma[k];
            }
        }
        output.push(Sample {
            time_s,
            strain,
            stress_pa: tensor(stress),
            stress_error_pa: [2e-8; 6],
        });
    }
    output
}

fn original_law(terms: &[Term], strengths: &[f64]) -> TensorProny {
    let mut eq = [[0.0; 6]; 6];
    let mut kv = [[0.0; 6]; 6];
    let mut branches = Vec::new();
    for (term, &strength) in terms.iter().zip(strengths) {
        match term {
            Term::Maxwell {
                gram_rows,
                relaxation_time_s,
            } => branches.push(
                TensorPronyBranch::from_gram_rows(
                    gram_rows.map(|r| r.map(|v| strength.sqrt() * v)),
                    *relaxation_time_s,
                )
                .unwrap(),
            ),
            Term::Equilibrium { .. } | Term::Kelvin { .. } => {
                let c = gram(&factor(term));
                let out = if matches!(term, Term::Equilibrium { .. }) {
                    &mut eq
                } else {
                    &mut kv
                };
                for i in 0..6 {
                    for j in 0..6 {
                        out[i][j] += strength * c[i][j];
                    }
                }
            }
        }
    }
    TensorProny::new_mandel(eq, branches, 0.02)
        .unwrap()
        .with_kelvin_viscosity_mandel(kv)
        .unwrap()
}

// Continuous original material heat integral, independently integrated with
// composite Simpson quadrature; Kelvin and every tensor Maxwell arm stay separate.
fn heat_integral(terms: &[Term], strengths: &[f64], protocol: &[(f64, Voigt)]) -> f64 {
    let mut memory = vec![[0.0; 6]; terms.len()];
    let mut total = 0.0;
    for segment in protocol.windows(2) {
        let dt = segment[1].0 - segment[0].0;
        let a = mandel(segment[0].1);
        let b = mandel(segment[1].1);
        let rate: Voigt = core::array::from_fn(|k| (b[k] - a[k]) / dt);
        for (j, (term, &strength)) in terms.iter().zip(strengths).enumerate() {
            let c = gram(&factor(term));
            match term {
                Term::Equilibrium { .. } => {}
                Term::Kelvin { .. } => total += dt * strength * dot(rate, apply(&c, rate)),
                Term::Maxwell {
                    relaxation_time_s: tau,
                    ..
                } => {
                    let start = memory[j];
                    let evaluate = |t: f64| {
                        let decay = (-t / tau).exp();
                        let ramp = -(-t / tau).exp_m1();
                        core::array::from_fn(|k| decay * start[k] + tau * ramp * rate[k])
                    };
                    let mut integral = 0.0;
                    for k in 0..=256 {
                        let q = evaluate(dt * f64::from(k) / 256.0);
                        let weight = if k == 0 || k == 256 {
                            1.0
                        } else if k % 2 == 0 {
                            2.0
                        } else {
                            4.0
                        };
                        integral += weight * dot(q, apply(&c, q));
                    }
                    total += strength * dt * integral / (768.0 * tau);
                    memory[j] = evaluate(dt);
                }
            }
        }
    }
    total
}

#[test]
fn g1_joint_tensor_fit_recovers_coupled_material_and_independent_held_out_heat() {
    let terms = terms();
    let strengths = strengths();
    let first = observations(&terms, &strengths, &protocol(0.0, true));
    let second = observations(&terms, &strengths, &protocol(1.37, true));
    let result = fit(&[&first, &second], &terms, 0.02, budget()).unwrap();
    for (&actual, &expected) in result.strengths.iter().zip(&strengths).take(5) {
        relative(actual, expected, 3e-9);
    }
    assert!(result.strengths[5] < 1e-6);
    assert_eq!(result.model.branches().len(), 3);
    assert!(result.model.branches()[0].stiffness_mandel()[0][3] < 0.0);
    assert!(result.max_error_ratio <= 1.0 && result.rms_error_ratio <= result.max_error_ratio);
    assert_eq!(result.predicted_stress_pa[0].len(), first.len());
    assert_eq!(result.predicted_stress_pa[1].len(), second.len());
    let held_out = protocol(2.71, false);
    let original = original_law(&terms, &strengths);
    let mut expected = original.history();
    let mut actual = result.model.history();
    for pair in held_out.windows(2) {
        let dt = pair[1].0 - pair[0].0;
        let truth = expected.advance(pair[1].1, dt, || false).unwrap();
        let replay = actual.advance(pair[1].1, dt, || false).unwrap();
        for k in 0..6 {
            assert!((truth.stress[k] - replay.stress[k]).abs() < 2e-8);
        }
        assert!(replay.dissipated >= 0.0);
        assert!(replay.balance_residual.abs() < 2e-14);
    }
    relative(
        actual.dissipated(),
        heat_integral(&terms, &strengths, &held_out),
        2e-8,
    );
    relative(actual.work(), expected.work(), 3e-9);
    relative(actual.stored_energy(), expected.stored_energy(), 3e-9);
}

#[test]
fn g1_separate_experiments_supply_missing_tensor_directions_and_reset_memory() {
    let identity =
        core::array::from_fn(|i| core::array::from_fn(|j| if i == j { 1.0 } else { 0.0 }));
    let terms = [
        Term::Equilibrium {
            gram_rows: identity,
        },
        Term::Maxwell {
            gram_rows: rank_one([1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
            relaxation_time_s: 0.7,
        },
        Term::Maxwell {
            gram_rows: rank_one([0.0, 1.0, 0.0, 0.0, 0.0, 0.0]),
            relaxation_time_s: 0.7,
        },
        Term::Kelvin {
            gram_rows: rank_one([0.0, 0.0, 0.0, 1.0, 0.0, 0.0]),
        },
    ];
    let strengths = [900.0, 370.0, 820.0, 27.0];
    let experiments: Vec<_> = [0, 1, 3]
        .into_iter()
        .map(|axis| {
            let p: Vec<_> = protocol(axis as f64, true)
                .iter()
                .map(|&(t, e)| {
                    let mut uniaxial = [0.0; 6];
                    uniaxial[axis] = e[axis];
                    (t, uniaxial)
                })
                .collect();
            observations(&terms, &strengths, &p)
        })
        .collect();
    assert!(matches!(
        fit(&[&experiments[0]], &terms, 0.02, budget()),
        Err(Error::UnexcitedTerm(2))
    ));
    let refs: Vec<_> = experiments.iter().map(Vec::as_slice).collect();
    let result = fit(&refs, &terms, 0.02, budget()).unwrap();
    for (actual, expected) in result.strengths.iter().zip(strengths) {
        relative(*actual, expected, 2e-10);
    }
    // Every first stress has only relaxed equilibrium, even when each initial
    // strain is nonzero and the immediately preceding experiment ended loaded.
    for (predicted, source) in result.predicted_stress_pa.iter().zip(&experiments) {
        for (&actual, expected) in predicted[0].iter().zip(source[0].stress_pa) {
            assert!((actual - expected).abs() < 2e-8);
        }
    }
}

fn rotation() -> Tangent6 {
    let s = 0.5_f64.sqrt();
    [
        [0.5, 0.5, 0.0, -s, 0.0, 0.0],
        [0.5, 0.5, 0.0, s, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.0, 0.0],
        [s, -s, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, s, s],
        [0.0, 0.0, 0.0, 0.0, -s, s],
    ]
}

#[test]
fn g3_joint_material_rotation_preserves_fitted_strengths_and_signed_shear() {
    let basis = terms();
    let source = observations(&basis, &strengths(), &protocol(0.23, true));
    let r = rotation();
    let rotated: Vec<_> = basis
        .iter()
        .map(|term| {
            let b = factor(term);
            let gram_rows = core::array::from_fn(|i| {
                core::array::from_fn(|j| (0..6).map(|k| b[i][k] * r[j][k]).sum())
            });
            match term {
                Term::Equilibrium { .. } => Term::Equilibrium { gram_rows },
                Term::Maxwell {
                    relaxation_time_s, ..
                } => Term::Maxwell {
                    gram_rows,
                    relaxation_time_s: *relaxation_time_s,
                },
                Term::Kelvin { .. } => Term::Kelvin { gram_rows },
            }
        })
        .collect();
    let samples: Vec<_> = source
        .iter()
        .map(|s| Sample {
            strain: tensor(apply(&r, mandel(s.strain))),
            stress_pa: tensor(apply(&r, mandel(s.stress_pa))),
            ..*s
        })
        .collect();
    let original = fit(&[&source], &basis, 0.02, budget()).unwrap();
    let actual = fit(&[&samples], &rotated, 0.02, budget()).unwrap();
    for j in 0..5 {
        relative(actual.strengths[j], original.strengths[j], 2e-9);
    }
}

#[test]
fn g3_stress_clock_and_authored_factor_normalizations_preserve_the_material() {
    let original = terms();
    let original_strengths = strengths();
    for (stress_scale, clock_scale) in [(1e-120, 1e-6), (1e120, 1e6)] {
        let factors = [0.01, 7.0, 0.3, 4.0, 0.02, 1.3];
        let terms: Vec<_> = original
            .iter()
            .enumerate()
            .map(|(j, term)| {
                let gram_rows = factor(term).map(|r| r.map(|x| x * factors[j]));
                match term {
                    Term::Equilibrium { .. } => Term::Equilibrium { gram_rows },
                    Term::Maxwell {
                        relaxation_time_s, ..
                    } => Term::Maxwell {
                        gram_rows,
                        relaxation_time_s: relaxation_time_s * clock_scale,
                    },
                    Term::Kelvin { .. } => Term::Kelvin { gram_rows },
                }
            })
            .collect();
        let expected: Vec<_> = original_strengths
            .iter()
            .enumerate()
            .map(|(j, s)| {
                s * stress_scale
                    * if matches!(original[j], Term::Kelvin { .. }) {
                        clock_scale
                    } else {
                        1.0
                    }
                    / (factors[j] * factors[j])
            })
            .collect();
        let p: Vec<_> = protocol(0.71, true)
            .iter()
            .map(|&(t, e)| (t * clock_scale, e))
            .collect();
        let mut samples = observations(&terms, &expected, &p);
        for s in &mut samples {
            s.stress_error_pa = [2e-8 * stress_scale; 6];
        }
        let actual = fit(&[&samples], &terms, 0.02, budget()).unwrap();
        for (&actual, expected) in actual.strengths.iter().zip(expected).take(5) {
            relative(actual, expected, 3e-9);
        }
    }
}

#[test]
fn g1_component_errors_weight_noise_and_replay_rejects_an_incompatible_row() {
    let terms = terms();
    let mut samples = observations(&terms, &strengths(), &protocol(0.17, true));
    samples[17].stress_error_pa[4] = 0.5;
    samples[17].stress_pa[4] += 0.3;
    let result = fit(&[&samples], &terms, 0.02, budget()).unwrap();
    assert!((result.max_error_ratio - 0.6).abs() < 1e-5);
    samples[17].stress_error_pa[4] = 2e-8;
    assert!(matches!(
        fit(&[&samples], &terms, 0.02, budget()),
        Err(Error::StressResidual { .. })
    ));
}

#[test]
fn g0_passivity_and_full_equilibrium_cannot_be_inferred_from_inadmissible_data() {
    let terms = terms();
    let mut samples = observations(&terms, &strengths(), &protocol(0.41, true));
    for sample in &mut samples {
        sample.stress_pa = sample.stress_pa.map(|s| -s);
    }
    assert!(fit(&[&samples], &terms, 0.02, budget()).is_err());
    let singular = [Term::Equilibrium {
        gram_rows: rank_one([1.0, -0.4, 0.0, 0.2, 0.0, 0.0]),
    }];
    let data = observations(&singular, &[1000.0], &protocol(0.0, true));
    assert!(matches!(
        fit(&[&data], &singular, 0.02, budget()),
        Err(Error::Material(_))
    ));
}

#[test]
#[allow(clippy::too_many_lines)] // One admission/cancellation/retry scenario.
fn g4_admission_capacities_and_late_cancellation_publish_no_partial_fit() {
    let terms = terms();
    let samples = observations(&terms, &strengths(), &protocol(0.53, true));
    assert!(matches!(
        fit(
            &[&samples],
            &terms,
            0.02,
            Budget {
                max_terms: 5,
                ..budget()
            }
        ),
        Err(Error::Capacity {
            quantity: "tensor terms",
            ..
        })
    ));
    assert!(matches!(
        fit(
            &[&samples],
            &terms,
            0.02,
            Budget {
                max_samples: samples.len() - 1,
                ..budget()
            }
        ),
        Err(Error::Capacity {
            quantity: "total observations",
            ..
        })
    ));
    assert!(matches!(
        fit(
            &[&samples, &samples],
            &terms,
            0.02,
            Budget {
                max_experiments: 1,
                ..budget()
            }
        ),
        Err(Error::Capacity {
            quantity: "experiments",
            ..
        })
    ));
    let too_many = vec![terms[0]; 65];
    assert!(matches!(
        fit(
            &[&samples],
            &too_many,
            0.02,
            Budget {
                max_terms: 100,
                ..budget()
            }
        ),
        Err(Error::Capacity {
            quantity: "tensor terms",
            limit: 64,
            ..
        })
    ));
    for invalid in [f64::NAN, 0.0, -1.0] {
        assert!(matches!(
            fit(&[&samples], &terms, invalid, budget()),
            Err(Error::Invalid(_))
        ));
    }
    let mut bad = samples.clone();
    bad[3].time_s = bad[2].time_s;
    assert!(matches!(
        fit(&[&bad], &terms, 0.02, budget()),
        Err(Error::Invalid(_))
    ));
    bad = samples.clone();
    bad[3].stress_error_pa[1] = 0.0;
    assert!(matches!(
        fit(&[&bad], &terms, 0.02, budget()),
        Err(Error::Invalid(_))
    ));
    let mut polls = 0;
    let baseline = fit_cancel(&[&samples], &terms, 0.02, budget(), || {
        polls += 1;
        false
    })
    .unwrap();
    for stop in [1, polls / 2, polls - 1] {
        let mut count = 0;
        assert!(matches!(
            fit_cancel(&[&samples], &terms, 0.02, budget(), || {
                count += 1;
                count >= stop
            }),
            Err(Error::Cancelled)
        ));
    }
    let retry = fit(&[&samples], &terms, 0.02, budget()).unwrap();
    assert_eq!(baseline.strengths, retry.strengths);
    assert_eq!(baseline.predicted_stress_pa, retry.predicted_stress_pa);
}
