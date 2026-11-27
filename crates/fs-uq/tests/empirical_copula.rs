//! End-to-end fixtures for empirical physical coordinates in the two copula executors.
use fs_blake3::{hash_bytes, hash_domain};
use fs_uq::{
    CorrelationModel, EmpiricalMarginal, GaussianCopulaExecution, GaussianCopulaQmcExecution,
    INVERSE_ECDF_RULE_ID, ParameterUncertainty, PropagationMethod, QmcConfig, QmcExecution,
    UqCheckpointError, UqControlError, UqExecution, UqPlan, UqStatus,
};

const METHODS: [PropagationMethod; 2] = [
    PropagationMethod::MonteCarlo,
    PropagationMethod::QuasiMonteCarlo,
];
const NET: QmcConfig = QmcConfig {
    replicates: 4,
    samples_per_replicate: 8,
};

fn plan(method: PropagationMethod, count: usize) -> UqPlan {
    UqPlan::new("load", method, count)
        .with_correlation(CorrelationModel::Independent)
        .with_compliance_threshold(10.0)
}

fn empirical(name: &str, values: Vec<f64>) -> ParameterUncertainty {
    let marginal = EmpiricalMarginal::new(values, INVERSE_ECDF_RULE_ID).unwrap();
    ParameterUncertainty::empirical(name, marginal, "kg")
}

fn matrix(rho: f64) -> Vec<Vec<f64>> {
    vec![vec![1.0, rho], vec![rho, 1.0]]
}

fn observe(plan: &UqPlan, matrix: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let mut rows = Vec::with_capacity(plan.budget_max_samples);
    let collect = |x: &[f64]| {
        rows.push(x.to_vec());
        Ok::<_, &str>(x[0])
    };
    let status = match plan.method {
        PropagationMethod::MonteCarlo => {
            GaussianCopulaExecution::new(plan, matrix)
                .unwrap()
                .advance(plan.budget_max_samples, || false, collect)
                .status
        }
        PropagationMethod::QuasiMonteCarlo => {
            GaussianCopulaQmcExecution::new(plan, matrix, NET)
                .unwrap()
                .advance(plan.budget_max_samples, || false, collect)
                .status
        }
        _ => unreachable!(),
    };
    assert_eq!(status, UqStatus::Complete);
    assert_eq!(rows.len(), plan.budget_max_samples);
    rows
}

#[test]
fn mc_marginal_cdfs_match_the_observed_masses_including_ties() {
    let source = [
        vec![-7.0, 0.0, 0.0, 0.0, 19.0],
        vec![1.0, 3.0, 5.0, 8.0, 8.0, 10.0],
    ];
    let count = 8192;
    let p = plan(PropagationMethod::MonteCarlo, count)
        .with_parameter(empirical("a", source[0].clone()))
        .with_parameter(empirical("b", source[1].clone()));
    let rows = observe(&p, &matrix(0.0));
    // Two-sided DKW at alpha=0.0005 per marginal (joint level >=99.9% by union bound).
    let epsilon = ((2.0_f64 / 0.0005).ln() / (2.0 * count as f64)).sqrt();
    for (j, sample) in source.iter().enumerate() {
        assert!(rows.iter().all(|row| sample.contains(&row[j])));
        for &x in sample {
            let expected = sample.partition_point(|&v| v <= x) as f64 / sample.len() as f64;
            let actual = rows.iter().filter(|row| row[j] <= x).count() as f64 / count as f64;
            assert!(
                (actual - expected).abs() <= epsilon,
                "coordinate {j} at {x}: {actual} vs {expected}"
            );
        }
    }
}

// Centered average ranks, including ties caused by repeated draws of a single order statistic.
fn centered_ranks(rows: &[Vec<f64>], j: usize) -> Vec<f64> {
    let mut order: Vec<usize> = (0..rows.len()).collect();
    order.sort_unstable_by(|&a, &b| rows[a][j].total_cmp(&rows[b][j]));
    let mut ranks = vec![0.0; rows.len()];
    let mut start = 0;
    while start < order.len() {
        let value = rows[order[start]][j];
        let end = start + order[start..].partition_point(|&i| rows[i][j] <= value);
        let rank = (start + end) as f64 * 0.5 - rows.len() as f64 * 0.5;
        for &i in &order[start..end] {
            ranks[i] = rank;
        }
        start = end;
    }
    ranks
}

#[test]
fn empirical_ranks_follow_the_declared_gaussian_copula() {
    let count = 8192;
    let rho = 0.7;
    let p = plan(PropagationMethod::MonteCarlo, count)
        .with_parameter(empirical("mass", (1..=8192).map(f64::from).collect()))
        .with_parameter(empirical(
            "square",
            (1..=8192).map(|x| f64::from(x).powi(2)).collect(),
        ));
    let rows = observe(&p, &matrix(rho));
    let a = centered_ranks(&rows, 0);
    let b = centered_ranks(&rows, 1);
    let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
    let measured = dot(&a, &b) / (dot(&a, &a) * dot(&b, &b)).sqrt();
    let target = 6.0 / core::f64::consts::PI * (rho / 2.0).asin();
    assert!(
        (measured - target).abs() < 0.03,
        "Spearman {measured}, expected {target}"
    );
}

#[test]
fn mixed_coordinates_and_matrix_permute_together() {
    // A full-rank matrix with distinct correlations; all three coordinate positions move.
    let matrix = [
        vec![1.0, 0.65, -0.4],
        vec![0.65, 1.0, 0.1],
        vec![-0.4, 0.1, 1.0],
    ];
    let rotated_matrix = [
        vec![1.0, -0.4, 0.1],
        vec![-0.4, 1.0, 0.65],
        vec![0.1, 0.65, 1.0],
    ];
    for method in METHODS {
        let mut p = plan(method, 32)
            .with_parameter(ParameterUncertainty::uniform("u", 0.0, 1.0, "1"))
            .with_parameter(empirical("short", vec![-1.0, 2.0]))
            .with_parameter(empirical("long", vec![5.0, 8.0, 9.0, 15.0]));
        let original = observe(&p, &matrix);
        p.parameters.rotate_right(1);
        for (row, rotated) in original.iter().zip(observe(&p, &rotated_matrix)) {
            assert!((0.0..=1.0).contains(&row[0]));
            assert!([-1.0, 2.0].contains(&row[1]));
            assert!([5.0, 8.0, 9.0, 15.0].contains(&row[2]));
            assert_eq!(rotated, [row[2], row[0], row[1]]);
        }
    }
}

#[test]
fn both_executors_report_empirical_sample_counts_and_refuse_mean_controls() {
    for method in METHODS {
        let p = plan(method, 32)
            .with_parameter(ParameterUncertainty::uniform("u", 0.0, 1.0, "1"))
            .with_parameter(empirical("lead", vec![2.0, 2.0, 7.0]));
        let m = matrix(0.4);
        let summarize = |(p, law): (&ParameterUncertainty, &EmpiricalMarginal)| {
            (p.name.clone(), law.sample_count(), law.observed_support())
        };
        let (metadata, refusal) = if method == PropagationMethod::MonteCarlo {
            let run = GaussianCopulaExecution::new(&p, &m).unwrap();
            (
                run.empirical_marginals().map(summarize).collect::<Vec<_>>(),
                run.freeze_linear_control_variate(&[0.0; 2]).unwrap_err(),
            )
        } else {
            let run = GaussianCopulaQmcExecution::new(&p, &m, NET).unwrap();
            (
                run.empirical_marginals().map(summarize).collect::<Vec<_>>(),
                run.freeze_linear_control_variate(&[0.0; 2]).unwrap_err(),
            )
        };
        assert_eq!(metadata, [("lead".to_owned(), 3, [2.0, 7.0])]);
        assert_eq!(refusal, UqControlError::EmpiricalMarginalUnsupported);
    }
}

#[test]
fn both_checkpoints_bind_empirical_values_length_units_and_order() {
    let m = matrix(0.25);
    let model = hash_domain("empirical-copula-test", b"load");
    for method in METHODS {
        let p = plan(method, 32)
            .with_parameter(empirical("a", vec![1.0, 5.0, 6.0]))
            .with_parameter(empirical("b", vec![2.0, 7.0]));
        let check = |plan: &UqPlan| -> Result<Vec<u8>, UqCheckpointError> {
            if method == PropagationMethod::MonteCarlo {
                GaussianCopulaExecution::new(plan, &m)
                    .unwrap()
                    .checkpoint(model)
            } else {
                GaussianCopulaQmcExecution::new(plan, &m, NET)
                    .unwrap()
                    .checkpoint(model)
            }
        };
        let bytes = check(&p).unwrap();
        let restore = |plan: &UqPlan| {
            if method == PropagationMethod::MonteCarlo {
                GaussianCopulaExecution::restore(plan, &m, model, &bytes).map(|_| ())
            } else {
                GaussianCopulaQmcExecution::restore(plan, &m, NET, model, &bytes).map(|_| ())
            }
        };
        restore(&p).unwrap();
        for change in 0..5 {
            let mut changed = p.clone();
            match change {
                0 => changed.parameters[0] = empirical("a", vec![1.0, 4.0, 6.0]),
                1 => changed.parameters[0] = empirical("a", vec![1.0, 5.0, 5.0, 6.0]),
                2 => changed.parameters[0].unit = "g".into(),
                3 => {
                    // Swap just the laws: the latent names and matrix remain identical.
                    let (a, b) = changed.parameters.split_at_mut(1);
                    core::mem::swap(&mut a[0].kind, &mut b[0].kind);
                }
                _ => {
                    // The symmetric matrix and canonical latent names also survive this change.
                    let (a, b) = changed.parameters.split_at_mut(1);
                    core::mem::swap(&mut a[0].name, &mut b[0].name);
                }
            }
            assert_ne!(check(&changed).unwrap(), bytes);
            assert!(matches!(
                restore(&changed),
                Err(UqCheckpointError::IdentityMismatch)
            ));
        }
        let mut missing_unit = p;
        for unit in ["", " \t"] {
            missing_unit.parameters[0].unit = unit.to_owned();
            assert!(if method == PropagationMethod::MonteCarlo {
                GaussianCopulaExecution::new(&missing_unit, &m).is_err()
            } else {
                GaussianCopulaQmcExecution::new(&missing_unit, &m, NET).is_err()
            });
        }
    }
}

#[test]
fn interrupted_mixed_run_resumes_the_exact_physical_suffix() {
    let matrix = matrix(-0.6);
    let model = hash_domain("empirical-copula-test", b"sum");
    let sum = |x: &[f64]| Ok::<_, &str>(x.iter().sum::<f64>());
    for method in METHODS {
        let p = plan(method, 32)
            .with_parameter(empirical("measured", vec![1.0, 3.0, 3.0, 8.0]))
            .with_parameter(ParameterUncertainty::uniform("factor", -3.0, 9.0, "1"));
        // The same assertions exercise two concrete types with distinct report types.
        macro_rules! check_replay {
            ($new:expr, $restore:expr) => {{
                let mut full = $new.unwrap();
                let mut expected = Vec::new();
                let report = full.advance(
                    32,
                    || false,
                    |x| {
                        expected.push(x.to_vec());
                        sum(x)
                    },
                );
                let mut partial = $new.unwrap();
                let mut observed = Vec::new();
                partial.advance(
                    13,
                    || false,
                    |x| {
                        observed.push(x.to_vec());
                        sum(x)
                    },
                );
                let interrupted = partial.advance_interruptible(
                    1,
                    || false,
                    |x| {
                        assert_eq!(x, expected[13]);
                        Ok::<_, &str>(None)
                    },
                );
                assert_eq!(interrupted.status, UqStatus::Cancelled);
                let bytes = partial.checkpoint(model).unwrap();
                let mut resumed = ($restore)(&bytes).unwrap();
                let actual = resumed.advance(
                    32,
                    || false,
                    |x| {
                        observed.push(x.to_vec());
                        sum(x)
                    },
                );
                assert_eq!(actual, report);
                assert_eq!(observed, expected);
                assert_eq!(resumed.checkpoint(model), full.checkpoint(model));
            }};
        }
        if method == PropagationMethod::MonteCarlo {
            check_replay!(
                GaussianCopulaExecution::new(&p, &matrix),
                |bytes: &[u8]| { GaussianCopulaExecution::restore(&p, &matrix, model, bytes) }
            );
        } else {
            check_replay!(
                GaussianCopulaQmcExecution::new(&p, &matrix, NET),
                |bytes: &[u8]| {
                    GaussianCopulaQmcExecution::restore(&p, &matrix, NET, model, bytes)
                }
            );
        }
    }
}

#[test]
fn ordinary_drivers_still_refuse_empirical_plans() {
    for method in METHODS {
        let mut p = plan(method, 32)
            .with_parameter(empirical("a", vec![1.0, 2.0]))
            .with_parameter(empirical("b", vec![3.0, 5.0]));
        for correlation in [
            CorrelationModel::Independent,
            CorrelationModel::Correlated {
                matrix: matrix(0.5),
            },
            CorrelationModel::JointGaussian {
                matrix: matrix(0.5),
            },
        ] {
            p.correlation = correlation;
            if method == PropagationMethod::MonteCarlo {
                assert!(UqExecution::new(&p).is_err());
            } else {
                assert!(QmcExecution::new(&p, NET).is_err());
            }
        }
    }
}

#[test]
fn uniform_only_reports_and_checkpoints_match_pre_empirical_goldens() {
    // Captured at 151162deb90c658291a67a8254ad9dc1a77a9da7. Prefix and final
    // checkpoint bytes cover the model binding and each individual observation bit.
    let m = matrix(-0.45);
    let model = hash_domain("empirical-copula-test", b"sum");
    let sum = |x: &[f64]| Ok::<_, &str>(x.iter().sum::<f64>());
    for method in METHODS {
        let p = plan(method, 32)
            .with_parameter(ParameterUncertainty::uniform("mass", -5.0, 3.0, "kg"))
            .with_parameter(ParameterUncertainty::uniform(
                "temperature",
                0.25,
                110.75,
                "K",
            ));
        let (prefix, complete, report, expected) = if method == PropagationMethod::MonteCarlo {
            let mut run = GaussianCopulaExecution::new(&p, &m).unwrap();
            run.advance(13, || false, sum);
            let prefix = run.checkpoint(model).unwrap();
            let mut restored = GaussianCopulaExecution::restore(&p, &m, model, &prefix).unwrap();
            let report = restored.advance(32, || false, sum);
            (
                prefix,
                restored.checkpoint(model).unwrap(),
                format!("{report:?}"),
                [
                    "2cfda4f2401cf101ae3e06e6f1cf8898515c7e10edaec51c0796fed2e716e7fe",
                    "3e3bdb1a1be29b756adf9848d729607dc8f1c0eb2ceb65d76b627f634395e838",
                    "5587c8e833d54f3a9863b125969200bceaf02fab27aa91c0af010918597dedd1",
                ],
            )
        } else {
            let mut run = GaussianCopulaQmcExecution::new(&p, &m, NET).unwrap();
            run.advance(13, || false, sum);
            let prefix = run.checkpoint(model).unwrap();
            let mut restored =
                GaussianCopulaQmcExecution::restore(&p, &m, NET, model, &prefix).unwrap();
            let report = restored.advance(32, || false, sum);
            (
                prefix,
                restored.checkpoint(model).unwrap(),
                format!("{report:?}"),
                [
                    "f3b0fed037fec5be0a8f074bca9aa2e9df132e7a5abb14ee2eb7f3560c4ba161",
                    "c679e875a55bda79ee8724382dba5cd3e5a9e3f45acbb20ed4c67cc69be545b7",
                    "4ed18e501273a44bd0d233923e3a34262fc9d427b6c95e9680287c6349c55345",
                ],
            )
        };
        let actual = [prefix.as_slice(), complete.as_slice(), report.as_bytes()]
            .map(|bytes| hash_bytes(bytes).to_hex());
        assert_eq!(actual, expected);
    }
}
