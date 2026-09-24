use fs_uq::{CorrelationModel, ParameterUncertainty, PropagationMethod, QmcConfig,
    QmcExecution, UqPlan, UqStatus};
use fs_rand::qmc::Sobol;

fn config() -> QmcConfig { QmcConfig { replicates: 3, samples_per_replicate: 128 } }
fn plan() -> UqPlan {
    let mut p = UqPlan::new("response", PropagationMethod::QuasiMonteCarlo, 384)
        .with_parameter(ParameterUncertainty::uniform("x", 0.0, 1.0, "1"))
        .with_correlation(CorrelationModel::Independent).with_compliance_threshold(0.5);
    p.seed = 83;
    p
}
fn value(x: &[f64]) -> Result<f64, &'static str> { Ok(x[0]) }

#[test]
fn exact_existing_net_points_include_the_origin_and_independent_replica_keys() {
    let p = plan();
    let mut run = QmcExecution::new(&p, config()).unwrap();
    let mut visited = Vec::new();
    let report = run.advance(usize::MAX, || false, |x| {
        visited.push(x[0]); value(x)
    });
    assert_eq!(report.status, UqStatus::Complete);
    assert_eq!(report.completed_replicates, 3);
    assert_eq!(report.samples_evaluated, 384);
    for r in 0..3 {
        let seed = p.seed.wrapping_add((r as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15));
        let sobol = Sobol::scrambled(1, seed);
        let mut bins = vec![false; 128];
        for i in 0..128 {
            let mut expected = [0.0]; sobol.point(i as u32, &mut expected);
            expected[0] += 1.0 / 8_589_934_592.0;
            assert_eq!(visited[r * 128 + i].to_bits(), expected[0].to_bits());
            assert!(expected[0] > 0.0 && expected[0] < 1.0);
            let bin = (expected[0] * 128.0) as usize;
            assert!(!bins[bin]); bins[bin] = true;
        }
        assert!(bins.into_iter().all(|b| b));
    }
    assert_ne!(&visited[..128], &visited[128..256]);
    let means: Vec<_> = visited.chunks_exact(128).map(|v| v.iter().sum::<f64>() / 128.0).collect();
    let mean = means.iter().sum::<f64>() / 3.0;
    let error = (means.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / 6.0).sqrt();
    let estimate = report.estimate.unwrap();
    assert!((estimate.mean - mean).abs() < 1e-14);
    assert!((estimate.standard_error.unwrap() - error).abs() < 1e-14);
    assert!((estimate.mean - 0.5).abs() < 1.0 / 128.0);
    // Every net has exactly half its points below the threshold. The within-net
    // Bernoulli variation is nonzero, but the between-replicate error is zero.
    let compliance = report.compliance.unwrap();
    assert_eq!(compliance.mean, 0.5);
    assert_eq!(compliance.standard_error, Some(0.0));
}

#[test]
fn incomplete_nets_are_paid_work_not_replica_estimates() {
    let mut run = QmcExecution::new(&plan(), config()).unwrap();
    let report = run.advance(127, || false, value);
    assert_eq!(report.samples_accepted, 127);
    assert!(report.estimate.is_none()); assert!(report.compliance.is_none());
    let report = run.advance(1, || false, value);
    assert_eq!(report.completed_replicates, 1);
    assert!(report.estimate.as_ref().unwrap().standard_error.is_none());
    let previous = report.estimate;
    let report = run.advance(71, || false, value);
    assert_eq!(report.estimate, previous);
    assert_eq!(report.samples_accepted, 199);
    assert_eq!(run.advance(0, || panic!("zero work polls nothing"), value), report);
    let report = run.advance(usize::MAX, || false, value);
    assert_eq!(report.completed_replicates, 3);
}

#[test]
fn split_clone_and_interruption_replay_the_exact_same_points_and_reports() {
    let p = plan();
    let mut full = QmcExecution::new(&p, config()).unwrap();
    full.advance(384, || false, value);
    let mut split = QmcExecution::new(&p, config()).unwrap();
    split.advance(133, || false, value);
    let work = split.evaluations_attempted();
    let mut interrupted_point = Vec::new();
    let report = split.advance_interruptible(50, || false, |x| {
        interrupted_point = x.to_vec(); Ok::<_, &str>(None)
    });
    assert_eq!(report.status, UqStatus::Cancelled);
    assert_eq!(split.evaluations_attempted(), work);
    let mut resumed = split.clone();
    resumed.advance(1, || false, |x| { assert_eq!(x, interrupted_point.as_slice()); value(x) });
    for chunk in [7, 113, 1, 900] { resumed.advance(chunk, || false, value); }
    assert_eq!(resumed.observations(), full.observations());
    assert_eq!(resumed.report(), full.report());
    let after = resumed.advance(100, || false, |_| -> Result<f64, &str> { panic!("complete is terminal") });
    assert_eq!(after, full.report());
}

#[test]
fn rejected_samples_hide_all_estimates_and_remain_terminal() {
    for result in [Err("physical model refused"), Ok(f64::NAN), Ok(f64::INFINITY)] {
        let mut run = QmcExecution::new(&plan(), config()).unwrap();
        run.advance(256, || false, value);
        let report = run.advance(10, || false, |_| result);
        assert_eq!(report.status, UqStatus::Refused);
        assert_eq!(report.samples_evaluated, 257);
        assert_eq!(report.samples_accepted, 256);
        assert!(report.estimate.is_none()); assert!(report.compliance.is_none());
        assert!(report.replicate_means.is_empty());
        assert!(report.rejection_reason.is_some());
        assert_eq!(run.advance(8, || false, |_| -> Result<f64, &str> { panic!("refused is terminal") }), report);
    }
}

#[test]
fn cancellation_before_work_and_finite_extreme_constant_values() {
    let mut run = QmcExecution::new(&plan(), config()).unwrap();
    let report = run.advance(100, || true, |_| -> Result<f64, &str> { panic!("pre-cancel") });
    assert_eq!(report.status, UqStatus::Cancelled);
    assert_eq!(report.samples_evaluated, 0);
    let report = run.advance(384, || false, |_| Ok::<_, &str>(f64::MAX));
    assert_eq!(report.status, UqStatus::Complete);
    let estimate = report.estimate.unwrap();
    assert_eq!(estimate.mean, f64::MAX); assert_eq!(estimate.standard_error, Some(0.0));
}

#[test]
fn gaussian_marginals_and_singular_joint_gaussians_use_the_existing_transform() {
    let config = QmcConfig { replicates: 4, samples_per_replicate: 1024 };
    let p = UqPlan::new("joint", PropagationMethod::QuasiMonteCarlo, 4096)
        .with_parameter(ParameterUncertainty::gaussian("a", 2.0, 3.0, "m"))
        .with_parameter(ParameterUncertainty::gaussian("b", -1.0, 2.0, "m"))
        .with_correlation(CorrelationModel::JointGaussian { matrix: vec![vec![1.0, 1.0], vec![1.0, 1.0]] });
    let mut run = QmcExecution::new(&p, config).unwrap();
    let report = run.advance(4096, || false, |x| {
        let a = (x[0] - 2.0) / 3.0;
        let b = (x[1] + 1.0) / 2.0;
        assert!((a - b).abs() < 1e-13);
        Ok::<_, &str>(a * b)
    });
    assert_eq!(report.status, UqStatus::Complete);
    assert!((report.estimate.unwrap().mean - 1.0).abs() < 0.04);
}

#[test]
fn every_layout_and_joint_measure_is_admitted_before_evaluation() {
    for c in [QmcConfig { replicates: 1, samples_per_replicate: 384 },
              QmcConfig { replicates: 3, samples_per_replicate: 127 },
              QmcConfig { replicates: 3, samples_per_replicate: 64 },
              QmcConfig { replicates: usize::MAX, samples_per_replicate: 128 }] {
        assert!(QmcExecution::new(&plan(), c).is_err());
    }
    let mut p = plan(); p.method = PropagationMethod::MonteCarlo;
    assert!(QmcExecution::new(&p, config()).is_err());
    let mut p = plan(); p.parameters[0].unit.clear();
    assert!(QmcExecution::new(&p, config()).is_err());
    let mut p = plan(); p.parameters[0] = ParameterUncertainty::interval("x", 0.0, 1.0, "m");
    assert!(QmcExecution::new(&p, config()).is_err());
    let mut p = plan();
    p.parameters.push(ParameterUncertainty::uniform("y", 0.0, 1.0, "1"));
    p.correlation = CorrelationModel::Unknown;
    assert!(QmcExecution::new(&p, config()).is_err());
    p.correlation = CorrelationModel::JointGaussian { matrix: vec![vec![1.0, 0.0], vec![0.0, 1.0]] };
    assert!(QmcExecution::new(&p, config()).is_err());
    let mut p = plan();
    for i in 1..11 { p.parameters.push(ParameterUncertainty::uniform(format!("x{i}"), 0.0, 1.0, "1")); }
    assert!(QmcExecution::new(&p, config()).is_err());
}

#[test]
fn durable_recovery_retains_an_incomplete_net_and_retries_the_interrupted_point() {
    let p = plan();
    let model = fs_blake3::hash_domain("qmc-test-model", b"response = x");
    let mut full = QmcExecution::new(&p, config()).unwrap();
    full.advance(384, || false, value);
    let mut split = QmcExecution::new(&p, config()).unwrap();
    split.advance(133, || false, value);
    let mut interrupted = Vec::new();
    split.advance_interruptible(1, || false, |x| {
        interrupted = x.to_vec(); Ok::<_, &str>(None)
    });
    let bytes = split.checkpoint(model).unwrap();
    let mut resumed = QmcExecution::restore(&p, config(), model, &bytes).unwrap();
    assert_eq!(resumed.report(), split.report());
    assert_eq!(resumed.observations().len(), 133);
    assert_eq!(resumed.report().completed_replicates, 1);
    resumed.advance(1, || false, |x| { assert_eq!(x, interrupted.as_slice()); value(x) });
    resumed.advance(384, || false, value);
    assert_eq!(resumed.report(), full.report());
    assert_eq!(resumed.checkpoint(model).unwrap(), full.checkpoint(model).unwrap());
    let mut terminal = QmcExecution::restore(&p, config(), model, &resumed.checkpoint(model).unwrap()).unwrap();
    assert_eq!(terminal.advance(1, || panic!("complete"), |_| -> Result<f64, &str> { panic!("complete") }), full.report());
}

#[test]
fn durable_qmc_identity_binds_layout_plan_model_and_sampler() {
    let p = plan();
    let model = fs_blake3::hash_domain("qmc-test-model", b"response = x");
    let mut execution = QmcExecution::new(&p, config()).unwrap();
    execution.advance(2, || false, value);
    let bytes = execution.checkpoint(model).unwrap();
    // Same total sample count with different net grouping changes EVERY point
    // after the first net and the units of the error estimate; reject it.
    let alternate = QmcConfig { replicates: 6, samples_per_replicate: 64 };
    assert_eq!(QmcExecution::restore(&p, alternate, model, &bytes).unwrap_err(), fs_uq::UqCheckpointError::IdentityMismatch);
    let mut changed = p.clone(); changed.seed += 1;
    assert_eq!(QmcExecution::restore(&changed, config(), model, &bytes).unwrap_err(), fs_uq::UqCheckpointError::IdentityMismatch);
    assert_eq!(QmcExecution::restore(&p, config(), fs_blake3::ContentHash([0; 32]), &bytes).unwrap_err(), fs_uq::UqCheckpointError::IdentityMismatch);
    let mut mc_plan = p.clone(); mc_plan.method = PropagationMethod::MonteCarlo;
    let mc = fs_uq::UqExecution::new(&mc_plan).unwrap().checkpoint(model).unwrap();
    assert!(QmcExecution::restore(&p, config(), model, &mc).is_err());
    assert!(fs_uq::UqExecution::restore(&mc_plan, model, &bytes).is_err());
    execution.advance(1, || false, |_| Err::<f64, _>("solver refused"));
    assert_eq!(execution.checkpoint(model), Err(fs_uq::UqCheckpointError::NotResumable));
}

#[test]
fn durable_qmc_rejects_corruption_and_preserves_exact_finite_bits() {
    let p = plan();
    let model = fs_blake3::hash_domain("qmc-test-model", b"response = x");
    let mut execution = QmcExecution::new(&p, config()).unwrap();
    execution.advance(1, || false, |_| Ok::<_, &str>(-0.0));
    let bytes = execution.checkpoint(model).unwrap();
    let restored = QmcExecution::restore(&p, config(), model, &bytes).unwrap();
    assert_eq!(restored.observations()[0].to_bits(), (-0.0_f64).to_bits());
    for index in [0, 8, 40, 41, 49, bytes.len() - 1] {
        let mut changed = bytes.clone(); changed[index] ^= 0x80;
        assert!(QmcExecution::restore(&p, config(), model, &changed).is_err());
    }
    let mut nonfinite = bytes;
    nonfinite[49..57].copy_from_slice(&f64::NAN.to_bits().to_le_bytes());
    let checksum_start = nonfinite.len() - 32;
    let checksum = fs_blake3::hash_domain("org.frankensim.uq.qmc.checkpoint.v1", &nonfinite[..checksum_start]);
    nonfinite[checksum_start..].copy_from_slice(&checksum.0);
    assert!(QmcExecution::restore(&p, config(), model, &nonfinite).is_err());
}
