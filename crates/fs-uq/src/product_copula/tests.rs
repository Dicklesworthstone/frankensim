use super::*;
use crate::{PropagationMethod, UqStatus};

fn plan(n: usize) -> UqPlan {
    UqPlan::new("temperature", PropagationMethod::MonteCarlo, n)
        .with_correlation(CorrelationModel::Independent)
        .with_parameter(ParameterUncertainty::uniform("power", 0.0, 1.0, "W"))
        .with_parameter(ParameterUncertainty::uniform("inlet", 0.0, 1.0, "K"))
        .with_compliance_threshold(1.0)
}
fn matrix(rho: f64) -> Vec<Vec<f64>> { vec![vec![1.0, rho], vec![rho, 1.0]] }
fn model() -> ContentHash { fs_blake3::hash_domain("test-copula-model", b"x0+x1") }
fn sum(x: &[f64]) -> Result<f64, &'static str> { Ok(x[0] + x[1]) }

#[test]
fn normal_cdf_matches_independent_reference_and_preserves_support() {
    // Standard normal CDF values from independent high-precision evaluation.
    for (z, expected) in [(-8.0, 6.220960574271784e-16), (-3.0, 0.0013498980316300945),
        (-1.0, 0.15865525393145705), (0.0, 0.5), (0.25, 0.5987063256829237),
        (1.0, 0.8413447460685429), (3.0, 0.9986501019683699), (8.0, 0.9999999999999993)] {
        let value = physical_parameters(&plan(2).parameters[..1], &[z]).unwrap()[0];
        assert!((value - expected).abs() <= 2e-14 * expected.max(1e-15), "z={z}, value={value}");
    }
    let parameters = [ParameterUncertainty::uniform("wide", -f64::MAX, f64::MAX, "W"),
        ParameterUncertainty::uniform("constant", -0.0, 0.0, "K")];
    for z in [-f64::MAX, -100.0, -8.0, 0.0, 8.0, 100.0, f64::MAX] {
        let values = physical_parameters(&parameters, &[z, z]).unwrap();
        assert!(values[0].is_finite());
        assert_eq!(values[1].to_bits(), (-0.0_f64).to_bits());
    }
    for latent in [vec![0.0], vec![f64::NAN, 0.0], vec![0.0, f64::INFINITY]] {
        assert!(physical_parameters(&parameters, &latent).is_err());
    }
}

#[test]
fn perfect_latent_dependence_is_not_replaced_by_independent_draws() {
    for rho in [-1.0, 1.0] {
        let mut execution = GaussianCopulaExecution::new(&plan(64), &matrix(rho)).unwrap();
        let report = execution.advance(64, || false, |x| {
            if rho == 1.0 { assert_eq!(x[0].to_bits(), x[1].to_bits()); }
            else { assert!((x[0] + x[1] - 1.0).abs() < 3e-16); }
            sum(x)
        });
        assert_eq!(report.status, UqStatus::Complete);
    }
}

#[test]
fn marginal_moments_and_physical_uniform_correlation_follow_the_declared_copula() {
    let count = 8192;
    let mut p = plan(count); p.seed = 987654;
    let mut execution = GaussianCopulaExecution::new(&p, &matrix(0.8)).unwrap();
    let mut sums = [0.0; 5];
    let report = execution.advance(count, || false, |x| {
        for (total, value) in sums.iter_mut().zip([x[0], x[1], x[0]*x[0], x[1]*x[1], x[0]*x[1]]) {
            *total += value;
        }
        sum(x)
    });
    assert_eq!(report.status, UqStatus::Complete);
    for value in &mut sums { *value /= count as f64; }
    let a = sums[2] - sums[0]*sums[0]; let b = sums[3] - sums[1]*sums[1];
    assert!((sums[0]-0.5).abs() < 0.02 && (sums[1]-0.5).abs() < 0.02);
    assert!((a-1.0/12.0).abs() < 0.005 && (b-1.0/12.0).abs() < 0.005);
    let correlation = (sums[4]-sums[0]*sums[1]) / (a*b).sqrt();
    // 6/pi * asin(0.8/2), not the latent Pearson coefficient 0.8.
    assert!((correlation - 0.7859392826067277).abs() < 0.025, "{correlation}");
}

#[test]
fn malformed_joint_laws_refuse_before_any_model_work() {
    for m in [vec![], vec![vec![1.0]], vec![vec![1.0, 0.5], vec![0.4, 1.0]],
        matrix(f64::NAN), matrix(1.01), vec![vec![0.9, 0.0], vec![0.0, 1.0]]] {
        assert!(GaussianCopulaExecution::new(&plan(8), &m).is_err());
    }
    let mut p = plan(8);
    p.parameters.push(ParameterUncertainty::uniform("third", 0.0, 1.0, "1"));
    let non_psd = vec![vec![1.0, 0.9, 0.9], vec![0.9, 1.0, -0.9], vec![0.9, -0.9, 1.0]];
    assert!(GaussianCopulaExecution::new(&p, &non_psd).is_err());
    for kind in [UncertaintyKind::AleatoryGaussian { mean: 0.0, std_dev: 1.0 },
        UncertaintyKind::AleatoryUniform { lo: 2.0, hi: 1.0 },
        UncertaintyKind::EpistemicInterval { lo: 0.0, hi: 1.0 }] {
        let mut p = plan(8); p.parameters[0].kind = kind;
        assert!(GaussianCopulaExecution::new(&p, &matrix(0.0)).is_err());
    }
    let mut p = plan(8); p.correlation = CorrelationModel::Unknown;
    assert!(GaussianCopulaExecution::new(&p, &matrix(0.0)).is_err());
    let mut p = plan(8); p.method = PropagationMethod::QuasiMonteCarlo;
    assert!(GaussianCopulaExecution::new(&p, &matrix(0.0)).is_err());
}

#[test]
fn all_chunk_boundaries_restore_exact_physical_draws_and_results() {
    let p = plan(19); let m = matrix(0.7);
    let mut full = GaussianCopulaExecution::new(&p, &m).unwrap();
    let mut trace = Vec::new();
    let expected = full.advance(19, || false, |x| { trace.push(x.to_vec()); sum(x) });
    for split in 0..=19 {
        let mut actual = GaussianCopulaExecution::new(&p, &m).unwrap();
        let mut observed = Vec::new();
        actual.advance(split, || false, |x| { observed.push(x.to_vec()); sum(x) });
        let bytes = actual.checkpoint(model()).unwrap();
        let mut restored = GaussianCopulaExecution::restore(&p, &m, model(), &bytes).unwrap();
        let result = restored.advance(19, || false, |x| { observed.push(x.to_vec()); sum(x) });
        assert_eq!(result, expected);
        assert_eq!(observed, trace);
        assert_eq!(restored.checkpoint(model()).unwrap(), full.checkpoint(model()).unwrap());
    }
}

#[test]
fn physical_supports_units_matrix_and_model_are_bound_even_with_unchanged_latent_marginals() {
    let p = plan(8); let m = matrix(0.2);
    let mut execution = GaussianCopulaExecution::new(&p, &m).unwrap();
    execution.advance(2, || false, sum);
    let bytes = execution.checkpoint(model()).unwrap();
    let mut changed = p.clone(); changed.parameters[0].kind = UncertaintyKind::AleatoryUniform { lo: 0.0, hi: 2.0 };
    assert!(matches!(GaussianCopulaExecution::restore(&changed, &m, model(), &bytes), Err(UqCheckpointError::IdentityMismatch)));
    let mut changed = p.clone(); changed.parameters[0].unit = "kW".into();
    assert!(matches!(GaussianCopulaExecution::restore(&changed, &m, model(), &bytes), Err(UqCheckpointError::IdentityMismatch)));
    assert!(GaussianCopulaExecution::restore(&p, &matrix(0.3), model(), &bytes).is_err());
    let other = fs_blake3::hash_domain("other", b"model");
    assert!(GaussianCopulaExecution::restore(&p, &m, other, &bytes).is_err());
    for length in [0, 8, bytes.len()-1] { assert!(GaussianCopulaExecution::restore(&p, &m, model(), &bytes[..length]).is_err()); }
    let mut corrupt = bytes; let last = corrupt.len()-1; corrupt[last] ^= 1;
    assert!(GaussianCopulaExecution::restore(&p, &m, model(), &corrupt).is_err());
}

#[test]
fn interrupted_draw_retries_but_rejected_draw_remains_terminal() {
    let p = plan(8); let m = matrix(-0.4);
    let mut execution = GaussianCopulaExecution::new(&p, &m).unwrap();
    execution.advance(2, || false, sum);
    let mut interrupted = Vec::new();
    let report = execution.advance_interruptible(1, || false, |x| {
        interrupted = x.to_vec(); Ok::<_, &str>(None)
    });
    assert_eq!(report.status, UqStatus::Cancelled);
    let bytes = execution.checkpoint(model()).unwrap();
    let mut restored = GaussianCopulaExecution::restore(&p, &m, model(), &bytes).unwrap();
    restored.advance(1, || false, |x| { assert_eq!(x, interrupted.as_slice()); sum(x) });
    let report = restored.advance(1, || false, |_| Err::<f64, _>("physical refusal"));
    assert_eq!(report.status, UqStatus::Refused);
    assert!(report.rejection_reason.unwrap().contains("physical refusal"));
    assert!(matches!(restored.checkpoint(model()), Err(UqCheckpointError::NotResumable)));
    restored.advance(10, || false, |_| -> Result<f64, &str> { panic!("terminal refusal must not resume") });
}

#[test]
fn zero_budget_and_precancelled_execution_do_not_call_the_model() {
    let mut execution = GaussianCopulaExecution::new(&plan(8), &matrix(0.5)).unwrap();
    let report = execution.advance(0, || false, |_| -> Result<f64, &str> { panic!("zero allowance") });
    assert_eq!(report.samples_evaluated, 0);
    let report = execution.advance(8, || true, |_| -> Result<f64, &str> { panic!("cancelled") });
    assert_eq!(report.status, UqStatus::Cancelled);
    assert!(execution.monte_carlo().observations().is_empty());
}
