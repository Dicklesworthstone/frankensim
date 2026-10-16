use super::*;
use crate::{CorrelationModel, PropagationMethod, UqStatus};

fn plan() -> UqPlan {
    UqPlan::new("temperature", PropagationMethod::QuasiMonteCarlo, 32)
        .with_parameter(ParameterUncertainty::uniform("power", 4.0, 6.0, "W"))
        .with_parameter(ParameterUncertainty::uniform("inlet", 290.0, 310.0, "K"))
        .with_correlation(CorrelationModel::Independent)
        .with_compliance_threshold(305.0)
}
fn layout() -> QmcConfig { QmcConfig { replicates: 4, samples_per_replicate: 8 } }
fn matrix(rho: f64) -> Vec<Vec<f64>> { vec![vec![1.0, rho], vec![rho, 1.0]] }
fn model() -> ContentHash { fs_blake3::hash_domain("copula-qmc-test", b"power+inlet") }
fn value(x: &[f64]) -> Result<f64, &'static str> { Ok(x[0] + x[1]) }

#[test]
fn singular_dependent_inputs_keep_physical_support_and_units() {
    for rho in [-1.0, 1.0] {
        let mut run = GaussianCopulaQmcExecution::new(&plan(), &matrix(rho), layout()).unwrap();
        assert_eq!(run.marginals()[0].unit, "W");
        let report = run.advance(32, || false, |x| {
            assert!((4.0..=6.0).contains(&x[0]) && (290.0..=310.0).contains(&x[1]));
            let (a, b) = ((x[0] - 4.0) / 2.0, (x[1] - 290.0) / 20.0);
            assert!((if rho > 0.0 { a - b } else { a + b - 1.0 }).abs() < 4e-15);
            value(x)
        });
        assert_eq!(report.status, UqStatus::Complete);
        assert_eq!(report.completed_replicates, 4);
        assert!(report.estimate.unwrap().standard_error.is_some());
    }
}

#[test]
fn physical_qoi_statistics_are_not_statistics_of_latent_normals() {
    let mut p = plan();
    p.parameters[0] = ParameterUncertainty::uniform("power", -0.0, -0.0, "W");
    p.parameters[1] = ParameterUncertainty::uniform("inlet", 300.0, 300.0, "K");
    let mut run = GaussianCopulaQmcExecution::new(&p, &matrix(0.75), layout()).unwrap();
    let report = run.advance(32, || false, |x| {
        assert_eq!(x[0].to_bits(), (-0.0_f64).to_bits()); value(x)
    });
    assert_eq!(report.estimate.unwrap().mean, 300.0);
    assert_eq!(report.compliance.unwrap().mean, 1.0);
}

#[test]
fn interrupted_net_resumes_exact_physical_vector_and_replicate_statistics() {
    let p = plan(); let r = matrix(0.6);
    let mut full = GaussianCopulaQmcExecution::new(&p, &r, layout()).unwrap();
    let mut expected = Vec::new();
    full.advance(32, || false, |x| { expected.push(x.to_vec()); value(x) });
    let mut split = GaussianCopulaQmcExecution::new(&p, &r, layout()).unwrap();
    let mut actual = Vec::new();
    split.advance(11, || false, |x| { actual.push(x.to_vec()); value(x) });
    assert_eq!(split.report().completed_replicates, 1);
    assert!(split.report().estimate.unwrap().standard_error.is_none());
    split.advance_interruptible(1, || false, |x| {
        assert_eq!(x, expected[11].as_slice()); Ok::<_, &str>(None)
    });
    assert_eq!(split.evaluations_attempted(), 11);
    let bytes = split.checkpoint(model()).unwrap();
    let mut resumed = GaussianCopulaQmcExecution::restore(&p, &r, layout(), model(), &bytes).unwrap();
    resumed.advance(32, || false, |x| { actual.push(x.to_vec()); value(x) });
    assert_eq!(actual, expected);
    assert_eq!(resumed.observations(), full.observations());
    assert_eq!(resumed.report(), full.report());
}

#[test]
fn checkpoint_binds_physical_support_units_joint_law_and_layout() {
    let p = plan(); let r = matrix(0.5);
    let bytes = GaussianCopulaQmcExecution::new(&p, &r, layout()).unwrap().checkpoint(model()).unwrap();
    for change in 0..3 {
        let mut different = p.clone();
        match change {
            0 => different.parameters[0] = ParameterUncertainty::uniform("power", 3.0, 6.0, "W"),
            1 => different.parameters[0].unit = "kW".into(),
            _ => different.seed += 1,
        }
        assert!(GaussianCopulaQmcExecution::restore(&different, &r, layout(), model(), &bytes).is_err());
    }
    assert!(GaussianCopulaQmcExecution::restore(&p, &matrix(-0.5), layout(), model(), &bytes).is_err());
    assert!(GaussianCopulaQmcExecution::restore(&p, &r,
        QmcConfig { replicates: 2, samples_per_replicate: 16 }, model(), &bytes).is_err());
}

#[test]
fn numerical_or_physical_refusal_is_terminal_not_a_redraw() {
    for nonfinite in [false, true] {
        let mut run = GaussianCopulaQmcExecution::new(&plan(), &matrix(0.5), layout()).unwrap();
        run.advance(3, || false, value);
        let report = run.advance(32, || false, |_| {
            if nonfinite { Ok(f64::NAN) } else { Err("physical solve failed") }
        });
        assert_eq!(report.status, UqStatus::Refused);
        assert!(report.estimate.is_none() && report.compliance.is_none());
        assert_eq!(report.samples_evaluated, 4);
        run.advance(32, || false, |_| -> Result<f64, &str> { panic!("terminal execution") });
        assert!(run.checkpoint(model()).is_err());
    }
}

#[test]
fn admission_and_zero_budget_do_not_run_physics() {
    let mut p = plan();
    p.method = PropagationMethod::MonteCarlo;
    assert!(GaussianCopulaQmcExecution::new(&p, &matrix(0.0), layout()).is_err());
    p = plan(); p.correlation = CorrelationModel::Unknown;
    assert!(GaussianCopulaQmcExecution::new(&p, &matrix(0.0), layout()).is_err());
    assert!(GaussianCopulaQmcExecution::new(&plan(), &matrix(1.1), layout()).is_err());
    assert!(GaussianCopulaQmcExecution::new(&plan(), &matrix(f64::NAN), layout()).is_err());
    assert!(GaussianCopulaQmcExecution::new(&plan(), &[vec![1.0]], layout()).is_err());
    let mut run = GaussianCopulaQmcExecution::new(&plan(), &matrix(0.0), layout()).unwrap();
    assert_eq!(run.advance(0, || false,
        |_| -> Result<f64, &str> { panic!("no work") }).samples_evaluated, 0);
    assert_eq!(run.advance(1, || true,
        |_| -> Result<f64, &str> { panic!("cancelled") }).status, UqStatus::Cancelled);
}
