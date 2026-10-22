use super::*;
use crate::{CorrelationModel, LhsPlan, ParameterUncertainty};

fn plan() -> LhsPlan {
    LhsPlan {
        target_qoi: "temperature".into(),
        compliance_threshold: Some(28.0),
        parameters: vec![
            ParameterUncertainty::uniform("power", -2.0, 4.0, "W"),
            ParameterUncertainty::gaussian("ambient", 3.0, 0.2, "K"),
        ],
        correlation: CorrelationModel::Independent,
        replicates: 3,
        samples_per_replicate: 5,
        seed: 59,
    }
}

fn affine(x: &[f64]) -> Result<f64, &'static str> {
    Ok(17.0 + 2.0 * x[0] + 3.0 * x[1])
}

#[test]
fn g1_affine_mean_control_preserves_physics_compliance_and_checkpoint_bytes() {
    let mut run = LhsExecution::new(&plan()).unwrap();
    let control = run.freeze_linear_control_variate(&[2.0, 3.0]).unwrap();
    assert_eq!(control.parameter_means(), [1.0, 3.0]);
    assert_eq!(control.gradient(), [2.0, 3.0]);
    let report = run.advance(15, || false, affine);
    let identity = fs_blake3::hash_domain("lhs-control-test", b"affine");
    let checkpoint = run.checkpoint(identity).unwrap();
    let result = run.assess_linear_control_variate(&control).unwrap();
    assert_eq!(result.raw, report.estimate);
    assert_eq!(result.raw_replicate_means, report.replicate_means);
    assert_eq!(result.samples_in_estimate, 15);
    let controlled = result.controlled.unwrap();
    assert!((controlled.mean - 28.0).abs() < 2e-14);
    assert!(controlled.standard_error.unwrap() < 2e-14);
    assert!(result.variance_ratio.unwrap() < 1e-20);
    assert_eq!(run.report(), report);
    assert_eq!(run.checkpoint(identity).unwrap(), checkpoint);
}

#[test]
fn g4_g5_partial_replicates_resume_without_invented_observations_or_recalibration() {
    let p = plan();
    let mut whole = LhsExecution::new(&p).unwrap();
    let control = whole.freeze_linear_control_variate(&[2.0, 3.0]).unwrap();
    whole.advance(15, || false, affine);
    let expected = whole.assess_linear_control_variate(&control).unwrap();
    let mut split = LhsExecution::new(&p).unwrap();
    split.advance(4, || false, affine);
    let partial = split.assess_linear_control_variate(&control).unwrap();
    assert_eq!(partial.samples_accepted, 4);
    assert_eq!(partial.samples_in_estimate, 0);
    assert!(partial.raw.is_none() && partial.controlled.is_none());
    split.advance(3, || false, affine);
    let one = split.assess_linear_control_variate(&control).unwrap();
    assert_eq!(one.samples_accepted, 7);
    assert_eq!(one.samples_in_estimate, 5);
    assert_eq!(one.completed_replicates, 1);
    assert_eq!(one.controlled.unwrap().standard_error, None);
    let identity = fs_blake3::hash_domain("lhs-control-test", b"affine");
    let bytes = split.checkpoint(identity).unwrap();
    let mut polls = 0;
    assert_eq!(
        split.assess_linear_control_variate_interruptible(&control, || {
            polls += 1;
            polls == 5
        }),
        Err(UqControlError::Cancelled)
    );
    assert_eq!(split.checkpoint(identity).unwrap(), bytes);
    let mut restored = LhsExecution::restore(&p, identity, &bytes).unwrap();
    let mut calls = 0;
    restored.advance(
        15,
        || false,
        |x| {
            calls += 1;
            affine(x)
        },
    );
    assert_eq!(calls, 8);
    assert_eq!(
        restored.assess_linear_control_variate(&control).unwrap(),
        expected
    );
    assert_eq!(
        restored.checkpoint(identity).unwrap(),
        whole.checkpoint(identity).unwrap()
    );
}

#[test]
fn g0_controls_bind_layout_physical_law_and_preobservation_state() {
    let p = plan();
    let mut run = LhsExecution::new(&p).unwrap();
    for gradient in [&[][..], &[2.0][..], &[2.0, f64::NAN][..]] {
        assert_eq!(
            run.freeze_linear_control_variate(gradient).unwrap_err(),
            UqControlError::InvalidGradient
        );
    }
    let control = run.freeze_linear_control_variate(&[2.0, 3.0]).unwrap();
    for variant in 0..5 {
        let mut changed = p.clone();
        match variant {
            0 => changed.seed += 1,
            1 => {
                changed.replicates = 5;
                changed.samples_per_replicate = 3;
            }
            2 => changed.parameters[0] = ParameterUncertainty::uniform("power", -3.0, 4.0, "W"),
            3 => changed.parameters[0] = ParameterUncertainty::uniform("power", -2.0, 4.0, "kW"),
            _ => changed.compliance_threshold = Some(29.0),
        }
        assert_eq!(
            LhsExecution::new(&changed)
                .unwrap()
                .assess_linear_control_variate(&control),
            Err(UqControlError::PlanMismatch)
        );
    }
    run.advance(1, || false, affine);
    assert_eq!(
        run.freeze_linear_control_variate(&[2.0, 3.0]).unwrap_err(),
        UqControlError::AlreadySampled
    );
    run.advance(1, || false, |_| Err::<f64, _>("physical refusal"));
    assert_eq!(
        run.freeze_linear_control_variate(&[2.0, 3.0]).unwrap_err(),
        UqControlError::RefusedExecution
    );
    assert_eq!(
        run.assess_linear_control_variate(&control),
        Err(UqControlError::RefusedExecution)
    );
}

#[test]
fn g3_wrong_fixed_coefficients_can_increase_variance_and_zero_control_is_identity() {
    let mut run = LhsExecution::new(&plan()).unwrap();
    let zero = run.freeze_linear_control_variate(&[0.0, 0.0]).unwrap();
    let wrong = run.freeze_linear_control_variate(&[-2.0, -3.0]).unwrap();
    run.advance(15, || false, affine);
    let unchanged = run.assess_linear_control_variate(&zero).unwrap();
    assert_eq!(unchanged.raw, unchanged.controlled);
    assert_eq!(
        unchanged.raw_replicate_means,
        unchanged.controlled_replicate_means
    );
    let worse = run.assess_linear_control_variate(&wrong).unwrap();
    assert!((worse.variance_ratio.unwrap() - 4.0).abs() < 1e-10);
}
