use super::*;
use fs_assimilate::groups::stream::design::SharedNoiseProbe;

fn campaign(cx: &Cx<'_>) -> (SharedNoiseModel, SharedNoiseStream) {
    let model = offset(cx, "fixed-fit", 4.0);
    let mut stream = SharedNoiseStream::new(&Belief::scalar(0.0, 10.0).unwrap(), &model, 64, cx).unwrap();
    for i in 0..16 {
        stream.assimilate(&batch(&model, vec![record(&model, &format!("past-{i}"), vec![1.0])], cx), cx).unwrap();
    }
    (model, stream)
}

fn probe(model: &SharedNoiseModel, row: &str, h: f64, j: f64, noise: f64) -> SharedNoiseProbe {
    SharedNoiseProbe::new(model, RowId::try_new(row).unwrap(), vec![h], vec![j], noise).unwrap()
}

#[test]
fn planning_prefers_a_probe_that_breaks_the_retained_calibration_ambiguity() {
    with_budget(Budget::INFINITE, |cx| {
        let (model, stream) = campaign(cx);
        let before = stream.clone();
        let probes = [probe(&model, "repeat-direction", 1.0, 1.0, 0.25), probe(&model, "independent-reference", 1.0, 0.0, 0.25), probe(&model, "opposite-direction", 1.0, -1.0, 0.25)];
        let scores = stream.rank_probes(&[1.0], &probes, cx).unwrap();
        assert_eq!(scores[0].row().as_str(), "opposite-direction");
        assert_eq!(scores[1].row().as_str(), "independent-reference");
        assert_eq!(scores[2].row().as_str(), "repeat-direction");
        assert!(scores[0].variance_reduction() > 100.0 * scores[2].variance_reduction());
        assert!(scores.iter().all(|score| score.history_id() == stream.history_id()));
        assert_eq!(stream, before);
    });
}

#[test]
fn predicted_target_covariance_matches_acquisition_for_every_actual_outcome() {
    with_budget(Budget::INFINITE, |cx| {
        let (model, stream) = campaign(cx);
        for j in [-1.0, 0.0, 1.0] {
            let candidate = probe(&model, "next", 1.0, j, 0.25);
            let score = stream.rank_probes(&[1.0], &[candidate.clone()], cx).unwrap().remove(0);
            close(score.predicted_mean(), stream.joint().mean()[0] + j * stream.joint().mean()[1]);
            for measured in [-20.0, 1.0, 30.0] {
                let observation = candidate.record(&model, measured, "real-instrument").unwrap();
                let mut acquired = stream.clone();
                acquired.assimilate(&batch(&model, vec![observation], cx), cx).unwrap();
                close(acquired.joint().variance(0).unwrap(), score.posterior_target_variance());
                close(score.prior_target_variance() - score.posterior_target_variance(), score.variance_reduction());
            }
        }
    });
}

#[test]
fn target_direction_and_independent_noise_determine_the_ranking() {
    with_budget(Budget::INFINITE, |cx| {
        let model = offset(cx, "fit", 1.0);
        let prior = Belief::diagonal(vec![0.0, 0.0], &[4.0, 1.0], cx).unwrap();
        let stream = SharedNoiseStream::new(&prior, &model, 8, cx).unwrap();
        let make = |id: &str, h: Vec<f64>, noise| SharedNoiseProbe::new(&model, RowId::try_new(id).unwrap(), h, vec![0.0], noise).unwrap();
        let probes = [make("x-noisy", vec![1.0, 0.0], 100.0), make("y", vec![0.0, 1.0], 0.25), make("x-precise", vec![1.0, 0.0], 0.25)];
        let x = stream.rank_probes(&[1.0, 0.0], &probes, cx).unwrap();
        assert_eq!(x[0].row().as_str(), "x-precise");
        assert_eq!(x[1].row().as_str(), "x-noisy");
        assert_eq!(x[2].variance_reduction(), 0.0);
        let y = stream.rank_probes(&[0.0, 1.0], &probes, cx).unwrap();
        assert_eq!(y[0].row().as_str(), "y");
        close(x[0].posterior_target_variance(), 1.0 / (0.25 + 4.0));
    });
}

#[test]
fn sensor_unit_scaling_target_scaling_and_tie_order_are_consistent() {
    with_budget(Budget::INFINITE, |cx| {
        let (model, stream) = campaign(cx);
        let base = probe(&model, "base", 1.0, -1.0, 0.25);
        let scaled = probe(&model, "scaled", -7.0, 7.0, 49.0 * 0.25);
        let a = stream.rank_probes(&[1.0], &[base.clone()], cx).unwrap().remove(0);
        let b = stream.rank_probes(&[1.0], &[scaled], cx).unwrap().remove(0);
        close(a.variance_reduction(), b.variance_reduction());
        close(a.posterior_target_variance(), b.posterior_target_variance());
        close(-7.0 * a.predicted_mean(), b.predicted_mean());
        close(49.0 * a.innovation_variance(), b.innovation_variance());
        let target_scaled = stream.rank_probes(&[-3.0], &[base], cx).unwrap().remove(0);
        close(target_scaled.variance_reduction(), 9.0 * a.variance_reduction());
        close(target_scaled.posterior_target_variance(), 9.0 * a.posterior_target_variance());
        let left = probe(&model, "a", 1.0, 0.0, 0.25);
        let right = probe(&model, "z", 1.0, 0.0, 0.25);
        let order_a = stream.rank_probes(&[1.0], &[right.clone(), left.clone()], cx).unwrap();
        let order_b = stream.rank_probes(&[1.0], &[left, right], cx).unwrap();
        assert_eq!(order_a, order_b);
        assert_eq!(order_a[0].row().as_str(), "a");
    });
}

#[test]
fn near_exact_sensor_keeps_positive_conditional_variance_without_subtraction() {
    with_budget(Budget::INFINITE, |cx| {
        let model = offset(cx, "zero-error", 0.0);
        let stream = SharedNoiseStream::new(&Belief::scalar(0.0, 1.0).unwrap(), &model, 2, cx).unwrap();
        let sensor = probe(&model, "near-exact", 1.0, 0.0, 1e-20);
        let score = stream.rank_probes(&[1.0], &[sensor], cx).unwrap().remove(0);
        assert!(score.posterior_target_variance() > 0.0);
        assert!((score.posterior_target_variance() / 1e-20 - 1.0).abs() < 1e-14);
        // Subtraction of the rounded improvement would incorrectly give zero.
        assert_eq!(score.prior_target_variance() - score.variance_reduction(), 0.0);
        // An unrepresentable positive result must not become an exact-zero claim.
        let tiny_noise = probe(&model, "underflow", 1e154, 0.0, 1e-320);
        assert!(stream.rank_probes(&[1e-160], &[tiny_noise], cx).is_err());
        let constrained = SharedNoiseStream::new(&Belief::scalar(3.0, 0.0).unwrap(), &model, 2, cx).unwrap();
        let known = constrained.rank_probes(&[1.0], &[probe(&model, "known", 1.0, 1.0, 0.25)], cx).unwrap().remove(0);
        assert_eq!(known.variance_reduction(), 0.0);
        assert_eq!(known.posterior_target_variance(), 0.0);
        assert_eq!(known.predicted_mean(), 3.0);
    });
}

#[test]
fn invalid_or_already_acquired_probes_and_budget_refusals_publish_no_plan() {
    let (stream, candidate) = with_budget(Budget::INFINITE, |cx| {
        let (model, stream) = campaign(cx);
        let candidate = probe(&model, "new", 1.0, 1.0, 0.25);
        let prior = stream.clone();
        let other = offset(cx, "other-model", 4.0);
        assert!(matches!(stream.rank_probes(&[1.0], &[probe(&other, "new", 1.0, 1.0, 0.25)], cx), Err(AssimError::SharedSourceDeclaration { .. })));
        assert!(matches!(stream.rank_probes(&[1.0], &[probe(&model, "past-0", 1.0, 1.0, 0.25)], cx), Err(AssimError::DuplicateDatasetRow { .. })));
        assert!(matches!(stream.rank_probes(&[1.0], &[candidate.clone(), candidate.clone()], cx), Err(AssimError::DuplicateDatasetRow { .. })));
        for target in [vec![0.0], vec![f64::NAN], vec![1.0, 0.0]] {
            assert!(stream.rank_probes(&target, &[candidate.clone()], cx).is_err());
        }
        assert!(stream.rank_probes(&[1.0], &[], cx).is_err());
        assert!(stream.rank_probes(&[1.0], &vec![candidate.clone(); 257], cx).is_err());
        assert!(candidate.record(&other, 1.0, "instrument").is_err());
        assert!(candidate.record(&model, f64::NAN, "instrument").is_err());
        assert!(SharedNoiseProbe::new(&model, RowId::try_new("zero-noise").unwrap(), vec![1.0], vec![1.0], 0.0).is_err());
        assert_eq!(stream, prior);
        (stream, candidate)
    });
    for budget in [Budget::INFINITE.with_poll_quota(0), Budget::INFINITE.with_poll_quota(1), Budget::INFINITE.with_cost_quota(1)] {
        with_budget(budget, |cx| {
            assert!(stream.rank_probes(&[1.0], &[candidate.clone()], cx).is_err());
        });
    }
}
