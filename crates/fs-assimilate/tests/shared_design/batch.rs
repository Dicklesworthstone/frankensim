use super::*;
use fs_assimilate::groups::stream::design::SharedNoiseProbe;
use fs_assimilate::groups::stream::design::batch::{CostedProbe, ProbeBatchBudget, ProbeBatchStop};
use fs_assimilate::AssimError;

fn model(cx: &Cx<'_>, variance: f64) -> SharedNoiseModel {
    SharedNoiseModel::new("batch-calibration", vec!["offset".into()], vec![vec![variance]], cx).unwrap()
}
fn candidate(model: &SharedNoiseModel, row: &str, h: Vec<f64>, j: f64, noise: f64, cost: u32) -> CostedProbe {
    CostedProbe::new(SharedNoiseProbe::new(model, RowId::try_new(row).unwrap(), h, vec![j], noise).unwrap(), cost).unwrap()
}
fn target_variance(belief: &Belief, target: &[f64]) -> f64 {
    let mut sum = 0.0;
    for i in 0..target.len() { for j in 0..target.len() { sum += target[i] * belief.covariance()[i][j] * target[j]; } }
    sum
}

#[test]
fn joint_selection_rejects_redundant_top_rankings_and_matches_every_acquired_prefix() {
    with_budget(Budget::INFINITE, |cx| {
        let model = model(cx, 0.0);
        let prior = Belief::diagonal(vec![0.0, 0.0], &[2.0, 1.0], cx).unwrap();
        let mut stream = SharedNoiseStream::new(&prior, &model, 4, cx).unwrap();
        let candidates = vec![
            candidate(&model, "a", vec![1.0, 0.0], 0.0, 0.01, 1),
            candidate(&model, "b", vec![1.0, 0.0], 0.0, 0.01, 1),
            candidate(&model, "c", vec![0.0, 1.0], 0.0, 0.01, 1),
        ];
        let target = [1.0, 1.0];
        let independent = stream.rank_probes(&target, &candidates.iter().map(|c| c.probe().clone()).collect::<Vec<_>>(), cx).unwrap();
        assert_eq!(independent[0].row().as_str(), "a");
        assert_eq!(independent[1].row().as_str(), "b");
        let before = stream.clone();
        let plan = stream.plan_probe_batch(&target, &candidates, ProbeBatchBudget::new(2, 2).unwrap(), cx).unwrap();
        assert_eq!(stream, before);
        assert_eq!(plan.steps().iter().map(|s| s.candidate_index()).collect::<Vec<_>>(), vec![0, 2]);
        assert_eq!(plan.total_cost(), 2);
        assert_eq!(plan.stop(), ProbeBatchStop::ReadingLimit);
        close(plan.steps()[0].variance_reduction(), independent[0].variance_reduction());
        close(plan.posterior_target_variance(), 1.0 / (0.5 + 100.0) + 1.0 / (1.0 + 100.0));
        for (i, step) in plan.steps().iter().enumerate() {
            let reading = step.candidate().probe().record(&model, 3.0 - i as f64, "instrument").unwrap();
            stream.assimilate(&JacobianBatch::new(model.clone(), vec![reading], cx).unwrap(), cx).unwrap();
            close(step.posterior_target_variance(), target_variance(&stream.physical_belief(cx).unwrap(), &target));
        }
        assert!(plan.posterior_target_variance() < 0.02);
        assert!(1.0 + 1.0 / (0.5 + 200.0) > 1.0); // naive a,b selection
    });
}

#[test]
fn acquisition_cost_count_capacity_and_zero_information_are_explicit_stops() {
    with_budget(Budget::INFINITE, |cx| {
        let model = model(cx, 0.0);
        let prior = Belief::diagonal(vec![0.0, 0.0], &[2.0, 1.0], cx).unwrap();
        let stream = SharedNoiseStream::new(&prior, &model, 4, cx).unwrap();
        let candidates = vec![
            candidate(&model, "expensive", vec![1.0, 0.0], 0.0, 0.01, 3),
            candidate(&model, "cheap-x", vec![1.0, 0.0], 0.0, 0.25, 1),
            candidate(&model, "cheap-y", vec![0.0, 1.0], 0.0, 0.01, 1),
        ];
        let plan = stream.plan_probe_batch(&[1.0, 1.0], &candidates, ProbeBatchBudget::new(3, 2).unwrap(), cx).unwrap();
        assert_eq!(plan.steps().iter().map(|s| s.candidate_index()).collect::<Vec<_>>(), vec![1, 2]);
        assert_eq!(plan.stop(), ProbeBatchStop::CostLimit);
        assert_eq!(plan.total_cost(), 2);
        for (count, cost, stop) in [(0, 5, ProbeBatchStop::ReadingLimit), (3, 0, ProbeBatchStop::CostLimit)] {
            let empty = stream.plan_probe_batch(&[1.0, 1.0], &candidates, ProbeBatchBudget::new(count, cost).unwrap(), cx).unwrap();
            assert!(empty.steps().is_empty()); assert_eq!(empty.stop(), stop);
            assert_eq!(empty.prior_target_variance(), empty.posterior_target_variance());
        }
        let limited = SharedNoiseStream::new(&prior, &model, 1, cx).unwrap();
        let plan = limited.plan_probe_batch(&[1.0, 1.0], &candidates, ProbeBatchBudget::new(3, 10).unwrap(), cx).unwrap();
        assert_eq!(plan.steps().len(), 1); assert_eq!(plan.stop(), ProbeBatchStop::RowCapacity);
        let plan = stream.plan_probe_batch(&[1.0, 1.0], &candidates[..1], ProbeBatchBudget::new(3, 10).unwrap(), cx).unwrap();
        assert_eq!(plan.stop(), ProbeBatchStop::CandidatesExhausted);
        let exact = SharedNoiseStream::new(&Belief::diagonal(vec![1.0, 0.0], &[0.0, 1.0], cx).unwrap(), &model, 4, cx).unwrap();
        let plan = exact.plan_probe_batch(&[1.0, 0.0], &candidates, ProbeBatchBudget::new(3, 10).unwrap(), cx).unwrap();
        assert_eq!(plan.stop(), ProbeBatchStop::NoInformation);
        assert_eq!(plan.posterior_target_variance(), 0.0);
    });
}

#[test]
fn conditional_forecasts_retain_learned_shared_error_and_ignore_realized_values() {
    with_budget(Budget::INFINITE, |cx| {
        let model = model(cx, 4.0);
        let mut stream = SharedNoiseStream::new(&Belief::scalar(0.0, 10.0).unwrap(), &model, 16, cx).unwrap();
        for i in 0..8 {
            let probe = candidate(&model, &format!("seen-{i}"), vec![1.0], 1.0, 0.25, 1);
            let row = probe.probe().record(&model, 1.0, "instrument").unwrap();
            stream.assimilate(&JacobianBatch::new(model.clone(), vec![row], cx).unwrap(), cx).unwrap();
        }
        let candidates = vec![
            candidate(&model, "same", vec![1.0], 1.0, 0.25, 1),
            candidate(&model, "opposite", vec![1.0], -1.0, 0.25, 1),
            candidate(&model, "independent", vec![1.0], 0.0, 0.25, 1),
        ];
        let plan = stream.plan_probe_batch(&[1.0], &candidates, ProbeBatchBudget::new(2, 2).unwrap(), cx).unwrap();
        assert_eq!(plan.steps()[0].candidate_index(), 1);
        for value in [-3.0, 0.0, 8.0] {
            let mut actual = stream.clone();
            for step in plan.steps() {
                let row = step.candidate().probe().record(&model, value, "instrument").unwrap();
                actual.assimilate(&JacobianBatch::new(model.clone(), vec![row], cx).unwrap(), cx).unwrap();
                close(actual.physical_belief(cx).unwrap().variance(0).unwrap(), step.posterior_target_variance());
            }
        }
        assert_eq!(stream.accepted_rows(), 8);
        assert_eq!(plan.history_id(), stream.history_id());
    });
}

#[test]
fn selection_is_candidate_permutation_stable_and_respects_measurement_units() {
    with_budget(Budget::INFINITE, |cx| {
        let model = model(cx, 1.0);
        let prior = Belief::diagonal(vec![0.2, -0.3], &[2.0, 1.0], cx).unwrap();
        let stream = SharedNoiseStream::new(&prior, &model, 8, cx).unwrap();
        let base = vec![candidate(&model, "a", vec![1.0, 0.0], 1.0, 0.25, 1),
            candidate(&model, "b", vec![0.0, 1.0], -0.5, 0.5, 2),
            candidate(&model, "c", vec![1.0, 1.0], 0.25, 1.0, 2)];
        let budget = ProbeBatchBudget::new(2, 4).unwrap();
        let expected = stream.plan_probe_batch(&[1.0, 0.5], &base, budget, cx).unwrap();
        for scale in [-2.0, 0.5, 1.0, 16.0] {
            let probes: Vec<_> = base.iter().rev().map(|c| candidate(&model, c.probe().row().as_str(),
                c.probe().physical_operator().iter().map(|x| x * scale).collect(),
                c.probe().jacobian()[0] * scale, c.probe().noise_variance() * scale * scale, c.cost())).collect();
            let actual = stream.plan_probe_batch(&[1.0, 0.5], &probes, budget, cx).unwrap();
            assert_eq!(actual.steps().iter().map(|s| s.candidate().probe().row()).collect::<Vec<_>>(),
                expected.steps().iter().map(|s| s.candidate().probe().row()).collect::<Vec<_>>());
            close(actual.posterior_target_variance(), expected.posterior_target_variance());
        }
    });
}

#[test]
fn near_exact_reading_retains_positive_conditional_variance() {
    with_budget(Budget::INFINITE, |cx| {
        let model = model(cx, 0.0);
        let stream = SharedNoiseStream::new(&Belief::scalar(0.0, 1.0).unwrap(), &model, 4, cx).unwrap();
        let probes = vec![candidate(&model, "a", vec![1.0], 0.0, 1.0e-24, 1),
            candidate(&model, "b", vec![1.0], 0.0, 1.0e-24, 1)];
        let plan = stream.plan_probe_batch(&[1.0], &probes, ProbeBatchBudget::new(2, 2).unwrap(), cx).unwrap();
        assert!(plan.posterior_target_variance() > 0.0);
        assert!((plan.posterior_target_variance() / 5.0e-25 - 1.0).abs() < 1.0e-12);
    });
}

#[test]
fn whole_candidate_admission_and_compute_refusals_do_not_change_stream() {
    let (stream, probes) = with_budget(Budget::INFINITE, |cx| {
        let model = model(cx, 0.0);
        let stream = SharedNoiseStream::new(&Belief::scalar(0.0, 1.0).unwrap(), &model, 4, cx).unwrap();
        let probe = candidate(&model, "a", vec![1.0], 0.0, 0.25, 1);
        assert!(CostedProbe::new(probe.probe().clone(), 0).is_err());
        assert!(ProbeBatchBudget::new(257, 1).is_err());
        assert!(stream.plan_probe_batch(&[1.0], &[probe.clone(), probe.clone()], ProbeBatchBudget::new(1, 0).unwrap(), cx).is_err());
        let other = SharedNoiseModel::new("other", vec!["offset".into()], vec![vec![0.0]], cx).unwrap();
        let wrong = candidate(&other, "b", vec![1.0], 0.0, 0.25, 100);
        assert!(stream.plan_probe_batch(&[1.0], &[probe.clone(), wrong], ProbeBatchBudget::new(1, 1).unwrap(), cx).is_err());
        for target in [vec![], vec![1.0, 0.0], vec![0.0], vec![f64::NAN]] {
            assert!(stream.plan_probe_batch(&target, &[probe.clone()], ProbeBatchBudget::new(1, 1).unwrap(), cx).is_err());
        }
        (stream, vec![probe])
    });
    for budget in [Budget::INFINITE.with_poll_quota(0), Budget::INFINITE.with_poll_quota(2), Budget::INFINITE.with_cost_quota(1)] {
        with_budget(budget, |cx| {
            let before = stream.clone();
            assert!(matches!(stream.plan_probe_batch(&[1.0], &probes, ProbeBatchBudget::new(1, 1).unwrap(), cx),
                Err(AssimError::BudgetRefused(_)) | Err(AssimError::Cancelled { .. })));
            assert_eq!(stream, before);
        });
    }
}
