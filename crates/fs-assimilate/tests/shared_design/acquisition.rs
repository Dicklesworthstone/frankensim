use super::*;
use fs_assimilate::AssimError;
use fs_assimilate::groups::stream::design::SharedNoiseProbe;
use fs_assimilate::groups::stream::design::batch::{CostedProbe, ProbeBatchBudget, ProbeBatchPlan, ProbeReading};

fn fixture(cx: &Cx<'_>) -> (SharedNoiseStream, SharedNoiseModel, ProbeBatchPlan) {
    let model = SharedNoiseModel::new("fit", vec!["offset".into()], vec![vec![1.0]], cx).unwrap();
    let stream = SharedNoiseStream::new(&Belief::scalar(0.0, 2.0).unwrap(), &model, 8, cx).unwrap();
    let probes: Vec<_> = [("left", 1.0), ("right", -1.0)].into_iter().map(|(row, j)|
        CostedProbe::new(SharedNoiseProbe::new(&model, RowId::try_new(row).unwrap(), vec![1.0], vec![j], 0.25).unwrap(), 1).unwrap()
    ).collect();
    let plan = stream.plan_probe_batch(&[1.0], &probes, ProbeBatchBudget::new(2, 2).unwrap(), cx).unwrap();
    (stream, model, plan)
}
fn readings(plan: &ProbeBatchPlan) -> Vec<ProbeReading> {
    plan.steps().iter().enumerate().map(|(index, step)|
        ProbeReading::new(step.candidate().probe().row().clone(), 3.0 - 2.0 * index as f64, "instrument").unwrap()
    ).collect()
}

#[test]
fn real_row_mapping_is_order_independent_and_matches_normal_stream_updates() {
    with_budget(Budget::INFINITE, |cx| {
        let (mut stream, model, plan) = fixture(cx);
        let measured = readings(&plan);
        let mut reference = stream.clone();
        let rows = plan.steps().iter().zip(&measured).map(|(s, r)| s.candidate().probe().record(&model, r.value(), r.instrument()).unwrap()).collect();
        reference.assimilate(&JacobianBatch::new(model.clone(), rows, cx).unwrap(), cx).unwrap();
        let reversed: Vec<_> = measured.iter().cloned().rev().collect();
        plan.check_current(&stream).unwrap();
        let update = plan.assimilate_readings(&mut stream, &model, &reversed, cx).unwrap();
        assert_eq!(update.accepted_indices(), &[0, 1]); // selection order, not input order
        assert_eq!(stream, reference);
        close(stream.physical_belief(cx).unwrap().variance(0).unwrap(), plan.posterior_target_variance());
        let before = stream.clone();
        assert!(plan.check_current(&stream).is_err());
        assert!(plan.assimilate_readings(&mut stream, &model, &measured, cx).is_err());
        assert_eq!(stream, before);
    });
}

#[test]
fn missing_duplicate_unselected_and_foreign_model_readings_refuse_atomically() {
    with_budget(Budget::INFINITE, |cx| {
        let (mut stream, model, plan) = fixture(cx);
        let before = stream.clone();
        let correct = readings(&plan);
        let unknown = ProbeReading::new(RowId::try_new("unselected").unwrap(), 2.0, "instrument").unwrap();
        for wrong in [vec![], vec![correct[0].clone()], vec![correct[0].clone(), correct[0].clone()],
            vec![correct[0].clone(), unknown], vec![correct[0].clone(), correct[1].clone(), correct[0].clone()]] {
            assert!(plan.assimilate_readings(&mut stream, &model, &wrong, cx).is_err());
            assert_eq!(stream, before);
        }
        let other = SharedNoiseModel::new("other-fit", vec!["offset".into()], vec![vec![1.0]], cx).unwrap();
        assert!(plan.assimilate_readings(&mut stream, &other, &correct, cx).is_err());
        assert_eq!(stream, before);
        let mut small = SharedNoiseStream::new(&Belief::scalar(0.0, 2.0).unwrap(), &model, 1, cx).unwrap();
        assert_eq!(small.history_id(), plan.history_id());
        let small_before = small.clone();
        assert!(plan.assimilate_readings(&mut small, &model, &correct, cx).is_err());
        assert_eq!(small, small_before);
        assert!(ProbeReading::new(correct[0].row().clone(), f64::NAN, "instrument").is_err());
        assert!(ProbeReading::new(correct[0].row().clone(), 1.0, "").is_err());
    });
}

#[test]
fn unrelated_arrival_invalidates_plan_without_consuming_planned_rows() {
    with_budget(Budget::INFINITE, |cx| {
        let (mut stream, model, plan) = fixture(cx);
        let row = JacobianRecord::new(Observation::new(vec![1.0], 0.0, 0.25, "sensor").unwrap(),
            RowId::try_new("unrelated").unwrap(), &model, vec![1.0]).unwrap();
        stream.assimilate(&JacobianBatch::new(model.clone(), vec![row], cx).unwrap(), cx).unwrap();
        let before = stream.clone();
        assert!(plan.assimilate_readings(&mut stream, &model, &readings(&plan), cx).is_err());
        assert_eq!(stream, before);
        assert_eq!(stream.accepted_rows(), 1);
    });
}

#[test]
fn failed_numerical_update_and_cancelled_preparation_allow_exact_retry() {
    let (original, model, plan, correct) = with_budget(Budget::INFINITE, |cx| {
        let (stream, model, plan) = fixture(cx);
        let correct = readings(&plan);
        let huge = vec![correct[0].clone(), ProbeReading::new(correct[1].row().clone(), 1e200, "instrument").unwrap()];
        let mut failed = stream.clone();
        assert!(plan.assimilate_readings(&mut failed, &model, &huge, cx).is_err());
        assert_eq!(failed, stream);
        (stream, model, plan, correct)
    });
    for budget in [Budget::INFINITE.with_poll_quota(0), Budget::INFINITE.with_poll_quota(1), Budget::INFINITE.with_cost_quota(1)] {
        with_budget(budget, |cx| {
            let mut stream = original.clone();
            assert!(matches!(plan.assimilate_readings(&mut stream, &model, &correct, cx),
                Err(AssimError::BudgetRefused(_)) | Err(AssimError::Cancelled { .. })));
            assert_eq!(stream, original);
        });
    }
    with_budget(Budget::INFINITE, |cx| {
        let mut stream = original;
        plan.assimilate_readings(&mut stream, &model, &correct, cx).unwrap();
        assert_eq!(stream.accepted_rows(), 2);
        close(stream.physical_belief(cx).unwrap().variance(0).unwrap(), plan.posterior_target_variance());
    });
}

#[test]
fn empty_selection_has_a_checked_empty_acquisition() {
    with_budget(Budget::INFINITE, |cx| {
        let (mut stream, model, plan) = fixture(cx);
        let candidate = plan.steps()[0].candidate().clone();
        let empty = stream.plan_probe_batch(&[1.0], &[candidate], ProbeBatchBudget::new(0, 0).unwrap(), cx).unwrap();
        let before = stream.clone();
        let update = empty.assimilate_readings(&mut stream, &model, &[], cx).unwrap();
        assert!(update.accepted_indices().is_empty()); assert!(update.diagnostics().is_none());
        assert_eq!(stream, before);
        assert!(empty.assimilate_readings(&mut stream, &model, &readings(&plan), cx).is_err());
        assert_eq!(stream, before);
    });
}

#[test]
fn unselectable_candidates_do_not_trigger_unused_numerical_projections() {
    with_budget(Budget::INFINITE, |cx| {
        let (stream, model, previous) = fixture(cx);
        let giant = CostedProbe::new(SharedNoiseProbe::new(&model, RowId::try_new("unaffordable").unwrap(),
            vec![f64::MAX], vec![0.0], 0.25).unwrap(), 100).unwrap();
        let probes = vec![previous.steps()[0].candidate().clone(), giant];
        let plan = stream.plan_probe_batch(&[1.0], &probes, ProbeBatchBudget::new(2, 1).unwrap(), cx).unwrap();
        assert_eq!(plan.steps().len(), 1); assert_eq!(plan.total_cost(), 1);
        let empty = stream.plan_probe_batch(&[1.0], &probes, ProbeBatchBudget::new(0, 100).unwrap(), cx).unwrap();
        assert!(empty.steps().is_empty());
    });
}
