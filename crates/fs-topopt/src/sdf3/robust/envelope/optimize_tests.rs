use super::*;
use super::super::tests::fixture;
use crate::pipeline::LoadCase;
use crate::{SolveBudget, SolveProgress};

fn options() -> EnvelopeDesignOptions3 {
    EnvelopeDesignOptions3 { volume_limit: 0.5, ..Default::default() }
}
fn policy(study: &mut CutDensityStudy3<fs_cutfem::elastic3::CutElasticity3>,
    rho: &[f64], y: &[f64], z: &[f64], control: &mut SolveControl<'_>) -> LoadEnvelope3 {
    let e = study.evaluate_robust(rho, &[LoadCase { force: y, weight: 0.5 }, LoadCase { force: z, weight: 0.5 }],
        options().robustness, control).unwrap();
    LoadEnvelope3::new(e.eroded.objective.case_compliances, 0.05, 2).unwrap()
}

#[test]
fn real_worst_load_design_improves_and_replays_its_retained_physics() {
    let (mut study,y,z) = fixture(2.0);
    let rho = vec![0.3;study.cells()];
    let keys = study.operator().cell_keys(); let scales = study.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let p = policy(&mut study,&rho,&y,&z,&mut control);
    let forces = [&y[..],&z[..]];
    let mut session = EnvelopeDesignStudy3::new(&mut study,&forces,&rho,p.clone(),options(),&mut control).unwrap();
    assert!(session.gradient_check().numerical.passed());
    let initial = session.accepted().clone();
    let outcome = session.run(30).unwrap();
    assert!(session.history().len()>1);
    let best = session.best_feasible().unwrap().clone();
    assert!(best.evaluation.envelope.value < initial.evaluation.envelope.value);
    assert!(best.evaluation.envelope.maximum < initial.evaluation.envelope.maximum);
    assert_eq!(outcome.objective,session.accepted().evaluation.envelope.value);
    assert!(best.evaluation.physical.dilated_volume_fraction <= options().volume_limit*(1.0+options().optimizer.tolerance));
    drop(session);
    let fresh = study.evaluate_robust_envelope(&best.rho,&forces,options().robustness,&p,&mut control).unwrap();
    assert_eq!(fresh.envelope,best.evaluation.envelope);
    assert_eq!(fresh.physical.eroded.objective.displacements,best.evaluation.physical.eroded.objective.displacements);
    assert_eq!(study.operator().cell_keys(),keys); assert_eq!(study.operator().scales(),scales);
}

#[test]
fn segmented_runs_preserve_the_trajectory_and_cumulative_work() {
    let run = |chunks: &[usize]| {
        let (mut study,y,z) = fixture(1.0); let rho = vec![0.3;study.cells()];
        let mut poll = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
        let p = policy(&mut study,&rho,&y,&z,&mut control); let forces = [&y[..],&z[..]];
        let mut session = EnvelopeDesignStudy3::new(&mut study,&forces,&rho,p,options(),&mut control).unwrap();
        let initial_work = session.work();
        session.run(0).unwrap(); assert_eq!(initial_work,session.work());
        for &chunk in chunks {
            if session.run(chunk).unwrap().stop != ProjectedAlStop::IterationLimit { break; }
        }
        (session.history().to_vec(),session.accepted().rho.clone(),session.work(),session.optimizer_work())
    };
    let whole = run(&[12]); let split = run(&[5,7]);
    assert!(whole.0.len()>1); assert_eq!(whole,split);
}

#[test]
fn late_cancellation_preserves_the_accepted_point_and_its_incumbent() {
    let (mut study,y,z) = fixture(1.0); let rho = vec![0.3;study.cells()];
    let scales = study.operator().scales().to_vec();
    let mut poll = |p: SolveProgress| if p.stage=="sdf3-envelope-design-accepted" { ControlFlow::Break(()) } else { ControlFlow::Continue(()) };
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let p = policy(&mut study,&rho,&y,&z,&mut control); let forces = [&y[..],&z[..]];
    let mut session = EnvelopeDesignStudy3::new(&mut study,&forces,&rho,p,options(),&mut control).unwrap();
    assert!(matches!(session.run(12),Err(ProjectedAlError::Evaluation(EvaluationStop::Cancelled))));
    assert_eq!(session.history().len(),2);
    assert_eq!(session.accepted().rho,session.state.point());
    assert_eq!(session.history().last().unwrap().envelope,session.accepted().evaluation.envelope.value);
    assert!(session.best_feasible().is_some());
    drop(session); assert_eq!(study.operator().scales(),scales);
}

#[test]
fn evaluation_limit_does_not_renew_work_or_publish_a_trial() {
    let (mut study,y,z) = fixture(1.0); let rho = vec![0.3;study.cells()];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let p = policy(&mut study,&rho,&y,&z,&mut control); let forces = [&y[..],&z[..]];
    let mut opts = options(); opts.optimizer.max_evaluations = 1;
    let mut session = EnvelopeDesignStudy3::new(&mut study,&forces,&rho,p,opts,&mut control).unwrap();
    let work = session.work();
    assert_eq!(session.run(20).unwrap().stop,ProjectedAlStop::EvaluationLimit);
    assert_eq!(session.work(),work);
    assert_eq!(session.history().len(),1); assert_eq!(session.accepted().rho,rho);
}

#[test]
fn exhausted_physics_and_invalid_case_shapes_refuse_without_mutation() {
    let (mut study,y,z) = fixture(1.0); let rho = vec![0.3;study.cells()];
    let scales = study.operator().scales().to_vec();
    let p = LoadEnvelope3::new(vec![1000.0;2],0.05,2).unwrap();
    let forces = [&y[..],&z[..]];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: 1, ..Default::default() }, &mut poll);
    assert!(matches!(EnvelopeDesignStudy3::new(&mut study,&forces,&rho,p.clone(),options(),&mut control),
        Err(ProjectedAlError::Evaluation(EvaluationStop::TotalBudget {..}))));
    assert_eq!(control.work().linear_iterations,1);
    assert_eq!(study.operator().scales(),scales);
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    assert!(matches!(EnvelopeDesignStudy3::new(&mut study,&forces[..1],&rho,p,options(),&mut control),Err(ProjectedAlError::Invalid(_))));
    assert_eq!(control.work().linear_solves,0);
}
