use super::*;
use crate::sdf3::PhysicalRegion3;
use std::cell::Cell;

#[test]
fn infeasible_start_recovers_actual_dilated_volume_without_installing_optimizer_state() {
    let (study, y, z) = fixture(2.0);
    let mut regions = vec![PhysicalRegion3::Design; study.cells()];
    regions[0] = PhysicalRegion3::Solid;
    regions[1] = PhysicalRegion3::Void;
    let mut study = study.with_physical_regions(regions).unwrap();
    let raw = vec![0.8; study.cells()];
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    study.evaluate(&raw, &[LoadCase { force: &z, weight: 1.0 }], &mut control).unwrap();
    let scales = study.operator().scales().to_vec();
    let original = study.params();
    let p = policy(&mut study, &raw, &y, &z, &mut control);
    let forces = [&y[..], &z[..]];
    let mut session = EnvelopeDesignStudy3::new(&mut study, &forces, &raw, p.clone(), options(), &mut control).unwrap();
    assert!(session.best_feasible().is_none());
    let checkpoint = session.state.checkpoint();
    let history = session.history().to_vec();
    let recovered = session.recover_feasible().unwrap();
    assert_ne!(recovered.rho, raw);
    assert!(recovered.evaluation.physical.dilated_volume_fraction <= options().volume_limit);
    assert_eq!(session.state.checkpoint(), checkpoint);
    assert_eq!(session.history(), history);
    assert!(session.best_feasible().is_none());
    for field in [&recovered.evaluation.physical.eroded.projected_rho,
        &recovered.evaluation.physical.nominal_rho, &recovered.evaluation.physical.dilated_rho] {
        assert_eq!(field[0], 1.0); assert_eq!(field[1], 0.0);
    }
    drop(session);
    let fresh = study.evaluate_robust_envelope(&recovered.rho, &forces, options().robustness, &p, &mut control).unwrap();
    assert_eq!(recovered.evaluation.envelope, fresh.envelope);
    assert_eq!(recovered.evaluation.physical.eroded.objective.displacements, fresh.physical.eroded.objective.displacements);
    assert_eq!(study.operator().scales(), scales);
    assert_eq!(study.params().eta, original.eta);
    assert_eq!(study.params().beta, original.beta);
}

#[test]
fn interleaved_recovery_does_not_restart_the_following_al_segment() {
    let run = |recover: bool| {
        let (mut study, y, z) = fixture(2.0);
        let raw = vec![0.3; study.cells()];
        let mut poll = |_| ControlFlow::Continue(());
        let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
        let p = policy(&mut study, &raw, &y, &z, &mut control);
        let forces = [&y[..], &z[..]];
        let mut session = EnvelopeDesignStudy3::new(&mut study, &forces, &raw, p, options(), &mut control).unwrap();
        session.run(5).unwrap();
        if recover { session.recover_feasible().unwrap(); }
        session.run(7).unwrap();
        (session.state.checkpoint(), session.history().to_vec(),
            session.accepted().evaluation.physical.eroded.objective.displacements.clone(), session.work())
    };
    let plain = run(false); let recovered = run(true);
    assert_eq!(plain.0, recovered.0);
    assert_eq!(plain.1, recovered.1);
    assert_eq!(plain.2, recovered.2);
    assert!(recovered.3.linear_iterations > plain.3.linear_iterations);
}

#[test]
fn recovery_cancellation_preserves_cached_fields_model_and_nominal_stiffness() {
    for stage in ["sdf3-volume-restoration", "sdf3-envelope-adjoint", "sdf3-envelope-recovery-publish"] {
        let (mut study, y, z) = fixture(2.0);
        let raw = vec![0.8; study.cells()];
        let scales = study.operator().scales().to_vec();
        let original = study.params();
        let recovering = Cell::new(false);
        let mut poll = |p: SolveProgress| if recovering.get() && p.stage == stage {
            ControlFlow::Break(())
        } else { ControlFlow::Continue(()) };
        let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
        let forces = [&y[..], &z[..]];
        let p = LoadEnvelope3::new(vec![1000.0; 2], 0.05, 2).unwrap();
        let mut session = EnvelopeDesignStudy3::new(&mut study, &forces, &raw, p, options(), &mut control).unwrap();
        let checkpoint = session.state.checkpoint();
        let fields = session.accepted().evaluation.physical.eroded.objective.displacements.clone();
        recovering.set(true);
        assert!(matches!(session.recover_feasible(), Err(EvaluationStop::Cancelled)), "{stage}");
        assert_eq!(session.state.checkpoint(), checkpoint);
        assert_eq!(session.accepted().evaluation.physical.eroded.objective.displacements, fields);
        assert_eq!(session.history().len(), 1);
        assert!(session.best_feasible().is_none());
        drop(session);
        assert_eq!(study.operator().scales(), scales);
        assert_eq!(study.params().eta, original.eta);
    }
}

#[test]
fn recovery_cannot_renew_the_original_physical_iteration_budget() {
    let (mut study, y, z) = fixture(2.0);
    let raw = vec![0.8; study.cells()];
    let forces = [&y[..], &z[..]];
    let p = LoadEnvelope3::new(vec![1000.0; 2], 0.05, 2).unwrap();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let session = EnvelopeDesignStudy3::new(&mut study, &forces, &raw, p.clone(), options(), &mut control).unwrap();
    let budget = session.work().linear_iterations + 1;
    drop(session);
    let scales = study.operator().scales().to_vec();
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget { total_iterations: budget, ..Default::default() }, &mut poll);
    let mut session = EnvelopeDesignStudy3::new(&mut study, &forces, &raw, p, options(), &mut control).unwrap();
    let checkpoint = session.state.checkpoint();
    assert!(matches!(session.recover_feasible(), Err(EvaluationStop::TotalBudget { .. })));
    assert_eq!(session.work().linear_iterations, budget);
    assert_eq!(session.state.checkpoint(), checkpoint);
    drop(session);
    assert_eq!(study.operator().scales(), scales);
}

#[test]
fn a_rejected_declared_floor_does_not_publish_a_feasible_point() {
    let (mut study, y, z) = fixture(2.0);
    let raw = vec![0.9; study.cells()];
    let forces = [&y[..], &z[..]];
    let p = LoadEnvelope3::new(vec![1000.0; 2], 0.05, 2).unwrap();
    let opts = EnvelopeDesignOptions3 { density_floor: 0.8, ..options() };
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let mut session = EnvelopeDesignStudy3::new(&mut study, &forces, &raw, p, opts, &mut control).unwrap();
    let checkpoint = session.state.checkpoint();
    assert!(matches!(session.recover_feasible(), Err(EvaluationStop::Breakdown { stage: "sdf3-volume-restoration" })));
    assert_eq!(session.state.checkpoint(), checkpoint);
    assert!(session.best_feasible().is_none());
}
