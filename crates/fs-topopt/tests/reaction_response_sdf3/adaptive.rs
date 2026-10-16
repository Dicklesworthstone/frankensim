//! Real mixed-observation fitting through new-grid baseline and accepted steps.
use super::*;
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_topopt::sdf3::response::{ResponseEvaluation3, ResponseOptions3};
use fs_topopt::sdf3::response::refinement::{ReferenceResponseCase3, ReferenceResponseTarget3};

fn observe(p: [f64; 3]) -> [f64; 3] { [1.0 + p[1], 0.0, 0.0] }
fn tree() -> Octree3 {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    tree.refined(&[*tree.leaves().iter().next().unwrap()], || ControlFlow::Continue(())).unwrap()
}
fn build(tree: &Octree3) -> AdaptiveSolveSpace3 {
    let mut p = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut p).unwrap();
    let op = AdaptiveElasticity3::build_with_embedded_dirichlet(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(), tree, &Slab,
        &IsotropicElastic::new(1.0, 0.0, 1.0).unwrap(), &|_| false, &|_, _| true,
        ElasticityOptions3::default(), Default::default(), Default::default(), &mut q).unwrap();
    AdaptiveSolveSpace3::jacobi(op, 100_000_000)
}
fn candidate() -> CutDensityStudy3<AdaptiveSolveSpace3> {
    CutDensityStudy3::new(build(&Octree3::uniform(2, 4, 4096).unwrap()), 0.15,
        SimpParams { beta: 2.0, ..Default::default() })
}
fn reaction() -> ReactionTarget3<'static> {
    ReactionTarget3 { mode: &right, target: 0.001, scale: 0.002, weight: 1.0 }
}
fn target() -> ReferenceResponseTarget3<'static> {
    ReferenceResponseTarget3 { observation: ReferenceLoad3::body(&observe),
        target: 0.004, scale: 0.02, weight: 0.3 }
}
fn fields_equal(a: &ResponseEvaluation3, b: &ResponseEvaluation3) {
    assert_eq!(a.rho, b.rho);
    assert_eq!(a.projected_rho, b.projected_rho);
    assert_eq!(a.objective.to_bits(), b.objective.to_bits());
    assert_eq!(a.gradient, b.gradient);
    assert_eq!(a.volume_fraction.to_bits(), b.volume_fraction.to_bits());
    assert_eq!(a.volume_gradient, b.volume_gradient);
    assert_eq!(a.displacements, b.displacements);
    assert_eq!(a.adjoints, b.adjoints);
    assert_eq!(a.responses, b.responses);
    assert_eq!(a.reaction_responses, b.reaction_responses);
    assert_eq!(a.external_work, b.external_work);
    assert_eq!(a.augmented_rhs_work, b.augmented_rhs_work);
}
// Independently bind reference observations on the returned candidate and call
// the existing physical evaluator, not the new reference fitting entry point.
fn check_accepted(study: &mut CutDensityStudy3<AdaptiveSolveSpace3>, accepted: &ResponseEvaluation3,
    displacement: bool, options: ProjectedResponseOptions3) {
    let op = study.operator().elasticity();
    let before = op.scales().to_vec();
    let f = op.reference_load(ReferenceLoad3::default(), || ControlFlow::Continue(())).unwrap();
    let q = op.reference_load(target().observation, || ControlFlow::Continue(())).unwrap();
    let targets = [ResponseTarget3 { q: &q, target: target().target, scale: target().scale, weight: target().weight }];
    let case = [ResponseCase3 { force: &f, prescribed: Some(&motion), targets: if displacement { &targets } else { &[] } }];
    let r = [reaction()];
    let mut p = |_| ControlFlow::Continue(());
    let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let again = study.evaluate_responses_with_reactions(&accepted.rho, &case, &[&r], options.response, &mut c).unwrap();
    fields_equal(accepted, &again);
    assert_eq!(study.operator().elasticity().scales(), before);
}

#[test]
fn g3_reference_mixed_fit_is_the_same_optimizer_and_physics_as_explicit_nodal_fitting() {
    let (mut study, _, _) = fixture();
    let (mut direct, _, _) = fixture();
    let f = direct.operator().elasticity().reference_load(ReferenceLoad3::default(), || ControlFlow::Continue(())).unwrap();
    let q = direct.operator().elasticity().reference_load(target().observation, || ControlFlow::Continue(())).unwrap();
    let targets = [target()];
    let cases = [ReferenceResponseCase3 { load: ReferenceLoad3::default(), prescribed: Some(&motion), targets: &targets }; 2];
    let r = [reaction()]; let rows: [&[ReactionTarget3<'_>]; 2] = [&r, &r];
    let options = ProjectedResponseOptions3 { volume_cap: 0.6, ..Default::default() };
    let rho = vec![0.5; study.cells()];
    let mut pa = |_| ControlFlow::Continue(()); let mut ca = SolveControl::new(SolveBudget::default(), &mut pa);
    let fit = study.fit_reference_responses_with_reactions(&cases, &rows, &rho, options, 4, &mut ca).unwrap();
    assert!(fit.outcome.is_ok());
    let nodal_targets = [ResponseTarget3 { q: &q, target: targets[0].target, scale: targets[0].scale, weight: targets[0].weight }];
    let nodal = [ResponseCase3 { force: &f, prescribed: Some(&motion), targets: &nodal_targets }; 2];
    let mut pb = |_| ControlFlow::Continue(()); let mut cb = SolveControl::new(SolveBudget::default(), &mut pb);
    let mut session = ProjectedResponseStudy3::new_with_reactions(&mut direct, &nodal, &rows, &rho, options, &mut cb).unwrap();
    fields_equal(&fit.initial, session.accepted());
    session.run(4).unwrap();
    fields_equal(&fit.accepted, session.accepted());
    assert_eq!(fit.history, session.history());
    assert_eq!(fit.work, session.work());
    assert!(fit.history.len() > 1);
}

#[test]
fn g1_estimate_mark_and_refit_preserves_reaction_targets_and_rebuilds_the_fine_baseline() {
    let (mut source, _, _) = fixture();
    let targets = [target()];
    let cases = [ReferenceResponseCase3 { load: ReferenceLoad3::default(), prescribed: Some(&motion), targets: &targets }];
    let r = [reaction()]; let rows: [&[ReactionTarget3<'_>]; 1] = [&r];
    let options = ProjectedResponseOptions3 { volume_cap: 0.6, ..Default::default() };
    let rho: Vec<_> = (0..source.cells()).map(|i| 0.25 + 0.025 * i as f64).collect();
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let coarse = source.fit_reference_responses_with_reactions(&cases, &rows, &rho, options, 0, &mut c).unwrap();
    let before = source.operator().elasticity().scales().to_vec();
    let evidence = source.estimate_response_enrichment_with_reactions(
        build(&Octree3::uniform(2, 4, 4096).unwrap()), &coarse.accepted, &cases, &rows, Default::default(), &mut c).unwrap();
    let marks = evidence.mark(0.5, 2, || ControlFlow::Continue(())).unwrap();
    assert!(!marks.marked.is_empty());
    let fine_tree = tree().refined(&marks.marked, || ControlFlow::Continue(())).unwrap();
    let next = CutDensityStudy3::new(build(&fine_tree), 0.15, SimpParams { beta: 2.0, ..Default::default() });
    let transfer = AdaptiveTransfer3::new(source.operator().elasticity(), next.operator().elasticity(),
        2_000_000, || ControlFlow::Continue(())).unwrap();
    let inherited: Vec<_> = transfer.parents().iter().map(|&parent| coarse.accepted.rho[parent]).collect();
    drop(transfer);
    let work = c.work();
    let mut result = source.refit_reference_responses_with_reactions(next, &coarse.accepted,
        &cases, &rows, options, 6, 2_000_000, &mut c).unwrap();
    assert!(result.fit.outcome.is_ok());
    assert_eq!(result.inherited_rho, inherited);
    assert_eq!(result.fit.initial.rho, result.starting_rho);
    assert_eq!(result.fit.history[0].iteration, 0);
    assert!(result.study.cells() > source.cells());
    assert!(result.fit.initial.volume_fraction <= options.volume_cap);
    assert!(result.fit.history.len() > 1);
    assert!(c.work().linear_iterations > work.linear_iterations);
    assert_eq!(result.fit.work, c.work());
    assert_eq!(source.operator().elasticity().scales(), before);
    assert_eq!(coarse.accepted.rho, rho);
    check_accepted(&mut result.study, &result.fit.accepted, true, options);
    assert_eq!(source.operator().elasticity().scales(), before);
}

#[test]
fn g1_reaction_only_refit_restores_projected_volume_and_does_not_drop_the_force_objective() {
    let (mut source, _, _) = fixture();
    let cases = [ReferenceResponseCase3 { load: ReferenceLoad3::default(), prescribed: Some(&motion), targets: &[] }];
    let r = [reaction()]; let rows: [&[ReactionTarget3<'_>]; 1] = [&r];
    let options = ProjectedResponseOptions3 { volume_cap: 0.6, ..Default::default() };
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let rho = vec![0.8; source.cells()];
    let coarse = source.fit_reference_responses_with_reactions(&cases, &rows, &rho, options, 0, &mut c).unwrap();
    let before = source.operator().elasticity().scales().to_vec();
    let mut result = source.refit_reference_responses_with_reactions(candidate(), &coarse.accepted,
        &cases, &rows, options, 8, 2_000_000, &mut c).unwrap();
    assert!(result.fit.outcome.is_ok());
    assert!(result.inherited_rho.iter().all(|v| *v == 0.8));
    assert_ne!(result.inherited_rho, result.starting_rho);
    assert!(result.fit.initial.volume_fraction <= options.volume_cap);
    assert!(result.fit.accepted.responses[0].is_empty());
    assert_eq!(result.fit.accepted.reaction_responses[0].len(), 1);
    assert!(result.fit.history.iter().any(|row| row.objective < result.fit.initial.objective
        && row.constraint_violation <= options.optimizer.tolerance));
    check_accepted(&mut result.study, &result.fit.accepted, false, options);
    assert_eq!(source.operator().elasticity().scales(), before);
}

#[test]
fn g0_whole_reaction_family_is_admitted_before_reference_callbacks_or_solves() {
    let (mut study, _, _) = fixture();
    let called = Cell::new(0usize);
    let force = |_| { called.set(called.get() + 1); [0.0; 3] };
    let cases = [ReferenceResponseCase3 { load: ReferenceLoad3::body(&force), prescribed: Some(&motion), targets: &[] }; 2];
    let good = [reaction()]; let bad = [ReactionTarget3 { scale: f64::NAN, ..reaction() }];
    let rho = vec![0.5; study.cells()];
    let before = study.operator().elasticity().scales().to_vec();
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    for rows in [&[&good[..]][..], &[&good[..], &bad[..]][..]] {
        let work = c.work();
        assert!(study.fit_reference_responses_with_reactions(&cases, rows, &rho, Default::default(), 1, &mut c).is_err());
        assert_eq!(c.work(), work);
    }
    let options = ProjectedResponseOptions3 { response: ResponseOptions3 { max_observations: 1, ..Default::default() }, ..Default::default() };
    let work = c.work();
    assert!(study.fit_reference_responses_with_reactions(&cases, &[&good, &good], &rho, options, 1, &mut c).is_err());
    assert_eq!(called.get(), 0); assert_eq!(c.work(), work);
    assert_eq!(study.operator().elasticity().scales(), before);
}

#[test]
fn g4_refit_failures_preserve_source_and_late_cancellation_returns_aligned_fine_prefix() {
    let (mut source, _, _) = fixture();
    let cases = [ReferenceResponseCase3 { load: ReferenceLoad3::default(), prescribed: Some(&motion), targets: &[] }];
    let r = [reaction()]; let rows: [&[ReactionTarget3<'_>]; 1] = [&r];
    let options = ProjectedResponseOptions3 { volume_cap: 0.6, ..Default::default() };
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let coarse = source.fit_reference_responses_with_reactions(&cases, &rows, &vec![0.5; source.cells()], options, 0, &mut c).unwrap();
    let before = source.operator().elasticity().scales().to_vec();
    for stage in ["response-refit-transfer", "response-refit-baseline", "sdf3-response-reaction"] {
        let mut p = |s: SolveProgress| if s.stage == stage { ControlFlow::Break(()) } else { ControlFlow::Continue(()) };
        let mut c = SolveControl::new(SolveBudget::default(), &mut p);
        assert!(source.refit_reference_responses_with_reactions(candidate(), &coarse.accepted,
            &cases, &rows, options, 8, 2_000_000, &mut c).is_err());
        if stage == "sdf3-response-reaction" { assert!(c.work().linear_iterations > 0); }
        assert_eq!(source.operator().elasticity().scales(), before);
    }
    let reached = Cell::new(false);
    let mut p = |s: SolveProgress| if s.stage == "response-projected-accepted" {
        reached.set(true); ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let mut result = source.refit_reference_responses_with_reactions(candidate(), &coarse.accepted,
        &cases, &rows, options, 8, 2_000_000, &mut c).unwrap();
    assert!(reached.get()); assert!(result.fit.outcome.is_err()); assert!(result.fit.history.len() > 1);
    assert_eq!(result.fit.work, c.work());
    check_accepted(&mut result.study, &result.fit.accepted, false, options);
    assert_eq!(source.operator().elasticity().scales(), before);
    let mut p = |_| ControlFlow::Continue(());
    let mut c = SolveControl::new(SolveBudget { total_iterations: 1, ..Default::default() }, &mut p);
    assert!(source.refit_reference_responses_with_reactions(candidate(), &coarse.accepted,
        &cases, &rows, options, 8, 2_000_000, &mut c).is_err());
    assert_eq!(c.work().linear_iterations, 1);
    assert_eq!(source.operator().elasticity().scales(), before);
}

#[test]
fn g0_legacy_refit_cannot_silently_discard_reaction_observations() {
    let (mut source, _, _) = fixture();
    let targets = [target()];
    let cases = [ReferenceResponseCase3 { load: ReferenceLoad3::default(), prescribed: Some(&motion), targets: &targets }];
    let r = [reaction()]; let rows: [&[ReactionTarget3<'_>]; 1] = [&r];
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let fit = source.fit_reference_responses_with_reactions(&cases, &rows, &vec![0.5; source.cells()], Default::default(), 0, &mut c).unwrap();
    let before = source.operator().elasticity().scales().to_vec(); let work = c.work();
    assert!(source.refit_reference_responses(candidate(), &fit.accepted, &cases,
        Default::default(), 1, 2_000_000, &mut c).is_err());
    assert_eq!(c.work(), work);
    assert_eq!(source.operator().elasticity().scales(), before);
    let mut bad = fit.accepted.clone(); bad.projected_rho[0] += 0.01;
    assert!(source.refit_reference_responses_with_reactions(candidate(), &bad, &cases, &rows,
        Default::default(), 1, 2_000_000, &mut c).is_err());
    assert_eq!(source.operator().elasticity().scales(), before);
}

#[test]
fn g5_empty_reaction_reference_rows_preserve_fit_and_refit_arithmetic_and_work() {
    let (mut a, _, _) = fixture(); let (mut b, _, _) = fixture();
    let targets = [target()];
    let cases = [ReferenceResponseCase3 { load: ReferenceLoad3::default(), prescribed: Some(&motion), targets: &targets }];
    let rho = vec![0.5; a.cells()]; let options = ProjectedResponseOptions3::default();
    let mut pa = |_| ControlFlow::Continue(()); let mut ca = SolveControl::new(SolveBudget::default(), &mut pa);
    let mut pb = |_| ControlFlow::Continue(()); let mut cb = SolveControl::new(SolveBudget::default(), &mut pb);
    let fit_a = a.fit_reference_responses(&cases, &rho, options, 2, &mut ca).unwrap();
    let fit_b = b.fit_reference_responses_with_reactions(&cases, &[&[]], &rho, options, 2, &mut cb).unwrap();
    fields_equal(&fit_a.accepted, &fit_b.accepted);
    assert_eq!(fit_a.history, fit_b.history); assert_eq!(fit_a.work, fit_b.work);
    let next_a = a.refit_reference_responses(candidate(), &fit_a.accepted, &cases, options, 2, 2_000_000, &mut ca).unwrap();
    let next_b = b.refit_reference_responses_with_reactions(candidate(), &fit_b.accepted, &cases,
        &[&[]], options, 2, 2_000_000, &mut cb).unwrap();
    fields_equal(&next_a.fit.initial, &next_b.fit.initial);
    fields_equal(&next_a.fit.accepted, &next_b.fit.accepted);
    assert_eq!(next_a.fit.history, next_b.fit.history); assert_eq!(next_a.fit.work, next_b.fit.work);
}
