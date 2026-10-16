//! G0/G1/G3/G4/G5: real cut-cell mixed-response adaptive campaigns.
use std::{cell::Cell, ops::ControlFlow};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{adaptive::AdaptiveElasticity3, ElasticityOptions3};
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::elastic3::surface::ReferenceLoad3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{EvaluationStop, SimpParams, SolveBudget, SolveControl, SolveProgress};
use fs_topopt::sdf3::{CutDensityStudy3, PhysicalRegion3};
use fs_topopt::sdf3::response::{ReactionTarget3, ResponseEvaluation3};
use fs_topopt::sdf3::response::adaptive::*;
use fs_topopt::sdf3::response::refinement::{ReferenceResponseCase3, ReferenceResponseTarget3, ResponseRefinementOptions3};

type Study = CutDensityStudy3<AdaptiveSolveSpace3>;
struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 { (p[0]-0.17)*(p[0]-0.83) }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        let x = Interval::new(lo[0], hi[0]);
        (x-Interval::new(0.17, 0.17))*(x-Interval::new(0.83, 0.83))
    }
    fn derivative_enclose(&self, lo: [f64; 3], hi: [f64; 3], axis: HeightAxis) -> Interval {
        if axis == HeightAxis::X { Interval::new(2.0, 2.0)*Interval::new(lo[0], hi[0])-Interval::new(1.0, 1.0) }
        else { Interval::new(0.0, 0.0) }
    }
}
fn motion(p: [f64; 3], _: [f64; 3]) -> [f64; 3] { [0.02*(p[0]-0.17), 0.0, 0.0] }
fn right(_: [f64; 3], n: [f64; 3]) -> [f64; 3] { if n[0] > 0.0 { [1.0, 0.0, 0.0] } else { [0.0; 3] } }
fn observe(p: [f64; 3]) -> [f64; 3] { [1.0+p[1], 0.0, 0.0] }
fn tree() -> Octree3 {
    // Initial mixed leaves reach level two. Reserve level five for the final
    // probe after two selectively promoted grids, not just their fit geometry.
    let t = Octree3::uniform(1, 5, 4096).unwrap();
    t.refined(&[*t.leaves().iter().next().unwrap()], || ControlFlow::Continue(())).unwrap()
}
fn build(t: &Octree3) -> Result<AdaptiveSolveSpace3, &'static str> {
    let mut p = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut p).unwrap();
    let op = AdaptiveElasticity3::build_with_embedded_dirichlet(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(), t, &Slab,
        &IsotropicElastic::new(1.0, 0.0, 1.0).unwrap(), &|_| false, &|_, _| true,
        ElasticityOptions3::default(), Default::default(), Default::default(), &mut q).unwrap();
    Ok(AdaptiveSolveSpace3::jacobi(op, 100_000_000))
}
fn fixture() -> Study { CutDensityStudy3::new(build(&tree()).unwrap(), 0.15, SimpParams { beta: 2.0, ..Default::default() }) }
fn target() -> ReferenceResponseTarget3<'static> {
    ReferenceResponseTarget3 { observation: ReferenceLoad3::body(&observe), target: 0.004, scale: 0.02, weight: 0.3 }
}
fn reaction() -> ReactionTarget3<'static> { ReactionTarget3 { mode: &right, target: 0.001, scale: 0.002, weight: 1.0 } }
fn cases<'a>(targets: &'a [ReferenceResponseTarget3<'a>]) -> [ReferenceResponseCase3<'a>; 1] {
    [ReferenceResponseCase3 { load: ReferenceLoad3::default(), prescribed: Some(&motion), targets }]
}
fn rho(s: &Study) -> Vec<f64> { (0..s.cells()).map(|i| 0.25+0.025*i as f64).collect() }
fn options() -> AdaptiveResponseOptions3 {
    let mut o = AdaptiveResponseOptions3 { steps_per_grid: 0, max_refinements: 2,
        max_marks: 2, objective_tolerance: 0.0, response_tolerance: 0.0, volume_tolerance: 0.0,
        ..Default::default() };
    o.fit.volume_cap = 0.6;
    o
}
fn equal(a: &ResponseEvaluation3, b: &ResponseEvaluation3) {
    assert_eq!(a.rho, b.rho); assert_eq!(a.projected_rho, b.projected_rho);
    assert_eq!(a.objective.to_bits(), b.objective.to_bits());
    assert_eq!(a.responses, b.responses); assert_eq!(a.reaction_responses, b.reaction_responses);
    assert_eq!(a.displacements, b.displacements); assert_eq!(a.adjoints, b.adjoints);
    assert_eq!(a.gradient, b.gradient); assert_eq!(a.volume_gradient, b.volume_gradient);
}
fn aligned(s: &AdaptiveResponseStudy3<'_, AdaptiveSolveSpace3>) {
    let p = s.study().params();
    let scales: Vec<_> = s.accepted().projected_rho.iter().map(|r| {
        p.e_min+(1.0-p.e_min)*if *r == 0.0 { 0.0 } else { fs_math::det::pow(*r, p.penal) }
    }).collect();
    assert_eq!(s.study().operator().elasticity().scales(), scales);
    assert!(s.study().operator().elasticity().leaves().iter().all(|l| s.tree().leaves().contains(l)));
    assert!(s.accepted().volume_fraction <= 0.6+1e-6);
}

#[test]
fn g1_automatic_mixed_refit_matches_manual_physics_and_marks() {
    let mut a = fixture(); let mut b = fixture(); let raw = rho(&a);
    let targets = [target()]; let cases = cases(&targets); let rt = [reaction()]; let rows: [&[_]; 1] = [&rt];
    let o = options();
    let mut pa = |_| ControlFlow::Continue(()); let mut ca = SolveControl::new(SolveBudget::default(), &mut pa);
    let mut session = AdaptiveResponseStudy3::new(&mut a, tree(), &cases, Some(&rows), &raw, o, &mut ca).unwrap();
    assert_eq!(session.run(1, |t, _| build(t), &mut ca).unwrap(), AdaptiveResponseStop3::Paused);
    assert_eq!(session.rounds().len(), 1); aligned(&session);
    let mut pb = |_| ControlFlow::Continue(()); let mut cb = SolveControl::new(SolveBudget::default(), &mut pb);
    let fit = b.fit_reference_responses_with_reactions(&cases, &rows, &raw, o.fit, 0, &mut cb).unwrap();
    let source = fit.best_feasible.as_ref().unwrap();
    let probe = tree().refined(b.operator().elasticity().leaves(), || ControlFlow::Continue(())).unwrap();
    let evidence = b.estimate_response_enrichment_with_reactions(build(&probe).unwrap(), source, &cases, &rows,
        ResponseRefinementOptions3 { response: o.fit.response, ..Default::default() }, &mut cb).unwrap();
    let marks = evidence.mark(o.theta, o.max_marks, || ControlFlow::Continue(())).unwrap();
    assert_eq!(session.rounds()[0].marking.marked, marks.marked);
    let fine = tree().refined(&marks.marked, || ControlFlow::Continue(())).unwrap();
    let next = b.refit_reference_responses_with_reactions(CutDensityStudy3::new(build(&fine).unwrap(), 0.15, b.params()),
        source, &cases, &rows, o.fit, 0, o.max_transfer_terms, &mut cb).unwrap();
    equal(session.accepted(), next.fit.best_feasible.as_ref().unwrap());
    equal(&session.fit().initial, &next.fit.initial);
    assert_eq!(session.tree().leaves(), fine.leaves());
    assert!(session.probe().is_some(), "final evidence must describe the promoted grid");
    let loss = session.accepted().objective;
    assert!((session.probe().unwrap().coarse_objective-loss).abs() <= 1e-9*loss.abs().max(1e-12));
}

#[test]
fn g5_segmented_campaign_reuses_the_probe_and_preserves_physical_work() {
    let mut a = fixture(); let mut b = fixture(); let raw = rho(&a);
    let targets = [target()]; let cases = cases(&targets); let rt = [reaction()]; let rows: [&[_]; 1] = [&rt];
    let o = options();
    let mut pa = |_| ControlFlow::Continue(()); let mut ca = SolveControl::new(SolveBudget::default(), &mut pa);
    let mut pb = |_| ControlFlow::Continue(()); let mut cb = SolveControl::new(SolveBudget::default(), &mut pb);
    let mut a = AdaptiveResponseStudy3::new(&mut a, tree(), &cases, Some(&rows), &raw, o, &mut ca).unwrap();
    let mut b = AdaptiveResponseStudy3::new(&mut b, tree(), &cases, Some(&rows), &raw, o, &mut cb).unwrap();
    let mut na = 0; let mut nb = 0;
    let ra = a.run(2, |t, _| { na += 1; build(t) }, &mut ca).unwrap();
    b.run(1, |t, _| { nb += 1; build(t) }, &mut cb).unwrap();
    let rb = b.run(1, |t, _| { nb += 1; build(t) }, &mut cb).unwrap();
    assert_eq!(ra, rb); assert_eq!(na, nb); assert_eq!(a.rounds().len(), 2);
    equal(a.accepted(), b.accepted()); assert_eq!(a.tree().leaves(), b.tree().leaves());
    assert_eq!(a.measured_change().unwrap(), b.measured_change().unwrap());
    assert_eq!(ca.work(), cb.work()); aligned(&a); aligned(&b);
    let spent = cb.work();
    b.run(1, |_, _| -> Result<AdaptiveSolveSpace3, &'static str> { panic!("cached final probe must be reused") }, &mut cb).unwrap();
    assert_eq!(cb.work(), spent);
}

#[test]
fn g0_equal_squared_losses_cannot_hide_opposite_observation_errors() {
    let mut s = fixture(); let raw = rho(&s); let targets = [target()]; let first = cases(&targets);
    let rt = [reaction()]; let rows: [&[_]; 1] = [&rt]; let mut o = options(); o.max_refinements = 0;
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let mut provisional = AdaptiveResponseStudy3::new(&mut s, tree(), &first, Some(&rows), &raw, o, &mut c).unwrap();
    provisional.run(0, |t, _| build(t), &mut c).unwrap();
    let e = &provisional.probe().unwrap().cases[0];
    let targets = [ReferenceResponseTarget3 { target: 0.5*(e.coarse_responses[0]+e.fine_responses[0]), ..target() }];
    let rt = [ReactionTarget3 { target: 0.5*(e.coarse_reactions[0]+e.fine_reactions[0]), ..reaction() }];
    drop(provisional);
    let cases = cases(&targets); let rows: [&[_]; 1] = [&rt];
    o.objective_tolerance = 1e-6; o.response_tolerance = 1e-10; o.volume_tolerance = 1e-6;
    let mut session = AdaptiveResponseStudy3::new(&mut s, tree(), &cases, Some(&rows), &raw, o, &mut c).unwrap();
    let stop = session.run(0, |t, _| build(t), &mut c).unwrap();
    let d = session.measured_change().unwrap().unwrap();
    assert!(d.objective <= o.objective_tolerance);
    assert!(d.displacement > o.response_tolerance || d.reaction > o.response_tolerance);
    assert_eq!(stop, AdaptiveResponseStop3::RefinementLimit);
}

#[test]
fn g4_cancellation_on_either_side_of_promotion_keeps_an_aligned_resumable_source() {
    for stage in ["response-adapt-promote", "response-adapt-accepted"] {
        let mut s = fixture(); let raw = rho(&s); let targets = [target()]; let cases = cases(&targets);
        let rt = [reaction()]; let rows: [&[_]; 1] = [&rt]; let enabled = Cell::new(false);
        let mut p = |s: SolveProgress| if enabled.get() && s.stage == stage { ControlFlow::Break(()) } else { ControlFlow::Continue(()) };
        let mut c = SolveControl::new(SolveBudget::default(), &mut p);
        let mut session = AdaptiveResponseStudy3::new(&mut s, tree(), &cases, Some(&rows), &raw, options(), &mut c).unwrap();
        let before = session.accepted().clone(); let spent = c.work(); enabled.set(true);
        assert!(matches!(session.run(1, |t, _| build(t), &mut c),
            Err(AdaptiveResponseError3::Evaluation(EvaluationStop::Cancelled))));
        enabled.set(false); assert!(c.work().linear_iterations > spent.linear_iterations); aligned(&session);
        if stage == "response-adapt-promote" {
            assert!(session.rounds().is_empty()); equal(session.accepted(), &before);
            assert!(session.unpromoted().is_some());
            session.run(1, |t, _| build(t), &mut c).unwrap();
        } else {
            assert_eq!(session.rounds().len(), 1); assert!(session.probe().is_none());
            session.run(0, |t, _| build(t), &mut c).unwrap();
        }
        assert_eq!(session.rounds().len(), 1); aligned(&session);
    }
}

#[test]
fn g4_factory_failure_and_wrong_background_do_not_replace_the_source() {
    let mut s = fixture(); let raw = rho(&s); let targets = [target()]; let cases = cases(&targets);
    let rt = [reaction()]; let rows: [&[_]; 1] = [&rt];
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let mut session = AdaptiveResponseStudy3::new(&mut s, tree(), &cases, Some(&rows), &raw, options(), &mut c).unwrap();
    let before = session.accepted().clone(); let spent = c.work();
    assert!(matches!(session.run(1, |_, _| Err::<AdaptiveSolveSpace3, _>("geometry budget"), &mut c),
        Err(AdaptiveResponseError3::Build("geometry budget"))));
    assert!(matches!(session.run(1, |_, _| build(&tree()), &mut c), Err(AdaptiveResponseError3::Invalid(_))));
    assert_eq!(c.work(), spent); equal(session.accepted(), &before); aligned(&session);
}

#[test]
fn g1_protected_material_regions_follow_the_exact_parent_map() {
    let s = fixture(); let raw = rho(&s);
    let mut labels = vec![PhysicalRegion3::Design; s.cells()];
    labels[0] = PhysicalRegion3::Solid; labels[1] = PhysicalRegion3::Void;
    let mut s = s.with_physical_regions(labels.clone()).unwrap();
    let targets = [target()]; let cases = cases(&targets); let rt = [reaction()]; let rows: [&[_]; 1] = [&rt];
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let mut session = AdaptiveResponseStudy3::new(&mut s, tree(), &cases, Some(&rows), &raw, options(), &mut c).unwrap();
    session.run(1, |t, _| build(t), &mut c).unwrap(); assert_eq!(session.rounds().len(), 1);
    let old = build(&tree()).unwrap();
    let transfer = AdaptiveTransfer3::new(old.elasticity(), session.study().operator().elasticity(),
        2_000_000, || ControlFlow::Continue(())).unwrap();
    for ((&parent, label), physical) in transfer.parents().iter()
        .zip(session.study().physical_regions().unwrap()).zip(&session.accepted().projected_rho) {
        assert_eq!(*label, labels[parent]);
        match label { PhysicalRegion3::Solid => assert_eq!(*physical, 1.0),
            PhysicalRegion3::Void => assert_eq!(*physical, 0.0), PhysicalRegion3::Design => {} }
    }
    aligned(&session);
}

#[test]
fn g0_policy_and_whole_family_refusals_precede_reference_callbacks() {
    let mut s = fixture(); let raw = rho(&s); let calls = Cell::new(0);
    let observe = |p| { calls.set(calls.get()+1); observe(p) };
    let targets = [ReferenceResponseTarget3 { observation: ReferenceLoad3::body(&observe), ..target() }];
    let cases = cases(&targets); let rt = [reaction()]; let rows: [&[_]; 1] = [&rt];
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let bad = AdaptiveResponseOptions3 { response_tolerance: f64::NAN, ..options() };
    assert!(AdaptiveResponseStudy3::new(&mut s, tree(), &cases, Some(&rows), &raw, bad, &mut c).is_err());
    assert!(AdaptiveResponseStudy3::new(&mut s, tree(), &cases, Some(&[]), &raw, options(), &mut c).is_err());
    assert_eq!(calls.get(), 0); assert_eq!(c.work().linear_solves, 0);
}

#[test]
fn g3_empty_reaction_rows_preserve_displacement_only_fitting_and_estimation() {
    let mut a = fixture(); let mut b = fixture(); let raw = rho(&a); let targets = [target()]; let cases = cases(&targets);
    let mut o = options(); o.max_refinements = 0;
    let rows: [&[ReactionTarget3<'_>]; 1] = [&[]];
    let mut pa = |_| ControlFlow::Continue(()); let mut ca = SolveControl::new(SolveBudget::default(), &mut pa);
    let mut pb = |_| ControlFlow::Continue(()); let mut cb = SolveControl::new(SolveBudget::default(), &mut pb);
    let mut a = AdaptiveResponseStudy3::new(&mut a, tree(), &cases, None, &raw, o, &mut ca).unwrap();
    let mut b = AdaptiveResponseStudy3::new(&mut b, tree(), &cases, Some(&rows), &raw, o, &mut cb).unwrap();
    let ra = a.run(0, |t, _| build(t), &mut ca).unwrap(); let rb = b.run(0, |t, _| build(t), &mut cb).unwrap();
    assert_eq!(ra, rb); equal(a.accepted(), b.accepted());
    assert_eq!(a.measured_change().unwrap(), b.measured_change().unwrap()); assert_eq!(ca.work(), cb.work());
}

#[test]
fn g1_fitted_campaign_retains_the_same_best_feasible_endpoint_as_the_nodal_owner() {
    let mut s = fixture(); let raw = vec![0.5; s.cells()]; let targets = [target()]; let cases = cases(&targets);
    let rt = [reaction()]; let rows: [&[_]; 1] = [&rt];
    let o = AdaptiveResponseOptions3 { steps_per_grid: 8, ..options() };
    let mut p = |_| ControlFlow::Continue(()); let mut c = SolveControl::new(SolveBudget::default(), &mut p);
    let session = AdaptiveResponseStudy3::new(&mut s, tree(), &cases, Some(&rows), &raw, o, &mut c).unwrap();
    assert!(session.fit().outcome.is_ok());
    assert!(session.accepted().objective < session.fit().initial.objective);
    equal(session.accepted(), session.fit().best_feasible.as_ref().unwrap()); aligned(&session);
    let min = session.fit().history.iter().filter(|r| r.constraint_violation <= o.fit.optimizer.tolerance)
        .map(|r| r.objective).reduce(f64::min).unwrap();
    assert_eq!(session.accepted().objective.to_bits(), min.to_bits());
}
