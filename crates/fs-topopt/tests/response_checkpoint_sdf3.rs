//! Fresh-process-style response recovery on actual CutFEM equilibrium/adjoints.
use std::ops::ControlFlow;
use fs_ascent::projected_al::{ProjectedAlError, ProjectedAlOptions};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_cutfem::elastic3::{adaptive::AdaptiveElasticity3, ElasticityOptions3};
use fs_cutfem::elastic3::adaptive::enrichment::precondition::AdaptiveSolveSpace3;
use fs_cutfem::octree3::Octree3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_topopt::{SimpParams, SolveBudget, SolveControl, SolveProgress};
use fs_topopt::sdf3::CutDensityStudy3;
use fs_topopt::sdf3::response::{
    ProjectedResponseStudy3, ProjectedResponseOptions3, ResponseCase3, ResponseTarget3,
    ReactionTarget3,
};

struct Slab;
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 { (p[0] - 0.17) * (p[0] - 0.83) }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        let x = Interval::new(lo[0], hi[0]);
        (x - Interval::new(0.17, 0.17)) * (x - Interval::new(0.83, 0.83))
    }
    fn derivative_enclose(&self, lo: [f64; 3], hi: [f64; 3], axis: HeightAxis) -> Interval {
        if axis == HeightAxis::X {
            Interval::new(2.0, 2.0) * Interval::new(lo[0], hi[0]) - Interval::new(1.0, 1.0)
        } else { Interval::new(0.0, 0.0) }
    }
}
fn motion(p: [f64; 3], n: [f64; 3]) -> [f64; 3] {
    if n[0] > 0.0 { [0.02, 0.0, 0.0] } else { [0.0, 0.0, 0.007 * p[1]] }
}
fn fixture() -> (CutDensityStudy3<AdaptiveSolveSpace3>, Vec<f64>, Vec<f64>) {
    let mut poll = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(QuadratureOptions3::default(), &mut poll).unwrap();
    let op = AdaptiveElasticity3::build_with_embedded_dirichlet(
        HexCell::try_new([0.0; 3], [1.0; 3]).unwrap(), &Octree3::uniform(1, 4, 4096).unwrap(),
        &Slab, &IsotropicElastic::new(1.0, 0.3, 1.0).unwrap(), &|_| false, &|_, _| true,
        ElasticityOptions3::default(), Default::default(), Default::default(), &mut q,
    ).unwrap();
    let f = op.body_load(&|_| [0.002, 0.0, -0.003], || ControlFlow::Continue(())).unwrap();
    let q = op.body_load(&|p| [0.0, 0.0, 1.0 + p[0]], || ControlFlow::Continue(())).unwrap();
    (CutDensityStudy3::new(AdaptiveSolveSpace3::jacobi(op, 100_000_000), 0.15, SimpParams::default()), f, q)
}
fn options() -> ProjectedResponseOptions3 {
    ProjectedResponseOptions3 {
        density_floor: 0.1, objective_scale: 2.0,
        optimizer: ProjectedAlOptions { tolerance: 1e-9, ..Default::default() },
        ..Default::default()
    }
}

#[test]
fn restoration_preserves_the_search_and_rebuilds_independent_motion_fields() {
    let (mut whole, f, q) = fixture();
    let (mut split, _, _) = fixture();
    let other_force: Vec<_> = f.iter().map(|f| -0.5 * f).collect();
    let targets = [ResponseTarget3 { q: &q, target: 0.003, scale: 0.02, weight: 1.0 }];
    let cases = [
        ResponseCase3 { force: &f, prescribed: Some(&motion), targets: &targets },
        ResponseCase3 { force: &other_force, prescribed: Some(&motion), targets: &targets },
    ];
    let rho = vec![0.5; whole.cells()];
    let mut pa = |_| ControlFlow::Continue(());
    let mut ca = SolveControl::new(SolveBudget::default(), &mut pa);
    let mut pb = |_| ControlFlow::Continue(());
    let mut cb = SolveControl::new(SolveBudget::default(), &mut pb);
    let mut a = ProjectedResponseStudy3::new(&mut whole, &cases, &rho, options(), &mut ca).unwrap();
    let ra = a.run(6).unwrap();
    let mut b = ProjectedResponseStudy3::new(&mut split, &cases, &rho, options(), &mut cb).unwrap();
    b.run(2).unwrap();
    assert_eq!(b.optimizer_work().iterations, 2, "require real accepted updates before restart");
    let checkpoint = b.checkpoint();
    let cost = checkpoint.restoration_cost();
    let accepted = b.accepted().clone();
    drop(b);
    let (mut rebuilt, _, _) = fixture();
    let mut b = ProjectedResponseStudy3::restore(&mut rebuilt, &cases, checkpoint, options(), &mut cb).unwrap();
    assert_eq!(b.accepted().displacements, accepted.displacements);
    assert_eq!(b.accepted().adjoints, accepted.adjoints);
    assert_eq!(b.accepted().gradient, accepted.gradient);
    let before = b.work();
    b.run(0).unwrap();
    assert_eq!(b.work(), before, "zero-step inspection cannot repeat restoration physics");
    let rb = b.run(4).unwrap();
    assert_eq!(ra.stop, rb.stop);
    assert_eq!((ra.multiplier, ra.penalty), (rb.multiplier, rb.penalty));
    assert_eq!(a.point(), b.point());
    assert_eq!(a.history(), b.history());
    assert_eq!(a.accepted().displacements, b.accepted().displacements);
    assert_eq!(a.accepted().adjoints, b.accepted().adjoints);
    assert_ne!(b.accepted().displacements[0], b.accepted().displacements[1]);
    assert_eq!(b.evaluations(), a.evaluations() + cost);
    assert!(b.work().linear_iterations > a.work().linear_iterations);
    assert_eq!(b.checkpoint().restoration_evaluations, cost);
}

#[test]
fn an_infeasible_accepted_step_keeps_and_restores_a_distinct_feasible_incumbent() {
    let (mut study, f, _) = fixture();
    // A squared compliance observation prefers more stiffness. A weak initial
    // AL penalty deliberately permits a real volume-infeasible accepted step.
    let targets = [ResponseTarget3 { q: &f, target: 0.0, scale: 1e-6, weight: 1.0 }];
    let cases = [ResponseCase3 { force: &f, prescribed: None, targets: &targets }];
    let mut opts = options();
    opts.volume_cap = 0.4001;
    opts.optimizer.initial_penalty = 1e-8;
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let rho = vec![0.4; study.cells()];
    let mut design = ProjectedResponseStudy3::new(&mut study, &cases, &rho, opts, &mut control).unwrap();
    let baseline = design.best_feasible().unwrap().clone();
    design.run(1).unwrap();
    assert_eq!(design.optimizer_work().iterations, 1);
    assert!(design.constraint_violation() > opts.optimizer.tolerance);
    assert_eq!(design.best_feasible().unwrap().rho, baseline.rho);
    assert!(design.accepted().objective < baseline.objective);
    let checkpoint = design.checkpoint();
    assert_eq!(checkpoint.restoration_cost(), 2);
    let accepted = design.accepted().clone();
    drop(design);
    let (mut other, _, _) = fixture();
    let restored = ProjectedResponseStudy3::restore(&mut other, &cases, checkpoint, opts, &mut control).unwrap();
    assert_eq!(restored.accepted().displacements, accepted.displacements);
    assert_eq!(restored.best_feasible().unwrap().displacements, baseline.displacements);
    assert_eq!(restored.best_feasible().unwrap().gradient, baseline.gradient);
    assert_ne!(restored.best_feasible().unwrap().rho, restored.point());
    assert_eq!(restored.checkpoint().restoration_evaluations, 2);
}

#[test]
fn mutation_exhaustion_and_cancelled_readmission_cannot_install_unverified_scales() {
    let (mut study, f, q) = fixture();
    let targets = [ResponseTarget3 { q: &q, target: 0.003, scale: 0.02, weight: 1.0 }];
    let cases = [ResponseCase3 { force: &f, prescribed: Some(&motion), targets: &targets }];
    let mut opts = options();
    opts.volume_cap = 1.0;
    let mut poll = |_| ControlFlow::Continue(());
    let mut control = SolveControl::new(SolveBudget::default(), &mut poll);
    let rho = vec![0.5; study.cells()];
    let mut design = ProjectedResponseStudy3::new(&mut study, &cases, &rho, opts, &mut control).unwrap();
    design.run(2).unwrap();
    let saved = design.checkpoint();
    let original = format!("{saved:?}");
    drop(design);
    for mutation in 0..5 {
        let (mut rebuilt, _, _) = fixture();
        let incoming = rebuilt.operator().elasticity().scales().to_vec();
        let mut damaged = saved.clone();
        match mutation {
            0 => damaged.optimizer.sample.gradient[0] += 0.125,
            1 => damaged.history.last_mut().unwrap().objective += 1.0,
            2 => damaged.best_feasible_density = None,
            3 => damaged.history[0].iteration = 1,
            _ => damaged.restoration_evaluations = usize::MAX,
        }
        let mut poll = |_| ControlFlow::Continue(());
        let mut c = SolveControl::new(SolveBudget::default(), &mut poll);
        assert!(ProjectedResponseStudy3::restore(&mut rebuilt, &cases, damaged, opts, &mut c).is_err());
        assert_eq!(rebuilt.operator().elasticity().scales(), incoming);
    }
    let (mut rebuilt, _, _) = fixture();
    let incoming = rebuilt.operator().elasticity().scales().to_vec();
    let mut limited = opts;
    limited.optimizer.max_evaluations = saved.optimizer.work.evaluations + saved.restoration_cost() - 1;
    let mut poll = |_| ControlFlow::Continue(());
    let mut c = SolveControl::new(SolveBudget::default(), &mut poll);
    assert!(ProjectedResponseStudy3::restore(&mut rebuilt, &cases, saved.clone(), limited, &mut c).is_err());
    assert_eq!(c.work().linear_solves, 0, "reserve the complete restoration before solving");
    let mut poll = |p: SolveProgress| if p.stage == "response-projected-restore-publish" {
        ControlFlow::Break(())
    } else { ControlFlow::Continue(()) };
    let mut c = SolveControl::new(SolveBudget::default(), &mut poll);
    assert!(ProjectedResponseStudy3::restore(&mut rebuilt, &cases, saved.clone(), opts, &mut c).is_err());
    assert!(c.work().linear_solves > 0, "cancel after actual physics, not before it");
    assert_eq!(rebuilt.operator().elasticity().scales(), incoming);
    assert_eq!(format!("{saved:?}"), original);
}

#[test]
fn mixed_reaction_recovery_keeps_direct_material_derivatives_and_original_work() {
    let (mut study, f, q) = fixture();
    let targets = [ResponseTarget3 { q: &q, target: 0.003, scale: 0.02, weight: 0.7 }];
    let cases = [ResponseCase3 { force: &f, prescribed: Some(&motion), targets: &targets }];
    let mode = |_: [f64; 3], n: [f64; 3]| if n[0] > 0.0 { [1.0, 0.0, 0.0] } else { [0.0; 3] };
    let targets = [ReactionTarget3 { mode: &mode, target: 0.01, scale: 0.1, weight: 0.3 }];
    let reactions: [&[ReactionTarget3<'_>]; 1] = [&targets];
    let mut poll = |_| ControlFlow::Continue(());
    let mut c = SolveControl::new(SolveBudget::default(), &mut poll);
    let rho = vec![0.5; study.cells()];
    let mut design = ProjectedResponseStudy3::new_with_reactions(&mut study, &cases, &reactions,
        &rho, options(), &mut c).unwrap();
    design.run(2).unwrap();
    assert!(design.optimizer_work().iterations > 0);
    let saved = design.checkpoint();
    let accepted = design.accepted().clone();
    drop(design);
    let (mut rebuilt, _, _) = fixture();
    let restored = ProjectedResponseStudy3::restore_with_reactions(&mut rebuilt, &cases, &reactions,
        saved, options(), &mut c).unwrap();
    assert_eq!(restored.accepted().reaction_responses, accepted.reaction_responses);
    assert_eq!(restored.accepted().gradient, accepted.gradient);
    assert_eq!(restored.accepted().adjoints, accepted.adjoints);
    drop(restored);
    let changed = [ReactionTarget3 { target: 0.02, ..targets[0] }];
    let changed_family: [&[ReactionTarget3<'_>]; 1] = [&changed];
    // A displacement-only reconstruction must not admit a reaction checkpoint.
    let mut original = ProjectedResponseStudy3::new_with_reactions(&mut rebuilt, &cases, &reactions,
        &rho, options(), &mut c).unwrap();
    original.run(1).unwrap();
    let checkpoint = original.checkpoint();
    drop(original);
    assert!(matches!(ProjectedResponseStudy3::restore_with_reactions(&mut rebuilt, &cases, &changed_family,
        checkpoint, options(), &mut c), Err(ProjectedAlError::Invalid(_))));
}
