//! Fine-grid refits retain exact accepted optimizer state across interruption.
use super::*;
use fs_ascent::projected_al::ProjectedAlError;
use fs_topopt::sdf3::response::limits::refinement::{
    ReferenceBoundedResponse3, ReferenceResponseLimit3,
};
use fs_topopt::sdf3::robust::envelope::affine_load_space::constraints::{
    AffineLoadSpaceConstraints3, ReferenceAffineLoadSpaceConstraints3,
};
use fs_topopt::sdf3::robust::envelope::affine_load_space::refit::AffineLoadSpaceRefitError3;
use fs_topopt::sdf3::robust::envelope::affine_load_space::response::reference::ReferenceAffineResponseLimits3;
use fs_topopt::sdf3::robust::envelope::affine_load_space::stress::AffineStressLimit3;
use fs_topopt::sdf3::stress::StressOptions3;
use std::cell::Cell;

fn observation(p: [f64; 3]) -> [f64; 3] {
    [0.0, 1.0 + p[0], 0.4]
}

fn limits() -> [ReferenceResponseLimit3<'static>; 2] {
    let response = ReferenceBoundedResponse3::Displacement(ReferenceLoad3::body(&observation));
    [
        ReferenceResponseLimit3 {
            response,
            sense: ResponseLimitSense3::AtMost,
            bound: 1.0,
            scale: 1.0,
        },
        ReferenceResponseLimit3 {
            response,
            sense: ResponseLimitSense3::AtLeast,
            bound: -1.0,
            scale: 1.0,
        },
    ]
}

fn constraints<'a>(
    limits: &'a [ReferenceResponseLimit3<'a>],
) -> ReferenceAffineLoadSpaceConstraints3<'a> {
    ReferenceAffineLoadSpaceConstraints3 {
        responses: Some(ReferenceAffineResponseLimits3 {
            limits,
            options: AffineResponseOptions3::default(),
        }),
        stress: Some(AffineStressLimit3 {
            limit: 100.0,
            options: StressOptions3 {
                relaxation_power: 2.0,
                aggregation_power: 2.0,
                ..Default::default()
            },
        }),
        sharpness: 40.0,
    }
}

fn options() -> AffineLoadSpaceRefitOptions3 {
    let mut design = AffineLoadSpaceDesignOptions3 {
        density_floor: 0.6,
        ..Default::default()
    };
    design.optimizer.max_spectral_step = 0.1;
    AffineLoadSpaceRefitOptions3 {
        design,
        accepted_steps: 2,
        ..Default::default()
    }
}

#[test]
fn refitted_candidate_checkpoint_continues_exactly_on_a_rebuilt_native_fine_owner() {
    let tree = Octree3::uniform(1, 3, 512).unwrap();
    let source = adaptive_material_study(&tree);
    let rho: Vec<_> = (0..source.cells())
        .map(|i| 0.70 + 0.01 * i as f64)
        .collect();
    let fine_tree = tree
        .refined(&[source.operator().elasticity().leaves()[2]], poll)
        .unwrap();
    let fine = adaptive_on_tree(&fine_tree);
    let parents = AdaptiveTransfer3::new(
        source.operator().elasticity(),
        fine.elasticity(),
        options().max_transfer_terms,
        poll,
    )
    .unwrap()
    .parents()
    .to_vec();
    let expected_rho: Vec<_> = parents.iter().map(|&parent| rho[parent]).collect();
    let expected_regions: Vec<_> = parents
        .iter()
        .map(|&parent| source.physical_regions().unwrap()[parent])
        .collect();
    let source_scales = source.operator().elasticity().scales().to_vec();
    let source_fingerprint = source.material_body_load_fingerprint();
    let source_regions = source.physical_regions().unwrap().to_vec();
    let laws = design_laws();
    let dead = ReferenceLoad3::body(&refinement_dead);
    let cap = AffineLoadSpace3::new(laws.len(), 1.0).unwrap();
    let limits = limits();
    let requirements = constraints(&limits);
    let mut reference_poll = |_| poll();
    let mut reference_control = SolveControl::new(SolveBudget::default(), &mut reference_poll);
    let seed = source
        .prepare_affine_load_space_refit_with_constraints(
            fine,
            &rho,
            laws.len(),
            &cap,
            mass_goal(),
            options(),
            requirements,
            &mut reference_control,
        )
        .unwrap();
    assert_eq!(seed.inherited_rho, expected_rho);
    assert_eq!(seed.study.physical_regions().unwrap(), expected_regions);
    assert!(seed.study.can_refine_material_body_load());
    assert_ne!(
        seed.study.material_body_load_fingerprint(),
        source_fingerprint
    );
    assert_eq!(reference_control.work().linear_solves, 0);
    let expected = seed
        .study
        .solve_reference_affine_load_space_with_constraints(
            &seed.inherited_rho,
            dead,
            &laws,
            &cap,
            mass_goal(),
            options(),
            requirements,
            &mut reference_control,
        )
        .unwrap();
    assert!(expected.completed_and_rechecked(), "{:?}", expected.outcome);
    let expected_checkpoint = expected.checkpoint.as_ref().unwrap();
    assert_eq!(expected_checkpoint.optimizer.work.iterations, 2);

    // Stop after the native owner has retained its FIRST accepted update.
    // The error must not throw away its multiplier, spectral step or archive.
    let armed = Cell::new(true);
    let mut split_poll = |progress: SolveProgress| {
        if progress.stage == "sdf3-affine-design-accepted" && armed.replace(false) {
            ControlFlow::Break(())
        } else {
            poll()
        }
    };
    let mut split_control = SolveControl::new(SolveBudget::default(), &mut split_poll);
    let prefix = source
        .refit_reference_affine_load_space_with_constraints(
            adaptive_on_tree(&fine_tree),
            &rho,
            dead,
            &laws,
            &cap,
            mass_goal(),
            options(),
            requirements,
            &mut split_control,
        )
        .unwrap();
    assert!(matches!(
        prefix.outcome,
        Err(ProjectedAlError::Evaluation(EvaluationStop::Cancelled))
    ));
    assert!(!prefix.completed_and_rechecked());
    let mut checkpoint = prefix.checkpoint.as_ref().unwrap().clone();
    assert_eq!(checkpoint.optimizer.work.iterations, 1);
    assert_eq!(checkpoint.history, prefix.history);
    assert_eq!(
        checkpoint.optimizer.point,
        prefix.accepted.as_ref().unwrap().rho
    );
    let initial_audit = format!("{:?}", checkpoint.gradient_check);
    assert_eq!(
        initial_audit,
        format!("{:?}", expected_checkpoint.gradient_check)
    );
    let restore_cost = checkpoint.restoration_cost();
    let before_restore = split_control.work();

    let mut fresh = source
        .prepare_affine_load_space_refit_with_constraints(
            adaptive_on_tree(&fine_tree),
            &rho,
            laws.len(),
            &cap,
            mass_goal(),
            options(),
            requirements,
            &mut split_control,
        )
        .unwrap();
    assert_eq!(fresh.inherited_rho, expected_rho);
    assert_eq!(
        fresh.study.material_body_load_fingerprint(),
        prefix.study.material_body_load_fingerprint()
    );
    let native = fresh.study.operator().elasticity();
    let dead_force = native.reference_load(dead, poll).unwrap();
    let forces: Vec<_> = laws
        .iter()
        .map(|&law| native.reference_load(law, poll).unwrap())
        .collect();
    let basis: Vec<_> = forces.iter().map(Vec::as_slice).collect();
    let q = native
        .reference_load(ReferenceLoad3::body(&observation), poll)
        .unwrap();
    let nodal: Vec<_> = limits
        .iter()
        .map(|limit| ResponseLimit3 {
            response: BoundedResponse3::Displacement(&q),
            sense: limit.sense,
            bound: limit.bound,
            scale: limit.scale,
        })
        .collect();
    let constraints = AffineLoadSpaceConstraints3 {
        responses: Some(AffineResponseLimits3 {
            limits: &nodal,
            options: requirements.responses.unwrap().options,
        }),
        stress: requirements.stress,
        sharpness: requirements.sharpness,
    };
    let mut restored = AffineLoadSpaceDesign3::restore_with_constraints(
        &mut fresh.study,
        &dead_force,
        &basis,
        &mut checkpoint,
        cap,
        mass_goal(),
        options().design,
        constraints,
        &mut split_control,
    )
    .unwrap();
    assert!(restored.work().linear_iterations > before_restore.linear_iterations);
    restored.run(1).unwrap();
    let actual = restored.checkpoint();
    let mut optimizer = expected_checkpoint.optimizer.clone();
    optimizer.work.evaluations += restore_cost;
    assert_eq!(actual.optimizer, optimizer);
    assert_eq!(actual.restoration_evaluations, restore_cost);
    assert_eq!(actual.history, expected_checkpoint.history);
    assert_eq!(
        actual.best_feasible_density,
        expected_checkpoint.best_feasible_density
    );
    assert_eq!(actual.binding, expected_checkpoint.binding);
    assert_eq!(format!("{:?}", actual.gradient_check), initial_audit);
    assert_eq!(
        restored.accepted().evaluation.worst_displacement,
        expected
            .accepted
            .as_ref()
            .unwrap()
            .evaluation
            .worst_displacement
    );
    assert_eq!(
        restored.accepted().evaluation.basis_displacements,
        expected
            .accepted
            .as_ref()
            .unwrap()
            .evaluation
            .basis_displacements
    );
    let checked = restored.recheck_best().unwrap().unwrap();
    let expected_best = expected.recheck.as_ref().unwrap().as_ref().unwrap();
    assert_eq!(checked.rho, expected_best.rho);
    assert_eq!(
        checked.evaluation.worst_displacement,
        expected_best.evaluation.worst_displacement
    );
    assert!(
        checked
            .evaluation
            .responses
            .as_ref()
            .unwrap()
            .physically_feasible()
    );
    assert!(
        checked
            .evaluation
            .stress
            .as_ref()
            .unwrap()
            .physically_feasible()
    );
    assert_eq!(source.operator().elasticity().scales(), source_scales);
    assert_eq!(source.material_body_load_fingerprint(), source_fingerprint);
    assert_eq!(source.physical_regions().unwrap(), source_regions);
}

#[test]
fn failed_refit_preparation_or_recheck_keeps_the_source_and_any_admitted_checkpoint() {
    let tree = Octree3::uniform(1, 3, 512).unwrap();
    let source = adaptive_material_study(&tree);
    let rho = vec![0.7; source.cells()];
    let fine_tree = tree
        .refined(&[source.operator().elasticity().leaves()[2]], poll)
        .unwrap();
    let scales = source.operator().elasticity().scales().to_vec();
    let fingerprint = source.material_body_load_fingerprint();
    let regions = source.physical_regions().unwrap().to_vec();
    let laws = design_laws();
    let limits = limits();
    let requirements = constraints(&limits);
    let cap = AffineLoadSpace3::new(laws.len(), 1.0).unwrap();

    let mut stopped = |progress: SolveProgress| {
        if progress.stage == "sdf3-affine-refit-transfer" {
            ControlFlow::Break(())
        } else {
            poll()
        }
    };
    let mut control = SolveControl::new(SolveBudget::default(), &mut stopped);
    let refused = source.prepare_affine_load_space_refit_with_constraints(
        adaptive_on_tree(&fine_tree),
        &rho,
        laws.len(),
        &cap,
        mass_goal(),
        options(),
        requirements,
        &mut control,
    );
    assert!(matches!(
        refused,
        Err(AffineLoadSpaceRefitError3::Goal(
            GoalRefinementError3::Evaluation(EvaluationStop::Cancelled)
        ))
    ));
    assert_eq!(control.work().linear_solves, 0);

    let mut stopped = |progress: SolveProgress| {
        if progress.stage == "sdf3-affine-design-recheck" {
            ControlFlow::Break(())
        } else {
            poll()
        }
    };
    let mut control = SolveControl::new(SolveBudget::default(), &mut stopped);
    let retained = source
        .refit_reference_affine_load_space_with_constraints(
            adaptive_on_tree(&fine_tree),
            &rho,
            ReferenceLoad3::body(&refinement_dead),
            &laws,
            &cap,
            mass_goal(),
            options(),
            requirements,
            &mut control,
        )
        .unwrap();
    assert!(retained.outcome.is_ok());
    assert!(matches!(
        retained.recheck,
        Some(Err(EvaluationStop::Cancelled))
    ));
    assert!(!retained.completed_and_rechecked());
    assert!(retained.gradient_check.as_ref().unwrap().passed());
    let checkpoint = retained.checkpoint.as_ref().unwrap();
    assert_eq!(checkpoint.optimizer.work.iterations, 2);
    assert_eq!(
        checkpoint.optimizer.point,
        retained.accepted.as_ref().unwrap().rho
    );
    assert_eq!(checkpoint.history, retained.history);
    assert_eq!(
        checkpoint.best_feasible_density,
        retained.best.as_ref().map(|point| point.rho.clone())
    );
    assert!(retained.work.linear_iterations > 0);
    assert_eq!(source.operator().elasticity().scales(), scales);
    assert_eq!(source.material_body_load_fingerprint(), fingerprint);
    assert_eq!(source.physical_regions().unwrap(), regions);
}
