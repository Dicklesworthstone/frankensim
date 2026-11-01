//! Actual new-grid optimization must preserve sources and improve its own baseline.
use super::*;
use fs_topopt::sdf3::modal::design::{
    Sdf3EigenfrequencyDesignPoint3, Sdf3EigenfrequencyGoal3,
};
use fs_topopt::sdf3::modal::refit::{
    Sdf3EigenfrequencyRefit3, Sdf3EigenfrequencyRefitOptions3,
};
use fs_topopt::sdf3::response::limits::design::MinimumMaterialObjective3;

fn settings(reference: f64) -> Sdf3EigenfrequencyRefitOptions3 {
    Sdf3EigenfrequencyRefitOptions3 {
        design: design_settings(reference),
        accepted_steps: 2,
        ..Default::default()
    }
}

pub(super) fn point(result: &Sdf3EigenfrequencyRefit3<AdaptiveSolveSpace3>)
    -> &Sdf3EigenfrequencyDesignPoint3
{
    result.recheck.as_ref().unwrap().as_ref().unwrap()
}

pub(super) fn assert_physical(result: &Sdf3EigenfrequencyRefit3<AdaptiveSolveSpace3>) {
    assert!(result.completed_and_rechecked());
    assert!(result.gradient_check.as_ref().unwrap().passed());
    let evaluation = &point(result).evaluation;
    match result.goal {
        Sdf3EigenfrequencyGoal3::MaximumFrequency { volume_limit } => {
            assert!(evaluation.volume_fraction <= volume_limit);
        }
        Sdf3EigenfrequencyGoal3::MinimumMaterial { .. } => {
            assert!(evaluation.objective.aggregate >= result.options.design.eigenvalue_reference);
            assert!(evaluation.objective.modal.values[0] >= result.options.design.eigenvalue_reference);
        }
    }
    result.study.operator().validate_material_mass(&result.mass).unwrap();
}

pub(super) fn seed(
    minimum_mass: bool,
    control: &mut SolveControl<'_>,
) -> (Octree3, Sdf3EigenfrequencyRefit3<AdaptiveSolveSpace3>) {
    let tree = Octree3::uniform(0, 3, 512).unwrap();
    let operator = build(&tree, 3.0, &mut poll).unwrap();
    let mass = operator.prepare_material_mass(100_000, poll).unwrap();
    let mut study = CutDensityStudy3::new(operator, 0.2, params())
        .with_material_acceleration(|_| [0.0, 0.0, -9.81], 100_000, poll)
        .unwrap();
    let baseline = study.evaluate_eigenfrequency(&mass, &[0.72], frequency(), control).unwrap();
    assert!(baseline.objective.aggregate > 0.0);
    let goal = if minimum_mass {
        Sdf3EigenfrequencyGoal3::MinimumMaterial {
            objective: MinimumMaterialObjective3::MassKg {
                scale_kg: baseline.physical_mass_kg,
            },
            volume_scale: 1.0,
        }
    } else {
        Sdf3EigenfrequencyGoal3::MaximumFrequency { volume_limit: 0.98 }
    };
    let reference = baseline.objective.aggregate * if minimum_mass { 0.05 } else { 1.0 };
    let result = study.solve_eigenfrequency_design(&[0.72], goal, settings(reference), control).unwrap();
    assert_physical(&result);
    (tree, result)
}

#[test]
fn refits_use_new_native_mass_and_improve_their_own_grid_for_both_goals() {
    for minimum_mass in [false, true] {
        let mut callback = |_| poll();
        let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
        let (tree, source) = seed(minimum_mass, &mut control);
        let marks: Vec<_> = tree.leaves().iter().copied().collect();
        let fine_tree = tree.refined(&marks, poll).unwrap();
        let fine = build(&fine_tree, 3.0, &mut poll).unwrap();
        let source_n = source.study.operator().n();
        let source_scales = source.study.operator().scales().to_vec();
        let source_body = source.study.material_body_load_fingerprint();
        let rho = point(&source).rho.clone();
        let result = source.study.refit_eigenfrequency_design(
            fine, &rho, source.goal, source.options, &mut control,
        ).unwrap();
        assert_physical(&result);
        assert!(result.study.cells() > source.study.cells());
        assert!(result.study.operator().n() > source_n);
        assert_eq!(result.inherited_rho, vec![rho[0]; result.study.cells()]);
        assert!(result.study.operator().validate_material_mass(&source.mass).is_err());
        assert_eq!(source.study.operator().scales(), source_scales);
        assert_eq!(source.study.operator().n(), source_n);
        assert_eq!(source.study.material_body_load_fingerprint(), source_body);
        assert!(source_body.is_some());
        assert!(result.study.material_body_load_fingerprint().is_some());
        assert_ne!(result.study.material_body_load_fingerprint(), source_body);
        assert!(result.study.can_refine_material_body_load());
        let before = &result.initial.as_ref().unwrap().evaluation;
        let after = &point(&result).evaluation;
        if minimum_mass {
            assert!(after.physical_mass_kg < before.physical_mass_kg * (1.0 - 1e-7));
        } else {
            assert!(after.objective.aggregate > before.objective.aggregate * (1.0 + 1e-7));
        }
    }
}
