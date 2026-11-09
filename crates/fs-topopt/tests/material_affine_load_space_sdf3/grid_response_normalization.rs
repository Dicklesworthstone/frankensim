//! Distinct native-grid observation laws retain one signed physical policy.
use super::*;
use fs_topopt::SolveWork;
use fs_topopt::sdf3::response::limits::refinement::{
    ReferenceBoundedResponse3, ReferenceResponseLimit3,
};
use fs_topopt::sdf3::robust::envelope::affine_load_space::constraints::ReferenceAffineLoadSpaceConstraints3;
use fs_topopt::sdf3::robust::envelope::affine_load_space::response::reference::ReferenceAffineResponseLimits3;
use fs_topopt::sdf3::robust::envelope::affine_load_space::stress::AffineStressLimit3;
use fs_topopt::sdf3::stress::StressOptions3;
use std::cell::Cell;

fn limits(observation: &dyn Fn([f64; 3]) -> [f64; 3]) -> [ReferenceResponseLimit3<'_>; 2] {
    let response = ReferenceBoundedResponse3::Displacement(ReferenceLoad3::body(observation));
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

fn requirements<'a>(
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

#[test]
fn grid_specific_affine_response_normalizers_keep_native_joint_physics_and_shared_primals() {
    let tree = Octree3::uniform(1, 3, 512).unwrap();
    let mut study = adaptive_material_study(&tree);
    let fine = adaptive_operator(2);
    let rho: Vec<_> = (0..study.cells()).map(|i| 0.7 + 0.01 * i as f64).collect();
    let scales = study.operator().elasticity().scales().to_vec();
    let fingerprint = study.material_body_load_fingerprint();
    let coarse_volume: f64 = study.operator().elasticity().volumes().iter().sum();
    let fine_volume: f64 = fine.elasticity().volumes().iter().sum();
    // Distinct declared normalizers deliberately isolate per-grid binding even
    // though this planar fixture happens to integrate the same physical volume.
    // Reusing the coarse callback would double the intended fine observation.
    let coarse_observation = |_: [f64; 3]| [0.0, 0.0, 1.0 / coarse_volume];
    let fine_observation = |_: [f64; 3]| [0.0, 0.0, 1.0 / (2.0 * fine_volume)];
    let coarse_limits = limits(&coarse_observation);
    let fine_limits = limits(&fine_observation);
    let coarse_constraints = requirements(&coarse_limits);
    let fine_constraints = requirements(&fine_limits);
    let loads = design_laws();
    let cap = AffineLoadSpace3::new(loads.len(), 1.0).unwrap();
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let report = study
        .estimate_affine_load_space_enrichment_for_goal_with_grid_constraints(
            fine,
            &rho,
            ReferenceLoad3::body(&refinement_dead),
            &loads,
            RobustOptions3::default(),
            &cap,
            mass_goal(),
            coarse_constraints,
            fine_constraints,
            AffineLoadSpaceRefinementOptions3::default(),
            &mut control,
        )
        .unwrap();
    assert_eq!(report.rho, rho);
    assert_eq!(report.constraints, coarse_constraints.policy());
    assert_eq!(report.work, control.work());
    assert!(report.work.linear_solves <= 1 + 2 * (loads.len() + 1));
    assert_eq!(study.operator().elasticity().scales(), scales);
    assert_eq!(study.material_body_load_fingerprint(), fingerprint);
    let responses = report.responses.as_ref().unwrap();
    let stress = report.stress.as_ref().unwrap();
    assert!(responses.coarse.adjoints.is_empty() && responses.fine.adjoints.is_empty());
    assert!(responses.coarse.gradient.is_none() && responses.fine.gradient.is_none());
    assert!(stress.coarse.adjoints.is_empty() && stress.fine.adjoints.is_empty());
    assert!(stress.coarse.gradient.is_none() && stress.fine.gradient.is_none());
    assert_eq!(stress.coarse.limit, coarse_constraints.stress.unwrap());
    assert_eq!(stress.fine.limit, fine_constraints.stress.unwrap());
    assert_eq!(stress.changes.len(), stress.fine.point_count);
    assert!(report.inventory.is_some());

    let native_fine = adaptive_operator(2);
    for (native, grid, response, observation) in [
        (
            study.operator().elasticity(),
            &report.coarse,
            &responses.coarse,
            &coarse_observation as &dyn Fn([f64; 3]) -> [f64; 3],
        ),
        (
            native_fine.elasticity(),
            &report.fine,
            &responses.fine,
            &fine_observation as &dyn Fn([f64; 3]) -> [f64; 3],
        ),
    ] {
        let q = native
            .reference_load(ReferenceLoad3::body(observation), poll)
            .unwrap();
        let baseline = dot(&q, &grid.offset_displacement);
        let basis: Vec<_> = grid
            .basis_displacements
            .iter()
            .map(|u| dot(&q, u))
            .collect();
        let radius = dot(&basis, &basis).sqrt();
        assert!(radius > 0.0);
        for (i, witness) in response.witnesses.iter().enumerate() {
            close(witness.baseline, baseline, 1e-9);
            close_vector(&witness.basis_responses, &basis);
            close(
                witness.worst_value,
                baseline + if i == 0 { radius } else { -radius },
                1e-9,
            );
            close(witness.witness_value, witness.worst_value, 1e-9);
        }
    }
    let wrong_fine_q = native_fine
        .elasticity()
        .reference_load(ReferenceLoad3::body(&coarse_observation), poll)
        .unwrap();
    let wrong_baseline = dot(&wrong_fine_q, &report.fine.offset_displacement);
    let actual_baseline = responses.fine.witnesses[0].baseline;
    assert!(actual_baseline.abs() > 1e-10);
    assert!((wrong_baseline - actual_baseline).abs() > 0.5 * actual_baseline.abs());
    for ((coarse, fine), change) in responses
        .coarse
        .witnesses
        .iter()
        .zip(&responses.fine.witnesses)
        .zip(&responses.normalized_model_changes)
    {
        let basis_change = coarse
            .basis_responses
            .iter()
            .zip(&fine.basis_responses)
            .map(|(a, b)| (b - a).powi(2))
            .sum::<f64>()
            .sqrt();
        close(
            *change,
            ((fine.baseline - coarse.baseline).abs() + basis_change) / coarse.scale,
            1e-9,
        );
    }
    let joint = report.joint_constraints().unwrap();
    for (i, (grid, response, stress)) in [
        (&report.coarse, &responses.coarse, &stress.coarse),
        (&report.fine, &responses.fine, &stress.fine),
    ]
    .into_iter()
    .enumerate()
    {
        let values = [grid.maximum - 1.0, response.constraint, stress.constraint];
        let maximum = values.into_iter().fold(f64::NEG_INFINITY, f64::max);
        let expected = maximum
            + values
                .iter()
                .map(|v| (report.constraints.sharpness * (v - maximum)).exp())
                .sum::<f64>()
                .ln()
                / report.constraints.sharpness;
        close(joint[i], expected, 1e-9);
    }
}

#[test]
fn incompatible_affine_grid_response_families_refuse_before_callbacks_or_solver_work() {
    let tree = Octree3::uniform(1, 3, 512).unwrap();
    let mut study = adaptive_material_study(&tree);
    let rho = vec![0.7; study.cells()];
    let scales = study.operator().elasticity().scales().to_vec();
    let fingerprint = study.material_body_load_fingerprint();
    for mismatch in [
        "bound", "scale", "sense", "kind", "options", "stress", "joint",
    ] {
        let invoked = Cell::new(0);
        let observation = |_: [f64; 3]| {
            invoked.set(invoked.get() + 1);
            [0.0, 0.0, 1.0]
        };
        let reaction = |_: [f64; 3], _: [f64; 3]| [0.0, 0.0, 1.0];
        let coarse_limits = limits(&observation);
        let mut fine_limits = coarse_limits;
        match mismatch {
            "bound" => fine_limits[0].bound = 2.0,
            "scale" => fine_limits[0].scale = 2.0,
            "sense" => fine_limits[0].sense = ResponseLimitSense3::AtLeast,
            "kind" => fine_limits[0].response = ReferenceBoundedResponse3::Reaction(&reaction),
            _ => {}
        }
        let coarse = requirements(&coarse_limits);
        let mut fine = requirements(&fine_limits);
        match mismatch {
            "options" => fine.responses.as_mut().unwrap().options.sharpness = 50.0,
            "stress" => fine.stress.as_mut().unwrap().limit = 50.0,
            "joint" => fine.sharpness = 50.0,
            _ => {}
        }
        let checkpoints = Cell::new(0);
        let mut callback = |_| {
            checkpoints.set(checkpoints.get() + 1);
            poll()
        };
        let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
        let result = study.estimate_affine_load_space_enrichment_for_goal_with_grid_constraints(
            adaptive_operator(2),
            &rho,
            ReferenceLoad3::default(),
            &design_laws(),
            RobustOptions3::default(),
            &policy(),
            mass_goal(),
            coarse,
            fine,
            AffineLoadSpaceRefinementOptions3::default(),
            &mut control,
        );
        assert!(
            matches!(
                result,
                Err(GoalRefinementError3::Invalid(
                    "coarse/fine affine constraint declarations differ"
                ))
            ),
            "mismatch: {mismatch}"
        );
        assert_eq!(invoked.get(), 0, "mismatch: {mismatch}");
        assert_eq!(checkpoints.get(), 0, "mismatch: {mismatch}");
        assert_eq!(control.work(), SolveWork::default(), "mismatch: {mismatch}");
        assert_eq!(study.operator().elasticity().scales(), scales);
        assert_eq!(study.material_body_load_fingerprint(), fingerprint);
    }
}
