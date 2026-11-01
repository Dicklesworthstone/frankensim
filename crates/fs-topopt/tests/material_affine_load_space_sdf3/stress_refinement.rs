//! G1/G3: adaptive strength, response and mass remain one physical problem.
use super::*;
use fs_topopt::sdf3::response::limits::refinement::{
    ReferenceBoundedResponse3, ReferenceResponseLimit3,
};
use fs_topopt::sdf3::robust::envelope::affine_load_space::constraints::{
    AffineLoadSpaceConstraints3, ReferenceAffineLoadSpaceConstraints3,
};
use fs_topopt::sdf3::robust::envelope::affine_load_space::response::reference::ReferenceAffineResponseLimits3;
use fs_topopt::sdf3::robust::envelope::affine_load_space::stress::AffineStressLimit3;
use fs_topopt::sdf3::stress::StressOptions3;
use std::cell::Cell;

fn observation(p: [f64; 3]) -> [f64; 3] {
    [0.0, (1.0 + 2.0 * p[0] * p[0]) * (1.0 + p[1]), 0.4 * p[2]]
}

fn refit_options(accepted_steps: usize) -> AffineLoadSpaceRefitOptions3 {
    let mut design = AffineLoadSpaceDesignOptions3 {
        density_floor: 0.6,
        ..Default::default()
    };
    // Leave room for a measurable accepted update on both actual grids.
    design.optimizer.max_spectral_step = 0.1;
    AffineLoadSpaceRefitOptions3 {
        design,
        accepted_steps,
        ..Default::default()
    }
}

fn calibrated_limits(
    study: &mut CutDensityStudy3<AdaptiveSolveSpace3>,
    rho: &[f64],
    dead: ReferenceLoad3<'_>,
    loads: &[ReferenceLoad3<'_>],
    control: &mut SolveControl<'_>,
) -> (
    AffineLoadSpace3,
    AffineStressLimit3,
    [ReferenceResponseLimit3<'static>; 2],
) {
    let forces: Vec<_> = std::iter::once(dead)
        .chain(loads.iter().copied())
        .map(|law| {
            study
                .operator()
                .elasticity()
                .reference_load(law, poll)
                .unwrap()
        })
        .collect();
    let basis: Vec<_> = forces[1..].iter().map(Vec::as_slice).collect();
    let mut strength = AffineStressLimit3 {
        limit: 1.0,
        options: StressOptions3 {
            relaxation_power: 2.0,
            aggregation_power: 2.0,
            ..Default::default()
        },
    };
    let baseline = study
        .evaluate_robust_affine_constraints(
            rho,
            &forces[0],
            &basis,
            RobustOptions3::default(),
            &policy(),
            AffineLoadSpaceConstraints3 {
                stress: Some(strength),
                ..Default::default()
            },
            control,
        )
        .unwrap();
    let q = study
        .operator()
        .elasticity()
        .reference_load(ReferenceLoad3::body(&observation), poll)
        .unwrap();
    let responses: Vec<_> = baseline
        .basis_displacements
        .iter()
        .map(|u| dot(&q, u))
        .collect();
    let magnitude =
        dot(&q, &baseline.offset_displacement).abs() + dot(&responses, &responses).sqrt();
    assert!(magnitude > 0.0);
    let scale = 8.0 * magnitude;
    let response = ReferenceBoundedResponse3::Displacement(ReferenceLoad3::body(&observation));
    let limits = [
        ReferenceResponseLimit3 {
            response,
            sense: ResponseLimitSense3::AtMost,
            bound: scale,
            scale,
        },
        ReferenceResponseLimit3 {
            response,
            sense: ResponseLimitSense3::AtLeast,
            bound: -scale,
            scale,
        },
    ];
    // The p-norm grows with sample count. Keep all physical caps loose for this
    // bounded campaign, while native normalization resolves the joint FD audit.
    strength.limit = 16.0 * baseline.stress.as_ref().unwrap().aggregate;
    let cap = AffineLoadSpace3::new(loads.len(), 16.0 * baseline.worst_compliance).unwrap();
    (cap, strength, limits)
}

fn von_mises(s: &[f64; 6]) -> f64 {
    (0.5 * ((s[0] - s[1]).powi(2) + (s[1] - s[2]).powi(2) + (s[2] - s[0]).powi(2))
        + 3.0 * (s[3].powi(2) + s[4].powi(2) + s[5].powi(2)))
    .sqrt()
}

#[test]
fn affine_strength_campaign_rechecks_joint_physics_and_all_matched_stress_directions() {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let mut study = adaptive_material_study(&tree);
    let rho = vec![0.7; study.cells()];
    let old_dimension = study.operator().n();
    let old_fingerprint = study.material_body_load_fingerprint();
    let dead = ReferenceLoad3::body(&refinement_dead);
    let loads = design_laws();
    let stress_work = RefCell::new(Vec::new());
    let mut callback = |progress: SolveProgress| {
        if matches!(
            progress.stage,
            "sdf3-affine-stress-refinement-start" | "sdf3-affine-stress-refinement-publish"
        ) {
            stress_work
                .borrow_mut()
                .push((progress.stage, progress.work));
        }
        poll()
    };
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let (cap, strength, limits) = calibrated_limits(&mut study, &rho, dead, &loads, &mut control);
    let requirements = ReferenceAffineLoadSpaceConstraints3 {
        responses: Some(ReferenceAffineResponseLimits3 {
            limits: &limits,
            options: AffineResponseOptions3::default(),
        }),
        stress: Some(strength),
        ..Default::default()
    };
    let initial = study
        .solve_reference_affine_load_space_with_constraints(
            &rho,
            dead,
            &loads,
            &cap,
            mass_goal(),
            refit_options(1),
            requirements,
            &mut control,
        )
        .unwrap();
    assert!(initial.completed_and_rechecked(), "{:?}", initial.outcome);
    assert_eq!(
        initial.gradient_check.as_ref().unwrap().constraints,
        requirements.policy()
    );
    let builds = RefCell::new(Vec::new());
    let mut campaign = AffineLoadSpaceCampaign3::new_with_constraints(
        initial,
        tree,
        dead,
        &loads,
        requirements,
        AffineLoadSpaceCampaignOptions3 {
            max_refits: 1,
            max_candidate_attempts: 1,
            change_tolerance: 0.0,
            response_change_tolerance: 1e6,
            stress_change_tolerance: 1e6,
            volume_tolerance: 1e-6,
            objective_tolerance: 1e-6,
            marking: AffineLoadSpaceMarking3::FrozenWitness {
                theta: 0.5,
                max_marks: 1,
            },
            ..Default::default()
        },
        &mut control,
        |tree, kind, checkpoint| {
            builds.borrow_mut().push(kind);
            adaptive_on_tree_controlled(tree, checkpoint)
        },
    )
    .unwrap();
    assert_eq!(
        campaign.run(0).unwrap(),
        AffineLoadSpaceCampaignStop3::Paused
    );
    assert!(campaign.current_probe().unwrap().stress.is_some());
    assert_eq!(
        campaign.run(1).unwrap(),
        AffineLoadSpaceCampaignStop3::RefinementLimit
    );
    assert_eq!(
        builds.borrow().as_slice(),
        &[
            AffineLoadSpaceBuild3::Probe,
            AffineLoadSpaceBuild3::Refit,
            AffineLoadSpaceBuild3::Probe,
        ]
    );
    assert_eq!(campaign.grids().len(), 2);
    assert_eq!(campaign.probes().len(), 2);
    let current = campaign.current();
    assert!(current.completed_and_rechecked());
    assert!(current.study.operator().n() > old_dimension);
    assert_ne!(
        current.study.material_body_load_fingerprint(),
        old_fingerprint
    );
    assert_eq!(current.goal, mass_goal());
    assert_eq!(current.policy, cap);
    assert_eq!(
        current.gradient_check.as_ref().unwrap().constraints,
        requirements.policy()
    );
    let checked = current.recheck.as_ref().unwrap().as_ref().unwrap();
    let initial = current.initial.as_ref().unwrap();
    assert_eq!(initial.rho, current.inherited_rho);
    assert!(
        checked.evaluation.physical.dilated_mass_kg.unwrap()
            < initial.evaluation.physical.dilated_mass_kg.unwrap()
    );
    for point in [
        initial,
        current.accepted.as_ref().unwrap(),
        current.best.as_ref().unwrap(),
        checked,
    ] {
        assert_eq!(point.evaluation.constraints, requirements.policy());
        assert_eq!(point.evaluation.stress.as_ref().unwrap().limit, strength);
        assert!(point.evaluation.responses.is_some());
    }
    assert!(
        checked
            .evaluation
            .stress
            .as_ref()
            .unwrap()
            .physically_feasible()
    );
    assert!(
        checked
            .evaluation
            .responses
            .as_ref()
            .unwrap()
            .physically_feasible()
    );
    assert!(checked.evaluation.worst_compliance <= cap.reference());
    let probe = campaign.current_probe().unwrap();
    assert_eq!(probe.rho, checked.rho);
    assert_eq!(probe.constraints, requirements.policy());
    let stress = probe.stress.as_ref().unwrap();
    assert_eq!(stress.coarse.limit, strength);
    assert_eq!(stress.fine.limit, strength);
    assert!(stress.stable_and_feasible(1e6));
    assert!(!stress.stable_and_feasible(0.0));
    assert!(stress.coarse.adjoints.is_empty() && stress.fine.adjoints.is_empty());
    assert!(stress.coarse.gradient.is_none() && stress.fine.gradient.is_none());
    close(
        stress.coarse.aggregate,
        checked.evaluation.stress.as_ref().unwrap().aggregate,
        1e-9,
    );
    close(
        stress.normalized_aggregate_change,
        (stress.fine.aggregate - stress.coarse.aggregate).abs() / strength.limit,
        1e-9,
    );
    let work = stress_work.borrow();
    assert_eq!(work.len(), 4);
    for pair in work.chunks_exact(2) {
        assert_eq!(pair[0].0, "sdf3-affine-stress-refinement-start");
        assert_eq!(pair[1].0, "sdf3-affine-stress-refinement-publish");
        assert_eq!(pair[0].1, pair[1].1);
    }

    let probe_tree = campaign
        .tree()
        .refined(current.study.operator().elasticity().leaves(), poll)
        .unwrap();
    let mut coarse = adaptive_on_tree(campaign.tree());
    let mut fine = adaptive_on_tree(&probe_tree);
    coarse.set_scales(&probe.coarse_eroded_scales).unwrap();
    fine.set_scales(&probe.fine_eroded_scales).unwrap();
    let transfer =
        AdaptiveTransfer3::new(coarse.elasticity(), fine.elasticity(), 2_000_000, poll).unwrap();
    for (child, &parent) in transfer.parents().iter().enumerate() {
        assert_eq!(
            probe.fine_eroded_rho[child],
            probe.coarse_eroded_rho[parent]
        );
        assert_eq!(
            probe.fine_eroded_scales[child],
            probe.coarse_eroded_scales[parent]
        );
    }
    // An independent native check: interpolate and subtract displacement fields,
    // then sample their small, admissible strain differences with the fine tensor.
    // Production instead admits both complete physical families before differencing.
    let differences: Vec<_> = std::iter::once(&probe.coarse.offset_displacement)
        .chain(&probe.coarse.basis_displacements)
        .zip(
            std::iter::once(&probe.fine.offset_displacement).chain(&probe.fine.basis_displacements),
        )
        .map(|(a, b)| {
            let prolonged = transfer.prolongate(a, poll).unwrap();
            let delta: Vec<_> = b.iter().zip(prolonged).map(|(f, c)| f - c).collect();
            fine.elasticity()
                .bulk_stress(&delta, strength.options.max_points, poll)
                .unwrap()
        })
        .collect();
    assert!(
        differences
            .iter()
            .all(|points| points.len() == stress.changes.len())
    );
    assert_eq!(stress.changes.len(), stress.fine.point_count);
    let mut maximum = 0.0_f64;
    let mut reference_maximum = 0.0_f64;
    for (index, change) in stress.changes.iter().enumerate() {
        let baseline = &differences[0][index];
        assert_eq!(change.cell, baseline.cell);
        assert_eq!(change.position, baseline.position);
        assert_eq!(change.weight, baseline.weight);
        let reference = (von_mises(&baseline.reference_stress)
            + differences[1..]
                .iter()
                .map(|d| von_mises(&d[index].reference_stress).powi(2))
                .sum::<f64>()
                .sqrt())
            / strength.limit;
        let physical = (von_mises(&baseline.stress)
            + differences[1..]
                .iter()
                .map(|d| von_mises(&d[index].stress).powi(2))
                .sum::<f64>()
                .sqrt())
            / strength.limit;
        let relaxed =
            reference * probe.fine_eroded_rho[change.cell].powf(strength.options.relaxation_power);
        close(change.normalized_reference_change, reference, 1e-7);
        close(change.normalized_relaxed_change, relaxed, 1e-7);
        close(change.normalized_physical_change, physical, 1e-7);
        reference_maximum = reference_maximum.max(reference);
        maximum = maximum.max(relaxed).max(physical);
        for j in 0..64 {
            let angle = std::f64::consts::TAU * j as f64 / 64.0;
            let tensor = std::array::from_fn(|k| {
                baseline.reference_stress[k]
                    + angle.cos() * differences[1][index].reference_stress[k]
                    + angle.sin() * differences[2][index].reference_stress[k]
            });
            assert!(
                von_mises(&tensor) / strength.limit
                    <= change.normalized_reference_change * (1.0 + 1e-8) + 1e-14
            );
        }
    }
    close(stress.maximum_reference_change, reference_maximum, 1e-7);
    close(stress.maximum_model_change, maximum, 1e-7);
    assert!(maximum > 0.0);
    let responses = probe.responses.as_ref().unwrap();
    assert!(responses.coarse.physically_feasible() && responses.fine.physically_feasible());
    for (native, fields, response) in [
        (&coarse, &probe.coarse, &responses.coarse),
        (&fine, &probe.fine, &responses.fine),
    ] {
        let q = native
            .elasticity()
            .reference_load(ReferenceLoad3::body(&observation), poll)
            .unwrap();
        let columns: Vec<_> = fields
            .basis_displacements
            .iter()
            .map(|u| dot(&q, u))
            .collect();
        for witness in &response.witnesses {
            close(witness.baseline, dot(&q, &fields.offset_displacement), 1e-9);
            close_vector(&witness.basis_responses, &columns);
        }
    }
    let inventory = probe.inventory.as_ref().unwrap();
    assert_eq!(
        inventory.objective,
        MinimumMaterialObjective3::MassKg { scale_kg: 2.0 }
    );
    close(
        inventory.coarse_mass_kg.unwrap(),
        checked.evaluation.physical.dilated_mass_kg.unwrap(),
        1e-12,
    );
    let kg: f64 = fine
        .elasticity()
        .reference_cell_masses_kg()
        .unwrap()
        .iter()
        .zip(transfer.parents())
        .map(|(m, &p)| m * checked.evaluation.physical.dilated_rho[p])
        .sum();
    close(inventory.fine_mass_kg.unwrap(), kg, 1e-12);
    assert!(campaign.pending().is_none() && campaign.rejected().is_none());
}

#[test]
fn legacy_affine_campaign_constructors_cannot_drop_an_optimized_strength_cap() {
    for preserve_responses in [false, true] {
        let tree = Octree3::uniform(1, 4, 4096).unwrap();
        let mut study = adaptive_material_study(&tree);
        let rho = vec![0.7; study.cells()];
        let dead = ReferenceLoad3::body(&refinement_dead);
        let loads = design_laws();
        let mut callback = |_| poll();
        let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
        let (cap, strength, limits) =
            calibrated_limits(&mut study, &rho, dead, &loads, &mut control);
        let responses = ReferenceAffineResponseLimits3 {
            limits: &limits,
            options: AffineResponseOptions3::default(),
        };
        let requirements = ReferenceAffineLoadSpaceConstraints3 {
            responses: Some(responses),
            stress: Some(strength),
            ..Default::default()
        };
        let initial = study
            .solve_reference_affine_load_space_with_constraints(
                &rho,
                dead,
                &loads,
                &cap,
                mass_goal(),
                refit_options(0),
                requirements,
                &mut control,
            )
            .unwrap();
        assert!(initial.completed_and_rechecked(), "{:?}", initial.outcome);
        let checked = initial.recheck.as_ref().unwrap().as_ref().unwrap();
        assert!(
            checked
                .evaluation
                .stress
                .as_ref()
                .unwrap()
                .physically_feasible()
        );
        let work = control.work();
        let builds = Cell::new(0);
        let build = |tree: &Octree3,
                     _: AffineLoadSpaceBuild3,
                     checkpoint: &mut dyn FnMut() -> ControlFlow<()>| {
            builds.set(builds.get() + 1);
            adaptive_on_tree_controlled(tree, checkpoint)
        };
        let refused = if preserve_responses {
            AffineLoadSpaceCampaign3::new_with_response_limits(
                initial,
                tree,
                dead,
                &loads,
                responses,
                Default::default(),
                &mut control,
                build,
            )
        } else {
            AffineLoadSpaceCampaign3::new(
                initial,
                tree,
                dead,
                &loads,
                Default::default(),
                &mut control,
                build,
            )
        };
        assert!(refused.is_err());
        drop(refused);
        assert_eq!(builds.get(), 0);
        assert_eq!(control.work(), work);
    }
}

#[test]
fn affine_compliance_campaign_keeps_volume_response_and_strength_limits_on_every_grid() {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let mut study = adaptive_material_study(&tree);
    let rho = vec![0.7; study.cells()];
    let original_dimension = study.operator().n();
    let dead = ReferenceLoad3::body(&refinement_dead);
    let loads = design_laws();
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let (loose_cap, strength, limits) =
        calibrated_limits(&mut study, &rho, dead, &loads, &mut control);
    // This is an objective normalization, deliberately BELOW actual compliance.
    // Reusing the minimum-material cap predicate would falsely reject the run.
    let normalization = AffineLoadSpace3::new(loads.len(), loose_cap.reference() / 64.0).unwrap();
    let volume_limit = 0.95;
    let goal = AffineLoadSpaceDesignGoal3::WorstCompliance { volume_limit };
    let requirements = ReferenceAffineLoadSpaceConstraints3 {
        responses: Some(ReferenceAffineResponseLimits3 {
            limits: &limits,
            options: AffineResponseOptions3::default(),
        }),
        stress: Some(strength),
        ..Default::default()
    };
    let mut options = refit_options(1);
    // Keep a measurable descent step inside the original volume constraint on
    // both grids. Accepted-step completion is not a stationarity assertion.
    options.design.optimizer.max_spectral_step = 0.001;
    let initial = study
        .solve_reference_affine_load_space_with_constraints(
            &rho,
            dead,
            &loads,
            &normalization,
            goal,
            options,
            requirements,
            &mut control,
        )
        .unwrap();
    assert!(initial.completed_and_rechecked(), "{:?}", initial.outcome);
    let checked = initial.recheck.as_ref().unwrap().as_ref().unwrap();
    assert!(checked.evaluation.maximum > 1.0);
    assert!(
        checked.evaluation.worst_compliance
            < initial
                .initial
                .as_ref()
                .unwrap()
                .evaluation
                .worst_compliance
    );
    let first_objective = checked.evaluation.maximum;
    let builds = RefCell::new(Vec::new());
    let mut campaign = AffineLoadSpaceCampaign3::new_with_constraints(
        initial,
        tree,
        dead,
        &loads,
        requirements,
        AffineLoadSpaceCampaignOptions3 {
            max_refits: 1,
            max_candidate_attempts: 1,
            // A genuine nonzero model change forces one native refit.
            change_tolerance: 0.0,
            response_change_tolerance: 1e6,
            stress_change_tolerance: 1e6,
            volume_tolerance: 1e-6,
            objective_tolerance: 1e6,
            ..Default::default()
        },
        &mut control,
        |tree, kind, checkpoint| {
            builds.borrow_mut().push(kind);
            adaptive_on_tree_controlled(tree, checkpoint)
        },
    )
    .unwrap();
    assert_eq!(
        campaign.run(0).unwrap(),
        AffineLoadSpaceCampaignStop3::Paused
    );
    let first_probe = campaign.current_probe().unwrap();
    assert_eq!(first_probe.goal, Some(goal));
    assert!(first_probe.coarse.maximum > 1.0 && first_probe.fine.maximum > 1.0);
    assert!(first_probe.stable_and_feasible_with_constraint_tolerances(1e6, 1e6, 1e6));
    assert!(first_probe.normalized_model_change > 0.0);
    assert_eq!(
        campaign.run(1).unwrap(),
        AffineLoadSpaceCampaignStop3::RefinementLimit
    );
    assert_eq!(
        builds.borrow().as_slice(),
        &[
            AffineLoadSpaceBuild3::Probe,
            AffineLoadSpaceBuild3::Refit,
            AffineLoadSpaceBuild3::Probe,
        ]
    );
    assert_eq!(campaign.grids().len(), 2);
    assert_eq!(campaign.probes().len(), 2);
    close(campaign.grids()[0].objective, first_objective, 1e-10);
    for grid in campaign.grids() {
        close(
            grid.objective,
            grid.worst_compliance / normalization.reference(),
            1e-9,
        );
        assert!(grid.objective > 1.0);
        assert!(grid.dilated_volume_fraction <= volume_limit);
    }

    let current = campaign.current();
    assert!(current.completed_and_rechecked(), "{:?}", current.outcome);
    assert!(current.study.operator().n() > original_dimension);
    assert_eq!(current.goal, goal);
    assert_eq!(current.policy, normalization);
    let audit = current.gradient_check.as_ref().unwrap();
    assert!(audit.passed());
    assert_eq!(audit.goal, goal);
    assert_eq!(audit.constraints, requirements.policy());
    let inherited = current.initial.as_ref().unwrap();
    let checked = current.recheck.as_ref().unwrap().as_ref().unwrap();
    // Compare the fine endpoint with its own freshly solved starting design.
    // Different grids need not have ordered objective values.
    assert!(checked.evaluation.worst_compliance < inherited.evaluation.worst_compliance);
    assert!(checked.evaluation.maximum > 1.0);
    assert_eq!(inherited.rho, current.inherited_rho);
    for point in [
        inherited,
        current.accepted.as_ref().unwrap(),
        current.best.as_ref().unwrap(),
        checked,
    ] {
        let evaluation = &point.evaluation;
        assert_eq!(evaluation.constraints, requirements.policy());
        assert!(evaluation.physical.dilated_volume_fraction <= volume_limit);
        let stress = evaluation.stress.as_ref().unwrap();
        assert_eq!(stress.limit, strength);
        assert!(stress.physically_feasible());
        let responses = evaluation.responses.as_ref().unwrap();
        assert!(responses.physically_feasible());
        for (witness, limit) in responses.witnesses.iter().zip(&limits) {
            assert_eq!(witness.bound, limit.bound);
            assert_eq!(witness.sense, limit.sense);
            assert_eq!(witness.scale, limit.scale);
        }
    }
    let probe = campaign.current_probe().unwrap();
    assert_eq!(probe.goal, Some(goal));
    assert_eq!(probe.rho, checked.rho);
    assert_eq!(probe.constraints, requirements.policy());
    assert!(probe.coarse.maximum > 1.0 && probe.fine.maximum > 1.0);
    assert!(probe.stable_and_feasible_with_constraint_tolerances(1e6, 1e6, 1e6));
    let inventory = probe.inventory.as_ref().unwrap();
    assert_eq!(
        inventory.objective,
        MinimumMaterialObjective3::VolumeFraction
    );
    assert_eq!(inventory.volume_scale, volume_limit);
    close(
        inventory.coarse_volume_fraction,
        checked.evaluation.physical.dilated_volume_fraction,
        1e-12,
    );
    assert!(inventory.coarse_volume_fraction <= volume_limit);
    assert!(inventory.fine_volume_fraction <= volume_limit);
    close(
        inventory.coarse_objective,
        inventory.coarse_volume_fraction / volume_limit,
        1e-12,
    );
    close(
        probe.objective_change().unwrap(),
        (probe.fine.maximum - probe.coarse.maximum).abs(),
        1e-12,
    );
    close(
        campaign.probes().last().unwrap().objective_change,
        probe.objective_change().unwrap(),
        1e-12,
    );
    let stress = probe.stress.as_ref().unwrap();
    assert_eq!(stress.coarse.limit, strength);
    assert_eq!(stress.fine.limit, strength);
    assert!(stress.coarse.physically_feasible() && stress.fine.physically_feasible());
    let responses = probe.responses.as_ref().unwrap();
    assert!(responses.coarse.physically_feasible() && responses.fine.physically_feasible());

    // Predicate regression on an otherwise admitted native report: no change
    // tolerance can excuse even a small excess of the original fine volume cap.
    let mut excessive = probe.clone();
    let inventory = excessive.inventory.as_mut().unwrap();
    inventory.fine_volume_fraction = volume_limit + 1e-8;
    inventory.fine_objective = inventory.fine_volume_fraction / volume_limit;
    inventory.volume_change =
        (inventory.fine_volume_fraction - inventory.coarse_volume_fraction).abs();
    inventory.objective_change = (inventory.fine_objective - inventory.coarse_objective).abs();
    assert!(!excessive.stable_and_feasible_with_constraint_tolerances(1e6, 1e6, 1e6));
    assert!(campaign.pending().is_none() && campaign.rejected().is_none());
}
