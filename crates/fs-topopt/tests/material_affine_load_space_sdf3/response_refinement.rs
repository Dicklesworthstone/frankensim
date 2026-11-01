//! G1/G3: every adaptive grid binds the same physical response laws anew.
use super::*;
use fs_topopt::sdf3::response::limits::refinement::{
    ReferenceBoundedResponse3, ReferenceResponseLimit3,
};
use fs_topopt::sdf3::robust::envelope::affine_load_space::response::reference::ReferenceAffineResponseLimits3;
use std::cell::Cell;

fn observation(p: [f64; 3]) -> [f64; 3] {
    [0.0, (1.0 + 2.0 * p[0] * p[0]) * (1.0 + p[1]), 0.4 * p[2]]
}

fn calibrated_limits(
    study: &mut CutDensityStudy3<AdaptiveSolveSpace3>,
    rho: &[f64],
    dead: ReferenceLoad3<'_>,
    loads: &[ReferenceLoad3<'_>],
    control: &mut SolveControl<'_>,
) -> (AffineLoadSpace3, [ReferenceResponseLimit3<'static>; 2]) {
    let baseline = reference_baseline(study, rho, dead, loads, control);
    let q = study
        .operator()
        .elasticity()
        .reference_load(ReferenceLoad3::body(&observation), poll)
        .unwrap();
    let columns: Vec<_> = baseline
        .basis_displacements
        .iter()
        .map(|u| dot(&q, u))
        .collect();
    let magnitude = dot(&q, &baseline.offset_displacement).abs() + dot(&columns, &columns).sqrt();
    assert!(magnitude.is_finite() && magnitude > 0.0);
    // Loose physical limits isolate refinement wiring while native scales keep
    // the joint constraint's finite differences resolvable for these small loads.
    let scale = 8.0 * magnitude;
    let cap = AffineLoadSpace3::new(loads.len(), 16.0 * baseline.worst_compliance).unwrap();
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
    (cap, limits)
}

#[test]
fn affine_response_campaign_reintegrates_observations_and_rechecks_the_promoted_kg_design() {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let mut study = adaptive_material_study(&tree);
    let rho = vec![0.7; study.cells()];
    let old_dimension = study.operator().n();
    let old_fingerprint = study.material_body_load_fingerprint();
    let dead = ReferenceLoad3::body(&refinement_dead);
    let loads = design_laws();
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let (cap, limits) = calibrated_limits(&mut study, &rho, dead, &loads, &mut control);
    let family = ReferenceAffineResponseLimits3 {
        limits: &limits,
        options: AffineResponseOptions3::default(),
    };
    let initial = study
        .solve_reference_affine_load_space_with_responses(
            &rho,
            dead,
            &loads,
            &cap,
            mass_goal(),
            AffineLoadSpaceRefitOptions3 {
                accepted_steps: 1,
                ..Default::default()
            },
            family,
            &mut control,
        )
        .unwrap();
    assert!(initial.completed_and_rechecked(), "{:?}", initial.outcome);
    assert_eq!(
        initial
            .gradient_check
            .as_ref()
            .unwrap()
            .constraints
            .responses,
        Some(family.options)
    );
    let builds = RefCell::new(Vec::new());
    let mut campaign = AffineLoadSpaceCampaign3::new_with_response_limits(
        initial,
        tree,
        dead,
        &loads,
        family,
        AffineLoadSpaceCampaignOptions3 {
            max_refits: 1,
            max_candidate_attempts: 1,
            change_tolerance: 0.0,
            response_change_tolerance: 1e6,
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
    assert!(campaign.current_probe().unwrap().responses.is_some());
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
        current
            .gradient_check
            .as_ref()
            .unwrap()
            .constraints
            .responses,
        Some(family.options)
    );
    let checked = current.recheck.as_ref().unwrap().as_ref().unwrap();
    let initial = current.initial.as_ref().unwrap();
    assert_eq!(initial.rho, current.inherited_rho);
    assert!(
        checked.evaluation.physical.dilated_mass_kg.unwrap()
            < initial.evaluation.physical.dilated_mass_kg.unwrap()
    );
    let response = checked.evaluation.responses.as_ref().unwrap();
    assert!(response.physically_feasible());
    assert_eq!(response.options, family.options);
    assert!(checked.evaluation.worst_compliance <= cap.reference());

    // Independently integrate the same observation law in the NEW native DOF
    // ordering, then recover its entire response ball from actual solved fields.
    let native = adaptive_on_tree(campaign.tree());
    let q = native
        .elasticity()
        .reference_load(ReferenceLoad3::body(&observation), poll)
        .unwrap();
    assert_eq!(q.len(), checked.evaluation.offset_displacement.len());
    assert!(q.len() > old_dimension);
    let baseline = dot(&q, &checked.evaluation.offset_displacement);
    let columns: Vec<_> = checked
        .evaluation
        .basis_displacements
        .iter()
        .map(|u| dot(&q, u))
        .collect();
    let radius = dot(&columns, &columns).sqrt();
    for (i, witness) in response.witnesses.iter().enumerate() {
        close(witness.baseline, baseline, 1e-9);
        close_vector(&witness.basis_responses, &columns);
        let sign = if i == 0 { 1.0 } else { -1.0 };
        close(witness.worst_value, baseline + sign * radius, 1e-9);
        assert_eq!(witness.sense, limits[i].sense);
        assert_eq!(witness.bound, limits[i].bound);
        assert_eq!(witness.scale, limits[i].scale);
    }
    let probe = campaign.current_probe().unwrap();
    assert_eq!(probe.rho, checked.rho);
    let responses = probe.responses.as_ref().unwrap();
    assert_eq!(responses.normalized_model_changes.len(), limits.len());
    close(responses.coarse.constraint, response.constraint, 1e-9);
    assert!(responses.coarse.physically_feasible());
    assert!(responses.fine.physically_feasible());
    let probe_tree = campaign
        .tree()
        .refined(current.study.operator().elasticity().leaves(), poll)
        .unwrap();
    let probe_native = adaptive_on_tree(&probe_tree);
    let probe_q = probe_native
        .elasticity()
        .reference_load(ReferenceLoad3::body(&observation), poll)
        .unwrap();
    assert_eq!(probe_q.len(), probe.fine.offset_displacement.len());
    let probe_baseline = dot(&probe_q, &probe.fine.offset_displacement);
    let probe_columns: Vec<_> = probe
        .fine
        .basis_displacements
        .iter()
        .map(|u| dot(&probe_q, u))
        .collect();
    for witness in &responses.fine.witnesses {
        close(witness.baseline, probe_baseline, 1e-9);
        close_vector(&witness.basis_responses, &probe_columns);
    }
    let mut largest = 0.0_f64;
    for ((coarse, fine), &change) in responses
        .coarse
        .witnesses
        .iter()
        .zip(&responses.fine.witnesses)
        .zip(&responses.normalized_model_changes)
    {
        let difference: Vec<_> = fine
            .basis_responses
            .iter()
            .zip(&coarse.basis_responses)
            .map(|(f, c)| f - c)
            .collect();
        let full_change = ((fine.baseline - coarse.baseline).abs()
            + dot(&difference, &difference).sqrt())
            / coarse.scale;
        close(change, full_change, 1e-9);
        largest = largest.max(change);
    }
    close(responses.maximum_model_change, largest, 1e-12);
    assert!(responses.maximum_model_change > 0.0);
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
    assert!(inventory.fine_mass_kg.unwrap() > 0.0);
    assert!(campaign.pending().is_none());
    assert!(campaign.rejected().is_none());
}

#[test]
fn affine_campaign_cannot_drop_response_constraints_when_accepting_a_rechecked_design() {
    let tree = Octree3::uniform(1, 4, 4096).unwrap();
    let mut study = adaptive_material_study(&tree);
    let rho = vec![0.7; study.cells()];
    let dead = ReferenceLoad3::body(&refinement_dead);
    let loads = design_laws();
    let mut callback = |_| poll();
    let mut control = SolveControl::new(SolveBudget::default(), &mut callback);
    let (cap, limits) = calibrated_limits(&mut study, &rho, dead, &loads, &mut control);
    let family = ReferenceAffineResponseLimits3 {
        limits: &limits,
        options: AffineResponseOptions3::default(),
    };
    let initial = study
        .solve_reference_affine_load_space_with_responses(
            &rho,
            dead,
            &loads,
            &cap,
            mass_goal(),
            AffineLoadSpaceRefitOptions3 {
                accepted_steps: 0,
                ..Default::default()
            },
            family,
            &mut control,
        )
        .unwrap();
    assert!(initial.completed_and_rechecked(), "{:?}", initial.outcome);
    assert!(
        initial
            .recheck
            .as_ref()
            .unwrap()
            .as_ref()
            .unwrap()
            .evaluation
            .responses
            .is_some()
    );
    let work = control.work();
    let builds = Cell::new(0);
    let refused = AffineLoadSpaceCampaign3::new(
        initial,
        tree,
        dead,
        &loads,
        AffineLoadSpaceCampaignOptions3::default(),
        &mut control,
        |tree, _, checkpoint| {
            builds.set(builds.get() + 1);
            adaptive_on_tree_controlled(tree, checkpoint)
        },
    );
    assert!(refused.is_err());
    drop(refused);
    assert_eq!(builds.get(), 0);
    assert_eq!(control.work(), work);
}
