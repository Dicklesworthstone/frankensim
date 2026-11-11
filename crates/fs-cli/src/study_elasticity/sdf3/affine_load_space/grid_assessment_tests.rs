use super::*;
use fs_cutfem::elastic3::adaptive::enrichment::AdaptiveTransfer3;
use fs_topopt::SolveProgress;
use fs_topopt::sdf3::response::limits::{BoundedResponse3, ResponseLimitSense3};
use fs_topopt::sdf3::robust::envelope::affine_load_space::refinement::AffineLoadSpaceRefinement3;
use spec::responses::Observation;

const FIXTURE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-3d-affine-grid-assessment.fsim"
));

fn fixture() -> Spec {
    spec::parse(FIXTURE).unwrap()
}

fn compute(spec: &Spec) -> Computation {
    compute_observed(spec, &CancelGate::new_clock_free(), 0, None, 0.0, |_, _| {
        Ok(())
    })
    .unwrap()
}

fn coarse(spec: &Spec) -> Computation {
    let mut source = spec.clone();
    source.affine_grid_assessment = None;
    source.max_level = source.level;
    compute(&source)
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 2e-7 * a.abs().max(b.abs()).max(1e-12),
        "{a:e} != {b:e}"
    );
}

// A basis response of a limit is compared on that limit's own response scale. `radius` is the norm of all of the limit's basis responses, so it bounds every one of them: a direction
// whose response is at rounding level beside its neighbours' is read against the size of its neighbours (here 0.36) and not against the unitless 1e-12 of `near`, which no response of any size
// beats. In the rows that failed, the stored and the re-measured response differ by 2.5 and 4.2 eps of the radius (a reaction, so no sum of absolute terms exists to scale by).
fn near_to_radius(a: f64, b: f64, radius: f64) {
    assert!(
        (a - b).abs() <= 2e-7 * a.abs().max(b.abs()).max(radius),
        "{a:e} != {b:e} on the response scale {radius:e}"
    );
}

fn difference(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt()
}

fn sign(sense: ResponseLimitSense3) -> f64 {
    if sense == ResponseLimitSense3::AtMost {
        1.0
    } else {
        -1.0
    }
}

fn checked(run: &Computation) -> &AffineLoadSpaceRefinement3 {
    let outcome = run.state.assessment.as_ref().unwrap();
    assert_eq!(outcome.status, "checked", "{:?}", outcome.reason);
    outcome.report.as_ref().unwrap()
}

fn loose(mut spec: Spec, report: &AffineLoadSpaceRefinement3) -> Spec {
    let config = spec.affine_load_space.as_mut().unwrap();
    config.cap_j = 2.0
        * report
            .coarse
            .worst_compliance
            .max(report.fine.worst_compliance);
    let stress = report.stress.as_ref().unwrap();
    config.stress.as_mut().unwrap().limit =
        2.0 * stress.coarse.aggregate.max(stress.fine.aggregate);
    let responses = report.responses.as_ref().unwrap();
    for ((limit, a), b) in spec
        .affine_responses
        .as_mut()
        .unwrap()
        .limits
        .iter_mut()
        .zip(&responses.coarse.witnesses)
        .zip(&responses.fine.witnesses)
    {
        let scale = a
            .radius
            .max(b.radius)
            .max((a.baseline - b.baseline).abs())
            .max(1e-9);
        limit.bound = if limit.sense == ResponseLimitSense3::AtMost {
            a.worst_value.max(b.worst_value) + 2.0 * scale
        } else {
            a.worst_value.min(b.worst_value) - 2.0 * scale
        };
        limit.scale = scale;
    }
    let p = spec.affine_grid_assessment.as_mut().unwrap();
    p.compliance_tolerance = 1e6;
    p.response_tolerance = 1e6;
    p.stress_tolerance = 1e6;
    p.inventory_tolerance = 1e6;
    spec
}

fn fine_operator(spec: &Spec) -> AdaptiveElasticity3 {
    let tree = Octree3::uniform((spec.level + 1) as u8, spec.max_level as u8, spec.leaves).unwrap();
    let domain = PhysicalDomain {
        bounds: spec.bounds,
        height: spec.height,
        curvature: spec.curvature,
    };
    let mut poll = |_| ControlFlow::Continue(());
    let mut q = QuadratureControl3::new(
        QuadratureOptions3 {
            max_boxes: spec.boxes,
            max_points: spec.points,
            ..Default::default()
        },
        &mut poll,
    )
    .unwrap();
    loading::build_operator(
        spec,
        HexCell::try_new(spec.bounds.0, spec.bounds.1).unwrap(),
        &tree,
        &domain,
        &loading::material(spec).unwrap(),
        &mut q,
    )
    .unwrap()
}

fn attach(
    spec: &Spec,
    run: &mut Computation,
    budget: SolveBudget,
    mut callback: impl FnMut(SolveProgress) -> ControlFlow<()>,
    gate: &CancelGate,
) {
    let tree = Octree3::uniform(spec.level as u8, spec.max_level as u8, spec.leaves).unwrap();
    let mut qp = |_| ControlFlow::Continue(());
    let mut quadrature = QuadratureControl3::new(
        QuadratureOptions3 {
            max_boxes: spec.boxes - run.state.spent.geometry.boxes,
            max_points: spec.points - run.state.spent.geometry.points,
            ..Default::default()
        },
        &mut qp,
    )
    .unwrap();
    let mut control = SolveControl::new(budget, &mut callback);
    grid_assessment::attach(
        spec.affine_grid_assessment.unwrap(),
        spec,
        &tree,
        &mut run.study,
        &mut run.state,
        &mut quadrature,
        &mut control,
        gate,
        || {
            if gate.is_requested() {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        },
    );
    let outcome = run.state.assessment.as_ref().unwrap();
    assert_eq!(outcome.work, control.work());
    assert_eq!(outcome.geometry, quadrature.work());
}

fn report(ledger: &Ledger, retained: &Loaded) -> JsonValue {
    let bytes = linked(ledger, &retained.value, "report_json", "study-report-json").unwrap();
    JsonValue::parse(std::str::from_utf8(&bytes).unwrap()).unwrap()
}

#[test]
fn g1_affine_grid_assessment_rebuilds_full_models_and_keeps_fine_physical_caps_strict() {
    let spec = fixture();
    let mut prefix = None;
    let initial = compute_observed(
        &spec,
        &CancelGate::new_clock_free(),
        0,
        None,
        0.0,
        |_, state| {
            assert!(state.assessment.is_none());
            prefix = Some(state.spent);
            Ok(())
        },
    )
    .unwrap();
    let r = checked(&initial);
    let a = initial.state.assessment.as_ref().unwrap();
    assert_eq!(r.rho, initial.state.selected().rho);
    assert_eq!(
        r.constraints,
        spec.affine_load_space.unwrap().constraint_policy()
    );
    assert_eq!(r.goal, Some(GOAL));
    assert!(a.fine_cells > a.coarse_cells && a.fine_dofs > a.coarse_dofs);
    assert!(a.work.linear_iterations > 0);
    assert!(a.work.linear_solves <= 1 + 2 * (spec.loads.len() + 1));
    let prefix = prefix.unwrap();
    assert_eq!(
        initial.state.spent.linear.linear_iterations - prefix.linear.linear_iterations,
        a.work.linear_iterations
    );
    assert_eq!(
        initial.state.spent.linear.linear_solves - prefix.linear.linear_solves,
        a.work.linear_solves
    );
    assert_eq!(
        initial.state.spent.geometry.points - prefix.geometry.points,
        a.geometry.points
    );
    near(
        r.normalized_model_change,
        difference(&r.coarse.gram, &r.fine.gram)
            + 2.0 * difference(&r.coarse.linear, &r.fine.linear)
            + (r.fine.offset - r.coarse.offset).abs(),
    );
    assert!(r.normalized_model_change > 0.0);
    near(
        r.normalized_maximum_change,
        r.fine.maximum - r.coarse.maximum,
    );
    near(
        r.normalized_witness_change,
        (r.fine.worst_compliance - r.coarse.worst_compliance)
            / spec.affine_load_space.unwrap().cap_j,
    );
    near(
        r.coarse.worst_compliance,
        initial.state.selected().evaluation.worst_compliance,
    );

    let mut fine = fine_operator(&spec);
    let parents = AdaptiveTransfer3::new(
        initial.study.operator().elasticity(),
        &fine,
        spec.affine_grid_assessment.unwrap().max_transfer_terms,
        || ControlFlow::Continue(()),
    )
    .unwrap()
    .parents()
    .to_vec();
    for (child, parent) in parents.into_iter().enumerate() {
        assert_eq!(r.fine_eroded_rho[child], r.coarse_eroded_rho[parent]);
        assert_eq!(r.fine_eroded_scales[child], r.coarse_eroded_scales[parent]);
    }
    fine.set_scales(&r.fine_eroded_scales).unwrap();
    let responses = r.responses.as_ref().unwrap();
    let prepared = response_binding::assemble(&spec, &fine, || ControlFlow::Continue(())).unwrap();
    response_binding::with_limits(prepared.as_ref(), |family| {
        for (limit, witness) in family.unwrap().limits.iter().zip(&responses.fine.witnesses) {
            let measure = |u: &[f64]| match limit.response {
                BoundedResponse3::Displacement(q) => dot(q, u),
                BoundedResponse3::Reaction(mode) => {
                    fine.embedded_reaction(u, None, mode, || ControlFlow::Continue(()))
                        .unwrap()
                        .value
                }
            };
            near(witness.baseline, measure(&r.fine.offset_displacement));
            for (u, value) in r
                .fine
                .basis_displacements
                .iter()
                .zip(&witness.basis_responses)
            {
                near_to_radius(*value, measure(u), witness.radius);
            }
        }
    });
    let volumes = fine.volumes();
    for (index, declared) in spec
        .affine_responses
        .as_ref()
        .unwrap()
        .limits
        .iter()
        .enumerate()
    {
        match declared.observation {
            Observation::MeanDisplacement { x_fraction, .. } => {
                let expected: f64 = fine
                    .leaves()
                    .iter()
                    .zip(&volumes)
                    .filter_map(|(leaf, volume)| {
                        let left = f64::from(leaf.index()[0]) / f64::from(1_u32 << leaf.level());
                        (x_fraction[0] <= left && left < x_fraction[1]).then_some(*volume)
                    })
                    .sum();
                near(a.fine_mean_volumes[index].unwrap(), expected);
                assert!(a.coarse_mean_volumes[index].unwrap() > 0.0);
            }
            _ => {
                assert_eq!(a.coarse_mean_volumes[index], None);
                assert_eq!(a.fine_mean_volumes[index], None);
            }
        }
        let coarse = &responses.coarse.witnesses[index];
        let fine = &responses.fine.witnesses[index];
        near(
            responses.normalized_model_changes[index],
            ((fine.baseline - coarse.baseline).abs()
                + difference(&coarse.basis_responses, &fine.basis_responses))
                / coarse.scale,
        );
    }
    let stress = r.stress.as_ref().unwrap();
    assert_eq!(stress.changes.len(), stress.fine.point_count);
    assert!(stress.coarse.adjoints.is_empty() && stress.fine.adjoints.is_empty());
    near(
        stress.normalized_aggregate_change,
        (stress.fine.aggregate - stress.coarse.aggregate).abs() / stress.fine.limit.limit,
    );

    let loose_spec = loose(spec, r);
    let admitted = compute(&loose_spec);
    assert!(
        admitted
            .state
            .assessment
            .as_ref()
            .unwrap()
            .passed(loose_spec.affine_grid_assessment.unwrap(), loose_spec.affine_load_space.unwrap().goal)
    );
    let mut strict_stress = loose_spec.clone();
    strict_stress
        .affine_load_space
        .as_mut()
        .unwrap()
        .stress
        .as_mut()
        .unwrap()
        .limit = 0.5 * stress.fine.sampled_physical_max;
    let rejected = compute(&strict_stress);
    assert!(
        !checked(&rejected)
            .stress
            .as_ref()
            .unwrap()
            .fine
            .physically_feasible()
    );
    assert!(
        !rejected
            .state
            .assessment
            .as_ref()
            .unwrap()
            .passed(strict_stress.affine_grid_assessment.unwrap(), strict_stress.affine_load_space.unwrap().goal)
    );
    let mut strict_response = loose_spec.clone();
    let witness = &responses.fine.witnesses[0];
    let limit = &mut strict_response.affine_responses.as_mut().unwrap().limits[0];
    limit.bound = witness.baseline + sign(limit.sense) * 0.5 * witness.radius;
    limit.scale = witness.radius;
    let rejected = compute(&strict_response);
    let checked = checked(&rejected);
    near(
        checked.responses.as_ref().unwrap().fine.witnesses[0].residual,
        0.5,
    );
    assert!(checked.stress.as_ref().unwrap().fine.physically_feasible());
    assert!(
        !rejected
            .state
            .assessment
            .as_ref()
            .unwrap()
            .passed(strict_response.affine_grid_assessment.unwrap(), strict_response.affine_load_space.unwrap().goal)
    );
}

#[test]
fn g3_g4_affine_grid_assessment_selects_the_incumbent_and_retains_partial_work_on_stops() {
    let baseline = compute(&fixture());
    let spec = loose(fixture(), checked(&baseline));
    let incumbent = coarse(&spec).state.accepted;
    assert!(spec.affine_load_space.unwrap().feasible(&incumbent));
    let mut dense = spec.clone();
    dense.density = 0.85;
    let mut run = coarse(&dense);
    assert_ne!(run.state.accepted.rho, incumbent.rho);
    assert!(
        incumbent.evaluation.physical.dilated_volume_fraction
            < run
                .state
                .accepted
                .evaluation
                .physical
                .dilated_volume_fraction
    );
    run.state.best = Some(incumbent);
    let accepted = run.state.accepted.rho.clone();
    let selected = run.state.selected().rho.clone();
    let checkpoint = run.state.checkpoint.optimizer.clone();
    let history = run.state.history.clone();
    let scales = run.study.operator().elasticity().scales().to_vec();
    let gate = CancelGate::new_clock_free();
    attach(
        &spec,
        &mut run,
        SolveBudget::default(),
        |_| ControlFlow::Continue(()),
        &gate,
    );
    assert_eq!(checked(&run).rho, selected);
    assert_ne!(checked(&run).rho, accepted);
    assert_eq!(run.state.checkpoint.optimizer, checkpoint);
    assert_eq!(run.state.history, history);

    // These independent attempts start from the same retained optimizer slice.
    // A one-iteration remaining allowance cannot be renewed by the fine probe.
    attach(
        &spec,
        &mut run,
        SolveBudget {
            total_iterations: 1,
            per_solve_iterations: spec.per_solve,
        },
        |_| ControlFlow::Continue(()),
        &gate,
    );
    let stopped = run.state.assessment.as_ref().unwrap();
    assert_eq!(stopped.status, "budget-exhausted", "{:?}", stopped.reason);
    assert!(stopped.report.is_none());
    assert_eq!(stopped.work.linear_iterations, 1);
    assert!(stopped.geometry.points > 0);
    assert_eq!(run.study.operator().elasticity().scales(), scales);

    let gate = CancelGate::new_clock_free();
    attach(
        &spec,
        &mut run,
        SolveBudget::default(),
        |progress| {
            if progress.stage == "sdf3-affine-goal-publish" {
                gate.request();
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        },
        &gate,
    );
    let stopped = run.state.assessment.as_ref().unwrap();
    assert_eq!(stopped.status, "cancelled");
    assert!(stopped.report.is_none());
    assert!(stopped.work.linear_iterations > 0 && stopped.geometry.points > 0);
    assert_eq!(run.state.accepted.rho, accepted);
    assert_eq!(run.state.selected().rho, selected);
    assert_eq!(run.state.checkpoint.optimizer, checkpoint);
    assert_eq!(run.state.history, history);
    assert_eq!(run.study.operator().elasticity().scales(), scales);

    // Calibrate the coarse map's actual term allowance through the native
    // admission boundary; no guessed stencil size or Debug-string parsing.
    let law = spec.material_acceleration.unwrap();
    let operator = run.study.operator().elasticity();
    let prepare = |terms| {
        operator.prepare_material_body_load(&|_| law.acceleration_m_s2, terms, || {
            ControlFlow::Continue(())
        })
    };
    let (mut lower, mut upper) = (1, law.max_terms);
    assert!(prepare(upper).is_ok());
    while lower < upper {
        let midpoint = lower + (upper - lower) / 2;
        match prepare(midpoint) {
            Ok(_) => upper = midpoint,
            Err(ElasticityError3::Invalid("material body-load term allowance exhausted")) => {
                lower = midpoint + 1
            }
            other => panic!("unexpected native body-map admission: {other:?}"),
        }
    }
    let fine = fine_operator(&spec);
    assert!(matches!(
        fine.prepare_material_body_load(&|_| law.acceleration_m_s2, lower, || {
            ControlFlow::Continue(())
        }),
        Err(ElasticityError3::Invalid(
            "material body-load term allowance exhausted"
        ))
    ));
    let mut limited = spec.clone();
    limited.material_acceleration.as_mut().unwrap().max_terms = lower;
    let bounded = compute(&limited);
    let stopped = bounded.state.assessment.as_ref().unwrap();
    assert_eq!(stopped.status, "budget-exhausted", "{:?}", stopped.reason);
    assert!(
        stopped
            .reason
            .as_ref()
            .unwrap()
            .contains("sdf3-material-body-load-budget")
    );
    assert!(stopped.report.is_none());
    assert!(stopped.geometry.points > 0 && stopped.work.linear_solves > 0);
    assert!(bounded.state.gradient.passed());
    assert_eq!(bounded.state.optimizer.iterations, 0);
    assert_eq!(
        bounded.state.accepted.rho,
        bounded.state.checkpoint.optimizer.point
    );
}

#[test]
fn g3_g5_affine_resume_finishes_a_stopped_grid_probe_without_new_optimizer_steps() {
    let spec = spec::parse(&FIXTURE.replace(":max-updates 80", ":max-updates 1")).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let gate = CancelGate::new_clock_free();
    let producer = checkpoint::producer_identity(&gate).unwrap();
    let mut previous = None;
    let mut pending = None;
    let run = compute_observed(&spec, &gate, 1, None, 0.0, |study, state| {
        let retained = output::persist(&spec, &ledger, study, state, producer, previous, 1)?;
        let retained = load(&ledger, &retained.pointer)?;
        previous = Some(retained.hash);
        if state.optimizer.iterations == spec.updates && state.optimization_status != "checkpointed"
        {
            pending = Some(retained);
            gate.request();
        }
        Ok(())
    })
    .unwrap();
    assert_eq!(run.state.optimizer.iterations, 1);
    assert_eq!(run.state.assessment.as_ref().unwrap().status, "cancelled");
    assert_ne!(run.state.optimization_status, "cancelled");
    let pending = pending.unwrap();
    let pending_report = report(&ledger, &pending);
    assert_eq!(
        pending_report
            .get("grid_assessment")
            .unwrap()
            .str_field("status"),
        Some("pending")
    );
    let stopped = output::persist(
        &spec, &ledger, &run.study, &run.state, producer, previous, 1,
    )
    .unwrap();
    let old = load(&ledger, &stopped.pointer).unwrap();
    let before = report(&ledger, &old);
    let resumed = resume(&ledger, &old, Some(100), &CancelGate::new_clock_free()).unwrap();
    let resumed = load(&ledger, &resumed.pointer).unwrap();
    let after = report(&ledger, &resumed);
    assert_eq!(integer(&resumed.value, "iterations_completed").unwrap(), 1);
    assert_eq!(
        after.str_field("optimization_status"),
        before.str_field("optimization_status")
    );
    assert_eq!(after.get("gradient_check"), before.get("gradient_check"));
    assert_eq!(
        linked(&ledger, &resumed.value, "iterations", "study-iterations").unwrap(),
        linked(&ledger, &old.value, "iterations", "study-iterations").unwrap()
    );
    let assessment = after.get("grid_assessment").unwrap();
    assert_eq!(assessment.str_field("status"), Some("checked"));
    assert_eq!(
        assessment.str_field("assessed_design"),
        after.str_field("selected_design")
    );
    assert!(
        assessment
            .get("result")
            .unwrap()
            .get("responses")
            .unwrap()
            .get("fine")
            .is_some()
    );
    assert!(
        assessment
            .get("result")
            .unwrap()
            .get("strength")
            .unwrap()
            .get("fine")
            .is_some()
    );
    let old_optimizer = old.value.get("optimizer_work").unwrap();
    let new_optimizer = resumed.value.get("optimizer_work").unwrap();
    let restorations = integer(new_optimizer, "restoration_evaluations").unwrap()
        - integer(old_optimizer, "restoration_evaluations").unwrap();
    assert!((1..=2).contains(&restorations));
    assert_eq!(
        integer(new_optimizer, "evaluations").unwrap(),
        integer(old_optimizer, "evaluations").unwrap() + restorations
    );
    let old_work = Spent::read(&old.value).unwrap();
    let new_work = Spent::read(&resumed.value).unwrap();
    let probe_work = assessment.get("work").unwrap();
    assert!(
        new_work.linear.linear_iterations
            >= old_work.linear.linear_iterations
                + integer(probe_work, "linear_iterations").unwrap()
    );
    assert!(new_work.geometry.points > old_work.geometry.points);
    new_work.validate(&spec).unwrap();
    let unchanged = resume(&ledger, &resumed, Some(100), &CancelGate::new_clock_free()).unwrap();
    assert_eq!(
        unchanged.pointer,
        format!("study-{}", resumed.hash.to_hex())
    );
    assert_eq!(unchanged.receipt, resumed.bytes);
    assert_eq!(load(&ledger, &stopped.pointer).unwrap().bytes, old.bytes);
    assert_eq!(
        load(&ledger, &format!("study-{}", pending.hash.to_hex()))
            .unwrap()
            .bytes,
        pending.bytes
    );
}
