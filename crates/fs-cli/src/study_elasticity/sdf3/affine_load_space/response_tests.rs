use super::*;
use fs_topopt::sdf3::response::limits::{BoundedResponse3, ResponseLimitSense3};
use spec::responses::Observation;

const FIXTURE: &str = include_str!(concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../examples/marquee/bracket-3d-affine-responses.fsim"
));

fn fixture() -> Spec {
    spec::parse(FIXTURE).unwrap()
}

fn compute(spec: &Spec, updates: usize) -> Computation {
    compute_observed(
        spec,
        &CancelGate::new_clock_free(),
        updates,
        None,
        0.0,
        |_, _| Ok(()),
    )
    .unwrap()
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

fn sign(sense: ResponseLimitSense3) -> f64 {
    match sense {
        ResponseLimitSense3::AtMost => 1.0,
        ResponseLimitSense3::AtLeast => -1.0,
    }
}

fn native_operator(spec: &Spec) -> AdaptiveElasticity3 {
    assert!(spec.physical);
    let tree = Octree3::uniform(spec.level as u8, spec.max_level as u8, spec.leaves).unwrap();
    let domain = PhysicalDomain {
        bounds: spec.bounds,
        height: spec.height,
        curvature: spec.curvature,
    };
    let mut poll = |_| ControlFlow::Continue(());
    let mut quadrature = QuadratureControl3::new(
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
        &mut quadrature,
    )
    .unwrap()
}

fn measure(operator: &AdaptiveElasticity3, response: BoundedResponse3<'_>, u: &[f64]) -> f64 {
    match response {
        BoundedResponse3::Displacement(q) => dot(q, u),
        BoundedResponse3::Reaction(mode) => {
            operator
                .embedded_reaction(u, None, mode, || ControlFlow::Continue(()))
                .unwrap()
                .value
        }
    }
}

fn report(ledger: &Ledger, retained: &Loaded) -> JsonValue {
    let bytes = linked(ledger, &retained.value, "report_json", "study-report-json").unwrap();
    JsonValue::parse(std::str::from_utf8(&bytes).unwrap()).unwrap()
}

#[test]
fn g1_affine_response_windows_use_native_displacement_force_and_moment_witnesses() {
    let spec = fixture();
    let run = compute(&spec, 0);
    let point = &run.state.accepted;
    let e = &point.evaluation;
    let response = e.responses.as_ref().unwrap();
    let authored = spec.affine_responses.as_ref().unwrap();
    assert!(run.state.gradient.passed());
    assert_eq!(
        run.state.gradient.constraints,
        spec.affine_load_space.unwrap().constraint_policy()
    );
    assert!(
        response
            .gradient
            .as_ref()
            .is_some_and(|g| g.len() == point.rho.len())
    );
    assert!(e.stress.as_ref().unwrap().gradient_admissible());
    assert_eq!(response.witnesses.len(), authored.limits.len());
    assert!(response.adjoints.len() <= authored.limits.len().min(spec.loads.len() + 1));
    near(
        run.state.gradient.baseline_constraint,
        point.constraint(GOAL).unwrap(),
    );
    assert!(point.constraint(GOAL).unwrap() >= response.constraint);

    // Rebuild the declared native operator, then use the actual retained eroded
    // material. Reaction observations cannot use the restored unit-stiffness K.
    let mut operator = native_operator(&spec);
    operator.set_scales(&e.physical.eroded_scales).unwrap();
    let prepared =
        response_binding::assemble(&spec, &operator, || ControlFlow::Continue(())).unwrap();
    let mut independent_adversary = false;
    response_binding::with_limits(prepared.as_ref(), |family| {
        let family = family.unwrap();
        for ((limit, declared), witness) in family
            .limits
            .iter()
            .zip(&authored.limits)
            .zip(&response.witnesses)
        {
            if let (Observation::MeanDisplacement { axis, .. }, BoundedResponse3::Displacement(q)) =
                (declared.observation, limit.response)
            {
                // Embedded support has no masked box nodes. A unit translation
                // must have mean one, not the unnormalized patch volume.
                assert!(operator.fixed().iter().all(|fixed| !fixed));
                near(q.chunks_exact(3).map(|node| node[axis]).sum(), 1.0);
                assert!(q.chunks_exact(3).all(|node| {
                    node.iter()
                        .enumerate()
                        .all(|(j, value)| j == axis || *value == 0.0)
                }));
            }
            let baseline = measure(&operator, limit.response, &e.offset_displacement);
            let basis: Vec<_> = e
                .basis_displacements
                .iter()
                .map(|u| measure(&operator, limit.response, u))
                .collect();
            let radius = dot(&basis, &basis).sqrt();
            assert!(radius > 0.0);
            near(witness.baseline, baseline);
            near(witness.radius, radius);
            for ((actual, expected), coefficient) in witness
                .basis_responses
                .iter()
                .zip(&basis)
                .zip(&witness.coefficients)
            {
                near_to_radius(*actual, *expected, radius);
                // Coefficients are components of a unit vector (asserted below), so their scale is 1.
                near_to_radius(*coefficient, sign(limit.sense) * expected / radius, 1.0);
            }
            near(dot(&witness.coefficients, &witness.coefficients), 1.0);
            near(witness.worst_value, baseline + sign(limit.sense) * radius);
            near(
                witness.residual,
                sign(limit.sense) * (witness.worst_value - limit.bound) / limit.scale,
            );
            let combined: Vec<_> = (0..e.offset_displacement.len())
                .map(|i| {
                    e.offset_displacement[i]
                        + e.basis_displacements
                            .iter()
                            .zip(&witness.coefficients)
                            .map(|(u, coefficient)| u[i] * coefficient)
                            .sum::<f64>()
                })
                .collect();
            near(
                measure(&operator, limit.response, &combined),
                witness.worst_value,
            );
            if !matches!(declared.observation, Observation::MeanDisplacement { .. }) {
                // Check authored axes and moment handedness independently of
                // the binding helper's virtual-mode closure.
                let mode = |p: [f64; 3], _: [f64; 3]| match declared.observation {
                    Observation::SupportForce { axis } => {
                        std::array::from_fn(|i| if i == axis { 1.0 } else { 0.0 })
                    }
                    Observation::SupportMoment { axis, about_m } => {
                        let r: [f64; 3] = std::array::from_fn(|i| p[i] - about_m[i]);
                        match axis {
                            0 => [0.0, -r[2], r[1]],
                            1 => [r[2], 0.0, -r[0]],
                            2 => [-r[1], r[0], 0.0],
                            _ => unreachable!("validated axis"),
                        }
                    }
                    Observation::MeanDisplacement { .. } => unreachable!(),
                };
                let direct = operator
                    .embedded_reaction(&combined, None, &mode, || ControlFlow::Continue(()))
                    .unwrap();
                near(direct.value, witness.worst_value);
            }
            near(witness.witness_value, witness.worst_value);
            independent_adversary |= witness
                .coefficients
                .iter()
                .zip(&e.worst_coefficients)
                .any(|(response, compliance)| (response - compliance).abs() > 1e-4);
        }
    });
    assert!(
        independent_adversary,
        "response bounds need their own adverse load directions"
    );

    // Strong support removes displacement DOFs, not physical averaging volume.
    // Renormalizing by the surviving load-vector weights would give a wrong
    // mean whenever the observation patch includes clamped material.
    let mut clamped = spec.clone();
    clamped.fixed = FixedFace::Left;
    let family = clamped.affine_responses.as_mut().unwrap();
    family
        .limits
        .retain(|limit| matches!(limit.observation, Observation::MeanDisplacement { .. }));
    for limit in &mut family.limits {
        limit.observation = Observation::MeanDisplacement {
            axis: 2,
            x_fraction: [0.0, 1.0],
        };
    }
    let operator = native_operator(&clamped);
    let volume: f64 = operator.volumes().iter().sum();
    let unnormalized = operator
        .reference_load(ReferenceLoad3::body(&|_| [0.0, 0.0, 1.0]), || {
            ControlFlow::Continue(())
        })
        .unwrap();
    let prepared =
        response_binding::assemble(&clamped, &operator, || ControlFlow::Continue(())).unwrap();
    response_binding::with_limits(prepared.as_ref(), |family| {
        for limit in family.unwrap().limits {
            let BoundedResponse3::Displacement(q) = limit.response else {
                panic!("displacement mean")
            };
            for (actual, raw) in q.iter().zip(&unnormalized) {
                near(*actual, raw / volume);
            }
            assert!(q.chunks_exact(3).map(|node| node[2]).sum::<f64>() < 1.0);
        }
    });
}

#[test]
fn g1_g5_affine_response_bounds_remain_independent_of_compliance_and_strength_caps() {
    let mut spec = fixture();
    let initial = compute(&spec, 0);
    let e = &initial.state.accepted.evaluation;
    let responses = e.responses.as_ref().unwrap();
    let config = spec.affine_load_space.as_mut().unwrap();
    config.cap_j = 2.0 * e.worst_compliance;
    config.stress.as_mut().unwrap().limit = 2.0 * e.stress.as_ref().unwrap().aggregate;
    for (limit, witness) in spec
        .affine_responses
        .as_mut()
        .unwrap()
        .limits
        .iter_mut()
        .zip(&responses.witnesses)
    {
        assert!(witness.radius > 0.0);
        limit.bound = witness.baseline + sign(limit.sense) * 2.0 * witness.radius;
        limit.scale = witness.radius;
    }
    let loose = compute(&spec, 4);
    let best = loose
        .state
        .best
        .as_ref()
        .expect("all calibrated caps are initially feasible");
    assert!(spec.affine_load_space.unwrap().feasible(best));
    assert!(
        best.evaluation
            .responses
            .as_ref()
            .unwrap()
            .physically_feasible()
    );
    assert!(
        best.evaluation
            .stress
            .as_ref()
            .unwrap()
            .physically_feasible()
    );
    assert!(best.evaluation.physical.dilated_volume_fraction < e.physical.dilated_volume_fraction);
    assert!(loose.state.gradient.passed());

    for sense in [ResponseLimitSense3::AtMost, ResponseLimitSense3::AtLeast] {
        let mut strict = spec.clone();
        let limits = &mut strict.affine_responses.as_mut().unwrap().limits;
        let index = limits
            .iter()
            .position(|limit| {
                limit.sense == sense
                    && matches!(limit.observation, Observation::MeanDisplacement { .. })
            })
            .unwrap();
        let baseline = &responses.witnesses[index];
        limits[index].bound = baseline.baseline + sign(sense) * 0.5 * baseline.radius;
        let rejected = compute(&strict, 0);
        let point = &rejected.state.accepted;
        let evaluation = &point.evaluation;
        let response = evaluation.responses.as_ref().unwrap();
        near(response.witnesses[index].residual, 0.5);
        assert!(!response.physically_feasible());
        assert!(evaluation.worst_compliance < strict.affine_load_space.unwrap().cap_j);
        assert!(evaluation.stress.as_ref().unwrap().physically_feasible());
        assert!(point.constraint(GOAL).unwrap() >= response.constraint);
        assert!(point.constraint(GOAL).unwrap() > 0.0);
        assert!(!strict.affine_load_space.unwrap().feasible(point));
        assert!(rejected.state.best.is_none());
        assert!(rejected.state.gradient.passed());
        assert_ne!(rejected.state.status, "completed");
    }
}

#[test]
fn g3_g5_affine_response_resume_preserves_joint_policy_and_original_signed_evidence() {
    let spec = fixture();
    let ledger = Ledger::open(":memory:").unwrap();
    let gate = CancelGate::new_clock_free();
    let prefix = drive(&spec, &ledger, Some(0), &gate).unwrap();
    let old = load(&ledger, &prefix.pointer).unwrap();
    let before = report(&ledger, &old);
    let resumed = resume(&ledger, &old, Some(1), &gate).unwrap();
    let resumed = load(&ledger, &resumed.pointer).unwrap();
    let after = report(&ledger, &resumed);
    for key in ["response_constraints", "strength_constraint"] {
        assert_eq!(
            before.get("model").unwrap().get(key),
            after.get("model").unwrap().get(key)
        );
    }
    assert_eq!(before.get("gradient_check"), after.get("gradient_check"));
    let policy = after
        .get("gradient_check")
        .unwrap()
        .get("constraint_policy")
        .unwrap();
    assert_eq!(policy.str_field("kind"), Some("joint-affine-constraints"));
    assert!(policy.get("responses").is_some());
    assert!(policy.get("strength").is_some());
    let responses = after.get("selected").unwrap().get("responses").unwrap();
    let witnesses = responses.get("limits").unwrap().as_array().unwrap();
    let declared = &spec.affine_responses.as_ref().unwrap().limits;
    assert_eq!(witnesses.len(), declared.len());
    for (witness, limit) in witnesses.iter().zip(declared) {
        assert_eq!(witness.str_field("name"), Some(limit.name.as_str()));
        assert_eq!(witness.str_field("kind"), Some(limit.observation.kind()));
        assert_eq!(witness.str_field("unit"), Some(limit.observation.unit()));
        near(witness.get("bound").unwrap().as_f64().unwrap(), limit.bound);
        assert_eq!(
            witness
                .get("coefficients")
                .unwrap()
                .as_array()
                .unwrap()
                .len(),
            spec.loads.len()
        );
    }
    assert_eq!(integer(&resumed.value, "iterations_completed").unwrap(), 1);
    assert_eq!(
        integer(
            resumed.value.get("optimizer_work").unwrap(),
            "restoration_evaluations"
        )
        .unwrap(),
        1
    );
    assert!(
        Spent::read(&resumed.value).unwrap().linear.linear_solves
            > Spent::read(&old.value).unwrap().linear.linear_solves
    );
    let design = linked(&ledger, &resumed.value, "design", "study-design").unwrap();
    let design = JsonValue::parse(std::str::from_utf8(&design).unwrap()).unwrap();
    assert_eq!(
        design.get("response_constrained"),
        Some(&JsonValue::Bool(true))
    );
    assert_eq!(
        design
            .get("selected")
            .unwrap()
            .get("responses")
            .unwrap()
            .get("witnesses")
            .unwrap()
            .as_array()
            .unwrap()
            .len(),
        declared.len()
    );

    for source in [
        FIXTURE.replace(
            ":name \"tip-z-lower\" :axis z",
            ":name \"tip-z-lower\" :axis y",
        ),
        FIXTURE.replace(":bound-m -0.00001953125", ":bound-m -0.000029296875"),
        FIXTURE.replace(":penalty 32.0", ":penalty 64.0"),
        FIXTURE.replace(":about-m (0.0 0.0 0.0)", ":about-m (0.25 0.0 0.0)"),
    ] {
        let changed = spec::parse(&source).unwrap();
        assert_ne!(changed.id, spec.id);
        let error = resume::drive(&changed, &ledger, Some(1), &gate, Some(&old)).unwrap_err();
        assert_eq!(error.code, "cli-study-sdf3-affine-load-space-checkpoint");
        assert_eq!(error.exit, exit::REFUSED);
    }
    assert_eq!(load(&ledger, &prefix.pointer).unwrap().bytes, old.bytes);
}
