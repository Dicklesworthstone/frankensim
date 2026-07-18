//! Battery for geometry program synthesis (fs-shapeprog). Covers DSL
//! round-trip, SDF semantics, the load-bearing rewrite-safety property (SDF
//! preserved within the declared certificate, checked by sampling),
//! canonicalization/dedup, and seeded shape-grammar derivation.

use fs_shapeprog::{
    BoundOperation, Certificate, Geom, ParseError, RewritePathStep, Simplified, SimplifyRefusal,
    linear_repeat, max_sdf_discrepancy, parse, simplify, stochastic_repeat,
};

fn grid() -> Vec<[f64; 3]> {
    let mut pts = Vec::new();
    for i in -3..=3 {
        for j in -3..=3 {
            for k in -3..=3 {
                pts.push([f64::from(i), f64::from(j), f64::from(k)]);
            }
        }
    }
    pts
}

fn next_up(value: f64) -> f64 {
    assert!(value.is_finite() && value >= 0.0);
    if value == 0.0 {
        f64::from_bits(1)
    } else {
        f64::from_bits(value.to_bits() + 1)
    }
}

fn next_down(value: f64) -> f64 {
    assert!(value.is_finite() && value > 0.0);
    f64::from_bits(value.to_bits() - 1)
}

fn assert_certificate_holds(
    label: &str,
    original: &Geom,
    simplified: &Simplified,
    samples: &[[f64; 3]],
) -> f64 {
    assert!(
        !simplified.is_refused(),
        "{label}: unexpected transactional refusal: {simplified:#?}"
    );
    assert!(
        simplified.max_error.is_finite() && simplified.max_error >= 0.0,
        "{label}: successful certificate must be finite and nonnegative: {simplified:#?}"
    );
    let actual = max_sdf_discrepancy(original, &simplified.program, samples);
    assert!(
        actual <= simplified.max_error,
        "{label}: actual discrepancy {actual:e} exceeds bound {:e}; result={simplified:#?}; original={original:#?}",
        simplified.max_error
    );
    actual
}

#[test]
fn the_dsl_round_trips() {
    let g = Geom::sphere(2.0)
        .offset(1.0)
        .translate([1.0, 2.0, 3.0])
        .union(Geom::cube(3.0));
    assert_eq!(parse(&g.to_sexpr()), Ok(g.clone()));
    // print/parse/print is stable.
    assert_eq!(parse(&g.to_sexpr()).unwrap().to_sexpr(), g.to_sexpr());
}

#[test]
fn sdf_semantics_are_correct() {
    let s = Geom::sphere(2.0);
    assert!((s.sdf([0.0, 0.0, 0.0]) + 2.0).abs() < 1e-12); // centre: -radius
    assert!((s.sdf([3.0, 0.0, 0.0]) - 1.0).abs() < 1e-12); // 3-2 = 1
    // union is the min of the two SDFs.
    let u = Geom::sphere(1.0).union(Geom::sphere(1.0).translate([5.0, 0.0, 0.0]));
    assert!((u.sdf([5.0, 0.0, 0.0]) + 1.0).abs() < 1e-12);
}

#[test]
fn consecutive_offsets_preserve_the_floating_evaluation_order() {
    // The real-arithmetic identity would merge these radii, but the interpreter
    // performs two rounded subtractions. Without a bit-equivalence proof the
    // only exact simplification is to preserve the two nodes.
    let pre = Geom::sphere(2.0).offset(1.0).offset(2.0);
    let out = simplify(&pre, 1e-6);
    assert_eq!(out.program, pre);
    assert!(
        out.rewrites.is_empty(),
        "no reassociation is authorized: {out:#?}"
    );
    assert_eq!(out.max_error, 0.0);
    assert_certificate_holds("preserved offset order", &pre, &out, &grid());
}

#[test]
fn identity_and_distribution_rewrites_hold() {
    // union(a, empty) -> a.
    let pre = Geom::sphere(2.0).union(Geom::Empty);
    let out = simplify(&pre, 1e-6);
    assert_eq!(out.program, Geom::sphere(2.0));
    assert!(max_sdf_discrepancy(&pre, &out.program, &grid()) < 1e-9);
    // translate(union(a,b), t) -> union(translate a, translate b).
    let pre2 = Geom::sphere(1.0)
        .union(Geom::cube(1.0))
        .translate([2.0, 0.0, 0.0]);
    let out2 = simplify(&pre2, 1e-6);
    assert!(
        out2.rewrites
            .iter()
            .any(|r| r.rule == "translate-distributes")
    );
    assert!(max_sdf_discrepancy(&pre2, &out2.program, &grid()) < 1e-9);
}

#[test]
fn a_certified_approximate_rewrite_stays_within_its_bound() {
    // Correctly rounded subtraction can move to the neighbouring lattice point
    // by more than |r|, so the context-free global envelope is 2|r|.
    let pre = Geom::sphere(2.0).offset(0.001);
    let out = simplify(&pre, 0.01);
    assert_eq!(out.program, Geom::sphere(2.0));
    let bound = out
        .rewrites
        .iter()
        .find(|r| r.rule == "drop-tiny-offset")
        .map(|r| r.certificate);
    assert!(matches!(bound, Some(Certificate::Approximate { .. })));
    assert_eq!(out.max_error.to_bits(), 0.002_f64.to_bits());
    assert_certificate_holds("single tiny offset", &pre, &out, &grid());
}

#[test]
fn sequential_tiny_offsets_add_their_outward_bounds() {
    // G0 regression: sequential lossy rewrites lie on the same evaluation path
    // and therefore ADD. A global maximum would forge a certificate.
    let pre = Geom::sphere(1.0).offset(0.006).offset(0.006);
    let out = simplify(&pre, 0.01);
    assert_eq!(out.program, Geom::sphere(1.0));
    assert_eq!(out.max_error.to_bits(), 0.024_f64.to_bits());
    assert_eq!(out.first_lossy_rewrite, Some(0));
    assert_eq!(
        out.rewrites.len(),
        2,
        "both drops remain replayable: {out:#?}"
    );
    assert_eq!(out.rewrites[0].path, Vec::<RewritePathStep>::new());
    assert_eq!(out.rewrites[1].path, vec![RewritePathStep::OffsetChild]);
    assert_eq!(
        out.rewrites[0].accumulated_bound.to_bits(),
        out.max_error.to_bits(),
        "the outer trace entry must replay the complete sequential subtree"
    );
    assert_eq!(
        out.rewrites[1].accumulated_bound.to_bits(),
        0.012_f64.to_bits(),
        "the inner trace entry must retain its local subtree envelope"
    );
    assert_certificate_holds("sequential offset stack", &pre, &out, &grid());
}

#[test]
fn sequential_addition_rounds_outward_when_a_tiny_term_is_absorbed() {
    // The second local envelope is far below one ulp of the first. Ordinary
    // round-to-nearest addition would silently return 0.5; the certificate must
    // advance to the next representable upper bound instead.
    let min_subnormal = f64::from_bits(1);
    let pre = Geom::sphere(1.0).offset(min_subnormal).offset(0.25);
    let out = simplify(&pre, 1.0);
    assert_eq!(out.program, Geom::sphere(1.0));
    assert_eq!(out.max_error.to_bits(), next_up(0.5).to_bits());
    assert_certificate_holds("outward sequential addition", &pre, &out, &grid());
}

#[test]
fn sequential_loss_across_a_boolean_constructor_is_not_branch_loss() {
    // Mandatory G0 witness from frankensim-shapeprog-compositional-error-rhn9w:
    // the inner and outer 0.006 offsets are sequential on the active origin
    // branch even though a Union separates their AST nodes. The historical
    // global-max report was 0.006 while the actual discrepancy is about 0.012.
    let pre = Geom::sphere(1.0)
        .offset(0.006)
        .union(Geom::sphere(1.0).translate([100.0, 0.0, 0.0]))
        .offset(0.006);
    let out = simplify(&pre, 0.01);
    let origin = [[0.0, 0.0, 0.0]];
    let actual = assert_certificate_holds("sequential union witness", &pre, &out, &origin);

    assert!(
        actual > 0.006,
        "witness must exercise the old false certificate, actual={actual:e}; {out:#?}"
    );
    assert_eq!(out.max_error.to_bits(), 0.024_f64.to_bits());
    assert_eq!(out.first_lossy_rewrite, Some(0));
    assert_eq!(out.rewrites[0].path, Vec::<RewritePathStep>::new());
    assert_eq!(
        out.rewrites[1].path,
        vec![RewritePathStep::OffsetChild, RewritePathStep::UnionLeft]
    );
    assert_eq!(
        out.pass_bounds.len(),
        1,
        "both sequential drops occur in one pass"
    );
    assert_eq!(out.pass_bounds[0].bound.to_bits(), out.max_error.to_bits());
}

#[test]
fn catastrophic_offset_cancellation_is_never_labeled_exact() {
    // Mandatory G0 witness: two rounded subtractions are not equivalent to one
    // subtraction by the rounded parameter sum. At the origin the reassociated
    // program differs by exactly one while the preserved program is bit-stable.
    for pre in [
        Geom::sphere(1.0).offset(f64::MAX).offset(-f64::MAX),
        Geom::sphere(1.0).offset(-f64::MAX).offset(f64::MAX),
    ] {
        let out = simplify(&pre, 0.01);
        assert_eq!(
            out.program, pre,
            "offset evaluation order must be preserved"
        );
        assert!(
            out.rewrites.is_empty(),
            "no offset-compose claim is legal: {out:#?}"
        );
        assert_eq!(out.max_error, 0.0);
        assert_eq!(
            out.program.sdf([0.0, 0.0, 0.0]).to_bits(),
            pre.sdf([0.0, 0.0, 0.0]).to_bits()
        );

        let reassociated = Geom::sphere(1.0).offset(0.0);
        assert_eq!(
            max_sdf_discrepancy(&pre, &reassociated, &[[0.0, 0.0, 0.0]]),
            1.0,
            "the exact witness must remain sensitive to reassociation"
        );
    }
}

#[test]
fn rounded_subtraction_threshold_requires_the_two_radius_envelope() {
    // Around x=1, the predecessor gap is 2^-53. A radius just above half that
    // gap crosses the rounding-cell boundary: the computed SDF moves by one
    // full predecessor gap, which is strictly greater than |radius|.
    let half_predecessor_gap = f64::from_bits((969_u64) << 52); // 2^-54
    let radius = next_up(half_predecessor_gap);
    let pre = Geom::sphere(0.0).offset(radius);
    let out = simplify(&pre, next_up(radius));
    let actual = assert_certificate_holds(
        "rounded subtraction threshold",
        &pre,
        &out,
        &[[1.0, 0.0, 0.0]],
    );

    assert!(
        actual > radius,
        "the witness must falsify the historical |radius| certificate: radius={radius:e}, actual={actual:e}"
    );
    assert_eq!(out.max_error.to_bits(), (radius * 2.0).to_bits());
}

#[test]
fn independent_boolean_branches_use_max_not_sum() {
    // G3: min/max/difference are 1-Lipschitz in the infinity norm over their
    // independently simplified branches. The two local envelopes are 0.006
    // and 0.008, so the root bound is 0.008 rather than their sum.
    let left = Geom::sphere(1.0).offset(0.003);
    let right = Geom::cube(1.0).offset(-0.004).translate([10.0, 0.0, 0.0]);
    let samples = [[0.0, 0.0, 0.0], [10.0, 0.0, 0.0]];

    let union = left.clone().union(right.clone());
    let union_out = simplify(&union, 0.01);
    assert_eq!(union_out.max_error.to_bits(), 0.008_f64.to_bits());
    assert_certificate_holds("union alternative", &union, &union_out, &samples);

    let intersect = Geom::Intersect(Box::new(left.clone()), Box::new(right.clone()));
    let intersect_out = simplify(&intersect, 0.01);
    assert_eq!(intersect_out.max_error.to_bits(), 0.008_f64.to_bits());
    assert_certificate_holds(
        "intersection alternative",
        &intersect,
        &intersect_out,
        &samples,
    );

    let difference = Geom::Difference(Box::new(left), Box::new(right));
    let difference_out = simplify(&difference, 0.01);
    assert_eq!(difference_out.max_error.to_bits(), 0.008_f64.to_bits());
    assert_certificate_holds(
        "difference alternative",
        &difference,
        &difference_out,
        &samples,
    );
}

#[test]
fn exact_identities_transport_or_discharge_child_error_by_context() {
    let lossy = Geom::sphere(1.0).offset(0.003);

    // Removing an Empty union alternative preserves the live child's error;
    // the exact wrapper rewrite must not reset that inherited envelope.
    let live_union = Geom::Union(Box::new(Geom::Empty), Box::new(lossy.clone()));
    let live_out = simplify(&live_union, 0.01);
    assert_eq!(live_out.program, Geom::sphere(1.0));
    assert_eq!(live_out.max_error.to_bits(), 0.006_f64.to_bits());
    assert!(live_out.rewrites.iter().any(|rewrite| {
        rewrite.rule == "union-identity"
            && rewrite.accumulated_bound.to_bits() == 0.006_f64.to_bits()
    }));
    assert_certificate_holds(
        "identity transports live child error",
        &live_union,
        &live_out,
        &grid(),
    );

    // An Empty intersection makes the other branch semantically irrelevant.
    // Skipping its lossy rewrite is exact and contributes no phantom error.
    let empty_intersection = Geom::Intersect(Box::new(Geom::Empty), Box::new(lossy));
    let empty_out = simplify(&empty_intersection, 0.01);
    assert_eq!(empty_out.program, Geom::Empty);
    assert_eq!(empty_out.max_error, 0.0);
    assert!(
        empty_out
            .rewrites
            .iter()
            .all(|rewrite| rewrite.certificate == Certificate::Exact),
        "irrelevant branch must not create a lossy certificate: {empty_out:#?}"
    );
    assert_certificate_holds(
        "identity discharges irrelevant child error",
        &empty_intersection,
        &empty_out,
        &grid(),
    );
}

#[test]
fn branch_permutation_preserves_the_bound_and_canonical_result() {
    // G3: commutative branch order changes deterministic paths but not the
    // theorem. Both independent alternatives have the same root max bound.
    let a = Geom::sphere(1.0).offset(0.003);
    let b = Geom::cube(2.0).offset(-0.004);
    let ab = simplify(&a.clone().union(b.clone()), 0.01);
    let ba = simplify(&b.union(a), 0.01);

    assert_eq!(ab.max_error.to_bits(), ba.max_error.to_bits());
    assert_eq!(ab.program.canonical(), ba.program.canonical());
    assert_eq!(ab.rewrites.len(), ba.rewrites.len());
    assert!(
        ab.rewrites
            .iter()
            .chain(&ba.rewrites)
            .all(|rewrite| matches!(rewrite.certificate, Certificate::Approximate { .. }))
    );
}

#[test]
fn unary_contexts_apply_their_declared_affine_and_rounding_envelopes() {
    let inner = Geom::sphere(1.0).offset(0.003);

    // Translation evaluates both programs at the same transformed point and
    // therefore carries the child bound with factor one.
    let translated = inner.clone().translate([2.0, -3.0, 4.0]);
    let translated_out = simplify(&translated, 0.01);
    assert_eq!(translated_out.max_error.to_bits(), 0.006_f64.to_bits());
    assert_certificate_holds(
        "translation factor one",
        &translated,
        &translated_out,
        &grid(),
    );

    // A retained rounded Offset is not globally Lipschitz on the floating
    // lattice. Its sound range-free envelope is E + 2|radius|: factor one for
    // the real affine map plus a nearest-rounding envelope at both endpoints.
    let retained_offset = inner.offset(1.0);
    let offset_out = simplify(&retained_offset, 0.01);
    assert!(
        offset_out.max_error >= 2.006,
        "retained offset must include both radius-rounding envelopes: {offset_out:#?}"
    );
    assert_certificate_holds(
        "retained offset affine plus rounding",
        &retained_offset,
        &offset_out,
        &grid(),
    );
}

#[test]
fn retained_offset_does_not_assume_a_false_global_lipschitz_factor() {
    // G0/G3 rounding-cell witness. The simplified child is just above the
    // midpoint between 1 and next_up(1); the original child is just below it.
    // Adding the same representable 1.0 therefore expands a ~2^-104 child
    // perturbation into a 2^-52 output jump. No constant factor-two child-bound
    // rule could certify this. The range-free E + 2|radius| envelope can.
    let midpoint = f64::from_bits((970_u64) << 52); // 2^-53
    let child_ulp = f64::from_bits((918_u64) << 52); // 2^-105
    let base = next_up(midpoint);
    let inner_radius = child_ulp * 2.0;
    let pre = Geom::sphere(-base).offset(inner_radius).offset(-1.0);
    let out = simplify(&pre, next_up(inner_radius));
    let actual = assert_certificate_holds(
        "retained-offset rounding-cell expansion",
        &pre,
        &out,
        &[[0.0, 0.0, 0.0]],
    );

    let expanded_gap = f64::from_bits((971_u64) << 52); // 2^-52
    assert_eq!(
        actual.to_bits(),
        expanded_gap.to_bits(),
        "witness must cross the binary64 rounding midpoint: {out:#?}"
    );
    assert!(
        actual > inner_radius * 4.0,
        "the witness must falsify child-bound factor-two propagation"
    );
    assert!(
        out.max_error > 2.0,
        "the sound range-free envelope must retain the parent rounding term"
    );
}

#[test]
fn threshold_neighbours_are_monotone_in_tolerance() {
    let radius = 0.01_f64;
    let tolerances = [0.0, next_down(radius), radius, next_up(radius), 1.0];
    let input = Geom::sphere(1.0).offset(radius);
    let outputs: Vec<_> = tolerances
        .iter()
        .map(|&tol| simplify(&input, tol))
        .collect();

    for window in outputs.windows(2) {
        assert!(
            window[1].program.size() <= window[0].program.size(),
            "larger tolerance must not make the program larger: {window:#?}"
        );
    }
    assert_eq!(outputs[0].program, input);
    assert_eq!(outputs[1].program, input);
    assert_eq!(
        outputs[2].program, input,
        "the admission predicate is strict"
    );
    assert_eq!(outputs[3].program, Geom::sphere(1.0));
    assert_eq!(outputs[4].program, Geom::sphere(1.0));
}

#[test]
fn signed_zero_and_subnormal_offset_bounds_are_finite_and_replayable() {
    let min_subnormal = f64::from_bits(1);
    let tolerance = f64::from_bits(2);
    for radius in [0.0, -0.0, min_subnormal, -min_subnormal] {
        let pre = Geom::sphere(1.0).offset(radius);
        let out = simplify(&pre, tolerance);
        assert_eq!(out.program, Geom::sphere(1.0));
        let expected = radius.abs() * 2.0;
        assert_eq!(
            out.max_error.to_bits(),
            expected.to_bits(),
            "signed/subnormal case radius={radius:?}: {out:#?}"
        );
        assert_certificate_holds("signed/subnormal drop", &pre, &out, &grid());
    }
}

#[test]
fn invalid_or_unrepresentable_bound_inputs_refuse_transactionally() {
    let ordinary = Geom::sphere(1.0).offset(0.001);
    for tolerance in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -1.0] {
        let out = simplify(&ordinary, tolerance);
        assert_eq!(out.program, ordinary);
        assert_eq!(out.max_error, 0.0);
        assert!(out.rewrites.is_empty());
        assert_eq!(
            out.refusal,
            Some(SimplifyRefusal::InvalidTolerance {
                tolerance_bits: tolerance.to_bits(),
            }),
            "invalid tolerance diagnostics must preserve exact IEEE bits"
        );
    }

    let negative_zero_tolerance = simplify(&ordinary, -0.0);
    assert!(
        !negative_zero_tolerance.is_refused(),
        "negative zero is a finite nonnegative threshold"
    );
    assert_eq!(negative_zero_tolerance.program, ordinary);

    let invalid_program = Geom::sphere(f64::NAN).offset(0.001);
    let invalid_out = simplify(&invalid_program, 0.01);
    assert_eq!(invalid_out.program.to_sexpr(), invalid_program.to_sexpr());
    assert_eq!(invalid_out.refusal, Some(SimplifyRefusal::NonFiniteProgram));

    // The radius is admitted by the tolerance, but its mandatory 2|r| local
    // envelope overflows. The transaction must return the untouched program,
    // not a partial simplification with a flattering finite bound.
    let huge_radius = f64::MAX / 1.5;
    assert!(huge_radius.is_finite() && huge_radius < f64::MAX);
    let overflow = Geom::sphere(1.0).offset(huge_radius);
    let overflow_out = simplify(&overflow, f64::MAX);
    assert_eq!(overflow_out.program, overflow);
    assert_eq!(overflow_out.max_error, 0.0);
    assert!(overflow_out.rewrites.is_empty());
    assert!(matches!(
        overflow_out.refusal,
        Some(SimplifyRefusal::UnrepresentableBound {
            pass: 0,
            ref path,
            operation: BoundOperation::LocalOffset,
        }) if path.is_empty()
    ));

    // Each local envelope is finite, but their sequential outward sum is not.
    let component_radius = f64::MAX / 3.0;
    let sequential_overflow = Geom::sphere(1.0)
        .offset(component_radius)
        .offset(component_radius);
    let sequential_out = simplify(&sequential_overflow, f64::MAX);
    assert_eq!(sequential_out.program, sequential_overflow);
    assert!(matches!(
        sequential_out.refusal,
        Some(SimplifyRefusal::UnrepresentableBound {
            pass: 0,
            ref path,
            operation: BoundOperation::Sequential,
        }) if path.is_empty()
    ));

    // The inner local envelope is finite, but a retained parent Offset's
    // mandatory 2|radius| rounding envelope is not. This is a scaling-context
    // refusal, not a local-drop refusal.
    let scale_overflow = Geom::sphere(1.0).offset(component_radius).offset(f64::MAX);
    let scale_out = simplify(&scale_overflow, f64::MAX);
    assert_eq!(scale_out.program, scale_overflow);
    assert!(matches!(
        scale_out.refusal,
        Some(SimplifyRefusal::UnrepresentableBound {
            pass: 0,
            ref path,
            operation: BoundOperation::Scale,
        }) if path.is_empty()
    ));
}

#[test]
fn successful_simplification_is_deterministic_and_idempotent() {
    let input = Geom::sphere(1.0)
        .offset(0.003)
        .union(Geom::cube(1.0).offset(-0.004))
        .translate([2.0, 0.0, 0.0]);
    let first = simplify(&input, 0.01);
    let replay = simplify(&input, 0.01);
    assert_eq!(
        first, replay,
        "rewrite trace and bounds must replay bit-for-bit"
    );
    assert_certificate_holds("deterministic first pass", &input, &first, &grid());

    let fixed = simplify(&first.program, 0.01);
    assert!(
        !fixed.is_refused(),
        "fixed point unexpectedly refused: {fixed:#?}"
    );
    assert_eq!(fixed.program, first.program);
    assert!(
        fixed.rewrites.is_empty(),
        "a true fixed point applies no rewrites"
    );
    assert!(fixed.pass_bounds.is_empty());
    assert_eq!(fixed.max_error, 0.0);
}

#[test]
fn pass_budget_exhaustion_refuses_instead_of_returning_a_partial_program() {
    // Translation distribution advances one union level per pass. A depth
    // beyond the deterministic budget must roll back rather than call the
    // partially distributed tree a fixed point.
    let mut deep = Geom::sphere(1.0);
    for _ in 0..66 {
        deep = deep.union(Geom::cube(1.0));
    }
    let input = deep.translate([1.0, 0.0, 0.0]);
    let out = simplify(&input, 0.0);
    assert_eq!(out.program, input);
    assert_eq!(out.max_error, 0.0);
    assert!(out.rewrites.is_empty());
    assert_eq!(
        out.refusal,
        Some(SimplifyRefusal::PassLimitExceeded { limit: 64 })
    );
}

#[test]
fn canonicalization_deduplicates_commutative_programs() {
    let a = Geom::sphere(1.0);
    let b = Geom::cube(2.0);
    let ab = a.clone().union(b.clone());
    let ba = b.union(a);
    // union is commutative -> the same content hash (archive/ledger dedup).
    assert_eq!(ab.canonical_hash(), ba.canonical_hash());
    assert_eq!(ab.canonical(), ba.canonical());
}

#[test]
fn grammar_derivations_are_reproducible_from_seeds() {
    let unit = Geom::sphere(0.4);
    let three = linear_repeat(&unit, 3, [1.0, 0.0, 0.0]);
    // three unit spheres -> each centre is inside.
    for i in 0..3 {
        assert!(three.sdf([f64::from(i), 0.0, 0.0]) < 0.0);
    }
    // the same seed yields the same derivation; different seeds may differ.
    let s1 = stochastic_repeat(&unit, 5, [1.0, 0.0, 0.0], 42);
    let s1b = stochastic_repeat(&unit, 5, [1.0, 0.0, 0.0], 42);
    assert_eq!(s1, s1b);
    assert_eq!(s1.canonical_hash(), s1b.canonical_hash());
}

#[test]
fn malformed_programs_are_rejected() {
    assert!(matches!(parse("("), Err(ParseError::UnexpectedEnd)));
    assert!(matches!(
        parse("(banana 3)"),
        Err(ParseError::Unexpected(_))
    ));
    assert!(matches!(parse("(sphere x)"), Err(ParseError::BadNumber(_))));
    // a truncated primitive reads no number -> unexpected end.
    assert!(matches!(parse("(sphere"), Err(ParseError::UnexpectedEnd)));
}

#[test]
fn non_finite_numeric_atoms_are_rejected_in_every_grammar_slot() {
    for atom in ["NaN", "inf", "-inf", "1e999", "-1e999"] {
        for source in [
            format!("(sphere {atom})"),
            format!("(cube {atom})"),
            format!("(offset (sphere 1) {atom})"),
            format!("(translate (sphere 1) {atom} 0 0)"),
            format!("(translate (sphere 1) 0 {atom} 0)"),
            format!("(translate (sphere 1) 0 0 {atom})"),
        ] {
            assert_eq!(
                parse(&source),
                Err(ParseError::BadNumber(atom.to_string())),
                "source {source} must fail closed"
            );
        }
    }

    let parsed = parse("(sphere -0.0)").expect("signed zero remains finite");
    let Geom::Primitive { size, .. } = parsed else {
        panic!("sphere parser must return a primitive");
    };
    assert_eq!(size.to_bits(), (-0.0_f64).to_bits());
}

#[test]
fn discrepancy_refuses_invalid_or_unrepresentable_evidence() {
    let one = Geom::sphere(1.0);
    let two = Geom::sphere(2.0);
    let sample = [[0.0, 0.0, 0.0]];
    assert_eq!(max_sdf_discrepancy(&one, &two, &sample), 1.0);
    assert!(max_sdf_discrepancy(&one, &two, &[]).is_infinite());

    for invalid_sample in [
        [[f64::NAN, 0.0, 0.0]],
        [[f64::INFINITY, 0.0, 0.0]],
        [[f64::NEG_INFINITY, 0.0, 0.0]],
    ] {
        assert!(max_sdf_discrepancy(&one, &one, &invalid_sample).is_infinite());
    }

    for invalid_program in [
        Geom::sphere(f64::NAN),
        Geom::cube(f64::INFINITY),
        Geom::sphere(1.0).offset(f64::NEG_INFINITY),
        Geom::sphere(1.0).translate([0.0, f64::NAN, 0.0]),
    ] {
        assert!(max_sdf_discrepancy(&invalid_program, &invalid_program, &sample).is_infinite());
    }

    let extreme_sample = [[f64::MAX, f64::MAX, f64::MAX]];
    assert!(max_sdf_discrepancy(&one, &one, &extreme_sample).is_infinite());
    assert!(
        max_sdf_discrepancy(&Geom::sphere(-f64::MAX), &Geom::sphere(f64::MAX), &sample)
            .is_infinite()
    );

    // The exact difference 1 + min_subnormal rounds to 1 in binary64. Evidence
    // must advance outward instead of silently understating that discrepancy.
    let min_subnormal = f64::from_bits(1);
    assert_eq!(
        max_sdf_discrepancy(&Geom::sphere(-1.0), &Geom::sphere(min_subnormal), &sample,).to_bits(),
        next_up(1.0).to_bits()
    );

    // Near cancellation between adjacent floats is exact by the subtraction
    // geometry: the evidence result is exactly one ulp, not a gratuitously
    // widened neighbour.
    let one_ulp_at_one = next_up(1.0) - 1.0;
    assert_eq!(
        max_sdf_discrepancy(&Geom::sphere(-next_up(1.0)), &Geom::sphere(-1.0), &sample,).to_bits(),
        one_ulp_at_one.to_bits()
    );

    // A finite Union root must not launder a non-finite inactive branch into
    // apparently valid evidence. The Boolean max-bound theorem requires every
    // non-structural branch value, not merely the selected root value, finite.
    let masked_non_finite =
        Geom::sphere(1.0).union(Geom::sphere(1.0).translate([f64::MAX, f64::MAX, f64::MAX]));
    assert_eq!(masked_non_finite.sdf([0.0, 0.0, 0.0]), -1.0);
    assert!(
        max_sdf_discrepancy(&masked_non_finite, &masked_non_finite, &[[0.0, 0.0, 0.0]],)
            .is_infinite(),
        "a finite selected root cannot certify a masked overflowing branch"
    );

    for structural_empty in [
        Geom::Empty,
        Geom::Union(Box::new(Geom::Empty), Box::new(Geom::Empty)),
        Geom::Intersect(Box::new(Geom::Empty), Box::new(one.clone())),
        Geom::Difference(Box::new(Geom::Empty), Box::new(one)),
        Geom::Empty.offset(1.0),
        Geom::Empty.translate([1.0, 2.0, 3.0]),
        Geom::Intersect(
            Box::new(Geom::Empty),
            Box::new(Geom::sphere(1.0).translate([f64::MAX, f64::MAX, f64::MAX])),
        ),
    ] {
        assert_eq!(
            max_sdf_discrepancy(&structural_empty, &Geom::Empty, &sample),
            0.0
        );
    }
}
