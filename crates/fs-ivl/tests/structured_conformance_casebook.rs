//! Cheap structured per-crate evidence for fs-ivl (bead 6ys.18.10).
//!
//! These fixed public-API fixtures register representative interval/affine,
//! Taylor/root, and exact-predicate laws with the shared Casebook harness.
//! They complement rather than replace the retained random-DAG/property,
//! Taylor-convergence/certified-root, adversarial-predicate, aggregate-golden,
//! performance, and dual-ISA batteries.

use core::fmt::Write as _;
use fs_casebook::{CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, ToleranceSpec, fnv1a64};
use fs_ivl::{
    AffineCtx, Interval, MAX_TAYLOR_ORDER, RootSearchConfig, RootSearchError, RootSearchReport,
    Sign, Stage, TaylorModel1, TaylorModelError, incircle, newton_roots_bounded, orient2d_sos,
    orient2d_with_stage,
};

const SUITE: &str = "fs-ivl/structured-conformance-v1";
const AFFINE_ABSOLUTE_TOLERANCE: f64 = 1.0e-13;
const PLAIN_MINIMUM_WIDTH: f64 = 1.0;
const WHOLE_BITS: [u64; 2] = [0xfff0_0000_0000_0000, 0x7ff0_0000_0000_0000];

const CCW_POINTS: [[f64; 2]; 3] = [[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]];
const COLLINEAR_POINTS: [[f64; 2]; 3] = [[0.0, 0.0], [1.0, 1.0], [2.0, 2.0]];
const TIE_POINT: [f64; 2] = [1.0, 1.0];
const CIRCLE_POINTS: [[f64; 2]; 4] = [[5.0, 0.0], [3.0, 4.0], [-3.0, 4.0], [0.0, -5.0]];
const PREDICATE_REFERENCES: [Sign; 6] = [
    Sign::Positive,
    Sign::Negative,
    Sign::Zero,
    Sign::Positive,
    Sign::Negative,
    Sign::Zero,
];

fn push_u64(bytes: &mut Vec<u8>, value: u64) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn push_len(bytes: &mut Vec<u8>, value: usize) {
    push_u64(
        bytes,
        u64::try_from(value).expect("conformance fixture lengths fit u64"),
    );
}

fn push_text(bytes: &mut Vec<u8>, value: &str) {
    push_len(bytes, value.len());
    bytes.extend_from_slice(value.as_bytes());
}

fn push_f64s(bytes: &mut Vec<u8>, values: &[f64]) {
    push_len(bytes, values.len());
    for value in values {
        push_u64(bytes, value.to_bits());
    }
}

fn push_nested(bytes: &mut Vec<u8>, label: &str, frame: &[u8]) {
    push_text(bytes, label);
    push_len(bytes, frame.len());
    bytes.extend_from_slice(frame);
}

fn hex_bytes(bytes: &[u8]) -> String {
    let mut encoded = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        write!(&mut encoded, "{byte:02x}").expect("writing to String cannot fail");
    }
    encoded
}

const fn sign_name(sign: Sign) -> &'static str {
    match sign {
        Sign::Negative => "negative",
        Sign::Zero => "zero",
        Sign::Positive => "positive",
    }
}

const fn stage_name(stage: Stage) -> &'static str {
    match stage {
        Stage::Filtered => "filtered",
        Stage::Adaptive => "adaptive",
        Stage::Exact => "exact",
    }
}

fn core_inputs() -> Vec<u8> {
    let mut bytes = b"fs-ivl:interval-affine-policy:v1".to_vec();
    push_text(
        &mut bytes,
        "Interval::{new,div,sub,width,WHOLE}+AffineCtx::from_interval+Affine::to_interval",
    );
    push_text(&mut bytes, "source-interval");
    push_f64s(&mut bytes, &[1.0, 2.0]);
    push_text(&mut bytes, "zero-containing-divisor");
    push_f64s(&mut bytes, &[-1.0, 1.0]);
    push_text(
        &mut bytes,
        "division-policy=zero-containing-divisor-yields-WHOLE:v1",
    );
    push_text(&mut bytes, "literal-WHOLE-endpoint-bits");
    for endpoint in WHOLE_BITS {
        push_u64(&mut bytes, endpoint);
    }
    push_text(&mut bytes, "correlation-policy=same-context-same-symbol:v1");
    push_u64(&mut bytes, AFFINE_ABSOLUTE_TOLERANCE.to_bits());
    push_u64(&mut bytes, PLAIN_MINIMUM_WIDTH.to_bits());
    bytes
}

fn taylor_root_inputs() -> Vec<u8> {
    let mut bytes = b"fs-ivl:taylor-root-bounded-policy:v1".to_vec();
    push_text(
        &mut bytes,
        "TaylorModel1::{variable,constant}+newton_roots_bounded+RootSearchConfig",
    );
    push_text(&mut bytes, "domain");
    push_f64s(&mut bytes, &[-1.0, 1.0]);
    push_text(&mut bytes, "admitted-order-zero-constant-at-zero");
    push_f64s(&mut bytes, &[0.25, 0.0]);
    push_u64(&mut bytes, 0);
    push_u64(&mut bytes, 0x3fcf_ffff_ffff_fffe);
    push_u64(&mut bytes, 0x3fd0_0000_0000_0002);
    push_text(&mut bytes, "variable-order-too-small");
    push_u64(&mut bytes, 0);
    push_text(&mut bytes, "maximum-order");
    push_u64(
        &mut bytes,
        u64::try_from(MAX_TAYLOR_ORDER).expect("Taylor order fits u64"),
    );
    push_text(&mut bytes, "maximum-order-is-admitted");
    push_text(&mut bytes, "rejected-order");
    push_u64(
        &mut bytes,
        u64::try_from(MAX_TAYLOR_ORDER + 1).expect("Taylor order fits u64"),
    );
    push_text(&mut bytes, "non-finite-constant=canonical-quiet-nan");
    push_u64(&mut bytes, 0x7ff8_0000_0000_0000);
    push_text(&mut bytes, "empty-root-budget");
    push_u64(&mut bytes, 0);
    push_text(&mut bytes, "exhaustion-root-budget");
    push_u64(&mut bytes, 1);
    push_u64(&mut bytes, f64::MIN_POSITIVE.to_bits());
    push_text(
        &mut bytes,
        "zero-function-and-derivative;incomplete=>Possible:v1",
    );
    bytes
}

fn push_points(bytes: &mut Vec<u8>, label: &str, points: &[[f64; 2]]) {
    push_text(bytes, label);
    push_len(bytes, points.len());
    for point in points {
        push_f64s(bytes, point);
    }
}

fn predicate_inputs() -> Vec<u8> {
    let mut bytes = b"fs-ivl:exact-predicate-known-answers:v1".to_vec();
    push_text(
        &mut bytes,
        "orient2d_with_stage+orient2d_sos+incircle+Sign::flip",
    );
    push_points(&mut bytes, "counterclockwise", &CCW_POINTS);
    push_points(
        &mut bytes,
        "clockwise-swap",
        &[CCW_POINTS[1], CCW_POINTS[0], CCW_POINTS[2]],
    );
    push_points(&mut bytes, "collinear", &COLLINEAR_POINTS);
    push_points(&mut bytes, "coincident-sos", &[TIE_POINT; 3]);
    push_text(&mut bytes, "sos-index-orders");
    for indices in [[0_u64, 1, 2], [1, 0, 2]] {
        for index in indices {
            push_u64(&mut bytes, index);
        }
    }
    push_points(&mut bytes, "cocircular", &CIRCLE_POINTS);
    push_text(&mut bytes, "expected-signs");
    push_len(&mut bytes, PREDICATE_REFERENCES.len());
    for sign in PREDICATE_REFERENCES {
        push_text(&mut bytes, sign_name(sign));
    }
    push_text(&mut bytes, "expected-general-position-stage");
    push_text(&mut bytes, stage_name(Stage::Filtered));
    bytes
}

fn interval_bits(interval: Interval) -> [u64; 2] {
    [interval.lo().to_bits(), interval.hi().to_bits()]
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CoreMeasurement {
    division: [u64; 2],
    affine_difference: [u64; 2],
    affine_width: u64,
    affine_contains_zero: bool,
    plain_width: u64,
}

fn measure_core() -> CoreMeasurement {
    let source = Interval::new(1.0, 2.0);
    let divisor = Interval::new(-1.0, 1.0);
    let division = source / divisor;
    let mut context = AffineCtx::new();
    let affine = context.from_interval(source);
    let affine_difference = (&affine - &affine).to_interval();
    let plain_difference = source - source;
    CoreMeasurement {
        division: interval_bits(division),
        affine_difference: interval_bits(affine_difference),
        affine_width: affine_difference.width().to_bits(),
        affine_contains_zero: affine_difference.contains(0.0),
        plain_width: plain_difference.width().to_bits(),
    }
}

fn core_outcome(input_frame: &[u8]) -> CaseOutcome {
    let run = measure_core();
    let replay = measure_core();
    let inputs_hex = hex_bytes(input_frame);
    if run != replay {
        return CaseOutcome::fail(format!(
            "stage=same-run-replay; first={run:?}; second={replay:?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-ivl/CONTRACT.md#determinism-class");
    }
    let whole = interval_bits(Interval::WHOLE);
    if whole != WHOLE_BITS {
        return CaseOutcome::fail(format!(
            "stage=whole-constant; computed={whole:016x?}; literal={WHOLE_BITS:016x?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-ivl/CONTRACT.md#public-types-and-semantics");
    }
    if run.division != WHOLE_BITS {
        return CaseOutcome::fail(format!(
            "stage=zero-divisor-policy; computed={:016x?}; reference={WHOLE_BITS:016x?}; inputs_hex={inputs_hex}",
            run.division,
        ))
        .with_evidence("crates/fs-ivl/CONTRACT.md#error-model");
    }
    let affine_width = f64::from_bits(run.affine_width);
    if !run.affine_contains_zero
        || !affine_width.is_finite()
        || affine_width > AFFINE_ABSOLUTE_TOLERANCE
    {
        return CaseOutcome::fail(format!(
            "stage=affine-correlation; interval={:016x?}; contains_zero={}; width={affine_width}; tolerance={AFFINE_ABSOLUTE_TOLERANCE}; inputs_hex={inputs_hex}",
            run.affine_difference, run.affine_contains_zero,
        ))
        .with_evidence("crates/fs-ivl/CONTRACT.md#invariants");
    }
    let plain_width = f64::from_bits(run.plain_width);
    if !plain_width.is_finite() || plain_width < PLAIN_MINIMUM_WIDTH {
        return CaseOutcome::fail(format!(
            "stage=plain-dependency-reference; width={plain_width}; minimum={PLAIN_MINIMUM_WIDTH}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-ivl/CONTRACT.md#invariants");
    }
    CaseOutcome::pass(
        "zero_divisor=WHOLE; affine_x_minus_x_contains_zero; affine_width<=1e-13; plain_width>=1; same_run_bits=identical",
    )
    .with_evidence("crates/fs-ivl/CONTRACT.md#error-model")
    .with_evidence("crates/fs-ivl/CONTRACT.md#invariants")
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ConstantReceipt {
    order: usize,
    remainder: [u64; 2],
    at_zero: [u64; 2],
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RootReceipt {
    boxes_examined: usize,
    complete: bool,
    root_count: usize,
    certified_count: usize,
    first_root: Option<[u64; 2]>,
}

fn root_receipt(report: RootSearchReport) -> RootReceipt {
    let RootSearchReport {
        roots,
        boxes_examined,
        complete,
        ..
    } = report;
    RootReceipt {
        boxes_examined,
        complete,
        root_count: roots.len(),
        certified_count: roots.iter().filter(|root| root.is_certified()).count(),
        first_root: roots.first().map(|root| interval_bits(root.interval())),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct TaylorRootMeasurement {
    constant: Result<ConstantReceipt, TaylorModelError>,
    maximum_constant: Result<usize, TaylorModelError>,
    variable_zero: Result<usize, TaylorModelError>,
    order_too_large: Result<usize, TaylorModelError>,
    non_finite_constant: Result<usize, TaylorModelError>,
    empty_budget: Result<RootReceipt, RootSearchError>,
    exhausted_budget: Result<RootReceipt, RootSearchError>,
}

fn measure_taylor_root() -> TaylorRootMeasurement {
    let domain = Interval::new(-1.0, 1.0);
    let constant = TaylorModel1::constant(0.25, domain, 0).map(|constant| ConstantReceipt {
        order: constant.order(),
        remainder: interval_bits(constant.remainder()),
        at_zero: interval_bits(constant.eval_interval(Interval::point(0.0))),
    });
    let maximum_constant =
        TaylorModel1::constant(0.0, domain, MAX_TAYLOR_ORDER).map(|model| model.order());
    let variable_zero = TaylorModel1::variable(domain, 0).map(|model| model.order());
    let order_too_large =
        TaylorModel1::constant(0.0, domain, MAX_TAYLOR_ORDER + 1).map(|model| model.order());
    let non_finite_constant =
        TaylorModel1::constant(f64::NAN, domain, 0).map(|model| model.order());
    let zero = |_x: Interval| Interval::point(0.0);
    let empty_budget = newton_roots_bounded(
        &zero,
        &zero,
        domain,
        RootSearchConfig {
            min_width: f64::MIN_POSITIVE,
            max_boxes: 0,
        },
    )
    .map(root_receipt);
    let exhausted_budget = newton_roots_bounded(
        &zero,
        &zero,
        domain,
        RootSearchConfig {
            min_width: f64::MIN_POSITIVE,
            max_boxes: 1,
        },
    )
    .map(root_receipt);
    TaylorRootMeasurement {
        constant,
        maximum_constant,
        variable_zero,
        order_too_large,
        non_finite_constant,
        empty_budget,
        exhausted_budget,
    }
}

fn taylor_root_outcome(input_frame: &[u8]) -> CaseOutcome {
    let run = measure_taylor_root();
    let replay = measure_taylor_root();
    let inputs_hex = hex_bytes(input_frame);
    if run != replay {
        return CaseOutcome::fail(format!(
            "stage=same-run-replay; first={run:?}; second={replay:?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-ivl/CONTRACT.md#determinism-class");
    }
    let expected = TaylorRootMeasurement {
        constant: Ok(ConstantReceipt {
            order: 0,
            remainder: [0, 0],
            at_zero: [0x3fcf_ffff_ffff_fffe, 0x3fd0_0000_0000_0002],
        }),
        maximum_constant: Ok(MAX_TAYLOR_ORDER),
        variable_zero: Err(TaylorModelError::VariableOrderTooSmall {
            requested: 0,
            minimum: 1,
        }),
        order_too_large: Err(TaylorModelError::OrderTooLarge {
            requested: MAX_TAYLOR_ORDER + 1,
            maximum: MAX_TAYLOR_ORDER,
        }),
        non_finite_constant: Err(TaylorModelError::NonFiniteConstant),
        empty_budget: Err(RootSearchError::EmptyBudget),
        exhausted_budget: Ok(RootReceipt {
            boxes_examined: 1,
            complete: false,
            root_count: 1,
            certified_count: 0,
            first_root: Some(interval_bits(Interval::new(-1.0, 1.0))),
        }),
    };
    if run != expected {
        return CaseOutcome::fail(format!(
            "stage=taylor-root-policy; computed={run:?}; reference={expected:?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-ivl/CONTRACT.md#error-model")
        .with_evidence("crates/fs-ivl/CONTRACT.md#cancellation-behavior");
    }
    CaseOutcome::pass(
        "taylor_constant_order0=admitted; taylor_maximum_order=admitted; taylor_refusals=variable-order,above-maximum-order,non-finite; root_empty_budget=refused; one_box_examined=1; complete=false; roots=1-Possible; same_run=identical",
    )
    .with_evidence("crates/fs-ivl/CONTRACT.md#error-model")
    .with_evidence("crates/fs-ivl/CONTRACT.md#cancellation-behavior")
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PredicateMeasurement {
    signs: [Sign; 6],
    general_position_stages: [Stage; 2],
}

fn measure_predicates() -> PredicateMeasurement {
    let (ccw, ccw_stage) = orient2d_with_stage(CCW_POINTS[0], CCW_POINTS[1], CCW_POINTS[2]);
    let (cw, cw_stage) = orient2d_with_stage(CCW_POINTS[1], CCW_POINTS[0], CCW_POINTS[2]);
    let collinear = orient2d_with_stage(
        COLLINEAR_POINTS[0],
        COLLINEAR_POINTS[1],
        COLLINEAR_POINTS[2],
    )
    .0;
    let sos = orient2d_sos(TIE_POINT, TIE_POINT, TIE_POINT, 0, 1, 2);
    let sos_swapped = orient2d_sos(TIE_POINT, TIE_POINT, TIE_POINT, 1, 0, 2);
    let circle = incircle(
        CIRCLE_POINTS[0],
        CIRCLE_POINTS[1],
        CIRCLE_POINTS[2],
        CIRCLE_POINTS[3],
    );
    PredicateMeasurement {
        signs: [ccw, cw, collinear, sos, sos_swapped, circle],
        general_position_stages: [ccw_stage, cw_stage],
    }
}

#[derive(Debug, Clone, Copy)]
struct Corruption {
    seed: u64,
    vector: usize,
}

fn predicate_outcome(
    reference: [Sign; 6],
    corruption: Option<Corruption>,
    input_frame: &[u8],
) -> CaseOutcome {
    let run = measure_predicates();
    let replay = measure_predicates();
    let inputs_hex = hex_bytes(input_frame);
    let context = corruption.map_or_else(
        || "mode=canonical".to_owned(),
        |corruption| {
            format!(
                "seed=0x{:016x}; vector={}",
                corruption.seed, corruption.vector
            )
        },
    );
    if run != replay {
        return CaseOutcome::fail(format!(
            "{context}; stage=same-run-replay; first={run:?}; second={replay:?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-ivl/CONTRACT.md#determinism-class");
    }
    if run.general_position_stages != [Stage::Filtered, Stage::Filtered] {
        return CaseOutcome::fail(format!(
            "{context}; stage=general-position-filter; computed={:?}; reference=[Filtered, Filtered]; inputs_hex={inputs_hex}",
            run.general_position_stages,
        ))
        .with_evidence("crates/fs-ivl/CONTRACT.md#public-types-and-semantics");
    }
    if run.signs != reference {
        return CaseOutcome::fail(format!(
            "{context}; stage=predicate-known-answers; computed={:?}; reference={reference:?}; inputs_hex={inputs_hex}",
            run.signs,
        ))
        .with_evidence("crates/fs-ivl/CONTRACT.md#public-types-and-semantics");
    }
    if run.signs[1] != run.signs[0].flip() || run.signs[4] != run.signs[3].flip() {
        return CaseOutcome::fail(format!(
            "{context}; stage=predicate-antisymmetry; signs={:?}; inputs_hex={inputs_hex}",
            run.signs,
        ))
        .with_evidence("crates/fs-ivl/CONTRACT.md#public-types-and-semantics");
    }
    CaseOutcome::pass(
        "orient2d=Positive/Negative/Zero; general_position=Filtered; sos_tie=Positive/Negative; incircle=Zero; antisymmetry=exact; same_run=identical",
    )
    .with_evidence("crates/fs-ivl/CONTRACT.md#public-types-and-semantics")
}

#[test]
fn structured_casebook_emits_replay_complete_green_records() {
    assert_eq!(CASEBOOK_RECORD_VERSION, 1);
    let core_frame = core_inputs();
    let taylor_frame = taylor_root_inputs();
    let predicate_frame = predicate_inputs();
    let core_digest = fnv1a64(&core_frame);
    let taylor_digest = fnv1a64(&taylor_frame);
    let predicate_digest = fnv1a64(&predicate_frame);
    assert_eq!(core_digest, 0x6993_9f28_0ebc_0702);
    assert_eq!(taylor_digest, 0xe6f2_d8ff_51ad_3676);
    assert_eq!(predicate_digest, 0x0162_cccc_05db_6ba0);

    let report = Suite::new(SUITE)
        .case(
            "interval-affine-fail-closed-policy",
            core_digest,
            ToleranceSpec::AbsoluteLe(AFFINE_ABSOLUTE_TOLERANCE),
            move || core_outcome(&core_frame),
        )
        .case(
            "taylor-root-bounded-policy",
            taylor_digest,
            ToleranceSpec::Exact,
            move || taylor_root_outcome(&taylor_frame),
        )
        .case(
            "exact-predicate-known-answers",
            predicate_digest,
            ToleranceSpec::Exact,
            move || predicate_outcome(PREDICATE_REFERENCES, None, &predicate_frame),
        )
        .run();

    report.assert_green();
    assert_eq!(
        report
            .records
            .iter()
            .map(|record| record.case.as_str())
            .collect::<Vec<_>>(),
        [
            "interval-affine-fail-closed-policy",
            "taylor-root-bounded-policy",
            "exact-predicate-known-answers",
        ]
    );
    assert_eq!(
        report.records[0].json_line(),
        format!(
            concat!(
                "{{\"casebook\":{},\"suite\":\"fs-ivl/structured-conformance-v1\",",
                "\"case\":\"interval-affine-fail-closed-policy\",\"inputs_digest\":\"69939f280ebc0702\",",
                "\"tolerance\":\"abs<=1e-13\",\"pass\":true,",
                "\"details\":\"zero_divisor=WHOLE; affine_x_minus_x_contains_zero; affine_width<=1e-13; plain_width>=1; same_run_bits=identical\",",
                "\"evidence\":[\"crates/fs-ivl/CONTRACT.md#error-model\",",
                "\"crates/fs-ivl/CONTRACT.md#invariants\"]}}"
            ),
            CASEBOOK_RECORD_VERSION,
        ),
        "the structured record schema and field order are contract"
    );
}

#[test]
#[allow(clippy::too_many_lines)] // Keep the corruption frame and every refusal assertion together.
fn disclosed_seeded_corruption_is_replay_identical_and_turns_the_suite_red() {
    const CORRUPTION_SEED: u64 = 0xF51E_0000;
    let vector = (CORRUPTION_SEED & 0x7) as usize;
    assert_eq!(vector, 0);
    let canonical = PREDICATE_REFERENCES;
    let mut corrupted = canonical;
    corrupted[vector] = corrupted[vector].flip();
    assert_eq!(corrupted[vector], Sign::Negative);
    let corruption = Corruption {
        seed: CORRUPTION_SEED,
        vector,
    };

    let predicate_frame = predicate_inputs();
    let mut inputs = b"fs-ivl:seeded-predicate-reference-corruption:v1".to_vec();
    push_u64(&mut inputs, CORRUPTION_SEED);
    push_len(&mut inputs, vector);
    push_nested(&mut inputs, "nested-predicate-frame", &predicate_frame);
    push_text(&mut inputs, "canonical-signs");
    for sign in canonical {
        push_text(&mut inputs, sign_name(sign));
    }
    push_text(&mut inputs, "corrupted-signs");
    for sign in corrupted {
        push_text(&mut inputs, sign_name(sign));
    }
    let inputs_digest = fnv1a64(&inputs);
    assert_eq!(inputs_digest, 0x3854_abe5_7cb6_7750);

    let make_report = || {
        let input_frame = inputs.clone();
        Suite::new(SUITE)
            .case(
                "seeded-predicate-reference-corruption",
                inputs_digest,
                ToleranceSpec::Exact,
                move || predicate_outcome(corrupted, Some(corruption), &input_frame),
            )
            .run()
    };
    let first = make_report();
    let replay = make_report();

    assert!(!first.all_passed());
    assert!(!replay.all_passed());
    let first_failures = first.failures();
    let replay_failures = replay.failures();
    let [first_failure] = first_failures.as_slice() else {
        panic!("the seeded corruption must produce exactly one structured failure");
    };
    let [replay_failure] = replay_failures.as_slice() else {
        panic!("the replayed corruption must produce exactly one structured failure");
    };
    assert_eq!(first_failure.case, "seeded-predicate-reference-corruption");
    assert_eq!(first_failure.inputs_digest, "3854abe57cb67750");
    assert_eq!(first_failure.json_line(), replay_failure.json_line());
    assert!(
        first_failure
            .details
            .contains(&format!("seed=0x{CORRUPTION_SEED:016x}; vector={vector}"))
    );
    assert!(
        first_failure
            .details
            .contains("stage=predicate-known-answers")
    );
    assert!(first_failure.details.contains("computed=[Positive"));
    assert!(first_failure.details.contains("reference=[Negative"));
    assert!(first_failure.details.contains("inputs_hex="));
    assert!(
        first_failure
            .json_line()
            .contains("\"tolerance\":\"exact\",\"pass\":false")
    );

    let panic = std::panic::catch_unwind(|| first.assert_green())
        .expect_err("the merge-gate assertion must reject the seeded failure");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("casebook panic carries text");
    assert!(message.contains("seeded-predicate-reference-corruption"));
    assert!(message.contains(&format!("seed=0x{CORRUPTION_SEED:016x}")));
}
