//! Structured BEDROCK conformance records for fs-ad.
//!
//! The crate's source-level batteries retain broad primal-fidelity, gradient,
//! IFT, checkpoint, spill, and optional tape-bridge coverage. These bounded
//! cases expose exact forward-AD and treeverse policies through fs-casebook so
//! a central run can diagnose and replay failures from structured records
//! (Gauntlet G0).

use fs_ad::{
    Real, checkpointed_adjoint, full_adjoint, gradient, jvp, min_budget, second_directional,
};
use fs_casebook::{CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, ToleranceSpec, fnv1a64};

const SUITE: &str = "fs-ad/bedrock-conformance-v1";
const POLYNOMIAL_FORMULA: &str = "(((x*x)+((3*x)*y))+((2*y)*y)):real-v1";
const POLYNOMIAL_POINT: [f64; 2] = [1.5, -0.5];
const POLYNOMIAL_DIRECTION: [f64; 2] = [1.0, 2.0];
const POLYNOMIAL_OUTPUTS: [(&str, f64); 8] = [
    ("gradient.value", 0.5),
    ("gradient.dx", 1.5),
    ("gradient.dy", 2.5),
    ("jvp.value", 0.5),
    ("jvp.directional", 6.5),
    ("second_directional.value", 0.5),
    ("second_directional.first", 6.5),
    ("second_directional.second", 30.0),
];
const REVOLVE_STEPS: usize = 100;
const REVOLVE_BUDGET: usize = 8;
const REVOLVE_LOG2_BOUND: u64 = 7;
const REVOLVE_MIN_FORWARD: u64 = 1;
const REVOLVE_MAX_FORWARD: u64 = 700;
const REVOLVE_MIN_PEAK: usize = 1;
const REFUSAL_STEPS: usize = 64;
const REFUSAL_BUDGET: usize = 2;
const REVOLVE_X0: f64 = 0.3;
const REVOLVE_THETA: f64 = 0.8;

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

fn polynomial<T: Real>([x, y]: [T; 2]) -> T {
    x * x + T::from_f64(3.0) * x * y + T::from_f64(2.0) * y * y
}

fn analytic_inputs() -> Vec<u8> {
    let mut bytes = b"fs-ad:dual-polynomial-known-answers:v1".to_vec();
    push_text(&mut bytes, "gradient+jvp+second_directional");
    push_text(&mut bytes, POLYNOMIAL_FORMULA);
    push_text(&mut bytes, "point[x,y]");
    push_f64s(&mut bytes, &POLYNOMIAL_POINT);
    push_text(&mut bytes, "direction[x,y]");
    push_f64s(&mut bytes, &POLYNOMIAL_DIRECTION);
    for (output, reference) in POLYNOMIAL_OUTPUTS {
        push_text(&mut bytes, output);
        push_u64(&mut bytes, reference.to_bits());
    }
    bytes
}

fn boundary_inputs() -> Vec<u8> {
    let mut bytes = b"fs-ad:dual-boundary-policy:v1".to_vec();
    for (operation, point, expected_value, expected_derivative) in [
        ("abs", 0.0_f64, 0.0_f64, 0.0_f64),
        ("sqrt", 0.0, 0.0, f64::INFINITY),
        ("asin", 1.0, core::f64::consts::FRAC_PI_2, f64::INFINITY),
        ("acos", -1.0, core::f64::consts::PI, f64::NEG_INFINITY),
    ] {
        push_text(&mut bytes, operation);
        push_text(&mut bytes, "point-f64");
        push_u64(&mut bytes, point.to_bits());
        push_text(&mut bytes, "expected-value");
        push_u64(&mut bytes, expected_value.to_bits());
        push_text(&mut bytes, "expected-derivative");
        push_u64(&mut bytes, expected_derivative.to_bits());
    }
    push_text(&mut bytes, "powi");
    push_text(&mut bytes, "point-f64");
    push_u64(&mut bytes, 1.0_f64.to_bits());
    push_text(&mut bytes, "exponent-i32-le");
    bytes.extend_from_slice(&i32::MIN.to_le_bytes());
    push_text(&mut bytes, "expected-value");
    push_u64(&mut bytes, 1.0_f64.to_bits());
    push_text(&mut bytes, "expected-derivative");
    push_u64(&mut bytes, f64::from(i32::MIN).to_bits());
    bytes
}

fn revolve_inputs() -> Vec<u8> {
    let mut bytes = b"fs-ad:treeverse-bitwise-budget-policy:v1".to_vec();
    for operation in ["full_adjoint", "checkpointed_adjoint", "min_budget"] {
        push_text(&mut bytes, operation);
    }
    push_text(&mut bytes, "forward:(-0.1).mul_add(((x*x)*x)-theta,x):v1");
    push_text(
        &mut bytes,
        "reverse:((((-0.3*x)*x)+1.0)*xbar,0.1.mul_add(xbar,thetabar)):v1",
    );
    for (field, value) in [
        ("steps", REVOLVE_STEPS),
        ("snapshot-budget", REVOLVE_BUDGET),
        ("refusal-steps", REFUSAL_STEPS),
        ("refusal-budget", REFUSAL_BUDGET),
    ] {
        push_text(&mut bytes, field);
        push_len(&mut bytes, value);
    }
    push_text(&mut bytes, "x0+theta+seed[xbar,thetabar]");
    push_f64s(&mut bytes, &[REVOLVE_X0, REVOLVE_THETA, 1.0, 0.0]);
    push_text(&mut bytes, "ceil-log2-steps");
    push_u64(&mut bytes, REVOLVE_LOG2_BOUND);
    push_text(&mut bytes, "min-forward-steps");
    push_u64(&mut bytes, REVOLVE_MIN_FORWARD);
    push_text(&mut bytes, "max-forward-steps");
    push_u64(&mut bytes, REVOLVE_MAX_FORWARD);
    push_text(&mut bytes, "min-peak-snapshots");
    push_len(&mut bytes, REVOLVE_MIN_PEAK);
    push_text(&mut bytes, "max-peak-snapshots");
    push_len(&mut bytes, REVOLVE_BUDGET);
    push_text(&mut bytes, "required-panic-fragment");
    push_text(&mut bytes, "snapshot budget 2 < required");
    bytes
}

fn polynomial_outcome() -> CaseOutcome {
    let point = POLYNOMIAL_POINT;
    let direction = POLYNOMIAL_DIRECTION;
    let (value, grad) = gradient(point, polynomial);
    let (jvp_value, directional) = jvp(point, direction, polynomial);
    let (second_value, first_directional, second) =
        second_directional(point, direction, polynomial);
    let computed = [
        value,
        grad[0],
        grad[1],
        jvp_value,
        directional,
        second_value,
        first_directional,
        second,
    ];
    for (index, (computed, (output, reference))) in
        computed.into_iter().zip(POLYNOMIAL_OUTPUTS).enumerate()
    {
        if computed.to_bits() != reference.to_bits() {
            return CaseOutcome::fail(format!(
                "operation=gradient+jvp+second_directional; formula={POLYNOMIAL_FORMULA}; point=[1.5,-0.5]; direction=[1,2]; output={output}; output_index={index}; computed_bits=0x{:016x}; reference_bits=0x{:016x}; computed={computed}; reference={reference}",
                computed.to_bits(),
                reference.to_bits(),
            ))
            .with_evidence("crates/fs-ad/CONTRACT.md#public-types-and-semantics")
            .with_evidence("crates/fs-ad/CONTRACT.md#invariants");
        }
    }
    CaseOutcome::pass(format!(
        "operation=gradient+jvp+second_directional; formula={POLYNOMIAL_FORMULA}; point=[1.5,-0.5]; direction=[1,2]; value=0.5; gradient=[1.5,2.5]; jvp=6.5; second_directional=30"
    ))
    .with_evidence("crates/fs-ad/CONTRACT.md#public-types-and-semantics")
    .with_evidence("crates/fs-ad/CONTRACT.md#invariants")
}

fn boundary_outcome() -> CaseOutcome {
    let (abs_value, abs_grad) = gradient([0.0], |[x]| x.abs());
    let (sqrt_value, sqrt_grad) = gradient([0.0], |[x]| x.sqrt());
    let (asin_value, asin_grad) = gradient([1.0], |[x]| x.asin());
    let (acos_value, acos_grad) = gradient([-1.0], |[x]| x.acos());
    let (powi_value, powi_grad) = gradient([1.0], |[x]| x.powi(i32::MIN));
    let results = [
        ("abs", "point=+0", "value", abs_value, 0.0),
        ("abs", "point=+0", "derivative", abs_grad[0], 0.0),
        ("sqrt", "point=+0", "value", sqrt_value, 0.0),
        (
            "sqrt",
            "point=+0",
            "derivative",
            sqrt_grad[0],
            f64::INFINITY,
        ),
        (
            "asin",
            "point=1",
            "value",
            asin_value,
            core::f64::consts::FRAC_PI_2,
        ),
        ("asin", "point=1", "derivative", asin_grad[0], f64::INFINITY),
        (
            "acos",
            "point=-1",
            "value",
            acos_value,
            core::f64::consts::PI,
        ),
        (
            "acos",
            "point=-1",
            "derivative",
            acos_grad[0],
            f64::NEG_INFINITY,
        ),
        (
            "powi",
            "point=1; exponent=-2147483648i32",
            "value",
            powi_value,
            1.0,
        ),
        (
            "powi",
            "point=1; exponent=-2147483648i32",
            "derivative",
            powi_grad[0],
            f64::from(i32::MIN),
        ),
    ];
    for (index, (operation, input, component, computed, reference)) in
        results.into_iter().enumerate()
    {
        if computed.to_bits() != reference.to_bits() {
            return CaseOutcome::fail(format!(
                "operation={operation}; {input}; component={component}; output_index={index}; computed_bits=0x{:016x}; reference_bits=0x{:016x}; computed={computed}; reference={reference}",
                computed.to_bits(),
                reference.to_bits(),
            ))
            .with_evidence("crates/fs-ad/CONTRACT.md#public-types-and-semantics")
            .with_evidence("crates/fs-ad/CONTRACT.md#invariants");
        }
    }
    CaseOutcome::pass(
        "abs_prime_zero=+0; sqrt_prime_zero=+inf; asin_prime_one=+inf; acos_prime_minus_one=-inf; powi_i32_min_prime=-2147483648",
    )
    .with_evidence("crates/fs-ad/CONTRACT.md#public-types-and-semantics")
    .with_evidence("crates/fs-ad/CONTRACT.md#invariants")
}

fn forward(theta: f64) -> impl Fn(usize, &f64) -> f64 {
    move |_index, &x| (-0.1_f64).mul_add(x * x * x - theta, x)
}

fn reverse(_theta: f64) -> impl Fn(usize, &f64, (f64, f64)) -> (f64, f64) {
    move |_index, &x, (xbar, thetabar)| {
        (
            (-0.3_f64 * x * x + 1.0) * xbar,
            0.1_f64.mul_add(xbar, thetabar),
        )
    }
}

fn revolve_outcome() -> CaseOutcome {
    let forward = forward(REVOLVE_THETA);
    let reverse = reverse(REVOLVE_THETA);
    let budget = min_budget(REVOLVE_STEPS);
    if budget != REVOLVE_BUDGET {
        return CaseOutcome::fail(format!(
            "steps={REVOLVE_STEPS}; computed_min_budget={budget}; expected_min_budget={REVOLVE_BUDGET}"
        ))
        .with_evidence("crates/fs-ad/CONTRACT.md#public-types-and-semantics");
    }
    let full = full_adjoint(&REVOLVE_X0, REVOLVE_STEPS, &forward, &reverse, (1.0, 0.0));
    let (checkpointed, stats) = checkpointed_adjoint(
        &REVOLVE_X0,
        REVOLVE_STEPS,
        budget,
        &forward,
        &reverse,
        (1.0, 0.0),
    );
    let full_bits = [full.0.to_bits(), full.1.to_bits()];
    let checkpointed_bits = [checkpointed.0.to_bits(), checkpointed.1.to_bits()];
    if checkpointed_bits != full_bits {
        return CaseOutcome::fail(format!(
            "steps={REVOLVE_STEPS}; budget={budget}; full_bits={full_bits:016x?}; checkpointed_bits={checkpointed_bits:016x?}; forward_steps={}; peak_snapshots={}",
            stats.forward_steps,
            stats.peak_snapshots,
        ))
        .with_evidence("crates/fs-ad/CONTRACT.md#public-types-and-semantics");
    }

    let max_forward =
        u64::try_from(REVOLVE_STEPS).expect("step count fits u64") * REVOLVE_LOG2_BOUND;
    if max_forward != REVOLVE_MAX_FORWARD {
        return CaseOutcome::fail(format!(
            "steps={REVOLVE_STEPS}; log2_bound={REVOLVE_LOG2_BOUND}; computed_max_forward={max_forward}; expected_max_forward={REVOLVE_MAX_FORWARD}"
        ))
        .with_evidence("crates/fs-ad/CONTRACT.md#public-types-and-semantics");
    }
    if stats.forward_steps < REVOLVE_MIN_FORWARD
        || stats.forward_steps > REVOLVE_MAX_FORWARD
        || stats.peak_snapshots < REVOLVE_MIN_PEAK
        || stats.peak_snapshots > REVOLVE_BUDGET
    {
        return CaseOutcome::fail(format!(
            "steps={REVOLVE_STEPS}; budget={budget}; log2_bound={REVOLVE_LOG2_BOUND}; min_forward={REVOLVE_MIN_FORWARD}; forward_steps={}; max_forward={REVOLVE_MAX_FORWARD}; min_peak={REVOLVE_MIN_PEAK}; peak_snapshots={}; max_peak={REVOLVE_BUDGET}; full_bits={full_bits:016x?}",
            stats.forward_steps,
            stats.peak_snapshots,
        ))
        .with_evidence("crates/fs-ad/CONTRACT.md#public-types-and-semantics");
    }

    let refusal_forward = forward(REVOLVE_THETA);
    let refusal_reverse = reverse(REVOLVE_THETA);
    let refused = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
        checkpointed_adjoint(
            &REVOLVE_X0,
            REFUSAL_STEPS,
            REFUSAL_BUDGET,
            &refusal_forward,
            &refusal_reverse,
            (1.0, 0.0),
        )
    }));
    let panic = match refused {
        Err(panic) => panic,
        Ok((bar, refusal_stats)) => {
            return CaseOutcome::fail(format!(
                "refusal_steps=64; refusal_budget=2; unexpected_bar_bits=[0x{:016x},0x{:016x}]; unexpected_forward_steps={}; unexpected_peak_snapshots={}",
                bar.0.to_bits(),
                bar.1.to_bits(),
                refusal_stats.forward_steps,
                refusal_stats.peak_snapshots,
            ))
            .with_evidence("crates/fs-ad/CONTRACT.md#public-types-and-semantics");
        }
    };
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .unwrap_or("non-string panic payload");
    if !message.contains("snapshot budget 2 < required") || !message.contains("64 steps") {
        return CaseOutcome::fail(format!(
            "refusal_steps=64; refusal_budget=2; panic={message:?}; expected_fragments=[\"snapshot budget 2 < required\",\"64 steps\"]"
        ))
        .with_evidence("crates/fs-ad/CONTRACT.md#public-types-and-semantics");
    }

    CaseOutcome::pass(format!(
        "steps={REVOLVE_STEPS}; budget={budget}; log2_bound={REVOLVE_LOG2_BOUND}; forward_steps={}; peak_snapshots={}; full_bits={full_bits:016x?}; insufficient_budget_refused=true",
        stats.forward_steps,
        stats.peak_snapshots,
    ))
    .with_evidence("crates/fs-ad/CONTRACT.md#public-types-and-semantics")
}

#[test]
fn bedrock_casebook_suite_emits_replay_complete_green_records() {
    let analytic_digest = fnv1a64(&analytic_inputs());
    let boundary_digest = fnv1a64(&boundary_inputs());
    let revolve_digest = fnv1a64(&revolve_inputs());
    assert_eq!(analytic_digest, 0xb978_99e3_2eea_35e6);
    assert_eq!(boundary_digest, 0x62a4_f4f0_1c0b_e3c9);
    assert_eq!(revolve_digest, 0x4c0b_c870_5b1b_a371);

    let report = Suite::new(SUITE)
        .case(
            "dual-polynomial-known-answers",
            analytic_digest,
            ToleranceSpec::Exact,
            polynomial_outcome,
        )
        .case(
            "dual-boundary-policy",
            boundary_digest,
            ToleranceSpec::Exact,
            boundary_outcome,
        )
        .case(
            "treeverse-bitwise-budget-policy",
            revolve_digest,
            ToleranceSpec::Structural,
            revolve_outcome,
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
            "dual-polynomial-known-answers",
            "dual-boundary-policy",
            "treeverse-bitwise-budget-policy",
        ]
    );
    assert_eq!(
        report.records[0].json_line(),
        format!(
            concat!(
                "{{\"casebook\":{},\"suite\":\"fs-ad/bedrock-conformance-v1\",",
                "\"case\":\"dual-polynomial-known-answers\",\"inputs_digest\":\"b97899e32eea35e6\",",
                "\"tolerance\":\"exact\",\"pass\":true,",
                "\"details\":\"operation=gradient+jvp+second_directional; formula=(((x*x)+((3*x)*y))+((2*y)*y)):real-v1; point=[1.5,-0.5]; direction=[1,2]; value=0.5; gradient=[1.5,2.5]; jvp=6.5; second_directional=30\",",
                "\"evidence\":[\"crates/fs-ad/CONTRACT.md#public-types-and-semantics\",",
                "\"crates/fs-ad/CONTRACT.md#invariants\"]}}"
            ),
            CASEBOOK_RECORD_VERSION,
        ),
        "the structured analytic record schema and field order are contract"
    );
}

#[test]
fn disclosed_seeded_corruption_turns_the_casebook_suite_red() {
    const CORRUPTION_SEED: u64 = 0xF5AD_0001;
    let component = (CORRUPTION_SEED as usize) % 2;
    let bit = ((CORRUPTION_SEED >> 1) % 64) as u32;
    assert_eq!(component, 1);
    assert_eq!(bit, 0);
    let canonical = [1.5_f64.to_bits(), 2.5_f64.to_bits()];
    let mut corrupted = canonical;
    corrupted[component] ^= 1_u64 << bit;

    let analytic_inputs = analytic_inputs();
    let mut inputs = b"fs-ad:seeded-polynomial-gradient-corruption:v1".to_vec();
    push_u64(&mut inputs, CORRUPTION_SEED);
    push_len(&mut inputs, component);
    push_u64(&mut inputs, u64::from(bit));
    push_len(&mut inputs, analytic_inputs.len());
    inputs.extend_from_slice(&analytic_inputs);
    for value in canonical.into_iter().chain(corrupted) {
        push_u64(&mut inputs, value);
    }
    let inputs_digest = fnv1a64(&inputs);
    assert_eq!(inputs_digest, 0xbdcc_e701_1112_1b5b);

    let report = Suite::new(SUITE)
        .case(
            "seeded-polynomial-gradient-corruption",
            inputs_digest,
            ToleranceSpec::Exact,
            move || {
                let (_, computed) = gradient([1.5, -0.5], polynomial);
                let computed = computed.map(f64::to_bits);
                if computed == corrupted {
                    CaseOutcome::pass("seeded corruption was not detected")
                } else {
                    CaseOutcome::fail(format!(
                        "seed=0x{CORRUPTION_SEED:016x}; operation=gradient; formula={POLYNOMIAL_FORMULA}; point=[1.5,-0.5]; component={component}; component_name=gradient.dy; bit={bit}; computed={computed:016x?}; canonical={canonical:016x?}; corrupted={corrupted:016x?}"
                    ))
                    .with_evidence("crates/fs-ad/tests/conformance.rs#seeded-corruption")
                }
            },
        )
        .run();

    assert!(
        !report.all_passed(),
        "the deliberately corrupted oracle must turn red"
    );
    let failures = report.failures();
    let [failure] = failures.as_slice() else {
        panic!("the seeded corruption must produce exactly one structured failure");
    };
    assert_eq!(failure.case, "seeded-polynomial-gradient-corruption");
    assert_eq!(failure.inputs_digest, "bdcce70111121b5b");
    assert!(
        failure
            .details
            .contains(&format!("seed=0x{CORRUPTION_SEED:016x}"))
    );
    assert!(failure.details.contains(&format!("component={component}")));
    assert!(
        failure
            .details
            .contains(&format!("component_name=gradient.dy; bit={bit}"))
    );
    assert!(failure.details.contains("computed=["));
    assert!(failure.details.contains("canonical=["));
    assert!(failure.details.contains("corrupted=["));
    let line = failure.json_line();
    assert!(line.contains("\"tolerance\":\"exact\",\"pass\":false"));

    let panic = std::panic::catch_unwind(|| report.assert_green())
        .expect_err("the merge-gate assertion must reject the seeded failure");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("casebook panic carries text");
    assert!(message.contains("seeded-polynomial-gradient-corruption"));
    assert!(message.contains(&format!("seed=0x{CORRUPTION_SEED:016x}")));
}
