//! Structured BEDROCK conformance records for fs-cheb.
//!
//! The existing batteries retain the adaptive, root, budget, variant,
//! Orr-Sommerfeld, and admitted aggregate-golden coverage. These bounded G0
//! cases are the cheap PR subset: direct-coefficient spectral semantics,
//! checkpoint refusal, and the smallest exact collocation matrix, emitted as
//! replay-complete fs-casebook records.

use fs_casebook::{CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, ToleranceSpec, fnv1a64};
use fs_cheb::{Cheb1, diff_matrix, lobatto_points};

const SUITE: &str = "fs-cheb/bedrock-conformance-v1";
const CLENSHAW_CONVENTION: &str =
    "reverse(c1..):b0=(2*t).mul_add(b1,c-b2);value=t.mul_add(b1,0.5.mul_add(c0,-b2)):v1";
const DERIVATIVE_CONVENTION: &str =
    "reverse(k=1..n):d[k-1]=(2*k).mul_add(c[k],d[k+1]);scale=2/(b-a):v1";
const INTEGRAL_CONVENTION: &str = "half_width*(c0+sum(even k>=2,2*c[k]/(1-k*k))):stable-sum-v1";
const DOMAIN: [f64; 2] = [-1.0, 1.0];
const COEFFICIENTS: [f64; 4] = [2.0, 3.0, 0.0, 5.0];
const EVAL_POINTS: [f64; 4] = [-1.0, 0.0, 0.5, 1.0];
const KNOWN_OUTPUTS: [(&str, f64); 9] = [
    ("eval[-1]", -7.0),
    ("eval[0]", 1.0),
    ("eval[0.5]", -2.5),
    ("eval[1]", 9.0),
    ("differentiate.coeff[0]", 36.0),
    ("differentiate.coeff[1]", 0.0),
    ("differentiate.coeff[2]", 30.0),
    ("differentiate.eval[0.5]", 3.0),
    ("integral[-1,1]", 2.0),
];
const CHECKPOINT_POINT: f64 = 0.25;
const CHECKPOINT_VALUE: f64 = -1.6875;
const EXPECTED_POLLS: usize = 3;
const REFUSAL_POLL: usize = 2;

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

fn push_u64s(bytes: &mut Vec<u8>, values: &[u64]) {
    push_len(bytes, values.len());
    for value in values {
        push_u64(bytes, *value);
    }
}

fn fixture() -> Cheb1 {
    Cheb1::from_coeffs(DOMAIN[0], DOMAIN[1], COEFFICIENTS.to_vec())
}

fn known_answer_inputs() -> Vec<u8> {
    let mut bytes = b"fs-cheb:clenshaw-calculus-known-answers:v1".to_vec();
    push_text(&mut bytes, "Cheb1::from_coeffs+eval+differentiate+integral");
    push_text(&mut bytes, CLENSHAW_CONVENTION);
    push_text(&mut bytes, DERIVATIVE_CONVENTION);
    push_text(&mut bytes, INTEGRAL_CONVENTION);
    push_text(&mut bytes, "domain[a,b]");
    push_f64s(&mut bytes, &DOMAIN);
    push_text(&mut bytes, "stored-coefficients[c0,c1,c2,c3]");
    push_f64s(&mut bytes, &COEFFICIENTS);
    push_text(&mut bytes, "expected-degree");
    push_len(&mut bytes, 3);
    push_text(&mut bytes, "evaluation-points");
    push_f64s(&mut bytes, &EVAL_POINTS);
    for (output, reference) in KNOWN_OUTPUTS {
        push_text(&mut bytes, output);
        push_u64(&mut bytes, reference.to_bits());
    }
    bytes
}

fn checkpoint_inputs() -> Vec<u8> {
    let known = known_answer_inputs();
    let mut bytes = b"fs-cheb:eval-checkpoint-refusal-policy:v1".to_vec();
    push_text(&mut bytes, "Cheb1::eval_with_checkpoint");
    push_text(&mut bytes, "nested-known-answer-frame");
    push_len(&mut bytes, known.len());
    bytes.extend_from_slice(&known);
    push_text(&mut bytes, "point");
    push_u64(&mut bytes, CHECKPOINT_POINT.to_bits());
    push_text(&mut bytes, "expected-value");
    push_u64(&mut bytes, CHECKPOINT_VALUE.to_bits());
    push_text(&mut bytes, "expected-success-polls");
    push_len(&mut bytes, EXPECTED_POLLS);
    push_text(&mut bytes, "refusal-poll");
    push_len(&mut bytes, REFUSAL_POLL);
    push_text(&mut bytes, "CheckpointStop::RefusedAt(2)");
    bytes
}

fn lobatto_inputs() -> Vec<u8> {
    let mut bytes = b"fs-cheb:lobatto-negative-sum-kat:v1".to_vec();
    push_text(&mut bytes, "lobatto_points+diff_matrix");
    push_text(&mut bytes, "n=1; descending-points; row-major-matrix");
    push_text(
        &mut bytes,
        "offdiag=(c_i/c_j)/(x_i-x_j);diagonal=-sum(offdiag in ascending j):v1",
    );
    push_len(&mut bytes, 1);
    push_text(&mut bytes, "expected-points");
    push_f64s(&mut bytes, &[1.0, -1.0]);
    push_text(&mut bytes, "expected-row-major-diff-matrix");
    push_f64s(&mut bytes, &[0.5, -0.5, 0.5, -0.5]);
    push_text(&mut bytes, "expected-negative-sum-row-zeros");
    push_u64s(&mut bytes, &[0.0_f64.to_bits(), 0.0_f64.to_bits()]);
    bytes
}

fn known_answer_outcome() -> CaseOutcome {
    let cheb = fixture();
    if cheb.degree() != 3 {
        return CaseOutcome::fail(format!(
            "operation=Cheb1::degree; clenshaw={CLENSHAW_CONVENTION}; differentiate={DERIVATIVE_CONVENTION}; integral={INTEGRAL_CONVENTION}; domain=[-1,1]; stored_coeffs=[2,3,0,5]; computed={}; reference=3",
            cheb.degree(),
        ))
        .with_evidence("crates/fs-cheb/CONTRACT.md#public-types-and-semantics");
    }

    let derivative = cheb.differentiate();
    if derivative.coeffs().len() != 3 {
        return CaseOutcome::fail(format!(
            "operation=Cheb1::differentiate; clenshaw={CLENSHAW_CONVENTION}; differentiate={DERIVATIVE_CONVENTION}; integral={INTEGRAL_CONVENTION}; domain=[-1,1]; stored_coeffs=[2,3,0,5]; computed_coeff_count={}; reference_coeff_count=3",
            derivative.coeffs().len(),
        ))
        .with_evidence("crates/fs-cheb/CONTRACT.md#public-types-and-semantics");
    }

    let computed = [
        cheb.eval(EVAL_POINTS[0]),
        cheb.eval(EVAL_POINTS[1]),
        cheb.eval(EVAL_POINTS[2]),
        cheb.eval(EVAL_POINTS[3]),
        derivative.coeffs()[0],
        derivative.coeffs()[1],
        derivative.coeffs()[2],
        derivative.eval(0.5),
        cheb.integral(),
    ];
    for (index, (computed, (output, reference))) in
        computed.into_iter().zip(KNOWN_OUTPUTS).enumerate()
    {
        if computed.to_bits() != reference.to_bits() {
            return CaseOutcome::fail(format!(
                "operation=Cheb1::from_coeffs+eval+differentiate+integral; clenshaw={CLENSHAW_CONVENTION}; differentiate={DERIVATIVE_CONVENTION}; integral={INTEGRAL_CONVENTION}; domain=[-1,1]; stored_coeffs=[2,3,0,5]; output={output}; output_index={index}; computed_bits=0x{:016x}; reference_bits=0x{:016x}; computed={computed}; reference={reference}",
                computed.to_bits(),
                reference.to_bits(),
            ))
            .with_evidence("crates/fs-cheb/CONTRACT.md#public-types-and-semantics")
            .with_evidence("crates/fs-cheb/CONTRACT.md#invariants");
        }
    }

    CaseOutcome::pass(
        "domain=[-1,1]; stored_coeffs=[2,3,0,5]; convention=half-c0; eval=[-7,1,-2.5,9]; derivative_coeffs=[36,0,30]; derivative_at_0.5=3; integral=2",
    )
    .with_evidence("crates/fs-cheb/CONTRACT.md#public-types-and-semantics")
    .with_evidence("crates/fs-cheb/CONTRACT.md#invariants")
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CheckpointStop {
    RefusedAt(usize),
}

fn checkpoint_outcome() -> CaseOutcome {
    let cheb = fixture();
    let baseline = cheb.eval(CHECKPOINT_POINT);
    let mut success_polls = 0_usize;
    let checked: Result<f64, &'static str> =
        cheb.eval_with_checkpoint(CHECKPOINT_POINT, &mut || {
            success_polls += 1;
            Ok(())
        });
    let checked = match checked {
        Ok(value) => value,
        Err(reason) => {
            return CaseOutcome::fail(format!(
                "operation=Cheb1::eval_with_checkpoint; point={CHECKPOINT_POINT}; unexpected_success_error={reason}; success_polls={success_polls}; expected_polls={EXPECTED_POLLS}"
            ))
            .with_evidence("crates/fs-cheb/CONTRACT.md#cancellation-behavior");
        }
    };
    if checked.to_bits() != CHECKPOINT_VALUE.to_bits()
        || checked.to_bits() != baseline.to_bits()
        || success_polls != EXPECTED_POLLS
    {
        return CaseOutcome::fail(format!(
            "operation=Cheb1::eval_with_checkpoint; clenshaw={CLENSHAW_CONVENTION}; point={CHECKPOINT_POINT}; checked_bits=0x{:016x}; baseline_bits=0x{:016x}; reference_bits=0x{:016x}; checked={checked}; reference={CHECKPOINT_VALUE}; success_polls={success_polls}; expected_polls={EXPECTED_POLLS}",
            checked.to_bits(),
            baseline.to_bits(),
            CHECKPOINT_VALUE.to_bits(),
        ))
        .with_evidence("crates/fs-cheb/CONTRACT.md#public-types-and-semantics")
        .with_evidence("crates/fs-cheb/CONTRACT.md#cancellation-behavior");
    }

    let mut refusal_polls = 0_usize;
    let refused = cheb.eval_with_checkpoint(CHECKPOINT_POINT, &mut || {
        refusal_polls += 1;
        if refusal_polls == REFUSAL_POLL {
            Err(CheckpointStop::RefusedAt(refusal_polls))
        } else {
            Ok(())
        }
    });
    if refused != Err(CheckpointStop::RefusedAt(REFUSAL_POLL)) || refusal_polls != REFUSAL_POLL {
        return CaseOutcome::fail(format!(
            "operation=Cheb1::eval_with_checkpoint; point={CHECKPOINT_POINT}; refusal_result={refused:?}; refusal_polls={refusal_polls}; expected_result=Err(CheckpointStop::RefusedAt({REFUSAL_POLL})); expected_refusal_polls={REFUSAL_POLL}"
        ))
        .with_evidence("crates/fs-cheb/CONTRACT.md#cancellation-behavior");
    }

    CaseOutcome::pass(format!(
        "point={CHECKPOINT_POINT}; value={CHECKPOINT_VALUE}; success_polls={success_polls}; refusal_poll={refusal_polls}; typed_refusal=CheckpointStop::RefusedAt({REFUSAL_POLL}); partial_value_published=false"
    ))
    .with_evidence("crates/fs-cheb/CONTRACT.md#public-types-and-semantics")
    .with_evidence("crates/fs-cheb/CONTRACT.md#cancellation-behavior")
}

fn lobatto_outcome() -> CaseOutcome {
    let points = lobatto_points(1);
    let matrix = diff_matrix(1);
    let expected_points: [f64; 2] = [1.0, -1.0];
    let expected_matrix: [f64; 4] = [0.5, -0.5, 0.5, -0.5];
    if points.len() != expected_points.len() || matrix.len() != expected_matrix.len() {
        return CaseOutcome::fail(format!(
            "operation=lobatto_points+diff_matrix; n=1; computed_point_count={}; reference_point_count={}; computed_matrix_count={}; reference_matrix_count={}",
            points.len(),
            expected_points.len(),
            matrix.len(),
            expected_matrix.len(),
        ))
        .with_evidence("crates/fs-cheb/CONTRACT.md#invariants");
    }
    for (index, (&computed, reference)) in points.iter().zip(expected_points).enumerate() {
        if computed.to_bits() != reference.to_bits() {
            return CaseOutcome::fail(format!(
                "operation=lobatto_points; n=1; ordering=descending; point_index={index}; computed_bits=0x{:016x}; reference_bits=0x{:016x}; computed={computed}; reference={reference}",
                computed.to_bits(),
                reference.to_bits(),
            ))
            .with_evidence("crates/fs-cheb/CONTRACT.md#public-types-and-semantics");
        }
    }
    for (index, (&computed, reference)) in matrix.iter().zip(expected_matrix).enumerate() {
        if computed.to_bits() != reference.to_bits() {
            return CaseOutcome::fail(format!(
                "operation=diff_matrix; n=1; ordering=row-major; matrix_index={index}; row={}; column={}; computed_bits=0x{:016x}; reference_bits=0x{:016x}; computed={computed}; reference={reference}",
                index / 2,
                index % 2,
                computed.to_bits(),
                reference.to_bits(),
            ))
            .with_evidence("crates/fs-cheb/CONTRACT.md#public-types-and-semantics")
            .with_evidence("crates/fs-cheb/CONTRACT.md#invariants");
        }
    }
    let matrix_bits = [
        matrix[0].to_bits(),
        matrix[1].to_bits(),
        matrix[2].to_bits(),
        matrix[3].to_bits(),
    ];
    let row_zeros = [matrix[1] + matrix[0], matrix[2] + matrix[3]];
    for (row, zero) in row_zeros.into_iter().enumerate() {
        if zero.to_bits() != 0.0_f64.to_bits() {
            return CaseOutcome::fail(format!(
                "operation=diff_matrix-negative-sum; n=1; row={row}; computed_zero_bits=0x{:016x}; reference_zero_bits=0x{:016x}; row_major_matrix_bits={matrix_bits:016x?}",
                zero.to_bits(),
                0.0_f64.to_bits(),
            ))
            .with_evidence("crates/fs-cheb/CONTRACT.md#invariants");
        }
    }

    CaseOutcome::pass(
        "n=1; descending_points=[1,-1]; row_major_diff_matrix=[0.5,-0.5,0.5,-0.5]; negative_sum_row_zero_bits=[0000000000000000,0000000000000000]",
    )
    .with_evidence("crates/fs-cheb/CONTRACT.md#public-types-and-semantics")
    .with_evidence("crates/fs-cheb/CONTRACT.md#invariants")
}

#[test]
fn bedrock_casebook_suite_emits_replay_complete_green_records() {
    let known_digest = fnv1a64(&known_answer_inputs());
    let checkpoint_digest = fnv1a64(&checkpoint_inputs());
    let lobatto_digest = fnv1a64(&lobatto_inputs());
    assert_eq!(known_digest, 0x117a_6c1c_8ac1_b144);
    assert_eq!(checkpoint_digest, 0xea1c_d25c_b116_3230);
    assert_eq!(lobatto_digest, 0x2f8f_d300_ac3b_af3e);

    let report = Suite::new(SUITE)
        .case(
            "clenshaw-calculus-known-answers",
            known_digest,
            ToleranceSpec::Exact,
            known_answer_outcome,
        )
        .case(
            "eval-checkpoint-refusal-policy",
            checkpoint_digest,
            ToleranceSpec::Structural,
            checkpoint_outcome,
        )
        .case(
            "lobatto-negative-sum-kat",
            lobatto_digest,
            ToleranceSpec::Exact,
            lobatto_outcome,
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
            "clenshaw-calculus-known-answers",
            "eval-checkpoint-refusal-policy",
            "lobatto-negative-sum-kat",
        ]
    );
    assert_eq!(
        report.records[0].json_line(),
        format!(
            concat!(
                "{{\"casebook\":{},\"suite\":\"fs-cheb/bedrock-conformance-v1\",",
                "\"case\":\"clenshaw-calculus-known-answers\",\"inputs_digest\":\"117a6c1c8ac1b144\",",
                "\"tolerance\":\"exact\",\"pass\":true,",
                "\"details\":\"domain=[-1,1]; stored_coeffs=[2,3,0,5]; convention=half-c0; eval=[-7,1,-2.5,9]; derivative_coeffs=[36,0,30]; derivative_at_0.5=3; integral=2\",",
                "\"evidence\":[\"crates/fs-cheb/CONTRACT.md#public-types-and-semantics\",",
                "\"crates/fs-cheb/CONTRACT.md#invariants\"]}}"
            ),
            CASEBOOK_RECORD_VERSION,
        ),
        "the structured spectral-core record schema and field order are contract"
    );
}

#[test]
fn disclosed_seeded_corruption_turns_the_casebook_suite_red() {
    const CORRUPTION_SEED: u64 = 0xF5CE_0001;
    let component = (CORRUPTION_SEED & 0x3) as usize;
    let bit = CORRUPTION_SEED.trailing_zeros();
    assert_eq!(component, 1);
    assert_eq!(bit, 0);
    let canonical = [
        (-7.0_f64).to_bits(),
        1.0_f64.to_bits(),
        (-2.5_f64).to_bits(),
        9.0_f64.to_bits(),
    ];
    let mut corrupted = canonical;
    corrupted[component] ^= 1_u64 << bit;
    let component_name = KNOWN_OUTPUTS[component].0;

    let known = known_answer_inputs();
    let mut inputs = b"fs-cheb:seeded-known-answer-corruption:v1".to_vec();
    push_u64(&mut inputs, CORRUPTION_SEED);
    push_len(&mut inputs, component);
    push_u64(&mut inputs, u64::from(bit));
    push_text(&mut inputs, "nested-known-answer-frame");
    push_len(&mut inputs, known.len());
    inputs.extend_from_slice(&known);
    push_text(&mut inputs, "canonical-eval-bits[-1,0,0.5,1]");
    push_u64s(&mut inputs, &canonical);
    push_text(&mut inputs, "corrupted-eval-bits[-1,0,0.5,1]");
    push_u64s(&mut inputs, &corrupted);
    let inputs_digest = fnv1a64(&inputs);
    assert_eq!(inputs_digest, 0x0b15_e93c_2468_1f16);

    let report = Suite::new(SUITE)
        .case(
            "seeded-known-answer-corruption",
            inputs_digest,
            ToleranceSpec::Exact,
            move || {
                let cheb = fixture();
                let computed = EVAL_POINTS.map(|point| cheb.eval(point).to_bits());
                if computed == corrupted {
                    CaseOutcome::pass("seeded corruption was not detected")
                } else {
                    CaseOutcome::fail(format!(
                        "seed=0x{CORRUPTION_SEED:016x}; operation=Cheb1::eval; clenshaw={CLENSHAW_CONVENTION}; domain=[-1,1]; stored_coeffs=[2,3,0,5]; eval_points=[-1,0,0.5,1]; component={component}; component_name={component_name}; bit={bit}; computed={computed:016x?}; canonical={canonical:016x?}; corrupted={corrupted:016x?}"
                    ))
                    .with_evidence("crates/fs-cheb/tests/conformance.rs#seeded-corruption")
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
    assert_eq!(failure.case, "seeded-known-answer-corruption");
    assert_eq!(failure.inputs_digest, "0b15e93c24681f16");
    assert!(
        failure
            .details
            .contains(&format!("seed=0x{CORRUPTION_SEED:016x}"))
    );
    assert!(failure.details.contains(&format!("component={component}")));
    assert!(
        failure
            .details
            .contains(&format!("component_name={component_name}; bit={bit}"))
    );
    assert!(failure.details.contains("computed=["));
    assert!(failure.details.contains("canonical=["));
    assert!(failure.details.contains("corrupted=["));
    assert!(
        failure
            .json_line()
            .contains("\"tolerance\":\"exact\",\"pass\":false")
    );

    let panic = std::panic::catch_unwind(|| report.assert_green())
        .expect_err("the merge-gate assertion must reject the seeded failure");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("casebook panic carries text");
    assert!(message.contains("seeded-known-answer-corruption"));
    assert!(message.contains(&format!("seed=0x{CORRUPTION_SEED:016x}")));
}
