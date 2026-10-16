//! Structured SUBSTRATE/BEDROCK conformance records for fs-math.
//!
//! The crate's larger numerical batteries retain broad ULP and algebraic
//! coverage. These cases expose three load-bearing deterministic laws through
//! the shared fs-casebook schema so a failing central run is replayable from
//! its JSON line alone (Gauntlet G0/G5 instrumentation).

use core::fmt::Write as _;
use fs_casebook::{CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, ToleranceSpec, fnv1a64};
use fs_math::eft::{quick_two_sum, two_prod, two_sum};
use fs_math::{canonical_nan, det, next_down, next_up, nudge_out, ulp_distance};

const SUITE: &str = "fs-math/bedrock-conformance-v1";
const CORE_SEED: u64 = 0x0DD_BA11;
const CORE_ITERATIONS: u64 = 25_000;
const CORE_GOLDEN: u64 = 0xeb79_cab7_a016_43e5;
const LCG_MULTIPLIER: u64 = 6_364_136_223_846_793_005;
const LCG_INCREMENT: u64 = 1_442_695_040_888_963_407;

#[derive(Clone, Copy)]
struct NudgeVector {
    input: u64,
    down: u64,
    up: u64,
    same_sign_span: Option<u64>,
}

const NUDGE_VECTORS: [NudgeVector; 4] = [
    NudgeVector {
        input: 0x0000_0000_0000_0000,
        down: 0x8000_0000_0000_0001,
        up: 0x0000_0000_0000_0001,
        same_sign_span: None,
    },
    NudgeVector {
        input: 0x3ff0_0000_0000_0000,
        down: 0x3fef_ffff_ffff_ffff,
        up: 0x3ff0_0000_0000_0001,
        same_sign_span: Some(2),
    },
    NudgeVector {
        input: 0xbff0_0000_0000_0000,
        down: 0xbff0_0000_0000_0001,
        up: 0xbfef_ffff_ffff_ffff,
        same_sign_span: Some(2),
    },
    NudgeVector {
        input: 0x0010_0000_0000_0000,
        down: 0x000f_ffff_ffff_ffff,
        up: 0x0010_0000_0000_0001,
        same_sign_span: Some(2),
    },
];

#[derive(Clone, Copy)]
struct EftVector {
    a: u64,
    b: u64,
    rounded: u64,
    error: u64,
}

const SUM_VECTORS: [EftVector; 1] = [EftVector {
    a: 0x4341_c379_37e0_8000,
    b: 0x3ff0_0000_0000_0000,
    rounded: 0x4341_c379_37e0_8000,
    error: 0x3ff0_0000_0000_0000,
}];

const PRODUCT_VECTORS: [EftVector; 2] = [
    EftVector {
        a: 0x4008_0000_0000_0000,
        b: 0x3fe0_0000_0000_0000,
        rounded: 0x3ff8_0000_0000_0000,
        error: 0x0000_0000_0000_0000,
    },
    EftVector {
        a: 0x3ff0_0000_0000_0001,
        b: 0x3ff0_0000_0000_0001,
        rounded: 0x3ff0_0000_0000_0002,
        error: 0x3970_0000_0000_0000,
    },
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

fn hex_bytes(bytes: &[u8]) -> String {
    let mut encoded = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        write!(&mut encoded, "{byte:02x}").expect("writing to String cannot fail");
    }
    encoded
}

fn core_inputs() -> Vec<u8> {
    let mut bytes = b"fs-math:strict-core-cross-isa:v1".to_vec();
    for value in [
        CORE_SEED,
        LCG_MULTIPLIER,
        LCG_INCREMENT,
        11,
        1_u64 << 53,
        CORE_ITERATIONS,
    ] {
        push_u64(&mut bytes, value);
    }
    let functions = ["exp", "expm1", "ln", "sin", "cos", "tanh", "sqrt"];
    push_len(&mut bytes, functions.len());
    for function in functions {
        push_text(&mut bytes, function);
    }
    let scales = [0.5_f64, 1.0e4, 0.05, 1.0e-9, 0.01, 1.0e6];
    push_len(&mut bytes, scales.len());
    for scale in scales {
        push_u64(&mut bytes, scale.to_bits());
    }
    push_u64(&mut bytes, CORE_GOLDEN);
    bytes
}

fn policy_inputs() -> Vec<u8> {
    let mut bytes = b"fs-math:ieee-policy-and-nudge:v1".to_vec();
    push_u64(&mut bytes, 0x7ff8_0000_0000_0000);
    push_len(&mut bytes, NUDGE_VECTORS.len());
    for vector in NUDGE_VECTORS {
        for value in [vector.input, vector.down, vector.up] {
            push_u64(&mut bytes, value);
        }
        bytes.push(u8::from(vector.same_sign_span.is_some()));
        if let Some(span) = vector.same_sign_span {
            push_u64(&mut bytes, span);
        }
    }
    bytes
}

fn eft_inputs() -> Vec<u8> {
    let mut bytes = b"fs-math:eft-known-answers:v1".to_vec();
    push_len(&mut bytes, SUM_VECTORS.len());
    for vector in SUM_VECTORS {
        for value in [vector.a, vector.b, vector.rounded, vector.error] {
            push_u64(&mut bytes, value);
        }
    }
    push_len(&mut bytes, PRODUCT_VECTORS.len());
    for vector in PRODUCT_VECTORS {
        for value in [vector.a, vector.b, vector.rounded, vector.error] {
            push_u64(&mut bytes, value);
        }
    }
    bytes
}

fn seeded_corruption_inputs(seed: u64, bit: u32, corrupted: u64) -> Vec<u8> {
    let core = core_inputs();
    let mut bytes = b"fs-math:seeded-core-golden-corruption:v1".to_vec();
    push_u64(&mut bytes, seed);
    push_u64(&mut bytes, u64::from(bit));
    push_len(&mut bytes, core.len());
    bytes.extend_from_slice(&core);
    push_u64(&mut bytes, CORE_GOLDEN);
    push_u64(&mut bytes, corrupted);
    bytes
}

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(LCG_MULTIPLIER)
        .wrapping_add(LCG_INCREMENT);
    ((*seed >> 11) as f64) / ((1_u64 << 53) as f64)
}

fn strict_core_fingerprint() -> u64 {
    let mut outputs = Vec::new();
    let mut seed = CORE_SEED;
    for _ in 0..CORE_ITERATIONS {
        let t = lcg(&mut seed);
        let x = (t - 0.5) * 1.0e4;
        for value in [
            det::exp(x * 0.05),
            det::expm1(x * 0.05),
            det::ln(t + 1.0e-9),
            det::sin(x),
            det::cos(x),
            det::tanh(x * 0.01),
            det::sqrt(t * 1.0e6),
        ] {
            outputs.extend_from_slice(&value.to_bits().to_le_bytes());
        }
    }
    fnv1a64(&outputs)
}

fn core_golden_outcome() -> CaseOutcome {
    let inputs_hex = hex_bytes(&core_inputs());
    let computed = strict_core_fingerprint();
    if computed != CORE_GOLDEN {
        return CaseOutcome::fail(format!(
            "seed=0x{CORE_SEED:016x}; iterations={CORE_ITERATIONS}; functions=7; computed=0x{computed:016x}; reference=0x{CORE_GOLDEN:016x}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-math/CONTRACT.md#determinism-class");
    }
    CaseOutcome::pass(format!(
        "seed=0x{CORE_SEED:016x}; iterations={CORE_ITERATIONS}; functions=7; fingerprint=0x{computed:016x}"
    ))
    .with_evidence("crates/fs-math/CONTRACT.md#determinism-class")
}

fn policy_outcome() -> CaseOutcome {
    let inputs_hex = hex_bytes(&policy_inputs());
    let canonical = canonical_nan().to_bits();
    if canonical != 0x7ff8_0000_0000_0000 {
        return CaseOutcome::fail(format!(
            "operation=canonical_nan; computed=0x{canonical:016x}; reference=0x7ff8000000000000; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-math/CONTRACT.md#invariants");
    }
    for (index, vector) in NUDGE_VECTORS.into_iter().enumerate() {
        let input = f64::from_bits(vector.input);
        let direct = (next_down(input).to_bits(), next_up(input).to_bits());
        let nudged = nudge_out(input);
        let nudged_bits = (nudged.0.to_bits(), nudged.1.to_bits());
        let span = vector
            .same_sign_span
            .map(|_| ulp_distance(nudged.0, nudged.1));
        if direct != (vector.down, vector.up)
            || nudged_bits != (vector.down, vector.up)
            || span != vector.same_sign_span
        {
            return CaseOutcome::fail(format!(
                "vector={index}; input_bits=0x{:016x}; direct_bits=[0x{:016x},0x{:016x}]; nudge_bits=[0x{:016x},0x{:016x}]; expected_bits=[0x{:016x},0x{:016x}]; same_sign_span={span:?}; expected_same_sign_span={:?}; inputs_hex={inputs_hex}",
                vector.input,
                direct.0,
                direct.1,
                nudged_bits.0,
                nudged_bits.1,
                vector.down,
                vector.up,
                vector.same_sign_span,
            ))
            .with_evidence("crates/fs-math/CONTRACT.md#invariants");
        }
    }
    CaseOutcome::pass(
        "canonical_nan=0x7ff8000000000000; nudge_vectors=4; same_sign_span_vectors=3; exact_bit_mismatches=0",
    )
    .with_evidence("crates/fs-math/CONTRACT.md#invariants")
}

fn eft_outcome() -> CaseOutcome {
    let inputs_hex = hex_bytes(&eft_inputs());
    for (index, vector) in SUM_VECTORS.into_iter().enumerate() {
        let (a, b) = (f64::from_bits(vector.a), f64::from_bits(vector.b));
        let sum = two_sum(a, b);
        let quick = quick_two_sum(a, b);
        let observed = (
            sum.0.to_bits(),
            sum.1.to_bits(),
            quick.0.to_bits(),
            quick.1.to_bits(),
        );
        let expected = (vector.rounded, vector.error, vector.rounded, vector.error);
        if observed != expected {
            return CaseOutcome::fail(format!(
                "sum_vector={index}; a_bits=0x{:016x}; b_bits=0x{:016x}; observed=[0x{:016x},0x{:016x},0x{:016x},0x{:016x}]; expected=[0x{:016x},0x{:016x},0x{:016x},0x{:016x}]; inputs_hex={inputs_hex}",
                vector.a,
                vector.b,
                observed.0,
                observed.1,
                observed.2,
                observed.3,
                expected.0,
                expected.1,
                expected.2,
                expected.3,
            ))
            .with_evidence("crates/fs-math/CONTRACT.md#public-types-and-semantics");
        }
    }
    for (index, vector) in PRODUCT_VECTORS.into_iter().enumerate() {
        let product = two_prod(f64::from_bits(vector.a), f64::from_bits(vector.b));
        let observed = (product.0.to_bits(), product.1.to_bits());
        let expected = (vector.rounded, vector.error);
        if observed != expected {
            return CaseOutcome::fail(format!(
                "product_vector={index}; a_bits=0x{:016x}; b_bits=0x{:016x}; observed=[0x{:016x},0x{:016x}]; expected=[0x{:016x},0x{:016x}]; inputs_hex={inputs_hex}",
                vector.a,
                vector.b,
                observed.0,
                observed.1,
                expected.0,
                expected.1,
            ))
            .with_evidence("crates/fs-math/CONTRACT.md#public-types-and-semantics");
        }
    }
    CaseOutcome::pass("sum_vectors=1; product_vectors=2; exact_bit_mismatches=0")
        .with_evidence("crates/fs-math/CONTRACT.md#public-types-and-semantics")
}

#[test]
fn bedrock_casebook_suite_emits_replay_complete_green_records() {
    assert_eq!(CASEBOOK_RECORD_VERSION, 1);
    let core_digest = fnv1a64(&core_inputs());
    let policy_digest = fnv1a64(&policy_inputs());
    let eft_digest = fnv1a64(&eft_inputs());
    assert_eq!(core_digest, 0x3526_9174_11ef_9b26);
    assert_eq!(policy_digest, 0x876c_1adc_d949_ee59);
    assert_eq!(eft_digest, 0xde4e_84b6_88f7_5a88);

    let report = Suite::new(SUITE)
        .case(
            "strict-core-cross-isa-golden",
            core_digest,
            ToleranceSpec::Exact,
            core_golden_outcome,
        )
        .case(
            "ieee-policy-and-nudge-known-answers",
            policy_digest,
            ToleranceSpec::Exact,
            policy_outcome,
        )
        .case(
            "eft-known-answers",
            eft_digest,
            ToleranceSpec::Exact,
            eft_outcome,
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
            "strict-core-cross-isa-golden",
            "ieee-policy-and-nudge-known-answers",
            "eft-known-answers",
        ]
    );
    assert_eq!(
        report.records[0].json_line(),
        format!(
            concat!(
                "{{\"casebook\":{},\"suite\":\"fs-math/bedrock-conformance-v1\",",
                "\"case\":\"strict-core-cross-isa-golden\",",
                "\"inputs_digest\":\"3526917411ef9b26\",\"tolerance\":\"exact\",",
                "\"pass\":true,\"details\":\"seed=0x0000000000ddba11; iterations=25000; functions=7; fingerprint=0xeb79cab7a01643e5\",",
                "\"evidence\":[\"crates/fs-math/CONTRACT.md#determinism-class\"]}}"
            ),
            CASEBOOK_RECORD_VERSION,
        ),
        "the structured record schema and field order are contract"
    );
}

#[test]
fn disclosed_seeded_corruption_turns_the_casebook_suite_red() {
    const CORRUPTION_SEED: u64 = 0xF5A7_0001;
    let bit = CORRUPTION_SEED.trailing_zeros();
    assert_eq!(bit, 0);
    let corrupted = CORE_GOLDEN ^ (1_u64 << bit);
    assert_eq!(corrupted, 0xeb79_cab7_a016_43e4);
    let inputs = seeded_corruption_inputs(CORRUPTION_SEED, bit, corrupted);
    let inputs_hex = hex_bytes(&inputs);
    let inputs_digest = fnv1a64(&inputs);
    assert_eq!(inputs_digest, 0xe1f9_4975_c093_f753);

    let report = Suite::new(SUITE)
        .case(
            "seeded-core-golden-corruption",
            inputs_digest,
            ToleranceSpec::Exact,
            move || {
                let computed = strict_core_fingerprint();
                if computed == corrupted {
                    CaseOutcome::pass("seeded core corruption escaped detection")
                } else {
                    CaseOutcome::fail(format!(
                        "seed=0x{CORRUPTION_SEED:016x}; bit={bit}; computed=0x{computed:016x}; canonical=0x{CORE_GOLDEN:016x}; corrupted=0x{corrupted:016x}; inputs_hex={inputs_hex}"
                    ))
                    .with_evidence(
                        "crates/fs-math/tests/conformance.rs#disclosed-seeded-corruption",
                    )
                }
            },
        )
        .run();

    assert!(!report.all_passed(), "the corrupted oracle must turn red");
    let failures = report.failures();
    let [failure] = failures.as_slice() else {
        panic!("the seeded corruption must produce exactly one structured failure");
    };
    assert_eq!(failure.case, "seeded-core-golden-corruption");
    assert_eq!(failure.inputs_digest, "e1f94975c093f753");
    assert!(
        failure
            .details
            .contains(&format!("seed=0x{CORRUPTION_SEED:016x}"))
    );
    assert!(failure.details.contains("bit=0"));
    assert!(failure.details.contains("computed=0xeb79cab7a01643e5"));
    assert!(failure.details.contains("canonical=0xeb79cab7a01643e5"));
    assert!(failure.details.contains("corrupted=0xeb79cab7a01643e4"));
    assert!(failure.details.contains("inputs_hex="));
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
    assert!(message.contains("seeded-core-golden-corruption"));
    assert!(message.contains(&format!("seed=0x{CORRUPTION_SEED:016x}")));
}
