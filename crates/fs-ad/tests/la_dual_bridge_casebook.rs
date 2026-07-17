//! Structured cross-crate evidence for the fs-ad -> fs-la scalar bridge.
//!
//! The fixed fixtures prove that packed forward-dual lanes execute inside the
//! public fs-la reference GEMM, including exact-zero and transactional shape
//! policy. They do not claim that the optimized f64 microkernel is generic,
//! SIMD throughput for duals, or fresh cross-ISA execution.

use core::fmt::Write as _;

use fs_ad::{Dual, Dual64};
use fs_casebook::{CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, ToleranceSpec, fnv1a64};
use fs_la::{GEMM_SCALAR_SEMANTICS_VERSION, GemmShapeError, gemm_f64, gemm_scalar_checked};

const SUITE: &str = "bedrock/fs-ad-fs-la-dual-gemm-v1";
const MATRIX_DIMENSION: usize = 2;
const POISON_BITS: u64 = 0x7ff8_0000_0000_0000;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct KatMeasurement {
    packed: [[u64; 3]; 4],
    x_single: [[u64; 2]; 4],
    y_single: [[u64; 2]; 4],
    scalar: [u64; 4],
}

const EXPECTED_KAT: KatMeasurement = KatMeasurement {
    packed: [
        [0x4014_0000_0000_0000, 0x4010_0000_0000_0000, 0],
        [0, 0xc000_0000_0000_0000, 0],
        [0xc004_0000_0000_0000, 0, 0x3fe0_0000_0000_0000],
        [0x4014_0000_0000_0000, 0, 0x3ff0_0000_0000_0000],
    ],
    x_single: [
        [0x4014_0000_0000_0000, 0x4010_0000_0000_0000],
        [0, 0xc000_0000_0000_0000],
        [0xc004_0000_0000_0000, 0],
        [0x4014_0000_0000_0000, 0],
    ],
    y_single: [
        [0x4014_0000_0000_0000, 0],
        [0, 0],
        [0xc004_0000_0000_0000, 0x3fe0_0000_0000_0000],
        [0x4014_0000_0000_0000, 0x3ff0_0000_0000_0000],
    ],
    scalar: [
        0x4014_0000_0000_0000,
        0,
        0xc004_0000_0000_0000,
        0x4014_0000_0000_0000,
    ],
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PolicyMeasurement {
    zero_primal_alpha: [u64; 2],
    exact_zero_alpha: [u64; 2],
    nested_square: [u64; 4],
    shape_result: Result<(), GemmShapeError>,
    output_after_refusal: [u64; 2],
}

const EXPECTED_POLICY: PolicyMeasurement = PolicyMeasurement {
    zero_primal_alpha: [0, 0x4018_0000_0000_0000],
    exact_zero_alpha: [0x4018_0000_0000_0000, 0],
    nested_square: [
        0x4010_0000_0000_0000,
        0x4010_0000_0000_0000,
        0x4010_0000_0000_0000,
        0x4000_0000_0000_0000,
    ],
    shape_result: Err(GemmShapeError::LengthMismatch {
        operand: "a",
        expected: 2,
        actual: 1,
    }),
    output_after_refusal: [0x401c_0000_0000_0000, 0x4026_0000_0000_0000],
};

#[derive(Debug, Clone, Copy)]
struct Corruption {
    seed: u64,
    output: usize,
    derivative_lane: usize,
}

fn push_u64(bytes: &mut Vec<u8>, value: u64) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn push_len(bytes: &mut Vec<u8>, value: usize) {
    push_u64(
        bytes,
        u64::try_from(value).expect("bridge fixture lengths fit u64"),
    );
}

fn push_text(bytes: &mut Vec<u8>, value: &str) {
    push_len(bytes, value.len());
    bytes.extend_from_slice(value.as_bytes());
}

fn push_f64_bits(bytes: &mut Vec<u8>, values: &[u64]) {
    push_len(bytes, values.len());
    for value in values {
        push_u64(bytes, *value);
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

fn push_kat_measurement(bytes: &mut Vec<u8>, measurement: KatMeasurement) {
    for output in measurement.packed {
        push_f64_bits(bytes, &output);
    }
    for output in measurement.x_single {
        push_f64_bits(bytes, &output);
    }
    for output in measurement.y_single {
        push_f64_bits(bytes, &output);
    }
    push_f64_bits(bytes, &measurement.scalar);
}

fn kat_inputs() -> Vec<u8> {
    let mut bytes = b"bedrock:fs-ad-fs-la:packed-dual-gemm:v1".to_vec();
    push_text(
        &mut bytes,
        "fs_la::gemm_scalar_checked<fs_ad::Dual64<N>>+fs_la::gemm_f64",
    );
    push_u64(&mut bytes, u64::from(GEMM_SCALAR_SEMANTICS_VERSION));
    push_text(&mut bytes, "shape=m2-n2-k2-row-major");
    push_text(&mut bytes, "packed-a=re,dx,dy");
    for input in [
        [1.0_f64.to_bits(), 1.0_f64.to_bits(), 0],
        [2.0_f64.to_bits(), 0, 0],
        [(-1.0_f64).to_bits(), 0, 0],
        [3.0_f64.to_bits(), 0, 1.0_f64.to_bits()],
    ] {
        push_f64_bits(&mut bytes, &input);
    }
    push_text(&mut bytes, "b-constants");
    push_f64_bits(
        &mut bytes,
        &[
            4.0_f64.to_bits(),
            (-2.0_f64).to_bits(),
            0.5_f64.to_bits(),
            1.0_f64.to_bits(),
        ],
    );
    push_text(&mut bytes, "alpha=(one,zero,zero)");
    push_f64_bits(&mut bytes, &[1.0_f64.to_bits(), 0, 0]);
    push_text(&mut bytes, "beta=(zero,zero,zero)");
    push_f64_bits(&mut bytes, &[0, 0, 0]);
    push_text(&mut bytes, "c=4x(re,dx,dy)=canonical-qnan");
    push_u64(&mut bytes, 4);
    push_f64_bits(&mut bytes, &[POISON_BITS; 3]);
    push_text(
        &mut bytes,
        "expected=packed+x-single+y-single+optimized-f64",
    );
    push_kat_measurement(&mut bytes, EXPECTED_KAT);
    push_text(
        &mut bytes,
        "fixed-i-j-k;packed-lanes-equal-single-runs;primal-equals-f64:v1",
    );
    bytes
}

fn policy_inputs() -> Vec<u8> {
    let mut bytes = b"bedrock:fs-ad-fs-la:dual-gemm-policy:v1".to_vec();
    push_text(
        &mut bytes,
        "GemmScalar::is_exact_zero+gemm_scalar_checked+GemmShapeError",
    );
    push_u64(&mut bytes, u64::from(GEMM_SCALAR_SEMANTICS_VERSION));
    push_text(
        &mut bytes,
        "zero-primal-alpha=(re=0,eps=1);a=2;b=3;beta=0;c=qnan",
    );
    push_f64_bits(
        &mut bytes,
        &[
            0,
            1.0_f64.to_bits(),
            2.0_f64.to_bits(),
            3.0_f64.to_bits(),
            POISON_BITS,
        ],
    );
    push_text(&mut bytes, "exact-zero-alpha;poison-a-b;beta=2;c=(3,0)");
    push_f64_bits(
        &mut bytes,
        &[
            0,
            0,
            POISON_BITS,
            POISON_BITS,
            2.0_f64.to_bits(),
            3.0_f64.to_bits(),
        ],
    );
    push_text(
        &mut bytes,
        "nested-dual-square:x=(re=(2,1),eps=[(1,0)]);expected=(4,4,4,2)",
    );
    push_f64_bits(
        &mut bytes,
        &[
            2.0_f64.to_bits(),
            1.0_f64.to_bits(),
            1.0_f64.to_bits(),
            0,
            4.0_f64.to_bits(),
            4.0_f64.to_bits(),
            4.0_f64.to_bits(),
            2.0_f64.to_bits(),
        ],
    );
    push_text(
        &mut bytes,
        "malformed=m1-n1-k2;a-len1;b-len2;c=(7,11);expected-a-length-mismatch",
    );
    push_f64_bits(
        &mut bytes,
        &[
            1.0_f64.to_bits(),
            1.0_f64.to_bits(),
            1.0_f64.to_bits(),
            7.0_f64.to_bits(),
            11.0_f64.to_bits(),
        ],
    );
    push_text(
        &mut bytes,
        "policy=full-scalar-zero;alpha-no-read;beta-no-read;shape-preflight-transactional:v1",
    );
    bytes
}

fn dual2(re: f64, dx: f64, dy: f64) -> Dual64<2> {
    Dual { re, eps: [dx, dy] }
}

fn dual1(re: f64, derivative: f64) -> Dual64<1> {
    Dual {
        re,
        eps: [derivative],
    }
}

fn dual2_bits(value: Dual64<2>) -> [u64; 3] {
    [
        value.re.to_bits(),
        value.eps[0].to_bits(),
        value.eps[1].to_bits(),
    ]
}

fn dual1_bits(value: Dual64<1>) -> [u64; 2] {
    [value.re.to_bits(), value.eps[0].to_bits()]
}

fn nested_bits(value: Dual<Dual64<1>, 1>) -> [u64; 4] {
    [
        value.re.re.to_bits(),
        value.re.eps[0].to_bits(),
        value.eps[0].re.to_bits(),
        value.eps[0].eps[0].to_bits(),
    ]
}

fn measure_kat() -> Result<KatMeasurement, GemmShapeError> {
    let packed_a = [
        dual2(1.0, 1.0, 0.0),
        dual2(2.0, 0.0, 0.0),
        dual2(-1.0, 0.0, 0.0),
        dual2(3.0, 0.0, 1.0),
    ];
    let packed_b = [4.0, -2.0, 0.5, 1.0].map(Dual64::<2>::constant);
    let poison_value = f64::from_bits(POISON_BITS);
    let poison2 = dual2(poison_value, poison_value, poison_value);
    let mut packed = [poison2; 4];
    gemm_scalar_checked(
        MATRIX_DIMENSION,
        MATRIX_DIMENSION,
        MATRIX_DIMENSION,
        Dual64::<2>::constant(1.0),
        &packed_a,
        &packed_b,
        Dual64::<2>::constant(0.0),
        &mut packed,
    )?;

    let x_a = [
        dual1(1.0, 1.0),
        dual1(2.0, 0.0),
        dual1(-1.0, 0.0),
        dual1(3.0, 0.0),
    ];
    let y_a = [
        dual1(1.0, 0.0),
        dual1(2.0, 0.0),
        dual1(-1.0, 0.0),
        dual1(3.0, 1.0),
    ];
    let single_b = [4.0, -2.0, 0.5, 1.0].map(Dual64::<1>::constant);
    let poison1 = dual1(poison_value, poison_value);
    let mut x_single = [poison1; 4];
    let mut y_single = [poison1; 4];
    for (a, output) in [(&x_a, &mut x_single), (&y_a, &mut y_single)] {
        gemm_scalar_checked(
            MATRIX_DIMENSION,
            MATRIX_DIMENSION,
            MATRIX_DIMENSION,
            Dual64::<1>::constant(1.0),
            a,
            &single_b,
            Dual64::<1>::constant(0.0),
            output,
        )?;
    }

    let scalar_a = [1.0, 2.0, -1.0, 3.0];
    let scalar_b = [4.0, -2.0, 0.5, 1.0];
    let mut scalar = [f64::from_bits(POISON_BITS); 4];
    gemm_f64(2, 2, 2, 1.0, &scalar_a, &scalar_b, 0.0, &mut scalar);

    Ok(KatMeasurement {
        packed: packed.map(dual2_bits),
        x_single: x_single.map(dual1_bits),
        y_single: y_single.map(dual1_bits),
        scalar: scalar.map(f64::to_bits),
    })
}

fn measure_policy() -> Result<PolicyMeasurement, GemmShapeError> {
    let poison_value = f64::from_bits(POISON_BITS);
    let poison = dual1(poison_value, poison_value);
    let mut derivative_output = [poison];
    gemm_scalar_checked(
        1,
        1,
        1,
        dual1(0.0, 1.0),
        &[Dual64::<1>::constant(2.0)],
        &[Dual64::<1>::constant(3.0)],
        Dual64::<1>::constant(0.0),
        &mut derivative_output,
    )?;

    let mut exact_zero_output = [Dual64::<1>::constant(3.0)];
    gemm_scalar_checked(
        1,
        1,
        1,
        Dual64::<1>::constant(0.0),
        &[poison],
        &[poison],
        Dual64::<1>::constant(2.0),
        &mut exact_zero_output,
    )?;

    let nested_x = Dual {
        re: dual1(2.0, 1.0),
        eps: [dual1(1.0, 0.0)],
    };
    let nested_zero = Dual::<Dual64<1>, 1>::constant(Dual64::<1>::constant(0.0));
    let nested_one = Dual::<Dual64<1>, 1>::constant(Dual64::<1>::constant(1.0));
    let mut nested_output = [nested_zero];
    gemm_scalar_checked(
        1,
        1,
        1,
        nested_one,
        &[nested_x],
        &[nested_x],
        nested_zero,
        &mut nested_output,
    )?;

    let mut refusal_output = [dual1(7.0, 11.0)];
    let shape_result = gemm_scalar_checked(
        1,
        1,
        2,
        Dual64::<1>::constant(1.0),
        &[Dual64::<1>::constant(1.0)],
        &[Dual64::<1>::constant(1.0); 2],
        Dual64::<1>::constant(0.0),
        &mut refusal_output,
    );

    Ok(PolicyMeasurement {
        zero_primal_alpha: dual1_bits(derivative_output[0]),
        exact_zero_alpha: dual1_bits(exact_zero_output[0]),
        nested_square: nested_bits(nested_output[0]),
        shape_result,
        output_after_refusal: dual1_bits(refusal_output[0]),
    })
}

fn kat_outcome(
    reference: KatMeasurement,
    corruption: Option<Corruption>,
    input_frame: &[u8],
) -> CaseOutcome {
    let run = measure_kat();
    let replay = measure_kat();
    let inputs_hex = hex_bytes(input_frame);
    let context = corruption.map_or_else(
        || "mode=canonical".to_owned(),
        |corruption| {
            format!(
                "seed=0x{:016x}; output={}; derivative_lane={}",
                corruption.seed, corruption.output, corruption.derivative_lane
            )
        },
    );
    if run != replay {
        return CaseOutcome::fail(format!(
            "{context}; stage=same-run-replay; first={run:?}; second={replay:?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-la/CONTRACT.md#determinism-class")
        .with_evidence("crates/fs-ad/CONTRACT.md#determinism-class");
    }
    let Ok(measurement) = run else {
        return CaseOutcome::fail(format!(
            "{context}; stage=packed-dual-gemm-admission; result={run:?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-la/CONTRACT.md#error-model");
    };
    if measurement != reference {
        return CaseOutcome::fail(format!(
            "{context}; stage=packed-dual-known-answer; computed={measurement:?}; reference={reference:?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-la/CONTRACT.md#public-types-and-semantics")
        .with_evidence("crates/fs-ad/CONTRACT.md#invariants");
    }
    for output in 0..4 {
        if measurement.packed[output][0] != measurement.scalar[output]
            || measurement.packed[output][1] != measurement.x_single[output][1]
            || measurement.packed[output][2] != measurement.y_single[output][1]
        {
            return CaseOutcome::fail(format!(
                "{context}; stage=packed-single-primal-equivalence; output={output}; measurement={measurement:?}; inputs_hex={inputs_hex}"
            ))
            .with_evidence("crates/fs-ad/CONTRACT.md#invariants");
        }
    }
    CaseOutcome::pass(
        "dual2_gemm=exact; packed_lanes=two_dual1_runs; primal=optimized_f64; beta_zero_poison=overwritten; same_run=identical",
    )
    .with_evidence("crates/fs-la/CONTRACT.md#public-types-and-semantics")
    .with_evidence("crates/fs-ad/CONTRACT.md#invariants")
}

fn policy_outcome(input_frame: &[u8]) -> CaseOutcome {
    let run = measure_policy();
    let replay = measure_policy();
    let inputs_hex = hex_bytes(input_frame);
    if run != replay {
        return CaseOutcome::fail(format!(
            "stage=same-run-replay; first={run:?}; second={replay:?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-la/CONTRACT.md#determinism-class");
    }
    let Ok(measurement) = run else {
        return CaseOutcome::fail(format!(
            "stage=valid-policy-probe-refused; result={run:?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-la/CONTRACT.md#error-model");
    };
    if measurement != EXPECTED_POLICY {
        return CaseOutcome::fail(format!(
            "stage=dual-zero-or-shape-policy; computed={measurement:?}; reference={EXPECTED_POLICY:?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-la/CONTRACT.md#error-model")
        .with_evidence("crates/fs-ad/CONTRACT.md#public-types-and-semantics");
    }
    CaseOutcome::pass(
        "zero_primal_nonzero_lane=executed; exact_zero_alpha=no_operand_read; exact_zero_beta=no_output_read; nested_dual=exact; malformed_a=typed_refusal; c_after_refusal=unchanged; same_run=identical",
    )
    .with_evidence("crates/fs-la/CONTRACT.md#error-model")
    .with_evidence("crates/fs-ad/CONTRACT.md#public-types-and-semantics")
}

#[test]
fn packed_dual_fs_la_bridge_emits_replay_complete_green_records() {
    assert_eq!(CASEBOOK_RECORD_VERSION, 1);
    assert_eq!(GEMM_SCALAR_SEMANTICS_VERSION, 1);
    let kat_frame = kat_inputs();
    let policy_frame = policy_inputs();
    let kat_digest = fnv1a64(&kat_frame);
    let policy_digest = fnv1a64(&policy_frame);
    assert_eq!(kat_digest, 0x176c_af7d_7086_7033);
    assert_eq!(policy_digest, 0x8ec4_5621_dec8_d560);

    let report = Suite::new(SUITE)
        .case(
            "packed-dual-gemm-known-answer",
            kat_digest,
            ToleranceSpec::Exact,
            move || kat_outcome(EXPECTED_KAT, None, &kat_frame),
        )
        .case(
            "dual-zero-and-shape-policy",
            policy_digest,
            ToleranceSpec::Exact,
            move || policy_outcome(&policy_frame),
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
            "packed-dual-gemm-known-answer",
            "dual-zero-and-shape-policy"
        ]
    );
    assert_eq!(
        report.records[0].json_line(),
        format!(
            concat!(
                "{{\"casebook\":{},\"suite\":\"bedrock/fs-ad-fs-la-dual-gemm-v1\",",
                "\"case\":\"packed-dual-gemm-known-answer\",\"inputs_digest\":\"176caf7d70867033\",",
                "\"tolerance\":\"exact\",\"pass\":true,",
                "\"details\":\"dual2_gemm=exact; packed_lanes=two_dual1_runs; primal=optimized_f64; beta_zero_poison=overwritten; same_run=identical\",",
                "\"evidence\":[\"crates/fs-la/CONTRACT.md#public-types-and-semantics\",",
                "\"crates/fs-ad/CONTRACT.md#invariants\"]}}"
            ),
            CASEBOOK_RECORD_VERSION,
        ),
        "the bridge record schema and field order are contract"
    );
}

#[test]
fn disclosed_seeded_bridge_reference_corruption_turns_the_suite_red() {
    const CORRUPTION_SEED: u64 = 0xF5A1_0000;
    let output = (CORRUPTION_SEED & 0x3) as usize;
    let derivative_lane = ((CORRUPTION_SEED >> 2) & 0x1) as usize;
    assert_eq!((output, derivative_lane), (0, 0));
    let corruption = Corruption {
        seed: CORRUPTION_SEED,
        output,
        derivative_lane,
    };
    let canonical = EXPECTED_KAT;
    let mut corrupted = canonical;
    corrupted.packed[output][derivative_lane + 1] ^= 1;
    assert_eq!(corrupted.packed[0][1], EXPECTED_KAT.packed[0][1] ^ 1);

    let kat_frame = kat_inputs();
    let mut inputs = b"bedrock:fs-ad-fs-la:seeded-dual-reference-corruption:v1".to_vec();
    push_u64(&mut inputs, CORRUPTION_SEED);
    push_len(&mut inputs, output);
    push_len(&mut inputs, derivative_lane);
    push_nested(&mut inputs, "nested-packed-dual-kat", &kat_frame);
    push_text(&mut inputs, "canonical-reference");
    push_kat_measurement(&mut inputs, canonical);
    push_text(&mut inputs, "corrupted-reference");
    push_kat_measurement(&mut inputs, corrupted);
    let inputs_digest = fnv1a64(&inputs);
    assert_eq!(inputs_digest, 0x1376_208c_ab99_b80b);

    let make_report = || {
        let input_frame = inputs.clone();
        Suite::new(SUITE)
            .case(
                "seeded-packed-dual-reference-corruption",
                inputs_digest,
                ToleranceSpec::Exact,
                move || kat_outcome(corrupted, Some(corruption), &input_frame),
            )
            .run()
    };
    let first = make_report();
    let replay = make_report();
    let first_failures = first.failures();
    let replay_failures = replay.failures();
    let [first_failure] = first_failures.as_slice() else {
        panic!("the disclosed corruption must produce exactly one failure");
    };
    let [replay_failure] = replay_failures.as_slice() else {
        panic!("the replayed corruption must produce exactly one failure");
    };
    assert_eq!(first_failure.json_line(), replay_failure.json_line());
    assert_eq!(first_failure.inputs_digest, "1376208cab99b80b");
    assert!(
        first_failure
            .details
            .contains("stage=packed-dual-known-answer")
    );
    assert!(
        first_failure
            .details
            .contains(&format!("seed=0x{CORRUPTION_SEED:016x}"))
    );
    assert!(first_failure.details.contains("inputs_hex="));
    assert!(
        first_failure
            .json_line()
            .contains("\"tolerance\":\"exact\",\"pass\":false")
    );

    let panic = std::panic::catch_unwind(|| first.assert_green())
        .expect_err("the merge gate must reject the disclosed corruption");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("casebook panic carries text");
    assert!(message.contains("seeded-packed-dual-reference-corruption"));
    assert!(message.contains(&format!("seed=0x{CORRUPTION_SEED:016x}")));
}
