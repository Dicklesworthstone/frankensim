//! G5 composition evidence for logical Philox inputs flowing through GEMM.
//!
//! The fixture is deliberately small but enters fs-la's real scoped-thread
//! path for every worker request above one. Matrix values are keyed only by
//! logical coordinates, then materialized under several simulated generation
//! partitions and traversal orders before serial and parallel GEMM are
//! compared bit-for-bit.
//!
//! This is same-build, same-ISA evidence for one disclosed finite fixture. It
//! is not fresh cross-ISA execution proof and makes no randomness-quality,
//! general-shape, performance, placement, NUMA, cancellation, or drain claim.

use core::fmt::Write as _;
use std::panic::{AssertUnwindSafe, catch_unwind};

use fs_casebook::{CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, ToleranceSpec, fnv1a64};
use fs_la::gemm::GEMM_BIT_SEMANTICS_VERSION;
use fs_la::{
    GEMM_IMPLEMENTATION_VERSION, VERSION as FS_LA_VERSION, gemm_execution_tier, gemm_f64,
    gemm_f64_parallel,
};
use fs_rand::{
    STREAM_POSITION_IDENTITY_DOMAIN, STREAM_SEMANTICS_VERSION, Stream, StreamKey,
    VERSION as FS_RAND_VERSION,
};

const SUITE: &str = "bedrock/fs-rand-to-fs-la-gemm-replay-v1";
const M: usize = 257;
const N: usize = 7;
const K: usize = 9;
const ALPHA: f64 = 1.25;
const BETA: f64 = -0.5;
const ROOT_SEED: u64 = 0x6A5E_5EED_2026_0717;
const A_KERNEL: u32 = 0x4745_4D41;
const B_KERNEL: u32 = 0x4745_4D42;
const C_KERNEL: u32 = 0x4745_4D43;
const PARALLEL_WORKERS: [usize; 5] = [1, 2, 3, 5, 8];
const RED_SEED: u64 = 0x6A5E_C0DE_0000_0011;
const GREEN_FRAME_LEN: usize = 34_572;
const GREEN_FRAME_DIGEST: u64 = 0xb490_a105_e7d5_32d4;
const RED_FRAME_LEN: usize = 35_177;
const RED_FRAME_DIGEST: u64 = 0x917d_cd5b_1fb6_7a25;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct GenerationPlan {
    partitions: usize,
    reverse_partitions: bool,
}

const GENERATION_PLANS: [GenerationPlan; 9] = [
    GenerationPlan {
        partitions: 1,
        reverse_partitions: false,
    },
    GenerationPlan {
        partitions: 2,
        reverse_partitions: false,
    },
    GenerationPlan {
        partitions: 2,
        reverse_partitions: true,
    },
    GenerationPlan {
        partitions: 3,
        reverse_partitions: false,
    },
    GenerationPlan {
        partitions: 3,
        reverse_partitions: true,
    },
    GenerationPlan {
        partitions: 5,
        reverse_partitions: false,
    },
    GenerationPlan {
        partitions: 5,
        reverse_partitions: true,
    },
    GenerationPlan {
        partitions: 8,
        reverse_partitions: false,
    },
    GenerationPlan {
        partitions: 8,
        reverse_partitions: true,
    },
];

#[derive(Debug, Clone, PartialEq, Eq)]
struct FixtureBits {
    a: Vec<u64>,
    b: Vec<u64>,
    c: Vec<u64>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Evaluation {
    fixture_digest: u64,
    serial: Vec<u64>,
    parallel: Vec<(usize, Vec<u64>)>,
}

fn push_u64(bytes: &mut Vec<u8>, value: u64) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn push_len(bytes: &mut Vec<u8>, value: usize) {
    push_u64(
        bytes,
        u64::try_from(value).expect("Casebook fixture lengths fit u64"),
    );
}

fn push_text(bytes: &mut Vec<u8>, value: &str) {
    push_len(bytes, value.len());
    bytes.extend_from_slice(value.as_bytes());
}

fn push_bits(bytes: &mut Vec<u8>, values: &[u64]) {
    push_len(bytes, values.len());
    for &value in values {
        push_u64(bytes, value);
    }
}

fn push_nested(bytes: &mut Vec<u8>, label: &str, nested: &[u8]) {
    push_text(bytes, label);
    push_len(bytes, nested.len());
    bytes.extend_from_slice(nested);
}

fn hex_bytes(bytes: &[u8]) -> String {
    let mut encoded = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        write!(&mut encoded, "{byte:02x}").expect("writing to String cannot fail");
    }
    encoded
}

fn digest_bits(parts: &[&[u64]]) -> u64 {
    let capacity = parts
        .iter()
        .map(|part| part.len())
        .sum::<usize>()
        .saturating_mul(8);
    let mut bytes = Vec::with_capacity(capacity);
    for part in parts {
        for &value in *part {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
    }
    fnv1a64(&bytes)
}

fn generated_value_bits(kernel: u32, row: usize, column: usize) -> u64 {
    let key = StreamKey {
        seed: ROOT_SEED,
        kernel,
        tile: u32::try_from(row).expect("fixture row fits the logical tile field"),
    };
    let block = Stream::at(
        key,
        u64::try_from(column).expect("fixture column fits the stream index"),
    );
    let word = (u64::from(block[1]) << 32) | u64::from(block[0]);
    let top53 = word >> 11;
    // Exact dyadic ladder in [-4, 4): integer conversion is exact through
    // 2^53, centering is exact, and scaling is a power of two.
    ((top53 as f64 - 4_503_599_627_370_496.0) * 8.881_784_197_001_252e-16).to_bits()
}

fn materialize_matrix(rows: usize, columns: usize, kernel: u32, plan: GenerationPlan) -> Vec<u64> {
    assert!(
        plan.partitions > 0,
        "generation partitions must be positive"
    );
    let mut values = vec![0; rows * columns];
    let mut fill_partition = |partition: usize| {
        for flat in (partition..values.len()).step_by(plan.partitions) {
            values[flat] = generated_value_bits(kernel, flat / columns, flat % columns);
        }
    };
    if plan.reverse_partitions {
        for partition in (0..plan.partitions).rev() {
            fill_partition(partition);
        }
    } else {
        for partition in 0..plan.partitions {
            fill_partition(partition);
        }
    }
    values
}

fn materialize_fixture(plan: GenerationPlan) -> FixtureBits {
    FixtureBits {
        a: materialize_matrix(M, K, A_KERNEL, plan),
        b: materialize_matrix(K, N, B_KERNEL, plan),
        c: materialize_matrix(M, N, C_KERNEL, plan),
    }
}

fn first_mismatch(left: &[u64], right: &[u64]) -> Option<(usize, u64, u64)> {
    left.iter()
        .zip(right)
        .enumerate()
        .find_map(|(index, (&lhs, &rhs))| (lhs != rhs).then_some((index, lhs, rhs)))
        .or_else(|| (left.len() != right.len()).then_some((left.len().min(right.len()), 0, 0)))
}

fn validate_generation() -> Result<FixtureBits, String> {
    let canonical = materialize_fixture(GENERATION_PLANS[0]);
    for plan in GENERATION_PLANS {
        let candidate = materialize_fixture(plan);
        for (matrix, reference, actual) in [
            ("A", canonical.a.as_slice(), candidate.a.as_slice()),
            ("B", canonical.b.as_slice(), candidate.b.as_slice()),
            ("C", canonical.c.as_slice(), candidate.c.as_slice()),
        ] {
            if let Some((index, expected, computed)) = first_mismatch(reference, actual) {
                return Err(format!(
                    "stage=logical-generation-replay; matrix={matrix}; partitions={}; reverse={}; index={index}; canonical_bits=0x{expected:016x}; computed_bits=0x{computed:016x}",
                    plan.partitions, plan.reverse_partitions,
                ));
            }
        }
    }
    Ok(canonical)
}

fn common_frame_prefix(domain: &[u8]) -> Vec<u8> {
    let mut bytes = domain.to_vec();
    push_text(&mut bytes, "encoding");
    push_text(
        &mut bytes,
        "length-prefixed-little-endian-u64-and-f64-bits:v1",
    );
    push_text(&mut bytes, "casebook-record-version");
    push_u64(&mut bytes, u64::from(CASEBOOK_RECORD_VERSION));
    push_text(&mut bytes, "fs-la-version");
    push_text(&mut bytes, FS_LA_VERSION);
    push_text(&mut bytes, "fs-rand-version");
    push_text(&mut bytes, FS_RAND_VERSION);
    push_text(&mut bytes, "gemm-bit-semantics-version");
    push_u64(&mut bytes, u64::from(GEMM_BIT_SEMANTICS_VERSION));
    push_text(&mut bytes, "gemm-implementation-version");
    push_u64(&mut bytes, u64::from(GEMM_IMPLEMENTATION_VERSION));
    push_text(&mut bytes, "stream-semantics-version");
    push_u64(&mut bytes, u64::from(STREAM_SEMANTICS_VERSION));
    push_text(&mut bytes, "stream-position-identity-domain");
    push_text(&mut bytes, STREAM_POSITION_IDENTITY_DOMAIN);
    bytes
}

fn green_inputs(fixture: &FixtureBits) -> Vec<u8> {
    let mut bytes = common_frame_prefix(b"bedrock:fs-rand-to-fs-la-gemm-replay:v1");
    push_text(&mut bytes, "logical-coordinate-policy");
    push_text(
        &mut bytes,
        "matrix-kernel;tile=row-u32;index=column-u64;Stream::at;word=(block[1]<<32)|block[0]",
    );
    push_text(&mut bytes, "finite-f64-policy");
    push_text(
        &mut bytes,
        "top53=word>>11;value=(f64(top53)-2^52)*2^-50;exact-dyadic:[-4,4):v1",
    );
    push_text(&mut bytes, "root-seed");
    push_u64(&mut bytes, ROOT_SEED);
    push_text(&mut bytes, "matrix-kernels-a-b-c");
    for kernel in [A_KERNEL, B_KERNEL, C_KERNEL] {
        push_u64(&mut bytes, u64::from(kernel));
    }
    push_text(&mut bytes, "row-major-dimensions-m-n-k");
    for dimension in [M, N, K] {
        push_len(&mut bytes, dimension);
    }
    push_text(&mut bytes, "alpha-beta-bits");
    push_u64(&mut bytes, ALPHA.to_bits());
    push_u64(&mut bytes, BETA.to_bits());
    push_text(&mut bytes, "generation-plans");
    push_len(&mut bytes, GENERATION_PLANS.len());
    for plan in GENERATION_PLANS {
        push_len(&mut bytes, plan.partitions);
        push_u64(&mut bytes, u64::from(u8::from(plan.reverse_partitions)));
    }
    push_text(&mut bytes, "parallel-worker-requests");
    push_len(&mut bytes, PARALLEL_WORKERS.len());
    for workers in PARALLEL_WORKERS {
        push_len(&mut bytes, workers);
    }
    push_text(&mut bytes, "a-input-bits");
    push_bits(&mut bytes, &fixture.a);
    push_text(&mut bytes, "b-input-bits");
    push_bits(&mut bytes, &fixture.b);
    push_text(&mut bytes, "c-input-bits");
    push_bits(&mut bytes, &fixture.c);
    bytes
}

fn corruption_coordinates() -> (usize, u32) {
    let output_count = u64::try_from(M * N).expect("fixture output count fits u64");
    let output = usize::try_from(RED_SEED % output_count).expect("derived output index fits usize");
    let bit = u32::try_from((RED_SEED >> 16) % 52).expect("derived mantissa bit fits u32");
    (output, bit)
}

fn red_inputs(green: &[u8]) -> Vec<u8> {
    let (output, bit) = corruption_coordinates();
    let mut bytes = common_frame_prefix(b"bedrock:fs-rand-to-fs-la-gemm-red:v1");
    push_nested(&mut bytes, "nested-green-input-frame", green);
    push_text(&mut bytes, "corruption-seed-output-mantissa-bit");
    push_u64(&mut bytes, RED_SEED);
    push_len(&mut bytes, output);
    push_u64(&mut bytes, u64::from(bit));
    push_text(
        &mut bytes,
        "policy=flip-one-derived-mantissa-bit-in-serial-reference:v1",
    );
    bytes
}

fn panic_message(payload: &(dyn core::any::Any + Send)) -> String {
    payload
        .downcast_ref::<String>()
        .cloned()
        .or_else(|| {
            payload
                .downcast_ref::<&str>()
                .map(|text| (*text).to_owned())
        })
        .unwrap_or_else(|| "non-text panic payload".to_owned())
}

fn f64_values(bits: &[u64]) -> Vec<f64> {
    bits.iter().map(|&value| f64::from_bits(value)).collect()
}

fn output_bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|value| value.to_bits()).collect()
}

fn execute(fixture: &FixtureBits) -> Result<Evaluation, String> {
    let a = f64_values(&fixture.a);
    let b = f64_values(&fixture.b);
    let c = f64_values(&fixture.c);
    let mut serial_values = c.clone();
    catch_unwind(AssertUnwindSafe(|| {
        gemm_f64(M, N, K, ALPHA, &a, &b, BETA, &mut serial_values);
    }))
    .map_err(|payload| format!("stage=serial-gemm; panic={}", panic_message(&*payload)))?;
    let serial = output_bits(&serial_values);

    let mut parallel = Vec::with_capacity(PARALLEL_WORKERS.len());
    for workers in PARALLEL_WORKERS {
        let mut values = c.clone();
        catch_unwind(AssertUnwindSafe(|| {
            gemm_f64_parallel(M, N, K, ALPHA, &a, &b, BETA, &mut values, workers);
        }))
        .map_err(|payload| {
            format!(
                "stage=parallel-gemm; workers={workers}; panic={}",
                panic_message(&*payload)
            )
        })?;
        let bits = output_bits(&values);
        if let Some((index, expected, computed)) = first_mismatch(&serial, &bits) {
            return Err(format!(
                "stage=parallel-bit-identity; workers={workers}; index={index}; serial_bits=0x{expected:016x}; parallel_bits=0x{computed:016x}"
            ));
        }
        parallel.push((workers, bits));
    }

    Ok(Evaluation {
        fixture_digest: digest_bits(&[
            fixture.a.as_slice(),
            fixture.b.as_slice(),
            fixture.c.as_slice(),
        ]),
        serial,
        parallel,
    })
}

fn evaluate() -> Result<Evaluation, String> {
    execute(&validate_generation()?)
}

fn green_outcome(input_frame: &[u8]) -> CaseOutcome {
    let inputs_hex = hex_bytes(input_frame);
    let first = match evaluate() {
        Ok(evidence) => evidence,
        Err(error) => {
            return CaseOutcome::fail(format!("{error}; inputs_hex={inputs_hex}"))
                .with_evidence("crates/fs-la/CONTRACT.md#determinism-class")
                .with_evidence("crates/fs-rand/CONTRACT.md#determinism-class");
        }
    };
    let replay = match evaluate() {
        Ok(evidence) => evidence,
        Err(error) => {
            return CaseOutcome::fail(format!(
                "stage=same-run-replay; replay_error={error}; inputs_hex={inputs_hex}"
            ))
            .with_evidence("crates/fs-la/CONTRACT.md#determinism-class");
        }
    };
    if first != replay {
        return CaseOutcome::fail(format!(
            "stage=same-run-replay; first={first:016x?}; replay={replay:016x?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-la/CONTRACT.md#determinism-class");
    }
    CaseOutcome::pass(format!(
        "shape={M}x{N}x{K}; generation_plans={}; worker_requests=1,2,3,5,8; threaded_worker_requests=2,3,5,8; fixture_digest={:016x}; output_digest={:016x}; execution_tier={}; same_run=identical",
        GENERATION_PLANS.len(),
        first.fixture_digest,
        digest_bits(&[first.serial.as_slice()]),
        gemm_execution_tier(),
    ))
    .with_evidence("crates/fs-la/CONTRACT.md#determinism-class")
    .with_evidence("crates/fs-rand/CONTRACT.md#determinism-class")
}

fn red_outcome(input_frame: &[u8]) -> CaseOutcome {
    let inputs_hex = hex_bytes(input_frame);
    let first = match evaluate() {
        Ok(evidence) => evidence,
        Err(error) => {
            return CaseOutcome::fail(format!(
                "stage=red-prerequisite; error={error}; inputs_hex={inputs_hex}"
            ));
        }
    };
    let replay = match evaluate() {
        Ok(evidence) => evidence,
        Err(error) => {
            return CaseOutcome::fail(format!(
                "stage=red-replay-prerequisite; error={error}; inputs_hex={inputs_hex}"
            ));
        }
    };
    if first != replay {
        return CaseOutcome::fail(format!(
            "stage=red-same-run-replay; first={first:016x?}; replay={replay:016x?}; inputs_hex={inputs_hex}"
        ));
    }

    let (output, bit) = corruption_coordinates();
    let canonical = first.serial[output];
    let corrupted = canonical ^ (1_u64 << bit);
    if canonical == corrupted {
        return CaseOutcome::pass("derived corruption did not move the reference");
    }
    CaseOutcome::fail(format!(
        "stage=seeded-output-reference-corruption; seed=0x{RED_SEED:016x}; output={output}; bit={bit}; actual_bits=0x{canonical:016x}; canonical_bits=0x{canonical:016x}; corrupted_bits=0x{corrupted:016x}; inputs_hex={inputs_hex}"
    ))
    .with_evidence("crates/fs-la/tests/rand_gemm_replay_casebook.rs#seeded-corruption")
}

#[test]
fn philox_to_parallel_gemm_casebook_emits_replay_complete_green_record() {
    assert_eq!(CASEBOOK_RECORD_VERSION, 1);
    assert_eq!(GEMM_BIT_SEMANTICS_VERSION, 1);
    assert_eq!(STREAM_SEMANTICS_VERSION, 1);
    assert_eq!(GEMM_IMPLEMENTATION_VERSION, 4);
    let fixture = materialize_fixture(GENERATION_PLANS[0]);
    let inputs = green_inputs(&fixture);
    let inputs_digest = fnv1a64(&inputs);
    assert_eq!(
        (inputs.len(), inputs_digest),
        (GREEN_FRAME_LEN, GREEN_FRAME_DIGEST)
    );

    let report = Suite::new(SUITE)
        .case(
            "philox-logical-input-parallel-gemm-bit-replay",
            inputs_digest,
            ToleranceSpec::Exact,
            move || green_outcome(&inputs),
        )
        .run();
    report.assert_green();
    let [record] = report.records.as_slice() else {
        panic!("the composition suite must emit exactly one record");
    };
    assert_eq!(record.case, "philox-logical-input-parallel-gemm-bit-replay");
    assert!(record.details.contains("threaded_worker_requests=2,3,5,8"));
    assert!(record.details.contains("same_run=identical"));
}

#[test]
fn disclosed_seeded_output_reference_corruption_turns_suite_red() {
    let fixture = materialize_fixture(GENERATION_PLANS[0]);
    let green = green_inputs(&fixture);
    let inputs = red_inputs(&green);
    let inputs_digest = fnv1a64(&inputs);
    assert_eq!(
        (inputs.len(), inputs_digest),
        (RED_FRAME_LEN, RED_FRAME_DIGEST)
    );
    let (output, bit) = corruption_coordinates();

    let make_report = || {
        let input_frame = inputs.clone();
        Suite::new(SUITE)
            .case(
                "seeded-parallel-gemm-output-reference-corruption",
                inputs_digest,
                ToleranceSpec::Exact,
                move || red_outcome(&input_frame),
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
    assert!(
        first_failure
            .details
            .contains("stage=seeded-output-reference-corruption")
    );
    assert!(
        first_failure
            .details
            .contains(&format!("seed=0x{RED_SEED:016x}"))
    );
    assert!(first_failure.details.contains(&format!("output={output}")));
    assert!(first_failure.details.contains(&format!("bit={bit}")));
    assert!(first_failure.details.contains("inputs_hex="));
    assert!(
        first_failure
            .json_line()
            .contains("\"tolerance\":\"exact\",\"pass\":false")
    );
    let panic = catch_unwind(|| first.assert_green())
        .expect_err("the Casebook merge gate must reject the disclosed corruption");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("Casebook panic carries text");
    assert!(message.contains("seeded-parallel-gemm-output-reference-corruption"));
    assert!(message.contains(&format!("seed=0x{RED_SEED:016x}")));
}
