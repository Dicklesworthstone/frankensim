//! Structured BEDROCK conformance records for fs-fft.
//!
//! The larger FFT battery retains broad oracle, theorem, cancellation, and
//! performance coverage. These cases expose the load-bearing bit contracts
//! through fs-casebook so a central failure carries a canonical reproducer.

use core::fmt::Write as _;
use fs_casebook::{CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, ToleranceSpec, fnv1a64};
use fs_exec::{CancelGate, PoolConfig, TilePool};
use fs_fft::{C64, Fft, FftNd, TRANSFORM_BIT_SEMANTICS_VERSION};
use fs_substrate::affinity::CcdTopology;

const SUITE: &str = "fs-fft/bedrock-conformance-v1";
const LCG_MULTIPLIER: u64 = 6_364_136_223_846_793_005;
const LCG_INCREMENT: u64 = 1_442_695_040_888_963_407;
const STAGE_SEED: u64 = 0xD15C;
const STAGE_GOLDEN: u64 = 0x22dd_b617_266e_a792;
const ND_SEED: u64 = 0xD1D5;
const POOL_SEED: u64 = 0xFD1D;
const POOL_TOPOLOGY: CcdTopology = CcdTopology {
    ccds: 2,
    cores_per_ccd: 8,
};
const WORKERS: [usize; 4] = [1, 2, 3, 7];
const ND_SHAPES: [&[usize]; 5] = [&[8, 16], &[4, 8, 2], &[2, 2, 2, 4], &[1, 16, 1, 4], &[32]];

const KAT_INPUT: [[u64; 2]; 2] = [
    [0x3ff0_0000_0000_0000, 0x4000_0000_0000_0000],
    [0x4008_0000_0000_0000, 0x4010_0000_0000_0000],
];
const KAT_FORWARD: [[u64; 2]; 2] = [
    [0x4010_0000_0000_0000, 0x4018_0000_0000_0000],
    [0xc000_0000_0000_0000, 0xc000_0000_0000_0000],
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

fn hex_bytes(bytes: &[u8]) -> String {
    let mut encoded = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        write!(&mut encoded, "{byte:02x}").expect("writing to String cannot fail");
    }
    encoded
}

fn push_complex_vectors(bytes: &mut Vec<u8>, vectors: &[[u64; 2]]) {
    push_len(bytes, vectors.len());
    for [re, im] in vectors {
        push_u64(bytes, *re);
        push_u64(bytes, *im);
    }
}

fn kat_inputs() -> Vec<u8> {
    let mut bytes = b"fs-fft:n2-complex-kat-roundtrip:v1".to_vec();
    push_u64(&mut bytes, u64::from(TRANSFORM_BIT_SEMANTICS_VERSION));
    push_u64(&mut bytes, 2);
    push_complex_vectors(&mut bytes, &KAT_INPUT);
    push_complex_vectors(&mut bytes, &KAT_FORWARD);
    push_complex_vectors(&mut bytes, &KAT_INPUT);
    bytes
}

fn stage_inputs() -> Vec<u8> {
    let mut bytes = b"fs-fft:stage-golden:v1".to_vec();
    for value in [
        u64::from(TRANSFORM_BIT_SEMANTICS_VERSION),
        128,
        16,
        STAGE_SEED,
        LCG_MULTIPLIER,
        LCG_INCREMENT,
        11,
        1_u64 << 53,
        0.5_f64.to_bits(),
        STAGE_GOLDEN,
    ] {
        push_u64(&mut bytes, value);
    }
    bytes
}

fn pooled_inputs() -> Vec<u8> {
    let mut bytes = b"fs-fft:nd-pooled-serial-bitwise:v1".to_vec();
    for value in [
        u64::from(TRANSFORM_BIT_SEMANTICS_VERSION),
        ND_SEED,
        LCG_MULTIPLIER,
        LCG_INCREMENT,
        11,
        1_u64 << 53,
        0.5_f64.to_bits(),
    ] {
        push_u64(&mut bytes, value);
    }
    push_len(&mut bytes, ND_SHAPES.len());
    for shape in ND_SHAPES {
        push_len(&mut bytes, shape.len());
        for dimension in shape {
            push_u64(
                &mut bytes,
                u64::try_from(*dimension).expect("FFT fixture dimensions fit u64"),
            );
        }
    }
    push_len(&mut bytes, WORKERS.len());
    for workers in WORKERS {
        push_u64(
            &mut bytes,
            u64::try_from(workers).expect("worker counts fit u64"),
        );
    }
    push_u64(&mut bytes, POOL_SEED);
    push_u64(&mut bytes, u64::from(POOL_TOPOLOGY.ccds));
    push_u64(&mut bytes, u64::from(POOL_TOPOLOGY.cores_per_ccd));
    push_u64(&mut bytes, 2);
    push_u64(&mut bytes, 0); // forward
    push_u64(&mut bytes, 1); // inverse
    bytes
}

fn seeded_corruption_inputs(seed: u64, bit: u32, corrupted: u64) -> Vec<u8> {
    let stage = stage_inputs();
    let mut bytes = b"fs-fft:seeded-stage-golden-corruption:v1".to_vec();
    push_u64(&mut bytes, seed);
    push_u64(&mut bytes, u64::from(bit));
    push_len(&mut bytes, stage.len());
    bytes.extend_from_slice(&stage);
    push_u64(&mut bytes, STAGE_GOLDEN);
    push_u64(&mut bytes, corrupted);
    bytes
}

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(LCG_MULTIPLIER)
        .wrapping_add(LCG_INCREMENT);
    ((*seed >> 11) as f64) / ((1_u64 << 53) as f64) - 0.5
}

fn complex_from_bits([re, im]: [u64; 2]) -> C64 {
    C64::new(f64::from_bits(re), f64::from_bits(im))
}

fn complex_bits(value: C64) -> [u64; 2] {
    [value.re.to_bits(), value.im.to_bits()]
}

fn first_mismatch(actual: &[C64], reference: &[C64]) -> Option<(usize, [u64; 2], [u64; 2])> {
    actual
        .iter()
        .zip(reference)
        .enumerate()
        .find_map(|(index, (actual, reference))| {
            let actual = complex_bits(*actual);
            let reference = complex_bits(*reference);
            (actual != reference).then_some((index, actual, reference))
        })
}

fn first_expected_mismatch(
    actual: &[C64],
    reference: &[[u64; 2]],
) -> Option<(usize, [u64; 2], [u64; 2])> {
    actual
        .iter()
        .zip(reference)
        .enumerate()
        .find_map(|(index, (actual, reference))| {
            let actual = complex_bits(*actual);
            (actual != *reference).then_some((index, actual, *reference))
        })
}

fn fixture(total: usize) -> Vec<C64> {
    let mut seed = ND_SEED;
    (0..total)
        .map(|_| C64::new(lcg(&mut seed), lcg(&mut seed)))
        .collect()
}

fn stage_fingerprint() -> u64 {
    let plan = Fft::new(128);
    let mut scratch = vec![C64::default(); 128];
    let mut seed = STAGE_SEED;
    let mut outputs = Vec::new();
    for _ in 0..16 {
        let mut data: Vec<C64> = (0..128)
            .map(|_| C64::new(lcg(&mut seed), lcg(&mut seed)))
            .collect();
        plan.forward(&mut data, &mut scratch);
        for value in data {
            outputs.extend_from_slice(&value.re.to_bits().to_le_bytes());
            outputs.extend_from_slice(&value.im.to_bits().to_le_bytes());
        }
    }
    fnv1a64(&outputs)
}

fn kat_outcome() -> CaseOutcome {
    let inputs_hex = hex_bytes(&kat_inputs());
    let plan = Fft::new(2);
    let mut scratch = [C64::default(); 2];
    let mut data = KAT_INPUT.map(complex_from_bits);
    plan.forward(&mut data, &mut scratch);
    if let Some((index, actual, reference)) = first_expected_mismatch(&data, &KAT_FORWARD) {
        return CaseOutcome::fail(format!(
            "direction=forward; n=2; index={index}; actual_bits={actual:016x?}; reference_bits={reference:016x?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-fft/CONTRACT.md#invariants");
    }
    plan.inverse(&mut data, &mut scratch);
    if let Some((index, actual, reference)) = first_expected_mismatch(&data, &KAT_INPUT) {
        return CaseOutcome::fail(format!(
            "direction=inverse-roundtrip; n=2; index={index}; actual_bits={actual:016x?}; reference_bits={reference:016x?}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-fft/CONTRACT.md#invariants");
    }
    CaseOutcome::pass("n=2; forward_exact=true; inverse_roundtrip_exact=true")
        .with_evidence("crates/fs-fft/CONTRACT.md#invariants")
}

fn stage_outcome() -> CaseOutcome {
    let inputs_hex = hex_bytes(&stage_inputs());
    let computed = stage_fingerprint();
    if computed != STAGE_GOLDEN {
        return CaseOutcome::fail(format!(
            "n=128; batches=16; seed=0x{STAGE_SEED:016x}; computed=0x{computed:016x}; reference=0x{STAGE_GOLDEN:016x}; inputs_hex={inputs_hex}"
        ))
        .with_evidence("crates/fs-fft/CONTRACT.md#determinism-class");
    }
    CaseOutcome::pass(format!(
        "n=128; batches=16; seed=0x{STAGE_SEED:016x}; fingerprint=0x{computed:016x}"
    ))
    .with_evidence("crates/fs-fft/CONTRACT.md#determinism-class")
}

fn pooled_outcome() -> CaseOutcome {
    let inputs_hex = hex_bytes(&pooled_inputs());
    for shape in ND_SHAPES {
        let plan = FftNd::new(shape);
        let source = fixture(plan.total());
        let mut forward_reference = source.clone();
        plan.forward(&mut forward_reference);
        let mut inverse_reference = forward_reference.clone();
        plan.inverse(&mut inverse_reference);
        for workers in WORKERS {
            let pool = TilePool::new(PoolConfig::new(workers, POOL_TOPOLOGY, POOL_SEED));
            let gate = CancelGate::new();
            let mut forward = source.clone();
            if let Err(error) = plan.forward_pooled(&mut forward, &pool, &gate) {
                return CaseOutcome::fail(format!(
                    "direction=forward; shape={shape:?}; workers={workers}; run_error={error}; inputs_hex={inputs_hex}"
                ))
                .with_evidence("crates/fs-fft/CONTRACT.md#invariants");
            }
            if let Some((index, actual, reference)) = first_mismatch(&forward, &forward_reference) {
                return CaseOutcome::fail(format!(
                    "direction=forward; shape={shape:?}; workers={workers}; index={index}; actual_bits={actual:016x?}; reference_bits={reference:016x?}; inputs_hex={inputs_hex}"
                ))
                .with_evidence("crates/fs-fft/CONTRACT.md#invariants");
            }
            let mut inverse = forward;
            if let Err(error) = plan.inverse_pooled(&mut inverse, &pool, &gate) {
                return CaseOutcome::fail(format!(
                    "direction=inverse; shape={shape:?}; workers={workers}; run_error={error}; inputs_hex={inputs_hex}"
                ))
                .with_evidence("crates/fs-fft/CONTRACT.md#invariants");
            }
            if let Some((index, actual, reference)) = first_mismatch(&inverse, &inverse_reference) {
                return CaseOutcome::fail(format!(
                    "direction=inverse; shape={shape:?}; workers={workers}; index={index}; actual_bits={actual:016x?}; reference_bits={reference:016x?}; inputs_hex={inputs_hex}"
                ))
                .with_evidence("crates/fs-fft/CONTRACT.md#invariants");
            }
        }
    }
    CaseOutcome::pass("shapes=5; worker_counts=4; directions=2; serial_pooled_bit_mismatches=0")
        .with_evidence("crates/fs-fft/CONTRACT.md#invariants")
}

#[test]
fn bedrock_casebook_suite_emits_replay_complete_green_records() {
    assert_eq!(CASEBOOK_RECORD_VERSION, 1);
    let kat_digest = fnv1a64(&kat_inputs());
    let stage_digest = fnv1a64(&stage_inputs());
    let pooled_digest = fnv1a64(&pooled_inputs());
    assert_eq!(kat_digest, 0x2d68_1fb6_9ef4_2581);
    assert_eq!(stage_digest, 0xe02b_6506_d652_a5ae);
    assert_eq!(pooled_digest, 0x36d8_d637_d73d_2496);

    let report = Suite::new(SUITE)
        .case(
            "n2-complex-kat-roundtrip",
            kat_digest,
            ToleranceSpec::Exact,
            kat_outcome,
        )
        .case(
            "stage-cross-isa-golden",
            stage_digest,
            ToleranceSpec::Exact,
            stage_outcome,
        )
        .case(
            "nd-pooled-serial-bitwise",
            pooled_digest,
            ToleranceSpec::Exact,
            pooled_outcome,
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
            "n2-complex-kat-roundtrip",
            "stage-cross-isa-golden",
            "nd-pooled-serial-bitwise",
        ]
    );
    assert_eq!(
        report.records[0].json_line(),
        format!(
            concat!(
                "{{\"casebook\":{},\"suite\":\"fs-fft/bedrock-conformance-v1\",",
                "\"case\":\"n2-complex-kat-roundtrip\",",
                "\"inputs_digest\":\"2d681fb69ef42581\",\"tolerance\":\"exact\",",
                "\"pass\":true,\"details\":\"n=2; forward_exact=true; inverse_roundtrip_exact=true\",",
                "\"evidence\":[\"crates/fs-fft/CONTRACT.md#invariants\"]}}"
            ),
            CASEBOOK_RECORD_VERSION,
        ),
        "the structured record schema and field order are contract"
    );
}

#[test]
fn disclosed_seeded_corruption_turns_the_casebook_suite_red() {
    const CORRUPTION_SEED: u64 = 0xF5FF_0001;
    let bit = CORRUPTION_SEED.trailing_zeros();
    assert_eq!(bit, 0);
    let corrupted = STAGE_GOLDEN ^ (1_u64 << bit);
    assert_eq!(corrupted, 0x22dd_b617_266e_a793);
    let inputs = seeded_corruption_inputs(CORRUPTION_SEED, bit, corrupted);
    let inputs_hex = hex_bytes(&inputs);
    let inputs_digest = fnv1a64(&inputs);
    assert_eq!(inputs_digest, 0x428c_5b33_ad58_f30e);

    let report = Suite::new(SUITE)
        .case(
            "seeded-stage-golden-corruption",
            inputs_digest,
            ToleranceSpec::Exact,
            move || {
                let computed = stage_fingerprint();
                if computed == corrupted {
                    CaseOutcome::pass("seeded stage corruption escaped detection")
                } else {
                    CaseOutcome::fail(format!(
                        "seed=0x{CORRUPTION_SEED:016x}; bit={bit}; computed=0x{computed:016x}; canonical=0x{STAGE_GOLDEN:016x}; corrupted=0x{corrupted:016x}; inputs_hex={inputs_hex}"
                    ))
                    .with_evidence(
                        "crates/fs-fft/tests/conformance.rs#disclosed-seeded-corruption",
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
    assert_eq!(failure.case, "seeded-stage-golden-corruption");
    assert_eq!(failure.inputs_digest, "428c5b33ad58f30e");
    assert!(
        failure
            .details
            .contains(&format!("seed=0x{CORRUPTION_SEED:016x}"))
    );
    assert!(failure.details.contains("bit=0"));
    assert!(failure.details.contains("computed=0x22ddb617266ea792"));
    assert!(failure.details.contains("canonical=0x22ddb617266ea792"));
    assert!(failure.details.contains("corrupted=0x22ddb617266ea793"));
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
    assert!(message.contains("seeded-stage-golden-corruption"));
    assert!(message.contains(&format!("seed=0x{CORRUPTION_SEED:016x}")));
}
