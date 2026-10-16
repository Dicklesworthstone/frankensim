//! Structured BEDROCK conformance records for fs-rand.
//!
//! These cases do not replace the crate's larger statistical and replay
//! batteries. They expose the load-bearing deterministic laws through the
//! shared fs-casebook schema so a failed central run is reproducible from its
//! JSON line alone (Gauntlet G0/G5 instrumentation).

use fs_casebook::{CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, ToleranceSpec, fnv1a64};
use fs_rand::philox::philox4x32_10;
use fs_rand::{STREAM_CHECKPOINT_VERSION, STREAM_SEMANTICS_VERSION, Stream, StreamKey};

const SUITE: &str = "fs-rand/bedrock-conformance-v1";
const STREAM_KEY: StreamKey = StreamKey {
    seed: 0x5EED_0001_DEAD_BEEF,
    kernel: 7,
    tile: 42,
};
const KATS: [([u32; 4], [u32; 2], [u32; 4]); 3] = [
    (
        [0, 0, 0, 0],
        [0, 0],
        [0x6627_e8d5, 0xe169_c58d, 0xbc57_ac4c, 0x9b00_dbd8],
    ),
    (
        [u32::MAX; 4],
        [u32::MAX; 2],
        [0x408f_276d, 0x41c8_3b0e, 0xa20b_c7c6, 0x6d54_51fd],
    ),
    (
        [0x243f_6a88, 0x85a3_08d3, 0x1319_8a2e, 0x0370_7344],
        [0xa409_3822, 0x299f_31d0],
        [0xd16c_fe09, 0x94fd_cceb, 0x5001_e420, 0x2412_6ea1],
    ),
];

fn push_key(bytes: &mut Vec<u8>, key: StreamKey) {
    bytes.extend_from_slice(&key.seed.to_le_bytes());
    bytes.extend_from_slice(&key.kernel.to_le_bytes());
    bytes.extend_from_slice(&key.tile.to_le_bytes());
}

fn kat_inputs() -> Vec<u8> {
    let mut bytes = b"fs-rand:philox4x32-10:kats:v1".to_vec();
    for (counter, key, reference) in KATS {
        for word in counter {
            bytes.extend_from_slice(&word.to_le_bytes());
        }
        for word in key {
            bytes.extend_from_slice(&word.to_le_bytes());
        }
        for word in reference {
            bytes.extend_from_slice(&word.to_le_bytes());
        }
    }
    bytes
}

fn stream_inputs(case: &[u8], draws: u64) -> Vec<u8> {
    let mut bytes = b"fs-rand:stream-conformance:v1".to_vec();
    let case_len = u64::try_from(case.len()).expect("case labels fit the u64 input frame");
    bytes.extend_from_slice(&case_len.to_le_bytes());
    bytes.extend_from_slice(case);
    bytes.extend_from_slice(&STREAM_SEMANTICS_VERSION.to_le_bytes());
    bytes.extend_from_slice(&STREAM_CHECKPOINT_VERSION.to_le_bytes());
    push_key(&mut bytes, STREAM_KEY);
    bytes.extend_from_slice(&draws.to_le_bytes());
    bytes
}

fn philox_kat_outcome() -> CaseOutcome {
    for (vector, (counter, key, reference)) in KATS.into_iter().enumerate() {
        let computed = philox4x32_10(counter, key);
        if computed != reference {
            return CaseOutcome::fail(format!(
                "vector={vector}; counter={counter:08x?}; key={key:08x?}; computed={computed:08x?}; reference={reference:08x?}"
            ))
            .with_evidence("Random123/philox4x32-10/kat_vectors");
        }
    }
    CaseOutcome::pass("vectors=3; mismatches=0")
        .with_evidence("Random123/philox4x32-10/kat_vectors")
}

fn random_access_outcome() -> CaseOutcome {
    const DRAWS: u64 = 64;
    let mut stream = STREAM_KEY.stream();
    for index in 0..DRAWS {
        let sequential = stream.next_u64();
        let block = Stream::at(STREAM_KEY, index);
        let random_access = (u64::from(block[1]) << 32) | u64::from(block[0]);
        if random_access != sequential {
            return CaseOutcome::fail(format!(
                "stream_semantics_version={STREAM_SEMANTICS_VERSION}; checkpoint_version={STREAM_CHECKPOINT_VERSION}; seed=0x5eed0001deadbeef; kernel=7; tile=42; draws=64; index={index}; random_access=0x{random_access:016x}; sequential=0x{sequential:016x}"
            ))
            .with_evidence("crates/fs-rand/CONTRACT.md#invariants");
        }
    }
    CaseOutcome::pass("seed=0x5eed0001deadbeef; kernel=7; tile=42; draws=64; mismatches=0")
        .with_evidence("crates/fs-rand/CONTRACT.md#invariants")
}

fn checkpoint_resume_outcome() -> CaseOutcome {
    const WARMUP: usize = 32;
    const TAIL: usize = 16;
    let mut stream = STREAM_KEY.stream();
    for _ in 0..WARMUP {
        let _ = stream.next_u64();
    }
    let checkpoint = stream.checkpoint();
    let retained = checkpoint.to_canonical_le_bytes();
    let reference: Vec<u64> = (0..TAIL).map(|_| stream.next_u64()).collect();
    let mut resumed = match Stream::resume_retained(&retained) {
        Ok(resumed) => resumed,
        Err(error) => {
            return CaseOutcome::fail(format!(
                "checkpoint_version={}; stream_semantics_version={}; seed=0x5eed0001deadbeef; kernel=7; tile=42; checkpoint_index={}; warmup=32; replay_tail=16; checkpoint_frame={retained:02x?}; resume_refusal={error}",
                checkpoint.checkpoint_version,
                checkpoint.stream_semantics_version,
                checkpoint.index,
            ))
            .with_evidence("crates/fs-rand/CONTRACT.md#public-types-and-semantics");
        }
    };
    for (offset, expected) in reference.into_iter().enumerate() {
        let computed = resumed.next_u64();
        if computed != expected {
            return CaseOutcome::fail(format!(
                "checkpoint_version={}; stream_semantics_version={}; seed=0x5eed0001deadbeef; kernel=7; tile=42; checkpoint_index={}; warmup=32; replay_tail=16; checkpoint_frame={retained:02x?}; tail_offset={offset}; computed=0x{computed:016x}; reference=0x{expected:016x}",
                checkpoint.checkpoint_version,
                checkpoint.stream_semantics_version,
                checkpoint.index,
            ))
            .with_evidence("crates/fs-rand/CONTRACT.md#public-types-and-semantics");
        }
    }
    CaseOutcome::pass(
        "seed=0x5eed0001deadbeef; kernel=7; tile=42; warmup=32; replay_tail=16; mismatches=0",
    )
    .with_evidence("crates/fs-rand/CONTRACT.md#public-types-and-semantics")
}

#[test]
fn bedrock_casebook_suite_emits_replay_complete_green_records() {
    let kat_digest = fnv1a64(&kat_inputs());
    let random_access_digest = fnv1a64(&stream_inputs(b"random-access", 64));
    let checkpoint_digest = fnv1a64(&stream_inputs(b"checkpoint-resume:32+16", 48));
    assert_eq!(kat_digest, 0x947b_c3da_1b13_24ca);
    assert_eq!(random_access_digest, 0x2048_001f_faca_223b);
    assert_eq!(checkpoint_digest, 0x1d9e_2e49_ccc0_67d0);
    let report = Suite::new(SUITE)
        .case(
            "philox-random123-kat",
            kat_digest,
            ToleranceSpec::Exact,
            philox_kat_outcome,
        )
        .case(
            "stream-random-access-sequential",
            random_access_digest,
            ToleranceSpec::Exact,
            random_access_outcome,
        )
        .case(
            "checkpoint-resume-tail",
            checkpoint_digest,
            ToleranceSpec::Exact,
            checkpoint_resume_outcome,
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
            "philox-random123-kat",
            "stream-random-access-sequential",
            "checkpoint-resume-tail",
        ]
    );
    assert_eq!(
        report.records[0].json_line(),
        format!(
            concat!(
                "{{\"casebook\":{},\"suite\":\"fs-rand/bedrock-conformance-v1\",",
                "\"case\":\"philox-random123-kat\",\"inputs_digest\":\"947bc3da1b1324ca\",",
                "\"tolerance\":\"exact\",\"pass\":true,",
                "\"details\":\"vectors=3; mismatches=0\",",
                "\"evidence\":[\"Random123/philox4x32-10/kat_vectors\"]}}"
            ),
            CASEBOOK_RECORD_VERSION,
        ),
        "the structured record schema and field order are contract"
    );
}

#[test]
fn disclosed_seeded_corruption_turns_the_casebook_suite_red() {
    const CORRUPTION_SEED: u64 = 0xF5AA_0001;
    let (counter, key, reference) = KATS[0];
    let word = (CORRUPTION_SEED as usize) % reference.len();
    let bit = ((CORRUPTION_SEED.rotate_left(17) ^ CORRUPTION_SEED.rotate_right(7)) % 32) as u32;
    assert_eq!(word, 1);
    assert_eq!(bit, 0);
    let mut corrupted_reference = reference;
    corrupted_reference[word] ^= 1_u32 << bit;

    let mut inputs = b"fs-rand:seeded-kat-corruption:v1".to_vec();
    inputs.extend_from_slice(&CORRUPTION_SEED.to_le_bytes());
    for value in counter.into_iter().chain(key) {
        inputs.extend_from_slice(&value.to_le_bytes());
    }
    for value in corrupted_reference {
        inputs.extend_from_slice(&value.to_le_bytes());
    }
    let inputs_digest = fnv1a64(&inputs);
    assert_eq!(inputs_digest, 0x5559_482e_52cc_f742);
    let report = Suite::new(SUITE)
        .case(
            "seeded-philox-reference-corruption",
            inputs_digest,
            ToleranceSpec::Exact,
            move || {
                let computed = philox4x32_10(counter, key);
                if computed == corrupted_reference {
                    CaseOutcome::pass("seeded corruption was not detected")
                } else {
                    CaseOutcome::fail(format!(
                        "seed=0x{CORRUPTION_SEED:016x}; vector=0; word={word}; bit={bit}; counter={counter:08x?}; key={key:08x?}; computed={computed:08x?}; corrupted_reference={corrupted_reference:08x?}"
                    ))
                    .with_evidence("crates/fs-rand/tests/conformance.rs#seeded-corruption")
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
    assert_eq!(failure.case, "seeded-philox-reference-corruption");
    assert_eq!(failure.inputs_digest, "5559482e52ccf742");
    assert!(
        failure
            .details
            .contains(&format!("seed=0x{CORRUPTION_SEED:016x}"))
    );
    assert!(failure.details.contains(&format!("word={word}; bit={bit}")));
    assert!(failure.details.contains("vector=0"));
    assert!(failure.details.contains(&format!("counter={counter:08x?}")));
    assert!(failure.details.contains(&format!("key={key:08x?}")));
    let line = failure.json_line();
    assert!(line.contains("\"tolerance\":\"exact\",\"pass\":false"));
    assert!(line.contains("computed=["));
    assert!(line.contains("corrupted_reference=["));

    let panic = std::panic::catch_unwind(|| report.assert_green())
        .expect_err("the merge-gate assertion must reject the seeded failure");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("casebook panic carries text");
    assert!(message.contains("seeded-philox-reference-corruption"));
    assert!(message.contains(&format!("seed=0x{CORRUPTION_SEED:016x}")));
}
