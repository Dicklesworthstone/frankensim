//! Structured BEDROCK conformance records for fs-la.
//!
//! The crate's larger numerical batteries retain broad shape, factorization,
//! cancellation, property, and performance coverage. These cases expose three
//! load-bearing deterministic policies through fs-casebook so a central run can
//! diagnose and replay failures from structured records (Gauntlet G0).

use fs_casebook::{CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, ToleranceSpec, fnv1a64};
use fs_la::factor::{FACTOR_BIT_SEMANTICS_VERSION, FactorError, cholesky, lu};
use fs_la::gemm::GEMM_BIT_SEMANTICS_VERSION;
use fs_la::{gemm_execution_tier, gemm_f64};

const SUITE: &str = "fs-la/bedrock-conformance-v1";
const LCG_MULTIPLIER: u64 = 6_364_136_223_846_793_005;
const LCG_INCREMENT: u64 = 1_442_695_040_888_963_407;
const GEMM_GOLDEN: u64 = 0x1d7a_a3c6_b631_7ef0;
const GOLDEN_M: usize = 48;
const GOLDEN_N: usize = 36;
const GOLDEN_K: usize = 300;
const GOLDEN_A_SEED: u64 = 0x60;
const GOLDEN_B_SEED: u64 = 0x61;

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

fn lcg(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(LCG_MULTIPLIER)
        .wrapping_add(LCG_INCREMENT);
    ((*seed >> 11) as f64) / ((1_u64 << 53) as f64) - 0.5
}

fn random_matrix(rows: usize, columns: usize, seed: u64) -> Vec<f64> {
    let mut state = seed;
    (0..rows * columns).map(|_| lcg(&mut state)).collect()
}

fn output_fingerprint(values: &[f64]) -> u64 {
    let mut bytes = Vec::with_capacity(core::mem::size_of_val(values));
    for value in values {
        bytes.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    fnv1a64(&bytes)
}

fn gemm_golden_inputs() -> Vec<u8> {
    let mut bytes = b"fs-la:gemm-cross-isa-golden:v1".to_vec();
    push_text(&mut bytes, "gemm_f64");
    push_text(&mut bytes, "lcg-update-then-sample-v1");
    push_text(&mut bytes, "fnv1a64-f64-le-bits");
    for value in [
        u64::from(GEMM_BIT_SEMANTICS_VERSION),
        GOLDEN_M as u64,
        GOLDEN_N as u64,
        GOLDEN_K as u64,
        GOLDEN_A_SEED,
        GOLDEN_B_SEED,
        LCG_MULTIPLIER,
        LCG_INCREMENT,
        11,
        1_u64 << 53,
    ] {
        push_u64(&mut bytes, value);
    }
    for value in [-0.5_f64, 1.25, 0.0] {
        push_u64(&mut bytes, value.to_bits());
    }
    push_u64(&mut bytes, GEMM_GOLDEN);
    bytes
}

fn no_read_inputs() -> Vec<u8> {
    let mut bytes = b"fs-la:gemm-zero-no-read:v1".to_vec();
    push_u64(&mut bytes, u64::from(GEMM_BIT_SEMANTICS_VERSION));
    push_text(&mut bytes, "gemm_f64");

    push_text(&mut bytes, "beta-zero-overwrite");
    for value in [2_u64, 2, 3] {
        push_u64(&mut bytes, value);
    }
    for value in [1.0_f64, 0.0] {
        push_u64(&mut bytes, value.to_bits());
    }
    push_f64s(&mut bytes, &[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);
    push_f64s(&mut bytes, &[7.0, 8.0, 9.0, 10.0, 11.0, 12.0]);
    push_f64s(&mut bytes, &[f64::NAN; 4]);
    push_f64s(&mut bytes, &[58.0, 64.0, 139.0, 154.0]);

    push_text(&mut bytes, "alpha-zero-no-read");
    for value in [2_u64, 2, 3] {
        push_u64(&mut bytes, value);
    }
    for value in [0.0_f64, -0.5] {
        push_u64(&mut bytes, value.to_bits());
    }
    push_f64s(&mut bytes, &[f64::NAN; 6]);
    push_f64s(&mut bytes, &[f64::INFINITY; 6]);
    push_f64s(&mut bytes, &[2.0, 4.0, 6.0, 8.0]);
    push_f64s(&mut bytes, &[-1.0, -2.0, -3.0, -4.0]);
    bytes
}

fn factor_policy_inputs() -> Vec<u8> {
    let mut bytes = b"fs-la:factor-policy-refusals:v1".to_vec();
    push_u64(&mut bytes, u64::from(FACTOR_BIT_SEMANTICS_VERSION));

    push_text(&mut bytes, "lu-lowest-index-tie");
    push_len(&mut bytes, 2);
    push_f64s(&mut bytes, &[1.0, 2.0, -1.0, 4.0]);
    push_len(&mut bytes, 2);
    for index in [0_usize, 1] {
        push_len(&mut bytes, index);
    }

    push_text(&mut bytes, "lu-singular-refusal");
    push_len(&mut bytes, 2);
    push_f64s(&mut bytes, &[1.0, 2.0, 2.0, 4.0]);
    push_text(&mut bytes, "Singular");
    push_len(&mut bytes, 1);

    push_text(&mut bytes, "cholesky-not-spd-refusal");
    push_len(&mut bytes, 2);
    push_f64s(&mut bytes, &[4.0, 10.0, 10.0, 4.0]);
    push_text(&mut bytes, "NotSpd");
    push_len(&mut bytes, 1);
    bytes
}

fn first_bit_mismatch(computed: &[f64], reference: &[f64]) -> Option<(usize, u64, u64)> {
    computed
        .iter()
        .zip(reference)
        .enumerate()
        .find_map(|(index, (computed, reference))| {
            (computed.to_bits() != reference.to_bits()).then_some((
                index,
                computed.to_bits(),
                reference.to_bits(),
            ))
        })
}

fn gemm_golden_outcome() -> CaseOutcome {
    let a = random_matrix(GOLDEN_M, GOLDEN_K, GOLDEN_A_SEED);
    let b = random_matrix(GOLDEN_K, GOLDEN_N, GOLDEN_B_SEED);
    let mut c = vec![0.0; GOLDEN_M * GOLDEN_N];
    gemm_f64(GOLDEN_M, GOLDEN_N, GOLDEN_K, 1.25, &a, &b, 0.0, &mut c);
    let computed = output_fingerprint(&c);
    if computed != GEMM_GOLDEN {
        return CaseOutcome::fail(format!(
            "bit_semantics_version={GEMM_BIT_SEMANTICS_VERSION}; execution_tier={}; m={GOLDEN_M}; n={GOLDEN_N}; k={GOLDEN_K}; a_seed=0x{GOLDEN_A_SEED:016x}; b_seed=0x{GOLDEN_B_SEED:016x}; alpha=1.25; beta=0; computed=0x{computed:016x}; reference=0x{GEMM_GOLDEN:016x}",
            gemm_execution_tier(),
        ))
        .with_evidence("crates/fs-la/CONTRACT.md#invariants");
    }
    CaseOutcome::pass(format!(
        "bit_semantics_version={GEMM_BIT_SEMANTICS_VERSION}; execution_tier={}; m={GOLDEN_M}; n={GOLDEN_N}; k={GOLDEN_K}; fingerprint=0x{computed:016x}",
        gemm_execution_tier(),
    ))
    .with_evidence("crates/fs-la/CONTRACT.md#invariants")
}

fn zero_no_read_outcome() -> CaseOutcome {
    let a = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
    let b = [7.0, 8.0, 9.0, 10.0, 11.0, 12.0];
    let mut overwrite = [f64::NAN; 4];
    gemm_f64(2, 2, 3, 1.0, &a, &b, 0.0, &mut overwrite);
    let overwrite_reference = [58.0, 64.0, 139.0, 154.0];
    if let Some((index, computed, reference)) = first_bit_mismatch(&overwrite, &overwrite_reference)
    {
        return CaseOutcome::fail(format!(
            "subcase=beta-zero-overwrite; index={index}; computed_bits=0x{computed:016x}; reference_bits=0x{reference:016x}; initial_c=canonical_nan"
        ))
        .with_evidence("crates/fs-la/CONTRACT.md#public-types-and-semantics");
    }

    let poisoned_a = [f64::NAN; 6];
    let poisoned_b = [f64::INFINITY; 6];
    let mut beta_only = [2.0, 4.0, 6.0, 8.0];
    gemm_f64(2, 2, 3, 0.0, &poisoned_a, &poisoned_b, -0.5, &mut beta_only);
    let beta_only_reference = [-1.0, -2.0, -3.0, -4.0];
    if let Some((index, computed, reference)) = first_bit_mismatch(&beta_only, &beta_only_reference)
    {
        return CaseOutcome::fail(format!(
            "subcase=alpha-zero-no-read; index={index}; computed_bits=0x{computed:016x}; reference_bits=0x{reference:016x}; a=canonical_nan; b=positive_infinity; beta=-0.5"
        ))
        .with_evidence("crates/fs-la/CONTRACT.md#public-types-and-semantics");
    }

    CaseOutcome::pass("beta_zero_nan_overwrite=exact; alpha_zero_poison_no_read=exact; fixtures=2")
        .with_evidence("crates/fs-la/CONTRACT.md#public-types-and-semantics")
}

fn factor_policy_outcome() -> CaseOutcome {
    let tied = match lu(&[1.0, 2.0, -1.0, 4.0], 2) {
        Ok(factor) => factor,
        Err(error) => {
            return CaseOutcome::fail(format!(
                "subcase=lu-lowest-index-tie; unexpected_error={error:?}"
            ))
            .with_evidence("crates/fs-la/CONTRACT.md#invariants");
        }
    };
    if tied.perm() != [0, 1] {
        return CaseOutcome::fail(format!(
            "subcase=lu-lowest-index-tie; computed_perm={:?}; reference_perm=[0, 1]",
            tied.perm(),
        ))
        .with_evidence("crates/fs-la/CONTRACT.md#invariants");
    }

    match lu(&[1.0, 2.0, 2.0, 4.0], 2) {
        Err(FactorError::Singular { index: 1 }) => {}
        Err(error) => {
            return CaseOutcome::fail(format!(
                "subcase=lu-singular-refusal; computed={error:?}; reference=Singular {{ index: 1 }}"
            ))
            .with_evidence("crates/fs-la/CONTRACT.md#error-model");
        }
        Ok(factor) => {
            return CaseOutcome::fail(format!(
                "subcase=lu-singular-refusal; unexpected_success_perm={:?}; reference=Singular {{ index: 1 }}",
                factor.perm(),
            ))
            .with_evidence("crates/fs-la/CONTRACT.md#error-model");
        }
    }

    match cholesky(&[4.0, 10.0, 10.0, 4.0], 2) {
        Err(FactorError::NotSpd { index: 1 }) => {}
        Err(error) => {
            return CaseOutcome::fail(format!(
                "subcase=cholesky-not-spd-refusal; computed={error:?}; reference=NotSpd {{ index: 1 }}"
            ))
            .with_evidence("crates/fs-la/CONTRACT.md#error-model");
        }
        Ok(_) => {
            return CaseOutcome::fail(
                "subcase=cholesky-not-spd-refusal; unexpected_success=true; reference=NotSpd { index: 1 }",
            )
            .with_evidence("crates/fs-la/CONTRACT.md#error-model");
        }
    }

    CaseOutcome::pass("lu_tie_perm=[0, 1]; singular_index=1; not_spd_index=1")
        .with_evidence("crates/fs-la/CONTRACT.md#invariants")
        .with_evidence("crates/fs-la/CONTRACT.md#error-model")
}

#[test]
fn bedrock_casebook_suite_emits_replay_complete_green_records() {
    let golden_digest = fnv1a64(&gemm_golden_inputs());
    let no_read_digest = fnv1a64(&no_read_inputs());
    let factor_digest = fnv1a64(&factor_policy_inputs());
    assert_eq!(golden_digest, 0xa1e0_32e1_8d72_3b70);
    assert_eq!(no_read_digest, 0xa897_39ad_64c6_1fc4);
    assert_eq!(factor_digest, 0xd7d2_c9d2_6943_e940);

    let report = Suite::new(SUITE)
        .case(
            "gemm-cross-isa-golden",
            golden_digest,
            ToleranceSpec::Exact,
            gemm_golden_outcome,
        )
        .case(
            "gemm-zero-no-read",
            no_read_digest,
            ToleranceSpec::Exact,
            zero_no_read_outcome,
        )
        .case(
            "factorization-policy-and-refusals",
            factor_digest,
            ToleranceSpec::Structural,
            factor_policy_outcome,
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
            "gemm-cross-isa-golden",
            "gemm-zero-no-read",
            "factorization-policy-and-refusals",
        ]
    );
    assert_eq!(
        report.records[2].json_line(),
        format!(
            concat!(
                "{{\"casebook\":{},\"suite\":\"fs-la/bedrock-conformance-v1\",",
                "\"case\":\"factorization-policy-and-refusals\",\"inputs_digest\":\"d7d2c9d26943e940\",",
                "\"tolerance\":\"structural\",\"pass\":true,",
                "\"details\":\"lu_tie_perm=[0, 1]; singular_index=1; not_spd_index=1\",",
                "\"evidence\":[\"crates/fs-la/CONTRACT.md#invariants\",",
                "\"crates/fs-la/CONTRACT.md#error-model\"]}}"
            ),
            CASEBOOK_RECORD_VERSION,
        ),
        "the structured factor-policy record schema and field order are contract"
    );
}

#[test]
fn disclosed_seeded_corruption_turns_the_casebook_suite_red() {
    const CORRUPTION_SEED: u64 = 0xF51A_0001;
    let bit = ((CORRUPTION_SEED >> 1) % 64) as u32;
    assert_eq!(bit, 0);
    let corrupted_golden = GEMM_GOLDEN ^ (1_u64 << bit);

    let golden_inputs = gemm_golden_inputs();
    let mut inputs = b"fs-la:seeded-gemm-golden-corruption:v1".to_vec();
    push_u64(&mut inputs, CORRUPTION_SEED);
    push_u64(&mut inputs, u64::from(bit));
    push_len(&mut inputs, golden_inputs.len());
    inputs.extend_from_slice(&golden_inputs);
    push_u64(&mut inputs, GEMM_GOLDEN);
    push_u64(&mut inputs, corrupted_golden);
    let inputs_digest = fnv1a64(&inputs);
    assert_eq!(inputs_digest, 0x8168_3ae6_ea29_894e);

    let report = Suite::new(SUITE)
        .case(
            "seeded-gemm-golden-corruption",
            inputs_digest,
            ToleranceSpec::Exact,
            move || {
                let a = random_matrix(GOLDEN_M, GOLDEN_K, GOLDEN_A_SEED);
                let b = random_matrix(GOLDEN_K, GOLDEN_N, GOLDEN_B_SEED);
                let mut c = vec![0.0; GOLDEN_M * GOLDEN_N];
                gemm_f64(
                    GOLDEN_M, GOLDEN_N, GOLDEN_K, 1.25, &a, &b, 0.0, &mut c,
                );
                let computed = output_fingerprint(&c);
                if computed == corrupted_golden {
                    CaseOutcome::pass("seeded corruption was not detected")
                } else {
                    CaseOutcome::fail(format!(
                        "seed=0x{CORRUPTION_SEED:016x}; bit={bit}; bit_semantics_version={GEMM_BIT_SEMANTICS_VERSION}; execution_tier={}; m={GOLDEN_M}; n={GOLDEN_N}; k={GOLDEN_K}; a_seed=0x{GOLDEN_A_SEED:016x}; b_seed=0x{GOLDEN_B_SEED:016x}; alpha=1.25; beta=0; computed=0x{computed:016x}; canonical=0x{GEMM_GOLDEN:016x}; corrupted=0x{corrupted_golden:016x}",
                        gemm_execution_tier(),
                    ))
                    .with_evidence("crates/fs-la/tests/conformance.rs#seeded-corruption")
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
    assert_eq!(failure.case, "seeded-gemm-golden-corruption");
    assert_eq!(failure.inputs_digest, "81683ae6ea29894e");
    assert!(
        failure
            .details
            .contains(&format!("seed=0x{CORRUPTION_SEED:016x}"))
    );
    assert!(failure.details.contains(&format!("bit={bit}")));
    assert!(failure.details.contains("computed=0x"));
    assert!(failure.details.contains("canonical=0x1d7aa3c6b6317ef0"));
    assert!(failure.details.contains("corrupted=0x1d7aa3c6b6317ef1"));
    let line = failure.json_line();
    assert!(line.contains("\"tolerance\":\"exact\",\"pass\":false"));

    let panic = std::panic::catch_unwind(|| report.assert_green())
        .expect_err("the merge-gate assertion must reject the seeded failure");
    let message = panic
        .downcast_ref::<String>()
        .map(String::as_str)
        .or_else(|| panic.downcast_ref::<&str>().copied())
        .expect("casebook panic carries text");
    assert!(message.contains("seeded-gemm-golden-corruption"));
    assert!(message.contains(&format!("seed=0x{CORRUPTION_SEED:016x}")));
}
