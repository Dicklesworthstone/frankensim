//! Structured SUBSTRATE conformance records for fs-simd.
//!
//! The crate's source-level batteries retain broad tail, exceptional-value,
//! geometry-refusal, and capsule-equivalence coverage. These cases expose the
//! load-bearing portable laws through fs-casebook so a central run can diagnose
//! and replay a failure from its structured record (Gauntlet G0).

use fs_casebook::{CASEBOOK_RECORD_VERSION, CaseOutcome, Suite, ToleranceSpec, fnv1a64};
use fs_simd::{mk8x4_f64_tier_for, ops, scalar};
use fs_substrate::SimdTier;

const SUITE: &str = "fs-simd/substrate-conformance-v1";
const LCG_MULTIPLIER: u64 = 6_364_136_223_846_793_005;
const LCG_INCREMENT: u64 = 1_442_695_040_888_963_407;
const VECTOR_LENGTHS: [usize; 10] = [0, 1, 3, 4, 7, 8, 9, 17, 31, 65];
const VECTOR_SEEDS: [u64; 3] = [1, 42, 0xDEAD];
const TILE_DEPTHS: [usize; 4] = [0, 1, 3, 7];
const TILE_SEEDS: [u64; 2] = [7, 0xC0DE];
const SELECTOR_ROWS: [(SimdTier, bool, bool, SimdTier); 16] = [
    (SimdTier::Scalar, false, false, SimdTier::Scalar),
    (SimdTier::Scalar, false, true, SimdTier::Scalar),
    (SimdTier::Scalar, true, false, SimdTier::Scalar),
    (SimdTier::Scalar, true, true, SimdTier::Scalar),
    (SimdTier::Neon, false, false, SimdTier::Neon),
    (SimdTier::Neon, false, true, SimdTier::Neon),
    (SimdTier::Neon, true, false, SimdTier::Neon),
    (SimdTier::Neon, true, true, SimdTier::Neon),
    (SimdTier::Avx2, false, false, SimdTier::Scalar),
    (SimdTier::Avx2, false, true, SimdTier::Scalar),
    (SimdTier::Avx2, true, false, SimdTier::Scalar),
    (SimdTier::Avx2, true, true, SimdTier::Avx2),
    (SimdTier::Avx512, false, false, SimdTier::Scalar),
    (SimdTier::Avx512, false, true, SimdTier::Scalar),
    (SimdTier::Avx512, true, false, SimdTier::Scalar),
    (SimdTier::Avx512, true, true, SimdTier::Avx2),
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

fn push_f64s(bytes: &mut Vec<u8>, values: &[f64]) {
    push_len(bytes, values.len());
    for value in values {
        push_u64(bytes, value.to_bits());
    }
}

fn tier_code(tier: SimdTier) -> u8 {
    match tier {
        SimdTier::Scalar => 0,
        SimdTier::Neon => 1,
        SimdTier::Avx2 => 2,
        SimdTier::Avx512 => 3,
    }
}

fn primitive_inputs() -> Vec<u8> {
    let x = [1.0, -2.0, 4.0, 0.5];
    let y = [3.0, 5.0, -1.0, 2.0];
    let c = [-1.0, 2.0, 3.0, -4.0];
    let mut bytes = b"fs-simd:primitive-known-answers:v1".to_vec();
    push_f64s(&mut bytes, &x);
    push_f64s(&mut bytes, &y);
    push_f64s(&mut bytes, &c);
    for scalar in [2.0_f64, -0.5, 1.0] {
        push_u64(&mut bytes, scalar.to_bits());
    }
    for expected in [
        [5.0, 1.0, 7.0, 3.0],
        [-0.5, 1.0, -2.0, -0.25],
        [3.0, -10.0, -4.0, 1.0],
        [2.0, -8.0, -1.0, -3.0],
        [4.0, -9.0, -3.0, 2.0],
    ] {
        push_f64s(&mut bytes, &expected);
    }
    for expected in [-10.0_f64, 3.5] {
        push_u64(&mut bytes, expected.to_bits());
    }
    bytes
}

fn selector_inputs() -> Vec<u8> {
    let mut bytes = b"fs-simd:mk8x4-tier-selector:v1".to_vec();
    for (global, avx2, fma, expected) in SELECTOR_ROWS {
        bytes.extend_from_slice(&[
            tier_code(global),
            u8::from(avx2),
            u8::from(fma),
            tier_code(expected),
        ]);
    }
    bytes
}

fn dispatch_equivalence_inputs() -> Vec<u8> {
    let mut bytes = b"fs-simd:bounded-dispatch-table-scalar-equivalence:v1".to_vec();
    for value in [
        LCG_MULTIPLIER,
        LCG_INCREMENT,
        1,           // Seed oddness mask.
        11,          // Mantissa shift.
        1_u64 << 53, // Unit-interval denominator.
        7,           // Index scale modulus.
        1,           // Index scale bias.
        0x9E37,      // y-vector seed salt.
        0xC6A4,      // z-vector seed salt.
        0xA5A5,      // B-panel seed salt.
        0x5A5A,      // Accumulator seed salt.
        8,           // A-panel rows.
        4,           // B-panel columns.
        32,          // Initial accumulator elements.
    ] {
        push_u64(&mut bytes, value);
    }
    for value in [-0.5_f64, 1.25, -0.75] {
        push_u64(&mut bytes, value.to_bits());
    }
    push_len(&mut bytes, VECTOR_LENGTHS.len());
    for len in VECTOR_LENGTHS {
        push_len(&mut bytes, len);
    }
    push_len(&mut bytes, VECTOR_SEEDS.len());
    for seed in VECTOR_SEEDS {
        push_u64(&mut bytes, seed);
    }
    push_len(&mut bytes, TILE_DEPTHS.len());
    for depth in TILE_DEPTHS {
        push_len(&mut bytes, depth);
    }
    push_len(&mut bytes, TILE_SEEDS.len());
    for seed in TILE_SEEDS {
        push_u64(&mut bytes, seed);
    }
    for operation in ["axpy", "scale", "mul_elem", "fma3", "fmacc", "mk8x4_f64"] {
        push_len(&mut bytes, operation.len());
        bytes.extend_from_slice(operation.as_bytes());
    }
    bytes
}

fn generated_values(len: usize, seed: u64) -> Vec<f64> {
    let mut state = seed | 1;
    (0..len)
        .map(|index| {
            state = state
                .wrapping_mul(LCG_MULTIPLIER)
                .wrapping_add(LCG_INCREMENT);
            let unit = ((state >> 11) as f64) / ((1_u64 << 53) as f64);
            (unit - 0.5) * ((index % 7 + 1) as f64)
        })
        .collect()
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

fn exact_vector_outcome() -> CaseOutcome {
    let table = ops();
    let x = [1.0, -2.0, 4.0, 0.5];
    let y = [3.0, 5.0, -1.0, 2.0];
    let c = [-1.0, 2.0, 3.0, -4.0];

    let mut axpy = y;
    (table.axpy)(2.0, &x, &mut axpy);
    let mut scaled = x;
    (table.scale)(-0.5, &mut scaled);
    let mut product = [0.0; 4];
    (table.mul_elem)(&x, &y, &mut product);
    let mut fused = [0.0; 4];
    (table.fma3)(&x, &y, &c, &mut fused);
    let mut accumulated = [1.0; 4];
    (table.fmacc)(&x, &y, &mut accumulated);

    let vectors = [
        ("axpy", axpy.as_slice(), [5.0, 1.0, 7.0, 3.0]),
        ("scale", scaled.as_slice(), [-0.5, 1.0, -2.0, -0.25]),
        ("mul_elem", product.as_slice(), [3.0, -10.0, -4.0, 1.0]),
        ("fma3", fused.as_slice(), [2.0, -8.0, -1.0, -3.0]),
        ("fmacc", accumulated.as_slice(), [4.0, -9.0, -3.0, 2.0]),
    ];
    for (operation, computed, reference) in vectors {
        if let Some((index, computed_bits, reference_bits)) =
            first_bit_mismatch(computed, &reference)
        {
            return CaseOutcome::fail(format!(
                "dispatch_tier={}; operation={operation}; index={index}; computed_bits=0x{computed_bits:016x}; reference_bits=0x{reference_bits:016x}",
                table.tier.name(),
            ))
            .with_evidence("crates/fs-simd/CONTRACT.md#invariants");
        }
    }

    let dot = (table.dot)(&x, &y);
    let sum = (table.sum)(&x);
    if dot.to_bits() != (-10.0_f64).to_bits() || sum.to_bits() != 3.5_f64.to_bits() {
        return CaseOutcome::fail(format!(
            "dispatch_tier={}; operation=reductions; dot_bits=0x{:016x}; dot_reference_bits=0x{:016x}; sum_bits=0x{:016x}; sum_reference_bits=0x{:016x}",
            table.tier.name(),
            dot.to_bits(),
            (-10.0_f64).to_bits(),
            sum.to_bits(),
            3.5_f64.to_bits(),
        ))
        .with_evidence("crates/fs-simd/CONTRACT.md#invariants");
    }

    CaseOutcome::pass(format!(
        "dispatch_tier={}; vector_operations=5; reductions=2; bit_mismatches=0",
        table.tier.name(),
    ))
    .with_evidence("crates/fs-simd/CONTRACT.md#invariants")
}

fn selector_outcome() -> CaseOutcome {
    for (row, (global, avx2, fma, expected)) in SELECTOR_ROWS.into_iter().enumerate() {
        let computed = mk8x4_f64_tier_for(global, avx2, fma);
        if computed != expected {
            return CaseOutcome::fail(format!(
                "row={row}; global={}; avx2_available={avx2}; fma_available={fma}; computed={}; reference={}",
                global.name(),
                computed.name(),
                expected.name(),
            ))
            .with_evidence("crates/fs-simd/CONTRACT.md#public-types-and-semantics");
        }
    }
    CaseOutcome::pass(format!(
        "truth_table_rows={}; mismatches=0",
        SELECTOR_ROWS.len()
    ))
    .with_evidence("crates/fs-simd/CONTRACT.md#public-types-and-semantics")
}

fn compare_operation(
    dispatch_tier: SimdTier,
    operation: &str,
    len: usize,
    seed: u64,
    computed: &[f64],
    reference: &[f64],
) -> Option<CaseOutcome> {
    first_bit_mismatch(computed, reference).map(|(index, computed_bits, reference_bits)| {
        CaseOutcome::fail(format!(
            "dispatch_tier={}; operation={operation}; len={len}; seed=0x{seed:016x}; index={index}; computed_bits=0x{computed_bits:016x}; reference_bits=0x{reference_bits:016x}",
            dispatch_tier.name(),
        ))
        .with_evidence("crates/fs-simd/CONTRACT.md#determinism-class")
    })
}

#[allow(clippy::too_many_lines)] // one compact record exercises every bitwise dispatch facade
fn dispatch_equivalence_outcome() -> CaseOutcome {
    let table = ops();
    let dispatch_tier = table.tier;
    for len in VECTOR_LENGTHS {
        for seed in VECTOR_SEEDS {
            let x = generated_values(len, seed);
            let y = generated_values(len, seed ^ 0x9E37);
            let z = generated_values(len, seed ^ 0xC6A4);

            let mut computed = y.clone();
            let mut reference = y.clone();
            (table.axpy)(1.25, &x, &mut computed);
            scalar::axpy(1.25, &x, &mut reference);
            if let Some(failure) =
                compare_operation(dispatch_tier, "axpy", len, seed, &computed, &reference)
            {
                return failure;
            }

            let mut computed = x.clone();
            let mut reference = x.clone();
            (table.scale)(-0.75, &mut computed);
            scalar::scale(-0.75, &mut reference);
            if let Some(failure) =
                compare_operation(dispatch_tier, "scale", len, seed, &computed, &reference)
            {
                return failure;
            }

            let mut computed = vec![0.0; len];
            let mut reference = vec![0.0; len];
            (table.mul_elem)(&x, &y, &mut computed);
            scalar::mul_elem(&x, &y, &mut reference);
            if let Some(failure) =
                compare_operation(dispatch_tier, "mul_elem", len, seed, &computed, &reference)
            {
                return failure;
            }

            (table.fma3)(&x, &y, &z, &mut computed);
            scalar::fma3(&x, &y, &z, &mut reference);
            if let Some(failure) =
                compare_operation(dispatch_tier, "fma3", len, seed, &computed, &reference)
            {
                return failure;
            }

            computed.copy_from_slice(&z);
            reference.copy_from_slice(&z);
            (table.fmacc)(&x, &y, &mut computed);
            scalar::fmacc(&x, &y, &mut reference);
            if let Some(failure) =
                compare_operation(dispatch_tier, "fmacc", len, seed, &computed, &reference)
            {
                return failure;
            }
        }
    }

    for depth in TILE_DEPTHS {
        for seed in TILE_SEEDS {
            let a = generated_values(depth * 8, seed);
            let b = generated_values(depth * 4, seed ^ 0xA5A5);
            let initial = generated_values(32, seed ^ 0x5A5A);
            let mut computed = [[0.0; 4]; 8];
            let mut reference = [[0.0; 4]; 8];
            for ((computed_row, reference_row), initial_row) in computed
                .iter_mut()
                .zip(&mut reference)
                .zip(initial.chunks_exact(4))
            {
                computed_row.copy_from_slice(initial_row);
                reference_row.copy_from_slice(initial_row);
            }
            (table.mk8x4_f64)(&a, &b, depth, &mut computed);
            scalar::mk8x4_f64(&a, &b, depth, &mut reference);
            for (row, (computed_row, reference_row)) in computed.iter().zip(&reference).enumerate()
            {
                if let Some((column, computed_bits, reference_bits)) =
                    first_bit_mismatch(computed_row, reference_row)
                {
                    return CaseOutcome::fail(format!(
                        "dispatch_tier={}; mk8x4_effective_tier={}; operation=mk8x4_f64; depth={depth}; seed=0x{seed:016x}; row={row}; column={column}; computed_bits=0x{computed_bits:016x}; reference_bits=0x{reference_bits:016x}",
                        dispatch_tier.name(),
                        table.mk8x4_f64_tier.name(),
                    ))
                    .with_evidence("crates/fs-simd/CONTRACT.md#determinism-class");
                }
            }
        }
    }

    CaseOutcome::pass(format!(
        "dispatch_tier={}; mk8x4_effective_tier={}; vector_fixtures=30; vector_operations=5; tile_fixtures=8; bit_mismatches=0",
        dispatch_tier.name(),
        table.mk8x4_f64_tier.name(),
    ))
    .with_evidence("crates/fs-simd/CONTRACT.md#determinism-class")
}

#[test]
fn substrate_casebook_suite_emits_replay_complete_green_records() {
    let primitive_digest = fnv1a64(&primitive_inputs());
    let selector_digest = fnv1a64(&selector_inputs());
    let dispatch_digest = fnv1a64(&dispatch_equivalence_inputs());
    assert_eq!(primitive_digest, 0x5e67_d0e4_131e_ad10);
    assert_eq!(selector_digest, 0xb1a1_3db7_757a_39b4);
    assert_eq!(dispatch_digest, 0xe824_3aa2_ef9d_116e);

    let report = Suite::new(SUITE)
        .case(
            "primitive-known-answers",
            primitive_digest,
            ToleranceSpec::Exact,
            exact_vector_outcome,
        )
        .case(
            "mk8x4-tier-admission-truth-table",
            selector_digest,
            ToleranceSpec::Structural,
            selector_outcome,
        )
        .case(
            "bounded-dispatch-table-scalar-bit-equivalence",
            dispatch_digest,
            ToleranceSpec::Exact,
            dispatch_equivalence_outcome,
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
            "primitive-known-answers",
            "mk8x4-tier-admission-truth-table",
            "bounded-dispatch-table-scalar-bit-equivalence",
        ]
    );
    assert_eq!(
        report.records[1].json_line(),
        format!(
            concat!(
                "{{\"casebook\":{},\"suite\":\"fs-simd/substrate-conformance-v1\",",
                "\"case\":\"mk8x4-tier-admission-truth-table\",\"inputs_digest\":\"b1a13db7757a39b4\",",
                "\"tolerance\":\"structural\",\"pass\":true,",
                "\"details\":\"truth_table_rows=16; mismatches=0\",",
                "\"evidence\":[\"crates/fs-simd/CONTRACT.md#public-types-and-semantics\"]}}"
            ),
            CASEBOOK_RECORD_VERSION,
        ),
        "the structured selector record schema and field order are contract"
    );
}

#[test]
fn disclosed_seeded_corruption_turns_the_casebook_suite_red() {
    const CORRUPTION_SEED: u64 = 0xF51D_0001;
    let x = [1.0, -2.0, 4.0, 0.5];
    let mut computed = [3.0, 5.0, -1.0, 2.0];
    (ops().axpy)(2.0, &x, &mut computed);
    let mut corrupted_reference = [5.0_f64, 1.0, 7.0, 3.0].map(f64::to_bits);
    let word = (CORRUPTION_SEED as usize) % corrupted_reference.len();
    let bit = ((CORRUPTION_SEED.rotate_left(13) ^ CORRUPTION_SEED.rotate_right(11)) % 64) as u32;
    corrupted_reference[word] ^= 1_u64 << bit;
    assert_eq!(word, 1);
    assert_eq!(bit, 32);

    let mut inputs = b"fs-simd:seeded-known-answer-corruption:v1".to_vec();
    push_u64(&mut inputs, CORRUPTION_SEED);
    push_len(&mut inputs, "axpy".len());
    inputs.extend_from_slice(b"axpy");
    push_u64(&mut inputs, 2.0_f64.to_bits());
    push_len(&mut inputs, word);
    push_u64(&mut inputs, u64::from(bit));
    push_f64s(&mut inputs, &x);
    push_f64s(&mut inputs, &[3.0, 5.0, -1.0, 2.0]);
    for bits in corrupted_reference {
        push_u64(&mut inputs, bits);
    }
    let inputs_digest = fnv1a64(&inputs);
    assert_eq!(inputs_digest, 0x3450_7883_6384_1a85);

    let report = Suite::new(SUITE)
        .case(
            "seeded-primitive-reference-corruption",
            inputs_digest,
            ToleranceSpec::Exact,
            move || {
                let observed = computed.map(f64::to_bits);
                if observed == corrupted_reference {
                    CaseOutcome::pass("seeded corruption was not detected")
                } else {
                    CaseOutcome::fail(format!(
                        "seed=0x{CORRUPTION_SEED:016x}; operation=axpy; word={word}; bit={bit}; computed={observed:016x?}; corrupted_reference={corrupted_reference:016x?}"
                    ))
                    .with_evidence("crates/fs-simd/tests/conformance.rs#seeded-corruption")
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
    assert_eq!(failure.case, "seeded-primitive-reference-corruption");
    assert_eq!(failure.inputs_digest, "3450788363841a85");
    assert!(
        failure
            .details
            .contains(&format!("seed=0x{CORRUPTION_SEED:016x}"))
    );
    assert!(failure.details.contains("operation=axpy"));
    assert!(failure.details.contains(&format!("word={word}; bit={bit}")));
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
    assert!(message.contains("seeded-primitive-reference-corruption"));
    assert!(message.contains(&format!("seed=0x{CORRUPTION_SEED:016x}")));
}
