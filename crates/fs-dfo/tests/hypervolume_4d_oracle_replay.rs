//! G0/G3 independent-oracle replay at the exact hypervolume ceiling.
//!
//! Production `hypervolume` uses recursive objective-axis slicing. This target
//! checks its documented four-objective ceiling with a separately implemented
//! union-of-boxes inclusion-exclusion oracle over small dyadic fronts. It also
//! exercises public filtering policy, affine covariance, same-build replay,
//! and a disclosed one-bit input mutation that stale and frozen retained
//! receipts refuse.
//! The filtering cases replay the malformed/outside-point ignore behavior
//! already recorded by `fs-dfo/CONTRACT.md` and `moo_battery`; they do not
//! extend that behavior to non-finite data.
//!
//! This target makes no claim about exactness above four objectives, non-finite
//! inputs or references, large-front complexity, optimizer or archive quality,
//! Monte Carlo accuracy, authenticated admission, cross-ISA equality,
//! cancellation, or performance.

#![deny(unsafe_code)]

use fs_dfo::hypervolume;
use fs_rand::StreamKey;

const DIMENSION: usize = 4;
const REFERENCE: [f64; DIMENSION] = [1.0; DIMENSION];
const FIXTURE_SEEDS: [u64; 6] = [3, 17, 41, 97, 257, 65_537];
const FRONT_KERNEL: u32 = 0x4D48;
const MUTATION_COORDINATE: usize = 2;
const MUTATION_MANTISSA_BIT: u32 = 24;
const GOLDEN_ORIGINAL_INPUT_DIGEST: u64 = 0x8365_3e79_2771_b19e;
const GOLDEN_MUTANT_INPUT_DIGEST: u64 = 0xb280_e021_948d_de21;
const GOLDEN_ORIGINAL_VOLUME_BITS: u64 = 0x3fb0_0000_0000_0000;
const GOLDEN_MUTANT_VOLUME_BITS: u64 = 0x3faf_ffff_fe00_0000;
const GOLDEN_RED_LINE: &str = "{\"suite\":\"fs-dfo/hypervolume-4d-oracle-v1\",\"case\":\"fixed-front-bit-mutation\",\"verdict\":\"red\",\"target\":\"front[0][2]\",\"mantissa_bit\":24,\"before\":\"0x3fe0000000000000\",\"after\":\"0x3fe0000001000000\",\"input_before\":\"0x83653e792771b19e\",\"input_after\":\"0xb280e021948dde21\",\"volume_before\":\"0x3fb0000000000000\",\"volume_after\":\"0x3faffffffe000000\",\"stale_gate\":\"InputIdentityMismatch\",\"retained_gate\":\"RetainedVolumeMismatch\"}";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct HvReceipt {
    input_digest: u64,
    volume_bits: u64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ReceiptError {
    InputIdentityMismatch { expected: u64, found: u64 },
    RetainedVolumeMismatch { expected: u64, found: u64 },
}

impl HvReceipt {
    fn seal(front: &[Vec<f64>], reference: &[f64; DIMENSION]) -> Self {
        Self {
            input_digest: fnv1a64(&canonical_input_bytes(front, reference)),
            volume_bits: hypervolume(front, reference).to_bits(),
        }
    }

    fn validate_input(
        &self,
        front: &[Vec<f64>],
        reference: &[f64; DIMENSION],
    ) -> Result<(), ReceiptError> {
        let found = fnv1a64(&canonical_input_bytes(front, reference));
        if found == self.input_digest {
            Ok(())
        } else {
            Err(ReceiptError::InputIdentityMismatch {
                expected: self.input_digest,
                found,
            })
        }
    }

    fn admit_against(&self, retained_volume_bits: u64) -> Result<(), ReceiptError> {
        if self.volume_bits == retained_volume_bits {
            Ok(())
        } else {
            Err(ReceiptError::RetainedVolumeMismatch {
                expected: retained_volume_bits,
                found: self.volume_bits,
            })
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Mutation {
    coordinate: usize,
    mantissa_bit: u32,
    before: u64,
    after: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct MutationRun {
    mutation: Mutation,
    original: HvReceipt,
    mutant: HvReceipt,
    stale_error: ReceiptError,
    retained_error: ReceiptError,
    oracle_bits: u64,
    red_line: String,
}

fn push_u64(bytes: &mut Vec<u8>, value: u64) {
    bytes.extend_from_slice(&value.to_le_bytes());
}

fn push_len(bytes: &mut Vec<u8>, value: usize) {
    push_u64(
        bytes,
        u64::try_from(value).expect("fixture cardinality fits u64"),
    );
}

fn canonical_input_bytes(front: &[Vec<f64>], reference: &[f64; DIMENSION]) -> Vec<u8> {
    let mut bytes = b"fs-dfo-exact-hypervolume-4d-input-v1".to_vec();
    push_len(&mut bytes, reference.len());
    for value in reference {
        push_u64(&mut bytes, value.to_bits());
    }
    push_len(&mut bytes, front.len());
    for point in front {
        push_len(&mut bytes, point.len());
        for value in point {
            push_u64(&mut bytes, value.to_bits());
        }
    }
    bytes
}

fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut hash = 0xcbf2_9ce4_8422_2325_u64;
    for byte in bytes {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

/// Independent union-of-boxes inclusion-exclusion oracle. Production instead
/// recursively slices on the final objective, so neither control flow nor
/// accumulation order is shared.
fn inclusion_exclusion_oracle(front: &[Vec<f64>], reference: &[f64; DIMENSION]) -> f64 {
    assert!(reference.iter().all(|value| value.is_finite()));
    assert!(front.iter().all(|point| {
        point.len() == DIMENSION
            && point.iter().all(|value| value.is_finite())
            && point
                .iter()
                .zip(reference)
                .all(|(value, limit)| value < limit)
    }));
    if front.is_empty() {
        return 0.0;
    }
    let shift = u32::try_from(front.len()).expect("small oracle front fits u32");
    let subset_count = 1usize
        .checked_shl(shift)
        .expect("small oracle front has a representable power set");
    let mut union = 0.0;
    for subset in 1..subset_count {
        let mut lower = [f64::NEG_INFINITY; DIMENSION];
        for (point_index, point) in front.iter().enumerate() {
            if subset & (1usize << point_index) != 0 {
                for (lower_value, point_value) in lower.iter_mut().zip(point) {
                    *lower_value = lower_value.max(*point_value);
                }
            }
        }
        let intersection = reference
            .iter()
            .zip(&lower)
            .map(|(limit, lower)| limit - lower)
            .product::<f64>();
        if subset.count_ones().is_multiple_of(2) {
            union -= intersection;
        } else {
            union += intersection;
        }
    }
    union
}

fn assert_exact_oracle(front: &[Vec<f64>], reference: &[f64; DIMENSION]) -> u64 {
    let production = hypervolume(front, reference);
    let oracle = inclusion_exclusion_oracle(front, reference);
    assert!(production.is_finite() && oracle.is_finite());
    assert_eq!(
        production.to_bits(),
        oracle.to_bits(),
        "production={production:.17e}; oracle={oracle:.17e}; front={front:?}; reference={reference:?}"
    );
    production.to_bits()
}

fn seeded_front(seed: u64, points: usize) -> Vec<Vec<f64>> {
    let tile = u32::try_from(points).expect("fixture point count fits u32");
    let mut stream = StreamKey {
        seed,
        kernel: FRONT_KERNEL,
        tile,
    }
    .stream();
    (0..points)
        .map(|_| {
            (0..DIMENSION)
                .map(|_| {
                    let numerator = u32::try_from(stream.next_below(15) + 1)
                        .expect("dyadic numerator fits u32");
                    f64::from(numerator) / 16.0
                })
                .collect()
        })
        .collect()
}

fn affine_transform(
    front: &[Vec<f64>],
    reference: &[f64; DIMENSION],
) -> (Vec<Vec<f64>>, [f64; DIMENSION]) {
    let scales = [2.0, 0.5, 4.0, 0.25];
    let shifts = [-1.0, 2.0, -3.0, 0.5];
    let transformed_front = front
        .iter()
        .map(|point| {
            point
                .iter()
                .enumerate()
                .map(|(objective, value)| shifts[objective] + scales[objective] * value)
                .collect()
        })
        .collect();
    let transformed_reference = core::array::from_fn(|objective| {
        shifts[objective] + scales[objective] * reference[objective]
    });
    (transformed_front, transformed_reference)
}

fn fixed_mutation_run() -> MutationRun {
    // The three untouched half-widths make every post-mutation product an
    // exact power-of-two scaling, so the independent and production
    // accumulation orders must agree bit for bit.
    let front = [vec![0.5; DIMENSION]];
    let original = HvReceipt::seal(&front, &REFERENCE);
    assert_eq!(
        original,
        HvReceipt {
            input_digest: GOLDEN_ORIGINAL_INPUT_DIGEST,
            volume_bits: GOLDEN_ORIGINAL_VOLUME_BITS,
        }
    );
    assert_eq!(
        original.volume_bits,
        inclusion_exclusion_oracle(&front, &REFERENCE).to_bits()
    );

    let mut mutated_front = front.clone();
    let before = mutated_front[0][MUTATION_COORDINATE].to_bits();
    let after = before ^ (1u64 << MUTATION_MANTISSA_BIT);
    mutated_front[0][MUTATION_COORDINATE] = f64::from_bits(after);
    assert!(mutated_front[0][MUTATION_COORDINATE].is_finite());
    assert!(mutated_front[0][MUTATION_COORDINATE] < REFERENCE[MUTATION_COORDINATE]);

    let stale_error = original
        .validate_input(&mutated_front, &REFERENCE)
        .expect_err("one-bit front mutation must miss stale input identity");
    let mutant = HvReceipt::seal(&mutated_front, &REFERENCE);
    assert_eq!(
        mutant,
        HvReceipt {
            input_digest: GOLDEN_MUTANT_INPUT_DIGEST,
            volume_bits: GOLDEN_MUTANT_VOLUME_BITS,
        }
    );
    let retained_error = mutant
        .admit_against(original.volume_bits)
        .expect_err("one-bit front mutation must miss retained volume");
    let oracle_bits = inclusion_exclusion_oracle(&mutated_front, &REFERENCE).to_bits();
    let red_line = format!(
        "{{\"suite\":\"fs-dfo/hypervolume-4d-oracle-v1\",\"case\":\"fixed-front-bit-mutation\",\"verdict\":\"red\",\"target\":\"front[0][{MUTATION_COORDINATE}]\",\"mantissa_bit\":{MUTATION_MANTISSA_BIT},\"before\":\"0x{before:016x}\",\"after\":\"0x{after:016x}\",\"input_before\":\"0x{:016x}\",\"input_after\":\"0x{:016x}\",\"volume_before\":\"0x{:016x}\",\"volume_after\":\"0x{:016x}\",\"stale_gate\":\"InputIdentityMismatch\",\"retained_gate\":\"RetainedVolumeMismatch\"}}",
        original.input_digest, mutant.input_digest, original.volume_bits, mutant.volume_bits,
    );
    assert_eq!(red_line, GOLDEN_RED_LINE);
    MutationRun {
        mutation: Mutation {
            coordinate: MUTATION_COORDINATE,
            mantissa_bit: MUTATION_MANTISSA_BIT,
            before,
            after,
        },
        original,
        mutant,
        stale_error,
        retained_error,
        oracle_bits,
        red_line,
    }
}

#[test]
fn g0_four_dimensional_closed_forms_match_independent_oracle() {
    let singleton = [vec![0.5; DIMENSION]];
    assert_eq!(
        assert_exact_oracle(&singleton, &REFERENCE),
        0.0625f64.to_bits()
    );

    let pair = [vec![0.25, 0.5, 0.75, 0.5], vec![0.5, 0.25, 0.5, 0.75]];
    assert_eq!(
        assert_exact_oracle(&pair, &REFERENCE),
        0.078125f64.to_bits()
    );
}

#[test]
fn g3_seeded_fronts_preserve_union_semantics_under_public_metamorphisms() {
    for (seed_index, seed) in FIXTURE_SEEDS.into_iter().enumerate() {
        let front = seeded_front(seed, seed_index + 2);
        let baseline = assert_exact_oracle(&front, &REFERENCE);
        assert_eq!(hypervolume(&front, &REFERENCE).to_bits(), baseline);

        let mut permuted = front.clone();
        permuted.reverse();
        assert_eq!(assert_exact_oracle(&permuted, &REFERENCE), baseline);

        let mut admitted_augmentation = front.clone();
        admitted_augmentation.push(front[0].clone());
        admitted_augmentation.push(front[0].iter().map(|value| 0.5 * (*value + 1.0)).collect());
        assert_eq!(
            assert_exact_oracle(&admitted_augmentation, &REFERENCE),
            baseline
        );

        let mut mixed_public_input = admitted_augmentation;
        mixed_public_input.push(vec![0.25, 0.5, 0.75]);
        mixed_public_input.push(vec![1.125, 0.0, 0.0, 0.0]);
        assert_eq!(
            hypervolume(&mixed_public_input, &REFERENCE).to_bits(),
            baseline,
            "recorded malformed and finite outside-reference points must be ignored"
        );

        let (transformed, transformed_reference) = affine_transform(&front, &REFERENCE);
        assert_eq!(
            assert_exact_oracle(&transformed, &transformed_reference),
            baseline,
            "the per-axis power-of-two Jacobian is exactly one"
        );
    }
}

#[test]
fn g3_same_build_fixed_inputs_repeat_every_four_dimensional_volume_bit() {
    for (seed_index, seed) in FIXTURE_SEEDS.into_iter().enumerate() {
        let front = seeded_front(seed, seed_index + 2);
        let first = HvReceipt::seal(&front, &REFERENCE);
        let second = HvReceipt::seal(&front, &REFERENCE);
        assert_eq!(first, second);
        assert_eq!(first.volume_bits, assert_exact_oracle(&front, &REFERENCE));
        assert!(first.validate_input(&front, &REFERENCE).is_ok());
        assert!(first.admit_against(first.volume_bits).is_ok());
    }
}

#[test]
fn g3_fixed_one_bit_front_mutation_matches_frozen_red_receipt() {
    let first = fixed_mutation_run();
    let second = fixed_mutation_run();

    assert_eq!(first, second);
    assert_ne!(first.original.input_digest, first.mutant.input_digest);
    assert_ne!(first.original.volume_bits, first.mutant.volume_bits);
    assert_eq!(first.mutant.volume_bits, first.oracle_bits);
    assert_eq!(first.original.input_digest, GOLDEN_ORIGINAL_INPUT_DIGEST);
    assert_eq!(first.mutant.input_digest, GOLDEN_MUTANT_INPUT_DIGEST);
    assert_eq!(first.original.volume_bits, GOLDEN_ORIGINAL_VOLUME_BITS);
    assert_eq!(first.mutant.volume_bits, GOLDEN_MUTANT_VOLUME_BITS);
    assert_eq!(
        first.mutation.before ^ first.mutation.after,
        1u64 << first.mutation.mantissa_bit
    );
    assert!(matches!(
        first.stale_error,
        ReceiptError::InputIdentityMismatch { expected, found }
            if expected == first.original.input_digest && found == first.mutant.input_digest
    ));
    assert!(matches!(
        first.retained_error,
        ReceiptError::RetainedVolumeMismatch { expected, found }
            if expected == first.original.volume_bits && found == first.mutant.volume_bits
    ));
    assert_eq!(first.red_line, GOLDEN_RED_LINE);
    assert!(
        first
            .red_line
            .contains(&format!("front[0][{}]", first.mutation.coordinate))
    );
    assert!(first.red_line.contains("InputIdentityMismatch"));
    assert!(first.red_line.contains("RetainedVolumeMismatch"));
}
