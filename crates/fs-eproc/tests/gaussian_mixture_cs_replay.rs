//! G0/G3/G5 replay receipt for [`fs_eproc::GaussianMixtureCs`].
//!
//! A fixed Philox coordinate rotates and orients a bank of complementary
//! dyadic observations. The production process is inspected after every
//! observation through its complete public surface: length, empty state,
//! interval, and a fixed symmetric lattice of null-mean e-values. A
//! test-local oracle accumulates the dyadic numerators as integers and uses a
//! separately arranged log-domain mixture formula, so it does not trust the
//! process's private sufficient statistics or either composite public method.
//!
//! This finite fixture is replay and certifier plumbing, not a general
//! coverage proof. It makes no all-law/configuration/seed/horizon,
//! adversarial-stopping, optimizer-integration, cross-process/cross-ISA,
//! cancellation, persistence/authentication, private-state, or performance
//! claim.

use fs_eproc::GaussianMixtureCs;
use fs_obs::ident::{IdentityBuilder, ReplayIdentity};
use fs_obs::{Emitter, Event, EventKind, Severity};
use fs_rand::StreamKey;
use std::panic::catch_unwind;

const SUITE: &str = "fs-eproc/gaussian-mixture-cs-replay-v1";
const CASE: &str = "complementary-dyadic-full-trajectory";
const RED_CASE: &str = "seeded-checkpoint-radius-corruption";

const INPUT_SEED: u64 = 0x6A75_5C50_2026_0718;
const INPUT_KERNEL: u32 = 0x6C51;
const INPUT_TILE: u32 = 0;
const INPUT_DRAWS: u64 = 2;

const MUTATION_SEED: u64 = 0xB17F_11F0_2026_0718;
const MUTATION_KERNEL: u32 = 0x6C52;
const MUTATION_TILE: u32 = 0;
const MUTATION_DRAWS: u64 = 2;

const SIGMA: f64 = 0.5;
const RHO: f64 = 8.0;
const ALPHA: f64 = 0.05;
const TRUE_MEAN_NUMERATOR: u64 = 8;
const DYADIC_DENOMINATOR: u64 = 16;
const ORACLE_ULPS: u64 = 128;
const BOUNDARY_REL_TOLERANCE: f64 = 2.0e-12;

const COMPLEMENTARY_PAIRS: [[u64; 2]; 8] = [
    [0, 16],
    [1, 15],
    [2, 14],
    [3, 13],
    [4, 12],
    [5, 11],
    [6, 10],
    [7, 9],
];

// Symmetric about 8/16 = 1/2. At each completed complementary pair the
// running center is exactly one half, making the symmetry comparison exact.
const NULL_NUMERATORS: [u64; 5] = [4, 6, 8, 10, 12];

#[derive(Debug, Clone, PartialEq, Eq)]
struct FixtureSchedule {
    rotation: usize,
    orientation_word: u64,
    draws: u64,
    observations: Vec<u64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct InitialState {
    len: u64,
    is_empty: bool,
    interval_is_none: bool,
    e_value_bits: [u64; NULL_NUMERATORS.len()],
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Checkpoint {
    ordinal: u64,
    observation_numerator: u64,
    observation_bits: u64,
    prefix_sum_numerator: u64,
    len: u64,
    is_empty: bool,
    center_bits: u64,
    radius_bits: u64,
    fixed_e_value_bits: [u64; NULL_NUMERATORS.len()],
    boundary_e_value_bits: [u64; 2],
    inner_e_value_bits: [u64; 2],
    outer_e_value_bits: [u64; 2],
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRecord {
    rotation: usize,
    orientation_word: u64,
    input_draws: u64,
    observations: Vec<u64>,
    initial: InitialState,
    checkpoints: Vec<Checkpoint>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct SealedStudy {
    record: StudyRecord,
    identity: ReplayIdentity,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AdmissionError {
    PayloadIdentityMismatch { declared: u64, computed: u64 },
    ReferenceIdentityMismatch { expected: u64, found: u64 },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Mutation {
    seed: u64,
    kernel: u32,
    tile: u32,
    selector_word: u64,
    bit_selector_word: u64,
    draws: u64,
    checkpoint: usize,
    mantissa_bit: u32,
    before: u64,
    after: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Corruption {
    mutant: SealedStudy,
    mutation: Mutation,
    stale_error: AdmissionError,
    reference_error: AdmissionError,
    semantic_mismatch: String,
}

#[derive(Debug, Clone, Copy)]
struct OracleCheckpoint {
    center: f64,
    radius: f64,
    fixed_e_values: [f64; NULL_NUMERATORS.len()],
}

fn usize_u64(value: usize) -> u64 {
    u64::try_from(value).expect("fixture sizes fit u64")
}

fn dyadic(numerator: u64) -> f64 {
    numerator as f64 / DYADIC_DENOMINATOR as f64
}

fn fixed_nulls() -> [f64; NULL_NUMERATORS.len()] {
    NULL_NUMERATORS.map(dyadic)
}

fn fixture_schedule() -> FixtureSchedule {
    let mut stream = StreamKey {
        seed: INPUT_SEED,
        kernel: INPUT_KERNEL,
        tile: INPUT_TILE,
    }
    .stream();
    let rotation = usize::try_from(
        stream.next_below(u64::try_from(COMPLEMENTARY_PAIRS.len()).expect("eight pairs")),
    )
    .expect("rotation fits usize");
    let orientation_word = stream.next_u64();
    let draws = stream.index();
    assert_eq!(
        draws, INPUT_DRAWS,
        "power-of-two rotation plus orientation must consume two draws"
    );
    let mut observations = Vec::with_capacity(2 * COMPLEMENTARY_PAIRS.len());
    for step in 0..COMPLEMENTARY_PAIRS.len() {
        let pair = COMPLEMENTARY_PAIRS[(rotation + step) % COMPLEMENTARY_PAIRS.len()];
        let reverse = ((orientation_word >> step) & 1) != 0;
        if reverse {
            observations.extend([pair[1], pair[0]]);
        } else {
            observations.extend(pair);
        }
    }
    FixtureSchedule {
        rotation,
        orientation_word,
        draws,
        observations,
    }
}

fn config_identity() -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-eproc-gaussian-mixture-cs-config-v1")
        .str("statistical-object", "Robbins-normal-mixture-confidence-sequence")
        .str("observation-units", "dimensionless-bounded-score")
        .str("mean-units", "dimensionless-bounded-score")
        .f64_bits("support-lower", 0.0)
        .f64_bits("support-upper", 1.0)
        .str("sub-Gaussian-justification", "Hoeffding sigma one-half for support zero-to-one")
        .f64_bits("sigma", SIGMA)
        .f64_bits("rho", RHO)
        .f64_bits("alpha", ALPHA)
        .u64("dyadic-denominator", DYADIC_DENOMINATOR)
        .u64("true-mean-numerator", TRUE_MEAN_NUMERATOR)
        .u64("input-seed", INPUT_SEED)
        .u64("input-kernel", u64::from(INPUT_KERNEL))
        .u64("input-tile", u64::from(INPUT_TILE))
        .u64("input-draws", INPUT_DRAWS)
        .u64(
            "stream-semantics-version",
            u64::from(fs_rand::STREAM_SEMANTICS_VERSION),
        )
        .str(
            "schedule",
            "next-below-eight rotation then one orientation bit per complementary pair",
        )
        .str(
            "oracle",
            "integer-dyadic sufficient statistics plus factored-radius and log-e formulas",
        )
        .u64("oracle-ulp-tolerance", ORACLE_ULPS)
        .f64_bits("boundary-relative-tolerance", BOUNDARY_REL_TOLERANCE)
        .u64("mutation-seed", MUTATION_SEED)
        .u64("mutation-kernel", u64::from(MUTATION_KERNEL))
        .u64("mutation-tile", u64::from(MUTATION_TILE))
        .u64("mutation-draws", MUTATION_DRAWS)
        .str("mutation-target", "one retained checkpoint radius significand bit")
        .u64("mutation-minimum-mantissa-bit", 12)
        .u64("mutation-exclusive-maximum-mantissa-bit", 24)
        .str("fs-eproc-version", fs_eproc::VERSION)
        .str("fs-math-version", fs_math::VERSION)
        .str("fs-obs-version", fs_obs::VERSION)
        .str("fs-rand-version", fs_rand::VERSION)
        .str(
            "no-claim",
            "finite-fixture-only-no-general-coverage-all-law-config-seed-horizon-adversarial-stopping-optimizer-cross-process-cross-ISA-Cx-persistence-auth-private-state-performance-claim",
        );
    for pair in COMPLEMENTARY_PAIRS {
        builder = builder
            .u64("pair-low-numerator", pair[0])
            .u64("pair-high-numerator", pair[1]);
    }
    for numerator in NULL_NUMERATORS {
        builder = builder.u64("fixed-null-numerator", numerator);
    }
    builder.finish()
}

fn result_identity(config: &ReplayIdentity, record: &StudyRecord) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-eproc-gaussian-mixture-cs-result-v1")
        .child("config", config)
        .u64("rotation", usize_u64(record.rotation))
        .u64("orientation-word", record.orientation_word)
        .u64("input-draws", record.input_draws)
        .u64("initial-len", record.initial.len)
        .flag("initial-is-empty", record.initial.is_empty)
        .flag("initial-interval-is-none", record.initial.interval_is_none);
    for bits in record.initial.e_value_bits {
        builder = builder.u64("initial-e-value-bits", bits);
    }
    for numerator in &record.observations {
        builder = builder.u64("observation-numerator", *numerator);
    }
    for checkpoint in &record.checkpoints {
        builder = builder
            .u64("checkpoint-ordinal", checkpoint.ordinal)
            .u64(
                "checkpoint-observation-numerator",
                checkpoint.observation_numerator,
            )
            .u64("checkpoint-observation-bits", checkpoint.observation_bits)
            .u64(
                "checkpoint-prefix-sum-numerator",
                checkpoint.prefix_sum_numerator,
            )
            .u64("checkpoint-len", checkpoint.len)
            .flag("checkpoint-is-empty", checkpoint.is_empty)
            .u64("checkpoint-center-bits", checkpoint.center_bits)
            .u64("checkpoint-radius-bits", checkpoint.radius_bits);
        for bits in checkpoint.fixed_e_value_bits {
            builder = builder.u64("checkpoint-fixed-e-value-bits", bits);
        }
        for bits in checkpoint.boundary_e_value_bits {
            builder = builder.u64("checkpoint-boundary-e-value-bits", bits);
        }
        for bits in checkpoint.inner_e_value_bits {
            builder = builder.u64("checkpoint-inner-e-value-bits", bits);
        }
        for bits in checkpoint.outer_e_value_bits {
            builder = builder.u64("checkpoint-outer-e-value-bits", bits);
        }
    }
    builder.finish()
}

impl SealedStudy {
    fn seal(config: &ReplayIdentity, record: StudyRecord) -> Self {
        let identity = result_identity(config, &record);
        Self { record, identity }
    }

    fn validate_payload(&self, config: &ReplayIdentity) -> Result<(), AdmissionError> {
        let computed = result_identity(config, &self.record);
        if computed == self.identity {
            Ok(())
        } else {
            Err(AdmissionError::PayloadIdentityMismatch {
                declared: self.identity.root(),
                computed: computed.root(),
            })
        }
    }

    fn admit_against(
        &self,
        config: &ReplayIdentity,
        reference: &ReplayIdentity,
    ) -> Result<(), AdmissionError> {
        self.validate_payload(config)?;
        if &self.identity == reference {
            Ok(())
        } else {
            Err(AdmissionError::ReferenceIdentityMismatch {
                expected: reference.root(),
                found: self.identity.root(),
            })
        }
    }
}

fn oracle_log_e_value(n: u64, sum_numerator: u64, null_mean: f64) -> f64 {
    let t = n as f64;
    let sum = dyadic(sum_numerator);
    let variance_clock = t.mul_add(SIGMA * SIGMA, RHO);
    let centered_sum = (-null_mean).mul_add(t, sum);
    0.5 * (fs_math::det::ln(RHO) - fs_math::det::ln(variance_clock))
        + centered_sum * centered_sum / (2.0 * variance_clock)
}

fn oracle_checkpoint(n: u64, sum_numerator: u64) -> OracleCheckpoint {
    let t = n as f64;
    let center = dyadic(sum_numerator) / t;
    let variance_clock = t.mul_add(SIGMA * SIGMA, RHO);
    let log_boundary = fs_math::det::ln(variance_clock / RHO) + 2.0 * fs_math::det::ln(1.0 / ALPHA);
    // Deliberately factored differently from the production expression.
    let radius = fs_math::det::sqrt(variance_clock) * fs_math::det::sqrt(log_boundary) / t;
    let fixed_e_values =
        fixed_nulls().map(|null| fs_math::det::exp(oracle_log_e_value(n, sum_numerator, null)));
    OracleCheckpoint {
        center,
        radius,
        fixed_e_values,
    }
}

fn query_pair(process: &GaussianMixtureCs, center: f64, offset: f64) -> [u64; 2] {
    [
        process.e_value_for(center - offset).to_bits(),
        process.e_value_for(center + offset).to_bits(),
    ]
}

fn run_study(config: &ReplayIdentity) -> SealedStudy {
    let schedule = fixture_schedule();
    let nulls = fixed_nulls();
    let mut process = GaussianMixtureCs::new(SIGMA, RHO, ALPHA);
    let initial = InitialState {
        len: process.len(),
        is_empty: process.is_empty(),
        interval_is_none: process.interval().is_none(),
        e_value_bits: nulls.map(|null| process.e_value_for(null).to_bits()),
    };

    let mut prefix_sum_numerator = 0u64;
    let mut checkpoints = Vec::with_capacity(schedule.observations.len());
    for (index, &numerator) in schedule.observations.iter().enumerate() {
        let observation = dyadic(numerator);
        process.observe(observation);
        prefix_sum_numerator = prefix_sum_numerator
            .checked_add(numerator)
            .expect("fixture numerator sum fits u64");
        let (center, radius) = process
            .interval()
            .expect("one observation produces an interval");
        checkpoints.push(Checkpoint {
            ordinal: usize_u64(index + 1),
            observation_numerator: numerator,
            observation_bits: observation.to_bits(),
            prefix_sum_numerator,
            len: process.len(),
            is_empty: process.is_empty(),
            center_bits: center.to_bits(),
            radius_bits: radius.to_bits(),
            fixed_e_value_bits: nulls.map(|null| process.e_value_for(null).to_bits()),
            boundary_e_value_bits: query_pair(&process, center, radius),
            inner_e_value_bits: query_pair(&process, center, 0.5 * radius),
            outer_e_value_bits: query_pair(&process, center, 1.5 * radius),
        });
    }

    SealedStudy::seal(
        config,
        StudyRecord {
            rotation: schedule.rotation,
            orientation_word: schedule.orientation_word,
            input_draws: schedule.draws,
            observations: schedule.observations,
            initial,
            checkpoints,
        },
    )
}

fn within_ulps(actual: f64, expected: f64, tolerance: u64) -> bool {
    actual.is_finite()
        && expected.is_finite()
        && actual.is_sign_negative() == expected.is_sign_negative()
        && actual.to_bits().abs_diff(expected.to_bits()) <= tolerance
}

fn relative_close(actual: f64, expected: f64, tolerance: f64) -> bool {
    actual.is_finite()
        && expected.is_finite()
        && (actual - expected).abs() <= tolerance * actual.abs().max(expected.abs()).max(1.0)
}

fn decode_pair(bits: [u64; 2]) -> [f64; 2] {
    bits.map(f64::from_bits)
}

#[allow(clippy::too_many_lines)]
fn semantic_mismatch(record: &StudyRecord) -> Option<String> {
    let expected_schedule = fixture_schedule();
    if record.rotation != expected_schedule.rotation {
        return Some("rotation".to_string());
    }
    if record.orientation_word != expected_schedule.orientation_word {
        return Some("orientation_word".to_string());
    }
    if record.input_draws != expected_schedule.draws {
        return Some("input_draws".to_string());
    }
    if record.observations != expected_schedule.observations {
        return Some("observations".to_string());
    }
    if record.initial.len != 0 {
        return Some("initial.len".to_string());
    }
    if !record.initial.is_empty {
        return Some("initial.is_empty".to_string());
    }
    if !record.initial.interval_is_none {
        return Some("initial.interval".to_string());
    }
    if record.initial.e_value_bits != [1.0f64.to_bits(); NULL_NUMERATORS.len()] {
        return Some("initial.e_values".to_string());
    }
    if record.checkpoints.len() != record.observations.len() {
        return Some("checkpoint_count".to_string());
    }

    let threshold = 1.0 / ALPHA;
    let mut prefix_sum_numerator = 0u64;
    for (index, checkpoint) in record.checkpoints.iter().enumerate() {
        let path = |field: &str| format!("checkpoints[{index}].{field}");
        let expected_numerator = record.observations[index];
        prefix_sum_numerator = prefix_sum_numerator
            .checked_add(expected_numerator)
            .expect("fixture numerator sum fits u64");
        let expected_len = usize_u64(index + 1);
        if checkpoint.ordinal != expected_len {
            return Some(path("ordinal"));
        }
        if checkpoint.observation_numerator != expected_numerator {
            return Some(path("observation_numerator"));
        }
        if checkpoint.observation_bits != dyadic(expected_numerator).to_bits() {
            return Some(path("observation_bits"));
        }
        if checkpoint.prefix_sum_numerator != prefix_sum_numerator {
            return Some(path("prefix_sum_numerator"));
        }
        if checkpoint.len != expected_len {
            return Some(path("len"));
        }
        if checkpoint.is_empty {
            return Some(path("is_empty"));
        }

        let oracle = oracle_checkpoint(expected_len, prefix_sum_numerator);
        let center = f64::from_bits(checkpoint.center_bits);
        let radius = f64::from_bits(checkpoint.radius_bits);
        if center.to_bits() != oracle.center.to_bits() {
            return Some(path("center"));
        }
        if !within_ulps(radius, oracle.radius, ORACLE_ULPS) || radius <= 0.0 {
            return Some(path("radius"));
        }
        for (null_index, (&actual_bits, expected)) in checkpoint
            .fixed_e_value_bits
            .iter()
            .zip(oracle.fixed_e_values)
            .enumerate()
        {
            if !within_ulps(f64::from_bits(actual_bits), expected, ORACLE_ULPS) {
                return Some(format!("checkpoints[{index}].fixed_e_values[{null_index}]"));
            }
        }

        let boundary = decode_pair(checkpoint.boundary_e_value_bits);
        if !boundary
            .into_iter()
            .all(|value| relative_close(value, threshold, BOUNDARY_REL_TOLERANCE))
        {
            return Some(path("boundary_e_values"));
        }
        let inner = decode_pair(checkpoint.inner_e_value_bits);
        if !inner.into_iter().all(|value| value < threshold) {
            return Some(path("inner_e_values"));
        }
        let outer = decode_pair(checkpoint.outer_e_value_bits);
        if !outer.into_iter().all(|value| value > threshold) {
            return Some(path("outer_e_values"));
        }

        if expected_len.is_multiple_of(2) {
            let values = checkpoint.fixed_e_value_bits;
            if center.to_bits() != dyadic(TRUE_MEAN_NUMERATOR).to_bits() {
                return Some(path("completed_pair_center"));
            }
            if values[0] != values[4] || values[1] != values[3] {
                return Some(path("fixed_lattice_symmetry"));
            }
            let decoded = values.map(f64::from_bits);
            if !(decoded[2] <= decoded[1]
                && decoded[1] <= decoded[0]
                && decoded[2] <= decoded[3]
                && decoded[3] <= decoded[4])
            {
                return Some(path("fixed_lattice_monotonicity"));
            }
        }
    }
    None
}

fn exact_radius_bit_delta(
    reference: &SealedStudy,
    mutant: &SealedStudy,
    mutation: Mutation,
) -> bool {
    if mutation.before ^ mutation.after != 1u64 << mutation.mantissa_bit {
        return false;
    }
    let mut expected = reference.record.clone();
    let Some(checkpoint) = expected.checkpoints.get_mut(mutation.checkpoint) else {
        return false;
    };
    if checkpoint.radius_bits != mutation.before {
        return false;
    }
    checkpoint.radius_bits = mutation.after;
    expected == mutant.record
}

fn seeded_corruption(config: &ReplayIdentity, reference: &SealedStudy) -> Corruption {
    let mut selector = StreamKey {
        seed: MUTATION_SEED,
        kernel: MUTATION_KERNEL,
        tile: MUTATION_TILE,
    }
    .stream();
    let selector_word = selector.next_u64();
    let bit_selector_word = selector.next_u64();
    let draws = selector.index();
    assert_eq!(
        draws, MUTATION_DRAWS,
        "checkpoint and bit selection must consume exactly two draws"
    );
    let checkpoint = usize::try_from(
        selector_word % u64::try_from(reference.record.checkpoints.len()).expect("short fixture"),
    )
    .expect("checkpoint index fits usize");
    // Stay in the low significand: the mutated radius remains finite,
    // positive, and close enough that only the retained semantic gate—not a
    // coarse finiteness check—can catch it.
    let mantissa_bit = 12 + u32::try_from(bit_selector_word % 12).expect("bit index fits u32");
    let before = reference.record.checkpoints[checkpoint].radius_bits;
    let after = before ^ (1u64 << mantissa_bit);
    assert!(f64::from_bits(after).is_finite() && f64::from_bits(after) > 0.0);

    let mutation = Mutation {
        seed: MUTATION_SEED,
        kernel: MUTATION_KERNEL,
        tile: MUTATION_TILE,
        selector_word,
        bit_selector_word,
        draws,
        checkpoint,
        mantissa_bit,
        before,
        after,
    };
    let mut stale = reference.clone();
    stale.record.checkpoints[checkpoint].radius_bits = after;
    let stale_error = stale
        .validate_payload(config)
        .expect_err("unsealed checkpoint mutation must fail payload identity");
    let mutant = SealedStudy::seal(config, stale.record);
    let reference_error = mutant
        .admit_against(config, &reference.identity)
        .expect_err("resealed checkpoint mutation must miss retained reference");
    let semantic_mismatch = semantic_mismatch(&mutant.record)
        .expect("resealed checkpoint mutation must fail the independent oracle");
    Corruption {
        mutant,
        mutation,
        stale_error,
        reference_error,
        semantic_mismatch,
    }
}

fn red_event(reference: &SealedStudy, corruption: &Corruption) -> Event {
    let mutation = corruption.mutation;
    let detail = format!(
        "input_seed=0x{INPUT_SEED:016x}; corruption_seed=0x{:016x}; kernel=0x{:04x}; \
         tile={}; selector_word=0x{:016x}; bit_selector_word=0x{:016x}; draws={}; \
         target=checkpoints[{}].radius; mantissa_bit={}; before=0x{:016x}; \
         after=0x{:016x}; reference=0x{:016x}; mutant=0x{:016x}; stale_gate={:?}; \
         reference_gate={:?}; semantic_gate={}",
        mutation.seed,
        mutation.kernel,
        mutation.tile,
        mutation.selector_word,
        mutation.bit_selector_word,
        mutation.draws,
        mutation.checkpoint,
        mutation.mantissa_bit,
        mutation.before,
        mutation.after,
        reference.identity.root(),
        corruption.mutant.identity.root(),
        corruption.stale_error,
        corruption.reference_error,
        corruption.semantic_mismatch,
    );
    let mut emitter = Emitter::new(SUITE, RED_CASE);
    let event = emitter.emit(
        Severity::Error,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: RED_CASE.to_string(),
            pass: false,
            detail,
            seed: MUTATION_SEED,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("seeded red record must be replayable");
    fs_obs::validate_line(&event.to_jsonl()).expect("seeded red record must be wire-valid");
    event
}

fn assert_mergeable(event: &Event) {
    if let EventKind::ConformanceCase {
        case,
        pass: false,
        detail,
        ..
    } = &event.kind
    {
        panic!("merge refused by {case}: {detail}");
    }
}

fn panic_message(payload: &(dyn core::any::Any + Send)) -> String {
    payload
        .downcast_ref::<String>()
        .cloned()
        .or_else(|| {
            payload
                .downcast_ref::<&str>()
                .map(|text| (*text).to_string())
        })
        .unwrap_or_else(|| "non-string panic".to_string())
}

fn emit_case(emitter: &mut Emitter, case: &str, detail: String) {
    let event = emitter.emit(
        Severity::Info,
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: case.to_string(),
            pass: true,
            detail,
            seed: INPUT_SEED,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("green verdict must satisfy replay lint");
    fs_obs::validate_line(&event.to_jsonl()).expect("green verdict must be wire-valid");
    let line = event.to_jsonl();
    println!("{line}");
}

fn emit_receipt(
    config: &ReplayIdentity,
    result: &ReplayIdentity,
    record: &StudyRecord,
    corruption: &Corruption,
    red: &Event,
) {
    let mut emitter = Emitter::new(SUITE, CASE);
    let receipt = emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "gaussian-mixture-cs-full-trajectory-replay".to_string(),
            json: format!(
                "{{\"config_identity\":\"{}\",\"result_identity\":\"{}\",\
                 \"input_seed\":{},\"input_kernel\":{},\"input_tile\":{},\
                 \"input_draws\":{},\"stream_semantics_version\":{},\
                 \"rotation\":{},\"orientation_word\":\"0x{:016x}\",\
                 \"observations\":{},\"checkpoints\":{},\
                 \"mutation_seed\":{},\"mutation_checkpoint\":{},\
                 \"mutation_bit\":{},\"mutation_before\":\"0x{:016x}\",\
                 \"mutation_after\":\"0x{:016x}\",\"mutant_identity\":\"{}\",\
                 \"red_event_identity\":\"0x{:016x}\",\
                 \"scope\":\"finite same-process complete-public-surface fixture\"}}",
                config.hex(),
                result.hex(),
                INPUT_SEED,
                INPUT_KERNEL,
                INPUT_TILE,
                INPUT_DRAWS,
                fs_rand::STREAM_SEMANTICS_VERSION,
                record.rotation,
                record.orientation_word,
                record.observations.len(),
                record.checkpoints.len(),
                corruption.mutation.seed,
                corruption.mutation.checkpoint,
                corruption.mutation.mantissa_bit,
                corruption.mutation.before,
                corruption.mutation.after,
                corruption.mutant.identity.hex(),
                red.content_hash(),
            ),
        },
        None,
    );
    fs_obs::validate_line(&receipt.to_jsonl()).expect("replay receipt must be wire-valid");
    let line = receipt.to_jsonl();
    println!("{line}");
    emit_case(
        &mut emitter,
        "integer-oracle-and-mixture-duality",
        format!(
            "{} checkpoints match the integer/log-domain oracle; radius/e-values <= {ORACLE_ULPS} ulp; boundary relative tolerance={BOUNDARY_REL_TOLERANCE:e}",
            record.checkpoints.len(),
        ),
    );
    emit_case(
        &mut emitter,
        "full-public-trajectory-replay",
        format!(
            "config={}; result={}; observations={}; checkpoints={}",
            config.hex(),
            result.hex(),
            record.observations.len(),
            record.checkpoints.len(),
        ),
    );
    emit_case(
        &mut emitter,
        "seeded-corruption-refused",
        format!(
            "seed=0x{:016x}; selector_word=0x{:016x}; bit_selector_word=0x{:016x}; checkpoint={}; mantissa_bit={}; before=0x{:016x}; after=0x{:016x}; stale={:?}; reference={:?}; semantic={}",
            corruption.mutation.seed,
            corruption.mutation.selector_word,
            corruption.mutation.bit_selector_word,
            corruption.mutation.checkpoint,
            corruption.mutation.mantissa_bit,
            corruption.mutation.before,
            corruption.mutation.after,
            corruption.stale_error,
            corruption.reference_error,
            corruption.semantic_mismatch,
        ),
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn gaussian_mixture_cs_full_trajectory_replays_and_seeded_failure_is_refused() {
    let config = config_identity();
    let replayed_config = config_identity();
    assert_eq!(
        config, replayed_config,
        "configuration identity must replay"
    );

    let original = run_study(&config);
    let replay = run_study(&replayed_config);
    assert_eq!(
        semantic_mismatch(&original.record),
        None,
        "original trajectory must satisfy the independent oracle"
    );
    assert_eq!(
        semantic_mismatch(&replay.record),
        None,
        "replayed trajectory must satisfy the independent oracle"
    );
    assert_eq!(original, replay, "complete public trajectory must replay");
    original
        .validate_payload(&config)
        .expect("original payload identity must validate");
    replay
        .admit_against(&replayed_config, &original.identity)
        .expect("replayed payload must match the retained reference");

    let first = seeded_corruption(&config, &original);
    let second = seeded_corruption(&replayed_config, &replay);
    assert_eq!(first, second, "seeded corruption must replay exactly");
    assert!(exact_radius_bit_delta(
        &original,
        &first.mutant,
        first.mutation
    ));
    assert_eq!(
        first.semantic_mismatch,
        format!("checkpoints[{}].radius", first.mutation.checkpoint)
    );
    assert!(matches!(
        first.stale_error,
        AdmissionError::PayloadIdentityMismatch { .. }
    ));
    assert!(matches!(
        first.reference_error,
        AdmissionError::ReferenceIdentityMismatch { .. }
    ));
    assert_eq!(
        first.stale_error,
        AdmissionError::PayloadIdentityMismatch {
            declared: original.identity.root(),
            computed: first.mutant.identity.root(),
        }
    );
    assert_eq!(
        first.reference_error,
        AdmissionError::ReferenceIdentityMismatch {
            expected: original.identity.root(),
            found: first.mutant.identity.root(),
        }
    );
    first
        .mutant
        .validate_payload(&config)
        .expect("resealed mutant payload must be internally self-consistent");

    let first_red = red_event(&original, &first);
    let second_red = red_event(&replay, &second);
    assert_eq!(
        first_red.to_jsonl(),
        second_red.to_jsonl(),
        "seeded red evidence must be byte-stable"
    );
    assert_eq!(first_red.content_identity(), second_red.content_identity());
    let merge_panic = catch_unwind(|| assert_mergeable(&first_red))
        .expect_err("test-local merge gate must refuse the seeded corruption");
    let merge_message = panic_message(merge_panic.as_ref());
    assert!(merge_message.contains(RED_CASE));
    assert!(merge_message.contains("PayloadIdentityMismatch"));
    assert!(merge_message.contains("ReferenceIdentityMismatch"));
    assert!(merge_message.contains(&first.semantic_mismatch));

    emit_receipt(
        &config,
        &original.identity,
        &original.record,
        &first,
        &first_red,
    );
}
