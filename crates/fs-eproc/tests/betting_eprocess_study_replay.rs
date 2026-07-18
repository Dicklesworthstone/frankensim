//! G0/G3/G5 full-trajectory replay for [`fs_eproc::BettingEProcess`].
//!
//! The production process consumes one fixed bounded trajectory while a
//! test-local shadow independently reconstructs reference-model predictable
//! bets, sufficient statistics, strict-log wealth, and a fixed alpha lattice.
//! Every pre-observation production state is also cloned and fed a non-null
//! probe, making the otherwise wealth-neutral observation at the null unable
//! to conceal a production-bet change. Shadow-only fields are annotations, not
//! claims of production private-state introspection. The fixture exercises
//! clone/resume and verifies that refused observations leave cloned public
//! state untouched. A separately
//! seeded one-bit checkpoint corruption is retained as stable red evidence.
//! Its `StreamKey` is evidence-generator provenance only; it does not randomize
//! the statistical fixture or strengthen any validity claim.
//!
//! This finite deterministic fixture does not claim general Ville validity,
//! type-I error or power, all-null or stopping-time laws, PairwiseRace or
//! optimizer equivalence, arbitrary horizons, private-state serialization,
//! cross-ISA equality, cancellation/Cx behavior, authenticated persistence,
//! or performance.

#![deny(unsafe_code)]

use fs_eproc::BettingEProcess;
use fs_obs::ident::{IdentityBuilder, ReplayIdentity};
use fs_obs::{Emitter, Event, EventKind, Severity};
use fs_rand::StreamKey;
use std::panic::{AssertUnwindSafe, catch_unwind};

const SUITE: &str = "fs-eproc/betting-eprocess-study-replay-v2";
const CASE: &str = "predictable-bet-full-trajectory";
const RED_CASE: &str = "seeded-checkpoint-log-e-corruption";

const NULL_MEAN: f64 = 0.5;
const AGGRESSIVENESS: f64 = 0.5;
const VARIANCE_FLOOR: f64 = 1.0e-4;
const PROBE_OBSERVATION: f64 = 1.0;
const PREFIX: [f64; 7] = [1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0];
const ONES: usize = 24;
const SUFFIX: [f64; 3] = [0.5, 0.0, 1.0];
const ALPHAS: [f64; 10] = [
    0.5,
    0.25,
    0.125,
    0.0625,
    0.03125,
    0.015625,
    0.0078125,
    0.00390625,
    0.001953125,
    0.0009765625,
];

const MUTATION_SEED: u64 = 0xB17E_2142_2026_0718;
const MUTATION_KERNEL: u32 = 0xE521;
const MUTATION_TILE: u32 = 0;
const MUTATION_DRAWS: u64 = 2;
const MUTATION_CANDIDATES: [usize; 10] = [1, 2, 3, 7, 15, 23, 30, 31, 32, 33];
const LOW_MANTISSA_START: u32 = 8;
const LOW_MANTISSA_WIDTH: u32 = 16;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BetClass {
    ZeroClamped,
    Interior,
    CapClamped,
}

impl BetClass {
    const fn label(self) -> &'static str {
        match self {
            Self::ZeroClamped => "zero-clamped",
            Self::Interior => "interior",
            Self::CapClamped => "cap-clamped",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct InitialPublicState {
    len: u64,
    is_empty: bool,
    log_e_bits: u64,
    e_value_bits: u64,
    rejection_decisions: [bool; ALPHAS.len()],
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Checkpoint {
    ordinal: u64,
    observation_bits: u64,
    /// Independent reference-model annotations, never production introspection.
    shadow: ShadowStep,
    /// Public outputs from the production trajectory observation.
    production: PublicStep,
    /// Public outputs from a clone of the pre-observation production state
    /// after feeding [`PROBE_OBSERVATION`]. This makes its next bet observable.
    production_probe: PublicStep,
    first_crossing_ordinals: [Option<u64>; ALPHAS.len()],
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PublicStep {
    returned_log_e_bits: u64,
    accessor_log_e_bits: u64,
    e_value_bits: u64,
    len: u64,
    is_empty: bool,
    rejection_decisions: [bool; ALPHAS.len()],
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RetainedGolden {
    config_root: u64,
    result_root: u64,
    mutant_root: u64,
    red_event_content_hash: u64,
    first_crossing_ordinals: [Option<u64>; ALPHAS.len()],
}

// Re-pinned from the first centralized v2 runtime after this static repair.
// The all-sentinel value deliberately keeps the code-first snapshot from being
// mistaken for retained semantic proof before that batch run occurs.
const EXPECTED_RETAINED_GOLDEN: RetainedGolden = RetainedGolden {
    config_root: u64::MAX,
    result_root: u64::MAX,
    mutant_root: u64::MAX,
    red_event_content_hash: u64::MAX,
    first_crossing_ordinals: [None; ALPHAS.len()],
};

#[derive(Debug, Clone, PartialEq, Eq)]
struct StudyRecord {
    observations: Vec<u64>,
    initial: InitialPublicState,
    checkpoints: Vec<Checkpoint>,
    first_crossing_ordinals: [Option<u64>; ALPHAS.len()],
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
struct MutationCoordinates {
    selector_word: u64,
    bit_selector_word: u64,
    draws: u64,
    checkpoint: usize,
    mantissa_bit: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Mutation {
    coordinates: MutationCoordinates,
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

#[derive(Debug, Clone)]
struct ShadowEProcess {
    null_mean: f64,
    aggressiveness: f64,
    n: u64,
    sum: f64,
    sum_sq: f64,
    log_wealth: f64,
    observation_bits: Vec<u64>,
    first_crossing_ordinals: [Option<u64>; ALPHAS.len()],
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct ShadowSnapshot {
    null_mean_bits: u64,
    aggressiveness_bits: u64,
    n: u64,
    sum_bits: u64,
    sum_sq_bits: u64,
    log_wealth_bits: u64,
    observation_bits: Vec<u64>,
    first_crossing_ordinals: [Option<u64>; ALPHAS.len()],
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ShadowStep {
    n_before: u64,
    sum_before_bits: u64,
    sum_sq_before_bits: u64,
    regularized_mean_bits: u64,
    regularized_second_moment_bits: u64,
    unfloored_variance_bits: u64,
    variance_bits: u64,
    raw_lambda_bits: u64,
    lambda_cap_bits: u64,
    clipped_lambda_bits: u64,
    bet_class: BetClass,
    wealth_factor_bits: u64,
    log_increment_bits: u64,
    log_wealth_before_bits: u64,
    n_after: u64,
    sum_after_bits: u64,
    sum_sq_after_bits: u64,
    log_wealth_after_bits: u64,
    e_value_bits: u64,
    rejection_decisions: [bool; ALPHAS.len()],
    first_crossing_ordinals: [Option<u64>; ALPHAS.len()],
}

impl ShadowEProcess {
    fn new(null_mean: f64) -> Self {
        Self {
            null_mean,
            aggressiveness: AGGRESSIVENESS,
            n: 0,
            sum: 0.0,
            sum_sq: 0.0,
            log_wealth: 0.0,
            observation_bits: Vec::new(),
            first_crossing_ordinals: [None; ALPHAS.len()],
        }
    }

    fn snapshot(&self) -> ShadowSnapshot {
        ShadowSnapshot {
            null_mean_bits: self.null_mean.to_bits(),
            aggressiveness_bits: self.aggressiveness.to_bits(),
            n: self.n,
            sum_bits: self.sum.to_bits(),
            sum_sq_bits: self.sum_sq.to_bits(),
            log_wealth_bits: self.log_wealth.to_bits(),
            observation_bits: self.observation_bits.clone(),
            first_crossing_ordinals: self.first_crossing_ordinals,
        }
    }

    fn observe(&mut self, observation: f64) -> ShadowStep {
        assert!((0.0..=1.0).contains(&observation));
        let n_before = self.n;
        let sum_before = self.sum;
        let sum_sq_before = self.sum_sq;
        let log_wealth_before = self.log_wealth;
        let n = n_before as f64;
        let regularized_mean = (sum_before + self.null_mean) / (n + 1.0);
        let regularized_second_moment =
            (sum_sq_before + self.null_mean * self.null_mean) / (n + 1.0);
        let unfloored_variance = regularized_second_moment - regularized_mean * regularized_mean;
        let variance = unfloored_variance.max(VARIANCE_FLOOR);
        let raw_lambda = (regularized_mean - self.null_mean) / variance;
        let lambda_cap = self.aggressiveness / self.null_mean.max(1.0 - self.null_mean);
        let clipped_lambda = raw_lambda.clamp(0.0, lambda_cap);
        let bet_class = if raw_lambda <= 0.0 {
            BetClass::ZeroClamped
        } else if raw_lambda >= lambda_cap {
            BetClass::CapClamped
        } else {
            BetClass::Interior
        };
        let wealth_factor = clipped_lambda.mul_add(observation - self.null_mean, 1.0);
        let log_increment = fs_math::det::ln(wealth_factor);
        self.log_wealth += log_increment;
        self.n += 1;
        self.sum += observation;
        self.sum_sq += observation * observation;
        self.observation_bits.push(observation.to_bits());

        let rejection_decisions = ALPHAS.map(|alpha| self.log_wealth >= -fs_math::det::ln(alpha));
        for (index, &rejects) in rejection_decisions.iter().enumerate() {
            if rejects && self.first_crossing_ordinals[index].is_none() {
                self.first_crossing_ordinals[index] = Some(self.n);
            }
        }

        ShadowStep {
            n_before,
            sum_before_bits: sum_before.to_bits(),
            sum_sq_before_bits: sum_sq_before.to_bits(),
            regularized_mean_bits: regularized_mean.to_bits(),
            regularized_second_moment_bits: regularized_second_moment.to_bits(),
            unfloored_variance_bits: unfloored_variance.to_bits(),
            variance_bits: variance.to_bits(),
            raw_lambda_bits: raw_lambda.to_bits(),
            lambda_cap_bits: lambda_cap.to_bits(),
            clipped_lambda_bits: clipped_lambda.to_bits(),
            bet_class,
            wealth_factor_bits: wealth_factor.to_bits(),
            log_increment_bits: log_increment.to_bits(),
            log_wealth_before_bits: log_wealth_before.to_bits(),
            n_after: self.n,
            sum_after_bits: self.sum.to_bits(),
            sum_sq_after_bits: self.sum_sq.to_bits(),
            log_wealth_after_bits: self.log_wealth.to_bits(),
            e_value_bits: fs_math::det::exp(self.log_wealth).to_bits(),
            rejection_decisions,
            first_crossing_ordinals: self.first_crossing_ordinals,
        }
    }
}

fn usize_u64(value: usize) -> u64 {
    u64::try_from(value).expect("fixture cardinality fits u64")
}

fn observations() -> Vec<f64> {
    let mut values = Vec::with_capacity(PREFIX.len() + ONES + SUFFIX.len());
    values.extend(PREFIX);
    values.extend(std::iter::repeat_n(1.0, ONES));
    values.extend(SUFFIX);
    values
}

fn public_state(process: &BettingEProcess) -> InitialPublicState {
    InitialPublicState {
        len: process.len(),
        is_empty: process.is_empty(),
        log_e_bits: process.log_e_value().to_bits(),
        e_value_bits: process.e_value().to_bits(),
        rejection_decisions: ALPHAS.map(|alpha| process.rejects_at(alpha)),
    }
}

fn mutation_coordinates() -> MutationCoordinates {
    let mut selector = StreamKey {
        seed: MUTATION_SEED,
        kernel: MUTATION_KERNEL,
        tile: MUTATION_TILE,
    }
    .stream();
    let selector_word = selector.next_u64();
    let bit_selector_word = selector.next_u64();
    let draws = selector.index();
    assert_eq!(draws, MUTATION_DRAWS);
    let candidate = usize::try_from(selector_word % usize_u64(MUTATION_CANDIDATES.len()))
        .expect("candidate index fits usize");
    let mantissa_bit = LOW_MANTISSA_START
        + u32::try_from(bit_selector_word % u64::from(LOW_MANTISSA_WIDTH))
            .expect("mantissa selector fits u32");
    MutationCoordinates {
        selector_word,
        bit_selector_word,
        draws,
        checkpoint: MUTATION_CANDIDATES[candidate],
        mantissa_bit,
    }
}

fn config_identity(coordinates: MutationCoordinates) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-eproc-betting-eprocess-config-v2")
        .str("statistical-object", "one-sided-predictable-betting-e-process")
        .str("observation-units", "dimensionless-bounded-score")
        .f64_bits("support-lower", 0.0)
        .f64_bits("support-upper", 1.0)
        .f64_bits("null-mean", NULL_MEAN)
        .f64_bits("aggressiveness", AGGRESSIVENESS)
        .f64_bits("variance-floor", VARIANCE_FLOOR)
        .f64_bits("production-probe-observation", PROBE_OBSERVATION)
        .str("regularization", "one-pseudo-observation-at-null")
        .str("lambda-policy", "raw=(regularized-mean-null)/floored-variance")
        .str("lambda-clamp", "zero-to-aggressiveness-over-max-null-complement")
        .u64("prefix-length", usize_u64(PREFIX.len()))
        .u64("one-run-length", usize_u64(ONES))
        .u64("suffix-length", usize_u64(SUFFIX.len()))
        .u64(
            "stream-semantics-version",
            u64::from(fs_rand::STREAM_SEMANTICS_VERSION),
        )
        .u64("mutation-seed", MUTATION_SEED)
        .u64("mutation-kernel", u64::from(MUTATION_KERNEL))
        .u64("mutation-tile", u64::from(MUTATION_TILE))
        .u64("mutation-draws", coordinates.draws)
        .u64("mutation-selector-word", coordinates.selector_word)
        .u64("mutation-bit-selector-word", coordinates.bit_selector_word)
        .u64("mutation-checkpoint-index", usize_u64(coordinates.checkpoint))
        .u64("mutation-mantissa-bit", u64::from(coordinates.mantissa_bit))
        .u64("low-mantissa-start", u64::from(LOW_MANTISSA_START))
        .u64("low-mantissa-width", u64::from(LOW_MANTISSA_WIDTH))
        .str("mutation-rng-role", "evidence-generator-provenance-only")
        .str("fs-eproc-version", fs_eproc::VERSION)
        .str("fs-math-version", fs_math::VERSION)
        .str("fs-obs-version", fs_obs::VERSION)
        .str("fs-rand-version", fs_rand::VERSION)
        .str(
            "no-claim",
            "finite-fixture-only-no-general-Ville-type-I-power-all-null-stopping-law-PairwiseRace-optimizer-arbitrary-horizon-private-state-cross-ISA-Cx-cancellation-authenticated-persistence-performance-claim",
        );
    for value in observations() {
        builder = builder.f64_bits("ordered-observation", value);
    }
    for alpha in ALPHAS {
        builder = builder.f64_bits("alpha-lattice", alpha);
    }
    for candidate in MUTATION_CANDIDATES {
        builder = builder.u64("mutation-candidate-index", usize_u64(candidate));
    }
    builder.finish()
}

fn result_identity(config: &ReplayIdentity, record: &StudyRecord) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-eproc-betting-eprocess-result-v2")
        .child("config", config)
        .u64("observation-count", usize_u64(record.observations.len()))
        .u64("checkpoint-count", usize_u64(record.checkpoints.len()))
        .u64("alpha-count", usize_u64(ALPHAS.len()))
        .u64("initial-len", record.initial.len)
        .flag("initial-is-empty", record.initial.is_empty)
        .u64("initial-log-e-bits", record.initial.log_e_bits)
        .u64("initial-e-value-bits", record.initial.e_value_bits);
    for decision in record.initial.rejection_decisions {
        builder = builder.flag("initial-rejection-decision", decision);
    }
    for observation_bits in &record.observations {
        builder = builder.u64("ordered-observation-bits", *observation_bits);
    }
    for checkpoint in &record.checkpoints {
        builder = builder
            .u64("checkpoint-ordinal", checkpoint.ordinal)
            .u64("checkpoint-observation-bits", checkpoint.observation_bits)
            .u64("shadow-n-before", checkpoint.shadow.n_before)
            .u64("shadow-sum-before-bits", checkpoint.shadow.sum_before_bits)
            .u64(
                "shadow-sum-sq-before-bits",
                checkpoint.shadow.sum_sq_before_bits,
            )
            .u64(
                "shadow-regularized-mean-bits",
                checkpoint.shadow.regularized_mean_bits,
            )
            .u64(
                "shadow-regularized-second-moment-bits",
                checkpoint.shadow.regularized_second_moment_bits,
            )
            .u64(
                "shadow-unfloored-variance-bits",
                checkpoint.shadow.unfloored_variance_bits,
            )
            .u64("shadow-variance-bits", checkpoint.shadow.variance_bits)
            .u64("shadow-raw-lambda-bits", checkpoint.shadow.raw_lambda_bits)
            .u64("shadow-lambda-cap-bits", checkpoint.shadow.lambda_cap_bits)
            .u64(
                "shadow-clipped-lambda-bits",
                checkpoint.shadow.clipped_lambda_bits,
            )
            .str("shadow-bet-class", checkpoint.shadow.bet_class.label())
            .u64(
                "shadow-wealth-factor-bits",
                checkpoint.shadow.wealth_factor_bits,
            )
            .u64(
                "shadow-log-increment-bits",
                checkpoint.shadow.log_increment_bits,
            )
            .u64(
                "shadow-log-wealth-before-bits",
                checkpoint.shadow.log_wealth_before_bits,
            )
            .u64("shadow-n-after", checkpoint.shadow.n_after)
            .u64("shadow-sum-after-bits", checkpoint.shadow.sum_after_bits)
            .u64(
                "shadow-sum-sq-after-bits",
                checkpoint.shadow.sum_sq_after_bits,
            )
            .u64(
                "shadow-log-wealth-after-bits",
                checkpoint.shadow.log_wealth_after_bits,
            )
            .u64("shadow-e-value-bits", checkpoint.shadow.e_value_bits)
            .u64(
                "production-returned-log-e-bits",
                checkpoint.production.returned_log_e_bits,
            )
            .u64(
                "production-accessor-log-e-bits",
                checkpoint.production.accessor_log_e_bits,
            )
            .u64(
                "production-e-value-bits",
                checkpoint.production.e_value_bits,
            )
            .u64("production-len", checkpoint.production.len)
            .flag("production-is-empty", checkpoint.production.is_empty)
            .u64(
                "production-probe-returned-log-e-bits",
                checkpoint.production_probe.returned_log_e_bits,
            )
            .u64(
                "production-probe-accessor-log-e-bits",
                checkpoint.production_probe.accessor_log_e_bits,
            )
            .u64(
                "production-probe-e-value-bits",
                checkpoint.production_probe.e_value_bits,
            )
            .u64("production-probe-len", checkpoint.production_probe.len)
            .flag(
                "production-probe-is-empty",
                checkpoint.production_probe.is_empty,
            );
        for decision in checkpoint.shadow.rejection_decisions {
            builder = builder.flag("shadow-rejection-decision", decision);
        }
        for crossing in checkpoint.shadow.first_crossing_ordinals {
            builder = builder
                .flag("shadow-first-crossing-present", crossing.is_some())
                .u64("shadow-first-crossing-ordinal", crossing.unwrap_or(0));
        }
        for decision in checkpoint.production.rejection_decisions {
            builder = builder.flag("production-rejection-decision", decision);
        }
        for decision in checkpoint.production_probe.rejection_decisions {
            builder = builder.flag("production-probe-rejection-decision", decision);
        }
        for crossing in checkpoint.first_crossing_ordinals {
            builder = builder
                .flag("checkpoint-first-crossing-present", crossing.is_some())
                .u64("checkpoint-first-crossing-ordinal", crossing.unwrap_or(0));
        }
    }
    for crossing in record.first_crossing_ordinals {
        builder = builder
            .flag("result-first-crossing-present", crossing.is_some())
            .u64("result-first-crossing-ordinal", crossing.unwrap_or(0));
    }
    builder.finish()
}

impl SealedStudy {
    fn seal(config: &ReplayIdentity, record: StudyRecord) -> Self {
        Self {
            identity: result_identity(config, &record),
            record,
        }
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

fn public_step(process: &mut BettingEProcess, observation: f64) -> PublicStep {
    let returned_log_e = process.observe(observation);
    PublicStep {
        returned_log_e_bits: returned_log_e.to_bits(),
        accessor_log_e_bits: process.log_e_value().to_bits(),
        e_value_bits: process.e_value().to_bits(),
        len: process.len(),
        is_empty: process.is_empty(),
        rejection_decisions: ALPHAS.map(|alpha| process.rejects_at(alpha)),
    }
}

fn shadow_public_step(step: ShadowStep) -> PublicStep {
    PublicStep {
        returned_log_e_bits: step.log_wealth_after_bits,
        accessor_log_e_bits: step.log_wealth_after_bits,
        e_value_bits: step.e_value_bits,
        len: step.n_after,
        is_empty: false,
        rejection_decisions: step.rejection_decisions,
    }
}

fn observe_checkpoint(
    process: &mut BettingEProcess,
    shadow: &mut ShadowEProcess,
    ordinal: usize,
    observation: f64,
    first_crossings: &mut [Option<u64>; ALPHAS.len()],
) -> Checkpoint {
    let mut shadow_probe = shadow.clone();
    let expected_probe = shadow_probe.observe(PROBE_OBSERVATION);
    let mut process_probe = process.clone();
    let production_probe = public_step(&mut process_probe, PROBE_OBSERVATION);
    assert_eq!(
        production_probe,
        shadow_public_step(expected_probe),
        "non-null clone probe must expose the production next bet",
    );

    let expected = shadow.observe(observation);
    let production = public_step(process, observation);
    assert_eq!(
        production,
        shadow_public_step(expected),
        "production public trajectory must match the shadow reference",
    );
    let ordinal_u64 = usize_u64(ordinal);
    for (index, &rejects) in production.rejection_decisions.iter().enumerate() {
        if rejects && first_crossings[index].is_none() {
            first_crossings[index] = Some(ordinal_u64);
        }
    }

    assert_eq!(*first_crossings, expected.first_crossing_ordinals);

    Checkpoint {
        ordinal: ordinal_u64,
        observation_bits: observation.to_bits(),
        shadow: expected,
        production,
        production_probe,
        first_crossing_ordinals: *first_crossings,
    }
}

fn run_study(config: &ReplayIdentity) -> SealedStudy {
    let values = observations();
    let mut process = BettingEProcess::new(NULL_MEAN);
    let mut shadow = ShadowEProcess::new(NULL_MEAN);
    let initial = public_state(&process);
    assert_eq!(initial.len, 0);
    assert!(initial.is_empty);
    assert_eq!(initial.log_e_bits, 0.0f64.to_bits());
    assert_eq!(initial.e_value_bits, 1.0f64.to_bits());
    assert_eq!(initial.rejection_decisions, [false; ALPHAS.len()]);

    let mut first_crossings = [None; ALPHAS.len()];
    let checkpoints = values
        .iter()
        .enumerate()
        .map(|(index, &observation)| {
            observe_checkpoint(
                &mut process,
                &mut shadow,
                index + 1,
                observation,
                &mut first_crossings,
            )
        })
        .collect();
    assert_eq!(
        shadow.observation_bits,
        values.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
    );

    SealedStudy::seal(
        config,
        StudyRecord {
            observations: values.iter().map(|value| value.to_bits()).collect(),
            initial,
            checkpoints,
            first_crossing_ordinals: first_crossings,
        },
    )
}

fn checkpoint_mismatch(index: usize, actual: &Checkpoint, expected: &Checkpoint) -> Option<String> {
    macro_rules! check_top {
        ($field:ident) => {
            if actual.$field != expected.$field {
                return Some(format!("checkpoints[{index}].{}", stringify!($field)));
            }
        };
    }
    macro_rules! check_shadow {
        ($field:ident) => {
            if actual.shadow.$field != expected.shadow.$field {
                return Some(format!(
                    "checkpoints[{index}].shadow.{}",
                    stringify!($field),
                ));
            }
        };
    }
    macro_rules! check_public {
        ($which:ident, $field:ident) => {
            if actual.$which.$field != expected.$which.$field {
                return Some(format!(
                    "checkpoints[{index}].{}.{}",
                    stringify!($which),
                    stringify!($field),
                ));
            }
        };
    }
    check_top!(ordinal);
    check_top!(observation_bits);
    check_shadow!(n_before);
    check_shadow!(sum_before_bits);
    check_shadow!(sum_sq_before_bits);
    check_shadow!(regularized_mean_bits);
    check_shadow!(regularized_second_moment_bits);
    check_shadow!(unfloored_variance_bits);
    check_shadow!(variance_bits);
    check_shadow!(raw_lambda_bits);
    check_shadow!(lambda_cap_bits);
    check_shadow!(clipped_lambda_bits);
    check_shadow!(bet_class);
    check_shadow!(wealth_factor_bits);
    check_shadow!(log_increment_bits);
    check_shadow!(log_wealth_before_bits);
    check_shadow!(n_after);
    check_shadow!(sum_after_bits);
    check_shadow!(sum_sq_after_bits);
    check_shadow!(log_wealth_after_bits);
    check_shadow!(e_value_bits);
    check_shadow!(rejection_decisions);
    check_shadow!(first_crossing_ordinals);
    check_public!(production, returned_log_e_bits);
    check_public!(production, accessor_log_e_bits);
    check_public!(production, e_value_bits);
    check_public!(production, len);
    check_public!(production, is_empty);
    check_public!(production, rejection_decisions);
    check_public!(production_probe, returned_log_e_bits);
    check_public!(production_probe, accessor_log_e_bits);
    check_public!(production_probe, e_value_bits);
    check_public!(production_probe, len);
    check_public!(production_probe, is_empty);
    check_public!(production_probe, rejection_decisions);
    check_top!(first_crossing_ordinals);
    None
}

#[allow(clippy::too_many_lines)]
fn semantic_mismatch(record: &StudyRecord) -> Option<String> {
    let expected_observations = observations();
    let expected_bits: Vec<_> = expected_observations
        .iter()
        .map(|value| value.to_bits())
        .collect();
    if record.observations != expected_bits {
        return Some("observations".to_string());
    }
    let expected_initial = InitialPublicState {
        len: 0,
        is_empty: true,
        log_e_bits: 0.0f64.to_bits(),
        e_value_bits: 1.0f64.to_bits(),
        rejection_decisions: [false; ALPHAS.len()],
    };
    if record.initial != expected_initial {
        return Some("initial_public_state".to_string());
    }
    if record.checkpoints.len() != expected_observations.len() {
        return Some("checkpoint_count".to_string());
    }

    let mut shadow = ShadowEProcess::new(NULL_MEAN);
    for (index, (&observation, actual)) in expected_observations
        .iter()
        .zip(&record.checkpoints)
        .enumerate()
    {
        let mut shadow_probe = shadow.clone();
        let expected_probe = shadow_probe.observe(PROBE_OBSERVATION);
        let expected_step = shadow.observe(observation);
        let expected = Checkpoint {
            ordinal: usize_u64(index + 1),
            observation_bits: observation.to_bits(),
            shadow: expected_step,
            production: shadow_public_step(expected_step),
            production_probe: shadow_public_step(expected_probe),
            first_crossing_ordinals: expected_step.first_crossing_ordinals,
        };
        if let Some(mismatch) = checkpoint_mismatch(index, actual, &expected) {
            return Some(mismatch);
        }
    }
    if record.first_crossing_ordinals != shadow.first_crossing_ordinals {
        return Some("first_crossing_ordinals".to_string());
    }
    None
}

fn assert_clone_resume_matches(record: &StudyRecord) {
    let values = observations();
    let mut process = BettingEProcess::new(NULL_MEAN);
    let mut shadow = ShadowEProcess::new(NULL_MEAN);
    let mut first_crossings = [None; ALPHAS.len()];

    for (index, &observation) in values.iter().take(PREFIX.len()).enumerate() {
        let checkpoint = observe_checkpoint(
            &mut process,
            &mut shadow,
            index + 1,
            observation,
            &mut first_crossings,
        );
        assert_eq!(checkpoint, record.checkpoints[index]);
    }

    let mut resumed_process = process.clone();
    let mut resumed_shadow = shadow.clone();
    let mut resumed_crossings = first_crossings;
    for (index, &observation) in values.iter().enumerate().skip(PREFIX.len()) {
        let uninterrupted = observe_checkpoint(
            &mut process,
            &mut shadow,
            index + 1,
            observation,
            &mut first_crossings,
        );
        let resumed = observe_checkpoint(
            &mut resumed_process,
            &mut resumed_shadow,
            index + 1,
            observation,
            &mut resumed_crossings,
        );
        assert_eq!(uninterrupted, record.checkpoints[index]);
        assert_eq!(resumed, uninterrupted);
    }
    assert_eq!(public_state(&resumed_process), public_state(&process));
    assert_eq!(resumed_shadow.snapshot(), shadow.snapshot());
    assert_eq!(resumed_crossings, first_crossings);
}

fn assert_invalid_observations_preserve_public_state() {
    let mut process = BettingEProcess::new(NULL_MEAN);
    for observation in PREFIX {
        let _ = process.observe(observation);
    }
    let invalid = [
        f64::NAN,
        f64::from_bits(0x8000_0000_0000_0001),
        f64::from_bits(1.0f64.to_bits() + 1),
    ];
    for observation in invalid {
        let mut clone = process.clone();
        let before = public_state(&clone);
        let panic = catch_unwind(AssertUnwindSafe(|| clone.observe(observation)));
        assert!(panic.is_err(), "invalid observation must be refused");
        assert_eq!(
            public_state(&clone),
            before,
            "refused observation changed public state: 0x{:016x}",
            observation.to_bits(),
        );
    }
}

fn assert_trajectory_semantics(record: &StudyRecord) {
    let classes: Vec<_> = record
        .checkpoints
        .iter()
        .map(|checkpoint| checkpoint.shadow.bet_class)
        .collect();
    assert!(classes.contains(&BetClass::ZeroClamped));
    assert!(classes.contains(&BetClass::Interior));
    assert!(classes.contains(&BetClass::CapClamped));
    assert_eq!(
        record.checkpoints[0].shadow.bet_class,
        BetClass::ZeroClamped
    );
    assert_eq!(record.checkpoints[1].shadow.bet_class, BetClass::CapClamped);
    assert_eq!(record.checkpoints[3].shadow.bet_class, BetClass::Interior);
    assert_eq!(
        record.checkpoints[0].shadow.unfloored_variance_bits,
        0.0f64.to_bits()
    );
    assert_eq!(
        record.checkpoints[0].shadow.variance_bits,
        VARIANCE_FLOOR.to_bits()
    );
    assert_eq!(
        record.checkpoints[0].shadow.lambda_cap_bits,
        1.0f64.to_bits()
    );

    // After [1, 1, 0, 0], the regularized mean lands exactly on the null:
    // this is the exact-zero boundary, not the negative one-sided branch.
    let exact_zero = record.checkpoints[4];
    assert_eq!(exact_zero.shadow.raw_lambda_bits, 0.0f64.to_bits());
    assert_eq!(exact_zero.shadow.clipped_lambda_bits, 0.0f64.to_bits());
    assert_eq!(exact_zero.shadow.bet_class, BetClass::ZeroClamped);
    // After the fifth observation (another zero), the raw predictable bet is
    // strictly negative and the one-sided policy must clamp it to positive 0.
    let negative_raw = record.checkpoints[5];
    assert!(f64::from_bits(negative_raw.shadow.raw_lambda_bits) < 0.0);
    assert_eq!(negative_raw.shadow.clipped_lambda_bits, 0.0f64.to_bits());
    assert_eq!(negative_raw.shadow.bet_class, BetClass::ZeroClamped);

    // The post-prefix run first recovers through a genuine interior bet, then
    // reaches and stays at the configured cap as evidence accumulates.
    let interior_recovery = record.checkpoints[PREFIX.len()];
    assert_eq!(interior_recovery.shadow.bet_class, BetClass::Interior);
    assert!(f64::from_bits(interior_recovery.shadow.clipped_lambda_bits) > 0.0);
    assert!(
        f64::from_bits(interior_recovery.shadow.clipped_lambda_bits)
            < f64::from_bits(interior_recovery.shadow.lambda_cap_bits)
    );
    let cap_recovery = record.checkpoints[PREFIX.len() + 4];
    assert_eq!(cap_recovery.shadow.bet_class, BetClass::CapClamped);
    assert_eq!(
        cap_recovery.shadow.clipped_lambda_bits,
        cap_recovery.shadow.lambda_cap_bits
    );

    let neutral = PREFIX.len() + ONES;
    let before_neutral = record.checkpoints[neutral - 1];
    let after_neutral = record.checkpoints[neutral];
    let after_loss = record.checkpoints[neutral + 1];
    let after_recovery = record.checkpoints[neutral + 2];
    assert_eq!(after_neutral.observation_bits, 0.5f64.to_bits());
    assert_eq!(after_neutral.shadow.wealth_factor_bits, 1.0f64.to_bits());
    assert_eq!(after_neutral.shadow.log_increment_bits, 0.0f64.to_bits());
    assert_eq!(
        after_neutral.production.returned_log_e_bits, before_neutral.production.returned_log_e_bits,
        "the null-centered observation must be wealth-neutral",
    );
    assert_ne!(
        after_neutral.production_probe.returned_log_e_bits,
        after_neutral.production.returned_log_e_bits,
        "the non-null clone probe must make the neutral checkpoint bet observable",
    );
    assert!(
        f64::from_bits(after_loss.production.returned_log_e_bits)
            < f64::from_bits(after_neutral.production.returned_log_e_bits),
        "zero outcome must decrease cap-bet wealth",
    );
    assert!(
        f64::from_bits(after_recovery.production.returned_log_e_bits)
            > f64::from_bits(after_loss.production.returned_log_e_bits),
        "final one outcome must recover wealth after the loss",
    );

    assert!(
        record.first_crossing_ordinals.iter().all(Option::is_some),
        "the 24-one run must cross every fixed alpha threshold",
    );
    for pair in record.first_crossing_ordinals.windows(2) {
        assert!(pair[0] <= pair[1], "stricter alpha crossed earlier");
    }
    for (alpha_index, crossing) in record.first_crossing_ordinals.iter().enumerate() {
        let ordinal = usize::try_from(crossing.expect("all thresholds cross"))
            .expect("crossing ordinal fits usize");
        assert!(
            record.checkpoints[ordinal - 1]
                .production
                .rejection_decisions[alpha_index]
        );
        if ordinal > 1 {
            assert!(
                !record.checkpoints[ordinal - 2]
                    .production
                    .rejection_decisions[alpha_index]
            );
        }
    }
}

fn exact_one_checkpoint_bit_delta(
    reference: &SealedStudy,
    mutant: &SealedStudy,
    mutation: Mutation,
) -> bool {
    if mutation.before ^ mutation.after != 1u64 << mutation.coordinates.mantissa_bit {
        return false;
    }
    let mut expected = reference.record.clone();
    let Some(checkpoint) = expected
        .checkpoints
        .get_mut(mutation.coordinates.checkpoint)
    else {
        return false;
    };
    if checkpoint.production.returned_log_e_bits != mutation.before {
        return false;
    }
    checkpoint.production.returned_log_e_bits = mutation.after;
    expected == mutant.record
}

fn seeded_corruption(
    config: &ReplayIdentity,
    reference: &SealedStudy,
    coordinates: MutationCoordinates,
) -> Corruption {
    assert_eq!(coordinates, mutation_coordinates());
    let before = reference.record.checkpoints[coordinates.checkpoint]
        .production
        .returned_log_e_bits;
    let before_value = f64::from_bits(before);
    assert!(before_value.is_finite() && before_value != 0.0);
    let after = before ^ (1u64 << coordinates.mantissa_bit);
    let after_value = f64::from_bits(after);
    assert!(after_value.is_finite() && after_value != 0.0);
    let mutation = Mutation {
        coordinates,
        before,
        after,
    };

    let mut stale = reference.clone();
    stale.record.checkpoints[coordinates.checkpoint]
        .production
        .returned_log_e_bits = after;
    let stale_error = stale
        .validate_payload(config)
        .expect_err("unsealed log-e corruption must fail payload identity");
    let mutant = SealedStudy::seal(config, stale.record);
    let reference_error = mutant
        .admit_against(config, &reference.identity)
        .expect_err("resealed log-e corruption must miss retained reference");
    let semantic_mismatch = semantic_mismatch(&mutant.record)
        .expect("resealed log-e corruption must fail independent shadow replay");

    Corruption {
        mutant,
        mutation,
        stale_error,
        reference_error,
        semantic_mismatch,
    }
}

fn red_event(reference: &SealedStudy, corruption: &Corruption) -> Event {
    let coordinates = corruption.mutation.coordinates;
    let detail = format!(
        "corruption_seed=0x{MUTATION_SEED:016x}; kernel=0x{MUTATION_KERNEL:04x}; \
         tile={MUTATION_TILE}; selector_word=0x{:016x}; bit_selector_word=0x{:016x}; \
         draws={}; target=checkpoints[{}].production.returned_log_e_bits; mantissa_bit={}; \
         before=0x{:016x}; after=0x{:016x}; reference={}; mutant={}; \
         stale_gate={:?}; reference_gate={:?}; first_semantic_mismatch={}",
        coordinates.selector_word,
        coordinates.bit_selector_word,
        coordinates.draws,
        coordinates.checkpoint,
        coordinates.mantissa_bit,
        corruption.mutation.before,
        corruption.mutation.after,
        reference.identity.hex(),
        corruption.mutant.identity.hex(),
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
            seed: 0,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("green verdict must satisfy replay lint");
    fs_obs::validate_line(&event.to_jsonl()).expect("green verdict must be wire-valid");
    println!("{}", event.to_jsonl());
}

fn emit_receipt(
    config: &ReplayIdentity,
    study: &SealedStudy,
    corruption: &Corruption,
    red: &Event,
) {
    let mut emitter = Emitter::new(SUITE, CASE);
    let receipt = emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "betting-eprocess-full-trajectory-replay".to_string(),
            json: format!(
                "{{\"config_identity\":\"{}\",\"result_identity\":\"{}\",\
                 \"observations\":{},\"checkpoints\":{},\"alpha_levels\":{},\
                 \"mutation_seed\":{},\"mutation_checkpoint\":{},\
                 \"mutation_bit\":{},\"mutation_before\":\"0x{:016x}\",\
                 \"mutation_after\":\"0x{:016x}\",\"mutant_identity\":\"{}\",\
                 \"red_event_identity\":\"0x{:016x}\",\
                 \"scope\":\"finite production-public trajectory plus non-null clone probes with explicitly shadow-only internal annotations\"}}",
                config.hex(),
                study.identity.hex(),
                study.record.observations.len(),
                study.record.checkpoints.len(),
                ALPHAS.len(),
                MUTATION_SEED,
                corruption.mutation.coordinates.checkpoint,
                corruption.mutation.coordinates.mantissa_bit,
                corruption.mutation.before,
                corruption.mutation.after,
                corruption.mutant.identity.hex(),
                red.content_hash(),
            ),
        },
        None,
    );
    fs_obs::validate_line(&receipt.to_jsonl()).expect("replay receipt must be wire-valid");
    println!("{}", receipt.to_jsonl());
    emit_case(
        &mut emitter,
        "predictable-bet-shadow-oracle",
        format!(
            "{} checkpoints matched production public outputs and non-null clone probes against shadow-only reference moments, bets, strict-log wealth, and alpha crossings",
            study.record.checkpoints.len(),
        ),
    );
    emit_case(
        &mut emitter,
        "clone-resume-and-invalid-state-preservation",
        "prefix clones replayed the uninterrupted suffix; NaN and one-ULP support violations preserved all public state".to_string(),
    );
    emit_case(
        &mut emitter,
        "seeded-corruption-refused",
        format!(
            "seed=0x{MUTATION_SEED:016x}; checkpoint={}; mantissa_bit={}; stale={:?}; reference={:?}; semantic={}",
            corruption.mutation.coordinates.checkpoint,
            corruption.mutation.coordinates.mantissa_bit,
            corruption.stale_error,
            corruption.reference_error,
            corruption.semantic_mismatch,
        ),
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn betting_eprocess_full_trajectory_replays_and_seeded_failure_is_refused() {
    let coordinates = mutation_coordinates();
    let config = config_identity(coordinates);
    let replayed_config = config_identity(mutation_coordinates());
    assert_eq!(config, replayed_config, "fixture identity must replay");
    let mut draw_count_mutant = coordinates;
    draw_count_mutant.draws ^= 1;
    assert_ne!(
        config_identity(draw_count_mutant),
        config,
        "fixture identity must bind the supplied mutation selector draw count",
    );

    let original = run_study(&config);
    let replay = run_study(&replayed_config);
    assert_eq!(
        semantic_mismatch(&original.record),
        None,
        "reference trajectory must pass independent shadow replay",
    );
    assert_eq!(
        semantic_mismatch(&replay.record),
        None,
        "replayed trajectory must pass independent shadow replay",
    );
    assert_eq!(
        original, replay,
        "complete trajectory and identity must replay"
    );
    original
        .validate_payload(&config)
        .expect("reference payload identity must validate");
    replay
        .admit_against(&replayed_config, &original.identity)
        .expect("independent rerun must match retained reference");

    assert_trajectory_semantics(&original.record);
    assert_clone_resume_matches(&original.record);
    assert_invalid_observations_preserve_public_state();

    // Presence is semantic independently of the ordinal value: an impossible
    // Some(0) must not collide with a genuine not-yet-crossed None slot.
    let mut crossing_presence_mutant = original.record.clone();
    assert_eq!(
        crossing_presence_mutant.checkpoints[0].first_crossing_ordinals[0],
        None
    );
    crossing_presence_mutant.checkpoints[0].first_crossing_ordinals[0] = Some(0);
    assert_ne!(
        result_identity(&config, &crossing_presence_mutant),
        original.identity,
        "result identity must distinguish absent crossings from ordinal zero",
    );

    // The checkpoint at the null has no trajectory wealth sensitivity to its
    // bet. Its pre-state clone probe must therefore be both identity-bound and
    // independently replay-checked as production output.
    let neutral_checkpoint = PREFIX.len() + ONES;
    let mut probe_mutant = original.record.clone();
    probe_mutant.checkpoints[neutral_checkpoint]
        .production_probe
        .returned_log_e_bits ^= 1;
    assert_ne!(
        result_identity(&config, &probe_mutant),
        original.identity,
        "production probe output must participate in result identity",
    );
    assert_eq!(
        semantic_mismatch(&probe_mutant),
        Some(format!(
            "checkpoints[{neutral_checkpoint}].production_probe.returned_log_e_bits"
        )),
        "shadow replay must localize production-probe drift at the neutral checkpoint",
    );

    let first = seeded_corruption(&config, &original, coordinates);
    let second = seeded_corruption(&replayed_config, &replay, mutation_coordinates());
    assert_eq!(first, second, "seeded corruption must replay exactly");
    assert!(exact_one_checkpoint_bit_delta(
        &original,
        &first.mutant,
        first.mutation,
    ));
    assert_eq!(
        first.semantic_mismatch,
        format!(
            "checkpoints[{}].production.returned_log_e_bits",
            coordinates.checkpoint,
        ),
        "sufficient statistics must remain green until the corrupted wealth field",
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
        },
    );
    assert_eq!(
        first.reference_error,
        AdmissionError::ReferenceIdentityMismatch {
            expected: original.identity.root(),
            found: first.mutant.identity.root(),
        },
    );
    first
        .mutant
        .validate_payload(&config)
        .expect("resealed mutant must be internally identity-consistent");

    let first_red = red_event(&original, &first);
    let second_red = red_event(&replay, &second);
    assert_eq!(
        first_red.to_jsonl(),
        second_red.to_jsonl(),
        "seeded red evidence must be byte-stable",
    );
    assert_eq!(first_red.content_identity(), second_red.content_identity());
    let merge_panic = catch_unwind(|| assert_mergeable(&first_red))
        .expect_err("test-local merge gate must refuse the seeded red receipt");
    let merge_message = panic_message(merge_panic.as_ref());
    assert!(merge_message.contains(RED_CASE));
    assert!(merge_message.contains("PayloadIdentityMismatch"));
    assert!(merge_message.contains("ReferenceIdentityMismatch"));
    assert!(merge_message.contains(&first.semantic_mismatch));

    let retained = RetainedGolden {
        config_root: config.root(),
        result_root: original.identity.root(),
        mutant_root: first.mutant.identity.root(),
        red_event_content_hash: first_red.content_hash(),
        first_crossing_ordinals: original.record.first_crossing_ordinals,
    };
    assert_eq!(
        retained, EXPECTED_RETAINED_GOLDEN,
        "v2 retained semantic golden must be explicitly re-pinned from the centralized runtime; actual={retained:#?}",
    );

    emit_receipt(&config, &original, &first, &first_red);
}
