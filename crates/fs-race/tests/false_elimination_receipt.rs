//! Retained false-elimination calibration for the production e-race.
//!
//! This battery runs two fixed-seed Bernoulli campaigns through
//! `race_field`: an exchangeable global null, where any elimination is
//! false, and a unique-best field, where eliminating candidate zero is
//! false.  Random access Philox keys make replay independent of survivor
//! schedules.  The campaign is an empirical certifier at one recorded
//! seed and finite budget; it is not a replacement for the e-process
//! proof in the crate contract.

use fs_exec::KillRegistry;
use fs_obs::ident::{IdentityBuilder, ReplayIdentity};
use fs_obs::{Emitter, EventKind, Severity};
use fs_race::{LossSpan, RaceOutcome, RaceSettings, race_field};
use fs_rand::{Stream, StreamKey};

const SUITE: &str = "fs-race-false-elimination";
const CASE: &str = "bernoulli-calibration-v1";
const INPUT_SEED: u64 = 0xF51A_EACE_0000_0001;
const GLOBAL_NULL_KERNEL: u32 = 0xE701;
const UNIQUE_BEST_KERNEL: u32 = 0xE702;

const CANDIDATES: usize = 8;
const REPLAYS: usize = 512;
const MAX_ROUNDS: u32 = 256;
const MIN_ROUNDS: u32 = 8;
const ALPHA_NUMERATOR: u64 = 1;
const ALPHA_DENOMINATOR: u64 = 16;
const ALPHA: f64 = ALPHA_NUMERATOR as f64 / ALPHA_DENOMINATOR as f64;
const BUCKETS: u32 = 16;
const NULL_LOSS_NUMERATOR: u32 = 8;
const BEST_LOSS_NUMERATOR: u32 = 7;
const UNIQUE_BEST_CANDIDATE: usize = 0;
const BUCKET_MASK: u32 = BUCKETS - 1;
const MAX_FALSE_EVENTS: usize = 54;
#[allow(clippy::cast_lossless)] // `u64::from(MAX_ROUNDS)` is not const.
const FIXED_EVALUATIONS_PER_RACE: u64 = CANDIDATES as u64 * MAX_ROUNDS as u64;
const IID_BOUNDARY_BINOMIAL_TAIL_FORMULA: &str =
    "sum[k=maximum-false-events+1..replays] binomial-pmf(k;replays,alpha)";
const IID_BOUNDARY_BINOMIAL_TAIL_UPPER_BOUND: &str =
    "0.000078342973178362265874595602263698505767237360887";
const TWO_GATE_UNION_UPPER_BOUND: &str = "0.000156685946356724531749191204527397011534474721774";

const _: () = assert!(BUCKETS.is_power_of_two());
const _: () = assert!(BEST_LOSS_NUMERATOR <= NULL_LOSS_NUMERATOR);
const _: () = assert!(NULL_LOSS_NUMERATOR <= BUCKETS);
const _: () = assert!(UNIQUE_BEST_CANDIDATE < CANDIDATES);
const _: () = assert!(
    REPLAYS == 512 && ALPHA_NUMERATOR == 1 && ALPHA_DENOMINATOR == 16 && MAX_FALSE_EVENTS == 54
);

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Law {
    GlobalNull,
    UniqueBest,
}

impl Law {
    const fn name(self) -> &'static str {
        match self {
            Self::GlobalNull => "global-null-bernoulli-half",
            Self::UniqueBest => "unique-best-seven-sixteenths-vs-half",
        }
    }

    const fn kernel(self) -> u32 {
        match self {
            Self::GlobalNull => GLOBAL_NULL_KERNEL,
            Self::UniqueBest => UNIQUE_BEST_KERNEL,
        }
    }

    const fn loss_numerator(self, candidate: usize) -> u32 {
        match self {
            Self::GlobalNull => NULL_LOSS_NUMERATOR,
            Self::UniqueBest => {
                if candidate == UNIQUE_BEST_CANDIDATE {
                    BEST_LOSS_NUMERATOR
                } else {
                    NULL_LOSS_NUMERATOR
                }
            }
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct RunRecord {
    survivors: Vec<usize>,
    eliminated: Vec<(u32, usize)>,
    winner: usize,
    evaluations_used: u64,
    closure_calls: u64,
    fixed_n_equivalent: u64,
    rounds: u32,
    loss_span_bits: u64,
}

impl RunRecord {
    fn from_outcome(outcome: RaceOutcome, closure_calls: u64) -> Self {
        Self {
            survivors: outcome.survivors,
            eliminated: outcome.eliminated,
            winner: outcome.winner,
            evaluations_used: outcome.evaluations_used,
            closure_calls,
            fixed_n_equivalent: outcome.fixed_n_equivalent,
            rounds: outcome.rounds,
            loss_span_bits: outcome.loss_span.get().to_bits(),
        }
    }

    fn scheduled_evaluations(&self) -> Option<u64> {
        let eliminated = self
            .eliminated
            .iter()
            .try_fold(0u64, |total, &(round, _)| {
                total.checked_add(u64::from(round))
            })?;
        let survivors = usize_u64(self.survivors.len()).checked_mul(u64::from(self.rounds))?;
        eliminated.checked_add(survivors)
    }

    fn accounting_pass(&self) -> bool {
        let mut membership = [0u8; CANDIDATES];
        let survivors_valid = self
            .survivors
            .iter()
            .all(|&candidate| mark_candidate(&mut membership, candidate));
        let eliminations_valid = self.eliminated.iter().all(|&(round, candidate)| {
            (MIN_ROUNDS..=self.rounds).contains(&round)
                && mark_candidate(&mut membership, candidate)
        });

        self.closure_calls == self.evaluations_used
            && self.scheduled_evaluations() == Some(self.evaluations_used)
            && self.evaluations_used <= self.fixed_n_equivalent
            && self.fixed_n_equivalent == FIXED_EVALUATIONS_PER_RACE
            && survivors_valid
            && eliminations_valid
            && membership.iter().all(|&count| count == 1)
            && self.survivors.contains(&self.winner)
            && (MIN_ROUNDS..=MAX_ROUNDS).contains(&self.rounds)
            && self.loss_span_bits == 1.0f64.to_bits()
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct Campaign {
    law: Law,
    runs: Vec<RunRecord>,
}

impl Campaign {
    fn false_events(&self) -> usize {
        match self.law {
            Law::GlobalNull => self
                .runs
                .iter()
                .filter(|run| !run.eliminated.is_empty())
                .count(),
            Law::UniqueBest => self
                .runs
                .iter()
                .filter(|run| {
                    run.eliminated
                        .iter()
                        .any(|&(_, candidate)| candidate == UNIQUE_BEST_CANDIDATE)
                })
                .count(),
        }
    }

    fn total_evaluations(&self) -> u64 {
        self.runs.iter().map(|run| run.evaluations_used).sum()
    }

    fn total_fixed_evaluations(&self) -> u64 {
        self.runs.iter().map(|run| run.fixed_n_equivalent).sum()
    }

    fn winners_at(&self, candidate: usize) -> usize {
        self.runs
            .iter()
            .filter(|run| run.winner == candidate)
            .count()
    }

    fn calibration_pass(&self) -> bool {
        self.false_events() <= MAX_FALSE_EVENTS
    }

    fn false_event_rate(&self) -> f64 {
        self.false_events() as f64 / self.runs.len() as f64
    }

    fn evaluation_fraction(&self) -> f64 {
        self.total_evaluations() as f64 / self.total_fixed_evaluations() as f64
    }
}

#[allow(clippy::struct_excessive_bools)]
#[derive(Debug, Clone, Copy)]
struct Verdicts {
    global_null: bool,
    unique_best: bool,
    accounting: bool,
    replay: bool,
}

impl Verdicts {
    fn pass(self) -> bool {
        self.global_null && self.unique_best && self.accounting && self.replay
    }
}

fn usize_u64(value: usize) -> u64 {
    u64::try_from(value).expect("false-elimination fixture cardinality fits u64")
}

fn mark_candidate(membership: &mut [u8; CANDIDATES], candidate: usize) -> bool {
    let Some(slot) = membership.get_mut(candidate) else {
        return false;
    };
    if *slot != 0 {
        return false;
    }
    *slot = 1;
    true
}

fn registered_kills() -> KillRegistry {
    let kills = KillRegistry::new();
    for candidate in 0..CANDIDATES {
        let _ = kills.register(usize_u64(candidate));
    }
    kills
}

fn tile_for(replay: usize, candidate: usize) -> u32 {
    let tile = replay
        .checked_mul(CANDIDATES)
        .and_then(|base| base.checked_add(candidate))
        .expect("fixture tile arithmetic");
    u32::try_from(tile).expect("fixture tile fits the fs-rand identity slot")
}

fn bernoulli_loss(law: Law, replay: usize, candidate: usize, round: u64) -> f64 {
    let key = StreamKey {
        seed: INPUT_SEED,
        kernel: law.kernel(),
        tile: tile_for(replay, candidate),
    };
    let bucket = Stream::at(key, round)[0] & BUCKET_MASK;
    f64::from(u8::from(bucket < law.loss_numerator(candidate)))
}

fn settings() -> RaceSettings {
    RaceSettings {
        alpha: ALPHA,
        max_rounds: MAX_ROUNDS,
        min_rounds: MIN_ROUNDS,
        loss_span: LossSpan::ONE,
    }
}

fn run_campaign(law: Law) -> Campaign {
    let runs = (0..REPLAYS)
        .map(|replay| {
            let kills = registered_kills();
            let mut closure_calls = 0u64;
            let outcome = {
                let mut loss = |candidate: usize, round: u64| {
                    closure_calls = closure_calls
                        .checked_add(1)
                        .expect("fixture closure-call counter fits u64");
                    bernoulli_loss(law, replay, candidate, round)
                };
                race_field(&mut loss, CANDIDATES, settings(), &kills)
                    .expect("Bernoulli losses stay in the declared unit span")
            };
            RunRecord::from_outcome(outcome, closure_calls)
        })
        .collect();
    Campaign { law, runs }
}

fn config_identity() -> ReplayIdentity {
    IdentityBuilder::new("fs-race-false-elimination-config-v1")
        .str("units", "dimensionless-binary-loss")
        .str("global-null-law", Law::GlobalNull.name())
        .str("unique-best-law", Law::UniqueBest.name())
        .u64("input-seed", INPUT_SEED)
        .u64("global-null-kernel", u64::from(GLOBAL_NULL_KERNEL))
        .u64("unique-best-kernel", u64::from(UNIQUE_BEST_KERNEL))
        .u64("unique-best-candidate", usize_u64(UNIQUE_BEST_CANDIDATE))
        .str("tile-rule", "replay-index*candidate-count+candidate-index")
        .str("draw-index", "race-round-zero-based")
        .str("draw-word", "philox-word-zero-and-bucket-mask")
        .u64(
            "stream-semantics-version",
            u64::from(fs_rand::STREAM_SEMANTICS_VERSION),
        )
        .u64("candidate-count", usize_u64(CANDIDATES))
        .u64("replays-per-law", usize_u64(REPLAYS))
        .u64("maximum-rounds", u64::from(MAX_ROUNDS))
        .u64("minimum-rounds", u64::from(MIN_ROUNDS))
        .u64("alpha-numerator", ALPHA_NUMERATOR)
        .u64("alpha-denominator", ALPHA_DENOMINATOR)
        .f64_bits("alpha", ALPHA)
        .f64_bits("loss-span", LossSpan::ONE.get())
        .u64("bernoulli-buckets", u64::from(BUCKETS))
        .u64("bernoulli-bucket-mask", u64::from(BUCKET_MASK))
        .u64("global-null-loss-numerator", u64::from(NULL_LOSS_NUMERATOR))
        .u64("unique-best-loss-numerator", u64::from(BEST_LOSS_NUMERATOR))
        .u64("maximum-false-events", usize_u64(MAX_FALSE_EVENTS))
        .str(
            "calibration-envelope",
            "ceil(replays*alpha+4*sqrt(replays*alpha*(1-alpha)))",
        )
        .str(
            "iid-boundary-binomial-upper-tail-formula",
            IID_BOUNDARY_BINOMIAL_TAIL_FORMULA,
        )
        .str(
            "iid-boundary-binomial-upper-tail-upper-bound",
            IID_BOUNDARY_BINOMIAL_TAIL_UPPER_BOUND,
        )
        .str("two-gate-union-upper-bound", TWO_GATE_UNION_UPPER_BOUND)
        .u64("fixed-evaluations-per-race", FIXED_EVALUATIONS_PER_RACE)
        .str(
            "capabilities",
            "safe-rust;production-race-field;keyed-random-access-philox",
        )
        .str("execution-context", "synchronous-direct-test-no-Cx")
        .str("fs-race-version", fs_race::VERSION)
        .str("fs-eproc-version", fs_eproc::VERSION)
        .str("fs-exec-version", fs_exec::VERSION)
        .str("fs-rand-version", fs_rand::VERSION)
        .str("fs-obs-version", fs_obs::VERSION)
        .finish()
}

fn campaign_identity(config: &ReplayIdentity, campaign: &Campaign) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-race-false-elimination-campaign-v1")
        .child("config", config)
        .str("law", campaign.law.name())
        .u64("run-count", usize_u64(campaign.runs.len()));
    for (replay, run) in campaign.runs.iter().enumerate() {
        builder = builder
            .u64("replay", usize_u64(replay))
            .u64("winner", usize_u64(run.winner))
            .u64("evaluations-used", run.evaluations_used)
            .u64("closure-calls", run.closure_calls)
            .u64("fixed-n-equivalent", run.fixed_n_equivalent)
            .u64("rounds", u64::from(run.rounds))
            .u64("loss-span-bits", run.loss_span_bits)
            .u64("survivor-count", usize_u64(run.survivors.len()));
        for &survivor in &run.survivors {
            builder = builder.u64("survivor", usize_u64(survivor));
        }
        builder = builder.u64("elimination-count", usize_u64(run.eliminated.len()));
        for &(round, candidate) in &run.eliminated {
            builder = builder
                .u64("elimination-round", u64::from(round))
                .u64("eliminated-candidate", usize_u64(candidate));
        }
    }
    builder.finish()
}

fn first_campaign_mismatch(left: &Campaign, right: &Campaign) -> Option<String> {
    if left.law != right.law {
        return Some(format!("law:{:?}!={:?}", left.law, right.law));
    }
    if left.runs.len() != right.runs.len() {
        return Some(format!(
            "run-count:{}!={}",
            left.runs.len(),
            right.runs.len()
        ));
    }
    left.runs
        .iter()
        .zip(&right.runs)
        .enumerate()
        .find_map(|(replay, (a, b))| {
            (a != b).then(|| format!("replay[{replay}]:left={a:?};right={b:?}"))
        })
}

fn campaign_accounting_mismatch(label: &str, campaign: &Campaign) -> Option<String> {
    if campaign.runs.len() != REPLAYS {
        return Some(format!(
            "{label}:run-count={}!=expected-{REPLAYS}",
            campaign.runs.len()
        ));
    }
    campaign.runs.iter().enumerate().find_map(|(replay, run)| {
        (!run.accounting_pass()).then(|| {
            format!(
                "{label}[{replay}]:scheduled={:?};record={run:?}",
                run.scheduled_evaluations()
            )
        })
    })
}

fn first_accounting_mismatch(
    global: &Campaign,
    best: &Campaign,
    global_replayed: &Campaign,
    best_replayed: &Campaign,
) -> Option<String> {
    [
        ("global-null", global),
        ("unique-best", best),
        ("global-null-replay", global_replayed),
        ("unique-best-replay", best_replayed),
    ]
    .into_iter()
    .find_map(|(label, campaign)| campaign_accounting_mismatch(label, campaign))
}

#[allow(clippy::too_many_arguments)]
fn result_identity(
    config: &ReplayIdentity,
    global: &ReplayIdentity,
    global_replay: &ReplayIdentity,
    best: &ReplayIdentity,
    best_replay: &ReplayIdentity,
    global_campaign: &Campaign,
    best_campaign: &Campaign,
    global_mismatch: Option<&str>,
    best_mismatch: Option<&str>,
    accounting_mismatch: Option<&str>,
    verdicts: Verdicts,
) -> ReplayIdentity {
    IdentityBuilder::new("fs-race-false-elimination-result-v1")
        .child("config", config)
        .child("global-null", global)
        .child("global-null-replay", global_replay)
        .child("unique-best", best)
        .child("unique-best-replay", best_replay)
        .u64(
            "global-null-false-events",
            usize_u64(global_campaign.false_events()),
        )
        .u64(
            "unique-best-false-events",
            usize_u64(best_campaign.false_events()),
        )
        .u64(
            "global-null-total-evaluations",
            global_campaign.total_evaluations(),
        )
        .u64(
            "unique-best-total-evaluations",
            best_campaign.total_evaluations(),
        )
        .u64(
            "global-null-designated-candidate-winner-count",
            usize_u64(global_campaign.winners_at(UNIQUE_BEST_CANDIDATE)),
        )
        .u64(
            "unique-best-winner-count",
            usize_u64(best_campaign.winners_at(UNIQUE_BEST_CANDIDATE)),
        )
        .flag("global-null-pass", verdicts.global_null)
        .flag("unique-best-pass", verdicts.unique_best)
        .flag("accounting-pass", verdicts.accounting)
        .flag("replay-pass", verdicts.replay)
        .str(
            "global-null-first-replay-mismatch",
            global_mismatch.unwrap_or("none"),
        )
        .str(
            "unique-best-first-replay-mismatch",
            best_mismatch.unwrap_or("none"),
        )
        .str(
            "first-accounting-mismatch",
            accounting_mismatch.unwrap_or("none"),
        )
        .finish()
}

fn optional_json_string(value: Option<&str>) -> String {
    value.map_or_else(
        || "null".to_string(),
        |value| {
            let mut escaped = String::with_capacity(value.len() + 2);
            escaped.push('"');
            for character in value.chars() {
                match character {
                    '"' => escaped.push_str("\\\""),
                    '\\' => escaped.push_str("\\\\"),
                    '\n' => escaped.push_str("\\n"),
                    '\r' => escaped.push_str("\\r"),
                    '\t' => escaped.push_str("\\t"),
                    other => escaped.push(other),
                }
            }
            escaped.push('"');
            escaped
        },
    )
}

fn emit_case(emitter: &mut Emitter, case: &str, pass: bool, detail: String) {
    let event = emitter.emit(
        if pass {
            Severity::Info
        } else {
            Severity::Error
        },
        EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: case.to_string(),
            pass,
            detail,
            seed: INPUT_SEED,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("false-elimination verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("false-elimination verdict must use the fs-obs schema");
    println!("{line}");
}

fn emit_benchmark(emitter: &mut Emitter, metric: &str, value: f64) {
    let event = emitter.emit(
        Severity::Info,
        EventKind::BenchmarkResult {
            kernel: CASE.to_string(),
            metric: metric.to_string(),
            value,
            machine: 0,
        },
        None,
    );
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("false-elimination row must use the fs-obs schema");
    println!("{line}");
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn emit_receipt(
    emitter: &mut Emitter,
    config: &ReplayIdentity,
    result: &ReplayIdentity,
    global: &ReplayIdentity,
    global_replay: &ReplayIdentity,
    best: &ReplayIdentity,
    best_replay: &ReplayIdentity,
    global_campaign: &Campaign,
    best_campaign: &Campaign,
    global_mismatch: Option<&str>,
    best_mismatch: Option<&str>,
    accounting_mismatch: Option<&str>,
    verdicts: Verdicts,
) {
    let global_mismatch_json = optional_json_string(global_mismatch);
    let best_mismatch_json = optional_json_string(best_mismatch);
    let accounting_mismatch_json = optional_json_string(accounting_mismatch);
    let event = emitter.emit(
        if verdicts.pass() {
            Severity::Info
        } else {
            Severity::Error
        },
        EventKind::Custom {
            name: "e-race-false-elimination-calibration".to_string(),
            json: format!(
                "{{\"config_identity\":\"{}\",\"result_identity\":\"{}\",\
                 \"global_identity\":\"{}\",\"global_replay_identity\":\"{}\",\
                 \"unique_best_identity\":\"{}\",\
                 \"unique_best_replay_identity\":\"{}\",\
                 \"units\":\"dimensionless-binary-loss\",\
                 \"input_seed\":{INPUT_SEED},\"stream_semantics_version\":{},\
                 \"kernels\":{{\"global_null\":{GLOBAL_NULL_KERNEL},\
                 \"unique_best\":{UNIQUE_BEST_KERNEL}}},\
                 \"unique_best_candidate\":{UNIQUE_BEST_CANDIDATE},\
                 \"tile_rule\":\
                 \"replay-index*candidate-count+candidate-index\",\
                 \"draw_index\":\"race-round-zero-based\",\
                 \"draw_word\":\"philox-word-zero-and-bucket-mask\",\
                 \"candidate_count\":{CANDIDATES},\"replays_per_law\":{REPLAYS},\
                 \"alpha_numerator\":{ALPHA_NUMERATOR},\
                 \"alpha_denominator\":{ALPHA_DENOMINATOR},\
                 \"alpha_bits\":\"0x{:016x}\",\"max_rounds\":{MAX_ROUNDS},\
                 \"min_rounds\":{MIN_ROUNDS},\"loss_span\":1.0,\
                 \"bernoulli_buckets\":{BUCKETS},\
                 \"bernoulli_bucket_mask\":{BUCKET_MASK},\
                 \"null_loss_numerator\":{NULL_LOSS_NUMERATOR},\
                 \"best_loss_numerator\":{BEST_LOSS_NUMERATOR},\
                 \"maximum_false_events\":{MAX_FALSE_EVENTS},\
                 \"iid_boundary_binomial_tail_formula\":\
                 \"{IID_BOUNDARY_BINOMIAL_TAIL_FORMULA}\",\
                 \"iid_boundary_binomial_tail_upper_bound\":\
                 \"{IID_BOUNDARY_BINOMIAL_TAIL_UPPER_BOUND}\",\
                 \"two_gate_union_upper_bound\":\
                 \"{TWO_GATE_UNION_UPPER_BOUND}\",\
                 \"global_null\":{{\"false_events\":{},\"rate\":{},\
                 \"total_evaluations\":{},\"fixed_evaluations\":{},\
                 \"designated_candidate_winner_count\":{},\
                 \"first_replay_mismatch\":{global_mismatch_json},\
                 \"pass\":{}}},\"unique_best\":{{\"false_events\":{},\
                 \"rate\":{},\"total_evaluations\":{},\"fixed_evaluations\":{},\
                 \"best_candidate_winner_count\":{},\
                 \"first_replay_mismatch\":{best_mismatch_json},\
                 \"pass\":{}}},\"accounting_pass\":{},\
                 \"first_accounting_mismatch\":{accounting_mismatch_json},\
                 \"replay_pass\":{},\
                 \"versions\":{{\"fs_race\":\"{}\",\"fs_eproc\":\"{}\",\
                 \"fs_exec\":\"{}\",\"fs_rand\":\"{}\",\"fs_obs\":\"{}\"}},\
                 \"no_claims\":[\"RNG-independence-or-randomness\",\
                 \"seeds-beyond-the-recorded-corpus\",\"all-laws\",\
                 \"all-dependence-structures\",\"all-candidate-counts\",\
                 \"all-alphas\",\"all-horizons\",\
                 \"general-FWER-outside-global-null\",\
                 \"unique-best-winner-error\",\"power-or-regret\",\
                 \"performance-or-savings\",\"cross-ISA-execution\",\
                 \"Cx-or-cancellation-latency\"],\"pass\":{}}}",
                config.hex(),
                result.hex(),
                global.hex(),
                global_replay.hex(),
                best.hex(),
                best_replay.hex(),
                fs_rand::STREAM_SEMANTICS_VERSION,
                ALPHA.to_bits(),
                global_campaign.false_events(),
                global_campaign.false_event_rate(),
                global_campaign.total_evaluations(),
                global_campaign.total_fixed_evaluations(),
                global_campaign.winners_at(UNIQUE_BEST_CANDIDATE),
                verdicts.global_null,
                best_campaign.false_events(),
                best_campaign.false_event_rate(),
                best_campaign.total_evaluations(),
                best_campaign.total_fixed_evaluations(),
                best_campaign.winners_at(UNIQUE_BEST_CANDIDATE),
                verdicts.unique_best,
                verdicts.accounting,
                verdicts.replay,
                fs_race::VERSION,
                fs_eproc::VERSION,
                fs_exec::VERSION,
                fs_rand::VERSION,
                fs_obs::VERSION,
                verdicts.pass(),
            ),
        },
        None,
    );
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("false-elimination receipt must use the fs-obs schema");
    println!("{line}");
}

#[test]
#[allow(clippy::too_many_lines)]
fn production_race_retains_false_elimination_calibration_and_replay() {
    let config = config_identity();
    let global = run_campaign(Law::GlobalNull);
    let best = run_campaign(Law::UniqueBest);
    let global_replayed = run_campaign(Law::GlobalNull);
    let best_replayed = run_campaign(Law::UniqueBest);

    let global_identity = campaign_identity(&config, &global);
    let global_replay_identity = campaign_identity(&config, &global_replayed);
    let best_identity = campaign_identity(&config, &best);
    let best_replay_identity = campaign_identity(&config, &best_replayed);
    let global_mismatch = first_campaign_mismatch(&global, &global_replayed);
    let best_mismatch = first_campaign_mismatch(&best, &best_replayed);
    let accounting_mismatch =
        first_accounting_mismatch(&global, &best, &global_replayed, &best_replayed);

    let verdicts = Verdicts {
        global_null: global.calibration_pass(),
        unique_best: best.calibration_pass(),
        accounting: accounting_mismatch.is_none(),
        replay: global_mismatch.is_none()
            && best_mismatch.is_none()
            && global_identity.root() == global_replay_identity.root()
            && best_identity.root() == best_replay_identity.root(),
    };
    let result = result_identity(
        &config,
        &global_identity,
        &global_replay_identity,
        &best_identity,
        &best_replay_identity,
        &global,
        &best,
        global_mismatch.as_deref(),
        best_mismatch.as_deref(),
        accounting_mismatch.as_deref(),
        verdicts,
    );

    let mut emitter = Emitter::new(SUITE, CASE);
    emit_receipt(
        &mut emitter,
        &config,
        &result,
        &global_identity,
        &global_replay_identity,
        &best_identity,
        &best_replay_identity,
        &global,
        &best,
        global_mismatch.as_deref(),
        best_mismatch.as_deref(),
        accounting_mismatch.as_deref(),
        verdicts,
    );
    emit_benchmark(
        &mut emitter,
        "global_null_any_elimination_rate",
        global.false_event_rate(),
    );
    emit_benchmark(
        &mut emitter,
        "unique_best_false_elimination_rate",
        best.false_event_rate(),
    );
    emit_benchmark(
        &mut emitter,
        "global_null_evaluation_fraction",
        global.evaluation_fraction(),
    );
    emit_benchmark(
        &mut emitter,
        "unique_best_evaluation_fraction",
        best.evaluation_fraction(),
    );
    emit_case(
        &mut emitter,
        "global-null-any-elimination",
        verdicts.global_null,
        format!(
            "config={}; result={}; false events={}/{REPLAYS} (rate={:.6}) \
             <= {MAX_FALSE_EVENTS}; \
             campaign={}; total evaluations={}/{}",
            config.hex(),
            result.hex(),
            global.false_events(),
            global.false_event_rate(),
            global_identity.hex(),
            global.total_evaluations(),
            global.total_fixed_evaluations(),
        ),
    );
    emit_case(
        &mut emitter,
        "unique-best-false-elimination",
        verdicts.unique_best,
        format!(
            "config={}; result={}; best-candidate false eliminations={}/{REPLAYS} \
             (rate={:.6}) <= {MAX_FALSE_EVENTS}; campaign={}; \
             best-candidate winners={}/{REPLAYS}",
            config.hex(),
            result.hex(),
            best.false_events(),
            best.false_event_rate(),
            best_identity.hex(),
            best.winners_at(UNIQUE_BEST_CANDIDATE),
        ),
    );
    emit_case(
        &mut emitter,
        "closure-and-budget-accounting",
        verdicts.accounting,
        format!(
            "config={}; result={}; all four campaigns retain closure_calls == \
             evaluations_used == independently scheduled evaluations <= fixed_n; \
             fixed evaluations/race={FIXED_EVALUATIONS_PER_RACE}; \
             first mismatch={accounting_mismatch:?}",
            config.hex(),
            result.hex(),
        ),
    );
    emit_case(
        &mut emitter,
        "full-campaign-bitwise-replay",
        verdicts.replay,
        format!(
            "config={}; result={}; global roots={}/{} mismatch={global_mismatch:?}; \
             unique-best roots={}/{} mismatch={best_mismatch:?}",
            config.hex(),
            result.hex(),
            global_identity.hex(),
            global_replay_identity.hex(),
            best_identity.hex(),
            best_replay_identity.hex(),
        ),
    );

    assert!(
        verdicts.pass(),
        "false-elimination campaign failed: verdicts={verdicts:?}; \
         global false={}; best false={}; global mismatch={global_mismatch:?}; \
         best mismatch={best_mismatch:?}; accounting mismatch={accounting_mismatch:?}",
        global.false_events(),
        best.false_events(),
    );
}
