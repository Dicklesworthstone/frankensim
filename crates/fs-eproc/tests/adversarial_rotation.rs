//! Certifier-of-certifier trial for [`fs_eproc::BettingEProcess`].
//!
//! An adversary chooses the next bounded null law from the observed history and
//! stops at data-dependent times. Every available law still has conditional
//! mean exactly one half, so the one-sided null remains valid. The finite,
//! seeded campaign checks that threshold crossings stay inside a conservative
//! binomial envelope and retains enough schedule/result provenance to replay
//! the empirical claim.

use fs_eproc::BettingEProcess;
use fs_obs::ident::{IdentityBuilder, ReplayIdentity};
use fs_obs::{Emitter, EventKind, Severity};
use fs_rand::{Stream, StreamKey};

const SUITE: &str = "fs-eproc/adversarial-rotation";
const CASE: &str = "null-law-rotation";
const INPUT_SEED: u64 = 0xAD0E_2026_0717_0021;
const RNG_KERNEL: u32 = 1;
const SIMS_PER_STRATEGY: u32 = 1_024;
const HORIZON: u32 = 1_024;
const ALPHA: f64 = 0.05;
// Under a level-0.05 null, E[crossings] <= 51.2 and the binomial standard
// deviation is <= 6.98. Eighty is more than four standard deviations above
// that boundary, keeping this fixed empirical gate conservative.
const MAX_CROSSINGS_PER_STRATEGY: u32 = 80;
const WEALTH_NEAR_PEAK: f64 = 0.25;
const WEALTH_LOW: f64 = -2.0;
const OUTCOME_HIGH: f64 = 0.75;
const OUTCOME_LOW: f64 = 0.25;
const PEAK_NEAR: f64 = 0.5;
const WEALTH_STOP_AFTER: u32 = 128;
const WEALTH_STOP_FLOOR: f64 = -4.0;
const OUTCOME_STOP_AFTER: u32 = 64;
const OUTCOME_HIGH_RUN: u32 = 8;
const RETREAT_STOP_AFTER: u32 = 128;
const RETREAT_LOG_UNITS: f64 = 2.0;
const INITIAL_PREVIOUS: f64 = 0.5;
const INITIAL_LOG_E: f64 = 0.0;
const INITIAL_PEAK_LOG_E: f64 = 0.0;
const INITIAL_HIGH_RUN: u32 = 0;
const HIGH_OUTCOME_THRESHOLD: f64 = 0.5;

#[derive(Clone, Copy)]
struct NullLaw {
    name: &'static str,
    low_numerator: u64,
    high_numerator: u64,
    value_denominator: u64,
    buckets: u64,
    high_buckets: u64,
}

const NULL_LAWS: [NullLaw; 4] = [
    NullLaw {
        name: "extremal",
        low_numerator: 0,
        high_numerator: 1,
        value_denominator: 1,
        buckets: 2,
        high_buckets: 1,
    },
    NullLaw {
        name: "narrow",
        low_numerator: 1,
        high_numerator: 3,
        value_denominator: 4,
        buckets: 2,
        high_buckets: 1,
    },
    NullLaw {
        name: "wide",
        low_numerator: 1,
        high_numerator: 15,
        value_denominator: 16,
        buckets: 2,
        high_buckets: 1,
    },
    NullLaw {
        name: "rare-high",
        low_numerator: 3,
        high_numerator: 7,
        value_denominator: 8,
        buckets: 4,
        high_buckets: 1,
    },
];

impl NullLaw {
    fn is_exact_bounded_null(self) -> bool {
        if self.value_denominator == 0
            || !self.value_denominator.is_power_of_two()
            || self.value_denominator > (1u64 << 53)
            || self.buckets == 0
            || !self.buckets.is_power_of_two()
            || self.high_buckets > self.buckets
            || self.low_numerator > self.value_denominator
            || self.high_numerator > self.value_denominator
        {
            return false;
        }
        let low_buckets = self.buckets - self.high_buckets;
        let weighted_numerator = u128::from(low_buckets) * u128::from(self.low_numerator)
            + u128::from(self.high_buckets) * u128::from(self.high_numerator);
        2 * weighted_numerator == u128::from(self.buckets) * u128::from(self.value_denominator)
    }
}

#[derive(Clone, Copy)]
enum Strategy {
    Cyclic,
    WealthChasing,
    OutcomeReactive,
    PeakRetreat,
}

const STRATEGIES: [Strategy; 4] = [
    Strategy::Cyclic,
    Strategy::WealthChasing,
    Strategy::OutcomeReactive,
    Strategy::PeakRetreat,
];

impl Strategy {
    const fn name(self) -> &'static str {
        match self {
            Self::Cyclic => "cyclic",
            Self::WealthChasing => "wealth-chasing",
            Self::OutcomeReactive => "outcome-reactive",
            Self::PeakRetreat => "peak-retreat",
        }
    }

    const fn stop_rule(self) -> &'static str {
        match self {
            Self::Cyclic => "threshold-or-horizon",
            Self::WealthChasing => "threshold-or-log-wealth-below-minus-four-after-128",
            Self::OutcomeReactive => "threshold-or-eight-high-outcomes-after-64",
            Self::PeakRetreat => "threshold-or-two-log-units-below-peak-after-128",
        }
    }

    const fn rotation_rule(self) -> &'static str {
        match self {
            Self::Cyclic => "law=(simulation+step)&3",
            Self::WealthChasing => {
                "if log_e>=peak-wealth-near-peak:0; else if log_e<wealth-low:1; else if previous>.5:2; else:3"
            }
            Self::OutcomeReactive => {
                "if previous>=outcome-high:0; else if previous<=outcome-low:3; else if even-step:1; else:2"
            }
            Self::PeakRetreat => "if peak-log_e<=peak-near:2; else if even-step:1; else:3",
        }
    }

    fn law(self, sim: u32, step: u32, log_e: f64, peak_log_e: f64, previous: f64) -> usize {
        match self {
            Self::Cyclic => {
                usize::try_from(sim.wrapping_add(step) & 3).expect("two-bit law index fits usize")
            }
            Self::WealthChasing => {
                if log_e >= peak_log_e - WEALTH_NEAR_PEAK {
                    0
                } else if log_e < WEALTH_LOW {
                    1
                } else if previous > 0.5 {
                    2
                } else {
                    3
                }
            }
            Self::OutcomeReactive => {
                if previous >= OUTCOME_HIGH {
                    0
                } else if previous <= OUTCOME_LOW {
                    3
                } else if step & 1 == 0 {
                    1
                } else {
                    2
                }
            }
            Self::PeakRetreat => {
                if peak_log_e - log_e <= PEAK_NEAR {
                    2
                } else if step & 1 == 0 {
                    1
                } else {
                    3
                }
            }
        }
    }

    fn policy_stop(self, observations: u32, log_e: f64, peak_log_e: f64, high_run: u32) -> bool {
        match self {
            Self::Cyclic => false,
            Self::WealthChasing => observations >= WEALTH_STOP_AFTER && log_e <= WEALTH_STOP_FLOOR,
            Self::OutcomeReactive => {
                observations >= OUTCOME_STOP_AFTER && high_run >= OUTCOME_HIGH_RUN
            }
            Self::PeakRetreat => {
                observations >= RETREAT_STOP_AFTER && peak_log_e - log_e >= RETREAT_LOG_UNITS
            }
        }
    }
}

struct StrategyResult {
    strategy: Strategy,
    strategy_index: u32,
    crossings: u32,
    policy_stops: u32,
    horizon_stops: u32,
    total_observations: u64,
    max_stop: u64,
    law_counts: [u64; 4],
}

impl StrategyResult {
    fn pass(&self) -> bool {
        self.crossings <= MAX_CROSSINGS_PER_STRATEGY
            && self.crossings + self.policy_stops + self.horizon_stops == SIMS_PER_STRATEGY
            && self.law_counts.iter().sum::<u64>() == self.total_observations
            && self.total_observations >= u64::from(SIMS_PER_STRATEGY)
            && self.total_observations <= u64::from(SIMS_PER_STRATEGY) * u64::from(HORIZON)
            && self.max_stop <= u64::from(HORIZON)
    }

    fn crossing_rate(&self) -> f64 {
        f64::from(self.crossings) / f64::from(SIMS_PER_STRATEGY)
    }

    fn mean_stop(&self) -> f64 {
        self.total_observations as f64 / f64::from(SIMS_PER_STRATEGY)
    }

    fn first_tile(&self) -> u32 {
        self.strategy_index * SIMS_PER_STRATEGY
    }

    fn last_tile(&self) -> u32 {
        self.first_tile() + SIMS_PER_STRATEGY - 1
    }
}

/// Draw from one mechanically checked dyadic null law. Selection happens
/// before this draw, from prior history only.
fn draw_null(stream: &mut Stream, law: usize) -> f64 {
    let law = NULL_LAWS
        .get(law)
        .expect("the strategy selects one of four null laws");
    let high = stream.next_below(law.buckets) < law.high_buckets;
    let numerator = if high {
        law.high_numerator
    } else {
        law.low_numerator
    };
    numerator as f64 / law.value_denominator as f64
}

fn run_strategy(strategy: Strategy, strategy_index: u32) -> StrategyResult {
    let mut result = StrategyResult {
        strategy,
        strategy_index,
        crossings: 0,
        policy_stops: 0,
        horizon_stops: 0,
        total_observations: 0,
        max_stop: 0,
        law_counts: [0; 4],
    };

    for sim in 0..SIMS_PER_STRATEGY {
        let tile = strategy_index * SIMS_PER_STRATEGY + sim;
        let mut stream = StreamKey {
            seed: INPUT_SEED,
            kernel: RNG_KERNEL,
            tile,
        }
        .stream();
        let mut process = BettingEProcess::new(0.5);
        let mut previous = INITIAL_PREVIOUS;
        let mut peak_log_e = INITIAL_PEAK_LOG_E;
        let mut high_run = INITIAL_HIGH_RUN;
        let mut policy_stopped = false;
        let mut crossed = false;

        for step in 0..HORIZON {
            // The adversary sees only F_(t-1) here; the current draw has not
            // happened, so choosing among mean-one-half laws is predictable.
            let law = strategy.law(sim, step, process.log_e_value(), peak_log_e, previous);
            result.law_counts[law] += 1;
            let value = draw_null(&mut stream, law);
            let log_e = process.observe(value);
            previous = value;
            high_run = if value > HIGH_OUTCOME_THRESHOLD {
                high_run + 1
            } else {
                0
            };
            peak_log_e = peak_log_e.max(log_e);

            if process.rejects_at(ALPHA) {
                crossed = true;
                break;
            }
            if strategy.policy_stop(step + 1, log_e, peak_log_e, high_run) {
                policy_stopped = true;
                break;
            }
        }

        result.crossings += u32::from(crossed);
        result.policy_stops += u32::from(policy_stopped);
        result.horizon_stops += u32::from(!crossed && !policy_stopped);
        result.total_observations += process.len();
        result.max_stop = result.max_stop.max(process.len());
    }
    result
}

fn campaign_identity() -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-eproc-adversarial-null-rotation-config-v1")
        .str("null", "conditional-mean-at-most-one-half")
        .f64_bits("null-mean", 0.5)
        .f64_bits("alpha", ALPHA)
        .u64("simulations-per-strategy", u64::from(SIMS_PER_STRATEGY))
        .u64("horizon", u64::from(HORIZON))
        .u64(
            "maximum-crossings-per-strategy",
            u64::from(MAX_CROSSINGS_PER_STRATEGY),
        )
        .u64("input-seed", INPUT_SEED)
        .u64("rng-kernel", u64::from(RNG_KERNEL))
        .u64(
            "stream-semantics-version",
            u64::from(fs_rand::STREAM_SEMANTICS_VERSION),
        )
        .str("tile-layout", "strategy-index*tile-stride+simulation")
        .u64("tile-stride", u64::from(SIMS_PER_STRATEGY))
        .str("draw-method", "Stream::next_below(power-of-two-buckets)-v1")
        .u64("draws-per-observation", 1)
        .u64("stream-index-start", 0)
        .u64("stream-index-exclusive-limit", u64::from(HORIZON))
        .flag("law-selected-before-draw", true)
        .f64_bits("wealth-near-peak", WEALTH_NEAR_PEAK)
        .f64_bits("wealth-low", WEALTH_LOW)
        .f64_bits("outcome-high", OUTCOME_HIGH)
        .f64_bits("outcome-low", OUTCOME_LOW)
        .f64_bits("peak-near", PEAK_NEAR)
        .u64("wealth-stop-after", u64::from(WEALTH_STOP_AFTER))
        .f64_bits("wealth-stop-floor", WEALTH_STOP_FLOOR)
        .u64("outcome-stop-after", u64::from(OUTCOME_STOP_AFTER))
        .u64("outcome-high-run", u64::from(OUTCOME_HIGH_RUN))
        .u64("retreat-stop-after", u64::from(RETREAT_STOP_AFTER))
        .f64_bits("retreat-log-units", RETREAT_LOG_UNITS)
        .f64_bits("initial-previous", INITIAL_PREVIOUS)
        .f64_bits("initial-log-e", INITIAL_LOG_E)
        .f64_bits("initial-peak-log-e", INITIAL_PEAK_LOG_E)
        .u64("initial-high-run", u64::from(INITIAL_HIGH_RUN))
        .f64_bits("high-outcome-threshold", HIGH_OUTCOME_THRESHOLD)
        .str(
            "observation-order",
            "select-law-pre-draw>draw>observe>update-history>threshold-stop>policy-stop>horizon",
        )
        .str(
            "high-run-update",
            "if value>high-outcome-threshold then previous+1 else 0",
        )
        .str("fs-eproc-version", fs_eproc::VERSION)
        .str("fs-rand-version", fs_rand::VERSION);
    for law in NULL_LAWS {
        builder = builder
            .str("law", law.name)
            .u64("low-numerator", law.low_numerator)
            .u64("high-numerator", law.high_numerator)
            .u64("value-denominator", law.value_denominator)
            .u64("buckets", law.buckets)
            .u64("high-buckets", law.high_buckets);
    }
    for strategy in STRATEGIES {
        builder = builder
            .str("strategy", strategy.name())
            .str("rotation-rule", strategy.rotation_rule())
            .str("stopping-rule", strategy.stop_rule());
    }
    builder.finish()
}

fn result_identity(campaign: &ReplayIdentity, results: &[StrategyResult]) -> ReplayIdentity {
    let mut builder = IdentityBuilder::new("fs-eproc-adversarial-null-rotation-result-v1")
        .child("campaign", campaign);
    for result in results {
        builder = builder
            .str("strategy", result.strategy.name())
            .u64("crossings", u64::from(result.crossings))
            .u64("policy-stops", u64::from(result.policy_stops))
            .u64("horizon-stops", u64::from(result.horizon_stops))
            .u64("total-observations", result.total_observations)
            .u64("maximum-stop", result.max_stop);
        for &count in &result.law_counts {
            builder = builder.u64("law-count", count);
        }
    }
    builder.finish()
}

fn law_counts_json(counts: &[u64; 4]) -> String {
    format!("[{},{},{},{}]", counts[0], counts[1], counts[2], counts[3])
}

fn null_laws_json() -> String {
    NULL_LAWS
        .iter()
        .map(|law| {
            format!(
                "{{\"name\":\"{}\",\"low_numerator\":{},\"high_numerator\":{},\
                 \"value_denominator\":{},\"buckets\":{},\"high_buckets\":{},\
                 \"exact_bounded_null\":{}}}",
                law.name,
                law.low_numerator,
                law.high_numerator,
                law.value_denominator,
                law.buckets,
                law.high_buckets,
                law.is_exact_bounded_null(),
            )
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn strategy_json(result: &StrategyResult, laws_valid: bool) -> String {
    format!(
        "{{\"strategy\":\"{}\",\"rotation_rule\":\"{}\",\"stop_rule\":\"{}\",\
         \"tile_first\":{},\"tile_last\":{},\"crossings\":{},\
         \"crossing_rate\":{},\"maximum_crossings\":{},\
         \"policy_stops\":{},\"horizon_stops\":{},\
         \"total_observations\":{},\"mean_stop\":{},\"maximum_stop\":{},\
         \"law_counts\":{},\"pass\":{}}}",
        result.strategy.name(),
        result.strategy.rotation_rule(),
        result.strategy.stop_rule(),
        result.first_tile(),
        result.last_tile(),
        result.crossings,
        result.crossing_rate(),
        MAX_CROSSINGS_PER_STRATEGY,
        result.policy_stops,
        result.horizon_stops,
        result.total_observations,
        result.mean_stop(),
        result.max_stop,
        law_counts_json(&result.law_counts),
        laws_valid && result.pass(),
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
    fs_obs::lint_failure_record(&event).expect("rotation verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("rotation verdict must use the fs-obs wire schema");
    println!("{line}");
}

fn aggregate_law_counts(results: &[StrategyResult]) -> [u64; 4] {
    let mut aggregate = [0u64; 4];
    for row in results {
        for (total, count) in aggregate.iter_mut().zip(row.law_counts) {
            *total += count;
        }
    }
    aggregate
}

fn emitter_with_receipt(
    campaign: &ReplayIdentity,
    result: &ReplayIdentity,
    results: &[StrategyResult],
    aggregate_law_counts: &[u64; 4],
    laws_valid: bool,
    pass: bool,
    first_failure: Option<&str>,
) -> Emitter {
    let rows_json = results
        .iter()
        .map(|row| strategy_json(row, laws_valid))
        .collect::<Vec<_>>()
        .join(",");
    let first_failure_json =
        first_failure.map_or_else(|| "null".to_string(), |failure| format!("\"{failure}\""));
    let mut emitter = Emitter::new(SUITE, CASE);
    let receipt = emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "e-process-adversarial-null-rotation".to_string(),
            json: format!(
                "{{\"campaign_identity\":\"{}\",\"result_identity\":\"{}\",\
                 \"input_seed\":{INPUT_SEED},\"rng_kernel\":{RNG_KERNEL},\
                 \"stream_semantics_version\":{},\"simulations_per_strategy\":{},\
                 \"tile_stride\":{},\"draw_method\":\"next_below-power-of-two-v1\",\
                 \"draws_per_observation\":1,\"stream_index_start\":0,\
                 \"horizon\":{HORIZON},\"alpha_bits\":\"0x{:016x}\",\
                 \"maximum_crossings_per_strategy\":{},\"null_laws\":[{}],\
                 \"laws_valid\":{laws_valid},\"aggregate_law_counts\":{},\
                 \"strategies\":[{rows_json}],\"first_failure\":{first_failure_json},\
                 \"pass\":{pass}}}",
                campaign.hex(),
                result.hex(),
                fs_rand::STREAM_SEMANTICS_VERSION,
                SIMS_PER_STRATEGY,
                SIMS_PER_STRATEGY,
                ALPHA.to_bits(),
                MAX_CROSSINGS_PER_STRATEGY,
                null_laws_json(),
                law_counts_json(aggregate_law_counts),
            ),
        },
        None,
    );
    let line = receipt.to_jsonl();
    fs_obs::validate_line(&line).expect("rotation receipt must use the fs-obs wire schema");
    println!("{line}");
    emitter
}

#[test]
fn adversarial_null_rotation_preserves_anytime_validity() {
    let campaign = campaign_identity();
    let results: Vec<_> = STRATEGIES
        .into_iter()
        .enumerate()
        .map(|(index, strategy)| {
            run_strategy(
                strategy,
                u32::try_from(index).expect("strategy index fits u32"),
            )
        })
        .collect();
    let result = result_identity(&campaign, &results);
    let aggregate_law_counts = aggregate_law_counts(&results);
    let laws_valid = NULL_LAWS.iter().all(|law| law.is_exact_bounded_null());
    let all_laws_exercised = aggregate_law_counts.iter().all(|&count| count > 0);
    let pass = laws_valid && all_laws_exercised && results.iter().all(StrategyResult::pass);
    let first_failure = if !laws_valid {
        Some("null-law-definition")
    } else if !all_laws_exercised {
        Some("law-coverage")
    } else {
        results
            .iter()
            .find(|row| !row.pass())
            .map(|row| row.strategy.name())
    };
    let mut emitter = emitter_with_receipt(
        &campaign,
        &result,
        &results,
        &aggregate_law_counts,
        laws_valid,
        pass,
        first_failure,
    );

    for row in &results {
        let row_pass = laws_valid && row.pass();
        emit_case(
            &mut emitter,
            &format!("{CASE}/{}", row.strategy.name()),
            row_pass,
            format!(
                "campaign={}; result={}; crossings={}/{} (rate={:.6}, max={}); \
                 stops threshold/policy/horizon={}/{}/{}; observations={} (mean={:.3}, max={}); \
                 laws_valid={laws_valid}; laws={:?}",
                campaign.hex(),
                result.hex(),
                row.crossings,
                SIMS_PER_STRATEGY,
                row.crossing_rate(),
                MAX_CROSSINGS_PER_STRATEGY,
                row.crossings,
                row.policy_stops,
                row.horizon_stops,
                row.total_observations,
                row.mean_stop(),
                row.max_stop,
                row.law_counts,
            ),
        );
    }
    emit_case(
        &mut emitter,
        CASE,
        pass,
        format!(
            "campaign={}; result={}; four predictable rotations x {} trials x horizon {}; \
             alpha={ALPHA}; max crossings/strategy={MAX_CROSSINGS_PER_STRATEGY}; \
             laws_valid={laws_valid}; aggregate laws={aggregate_law_counts:?}; \
             first_failure={first_failure:?}",
            campaign.hex(),
            result.hex(),
            SIMS_PER_STRATEGY,
            HORIZON,
        ),
    );
    assert!(
        pass,
        "adversarial null-rotation campaign failed: {first_failure:?}"
    );
}
