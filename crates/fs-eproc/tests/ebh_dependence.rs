//! G0 certifier-of-certifier battery for e-BH under strong dependence.
//!
//! Every true-null e-value in the finite laws below has expectation exactly
//! one.  The laws then make thousands of hypotheses perfectly dependent in
//! blocks, mutually exclusive across blocks, or perfectly dependent across
//! the entire null family.  Exhaustive phase enumeration avoids Monte Carlo
//! slack while exercising the production step-up implementation and its
//! deterministic recovery of original hypothesis indices.

use fs_eproc::e_benjamini_hochberg;
use fs_math::det;
use fs_obs::ident::{IdentityBuilder, ReplayIdentity};
use fs_obs::{Emitter, EventKind, Severity};
use fs_rand::StreamKey;

const SUITE: &str = "fs-eproc/e-bh-dependence";
const CASE: &str = "strong-dependence";
const INPUT_SEED: u64 = 0xEB14_2026_0717_0022;
const RNG_KERNEL: u32 = 1;
const ROTATION_TILE: u32 = 0;
const FAMILY_SIZE: usize = 4_096;
const ALPHA_NUMERATOR: usize = 1;
const ALPHA_DENOMINATOR: usize = 16;
const ALPHA: f64 = ALPHA_NUMERATOR as f64 / ALPHA_DENOMINATOR as f64;
const ZERO_LOG_E: f64 = f64::NEG_INFINITY;

const GLOBAL_BLOCKS: usize = 64;
const GLOBAL_BLOCK_SIZE: usize = 64;
const GLOBAL_PHASES: usize = 1_024;
const GLOBAL_ACTIVE_E: u64 = 1_024;

const PERFECT_ALTERNATIVES: usize = 20;
const PERFECT_NULLS: usize = FAMILY_SIZE - PERFECT_ALTERNATIVES;
const PERFECT_PHASES: usize = 17;
const PERFECT_NULL_E: u64 = 17;

const MUTANT_ALTERNATIVES: usize = 64;
const MUTANT_NULLS: usize = FAMILY_SIZE - MUTANT_ALTERNATIVES;
const MUTANT_BLOCKS: usize = 63;
const MUTANT_BLOCK_SIZE: usize = 64;
const MUTANT_PHASES: usize = MUTANT_BLOCKS;
const MUTANT_NULL_E: u64 = 63;
const ALTERNATIVE_E: u64 = 4_096;
const MUTANT_NAME: &str = "missing-family-size-factor";

#[derive(Debug)]
struct GlobalResult {
    active_phases: usize,
    inactive_phases: usize,
    total_rejections: usize,
    first_shape_failure: Option<usize>,
    marginal_means_exact: bool,
}

impl GlobalResult {
    fn pass(&self) -> bool {
        self.active_phases == GLOBAL_BLOCKS
            && self.inactive_phases == GLOBAL_PHASES - GLOBAL_BLOCKS
            && self.total_rejections == GLOBAL_BLOCKS * GLOBAL_BLOCK_SIZE
            && self.first_shape_failure.is_none()
            && self.marginal_means_exact
            && self.active_phases * ALPHA_DENOMINATOR == GLOBAL_PHASES * ALPHA_NUMERATOR
    }
}

#[derive(Debug)]
struct PerfectResult {
    shock_rejections: usize,
    shock_false_rejections: usize,
    ordinary_rejections: usize,
    ordinary_false_rejections: usize,
    first_shape_failure: Option<usize>,
    marginal_means_exact: bool,
}

impl PerfectResult {
    fn pass(&self) -> bool {
        self.shock_rejections == FAMILY_SIZE
            && self.shock_false_rejections == PERFECT_NULLS
            && self.ordinary_rejections == (PERFECT_PHASES - 1) * PERFECT_ALTERNATIVES
            && self.ordinary_false_rejections == 0
            && self.first_shape_failure.is_none()
            && self.marginal_means_exact
            && PERFECT_NULLS * ALPHA_DENOMINATOR <= FAMILY_SIZE * PERFECT_PHASES * ALPHA_NUMERATOR
    }
}

#[derive(Debug)]
struct MutantResult {
    correct_rejections: usize,
    correct_false_rejections: usize,
    mutant_rejections: usize,
    mutant_false_rejections: usize,
    first_shape_failure: Option<usize>,
    marginal_means_exact: bool,
}

impl MutantResult {
    fn pass(&self) -> bool {
        self.correct_rejections == MUTANT_PHASES * MUTANT_ALTERNATIVES
            && self.correct_false_rejections == 0
            && self.mutant_rejections == MUTANT_PHASES * (MUTANT_ALTERNATIVES + MUTANT_BLOCK_SIZE)
            && self.mutant_false_rejections == MUTANT_PHASES * MUTANT_BLOCK_SIZE
            && self.first_shape_failure.is_none()
            && self.marginal_means_exact
    }
}

fn input_rotation() -> usize {
    let mut stream = StreamKey {
        seed: INPUT_SEED,
        kernel: RNG_KERNEL,
        tile: ROTATION_TILE,
    }
    .stream();
    let rotation = stream.next_below(FAMILY_SIZE as u64);
    assert_eq!(
        stream.index(),
        1,
        "power-of-two rotation draw must consume exactly one Philox position"
    );
    usize::try_from(rotation).expect("rotation is smaller than the family size")
}

fn hypothesis_at(position: usize, rotation: usize) -> usize {
    (position + rotation) & (FAMILY_SIZE - 1)
}

fn sorted_positions(range: core::ops::Range<usize>, rotation: usize) -> Vec<usize> {
    let mut indices: Vec<_> = range
        .map(|position| hypothesis_at(position, rotation))
        .collect();
    indices.sort_unstable();
    indices
}

fn mark_positions(log_e: &mut [f64], range: core::ops::Range<usize>, rotation: usize, value: f64) {
    for position in range {
        log_e[hypothesis_at(position, rotation)] = value;
    }
}

fn global_null_trial(rotation: usize) -> GlobalResult {
    let active_log_e = det::ln(GLOBAL_ACTIVE_E as f64);
    let mut log_e = vec![ZERO_LOG_E; FAMILY_SIZE];
    let mut activation_counts = vec![0u16; FAMILY_SIZE];
    let mut active_phases = 0;
    let mut inactive_phases = 0;
    let mut total_rejections = 0;
    let mut first_shape_failure = None;

    for phase in 0..GLOBAL_PHASES {
        log_e.fill(ZERO_LOG_E);
        let expected = if phase < GLOBAL_BLOCKS {
            let start = phase * GLOBAL_BLOCK_SIZE;
            let end = start + GLOBAL_BLOCK_SIZE;
            mark_positions(&mut log_e, start..end, rotation, active_log_e);
            for position in start..end {
                activation_counts[hypothesis_at(position, rotation)] += 1;
            }
            sorted_positions(start..end, rotation)
        } else {
            Vec::new()
        };
        let rejected = e_benjamini_hochberg(&log_e, ALPHA);
        if rejected.is_empty() {
            inactive_phases += 1;
        } else {
            active_phases += 1;
        }
        total_rejections += rejected.len();
        if rejected != expected && first_shape_failure.is_none() {
            first_shape_failure = Some(phase);
        }
    }

    GlobalResult {
        active_phases,
        inactive_phases,
        total_rejections,
        first_shape_failure,
        marginal_means_exact: activation_counts.iter().all(|&count| count == 1)
            && GLOBAL_ACTIVE_E == GLOBAL_PHASES as u64,
    }
}

fn perfect_dependence_trial(rotation: usize) -> PerfectResult {
    let alternative_log_e = det::ln(ALTERNATIVE_E as f64);
    let null_log_e = det::ln(PERFECT_NULL_E as f64);
    let alternatives = sorted_positions(0..PERFECT_ALTERNATIVES, rotation);
    let mut is_alternative = vec![false; FAMILY_SIZE];
    for &hypothesis in &alternatives {
        is_alternative[hypothesis] = true;
    }
    let mut activation_counts = vec![0u8; FAMILY_SIZE];
    let mut log_e = vec![ZERO_LOG_E; FAMILY_SIZE];
    let mut result = PerfectResult {
        shock_rejections: 0,
        shock_false_rejections: 0,
        ordinary_rejections: 0,
        ordinary_false_rejections: 0,
        first_shape_failure: None,
        marginal_means_exact: false,
    };

    for phase in 0..PERFECT_PHASES {
        log_e.fill(ZERO_LOG_E);
        mark_positions(
            &mut log_e,
            0..PERFECT_ALTERNATIVES,
            rotation,
            alternative_log_e,
        );
        let shock = phase == 0;
        if shock {
            mark_positions(
                &mut log_e,
                PERFECT_ALTERNATIVES..FAMILY_SIZE,
                rotation,
                null_log_e,
            );
            for position in PERFECT_ALTERNATIVES..FAMILY_SIZE {
                activation_counts[hypothesis_at(position, rotation)] += 1;
            }
        }
        let rejected = e_benjamini_hochberg(&log_e, ALPHA);
        let expected = if shock {
            (0..FAMILY_SIZE).collect()
        } else {
            alternatives.clone()
        };
        let false_rejections = rejected
            .iter()
            .filter(|&&hypothesis| !is_alternative[hypothesis])
            .count();
        if shock {
            result.shock_rejections = rejected.len();
            result.shock_false_rejections = false_rejections;
        } else {
            result.ordinary_rejections += rejected.len();
            result.ordinary_false_rejections += false_rejections;
        }
        if rejected != expected && result.first_shape_failure.is_none() {
            result.first_shape_failure = Some(phase);
        }
    }
    result.marginal_means_exact = activation_counts
        .iter()
        .enumerate()
        .all(|(hypothesis, &count)| count == u8::from(!is_alternative[hypothesis]))
        && PERFECT_NULL_E == PERFECT_PHASES as u64;
    result
}

/// Deliberately broken local comparator: omitting `m` from `m / (alpha*k)`
/// makes the threshold collapse as the family grows.  This mutant is never
/// used by production code; the battery requires the exact strong-dependence
/// fixture to distinguish it from the shipped implementation.
fn missing_family_size_e_bh(log_e: &[f64], alpha: f64) -> Vec<usize> {
    let mut order: Vec<usize> = (0..log_e.len()).collect();
    order.sort_by(|&a, &b| log_e[b].total_cmp(&log_e[a]).then(a.cmp(&b)));
    let mut k_hat = 0;
    for (rank0, &hypothesis) in order.iter().enumerate() {
        let k = rank0 + 1;
        let broken_threshold = det::ln(1.0 / (alpha * k as f64));
        if log_e[hypothesis] >= broken_threshold {
            k_hat = k;
        }
    }
    order.truncate(k_hat);
    order.sort_unstable();
    order
}

fn mutant_trial(rotation: usize) -> MutantResult {
    let alternative_log_e = det::ln(ALTERNATIVE_E as f64);
    let null_log_e = det::ln(MUTANT_NULL_E as f64);
    let alternatives = sorted_positions(0..MUTANT_ALTERNATIVES, rotation);
    let mut is_alternative = vec![false; FAMILY_SIZE];
    for &hypothesis in &alternatives {
        is_alternative[hypothesis] = true;
    }
    let mut activation_counts = vec![0u8; FAMILY_SIZE];
    let mut result = MutantResult {
        correct_rejections: 0,
        correct_false_rejections: 0,
        mutant_rejections: 0,
        mutant_false_rejections: 0,
        first_shape_failure: None,
        marginal_means_exact: false,
    };
    let mut log_e = vec![ZERO_LOG_E; FAMILY_SIZE];

    for phase in 0..MUTANT_PHASES {
        log_e.fill(ZERO_LOG_E);
        mark_positions(
            &mut log_e,
            0..MUTANT_ALTERNATIVES,
            rotation,
            alternative_log_e,
        );
        let start = MUTANT_ALTERNATIVES + phase * MUTANT_BLOCK_SIZE;
        let end = start + MUTANT_BLOCK_SIZE;
        mark_positions(&mut log_e, start..end, rotation, null_log_e);
        for position in start..end {
            activation_counts[hypothesis_at(position, rotation)] += 1;
        }

        let correct = e_benjamini_hochberg(&log_e, ALPHA);
        let broken = missing_family_size_e_bh(&log_e, ALPHA);
        let mut expected_broken = alternatives.clone();
        expected_broken.extend(sorted_positions(start..end, rotation));
        expected_broken.sort_unstable();
        let correct_false = correct
            .iter()
            .filter(|&&hypothesis| !is_alternative[hypothesis])
            .count();
        let broken_false = broken
            .iter()
            .filter(|&&hypothesis| !is_alternative[hypothesis])
            .count();
        result.correct_rejections += correct.len();
        result.correct_false_rejections += correct_false;
        result.mutant_rejections += broken.len();
        result.mutant_false_rejections += broken_false;
        if (correct != alternatives || broken != expected_broken)
            && result.first_shape_failure.is_none()
        {
            result.first_shape_failure = Some(phase);
        }
    }
    result.marginal_means_exact = activation_counts
        .iter()
        .enumerate()
        .all(|(hypothesis, &count)| count == u8::from(!is_alternative[hypothesis]))
        && MUTANT_NULL_E == MUTANT_PHASES as u64;
    result
}

fn campaign_identity(rotation: usize) -> ReplayIdentity {
    IdentityBuilder::new("fs-eproc-e-bh-strong-dependence-config-v1")
        .u64("family-size", FAMILY_SIZE as u64)
        .f64_bits("alpha", ALPHA)
        .u64("alpha-numerator", ALPHA_NUMERATOR as u64)
        .u64("alpha-denominator", ALPHA_DENOMINATOR as u64)
        .f64_bits("zero-log-e", ZERO_LOG_E)
        .u64("input-seed", INPUT_SEED)
        .u64("rng-kernel", u64::from(RNG_KERNEL))
        .u64("rotation-tile", u64::from(ROTATION_TILE))
        .u64("rotation-draw-index", 0)
        .u64("rotation", rotation as u64)
        .str("rotation-method", "next_below(4096)-then-add-modulo-4096")
        .u64(
            "stream-semantics-version",
            u64::from(fs_rand::STREAM_SEMANTICS_VERSION),
        )
        .str("global-law", "64-of-1024 exclusive block phases")
        .str(
            "global-phase-order",
            "phases-0-through-63-activate-corresponding-block;64-through-1023-inactive",
        )
        .u64("global-blocks", GLOBAL_BLOCKS as u64)
        .u64("global-block-size", GLOBAL_BLOCK_SIZE as u64)
        .u64("global-phases", GLOBAL_PHASES as u64)
        .f64_bits("global-active-e", GLOBAL_ACTIVE_E as f64)
        .str("perfect-law", "one all-null shock among 17 phases")
        .str(
            "perfect-phase-order",
            "phase-0-shock;phases-1-through-16-null-zero",
        )
        .u64("perfect-alternatives", PERFECT_ALTERNATIVES as u64)
        .u64("perfect-nulls", PERFECT_NULLS as u64)
        .u64("perfect-phases", PERFECT_PHASES as u64)
        .f64_bits("perfect-null-e", PERFECT_NULL_E as f64)
        .str("mutant-law", "63 mutually exclusive 64-null blocks")
        .str(
            "mutant-phase-order",
            "phase-index-activates-corresponding-null-block",
        )
        .u64("mutant-alternatives", MUTANT_ALTERNATIVES as u64)
        .u64("mutant-nulls", MUTANT_NULLS as u64)
        .u64("mutant-blocks", MUTANT_BLOCKS as u64)
        .u64("mutant-block-size", MUTANT_BLOCK_SIZE as u64)
        .u64("mutant-phases", MUTANT_PHASES as u64)
        .f64_bits("mutant-null-e", MUTANT_NULL_E as f64)
        .f64_bits("alternative-e", ALTERNATIVE_E as f64)
        .str("mutant", MUTANT_NAME)
        .str("mutant-threshold", "1/(alpha*k)")
        .str("production-threshold", "family-size/(alpha*k)")
        .str(
            "fdr-definition",
            "uniform-phase-mean(false-rejections/max(rejections,1))",
        )
        .str("zero-e-representation", "log-e=negative-infinity")
        .str("log-transform", "fs-math::det::ln")
        .str("fs-eproc-version", fs_eproc::VERSION)
        .str("fs-math-version", fs_math::VERSION)
        .str("fs-rand-version", fs_rand::VERSION)
        .str("fs-obs-version", fs_obs::VERSION)
        .finish()
}

fn result_identity(
    campaign: &ReplayIdentity,
    global: &GlobalResult,
    perfect: &PerfectResult,
    mutant: &MutantResult,
) -> ReplayIdentity {
    IdentityBuilder::new("fs-eproc-e-bh-strong-dependence-result-v1")
        .child("campaign", campaign)
        .flag("global-pass", global.pass())
        .u64("global-active-phases", global.active_phases as u64)
        .u64("global-inactive-phases", global.inactive_phases as u64)
        .u64("global-total-rejections", global.total_rejections as u64)
        .flag("global-marginal-means-exact", global.marginal_means_exact)
        .flag(
            "global-shape-failure-present",
            global.first_shape_failure.is_some(),
        )
        .u64(
            "global-first-shape-failure-phase",
            global.first_shape_failure.map_or(0, |phase| phase as u64),
        )
        .flag("perfect-pass", perfect.pass())
        .u64("perfect-shock-rejections", perfect.shock_rejections as u64)
        .u64(
            "perfect-shock-false-rejections",
            perfect.shock_false_rejections as u64,
        )
        .u64(
            "perfect-ordinary-rejections",
            perfect.ordinary_rejections as u64,
        )
        .u64(
            "perfect-ordinary-false-rejections",
            perfect.ordinary_false_rejections as u64,
        )
        .flag("perfect-marginal-means-exact", perfect.marginal_means_exact)
        .flag(
            "perfect-shape-failure-present",
            perfect.first_shape_failure.is_some(),
        )
        .u64(
            "perfect-first-shape-failure-phase",
            perfect.first_shape_failure.map_or(0, |phase| phase as u64),
        )
        .flag("mutant-caught", mutant.pass())
        .u64("correct-rejections", mutant.correct_rejections as u64)
        .u64(
            "correct-false-rejections",
            mutant.correct_false_rejections as u64,
        )
        .u64("mutant-rejections", mutant.mutant_rejections as u64)
        .u64(
            "mutant-false-rejections",
            mutant.mutant_false_rejections as u64,
        )
        .flag("mutant-marginal-means-exact", mutant.marginal_means_exact)
        .flag(
            "mutant-shape-failure-present",
            mutant.first_shape_failure.is_some(),
        )
        .u64(
            "mutant-first-shape-failure-phase",
            mutant.first_shape_failure.map_or(0, |phase| phase as u64),
        )
        .finish()
}

fn optional_phase_json(phase: Option<usize>) -> String {
    phase.map_or_else(|| "null".to_string(), |value| value.to_string())
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
    fs_obs::lint_failure_record(&event).expect("e-BH verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("e-BH verdict must use the fs-obs wire schema");
    println!("{line}");
}

fn emitter_with_receipt(
    campaign: &ReplayIdentity,
    result: &ReplayIdentity,
    rotation: usize,
    global: &GlobalResult,
    perfect: &PerfectResult,
    mutant: &MutantResult,
) -> Emitter {
    let pass = global.pass() && perfect.pass() && mutant.pass();
    let mut emitter = Emitter::new(SUITE, CASE);
    let receipt = emitter.emit(
        Severity::Info,
        EventKind::Custom {
            name: "e-bh-strong-dependence-certifier".to_string(),
            json: format!(
                "{{\"campaign_identity\":\"{}\",\"result_identity\":\"{}\",\
                 \"input_seed\":{INPUT_SEED},\"rng_kernel\":{RNG_KERNEL},\
                 \"rotation_tile\":{ROTATION_TILE},\"rotation_draw_index\":0,\
                 \"rotation\":{rotation},\"family_size\":{FAMILY_SIZE},\
                 \"alpha_bits\":\"0x{:016x}\",\"alpha_numerator\":{ALPHA_NUMERATOR},\
                 \"alpha_denominator\":{ALPHA_DENOMINATOR},\
                 \"zero_log_e\":\"negative-infinity\",\"alternative_e\":{ALTERNATIVE_E},\
                 \"stream_semantics_version\":{},\"rotation_method\":\
                 \"next_below(4096)-then-add-modulo-4096\",\"rotation_draw_count\":1,\
                 \"phase_weighting\":\"uniform-exhaustive\",\
                 \"log_transform\":\"fs-math::det::ln\",\
                 \"production_threshold\":\"family-size/(alpha*k)\",\
                 \"mutant_threshold\":\"1/(alpha*k)\",\
                 \"versions\":{{\"fs_eproc\":\"{}\",\"fs_math\":\"{}\",\
                 \"fs_rand\":\"{}\",\"fs_obs\":\"{}\"}},\
                 \"global\":{{\"phase_order\":\"0..63 active block;64..1023 inactive\",\
                 \"blocks\":{GLOBAL_BLOCKS},\"block_size\":{GLOBAL_BLOCK_SIZE},\
                 \"phases\":{GLOBAL_PHASES},\"active_e\":{GLOBAL_ACTIVE_E},\
                 \"active_phases\":{},\"inactive_phases\":{},\"total_rejections\":{},\
                 \"marginal_means_exact\":{},\"first_shape_failure\":{},\"pass\":{}}},\
                 \"perfect\":{{\"phase_order\":\"phase 0 shock;phases 1..16 null zero\",\
                 \"alternatives\":{PERFECT_ALTERNATIVES},\
                 \"nulls\":{PERFECT_NULLS},\"phases\":{PERFECT_PHASES},\
                 \"null_e\":{PERFECT_NULL_E},\"shock_rejections\":{},\
                 \"shock_false_rejections\":{},\"ordinary_rejections\":{},\
                 \"ordinary_false_rejections\":{},\"fdr_numerator\":{PERFECT_NULLS},\
                 \"fdr_denominator\":{},\"marginal_means_exact\":{},\
                 \"first_shape_failure\":{},\"pass\":{}}},\
                 \"mutant\":{{\"name\":\"{MUTANT_NAME}\",\
                 \"phase_order\":\"phase index activates corresponding null block\",\
                 \"alternatives\":{MUTANT_ALTERNATIVES},\"nulls\":{MUTANT_NULLS},\
                 \"blocks\":{MUTANT_BLOCKS},\"block_size\":{MUTANT_BLOCK_SIZE},\
                 \"phases\":{MUTANT_PHASES},\
                 \"null_e\":{MUTANT_NULL_E},\"correct_rejections\":{},\
                 \"correct_false_rejections\":{},\"mutant_rejections\":{},\
                 \"mutant_false_rejections\":{},\"marginal_means_exact\":{},\
                 \"first_shape_failure\":{},\"caught\":{}}},\"pass\":{pass}}}",
                campaign.hex(),
                result.hex(),
                ALPHA.to_bits(),
                fs_rand::STREAM_SEMANTICS_VERSION,
                fs_eproc::VERSION,
                fs_math::VERSION,
                fs_rand::VERSION,
                fs_obs::VERSION,
                global.active_phases,
                global.inactive_phases,
                global.total_rejections,
                global.marginal_means_exact,
                optional_phase_json(global.first_shape_failure),
                global.pass(),
                perfect.shock_rejections,
                perfect.shock_false_rejections,
                perfect.ordinary_rejections,
                perfect.ordinary_false_rejections,
                FAMILY_SIZE * PERFECT_PHASES,
                perfect.marginal_means_exact,
                optional_phase_json(perfect.first_shape_failure),
                perfect.pass(),
                mutant.correct_rejections,
                mutant.correct_false_rejections,
                mutant.mutant_rejections,
                mutant.mutant_false_rejections,
                mutant.marginal_means_exact,
                optional_phase_json(mutant.first_shape_failure),
                mutant.pass(),
            ),
        },
        None,
    );
    let line = receipt.to_jsonl();
    fs_obs::validate_line(&line).expect("e-BH receipt must use the fs-obs wire schema");
    println!("{line}");
    emitter
}

#[test]
fn e_bh_controls_exact_fdr_under_strong_dependence_and_catches_mutant() {
    let rotation = input_rotation();
    let global = global_null_trial(rotation);
    let perfect = perfect_dependence_trial(rotation);
    let mutant = mutant_trial(rotation);
    let campaign = campaign_identity(rotation);
    let result = result_identity(&campaign, &global, &perfect, &mutant);
    let mut emitter =
        emitter_with_receipt(&campaign, &result, rotation, &global, &perfect, &mutant);

    emit_case(
        &mut emitter,
        "global-null-block-fwer-fdr",
        global.pass(),
        format!(
            "campaign={}; result={}; 4096 nulls in 64 exclusive blocks; \
             any-rejection phases={}/{} = alpha={}/{}; total rejections={}; \
             marginal means exact={}; first shape failure={:?}",
            campaign.hex(),
            result.hex(),
            global.active_phases,
            GLOBAL_PHASES,
            ALPHA_NUMERATOR,
            ALPHA_DENOMINATOR,
            global.total_rejections,
            global.marginal_means_exact,
            global.first_shape_failure,
        ),
    );
    emit_case(
        &mut emitter,
        "mixed-perfect-dependence-fdr",
        perfect.pass(),
        format!(
            "campaign={}; result={}; exact FDR={}/({}*{}) <= {}/{}; \
             shock false/rejections={}/{}; ordinary false/rejections={}/{}; \
             marginal means exact={}; first shape failure={:?}",
            campaign.hex(),
            result.hex(),
            PERFECT_NULLS,
            FAMILY_SIZE,
            PERFECT_PHASES,
            ALPHA_NUMERATOR,
            ALPHA_DENOMINATOR,
            perfect.shock_false_rejections,
            perfect.shock_rejections,
            perfect.ordinary_false_rejections,
            perfect.ordinary_rejections,
            perfect.marginal_means_exact,
            perfect.first_shape_failure,
        ),
    );
    emit_case(
        &mut emitter,
        "missing-family-size-mutant",
        mutant.pass(),
        format!(
            "campaign={}; result={}; shipped false/rejections={}/{}; mutant \
             false/rejections={}/{} (FDP=1/2 in every phase); marginal means exact={}; \
             first shape failure={:?}",
            campaign.hex(),
            result.hex(),
            mutant.correct_false_rejections,
            mutant.correct_rejections,
            mutant.mutant_false_rejections,
            mutant.mutant_rejections,
            mutant.marginal_means_exact,
            mutant.first_shape_failure,
        ),
    );
    let pass = global.pass() && perfect.pass() && mutant.pass();
    emit_case(
        &mut emitter,
        CASE,
        pass,
        format!(
            "campaign={}; result={}; exhaustive finite laws over 4096 hypotheses; \
             global={} perfect={} mutant-caught={}; rotation={rotation}",
            campaign.hex(),
            result.hex(),
            global.pass(),
            perfect.pass(),
            mutant.pass(),
        ),
    );
    assert!(
        pass,
        "e-BH strong-dependence campaign failed: global={global:?}; \
         perfect={perfect:?}; mutant={mutant:?}"
    );
}
