//! Bead 7tv.21.17: statistical validity of optimizer comparisons.
//!
//! fs-eproc's `PairwiseRace` supplies anytime-valid e-process racing; this
//! battery proves the guarantee HOLDS on real optimizer loss streams:
//!
//! - NULL: races between IDENTICAL DE configurations (differing only in
//!   seed) monitored continuously under ADVERSARIAL stopping (a race
//!   counts as a false elimination if it rejects at ANY observation) must
//!   keep the either-direction elimination rate within the 2α e-process
//!   bound plus 3σ binomial slack.
//! - POWER: a budget-starved arm must be eliminated within the
//!   observation budget in at least 90% of races.
//! - REPLAY: identical seeds reproduce the exact log-e trajectory
//!   bit-for-bit.
//!
//! Losses are min(ln(1 + f_final), 2) on canonical Rastrigin-2 — a
//! documented monotone-capped transform: the cap tightens the loss span
//! (better bets per observation), reads everything above f ≈ 6.4 as
//! equally bad, and is symmetric across arms, so null exchangeability is
//! untouched.

use fs_eproc::{LossSpan, PairwiseRace};
use fsci_opt::{DifferentialEvolutionOptions, differential_evolution};

const SUITE: &str = "fs-ascent/eproc-race-validity";

fn verdict(name: &str, pass: bool, details: &str) {
    // Canonical fs-obs verdict (51bbx suite adoption): linted, wire-validated,
    // emitted BEFORE the assertion so failures reproduce from the log alone.
    let mut emitter = fs_obs::Emitter::new(SUITE, name);
    let event = emitter.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: name.to_string(),
            pass,
            detail: details.to_string(),
            seed: 0, // deterministic fixed-input batteries; embedded seeds live in detail
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("verdict must satisfy the failure-record lint");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("verdict must use the fs-obs wire schema");
    println!("{line}");
    assert!(pass, "{name}: {details}");
}

fn rastrigin2(x: &[f64]) -> f64 {
    20.0 + x
        .iter()
        .map(|&v| v * v - 10.0 * (2.0 * core::f64::consts::PI * v).cos())
        .sum::<f64>()
}

/// One DE run's loss under the documented transform.
fn de_loss(seed: u64, maxiter: usize) -> f64 {
    let res = differential_evolution(
        rastrigin2,
        &[(-5.12, 5.12), (-5.12, 5.12)],
        DifferentialEvolutionOptions {
            seed: Some(seed),
            maxiter,
            ..DifferentialEvolutionOptions::default()
        },
    )
    .expect("DE runs");
    (1.0 + rastrigin2(&res.x)).ln().min(LOSS_CAP)
}

const ALPHA: f64 = 0.05;
const LOSS_CAP: f64 = 2.0;
const OBS_PER_RACE: usize = 25;
const DEFAULT_MAXITER: usize = 1000;

/// Run one race under adversarial stopping: reject the moment either
/// direction's e-value crosses 1/α at ANY observation. Returns
/// (a_eliminated_b, b_eliminated_a, final log-e values).
fn adversarial_race(seed_base: u64, maxiter_a: usize, maxiter_b: usize) -> (bool, bool, f64, f64) {
    let span = LossSpan::new(LOSS_CAP).expect("valid span");
    let mut forward = PairwiseRace::new(span);
    let mut reverse = PairwiseRace::new(span);
    let mut a_wins = false;
    let mut b_wins = false;
    for obs in 0..OBS_PER_RACE {
        let loss_a = de_loss(seed_base + 2 * obs as u64, maxiter_a);
        let loss_b = de_loss(seed_base + 2 * obs as u64 + 1, maxiter_b);
        forward.observe(loss_a, loss_b).expect("in-span losses");
        reverse.observe(loss_b, loss_a).expect("in-span losses");
        a_wins |= forward.a_beats_b(ALPHA);
        b_wins |= reverse.a_beats_b(ALPHA);
    }
    (a_wins, b_wins, forward.log_e_value(), reverse.log_e_value())
}

#[test]
fn null_races_hold_the_anytime_false_elimination_bound() {
    // M identical-config races. The e-process guarantee is
    // P(sup E >= 1/alpha) <= alpha per direction, so either-direction
    // eliminations are bounded by 2*alpha even under adversarial
    // stopping. Gate: rate <= 2*alpha + 3*sigma (binomial), seeds fixed
    // so the gate is deterministic.
    const M: usize = 100;
    let mut eliminations = 0usize;
    for race in 0..M {
        let (a, b, _, _) = adversarial_race(
            10_000 + race as u64 * 1_000,
            DEFAULT_MAXITER,
            DEFAULT_MAXITER,
        );
        if a || b {
            eliminations += 1;
        }
    }
    let bound_rate = 2.0 * ALPHA;
    let slack = 3.0 * (bound_rate * (1.0 - bound_rate) / M as f64).sqrt();
    let gate = ((bound_rate + slack) * M as f64).floor() as usize;
    assert!(
        eliminations <= gate,
        "null false eliminations {eliminations}/{M} exceed the anytime bound gate {gate}"
    );
    verdict(
        "7tv21-eproc-null",
        true,
        &format!(
            "{eliminations}/{M} null races eliminated under adversarial stopping \
             (anytime gate {gate} = 2a+3s at a={ALPHA})"
        ),
    );
}

#[test]
fn budget_starved_arms_are_eliminated_with_power() {
    // Arm B gets maxiter=2 (budget-starved: DE barely improves on its
    // initial population), a genuinely worse optimizer configuration.
    const M: usize = 20;
    let mut eliminated = 0usize;
    let mut wrong_direction = 0usize;
    for race in 0..M {
        let (a, b, _, _) = adversarial_race(500_000 + race as u64 * 1_000, DEFAULT_MAXITER, 2);
        if a {
            eliminated += 1;
        }
        if b {
            wrong_direction += 1;
        }
    }
    assert_eq!(
        wrong_direction, 0,
        "the starved arm must never eliminate the healthy arm"
    );
    assert!(
        eliminated * 10 >= M * 9,
        "power too low: {eliminated}/{M} starved arms eliminated"
    );
    verdict(
        "7tv21-eproc-power",
        true,
        &format!(
            "{eliminated}/{M} budget-starved arms eliminated within {OBS_PER_RACE} \
             observations (gate 90%); 0 wrong-direction eliminations"
        ),
    );
}

#[test]
fn race_trajectories_replay_bit_for_bit() {
    let first = adversarial_race(777_000, DEFAULT_MAXITER, 2);
    let second = adversarial_race(777_000, DEFAULT_MAXITER, 2);
    assert_eq!(
        first.2.to_bits(),
        second.2.to_bits(),
        "forward log-e trajectory must replay exactly"
    );
    assert_eq!(
        first.3.to_bits(),
        second.3.to_bits(),
        "reverse log-e trajectory must replay exactly"
    );
    assert_eq!((first.0, first.1), (second.0, second.1));
    verdict(
        "7tv21-eproc-replay",
        true,
        "identical seeds reproduce both directions' log-e values bit-for-bit",
    );
}
