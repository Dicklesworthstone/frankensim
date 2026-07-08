//! fs-race conformance battery (bead 7tv.7): bitwise replay, ground-
//! truth domination, the FALSE-ELIMINATION calibration study (the
//! anytime-validity claim checked empirically against α), the measured
//! savings payoff on separated AND inseparable fields (the falsifiable
//! [M] claim), kill-registry wiring, and successive-halving brackets.

use fs_exec::KillRegistry;
use fs_race::{RaceSettings, race_field, successive_halving};

fn verdict(name: &str, pass: bool, details: &str) {
    println!("{{\"test\":\"{name}\",\"pass\":{pass},\"details\":\"{details}\"}}");
    assert!(pass, "{name}: {details}");
}

/// Deterministic noise: hash (seed, candidate, obs) → roughly N(0, 1)
/// via a 12-uniform sum (Irwin–Hall). Pure function — the racing
/// determinism contract.
fn noise(seed: u64, candidate: usize, obs: u64) -> f64 {
    let mut acc = 0.0f64;
    let mut h = seed
        .wrapping_mul(0x9e37_79b9_7f4a_7c15)
        .wrapping_add((candidate as u64) << 32)
        .wrapping_add(obs);
    for _ in 0..12 {
        h ^= h >> 33;
        h = h.wrapping_mul(0xff51_afd7_ed55_8ccd);
        h ^= h >> 33;
        #[allow(clippy::cast_precision_loss)]
        {
            acc += (h >> 11) as f64 / (1u64 << 53) as f64;
        }
    }
    acc - 6.0
}

/// race-001: bitwise replay — identical seeds give identical
/// elimination sequences, winners, and counters.
#[test]
fn race_001_replay() {
    let mus = [0.0f64, 0.5, 0.8, 1.2, 0.4, 0.9];
    let run = || {
        let kills = KillRegistry::new();
        let mut loss = |i: usize, t: u64| mus[i] + noise(0xACE, i, t);
        race_field(&mut loss, mus.len(), RaceSettings::default(), &kills)
    };
    let a = run();
    let b = run();
    verdict(
        "race-001-replay",
        a.eliminated == b.eliminated
            && a.winner == b.winner
            && a.evaluations_used == b.evaluations_used
            && a.rounds == b.rounds,
        &format!(
            "identical replays: winner {}, {} eliminations, {} evals, {} rounds",
            a.winner,
            a.eliminated.len(),
            a.evaluations_used,
            a.rounds
        ),
    );
}

/// race-002: ground-truth domination — on a well-separated field the
/// true best wins and every dominated candidate is eliminated before
/// the budget.
#[test]
fn race_002_domination() {
    let mus = [0.0f64, 1.0, 1.5, 2.0, 1.2, 0.9, 1.7, 1.3];
    let kills = KillRegistry::new();
    let mut loss = |i: usize, t: u64| mus[i] + noise(0xD0D0, i, t);
    let out = race_field(&mut loss, mus.len(), RaceSettings::default(), &kills);
    verdict(
        "race-002-domination",
        out.winner == 0 && out.survivors == vec![0] && out.eliminated.len() == 7,
        &format!(
            "true best (0) wins; eliminations {:?}; rounds {}",
            out.eliminated, out.rounds
        ),
    );
}

/// race-003: FALSE-ELIMINATION CALIBRATION — across 200 seeded
/// replays with a genuinely-best candidate, the true best is
/// eliminated no more often than α plus binomial slack (the
/// anytime-validity acceptance criterion, checked empirically).
#[test]
fn race_003_calibration() {
    let mus = [0.0f64, 0.35, 0.35, 0.5, 0.5, 0.65];
    let alpha = 0.05;
    let replays = 200u64;
    let mut false_elims = 0u32;
    for seed in 0..replays {
        let kills = KillRegistry::new();
        let mut loss = |i: usize, t: u64| mus[i] + noise(seed.wrapping_mul(0x5DEECE66D), i, t);
        let settings = RaceSettings {
            alpha,
            max_rounds: 300,
            min_rounds: 8,
        };
        let out = race_field(&mut loss, mus.len(), settings, &kills);
        if out.eliminated.iter().any(|&(_, c)| c == 0) {
            false_elims += 1;
        }
    }
    // Binomial 3σ slack around α·R.
    let expect = alpha * replays as f64;
    let slack = 3.0 * (replays as f64 * alpha * (1.0 - alpha)).sqrt();
    verdict(
        "race-003-calibration",
        f64::from(false_elims) <= expect + slack,
        &format!(
            "true best eliminated in {false_elims}/{replays} replays (alpha budget {expect:.1} + 3sigma {slack:.1}) — anytime validity holds empirically"
        ),
    );
}

/// race-004: the MEASURED payoff — a separated field saves ≥ 2× vs
/// fixed-N (the falsifiable Bet 8 claim, gated), while an INSEPARABLE
/// field runs to budget with no fake savings and the elimination
/// machinery stays quiet (α-controlled).
#[test]
fn race_004_savings() {
    let mus = [0.0f64, 1.0, 1.5, 2.0, 1.2, 0.9, 1.7, 1.3];
    let kills = KillRegistry::new();
    let mut loss = |i: usize, t: u64| mus[i] + noise(0x5A7E, i, t);
    let out = race_field(&mut loss, mus.len(), RaceSettings::default(), &kills);
    verdict(
        "race-004-separated-savings",
        out.savings() >= 2.0,
        &format!(
            "LEDGER separated field: {} evals vs fixed-N {} — savings {:.1}x (claimed 2-5x)",
            out.evaluations_used,
            out.fixed_n_equivalent,
            out.savings()
        ),
    );
    // Inseparable field: all equal means.
    let kills2 = KillRegistry::new();
    let mut loss2 = |i: usize, t: u64| noise(0xE0_01, i, t);
    let out2 = race_field(&mut loss2, 6, RaceSettings::default(), &kills2);
    verdict(
        "race-004-inseparable-honest",
        out2.savings() < 1.5 && out2.eliminated.len() <= 1,
        &format!(
            "LEDGER inseparable field: savings {:.2}x (no fake payoff), {} eliminations (alpha-controlled)",
            out2.savings(),
            out2.eliminated.len()
        ),
    );
}

/// race-005: kill wiring — eliminated candidates' registered gates
/// actually fire; survivors' gates stay clean.
#[test]
fn race_005_kill_wiring() {
    let mus = [0.0f64, 1.5, 2.0, 1.8];
    let kills = KillRegistry::new();
    let gates: Vec<_> = (0..mus.len()).map(|i| kills.register(i as u64)).collect();
    let mut loss = |i: usize, t: u64| mus[i] + noise(0x1 << 20, i, t);
    let out = race_field(&mut loss, mus.len(), RaceSettings::default(), &kills);
    let mut wiring_ok = true;
    for (i, gate) in gates.iter().enumerate() {
        let should_fire = out.eliminated.iter().any(|&(_, c)| c == i);
        if gate.is_requested() != should_fire {
            wiring_ok = false;
        }
    }
    verdict(
        "race-005-kill-wiring",
        wiring_ok && !out.eliminated.is_empty(),
        &format!(
            "gates fired exactly for the {} eliminated candidates; survivors clean",
            out.eliminated.len()
        ),
    );
}

/// race-006: successive-halving bracket — the true best survives all
/// brackets, the ledger records the halving schedule, and evaluations
/// beat fixed-N (rank-based semantics, NOT the e-guarantee —
/// documented).
#[test]
fn race_006_successive_halving() {
    let mus = [0.0f64, 0.6, 0.9, 1.2, 0.7, 1.1, 0.8, 1.4];
    let kills = KillRegistry::new();
    let mut loss = |i: usize, t: u64| mus[i] + noise(0x5_60, i, t);
    let ledger = successive_halving(&mut loss, mus.len(), 16, 2, &kills);
    let halves: Vec<usize> = ledger.brackets.iter().map(|&(_, _, after)| after).collect();
    verdict(
        "race-006-successive-halving",
        ledger.winner == 0
            && halves.windows(2).all(|w| w[1] < w[0] || w[0] == 1)
            && ledger.evaluations_used < ledger.fixed_n_equivalent,
        &format!(
            "winner {}; brackets {:?}; {} evals vs fixed-N {}",
            ledger.winner, ledger.brackets, ledger.evaluations_used, ledger.fixed_n_equivalent
        ),
    );
}
