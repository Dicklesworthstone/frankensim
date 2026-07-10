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

fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// Deterministic noise: splitmix64 keyed by the (seed, candidate, obs)
/// tuple → roughly N(0, 1) via a 12-uniform sum (Irwin–Hall). Pure
/// function — the racing determinism contract.
///
/// REPLACED under bead 7tv.7.1: the previous additive-hash version was
/// measurably NOT a null fixture — lag-1 autocorrelation ≈ 0.10 and
/// per-(seed, candidate) persistent offsets (mean² 0.047 vs the 0.017
/// an iid stream gives), i.e. candidates genuinely differed within one
/// seed. The E[e] ≤ 1 certifier below CAUGHT that (its first catch was
/// the fixture itself), and the old "calibration passes" evidence was
/// partly an artifact of it.
fn noise(seed: u64, candidate: usize, obs: u64) -> f64 {
    let mut state = mix64(
        seed ^ mix64(
            (candidate as u64)
                .wrapping_mul(0x9e37_79b9_7f4a_7c15)
                .wrapping_add(0xd1b5_4a32_d192_ed03),
        ),
    );
    state = mix64(state ^ mix64(obs.wrapping_mul(0x2545_f491_4f6c_dd1d)));
    let mut acc = 0.0f64;
    for _ in 0..12 {
        state = state.wrapping_add(0x9e37_79b9_7f4a_7c15);
        #[allow(clippy::cast_precision_loss)]
        {
            acc += (mix64(state) >> 11) as f64 / (1u64 << 53) as f64;
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

/// race-007: GLOBAL-NULL calibration (bead 7tv.7.1) — every candidate
/// identical, so ANY elimination is false. With adaptive elimination
/// and optional stopping in play (the race itself), the fixed-family
/// mixture evidence must keep the any-elimination rate within α plus
/// binomial slack. (Under the global null, e-BH's FDR equals the
/// probability of at least one rejection.)
#[test]
fn race_007_global_null() {
    let alpha = 0.05;
    let replays = 200u64;
    let mut any_elims = 0u32;
    for seed in 0..replays {
        let kills = KillRegistry::new();
        let mut loss = |i: usize, t: u64| noise(seed.wrapping_mul(0xA5A5_1234), i, t);
        let settings = RaceSettings {
            alpha,
            max_rounds: 300,
            min_rounds: 8,
        };
        let out = race_field(&mut loss, 6, settings, &kills);
        if !out.eliminated.is_empty() {
            any_elims += 1;
        }
    }
    let expect = alpha * replays as f64;
    let slack = 3.0 * (replays as f64 * alpha * (1.0 - alpha)).sqrt();
    verdict(
        "race-007-global-null",
        f64::from(any_elims) <= expect + slack,
        &format!(
            "eliminations under the global null in {any_elims}/{replays} replays \
             (alpha budget {expect:.1} + 3sigma {slack:.1})"
        ),
    );
}

/// race-008: the CERTIFIER catches the invalid max (bead 7tv.7.1).
/// Direct e-value validity check at a fixed round under the global
/// null: a genuine e-value satisfies E[e] ≤ 1 (Markov), and clamping
/// only lowers the mean, so the clamped sample mean must sit at or
/// below 1 within noise. The FIXED-FAMILY MIXTURE (shipped) passes;
/// the former SELECTIVE MAX — rebuilt here verbatim as the
/// deliberately invalid reference — is caught inflating E[e] well
/// above 1. This is the test that would have flagged the audited
/// construction.
#[test]
fn race_008_certifier_catches_max() {
    use fs_eproc::{PairwiseRace, combine_average};
    let (n, rounds, replays) = (6usize, 60u64, 300u64);
    let clamp = 20.0f64; // e-values capped — tames the tail, bias is downward only
    // Per-SEED means are the independent unit (the 6 per-candidate
    // values inside one seed share observations); empirical SEs over
    // seeds give the margins.
    let (mut mix_means, mut max_means) = (Vec::new(), Vec::new());
    for seed in 0..replays {
        let mut races: Vec<PairwiseRace> = (0..n * n).map(|_| PairwiseRace::new()).collect();
        for t in 0..rounds {
            let obs: Vec<f64> = (0..n)
                .map(|i| noise(seed.wrapping_mul(0x00BE_EF77), i, t))
                .collect();
            for i in 0..n {
                for j in (i + 1)..n {
                    races[i * n + j].observe(obs[i], obs[j]);
                    races[j * n + i].observe(obs[j], obs[i]);
                }
            }
        }
        let (mut smix, mut smax) = (0.0f64, 0.0f64);
        for i in 0..n {
            let family: Vec<f64> = (0..n)
                .filter(|&j| j != i)
                .map(|j| races[j * n + i].log_e_value())
                .collect();
            smix += combine_average(&family).exp().min(clamp);
            smax += family
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max)
                .exp()
                .min(clamp);
        }
        mix_means.push(smix / n as f64);
        max_means.push(smax / n as f64);
    }
    let stats = |v: &[f64]| {
        let m = v.iter().sum::<f64>() / v.len() as f64;
        let var = v.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / (v.len() - 1) as f64;
        (m, (var / v.len() as f64).sqrt())
    };
    let (mean_mix, se_mix) = stats(&mix_means);
    let (mean_max, se_max) = stats(&max_means);
    verdict(
        "race-008-certifier",
        mean_mix <= 1.0 + 3.0 * se_mix && mean_max > 1.0 + 3.0 * se_max,
        &format!(
            "E[e] certifier at round {rounds}: mixture {mean_mix:.3} (+-3SE {:.3} — a valid \
             e-value), selective max {mean_max:.3} (+-3SE {:.3} — CAUGHT exceeding an \
             e-value's Markov budget of 1)",
            3.0 * se_mix,
            3.0 * se_max
        ),
    );
}

/// race-009: non-finite losses are rejected STRUCTURALLY — the NaN
/// candidate is condemned that round (invalid + eliminated + kill
/// fired), its poison never reaches means or e-processes, no panic
/// anywhere, and the race's winner comes from the valid field. Same
/// contract in successive halving.
#[test]
fn race_009_nonfinite_structural() {
    let mus = [0.0f64, 0.6, 0.9, 1.2];
    let kills = KillRegistry::new();
    for i in 0..4u64 {
        let _ = kills.register(i);
    }
    let mut loss = |i: usize, t: u64| {
        if i == 2 && t == 3 {
            f64::NAN
        } else {
            mus[i] + noise(0xF1A5, i, t)
        }
    };
    let out = race_field(&mut loss, 4, RaceSettings::default(), &kills);
    let sh_kills = KillRegistry::new();
    let mut loss2 = |i: usize, t: u64| {
        if i == 1 && t == 2 {
            f64::INFINITY
        } else {
            mus[i] + noise(0xF1A6, i, t)
        }
    };
    let ledger = successive_halving(&mut loss2, 4, 8, 2, &sh_kills);
    verdict(
        "race-009-nonfinite",
        out.invalid == vec![(4, 2)]
            && out.eliminated.contains(&(4, 2))
            && out.winner != 2
            && ledger.invalid == vec![(3, 1)]
            && ledger.winner != 1,
        &format!(
            "race condemned candidate 2 at round 4 structurally (invalid {:?}, winner {}); \
             halving condemned candidate 1 at round 3 (invalid {:?}, winner {})",
            out.invalid, out.winner, ledger.invalid, ledger.winner
        ),
    );
}
