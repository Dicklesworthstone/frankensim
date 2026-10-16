//! E3.3b-i battery (bead wf-root-guzez.4.6.1) — the V-17 core, executed:
//! exact-discrete law vs the closed form, stationarity preservation
//! (variance flat over the path, no spin-up), the SUFFIX property (windows
//! of one physical path), checkpoint-addressed reconstruction equivalence
//! (bitwise), the 39-minute recurrence battery, caps, golden.
//! Repro: cargo test -p fs-atmo --test ou_battery

use fs_atmo::TICK_HZ;
use fs_atmo::ou::{MAX_OU_MODES, OuMode, STATIONARY_ANCHOR_TICK, StationaryOuPath};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-atmo-ou\",\"case\":\"{case}\",{payload}}}");
}

fn modes() -> Vec<OuMode> {
    (0..8)
        .map(|i| {
            OuMode::from_correlation_time(1.0 + 0.25 * f64::from(i), 0.5 + 0.4 * f64::from(i))
                .unwrap()
        })
        .collect()
}

#[test]
fn exact_discrete_law_matches_closed_form() {
    // ρ and σ_innov must satisfy the exact-discrete relations, and the
    // two-step composition must equal the one-double-step law:
    // ρ(2Δt) = ρ(Δt)², σ²(2Δt) = σ_st²(1 − ρ⁴).
    let m = OuMode::from_correlation_time(2.0, 1.3).unwrap();
    let rho_hand = (-1.0f64 / (TICK_HZ * 1.3)).exp();
    assert!(
        (m.rho - rho_hand).abs() < 1e-15,
        "rho {} vs {rho_hand}",
        m.rho
    );
    let innov_hand = 2.0 * (1.0 - rho_hand * rho_hand).sqrt();
    assert!((m.sigma_innov - innov_hand).abs() < 1e-14);
    // Variance recursion: propagating the stationary variance through one
    // exact step returns it EXACTLY: ρ²σ_st² + σ_innov² = σ_st².
    let recursed = m.rho * m.rho * 4.0 + m.sigma_innov * m.sigma_innov;
    assert!(
        (recursed - 4.0).abs() < 1e-12,
        "stationary fixed point violated: {recursed}"
    );
    jlog(
        "exact-law",
        &format!("\"rho\":{},\"fixed_point\":{recursed}", m.rho),
    );
}

#[test]
fn stationarity_no_spinup() {
    // The anchor draw is EXACTLY stationary, so the sample variance over
    // early and late windows of one long path must agree (no transient).
    // Fixed seed ⇒ deterministic numbers; the ±35% window band is a smoke
    // bound on ~4.6k-sample estimates, NOT a statistical claim (V-04b2
    // owns real statistics).
    let mut path = StationaryOuPath::stationary_at_anchor(11, modes()).unwrap();
    let mut early = Vec::new();
    let mut late = Vec::new();
    for tick in (STATIONARY_ANCHOR_TICK + 1)..=6000 {
        path.advance_to(tick).unwrap();
        let a0 = path.amplitudes()[0];
        if tick < STATIONARY_ANCHOR_TICK + 2400 {
            early.push(a0);
        } else if tick > 3600 {
            late.push(a0);
        }
    }
    let var = |xs: &[f64]| -> f64 {
        let mean = xs.iter().sum::<f64>() / xs.len() as f64;
        xs.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / xs.len() as f64
    };
    let sigma2 = modes()[0].sigma_stationary * modes()[0].sigma_stationary;
    let (ve, vl) = (var(&early), var(&late));
    for (name, v) in [("early", ve), ("late", vl)] {
        assert!(
            (v / sigma2 - 1.0).abs() < 0.35,
            "{name} window variance {v} vs stationary {sigma2} — spin-up transient?"
        );
    }
    jlog(
        "stationarity",
        &format!(
            "\"early_ratio\":{},\"late_ratio\":{}",
            ve / sigma2,
            vl / sigma2
        ),
    );
}

#[test]
fn suffix_property_one_physical_path() {
    // Two prerolls that both start at the FIXED anchor but stop sampling
    // at different points are the SAME path: the amplitude at tick 0 is
    // bit-identical whether we advanced straight there or paused at −1200
    // (did other work) and continued. CRN-preserving prelaunch depends on
    // exactly this.
    let mut straight = StationaryOuPath::stationary_at_anchor(23, modes()).unwrap();
    straight.advance_to(0).unwrap();
    let mut paused = StationaryOuPath::stationary_at_anchor(23, modes()).unwrap();
    paused.advance_to(-1200).unwrap();
    paused.advance_to(-77).unwrap();
    paused.advance_to(0).unwrap();
    for (a, b) in straight.amplitudes().iter().zip(paused.amplitudes()) {
        assert_eq!(
            a.to_bits(),
            b.to_bits(),
            "windows must be suffixes of ONE path"
        );
    }
    jlog("suffix", "\"one_physical_path\":true");
}

#[test]
fn checkpoint_reconstruction_is_bitwise() {
    // V-17 checkpoint equivalence: run to tick 500, checkpoint, run the
    // original to 2000; resume the checkpoint and run to 2000 — bitwise
    // identical amplitudes (innovations are counter-addressed, so the
    // resumed path consumes the same pure function).
    let mut original = StationaryOuPath::stationary_at_anchor(37, modes()).unwrap();
    original.advance_to(500).unwrap();
    let cp = original.checkpoint();
    original.advance_to(2000).unwrap();
    let mut resumed = StationaryOuPath::resume(cp, modes()).unwrap();
    resumed.advance_to(2000).unwrap();
    for (a, b) in original.amplitudes().iter().zip(resumed.amplitudes()) {
        assert_eq!(
            a.to_bits(),
            b.to_bits(),
            "resume must reproduce the path bitwise"
        );
    }
    // Checkpoint refusals: mode-count mismatch and pre-anchor tick.
    let mut bad = original.checkpoint();
    bad.amplitudes.pop();
    assert_eq!(
        StationaryOuPath::resume(bad, modes()).unwrap_err().code,
        "ou-checkpoint-invalid"
    );
    let mut early = original.checkpoint();
    early.tick = STATIONARY_ANCHOR_TICK - 1;
    assert_eq!(
        StationaryOuPath::resume(early, modes()).unwrap_err().code,
        "ou-checkpoint-invalid"
    );
    jlog("checkpoint", "\"resume\":\"bitwise\"");
}

#[test]
fn recurrence_over_the_39_minute_scenario() {
    // The 1905 Oct-5 flight: 39 min 23.8 s ≈ 283,656 ticks at 120 Hz. The
    // path must advance the full span (in checkpointed chunks under the
    // span cap) with every amplitude finite and the running variance sane.
    let mut path = StationaryOuPath::stationary_at_anchor(51, modes()).unwrap();
    let total: i64 = (39.0f64 * 60.0 + 23.8).ceil() as i64 * 120;
    let mut sum_sq = 0.0f64;
    let mut count = 0u64;
    let chunk = 40_000i64;
    let mut t = STATIONARY_ANCHOR_TICK;
    while t < total {
        t = (t + chunk).min(total);
        path.advance_to(t).unwrap();
        for a in path.amplitudes() {
            assert!(a.is_finite(), "amplitude diverged at tick {t}");
            sum_sq += a * a;
            count += 1;
        }
    }
    let sigma2 = modes()[0].sigma_stationary * modes()[0].sigma_stationary;
    let mean_sq = sum_sq / count as f64;
    // Coarse sanity only (8 sparse samples per chunk): within a decade.
    assert!(
        mean_sq > sigma2 / 10.0 && mean_sq < sigma2 * 40.0,
        "checkpoint-sampled mean-square {mean_sq} implausible vs sigma^2 {sigma2}"
    );
    assert_eq!(path.tick(), total);
    jlog(
        "recurrence",
        &format!(
            "\"ticks\":{},\"sampled_mean_sq\":{mean_sq}",
            total - STATIONARY_ANCHOR_TICK
        ),
    );
}

#[test]
fn refusals_at_cap_and_cap_plus_one() {
    let mode = OuMode::from_correlation_time(1.0, 1.0).unwrap();
    // Mode caps.
    assert!(StationaryOuPath::stationary_at_anchor(1, vec![mode; MAX_OU_MODES]).is_ok());
    assert_eq!(
        StationaryOuPath::stationary_at_anchor(1, vec![mode; MAX_OU_MODES + 1])
            .unwrap_err()
            .code,
        "mode-count-invalid"
    );
    assert_eq!(
        StationaryOuPath::stationary_at_anchor(1, vec![])
            .unwrap_err()
            .code,
        "mode-count-invalid"
    );
    // Backwards advance refuses.
    let mut path = StationaryOuPath::stationary_at_anchor(1, vec![mode]).unwrap();
    path.advance_to(10).unwrap();
    assert_eq!(path.advance_to(9).unwrap_err().code, "ou-advance-backwards");
    // Span cap at cap AND cap+1.
    let mut fresh = StationaryOuPath::stationary_at_anchor(2, vec![mode]).unwrap();
    assert!(
        fresh
            .advance_to(STATIONARY_ANCHOR_TICK + fs_atmo::ou::MAX_ADVANCE_TICKS)
            .is_ok()
    );
    let mut fresh2 = StationaryOuPath::stationary_at_anchor(3, vec![mode]).unwrap();
    assert_eq!(
        fresh2
            .advance_to(STATIONARY_ANCHOR_TICK + fs_atmo::ou::MAX_ADVANCE_TICKS + 1)
            .unwrap_err()
            .code,
        "ou-advance-span-exceeded"
    );
    // Parameter refusals.
    assert_eq!(
        OuMode::from_correlation_time(1.0, 0.0).unwrap_err().code,
        "ou-params-invalid"
    );
    jlog(
        "refusals",
        "\"gates\":\"modes cap/cap+1/0, backwards, span cap/cap+1, params\"",
    );
}

#[test]
fn ou_golden_digest() {
    // Amplitudes at ticks {−3840, −1000, 0, 5000} for the standard path:
    // exact-bit golden (measure-then-pin; golden-bump protocol).
    let mut path = StationaryOuPath::stationary_at_anchor(42, modes()).unwrap();
    let mut payload = Vec::new();
    for target in [STATIONARY_ANCHOR_TICK, -1000, 0, 5000] {
        path.advance_to(target).unwrap();
        for a in path.amplitudes() {
            payload.extend_from_slice(&a.to_bits().to_le_bytes());
        }
    }
    let digest = fs_blake3::hash_domain("org.frankensim.fs-atmo.ou-golden.v1", &payload).to_hex();
    jlog("golden", &format!("\"digest\":\"{digest}\""));
    assert_eq!(
        digest, "5b2208c63c79e67c1b34f3adeb0cd14e9c53f31c7d3bf5375e5c69fe29656ee5",
        "OU golden moved — determinism regression or an intentional \
         transition change requiring the golden-bump protocol"
    );
}
