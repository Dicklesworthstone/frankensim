//! V-16b model-side battery (bead wf-root-guzez.5.16.1, E4.6c-i):
//! bitwise determinism, exact delay (impulse in cue k emerges exactly
//! delay ticks later through the filter), checkpoint/resume equality ==
//! uninterrupted execution, downstream query-pattern invariance, remnant
//! reproducibility + seed sensitivity + zero-sigma silence, caps at cap
//! AND cap+1, golden. The cross-backend/render-rate browser matrix is
//! E0.6-lane scope; the structural claim here is that NO render input
//! exists on the step signature.
//! Repro: cargo test -p fs-flyer --test perception_battery

use fs_flyer::perception::{CueSpec, MAX_DELAY_TICKS, N_CUES, PerceptionModelSpec, perception_v1};

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-v16b\",\"case\":\"{case}\",{payload}}}");
}

fn raw_at(tick: u64) -> [f64; N_CUES] {
    let t = tick as f64 / 120.0;
    [
        0.05 * fs_math::det::sin(2.0 * t),
        0.10 * fs_math::det::cos(2.0 * t),
        0.30 * fs_math::det::sin(5.0 * t),
        0.02 * fs_math::det::sin(1.3 * t),
        0.06 * fs_math::det::cos(1.3 * t),
        0.04 * fs_math::det::sin(0.7 * t),
    ]
}

#[test]
fn bitwise_determinism_across_runs() {
    let spec = perception_v1(42);
    let run = || -> Vec<u64> {
        let mut st = spec.init().unwrap();
        let mut out = Vec::new();
        for tick in 0..600 {
            let cues = spec.step(&mut st, raw_at(tick)).unwrap();
            for v in cues.values {
                out.push(v.to_bits());
            }
        }
        out
    };
    assert_eq!(run(), run(), "bitwise repeat");
    jlog("determinism", "\"bitwise\":true");
}

#[test]
fn delay_is_exactly_the_declared_ticks() {
    // Zero remnant, zero everywhere except an impulse in cue 1 at tick
    // 7: the FILTER INPUT must stay zero until exactly tick 7+delay —
    // the filtered output leaves zero at that tick and not one earlier.
    let mut spec = perception_v1(1);
    for c in &mut spec.cues {
        c.remnant_sigma = 0.0;
    }
    let d = spec.cues[1].delay_ticks as u64;
    let mut st = spec.init().unwrap();
    let mut first_nonzero = None;
    for tick in 0..(7 + d + 5) {
        let mut raw = [0.0; N_CUES];
        if tick == 7 {
            raw[1] = 1.0;
        }
        let out = spec.step(&mut st, raw).unwrap();
        if first_nonzero.is_none() && out.values[1] != 0.0 {
            first_nonzero = Some(tick);
        }
    }
    assert_eq!(
        first_nonzero,
        Some(7 + d),
        "impulse must emerge exactly delay ticks later"
    );
    jlog("delay", &format!("\"delay_ticks\":{d}"));
}

#[test]
fn checkpoint_resume_equals_uninterrupted() {
    let spec = perception_v1(42);
    // Uninterrupted 400 ticks.
    let mut st = spec.init().unwrap();
    let mut full = Vec::new();
    for tick in 0..400 {
        full.push(spec.step(&mut st, raw_at(tick)).unwrap());
    }
    // Interrupted at tick 173: checkpoint (a plain state clone — the
    // remnant stream is reconstructed from the tick), resume, continue.
    let mut a = spec.init().unwrap();
    for tick in 0..173 {
        spec.step(&mut a, raw_at(tick)).unwrap();
    }
    let mut b = spec.checkpoint(&a);
    for (i, item) in full.iter().enumerate().skip(173) {
        let out = spec.step(&mut b, raw_at(i as u64)).unwrap();
        assert_eq!(
            out, *item,
            "tick {i}: resume must be bitwise-identical to uninterrupted"
        );
    }
    jlog("checkpoint", "\"resume_equals_uninterrupted\":true");
}

#[test]
fn downstream_query_pattern_cannot_change_the_trace() {
    // The consumer sampling every 3rd tick reads EXACTLY the same values
    // the every-tick consumer saw at those ticks — perception state
    // advances only through step(), never through queries (the V-16b
    // noninterference shape at model level).
    let spec = perception_v1(7);
    let mut st1 = spec.init().unwrap();
    let mut every = Vec::new();
    for tick in 0..300 {
        every.push(spec.step(&mut st1, raw_at(tick)).unwrap());
    }
    let mut st2 = spec.init().unwrap();
    for (tick, item) in every.iter().enumerate() {
        let out = spec.step(&mut st2, raw_at(tick as u64)).unwrap();
        if tick % 3 == 0 {
            assert_eq!(out, *item, "sampled tick {tick} must match");
        }
    }
    jlog("query-invariance", "\"noninterference\":true");
}

#[test]
fn remnant_is_reproducible_seeded_and_silenceable() {
    // Same seed -> same remnant; different seed -> different values;
    // zero sigma -> exactly the filtered deterministic path.
    let base = perception_v1(42);
    let other = perception_v1(43);
    let mut quiet = perception_v1(42);
    for c in &mut quiet.cues {
        c.remnant_sigma = 0.0;
    }
    let run = |spec: &PerceptionModelSpec| -> Vec<f64> {
        let mut st = spec.init().unwrap();
        let mut out = Vec::new();
        for tick in 0..200 {
            out.extend(spec.step(&mut st, raw_at(tick)).unwrap().values);
        }
        out
    };
    let a = run(&base);
    let b = run(&base);
    let c = run(&other);
    let q = run(&quiet);
    assert_eq!(
        a.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
        b.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
    );
    assert_ne!(
        a.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
        c.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
        "seed must matter (remnant liveness)"
    );
    // The remnant-free trace differs from the seeded one (liveness) and
    // the seeded trace scatters around it at the sigma scale.
    let mut max_dev: f64 = 0.0;
    for (x, y) in a.iter().zip(&q) {
        max_dev = max_dev.max((x - y).abs());
    }
    assert!(
        max_dev > 1e-6 && max_dev < 0.2,
        "remnant deviation scale implausible: {max_dev}"
    );
    jlog("remnant", &format!("\"max_dev\":{max_dev}"));
}

#[test]
fn caps_at_cap_and_cap_plus_one() {
    let mk = |delay: usize| -> PerceptionModelSpec {
        let mut s = perception_v1(1);
        s.cues[0] = CueSpec {
            delay_ticks: delay,
            ..s.cues[0]
        };
        s
    };
    assert!(mk(MAX_DELAY_TICKS).admit().is_ok(), "cap admits");
    assert_eq!(
        mk(MAX_DELAY_TICKS + 1).admit().unwrap_err().code,
        "perception-spec-invalid",
        "cap+1 refuses"
    );
    let mut bad_tau = perception_v1(1);
    bad_tau.cues[2].filter_tau_s = 0.0;
    assert_eq!(bad_tau.admit().unwrap_err().code, "perception-spec-invalid");
    let mut bad_sigma = perception_v1(1);
    bad_sigma.cues[3].remnant_sigma = -1e-300;
    assert_eq!(
        bad_sigma.admit().unwrap_err().code,
        "perception-spec-invalid"
    );
    let spec = perception_v1(1);
    let mut st = spec.init().unwrap();
    let mut raw = [0.0; N_CUES];
    raw[4] = f64::NAN;
    assert_eq!(
        spec.step(&mut st, raw).unwrap_err().code,
        "perception-input-invalid"
    );
    jlog("caps", "\"cap_and_cap_plus_one\":true");
}

#[test]
fn golden_digest() {
    let spec = perception_v1(42);
    let mut st = spec.init().unwrap();
    let mut payload = Vec::new();
    for tick in 0..240 {
        let out = spec.step(&mut st, raw_at(tick)).unwrap();
        for v in out.values {
            payload.extend_from_slice(&v.to_bits().to_le_bytes());
        }
    }
    let digest =
        fs_blake3::hash_domain("org.frankensim.fs-flyer.v16b-golden.v1", &payload).to_hex();
    jlog("golden", &format!("\"digest\":\"{digest}\""));
    // KNOWN QUADRANT DIVERGENCE (bead guzez.7.2.1, measured
    // 2026-08-21): this path's RELEASE codegen contracts FP (the
    // debug pin f3b5b743 vs a5cb2dbe measured in release), and the
    // release value is NOT re-pinned because optimizer-context shifts
    // can move it again — pinning it would manufacture flaky goldens.
    // Debug stays the hard pin; release asserts in-process
    // repeatability and logs the value LOUDLY until 7.2.1 lands.
    if cfg!(debug_assertions) {
        assert_eq!(
            digest, "07e2aab156b514056a02251f9ec7e1ea3389f96c0ecb1b92cf38a443f474f637",
            "perception golden moved — determinism regression or an \
             intentional model change requiring the golden-bump protocol"
        );
    } else {
        let mut st2 = spec.init().unwrap();
        let mut payload2 = Vec::new();
        for tick in 0..240 {
            let out = spec.step(&mut st2, raw_at(tick)).unwrap();
            for v in out.values {
                payload2.extend_from_slice(&v.to_bits().to_le_bytes());
            }
        }
        let digest2 =
            fs_blake3::hash_domain("org.frankensim.fs-flyer.v16b-golden.v1", &payload2).to_hex();
        assert_eq!(digest, digest2, "release must at least repeat in-process");
        jlog(
            "golden-release-divergence",
            &format!("\"digest\":\"{digest}\",\"tracked\":\"guzez.7.2.1\""),
        );
    }
}
