//! E10.2a battery (bead wf-root-guzez.11.3): the H-07 inference
//! engine on PINNED synthetic-truth fixtures. Posterior recovers the
//! planted truth (per-parameter oracles); four LOFO folds each cover
//! their left-out flight (posterior-predictive z per fold, never a
//! pooled score); the diagnostics gate SIGNS healthy artifacts and
//! the unconverged falsifier is EXECUTED (frozen chains, huge
//! dispersion -> unsigned -> require_signed refuses); prior
//! sensitivity REPORTED; deterministic manifests bitwise twice; caps
//! at cap AND cap+1.
//! Repro: cargo test -p fs-flyer --test hinference_battery

use fs_flyer::hinference::{
    InferenceContractV1, MAX_OBS, MAX_PARAMS, Observation, RHAT_SIGNING_GATE, SamplerConfig,
    condition_and_sign, require_signed,
};
use fs_rand::StreamKey;

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-hinference\",\"case\":\"{case}\",{payload}}}");
}

/// Planted truth for the synthetic-recovery fixture.
const TRUTH: [f64; 2] = [2.0, -0.7];
const OBS_SD: f64 = 0.05;

fn contract() -> InferenceContractV1 {
    InferenceContractV1 {
        param_names: vec!["intercept", "slope"],
        prior_mean: vec![0.0, 0.0],
        prior_sd: vec![5.0, 5.0],
        obs_sd: OBS_SD,
    }
}

/// Four synthetic "flights", three observations each, philox noise
/// at the pinned seed (never wall-clock, never thread order).
fn observations() -> Vec<Observation> {
    let mut s = StreamKey {
        seed: 1903,
        kernel: 0x4849_4e46, // "HINF"
        tile: 0,
    }
    .stream();
    let mut obs = Vec::new();
    for case in 0..4u32 {
        for k in 0..3 {
            let x = -1.0 + 0.6 * f64::from(case) + 0.2 * f64::from(k);
            let y = TRUTH[0] + TRUTH[1] * x + OBS_SD * s.next_normal();
            obs.push(Observation { case, x, y });
        }
    }
    obs
}

fn cfg() -> SamplerConfig {
    SamplerConfig {
        seed: 17,
        n_chains: 4,
        n_samples: 4000,
        n_warmup: 1000,
        proposal_frac: 0.02,
        start_spread: 1.0,
    }
}

#[test]
fn synthetic_truth_recovery_and_signing() {
    let art = condition_and_sign(&contract(), &observations(), &cfg()).unwrap();
    // Per-parameter recovery oracles: truth within 4 posterior sds
    // AND within an absolute band (the pinned fixture is informative).
    for (i, truth) in TRUTH.iter().enumerate() {
        let m = art.full_fit.posterior_mean[i];
        let sd = art.full_fit.posterior_sd[i];
        assert!(
            (m - truth).abs() < 4.0 * sd.max(1e-6),
            "param {i}: {m} vs {truth} (sd {sd})"
        );
        assert!(
            (m - truth).abs() < 0.1,
            "param {i} absolute: {m} vs {truth}"
        );
    }
    // Healthy diagnostics -> SIGNED.
    assert!(
        art.worst_rhat < RHAT_SIGNING_GATE,
        "rhat {}",
        art.worst_rhat
    );
    let sig = require_signed(&art).unwrap();
    assert_eq!(sig.len(), 64);
    // Prior sensitivity is REPORTED and small for informative data.
    for (i, shift) in art.prior_sensitivity_shift.iter().enumerate() {
        assert!(shift.abs() < 0.05, "prior shift {i}: {shift}");
    }
    // Deterministic manifests: bitwise twice.
    let again = condition_and_sign(&contract(), &observations(), &cfg()).unwrap();
    assert_eq!(
        again.manifest_digest, art.manifest_digest,
        "bit-identical twice"
    );
    jlog(
        "recovery",
        &format!(
            "\"mean\":[{},{}],\"worst_rhat\":{},\"manifest\":\"{}\"",
            art.full_fit.posterior_mean[0],
            art.full_fit.posterior_mean[1],
            art.worst_rhat,
            art.manifest_digest
        ),
    );
}

#[test]
fn lofo_folds_cover_their_left_out_flights_per_fold() {
    let obs = observations();
    let art = condition_and_sign(&contract(), &obs, &cfg()).unwrap();
    assert_eq!(
        art.lofo_fits.len(),
        4,
        "four LOFO folds (the Dec-17 protocol)"
    );
    for fit in &art.lofo_fits {
        let case = fit.left_out_case.expect("a LOFO fold names its fold");
        // Posterior-predictive coverage of EVERY left-out observation
        // in this fold (per-item, never pooled).
        for o in obs.iter().filter(|o| o.case == case) {
            let pred = fit.posterior_mean[0] + fit.posterior_mean[1] * o.x;
            // Predictive sd: obs noise + linearized parameter noise.
            let var = OBS_SD * OBS_SD
                + fit.posterior_sd[0] * fit.posterior_sd[0]
                + (o.x * fit.posterior_sd[1]) * (o.x * fit.posterior_sd[1]);
            let z = (o.y - pred).abs() / var.sqrt();
            assert!(z < 4.0, "fold {case}, x {}: z {z}", o.x);
        }
    }
    // Fold digests are distinct (each fold genuinely refit).
    let mut digests: Vec<&str> = art
        .lofo_fits
        .iter()
        .map(|f| f.samples_digest.as_str())
        .collect();
    digests.push(&art.full_fit.samples_digest);
    let before = digests.len();
    digests.sort_unstable();
    digests.dedup();
    assert_eq!(digests.len(), before, "every fit has its own samples");
    jlog("lofo", "\"folds\":4");
}

#[test]
fn unconverged_falsifier_blocks_signing() {
    // FALSIFIER: frozen chains (zero proposal scale) started far
    // apart never mix — split-R-hat explodes and signing refuses.
    let frozen = SamplerConfig {
        proposal_frac: 0.0,
        start_spread: 50.0,
        n_samples: 200,
        n_warmup: 0,
        ..cfg()
    };
    let art = condition_and_sign(&contract(), &observations(), &frozen).unwrap();
    assert!(
        art.worst_rhat > RHAT_SIGNING_GATE,
        "frozen chains must fail diagnostics: {}",
        art.worst_rhat
    );
    assert!(art.signature.is_none(), "unsigned");
    let err = require_signed(&art).unwrap_err();
    assert_eq!(err.code, "inference-diagnostics-failed");
    jlog(
        "falsifier",
        &format!("\"rhat\":{},\"code\":\"{}\"", art.worst_rhat, err.code),
    );
}

#[test]
fn contract_freeze_and_caps() {
    // The contract digest moves when a prior moves (frozen contract:
    // any post-hoc prior edit is visible to every consumer).
    let base = contract().admit().unwrap();
    let mut widened = contract();
    widened.prior_sd[0] = 6.0;
    assert_ne!(
        base,
        widened.admit().unwrap(),
        "prior edits move the digest"
    );
    // Param caps: 8 admits, 9 refuses.
    let mk = |n: usize| InferenceContractV1 {
        param_names: (0..n).map(|_| "p").collect(),
        prior_mean: vec![0.0; n],
        prior_sd: vec![1.0; n],
        obs_sd: 0.1,
    };
    assert!(mk(MAX_PARAMS).admit().is_ok(), "AT cap");
    assert_eq!(
        mk(MAX_PARAMS + 1).admit().unwrap_err().code,
        "inference-contract-invalid"
    );
    assert_eq!(
        mk(0).admit().unwrap_err().code,
        "inference-contract-invalid"
    );
    // Sampler caps.
    let bad = SamplerConfig {
        n_chains: 9,
        ..cfg()
    };
    assert_eq!(
        condition_and_sign(&contract(), &observations(), &bad)
            .unwrap_err()
            .code,
        "inference-sampler-invalid"
    );
    // Observation caps: 64 admits, 65 refuses; single-case refuses
    // (LOFO needs folds).
    let mut many: Vec<Observation> = (0..MAX_OBS)
        .map(|i| Observation {
            case: (i % 4) as u32,
            x: i as f64 * 0.01,
            y: 2.0,
        })
        .collect();
    let quick = SamplerConfig {
        n_samples: 50,
        n_warmup: 10,
        ..cfg()
    };
    assert!(
        condition_and_sign(&contract(), &many, &quick).is_ok(),
        "AT cap"
    );
    many.push(Observation {
        case: 0,
        x: 0.0,
        y: 0.0,
    });
    assert_eq!(
        condition_and_sign(&contract(), &many, &quick)
            .unwrap_err()
            .code,
        "inference-obs-invalid"
    );
    let single = vec![
        Observation {
            case: 7,
            x: 0.0,
            y: 1.0
        };
        3
    ];
    assert_eq!(
        condition_and_sign(&contract(), &single, &quick)
            .unwrap_err()
            .code,
        "inference-obs-invalid"
    );
    jlog(
        "caps",
        &format!("\"max_params\":{MAX_PARAMS},\"max_obs\":{MAX_OBS}"),
    );
}
