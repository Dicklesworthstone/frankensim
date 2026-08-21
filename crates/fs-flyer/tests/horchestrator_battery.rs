//! E10.2c battery (bead wf-root-guzez.11.5): the stage pipeline on
//! the pinned synthetic fixture. The held-out-flight-never-scores-
//! itself invariant is STRUCTURALLY tested (full-data fit and
//! mismatched fold both refuse); anti-vacuity baselines (V-19) score
//! MATERIALLY worse per flight per baseline; unsigned conditioning
//! refuses; missing baselines/receipt refuse; caps at cap AND cap+1;
//! deterministic artifact digest.
//! Repro: cargo test -p fs-flyer --test horchestrator_battery

use fs_flyer::hinference::{InferenceContractV1, Observation, SamplerConfig, condition_and_sign};
use fs_flyer::horchestrator::{
    DeficientBaseline, MAX_BASELINES, orchestrate_and_score, score_held_out_flight,
};
use fs_rand::StreamKey;

fn jlog(case: &str, payload: &str) {
    println!("{{\"suite\":\"fs-flyer-horchestrator\",\"case\":\"{case}\",{payload}}}");
}

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

fn observations() -> Vec<Observation> {
    let mut s = StreamKey {
        seed: 1903,
        kernel: 0x4849_4e46,
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

fn baselines() -> Vec<DeficientBaseline> {
    vec![
        DeficientBaseline {
            label: "constant-zero-no-physics",
            mean: 0.0,
            sd: 1.0,
        },
        DeficientBaseline {
            label: "grand-mean-wide",
            mean: 2.0,
            sd: 2.0,
        },
    ]
}

const CAMPAIGN_DIGEST: &str = "a03e3ace4941c5a1d8a6c0152c697db34a5cee3220b432a1eaeb8248e9649b81";

#[test]
fn pipeline_scores_and_beats_the_deficient_baselines() {
    let art = orchestrate_and_score(
        &contract(),
        &observations(),
        &cfg(),
        CAMPAIGN_DIGEST,
        &baselines(),
    )
    .unwrap();
    assert_eq!(art.flight_scores.len(), 4, "four held-out flights");
    assert_eq!(art.campaign_merge_digest, CAMPAIGN_DIGEST);
    // V-19: the model scores MATERIALLY better than every deficient
    // baseline on EVERY flight (per-flight, per-baseline — never a
    // pooled margin only).
    for fsr in &art.flight_scores {
        assert!(fsr.model_log_score.is_finite());
        for (bi, bl) in fsr.baseline_log_scores.iter().enumerate() {
            assert!(
                fsr.model_log_score > bl + 0.5,
                "flight {} vs baseline {bi}: {} <= {} + 0.5",
                fsr.flight,
                fsr.model_log_score,
                bl
            );
        }
    }
    assert!(art.worst_margin > 0.5, "worst margin {}", art.worst_margin);
    // The prior-predictive stage is recorded and worse than held-out
    // conditioning (learning happened).
    let mean_model = art
        .flight_scores
        .iter()
        .map(|f| f.model_log_score)
        .sum::<f64>()
        / 4.0;
    assert!(
        art.prior_predictive_log_score < mean_model,
        "prior {} vs conditioned {mean_model}",
        art.prior_predictive_log_score
    );
    // Joint region: truth inside the 2-sd box per parameter.
    for (i, (lo, hi)) in art.joint_region.iter().enumerate() {
        assert!(TRUTH[i] > *lo && TRUTH[i] < *hi, "param {i}: [{lo}, {hi}]");
    }
    // Determinism.
    let again = orchestrate_and_score(
        &contract(),
        &observations(),
        &cfg(),
        CAMPAIGN_DIGEST,
        &baselines(),
    )
    .unwrap();
    assert_eq!(
        again.artifact_digest, art.artifact_digest,
        "bit-identical twice"
    );
    jlog(
        "pipeline",
        &format!(
            "\"worst_margin\":{},\"prior\":{},\"digest\":\"{}\"",
            art.worst_margin, art.prior_predictive_log_score, art.artifact_digest
        ),
    );
}

#[test]
fn holdout_invariant_is_structural() {
    let obs = observations();
    let engine_art = condition_and_sign(&contract(), &obs, &cfg()).unwrap();
    // The FULL-DATA fit (fold None) may NEVER score a flight.
    let err = score_held_out_flight(&engine_art.full_fit, 0, &obs, OBS_SD).unwrap_err();
    assert_eq!(err.code, "holdout-leakage");
    // A fold may only score ITS OWN left-out flight.
    let fold0 = engine_art
        .lofo_fits
        .iter()
        .find(|f| f.left_out_case == Some(0))
        .unwrap();
    assert!(score_held_out_flight(fold0, 0, &obs, OBS_SD).is_ok());
    let err = score_held_out_flight(fold0, 1, &obs, OBS_SD).unwrap_err();
    assert_eq!(err.code, "holdout-leakage", "mismatched fold refuses");
    // A flight with no observations refuses honestly.
    let err = score_held_out_flight(fold0, 0, &[], OBS_SD).unwrap_err();
    assert_eq!(err.code, "holdout-flight-empty");
    jlog("holdout-invariant", "\"structural\":true");
}

#[test]
fn unsigned_conditioning_refuses_the_whole_pipeline() {
    // Frozen chains -> diagnostics fail -> the orchestrator refuses
    // (the E10.2a gate propagates; no artifact is minted).
    let frozen = SamplerConfig {
        proposal_frac: 0.0,
        start_spread: 50.0,
        n_samples: 200,
        n_warmup: 0,
        ..cfg()
    };
    let err = orchestrate_and_score(
        &contract(),
        &observations(),
        &frozen,
        CAMPAIGN_DIGEST,
        &baselines(),
    )
    .unwrap_err();
    assert_eq!(err.code, "inference-diagnostics-failed");
    jlog("unsigned", &format!("\"code\":\"{}\"", err.code));
}

#[test]
fn baseline_and_receipt_refusals_with_caps() {
    // V-19: no baselines, no claim.
    let err = orchestrate_and_score(&contract(), &observations(), &cfg(), CAMPAIGN_DIGEST, &[])
        .unwrap_err();
    assert_eq!(err.code, "anti-vacuity-baseline-missing");
    // Caps: 8 baselines admits, 9 refuses.
    let mk = |n: usize| {
        (0..n)
            .map(|_| DeficientBaseline {
                label: "b",
                mean: 0.0,
                sd: 1.0,
            })
            .collect::<Vec<_>>()
    };
    assert!(
        orchestrate_and_score(
            &contract(),
            &observations(),
            &cfg(),
            CAMPAIGN_DIGEST,
            &mk(MAX_BASELINES)
        )
        .is_ok(),
        "AT cap"
    );
    assert_eq!(
        orchestrate_and_score(
            &contract(),
            &observations(),
            &cfg(),
            CAMPAIGN_DIGEST,
            &mk(MAX_BASELINES + 1)
        )
        .unwrap_err()
        .code,
        "anti-vacuity-baseline-missing"
    );
    // A score without an executed campaign receipt refuses.
    let err =
        orchestrate_and_score(&contract(), &observations(), &cfg(), "", &baselines()).unwrap_err();
    assert_eq!(err.code, "campaign-receipt-missing");
    jlog("refusals", &format!("\"max_baselines\":{MAX_BASELINES}"));
}
