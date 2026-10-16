//! E10.2c — the H-07 stage ORCHESTRATOR + scorer (bead
//! wf-root-guzez.11.5): prior-predictive stage → LOFO conditioning
//! fits (SIGNED by the E10.2a engine or nothing) → four held-out
//! predictive scores → full-data diagnostic →
//! HistoricalCampaignScoreArtifactV1 referencing the E10.2b campaign
//! receipts.
//!
//! Laws:
//!   - the HELD-OUT flight never scores itself: the scoring path
//!     STRUCTURALLY refuses any fit whose left-out fold is not
//!     exactly the flight being scored (leakage is a refusal, not a
//!     footnote);
//!   - only SIGNED conditioning artifacts are consumed (the E10.2a
//!     diagnostics gate propagates);
//!   - anti-vacuity baselines (V-19) are scored with the SAME proper
//!     score and their margin is carried in the artifact — a model
//!     that cannot materially beat a deficient baseline has no claim.

use crate::hinference::{
    ConditioningArtifactV1, ConditioningFit, InferenceContractV1, Observation, SamplerConfig,
    condition_and_sign, require_signed,
};
use crate::{Refusal, refuse};
use fs_blake3::hash_domain;
use fs_math::det;
use fs_rand::StreamKey;

/// Baseline cap.
pub const MAX_BASELINES: usize = 8;

/// A deficient anti-vacuity baseline (V-19): a fixed predictive mean
/// + sd applied to every observation, no physics, no conditioning.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DeficientBaseline {
    /// Baseline label.
    pub label: &'static str,
    /// Constant predictive mean.
    pub mean: f64,
    /// Constant predictive sd.
    pub sd: f64,
}

/// Proper Gaussian log score of one observation under (mean, var).
fn log_score(y: f64, mean: f64, var: f64) -> f64 {
    let v = var.max(1e-30);
    let z = y - mean;
    -0.5 * (det::ln(2.0 * core::f64::consts::PI * v) + z * z / v)
}

/// Score ONE flight's observations under a LOFO fit. The structural
/// holdout invariant lives here.
///
/// # Errors
/// `holdout-leakage` (the fit's left-out fold is not exactly this
/// flight — a full-data fit or a mismatched fold NEVER scores it).
pub fn score_held_out_flight(
    fit: &ConditioningFit,
    flight: u32,
    obs: &[Observation],
    obs_sd: f64,
) -> Result<f64, Refusal> {
    if fit.left_out_case != Some(flight) {
        return Err(refuse(
            "holdout-leakage",
            format!(
                "flight {flight} scored by a fit whose fold is {:?}",
                fit.left_out_case
            ),
            "a held-out flight is scored ONLY by the fit that left it out",
        ));
    }
    let mine: Vec<&Observation> = obs.iter().filter(|o| o.case == flight).collect();
    if mine.is_empty() {
        return Err(refuse(
            "holdout-flight-empty",
            format!("flight {flight} has no observations"),
            "score only flights that exist",
        ));
    }
    let mut total = 0.0;
    for o in &mine {
        // Posterior-predictive: linearized parameter variance + noise.
        let mut mean = 0.0;
        let mut var = obs_sd * obs_sd;
        let mut xp = 1.0;
        for (i, m) in fit.posterior_mean.iter().enumerate() {
            mean += m * xp;
            var += (fit.posterior_sd[i] * xp) * (fit.posterior_sd[i] * xp);
            xp *= o.x;
        }
        total += log_score(o.y, mean, var);
    }
    Ok(total / mine.len() as f64)
}

/// One flight's held-out score row.
#[derive(Clone, Debug, PartialEq)]
pub struct FlightScore {
    /// Flight index.
    pub flight: u32,
    /// Mean proper log score over the flight's observations.
    pub model_log_score: f64,
    /// Per-baseline mean log scores (same observations, same score).
    pub baseline_log_scores: Vec<f64>,
}

/// The final score artifact (Round-5: scorer-owned ids live HERE,
/// never on the campaign receipt).
#[derive(Clone, Debug, PartialEq)]
pub struct HistoricalCampaignScoreArtifactV1 {
    /// Schema id.
    pub schema: &'static str,
    /// The E10.2b campaign receipt this scores (merge digest).
    pub campaign_merge_digest: String,
    /// The signed conditioning manifest it consumed.
    pub conditioning_manifest_digest: String,
    /// Prior-predictive mean log score (stage 1, recorded).
    pub prior_predictive_log_score: f64,
    /// Held-out rows, flight order.
    pub flight_scores: Vec<FlightScore>,
    /// Baseline labels, column order for the rows above.
    pub baseline_labels: Vec<&'static str>,
    /// Worst (model − baseline) margin across flights × baselines.
    pub worst_margin: f64,
    /// Joint 2-sd posterior box from the full-data fit (lo, hi) per
    /// parameter — the browser artifact's region.
    pub joint_region: Vec<(f64, f64)>,
    /// Prior-sensitivity shifts (passthrough from the engine).
    pub prior_sensitivity_shift: Vec<f64>,
    /// Artifact digest.
    pub artifact_digest: String,
}

/// Run the full H-07 stage pipeline.
///
/// # Errors
/// Engine refusals propagate (incl. `inference-diagnostics-failed`
/// when the conditioning artifact is unsigned);
/// `anti-vacuity-baseline-missing` (V-19: no baselines, no claim;
/// caps at cap AND cap+1); `campaign-receipt-missing` (an empty
/// merge digest — scores must reference executed campaigns).
pub fn orchestrate_and_score(
    contract: &InferenceContractV1,
    obs: &[Observation],
    cfg: &SamplerConfig,
    campaign_merge_digest: &str,
    baselines: &[DeficientBaseline],
) -> Result<HistoricalCampaignScoreArtifactV1, Refusal> {
    if campaign_merge_digest.len() != 64 {
        return Err(refuse(
            "campaign-receipt-missing",
            "scores must reference an executed E10.2b campaign receipt".into(),
            "run the campaign first; its merge digest is the reference",
        ));
    }
    if baselines.is_empty() || baselines.len() > MAX_BASELINES {
        return Err(refuse(
            "anti-vacuity-baseline-missing",
            format!("{} baselines outside [1, {MAX_BASELINES}]", baselines.len()),
            "V-19: a score without a deficient baseline is vacuous",
        ));
    }
    if baselines
        .iter()
        .any(|b| !(b.mean.is_finite() && b.sd.is_finite() && b.sd > 0.0))
    {
        return Err(refuse(
            "anti-vacuity-baseline-missing",
            "non-finite baseline".into(),
            "finite baseline mean, positive sd",
        ));
    }
    // Stage 1: prior-predictive score (deterministic prior draws).
    let contract_digest = contract.admit()?;
    let mut s = StreamKey {
        seed: cfg.seed,
        kernel: 0x4850_5250, // "HPRP"
        tile: 0,
    }
    .stream();
    let mut prior_score = 0.0;
    for o in obs {
        // 64 deterministic prior draws -> predictive moments.
        let mut sum = 0.0;
        let mut sum2 = 0.0;
        for _ in 0..64 {
            let mut mean = 0.0;
            let mut xp = 1.0;
            for i in 0..contract.param_names.len() {
                let theta = contract.prior_mean[i] + contract.prior_sd[i] * s.next_normal();
                mean += theta * xp;
                xp *= o.x;
            }
            sum += mean;
            sum2 += mean * mean;
        }
        let m = sum / 64.0;
        let v = (sum2 / 64.0 - m * m).max(0.0) + contract.obs_sd * contract.obs_sd;
        prior_score += log_score(o.y, m, v);
    }
    let prior_predictive_log_score = prior_score / obs.len().max(1) as f64;
    // Stages 2 + 4: the engine's conditioning protocol; SIGNED or
    // nothing.
    let artifact: ConditioningArtifactV1 = condition_and_sign(contract, obs, cfg)?;
    require_signed(&artifact)?;
    // Stage 3: held-out predictive scores, one fold per flight, the
    // structural invariant enforced by score_held_out_flight.
    let mut flight_scores = Vec::new();
    let mut worst_margin = f64::INFINITY;
    for fit in &artifact.lofo_fits {
        let flight = fit.left_out_case.expect("LOFO fold names its flight");
        let model = score_held_out_flight(fit, flight, obs, contract.obs_sd)?;
        let mut baseline_scores = Vec::with_capacity(baselines.len());
        for b in baselines {
            let mine: Vec<&Observation> = obs.iter().filter(|o| o.case == flight).collect();
            let bl = mine
                .iter()
                .map(|o| log_score(o.y, b.mean, b.sd * b.sd))
                .sum::<f64>()
                / mine.len() as f64;
            worst_margin = worst_margin.min(model - bl);
            baseline_scores.push(bl);
        }
        flight_scores.push(FlightScore {
            flight,
            model_log_score: model,
            baseline_log_scores: baseline_scores,
        });
    }
    // Stage 5: the artifact.
    let joint_region: Vec<(f64, f64)> = artifact
        .full_fit
        .posterior_mean
        .iter()
        .zip(artifact.full_fit.posterior_sd.iter())
        .map(|(m, sd)| (m - 2.0 * sd, m + 2.0 * sd))
        .collect();
    let mut b = contract_digest.as_bytes().to_vec();
    b.extend_from_slice(campaign_merge_digest.as_bytes());
    b.extend_from_slice(artifact.manifest_digest.as_bytes());
    b.extend_from_slice(&prior_predictive_log_score.to_bits().to_le_bytes());
    for fsr in &flight_scores {
        b.extend_from_slice(&fsr.flight.to_le_bytes());
        b.extend_from_slice(&fsr.model_log_score.to_bits().to_le_bytes());
        for v in &fsr.baseline_log_scores {
            b.extend_from_slice(&v.to_bits().to_le_bytes());
        }
    }
    let artifact_digest = hash_domain("org.frankensim.wf.h-campaign-score.v1", &b).to_hex();
    Ok(HistoricalCampaignScoreArtifactV1 {
        schema: "org.frankensim.wf.h-campaign-score.v1",
        campaign_merge_digest: campaign_merge_digest.to_string(),
        conditioning_manifest_digest: artifact.manifest_digest,
        prior_predictive_log_score,
        flight_scores,
        baseline_labels: baselines.iter().map(|b| b.label).collect(),
        worst_margin,
        joint_region,
        prior_sensitivity_shift: artifact.prior_sensitivity_shift,
        artifact_digest,
    })
}
