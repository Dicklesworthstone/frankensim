//! E10.2a — the H-07 Bayesian INFERENCE ENGINE (bead
//! wf-root-guzez.11.3, de-circularized per Round-4): frozen
//! prior/likelihood contract, deterministic multi-chain sampling
//! (philox-addressed random-walk Metropolis), four LOFO conditioning
//! fits + the full-data diagnostic fit, prior-sensitivity probe,
//! sampler diagnostics (split-R̂ class) that GATE artifact signing,
//! and deterministic inference manifests.
//!
//! This engine owns NO final predictive scores — scoring is E10.2c's
//! HistoricalCampaignScoreArtifactV1, computed by a different owner
//! against artifacts this engine SIGNS. There is no score API here.

use crate::{Refusal, refuse};
use fs_blake3::hash_domain;
use fs_math::det;
use fs_rand::StreamKey;

/// Parameter-dimension cap.
pub const MAX_PARAMS: usize = 8;
/// Chain cap.
pub const MAX_CHAINS: usize = 8;
/// Per-chain sample cap.
pub const MAX_SAMPLES: usize = 65_536;
/// Observation cap.
pub const MAX_OBS: usize = 64;
/// The split-R̂ signing gate.
pub const RHAT_SIGNING_GATE: f64 = 1.05;

/// The frozen inference contract: priors + likelihood noise, hashed.
/// The digest is minted at freeze time; every fit binds to it.
#[derive(Clone, Debug, PartialEq)]
pub struct InferenceContractV1 {
    /// Parameter names (frozen order).
    pub param_names: Vec<&'static str>,
    /// Gaussian prior means.
    pub prior_mean: Vec<f64>,
    /// Gaussian prior standard deviations.
    pub prior_sd: Vec<f64>,
    /// Gaussian observation noise sd.
    pub obs_sd: f64,
}

impl InferenceContractV1 {
    /// Admit + digest the contract.
    ///
    /// # Errors
    /// `inference-contract-invalid` (dims outside [1, 8] — AT the cap
    /// admits; non-finite/non-positive scales; length mismatch).
    pub fn admit(&self) -> Result<String, Refusal> {
        let n = self.param_names.len();
        if n == 0
            || n > MAX_PARAMS
            || self.prior_mean.len() != n
            || self.prior_sd.len() != n
            || self.prior_sd.iter().any(|s| !(s.is_finite() && *s > 0.0))
            || self.prior_mean.iter().any(|m| !m.is_finite())
            || !(self.obs_sd.is_finite() && self.obs_sd > 0.0)
        {
            return Err(refuse(
                "inference-contract-invalid",
                format!("{n} params"),
                "1..=8 params, finite means, positive scales",
            ));
        }
        let mut b = Vec::new();
        for name in &self.param_names {
            b.extend_from_slice(name.as_bytes());
        }
        for v in self
            .prior_mean
            .iter()
            .chain(self.prior_sd.iter())
            .chain([self.obs_sd].iter())
        {
            b.extend_from_slice(&v.to_bits().to_le_bytes());
        }
        Ok(hash_domain("org.frankensim.wf.h07-contract.v1", &b).to_hex())
    }
}

/// One observation: an input point and the observed value.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Observation {
    /// Case id (LOFO folds leave one of these out).
    pub case: u32,
    /// Input.
    pub x: f64,
    /// Observed value.
    pub y: f64,
}

/// A fitted posterior (one conditioning fit).
#[derive(Clone, Debug, PartialEq)]
pub struct ConditioningFit {
    /// Fold label: `None` = full-data diagnostic fit; `Some(case)` =
    /// that case was LEFT OUT.
    pub left_out_case: Option<u32>,
    /// Posterior means.
    pub posterior_mean: Vec<f64>,
    /// Posterior standard deviations.
    pub posterior_sd: Vec<f64>,
    /// Split-R̂ per parameter (worst across parameters gates signing).
    pub rhat: Vec<f64>,
    /// Mean acceptance rate across chains.
    pub accept_rate: f64,
    /// Digest over the fit's samples (the sample identity).
    pub samples_digest: String,
}

impl ConditioningFit {
    /// Worst R̂ across parameters.
    #[must_use]
    pub fn worst_rhat(&self) -> f64 {
        self.rhat.iter().copied().fold(1.0, f64::max)
    }
}

/// Sampler configuration (deterministic; enters the manifest).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SamplerConfig {
    /// Root seed (philox).
    pub seed: u64,
    /// Chains.
    pub n_chains: usize,
    /// Post-warmup samples per chain.
    pub n_samples: usize,
    /// Warmup steps per chain.
    pub n_warmup: usize,
    /// Random-walk proposal scale (fraction of prior sd).
    pub proposal_frac: f64,
    /// Chain-start dispersion (prior sds from the prior mean).
    pub start_spread: f64,
}

impl SamplerConfig {
    fn admit(&self) -> Result<(), Refusal> {
        if (1..=MAX_CHAINS).contains(&self.n_chains)
            && (2..=MAX_SAMPLES).contains(&self.n_samples)
            && self.n_warmup <= MAX_SAMPLES
            && self.proposal_frac.is_finite()
            && self.proposal_frac >= 0.0
            && self.start_spread.is_finite()
        {
            Ok(())
        } else {
            Err(refuse(
                "inference-sampler-invalid",
                format!("{self:?}"),
                "chains 1..=8, samples 2..=65536",
            ))
        }
    }
}

/// The forward model: linear-in-x with the contract's parameters
/// (theta[0] + theta[1]*x + ...) — the H-07 tier-0 observable map.
/// Higher-order parameters multiply successive powers of x.
fn forward(theta: &[f64], x: f64) -> f64 {
    let mut y = 0.0;
    let mut xp = 1.0;
    for t in theta {
        y += t * xp;
        xp *= x;
    }
    y
}

fn log_post(c: &InferenceContractV1, obs: &[Observation], theta: &[f64]) -> f64 {
    let mut lp = 0.0;
    for (i, t) in theta.iter().enumerate() {
        let z = (t - c.prior_mean[i]) / c.prior_sd[i];
        lp -= 0.5 * z * z;
    }
    for o in obs {
        let z = (o.y - forward(theta, o.x)) / c.obs_sd;
        lp -= 0.5 * z * z;
    }
    lp
}

fn run_fit(
    contract: &InferenceContractV1,
    obs: &[Observation],
    cfg: &SamplerConfig,
    left_out_case: Option<u32>,
    kernel_salt: u32,
) -> ConditioningFit {
    let n = contract.param_names.len();
    let n_chains = cfg.n_chains;
    // Per-chain sample storage (post-warmup), parameter-major.
    let mut chains: Vec<Vec<Vec<f64>>> = vec![vec![Vec::with_capacity(cfg.n_samples); n]; n_chains];
    let mut accepts = 0u64;
    let mut proposals = 0u64;
    let mut digest_bytes = Vec::new();
    for chain in 0..n_chains {
        let mut s = StreamKey {
            seed: cfg.seed,
            kernel: 0x4830_3700 | kernel_salt, // "H07" + fold salt
            tile: chain as u32,
        }
        .stream();
        // Overdispersed deterministic starts.
        let mut theta: Vec<f64> = (0..n)
            .map(|i| {
                contract.prior_mean[i] + cfg.start_spread * contract.prior_sd[i] * s.next_normal()
            })
            .collect();
        let mut lp = log_post(contract, obs, &theta);
        for step in 0..(cfg.n_warmup + cfg.n_samples) {
            let mut prop = theta.clone();
            for (i, p) in prop.iter_mut().enumerate() {
                *p += cfg.proposal_frac * contract.prior_sd[i] * s.next_normal();
            }
            let lp_prop = log_post(contract, obs, &prop);
            proposals += 1;
            let u = s.next_f64().max(1e-300);
            if lp_prop - lp > det::ln(u) {
                theta = prop;
                lp = lp_prop;
                accepts += 1;
            }
            if step >= cfg.n_warmup {
                for i in 0..n {
                    chains[chain][i].push(theta[i]);
                    digest_bytes.extend_from_slice(&theta[i].to_bits().to_le_bytes());
                }
            }
        }
    }
    // Posterior moments + split-R̂ per parameter.
    let mut posterior_mean = Vec::with_capacity(n);
    let mut posterior_sd = Vec::with_capacity(n);
    let mut rhat = Vec::with_capacity(n);
    for i in 0..n {
        let all: Vec<f64> = chains.iter().flat_map(|c| c[i].iter().copied()).collect();
        let m = all.iter().sum::<f64>() / all.len() as f64;
        let v = all.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / (all.len() - 1) as f64;
        posterior_mean.push(m);
        posterior_sd.push(det::sqrt(v.max(0.0)));
        rhat.push(split_rhat(&chains, i));
    }
    ConditioningFit {
        left_out_case,
        posterior_mean,
        posterior_sd,
        rhat,
        accept_rate: accepts as f64 / proposals.max(1) as f64,
        samples_digest: hash_domain("org.frankensim.wf.h07-samples.v1", &digest_bytes).to_hex(),
    }
}

/// Split-R̂ over the chains for parameter i (each chain halved).
fn split_rhat(chains: &[Vec<Vec<f64>>], i: usize) -> f64 {
    let mut halves: Vec<&[f64]> = Vec::new();
    for c in chains {
        let s = &c[i];
        let mid = s.len() / 2;
        halves.push(&s[..mid]);
        halves.push(&s[mid..]);
    }
    let m = halves.len() as f64;
    let n = halves.iter().map(|h| h.len()).min().unwrap_or(0) as f64;
    if n < 2.0 {
        return f64::INFINITY;
    }
    let means: Vec<f64> = halves
        .iter()
        .map(|h| h.iter().sum::<f64>() / h.len() as f64)
        .collect();
    let grand = means.iter().sum::<f64>() / m;
    let b = n / (m - 1.0) * means.iter().map(|x| (x - grand) * (x - grand)).sum::<f64>();
    let w = halves
        .iter()
        .zip(means.iter())
        .map(|(h, mu)| h.iter().map(|x| (x - mu) * (x - mu)).sum::<f64>() / (h.len() as f64 - 1.0))
        .sum::<f64>()
        / m;
    if w <= 0.0 {
        return f64::INFINITY;
    }
    det::sqrt(((n - 1.0) / n * w + b / n) / w)
}

/// The signed conditioning-artifact bundle: contract digest, the
/// full-data fit, the LOFO folds, the prior-sensitivity probe, and
/// the manifest digest. Signing is GATED on diagnostics.
#[derive(Clone, Debug, PartialEq)]
pub struct ConditioningArtifactV1 {
    /// Schema id.
    pub schema: &'static str,
    /// The frozen contract's digest.
    pub contract_digest: String,
    /// Full-data diagnostic fit.
    pub full_fit: ConditioningFit,
    /// One fit per left-out case (the LOFO protocol).
    pub lofo_fits: Vec<ConditioningFit>,
    /// Posterior-mean shift under the widened prior (per parameter,
    /// REPORTED — prior sensitivity is data, not a pass/fail).
    pub prior_sensitivity_shift: Vec<f64>,
    /// Worst split-R̂ across every fit.
    pub worst_rhat: f64,
    /// The deterministic inference manifest digest.
    pub manifest_digest: String,
    /// Signature (present ONLY when diagnostics passed the gate).
    pub signature: Option<String>,
}

/// Run the full H-07 conditioning protocol and sign if diagnostics
/// permit.
///
/// # Errors
/// Contract/sampler refusals; `inference-obs-invalid` (0 or beyond
/// the cap, non-finite, or fewer distinct cases than 2 — LOFO needs
/// folds).
pub fn condition_and_sign(
    contract: &InferenceContractV1,
    obs: &[Observation],
    cfg: &SamplerConfig,
) -> Result<ConditioningArtifactV1, Refusal> {
    let contract_digest = contract.admit()?;
    cfg.admit()?;
    if obs.is_empty()
        || obs.len() > MAX_OBS
        || obs.iter().any(|o| !(o.x.is_finite() && o.y.is_finite()))
    {
        return Err(refuse(
            "inference-obs-invalid",
            format!("{} observations", obs.len()),
            "1..=64 finite observations",
        ));
    }
    let mut cases: Vec<u32> = obs.iter().map(|o| o.case).collect();
    cases.sort_unstable();
    cases.dedup();
    if cases.len() < 2 {
        return Err(refuse(
            "inference-obs-invalid",
            "LOFO needs at least two distinct cases".into(),
            "the Dec-17 protocol uses four flights",
        ));
    }
    let full_fit = run_fit(contract, obs, cfg, None, 0);
    let mut lofo_fits = Vec::with_capacity(cases.len());
    for (k, case) in cases.iter().enumerate() {
        let fold: Vec<Observation> = obs.iter().copied().filter(|o| o.case != *case).collect();
        lofo_fits.push(run_fit(contract, &fold, cfg, Some(*case), 1 + k as u32));
    }
    // Prior sensitivity: widen every prior sd 2x, refit full data.
    let widened = InferenceContractV1 {
        prior_sd: contract.prior_sd.iter().map(|s| 2.0 * s).collect(),
        ..contract.clone()
    };
    let wide_fit = run_fit(&widened, obs, cfg, None, 99);
    let prior_sensitivity_shift: Vec<f64> = full_fit
        .posterior_mean
        .iter()
        .zip(wide_fit.posterior_mean.iter())
        .map(|(a, b)| b - a)
        .collect();
    let worst_rhat = lofo_fits
        .iter()
        .chain([&full_fit])
        .map(ConditioningFit::worst_rhat)
        .fold(0.0f64, f64::max);
    let mut b = contract_digest.as_bytes().to_vec();
    b.extend_from_slice(full_fit.samples_digest.as_bytes());
    for f in &lofo_fits {
        b.extend_from_slice(f.samples_digest.as_bytes());
    }
    for v in &prior_sensitivity_shift {
        b.extend_from_slice(&v.to_bits().to_le_bytes());
    }
    b.extend_from_slice(&worst_rhat.to_bits().to_le_bytes());
    let manifest_digest = hash_domain("org.frankensim.wf.h07-manifest.v1", &b).to_hex();
    let signature = if worst_rhat <= RHAT_SIGNING_GATE {
        Some(
            hash_domain(
                "org.frankensim.wf.h07-signature.v1",
                manifest_digest.as_bytes(),
            )
            .to_hex(),
        )
    } else {
        None
    };
    Ok(ConditioningArtifactV1 {
        schema: "org.frankensim.wf.h07-conditioning.v1",
        contract_digest,
        full_fit,
        lofo_fits,
        prior_sensitivity_shift,
        worst_rhat,
        manifest_digest,
        signature,
    })
}

/// The signing gate as its own callable (consumers verify before
/// trusting an artifact).
///
/// # Errors
/// `inference-diagnostics-failed` (worst R̂ beyond the gate — the
/// artifact exists but is UNSIGNED and may not be consumed).
pub fn require_signed(artifact: &ConditioningArtifactV1) -> Result<&str, Refusal> {
    match &artifact.signature {
        Some(sig) if artifact.worst_rhat <= RHAT_SIGNING_GATE => Ok(sig),
        _ => Err(refuse(
            "inference-diagnostics-failed",
            format!(
                "worst R-hat {} beyond the {RHAT_SIGNING_GATE} signing gate",
                artifact.worst_rhat
            ),
            "longer chains or better mixing; unsigned artifacts are never consumed",
        )),
    }
}
