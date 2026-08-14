//! Target-inaccessible ensemble executor with complete sample accounting
//! (bead frankensim-jmh21.2, core slice).
//!
//! Executes one deterministic seeded ensemble against an admitted
//! [`PredictionExecutionInput`] and produces the UNSEALED payload that the
//! bundle sealer later commits. Three properties hold by construction:
//!
//! - **Target inaccessibility**: the model callback receives ONLY logical
//!   sample coordinates and coordinate-derived seeds. No executor API
//!   accepts, stores, or forwards target outcomes; process separation is a
//!   type-level fact, not a convention.
//! - **Complete accounting**: the executor derives accounting FROM the
//!   retained per-sample outcomes; there is no path by which a caller
//!   supplies (or edits) denominators, so a dropped failure is
//!   unrepresentable rather than forbidden.
//! - **Coordinate-derived determinism**: every per-sample seed is a pure
//!   function of (stream declaration, sample index) through a versioned
//!   hash domain. Worker identity, execution order, and wall clock never
//!   reach a seed, so replay at any concurrency is bit-identical.
//!
//! No-claims: target inaccessibility protects process separation only. It
//! cannot prove the model, the uncertainty distribution, or the physical
//! prediction correct, and a sealed output of a completed run remains
//! exactly as scoreable-or-not as its referenced artifacts.

use std::collections::BTreeMap;

use fs_evidence::prediction_bundle::{
    PredictionBundleError, PredictionExecutionInput, SampleAccounting,
};
use fs_exec::Cx;

/// Versioned domain for coordinate-derived sample seeds.
pub const SAMPLE_SEED_DOMAIN: &str = "org.frankensim.fs-session.prediction-sample-seed.v1";

/// Hard ceiling on ensemble size (admission refuses beyond it).
pub const MAX_ENSEMBLE_SAMPLES: u64 = 1 << 20;

/// Logical coordinates of one sample: everything the model may see.
///
/// There is deliberately no field through which a target outcome, an
/// observation, or another sample's result could travel.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SampleCoordinates {
    /// Zero-based sample index within the requested ensemble.
    pub sample_index: u64,
    /// The admitted model rung this run executes.
    pub rung: String,
}

/// Coordinate-derived seeds, one per declared random stream.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SampleSeeds {
    seeds: BTreeMap<String, u64>,
}

impl SampleSeeds {
    /// Seed for a declared stream, if the input declared it.
    #[must_use]
    pub fn stream(&self, name: &str) -> Option<u64> {
        self.seeds.get(name).copied()
    }

    /// Declared stream names in canonical order.
    pub fn stream_names(&self) -> impl Iterator<Item = &str> {
        self.seeds.keys().map(String::as_str)
    }
}

/// What one sample produced. The executor retains every variant verbatim.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SampleOutcome {
    /// The model produced its per-sample contribution (content hashes of
    /// whatever artifacts it wrote; the executor does not interpret them).
    Succeeded {
        /// Digests of the sample's produced artifact bytes.
        artifact_hashes: Vec<fs_blake3::ContentHash>,
    },
    /// Declared-policy refusal (e.g. applicability outside the admitted
    /// domain under `ApplicabilityPolicy::Refuse`).
    Refused {
        /// Stable machine rule.
        rule: String,
    },
    /// Numerical or resource failure outside declared policy.
    Failed {
        /// Stable machine rule.
        rule: String,
    },
    /// Execution was cancelled before this sample ran (drain marker).
    Cancelled,
}

/// Typed refusals of the executor's admission and finalization boundary.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExecutorRefusal {
    /// Stable machine slug.
    pub rule: &'static str,
    /// Human diagnosis.
    pub detail: String,
}

impl core::fmt::Display for ExecutorRefusal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "{}: {}", self.rule, self.detail)
    }
}

impl std::error::Error for ExecutorRefusal {}

fn refuse(rule: &'static str, detail: impl Into<String>) -> ExecutorRefusal {
    ExecutorRefusal {
        rule,
        detail: detail.into(),
    }
}

/// Terminal disposition of one ensemble run.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunDisposition {
    /// Every requested sample reached a non-cancelled outcome.
    Completed,
    /// Cancellation was observed; every unexecuted sample carries the
    /// [`SampleOutcome::Cancelled`] drain marker.
    Cancelled {
        /// Index of the first sample that did NOT execute.
        drained_from: u64,
    },
}

/// The unsealed run payload: retained outcomes plus derived accounting.
///
/// Constructed only by [`execute_ensemble`]; the outcome vector and the
/// accounting can never disagree because the accounting is a projection.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EnsembleRun {
    input_root: fs_blake3::ContentHash,
    rung: String,
    outcomes: Vec<SampleOutcome>,
    disposition: RunDisposition,
}

impl EnsembleRun {
    /// Every retained outcome, in sample order.
    #[must_use]
    pub fn outcomes(&self) -> &[SampleOutcome] {
        &self.outcomes
    }

    /// Terminal disposition.
    #[must_use]
    pub const fn disposition(&self) -> RunDisposition {
        self.disposition
    }

    /// Identity of the sealed input this run executed.
    #[must_use]
    pub const fn input_root(&self) -> fs_blake3::ContentHash {
        self.input_root
    }

    /// Project the exact output-bundle accounting from the retained
    /// outcomes. A cancelled run has NO accounting: partial denominators
    /// are partial authority, and the sealer must never see them.
    ///
    /// # Errors
    /// Refuses on a cancelled run.
    pub fn accounting(&self) -> Result<SampleAccounting, ExecutorRefusal> {
        if let RunDisposition::Cancelled { drained_from } = self.disposition {
            return Err(refuse(
                "prediction-executor-cancelled-unscoreable",
                format!(
                    "run drained at sample {drained_from}; a cancelled ensemble \
                     has no denominators to publish"
                ),
            ));
        }
        let mut accounting = SampleAccounting {
            requested: self.outcomes.len() as u64,
            succeeded: 0,
            refused: 0,
            failed: 0,
        };
        for outcome in &self.outcomes {
            match outcome {
                SampleOutcome::Succeeded { .. } => accounting.succeeded += 1,
                SampleOutcome::Refused { .. } => accounting.refused += 1,
                SampleOutcome::Failed { .. } => accounting.failed += 1,
                SampleOutcome::Cancelled => unreachable!("completed runs hold no drain markers"),
            }
        }
        Ok(accounting)
    }
}

/// Derive the seeds for one sample from logical coordinates only.
///
/// Pure function of the input's stream declarations and the sample index:
/// `hash_domain(SAMPLE_SEED_DOMAIN, stream_domain ‖ stream_seed ‖ index)`
/// truncated to eight little-endian bytes per stream.
#[must_use]
pub fn sample_seeds(input: &PredictionExecutionInput, sample_index: u64) -> SampleSeeds {
    let mut seeds = BTreeMap::new();
    for stream in input.random_streams() {
        let mut payload = Vec::new();
        payload.extend_from_slice(stream.seed_domain.as_bytes());
        payload.push(0);
        payload.extend_from_slice(&stream.seed.to_le_bytes());
        payload.extend_from_slice(&sample_index.to_le_bytes());
        let digest = fs_blake3::hash_domain(SAMPLE_SEED_DOMAIN, &payload);
        let bytes: [u8; 8] = digest.as_bytes()[..8].try_into().expect("8 bytes");
        seeds.insert(stream.name.clone(), u64::from_le_bytes(bytes));
    }
    SampleSeeds { seeds }
}

/// Execute one deterministic ensemble.
///
/// `model` is called once per sample with coordinates and seeds only. Its
/// outcome is retained verbatim; the executor never edits, drops, or
/// reorders outcomes. Cancellation is polled BEFORE each sample; on
/// cancellation every unexecuted sample gets the drain marker and the run
/// finalizes with the cancelled disposition.
///
/// # Errors
/// Admission refusals: zero or over-cap `requested`, a `rung` outside the
/// input's admitted set (silent rung substitution is a refusal, never a
/// fallback), and a model returning the reserved [`SampleOutcome::Cancelled`]
/// variant (that marker is the executor's alone).
pub fn execute_ensemble<M>(
    cx: &Cx<'_>,
    input: &PredictionExecutionInput,
    rung: &str,
    requested: u64,
    mut model: M,
) -> Result<EnsembleRun, ExecutorRefusal>
where
    M: FnMut(&SampleCoordinates, &SampleSeeds) -> SampleOutcome,
{
    if requested == 0 || requested > MAX_ENSEMBLE_SAMPLES {
        return Err(refuse(
            "prediction-executor-ensemble-bounds",
            format!("requested must lie in 1..={MAX_ENSEMBLE_SAMPLES}, got {requested}"),
        ));
    }
    if !input
        .model_rungs()
        .allowed_rungs
        .iter()
        .any(|allowed| allowed == rung)
    {
        return Err(refuse(
            "prediction-executor-rung-not-admitted",
            format!(
                "rung {rung:?} is not in the input's admitted set; substituting \
                 another rung silently is forbidden"
            ),
        ));
    }
    let input_root = input.identity().map_err(|error: PredictionBundleError| {
        refuse(
            "prediction-executor-input-identity",
            format!("cannot derive the input root: {error}"),
        )
    })?;

    let mut outcomes = Vec::with_capacity(usize::try_from(requested).unwrap_or(0));
    let mut disposition = RunDisposition::Completed;
    for sample_index in 0..requested {
        if cx.checkpoint().is_err() {
            // Drain: mark every unexecuted sample, finalize honestly.
            for _ in sample_index..requested {
                outcomes.push(SampleOutcome::Cancelled);
            }
            disposition = RunDisposition::Cancelled {
                drained_from: sample_index,
            };
            break;
        }
        let coordinates = SampleCoordinates {
            sample_index,
            rung: rung.to_string(),
        };
        let seeds = sample_seeds(input, sample_index);
        let outcome = model(&coordinates, &seeds);
        if outcome == SampleOutcome::Cancelled {
            return Err(refuse(
                "prediction-executor-reserved-outcome",
                "the Cancelled drain marker is the executor's alone; a model \
                 refusing work must return Refused or Failed with a rule",
            ));
        }
        outcomes.push(outcome);
    }
    Ok(EnsembleRun {
        input_root,
        rung: rung.to_string(),
        outcomes,
        disposition,
    })
}
