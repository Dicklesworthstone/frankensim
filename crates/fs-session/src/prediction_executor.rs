//! Target-inaccessible ensemble executor with complete sample accounting
//! (bead frankensim-jmh21.2).
//!
//! Executes one deterministic seeded ensemble against an admitted
//! [`PredictionExecutionInput`] and produces the UNSEALED payload that the
//! bundle sealer later commits. Five properties hold by construction:
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
//! - **Explicit capability admission**: filesystem and network access are
//!   granted per run through [`ExecutionCapabilities`]. The grant binds
//!   checkpoint lineage and run-log identity, and resume refuses any
//!   grant other than the prefix's own, so a prefix cannot be laundered
//!   through a broader (or narrower) environment.
//! - **Adaptive diagnostics only**: the adaptive driver exposes the
//!   RETAINED prediction-side outcome prefix to the model — successes,
//!   refusals, and failures alike, so deletion is unrepresentable — and
//!   nothing else. No target type exists in this crate and the rung is
//!   fixed by the executor, so target-informed sampling and rung
//!   substitution have no API to reach.
//!
//! Forks are honest by construction: [`fork_ensemble`] records the parent
//! checkpoint identity on the child run, so a forked continuation can
//! never pose as the plain ensemble.
//!
//! No-claims: target inaccessibility protects process separation only. It
//! cannot prove the model, the uncertainty distribution, or the physical
//! prediction correct, and a sealed output of a completed run remains
//! exactly as scoreable-or-not as its referenced artifacts. Capability
//! admission is a declaration-and-audit boundary inside one process: safe
//! Rust cannot sandbox a closure against `std::fs`, so the grant makes
//! undeclared access a lineage-visible contract violation, not a physical
//! impossibility.

use std::collections::BTreeMap;

use fs_evidence::prediction_bundle::{
    PredictionBundleError, PredictionExecutionInput, SampleAccounting,
};
use fs_exec::Cx;

/// Versioned domain for coordinate-derived sample seeds.
pub const SAMPLE_SEED_DOMAIN: &str = "org.frankensim.fs-session.prediction-sample-seed.v1";

/// Hard ceiling on ensemble size (admission refuses beyond it).
pub const MAX_ENSEMBLE_SAMPLES: u64 = 1 << 20;

/// Versioned domain for executor capability-set identities.
pub const CAPABILITY_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-session.prediction-executor-capabilities.v1";

/// The execution environment granted to one ensemble run.
///
/// Compute is implicit — an executor that may not compute is absurd — so
/// the explicit axes are filesystem and network access. The grant is
/// recorded on the run, binds checkpoint lineage and run-log identity,
/// and resume requires the exact grant the executed prefix ran under.
///
/// No-claim: this is a declaration-and-audit boundary, not a sandbox.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExecutionCapabilities {
    filesystem: bool,
    network: bool,
}

impl ExecutionCapabilities {
    /// Compute only: no filesystem, no network.
    #[must_use]
    pub const fn compute_only() -> Self {
        Self {
            filesystem: false,
            network: false,
        }
    }

    /// Grant filesystem access.
    #[must_use]
    pub const fn granting_filesystem(mut self) -> Self {
        self.filesystem = true;
        self
    }

    /// Grant network access.
    #[must_use]
    pub const fn granting_network(mut self) -> Self {
        self.network = true;
        self
    }

    /// Whether filesystem access is admitted.
    #[must_use]
    pub const fn admits_filesystem(self) -> bool {
        self.filesystem
    }

    /// Whether network access is admitted.
    #[must_use]
    pub const fn admits_network(self) -> bool {
        self.network
    }

    /// Content identity of the capability set.
    #[must_use]
    pub fn identity(self) -> fs_blake3::ContentHash {
        fs_blake3::hash_domain(
            CAPABILITY_IDENTITY_DOMAIN,
            &[self.filesystem as u8, self.network as u8],
        )
    }
}

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
    /// The execution environment granted to this run.
    pub capabilities: ExecutionCapabilities,
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
/// Constructed only by [`execute_ensemble`], [`execute_ensemble_adaptive`],
/// [`resume_ensemble`], or [`fork_ensemble`]; the outcome vector and the
/// accounting can never disagree because the accounting is a projection.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EnsembleRun {
    input_root: fs_blake3::ContentHash,
    rung: String,
    outcomes: Vec<SampleOutcome>,
    disposition: RunDisposition,
    capabilities: ExecutionCapabilities,
    fork_parent: Option<fs_blake3::ContentHash>,
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

    /// The capability grant this run executed under.
    #[must_use]
    pub const fn capabilities(&self) -> ExecutionCapabilities {
        self.capabilities
    }

    /// Parent checkpoint identity for a forked continuation; `None` for a
    /// plain or resumed run.
    #[must_use]
    pub const fn fork_parent(&self) -> Option<fs_blake3::ContentHash> {
        self.fork_parent
    }

    /// Extract the resumable checkpoint from a cancelled run: the executed
    /// prefix, bound to input root, rung, and capability grant. A completed
    /// run has nothing to resume and refuses.
    ///
    /// # Errors
    /// Refuses on a completed run.
    pub fn checkpoint(&self) -> Result<EnsembleCheckpoint, ExecutorRefusal> {
        let RunDisposition::Cancelled { drained_from } = self.disposition else {
            return Err(refuse(
                "prediction-executor-nothing-to-resume",
                "a completed run has no resumable remainder",
            ));
        };
        Ok(EnsembleCheckpoint {
            input_root: self.input_root,
            rung: self.rung.clone(),
            executed: self.outcomes[..usize::try_from(drained_from).expect("bounded by cap")]
                .to_vec(),
            capabilities: self.capabilities,
        })
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
/// The `capabilities` grant is recorded on the run, its checkpoint, and
/// its log; see [`ExecutionCapabilities`].
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
    capabilities: ExecutionCapabilities,
    mut model: M,
) -> Result<EnsembleRun, ExecutorRefusal>
where
    M: FnMut(&SampleCoordinates, &SampleSeeds) -> SampleOutcome,
{
    execute_ensemble_adaptive(
        cx,
        input,
        rung,
        requested,
        capabilities,
        |coordinates, seeds, _| model(coordinates, seeds),
    )
}

/// Execute one deterministic ADAPTIVE ensemble.
///
/// Identical contract to [`execute_ensemble`] except the model also sees
/// the RETAINED prefix of prior outcomes — prediction-side diagnostics
/// only: successes, refusals, and failures alike, verbatim and in order.
/// Because the view IS the retained vector, failure deletion is
/// unrepresentable; because no target type exists in this crate and the
/// rung is fixed by the executor, target-informed sampling and rung
/// substitution have no API to reach.
pub fn execute_ensemble_adaptive<M>(
    cx: &Cx<'_>,
    input: &PredictionExecutionInput,
    rung: &str,
    requested: u64,
    capabilities: ExecutionCapabilities,
    mut model: M,
) -> Result<EnsembleRun, ExecutorRefusal>
where
    M: FnMut(&SampleCoordinates, &SampleSeeds, &[SampleOutcome]) -> SampleOutcome,
{
    let input_root = admit_run(input, rung, requested)?;
    drive_ensemble(
        cx,
        input,
        input_root,
        rung,
        requested,
        Vec::new(),
        capabilities,
        None,
        &mut model,
    )
}

/// Shared admission for one ensemble invocation: bounds, rung admission,
/// and the sealed input root.
///
/// # Errors
/// Bounds, rung-admission, and input-identity refusals verbatim.
fn admit_run(
    input: &PredictionExecutionInput,
    rung: &str,
    requested: u64,
) -> Result<fs_blake3::ContentHash, ExecutorRefusal> {
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
    input.identity().map_err(|error: PredictionBundleError| {
        refuse(
            "prediction-executor-input-identity",
            format!("cannot derive the input root: {error}"),
        )
    })
}

/// The single execution driver behind plain, adaptive, resumed, and forked
/// runs: one loop, one drain discipline, one reserved-marker rule.
///
/// # Errors
/// Only a model returning the reserved [`SampleOutcome::Cancelled`] marker
/// refuses here; everything else was admitted upstream.
#[expect(
    clippy::too_many_arguments,
    reason = "one driver keeps the four public entries from drifting apart; \
              each argument is a distinct authority the callers owe"
)]
fn drive_ensemble<M>(
    cx: &Cx<'_>,
    input: &PredictionExecutionInput,
    input_root: fs_blake3::ContentHash,
    rung: &str,
    requested: u64,
    prefix: Vec<SampleOutcome>,
    capabilities: ExecutionCapabilities,
    fork_parent: Option<fs_blake3::ContentHash>,
    model: &mut M,
) -> Result<EnsembleRun, ExecutorRefusal>
where
    M: FnMut(&SampleCoordinates, &SampleSeeds, &[SampleOutcome]) -> SampleOutcome,
{
    let mut outcomes = prefix;
    let mut disposition = RunDisposition::Completed;
    for sample_index in outcomes.len() as u64..requested {
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
            capabilities,
        };
        let seeds = sample_seeds(input, sample_index);
        let outcome = model(&coordinates, &seeds, outcomes.as_slice());
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
        capabilities,
        fork_parent,
    })
}

/// The executed prefix of a cancelled run, bound to its input root, rung,
/// and capability grant so resume lineage is attestable:
/// [`EnsembleCheckpoint::identity`] moves if ANY retained outcome, the
/// rung, the grant, or the input root changes.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EnsembleCheckpoint {
    input_root: fs_blake3::ContentHash,
    rung: String,
    executed: Vec<SampleOutcome>,
    capabilities: ExecutionCapabilities,
}

/// Versioned domain for checkpoint lineage identities. v2 binds the
/// capability grant into the lineage.
pub const CHECKPOINT_IDENTITY_DOMAIN: &str = "org.frankensim.fs-session.prediction-checkpoint.v2";

impl EnsembleCheckpoint {
    /// Number of executed samples in the prefix.
    #[must_use]
    pub fn executed_len(&self) -> u64 {
        self.executed.len() as u64
    }

    /// The capability grant the prefix was executed under.
    #[must_use]
    pub const fn capabilities(&self) -> ExecutionCapabilities {
        self.capabilities
    }

    /// Lineage identity binding input root, rung, capability grant, and
    /// every retained prefix outcome byte-exactly.
    #[must_use]
    pub fn identity(&self) -> fs_blake3::ContentHash {
        let mut payload = Vec::new();
        payload.extend_from_slice(self.input_root.as_bytes());
        payload.extend_from_slice(self.rung.as_bytes());
        payload.push(0);
        for outcome in &self.executed {
            match outcome {
                SampleOutcome::Succeeded { artifact_hashes } => {
                    payload.push(1);
                    payload.extend_from_slice(&(artifact_hashes.len() as u64).to_le_bytes());
                    for hash in artifact_hashes {
                        payload.extend_from_slice(hash.as_bytes());
                    }
                }
                SampleOutcome::Refused { rule } => {
                    payload.push(2);
                    payload.extend_from_slice(&(rule.len() as u64).to_le_bytes());
                    payload.extend_from_slice(rule.as_bytes());
                }
                SampleOutcome::Failed { rule } => {
                    payload.push(3);
                    payload.extend_from_slice(&(rule.len() as u64).to_le_bytes());
                    payload.extend_from_slice(rule.as_bytes());
                }
                SampleOutcome::Cancelled => payload.push(4),
            }
        }
        payload.extend_from_slice(self.capabilities.identity().as_bytes());
        fs_blake3::hash_domain(CHECKPOINT_IDENTITY_DOMAIN, &payload)
    }
}

/// Resume a cancelled ensemble from its checkpoint under the SAME
/// capability grant.
///
/// The retained prefix is trusted verbatim (its lineage identity is the
/// caller's attestation surface); execution continues from the first
/// unexecuted sample with the SAME coordinate-derived seeds, so an
/// uninterrupted run and a resumed run are bit-identical. Resuming is
/// continuing THE run: any grant other than the prefix's own refuses
/// rather than laundering the prefix through a different environment.
/// To continue under a DIFFERENT declared environment, use
/// [`fork_ensemble`], which records the branch honestly.
///
/// # Errors
/// Admission refusals plus: a checkpoint whose input root differs from
/// this input, whose rung differs from the requested rung, a capability
/// grant differing from the checkpoint's, whose prefix contains a drain
/// marker, or whose prefix is not shorter than `requested`.
pub fn resume_ensemble<M>(
    cx: &Cx<'_>,
    input: &PredictionExecutionInput,
    checkpoint: &EnsembleCheckpoint,
    requested: u64,
    capabilities: ExecutionCapabilities,
    mut model: M,
) -> Result<EnsembleRun, ExecutorRefusal>
where
    M: FnMut(&SampleCoordinates, &SampleSeeds) -> SampleOutcome,
{
    let input_root = admit_run(input, &checkpoint.rung, requested)?;
    if checkpoint.input_root != input_root {
        return Err(refuse(
            "prediction-executor-foreign-checkpoint",
            "checkpoint binds a different input root; resuming it here would              splice two ensembles",
        ));
    }
    if checkpoint.capabilities != capabilities {
        return Err(refuse(
            "prediction-executor-capability-mismatch",
            "resume must continue under the exact grant the executed prefix \
             ran under; narrowing or broadening both misrepresent it",
        ));
    }
    if checkpoint.executed.len() as u64 >= requested {
        return Err(refuse(
            "prediction-executor-checkpoint-bounds",
            format!(
                "checkpoint holds {} executed sample(s), not a strict prefix of {requested}",
                checkpoint.executed.len()
            ),
        ));
    }
    if checkpoint.executed.contains(&SampleOutcome::Cancelled) {
        return Err(refuse(
            "prediction-executor-checkpoint-bounds",
            "a checkpoint prefix holds executed outcomes only; a drain marker              inside it means the checkpoint was cut past the cancellation point",
        ));
    }
    drive_ensemble(
        cx,
        input,
        input_root,
        &checkpoint.rung,
        requested,
        checkpoint.executed.clone(),
        capabilities,
        None,
        &mut |coordinates, seeds, _history| model(coordinates, seeds),
    )
}

/// Fork a cancelled ensemble from its checkpoint into a NEW lineage.
///
/// Structural validation matches [`resume_ensemble`] (same input root,
/// admitted rung, strict prefix, no drain markers), but two things differ
/// BY CONSTRUCTION: any capability set may be declared for the
/// continuation because it is RECORDED on the child rather than inherited
/// silently, and the child run carries the parent checkpoint identity in
/// [`EnsembleRun::fork_parent`] — a forked continuation can never be
/// substituted for, or scored as, the plain uninterrupted ensemble.
///
/// # Errors
/// Same structural refusals as resume.
pub fn fork_ensemble<M>(
    cx: &Cx<'_>,
    input: &PredictionExecutionInput,
    checkpoint: &EnsembleCheckpoint,
    requested: u64,
    capabilities: ExecutionCapabilities,
    mut model: M,
) -> Result<EnsembleRun, ExecutorRefusal>
where
    M: FnMut(&SampleCoordinates, &SampleSeeds) -> SampleOutcome,
{
    let input_root = admit_run(input, &checkpoint.rung, requested)?;
    if checkpoint.input_root != input_root {
        return Err(refuse(
            "prediction-executor-foreign-checkpoint",
            "checkpoint binds a different input root; forking it here would              splice two ensembles",
        ));
    }
    if checkpoint.executed.len() as u64 >= requested {
        return Err(refuse(
            "prediction-executor-checkpoint-bounds",
            format!(
                "checkpoint holds {} executed sample(s), not a strict prefix of {requested}",
                checkpoint.executed.len()
            ),
        ));
    }
    if checkpoint.executed.contains(&SampleOutcome::Cancelled) {
        return Err(refuse(
            "prediction-executor-checkpoint-bounds",
            "a checkpoint prefix holds executed outcomes only; a drain marker              inside it means the checkpoint was cut past the cancellation point",
        ));
    }
    drive_ensemble(
        cx,
        input,
        input_root,
        &checkpoint.rung,
        requested,
        checkpoint.executed.clone(),
        capabilities,
        Some(checkpoint.identity()),
        &mut |coordinates, seeds, _history| model(coordinates, seeds),
    )
}

/// Versioned domain for run-log content addresses. v2 adds the capability
/// grant identity and the fork-parent binding.
pub const RUN_LOG_IDENTITY_DOMAIN: &str = "org.frankensim.fs-session.prediction-run-log.v2";
/// Run-log schema identity. v2 adds the capability grant identity and the
/// fork-parent binding.
pub const RUN_LOG_SCHEMA: &str = "frankensim.fs-session.prediction-run-log.v2";

/// Bounded, deterministic, content-addressed log of one ensemble run.
///
/// Redaction is BY CONSTRUCTION: the schema has no field for wall-clock
/// time, process identity, host identity, worker identity, or filesystem
/// paths, so none can leak into the buffered or hashed bytes. Two replays
/// of the same run produce byte-identical logs (the determinism test
/// executes this claim). v2 additionally binds the capability grant
/// identity and — for forks — the parent checkpoint identity, so runs
/// that differ only in granted environment or lineage are attestably
/// different logs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EnsembleRunLog {
    input_root: fs_blake3::ContentHash,
    rung: String,
    requested: u64,
    /// One class byte per sample: 1 succeeded, 2 refused, 3 failed,
    /// 4 cancelled-drain. Bounded by [`MAX_ENSEMBLE_SAMPLES`].
    outcome_classes: Vec<u8>,
    /// Deduplicated refusal/failure rules with occurrence counts, in
    /// canonical (rule, class) order.
    rule_counts: Vec<(String, u8, u64)>,
    /// First non-succeeded sample index, the log's "first divergence".
    first_divergence: Option<u64>,
    disposition: RunDisposition,
    /// Content identity of the capability grant in force.
    capabilities_identity: fs_blake3::ContentHash,
    /// Parent checkpoint identity for a forked continuation.
    fork_parent: Option<fs_blake3::ContentHash>,
}

impl EnsembleRunLog {
    /// Canonical log bytes (schema-versioned, deterministic).
    #[must_use]
    pub fn canonical_bytes(&self) -> Vec<u8> {
        let mut bytes = Vec::new();
        bytes.extend_from_slice(RUN_LOG_SCHEMA.as_bytes());
        bytes.push(0);
        bytes.extend_from_slice(self.input_root.as_bytes());
        bytes.extend_from_slice(self.rung.as_bytes());
        bytes.push(0);
        bytes.extend_from_slice(&self.requested.to_le_bytes());
        bytes.extend_from_slice(&self.outcome_classes);
        for (rule, class, count) in &self.rule_counts {
            bytes.push(*class);
            bytes.extend_from_slice(&count.to_le_bytes());
            bytes.extend_from_slice(rule.as_bytes());
            bytes.push(0);
        }
        match self.first_divergence {
            None => bytes.push(0),
            Some(index) => {
                bytes.push(1);
                bytes.extend_from_slice(&index.to_le_bytes());
            }
        }
        match self.disposition {
            RunDisposition::Completed => bytes.push(0),
            RunDisposition::Cancelled { drained_from } => {
                bytes.push(1);
                bytes.extend_from_slice(&drained_from.to_le_bytes());
            }
        }
        bytes.extend_from_slice(self.capabilities_identity.as_bytes());
        match self.fork_parent {
            None => bytes.push(0),
            Some(parent) => {
                bytes.push(1);
                bytes.extend_from_slice(parent.as_bytes());
            }
        }
        bytes
    }

    /// Content address of the canonical log bytes.
    #[must_use]
    pub fn identity(&self) -> fs_blake3::ContentHash {
        fs_blake3::hash_domain(RUN_LOG_IDENTITY_DOMAIN, &self.canonical_bytes())
    }

    /// First non-succeeded sample index.
    #[must_use]
    pub const fn first_divergence(&self) -> Option<u64> {
        self.first_divergence
    }

    /// Repository-relative reproduction command: replays this exact run
    /// through the executor battery harness. Contains no absolute paths.
    #[must_use]
    pub fn reproduction_command(&self) -> String {
        format!(
            "cargo test -p fs-session --test prediction_executor -- replay              # input_root={} rung={} requested={}",
            self.input_root.to_hex(),
            self.rung,
            self.requested
        )
    }
}

impl EnsembleRun {
    /// Project the bounded deterministic run log from retained outcomes.
    #[must_use]
    pub fn log(&self) -> EnsembleRunLog {
        let mut outcome_classes = Vec::with_capacity(self.outcomes.len());
        let mut rules: std::collections::BTreeMap<(String, u8), u64> =
            std::collections::BTreeMap::new();
        let mut first_divergence = None;
        for (index, outcome) in self.outcomes.iter().enumerate() {
            let class = match outcome {
                SampleOutcome::Succeeded { .. } => 1u8,
                SampleOutcome::Refused { rule } => {
                    *rules.entry((rule.clone(), 2)).or_insert(0) += 1;
                    2
                }
                SampleOutcome::Failed { rule } => {
                    *rules.entry((rule.clone(), 3)).or_insert(0) += 1;
                    3
                }
                SampleOutcome::Cancelled => 4,
            };
            if class != 1 && first_divergence.is_none() {
                first_divergence = Some(index as u64);
            }
            outcome_classes.push(class);
        }
        EnsembleRunLog {
            input_root: self.input_root,
            rung: self.rung.clone(),
            requested: self.outcomes.len() as u64,
            outcome_classes,
            rule_counts: rules
                .into_iter()
                .map(|((rule, class), count)| (rule, class, count))
                .collect(),
            first_divergence,
            disposition: self.disposition,
            capabilities_identity: self.capabilities.identity(),
            fork_parent: self.fork_parent,
        }
    }
}
