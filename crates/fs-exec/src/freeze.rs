//! Typed single-use snapshot freeze transactions (bead frankensim-sj31i.52.5.4.1).
//!
//! [`crate::cx::DrainFinalizeReport`] proves that cancellation was requested
//! and every registered worker drained, but the report is `Copy` and its
//! semantic identity contains only version/run/worker counts. Any path that
//! lets a caller *pair* such a report with caller-chosen pause labels and a
//! state value provides caller-relative consistency only: a copied report can
//! be bound to arbitrary labels later, and nothing proves the serialized state
//! was atomically frozen by the session owning the paused run.
//!
//! This module closes that gap with ownership instead of flags:
//!
//! 1. The session owner holds one private [`SnapshotFreezeRegistry`] per
//!    (owner, session, solver-instance generation). The registry mints and
//!    burns every permit and performs exactly one transaction in its life.
//! 2. [`SnapshotFreezeRegistry::begin_freeze`] closes mutation admission;
//!    every later transaction on this registry refuses forever.
//! 3. [`SnapshotFreezeRequest::freeze`] drains the exact run through
//!    executor-owned [`DrainTracker::finalize`] — callers never supply the
//!    report — immutably captures the state, and encodes the canonical payload
//!    **exactly once**. Payload bytes live inside the minted
//!    [`SnapshotFreezePermit`]; sealing consumes those exact bytes, so no
//!    TOCTOU window exists between freeze commitment and envelope encoding.
//! 4. [`SnapshotFreezePermit`] is linear: no `Clone`, no `Copy`, no
//!    serialization, private fields, no public constructor outside minting.
//!    Sealing burns the permit nonce in the registry before any work
//!    (burn-before-call). Panic, cancellation, and every seal refusal poison
//!    the registry terminal: the same identity can never publish another
//!    state afterwards. There is no silent retry; a fresh freeze requires a
//!    new solver-instance generation and therefore a new registry.
//! 5. The resulting [`SnapshotFreezeReceipt`] states precisely what happened:
//!    a live permit is [`FreezeDisposition::BytesPrepared`], a sealed outcome
//!    is [`FreezeDisposition::BytesCommitted`]. Publication and activation
//!    remain downstream obligations recorded elsewhere; this module claims
//!    neither.
//!
//! # Uniqueness contract
//!
//! [`FreezeOwnerBinding::instance_nonce`] must be unique per constructed
//! registry (the session layer persists it next to the solver-instance
//! identity). The full-width registry identity derives from it; two live
//! registries must never share one nonce value.
//!
//! # No-claim boundary
//!
//! A committed receipt proves: the registry owning (owner, session,
//! solver-instance generation) admitted exactly one transaction, the named run
//! drained through the executor drain tracker, the exact captured state was
//! encoded once, and the sealed envelope embeds those identities. It does not
//! authenticate the owner, prove wall-clock atomicity beyond Rust moves, or
//! validate physics. Header declarations remain replayable data and can never
//! reconstruct a permit.
//!
//! # Compile-fail misuse battery
//!
//! Forged construction and duplication are compile errors:
//!
//! ```compile_fail,E0277
//! use fs_exec::freeze::SnapshotFreezePermit;
//!
//! fn duplicate<S>(permit: &SnapshotFreezePermit<'_, S>) -> SnapshotFreezePermit<'_, S> {
//!     permit.clone() // ERROR: `Clone` is deliberately not implemented.
//! }
//! ```
//!
//! ```compile_fail,E0451
//! use fs_exec::freeze::SnapshotFreezePermit;
//!
//! // Private fields refuse struct literals outside the minting path.
//! let forged: SnapshotFreezePermit<'_, u8> = SnapshotFreezePermit {
//!     core: unimplemented!(),
//! };
//! ```

use core::fmt;

use crate::cx::{DrainFinalizeError, DrainFinalizeReport, DrainTracker, RunId};
use crate::solver::snapshot_v2::{
    ExpectedResumeContextV2, PausedSnapshotBoundaryV2, SealedSnapshotV2, SnapshotAlgorithmIdV2,
    SnapshotBudgetStateIdV2, SnapshotDeterminismV2, SnapshotExecutionFingerprintIdV2,
    SnapshotLimitsV2, SnapshotPauseRequestIdV2, SnapshotProblemIdV2, SnapshotProvenanceIdV2,
    SnapshotRngCounterIdV2, SnapshotV2Error,
};
use crate::solver::{SolverStateV2, snapshot_v2};

/// Identity domain of one freeze registry instance.
pub const FREEZE_REGISTRY_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-exec.solver-freeze-registry.v2";
/// Identity domain of one single-use freeze/burn nonce.
pub const FREEZE_NONCE_IDENTITY_DOMAIN: &str = "org.frankensim.fs-exec.solver-freeze-nonce.v2";
/// Identity domain of the pre-seal transaction commitment.
pub const FREEZE_COMMITMENT_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-exec.solver-freeze-commitment.v2";
/// Identity domain of committed freeze receipts.
pub const FREEZE_RECEIPT_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-exec.solver-freeze-receipt.v2";

/// Owner/session/solver-instance binding declared by the session owner.
///
/// All fields are explicit; the registry derives its full-width identity from
/// them plus a nonzero caller-persisted instance nonce.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FreezeOwnerBinding {
    /// Opaque full-width owner identity (e.g. derived from governor id).
    pub owner: [u8; 32],
    /// Session the paused run belongs to.
    pub session: u64,
    /// Solver-instance generation. One registry serves exactly one generation.
    pub solver_instance_generation: u64,
    /// Caller-persisted nonzero nonce distinguishing registries that would
    /// otherwise share a binding. Unique per constructed registry.
    pub instance_nonce: [u8; 32],
}

impl FreezeOwnerBinding {
    fn registry_identity(&self) -> Result<[u8; 32], SnapshotFreezeError> {
        if self.instance_nonce == [0_u8; 32] {
            return Err(SnapshotFreezeError::ZeroInstanceNonce);
        }
        let mut preimage = Vec::with_capacity(80);
        preimage.extend_from_slice(&self.owner);
        preimage.extend_from_slice(&self.session.to_le_bytes());
        preimage.extend_from_slice(&self.solver_instance_generation.to_le_bytes());
        preimage.extend_from_slice(&self.instance_nonce);
        Ok(*fs_blake3::hash_domain(FREEZE_REGISTRY_IDENTITY_DOMAIN, &preimage).as_bytes())
    }
}

/// Single-use freeze/burn nonce bound to exactly one transaction attempt.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct FreezeNonce([u8; 32]);

impl FreezeNonce {
    /// Exact domain-separated nonce identity.
    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; 32] {
        &self.0
    }
}

/// Explicit resume-semantics inputs required to build the expected context.
///
/// Typed declarations owned by the session/solver layer; the freeze path never
/// infers them from candidate bytes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FreezeResumeInputs {
    /// Algorithm family identity.
    pub algorithm: SnapshotAlgorithmIdV2,
    /// Algorithm implementation/semantic version.
    pub algorithm_version: u64,
    /// Semantic problem identity.
    pub problem: SnapshotProblemIdV2,
    /// RNG streams/counters/stochastic cursor identity at the pause boundary.
    pub rng_counter: SnapshotRngCounterIdV2,
    /// Declared determinism contract.
    pub determinism: SnapshotDeterminismV2,
    /// ISA/numeric/dispatch fingerprint required for replay.
    pub execution_fingerprint: SnapshotExecutionFingerprintIdV2,
    /// Remaining/consumed budget-state identity at the pause boundary.
    pub budget_state: SnapshotBudgetStateIdV2,
    /// Complete run/ledger provenance-context identity.
    pub provenance: SnapshotProvenanceIdV2,
}

/// Labels the session owner declares alongside its pause request.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FreezeBoundaryLabels {
    /// Pause-request identity previously minted by the session layer.
    pub pause_request: SnapshotPauseRequestIdV2,
    /// Session gate generation at which the old run stopped.
    pub gate_generation: u64,
}

/// Terminal phase of the single transaction slot inside a registry.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FreezePhase {
    /// Fresh registry; the transaction may still begin.
    Armed,
    /// Admission closed; state capture has not succeeded yet.
    Admitted,
    /// State captured and encoded once; the permit exists somewhere.
    Frozen,
    /// Burn-before-call recorded; sealing work in flight or just finished.
    Burning,
    /// Sealed bytes exist; the nonce is spent.
    Committed,
    /// Terminal failure. No further transaction can ever succeed.
    Poisoned { reason: FreezePoisonReason },
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FreezePoisonReason {
    DroppedBeforeCapture,
    CaptureRefused,
    DroppedBeforeSeal,
    SealRefused,
    SealPanic,
}

impl FreezePoisonReason {
    const fn describe(self) -> &'static str {
        match self {
            Self::DroppedBeforeCapture => "request dropped before state capture",
            Self::CaptureRefused => "state capture refused",
            Self::DroppedBeforeSeal => "permit dropped before sealing",
            Self::SealRefused => "envelope seal refused",
            Self::SealPanic => "panic unwound through sealing",
        }
    }
}

struct RegistrySlot {
    phase: FreezePhase,
    nonce: Option<FreezeNonce>,
    binding: Option<FreezeOwnerBinding>,
    labels: Option<FreezeBoundaryLabels>,
}

/// Owner-held authority that mints and burns snapshot freeze permits.
///
/// The session owner keeps this value private. Whoever can reach it can
/// authorize freezes for exactly one (owner, session, solver-instance
/// generation); nobody else can obtain a permit because every minting path
/// starts here and each registry performs exactly one transaction.
pub struct SnapshotFreezeRegistry {
    id: [u8; 32],
    binding: FreezeOwnerBinding,
    slot: std::sync::Mutex<RegistrySlot>,
}

impl fmt::Debug for SnapshotFreezeRegistry {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("SnapshotFreezeRegistry")
            .field("binding", &self.binding)
            .finish_non_exhaustive()
    }
}

impl SnapshotFreezeRegistry {
    /// Construct the one-shot registry for `binding`.
    ///
    /// # Errors
    /// [`SnapshotFreezeError::ZeroInstanceNonce`] on a zero instance nonce.
    pub fn new(binding: FreezeOwnerBinding) -> Result<Self, SnapshotFreezeError> {
        let id = binding.registry_identity()?;
        Ok(Self {
            id,
            binding,
            slot: std::sync::Mutex::new(RegistrySlot {
                phase: FreezePhase::Armed,
                nonce: None,
                binding: None,
                labels: None,
            }),
        })
    }

    /// Full-width registry identity.
    #[must_use]
    pub fn id(&self) -> &[u8; 32] {
        &self.id
    }

    /// The binding this registry was constructed for.
    #[must_use]
    pub const fn binding(&self) -> &FreezeOwnerBinding {
        &self.binding
    }

    /// Begin the single freeze transaction and close mutation admission.
    ///
    /// After this call the registry refuses every further transaction forever,
    /// whether or not capture later succeeds. Dropping the returned request
    /// without freezing poisons the registry.
    ///
    /// # Errors
    /// [`SnapshotFreezeError::TransactionAlreadyBegan`] unless freshly armed;
    /// [`SnapshotFreezeError::PoisonedTerminal`] after poisoning.
    pub fn begin_freeze(&self) -> Result<SnapshotFreezeRequest<'_>, SnapshotFreezeError> {
        let mut slot = self
            .slot
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        match slot.phase {
            FreezePhase::Armed => {
                slot.phase = FreezePhase::Admitted;
                Ok(SnapshotFreezeRequest {
                    registry: self,
                    live: true,
                })
            }
            FreezePhase::Poisoned { reason } => {
                Err(SnapshotFreezeError::PoisonedTerminal(reason.describe()))
            }
            _ => Err(SnapshotFreezeError::TransactionAlreadyBegan),
        }
    }

    fn poison(&self, reason: FreezePoisonReason) {
        if let Ok(mut slot) = self.slot.lock() {
            if !matches!(slot.phase, FreezePhase::Committed | FreezePhase::Poisoned { .. }) {
                slot.phase = FreezePhase::Poisoned { reason };
            }
        }
    }

    /// Burn-before-call: only the exact frozen record accepts it; forged or
    /// foreign permits are refused WITHOUT burning so the legitimate holder
    /// can still seal.
    fn burn_for_seal(&self, core: &PermitCore) -> Result<(), SnapshotSealError> {
        let mut slot = self
            .slot
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        match slot.phase {
            FreezePhase::Frozen => {}
            FreezePhase::Burning | FreezePhase::Committed => {
                return Err(SnapshotSealError::AlreadySpent);
            }
            FreezePhase::Poisoned { reason } => {
                return Err(SnapshotSealError::PoisonedTerminal(reason.describe()));
            }
            _ => return Err(SnapshotSealError::WrongPhase { expected: "frozen" }),
        }
        let recorded_ok = match (slot.nonce, slot.binding, slot.labels) {
            (Some(nonce), Some(binding), Some(labels)) => {
                nonce == core.nonce
                    && binding == core.binding
                    && labels == core.labels
                    && self.id == core.registry_id
            }
            _ => false,
        };
        if !recorded_ok {
            return Err(SnapshotSealError::RegistryMismatch);
        }
        slot.phase = FreezePhase::Burning;
        Ok(())
    }

    fn finish_committed(&self) {
        if let Ok(mut slot) = self.slot.lock() {
            if slot.phase == FreezePhase::Burning {
                slot.phase = FreezePhase::Committed;
            }
        }
    }

    fn record_capture(&self, core: &PermitCore) {
        if let Ok(mut slot) = self.slot.lock() {
            debug_assert_eq!(slot.phase, FreezePhase::Admitted);
            slot.phase = FreezePhase::Frozen;
            slot.nonce = Some(core.nonce);
            slot.binding = Some(core.binding);
            slot.labels = Some(core.labels);
        }
    }
}


/// Open transaction after [`SnapshotFreezeRegistry::begin_freeze`].
///
/// Mutation admission is closed. Dropping this value without capturing a state
/// poisons the registry: abandoned admissions cannot be silently reused.
#[derive(Debug)]
#[must_use = "complete the freeze or drop it knowingly; unconsumed requests poison their registry"]
pub struct SnapshotFreezeRequest<'reg> {
    registry: &'reg SnapshotFreezeRegistry,
    live: bool,
}

impl<'reg> SnapshotFreezeRequest<'reg> {
    /// Capture `state`, drain the exact run, and encode the canonical payload
    /// exactly once.
    ///
    /// The drain report comes from [`DrainTracker::finalize`] itself; callers
    /// cannot pair a copied report with these labels. On success the registry
    /// holds the freeze record and the returned permit is the only path to
    /// sealing.
    ///
    /// # Errors
    /// Every refusal poisons the registry terminal (fail closed): drain
    /// refusals as [`SnapshotFreezeError::Drain`], encode refusals as
    /// [`SnapshotFreezeError::Encode`].
    pub fn freeze<S>(
        mut self,
        state: S,
        drain: &DrainTracker<'_>,
        labels: FreezeBoundaryLabels,
        inputs: FreezeResumeInputs,
        limits: SnapshotLimitsV2,
        encode_cancellation: &mut dyn fs_blake3::identity::CancellationProbe,
    ) -> Result<SnapshotFreezePermit<'reg, S>, SnapshotFreezeError>
    where
        S: SolverStateV2,
    {
        self.live = false;
        // Drain first: the report must exist before any state byte is read.
        let report = drain.finalize().map_err(|error| {
            self.registry.poison(FreezePoisonReason::CaptureRefused);
            SnapshotFreezeError::Drain(error)
        })?;
        let nonce = derive_nonce(&self.registry.id, &labels);
        // Commitment covers every transaction identity plus the exact
        // drain-report hash.
        let mut commitment_input = Vec::with_capacity(112);
        commitment_input.extend_from_slice(&self.registry.id);
        commitment_input.extend_from_slice(nonce.as_bytes());
        commitment_input.extend_from_slice(labels.pause_request.as_bytes());
        commitment_input.extend_from_slice(&labels.gate_generation.to_le_bytes());
        commitment_input.extend_from_slice(report.content_hash().as_bytes());
        let commitment =
            *fs_blake3::hash_domain(FREEZE_COMMITMENT_IDENTITY_DOMAIN, &commitment_input)
                .as_bytes();
        let core = PermitCore {
            registry_id: self.registry.id,
            binding: self.registry.binding,
            nonce,
            labels,
            inputs,
            report,
            commitment,
        };
        // Encode exactly once; the payload lives inside the permit from now
        // on, so sealing can never observe a mutated state.
        let payload = match snapshot_v2::encode_state_payload::<S>(
            limits,
            encode_cancellation,
            &state,
        ) {
            Ok(payload) => payload,
            Err(error) => {
                self.registry.poison(FreezePoisonReason::CaptureRefused);
                return Err(SnapshotFreezeError::Encode(error));
            }
        };
        self.registry.record_capture(&core);
        Ok(SnapshotFreezePermit {
            registry: self.registry,
            core,
            payload,
            state: Some(state),
            limits_used: limits,
        })
    }
}

impl Drop for SnapshotFreezeRequest<'_> {
    fn drop(&mut self) {
        if self.live {
            self.registry.poison(FreezePoisonReason::DroppedBeforeCapture);
        }
    }
}

/// Private shared permit data. Fully `Copy` so sealing can consume it without
/// partial-move gymnastics.
#[derive(Clone, Copy)]
struct PermitCore {
    registry_id: [u8; 32],
    binding: FreezeOwnerBinding,
    nonce: FreezeNonce,
    labels: FreezeBoundaryLabels,
    inputs: FreezeResumeInputs,
    report: DrainFinalizeReport,
    commitment: [u8; 32],
}

impl fmt::Debug for PermitCore {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("PermitCore")
            .field("binding", &self.binding)
            .field("nonce", &self.nonce)
            .field("drained_run", &self.report.run().0)
            .finish_non_exhaustive()
    }
}

fn derive_nonce(registry_id: &[u8; 32], labels: &FreezeBoundaryLabels) -> FreezeNonce {
    let mut preimage = Vec::with_capacity(72);
    preimage.extend_from_slice(registry_id);
    preimage.extend_from_slice(labels.pause_request.as_bytes());
    preimage.extend_from_slice(&labels.gate_generation.to_le_bytes());
    FreezeNonce(*fs_blake3::hash_domain(FREEZE_NONCE_IDENTITY_DOMAIN, &preimage).as_bytes())
}

/// Single-use authority to seal exactly the frozen state bytes.
///
/// Minted only by [`SnapshotFreezeRequest::freeze`]. Linear: moving it is
/// allowed; cloning, copying, serializing, or rebuilding it is not. Sealing
/// burns the nonce in the minting registry before any work and poisons the
/// registry on panic, cancellation, or refusal.
#[must_use = "an unsealed frozen state poisons its registry when dropped"]
pub struct SnapshotFreezePermit<'reg, S> {
    registry: &'reg SnapshotFreezeRegistry,
    core: PermitCore,
    payload: Vec<u8>,
    state: Option<S>,
    limits_used: SnapshotLimitsV2,
}

impl<S> fmt::Debug for SnapshotFreezePermit<'_, S> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("SnapshotFreezePermit")
            .field("core", &self.core)
            .field("payload_bytes", &self.payload.len())
            .finish_non_exhaustive()
    }
}

impl<S> SnapshotFreezePermit<'_, S> {
    /// Pre-seal disposition: canonical payload bytes exist and are bound into
    /// the registry record, but no envelope exists yet.
    #[must_use]
    pub fn disposition(&self) -> FreezeDisposition {
        FreezeDisposition::BytesPrepared
    }

    /// Domain-separated commitment over the transaction identities and the
    /// executor drain-report hash.
    #[must_use]
    pub const fn payload_commitment(&self) -> &[u8; 32] {
        &self.core.commitment
    }

    /// The logical run whose workers drained before this capture.
    #[must_use]
    pub const fn drained_run(&self) -> RunId {
        self.core.report.run()
    }
}

impl<S: SolverStateV2> SnapshotFreezePermit<'_, S> {
    /// Seal the frozen payload into its final envelope, consuming the permit.
    ///
    /// Burn-before-call: the registry records the burning phase before any
    /// sealing work, so a crashed or duplicated sealer can never admit a
    /// second attempt under the same identity. Panic, cancellation, and every
    /// seal refusal poison the registry terminal.
    ///
    /// # Errors
    /// [`SnapshotSealError::Envelope`] on structural/refusal failure;
    /// [`SnapshotSealError::RegistryMismatch`] for a permit that does not
    /// match the registry's frozen record; [`SnapshotSealError::WrongPhase`],
    /// [`SnapshotSealError::AlreadySpent`],
    /// [`SnapshotSealError::PoisonedTerminal`] for stale/replayed/duplicate
    /// burns.
    pub fn seal<C>(self, seal_cancellation: C) -> Result<CommittedFreeze<S>, SnapshotSealError>
    where
        C: fs_blake3::identity::CancellationProbe,
    {
        self.registry.burn_for_seal(&self.core)?;
        // From here every exit path terminates the slot explicitly; suppress
        // Drop so it cannot double-poison after ownership moved out.
        let mut this = std::mem::ManuallyDrop::new(self);
        let core = this.core;
        let payload = core::mem::take(&mut this.payload);
        let state = this.state.take();
        let result = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            seal_registered_payload::<S, C>(&core, &payload, this.limits_used, seal_cancellation)
        }));
        match result {
            Ok(Ok(sealed)) => {
                let receipt = build_receipt(&core, &sealed);
                this.registry.finish_committed();
                Ok(CommittedFreeze {
                    state,
                    sealed,
                    receipt,
                })
            }
            Ok(Err(error)) => {
                this.registry.poison(FreezePoisonReason::SealRefused);
                Err(SnapshotSealError::Envelope(error))
            }
            Err(panic_payload) => {
                this.registry.poison(FreezePoisonReason::SealPanic);
                std::panic::resume_unwind(panic_payload);
            }
        }
    }

    /// Test-only duplicate construction used to exercise forged-duplicate and
    /// race refusals. Never compiled outside unit tests; not part of any
    /// public contract.
    #[cfg(test)]
    #[doc(hidden)]
    pub fn test_duplicate_for_refusal_probes(&self) -> Self {
        Self {
            registry: self.registry,
            core: self.core,
            payload: self.payload.clone(),
            state: None,
            limits_used: self.limits_used,
        }
    }

    /// Test-only relabeling used to prove label substitution fails closed.
    #[cfg(test)]
    #[doc(hidden)]
    pub fn test_with_labels(&self, labels: FreezeBoundaryLabels) -> Self {
        let mut core = self.core;
        core.labels = labels;
        core.nonce = derive_nonce(&core.registry_id, &labels);
        Self {
            registry: self.registry,
            core,
            payload: self.payload.clone(),
            state: None,
            limits_used: self.limits_used,
        }
    }

    /// Test-only view of the minting registry identity.
    #[cfg(test)]
    #[doc(hidden)]
    pub fn test_registry_id(&self) -> &[u8; 32] {
        &self.core.registry_id
    }
}

impl<S> Drop for SnapshotFreezePermit<'_, S> {
    fn drop(&mut self) {
        // Reaching Drop un-sealed means the permit was abandoned. Poisoning is
        // idempotent: already-committed or already-poisoned slots stay put.
        self.registry.poison(FreezePoisonReason::DroppedBeforeSeal);
    }
}

/// Outcome of a successful [`SnapshotFreezePermit::seal`].
///
/// The frozen state stays attached so the paused owner retains its
/// authoritative in-memory copy until activation decides otherwise; it is not
/// re-encoded.
pub struct CommittedFreeze<S> {
    state: Option<S>,
    sealed: SealedSnapshotV2,
    receipt: SnapshotFreezeReceipt,
}

impl<S> fmt::Debug for CommittedFreeze<S> {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("CommittedFreeze")
            .field("receipt", &self.receipt)
            .finish_non_exhaustive()
    }
}

impl<S: SolverStateV2> CommittedFreeze<S> {
    /// Exact sealed envelope bytes.
    #[must_use]
    pub fn bytes(&self) -> &[u8] {
        self.sealed.bytes()
    }

    /// Typed committed receipt.
    #[must_use]
    pub const fn receipt(&self) -> &SnapshotFreezeReceipt {
        &self.receipt
    }

    /// Take the frozen state out of the committed bundle.
    pub fn take_state(&mut self) -> Option<S> {
        self.state.take()
    }
}

/// Precise statement of how far a freeze artifact progressed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FreezeDisposition {
    /// Canonical payload encoded and bound into the registry record; no
    /// envelope exists yet.
    BytesPrepared,
    /// Complete sealed envelope bytes exist; publication/activation remain
    /// downstream obligations recorded elsewhere.
    BytesCommitted,
}

/// Typed receipt for a completed freeze transaction.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SnapshotFreezeReceipt {
    registry_id: [u8; 32],
    binding: FreezeOwnerBinding,
    nonce: FreezeNonce,
    labels: FreezeBoundaryLabels,
    drained_run: u64,
    registered_workers: u64,
    drained_workers: u64,
    drain_report_hash: [u8; 32],
    payload_commitment: [u8; 32],
    sealed_content_id: [u8; 32],
    disposition: FreezeDisposition,
    identity: [u8; 32],
}

impl SnapshotFreezeReceipt {
    /// Domain-separated receipt identity over every semantic field.
    #[must_use]
    pub const fn identity(&self) -> &[u8; 32] {
        &self.identity
    }

    /// Registry identity that minted and burned the permit.
    #[must_use]
    pub const fn registry_id(&self) -> &[u8; 32] {
        &self.registry_id
    }

    /// Owner binding of the freeze.
    #[must_use]
    pub const fn binding(&self) -> &FreezeOwnerBinding {
        &self.binding
    }

    /// Burned single-use nonce.
    #[must_use]
    pub const fn nonce(&self) -> &FreezeNonce {
        &self.nonce
    }

    /// Pause labels bound at capture time.
    #[must_use]
    pub const fn labels(&self) -> &FreezeBoundaryLabels {
        &self.labels
    }

    /// Logical run whose workers drained.
    #[must_use]
    pub const fn drained_run(&self) -> u64 {
        self.drained_run
    }

    /// Worker counts witnessed by the executor drain tracker.
    #[must_use]
    pub const fn worker_counts(&self) -> (u64, u64) {
        (self.registered_workers, self.drained_workers)
    }

    /// Executor drain-report content identity.
    #[must_use]
    pub const fn drain_report_hash(&self) -> &[u8; 32] {
        &self.drain_report_hash
    }

    /// Commitment covering transaction identities and the drain report hash.
    #[must_use]
    pub const fn payload_commitment(&self) -> &[u8; 32] {
        &self.payload_commitment
    }

    /// Content identity of the sealed envelope.
    #[must_use]
    pub const fn sealed_content_id(&self) -> &[u8; 32] {
        &self.sealed_content_id
    }

    /// Precise progress statement of this receipt.
    #[must_use]
    pub const fn disposition(&self) -> FreezeDisposition {
        self.disposition
    }
}

#[cfg(test)]
fn test_limits() -> SnapshotLimitsV2 {
    snapshot_v2::SnapshotLimitsV2::new(
        1 << 20,
        64,
        fs_blake3::identity::CanonicalLimits::new(16_384, 4_096, 32, 32, 64),
        4_096,
        1 << 20,
        64,
    )
}

/// Refusal while driving the transaction up to permit minting.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SnapshotFreezeError {
    /// The caller supplied the all-zero instance nonce.
    ZeroInstanceNonce,
    /// This registry already began (and closed) its single transaction.
    TransactionAlreadyBegan,
    /// The registry terminal was poisoned by an earlier failure.
    PoisonedTerminal(&'static str),
    /// The drain tracker refused to attest a fully drained run.
    Drain(DrainFinalizeError),
    /// Payload encoding refused.
    Encode(SnapshotV2Error),
}

impl From<DrainFinalizeError> for SnapshotFreezeError {
    fn from(value: DrainFinalizeError) -> Self {
        Self::Drain(value)
    }
}

impl fmt::Display for SnapshotFreezeError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroInstanceNonce => write!(formatter, "freeze instance nonce must be nonzero"),
            Self::TransactionAlreadyBegan => {
                write!(formatter, "freeze registry already ran its single transaction")
            }
            Self::PoisonedTerminal(reason) => {
                write!(formatter, "freeze terminal poisoned: {reason}")
            }
            Self::Drain(error) => write!(formatter, "drain finalize refused: {error}"),
            Self::Encode(error) => write!(formatter, "payload encode refused: {error}"),
        }
    }
}

impl core::error::Error for SnapshotFreezeError {}

/// Refusal or misuse detected while sealing a permit.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SnapshotSealError {
    /// Structural or resource refusal from the v2 envelope machinery.
    Envelope(SnapshotV2Error),
    /// The permit did not originate from this registry's frozen record.
    RegistryMismatch,
    /// The slot was not in the phase required at burn time.
    WrongPhase {
        /// Phase name expected at burn time.
        expected: &'static str,
    },
    /// The nonce was already burned: replay of a consumed permit.
    AlreadySpent,
    /// The registry terminal was poisoned by an earlier failure.
    PoisonedTerminal(&'static str),
    /// The minting registry mutex was poisoned by an unwind race.
    RegistryGone,
}

impl fmt::Display for SnapshotSealError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Envelope(error) => write!(formatter, "envelope seal refused: {error}"),
            Self::RegistryMismatch => {
                write!(formatter, "permit does not match the registry's frozen record")
            }
            Self::WrongPhase { expected } => {
                write!(formatter, "slot phase was not `{expected}` at burn time")
            }
            Self::AlreadySpent => write!(formatter, "freeze nonce already burned"),
            Self::PoisonedTerminal(reason) => {
                write!(formatter, "freeze terminal poisoned: {reason}")
            }
            Self::RegistryGone => write!(formatter, "minting registry lock failed"),
        }
    }
}

impl core::error::Error for SnapshotSealError {}

struct ProbeRef<'a>(&'a mut dyn fs_blake3::identity::CancellationProbe);

impl fs_blake3::identity::CancellationProbe for ProbeRef<'_> {
    fn is_cancelled(&mut self) -> bool {
        self.0.is_cancelled()
    }
}

fn seal_registered_payload<S: SolverStateV2, C: fs_blake3::identity::CancellationProbe>(
    core: &PermitCore,
    payload: &[u8],
    limits: SnapshotLimitsV2,
    cancellation: C,
) -> Result<SealedSnapshotV2, SnapshotV2Error> {
    let boundary = PausedSnapshotBoundaryV2::from_drain_report(
        core.report,
        core.labels.pause_request,
        core.labels.gate_generation,
    );
    let expected = ExpectedResumeContextV2::for_paused_state::<S>(
        core.inputs.algorithm,
        core.inputs.algorithm_version,
        core.inputs.problem,
        core.inputs.rng_counter,
        core.inputs.determinism,
        core.inputs.execution_fingerprint,
        core.inputs.budget_state,
        core.inputs.provenance,
        boundary,
    );
    let mut bytes = Vec::with_capacity(payload.len());
    bytes.extend_from_slice(payload);
    let mut cancellation = cancellation;
    let probe = ProbeRef(&mut cancellation);
    snapshot_v2::seal_encoded_payload(bytes, &expected, limits, probe)
}

fn build_receipt(core: &PermitCore, sealed: &SealedSnapshotV2) -> SnapshotFreezeReceipt {
    let receipt = SnapshotFreezeReceipt {
        registry_id: core.registry_id,
        binding: core.binding,
        nonce: core.nonce,
        labels: core.labels,
        drained_run: core.report.run().0,
        registered_workers: core.report.registered_workers(),
        drained_workers: core.report.drained_workers(),
        drain_report_hash: *core.report.content_hash().as_bytes(),
        payload_commitment: core.commitment,
        sealed_content_id: *sealed.content_id().as_bytes(),
        disposition: FreezeDisposition::BytesCommitted,
        identity: [0_u8; 32],
    };
    let identity_input = receipt_preimage(&receipt);
    let identity =
        *fs_blake3::hash_domain(FREEZE_RECEIPT_IDENTITY_DOMAIN, &identity_input).as_bytes();
    SnapshotFreezeReceipt {
        identity,
        ..receipt
    }
}

fn receipt_preimage(receipt: &SnapshotFreezeReceipt) -> Vec<u8> {
    let mut out = Vec::with_capacity(288);
    out.extend_from_slice(&receipt.registry_id);
    out.extend_from_slice(&receipt.binding.owner);
    out.extend_from_slice(&receipt.binding.session.to_le_bytes());
    out.extend_from_slice(&receipt.binding.solver_instance_generation.to_le_bytes());
    out.extend_from_slice(&receipt.binding.instance_nonce);
    out.extend_from_slice(receipt.nonce.as_bytes());
    out.extend_from_slice(receipt.labels.pause_request.as_bytes());
    out.extend_from_slice(&receipt.labels.gate_generation.to_le_bytes());
    out.extend_from_slice(&receipt.drained_run.to_le_bytes());
    out.extend_from_slice(&receipt.registered_workers.to_le_bytes());
    out.extend_from_slice(&receipt.drained_workers.to_le_bytes());
    out.extend_from_slice(&receipt.drain_report_hash);
    out.extend_from_slice(&receipt.payload_commitment);
    out.extend_from_slice(&receipt.sealed_content_id);
    out.push(match receipt.disposition {
        FreezeDisposition::BytesPrepared => 1_u8,
        FreezeDisposition::BytesCommitted => 2_u8,
    });
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cx::{CancelGate, ExecMode};
    use crate::solver::snapshot_v2::{
        SnapshotStateCodecIdV2, SnapshotStateSchemaIdV2, SnapshotStateTypeIdV2,
    };
    use std::sync::Arc;
    use std::sync::Barrier;

    const fn id32(tag: u8) -> [u8; 32] {
        [tag; 32]
    }

    struct TestState {
        marker: u64,
    }

    impl SolverStateV2 for TestState {
        const STATE_TYPE_ID_V2: SnapshotStateTypeIdV2 =
            SnapshotStateTypeIdV2::from_bytes(id32(0xA1));
        const STATE_SCHEMA_ID_V2: SnapshotStateSchemaIdV2 =
            SnapshotStateSchemaIdV2::from_bytes(id32(0xA2));
        const STATE_CODEC_ID_V2: SnapshotStateCodecIdV2 =
            SnapshotStateCodecIdV2::from_bytes(id32(0xA3));
        const STATE_CODEC_VERSION_V2: u32 = 1;

        fn encode_v2(
            &self,
            encoder: &mut snapshot_v2::SnapshotEncoderV2<'_>,
        ) -> Result<(), SnapshotV2Error> {
            encoder.put_u64(self.marker)
        }

        fn decode_v2(
            decoder: &mut snapshot_v2::SnapshotDecoderV2<'_, '_>,
        ) -> Result<Self, SnapshotV2Error> {
            Ok(Self {
                marker: decoder.get_u64()?,
            })
        }
    }

    /// Test probe: optional cancellation, optional seeded panic on first poll.
    #[derive(Clone)]
    struct TestProbe {
        cancelled: bool,
        panic_on_poll: bool,
    }

    impl TestProbe {
        fn ok() -> Self {
            Self {
                cancelled: false,
                panic_on_poll: false,
            }
        }
    }

    impl fs_blake3::identity::CancellationProbe for TestProbe {
        fn is_cancelled(&mut self) -> bool {
            if self.panic_on_poll {
                panic!("seeded seal unwind");
            }
            self.cancelled
        }
    }

    fn owner_binding(nonce_tag: u8) -> FreezeOwnerBinding {
        FreezeOwnerBinding {
            owner: id32(0x10),
            session: 7,
            solver_instance_generation: 3,
            instance_nonce: id32(nonce_tag),
        }
    }

    fn resume_inputs() -> FreezeResumeInputs {
        FreezeResumeInputs {
            algorithm: snapshot_v2::SnapshotAlgorithmIdV2::from_bytes(id32(0xB1)),
            algorithm_version: 11,
            problem: snapshot_v2::SnapshotProblemIdV2::from_bytes(id32(0xB2)),
            rng_counter: snapshot_v2::SnapshotRngCounterIdV2::from_bytes(id32(0xB3)),
            determinism: snapshot_v2::SnapshotDeterminismV2::Deterministic,
            execution_fingerprint: snapshot_v2::SnapshotExecutionFingerprintIdV2::from_bytes(
                id32(0xB4),
            ),
            budget_state: snapshot_v2::SnapshotBudgetStateIdV2::from_bytes(id32(0xB5)),
            provenance: snapshot_v2::SnapshotProvenanceIdV2::from_bytes(id32(0xB6)),
        }
    }

    fn boundary_labels(tag: u8) -> FreezeBoundaryLabels {
        FreezeBoundaryLabels {
            pause_request: snapshot_v2::SnapshotPauseRequestIdV2::from_bytes(id32(tag)),
            gate_generation: 42,
        }
    }

    fn test_limits() -> SnapshotLimitsV2 {
        super::test_limits()
    }

    /// Gate requested + one worker registered and drained.
    fn drained_tracker<'a>(gate: &'a CancelGate) -> DrainTracker<'a> {
        gate.request();
        let tracker = DrainTracker::new(RunId(77), gate);
        drop(tracker.register_worker().expect("register worker"));
        tracker
    }

    fn mint_permit<'a>(
        registry: &'a SnapshotFreezeRegistry,
        gate: &'a CancelGate,
    ) -> Result<SnapshotFreezePermit<'a, TestState>, SnapshotFreezeError> {
        let request = registry.begin_freeze()?;
        let tracker = drained_tracker(gate);
        request.freeze(
            TestState { marker: 0x5EED },
            &tracker,
            boundary_labels(0xC1),
            resume_inputs(),
            test_limits(),
            &mut TestProbe::ok(),
        )
    }

    #[test]
    fn zero_instance_nonce_is_refused() {
        let mut binding = owner_binding(0x01);
        binding.instance_nonce = [0_u8; 32];
        assert_eq!(
            SnapshotFreezeRegistry::new(binding).unwrap_err(),
            SnapshotFreezeError::ZeroInstanceNonce
        );
    }

    #[test]
    fn committed_receipt_binds_every_identity() {
        let registry = SnapshotFreezeRegistry::new(owner_binding(0x01)).unwrap();
        let gate = CancelGate::new();
        let permit = mint_permit(&registry, &gate).expect("freeze");
        assert_eq!(permit.disposition(), FreezeDisposition::BytesPrepared);
        assert_eq!(permit.drained_run(), RunId(77));

        let committed = permit.seal(TestProbe::ok()).expect("seal");
        let receipt = committed.receipt();
        assert_eq!(receipt.disposition(), FreezeDisposition::BytesCommitted);
        assert_eq!(receipt.registry_id(), registry.id());
        assert_eq!(*receipt.binding(), owner_binding(0x01));
        assert_eq!(receipt.labels(), &boundary_labels(0xC1));
        assert_eq!(receipt.drained_run(), 77);
        assert_eq!(receipt.worker_counts(), (1, 1));
        assert_ne!(receipt.identity(), &[0_u8; 32]);
        assert_ne!(receipt.sealed_content_id(), &[0_u8; 32]);
        // Deterministic recompute: same inputs produce the same receipt id.
        let identity_input = receipt_preimage(receipt);
        assert_eq!(
            *fs_blake3::hash_domain(FREEZE_RECEIPT_IDENTITY_DOMAIN, &identity_input).as_bytes(),
            *receipt.identity()
        );
    }

    #[test]
    fn second_begin_is_refused_after_commit() {
        let registry = SnapshotFreezeRegistry::new(owner_binding(0x02)).unwrap();
        let gate = CancelGate::new();
        let permit = mint_permit(&registry, &gate).unwrap();
        let _committed = permit.seal(TestProbe::ok()).unwrap();
        assert_eq!(
            registry.begin_freeze().unwrap_err(),
            SnapshotFreezeError::TransactionAlreadyBegan
        );
    }

    #[test]
    fn unconsumed_request_drop_poisons_terminal() {
        let registry = SnapshotFreezeRegistry::new(owner_binding(0x03)).unwrap();
        let first = registry.begin_freeze().unwrap();
        assert_eq!(
            registry.begin_freeze().unwrap_err(),
            SnapshotFreezeError::TransactionAlreadyBegan
        );
        drop(first);
        assert!(matches!(
            registry.begin_freeze().unwrap_err(),
            SnapshotFreezeError::PoisonedTerminal(_)
        ));
    }

    #[test]
    fn drain_refusal_poisons_terminal() {
        let registry = SnapshotFreezeRegistry::new(owner_binding(0x04)).unwrap();
        let gate = CancelGate::new();
        // Gate NOT requested: finalize must refuse and poison.
        let tracker = DrainTracker::new(RunId(5), &gate);
        drop(tracker.register_worker().unwrap());
        let request = registry.begin_freeze().unwrap();
        let error = request
            .freeze(
                TestState { marker: 1 },
                &tracker,
                boundary_labels(0xC2),
                resume_inputs(),
                test_limits(),
                &mut TestProbe::ok(),
            )
            .unwrap_err();
        assert!(matches!(error, SnapshotFreezeError::Drain(_)));
        assert!(matches!(
            registry.begin_freeze().unwrap_err(),
            SnapshotFreezeError::PoisonedTerminal(_)
        ));
    }

    #[test]
    fn cancelled_encode_poisons_terminal() {
        let registry = SnapshotFreezeRegistry::new(owner_binding(0x05)).unwrap();
        let gate = CancelGate::new();
        let request = registry.begin_freeze().unwrap();
        let tracker = drained_tracker(&gate);
        let error = request
            .freeze(
                TestState { marker: 2 },
                &tracker,
                boundary_labels(0xC3),
                resume_inputs(),
                test_limits(),
                &mut TestProbe {
                    cancelled: true,
                    panic_on_poll: false,
                },
            )
            .unwrap_err();
        assert!(matches!(error, SnapshotFreezeError::Encode(_)));
        assert!(matches!(
            registry.begin_freeze().unwrap_err(),
            SnapshotFreezeError::PoisonedTerminal(_)
        ));
    }

    #[test]
    fn seal_panic_poisons_terminal_and_propagates() {
        let registry = SnapshotFreezeRegistry::new(owner_binding(0x06)).unwrap();
        let gate = CancelGate::new();
        let permit = mint_permit(&registry, &gate).unwrap();
        let duplicate = permit.test_duplicate_for_refusal_probes();
        let unwound = std::panic::catch_unwind(std::panic::AssertUnwindSafe(move || {
            let _ = permit.seal(TestProbe {
                cancelled: false,
                panic_on_poll: true,
            });
        }));
        assert!(unwound.is_err(), "panic must propagate after poisoning");
        // Duplicate sealing afterwards hits the poisoned terminal.
        assert!(matches!(
            duplicate.seal(TestProbe::ok()).unwrap_err(),
            SnapshotSealError::PoisonedTerminal(_)
        ));
    }

    #[test]
    fn tampered_labels_fail_closed_without_burning() {
        let registry = SnapshotFreezeRegistry::new(owner_binding(0x07)).unwrap();
        let gate = CancelGate::new();
        let permit = mint_permit(&registry, &gate).unwrap();
        let tampered = permit.test_with_labels(boundary_labels(0xFF));
        assert!(matches!(
            tampered.seal(TestProbe::ok()).unwrap_err(),
            SnapshotSealError::RegistryMismatch
        ));
        // The legitimate permit still seals: nothing was burned.
        let committed = permit.seal(TestProbe::ok()).expect("legitimate seal");
        assert_eq!(
            committed.receipt().disposition(),
            FreezeDisposition::BytesCommitted
        );
    }

    #[test]
    fn replay_after_commit_refuses_as_already_spent() {
        let registry = SnapshotFreezeRegistry::new(owner_binding(0x08)).unwrap();
        let gate = CancelGate::new();
        let permit = mint_permit(&registry, &gate).unwrap();
        let replay = permit.test_duplicate_for_refusal_probes();
        let _committed = permit.seal(TestProbe::ok()).unwrap();
        assert_eq!(
            replay.seal(TestProbe::ok()).unwrap_err(),
            SnapshotSealError::AlreadySpent
        );
    }

    fn duplicate_seals_after_commit_refuse_as_spent() {
        let registry = SnapshotFreezeRegistry::new(owner_binding(0x09)).unwrap();
        let gate = CancelGate::new();
        let permit = mint_permit(&registry, &gate).unwrap();
        let first_duplicate = permit.test_duplicate_for_refusal_probes();
        let second_duplicate = permit.test_duplicate_for_refusal_probes();
        let _committed = permit.seal(TestProbe::ok()).unwrap();
        assert_eq!(
            first_duplicate.seal(TestProbe::ok()).unwrap_err(),
            SnapshotSealError::AlreadySpent
        );
        assert_eq!(
            second_duplicate.seal(TestProbe::ok()).unwrap_err(),
            SnapshotSealError::AlreadySpent
        );
    }

    #[test]
    fn concurrent_begins_admit_exactly_one_transaction() {
        let registry = Arc::new(SnapshotFreezeRegistry::new(owner_binding(0x0A)).unwrap());
        let barrier = Arc::new(Barrier::new(8));
        let handles: Vec<_> = (0..8)
            .map(|_| {
                std::thread::spawn(move || {
                    barrier.wait();
                    let admitted = registry.begin_freeze().is_ok();
                    // Dropping a winning request here poisons the terminal;
                    // admission counting stays exact either way.
                    admitted
                })
            })
            .collect();
        let mut winners = 0;
        for handle in handles {
            if handle.join().unwrap().is_ok() {
                winners += 1;
            }
        }
        assert_eq!(winners, 1);
    }

    #[test]
    fn exec_mode_remains_importable_for_callers() {
        // Guard against accidental removal of re-exports used by sessions.
        let _mode = ExecMode::Deterministic;
    }
}
