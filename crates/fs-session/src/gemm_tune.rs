//! The production GEMM autotune loop (bead yqug): measure → cache →
//! model → dispatch, closed end-to-end.
//!
//! [`gemm_f64_session`] is the production consumer the tuner was built
//! for: it resolves an MC/NC [`GemmBlockPlan`] for the caller's shape
//! class (pins beat cached rows beat the documented cold-start default),
//! runs a BOUNDED candidate sweep when the machine is cold, records the
//! ranked wall-time evidence as a tune row, applies the caller's explicit
//! cache policy, and dispatches
//! `fs_la::gemm_f64_parallel_with_pool` with the selected plan.
//!
//! Honesty boundaries, in the fs-exec tuner's division:
//! - The KERNEL KEY embeds fs-la's `GEMM_BIT_SEMANTICS_VERSION`, so rows
//!   measured under a different accumulation contract can never match a
//!   lookup (semantic filtering by construction). Rows are additionally
//!   bound to the exact probe dims, requested/normalized thread budget,
//!   resolved ISA tier, placement policy, and implementation version, then
//!   the exact compiler/profile/codegen build fingerprint, then machine-
//!   fingerprint-keyed. The ledger read path refuses stale,
//!   differently scoped, non-canonical, and params/body-disagreeing rows.
//! - MC/NC are BIT-NEUTRAL by fs-la's determinism contract, and the
//!   sweep ENFORCES that: every repeat of every effective candidate is
//!   compared word-for-word with the first output, else the loop fails
//!   closed with [`GemmTuneError::BitDrift`] and records nothing. KC is part
//!   of the bit contract and is NOT in this loop. The resolved SIMD tier is
//!   bit-neutral but remains performance identity.
//! - The "cost model" is declared and minimal: argmin of the per-
//!   candidate MINIMUM wall time, ties to the earlier candidate in
//!   lattice order — a recorded selection rule, never a statistical
//!   confidence claim.
//!
//! Determinism class: dispatch results are bit-identical to serial
//! `gemm_f64` for every plan the loop can select (enforced by the sweep
//! and gated in tests); WHICH plan wins is wall-clock-dependent by
//! nature and travels as evidence + a pinnable decision, never inside
//! numeric results.

use fs_exec::{
    CancelGate, GEMM_KERNEL_PREFIX, GemmBlockPlan, GemmExecutionIdentity, GemmTuneKey,
    PreparedGemmDecision, PreparedGemmRow, TilePool, TuneError, TuneEvidence, TuneObservation,
    TuneSource, Tuner,
};
use fs_ledger::Ledger;

/// The bounded sweep lattice: up to 4 × 2 candidates, lattice order
/// (mc-major ascending). Candidates that clamp to an identical effective
/// `(mc, nc)` pair are deduplicated before measurement. Chosen around the
/// measured xlvx s5 landscape: thin bands won both reference machines; the
/// extremes document the neighborhood.
const SWEEP_MC: [usize; 4] = [16, 32, 64, 128];
const SWEEP_NC_CAP: [usize; 2] = [512, 2048];

/// Probe M/K dims are capped so a cold-start sweep stays bounded (seconds,
/// not minutes) even when the caller's problem is huge. N has a separate
/// cap: it must extend beyond the smaller NC candidate or that axis is never
/// measured at all.
const PROBE_MK_DIM_CAP: usize = 512;
const PROBE_N_DIM_CAP: usize = 2048;

/// Wall-time samples per candidate (min-of ranking, all survive in the
/// evidence row).
const SWEEP_SAMPLES: usize = 3;

/// Durable identity of the autotune producer algorithm: candidate lattice,
/// probe dimensions/sample policy, ranking, and plan-to-dispatch mapping. Any
/// semantic change to those choices must bump this value before old rows may be
/// considered compatible.
pub const GEMM_TUNER_SCHEMA_VERSION: u32 = 1;

/// Logical stream seed for the compatibility session pool. GEMM itself is
/// non-stochastic; caller-owned pools retain their study seed for Cx identity.
const SESSION_GEMM_POOL_SEED: u64 = 0x4653_2D53_4553_534E;

const GEMM_SWEEP_RUN_DOMAIN: &str = "org.frankensim.fs-session.gemm-sweep-run.v1";

/// Globally unique BLAKE3 derive-key context for canonical GEMM tune-row
/// receipts.
pub const GEMM_TUNE_ROW_RECEIPT_DOMAIN: &str = "org.frankensim.fs-session.gemm-tune-row-receipt.v2";

fn push_json_string(out: &mut String, value: &str) {
    use core::fmt::Write as _;

    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c.is_control() => {
                let _ = write!(out, "\\u{:04x}", u32::from(c));
            }
            c => out.push(c),
        }
    }
    out.push('"');
}

/// A structured autotune-loop failure. Every variant fails closed: sweep
/// failures record no row and nothing dispatches under unvalidated blocking.
/// A cancellation during the final dispatch may retain the already validated
/// measured row, but records no successful decision and does not commit `C`.
#[derive(Debug)]
pub enum GemmTuneError {
    /// The cancel gate was requested. Compute may have completed in private
    /// staging, but the caller's output was not committed.
    Cancelled {
        /// The caller-visible envelope in force.
        limit_bytes: u64,
        /// Largest session-owned logical reservation concurrency reached.
        peak_used_bytes: u128,
        /// Drained numerical-run report when cancellation was returned by
        /// fs-la. `None` means the gate was observed between dispatch calls;
        /// earlier completed probes may still contribute to the peak.
        report: Option<Box<fs_la::GemmRunReport>>,
    },
    /// Tuner-side refusal (invalid pin, evidence, or adoption).
    Tune(TuneError),
    /// Ledger cache I/O failed (the loop does not guess around storage).
    Ledger(String),
    /// Two sweep candidates produced different output bits: the
    /// bit-neutrality contract is broken and NO plan may be selected.
    BitDrift {
        /// Canonical params of the candidate that diverged.
        candidate: String,
        /// One-based repeat whose exact output bits diverged.
        repeat: usize,
    },
    /// The GEMM path refused at its memory boundary (wf9.15): the plan
    /// exceeded the caller envelope or an allocator declined a reservation.
    /// Output is not committed; the retained report may contain drained
    /// private progress from panels that completed before refusal.
    MemoryRefused {
        /// Which reservation was refused.
        what: &'static str,
        /// Bytes the refused reservation asked for.
        requested_bytes: u128,
        /// The envelope in force.
        limit_bytes: u64,
        /// Largest session-owned logical reservation concurrency reached.
        peak_used_bytes: u128,
        /// Drained numerical-run report when refusal occurred inside fs-la.
        report: Option<Box<fs_la::GemmRunReport>>,
    },
    /// Checked arithmetic could not represent the session or fs-la memory plan.
    MemoryPlanOverflow {
        /// Component whose arithmetic overflowed.
        what: &'static str,
        /// The envelope in force.
        limit_bytes: u64,
    },
    /// The TilePool contained a tile panic or executor invariant refusal.
    Executor {
        /// Structured fs-exec outcome with logical tile provenance.
        error: fs_exec::RunError,
        /// The caller-visible envelope in force.
        limit_bytes: u64,
        /// Largest session-owned logical reservation concurrency reached.
        peak_used_bytes: u128,
        /// Drained numerical-run report, including memory and tile progress.
        report: Box<fs_la::GemmRunReport>,
    },
}

impl core::fmt::Display for GemmTuneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Cancelled {
                limit_bytes,
                peak_used_bytes,
                report,
            } => write!(
                f,
                "gemm work cancelled after {}/{} compute tiles under a {limit_bytes}-byte envelope \
                 after reaching {peak_used_bytes} logical bytes; output not committed",
                report.as_ref().map_or(0, |run| run.completed_tiles),
                report.as_ref().map_or(0, |run| run.total_tiles),
            ),
            Self::Tune(e) => write!(f, "gemm autotune: {e}"),
            Self::Ledger(detail) => write!(f, "gemm autotune ledger cache: {detail}"),
            Self::MemoryRefused {
                what,
                requested_bytes,
                limit_bytes,
                peak_used_bytes,
                report,
            } => write!(
                f,
                "gemm autotune memory refused at {what}: {requested_bytes} bytes against a \
                 {limit_bytes}-byte envelope after reaching {peak_used_bytes} logical bytes and \
                 {}/{} compute tiles; output not committed",
                report.as_ref().map_or(0, |run| run.completed_tiles),
                report.as_ref().map_or(0, |run| run.total_tiles),
            ),
            Self::MemoryPlanOverflow { what, limit_bytes } => write!(
                f,
                "gemm autotune memory-plan arithmetic overflowed at {what} under the \
                 {limit_bytes}-byte envelope; output not touched"
            ),
            Self::BitDrift { candidate, repeat } => write!(
                f,
                "gemm autotune: candidate {candidate} repeat {repeat} broke the MC/NC bit-neutrality contract"
            ),
            Self::Executor {
                error,
                limit_bytes,
                peak_used_bytes,
                report,
            } => write!(
                f,
                "gemm executor failed after {}/{} compute tiles under a {limit_bytes}-byte \
                 envelope after reaching {peak_used_bytes} logical bytes: {error}",
                report.completed_tiles, report.total_tiles,
            ),
        }
    }
}

impl core::error::Error for GemmTuneError {}

impl From<TuneError> for GemmTuneError {
    fn from(e: TuneError) -> Self {
        Self::Tune(e)
    }
}

impl From<fs_la::GemmCancelled> for GemmTuneError {
    fn from(cancelled: fs_la::GemmCancelled) -> Self {
        let limit_bytes = cancelled.report.memory.limit_bytes;
        let peak_used_bytes = cancelled.report.memory.peak_used_bytes;
        Self::Cancelled {
            limit_bytes,
            peak_used_bytes,
            report: Some(Box::new(cancelled.report)),
        }
    }
}

impl From<fs_la::GemmRunError> for GemmTuneError {
    fn from(error: fs_la::GemmRunError) -> Self {
        match error {
            fs_la::GemmRunError::Cancelled(cancelled) => Self::from(cancelled),
            fs_la::GemmRunError::Executor { error, report } => {
                let limit_bytes = report.memory.limit_bytes;
                let peak_used_bytes = report.memory.peak_used_bytes;
                Self::Executor {
                    error,
                    limit_bytes,
                    peak_used_bytes,
                    report: Box::new(report),
                }
            }
            fs_la::GemmRunError::MemoryRefused {
                what,
                requested_bytes,
                limit_bytes,
                report,
            } => Self::MemoryRefused {
                what,
                requested_bytes,
                limit_bytes,
                peak_used_bytes: report.memory.peak_used_bytes,
                report: Some(Box::new(report)),
            },
            fs_la::GemmRunError::MemoryPlanOverflow { what, limit_bytes } => {
                Self::MemoryPlanOverflow { what, limit_bytes }
            }
        }
    }
}

fn gemm_error_with_session_memory(
    error: fs_la::GemmRunError,
    envelope: fs_la::GemmMemoryEnvelope,
    session_bytes: u128,
) -> GemmTuneError {
    match GemmTuneError::from(error) {
        GemmTuneError::MemoryRefused {
            what,
            requested_bytes,
            peak_used_bytes,
            report,
            ..
        } => match session_bytes.checked_add(peak_used_bytes) {
            Some(peak_used_bytes) => GemmTuneError::MemoryRefused {
                what,
                requested_bytes,
                limit_bytes: envelope.limit_bytes,
                peak_used_bytes,
                report,
            },
            None => GemmTuneError::MemoryPlanOverflow {
                what: "session-plus-gemm-peak",
                limit_bytes: envelope.limit_bytes,
            },
        },
        GemmTuneError::Cancelled {
            peak_used_bytes,
            report,
            ..
        } => match session_bytes.checked_add(peak_used_bytes) {
            Some(peak_used_bytes) => GemmTuneError::Cancelled {
                limit_bytes: envelope.limit_bytes,
                peak_used_bytes,
                report,
            },
            None => GemmTuneError::MemoryPlanOverflow {
                what: "session-plus-gemm-peak",
                limit_bytes: envelope.limit_bytes,
            },
        },
        GemmTuneError::Executor {
            error,
            peak_used_bytes,
            report,
            ..
        } => match session_bytes.checked_add(peak_used_bytes) {
            Some(peak_used_bytes) => GemmTuneError::Executor {
                error,
                limit_bytes: envelope.limit_bytes,
                peak_used_bytes,
                report,
            },
            None => GemmTuneError::MemoryPlanOverflow {
                what: "session-plus-gemm-peak",
                limit_bytes: envelope.limit_bytes,
            },
        },
        other => other,
    }
}

fn cancelled_before_compute(envelope: fs_la::GemmMemoryEnvelope) -> GemmTuneError {
    GemmTuneError::Cancelled {
        limit_bytes: envelope.limit_bytes,
        peak_used_bytes: 0,
        report: None,
    }
}

fn cancelled_with_live_probe_memory(
    envelope: fs_la::GemmMemoryEnvelope,
    session_bytes: u128,
    numerical_peak: u128,
) -> GemmTuneError {
    let Some(peak_used_bytes) = session_bytes.checked_add(numerical_peak) else {
        return GemmTuneError::MemoryPlanOverflow {
            what: "session-plus-gemm-peak",
            limit_bytes: envelope.limit_bytes,
        };
    };
    GemmTuneError::Cancelled {
        limit_bytes: envelope.limit_bytes,
        peak_used_bytes,
        report: None,
    }
}

/// The receipt for one autotuned dispatch: what ran, under which plan,
/// and where the plan came from. A study records this; replay pins it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GemmDispatch {
    /// Exact scoped kernel key: numerical version plus complete execution
    /// identity.
    pub kernel: String,
    /// Shape class the plan was resolved for.
    pub shape_class: String,
    /// The MC/NC plan that dispatched.
    pub plan: GemmBlockPlan,
    /// Plan provenance (pinned / tuned / cold-start).
    pub source: TuneSource,
    /// True when this call ran the measurement sweep (cold cache).
    pub swept: bool,
    /// Sealed newly measured row. Read-only cache users can retain this
    /// process-locally and publish it only after their enclosing run passes
    /// admission. `None` when no sweep ran.
    pub new_tune_row: Option<ValidatedGemmTuneRow>,
    /// Sealed row adopted or measured during this call. Callers that need a
    /// citable execution receipt retain this identity across later warm-cache
    /// dispatches. `None` when this call reused an already local row or bypassed
    /// tuning.
    pub validated_tune_row: Option<ValidatedGemmTuneRow>,
    /// Final production execution receipt. Its `pool_runs` prove the selected
    /// plan traversed the caller's TilePool rather than a detached thread path.
    pub run: fs_la::GemmRunReport,
}

impl GemmDispatch {
    /// Deterministic execution facts suitable for replay and evidence binding.
    /// Scheduling measurements (steals, worker distribution, and cancellation
    /// latency) deliberately remain outside this identity.
    #[must_use]
    pub fn execution_receipt(&self) -> GemmExecutionReceipt {
        GemmExecutionReceipt::from_report(&self.run)
    }
}

/// Stable facts for one TilePool traversal of an NC/KC panel.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GemmPanelReceipt {
    /// Stable TileKernel name.
    pub kernel: String,
    /// Deterministic/fast execution mode recorded by the pool.
    pub mode: String,
    /// Deterministic NC/KC panel ordinal used as the fs-exec declared run.
    pub declared_run: u64,
    /// Logical M-band tiles completed.
    pub completed: u64,
    /// Logical M-band tiles planned.
    pub total: u64,
}

/// Deterministic production-path receipt for a GEMM dispatch.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GemmExecutionReceipt {
    /// Caller-ledgered identity of the final production dispatch.
    pub declared_run: u64,
    /// Completed bounded GEMM microtiles.
    pub completed_tiles: usize,
    /// Total bounded GEMM microtiles.
    pub total_tiles: usize,
    /// Declared logical-memory plan. Schedule-observed peak and refusal bytes
    /// remain in the full run report and are excluded from replay identity.
    pub memory: GemmMemoryReceipt,
    /// Ordered NC/KC panel traversals through the caller's TilePool.
    pub panels: Vec<GemmPanelReceipt>,
}

/// Identity-stable projection of [`fs_la::GemmMemoryReport`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GemmMemoryReceipt {
    /// Caller-declared memory ceiling.
    pub limit_bytes: u64,
    /// Transactional C staging bytes.
    pub staging_bytes: u128,
    /// Shared B-pack bytes.
    pub b_pack_bytes: u128,
    /// M-band metadata bytes.
    pub band_metadata_bytes: u128,
    /// fs-la panel-receipt vector bytes.
    pub pool_run_bytes: u128,
    /// Fresh arena reservation per active worker.
    pub arena_bytes_per_worker: u64,
    /// Planned maximum active arena workers.
    pub active_arena_workers: usize,
    /// Planned active arena bytes.
    pub arena_bytes: u128,
    /// Checked fs-la-owned plan total.
    pub requested_bytes: u128,
}

impl From<fs_la::GemmMemoryReport> for GemmMemoryReceipt {
    fn from(report: fs_la::GemmMemoryReport) -> Self {
        Self {
            limit_bytes: report.limit_bytes,
            staging_bytes: report.staging_bytes,
            b_pack_bytes: report.b_pack_bytes,
            band_metadata_bytes: report.band_metadata_bytes,
            pool_run_bytes: report.pool_run_bytes,
            arena_bytes_per_worker: report.arena_bytes_per_worker,
            active_arena_workers: report.active_arena_workers,
            arena_bytes: report.arena_bytes,
            requested_bytes: report.requested_bytes,
        }
    }
}

impl GemmExecutionReceipt {
    /// Project a successful numerical run report onto identity-stable fields.
    /// Error paths retain their full report in [`GemmTuneError`] instead.
    #[must_use]
    pub fn from_report(report: &fs_la::GemmRunReport) -> Self {
        Self {
            declared_run: report.declared_run.0,
            completed_tiles: report.completed_tiles,
            total_tiles: report.total_tiles,
            memory: report.memory.into(),
            panels: report
                .pool_runs
                .iter()
                .map(|panel| GemmPanelReceipt {
                    kernel: panel.kernel.to_string(),
                    mode: panel.mode.to_string(),
                    declared_run: panel.declared_run.0,
                    completed: panel.completed,
                    total: panel.total,
                })
                .collect(),
        }
    }

    /// Whether every panel completed and carries the exact child RunId derived
    /// from this receipt's declared operation identity and panel ordinal.
    #[must_use]
    pub fn is_complete(&self) -> bool {
        if self.completed_tiles != self.total_tiles {
            return false;
        }
        if self.total_tiles == 0 {
            return self.panels.is_empty();
        }
        !self.panels.is_empty()
            && self.panels.iter().enumerate().all(|(ordinal, panel)| {
                let Ok(ordinal) = u64::try_from(ordinal) else {
                    return false;
                };
                panel.declared_run
                    == fs_la::gemm_panel_run_id(fs_exec::RunId(self.declared_run), ordinal).0
                    && !panel.kernel.is_empty()
                    && !panel.mode.is_empty()
                    && panel.total > 0
                    && panel.completed == panel.total
            })
    }
}

/// Explicit access policy for the durable GEMM tune cache.
///
/// A read-only caller may adopt a previously admitted row but cannot publish a
/// new measurement during speculative or not-yet-admitted work. Newly measured
/// rows remain available through [`GemmDispatch::new_tune_row`].
#[derive(Clone, Copy)]
pub enum GemmTuneCache<'a> {
    /// Do not read or write durable tuning state.
    Disabled,
    /// Adopt validated rows, but never write the ledger.
    ReadOnly(&'a Ledger),
    /// Adopt validated rows and persist a newly measured row before local
    /// installation.
    ReadWrite(&'a Ledger),
}

impl<'a> GemmTuneCache<'a> {
    fn reader(self) -> Option<&'a Ledger> {
        match self {
            Self::Disabled => None,
            Self::ReadOnly(ledger) | Self::ReadWrite(ledger) => Some(ledger),
        }
    }
}

/// A validated tune row that can be published after a wider admission gate.
///
/// Fields are private: callers can neither forge nor alter the scoped kernel,
/// shape, machine fingerprint, selected parameters, or measured evidence.
/// Instances are created only after fs-exec validates the complete row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValidatedGemmTuneRow {
    kernel: String,
    shape_class: String,
    machine: [u8; 8],
    params: String,
    measured: String,
    memory_limit_bytes: u64,
    probe_buffer_bytes: u128,
}

impl ValidatedGemmTuneRow {
    fn from_prepared(prepared: &PreparedGemmRow, machine: u64) -> Result<Self, GemmTuneError> {
        let execution = prepared.key().execution();
        let probe_buffer_bytes = probe_buffer_bytes_for_dims(execution.probe_dims()).ok_or(
            GemmTuneError::MemoryPlanOverflow {
                what: "tune-probe-buffers",
                limit_bytes: execution.memory_limit_bytes(),
            },
        )?;
        Ok(Self {
            kernel: prepared.key().kernel().to_string(),
            shape_class: prepared.key().shape_class().to_string(),
            machine: machine.to_le_bytes(),
            params: prepared.params_json(),
            measured: prepared.row_json(),
            memory_limit_bytes: execution.memory_limit_bytes(),
            probe_buffer_bytes,
        })
    }

    /// Canonical JSON preimage for this exact tune-table tuple.
    #[must_use]
    pub fn receipt_json(&self) -> String {
        use core::fmt::Write as _;

        let mut out = String::new();
        out.push_str("{\"kernel\":");
        push_json_string(&mut out, &self.kernel);
        out.push_str(",\"shape_class\":");
        push_json_string(&mut out, &self.shape_class);
        let _ = write!(
            out,
            ",\"machine\":\"{:016x}\",\"params\":",
            u64::from_le_bytes(self.machine)
        );
        out.push_str(&self.params);
        out.push_str(",\"measured\":");
        out.push_str(&self.measured);
        let _ = write!(
            out,
            ",\"memory_limit_bytes\":{},\"probe_buffer_bytes\":{}",
            self.memory_limit_bytes, self.probe_buffer_bytes
        );
        out.push('}');
        out
    }

    /// Domain-separated identity of [`Self::receipt_json`]. It is stable for a
    /// freshly measured row and the same row adopted later.
    #[must_use]
    pub fn receipt_identity(&self) -> fs_ledger::ContentHash {
        fs_blake3::hash_domain(GEMM_TUNE_ROW_RECEIPT_DOMAIN, self.receipt_json().as_bytes())
    }

    /// Whether this sealed row is the exact evidence behind one dispatched
    /// decision.
    #[must_use]
    pub fn matches_decision(
        &self,
        scoped_kernel: &str,
        shape_class: &str,
        machine: u64,
        canonical_plan: &str,
    ) -> bool {
        self.kernel == scoped_kernel
            && self.shape_class == shape_class
            && self.machine == machine.to_le_bytes()
            && self.params == format!("\"{canonical_plan}\"")
    }

    /// Whether a ledger query returned this exact sealed tuple.
    #[must_use]
    pub fn matches_ledger_row(&self, row: &fs_ledger::TuneRow) -> bool {
        self.kernel == row.kernel
            && self.shape_class == row.shape_class
            && self.machine.as_slice() == row.machine
            && self.params == row.params
            && self.measured == row.measured
    }

    /// Publish this already validated row without replacing a different row,
    /// preserving the caller's ledger transaction when one is active.
    ///
    /// # Errors
    /// Propagates the original ledger diagnostic.
    pub fn publish_to_ledger(&self, ledger: &Ledger) -> Result<(), fs_ledger::LedgerError> {
        self.publish_if_absent_or_identical(ledger)
    }

    /// Insert this row into an evidence ledger without replacing a different
    /// tune decision already stored under the same key. An identical row is an
    /// idempotent success; a conflict fails closed.
    ///
    /// # Errors
    /// Propagates ledger failures or returns [`fs_ledger::LedgerError::Invalid`]
    /// when the destination key contains a different tuple.
    pub fn publish_if_absent_or_identical(
        &self,
        ledger: &Ledger,
    ) -> Result<(), fs_ledger::LedgerError> {
        ledger.tune_put_if_absent(
            &self.kernel,
            &self.shape_class,
            &self.machine,
            &self.params,
            &self.measured,
        )?;
        let stored = ledger
            .tune_get(&self.kernel, &self.shape_class, &self.machine)?
            .ok_or_else(|| fs_ledger::LedgerError::Invalid {
                field: "tune".to_string(),
                problem: "insert-if-absent returned without a stored tune row".to_string(),
            })?;
        if self.matches_ledger_row(&stored) {
            Ok(())
        } else {
            Err(fs_ledger::LedgerError::Invalid {
                field: "tune".to_string(),
                problem: format!(
                    "refusing to replace a conflicting tune row for kernel {:?}, shape {:?}",
                    self.kernel, self.shape_class
                ),
            })
        }
    }

    /// Persist this already validated row to a durable tune ledger.
    ///
    /// # Errors
    /// Returns [`GemmTuneError::Ledger`] when the ledger refuses the write.
    pub fn persist(&self, ledger: &Ledger) -> Result<(), GemmTuneError> {
        self.publish_to_ledger(ledger)
            .map_err(|error| GemmTuneError::Ledger(error.to_string()))
    }

    /// Install this sealed row as the current mutable cache decision.
    ///
    /// Unlike [`ValidatedGemmTuneRow::publish_if_absent_or_identical`], this
    /// method deliberately replaces a stale or malformed row under the same
    /// cache key. The `tune` table is a dispatch cache, not an append-only
    /// evidence history; citable benchmark publication uses the insert-only
    /// method and content-addressed ledger artifacts instead.
    ///
    /// # Errors
    /// Returns [`GemmTuneError::Ledger`] when the upsert or exact read-back
    /// verification fails.
    pub fn replace_cache_row(&self, ledger: &Ledger) -> Result<(), GemmTuneError> {
        ledger
            .tune_put(
                &self.kernel,
                &self.shape_class,
                &self.machine,
                &self.params,
                &self.measured,
            )
            .map_err(|error| GemmTuneError::Ledger(error.to_string()))?;
        let stored = ledger
            .tune_get(&self.kernel, &self.shape_class, &self.machine)
            .map_err(|error| GemmTuneError::Ledger(error.to_string()))?
            .ok_or_else(|| {
                GemmTuneError::Ledger(
                    "cache upsert returned without a stored GEMM tune row".to_string(),
                )
            })?;
        if self.matches_ledger_row(&stored) {
            Ok(())
        } else {
            Err(GemmTuneError::Ledger(format!(
                "cache read-back disagrees with the sealed GEMM row for kernel {:?}, shape {:?}",
                self.kernel, self.shape_class
            )))
        }
    }
}

/// The kernel key for this build's GEMM accumulation contract.
#[must_use]
pub fn gemm_kernel_key() -> String {
    format!(
        "{GEMM_KERNEL_PREFIX}{}",
        fs_la::gemm::GEMM_BIT_SEMANTICS_VERSION
    )
}

/// Bucket one extent to its shape-class quantum (next power of two,
/// clamped to [8, 65536]).
fn bucket(extent: usize) -> usize {
    extent.clamp(8, 65_536).next_power_of_two()
}

/// The shape class for an (m, n, k) problem: power-of-two buckets. Exact
/// measured probe dims remain in [`GemmTuneKey`], so a bucket never erases
/// the context that produced a row.
#[must_use]
pub fn gemm_shape_class(m: usize, n: usize, k: usize) -> String {
    format!("m{}-n{}-k{}", bucket(m), bucket(n), bucket(k))
}

fn probe_dims(m: usize, n: usize, k: usize) -> [usize; 3] {
    [
        m.clamp(1, PROBE_MK_DIM_CAP),
        n.clamp(1, PROBE_N_DIM_CAP),
        k.clamp(1, PROBE_MK_DIM_CAP),
    ]
}

fn probe_buffer_bytes_for_dims([m, n, k]: [u64; 3]) -> Option<u128> {
    let m = u128::from(m);
    let n = u128::from(n);
    let k = u128::from(k);
    let elements = m
        .checked_mul(k)?
        .checked_add(k.checked_mul(n)?)?
        .checked_add(m.checked_mul(n)?.checked_mul(2)?)?;
    elements.checked_mul(core::mem::size_of::<u64>() as u128)
}

/// Construct the exact persistent tuning identity for this invocation.
/// Studies normally replay the recorded decision key directly; exposing this
/// constructor also lets admission and diagnostics explain why two calls do
/// or do not share evidence.
///
/// # Errors
/// [`GemmTuneError::Tune`] if a dimension or implementation identity cannot
/// be represented canonically.
pub fn gemm_tune_key(
    threads: usize,
    m: usize,
    n: usize,
    k: usize,
) -> Result<GemmTuneKey, GemmTuneError> {
    gemm_tune_key_budgeted(threads, m, n, k, fs_la::GemmMemoryEnvelope::UNBOUNDED)
}

/// Construct the persistent tuning identity under an explicit memory envelope.
/// Otherwise-identical calls with different ceilings cannot share rows or pins.
///
/// # Errors
/// As [`gemm_tune_key`].
pub fn gemm_tune_key_budgeted(
    threads: usize,
    m: usize,
    n: usize,
    k: usize,
    envelope: fs_la::GemmMemoryEnvelope,
) -> Result<GemmTuneKey, GemmTuneError> {
    let pool = TilePool::for_host(threads, SESSION_GEMM_POOL_SEED);
    gemm_tune_key_with_pool_budgeted(&pool, m, n, k, envelope)
}

/// Construct the persistent tuning identity from the TilePool that will
/// actually execute the measured and selected plans.
///
/// # Errors
/// [`GemmTuneError::Tune`] if a pool dimension or identity component cannot
/// be represented canonically.
pub fn gemm_tune_key_with_pool(
    pool: &TilePool,
    m: usize,
    n: usize,
    k: usize,
) -> Result<GemmTuneKey, GemmTuneError> {
    gemm_tune_key_with_pool_budgeted(pool, m, n, k, fs_la::GemmMemoryEnvelope::UNBOUNDED)
}

/// Construct the persistent tuning identity from the executing pool and an
/// explicit memory envelope.
///
/// # Errors
/// As [`gemm_tune_key_with_pool`].
pub fn gemm_tune_key_with_pool_budgeted(
    pool: &TilePool,
    m: usize,
    n: usize,
    k: usize,
    envelope: fs_la::GemmMemoryEnvelope,
) -> Result<GemmTuneKey, GemmTuneError> {
    gemm_tune_key_for_execution(
        pool.workers(),
        pool.workers(),
        envelope.limit_bytes,
        &pool.placement_identity(),
        m,
        n,
        k,
    )
}

fn gemm_tune_key_for_execution(
    requested_threads: usize,
    thread_budget: usize,
    memory_limit_bytes: u64,
    placement: &str,
    m: usize,
    n: usize,
    k: usize,
) -> Result<GemmTuneKey, GemmTuneError> {
    gemm_tune_key_for_execution_schema(
        requested_threads,
        thread_budget,
        memory_limit_bytes,
        placement,
        m,
        n,
        k,
        GEMM_TUNER_SCHEMA_VERSION,
    )
}

#[allow(clippy::too_many_arguments)]
fn gemm_tune_key_for_execution_schema(
    requested_threads: usize,
    thread_budget: usize,
    memory_limit_bytes: u64,
    placement: &str,
    m: usize,
    n: usize,
    k: usize,
    tuner_schema: u32,
) -> Result<GemmTuneKey, GemmTuneError> {
    debug_assert!(tuner_schema > 0);
    let implementation = format!(
        "fs-la-{}-gemm-v{}-fs-session-tuner-v{tuner_schema}",
        fs_la::VERSION,
        fs_la::GEMM_IMPLEMENTATION_VERSION
    );
    let execution = GemmExecutionIdentity::new(
        requested_threads,
        thread_budget,
        memory_limit_bytes,
        probe_dims(m, n, k),
        fs_la::gemm_execution_tier(),
        placement,
        implementation,
        fs_la::gemm_build_identity(),
    )?;
    Ok(GemmTuneKey::new(
        gemm_kernel_key(),
        gemm_shape_class(m, n, k),
        execution,
    )?)
}

#[track_caller]
fn checked_product(label: &str, lhs: usize, rhs: usize) -> usize {
    lhs.checked_mul(rhs)
        .unwrap_or_else(|| panic!("{label} extent overflow: {lhs} * {rhs}"))
}

fn try_filled_buffer<T: Copy>(
    len: usize,
    value: T,
    what: &'static str,
    envelope: fs_la::GemmMemoryEnvelope,
    peak_used_bytes: u128,
) -> Result<Vec<T>, GemmTuneError> {
    let requested_bytes = (len as u128)
        .checked_mul(core::mem::size_of::<T>() as u128)
        .ok_or(GemmTuneError::MemoryPlanOverflow {
            what,
            limit_bytes: envelope.limit_bytes,
        })?;
    let mut values = Vec::new();
    values
        .try_reserve_exact(len)
        .map_err(|_| GemmTuneError::MemoryRefused {
            what,
            requested_bytes,
            limit_bytes: envelope.limit_bytes,
            peak_used_bytes,
            report: None,
        })?;
    values.resize(len, value);
    Ok(values)
}

/// Mirror fs-la's public contiguous-slice precondition before consulting or
/// mutating tuning state. fs-la validates again at the execution boundary;
/// this ordering is the session-level no-phantom-row guarantee.
#[track_caller]
fn assert_contiguous_shapes(m: usize, n: usize, k: usize, a: &[f64], b: &[f64], c: &[f64]) {
    let a_len = checked_product("a", m, k);
    let b_len = checked_product("b", k, n);
    let c_len = checked_product("c", m, n);
    assert_eq!(a.len(), a_len, "a must be m*k = {a_len}");
    assert_eq!(b.len(), b_len, "b must be k*n = {b_len}");
    assert_eq!(c.len(), c_len, "c must be m*n = {c_len}");
}

/// Deterministic probe fill (splitmix64 bits folded to [-0.5, 0.5)):
/// integer-only, so probe inputs are bit-identical on every ISA.
fn probe_fill(buf: &mut [f64], salt: u64) {
    for (i, slot) in buf.iter_mut().enumerate() {
        let mut z = (i as u64)
            .wrapping_add(salt)
            .wrapping_add(0x9E37_79B9_7F4A_7C15);
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        // 53 mantissa bits → [0, 1), then center.
        *slot = (z >> 11) as f64 / 9_007_199_254_740_992.0 - 0.5;
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SweepCandidate {
    plan: GemmBlockPlan,
    effective_mc: usize,
    effective_nc: usize,
}

#[derive(Debug)]
struct SweepResult {
    winner: GemmBlockPlan,
    evidence: TuneEvidence,
}

/// Build the lattice the kernel will ACTUALLY execute. Nominal plans that
/// collapse to the same clamped `(mc, nc)` pair are measured only once.
fn effective_sweep_candidates(pm: usize, pn: usize) -> Result<Vec<SweepCandidate>, GemmTuneError> {
    let mut seen = std::collections::BTreeSet::new();
    let mut candidates = Vec::with_capacity(SWEEP_MC.len() * SWEEP_NC_CAP.len());
    for (mc, nc_cap) in SWEEP_MC
        .iter()
        .flat_map(|&mc| SWEEP_NC_CAP.iter().map(move |&nc| (mc, nc)))
    {
        let plan = GemmBlockPlan::new(mc, nc_cap)?;
        // These are the clamps applied by fs-la's packed parallel engine.
        let effective_mc = plan.mc.max(8).min(pm.max(8));
        let effective_nc = pn.min(plan.nc_cap).max(4);
        if seen.insert((effective_mc, effective_nc)) {
            candidates.push(SweepCandidate {
                plan,
                effective_mc,
                effective_nc,
            });
        }
    }
    Ok(candidates)
}

/// Measure candidate executions supplied by `run`. Keeping this core
/// injectable lets the Gauntlet force drift in each repeat and cache faults
/// without adding test behavior to the production GEMM implementation.
fn measure_candidates<R>(
    gate: &CancelGate,
    candidates: &[SweepCandidate],
    output_len: usize,
    envelope: fs_la::GemmMemoryEnvelope,
    base_used_bytes: u128,
    mut run: R,
) -> Result<SweepResult, GemmTuneError>
where
    R: FnMut(&SweepCandidate, &mut [f64]) -> Result<u128, GemmTuneError>,
{
    let output_bytes = (output_len as u128)
        .checked_mul(core::mem::size_of::<u64>() as u128)
        .ok_or(GemmTuneError::MemoryPlanOverflow {
            what: "tune-probe-output",
            limit_bytes: envelope.limit_bytes,
        })?;
    let first_output_peak =
        base_used_bytes
            .checked_add(output_bytes)
            .ok_or(GemmTuneError::MemoryPlanOverflow {
                what: "tune-probe-output-peak",
                limit_bytes: envelope.limit_bytes,
            })?;
    let live_probe_bytes =
        first_output_peak
            .checked_add(output_bytes)
            .ok_or(GemmTuneError::MemoryPlanOverflow {
                what: "tune-probe-reference-peak",
                limit_bytes: envelope.limit_bytes,
            })?;
    let mut c = try_filled_buffer(
        output_len,
        0.0_f64,
        "tune-probe-c",
        envelope,
        base_used_bytes,
    )?;
    let mut reference_bits = try_filled_buffer(
        output_len,
        0_u64,
        "tune-probe-reference",
        envelope,
        first_output_peak,
    )?;
    let mut observations = Vec::with_capacity(candidates.len());
    let mut ranked: Vec<(u64, usize, GemmBlockPlan)> = Vec::with_capacity(candidates.len());
    let mut reference_initialized = false;
    let mut numerical_peak = 0_u128;
    for (index, candidate) in candidates.iter().enumerate() {
        if gate.is_requested() {
            return Err(cancelled_with_live_probe_memory(
                envelope,
                live_probe_bytes,
                numerical_peak,
            ));
        }
        let mut samples_ns = Vec::with_capacity(SWEEP_SAMPLES);
        for repeat in 1..=SWEEP_SAMPLES {
            if gate.is_requested() {
                return Err(cancelled_with_live_probe_memory(
                    envelope,
                    live_probe_bytes,
                    numerical_peak,
                ));
            }
            c.fill(0.0);
            let t0 = std::time::Instant::now();
            numerical_peak = numerical_peak.max(run(candidate, &mut c)?);
            let ns = u64::try_from(t0.elapsed().as_nanos()).unwrap_or(u64::MAX);
            samples_ns.push(ns.max(1));
            if gate.is_requested() {
                return Err(cancelled_with_live_probe_memory(
                    envelope,
                    live_probe_bytes,
                    numerical_peak,
                ));
            }

            // Compare every output word directly. A fixed-width digest is
            // not a proof of bit-neutrality and would also hide which repeat
            // drifted. `to_bits` intentionally distinguishes signed zero and
            // every NaN payload.
            if !reference_initialized {
                for (dst, value) in reference_bits.iter_mut().zip(&c) {
                    *dst = value.to_bits();
                }
                reference_initialized = true;
            } else if !reference_bits
                .iter()
                .zip(&c)
                .all(|(&expected, value)| expected == value.to_bits())
            {
                return Err(GemmTuneError::BitDrift {
                    candidate: candidate.plan.canonical(),
                    repeat,
                });
            }
        }
        let best = samples_ns.iter().copied().min().unwrap_or(u64::MAX);
        ranked.push((best, index, candidate.plan));
        observations.push(TuneObservation::wall_time(
            candidate.plan.canonical(),
            samples_ns,
        )?);
    }
    ranked.sort_unstable_by_key(|&(ns, index, _)| (ns, index));
    let winner = ranked
        .first()
        .map(|entry| entry.2)
        .ok_or_else(|| TuneError {
            detail: "the effective GEMM candidate lattice is empty".to_string(),
        })?;
    let evidence = TuneEvidence::ranked_wall_times(observations)?;
    Ok(SweepResult { winner, evidence })
}

/// Run the bounded candidate sweep for one exact probe. This function only
/// measures and validates; its caller persists first and commits the tuner
/// row second so a cache failure cannot leave a phantom in-memory success.
fn run_sweep(
    gate: &CancelGate,
    pool: &TilePool,
    declared_run: fs_exec::RunId,
    m: usize,
    n: usize,
    k: usize,
    envelope: fs_la::GemmMemoryEnvelope,
) -> Result<SweepResult, GemmTuneError> {
    // Probe at the CALLER's dims (capped): the oracle lane showed that
    // probing at the class's power-of-two bucket flips winners — at
    // m = 320 the band count under each mc differs from m = 512, and
    // band balance decides the ranking. The row retains the bucketed shape
    // class, but the exact capped probe is also part of the scoped key so a
    // neighboring caller cannot silently inherit different evidence.
    let [pm, pn, pk] = probe_dims(m, n, k);
    let probe_dims_u64 = [pm, pn, pk]
        .map(|extent| u64::try_from(extent).expect("capped GEMM probe dimensions fit u64"));
    let probe_buffer_bytes =
        probe_buffer_bytes_for_dims(probe_dims_u64).ok_or(GemmTuneError::MemoryPlanOverflow {
            what: "tune-probe-buffers",
            limit_bytes: envelope.limit_bytes,
        })?;
    if probe_buffer_bytes > u128::from(envelope.limit_bytes) {
        return Err(GemmTuneError::MemoryRefused {
            what: "tune-probe-envelope",
            requested_bytes: probe_buffer_bytes,
            limit_bytes: envelope.limit_bytes,
            peak_used_bytes: 0,
            report: None,
        });
    }
    let child_limit_bytes = if envelope == fs_la::GemmMemoryEnvelope::UNBOUNDED {
        u64::MAX
    } else {
        u64::try_from(u128::from(envelope.limit_bytes) - probe_buffer_bytes)
            .expect("bounded probe preflight leaves a u64 child envelope")
    };
    let child_envelope = fs_la::GemmMemoryEnvelope {
        limit_bytes: child_limit_bytes,
    };

    let a_len = checked_product("tune probe A", pm, pk);
    let b_len = checked_product("tune probe B", pk, pn);
    let a_bytes = (a_len as u128)
        .checked_mul(core::mem::size_of::<f64>() as u128)
        .ok_or(GemmTuneError::MemoryPlanOverflow {
            what: "tune-probe-a",
            limit_bytes: envelope.limit_bytes,
        })?;
    let b_bytes = (b_len as u128)
        .checked_mul(core::mem::size_of::<f64>() as u128)
        .ok_or(GemmTuneError::MemoryPlanOverflow {
            what: "tune-probe-b",
            limit_bytes: envelope.limit_bytes,
        })?;
    let ab_bytes = a_bytes
        .checked_add(b_bytes)
        .ok_or(GemmTuneError::MemoryPlanOverflow {
            what: "tune-probe-a-plus-b",
            limit_bytes: envelope.limit_bytes,
        })?;
    let mut a = try_filled_buffer(a_len, 0.0_f64, "tune-probe-a", envelope, 0)?;
    let mut b = try_filled_buffer(b_len, 0.0_f64, "tune-probe-b", envelope, a_bytes)?;
    probe_fill(&mut a, 0xA);
    probe_fill(&mut b, 0xB);
    let candidates = effective_sweep_candidates(pm, pn)?;
    let mut sweep_ordinal = 0_u64;
    measure_candidates(
        gate,
        &candidates,
        checked_product("tune probe C", pm, pn),
        envelope,
        ab_bytes,
        |candidate, c| {
            let sweep_run = declared_run.derive(GEMM_SWEEP_RUN_DOMAIN, sweep_ordinal);
            sweep_ordinal = sweep_ordinal.checked_add(1).ok_or_else(|| TuneError {
                detail: "GEMM sweep run ordinal exhausted".to_string(),
            })?;
            fs_la::gemm_f64_parallel_with_pool_budgeted(
                pm,
                pn,
                pk,
                1.0,
                &a,
                &b,
                0.0,
                c,
                pool,
                candidate.effective_mc,
                candidate.effective_nc,
                gate,
                sweep_run,
                child_envelope,
            )
            .map(|report| report.memory.peak_used_bytes)
            .map_err(|error| gemm_error_with_session_memory(error, envelope, probe_buffer_bytes))
        },
    )
}

/// Persist a validated measured row before installing it in the process-local
/// tuner. `persist` is injectable so the failure-atomic boundary is directly
/// testable without corrupting a real ledger connection.
fn install_sweep_row<P>(
    tuner: &mut Tuner,
    key: &GemmTuneKey,
    sweep: SweepResult,
    persist: P,
) -> Result<(GemmBlockPlan, ValidatedGemmTuneRow), GemmTuneError>
where
    P: FnOnce(&ValidatedGemmTuneRow) -> Result<(), GemmTuneError>,
{
    let prepared = tuner.prepare_gemm_row(key, sweep.winner, sweep.evidence)?;
    let validated = ValidatedGemmTuneRow::from_prepared(&prepared, tuner.machine())?;
    persist(&validated)?;
    let winner = sweep.winner;
    tuner.commit_gemm_row(prepared)?;
    Ok((winner, validated))
}

fn adopt_cached_row(
    tuner: &mut Tuner,
    key: &GemmTuneKey,
    params: &str,
    measured: &str,
) -> Result<Option<ValidatedGemmTuneRow>, GemmTuneError> {
    let Ok(prepared) = tuner.prepare_adopt_gemm_row_json(key, measured) else {
        return Ok(None);
    };
    if params != prepared.params_json() {
        return Ok(None);
    }
    let validated = ValidatedGemmTuneRow::from_prepared(&prepared, tuner.machine())?;
    tuner.commit_gemm_row(prepared)?;
    Ok(Some(validated))
}

fn execute_prepared_decision<R, F>(
    tuner: &mut Tuner,
    decision: PreparedGemmDecision,
    run: F,
) -> Result<(GemmBlockPlan, TuneSource, R), GemmTuneError>
where
    F: FnOnce(GemmBlockPlan) -> Result<R, GemmTuneError>,
{
    let plan = decision.plan();
    let source = decision.source();
    let output = run(plan)?;
    // Exclusive access to `tuner` spans prepare -> run -> commit, so no
    // applicable pin/row can change and make this prepared decision stale.
    tuner
        .commit_gemm_decision(decision)
        .expect("exclusive tuner borrow preserves a prepared GEMM decision");
    Ok((plan, source, output))
}

/// The production autotuned f64 GEMM: `c = alpha·a·b + beta·c` through
/// the measure → cache → model → dispatch loop.
///
/// Resolution order after shape and cancellation preflight: a pinned plan
/// dispatches without measurement; else an exact cached row (in the tuner,
/// seeded from the cache when permitted); else the bounded sweep measures,
/// applies the explicit write policy, commits it locally, and dispatches.
/// Serial, small-M, and no-product calls bypass tuning entirely.
///
/// # Errors
/// [`GemmTuneError`] — cancellation, tuner refusals, ledger I/O, or a
/// bit-neutrality violation. On every returned error, `c` retains its exact
/// original bits. Cancellable GEMM computes in private staging, drains its
/// workers, and commits only after its final poll.
///
/// # Panics
/// Inherits fs-la's structured shape panics for mismatched slice
/// lengths.
#[allow(clippy::too_many_arguments)] // BLAS-shape signature + orchestration handles
pub fn gemm_f64_session(
    tuner: &mut Tuner,
    cache: GemmTuneCache<'_>,
    gate: &CancelGate,
    threads: usize,
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    beta: f64,
    c: &mut [f64],
) -> Result<GemmDispatch, GemmTuneError> {
    gemm_f64_session_budgeted(
        tuner,
        cache,
        gate,
        threads,
        m,
        n,
        k,
        alpha,
        a,
        b,
        beta,
        c,
        fs_la::GemmMemoryEnvelope::UNBOUNDED,
    )
}

/// As [`gemm_f64_session`], under an explicit memory envelope bound into tune
/// identity and every numerical dispatch.
///
/// # Errors
/// As [`gemm_f64_session`], plus structured memory refusal.
///
/// # Panics
/// As [`gemm_f64_session`].
#[allow(clippy::too_many_arguments)]
pub fn gemm_f64_session_budgeted(
    tuner: &mut Tuner,
    cache: GemmTuneCache<'_>,
    gate: &CancelGate,
    threads: usize,
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    beta: f64,
    c: &mut [f64],
    envelope: fs_la::GemmMemoryEnvelope,
) -> Result<GemmDispatch, GemmTuneError> {
    let pool = TilePool::for_host(threads, SESSION_GEMM_POOL_SEED);
    gemm_f64_session_with_pool_budgeted(
        tuner, cache, &pool, gate, m, n, k, alpha, a, b, beta, c, envelope,
    )
}

/// The production autotuned f64 GEMM on a caller-owned, reusable TilePool.
/// The same pool executes every sweep candidate and the selected plan; its
/// normalized worker budget and placement policy are bound into the tune key.
///
/// # Errors
/// As [`gemm_f64_session`], plus a structured executor failure if TilePool
/// contains a tile panic or detects an incomplete traversal.
///
/// # Panics
/// Inherits fs-la's structured shape panics for mismatched slice lengths.
#[allow(clippy::too_many_arguments)]
pub fn gemm_f64_session_with_pool(
    tuner: &mut Tuner,
    cache: GemmTuneCache<'_>,
    pool: &TilePool,
    gate: &CancelGate,
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    beta: f64,
    c: &mut [f64],
) -> Result<GemmDispatch, GemmTuneError> {
    gemm_f64_session_with_pool_budgeted(
        tuner,
        cache,
        pool,
        gate,
        m,
        n,
        k,
        alpha,
        a,
        b,
        beta,
        c,
        fs_la::GemmMemoryEnvelope::UNBOUNDED,
    )
}

/// As [`gemm_f64_session_with_pool`], under an explicit memory envelope.
///
/// # Errors
/// As [`gemm_f64_session_with_pool`], plus structured memory refusal.
///
/// # Panics
/// As [`gemm_f64_session_with_pool`].
#[allow(clippy::too_many_arguments)]
pub fn gemm_f64_session_with_pool_budgeted(
    tuner: &mut Tuner,
    cache: GemmTuneCache<'_>,
    pool: &TilePool,
    gate: &CancelGate,
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    beta: f64,
    c: &mut [f64],
    envelope: fs_la::GemmMemoryEnvelope,
) -> Result<GemmDispatch, GemmTuneError> {
    gemm_f64_session_with_pool_declared_budgeted(
        tuner,
        cache,
        pool,
        gate,
        fs_exec::RunId::default(),
        m,
        n,
        k,
        alpha,
        a,
        b,
        beta,
        c,
        envelope,
    )
}

/// As [`gemm_f64_session_with_pool`], with the caller-ledgered identity of the
/// final production dispatch. Sweep repetitions receive separate
/// domain-derived children and cannot collide with the final run's tile
/// streams.
///
/// # Errors
/// As [`gemm_f64_session_with_pool`].
///
/// # Panics
/// As [`gemm_f64_session_with_pool`].
#[allow(clippy::too_many_arguments)]
pub fn gemm_f64_session_with_pool_declared(
    tuner: &mut Tuner,
    cache: GemmTuneCache<'_>,
    pool: &TilePool,
    gate: &CancelGate,
    declared_run: fs_exec::RunId,
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    beta: f64,
    c: &mut [f64],
) -> Result<GemmDispatch, GemmTuneError> {
    gemm_f64_session_with_pool_declared_budgeted(
        tuner,
        cache,
        pool,
        gate,
        declared_run,
        m,
        n,
        k,
        alpha,
        a,
        b,
        beta,
        c,
        fs_la::GemmMemoryEnvelope::UNBOUNDED,
    )
}

/// As [`gemm_f64_session_with_pool_declared`], under an explicit memory
/// envelope bound into tune identity, sweep admission, and final dispatch.
///
/// # Errors
/// As [`gemm_f64_session_with_pool_declared`], plus structured memory refusal.
///
/// # Panics
/// As [`gemm_f64_session_with_pool_declared`].
#[allow(clippy::too_many_arguments)]
pub fn gemm_f64_session_with_pool_declared_budgeted(
    tuner: &mut Tuner,
    cache: GemmTuneCache<'_>,
    pool: &TilePool,
    gate: &CancelGate,
    declared_run: fs_exec::RunId,
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    beta: f64,
    c: &mut [f64],
    envelope: fs_la::GemmMemoryEnvelope,
) -> Result<GemmDispatch, GemmTuneError> {
    // Public slice/extent preconditions are checked before tier resolution,
    // cache reads, sweeps, rows, or decisions. Invalid work cannot poison the
    // tuning state and `c` is still untouched when this panics.
    assert_contiguous_shapes(m, n, k, a, b, c);
    if gate.is_requested() {
        return Err(cancelled_before_compute(envelope));
    }

    let key = gemm_tune_key_with_pool_budgeted(pool, m, n, k, envelope)?;
    let kernel = key.kernel().to_string();
    let shape_class = gemm_shape_class(m, n, k);
    let mut swept = false;
    let mut new_tune_row = None;
    let mut validated_tune_row = None;

    // No product, one-thread, and small-M routes do not have a meaningful
    // production MC/NC choice. Dispatch them cancellation-correctly under the
    // documented cold plan without reading or mutating tune state.
    if !fs_la::gemm_tuning_is_effective(m, n, k, alpha, pool.workers()) {
        let plan = GemmBlockPlan::COLD_START;
        let run = fs_la::gemm_f64_parallel_with_pool_budgeted(
            m,
            n,
            k,
            alpha,
            a,
            b,
            beta,
            c,
            pool,
            plan.mc,
            n.min(plan.nc_cap).max(1),
            gate,
            declared_run,
            envelope,
        )
        .map_err(GemmTuneError::from)?;
        return Ok(GemmDispatch {
            kernel,
            shape_class,
            plan,
            source: TuneSource::ColdStart,
            swept,
            new_tune_row,
            validated_tune_row,
            run,
        });
    }

    if !tuner.has_gemm_pin(&key) && !tuner.has_gemm_row(&key) {
        // Cache tier: try the ledger before measuring. Stale
        // (other-machine) or non-canonical rows are refused by
        // prepare_adopt_gemm_row_json and we fall through to a fresh sweep.
        // The ledger's separate params column must agree byte-for-byte with
        // the validated row body before either is allowed into the tuner.
        if let Some(ledger) = cache.reader() {
            let cached = ledger
                .tune_get(
                    key.kernel(),
                    key.shape_class(),
                    &tuner.machine().to_le_bytes(),
                )
                .map_err(|e| GemmTuneError::Ledger(e.to_string()))?;
            if let Some(row) = cached {
                validated_tune_row = adopt_cached_row(tuner, &key, &row.params, &row.measured)?;
            }
        }
        if validated_tune_row.is_none() {
            let sweep = run_sweep(gate, pool, declared_run, m, n, k, envelope)?;
            swept = true;
            let (_, validated) = install_sweep_row(tuner, &key, sweep, |row| match cache {
                GemmTuneCache::ReadWrite(ledger) => row.replace_cache_row(ledger),
                GemmTuneCache::Disabled | GemmTuneCache::ReadOnly(_) => Ok(()),
            })?;
            new_tune_row = Some(validated.clone());
            validated_tune_row = Some(validated);
        }
    }

    if gate.is_requested() {
        return Err(cancelled_before_compute(envelope));
    }
    let decision = tuner.prepare_gemm_decision(&key);
    let (plan, source, run) = execute_prepared_decision(tuner, decision, |plan| {
        fs_la::gemm_f64_parallel_with_pool_budgeted(
            m,
            n,
            k,
            alpha,
            a,
            b,
            beta,
            c,
            pool,
            plan.mc,
            n.min(plan.nc_cap).max(1),
            gate,
            declared_run,
            envelope,
        )
        .map_err(GemmTuneError::from)
    })?;
    Ok(GemmDispatch {
        kernel,
        shape_class,
        plan,
        source,
        swept,
        new_tune_row,
        validated_tune_row,
        run,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tuner_schema_bump_separates_durable_keys() {
        let v1 =
            gemm_tune_key_for_execution_schema(4, 4, u64::MAX, "test-placement", 512, 640, 512, 1)
                .expect("schema-v1 key");
        let v2 =
            gemm_tune_key_for_execution_schema(4, 4, u64::MAX, "test-placement", 512, 640, 512, 2)
                .expect("schema-v2 key");

        assert_ne!(v1.kernel(), v2.kernel());
        assert_eq!(v1.shape_class(), v2.shape_class());
        assert!(v1.kernel().contains("fs-session-tuner-v1"));
        assert!(v2.kernel().contains("fs-session-tuner-v2"));
    }

    #[test]
    fn execution_receipt_excludes_schedule_measurements() {
        let operation_run = fs_exec::RunId(7);
        let no_product = fs_la::GemmRunReport {
            declared_run: operation_run,
            completed_tiles: 0,
            total_tiles: 0,
            pool_runs: Vec::new(),
            memory: fs_la::GemmMemoryReport::default(),
        };
        assert!(
            GemmExecutionReceipt::from_report(&no_product).is_complete(),
            "a successful no-product dispatch is complete without panel traversals"
        );
        let base_panel = fs_exec::RunReport {
            kernel: "fs-la/gemm-f64-m-band-v1",
            mode: "deterministic",
            declared_run: fs_la::gemm_panel_run_id(operation_run, 0),
            completed: 4,
            total: 4,
            steals: 0,
            cross_ccd_steals: 0,
            cancel_latencies_ns: Vec::new(),
            tiles_by_worker: vec![2, 2],
        };
        let first = fs_la::GemmRunReport {
            declared_run: operation_run,
            completed_tiles: 32,
            total_tiles: 32,
            pool_runs: vec![base_panel.clone()],
            memory: fs_la::GemmMemoryReport::default(),
        };
        let mut noisy_panel = base_panel;
        noisy_panel.steals = 99;
        noisy_panel.cross_ccd_steals = 17;
        noisy_panel.cancel_latencies_ns = vec![3, 5, 8];
        noisy_panel.tiles_by_worker = vec![4, 0];
        let mut second = fs_la::GemmRunReport {
            declared_run: operation_run,
            completed_tiles: 32,
            total_tiles: 32,
            pool_runs: vec![noisy_panel],
            memory: fs_la::GemmMemoryReport::default(),
        };
        second.memory.peak_used_bytes = 999;
        second.memory.refused_bytes = 17;
        assert_eq!(
            GemmExecutionReceipt::from_report(&first),
            GemmExecutionReceipt::from_report(&second),
            "steal, latency, and worker-distribution envelopes are not replay identity"
        );
        assert!(GemmExecutionReceipt::from_report(&first).is_complete());
        let mut different_memory_plan = second.clone();
        different_memory_plan.memory.limit_bytes = 1 << 20;
        assert_ne!(
            GemmExecutionReceipt::from_report(&first),
            GemmExecutionReceipt::from_report(&different_memory_plan),
            "the declared memory plan is replay identity"
        );
        let mut different_run = first;
        different_run.pool_runs[0].declared_run = fs_exec::RunId(1);
        assert_ne!(
            GemmExecutionReceipt::from_report(&different_run),
            GemmExecutionReceipt::from_report(&second),
            "declared logical run is part of replay identity"
        );
    }

    fn synthetic_sweep() -> SweepResult {
        let winner = GemmBlockPlan::new(16, 512).expect("winner plan");
        let runner_up = GemmBlockPlan::new(32, 512).expect("runner-up plan");
        let evidence = TuneEvidence::ranked_wall_times(vec![
            TuneObservation::wall_time(winner.canonical(), vec![10, 11, 12])
                .expect("winner evidence"),
            TuneObservation::wall_time(runner_up.canonical(), vec![20, 21, 22])
                .expect("runner-up evidence"),
        ])
        .expect("ranked evidence");
        SweepResult { winner, evidence }
    }

    #[test]
    fn exact_bits_gate_catches_drift_in_every_repeat() {
        let candidates = effective_sweep_candidates(320, 2048).expect("candidate lattice");
        assert!(candidates.len() >= 2);
        for drift_repeat in 1..=SWEEP_SAMPLES {
            let mut call = 0usize;
            let error = measure_candidates(
                &CancelGate::new(),
                &candidates,
                2,
                fs_la::GemmMemoryEnvelope::UNBOUNDED,
                0,
                |_, c| {
                    let candidate = call / SWEEP_SAMPLES;
                    let repeat = call % SWEEP_SAMPLES + 1;
                    call += 1;
                    c[0] = 0.0;
                    c[1] = f64::from_bits(0x7ff8_0000_0000_0001);
                    if candidate == 1 && repeat == drift_repeat {
                        // Both changes are invisible to ordinary floating-point
                        // equality: signed zero compares equal and NaNs compare
                        // unequal regardless of payload. The contract is bits.
                        c[0] = -0.0;
                        c[1] = f64::from_bits(0x7ff8_0000_0000_0002);
                    }
                    Ok(0)
                },
            )
            .expect_err("the injected repeat must fail closed");
            assert!(
                matches!(
                    error,
                    GemmTuneError::BitDrift {
                        repeat,
                        ..
                    } if repeat == drift_repeat
                ),
                "repeat {drift_repeat}: {error}"
            );
        }
    }

    #[test]
    fn effective_candidate_lattice_is_unique_and_exercises_nc() {
        let narrow = effective_sweep_candidates(320, 288).expect("narrow lattice");
        assert_eq!(narrow.len(), SWEEP_MC.len());
        let narrow_pairs: std::collections::BTreeSet<_> = narrow
            .iter()
            .map(|candidate| (candidate.effective_mc, candidate.effective_nc))
            .collect();
        assert_eq!(narrow_pairs.len(), narrow.len());

        let wide = effective_sweep_candidates(320, 2048).expect("wide lattice");
        let wide_pairs: std::collections::BTreeSet<_> = wide
            .iter()
            .map(|candidate| (candidate.effective_mc, candidate.effective_nc))
            .collect();
        assert_eq!(wide_pairs.len(), wide.len());
        assert_eq!(
            wide.iter()
                .map(|candidate| candidate.effective_nc)
                .collect::<std::collections::BTreeSet<_>>(),
            std::collections::BTreeSet::from([512, 2048]),
            "n > 512 must measure both a multi-panel NC=512 execution and the wider panel"
        );

        let mut executed = Vec::new();
        measure_candidates(
            &CancelGate::new(),
            &wide,
            1,
            fs_la::GemmMemoryEnvelope::UNBOUNDED,
            0,
            |candidate, c| {
                executed.push((candidate.effective_mc, candidate.effective_nc));
                c[0] = 1.0;
                Ok(0)
            },
        )
        .expect("synthetic sweep");
        for pair in wide_pairs {
            assert_eq!(
                executed
                    .iter()
                    .filter(|&&observed| observed == pair)
                    .count(),
                SWEEP_SAMPLES,
                "each unique effective pair runs every repeat"
            );
        }
    }

    #[test]
    fn cancellation_between_repeats_returns_no_partial_evidence() {
        let candidates = effective_sweep_candidates(320, 2048).expect("candidate lattice");
        let gate = CancelGate::new();
        let error = measure_candidates(
            &gate,
            &candidates,
            1,
            fs_la::GemmMemoryEnvelope::UNBOUNDED,
            0,
            |_, c| {
                c[0] = 1.0;
                gate.request();
                Ok(0)
            },
        )
        .expect_err("the post-repeat poll must observe cancellation");
        assert!(matches!(
            error,
            GemmTuneError::Cancelled {
                peak_used_bytes: 16,
                report: None,
                ..
            }
        ));
    }

    #[test]
    fn cache_persistence_failure_is_atomic_and_retryable() {
        let key = gemm_tune_key(4, 320, 288, 300).expect("key");
        let mut tuner = Tuner::cold(0xAA55);
        let error = install_sweep_row(&mut tuner, &key, synthetic_sweep(), |_| {
            Err(GemmTuneError::Ledger("injected write failure".to_string()))
        })
        .expect_err("faulted cache write");
        assert!(matches!(error, GemmTuneError::Ledger(_)));
        assert!(!tuner.has_gemm_row(&key));
        assert!(tuner.decisions().is_empty());

        let (winner, _) = install_sweep_row(&mut tuner, &key, synthetic_sweep(), |_| Ok(()))
            .expect("retry installs the row");
        assert_eq!(winner, GemmBlockPlan::new(16, 512).expect("plan"));
        assert!(tuner.has_gemm_row(&key));
    }

    #[test]
    fn cached_params_and_body_must_agree_before_adoption() {
        let key = gemm_tune_key(4, 320, 288, 300).expect("key");
        let mut producer = Tuner::cold(0xAA55);
        let mut params = String::new();
        let mut measured = String::new();
        let mut sealed = None;
        install_sweep_row(&mut producer, &key, synthetic_sweep(), |validated| {
            params.clone_from(&validated.params);
            measured.clone_from(&validated.measured);
            sealed = Some(validated.clone());
            Ok(())
        })
        .expect("produce cached row");
        let sealed = sealed.expect("sealed row");

        let mut consumer = Tuner::cold(0xAA55);
        assert!(
            adopt_cached_row(&mut consumer, &key, "\"mc=32,nc-cap=512\"", &measured)
                .expect("mismatch is a cache miss")
                .is_none()
        );
        assert!(!consumer.has_gemm_row(&key));
        let adopted = adopt_cached_row(&mut consumer, &key, &params, &measured)
            .expect("adopt")
            .expect("validated adopted row");
        assert_eq!(adopted.receipt_identity(), sealed.receipt_identity());
        assert!(adopted.matches_decision(
            key.kernel(),
            key.shape_class(),
            0xAA55,
            "mc=16,nc-cap=512"
        ));
        assert!(consumer.has_gemm_row(&key));

        let original_identity = sealed.receipt_identity();
        let mut field_tampers = Vec::new();
        let mut tampered = sealed.clone();
        tampered.kernel.push('x');
        field_tampers.push(tampered);
        let mut tampered = sealed.clone();
        tampered.shape_class.push('x');
        field_tampers.push(tampered);
        let mut tampered = sealed.clone();
        tampered.machine[0] ^= 1;
        field_tampers.push(tampered);
        let mut tampered = sealed.clone();
        tampered.params.push(' ');
        field_tampers.push(tampered);
        let mut tampered = sealed.clone();
        tampered.measured.push(' ');
        field_tampers.push(tampered);
        let mut tampered = sealed.clone();
        tampered.memory_limit_bytes ^= 1;
        field_tampers.push(tampered);
        let mut tampered = sealed.clone();
        tampered.probe_buffer_bytes ^= 1;
        field_tampers.push(tampered);
        assert!(
            field_tampers
                .iter()
                .all(|tampered| tampered.receipt_identity() != original_identity),
            "every ledger tuple field must participate in the derive-key identity"
        );

        let other_probe = gemm_tune_key(4, 320, 289, 300).expect("other key");
        let mut wrong_context = Tuner::cold(0xAA55);
        assert!(
            adopt_cached_row(&mut wrong_context, &other_probe, &params, &measured)
                .expect("wrong context is a cache miss")
                .is_none()
        );
        assert!(!wrong_context.has_gemm_row(&other_probe));
    }

    #[test]
    fn cancelled_dispatch_preserves_progress_but_records_no_success_decision() {
        let key = gemm_tune_key(4, 320, 288, 300).expect("key");
        let mut tuner = Tuner::cold(0xAA55);
        tuner
            .pin_gemm_blocking(&key, GemmBlockPlan::COLD_START)
            .expect("pin");
        let decision = tuner.prepare_gemm_decision(&key);
        let error = execute_prepared_decision(&mut tuner, decision, |_| {
            Err::<(), _>(GemmTuneError::from(fs_la::GemmCancelled {
                report: fs_la::GemmRunReport {
                    declared_run: fs_exec::RunId(9),
                    completed_tiles: 7,
                    total_tiles: 19,
                    pool_runs: Vec::new(),
                    memory: fs_la::GemmMemoryReport {
                        limit_bytes: 1_024,
                        requested_bytes: 512,
                        peak_used_bytes: 384,
                        ..fs_la::GemmMemoryReport::default()
                    },
                },
            }))
        })
        .expect_err("cancelled producer");
        let GemmTuneError::Cancelled {
            limit_bytes,
            peak_used_bytes,
            report: Some(report),
        } = error
        else {
            panic!("expected retained cancellation report");
        };
        assert_eq!(limit_bytes, 1_024);
        assert_eq!(peak_used_bytes, 384);
        assert_eq!(report.completed_tiles, 7);
        assert_eq!(report.total_tiles, 19);
        assert_eq!(report.memory.requested_bytes, 512);
        assert_eq!(report.memory.peak_used_bytes, 384);
        assert!(tuner.decisions().is_empty());
        assert!(tuner.has_gemm_pin(&key));
    }

    #[test]
    fn executor_failure_retains_full_memory_and_progress_report() {
        let report = fs_la::GemmRunReport {
            declared_run: fs_exec::RunId(12),
            completed_tiles: 3,
            total_tiles: 11,
            pool_runs: Vec::new(),
            memory: fs_la::GemmMemoryReport {
                limit_bytes: 2_048,
                requested_bytes: 1_024,
                peak_used_bytes: 768,
                ..fs_la::GemmMemoryReport::default()
            },
        };
        let error = GemmTuneError::from(fs_la::GemmRunError::Executor {
            error: fs_exec::RunError::Incomplete {
                kernel: "fixture",
                tile: 4,
            },
            report,
        });
        let GemmTuneError::Executor {
            limit_bytes,
            peak_used_bytes,
            report,
            ..
        } = error
        else {
            panic!("expected retained executor report");
        };
        assert_eq!(limit_bytes, 2_048);
        assert_eq!(peak_used_bytes, 768);
        assert_eq!(report.completed_tiles, 3);
        assert_eq!(report.total_tiles, 11);
        assert_eq!(report.memory.requested_bytes, 1_024);
        assert_eq!(report.memory.peak_used_bytes, 768);
    }

    #[test]
    fn error_reports_are_indirected_off_the_success_path() {
        assert!(
            core::mem::size_of::<GemmTuneError>() <= 128,
            "GemmTuneError grew to {} bytes; keep full reports behind indirection",
            core::mem::size_of::<GemmTuneError>()
        );
    }
}
