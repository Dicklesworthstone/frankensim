//! The production GEMM autotune loop (bead yqug): measure → cache →
//! model → dispatch, closed end-to-end.
//!
//! [`gemm_f64_session`] is the production consumer the tuner was built
//! for: it resolves an MC/NC [`GemmBlockPlan`] for the caller's shape
//! class (pins beat cached rows beat the documented cold-start default),
//! runs a BOUNDED candidate sweep when the machine is cold, records the
//! ranked wall-time evidence as a tune row, writes it through to the
//! ledger `tune` table, and dispatches
//! `fs_la::gemm_f64_parallel_with_cancel` with the selected plan.
//!
//! Honesty boundaries, in the fs-exec tuner's division:
//! - The KERNEL KEY embeds fs-la's `GEMM_BIT_SEMANTICS_VERSION`, so rows
//!   measured under a different accumulation contract can never match a
//!   lookup (semantic filtering by construction). Rows are additionally
//!   bound to the exact probe dims, requested/normalized thread budget,
//!   resolved ISA tier, placement policy, and implementation version, then
//!   machine-fingerprint-keyed. The ledger read path refuses stale,
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
    PreparedGemmDecision, PreparedGemmRow, TuneError, TuneEvidence, TuneObservation, TuneSource,
    Tuner,
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

/// A structured autotune-loop failure. Every variant fails closed: sweep
/// failures record no row and nothing dispatches under unvalidated blocking.
/// A cancellation during the final dispatch may retain the already validated
/// measured row, but records no successful decision and does not commit `C`.
#[derive(Debug)]
pub enum GemmTuneError {
    /// The cancel gate was requested. Compute may have completed in private
    /// staging, but the caller's output was not committed.
    Cancelled {
        /// Completed bounded compute tiles before the request was observed.
        completed_tiles: usize,
        /// Total bounded compute tiles in the interrupted dispatch.
        total_tiles: usize,
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
}

impl core::fmt::Display for GemmTuneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Cancelled {
                completed_tiles,
                total_tiles,
            } => write!(
                f,
                "gemm work cancelled after {completed_tiles}/{total_tiles} compute tiles; output not committed"
            ),
            Self::Tune(e) => write!(f, "gemm autotune: {e}"),
            Self::Ledger(detail) => write!(f, "gemm autotune ledger cache: {detail}"),
            Self::BitDrift { candidate, repeat } => write!(
                f,
                "gemm autotune: candidate {candidate} repeat {repeat} broke the MC/NC bit-neutrality contract"
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
        Self::Cancelled {
            completed_tiles: cancelled.report.completed_tiles,
            total_tiles: cancelled.report.total_tiles,
        }
    }
}

fn cancelled_before_compute() -> GemmTuneError {
    GemmTuneError::Cancelled {
        completed_tiles: 0,
        total_tiles: 0,
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
    let implementation = format!(
        "fs-la-{}-gemm-v{}",
        fs_la::VERSION,
        fs_la::GEMM_IMPLEMENTATION_VERSION
    );
    let execution = GemmExecutionIdentity::new(
        threads,
        threads.max(1),
        probe_dims(m, n, k),
        fs_la::gemm_execution_tier(),
        fs_la::GEMM_PARALLEL_IMPLEMENTATION,
        implementation,
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
    mut run: R,
) -> Result<SweepResult, GemmTuneError>
where
    R: FnMut(&SweepCandidate, &mut [f64]) -> Result<(), GemmTuneError>,
{
    let mut c = vec![0.0f64; output_len];
    let mut observations = Vec::with_capacity(candidates.len());
    let mut ranked: Vec<(u64, usize, GemmBlockPlan)> = Vec::with_capacity(candidates.len());
    let mut reference_bits: Option<Vec<u64>> = None;
    for (index, candidate) in candidates.iter().enumerate() {
        if gate.is_requested() {
            return Err(cancelled_before_compute());
        }
        let mut samples_ns = Vec::with_capacity(SWEEP_SAMPLES);
        for repeat in 1..=SWEEP_SAMPLES {
            if gate.is_requested() {
                return Err(cancelled_before_compute());
            }
            c.fill(0.0);
            let t0 = std::time::Instant::now();
            run(candidate, &mut c)?;
            let ns = u64::try_from(t0.elapsed().as_nanos()).unwrap_or(u64::MAX);
            samples_ns.push(ns.max(1));
            if gate.is_requested() {
                return Err(cancelled_before_compute());
            }

            // Compare every output word directly. A fixed-width digest is
            // not a proof of bit-neutrality and would also hide which repeat
            // drifted. `to_bits` intentionally distinguishes signed zero and
            // every NaN payload.
            let bits: Vec<u64> = c.iter().map(|value| value.to_bits()).collect();
            match &reference_bits {
                None => reference_bits = Some(bits),
                Some(expected) if bits != *expected => {
                    return Err(GemmTuneError::BitDrift {
                        candidate: candidate.plan.canonical(),
                        repeat,
                    });
                }
                Some(_) => {}
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
    threads: usize,
    m: usize,
    n: usize,
    k: usize,
) -> Result<SweepResult, GemmTuneError> {
    // Probe at the CALLER's dims (capped): the oracle lane showed that
    // probing at the class's power-of-two bucket flips winners — at
    // m = 320 the band count under each mc differs from m = 512, and
    // band balance decides the ranking. The row retains the bucketed shape
    // class, but the exact capped probe is also part of the scoped key so a
    // neighboring caller cannot silently inherit different evidence.
    let [pm, pn, pk] = probe_dims(m, n, k);
    let mut a = vec![0.0f64; pm * pk];
    let mut b = vec![0.0f64; pk * pn];
    probe_fill(&mut a, 0xA);
    probe_fill(&mut b, 0xB);
    let candidates = effective_sweep_candidates(pm, pn)?;
    measure_candidates(gate, &candidates, pm * pn, |candidate, c| {
        fs_la::gemm_f64_parallel_with_cancel(
            pm,
            pn,
            pk,
            1.0,
            &a,
            &b,
            0.0,
            c,
            threads,
            candidate.effective_mc,
            candidate.effective_nc,
            gate,
        )
        .map(|_| ())
        .map_err(GemmTuneError::from)
    })
}

/// Persist a validated measured row before installing it in the process-local
/// tuner. `persist` is injectable so the failure-atomic boundary is directly
/// testable without corrupting a real ledger connection.
fn install_sweep_row<P>(
    tuner: &mut Tuner,
    key: &GemmTuneKey,
    sweep: SweepResult,
    persist: P,
) -> Result<GemmBlockPlan, GemmTuneError>
where
    P: FnOnce(&PreparedGemmRow) -> Result<(), GemmTuneError>,
{
    let prepared = tuner.prepare_gemm_row(key, sweep.winner, sweep.evidence)?;
    persist(&prepared)?;
    let winner = sweep.winner;
    tuner.commit_gemm_row(prepared)?;
    Ok(winner)
}

fn adopt_cached_row(
    tuner: &mut Tuner,
    key: &GemmTuneKey,
    params: &str,
    measured: &str,
) -> Result<bool, GemmTuneError> {
    let Ok(prepared) = tuner.prepare_adopt_gemm_row_json(key, measured) else {
        return Ok(false);
    };
    if params != prepared.params_json() {
        return Ok(false);
    }
    tuner.commit_gemm_row(prepared)?;
    Ok(true)
}

fn execute_prepared_decision<R>(
    tuner: &mut Tuner,
    decision: PreparedGemmDecision,
    run: R,
) -> Result<(GemmBlockPlan, TuneSource), GemmTuneError>
where
    R: FnOnce(GemmBlockPlan) -> Result<(), GemmTuneError>,
{
    let plan = decision.plan();
    let source = decision.source();
    run(plan)?;
    // Exclusive access to `tuner` spans prepare -> run -> commit, so no
    // applicable pin/row can change and make this prepared decision stale.
    tuner
        .commit_gemm_decision(decision)
        .expect("exclusive tuner borrow preserves a prepared GEMM decision");
    Ok((plan, source))
}

/// The production autotuned f64 GEMM: `c = alpha·a·b + beta·c` through
/// the measure → cache → model → dispatch loop.
///
/// Resolution order after shape and cancellation preflight: a pinned plan
/// dispatches without measurement; else an exact cached row (in the tuner,
/// seeded from `ledger` when supplied); else the bounded sweep measures,
/// persists a prepared row, commits it locally, and dispatches. Serial,
/// small-M, and no-product calls bypass tuning entirely.
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
    ledger: Option<&Ledger>,
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
    // Public slice/extent preconditions are checked before tier resolution,
    // cache reads, sweeps, rows, or decisions. Invalid work cannot poison the
    // tuning state and `c` is still untouched when this panics.
    assert_contiguous_shapes(m, n, k, a, b, c);
    if gate.is_requested() {
        return Err(cancelled_before_compute());
    }

    let key = gemm_tune_key(threads, m, n, k)?;
    let kernel = key.kernel().to_string();
    let shape_class = gemm_shape_class(m, n, k);
    let mut swept = false;

    // No product, one-thread, and small-M routes do not have a meaningful
    // production MC/NC choice. Dispatch them cancellation-correctly under the
    // documented cold plan without reading or mutating tune state.
    if !fs_la::gemm_tuning_is_effective(m, n, k, alpha, threads) {
        let plan = GemmBlockPlan::COLD_START;
        fs_la::gemm_f64_parallel_with_cancel(
            m,
            n,
            k,
            alpha,
            a,
            b,
            beta,
            c,
            threads,
            plan.mc,
            n.min(plan.nc_cap).max(1),
            gate,
        )
        .map_err(GemmTuneError::from)?;
        return Ok(GemmDispatch {
            kernel,
            shape_class,
            plan,
            source: TuneSource::ColdStart,
            swept,
        });
    }

    if !tuner.has_gemm_pin(&key) && !tuner.has_gemm_row(&key) {
        // Cache tier: try the ledger before measuring. Stale
        // (other-machine) or non-canonical rows are refused by
        // prepare_adopt_gemm_row_json and we fall through to a fresh sweep.
        // The ledger's separate params column must agree byte-for-byte with
        // the validated row body before either is allowed into the tuner.
        let mut seeded = false;
        if let Some(ledger) = ledger {
            let cached = ledger
                .tune_get(
                    key.kernel(),
                    key.shape_class(),
                    &tuner.machine().to_le_bytes(),
                )
                .map_err(|e| GemmTuneError::Ledger(e.to_string()))?;
            if let Some(row) = cached {
                seeded = adopt_cached_row(tuner, &key, &row.params, &row.measured)?;
            }
        }
        if !seeded {
            let sweep = run_sweep(gate, threads, m, n, k)?;
            swept = true;
            let machine = tuner.machine().to_le_bytes();
            install_sweep_row(tuner, &key, sweep, |prepared| {
                let Some(ledger) = ledger else {
                    return Ok(());
                };
                ledger
                    .tune_put(
                        key.kernel(),
                        key.shape_class(),
                        &machine,
                        &prepared.params_json(),
                        &prepared.row_json(),
                    )
                    .map_err(|e| GemmTuneError::Ledger(e.to_string()))
            })?;
        }
    }

    if gate.is_requested() {
        return Err(cancelled_before_compute());
    }
    let decision = tuner.prepare_gemm_decision(&key);
    let (plan, source) = execute_prepared_decision(tuner, decision, |plan| {
        fs_la::gemm_f64_parallel_with_cancel(
            m,
            n,
            k,
            alpha,
            a,
            b,
            beta,
            c,
            threads,
            plan.mc,
            n.min(plan.nc_cap).max(1),
            gate,
        )
        .map(|_| ())
        .map_err(GemmTuneError::from)
    })?;
    Ok(GemmDispatch {
        kernel,
        shape_class,
        plan,
        source,
        swept,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

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
            let error = measure_candidates(&CancelGate::new(), &candidates, 2, |_, c| {
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
                Ok(())
            })
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
        measure_candidates(&CancelGate::new(), &wide, 1, |candidate, c| {
            executed.push((candidate.effective_mc, candidate.effective_nc));
            c[0] = 1.0;
            Ok(())
        })
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
        let error = measure_candidates(&gate, &candidates, 1, |_, c| {
            c[0] = 1.0;
            gate.request();
            Ok(())
        })
        .expect_err("the post-repeat poll must observe cancellation");
        assert!(matches!(
            error,
            GemmTuneError::Cancelled {
                completed_tiles: 0,
                total_tiles: 0
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

        let winner = install_sweep_row(&mut tuner, &key, synthetic_sweep(), |_| Ok(()))
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
        install_sweep_row(&mut producer, &key, synthetic_sweep(), |prepared| {
            params = prepared.params_json();
            measured = prepared.row_json();
            Ok(())
        })
        .expect("produce cached row");

        let mut consumer = Tuner::cold(0xAA55);
        assert!(
            !adopt_cached_row(&mut consumer, &key, "\"mc=32,nc-cap=512\"", &measured)
                .expect("mismatch is a cache miss")
        );
        assert!(!consumer.has_gemm_row(&key));
        assert!(adopt_cached_row(&mut consumer, &key, &params, &measured).expect("adopt"));
        assert!(consumer.has_gemm_row(&key));

        let other_probe = gemm_tune_key(4, 320, 289, 300).expect("other key");
        let mut wrong_context = Tuner::cold(0xAA55);
        assert!(
            !adopt_cached_row(&mut wrong_context, &other_probe, &params, &measured)
                .expect("wrong context is a cache miss")
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
            Err(GemmTuneError::from(fs_la::GemmCancelled {
                report: fs_la::GemmRunReport {
                    completed_tiles: 7,
                    total_tiles: 19,
                },
            }))
        })
        .expect_err("cancelled producer");
        assert!(matches!(
            error,
            GemmTuneError::Cancelled {
                completed_tiles: 7,
                total_tiles: 19
            }
        ));
        assert!(tuner.decisions().is_empty());
        assert!(tuner.has_gemm_pin(&key));
    }
}
