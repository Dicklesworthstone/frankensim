//! The production GEMM autotune loop (bead yqug): measure → cache →
//! model → dispatch, closed end-to-end.
//!
//! [`gemm_f64_session`] is the production consumer the tuner was built
//! for: it resolves an MC/NC [`GemmBlockPlan`] for the caller's shape
//! class (pins beat cached rows beat the documented cold-start default),
//! runs a BOUNDED candidate sweep when the machine is cold, records the
//! ranked wall-time evidence as a tune row, writes it through to the
//! ledger `tune` table, and dispatches `fs_la::gemm_f64_parallel_with`
//! with the selected plan.
//!
//! Honesty boundaries, in the fs-exec tuner's division:
//! - The KERNEL KEY embeds fs-la's `GEMM_BIT_SEMANTICS_VERSION`, so rows
//!   measured under a different accumulation contract can never match a
//!   lookup (semantic filtering by construction). Rows are additionally
//!   machine-fingerprint-keyed; the ledger read path refuses stale
//!   (other-machine) and non-canonical rows instead of adopting them.
//! - MC/NC are BIT-NEUTRAL by fs-la's determinism contract, and the
//!   sweep ENFORCES that: every candidate's output must be bitwise
//!   identical to the first candidate's, else the loop fails closed
//!   with [`GemmTuneError::BitDrift`] and records nothing. KC and the
//!   SIMD tier are part of the bit contract and are NOT in this loop.
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
    CancelGate, GEMM_KERNEL_PREFIX, GemmBlockPlan, TuneError, TuneEvidence, TuneObservation,
    TuneSource, Tuner,
};
use fs_ledger::Ledger;

/// The bounded sweep lattice: 4 × 2 candidates, lattice order (mc-major
/// ascending). Chosen around the measured xlvx s5 landscape: thin bands
/// won both reference machines; the extremes document the neighborhood.
const SWEEP_MC: [usize; 4] = [16, 32, 64, 128];
const SWEEP_NC_CAP: [usize; 2] = [512, 2048];

/// Probe dims are capped so a cold-start sweep stays bounded (seconds,
/// not minutes) even when the caller's problem is huge.
const PROBE_DIM_CAP: usize = 512;

/// Wall-time samples per candidate (min-of ranking, all survive in the
/// evidence row).
const SWEEP_SAMPLES: usize = 3;

/// A structured autotune-loop failure. Every variant fails closed: no
/// tune row is recorded and nothing is dispatched with unvalidated
/// blocking.
#[derive(Debug)]
pub enum GemmTuneError {
    /// The cancel gate was requested during the sweep.
    Cancelled,
    /// Tuner-side refusal (invalid pin, evidence, or adoption).
    Tune(TuneError),
    /// Ledger cache I/O failed (the loop does not guess around storage).
    Ledger(String),
    /// Two sweep candidates produced different output bits: the
    /// bit-neutrality contract is broken and NO plan may be selected.
    BitDrift {
        /// Canonical params of the candidate that diverged.
        candidate: String,
    },
}

impl core::fmt::Display for GemmTuneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Cancelled => write!(f, "gemm autotune sweep cancelled"),
            Self::Tune(e) => write!(f, "gemm autotune: {e}"),
            Self::Ledger(detail) => write!(f, "gemm autotune ledger cache: {detail}"),
            Self::BitDrift { candidate } => write!(
                f,
                "gemm autotune: candidate {candidate} broke the MC/NC bit-neutrality contract"
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

/// The receipt for one autotuned dispatch: what ran, under which plan,
/// and where the plan came from. A study records this; replay pins it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GemmDispatch {
    /// Kernel key (embeds fs-la's GEMM bit-semantics version).
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

/// The shape class for an (m, n, k) problem: power-of-two buckets, so
/// nearby problems share rows and the class count stays bounded.
#[must_use]
pub fn gemm_shape_class(m: usize, n: usize, k: usize) -> String {
    format!("m{}-n{}-k{}", bucket(m), bucket(n), bucket(k))
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

/// FNV-1a over the output bits: the sweep's bit-neutrality witness.
fn bits_hash(c: &[f64]) -> u64 {
    let mut h = 0xcbf2_9ce4_8422_2325_u64;
    for value in c {
        for byte in value.to_bits().to_le_bytes() {
            h ^= u64::from(byte);
            h = h.wrapping_mul(0x100_0000_01b3);
        }
    }
    h
}

/// Run the bounded candidate sweep for one shape class and record the
/// winning row (ranked wall-time evidence). Returns the winning plan.
fn run_sweep(
    tuner: &mut Tuner,
    kernel: &str,
    shape_class: &str,
    gate: &CancelGate,
    threads: usize,
    m: usize,
    n: usize,
    k: usize,
) -> Result<GemmBlockPlan, GemmTuneError> {
    // Probe at the CALLER's dims (capped): the oracle lane showed that
    // probing at the class's power-of-two bucket flips winners — at
    // m = 320 the band count under each mc differs from m = 512, and
    // band balance decides the ranking. The row is still keyed by the
    // shared shape class; its evidence records the probe that measured
    // it (first-measurer-wins within a class, honestly labeled).
    let pm = m.clamp(1, PROBE_DIM_CAP);
    let pn = n.clamp(1, PROBE_DIM_CAP);
    let pk = k.clamp(1, PROBE_DIM_CAP);
    let mut a = vec![0.0f64; pm * pk];
    let mut b = vec![0.0f64; pk * pn];
    probe_fill(&mut a, 0xA);
    probe_fill(&mut b, 0xB);
    let mut c = vec![0.0f64; pm * pn];

    let mut observations = Vec::with_capacity(SWEEP_MC.len() * SWEEP_NC_CAP.len());
    let mut ranked: Vec<(u64, usize, GemmBlockPlan)> = Vec::new();
    let mut reference_bits: Option<u64> = None;
    for (index, (mc, nc_cap)) in SWEEP_MC
        .iter()
        .flat_map(|&mc| SWEEP_NC_CAP.iter().map(move |&nc| (mc, nc)))
        .enumerate()
    {
        if gate.is_requested() {
            return Err(GemmTuneError::Cancelled);
        }
        let plan = GemmBlockPlan::new(mc, nc_cap)?;
        let mut samples_ns = Vec::with_capacity(SWEEP_SAMPLES);
        for _ in 0..SWEEP_SAMPLES {
            c.fill(0.0);
            let t0 = std::time::Instant::now();
            fs_la::gemm_f64_parallel_with(
                pm,
                pn,
                pk,
                1.0,
                &a,
                &b,
                0.0,
                &mut c,
                threads,
                plan.mc,
                pn.min(plan.nc_cap).max(1),
            );
            let ns = u64::try_from(t0.elapsed().as_nanos()).unwrap_or(u64::MAX);
            samples_ns.push(ns.max(1));
        }
        // Bit-neutrality enforcement: every candidate must reproduce the
        // first candidate's bits exactly, or the loop refuses to select.
        let bits = bits_hash(&c);
        match reference_bits {
            None => reference_bits = Some(bits),
            Some(expected) if bits != expected => {
                return Err(GemmTuneError::BitDrift {
                    candidate: plan.canonical(),
                });
            }
            Some(_) => {}
        }
        let best = samples_ns.iter().copied().min().unwrap_or(u64::MAX);
        ranked.push((best, index, plan));
        observations.push(TuneObservation::wall_time(plan.canonical(), samples_ns)?);
    }
    // The declared selection rule: argmin of minimum wall time, ties to
    // the earlier lattice index.
    ranked.sort_unstable_by_key(|&(ns, index, _)| (ns, index));
    let winner = ranked[0].2;
    let evidence = TuneEvidence::ranked_wall_times(observations)?;
    tuner.record_gemm_row(kernel, shape_class, winner, evidence)?;
    Ok(winner)
}

/// The production autotuned f64 GEMM: `c = alpha·a·b + beta·c` through
/// the measure → cache → model → dispatch loop.
///
/// Resolution order: a pinned plan dispatches immediately (replay
/// fidelity — no measurement on a pinned path); else a cached row for
/// this kernel key × shape class (in the tuner, seeded from `ledger`'s
/// `tune` table when supplied); else the bounded sweep measures, records
/// evidence, writes through to the ledger, and the tuned row dispatches.
///
/// # Errors
/// [`GemmTuneError`] — cancelled sweeps, tuner refusals, ledger I/O, or
/// a bit-neutrality violation. On error `c` holds either its original
/// contents (pre-sweep failures) or the correct GEMM result is never
/// partially written: dispatch happens only after a plan is selected.
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
    let kernel = gemm_kernel_key();
    let shape_class = gemm_shape_class(m, n, k);
    let mut swept = false;

    if !tuner.has_pin(&kernel) && !tuner.has_row(&kernel, &shape_class) {
        // Cache tier: try the ledger before measuring. Stale
        // (other-machine) or non-canonical rows are refused by
        // adopt_row_json and we fall through to a fresh sweep.
        let mut seeded = false;
        if let Some(ledger) = ledger {
            let cached = ledger
                .tune_get(&kernel, &shape_class, &tuner.machine().to_le_bytes())
                .map_err(|e| GemmTuneError::Ledger(e.to_string()))?;
            if let Some(row) = cached {
                seeded = tuner.adopt_row_json(&row.measured).is_ok();
            }
        }
        if !seeded {
            run_sweep(tuner, &kernel, &shape_class, gate, threads, m, n, k)?;
            swept = true;
            if let Some(ledger) = ledger {
                let line = tuner
                    .row_json(&kernel, &shape_class)
                    .expect("the sweep just recorded this row");
                ledger
                    .tune_put(
                        &kernel,
                        &shape_class,
                        &tuner.machine().to_le_bytes(),
                        "\"gemm-block-plan\"",
                        &line,
                    )
                    .map_err(|e| GemmTuneError::Ledger(e.to_string()))?;
            }
        }
    }

    let (plan, source) = tuner.gemm_blocking_for(&kernel, &shape_class);
    fs_la::gemm_f64_parallel_with(
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
    );
    Ok(GemmDispatch {
        kernel,
        shape_class,
        plan,
        source,
        swept,
    })
}
