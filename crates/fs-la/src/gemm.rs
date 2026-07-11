//! BLIS-style GEMM (plan §6.1): C = α·A·B + β·C, row-major, with cache
//! blocking and panel packing. v1 is the CORRECTNESS + ARCHITECTURE layer
//! in safe Rust: the register-tiled microkernel accumulates in local
//! arrays with fused `mul_add` (auto-vectorizes respectably); the
//! arch-specific fs-simd capsule microkernels and the autotuned blocking
//! sweep are the recorded perf follow-up (gated on the autotuner bead).
//!
//! DETERMINISM CONTRACT: accumulation order is k-ascending within KC
//! chunks, with per-chunk register partials folded into C in chunk order.
//! Hence KC is PART of the bit contract (changing it legitimately changes
//! bits → golden bump with justification); MC/NC are bit-neutral (pure
//! m/n tiling — the fact the future parallel lane relies on). Everything
//! is fixed-order +/×/mul_add: cross-ISA bit-deterministic by
//! construction, golden-hashed in tests.

/// GEMM BIT-SEMANTICS VERSION (bead y4pt): bump on ANY change to the
/// accumulation order contract — KC retuning, k-order changes, or a
/// microkernel that stops being bitwise-equal to the scalar twin.
/// Downstream goldens pin this in golden-couplings.json;
/// `cargo run -p xtask -- check-goldens` fails on drift until they are
/// deliberately re-frozen.
pub const GEMM_BIT_SEMANTICS_VERSION: u32 = 1;

/// Micro-tile rows (A panel height). Pre-autotuner default.
const MR: usize = 8;
/// Micro-tile cols (B panel width). Pre-autotuner default.
const NR: usize = 4;
/// K blocking — PART OF THE BIT CONTRACT (see module docs).
const KC: usize = 256;
/// M blocking (bit-neutral).
const MC: usize = 128;
/// N blocking (bit-neutral).
const NC: usize = 512;

/// Version of the production f64 GEMM implementation and scheduling
/// surface. This is separate from [`GEMM_BIT_SEMANTICS_VERSION`]: a
/// bit-neutral scheduling or cancellation change bumps this identity while
/// leaving the numerical contract alone.
pub const GEMM_IMPLEMENTATION_VERSION: u32 = 4;

/// Domain used to derive each NC/KC panel's child run identity from the
/// caller-ledgered GEMM operation run.
pub const GEMM_PANEL_RUN_DOMAIN: &str = "org.frankensim.fs-la.gemm-panel-run.v1";

/// Deterministic child run for one NC/KC panel of a declared GEMM operation.
#[must_use]
pub fn gemm_panel_run_id(operation: fs_exec::RunId, panel_ordinal: u64) -> fs_exec::RunId {
    operation.derive(GEMM_PANEL_RUN_DOMAIN, panel_ordinal)
}

/// BLAKE3 fingerprint of the compiler, Cargo profile/codegen inputs, target,
/// explicit Rust flags, workspace manifests, the bounded GEMM execution source
/// closure, and optional operator-supplied `FRANKENSIM_GEMM_CODEGEN_ID` salt
/// for this build.
///
/// This is performance identity rather than numerical identity: two binaries
/// can preserve [`GEMM_BIT_SEMANTICS_VERSION`] while requiring independent
/// tune rows because their generated code differs.
pub const GEMM_BUILD_FINGERPRINT: &str = env!("FS_LA_GEMM_BUILD_FINGERPRINT");

/// Maximum arithmetic work in one cancellable GEMM compute quantum. Packing
/// and beta staging use smaller fixed quanta; the largest poll interval is one
/// `MR x NR x KC` microtile plus its alpha/write-back FMA per output.
pub const GEMM_MAX_FMAS_BETWEEN_POLLS: usize = MR * NR * (KC + 1);

/// Number of output elements staged between cancellation polls.
const C_STAGE_TILE_ELEMENTS: usize = 4096;

/// Explicit operation memory envelope for the pool GEMM path (bead
/// wf9.15): the byte ceiling the caller grants root orchestration.
/// Preflight-checked BEFORE any allocation or C mutation; the default
/// is unbounded (existing callers unchanged).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GemmMemoryEnvelope {
    /// Ceiling over staging + shared pack + band metadata + per-worker
    /// arena panels, in bytes.
    pub limit_bytes: u64,
}

impl GemmMemoryEnvelope {
    /// No caller-imposed ceiling.
    pub const UNBOUNDED: GemmMemoryEnvelope = GemmMemoryEnvelope {
        limit_bytes: u64::MAX,
    };
}

impl Default for GemmMemoryEnvelope {
    fn default() -> Self {
        Self::UNBOUNDED
    }
}

/// Deterministic logical-memory accounting for one pool GEMM (wf9.15).
///
/// Component and requested bytes describe the checked admission plan. Peak-used
/// bytes are the largest fs-la-owned logical reservation concurrency entered:
/// root buffers successfully reserved plus tile arena reservations whose
/// allocation attempt had begun. Counting the attempt before the allocator call
/// makes the high-water mark conservative even when one attempt refuses. This
/// is not process RSS and deliberately excludes generic TilePool worker, deque,
/// stack, and receipt internals (tracked separately by wf9.16).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub struct GemmMemoryReport {
    /// The envelope in force (u64::MAX = unbounded).
    pub limit_bytes: u64,
    /// Transactional C staging bytes.
    pub staging_bytes: u128,
    /// Shared packed-B panel bytes.
    pub b_pack_bytes: u128,
    /// Reusable band-metadata bytes (mutex-guarded band slots).
    pub band_metadata_bytes: u128,
    /// Capacity reserved for fs-la's ordered panel-run receipt vector.
    pub pool_run_bytes: u128,
    /// Fresh-arena reservation for one A micro-panel.
    pub arena_bytes_per_worker: u64,
    /// Maximum number of workers that can be active for one M-band traversal.
    pub active_arena_workers: usize,
    /// Per-active-worker A micro-panel arena reservations.
    pub arena_bytes: u128,
    /// Preflight total (checked sum of the above).
    pub requested_bytes: u128,
    /// Largest fs-la-owned logical live set reached by this attempt.
    pub peak_used_bytes: u128,
    /// Reservation bytes rejected at the failure boundary, or zero on success
    /// and cancellation.
    pub refused_bytes: u128,
}

/// Structured progress for a cancellation-aware GEMM dispatch. A successful
/// return always has `completed_tiles == total_tiles`; a cancelled return may
/// contain completed work, but that work exists only in private staging.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GemmRunReport {
    /// Caller-ledgered identity of the complete GEMM operation.
    pub declared_run: fs_exec::RunId,
    /// Fully written `MR x NR x KC` compute tiles.
    pub completed_tiles: usize,
    /// Total compute tiles required by this dispatch (zero for beta-only work).
    pub total_tiles: usize,
    /// One real fs-exec traversal receipt per dispatched NC/KC panel.
    /// Empty only when the operation has no arithmetic product or cancellation
    /// happened before the first M-band dispatch.
    pub pool_runs: Vec<fs_exec::RunReport>,
    /// The memory plan the operation was admitted under (wf9.15):
    /// requested/limit bytes for staging, pack, metadata, and arenas.
    pub memory: GemmMemoryReport,
}

/// Cancellation observed at a bounded GEMM poll point after all scoped
/// workers drained. The caller's output remains bitwise unchanged.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GemmCancelled {
    /// Work completed in private staging before the request was observed.
    pub report: Box<GemmRunReport>,
}

impl core::fmt::Display for GemmCancelled {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "gemm cancelled after {}/{} compute tiles; output not committed",
            self.report.completed_tiles, self.report.total_tiles
        )
    }
}

impl core::error::Error for GemmCancelled {}

/// Failure from the caller-owned TilePool GEMM path.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GemmRunError {
    /// The gate was observed at a bounded poll point after every worker and
    /// arena scope drained. Caller-visible C is unchanged.
    Cancelled(GemmCancelled),
    /// Tile panic or executor invariant refusal, with tile provenance from
    /// fs-exec. Caller-visible C is unchanged.
    Executor {
        /// Structured fs-exec failure with logical tile provenance.
        error: fs_exec::RunError,
        /// Compute and traversal progress accumulated before refusal.
        report: Box<GemmRunReport>,
    },
    /// Refused at the memory boundary (wf9.15): the preflight plan exceeded
    /// the caller's envelope, or an fs-la-owned reservation was declined. Any
    /// dispatched panel was drained by the pool's scoped join, and
    /// caller-visible C is bitwise unchanged.
    MemoryRefused {
        /// Which reservation was refused.
        what: &'static str,
        /// Bytes the refused reservation asked for.
        requested_bytes: u128,
        /// The envelope in force.
        limit_bytes: u64,
        /// Progress and logical-memory accounting at the drained boundary.
        report: Box<GemmRunReport>,
    },
    /// Arithmetic needed to represent the logical memory plan exceeded u128.
    /// No allocation, dispatch, or caller-visible mutation occurred.
    MemoryPlanOverflow {
        /// Component whose checked arithmetic overflowed.
        what: &'static str,
        /// The envelope in force.
        limit_bytes: u64,
    },
}

impl core::fmt::Display for GemmRunError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Cancelled(error) => error.fmt(f),
            Self::Executor { error, .. } => write!(f, "gemm executor: {error}"),
            Self::MemoryRefused {
                what,
                requested_bytes,
                limit_bytes,
                ..
            } => write!(
                f,
                "gemm memory refused at {what}: {requested_bytes} bytes requested against a \
                 {limit_bytes}-byte envelope; output not committed; raise the envelope or shrink \
                 the operation"
            ),
            Self::MemoryPlanOverflow { what, limit_bytes } => write!(
                f,
                "gemm memory-plan arithmetic overflowed at {what} under the \
                 {limit_bytes}-byte envelope; output not touched"
            ),
        }
    }
}

impl core::error::Error for GemmRunError {
    fn source(&self) -> Option<&(dyn core::error::Error + 'static)> {
        match self {
            Self::Cancelled(error) => Some(error),
            Self::Executor { error, .. } => Some(error),
            Self::MemoryRefused { .. } | Self::MemoryPlanOverflow { .. } => None,
        }
    }
}

/// The SIMD tier ACTUALLY selected for GEMM's 8x4 microkernel.
///
/// This is operation-specific rather than the process-wide maximum: an
/// AVX-512-capable host currently reports `avx2` because GEMM's audited x86
/// capsule uses AVX2/FMA. It remains truthful under Miri, where the table
/// deliberately routes to the scalar implementation.
#[must_use]
pub fn gemm_execution_tier() -> &'static str {
    fs_simd::mk8x4_f64_tier().name()
}

/// Exact codegen/build fingerprint carried by production GEMM tune keys.
#[must_use]
pub const fn gemm_build_identity() -> &'static str {
    GEMM_BUILD_FINGERPRINT
}

/// Whether MC/NC tuning can affect the legacy production parallel route.
/// Small and single-thread calls deliberately bypass tuning because thread
/// dispatch overhead dominates and [`gemm_f64_parallel`] routes them through
/// the serial kernel. No-op products likewise have no blocking decision to
/// measure.
#[must_use]
pub fn gemm_tuning_is_effective(m: usize, n: usize, k: usize, alpha: f64, threads: usize) -> bool {
    threads > 1 && m >= 2 * MC && n != 0 && k != 0 && alpha != 0.0
}

#[track_caller]
fn checked_product(label: &str, lhs: usize, rhs: usize) -> usize {
    lhs.checked_mul(rhs)
        .unwrap_or_else(|| panic!("{label} extent overflow: {lhs} * {rhs}"))
}

#[track_caller]
fn assert_contiguous_shapes<T, U>(m: usize, n: usize, k: usize, a: &[T], b: &[T], c: &[U]) {
    let a_len = checked_product("a", m, k);
    let b_len = checked_product("b", k, n);
    let c_len = checked_product("c", m, n);
    assert_eq!(a.len(), a_len, "a must be m*k = {a_len}");
    assert_eq!(b.len(), b_len, "b must be k*n = {b_len}");
    assert_eq!(c.len(), c_len, "c must be m*n = {c_len}");
}

#[track_caller]
fn assert_view_shape(name: &str, len: usize, rows: usize, cols: usize, ld: usize) {
    assert!(ld >= cols.max(1), "{name}: ld {ld} < view cols {cols}");
    if rows == 0 {
        return;
    }
    let row_offset = (rows - 1)
        .checked_mul(ld)
        .unwrap_or_else(|| panic!("{name}: row-stride extent overflow"));
    let need = row_offset
        .checked_add(cols)
        .unwrap_or_else(|| panic!("{name}: view extent overflow"));
    assert!(len >= need, "{name}: slice len {len} < view need {need}");
}

#[track_caller]
fn checked_round_up(label: &str, value: usize, quantum: usize) -> usize {
    debug_assert!(quantum > 0);
    let remainder = value % quantum;
    if remainder == 0 {
        value
    } else {
        value
            .checked_add(quantum - remainder)
            .unwrap_or_else(|| panic!("{label} extent overflow"))
    }
}

/// f64 GEMM: `c[m×n] = alpha · a[m×k] · b[k×n] + beta · c`, row-major
/// contiguous slices. β = 0 OVERWRITES c (existing NaN/garbage in c is
/// ignored — the BLAS convention callers expect for uninitialized output).
///
/// # Panics
/// Structured panics on slice-length mismatches.
#[allow(clippy::too_many_arguments)] // BLAS-shape signature: m,n,k,alpha,a,b,beta,c
pub fn gemm_f64(
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    beta: f64,
    c: &mut [f64],
) {
    assert_contiguous_shapes(m, n, k, a, b, c);
    // β pass first (once, before any KC chunk): scale or overwrite.
    scale_c(c, beta);
    if m == 0 || n == 0 || alpha == 0.0 {
        return;
    }
    if k == 0 {
        return; // C = beta*C only (already applied).
    }
    let mut a_pack = vec![0.0f64; MC * KC];
    let mut b_pack = vec![0.0f64; KC * NC];
    // Loop nest (BLIS order): NC → KC → MC → NR → MR → K.
    let mut jc = 0;
    while jc < n {
        let nc = NC.min(n - jc);
        let mut pc = 0;
        while pc < k {
            let kc = KC.min(k - pc);
            pack_b(&mut b_pack, b, n, pc, jc, kc, nc);
            let mut ic = 0;
            while ic < m {
                let mc = MC.min(m - ic);
                pack_a(&mut a_pack, a, k, ic, pc, mc, kc);
                macro_kernel(&a_pack, &b_pack, c, m, n, ic, jc, mc, nc, kc, alpha);
                ic += MC;
            }
            pc += KC;
        }
        jc += NC;
    }
    let _ = m; // (m used above; silences pedantic when MC >= m)
}

/// PARALLEL GEMM, shared-B design (bead xlvx item 3, v2): the packed
/// B panel for each (jc, pc) chunk is built ONCE and SHARED read-only
/// across threads, which then split the MC loop — each thread packs
/// its own A block and owns a disjoint contiguous C row band. The v1
/// row-band design (each thread running the whole loop nest) was
/// MEASURED to repack the entire B per thread and topped out at 0.107
/// of the all-core axis on a 64-thread Threadripper; sharing the pack
/// is the standard BLIS parallelization. BITWISE-FREE: the per-element
/// accumulation order (jc/pc chunk order, k order within) is exactly
/// the serial kernel's — gated across thread counts, no golden bump
/// (xdgf recorded fact (b)).
#[allow(clippy::too_many_arguments)]
pub fn gemm_f64_parallel(
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    beta: f64,
    c: &mut [f64],
    threads: usize,
) {
    let t = threads.max(1);
    gemm_f64_parallel_with(
        m,
        n,
        k,
        alpha,
        a,
        b,
        beta,
        c,
        t,
        MC_PAR,
        n.clamp(NR, NC_PAR_CAP),
    );
}

/// Parallel-path M blocking, MEASURED (xlvx s5 sweep at n = 2048 on a
/// 14t M4 Pro and an idle 128t 5995WX): thin mc = 32 bands won BOTH
/// machines (213 / 386 GFLOP/s vs 159 / 201 at the serial 128/512
/// defaults). The serial MC = 128 caps parallelism at m/128 bands
/// (16 threads matched 128 on the 5995WX before this); the opposite
/// extreme — bands ~= 3 per worker, mc = 8 on 128t — measured 94:
/// per-band pack/dispatch overhead swamps the extra workers.
const MC_PAR: usize = 32;
/// Parallel-path N blocking: nc = n (one A pass, one scope barrier per
/// KC chunk) dominated every mc row on both machines — the sweep was
/// monotone in nc. Capped so b_pack stays L3-resident (KC×2048×8 =
/// 4 MB) for huge n.
const NC_PAR_CAP: usize = 2048;

/// The tunable parallel engine behind [`gemm_f64_parallel`]: explicit
/// `mc_q` (band height) and `nc_q` (B-panel width) blocking. Both are
/// BIT-NEUTRAL (module docs): per-element accumulation stays jc/pc
/// chunk order with k ascending regardless of the m/n tiling — gated
/// in gemm_suite across an (mc, nc) grid. Public for the autotune
/// sweep lane; library callers want [`gemm_f64_parallel`].
#[doc(hidden)]
#[allow(clippy::too_many_arguments)]
pub fn gemm_f64_parallel_with(
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    beta: f64,
    c: &mut [f64],
    threads: usize,
    mc_q: usize,
    nc_q: usize,
) {
    assert_contiguous_shapes(m, n, k, a, b, c);
    let t = threads.max(1);
    // Values wider than the problem are equivalent to one logical block;
    // clamp them before sizing packs so an untrusted tune row cannot turn a
    // small multiplication into an unbounded allocation request.
    let mc_q = mc_q.max(MR).min(m.max(MR));
    let nc_q = nc_q.max(NR).min(n.max(NR));
    if t == 1 || m < 2 * MC {
        gemm_f64(m, n, k, alpha, a, b, beta, c);
        return;
    }
    if m == 0 || n == 0 || alpha == 0.0 || k == 0 {
        scale_c(c, beta);
        return;
    }
    let a_pack_rows = checked_round_up("parallel A pack rows", mc_q, MR);
    let b_pack_cols = checked_round_up("parallel B pack columns", nc_q, NR);
    let a_pack_len = checked_product("parallel A pack", a_pack_rows, KC);
    let b_pack_len = checked_product("parallel B pack", KC, b_pack_cols);
    let band_len = checked_product("parallel C band", mc_q, n);
    scale_c(c, beta);
    let mut b_pack = vec![0.0f64; b_pack_len];
    let mut jc = 0;
    while jc < n {
        let nc = nc_q.min(n - jc);
        let mut pc = 0;
        while pc < k {
            let kc = KC.min(k - pc);
            pack_b(&mut b_pack, b, n, pc, jc, kc, nc);
            let bp: &[f64] = &b_pack;
            // WORK-STEALING band dispenser (safe Rust, no capsule):
            // mc_q-row C bands behind a Mutex-guarded iterator; threads
            // pull the next band as they finish, so slow cores take
            // fewer (equal static shares let heterogeneous E-cores
            // drag the whole chunk — measured on M4 Pro). Bitwise
            // invariant: a band's content is a pure function of the
            // band, never of which thread computed it or in what
            // order; the lock guards ASSIGNMENT only.
            let dispenser = std::sync::Mutex::new(c.chunks_mut(band_len).enumerate());
            // Never spawn more workers than bands: excess threads only
            // lock, see None, and exit — 64 spawns for 4-16 bands
            // measured 2-9x slower than v2 on the 64-thread ts1.
            let workers = t.min(m.div_ceil(mc_q));
            std::thread::scope(|scope| {
                for _ in 0..workers {
                    let disp = &dispenser;
                    scope.spawn(move || {
                        let mut a_pack = vec![0.0f64; a_pack_len];
                        loop {
                            let next = disp.lock().expect("dispenser lock").next();
                            let Some((bi, band)) = next else { break };
                            let ic = bi * mc_q;
                            let mc = mc_q.min(m - ic);
                            pack_a(&mut a_pack, a, k, ic, pc, mc, kc);
                            // Band-local rows (offset 0); ld stays n.
                            macro_kernel(&a_pack, bp, band, m, n, 0, jc, mc, nc, kc, alpha);
                        }
                    });
                }
            });
            pc += KC;
        }
        jc += nc_q;
    }
}

/// Cancellation-aware tunable parallel GEMM.
///
/// This is the request -> drain -> finalize path for session/solver
/// orchestration. Computation happens in private staging. A cancellation
/// request stops new work, every scoped worker finishes at most its current
/// bounded packing panel or `MR x NR x KC` microtile, and the scope is joined
/// before [`GemmCancelled`] is returned. Therefore `Err` leaves `c` bitwise
/// unchanged. On success, the final gate poll is the FINALIZATION CUTOFF;
/// committing the completed staging buffer is non-cancellable and a request
/// arriving during that copy belongs to the caller's next operation.
///
/// Unlike [`gemm_f64_parallel_with`], this entry point does not route small or
/// single-thread problems through the non-cancellable serial facade. The
/// supplied MC/NC quanta remain effective so sweep and dispatch execute the
/// same kernel. Orchestrators should use [`gemm_tuning_is_effective`] to skip
/// meaningless wall-time sweeps while still dispatching here with the cold
/// plan.
///
/// # Errors
/// [`GemmRunError`] after observing `gate` or a contained executor failure and
/// draining every worker. The error reports private progress; `c` is unchanged.
///
/// # Panics
/// Structured panics on slice-length or extent mismatches, before `c` can be
/// mutated.
#[allow(clippy::too_many_arguments)]
pub fn gemm_f64_parallel_with_cancel(
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    beta: f64,
    c: &mut [f64],
    threads: usize,
    mc_q: usize,
    nc_q: usize,
    gate: &fs_exec::CancelGate,
) -> Result<GemmRunReport, GemmRunError> {
    let pool = fs_exec::TilePool::for_host(threads, 0x4653_2D4C_412D_474D);
    gemm_f64_parallel_with_pool(m, n, k, alpha, a, b, beta, c, &pool, mc_q, nc_q, gate)
}

/// Cancellation-aware GEMM on a caller-owned, reusable fs-exec pool.
///
/// Each NC/KC panel packs B once. Its disjoint MC row bands are then logical
/// [`fs_exec::TileKernel`] tiles scheduled by `pool`, so the pool's worker
/// budget, topology, quantum weights, pinning policy, `Cx` budget, stream
/// identity, and tile-scoped arenas are the execution path rather than tune
/// metadata detached from the kernel. The caller can reuse one pool across a
/// session; current fs-exec workers are still joined `std::thread` scopes per
/// run (see its contract no-claim).
///
/// # Errors
/// [`GemmRunError::Cancelled`] after a request -> drain -> finalize refusal, or
/// [`GemmRunError::Executor`] for a contained tile/executor failure. Both leave
/// `c` bitwise unchanged.
///
/// # Panics
/// Structured panics on slice-length or extent mismatches before `c` mutation.
#[allow(clippy::too_many_arguments)]
pub fn gemm_f64_parallel_with_pool(
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    beta: f64,
    c: &mut [f64],
    pool: &fs_exec::TilePool,
    mc_q: usize,
    nc_q: usize,
    gate: &fs_exec::CancelGate,
) -> Result<GemmRunReport, GemmRunError> {
    gemm_f64_parallel_with_pool_declared(
        m,
        n,
        k,
        alpha,
        a,
        b,
        beta,
        c,
        pool,
        mc_q,
        nc_q,
        gate,
        fs_exec::RunId::default(),
    )
}

/// As [`gemm_f64_parallel_with_pool`], with an explicit caller-ledgered
/// operation identity. Every NC/KC panel receives a domain-separated child
/// [`fs_exec::RunId`], so distinct operations cannot reuse corresponding tile
/// stream identities accidentally.
///
/// # Errors
/// As [`gemm_f64_parallel_with_pool`].
///
/// # Panics
/// As [`gemm_f64_parallel_with_pool`].
#[allow(clippy::too_many_arguments)]
pub fn gemm_f64_parallel_with_pool_declared(
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    beta: f64,
    c: &mut [f64],
    pool: &fs_exec::TilePool,
    mc_q: usize,
    nc_q: usize,
    gate: &fs_exec::CancelGate,
    declared_run: fs_exec::RunId,
) -> Result<GemmRunReport, GemmRunError> {
    gemm_f64_parallel_with_pool_and_poll(
        m,
        n,
        k,
        alpha,
        a,
        b,
        beta,
        c,
        pool,
        mc_q,
        nc_q,
        gate,
        declared_run,
        GemmMemoryEnvelope::UNBOUNDED,
        &|| gate.is_requested(),
    )
}

/// As [`gemm_f64_parallel_with_pool_declared`], with an explicit
/// operation MEMORY ENVELOPE (bead wf9.15): the full reservation plan
/// (C staging + shared B pack + band metadata + per-worker arena
/// panels) is preflight-checked against `envelope` BEFORE any
/// allocation or C mutation, every root reservation is fallible, and
/// the admitted plan is ledgered in [`GemmRunReport::memory`].
///
/// # Errors
/// As [`gemm_f64_parallel_with_pool`], plus
/// [`GemmRunError::MemoryRefused`] when the plan exceeds the envelope
/// or the allocator declines a reservation — C is bitwise unchanged
/// and any dispatched panels were drained.
///
/// # Panics
/// As [`gemm_f64_parallel_with_pool`].
#[allow(clippy::too_many_arguments)]
pub fn gemm_f64_parallel_with_pool_budgeted(
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    beta: f64,
    c: &mut [f64],
    pool: &fs_exec::TilePool,
    mc_q: usize,
    nc_q: usize,
    gate: &fs_exec::CancelGate,
    declared_run: fs_exec::RunId,
    envelope: GemmMemoryEnvelope,
) -> Result<GemmRunReport, GemmRunError> {
    gemm_f64_parallel_with_pool_and_poll(
        m,
        n,
        k,
        alpha,
        a,
        b,
        beta,
        c,
        pool,
        mc_q,
        nc_q,
        gate,
        declared_run,
        envelope,
        &|| gate.is_requested(),
    )
}

#[track_caller]
fn segmented_micro_tiles(extent: usize, block: usize, micro: usize, label: &str) -> usize {
    debug_assert!(block > 0 && micro > 0);
    let full_blocks = extent / block;
    let tail = extent % block;
    let per_full = block.div_ceil(micro);
    full_blocks
        .checked_mul(per_full)
        .and_then(|count| count.checked_add(tail.div_ceil(micro)))
        .unwrap_or_else(|| panic!("{label} tile-count overflow"))
}

#[track_caller]
fn cancellable_tile_count(m: usize, n: usize, k: usize, mc_q: usize, nc_q: usize) -> usize {
    let mt = segmented_micro_tiles(m, mc_q, MR, "parallel M");
    let nt = segmented_micro_tiles(n, nc_q, NR, "parallel N");
    let kt = k.div_ceil(KC);
    checked_product("parallel MN tile count", mt, nt)
        .checked_mul(kt)
        .unwrap_or_else(|| panic!("parallel MNK tile-count overflow"))
}

fn report(
    completed: &std::sync::atomic::AtomicUsize,
    total_tiles: usize,
    pool_runs: Vec<fs_exec::RunReport>,
    declared_run: fs_exec::RunId,
    memory: GemmMemoryReport,
) -> GemmRunReport {
    GemmRunReport {
        declared_run,
        completed_tiles: completed.load(std::sync::atomic::Ordering::Acquire),
        total_tiles,
        pool_runs,
        memory,
    }
}

fn cancelled(
    completed: &std::sync::atomic::AtomicUsize,
    total_tiles: usize,
    pool_runs: Vec<fs_exec::RunReport>,
    declared_run: fs_exec::RunId,
    memory: GemmMemoryReport,
) -> GemmCancelled {
    GemmCancelled {
        report: Box::new(report(
            completed,
            total_tiles,
            pool_runs,
            declared_run,
            memory,
        )),
    }
}

fn alloc_error_requested_bytes(error: &fs_alloc::AllocError) -> u128 {
    match error {
        fs_alloc::AllocError::Exhausted {
            requested_bytes, ..
        }
        | fs_alloc::AllocError::OutOfMemory {
            requested_bytes, ..
        } => *requested_bytes as u128,
        fs_alloc::AllocError::LayoutOverflow {
            len, elem_bytes, ..
        } => (*len as u128) * (*elem_bytes as u128),
        fs_alloc::AllocError::ReservationOverflow {
            base_bytes,
            additional_bytes,
            ..
        } => (*base_bytes as u128) + (*additional_bytes as u128),
    }
}

struct GemmMemoryPlan {
    report: GemmMemoryReport,
    b_pack_len: Option<usize>,
    band_count: usize,
    panel_count: Option<usize>,
}

fn checked_memory_product(what: &'static str, lhs: u128, rhs: u128) -> Result<u128, &'static str> {
    lhs.checked_mul(rhs).ok_or(what)
}

fn checked_memory_sum(
    what: &'static str,
    values: impl IntoIterator<Item = u128>,
) -> Result<u128, &'static str> {
    values
        .into_iter()
        .try_fold(0_u128, |sum, value| sum.checked_add(value).ok_or(what))
}

/// Preflight byte plan (wf9.15): every fs-la-owned reservation is computed
/// with checked arithmetic before allocation or C mutation. No-product calls
/// require only transactional C staging. Product calls size arenas from the
/// pool's exact fresh-arena reservation and the active M-band worker count.
#[allow(clippy::too_many_arguments)]
fn preflight_memory(
    c_len: usize,
    m: usize,
    n: usize,
    k: usize,
    has_product: bool,
    mc_q: usize,
    nc_q: usize,
    workers: usize,
    arena_bytes_per_worker: u64,
    envelope: GemmMemoryEnvelope,
) -> Result<GemmMemoryPlan, &'static str> {
    let f64_bytes = core::mem::size_of::<f64>() as u128;
    let staging_bytes = checked_memory_product("c-staging", c_len as u128, f64_bytes)?;
    if !has_product {
        return Ok(GemmMemoryPlan {
            report: GemmMemoryReport {
                limit_bytes: envelope.limit_bytes,
                staging_bytes,
                requested_bytes: staging_bytes,
                ..GemmMemoryReport::default()
            },
            b_pack_len: Some(0),
            band_count: 0,
            panel_count: Some(0),
        });
    }

    let b_pack_cols = (nc_q as u128)
        .checked_add((NR - 1) as u128)
        .ok_or("b-pack-columns")?
        / NR as u128
        * NR as u128;
    let b_pack_len_u128 = checked_memory_product("b-pack-elements", KC as u128, b_pack_cols)?;
    let b_pack_bytes = checked_memory_product("b-pack-bytes", b_pack_len_u128, f64_bytes)?;
    let b_pack_len = usize::try_from(b_pack_len_u128).ok();

    let band_count = m.div_ceil(mc_q);
    let band_metadata_bytes = checked_memory_product(
        "band-metadata",
        band_count as u128,
        core::mem::size_of::<std::sync::Mutex<&mut [f64]>>() as u128,
    )?;
    let panel_count_u128 = checked_memory_product(
        "panel-run-count",
        n.div_ceil(nc_q) as u128,
        k.div_ceil(KC) as u128,
    )?;
    let pool_run_bytes = checked_memory_product(
        "panel-run-receipts",
        panel_count_u128,
        core::mem::size_of::<fs_exec::RunReport>() as u128,
    )?;
    let panel_count = usize::try_from(panel_count_u128).ok();
    let active_arena_workers = workers.min(band_count);
    let arena_bytes = checked_memory_product(
        "a-pack-arenas",
        active_arena_workers as u128,
        u128::from(arena_bytes_per_worker),
    )?;
    let requested_bytes = checked_memory_sum(
        "gemm-memory-total",
        [
            staging_bytes,
            b_pack_bytes,
            band_metadata_bytes,
            pool_run_bytes,
            arena_bytes,
        ],
    )?;
    Ok(GemmMemoryPlan {
        report: GemmMemoryReport {
            limit_bytes: envelope.limit_bytes,
            staging_bytes,
            b_pack_bytes,
            band_metadata_bytes,
            pool_run_bytes,
            arena_bytes_per_worker,
            active_arena_workers,
            arena_bytes,
            requested_bytes,
            peak_used_bytes: 0,
            refused_bytes: 0,
        },
        b_pack_len,
        band_count,
        panel_count,
    })
}

#[derive(Default)]
struct ArenaUseTracker {
    current: std::sync::atomic::AtomicUsize,
    peak: std::sync::atomic::AtomicUsize,
}

impl ArenaUseTracker {
    fn enter(&self) -> ArenaUseGuard<'_> {
        let current = self
            .current
            .fetch_add(1, std::sync::atomic::Ordering::AcqRel)
            + 1;
        self.peak
            .fetch_max(current, std::sync::atomic::Ordering::AcqRel);
        ArenaUseGuard { tracker: self }
    }

    fn peak(&self) -> usize {
        self.peak.load(std::sync::atomic::Ordering::Acquire)
    }
}

struct ArenaUseGuard<'a> {
    tracker: &'a ArenaUseTracker,
}

impl Drop for ArenaUseGuard<'_> {
    fn drop(&mut self) {
        self.tracker
            .current
            .fetch_sub(1, std::sync::atomic::Ordering::AcqRel);
    }
}

fn refresh_peak_memory(
    memory: &mut GemmMemoryReport,
    root_used_bytes: u128,
    arena_use: &ArenaUseTracker,
) {
    let live =
        root_used_bytes + (arena_use.peak() as u128) * u128::from(memory.arena_bytes_per_worker);
    memory.peak_used_bytes = memory.peak_used_bytes.max(live);
}

#[allow(clippy::too_many_arguments)]
fn memory_refused(
    what: &'static str,
    refused_bytes: u128,
    completed: &std::sync::atomic::AtomicUsize,
    total_tiles: usize,
    pool_runs: Vec<fs_exec::RunReport>,
    declared_run: fs_exec::RunId,
    mut memory: GemmMemoryReport,
    root_used_bytes: u128,
    arena_use: &ArenaUseTracker,
) -> GemmRunError {
    refresh_peak_memory(&mut memory, root_used_bytes, arena_use);
    memory.refused_bytes = refused_bytes;
    GemmRunError::MemoryRefused {
        what,
        requested_bytes: refused_bytes,
        limit_bytes: memory.limit_bytes,
        report: Box::new(report(
            completed,
            total_tiles,
            pool_runs,
            declared_run,
            memory,
        )),
    }
}

struct GemmBandKernel<'kernel, 'staged, 'shared, P> {
    bands: &'kernel [std::sync::Mutex<&'staged mut [f64]>],
    m: usize,
    n: usize,
    k: usize,
    ic_quantum: usize,
    pc: usize,
    jc: usize,
    kc: usize,
    nc: usize,
    alpha: f64,
    a: &'shared [f64],
    b_pack: &'kernel [f64],
    completed: &'shared std::sync::atomic::AtomicUsize,
    arena_bytes_per_worker: usize,
    arena_use: &'shared ArenaUseTracker,
    poll: &'shared P,
}

impl<P> fs_exec::TileKernel for GemmBandKernel<'_, '_, '_, P>
where
    P: Fn() -> bool + Sync,
{
    type Out = ();

    fn tiles(&self) -> fs_exec::TilePlan {
        fs_exec::TilePlan::new(
            "fs-la/gemm-f64-m-band-v1",
            u64::try_from(self.bands.len()).expect("GEMM M-band count exceeds u64"),
        )
    }

    fn run(
        &self,
        tile: u64,
        cx: &fs_exec::Cx<'_>,
    ) -> core::ops::ControlFlow<fs_exec::Cancelled, ()> {
        if cx.checkpoint().is_err() || (self.poll)() {
            return core::ops::ControlFlow::Break(fs_exec::Cancelled);
        }
        let tile = usize::try_from(tile).expect("GEMM tile index exceeds usize");
        let ic = tile
            .checked_mul(self.ic_quantum)
            .expect("GEMM M-band offset overflow");
        let mc = self.ic_quantum.min(self.m - ic);
        let mut band = self.bands[tile].lock().expect("GEMM C-band lock poisoned");

        // The run boundary must carry a finite quota equal to one fresh arena
        // reservation. This executable check makes an accidental regression to
        // TilePool::run_declared (Budget::INFINITE) fail closed in the same typed
        // memory channel as an allocator refusal.
        let Some(cost_quota) = cx.budget().remaining_cost() else {
            let error = fs_alloc::AllocError::Exhausted {
                site: "fs-la/gemm-a-micro-panel-budget",
                requested_bytes: self.arena_bytes_per_worker,
                reserved_bytes: 0,
                limit_bytes: 0,
            };
            return core::ops::ControlFlow::Break(
                cx.refuse(fs_exec::TileFailure::Allocation(error)),
            );
        };
        let expected_cost = u64::try_from(self.arena_bytes_per_worker)
            .expect("the preflight rejected arena reservations above u64");
        if cost_quota < expected_cost {
            let error = fs_alloc::AllocError::Exhausted {
                site: "fs-la/gemm-a-micro-panel-budget",
                requested_bytes: self.arena_bytes_per_worker,
                reserved_bytes: 0,
                limit_bytes: usize::try_from(cost_quota)
                    .expect("quota below a usize-derived arena reservation fits usize"),
            };
            return core::ops::ControlFlow::Break(
                cx.refuse(fs_exec::TileFailure::Allocation(error)),
            );
        }

        // One bounded A micro-panel lives in the Cx arena. It is reclaimed
        // when this tile completes, cancels, refuses, or panics; no packing
        // allocation can escape the executor scope.
        let arena_use = self.arena_use.enter();
        let a_pack = match cx.arena().alloc_slice_fill(
            fs_alloc::Site::named("fs-la/gemm-a-micro-panel"),
            MR * KC,
            0.0,
        ) {
            Ok(pack) => pack,
            Err(error) => {
                return core::ops::ControlFlow::Break(
                    cx.refuse(fs_exec::TileFailure::Allocation(error)),
                );
            }
        };
        let _arena_use = arena_use;
        let poll = || cx.checkpoint().is_err() || (self.poll)();
        let mut p = 0;
        while p < mc {
            let rows = MR.min(mc - p);
            if !pack_a_with_poll(
                a_pack,
                self.a,
                self.k,
                ic + p,
                self.pc,
                rows,
                self.kc,
                &poll,
            ) || !macro_kernel_with_poll(
                a_pack,
                self.b_pack,
                &mut band[p * self.n..],
                self.n,
                self.jc,
                rows,
                self.nc,
                self.kc,
                self.alpha,
                self.completed,
                &poll,
            ) {
                return core::ops::ControlFlow::Break(fs_exec::Cancelled);
            }
            p += MR;
        }
        core::ops::ControlFlow::Continue(())
    }
}

/// Internal generic poll seam. Production passes a monotonic CancelGate poll;
/// the generic form lets G4 deterministically inject a request after a known
/// number of boundaries without timing sleeps.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn gemm_f64_parallel_with_pool_and_poll<P>(
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    b: &[f64],
    beta: f64,
    c: &mut [f64],
    pool: &fs_exec::TilePool,
    mc_q: usize,
    nc_q: usize,
    gate: &fs_exec::CancelGate,
    declared_run: fs_exec::RunId,
    envelope: GemmMemoryEnvelope,
    poll: &P,
) -> Result<GemmRunReport, GemmRunError>
where
    P: Fn() -> bool + Sync,
{
    // Shape rejection precedes every allocation and every possible mutation.
    assert_contiguous_shapes(m, n, k, a, b, c);
    let mc_q = mc_q.max(MR).min(m.max(MR));
    let nc_q = nc_q.max(NR).min(n.max(NR));
    let has_product = m != 0 && n != 0 && k != 0 && alpha != 0.0;
    let total_tiles = if has_product {
        cancellable_tile_count(m, n, k, mc_q, nc_q)
    } else {
        0
    };
    let completed = std::sync::atomic::AtomicUsize::new(0);
    let arena_use = ArenaUseTracker::default();
    // MEMORY PREFLIGHT (wf9.15): the full reservation plan is checked
    // against the caller's envelope BEFORE anything is allocated or
    // mutated — a refusal here has touched nothing.
    let arena_reservation = if has_product {
        match pool.arena_pool().reservation_bytes_for_slice::<f64>(
            fs_alloc::Site::named("fs-la/gemm-a-micro-panel"),
            MR * KC,
        ) {
            Ok(bytes) => bytes,
            Err(error) => {
                let refused_bytes = alloc_error_requested_bytes(&error);
                let memory = GemmMemoryReport {
                    limit_bytes: envelope.limit_bytes,
                    requested_bytes: refused_bytes,
                    refused_bytes,
                    ..GemmMemoryReport::default()
                };
                return Err(GemmRunError::MemoryRefused {
                    what: "a-pack-arena-plan",
                    requested_bytes: refused_bytes,
                    limit_bytes: envelope.limit_bytes,
                    report: Box::new(report(
                        &completed,
                        total_tiles,
                        Vec::new(),
                        declared_run,
                        memory,
                    )),
                });
            }
        }
    } else {
        0
    };
    let arena_reservation_u64 =
        u64::try_from(arena_reservation).map_err(|_| GemmRunError::MemoryPlanOverflow {
            what: "a-pack-arena-reservation",
            limit_bytes: envelope.limit_bytes,
        })?;
    let plan = preflight_memory(
        c.len(),
        m,
        n,
        k,
        has_product,
        mc_q,
        nc_q,
        pool.workers(),
        arena_reservation_u64,
        envelope,
    )
    .map_err(|what| GemmRunError::MemoryPlanOverflow {
        what,
        limit_bytes: envelope.limit_bytes,
    })?;
    let mut memory = plan.report;
    if memory.requested_bytes > u128::from(memory.limit_bytes) {
        memory.refused_bytes = memory.requested_bytes;
        return Err(GemmRunError::MemoryRefused {
            what: "preflight-envelope",
            requested_bytes: memory.requested_bytes,
            limit_bytes: memory.limit_bytes,
            report: Box::new(report(
                &completed,
                total_tiles,
                Vec::new(),
                declared_run,
                memory,
            )),
        });
    }
    let Some(b_pack_len) = plan.b_pack_len else {
        return Err(memory_refused(
            "b-pack-layout",
            memory.b_pack_bytes,
            &completed,
            total_tiles,
            Vec::new(),
            declared_run,
            memory,
            0,
            &arena_use,
        ));
    };
    let Some(panel_count) = plan.panel_count else {
        return Err(memory_refused(
            "panel-run-layout",
            memory.pool_run_bytes,
            &completed,
            total_tiles,
            Vec::new(),
            declared_run,
            memory,
            0,
            &arena_use,
        ));
    };
    if poll() {
        return Err(GemmRunError::Cancelled(cancelled(
            &completed,
            total_tiles,
            Vec::new(),
            declared_run,
            memory,
        )));
    }

    let mut root_used_bytes = 0_u128;
    let mut pool_runs = Vec::new();
    if pool_runs.try_reserve_exact(panel_count).is_err() {
        return Err(memory_refused(
            "panel-run-receipts",
            memory.pool_run_bytes,
            &completed,
            total_tiles,
            pool_runs,
            declared_run,
            memory,
            root_used_bytes,
            &arena_use,
        ));
    }
    root_used_bytes += memory.pool_run_bytes;
    refresh_peak_memory(&mut memory, root_used_bytes, &arena_use);

    // Transactional staging is the no-torn-C boundary. Capacity reservation
    // itself is not a poll point, but initialization/copying is chunked under
    // the gate so beta=0 cannot hide an unbounded zero-fill.
    let mut staged = match stage_beta(c, beta, poll) {
        Ok(staged) => staged,
        Err(StageAbort::AllocRefused) => {
            return Err(memory_refused(
                "c-staging",
                memory.staging_bytes,
                &completed,
                total_tiles,
                pool_runs,
                declared_run,
                memory,
                root_used_bytes,
                &arena_use,
            ));
        }
        Err(StageAbort::Cancelled) => {
            root_used_bytes += memory.staging_bytes;
            refresh_peak_memory(&mut memory, root_used_bytes, &arena_use);
            return Err(GemmRunError::Cancelled(cancelled(
                &completed,
                total_tiles,
                pool_runs,
                declared_run,
                memory,
            )));
        }
    };
    root_used_bytes += memory.staging_bytes;
    refresh_peak_memory(&mut memory, root_used_bytes, &arena_use);
    if !has_product {
        if poll() {
            return Err(GemmRunError::Cancelled(cancelled(
                &completed,
                total_tiles,
                pool_runs,
                declared_run,
                memory,
            )));
        }
        c.copy_from_slice(&staged);
        return Ok(report(
            &completed,
            total_tiles,
            pool_runs,
            declared_run,
            memory,
        ));
    }

    let band_len = checked_product("parallel C band", mc_q, n);
    let mut b_pack = match zeroed_with_poll(b_pack_len, poll) {
        Ok(pack) => pack,
        Err(StageAbort::AllocRefused) => {
            return Err(memory_refused(
                "b-pack",
                memory.b_pack_bytes,
                &completed,
                total_tiles,
                pool_runs,
                declared_run,
                memory,
                root_used_bytes,
                &arena_use,
            ));
        }
        Err(StageAbort::Cancelled) => {
            root_used_bytes += memory.b_pack_bytes;
            refresh_peak_memory(&mut memory, root_used_bytes, &arena_use);
            return Err(GemmRunError::Cancelled(cancelled(
                &completed,
                total_tiles,
                pool_runs,
                declared_run,
                memory,
            )));
        }
    };
    root_used_bytes += memory.b_pack_bytes;
    refresh_peak_memory(&mut memory, root_used_bytes, &arena_use);

    let mut bands: Vec<std::sync::Mutex<&mut [f64]>> = Vec::new();
    if bands.try_reserve_exact(plan.band_count).is_err() {
        return Err(memory_refused(
            "band-metadata",
            memory.band_metadata_bytes,
            &completed,
            total_tiles,
            pool_runs,
            declared_run,
            memory,
            root_used_bytes,
            &arena_use,
        ));
    }
    bands.extend(staged.chunks_mut(band_len).map(std::sync::Mutex::new));
    debug_assert_eq!(bands.len(), plan.band_count);
    root_used_bytes += memory.band_metadata_bytes;
    refresh_peak_memory(&mut memory, root_used_bytes, &arena_use);

    let mut jc = 0;
    while jc < n {
        let nc = nc_q.min(n - jc);
        let mut pc = 0;
        while pc < k {
            let kc = KC.min(k - pc);
            if !pack_b_with_poll(&mut b_pack, b, n, pc, jc, kc, nc, poll) {
                return Err(GemmRunError::Cancelled(cancelled(
                    &completed,
                    total_tiles,
                    pool_runs,
                    declared_run,
                    memory,
                )));
            }
            let kernel = GemmBandKernel {
                bands: &bands,
                m,
                n,
                k,
                ic_quantum: mc_q,
                pc,
                jc,
                kc,
                nc,
                alpha,
                a,
                b_pack: &b_pack,
                completed: &completed,
                arena_bytes_per_worker: arena_reservation,
                arena_use: &arena_use,
                poll,
            };
            let panel_ordinal = u64::try_from(pool_runs.len())
                .expect("supported Rust targets have at most 64-bit usize");
            let panel_run = gemm_panel_run_id(declared_run, panel_ordinal);
            let budget = fs_exec::Budget::new().with_cost_quota(memory.arena_bytes_per_worker);
            let (outcome, pool_report) =
                pool.run_declared_budgeted(&kernel, gate, panel_run, budget);
            pool_runs.push(pool_report);
            refresh_peak_memory(&mut memory, root_used_bytes, &arena_use);
            match outcome {
                Ok(()) => {}
                Err(fs_exec::RunError::Cancelled { .. }) => {
                    return Err(GemmRunError::Cancelled(cancelled(
                        &completed,
                        total_tiles,
                        pool_runs,
                        declared_run,
                        memory,
                    )));
                }
                Err(fs_exec::RunError::TileFailed {
                    failure: fs_exec::TileFailure::Allocation(error),
                    ..
                }) => {
                    let refused_bytes = alloc_error_requested_bytes(&error);
                    return Err(memory_refused(
                        "a-pack-arena",
                        refused_bytes,
                        &completed,
                        total_tiles,
                        pool_runs,
                        declared_run,
                        memory,
                        root_used_bytes,
                        &arena_use,
                    ));
                }
                Err(error) => {
                    return Err(GemmRunError::Executor {
                        error,
                        report: Box::new(report(
                            &completed,
                            total_tiles,
                            pool_runs,
                            declared_run,
                            memory,
                        )),
                    });
                }
            }
            if poll() {
                return Err(GemmRunError::Cancelled(cancelled(
                    &completed,
                    total_tiles,
                    pool_runs,
                    declared_run,
                    memory,
                )));
            }
            pc += KC;
        }
        jc += nc_q;
    }

    drop(bands);
    refresh_peak_memory(&mut memory, root_used_bytes, &arena_use);
    let final_report = report(&completed, total_tiles, pool_runs, declared_run, memory);
    debug_assert_eq!(final_report.completed_tiles, final_report.total_tiles);
    if poll() {
        return Err(GemmRunError::Cancelled(GemmCancelled {
            report: Box::new(final_report),
        }));
    }
    // FINALIZE: &mut C excludes safe concurrent observers. Once the final poll
    // passes, copy is deliberately non-cancellable and the call returns a
    // complete result even if a later request races this commit.
    c.copy_from_slice(&staged);
    Ok(final_report)
}

/// Why a root-orchestration reservation stopped (wf9.15).
enum StageAbort {
    /// The gate/poll tripped mid-initialization.
    Cancelled,
    /// The allocator declined the reservation (never a process abort).
    AllocRefused,
}

/// Fallible, poll-chunked zeroed buffer: capacity via try_reserve_exact
/// (allocator refusal is a STRUCTURED outcome, not an abort), fill
/// chunked under the gate.
fn zeroed_with_poll<P>(len: usize, poll: &P) -> Result<Vec<f64>, StageAbort>
where
    P: Fn() -> bool,
{
    let mut values = Vec::new();
    values
        .try_reserve_exact(len)
        .map_err(|_| StageAbort::AllocRefused)?;
    while values.len() < len {
        if poll() {
            return Err(StageAbort::Cancelled);
        }
        let next = (len - values.len()).min(C_STAGE_TILE_ELEMENTS);
        values.resize(values.len() + next, 0.0);
    }
    if poll() {
        return Err(StageAbort::Cancelled);
    }
    Ok(values)
}

fn stage_beta<P>(source: &[f64], beta: f64, poll: &P) -> Result<Vec<f64>, StageAbort>
where
    P: Fn() -> bool,
{
    let mut staged = Vec::new();
    staged
        .try_reserve_exact(source.len())
        .map_err(|_| StageAbort::AllocRefused)?;
    for src in source.chunks(C_STAGE_TILE_ELEMENTS) {
        if poll() {
            return Err(StageAbort::Cancelled);
        }
        if beta == 0.0 {
            staged.resize(staged.len() + src.len(), 0.0);
        } else if beta.to_bits() == 1.0f64.to_bits() {
            staged.extend_from_slice(src);
        } else {
            staged.extend(src.iter().map(|&value| value * beta));
        }
    }
    if poll() {
        return Err(StageAbort::Cancelled);
    }
    Ok(staged)
}

#[allow(clippy::too_many_arguments)]
fn pack_a_with_poll<P>(
    dst: &mut [f64],
    a: &[f64],
    lda: usize,
    ic: usize,
    pc: usize,
    mc: usize,
    kc: usize,
    poll: &P,
) -> bool
where
    P: Fn() -> bool,
{
    let mut w = 0;
    let mut p = 0;
    while p < mc {
        if poll() {
            return false;
        }
        let rows = MR.min(mc - p);
        for kk in 0..kc {
            for r in 0..MR {
                dst[w] = if r < rows {
                    a[(ic + p + r) * lda + pc + kk]
                } else {
                    0.0
                };
                w += 1;
            }
        }
        p += MR;
    }
    !poll()
}

#[allow(clippy::too_many_arguments)]
fn pack_b_with_poll<P>(
    dst: &mut [f64],
    b: &[f64],
    ldb: usize,
    pc: usize,
    jc: usize,
    kc: usize,
    nc: usize,
    poll: &P,
) -> bool
where
    P: Fn() -> bool,
{
    let mut w = 0;
    let mut q = 0;
    while q < nc {
        if poll() {
            return false;
        }
        let cols = NR.min(nc - q);
        for kk in 0..kc {
            for s in 0..NR {
                dst[w] = if s < cols {
                    b[(pc + kk) * ldb + jc + q + s]
                } else {
                    0.0
                };
                w += 1;
            }
        }
        q += NR;
    }
    !poll()
}

#[allow(clippy::too_many_arguments)]
fn macro_kernel_with_poll<P>(
    a_pack: &[f64],
    b_pack: &[f64],
    c: &mut [f64],
    n: usize,
    jc: usize,
    mc: usize,
    nc: usize,
    kc: usize,
    alpha: f64,
    completed: &std::sync::atomic::AtomicUsize,
    poll: &P,
) -> bool
where
    P: Fn() -> bool,
{
    let mut p = 0;
    while p < mc {
        let rows = MR.min(mc - p);
        let a_panel = &a_pack[(p / MR) * MR * kc..][..MR * kc];
        let mut q = 0;
        while q < nc {
            if poll() {
                return false;
            }
            let cols = NR.min(nc - q);
            let b_panel = &b_pack[(q / NR) * NR * kc..][..NR * kc];
            let mut acc = [[0.0f64; NR]; MR];
            (fs_simd::ops().mk8x4_f64)(a_panel, b_panel, kc, &mut acc);
            for (r, accr) in acc.iter().enumerate().take(rows) {
                let crow = (p + r) * n + jc + q;
                for (s, &value) in accr.iter().enumerate().take(cols) {
                    c[crow + s] = alpha.mul_add(value, c[crow + s]);
                }
            }
            completed.fetch_add(1, std::sync::atomic::Ordering::Release);
            q += NR;
        }
        p += MR;
    }
    !poll()
}

/// Operand orientation for the op-form GEMM entry points.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Trans {
    /// Use the operand as stored.
    N,
    /// Use the operand transposed.
    T,
}

/// f64 GEMM over TRANSPOSED/STRIDED operand views (xlvx s7):
/// `C = α·op(A)·op(B) + β·C` where `op(X)` is `X` or `Xᵀ`, each operand
/// carrying its own leading dimension (row stride of the STORED
/// matrix), so submatrix views compute without copies.
///
/// BIT CONTRACT: op() and the leading dimensions are absorbed entirely
/// by pack addressing — the packed panels are byte-identical to the
/// contiguous non-transposed equivalent, hence the OUTPUT IS BITWISE
/// [`gemm_f64`] on materialized operands (gated in-module). Rows of C
/// outside the m×n view are never touched.
///
/// # Panics
/// Structured panics when a leading dimension is too small or a slice
/// cannot hold its view.
#[allow(clippy::too_many_arguments)] // BLAS-shape signature
pub fn gemm_f64_op(
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f64],
    lda: usize,
    ta: Trans,
    b: &[f64],
    ldb: usize,
    tb: Trans,
    beta: f64,
    c: &mut [f64],
    ldc: usize,
) {
    match ta {
        Trans::N => assert_view_shape("a", a.len(), m, k, lda),
        Trans::T => assert_view_shape("a", a.len(), k, m, lda),
    }
    match tb {
        Trans::N => assert_view_shape("b", b.len(), k, n, ldb),
        Trans::T => assert_view_shape("b", b.len(), n, k, ldb),
    }
    assert_view_shape("c", c.len(), m, n, ldc);
    // β pass over the VIEW only (op-form C may be a submatrix).
    if beta == 0.0 {
        for row in c.chunks_mut(ldc).take(m) {
            row[..n].fill(0.0);
        }
    } else if beta.to_bits() != 1.0f64.to_bits() {
        for row in c.chunks_mut(ldc).take(m) {
            for v in &mut row[..n] {
                *v *= beta;
            }
        }
    }
    if m == 0 || n == 0 || alpha == 0.0 || k == 0 {
        return;
    }
    let mut a_pack = vec![0.0f64; MC * KC];
    let mut b_pack = vec![0.0f64; KC * NC];
    let mut jc = 0;
    while jc < n {
        let nc = NC.min(n - jc);
        let mut pc = 0;
        while pc < k {
            let kc = KC.min(k - pc);
            pack_b_op(&mut b_pack, b, ldb, tb, pc, jc, kc, nc);
            let mut ic = 0;
            while ic < m {
                let mc = MC.min(m - ic);
                pack_a_op(&mut a_pack, a, lda, ta, ic, pc, mc, kc);
                macro_kernel(&a_pack, &b_pack, c, m, ldc, ic, jc, mc, nc, kc, alpha);
                ic += MC;
            }
            pc += KC;
        }
        jc += NC;
    }
}

/// Op-form A packer: element op(A)[ic+p+r, pc+kk] addressed through
/// (lda, ta) — N reads a[row·lda + col], T reads a[col·lda + row]. The
/// packed layout (and bytes) are exactly [`pack_a`]'s.
#[allow(clippy::too_many_arguments)]
fn pack_a_op(
    dst: &mut [f64],
    a: &[f64],
    lda: usize,
    ta: Trans,
    ic: usize,
    pc: usize,
    mc: usize,
    kc: usize,
) {
    let mut w = 0;
    let mut p = 0;
    while p < mc {
        let rows = MR.min(mc - p);
        for kk in 0..kc {
            for r in 0..MR {
                dst[w] = if r < rows {
                    match ta {
                        Trans::N => a[(ic + p + r) * lda + pc + kk],
                        Trans::T => a[(pc + kk) * lda + ic + p + r],
                    }
                } else {
                    0.0
                };
                w += 1;
            }
        }
        p += MR;
    }
}

/// Op-form B packer: element op(B)[pc+kk, jc+q+s] addressed through
/// (ldb, tb); packed bytes are exactly [`pack_b`]'s.
#[allow(clippy::too_many_arguments)]
fn pack_b_op(
    dst: &mut [f64],
    b: &[f64],
    ldb: usize,
    tb: Trans,
    pc: usize,
    jc: usize,
    kc: usize,
    nc: usize,
) {
    let mut w = 0;
    let mut q = 0;
    while q < nc {
        let cols = NR.min(nc - q);
        for kk in 0..kc {
            for s in 0..NR {
                dst[w] = if s < cols {
                    match tb {
                        Trans::N => b[(pc + kk) * ldb + jc + q + s],
                        Trans::T => b[(jc + q + s) * ldb + pc + kk],
                    }
                } else {
                    0.0
                };
                w += 1;
            }
        }
        q += NR;
    }
}

/// β application with BLAS overwrite semantics for β = 0.
fn scale_c(c: &mut [f64], beta: f64) {
    if beta == 0.0 {
        c.fill(0.0);
    } else if beta.to_bits() != 1.0f64.to_bits() {
        for v in c.iter_mut() {
            *v *= beta;
        }
    }
}

/// Pack an mc×kc block of A (row-major, ld = k) into MR-row micro-panels:
/// panel p holds rows [p·MR, p·MR+MR) column-major-within-panel
/// (k-index fastest across the MR lanes). Short tail rows are zero-padded —
/// zero lanes contribute exact +0.0 products which never reach C (tail
/// handling masks them on write-back).
fn pack_a(dst: &mut [f64], a: &[f64], lda: usize, ic: usize, pc: usize, mc: usize, kc: usize) {
    let mut w = 0;
    let mut p = 0;
    while p < mc {
        let rows = MR.min(mc - p);
        for kk in 0..kc {
            for r in 0..MR {
                dst[w] = if r < rows {
                    a[(ic + p + r) * lda + pc + kk]
                } else {
                    0.0
                };
                w += 1;
            }
        }
        p += MR;
    }
}

/// Pack a kc×nc block of B (row-major, ld = n) into NR-column micro-panels
/// (k-index outer, NR lanes inner), zero-padded tails.
fn pack_b(dst: &mut [f64], b: &[f64], ldb: usize, pc: usize, jc: usize, kc: usize, nc: usize) {
    let mut w = 0;
    let mut q = 0;
    while q < nc {
        let cols = NR.min(nc - q);
        for kk in 0..kc {
            for s in 0..NR {
                dst[w] = if s < cols {
                    b[(pc + kk) * ldb + jc + q + s]
                } else {
                    0.0
                };
                w += 1;
            }
        }
        q += NR;
    }
}

/// The macro kernel: sweep micro-tiles of the packed panels.
#[allow(clippy::too_many_arguments)]
fn macro_kernel(
    a_pack: &[f64],
    b_pack: &[f64],
    c: &mut [f64],
    _m: usize,
    n: usize,
    ic: usize,
    jc: usize,
    mc: usize,
    nc: usize,
    kc: usize,
    alpha: f64,
) {
    let mut p = 0;
    while p < mc {
        let rows = MR.min(mc - p);
        let a_panel = &a_pack[(p / MR) * MR * kc..][..MR * kc];
        let mut q = 0;
        while q < nc {
            let cols = NR.min(nc - q);
            let b_panel = &b_pack[(q / NR) * NR * kc..][..NR * kc];
            // Register-tiled microkernel: MR×NR accumulators, k
            // ascending — through the fs-simd dispatch table (bead
            // xdgf). The NEON capsule is per-element bitwise-identical
            // to the scalar twin (which IS the former inline loop), so
            // the golden hash is tier-invariant.
            let mut acc = [[0.0f64; NR]; MR];
            (fs_simd::ops().mk8x4_f64)(a_panel, b_panel, kc, &mut acc);
            // Write-back with α, masking padded tail lanes.
            for (r, accr) in acc.iter().enumerate().take(rows) {
                let crow = (ic + p + r) * n + jc + q;
                for (s, &av) in accr.iter().enumerate().take(cols) {
                    c[crow + s] = alpha.mul_add(av, c[crow + s]);
                }
            }
            q += NR;
        }
        p += MR;
    }
}

/// f32 GEMM, same structure and contract (KC shared).
#[allow(clippy::too_many_arguments)] // BLAS-shape signature: m,n,k,alpha,a,b,beta,c
pub fn gemm_f32(
    m: usize,
    n: usize,
    k: usize,
    alpha: f32,
    a: &[f32],
    b: &[f32],
    beta: f32,
    c: &mut [f32],
) {
    assert_contiguous_shapes(m, n, k, a, b, c);
    if beta == 0.0 {
        c.fill(0.0);
    } else if beta.to_bits() != 1.0f32.to_bits() {
        for v in c.iter_mut() {
            *v *= beta;
        }
    }
    if m == 0 || n == 0 || k == 0 || alpha == 0.0 {
        return;
    }
    // PACKED path (xlvx s6): same BLIS nest as f64. Bitwise identical
    // to the former naive-chunked loop — per-element accumulation is
    // still k-ascending within each KC chunk with chunk partials folded
    // into C in chunk order; packing changes layout, never arithmetic
    // (gated vs the naive-chunked oracle in-module).
    let mut a_pack = vec![0.0f32; MC * KC];
    let mut b_pack = vec![0.0f32; KC * NC];
    let mut jc = 0;
    while jc < n {
        let nc = NC.min(n - jc);
        let mut pc = 0;
        while pc < k {
            let kc = KC.min(k - pc);
            pack_b_f32(&mut b_pack, b, n, pc, jc, kc, nc);
            let mut ic = 0;
            while ic < m {
                let mc = MC.min(m - ic);
                pack_a_f32(&mut a_pack, a, k, ic, pc, mc, kc);
                macro_kernel_f32(&a_pack, &b_pack, c, n, ic, jc, mc, nc, kc, alpha);
                ic += MC;
            }
            pc += KC;
        }
        jc += NC;
    }
}

/// f32 twin of [`pack_a`]: MR-row micro-panels, zero-padded tails.
fn pack_a_f32(dst: &mut [f32], a: &[f32], lda: usize, ic: usize, pc: usize, mc: usize, kc: usize) {
    let mut w = 0;
    let mut p = 0;
    while p < mc {
        let rows = MR.min(mc - p);
        for kk in 0..kc {
            for r in 0..MR {
                dst[w] = if r < rows {
                    a[(ic + p + r) * lda + pc + kk]
                } else {
                    0.0
                };
                w += 1;
            }
        }
        p += MR;
    }
}

/// f32 twin of [`pack_b`]: NR-column micro-panels, zero-padded tails.
fn pack_b_f32(dst: &mut [f32], b: &[f32], ldb: usize, pc: usize, jc: usize, kc: usize, nc: usize) {
    let mut w = 0;
    let mut q = 0;
    while q < nc {
        let cols = NR.min(nc - q);
        for kk in 0..kc {
            for s in 0..NR {
                dst[w] = if s < cols {
                    b[(pc + kk) * ldb + jc + q + s]
                } else {
                    0.0
                };
                w += 1;
            }
        }
        q += NR;
    }
}

/// f32 macro kernel: scalar MR×NR register tile, k ascending — the
/// fs-simd f32 capsule microkernel is a recorded follow-up; this scalar
/// twin fixes the bit contract it will have to match.
#[allow(clippy::too_many_arguments)]
fn macro_kernel_f32(
    a_pack: &[f32],
    b_pack: &[f32],
    c: &mut [f32],
    n: usize,
    ic: usize,
    jc: usize,
    mc: usize,
    nc: usize,
    kc: usize,
    alpha: f32,
) {
    let mut p = 0;
    while p < mc {
        let rows = MR.min(mc - p);
        let a_panel = &a_pack[(p / MR) * MR * kc..][..MR * kc];
        let mut q = 0;
        while q < nc {
            let cols = NR.min(nc - q);
            let b_panel = &b_pack[(q / NR) * NR * kc..][..NR * kc];
            let mut acc = [[0.0f32; NR]; MR];
            for kk in 0..kc {
                for (r, accr) in acc.iter_mut().enumerate() {
                    let av = a_panel[kk * MR + r];
                    for (s, slot) in accr.iter_mut().enumerate() {
                        *slot = av.mul_add(b_panel[kk * NR + s], *slot);
                    }
                }
            }
            for (r, accr) in acc.iter().enumerate().take(rows) {
                let crow = (ic + p + r) * n + jc + q;
                for (s, &av) in accr.iter().enumerate().take(cols) {
                    c[crow + s] = alpha.mul_add(av, c[crow + s]);
                }
            }
            q += NR;
        }
        p += MR;
    }
}

/// Mixed-precision GEMM: f32 STORAGE, f64 ACCUMULATION — the
/// bandwidth-vs-accuracy compromise used throughout the plan (§6.1). Each
/// f32 element is widened exactly (f32→f64 is exact); all arithmetic is
/// f64 mul_add in the same k-ascending order.
#[allow(clippy::too_many_arguments)] // BLAS-shape signature: m,n,k,alpha,a,b,beta,c
pub fn gemm_mixed(
    m: usize,
    n: usize,
    k: usize,
    alpha: f64,
    a: &[f32],
    b: &[f32],
    beta: f64,
    c: &mut [f64],
) {
    assert_contiguous_shapes(m, n, k, a, b, c);
    scale_c(c, beta);
    if m == 0 || n == 0 || k == 0 || alpha == 0.0 {
        return;
    }
    // PACKED path (xlvx s6): panels stay f32 in memory (the bandwidth
    // saving is the point of mixed) and widen exactly at the multiply.
    // Bitwise identical to the former naive-chunked loop — same
    // per-element f64 k-ascending order (gated in-module).
    let mut a_pack = vec![0.0f32; MC * KC];
    let mut b_pack = vec![0.0f32; KC * NC];
    let mut jc = 0;
    while jc < n {
        let nc = NC.min(n - jc);
        let mut pc = 0;
        while pc < k {
            let kc = KC.min(k - pc);
            pack_b_f32(&mut b_pack, b, n, pc, jc, kc, nc);
            let mut ic = 0;
            while ic < m {
                let mc = MC.min(m - ic);
                pack_a_f32(&mut a_pack, a, k, ic, pc, mc, kc);
                macro_kernel_mixed(&a_pack, &b_pack, c, n, ic, jc, mc, nc, kc, alpha);
                ic += MC;
            }
            pc += KC;
        }
        jc += NC;
    }
}

/// Mixed macro kernel: f32 panels, f64 MR×NR accumulators; each lane
/// widens exactly (f32→f64) at the multiply, k ascending.
#[allow(clippy::too_many_arguments)]
fn macro_kernel_mixed(
    a_pack: &[f32],
    b_pack: &[f32],
    c: &mut [f64],
    n: usize,
    ic: usize,
    jc: usize,
    mc: usize,
    nc: usize,
    kc: usize,
    alpha: f64,
) {
    let mut p = 0;
    while p < mc {
        let rows = MR.min(mc - p);
        let a_panel = &a_pack[(p / MR) * MR * kc..][..MR * kc];
        let mut q = 0;
        while q < nc {
            let cols = NR.min(nc - q);
            let b_panel = &b_pack[(q / NR) * NR * kc..][..NR * kc];
            let mut acc = [[0.0f64; NR]; MR];
            for kk in 0..kc {
                for (r, accr) in acc.iter_mut().enumerate() {
                    let av = f64::from(a_panel[kk * MR + r]);
                    for (s, slot) in accr.iter_mut().enumerate() {
                        *slot = av.mul_add(f64::from(b_panel[kk * NR + s]), *slot);
                    }
                }
            }
            for (r, accr) in acc.iter().enumerate().take(rows) {
                let crow = (ic + p + r) * n + jc + q;
                for (s, &av) in accr.iter().enumerate().take(cols) {
                    c[crow + s] = alpha.mul_add(av, c[crow + s]);
                }
            }
            q += NR;
        }
        p += MR;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn root_reservations_refuse_structurally_never_abort() {
        // wf9.15: an absurd reservation is a STRUCTURED refusal from
        // try_reserve_exact (capacity overflow / allocator refusal) —
        // the process must not abort. This is the fault-injection seam
        // for the in-API MemoryRefused paths, which route through the
        // same helpers.
        assert!(matches!(
            zeroed_with_poll(usize::MAX / 16, &|| false),
            Err(StageAbort::AllocRefused)
        ));
        // And the poll path still cancels cleanly.
        assert!(matches!(
            zeroed_with_poll(1 << 20, &|| true),
            Err(StageAbort::Cancelled)
        ));
    }

    #[test]
    fn build_fingerprint_is_stable_and_canonical_for_this_binary() {
        let first = gemm_build_identity();
        let second = gemm_build_identity();
        assert_eq!(first, second, "one build has one stable codegen identity");
        assert_eq!(first.len(), 64, "BLAKE3 identity is full width");
        assert!(
            first
                .bytes()
                .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte)),
            "build identity must be canonical lowercase hex"
        );
    }

    fn lcg(seed: &mut u64) -> f64 {
        *seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((*seed >> 11) as f64) / (1u64 << 53) as f64 - 0.5
    }

    /// The oracle: naive triple loop with the SAME KC chunking and fused
    /// arithmetic — bitwise-comparable; plus a plain tolerance oracle.
    #[allow(clippy::too_many_arguments)] // mirrors the BLAS-shape signature
    fn naive_chunked(
        m: usize,
        n: usize,
        k: usize,
        alpha: f64,
        a: &[f64],
        b: &[f64],
        beta: f64,
        c0: &[f64],
    ) -> Vec<f64> {
        let mut c: Vec<f64> = if beta == 0.0 {
            vec![0.0; m * n]
        } else {
            c0.iter()
                .map(|&v| {
                    if beta.to_bits() == 1.0f64.to_bits() {
                        v
                    } else {
                        v * beta
                    }
                })
                .collect()
        };
        let mut pc = 0;
        while pc < k {
            let kc = KC.min(k - pc);
            for i in 0..m {
                for j in 0..n {
                    let mut acc = 0.0f64;
                    for kk in 0..kc {
                        acc = a[i * k + pc + kk].mul_add(b[(pc + kk) * n + j], acc);
                    }
                    c[i * n + j] = alpha.mul_add(acc, c[i * n + j]);
                }
            }
            pc += KC;
        }
        c
    }

    fn rand_mat(rows: usize, cols: usize, seed: u64) -> Vec<f64> {
        let mut s = seed;
        (0..rows * cols).map(|_| lcg(&mut s)).collect()
    }

    #[test]
    fn matches_oracle_bitwise_across_shape_sweep() {
        // The packed/blocked path must be BIT-IDENTICAL to the same-order
        // naive path: packing must not change arithmetic, only layout.
        let shapes = [
            (1usize, 1usize, 1usize),
            (1, 7, 3),
            (5, 1, 9),
            (8, 4, 256),   // exactly one micro-tile, one KC chunk
            (9, 5, 257),   // tails in every dimension
            (33, 17, 300), // KC chunking engaged
            (64, 64, 64),
            (3, 200, 2), // wide
            (200, 3, 2), // tall-skinny
        ];
        for (idx, &(m, n, k)) in shapes.iter().enumerate() {
            let a = rand_mat(m, k, 0xA + idx as u64);
            let b = rand_mat(k, n, 0xB + idx as u64);
            let c0 = rand_mat(m, n, 0xC + idx as u64);
            for (alpha, beta) in [(1.0, 0.0), (2.5, 1.0), (-0.75, 0.5)] {
                let mut c = c0.clone();
                gemm_f64(m, n, k, alpha, &a, &b, beta, &mut c);
                let want = naive_chunked(m, n, k, alpha, &a, &b, beta, &c0);
                for (i, (&got, &w)) in c.iter().zip(&want).enumerate() {
                    assert_eq!(
                        got.to_bits(),
                        w.to_bits(),
                        "({m}x{n}x{k}) alpha={alpha} beta={beta} at {i}: {got} vs {w}"
                    );
                }
            }
        }
        println!(
            "{{\"suite\":\"fs-la\",\"case\":\"gemm-oracle\",\"verdict\":\"pass\",\"detail\":\"9 shapes x 3 (alpha,beta) bitwise vs same-order oracle\"}}"
        );
    }

    #[test]
    fn degenerate_and_beta_semantics() {
        // k = 0: C = beta*C, and beta = 0 must OVERWRITE garbage/NaN.
        let mut c = vec![f64::NAN, 3.0, -2.0, 1.0];
        gemm_f64(2, 2, 0, 1.0, &[], &[], 0.0, &mut c);
        assert!(
            c.iter().all(|&v| v == 0.0),
            "beta=0 must overwrite NaN: {c:?}"
        );
        let mut c2 = vec![1.0, 2.0, 3.0, 4.0];
        gemm_f64(2, 2, 0, 1.0, &[], &[], 2.0, &mut c2);
        assert_eq!(c2, vec![2.0, 4.0, 6.0, 8.0]);
        // m or n zero: no-op, no panic.
        let mut empty: Vec<f64> = vec![];
        gemm_f64(0, 5, 3, 1.0, &[], &rand_mat(3, 5, 1), 0.0, &mut empty);
        // alpha = 0 leaves beta*C.
        let a = rand_mat(2, 3, 2);
        let b = rand_mat(3, 2, 3);
        let mut c3 = vec![1.0; 4];
        gemm_f64(2, 2, 3, 0.0, &a, &b, 1.0, &mut c3);
        assert_eq!(c3, vec![1.0; 4]);
    }

    #[test]
    fn transpose_identity_and_submatrix_consistency() {
        let (m, n, k) = (24usize, 18usize, 40usize);
        let a = rand_mat(m, k, 7);
        let b = rand_mat(k, n, 8);
        // (A·B)ᵀ == Bᵀ·Aᵀ within tight tolerance (orders differ → not bitwise).
        let mut ab = vec![0.0; m * n];
        gemm_f64(m, n, k, 1.0, &a, &b, 0.0, &mut ab);
        let at: Vec<f64> = (0..k * m).map(|i| a[(i % m) * k + i / m]).collect();
        let bt: Vec<f64> = (0..n * k).map(|i| b[(i % k) * n + i / k]).collect();
        let mut btat = vec![0.0; n * m];
        gemm_f64(n, m, k, 1.0, &bt, &at, 0.0, &mut btat);
        for i in 0..m {
            for j in 0..n {
                let x = ab[i * n + j];
                let y = btat[j * m + i];
                assert!(
                    (x - y).abs() <= 1e-13 * x.abs().max(1.0),
                    "transpose identity at ({i},{j}): {x} vs {y}"
                );
            }
        }
        // Submatrix consistency: the top-left quadrant of C equals the GEMM
        // of the corresponding A rows with B (exact: row tiling is
        // bit-neutral, same k order).
        let m2 = m / 2;
        let a_top = &a[..m2 * k];
        let mut c_top = vec![0.0; m2 * n];
        gemm_f64(m2, n, k, 1.0, a_top, &b, 0.0, &mut c_top);
        for i in 0..m2 * n {
            assert_eq!(
                c_top[i].to_bits(),
                ab[i].to_bits(),
                "row-tiling changed bits at {i}"
            );
        }
    }

    /// f32 twin of `naive_chunked` — the pre-s6 unpacked implementation,
    /// kept as the bit oracle for the packed f32 path.
    fn naive_chunked_f32(
        m: usize,
        n: usize,
        k: usize,
        alpha: f32,
        a: &[f32],
        b: &[f32],
        c: &mut [f32],
    ) {
        let mut pc = 0;
        while pc < k {
            let kc = KC.min(k - pc);
            for i in 0..m {
                for j in 0..n {
                    let mut acc = 0.0f32;
                    for kk in 0..kc {
                        acc = a[i * k + pc + kk].mul_add(b[(pc + kk) * n + j], acc);
                    }
                    c[i * n + j] = alpha.mul_add(acc, c[i * n + j]);
                }
            }
            pc += KC;
        }
    }

    #[test]
    fn f32_and_mixed_paths() {
        // Shape sweep: tails in every dimension, multi-tile, KC chunking.
        for (idx, &(m, n, k)) in [(17usize, 13usize, 129usize), (9, 5, 257), (33, 17, 300)]
            .iter()
            .enumerate()
        {
            let mut s = 0x32_u64 + idx as u64;
            let af: Vec<f32> = (0..m * k).map(|_| lcg(&mut s) as f32).collect();
            let bf: Vec<f32> = (0..k * n).map(|_| lcg(&mut s) as f32).collect();
            // Mixed vs full-f64 reference on the widened inputs: mixed IS the
            // f64 computation on exactly-widened values — bitwise equal.
            let ad: Vec<f64> = af.iter().map(|&v| f64::from(v)).collect();
            let bd: Vec<f64> = bf.iter().map(|&v| f64::from(v)).collect();
            let mut c_mixed = vec![0.0f64; m * n];
            gemm_mixed(m, n, k, 1.0, &af, &bf, 0.0, &mut c_mixed);
            let c_ref = naive_chunked(m, n, k, 1.0, &ad, &bd, 0.0, &vec![0.0; m * n]);
            for i in 0..m * n {
                assert_eq!(
                    c_mixed[i].to_bits(),
                    c_ref[i].to_bits(),
                    "mixed != widened f64 at {i} ({m}x{n}x{k})"
                );
            }
            // Packed f32 vs the naive-chunked f32 oracle: BITWISE (packing
            // is layout, not arithmetic — s6 contract).
            let mut c32 = vec![0.0f32; m * n];
            gemm_f32(m, n, k, 1.25, &af, &bf, 0.0, &mut c32);
            let mut c32_ref = vec![0.0f32; m * n];
            naive_chunked_f32(m, n, k, 1.25, &af, &bf, &mut c32_ref);
            for i in 0..m * n {
                assert_eq!(
                    c32[i].to_bits(),
                    c32_ref[i].to_bits(),
                    "packed f32 != naive-chunked oracle at {i} ({m}x{n}x{k})"
                );
            }
            // And the accuracy envelope vs f64 still holds.
            let mut c32a = vec![0.0f32; m * n];
            gemm_f32(m, n, k, 1.0, &af, &bf, 0.0, &mut c32a);
            for i in 0..m * n {
                let err = (f64::from(c32a[i]) - c_ref[i]).abs();
                assert!(
                    err <= 1e-4 * c_ref[i].abs().max(1.0),
                    "f32 path error {err} at {i} ({m}x{n}x{k})"
                );
            }
        }
        println!(
            "{{\"suite\":\"fs-la\",\"case\":\"gemm-precisions\",\"verdict\":\"pass\",\"detail\":\"3 shapes: mixed == widened-f64 bitwise; packed f32 == naive-chunked oracle bitwise; f32 within 1e-4 of f64\"}}"
        );
    }

    #[test]
    fn op_forms_bitwise_vs_materialized() {
        // The op-form contract: packing absorbs op()/ld, so every
        // combination is BITWISE the plain gemm_f64 on materialized
        // operands. Shapes with tails in all dims + multi-tile.
        for &(m, n, k) in &[(9usize, 5usize, 257usize), (33, 17, 300), (24, 18, 40)] {
            let a = rand_mat(m, k, 0x70);
            let b = rand_mat(k, n, 0x71);
            let c0 = rand_mat(m, n, 0x72);
            let mut want = c0.clone();
            gemm_f64(m, n, k, 1.25, &a, &b, 0.5, &mut want);
            let at: Vec<f64> = (0..k * m).map(|i| a[(i % m) * k + i / m]).collect();
            let bt: Vec<f64> = (0..n * k).map(|i| b[(i % k) * n + i / k]).collect();
            for (ta, tb) in [
                (Trans::N, Trans::N),
                (Trans::T, Trans::N),
                (Trans::N, Trans::T),
                (Trans::T, Trans::T),
            ] {
                let (av, lda) = match ta {
                    Trans::N => (&a, k),
                    Trans::T => (&at, m),
                };
                let (bv, ldb) = match tb {
                    Trans::N => (&b, n),
                    Trans::T => (&bt, k),
                };
                let mut c = c0.clone();
                gemm_f64_op(m, n, k, 1.25, av, lda, ta, bv, ldb, tb, 0.5, &mut c, n);
                for i in 0..m * n {
                    assert_eq!(
                        c[i].to_bits(),
                        want[i].to_bits(),
                        "op({ta:?},{tb:?}) != gemm_f64 at {i} ({m}x{n}x{k})"
                    );
                }
            }
        }
        println!(
            "{{\"suite\":\"fs-la\",\"case\":\"gemm-op-forms\",\"verdict\":\"pass\",\"detail\":\"NN/TN/NT/TT bitwise == materialized gemm_f64 over 3 shapes\"}}"
        );
    }

    #[test]
    fn strided_views_bitwise_and_untouched_outside() {
        // Submatrix views (ld > view cols) compute bitwise-identically
        // to contiguous copies, and C outside the view is UNTOUCHED.
        let (m, n, k) = (13usize, 11usize, 70usize);
        let (lda, ldb, ldc) = (k + 5, n + 3, n + 7);
        let mut s = 0x77_u64;
        let a_buf: Vec<f64> = (0..m * lda).map(|_| lcg(&mut s)).collect();
        let b_buf: Vec<f64> = (0..k * ldb).map(|_| lcg(&mut s)).collect();
        let c_buf0: Vec<f64> = (0..m * ldc).map(|_| lcg(&mut s)).collect();
        // Contiguous copies of the views.
        let a: Vec<f64> = (0..m * k).map(|i| a_buf[(i / k) * lda + i % k]).collect();
        let b: Vec<f64> = (0..k * n).map(|i| b_buf[(i / n) * ldb + i % n]).collect();
        let c0: Vec<f64> = (0..m * n).map(|i| c_buf0[(i / n) * ldc + i % n]).collect();
        let mut want = c0.clone();
        gemm_f64(m, n, k, -0.75, &a, &b, 0.5, &mut want);
        let mut c_buf = c_buf0.clone();
        gemm_f64_op(
            m,
            n,
            k,
            -0.75,
            &a_buf,
            lda,
            Trans::N,
            &b_buf,
            ldb,
            Trans::N,
            0.5,
            &mut c_buf,
            ldc,
        );
        for i in 0..m {
            for j in 0..ldc {
                let got = c_buf[i * ldc + j];
                if j < n {
                    assert_eq!(
                        got.to_bits(),
                        want[i * n + j].to_bits(),
                        "strided view != contiguous at ({i},{j})"
                    );
                } else {
                    assert_eq!(
                        got.to_bits(),
                        c_buf0[i * ldc + j].to_bits(),
                        "C outside the view was touched at ({i},{j})"
                    );
                }
            }
        }
        println!(
            "{{\"suite\":\"fs-la\",\"case\":\"gemm-strided\",\"verdict\":\"pass\",\"detail\":\"lda/ldb/ldc views bitwise == contiguous; outside-view C untouched\"}}"
        );
    }

    #[test]
    fn deterministic_golden_hash() {
        let (m, n, k) = (48usize, 36usize, 300usize);
        let a = rand_mat(m, k, 0x60);
        let b = rand_mat(k, n, 0x61);
        let run = || {
            let mut c = vec![0.0; m * n];
            gemm_f64(m, n, k, 1.25, &a, &b, 0.0, &mut c);
            let mut acc: u64 = 0xcbf2_9ce4_8422_2325;
            for v in &c {
                for byte in v.to_bits().to_le_bytes() {
                    acc ^= u64::from(byte);
                    acc = acc.wrapping_mul(0x0000_0100_0000_01b3);
                }
            }
            acc
        };
        let h = run();
        assert_eq!(h, run(), "same inputs must give identical bits");
        println!(
            "{{\"suite\":\"fs-la\",\"case\":\"gemm-golden\",\"verdict\":\"info\",\"detail\":\"{h:#018x}\"}}"
        );
        assert_eq!(
            h, GOLDEN_HASH,
            "GEMM output bits changed: {h:#018x} vs {GOLDEN_HASH:#018x} — KC is part of the \
             bit contract; bump only with semantic justification"
        );
    }

    /// G4: inject cancellation after a deterministic number of real packing /
    /// compute boundaries. This avoids a sleep-based race while proving that
    /// a mid-dispatch request drains with partial PRIVATE progress and never
    /// exposes a partially scaled or multiplied C.
    #[test]
    fn cancellable_gemm_mid_dispatch_drains_transactionally() {
        use std::sync::atomic::{AtomicUsize, Ordering};

        const CANCEL_AFTER: usize = 80;

        let (m, n, k) = (263usize, 37usize, 257usize);
        let a = rand_mat(m, k, 0xCA11);
        let b = rand_mat(k, n, 0xCE11);
        let original = rand_mat(m, n, 0xC0DE);
        let mut c = original.clone();
        let polls = AtomicUsize::new(0);
        let poll = || polls.fetch_add(1, Ordering::SeqCst) >= CANCEL_AFTER;
        let pool = fs_exec::TilePool::for_host(3, 0xCA11);
        let gate = fs_exec::CancelGate::new();

        let error = match gemm_f64_parallel_with_pool_and_poll(
            m,
            n,
            k,
            1.25,
            &a,
            &b,
            0.5,
            &mut c,
            &pool,
            32,
            37,
            &gate,
            fs_exec::RunId(41),
            GemmMemoryEnvelope::UNBOUNDED,
            &poll,
        ) {
            Err(GemmRunError::Cancelled(error)) => error,
            other => panic!("expected structured cancellation, got {other:?}"),
        };
        assert!(
            error.report.completed_tiles > 0,
            "injection must happen after real compute, not during setup: {error:?}"
        );
        assert!(
            error.report.completed_tiles < error.report.total_tiles,
            "injection must interrupt before finalize: {error:?}"
        );
        assert!(
            c.iter()
                .zip(&original)
                .all(|(got, before)| got.to_bits() == before.to_bits()),
            "cancelled transactional GEMM changed caller-visible C"
        );
        assert!(
            error
                .report
                .pool_runs
                .iter()
                .any(|run| run.kernel == "fs-la/gemm-f64-m-band-v1" && run.total > 0),
            "cancellation receipt must prove traversal through the real TilePool"
        );
        assert!(
            pool.arena_pool().stats().quiescent(),
            "cancelled GEMM must fully drain every Cx arena"
        );
        // Every worker exits on its first true poll; only the scope-finalizer
        // may add another observation. This guards accidental unbounded work
        // after a request as the loop nest evolves.
        assert!(
            polls.load(Ordering::SeqCst) <= CANCEL_AFTER + 2 * 3 + 4,
            "workers kept polling/working after cancellation"
        );
        assert_eq!(GEMM_MAX_FMAS_BETWEEN_POLLS, 8 * 4 * 257);
    }

    /// G0/G5: the cancellation-capable success path commits exactly once and
    /// remains bitwise the established serial accumulation contract.
    #[test]
    fn cancellable_gemm_success_is_bitwise_across_pool_placement() {
        let (m, n, k) = (263usize, 37usize, 257usize);
        let a = rand_mat(m, k, 0xA11);
        let b = rand_mat(k, n, 0xB11);
        let original = rand_mat(m, n, 0xC11);
        let mut expected = original.clone();
        gemm_f64(m, n, k, 1.25, &a, &b, 0.5, &mut expected);

        let parallelism =
            std::thread::available_parallelism().map_or(2, core::num::NonZeroUsize::get);
        let mut configs = vec![
            ("one", fs_exec::PoolConfig::for_host(1, 0xA11)),
            ("two", fs_exec::PoolConfig::for_host(2, 0xA11)),
            (
                "available-parallelism",
                fs_exec::PoolConfig::for_host(parallelism, 0xA11),
            ),
        ];
        let mut pinned = fs_exec::PoolConfig::for_host(parallelism, 0xA11);
        pinned.pin_groups = vec![vec![9999], vec![0]];
        configs.push(("advisory-pinned", pinned));

        for (placement, config) in configs {
            let pool = fs_exec::TilePool::new(config);
            let gate = fs_exec::CancelGate::new();
            let mut actual = original.clone();
            let report = gemm_f64_parallel_with_pool(
                m,
                n,
                k,
                1.25,
                &a,
                &b,
                0.5,
                &mut actual,
                &pool,
                32,
                37,
                &gate,
            )
            .expect("unrequested dispatch");
            assert_eq!(report.completed_tiles, report.total_tiles);
            assert!(report.total_tiles > 0);
            assert!(!report.pool_runs.is_empty());
            assert!(report.pool_runs.iter().all(|run| {
                run.kernel == "fs-la/gemm-f64-m-band-v1"
                    && run.completed == run.total
                    && run.total > 0
            }));
            assert_eq!(
                report
                    .pool_runs
                    .iter()
                    .map(|run| run.declared_run.0)
                    .collect::<Vec<_>>(),
                (0..u64::try_from(report.pool_runs.len()).expect("panel count fits u64"))
                    .map(|ordinal| gemm_panel_run_id(fs_exec::RunId::default(), ordinal).0)
                    .collect::<Vec<_>>(),
                "NC/KC panel stream identities must be distinct and deterministic"
            );
            assert!(
                actual
                    .iter()
                    .zip(&expected)
                    .all(|(got, want)| got.to_bits() == want.to_bits()),
                "cancellable GEMM diverged under {placement}"
            );
            assert!(pool.arena_pool().stats().quiescent(), "{placement}");
        }
    }

    /// G4: an already-requested gate is refused after shape validation and
    /// before staging allocation or mutation.
    #[test]
    fn cancellable_gemm_pre_requested_gate_leaves_c_untouched() {
        let (m, n, k) = (2usize, 3usize, 4usize);
        let a = rand_mat(m, k, 1);
        let b = rand_mat(k, n, 2);
        let original = rand_mat(m, n, 3);
        let mut c = original.clone();
        let gate = fs_exec::CancelGate::new();
        gate.request();
        let error =
            gemm_f64_parallel_with_cancel(m, n, k, 1.0, &a, &b, 0.0, &mut c, 2, 32, 128, &gate)
                .expect_err("pre-requested gate");
        let GemmRunError::Cancelled(error) = error else {
            panic!("pre-requested gate returned an executor failure");
        };
        assert_eq!(error.report.completed_tiles, 0);
        assert!(error.report.total_tiles > 0);
        assert!(
            c.iter()
                .zip(&original)
                .all(|(got, before)| got.to_bits() == before.to_bits())
        );
    }

    #[test]
    fn caller_pool_arena_refusal_is_structured_and_transactional() {
        let (m, n, k) = (16usize, 7usize, 9usize);
        let a = rand_mat(m, k, 0xA110);
        let b = rand_mat(k, n, 0xB110);
        let original = rand_mat(m, n, 0xC110);
        let mut c = original.clone();
        let mut config = fs_exec::PoolConfig::for_host(2, 0xA110);
        config.arena.limit_bytes = Some(0);
        let pool = fs_exec::TilePool::new(config);
        let gate = fs_exec::CancelGate::new();

        let error =
            gemm_f64_parallel_with_pool(m, n, k, 1.0, &a, &b, 0.5, &mut c, &pool, 8, 7, &gate)
                .expect_err("the zero-byte arena budget must refuse A packing");
        match error {
            GemmRunError::MemoryRefused {
                what,
                requested_bytes,
                report,
                ..
            } => {
                assert_eq!(what, "a-pack-arena");
                assert!(requested_bytes > 0);
                assert_eq!(report.completed_tiles, 0);
                assert_eq!(report.pool_runs.len(), 1);
                assert_eq!(report.memory.refused_bytes, requested_bytes);
                assert!(report.memory.peak_used_bytes > 0);
            }
            other => panic!("expected typed arena refusal, got {other:?}"),
        }
        assert!(
            gate.is_requested(),
            "tile containment drains through the gate"
        );
        assert!(
            c.iter()
                .zip(&original)
                .all(|(got, before)| got.to_bits() == before.to_bits()),
            "executor refusal must not expose staged output"
        );
        assert!(pool.arena_pool().stats().quiescent());
    }

    #[test]
    fn arena_refusal_after_completed_panel_preserves_progress_and_c() {
        let (m, n, k) = (8usize, 4usize, KC + 1);
        let a = rand_mat(m, k, 0xA220);
        let b = rand_mat(k, n, 0xB220);
        let original = rand_mat(m, n, 0xC220);
        let mut c = original.clone();
        let mut config = fs_exec::PoolConfig::for_host(1, 0xA220);
        config.arena.chunk_bytes = 1 << 20;
        config.arena.limit_bytes = Some(1 << 20);
        let pool = fs_exec::TilePool::new(config);
        let gate = fs_exec::CancelGate::new();
        let holder = std::sync::Mutex::new(None::<fs_alloc::Arena>);
        let poll = || {
            let stats = pool.arena_pool().stats();
            if stats.chunks_created > 0 && stats.arenas_live == 0 {
                let mut held = holder.lock().expect("holder lock");
                if held.is_none() {
                    let arena = pool.arena_pool().arena();
                    arena
                        .alloc_slice_fill(fs_alloc::Site::named("fs-la/test-held-arena"), 1, 0_u8)
                        .expect("recycle the first panel's chunk");
                    *held = Some(arena);
                }
            }
            false
        };

        let error = gemm_f64_parallel_with_pool_and_poll(
            m,
            n,
            k,
            1.0,
            &a,
            &b,
            0.5,
            &mut c,
            &pool,
            8,
            4,
            &gate,
            fs_exec::RunId(22),
            GemmMemoryEnvelope::UNBOUNDED,
            &poll,
        )
        .expect_err("the held recycled chunk must refuse the second panel arena");
        let GemmRunError::MemoryRefused {
            what,
            requested_bytes,
            report,
            ..
        } = error
        else {
            panic!("expected post-progress memory refusal, got {error:?}");
        };
        assert_eq!(what, "a-pack-arena");
        assert_eq!(report.completed_tiles, 1);
        assert_eq!(report.total_tiles, 2);
        assert_eq!(report.pool_runs.len(), 2);
        assert_eq!(report.memory.refused_bytes, requested_bytes);
        assert!(report.memory.peak_used_bytes <= report.memory.requested_bytes);
        assert!(
            c.iter()
                .zip(&original)
                .all(|(got, before)| got.to_bits() == before.to_bits())
        );
        drop(holder.lock().expect("holder lock").take());
        assert!(pool.arena_pool().stats().quiescent());
    }

    #[test]
    fn memory_planning_is_checked_and_no_product_ignores_huge_operands() {
        assert_eq!(
            checked_memory_sum("test-overflow", [u128::MAX, 1]),
            Err("test-overflow")
        );
        let no_product = preflight_memory(
            17,
            usize::MAX,
            usize::MAX,
            usize::MAX,
            false,
            8,
            4,
            usize::MAX,
            0,
            GemmMemoryEnvelope { limit_bytes: 136 },
        )
        .expect("no-product planning does not inspect unused operand extents");
        assert_eq!(no_product.report.requested_bytes, 136);
        assert_eq!(no_product.report.b_pack_bytes, 0);
        assert_eq!(no_product.report.arena_bytes, 0);

        let huge_product = preflight_memory(
            0,
            1,
            usize::MAX,
            1,
            true,
            8,
            usize::MAX,
            1,
            1,
            GemmMemoryEnvelope::UNBOUNDED,
        )
        .expect("u128 plan represents the exact oversized layout");
        assert!(huge_product.b_pack_len.is_none());
        assert!(huge_product.report.b_pack_bytes > u128::from(u64::MAX));
    }

    /// Recorded on aarch64-apple (M4 Pro); must match on x86-64 (trj).
    const GOLDEN_HASH: u64 = 0x1d7a_a3c6_b631_7ef0;
}
