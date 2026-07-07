//! fs-tilelang — the safe tile-kernel DSL runtime (plan patch Rev C).
//! Layer: L0 SUBSTRATE.
//!
//! The `kernel!` macro (re-exported from fs-tilelang-macros) lowers
//! ONE restricted kernel body into: a scalar reference variant, a
//! lane-shaped variant (chunked loops the autovectorizer maps onto the
//! resolved SIMD tier — per-element arithmetic is IDENTICAL, so the
//! variants are bitwise-equal by construction), kernel METADATA
//! (arithmetic intensity for the roofline harness and autotuner — P6:
//! every kernel ships its intensity analysis), and generated
//! G0 tier-equivalence + G5 determinism twin tests.
//!
//! This crate is the runtime the generated code targets: metadata
//! types, lane-width resolution (once, via fs-substrate dispatch —
//! never in hot loops), and deterministic/fast reduction combiners.

pub use fs_tilelang_macros::kernel;

/// Determinism class a kernel declares (part of its metadata).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DeterminismClass {
    /// Bitwise identical across tiers, batch shapes, and runs.
    BitwiseAllTiers,
    /// Deterministic per tier; reductions envelope-bounded across
    /// tiers (the fs-simd reduction class).
    PerTier,
}

/// Reduction flavor: the DETERMINISTIC variant combines fixed-width
/// chunk partials in index order (a fixed-shape tree keyed by logical
/// position — never by worker); the FAST variant may reassociate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReductionKind {
    /// No reduction output.
    None,
    /// Fixed-shape deterministic sum.
    DeterministicSum,
    /// Reassociation-permitted sum (bit-pattern NOT part of any
    /// contract; must agree with the deterministic variant within an
    /// envelope).
    FastSum,
}

/// Static per-kernel metadata emitted by the macro.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KernelMeta {
    /// Kernel name (the `kernel!` declaration name).
    pub name: &'static str,
    /// Floating-point operations per processed element (macro-time
    /// count of arithmetic operators and mul_add calls in the body;
    /// mul_add counts as 2).
    pub flops_per_elem: u32,
    /// Bytes moved per processed element (8 per f64 read/write
    /// buffer, 4 per u32 index buffer).
    pub bytes_per_elem: u32,
    /// Declared halo (elements skipped at each end for stencils).
    pub halo: u32,
    /// Reduction flavor.
    pub reduction: ReductionKind,
    /// Determinism class.
    pub determinism: DeterminismClass,
}

impl KernelMeta {
    /// Arithmetic intensity in FLOP/byte (the roofline x-axis).
    #[must_use]
    pub fn intensity(&self) -> f64 {
        f64::from(self.flops_per_elem) / f64::from(self.bytes_per_elem.max(1))
    }

    /// One JSON metadata line (roofline/autotuner food, ledger-ready).
    #[must_use]
    pub fn descr(&self) -> String {
        format!(
            "{{\"kernel\":\"{}\",\"flops_per_elem\":{},\"bytes_per_elem\":{},\
             \"intensity\":{:.4},\"halo\":{},\"reduction\":\"{:?}\",\"determinism\":\"{:?}\"}}",
            self.name,
            self.flops_per_elem,
            self.bytes_per_elem,
            self.intensity(),
            self.halo,
            self.reduction,
            self.determinism
        )
    }
}

/// Lane width for the RESOLVED SIMD tier (elements of f64 per lane
/// group): Scalar = 1, NEON = 2, AVX2 = 4, AVX-512 = 8. Resolved once
/// per call site through fs-substrate's cached dispatch — generated
/// kernels hoist this out of their loops.
#[must_use]
pub fn lane_width() -> usize {
    match fs_substrate::dispatch_tier() {
        fs_substrate::SimdTier::Scalar => 1,
        fs_substrate::SimdTier::Neon => 2,
        fs_substrate::SimdTier::Avx2 => 4,
        fs_substrate::SimdTier::Avx512 => 8,
    }
}

/// Chunk quantum for deterministic reductions: partials are formed
/// over fixed 64-element chunks and combined in index order. The
/// shape is a function of LENGTH ONLY — never of tier or thread — so
/// deterministic-sum results are bitwise identical everywhere.
pub const REDUCTION_CHUNK: usize = 64;

/// Fixed-shape deterministic sum: per-chunk sequential partials
/// combined in chunk-index order.
#[must_use]
pub fn deterministic_sum(values: &[f64]) -> f64 {
    let mut total = 0.0f64;
    for chunk in values.chunks(REDUCTION_CHUNK) {
        let mut partial = 0.0f64;
        for &v in chunk {
            partial += v;
        }
        total += partial;
    }
    total
}

/// Reassociation-permitted sum (currently a straight fold; the
/// contract permits future lane-parallel reassociation, which is why
/// its bit pattern is NOT part of any golden).
#[must_use]
pub fn fast_sum(values: &[f64]) -> f64 {
    values.iter().sum()
}

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
