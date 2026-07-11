//! The autotuner (plan §5.5): first-run and periodically refreshed
//! calibration so that NOTHING performance-critical is hardcoded to either
//! microarchitecture.
//!
//! Division of honesty: measurements are wall-clock and jittery — they live
//! in tune ROWS and calibration REPORTS (envelope-class, like fs-substrate
//! bandwidth). DECISIONS derived from them (tile edges, schedule kind) are
//! recorded values a study can PIN, and replay uses the recorded plans,
//! never re-tuned ones — that is what keeps deterministic replays faithful
//! on a different (or re-calibrated) machine.
//!
//! Persistence: a strict JSON-lines file store keyed by kernel × shape-class
//! × machine fingerprint. Evidence schema v1 distinguishes wall-time
//! samples, work counters, and scaled throughput instead of laundering every
//! integer as nanoseconds. Rows keyed to a DIFFERENT fingerprint are stale by
//! definition and ignored on load.

use crate::cx::{CancelGate, Cx};
use crate::kernel::{TileKernel, TilePlan};
use crate::pool::{PoolConfig, TilePool};
use core::fmt;
use core::ops::ControlFlow;
use fs_substrate::affinity::CcdTopology;
use fs_substrate::tile::TileEdge;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;
use std::io::Read as _;

const MAX_TUNE_STORE_BYTES: usize = 16 * 1024 * 1024;
const MAX_TUNE_ROW_BYTES: usize = 1024 * 1024;
const MAX_TUNE_STRING_BYTES: usize = 64 * 1024;
const MAX_TUNE_OBSERVATIONS: usize = 4096;
const MAX_WALL_TIME_SAMPLES: usize = 4096;
const MAX_TUNE_ROW_WALL_TIME_SAMPLES: usize = 4096;
const MAX_GEMM_IDENTITY_COMPONENT_BYTES: usize = 256;
const MAX_RETAINED_TUNING_DECISIONS: usize = 4096;
const MAX_RETAINED_TUNING_DECISION_BYTES: usize = 1024 * 1024;

/// Schedule polymorphism (plan §5.1 consequence 2): the same algorithm
/// ships a bandwidth-rich schedule (fewer, fatter, streaming-friendly
/// tiles) and a bandwidth-starved one (aggressive blocking,
/// recomputation-over-reload); the tuner selects per machine.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScheduleKind {
    /// Plenty of bandwidth per core (Apple unified-memory class).
    BandwidthRich,
    /// Bandwidth-starved cores (high-core-count x86 class).
    BandwidthStarved,
}

impl ScheduleKind {
    /// Stable name for rows/logs.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            ScheduleKind::BandwidthRich => "bandwidth-rich",
            ScheduleKind::BandwidthStarved => "bandwidth-starved",
        }
    }
}

/// Kernel-key prefix for the parallel-GEMM blocking lane. Consumers
/// append the PRODUCER's bit-semantics version (fs-la's
/// `GEMM_BIT_SEMANTICS_VERSION`) so rows measured under a different
/// accumulation contract can never match a lookup — semantic staleness
/// is filtered by key construction, not by trust.
pub const GEMM_KERNEL_PREFIX: &str = "gemm-f64-parallel/bits-v";

/// An MC/NC blocking plan for the parallel GEMM lane. Both members are
/// BIT-NEUTRAL by the producer's determinism contract (pure m/n tiling);
/// KC and the SIMD tier are part of the bit contract and stay OUTSIDE
/// this tuning loop.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GemmBlockPlan {
    /// Band height for the parallel M loop (rows per work-stealing band).
    pub mc: usize,
    /// Cap on the packed-B panel width (the engine takes `min(n, nc_cap)`).
    pub nc_cap: usize,
}

impl GemmBlockPlan {
    /// The documented cold-start default: the xlvx s5 sweep winner on
    /// both reference machines (thin mc = 32 bands, nc = n capped so the
    /// B pack stays L3-resident).
    pub const COLD_START: GemmBlockPlan = GemmBlockPlan {
        mc: 32,
        nc_cap: 2048,
    };

    /// A validated plan: `mc` a multiple of 8 in `[8, 1024]`, `nc_cap` a
    /// multiple of 128 in `[128, 8192]` — the bounded candidate lattice.
    ///
    /// # Errors
    /// [`TuneError`] outside the lattice (fail closed: an untrusted row
    /// never becomes an unbounded pack allocation).
    pub fn new(mc: usize, nc_cap: usize) -> Result<Self, TuneError> {
        if !(8..=1024).contains(&mc) || !mc.is_multiple_of(8) {
            return Err(TuneError {
                detail: format!("gemm mc {mc} outside the lattice (multiple of 8 in [8, 1024])"),
            });
        }
        if !(128..=8192).contains(&nc_cap) || !nc_cap.is_multiple_of(128) {
            return Err(TuneError {
                detail: format!(
                    "gemm nc_cap {nc_cap} outside the lattice (multiple of 128 in [128, 8192])"
                ),
            });
        }
        Ok(GemmBlockPlan { mc, nc_cap })
    }

    /// Canonical params form (`mc=32,nc-cap=2048`).
    #[must_use]
    pub fn canonical(&self) -> String {
        format!("mc={},nc-cap={}", self.mc, self.nc_cap)
    }

    /// Parse the canonical form; `None` for anything else (fail closed).
    #[must_use]
    pub fn parse(params: &str) -> Option<Self> {
        let rest = params.strip_prefix("mc=")?;
        let (mc, nc) = rest.split_once(",nc-cap=")?;
        let mc: usize = mc.parse().ok()?;
        let nc_cap: usize = nc.parse().ok()?;
        let plan = GemmBlockPlan::new(mc, nc_cap).ok()?;
        // Round-trip discipline: only canonical spellings are pinnable.
        (plan.canonical() == params).then_some(plan)
    }
}

/// Execution dimensions that can change which GEMM blocking plan wins.
///
/// Every field participates in the persistent tune key. This prevents a row
/// measured with one thread count, memory envelope, probe, ISA tier, placement
/// policy, or implementation/build from silently dispatching under another
/// configuration.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GemmExecutionIdentity {
    requested_threads: u64,
    thread_budget: u64,
    memory_limit_bytes: u64,
    probe_dims: [u64; 3],
    isa_tier: String,
    placement: String,
    implementation: String,
    build: String,
}

impl GemmExecutionIdentity {
    /// Construct a validated execution identity.
    ///
    /// `requested_threads` may be zero when the producer defines zero as an
    /// automatic/default request, but `thread_budget` (the producer-normalized
    /// maximum, not a claim about candidate-dependent spawned workers) and
    /// every probe extent must be non-zero. Text components use a deliberately
    /// narrow canonical alphabet so the scoped key has one unambiguous
    /// spelling.
    ///
    /// # Errors
    /// Returns [`TuneError`] for unrepresentable dimensions, zero effective
    /// dimensions, or non-canonical text components.
    pub fn new(
        requested_threads: usize,
        thread_budget: usize,
        memory_limit_bytes: u64,
        probe_dims: [usize; 3],
        isa_tier: impl Into<String>,
        placement: impl Into<String>,
        implementation: impl Into<String>,
        build: impl Into<String>,
    ) -> Result<Self, TuneError> {
        let requested_threads = u64::try_from(requested_threads).map_err(|_| TuneError {
            detail: "requested GEMM thread count is not representable as u64".to_string(),
        })?;
        let thread_budget = u64::try_from(thread_budget).map_err(|_| TuneError {
            detail: "normalized GEMM thread budget is not representable as u64".to_string(),
        })?;
        if thread_budget == 0 {
            return Err(TuneError {
                detail: "normalized GEMM thread budget must be non-zero".to_string(),
            });
        }
        let mut canonical_probe = [0_u64; 3];
        for (axis, (dst, extent)) in canonical_probe.iter_mut().zip(probe_dims).enumerate() {
            *dst = u64::try_from(extent).map_err(|_| TuneError {
                detail: format!("GEMM probe extent {axis} is not representable as u64"),
            })?;
            if *dst == 0 {
                return Err(TuneError {
                    detail: format!("GEMM probe extent {axis} must be non-zero"),
                });
            }
        }
        let isa_tier = isa_tier.into();
        let placement = placement.into();
        let implementation = implementation.into();
        let build = build.into();
        require_gemm_identity_component("ISA tier", &isa_tier)?;
        require_gemm_identity_component("placement", &placement)?;
        require_gemm_identity_component("implementation", &implementation)?;
        require_gemm_identity_component("build", &build)?;
        Ok(Self {
            requested_threads,
            thread_budget,
            memory_limit_bytes,
            probe_dims: canonical_probe,
            isa_tier,
            placement,
            implementation,
            build,
        })
    }

    /// Requested thread count, before producer normalization.
    #[must_use]
    pub const fn requested_threads(&self) -> u64 {
        self.requested_threads
    }

    /// Producer-normalized maximum thread budget. Actual spawned workers may
    /// be lower and candidate-dependent.
    #[must_use]
    pub const fn thread_budget(&self) -> u64 {
        self.thread_budget
    }

    /// Caller-declared memory ceiling for the measured and selected GEMM.
    /// `u64::MAX` is the explicit legacy/unbounded class, not an omitted field.
    #[must_use]
    pub const fn memory_limit_bytes(&self) -> u64 {
        self.memory_limit_bytes
    }

    /// Exact dimensions of the measured probe.
    #[must_use]
    pub const fn probe_dims(&self) -> [u64; 3] {
        self.probe_dims
    }

    /// Producer-resolved ISA tier.
    #[must_use]
    pub fn isa_tier(&self) -> &str {
        &self.isa_tier
    }

    /// Producer-resolved placement policy.
    #[must_use]
    pub fn placement(&self) -> &str {
        &self.placement
    }

    /// Producer implementation identity/version.
    #[must_use]
    pub fn implementation(&self) -> &str {
        &self.implementation
    }

    /// Producer codegen/build identity. Durable rows never cross this seam.
    #[must_use]
    pub fn build(&self) -> &str {
        &self.build
    }
}

/// Canonical GEMM tuning key. The storage kernel embeds the shape class and
/// every execution dimension, so row lookup, pin lookup, ledger lookup, and
/// recorded decisions all use the same identity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GemmTuneKey {
    base_kernel: String,
    shape_class: String,
    execution: GemmExecutionIdentity,
    scoped_kernel: String,
}

impl GemmTuneKey {
    /// Build a fully scoped key from a bit-semantics kernel, shape class, and
    /// execution identity.
    ///
    /// # Errors
    /// Returns [`TuneError`] unless `base_kernel` is exactly
    /// `gemm-f64-parallel/bits-vN` with a canonical decimal `N`, and the shape
    /// class is a canonical identity component.
    pub fn new(
        base_kernel: impl Into<String>,
        shape_class: impl Into<String>,
        execution: GemmExecutionIdentity,
    ) -> Result<Self, TuneError> {
        let base_kernel = base_kernel.into();
        if parse_gemm_base_version(&base_kernel).is_none() {
            return Err(TuneError {
                detail: format!(
                    "GEMM base kernel must be {GEMM_KERNEL_PREFIX}N with a canonical decimal version, got {base_kernel:?}"
                ),
            });
        }
        let shape_class = shape_class.into();
        require_gemm_identity_component("shape class", &shape_class)?;
        let scoped_kernel = format!(
            "{base_kernel}/tune-v3/shape={shape_class}/requested={}/thread-budget={}/memory-limit={}/probe={}x{}x{}/tier={}/placement={}/implementation={}/build={}",
            execution.requested_threads,
            execution.thread_budget,
            execution.memory_limit_bytes,
            execution.probe_dims[0],
            execution.probe_dims[1],
            execution.probe_dims[2],
            execution.isa_tier,
            execution.placement,
            execution.implementation,
            execution.build,
        );
        if scoped_kernel.len() > MAX_TUNE_STRING_BYTES {
            return Err(TuneError {
                detail: format!(
                    "scoped GEMM kernel exceeds the {MAX_TUNE_STRING_BYTES}-byte limit"
                ),
            });
        }
        debug_assert_eq!(
            gemm_shape_from_scoped_kernel(&scoped_kernel),
            Some(shape_class.as_str())
        );
        Ok(Self {
            base_kernel,
            shape_class,
            execution,
            scoped_kernel,
        })
    }

    /// Fully scoped kernel used by rows, pins, ledger entries, and decisions.
    #[must_use]
    pub fn kernel(&self) -> &str {
        &self.scoped_kernel
    }

    /// Producer bit-semantics kernel before execution scoping.
    #[must_use]
    pub fn base_kernel(&self) -> &str {
        &self.base_kernel
    }

    /// Shape class redundantly bound in both the scoped kernel and row field.
    #[must_use]
    pub fn shape_class(&self) -> &str {
        &self.shape_class
    }

    /// Execution identity carried by this key.
    #[must_use]
    pub const fn execution(&self) -> &GemmExecutionIdentity {
        &self.execution
    }
}

fn require_gemm_identity_component(label: &str, value: &str) -> Result<(), TuneError> {
    if value.is_empty()
        || value.len() > MAX_GEMM_IDENTITY_COMPONENT_BYTES
        || !value
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'-' | b'.' | b'_'))
    {
        return Err(TuneError {
            detail: format!(
                "GEMM {label} must be 1..={MAX_GEMM_IDENTITY_COMPONENT_BYTES} ASCII [A-Za-z0-9._-] bytes, got {value:?}"
            ),
        });
    }
    Ok(())
}

fn require_decision_kernel(kernel: &str) -> Result<(), TuneError> {
    if kernel.trim().is_empty() {
        return Err(TuneError {
            detail: "a tuning-decision kernel identity must be non-blank".to_string(),
        });
    }
    if kernel.len() > MAX_TUNE_STRING_BYTES {
        return Err(TuneError {
            detail: format!(
                "tuning-decision kernel identity exceeds the {MAX_TUNE_STRING_BYTES}-byte limit"
            ),
        });
    }
    Ok(())
}

fn parse_canonical_u64(text: &str) -> Option<u64> {
    let value = text.parse::<u64>().ok()?;
    (value.to_string() == text).then_some(value)
}

fn parse_gemm_base_version(kernel: &str) -> Option<u64> {
    parse_canonical_u64(kernel.strip_prefix(GEMM_KERNEL_PREFIX)?)
}

fn gemm_shape_from_scoped_kernel(kernel: &str) -> Option<&str> {
    let (base, scope) = kernel.split_once("/tune-v3/")?;
    parse_gemm_base_version(base)?;
    let mut parts = scope.split('/');
    let shape = parts.next()?.strip_prefix("shape=")?;
    let requested = parts.next()?.strip_prefix("requested=")?;
    let thread_budget = parts.next()?.strip_prefix("thread-budget=")?;
    let memory_limit = parts.next()?.strip_prefix("memory-limit=")?;
    let probe = parts.next()?.strip_prefix("probe=")?;
    let tier = parts.next()?.strip_prefix("tier=")?;
    let placement = parts.next()?.strip_prefix("placement=")?;
    let implementation = parts.next()?.strip_prefix("implementation=")?;
    let build = parts.next()?.strip_prefix("build=")?;
    if parts.next().is_some()
        || parse_canonical_u64(requested).is_none()
        || parse_canonical_u64(thread_budget)? == 0
        || parse_canonical_u64(memory_limit).is_none()
        || require_gemm_identity_component("shape class", shape).is_err()
        || require_gemm_identity_component("ISA tier", tier).is_err()
        || require_gemm_identity_component("placement", placement).is_err()
        || require_gemm_identity_component("implementation", implementation).is_err()
        || require_gemm_identity_component("build", build).is_err()
    {
        return None;
    }
    let mut probe_parts = probe.split('x');
    for _ in 0..3 {
        if parse_canonical_u64(probe_parts.next()?)? == 0 {
            return None;
        }
    }
    probe_parts.next().is_none().then_some(shape)
}

/// Where a decision came from — part of every recorded decision.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TuneSource {
    /// A pinned (replayed) decision — highest precedence.
    Pinned,
    /// Derived from this machine's calibration rows.
    Tuned,
    /// No applicable row: the documented cold-start default.
    ColdStart,
}

impl TuneSource {
    fn name(self) -> &'static str {
        match self {
            TuneSource::Pinned => "pinned",
            TuneSource::Tuned => "tuned",
            TuneSource::ColdStart => "cold-start",
        }
    }
}

/// Summary derived from all wall-time samples for one candidate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WallTimeSummary {
    /// Smallest observed wall time.
    pub minimum_ns: u64,
    /// Largest observed wall time.
    pub maximum_ns: u64,
}

/// Unit carried by a work-counter observation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum WorkUnit {
    /// Tiles completed by one worker during a probe.
    CompletedTiles,
    /// Steal operations reported by the pool during a probe.
    Steals,
}

impl WorkUnit {
    const fn name(self) -> &'static str {
        match self {
            Self::CompletedTiles => "completed-tiles",
            Self::Steals => "steals",
        }
    }

    fn parse(name: &str) -> Option<Self> {
        match name {
            "completed-tiles" => Some(Self::CompletedTiles),
            "steals" => Some(Self::Steals),
            _ => None,
        }
    }
}

/// Unit carried by a throughput observation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ThroughputUnit {
    /// Decimal gigabytes per second, stored in thousandths of the unit.
    GigabytesPerSecond,
}

impl ThroughputUnit {
    const fn name(self) -> &'static str {
        match self {
            Self::GigabytesPerSecond => "gigabytes-per-second",
        }
    }

    fn parse(name: &str) -> Option<Self> {
        match name {
            "gigabytes-per-second" => Some(Self::GigabytesPerSecond),
            _ => None,
        }
    }
}

/// One typed tuning observation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TuneObservation {
    /// Exact wall-time samples and their revalidated extrema.
    WallTime {
        /// Candidate or run label.
        label: String,
        /// Every observed duration, without floating-point conversion.
        samples_ns: Vec<u64>,
        /// Summary derived from `samples_ns`.
        summary: WallTimeSummary,
    },
    /// A dimensioned integral work counter.
    WorkCount {
        /// Worker or counter label.
        label: String,
        /// Counter unit.
        unit: WorkUnit,
        /// Exact counter value.
        count: u64,
    },
    /// A dimensioned throughput measurement.
    Throughput {
        /// Throughput axis label.
        label: String,
        /// Throughput unit.
        unit: ThroughputUnit,
        /// Thousandths of `unit`, stored exactly.
        milli_units: u64,
    },
}

impl TuneObservation {
    /// Construct a wall-time observation and derive its summary.
    ///
    /// # Errors
    /// Returns [`TuneError`] for a blank label, an empty sample set, or any
    /// zero-duration sample.
    pub fn wall_time(label: impl Into<String>, samples_ns: Vec<u64>) -> Result<Self, TuneError> {
        let label = label.into();
        require_observation_label(&label)?;
        if samples_ns.len() > MAX_WALL_TIME_SAMPLES {
            return Err(TuneError {
                detail: format!(
                    "wall-time observation {label:?} exceeds the {MAX_WALL_TIME_SAMPLES}-sample limit"
                ),
            });
        }
        let minimum_ns = samples_ns.iter().copied().min().ok_or_else(|| TuneError {
            detail: format!("wall-time observation {label:?} has no samples"),
        })?;
        if minimum_ns == 0 {
            return Err(TuneError {
                detail: format!("wall-time observation {label:?} contains a zero-duration sample"),
            });
        }
        let maximum_ns = samples_ns
            .iter()
            .copied()
            .max()
            .expect("non-empty samples have a maximum");
        Ok(Self::WallTime {
            label,
            samples_ns,
            summary: WallTimeSummary {
                minimum_ns,
                maximum_ns,
            },
        })
    }

    /// Construct a dimensioned work-counter observation.
    ///
    /// # Errors
    /// Returns [`TuneError`] for a blank label.
    pub fn work_count(
        label: impl Into<String>,
        unit: WorkUnit,
        count: u64,
    ) -> Result<Self, TuneError> {
        let label = label.into();
        require_observation_label(&label)?;
        Ok(Self::WorkCount { label, unit, count })
    }

    /// Construct a dimensioned, exactly scaled throughput observation.
    ///
    /// # Errors
    /// Returns [`TuneError`] for a blank label.
    pub fn throughput(
        label: impl Into<String>,
        unit: ThroughputUnit,
        milli_units: u64,
    ) -> Result<Self, TuneError> {
        let label = label.into();
        require_observation_label(&label)?;
        Ok(Self::Throughput {
            label,
            unit,
            milli_units,
        })
    }

    fn label(&self) -> &str {
        match self {
            Self::WallTime { label, .. }
            | Self::WorkCount { label, .. }
            | Self::Throughput { label, .. } => label,
        }
    }

    fn wall_minimum(&self) -> Option<u64> {
        match self {
            Self::WallTime { summary, .. } => Some(summary.minimum_ns),
            Self::WorkCount { .. } | Self::Throughput { .. } => None,
        }
    }

    fn validate(&self) -> Result<(), TuneError> {
        require_observation_label(self.label())?;
        if let Self::WallTime {
            label,
            samples_ns,
            summary,
        } = self
        {
            let expected = Self::wall_time(label.clone(), samples_ns.clone())?;
            if &expected != self {
                return Err(TuneError {
                    detail: format!(
                        "wall-time summary for {label:?} does not match its exact samples"
                    ),
                });
            }
            if summary.minimum_ns > summary.maximum_ns {
                return Err(TuneError {
                    detail: format!("wall-time extrema for {label:?} are reversed"),
                });
            }
        }
        Ok(())
    }
}

fn require_observation_label(label: &str) -> Result<(), TuneError> {
    if label.trim().is_empty() {
        return Err(TuneError {
            detail: "tune observation labels must be non-blank".to_string(),
        });
    }
    if label.len() > MAX_TUNE_STRING_BYTES {
        return Err(TuneError {
            detail: format!(
                "tune observation label exceeds the {MAX_TUNE_STRING_BYTES}-byte limit"
            ),
        });
    }
    Ok(())
}

/// Versioned evidence carried by one tune row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TuneEvidence {
    observations: Vec<TuneObservation>,
    candidate_separation_ppm: Option<u32>,
}

impl TuneEvidence {
    /// Canonical serialized evidence schema.
    pub const VERSION: u32 = 1;

    /// Validate observations without asserting that they are ranked candidates.
    ///
    /// # Errors
    /// Returns [`TuneError`] for no observations, invalid observation
    /// summaries, blank or duplicate labels, or an aggregate wall-time sample
    /// count beyond the per-row resource budget.
    pub fn new(observations: Vec<TuneObservation>) -> Result<Self, TuneError> {
        validate_observations(&observations)?;
        Ok(Self {
            observations,
            candidate_separation_ppm: None,
        })
    }

    /// Validate a ranked wall-time candidate set and derive its separation.
    ///
    /// This explicit constructor is the only way to assert that multiple
    /// timing observations are comparable candidates. The separation is the
    /// relative fastest-to-runner-up gap, not a confidence or uncertainty
    /// estimate.
    ///
    /// # Errors
    /// Returns [`TuneError`] unless at least two valid observations are
    /// present, every observation is wall time, and their aggregate sample
    /// count is within the per-row resource budget.
    pub fn ranked_wall_times(observations: Vec<TuneObservation>) -> Result<Self, TuneError> {
        validate_observations(&observations)?;
        let candidate_separation_ppm =
            candidate_separation_ppm(&observations).ok_or_else(|| TuneError {
                detail: "ranked timing evidence requires at least two wall-time candidates"
                    .to_string(),
            })?;
        Ok(Self {
            observations,
            candidate_separation_ppm: Some(candidate_separation_ppm),
        })
    }

    /// Typed observations in canonical insertion order.
    #[must_use]
    pub fn observations(&self) -> &[TuneObservation] {
        &self.observations
    }

    /// Descriptive fastest-to-runner-up wall-time gap in parts per million.
    ///
    /// `None` means the row makes no ranked-candidate claim. The value is
    /// never a statistical confidence claim.
    #[must_use]
    pub const fn candidate_separation_ppm(&self) -> Option<u32> {
        self.candidate_separation_ppm
    }

    fn validate(&self) -> Result<(), TuneError> {
        validate_observations(&self.observations)?;
        if let Some(stored) = self.candidate_separation_ppm {
            let expected =
                candidate_separation_ppm(&self.observations).ok_or_else(|| TuneError {
                    detail: "candidate separation requires comparable wall-time observations"
                        .to_string(),
                })?;
            if expected != stored {
                return Err(TuneError {
                    detail: "candidate separation does not match wall-time observations"
                        .to_string(),
                });
            }
        }
        Ok(())
    }

    fn write_json<W: fmt::Write>(&self, out: &mut W) -> fmt::Result {
        write!(out, "{{\"version\":{},\"observations\":[", Self::VERSION)?;
        for (index, observation) in self.observations.iter().enumerate() {
            if index > 0 {
                out.write_char(',')?;
            }
            match observation {
                TuneObservation::WallTime {
                    label,
                    samples_ns,
                    summary,
                } => {
                    out.write_str("{\"kind\":\"wall-time\",\"label\":")?;
                    write_json_string(out, label)?;
                    out.write_str(",\"samples_ns\":[")?;
                    for (sample_index, sample) in samples_ns.iter().enumerate() {
                        if sample_index > 0 {
                            out.write_char(',')?;
                        }
                        write!(out, "{sample}")?;
                    }
                    write!(
                        out,
                        "],\"summary\":{{\"minimum_ns\":{},\"maximum_ns\":{}}}}}",
                        summary.minimum_ns, summary.maximum_ns
                    )?;
                }
                TuneObservation::WorkCount { label, unit, count } => {
                    out.write_str("{\"kind\":\"work-count\",\"label\":")?;
                    write_json_string(out, label)?;
                    out.write_str(",\"unit\":")?;
                    write_json_string(out, unit.name())?;
                    write!(out, ",\"count\":{count}}}")?;
                }
                TuneObservation::Throughput {
                    label,
                    unit,
                    milli_units,
                } => {
                    out.write_str("{\"kind\":\"throughput\",\"label\":")?;
                    write_json_string(out, label)?;
                    out.write_str(",\"unit\":")?;
                    write_json_string(out, unit.name())?;
                    write!(out, ",\"milli_units\":{milli_units}}}")?;
                }
            }
        }
        out.write_str("],\"candidate_separation_ppm\":")?;
        if let Some(separation) = self.candidate_separation_ppm {
            write!(out, "{separation}")?;
        } else {
            out.write_str("null")?;
        }
        out.write_char('}')
    }
}

fn validate_observations(observations: &[TuneObservation]) -> Result<(), TuneError> {
    if observations.is_empty() {
        return Err(TuneError {
            detail: "tune evidence must contain at least one observation".to_string(),
        });
    }
    if observations.len() > MAX_TUNE_OBSERVATIONS {
        return Err(TuneError {
            detail: format!("tune evidence exceeds the {MAX_TUNE_OBSERVATIONS}-observation limit"),
        });
    }
    let mut labels = BTreeSet::new();
    let mut total_wall_time_samples = 0_usize;
    for observation in observations {
        observation.validate()?;
        if let TuneObservation::WallTime { samples_ns, .. } = observation {
            total_wall_time_samples = total_wall_time_samples
                .checked_add(samples_ns.len())
                .filter(|&total| total <= MAX_TUNE_ROW_WALL_TIME_SAMPLES)
                .ok_or_else(|| TuneError {
                    detail: format!(
                        "tune evidence exceeds the {MAX_TUNE_ROW_WALL_TIME_SAMPLES}-sample aggregate wall-time limit"
                    ),
                })?;
        }
        if !labels.insert(observation.label()) {
            return Err(TuneError {
                detail: format!("duplicate tune observation label {:?}", observation.label()),
            });
        }
    }
    Ok(())
}

fn candidate_separation_ppm(observations: &[TuneObservation]) -> Option<u32> {
    if observations.len() < 2
        || observations
            .iter()
            .any(|item| item.wall_minimum().is_none())
    {
        return None;
    }
    let mut minima: Vec<(usize, u64)> = observations
        .iter()
        .enumerate()
        .map(|(index, observation)| {
            (
                index,
                observation
                    .wall_minimum()
                    .expect("all observations were checked as wall time"),
            )
        })
        .collect();
    minima.sort_unstable_by_key(|&(index, ns)| (ns, index));
    let best = minima[0].1;
    let runner_up = minima[1].1;
    if runner_up == 0 {
        return Some(0);
    }
    let delta = u128::from(runner_up - best);
    let ppm = delta * 1_000_000 / u128::from(runner_up);
    Some(u32::try_from(ppm).expect("a relative separation is at most one million ppm"))
}

fn gemm_evidence_argmin(evidence: &TuneEvidence) -> Result<GemmBlockPlan, TuneError> {
    if evidence.candidate_separation_ppm().is_none() {
        return Err(TuneError {
            detail: "a GEMM selection row requires ranked wall-time candidate evidence".to_string(),
        });
    }
    let mut winner: Option<(u64, usize, GemmBlockPlan)> = None;
    for (index, observation) in evidence.observations().iter().enumerate() {
        let plan = GemmBlockPlan::parse(observation.label()).ok_or_else(|| TuneError {
            detail: format!(
                "GEMM candidate label {:?} is not a canonical blocking plan",
                observation.label()
            ),
        })?;
        let minimum_ns = observation.wall_minimum().ok_or_else(|| TuneError {
            detail: "GEMM ranked evidence must contain only wall-time candidates".to_string(),
        })?;
        let candidate = (minimum_ns, index, plan);
        if winner
            .as_ref()
            .is_none_or(|&(best_ns, best_index, _)| (minimum_ns, index) < (best_ns, best_index))
        {
            winner = Some(candidate);
        }
    }
    winner.map(|(_, _, plan)| plan).ok_or_else(|| TuneError {
        detail: "GEMM ranked evidence contains no candidates".to_string(),
    })
}

fn validate_gemm_selection(
    selected: GemmBlockPlan,
    evidence: &TuneEvidence,
) -> Result<(), TuneError> {
    let argmin = gemm_evidence_argmin(evidence)?;
    if selected != argmin {
        return Err(TuneError {
            detail: format!(
                "selected GEMM plan {} is not evidence argmin {}",
                selected.canonical(),
                argmin.canonical()
            ),
        });
    }
    Ok(())
}

fn validate_scoped_gemm_row(row: &TuneRow) -> Result<GemmBlockPlan, TuneError> {
    let embedded_shape = gemm_shape_from_scoped_kernel(&row.kernel).ok_or_else(|| TuneError {
        detail: format!(
            "GEMM row kernel {:?} is not a canonical scoped key",
            row.kernel
        ),
    })?;
    if embedded_shape != row.shape_class {
        return Err(TuneError {
            detail: format!(
                "GEMM row shape {:?} disagrees with scoped kernel shape {embedded_shape:?}",
                row.shape_class
            ),
        });
    }
    let selected = GemmBlockPlan::parse(&row.params).ok_or_else(|| TuneError {
        detail: format!(
            "GEMM row params {:?} are not a canonical blocking plan",
            row.params
        ),
    })?;
    validate_gemm_selection(selected, &row.evidence)?;
    Ok(selected)
}

/// One tune-table row (strict file-store edition).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TuneRow {
    /// Kernel identity (e.g. "stencil7-f32").
    pub kernel: String,
    /// Shape class (e.g. "48c-cube").
    pub shape_class: String,
    /// Machine fingerprint the measurements belong to.
    pub machine: u64,
    /// Chosen parameter, canonical form (e.g. "edge=8", "schedule=bandwidth-rich").
    pub params: String,
    /// Typed, versioned measurement evidence.
    pub evidence: TuneEvidence,
    /// Recalibration counter (idempotence witness).
    pub refresh: u32,
}

impl TuneRow {
    fn write_json<W: fmt::Write>(&self, out: &mut W) -> fmt::Result {
        out.write_str("{\"kernel\":")?;
        write_json_string(out, &self.kernel)?;
        out.write_str(",\"shape_class\":")?;
        write_json_string(out, &self.shape_class)?;
        write!(out, ",\"machine\":\"{:016x}\",\"params\":", self.machine)?;
        write_json_string(out, &self.params)?;
        out.write_str(",\"evidence\":")?;
        self.evidence.write_json(out)?;
        write!(out, ",\"refresh\":{}}}", self.refresh)
    }

    /// Canonical JSON-line (deterministic field order).
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut s = String::with_capacity(160);
        self.write_json(&mut s).expect("String writes cannot fail");
        s
    }

    fn validated_generated_json(&self) -> Result<String, TuneError> {
        self.evidence.validate()?;
        let mut bounded = BoundedJsonWriter::new(MAX_TUNE_ROW_BYTES);
        self.write_json(&mut bounded).map_err(|_| TuneError {
            detail: format!(
                "generated tune row exceeds the {MAX_TUNE_ROW_BYTES}-byte canonical row limit"
            ),
        })?;
        let json = bounded.finish();
        let reparsed = parse_row(&json).ok_or_else(|| TuneError {
            detail: "generated tune row is outside the canonical parser domain".to_string(),
        })?;
        if reparsed != *self {
            return Err(TuneError {
                detail: "generated tune row does not preserve its fields through canonical parsing"
                    .to_string(),
            });
        }
        Ok(json)
    }
}

/// A validated GEMM row awaiting failure-atomic installation.
///
/// The fields are private so callers cannot change the key, machine,
/// selection, evidence, or refresh witness after validation. A session can
/// persist [`PreparedGemmRow::row_json`] and [`PreparedGemmRow::params_json`]
/// before handing the value to [`Tuner::commit_gemm_row`].
#[derive(Debug)]
pub struct PreparedGemmRow {
    key: GemmTuneKey,
    row: TuneRow,
    expected_refresh: Option<u32>,
}

/// A resolved GEMM plan whose decision receipt has not yet been recorded.
/// Sessions commit it only after the cancellable producer reports success.
#[derive(Debug, Clone)]
pub struct PreparedGemmDecision {
    key: GemmTuneKey,
    plan: GemmBlockPlan,
    source: TuneSource,
}

impl PreparedGemmDecision {
    /// Resolved blocking plan.
    #[must_use]
    pub const fn plan(&self) -> GemmBlockPlan {
        self.plan
    }

    /// Provenance of the resolved plan.
    #[must_use]
    pub const fn source(&self) -> TuneSource {
        self.source
    }

    /// Exact scoped key that was resolved.
    #[must_use]
    pub const fn key(&self) -> &GemmTuneKey {
        &self.key
    }
}

impl PreparedGemmRow {
    /// Canonical selected-plan JSON for the ledger's separate `params` field.
    #[must_use]
    pub fn params_json(&self) -> String {
        let mut out = String::new();
        write_json_string(&mut out, &self.row.params).expect("String writes cannot fail");
        out
    }

    /// Canonical tune-row JSON for the ledger's measured body.
    #[must_use]
    pub fn row_json(&self) -> String {
        self.row.to_json()
    }

    /// Fully scoped key this row was validated against.
    #[must_use]
    pub const fn key(&self) -> &GemmTuneKey {
        &self.key
    }
}

struct BoundedJsonWriter {
    out: String,
    limit: usize,
}

impl BoundedJsonWriter {
    fn new(limit: usize) -> Self {
        Self {
            out: String::with_capacity(limit.min(1024)),
            limit,
        }
    }

    fn finish(self) -> String {
        self.out
    }
}

impl fmt::Write for BoundedJsonWriter {
    fn write_str(&mut self, value: &str) -> fmt::Result {
        let Some(new_len) = self.out.len().checked_add(value.len()) else {
            return Err(fmt::Error);
        };
        if new_len > self.limit {
            return Err(fmt::Error);
        }
        self.out.push_str(value);
        Ok(())
    }
}

fn write_json_string<W: fmt::Write>(out: &mut W, value: &str) -> fmt::Result {
    out.write_char('"')?;
    for ch in value.chars() {
        match ch {
            '"' => out.write_str("\\\"")?,
            '\\' => out.write_str("\\\\")?,
            '\n' => out.write_str("\\n")?,
            '\r' => out.write_str("\\r")?,
            '\t' => out.write_str("\\t")?,
            c if c.is_control() => {
                write!(out, "\\u{:04x}", c as u32)?;
            }
            c => out.write_char(c)?,
        }
    }
    out.write_char('"')
}

/// Structured tune-store failure (Decalogue P10).
#[derive(Debug)]
pub struct TuneError {
    /// What went wrong, with the path.
    pub detail: String,
}

impl fmt::Display for TuneError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "tuner error: {}", self.detail)
    }
}

impl core::error::Error for TuneError {}

/// A recorded decision (what a study pins for replay fidelity).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TuningDecision {
    /// Kernel the decision applies to.
    pub kernel: String,
    /// Chosen parameter, canonical form.
    pub params: String,
    /// Provenance of the choice.
    pub source: &'static str,
}

impl TuningDecision {
    /// Canonical JSON object.
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut out = String::from("{\"kernel\":");
        write_json_string(&mut out, &self.kernel).expect("String writes cannot fail");
        out.push_str(",\"params\":");
        write_json_string(&mut out, &self.params).expect("String writes cannot fail");
        out.push_str(",\"source\":");
        write_json_string(&mut out, self.source).expect("String writes cannot fail");
        out.push('}');
        out
    }

    fn retained_bytes(&self) -> Option<usize> {
        self.kernel.len().checked_add(self.params.len())
    }
}

/// Borrowed metadata for the tuner's bounded in-memory decision window.
///
/// [`TuningDecisionHistory::is_complete`] is false after any oldest-prefix
/// eviction. Callers must not treat an incomplete window as the full replay
/// record; production dispatch receipts belong in the Design Ledger.
#[derive(Debug, Clone, Copy)]
pub struct TuningDecisionHistory<'a> {
    decisions: &'a [TuningDecision],
    evicted: usize,
    retained_bytes: usize,
}

impl<'a> TuningDecisionHistory<'a> {
    /// Decisions still retained, in original dispatch order.
    #[must_use]
    pub const fn decisions(self) -> &'a [TuningDecision] {
        self.decisions
    }

    /// Number of older decisions evicted from this tuner.
    #[must_use]
    pub const fn evicted(self) -> usize {
        self.evicted
    }

    /// Owned string payload retained by the current window.
    #[must_use]
    pub const fn retained_bytes(self) -> usize {
        self.retained_bytes
    }

    /// True only when the retained slice still starts at the tuner's first
    /// recorded decision.
    #[must_use]
    pub const fn is_complete(self) -> bool {
        self.evicted == 0
    }

    /// Total decisions recorded by this tuner, saturating only after an
    /// unrepresentable process-lifetime count.
    #[must_use]
    pub fn total_recorded(self) -> usize {
        self.evicted.saturating_add(self.decisions.len())
    }
}

/// The per-machine tuner: tune rows + pins + decision log.
#[derive(Debug)]
pub struct Tuner {
    fingerprint: u64,
    rows: BTreeMap<(String, String), TuneRow>,
    pins: BTreeMap<String, PinnedParam>,
    decisions: Vec<TuningDecision>,
    decision_bytes: usize,
    evicted_decisions: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PinnedParam {
    TileEdge(TileEdge),
    Schedule(ScheduleKind),
    GemmBlocking(GemmBlockPlan),
}

impl PinnedParam {
    fn parse(kernel: &str, params: &str) -> Option<Self> {
        if kernel == SCHEDULE_KERNEL {
            return match params {
                "schedule=bandwidth-rich" => Some(Self::Schedule(ScheduleKind::BandwidthRich)),
                "schedule=bandwidth-starved" => {
                    Some(Self::Schedule(ScheduleKind::BandwidthStarved))
                }
                _ => None,
            };
        }
        if kernel.starts_with(GEMM_KERNEL_PREFIX) {
            return gemm_shape_from_scoped_kernel(kernel)
                .and_then(|_| GemmBlockPlan::parse(params))
                .map(Self::GemmBlocking);
        }
        match params {
            "edge=4" => Some(Self::TileEdge(TileEdge::E4)),
            "edge=8" => Some(Self::TileEdge(TileEdge::E8)),
            "edge=16" => Some(Self::TileEdge(TileEdge::E16)),
            _ => None,
        }
    }

    fn canonical(self) -> String {
        match self {
            Self::TileEdge(edge) => format!("edge={}", edge.cells()),
            Self::Schedule(kind) => format!("schedule={}", kind.name()),
            Self::GemmBlocking(plan) => plan.canonical(),
        }
    }
}

/// The reference stencil kernel used for tile-edge calibration (a real
/// tiled workload through the real pool, not a synthetic loop).
struct CalStencil {
    field: fs_substrate::field::TiledField<f32>,
}

impl TileKernel for CalStencil {
    type Out = f64;

    fn tiles(&self) -> TilePlan {
        TilePlan::new("tune/stencil7", self.field.grid().tile_count())
    }

    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<crate::Cancelled, f64> {
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(crate::Cancelled);
        }
        // Visit this tile's cells in z-order rank order; 7-point Laplacian
        // with clamped ghosts via the boundary API.
        let zorder = self.field.grid().iter_zorder();
        let t = zorder[tile as usize];
        let e = i64::from(self.field.grid().edge().cells());
        let base = [i64::from(t.x) * e, i64::from(t.y) * e, i64::from(t.z) * e];
        let bc = fs_substrate::field::Boundary::Clamp;
        let mut acc = 0.0f64;
        for lz in 0..e {
            for ly in 0..e {
                for lx in 0..e {
                    let c = [base[0] + lx, base[1] + ly, base[2] + lz];
                    let mid = f64::from(self.field.get_bc(c, bc));
                    let l = f64::from(self.field.get_bc([c[0] - 1, c[1], c[2]], bc))
                        + f64::from(self.field.get_bc([c[0] + 1, c[1], c[2]], bc))
                        + f64::from(self.field.get_bc([c[0], c[1] - 1, c[2]], bc))
                        + f64::from(self.field.get_bc([c[0], c[1] + 1, c[2]], bc))
                        + f64::from(self.field.get_bc([c[0], c[1], c[2] - 1], bc))
                        + f64::from(self.field.get_bc([c[0], c[1], c[2] + 1], bc))
                        - 6.0 * mid;
                    acc += l;
                }
            }
        }
        ControlFlow::Continue(acc)
    }
}

const STENCIL_KERNEL: &str = "stencil7-f32";
const SCHEDULE_KERNEL: &str = "schedule";
const CLASS_WORK_KERNEL: &str = "class-work-counts";
const SHAPE_DEFAULT: &str = "48c-cube";

impl Tuner {
    /// Cold tuner for a machine (no rows; defaults answer everything).
    #[must_use]
    pub fn cold(fingerprint: u64) -> Self {
        Tuner {
            fingerprint,
            rows: BTreeMap::new(),
            pins: BTreeMap::new(),
            decisions: Vec::new(),
            decision_bytes: 0,
            evicted_decisions: 0,
        }
    }

    /// True when no calibration rows exist for this machine.
    #[must_use]
    pub fn needs_calibration(&self) -> bool {
        self.rows.is_empty()
    }

    /// Pin a typed tile edge for a kernel.
    ///
    /// # Errors
    /// Returns [`TuneError`] when `kernel` is the reserved schedule key.
    pub fn pin_tile_edge(
        &mut self,
        kernel: impl Into<String>,
        edge: TileEdge,
    ) -> Result<(), TuneError> {
        let kernel = kernel.into();
        require_decision_kernel(&kernel)?;
        if kernel == SCHEDULE_KERNEL || kernel.starts_with(GEMM_KERNEL_PREFIX) {
            return Err(TuneError {
                detail: if kernel == SCHEDULE_KERNEL {
                    "the reserved schedule key requires pin_schedule".to_string()
                } else {
                    "a GEMM-scoped key requires pin_gemm_blocking".to_string()
                },
            });
        }
        self.pins.insert(kernel, PinnedParam::TileEdge(edge));
        Ok(())
    }

    /// Pin the typed schedule decision.
    pub fn pin_schedule(&mut self, kind: ScheduleKind) {
        self.pins
            .insert(SCHEDULE_KERNEL.to_string(), PinnedParam::Schedule(kind));
    }

    /// The machine fingerprint this tuner's rows belong to.
    #[must_use]
    pub fn machine(&self) -> u64 {
        self.fingerprint
    }

    /// True when a pin exists for `kernel` (any param kind).
    #[must_use]
    pub fn has_pin(&self, kernel: &str) -> bool {
        self.pins.contains_key(kernel)
    }

    /// True when a tuned row exists for `kernel` × `shape_class`.
    #[must_use]
    pub fn has_row(&self, kernel: &str, shape_class: &str) -> bool {
        self.rows
            .contains_key(&(kernel.to_string(), shape_class.to_string()))
    }

    /// True when a typed GEMM pin exists for this exact shape and execution
    /// identity.
    #[must_use]
    pub fn has_gemm_pin(&self, key: &GemmTuneKey) -> bool {
        matches!(
            self.pins.get(key.kernel()),
            Some(PinnedParam::GemmBlocking(_))
        )
    }

    /// True when a validated GEMM row exists for this exact shape and
    /// execution identity.
    #[must_use]
    pub fn has_gemm_row(&self, key: &GemmTuneKey) -> bool {
        self.rows
            .contains_key(&(key.kernel().to_string(), key.shape_class().to_string()))
    }

    /// Pin a typed GEMM blocking plan for one exact scoped key.
    ///
    /// # Errors
    /// [`TuneError`] if the typed key's internal canonical binding is invalid.
    pub fn pin_gemm_blocking(
        &mut self,
        key: &GemmTuneKey,
        plan: GemmBlockPlan,
    ) -> Result<(), TuneError> {
        if gemm_shape_from_scoped_kernel(key.kernel()) != Some(key.shape_class()) {
            return Err(TuneError {
                detail: "typed GEMM key has an invalid internal shape binding".to_string(),
            });
        }
        self.pins
            .insert(key.kernel().to_string(), PinnedParam::GemmBlocking(plan));
        Ok(())
    }

    /// Resolve a GEMM plan without recording that it executed. Pins beat
    /// tuned rows, which beat [`GemmBlockPlan::COLD_START`]. The caller must
    /// explicitly commit the returned decision after successful dispatch.
    #[must_use]
    pub fn prepare_gemm_decision(&self, key: &GemmTuneKey) -> PreparedGemmDecision {
        let (plan, source) = self.resolve_gemm(key);
        PreparedGemmDecision {
            key: key.clone(),
            plan,
            source,
        }
    }

    /// Record a prepared GEMM decision after successful dispatch.
    ///
    /// # Errors
    /// [`TuneError`] if the applicable pin/row changed after preparation.
    pub fn commit_gemm_decision(
        &mut self,
        prepared: PreparedGemmDecision,
    ) -> Result<(), TuneError> {
        let PreparedGemmDecision { key, plan, source } = prepared;
        if self.resolve_gemm(&key) != (plan, source) {
            return Err(TuneError {
                detail: "prepared GEMM decision is stale because tuner state changed before commit"
                    .to_string(),
            });
        }
        self.record_decision(TuningDecision {
            kernel: key.kernel().to_string(),
            params: plan.canonical(),
            source: source.name(),
        })
    }

    fn resolve_gemm(&self, key: &GemmTuneKey) -> (GemmBlockPlan, TuneSource) {
        if let Some(PinnedParam::GemmBlocking(plan)) = self.pins.get(key.kernel()) {
            (*plan, TuneSource::Pinned)
        } else if let Some(plan) = self
            .rows
            .get(&(key.kernel().to_string(), key.shape_class().to_string()))
            .and_then(|row| GemmBlockPlan::parse(&row.params))
        {
            (plan, TuneSource::Tuned)
        } else {
            (GemmBlockPlan::COLD_START, TuneSource::ColdStart)
        }
    }

    /// Validate a measured GEMM selection without mutating this tuner.
    ///
    /// The returned value can be persisted first and committed only after the
    /// external cache write succeeds, so a ledger failure cannot leave a
    /// process-local row that was never durably recorded.
    ///
    /// # Errors
    /// [`TuneError`] when evidence is invalid, candidate labels are not
    /// canonical plans, the selected plan is not the deterministic evidence
    /// argmin, the generated row exceeds persistence bounds or cannot be
    /// reparsed canonically, or the refresh counter is exhausted.
    pub fn prepare_gemm_row(
        &self,
        key: &GemmTuneKey,
        plan: GemmBlockPlan,
        evidence: TuneEvidence,
    ) -> Result<PreparedGemmRow, TuneError> {
        evidence.validate()?;
        validate_gemm_selection(plan, &evidence)?;
        let map_key = (key.kernel().to_string(), key.shape_class().to_string());
        let expected_refresh = self.rows.get(&map_key).map(|row| row.refresh);
        let refresh = expected_refresh.map_or(Ok(1), |refresh| {
            refresh.checked_add(1).ok_or_else(|| TuneError {
                detail: format!(
                    "refresh counter exhausted for kernel {:?} shape {:?}",
                    key.kernel(),
                    key.shape_class()
                ),
            })
        })?;
        let row = TuneRow {
            kernel: key.kernel().to_string(),
            shape_class: key.shape_class().to_string(),
            machine: self.fingerprint,
            params: plan.canonical(),
            evidence,
            refresh,
        };
        row.validated_generated_json()?;
        Ok(PreparedGemmRow {
            key: key.clone(),
            row,
            expected_refresh,
        })
    }

    /// Parse and validate a cached GEMM row against the exact requested key
    /// without mutating this tuner.
    ///
    /// This binds the embedded kernel, shape class, machine fingerprint,
    /// canonical plan, ranked evidence, and evidence argmin before a session
    /// compares the ledger's separate `params` field and commits the row.
    ///
    /// # Errors
    /// [`TuneError`] on any mismatch or when the requested key already has a
    /// local row (cache adoption is a seeding operation, not a rollback path).
    pub fn prepare_adopt_gemm_row_json(
        &self,
        key: &GemmTuneKey,
        line: &str,
    ) -> Result<PreparedGemmRow, TuneError> {
        let row = parse_row(line).ok_or_else(|| TuneError {
            detail: "external GEMM tune row is not a canonical validated row line".to_string(),
        })?;
        if row.machine != self.fingerprint {
            return Err(TuneError {
                detail: format!(
                    "external GEMM tune row is stale: machine {:016x} != this machine {:016x}",
                    row.machine, self.fingerprint
                ),
            });
        }
        if row.kernel != key.kernel() || row.shape_class != key.shape_class() {
            return Err(TuneError {
                detail: format!(
                    "external GEMM tune row key mismatch: embedded {:?} × {:?}, requested {:?} × {:?}",
                    row.kernel,
                    row.shape_class,
                    key.kernel(),
                    key.shape_class()
                ),
            });
        }
        validate_scoped_gemm_row(&row)?;
        let map_key = (key.kernel().to_string(), key.shape_class().to_string());
        if self.rows.contains_key(&map_key) {
            return Err(TuneError {
                detail: "external GEMM tune-row adoption cannot replace a local row".to_string(),
            });
        }
        Ok(PreparedGemmRow {
            key: key.clone(),
            row,
            expected_refresh: None,
        })
    }

    /// Install a previously validated GEMM row if tuner state still matches
    /// the state observed during preparation.
    ///
    /// # Errors
    /// [`TuneError`] when the prepared row belongs to another tuner, is no
    /// longer canonically bounded, or the row changed between preparation and
    /// commit.
    pub fn commit_gemm_row(&mut self, prepared: PreparedGemmRow) -> Result<(), TuneError> {
        if prepared.row.machine != self.fingerprint
            || prepared.row.kernel != prepared.key.kernel()
            || prepared.row.shape_class != prepared.key.shape_class()
        {
            return Err(TuneError {
                detail: "prepared GEMM row does not belong to this tuner/key".to_string(),
            });
        }
        prepared.row.validated_generated_json()?;
        validate_scoped_gemm_row(&prepared.row)?;
        let map_key = (
            prepared.key.kernel().to_string(),
            prepared.key.shape_class().to_string(),
        );
        let actual_refresh = self.rows.get(&map_key).map(|row| row.refresh);
        if actual_refresh != prepared.expected_refresh {
            return Err(TuneError {
                detail: "prepared GEMM row is stale because tuner state changed before commit"
                    .to_string(),
            });
        }
        self.rows.insert(map_key, prepared.row);
        Ok(())
    }

    /// Record a measured GEMM sweep row: the winning plan plus its RANKED
    /// wall-time candidate evidence. Same key replaces the row and
    /// increments `refresh` (idempotence witness).
    ///
    /// # Errors
    /// [`TuneError`] when evidence makes no ranked-candidate claim, contains
    /// non-canonical plan labels, or selects anything other than its argmin.
    pub fn record_gemm_row(
        &mut self,
        key: &GemmTuneKey,
        plan: GemmBlockPlan,
        evidence: TuneEvidence,
    ) -> Result<(), TuneError> {
        let prepared = self.prepare_gemm_row(key, plan, evidence)?;
        self.commit_gemm_row(prepared)
    }

    /// Canonical JSON for a validated GEMM row under this exact key.
    #[must_use]
    pub fn gemm_row_json(&self, key: &GemmTuneKey) -> Option<String> {
        self.rows
            .get(&(key.kernel().to_string(), key.shape_class().to_string()))
            .map(TuneRow::to_json)
    }

    /// The canonical JSON line for one general row (what an external cache
    /// stores and [`Tuner::adopt_row_json`] reads). GEMM consumers use
    /// [`Tuner::gemm_row_json`] and expected-key adoption instead.
    #[must_use]
    pub fn row_json(&self, kernel: &str, shape_class: &str) -> Option<String> {
        self.rows
            .get(&(kernel.to_string(), shape_class.to_string()))
            .map(TuneRow::to_json)
    }

    /// Adopt one canonical JSON tune-row line from an external cache (the
    /// ledger `tune` table). Fail-closed filters: the line must parse as
    /// a canonical row, its evidence must re-validate, and its machine
    /// fingerprint must equal THIS tuner's — a stale (other-machine) or
    /// tampered row is refused, never silently adopted. The stored
    /// `refresh` counter is preserved.
    ///
    /// # Errors
    /// [`TuneError`] naming the refusal.
    pub fn adopt_row_json(&mut self, line: &str) -> Result<(), TuneError> {
        let row = parse_row(line).ok_or_else(|| TuneError {
            detail: "external tune row is not a canonical row line".to_string(),
        })?;
        if row.kernel.starts_with(GEMM_KERNEL_PREFIX) {
            return Err(TuneError {
                detail: "GEMM rows require prepare_adopt_gemm_row_json with an expected scoped key"
                    .to_string(),
            });
        }
        if row.machine != self.fingerprint {
            return Err(TuneError {
                detail: format!(
                    "external tune row is stale: machine {:016x} != this machine {:016x}",
                    row.machine, self.fingerprint
                ),
            });
        }
        self.rows
            .insert((row.kernel.clone(), row.shape_class.clone()), row);
        Ok(())
    }

    /// Validate and install one GEMM cache row against an exact requested key.
    /// Sessions that also persist a separate params column should instead call
    /// [`Tuner::prepare_adopt_gemm_row_json`], compare that column with
    /// [`PreparedGemmRow::params_json`], then call
    /// [`Tuner::commit_gemm_row`].
    ///
    /// # Errors
    /// [`TuneError`] naming any semantic or identity mismatch.
    pub fn adopt_gemm_row_json(&mut self, key: &GemmTuneKey, line: &str) -> Result<(), TuneError> {
        let prepared = self.prepare_adopt_gemm_row_json(key, line)?;
        self.commit_gemm_row(prepared)
    }

    /// Validate and pin canonical params from a recorded decision.
    ///
    /// The accepted replay forms are exactly `edge={4,8,16}` for tile
    /// kernels, `schedule={bandwidth-rich,bandwidth-starved}` for the schedule
    /// key, and a canonical bounded GEMM plan for a canonical scoped GEMM key.
    /// Invalid spellings fail closed instead of resolving to a default that
    /// would be falsely recorded as pinned.
    ///
    /// # Errors
    /// Returns [`TuneError`] when `params` is not canonical for `kernel`.
    pub fn pin(
        &mut self,
        kernel: impl Into<String>,
        params: impl Into<String>,
    ) -> Result<(), TuneError> {
        let kernel = kernel.into();
        let params = params.into();
        require_decision_kernel(&kernel)?;
        let parsed = PinnedParam::parse(&kernel, &params).ok_or_else(|| TuneError {
            detail: format!("invalid pinned params {params:?} for kernel {kernel:?}"),
        })?;
        self.pins.insert(kernel, parsed);
        Ok(())
    }

    /// The currently retained decision window in dispatch order.
    ///
    /// This compatibility view is not, by itself, proof that the prefix is
    /// complete. Long-running callers must check [`Tuner::decision_history`]
    /// before using it as replay input and persist production dispatch receipts
    /// to the Design Ledger.
    #[must_use]
    pub fn decisions(&self) -> &[TuningDecision] {
        &self.decisions
    }

    /// Metadata for the bounded decision window.
    #[must_use]
    pub fn decision_history(&self) -> TuningDecisionHistory<'_> {
        TuningDecisionHistory {
            decisions: &self.decisions,
            evicted: self.evicted_decisions,
            retained_bytes: self.decision_bytes,
        }
    }

    /// Resolve a tile edge for `kernel`: pins beat tuned rows beat the
    /// cold-start default (8³ — plan §5.3). The decision is recorded.
    ///
    /// # Errors
    /// [`TuneError`] when the kernel identity is blank or exceeds the bounded
    /// tune-string domain.
    pub fn tile_edge_for(&mut self, kernel: &str) -> Result<(TileEdge, TuneSource), TuneError> {
        require_decision_kernel(kernel)?;
        let (params, source) = self.resolve(kernel)?;
        let edge = match params.as_str() {
            "edge=4" => TileEdge::E4,
            "edge=16" => TileEdge::E16,
            _ => TileEdge::E8,
        };
        Ok((edge, source))
    }

    /// Resolve the schedule kind: pins beat tuned beat the cold-start
    /// default (bandwidth-rich — harmless on starved machines, just less
    /// blocked). The decision is recorded.
    pub fn schedule(&mut self) -> (ScheduleKind, TuneSource) {
        let (params, source) = self
            .resolve(SCHEDULE_KERNEL)
            .expect("the reserved schedule decision is statically bounded");
        let kind = if params == "schedule=bandwidth-starved" {
            ScheduleKind::BandwidthStarved
        } else {
            ScheduleKind::BandwidthRich
        };
        (kind, source)
    }

    fn resolve(&mut self, kernel: &str) -> Result<(String, TuneSource), TuneError> {
        let (params, source) = if let Some(p) = self.pins.get(kernel) {
            (p.canonical(), TuneSource::Pinned)
        } else if let Some(row) = self
            .rows
            .get(&(kernel.to_string(), SHAPE_DEFAULT.to_string()))
        {
            (row.params.clone(), TuneSource::Tuned)
        } else {
            let default = if kernel == SCHEDULE_KERNEL {
                "schedule=bandwidth-rich"
            } else {
                "edge=8"
            };
            (default.to_string(), TuneSource::ColdStart)
        };
        self.record_decision(TuningDecision {
            kernel: kernel.to_string(),
            params: params.clone(),
            source: source.name(),
        })?;
        Ok((params, source))
    }

    fn record_decision(&mut self, decision: TuningDecision) -> Result<(), TuneError> {
        let bytes = decision.retained_bytes().ok_or_else(|| TuneError {
            detail: "tuning decision payload length overflowed usize".to_string(),
        })?;
        if bytes > MAX_RETAINED_TUNING_DECISION_BYTES {
            return Err(TuneError {
                detail: format!(
                    "tuning decision exceeds the {MAX_RETAINED_TUNING_DECISION_BYTES}-byte retained-payload limit"
                ),
            });
        }

        let exceeds_count = self.decisions.len() >= MAX_RETAINED_TUNING_DECISIONS;
        let exceeds_bytes = self
            .decision_bytes
            .checked_add(bytes)
            .is_none_or(|total| total > MAX_RETAINED_TUNING_DECISION_BYTES);
        if exceeds_count || exceeds_bytes {
            // Evict to half capacity in one deterministic oldest-prefix batch.
            // This keeps append cost amortized rather than shifting a full Vec
            // on every dispatch once the window reaches capacity.
            let target_count = MAX_RETAINED_TUNING_DECISIONS / 2;
            let target_bytes = MAX_RETAINED_TUNING_DECISION_BYTES / 2;
            let mut remove = 0_usize;
            let mut remaining_bytes = self.decision_bytes;
            while remove < self.decisions.len()
                && (self.decisions.len() - remove > target_count
                    || remaining_bytes
                        .checked_add(bytes)
                        .is_none_or(|total| total > target_bytes))
            {
                let old_bytes = self.decisions[remove]
                    .retained_bytes()
                    .expect("retained decision lengths were checked on insertion");
                remaining_bytes = remaining_bytes
                    .checked_sub(old_bytes)
                    .expect("retained decision byte accounting is exact");
                remove += 1;
            }
            self.decisions.drain(..remove);
            self.decision_bytes = remaining_bytes;
            self.evicted_decisions = self.evicted_decisions.saturating_add(remove);
        }

        self.decision_bytes = self
            .decision_bytes
            .checked_add(bytes)
            .filter(|&total| total <= MAX_RETAINED_TUNING_DECISION_BYTES)
            .ok_or_else(|| TuneError {
                detail: "tuning decision window failed its byte bound after eviction".to_string(),
            })?;
        if self.decisions.len() >= MAX_RETAINED_TUNING_DECISIONS {
            return Err(TuneError {
                detail: "tuning decision window failed its count bound after eviction".to_string(),
            });
        }
        self.decisions.push(decision);
        Ok(())
    }

    /// Run the calibration pass: stencil-edge sweep through the REAL tile
    /// pool, reduction-cost measurement, steal-cost measurement, and
    /// schedule selection from the probe's measured bandwidth. Idempotent:
    /// same fingerprint keys are replaced, `refresh` increments. Returns
    /// the machine calibration report (canonical JSON, ledger-bound).
    ///
    /// # Errors
    /// Returns [`TuneError`] without measuring or mutating rows when the
    /// probe's stable fingerprint differs from this tuner's machine key.
    pub fn calibrate(
        &mut self,
        probe: &fs_substrate::CapabilityProbe,
    ) -> Result<String, TuneError> {
        self.require_matching_probe(probe)?;
        self.require_refresh_capacity()?;
        let workers = (probe.logical_cpus as usize).clamp(1, 8);
        let topo = CcdTopology::from_probe(probe);
        let (schedule, throughput_milli) = schedule_measurement(probe)?;
        self.calibrate_stencil(workers, topo)?;
        self.calibrate_reduction()?;
        self.calibrate_steals(topo)?;
        self.class_work_count_row(probe, topo)?;
        self.insert_schedule_row(schedule, throughput_milli)?;
        Ok(self.calibration_report())
    }

    fn calibrate_stencil(&mut self, workers: usize, topo: CcdTopology) -> Result<(), TuneError> {
        // Tile-edge sweep: minimum-of-2 per edge, argmin wins; ties break to
        // the LOWER candidate index (deterministic tie law). Both exact
        // samples survive in the evidence row; the summary is re-derived on
        // load rather than trusted as an independent claim.
        let mut measured = Vec::new();
        for edge in [TileEdge::E4, TileEdge::E8, TileEdge::E16] {
            let grid = fs_substrate::tile::TileGrid::new([48, 48, 48], edge)
                .expect("48-cube fits every edge");
            let mut field = fs_substrate::field::TiledField::new(grid, 0.0f32);
            let dims = field.grid().cell_dims();
            for z in 0..dims[2] {
                for y in 0..dims[1] {
                    for x in 0..dims[0] {
                        field.set([x, y, z], ((x * 31 + y * 7 + z) % 97) as f32 / 97.0);
                    }
                }
            }
            let kernel = CalStencil { field };
            let pool = TilePool::new(PoolConfig::new(workers, topo, 0x7C4E));
            let mut samples_ns = Vec::with_capacity(2);
            for _ in 0..2 {
                let gate = CancelGate::new();
                let t0 = std::time::Instant::now();
                let (r, _) = pool.run_with_gate(&kernel, &gate);
                let ns = duration_ns(t0.elapsed(), "stencil calibration")?;
                r.expect("calibration stencil runs");
                samples_ns.push(ns);
            }
            measured.push(TuneObservation::wall_time(
                format!("edge={}", edge.cells()),
                samples_ns,
            )?);
        }
        let winner = measured
            .iter()
            .enumerate()
            .min_by_key(|(index, observation)| {
                (
                    observation
                        .wall_minimum()
                        .expect("stencil candidates are wall-time observations"),
                    *index,
                )
            })
            .map(|(_, observation)| observation.label().to_string())
            .expect("three candidates");
        self.insert_row(
            STENCIL_KERNEL,
            SHAPE_DEFAULT,
            winner,
            TuneEvidence::ranked_wall_times(measured)?,
        )
    }

    fn calibrate_reduction(&mut self) -> Result<(), TuneError> {
        // Reduction cost: deterministic compensated sum over 100k terms.
        let xs: Vec<f64> = (0..100_000).map(|i| 1.0 / f64::from(i + 1)).collect();
        let t0 = std::time::Instant::now();
        let s = crate::reduce::det_sum(&xs);
        let red_ns = duration_ns(t0.elapsed(), "deterministic reduction calibration")?;
        assert!(s.is_finite());
        self.insert_row(
            "det-sum-f64",
            "100k",
            "block=256".to_string(),
            TuneEvidence::new(vec![TuneObservation::wall_time("block=256", vec![red_ns])?])?,
        )
    }

    fn calibrate_steals(&mut self, topo: CcdTopology) -> Result<(), TuneError> {
        // Steal cost: an imbalanced 2-worker run forces steals; record the
        // exact run duration and steal counter as different units.
        let pool = TilePool::new(PoolConfig {
            quantum_weights: vec![15, 1],
            ..PoolConfig::new(2, topo, 0x7C4E)
        });
        let gate = CancelGate::new();
        let t0 = std::time::Instant::now();
        let (r, report) = pool.run_with_gate(&StealProbe, &gate);
        let steal_ns = duration_ns(t0.elapsed(), "steal calibration")?;
        r.expect("steal probe runs");
        self.insert_row(
            "steal-probe",
            "2w-imbalanced",
            format!("steals={}", report.steals),
            TuneEvidence::new(vec![
                TuneObservation::wall_time("run", vec![steal_ns])?,
                TuneObservation::work_count("steals", WorkUnit::Steals, report.steals)?,
            ])?,
        )
    }

    fn insert_schedule_row(
        &mut self,
        kind: ScheduleKind,
        throughput_milli: u64,
    ) -> Result<(), TuneError> {
        self.insert_row(
            SCHEDULE_KERNEL,
            SHAPE_DEFAULT,
            format!("schedule={}", kind.name()),
            TuneEvidence::new(vec![TuneObservation::throughput(
                "per-core",
                ThroughputUnit::GigabytesPerSecond,
                throughput_milli,
            )?])?,
        )
    }

    fn calibration_report(&self) -> String {
        let mut s = String::with_capacity(512);
        let _ = write!(
            s,
            "{{\"machine\":\"{:016x}\",\"evidence_version\":{},\"rows\":[",
            self.fingerprint,
            TuneEvidence::VERSION
        );
        for (i, row) in self.rows.values().enumerate() {
            if i > 0 {
                s.push(',');
            }
            s.push_str(&row.to_json());
        }
        s.push_str("]}");
        s
    }

    fn require_matching_probe(
        &self,
        probe: &fs_substrate::CapabilityProbe,
    ) -> Result<(), TuneError> {
        let probe_fingerprint = probe.fingerprint();
        if probe_fingerprint == self.fingerprint {
            return Ok(());
        }
        Err(TuneError {
            detail: format!(
                "calibration probe fingerprint {probe_fingerprint:016x} does not match tuner \
                 fingerprint {:016x}",
                self.fingerprint
            ),
        })
    }

    fn require_refresh_capacity(&self) -> Result<(), TuneError> {
        for (kernel, shape) in [
            (STENCIL_KERNEL, SHAPE_DEFAULT),
            ("det-sum-f64", "100k"),
            ("steal-probe", "2w-imbalanced"),
            (CLASS_WORK_KERNEL, SHAPE_DEFAULT),
            (SCHEDULE_KERNEL, SHAPE_DEFAULT),
        ] {
            if self
                .rows
                .get(&(kernel.to_string(), shape.to_string()))
                .is_some_and(|row| row.refresh == u32::MAX)
            {
                return Err(TuneError {
                    detail: format!(
                        "refresh counter exhausted for kernel {kernel:?} shape {shape:?}"
                    ),
                });
            }
        }
        Ok(())
    }

    /// Per-worker work-count row (fz2.2): FULL parallelism, uniform tiles.
    /// Counts are an observed allocation signal, not a throughput rate because
    /// the row does not carry a per-worker elapsed-time denominator. Records
    /// the sorted distribution and fast:slow count ratio; the quantum_weights
    /// verification lives in the machine A/B lane.
    fn class_work_count_row(
        &mut self,
        probe: &fs_substrate::CapabilityProbe,
        topo: CcdTopology,
    ) -> Result<(), TuneError> {
        let full = (probe.logical_cpus as usize).max(1);
        let pool = TilePool::new(PoolConfig::new(full, topo, 0x7C4E));
        let gate = CancelGate::new();
        let (r, report) = pool.run_with_gate(&ClassProbe, &gate);
        r.expect("class probe runs");
        let mut counts = report.tiles_by_worker.clone();
        counts.sort_unstable_by(|a, b| b.cmp(a));
        let fast = counts.first().copied().unwrap_or(0).max(1);
        let slow = counts.last().copied().unwrap_or(0).max(1);
        let ratio_x1000 = u64::try_from(u128::from(slow) * 1000 / u128::from(fast))
            .expect("a completed-tile ratio is representable as u64");
        self.insert_row(
            CLASS_WORK_KERNEL,
            SHAPE_DEFAULT,
            format!("workers={full};slow_over_fast_x1000={ratio_x1000}"),
            TuneEvidence::new(
                counts
                    .iter()
                    .enumerate()
                    .map(|(i, &count)| {
                        TuneObservation::work_count(
                            format!("worker-{i}"),
                            WorkUnit::CompletedTiles,
                            count,
                        )
                    })
                    .collect::<Result<Vec<_>, _>>()?,
            )?,
        )
    }

    fn insert_row(
        &mut self,
        kernel: &str,
        shape_class: &str,
        params: String,
        evidence: TuneEvidence,
    ) -> Result<(), TuneError> {
        evidence.validate()?;
        let key = (kernel.to_string(), shape_class.to_string());
        let refresh = self.rows.get(&key).map_or(Ok(1), |row| {
            row.refresh.checked_add(1).ok_or_else(|| TuneError {
                detail: format!(
                    "refresh counter exhausted for kernel {kernel:?} shape {shape_class:?}"
                ),
            })
        })?;
        let row = TuneRow {
            kernel: kernel.to_string(),
            shape_class: shape_class.to_string(),
            machine: self.fingerprint,
            params,
            evidence,
            refresh,
        };
        row.validated_generated_json()?;
        self.rows.insert(key, row);
        Ok(())
    }

    /// Persist the table (JSON-lines, one row per line — the ledger `tune`
    /// schema's file edition).
    ///
    /// # Errors
    /// [`TuneError`] when any row is no longer canonically bounded, the total
    /// store exceeds its persistence budget, or the final write fails.
    pub fn save(&self, path: &std::path::Path) -> Result<(), TuneError> {
        let mut out = String::new();
        for row in self.rows.values() {
            out.push_str(&row.validated_generated_json()?);
            out.push('\n');
            if out.len() > MAX_TUNE_STORE_BYTES {
                return Err(TuneError {
                    detail: format!(
                        "tune store exceeds the {MAX_TUNE_STORE_BYTES}-byte persistence limit"
                    ),
                });
            }
        }
        std::fs::write(path, out).map_err(|e| TuneError {
            detail: format!("writing {}: {e}", path.display()),
        })
    }

    /// Load a table, KEEPING only rows for `fingerprint` (rows from other
    /// machines are stale by definition — fingerprint drift detection). A
    /// missing file is a cold start, not an error.
    ///
    /// # Errors
    /// [`TuneError`] on unreadable/corrupt stores.
    pub fn load(path: &std::path::Path, fingerprint: u64) -> Result<Self, TuneError> {
        let mut tuner = Tuner::cold(fingerprint);
        let file = match std::fs::File::open(path) {
            Ok(file) => file,
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(tuner),
            Err(e) => {
                return Err(TuneError {
                    detail: format!("reading {}: {e}", path.display()),
                });
            }
        };
        let mut bytes = Vec::new();
        let read_limit = u64::try_from(MAX_TUNE_STORE_BYTES + 1)
            .expect("the fixed tune-store limit is representable as u64");
        file.take(read_limit)
            .read_to_end(&mut bytes)
            .map_err(|error| TuneError {
                detail: format!("reading {}: {error}", path.display()),
            })?;
        if bytes.len() > MAX_TUNE_STORE_BYTES {
            return Err(TuneError {
                detail: format!(
                    "tune store {} exceeds the {MAX_TUNE_STORE_BYTES}-byte limit; remove it and recalibrate",
                    path.display()
                ),
            });
        }
        let text = String::from_utf8(bytes).map_err(|_| TuneError {
            detail: format!(
                "tune store {} is not UTF-8; remove it and recalibrate",
                path.display()
            ),
        })?;
        for (lineno, line) in text.lines().enumerate() {
            if line.len() > MAX_TUNE_ROW_BYTES {
                return Err(TuneError {
                    detail: format!(
                        "tune row at {}:{} exceeds the {MAX_TUNE_ROW_BYTES}-byte limit; remove the store and recalibrate",
                        path.display(),
                        lineno + 1
                    ),
                });
            }
            let row = parse_row(line).ok_or_else(|| TuneError {
                detail: format!(
                    "corrupt row at {}:{}; remove the store and recalibrate",
                    path.display(),
                    lineno + 1
                ),
            })?;
            if row.machine == fingerprint {
                let key = (row.kernel.clone(), row.shape_class.clone());
                if tuner.rows.insert(key, row).is_some() {
                    return Err(TuneError {
                        detail: format!(
                            "duplicate tune row for the selected fingerprint at {}:{}; remove the store and recalibrate",
                            path.display(),
                            lineno + 1
                        ),
                    });
                }
            }
        }
        Ok(tuner)
    }
}

/// Uniform CPU-bound tiles for the per-worker work-distribution probe (fz2.2).
/// Every tile has the same loop budget, so completed-tile counts are comparable
/// work counters. Stealing and per-worker active time prevent interpreting
/// those counters as standalone throughput rates.
struct ClassProbe;

impl TileKernel for ClassProbe {
    type Out = u64;

    fn tiles(&self) -> TilePlan {
        TilePlan::new("tune/class-probe", 512)
    }

    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<crate::Cancelled, u64> {
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(crate::Cancelled);
        }
        let mut acc = (tile as f64).mul_add(1e-9, 1.0);
        for _ in 0..1_000_000 {
            acc = acc.mul_add(0.999_999_9, 1.0e-9);
        }
        ControlFlow::Continue(acc.to_bits() & 1)
    }
}

/// Tiny imbalance workload for the steal-cost probe.
struct StealProbe;

impl TileKernel for StealProbe {
    type Out = u64;

    fn tiles(&self) -> TilePlan {
        TilePlan::new("tune/steal-probe", 4096)
    }

    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<crate::Cancelled, u64> {
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(crate::Cancelled);
        }
        let mut acc = tile;
        for i in 0..100 {
            acc = acc.wrapping_mul(6364136223846793005).wrapping_add(i);
        }
        ControlFlow::Continue(acc & 1)
    }
}

fn duration_ns(duration: std::time::Duration, context: &str) -> Result<u64, TuneError> {
    let nanoseconds = u64::try_from(duration.as_nanos()).map_err(|_| TuneError {
        detail: format!("{context} duration exceeds the tune evidence u64 nanosecond range"),
    })?;
    Ok(nanoseconds.max(1))
}

fn schedule_measurement(
    probe: &fs_substrate::CapabilityProbe,
) -> Result<(ScheduleKind, u64), TuneError> {
    if !probe.measured.all_core_gbs.is_finite() || probe.measured.all_core_gbs < 0.0 {
        return Err(TuneError {
            detail: format!(
                "aggregate throughput {:?} GB/s is not a finite non-negative measurement",
                probe.measured.all_core_gbs
            ),
        });
    }
    let logical_cpus = probe.logical_cpus as usize;
    if logical_cpus > MAX_TUNE_OBSERVATIONS {
        return Err(TuneError {
            detail: format!(
                "probe reports {logical_cpus} logical CPUs, exceeding the {MAX_TUNE_OBSERVATIONS}-worker evidence limit"
            ),
        });
    }
    if logical_cpus == 0 && probe.measured.all_core_gbs > 0.0 {
        return Err(TuneError {
            detail: "positive aggregate throughput cannot be attributed to zero logical CPUs"
                .to_string(),
        });
    }
    let per_core = if probe.logical_cpus > 0 {
        probe.measured.all_core_gbs / f64::from(probe.logical_cpus)
    } else {
        0.0
    };
    let throughput_milli = throughput_milli_units(per_core)?;
    let kind = if probe.measured.all_core_gbs > 0.0 && per_core < 8.0 {
        ScheduleKind::BandwidthStarved
    } else {
        ScheduleKind::BandwidthRich
    };
    Ok((kind, throughput_milli))
}

fn throughput_milli_units(value: f64) -> Result<u64, TuneError> {
    let scaled = value * 1_000.0;
    if !value.is_finite() || value < 0.0 || !scaled.is_finite() || scaled >= u64::MAX as f64 {
        return Err(TuneError {
            detail: format!(
                "per-core throughput {value:?} GB/s cannot be represented as exact milli-units"
            ),
        });
    }
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    Ok(scaled.round() as u64)
}

/// Exact parser for the canonical row writer above. Field reordering,
/// duplicate fields, non-canonical numbers, and trailing content are
/// corruption rather than dialects of the tune-store schema.
fn parse_row(line: &str) -> Option<TuneRow> {
    if line.len() > MAX_TUNE_ROW_BYTES {
        return None;
    }
    let mut parser = RowParser { rest: line };
    parser.take("{\"kernel\":")?;
    let kernel = parser.string()?;
    parser.take(",\"shape_class\":")?;
    let shape_class = parser.string()?;
    parser.take(",\"machine\":")?;
    let machine_hex = parser.string()?;
    if machine_hex.len() != 16
        || !machine_hex
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    {
        return None;
    }
    let machine = u64::from_str_radix(&machine_hex, 16).ok()?;
    parser.take(",\"params\":")?;
    let params = parser.string()?;
    parser.take(",\"evidence\":")?;
    let evidence = parse_evidence(&mut parser)?;
    parser.take(",\"refresh\":")?;
    let refresh = u32::try_from(parser.canonical_u64()?).ok()?;
    parser.take("}")?;

    if !parser.rest.is_empty()
        || kernel.trim().is_empty()
        || shape_class.trim().is_empty()
        || params.trim().is_empty()
        || refresh == 0
        || (kernel == STENCIL_KERNEL && PinnedParam::parse(&kernel, &params).is_none())
        || (kernel == SCHEDULE_KERNEL && PinnedParam::parse(&kernel, &params).is_none())
    {
        return None;
    }

    let row = TuneRow {
        kernel,
        shape_class,
        machine,
        params,
        evidence,
        refresh,
    };
    if row.kernel.starts_with(GEMM_KERNEL_PREFIX) && validate_scoped_gemm_row(&row).is_err() {
        return None;
    }
    (row.to_json() == line).then_some(row)
}

fn parse_evidence(parser: &mut RowParser<'_>) -> Option<TuneEvidence> {
    parser.take("{\"version\":")?;
    let version = u32::try_from(parser.canonical_u64()?).ok()?;
    if version != TuneEvidence::VERSION {
        return None;
    }
    parser.take(",\"observations\":[")?;
    let mut observations = Vec::new();
    loop {
        if observations.len() == MAX_TUNE_OBSERVATIONS {
            return None;
        }
        observations.push(parse_observation(parser)?);
        if parser.rest.starts_with(',') {
            parser.take(",")?;
        } else {
            break;
        }
    }
    parser.take("],\"candidate_separation_ppm\":")?;
    let stored_separation = if parser.rest.starts_with("null") {
        parser.take("null")?;
        None
    } else {
        Some(u32::try_from(parser.canonical_u64()?).ok()?)
    };
    parser.take("}")?;

    let evidence = if stored_separation.is_some() {
        TuneEvidence::ranked_wall_times(observations).ok()?
    } else {
        TuneEvidence::new(observations).ok()?
    };
    (evidence.candidate_separation_ppm == stored_separation).then_some(evidence)
}

fn parse_observation(parser: &mut RowParser<'_>) -> Option<TuneObservation> {
    parser.take("{\"kind\":")?;
    let kind = parser.string()?;
    parser.take(",\"label\":")?;
    let label = parser.string()?;
    match kind.as_str() {
        "wall-time" => {
            parser.take(",\"samples_ns\":[")?;
            let mut samples_ns = Vec::new();
            loop {
                if samples_ns.len() == MAX_WALL_TIME_SAMPLES {
                    return None;
                }
                samples_ns.push(parser.canonical_u64()?);
                if parser.rest.starts_with(',') {
                    parser.take(",")?;
                } else {
                    break;
                }
            }
            parser.take("],\"summary\":{\"minimum_ns\":")?;
            let minimum_ns = parser.canonical_u64()?;
            parser.take(",\"maximum_ns\":")?;
            let maximum_ns = parser.canonical_u64()?;
            parser.take("}}")?;
            let observation = TuneObservation::wall_time(label, samples_ns).ok()?;
            match &observation {
                TuneObservation::WallTime { summary, .. }
                    if summary.minimum_ns == minimum_ns && summary.maximum_ns == maximum_ns =>
                {
                    Some(observation)
                }
                _ => None,
            }
        }
        "work-count" => {
            parser.take(",\"unit\":")?;
            let unit = WorkUnit::parse(&parser.string()?)?;
            parser.take(",\"count\":")?;
            let count = parser.canonical_u64()?;
            parser.take("}")?;
            TuneObservation::work_count(label, unit, count).ok()
        }
        "throughput" => {
            parser.take(",\"unit\":")?;
            let unit = ThroughputUnit::parse(&parser.string()?)?;
            parser.take(",\"milli_units\":")?;
            let milli_units = parser.canonical_u64()?;
            parser.take("}")?;
            TuneObservation::throughput(label, unit, milli_units).ok()
        }
        _ => None,
    }
}

struct RowParser<'a> {
    rest: &'a str,
}

impl RowParser<'_> {
    fn take(&mut self, expected: &str) -> Option<()> {
        self.rest = self.rest.strip_prefix(expected)?;
        Some(())
    }

    fn string(&mut self) -> Option<String> {
        self.take("\"")?;
        let mut out = String::new();
        loop {
            let ch = self.rest.chars().next()?;
            self.rest = &self.rest[ch.len_utf8()..];
            match ch {
                '"' => return Some(out),
                '\\' => {
                    let escaped = self.rest.chars().next()?;
                    self.rest = &self.rest[escaped.len_utf8()..];
                    let decoded = match escaped {
                        '"' => '"',
                        '\\' => '\\',
                        'n' => '\n',
                        'r' => '\r',
                        't' => '\t',
                        'u' => {
                            let hex = self.rest.get(..4)?;
                            if !hex.bytes().all(|byte| byte.is_ascii_hexdigit()) {
                                return None;
                            }
                            let code = u32::from_str_radix(hex, 16).ok()?;
                            self.rest = &self.rest[4..];
                            char::from_u32(code)?
                        }
                        _ => return None,
                    };
                    push_bounded_char(&mut out, decoded)?;
                }
                c if c.is_control() => return None,
                c => push_bounded_char(&mut out, c)?,
            }
        }
    }

    fn canonical_u64(&mut self) -> Option<u64> {
        let len = self.rest.bytes().take_while(u8::is_ascii_digit).count();
        if len == 0 {
            return None;
        }
        let token = &self.rest[..len];
        if token.len() > 1 && token.starts_with('0') {
            return None;
        }
        self.rest = &self.rest[len..];
        token.parse().ok()
    }
}

fn push_bounded_char(out: &mut String, ch: char) -> Option<()> {
    if out.len().checked_add(ch.len_utf8())? > MAX_TUNE_STRING_BYTES {
        return None;
    }
    out.push(ch);
    Some(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn gemm_identity(
        requested_threads: usize,
        thread_budget: usize,
        probe_dims: [usize; 3],
        tier: &str,
        placement: &str,
        implementation: &str,
        build: &str,
    ) -> GemmExecutionIdentity {
        GemmExecutionIdentity::new(
            requested_threads,
            thread_budget,
            u64::MAX,
            probe_dims,
            tier,
            placement,
            implementation,
            build,
        )
        .expect("test execution identity")
    }

    fn gemm_key_with(identity: GemmExecutionIdentity) -> GemmTuneKey {
        GemmTuneKey::new(format!("{GEMM_KERNEL_PREFIX}1"), "m64-n128-k32", identity)
            .expect("test GEMM key")
    }

    fn gemm_key() -> GemmTuneKey {
        gemm_key_with(gemm_identity(
            4,
            4,
            [64, 128, 32],
            "avx2-fma",
            "unpinned",
            "fs-la-parallel-v1",
            "release-opt3-cgu1-build-a",
        ))
    }

    fn ranked_gemm_evidence(
        first: (GemmBlockPlan, u64),
        second: (GemmBlockPlan, u64),
    ) -> TuneEvidence {
        TuneEvidence::ranked_wall_times(vec![
            TuneObservation::wall_time(first.0.canonical(), vec![first.1]).expect("first timing"),
            TuneObservation::wall_time(second.0.canonical(), vec![second.1])
                .expect("second timing"),
        ])
        .expect("ranked GEMM evidence")
    }

    fn work_count_evidence_with_label_lengths(label_lengths: &[usize]) -> TuneEvidence {
        TuneEvidence::new(
            label_lengths
                .iter()
                .enumerate()
                .map(|(index, &length)| {
                    let prefix = format!("work-{index}:");
                    assert!(length >= prefix.len(), "test label length is too short");
                    TuneObservation::work_count(
                        format!("{prefix}{}", "x".repeat(length - prefix.len())),
                        WorkUnit::Steals,
                        index as u64,
                    )
                })
                .collect::<Result<Vec<_>, _>>()
                .expect("bounded work-count observations"),
        )
        .expect("bounded work-count evidence")
    }

    fn work_count_row(evidence: TuneEvidence) -> TuneRow {
        TuneRow {
            kernel: "aggregate-row".to_string(),
            shape_class: "shape".to_string(),
            machine: 0xAC,
            params: "mode=bounded".to_string(),
            evidence,
            refresh: 1,
        }
    }

    #[test]
    fn gemm_scoped_identity_has_one_canonical_spelling() {
        let key = gemm_key();
        assert!(key.kernel().contains("/tune-v3/"));
        assert!(key.kernel().contains("/requested=4/thread-budget=4/"));
        assert!(key.kernel().contains("/memory-limit=18446744073709551615/"));
        assert!(key.kernel().contains("/build=release-opt3-cgu1-build-a"));
        assert_eq!(
            gemm_shape_from_scoped_kernel(key.kernel()),
            Some(key.shape_class())
        );
        let mut bounded_identity = key.execution().clone();
        bounded_identity.memory_limit_bytes = 1 << 20;
        let bounded = gemm_key_with(bounded_identity);
        assert_ne!(bounded.kernel(), key.kernel());
        assert_eq!(bounded.execution().memory_limit_bytes(), 1 << 20);
        assert!(bounded.kernel().contains("/memory-limit=1048576/"));
        let legacy = key.kernel().replacen("/tune-v3/", "/tune-v2/", 1);
        assert_eq!(
            gemm_shape_from_scoped_kernel(&legacy),
            None,
            "legacy keys without an explicit memory seam fail closed"
        );
        assert!(
            Tuner::cold(1)
                .pin(legacy, GemmBlockPlan::COLD_START.canonical())
                .is_err(),
            "legacy durable pins cannot bypass the current schema"
        );
        assert!(
            GemmTuneKey::new(
                format!("{GEMM_KERNEL_PREFIX}01"),
                key.shape_class(),
                key.execution().clone(),
            )
            .is_err(),
            "the bit-semantics version is canonical decimal"
        );
        assert!(
            GemmExecutionIdentity::new(
                4,
                0,
                u64::MAX,
                [64, 128, 32],
                "avx2",
                "unpinned",
                "v1",
                "release-a",
            )
            .is_err()
        );
        assert!(
            GemmExecutionIdentity::new(
                4,
                4,
                u64::MAX,
                [64, 0, 32],
                "avx2",
                "unpinned",
                "v1",
                "release-a",
            )
            .is_err()
        );
        assert!(
            GemmExecutionIdentity::new(
                4,
                4,
                u64::MAX,
                [64, 128, 32],
                "avx/2",
                "unpinned",
                "v1",
                "release-a",
            )
            .is_err()
        );
        assert!(
            GemmExecutionIdentity::new(
                4,
                4,
                u64::MAX,
                [64, 128, 32],
                "avx2",
                "unpinned",
                "v1",
                "release/a",
            )
            .is_err()
        );
    }

    #[test]
    fn selected_plan_must_equal_evidence_argmin() {
        let key = gemm_key();
        let fast = GemmBlockPlan::new(32, 2048).expect("fast plan");
        let slow = GemmBlockPlan::new(64, 1024).expect("slow plan");
        let evidence = ranked_gemm_evidence((fast, 10), (slow, 20));
        let tuner = Tuner::cold(0xA1);
        let error = tuner
            .prepare_gemm_row(&key, slow, evidence)
            .expect_err("non-argmin selection must fail");
        assert!(error.detail.contains("not evidence argmin"), "{error}");

        let tied = ranked_gemm_evidence((fast, 10), (slow, 10));
        assert!(tuner.prepare_gemm_row(&key, fast, tied.clone()).is_ok());
        assert!(
            tuner.prepare_gemm_row(&key, slow, tied).is_err(),
            "ranked evidence ties resolve to the earliest candidate"
        );
    }

    #[test]
    fn gemm_import_rejects_invalid_params() {
        let key = gemm_key();
        let selected = GemmBlockPlan::new(32, 2048).expect("selected");
        let other = GemmBlockPlan::new(64, 1024).expect("other");
        let tuner = Tuner::cold(0xA2);
        let valid = tuner
            .prepare_gemm_row(
                &key,
                selected,
                ranked_gemm_evidence((selected, 10), (other, 20)),
            )
            .expect("valid prepared row")
            .row_json();
        let invalid = valid.replace(
            "\"params\":\"mc=32,nc-cap=2048\"",
            "\"params\":\"mc=33,nc-cap=2048\"",
        );
        assert!(
            tuner.prepare_adopt_gemm_row_json(&key, &invalid).is_err(),
            "off-lattice cache params must fail closed"
        );
    }

    #[test]
    fn gemm_import_requires_ranked_evidence() {
        let key = gemm_key();
        let selected = GemmBlockPlan::new(32, 2048).expect("selected");
        let other = GemmBlockPlan::new(64, 1024).expect("other");
        let row = TuneRow {
            kernel: key.kernel().to_string(),
            shape_class: key.shape_class().to_string(),
            machine: 0xA3,
            params: selected.canonical(),
            evidence: TuneEvidence::new(vec![
                TuneObservation::wall_time(selected.canonical(), vec![10]).expect("timing"),
                TuneObservation::wall_time(other.canonical(), vec![20]).expect("timing"),
            ])
            .expect("unranked evidence"),
            refresh: 1,
        };
        let tuner = Tuner::cold(0xA3);
        assert!(
            tuner
                .prepare_adopt_gemm_row_json(&key, &row.to_json())
                .is_err(),
            "multiple timings are not ranked evidence unless explicitly declared"
        );
    }

    #[test]
    fn foreign_embedded_key_is_refused_and_remeasured() {
        let requested = gemm_key();
        let foreign = GemmTuneKey::new(
            requested.base_kernel(),
            "m128-n128-k32",
            requested.execution().clone(),
        )
        .expect("foreign key");
        let selected = GemmBlockPlan::new(32, 2048).expect("selected");
        let other = GemmBlockPlan::new(64, 1024).expect("other");
        let producer = Tuner::cold(0xA4);
        let foreign_json = producer
            .prepare_gemm_row(
                &foreign,
                selected,
                ranked_gemm_evidence((selected, 10), (other, 20)),
            )
            .expect("foreign row")
            .row_json();
        let mut consumer = Tuner::cold(0xA4);
        assert!(
            consumer
                .adopt_gemm_row_json(&requested, &foreign_json)
                .is_err(),
            "cache body must bind the exact requested key and shape"
        );
        assert!(!consumer.has_gemm_row(&requested));
        assert_eq!(
            consumer.prepare_gemm_decision(&requested).source(),
            TuneSource::ColdStart,
            "a caller must remeasure after the refused cache seed"
        );
    }

    #[test]
    fn gemm_execution_identity_dimensions_do_not_cross_reuse() {
        let baseline = gemm_key();
        let variants = [
            gemm_key_with(gemm_identity(
                8,
                4,
                [64, 128, 32],
                "avx2-fma",
                "unpinned",
                "fs-la-parallel-v1",
                "release-opt3-cgu1-build-a",
            )),
            gemm_key_with(gemm_identity(
                4,
                2,
                [64, 128, 32],
                "avx2-fma",
                "unpinned",
                "fs-la-parallel-v1",
                "release-opt3-cgu1-build-a",
            )),
            gemm_key_with(gemm_identity(
                4,
                4,
                [63, 128, 32],
                "avx2-fma",
                "unpinned",
                "fs-la-parallel-v1",
                "release-opt3-cgu1-build-a",
            )),
            gemm_key_with(gemm_identity(
                4,
                4,
                [64, 128, 32],
                "scalar",
                "unpinned",
                "fs-la-parallel-v1",
                "release-opt3-cgu1-build-a",
            )),
            gemm_key_with(gemm_identity(
                4,
                4,
                [64, 128, 32],
                "avx2-fma",
                "ccd-pinned",
                "fs-la-parallel-v1",
                "release-opt3-cgu1-build-a",
            )),
            gemm_key_with(gemm_identity(
                4,
                4,
                [64, 128, 32],
                "avx2-fma",
                "unpinned",
                "fs-la-parallel-v2",
                "release-opt3-cgu1-build-a",
            )),
            gemm_key_with(gemm_identity(
                4,
                4,
                [64, 128, 32],
                "avx2-fma",
                "unpinned",
                "fs-la-parallel-v1",
                "debug-opt0-cgu256-build-a",
            )),
        ];
        let selected = GemmBlockPlan::new(64, 1024).expect("selected");
        let other = GemmBlockPlan::new(32, 2048).expect("other");
        let mut tuner = Tuner::cold(0xA5);
        tuner
            .record_gemm_row(
                &baseline,
                selected,
                ranked_gemm_evidence((selected, 10), (other, 20)),
            )
            .expect("baseline row");
        tuner
            .pin_gemm_blocking(&baseline, selected)
            .expect("baseline pin");
        for variant in variants {
            assert_ne!(baseline.kernel(), variant.kernel());
            assert!(!tuner.has_gemm_row(&variant));
            assert!(!tuner.has_gemm_pin(&variant));
            assert_eq!(
                tuner.prepare_gemm_decision(&variant).source(),
                TuneSource::ColdStart
            );
        }
    }

    #[test]
    fn gemm_build_identity_is_stable_and_cache_rows_cannot_cross() {
        let baseline = gemm_key();
        let same_build = gemm_key();
        assert_eq!(
            baseline, same_build,
            "the same explicit build identity produces one stable key"
        );

        assert_eq!(baseline.execution().build(), "release-opt3-cgu1-build-a");

        let selected = GemmBlockPlan::new(64, 1024).expect("selected");
        let other = GemmBlockPlan::new(32, 2048).expect("other");
        let producer = Tuner::cold(0x0B01_7D1D);
        let row = producer
            .prepare_gemm_row(
                &baseline,
                selected,
                ranked_gemm_evidence((selected, 10), (other, 20)),
            )
            .expect("baseline row")
            .row_json();

        for foreign_identity in ["release-opt2-cgu1-build-a", "release-opt3-cgu1-build-b"] {
            let foreign_build = gemm_key_with(gemm_identity(
                4,
                4,
                [64, 128, 32],
                "avx2-fma",
                "unpinned",
                "fs-la-parallel-v1",
                foreign_identity,
            ));
            assert_ne!(baseline.kernel(), foreign_build.kernel());

            let mut consumer = Tuner::cold(0x0B01_7D1D);
            assert!(
                consumer.adopt_gemm_row_json(&foreign_build, &row).is_err(),
                "a durable row must not cross codegen/build identity {foreign_identity}"
            );
            assert!(!consumer.has_gemm_row(&foreign_build));
            assert_eq!(
                consumer.prepare_gemm_decision(&foreign_build).source(),
                TuneSource::ColdStart
            );
        }
    }

    #[test]
    fn recorded_gemm_key_replays_only_the_exact_scope() {
        let key = gemm_key();
        let neighbor = GemmTuneKey::new(key.base_kernel(), "m64-n256-k32", key.execution().clone())
            .expect("neighbor key");
        let plan = GemmBlockPlan::new(64, 1024).expect("plan");
        let mut live = Tuner::cold(0xA6);
        live.pin_gemm_blocking(&key, plan).expect("pin");
        let live_decision = live.prepare_gemm_decision(&key);
        assert_eq!(live_decision.plan(), plan);
        live.commit_gemm_decision(live_decision)
            .expect("successful live dispatch");
        let recorded = live.decisions()[0].clone();
        assert_eq!(recorded.kernel, key.kernel());

        let mut replay = Tuner::cold(0xFFFF);
        replay
            .pin(recorded.kernel, recorded.params)
            .expect("recorded scoped key is canonical");
        let replay_decision = replay.prepare_gemm_decision(&key);
        assert_eq!(
            (replay_decision.plan(), replay_decision.source()),
            (plan, TuneSource::Pinned)
        );
        replay
            .commit_gemm_decision(replay_decision)
            .expect("successful replay dispatch");
        assert_eq!(
            replay.prepare_gemm_decision(&neighbor).source(),
            TuneSource::ColdStart,
            "a scoped replay pin cannot leak into a neighboring shape"
        );
    }

    #[test]
    fn typed_pin_apis_reject_parameter_family_collisions() {
        let key = gemm_key();
        let mut tuner = Tuner::cold(0xA7);
        assert!(tuner.pin_tile_edge(key.kernel(), TileEdge::E4).is_err());
        assert!(tuner.pin(key.kernel(), "edge=4").is_err());
        assert!(
            tuner
                .pin(format!("{GEMM_KERNEL_PREFIX}1"), "edge=4")
                .is_err(),
            "an unscoped GEMM-family key cannot fall through to tile-edge parsing"
        );
        assert!(
            tuner
                .pin(STENCIL_KERNEL, GemmBlockPlan::COLD_START.canonical())
                .is_err()
        );
        tuner
            .pin_gemm_blocking(&key, GemmBlockPlan::COLD_START)
            .expect("typed GEMM pin");
        assert!(tuner.has_gemm_pin(&key));
    }

    #[test]
    fn prepared_gemm_rows_and_decisions_commit_explicitly() {
        let key = gemm_key();
        let selected = GemmBlockPlan::new(32, 2048).expect("selected");
        let other = GemmBlockPlan::new(64, 1024).expect("other");
        let mut tuner = Tuner::cold(0xA8);
        let evidence = ranked_gemm_evidence((selected, 10), (other, 20));
        let prepared = tuner
            .prepare_gemm_row(&key, selected, evidence.clone())
            .expect("prepared");
        let concurrently_stale = tuner
            .prepare_gemm_row(&key, selected, evidence)
            .expect("second preparation observes the same state");
        assert_eq!(prepared.params_json(), "\"mc=32,nc-cap=2048\"");
        assert!(!tuner.has_gemm_row(&key));
        tuner.commit_gemm_row(prepared).expect("commit row");
        assert!(tuner.has_gemm_row(&key));
        assert!(
            tuner.commit_gemm_row(concurrently_stale).is_err(),
            "a prepared row cannot overwrite state changed after preparation"
        );

        let decision = tuner.prepare_gemm_decision(&key);
        assert!(tuner.decisions().is_empty());
        tuner
            .commit_gemm_decision(decision)
            .expect("commit successful dispatch decision");
        assert_eq!(tuner.decisions().len(), 1);

        let neighbor =
            GemmTuneKey::new(key.base_kernel(), "m128-n128-k32", key.execution().clone())
                .expect("neighbor");
        let stale_decision = tuner.prepare_gemm_decision(&neighbor);
        tuner
            .pin_gemm_blocking(&neighbor, other)
            .expect("state changes after resolution");
        assert!(tuner.commit_gemm_decision(stale_decision).is_err());
        assert_eq!(
            tuner.decisions().len(),
            1,
            "a stale decision cannot fabricate a dispatch receipt"
        );
    }

    #[test]
    fn prepared_gemm_rows_are_bounded_writer_parser_adoption_fixed_points() {
        let key = gemm_key();
        let selected = GemmBlockPlan::new(32, 2048).expect("selected");
        let other = GemmBlockPlan::new(64, 1024).expect("other");
        let producer = Tuner::cold(0xAD);
        let prepared = producer
            .prepare_gemm_row(
                &key,
                selected,
                ranked_gemm_evidence((selected, 10), (other, 20)),
            )
            .expect("generated row");
        let json = prepared.row_json();
        assert!(json.len() <= MAX_TUNE_ROW_BYTES);
        assert_eq!(parse_row(&json).as_ref(), Some(&prepared.row));

        let mut consumer = Tuner::cold(0xAD);
        let adopted = consumer
            .prepare_adopt_gemm_row_json(&key, &json)
            .expect("writer output must be adoptable under the same key");
        assert_eq!(adopted.row_json(), json);
        consumer
            .commit_gemm_row(adopted)
            .expect("canonical adoption commit");
        assert_eq!(consumer.gemm_row_json(&key).as_deref(), Some(json.as_str()));
    }

    #[test]
    fn cold_start_defaults_are_documented_values_and_recorded() {
        let mut t = Tuner::cold(0xF1);
        assert!(t.needs_calibration());
        let (edge, src) = t.tile_edge_for("anything").expect("bounded kernel");
        assert_eq!((edge, src), (TileEdge::E8, TuneSource::ColdStart));
        let (kind, src) = t.schedule();
        assert_eq!(
            (kind, src),
            (ScheduleKind::BandwidthRich, TuneSource::ColdStart)
        );
        assert_eq!(t.decisions().len(), 2);
        assert!(t.decisions()[0].to_json().contains("cold-start"));
    }

    #[test]
    fn decision_history_is_bounded_and_eviction_is_explicit() {
        let mut tuner = Tuner::cold(0xF2);
        let calls = MAX_RETAINED_TUNING_DECISIONS + 1;
        for _ in 0..calls {
            assert_eq!(
                tuner.schedule(),
                (ScheduleKind::BandwidthRich, TuneSource::ColdStart)
            );
        }

        let history = tuner.decision_history();
        assert!(!history.is_complete());
        assert!(history.evicted() > 0);
        assert_eq!(history.total_recorded(), calls);
        assert!(history.decisions().len() <= MAX_RETAINED_TUNING_DECISIONS);
        assert!(history.retained_bytes() <= MAX_RETAINED_TUNING_DECISION_BYTES);
    }

    #[test]
    fn oversized_decision_kernel_is_refused_without_mutating_history() {
        let mut tuner = Tuner::cold(0xF3);
        let oversized = "k".repeat(MAX_TUNE_STRING_BYTES + 1);
        let error = tuner
            .tile_edge_for(&oversized)
            .expect_err("oversized decision identity must fail closed");
        assert!(error.detail.contains("kernel identity exceeds"), "{error}");
        assert!(tuner.decision_history().is_complete());
        assert!(tuner.decisions().is_empty());
    }

    #[test]
    fn schedule_selection_follows_the_measured_axis_both_ways() {
        // §5.1 consequence 2, end to end: the SAME calibration pass
        // picks a DIFFERENT schedule per machine class, driven by the
        // measured bandwidth axis — synthetic probes make the doctrine
        // testable deterministically on any host.
        let mut rich = fs_substrate::CapabilityProbe::topology_only();
        rich.logical_cpus = 14;
        rich.measured.all_core_gbs = 273.0; // M4 Pro class: ~19.5 GB/s/core
        let mut t_rich = Tuner::cold(rich.fingerprint());
        t_rich.calibrate(&rich).expect("matching rich probe");
        assert_eq!(t_rich.schedule().0, ScheduleKind::BandwidthRich);

        let mut starved = fs_substrate::CapabilityProbe::topology_only();
        starved.logical_cpus = 128;
        starved.measured.all_core_gbs = 120.0; // 5995WX class: ~0.94 GB/s/core
        let mut t_starved = Tuner::cold(starved.fingerprint());
        t_starved
            .calibrate(&starved)
            .expect("matching starved probe");
        assert_eq!(t_starved.schedule().0, ScheduleKind::BandwidthStarved);

        // The per-worker work-count row landed with the sorted distribution.
        let row = t_starved
            .rows
            .get(&(CLASS_WORK_KERNEL.to_string(), SHAPE_DEFAULT.to_string()))
            .expect("class row recorded");
        assert!(row.params.contains("slow_over_fast_x1000="));
    }

    #[test]
    fn pins_beat_everything_and_replay_reproducibly() {
        let mut t = Tuner::cold(0xF1);
        t.pin(STENCIL_KERNEL, "edge=16").expect("canonical pin");
        let (edge, src) = t.tile_edge_for(STENCIL_KERNEL).expect("bounded kernel");
        assert_eq!((edge, src), (TileEdge::E16, TuneSource::Pinned));
        // Replay: a second tuner pinned from the recorded decision agrees.
        let recorded = t.decisions()[0].clone();
        let mut replay = Tuner::cold(0xDEAD_BEEF); // different machine!
        replay
            .pin(recorded.kernel.clone(), recorded.params.clone())
            .expect("recorded pin is canonical");
        let (edge2, src2) = replay
            .tile_edge_for(STENCIL_KERNEL)
            .expect("bounded kernel");
        assert_eq!((edge2, src2), (TileEdge::E16, TuneSource::Pinned));
    }

    #[test]
    fn invalid_pins_fail_closed_instead_of_becoming_pinned_defaults() {
        let mut t = Tuner::cold(0xF1);
        assert!(t.pin(" ", "edge=4").is_err());
        assert!(t.pin_tile_edge("", TileEdge::E4).is_err());
        assert!(t.pin(STENCIL_KERNEL, "edge=32").is_err());
        assert!(t.pin(SCHEDULE_KERNEL, "schedule=rich").is_err());
        assert!(t.pin_tile_edge(SCHEDULE_KERNEL, TileEdge::E4).is_err());

        let (edge, source) = t.tile_edge_for(STENCIL_KERNEL).expect("bounded kernel");
        assert_eq!((edge, source), (TileEdge::E8, TuneSource::ColdStart));
        t.pin_tile_edge(STENCIL_KERNEL, TileEdge::E4)
            .expect("typed tile pin");
        assert_eq!(
            t.tile_edge_for(STENCIL_KERNEL).expect("bounded kernel"),
            (TileEdge::E4, TuneSource::Pinned)
        );
        t.pin_schedule(ScheduleKind::BandwidthStarved);
        assert_eq!(
            t.schedule(),
            (ScheduleKind::BandwidthStarved, TuneSource::Pinned)
        );
    }

    #[test]
    fn calibration_rejects_a_foreign_probe_before_mutating_rows() {
        let probe = fs_substrate::CapabilityProbe::topology_only();
        let mut t = Tuner::cold(probe.fingerprint() ^ 1);
        let err = t.calibrate(&probe).expect_err("foreign probe must fail");
        assert!(err.detail.contains("does not match"), "{err}");
        assert!(t.needs_calibration(), "failed calibration added no rows");
        assert!(t.decisions().is_empty());
    }

    #[test]
    fn row_json_round_trips_and_foreign_fingerprints_are_stale() {
        let dir = std::env::temp_dir().join("fs-exec-tune-test");
        std::fs::create_dir_all(&dir).expect("tmp dir");
        let path = dir.join("tune-roundtrip.jsonl");
        let mut t = Tuner::cold(0xAB);
        t.insert_row(
            "k1",
            "s1",
            "edge=4".to_string(),
            TuneEvidence::new(vec![
                TuneObservation::wall_time("edge=4", vec![100, 110]).expect("timing"),
                TuneObservation::wall_time("edge=8", vec![200, 210]).expect("timing"),
            ])
            .expect("evidence"),
        )
        .expect("first row insert");
        t.save(&path).expect("save");
        let same = Tuner::load(&path, 0xAB).expect("load");
        assert_eq!(
            same.rows.values().next().expect("row").params,
            "edge=4",
            "round trip"
        );
        let other = Tuner::load(&path, 0xCD).expect("load other machine");
        assert!(
            other.needs_calibration(),
            "foreign-fingerprint rows are stale and ignored"
        );
        let missing = Tuner::load(&dir.join("nope.jsonl"), 0xAB).expect("missing file");
        assert!(missing.needs_calibration(), "missing store is a cold start");
        let err = {
            std::fs::write(&path, "not a row\n").expect("write garbage");
            Tuner::load(&path, 0xAB).expect_err("corrupt store")
        };
        assert!(err.to_string().contains("recalibrate"), "{err}");
    }

    #[test]
    fn loader_rejects_duplicate_rows_for_the_selected_fingerprint() {
        let dir = std::env::temp_dir().join("fs-exec-tune-duplicate-test");
        std::fs::create_dir_all(&dir).expect("tmp dir");
        let path = dir.join("duplicate.jsonl");
        let row = TuneRow {
            kernel: "kernel".to_string(),
            shape_class: "shape".to_string(),
            machine: 0xAB,
            params: "edge=4".to_string(),
            evidence: TuneEvidence::new(vec![
                TuneObservation::wall_time("candidate", vec![1]).expect("timing"),
            ])
            .expect("evidence"),
            refresh: 1,
        }
        .to_json();
        std::fs::write(&path, format!("{row}\n{row}\n")).expect("duplicate store");
        let error = Tuner::load(&path, 0xAB).expect_err("duplicate selected rows must fail");
        assert!(error.detail.contains("duplicate tune row"), "{error}");

        let foreign = Tuner::load(&path, 0xCD).expect("foreign duplicate rows are stale");
        assert!(foreign.needs_calibration());
    }

    #[test]
    fn candidate_separation_is_descriptive_and_absent_for_a_single_observation() {
        let singleton = TuneEvidence::new(vec![
            TuneObservation::wall_time("a", vec![100]).expect("timing"),
        ])
        .expect("evidence");
        assert_eq!(
            singleton.candidate_separation_ppm(),
            None,
            "one observation cannot imply confidence or candidate separation"
        );

        let near_tie = TuneEvidence::ranked_wall_times(vec![
            TuneObservation::wall_time("a", vec![100]).expect("timing"),
            TuneObservation::wall_time("b", vec![101]).expect("timing"),
        ])
        .expect("evidence");
        let clear = TuneEvidence::ranked_wall_times(vec![
            TuneObservation::wall_time("a", vec![100]).expect("timing"),
            TuneObservation::wall_time("b", vec![400]).expect("timing"),
        ])
        .expect("evidence");
        assert_eq!(near_tie.candidate_separation_ppm(), Some(9_900));
        assert_eq!(clear.candidate_separation_ppm(), Some(750_000));
        let ranked_row = TuneRow {
            kernel: "ranked".to_string(),
            shape_class: "shape".to_string(),
            machine: 1,
            params: "winner=a".to_string(),
            evidence: clear.clone(),
            refresh: 1,
        }
        .to_json();
        assert!(parse_row(&ranked_row).is_some());
        assert!(
            parse_row(&ranked_row.replace(
                "\"candidate_separation_ppm\":750000",
                "\"candidate_separation_ppm\":749999",
            ))
            .is_none(),
            "stored candidate separation must be re-derived"
        );

        let unranked = TuneEvidence::new(vec![
            TuneObservation::wall_time("phase-a", vec![100]).expect("timing"),
            TuneObservation::wall_time("phase-b", vec![400]).expect("timing"),
        ])
        .expect("unranked evidence");
        assert_eq!(
            unranked.candidate_separation_ppm(),
            None,
            "multiple durations are not candidates unless declared as such"
        );

        let mixed = TuneEvidence::new(vec![
            TuneObservation::wall_time("run", vec![100]).expect("timing"),
            TuneObservation::work_count("steals", WorkUnit::Steals, 7).expect("counter"),
        ])
        .expect("mixed evidence");
        assert_eq!(mixed.candidate_separation_ppm(), None);
        assert!(TuneEvidence::ranked_wall_times(mixed.observations().to_vec()).is_err());
    }

    #[test]
    fn row_parser_rejects_noncanonical_or_out_of_domain_fields() {
        let row = TuneRow {
            kernel: STENCIL_KERNEL.to_string(),
            shape_class: SHAPE_DEFAULT.to_string(),
            machine: 0xAB,
            params: "edge=4".to_string(),
            evidence: TuneEvidence::new(vec![
                TuneObservation::wall_time("edge=4", vec![100]).expect("timing"),
            ])
            .expect("evidence"),
            refresh: 1,
        }
        .to_json();
        assert!(parse_row(&row).is_some(), "writer output must parse");

        for corrupt in [
            format!("{row} trailing"),
            row.replace("\"version\":1", "\"version\":2"),
            row.replace("\"minimum_ns\":100", "\"minimum_ns\":99"),
            row.replace("\"maximum_ns\":100", "\"maximum_ns\":101"),
            row.replace(
                "\"candidate_separation_ppm\":null",
                "\"candidate_separation_ppm\":1000000",
            ),
            row.replace("\"samples_ns\":[100]", "\"samples_ns\":[]"),
            row.replace("\"samples_ns\":[100]", "\"samples_ns\":[0100]"),
            row.replace("\"kind\":\"wall-time\"", "\"kind\":\"nanoseconds\""),
            row.replace("\"refresh\":1", "\"refresh\":1.5"),
            row.replace("\"refresh\":1", "\"refresh\":-1"),
            row.replace("\"refresh\":1", "\"refresh\":0"),
            row.replace("\"refresh\":1", "\"refresh\":4294967296"),
            row.replace("\"params\":\"edge=4\"", "\"params\":\"edge=32\""),
        ] {
            assert!(
                parse_row(&corrupt).is_none(),
                "accepted corrupt row: {corrupt}"
            );
        }

        let escaped = TuneRow {
            kernel: "custom\"kernel".to_string(),
            shape_class: "line\nbreak".to_string(),
            machine: 0xCD,
            params: "edge=4".to_string(),
            evidence: TuneEvidence::new(vec![
                TuneObservation::wall_time("candidate\\path", vec![7]).expect("timing"),
            ])
            .expect("evidence"),
            refresh: 1,
        }
        .to_json();
        let parsed = parse_row(&escaped).expect("canonical escaped strings round trip");
        assert_eq!(parsed.kernel, "custom\"kernel");
        assert_eq!(parsed.shape_class, "line\nbreak");
        assert_eq!(parsed.evidence.observations()[0].label(), "candidate\\path");
    }

    #[test]
    fn row_parser_preserves_full_width_typed_integers() {
        let row = TuneRow {
            kernel: "typed-evidence".to_string(),
            shape_class: "full-width".to_string(),
            machine: u64::MAX,
            params: "mode=probe".to_string(),
            evidence: TuneEvidence::new(vec![
                TuneObservation::wall_time("wall", vec![u64::MAX]).expect("timing"),
                TuneObservation::work_count("work", WorkUnit::CompletedTiles, u64::MAX)
                    .expect("work"),
                TuneObservation::throughput(
                    "bandwidth",
                    ThroughputUnit::GigabytesPerSecond,
                    u64::MAX,
                )
                .expect("throughput"),
            ])
            .expect("evidence"),
            refresh: u32::MAX,
        };
        let json = row.to_json();
        let parsed = parse_row(&json).expect("canonical full-width row");
        assert_eq!(parsed, row);
        assert!(json.matches("18446744073709551615").count() >= 4);

        for corrupt in [
            json.replace("\"unit\":\"completed-tiles\"", "\"unit\":\"tiles\""),
            json.replace(
                "\"unit\":\"gigabytes-per-second\"",
                "\"unit\":\"bytes-per-second\"",
            ),
            json.replace("18446744073709551615", "18446744073709551616"),
        ] {
            assert!(
                parse_row(&corrupt).is_none(),
                "accepted corrupt row: {corrupt}"
            );
        }
    }

    #[test]
    fn evidence_rejects_blank_duplicate_and_inconsistent_observations() {
        assert!(TuneObservation::wall_time(" ", vec![1]).is_err());
        assert!(TuneObservation::wall_time("empty", Vec::new()).is_err());
        assert!(TuneObservation::wall_time("zero", vec![0]).is_err());
        assert!(TuneObservation::wall_time("mixed-zero", vec![10, 0, 20]).is_err());
        let duplicate = TuneEvidence::new(vec![
            TuneObservation::wall_time("same", vec![1]).expect("timing"),
            TuneObservation::work_count("same", WorkUnit::Steals, 2).expect("work"),
        ]);
        assert!(duplicate.is_err());

        let inconsistent = TuneObservation::WallTime {
            label: "candidate".to_string(),
            samples_ns: vec![10, 20],
            summary: WallTimeSummary {
                minimum_ns: 9,
                maximum_ns: 20,
            },
        };
        assert!(TuneEvidence::new(vec![inconsistent]).is_err());
    }

    #[test]
    fn canonical_ranked_row_adoption_rejects_zero_wall_time() {
        let row = TuneRow {
            kernel: "hostile-ranked-kernel".to_string(),
            shape_class: "shape".to_string(),
            machine: 0xAB,
            params: "candidate-a".to_string(),
            evidence: TuneEvidence::ranked_wall_times(vec![
                TuneObservation::wall_time("candidate-a", vec![1]).expect("first timing"),
                TuneObservation::wall_time("candidate-b", vec![2]).expect("second timing"),
            ])
            .expect("ranked evidence"),
            refresh: 1,
        }
        .to_json();
        let hostile = row.replace(
            "\"samples_ns\":[1],\"summary\":{\"minimum_ns\":1,\"maximum_ns\":1}",
            "\"samples_ns\":[0],\"summary\":{\"minimum_ns\":0,\"maximum_ns\":0}",
        );
        assert_ne!(
            hostile, row,
            "hostile fixture must replace the fastest sample"
        );
        let hostile = hostile.replace(
            "\"candidate_separation_ppm\":500000",
            "\"candidate_separation_ppm\":1000000",
        );
        assert!(
            hostile.contains("\"candidate_separation_ppm\":1000000"),
            "hostile fixture must carry the separation derived from zero and two"
        );
        assert!(
            parse_row(&hostile).is_none(),
            "a structurally canonical zero-duration row must fail closed"
        );

        let mut tuner = Tuner::cold(0xAB);
        assert!(tuner.adopt_row_json(&hostile).is_err());
        assert!(
            tuner.row_json("hostile-ranked-kernel", "shape").is_none(),
            "failed adoption must not mutate tuner state"
        );
    }

    #[test]
    fn numeric_conversions_fail_closed_without_partial_calibration() {
        assert_eq!(throughput_milli_units(1.234).expect("scaled"), 1234);
        assert_eq!(
            duration_ns(std::time::Duration::ZERO, "test").expect("measurement floor"),
            1
        );
        for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, -0.001, f64::MAX] {
            assert!(
                throughput_milli_units(invalid).is_err(),
                "accepted invalid throughput {invalid:?}"
            );
        }
        let enormous = std::time::Duration::new(u64::MAX, 999_999_999);
        assert!(duration_ns(enormous, "test").is_err());

        let mut probe = fs_substrate::CapabilityProbe::topology_only();
        probe.measured.all_core_gbs = f64::NAN;
        let mut tuner = Tuner::cold(probe.fingerprint());
        assert!(tuner.calibrate(&probe).is_err());
        assert!(
            tuner.needs_calibration(),
            "invalid numeric evidence must fail before any row mutation"
        );

        let mut unattributed = fs_substrate::CapabilityProbe::topology_only();
        unattributed.logical_cpus = 0;
        unattributed.measured.all_core_gbs = 1.0;
        let mut tuner = Tuner::cold(unattributed.fingerprint());
        assert!(tuner.calibrate(&unattributed).is_err());
        assert!(tuner.needs_calibration());

        let mut oversized = fs_substrate::CapabilityProbe::topology_only();
        oversized.logical_cpus = (MAX_TUNE_OBSERVATIONS + 1) as u32;
        let mut tuner = Tuner::cold(oversized.fingerprint());
        assert!(tuner.calibrate(&oversized).is_err());
        assert!(tuner.needs_calibration());
    }

    #[test]
    fn evidence_and_parser_resource_limits_fail_closed() {
        assert!(
            TuneObservation::wall_time("oversampled", vec![1; MAX_WALL_TIME_SAMPLES + 1]).is_err()
        );
        let observation =
            TuneObservation::work_count("work", WorkUnit::Steals, 1).expect("bounded observation");
        assert!(TuneEvidence::new(vec![observation; MAX_TUNE_OBSERVATIONS + 1]).is_err());
        assert!(
            TuneObservation::work_count(
                "x".repeat(MAX_TUNE_STRING_BYTES + 1),
                WorkUnit::Steals,
                1,
            )
            .is_err()
        );
        assert!(parse_row(&"x".repeat(MAX_TUNE_ROW_BYTES + 1)).is_none());

        let bounded = TuneRow {
            kernel: "kernel".to_string(),
            shape_class: "shape".to_string(),
            machine: 1,
            params: "edge=4".to_string(),
            evidence: TuneEvidence::new(vec![
                TuneObservation::wall_time("candidate", vec![1]).expect("timing"),
            ])
            .expect("evidence"),
            refresh: 1,
        }
        .to_json();
        let mut empty_observations = bounded.clone();
        let start_marker = "\"observations\":[";
        let start =
            empty_observations.find(start_marker).expect("observations") + start_marker.len();
        let end = start
            + empty_observations[start..]
                .find("],\"candidate_separation_ppm\"")
                .expect("observation end");
        empty_observations.replace_range(start..end, "");
        assert!(parse_row(&empty_observations).is_none());

        let too_many_samples = bounded.replace(
            "\"samples_ns\":[1]",
            &format!(
                "\"samples_ns\":[{}]",
                vec!["1"; MAX_WALL_TIME_SAMPLES + 1].join(",")
            ),
        );
        assert!(parse_row(&too_many_samples).is_none());
    }

    #[test]
    fn evidence_enforces_the_aggregate_wall_time_sample_budget() {
        let first = MAX_TUNE_ROW_WALL_TIME_SAMPLES / 2;
        let second = MAX_TUNE_ROW_WALL_TIME_SAMPLES - first;
        let boundary = TuneEvidence::new(vec![
            TuneObservation::wall_time("first", vec![1; first]).expect("first samples"),
            TuneObservation::wall_time("second", vec![2; second]).expect("second samples"),
        ])
        .expect("the exact aggregate boundary is accepted");
        assert_eq!(
            boundary
                .observations()
                .iter()
                .map(|observation| match observation {
                    TuneObservation::WallTime { samples_ns, .. } => samples_ns.len(),
                    TuneObservation::WorkCount { .. } | TuneObservation::Throughput { .. } => 0,
                })
                .sum::<usize>(),
            MAX_TUNE_ROW_WALL_TIME_SAMPLES
        );

        let error = TuneEvidence::new(vec![
            TuneObservation::wall_time("first", vec![1; first]).expect("first samples"),
            TuneObservation::wall_time("second", vec![2; second + 1]).expect("second samples"),
        ])
        .expect_err("one sample beyond the aggregate budget must fail");
        assert!(
            error.detail.contains("aggregate wall-time limit"),
            "{error}"
        );
    }

    #[test]
    fn generated_rows_accept_the_exact_byte_boundary_and_refuse_one_more_byte() {
        let mut label_lengths = vec![MAX_TUNE_STRING_BYTES; 15];
        let shortest_last = "work-15:".len();
        label_lengths.push(shortest_last);
        let shortest_row = work_count_row(work_count_evidence_with_label_lengths(&label_lengths));
        let remaining = MAX_TUNE_ROW_BYTES
            .checked_sub(shortest_row.to_json().len())
            .expect("fifteen maximum labels leave room for the final observation");
        let exact_last = shortest_last + remaining;
        assert!(
            exact_last < MAX_TUNE_STRING_BYTES,
            "the one-byte overflow fixture must remain inside the per-string limit"
        );
        *label_lengths.last_mut().expect("last label") = exact_last;

        let boundary_evidence = work_count_evidence_with_label_lengths(&label_lengths);
        let boundary_row = work_count_row(boundary_evidence.clone());
        let boundary_json = boundary_row
            .validated_generated_json()
            .expect("the exact canonical row boundary is accepted");
        assert_eq!(boundary_json.len(), MAX_TUNE_ROW_BYTES);
        assert_eq!(parse_row(&boundary_json), Some(boundary_row));

        let mut tuner = Tuner::cold(0xAC);
        tuner
            .insert_row(
                "aggregate-row",
                "shape",
                "mode=bounded".to_string(),
                boundary_evidence,
            )
            .expect("the generated-row insertion boundary is accepted");
        assert_eq!(
            tuner
                .row_json("aggregate-row", "shape")
                .expect("inserted boundary row")
                .len(),
            MAX_TUNE_ROW_BYTES
        );

        *label_lengths.last_mut().expect("last label") += 1;
        let oversized_evidence = work_count_evidence_with_label_lengths(&label_lengths);
        let mut rejected = Tuner::cold(0xAC);
        let error = rejected
            .insert_row(
                "aggregate-row",
                "shape",
                "mode=bounded".to_string(),
                oversized_evidence.clone(),
            )
            .expect_err("one byte beyond the canonical row limit must fail");
        assert!(error.detail.contains("canonical row limit"), "{error}");
        assert!(rejected.row_json("aggregate-row", "shape").is_none());

        let mut poisoned = Tuner::cold(0xAC);
        let oversized_row = work_count_row(oversized_evidence);
        poisoned.rows.insert(
            ("aggregate-row".to_string(), "shape".to_string()),
            oversized_row,
        );
        let path = std::env::temp_dir().join(format!(
            "fs-exec-oversized-generated-row-{}.jsonl",
            std::process::id()
        ));
        let error = poisoned
            .save(&path)
            .expect_err("persistence revalidates even internally injected rows");
        assert!(error.detail.contains("canonical row limit"), "{error}");
    }

    #[test]
    fn parser_bounds_observation_and_decoded_string_growth() {
        let row = TuneRow {
            kernel: "kernel".to_string(),
            shape_class: "shape".to_string(),
            machine: 1,
            params: "edge=4".to_string(),
            evidence: TuneEvidence::new(vec![
                TuneObservation::work_count("work", WorkUnit::Steals, 1).expect("work"),
            ])
            .expect("evidence"),
            refresh: 1,
        }
        .to_json();
        let oversized_string = row.replace(
            "\"kernel\":\"kernel\"",
            &format!("\"kernel\":\"{}\"", "k".repeat(MAX_TUNE_STRING_BYTES + 1)),
        );
        assert!(parse_row(&oversized_string).is_none());

        let observations = (0..=MAX_TUNE_OBSERVATIONS)
            .map(|index| {
                format!(
                    "{{\"kind\":\"work-count\",\"label\":\"w{index}\",\"unit\":\"steals\",\"count\":1}}"
                )
            })
            .collect::<Vec<_>>()
            .join(",");
        let start_marker = "\"observations\":[";
        let start = row.find(start_marker).expect("observations") + start_marker.len();
        let end = start
            + row[start..]
                .find("],\"candidate_separation_ppm\"")
                .expect("observation end");
        let mut too_many_observations = row;
        too_many_observations.replace_range(start..end, &observations);
        assert!(parse_row(&too_many_observations).is_none());
    }

    #[test]
    fn loader_bounds_total_store_before_parsing() {
        let dir = std::env::temp_dir().join("fs-exec-tune-store-bound-test");
        std::fs::create_dir_all(&dir).expect("tmp dir");
        let path = dir.join("oversized.jsonl");
        std::fs::write(&path, vec![b'x'; MAX_TUNE_STORE_BYTES + 1]).expect("oversized store");
        let error = Tuner::load(&path, 0).expect_err("oversized store must fail");
        assert!(error.detail.contains("exceeds"), "{error}");
    }

    #[test]
    fn exhausted_refresh_counter_refuses_before_calibration() {
        let mut tuner = Tuner::cold(0xAB);
        tuner.rows.insert(
            (STENCIL_KERNEL.to_string(), SHAPE_DEFAULT.to_string()),
            TuneRow {
                kernel: STENCIL_KERNEL.to_string(),
                shape_class: SHAPE_DEFAULT.to_string(),
                machine: 0xAB,
                params: "edge=4".to_string(),
                evidence: TuneEvidence::new(vec![
                    TuneObservation::wall_time("edge=4", vec![1]).expect("timing"),
                ])
                .expect("evidence"),
                refresh: u32::MAX,
            },
        );
        let err = tuner
            .require_refresh_capacity()
            .expect_err("counter exhaustion must fail closed");
        assert!(err.detail.contains("counter exhausted"), "{err}");
    }
}
