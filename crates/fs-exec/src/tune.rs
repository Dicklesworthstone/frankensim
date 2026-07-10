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
    /// Returns [`TuneError`] for a blank label or an empty sample set.
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
    /// summaries, blank labels, or duplicate labels.
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
    /// present and every observation is wall time.
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

    fn push_json(&self, out: &mut String) {
        let _ = write!(out, "{{\"version\":{},\"observations\":[", Self::VERSION);
        for (index, observation) in self.observations.iter().enumerate() {
            if index > 0 {
                out.push(',');
            }
            match observation {
                TuneObservation::WallTime {
                    label,
                    samples_ns,
                    summary,
                } => {
                    out.push_str("{\"kind\":\"wall-time\",\"label\":");
                    push_json_string(out, label);
                    out.push_str(",\"samples_ns\":[");
                    for (sample_index, sample) in samples_ns.iter().enumerate() {
                        if sample_index > 0 {
                            out.push(',');
                        }
                        let _ = write!(out, "{sample}");
                    }
                    let _ = write!(
                        out,
                        "],\"summary\":{{\"minimum_ns\":{},\"maximum_ns\":{}}}}}",
                        summary.minimum_ns, summary.maximum_ns
                    );
                }
                TuneObservation::WorkCount { label, unit, count } => {
                    out.push_str("{\"kind\":\"work-count\",\"label\":");
                    push_json_string(out, label);
                    out.push_str(",\"unit\":");
                    push_json_string(out, unit.name());
                    let _ = write!(out, ",\"count\":{count}}}");
                }
                TuneObservation::Throughput {
                    label,
                    unit,
                    milli_units,
                } => {
                    out.push_str("{\"kind\":\"throughput\",\"label\":");
                    push_json_string(out, label);
                    out.push_str(",\"unit\":");
                    push_json_string(out, unit.name());
                    let _ = write!(out, ",\"milli_units\":{milli_units}}}");
                }
            }
        }
        out.push_str("],\"candidate_separation_ppm\":");
        if let Some(separation) = self.candidate_separation_ppm {
            let _ = write!(out, "{separation}");
        } else {
            out.push_str("null");
        }
        out.push('}');
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
    for observation in observations {
        observation.validate()?;
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
    /// Canonical JSON-line (deterministic field order).
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut s = String::with_capacity(160);
        s.push_str("{\"kernel\":");
        push_json_string(&mut s, &self.kernel);
        s.push_str(",\"shape_class\":");
        push_json_string(&mut s, &self.shape_class);
        let _ = write!(s, ",\"machine\":\"{:016x}\",\"params\":", self.machine);
        push_json_string(&mut s, &self.params);
        s.push_str(",\"evidence\":");
        self.evidence.push_json(&mut s);
        let _ = write!(s, ",\"refresh\":{}}}", self.refresh);
        s
    }
}

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
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out.push('"');
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
        push_json_string(&mut out, &self.kernel);
        out.push_str(",\"params\":");
        push_json_string(&mut out, &self.params);
        out.push_str(",\"source\":");
        push_json_string(&mut out, self.source);
        out.push('}');
        out
    }
}

/// The per-machine tuner: tune rows + pins + decision log.
#[derive(Debug)]
pub struct Tuner {
    fingerprint: u64,
    rows: BTreeMap<(String, String), TuneRow>,
    pins: BTreeMap<String, PinnedParam>,
    decisions: Vec<TuningDecision>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PinnedParam {
    TileEdge(TileEdge),
    Schedule(ScheduleKind),
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
        if kernel.trim().is_empty() {
            return Err(TuneError {
                detail: "a pinned kernel identity must be non-blank".to_string(),
            });
        }
        if kernel == SCHEDULE_KERNEL {
            return Err(TuneError {
                detail: "the reserved schedule key requires pin_schedule".to_string(),
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

    /// Validate and pin canonical params from a recorded decision.
    ///
    /// The accepted replay forms are exactly `edge={4,8,16}` for tile
    /// kernels and `schedule={bandwidth-rich,bandwidth-starved}` for the
    /// schedule key. Invalid spellings fail closed instead of resolving to a
    /// default that would be falsely recorded as pinned.
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
        if kernel.trim().is_empty() {
            return Err(TuneError {
                detail: "a pinned kernel identity must be non-blank".to_string(),
            });
        }
        let parsed = PinnedParam::parse(&kernel, &params).ok_or_else(|| TuneError {
            detail: format!("invalid pinned params {params:?} for kernel {kernel:?}"),
        })?;
        self.pins.insert(kernel, parsed);
        Ok(())
    }

    /// The decisions handed out so far (a study records these; replaying
    /// the study pins them).
    #[must_use]
    pub fn decisions(&self) -> &[TuningDecision] {
        &self.decisions
    }

    /// Resolve a tile edge for `kernel`: pins beat tuned rows beat the
    /// cold-start default (8³ — plan §5.3). The decision is recorded.
    pub fn tile_edge_for(&mut self, kernel: &str) -> (TileEdge, TuneSource) {
        let (params, source) = self.resolve(kernel);
        let edge = match params.as_str() {
            "edge=4" => TileEdge::E4,
            "edge=16" => TileEdge::E16,
            _ => TileEdge::E8,
        };
        (edge, source)
    }

    /// Resolve the schedule kind: pins beat tuned beat the cold-start
    /// default (bandwidth-rich — harmless on starved machines, just less
    /// blocked). The decision is recorded.
    pub fn schedule(&mut self) -> (ScheduleKind, TuneSource) {
        let (params, source) = self.resolve(SCHEDULE_KERNEL);
        let kind = if params == "schedule=bandwidth-starved" {
            ScheduleKind::BandwidthStarved
        } else {
            ScheduleKind::BandwidthRich
        };
        (kind, source)
    }

    fn resolve(&mut self, kernel: &str) -> (String, TuneSource) {
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
        self.decisions.push(TuningDecision {
            kernel: kernel.to_string(),
            params: params.clone(),
            source: source.name(),
        });
        (params, source)
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
        self.rows.insert(
            key,
            TuneRow {
                kernel: kernel.to_string(),
                shape_class: shape_class.to_string(),
                machine: self.fingerprint,
                params,
                evidence,
                refresh,
            },
        );
        Ok(())
    }

    /// Persist the table (JSON-lines, one row per line — the ledger `tune`
    /// schema's file edition).
    ///
    /// # Errors
    /// [`TuneError`] on I/O failure.
    pub fn save(&self, path: &std::path::Path) -> Result<(), TuneError> {
        let mut out = String::new();
        for row in self.rows.values() {
            out.push_str(&row.to_json());
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
    u64::try_from(duration.as_nanos()).map_err(|_| TuneError {
        detail: format!("{context} duration exceeds the tune evidence u64 nanosecond range"),
    })
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

    #[test]
    fn cold_start_defaults_are_documented_values_and_recorded() {
        let mut t = Tuner::cold(0xF1);
        assert!(t.needs_calibration());
        let (edge, src) = t.tile_edge_for("anything");
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
        let (edge, src) = t.tile_edge_for(STENCIL_KERNEL);
        assert_eq!((edge, src), (TileEdge::E16, TuneSource::Pinned));
        // Replay: a second tuner pinned from the recorded decision agrees.
        let recorded = t.decisions()[0].clone();
        let mut replay = Tuner::cold(0xDEAD_BEEF); // different machine!
        replay
            .pin(recorded.kernel.clone(), recorded.params.clone())
            .expect("recorded pin is canonical");
        let (edge2, src2) = replay.tile_edge_for(STENCIL_KERNEL);
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

        let (edge, source) = t.tile_edge_for(STENCIL_KERNEL);
        assert_eq!((edge, source), (TileEdge::E8, TuneSource::ColdStart));
        t.pin_tile_edge(STENCIL_KERNEL, TileEdge::E4)
            .expect("typed tile pin");
        assert_eq!(
            t.tile_edge_for(STENCIL_KERNEL),
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
    fn numeric_conversions_fail_closed_without_partial_calibration() {
        assert_eq!(throughput_milli_units(1.234).expect("scaled"), 1234);
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
