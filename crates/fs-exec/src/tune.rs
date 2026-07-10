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
//! Persistence: a JSON-lines file store whose row shape mirrors the ledger
//! `tune` table (kernel × shape-class × machine fingerprint → params +
//! measured + confidence); migrating to fs-ledger is a rename, not a
//! rewrite. Rows keyed to a DIFFERENT fingerprint are stale by definition
//! and ignored on load.

use crate::cx::{CancelGate, Cx};
use crate::kernel::{TileKernel, TilePlan};
use crate::pool::{PoolConfig, TilePool};
use core::fmt;
use core::ops::ControlFlow;
use fs_substrate::affinity::CcdTopology;
use fs_substrate::tile::TileEdge;
use std::collections::BTreeMap;
use std::fmt::Write as _;

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

/// One tune-table row (the ledger `tune` schema, file-store edition).
#[derive(Debug, Clone, PartialEq)]
pub struct TuneRow {
    /// Kernel identity (e.g. "stencil7-f32").
    pub kernel: String,
    /// Shape class (e.g. "48c-cube").
    pub shape_class: String,
    /// Machine fingerprint the measurements belong to.
    pub machine: u64,
    /// Chosen parameter, canonical form (e.g. "edge=8", "schedule=bandwidth-rich").
    pub params: String,
    /// Measured evidence: `(candidate, best-of-repeats ns)` pairs.
    pub measured_ns: Vec<(String, u64)>,
    /// Agreement between repeats, 0..=1 (1 = perfectly stable).
    pub confidence: f64,
    /// Recalibration counter (idempotence witness).
    pub refresh: u32,
}

impl TuneRow {
    /// Canonical JSON-line (deterministic field order).
    #[must_use]
    pub fn to_json(&self) -> String {
        let mut s = String::with_capacity(160);
        let _ = write!(
            s,
            "{{\"kernel\":\"{}\",\"shape_class\":\"{}\",\"machine\":\"{:016x}\",\
             \"params\":\"{}\",\"measured_ns\":[",
            self.kernel, self.shape_class, self.machine, self.params
        );
        for (i, (name, ns)) in self.measured_ns.iter().enumerate() {
            if i > 0 {
                s.push(',');
            }
            let _ = write!(s, "{{\"candidate\":\"{name}\",\"ns\":{ns}}}");
        }
        let _ = write!(
            s,
            "],\"confidence\":{:.3},\"refresh\":{}}}",
            self.confidence, self.refresh
        );
        s
    }
}

/// Structured tune-store failure (Decalogue P10).
#[derive(Debug)]
pub struct TuneError {
    /// What went wrong, with the path.
    pub detail: String,
}

impl fmt::Display for TuneError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "tune store error: {}; calibration data is reproducible — delete the store and \
             recalibrate if it is corrupt",
            self.detail
        )
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
        format!(
            "{{\"kernel\":\"{}\",\"params\":\"{}\",\"source\":\"{}\"}}",
            self.kernel, self.params, self.source
        )
    }
}

/// The per-machine tuner: tune rows + pins + decision log.
#[derive(Debug)]
pub struct Tuner {
    fingerprint: u64,
    rows: BTreeMap<(String, String), TuneRow>,
    pins: BTreeMap<String, String>,
    decisions: Vec<TuningDecision>,
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

    /// Pin a kernel's params (replay path: recorded plans, never re-tuned).
    pub fn pin(&mut self, kernel: impl Into<String>, params: impl Into<String>) {
        self.pins.insert(kernel.into(), params.into());
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
            (p.clone(), TuneSource::Pinned)
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
    pub fn calibrate(&mut self, probe: &fs_substrate::CapabilityProbe) -> String {
        let workers = (probe.logical_cpus as usize).clamp(1, 8);
        let topo = CcdTopology::from_probe(probe);
        // Tile-edge sweep: best-of-2 per edge, argmin wins; ties break to
        // the LOWER candidate index (deterministic tie law).
        let mut measured: Vec<(String, u64)> = Vec::new();
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
            let mut best = u64::MAX;
            for _ in 0..2 {
                let gate = CancelGate::new();
                let t0 = std::time::Instant::now();
                let (r, _) = pool.run_with_gate(&kernel, &gate);
                let ns = t0.elapsed().as_nanos() as u64;
                r.expect("calibration stencil runs");
                best = best.min(ns);
            }
            measured.push((format!("edge={}", edge.cells()), best));
        }
        let winner = measured
            .iter()
            .enumerate()
            .min_by_key(|(i, (_, ns))| (*ns, *i))
            .map(|(_, (name, _))| name.clone())
            .expect("three candidates");
        let confidence = confidence_from(&measured);
        self.insert_row(STENCIL_KERNEL, SHAPE_DEFAULT, winner, measured);

        // Reduction cost: deterministic compensated sum over 100k terms.
        let xs: Vec<f64> = (0..100_000).map(|i| 1.0 / f64::from(i + 1)).collect();
        let t0 = std::time::Instant::now();
        let s = crate::reduce::det_sum(&xs);
        let red_ns = t0.elapsed().as_nanos() as u64;
        assert!(s.is_finite());
        self.insert_row(
            "det-sum-f64",
            "100k",
            "block=256".to_string(),
            vec![("block=256".to_string(), red_ns)],
        );

        // Steal cost: an imbalanced 2-worker run forces steals; record the
        // per-steal wall cost estimate.
        let pool = TilePool::new(PoolConfig {
            quantum_weights: vec![15, 1],
            ..PoolConfig::new(2, topo, 0x7C4E)
        });
        let gate = CancelGate::new();
        let t0 = std::time::Instant::now();
        let (r, report) = pool.run_with_gate(&StealProbe, &gate);
        let steal_ns = t0.elapsed().as_nanos() as u64;
        r.expect("steal probe runs");
        self.insert_row(
            "steal-probe",
            "2w-imbalanced",
            format!("steals={}", report.steals),
            vec![("wall".to_string(), steal_ns)],
        );

        // Per-class throughput (fz2.2): the measured class signal.
        self.class_throughput_row(probe, topo);

        // Schedule: measured bandwidth per logical core decides (the §5.1
        // consequence-2 doctrine); zero measurement -> rich default,
        // recorded as such.
        let per_core = if probe.logical_cpus > 0 {
            probe.measured.all_core_gbs / f64::from(probe.logical_cpus)
        } else {
            0.0
        };
        let kind = if probe.measured.all_core_gbs > 0.0 && per_core < 8.0 {
            ScheduleKind::BandwidthStarved
        } else {
            ScheduleKind::BandwidthRich
        };
        self.insert_row(
            SCHEDULE_KERNEL,
            SHAPE_DEFAULT,
            format!("schedule={}", kind.name()),
            vec![("per-core-gbs-x1000".to_string(), (per_core * 1000.0) as u64)],
        );

        let mut s = String::with_capacity(512);
        let _ = write!(
            s,
            "{{\"machine\":\"{:016x}\",\"confidence\":{confidence:.3},\"rows\":[",
            self.fingerprint
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

    /// Per-class throughput row (fz2.2): FULL parallelism, uniform
    /// tiles; per-worker completed counts are the measured class
    /// signal. Records the sorted distribution and the fast:slow
    /// ratio — the quantum_weights VERIFICATION (weighted initial
    /// shares vs pure stealing) lives in the machine A/B lane.
    fn class_throughput_row(&mut self, probe: &fs_substrate::CapabilityProbe, topo: CcdTopology) {
        let full = (probe.logical_cpus as usize).max(1);
        let pool = TilePool::new(PoolConfig::new(full, topo, 0x7C4E));
        let gate = CancelGate::new();
        let (r, report) = pool.run_with_gate(&ClassProbe, &gate);
        r.expect("class probe runs");
        let mut counts = report.tiles_by_worker.clone();
        counts.sort_unstable_by(|a, b| b.cmp(a));
        let fast = counts.first().copied().unwrap_or(0).max(1);
        let slow = counts.last().copied().unwrap_or(0).max(1);
        #[allow(clippy::cast_precision_loss)]
        let ratio_x1000 = (slow as f64 / fast as f64 * 1000.0) as u64;
        self.insert_row(
            "class-throughput",
            SHAPE_DEFAULT,
            format!("workers={full};slow_over_fast_x1000={ratio_x1000}"),
            counts
                .iter()
                .enumerate()
                .map(|(i, &c)| (format!("w{i}"), c))
                .collect(),
        );
    }

    fn insert_row(
        &mut self,
        kernel: &str,
        shape_class: &str,
        params: String,
        measured_ns: Vec<(String, u64)>,
    ) {
        let key = (kernel.to_string(), shape_class.to_string());
        let refresh = self.rows.get(&key).map_or(1, |r| r.refresh + 1);
        let confidence = confidence_from(&measured_ns);
        self.rows.insert(
            key,
            TuneRow {
                kernel: kernel.to_string(),
                shape_class: shape_class.to_string(),
                machine: self.fingerprint,
                params,
                measured_ns,
                confidence,
                refresh,
            },
        );
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
        let text = match std::fs::read_to_string(path) {
            Ok(t) => t,
            Err(e) if e.kind() == std::io::ErrorKind::NotFound => return Ok(tuner),
            Err(e) => {
                return Err(TuneError {
                    detail: format!("reading {}: {e}", path.display()),
                });
            }
        };
        for (lineno, line) in text.lines().enumerate() {
            let row = parse_row(line).ok_or_else(|| TuneError {
                detail: format!("corrupt row at {}:{}", path.display(), lineno + 1),
            })?;
            if row.machine == fingerprint {
                tuner
                    .rows
                    .insert((row.kernel.clone(), row.shape_class.clone()), row);
            }
        }
        Ok(tuner)
    }
}

/// Uniform CPU-bound tiles for the per-class throughput probe (fz2.2):
/// every tile costs the same FLOP budget, so per-worker completed-tile
/// counts under work-stealing ARE the measured per-class throughput
/// (heterogeneous slow-class workers complete measurably fewer).
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

/// Repeat-agreement confidence: 1 - spread/best over candidates' best
/// times, clamped to [0, 1] (1 candidate → 1.0 by convention).
fn confidence_from(measured: &[(String, u64)]) -> f64 {
    let best = measured.iter().map(|&(_, ns)| ns).min().unwrap_or(1).max(1);
    let worst = measured.iter().map(|&(_, ns)| ns).max().unwrap_or(1);
    if measured.len() < 2 {
        1.0
    } else {
        // Wide spreads mean the choice mattered AND the ranking is clear;
        // near-ties mean low decision confidence.
        (1.0 - (best as f64) / (worst as f64)).clamp(0.0, 1.0)
    }
}

/// Minimal strict parser for our own row writer (deviation = corruption,
/// same doctrine as fs-obs).
fn parse_row(line: &str) -> Option<TuneRow> {
    let field = |key: &str| -> Option<String> {
        let pat = format!("\"{key}\":\"");
        let start = line.find(&pat)? + pat.len();
        let end = line[start..].find('"')? + start;
        Some(line[start..end].to_string())
    };
    let num_field = |key: &str| -> Option<f64> {
        let pat = format!("\"{key}\":");
        let start = line.find(&pat)? + pat.len();
        let end = line[start..]
            .find([',', '}'])
            .map_or(line.len(), |e| e + start);
        line[start..end].trim().parse().ok()
    };
    let mut measured_ns = Vec::new();
    let mut rest = line;
    while let Some(p) = rest.find("{\"candidate\":\"") {
        let after = &rest[p + 14..];
        let name_end = after.find('"')?;
        let name = after[..name_end].to_string();
        let ns_pat = "\"ns\":";
        let ns_start = after.find(ns_pat)? + ns_pat.len();
        let ns_end = after[ns_start..].find('}')? + ns_start;
        let ns: u64 = after[ns_start..ns_end].trim().parse().ok()?;
        measured_ns.push((name, ns));
        rest = &after[ns_end..];
    }
    Some(TuneRow {
        kernel: field("kernel")?,
        shape_class: field("shape_class")?,
        machine: u64::from_str_radix(&field("machine")?, 16).ok()?,
        params: field("params")?,
        measured_ns,
        confidence: num_field("confidence")?,
        refresh: num_field("refresh")? as u32,
    })
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
        let mut t_rich = Tuner::cold(0xA);
        t_rich.calibrate(&rich);
        assert_eq!(t_rich.schedule().0, ScheduleKind::BandwidthRich);

        let mut starved = fs_substrate::CapabilityProbe::topology_only();
        starved.logical_cpus = 128;
        starved.measured.all_core_gbs = 120.0; // 5995WX class: ~0.94 GB/s/core
        let mut t_starved = Tuner::cold(0xB);
        t_starved.calibrate(&starved);
        assert_eq!(t_starved.schedule().0, ScheduleKind::BandwidthStarved);

        // The class-throughput row landed with the sorted distribution.
        let row = t_starved
            .rows
            .get(&("class-throughput".to_string(), SHAPE_DEFAULT.to_string()))
            .expect("class row recorded");
        assert!(row.params.contains("slow_over_fast_x1000="));
    }

    #[test]
    fn pins_beat_everything_and_replay_reproducibly() {
        let mut t = Tuner::cold(0xF1);
        t.pin(STENCIL_KERNEL, "edge=16");
        let (edge, src) = t.tile_edge_for(STENCIL_KERNEL);
        assert_eq!((edge, src), (TileEdge::E16, TuneSource::Pinned));
        // Replay: a second tuner pinned from the recorded decision agrees.
        let recorded = t.decisions()[0].clone();
        let mut replay = Tuner::cold(0xDEAD_BEEF); // different machine!
        replay.pin(recorded.kernel.clone(), recorded.params.clone());
        let (edge2, src2) = replay.tile_edge_for(STENCIL_KERNEL);
        assert_eq!((edge2, src2), (TileEdge::E16, TuneSource::Pinned));
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
            vec![("edge=4".to_string(), 100), ("edge=8".to_string(), 200)],
        );
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
    fn confidence_reflects_candidate_spread() {
        assert!((confidence_from(&[("a".into(), 100)]) - 1.0).abs() < 1e-12);
        let near_tie = confidence_from(&[("a".into(), 100), ("b".into(), 101)]);
        let clear = confidence_from(&[("a".into(), 100), ("b".into(), 400)]);
        assert!(near_tie < 0.05, "near ties are low confidence: {near_tie}");
        assert!(clear > 0.7, "clear rankings are high confidence: {clear}");
    }
}
