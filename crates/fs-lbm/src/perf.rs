//! D3Q19 sparse-sweep performance evidence model (bead 712t).
//!
//! This module is the dependency-neutral semantic core beneath the ignored
//! measurement lane. It fixes the logical byte/FLOP model, 128³ dense-active
//! and ten-percent-active workload shapes, raw plan targets, evidence classes,
//! gate semantics, and deterministic max-plus attribution before timing code
//! can emit rows. Keeping this algebra in production code makes the expensive
//! runner a thin observer rather than an independent source of truth.
//!
//! No-claims: the model performs no timing, hardware probing, admission, or
//! baseline authentication. A report-only row is never citable. The raw plan
//! target is informational until both reference-ISA campaigns retain admitted
//! measurements and a separately authorized anti-collapse floor.

use core::fmt;
use std::collections::{BTreeMap, BTreeSet};

/// Version of the D3Q19 sparse-sweep traffic and receipt semantics.
pub const D3Q19_PERF_MODEL_VERSION: &str = "d3q19-sparse-sweep-v2-local-halo";
/// Parts per million used for exact occupancy/share ratios.
pub const RATIO_PPM: u32 = 1_000_000;
/// Edge length of the current sparse tile.
pub const SPARSE_TILE_EDGE: usize = 4;
/// Cells in one 4×4×4 tile.
pub const SPARSE_TILE_CELLS: usize = SPARSE_TILE_EDGE * SPARSE_TILE_EDGE * SPARSE_TILE_EDGE;
/// D3Q19 population count.
pub const D3Q19_DISTRIBUTIONS: usize = 19;
/// Total distribution links visited by one 4³ tile update.
pub const D3Q19_LINKS_PER_TILE: usize = D3Q19_DISTRIBUTIONS * SPARSE_TILE_CELLS;
/// Same-tile links served without a sparse key/slot lookup.
///
/// D3Q19 contributes 1,216 links per 4³ tile. The six face directions cross
/// one 4×4 plane each (96 links), and the twelve edge directions cross the
/// union of two signed planes (12 × (16 + 16 - 4) = 336 links). The remaining
/// 784 links are local.
pub const D3Q19_LOCAL_LINKS_PER_TILE: usize = 784;
/// Cross-tile/domain links that enter the sparse key/slot lookup path.
pub const D3Q19_HALO_LINKS_PER_TILE: usize = 432;
/// Bytes in one f64 distribution value.
pub const DISTRIBUTION_BYTES: usize = core::mem::size_of::<f64>();
/// BGK density/momentum reductions plus forced-velocity reconstruction.
pub const BGK_MACRO_VELOCITY_FLOPS: u16 = 142;
/// D3Q19 equilibrium construction including `u²` and all 19 directions.
pub const BGK_EQUILIBRIUM_FLOPS: u16 = 271;
/// Generic three-axis Guo projection plus BGK relaxation for all directions.
pub const BGK_FORCE_RELAX_FLOPS: u16 = 611;
/// Source-counted FLOPs per active BGK lattice update.
pub const BGK_FLOPS_PER_CELL: u16 =
    BGK_MACRO_VELOCITY_FLOPS + BGK_EQUILIBRIUM_FLOPS + BGK_FORCE_RELAX_FLOPS;
/// Maximum timed tasks retained in one attribution receipt.
pub const MAX_CRITICAL_TASKS: usize = 65_536;
/// Maximum precedence edges retained in one attribution receipt.
pub const MAX_CRITICAL_EDGES: usize = 262_144;
/// Maximum timed repetitions retained in one GLUP/s row.
pub const MAX_PERF_REPETITIONS: usize = 64;
const MAX_PREDECESSORS_PER_TASK: usize = 64;

/// Stable refusal from workload, receipt, or max-plus validation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PerfModelError {
    /// Dimensions are zero or not divisible by the frozen tile edge.
    InvalidDimensions([usize; 3]),
    /// Thread count is zero or contradicts the selected lane.
    InvalidWorkers {
        /// Selected threading lane.
        threading: ThreadingClass,
        /// Offered worker count.
        workers: usize,
    },
    /// Checked shape/byte/time arithmetic overflowed.
    ArithmeticOverflow(&'static str),
    /// A collection exceeded a declared resource cap.
    ResourceLimit {
        /// Bounded collection.
        resource: &'static str,
        /// Configured cap.
        limit: usize,
        /// Offered length.
        observed: usize,
    },
    /// Task identity appeared more than once.
    DuplicateTask(u32),
    /// One task repeated the same predecessor.
    DuplicatePredecessor {
        /// Task containing the duplicate.
        task: u32,
        /// Repeated predecessor.
        predecessor: u32,
    },
    /// A precedence edge names an absent task.
    MissingPredecessor {
        /// Task owning the edge.
        task: u32,
        /// Missing predecessor id.
        predecessor: u32,
    },
    /// A task has a zero wall sample.
    ZeroWallSample(u32),
    /// Precedence edges contain a cycle.
    CyclicTaskGraph,
    /// A row field is malformed or contradicts its evidence class.
    InvalidReceipt(&'static str),
}

impl fmt::Display for PerfModelError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidDimensions(dims) => write!(
                f,
                "D3Q19 perf dims {dims:?} must be positive multiples of {SPARSE_TILE_EDGE}"
            ),
            Self::InvalidWorkers { threading, workers } => write!(
                f,
                "{threading:?} D3Q19 perf lane cannot use {workers} workers"
            ),
            Self::ArithmeticOverflow(operation) => {
                write!(f, "D3Q19 perf arithmetic overflow during {operation}")
            }
            Self::ResourceLimit {
                resource,
                limit,
                observed,
            } => write!(
                f,
                "D3Q19 perf {resource} cap {limit} exceeded by {observed}"
            ),
            Self::DuplicateTask(task) => write!(f, "duplicate timing task id {task}"),
            Self::DuplicatePredecessor { task, predecessor } => {
                write!(f, "timing task {task} repeats predecessor {predecessor}")
            }
            Self::MissingPredecessor { task, predecessor } => write!(
                f,
                "timing task {task} names absent predecessor {predecessor}"
            ),
            Self::ZeroWallSample(task) => write!(f, "timing task {task} has zero wall_ns"),
            Self::CyclicTaskGraph => f.write_str("timing task graph contains a cycle"),
            Self::InvalidReceipt(field) => {
                write!(f, "D3Q19 performance receipt has invalid {field}")
            }
        }
    }
}

impl std::error::Error for PerfModelError {}

/// Reference target family from the plan.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ReferenceIsa {
    /// Apple M-class aarch64 host (1.0 GLUP/s plan target).
    AppleMClass,
    /// Threadripper/EPYC-class x86-64 host (0.6 GLUP/s plan target).
    ThreadripperClass,
    /// Any host outside the two target families; no plan target applies.
    Other,
}

impl ReferenceIsa {
    /// Informational raw-throughput target. It is never the initial floor.
    #[must_use]
    pub const fn plan_target_glups(self) -> Option<f64> {
        match self {
            Self::AppleMClass => Some(1.0),
            Self::ThreadripperClass => Some(0.6),
            Self::Other => None,
        }
    }

    /// Stable receipt identity.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::AppleMClass => "apple-m-class",
            Self::ThreadripperClass => "threadripper-class",
            Self::Other => "other",
        }
    }
}

/// Frozen occupancy shape for the memory-resident sweep.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum OccupancyClass {
    /// Every tile in the 128³ domain is active.
    DenseActive,
    /// Smallest whole-tile set covering at least ten percent of the domain.
    SparseTenPercent,
}

impl OccupancyClass {
    /// Exact requested active fraction in parts per million.
    #[must_use]
    pub const fn active_fraction_ppm(self) -> u32 {
        match self {
            Self::DenseActive => RATIO_PPM,
            Self::SparseTenPercent => 100_000,
        }
    }

    /// Stable receipt identity.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::DenseActive => "dense-active",
            Self::SparseTenPercent => "sparse-ten-percent",
        }
    }
}

/// Threading axis of one throughput row.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ThreadingClass {
    /// Serial reference sweep with exactly one worker.
    SingleThread,
    /// Pooled sweep against aggregate machine axes.
    AllCore,
}

impl ThreadingClass {
    /// Stable receipt identity.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::SingleThread => "single-thread",
            Self::AllCore => "all-core",
        }
    }
}

/// Exact workload geometry of one D3Q19 measurement row.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LaneShape {
    /// Cell dimensions, each a multiple of four.
    pub dims: [usize; 3],
    /// Active-tile occupancy class.
    pub occupancy: OccupancyClass,
    /// Single-thread or all-core lane.
    pub threading: ThreadingClass,
    /// Exact worker count retained in the row.
    pub workers: usize,
}

impl LaneShape {
    /// Canonical 128³ memory-resident shape.
    pub fn memory_resident(
        occupancy: OccupancyClass,
        threading: ThreadingClass,
        workers: usize,
    ) -> Result<Self, PerfModelError> {
        let shape = Self {
            dims: [128; 3],
            occupancy,
            threading,
            workers,
        };
        shape.validate()?;
        Ok(shape)
    }

    /// Validate tiling and worker semantics.
    pub fn validate(self) -> Result<(), PerfModelError> {
        if self
            .dims
            .iter()
            .any(|&dim| dim == 0 || dim % SPARSE_TILE_EDGE != 0)
        {
            return Err(PerfModelError::InvalidDimensions(self.dims));
        }
        if self.workers == 0
            || (self.threading == ThreadingClass::SingleThread && self.workers != 1)
        {
            return Err(PerfModelError::InvalidWorkers {
                threading: self.threading,
                workers: self.workers,
            });
        }
        Ok(())
    }

    /// Total whole tiles in the declared domain.
    pub fn total_tiles(self) -> Result<usize, PerfModelError> {
        self.validate()?;
        self.dims
            .iter()
            .try_fold(1usize, |product, dim| {
                product.checked_mul(*dim / SPARSE_TILE_EDGE)
            })
            .ok_or(PerfModelError::ArithmeticOverflow("total tiles"))
    }

    /// Active whole tiles. Sparse occupancy rounds upward so the retained
    /// workload never silently undershoots ten percent.
    pub fn active_tiles(self) -> Result<usize, PerfModelError> {
        let tiles = self.total_tiles()?;
        let numerator = tiles
            .checked_mul(self.occupancy.active_fraction_ppm() as usize)
            .ok_or(PerfModelError::ArithmeticOverflow("active tile numerator"))?;
        numerator
            .checked_add(RATIO_PPM as usize - 1)
            .map(|rounded| rounded / RATIO_PPM as usize)
            .ok_or(PerfModelError::ArithmeticOverflow("active tile rounding"))
    }

    /// Active lattice cells represented by the whole-tile workload.
    pub fn active_cells(self) -> Result<usize, PerfModelError> {
        self.active_tiles()?
            .checked_mul(SPARSE_TILE_CELLS)
            .ok_or(PerfModelError::ArithmeticOverflow("active cells"))
    }

    /// Three population buffers (published, collided, transactional next).
    pub fn allocated_population_bytes(self) -> Result<usize, PerfModelError> {
        self.active_cells()?
            .checked_mul(D3Q19_DISTRIBUTIONS)
            .and_then(|values| values.checked_mul(DISTRIBUTION_BYTES))
            .and_then(|one_buffer| one_buffer.checked_mul(3))
            .ok_or(PerfModelError::ArithmeticOverflow(
                "allocated population bytes",
            ))
    }
}

/// Versioned logical traffic model for one two-pass sparse sweep.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct D3q19TrafficModel {
    /// Distribution loads per cell: collide input plus stream source.
    pub distribution_reads: u8,
    /// Distribution stores per cell: collided plus transactional output.
    pub distribution_writes: u8,
    /// Logical key/slot bytes amortized per active tile.
    pub tile_metadata_bytes: u16,
    /// Logical sparse source-lookup bytes per cross-tile/domain link.
    pub lookup_bytes_per_halo_link: u8,
    /// Versioned scalar FLOP count for collision/equilibrium/macros.
    pub flops_per_cell: u16,
}

impl Default for D3q19TrafficModel {
    fn default() -> Self {
        Self {
            distribution_reads: 2,
            distribution_writes: 2,
            tile_metadata_bytes: 16,
            lookup_bytes_per_halo_link: 16,
            flops_per_cell: BGK_FLOPS_PER_CELL,
        }
    }
}

impl D3q19TrafficModel {
    /// Population payload traffic (19 f64 values × two reads × two writes).
    #[must_use]
    pub fn population_bytes_per_cell(self) -> f64 {
        f64::from(self.distribution_reads + self.distribution_writes)
            * D3Q19_DISTRIBUTIONS as f64
            * DISTRIBUTION_BYTES as f64
    }

    /// Logical sparse-index traffic for cross-tile/domain links, including
    /// amortized active-tile key/slot identity. Same-tile pulls use direct
    /// lane addressing and are deliberately excluded.
    #[must_use]
    pub fn sparse_overhead_bytes_per_cell(self) -> f64 {
        f64::from(self.lookup_bytes_per_halo_link) * D3Q19_HALO_LINKS_PER_TILE as f64
            / SPARSE_TILE_CELLS as f64
            + f64::from(self.tile_metadata_bytes) / SPARSE_TILE_CELLS as f64
    }

    /// Total versioned logical bytes per active lattice update.
    #[must_use]
    pub fn bytes_per_cell(self) -> f64 {
        self.population_bytes_per_cell() + self.sparse_overhead_bytes_per_cell()
    }

    /// FLOP/byte arithmetic intensity used to identify the binding roof.
    #[must_use]
    pub fn arithmetic_intensity(self) -> f64 {
        f64::from(self.flops_per_cell) / self.bytes_per_cell()
    }

    /// Canonical JSON fragment for lane headers and retained receipts.
    #[must_use]
    pub fn receipt_json(self) -> String {
        format!(
            "{{\"version\":\"{D3Q19_PERF_MODEL_VERSION}\",\"q\":{D3Q19_DISTRIBUTIONS},\"distribution_bytes\":{DISTRIBUTION_BYTES},\"distribution_reads\":{},\"distribution_writes\":{},\"tile_metadata_bytes\":{},\"lookup_bytes_per_halo_link\":{},\"links_per_tile\":{D3Q19_LINKS_PER_TILE},\"local_links_per_tile\":{D3Q19_LOCAL_LINKS_PER_TILE},\"halo_links_per_tile\":{D3Q19_HALO_LINKS_PER_TILE},\"population_bytes_per_cell\":{:.6},\"sparse_overhead_bytes_per_cell\":{:.6},\"bytes_per_cell\":{:.6},\"flops_per_cell\":{},\"flops_per_byte\":{:.9}}}",
            self.distribution_reads,
            self.distribution_writes,
            self.tile_metadata_bytes,
            self.lookup_bytes_per_halo_link,
            self.population_bytes_per_cell(),
            self.sparse_overhead_bytes_per_cell(),
            self.bytes_per_cell(),
            self.flops_per_cell,
            self.arithmetic_intensity(),
        )
    }
}

/// Timed kernel class used for bottleneck attribution.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum KernelClass {
    /// Activation/deactivation and active-set maintenance.
    Activation = 0,
    /// Per-cell collision/equilibrium work.
    Collide = 1,
    /// Inter-tile boundary/halo exchange or lookup.
    Halo = 2,
    /// Pull-stream population movement.
    Stream = 3,
}

impl KernelClass {
    /// Stable receipt identity.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Activation => "activation",
            Self::Collide => "collide",
            Self::Halo => "halo",
            Self::Stream => "stream",
        }
    }

    const fn index(self) -> usize {
        self as usize
    }
}

/// One bounded wall sample and its task-DAG predecessors.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TaskSample {
    /// Stable task id within one repetition.
    pub id: u32,
    /// Kernel class executed by the task.
    pub class: KernelClass,
    /// Positive wall duration in nanoseconds.
    pub wall_ns: u64,
    /// Predecessor task ids; input order is not semantic.
    pub predecessors: Vec<u32>,
}

/// Deterministic max-plus critical-path attribution for one repetition.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CriticalPathAttribution {
    /// Source-to-sink task ids on the selected critical path.
    pub path: Vec<u32>,
    /// Max-plus makespan in nanoseconds.
    pub makespan_ns: u64,
    /// Per-class wall totals along the selected path in enum order.
    pub class_wall_ns: [u64; 4],
    /// Strictly deterministic dominant class (enum order breaks exact ties).
    pub dominant_class: KernelClass,
    /// Dominant-class share of critical-path wall in parts per million.
    pub dominant_share_ppm: u32,
}

impl CriticalPathAttribution {
    /// Validate internal max-plus totals and deterministic dominant-class math.
    pub fn validate(&self) -> Result<(), PerfModelError> {
        if self.path.is_empty() || self.makespan_ns == 0 {
            return Err(PerfModelError::InvalidReceipt("critical path"));
        }
        let unique: BTreeSet<u32> = self.path.iter().copied().collect();
        if unique.len() != self.path.len() {
            return Err(PerfModelError::InvalidReceipt("critical path cycle"));
        }
        let total = self
            .class_wall_ns
            .iter()
            .try_fold(0u64, |sum, wall| sum.checked_add(*wall));
        if total != Some(self.makespan_ns) {
            return Err(PerfModelError::InvalidReceipt("critical path class totals"));
        }
        let expected_class = [
            KernelClass::Activation,
            KernelClass::Collide,
            KernelClass::Halo,
            KernelClass::Stream,
        ]
        .into_iter()
        .max_by(|left, right| {
            self.class_wall_ns[left.index()]
                .cmp(&self.class_wall_ns[right.index()])
                .then_with(|| right.cmp(left))
        })
        .ok_or(PerfModelError::InvalidReceipt("kernel class inventory"))?;
        let share = (u128::from(self.class_wall_ns[expected_class.index()])
            * u128::from(RATIO_PPM)
            + u128::from(self.makespan_ns / 2))
            / u128::from(self.makespan_ns);
        if self.dominant_class != expected_class || u128::from(self.dominant_share_ppm) != share {
            return Err(PerfModelError::InvalidReceipt(
                "critical path dominant class/share",
            ));
        }
        Ok(())
    }

    /// Stable signature whose equality across repetitions is required in
    /// deterministic mode. Durations may vary; topology/class must not.
    #[must_use]
    pub fn stable_signature(&self) -> (&[u32], KernelClass) {
        (&self.path, self.dominant_class)
    }

    /// Canonical JSON fragment retained by measurement rows.
    #[must_use]
    pub fn receipt_json(&self) -> String {
        let path = self
            .path
            .iter()
            .map(u32::to_string)
            .collect::<Vec<_>>()
            .join(",");
        format!(
            "{{\"path\":[{path}],\"makespan_ns\":{},\"class_wall_ns\":{{\"activation\":{},\"collide\":{},\"halo\":{},\"stream\":{}}},\"dominant_class\":\"{}\",\"dominant_share_ppm\":{}}}",
            self.makespan_ns,
            self.class_wall_ns[KernelClass::Activation.index()],
            self.class_wall_ns[KernelClass::Collide.index()],
            self.class_wall_ns[KernelClass::Halo.index()],
            self.class_wall_ns[KernelClass::Stream.index()],
            self.dominant_class.as_str(),
            self.dominant_share_ppm,
        )
    }
}

/// Analyze a bounded task DAG in the max-plus semiring. Input enumeration and
/// predecessor enumeration do not affect the result; smaller task ids break
/// exact path ties and enum order breaks exact class-total ties.
pub fn attribute_critical_path(
    samples: &[TaskSample],
) -> Result<CriticalPathAttribution, PerfModelError> {
    if samples.is_empty() {
        return Err(PerfModelError::InvalidReceipt("empty task graph"));
    }
    if samples.len() > MAX_CRITICAL_TASKS {
        return Err(PerfModelError::ResourceLimit {
            resource: "critical tasks",
            limit: MAX_CRITICAL_TASKS,
            observed: samples.len(),
        });
    }
    let mut tasks = samples.to_vec();
    tasks.sort_by_key(|task| task.id);
    let mut index = BTreeMap::new();
    for (slot, task) in tasks.iter().enumerate() {
        if task.wall_ns == 0 {
            return Err(PerfModelError::ZeroWallSample(task.id));
        }
        if index.insert(task.id, slot).is_some() {
            return Err(PerfModelError::DuplicateTask(task.id));
        }
        if task.predecessors.len() > MAX_PREDECESSORS_PER_TASK {
            return Err(PerfModelError::ResourceLimit {
                resource: "predecessors per task",
                limit: MAX_PREDECESSORS_PER_TASK,
                observed: task.predecessors.len(),
            });
        }
    }
    let edge_count = tasks.iter().try_fold(0usize, |count, task| {
        count.checked_add(task.predecessors.len())
    });
    let Some(edge_count) = edge_count else {
        return Err(PerfModelError::ArithmeticOverflow("edge count"));
    };
    if edge_count > MAX_CRITICAL_EDGES {
        return Err(PerfModelError::ResourceLimit {
            resource: "critical edges",
            limit: MAX_CRITICAL_EDGES,
            observed: edge_count,
        });
    }

    let mut predecessors = vec![Vec::new(); tasks.len()];
    let mut successors = vec![Vec::new(); tasks.len()];
    let mut indegree = vec![0usize; tasks.len()];
    for (slot, task) in tasks.iter().enumerate() {
        let mut unique = BTreeSet::new();
        for &predecessor in &task.predecessors {
            if !unique.insert(predecessor) {
                return Err(PerfModelError::DuplicatePredecessor {
                    task: task.id,
                    predecessor,
                });
            }
            let Some(&predecessor_slot) = index.get(&predecessor) else {
                return Err(PerfModelError::MissingPredecessor {
                    task: task.id,
                    predecessor,
                });
            };
            predecessors[slot].push(predecessor_slot);
            successors[predecessor_slot].push(slot);
            indegree[slot] += 1;
        }
    }
    for next in &mut successors {
        next.sort_by_key(|&slot| tasks[slot].id);
    }

    let mut ready: BTreeSet<(u32, usize)> = tasks
        .iter()
        .enumerate()
        .filter(|(slot, _)| indegree[*slot] == 0)
        .map(|(slot, task)| (task.id, slot))
        .collect();
    let mut finish = vec![0u64; tasks.len()];
    let mut best_predecessor = vec![None; tasks.len()];
    let mut visited = 0usize;
    while let Some(&(id, slot)) = ready.first() {
        ready.remove(&(id, slot));
        visited += 1;
        let best = predecessors[slot].iter().copied().max_by(|&left, &right| {
            finish[left]
                .cmp(&finish[right])
                .then_with(|| tasks[right].id.cmp(&tasks[left].id))
        });
        let predecessor_finish = best.map_or(0, |predecessor| finish[predecessor]);
        finish[slot] = predecessor_finish.checked_add(tasks[slot].wall_ns).ok_or(
            PerfModelError::ArithmeticOverflow("critical path accumulation"),
        )?;
        best_predecessor[slot] = best;
        for &next in &successors[slot] {
            indegree[next] -= 1;
            if indegree[next] == 0 {
                ready.insert((tasks[next].id, next));
            }
        }
    }
    if visited != tasks.len() {
        return Err(PerfModelError::CyclicTaskGraph);
    }
    let sink = (0..tasks.len())
        .max_by(|&left, &right| {
            finish[left]
                .cmp(&finish[right])
                .then_with(|| tasks[right].id.cmp(&tasks[left].id))
        })
        .ok_or(PerfModelError::InvalidReceipt("empty task graph"))?;
    let mut path_slots = vec![sink];
    let mut current = sink;
    while let Some(predecessor) = best_predecessor[current] {
        path_slots.push(predecessor);
        current = predecessor;
    }
    path_slots.reverse();

    let mut class_wall_ns = [0u64; 4];
    for &slot in &path_slots {
        let bucket = &mut class_wall_ns[tasks[slot].class.index()];
        *bucket = bucket
            .checked_add(tasks[slot].wall_ns)
            .ok_or(PerfModelError::ArithmeticOverflow("class wall total"))?;
    }
    let dominant_class = [
        KernelClass::Activation,
        KernelClass::Collide,
        KernelClass::Halo,
        KernelClass::Stream,
    ]
    .into_iter()
    .max_by(|left, right| {
        class_wall_ns[left.index()]
            .cmp(&class_wall_ns[right.index()])
            .then_with(|| right.cmp(left))
    })
    .ok_or(PerfModelError::InvalidReceipt("kernel class inventory"))?;
    let makespan_ns = finish[sink];
    let share = (u128::from(class_wall_ns[dominant_class.index()]) * u128::from(RATIO_PPM)
        + u128::from(makespan_ns / 2))
        / u128::from(makespan_ns);
    let dominant_share_ppm =
        u32::try_from(share).map_err(|_| PerfModelError::ArithmeticOverflow("dominant share"))?;
    Ok(CriticalPathAttribution {
        path: path_slots.iter().map(|&slot| tasks[slot].id).collect(),
        makespan_ns,
        class_wall_ns,
        dominant_class,
        dominant_share_ppm,
    })
}

/// Whether all repetitions retain the same deterministic path/class signature.
#[must_use]
pub fn critical_path_is_stable(repetitions: &[CriticalPathAttribution]) -> bool {
    let Some(first) = repetitions.first() else {
        return false;
    };
    repetitions
        .iter()
        .skip(1)
        .all(|row| row.stable_signature() == first.stable_signature())
}

/// Evidence authority of one measured row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EvidenceClass {
    /// Attested baseline/admission and durable receipt permit a positive gate.
    Citable {
        /// Lowercase 32-byte admission receipt digest.
        admission_receipt: String,
    },
    /// Measurement is retained but cannot support a positive/negative gate.
    ReportOnly {
        /// Exact configuration/admission refusal.
        reason: String,
    },
    /// Host/axes contamination invalidates the measurement population.
    EnvironmentInvalid {
        /// Exact environment-invalid reason.
        reason: String,
    },
}

impl EvidenceClass {
    /// Stable evidence-class name.
    #[must_use]
    pub const fn as_str(&self) -> &'static str {
        match self {
            Self::Citable { .. } => "citable",
            Self::ReportOnly { .. } => "report_only",
            Self::EnvironmentInvalid { .. } => "environment_invalid",
        }
    }

    fn validate(&self) -> Result<(), PerfModelError> {
        match self {
            Self::Citable { admission_receipt } => {
                if admission_receipt.len() != 64
                    || !admission_receipt
                        .bytes()
                        .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
                {
                    return Err(PerfModelError::InvalidReceipt("citable admission receipt"));
                }
            }
            Self::ReportOnly { reason } | Self::EnvironmentInvalid { reason } => {
                if reason.trim().is_empty() || reason.len() > 4096 {
                    return Err(PerfModelError::InvalidReceipt("evidence reason"));
                }
            }
        }
        Ok(())
    }
}

/// Gate result. The plan target is reported independently and never becomes
/// the initial anti-collapse floor by implication.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PerfGateVerdict {
    /// Admitted row met its separately authorized floor.
    FloorMet,
    /// Admitted row missed its separately authorized floor.
    FloorMiss,
    /// Row is measured but not citable.
    ReportOnly,
    /// Environment invalidates the row rather than passing/failing it.
    EnvironmentInvalid,
}

impl PerfGateVerdict {
    /// Stable receipt identity.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::FloorMet => "floor_met",
            Self::FloorMiss => "floor_miss",
            Self::ReportOnly => "report_only",
            Self::EnvironmentInvalid => "environment_invalid",
        }
    }
}

/// One fully typed D3Q19 throughput row emitted by the measurement driver.
#[derive(Debug, Clone, PartialEq)]
pub struct D3q19PerfRow {
    /// Reference target family.
    pub reference_isa: ReferenceIsa,
    /// Exact workload/thread shape.
    pub shape: LaneShape,
    /// Measured billions of active lattice updates per second.
    pub glups: f64,
    /// Interquartile timing dispersion in parts per million.
    pub dispersion_ppm: u32,
    /// Separately authorized anti-collapse floor, absent before calibration.
    pub floor_glups: Option<f64>,
    /// Admission/evidence authority.
    pub evidence: EvidenceClass,
    /// Exact pool/tile placement identity.
    pub placement_identity: String,
    /// Tropical attribution for every retained timed repetition.
    pub critical_paths: Vec<CriticalPathAttribution>,
}

impl D3q19PerfRow {
    /// Validate all bounded axes and citable-floor requirements.
    pub fn validate(&self) -> Result<(), PerfModelError> {
        self.shape.validate()?;
        if !self.glups.is_finite() || self.glups <= 0.0 {
            return Err(PerfModelError::InvalidReceipt("glups"));
        }
        if self.dispersion_ppm > RATIO_PPM {
            return Err(PerfModelError::InvalidReceipt("dispersion_ppm"));
        }
        if self.critical_paths.is_empty() || self.critical_paths.len() > MAX_PERF_REPETITIONS {
            return Err(PerfModelError::ResourceLimit {
                resource: "performance repetitions",
                limit: MAX_PERF_REPETITIONS,
                observed: self.critical_paths.len(),
            });
        }
        if self
            .floor_glups
            .is_some_and(|floor| !floor.is_finite() || floor <= 0.0)
        {
            return Err(PerfModelError::InvalidReceipt("floor_glups"));
        }
        if self.placement_identity.trim().is_empty() || self.placement_identity.len() > 4096 {
            return Err(PerfModelError::InvalidReceipt("placement identity"));
        }
        self.evidence.validate()?;
        for critical_path in &self.critical_paths {
            critical_path.validate()?;
        }
        if matches!(&self.evidence, EvidenceClass::Citable { .. })
            && (self.floor_glups.is_none()
                || self.critical_paths.len() < 3
                || !self.critical_path_stable())
        {
            return Err(PerfModelError::InvalidReceipt(
                "citable floor or stable critical path",
            ));
        }
        Ok(())
    }

    /// Report-only comparison against the raw plan target.
    #[must_use]
    pub fn plan_target_met(&self) -> Option<bool> {
        self.reference_isa
            .plan_target_glups()
            .map(|target| self.glups >= target)
    }

    /// Whether every retained repetition has the same path/class signature.
    #[must_use]
    pub fn critical_path_stable(&self) -> bool {
        critical_path_is_stable(&self.critical_paths)
    }

    /// Gate against only a separately authorized floor and only on citable
    /// evidence. Environment-invalid evidence is neither pass nor fail.
    pub fn gate_verdict(&self) -> Result<PerfGateVerdict, PerfModelError> {
        self.validate()?;
        Ok(match &self.evidence {
            EvidenceClass::Citable { .. } => {
                let Some(floor) = self.floor_glups else {
                    return Err(PerfModelError::InvalidReceipt("citable floor"));
                };
                if self.glups >= floor {
                    PerfGateVerdict::FloorMet
                } else {
                    PerfGateVerdict::FloorMiss
                }
            }
            EvidenceClass::ReportOnly { .. } => PerfGateVerdict::ReportOnly,
            EvidenceClass::EnvironmentInvalid { .. } => PerfGateVerdict::EnvironmentInvalid,
        })
    }

    /// Canonical bounded JSONL row. Human reporting must project this same
    /// typed record rather than recomputing truth from wall samples.
    pub fn receipt_json(&self) -> Result<String, PerfModelError> {
        self.validate()?;
        let target = self
            .reference_isa
            .plan_target_glups()
            .map_or_else(|| "null".to_owned(), |value| format!("{value:.6}"));
        let target_met = self
            .plan_target_met()
            .map_or_else(|| "null".to_owned(), |value| value.to_string());
        let floor = self
            .floor_glups
            .map_or_else(|| "null".to_owned(), |value| format!("{value:.6}"));
        let (admission_receipt, reason) = match &self.evidence {
            EvidenceClass::Citable { admission_receipt } => {
                (format!("\"{admission_receipt}\""), "null".to_owned())
            }
            EvidenceClass::ReportOnly { reason } | EvidenceClass::EnvironmentInvalid { reason } => {
                ("null".to_owned(), format!("\"{}\"", json_escape(reason)))
            }
        };
        let representative = self
            .critical_paths
            .first()
            .ok_or(PerfModelError::InvalidReceipt(
                "representative critical path",
            ))?;
        Ok(format!(
            "{{\"metric\":\"lbm-d3q19-sweep\",\"model_version\":\"{D3Q19_PERF_MODEL_VERSION}\",\"reference_isa\":\"{}\",\"dims\":[{},{},{}],\"occupancy\":\"{}\",\"active_tiles\":{},\"active_cells\":{},\"threading\":\"{}\",\"workers\":{},\"glups\":{:.9},\"dispersion_ppm\":{},\"repetitions\":{},\"plan_target_glups\":{target},\"plan_target_met\":{target_met},\"floor_glups\":{floor},\"evidence_class\":\"{}\",\"gate_verdict\":\"{}\",\"admission_receipt\":{admission_receipt},\"reason\":{reason},\"placement_identity\":\"{}\",\"critical_path_stable\":{},\"critical_path\":{}}}",
            self.reference_isa.as_str(),
            self.shape.dims[0],
            self.shape.dims[1],
            self.shape.dims[2],
            self.shape.occupancy.as_str(),
            self.shape.active_tiles()?,
            self.shape.active_cells()?,
            self.shape.threading.as_str(),
            self.shape.workers,
            self.glups,
            self.dispersion_ppm,
            self.critical_paths.len(),
            self.evidence.as_str(),
            self.gate_verdict()?.as_str(),
            json_escape(&self.placement_identity),
            self.critical_path_stable(),
            representative.receipt_json(),
        ))
    }
}

fn json_escape(value: &str) -> String {
    use core::fmt::Write as _;

    let mut escaped = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            control if control.is_control() => {
                let _ = write!(escaped, "\\u{:04x}", u32::from(control));
            }
            other => escaped.push(other),
        }
    }
    escaped
}
