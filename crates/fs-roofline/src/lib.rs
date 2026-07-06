//! fs-roofline: the roofline harness (plan §14; Decalogue P6).
//!
//! Performance claims as FALSIFIABLE targets: every registered kernel is
//! benchmarked against its arithmetic-intensity-derived limit on the actual
//! machine — measured axes, never spec-sheet numbers — with dispersion
//! reported and results ledgered under the machine fingerprint. "A target
//! that was never re-measured is a lie waiting to happen."
//!
//! Layer: L6 (consumes fs-substrate probes, fs-simd primitives, and writes
//! fs-ledger records). Reporting-only in v0: attainment verdicts inform;
//! gating bands belong to nightly runs on ledgered reference machines.

pub mod axes;
pub mod kernels;
pub mod stats;

pub use axes::MachineAxes;

use fs_ledger::{EventRow, FiveExplicits, Ledger, LedgerError, OpOutcome, now_wall_ns};

/// Crate version (compile-time stamp).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Shape-class key under which roofline rows land in the ledger `tune` table.
pub const TUNE_SHAPE_CLASS: &str = "roofline-v1";

/// Which machine axis a kernel is measured against.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Threading {
    /// One thread: per-core bandwidth/compute axes.
    SingleThread,
    /// All logical cores: aggregate axes.
    AllCore,
}

/// Static description of a benchmarkable kernel: identity plus the
/// arithmetic-intensity model that derives its machine-specific limit.
#[derive(Debug, Clone, Copy)]
pub struct KernelSpec {
    /// Registry name (ledger key; kebab-case).
    pub name: &'static str,
    /// Kernel version (bumped when the implementation changes — attainment
    /// history is only comparable within one version).
    pub version: &'static str,
    /// Bytes moved to/from memory per element processed.
    pub bytes_per_elem: f64,
    /// Floating-point operations per element processed.
    pub flops_per_elem: f64,
    /// Measurement threading model.
    pub threading: Threading,
    /// Target as a fraction of the roofline limit (e.g. 0.85 for "≥85% of
    /// STREAM"). `None` = report-only, no band claimed.
    pub target_fraction: Option<f64>,
}

/// One benchmarkable kernel: owns its buffers; `run_once` is the timed unit.
pub trait RooflineKernel {
    /// The kernel's spec (identity + intensity model + target).
    fn spec(&self) -> KernelSpec;
    /// Elements processed per `run_once` call.
    fn elements(&self) -> usize;
    /// Execute one timed repetition.
    fn run_once(&mut self);
}

/// Which side of the roofline binds the limit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RoofSide {
    /// Memory-bandwidth-bound at this intensity.
    Bandwidth,
    /// Compute-bound at this intensity.
    Compute,
}

/// Verdict against the kernel's declared band.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Verdict {
    /// Attainment ≥ target fraction.
    WithinBand,
    /// Attainment < target fraction.
    BelowBand,
    /// No target declared: report-only.
    NoTarget,
}

impl Verdict {
    /// Stable lowercase name for logs and ledger rows.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Verdict::WithinBand => "within_band",
            Verdict::BelowBand => "below_band",
            Verdict::NoTarget => "no_target",
        }
    }
}

/// One kernel's measured attainment against the machine roofline.
#[derive(Debug, Clone)]
pub struct Attainment {
    /// Kernel name.
    pub kernel: String,
    /// Kernel version.
    pub version: String,
    /// Median elements/second across repetitions.
    pub elems_per_sec: f64,
    /// Achieved memory traffic, GB/s.
    pub achieved_gbs: f64,
    /// Achieved compute, GFLOP/s.
    pub achieved_gflops: f64,
    /// Roofline limit in elements/second for this machine + intensity.
    pub limit_elems_per_sec: f64,
    /// Which axis binds.
    pub roof: RoofSide,
    /// `elems_per_sec / limit_elems_per_sec` (1.0 = at the roof).
    pub attainment: f64,
    /// Relative interquartile dispersion of the repetition times
    /// ((p75 − p25) / median): a benchmark without variance bars is
    /// folklore.
    pub dispersion: f64,
    /// Repetitions measured (after warmup).
    pub reps: usize,
    /// Verdict against the declared band.
    pub verdict: Verdict,
}

impl Attainment {
    /// One JSON line for logs/agents (stable field order).
    #[must_use]
    pub fn to_jsonl(&self) -> String {
        format!(
            "{{\"kernel\":\"{}\",\"version\":\"{}\",\"elems_per_sec\":{:.3e},\
             \"gbs\":{:.3},\"gflops\":{:.3},\"limit_elems_per_sec\":{:.3e},\
             \"roof\":\"{}\",\"attainment\":{:.4},\"dispersion\":{:.4},\
             \"reps\":{},\"verdict\":\"{}\"}}",
            self.kernel,
            self.version,
            self.elems_per_sec,
            self.achieved_gbs,
            self.achieved_gflops,
            self.limit_elems_per_sec,
            match self.roof {
                RoofSide::Bandwidth => "bandwidth",
                RoofSide::Compute => "compute",
            },
            self.attainment,
            self.dispersion,
            self.reps,
            self.verdict.name(),
        )
    }
}

/// Compute attainment from a measured rate and the machine axes. Pure
/// arithmetic — meta-tested against hand calculations (`rf_002`).
#[must_use]
pub fn attainment_for(spec: &KernelSpec, elems_per_sec: f64, axes: &MachineAxes) -> Attainment {
    attainment_with_dispersion(spec, elems_per_sec, 0.0, 0, axes)
}

/// [`attainment_for`] with measured dispersion and repetition count.
#[must_use]
pub fn attainment_with_dispersion(
    spec: &KernelSpec,
    elems_per_sec: f64,
    dispersion: f64,
    reps: usize,
    axes: &MachineAxes,
) -> Attainment {
    let (bandwidth_gbs, peak_gflops) = match spec.threading {
        Threading::SingleThread => (axes.bandwidth_single_gbs, axes.peak_single_gflops),
        Threading::AllCore => (axes.bandwidth_all_core_gbs, axes.peak_all_core_gflops),
    };
    // Limits in elements/second on each axis; +inf when the kernel does not
    // exercise an axis (zero bytes or zero flops per element).
    let bw_limit = if spec.bytes_per_elem > 0.0 {
        bandwidth_gbs * 1e9 / spec.bytes_per_elem
    } else {
        f64::INFINITY
    };
    let comp_limit = if spec.flops_per_elem > 0.0 {
        peak_gflops * 1e9 / spec.flops_per_elem
    } else {
        f64::INFINITY
    };
    let (limit, roof) = if bw_limit <= comp_limit {
        (bw_limit, RoofSide::Bandwidth)
    } else {
        (comp_limit, RoofSide::Compute)
    };
    let attainment = if limit.is_finite() && limit > 0.0 {
        elems_per_sec / limit
    } else {
        0.0
    };
    let verdict = match spec.target_fraction {
        None => Verdict::NoTarget,
        Some(t) if attainment >= t => Verdict::WithinBand,
        Some(_) => Verdict::BelowBand,
    };
    Attainment {
        kernel: spec.name.to_string(),
        version: spec.version.to_string(),
        elems_per_sec,
        achieved_gbs: elems_per_sec * spec.bytes_per_elem / 1e9,
        achieved_gflops: elems_per_sec * spec.flops_per_elem / 1e9,
        limit_elems_per_sec: limit,
        roof,
        attainment,
        dispersion,
        reps,
        verdict,
    }
}

/// Measure one kernel (warmup + repetitions) and compute its attainment.
pub fn measure(
    kernel: &mut dyn RooflineKernel,
    warmup: usize,
    reps: usize,
    axes: &MachineAxes,
) -> Attainment {
    let spec = kernel.spec();
    let elems = kernel.elements() as f64;
    let sample = stats::time_reps(&mut || kernel.run_once(), warmup, reps);
    let elems_per_sec = if sample.median > 0.0 {
        elems / sample.median
    } else {
        0.0
    };
    attainment_with_dispersion(&spec, elems_per_sec, sample.dispersion, reps, axes)
}

/// Run every kernel in the registry.
pub fn run_registry(
    registry: &mut [Box<dyn RooflineKernel>],
    warmup: usize,
    reps: usize,
    axes: &MachineAxes,
) -> Vec<Attainment> {
    registry
        .iter_mut()
        .map(|k| measure(k.as_mut(), warmup, reps, axes))
        .collect()
}

// ---------------------------------------------------------------------------
// §14.1 target table as data
// ---------------------------------------------------------------------------

/// One row of the plan §14.1 target table. `landed = false` rows are
/// visible from day one so nothing is silently uncovered — they flip as the
/// owning kernels register.
#[derive(Debug, Clone, Copy)]
pub struct TargetRow {
    /// Kernel family name.
    pub kernel: &'static str,
    /// What the target means (unit and roof context).
    pub statement: &'static str,
    /// Whether an implementation is registered in this harness yet.
    pub landed: bool,
}

/// The §14.1 table. Statements, not claims: every value must be re-measured
/// on a fingerprinted machine before anyone may cite it.
pub const SECTION_14_1_TARGETS: &[TargetRow] = &[
    TargetRow {
        kernel: "lbm-d3q19-stream-collide",
        statement: "≥1.0 GLUP/s (M-class) / ≥0.6 GLUP/s (TR-class), bandwidth-bound",
        landed: false,
    },
    TargetRow {
        kernel: "gemm-f64",
        statement: "≥75% of measured peak FLOPs for the selected SIMD tier",
        landed: false,
    },
    TargetRow {
        kernel: "spmv-sell-c-sigma",
        statement: "≥85% of measured STREAM-class bandwidth",
        landed: false,
    },
    TargetRow {
        kernel: "feec-apply-p4",
        statement: "≥30% of peak FLOPs, sum-factorized",
        landed: false,
    },
    TargetRow {
        kernel: "batched-small-dense",
        statement: "≥60% of peak FLOPs, SIMD-across-elements",
        landed: false,
    },
    TargetRow {
        kernel: "fft-3d-pencil",
        statement: "≥40% of the memory-bound limit",
        landed: false,
    },
    TargetRow {
        kernel: "sdf-primary-rays",
        statement: "≥80 Mray/s (M-class) / ≥120 Mray/s (TR-class)",
        landed: false,
    },
];

// ---------------------------------------------------------------------------
// Ledger integration and staleness
// ---------------------------------------------------------------------------

/// Record a harness run in the ledger: one op (frozen Five Explicits),
/// per-kernel metric rows, a `benchmark_result` event per kernel, and tune
/// rows keyed by machine fingerprint.
///
/// # Errors
/// Ledger errors propagate; the op is finished `Error` when any kernel row
/// fails to record.
pub fn record_run(
    ledger: &Ledger,
    axes: &MachineAxes,
    results: &[Attainment],
) -> Result<i64, LedgerError> {
    let versions = format!(
        "{{\"frankensim\":\"{}\",\"fs-roofline\":\"{VERSION}\"}}",
        std::env::var("GITHUB_SHA").unwrap_or_else(|_| "local".to_string())
    );
    let explicits = FiveExplicits {
        seed: b"roofline",
        versions: &versions,
        budget: "{\"wall_s\":600}",
        capability: "{\"ops\":[\"perf.roofline\"]}",
    };
    let ir = format!(
        "{{\"op\":\"perf.roofline\",\"kernels\":{},\"fingerprint\":\"{:016x}\"}}",
        results.len(),
        axes.fingerprint
    );
    let op = ledger.begin_op(Some(b"roofline"), &ir, &explicits, now_wall_ns())?;
    let fp_bytes = axes.fingerprint.to_le_bytes();
    for r in results {
        ledger.record_metric(
            op,
            0,
            &format!("{}.elems_per_sec", r.kernel),
            r.elems_per_sec,
        )?;
        ledger.record_metric(op, 0, &format!("{}.attainment", r.kernel), r.attainment)?;
        ledger.record_metric(op, 0, &format!("{}.dispersion", r.kernel), r.dispersion)?;
        ledger.append_event(&EventRow {
            session: Some(b"roofline"),
            t: 0,
            kind: "benchmark_result",
            payload: Some(&r.to_jsonl()),
        })?;
        ledger.tune_put(
            &r.kernel,
            TUNE_SHAPE_CLASS,
            &fp_bytes,
            &format!("{{\"version\":\"{}\",\"reps\":{}}}", r.version, r.reps),
            &r.to_jsonl(),
        )?;
    }
    ledger.finish_op(op, OpOutcome::Ok, None, now_wall_ns())?;
    Ok(op)
}

/// Staleness state of one kernel's ledgered attainment on this machine.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Staleness {
    /// A row exists for the current machine fingerprint.
    Fresh,
    /// Rows exist, but none for the current fingerprint — the machine
    /// drifted and every cited number is stale until re-measured.
    FingerprintDrift,
    /// No roofline rows at all: never measured.
    NeverMeasured,
}

/// Check one kernel's staleness against the current fingerprint.
///
/// # Errors
/// Ledger errors propagate.
pub fn staleness(
    ledger: &Ledger,
    kernel: &str,
    current_fingerprint: u64,
) -> Result<Staleness, LedgerError> {
    let rows = ledger.tune_rows(kernel)?;
    let roofline_rows: Vec<_> = rows
        .iter()
        .filter(|r| r.shape_class == TUNE_SHAPE_CLASS)
        .collect();
    if roofline_rows.is_empty() {
        return Ok(Staleness::NeverMeasured);
    }
    let fp = current_fingerprint.to_le_bytes();
    if roofline_rows.iter().any(|r| r.machine == fp) {
        Ok(Staleness::Fresh)
    } else {
        Ok(Staleness::FingerprintDrift)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn synthetic_axes() -> MachineAxes {
        MachineAxes {
            fingerprint: 0xABCD,
            cpu_brand: "synthetic".to_string(),
            logical_cpus: 8,
            bandwidth_single_gbs: 100.0,
            bandwidth_all_core_gbs: 400.0,
            peak_single_gflops: 50.0,
            peak_all_core_gflops: 300.0,
        }
    }

    #[test]
    fn attainment_matches_hand_calculation_bandwidth_bound() {
        // axpy: 24 B/elem, 2 flop/elem on axes (100 GB/s, 50 GFLOP/s):
        // bw limit = 100e9/24 = 4.1667e9 elem/s; compute = 50e9/2 = 25e9.
        // Bandwidth binds. At 2.0833e9 elem/s attainment = 0.5 exactly.
        let spec = KernelSpec {
            name: "axpy",
            version: "1",
            bytes_per_elem: 24.0,
            flops_per_elem: 2.0,
            threading: Threading::SingleThread,
            target_fraction: Some(0.6),
        };
        let a = attainment_for(&spec, 100.0e9 / 24.0 / 2.0, &synthetic_axes());
        assert_eq!(a.roof, RoofSide::Bandwidth);
        assert!((a.attainment - 0.5).abs() < 1e-12, "got {}", a.attainment);
        assert!((a.achieved_gbs - 50.0).abs() < 1e-9);
        assert_eq!(a.verdict, Verdict::BelowBand);
    }

    #[test]
    fn attainment_matches_hand_calculation_compute_bound() {
        // High-intensity kernel: 1 B/elem, 100 flop/elem.
        // bw limit = 100e9 elem/s; compute = 50e9/100 = 0.5e9 → compute binds.
        let spec = KernelSpec {
            name: "dense",
            version: "1",
            bytes_per_elem: 1.0,
            flops_per_elem: 100.0,
            threading: Threading::SingleThread,
            target_fraction: Some(0.5),
        };
        let a = attainment_for(&spec, 0.4e9, &synthetic_axes());
        assert_eq!(a.roof, RoofSide::Compute);
        assert!((a.attainment - 0.8).abs() < 1e-12);
        assert_eq!(a.verdict, Verdict::WithinBand);
        // All-core axes flip the limit.
        let all = KernelSpec {
            threading: Threading::AllCore,
            ..spec
        };
        let b = attainment_for(&all, 0.4e9, &synthetic_axes());
        assert!((b.limit_elems_per_sec - 3.0e9).abs() < 1.0);
    }

    #[test]
    fn no_target_reports_without_verdict() {
        let spec = KernelSpec {
            name: "probe",
            version: "1",
            bytes_per_elem: 8.0,
            flops_per_elem: 1.0,
            threading: Threading::SingleThread,
            target_fraction: None,
        };
        let a = attainment_for(&spec, 1.0e9, &synthetic_axes());
        assert_eq!(a.verdict, Verdict::NoTarget);
        assert!(a.to_jsonl().contains("\"verdict\":\"no_target\""));
    }

    #[test]
    fn section_14_1_table_is_complete_and_honest() {
        assert_eq!(SECTION_14_1_TARGETS.len(), 7, "all §14.1 families present");
        // Nothing may claim to be landed until its kernel registers here.
        for row in SECTION_14_1_TARGETS {
            assert!(
                !row.landed,
                "{} claims landed without a registered kernel",
                row.kernel
            );
        }
    }
}
