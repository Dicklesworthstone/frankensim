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

pub mod regress;

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
    /// The MEASUREMENT ENVIRONMENT is invalid: the probed axes fail
    /// absolute plausibility floors, or the kernel "beat" its roofline
    /// by an impossible margin (stale/contention-crushed axes). A gate
    /// must never pass — or fail — on a machine that was useless while
    /// it was measured (bead 1n61: a load-68 window collapsed both the
    /// STREAM probe and the kernel ~1000× together, and the RATIO
    /// self-normalized to a vacuous within_band).
    EnvironmentInvalid,
}

impl Verdict {
    /// Stable lowercase name for logs and ledger rows.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Verdict::WithinBand => "within_band",
            Verdict::BelowBand => "below_band",
            Verdict::NoTarget => "no_target",
            Verdict::EnvironmentInvalid => "environment_invalid",
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
    /// Why this row cannot support a verdict. Present exactly when
    /// `verdict == EnvironmentInvalid`.
    pub invalid_reason: Option<String>,
}

impl Attainment {
    /// One JSON line for logs/agents (stable field order).
    #[must_use]
    pub fn to_jsonl(&self) -> String {
        let invalid_reason = self.invalid_reason.as_ref().map_or_else(
            || "null".to_string(),
            |reason| format!("\"{}\"", json_escape(reason)),
        );
        format!(
            "{{\"kernel\":\"{}\",\"version\":\"{}\",\"elems_per_sec\":{:.3e},\
             \"gbs\":{:.3},\"gflops\":{:.3},\"limit_elems_per_sec\":{:.3e},\
             \"roof\":\"{}\",\"attainment\":{:.4},\"dispersion\":{:.4},\
             \"reps\":{},\"verdict\":\"{}\",\"invalid_reason\":{}}}",
            json_escape(&self.kernel),
            json_escape(&self.version),
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
            invalid_reason,
        )
    }
}

fn json_escape(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            c if c.is_control() => {
                use std::fmt::Write as _;
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out
}

fn spec_error(spec: &KernelSpec) -> Option<&'static str> {
    if spec.name.trim().is_empty() || spec.version.trim().is_empty() {
        return Some("kernel name and version must be non-empty");
    }
    if !spec.bytes_per_elem.is_finite()
        || !spec.flops_per_elem.is_finite()
        || spec.bytes_per_elem < 0.0
        || spec.flops_per_elem < 0.0
        || (spec.bytes_per_elem == 0.0 && spec.flops_per_elem == 0.0)
    {
        return Some("kernel intensity must be finite, non-negative, and exercise an axis");
    }
    if let Some(target) = spec.target_fraction
        && (!target.is_finite() || target <= 0.0 || target > 1.0)
    {
        return Some("target fraction must be finite and in (0, 1]");
    }
    None
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
    let spec_error = spec_error(spec);
    let measurement_error = if !elems_per_sec.is_finite() || elems_per_sec < 0.0 {
        Some("measured element rate is non-finite or negative")
    } else if !dispersion.is_finite() || dispersion < 0.0 {
        Some("measured dispersion is non-finite or negative")
    } else {
        None
    };
    let safe_rate = if measurement_error.is_none() {
        elems_per_sec
    } else {
        0.0
    };
    let safe_bytes = if spec.bytes_per_elem.is_finite() && spec.bytes_per_elem >= 0.0 {
        spec.bytes_per_elem
    } else {
        0.0
    };
    let safe_flops = if spec.flops_per_elem.is_finite() && spec.flops_per_elem >= 0.0 {
        spec.flops_per_elem
    } else {
        0.0
    };
    let (bandwidth_gbs, peak_gflops) = match spec.threading {
        Threading::SingleThread => (axes.bandwidth_single_gbs, axes.peak_single_gflops),
        Threading::AllCore => (axes.bandwidth_all_core_gbs, axes.peak_all_core_gflops),
    };
    // Limits in elements/second on each axis; +inf when the kernel does not
    // exercise an axis (zero bytes or zero flops per element).
    let bw_limit = if safe_bytes > 0.0 && bandwidth_gbs.is_finite() && bandwidth_gbs > 0.0 {
        bandwidth_gbs * 1e9 / safe_bytes
    } else {
        f64::INFINITY
    };
    let comp_limit = if safe_flops > 0.0 && peak_gflops.is_finite() && peak_gflops > 0.0 {
        peak_gflops * 1e9 / safe_flops
    } else {
        f64::INFINITY
    };
    let (limit, roof) = if bw_limit <= comp_limit {
        (bw_limit, RoofSide::Bandwidth)
    } else {
        (comp_limit, RoofSide::Compute)
    };
    let raw_attainment = if limit.is_finite() && limit > 0.0 {
        safe_rate / limit
    } else {
        0.0
    };
    let raw_achieved_gbs = safe_rate * safe_bytes / 1e9;
    let raw_achieved_gflops = safe_rate * safe_flops / 1e9;
    let derived_error = if !raw_attainment.is_finite()
        || !raw_achieved_gbs.is_finite()
        || !raw_achieved_gflops.is_finite()
    {
        Some("derived roofline quantities overflowed or became non-finite")
    } else {
        None
    };
    let attainment = if raw_attainment.is_finite() {
        raw_attainment
    } else {
        0.0
    };
    // Environment validity BEFORE band comparison (bead 1n61): the
    // ratio is meaningless when the axes are implausible, and an
    // attainment materially above 1 means the kernel outran its own
    // roofline — the axes were probed under different (crushed)
    // conditions than the kernel run. Refuse to gate either way.
    let invalid_reason = axes
        .plausibility_error()
        .or(spec_error)
        .or(measurement_error)
        .or(derived_error)
        .or_else(|| (attainment > 1.5).then_some("attainment exceeds the credible roofline band"));
    let verdict = if invalid_reason.is_some() {
        Verdict::EnvironmentInvalid
    } else {
        match spec.target_fraction {
            None => Verdict::NoTarget,
            Some(t) if attainment >= t => Verdict::WithinBand,
            Some(_) => Verdict::BelowBand,
        }
    };
    Attainment {
        kernel: spec.name.to_string(),
        version: spec.version.to_string(),
        elems_per_sec: safe_rate,
        achieved_gbs: if raw_achieved_gbs.is_finite() {
            raw_achieved_gbs
        } else {
            0.0
        },
        achieved_gflops: if raw_achieved_gflops.is_finite() {
            raw_achieved_gflops
        } else {
            0.0
        },
        limit_elems_per_sec: if limit.is_finite() { limit } else { 0.0 },
        roof,
        attainment,
        dispersion: if dispersion.is_finite() && dispersion >= 0.0 {
            dispersion
        } else {
            0.0
        },
        reps,
        verdict,
        invalid_reason: invalid_reason.map(str::to_string),
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
    attainment_with_dispersion(&spec, elems_per_sec, sample.dispersion, reps.max(1), axes)
}

/// Run every kernel in the registry.
pub fn run_registry(
    registry: &mut [Box<dyn RooflineKernel>],
    warmup: usize,
    reps: usize,
    axes: &MachineAxes,
) -> Vec<Attainment> {
    let mut results: Vec<_> = registry
        .iter_mut()
        .map(|k| measure(k.as_mut(), warmup, reps, axes))
        .collect();
    poison_invalid_run(&mut results);
    results
}

fn poison_invalid_run(results: &mut [Attainment]) {
    let Some(origin) = results
        .iter()
        .find(|result| result.verdict == Verdict::EnvironmentInvalid)
    else {
        return;
    };
    let reason = format!(
        "registry invalidated by {}: {}",
        origin.kernel,
        origin
            .invalid_reason
            .as_deref()
            .unwrap_or("invalid evidence")
    );
    for result in results {
        result.verdict = Verdict::EnvironmentInvalid;
        result.invalid_reason = Some(reason.clone());
    }
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
    let run_valid = results
        .iter()
        .all(|result| result.verdict != Verdict::EnvironmentInvalid);
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
        if run_valid {
            ledger.tune_put(
                &r.kernel,
                TUNE_SHAPE_CLASS,
                &fp_bytes,
                &format!(
                    "{{\"version\":\"{}\",\"reps\":{}}}",
                    json_escape(&r.version),
                    r.reps
                ),
                &r.to_jsonl(),
            )?;
        }
    }
    if run_valid {
        ledger.finish_op(op, OpOutcome::Ok, None, now_wall_ns())?;
    } else {
        ledger.append_event(&EventRow {
            session: Some(b"roofline"),
            t: 0,
            kind: "roofline_run_invalid",
            payload: Some(
                "{\"code\":\"environment_invalid\",\"effect\":\"no_tune_rows_published\"}",
            ),
        })?;
        ledger.finish_op(
            op,
            OpOutcome::Error,
            Some(
                "{\"code\":\"roofline_environment_invalid\",\"action\":\"rerun_on_a_quiet_reference_machine\"}",
            ),
            now_wall_ns(),
        )?;
    }
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
    fn crushed_axes_cannot_gate_vacuously() {
        // Bead 1n61 counterexample replay: on a load-68 host both the
        // STREAM probe and the FFT kernel collapsed ~1000× together
        // (axes 0.2 GB/s, kernel 0.156 GB/s on a ~200 GB/s machine) and
        // the RATIO self-normalized to 0.89 = a vacuous within_band.
        let crushed = MachineAxes {
            fingerprint: 0xDEAD,
            cpu_brand: "crushed".to_string(),
            logical_cpus: 128,
            bandwidth_single_gbs: 0.2,
            bandwidth_all_core_gbs: 0.4,
            peak_single_gflops: 0.05,
            peak_all_core_gflops: 0.4,
        };
        assert!(!crushed.plausible());
        let spec = KernelSpec {
            name: "fft-roundtrip",
            version: "1n61",
            bytes_per_elem: 672.0,
            flops_per_elem: 172.0,
            threading: Threading::SingleThread,
            target_fraction: Some(0.40),
        };
        // 0.156 GB/s effective — 89% of the crushed axis.
        let a = attainment_for(&spec, 0.156e9 / 672.0, &crushed);
        assert_eq!(
            a.verdict,
            Verdict::EnvironmentInvalid,
            "a crushed environment must refuse to gate (got attainment {:.3})",
            a.attainment
        );
    }

    #[test]
    fn over_roof_attainment_poisons_the_gate() {
        // Healthy axes but the kernel 'beats' its roofline by 2× — the
        // axes are stale relative to the run; refuse.
        let spec = KernelSpec {
            name: "axpy",
            version: "1n61",
            bytes_per_elem: 24.0,
            flops_per_elem: 2.0,
            threading: Threading::SingleThread,
            target_fraction: Some(0.6),
        };
        let a = attainment_for(&spec, 2.0 * 100.0e9 / 24.0, &synthetic_axes());
        assert_eq!(a.verdict, Verdict::EnvironmentInvalid);
        // Slightly over 1 (measurement jitter) still gates normally.
        let b = attainment_for(&spec, 1.2 * 100.0e9 / 24.0, &synthetic_axes());
        assert_eq!(b.verdict, Verdict::WithinBand);
    }

    #[test]
    fn invalid_numeric_inputs_fail_closed_and_remain_json() {
        let base = KernelSpec {
            name: "probe\"escaped",
            version: "1",
            bytes_per_elem: 8.0,
            flops_per_elem: 1.0,
            threading: Threading::SingleThread,
            target_fraction: Some(0.5),
        };
        let bad_rate = attainment_with_dispersion(&base, f64::NAN, 0.0, 3, &synthetic_axes());
        assert_eq!(bad_rate.verdict, Verdict::EnvironmentInvalid);
        assert!(bad_rate.elems_per_sec.is_finite());
        assert!(bad_rate.to_jsonl().contains("probe\\\"escaped"));
        assert!(!bad_rate.to_jsonl().contains("NaN"));

        let bad_dispersion =
            attainment_with_dispersion(&base, 1.0, f64::INFINITY, 3, &synthetic_axes());
        assert_eq!(bad_dispersion.verdict, Verdict::EnvironmentInvalid);
        let bad_target = KernelSpec {
            target_fraction: Some(f64::NAN),
            ..base
        };
        assert_eq!(
            attainment_for(&bad_target, 1.0, &synthetic_axes()).verdict,
            Verdict::EnvironmentInvalid
        );
    }

    #[test]
    fn one_invalid_row_poisons_every_registry_verdict() {
        let normal = KernelSpec {
            name: "normal",
            version: "1",
            bytes_per_elem: 24.0,
            flops_per_elem: 2.0,
            threading: Threading::SingleThread,
            target_fraction: Some(0.5),
        };
        let impossible = KernelSpec {
            name: "impossible",
            ..normal
        };
        let mut results = vec![
            attainment_for(&normal, 100.0e9 / 24.0 * 0.8, &synthetic_axes()),
            attainment_for(&impossible, 100.0e9 / 24.0 * 2.0, &synthetic_axes()),
        ];
        assert_eq!(results[0].verdict, Verdict::WithinBand);
        assert_eq!(results[1].verdict, Verdict::EnvironmentInvalid);
        poison_invalid_run(&mut results);
        assert!(
            results
                .iter()
                .all(|row| row.verdict == Verdict::EnvironmentInvalid)
        );
        assert!(results.iter().all(|row| {
            row.invalid_reason
                .as_deref()
                .is_some_and(|r| r.contains("impossible"))
        }));
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
