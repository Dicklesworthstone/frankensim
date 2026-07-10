//! Machine roofline axes: measured bandwidth and measured peak FLOPs.
//!
//! NEVER spec-sheet numbers (plan §14.1). Bandwidth comes from
//! fs-substrate's STREAM-triad probe; peak FLOPs from an in-house fused
//! multiply-add chain microbench (independent accumulator lanes so the
//! autovectorizer can fill the SIMD units). The compute axis is therefore
//! "compiler-achievable peak" — conservative on tiers the autovectorizer
//! misses, which is the honest direction for a limit that divides other
//! kernels' attainment.

use fs_substrate::CapabilityProbe;

/// Lowest credible single-thread STREAM bandwidth on FrankenSim's
/// post-2015 reference-machine families. Values below this are evidence
/// that the probe ran in a crushed or otherwise unusable environment.
pub const MIN_REFERENCE_SINGLE_BANDWIDTH_GBS: f64 = 5.0;

/// Lowest credible single-thread fused throughput on FrankenSim's
/// post-2015 reference-machine families.
pub const MIN_REFERENCE_SINGLE_PEAK_GFLOPS: f64 = 5.0;

/// Maximum relative change allowed between the pre-run and post-run probes.
pub const MAX_AXIS_REPROBE_DRIFT: f64 = 0.25;

/// The roofline's axes for one machine, all measured.
#[derive(Debug, Clone)]
pub struct MachineAxes {
    /// fs-substrate topology fingerprint (bandwidth excluded by design).
    pub fingerprint: u64,
    /// CPU brand string (report context).
    pub cpu_brand: String,
    /// Logical CPU count.
    pub logical_cpus: u32,
    /// Measured single-thread STREAM-triad bandwidth, GB/s.
    pub bandwidth_single_gbs: f64,
    /// Measured all-core STREAM-triad bandwidth, GB/s.
    pub bandwidth_all_core_gbs: f64,
    /// Measured single-thread FMA throughput, GFLOP/s.
    pub peak_single_gflops: f64,
    /// Measured all-core FMA throughput, GFLOP/s.
    pub peak_all_core_gflops: f64,
}

impl MachineAxes {
    /// Probe the current machine (topology + bandwidth via fs-substrate,
    /// peak FLOPs via the FMA microbench). Takes a few hundred ms by design.
    #[must_use]
    pub fn probe() -> MachineAxes {
        let probe = CapabilityProbe::run();
        let peak_single = fma_gflops_single();
        let peak_all = fma_gflops_all_core(probe.logical_cpus as usize);
        MachineAxes {
            fingerprint: probe.fingerprint(),
            cpu_brand: probe.cpu_brand.clone(),
            logical_cpus: probe.logical_cpus,
            bandwidth_single_gbs: probe.measured.single_thread_gbs,
            bandwidth_all_core_gbs: probe.measured.all_core_gbs,
            peak_single_gflops: peak_single,
            peak_all_core_gflops: peak_all,
        }
    }

    /// Explain why these axes cannot support an attainment verdict.
    ///
    /// Absolute floors catch the bead-1n61 counterexample where both the
    /// axis and kernel collapsed together. Finite/positive checks keep
    /// malformed or overflowed measurements out of JSON and the ledger.
    /// All-core measurements may be noisy, but falling below half the
    /// single-thread value means the aggregate probe was not usable.
    #[must_use]
    pub fn plausibility_error(&self) -> Option<&'static str> {
        if self.logical_cpus == 0 {
            return Some("logical CPU count is zero");
        }
        let axes = [
            self.bandwidth_single_gbs,
            self.bandwidth_all_core_gbs,
            self.peak_single_gflops,
            self.peak_all_core_gflops,
        ];
        if axes.iter().any(|value| !value.is_finite() || *value <= 0.0) {
            return Some("one or more measured axes are non-finite or non-positive");
        }
        if self.bandwidth_single_gbs < MIN_REFERENCE_SINGLE_BANDWIDTH_GBS {
            return Some("single-thread bandwidth is below the reference-family floor");
        }
        if self.peak_single_gflops < MIN_REFERENCE_SINGLE_PEAK_GFLOPS {
            return Some("single-thread FMA throughput is below the reference-family floor");
        }
        if self.bandwidth_all_core_gbs < self.bandwidth_single_gbs * 0.5 {
            return Some("all-core bandwidth is less than half the single-thread axis");
        }
        if self.peak_all_core_gflops < self.peak_single_gflops * 0.5 {
            return Some("all-core FMA throughput is less than half the single-thread axis");
        }
        None
    }

    /// Whether these axes can support an attainment verdict.
    #[must_use]
    pub fn plausible(&self) -> bool {
        self.plausibility_error().is_none()
    }

    /// Explain why a post-run probe does not corroborate this pre-run probe.
    #[must_use]
    pub fn reprobe_error(&self, after: &Self) -> Option<&'static str> {
        if self.plausibility_error().is_some() || after.plausibility_error().is_some() {
            return Some("pre-run or post-run axes are implausible");
        }
        if self.fingerprint != after.fingerprint || self.logical_cpus != after.logical_cpus {
            return Some("machine identity changed between axis probes");
        }
        let pairs = [
            (self.bandwidth_single_gbs, after.bandwidth_single_gbs),
            (self.bandwidth_all_core_gbs, after.bandwidth_all_core_gbs),
            (self.peak_single_gflops, after.peak_single_gflops),
            (self.peak_all_core_gflops, after.peak_all_core_gflops),
        ];
        if pairs.iter().any(|&(before, after)| {
            (before - after).abs() / before.abs().max(after.abs()) > MAX_AXIS_REPROBE_DRIFT
        }) {
            return Some("pre-run and post-run axes disagree beyond the drift band");
        }
        None
    }

    /// One JSON line describing the axes (report/ledger context).
    #[must_use]
    pub fn to_jsonl(&self) -> String {
        format!(
            "{{\"fingerprint\":\"{:016x}\",\"cpu\":\"{}\",\"logical_cpus\":{},\
             \"bandwidth_single_gbs\":{},\"bandwidth_all_core_gbs\":{},\
             \"peak_single_gflops\":{},\"peak_all_core_gflops\":{}}}",
            self.fingerprint,
            json_escape(&self.cpu_brand),
            self.logical_cpus,
            json_axis(self.bandwidth_single_gbs),
            json_axis(self.bandwidth_all_core_gbs),
            json_axis(self.peak_single_gflops),
            json_axis(self.peak_all_core_gflops),
        )
    }
}

fn json_axis(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.2}")
    } else {
        "null".to_string()
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

/// Independent FMA accumulator lanes per chain step — sized to the
/// architecture's REGISTER FILE, not just its issue width: 64 f64 fill
/// 8 zmm (AVX-512) or 16 ymm (AVX2) exactly, but need 32 NEON vregs —
/// all of them — so on aarch64 the accumulators spilled every step and
/// the "peak" axis read ~25% low (and unstably). 48 lanes (24 vregs)
/// sit on the measured M4 plateau (32/48 lanes within 0.2%; 64 lanes
/// 26% lower — ledgered in bead xdgf).
#[cfg(target_arch = "aarch64")]
const LANES: usize = 48;
#[cfg(not(target_arch = "aarch64"))]
const LANES: usize = 64;
/// Chain steps per timed pass.
const STEPS: usize = 4096;
/// Timed passes per measurement; best-of keeps thermal/scheduler noise from
/// deflating the axis (an axis too low would inflate everyone's attainment).
const PASSES: usize = 5;

use fma_kernel::fma_pass;

/// The FMA chain kernel (own registered capsule — the x86 body needs
/// `target_feature`; everything else in this file is safe).
#[path = "fma_kernel.rs"]
mod fma_kernel;

/// Best-of-passes single-thread FMA throughput in GFLOP/s.
///
/// Each timed sample REPEATS the chain until it spans ≥ 5 ms of wall
/// clock: a single 4096-step pass is only a few microseconds, which is
/// inside the frequency-ramp/scheduler-noise floor on modern cores and
/// made the axis wander tens of percent between probes — an unstable
/// denominator that can push honest kernels past attainment 1.0.
#[must_use]
pub fn fma_gflops_single() -> f64 {
    let mut acc = [1.0f64; LANES];
    // Multiplier chosen so accumulators orbit without overflow/denormal.
    let (m, a) = (0.999_999_9, 1.0e-9);
    let flops_per_pass = (2 * LANES * STEPS) as f64;
    // Calibrate the repeat count to reach the 5 ms sample floor.
    let mut reps = 1usize;
    loop {
        let start = std::time::Instant::now();
        for _ in 0..reps {
            fma_pass(&mut acc, m, a);
        }
        if start.elapsed().as_secs_f64() >= 5e-3 || reps >= 1 << 20 {
            break;
        }
        reps *= 2;
    }
    let mut best = 0.0f64;
    for _ in 0..PASSES {
        let start = std::time::Instant::now();
        for _ in 0..reps {
            fma_pass(&mut acc, m, a);
        }
        let dt = start.elapsed().as_secs_f64();
        if dt > 0.0 {
            best = best.max(flops_per_pass * reps as f64 / dt / 1e9);
        }
    }
    std::hint::black_box(acc[LANES / 2]);
    best
}

/// All-core FMA throughput: every logical CPU runs the single-thread bench
/// concurrently; the sum is the aggregate compute axis.
#[must_use]
pub fn fma_gflops_all_core(logical_cpus: usize) -> f64 {
    let threads = logical_cpus.max(1);
    std::thread::scope(|s| {
        let handles: Vec<_> = (0..threads).map(|_| s.spawn(fma_gflops_single)).collect();
        handles.into_iter().map(|h| h.join().unwrap_or(0.0)).sum()
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fma_bench_reports_positive_throughput() {
        let gflops = fma_gflops_single();
        assert!(
            gflops > 0.05,
            "implausibly low FMA throughput: {gflops} GFLOP/s"
        );
    }

    #[test]
    fn probe_axes_are_positive_and_fingerprint_stable() {
        let a = MachineAxes::probe();
        assert!(a.bandwidth_single_gbs > 0.0);
        assert!(a.bandwidth_all_core_gbs > 0.0);
        assert!(a.peak_single_gflops > 0.0);
        assert!(a.peak_all_core_gflops >= a.peak_single_gflops * 0.5);
        // Fingerprint covers topology only: re-probing the same machine in
        // the same process must agree.
        let b = MachineAxes::probe();
        assert_eq!(a.fingerprint, b.fingerprint);
        assert!(a.to_jsonl().contains("fingerprint"));
    }

    #[test]
    fn plausibility_checks_every_axis_and_topology_field() {
        let healthy = MachineAxes {
            fingerprint: 1,
            cpu_brand: "synthetic".to_string(),
            logical_cpus: 8,
            bandwidth_single_gbs: 50.0,
            bandwidth_all_core_gbs: 200.0,
            peak_single_gflops: 40.0,
            peak_all_core_gflops: 240.0,
        };
        assert!(healthy.plausible());
        for field in 0..4 {
            let mut poisoned = healthy.clone();
            match field {
                0 => poisoned.bandwidth_single_gbs = f64::NAN,
                1 => poisoned.bandwidth_all_core_gbs = f64::INFINITY,
                2 => poisoned.peak_single_gflops = -1.0,
                _ => poisoned.peak_all_core_gflops = 0.0,
            }
            assert!(
                !poisoned.plausible(),
                "axis field {field} escaped validation"
            );
            let json = poisoned.to_jsonl();
            assert!(!json.contains("NaN") && !json.contains("inf"));
        }
        let mut no_cpus = healthy;
        no_cpus.logical_cpus = 0;
        assert!(!no_cpus.plausible());
    }

    #[test]
    fn post_run_probe_must_match_machine_and_axis_band() {
        let before = MachineAxes {
            fingerprint: 1,
            cpu_brand: "synthetic".to_string(),
            logical_cpus: 8,
            bandwidth_single_gbs: 50.0,
            bandwidth_all_core_gbs: 200.0,
            peak_single_gflops: 40.0,
            peak_all_core_gflops: 240.0,
        };
        assert!(before.reprobe_error(&before).is_none());
        let mut drifted = before.clone();
        drifted.bandwidth_single_gbs = 20.0;
        assert!(before.reprobe_error(&drifted).is_some());
        let mut other = before.clone();
        other.fingerprint = 2;
        assert!(before.reprobe_error(&other).is_some());
    }
}
