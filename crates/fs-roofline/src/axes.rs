//! Machine roofline axes: measured bandwidth and measured peak FLOPs.
//!
//! NEVER spec-sheet numbers (plan §14.1). Bandwidth comes from
//! fs-substrate's STREAM-triad probe; peak FLOPs from an in-house fused
//! multiply-add chain microbench (independent accumulator lanes so the
//! autovectorizer can fill the SIMD units). The compute axis is therefore
//! "compiler-achievable peak" — conservative on tiers the autovectorizer
//! misses, which is the honest direction for a limit that divides other
//! kernels' attainment.

#![allow(unsafe_code)]

use fs_substrate::CapabilityProbe;

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

    /// One JSON line describing the axes (report/ledger context).
    #[must_use]
    pub fn to_jsonl(&self) -> String {
        format!(
            "{{\"fingerprint\":\"{:016x}\",\"cpu\":\"{}\",\"logical_cpus\":{},\
             \"bandwidth_single_gbs\":{:.2},\"bandwidth_all_core_gbs\":{:.2},\
             \"peak_single_gflops\":{:.2},\"peak_all_core_gflops\":{:.2}}}",
            self.fingerprint,
            self.cpu_brand,
            self.logical_cpus,
            self.bandwidth_single_gbs,
            self.bandwidth_all_core_gbs,
            self.peak_single_gflops,
            self.peak_all_core_gflops,
        )
    }
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

fn fma_pass_portable(acc: &mut [f64; LANES], m: f64, a: f64) {
    for _ in 0..STEPS {
        for lane in acc.iter_mut() {
            *lane = lane.mul_add(m, a);
        }
    }
}

/// x86 body with REAL fused-multiply-add codegen: the baseline x86-64
/// target has no compile-time FMA, so `f64::mul_add` lowers to a libm
/// CALL and the "peak" axis measured call overhead — 1.0 GFLOP/s on a
/// Threadripper whose GEMM then read attainment 40× (found by the
/// xlvx x86 attainment row). Runtime-detected, capsule-registered.
///
/// # Safety
/// Requires avx+fma, verified by the dispatcher immediately before the
/// call. The body is pure safe arithmetic on a caller-owned array.
#[cfg(target_arch = "x86_64")]
#[target_feature(enable = "avx,fma")]
unsafe fn fma_pass_x86(acc: &mut [f64; LANES], m: f64, a: f64) {
    for _ in 0..STEPS {
        for lane in acc.iter_mut() {
            *lane = lane.mul_add(m, a);
        }
    }
}

fn fma_pass(acc: &mut [f64; LANES], m: f64, a: f64) {
    #[cfg(target_arch = "x86_64")]
    if std::arch::is_x86_feature_detected!("avx") && std::arch::is_x86_feature_detected!("fma") {
        // SAFETY: features verified on this CPU immediately above.
        return unsafe { fma_pass_x86(acc, m, a) };
    }
    fma_pass_portable(acc, m, a);
}

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
}
