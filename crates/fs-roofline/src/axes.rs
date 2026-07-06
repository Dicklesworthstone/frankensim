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

/// Independent FMA accumulator lanes per chain step (wide enough to fill
/// 512-bit units with room for dual issue).
const LANES: usize = 64;
/// Chain steps per timed pass.
const STEPS: usize = 4096;
/// Timed passes per measurement; best-of keeps thermal/scheduler noise from
/// deflating the axis (an axis too low would inflate everyone's attainment).
const PASSES: usize = 5;

fn fma_pass(acc: &mut [f64; LANES], m: f64, a: f64) {
    for _ in 0..STEPS {
        for lane in acc.iter_mut() {
            *lane = lane.mul_add(m, a);
        }
    }
}

/// Best-of-passes single-thread FMA throughput in GFLOP/s.
#[must_use]
pub fn fma_gflops_single() -> f64 {
    let mut acc = [1.0f64; LANES];
    // Multiplier chosen so accumulators orbit without overflow/denormal.
    let (m, a) = (0.999_999_9, 1.0e-9);
    let flops = (2 * LANES * STEPS) as f64;
    let mut best = 0.0f64;
    for _ in 0..PASSES {
        let start = std::time::Instant::now();
        fma_pass(&mut acc, m, a);
        let dt = start.elapsed().as_secs_f64();
        if dt > 0.0 {
            best = best.max(flops / dt / 1e9);
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
