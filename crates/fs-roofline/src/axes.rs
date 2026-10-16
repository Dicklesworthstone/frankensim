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
/// Samples per probe. The single-thread axis selects their best rate; the
/// all-core axis executes all of them inside one cumulative common window.
const PASSES: usize = 5;
const MIN_FMA_SAMPLE_SECONDS: f64 = 5e-3;
const MAX_FMA_REPETITIONS: usize = 1 << 20;
const MAX_FMA_WORKERS: usize = 4096;
const FMA_MULTIPLIER: f64 = 0.999_999_9;
const FMA_ADDEND: f64 = 1.0e-9;

use fma_kernel::fma_pass;

/// The FMA chain kernel (own registered capsule — the x86 body needs
/// `target_feature`; everything else in this file is safe).
#[path = "fma_kernel.rs"]
mod fma_kernel;

fn run_fma_repetitions(acc: &mut [f64; LANES], repetitions: usize) {
    for _ in 0..repetitions {
        fma_pass(acc, FMA_MULTIPLIER, FMA_ADDEND);
    }
}

fn calibrate_fma_repetitions(acc: &mut [f64; LANES]) -> usize {
    let mut repetitions = 1usize;
    loop {
        let start = std::time::Instant::now();
        run_fma_repetitions(acc, repetitions);
        if start.elapsed().as_secs_f64() >= MIN_FMA_SAMPLE_SECONDS
            || repetitions >= MAX_FMA_REPETITIONS
        {
            return repetitions;
        }
        repetitions *= 2;
    }
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
    let flops_per_pass = (2 * LANES * STEPS) as f64;
    let repetitions = calibrate_fma_repetitions(&mut acc);
    let mut best = 0.0f64;
    for _ in 0..PASSES {
        let start = std::time::Instant::now();
        run_fma_repetitions(&mut acc, repetitions);
        let dt = start.elapsed().as_secs_f64();
        if dt > 0.0 {
            best = best.max(flops_per_pass * repetitions as f64 / dt / 1e9);
        }
    }
    std::hint::black_box(acc);
    best
}

fn synchronized_worker_rate<F>(threads: usize, work_per_worker: f64, work: F) -> f64
where
    F: Fn(usize) + Sync,
{
    synchronized_worker_rate_with_spawn_limit(threads, work_per_worker, usize::MAX, work)
}

#[derive(Default)]
struct WorkerStartGate {
    ready: usize,
    released: bool,
    aborted: bool,
}

fn synchronized_worker_rate_with_spawn_limit<F>(
    threads: usize,
    work_per_worker: f64,
    spawn_limit: usize,
    work: F,
) -> f64
where
    F: Fn(usize) + Sync,
{
    let threads = threads.max(1);
    if threads > MAX_FMA_WORKERS || !work_per_worker.is_finite() || work_per_worker <= 0.0 {
        return f64::NAN;
    }

    let start_gate = (
        std::sync::Mutex::new(WorkerStartGate::default()),
        std::sync::Condvar::new(),
        std::sync::Condvar::new(),
    );
    let (elapsed, worker_failed) = std::thread::scope(|scope| {
        let mut handles = Vec::with_capacity(threads);
        let mut spawn_failed = false;
        for worker in 0..threads {
            if worker >= spawn_limit {
                spawn_failed = true;
                break;
            }
            let start_gate = &start_gate;
            let work = &work;
            match std::thread::Builder::new().spawn_scoped(scope, move || {
                let (state_lock, ready_changed, start_changed) = start_gate;
                let mut state = state_lock
                    .lock()
                    .unwrap_or_else(std::sync::PoisonError::into_inner);
                state.ready += 1;
                ready_changed.notify_one();
                while !state.released && !state.aborted {
                    state = start_changed
                        .wait(state)
                        .unwrap_or_else(std::sync::PoisonError::into_inner);
                }
                let aborted = state.aborted;
                drop(state);
                if !aborted {
                    work(worker);
                }
            }) {
                Ok(handle) => handles.push(handle),
                Err(_) => {
                    spawn_failed = true;
                    break;
                }
            }
        }

        let (state_lock, ready_changed, start_changed) = &start_gate;
        let mut state = state_lock
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner);
        let start = if spawn_failed {
            state.aborted = true;
            start_changed.notify_all();
            None
        } else {
            while state.ready < threads {
                state = ready_changed
                    .wait(state)
                    .unwrap_or_else(std::sync::PoisonError::into_inner);
            }
            let start = std::time::Instant::now();
            state.released = true;
            start_changed.notify_all();
            Some(start)
        };
        drop(state);

        let mut worker_failed = spawn_failed;
        for handle in handles {
            worker_failed |= handle.join().is_err();
        }
        (
            start.map_or(f64::NAN, |start| start.elapsed().as_secs_f64()),
            worker_failed,
        )
    });

    let total_work = work_per_worker * threads as f64;
    if worker_failed || !elapsed.is_finite() || elapsed <= 0.0 || !total_work.is_finite() {
        return f64::NAN;
    }
    let rate = total_work / elapsed;
    if rate.is_finite() && rate > 0.0 {
        rate
    } else {
        f64::NAN
    }
}

/// All-core FMA throughput from one synchronized wall-clock window.
///
/// A single caller-side calibration supplies the identical workload to every
/// logical worker. An abortable start gate releases the complete worker wave;
/// any creation failure or worker unwind invalidates the whole axis rather than
/// silently turning it into a partial-machine observation. Failure is exposed
/// as `NaN`; [`MachineAxes::plausibility_error`] rejects that sentinel and JSON
/// reports the unavailable axis as `null`. Requests above 4096 workers are
/// refused through the same sentinel before any worker is created.
#[must_use]
pub fn fma_gflops_all_core(logical_cpus: usize) -> f64 {
    let threads = logical_cpus.max(1);
    if threads > MAX_FMA_WORKERS {
        return f64::NAN;
    }
    let mut calibration_acc = [1.0f64; LANES];
    let repetitions = calibrate_fma_repetitions(&mut calibration_acc);
    std::hint::black_box(calibration_acc);
    let flops_per_worker = (2 * LANES * STEPS) as f64 * repetitions as f64 * PASSES as f64;

    synchronized_worker_rate(threads, flops_per_worker, |_| {
        let mut acc = [1.0f64; LANES];
        for _ in 0..PASSES {
            run_fma_repetitions(&mut acc, repetitions);
        }
        std::hint::black_box(acc);
    }) / 1e9
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
    fn all_core_rate_uses_the_staggered_workers_common_window() {
        // G3: summing independent 5/15/45 ms worker rates would imply a
        // much shorter aggregate window than the slowest concurrent worker.
        let worker_count = 3usize;
        let visits: [std::sync::atomic::AtomicUsize; 3] =
            std::array::from_fn(|_| std::sync::atomic::AtomicUsize::new(0));
        let rate = synchronized_worker_rate(worker_count, 1.0, |worker| {
            visits[worker].fetch_add(1, std::sync::atomic::Ordering::SeqCst);
            let delay_ms = [5, 15, 45][worker];
            std::thread::sleep(std::time::Duration::from_millis(delay_ms));
        });
        let inferred_common_seconds = worker_count as f64 / rate;
        assert!(
            inferred_common_seconds >= 0.035,
            "aggregate rate did not retain the slowest worker's common window: {inferred_common_seconds}s"
        );
        assert!(
            visits
                .iter()
                .all(|count| count.load(std::sync::atomic::Ordering::SeqCst) == 1),
            "each requested worker must execute the shared workload exactly once"
        );
    }

    #[test]
    fn all_core_rate_is_unavailable_after_any_worker_panics() {
        let drained = std::sync::atomic::AtomicUsize::new(0);
        let rate = synchronized_worker_rate(3, 1.0, |worker| {
            assert_ne!(worker, 1, "synthetic all-core worker failure");
            if worker == 2 {
                std::thread::sleep(std::time::Duration::from_millis(20));
            }
            drained.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
        });
        assert!(rate.is_nan(), "partial-machine rate escaped: {rate}");
        assert_eq!(drained.load(std::sync::atomic::Ordering::SeqCst), 2);
    }

    #[test]
    fn all_core_spawn_refusal_aborts_and_drains_the_waiting_wave() {
        let executed = std::sync::atomic::AtomicUsize::new(0);
        let rate = synchronized_worker_rate_with_spawn_limit(3, 1.0, 1, |_| {
            executed.fetch_add(1, std::sync::atomic::Ordering::SeqCst);
        });
        assert!(rate.is_nan(), "incomplete worker wave escaped: {rate}");
        assert_eq!(executed.load(std::sync::atomic::Ordering::SeqCst), 0);
        assert!(synchronized_worker_rate(MAX_FMA_WORKERS + 1, 1.0, |_| {}).is_nan());
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

        let mut unavailable_worker_wave = no_cpus;
        unavailable_worker_wave.logical_cpus = 8;
        unavailable_worker_wave.peak_all_core_gflops = f64::NAN;
        assert!(!unavailable_worker_wave.plausible());
        assert!(
            unavailable_worker_wave
                .to_jsonl()
                .contains("\"peak_all_core_gflops\":null")
        );
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
