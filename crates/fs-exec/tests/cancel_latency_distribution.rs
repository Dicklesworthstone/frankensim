//! Cancel-latency distribution measurement over >= 10,000 trials (q61wp.28).
//!
//! Measures request -> drain-complete latency distribution (p50/p90/p99/max)
//! for tile kernels under `Cx` polling at tile boundaries across warm worker threads.

use core::ops::ControlFlow;
use fs_exec::{
    CancelGate, Cancelled, Cx, PoolConfig, RunError, TileKernel, TilePlan, TilePool,
};
struct LatencyProbeKernel<'a> {
    gate: &'a CancelGate,
    trigger_tile: u64,
    ops_per_tile: usize,
    tile_count: u64,
}

impl TileKernel for LatencyProbeKernel<'_> {
    type Out = u64;

    fn tiles(&self) -> TilePlan {
        TilePlan::new("perf/cancel-latency-probe", self.tile_count)
    }

    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, u64> {
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(Cancelled);
        }
        if tile == self.trigger_tile {
            self.gate.request();
        }
        let mut acc = tile.wrapping_add(1);
        for i in 0..self.ops_per_tile {
            acc = acc.wrapping_mul(6364136223846793005).wrapping_add(i as u64);
        }
        std::hint::black_box(acc);
        ControlFlow::Continue(1)
    }
}

/// Simple LCG for deterministic pseudo-random triggers without extra dependencies.
struct SimpleLcg(u64);
impl SimpleLcg {
    fn new(seed: u64) -> Self {
        Self(seed)
    }
    fn next_u64(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.0
    }
    fn range(&mut self, min: u64, max: u64) -> u64 {
        min + (self.next_u64() >> 16) % (max - min + 1)
    }
}

/// Distribution statistics for cancel latency measurements.
#[derive(Debug, Clone)]
pub struct LatencyDistribution {
    /// Number of measurement trials.
    pub trials: usize,
    /// 50th percentile latency in nanoseconds.
    pub p50_ns: u64,
    /// 90th percentile latency in nanoseconds.
    pub p90_ns: u64,
    /// 95th percentile latency in nanoseconds.
    pub p95_ns: u64,
    /// 99th percentile latency in nanoseconds.
    pub p99_ns: u64,
    /// 99.9th percentile latency in nanoseconds.
    pub p999_ns: u64,
    /// Maximum observed latency in nanoseconds.
    pub max_ns: u64,
    /// Minimum observed latency in nanoseconds.
    pub min_ns: u64,
    /// Arithmetic mean latency in nanoseconds.
    pub mean_ns: u64,
    /// Operations per tile during measurement.
    pub tile_ops: usize,
}

impl LatencyDistribution {
    /// Compute percentile distribution from raw sample latencies.
    pub fn compute(mut samples: Vec<u64>, tile_ops: usize) -> Self {
        assert!(!samples.is_empty(), "samples must not be empty");
        samples.sort_unstable();
        let n = samples.len();
        let sum: u128 = samples.iter().map(|&s| s as u128).sum();
        let mean_ns = (sum / n as u128) as u64;

        Self {
            trials: n,
            p50_ns: samples[(n as f64 * 0.50) as usize],
            p90_ns: samples[(n as f64 * 0.90) as usize],
            p95_ns: samples[(n as f64 * 0.95) as usize],
            p99_ns: samples[((n as f64 * 0.99) as usize).min(n - 1)],
            p999_ns: samples[((n as f64 * 0.999) as usize).min(n - 1)],
            max_ns: samples[n - 1],
            min_ns: samples[0],
            mean_ns,
            tile_ops,
        }
    }
}

#[test]
fn measure_cancel_latency_distribution_10k_trials() {
    const TRIALS: usize = 10_000;
    const WORKERS: usize = 4;
    const TILE_OPS: usize = 500; // ~1-3 microseconds on modern cores

    let pool = TilePool::new(PoolConfig::for_host(WORKERS, 0xca_4c_e1));
    let mut rng = SimpleLcg::new(0x2026_0903);

    pool.with_parked_crew_local(|parked| {
        let mut drain_samples_ns = Vec::with_capacity(TRIALS);
        let mut worker_observed_p99_samples_ns = Vec::with_capacity(TRIALS);

        for _ in 0..TRIALS {
            let gate = CancelGate::new();
            let trigger_tile = rng.range(2, 10);
            let kernel = LatencyProbeKernel {
                gate: &gate,
                trigger_tile,
                ops_per_tile: TILE_OPS,
                tile_count: 64,
            };

            let (result, report) = parked.run_with_gate(&kernel, &gate);
            let drain_finished_ns = gate.now_ns();

            assert!(
                matches!(result, Err(RunError::Cancelled { .. })),
                "run must be cancelled"
            );

            if let Some(req_at) = gate.requested_at_ns() {
                let drain_latency = drain_finished_ns.saturating_sub(req_at);
                drain_samples_ns.push(drain_latency);
            }
            if let Some(p99_obs) = report.cancel_latency_p99_ns() {
                worker_observed_p99_samples_ns.push(p99_obs);
            }
        }

        let dist = LatencyDistribution::compute(drain_samples_ns, TILE_OPS);
        let obs_dist = LatencyDistribution::compute(worker_observed_p99_samples_ns, TILE_OPS);

        println!("=== CANCEL LATENCY OVER {} TRIALS (IN-BAND TRIGGER) ===", dist.trials);
        println!("Request -> Drain-Complete Latency:");
        println!("  min:    {:>8.2} µs ({} ns)", dist.min_ns as f64 / 1_000.0, dist.min_ns);
        println!("  p50:    {:>8.2} µs ({} ns)", dist.p50_ns as f64 / 1_000.0, dist.p50_ns);
        println!("  p90:    {:>8.2} µs ({} ns)", dist.p90_ns as f64 / 1_000.0, dist.p90_ns);
        println!("  p95:    {:>8.2} µs ({} ns)", dist.p95_ns as f64 / 1_000.0, dist.p95_ns);
        println!("  p99:    {:>8.2} µs ({} ns)", dist.p99_ns as f64 / 1_000.0, dist.p99_ns);
        println!("  p99.9:  {:>8.2} µs ({} ns)", dist.p999_ns as f64 / 1_000.0, dist.p999_ns);
        println!("  max:    {:>8.2} µs ({} ns)", dist.max_ns as f64 / 1_000.0, dist.max_ns);
        println!("  mean:   {:>8.2} µs ({} ns)", dist.mean_ns as f64 / 1_000.0, dist.mean_ns);
        println!("Worker-Observed p99 Latency Distribution:");
        println!("  p50:    {:>8.2} µs ({} ns)", obs_dist.p50_ns as f64 / 1_000.0, obs_dist.p50_ns);
        println!("  p99:    {:>8.2} µs ({} ns)", obs_dist.p99_ns as f64 / 1_000.0, obs_dist.p99_ns);
        println!("  max:    {:>8.2} µs ({} ns)", obs_dist.max_ns as f64 / 1_000.0, obs_dist.max_ns);
    });
}
