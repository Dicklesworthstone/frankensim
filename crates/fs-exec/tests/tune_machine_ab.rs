//! Machine calibration A/B lane (bead fz2.2 s3): the REAL probe +
//! calibration report for this machine (the ledgered tuning decisions),
//! plus the quantum-weights verification — class-shaped initial shares
//! vs pure work-stealing, measured. Report rows; the selection DOCTRINE
//! (same kernels, different winner per machine class) is unit-gated in
//! tune.rs with synthetic probes.
//! Run: `cargo test -p fs-exec --release --test tune_machine_ab -- --ignored --nocapture`

use core::ops::ControlFlow;
use fs_exec::{CancelGate, Cancelled, Cx, PoolConfig, TileKernel, TilePlan, TilePool, Tuner};
use fs_substrate::CapabilityProbe;
use fs_substrate::affinity::CcdTopology;
use std::time::Instant;

/// Uniform CPU-bound tiles (the same shape as the tuner's class probe).
struct UniformKernel;

impl TileKernel for UniformKernel {
    type Out = u64;

    fn tiles(&self) -> TilePlan {
        TilePlan::new("fz22/uniform", 1024)
    }

    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, u64> {
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(Cancelled);
        }
        let mut acc = (tile as f64).mul_add(1e-9, 1.0);
        for _ in 0..500_000 {
            acc = acc.mul_add(0.999_999_9, 1.0e-9);
        }
        ControlFlow::Continue(acc.to_bits() & 1)
    }
}

#[test]
#[ignore = "perf harness: run explicitly in release with --ignored"]
fn machine_calibration_and_weight_verification() {
    // --- The real calibration report for THIS machine (ledger row). ---
    let probe = CapabilityProbe::run();
    let mut tuner = Tuner::cold(probe.fingerprint());
    let report = tuner
        .calibrate(&probe)
        .expect("probe matches tuner fingerprint");
    println!("{report}");
    let (kind, _) = tuner.schedule();
    println!(
        "{{\"metric\":\"schedule-selected\",\"cpu\":\"{}\",\"schedule\":\"{}\"}}",
        probe.cpu_brand,
        kind.name()
    );

    // --- Quantum-weights verification (fz2.2): does shaping initial
    // shares from the measured class distribution beat pure stealing?
    // Measured honestly either way — a tie is the "stealing already
    // absorbs the asymmetry" verdict, not a failure. ---
    let workers = (probe.logical_cpus as usize).max(1);
    let topo = CcdTopology::from_probe(&probe);
    let gate = CancelGate::new();
    let probe_pool = TilePool::new(PoolConfig::new(workers, topo, 0xF22));
    let (r, class_report) = probe_pool.run_with_gate(&UniformKernel, &gate);
    r.expect("class measurement runs");
    let mut counts = class_report.tiles_by_worker.clone();
    counts.sort_unstable_by(|a, b| b.cmp(a));
    let fast = counts.first().copied().unwrap_or(1).max(1);
    // Weights 1..=4 scaled by measured per-worker throughput; the slow
    // tail gets proportionally smaller initial quanta.
    let weights: Vec<u32> = counts
        .iter()
        .map(|&c| u32::try_from((4 * c).div_ceil(fast).max(1)).unwrap_or(1))
        .collect();
    println!(
        "{{\"metric\":\"class-distribution\",\"workers\":{workers},\"tiles_sorted\":{counts:?},\"weights\":{weights:?}}}"
    );
    let best = |cfg: PoolConfig| -> (f64, u64, u64) {
        let pool = TilePool::new(cfg);
        let mut best = f64::INFINITY;
        let mut steals = 0;
        let mut out = 0;
        for _ in 0..3 {
            let gate = CancelGate::new();
            let t0 = Instant::now();
            let (r, rep) = pool.run_with_gate(&UniformKernel, &gate);
            out = r.expect("run");
            let dt = t0.elapsed().as_secs_f64();
            if dt < best {
                best = dt;
                steals = rep.steals;
            }
        }
        (best, steals, out)
    };
    let (t_eq, s_eq, out_eq) = best(PoolConfig::new(workers, topo, 0xF22));
    let (t_w, s_w, out_w) = best(PoolConfig {
        quantum_weights: weights,
        ..PoolConfig::new(workers, topo, 0xF22)
    });
    assert_eq!(out_eq, out_w, "P2: weights must never change bits");
    println!(
        "{{\"metric\":\"weight-ab\",\"equal_ms\":{:.1},\"equal_steals\":{s_eq},\
         \"weighted_ms\":{:.1},\"weighted_steals\":{s_w},\"speedup\":{:.3}}}",
        t_eq * 1e3,
        t_w * 1e3,
        t_eq / t_w,
    );
}
