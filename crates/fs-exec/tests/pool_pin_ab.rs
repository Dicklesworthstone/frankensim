//! Pool-level CCD-pinning A/B (bead fz2.2 s2): the TilePool running a
//! cache-resident kernel with and without measured-topology pinning.
//! Worker w's contiguous tile range maps to buffer w, so with one
//! worker per CCD each 24 MiB working set can live in its own L3
//! island — IF the workers stay put. Report rows; the bit-invariance
//! gate lives in the pool unit tests (P2: pinning is timing-only).
//! Run: `cargo test -p fs-exec --release --test pool_pin_ab -- --ignored --nocapture`

use core::ops::ControlFlow;
use fs_exec::{Cancelled, Cx, PoolConfig, TileKernel, TilePlan, TilePool};
use fs_substrate::affinity::{CcdTopology, measured_l3_groups};
use std::time::Instant;

struct StreamKernel {
    buffers: Vec<Vec<u64>>,
    tiles_per_buf: u64,
}

impl TileKernel for StreamKernel {
    type Out = u64;

    fn tiles(&self) -> TilePlan {
        TilePlan::new(
            "fz22/stream",
            self.buffers.len() as u64 * self.tiles_per_buf,
        )
    }

    fn run(&self, tile: u64, cx: &Cx<'_>) -> ControlFlow<Cancelled, u64> {
        if cx.checkpoint().is_err() {
            return ControlFlow::Break(Cancelled);
        }
        let buf = &self.buffers[(tile / self.tiles_per_buf) as usize];
        let mut acc = 0u64;
        for &v in buf {
            acc = acc.wrapping_add(v);
        }
        ControlFlow::Continue(acc & 0xFF)
    }
}

#[test]
#[ignore = "perf harness: run explicitly in release with --ignored"]
fn pool_pinning_ab() {
    let groups = measured_l3_groups();
    let g = groups.len().max(2);
    let words = (24 << 20) / 8;
    let kernel = StreamKernel {
        buffers: (0..g)
            .map(|k| (0..words).map(|i| (i as u64) ^ (k as u64)).collect())
            .collect(),
        tiles_per_buf: 16,
    };
    let topo = CcdTopology::from_l3_groups(&groups).unwrap_or(CcdTopology::APPLE_M_CLASS);
    let best = |pool: &TilePool| -> (f64, u64) {
        let mut best = f64::INFINITY;
        let mut out = 0;
        for _ in 0..3 {
            let t0 = Instant::now();
            out = pool.run(&kernel).expect("run");
            best = best.min(t0.elapsed().as_secs_f64());
        }
        (best, out)
    };
    let unpinned = TilePool::new(PoolConfig::new(g, topo, 0xF22));
    let pinned_cfg = PoolConfig::new(g, topo, 0xF22).with_measured_pinning();
    let pinned_active = !pinned_cfg.pin_groups.is_empty();
    let pinned = TilePool::new(pinned_cfg);
    let (t_free, out_free) = best(&unpinned);
    let (t_pin, out_pin) = best(&pinned);
    assert_eq!(out_free, out_pin, "P2: pinning must never change bits");
    let bytes = (g * words * 8) as f64 * f64::from(u32::try_from(kernel.tiles_per_buf).unwrap());
    println!(
        "{{\"metric\":\"pool-pin-ab\",\"workers\":{g},\"pinning\":{pinned_active},\
         \"unpinned_gbs\":{:.1},\"pinned_gbs\":{:.1},\"speedup\":{:.2}}}",
        bytes / t_free / 1e9,
        bytes / t_pin / 1e9,
        t_free / t_pin,
    );
}
