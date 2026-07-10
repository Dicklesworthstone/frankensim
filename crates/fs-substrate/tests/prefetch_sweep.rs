//! Prefetch-distance sweep (bead fz2.2): data-DEPENDENT gather walk —
//! the pattern hardware prefetchers cannot predict — with software
//! read-ahead at candidate distances. Reports GB/s per distance and
//! the per-machine winner; the winning distance is the value the
//! tuner records. Report rows (P2: prefetch changes timing, never
//! bits — asserted).
//! Run: `cargo test -p fs-substrate --release --test prefetch_sweep -- --ignored --nocapture`

use fs_substrate::prefetch::read_ahead;
use std::time::Instant;

#[test]
#[ignore = "perf harness: run explicitly in release with --ignored"]
fn gather_prefetch_distance_sweep() {
    // 256 MiB data, 8M pseudo-random gather indices (LCG, deterministic).
    let n = (256usize << 20) / 8;
    let data: Vec<u64> = (0..n)
        .map(|i| (i as u64).wrapping_mul(0x9E37_79B9))
        .collect();
    let m = 8usize << 20;
    let mut seed = 0x5EED_5EEDu64;
    let idx: Vec<u32> = (0..m)
        .map(|_| {
            seed = seed
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((seed >> 33) as usize % n) as u32
        })
        .collect();
    let walk = |dist: usize| -> (f64, u64) {
        let mut best = f64::INFINITY;
        let mut acc_out = 0u64;
        for _ in 0..3 {
            let t0 = Instant::now();
            let mut acc = 0u64;
            for i in 0..m {
                if dist > 0
                    && let Some(&pf) = idx.get(i + dist)
                {
                    read_ahead(&data, pf as usize);
                }
                acc = acc.wrapping_add(data[idx[i] as usize]);
            }
            std::hint::black_box(acc);
            acc_out = acc;
            best = best.min(t0.elapsed().as_secs_f64());
        }
        ((m * 8) as f64 / best / 1e9, acc_out)
    };
    let (base_gbs, base_acc) = walk(0);
    let mut rows = vec![(0usize, base_gbs)];
    let mut winner = (0usize, base_gbs);
    for dist in [2usize, 4, 8, 16, 32, 64, 128, 256] {
        let (gbs, acc) = walk(dist);
        assert_eq!(acc, base_acc, "P2: prefetch must never change bits");
        rows.push((dist, gbs));
        if gbs > winner.1 {
            winner = (dist, gbs);
        }
    }
    for (dist, gbs) in &rows {
        println!("{{\"metric\":\"prefetch-sweep\",\"distance\":{dist},\"gather_gbs\":{gbs:.2}}}");
    }
    println!(
        "{{\"metric\":\"prefetch-winner\",\"distance\":{},\"gather_gbs\":{:.2},\"over_no_prefetch\":{:.2}}}",
        winner.0,
        winner.1,
        winner.1 / base_gbs.max(1e-9)
    );
}
