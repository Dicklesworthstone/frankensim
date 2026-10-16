//! Prefetch-distance sweep (bead fz2.2): data-DEPENDENT gather walk —
//! the pattern hardware prefetchers cannot predict — with software
//! read-ahead at candidate distances. Reports GB/s per distance and
//! the per-machine winner; the winning distance is the value the
//! tuner records. Report rows (P2: prefetch changes timing, never
//! bits — asserted).
//! Run: `cargo test -p fs-substrate --release --test prefetch_sweep -- --ignored --nocapture`

use fs_substrate::prefetch::read_ahead;
use std::time::Instant;

const SUITE: &str = "fs-substrate/prefetch-sweep";
const GATHER_INPUT_SEED: u64 = 0x5EED_5EED;

fn measurement(identity: &str, name: &str, json: String) {
    let mut emitter = fs_obs::Emitter::new(SUITE, identity);
    let event = emitter.emit(
        fs_obs::Severity::Info,
        fs_obs::EventKind::Custom {
            name: name.to_string(),
            json,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("prefetch sweep row must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("prefetch sweep row must use the fs-obs wire schema");
    println!("{line}");
}

fn finite_json(value: f64, precision: usize) -> String {
    if value.is_finite() {
        format!("{value:.precision$}")
    } else {
        "null".to_string()
    }
}

#[test]
#[ignore = "perf harness: run explicitly in release with --ignored"]
fn gather_prefetch_distance_sweep() {
    // 256 MiB data, 8M pseudo-random gather indices (LCG, deterministic).
    let n = (256usize << 20) / 8;
    let data: Vec<u64> = (0..n)
        .map(|i| (i as u64).wrapping_mul(0x9E37_79B9))
        .collect();
    let m = 8usize << 20;
    let mut seed = GATHER_INPUT_SEED;
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
    let machine = fs_substrate::CapabilityProbe::topology_only().fingerprint();
    for (dist, gbs) in &rows {
        measurement(
            &format!("prefetch-sweep/distance-{dist}/measurement"),
            "prefetch-sweep",
            format!(
                "{{\"metric\":\"prefetch-sweep\",\"distance\":{dist},\"gather_gbs\":{},\
                 \"machine\":{machine},\"input_seed\":{GATHER_INPUT_SEED},\"trials\":3,\
                 \"timing_seed\":null,\"timing_replay\":false}}",
                finite_json(*gbs, 2),
            ),
        );
    }
    let over_no_prefetch = winner.1 / base_gbs.max(1e-9);
    measurement(
        "prefetch-winner/measurement",
        "prefetch-winner",
        format!(
            "{{\"metric\":\"prefetch-winner\",\"distance\":{},\"gather_gbs\":{},\
             \"over_no_prefetch\":{},\"machine\":{machine},\
             \"input_seed\":{GATHER_INPUT_SEED},\"trials\":3,\"timing_seed\":null,\
             \"timing_replay\":false}}",
            winner.0,
            finite_json(winner.1, 2),
            finite_json(over_no_prefetch, 2),
        ),
    );
}
