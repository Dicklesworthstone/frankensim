//! The wsbf ROOFLINE lane: SpMV effective bandwidth vs the measured
//! STREAM triad (fs-substrate), release profile, `--ignored` (perf
//! lanes run on demand / perf-CI cadence). Conventional accounting:
//! bytes = nnz·(8 val + idx bytes) + nrows·8 (y) + nrows·8 (row_ptr)
//! + ncols·8 (x once). Attainment is LEDGERED; the >=85% acceptance
//! gate is asserted for the SHARDED all-core kernel only (the bead's
//! criterion; serial single-thread numbers are reported as evidence).

use std::time::Instant;

use fs_sparse::{Coo, CsrCompact};

fn banded_matrix(nrows: usize, band: usize) -> fs_sparse::Csr {
    let mut coo = Coo::new(nrows, nrows);
    let mut seed = 0xBEEF_2026_u64;
    let mut lcg = move || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        seed
    };
    for r in 0..nrows {
        for k in 0..band {
            // Spread within a +-4*band window (index locality similar
            // to FEM stencils; defeats pure streaming but is honest).
            let off = (lcg() % (8 * band as u64)) as i64 - 4 * band as i64;
            let c = (r as i64 + off).clamp(0, nrows as i64 - 1) as usize;
            let v = ((lcg() >> 11) as f64) / (1u64 << 53) as f64 + 0.5;
            coo.push(r, c, v);
        }
    }
    coo.assemble()
}

#[test]
#[ignore = "perf lane: run in release on demand (mac + ts1); nightly cadence is fz2.4"]
fn wsbf_roofline() {
    // FS_SPARSE_THREADS overrides (heterogeneous-core machines: equal
    // nnz shards let E-cores drag the tail — pin to P-core count).
    let threads = std::env::var("FS_SPARSE_THREADS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(8, std::num::NonZero::get));
    let nrows = 4_000_000usize;
    let band = 8usize;
    let a = banded_matrix(nrows, band);
    let compact = CsrCompact::from_csr(&a).numa_localized(threads);
    let nnz = compact.nnz();
    let x: Vec<f64> = (0..nrows).map(|i| 0.5 + (i % 13) as f64 * 0.01).collect();
    let mut y = vec![0.0f64; nrows];
    let stream = fs_substrate::bandwidth::measure(threads);
    let bytes = |idx_bytes: usize| -> f64 {
        (nnz * (8 + idx_bytes) + nrows * 8 + nrows * 8 + nrows * 8) as f64
    };
    let time_best = |f: &mut dyn FnMut()| -> f64 {
        let mut best = f64::INFINITY;
        for _ in 0..3 {
            let t0 = Instant::now();
            f();
            best = best.min(t0.elapsed().as_secs_f64());
        }
        best
    };
    // Serial wide (usize idx), serial compact, sharded compact.
    let t_wide = time_best(&mut || a.spmv(&x, &mut y));
    let t_cmp = time_best(&mut || compact.spmv(&x, &mut y));
    let t_shard = time_best(&mut || compact.spmv_sharded(&x, &mut y, threads));
    std::hint::black_box(y[nrows / 2]);
    let g_wide = bytes(8) / t_wide / 1e9;
    let g_cmp = bytes(4) / t_cmp / 1e9;
    let g_shard = bytes(4) / t_shard / 1e9;
    let att_serial = g_cmp / stream.single_thread_gbs;
    let att_shard = g_shard / stream.all_core_gbs;
    println!(
        "{{\"metric\":\"wsbf-roofline\",\"nrows\":{nrows},\"nnz\":{nnz},\"threads\":{threads},\
         \"stream_single_gbs\":{:.1},\"stream_allcore_gbs\":{:.1},\
         \"spmv_wide_gbs\":{g_wide:.1},\"spmv_compact_gbs\":{g_cmp:.1},\"spmv_sharded_gbs\":{g_shard:.1},\
         \"attainment_serial\":{att_serial:.3},\"attainment_sharded\":{att_shard:.3}}}",
        stream.single_thread_gbs, stream.all_core_gbs
    );
    // The 85% acceptance GATE asserts under FS_SPARSE_ROOFLINE_GATE=1
    // (the perf-CI lanes / dedicated machines); ad-hoc runs on loaded
    // dev boxes REPORT — a hard gate there measures the neighbors'
    // builds, not the kernel (mac numbers swung 25% run-to-run while
    // the swarm compiled).
    if std::env::var("FS_SPARSE_ROOFLINE_GATE").as_deref() == Ok("1") {
        assert!(
            att_shard >= 0.85,
            "sharded SpMV attainment {att_shard:.3} below the 85% STREAM gate \
             (sharded {g_shard:.1} GB/s vs all-core triad {:.1} GB/s)",
            stream.all_core_gbs
        );
    }
}
