//! The wsbf ROOFLINE lane: SpMV effective bandwidth vs the measured
//! STREAM triad (fs-substrate), release profile, `--ignored` (perf
//! lanes run on demand / perf-CI cadence). Logical traffic accounting follows
//! the arrays the kernels actually dereference:
//! bytes = nnz·(8 val + idx bytes + 8 x[c]) plus each format's actual
//! row/chunk metadata and output accesses.
//! Attainment is LEDGERED; the >=85% acceptance
//!     gate is asserted for the best sharded kernel at the reported STREAM
//!     concurrency only (the bead's criterion; serial single-thread numbers
//!     are reported as evidence).

use std::time::Instant;

use fs_sparse::{Coo, Csr, CsrCompact, Sell};

/// Logical byte traffic for one CSR-style SpMV invocation.
///
/// The inner loop loads one value, one column index, and `x[c]` for each
/// stored nonzero. Both CSR bodies read two row bounds per output row. The
/// sharded compact path additionally performs one binary `partition_point`
/// search per internal shard boundary; its bounded row-pointer probes are
/// charged separately. This is logical source-array traffic, not a cache-miss
/// claim.
fn partition_point_row_ptr_probes(nrows: usize) -> u128 {
    if nrows == 0 {
        return 0;
    }
    // Slice partition_point delegates to binary_search_by. The current
    // standard-library algorithm deliberately makes its loop count depend
    // only on the slice length, then performs one final probe. This is
    // ceil(log2(nrows)) + 1 without overflowing at usize::MAX.
    u128::from(usize::BITS - (nrows - 1).leading_zeros()) + 1
}

fn sharded_worker_count(requested: usize, nonempty_rows: usize) -> usize {
    let host_parallelism = std::thread::available_parallelism().map_or(1, std::num::NonZero::get);
    requested
        .max(1)
        .min(nonempty_rows.max(1))
        .min(host_parallelism)
}

fn stream_concurrency_label(thread_override: bool) -> (&'static str, &'static str) {
    if thread_override {
        ("stream_selected_concurrency_gbs", "selected-concurrency")
    } else {
        ("stream_allcore_gbs", "all-core")
    }
}

fn csr_spmv_logical_bytes(
    nnz: usize,
    nrows: usize,
    index_bytes: usize,
    row_ptr_loads_per_row: usize,
    partition_searches: usize,
) -> f64 {
    let per_nnz = (2 * std::mem::size_of::<f64>() + index_bytes) as u128;
    let per_row =
        (row_ptr_loads_per_row * std::mem::size_of::<usize>() + std::mem::size_of::<f64>()) as u128;
    let partition_probes = (partition_searches as u128)
        * partition_point_row_ptr_probes(nrows)
        * std::mem::size_of::<usize>() as u128;
    (nnz as u128 * per_nnz + nrows as u128 * per_row + partition_probes) as f64
}

/// SELL metadata accesses that the chunk-major kernels perform for one matrix.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SellLogicalLayout {
    nrows: usize,
    nchunks: usize,
    row_len_loads: usize,
}

fn sell_logical_layout(a: &Csr, chunk_height: usize, sigma: usize) -> SellLogicalLayout {
    let sigma = sigma.max(chunk_height).div_ceil(chunk_height) * chunk_height;
    let mut row_lengths: Vec<usize> = (0..a.nrows()).map(|row| a.row(row).0.len()).collect();
    for window in row_lengths.chunks_mut(sigma) {
        window.sort_unstable_by(|left, right| right.cmp(left));
    }
    let nchunks = a.nrows().div_ceil(chunk_height);
    let mut row_len_loads = 0usize;
    for chunk in 0..nchunks {
        let row0 = chunk * chunk_height;
        let row1 = (row0 + chunk_height).min(a.nrows());
        let width = row_lengths[row0..row1].iter().copied().max().unwrap_or(0);
        row_len_loads += width * (row1 - row0);
    }
    SellLogicalLayout {
        nrows: a.nrows(),
        nchunks,
        row_len_loads,
    }
}

/// Logical traffic for the chunk-major SELL kernels.
///
/// Pad value/index slots are not counted because the kernel's row-length guard
/// does not dereference them. It does dereference `row_len` for every live lane
/// at every chunk width, two `chunk_ptr` entries per chunk, and one `perm` entry
/// per output. The sharded path additionally writes then reads its `(row, value)`
/// staging pair before writing `y`.
fn sell_spmv_logical_bytes(nnz: usize, layout: SellLogicalLayout, sharded: bool) -> f64 {
    let word = std::mem::size_of::<usize>() as u128;
    let value = std::mem::size_of::<f64>() as u128;
    let mut bytes = nnz as u128 * (value + word + value)
        + layout.row_len_loads as u128 * word
        + layout.nchunks as u128 * 2 * word
        + layout.nrows as u128 * (word + value);
    if sharded {
        bytes += layout.nrows as u128 * 2 * (word + value);
    }
    bytes as f64
}

/// Sharded attainment is associated with the better measured sharded kernel,
/// rather than implicitly with compact CSR.
fn best_sharded_kernel(compact_gbs: f64, sell_gbs: f64) -> (&'static str, f64) {
    if compact_gbs >= sell_gbs {
        ("compact_csr", compact_gbs)
    } else {
        ("sell_chunked", sell_gbs)
    }
}

fn banded_matrix(nrows: usize, band: usize) -> fs_sparse::Csr {
    let mut coo = Coo::new(nrows, nrows);
    let mut seed = 0xBEEF_2026_u64;
    let mut lcg = move || {
        seed = seed
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        seed
    };
    let nrows_i64 = i64::try_from(nrows).unwrap_or(i64::MAX);
    let band_i64 = i64::try_from(band).unwrap_or(i64::MAX / 8);
    let band_u64 = u64::try_from(band).unwrap_or(u64::MAX / 8);
    let window = 8_u64.saturating_mul(band_u64).max(1);
    let radius = band_i64.saturating_mul(4);
    for r in 0..nrows {
        let r_i64 = i64::try_from(r).unwrap_or(i64::MAX);
        for _ in 0..band {
            // Spread within a +-4*band window (index locality similar
            // to FEM stencils; defeats pure streaming but is honest).
            let draw = i64::try_from(lcg() % window).unwrap_or(i64::MAX);
            let off = draw - radius;
            let c = usize::try_from((r_i64 + off).clamp(0, nrows_i64 - 1)).unwrap_or(0);
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
    let requested_threads = std::env::var("FS_SPARSE_THREADS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .map(|threads| threads.max(1));
    let thread_override = requested_threads.is_some();
    let requested_threads = requested_threads
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(8, std::num::NonZero::get));
    let nrows = 4_000_000usize;
    let band = 8usize;
    // `banded_matrix` stores `band` entries in every row, so nrows is its
    // nonempty-row count. Use the same capped count for the compact kernel and
    // its STREAM denominator; the output's `threads` is actual concurrency.
    let threads = sharded_worker_count(requested_threads, nrows);
    let (stream_metric, stream_scope) = stream_concurrency_label(thread_override);
    let a = banded_matrix(nrows, band);
    let compact = CsrCompact::from_csr(&a).numa_localized(threads);
    let nnz = compact.nnz();
    let x: Vec<f64> = (0..nrows).map(|i| 0.5 + (i % 13) as f64 * 0.01).collect();
    let mut y = vec![0.0f64; nrows];
    let stream = fs_substrate::bandwidth::measure(threads);
    let time_best = |f: &mut dyn FnMut()| -> f64 {
        let mut best = f64::INFINITY;
        for _ in 0..3 {
            let t0 = Instant::now();
            f();
            best = best.min(t0.elapsed().as_secs_f64());
        }
        best
    };
    // Serial wide (usize idx), serial compact, sharded compact,
    // chunk-major SELL (serial + sharded).
    const SELL_C: usize = 8;
    const SELL_SIGMA: usize = 64;
    let sell =
        Sell::from_csr(&a, SELL_C, SELL_SIGMA).expect("fixed SELL geometry is representable");
    let sell_layout = sell_logical_layout(&a, SELL_C, SELL_SIGMA);
    let t_wide = time_best(&mut || a.spmv(&x, &mut y));
    let t_cmp = time_best(&mut || compact.spmv(&x, &mut y));
    let t_shard = time_best(&mut || compact.spmv_sharded(&x, &mut y, threads));
    let t_sell = time_best(&mut || sell.spmv_chunked(&x, &mut y));
    let t_sell_sh = time_best(&mut || sell.spmv_chunked_sharded(&x, &mut y, threads));
    std::hint::black_box(y[nrows / 2]);
    let g_wide =
        csr_spmv_logical_bytes(nnz, nrows, std::mem::size_of::<usize>(), 2, 0) / t_wide / 1e9;
    let g_cmp = csr_spmv_logical_bytes(nnz, nrows, std::mem::size_of::<u32>(), 2, 0) / t_cmp / 1e9;
    let g_shard = csr_spmv_logical_bytes(
        nnz,
        nrows,
        std::mem::size_of::<u32>(),
        2,
        threads.saturating_sub(1),
    ) / t_shard
        / 1e9;
    // SELL moves usize indices today (u32 SELL is follow-up) and has distinct
    // chunk metadata plus sharded staging traffic.
    let g_sell = sell_spmv_logical_bytes(nnz, sell_layout, false) / t_sell / 1e9;
    let g_sell_sh = sell_spmv_logical_bytes(nnz, sell_layout, true) / t_sell_sh / 1e9;
    let att_serial = g_cmp / stream.single_thread_gbs;
    let (sharded_kernel, sharded_gbs) = best_sharded_kernel(g_shard, g_sell_sh);
    let att_shard = sharded_gbs / stream.all_core_gbs;
    println!(
        "{{\"metric\":\"wsbf-roofline\",\"nrows\":{nrows},\"nnz\":{nnz},\"threads\":{threads},\
         \"stream_single_gbs\":{:.1},\"{stream_metric}\":{:.1},\
         \"spmv_wide_gbs\":{g_wide:.1},\"spmv_compact_gbs\":{g_cmp:.1},\"spmv_sharded_gbs\":{g_shard:.1},\
         \"sell_chunked_gbs\":{g_sell:.1},\"sell_sharded_gbs\":{g_sell_sh:.1},\
         \"attainment_serial\":{att_serial:.3},\"attainment_sharded_kernel\":\"{sharded_kernel}\",\
         \"attainment_sharded\":{att_shard:.3}}}",
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
             ({sharded_kernel} {sharded_gbs:.1} GB/s vs {stream_scope} STREAM {:.1} GB/s)",
            stream.all_core_gbs
        );
    }
}

#[test]
fn logical_traffic_counts_kernel_array_accesses() {
    let nnz = 3;
    let nrows = 2;
    let index_bytes = std::mem::size_of::<u32>();
    let word = std::mem::size_of::<usize>();
    let csr_expected = nnz * (8 + index_bytes + 8) + nrows * (2 * word + 8);
    let partition_probes = partition_point_row_ptr_probes(nrows) as usize;
    let shard_expected = csr_expected + partition_probes * word;
    let mut sell_coo = Coo::new(5, 5);
    sell_coo.push(0, 0, 1.0);
    sell_coo.push(0, 1, 1.0);
    sell_coo.push(1, 0, 1.0);
    // Rows 2 and 3 form a zero-width chunk; row 4 is a ragged final chunk.
    sell_coo.push(4, 4, 1.0);
    let sell_matrix = sell_coo.assemble();
    let sell = sell_logical_layout(&sell_matrix, 2, 2);
    assert_eq!(
        sell,
        SellLogicalLayout {
            nrows: 5,
            nchunks: 3,
            row_len_loads: 5,
        },
        "a zero-width chunk reads no row lengths and the final chunk checks only its live lane"
    );
    let sell_expected = 4 * (8 + word + 8) + 5 * word + 6 * word + 5 * (word + 8);
    let sell_sharded_expected = sell_expected + 5 * 2 * (word + 8);

    assert_eq!(
        csr_spmv_logical_bytes(nnz, nrows, index_bytes, 2, 0),
        csr_expected as f64
    );
    assert_eq!(
        csr_spmv_logical_bytes(nnz, nrows, index_bytes, 2, 1),
        shard_expected as f64
    );
    assert_eq!(
        csr_spmv_logical_bytes(nnz, nrows, index_bytes, 2, 2),
        (csr_expected + 2 * partition_probes * word) as f64
    );
    assert_eq!(partition_point_row_ptr_probes(0), 0);
    assert_eq!(partition_point_row_ptr_probes(1), 1);
    assert_eq!(partition_point_row_ptr_probes(2), 2);
    assert_eq!(partition_point_row_ptr_probes(3), 3);
    assert_eq!(partition_point_row_ptr_probes(4), 3);
    assert_eq!(
        sell_spmv_logical_bytes(sell_matrix.nnz(), sell, false),
        sell_expected as f64
    );
    assert_eq!(
        sell_spmv_logical_bytes(sell_matrix.nnz(), sell, true),
        sell_sharded_expected as f64
    );
    assert_eq!(best_sharded_kernel(91.0, 90.0), ("compact_csr", 91.0));
    assert_eq!(best_sharded_kernel(90.0, 91.0), ("sell_chunked", 91.0));
}

#[test]
fn roofline_stream_label_tracks_thread_override() {
    assert_eq!(
        stream_concurrency_label(false),
        ("stream_allcore_gbs", "all-core")
    );
    assert_eq!(
        stream_concurrency_label(true),
        ("stream_selected_concurrency_gbs", "selected-concurrency")
    );
}
