//! Sustained-bandwidth measurement: a STREAM-triad-style sweep
//! (`a[i] = b[i] + s * c[i]`) over buffers far larger than the last-level
//! cache, best-of-N repetitions, single-thread and all-core variants.
//!
//! These numbers are the roofline's memory axis (plan §14.1) and the
//! bandwidth-rich vs bandwidth-starved schedule decision input (§5.1
//! consequence 2: ~34 GB/s/core on M-class vs ~3.4 GB/s/core on 96-core
//! parts — the SAME kernel needs different schedules).
//!
//! Honest limits (no-claim until the executor exists): per-CORE-CLASS
//! bandwidth (P vs E pinning) needs QoS/affinity APIs outside safe std;
//! we measure single-thread and all-core aggregates and record per-class
//! CORE COUNTS from the topology probe. Timing uses `std::time::Instant`
//! (measurement only — never on a determinism-bearing path).

use super::Measured;

/// Bytes moved per triad element: read b, read c, write a (allocate-on-write
/// traffic varies by microarchitecture; 24 B/elem is the conventional
/// STREAM accounting and is what we report against).
const BYTES_PER_ELEM: f64 = 24.0;

/// Buffer length per array (64 MiB of f64 per array ≫ any current LLC).
const LEN: usize = 8 << 20;

/// Repetitions; best (max) is reported, standard STREAM practice.
const REPS: usize = 3;

/// Single-thread timed samples (best-of; more than the all-core REPS
/// because this axis divides single-thread attainments directly).
const SINGLE_SAMPLES: usize = 5;

/// Four AVX2 vectors (or eight Neon vectors) per source-level body. The fixed
/// body exposes independent loads/stores to the optimizer instead of relying
/// on target-specific vectorization of one indexed scalar statement.
const TRIAD_UNROLL: usize = 16;

fn triad(a: &mut [f64], b: &[f64], c: &[f64], s: f64) {
    debug_assert_eq!(a.len(), b.len());
    debug_assert_eq!(a.len(), c.len());

    let unrolled_len = a.len() / TRIAD_UNROLL * TRIAD_UNROLL;
    let (a_head, a_tail) = a.split_at_mut(unrolled_len);
    let (b_head, b_tail) = b.split_at(unrolled_len);
    let (c_head, c_tail) = c.split_at(unrolled_len);
    for ((a_chunk, b_chunk), c_chunk) in a_head
        .chunks_exact_mut(TRIAD_UNROLL)
        .zip(b_head.chunks_exact(TRIAD_UNROLL))
        .zip(c_head.chunks_exact(TRIAD_UNROLL))
    {
        a_chunk[0] = b_chunk[0] + s * c_chunk[0];
        a_chunk[1] = b_chunk[1] + s * c_chunk[1];
        a_chunk[2] = b_chunk[2] + s * c_chunk[2];
        a_chunk[3] = b_chunk[3] + s * c_chunk[3];
        a_chunk[4] = b_chunk[4] + s * c_chunk[4];
        a_chunk[5] = b_chunk[5] + s * c_chunk[5];
        a_chunk[6] = b_chunk[6] + s * c_chunk[6];
        a_chunk[7] = b_chunk[7] + s * c_chunk[7];
        a_chunk[8] = b_chunk[8] + s * c_chunk[8];
        a_chunk[9] = b_chunk[9] + s * c_chunk[9];
        a_chunk[10] = b_chunk[10] + s * c_chunk[10];
        a_chunk[11] = b_chunk[11] + s * c_chunk[11];
        a_chunk[12] = b_chunk[12] + s * c_chunk[12];
        a_chunk[13] = b_chunk[13] + s * c_chunk[13];
        a_chunk[14] = b_chunk[14] + s * c_chunk[14];
        a_chunk[15] = b_chunk[15] + s * c_chunk[15];
    }
    for ((out, left), right) in a_tail.iter_mut().zip(b_tail).zip(c_tail) {
        *out = *left + s * *right;
    }
}

fn triad_gbs_once(len: usize) -> f64 {
    let mut a = vec![0.0f64; len];
    let b = vec![1.0f64; len];
    let c = vec![2.0f64; len];
    // Touch `a` before timing: a fresh zeroed vec is lazily mapped, so
    // the first sweep would pay every page fault. Then CALIBRATE the
    // sweeps-per-sample to span ≥ 5 ms of wall clock: single ~2 ms
    // sweeps sat inside the scheduler-placement noise floor (one
    // E-core window halved the axis run-to-run — measured 55 vs 89
    // GB/s on the same M4 Pro, flipping downstream roofline verdicts).
    triad(&mut a, &b, &c, 0.5);
    let mut sweeps = 1usize;
    loop {
        let start = std::time::Instant::now();
        for k in 0..sweeps {
            triad(&mut a, &b, &c, 1.0 + k as f64);
        }
        if start.elapsed().as_secs_f64() >= 5e-3 || sweeps >= 1 << 12 {
            break;
        }
        sweeps *= 2;
    }
    let mut best = 0.0f64;
    for rep in 0..SINGLE_SAMPLES {
        let start = std::time::Instant::now();
        for k in 0..sweeps {
            triad(&mut a, &b, &c, 1.0 + (rep * sweeps + k) as f64);
        }
        let dt = start.elapsed().as_secs_f64();
        let gbs = (len as f64 * BYTES_PER_ELEM) * sweeps as f64 / dt / 1e9;
        best = best.max(gbs);
    }
    // Defeat dead-store elimination.
    std::hint::black_box(a[len / 2]);
    best
}

/// Measure single-thread and all-core sustained triad bandwidth.
#[must_use]
pub fn measure(logical_cpus: usize) -> Measured {
    let single = triad_gbs_once(LEN);
    // All-core: threads stream DISJOINT SLICES of one shared triple of
    // full-size arrays (classic parallel STREAM). MEASURED rejection of
    // the previous design (private per-thread buffers of LEN/threads):
    // at 14 threads each buffer fit in cache and the "bandwidth" read
    // 922 GB/s on an M4 Pro whose DRAM peak is ~273 — the roofline
    // denominator was measuring SLC, not memory (found by the wsbf
    // SpMV attainment gate reading 8%).
    let threads = logical_cpus.max(1);
    let mut a = vec![0.0f64; LEN];
    let mut b = vec![0.0f64; LEN];
    let mut c = vec![0.0f64; LEN];
    let chunk = LEN.div_ceil(threads);
    // FIRST-TOUCH: each timing thread initializes ITS chunk, so pages
    // land on the toucher's NUMA node. MEASURED rejection of serial
    // init: on a 64-core Threadripper (8-channel DDR4) all-core triad
    // read 30 GB/s — every page faulted on one node and 64 threads
    // queued on one memory controller.
    std::thread::scope(|s| {
        let mut ra = a.as_mut_slice();
        let mut rb = b.as_mut_slice();
        let mut rc = c.as_mut_slice();
        while !ra.is_empty() {
            let take = chunk.min(ra.len());
            let (ma, ta) = ra.split_at_mut(take);
            let (mb, tb) = rb.split_at_mut(take);
            let (mc, tc) = rc.split_at_mut(take);
            s.spawn(move || {
                ma.fill(0.0);
                mb.fill(1.0);
                mc.fill(2.0);
            });
            ra = ta;
            rb = tb;
            rc = tc;
        }
    });
    let (b, c) = (b, c);
    let mut best = 0.0f64;
    for rep in 0..REPS {
        let scale = 1.0 + rep as f64;
        let start = std::time::Instant::now();
        std::thread::scope(|s| {
            let mut rest = a.as_mut_slice();
            let mut off = 0usize;
            while !rest.is_empty() {
                let take = chunk.min(rest.len());
                let (mine, tail) = rest.split_at_mut(take);
                let (bs, cs) = (&b[off..off + take], &c[off..off + take]);
                s.spawn(move || triad(mine, bs, cs, scale));
                rest = tail;
                off += take;
            }
        });
        let dt = start.elapsed().as_secs_f64();
        best = best.max((LEN as f64 * BYTES_PER_ELEM) / dt / 1e9);
    }
    std::hint::black_box(a[LEN / 2]);
    Measured {
        single_thread_gbs: single,
        all_core_gbs: best,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn unrolled_triad_matches_scalar_reference_including_tails() {
        for len in [0, 1, 7, 15, 16, 17, 31, 32, 33, 257] {
            let b: Vec<f64> = (0..len).map(|i| i as f64 * 0.25 - 7.0).collect();
            let c: Vec<f64> = (0..len).map(|i| 11.0 - i as f64 * 0.125).collect();
            let scale = -1.25;
            let expected: Vec<f64> = b
                .iter()
                .zip(&c)
                .map(|(left, right)| *left + scale * *right)
                .collect();
            let mut actual = vec![f64::NAN; len];
            triad(&mut actual, &b, &c, scale);
            assert_eq!(actual, expected, "tail mismatch at len={len}");
        }
    }

    #[test]
    fn bandwidth_is_physically_plausible() {
        // Small, fast variant for CI: 8 MiB per array.
        let single = triad_gbs_once(1 << 20);
        // Wide bounds ON PURPOSE: this guards against ACCOUNTING bugs
        // (wrong BYTES_PER_ELEM, ms-vs-s confusion produce 1000x errors),
        // not against slow machines — a debug build on a loaded 128-thread
        // box legitimately measures under 1 GB/s (learned on trj, load 27).
        assert!(
            (0.01..=20_000.0).contains(&single),
            "single-thread triad {single} GB/s outside sanity bounds (accounting bug?)"
        );
    }

    #[test]
    fn all_core_measures_and_exceeds_zero() {
        let m = measure(2);
        assert!(m.single_thread_gbs > 0.0);
        assert!(m.all_core_gbs > 0.0);
        println!(
            "{{\"suite\":\"fs-substrate/bandwidth\",\"case\":\"triad\",\"verdict\":\"pass\",\"detail\":\"single={:.1}GB/s all={:.1}GB/s\"}}",
            m.single_thread_gbs, m.all_core_gbs
        );
    }
}
