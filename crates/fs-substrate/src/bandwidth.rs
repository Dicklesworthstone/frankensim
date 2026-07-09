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

fn triad(a: &mut [f64], b: &[f64], c: &[f64], s: f64) {
    for i in 0..a.len() {
        a[i] = b[i] + s * c[i];
    }
}

fn triad_gbs_once(len: usize) -> f64 {
    let mut a = vec![0.0f64; len];
    let b = vec![1.0f64; len];
    let c = vec![2.0f64; len];
    let mut best = 0.0f64;
    for rep in 0..REPS {
        let start = std::time::Instant::now();
        triad(&mut a, &b, &c, 1.0 + rep as f64);
        let dt = start.elapsed().as_secs_f64();
        let gbs = (len as f64 * BYTES_PER_ELEM) / dt / 1e9;
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
    let b = vec![1.0f64; LEN];
    let c = vec![2.0f64; LEN];
    let chunk = LEN.div_ceil(threads);
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
