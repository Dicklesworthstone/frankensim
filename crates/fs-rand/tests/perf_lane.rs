//! The fs-rand PERF LANE (bead 1za9): ziggurat-vs-Box-Muller and bulk
//! fill throughput, MEASURED with machine fingerprint (no "fast"
//! without a benchmark). Run explicitly in release:
//! `cargo test -p fs-rand --release --test perf_lane -- --ignored --nocapture`

use fs_rand::{Stream, StreamKey};

fn stream(seed: u64) -> Stream {
    StreamKey {
        seed,
        kernel: 0x9A2D,
        tile: 3,
    }
    .stream()
}

#[test]
#[ignore = "perf lane: run explicitly in release with --ignored"]
fn perf_normal_and_bulk() {
    let n = 10_000_000u64;
    // Box-Muller (the strict default).
    let mut st = stream(1);
    let t0 = std::time::Instant::now();
    let mut acc = 0.0f64;
    for _ in 0..n {
        acc += st.next_normal();
    }
    let bm = n as f64 / t0.elapsed().as_secs_f64();
    // Ziggurat (the fast-mode perf path).
    let mut st = stream(1);
    let t0 = std::time::Instant::now();
    for _ in 0..n {
        acc += st.next_normal_ziggurat();
    }
    let zig = n as f64 / t0.elapsed().as_secs_f64();
    // Bulk uniform fill (the SoA batch substrate).
    let mut st = stream(2);
    let mut buf = vec![0.0f64; 1 << 20];
    let t0 = std::time::Instant::now();
    for _ in 0..10 {
        st.fill_f64(&mut buf);
        acc += buf[0];
    }
    let bulk = f64::from(10 * (1 << 20)) / t0.elapsed().as_secs_f64();
    // Scalar uniform baseline for the bulk ratio.
    let mut st = stream(2);
    let t0 = std::time::Instant::now();
    for _ in 0..(10 * (1 << 20)) {
        acc += st.next_f64();
    }
    let scalar = f64::from(10 * (1 << 20)) / t0.elapsed().as_secs_f64();
    println!(
        "{{\"metric\":\"rand-perf\",\"box_muller_per_s\":{bm:.0},\"ziggurat_per_s\":{zig:.0},\
         \"zig_speedup\":{:.2},\"bulk_fill_per_s\":{bulk:.0},\"scalar_per_s\":{scalar:.0},\
         \"bulk_speedup\":{:.2},\"build\":\"release\",\"machine\":\"{}-{}\",\"sink\":{acc:.1}}}",
        zig / bm,
        bulk / scalar,
        std::env::consts::OS,
        std::env::consts::ARCH
    );
    assert!(
        zig > bm,
        "the ziggurat IS the perf path: {zig:.0} vs {bm:.0}"
    );
    // HONEST RESULT: the SoA batch is bitwise-equivalent but measured
    // ~parity with the scalar loop on this machine (0.96x) — the bulk
    // SPEEDUP claim awaits the hand-written NEON/AVX Philox capsule,
    // exactly the resource-gated no-claim the CONTRACT records. The
    // gate here only trips on pathological regression.
    assert!(
        bulk > 0.7 * scalar,
        "bulk fill stays within the parity band: {bulk:.0} vs {scalar:.0}"
    );
}
