//! The fs-mesh PERF LANE (bead uee3, item 5): large-point-count
//! Delaunay throughput, MEASURED and ledgered with a machine
//! fingerprint per the perf program (no "fast" without a benchmark).
//!
//! Run explicitly, in release:
//! `cargo test -p fs-mesh --release --test perf_lane -- --ignored --nocapture`
//! The 1e7 rung additionally wants `FS_MESH_PERF_FULL=1` (memory and
//! minutes); the default invocation measures 1e5 and 1e6.
//!
//! Correctness rides along: the exact structural audit must stay clean
//! at every rung (full insphere certification at the smallest rung
//! only — it is O(n·tets) and the trip-wire lives in conformance).

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::Point3;
use fs_mesh::delaunay;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 41,
                kernel_id: 77,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

/// Deterministic point cloud (splitmix), pre-generated so generation
/// cost stays out of the measurement.
fn cloud(n: usize, seed: u64) -> Vec<Point3> {
    let unit = |k: u64| -> f64 {
        let mut z = seed ^ 0x9e37_79b9_7f4a_7c15u64.wrapping_mul(k + 1);
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
        z ^= z >> 31;
        (z >> 11) as f64 / (1u64 << 53) as f64
    };
    (0..n)
        .map(|i| {
            let k = i as u64 * 3;
            Point3::new(unit(k), unit(k + 1), unit(k + 2))
        })
        .collect()
}

fn machine_fingerprint() -> String {
    let model = std::process::Command::new("sysctl")
        .args(["-n", "hw.model"])
        .output()
        .ok()
        .and_then(|o| String::from_utf8(o.stdout).ok())
        .map_or_else(|| "unknown".to_string(), |s| s.trim().to_string());
    format!(
        "{}-{}-{model}",
        std::env::consts::OS,
        std::env::consts::ARCH
    )
}

fn run_rung(n: usize, full_insphere_audit: bool) -> (f64, usize) {
    let points = cloud(n, 0xbead_5eed);
    let start = std::time::Instant::now();
    let tetra = with_cx(|cx| delaunay(&points, cx).expect("delaunay"));
    let secs = start.elapsed().as_secs_f64();
    let audit = tetra.audit(full_insphere_audit);
    assert!(audit.clean(), "the exact audit stays clean at n = {n}");
    let tets = tetra.tets().len();
    #[allow(clippy::cast_precision_loss)]
    let rate = n as f64 / secs;
    println!(
        "{{\"metric\":\"mesh-perf\",\"n\":{n},\"secs\":{secs:.3},\
         \"points_per_s\":{rate:.0},\"tets\":{tets},\
         \"insphere_audit\":{full_insphere_audit},\
         \"build\":\"{}\",\"machine\":\"{}\"}}",
        if cfg!(debug_assertions) {
            "debug"
        } else {
            "release"
        },
        machine_fingerprint()
    );
    (rate, tets)
}

#[test]
#[ignore = "perf lane: run explicitly in release with --ignored"]
fn perf_lane_ladder() {
    // Full insphere certification anchors at 1e4 (it is exact
    // arithmetic over all tets — the 1e5 attempt ran for minutes and
    // taught the calibration); the larger rungs ride the structural
    // exact audit, whose invariants conformance already trip-wires.
    let (_anchor, _) = run_rung(10_000, true);
    let (rate5, _) = run_rung(100_000, false);
    let (rate6, _) = run_rung(1_000_000, false);
    // Throughput should not collapse with scale (BRIO locality's job):
    // allow a generous 4x degradation from 1e5 to 1e6.
    assert!(
        rate6 > rate5 / 4.0,
        "throughput holds within 4x across a decade: {rate5:.0} -> {rate6:.0} pts/s"
    );
    // The 1e7 rung on request (memory + minutes).
    if std::env::var("FS_MESH_PERF_FULL").as_deref() == Ok("1") {
        let (rate7, _) = run_rung(10_000_000, false);
        assert!(
            rate7 > rate6 / 4.0,
            "the 1e7 rung holds within 4x of 1e6: {rate7:.0} vs {rate6:.0}"
        );
    } else {
        println!(
            "{{\"metric\":\"mesh-perf\",\"n\":10000000,\
             \"status\":\"SKIPPED (set FS_MESH_PERF_FULL=1)\"}}"
        );
    }
}
