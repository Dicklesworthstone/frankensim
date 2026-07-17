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

const SESSION: &str = "fs-mesh/perf-lane";
const CLOUD_INPUT_SEED: u64 = 0xbead_5eed;
const CX_EXECUTION_SEED: u64 = 41;
const CX_KERNEL_ID: u64 = 77;

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: CX_EXECUTION_SEED,
                kernel_id: CX_KERNEL_ID,
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

fn emit_custom(identity: &str, severity: fs_obs::Severity, name: &str, json: String) {
    let mut emitter = fs_obs::Emitter::new(SESSION, identity);
    let event = emitter.emit(
        severity,
        fs_obs::EventKind::Custom {
            name: name.to_string(),
            json,
        },
        None,
    );
    fs_obs::lint_failure_record(&event)
        .expect("mesh perf-lane row must satisfy the failure-record lint");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("mesh perf-lane row must use the fs-obs wire schema");
    println!("{line}");
}

fn finite_json(value: f64, precision: usize) -> String {
    if value.is_finite() {
        format!("{value:.precision$}")
    } else {
        "null".to_string()
    }
}

fn json_string(value: &str) -> String {
    use core::fmt::Write as _;

    let mut out = String::with_capacity(value.len() + 2);
    out.push('"');
    for ch in value.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\u{0008}' => out.push_str("\\b"),
            '\u{000c}' => out.push_str("\\f"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            ch if ch <= '\u{001f}' => {
                let _ = write!(out, "\\u{:04x}", u32::from(ch));
            }
            ch => out.push(ch),
        }
    }
    out.push('"');
    out
}

fn run_rung(n: usize, full_insphere_audit: bool) -> (f64, usize) {
    let points = cloud(n, CLOUD_INPUT_SEED);
    let start = std::time::Instant::now();
    let tetra = with_cx(|cx| delaunay(&points, cx).expect("delaunay"));
    let secs = start.elapsed().as_secs_f64();
    let audit = tetra.audit(full_insphere_audit);
    assert!(audit.clean(), "the exact audit stays clean at n = {n}");
    let tets = tetra.tets().len();
    #[allow(clippy::cast_precision_loss)]
    let rate = n as f64 / secs;
    let build = if cfg!(debug_assertions) {
        "debug"
    } else {
        "release"
    };
    let machine = machine_fingerprint();
    let machine_id = fs_obs::fnv1a64(machine.as_bytes());
    emit_custom(
        &format!("mesh-perf/n-{n}/measurement"),
        fs_obs::Severity::Info,
        "mesh-perf-rung",
        format!(
            "{{\"metric\":\"mesh-perf\",\"n\":{n},\"secs\":{},\
             \"points_per_s\":{},\"tets\":{tets},\
             \"insphere_audit\":{full_insphere_audit},\"build\":{build_json},\
             \"machine\":{machine_json},\"machine_fingerprint\":{machine_id},\
             \"machine_fingerprint_domain\":\"fnv1a64(os-arch-hw.model-label)\",\
             \"input_seed\":{CLOUD_INPUT_SEED},\
             \"cx_execution_seed\":{CX_EXECUTION_SEED},\
             \"cx_kernel_id\":{CX_KERNEL_ID},\"cx_tile\":0,\"cx_iteration\":0,\
             \"timing_replayable\":false,\"claim_scope\":\"this-run-this-machine\"}}",
            finite_json(secs, 3),
            finite_json(rate, 0),
            build_json = json_string(build),
            machine_json = json_string(&machine),
        ),
    );
    (rate, tets)
}

fn emit_full_rung_skip() {
    emit_custom(
        "mesh-perf/n-10000000/skip",
        fs_obs::Severity::Warn,
        "mesh-perf-skip",
        format!(
            "{{\"metric\":\"mesh-perf\",\"n\":10000000,\"status\":\"skipped\",\
             \"display\":\"SKIPPED (set FS_MESH_PERF_FULL=1)\",\
             \"required_env\":\"FS_MESH_PERF_FULL=1\",\"input_seed\":null,\
             \"configured_input_seed\":{CLOUD_INPUT_SEED},\"execution_stream\":null,\
             \"timing_replayable\":false}}"
        ),
    );
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
        emit_full_rung_skip();
    }
}
