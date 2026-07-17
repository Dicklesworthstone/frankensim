//! CCD-locality A/B harness (bead fz2.2): measured L3 topology, pin
//! verification, first-touch page-to-node audit, and the L3-island
//! effect quantified by timing proxy — spread (one worker per CCD,
//! working set inside its own 32 MiB island) vs packed (every worker
//! contending for ONE island) vs unpinned (OS placement).
//!
//! Report-only rows by design: locality changes TIMING, never bits
//! (P2), so the perf gates live in the roofline lanes; this harness
//! proves the MECHANISMS work and ledgers their measured effect.
//! Run: `cargo test -p fs-substrate --release --test ccd_ab -- --ignored --nocapture`

use fs_substrate::affinity::{CcdTopology, measured_l3_groups};
use fs_substrate::os_affinity::{OsAffinityError, current_cpu, page_nodes, pin_current_thread};
use std::time::Instant;

const SUITE: &str = "fs-substrate/ccd-ab";

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
    fs_obs::lint_failure_record(&event).expect("CCD A/B row must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("CCD A/B row must use the fs-obs wire schema");
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

fn option_json(value: Option<&str>) -> String {
    value.map_or_else(|| "null".to_string(), json_string)
}

fn node_counts_json(per_node: &std::collections::BTreeMap<i32, usize>) -> String {
    use core::fmt::Write as _;

    let mut out = String::from("{");
    for (index, (node, pages)) in per_node.iter().enumerate() {
        if index > 0 {
            out.push(',');
        }
        let _ = write!(out, "\"{node}\":{pages}");
    }
    out.push('}');
    out
}

/// Stream-sum over a buffer; returns (checksum, GB/s aggregate for one pass).
fn stream_pass(buf: &[u64]) -> u64 {
    let mut acc = 0u64;
    for &v in buf {
        acc = acc.wrapping_add(v);
    }
    acc
}

/// Aggregate GB/s for `threads` workers each streaming its own buffer
/// `reps` times, with a per-worker pin list (empty = unpinned).
fn measure(buffers: &mut [Vec<u64>], pins: &[Vec<u32>], reps: usize) -> f64 {
    let bytes_total = buffers.iter().map(|b| b.len() * 8 * reps).sum::<usize>();
    let t0 = Instant::now();
    std::thread::scope(|s| {
        for (buf, pin) in buffers.iter().zip(pins) {
            let pin = pin.clone();
            s.spawn(move || {
                if !pin.is_empty() {
                    let _ = pin_current_thread(&pin);
                }
                let mut acc = 0u64;
                for _ in 0..reps {
                    acc = acc.wrapping_add(stream_pass(buf));
                }
                std::hint::black_box(acc);
            });
        }
    });
    bytes_total as f64 / t0.elapsed().as_secs_f64() / 1e9
}

#[test]
#[ignore = "perf harness: run explicitly in release with --ignored"]
fn ccd_locality_ab() {
    // --- Topology: measured where possible, heuristic otherwise. ---
    let groups = measured_l3_groups();
    let probe = fs_substrate::CapabilityProbe::topology_only();
    let machine = probe.fingerprint();
    let (topo, source) = match CcdTopology::from_l3_groups(&groups) {
        Some(t) => (t, "measured-sysfs"),
        None => (CcdTopology::from_probe(&probe), "heuristic-probe"),
    };
    measurement(
        "ccd-topology/measurement",
        "ccd-topology",
        format!(
            "{{\"metric\":\"ccd-topology\",\"source\":{},\"groups\":{},\
             \"cores_per_group\":{},\"machine\":{machine},\"input_seed\":null}}",
            json_string(source),
            topo.ccds,
            topo.cores_per_ccd,
        ),
    );

    // --- Pin verification (structural: refusal or proof, never a no-op). ---
    match current_cpu() {
        Ok(_) => {
            let target = groups.first().map_or(0, |g| g[0]);
            pin_current_thread(&[target]).expect("pin to first core of group 0");
            let now = current_cpu().expect("getcpu");
            assert_eq!(now, target, "pinned thread must run on its target CPU");
            measurement(
                "pin-verify/measurement",
                "pin-verify",
                format!(
                    "{{\"metric\":\"pin-verify\",\"verdict\":\"pass\",\"cpu\":{now},\
                     \"machine\":{machine},\"input_seed\":null,\"os_schedule_replay\":false}}"
                ),
            );
        }
        Err(OsAffinityError::Unsupported(why)) => {
            measurement(
                "pin-verify/measurement",
                "pin-verify",
                format!(
                    "{{\"metric\":\"pin-verify\",\"verdict\":\"skip\",\"why\":{},\
                     \"machine\":{machine},\"input_seed\":null,\"os_schedule_replay\":false}}",
                    json_string(why),
                ),
            );
        }
        Err(e) => panic!("getcpu failed structurally: {e}"),
    }

    // --- First-touch page-to-node audit (documents NPS mode as configured). ---
    let mut audit_buf = vec![0u8; 64 << 20];
    for (i, b) in audit_buf.iter_mut().enumerate().step_by(4096) {
        *b = (i % 251) as u8; // fault every page in
    }
    // An unaligned base can leave the FINAL partial OS page untouched
    // by the strided loop (measured on ts2: one -ENOENT page) — fault
    // it explicitly.
    let len = audit_buf.len();
    audit_buf[len - 1] = 1;
    match page_nodes(&audit_buf, 4096) {
        Ok(nodes) => {
            let mut per_node = std::collections::BTreeMap::<i32, usize>::new();
            for &n in &nodes {
                *per_node.entry(n).or_default() += 1;
            }
            measurement(
                "first-touch-audit/measurement",
                "first-touch-audit",
                format!(
                    "{{\"metric\":\"first-touch-audit\",\"pages\":{},\"per_node\":{},\
                     \"machine\":{machine},\"input_seed\":null,\"os_placement_replay\":false}}",
                    nodes.len(),
                    node_counts_json(&per_node),
                ),
            );
            assert!(
                per_node.keys().all(|&n| n >= 0),
                "every touched page reports a real node"
            );
        }
        Err(OsAffinityError::Unsupported(why)) => {
            measurement(
                "first-touch-audit/measurement",
                "first-touch-audit",
                format!(
                    "{{\"metric\":\"first-touch-audit\",\"verdict\":\"skip\",\"why\":{},\
                     \"machine\":{machine},\"input_seed\":null,\"os_placement_replay\":false}}",
                    json_string(why),
                ),
            );
        }
        Err(e) => panic!("move_pages query failed structurally: {e}"),
    }
    drop(audit_buf);

    // --- L3-island A/B (needs >= 2 measured groups + working pinning). ---
    if groups.len() < 2 || current_cpu().is_err() {
        measurement(
            "l3-island-ab/measurement",
            "l3-island-ab",
            format!(
                "{{\"metric\":\"l3-island-ab\",\"verdict\":\"skip\",\
                 \"why\":\"needs >=2 measured L3 groups + pinning\",\"machine\":{machine},\
                 \"input_seed\":null,\"timing_seed\":null,\"timing_replay\":false}}"
            ),
        );
        return;
    }
    // One worker per group; per-worker working set ~24 MiB (inside one
    // 32 MiB island), streamed repeatedly.
    let g = groups.len();
    let words = (24 << 20) / 8;
    let mut buffers: Vec<Vec<u64>> = (0..g)
        .map(|k| (0..words).map(|i| (i as u64) ^ (k as u64)).collect())
        .collect();
    let reps = 8;
    // SPREAD: worker k owns group k's cores.
    let spread_pins: Vec<Vec<u32>> = groups.clone();
    // PACKED: every worker squeezed onto group 0's cores.
    let packed_pins: Vec<Vec<u32>> = (0..g).map(|_| groups[0].clone()).collect();
    // UNPINNED: OS placement.
    let free_pins: Vec<Vec<u32>> = (0..g).map(|_| Vec::new()).collect();
    // Warm + best-of-3 per configuration.
    let best = |buffers: &mut [Vec<u64>], pins: &[Vec<u32>]| -> f64 {
        let mut best = 0.0f64;
        for _ in 0..3 {
            best = best.max(measure(buffers, pins, reps));
        }
        best
    };
    let spread = best(&mut buffers, &spread_pins);
    let packed = best(&mut buffers, &packed_pins);
    let free = best(&mut buffers, &free_pins);
    let spread_over_packed = spread / packed.max(1e-9);
    let spread_over_unpinned = spread / free.max(1e-9);
    measurement(
        "l3-island-ab/measurement",
        "l3-island-ab",
        format!(
            "{{\"metric\":\"l3-island-ab\",\"groups\":{g},\"ws_mib_per_worker\":24,\
             \"spread_gbs\":{},\"packed_gbs\":{},\"unpinned_gbs\":{},\
             \"spread_over_packed\":{},\"spread_over_unpinned\":{},\
             \"machine\":{machine},\"input_seed\":null,\"trials\":3,\
             \"timing_seed\":null,\"timing_replay\":false}}",
            finite_json(spread, 1),
            finite_json(packed, 1),
            finite_json(free, 1),
            finite_json(spread_over_packed, 2),
            finite_json(spread_over_unpinned, 2),
        ),
    );
    assert!(
        spread > packed,
        "one working set per L3 island must beat {g} sets contending for one island \
         (spread {spread:.1} vs packed {packed:.1} GB/s)"
    );
}

/// Hugepage A/B (fz2.2): TLB-stressing strided walk over 1 GiB with
/// and without MADV_HUGEPAGE. Only meaningful where the THP mode is
/// `madvise` (ts2's configuration); under `always` the baseline is
/// already huge-backed and the row documents that instead. Advice is
/// applied to the freshly zero-mapped allocation BEFORE first write,
/// so the dirtying faults populate huge pages.
#[test]
#[ignore = "perf harness: run explicitly in release with --ignored"]
fn hugepage_ab() {
    use fs_substrate::affinity::thp_mode;
    use fs_substrate::os_affinity::advise_hugepages_words;

    let mode = thp_mode();
    let machine = fs_substrate::CapabilityProbe::topology_only().fingerprint();
    measurement(
        "thp-mode/measurement",
        "thp-mode",
        format!(
            "{{\"metric\":\"thp-mode\",\"mode\":{},\"machine\":{machine},\
             \"input_seed\":null}}",
            option_json(mode.as_deref()),
        ),
    );
    let words = (1usize << 30) / 8;
    let walk = |buf: &[u64]| -> f64 {
        // One touch per >4 KiB stride, wrapping: TLB misses dominate
        // on 4 KiB pages, shrink under 2 MiB pages.
        let n = buf.len();
        let mut best = f64::INFINITY;
        for _ in 0..3 {
            let t0 = Instant::now();
            let mut acc = 0u64;
            let mut i = 0usize;
            for _ in 0..n / 8 {
                acc = acc.wrapping_add(buf[i]);
                i = (i + 4099) % n;
            }
            std::hint::black_box(acc);
            best = best.min(t0.elapsed().as_secs_f64());
        }
        (n / 8 * 8) as f64 / best / 1e9
    };
    let fill = |buf: &mut [u64]| {
        for (i, v) in buf.iter_mut().enumerate() {
            *v = i as u64;
        }
    };
    // Baseline: zero-mapped, dirtied without advice, walked.
    let mut base = vec![0u64; words];
    fill(&mut base);
    let base_gbs = walk(&base);
    drop(base);
    // Advised: same shape, MADV_HUGEPAGE between map and dirty.
    let mut adv = vec![0u64; words];
    match advise_hugepages_words(&adv, 4096) {
        Ok(bytes) => {
            fill(&mut adv);
            let adv_gbs = walk(&adv);
            let speedup = adv_gbs / base_gbs.max(1e-9);
            measurement(
                "hugepage-ab/measurement",
                "hugepage-ab",
                format!(
                    "{{\"metric\":\"hugepage-ab\",\"advised_bytes\":{bytes},\
                     \"plain_gbs\":{},\"advised_gbs\":{},\"speedup\":{},\
                     \"machine\":{machine},\"input_seed\":null,\"trials\":3,\
                     \"timing_seed\":null,\"timing_replay\":false}}",
                    finite_json(base_gbs, 2),
                    finite_json(adv_gbs, 2),
                    finite_json(speedup, 3),
                ),
            );
        }
        Err(e) => {
            measurement(
                "hugepage-ab/measurement",
                "hugepage-ab",
                format!(
                    "{{\"metric\":\"hugepage-ab\",\"verdict\":\"skip\",\"why\":{},\
                     \"plain_gbs\":{},\"machine\":{machine},\"input_seed\":null,\
                     \"trials\":3,\"timing_seed\":null,\"timing_replay\":false}}",
                    json_string(&e.to_string()),
                    finite_json(base_gbs, 2),
                ),
            );
        }
    }
}
