//! Cross-environment determinism probe for the bore oracles (bead
//! `frankensim-uf7cw`). This binary reproduces the gb_003 bent-tube and
//! gb_004 cone-frustum scenarios from `tests/bore.rs` with byte-for-byte
//! identical geometry construction, seeds, kernel IDs, configuration, and
//! authored pass bands — then emits one stable machine-readable payload
//! line per case instead of relying on the test harness's assertion
//! outcomes alone.
//!
//! WHY an example and not the tests themselves: the b2can incident showed
//! golden outcomes flipping across environments while every individual
//! process still printed plausible logs. A flip cannot be diagnosed from
//! PASS/FAIL text; it needs bit-level numeric receipts (station tables,
//! pole-cloud digests, dependency fingerprints) captured identically on
//! every host, in every sibling-constellation mode, by ONE code path. The
//! tests remain the product authority; this probe is their determinism
//! witness.
//!
//! Two targets cannot share code silently, so any change to the gb_003 /
//! gb_004 scenarios in `tests/bore.rs` MUST be mirrored here (and vice
//! versa): geometry builders, tolerances, source IDs, seeds, kernel IDs,
//! and pass bands are duplicated ON PURPOSE and kept honest by this
//! comment contract.
//!
//! Gauntlet tier: G5 determinism audit (per-environment digests; the
//! comparing lane lives in xtask's oracle-determinism family).
//!
//! Output contract (stdout, strict JSONL):
//!   {"record":"header","schema":"fs-query-determinism-probe-v1", ...}
//!   {"record":"case","case":"gb-003","pass":true, ...}
//!   {"record":"case","case":"gb-004","pass":true, ...}

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_rep_frep::{
    AxisymmetricChart, BoolOp, BoolStyle, FrepBuilder, MeridianPoint, MeridianSegment,
};
use fs_rep_mesh::{DcOptions, Soup, dual_contour};
use std::f64::consts::PI;

use fs_geom::{Chart, Point3, Vec3};
use fs_query::{BoreConfig, VolumeClosure, extract_bore, medial_poles};

/// Mirrors `tests/bore.rs`: identical stream identity or the numerics are
/// not the numerics under audit.
const FIXED_INPUT_SEED: u64 = 0;
const EXECUTION_SEED: u64 = 0xB0E5;
const BORE_KERNEL_ID: u64 = 30;

/// FNV-1a 64-bit over raw bytes; matches the digest family already used
/// for `BoreExtraction::boundary_digest`.
pub(crate) fn fnv1a(bytes: &[u8]) -> u64 {
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for &byte in bytes {
        hash ^= u64::from(byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

/// Bit-exact hex encoding of one f64 — never a float-to-string rendering,
/// which can differ across libc/locale while the IEEE bits agree.
pub(crate) fn bits(value: f64) -> String {
    format!("{:016x}", value.to_bits())
}

/// Order-sensitive digest over a sequence of floats.
pub(crate) fn digest_seq<'a>(values: impl Iterator<Item = &'a f64>) -> String {
    let mut buf = Vec::new();
    for value in values {
        buf.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    format!("{:016x}", fnv1a(&buf))
}

pub(crate) fn digest_pairs(poles: &[(Point3, f64)]) -> String {
    let mut buf = Vec::new();
    for (point, radius) in poles {
        for component in [point.x, point.y, point.z] {
            buf.extend_from_slice(&component.to_bits().to_le_bytes());
        }
        buf.extend_from_slice(&radius.to_bits().to_le_bytes());
    }
    format!("{:016x}", fnv1a(&buf))
}

/// Same execution context recipe as the goldens: deterministic mode,
/// fixed input seed, single bore kernel stream.
fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = ArenaPool::new(ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: EXECUTION_SEED,
                kernel_id: BORE_KERNEL_ID,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn closure_label(closure: &VolumeClosure) -> &'static str {
    match closure {
        VolumeClosure::Certified { .. } => "Certified",
        VolumeClosure::Unavailable { .. } => "Unavailable",
        VolumeClosure::Skipped => "Skipped",
    }
}

fn emit(line: &str) {
    println!("{line}");
}

#[allow(clippy::too_many_arguments)]
fn emit_case(
    case: &str,
    pass: bool,
    extraction: &fs_query::BoreExtraction,
    poles_len: usize,
    poles_digest: &str,
    lambda: f64,
) {
    let areas = digest_seq(extraction.stations.iter().map(|s| &s.area_m2));
    let radii = digest_seq(extraction.stations.iter().map(|s| &s.equivalent_radius_m));
    let arcs = digest_seq(extraction.stations.iter().map(|s| &s.arc_length_m));
    emit(&format!(
        "{{\"record\":\"case\",\"case\":\"{case}\",\"pass\":{},\"stations\":{},\
         \"boundary_digest\":\"{:016x}\",\"areas_digest\":\"{areas}\",\
         \"radius_digest\":\"{radii}\",\"arc_digest\":\"{arcs}\",\
         \"pole_count\":{poles_len},\"poles_digest\":\"{poles_digest}\",\
         \"closure\":\"{}\",\"total_length\":\"{}\",\"axial_volume\":\"{}\",\
         \"thickness_skipped\":{},\"config_lambda\":\"{}\"}}",
        pass,
        extraction.stations.len(),
        extraction.boundary_digest,
        closure_label(&extraction.volume_closure),
        bits(extraction.total_length_m),
        bits(extraction.axial_volume_m3),
        extraction.thickness_skipped,
        bits(lambda),
    ));
}

struct ProbeOutcome {
    extraction: fs_query::BoreExtraction,
    poles_len: usize,
    poles_digest: String,
}

/// Raw `medial_poles` cloud for the SAME boundary soup the extraction
/// consumes (the bead names the pole-CLOUD digest explicitly; thinning
/// downstream of this point belongs to `extract_bore`, not the witness).
fn pole_cloud(
    chart: &dyn Chart,
    soup: &Soup,
    config: &BoreConfig,
    cx: &Cx<'_>,
) -> Result<(usize, String), String> {
    let poles = medial_poles(chart, soup, config.lambda, cx)
        .map_err(|error| format!("medial_poles refused: {error}"))?;
    Ok((poles.len(), digest_pairs(&poles)))
}

fn run_bent_tube() -> Result<ProbeOutcome, String> {
    with_cx(|cx| {
        let major = 1.0_f64;
        let minor = 0.25_f64;
        let mut builder = FrepBuilder::new();
        let torus = builder
            .torus(Point3::new(0.0, 0.0, 0.0), major, minor)
            .expect("torus");
        // Keep the quadrant x >= 0, y >= 0: a 90-degree bend (as gb_003).
        let hx = builder
            .half_space(Vec3::new(-1.0, 0.0, 0.0), 0.0)
            .expect("hx");
        let hy = builder
            .half_space(Vec3::new(0.0, -1.0, 0.0), 0.0)
            .expect("hy");
        let quad = builder
            .boolean(BoolOp::Intersect, BoolStyle::Hard, hx, hy)
            .expect("quad");
        let root = builder
            .boolean(BoolOp::Intersect, BoolStyle::Hard, torus, quad)
            .expect("root");
        let bend = builder.finish(root).expect("frep");
        let (soup, _) = dual_contour(&bend, DcOptions::sharp(0.08), cx).expect("dc");
        let config = BoreConfig::default();
        let cloud = pole_cloud(&bend, &soup, &config, cx)?;
        let extraction = extract_bore(&bend, &soup, &config, "test/bend90/R1-r0.25/v1", cx)
            .map_err(|error| format!("bend extracts: {error}"))?;
        Ok(ProbeOutcome {
            extraction,
            poles_len: cloud.0,
            poles_digest: cloud.1,
        })
    })
}

fn run_cone_frustum() -> Result<ProbeOutcome, String> {
    with_cx(|cx| {
        let r0 = 0.5_f64;
        let r1 = 0.3_f64;
        let length = 1.6_f64;
        let p = |radius: f64, axial: f64| MeridianPoint { radius, axial };
        let line = |a: MeridianPoint, b: MeridianPoint| MeridianSegment::Line { start: a, end: b };
        let chart = AxisymmetricChart::try_new(vec![
            line(p(0.0, 0.0), p(r0, 0.0)),
            line(p(r0, 0.0), p(r1, length)),
            line(p(r1, length), p(0.0, length)),
            line(p(0.0, length), p(0.0, 0.0)),
        ])
        .expect("frustum admits");
        let (soup, _) = dual_contour(&chart, DcOptions::sharp(0.08), cx).expect("dc");
        let config = BoreConfig::default();
        let cloud = pole_cloud(&chart, &soup, &config, cx)?;
        let extraction = extract_bore(&chart, &soup, &config, "test/frustum/r0.5-r0.3/v1", cx)
            .map_err(|error| format!("frustum extracts: {error}"))?;
        Ok(ProbeOutcome {
            extraction,
            poles_len: cloud.0,
            poles_digest: cloud.1,
        })
    })
}

/// Authored pass-band predicate signature shared by both probe cases.
type VerdictFn = fn(&fs_query::BoreExtraction) -> bool;

/// Authored pass bands, copied verbatim from gb_003 / gb_004 verdicts.
fn bent_tube_passes(extraction: &fs_query::BoreExtraction) -> bool {
    let truth = PI * 0.25 * 0.25;
    let near_face = |s: &fs_query::BoreStation| {
        s.arc_length_m < 2.0 * 0.25 || s.arc_length_m > extraction.total_length_m - 2.0 * 0.25
    };
    let interior_dev = extraction
        .stations
        .iter()
        .filter(|s| !near_face(s))
        .map(|s| (s.area_m2 - truth).abs() / truth)
        .fold(0.0_f64, f64::max);
    let end_dev = extraction
        .stations
        .iter()
        .filter(|s| near_face(s))
        .map(|s| (s.area_m2 - truth).abs() / truth)
        .fold(0.0_f64, f64::max);
    let analytic_length = 0.5 * PI * 1.0;
    let length_rel = (extraction.total_length_m - analytic_length).abs() / analytic_length;
    !extraction.closed_loop && interior_dev < 0.01 && end_dev < 0.06 && length_rel < 0.10
}

fn cone_frustum_passes(extraction: &fs_query::BoreExtraction) -> bool {
    let r0 = 0.5_f64;
    let r1 = 0.3_f64;
    let length = 1.6_f64;
    let worst_taper = extraction
        .stations
        .iter()
        .map(|s| {
            let z = s.center_m.z.clamp(0.0, length);
            let radius = r0 + (r1 - r0) * z / length;
            let truth = PI * radius * radius;
            (s.area_m2 - truth).abs() / truth
        })
        .fold(0.0_f64, f64::max);
    let radii: Vec<f64> = extraction
        .stations
        .iter()
        .map(|s| s.equivalent_radius_m)
        .collect();
    let ratio = radii[0] / radii[radii.len() - 1];
    let taper_seen = ratio.max(1.0 / ratio);
    let closure_honest = matches!(
        extraction.volume_closure,
        VolumeClosure::Skipped | VolumeClosure::Unavailable { .. }
    );
    !extraction.closed_loop && worst_taper < 0.04 && taper_seen > 1.4 && closure_honest
}

fn probe_error(case: &str, reason: &str) {
    emit(&format!(
        "{{\"record\":\"error\",\"case\":\"{case}\",\"reason\":{}}}",
        json_escape(reason),
    ));
}

fn json_escape(text: &str) -> String {
    let mut out = String::from("\"");
    for ch in text.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' | '\r' | '\t' => out.push(' '),
            c if (c as u32) < 0x20 => {}
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

fn main() {
    let environment = std::env::var("FS_DET_ENV").unwrap_or_else(|_| "unnamed".to_string());
    emit(&format!(
        "{{\"record\":\"header\",\"schema\":\"fs-query-determinism-probe-v1\",\
         \"env\":\"{}\",\"arch\":\"{}\",\"os\":\"{}\",\
         \"input_seed\":{FIXED_INPUT_SEED:#x},\"execution_seed\":{EXECUTION_SEED:#x},\
         \"kernel_id\":{BORE_KERNEL_ID}}}",
        json_escape(&environment),
        std::env::consts::ARCH,
        std::env::consts::OS,
    ));
    let cases: [(&str, Result<ProbeOutcome, String>, VerdictFn); 2] = [
        ("gb-003", run_bent_tube(), bent_tube_passes),
        ("gb-004", run_cone_frustum(), cone_frustum_passes),
    ];
    let mut failed = false;
    for (case, outcome, passes) in cases {
        match outcome {
            Ok(probe) => {
                let pass = passes(&probe.extraction);
                failed |= !pass;
                emit_case(
                    case,
                    pass,
                    &probe.extraction,
                    probe.poles_len,
                    &probe.poles_digest,
                    BoreConfig::default().lambda,
                );
            }
            Err(reason) => {
                failed = true;
                probe_error(case, &reason);
            }
        }
    }
    if failed {
        eprintln!("determinism_probe: at least one case failed or refused; see payload lines");
        std::process::exit(1);
    }
}
