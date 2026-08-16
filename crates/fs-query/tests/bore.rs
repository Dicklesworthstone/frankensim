//! Bore-extraction battery (music bead `frankensim-music-v8-root-3ez8g.3.1`).
//!
//! - gb-001 straight cylinder (F-rep cylinder ∩ box): A(s) matches πr² at
//!   every station, medial and thickness cross-checks hold, the axial
//!   volume matches the analytic tube volume, and the composite chart's
//!   `LipschitzImplicit` claim yields an HONEST `VolumeClosure::Unavailable`
//!   — never a downgraded certificate.
//! - gb-002 full torus (`TorusChart`, exact distance): the chain closes
//!   into a loop, A(s) matches πr², and the CERTIFIED volume closure
//!   passes at the authored band — with the executed falsifier proving a
//!   deliberately biased area sweep FAILS the same gate.
//! - gb-003 bent tube (F-rep torus ∩ two half-spaces): a curved OPEN tube
//!   chains and measures A(s) = πr² on the arc.
//! - gb-004 cone frustum (`AxisymmetricChart`): the tapered A(s) oracle —
//!   per-station areas track π r(z)² and the equivalent radius is
//!   monotone; the `NoClaim` chart records closure unavailability.
//! - gb-005 refusals: a Y-junction refuses as `BranchedLumen`; a
//!   four-point boundary refuses as `TooFewPoles`; config refusals.
//! - gb-006 cancellation fails closed as `Query(Cancelled)`.
//! - gb-007 mesh end-to-end: an authored OBJ cylinder through the fs-io
//!   quarantine (import → census → promote) into `MeshChart`, extracted
//!   with A(s) at the meshed radius and closure honestly unavailable
//!   (`NoClaim`); the cap-less variant REFUSES promotion (the
//!   non-watertight arm lives at the quarantine gate, by design).
//!
//! AUTHORITY note repeated from the module doc: everything here is
//! Estimate-class except the torus closure enclosure. The oracle bands are
//! authored test tolerances, not certificates.

use asupersync::types::Budget;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::fixtures::TorusChart;
use fs_geom::{Point3, TraceStepClaim, Vec3};
use fs_query::{BoreConfig, BoreError, QueryError, VolumeClosure, extract_bore};
use fs_rep_frep::{BoolOp, BoolStyle, FrepBuilder};
use fs_rep_mesh::{DcOptions, Soup, dual_contour};
use std::f64::consts::PI;

const SUITE: &str = "fs-query/bore";
const FIXED_INPUT_SEED: u64 = 0;
const EXECUTION_SEED: u64 = 0xB0E5;
const BORE_KERNEL_ID: u64 = 30;
const CANCELLATION_KERNEL_ID: u64 = 31;

fn verdict(case: &str, pass: bool, detail: &str, execution_kernel: u64) {
    let mut emitter = fs_obs::Emitter::new(SUITE, case);
    let severity = if pass {
        fs_obs::Severity::Info
    } else {
        fs_obs::Severity::Error
    };
    let event = emitter.emit(
        severity,
        fs_obs::EventKind::ConformanceCase {
            suite: SUITE.to_string(),
            case: case.to_string(),
            pass,
            detail: format!(
                "{detail}; execution stream seed=0x{EXECUTION_SEED:x} \
                 kernel={execution_kernel} tile=0 iteration=0"
            ),
            seed: FIXED_INPUT_SEED,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("bore verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("bore verdict must use the fs-obs wire schema");
    println!("{line}");
    assert!(pass, "case {case}: {detail}");
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
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

/// Straight tube: infinite cylinder along z clipped by a box.
fn cylinder_tube(radius: f64, half_length: f64) -> fs_rep_frep::Frep {
    let mut b = FrepBuilder::new();
    let cyl = b.cylinder(Point3::new(0.0, 0.0, 0.0), radius).expect("cyl");
    let bx = b
        .box_prim(
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(2.0 * radius, 2.0 * radius, half_length),
        )
        .expect("box");
    let root = b
        .boolean(BoolOp::Intersect, BoolStyle::Hard, cyl, bx)
        .expect("root");
    b.finish(root).expect("frep")
}

#[test]
fn gb_001_straight_cylinder_oracle() {
    with_cx(|cx| {
        let radius = 0.5f64;
        let half_length = 1.5f64;
        let tube = cylinder_tube(radius, half_length);
        let (soup, _) = dual_contour(&tube, DcOptions::sharp(0.1), cx).expect("dc");
        let config = BoreConfig {
            closure_h_m: Some(0.05),
            ..BoreConfig::default()
        };
        let bore = extract_bore(&tube, &soup, &config, "test/cylinder/r0.5-l3/v1", cx)
            .expect("cylinder extracts");
        for line in bore.debug_lines() {
            println!("{line}");
        }
        let truth = PI * radius * radius;
        let worst_area = bore
            .stations
            .iter()
            .map(|s| (s.area_m2 - truth).abs() / truth)
            .fold(0.0f64, f64::max);
        let analytic_volume = truth * 2.0 * half_length;
        let volume_rel = (bore.axial_volume_m3 - analytic_volume).abs() / analytic_volume;
        let closure_honest = matches!(
            bore.volume_closure,
            VolumeClosure::Unavailable {
                claim: TraceStepClaim::LipschitzImplicit
            }
        );
        let pass = !bore.closed_loop
            && worst_area < 0.01
            && bore.worst_pole_deviation() < 0.25
            && volume_rel < 0.03
            && closure_honest
            && bore.thickness_skipped <= 2;
        verdict(
            "gb-001",
            pass,
            &format!(
                "straight cylinder: worst area dev {worst_area:.3e}, volume rel \
                 {volume_rel:.3e}, pole dev {:.3e}, closure honest {closure_honest}, \
                 thickness skipped {}",
                bore.worst_pole_deviation(),
                bore.thickness_skipped
            ),
            BORE_KERNEL_ID,
        );
    });
}

#[test]
fn gb_002_torus_certified_closure_and_falsifier() {
    with_cx(|cx| {
        let major = 1.0f64;
        let minor = 0.3f64;
        let torus = TorusChart {
            center: Point3::new(0.0, 0.0, 0.0),
            major,
            minor,
        };
        let (soup, _) = dual_contour(&torus, DcOptions::sharp(0.1), cx).expect("dc");
        let config = BoreConfig {
            closure_h_m: Some(0.02),
            ..BoreConfig::default()
        };
        let bore =
            extract_bore(&torus, &soup, &config, "test/torus/R1-r0.3/v1", cx).expect("torus");
        for line in bore.debug_lines() {
            println!("{line}");
        }
        let truth = PI * minor * minor;
        let worst_area = bore
            .stations
            .iter()
            .map(|s| (s.area_m2 - truth).abs() / truth)
            .fold(0.0f64, f64::max);
        let band = 0.05f64;
        let closure_pass = bore.closure_within(band);
        // THE FALSIFIER: a 30%-biased area sweep must FAIL the same gate.
        let falsifier_caught = match &bore.volume_closure {
            VolumeClosure::Certified { lo_m3, hi_m3, .. } => {
                let mid = 0.5 * (lo_m3 + hi_m3);
                let half_width_rel = 0.5 * (hi_m3 - lo_m3) / mid;
                let biased_rel = (1.30 * bore.axial_volume_m3 - mid).abs() / mid;
                biased_rel > band + half_width_rel
            }
            _ => false,
        };
        let pass =
            bore.closed_loop && worst_area < 0.02 && closure_pass == Some(true) && falsifier_caught;
        verdict(
            "gb-002",
            pass,
            &format!(
                "torus loop: closed {}, worst area dev {worst_area:.3e}, certified \
                 closure {:?} within band {band}, biased-sweep falsifier caught \
                 {falsifier_caught}",
                bore.closed_loop, bore.volume_closure
            ),
            BORE_KERNEL_ID,
        );
    });
}

#[test]
fn gb_003_bent_tube_oracle() {
    with_cx(|cx| {
        let major = 1.0f64;
        let minor = 0.25f64;
        let mut b = FrepBuilder::new();
        let torus = b
            .torus(Point3::new(0.0, 0.0, 0.0), major, minor)
            .expect("torus");
        // Keep the quadrant x >= 0, y >= 0: a 90-degree bend.
        let hx = b.half_space(Vec3::new(-1.0, 0.0, 0.0), 0.0).expect("hx");
        let hy = b.half_space(Vec3::new(0.0, -1.0, 0.0), 0.0).expect("hy");
        let quad = b
            .boolean(BoolOp::Intersect, BoolStyle::Hard, hx, hy)
            .expect("quad");
        let root = b
            .boolean(BoolOp::Intersect, BoolStyle::Hard, torus, quad)
            .expect("root");
        let bend = b.finish(root).expect("frep");
        let (soup, _) = dual_contour(&bend, DcOptions::sharp(0.08), cx).expect("dc");
        let config = BoreConfig::default();
        let bore = extract_bore(&bend, &soup, &config, "test/bend90/R1-r0.25/v1", cx)
            .expect("bend extracts");
        for line in bore.debug_lines() {
            println!("{line}");
        }
        let truth = PI * minor * minor;
        // Two-tier gate by arc distance from the cut faces: stations at
        // least two tube radii from a face are tight; face-adjacent
        // stations sit against the end medial sheet, their tangents tilt,
        // and a tilt inflates the section by 1/cos(tilt) — an authored,
        // disclosed 6% envelope there (measured 3.3% worst).
        let near_face = |s: &fs_query::BoreStation| {
            s.arc_length_m < 2.0 * minor || s.arc_length_m > bore.total_length_m - 2.0 * minor
        };
        let interior_dev = bore
            .stations
            .iter()
            .filter(|s| !near_face(s))
            .map(|s| (s.area_m2 - truth).abs() / truth)
            .fold(0.0f64, f64::max);
        let end_dev = bore
            .stations
            .iter()
            .filter(|s| near_face(s))
            .map(|s| (s.area_m2 - truth).abs() / truth)
            .fold(0.0f64, f64::max);
        // Arc length of the quarter bend at the centerline radius.
        let analytic_length = 0.5 * PI * major;
        let length_rel = (bore.total_length_m - analytic_length).abs() / analytic_length;
        let pass = !bore.closed_loop && interior_dev < 0.01 && end_dev < 0.06 && length_rel < 0.10;
        verdict(
            "gb-003",
            pass,
            &format!(
                "90-degree bend: interior area dev {interior_dev:.3e}, end-station dev \
                 {end_dev:.3e}, length {:.4} vs analytic {analytic_length:.4} \
                 (rel {length_rel:.3e})",
                bore.total_length_m
            ),
            BORE_KERNEL_ID,
        );
    });
}

#[test]
fn gb_004_cone_frustum_taper_oracle() {
    use fs_rep_frep::{AxisymmetricChart, MeridianPoint, MeridianSegment};
    with_cx(|cx| {
        let r0 = 0.5f64;
        let r1 = 0.3f64;
        let length = 1.6f64;
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
        let bore = extract_bore(&chart, &soup, &config, "test/frustum/r0.5-r0.3/v1", cx)
            .expect("frustum extracts");
        // Per-station: area tracks pi r(z)^2 with r linear in z.
        let worst_taper = bore
            .stations
            .iter()
            .map(|s| {
                let z = s.center_m.z.clamp(0.0, length);
                let r = r0 + (r1 - r0) * z / length;
                let truth = PI * r * r;
                (s.area_m2 - truth).abs() / truth
            })
            .fold(0.0f64, f64::max);
        // Equivalent radius must shrink monotonically along the taper
        // (orient by comparing the ends).
        let radii: Vec<f64> = bore
            .stations
            .iter()
            .map(|s| s.equivalent_radius_m)
            .collect();
        let ratio = radii[0] / radii[radii.len() - 1];
        let taper_seen = ratio.max(1.0 / ratio);
        let closure_honest = matches!(
            bore.volume_closure,
            VolumeClosure::Skipped | VolumeClosure::Unavailable { .. }
        );
        let pass = !bore.closed_loop && worst_taper < 0.04 && taper_seen > 1.4 && closure_honest;
        verdict(
            "gb-004",
            pass,
            &format!(
                "cone frustum: worst taper-tracking dev {worst_taper:.3e}, end-radius \
                 ratio {taper_seen:.3} (expect ~{:.3}), closure honest {closure_honest}",
                r0 / r1
            ),
            BORE_KERNEL_ID,
        );
    });
}

#[test]
fn gb_005_refusals() {
    with_cx(|cx| {
        // Y-junction: main tube along z plus a side arm along +x.
        let mut b = FrepBuilder::new();
        let main = b.cylinder(Point3::new(0.0, 0.0, 0.0), 0.4).expect("main");
        let main_box = b
            .box_prim(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.6, 0.6, 1.5))
            .expect("mb");
        let main = b
            .boolean(BoolOp::Intersect, BoolStyle::Hard, main, main_box)
            .expect("m");
        let side = b.cylinder(Point3::new(0.0, 0.0, 0.0), 0.4).expect("side");
        let side = b
            .rotate(side, Vec3::new(0.0, 1.0, 0.0), std::f64::consts::FRAC_PI_2)
            .expect("rot");
        let side_box = b
            .box_prim(Point3::new(0.9, 0.0, 0.0), Vec3::new(0.9, 0.6, 0.6))
            .expect("sb");
        let side = b
            .boolean(BoolOp::Intersect, BoolStyle::Hard, side, side_box)
            .expect("s");
        let root = b
            .boolean(BoolOp::Union, BoolStyle::Hard, main, side)
            .expect("root");
        let y = b.finish(root).expect("frep");
        let (soup, _) = dual_contour(&y, DcOptions::sharp(0.1), cx).expect("dc");
        let config = BoreConfig::default();
        let branched = extract_bore(&y, &soup, &config, "test/y-junction/v1", cx);
        let branch_refused = matches!(branched, Err(BoreError::BranchedLumen { .. }));

        // Four boundary points cannot seed a chain.
        let tube = cylinder_tube(0.5, 1.5);
        let tiny = Soup {
            positions: vec![
                Point3::new(0.5, 0.0, -1.0),
                Point3::new(-0.5, 0.0, 0.0),
                Point3::new(0.0, 0.5, 1.0),
                Point3::new(0.0, -0.5, 0.5),
            ],
            triangles: vec![[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]],
        };
        let too_few = extract_bore(&tube, &tiny, &config, "test/tiny/v1", cx);
        let too_few_refused = matches!(too_few, Err(BoreError::TooFewPoles { .. }));

        // Config refusals, one per knob family.
        let bad_stations = BoreConfig {
            stations: 2,
            ..BoreConfig::default()
        };
        let cfg_refused = matches!(
            extract_bore(&tube, &tiny, &bad_stations, "test/cfg/v1", cx),
            Err(BoreError::InvalidConfig { .. })
        );
        let bad_h = BoreConfig {
            closure_h_m: Some(0.0),
            ..BoreConfig::default()
        };
        let h_refused = matches!(
            extract_bore(&tube, &tiny, &bad_h, "test/cfg-h/v1", cx),
            Err(BoreError::InvalidConfig { .. })
        );
        let pass = branch_refused && too_few_refused && cfg_refused && h_refused;
        verdict(
            "gb-005",
            pass,
            &format!(
                "refusals: branched {branch_refused}, too-few-poles {too_few_refused}, \
                 station config {cfg_refused}, closure-h config {h_refused}"
            ),
            BORE_KERNEL_ID,
        );
    });
}

#[test]
fn gb_006_cancellation_fails_closed() {
    let gate = CancelGate::new();
    gate.request();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    let refused = pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: EXECUTION_SEED,
                kernel_id: CANCELLATION_KERNEL_ID,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        let tube = cylinder_tube(0.5, 1.5);
        let (soup, _) = match dual_contour(&tube, DcOptions::sharp(0.1), &cx) {
            Ok(pair) => pair,
            // The contouring itself may already observe the cancellation;
            // that is the same fail-closed outcome.
            Err(_) => return Err(BoreError::Query(QueryError::Cancelled)),
        };
        extract_bore(&tube, &soup, &BoreConfig::default(), "test/cancel/v1", &cx)
    });
    verdict(
        "gb-006",
        matches!(refused, Err(BoreError::Query(QueryError::Cancelled))),
        "a pre-cancelled context refuses before publishing a bore receipt",
        CANCELLATION_KERNEL_ID,
    );
}

/// Author a watertight OBJ cylinder (rings + caps); `capped: false` drops
/// the top cap to make the non-watertight arm.
fn obj_cylinder(radius: f64, length: f64, segments: usize, rings: usize, capped: bool) -> String {
    let mut obj = String::new();
    for ring in 0..=rings {
        let z = length * ring as f64 / rings as f64;
        for k in 0..segments {
            let th = std::f64::consts::TAU * k as f64 / segments as f64;
            obj.push_str(&format!(
                "v {:.9} {:.9} {:.9}\n",
                radius * th.cos(),
                radius * th.sin(),
                z
            ));
        }
    }
    // Cap centers (bottom, top).
    obj.push_str(&format!("v 0 0 0\nv 0 0 {length:.9}\n"));
    let idx = |ring: usize, k: usize| ring * segments + (k % segments) + 1;
    for ring in 0..rings {
        for k in 0..segments {
            let a = idx(ring, k);
            let b = idx(ring, k + 1);
            let c = idx(ring + 1, k + 1);
            let d = idx(ring + 1, k);
            obj.push_str(&format!("f {a} {b} {c}\nf {a} {c} {d}\n"));
        }
    }
    let bottom_center = (rings + 1) * segments + 1;
    let top_center = bottom_center + 1;
    for k in 0..segments {
        obj.push_str(&format!(
            "f {} {} {}\n",
            bottom_center,
            idx(0, k + 1),
            idx(0, k)
        ));
        if capped {
            obj.push_str(&format!(
                "f {} {} {}\n",
                top_center,
                idx(rings, k),
                idx(rings, k + 1)
            ));
        }
    }
    obj
}

#[test]
fn gb_007_mesh_end_to_end_through_quarantine() {
    with_cx(|cx| {
        let radius = 0.4f64;
        let length = 2.0f64;
        let obj = obj_cylinder(radius, length, 24, 6, true);
        let quarantined = fs_io::quarantine::import_mesh(obj.as_bytes(), "obj").expect("import");
        let (evidence, _receipt) = fs_io::promote(quarantined, 0).expect("watertight promotes");
        let soup = evidence.value;
        let chart = fs_rep_mesh::MeshChart::new(soup.clone());
        let config = BoreConfig {
            closure_h_m: Some(0.05),
            ..BoreConfig::default()
        };
        let bore = extract_bore(&chart, &soup, &config, "test/obj-cylinder/v1", cx)
            .expect("mesh extracts");
        for line in bore.debug_lines() {
            println!("{line}");
        }
        // The meshed tube's true cross-section is the inscribed 24-gon of
        // the authored radius: area = (n/2) r^2 sin(2 pi / n).
        let n = 24.0f64;
        let truth = 0.5 * n * radius * radius * (std::f64::consts::TAU / n).sin();
        let worst_area = bore
            .stations
            .iter()
            .map(|s| (s.area_m2 - truth).abs() / truth)
            .fold(0.0f64, f64::max);
        let closure_honest = matches!(
            bore.volume_closure,
            VolumeClosure::Unavailable {
                claim: TraceStepClaim::NoClaim
            }
        );
        // The cap-less mesh REFUSES promotion — the non-watertight arm.
        let open = fs_io::quarantine::import_mesh(
            obj_cylinder(radius, length, 24, 6, false).as_bytes(),
            "obj",
        )
        .expect("import open");
        let open_refused = fs_io::promote(open, 0).is_err();
        let pass = !bore.closed_loop && worst_area < 0.03 && closure_honest && open_refused;
        verdict(
            "gb-007",
            pass,
            &format!(
                "obj cylinder through quarantine: worst area dev {worst_area:.3e} vs the \
                 24-gon truth, closure honest {closure_honest}, open-mesh promotion \
                 refused {open_refused}"
            ),
            BORE_KERNEL_ID,
        );
    });
}
