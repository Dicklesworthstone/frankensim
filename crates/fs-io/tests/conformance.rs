//! fs-io conformance (the wqd.25 bead). Acceptance: import→repair→
//! promote works on a dirty-mesh defect zoo with correct receipts;
//! unrepaired defects BLOCK promotion with actionable diagnostics;
//! export→import round trips agree within format precision; fuzzing
//! finds no panics; catalogs validate with helpful errors; container
//! exports (3MF/GLB/VTK) are structurally valid.
//! Aggregate outcomes use canonical fs-obs conformance records; assertions
//! reached before those outcomes remain ordinary Rust test diagnostics.

use fs_geom::Point3;
use fs_io::quarantine::import_mesh;
use fs_io::{ColumnKind, ColumnSpec, Schema, export_3mf, export_glb, export_vtk, promote};
use fs_rep_mesh::Soup;

const DETERMINISTIC_SEED: u64 = 0;
const IO_003_INPUT_SEED: u64 = 0x10_0003;

fn verdict(case: &str, pass: bool, detail: &str, seed: u64) {
    let mut emitter = fs_obs::Emitter::new("fs-io/conformance", case);
    let event = emitter.emit(
        if pass {
            fs_obs::Severity::Info
        } else {
            fs_obs::Severity::Error
        },
        fs_obs::EventKind::ConformanceCase {
            suite: "fs-io/conformance".to_string(),
            case: case.to_string(),
            pass,
            detail: detail.to_string(),
            seed,
        },
        None,
    );
    fs_obs::lint_failure_record(&event).expect("I/O verdict must be replayable");
    let line = event.to_jsonl();
    fs_obs::validate_line(&line).expect("I/O verdict must use the fs-obs wire schema");
    println!("{line}");
    assert!(pass, "case {case}: {detail}");
}

fn measurement(case: &str, name: &str, json: String) {
    let mut emitter = fs_obs::Emitter::new("fs-io/conformance", case);
    let line = emitter
        .emit(
            fs_obs::Severity::Info,
            fs_obs::EventKind::Custom {
                name: name.to_string(),
                json,
            },
            None,
        )
        .to_jsonl();
    fs_obs::validate_line(&line).expect("I/O measurement must use the fs-obs wire schema");
    println!("{line}");
}

fn lcg(seed: &mut u64) -> u64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *seed
}

/// A closed tetrahedron (watertight, manifold).
fn tetra() -> Soup {
    Soup {
        positions: vec![
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(0.5, 1.0, 0.0),
            Point3::new(0.5, 0.3, 1.0),
        ],
        triangles: vec![[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]],
    }
}

/// An icosahedron-ish closed mesh (more faces for round-trip checks).
fn octa() -> Soup {
    Soup {
        positions: vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(-1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, -1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(0.0, 0.0, -1.0),
        ],
        triangles: vec![
            [0, 2, 4],
            [2, 1, 4],
            [1, 3, 4],
            [3, 0, 4],
            [2, 0, 5],
            [1, 2, 5],
            [3, 1, 5],
            [0, 3, 5],
        ],
    }
}

fn assert_same_geometry(a: &Soup, b: &Soup, tol: f64, what: &str) {
    assert_eq!(
        a.triangles.len(),
        b.triangles.len(),
        "{what}: triangle count"
    );
    // Compare as centroid multisets (welding may renumber).
    let centroids = |s: &Soup| {
        let mut cs: Vec<[i64; 3]> = (0..s.triangles.len())
            .map(|t| {
                let [p, q, r] = s.tri(t);
                [
                    ((p.x + q.x + r.x) / 3.0 / tol).round() as i64,
                    ((p.y + q.y + r.y) / 3.0 / tol).round() as i64,
                    ((p.z + q.z + r.z) / 3.0 / tol).round() as i64,
                ]
            })
            .collect();
        cs.sort_unstable();
        cs
    };
    assert_eq!(centroids(a), centroids(b), "{what}: centroid multiset");
}

#[test]
fn io_001_round_trips_within_format_precision() {
    let mesh = octa();
    // STL binary: f32 precision.
    let stl = fs_io::stl::write_stl(&mesh);
    let back = fs_io::stl::read_stl(&stl).expect("stl parses");
    assert_same_geometry(&mesh, &back, 1e-6, "stl");
    // Determinism: same soup, same bytes.
    assert_eq!(stl, fs_io::stl::write_stl(&mesh), "stl bytes deterministic");
    // OBJ: f64 text round trip is exact.
    let obj = fs_io::obj::write_obj(&mesh);
    let back = fs_io::obj::read_obj(&obj).expect("obj parses");
    assert_same_geometry(&mesh, &back, 1e-14, "obj");
    assert_eq!(back.positions[0].x.to_bits(), mesh.positions[0].x.to_bits());
    // PLY ASCII: double round trip is exact.
    let ply = fs_io::ply::write_ply(&mesh);
    let back = fs_io::ply::read_ply(ply.as_bytes()).expect("ply parses");
    assert_same_geometry(&mesh, &back, 1e-14, "ply");
    // ASCII STL import (hand-authored fixture).
    let ascii = "solid t\n facet normal 0 0 1\n  outer loop\n   vertex 0 0 0\n   \
                 vertex 1 0 0\n   vertex 0 1 0\n  endloop\n endfacet\nendsolid t\n";
    let tri = fs_io::stl::read_stl(ascii.as_bytes()).expect("ascii stl");
    assert_eq!(tri.triangles.len(), 1);
    verdict(
        "io-001",
        true,
        "STL(f32)/OBJ(exact)/PLY(exact) round trips; deterministic bytes",
        DETERMINISTIC_SEED,
    );
}

#[test]
fn io_002_quarantine_repair_promote_on_the_defect_zoo() {
    // Dirty mesh: tetra + duplicate face + degenerate face + a hole
    // (remove one face) + an unreferenced vertex.
    let clean = tetra();
    let mut dirty = clean.clone();
    dirty.triangles.push(dirty.triangles[0]); // duplicate
    dirty.triangles.push([1, 1, 2]); // degenerate
    dirty.triangles.remove(2); // open a hole
    dirty.positions.push(Point3::new(9.0, 9.0, 9.0)); // unreferenced
    let obj_text = fs_io::obj::write_obj(&dirty);
    let q = import_mesh(obj_text.as_bytes(), "obj").expect("import");
    // The census names every defect class.
    let classes: Vec<&str> = q.defects.iter().map(|d| d.class).collect();
    for expect in [
        "duplicate-face",
        "degenerate-face",
        "unreferenced-vertex",
        "non-manifold-or-open",
    ] {
        assert!(
            classes.contains(&expect),
            "census must include {expect}: {classes:?}"
        );
    }
    // Promotion repairs everything and yields Evidence + a receipt.
    let (evidence, receipt) = promote(q, 16).expect("promotes after repair");
    assert_eq!(
        evidence.value.triangles.len(),
        4,
        "healed back to a closed tetra"
    );
    assert!(receipt.contains("\"trust\":\"promoted\""));
    assert!(
        receipt.contains("duplicate-face"),
        "receipt records the census"
    );
    measurement("io-002/promotion-receipt", "io-promotion-receipt", receipt);
    // Unrepairable: a hole larger than the fill budget BLOCKS promotion.
    let mut gaping = octa();
    gaping.triangles.truncate(3); // massive open boundary
    let obj2 = fs_io::obj::write_obj(&gaping);
    let q2 = import_mesh(obj2.as_bytes(), "obj").expect("import");
    let refusal = promote(q2, 0).expect_err("must refuse");
    assert!(
        refusal
            .blocking
            .iter()
            .any(|b| b.contains("non-manifold-or-open")),
        "refusal names the blocker: {:?}",
        refusal.blocking
    );
    assert!(
        refusal.fixes.iter().any(|f| f.contains("max_hole_edges")),
        "diagnostics must be actionable: {:?}",
        refusal.fixes
    );
    assert!(refusal.receipt_json.contains("\"trust\":\"refused\""));
    verdict(
        "io-002",
        true,
        "defect zoo censused, repaired, promoted with receipts; oversized hole refused \
         with actionable fixes",
        DETERMINISTIC_SEED,
    );
}

#[test]
fn io_003_fuzz_never_panics() {
    let mut seed = IO_003_INPUT_SEED;
    let stl = fs_io::stl::write_stl(&octa());
    let ply = fs_io::ply::write_ply(&octa());
    let obj = fs_io::obj::write_obj(&octa());
    let mut parsed = 0usize;
    for _ in 0..1500 {
        for base in [&stl[..], ply.as_bytes(), obj.as_bytes()] {
            let mut mutated = base.to_vec();
            for _ in 0..=(lcg(&mut seed) % 8) {
                let pos = (lcg(&mut seed) as usize) % mutated.len();
                mutated[pos] = (lcg(&mut seed) % 256) as u8;
            }
            for format in ["stl", "ply", "obj"] {
                if import_mesh(&mutated, format).is_ok() {
                    parsed += 1;
                }
            }
        }
    }
    // Truncation prefixes of every format.
    for base in [&stl[..], ply.as_bytes(), obj.as_bytes()] {
        for cut in (0..base.len()).step_by(7) {
            for format in ["stl", "ply", "obj"] {
                let _ = import_mesh(&base[..cut], format);
            }
        }
    }
    // Pure junk.
    for _ in 0..500 {
        let len = (lcg(&mut seed) % 200) as usize;
        let junk: Vec<u8> = (0..len).map(|_| (lcg(&mut seed) % 256) as u8).collect();
        for format in ["stl", "ply", "obj"] {
            let _ = import_mesh(&junk, format);
        }
    }
    measurement(
        "io-003/fuzz",
        "io-parser-fuzz",
        format!("{{\"seed\":{IO_003_INPUT_SEED},\"mutants\":13500,\"still_parse\":{parsed}}}"),
    );
    verdict(
        "io-003",
        true,
        "13.5k mutants + truncations + junk: structured results, no panics",
        IO_003_INPUT_SEED,
    );
}

#[test]
fn io_004_ply_face_lists_reject_non_integer_values() {
    let base = "ply\nformat ascii 1.0\nelement vertex 3\nproperty double x\nproperty double y\n\
                property double z\nelement face 1\nproperty list uchar uint vertex_indices\n\
                end_header\n0 0 0\n1 0 0\n0 1 0\n";
    for (tail, needle) in [
        ("-3 0 1 2\n", "list count"),
        ("3.5 0 1 2\n", "list count"),
        ("3 0 -1 2\n", "list item"),
        ("3 0 1.25 2\n", "list item"),
    ] {
        let got = fs_io::ply::read_ply(format!("{base}{tail}").as_bytes());
        assert!(
            matches!(&got, Err(fs_io::IoError::Malformed { .. })),
            "expected malformed PLY list value for {tail:?}, got {got:?}"
        );
        if let Err(fs_io::IoError::Malformed { what, .. }) = &got {
            assert!(what.contains(needle), "{tail:?} error was {what:?}");
        }
    }
    verdict(
        "io-004",
        true,
        "PLY face list counts and indices reject negative/fractional values",
        DETERMINISTIC_SEED,
    );
}

#[test]
fn io_005_catalog_schema_validation_teaches() {
    let schema = Schema {
        columns: vec![
            ColumnSpec {
                name: "section",
                kind: ColumnKind::Text,
                required: true,
            },
            ColumnSpec {
                name: "area_in2",
                kind: ColumnKind::Number { min: 0.0, max: 1e4 },
                required: true,
            },
            ColumnSpec {
                name: "ix_in4",
                kind: ColumnKind::Number { min: 0.0, max: 1e6 },
                required: true,
            },
        ],
    };
    // The AISC-flavored happy path (quoted field included).
    let csv = "section,area_in2,ix_in4\n\"W14x90\",26.5,999\nW12x65,19.1,533\n";
    let catalog = schema.parse_csv(csv).expect("valid catalog");
    assert_eq!(catalog.rows.len(), 2);
    assert!((catalog.numbers[0]["area_in2"] - 26.5).abs() < 1e-12);
    assert_eq!(catalog.rows[0]["section"], "W14x90");
    // Violations teach: row, column, offending text, expectation.
    let bad_number = schema.parse_csv("section,area_in2,ix_in4\nW1,abc,3\n");
    assert!(
        matches!(&bad_number, Err(fs_io::IoError::Schema { .. })),
        "expected a schema error, got {bad_number:?}"
    );
    if let Err(fs_io::IoError::Schema { row, column, what }) = &bad_number {
        assert_eq!((*row, column.as_str()), (1, "area_in2"));
        assert!(what.contains("abc"), "must name the offender: {what}");
    }
    let out_of_range = schema.parse_csv("section,area_in2,ix_in4\nW1,-3,3\n");
    assert!(matches!(out_of_range, Err(fs_io::IoError::Schema { .. })));
    let missing_col = schema.parse_csv("section,area_in2\nW1,3\n");
    assert!(
        matches!(&missing_col, Err(fs_io::IoError::Schema { .. })),
        "expected a schema error, got {missing_col:?}"
    );
    if let Err(fs_io::IoError::Schema { column, what, .. }) = &missing_col {
        assert_eq!(column, "ix_in4");
        assert!(what.contains("found:"), "lists what WAS found: {what}");
    }
    // JSON path, same schema.
    let json = r#"[{"section": "W14x90", "area_in2": 26.5, "ix_in4": 999}]"#;
    let jcat = schema.parse_json(json).expect("json catalog");
    assert!((jcat.numbers[0]["area_in2"] - 26.5).abs() < 1e-12);
    let json_whitespace = format!("{json} \n\r\t");
    assert_eq!(
        schema
            .parse_json(&json_whitespace)
            .expect("JSON whitespace may follow the top-level array"),
        jcat
    );
    let suffix = " \n\t";
    let trailing = format!("{json}{suffix}THIS_IS_NOT_JSON");
    let trailing_error = schema
        .parse_json(&trailing)
        .expect_err("trailing bytes must not be promoted as a catalog");
    assert!(
        matches!(
            &trailing_error,
            fs_io::IoError::Malformed { at, what }
                if *at == json.len() + suffix.len() && what.contains("trailing bytes")
        ),
        "unexpected trailing-byte diagnosis: {trailing_error:?}"
    );
    let non_json_whitespace = format!("{json}\u{000b}");
    assert!(matches!(
        schema.parse_json(&non_json_whitespace),
        Err(fs_io::IoError::Malformed { at, .. }) if at == json.len()
    ));
    // regression: a JSON string with multi-byte UTF-8 must round-trip, not be
    // split into Latin-1 chars ("café–90" was becoming "cafÃ©â\u{80}\u{93}90").
    let utf8 = r#"[{"section": "café–90", "area_in2": 26.5, "ix_in4": 999}]"#;
    let ucat = schema.parse_json(utf8).expect("utf-8 json catalog");
    assert_eq!(ucat.rows[0]["section"], "café–90");
    assert!(
        schema.parse_json("[{\"section\": []}]").is_err(),
        "nested JSON refused"
    );
    verdict(
        "io-005",
        true,
        "CSV+JSON catalogs validate; errors name row/column/offender",
        DETERMINISTIC_SEED,
    );
}

#[test]
fn io_006_container_exports_are_structurally_valid() {
    let mesh = tetra();
    // 3MF: ZIP magic, EOCD present, model XML findable, entry count 3.
    let pkg = export_3mf(&mesh);
    assert_eq!(&pkg[0..4], &[0x50, 0x4B, 0x03, 0x04], "zip local header");
    let eocd = pkg
        .windows(4)
        .rposition(|w| w == [0x50, 0x4B, 0x05, 0x06])
        .expect("EOCD present");
    let entries = u16::from_le_bytes([pkg[eocd + 10], pkg[eocd + 11]]);
    assert_eq!(entries, 3, "content-types + rels + model");
    let model_needle = b"3dmodel.model";
    assert!(
        pkg.windows(model_needle.len()).any(|w| w == model_needle),
        "model part present"
    );
    let xml_needle = b"<triangle v1=";
    assert!(pkg.windows(xml_needle.len()).any(|w| w == xml_needle));
    assert_eq!(pkg, export_3mf(&mesh), "3MF bytes deterministic");
    // GLB: header, chunk sizes consistent, JSON chunk parses shape-wise.
    let glb = export_glb(&mesh);
    assert_eq!(&glb[0..4], b"glTF");
    let total = u32::from_le_bytes([glb[8], glb[9], glb[10], glb[11]]) as usize;
    assert_eq!(total, glb.len(), "declared GLB length matches");
    let json_len = u32::from_le_bytes([glb[12], glb[13], glb[14], glb[15]]) as usize;
    let json = core::str::from_utf8(&glb[20..20 + json_len]).expect("JSON chunk is UTF-8");
    assert!(json.contains("\"POSITION\":0") && json.contains("\"version\":\"2.0\""));
    let bin_len_at = 20 + json_len;
    let bin_len =
        u32::from_le_bytes(glb[bin_len_at..bin_len_at + 4].try_into().expect("4")) as usize;
    assert_eq!(20 + json_len + 8 + bin_len, glb.len(), "chunk accounting");
    // VTK: counts + field section.
    let field: Vec<f64> = (0..mesh.positions.len()).map(|i| i as f64 * 0.5).collect();
    let vtk = export_vtk(&mesh, Some(("temperature", &field)));
    assert!(vtk.contains("POINTS 4 double"));
    assert!(vtk.contains("CELLS 4 16"));
    assert!(vtk.contains("SCALARS temperature double 1"));
    assert_eq!(
        vtk.matches('\n').count(),
        vtk.lines().count(),
        "newline-terminated records"
    );
    verdict(
        "io-006",
        true,
        "3MF zip structure, GLB chunk accounting, VTK sections all check out",
        DETERMINISTIC_SEED,
    );
}

#[test]
fn io_006b_3mf_central_directory_records_are_well_formed() {
    // Regression: each ZIP central-directory file header is a FIXED 46 bytes
    // before the name (this writer emits no extra/comment fields). A 4-byte
    // inflation (50-byte records) misaligns every following record's
    // PK\x01\x02 signature and the EOCD offsets, yielding a ZIP/OPC package no
    // conformant reader can open — invisible to a test that only checks the
    // EOCD magic and entry count. Walk the directory at the 46-byte stride and
    // require it to land exactly on the EOCD.
    let mesh = tetra();
    let pkg = export_3mf(&mesh);
    let eocd = pkg
        .windows(4)
        .rposition(|w| w == [0x50, 0x4B, 0x05, 0x06])
        .expect("EOCD present");
    let entries = u16::from_le_bytes([pkg[eocd + 10], pkg[eocd + 11]]) as usize;
    let dir_start = u32::from_le_bytes(pkg[eocd + 16..eocd + 20].try_into().expect("4")) as usize;
    let mut pos = dir_start;
    for record in 0..entries {
        assert_eq!(
            &pkg[pos..pos + 4],
            &[0x50, 0x4B, 0x01, 0x02],
            "central-directory record {record} must begin at the 46-byte-strided offset"
        );
        let name_len = u16::from_le_bytes([pkg[pos + 28], pkg[pos + 29]]) as usize;
        pos += 46 + name_len;
    }
    assert_eq!(
        pos, eocd,
        "the central directory must end exactly at the EOCD (records are 46+name bytes)"
    );
}

/// wqd.25.1 — PLY element order is the HEADER's to define: face-first
/// files import identically to vertex-first files (ASCII and binary),
/// out-of-range indices still refuse after the full parse with the
/// exact triangle ordinal, and the list/triangle resource bounds hold
/// regardless of order.
#[test]
fn ply_face_first_element_order_is_legal() {
    let vertex_first = "ply\nformat ascii 1.0\nelement vertex 3\n\
        property double x\nproperty double y\nproperty double z\n\
        element face 1\nproperty list uchar uint vertex_indices\nend_header\n\
        0 0 0\n1 0 0\n0 1 0\n3 0 1 2\n";
    let face_first = "ply\nformat ascii 1.0\n\
        element face 1\nproperty list uchar uint vertex_indices\n\
        element vertex 3\nproperty double x\nproperty double y\nproperty double z\n\
        end_header\n3 0 1 2\n0 0 0\n1 0 0\n0 1 0\n";
    let a = fs_io::ply::read_ply(vertex_first.as_bytes()).expect("vertex-first parses");
    let b = fs_io::ply::read_ply(face_first.as_bytes()).expect("face-first parses");
    assert_eq!(a.positions, b.positions);
    assert_eq!(a.triangles, b.triangles);
    // Binary face-first: same header discipline, little-endian body.
    let mut bin = Vec::new();
    bin.extend_from_slice(
        b"ply\nformat binary_little_endian 1.0\n\
          element face 1\nproperty list uchar uint vertex_indices\n\
          element vertex 3\nproperty float x\nproperty float y\nproperty float z\n\
          end_header\n",
    );
    bin.push(3u8);
    for i in [0u32, 1, 2] {
        bin.extend_from_slice(&i.to_le_bytes());
    }
    for v in [[0f32, 0., 0.], [1., 0., 0.], [0., 1., 0.]] {
        for c in v {
            bin.extend_from_slice(&c.to_le_bytes());
        }
    }
    let c = fs_io::ply::read_ply(&bin).expect("binary face-first parses");
    assert_eq!(c.triangles, a.triangles);
    // Out-of-range indices still refuse AFTER the full parse, naming
    // the offending triangle ordinal.
    let bad = face_first.replace("3 0 1 2", "3 0 1 9");
    match fs_io::ply::read_ply(bad.as_bytes()) {
        Err(fs_io::IoError::Malformed { at, what }) => {
            assert_eq!(at, 0, "first triangle is the offender");
            assert!(
                what.contains("index 9") && what.contains("3 vertices"),
                "{what}"
            );
        }
        other => panic!("expected out-of-range refusal, got {other:?}"),
    }
    // Resource bounds hold in face-first order too.
    let oversized = face_first.replace("3 0 1 2", "1025 0 1 2");
    assert!(matches!(
        fs_io::ply::read_ply(oversized.as_bytes()),
        Err(fs_io::IoError::ResourceBound { .. })
    ));
}
