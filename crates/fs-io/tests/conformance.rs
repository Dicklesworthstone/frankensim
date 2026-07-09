//! fs-io conformance (the wqd.25 bead). Acceptance: import→repair→
//! promote works on a dirty-mesh defect zoo with correct receipts;
//! unrepaired defects BLOCK promotion with actionable diagnostics;
//! export→import round trips agree within format precision; fuzzing
//! finds no panics; catalogs validate with helpful errors; container
//! exports (3MF/GLB/VTK) are structurally valid.

use fs_geom::Point3;
use fs_io::quarantine::import_mesh;
use fs_io::{ColumnKind, ColumnSpec, Schema, export_3mf, export_glb, export_vtk, promote};
use fs_rep_mesh::Soup;

fn verdict(case: &str, detail: &str) {
    println!(
        "{{\"suite\":\"fs-io/conformance\",\"case\":\"{case}\",\"verdict\":\"pass\",\
         \"detail\":\"{detail}\"}}"
    );
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
        "STL(f32)/OBJ(exact)/PLY(exact) round trips; deterministic bytes",
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
    println!("{receipt}");
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
        "defect zoo censused, repaired, promoted with receipts; oversized hole refused \
         with actionable fixes",
    );
}

#[test]
fn io_003_fuzz_never_panics() {
    let mut seed = 0x10_0003u64;
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
    println!(
        "{{\"suite\":\"fs-io/conformance\",\"metric\":\"fuzz\",\"mutants\":13500,\
         \"still_parse\":{parsed}}}"
    );
    verdict(
        "io-003",
        "13.5k mutants + truncations + junk: structured results, no panics",
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
        match fs_io::ply::read_ply(format!("{base}{tail}").as_bytes()) {
            Err(fs_io::IoError::Malformed { what, .. }) => {
                assert!(what.contains(needle), "{tail:?} error was {what:?}");
            }
            other => panic!("expected malformed PLY list value for {tail:?}, got {other:?}"),
        }
    }
    verdict(
        "io-004",
        "PLY face list counts and indices reject negative/fractional values",
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
    match bad_number {
        Err(fs_io::IoError::Schema { row, column, what }) => {
            assert_eq!((row, column.as_str()), (1, "area_in2"));
            assert!(what.contains("abc"), "must name the offender: {what}");
        }
        other => panic!("expected a schema error, got {other:?}"),
    }
    let out_of_range = schema.parse_csv("section,area_in2,ix_in4\nW1,-3,3\n");
    assert!(matches!(out_of_range, Err(fs_io::IoError::Schema { .. })));
    let missing_col = schema.parse_csv("section,area_in2\nW1,3\n");
    match missing_col {
        Err(fs_io::IoError::Schema { column, what, .. }) => {
            assert_eq!(column, "ix_in4");
            assert!(what.contains("found:"), "lists what WAS found: {what}");
        }
        other => panic!("expected a schema error, got {other:?}"),
    }
    // JSON path, same schema.
    let json = r#"[{"section": "W14x90", "area_in2": 26.5, "ix_in4": 999}]"#;
    let jcat = schema.parse_json(json).expect("json catalog");
    assert!((jcat.numbers[0]["area_in2"] - 26.5).abs() < 1e-12);
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
        "CSV+JSON catalogs validate; errors name row/column/offender",
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
        "3MF zip structure, GLB chunk accounting, VTK sections all check out",
    );
}
