//! Real-binary, no-mock end-to-end tests for bead `frankensim-vm286`:
//! - `frankensim mesh export <promoted-mesh-hash> <ledger.db> --format glb|3mf|vtk --out <path>`
//! - `frankensim table export <solution-hash> <ledger.db> --format csv|arrow --out <path>`

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

use fs_blake3::hash_bytes;
use fs_cli::exit;
use fs_ledger::{EdgeRole, Ledger, MAIN_BRANCH};

static COUNTER: AtomicUsize = AtomicUsize::new(0);

fn scratch_dir(tag: &str) -> PathBuf {
    let n = COUNTER.fetch_add(1, Ordering::Relaxed);
    let dir = std::env::temp_dir().join(format!("fs-cli-vm286-{tag}-{}-{n}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("scratch dir");
    dir
}

fn ref_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../../data/reference-project")
}

fn run_bin(args: &[&str], expected: u8) -> (String, String) {
    let out = Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(args)
        .output()
        .expect("spawn frankensim");
    let stdout = String::from_utf8(out.stdout).expect("UTF-8 stdout");
    let stderr = String::from_utf8(out.stderr).expect("UTF-8 stderr");
    assert_eq!(
        out.status.code(),
        Some(i32::from(expected)),
        "{args:?}\n{stdout}\n{stderr}"
    );
    (stdout, stderr)
}

fn jstr<'a>(json: &'a str, key: &str) -> &'a str {
    let needle = format!("\"{key}\":\"");
    json.split(&needle)
        .nth(1)
        .and_then(|r| r.split('"').next())
        .unwrap_or_else(|| panic!("missing `{key}` in {json}"))
}

fn ju64(json: &str, key: &str) -> u64 {
    let needle = format!("\"{key}\":");
    json.split(&needle)
        .nth(1)
        .and_then(|r| r.split([',', '}']).next())
        .and_then(|raw| raw.trim().parse::<u64>().ok())
        .unwrap_or_else(|| panic!("missing u64 `{key}` in {json}"))
}

fn find_artifact_by_kind(ledger_path: &Path, kind: &str) -> (String, Vec<u8>) {
    let ledger = Ledger::open(ledger_path.to_str().unwrap()).unwrap();
    for op in ledger.visible_op_ids(MAIN_BRANCH, None).unwrap() {
        for edge in ledger.op_artifact_edges_bounded(op, 64).unwrap().edges {
            if edge.role == EdgeRole::Out
                && ledger
                    .artifact_info(&edge.artifact)
                    .unwrap()
                    .is_some_and(|i| i.kind == kind)
            {
                let bytes = ledger.get_artifact(&edge.artifact).unwrap().unwrap();
                return (edge.artifact.to_hex(), bytes);
            }
        }
    }
    panic!("no artifact of kind `{kind}` in {}", ledger_path.display());
}

#[test]
#[allow(clippy::too_many_lines)]
fn binary_mesh_export_glb_3mf_vtk_matches_import_receipt_and_enforces_guards() {
    let dir = scratch_dir("mesh");
    let r = ref_dir();
    let db = dir.join("mesh.db");
    let db_s = db.to_str().unwrap();

    run_bin(
        &[
            "--json",
            "import",
            r.join("cooling-reference.fsim").to_str().unwrap(),
            r.join("plate.stl").to_str().unwrap(),
            db_s,
            "--unit",
            "m",
            "--max-hole-edges",
            "0",
        ],
        exit::SUCCESS,
    );

    let (mesh_hash, ply_bytes) = find_artifact_by_kind(&db, "geometry-mesh-ply");
    let (rcpt_hash, rcpt_bytes) = find_artifact_by_kind(&db, "geometry-import-receipt");
    let rcpt_json = String::from_utf8(rcpt_bytes).unwrap();
    let exp_verts = ju64(&rcpt_json, "vertices");
    let exp_tris = ju64(&rcpt_json, "triangles");
    let decoded = fs_io::ply::read_ply(&ply_bytes).expect("read_ply round-trip");
    assert_eq!(decoded.positions.len() as u64, exp_verts);
    assert_eq!(decoded.triangles.len() as u64, exp_tris);

    // 1. GLB (f32-preview)
    let glb_out = dir.join("plate.glb");
    let (gj, _) = run_bin(
        &[
            "--json",
            "mesh",
            "export",
            &mesh_hash,
            db_s,
            "--format",
            "glb",
            "--out",
            glb_out.to_str().unwrap(),
        ],
        exit::SUCCESS,
    );
    let glb_bytes = std::fs::read(&glb_out).unwrap();
    assert!(glb_bytes.starts_with(b"glTF"), "invalid GLB magic");
    assert_eq!(jstr(&gj, "command"), "mesh-export");
    assert_eq!(jstr(&gj, "precision"), "f32-preview");
    assert_eq!(jstr(&gj, "content_hash"), hash_bytes(&glb_bytes).to_hex());
    assert_eq!(
        (ju64(&gj, "vertices"), ju64(&gj, "triangles")),
        (exp_verts, exp_tris)
    );
    assert!(gj.contains("f32 preview"));

    // 2. 3MF (ZIP magic + 3dmodel.model)
    let tmf_out = dir.join("plate.3mf");
    let (mj, _) = run_bin(
        &[
            "--json",
            "mesh",
            "export",
            &mesh_hash,
            db_s,
            "--format",
            "3mf",
            "--out",
            tmf_out.to_str().unwrap(),
        ],
        exit::SUCCESS,
    );
    let tmf_bytes = std::fs::read(&tmf_out).unwrap();
    assert!(
        tmf_bytes.starts_with(b"PK\x03\x04"),
        "invalid 3MF ZIP magic"
    );
    assert!(tmf_bytes.windows(16).any(|w| w == b"3D/3dmodel.model"));
    assert_eq!(jstr(&mj, "precision"), "f64");
    assert_eq!(jstr(&mj, "content_hash"), hash_bytes(&tmf_bytes).to_hex());

    // 3. VTK ASCII (UNSTRUCTURED_GRID: POINTS + CELLS + CELL_TYPES)
    let vtk_out = dir.join("plate.vtk");
    let (vj, _) = run_bin(
        &[
            "--json",
            "mesh",
            "export",
            &mesh_hash,
            db_s,
            "--format",
            "vtk",
            "--out",
            vtk_out.to_str().unwrap(),
        ],
        exit::SUCCESS,
    );
    let vtk_text = std::fs::read_to_string(&vtk_out).unwrap();
    assert!(
        vtk_text.starts_with("# vtk DataFile")
            && vtk_text.contains(&format!("POINTS {exp_verts} double"))
            && vtk_text.contains(&format!("CELLS {exp_tris} "))
            && vtk_text.contains(&format!("CELL_TYPES {exp_tris}"))
    );
    assert_eq!(jstr(&vj, "precision"), "f64");
    assert_eq!(
        jstr(&vj, "content_hash"),
        hash_bytes(vtk_text.as_bytes()).to_hex()
    );

    // 4. Idempotent re-export to identical `--out` succeeds; differing `--out` refuses without clobber.
    run_bin(
        &[
            "--json",
            "mesh",
            "export",
            &mesh_hash,
            db_s,
            "--format",
            "glb",
            "--out",
            glb_out.to_str().unwrap(),
        ],
        exit::SUCCESS,
    );
    let conflict = dir.join("conflict.glb");
    std::fs::write(&conflict, b"do-not-overwrite").unwrap();
    let (_, ce) = run_bin(
        &[
            "--json",
            "mesh",
            "export",
            &mesh_hash,
            db_s,
            "--format",
            "glb",
            "--out",
            conflict.to_str().unwrap(),
        ],
        exit::REFUSED,
    );
    assert!(ce.contains("cli-export-output-conflict"));
    assert_eq!(std::fs::read(&conflict).unwrap(), b"do-not-overwrite");

    // 5. Wrong kind, unknown hash, invalid format, missing ledger (must not create DB).
    let chk_err = |args: &[&str], code: u8, needle: &str| {
        assert!(run_bin(args, code).1.contains(needle));
    };
    let (x_out, u_out, b_out, m_out) = (
        dir.join("x.glb"),
        dir.join("u.glb"),
        dir.join("b.obj"),
        dir.join("m.glb"),
    );
    chk_err(
        &[
            "--json",
            "mesh",
            "export",
            &rcpt_hash,
            db_s,
            "--format",
            "glb",
            "--out",
            x_out.to_str().unwrap(),
        ],
        exit::REFUSED,
        "cli-export-artifact-kind",
    );
    let unk = "c".repeat(64);
    chk_err(
        &[
            "--json",
            "mesh",
            "export",
            &unk,
            db_s,
            "--format",
            "glb",
            "--out",
            u_out.to_str().unwrap(),
        ],
        exit::REFUSED,
        "cli-export-unknown-artifact",
    );
    chk_err(
        &[
            "--json",
            "mesh",
            "export",
            &mesh_hash,
            db_s,
            "--format",
            "obj",
            "--out",
            b_out.to_str().unwrap(),
        ],
        exit::USAGE,
        "cli-mesh-export-usage",
    );
    let missing = dir.join("missing-mesh.db");
    chk_err(
        &[
            "--json",
            "mesh",
            "export",
            &mesh_hash,
            missing.to_str().unwrap(),
            "--format",
            "glb",
            "--out",
            m_out.to_str().unwrap(),
        ],
        exit::INPUT,
        "cli-export-ledger-missing",
    );
    assert!(
        !missing.exists(),
        "mesh export must never create a missing ledger"
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn binary_table_export_csv_and_arrow_match_retained_solution_and_enforce_guards() {
    let dir = scratch_dir("table");
    let r = ref_dir();
    let db = dir.join("solve.db");
    let db_s = db.to_str().unwrap();

    run_bin(
        &[
            "--json",
            "import",
            r.join("cooling-reference.fsim").to_str().unwrap(),
            r.join("plate.stl").to_str().unwrap(),
            db_s,
            "--unit",
            "m",
            "--max-hole-edges",
            "0",
        ],
        exit::SUCCESS,
    );
    run_bin(
        &[
            "--json",
            "run",
            r.join("cooling-reference.fsim").to_str().unwrap(),
            db_s,
            "--materials",
            r.join("aa6061.fsmcdpk").to_str().unwrap(),
        ],
        exit::SUCCESS,
    );

    let (sol_hash, sol_bytes) = find_artifact_by_kind(&db, "solve-conduction-solution");
    let (mesh_hash, _) = find_artifact_by_kind(&db, "geometry-mesh-ply");
    let sol_text = String::from_utf8(sol_bytes).unwrap();
    let temps: Vec<f64> = sol_text
        .split("\"temperature\":[")
        .nth(1)
        .and_then(|t| t.split(']').next())
        .expect("temperature array")
        .split(',')
        .map(|raw| raw.trim().parse::<f64>().expect("f64 temperature"))
        .collect();
    assert!(!temps.is_empty());

    // 1. CSV export: exact header, units, row count, fixed {:.10e} match, honest no_claim.
    let csv_out = dir.join("temps.csv");
    let (cj, _) = run_bin(
        &[
            "--json",
            "table",
            "export",
            &sol_hash,
            db_s,
            "--format",
            "csv",
            "--out",
            csv_out.to_str().unwrap(),
        ],
        exit::SUCCESS,
    );
    assert_eq!(jstr(&cj, "command"), "table-export");
    assert_eq!(ju64(&cj, "row_count") as usize, temps.len());
    assert_eq!(ju64(&cj, "column_count"), 2);
    assert!(cj.contains("{:.10e}"));
    let csv_text = std::fs::read_to_string(&csv_out).unwrap();
    assert_eq!(
        jstr(&cj, "content_hash"),
        hash_bytes(csv_text.as_bytes()).to_hex()
    );
    let lines: Vec<&str> = csv_text.lines().collect();
    assert_eq!(lines[0], "node_index,temperature");
    assert_eq!(lines[1], "1,K");
    assert_eq!(lines.len(), temps.len() + 2);
    for (idx, &t) in temps.iter().enumerate() {
        assert_eq!(lines[idx + 2], format!("{idx},{t:.10e}"));
    }

    // 2. Arrow binary export: ARROW1 magic + 8-byte EOS trailer + bit-exact IEEE-754 f64 bytes.
    let arrow_out = dir.join("temps.arrow");
    let (aj, _) = run_bin(
        &[
            "--json",
            "table",
            "export",
            &sol_hash,
            db_s,
            "--format",
            "arrow",
            "--out",
            arrow_out.to_str().unwrap(),
        ],
        exit::SUCCESS,
    );
    let arrow_bytes = std::fs::read(&arrow_out).unwrap();
    assert!(arrow_bytes.starts_with(b"ARROW1\0\0"));
    assert_eq!(
        &arrow_bytes[arrow_bytes.len() - 8..],
        &[0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0]
    );
    assert_eq!(jstr(&aj, "content_hash"), hash_bytes(&arrow_bytes).to_hex());
    assert_eq!(ju64(&aj, "row_count") as usize, temps.len());
    let data_end = arrow_bytes.len() - 8;
    let f64_tail = &arrow_bytes[data_end - temps.len() * 8..data_end];
    for (idx, &t) in temps.iter().enumerate() {
        let got = f64::from_le_bytes(f64_tail[idx * 8..(idx + 1) * 8].try_into().unwrap());
        assert_eq!(
            got.to_bits(),
            t.to_bits(),
            "bit-exact f64 mismatch at row {idx}"
        );
    }

    // 3. Text mode + idempotent re-export + no-clobber conflict + refusals.
    let (ts, _) = run_bin(
        &[
            "table",
            "export",
            &sol_hash,
            db_s,
            "--format",
            "csv",
            "--out",
            csv_out.to_str().unwrap(),
        ],
        exit::SUCCESS,
    );
    assert!(ts.contains("status=ok\ncommand=table-export\n") && ts.contains("row_count="));

    let conflict = dir.join("conflict.csv");
    std::fs::write(&conflict, "keep-me").unwrap();
    let run_tbl = |art: &str, db_p: &str, fmt: &str, out: &Path, code: u8| {
        run_bin(
            &[
                "--json",
                "table",
                "export",
                art,
                db_p,
                "--format",
                fmt,
                "--out",
                out.to_str().unwrap(),
            ],
            code,
        )
    };
    assert!(
        run_tbl(&sol_hash, db_s, "csv", &conflict, exit::REFUSED)
            .1
            .contains("cli-export-output-conflict")
    );
    assert_eq!(std::fs::read_to_string(&conflict).unwrap(), "keep-me");
    assert!(
        run_tbl(&mesh_hash, db_s, "csv", &dir.join("w.csv"), exit::REFUSED)
            .1
            .contains("cli-export-artifact-kind")
    );
    assert!(
        run_tbl(&sol_hash, db_s, "parquet", &dir.join("p.pq"), exit::USAGE)
            .1
            .contains("cli-table-export-usage")
    );
    let missing = dir.join("missing-table.db");
    assert!(
        run_tbl(
            &sol_hash,
            missing.to_str().unwrap(),
            "csv",
            &dir.join("m.csv"),
            exit::INPUT
        )
        .1
        .contains("cli-export-ledger-missing")
    );
    assert!(
        !missing.exists(),
        "table export must never create a missing ledger"
    );
}
