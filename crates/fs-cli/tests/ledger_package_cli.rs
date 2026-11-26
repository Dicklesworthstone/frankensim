//! Real-binary, no-mock end-to-end tests for bead `frankensim-dhqf7`:
//! `lineage`, `replay-audit`, and `package --coverage <standard>`.

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicUsize, Ordering};

use fs_cli::exit;
use fs_ledger::{EdgeRole, ExecMode, FiveExplicits, Ledger, MAIN_BRANCH, OpOutcome};
use fs_package::{CoverageStatus, EvidencePackage, Standard, verified_package_coverage};

const HS_RUN_ID: &str = "5aa11590e72d49ec218c9d5049c71ae19c1ce5e647efa083a08d20d9e21b0a91";
const HS_PKG_ROOT: &str = "1b42e3d70134c85c021e3e0442c75dd17e5d45a2c11b7ea5b750c533747a80b4";
static COUNTER: AtomicUsize = AtomicUsize::new(0);

fn scratch_dir(tag: &str) -> PathBuf {
    let n = COUNTER.fetch_add(1, Ordering::Relaxed);
    let dir = std::env::temp_dir().join(format!("fs-cli-dhqf7-{tag}-{}-{n}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("scratch dir");
    dir
}

fn repo_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
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

fn import_ref_plate(db: &Path) {
    let ref_dir = repo_root().join("data/reference-project");
    let fsim = ref_dir.join("cooling-reference.fsim");
    let stl = ref_dir.join("plate.stl");
    run_bin(
        &[
            "--json",
            "import",
            fsim.to_str().unwrap(),
            stl.to_str().unwrap(),
            db.to_str().unwrap(),
            "--unit",
            "m",
            "--max-hole-edges",
            "0",
        ],
        exit::SUCCESS,
    );
}

fn run_ref_cooling(db: &Path) -> String {
    let ref_dir = repo_root().join("data/reference-project");
    let fsim = ref_dir.join("cooling-reference.fsim");
    let pack = ref_dir.join("aa6061.fsmcdpk");
    import_ref_plate(db);
    let (rj, _) = run_bin(
        &[
            "--json",
            "run",
            fsim.to_str().unwrap(),
            db.to_str().unwrap(),
            "--materials",
            pack.to_str().unwrap(),
        ],
        exit::SUCCESS,
    );
    jstr(&rj, "run").to_string()
}

fn record_det_op(ledger: &Ledger, branch: i64, ir: &str, inp: &[u8], out: &[u8], t: i64) -> String {
    let fx = FiveExplicits {
        seed: &[7, 0, 0, 0],
        versions: "{\"lock\":\"2026-07\"}",
        budget: "{\"wall_s\":60}",
        capability: "{\"cores\":1}",
    };
    let in_art = ledger.put_artifact("input-bytes", inp, None).unwrap();
    let op = ledger
        .begin_op_on(branch, ExecMode::Deterministic, None, ir, &fx, t)
        .unwrap();
    ledger.link(op, &in_art.hash, EdgeRole::In).unwrap();
    let out_art = ledger.put_artifact("output-bytes", out, None).unwrap();
    ledger.link(op, &out_art.hash, EdgeRole::Out).unwrap();
    ledger.finish_op(op, OpOutcome::Ok, None, t + 10).unwrap();
    out_art.hash.to_hex()
}

#[test]
fn binary_lineage_walks_real_dag_truncates_depth_and_never_creates_missing_db() {
    let dir = scratch_dir("lineage");
    let db = dir.join("cooling.db");
    let db_s = db.to_str().unwrap();
    run_ref_cooling(&db);

    // Pick the latest op that links both input and output artifacts (multi-hop DAG).
    let target_hash = {
        let ledger = Ledger::open(db_s).unwrap();
        let ops = ledger.visible_op_ids(MAIN_BRANCH, None).unwrap();
        ops.into_iter()
            .rev()
            .find_map(|op| {
                let edges = ledger.op_artifact_edges_bounded(op, 64).unwrap().edges;
                let has_in = edges.iter().any(|e| e.role == EdgeRole::In);
                has_in.then(|| {
                    edges
                        .into_iter()
                        .find(|e| e.role == EdgeRole::Out)
                        .map(|e| e.artifact.to_hex())
                })?
            })
            .expect("op with both input and output edges")
    };

    let (fj, _) = run_bin(&["--json", "lineage", &target_hash, db_s], exit::SUCCESS);
    assert!(fj.contains("\"command\":\"lineage\"") && fj.contains("\"status\":\"ok\""));
    assert!(fj.contains("\"truncated\":false") && fj.contains("\"producer_count\":1"));
    assert!(fj.matches("\"produced_by\":[{").count() > 1);

    let (sj, _) = run_bin(
        &["--json", "lineage", &target_hash, db_s, "--depth", "1"],
        exit::SUCCESS,
    );
    assert!(sj.contains("\"max_depth\":1") && sj.contains("\"truncated\":true"));
    assert_eq!(sj.matches("\"produced_by\":[{").count(), 1);

    let (ts, _) = run_bin(&["lineage", &target_hash, db_s], exit::SUCCESS);
    assert!(ts.contains("status=ok\ncommand=lineage\n") && ts.contains("<- op "));

    let unk = "b".repeat(64);
    assert!(
        run_bin(&["--json", "lineage", &unk, db_s], exit::REFUSED)
            .1
            .contains("cli-lineage-unknown-artifact")
    );
    assert!(
        run_bin(&["--json", "lineage", "not-hex", db_s], exit::USAGE)
            .1
            .contains("cli-lineage-hash")
    );
    assert!(
        run_bin(
            &["--json", "lineage", &target_hash, db_s, "--depth", "0"],
            exit::USAGE
        )
        .1
        .contains("cli-lineage-usage")
    );

    let missing = dir.join("does-not-exist.db");
    assert!(
        run_bin(
            &["--json", "lineage", &target_hash, missing.to_str().unwrap()],
            exit::INPUT
        )
        .1
        .contains("cli-export-ledger-missing")
    );
    assert!(
        !missing.exists(),
        "lineage must never create a missing ledger file"
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn binary_replay_audit_clean_mismatch_branch_diff_and_no_file_creation() {
    let dir = scratch_dir("replay");
    let (db_a, db_b) = (dir.join("run_a.db"), dir.join("run_b.db"));
    import_ref_plate(&db_a);
    import_ref_plate(&db_b);

    let (cj, _) = run_bin(
        &[
            "--json",
            "replay-audit",
            db_a.to_str().unwrap(),
            db_b.to_str().unwrap(),
        ],
        exit::SUCCESS,
    );
    assert!(
        cj.contains("\"command\":\"replay-audit\"")
            && cj.contains("\"status\":\"clean\"")
            && cj.contains("\"clean\":true")
            && cj.contains("\"compared\":1")
            && cj.contains("\"deterministic_mismatches\":[]")
    );

    let (det_a, det_b) = (dir.join("det_a.db"), dir.join("det_b.db"));
    let (ha, hb) = {
        let la = Ledger::open(det_a.to_str().unwrap()).unwrap();
        let lb = Ledger::open(det_b.to_str().unwrap()).unwrap();
        (
            record_det_op(&la, MAIN_BRANCH, "{\"op\":\"step\"}", b"in1", b"out-A", 10),
            record_det_op(&lb, MAIN_BRANCH, "{\"op\":\"step\"}", b"in1", b"out-B", 10),
        )
    };
    let (dj, de) = run_bin(
        &[
            "--json",
            "replay-audit",
            det_a.to_str().unwrap(),
            det_b.to_str().unwrap(),
        ],
        exit::REFUSED,
    );
    assert!(
        dj.contains("\"status\":\"diverged\"")
            && dj.contains("\"clean\":false")
            && dj.contains(&ha)
            && dj.contains(&hb)
    );
    assert!(de.contains("cli-replay-audit-mismatch"));

    {
        let la = Ledger::open(det_a.to_str().unwrap()).unwrap();
        let alt = la.fork("alt", MAIN_BRANCH).unwrap();
        record_det_op(&la, alt, "{\"op\":\"extra\"}", b"in2", b"out-alt", 30);
    }
    let (bj, _) = run_bin(
        &[
            "--json",
            "replay-audit",
            det_a.to_str().unwrap(),
            det_a.to_str().unwrap(),
            "--branch-a",
            "main",
            "--branch-b",
            "alt",
        ],
        exit::REFUSED,
    );
    assert!(
        bj.contains("\"branch_diff\":{\"shared\":1,\"only_a\":[],\"only_b\":[2]}")
            && bj.contains("\"branch_b\":\"alt\"")
    );

    let (miss_a, miss_b) = (dir.join("missing_a.db"), dir.join("missing_b.db"));
    assert!(
        run_bin(
            &[
                "--json",
                "replay-audit",
                miss_a.to_str().unwrap(),
                db_b.to_str().unwrap()
            ],
            exit::INPUT
        )
        .1
        .contains("cli-export-ledger-missing")
    );
    assert!(
        !miss_a.exists(),
        "replay-audit must never create missing a.db"
    );
    assert!(
        run_bin(
            &[
                "--json",
                "replay-audit",
                db_a.to_str().unwrap(),
                miss_b.to_str().unwrap()
            ],
            exit::INPUT
        )
        .1
        .contains("cli-export-ledger-missing")
    );
    assert!(
        !miss_b.exists(),
        "replay-audit must never create missing b.db"
    );
    assert!(
        run_bin(
            &[
                "--json",
                "replay-audit",
                db_a.to_str().unwrap(),
                db_b.to_str().unwrap(),
                "--branch-a",
                "missing-br"
            ],
            exit::REFUSED
        )
        .1
        .contains("cli-replay-audit-branch")
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn binary_package_coverage_all_four_standards_default_unchanged_and_study_refused() {
    let root = repo_root();
    let fsim = root.join("examples/heatsink-fan/heatsink-fan.fsim");
    let stl = root.join("examples/heatsink-fan/heatsink.stl");
    let pack = root.join("data/reference-project/aa6061.fsmcdpk");
    let dir = scratch_dir("pkg_cov");
    let db = dir.join("heatsink.db");
    let db_s = db.to_str().unwrap();

    run_bin(
        &[
            "--json",
            "import",
            fsim.to_str().unwrap(),
            stl.to_str().unwrap(),
            db_s,
            "--unit",
            "m",
            "--max-hole-edges",
            "0",
        ],
        exit::SUCCESS,
    );
    let (rj, _) = run_bin(
        &[
            "--json",
            "run",
            fsim.to_str().unwrap(),
            db_s,
            "--materials",
            pack.to_str().unwrap(),
        ],
        exit::SUCCESS,
    );
    assert_eq!(jstr(&rj, "run"), HS_RUN_ID);

    let (bj, _) = run_bin(&["--json", "package", HS_RUN_ID, db_s], exit::SUCCESS);
    assert_eq!(jstr(&bj, "merkle_root"), HS_PKG_ROOT);
    assert!(bj.contains("\"checker\":\"pass\"") && !bj.contains("\"coverage\""));

    let pkg = EvidencePackage::from_json(&std::fs::read_to_string(jstr(&bj, "package")).unwrap())
        .unwrap();
    let verified = pkg
        .into_verified()
        .expect("heatsink package verifies under deny-all");

    let mut hashes = Vec::new();
    for (slug, std_enum) in [
        ("asme-vv-10", Standard::AsmeVvV10),
        ("asme-vv-20", Standard::AsmeVvV20),
        ("asme-vv-40", Standard::AsmeVvV40),
        ("faa-easa-cba", Standard::FaaEasaCbA),
    ] {
        let (j, _) = run_bin(
            &["--json", "package", HS_RUN_ID, db_s, "--coverage", slug],
            exit::SUCCESS,
        );
        let expected = verified_package_coverage(&verified, std_enum);
        let expected_hash = expected.decision_hash().to_hex();
        let (mut cov, mut abs, mut noc) = (0usize, 0usize, 0usize);
        for (concept, st, _) in expected.rows() {
            let st_slug = match st {
                CoverageStatus::Covered => {
                    cov += 1;
                    "covered"
                }
                CoverageStatus::MappedButAbsent => {
                    abs += 1;
                    "mapped-but-absent"
                }
                CoverageStatus::NoClaim => {
                    noc += 1;
                    "no-claim"
                }
            };
            assert!(j.contains(&format!(
                "\"concept\":\"{}\",\"status\":\"{st_slug}\"",
                concept.label()
            )));
        }
        assert_eq!(jstr(&j, "standard"), slug);
        assert_eq!(jstr(&j, "decision_hash"), expected_hash);
        assert!(j.contains("\"crosswalk_version\":6") && j.contains("\"package_format\":9"));
        assert!(j.contains(&format!(
            "\"summary\":{{\"covered\":{cov},\"mapped_but_absent\":{abs},\"no_claim\":{noc}}}"
        )));
        hashes.push(expected_hash);
    }
    hashes.sort();
    hashes.dedup();
    assert_eq!(hashes.len(), 4);

    let (ts, _) = run_bin(
        &["package", HS_RUN_ID, db_s, "--coverage", "asme-vv-20"],
        exit::SUCCESS,
    );
    assert!(
        ts.contains("coverage_standard=asme-vv-20\n") && ts.contains("coverage_decision_hash=")
    );
    assert!(
        run_bin(
            &[
                "--json",
                "package",
                HS_RUN_ID,
                db_s,
                "--coverage",
                "iso-9999"
            ],
            exit::USAGE
        )
        .1
        .contains("cli-package-coverage-usage")
    );
    assert!(
        run_bin(
            &[
                "--json",
                "package",
                "study-0000000000000000",
                db_s,
                "--coverage",
                "asme-vv-20"
            ],
            exit::REFUSED
        )
        .1
        .contains("cli-package-study-coverage-unsupported")
    );
}
