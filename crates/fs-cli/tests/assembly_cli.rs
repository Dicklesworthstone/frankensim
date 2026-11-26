//! Real-binary, no-mock end-to-end tests for `.11.7` Stage 3 and `.11.10`
//! Stage 1b: `frankensim run <assembly>`, `frankensim run --resume`,
//! `frankensim explain`, FrankenScript `(assembly.run …)` with `:hash` pin
//! enforcement, unchanged legacy `run <project>` and `run <program.fs>` paths,
//! heatsink bit-for-bit identity, and existing-table structured refusal codes.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};

use fs_blake3::hash_bytes;
use fs_cli::{
    AssemblyAnalysisMode, AssemblyCancelPlan, AssemblyStage, CardPackKind, CardPackSet,
    RawCardPack, exit, run_assembly,
};
use fs_ledger::Ledger;

const PINNED_HEATSINK_RUN_ID: &str =
    "5aa11590e72d49ec218c9d5049c71ae19c1ce5e647efa083a08d20d9e21b0a91";
const PINNED_HEATSINK_REPORT_HASH: &str =
    "5fe2c82237d11ea78a526a348f074c9e69311cbe9167e8e8fe92f4113c250399";
const PINNED_HEATSINK_PACKAGE_ROOT: &str =
    "1b42e3d70134c85c021e3e0442c75dd17e5d45a2c11b7ea5b750c533747a80b4";

const SLIDER_FSIM: &str = include_str!("../../../examples/assemblies/thermo_mech_slider.fsim");
const TWO_BODY_FSIM: &str = include_str!("../../../examples/assemblies/two_body_conduction.fsim");

static COUNTER: AtomicUsize = AtomicUsize::new(0);

fn scratch_dir(tag: &str) -> PathBuf {
    let n = COUNTER.fetch_add(1, Ordering::Relaxed);
    let dir = std::env::temp_dir().join(format!("fs-cli-s3-{tag}-{}-{n}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("scratch dir");
    dir
}

fn repo_path(rel: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rel)
}

fn bin_run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_frankensim"))
        .args(args)
        .output()
        .expect("spawn frankensim binary")
}

fn s(bytes: &[u8]) -> String {
    String::from_utf8(bytes.to_vec()).expect("UTF-8")
}

fn jstr<'a>(json: &'a str, key: &str) -> &'a str {
    let needle = format!("\"{key}\":\"");
    json.split(&needle)
        .nth(1)
        .and_then(|r| r.split('"').next())
        .unwrap_or_else(|| panic!("missing `{key}` in {json}"))
}

#[test]
#[allow(clippy::too_many_lines)]
fn binary_runs_three_assemblies_explains_and_resumes_interrupted_run() {
    let dir = scratch_dir("assemblies");
    let steel =
        repo_path("examples/assemblies/packs/steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk");
    let alum = repo_path("examples/assemblies/packs/aluminum-pure-thermal-nbs-c447.fsmcdpk");
    let atkins = repo_path("examples/system/atkins1964/al6061t6-lathe-0p30um-vacuum.fsintpk");
    let sp = repo_path("examples/assemblies/thermo_mech_slider.fsim");
    let tp = repo_path("examples/assemblies/two_body_conduction.fsim");
    let op = repo_path("examples/assemblies/forced_oscillator.fsim");
    let ledger_p = dir.join("assemblies.db");
    let (lstr, steel_s) = (ledger_p.to_str().unwrap(), steel.to_str().unwrap());

    // 1. Assembly 1 (transient) + explain
    let out1 = bin_run(&[
        "--json",
        "run",
        sp.to_str().unwrap(),
        lstr,
        "--materials",
        steel_s,
    ]);
    assert_eq!(
        out1.status.code(),
        Some(i32::from(exit::SUCCESS)),
        "{}",
        s(&out1.stderr)
    );
    let j1 = s(&out1.stdout);
    assert!(j1.contains("\"kind\":\"assembly\"") && j1.contains("\"status\":\"completed\""));
    assert!(j1.contains("\"stages_completed\":6") && j1.contains("\"checker\":\"pass\""));
    let (run1, plan1) = (jstr(&j1, "run"), jstr(&j1, "plan_hash"));
    assert_eq!(plan1.len(), 64);

    let exp1 = bin_run(&["--json", "explain", run1, lstr]);
    assert_eq!(
        exp1.status.code(),
        Some(i32::from(exit::SUCCESS)),
        "{}",
        s(&exp1.stderr)
    );
    let ej1 = s(&exp1.stdout);
    assert_eq!(jstr(&ej1, "plan_hash"), plan1);
    assert!(
        ej1.contains("frankensim.lowering_receipt.v1")
            && ej1.contains("frankensim.window_ledger.v1")
    );
    assert!(ej1.contains("\"causalization_audit\":{"));
    let ps = "\"probe_series\":{\"schema\":\"frankensim.cli.assembly-probe-series.v1\"";
    assert!(ej1.contains(ps) && ej1.contains("\"name\":\"slider_temp\""));
    assert!(
        ej1.contains("\"qoi_budget\":{\"schema\":\"frankensim.cli.assembly-qoi-budget.v1\"")
            && ej1.contains("\"id\":\"temp_limit\"")
            && ej1.contains("\"outcome\":\"indeterminate\"")
            && ej1.contains("unstated card uncertainty")
    );
    let txt1 = s(&bin_run(&["explain", run1, lstr]).stdout);
    assert!(
        txt1.contains("probe_series_present=true")
            && txt1.contains("qoi_budget_present=true")
            && txt1.contains("requirement=temp_limit")
            && txt1.contains("value=NO-DATA")
    );

    // 2. Real-binary `run --resume` completes an interrupted run on disk.
    let partial_p = dir.join("partial.db");
    let partial_lstr = partial_p.to_str().unwrap();
    let partial_run_id = {
        let ledger = Ledger::open(partial_lstr).unwrap();
        let packs = CardPackSet::admit(vec![RawCardPack {
            kind: CardPackKind::Material,
            source: steel_s.into(),
            bytes: std::fs::read(&steel).unwrap(),
            expect: None,
        }])
        .unwrap();
        let mut t = 0.0;
        let err = run_assembly(
            &ledger,
            SLIDER_FSIM.as_bytes(),
            &packs,
            AssemblyAnalysisMode::Transient,
            &mut || {
                t += 0.1;
                t
            },
            AssemblyCancelPlan {
                cancel_pre_publish: Some(AssemblyStage::Plan),
                cancel_after: None,
            },
        )
        .unwrap_err();
        assert_eq!((err.code, err.stage), ("cli-solve-cancelled", Some("plan")));
        err.run.expect("cancelled prefix must carry run id")
    };
    let res = bin_run(&["--json", "run", "--resume", &partial_run_id, partial_lstr]);
    assert_eq!(
        res.status.code(),
        Some(i32::from(exit::SUCCESS)),
        "{}",
        s(&res.stderr)
    );
    let rj = s(&res.stdout);
    assert!(rj.contains("\"status\":\"completed\"") && rj.contains("\"stages_completed\":6"));
    assert_eq!((jstr(&rj, "run"), jstr(&rj, "plan_hash")), (run1, plan1));

    // 3. Assembly 2 (two-body conduction with tracked NBS C447 Al + Atkins 1964 joint)
    let out2 = bin_run(&[
        "--json",
        "run",
        tp.to_str().unwrap(),
        lstr,
        "--materials",
        alum.to_str().unwrap(),
        "--interfaces",
        atkins.to_str().unwrap(),
    ]);
    assert_eq!(
        out2.status.code(),
        Some(i32::from(exit::SUCCESS)),
        "{}",
        s(&out2.stderr)
    );
    assert!(s(&out2.stdout).contains("\"status\":\"completed\""));

    // 4. Assembly 3 under modal, harmonic, and stationary modes
    let mut runs = Vec::new();
    for mode in ["modal", "harmonic", "stationary"] {
        let out = bin_run(&[
            "--json",
            "run",
            op.to_str().unwrap(),
            lstr,
            "--materials",
            steel_s,
            "--analysis",
            mode,
        ]);
        assert_eq!(
            out.status.code(),
            Some(i32::from(exit::SUCCESS)),
            "{mode}: {}",
            s(&out.stderr)
        );
        let j = s(&out.stdout);
        assert!(j.contains("\"status\":\"completed\"") && j.contains("\"checker\":\"pass\""));
        assert_eq!(jstr(&j, "analysis"), mode);
        runs.push(jstr(&j, "run").to_string());
    }
    assert!(runs[0] != runs[1] && runs[1] != runs[2]);
    let exp_s = bin_run(&["--json", "explain", &runs[2], lstr]);
    assert_eq!(exp_s.status.code(), Some(i32::from(exit::SUCCESS)));
    assert!(s(&exp_s.stdout).contains("frankensim.lowering_receipt.v1"));
}

#[test]
fn binary_frankenscript_assembly_run_matches_direct_cli_and_enforces_hash_pin() {
    let dir = scratch_dir("frankenscript_asm");
    let steel =
        repo_path("examples/assemblies/packs/steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk");
    let sp = repo_path("examples/assemblies/thermo_mech_slider.fsim");
    let (dl, pl) = (dir.join("direct.db"), dir.join("prog.db"));

    let direct = bin_run(&[
        "--json",
        "run",
        sp.to_str().unwrap(),
        dl.to_str().unwrap(),
        "--materials",
        steel.to_str().unwrap(),
    ]);
    assert_eq!(
        direct.status.code(),
        Some(i32::from(exit::SUCCESS)),
        "{}",
        s(&direct.stderr)
    );
    let dj = s(&direct.stdout);
    let (d_run, d_plan) = (jstr(&dj, "run"), jstr(&dj, "plan_hash"));

    let canonical =
        fs_ir::sexpr::print_checked(&fs_ir::sexpr::parse(SLIDER_FSIM).unwrap()).unwrap();
    let exact_hash = hash_bytes(canonical.as_bytes()).to_hex();
    let prog_ok = dir.join("run_slider.fs");
    std::fs::write(
        &prog_ok,
        format!(
            "(study \"slider_via_fs\" (seed 0x7) (versions (constellation :lock \"2026-07\")) (budget (wall 60s) (mem 64MiB)) (capability :cores 1 :mem 64MiB :wall 60s :ops (assembly.*)) (assembly.run \"{}\" :materials (\"{}\") :analysis \"transient\" :hash \"{exact_hash}\"))\n",
            sp.display(),
            steel.display(),
        ),
    )
    .unwrap();
    let via_fs = bin_run(&[
        "--json",
        "run",
        prog_ok.to_str().unwrap(),
        pl.to_str().unwrap(),
    ]);
    assert_eq!(
        via_fs.status.code(),
        Some(i32::from(exit::SUCCESS)),
        "{}",
        s(&via_fs.stderr)
    );
    let fj = s(&via_fs.stdout);
    assert!(fj.contains("\"verb\":\"assembly.run\"") && fj.contains("\"status\":\"completed\""));
    assert!(fj.contains(d_run) && fj.contains(d_plan), "mismatch: {fj}");

    let prog_bad = dir.join("bad_hash.fs");
    std::fs::write(
        &prog_bad,
        format!(
            "(study \"bad_hash_study\" (seed 0x7) (versions (constellation :lock \"2026-07\")) (budget (wall 60s) (mem 64MiB)) (capability :cores 1 :mem 64MiB :wall 60s :ops (assembly.*)) (assembly.run \"{}\" :materials (\"{}\") :hash \"{}\"))\n",
            sp.display(),
            steel.display(),
            "0".repeat(64)
        ),
    )
    .unwrap();
    let bad = bin_run(&[
        "--json",
        "run",
        prog_bad.to_str().unwrap(),
        pl.to_str().unwrap(),
    ]);
    assert_eq!(bad.status.code(), Some(i32::from(exit::REFUSED)));
    assert!(s(&bad.stderr).contains("frankenscript-assembly-hash"));
}

#[test]
fn binary_legacy_run_project_and_run_program_fs_paths_remain_unchanged() {
    let fsim = repo_path("examples/heatsink-fan/heatsink-fan.fsim");
    let stl = repo_path("examples/heatsink-fan/heatsink.stl");
    let pack = repo_path("data/reference-project/aa6061.fsmcdpk");

    // Old form 1: `frankensim run <project.fsim> <ledger>` (heatsink-fan identity).
    let dir1 = scratch_dir("legacy_project");
    let l1 = dir1.join("heatsink.db");
    let imp = bin_run(&[
        "--json",
        "import",
        fsim.to_str().unwrap(),
        stl.to_str().unwrap(),
        l1.to_str().unwrap(),
        "--unit",
        "m",
        "--max-hole-edges",
        "0",
    ]);
    assert_eq!(
        imp.status.code(),
        Some(i32::from(exit::SUCCESS)),
        "{}",
        s(&imp.stderr)
    );
    let out1 = bin_run(&[
        "--json",
        "run",
        fsim.to_str().unwrap(),
        l1.to_str().unwrap(),
        "--materials",
        pack.to_str().unwrap(),
    ]);
    assert_eq!(
        out1.status.code(),
        Some(i32::from(exit::SUCCESS)),
        "{}",
        s(&out1.stderr)
    );
    let j1 = s(&out1.stdout);
    assert!(j1.contains("\"stages_completed\":7") && j1.contains("\"checker\":\"pass\""));
    assert_eq!(jstr(&j1, "run"), PINNED_HEATSINK_RUN_ID);
    assert_eq!(
        jstr(&j1, "report_content_hash"),
        PINNED_HEATSINK_REPORT_HASH
    );
    assert_eq!(jstr(&j1, "package_root"), PINNED_HEATSINK_PACKAGE_ROOT);

    // Old form 2: `frankensim run <program.fs> <ledger>` cooling FrankenScript program.
    let dir2 = scratch_dir("legacy_program");
    let l2 = dir2.join("cooling_prog.db");
    let prog = dir2.join("heatsink.fs");
    std::fs::write(
        &prog,
        format!(
            "(study \"legacy_cooling_program\" (seed 0x7) (versions (constellation :lock \"2026-07-06\")) (budget (wall 60s) (mem 64MiB)) (capability :cores 1 :mem 64MiB :wall 60s :ops (cooling.*)) (let p (cooling.project \"{}\")) (cooling.import p :sources (\"{}\") :unit \"m\" :max-hole-edges 0) (cooling.run p :materials (\"{}\")))\n",
            fsim.display(), stl.display(), pack.display(),
        ),
    )
    .unwrap();
    let out2 = bin_run(&[
        "--json",
        "run",
        prog.to_str().unwrap(),
        l2.to_str().unwrap(),
    ]);
    assert_eq!(
        out2.status.code(),
        Some(i32::from(exit::SUCCESS)),
        "{}",
        s(&out2.stderr)
    );
    let j2 = s(&out2.stdout);
    assert!(j2.contains("\"verb\":\"cooling.run\"") && j2.contains("\"status\":\"completed\""));
    assert!(j2.contains(PINNED_HEATSINK_RUN_ID));
}

#[test]
fn binary_structured_refusals_use_existing_exit_codes() {
    let dir = scratch_dir("refusals");
    let l = dir.join("refusals.db");
    let lstr = l.to_str().unwrap();
    let sp = repo_path("examples/assemblies/thermo_mech_slider.fsim");
    let hs = repo_path("examples/heatsink-fan/heatsink-fan.fsim");
    let unk = "a".repeat(64);

    let e1 = bin_run(&["--json", "explain", &unk]);
    assert_eq!(e1.status.code(), Some(i32::from(exit::INPUT)));
    assert!(s(&e1.stderr).contains("cli-export-ledger-required"));

    let e2 = bin_run(&["--json", "run", sp.to_str().unwrap(), lstr]);
    assert_eq!(e2.status.code(), Some(i32::from(exit::REFUSED)));
    assert!(s(&e2.stderr).contains("cli-assembly-card-missing"));

    let e3 = bin_run(&["--json", "explain", &unk, lstr]);
    assert_eq!(e3.status.code(), Some(i32::from(exit::REFUSED)));
    assert!(s(&e3.stderr).contains("cli-solve-unknown-run"));

    let e4 = bin_run(&[
        "--json",
        "run",
        sp.to_str().unwrap(),
        lstr,
        "--analysis",
        "bogus",
    ]);
    assert_eq!(e4.status.code(), Some(i32::from(exit::USAGE)));
    assert!(s(&e4.stderr).contains("cli-run-analysis-usage"));

    let e5 = bin_run(&[
        "--json",
        "run",
        hs.to_str().unwrap(),
        lstr,
        "--analysis",
        "modal",
    ]);
    assert_eq!(e5.status.code(), Some(i32::from(exit::USAGE)));
    assert!(s(&e5.stderr).contains("cli-run-analysis-project"));
}

#[test]
fn binary_assembly_run_two_ledgers_pass_replay_audit_cleanly() {
    let dir = scratch_dir("assembly_replay");
    let steel =
        repo_path("examples/assemblies/packs/steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk");
    let sp = repo_path("examples/assemblies/thermo_mech_slider.fsim");
    let (db_a, db_b) = (dir.join("asm_a.db"), dir.join("asm_b.db"));
    for (idx, db) in [&db_a, &db_b].into_iter().enumerate() {
        if idx == 1 {
            std::thread::sleep(std::time::Duration::from_millis(25));
        }
        let out = bin_run(&[
            "--json",
            "run",
            sp.to_str().unwrap(),
            db.to_str().unwrap(),
            "--materials",
            steel.to_str().unwrap(),
        ]);
        assert_eq!(
            out.status.code(),
            Some(i32::from(exit::SUCCESS)),
            "{}",
            s(&out.stderr)
        );
        assert!(s(&out.stdout).contains("\"status\":\"completed\""));
    }
    let audit = bin_run(&[
        "--json",
        "replay-audit",
        db_a.to_str().unwrap(),
        db_b.to_str().unwrap(),
    ]);
    assert_eq!(
        audit.status.code(),
        Some(i32::from(exit::SUCCESS)),
        "{}",
        s(&audit.stderr)
    );
    let aj = s(&audit.stdout);
    for needle in [
        "\"status\":\"clean\"",
        "\"clean\":true",
        "\"compared\":7",
        "\"deterministic_mismatches\":[]",
        "\"structure_mismatch\":null",
    ] {
        assert!(aj.contains(needle), "missing {needle}: {aj}");
    }
    let ledger = Ledger::open(db_a.to_str().unwrap()).unwrap();
    let ops = ledger.visible_op_ids(fs_ledger::MAIN_BRANCH, None).unwrap();
    assert_eq!(ops.len(), 7);
    for (idx, &op) in ops[0..6].iter().enumerate() {
        let t_end = i64::try_from(idx).unwrap() * 2 + 1;
        let m = ledger.metrics_for_op(op).unwrap();
        assert!(
            m.iter()
                .any(|(t, n, _)| *t == t_end && n == "driver.consumed_wall_s")
                && m.iter()
                    .any(|(t, n, _)| *t == t_end && n == "driver.consumed_core_s")
        );
        assert!(
            m.iter()
                .all(|(t, n, _)| n.starts_with("driver.") || *t != t_end)
        );
    }
    let tm = ledger.metrics_for_op(ops[6]).unwrap();
    assert!(
        tm.iter()
            .any(|(t, n, _)| *t == 101 && n == "driver.consumed_wall_s")
            && tm
                .iter()
                .any(|(t, n, _)| *t == 101 && n == "driver.consumed_core_s")
    );
}

#[test]
fn binary_lowers_atkins_contact_resistance_and_refuses_missing_area() {
    let dir = scratch_dir("atkins_tcr");
    let alum = repo_path("examples/assemblies/packs/aluminum-pure-thermal-nbs-c447.fsmcdpk");
    let atk = repo_path("examples/system/atkins1964/al6061t6-lathe-0p30um-vacuum.fsintpk");
    let tp = repo_path("examples/assemblies/two_body_conduction.fsim");
    let (m, ic) = (alum.to_str().unwrap(), atk.to_str().unwrap());
    let db1 = dir.join("ok.db");
    let r = bin_run(&[
        "--json",
        "run",
        tp.to_str().unwrap(),
        db1.to_str().unwrap(),
        "--materials",
        m,
        "--interfaces",
        ic,
    ]);
    assert_eq!(r.status.code(), Some(i32::from(exit::SUCCESS)));
    let id = jstr(&s(&r.stdout), "run").to_owned();
    let ej = s(&bin_run(&["--json", "explain", &id, db1.to_str().unwrap()]).stdout);
    assert!(ej.contains("area-specific-thermal-contact-resistance") && ej.contains("7.909890527"));

    let bad_f = dir.join("missing_area.fsim");
    std::fs::write(&bad_f, TWO_BODY_FSIM.replace(":area 0.0025m^2 ", "")).unwrap();
    let db2 = dir.join("bad.db");
    let r2 = bin_run(&[
        "--json",
        "run",
        bad_f.to_str().unwrap(),
        db2.to_str().unwrap(),
        "--materials",
        m,
        "--interfaces",
        ic,
    ]);
    assert_eq!(r2.status.code(), Some(i32::from(exit::REFUSED)));
    assert!(s(&r2.stderr).contains("requires positive :area"));
}
