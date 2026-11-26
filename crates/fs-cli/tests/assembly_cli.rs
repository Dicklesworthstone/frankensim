//! Real-binary, no-mock end-to-end tests for `.11.7` Stage 3:
//! `frankensim run <assembly>`, `frankensim run --resume`, `frankensim explain`,
//! FrankenScript `(assembly.run …)` with `:hash` pin enforcement, unchanged
//! legacy `run <project>` and `run <program.fs>` paths, heatsink bit-for-bit
//! identity, and existing-table structured refusal exit codes.

use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};

use fs_blake3::hash_bytes;
use fs_cli::{
    AssemblyAnalysisMode, AssemblyCancelPlan, AssemblyStage, CardPackKind, CardPackSet,
    RawCardPack, exit, run_assembly,
};
use fs_evidence::ValidityDomain;
use fs_ledger::Ledger;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedInterfacePack,
    NormalizedMaterialCardPack, NormalizedPack, ObservationDataset, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_qty::Dims;

const PINNED_HEATSINK_RUN_ID: &str =
    "5aa11590e72d49ec218c9d5049c71ae19c1ce5e647efa083a08d20d9e21b0a91";
const PINNED_HEATSINK_REPORT_HASH: &str =
    "5fe2c82237d11ea78a526a348f074c9e69311cbe9167e8e8fe92f4113c250399";
const PINNED_HEATSINK_PACKAGE_ROOT: &str =
    "1b42e3d70134c85c021e3e0442c75dd17e5d45a2c11b7ea5b750c533747a80b4";

static COUNTER: AtomicUsize = AtomicUsize::new(0);

fn scratch_dir(tag: &str) -> PathBuf {
    let n = COUNTER.fetch_add(1, Ordering::Relaxed);
    let dir = std::env::temp_dir().join(format!("fs-cli-s3-{tag}-{}-{n}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).expect("scratch dir");
    dir
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

fn make_claims(tag: &str, entries: &[(&str, f64, Dims)]) -> NormalizedPack {
    let prov = |src: &str| Provenance {
        source: src.into(),
        license: "CC-BY-4.0".into(),
        artifact: Some(hash_bytes(src.as_bytes())),
    };
    let mut claims = ClaimSet::new();
    let obs = claims
        .register_observation(ObservationDataset {
            specimen: format!("coupon-{tag}"),
            method: "test-fixture".into(),
            artifact: hash_bytes(tag.as_bytes()),
            caveats: "test fixture".into(),
            provenance: prov(tag),
        })
        .unwrap();
    for &(name, value, dims) in entries {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 200.0, 600.0),
                uncertainty: UncertaintyModel::HalfWidth {
                    half_width: (value.abs() * 0.01).max(1.0e-6),
                    confidence: 0.95,
                },
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![obs],
                provenance: prov(name),
            })
            .unwrap();
    }
    NormalizedPack::new(
        format!("pack-{tag}"),
        "fs-cli-s3-test-v1",
        hash_bytes(tag.as_bytes()),
        "CC-BY-4.0",
        claims,
        Vec::new(),
        Vec::new(),
    )
    .unwrap()
}

fn write_mat(dir: &Path, chem: &str, rho: f64, cp: f64, youngs: f64) -> PathBuf {
    let norm = make_claims(
        chem,
        &[
            ("density", rho, Dims([-3, 1, 0, 0, 0, 0])),
            ("specific_heat_capacity", cp, Dims([2, 0, -2, -1, 0, 0])),
            ("youngs_modulus", youngs, Dims([-1, 1, -2, 0, 0, 0])),
        ],
    );
    let pack = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: chem.into(),
            phase: "wrought".into(),
            process: "annealed".into(),
            revision: 0,
        },
        norm,
    )
    .unwrap();
    let p = dir.join(format!("{chem}.fsmcdpk"));
    std::fs::write(&p, pack.to_bytes()).unwrap();
    p
}

fn write_tim(dir: &Path, g_w_k: f64) -> PathBuf {
    let norm = make_claims(
        "tim",
        &[("thermal_conductance", g_w_k, Dims([2, 1, -3, -1, 0, 0]))],
    );
    let surf = |chem: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: chem.into(),
            phase: "wrought".into(),
            process: "annealed".into(),
            revision: 0,
        },
        texture_frame: "ground".into(),
    };
    let pack = NormalizedInterfacePack::new(
        surf("AISI-304"),
        surf("AA6061"),
        SystemContext {
            medium: "dry-contact".into(),
            third_body: None,
            environment: "air".into(),
            history: "run-in".into(),
        },
        norm,
    )
    .unwrap();
    let p = dir.join("contact.fsintpk");
    std::fs::write(&p, pack.to_bytes()).unwrap();
    p
}

const SLIDER_FSIM: &str = r#"(assembly "thermo_mech_slider"
  (version 1) (units "SI") (seed 0x7)
  (versions (constellation "2026-07") (workspace "stage3-cli"))
  (budget (wall 60s) (mem 64MiB))
  (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*))
  (reference-state unloaded) (environment earth-lab-air)
  (time-plan :horizon 0.05s :dt 0.001s :sample-rate 1000Hz)
  (part "slider" :kind rigid-body :geometry (box 0.1m 0.05m 0.05m)
    :material "AISI-304" :query-temp 293.15K :physics (mechanics thermal) :fidelity lumped)
  (joint "suspension" :kind prismatic :parent root :child "slider"
    :axis (1.0 0.0 0.0) :damping 12N*s/m :stiffness 400N/m)
  (source "drive" :part "slider" :port "mechanics"
    :history (harmonic :amplitude 15N :frequency 3Hz :phase 0.25 :offset 5N))
  (probe "slider_temp" :kind thermocouple :target "slider" :sample-rate 1000Hz :time-constant 0.005s)
)"#;

const TWO_BODY_FSIM: &str = r#"(assembly "two_body_conduction"
  (version 1) (units "SI") (seed 0x11)
  (versions (constellation "2026-07") (workspace "stage3-cli"))
  (budget (wall 60s) (mem 64MiB))
  (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*))
  (reference-state assembled) (environment :ambient-temp 300K :ambient-pressure 101325Pa)
  (time-plan :horizon 0.1s :dt 0.005s :sample-rate 200Hz)
  (part "hot" :kind solid :geometry (box 0.05m 0.05m 0.05m)
    :material "AISI-304" :query-temp 300K :physics (thermal) :fidelity lumped)
  (part "cold" :kind solid :geometry (box 0.05m 0.05m 0.05m)
    :material "AA6061" :query-temp 300K :physics (thermal) :fidelity lumped)
  (interface "joint_tim" :kind thermal-contact :a "hot" :b "cold"
    :card "AISI-304" :family "thermal-conductance" :area 0.0025m^2)
  (source "heater" :part "hot" :port "thermal" :history (constant 500W))
  (probe "hot_temp" :kind thermocouple :target "hot" :sample-rate 200Hz)
)"#;

const FORCED_OSC_FSIM: &str = r#"(assembly "forced_oscillator"
  (version 1) (units "SI") (seed 0x22)
  (versions (constellation "2026-07") (workspace "stage3-cli"))
  (budget (wall 60s) (mem 64MiB))
  (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*))
  (reference-state unloaded) (environment earth-lab-air)
  (time-plan :horizon 0.05s :dt 0.001s :sample-rate 1000Hz)
  (part "mass" :kind rigid-body :geometry (box 0.1m 0.05m 0.05m)
    :material "AISI-304" :query-temp 293.15K :physics (mechanics) :fidelity lumped)
  (joint "spring_damper" :kind prismatic :parent root :child "mass"
    :axis (1.0 0.0 0.0) :damping 4N*s/m :stiffness 200N/m)
  (source "load" :part "mass" :port "mechanics" :history (constant 10N))
  (analysis-options :harmonic-frequencies (1.0 2.0 5.0) :modal-max 4)
)"#;

#[test]
#[allow(clippy::too_many_lines)]
fn binary_runs_three_assemblies_explains_and_resumes_interrupted_run() {
    let dir = scratch_dir("assemblies");
    let steel = write_mat(&dir, "AISI-304", 8000.0, 500.0, 193.0e9);
    let alum = write_mat(&dir, "AA6061", 2700.0, 896.0, 68.9e9);
    let tim = write_tim(&dir, 45.0);
    let (sp, tp, op) = (
        dir.join("slider.fsim"),
        dir.join("twobody.fsim"),
        dir.join("osc.fsim"),
    );
    std::fs::write(&sp, SLIDER_FSIM).unwrap();
    std::fs::write(&tp, TWO_BODY_FSIM).unwrap();
    std::fs::write(&op, FORCED_OSC_FSIM).unwrap();
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

    // 3. Assembly 2 (two-body conduction)
    let out2 = bin_run(&[
        "--json",
        "run",
        tp.to_str().unwrap(),
        lstr,
        "--materials",
        steel_s,
        "--materials",
        alum.to_str().unwrap(),
        "--interfaces",
        tim.to_str().unwrap(),
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
    let steel = write_mat(&dir, "AISI-304", 8000.0, 500.0, 193.0e9);
    let sp = dir.join("slider.fsim");
    std::fs::write(&sp, SLIDER_FSIM).unwrap();
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
            "(study \"slider_via_fs\" (seed 0x7) (versions (constellation :lock \"2026-07\")) (budget (wall 60s) (mem 64MiB)) (capability :cores 1 :mem 64MiB :wall 60s :ops (assembly.*)) (assembly.run \"slider.fsim\" :materials (\"AISI-304.fsmcdpk\") :analysis \"transient\" :hash \"{exact_hash}\"))\n"
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
            "(study \"bad_hash_study\" (seed 0x7) (versions (constellation :lock \"2026-07\")) (budget (wall 60s) (mem 64MiB)) (capability :cores 1 :mem 64MiB :wall 60s :ops (assembly.*)) (assembly.run \"slider.fsim\" :materials (\"AISI-304.fsmcdpk\") :hash \"{}\"))\n",
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
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let fsim = root.join("examples/heatsink-fan/heatsink-fan.fsim");
    let stl = root.join("examples/heatsink-fan/heatsink.stl");
    let pack = root.join("data/reference-project/aa6061.fsmcdpk");

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
    let sp = dir.join("slider.fsim");
    std::fs::write(&sp, SLIDER_FSIM).unwrap();
    let hs =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/heatsink-fan/heatsink-fan.fsim");
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
    let steel = write_mat(&dir, "AISI-304", 8000.0, 500.0, 193.0e9);
    let sp = dir.join("slider.fsim");
    std::fs::write(&sp, SLIDER_FSIM).unwrap();
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
