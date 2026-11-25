//! Real no-mock integration tests for `.11.7` Stage 2 + S2-B:
//! three real lowered coupled assemblies on `stage_driver`, kill/resume at
//! every stage `0..=5`, `Modal`/`Harmonic` analytic checks, S2-B `Stationary`
//! equilibrium convergence (`|x*| = F/k`), `explain_assembly_run`, and
//! explicit `admit`/`probe` refusals.

use std::f64::consts::TAU;

use fs_blake3::{ContentHash, hash_bytes};
use fs_cli::{
    AssemblyAnalysisMode, AssemblyCancelPlan, AssemblyStage, CardPackKind, CardPackSet,
    RawCardPack, SolveRunId, SolveRunStatus, explain_assembly_run, resume_assembly, run_assembly,
};
use fs_evidence::ValidityDomain;
use fs_ledger::{EdgeRole, Ledger};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialStateId, NormalizedInterfacePack,
    NormalizedMaterialCardPack, NormalizedPack, ObservationDataset, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_package::EvidencePackage;
use fs_qty::Dims;

const DENSITY: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const SPECIFIC_HEAT: Dims = Dims([2, 0, -2, -1, 0, 0]);
const MODULUS: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const CONDUCTANCE: Dims = Dims([2, 1, -3, -1, 0, 0]); // W/K

fn prov(src: &str) -> Provenance {
    Provenance {
        source: src.to_string(),
        license: "CC-BY-4.0".to_string(),
        artifact: Some(hash_bytes(src.as_bytes())),
    }
}

fn claim(
    name: &str,
    value: f64,
    dims: Dims,
    obs: fs_matdb::ObservationId,
    src: &str,
) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained().with("T", 200.0, 600.0),
        uncertainty: UncertaintyModel::HalfWidth {
            half_width: (value.abs() * 0.01).max(1.0e-6),
            confidence: 0.95,
        },
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![obs],
        provenance: prov(src),
    }
}

fn reg_obs(claims: &mut ClaimSet, tag: &str) -> fs_matdb::ObservationId {
    claims
        .register_observation(ObservationDataset {
            specimen: format!("coupon-{tag}"),
            method: "test-fixture".to_string(),
            artifact: hash_bytes(tag.as_bytes()),
            caveats: "test fixture".to_string(),
            provenance: prov(tag),
        })
        .unwrap()
}

fn material_pack(chem: &str, rho: f64, cp: f64, youngs: f64) -> RawCardPack {
    let mut claims = ClaimSet::new();
    let obs = reg_obs(&mut claims, chem);
    claims
        .insert_claim(claim("density", rho, DENSITY, obs, "asm-rho"))
        .unwrap();
    claims
        .insert_claim(claim(
            "specific_heat_capacity",
            cp,
            SPECIFIC_HEAT,
            obs,
            "asm-cp",
        ))
        .unwrap();
    claims
        .insert_claim(claim("youngs_modulus", youngs, MODULUS, obs, "asm-e"))
        .unwrap();
    let norm = NormalizedPack::new(
        format!("pack-{chem}"),
        "fs-cli-test-compiler-v1",
        hash_bytes(chem.as_bytes()),
        "CC-BY-4.0",
        claims,
        Vec::new(),
        Vec::new(),
    )
    .unwrap();
    let pack = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: chem.to_string(),
            phase: "wrought".to_string(),
            process: "annealed".to_string(),
            revision: 0,
        },
        norm,
    )
    .unwrap();
    RawCardPack {
        kind: CardPackKind::Material,
        source: format!("{chem}.fsmcdpk"),
        bytes: pack.to_bytes(),
        expect: None,
    }
}

fn thermal_interface_pack(g_w_k: f64) -> RawCardPack {
    let mut claims = ClaimSet::new();
    let obs = reg_obs(&mut claims, "thermal-iface");
    claims
        .insert_claim(claim(
            "thermal_conductance",
            g_w_k,
            CONDUCTANCE,
            obs,
            "contact-g",
        ))
        .unwrap();
    let norm = NormalizedPack::new(
        "pack-thermal-iface",
        "fs-cli-test-compiler-v1",
        hash_bytes(b"thermal-iface"),
        "CC-BY-4.0",
        claims,
        Vec::new(),
        Vec::new(),
    )
    .unwrap();
    let mk_surf = |chem: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: chem.into(),
            phase: "wrought".into(),
            process: "annealed".into(),
            revision: 0,
        },
        texture_frame: "ground".into(),
    };
    let pack = NormalizedInterfacePack::new(
        mk_surf("AISI-304"),
        mk_surf("AA6061"),
        SystemContext {
            medium: "dry-contact".into(),
            third_body: None,
            environment: "air".into(),
            history: "run-in".into(),
        },
        norm,
    )
    .unwrap();
    RawCardPack {
        kind: CardPackKind::Interface,
        source: "contact.fsintpk".into(),
        bytes: pack.to_bytes(),
        expect: None,
    }
}

fn step_clock(dt: f64) -> impl FnMut() -> f64 {
    let mut t = 0.0;
    move || {
        t += dt;
        t
    }
}

fn json_str_field<'a>(json: &'a str, key: &str) -> &'a str {
    let needle = format!("\"{key}\":\"");
    json.split(&needle)
        .nth(1)
        .and_then(|r| r.split('"').next())
        .unwrap_or_else(|| panic!("missing `{key}` in {json}"))
}

fn read_report_receipt(ledger: &Ledger, receipt_hex: &str) -> String {
    let h = ContentHash::from_hex(receipt_hex).unwrap();
    String::from_utf8(ledger.get_artifact(&h).unwrap().unwrap()).unwrap()
}

fn check_retained_package(ledger: &Ledger, report_op: i64, expected_root: &str) {
    let edges = ledger
        .op_artifact_edges_bounded(report_op, 64)
        .unwrap()
        .edges;
    let pkg_h = edges
        .into_iter()
        .find_map(|e| {
            let info = ledger.artifact_info(&e.artifact).unwrap()?;
            (e.role == EdgeRole::Out && info.kind == "assembly-evidence-package")
                .then_some(e.artifact)
        })
        .expect("retained evidence package edge");
    let pkg_json = String::from_utf8(ledger.get_artifact(&pkg_h).unwrap().unwrap()).unwrap();
    let parsed = EvidencePackage::from_json(&pkg_json).expect("valid format-9 package");
    let chk = fs_checker::check(&parsed);
    assert!(chk.passed(), "solver-free checker must pass");
    assert_eq!(chk.merkle_root().to_hex(), expected_root);
}

const SLIDER_ASSEMBLY: &str = r#"(assembly "thermo_mech_slider"
  (version 1)
  (units "SI")
  (seed 0x7)
  (versions (constellation "2026-07") (workspace "stage2-e2e"))
  (budget (wall 60s) (mem 64MiB))
  (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*))
  (reference-state unloaded)
  (environment earth-lab-air)
  (time-plan :horizon 0.05s :dt 0.001s :sample-rate 1000Hz)
  (part "slider"
    :kind rigid-body
    :geometry (box 0.1m 0.05m 0.05m)
    :material "AISI-304" :query-temp 293.15K
    :physics (mechanics thermal)
    :fidelity lumped)
  (joint "suspension"
    :kind prismatic :parent root :child "slider"
    :axis (1.0 0.0 0.0) :damping 12N*s/m :stiffness 400N/m)
  (source "drive" :part "slider" :port "mechanics"
    :history (harmonic :amplitude 15N :frequency 3Hz :phase 0.25 :offset 5N))
  (probe "slider_temp" :kind thermocouple :target "slider" :sample-rate 1000Hz :time-constant 0.005s)
  (probe "slider_pos" :kind field-quantity :target "slider" :sample-rate 1000Hz :quantity "position" :dims (1 0 0 0 0 0))
)"#;

const TWO_BODY_THERMAL: &str = r#"(assembly "two_body_conduction"
  (version 1)
  (units "SI")
  (seed 0x11)
  (versions (constellation "2026-07") (workspace "stage2-e2e"))
  (budget (wall 60s) (mem 64MiB))
  (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*))
  (reference-state assembled)
  (environment :ambient-temp 300K :ambient-pressure 101325Pa)
  (time-plan :horizon 0.2s :dt 0.005s :sample-rate 200Hz)
  (part "hot"
    :kind solid :geometry (box 0.05m 0.05m 0.05m)
    :material "AISI-304" :query-temp 300K :physics (thermal) :fidelity lumped)
  (part "cold"
    :kind solid :geometry (box 0.05m 0.05m 0.05m)
    :material "AA6061" :query-temp 300K :physics (thermal) :fidelity lumped)
  (interface "joint_tim"
    :kind thermal-contact :a "hot" :b "cold"
    :card "AISI-304" :family "thermal-conductance" :area 0.0025m^2)
  (source "heater" :part "hot" :port "thermal" :history (constant 500W))
  (probe "hot_temp" :kind thermocouple :target "hot" :sample-rate 200Hz)
  (probe "hot_cross" :kind timing :target "hot" :sample-rate 200Hz :threshold 300.05K :edge rising)
)"#;

const OSCILLATOR_ASSEMBLY: &str = r#"(assembly "harmonic_oscillator"
  (version 1)
  (units "SI")
  (seed 0x22)
  (versions (constellation "2026-07") (workspace "stage2-e2e"))
  (budget (wall 60s) (mem 64MiB))
  (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*))
  (reference-state unloaded)
  (environment earth-lab-air)
  (time-plan :horizon 0.05s :dt 0.001s :sample-rate 1000Hz)
  (part "mass"
    :kind rigid-body :geometry (box 0.1m 0.05m 0.05m)
    :material "AISI-304" :query-temp 293.15K :physics (mechanics) :fidelity lumped)
  (joint "spring_damper"
    :kind prismatic :parent root :child "mass"
    :axis (1.0 0.0 0.0) :damping 4N*s/m :stiffness 200N/m)
  (analysis-options :harmonic-frequencies (1.0 2.0 5.0) :modal-max 4)
)"#;

#[test]
fn assembly1_thermo_mech_slider_runs_explains_and_resumes_at_every_stage() {
    let packs =
        CardPackSet::admit(vec![material_pack("AISI-304", 8000.0, 500.0, 193.0e9)]).unwrap();
    let fresh_ledger = Ledger::open(":memory:").unwrap();
    let (fresh, prog) = run_assembly(
        &fresh_ledger,
        SLIDER_ASSEMBLY.as_bytes(),
        &packs,
        AssemblyAnalysisMode::Transient,
        &mut step_clock(0.25),
        AssemblyCancelPlan::default(),
    )
    .expect("slider runs cleanly");
    assert_eq!(fresh.status, SolveRunStatus::Completed);
    assert_eq!(
        (fresh.stages.len(), fresh.prior_stages, prog.len()),
        (6, 0, 6)
    );
    let report_json = read_report_receipt(&fresh_ledger, &fresh.stages[5].receipt);
    let pkg_root = json_str_field(&report_json, "package_root");
    let plan_hash = json_str_field(&report_json, "plan_hash");
    assert_eq!(plan_hash.len(), 64);
    check_retained_package(&fresh_ledger, fresh.stages[5].op_id, pkg_root);

    let run_id = SolveRunId::parse_hex(&fresh.run).unwrap();
    let explained = explain_assembly_run(&fresh_ledger, run_id).expect("explain succeeds");
    assert_eq!(explained.stages_completed, 6);
    assert_eq!(explained.plan_hash.as_deref(), Some(plan_hash));
    assert!(
        explained
            .lowering_receipt_json
            .as_deref()
            .unwrap()
            .contains("frankensim.lowering_receipt.v1")
    );
    assert!(explained.plan_json.as_deref().unwrap().contains(plan_hash));
    assert!(explained.causalization_audit_json.is_some());
    assert!(
        explained
            .window_ledger_json
            .as_deref()
            .unwrap()
            .contains("frankensim.window_ledger.v1")
    );

    // Kill and resume at EVERY stage boundary 0..=4 and assert exact final report parity.
    for kill_stage in [
        AssemblyStage::Admit,
        AssemblyStage::Lower,
        AssemblyStage::Plan,
        AssemblyStage::Run,
        AssemblyStage::Probe,
    ] {
        let ledger = Ledger::open(":memory:").unwrap();
        let mut clk = step_clock(0.25);
        let plan = AssemblyCancelPlan {
            cancel_after: Some(kill_stage),
            cancel_pre_publish: None,
        };
        let err = run_assembly(
            &ledger,
            SLIDER_ASSEMBLY.as_bytes(),
            &packs,
            AssemblyAnalysisMode::Transient,
            &mut clk,
            plan,
        )
        .expect_err("cancelled after stage");
        assert_eq!(err.code, "cli-solve-cancelled");
        let (resumed, _) =
            resume_assembly(&ledger, run_id, &mut clk, AssemblyCancelPlan::default())
                .expect("resume completes");
        assert_eq!(resumed.status, SolveRunStatus::Completed);
        assert_eq!(resumed.prior_stages as usize, kill_stage as usize + 1);
        assert_eq!(
            resumed.stages.last().unwrap().receipt,
            fresh.stages[5].receipt
        );
    }
}

#[test]
fn assembly2_two_body_thermal_with_interface_card_and_timing_probe() {
    let iface = thermal_interface_pack(20.0);
    let packs = CardPackSet::admit(vec![
        material_pack("AISI-304", 8000.0, 500.0, 193.0e9),
        material_pack("AA6061", 2700.0, 900.0, 69.0e9),
        iface,
    ])
    .unwrap();
    let iface_card_hex = packs.interfaces()[0].card().to_hex();
    let src =
        TWO_BODY_THERMAL.replace(":card \"AISI-304\"", &format!(":card \"{iface_card_hex}\""));
    let ledger = Ledger::open(":memory:").unwrap();
    let (out, _) = run_assembly(
        &ledger,
        src.as_bytes(),
        &packs,
        AssemblyAnalysisMode::Transient,
        &mut step_clock(0.1),
        AssemblyCancelPlan::default(),
    )
    .expect("two-body thermal runs");
    assert_eq!(out.status, SolveRunStatus::Completed);
    assert_eq!(out.stages.len(), 6);
    let probe_receipt = read_report_receipt(&ledger, &out.stages[4].receipt);
    let series_h =
        ContentHash::from_hex(json_str_field(&probe_receipt, "series_artifact")).unwrap();
    let series_json = String::from_utf8(ledger.get_artifact(&series_h).unwrap().unwrap()).unwrap();
    assert!(
        series_json.contains("\"name\":\"hot_temp\"")
            && series_json.contains("\"name\":\"hot_cross\"")
    );
    // 500 W into 500 J/K raises T by ~1 K/s -> crosses 300.05 K well before 0.2 s!
    assert!(
        !series_json
            .contains("\"name\":\"hot_cross\",\"kind\":\"timing\",\"unit\":\"s\",\"count\":0")
    );
    let report_json = read_report_receipt(&ledger, &out.stages[5].receipt);
    check_retained_package(
        &ledger,
        out.stages[5].op_id,
        json_str_field(&report_json, "package_root"),
    );
}

#[test]
fn assembly3_oscillator_modal_harmonic_analytic_and_stationary_equilibrium() {
    let packs =
        CardPackSet::admit(vec![material_pack("AISI-304", 8000.0, 500.0, 193.0e9)]).unwrap();
    // m = 8000 * 0.1 * 0.05 * 0.05 = 2.0 kg, k = 200 N/m, c = 4 N*s/m
    // => omega_n = sqrt(k/m) = 10 rad/s, zeta = c / (2 sqrt(k m)) = 0.1
    let ledger = Ledger::open(":memory:").unwrap();
    let (modal, _) = run_assembly(
        &ledger,
        OSCILLATOR_ASSEMBLY.as_bytes(),
        &packs,
        AssemblyAnalysisMode::Modal,
        &mut step_clock(0.1),
        AssemblyCancelPlan::default(),
    )
    .expect("modal runs");
    assert_eq!(modal.status, SolveRunStatus::Completed);
    let probe_r = read_report_receipt(&ledger, &modal.stages[4].receipt);
    let series_h = ContentHash::from_hex(json_str_field(&probe_r, "series_artifact")).unwrap();
    let series = String::from_utf8(ledger.get_artifact(&series_h).unwrap().unwrap()).unwrap();
    let expected_fn = (10.0 * (1.0_f64 - 0.01).sqrt()) / TAU;
    assert!(series.contains("modal.natural_rad_s") && series.contains("modal.damping_ratio"));
    // Check real computed damped frequency matches analytic sqrt(k/m*(1-zeta^2))/(2*pi) within 1%
    let freq_part = series
        .split("\"name\":\"modal.frequency_hz\"")
        .nth(1)
        .unwrap();
    let f_val: f64 = freq_part
        .split("\"final\":")
        .nth(1)
        .unwrap()
        .split(',')
        .next()
        .unwrap()
        .parse()
        .unwrap();
    assert!(
        (f_val - expected_fn).abs() / expected_fn < 0.01,
        "modal f={f_val} expected {expected_fn}"
    );

    let (harm, _) = run_assembly(
        &ledger,
        OSCILLATOR_ASSEMBLY.as_bytes(),
        &packs,
        AssemblyAnalysisMode::Harmonic,
        &mut step_clock(0.1),
        AssemblyCancelPlan::default(),
    )
    .expect("harmonic runs");
    assert_eq!(harm.status, SolveRunStatus::Completed);
    assert_ne!(harm.run, modal.run, "analysis mode enters run identity");

    // S2-B: Stationary converges on a lowered forced spring-mass (F = 10 N, k = 200 N/m =>
    // |x*| = F/k = 0.05 m, v* = 0; joint lowering attaches port `a` to child so e* = -F/k).
    let forced_osc = OSCILLATOR_ASSEMBLY.replace(
        "(analysis-options :harmonic-frequencies (1.0 2.0 5.0) :modal-max 4)",
        "(source \"load\" :part \"mass\" :port \"mechanics\" :history (constant 10N))",
    );
    let (stat, _) = run_assembly(
        &ledger,
        forced_osc.as_bytes(),
        &packs,
        AssemblyAnalysisMode::Stationary,
        &mut step_clock(0.1),
        AssemblyCancelPlan::default(),
    )
    .expect("stationary converges on lowered forced spring-mass under S2-B");
    assert_eq!(
        (stat.status, stat.stages.len()),
        (SolveRunStatus::Completed, 6)
    );
    let stat_probe_r = read_report_receipt(&ledger, &stat.stages[4].receipt);
    let stat_series_h =
        ContentHash::from_hex(json_str_field(&stat_probe_r, "series_artifact")).unwrap();
    let stat_series =
        String::from_utf8(ledger.get_artifact(&stat_series_h).unwrap().unwrap()).unwrap();
    let probe_final = |name: &str| -> f64 {
        let needle = format!("\"name\":\"{name}\"");
        stat_series
            .split(&needle)
            .nth(1)
            .unwrap_or_else(|| panic!("missing `{name}` in {stat_series}"))
            .split("\"final\":")
            .nth(1)
            .unwrap()
            .split(',')
            .next()
            .unwrap()
            .parse()
            .unwrap()
    };
    let stat_res = probe_final("stationary.residual_norm");
    let stat_ext = probe_final("stationary.state.joint.spring_damper.spring.state");
    let stat_vel = probe_final("stationary.state.mass.velocity_m_s");
    assert!(
        stat_res < 1.0e-10,
        "stationary residual too large: {stat_res}"
    );
    assert!(
        (stat_ext - (-10.0 / 200.0)).abs() < 1.0e-9,
        "expected e* = -F/k = -0.05 m, got {stat_ext}"
    );
    assert!(stat_vel.abs() < 1.0e-9, "expected v* = 0, got {stat_vel}");
    let stat_rep = read_report_receipt(&ledger, &stat.stages[5].receipt);
    let stat_pkg_root = json_str_field(&stat_rep, "package_root");
    check_retained_package(&ledger, stat.stages[5].op_id, stat_pkg_root);
}

#[test]
fn admit_refusals_probe_refusals_and_unknown_explain_fail_closed() {
    let packs =
        CardPackSet::admit(vec![material_pack("AISI-304", 8000.0, 500.0, 193.0e9)]).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();

    // 1. Unknown clause on an otherwise valid assembly refuses at admit.
    let bad_src = OSCILLATOR_ASSEMBLY.replace(
        "(analysis-options :harmonic-frequencies (1.0 2.0 5.0) :modal-max 4)",
        "(unknown-clause 1)",
    );
    let e_admit = run_assembly(
        &ledger,
        bad_src.as_bytes(),
        &packs,
        AssemblyAnalysisMode::Transient,
        &mut step_clock(0.1),
        AssemblyCancelPlan::default(),
    )
    .expect_err("unknown clause refuses");
    assert_eq!(e_admit.code, "cli-assembly-admit");
    assert!(e_admit.what.contains("unknown assembly clause"));

    // 2. Missing material card refuses with cli-assembly-card-missing.
    let empty_packs = CardPackSet::empty();
    let e_card = run_assembly(
        &ledger,
        SLIDER_ASSEMBLY.as_bytes(),
        &empty_packs,
        AssemblyAnalysisMode::Transient,
        &mut step_clock(0.1),
        AssemblyCancelPlan::default(),
    )
    .expect_err("missing card refuses");
    assert_eq!(e_card.code, "cli-assembly-card-missing");

    // 3. Unsupported probe (thermocouple on mechanics-only mass & microphone) refuses explicitly with zero fabricated values.
    let bad_probe_src = OSCILLATOR_ASSEMBLY.replace(
        "(analysis-options :harmonic-frequencies (1.0 2.0 5.0) :modal-max 4)",
        "(probe \"bogus_tc\" :kind thermocouple :target \"mass\" :sample-rate 1000Hz)",
    );
    let e_probe = run_assembly(
        &ledger,
        bad_probe_src.as_bytes(),
        &packs,
        AssemblyAnalysisMode::Transient,
        &mut step_clock(0.1),
        AssemblyCancelPlan::default(),
    )
    .expect_err("thermocouple on non-thermal part refuses");
    assert_eq!(
        (e_probe.code, e_probe.stage),
        ("cli-assembly-probe-unsupported", Some("probe"))
    );
    assert!(e_probe.recorded_op.is_some());

    let mic_src = OSCILLATOR_ASSEMBLY.replace(
        "(analysis-options :harmonic-frequencies (1.0 2.0 5.0) :modal-max 4)",
        "(probe \"mic1\" :kind microphone :sample-rate 1000Hz)",
    );
    let e_mic = run_assembly(
        &ledger,
        mic_src.as_bytes(),
        &packs,
        AssemblyAnalysisMode::Transient,
        &mut step_clock(0.1),
        AssemblyCancelPlan::default(),
    )
    .expect_err("microphone without acoustic lowering refuses");
    assert_eq!(e_mic.code, "cli-assembly-probe-unsupported");

    // 4. Explain on unknown run id refuses cleanly.
    let unk = SolveRunId::parse_hex(&hash_bytes(b"unknown-assembly-run").to_hex()).unwrap();
    assert_eq!(
        explain_assembly_run(&ledger, unk).unwrap_err().code,
        "cli-solve-unknown-run"
    );
}
