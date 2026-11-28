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
const CONDUCTIVITY: Dims = Dims([1, 1, -3, -1, 0, 0]); // W/(m K)

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

fn pack_from_claims(chem: &str, claims: ClaimSet) -> RawCardPack {
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
    pack_from_claims(chem, claims)
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
    let sides = ["AISI-304", "AA6061"];
    interface_pack("thermal-iface", sides, "dry-contact", claims)
}

/// An interface pack `pack-<tag>` of `claims` between two annealed wrought surfaces.
fn interface_pack(tag: &str, sides: [&str; 2], medium: &str, claims: ClaimSet) -> RawCardPack {
    let norm = NormalizedPack::new(
        format!("pack-{tag}"),
        "fs-cli-test-compiler-v1",
        hash_bytes(tag.as_bytes()),
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
        mk_surf(sides[0]),
        mk_surf(sides[1]),
        SystemContext {
            medium: medium.into(),
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

const HEATED_BODY: &str = r#"(assembly "horizon_heater"
  (version 1)
  (units "SI")
  (seed 0x31)
  (versions (constellation "2026-07") (workspace "horizon-regression"))
  (budget (wall 60s) (mem 64MiB))
  (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*))
  (reference-state assembled)
  (environment :ambient-temp 300K :ambient-pressure 101325Pa)
  (time-plan :horizon 0.01025s :dt 0.003s :sample-rate 333.3333333333333Hz)
  (part "body"
    :kind solid :geometry (box 0.1m 0.1m 0.1m)
    :material "fixture" :query-temp 300K :physics (thermal) :fidelity lumped)
  (source "heater" :part "body" :port "thermal" :history (constant 500W))
  (probe "temperature" :kind thermocouple :target "body" :sample-rate 333.3333333333333Hz)
)"#;

#[test]
fn transient_horizon_uses_accepted_time_and_caps_the_final_step() {
    // rho * volume * cp = 500 J/K, so 500 W gives exactly 1 K/s.
    let packs = CardPackSet::admit(vec![material_pack("fixture", 1000.0, 500.0, 1.0e9)]).unwrap();
    for (horizon, history) in [
        (0.01025, "(constant 500W)"), // rounding the count used to stop early
        (0.011, "(constant 500W)"),   // rounding up used to overshoot
        (0.01025, "(step :t 0.004s :before 500W :after 500W)"),
    ] {
        // Equal levels isolate source-knot negotiation from the heat law.
        let src = HEATED_BODY
            .replace("0.01025s", &format!("{horizon}s"))
            .replace("(constant 500W)", history);
        let ledger = Ledger::open(":memory:").unwrap();
        let (out, _) = run_assembly(
            &ledger,
            src.as_bytes(),
            &packs,
            AssemblyAnalysisMode::Transient,
            &mut step_clock(0.1),
            AssemblyCancelPlan::default(),
        )
        .expect("complete requested horizon");
        assert_eq!(out.status, SolveRunStatus::Completed);
        let receipt = read_report_receipt(&ledger, &out.stages[3].receipt);
        let trace = read_report_receipt(&ledger, json_str_field(&receipt, "trace_artifact"));
        let times: Vec<f64> = trace
            .split("\"times_s\":[")
            .nth(1)
            .unwrap()
            .split(']')
            .next()
            .unwrap()
            .split(',')
            .map(|v| v.parse().unwrap())
            .collect();
        assert_eq!(times[0], 0.0);
        assert_eq!(*times.last().unwrap(), horizon, "retained times {times:?}");
        assert!(
            times
                .windows(2)
                .all(|w| { w[1] > w[0] && w[1] - w[0] <= 0.003 * (1.0 + 8.0 * f64::EPSILON) })
        );
        if history.starts_with("(step") {
            assert!(times.len() > 5, "source knots must shorten accepted steps");
        } else {
            assert_eq!(times.len(), 5);
        }
        let final_temperature: f64 = trace
            .rsplit("\"temperature_k\":")
            .next()
            .unwrap()
            .split('}')
            .next()
            .unwrap()
            .parse()
            .unwrap();
        assert!((final_temperature - (300.0 + horizon)).abs() < 1.0e-11);
    }
}

/// Bead vr9rw: a transient probe is read at its declared `:sample-rate`, on the time plan's step.
#[test]
fn a_probe_reads_at_its_declared_rate_and_refuses_one_the_steps_cannot_hit() {
    let packs =
        CardPackSet::admit(vec![material_pack("AISI-304", 8000.0, 500.0, 193.0e9)]).unwrap();
    let extra = r#"  (probe "slow_temp" :kind thermocouple :target "slider" :sample-rate 200Hz :time-constant 0.005s)
  (probe "slow_pos" :kind field-quantity :target "slider" :sample-rate 200Hz :quantity "position" :dims (1 0 0 0 0 0))
"#;
    let marker = "  (probe \"slider_temp\"";
    let src = SLIDER_ASSEMBLY.replace(marker, &format!("{extra}{marker}"));
    assert_ne!(src, SLIDER_ASSEMBLY);
    let ledger = Ledger::open(":memory:").unwrap();
    let (out, _) = run_assembly(
        &ledger,
        src.as_bytes(),
        &packs,
        AssemblyAnalysisMode::Transient,
        &mut step_clock(0.1),
        AssemblyCancelPlan::default(),
    )
    .expect("200 Hz probes on a 1 ms plan run");
    assert_eq!(out.status, SolveRunStatus::Completed);
    let exp = explain_assembly_run(&ledger, SolveRunId::parse_hex(&out.run).unwrap()).unwrap();
    let pj = exp.probe_series_json.as_deref();
    // 0.05 s at 1 ms: 51 samples at 1000 Hz, every fifth of them (11) at 200 Hz, lag included.
    for (fast, slow) in [("slider_temp", "slow_temp"), ("slider_pos", "slow_pos")] {
        let (f, s) = (
            extract_probe_values(pj, fast),
            extract_probe_values(pj, slow),
        );
        assert_eq!((f.len(), s.len()), (51, 11), "{fast}");
        let every_fifth: Vec<f64> = f.iter().step_by(5).copied().collect();
        assert_eq!(s, every_fifth, "{slow}");
    }
    for (rate, words) in [
        ("3000Hz", "faster than the time plan's step rate"),
        ("300Hz", "not a whole number"),
    ] {
        let bad = SLIDER_ASSEMBLY.replace(
            "(probe \"slider_pos\" :kind field-quantity :target \"slider\" :sample-rate 1000Hz",
            &format!(
                "(probe \"slider_pos\" :kind field-quantity :target \"slider\" :sample-rate {rate}"
            ),
        );
        assert_ne!(bad, SLIDER_ASSEMBLY);
        let e = run_assembly(
            &Ledger::open(":memory:").unwrap(),
            bad.as_bytes(),
            &packs,
            AssemblyAnalysisMode::Transient,
            &mut step_clock(0.1),
            AssemblyCancelPlan::default(),
        )
        .expect_err("a rate the steps cannot hit refuses");
        assert_eq!(
            (e.code, e.stage),
            ("cli-assembly-probe-rate", Some("probe"))
        );
        assert!(
            e.what.contains("`slider_pos`") && e.what.contains(words),
            "{}",
            e.what
        );
    }
    // A source knot at 4 ms moves the 3 ms steps after it off the 6 ms grid of a 1/(2 dt) probe.
    let fixture = CardPackSet::admit(vec![material_pack("fixture", 1000.0, 500.0, 1.0e9)]).unwrap();
    let knotted = HEATED_BODY
        .replace(
            "(constant 500W)",
            "(step :t 0.004s :before 500W :after 500W)",
        )
        .replace(
            ":target \"body\" :sample-rate 333.3333333333333Hz",
            ":target \"body\" :sample-rate 166.66666666666666Hz",
        );
    assert!(knotted.contains("(step :t 0.004s") && knotted.contains("166.66666666666666Hz"));
    let e = run_assembly(
        &Ledger::open(":memory:").unwrap(),
        knotted.as_bytes(),
        &fixture,
        AssemblyAnalysisMode::Transient,
        &mut step_clock(0.1),
        AssemblyCancelPlan::default(),
    )
    .expect_err("a knot-shifted grid refuses");
    assert_eq!(e.code, "cli-assembly-probe-rate");
    assert!(e.what.contains("is not a step of this run"), "{}", e.what);
}

#[test]
fn transient_step_limit_refuses_instead_of_publishing_a_truncated_run() {
    let packs = CardPackSet::admit(vec![material_pack("fixture", 1000.0, 500.0, 1.0e9)]).unwrap();
    let src = HEATED_BODY.replace(":dt 0.003s", ":dt 0.00000001s");
    let ledger = Ledger::open(":memory:").unwrap();
    let refusal = run_assembly(
        &ledger,
        src.as_bytes(),
        &packs,
        AssemblyAnalysisMode::Transient,
        &mut step_clock(0.1),
        AssemblyCancelPlan::default(),
    )
    .expect_err("the whole horizon exceeds the accepted-step limit");
    assert_eq!(refusal.code, "cli-assembly-step-limit");
}

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
    // The ledger's lowering receipt is the one the lower stage's receipt names.
    let lower = &fresh.stages[1];
    assert_eq!(lower.stage, "lower");
    let stage_json = read_report_receipt(&fresh_ledger, &lower.receipt);
    let receipt = explained.lowering_receipt_json.as_deref().unwrap();
    let artifact = hash_bytes(receipt.as_bytes()).to_hex();
    assert_eq!(json_str_field(&stage_json, "lowering_artifact"), artifact);
    let receipt_hash = json_str_field(receipt, "receipt_hash");
    let named = json_str_field(&stage_json, "lowering_receipt_hash");
    assert_eq!(named, receipt_hash);
    let explained_hash = explained.lowering_receipt_hash.as_deref();
    assert_eq!(explained_hash, Some(receipt_hash));
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
        "(probe \"mic1\" :kind microphone :target \"mass\" :sample-rate 1000Hz)",
    );
    let e_mic = run_assembly(
        &ledger,
        mic_src.as_bytes(),
        &packs,
        AssemblyAnalysisMode::Transient,
        &mut step_clock(0.1),
        AssemblyCancelPlan::default(),
    )
    .expect_err("microphone targeting non-radiating mass refuses");
    assert_eq!(e_mic.code, "cli-assembly-lower");

    // 4. Explain on unknown run id refuses cleanly.
    let unk = SolveRunId::parse_hex(&hash_bytes(b"unknown-assembly-run").to_hex()).unwrap();
    assert_eq!(
        explain_assembly_run(&ledger, unk).unwrap_err().code,
        "cli-solve-unknown-run"
    );
}

#[test]
fn qoi_budget_measures_stated_card_corners_resumes_without_recompute_and_refuses_bad_dims_or_wall()
{
    let packs =
        CardPackSet::admit(vec![material_pack("AISI-304", 8000.0, 500.0, 193.0e9)]).unwrap();
    let req_src = SLIDER_ASSEMBLY.replace(
        "\n)",
        "\n  (requirement \"t_lim\" :probe \"slider_temp\" :stat max :relation at-most :limit 310K)\n)",
    );
    let ledger = Ledger::open(":memory:").unwrap();
    let e_cancel = run_assembly(
        &ledger,
        req_src.as_bytes(),
        &packs,
        AssemblyAnalysisMode::Transient,
        &mut step_clock(0.1),
        AssemblyCancelPlan {
            cancel_after: Some(AssemblyStage::Probe),
            cancel_pre_publish: None,
        },
    )
    .expect_err("cancel after probe");
    assert_eq!(
        (e_cancel.code, e_cancel.stage),
        ("cli-solve-cancelled", Some("report"))
    );
    let run_id = SolveRunId::parse_hex(e_cancel.run.as_deref().unwrap()).unwrap();
    let before_qb = explain_assembly_run(&ledger, run_id)
        .unwrap()
        .qoi_budget_json
        .unwrap();
    let (resumed, _) = resume_assembly(
        &ledger,
        run_id,
        &mut step_clock(0.1),
        AssemblyCancelPlan::default(),
    )
    .unwrap();
    assert_eq!(
        (resumed.status, resumed.prior_stages, resumed.stages.len()),
        (SolveRunStatus::Completed, 5, 1)
    );
    let after_qb = explain_assembly_run(&ledger, run_id)
        .unwrap()
        .qoi_budget_json
        .unwrap();
    assert_eq!(
        before_qb, after_qb,
        "resume restores retained qoi_budget without recomputing"
    );
    for needle in [
        "\"id\":\"t_lim\"",
        "\"outcome\":\"indeterminate\"",
        "\"kind\":\"parameters\",\"state\":\"measured\"",
        "\"kind\":\"discretization\",\"state\":\"measured\"",
        "\"kind\":\"roundoff\",\"state\":\"measured\"",
        "\"label\":\"param_corner_3\"",
        "\"wall_s\":",
        "conservative overlap with solver-algebraic",
    ] {
        assert!(
            after_qb.contains(needle),
            "missing `{needle}` in {after_qb}"
        );
    }
    assert!(!after_qb.contains(
        "\"kind\":\"parameters\",\"state\":\"measured\",\"value\":0.00000000000000000e0"
    ));
    let rep_rcpt = read_report_receipt(&ledger, &resumed.stages[0].receipt);
    let rep_json_h = ContentHash::from_hex(json_str_field(&rep_rcpt, "report_json")).unwrap();
    let rep_json = String::from_utf8(ledger.get_artifact(&rep_json_h).unwrap().unwrap()).unwrap();
    assert!(rep_json.contains("\"id\": \"t_lim\"") && rep_json.contains("QoI uncertainty budget"));

    let bad_dim = SLIDER_ASSEMBLY.replace(
        "\n)",
        "\n  (requirement \"bad_dim\" :probe \"slider_temp\" :stat max :relation at-most :limit 5m)\n)",
    );
    let e_dim = run_assembly(
        &ledger,
        bad_dim.as_bytes(),
        &packs,
        AssemblyAnalysisMode::Transient,
        &mut step_clock(0.1),
        AssemblyCancelPlan::default(),
    )
    .unwrap_err();
    assert_eq!(
        (e_dim.code, e_dim.stage),
        ("cli-assembly-requirement", Some("probe"))
    );

    let tight_wall = req_src
        .replace("(wall 60s)", "(wall 1s)")
        .replace(":wall 60s", ":wall 1s");
    let mut n = 0usize;
    let (w_out, _) = run_assembly(
        &ledger,
        tight_wall.as_bytes(),
        &packs,
        AssemblyAnalysisMode::Transient,
        &mut || {
            n += 1;
            if n <= 7 { 0.0 } else { 0.999_999_999 }
        },
        AssemblyCancelPlan::default(),
    )
    .unwrap();
    let w_qb = explain_assembly_run(&ledger, SolveRunId::parse_hex(&w_out.run).unwrap())
        .unwrap()
        .qoi_budget_json
        .unwrap();
    assert!(
        w_qb.contains("\"status\":\"budget_exceeded\"") && w_qb.contains("wall budget"),
        "{w_qb}"
    );
}

fn run_budget(packs: &CardPackSet, src: &str) -> String {
    let ledger = Ledger::open(":memory:").unwrap();
    let (out, _) = run_assembly(
        &ledger,
        src.as_bytes(),
        packs,
        AssemblyAnalysisMode::Transient,
        &mut step_clock(0.1),
        AssemblyCancelPlan::default(),
    )
    .unwrap();
    explain_assembly_run(&ledger, SolveRunId::parse_hex(&out.run).unwrap())
        .unwrap()
        .qoi_budget_json
        .unwrap()
}

fn claim_pack(
    chem: &str,
    prop: &str,
    value: f64,
    dims: Dims,
    unc: UncertaintyModel,
) -> RawCardPack {
    let mut claims = ClaimSet::new();
    let obs = reg_obs(&mut claims, chem);
    let mut c = claim(prop, value, dims, obs, chem);
    c.uncertainty = unc;
    claims.insert_claim(c).unwrap();
    pack_from_claims(chem, claims)
}

#[test]
fn shear_film_third_body_qoi_budget_handles_unstated_and_stated_viscosity() {
    let visc = Dims([-1, 1, -1, 0, 0, 0]);
    let zero = UncertaintyModel::HalfWidth {
        half_width: 0.0,
        confidence: 0.95,
    };
    let body = claim_pack("steel-exact", "density", 8000.0, DENSITY, zero);
    let src = r#"(assembly "shear_film_budget" (version 1) (units "SI") (seed 0x41)
  (versions (constellation "2026-07") (workspace "0h46c")) (budget (wall 60s) (mem 64MiB))
  (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*))
  (reference-state assembled) (environment earth-lab-air) (time-plan :horizon 0.05s :dt 0.005s :sample-rate 200Hz)
  (part "rotor" :kind rigid-body :geometry (box 0.1m 0.05m 0.05m) :material "steel-exact" :query-temp 300K :physics (mechanics) :fidelity lumped
    :initial (:position (0m 0m 0m) :orientation (1 0 0 0) :linear-velocity (1m/s 0m/s 0m/s) :angular-velocity (0rad/s 0rad/s 0rad/s)))
  (part "pad" :kind rigid-body :geometry (box 0.1m 0.05m 0.05m) :material "steel-exact" :query-temp 300K :physics (mechanics) :fidelity lumped
    :initial (:position (0m 0m 0m) :orientation (1 0 0 0) :linear-velocity (0m/s 0m/s 0m/s) :angular-velocity (0rad/s 0rad/s 0rad/s)))
  (interface "film" :kind shear-film :a "rotor" :b "pad" :fluid "test-oil" :gap 0.001m :fluid-at (T 300K) :area 0.01m^2)
  (probe "rotor_vel" :kind field-quantity :target "rotor" :sample-rate 200Hz :quantity "velocity-x" :dims (1 0 -1 0 0 0))
  (requirement "v_end" :probe "rotor_vel" :stat final :relation at-most :limit 0.95m/s))"#;
    let run_with = |unc| {
        let p = CardPackSet::admit(vec![
            body.clone(),
            claim_pack("test-oil", "dynamic_viscosity", 2.0, visc, unc),
        ])
        .unwrap();
        run_budget(&p, src)
    };
    let u = run_with(UncertaintyModel::Unstated);
    assert!(
        u.contains("\"kind\":\"parameters\",\"state\":\"no-data\"")
            && u.contains("dynamic_viscosity"),
        "{u}"
    );
    let s = run_with(UncertaintyModel::HalfWidth {
        half_width: 0.5,
        confidence: 0.95,
    });
    assert!(
        s.contains("\"kind\":\"parameters\",\"state\":\"measured\"")
            && s.contains("dynamic_viscosity"),
        "{s}"
    );
    assert!(
        !s.contains(
            "\"kind\":\"parameters\",\"state\":\"measured\",\"value\":0.00000000000000000e0"
        ),
        "{s}"
    );
}

#[test]
fn declared_envelopes_yield_determinate_pass_fail_and_indeterminate_on_drop() {
    let packs = CardPackSet::admit(vec![material_pack("fixture", 1000.0, 500.0, 1.0e9)]).unwrap();
    let g_env =
        "\n  (envelope geometry :part \"body\" :axis lx :half-width 0.001m :source \"ISO 2768-m\")";
    let rest = "\n  (envelope boundary-conditions :target \"heater\" :half-width 5W :source \"PSU 1%\")\
                \n  (envelope measurement :probe \"temperature\" :half-width 0.0005K :source \"SPRT cal\")\
                \n  (envelope model-form :probe \"temperature\" :half-width 0.0005K :source \"Biot<0.01\")";
    let case =
        |extra: &str| run_budget(&packs, &HEATED_BODY.replace("\n)", &format!("{extra}\n)")));
    let pass = case(&format!(
        "{g_env}{rest}\n  (requirement \"t_lim\" :probe \"temperature\" :stat final :relation at-most :limit 305K)"
    ));
    for k in [
        "\"outcome\":\"pass\"",
        "\"unmeasured_terms\":[]",
        "\"kind\":\"geometry\",\"state\":\"measured\"",
        "\"kind\":\"boundary-conditions\",\"state\":\"measured\"",
        "\"kind\":\"measurement\",\"state\":\"measured\"",
        "\"kind\":\"model-form\",\"state\":\"measured\"",
        "\"label\":\"geom_corner_1\"",
        "\"label\":\"bc_corner_1\"",
        "omits O(eps^2) cross-category terms",
        "\"wall_s\":",
    ] {
        assert!(pass.contains(k), "missing `{k}` in {pass}");
    }
    assert!(
        !pass.contains(
            "\"kind\":\"geometry\",\"state\":\"measured\",\"value\":0.00000000000000000e0"
        )
    );
    assert!(!pass.contains(
        "\"kind\":\"boundary-conditions\",\"state\":\"measured\",\"value\":0.00000000000000000e0"
    ));
    let fail = case(&format!(
        "{g_env}{rest}\n  (requirement \"t_lim\" :probe \"temperature\" :stat final :relation at-most :limit 300.001K)"
    ));
    assert!(
        fail.contains("\"outcome\":\"fail\"") && fail.contains("\"unmeasured_terms\":[]"),
        "{fail}"
    );
    let ind = case(&format!(
        "{rest}\n  (requirement \"t_lim\" :probe \"temperature\" :stat final :relation at-most :limit 305K)"
    ));
    assert!(
        ind.contains("\"outcome\":\"indeterminate\"")
            && ind.contains("\"unmeasured_terms\":[\"geometry\"]"),
        "{ind}"
    );
}

#[test]
fn srm1453_slab_stationary_matches_closed_form_and_refuses_transient() {
    let b = include_bytes!("../../../examples/assemblies/packs/eps-nist-srm1453-2013.fsmcdpk");
    let raw = RawCardPack {
        kind: CardPackKind::Material,
        source: "srm".into(),
        bytes: b.to_vec(),
        expect: None,
    };
    let p = CardPackSet::admit(vec![raw]).unwrap();
    let l = Ledger::open(":memory:").unwrap();
    let s = include_bytes!("../../../examples/assemblies/srm1453_slab.fsim");
    let (mut c, cp) = (step_clock(0.1), AssemblyCancelPlan::default());
    let mut r = |m| run_assembly(&l, s, &p, m, &mut c, cp);
    let e = r(AssemblyAnalysisMode::Transient).unwrap_err();
    assert_eq!(
        (e.code, e.what.contains("transient")),
        ("cli-assembly-fidelity-analysis", true)
    );
    let (out, _) = r(AssemblyAnalysisMode::Stationary).unwrap();
    let exp = explain_assembly_run(&l, SolveRunId::parse_hex(&out.run).unwrap()).unwrap();
    let ps = exp.probe_series_json.unwrap();
    let t = ps
        .split("\"slab_heat_rate\"")
        .nth(1)
        .unwrap()
        .split("\"final\":")
        .nth(1)
        .unwrap();
    let q: f64 = t.split(',').next().unwrap().parse().unwrap();
    let qr = (0.00111 - 0.0000424 * 41.0 + 0.000115 * 297.15) * (0.25 * 0.25 / 0.01343) * 20.0;
    assert!((q - qr).abs() / qr <= 1.0e-12 && ps.contains("stationary.residual_norm"));
    let qb = exp.qoi_budget_json.unwrap();
    for k in [
        "\"outcome\":\"pass\"",
        "\"unmeasured_terms\":[]",
        "\"kind\":\"discretization\",\"state\":\"negligible\"",
        "\"kind\":\"model-form\",\"state\":\"negligible\"",
        "covered by source",
    ] {
        assert!(qb.contains(k), "missing `{k}` in {qb}");
    }
    let s_str = std::str::from_utf8(s).unwrap();
    let run_src = |src: &str, packs: &CardPackSet, mode| {
        run_assembly(&l, src.as_bytes(), packs, mode, &mut step_clock(0.1), cp)
    };
    let q_of = |a: &str, b: &str| {
        let (o, _) = run_src(&s_str.replace(a, b), &p, AssemblyAnalysisMode::Stationary).unwrap();
        explain_assembly_run(&l, SolveRunId::parse_hex(&o.run).unwrap())
            .unwrap()
            .qoi_budget_json
            .unwrap()
    };
    assert!(q_of("3.35W", "2.95W").contains("\"outcome\":\"fail\""));
    assert!(q_of("3.35W", "3.20W").contains("\"outcome\":\"indeterminate\""));
    let drop_m = q_of(
        "(envelope measurement :probe \"slab_heat_rate\" :half-width 0.005W :source \"declared apparatus tolerance (not measured)\")",
        "",
    );
    assert!(
        drop_m.contains("\"outcome\":\"indeterminate\"")
            && drop_m.contains("\"unmeasured_terms\":[\"measurement\"]")
    );
    let assert_ref = |src: &str, packs: &CardPackSet, mode, needle: &str| {
        let w = run_src(src, packs, mode).unwrap_err().what;
        assert!(w.contains(needle), "expected `{needle}` in `{w}`");
    };
    assert_ref(
        &s_str.replace("0.01343m", "0.0254m"),
        &p,
        AssemblyAnalysisMode::Stationary,
        "[12.86, 13.88] mm",
    );
    assert_ref(
        &s_str.replace("307.15K", "309.15K"),
        &p,
        AssemblyAnalysisMode::Stationary,
        "20 K",
    );
    assert_ref(
        &s_str
            .replace("307.15K", "305.15K")
            .replace("287.15K", "285.15K"),
        &p,
        AssemblyAnalysisMode::Stationary,
        "mean-temperature",
    );
    let p_k = CardPackSet::admit(vec![material_pack("fixture", 1000.0, 500.0, 1.0e9)]).unwrap();
    let tr_cov = HEATED_BODY.replace(
        "\n)",
        "\n  (envelope model-form :probe \"temperature\" :covered-by-source \"NIST\")\n  (requirement \"r\" :probe \"temperature\" :stat final :relation at-most :limit 305K)\n)",
    );
    assert_ref(&tr_cov, &p_k, AssemblyAnalysisMode::Transient, "stationary");
    let p_non = CardPackSet::admit(vec![claim_pack(
        "other-eps",
        "thermal_conductivity",
        0.04,
        Dims([1, 1, -3, -1, 0, 0]),
        UncertaintyModel::RelativeHalfWidth {
            fraction: 0.01,
            confidence: 0.95,
        },
    )])
    .unwrap();
    let non_srm = s_str.replace("eps-nist-srm1453-2013", "other-eps").replace(
        ":query-at (mean-temperature 297.15K bulk-density 41kg/m^3 temperature-difference 20K\n               source-lot-srm1453 1 source-conditioning-protocol-satisfied 1\n               unstacked-unsliced-uncontaminated-board 1)",
        ":query-temp 297.15K",
    );
    assert_ref(
        &non_srm,
        &p_non,
        AssemblyAnalysisMode::Stationary,
        "other-eps",
    );
}

#[test]
fn stage3_mode_budgets_si_refusal_and_modal_param_sweep() {
    use AssemblyAnalysisMode as Am;
    let sb = include_bytes!(
        "../../../examples/assemblies/packs/steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk"
    );
    let raw = RawCardPack {
        kind: CardPackKind::Material,
        source: "steel".into(),
        bytes: sb.to_vec(),
        expect: None,
    };
    let (steel, l) = (
        CardPackSet::admit(vec![raw]).unwrap(),
        Ledger::open(":memory:").unwrap(),
    );
    let run = |s: &[u8], p: &CardPackSet, m| {
        run_assembly(
            &l,
            s,
            p,
            m,
            &mut step_clock(0.1),
            AssemblyCancelPlan::default(),
        )
    };
    let qb_of = |s: &[u8], p: &CardPackSet, m| {
        let id = SolveRunId::parse_hex(&run(s, p, m).unwrap().0.run).unwrap();
        explain_assembly_run(&l, id)
            .unwrap()
            .qoi_budget_json
            .unwrap()
    };
    let fo = include_bytes!("../../../examples/assemblies/forced_oscillator.fsim");
    for (mode, id, unit) in [
        (Am::Modal, "mode1_freq", "Hz"),
        (Am::Harmonic, "harm_peak", "s"),
        (Am::Stationary, "eq_extension", "m"),
    ] {
        let qb = qb_of(fo, &steel, mode);
        let (a, b) = (format!("\"id\":\"{id}\""), format!("\"unit\":\"{unit}\""));
        for k in [
            a.as_str(),
            b.as_str(),
            "\"outcome\":\"indeterminate\"",
            "\"kind\":\"roundoff\",\"state\":\"no-data\"",
            "\"kind\":\"solver-algebraic\",\"state\":\"measured\"",
            "\"kind\":\"discretization\",\"state\":\"measured\"",
        ] {
            assert!(qb.contains(k), "mode {mode:?} missing `{k}` in {qb}");
        }
    }
    let fo_s = std::str::from_utf8(fo).unwrap();
    let si = fo_s.replace(
        "(requirement \"eq_extension\" :analysis stationary :probe \"stationary.state.joint.spring_damper.spring.state\" :stat final :relation at-least :limit -0.06m)",
        "(requirement \"res_lim\" :analysis stationary :probe \"stationary.residual_norm\" :stat final :relation at-most :limit 1)",
    );
    for (src, pat) in [
        (si, "probe `stationary.residual_norm` has no declared unit"),
        (fo_s.replace("spring.state", "missing.state"), "missing"),
    ] {
        let w = run(src.as_bytes(), &steel, Am::Stationary)
            .unwrap_err()
            .what;
        assert!(w.contains(pat), "expected `{pat}` in `{w}`");
    }
    let u = UncertaintyModel::RelativeHalfWidth {
        fraction: 0.02,
        confidence: 0.95,
    };
    let uc = claim_pack("steel-unc", "density", 7800.0, Dims([-3, 1, 0, 0, 0, 0]), u);
    let unc = CardPackSet::admit(vec![uc]).unwrap();
    let fu = fo_s.replace("steel-music-wire-ligo-p1400229", "steel-unc");
    let qbu = qb_of(fu.as_bytes(), &unc, Am::Modal);
    let tail = qbu
        .split("\"kind\":\"parameters\",\"state\":\"measured\",\"value\":")
        .nth(1)
        .unwrap();
    let pv: f64 = tail.split(',').next().unwrap().parse().unwrap();
    let fd = |m: f64| {
        let (wn, z) = ((200.0 / m).sqrt(), 4.0 / (2.0 * (200.0 * m).sqrt()));
        wn * (1.0 - z * z).sqrt() / std::f64::consts::TAU
    };
    let ex = (fd(1.95 * 0.98) - fd(1.95)).abs();
    assert!((pv - ex).abs() / ex < 1.0e-6, "pv={pv} ex={ex}");
}

/// The rdh memory pair's bushing values, labelled as fixture values in each
/// claim's provenance: Kelvin viscosity 0.4 N s/m and one branch of 1.69 N/m
/// and 0.2 s, all on x, 0 on the other axes.
const FIXTURE: &str = "[fixture: rdh relaxation memory_pair, not a measured compound]";

/// A `relaxation-bushing` card of the fixture's values.
fn relaxation_pack() -> RawCardPack {
    let (damping, stiffness) = (Dims([0, 1, -1, 0, 0, 0]), Dims([0, 1, -2, 0, 0, 0]));
    let turning = [Dims([2, 1, -1, 0, 0, 0]), Dims([2, 1, -2, 0, 0, 0])];
    let time = Dims([0, 0, 1, 0, 0, 0]);
    let mut values = vec![("relaxation_time_0".to_owned(), 0.2, time)];
    for (axis, x) in ["x", "y", "z"].into_iter().enumerate() {
        let on_x = if axis == 0 { 1.0 } else { 0.0 };
        values.push((format!("damping_{x}"), 0.4 * on_x, damping));
        values.push((format!("stiffness_0_{x}"), 1.69 * on_x, stiffness));
        values.push((format!("damping_r{x}"), 0.0, turning[0]));
        values.push((format!("stiffness_0_r{x}"), 0.0, turning[1]));
    }
    let mut claims = ClaimSet::new();
    let obs = reg_obs(&mut claims, FIXTURE);
    for (name, value, dims) in values {
        let fixture = claim(&name, value, dims, obs, FIXTURE);
        claims.insert_claim(fixture).unwrap();
    }
    interface_pack("relaxation", ["steel"; 2], "relaxation-bushing", claims)
}

/// A 2 kg steel box at rest at the origin, room for a relaxation network's
/// dense preparation (340023768 multiply-adds a step for this pair).
fn steel_box(name: &str) -> String {
    format!(
        r#"(part "{name}" :kind rigid-body :geometry (box 0.1m 0.05m 0.05m) :material "steel"
    :query-temp 300K :physics (mechanics) :fidelity lumped :max-work-units 1000000000
    :initial (:position (0m 0m 0m) :orientation (1 0 0 0)
      :linear-velocity (0m/s 0m/s 0m/s) :angular-velocity (0rad/s 0rad/s 0rad/s)))"#
    )
}

/// A relaxation bushing `name` from `a` to `b` on `card`, its branch holding
/// the fixture's elastic displacement 0.6 / 1.3 m on x.
fn bushing(name: &str, (a, b): (&str, &str), card: &str) -> String {
    format!(
        r#"(interface "{name}" :kind mechanical-contact :a "{a}" :b "{b}"
    :card "{card}" :family "relaxation-bushing"
    :memory (displacement-0-x 0.4615384615384615m displacement-0-y 0m displacement-0-z 0m
      displacement-0-rx 0rad displacement-0-ry 0rad displacement-0-rz 0rad))"#
    )
}

/// Two boxes on a relaxation bushing `mount` whose `drive_x` holds `port`'s
/// source at 0.7 N, then `extra`.
fn held_bushing(card: &str, port: &str, extra: &str) -> String {
    let (rotor, pad) = (steel_box("rotor"), steel_box("pad"));
    let mount = bushing("mount", ("rotor", "pad"), card);
    format!(
        r#"(assembly "held_bushing" (version 1) (units "SI") (seed 0x43)
  (versions (constellation "2026-07") (workspace "11-2-codec")) (budget (wall 60s) (mem 64MiB))
  (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*))
  (reference-state assembled) (environment earth-lab-air)
  (time-plan :horizon 0.02s :dt 0.005s :sample-rate 200Hz)
  {rotor} {pad} {mount}
  (source "hold" :interface "mount" :port "{port}" :history (constant 0.7N))
  {extra}
  (probe "pad_vel" :kind field-quantity :target "pad" :sample-rate 200Hz :quantity "velocity-x"
    :dims (1 0 -1 0 0 0)))"#
    )
}

#[test]
fn a_held_relaxation_bushing_with_declared_memory_runs_from_its_fsim() {
    let steel = material_pack("steel", 8000.0, 500.0, 193.0e9);
    let packs = CardPackSet::admit(vec![steel, relaxation_pack()]).unwrap();
    let card = packs.interfaces()[0].card().to_hex();
    let run = |port: &str, extra: &str| {
        let ledger = Ledger::open(":memory:").unwrap();
        let src = held_bushing(&card, port, extra);
        let out = run_assembly(
            &ledger,
            src.as_bytes(),
            &packs,
            AssemblyAnalysisMode::Transient,
            &mut step_clock(0.1),
            AssemblyCancelPlan::default(),
        );
        // Boxed: a refusal is large.
        out.map(|(out, _)| (ledger, out)).map_err(Box::new)
    };
    let (ledger, out) = run("drive_x", "").expect("the held bushing runs");
    assert_eq!(out.status, SolveRunStatus::Completed);
    let explained =
        explain_assembly_run(&ledger, SolveRunId::parse_hex(&out.run).unwrap()).unwrap();
    let receipt = explained.lowering_receipt_json.as_deref().unwrap();
    for shown in [
        "\"card_family\":\"relaxation-bushing\"",
        "\"parameter\":\"state.displacement_0_x\"",
        "\"parameter\":\"state.displacement_0_rz\"",
        "{\"source_name\":\"hold\",\"target\":\"",
        ".drive_x\",\"created_nodes\":[\"relaxation-network\"]",
    ] {
        assert!(receipt.contains(shown), "{shown}: {receipt}");
    }
    // A cell for the temperature sources below to read.
    let cell = r#"(part "cell" :kind solid :geometry (box 0.1m 0.1m 0.1m) :material "steel"
    :query-temp 300K :physics (thermal) :fidelity lumped)"#;
    let free = cell.replace("\"cell\"", "\"free\"");
    let warm =
        |a: &str| format!("(interface \"warm\" :kind temperature-source :a \"{a}\" :b \"cell\")");
    let unjoined = format!("{cell} {free} {}", warm("free"));
    let mount2 = bushing("mount2", ("rotor", "pad2"), &card);
    let twice = format!("{cell} {} {mount2} {}", steel_box("pad2"), warm("rotor"));
    let one = "interface:warm: feature_a must be on a part exactly one bushing joins";
    for (port, extra, named) in [
        ("drive_w", "", "has no exposed port 'drive_w'"),
        ("drive_x", unjoined.as_str(), one),
        ("drive_x", twice.as_str(), one),
    ] {
        let why = run(port, extra).err().unwrap();
        assert_eq!(why.code, "cli-assembly-lower", "{named}: {}", why.what);
        assert!(why.what.contains(named), "{named}: {}", why.what);
    }
}

fn gas_cubic_modes(alpha: f64, w_iso2: f64, w_ad2: f64) -> (f64, f64) {
    let f = |x: f64| ((x + alpha) * x + w_ad2) * x + alpha * w_iso2;
    let df = |x: f64| (3.0 * x + 2.0 * alpha) * x + w_ad2;
    let (mut lo, mut hi) = (-alpha, 0.0);
    for _ in 0..80 {
        let m = 0.5 * (lo + hi);
        if f(m) > 0.0 {
            hi = m;
        } else {
            lo = m;
        }
    }
    let mut r = 0.5 * (lo + hi);
    for _ in 0..8 {
        r -= f(r) / df(r);
    }
    let (sig, wn2) = (0.5 * (alpha + r), w_ad2 + r * (alpha + r));
    ((wn2 - sig * sig).sqrt(), sig / wn2.sqrt())
}

fn gas_fsim(extra_phys: &str, extra_iface: &str, force_n: f64) -> String {
    format!(
        r#"(assembly "gas_piston_test" (version 1) (units "SI") (seed 0x41)
  (versions (constellation "2026-07") (workspace "test-v1")) (budget (wall 60s) (mem 64MiB))
  (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*)) (reference-state assembled)
  (environment :ambient-temp 300K :ambient-pressure 90000Pa)
  (time-plan :horizon 0.01s :dt 0.0005s :sample-rate 2000Hz)
  (part "piston" :kind rigid-body :geometry (box 0.1m 0.1m 0.1m) :material "steel-music-wire-ligo-p1400229" :query-temp 293.15K :physics (mechanics) :fidelity lumped)
  (part "gas" :kind cavity :geometry (cylinder 0.15m 0.20m) :material "N2" :query-at (pressure 90000Pa) :physics (fluid{extra_phys}) :fidelity lumped)
  (interface "face" :kind piston-face :a "piston" :b "gas"){extra_iface}
  (source "ambient_pressure" :part "piston" :port "mechanics" :history (constant {force_n:.17}N))
  (probe "gas_temp" :kind field-quantity :target "gas" :quantity "temperature" :dims (0 0 0 1 0 0) :sample-rate 2000Hz)
  (analysis-options :modal-max 4))"#
    )
}

#[test]
fn gas_piston_fsim_adiabatic_finite_g_and_transient_oracles() {
    use fs_qty::semantic::{QuantityKind as Qk, SemanticType, ValueForm};
    let dir =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/assemblies/packs");
    let raw = |s: &str, p| RawCardPack {
        kind: CardPackKind::Material,
        source: s.into(),
        bytes: std::fs::read(dir.join(p)).unwrap(),
        expect: None,
    };
    let n2 = raw("n2.fsmcdpk", "nitrogen-n2-mcbride2002-lemmon2004.fsmcdpk");
    let st = raw(
        "steel.fsmcdpk",
        "steel-music-wire-heptonstall-ligo-p1400229.fsmcdpk",
    );
    let g_wall = 50.0_f64;
    let card = NormalizedMaterialCardPack::from_bytes(&n2.bytes).unwrap();
    let packs = CardPackSet::admit(vec![n2, st, thermal_interface_pack(g_wall)]).unwrap();
    let (r, h, t0, p0) = (0.15_f64, 0.20_f64, 300.0_f64, 90_000.0_f64);
    let area = std::f64::consts::PI * r * r;
    let (v0, m_p, f_amb) = (area * h, 7.8, -p0 * area);
    let qs = |k| fs_qty::QuantitySpec::semantic(SemanticType::new(k, ValueForm::Static));
    let qp = fs_matdb::QueryPoint::new()
        .with_quantity("temperature", qs(Qk::AbsoluteTemperature), t0)
        .unwrap()
        .with_quantity("pressure", qs(Qk::Pressure), p0)
        .unwrap();
    let parcel = fs_gas::component::parcel_from_material_card(card.card(), qp, v0).unwrap();
    let gs = parcel.state().state();
    let cv = gs.specific_heat_cp - gs.specific_gas_constant;
    let w_iso2 = p0 * area * area / (m_p * v0);
    let w_ad2 = (gs.specific_heat_cp / cv) * w_iso2;
    let ledger = Ledger::open(":memory:").unwrap();
    let try_run = |src: &str, m| -> Result<String, String> {
        let p = AssemblyCancelPlan::default();
        let (o, _) = run_assembly(&ledger, src.as_bytes(), &packs, m, &mut step_clock(0.1), p)
            .map_err(|e| e.what)?;
        assert_eq!(o.status, SolveRunStatus::Completed);
        Ok(o.run)
    };
    let run_one = |src: &str, m| {
        let id = SolveRunId::parse_hex(&try_run(src, m).unwrap()).unwrap();
        explain_assembly_run(&ledger, id).unwrap()
    };
    let pvals = |json: &Option<String>, key: &str| -> Vec<f64> {
        let raw = json.as_deref().unwrap();
        let after = raw.split(&format!("\"name\":\"{key}\"")).nth(1).unwrap();
        let arr = after.split("\"values\":[").nth(1).unwrap();
        let vals = arr.split(']').next().unwrap();
        vals.split(',').map(|s| s.trim().parse().unwrap()).collect()
    };
    let ad_src = gas_fsim("", "", f_amb);
    let ad_m = run_one(&ad_src, AssemblyAnalysisMode::Modal);
    let ad_f = pvals(&ad_m.probe_series_json, "modal.frequency_hz");
    let ad_z = pvals(&ad_m.probe_series_json, "modal.damping_ratio");
    let f_osc = *ad_f.iter().max_by(|a, b| a.total_cmp(b)).unwrap();
    assert!((f_osc - w_ad2.sqrt() / TAU).abs() / (w_ad2.sqrt() / TAU) <= 1e-6);
    assert!(ad_z[2].abs() <= 1e-6);

    let wall = "\n  (part \"frame\" :kind driven-body :geometry (box 0.1m 0.1m 0.1m) :material \"steel-music-wire-ligo-p1400229\" :query-temp 293.15K :physics (mechanics) :fidelity lumped)\n  (interface \"wall\" :kind thermal-contact :a \"gas\" :b \"frame\" :card \"AISI-304\" :family \"thermal-conductance\")";
    let th_src = gas_fsim(" thermal", wall, f_amb);
    let th_m = run_one(&th_src, AssemblyAnalysisMode::Modal);
    let th_f = pvals(&th_m.probe_series_json, "modal.frequency_hz");
    let th_z = pvals(&th_m.probe_series_json, "modal.damping_ratio");
    let (wd_ex, zeta_ex) = gas_cubic_modes(g_wall / (parcel.mass_kg() * cv), w_iso2, w_ad2);
    let idx = th_f.iter().position(|&v| v > 1.0).unwrap();
    assert!((th_f[idx] - wd_ex / TAU).abs() / (wd_ex / TAU) <= 1e-6);
    assert!((th_z[idx] - zeta_ex).abs() / zeta_ex <= 1e-5);

    let tr_src = gas_fsim("", "", f_amb - 500.0);
    let tr = run_one(&tr_src, AssemblyAnalysisMode::Transient);
    let temps = pvals(&tr.probe_series_json, "gas_temp");
    let t_end = *temps.last().unwrap();
    assert!(temps.len() == 21 && (temps[0] - t0).abs() < 1e-12 && t_end > t0 + 0.01);
    let res = tr.window_ledger_json.as_deref().unwrap();
    let rv = res.split("\"energy_residual_j\":").nth(1).unwrap();
    let rv = rv.split([',', '}']).next().unwrap();
    assert!(rv.parse::<f64>().unwrap().abs() <= 1e-9);

    let bad = ad_src.replace(
        ":b \"gas\")",
        ":b \"gas\" :card \"AISI-304\" :family \"thermal-conductance\")",
    );
    let err = try_run(&bad, AssemblyAnalysisMode::Modal).unwrap_err();
    assert!(err.contains("piston face takes no card"));
}

fn dry_friction_pack(tag: &str, chem: &str, mu_s: Option<f64>, mu_k: f64) -> RawCardPack {
    let mut claims = ClaimSet::new();
    let obs = reg_obs(&mut claims, tag);
    let val = ValidityDomain::unconstrained()
        .with("temperature", 250.0, 500.0)
        .with("normal_force", 1.0, 50.0)
        .with("sliding_speed", 0.0, 5.0);
    let add = |cs: &mut ClaimSet, name: &str, v: f64| {
        let mut c = claim(name, v, Dims::NONE, obs, name);
        c.validity = val.clone();
        cs.insert_claim(c).unwrap();
    };
    if let Some(s) = mu_s {
        add(&mut claims, "static-friction-coefficient", s);
    }
    add(&mut claims, "kinetic-friction-coefficient", mu_k);
    interface_pack(tag, [chem, chem], "dry-air", claims)
}

fn extract_probe_values(json: Option<&str>, key: &str) -> Vec<f64> {
    let after = json
        .unwrap()
        .split(&format!("\"name\":\"{key}\""))
        .nth(1)
        .unwrap();
    let arr = after.split("\"values\":[").nth(1).unwrap();
    arr.split(']')
        .next()
        .unwrap()
        .split(',')
        .map(|s| s.trim().parse().unwrap())
        .collect()
}

fn dry_friction_packs() -> CardPackSet {
    CardPackSet::admit(vec![
        material_pack("unit-mat", 1.0, 100.0, 200.0e9),
        dry_friction_pack("dry-good", "S355-ground", Some(0.5), 0.3),
        dry_friction_pack("dry-kin-only", "AISI-52100", None, 0.3),
    ])
    .unwrap()
}

fn run_dry_fsim(ledger: &Ledger, packs: &CardPackSet, src: &str) -> Result<String, String> {
    let p = AssemblyCancelPlan::default();
    let (o, _) = run_assembly(
        ledger,
        src.as_bytes(),
        packs,
        AssemblyAnalysisMode::Transient,
        &mut step_clock(0.1),
        p,
    )
    .map_err(|e| e.what)?;
    assert_eq!(o.status, SolveRunStatus::Completed);
    Ok(o.run)
}

#[test]
fn dry_friction_fsim_stick_slip_and_refusals() {
    let (packs, ledger) = (dry_friction_packs(), Ledger::open(":memory:").unwrap());
    let stick_slip = r#"(assembly "dry_stick_slip" (version 1) (units "SI") (seed 0x43)
  (versions (constellation "2026-07") (workspace "test-v1")) (budget (wall 60s) (mem 64MiB))
  (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*)) (reference-state assembled)
  (environment :ambient-temp 300K :ambient-pressure 101325Pa)
  (time-plan :horizon 3.5s :dt 0.002s :sample-rate 500Hz)
  (part "block" :kind rigid-body :geometry (box 1m 1m 1m) :material "unit-mat" :query-temp 300K :physics (mechanics thermal) :fidelity lumped :initial (:position (0m 0m 0m) :orientation (1 0 0 0) :linear-velocity (0.1m/s 0m/s 0m/s) :angular-velocity (0rad/s 0rad/s 0rad/s) :temperature 300K))
  (part "belt" :kind driven-body :geometry (box 1m 1m 1m) :material "unit-mat" :query-temp 300K :physics (mechanics) :fidelity lumped)
  (joint "spring" :kind prismatic :parent root :child "block" :axis (1 0 0) :stiffness 100N/m)
  (interface "pad" :kind mechanical-contact :a "block" :b "belt" :card "S355-ground" :family "dry-friction" :normal-force 10N :heat-partition-a 0.6 :initial-slip 0m/s :query-at (sliding_speed 0.1m/s))
  (source "belt_v" :part "belt" :port "mechanics" :history (constant 0.1m/s))
  (probe "block_pos" :kind field-quantity :target "block" :quantity "position" :dims (1 0 0 0 0 0) :sample-rate 500Hz)
  (probe "block_vel" :kind field-quantity :target "block" :quantity "velocity" :dims (1 0 -1 0 0 0) :sample-rate 500Hz)
  (probe "block_temp" :kind field-quantity :target "block" :quantity "temperature" :dims (0 0 0 1 0 0) :sample-rate 500Hz)
  (requirement "temp_cap" :probe "block_temp" :stat max :relation at-most :limit 305K))"#;
    let id = SolveRunId::parse_hex(&run_dry_fsim(&ledger, &packs, stick_slip).unwrap()).unwrap();
    let exp = explain_assembly_run(&ledger, id).unwrap();
    let pj = exp.probe_series_json.as_deref();
    let (pos, vel, tmp) = (
        extract_probe_values(pj, "block_pos"),
        extract_probe_values(pj, "block_vel"),
        extract_probe_values(pj, "block_temp"),
    );
    let (dt, x_k) = (2.0e-3_f64, 0.03_f64);
    let mut crosses = Vec::new();
    for i in 1..pos.len() {
        if (vel[i] - 0.1).abs() < 1e-6 && pos[i - 1] < x_k && pos[i] >= x_k {
            let f = (x_k - pos[i - 1]) / (pos[i] - pos[i - 1]);
            crosses.push((f64::from(u32::try_from(i - 1).unwrap()) + f) * dt);
        }
    }
    let delta: f64 = 0.02;
    let exact_t = 2.0 * delta / 0.1 + 2.0 * (std::f64::consts::PI - (delta * 100.0).atan()) / 10.0;
    let meas_t = (crosses[crosses.len() - 1] - crosses[1])
        / f64::from(u32::try_from(crosses.len() - 2).unwrap());
    assert!(crosses.len() >= 4 && (meas_t - exact_t).abs() < dt && *tmp.last().unwrap() > 300.005);
    let rcpt = exp.lowering_receipt_json.as_deref().unwrap();
    assert!(
        rcpt.contains("static-friction-coefficient.held")
            && rcpt.contains("kinetic-friction-coefficient.held")
            && rcpt.contains("heat_partition_a")
            && rcpt.contains("heat_partition_ambient")
            && exp.qoi_budget_json.as_deref().unwrap().contains("temp_cap")
    );
    let wj = exp.window_ledger_json.as_deref().unwrap();
    let res: f64 = wj
        .split("\"energy_residual_j\":")
        .nth(1)
        .unwrap()
        .split([',', '}'])
        .next()
        .unwrap()
        .parse()
        .unwrap();
    assert!(res.abs() <= 1e-12, "res={res}");
    for (bad, want) in [
        (
            stick_slip.replace(":kind mechanical-contact", ":kind thermal-contact"),
            "needs a mechanical contact",
        ),
        (
            stick_slip.replace(
                ":query-at (sliding_speed 0.1m/s)",
                ":query-at (normal_force 20N)",
            ),
            "repeats `normal_force`",
        ),
        (
            stick_slip.replace(":initial-slip 0m/s", ":initial-slip 1m/s"),
            "initial_slip",
        ),
        (
            stick_slip.replace("S355-ground", "AISI-52100"),
            "static-friction-coefficient",
        ),
    ] {
        let msg = run_dry_fsim(&ledger, &packs, &bad).unwrap_err();
        assert!(msg.contains(want), "expected `{want}` in: {msg}");
    }
}

/// A material pack `pack-<chem>` that also claims a thermal conductivity, as an effusivity
/// heat partition needs (`e = sqrt(k rho cp)`).
fn conductive_pack(chem: &str, rho: f64, cp: f64, k: f64) -> RawCardPack {
    let mut claims = ClaimSet::new();
    let obs = reg_obs(&mut claims, chem);
    for (name, value, dims, src) in [
        ("density", rho, DENSITY, "asm-rho"),
        ("specific_heat_capacity", cp, SPECIFIC_HEAT, "asm-cp"),
        ("thermal_conductivity", k, CONDUCTIVITY, "asm-k"),
        ("youngs_modulus", 200.0e9, MODULUS, "asm-e"),
    ] {
        claims
            .insert_claim(claim(name, value, dims, obs, src))
            .unwrap();
    }
    pack_from_claims(chem, claims)
}

#[test]
fn dry_friction_fsim_effusivity_heat_partition() {
    let ((ra, ca, ka), (rb, cb, kb)) = ((8000.0, 500.0, 16.0), (7850.0, 486.0, 50.0));
    let packs = CardPackSet::admit(vec![
        conductive_pack("eff-a", ra, ca, ka),
        conductive_pack("eff-b", rb, cb, kb),
        dry_friction_pack("dry-good", "S355-ground", Some(0.5), 0.3),
    ])
    .unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let src = r#"(assembly "dry_effusivity" (version 1) (units "SI") (seed 0x45)
  (versions (constellation "2026-07") (workspace "test-v1")) (budget (wall 60s) (mem 64MiB))
  (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*)) (reference-state assembled)
  (environment :ambient-temp 300K :ambient-pressure 101325Pa)
  (time-plan :horizon 0.6s :dt 0.002s :sample-rate 500Hz)
  (part "a" :kind rigid-body :geometry (box 0.1m 0.1m 0.1m) :material "eff-a" :query-temp 300K :physics (mechanics thermal) :fidelity lumped :initial (:position (0m 0m 0m) :orientation (1 0 0 0) :linear-velocity (1m/s 0m/s 0m/s) :angular-velocity (0rad/s 0rad/s 0rad/s) :temperature 300K))
  (part "b" :kind rigid-body :geometry (box 0.1m 0.1m 0.1m) :material "eff-b" :query-temp 300K :physics (mechanics thermal) :fidelity lumped :initial (:position (0m 0m 0m) :orientation (1 0 0 0) :linear-velocity (0m/s 0m/s 0m/s) :angular-velocity (0rad/s 0rad/s 0rad/s) :temperature 300K))
  (interface "pad" :kind mechanical-contact :a "a" :b "b" :card "S355-ground" :family "dry-friction" :normal-force 40N :heat-partition-effusivity 1)
  (probe "a_temp" :kind field-quantity :target "a" :quantity "temperature" :dims (0 0 0 1 0 0) :sample-rate 500Hz)
  (probe "b_temp" :kind field-quantity :target "b" :quantity "temperature" :dims (0 0 0 1 0 0) :sample-rate 500Hz))"#;
    let id = SolveRunId::parse_hex(&run_dry_fsim(&ledger, &packs, src).unwrap()).unwrap();
    let exp = explain_assembly_run(&ledger, id).unwrap();
    let pj = exp.probe_series_json.as_deref();
    let (ta, tb) = (
        *extract_probe_values(pj, "a_temp").last().unwrap(),
        *extract_probe_values(pj, "b_temp").last().unwrap(),
    );
    // Heat stored per body is rho cp V (T - T0); the common volume cancels in the share.
    let (dua, dub) = (ra * ca * (ta - 300.0), rb * cb * (tb - 300.0));
    let (ea, eb) = ((ka * ra * ca).sqrt(), (kb * rb * cb).sqrt());
    let share = dua / (dua + dub);
    assert!(dua > 0.0 && dub > 0.0, "a {ta} K, b {tb} K");
    assert!((share - ea / (ea + eb)).abs() <= 1e-8, "share {share}");
    let rcpt = exp.lowering_receipt_json.as_deref().unwrap();
    assert!(rcpt.contains("model: effusivity partition") && rcpt.contains("heat_partition_b"));
    let both = src.replace(
        ":heat-partition-effusivity 1",
        ":heat-partition-effusivity 1 :heat-partition-a 0.5",
    );
    let msg = run_dry_fsim(&ledger, &packs, &both).unwrap_err();
    assert!(msg.contains("cannot set both"), "got: {msg}");
}

/// A material pack `pack-<chem>` with a density and the elastic pair a Hertz contact reads.
fn elastic_pack(chem: &str, rho: f64, young: f64, poisson: f64) -> RawCardPack {
    let mut claims = ClaimSet::new();
    let obs = reg_obs(&mut claims, chem);
    for (name, value, dims, src) in [
        ("density", rho, DENSITY, "asm-rho"),
        ("young_modulus", young, MODULUS, "asm-e"),
        ("poisson_ratio", poisson, Dims::NONE, "asm-nu"),
    ] {
        claims
            .insert_claim(claim(name, value, dims, obs, src))
            .unwrap();
    }
    pack_from_claims(chem, claims)
}

#[test]
fn hertz_contact_fsim_elastic_impact_matches_closed_form() {
    let (young, nu, rho, r, v0, dt) = (200.0e9_f64, 0.3_f64, 7850.0_f64, 0.01_f64, 0.5_f64, 5.0e-7);
    let packs = CardPackSet::admit(vec![elastic_pack("hz-steel", rho, young, nu)]).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let src = r#"(assembly "hertz_impact" (version 1) (units "SI") (seed 0x46)
  (versions (constellation "2026-07") (workspace "test-v1")) (budget (wall 60s) (mem 64MiB))
  (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*)) (reference-state assembled)
  (environment :ambient-temp 293.15K :ambient-pressure 101325Pa)
  (time-plan :horizon 0.00015s :dt 0.0000005s :sample-rate 2000000Hz)
  (part "ball" :kind rigid-body :geometry (sphere 0.01m) :material "hz-steel" :query-temp 293.15K :physics (mechanics) :fidelity lumped :initial (:position (0m 0m 0m) :orientation (1 0 0 0) :linear-velocity (-0.5m/s 0m/s 0m/s) :angular-velocity (0rad/s 0rad/s 0rad/s)))
  (part "anvil" :kind driven-body :geometry (box 0.1m 0.1m 0.02m) :material "hz-steel" :query-temp 293.15K :physics (mechanics) :fidelity lumped)
  (interface "hit" :kind hertz-contact :a "ball" :b "anvil" :maximum-force 10000N :maximum-penetration 0.001m)
  (source "hold" :part "anvil" :port "mechanics" :history (constant 0m/s))
  (probe "x" :kind field-quantity :target "ball" :quantity "position" :dims (1 0 0 0 0 0) :sample-rate 2000000Hz)
  (probe "v" :kind field-quantity :target "ball" :quantity "velocity" :dims (1 0 -1 0 0 0) :sample-rate 2000000Hz))"#;
    let id = SolveRunId::parse_hex(&run_dry_fsim(&ledger, &packs, src).unwrap()).unwrap();
    let exp = explain_assembly_run(&ledger, id).unwrap();
    let pj = exp.probe_series_json.as_deref();
    let (x, v) = (extract_probe_values(pj, "x"), extract_probe_values(pj, "v"));
    // The junction hands each member port minus its body's velocity, so the closure grows at
    // v_b - v_a: the ball (end a) above the anvil (end b) strikes it at -v0, and the closure is
    // minus the ball's position. Hertz impact (Johnson 1985, sec. 11.4): K = 4/3 E* sqrt(R), the
    // peak closure from m v0^2 / 2 = 2/5 K d^(5/2), the contact time 2 d / v0 times
    // int_0^1 (1 - s^(5/2))^(-1/2) ds.
    let m = rho * 4.0 / 3.0 * std::f64::consts::PI * r.powi(3);
    let k = 4.0 / 3.0 * young / (2.0 * (1.0 - nu * nu)) * r.sqrt();
    let d_max = (5.0 * m * v0 * v0 / (4.0 * k)).powf(0.4);
    let t_c = 2.943_275_184_324_705 * d_max / v0;
    let peak = -x.iter().copied().fold(f64::MAX, f64::min);
    let v_end = *v.last().unwrap();
    let last = v
        .iter()
        .rposition(|w| (w - v_end).abs() > 1e-9 * v0)
        .unwrap();
    let t_end = f64::from(u32::try_from(last + 1).unwrap()) * dt;
    assert!(
        ((peak - d_max) / d_max).abs() <= 1e-3,
        "peak {peak} m, Hertz {d_max} m"
    );
    assert!(((v_end - v0) / v0).abs() <= 1e-6, "rebound {v_end} m/s");
    assert!(
        (t_end - t_c).abs() <= 3.0 * dt,
        "contact {t_end} s, Hertz {t_c} s"
    );
    let rcpt = exp.lowering_receipt_json.as_deref().unwrap();
    assert!(rcpt.contains("effective_modulus") && rcpt.contains("effective_radius"));
    for (bad, want) in [
        (
            src.replace(":kind hertz-contact", ":kind hertz-contact :card \"x\""),
            "takes no card",
        ),
        (
            src.replace(":kind hertz-contact", ":kind heat-route"),
            ":maximum-force needs a hertz contact",
        ),
    ] {
        let msg = run_dry_fsim(&ledger, &packs, &bad).unwrap_err();
        assert!(msg.contains(want), "expected `{want}` in: {msg}");
    }
}

#[test]
fn dry_friction_fsim_reversing_table_belt() {
    let (packs, ledger) = (dry_friction_packs(), Ledger::open(":memory:").unwrap());
    let rev = r#"(assembly "dry_reversing" (version 1) (units "SI") (seed 0x44)
  (versions (constellation "2026-07") (workspace "test-v1")) (budget (wall 60s) (mem 64MiB))
  (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*)) (reference-state assembled)
  (environment :ambient-temp 300K :ambient-pressure 101325Pa)
  (time-plan :horizon 0.8s :dt 0.002s :sample-rate 500Hz)
  (part "block" :kind rigid-body :geometry (box 1m 1m 1m) :material "unit-mat" :query-temp 300K :physics (mechanics thermal) :fidelity lumped :initial (:position (0m 0m 0m) :orientation (1 0 0 0) :linear-velocity (0m/s 0m/s 0m/s) :angular-velocity (0rad/s 0rad/s 0rad/s) :temperature 300K))
  (part "belt" :kind driven-body :geometry (box 1m 1m 1m) :material "unit-mat" :query-temp 300K :physics (mechanics) :fidelity lumped)
  (interface "pad" :kind mechanical-contact :a "block" :b "belt" :card "S355-ground" :family "dry-friction" :normal-force 10N)
  (source "belt_v" :part "belt" :port "mechanics" :history (table (0s 0m/s) (0.2s 0.6m/s) (0.4s -0.6m/s) (0.8s -0.6m/s)))
  (probe "block_vel" :kind field-quantity :target "block" :quantity "velocity" :dims (1 0 -1 0 0 0) :sample-rate 500Hz)
  (probe "block_temp" :kind field-quantity :target "block" :quantity "temperature" :dims (0 0 0 1 0 0) :sample-rate 500Hz))"#;
    let id = SolveRunId::parse_hex(&run_dry_fsim(&ledger, &packs, rev).unwrap()).unwrap();
    let rexp = explain_assembly_run(&ledger, id).unwrap();
    let pj = rexp.probe_series_json.as_deref();
    let (rvel, rtmp) = (
        extract_probe_values(pj, "block_vel"),
        extract_probe_values(pj, "block_temp"),
    );
    assert!((rvel[100] - 0.6).abs() <= 1e-12 && (rtmp[100] - 300.0).abs() <= 1e-12);
    assert!((rvel.last().unwrap() + 0.6).abs() <= 1e-12 && *rtmp.last().unwrap() > 300.003);
    let bad = rev.replace(
        "(table (0s 0m/s) (0.2s 0.6m/s) (0.4s -0.6m/s) (0.8s -0.6m/s))",
        "(table (0s 0m/s))",
    );
    let msg = run_dry_fsim(&ledger, &packs, &bad).unwrap_err();
    assert!(msg.contains("at least two"), "got: {msg}");
}

/// A revolute steel box turned 0.6 rad about z and spinning at 10 rad/s, on a
/// spring of 2000 N m/rad, with a probe of `quantity` and `dims`.
fn rotor(quantity: &str, dims: &str) -> String {
    let (c, s) = (0.3_f64.cos(), 0.3_f64.sin());
    format!(
        r#"(assembly "rotor"
  (version 1) (units "SI") (seed 0x23)
  (versions (constellation "2026-07") (workspace "stage2-e2e"))
  (budget (wall 60s) (mem 64MiB))
  (capability :cores 1cores :mem 64MiB :wall 60s :ops (assembly.*))
  (reference-state unloaded) (environment earth-lab-air)
  (time-plan :horizon 0.05s :dt 0.001s :sample-rate 1000Hz)
  (part "rotor" :kind rigid-body :geometry (box 0.5m 0.4m 0.1m)
    :material "AISI-304" :query-temp 293.15K :physics (mechanics) :fidelity lumped
    :initial (:position (0m 0m 0m) :orientation ({c:?} 0 0 {s:?})
      :linear-velocity (0m/s 0m/s 0m/s) :angular-velocity (0rad/s 0rad/s 10rad/s)))
  (joint "pin" :kind revolute :parent root :child "rotor" :axis (0 0 1)
    :stiffness 2000N*m/rad)
  (probe "turn" :kind field-quantity :target "rotor" :sample-rate 1000Hz
    :quantity "{quantity}" :dims ({dims}))
)"#
    )
}

/// A turning body reports its turn (bead wkb9e): `angle` reads the swing
/// `theta0 + (w0 / omega) sin(n phi)`, `tan(phi / 2) = omega dt / 2`, with `I = rho V
/// (lx^2 + ly^2) / 12` from the card's density alone, within O4c's per-step bound;
/// `position` and `velocity` refuse on it, naming `angle`.
#[test]
fn a_turning_body_reports_its_angle_and_refuses_a_position() {
    let packs =
        CardPackSet::admit(vec![material_pack("AISI-304", 8000.0, 500.0, 193.0e9)]).unwrap();
    let ledger = Ledger::open(":memory:").unwrap();
    let (out, _) = run_assembly(
        &ledger,
        rotor("angle", "0 0 0 0 0 0").as_bytes(),
        &packs,
        AssemblyAnalysisMode::Transient,
        &mut step_clock(0.1),
        AssemblyCancelPlan::default(),
    )
    .expect("an angle probe on a revolute box runs");
    assert_eq!(out.status, SolveRunStatus::Completed);
    let exp = explain_assembly_run(&ledger, SolveRunId::parse_hex(&out.run).unwrap()).unwrap();
    let turn = extract_probe_values(exp.probe_series_json.as_deref(), "turn");
    let inertia = 8000.0_f64 * (0.5 * 0.4 * 0.1) * (0.5 * 0.5 + 0.4 * 0.4) / 12.0;
    let (theta0, w0, omega) = (0.6, 10.0, (2000.0 / inertia).sqrt());
    let phi = 2.0 * (omega * 0.001 / 2.0).atan();
    assert_eq!(turn.len(), 51);
    for (n, got) in (0_u32..).zip(&turn) {
        let k = f64::from(n);
        let want = theta0 + w0 / omega * (k * phi).sin();
        let tol = 16.0 * f64::EPSILON * (theta0 + w0 / omega) * k.max(1.0);
        assert!((got - want).abs() <= tol, "{n}: {got} {want}");
    }
    // Killed after the run, the resumed probe stage reads the turn from the retained trace.
    let killed = Ledger::open(":memory:").unwrap();
    let plan = AssemblyCancelPlan {
        cancel_after: Some(AssemblyStage::Run),
        cancel_pre_publish: None,
    };
    let src = rotor("angle", "0 0 0 0 0 0");
    let mut clk = step_clock(0.1);
    let err = run_assembly(
        &killed,
        src.as_bytes(),
        &packs,
        AssemblyAnalysisMode::Transient,
        &mut clk,
        plan,
    )
    .expect_err("cancelled after the run");
    assert_eq!(err.code, "cli-solve-cancelled");
    let run_id = SolveRunId::parse_hex(&out.run).unwrap();
    let (resumed, _) = resume_assembly(&killed, run_id, &mut clk, AssemblyCancelPlan::default())
        .expect("the probe stage resumes from the retained trace");
    assert_eq!(
        resumed.prior_stages as usize,
        AssemblyStage::Run as usize + 1
    );
    assert_eq!(
        resumed.stages.last().unwrap().receipt,
        out.stages.last().unwrap().receipt
    );
    for (quantity, dims) in [("position", "1 0 0 0 0 0"), ("velocity", "1 0 -1 0 0 0")] {
        let e = run_assembly(
            &Ledger::open(":memory:").unwrap(),
            rotor(quantity, dims).as_bytes(),
            &packs,
            AssemblyAnalysisMode::Transient,
            &mut step_clock(0.1),
            AssemblyCancelPlan::default(),
        )
        .expect_err("a position or velocity probe on a turning body refuses");
        assert_eq!(
            (e.code, e.stage),
            ("cli-assembly-probe-unsupported", Some("probe"))
        );
        assert!(e.what.contains("probe `angle`"), "{quantity}: {}", e.what);
    }
}
