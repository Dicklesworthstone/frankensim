//! IR-speaking conformance conformance (bead 6nb.8 slice 1): real
//! programs route through the REAL admission path, refusals become
//! structured failing records with the deterministic diagnosis, admitted
//! kernels compare content-addressed artifacts, and the canonical
//! program identities ride every record as the cross-agent negotiation
//! anchor.

use fs_casebook::ToleranceSpec;
use fs_ir::admission::{
    AdmissionContext, RegimePolicy, SealedSessionCapability, SessionCapability,
};
use fs_ir::conformance::{IrCase, artifact_hash, run_ir_suite};
use std::collections::BTreeMap;

const SPOUT: &str = r#"(study "spout-laminar-v3"
  (seed 0x5EED0001) (versions (constellation :lock "2026-07"))
  (budget (wall 2h) (mem 96GiB) (qoi-rel-error 2e-2))
  (let vessel (frep (revolve (cheb-profile "body.chb")) (fillet :edge lip :r 3mm)))
  (let lever  (xform.level-set-velocity vessel :band 12mm :dof 4096))
  (let pour   (flux.free-surface-lbm vessel
                (fluid :model (carreau :mu0 0.12Pa*s :n 0.8) :sigma 0.061N/m)
                (schedule :rate 0.5L/s :tilt (ramp 0deg 65deg 3s))))
  (let J (min (perturbation-growth pour :at lip :modes (1 .. 8))))
  (ascent.optimize J :over lever :method (lbfgs :m 17)
    :until (any (grad-norm 1e-5) (e-value 20) (budget-exhausted))
    :emit (pareto ledger report)))"#;

fn permissive_context() -> AdmissionContext<'static> {
    AdmissionContext {
        router: None,
        cost_freshness: None,
        chart_requirements: Vec::new(),
        cost_models: BTreeMap::new(),
        capability: Some(SealedSessionCapability::caller_declared(
            SessionCapability {
                ops: vec![
                    "flux.*".to_owned(),
                    "ascent.*".to_owned(),
                    "xform.*".to_owned(),
                ],
                cores: 1,
                mem_bytes: 1 << 30,
                wall_s: 3_600.0,
            },
        )),
        regime: None,
        regime_policy: RegimePolicy::Warn,
    }
}

fn denying_context() -> AdmissionContext<'static> {
    AdmissionContext {
        router: None,
        cost_freshness: None,
        chart_requirements: Vec::new(),
        cost_models: BTreeMap::new(),
        capability: Some(SealedSessionCapability::caller_declared(
            SessionCapability {
                ops: vec!["topo.*".to_owned()],
                cores: 1,
                mem_bytes: 1 << 30,
                wall_s: 3_600.0,
            },
        )),
        regime: None,
        regime_policy: RegimePolicy::Warn,
    }
}

fn spout_case(id: &'static str, artifact: &'static [u8]) -> IrCase {
    IrCase {
        id,
        program_sexpr: SPOUT.to_owned(),
        tolerance: ToleranceSpec::Structural,
        expected_artifact: artifact_hash(artifact),
        kernel: Box::new(move |report| {
            assert!(report.admitted, "kernel runs only after real admission");
            artifact.to_vec()
        }),
    }
}

#[test]
fn admitted_programs_run_kernels_and_match_artifacts() {
    let cx = permissive_context();
    let report = run_ir_suite(
        "ir-conformance-demo",
        vec![spout_case("ir-001-spout-artifact", b"pour-artifact-v1")],
        &cx,
    );
    report.assert_green();
    let record = &report.records[0];
    assert!(record.details.contains("admitted through the real path"));
    assert!(
        record.evidence.iter().any(|e| e.starts_with("ir-raw:")),
        "raw canonical identity rides the record"
    );
    assert!(
        record.evidence.iter().any(|e| e.starts_with("ir-lowered:")),
        "lowered canonical identity rides the record"
    );
}

#[test]
fn capability_refusal_is_a_structured_failing_record_not_a_side_door() {
    let cx = denying_context();
    let report = run_ir_suite(
        "ir-conformance-denied",
        vec![spout_case("ir-002-denied", b"never-produced")],
        &cx,
    );
    assert!(!report.all_passed(), "a refused program cannot read green");
    let record = &report.records[0];
    assert!(record.details.contains("REAL admission refused"));
    assert!(
        record.details.contains("REJECT"),
        "the deterministic diagnosis rides the record: {}",
        record.details
    );
    assert!(record.evidence.iter().any(|e| e.starts_with("ir-raw:")));
    assert!(
        !record.evidence.iter().any(|e| e.starts_with("ir-lowered:")),
        "no lowered authority identity exists for a refused program"
    );
}

#[test]
fn artifact_drift_fails_with_both_content_addresses() {
    let cx = permissive_context();
    let mut case = spout_case("ir-003-drift", b"expected-bytes");
    case.kernel = Box::new(|_| b"drifted-bytes".to_vec());
    let report = run_ir_suite("ir-conformance-drift", vec![case], &cx);
    assert!(!report.all_passed());
    let record = &report.records[0];
    assert!(record.details.contains("artifact drifted"));
    assert!(record.details.contains("produced"));
    assert!(record.details.contains("expected"));
}

#[test]
fn malformed_programs_refuse_before_admission() {
    let cx = permissive_context();
    let case = IrCase {
        id: "ir-004-malformed",
        program_sexpr: "(study \"broken".to_owned(),
        tolerance: ToleranceSpec::Structural,
        expected_artifact: artifact_hash(b"unused"),
        kernel: Box::new(|_| unreachable!("kernel must not run for a refused parse")),
    };
    let report = run_ir_suite("ir-conformance-malformed", vec![case], &cx);
    assert!(!report.all_passed());
    assert!(
        report.records[0]
            .details
            .contains("refused before admission")
    );
}

#[test]
fn identical_programs_negotiate_on_identical_canonical_identities() {
    let cx = permissive_context();
    let a = run_ir_suite(
        "ir-agent-a",
        vec![spout_case("ir-005-negotiate", b"artifact")],
        &cx,
    );
    let b = run_ir_suite(
        "ir-agent-b",
        vec![spout_case("ir-005-negotiate", b"artifact")],
        &cx,
    );
    let identities = |r: &fs_casebook::SuiteReport| -> Vec<String> {
        r.records[0]
            .evidence
            .iter()
            .filter(|e| e.starts_with("ir-"))
            .cloned()
            .collect()
    };
    assert_eq!(
        identities(&a),
        identities(&b),
        "two agents agree on the exact canonical program, not on prose"
    );
    assert_eq!(a.records[0].inputs_digest, b.records[0].inputs_digest);
}
