//! Machine-IR E0 PR-5 one-way domain-lowering conformance (G0/G3/G5).

#[path = "support/machine_stack.rs"]
mod machine_stack;

use core::num::NonZeroU64;

use fs_blake3::{ContentHash, hash_domain};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_ir::machine::lowering::*;
use fs_ir::machine::{ClockId, MachineElementKind};
use fs_ledger::{EdgeRole, FiveExplicits, Ledger, MAIN_BRANCH, OpOutcome};
use fs_package::{
    Claim, EvidencePackage, MAX_SEMANTIC_WITNESS_PAYLOAD_BYTES, Provenance, SemanticWitness,
};
use fs_qty::{Dims, QtyAny};
use fs_scenario::payload::{
    OrientationParity, Payload, PayloadId, PayloadMeta, QuantityContract, ReferenceSemantics,
    SampleSource, VectorPayload,
};
use fs_scenario::{
    BcKind, BcValue, BoundaryCondition, Environment, FrameId, LoadCase, Physics, Scenario,
    ValidationBudget,
};

const EXPLICITS: FiveExplicits<'static> = FiveExplicits {
    seed: b"machine-domain-lowering-test-v1",
    versions: r#"{"fs-ir":"machine-lowering-v1"}"#,
    budget: r#"{"artifacts":4,"bytes":262144}"#,
    capability: r#"{"ops":["machine-domain-lowering"]}"#,
};

fn nz(value: u64) -> NonZeroU64 {
    NonZeroU64::new(value).expect("test value is nonzero")
}

fn digest(label: &str) -> ContentHash {
    hash_domain(
        "org.frankensim.fs-ir.machine-lowering-test.v1",
        label.as_bytes(),
    )
}

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x5eed,
                kernel_id: 5,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn with_cancelled_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    gate.request();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 0x5eed,
                kernel_id: 5,
                tile: 0,
                iteration: 1,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn scenario(name: &str) -> Scenario {
    Scenario::new(name, 0x5eed, Environment::earth_lab())
}

fn typed_magnetic_scenario(name: &str, region: &str) -> Scenario {
    const MAGNETIC_VECTOR_POTENTIAL: Dims = Dims([1, 1, -2, 0, -1, 0]);

    let meta = PayloadMeta::new(
        QuantityContract::Dimensions(MAGNETIC_VECTOR_POTENTIAL),
        PayloadId::new("basis/world-cartesian").expect("canonical basis ID"),
        FrameId(0),
        OrientationParity::Even,
        ReferenceSemantics::Continuous,
    )
    .expect("valid typed-payload metadata");
    let payload = Payload::Vector(
        VectorPayload::new(
            meta,
            SampleSource::fixed(vec![
                QtyAny::new(1.0, MAGNETIC_VECTOR_POTENTIAL),
                QtyAny::new(2.0, MAGNETIC_VECTOR_POTENTIAL),
                QtyAny::new(3.0, MAGNETIC_VECTOR_POTENTIAL),
            ]),
        )
        .expect("valid magnetic vector payload"),
    );
    let mut scenario = scenario(name);
    scenario.base_bcs.push(BoundaryCondition {
        region: region.to_string(),
        physics: Physics::Magnetics,
        kind: BcKind::MagneticVectorPotential,
        value: Some(BcValue::Typed(payload)),
        compatibility: None,
        frame: 0,
    });
    scenario
}

fn artifact(hash_label: &str) -> MachineDomainArtifactRefV1 {
    MachineDomainArtifactRefV1::new(
        MachineDomainArtifactId::new("domain/behavior-semantics").expect("artifact id"),
        MachineDomainArtifactKindV1::External,
        "fixture/behavior-semantics",
        nz(1),
        digest(hash_label),
    )
    .expect("artifact ref")
}

fn law(hash_label: &str) -> MachineCrosswalkLawRefV1 {
    MachineCrosswalkLawRefV1::new("fixture/machine-scenario-law", nz(1), digest(hash_label))
        .expect("crosswalk law")
}

fn policy(hash_label: &str) -> MachineLoweringPolicyRefV1 {
    MachineLoweringPolicyRefV1::new("fixture/machine-lowering-policy", nz(1), digest(hash_label))
        .expect("lowering policy")
}

fn complete_draft(
    sources: &[MachineLoweringSourceV1],
    scenario: Scenario,
    artifact_hash: &str,
    law_hash: &str,
    policy_hash: &str,
) -> MachineDomainLoweringDraftV1 {
    let artifact = artifact(artifact_hash);
    let artifact_id = artifact.id().clone();
    let law = law(law_hash);
    let crosswalks = sources
        .iter()
        .cloned()
        .map(|source| {
            let mut targets = vec![MachineDomainTargetV1::Scenario(
                MachineScenarioLocatorV1::Root,
            )];
            if matches!(source, MachineLoweringSourceV1::Behavior(_)) {
                targets.push(
                    MachineDomainTargetV1::external(artifact_id.clone(), "behavior/root")
                        .expect("external target"),
                );
            }
            MachineDomainCrosswalkEntryV1::new(source, targets, law.clone())
                .expect("complete crosswalk row")
        })
        .collect();
    MachineDomainLoweringDraftV1 {
        scenario,
        external_artifacts: vec![artifact],
        crosswalks,
        policy: policy(policy_hash),
    }
}

fn record_payload(ledger: &Ledger, payload: &[u8]) {
    let operation = ledger
        .begin_op(
            None,
            r#"{"op":"machine-domain-lowering-v1"}"#,
            &EXPLICITS,
            10,
        )
        .expect("begin deterministic lowering op");
    let artifact = ledger
        .put_artifact(
            "fs-ir/machine-domain-lowering-v1",
            payload,
            Some(r#"{"schema_version":1}"#),
        )
        .expect("store lowering payload");
    ledger
        .link(operation, &artifact.hash, EdgeRole::Out)
        .expect("link lowering output");
    ledger
        .finish_op(operation, OpOutcome::Ok, None, 11)
        .expect("finish lowering op");
}

#[test]
#[allow(clippy::too_many_lines)]
fn g0_complete_one_way_projection_replays_and_preserves_all_identity_domains() {
    let (graph, behavior, assurance) = machine_stack::admitted_machine_stack();
    let sources = required_machine_lowering_sources_v1(&graph, &behavior, &assurance)
        .expect("exact stack source inventory");
    assert!(sources.len() > 12, "fixture should exercise durable IDs");
    assert!(sources.iter().any(|source| matches!(
        source,
        MachineLoweringSourceV1::Element(element)
            if element.kind() == MachineElementKind::StateSlot
    )));
    assert!(
        sources
            .iter()
            .any(|source| matches!(source, MachineLoweringSourceV1::Sensor(_)))
    );
    assert!(
        sources
            .iter()
            .any(|source| matches!(source, MachineLoweringSourceV1::AccountingWindow(_)))
    );

    let decision = with_cx(|cx| {
        admit_machine_domain_lowering_with_decision_v1(
            &graph,
            &behavior,
            &assurance,
            complete_draft(
                &sources,
                scenario("fixture-machine"),
                "artifact-a",
                "law-a",
                "policy-a",
            ),
            ValidationBudget::default(),
            cx,
        )
    });
    assert_eq!(decision.code(), "MachineDomainLoweringAdmitted");
    assert_eq!(decision.refusal_code(), None);
    let submitted = decision.submitted_counts();
    assert_eq!(submitted.external_artifacts(), 1);
    assert_eq!(submitted.crosswalks(), sources.len());
    assert_eq!(submitted.crosswalk_targets(), Some(sources.len() + 1));
    let admitted = decision.into_result().expect("complete projection admits");
    assert_eq!(admitted.graph(), graph.identity());
    assert_eq!(admitted.behavior(), behavior.identity());
    assert_eq!(admitted.assurance(), assurance.identity());
    assert_eq!(admitted.crosswalks().len(), sources.len());
    assert_eq!(admitted.external_artifacts().len(), 1);
    assert!(admitted.scenario().validate().is_empty());
    assert_eq!(
        fs_blake3::hash_bytes(admitted.canonical_scenario_ir().as_bytes()),
        admitted.scenario_hash()
    );
    let replay = admitted.verify_replay().expect("exact replay is stable");
    assert_eq!(replay.lowering(), admitted.identity());
    assert_eq!(replay.scenario_hash(), admitted.scenario_hash());
    assert_eq!(
        replay.portable_payload_hash(),
        admitted.portable_payload_hash()
    );

    assert!(
        admitted.portable_payload().len() <= MAX_SEMANTIC_WITNESS_PAYLOAD_BYTES,
        "this small fixture must fit the package witness's independent cap"
    );
    let witness = SemanticWitness::new(
        "fs-ir/machine-domain-lowering",
        MACHINE_DOMAIN_LOWERING_SCHEMA_VERSION_V1,
        admitted.portable_payload().to_vec(),
    );
    let witness_hash = witness.content_hash();
    let package = EvidencePackage::new(Provenance::new("fixture-commit", "fixture-lock"))
        .with_claim(Claim::from_portable_certificate(
            "machine-domain-lowering",
            "exact transport of an admitted Machine-domain projection",
            0.0,
            0.0,
            "fs-ir-test",
            witness,
        ));
    let package_root = package
        .verify_structural_integrity()
        .expect("package structural integrity");
    assert_eq!(
        package_root,
        package.try_merkle_root().expect("package Merkle root")
    );
    let json = package.to_json().expect("canonical package JSON");
    let decoded = EvidencePackage::from_json(&json).expect("package round trip");
    assert_eq!(decoded, package);
    assert_eq!(decoded.to_json().expect("reprint package"), json);
    let decoded_witness = decoded.declared_claims_unverified()[0]
        .declared_semantic_witness_unverified()
        .expect("portable witness retained");
    assert_eq!(decoded_witness.family(), "fs-ir/machine-domain-lowering");
    assert_eq!(
        decoded_witness.schema_version(),
        MACHINE_DOMAIN_LOWERING_SCHEMA_VERSION_V1
    );
    assert_eq!(
        decoded_witness.canonical_payload(),
        admitted.portable_payload()
    );
    assert_eq!(decoded_witness.content_hash(), witness_hash);
    let portable = parse_machine_domain_portable_payload_v1(decoded_witness.canonical_payload())
        .expect("independent portable payload decode");
    assert_eq!(portable.graph(), graph.identity());
    assert_eq!(portable.behavior(), behavior.identity());
    assert_eq!(portable.assurance(), assurance.identity());
    assert_eq!(portable.scenario_hash(), admitted.scenario_hash());
    assert_eq!(portable.manifest_bytes(), admitted.canonical_manifest());
    assert_eq!(
        portable.canonical_scenario_ir(),
        admitted.canonical_scenario_ir()
    );
    assert_eq!(portable.payload_hash(), admitted.portable_payload_hash());
    assert_eq!(portable.external_artifact_rows(), 1);
    assert_eq!(portable.crosswalk_rows(), sources.len());
    assert_ne!(witness_hash, package_root, "identity domains are distinct");
    assert_ne!(
        admitted.portable_payload_hash(),
        witness_hash,
        "plain artifact and semantic-witness hashes are distinct domains"
    );

    let mut wrong_outer_magic = admitted.portable_payload().to_vec();
    wrong_outer_magic[0] ^= 1;
    assert!(matches!(
        parse_machine_domain_portable_payload_v1(&wrong_outer_magic),
        Err(MachineDomainPortablePayloadError::Magic { field: "payload" })
    ));
    let mut trailing_payload = admitted.portable_payload().to_vec();
    trailing_payload.push(0);
    assert!(matches!(
        parse_machine_domain_portable_payload_v1(&trailing_payload),
        Err(MachineDomainPortablePayloadError::TrailingBytes { field: "payload" })
    ));
    let mut wrong_manifest_magic = admitted.portable_payload().to_vec();
    let manifest_offset = wrong_manifest_magic
        .windows(admitted.canonical_manifest().len())
        .position(|window| window == admitted.canonical_manifest())
        .expect("nested manifest bytes");
    wrong_manifest_magic[manifest_offset] ^= 1;
    assert!(matches!(
        parse_machine_domain_portable_payload_v1(&wrong_manifest_magic),
        Err(MachineDomainPortablePayloadError::Magic { field: "manifest" })
    ));

    let mut tampered_json = json.clone().into_bytes();
    let marker = b"\"payload_hex\":\"";
    let start = tampered_json
        .windows(marker.len())
        .position(|window| window == marker)
        .expect("payload marker")
        + marker.len();
    tampered_json[start] = if tampered_json[start] == b'0' {
        b'1'
    } else {
        b'0'
    };
    assert!(
        EvidencePackage::from_json(core::str::from_utf8(&tampered_json).expect("still UTF-8 JSON"))
            .is_err(),
        "payload tamper must fail declared-root or certificate binding"
    );

    let original_ledger = Ledger::open(":memory:").expect("open original ledger");
    let replay_ledger = Ledger::open(":memory:").expect("open replay ledger");
    record_payload(&original_ledger, admitted.portable_payload());
    record_payload(&replay_ledger, decoded_witness.canonical_payload());
    let verdict = original_ledger
        .replay_verdict(MAIN_BRANCH, &replay_ledger, MAIN_BRANCH)
        .expect("compare replay ledgers");
    assert!(
        verdict.is_replay_clean(),
        "exact payload replay: {verdict:?}"
    );

    let mut mutated_payload = admitted.portable_payload().to_vec();
    let last = mutated_payload.last_mut().expect("nonempty payload");
    *last ^= 1;
    assert!(
        parse_machine_domain_portable_payload_v1(&mutated_payload).is_err(),
        "standalone portable decoder must reject scenario/hash drift"
    );
    let divergent_ledger = Ledger::open(":memory:").expect("open divergent ledger");
    record_payload(&divergent_ledger, &mutated_payload);
    let divergent = original_ledger
        .replay_verdict(MAIN_BRANCH, &divergent_ledger, MAIN_BRANCH)
        .expect("compare divergent replay");
    assert!(!divergent.is_replay_clean());
    assert_eq!(divergent.deterministic_mismatches.len(), 1);
}

#[test]
fn g0_typed_scenario_payload_is_bound_into_l6_artifact_and_crosswalk() {
    const REGION: &str = "magnetic-interface";

    let (graph, behavior, assurance) = machine_stack::admitted_machine_stack();
    let sources = required_machine_lowering_sources_v1(&graph, &behavior, &assurance)
        .expect("exact stack source inventory");
    let expected_scenario = typed_magnetic_scenario("typed-machine", REGION);
    let mut draft = complete_draft(
        &sources,
        expected_scenario.clone(),
        "typed-artifact",
        "typed-law",
        "typed-policy",
    );

    let mapped_source = draft.crosswalks[0].source().clone();
    let mapped_law = draft.crosswalks[0].law().clone();
    let region_target =
        MachineDomainTargetV1::Scenario(MachineScenarioLocatorV1::Region(REGION.into()));
    let mut mapped_targets = draft.crosswalks[0].targets().to_vec();
    mapped_targets.push(region_target.clone());
    draft.crosswalks[0] =
        MachineDomainCrosswalkEntryV1::new(mapped_source.clone(), mapped_targets, mapped_law)
            .expect("typed scenario region is an explicit crosswalk target");

    let admitted = with_cx(|cx| {
        admit_machine_domain_lowering_v1(
            &graph,
            &behavior,
            &assurance,
            draft,
            ValidationBudget::default(),
            cx,
        )
    })
    .expect("typed scenario projection admits");

    assert_eq!(admitted.scenario(), &expected_scenario);
    assert!(
        admitted
            .canonical_scenario_ir()
            .contains("(typed :version 1 \"")
    );
    let mapped_row = admitted
        .crosswalks()
        .iter()
        .find(|entry| entry.source() == &mapped_source)
        .expect("mapped source row is retained");
    assert!(mapped_row.targets().contains(&region_target));

    let portable = parse_machine_domain_portable_payload_v1(admitted.portable_payload())
        .expect("typed scenario portable payload decodes");
    assert_eq!(portable.scenario_hash(), admitted.scenario_hash());
    assert_eq!(
        portable.canonical_scenario_ir(),
        admitted.canonical_scenario_ir()
    );
    let decoded = fs_scenario::ir::parse_ir(portable.canonical_scenario_ir())
        .expect("versioned typed scenario reparses");
    assert_eq!(
        decoded.source_version(),
        fs_scenario::ir::SCENARIO_IR_VERSION
    );
    assert!(decoded.migration().is_none());
    assert_eq!(decoded.scenario(), &expected_scenario);
    assert!(matches!(
        &decoded.scenario().base_bcs[0].value,
        Some(BcValue::Typed(Payload::Vector(_)))
    ));

    let replay = admitted
        .verify_replay()
        .expect("typed artifact replays exactly");
    assert_eq!(replay.lowering(), admitted.identity());
    assert_eq!(replay.scenario_hash(), admitted.scenario_hash());
}

#[test]
#[allow(clippy::too_many_lines)]
fn g3_crosswalk_completeness_foreign_sources_and_targets_refuse() {
    let (graph, behavior, assurance) = machine_stack::admitted_machine_stack();
    let sources = required_machine_lowering_sources_v1(&graph, &behavior, &assurance)
        .expect("source inventory");

    let mut missing = complete_draft(
        &sources,
        scenario("missing-row"),
        "artifact-a",
        "law-a",
        "policy-a",
    );
    let omitted = missing.crosswalks.pop().expect("fixture crosswalk");
    let decision = with_cx(|cx| {
        admit_machine_domain_lowering_with_decision_v1(
            &graph,
            &behavior,
            &assurance,
            missing,
            ValidationBudget::default(),
            cx,
        )
    });
    assert_eq!(decision.code(), "MachineDomainLoweringRefused");
    assert_eq!(
        decision.refusal_code(),
        Some("MachineLoweringMissingCrosswalk")
    );
    assert_eq!(decision.submitted_counts().crosswalks(), sources.len() - 1);
    let refusal = decision.into_result().expect_err("missing source refuses");
    assert_eq!(refusal.code(), "MachineLoweringMissingCrosswalk");
    assert_eq!(
        refusal,
        MachineDomainLoweringRefusal::MissingCrosswalk {
            source: omitted.source().clone()
        }
    );

    let mut duplicate = complete_draft(
        &sources,
        scenario("duplicate-row"),
        "artifact-a",
        "law-a",
        "policy-a",
    );
    duplicate.crosswalks.push(duplicate.crosswalks[0].clone());
    let refusal = with_cx(|cx| {
        admit_machine_domain_lowering_v1(
            &graph,
            &behavior,
            &assurance,
            duplicate,
            ValidationBudget::default(),
            cx,
        )
        .expect_err("duplicate source refuses")
    });
    assert_eq!(refusal.code(), "MachineLoweringDuplicateCrosswalk");

    let mut foreign = complete_draft(
        &sources,
        scenario("foreign-row"),
        "artifact-a",
        "law-a",
        "policy-a",
    );
    foreign.crosswalks.push(
        MachineDomainCrosswalkEntryV1::new(
            MachineLoweringSourceV1::Clock(
                ClockId::new("clock/foreign").expect("foreign clock ID"),
            ),
            vec![MachineDomainTargetV1::Scenario(
                MachineScenarioLocatorV1::Root,
            )],
            law("law-a"),
        )
        .expect("foreign row shape"),
    );
    let refusal = with_cx(|cx| {
        admit_machine_domain_lowering_v1(
            &graph,
            &behavior,
            &assurance,
            foreign,
            ValidationBudget::default(),
            cx,
        )
        .expect_err("foreign source refuses")
    });
    assert_eq!(refusal.code(), "MachineLoweringUnexpectedCrosswalk");

    let mut bad_target = complete_draft(
        &sources,
        scenario("bad-target"),
        "artifact-a",
        "law-a",
        "policy-a",
    );
    bad_target.crosswalks[0] = MachineDomainCrosswalkEntryV1::new(
        bad_target.crosswalks[0].source().clone(),
        vec![MachineDomainTargetV1::Scenario(
            MachineScenarioLocatorV1::Region("missing-region".into()),
        )],
        law("law-a"),
    )
    .expect("bad locator row shape");
    let refusal = with_cx(|cx| {
        admit_machine_domain_lowering_v1(
            &graph,
            &behavior,
            &assurance,
            bad_target,
            ValidationBudget::default(),
            cx,
        )
        .expect_err("missing scenario locator refuses")
    });
    assert_eq!(refusal.code(), "MachineLoweringScenarioTarget");

    let mut empty_selector = complete_draft(
        &sources,
        scenario("empty-selector"),
        "artifact-a",
        "law-a",
        "policy-a",
    );
    let artifact_id = empty_selector.external_artifacts[0].id().clone();
    let behavior_row = empty_selector
        .crosswalks
        .iter()
        .position(|entry| matches!(entry.source(), MachineLoweringSourceV1::Behavior(_)))
        .expect("behavior crosswalk row");
    let source = empty_selector.crosswalks[behavior_row].source().clone();
    empty_selector.crosswalks[behavior_row] = MachineDomainCrosswalkEntryV1::new(
        source,
        vec![MachineDomainTargetV1::ExternalArtifact {
            artifact: artifact_id.clone(),
            selector: String::new().into_boxed_str(),
        }],
        law("law-a"),
    )
    .expect("direct empty selector reaches admission boundary");
    let refusal = with_cx(|cx| {
        admit_machine_domain_lowering_v1(
            &graph,
            &behavior,
            &assurance,
            empty_selector,
            ValidationBudget::default(),
            cx,
        )
        .expect_err("empty external selector refuses")
    });
    assert_eq!(
        refusal,
        MachineDomainLoweringRefusal::EmptyExternalSelector {
            artifact: artifact_id
        }
    );
}

#[test]
#[allow(clippy::too_many_lines)]
fn g3_invalid_scenario_or_orphan_external_artifact_publishes_no_identity() {
    let (graph, behavior, assurance) = machine_stack::admitted_machine_stack();
    let sources = required_machine_lowering_sources_v1(&graph, &behavior, &assurance)
        .expect("source inventory");

    let mut invalid = complete_draft(
        &sources,
        scenario("invalid-scenario"),
        "artifact-a",
        "law-a",
        "policy-a",
    );
    invalid.scenario.environment.ambient_temperature.value = -1.0;
    let refusal = with_cx(|cx| {
        admit_machine_domain_lowering_v1(
            &graph,
            &behavior,
            &assurance,
            invalid,
            ValidationBudget::default(),
            cx,
        )
        .expect_err("invalid scenario refuses")
    });
    assert_eq!(refusal.code(), "MachineLoweringScenarioFindings");

    let mut orphan = complete_draft(
        &sources,
        scenario("orphan-artifact"),
        "artifact-a",
        "law-a",
        "policy-a",
    );
    orphan.external_artifacts.push(
        MachineDomainArtifactRefV1::new(
            MachineDomainArtifactId::new("domain/orphan").expect("orphan ID"),
            MachineDomainArtifactKindV1::Time,
            "fixture/orphan",
            nz(1),
            digest("orphan"),
        )
        .expect("orphan artifact"),
    );
    let refusal = with_cx(|cx| {
        admit_machine_domain_lowering_v1(
            &graph,
            &behavior,
            &assurance,
            orphan,
            ValidationBudget::default(),
            cx,
        )
        .expect_err("orphan artifact refuses")
    });
    assert_eq!(refusal.code(), "MachineLoweringOrphanExternalArtifact");

    let mut excessive_records = scenario("excessive-records");
    excessive_records.cases = (0..=MAX_MACHINE_DOMAIN_SCENARIO_RECORDS)
        .map(|index| LoadCase {
            name: format!("case-{index}"),
            bcs: Vec::new(),
        })
        .collect();
    let refusal = with_cx(|cx| {
        admit_machine_domain_lowering_v1(
            &graph,
            &behavior,
            &assurance,
            complete_draft(
                &sources,
                excessive_records,
                "artifact-a",
                "law-a",
                "policy-a",
            ),
            ValidationBudget::default(),
            cx,
        )
        .expect_err("hard scenario-record bound refuses before validation scan")
    });
    assert!(matches!(
        refusal,
        MachineDomainLoweringRefusal::ResourceLimit {
            resource: "scenario structural records",
            requested,
            limit: MAX_MACHINE_DOMAIN_SCENARIO_RECORDS,
        } if requested == MAX_MACHINE_DOMAIN_SCENARIO_RECORDS + 1
    ));
}

#[test]
fn g5_permutation_is_stable_and_representative_identity_inputs_move_output() {
    let (graph, behavior, assurance) = machine_stack::admitted_machine_stack();
    let sources = required_machine_lowering_sources_v1(&graph, &behavior, &assurance)
        .expect("source inventory");
    let admit = |draft| {
        with_cx(|cx| {
            admit_machine_domain_lowering_v1(
                &graph,
                &behavior,
                &assurance,
                draft,
                ValidationBudget::default(),
                cx,
            )
            .expect("fixture projection admits")
        })
    };

    let baseline = admit(complete_draft(
        &sources,
        scenario("identity-fixture"),
        "artifact-a",
        "law-a",
        "policy-a",
    ));
    let mut permuted = complete_draft(
        &sources,
        scenario("identity-fixture"),
        "artifact-a",
        "law-a",
        "policy-a",
    );
    permuted.crosswalks.reverse();
    permuted.external_artifacts.reverse();
    let permuted = admit(permuted);
    assert_eq!(baseline.identity(), permuted.identity());
    assert_eq!(baseline.canonical_manifest(), permuted.canonical_manifest());
    assert_eq!(baseline.portable_payload(), permuted.portable_payload());

    let artifact_changed = admit(complete_draft(
        &sources,
        scenario("identity-fixture"),
        "artifact-b",
        "law-a",
        "policy-a",
    ));
    let law_changed = admit(complete_draft(
        &sources,
        scenario("identity-fixture"),
        "artifact-a",
        "law-b",
        "policy-a",
    ));
    let policy_changed = admit(complete_draft(
        &sources,
        scenario("identity-fixture"),
        "artifact-a",
        "law-a",
        "policy-b",
    ));
    let scenario_changed = admit(complete_draft(
        &sources,
        scenario("identity-fixture-changed"),
        "artifact-a",
        "law-a",
        "policy-a",
    ));
    for changed in [
        artifact_changed,
        law_changed,
        policy_changed,
        scenario_changed,
    ] {
        assert_ne!(baseline.identity(), changed.identity());
        assert_ne!(
            baseline.portable_payload_hash(),
            changed.portable_payload_hash()
        );
    }
}

#[test]
fn g4_pre_cancelled_projection_refuses_before_publication() {
    let (graph, behavior, assurance) = machine_stack::admitted_machine_stack();
    let sources = required_machine_lowering_sources_v1(&graph, &behavior, &assurance)
        .expect("source inventory");
    let refusal = with_cancelled_cx(|cx| {
        admit_machine_domain_lowering_v1(
            &graph,
            &behavior,
            &assurance,
            complete_draft(
                &sources,
                scenario("cancelled"),
                "artifact-a",
                "law-a",
                "policy-a",
            ),
            ValidationBudget::default(),
            cx,
        )
        .expect_err("pre-cancelled projection refuses")
    });
    assert_eq!(
        refusal,
        MachineDomainLoweringRefusal::Cancelled { phase: "initial" }
    );
}
