//! VerificationManifest v1 identity-core conformance (bead i94v.7.1.1).
//!
//! The battery locks the bead's test plan: canonicalization and digest
//! stability under input reordering, collision-freedom and Unicode
//! distinctness, supersession immutability, source-authority conflict,
//! post-receipt source mutation, relation orientation and variance,
//! favorable relabeling, cycles, certified-equivalence SCCs, migration
//! classification, the per-field registry, and semantic parity of the
//! human/JSON/ledger projections.

use fs_blake3::{ContentHash, hash_domain};
use fs_vmanifest::v1::{
    ClaimId, ClaimKind, ClaimRelationReceipt, ClaimRevision, MANIFEST_RECORD_FIELDS,
    MANIFEST_V1_SCHEMA_VERSION, MAX_V1_TEXT_BYTES, QuantifierVariance, RelationKind,
    SourceAuthority, SourcePin, V1AdmissionLimits, V1Error, admit_graph, admit_graph_with_limits,
    classify_migration, resolve_authority,
};
use std::collections::BTreeSet;

fn revision(lineage: &str, statement: &str) -> ClaimRevision {
    ClaimRevision {
        claim: ClaimId::new(lineage).expect("lineage id admits"),
        kind: ClaimKind::Behavioral,
        statement: statement.to_owned(),
        quantifiers: "for all admitted fixtures".to_owned(),
        units_conventions: "SI, lattice units declared per case".to_owned(),
        hypotheses: "committed tree, deterministic mode".to_owned(),
        domain: "fixture-scale battery domain".to_owned(),
        surface: "fs-example::api/CONTRACT#claims".to_owned(),
        no_claim: "no production-scale claim".to_owned(),
        supersedes: None,
    }
}

fn revision_id(revision: &ClaimRevision) -> ContentHash {
    revision.revision_id().expect("test revision admits")
}

fn push_reference_field(bytes: &mut Vec<u8>, value: &str) {
    bytes.extend_from_slice(
        &u32::try_from(value.len())
            .expect("valid v1 field fits u32")
            .to_be_bytes(),
    );
    bytes.extend_from_slice(value.as_bytes());
}

fn reference_revision_id(revision: &ClaimRevision) -> ContentHash {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&MANIFEST_V1_SCHEMA_VERSION.to_be_bytes());
    push_reference_field(&mut bytes, revision.claim.as_str());
    bytes.push(match revision.kind {
        ClaimKind::Behavioral => 1,
        ClaimKind::QuantitativeBound => 2,
        ClaimKind::Determinism => 3,
        ClaimKind::Refusal => 4,
        ClaimKind::Theorem => 5,
        _ => panic!("test reference must be updated for a new claim kind"),
    });
    for value in [
        revision.statement.as_str(),
        revision.quantifiers.as_str(),
        revision.units_conventions.as_str(),
        revision.hypotheses.as_str(),
        revision.domain.as_str(),
        revision.surface.as_str(),
        revision.no_claim.as_str(),
    ] {
        push_reference_field(&mut bytes, value);
    }
    match revision.supersedes {
        None => bytes.push(0),
        Some(previous) => {
            bytes.push(1);
            bytes.extend_from_slice(previous.as_bytes());
        }
    }
    hash_domain("org.frankensim.fs-vmanifest.claim-revision.v1", &bytes)
}

fn reference_graph_digest(graph: &fs_vmanifest::v1::NormalizedGraph) -> ContentHash {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&MANIFEST_V1_SCHEMA_VERSION.to_be_bytes());
    bytes.extend_from_slice(
        &u32::try_from(graph.revisions().len())
            .expect("valid graph revision count fits u32")
            .to_be_bytes(),
    );
    for revision in graph.revisions() {
        bytes.extend_from_slice(revision.as_bytes());
    }
    bytes.extend_from_slice(
        &u32::try_from(graph.edges().len())
            .expect("valid graph receipt count fits u32")
            .to_be_bytes(),
    );
    for edge in graph.edges() {
        bytes.push(match edge.kind {
            RelationKind::Implication => 1,
            RelationKind::Refinement => 2,
            RelationKind::Restriction => 3,
            RelationKind::Counterexample => 4,
            RelationKind::CertifiedEquivalence => 5,
            _ => panic!("test reference must be updated for a new relation kind"),
        });
        bytes.extend_from_slice(edge.from.as_bytes());
        bytes.extend_from_slice(edge.to.as_bytes());
        push_reference_field(&mut bytes, &edge.checker);
        push_reference_field(&mut bytes, &edge.tcb);
        bytes.push(match edge.variance {
            QuantifierVariance::Preserved => 1,
            QuantifierVariance::Weakened => 2,
            QuantifierVariance::Strengthened => 3,
            _ => panic!("test reference must be updated for new quantifier variance"),
        });
        push_reference_field(&mut bytes, &edge.domain_note);
        bytes.extend_from_slice(&edge.policy_version.to_be_bytes());
    }
    for (member, representative) in graph.representatives() {
        bytes.extend_from_slice(member.as_bytes());
        bytes.extend_from_slice(representative.as_bytes());
    }
    hash_domain("org.frankensim.fs-vmanifest.manifest-graph.v1", &bytes)
}

fn assert_resource_refusal(
    error: &V1Error,
    quantity: &'static str,
    required: u64,
    admitted: u64,
    unit: &'static str,
) {
    assert_eq!(error.rule(), "v1-resource-limit");
    let resource = error
        .resource_refusal()
        .expect("resource refusal carries typed evidence");
    assert_eq!(resource.quantity, quantity);
    assert_eq!(resource.required, required);
    assert_eq!(resource.admitted, admitted);
    assert_eq!(resource.unit, unit);
    assert!(error.detail().contains(&required.to_string()));
    assert!(error.detail().contains(&admitted.to_string()));
    assert!(error.detail().contains(unit));
}

fn edge(
    kind: RelationKind,
    from: &ClaimRevision,
    to: &ClaimRevision,
    variance: QuantifierVariance,
) -> ClaimRelationReceipt {
    ClaimRelationReceipt {
        kind,
        from: revision_id(from),
        to: revision_id(to),
        checker: "fs-checker/relation-v1".to_owned(),
        tcb: "rustc+fs-blake3".to_owned(),
        variance,
        domain_note: "identical domain".to_owned(),
        policy_version: 1,
    }
}

fn reversed(mut receipt: ClaimRelationReceipt) -> ClaimRelationReceipt {
    core::mem::swap(&mut receipt.from, &mut receipt.to);
    receipt
}

#[test]
fn identical_content_is_idempotent_and_distinct_content_cannot_collide() {
    let a = revision("claim/energy-balance", "energy balances to 1e-12");
    let b = revision("claim/energy-balance", "energy balances to 1e-12");
    assert_eq!(
        revision_id(&a),
        revision_id(&b),
        "content addressing is idempotent"
    );

    let c = revision("claim/energy-balance", "energy balances to 1e-11");
    assert_ne!(
        revision_id(&a),
        revision_id(&c),
        "one character = new revision"
    );

    // Unicode distinctness: visually confusable statements are distinct
    // content, never folded.
    let latin = revision("claim/units", "tolerance is 5 µm");
    let confusable = revision("claim/units", "tolerance is 5 μm");
    assert_ne!(revision_id(&latin), revision_id(&confusable));
}

#[test]
fn supersession_appends_and_never_mutates_the_old_revision() {
    let old = revision("claim/energy-balance", "energy balances to 1e-11");
    let old_id = revision_id(&old);
    let mut new = revision("claim/energy-balance", "energy balances to 1e-12");
    new.supersedes = Some(old_id);
    let new_id = revision_id(&new);
    assert_ne!(old_id, new_id);
    // The old revision's identity is untouched by being superseded.
    assert_eq!(revision_id(&old), old_id);
    // Supersession participates in identity: the same statement WITHOUT
    // the supersession pointer is a different revision.
    let mut orphan = new.clone();
    orphan.supersedes = None;
    assert_ne!(revision_id(&orphan), new_id);
}

#[test]
fn authoritative_revision_identity_is_fallible_and_stream_matches_reference() {
    for length in [1usize, 63, 64, 65, 1023, 1024, MAX_V1_TEXT_BYTES] {
        let candidate = revision("claim/streaming", &"x".repeat(length));
        assert_eq!(
            revision_id(&candidate),
            reference_revision_id(&candidate),
            "streamed identity preserves the legacy preimage at length {length}"
        );
    }

    let mut superseding = revision("claim/streaming", "successor");
    superseding.supersedes = Some(revision_id(&revision("claim/streaming", "predecessor")));
    assert_eq!(
        revision_id(&superseding),
        reference_revision_id(&superseding)
    );

    let mut empty = revision("claim/invalid", "valid before mutation");
    empty.statement.clear();
    let error = empty
        .revision_id()
        .expect_err("empty raw draft cannot mint authority");
    assert_eq!(error.rule(), "v1-field-bounds");

    let mut oversized = revision("claim/invalid", "valid before mutation");
    oversized.statement = "x".repeat(MAX_V1_TEXT_BYTES + 1);
    let error = oversized
        .revision_id()
        .expect_err("oversized raw draft cannot mint authority");
    assert_eq!(error.rule(), "v1-field-bounds");
}

#[test]
fn revision_text_limits_count_utf8_bytes_not_scalar_values() {
    let exact = revision("claim/utf8", &"é".repeat(MAX_V1_TEXT_BYTES / 2));
    assert_eq!(exact.statement.len(), MAX_V1_TEXT_BYTES);
    exact.revision_id().expect("exact UTF-8 byte cap admits");

    let mut over = exact;
    over.statement.push('x');
    let error = over
        .revision_id()
        .expect_err("one UTF-8 byte over the cap refuses");
    assert_eq!(error.rule(), "v1-field-bounds");
    assert!(
        error
            .detail()
            .contains(&(MAX_V1_TEXT_BYTES + 1).to_string())
    );
}

#[test]
fn graph_resource_limits_admit_zero_and_exact_then_refuse_one_more() {
    let zero_limits = V1AdmissionLimits {
        max_revisions: 0,
        max_receipts: 0,
        max_semantic_bytes: 0,
        max_working_bytes: 0,
        max_projection_rows: 1,
        max_projection_bytes: 1024,
    };
    let empty = admit_graph_with_limits(&[], &[], zero_limits).expect("zero-sized graph admits");
    let empty_resources = empty.resource_envelope();
    assert_eq!(empty_resources.revision_count(), 0);
    assert_eq!(empty_resources.receipt_count(), 0);
    assert_eq!(empty_resources.semantic_bytes(), 0);
    assert_eq!(empty_resources.logical_working_bytes(), 0);
    assert_eq!(empty_resources.projection_rows(), 1);
    assert_eq!(empty_resources.projection_bytes_upper_bound(), 1024);
    assert_eq!(empty.digest(), reference_graph_digest(&empty));

    let a = revision("claim/a", "a holds");
    let one = admit_graph(&[a.clone()], &[]).expect("one revision admits");
    let one_resources = one.resource_envelope();
    let exact_one = V1AdmissionLimits {
        max_revisions: 1,
        max_receipts: 0,
        max_semantic_bytes: one_resources.semantic_bytes(),
        max_working_bytes: one_resources.logical_working_bytes(),
        max_projection_rows: one_resources.projection_rows(),
        max_projection_bytes: one_resources.projection_bytes_upper_bound(),
    };
    let exact = admit_graph_with_limits(&[a.clone()], &[], exact_one)
        .expect("equality with every resource limit admits");
    assert_eq!(exact.digest(), one.digest());
    assert_eq!(exact.resource_envelope(), one_resources);

    let b = revision("claim/b", "b holds");
    let error = admit_graph_with_limits(&[a.clone(), b.clone()], &[], exact_one)
        .expect_err("revision limit plus one refuses");
    assert_resource_refusal(&error, "revision count", 2, 1, "revisions");

    let receipt = edge(
        RelationKind::Implication,
        &a,
        &b,
        QuantifierVariance::Preserved,
    );
    let graph =
        admit_graph(&[a.clone(), b.clone()], &[receipt.clone()]).expect("one-receipt graph admits");
    let resources = graph.resource_envelope();
    let exact_limits = V1AdmissionLimits {
        max_revisions: resources.revision_count(),
        max_receipts: resources.receipt_count(),
        max_semantic_bytes: resources.semantic_bytes(),
        max_working_bytes: resources.logical_working_bytes(),
        max_projection_rows: resources.projection_rows(),
        max_projection_bytes: resources.projection_bytes_upper_bound(),
    };
    let exact = admit_graph_with_limits(&[a.clone(), b.clone()], &[receipt.clone()], exact_limits)
        .expect("complete exact graph envelope admits");
    assert_eq!(exact.digest(), graph.digest());
    assert_eq!(exact.resource_envelope(), resources);

    let mut parallel = receipt.clone();
    parallel.policy_version += 1;
    let error =
        admit_graph_with_limits(&[a.clone(), b.clone()], &[receipt, parallel], exact_limits)
            .expect_err("receipt limit plus one refuses");
    assert_resource_refusal(&error, "receipt count", 2, 1, "receipts");

    for (quantity, limits, required, admitted, unit) in [
        (
            "aggregate semantic bytes",
            V1AdmissionLimits {
                max_semantic_bytes: resources.semantic_bytes() - 1,
                ..exact_limits
            },
            resources.semantic_bytes(),
            resources.semantic_bytes() - 1,
            "UTF-8 bytes",
        ),
        (
            "logical graph working bytes",
            V1AdmissionLimits {
                max_working_bytes: resources.logical_working_bytes() - 1,
                ..exact_limits
            },
            resources.logical_working_bytes(),
            resources.logical_working_bytes() - 1,
            "logical bytes",
        ),
        (
            "projection row count",
            V1AdmissionLimits {
                max_projection_rows: resources.projection_rows() - 1,
                ..exact_limits
            },
            u64::from(resources.projection_rows()),
            u64::from(resources.projection_rows() - 1),
            "rows",
        ),
        (
            "projection byte upper bound",
            V1AdmissionLimits {
                max_projection_bytes: resources.projection_bytes_upper_bound() - 1,
                ..exact_limits
            },
            resources.projection_bytes_upper_bound(),
            resources.projection_bytes_upper_bound() - 1,
            "logical bytes",
        ),
    ] {
        let error = admit_graph_with_limits(
            &[a.clone(), b.clone()],
            &[edge(
                RelationKind::Implication,
                &a,
                &b,
                QuantifierVariance::Preserved,
            )],
            limits,
        )
        .expect_err("one unit below the required envelope refuses");
        assert_resource_refusal(&error, quantity, required, admitted, unit);
    }
}

#[test]
fn resource_preflight_precedes_hostile_payload_work_and_cannot_be_loosened() {
    let mut malformed = revision("claim/malformed", "initially valid");
    malformed.statement.clear();
    let other = revision("claim/other", "other holds");
    let error = admit_graph_with_limits(
        &[malformed, other],
        &[],
        V1AdmissionLimits {
            max_revisions: 1,
            ..V1AdmissionLimits::DEFAULT
        },
    )
    .expect_err("count refusal wins before malformed element scan");
    assert_resource_refusal(&error, "revision count", 2, 1, "revisions");

    let a = revision("claim/a", "a holds");
    let b = revision("claim/b", "b holds");
    let mut hostile = edge(
        RelationKind::Implication,
        &a,
        &b,
        QuantifierVariance::Preserved,
    );
    hostile.checker = "x".repeat(MAX_V1_TEXT_BYTES + 1);
    let mut second = hostile.clone();
    second.policy_version += 1;
    let error = admit_graph_with_limits(
        &[a.clone(), b.clone()],
        &[hostile.clone(), second],
        V1AdmissionLimits {
            max_receipts: 1,
            ..V1AdmissionLimits::DEFAULT
        },
    )
    .expect_err("receipt count refuses before hostile text scan or clone");
    assert_resource_refusal(&error, "receipt count", 2, 1, "receipts");

    let error = admit_graph_with_limits(
        &[a, b],
        &[hostile],
        V1AdmissionLimits {
            max_semantic_bytes: 1,
            ..V1AdmissionLimits::DEFAULT
        },
    )
    .expect_err("aggregate bytes refuse before per-field semantics");
    assert_eq!(error.rule(), "v1-resource-limit");
    assert_eq!(
        error.resource_refusal().expect("typed evidence").quantity,
        "aggregate semantic bytes"
    );

    let error = admit_graph_with_limits(
        &[],
        &[],
        V1AdmissionLimits {
            max_revisions: V1AdmissionLimits::DEFAULT.max_revisions + 1,
            ..V1AdmissionLimits::DEFAULT
        },
    )
    .expect_err("local policy cannot loosen the v1 protocol ceiling");
    assert_resource_refusal(
        &error,
        "configured revision ceiling",
        u64::from(V1AdmissionLimits::DEFAULT.max_revisions) + 1,
        u64::from(V1AdmissionLimits::DEFAULT.max_revisions),
        "revisions",
    );
}

#[test]
fn streamed_graph_digest_matches_reference_and_projection_envelope() {
    let a = revision("claim/a", "a holds");
    let b = revision("claim/b", "b holds");
    let c = revision("claim/c", "c holds");
    let mut implication = edge(
        RelationKind::Implication,
        &a,
        &b,
        QuantifierVariance::Preserved,
    );
    implication.checker = "checker\"\\\ncontrol".to_owned();
    let refinement = edge(
        RelationKind::Refinement,
        &b,
        &c,
        QuantifierVariance::Weakened,
    );
    let graph = admit_graph(&[c, a, b], &[refinement, implication]).expect("seed graph admits");

    assert_eq!(graph.digest(), reference_graph_digest(&graph));
    assert_eq!(graph.digest(), graph.digest(), "cached identity is stable");
    let resources = graph.resource_envelope();
    let human = graph.render_human();
    let json = graph.render_json_lines();
    let ledger = graph.render_ledger_rows();
    let json_bytes = json.iter().map(String::len).sum::<usize>() + json.len().saturating_sub(1);
    let ledger_bytes = ledger
        .iter()
        .map(|(kind, payload)| kind.len() + payload.len() + 1)
        .sum::<usize>();
    for actual in [human.len(), json_bytes, ledger_bytes] {
        assert!(
            u64::try_from(actual).expect("projection length fits u64")
                <= resources.projection_bytes_upper_bound(),
            "actual projection {actual} exceeds retained upper bound {}",
            resources.projection_bytes_upper_bound()
        );
    }
    assert_eq!(
        resources.projection_rows(),
        1 + resources.revision_count() + resources.receipt_count()
    );
}

#[test]
fn normalization_digest_is_stable_under_input_reordering() {
    let a = revision("claim/a", "a holds");
    let b = revision("claim/b", "b holds");
    let c = revision("claim/c", "c holds");
    let e1 = edge(
        RelationKind::Implication,
        &a,
        &b,
        QuantifierVariance::Preserved,
    );
    let e2 = edge(
        RelationKind::Refinement,
        &b,
        &c,
        QuantifierVariance::Weakened,
    );

    let g1 = admit_graph(
        &[a.clone(), b.clone(), c.clone()],
        &[e1.clone(), e2.clone()],
    )
    .expect("graph admits");
    let g2 = admit_graph(&[c, a, b], &[e2, e1]).expect("reordered graph admits");
    assert_eq!(
        g1.digest(),
        g2.digest(),
        "equivalent semantic manifests normalize to the same digest"
    );
    assert_eq!(g1.resource_envelope(), g2.resource_envelope());
    assert_eq!(g1.digest(), reference_graph_digest(&g1));
}

#[test]
fn certified_equivalence_normalizes_endpoint_orientation_everywhere() {
    let a = revision("claim/a", "a holds");
    let b = revision("claim/b", "b holds");
    let receipt = edge(
        RelationKind::CertifiedEquivalence,
        &a,
        &b,
        QuantifierVariance::Preserved,
    );
    let reverse = reversed(receipt.clone());
    let revisions = [a, b];
    let forward = admit_graph(&revisions, &[receipt]).expect("forward equivalence admits");
    let backward = admit_graph(&revisions, &[reverse]).expect("reverse equivalence admits");
    let minimum = revision_id(&revisions[0]).min(revision_id(&revisions[1]));
    let maximum = revision_id(&revisions[0]).max(revision_id(&revisions[1]));

    assert_eq!(forward.edges(), backward.edges());
    assert_eq!(forward.edges()[0].from, minimum);
    assert_eq!(forward.edges()[0].to, maximum);
    assert_eq!(forward.representatives(), backward.representatives());
    assert_eq!(forward.digest(), backward.digest());
    let human = forward.render_human();
    assert!(human.contains("<->"), "equivalence has no directed arrow");
    assert!(human.contains(&format!(
        "from={} <-> to={}",
        minimum.to_hex(),
        maximum.to_hex()
    )));
    assert_eq!(human, backward.render_human());
    let json = forward.render_json_lines();
    assert!(json.last().expect("edge row").contains(&format!(
        "\"from\":\"{}\",\"to\":\"{}\"",
        minimum.to_hex(),
        maximum.to_hex()
    )));
    assert_eq!(json, backward.render_json_lines());
    let ledger = forward.render_ledger_rows();
    assert!(ledger.last().expect("edge row").1.contains(&format!(
        "from={} to={}",
        minimum.to_hex(),
        maximum.to_hex()
    )));
    assert_eq!(ledger, backward.render_ledger_rows());
}

#[test]
fn certified_equivalence_normalizes_before_dangling_endpoint_diagnostics() {
    let a = revision("claim/a", "a holds");
    let b = revision("claim/b", "b holds");
    let receipt = edge(
        RelationKind::CertifiedEquivalence,
        &a,
        &b,
        QuantifierVariance::Preserved,
    );
    let forward = admit_graph(&[], &[receipt.clone()]).expect_err("both endpoints dangle");
    let backward = admit_graph(&[], &[reversed(receipt)]).expect_err("both endpoints dangle");

    assert_eq!(forward, backward, "orientation cannot leak through refusal");
    assert_eq!(forward.rule(), "v1-dangling-relation");
    assert!(
        forward
            .detail()
            .contains(&revision_id(&a).min(revision_id(&b)).to_hex())
    );
}

#[test]
fn certified_equivalence_refuses_directional_quantifier_variance() {
    let a = revision("claim/a", "a holds");
    let b = revision("claim/b", "b holds");
    for variance in [
        QuantifierVariance::Weakened,
        QuantifierVariance::Strengthened,
    ] {
        let receipt = edge(RelationKind::CertifiedEquivalence, &a, &b, variance);
        for input in [receipt.clone(), reversed(receipt)] {
            let error = admit_graph(&[a.clone(), b.clone()], &[input])
                .expect_err("directional variance cannot be normalized as equivalence");
            assert_eq!(error.rule(), "v1-equivalence-variance");
            assert!(error.detail().contains("requires preserved quantifiers"));
            assert!(error.ranked_fixes()[0].contains("directed relation kind"));
        }
    }
}

#[test]
fn reverse_equivalence_duplicates_refuse_after_normalization() {
    let a = revision("claim/a", "a holds");
    let b = revision("claim/b", "b holds");
    let receipt = edge(
        RelationKind::CertifiedEquivalence,
        &a,
        &b,
        QuantifierVariance::Preserved,
    );
    let error = admit_graph(&[a, b], &[receipt.clone(), reversed(receipt)])
        .expect_err("reverse duplicate must not mint additional authority");
    assert_eq!(error.rule(), "v1-duplicate-relation");
    assert!(error.detail().contains("canonical relation receipt"));
    assert!(error.ranked_fixes()[0].contains("reversing"));
}

#[test]
fn parallel_equivalence_evidence_remains_distinct_after_orientation_normalization() {
    let a = revision("claim/a", "a holds");
    let b = revision("claim/b", "b holds");
    let base = edge(
        RelationKind::CertifiedEquivalence,
        &a,
        &b,
        QuantifierVariance::Preserved,
    );
    let mut checker = reversed(base.clone());
    checker.checker = "fs-checker/relation-v2".to_owned();
    let mut tcb = base.clone();
    tcb.tcb = "rustc+fs-blake3+fs-ivl".to_owned();
    let mut policy = reversed(base.clone());
    policy.policy_version = 2;
    let mut domain = reversed(base.clone());
    domain.domain_note = "same domain certified by alternate partition".to_owned();
    let graph = admit_graph(&[a, b], &[policy, domain, tcb, checker, base])
        .expect("independent equivalence evidence remains distinct");

    assert_eq!(graph.edges().len(), 5);
    let endpoint = (graph.edges()[0].from, graph.edges()[0].to);
    assert!(
        graph
            .edges()
            .iter()
            .all(|receipt| (receipt.from, receipt.to) == endpoint)
    );
    assert_eq!(
        graph
            .edges()
            .iter()
            .map(|receipt| receipt.checker.as_str())
            .collect::<BTreeSet<_>>()
            .len(),
        2
    );
    assert_eq!(
        graph
            .edges()
            .iter()
            .map(|receipt| receipt.tcb.as_str())
            .collect::<BTreeSet<_>>()
            .len(),
        2
    );
    assert_eq!(
        graph
            .edges()
            .iter()
            .map(|receipt| receipt.policy_version)
            .collect::<BTreeSet<_>>()
            .len(),
        2
    );
    assert_eq!(
        graph
            .edges()
            .iter()
            .map(|receipt| receipt.domain_note.as_str())
            .collect::<BTreeSet<_>>()
            .len(),
        2
    );
    assert!(
        graph
            .edges()
            .iter()
            .all(|receipt| receipt.variance == QuantifierVariance::Preserved)
    );
}

#[test]
fn directed_relation_kinds_preserve_endpoint_orientation() {
    let a = revision("claim/a", "a holds");
    let b = revision("claim/b", "b holds");
    for kind in [
        RelationKind::Implication,
        RelationKind::Refinement,
        RelationKind::Restriction,
        RelationKind::Counterexample,
    ] {
        let receipt = edge(kind, &a, &b, QuantifierVariance::Preserved);
        let reverse = reversed(receipt.clone());
        let forward = admit_graph(&[a.clone(), b.clone()], &[receipt.clone()])
            .expect("directed relation admits");
        let backward = admit_graph(&[a.clone(), b.clone()], &[reverse.clone()])
            .expect("reverse directed relation admits");
        assert_eq!(forward.edges()[0].from, receipt.from);
        assert_eq!(forward.edges()[0].to, receipt.to);
        assert_eq!(backward.edges()[0].from, reverse.from);
        assert_eq!(backward.edges()[0].to, reverse.to);
        let forward_pair = format!(
            "from={} -> to={}",
            receipt.from.to_hex(),
            receipt.to.to_hex()
        );
        let backward_pair = format!(
            "from={} -> to={}",
            reverse.from.to_hex(),
            reverse.to.to_hex()
        );
        assert!(forward.render_human().contains(&forward_pair));
        assert!(backward.render_human().contains(&backward_pair));
        assert!(
            forward
                .render_json_lines()
                .last()
                .expect("edge row")
                .contains(&format!(
                    "\"from\":\"{}\",\"to\":\"{}\"",
                    receipt.from.to_hex(),
                    receipt.to.to_hex()
                ))
        );
        assert!(
            backward
                .render_json_lines()
                .last()
                .expect("edge row")
                .contains(&format!(
                    "\"from\":\"{}\",\"to\":\"{}\"",
                    reverse.from.to_hex(),
                    reverse.to.to_hex()
                ))
        );
        assert!(
            forward
                .render_ledger_rows()
                .last()
                .expect("edge row")
                .1
                .contains(&format!(
                    "from={} to={}",
                    receipt.from.to_hex(),
                    receipt.to.to_hex()
                ))
        );
        assert!(
            backward
                .render_ledger_rows()
                .last()
                .expect("edge row")
                .1
                .contains(&format!(
                    "from={} to={}",
                    reverse.from.to_hex(),
                    reverse.to.to_hex()
                ))
        );
        assert_ne!(
            forward.digest(),
            backward.digest(),
            "{kind:?} orientation remains identity-forming"
        );
        assert_ne!(forward.render_human(), backward.render_human());
        assert_ne!(forward.render_json_lines(), backward.render_json_lines());
        assert_ne!(forward.render_ledger_rows(), backward.render_ledger_rows());
    }
}

#[test]
fn admission_owns_and_normalizes_mutable_equivalence_drafts() {
    let a = revision("claim/a", "a holds");
    let b = revision("claim/b", "b holds");
    let mut draft = edge(
        RelationKind::Implication,
        &a,
        &b,
        QuantifierVariance::Preserved,
    );
    draft.kind = RelationKind::CertifiedEquivalence;
    if draft.from < draft.to {
        core::mem::swap(&mut draft.from, &mut draft.to);
    }
    let submitted = draft.clone();
    let graph = admit_graph(&[a, b], &[draft.clone()]).expect("kind-flipped draft admits");
    let admitted = graph.edges()[0].clone();

    assert_eq!(draft, submitted, "admission does not mutate caller state");
    assert!(admitted.from < admitted.to, "admitted endpoints normalize");

    draft.kind = RelationKind::Counterexample;
    core::mem::swap(&mut draft.from, &mut draft.to);
    draft.checker = "mutated-after-admission".to_owned();
    assert_eq!(
        graph.edges()[0],
        admitted,
        "later draft mutation cannot alter the admitted graph"
    );
}

#[test]
fn equivalence_subset_reversal_and_permutation_are_full_projection_invariants() {
    let a = revision("claim/a", "a holds");
    let b = revision("claim/b", "b holds");
    let c = revision("claim/c", "c holds");
    let receipts = [
        edge(
            RelationKind::CertifiedEquivalence,
            &a,
            &b,
            QuantifierVariance::Preserved,
        ),
        edge(
            RelationKind::CertifiedEquivalence,
            &b,
            &c,
            QuantifierVariance::Preserved,
        ),
        edge(
            RelationKind::CertifiedEquivalence,
            &a,
            &c,
            QuantifierVariance::Preserved,
        ),
    ];
    let revisions = [a, b, c];
    let reference = admit_graph(&revisions, &receipts).expect("reference graph admits");
    const PERMUTATIONS: [[usize; 3]; 6] = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];

    for reversed_mask in 0_u8..8 {
        let oriented: Vec<_> = receipts
            .iter()
            .cloned()
            .enumerate()
            .map(|(index, receipt)| {
                if reversed_mask & (1 << index) == 0 {
                    receipt
                } else {
                    reversed(receipt)
                }
            })
            .collect();
        for permutation in PERMUTATIONS {
            let input: Vec<_> = permutation
                .iter()
                .map(|&index| oriented[index].clone())
                .collect();
            let graph = admit_graph(
                &[
                    revisions[2].clone(),
                    revisions[0].clone(),
                    revisions[1].clone(),
                ],
                &input,
            )
            .expect("reoriented and permuted graph admits");
            assert_eq!(graph.edges(), reference.edges());
            assert_eq!(graph.representatives(), reference.representatives());
            assert_eq!(graph.resource_envelope(), reference.resource_envelope());
            assert_eq!(graph.digest(), reference.digest());
            assert_eq!(graph.render_human(), reference.render_human());
            assert_eq!(graph.render_json_lines(), reference.render_json_lines());
            assert_eq!(graph.render_ledger_rows(), reference.render_ledger_rows());
        }
    }
}

#[test]
fn favorable_relabeling_severs_nothing_silently() {
    // The attack the Background names: rename the surface, keep the old
    // receipt. The renamed claim is a NEW revision; the old receipt now
    // points at a revision absent from the manifest and the graph
    // refuses with ranked fixes instead of silently rebinding.
    let original = revision("claim/roof", "roofline attainment >= 0.4");
    let target = revision("claim/floor", "floor holds");
    let receipt = edge(
        RelationKind::Implication,
        &original,
        &target,
        QuantifierVariance::Preserved,
    );
    let mut relabeled = original.clone();
    relabeled.surface = "fs-example::renamed/CONTRACT#claims".to_owned();

    let err = admit_graph(&[relabeled, target], &[receipt]).expect_err("severed");
    assert_eq!(err.rule(), "v1-dangling-relation");
    assert!(
        !err.ranked_fixes().is_empty(),
        "diagnostics carry ranked fixes"
    );
}

#[test]
fn promotion_transfer_respects_kind_and_quantifier_variance() {
    let a = revision("claim/a", "a holds");
    let b = revision("claim/b", "b holds");
    let implication = edge(
        RelationKind::Implication,
        &a,
        &b,
        QuantifierVariance::Preserved,
    );
    assert!(implication.promotion_transfers());
    let weakened = edge(
        RelationKind::Refinement,
        &a,
        &b,
        QuantifierVariance::Weakened,
    );
    assert!(weakened.promotion_transfers());
    let counterexample = edge(
        RelationKind::Counterexample,
        &a,
        &b,
        QuantifierVariance::Preserved,
    );
    assert!(
        !counterexample.promotion_transfers(),
        "counterexamples transfer refutation, never promotion"
    );
    let strengthened = edge(
        RelationKind::Implication,
        &a,
        &b,
        QuantifierVariance::Strengthened,
    );
    assert!(
        !strengthened.promotion_transfers(),
        "quantifier strengthening along an edge is unsound for transfer"
    );
}

#[test]
fn directed_cycles_refuse_unless_certified_equivalent() {
    let a = revision("claim/a", "a holds");
    let b = revision("claim/b", "b holds");
    let ab = edge(
        RelationKind::Implication,
        &a,
        &b,
        QuantifierVariance::Preserved,
    );
    let ba = edge(
        RelationKind::Implication,
        &b,
        &a,
        QuantifierVariance::Preserved,
    );
    let err = admit_graph(&[a.clone(), b.clone()], &[ab.clone(), ba]).expect_err("cycle refuses");
    assert_eq!(err.rule(), "v1-relation-cycle");
    assert!(err.ranked_fixes().iter().any(|f| f.contains("equivalence")));

    // The same cycle under certified equivalence canonicalizes instead.
    let eq = edge(
        RelationKind::CertifiedEquivalence,
        &a,
        &b,
        QuantifierVariance::Preserved,
    );
    let graph = admit_graph(&[a.clone(), b.clone()], &[ab, eq]).expect("SCC admits");
    let ida = revision_id(&a);
    let idb = revision_id(&b);
    let representative = ida.min(idb);
    assert_eq!(graph.representative_of(&ida), Some(representative));
    assert_eq!(graph.representative_of(&idb), Some(representative));
    assert_eq!(
        graph.revisions().len(),
        2,
        "canonicalization never erases members"
    );
}

#[test]
fn counterexample_inside_an_equivalence_component_is_a_contradiction() {
    let a = revision("claim/a", "a holds");
    let b = revision("claim/b", "b holds");
    let eq = edge(
        RelationKind::CertifiedEquivalence,
        &a,
        &b,
        QuantifierVariance::Preserved,
    );
    let cx = edge(
        RelationKind::Counterexample,
        &a,
        &b,
        QuantifierVariance::Preserved,
    );
    let err = admit_graph(&[a, b], &[eq, cx]).expect_err("contradiction");
    assert_eq!(err.rule(), "v1-domain-contradiction");
}

#[test]
fn source_authority_conflicts_resolve_upward_or_refuse_with_fixes() {
    let snap = |byte: u8| fs_blake3::hash_domain("test-snapshot", &[byte]);
    let contract = SourcePin {
        source: "crates/fs-example/CONTRACT.md".to_owned(),
        authority: SourceAuthority::Contract,
        snapshot: snap(1),
    };
    let mutated = SourcePin {
        source: "crates/fs-example/CONTRACT.md".to_owned(),
        authority: SourceAuthority::FrozenSnapshot,
        snapshot: snap(2),
    };
    // Post-receipt source mutation: the frozen snapshot wins by minting
    // against ITS hash — the lower pin is never silently reinterpreted.
    let winner = resolve_authority(&contract, &mutated).expect("lattice resolves");
    assert_eq!(winner.snapshot, snap(2));

    let peer = SourcePin {
        source: "crates/fs-example/CONTRACT.md".to_owned(),
        authority: SourceAuthority::Contract,
        snapshot: snap(3),
    };
    let err = resolve_authority(&contract, &peer).expect_err("equal authority refuses");
    assert_eq!(err.rule(), "v1-authority-conflict");
    assert_eq!(err.ranked_fixes().len(), 2);

    assert!(SourceAuthority::FrozenSnapshot > SourceAuthority::BeadObligation);
    assert!(SourceAuthority::BeadObligation > SourceAuthority::Contract);
    assert!(SourceAuthority::Contract > SourceAuthority::TestSource);
    assert!(SourceAuthority::TestSource > SourceAuthority::GeneratedArtifact);
}

#[test]
fn migration_classification_is_additive_or_reported_lossy_never_silent() {
    let old: BTreeSet<String> = ["a", "b", "c"].iter().map(|s| (*s).to_owned()).collect();
    let grown: BTreeSet<String> = ["a", "b", "c", "d"]
        .iter()
        .map(|s| (*s).to_owned())
        .collect();
    let shrunk: BTreeSet<String> = ["a", "b"].iter().map(|s| (*s).to_owned()).collect();

    match classify_migration(&old, &grown, None).expect("additive admits") {
        fs_vmanifest::v1::Migration::Additive { added } => assert_eq!(added, vec!["d"]),
        other => panic!("additive misclassified: {other:?}"),
    }
    let err = classify_migration(&old, &shrunk, None).expect_err("silent loss refuses");
    assert_eq!(err.rule(), "v1-lossy-migration-undeclared");
    match classify_migration(&old, &shrunk, Some("field c merged into b per policy 2"))
        .expect("declared lossy admits")
    {
        fs_vmanifest::v1::Migration::Breaking { report } => {
            assert_eq!(report.dropped, vec!["c"]);
            assert!(report.reason.contains("policy 2"));
        }
        other => panic!("breaking misclassified: {other:?}"),
    }
}

#[test]
fn the_field_registry_gives_every_field_full_metadata() {
    assert_eq!(
        MANIFEST_RECORD_FIELDS.len(),
        22,
        "the v1 record binds exactly the declared field families"
    );
    let mut names = BTreeSet::new();
    for field in MANIFEST_RECORD_FIELDS {
        assert!(names.insert(field.name), "duplicate field {}", field.name);
        assert!(!field.units.is_empty());
        assert!(
            ["1", "0..1", "0..n", "1..n"].contains(&field.cardinality),
            "{} has malformed cardinality {}",
            field.name,
            field.cardinality
        );
        assert!(
            ["additive", "breaking-if-removed", "identity-forming"].contains(&field.migration),
            "{} has malformed migration semantics",
            field.name
        );
    }
    assert!(
        MANIFEST_RECORD_FIELDS
            .iter()
            .any(|f| f.migration == "identity-forming" && f.name == "claim-revision"),
        "the claim revision is identity-forming by definition"
    );
}

#[test]
fn projections_are_semantic_views_of_one_digest() {
    let a = revision("claim/a", "a holds");
    let b = revision("claim/b", "b holds");
    let graph = admit_graph(
        &[a.clone(), b.clone()],
        &[edge(
            RelationKind::Implication,
            &a,
            &b,
            QuantifierVariance::Preserved,
        )],
    )
    .expect("graph admits");
    let digest = graph.digest().to_hex();

    let human = graph.render_human();
    let json = graph.render_json_lines();
    let ledger = graph.render_ledger_rows();
    assert!(human.contains(digest.as_str()));
    assert!(json.iter().all(|record| record.contains(digest.as_str())));
    assert!(
        ledger
            .iter()
            .all(|(_, payload)| payload.contains(digest.as_str()))
    );
    // Header + every revision + every edge exist in the machine
    // projections. The human view likewise carries every revision and edge.
    assert_eq!(
        json.len(),
        1 + graph.revisions().len() + graph.edges().len()
    );
    assert_eq!(ledger.len(), json.len());
    assert_eq!(human.matches("  revision ").count(), 2);
    assert_eq!(human.matches("->").count(), 1);
    assert_ne!(
        human,
        json.join("\n"),
        "projections need not be byte-identical"
    );
}

#[test]
fn json_lines_are_an_exact_complete_canonical_projection() {
    let a = revision("claim/a", "a holds");
    let b = revision("claim/b", "b holds");
    let mut receipt = edge(
        RelationKind::Implication,
        &a,
        &b,
        QuantifierVariance::Weakened,
    );
    receipt.checker = "checker\"\\v2".to_owned();
    receipt.tcb = "rustc\nfs-blake3".to_owned();
    receipt.domain_note = "strict\tdomain".to_owned();
    receipt.policy_version = 17;
    let from = receipt.from.to_hex();
    let to = receipt.to.to_hex();
    let mut revision_ids = [revision_id(&a), revision_id(&b)];
    revision_ids.sort();
    let graph = admit_graph(&[a, b], &[receipt]).expect("graph admits");
    let digest = graph.digest().to_hex();

    let json = graph.render_json_lines();
    let mut expected = vec![format!(
        "{{\"vmanifest\":\"graph\",\"schema_version\":{},\"graph_digest\":\"{}\",\"ordinal\":0,\"revisions\":2,\"edges\":1}}",
        MANIFEST_V1_SCHEMA_VERSION, digest,
    )];
    for (index, revision) in revision_ids.iter().enumerate() {
        expected.push(format!(
            "{{\"vmanifest\":\"revision\",\"schema_version\":{},\"graph_digest\":\"{}\",\"ordinal\":{},\"revision\":\"{}\",\"representative\":\"{}\"}}",
            MANIFEST_V1_SCHEMA_VERSION,
            digest,
            index + 1,
            revision.to_hex(),
            revision.to_hex(),
        ));
    }
    expected.push(format!(
        "{{\"vmanifest\":\"edge\",\"schema_version\":{},\"graph_digest\":\"{}\",\"ordinal\":3,\"kind\":\"Implication\",\"from\":\"{}\",\"to\":\"{}\",\"checker\":\"checker\\\"\\\\v2\",\"tcb\":\"rustc\\nfs-blake3\",\"variance\":\"Weakened\",\"domain_note\":\"strict\\tdomain\",\"policy_version\":17}}",
        MANIFEST_V1_SCHEMA_VERSION,
        digest,
        from,
        to,
    ));
    assert_eq!(
        json, expected,
        "strict projection includes every normalized field in canonical order"
    );

    let human = graph.render_human();
    for required in [
        format!("schema_version={MANIFEST_V1_SCHEMA_VERSION}"),
        format!("graph_digest={digest}"),
        format!("revision={}", revision_ids[0].to_hex()),
        format!("revision={}", revision_ids[1].to_hex()),
        "kind=Implication".to_owned(),
        format!("from={from}"),
        format!("to={to}"),
        "checker=\"checker\\\"\\\\v2\"".to_owned(),
        "tcb=\"rustc\\nfs-blake3\"".to_owned(),
        "variance=Weakened".to_owned(),
        "domain_note=\"strict\\tdomain\"".to_owned(),
        "policy_version=17".to_owned(),
    ] {
        assert!(
            human.contains(required.as_str()),
            "human projection omitted {required}"
        );
    }
}

fn assert_receipt_field_is_bound_and_canonical(
    field: &str,
    revisions: &[ClaimRevision],
    left: ClaimRelationReceipt,
    right: ClaimRelationReceipt,
) {
    let left_only = admit_graph(revisions, &[left.clone()]).expect("left receipt admits");
    let right_only = admit_graph(revisions, &[right.clone()]).expect("right receipt admits");
    assert_ne!(
        left_only.digest(),
        right_only.digest(),
        "changing {field} must change graph identity"
    );

    let forward = admit_graph(revisions, &[left.clone(), right.clone()]).expect("pair admits");
    let reverse = admit_graph(revisions, &[right, left]).expect("permuted pair admits");
    assert_eq!(
        forward.edges(),
        reverse.edges(),
        "{field} participates in the canonical total order"
    );
    assert_eq!(
        forward.digest(),
        reverse.digest(),
        "receipt input order cannot change identity when {field} breaks the tie"
    );
    assert_eq!(forward.render_json_lines(), reverse.render_json_lines());
}

#[test]
fn every_receipt_identity_field_is_bound_and_breaks_canonical_ties() {
    let a = revision("claim/a", "a holds");
    let b = revision("claim/b", "b holds");
    let c = revision("claim/c", "c holds");
    let revisions = [a.clone(), b.clone(), c.clone()];
    let base = edge(
        RelationKind::Implication,
        &a,
        &b,
        QuantifierVariance::Preserved,
    );

    let mut changed = base.clone();
    changed.kind = RelationKind::Refinement;
    assert_receipt_field_is_bound_and_canonical("kind", &revisions, base.clone(), changed);

    let mut changed = base.clone();
    changed.from = revision_id(&c);
    assert_receipt_field_is_bound_and_canonical("from", &revisions, base.clone(), changed);

    let mut changed = base.clone();
    changed.to = revision_id(&c);
    assert_receipt_field_is_bound_and_canonical("to", &revisions, base.clone(), changed);

    let mut changed = base.clone();
    changed.checker = "fs-checker/relation-v2".to_owned();
    assert_receipt_field_is_bound_and_canonical("checker", &revisions, base.clone(), changed);

    let mut changed = base.clone();
    changed.tcb = "rustc+fs-blake3+fs-ivl".to_owned();
    assert_receipt_field_is_bound_and_canonical("tcb", &revisions, base.clone(), changed);

    let mut changed = base.clone();
    changed.variance = QuantifierVariance::Weakened;
    assert_receipt_field_is_bound_and_canonical("variance", &revisions, base.clone(), changed);

    let mut changed = base.clone();
    changed.domain_note = "strict subdomain".to_owned();
    assert_receipt_field_is_bound_and_canonical("domain_note", &revisions, base.clone(), changed);

    let mut changed = base.clone();
    changed.policy_version = 2;
    assert_receipt_field_is_bound_and_canonical("policy_version", &revisions, base, changed);
}

#[test]
fn exact_duplicate_receipts_refuse_instead_of_changing_cardinality() {
    let a = revision("claim/a", "a holds");
    let b = revision("claim/b", "b holds");
    let receipt = edge(
        RelationKind::Implication,
        &a,
        &b,
        QuantifierVariance::Preserved,
    );
    let err = admit_graph(&[a, b], &[receipt.clone(), receipt]).expect_err("duplicate refuses");
    assert_eq!(err.rule(), "v1-duplicate-relation");
    assert!(err.ranked_fixes()[0].contains("deduplicate"));
}

#[test]
fn all_receipt_text_fields_escape_every_c0_and_hostile_line_character() {
    let a = revision("claim/a", "a holds");
    let b = revision("claim/b", "b holds");
    let mut hostile = "prefix\"\\".to_owned();
    hostile.extend((0_u32..=31).map(|value| char::from_u32(value).expect("C0 scalar")));
    hostile.push('\u{007f}');
    hostile.push('\u{0085}');
    hostile.push('\u{009f}');
    hostile.push('\u{2028}');
    hostile.push('\u{2029}');
    hostile.push('é');
    let escaped = concat!(
        "prefix\\\"\\\\",
        "\\u0000\\u0001\\u0002\\u0003\\u0004\\u0005\\u0006\\u0007",
        "\\b\\t\\n\\u000b\\f\\r\\u000e\\u000f",
        "\\u0010\\u0011\\u0012\\u0013\\u0014\\u0015\\u0016\\u0017",
        "\\u0018\\u0019\\u001a\\u001b\\u001c\\u001d\\u001e\\u001f",
        "\\u007f\\u0085\\u009f\\u2028\\u2029é"
    );

    for field in ["checker", "tcb", "domain_note"] {
        let mut receipt = edge(
            RelationKind::Implication,
            &a,
            &b,
            QuantifierVariance::Preserved,
        );
        match field {
            "checker" => receipt.checker.clone_from(&hostile),
            "tcb" => receipt.tcb.clone_from(&hostile),
            "domain_note" => receipt.domain_note.clone_from(&hostile),
            _ => unreachable!(),
        }
        let graph = admit_graph(&[a.clone(), b.clone()], &[receipt]).expect("hostile text admits");
        let json = graph.render_json_lines();
        assert!(
            json.last()
                .expect("edge record")
                .contains(format!("\"{field}\":\"{escaped}\"").as_str()),
            "{field} must use the exact JSON escaping policy"
        );
        assert!(json.iter().all(|line| {
            !line
                .chars()
                .any(|ch| ch.is_control() || matches!(ch, '\u{2028}' | '\u{2029}'))
        }));

        let human = graph.render_human();
        assert_eq!(
            human.lines().count(),
            1 + graph.revisions().len() + graph.edges().len()
        );
        assert!(human.contains(escaped));
        assert!(human.lines().all(|line| {
            !line
                .chars()
                .any(|ch| ch.is_control() || matches!(ch, '\u{2028}' | '\u{2029}'))
        }));

        assert!(graph.render_ledger_rows().iter().all(|(_, payload)| {
            payload.lines().count() == 1
                && !payload
                    .chars()
                    .any(|ch| ch.is_control() || matches!(ch, '\u{2028}' | '\u{2029}'))
        }));
    }
}

#[test]
fn domain_note_is_admitted_under_the_same_text_bounds_as_checker_and_tcb() {
    let a = revision("claim/a", "a holds");
    let b = revision("claim/b", "b holds");
    let mut receipt = edge(
        RelationKind::Implication,
        &a,
        &b,
        QuantifierVariance::Preserved,
    );
    receipt.domain_note.clear();
    let err = admit_graph(&[a.clone(), b.clone()], &[receipt]).expect_err("empty domain note");
    assert_eq!(err.rule(), "v1-field-bounds");
    assert!(err.ranked_fixes()[0].contains("domain_note"));

    let mut receipt = edge(
        RelationKind::Implication,
        &a,
        &b,
        QuantifierVariance::Preserved,
    );
    receipt.domain_note = "x".repeat(4096);
    assert!(
        admit_graph(&[a.clone(), b.clone()], &[receipt]).is_ok(),
        "the documented inclusive maximum admits"
    );

    let mut receipt = edge(
        RelationKind::Implication,
        &a,
        &b,
        QuantifierVariance::Preserved,
    );
    receipt.domain_note = "x".repeat(4097);
    let err = admit_graph(&[a, b], &[receipt]).expect_err("oversized domain note");
    assert_eq!(err.rule(), "v1-field-bounds");
}

#[test]
fn malformed_identities_and_empty_fields_refuse_with_stable_rules() {
    assert_eq!(ClaimId::new("").expect_err("empty").rule(), "v1-id-bounds");
    assert_eq!(
        ClaimId::new("Upper/Case")
            .expect_err("case folds identity")
            .rule(),
        "v1-id-bounds"
    );
    let mut hollow = revision("claim/a", "a holds");
    hollow.no_claim = String::new();
    let err = admit_graph(&[hollow], &[]).expect_err("empty no-claim severs");
    assert_eq!(err.rule(), "v1-field-bounds");
    assert!(err.ranked_fixes()[0].contains("no_claim"));
}
