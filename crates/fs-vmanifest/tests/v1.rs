//! VerificationManifest v1 identity-core conformance (bead i94v.7.1.1).
//!
//! The battery locks the bead's test plan: canonicalization and digest
//! stability under input reordering, collision-freedom and Unicode
//! distinctness, supersession immutability, source-authority conflict,
//! post-receipt source mutation, relation orientation and variance,
//! favorable relabeling, cycles, certified-equivalence SCCs, migration
//! classification, the per-field registry, and semantic parity of the
//! human/JSON/ledger projections.

use fs_vmanifest::v1::{
    ClaimId, ClaimKind, ClaimRelationReceipt, ClaimRevision, MANIFEST_RECORD_FIELDS,
    QuantifierVariance, RelationKind, SourceAuthority, SourcePin, admit_graph, classify_migration,
    resolve_authority,
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

fn edge(
    kind: RelationKind,
    from: &ClaimRevision,
    to: &ClaimRevision,
    variance: QuantifierVariance,
) -> ClaimRelationReceipt {
    ClaimRelationReceipt {
        kind,
        from: from.revision_id(),
        to: to.revision_id(),
        checker: "fs-checker/relation-v1".to_owned(),
        tcb: "rustc+fs-blake3".to_owned(),
        variance,
        domain_note: "identical domain".to_owned(),
        policy_version: 1,
    }
}

#[test]
fn identical_content_is_idempotent_and_distinct_content_cannot_collide() {
    let a = revision("claim/energy-balance", "energy balances to 1e-12");
    let b = revision("claim/energy-balance", "energy balances to 1e-12");
    assert_eq!(
        a.revision_id(),
        b.revision_id(),
        "content addressing is idempotent"
    );

    let c = revision("claim/energy-balance", "energy balances to 1e-11");
    assert_ne!(
        a.revision_id(),
        c.revision_id(),
        "one character = new revision"
    );

    // Unicode distinctness: visually confusable statements are distinct
    // content, never folded.
    let latin = revision("claim/units", "tolerance is 5 µm");
    let confusable = revision("claim/units", "tolerance is 5 μm");
    assert_ne!(latin.revision_id(), confusable.revision_id());
}

#[test]
fn supersession_appends_and_never_mutates_the_old_revision() {
    let old = revision("claim/energy-balance", "energy balances to 1e-11");
    let old_id = old.revision_id();
    let mut new = revision("claim/energy-balance", "energy balances to 1e-12");
    new.supersedes = Some(old_id);
    let new_id = new.revision_id();
    assert_ne!(old_id, new_id);
    // The old revision's identity is untouched by being superseded.
    assert_eq!(old.revision_id(), old_id);
    // Supersession participates in identity: the same statement WITHOUT
    // the supersession pointer is a different revision.
    let mut orphan = new.clone();
    orphan.supersedes = None;
    assert_ne!(orphan.revision_id(), new_id);
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
    let ida = a.revision_id();
    let idb = b.revision_id();
    let representative = ida.min(idb);
    assert_eq!(graph.representative[&ida], representative);
    assert_eq!(graph.representative[&idb], representative);
    assert_eq!(
        graph.revisions.len(),
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
    assert!(human.contains(&digest));
    assert!(json[0].contains(&digest));
    assert_eq!(ledger[0].1, digest);
    // Same edge count in every projection (semantic parity, not
    // byte-identity).
    assert_eq!(json.len() - 1, 1);
    assert_eq!(ledger.len() - 1, 1);
    assert_eq!(human.matches("->").count(), 1);
    assert_ne!(
        human,
        json.join("\n"),
        "projections need not be byte-identical"
    );
}

#[test]
fn json_lines_escape_untrusted_checker_text() {
    let a = revision("claim/a", "a holds");
    let b = revision("claim/b", "b holds");
    let mut receipt = edge(
        RelationKind::Implication,
        &a,
        &b,
        QuantifierVariance::Preserved,
    );
    receipt.checker = "checker\"\\\u{0008}\u{000c}\n\r\t\u{0000}\u{001f}é".to_owned();
    let graph = admit_graph(&[a.clone(), b.clone()], &[receipt]).expect("graph admits");

    let json = graph.render_json_lines();
    assert_eq!(json.len(), 2, "one edge remains one JSON-lines record");
    assert_eq!(
        json[1],
        format!(
            "{{\"vmanifest\":\"edge\",\"kind\":\"Implication\",\"from\":\"{}\",\"to\":\"{}\",\"checker\":\"checker\\\"\\\\\\b\\f\\n\\r\\t\\u0000\\u001fé\",\"policy\":1}}",
            a.revision_id().to_hex(),
            b.revision_id().to_hex()
        )
    );
    assert_eq!(json[1].lines().count(), 1);
    assert!(
        json[1].chars().all(|ch| ch >= ' '),
        "JSON-lines output must not contain literal control characters"
    );
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
