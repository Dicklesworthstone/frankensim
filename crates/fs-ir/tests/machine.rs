//! Machine-IR PR-1 identity and lineage conformance battery.

use fs_blake3::identity::StrongIdentity;
use fs_ir::machine::{
    BodyId, ContactFeatureId, DependentBinding, DependentKind, LineageEvent, LineageRecord,
    LineageRefusal, LineageRelation, MAX_LINEAGE_DEPENDENTS, MAX_LINEAGE_ENDPOINTS,
    MAX_LINEAGE_RELATIONS, MAX_LINEAGE_TARGETS_PER_SOURCE, MachineElementId, MachineIdError,
    PortId, StateSlotId, SurfacePatchId, TerminalId,
};

#[test]
fn g0_entity_ids_are_repeatable_nominally_separated_and_human_auditable() {
    let body = BodyId::new("rotor/main-shaft").expect("canonical body id");
    let replay = BodyId::new("rotor/main-shaft").expect("canonical body replay");
    let patch = SurfacePatchId::new("rotor/main-shaft").expect("same key in another role");
    let contact = ContactFeatureId::new("rotor/main-shaft").expect("contact id");
    let terminal = TerminalId::new("rotor/main-shaft").expect("terminal id");
    let port = PortId::new("rotor/main-shaft").expect("port id");
    let state = StateSlotId::new("rotor/main-shaft").expect("state-slot id");

    assert_eq!(body, replay);
    assert_eq!(body.identity(), replay.identity());
    assert_eq!(body.canonical_key(), "rotor/main-shaft");
    assert_eq!(body.identity_receipt().field_count(), 1);
    assert_ne!(
        body.identity().as_bytes(),
        patch.identity().as_bytes(),
        "nominal roles must domain-separate identical text"
    );
    let role_digests = [
        *body.identity().as_bytes(),
        *patch.identity().as_bytes(),
        *contact.identity().as_bytes(),
        *terminal.identity().as_bytes(),
        *port.identity().as_bytes(),
        *state.identity().as_bytes(),
    ];
    for (left_index, left) in role_digests.iter().enumerate() {
        for right in role_digests.iter().skip(left_index + 1) {
            assert_ne!(
                left, right,
                "every nominal role requires a distinct static schema"
            );
        }
    }
    assert_ne!(
        MachineElementId::from(body),
        MachineElementId::from(patch),
        "erasing to the closed element enum must retain the nominal role"
    );
}

#[test]
fn g0_entity_key_grammar_refuses_alias_prone_or_index_only_spellings() {
    let empty = BodyId::new("").expect_err("empty key must refuse");
    assert_eq!(empty, MachineIdError::Empty { role: "body-id" });
    assert_eq!(empty.code(), "MachineIdEmpty");
    assert!(matches!(
        BodyId::new("Rotor/main"),
        Err(MachineIdError::InvalidSegmentStart {
            role: "body-id",
            segment: 0,
            at: 0,
            byte: b'R',
        })
    ));
    assert!(matches!(
        BodyId::new("rotor//main"),
        Err(MachineIdError::EmptySegment {
            role: "body-id",
            segment: 1,
        })
    ));
    assert!(matches!(
        BodyId::new("rotor/7"),
        Err(MachineIdError::InvalidSegmentStart {
            role: "body-id",
            segment: 1,
            byte: b'7',
            ..
        })
    ));
    assert_eq!(
        BodyId::new("rotor/main--shaft"),
        Err(MachineIdError::RepeatedSeparator {
            role: "body-id",
            at: 11,
        })
    );
    assert_eq!(
        BodyId::new("rotor/main-"),
        Err(MachineIdError::TrailingSeparator {
            role: "body-id",
            segment: 1,
        })
    );
    assert_eq!(
        BodyId::new("rotor/main_shaft"),
        Err(MachineIdError::InvalidByte {
            role: "body-id",
            at: 10,
            byte: b'_',
        })
    );
    assert!(matches!(
        BodyId::new(format!("a{}", "0".repeat(128))),
        Err(MachineIdError::TooLong { bytes: 129, .. })
    ));
    assert_eq!(
        BodyId::new(format!("a{}", "0".repeat(127)))
            .expect("the exact 128-byte boundary is admitted")
            .canonical_key()
            .len(),
        128
    );
}

#[test]
fn g3_unique_lineage_rebindings_are_order_canonical() {
    let old_body = MachineElementId::from(BodyId::new("rotor/body-v1").unwrap());
    let new_body = MachineElementId::from(BodyId::new("rotor/body-v2").unwrap());
    let old_patch = MachineElementId::from(SurfacePatchId::new("rotor/skin-v1").unwrap());
    let new_patch = MachineElementId::from(SurfacePatchId::new("rotor/skin-v2").unwrap());

    let body_relation = LineageRelation::new(old_body.clone(), vec![new_body.clone()]).unwrap();
    let patch_relation = LineageRelation::new(old_patch.clone(), vec![new_patch.clone()]).unwrap();
    let cache = DependentBinding::new(DependentKind::Cache, "stress/cache-1", old_body).unwrap();
    let winding =
        DependentBinding::new(DependentKind::Winding, "motor/winding-a", old_patch).unwrap();

    let first_decision = LineageRecord::admit_with_decision(
        LineageEvent::Remesh,
        vec![patch_relation.clone(), body_relation.clone()],
        vec![winding.clone(), cache.clone()],
    );
    assert_eq!(first_decision.code(), "LineageAdmitted");
    assert_eq!(first_decision.submitted_relation_count(), 2);
    assert_eq!(first_decision.submitted_dependent_count(), 2);
    let first = first_decision
        .into_result()
        .expect("one-to-one remesh is unambiguous");
    let replay = LineageRecord::admit(
        LineageEvent::Remesh,
        vec![body_relation, patch_relation],
        vec![cache, winding],
    )
    .expect("caller order is not semantic");

    assert_eq!(first, replay);
    assert_eq!(first.identity(), replay.identity());
    assert_eq!(first.identity_receipt().field_count(), 3);
    assert_eq!(first.rebindings().len(), 2);
    assert!(first.rebindings().iter().any(|binding| {
        binding.dependent().kind() == DependentKind::Cache && binding.target() == &new_body
    }));
    assert!(first.rebindings().iter().any(|binding| {
        binding.dependent().kind() == DependentKind::Winding && binding.target() == &new_patch
    }));
}

#[test]
fn g0_split_without_attachments_emits_lineage_but_grants_no_implicit_crosswalk() {
    let source = MachineElementId::from(SurfacePatchId::new("housing/seam").unwrap());
    let left = MachineElementId::from(SurfacePatchId::new("housing/seam-left").unwrap());
    let right = MachineElementId::from(SurfacePatchId::new("housing/seam-right").unwrap());
    let relation = LineageRelation::new(source, vec![right, left]).unwrap();

    let record = LineageRecord::admit(LineageEvent::Split, vec![relation], Vec::new())
        .expect("a split can be recorded when there is nothing to rebind");
    assert!(record.rebindings().is_empty());
    assert_eq!(record.relations()[0].targets().len(), 2);
}

#[test]
fn g0_ambiguous_split_refuses_and_enumerates_every_invalidated_dependent() {
    let source = MachineElementId::from(SurfacePatchId::new("housing/seam").unwrap());
    let left = MachineElementId::from(SurfacePatchId::new("housing/seam-left").unwrap());
    let right = MachineElementId::from(SurfacePatchId::new("housing/seam-right").unwrap());
    let relation = LineageRelation::new(source.clone(), vec![right, left]).unwrap();
    let mut dependents = vec![
        DependentBinding::new(DependentKind::Cache, "solver/cache-a", source.clone()).unwrap(),
        DependentBinding::new(DependentKind::Contact, "contact/pair-a", source.clone()).unwrap(),
        DependentBinding::new(DependentKind::Winding, "motor/winding-a", source.clone()).unwrap(),
        DependentBinding::new(DependentKind::Adjoint, "adjoint/tape-a", source).unwrap(),
    ];

    let decision = LineageRecord::admit_with_decision(
        LineageEvent::Split,
        vec![relation.clone()],
        dependents.clone(),
    );
    assert_eq!(decision.event(), LineageEvent::Split);
    assert_eq!(decision.submitted_relation_count(), 1);
    assert_eq!(decision.submitted_dependent_count(), 4);
    assert_eq!(decision.code(), "LineageAmbiguous");
    let refusal = decision
        .into_result()
        .expect_err("one source with two targets cannot silently rebind attachments");
    let invalidation = refusal.invalidation().expect("typed ambiguity payload");
    assert_eq!(invalidation.relations(), &[relation.clone()]);
    assert_eq!(invalidation.considered_dependents().len(), 4);
    assert_eq!(invalidation.ambiguous_relations(), &[relation.clone()]);
    assert_eq!(invalidation.invalidated_dependents().len(), 4);
    assert_eq!(invalidation.identity_receipt().field_count(), 5);
    assert_eq!(
        invalidation
            .invalidated_dependents()
            .iter()
            .map(DependentBinding::kind)
            .collect::<Vec<_>>(),
        vec![
            DependentKind::Cache,
            DependentKind::Contact,
            DependentKind::Winding,
            DependentKind::Adjoint,
        ]
    );

    dependents.reverse();
    let replay = LineageRecord::admit(LineageEvent::Split, vec![relation], dependents)
        .expect_err("reordered ambiguity still refuses");
    assert_eq!(
        invalidation.identity(),
        replay.invalidation().unwrap().identity(),
        "the invalidation receipt must be caller-order invariant"
    );
}

#[test]
fn g3_ambiguous_invalidation_identity_binds_the_complete_refused_event() {
    let seam = MachineElementId::from(SurfacePatchId::new("housing/seam").unwrap());
    let seam_left = MachineElementId::from(SurfacePatchId::new("housing/seam-left").unwrap());
    let seam_right = MachineElementId::from(SurfacePatchId::new("housing/seam-right").unwrap());
    let ambiguous = LineageRelation::new(seam.clone(), vec![seam_left, seam_right]).unwrap();
    let attachment = DependentBinding::new(DependentKind::Contact, "contact/seam", seam).unwrap();

    let body = MachineElementId::from(BodyId::new("housing/body-v1").unwrap());
    let body_v2 = MachineElementId::from(BodyId::new("housing/body-v2").unwrap());
    let body_v3 = MachineElementId::from(BodyId::new("housing/body-v3").unwrap());
    let unchanged_context = LineageRelation::new(body.clone(), vec![body_v2]).unwrap();
    let changed_context = LineageRelation::new(body.clone(), vec![body_v3]).unwrap();
    let considered_a =
        DependentBinding::new(DependentKind::Cache, "solver/cache-a", body.clone()).unwrap();
    let considered_b = DependentBinding::new(DependentKind::Cache, "solver/cache-b", body).unwrap();

    let first = LineageRecord::admit(
        LineageEvent::Remesh,
        vec![ambiguous.clone(), unchanged_context.clone()],
        vec![attachment.clone(), considered_a.clone()],
    )
    .expect_err("ambiguous attachment must fail closed");
    let changed_relation = LineageRecord::admit(
        LineageEvent::Remesh,
        vec![ambiguous.clone(), changed_context],
        vec![attachment.clone(), considered_a],
    )
    .expect_err("same ambiguity in a different event must fail closed");
    let changed_dependent = LineageRecord::admit(
        LineageEvent::Remesh,
        vec![ambiguous, unchanged_context],
        vec![attachment, considered_b],
    )
    .expect_err("considered dependent changes remain identity-bearing");

    let first = first.invalidation().expect("typed invalidation");
    let changed_relation = changed_relation.invalidation().expect("typed invalidation");
    let changed_dependent = changed_dependent
        .invalidation()
        .expect("typed invalidation");
    assert_eq!(first.relations().len(), 2);
    assert_eq!(first.considered_dependents().len(), 2);
    assert_eq!(first.ambiguous_relations().len(), 1);
    assert_eq!(first.invalidated_dependents().len(), 1);
    assert_ne!(
        first.identity(),
        changed_relation.identity(),
        "non-ambiguous relations remain part of the refused event identity"
    );
    assert_ne!(
        first.identity(),
        changed_dependent.identity(),
        "non-invalidated considered dependents remain identity-bearing inputs"
    );
}

#[test]
fn g0_lineage_shape_and_attachment_errors_name_the_offending_identity() {
    let body = MachineElementId::from(BodyId::new("rotor/body").unwrap());
    let body_next = MachineElementId::from(BodyId::new("rotor/body-next").unwrap());
    let patch = MachineElementId::from(SurfacePatchId::new("rotor/skin").unwrap());
    let patch_other = MachineElementId::from(SurfacePatchId::new("rotor/skin-other").unwrap());

    let forward = LineageRelation::new(body.clone(), vec![patch.clone(), patch_other.clone()])
        .expect_err("cross-role targets must refuse");
    let reverse = LineageRelation::new(body.clone(), vec![patch_other, patch.clone()])
        .expect_err("caller order cannot select a different diagnostic target");
    assert_eq!(forward, reverse);
    assert_eq!(forward.code(), "LineageTargetKindMismatch");
    assert!(matches!(
        forward,
        LineageRefusal::TargetKindMismatch { source, target }
            if source == body && target.kind() == fs_ir::machine::MachineElementKind::SurfacePatch
    ));

    let one_target = LineageRelation::new(body.clone(), vec![body_next]).unwrap();
    assert!(matches!(
        LineageRecord::admit(LineageEvent::Split, vec![one_target.clone()], Vec::new()),
        Err(LineageRefusal::EventShape {
            event: LineageEvent::Split,
            ..
        })
    ));

    let missing = MachineElementId::from(BodyId::new("rotor/missing").unwrap());
    let dependent =
        DependentBinding::new(DependentKind::Cache, "solver/cache-a", missing.clone()).unwrap();
    assert!(matches!(
        LineageRecord::admit(LineageEvent::Remesh, vec![one_target], vec![dependent]),
        Err(LineageRefusal::UnknownDependentSource { source, .. }) if source == missing
    ));
}

#[test]
fn g0_event_cardinality_laws_cover_merge_fracture_and_wear() {
    let source_a = MachineElementId::from(BodyId::new("assembly/body-a").unwrap());
    let source_b = MachineElementId::from(BodyId::new("assembly/body-b").unwrap());
    let merged = MachineElementId::from(BodyId::new("assembly/body-merged").unwrap());
    let other = MachineElementId::from(BodyId::new("assembly/body-other").unwrap());
    let relation_a = LineageRelation::new(source_a, vec![merged.clone()]).unwrap();
    let relation_b = LineageRelation::new(source_b.clone(), vec![merged]).unwrap();

    let merge = LineageRecord::admit(
        LineageEvent::Merge,
        vec![relation_b.clone(), relation_a.clone()],
        Vec::new(),
    )
    .expect("many sources with one shared target form a merge");
    assert_eq!(merge.event(), LineageEvent::Merge);
    assert_eq!(merge.relations().len(), 2);

    let mismatched_target = LineageRelation::new(source_b, vec![other]).unwrap();
    assert!(matches!(
        LineageRecord::admit(
            LineageEvent::Merge,
            vec![relation_a, mismatched_target],
            Vec::new(),
        ),
        Err(LineageRefusal::EventShape {
            event: LineageEvent::Merge,
            ..
        })
    ));

    let source = MachineElementId::from(SurfacePatchId::new("housing/crack").unwrap());
    let left = MachineElementId::from(SurfacePatchId::new("housing/crack-left").unwrap());
    let right = MachineElementId::from(SurfacePatchId::new("housing/crack-right").unwrap());
    let fracture = LineageRelation::new(source, vec![left, right]).unwrap();
    assert_eq!(
        LineageRecord::admit(LineageEvent::Fracture, vec![fracture.clone()], Vec::new(),)
            .expect("fracture has one source and multiple descendants")
            .event(),
        LineageEvent::Fracture
    );
    assert!(matches!(
        LineageRecord::admit(LineageEvent::Wear, vec![fracture], Vec::new()),
        Err(LineageRefusal::EventShape {
            event: LineageEvent::Wear,
            ..
        })
    ));
}

#[test]
fn g0_duplicate_and_public_count_limits_refuse_structurally() {
    let source = MachineElementId::from(BodyId::new("limits/source").unwrap());
    let target = MachineElementId::from(BodyId::new("limits/target").unwrap());

    let no_relations =
        LineageRecord::admit_with_decision(LineageEvent::Remesh, Vec::new(), Vec::new());
    assert_eq!(no_relations.code(), "LineageNoRelations");
    assert!(matches!(
        no_relations.into_result(),
        Err(LineageRefusal::NoRelations)
    ));
    let no_targets = LineageRelation::new(source.clone(), Vec::new())
        .expect_err("every source requires at least one explicit target");
    assert_eq!(no_targets.code(), "LineageNoTargets");

    assert!(matches!(
        LineageRelation::new(source.clone(), vec![target.clone(), target.clone()]),
        Err(LineageRefusal::DuplicateTarget { .. })
    ));
    assert!(matches!(
        LineageRelation::new(
            source.clone(),
            vec![target.clone(); MAX_LINEAGE_TARGETS_PER_SOURCE + 1],
        ),
        Err(LineageRefusal::TargetLimit { count, max, .. })
            if count == MAX_LINEAGE_TARGETS_PER_SOURCE + 1
                && max == MAX_LINEAGE_TARGETS_PER_SOURCE
    ));

    let relation = LineageRelation::new(source.clone(), vec![target]).unwrap();
    assert!(matches!(
        LineageRecord::admit(
            LineageEvent::Remesh,
            vec![relation.clone(), relation.clone()],
            Vec::new(),
        ),
        Err(LineageRefusal::DuplicateSource { .. })
    ));
    assert!(matches!(
        LineageRecord::admit(
            LineageEvent::Remesh,
            vec![relation.clone(); MAX_LINEAGE_RELATIONS + 1],
            Vec::new(),
        ),
        Err(LineageRefusal::RelationLimit { count, max })
            if count == MAX_LINEAGE_RELATIONS + 1 && max == MAX_LINEAGE_RELATIONS
    ));

    let dependent = DependentBinding::new(DependentKind::Cache, "limits/cache", source).unwrap();
    assert!(matches!(
        LineageRecord::admit(
            LineageEvent::Remesh,
            vec![relation.clone()],
            vec![dependent.clone(), dependent.clone()],
        ),
        Err(LineageRefusal::DuplicateDependent { .. })
    ));
    assert!(matches!(
        LineageRecord::admit(
            LineageEvent::Remesh,
            vec![relation],
            vec![dependent; MAX_LINEAGE_DEPENDENTS + 1],
        ),
        Err(LineageRefusal::DependentLimit { count, max })
            if count == MAX_LINEAGE_DEPENDENTS + 1 && max == MAX_LINEAGE_DEPENDENTS
    ));
}

#[test]
fn g0_maximum_endpoint_payload_fits_identity_envelope_and_plus_one_refuses() {
    let source_a = MachineElementId::from(BodyId::new("limits/source-a").unwrap());
    let source_b = MachineElementId::from(BodyId::new("limits/source-b").unwrap());
    let targets_a = (0..MAX_LINEAGE_TARGETS_PER_SOURCE)
        .map(|index| {
            MachineElementId::from(
                BodyId::new(format!("limits/target-a-{index}")).expect("bounded target id"),
            )
        })
        .collect();
    let targets_b = (0..MAX_LINEAGE_TARGETS_PER_SOURCE)
        .map(|index| {
            MachineElementId::from(
                BodyId::new(format!("limits/target-b-{index}")).expect("bounded target id"),
            )
        })
        .collect();
    let relation_a = LineageRelation::new(source_a, targets_a).unwrap();
    let relation_b = LineageRelation::new(source_b, targets_b).unwrap();

    let maximum = LineageRecord::admit(
        LineageEvent::Remesh,
        vec![relation_a.clone(), relation_b.clone()],
        Vec::new(),
    )
    .expect("the public maximum endpoint payload must fit the identity envelope");
    assert_eq!(
        maximum
            .relations()
            .iter()
            .map(|relation| relation.targets().len())
            .sum::<usize>(),
        MAX_LINEAGE_ENDPOINTS
    );

    let extra_source = MachineElementId::from(BodyId::new("limits/source-extra").unwrap());
    let extra_target = MachineElementId::from(BodyId::new("limits/target-extra").unwrap());
    let extra = LineageRelation::new(extra_source, vec![extra_target]).unwrap();
    assert!(matches!(
        LineageRecord::admit(
            LineageEvent::Remesh,
            vec![relation_a, relation_b, extra],
            Vec::new(),
        ),
        Err(LineageRefusal::EndpointLimit { count, max })
            if count == MAX_LINEAGE_ENDPOINTS + 1 && max == MAX_LINEAGE_ENDPOINTS
    ));
}

#[test]
fn g0_maximum_relation_and_dependent_invalidation_fits_identity_envelope() {
    let mut relations = Vec::with_capacity(MAX_LINEAGE_RELATIONS);
    let mut dependents = Vec::with_capacity(MAX_LINEAGE_DEPENDENTS);
    for index in 0..MAX_LINEAGE_RELATIONS {
        let source = MachineElementId::from(
            SurfacePatchId::new(format!("limits/patch-source-{index}")).expect("bounded source id"),
        );
        let left = MachineElementId::from(
            SurfacePatchId::new(format!("limits/patch-left-{index}"))
                .expect("bounded left target id"),
        );
        let right = MachineElementId::from(
            SurfacePatchId::new(format!("limits/patch-right-{index}"))
                .expect("bounded right target id"),
        );
        dependents.push(
            DependentBinding::new(
                DependentKind::Cache,
                format!("limits/cache-{index}"),
                source.clone(),
            )
            .expect("bounded dependent id"),
        );
        relations.push(LineageRelation::new(source, vec![left, right]).unwrap());
    }

    let refusal = LineageRecord::admit(LineageEvent::Remesh, relations, dependents)
        .expect_err("live dependents on maximum one-to-many relations must fail closed");
    let invalidation = refusal.invalidation().expect("typed invalidation");
    assert_eq!(invalidation.relations().len(), MAX_LINEAGE_RELATIONS);
    assert_eq!(
        invalidation.considered_dependents().len(),
        MAX_LINEAGE_DEPENDENTS
    );
    assert_eq!(
        invalidation.ambiguous_relations().len(),
        MAX_LINEAGE_RELATIONS
    );
    assert_eq!(
        invalidation.invalidated_dependents().len(),
        MAX_LINEAGE_DEPENDENTS
    );
    assert_eq!(
        invalidation.identity_receipt().collection_items(),
        2 * (MAX_LINEAGE_RELATIONS + MAX_LINEAGE_DEPENDENTS) as u64
    );
}

#[test]
fn g3_lineage_identity_moves_with_event_target_and_dependent_class() {
    let source = MachineElementId::from(BodyId::new("rotor/body-v1").unwrap());
    let target_a = MachineElementId::from(BodyId::new("rotor/body-v2").unwrap());
    let target_b = MachineElementId::from(BodyId::new("rotor/body-v3").unwrap());
    let relation_a = LineageRelation::new(source.clone(), vec![target_a]).unwrap();
    let relation_b = LineageRelation::new(source.clone(), vec![target_b]).unwrap();
    let cache =
        DependentBinding::new(DependentKind::Cache, "solver/derived-a", source.clone()).unwrap();
    let adjoint =
        DependentBinding::new(DependentKind::Adjoint, "solver/derived-a", source).unwrap();

    let remesh = LineageRecord::admit(
        LineageEvent::Remesh,
        vec![relation_a.clone()],
        vec![cache.clone()],
    )
    .unwrap();
    let wear = LineageRecord::admit(
        LineageEvent::Wear,
        vec![relation_a.clone()],
        vec![cache.clone()],
    )
    .unwrap();
    let changed_target =
        LineageRecord::admit(LineageEvent::Remesh, vec![relation_b], vec![cache]).unwrap();
    let changed_kind =
        LineageRecord::admit(LineageEvent::Remesh, vec![relation_a], vec![adjoint]).unwrap();

    assert_ne!(remesh.identity(), wear.identity());
    assert_ne!(remesh.identity(), changed_target.identity());
    assert_ne!(remesh.identity(), changed_kind.identity());
}
