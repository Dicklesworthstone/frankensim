//! Machine-IR graph-bound datum-system admission (Gauntlet G0/G3/G5).

use core::num::NonZeroU64;

use fs_ir::machine::manufacturing::datum_system::{
    DatumFeatureBindingV1, DatumFeatureIdV1, DatumFeatureTargetV1, DatumPrecedenceV1,
    DatumReferenceFrameIdV1, DatumReferenceFrameV1, MAX_MACHINE_DATUM_FEATURES_V1,
    MAX_MACHINE_DATUM_REFERENCE_FRAMES_V1, MachineDatumAdmissionErrorV1, MachineDatumSystemDraftV1,
};
use fs_ir::machine::{
    AdmittedMachineGraph, BodyId, ContactFeatureId, MachineGraphDraft, MaterialBinding,
    MaterialCardRef, MaterialTarget, ModelRef, SubsystemId, SubsystemSpec, SurfacePatchId,
};

fn nz(value: u64) -> NonZeroU64 {
    NonZeroU64::new(value).expect("fixture version is nonzero")
}

fn body(key: &str) -> BodyId {
    BodyId::new(key).expect("fixture body key is canonical")
}

fn patch(key: &str) -> SurfacePatchId {
    SurfacePatchId::new(key).expect("fixture surface key is canonical")
}

fn contact(key: &str) -> ContactFeatureId {
    ContactFeatureId::new(key).expect("fixture contact key is canonical")
}

fn datum(key: &str) -> DatumFeatureIdV1 {
    DatumFeatureIdV1::new(key).expect("fixture datum key is canonical")
}

fn frame_id(key: &str) -> DatumReferenceFrameIdV1 {
    DatumReferenceFrameIdV1::new(key).expect("fixture frame key is canonical")
}

fn material(body: BodyId, key: &str, byte: u8) -> MaterialBinding {
    MaterialBinding {
        target: MaterialTarget::Body(body),
        material: MaterialCardRef::new(key, nz(1), [byte; 32])
            .expect("fixture material is canonical"),
    }
}

fn admitted_graph(model_byte: u8) -> AdmittedMachineGraph {
    let part = body("body/part");
    let alternate = body("body/part-alternate");
    let other = body("body/other");
    MachineGraphDraft {
        clocks: Vec::new(),
        subsystems: vec![
            SubsystemSpec {
                id: SubsystemId::new("subsystem/part").expect("canonical subsystem"),
                model: ModelRef::new("models/datum-part", nz(1), [model_byte; 32])
                    .expect("canonical model"),
                bodies: vec![part.clone(), alternate.clone()],
                surface_patches: vec![
                    patch("surface/part/primary"),
                    patch("surface/part/secondary"),
                ],
                contact_features: vec![contact("contact/part/tertiary")],
                state_slots: Vec::new(),
            },
            SubsystemSpec {
                id: SubsystemId::new("subsystem/other").expect("canonical subsystem"),
                model: ModelRef::new("models/datum-other", nz(1), [0x52; 32])
                    .expect("canonical model"),
                bodies: vec![other.clone()],
                surface_patches: vec![patch("surface/other/primary")],
                contact_features: Vec::new(),
                state_slots: Vec::new(),
            },
        ],
        terminals: Vec::new(),
        ports: Vec::new(),
        relations: Vec::new(),
        materials: vec![
            material(part, "materials/part", 1),
            material(alternate, "materials/part-alternate", 2),
            material(other, "materials/other", 3),
        ],
        interfaces: Vec::new(),
    }
    .admit()
    .expect("datum fixture graph admits")
}

fn surface_binding(id: &str, declared_body: &str, surface: &str) -> DatumFeatureBindingV1 {
    DatumFeatureBindingV1::new(
        datum(id),
        body(declared_body),
        DatumFeatureTargetV1::SurfacePatch(patch(surface)),
    )
}

fn contact_binding(id: &str, declared_body: &str, feature: &str) -> DatumFeatureBindingV1 {
    DatumFeatureBindingV1::new(
        datum(id),
        body(declared_body),
        DatumFeatureTargetV1::ContactFeature(contact(feature)),
    )
}

fn reference_frame(
    id: &str,
    primary: &str,
    secondary: Option<&str>,
    tertiary: Option<&str>,
) -> DatumReferenceFrameV1 {
    DatumReferenceFrameV1::new(
        frame_id(id),
        datum(primary),
        secondary.map(datum),
        tertiary.map(datum),
    )
}

fn valid_draft() -> MachineDatumSystemDraftV1 {
    MachineDatumSystemDraftV1 {
        datum_features: vec![
            surface_binding("datum/a", "body/part", "surface/part/primary"),
            surface_binding("datum/b", "body/part", "surface/part/secondary"),
            contact_binding("datum/c", "body/part", "contact/part/tertiary"),
            surface_binding("datum/d", "body/other", "surface/other/primary"),
        ],
        reference_frames: vec![
            reference_frame(
                "datum-frame/part",
                "datum/a",
                Some("datum/b"),
                Some("datum/c"),
            ),
            reference_frame("datum-frame/other", "datum/d", None, None),
        ],
    }
}

fn single_draft(binding: DatumFeatureBindingV1, frame: &str) -> MachineDatumSystemDraftV1 {
    let feature = binding.id().clone();
    MachineDatumSystemDraftV1 {
        datum_features: vec![binding],
        reference_frames: vec![DatumReferenceFrameV1::new(
            frame_id(frame),
            feature,
            None,
            None,
        )],
    }
}

#[test]
fn md_001_datum_precedence_is_graph_bound_order_invariant_and_identity_complete() {
    let graph = admitted_graph(0x41);
    let baseline = valid_draft()
        .admit_against(&graph)
        .expect("valid datum catalog admits");

    let mut reordered_draft = valid_draft();
    reordered_draft.datum_features.reverse();
    reordered_draft.reference_frames.reverse();
    let reordered = reordered_draft
        .admit_against(&graph)
        .expect("caller collection order is non-semantic");

    assert_eq!(baseline.graph(), graph.identity());
    assert_eq!(baseline.identity(), reordered.identity());
    assert_eq!(
        baseline.identity_receipt().canonical_preimage(),
        reordered.identity_receipt().canonical_preimage()
    );
    assert_eq!(
        baseline
            .datum_features()
            .iter()
            .map(|binding| binding.id().canonical_key())
            .collect::<Vec<_>>(),
        ["datum/a", "datum/b", "datum/c", "datum/d"]
    );
    assert_eq!(
        baseline
            .reference_frames()
            .iter()
            .map(|frame| frame.id().canonical_key())
            .collect::<Vec<_>>(),
        ["datum-frame/other", "datum-frame/part"]
    );
    let part_frame = &baseline.reference_frames()[1];
    assert_eq!(part_frame.reference_count(), 3);
    assert_eq!(part_frame.primary().canonical_key(), "datum/a");
    assert_eq!(
        part_frame.secondary().expect("secondary").canonical_key(),
        "datum/b"
    );
    assert_eq!(
        part_frame.tertiary().expect("tertiary").canonical_key(),
        "datum/c"
    );
    assert_eq!(DatumPrecedenceV1::Primary.tag(), 1);
    assert_eq!(DatumPrecedenceV1::Secondary.name(), "secondary");

    let base_single = single_draft(
        surface_binding("datum/single", "body/part", "surface/part/primary"),
        "datum-frame/single",
    )
    .admit_against(&graph)
    .expect("single datum admits");
    for changed in [
        single_draft(
            surface_binding("datum/single-renamed", "body/part", "surface/part/primary"),
            "datum-frame/single",
        ),
        single_draft(
            surface_binding(
                "datum/single",
                "body/part-alternate",
                "surface/part/primary",
            ),
            "datum-frame/single",
        ),
        single_draft(
            contact_binding("datum/single", "body/part", "contact/part/tertiary"),
            "datum-frame/single",
        ),
        single_draft(
            surface_binding("datum/single", "body/part", "surface/part/primary"),
            "datum-frame/single-renamed",
        ),
    ] {
        let changed = changed
            .admit_against(&graph)
            .expect("identity mutation remains structurally admissible");
        assert_ne!(base_single.identity(), changed.identity());
    }
    let other_graph = admitted_graph(0x42);
    let other_graph_receipt = single_draft(
        surface_binding("datum/single", "body/part", "surface/part/primary"),
        "datum-frame/single",
    )
    .admit_against(&other_graph)
    .expect("same selectors admit against changed graph");
    assert_ne!(base_single.identity(), other_graph_receipt.identity());

    let swapped = MachineDatumSystemDraftV1 {
        datum_features: vec![
            surface_binding("datum/a", "body/part", "surface/part/primary"),
            surface_binding("datum/b", "body/part", "surface/part/secondary"),
        ],
        reference_frames: vec![reference_frame(
            "datum-frame/part",
            "datum/b",
            Some("datum/a"),
            None,
        )],
    }
    .admit_against(&graph)
    .expect("swapped precedence is explicit and structurally admissible");
    let unswapped = MachineDatumSystemDraftV1 {
        datum_features: vec![
            surface_binding("datum/a", "body/part", "surface/part/primary"),
            surface_binding("datum/b", "body/part", "surface/part/secondary"),
        ],
        reference_frames: vec![reference_frame(
            "datum-frame/part",
            "datum/a",
            Some("datum/b"),
            None,
        )],
    }
    .admit_against(&graph)
    .expect("unswapped precedence admits");
    assert_ne!(swapped.identity(), unswapped.identity());
}

#[test]
#[allow(clippy::too_many_lines)]
fn md_002_datum_admission_refuses_aliases_gaps_and_unprovable_ownership() {
    let graph = admitted_graph(0x41);

    assert_eq!(
        MachineDatumSystemDraftV1 {
            datum_features: Vec::new(),
            reference_frames: vec![reference_frame("datum-frame/part", "datum/a", None, None,)],
        }
        .admit_against(&graph),
        Err(MachineDatumAdmissionErrorV1::NoDatumFeatures)
    );
    assert_eq!(
        MachineDatumSystemDraftV1 {
            datum_features: vec![surface_binding(
                "datum/a",
                "body/part",
                "surface/part/primary",
            )],
            reference_frames: Vec::new(),
        }
        .admit_against(&graph),
        Err(MachineDatumAdmissionErrorV1::NoReferenceFrames)
    );

    let duplicate_id = MachineDatumSystemDraftV1 {
        datum_features: vec![
            surface_binding("datum/a", "body/part", "surface/part/primary"),
            surface_binding("datum/a", "body/part", "surface/part/secondary"),
        ],
        reference_frames: vec![reference_frame("datum-frame/part", "datum/a", None, None)],
    };
    assert_eq!(
        duplicate_id.admit_against(&graph),
        Err(MachineDatumAdmissionErrorV1::DuplicateDatumFeature {
            feature: datum("datum/a"),
        })
    );

    let duplicate_selector = MachineDatumSystemDraftV1 {
        datum_features: vec![
            surface_binding("datum/a", "body/part", "surface/part/primary"),
            surface_binding("datum/b", "body/part", "surface/part/primary"),
        ],
        reference_frames: vec![reference_frame(
            "datum-frame/part",
            "datum/a",
            Some("datum/b"),
            None,
        )],
    };
    assert_eq!(
        duplicate_selector.admit_against(&graph),
        Err(MachineDatumAdmissionErrorV1::DuplicateDatumSelector {
            first: datum("datum/a"),
            duplicate: datum("datum/b"),
        })
    );
    let duplicate_target_with_changed_body = MachineDatumSystemDraftV1 {
        datum_features: vec![
            surface_binding("datum/a", "body/part", "surface/part/primary"),
            surface_binding("datum/b", "body/part-alternate", "surface/part/primary"),
        ],
        reference_frames: vec![reference_frame(
            "datum-frame/part",
            "datum/a",
            Some("datum/b"),
            None,
        )],
    };
    assert_eq!(
        duplicate_target_with_changed_body.admit_against(&graph),
        Err(MachineDatumAdmissionErrorV1::DuplicateDatumSelector {
            first: datum("datum/a"),
            duplicate: datum("datum/b"),
        })
    );

    let unknown_body = single_draft(
        surface_binding("datum/a", "body/missing", "surface/part/primary"),
        "datum-frame/part",
    );
    assert_eq!(
        unknown_body.admit_against(&graph),
        Err(MachineDatumAdmissionErrorV1::UnknownBody {
            feature: datum("datum/a"),
            body: body("body/missing"),
        })
    );

    let unknown_target = single_draft(
        surface_binding("datum/a", "body/part", "surface/part/missing"),
        "datum-frame/part",
    );
    assert_eq!(
        unknown_target.admit_against(&graph),
        Err(MachineDatumAdmissionErrorV1::UnknownFeatureTarget {
            feature: datum("datum/a"),
            target: DatumFeatureTargetV1::SurfacePatch(patch("surface/part/missing")),
        })
    );

    let owner_mismatch = single_draft(
        surface_binding("datum/a", "body/part", "surface/other/primary"),
        "datum-frame/part",
    );
    assert!(matches!(
        owner_mismatch.admit_against(&graph),
        Err(MachineDatumAdmissionErrorV1::FeatureOwnerMismatch {
            feature,
            body,
            target: DatumFeatureTargetV1::SurfacePatch(target),
            body_owner,
            target_owner,
        }) if feature == datum("datum/a")
            && body == crate::body("body/part")
            && target == patch("surface/other/primary")
            && body_owner == SubsystemId::new("subsystem/part").unwrap()
            && target_owner == SubsystemId::new("subsystem/other").unwrap()
    ));

    let duplicate_frame = MachineDatumSystemDraftV1 {
        datum_features: vec![surface_binding(
            "datum/a",
            "body/part",
            "surface/part/primary",
        )],
        reference_frames: vec![
            reference_frame("datum-frame/part", "datum/a", None, None),
            reference_frame("datum-frame/part", "datum/a", None, None),
        ],
    };
    assert_eq!(
        duplicate_frame.admit_against(&graph),
        Err(MachineDatumAdmissionErrorV1::DuplicateReferenceFrame {
            frame: frame_id("datum-frame/part"),
        })
    );

    let tertiary_gap = MachineDatumSystemDraftV1 {
        datum_features: vec![
            surface_binding("datum/a", "body/part", "surface/part/primary"),
            contact_binding("datum/c", "body/part", "contact/part/tertiary"),
        ],
        reference_frames: vec![reference_frame(
            "datum-frame/part",
            "datum/a",
            None,
            Some("datum/c"),
        )],
    };
    assert_eq!(
        tertiary_gap.admit_against(&graph),
        Err(MachineDatumAdmissionErrorV1::TertiaryWithoutSecondary {
            frame: frame_id("datum-frame/part"),
        })
    );

    let missing_reference = MachineDatumSystemDraftV1 {
        datum_features: vec![surface_binding(
            "datum/a",
            "body/part",
            "surface/part/primary",
        )],
        reference_frames: vec![reference_frame(
            "datum-frame/part",
            "datum/missing",
            None,
            None,
        )],
    };
    assert_eq!(
        missing_reference.admit_against(&graph),
        Err(MachineDatumAdmissionErrorV1::MissingDatumReference {
            frame: frame_id("datum-frame/part"),
            precedence: DatumPrecedenceV1::Primary,
            feature: datum("datum/missing"),
        })
    );

    let repeated_reference = MachineDatumSystemDraftV1 {
        datum_features: vec![surface_binding(
            "datum/a",
            "body/part",
            "surface/part/primary",
        )],
        reference_frames: vec![reference_frame(
            "datum-frame/part",
            "datum/a",
            Some("datum/a"),
            None,
        )],
    };
    assert_eq!(
        repeated_reference.admit_against(&graph),
        Err(MachineDatumAdmissionErrorV1::RepeatedDatumReference {
            frame: frame_id("datum-frame/part"),
            feature: datum("datum/a"),
            first: DatumPrecedenceV1::Primary,
            repeated: DatumPrecedenceV1::Secondary,
        })
    );

    let mixed_body = MachineDatumSystemDraftV1 {
        datum_features: vec![
            surface_binding("datum/a", "body/part", "surface/part/primary"),
            surface_binding("datum/d", "body/other", "surface/other/primary"),
        ],
        reference_frames: vec![reference_frame(
            "datum-frame/mixed",
            "datum/a",
            Some("datum/d"),
            None,
        )],
    };
    assert_eq!(
        mixed_body.admit_against(&graph),
        Err(MachineDatumAdmissionErrorV1::MixedBodyReferenceFrame {
            frame: frame_id("datum-frame/mixed"),
            first_body: body("body/part"),
            conflicting_body: body("body/other"),
        })
    );

    let unused = MachineDatumSystemDraftV1 {
        datum_features: vec![
            surface_binding("datum/a", "body/part", "surface/part/primary"),
            surface_binding("datum/b", "body/part", "surface/part/secondary"),
        ],
        reference_frames: vec![reference_frame("datum-frame/part", "datum/a", None, None)],
    };
    assert_eq!(
        unused.admit_against(&graph),
        Err(MachineDatumAdmissionErrorV1::UnusedDatumFeature {
            feature: datum("datum/b"),
        })
    );

    assert_eq!(
        MachineDatumAdmissionErrorV1::NoDatumFeatures.code(),
        "MachineDatumNoFeatures"
    );
}

fn boundary_graph() -> AdmittedMachineGraph {
    let part = body("body/boundary");
    let surfaces = (0..MAX_MACHINE_DATUM_FEATURES_V1)
        .map(|index| patch(&format!("surface/boundary/p{index:04}")))
        .collect::<Vec<_>>();
    MachineGraphDraft {
        clocks: Vec::new(),
        subsystems: vec![SubsystemSpec {
            id: SubsystemId::new("subsystem/boundary").expect("canonical subsystem"),
            model: ModelRef::new("models/datum-boundary", nz(1), [0x61; 32])
                .expect("canonical model"),
            bodies: vec![part.clone()],
            surface_patches: surfaces,
            contact_features: Vec::new(),
            state_slots: Vec::new(),
        }],
        terminals: Vec::new(),
        ports: Vec::new(),
        relations: Vec::new(),
        materials: vec![material(part, "materials/boundary", 0x62)],
        interfaces: Vec::new(),
    }
    .admit()
    .expect("exact-cap datum graph admits")
}

fn boundary_draft() -> MachineDatumSystemDraftV1 {
    let datum_features = (0..MAX_MACHINE_DATUM_FEATURES_V1)
        .map(|index| {
            surface_binding(
                &format!("datum/f{index:04}"),
                "body/boundary",
                &format!("surface/boundary/p{index:04}"),
            )
        })
        .collect::<Vec<_>>();
    let reference_frames = (0..MAX_MACHINE_DATUM_REFERENCE_FRAMES_V1)
        .map(|index| {
            reference_frame(
                &format!("datum-frame/f{index:04}"),
                &format!("datum/f{:04}", 2 * index),
                Some(&format!("datum/f{:04}", 2 * index + 1)),
                None,
            )
        })
        .collect();
    MachineDatumSystemDraftV1 {
        datum_features,
        reference_frames,
    }
}

#[test]
fn md_003_exact_resource_caps_admit_and_one_over_refuses_before_deduplication() {
    let graph = boundary_graph();
    let exact = boundary_draft();
    let admitted = exact
        .clone()
        .admit_against(&graph)
        .expect("simultaneous exact feature/frame caps admit");
    assert_eq!(
        admitted.datum_features().len(),
        MAX_MACHINE_DATUM_FEATURES_V1
    );
    assert_eq!(
        admitted.reference_frames().len(),
        MAX_MACHINE_DATUM_REFERENCE_FRAMES_V1
    );

    let mut too_many_features = exact.clone();
    let repeated_feature = too_many_features.datum_features[0].clone();
    too_many_features.datum_features.push(repeated_feature);
    assert_eq!(
        too_many_features.admit_against(&graph),
        Err(MachineDatumAdmissionErrorV1::DatumFeatureLimit {
            actual: MAX_MACHINE_DATUM_FEATURES_V1 + 1,
            max: MAX_MACHINE_DATUM_FEATURES_V1,
        })
    );

    let mut too_many_frames = exact;
    let repeated_frame = too_many_frames.reference_frames[0].clone();
    too_many_frames.reference_frames.push(repeated_frame);
    assert_eq!(
        too_many_frames.admit_against(&graph),
        Err(MachineDatumAdmissionErrorV1::ReferenceFrameLimit {
            actual: MAX_MACHINE_DATUM_REFERENCE_FRAMES_V1 + 1,
            max: MAX_MACHINE_DATUM_REFERENCE_FRAMES_V1,
        })
    );
}

#[test]
fn md_004_identical_input_replays_the_complete_receipt() {
    let graph = admitted_graph(0x41);
    let first = valid_draft()
        .admit_against(&graph)
        .expect("first replay admits");
    let second = valid_draft()
        .admit_against(&graph)
        .expect("second replay admits");
    assert_eq!(first, second);
    assert_eq!(
        first.identity_receipt().canonical_preimage(),
        second.identity_receipt().canonical_preimage()
    );
    println!(
        "{{\"suite\":\"fs-ir/machine-manufacturing-datum\",\"case\":\"md-004\",\
         \"features\":{},\"frames\":{},\"replay\":true}}",
        first.datum_features().len(),
        first.reference_frames().len()
    );
}
