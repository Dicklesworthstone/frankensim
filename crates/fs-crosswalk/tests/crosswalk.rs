//! Battery for the regulatory crosswalk (addendum Proposal 12). Verifies the
//! table is complete (every concept × standard covered, no silent gaps), that
//! it is HONEST (explicit no-counterpart rows exist, not forced maps), the
//! per-concept / per-standard slices, specific lookups, and deterministic JSON.

use fs_crosswalk::{
    CROSSWALK_VERSION, Counterpart, PackageConcept, SUPPORTED_PACKAGE_FORMAT, Standard, audit,
    crosswalk, for_concept, for_standard, lookup, to_json,
};

#[test]
fn compatibility_versions_are_explicit() {
    assert_eq!(CROSSWALK_VERSION, 5);
    assert_eq!(SUPPORTED_PACKAGE_FORMAT, 8);
}

#[test]
fn the_table_covers_every_concept_by_every_standard() {
    assert_eq!(crosswalk().len(), 12 * 4);
    let a = audit();
    assert_eq!(a.expected, 48);
    assert!(a.ok(), "no silent gaps: {:?}", a.gaps);
    assert_eq!(a.mapped + a.no_counterpart, 48);
}

#[test]
fn the_crosswalk_is_honest_about_missing_counterparts() {
    let a = audit();
    // it does not force a false mapping: some fields have no named counterpart.
    assert!(a.no_counterpart > 0);
    // specifically, content-addressed integrity has no V&V-10 counterpart...
    assert!(matches!(
        lookup(PackageConcept::MerkleRoot, Standard::AsmeVvV10)
            .unwrap()
            .counterpart,
        Counterpart::NoCounterpart { .. }
    ));
    // ...while a certified bound DOES map to solution verification.
    assert!(matches!(
        lookup(PackageConcept::VerifiedColor, Standard::AsmeVvV10)
            .unwrap()
            .counterpart,
        Counterpart::Mapped { .. }
    ));
    // Schema-v6 origin traceability maps only where the vocabulary supports
    // it; an authenticated waiver is not laundered into scientific evidence.
    assert!(matches!(
        lookup(PackageConcept::ClaimOrigin, Standard::AsmeVvV10)
            .unwrap()
            .counterpart,
        Counterpart::NoCounterpart { .. }
    ));
    assert!(matches!(
        lookup(PackageConcept::ClaimOrigin, Standard::AsmeVvV40)
            .unwrap()
            .counterpart,
        Counterpart::Mapped { .. }
    ));
    assert!(matches!(
        lookup(PackageConcept::WaiverAuthorization, Standard::AsmeVvV40)
            .unwrap()
            .counterpart,
        Counterpart::NoCounterpart { .. }
    ));
    assert!(matches!(
        lookup(PackageConcept::WaiverAuthorization, Standard::FaaEasaCbA)
            .unwrap()
            .counterpart,
        Counterpart::Mapped { .. }
    ));
}

#[test]
fn every_concept_maps_across_all_four_standards() {
    for c in PackageConcept::ALL {
        assert_eq!(for_concept(c).len(), 4, "{:?}", c.label());
    }
}

#[test]
fn every_standard_covers_all_twelve_concepts() {
    for s in Standard::ALL {
        assert_eq!(for_standard(s).len(), 12, "{:?}", s.label());
        assert!(!s.full_name().is_empty());
    }
}

#[test]
fn certificate_rows_explicitly_cover_portable_witnesses_without_claiming_verification() {
    for standard in Standard::ALL {
        let entry = lookup(PackageConcept::Certificate, standard)
            .expect("every standard has one certificate row");
        let Counterpart::Mapped { note, .. } = entry.counterpart else {
            panic!("certificate unexpectedly lacks a counterpart for {standard:?}");
        };
        assert!(note.contains("portable witness"), "{standard:?}: {note}");
        assert!(
            note.contains("does not establish semantic re-verification"),
            "{standard:?}: {note}"
        );
    }
}

#[test]
fn a_validated_claim_maps_to_the_validation_metric() {
    let e = lookup(PackageConcept::ValidatedColor, Standard::AsmeVvV20).unwrap();
    assert!(matches!(
        e.counterpart,
        Counterpart::Mapped { clause, .. } if clause.contains("Validation")
    ));
}

#[test]
fn signature_mapping_does_not_claim_checker_release_admission() {
    for standard in Standard::ALL {
        let entry = lookup(PackageConcept::Signature, standard).expect("signature row");
        let note = match entry.counterpart {
            Counterpart::Mapped { note, .. } => note,
            Counterpart::NoCounterpart { .. } => "",
        };
        assert!(
            !note.is_empty()
                && (note.contains("checker admission") || note.contains("checker gate decision")),
            "signature relation must preserve the no-admission boundary: {note}"
        );
    }
}

#[test]
fn labels_are_unique() {
    let concept_labels: Vec<&str> = PackageConcept::ALL.iter().map(|c| c.label()).collect();
    let mut sorted = concept_labels.clone();
    sorted.sort_unstable();
    sorted.dedup();
    assert_eq!(sorted.len(), concept_labels.len());

    let std_labels: Vec<&str> = Standard::ALL.iter().map(|s| s.label()).collect();
    let mut ss = std_labels.clone();
    ss.sort_unstable();
    ss.dedup();
    assert_eq!(ss.len(), std_labels.len());
}

#[test]
fn json_is_well_formed_and_deterministic() {
    let j = to_json();
    assert_eq!(j, to_json());
    assert!(j.starts_with('{') && j.ends_with('}'));
    assert!(j.contains(&format!("\"version\":{CROSSWALK_VERSION}")));
    assert!(j.contains(&format!(
        "\"supported_package_format\":{SUPPORTED_PACKAGE_FORMAT}"
    )));
    assert_eq!(j.matches("\"concept\":").count(), 48);
    assert!(j.contains("verified-color") && j.contains("asme-vv-40"));
    assert!(j.contains("claim-origin") && j.contains("waiver-authorization"));
    assert!(j.contains("no_counterpart"));
    assert!(!j.contains(",,"));
}
