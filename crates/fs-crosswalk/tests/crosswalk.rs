//! Battery for the regulatory crosswalk (addendum Proposal 12). Verifies the
//! table is complete (every concept × standard covered, no silent gaps), that
//! it is HONEST (explicit no-counterpart rows exist, not forced maps), the
//! per-concept / per-standard slices, specific lookups, and deterministic JSON.

use fs_crosswalk::{
    CROSSWALK_VERSION, Counterpart, PackageConcept, Standard, audit, crosswalk, for_concept,
    for_standard, lookup, to_json,
};

#[test]
fn the_table_covers_every_concept_by_every_standard() {
    assert_eq!(crosswalk().len(), 10 * 4);
    let a = audit();
    assert_eq!(a.expected, 40);
    assert!(a.ok(), "no silent gaps: {:?}", a.gaps);
    assert_eq!(a.mapped + a.no_counterpart, 40);
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
}

#[test]
fn every_concept_maps_across_all_four_standards() {
    for c in PackageConcept::ALL {
        assert_eq!(for_concept(c).len(), 4, "{:?}", c.label());
    }
}

#[test]
fn every_standard_covers_all_ten_concepts() {
    for s in Standard::ALL {
        assert_eq!(for_standard(s).len(), 10, "{:?}", s.label());
        assert!(!s.full_name().is_empty());
    }
}

#[test]
fn a_validated_claim_maps_to_the_validation_metric() {
    let e = lookup(PackageConcept::ValidatedColor, Standard::AsmeVvV20).unwrap();
    match e.counterpart {
        Counterpart::Mapped { clause, .. } => assert!(clause.contains("Validation")),
        Counterpart::NoCounterpart { .. } => panic!("validated color should map in V&V 20"),
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
    assert_eq!(j.matches("\"concept\":").count(), 40);
    assert!(j.contains("verified-color") && j.contains("asme-vv-40"));
    assert!(j.contains("no_counterpart"));
    assert!(!j.contains(",,"));
}

/// qmao.6.1 — crosswalk coverage derives from fields ACTUALLY PRESENT
/// in a parsed package: a mapped concept with absent evidence is
/// "mapped but absent", never covered; unrepresentable concepts
/// (falsifier logs in schema v2) can never report covered.
#[test]
fn coverage_cannot_claim_absent_evidence() {
    use fs_crosswalk::{CoverageStatus, PackageConcept, Standard, package_coverage, package_presence};
    use fs_evidence::Color;
    use fs_package::{Claim, EvidencePackage, Provenance};
    // Unsigned, verified-only package: no validated claims, no regime
    // tags, no datasets, no signature.
    let pkg = EvidencePackage::new(Provenance::new("v", "lock")).with_claim(Claim::new(
        "c",
        "bounded",
        Color::Verified { lo: 0.0, hi: 1.0 },
    ));
    let presence = package_presence(&pkg);
    let get = |c: PackageConcept| presence.iter().find(|p| p.concept == c).unwrap();
    assert!(get(PackageConcept::VerifiedColor).present);
    assert!(!get(PackageConcept::ValidatedColor).present);
    assert!(!get(PackageConcept::RegimeTag).present);
    assert!(!get(PackageConcept::AnchoringDataset).present);
    assert!(!get(PackageConcept::Signature).present);
    assert!(
        !get(PackageConcept::FalsifierLog).present,
        "unrepresentable evidence can never read as present"
    );
    for standard in Standard::ALL {
        for (concept, status, _why) in package_coverage(&pkg, standard) {
            if matches!(status, CoverageStatus::Covered) {
                let p = presence.iter().find(|p| p.concept == concept).unwrap();
                assert!(p.present, "{concept:?} covered without evidence for {standard:?}");
            }
            if concept == PackageConcept::FalsifierLog {
                assert!(
                    !matches!(status, CoverageStatus::Covered),
                    "falsifier logs cannot be covered in v2 ({standard:?})"
                );
            }
        }
    }
}
