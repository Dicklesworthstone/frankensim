//! Package-grounded crosswalk coverage (bead qmao.6.1): which
//! regulatory concepts are EVIDENCED by the fields actually present in
//! a package — the static fs-crosswalk mapping alone can never produce
//! "covered".

// ---------------------------------------------------------------------------
// Coverage from ACTUAL package fields (bead qmao.6.1; moved here from
// fs-crosswalk under the layer rule — L6 may know UTIL, never the
// reverse): the static fs-crosswalk table maps concepts to standards; whether a concept is EVIDENCED must
// come from a parsed package, never from the mapping itself — a mapped
// concept with no evidence in hand is "mapped but absent", not covered.
// ---------------------------------------------------------------------------

use crate::EvidencePackage;
use fs_crosswalk::{PackageConcept, Standard, crosswalk};
use fs_evidence::{Color, ColorRank};

/// Whether one concept is actually evidenced by a package, with the
/// reason either way.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ConceptPresence {
    /// The concept.
    pub concept: PackageConcept,
    /// Is it evidenced by fields ACTUALLY PRESENT in this package?
    pub present: bool,
    /// Why / why not (teaching, deterministic).
    pub why: String,
}

fn anchoring_dataset_presence(pkg: &EvidencePackage) -> (bool, String) {
    let named = pkg
        .claims
        .iter()
        .filter(|claim| {
            matches!(&claim.color, Color::Validated { dataset, .. } if !dataset.trim().is_empty())
        })
        .count();
    let anchored: usize = pkg.claims.iter().map(|claim| claim.anchors.len()).sum();
    let matched = pkg
        .claims
        .iter()
        .filter(|claim| claim.has_matching_validated_anchor())
        .count();
    (
        matched > 0,
        format!(
            "{matched} validated claim(s) have matching content-hashed anchors; {named} named \
             dataset(s), {anchored} total anchor record(s) (schema v3)"
        ),
    )
}

fn signature_presence(pkg: &EvidencePackage) -> (bool, String) {
    match &pkg.signature {
        Some(_) => (
            false,
            "detached signature present but not authenticated; raw presence cannot establish \
             sign-off"
                .to_string(),
        ),
        None => (false, "no authenticated detached signature".to_string()),
    }
}

/// Judge every concept against the fields actually present in `pkg`.
#[must_use]
pub fn package_presence(pkg: &EvidencePackage) -> Vec<ConceptPresence> {
    if pkg.verify().is_err() {
        let why = "package verification failed; fail-closed coverage suppresses every concept"
            .to_string();
        return PackageConcept::ALL
            .iter()
            .map(|&concept| ConceptPresence {
                concept,
                present: false,
                why: why.clone(),
            })
            .collect();
    }

    let count = |rank: ColorRank| pkg.claims.iter().filter(|c| c.color.rank() == rank).count();
    let validated_with = |f: &dyn Fn(&Color) -> bool| {
        pkg.claims
            .iter()
            .filter(|c| matches!(c.color, Color::Validated { .. }) && f(&c.color))
            .count()
    };
    PackageConcept::ALL
        .iter()
        .map(|&concept| {
            let (present, why) = match concept {
                PackageConcept::VerifiedColor => {
                    let n = count(ColorRank::Verified);
                    (n > 0, format!("{n} verified claim(s) with interval certificates"))
                }
                PackageConcept::ValidatedColor => {
                    let n = count(ColorRank::Validated);
                    (n > 0, format!("{n} validated claim(s)"))
                }
                PackageConcept::EstimatedColor => {
                    let n = count(ColorRank::Estimated);
                    (n > 0, format!("{n} estimated claim(s)"))
                }
                PackageConcept::Certificate => {
                    let ok = !pkg.claims.is_empty();
                    (
                        ok,
                        if ok {
                            "every claim carries a complete, re-verified color payload (schema v3)"
                                .to_string()
                        } else {
                            "no claims to certify".to_string()
                        },
                    )
                }
                PackageConcept::FalsifierLog => {
                    let n: usize = pkg.claims.iter().map(|c| c.falsifiers.len()).sum();
                    (
                        n > 0,
                        if n > 0 {
                            format!("{n} falsifier record(s) attached to claims (schema v3)")
                        } else {
                            "no falsifier records attached — coverage cannot be claimed for \
                             absent evidence"
                                .to_string()
                        },
                    )
                }
                PackageConcept::RegimeTag => {
                    let n = validated_with(&|c| {
                        matches!(c, Color::Validated { regime, .. } if !regime.bounds().is_empty())
                    });
                    (n > 0, format!("{n} validated claim(s) with regime bounds"))
                }
                PackageConcept::AnchoringDataset => {
                    anchoring_dataset_presence(pkg)
                }
                PackageConcept::Provenance => {
                    let ok = !pkg.provenance.code_version.trim().is_empty()
                        && !pkg.provenance.constellation_lock.trim().is_empty();
                    (
                        ok,
                        if ok {
                            "code version + constellation lock present".to_string()
                        } else {
                            "provenance fields empty".to_string()
                        },
                    )
                }
                PackageConcept::MerkleRoot => {
                    (true, format!("content root {} recomputable", pkg.merkle_root()))
                }
                PackageConcept::Signature => signature_presence(pkg),
            };
            ConceptPresence {
                concept,
                present,
                why,
            }
        })
        .collect()
}

/// One row of the package-grounded coverage report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum CoverageStatus {
    /// Mapped to the standard AND evidenced by this package.
    Covered,
    /// Mapped to the standard but the evidence is ABSENT from this
    /// package — never reported as covered.
    MappedButAbsent,
    /// Deliberately unmapped for this standard (see the table's reason).
    NoClaim,
}

/// Coverage of `standard` derived from the fields actually present in
/// `pkg` (bead qmao.6.1): the static mapping alone can never produce
/// `Covered` — only the intersection with parsed package evidence can.
#[must_use]
pub fn package_coverage(
    pkg: &EvidencePackage,
    standard: Standard,
) -> Vec<(PackageConcept, CoverageStatus, String)> {
    let presence = package_presence(pkg);
    crosswalk()
        .iter()
        .filter(|e| e.standard == standard)
        .map(|e| {
            let p = presence
                .iter()
                .find(|pr| pr.concept == e.concept)
                .expect("every concept judged");
            let status = if !e.is_mapped() {
                CoverageStatus::NoClaim
            } else if p.present {
                CoverageStatus::Covered
            } else {
                CoverageStatus::MappedButAbsent
            };
            (e.concept, status, p.why.clone())
        })
        .collect()
}
