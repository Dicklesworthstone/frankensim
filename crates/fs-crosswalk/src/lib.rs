//! fs-crosswalk — the regulatory vocabulary crosswalk (plan addendum,
//! Proposal 12). Layer: UTIL (pure data + audit; no dependencies).
//!
//! An evidence package ([`fs-package`](../fs_package)) is only persuasive to an
//! auditor if it speaks the auditor's language. This crate is the machine-
//! readable CROSSWALK that maps every evidence-package field onto the
//! regulator's EXISTING vocabulary — ASME V&V 10 (computational solid
//! mechanics), V&V 20 (computational fluid dynamics & heat transfer), V&V 40
//! (credibility of in-silico medical-device models), and FAA/EASA
//! certification-by-analysis (CbA) — so the artifact lands in the standards'
//! own terms instead of demanding the regulator learn ours.
//!
//! It OWNS risk R9 (standards-body latency): because every field is either
//! mapped to a named clause or EXPLICITLY flagged as having no counterpart, the
//! package doubles as internal-QA and B2B due-diligence collateral, so the
//! investment is not stranded if standards bodies move slowly. [`audit`]
//! enforces the "no silent gap" rule — every `(concept, standard)` pair is
//! covered.
//!
//! The mapping is versioned ([`CROSSWALK_VERSION`]) alongside the package
//! format. NOTE: the mappings below are a first-party reading of the standards
//! for engineering use; they are not an official ASME/FAA/EASA determination.

/// The crosswalk version (moves with the evidence-package format).
pub const CROSSWALK_VERSION: u32 = 1;

/// An evidence-package field/concept (the FrankenSim side of the map).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PackageConcept {
    /// A verified-color claim (a-posteriori certified numerical error bound).
    VerifiedColor,
    /// A validated-color claim (matched to data within a regime).
    ValidatedColor,
    /// An estimated-color claim (a conjecture with dispersion).
    EstimatedColor,
    /// The raw certificate payload behind a claim.
    Certificate,
    /// A falsifier log (adversarial attempts / negative results).
    FalsifierLog,
    /// A regime / validity-domain tag.
    RegimeTag,
    /// The anchoring (experimental / reference) dataset.
    AnchoringDataset,
    /// Provenance (code version + dependency lockfile).
    Provenance,
    /// The Merkle content-address root.
    MerkleRoot,
    /// A detached attestation signature.
    Signature,
}

impl PackageConcept {
    /// Every concept, in order.
    pub const ALL: [PackageConcept; 10] = [
        PackageConcept::VerifiedColor,
        PackageConcept::ValidatedColor,
        PackageConcept::EstimatedColor,
        PackageConcept::Certificate,
        PackageConcept::FalsifierLog,
        PackageConcept::RegimeTag,
        PackageConcept::AnchoringDataset,
        PackageConcept::Provenance,
        PackageConcept::MerkleRoot,
        PackageConcept::Signature,
    ];

    /// A stable slug.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            PackageConcept::VerifiedColor => "verified-color",
            PackageConcept::ValidatedColor => "validated-color",
            PackageConcept::EstimatedColor => "estimated-color",
            PackageConcept::Certificate => "certificate",
            PackageConcept::FalsifierLog => "falsifier-log",
            PackageConcept::RegimeTag => "regime-tag",
            PackageConcept::AnchoringDataset => "anchoring-dataset",
            PackageConcept::Provenance => "provenance",
            PackageConcept::MerkleRoot => "merkle-root",
            PackageConcept::Signature => "signature",
        }
    }
}

/// A target certification standard (the regulator's side of the map).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Standard {
    /// ASME V&V 10 — computational solid mechanics.
    AsmeVvV10,
    /// ASME V&V 20 — computational fluid dynamics & heat transfer.
    AsmeVvV20,
    /// ASME V&V 40 — credibility of in-silico medical-device models.
    AsmeVvV40,
    /// FAA/EASA certification-by-analysis.
    FaaEasaCbA,
}

impl Standard {
    /// Every standard, in order.
    pub const ALL: [Standard; 4] = [
        Standard::AsmeVvV10,
        Standard::AsmeVvV20,
        Standard::AsmeVvV40,
        Standard::FaaEasaCbA,
    ];

    /// A stable slug.
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Standard::AsmeVvV10 => "asme-vv-10",
            Standard::AsmeVvV20 => "asme-vv-20",
            Standard::AsmeVvV40 => "asme-vv-40",
            Standard::FaaEasaCbA => "faa-easa-cba",
        }
    }

    /// The full standard name.
    #[must_use]
    pub const fn full_name(self) -> &'static str {
        match self {
            Standard::AsmeVvV10 => {
                "ASME V&V 10 (Verification & Validation in Computational Solid Mechanics)"
            }
            Standard::AsmeVvV20 => {
                "ASME V&V 20 (V&V in Computational Fluid Dynamics & Heat Transfer)"
            }
            Standard::AsmeVvV40 => {
                "ASME V&V 40 (Credibility of Computational Modeling for Medical Devices)"
            }
            Standard::FaaEasaCbA => "FAA/EASA Certification by Analysis",
        }
    }
}

/// Whether a concept has a named counterpart in a standard.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Counterpart {
    /// A named clause/concept in the standard.
    Mapped {
        /// The named clause / concept.
        clause: &'static str,
        /// How the mapping reads.
        note: &'static str,
    },
    /// No counterpart — flagged explicitly (honesty, not a silent gap).
    NoCounterpart {
        /// Why there is no counterpart.
        reason: &'static str,
    },
}

/// One crosswalk row: a `(concept, standard)` pair and its counterpart.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CrosswalkEntry {
    /// The evidence-package concept.
    pub concept: PackageConcept,
    /// The target standard.
    pub standard: Standard,
    /// The mapping (or explicit absence).
    pub counterpart: Counterpart,
}

impl CrosswalkEntry {
    /// Is this a positive mapping (vs an explicit no-counterpart)?
    #[must_use]
    pub const fn is_mapped(&self) -> bool {
        matches!(self.counterpart, Counterpart::Mapped { .. })
    }
}

const fn m(
    concept: PackageConcept,
    standard: Standard,
    clause: &'static str,
    note: &'static str,
) -> CrosswalkEntry {
    CrosswalkEntry {
        concept,
        standard,
        counterpart: Counterpart::Mapped { clause, note },
    }
}

const fn nc(concept: PackageConcept, standard: Standard, reason: &'static str) -> CrosswalkEntry {
    CrosswalkEntry {
        concept,
        standard,
        counterpart: Counterpart::NoCounterpart { reason },
    }
}

use PackageConcept::{
    AnchoringDataset, Certificate, EstimatedColor, FalsifierLog, MerkleRoot, Provenance, RegimeTag,
    Signature, ValidatedColor, VerifiedColor,
};
use Standard::{AsmeVvV10, AsmeVvV20, AsmeVvV40, FaaEasaCbA};

/// The canonical crosswalk: every `(concept, standard)` pair (10 × 4).
const CROSSWALK: [CrosswalkEntry; 40] = [
    // Verified color — certified numerical error bounds = SOLUTION VERIFICATION.
    m(
        VerifiedColor,
        AsmeVvV10,
        "Solution (calculation) verification",
        "an a-posteriori certified error bound is exactly numerical-error estimation",
    ),
    m(
        VerifiedColor,
        AsmeVvV20,
        "Numerical uncertainty u_num",
        "the certified bound quantifies u_num directly (a guaranteed, not estimated, bound)",
    ),
    m(
        VerifiedColor,
        AsmeVvV40,
        "Verification activities in the credibility assessment",
        "code + calculation verification evidence",
    ),
    m(
        VerifiedColor,
        FaaEasaCbA,
        "Numerical accuracy substantiation of the analysis method",
        "the certified bound is the accuracy substantiation",
    ),
    // Validated color — matched to data within a regime = VALIDATION.
    m(
        ValidatedColor,
        AsmeVvV10,
        "Validation (prediction vs experiment)",
        "a validated claim is a validation result with its referent data",
    ),
    m(
        ValidatedColor,
        AsmeVvV20,
        "Validation & the validation metric (E = S - D)",
        "regime + dataset provide S, D and the comparison",
    ),
    m(
        ValidatedColor,
        AsmeVvV40,
        "Validation within the context of use",
        "validated color = validated for its declared applicability",
    ),
    m(
        ValidatedColor,
        FaaEasaCbA,
        "Validation of the computational method against test evidence",
        "the CbA validation step",
    ),
    // Estimated color — a conjecture with dispersion = UNCERTAINTY.
    m(
        EstimatedColor,
        AsmeVvV10,
        "Uncertainty quantification (model/input uncertainty)",
        "an estimated claim carries model-form uncertainty, not a bound",
    ),
    m(
        EstimatedColor,
        AsmeVvV20,
        "Estimated input/model uncertainties (u_input)",
        "estimated color = uncertainty not reduced to a certified bound",
    ),
    m(
        EstimatedColor,
        AsmeVvV40,
        "Model uncertainty in the credibility assessment",
        "lowest credibility rung pending further evidence",
    ),
    m(
        EstimatedColor,
        FaaEasaCbA,
        "Analysis uncertainty / declared conservatism",
        "estimated results carry explicit conservatism assumptions",
    ),
    // Certificate — the raw bound payload.
    m(
        Certificate,
        AsmeVvV10,
        "Numerical error estimate (e.g. GCI)",
        "the certificate is the numerical-error estimate record",
    ),
    m(
        Certificate,
        AsmeVvV20,
        "u_num record (grid-convergence / residual basis)",
        "the raw interval is the u_num evidence",
    ),
    m(
        Certificate,
        AsmeVvV40,
        "Verification output evidence",
        "the certificate is the verification artifact",
    ),
    m(
        Certificate,
        FaaEasaCbA,
        "Documented accuracy metric",
        "the certificate is the reported accuracy figure",
    ),
    // Falsifier log — adversarial negative results.
    nc(
        FalsifierLog,
        AsmeVvV10,
        "V&V 10 names no adversarial-falsification artifact; sensitivity analysis is the nearest activity",
    ),
    nc(
        FalsifierLog,
        AsmeVvV20,
        "the validation-uncertainty framework does not name negative/falsification logs",
    ),
    m(
        FalsifierLog,
        AsmeVvV40,
        "Model-risk evidence within credibility (where the model fails)",
        "falsifier logs are direct model-risk evidence",
    ),
    m(
        FalsifierLog,
        FaaEasaCbA,
        "Documented limitations & applicability boundaries",
        "falsifier results bound the method's applicability",
    ),
    // Regime tag — the validity domain.
    m(
        RegimeTag,
        AsmeVvV10,
        "Domain (range) of validation",
        "the regime tag is the validated domain",
    ),
    m(
        RegimeTag,
        AsmeVvV20,
        "Domain of validation",
        "regime bounds define the validated domain",
    ),
    m(
        RegimeTag,
        AsmeVvV40,
        "Context of use (COU) / applicability",
        "the regime is the model's context of use",
    ),
    m(
        RegimeTag,
        FaaEasaCbA,
        "Validation envelope / applicability range",
        "regime = the certified analysis envelope",
    ),
    // Anchoring dataset — the reference/experimental data.
    m(
        AnchoringDataset,
        AsmeVvV10,
        "Experimental (validation referent) data",
        "the anchoring dataset is the validation referent",
    ),
    m(
        AnchoringDataset,
        AsmeVvV20,
        "Experimental data D and its uncertainty u_D",
        "the dataset provides D, u_D",
    ),
    m(
        AnchoringDataset,
        AsmeVvV40,
        "Comparator / reference dataset",
        "the anchoring dataset is the comparator",
    ),
    m(
        AnchoringDataset,
        FaaEasaCbA,
        "Substantiating test evidence",
        "the dataset is the substantiating test data",
    ),
    // Provenance — code version + lockfile.
    m(
        Provenance,
        AsmeVvV10,
        "Software quality / code-verification records",
        "provenance is the SQA/config record",
    ),
    m(
        Provenance,
        AsmeVvV20,
        "Simulation documentation & configuration control",
        "provenance is the configuration record",
    ),
    m(
        Provenance,
        AsmeVvV40,
        "Configuration management & traceability",
        "provenance is the traceability record",
    ),
    m(
        Provenance,
        FaaEasaCbA,
        "Tool qualification & configuration control records",
        "provenance supports tool qualification",
    ),
    // Merkle root — content-addressed integrity (a modern software concept).
    nc(
        MerkleRoot,
        AsmeVvV10,
        "content-addressed hash integrity is not a named V&V 10 concept (records integrity is procedural)",
    ),
    nc(
        MerkleRoot,
        AsmeVvV20,
        "no named data-integrity hash concept in V&V 20",
    ),
    nc(
        MerkleRoot,
        AsmeVvV40,
        "no named cryptographic-integrity concept; records management is procedural",
    ),
    nc(
        MerkleRoot,
        FaaEasaCbA,
        "hash-based data integrity is not a named CbA concept (data integrity is handled procedurally)",
    ),
    // Signature — attestation.
    m(
        Signature,
        AsmeVvV10,
        "Responsible-party review & approval",
        "a detached signature is the modern form of the responsible-engineer sign-off",
    ),
    m(
        Signature,
        AsmeVvV20,
        "Documented review/approval of results",
        "signature = authorized approval",
    ),
    m(
        Signature,
        AsmeVvV40,
        "Approval within the credibility assessment plan",
        "signature = plan sign-off",
    ),
    m(
        Signature,
        FaaEasaCbA,
        "Authorized sign-off / statement of compliance",
        "signature attests the means-of-compliance record",
    ),
];

/// The full crosswalk.
#[must_use]
pub fn crosswalk() -> &'static [CrosswalkEntry] {
    &CROSSWALK
}

/// All rows for one package concept.
#[must_use]
pub fn for_concept(concept: PackageConcept) -> Vec<&'static CrosswalkEntry> {
    CROSSWALK.iter().filter(|e| e.concept == concept).collect()
}

/// All rows for one standard.
#[must_use]
pub fn for_standard(standard: Standard) -> Vec<&'static CrosswalkEntry> {
    CROSSWALK
        .iter()
        .filter(|e| e.standard == standard)
        .collect()
}

/// The single row for a `(concept, standard)` pair, if present.
#[must_use]
pub fn lookup(concept: PackageConcept, standard: Standard) -> Option<&'static CrosswalkEntry> {
    CROSSWALK
        .iter()
        .find(|e| e.concept == concept && e.standard == standard)
}

/// The completeness audit of the crosswalk.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CrosswalkAudit {
    /// Total `(concept, standard)` pairs expected.
    pub expected: usize,
    /// Rows with a named mapping.
    pub mapped: usize,
    /// Rows explicitly flagged as having no counterpart.
    pub no_counterpart: usize,
    /// `(concept, standard)` pairs with NO row at all (the failure mode).
    pub gaps: Vec<(PackageConcept, Standard)>,
}

impl CrosswalkAudit {
    /// Is every pair covered (mapped or explicitly flagged)?
    #[must_use]
    pub fn ok(&self) -> bool {
        self.gaps.is_empty() && self.mapped + self.no_counterpart == self.expected
    }
}

/// Audit the crosswalk: every `(concept, standard)` pair must have exactly one
/// row (mapped or an explicit no-counterpart) — no silent gaps.
#[must_use]
pub fn audit() -> CrosswalkAudit {
    let mut mapped = 0;
    let mut no_counterpart = 0;
    let mut gaps = Vec::new();
    for &concept in &PackageConcept::ALL {
        for &standard in &Standard::ALL {
            match lookup(concept, standard) {
                Some(e) if e.is_mapped() => mapped += 1,
                Some(_) => no_counterpart += 1,
                None => gaps.push((concept, standard)),
            }
        }
    }
    CrosswalkAudit {
        expected: PackageConcept::ALL.len() * Standard::ALL.len(),
        mapped,
        no_counterpart,
        gaps,
    }
}

/// Emit the crosswalk as deterministic machine-readable JSON.
#[must_use]
pub fn to_json() -> String {
    use core::fmt::Write as _;
    let mut out = format!("{{\"version\":{CROSSWALK_VERSION},\"entries\":[");
    for (i, e) in CROSSWALK.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        let (kind, a, b) = match e.counterpart {
            Counterpart::Mapped { clause, note } => ("mapped", clause, note),
            Counterpart::NoCounterpart { reason } => ("no_counterpart", reason, ""),
        };
        write!(
            out,
            "{{\"concept\":\"{}\",\"standard\":\"{}\",\"kind\":\"{}\",\"clause\":\"{}\",\"note\":\"{}\"}}",
            e.concept.label(),
            e.standard.label(),
            kind,
            json_escape(a),
            json_escape(b),
        )
        .expect("write to String");
    }
    out.push_str("]}");
    out
}

/// Minimal JSON string escaping.
fn json_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            c => out.push(c),
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Coverage from ACTUAL package fields (bead qmao.6.1): the static table
// above maps concepts to standards; whether a concept is EVIDENCED must
// come from a parsed package, never from the mapping itself — a mapped
// concept with no evidence in hand is "mapped but absent", not covered.
// ---------------------------------------------------------------------------

use fs_evidence::{Color, ColorRank};
use fs_package::EvidencePackage;

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

/// Judge every concept against the fields actually present in `pkg`.
#[must_use]
pub fn package_presence(pkg: &EvidencePackage) -> Vec<ConceptPresence> {
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
                    let ok = !pkg.claims.is_empty() && pkg.verify().is_ok();
                    (
                        ok,
                        if ok {
                            "every claim carries its complete color payload (schema v2)"
                                .to_string()
                        } else {
                            "no claims, or a claim failed completeness verification".to_string()
                        },
                    )
                }
                PackageConcept::FalsifierLog => (
                    false,
                    "falsifier logs are not representable in package schema v2 — coverage \
                     cannot be claimed for absent evidence (recorded no-claim)"
                        .to_string(),
                ),
                PackageConcept::RegimeTag => {
                    let n = validated_with(&|c| {
                        matches!(c, Color::Validated { regime, .. } if !regime.bounds().is_empty())
                    });
                    (n > 0, format!("{n} validated claim(s) with regime bounds"))
                }
                PackageConcept::AnchoringDataset => {
                    let n = validated_with(&|c| {
                        matches!(c, Color::Validated { dataset, .. } if !dataset.trim().is_empty())
                    });
                    (n > 0, format!("{n} validated claim(s) with anchoring datasets"))
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
                    (true, format!("content root {:016x} recomputable", pkg.merkle_root()))
                }
                PackageConcept::Signature => match &pkg.signature {
                    Some(_) => (true, "detached signature present".to_string()),
                    None => (false, "no detached signature attached".to_string()),
                },
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
