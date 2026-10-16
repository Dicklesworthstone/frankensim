//! CAE Ecosystem format capability matrix and adapter ruling registry.
//!
//! Bead: `frankensim-extreal-program-f85xj.11.5`
//!
//! Under Decision Record ADPT-2026-07 (bead f85xj.11.1), external CAE formats
//! and legacy kernel interactions operate strictly through the quarantine boundary.
//! Foreign tools (OpenCASCADE, Gmsh binary, Exodus HDF5) ship as isolated out-of-process
//! adapters that mint Estimate-only authority; pure-core parsers mint Verified/Validated
//! authority upon passing full repair and validity checks.

/// Evidence and epistemic classification for CAE format tiers.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CaeEvidenceClass {
    /// Pure interval-certified numerics.
    Verified,
    /// Anchored to experimental or validation datasets in a declared regime.
    Validated,
    /// Estimated via out-of-process adapter, surrogate, or heuristic.
    Estimated,
    /// Explicitly waived claim.
    Waived,
}

/// Processing direction for CAE ecosystem formats.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CaeDirection {
    /// Reading external format into quarantined representation.
    Import,
    /// Emitting external format from simulation outputs.
    Export,
    /// Both import and export supported.
    Bidirectional,
}

/// Trust and quarantine status of a CAE integration tier.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CaeQuarantineStatus {
    /// Pure Rust, zero-dependency parser/writer inside the workspace.
    NativeCertified,
    /// Official out-of-process quarantined adapter binary (ADPT-2026-07).
    QuarantinedAdapter,
    /// Staged for future ruling or evaluation.
    StagedRuling,
}

/// Capability entry for a CAE format tier.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CaeCapabilityEntry {
    /// Format identifier name.
    pub format_id: &'static str,
    /// Standard file extensions.
    pub extensions: &'static [&'static str],
    /// Direction supported.
    pub direction: CaeDirection,
    /// Quarantine and dependency status.
    pub quarantine_status: CaeQuarantineStatus,
    /// Authoritative decision record ruling.
    pub decision_record: &'static str,
    /// Default evidence class minted on import/export.
    pub evidence_class: CaeEvidenceClass,
    /// Explicit scope and no-claim statement.
    pub no_claim_boundary: &'static str,
}

/// Authoritative registry of CAE ecosystem interchange tiers.
pub const CAE_CAPABILITY_MATRIX: &[CaeCapabilityEntry] = &[
    CaeCapabilityEntry {
        format_id: "Gmsh-MSH",
        extensions: &["msh"],
        direction: CaeDirection::Bidirectional,
        quarantine_status: CaeQuarantineStatus::NativeCertified,
        decision_record: "ADPT-2026-07 / fs-io::gmsh",
        evidence_class: CaeEvidenceClass::Verified,
        no_claim_boundary: "Parses and emits MSH 2.2/4.1 geometry and unstructured mesh topology; does not execute Gmsh scripting or meshing algorithms",
    },
    CaeCapabilityEntry {
        format_id: "Abaqus-INP",
        extensions: &["inp"],
        direction: CaeDirection::Import,
        quarantine_status: CaeQuarantineStatus::NativeCertified,
        decision_record: "ADPT-2026-07 / fs-io::inp_bdf",
        evidence_class: CaeEvidenceClass::Validated,
        no_claim_boundary: "Extracts thermal, solid geometry, sets, and boundary cards; unsupported dynamic/plasticity cards recorded in census; never claims Abaqus execution semantics",
    },
    CaeCapabilityEntry {
        format_id: "Nastran-BDF",
        extensions: &["bdf", "dat", "nas"],
        direction: CaeDirection::Import,
        quarantine_status: CaeQuarantineStatus::NativeCertified,
        decision_record: "ADPT-2026-07 / fs-io::inp_bdf",
        evidence_class: CaeEvidenceClass::Validated,
        no_claim_boundary: "Extracts GRID coordinates, solid elements, MAT4/5 thermal cards; never claims Nastran solver semantics",
    },
    CaeCapabilityEntry {
        format_id: "Arrow-IPC",
        extensions: &["arrow", "feather"],
        direction: CaeDirection::Export,
        quarantine_status: CaeQuarantineStatus::NativeCertified,
        decision_record: "ADPT-2026-07 / fs-io::tabular_export",
        evidence_class: CaeEvidenceClass::Verified,
        no_claim_boundary: "Emits zero-copy binary column buffers for tabular results and Monte Carlo samples with unit headers and BLAKE3 receipts",
    },
    CaeCapabilityEntry {
        format_id: "Exodus-II",
        extensions: &["exo", "e"],
        direction: CaeDirection::Export,
        quarantine_status: CaeQuarantineStatus::QuarantinedAdapter,
        decision_record: "ADPT-2026-07::EXODUS-ADAPTER",
        evidence_class: CaeEvidenceClass::Estimated,
        no_claim_boundary: "Routes through official quarantined adapter; output carries Estimate authority and adapter version in receipt",
    },
    CaeCapabilityEntry {
        format_id: "CGNS",
        extensions: &["cgns"],
        direction: CaeDirection::Export,
        quarantine_status: CaeQuarantineStatus::QuarantinedAdapter,
        decision_record: "ADPT-2026-07::CGNS-ADAPTER",
        evidence_class: CaeEvidenceClass::Estimated,
        no_claim_boundary: "Routes through official quarantined adapter; output carries Estimate authority and adapter version in receipt",
    },
    CaeCapabilityEntry {
        format_id: "OpenFOAM",
        extensions: &["foam"],
        direction: CaeDirection::Export,
        quarantine_status: CaeQuarantineStatus::StagedRuling,
        decision_record: "ADPT-2026-07::STAGED",
        evidence_class: CaeEvidenceClass::Estimated,
        no_claim_boundary: "Staged post-0.1 comparison workflow tier",
    },
    CaeCapabilityEntry {
        format_id: "FMI-FMU",
        extensions: &["fmu"],
        direction: CaeDirection::Export,
        quarantine_status: CaeQuarantineStatus::StagedRuling,
        decision_record: "ADPT-2026-07::STAGED",
        evidence_class: CaeEvidenceClass::Estimated,
        no_claim_boundary: "Staged post-0.1 system coupling tier",
    },
];
