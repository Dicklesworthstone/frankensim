//! Accelerator pilot doctrine and backend evidence-class registry (bead
//! `frankensim-extreal-program-f85xj.15.1`).
//!
//! This is descriptive governance data. It defines the evidence a conditional
//! `[M]` accelerator pilot would have to retain and the measured falsifier that
//! can refuse the pilot before any backend dependency or device implementation
//! is admitted. It does not provide an accelerator backend, authenticate
//! runtime evidence, or authorize a production dependency.

use core::fmt;

/// Version of the accelerator doctrine and backend evidence-class schema.
pub const ACCELERATOR_DOCTRINE_SCHEMA_VERSION: u16 = 1;

/// Bead that owns the doctrine.
pub const ACCELERATOR_DOCTRINE_BEAD: &str = "frankensim-extreal-program-f85xj.15.1";
/// Bead that must supply the end-to-end wall-time and energy profile.
pub const ACCELERATOR_PROFILE_BEAD: &str = "frankensim-extreal-program-f85xj.15.2";
/// Conditional single-kernel pilot governed by this doctrine.
pub const ACCELERATOR_PILOT_BEAD: &str = "frankensim-extreal-program-f85xj.15.3";
/// Production-dependency ruling required before a backend is selected.
pub const ACCELERATOR_DEPENDENCY_POLICY_BEAD: &str = "frankensim-extreal-program-f85xj.11.1";
/// Fixed-size moonshot displacement policy that governs any admitted pilot.
pub const ACCELERATOR_MOONSHOT_POLICY_BEAD: &str = "frankensim-extreal-program-f85xj.16.3";

/// Minimum aggregate end-to-end wall-time share of the profiled top three.
pub const MIN_TOP_THREE_WALL_SHARE_BPS: u16 = 5_000;
/// Minimum aggregate energy share where credible energy data is available.
pub const MIN_TOP_THREE_ENERGY_SHARE_BPS: u16 = 5_000;
/// Minimum end-to-end wall-time share of the kernel selected for a pilot.
pub const MIN_PILOT_KERNEL_WALL_SHARE_BPS: u16 = 1_500;

/// Stable no-authority boundary carried by rendered doctrine artifacts.
pub const ACCELERATOR_DOCTRINE_NO_CLAIM: &str = "governance schema only; no accelerator backend, dependency admission, device execution, speedup, energy saving, numerical equivalence, cancellation completion, or production authority is established";

/// Pilot ambition class.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AcceleratorAmbition {
    /// Research-grade work, off the default product path until proven.
    Moonshot,
}

impl AcceleratorAmbition {
    /// Stable plan tag.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::Moonshot => "[M]",
        }
    }
}

/// Whether an evidence field maps to a live or deliberately absent type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvidenceRecordStatus {
    /// The named type exists at the source locator.
    Existing,
    /// The named type is required but does not exist yet.
    ExplicitlyNew,
}

impl EvidenceRecordStatus {
    /// Stable machine code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::Existing => "existing",
            Self::ExplicitlyNew => "explicitly-new",
        }
    }
}

/// Named measured falsifier for the conditional accelerator pilot.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AcceleratorFalsifier {
    /// Observation that the profiling campaign must retain.
    pub observation: &'static str,
    /// Deterministic decision rule applied to that observation.
    pub decision_rule: &'static str,
    /// Terminal result when the rule refuses a pilot.
    pub refusal: &'static str,
}

/// Complete v1 accelerator governance record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AcceleratorDoctrine {
    /// Owning doctrine Bead.
    pub doctrine_bead: &'static str,
    /// Required profiling Bead.
    pub profile_bead: &'static str,
    /// Conditional pilot Bead.
    pub pilot_bead: &'static str,
    /// Required production-dependency admission ruling.
    pub dependency_policy_bead: &'static str,
    /// Required moonshot displacement policy.
    pub moonshot_policy_bead: &'static str,
    /// Pilot ambition class.
    pub ambition: AcceleratorAmbition,
    /// Preconditions that must all hold before a pilot may open.
    pub pilot_gate: &'static str,
    /// Named falsifier.
    pub falsifier: AcceleratorFalsifier,
    /// Whether the CPU path must remain a retained reference.
    pub permanent_cpu_reference: bool,
    /// Smallest device-work boundary at which cancellation must be observed.
    pub cancellation_boundary: &'static str,
    /// Minimum aggregate top-three wall-time share, in basis points.
    pub min_top_three_wall_share_bps: u16,
    /// Minimum aggregate top-three energy share, in basis points.
    pub min_top_three_energy_share_bps: u16,
    /// Minimum selected-kernel wall-time share, in basis points.
    pub min_pilot_kernel_wall_share_bps: u16,
    /// Explicit authority boundary.
    pub no_claim: &'static str,
}

/// One required field in a future per-run accelerator evidence record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct BackendEvidenceField {
    /// Canonical field id.
    pub id: &'static str,
    /// Human-readable evidence field.
    pub field: &'static str,
    /// Existing or proposed Rust record type.
    pub record_type: &'static str,
    /// Workspace crate responsible for the record.
    pub record_owner: &'static str,
    /// Honest current record maturity.
    pub status: EvidenceRecordStatus,
    /// Live implementation locator; absent for new records.
    pub source_locator: Option<&'static str>,
    /// What a pilot must retain.
    pub requirement: &'static str,
    /// What this field alone cannot establish.
    pub no_claim: &'static str,
}

/// One kernel family eligible for profiling, not preselected for a pilot.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AcceleratorCandidate {
    /// Canonical candidate id.
    pub id: &'static str,
    /// Candidate kernel family.
    pub kernel_family: &'static str,
    /// Current CPU implementation owner.
    pub crate_name: &'static str,
    /// Repository-relative CPU implementation locator.
    pub source_locator: &'static str,
    /// Why profiling may find accelerator leverage.
    pub suitability_hypothesis: &'static str,
    /// Known reason profiling may falsify that hypothesis.
    pub known_hazard: &'static str,
}

/// Canonical v1 doctrine.
pub const ACCELERATOR_DOCTRINE: AcceleratorDoctrine = AcceleratorDoctrine {
    doctrine_bead: ACCELERATOR_DOCTRINE_BEAD,
    profile_bead: ACCELERATOR_PROFILE_BEAD,
    pilot_bead: ACCELERATOR_PILOT_BEAD,
    dependency_policy_bead: ACCELERATOR_DEPENDENCY_POLICY_BEAD,
    moonshot_policy_bead: ACCELERATOR_MOONSHOT_POLICY_BEAD,
    ambition: AcceleratorAmbition::Moonshot,
    pilot_gate: "the end-to-end profile must satisfy the named falsifier, a separate production-dependency ruling must admit the chosen runtime, and the fixed-size moonshot portfolio must admit a displacement before one feature-gated kernel pilot opens",
    falsifier: AcceleratorFalsifier {
        observation: "profile representative workflows on both reference ISA families; rank kernel wall-time and, where credible platform data exists, energy share after retaining import, assembly, boundary setup, reporting, ledger I/O, transfer, and synchronization costs",
        decision_rule: "refuse when the top three kernels jointly account for less than 50% of end-to-end wall time, jointly account for less than 50% of measured energy where energy is available, or contain no transfer-and-synchronization-suitable candidate with at least 15% of end-to-end wall time",
        refusal: "close the conditional pilot as refused-with-evidence; retain the CPU path and do not admit an accelerator dependency",
    },
    permanent_cpu_reference: true,
    cancellation_boundary: "kernel-batch boundary followed by request, drain, and finalize evidence",
    min_top_three_wall_share_bps: MIN_TOP_THREE_WALL_SHARE_BPS,
    min_top_three_energy_share_bps: MIN_TOP_THREE_ENERGY_SHARE_BPS,
    min_pilot_kernel_wall_share_bps: MIN_PILOT_KERNEL_WALL_SHARE_BPS,
    no_claim: ACCELERATOR_DOCTRINE_NO_CLAIM,
};

/// Closed backend evidence-class table in canonical schema order.
pub const BACKEND_EVIDENCE_FIELDS: [BackendEvidenceField; 12] = [
    BackendEvidenceField {
        id: "AE-01",
        field: "device identity",
        record_type: "AcceleratorEnvironmentReceipt",
        record_owner: "fs-roofline",
        status: EvidenceRecordStatus::ExplicitlyNew,
        source_locator: None,
        requirement: "vendor, architecture, model, stable device identifier, memory topology, and capability fingerprint",
        no_claim: "a device name does not establish which binary ran or whether results are equivalent",
    },
    BackendEvidenceField {
        id: "AE-02",
        field: "driver and runtime version",
        record_type: "AcceleratorEnvironmentReceipt",
        record_owner: "fs-roofline",
        status: EvidenceRecordStatus::ExplicitlyNew,
        source_locator: None,
        requirement: "exact driver, userspace runtime, backend API, and enabled feature versions",
        no_claim: "version strings do not authenticate the loaded implementation",
    },
    BackendEvidenceField {
        id: "AE-03",
        field: "compiler and backend identity",
        record_type: "AcceleratorEnvironmentReceipt",
        record_owner: "fs-roofline",
        status: EvidenceRecordStatus::ExplicitlyNew,
        source_locator: None,
        requirement: "compiler, code-generation backend, flags, target features, build identity, and dependency decision receipt",
        no_claim: "compiler metadata does not prove semantic preservation or production admission",
    },
    BackendEvidenceField {
        id: "AE-04",
        field: "kernel source identity",
        record_type: "fs_blake3::ContentHash",
        record_owner: "fs-blake3",
        status: EvidenceRecordStatus::Existing,
        source_locator: Some("crates/fs-blake3/src/lib.rs"),
        requirement: "content identity of canonical kernel source and generated or embedded device binary",
        no_claim: "a content hash identifies bytes but does not authenticate their author or execution",
    },
    BackendEvidenceField {
        id: "AE-05",
        field: "kernel and workload identity",
        record_type: "fs_roofline::KernelSpec",
        record_owner: "fs-roofline",
        status: EvidenceRecordStatus::Existing,
        source_locator: Some("crates/fs-roofline/src/lib.rs"),
        requirement: "versioned kernel, dimensions, units, arithmetic-intensity model, dataset, phase, and denominator",
        no_claim: "a static kernel specification is not a timed production run",
    },
    BackendEvidenceField {
        id: "AE-06",
        field: "CPU machine identity and measured axes",
        record_type: "fs_roofline::MachineAxes",
        record_owner: "fs-roofline",
        status: EvidenceRecordStatus::Existing,
        source_locator: Some("crates/fs-roofline/src/axes.rs"),
        requirement: "CPU topology fingerprint plus measured bandwidth and compute axes before interpreting a comparison",
        no_claim: "CPU axes do not describe accelerator throughput, transfer cost, or energy",
    },
    BackendEvidenceField {
        id: "AE-07",
        field: "reduction and determinism policy",
        record_type: "AcceleratorReductionPolicyReceipt",
        record_owner: "fs-roofline",
        status: EvidenceRecordStatus::ExplicitlyNew,
        source_locator: None,
        requirement: "fixed-order reduction topology or a versioned tolerance policy, including accumulation type and tie breaking",
        no_claim: "declaring a policy does not prove that a device kernel followed it",
    },
    BackendEvidenceField {
        id: "AE-08",
        field: "CPU-versus-device numerical comparison",
        record_type: "AcceleratorEquivalenceReceipt",
        record_owner: "fs-evidence",
        status: EvidenceRecordStatus::ExplicitlyNew,
        source_locator: None,
        requirement: "per-QoI CPU and device values with units, numerical uncertainty, acceptance envelope, and corpus identity",
        no_claim: "no exact equivalence receipt exists today, and one future accepted metric would not establish equivalence outside its corpus, QoI, or envelope",
    },
    BackendEvidenceField {
        id: "AE-09",
        field: "permanent CPU reference",
        record_type: "fs_roofline::RecordedProductionRun and FreshProductionEvidence",
        record_owner: "fs-roofline",
        status: EvidenceRecordStatus::Existing,
        source_locator: Some("crates/fs-roofline/src/production.rs"),
        requirement: "retain the CPU implementation and exact comparison run; revalidate freshness whenever it is cited positively",
        no_claim: "a recorded operation is not fresh positive evidence until live authority revalidation succeeds",
    },
    BackendEvidenceField {
        id: "AE-10",
        field: "cancellation, drain, and finalize outcome",
        record_type: "fs_exec::DrainFinalizeReport",
        record_owner: "fs-exec",
        status: EvidenceRecordStatus::Existing,
        source_locator: Some("crates/fs-exec/src/cx.rs"),
        requirement: "executor-minted proof that every admitted kernel batch observed request, drained, and finalized",
        no_claim: "a caller-authored cancelled flag or dropped device handle is not drain evidence",
    },
    BackendEvidenceField {
        id: "AE-11",
        field: "phase wall-time and energy attribution",
        record_type: "PipelineAttributionReceipt",
        record_owner: "fs-roofline",
        status: EvidenceRecordStatus::ExplicitlyNew,
        source_locator: None,
        requirement: "end-to-end phase totals, kernel shares, transfer and synchronization costs, profiling overhead, energy provenance, and named gaps",
        no_claim: "a kernel microbenchmark cannot stand in for workflow-level speed or energy benefit",
    },
    BackendEvidenceField {
        id: "AE-12",
        field: "go or no-go decision",
        record_type: "AcceleratorPilotDecisionReceipt",
        record_owner: "fs-govern",
        status: EvidenceRecordStatus::ExplicitlyNew,
        source_locator: None,
        requirement: "bind the exact profile, thresholds, top-three ranking, suitability findings, dependency ruling, moonshot displacement, and terminal decision",
        no_claim: "a governance decision does not prove scientific correctness or future production fitness",
    },
];

/// Closed candidate-family table. Profiling, not this order, selects a pilot.
pub const ACCELERATOR_CANDIDATES: [AcceleratorCandidate; 5] = [
    AcceleratorCandidate {
        id: "AK-01",
        kernel_family: "D3Q19 LBM collide and stream",
        crate_name: "fs-lbm",
        source_locator: "crates/fs-lbm/src/d3q19/sparse.rs",
        suitability_hypothesis: "regular active-tile arithmetic may expose enough parallel work per launch",
        known_hazard: "halo traffic, sparse activity, boundary handling, and publication barriers may dominate",
    },
    AcceleratorCandidate {
        id: "AK-02",
        kernel_family: "sparse matrix-vector multiplication",
        crate_name: "fs-sparse",
        source_locator: "crates/fs-sparse/src/lib.rs",
        suitability_hypothesis: "large repeated sparse operators may benefit from device memory bandwidth",
        known_hazard: "irregular gathers and problem sizes may remain bandwidth-bound and transfer-dominated",
    },
    AcceleratorCandidate {
        id: "AK-03",
        kernel_family: "FFT and batched transforms",
        crate_name: "fs-fft",
        source_locator: "crates/fs-fft/src/lib.rs",
        suitability_hypothesis: "regular batched transforms can offer substantial parallel work",
        known_hazard: "pencil transposes, full-array passes, and host-device movement may erase kernel gains",
    },
    AcceleratorCandidate {
        id: "AK-04",
        kernel_family: "batched constitutive evaluation",
        crate_name: "fs-material",
        source_locator: "crates/fs-material/src/lib.rs",
        suitability_hypothesis: "many independent material points may form efficient batches",
        known_hazard: "no dedicated batch API exists today; branching laws, internal-state traffic, and small batches may underfill the device",
    },
    AcceleratorCandidate {
        id: "AK-05",
        kernel_family: "spectral path tracing",
        crate_name: "fs-render",
        source_locator: "crates/fs-render/src/tracer.rs",
        suitability_hypothesis: "independent pixel samples provide abundant parallelism",
        known_hazard: "divergent paths, deterministic accumulation, scene transfer, and feature-gated maturity may dominate",
    },
];

/// Fail-closed accelerator-doctrine validation error.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AcceleratorDoctrineError {
    /// The caller presented an unsupported schema version.
    SchemaVersion {
        /// Presented schema version.
        found: u16,
    },
    /// The closed evidence-field count changed.
    EvidenceCount {
        /// Presented field count.
        found: usize,
    },
    /// The closed candidate count changed.
    CandidateCount {
        /// Presented candidate count.
        found: usize,
    },
    /// A canonical evidence id changed or moved.
    EvidenceId {
        /// Zero-based row position.
        index: usize,
        /// Required id.
        expected: &'static str,
        /// Presented id.
        found: &'static str,
    },
    /// A canonical candidate id changed or moved.
    CandidateId {
        /// Zero-based row position.
        index: usize,
        /// Required id.
        expected: &'static str,
        /// Presented id.
        found: &'static str,
    },
    /// A required doctrine or row field was empty.
    EmptyField {
        /// Doctrine or row id.
        row: &'static str,
        /// Empty field.
        field: &'static str,
    },
    /// A workspace owner did not use the `fs-*` convention.
    InvalidCrateName {
        /// Row id.
        row: &'static str,
        /// Rejected crate name.
        crate_name: &'static str,
    },
    /// An existing record omitted its implementation locator.
    ExistingRecordMissingLocator {
        /// Evidence row id.
        row: &'static str,
    },
    /// A deliberately new record falsely named a live implementation.
    NewRecordHasLocator {
        /// Evidence row id.
        row: &'static str,
    },
    /// A source locator was empty, absolute, or escaped the repository.
    InvalidSourceLocator {
        /// Row id.
        row: &'static str,
        /// Rejected locator.
        locator: &'static str,
    },
    /// A basis-point threshold was invalid or internally inconsistent.
    InvalidThresholds,
    /// A non-moonshot pilot or removable CPU reference was presented.
    InvalidPilotBoundary,
    /// A valid-looking doctrine field changed without a schema bump.
    DoctrineDrift,
    /// A valid-looking evidence row changed without a schema bump.
    EvidenceRowDrift {
        /// Canonical row id.
        row: &'static str,
    },
    /// A valid-looking candidate row changed without a schema bump.
    CandidateRowDrift {
        /// Canonical row id.
        row: &'static str,
    },
}

impl fmt::Display for AcceleratorDoctrineError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for AcceleratorDoctrineError {}

const CANONICAL_EVIDENCE_IDS: [&str; 12] = [
    "AE-01", "AE-02", "AE-03", "AE-04", "AE-05", "AE-06", "AE-07", "AE-08", "AE-09", "AE-10",
    "AE-11", "AE-12",
];
const CANONICAL_CANDIDATE_IDS: [&str; 5] = ["AK-01", "AK-02", "AK-03", "AK-04", "AK-05"];

fn locator_is_valid(locator: &str) -> bool {
    !locator.trim().is_empty()
        && !locator.starts_with('/')
        && !locator.split('/').any(|component| component == "..")
}

/// Validate a candidate doctrine, evidence table, and candidate registry.
///
/// Mutation and migration tooling can use this to prove that a proposed v1
/// table fails before reaching documentation or a future runtime adapter.
pub fn validate_accelerator_doctrine(
    schema_version: u16,
    doctrine: &AcceleratorDoctrine,
    evidence_fields: &[BackendEvidenceField],
    candidates: &[AcceleratorCandidate],
) -> Result<(), AcceleratorDoctrineError> {
    if schema_version != ACCELERATOR_DOCTRINE_SCHEMA_VERSION {
        return Err(AcceleratorDoctrineError::SchemaVersion {
            found: schema_version,
        });
    }
    if evidence_fields.len() != BACKEND_EVIDENCE_FIELDS.len() {
        return Err(AcceleratorDoctrineError::EvidenceCount {
            found: evidence_fields.len(),
        });
    }
    if candidates.len() != ACCELERATOR_CANDIDATES.len() {
        return Err(AcceleratorDoctrineError::CandidateCount {
            found: candidates.len(),
        });
    }
    for (field, value) in [
        ("doctrine_bead", doctrine.doctrine_bead),
        ("profile_bead", doctrine.profile_bead),
        ("pilot_bead", doctrine.pilot_bead),
        ("dependency_policy_bead", doctrine.dependency_policy_bead),
        ("moonshot_policy_bead", doctrine.moonshot_policy_bead),
        ("pilot_gate", doctrine.pilot_gate),
        ("falsifier.observation", doctrine.falsifier.observation),
        ("falsifier.decision_rule", doctrine.falsifier.decision_rule),
        ("falsifier.refusal", doctrine.falsifier.refusal),
        ("cancellation_boundary", doctrine.cancellation_boundary),
        ("no_claim", doctrine.no_claim),
    ] {
        if value.trim().is_empty() {
            return Err(AcceleratorDoctrineError::EmptyField {
                row: "accelerator-doctrine",
                field,
            });
        }
    }
    if doctrine.ambition != AcceleratorAmbition::Moonshot || !doctrine.permanent_cpu_reference {
        return Err(AcceleratorDoctrineError::InvalidPilotBoundary);
    }
    let thresholds = [
        doctrine.min_top_three_wall_share_bps,
        doctrine.min_top_three_energy_share_bps,
        doctrine.min_pilot_kernel_wall_share_bps,
    ];
    if thresholds
        .iter()
        .any(|threshold| *threshold == 0 || *threshold > 10_000)
        || doctrine.min_pilot_kernel_wall_share_bps > doctrine.min_top_three_wall_share_bps
    {
        return Err(AcceleratorDoctrineError::InvalidThresholds);
    }

    for (index, row) in evidence_fields.iter().enumerate() {
        if row.id != CANONICAL_EVIDENCE_IDS[index] {
            return Err(AcceleratorDoctrineError::EvidenceId {
                index,
                expected: CANONICAL_EVIDENCE_IDS[index],
                found: row.id,
            });
        }
        for (field, value) in [
            ("field", row.field),
            ("record_type", row.record_type),
            ("record_owner", row.record_owner),
            ("requirement", row.requirement),
            ("no_claim", row.no_claim),
        ] {
            if value.trim().is_empty() {
                return Err(AcceleratorDoctrineError::EmptyField { row: row.id, field });
            }
        }
        if !row.record_owner.starts_with("fs-") {
            return Err(AcceleratorDoctrineError::InvalidCrateName {
                row: row.id,
                crate_name: row.record_owner,
            });
        }
        match (row.status, row.source_locator) {
            (EvidenceRecordStatus::Existing, None) => {
                return Err(AcceleratorDoctrineError::ExistingRecordMissingLocator { row: row.id });
            }
            (EvidenceRecordStatus::ExplicitlyNew, Some(_)) => {
                return Err(AcceleratorDoctrineError::NewRecordHasLocator { row: row.id });
            }
            (_, Some(locator)) if !locator_is_valid(locator) => {
                return Err(AcceleratorDoctrineError::InvalidSourceLocator {
                    row: row.id,
                    locator,
                });
            }
            _ => {}
        }
        if *row != BACKEND_EVIDENCE_FIELDS[index] {
            return Err(AcceleratorDoctrineError::EvidenceRowDrift { row: row.id });
        }
    }

    for (index, row) in candidates.iter().enumerate() {
        if row.id != CANONICAL_CANDIDATE_IDS[index] {
            return Err(AcceleratorDoctrineError::CandidateId {
                index,
                expected: CANONICAL_CANDIDATE_IDS[index],
                found: row.id,
            });
        }
        for (field, value) in [
            ("kernel_family", row.kernel_family),
            ("crate_name", row.crate_name),
            ("source_locator", row.source_locator),
            ("suitability_hypothesis", row.suitability_hypothesis),
            ("known_hazard", row.known_hazard),
        ] {
            if value.trim().is_empty() {
                return Err(AcceleratorDoctrineError::EmptyField { row: row.id, field });
            }
        }
        if !row.crate_name.starts_with("fs-") {
            return Err(AcceleratorDoctrineError::InvalidCrateName {
                row: row.id,
                crate_name: row.crate_name,
            });
        }
        if !locator_is_valid(row.source_locator) {
            return Err(AcceleratorDoctrineError::InvalidSourceLocator {
                row: row.id,
                locator: row.source_locator,
            });
        }
        if *row != ACCELERATOR_CANDIDATES[index] {
            return Err(AcceleratorDoctrineError::CandidateRowDrift { row: row.id });
        }
    }
    if *doctrine != ACCELERATOR_DOCTRINE {
        return Err(AcceleratorDoctrineError::DoctrineDrift);
    }
    Ok(())
}

/// Validate and return the canonical accelerator doctrine.
pub fn accelerator_doctrine() -> Result<&'static AcceleratorDoctrine, AcceleratorDoctrineError> {
    validate_accelerator_doctrine(
        ACCELERATOR_DOCTRINE_SCHEMA_VERSION,
        &ACCELERATOR_DOCTRINE,
        &BACKEND_EVIDENCE_FIELDS,
        &ACCELERATOR_CANDIDATES,
    )?;
    Ok(&ACCELERATOR_DOCTRINE)
}

fn push_json_string(output: &mut String, value: &str) {
    output.push('"');
    for ch in value.chars() {
        match ch {
            '"' => output.push_str("\\\""),
            '\\' => output.push_str("\\\\"),
            '\n' => output.push_str("\\n"),
            '\r' => output.push_str("\\r"),
            '\t' => output.push_str("\\t"),
            ch if ch <= '\u{1f}' => {
                use fmt::Write as _;
                let _ = write!(output, "\\u{:04x}", ch as u32);
            }
            ch => output.push(ch),
        }
    }
    output.push('"');
}

fn push_json_optional_string(output: &mut String, value: Option<&str>) {
    if let Some(value) = value {
        push_json_string(output, value);
    } else {
        output.push_str("null");
    }
}

/// Deterministic machine-readable accelerator doctrine.
pub fn accelerator_doctrine_json() -> Result<String, AcceleratorDoctrineError> {
    let doctrine = accelerator_doctrine()?;
    let mut output = String::from("{\"schema_version\":");
    output.push_str(&ACCELERATOR_DOCTRINE_SCHEMA_VERSION.to_string());
    for (name, value) in [
        ("doctrine_bead", doctrine.doctrine_bead),
        ("profile_bead", doctrine.profile_bead),
        ("pilot_bead", doctrine.pilot_bead),
        ("dependency_policy_bead", doctrine.dependency_policy_bead),
        ("moonshot_policy_bead", doctrine.moonshot_policy_bead),
        ("ambition", doctrine.ambition.code()),
        ("pilot_gate", doctrine.pilot_gate),
    ] {
        output.push_str(",\"");
        output.push_str(name);
        output.push_str("\":");
        push_json_string(&mut output, value);
    }
    output.push_str(",\"falsifier\":{\"observation\":");
    push_json_string(&mut output, doctrine.falsifier.observation);
    output.push_str(",\"decision_rule\":");
    push_json_string(&mut output, doctrine.falsifier.decision_rule);
    output.push_str(",\"refusal\":");
    push_json_string(&mut output, doctrine.falsifier.refusal);
    output.push_str("},\"permanent_cpu_reference\":");
    output.push_str(if doctrine.permanent_cpu_reference {
        "true"
    } else {
        "false"
    });
    output.push_str(",\"cancellation_boundary\":");
    push_json_string(&mut output, doctrine.cancellation_boundary);
    output.push_str(",\"thresholds_bps\":{\"top_three_wall\":");
    output.push_str(&doctrine.min_top_three_wall_share_bps.to_string());
    output.push_str(",\"top_three_energy\":");
    output.push_str(&doctrine.min_top_three_energy_share_bps.to_string());
    output.push_str(",\"pilot_kernel_wall\":");
    output.push_str(&doctrine.min_pilot_kernel_wall_share_bps.to_string());
    output.push_str("},\"no_claim\":");
    push_json_string(&mut output, doctrine.no_claim);
    output.push_str(",\"evidence_fields\":[");
    for (index, row) in BACKEND_EVIDENCE_FIELDS.iter().enumerate() {
        if index != 0 {
            output.push(',');
        }
        output.push_str("{\"id\":");
        push_json_string(&mut output, row.id);
        for (name, value) in [
            ("field", row.field),
            ("record_type", row.record_type),
            ("record_owner", row.record_owner),
            ("status", row.status.code()),
        ] {
            output.push_str(",\"");
            output.push_str(name);
            output.push_str("\":");
            push_json_string(&mut output, value);
        }
        output.push_str(",\"source_locator\":");
        push_json_optional_string(&mut output, row.source_locator);
        output.push_str(",\"requirement\":");
        push_json_string(&mut output, row.requirement);
        output.push_str(",\"no_claim\":");
        push_json_string(&mut output, row.no_claim);
        output.push('}');
    }
    output.push_str("],\"candidates\":[");
    for (index, row) in ACCELERATOR_CANDIDATES.iter().enumerate() {
        if index != 0 {
            output.push(',');
        }
        output.push_str("{\"id\":");
        push_json_string(&mut output, row.id);
        for (name, value) in [
            ("kernel_family", row.kernel_family),
            ("crate", row.crate_name),
            ("source_locator", row.source_locator),
            ("suitability_hypothesis", row.suitability_hypothesis),
            ("known_hazard", row.known_hazard),
        ] {
            output.push_str(",\"");
            output.push_str(name);
            output.push_str("\":");
            push_json_string(&mut output, value);
        }
        output.push('}');
    }
    output.push_str("]}");
    Ok(output)
}

/// Code-derived Markdown block embedded in `docs/ACCELERATOR_DOCTRINE.md`.
pub fn accelerator_doctrine_markdown() -> Result<String, AcceleratorDoctrineError> {
    let doctrine = accelerator_doctrine()?;
    let mut output = format!(
        "| Policy field | Canonical value |\n\
         | --- | --- |\n\
         | Ambition | `{}` |\n\
         | Profiling gate | `{}` |\n\
         | Conditional pilot | `{}` |\n\
         | Dependency ruling | `{}` |\n\
         | Moonshot displacement | `{}` |\n\
         | Thresholds | top-three wall: {:.1}%; top-three energy where measured: {:.1}%; selected kernel wall: {:.1}% |\n\
         | Cancellation | {} |\n\
         | Permanent CPU reference | required |\n\
         | Named falsifier | {} |\n\
         | Refusal | {} |\n\
         | No claim | {} |\n\n",
        doctrine.ambition.code(),
        doctrine.profile_bead,
        doctrine.pilot_bead,
        doctrine.dependency_policy_bead,
        doctrine.moonshot_policy_bead,
        f64::from(doctrine.min_top_three_wall_share_bps) / 100.0,
        f64::from(doctrine.min_top_three_energy_share_bps) / 100.0,
        f64::from(doctrine.min_pilot_kernel_wall_share_bps) / 100.0,
        doctrine.cancellation_boundary,
        doctrine.falsifier.decision_rule,
        doctrine.falsifier.refusal,
        doctrine.no_claim,
    );
    output.push_str(
        "| ID | Required field | Record mapping | Status | Requirement and boundary |\n\
         | --- | --- | --- | --- | --- |\n",
    );
    for row in BACKEND_EVIDENCE_FIELDS {
        output.push_str(&format!(
            "| `{}` | {} | `{}` in `{}`<br>{} | `{}` | Requirement: {}<br>No claim: {} |\n",
            row.id,
            row.field,
            row.record_type,
            row.record_owner,
            row.source_locator.unwrap_or("no implementation locator"),
            row.status.code(),
            row.requirement,
            row.no_claim,
        ));
    }
    output.push('\n');
    output.push_str(
        "| ID | Candidate family | CPU source | Profiling hypothesis | Known falsifier pressure |\n\
         | --- | --- | --- | --- | --- |\n",
    );
    for row in ACCELERATOR_CANDIDATES {
        output.push_str(&format!(
            "| `{}` | {} | `{}`<br>{} | {} | {} |\n",
            row.id,
            row.kernel_family,
            row.crate_name,
            row.source_locator,
            row.suitability_hypothesis,
            row.known_hazard,
        ));
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    const DOCTRINE_DOC: &str = include_str!("../../../docs/ACCELERATOR_DOCTRINE.md");
    const BLOCK_BEGIN: &str = "<!-- BEGIN CODE-DERIVED ACCELERATOR DOCTRINE -->";
    const BLOCK_END: &str = "<!-- END CODE-DERIVED ACCELERATOR DOCTRINE -->";

    #[test]
    fn g0_canonical_doctrine_is_complete_and_conditional() {
        let doctrine = accelerator_doctrine().expect("canonical doctrine validates");
        assert_eq!(doctrine.ambition.code(), "[M]");
        assert_eq!(doctrine.profile_bead, ACCELERATOR_PROFILE_BEAD);
        assert_eq!(doctrine.pilot_bead, ACCELERATOR_PILOT_BEAD);
        assert!(doctrine.permanent_cpu_reference);
        assert!(doctrine.pilot_gate.contains("must admit"));
        assert!(doctrine.falsifier.decision_rule.contains("top three"));
        assert!(doctrine.falsifier.refusal.contains("refused-with-evidence"));
        assert!(doctrine.cancellation_boundary.contains("drain"));
        assert_eq!(BACKEND_EVIDENCE_FIELDS.len(), 12);
        assert_eq!(ACCELERATOR_CANDIDATES.len(), 5);
    }

    #[test]
    fn g0_schema_count_and_identity_mutants_fail_closed() {
        assert_eq!(
            validate_accelerator_doctrine(
                ACCELERATOR_DOCTRINE_SCHEMA_VERSION + 1,
                &ACCELERATOR_DOCTRINE,
                &BACKEND_EVIDENCE_FIELDS,
                &ACCELERATOR_CANDIDATES,
            ),
            Err(AcceleratorDoctrineError::SchemaVersion {
                found: ACCELERATOR_DOCTRINE_SCHEMA_VERSION + 1,
            })
        );
        assert_eq!(
            validate_accelerator_doctrine(
                ACCELERATOR_DOCTRINE_SCHEMA_VERSION,
                &ACCELERATOR_DOCTRINE,
                &BACKEND_EVIDENCE_FIELDS[..11],
                &ACCELERATOR_CANDIDATES,
            ),
            Err(AcceleratorDoctrineError::EvidenceCount { found: 11 })
        );
        let mut evidence = BACKEND_EVIDENCE_FIELDS;
        evidence[0].id = "AE-00";
        assert_eq!(
            validate_accelerator_doctrine(
                ACCELERATOR_DOCTRINE_SCHEMA_VERSION,
                &ACCELERATOR_DOCTRINE,
                &evidence,
                &ACCELERATOR_CANDIDATES,
            ),
            Err(AcceleratorDoctrineError::EvidenceId {
                index: 0,
                expected: "AE-01",
                found: "AE-00",
            })
        );

        let mut candidates = ACCELERATOR_CANDIDATES;
        candidates[1].known_hazard = "valid-looking but semantically different hazard";
        assert_eq!(
            validate_accelerator_doctrine(
                ACCELERATOR_DOCTRINE_SCHEMA_VERSION,
                &ACCELERATOR_DOCTRINE,
                &BACKEND_EVIDENCE_FIELDS,
                &candidates,
            ),
            Err(AcceleratorDoctrineError::CandidateRowDrift { row: "AK-02" })
        );
    }

    #[test]
    fn g0_record_maturity_cannot_imply_missing_or_future_source() {
        let mut evidence = BACKEND_EVIDENCE_FIELDS;
        evidence[3].source_locator = None;
        assert_eq!(
            validate_accelerator_doctrine(
                ACCELERATOR_DOCTRINE_SCHEMA_VERSION,
                &ACCELERATOR_DOCTRINE,
                &evidence,
                &ACCELERATOR_CANDIDATES,
            ),
            Err(AcceleratorDoctrineError::ExistingRecordMissingLocator { row: "AE-04" })
        );
        let mut evidence = BACKEND_EVIDENCE_FIELDS;
        evidence[0].source_locator = Some("crates/fs-roofline/src/lib.rs");
        assert_eq!(
            validate_accelerator_doctrine(
                ACCELERATOR_DOCTRINE_SCHEMA_VERSION,
                &ACCELERATOR_DOCTRINE,
                &evidence,
                &ACCELERATOR_CANDIDATES,
            ),
            Err(AcceleratorDoctrineError::NewRecordHasLocator { row: "AE-01" })
        );
        let mut doctrine = ACCELERATOR_DOCTRINE;
        doctrine.permanent_cpu_reference = false;
        assert_eq!(
            validate_accelerator_doctrine(
                ACCELERATOR_DOCTRINE_SCHEMA_VERSION,
                &doctrine,
                &BACKEND_EVIDENCE_FIELDS,
                &ACCELERATOR_CANDIDATES,
            ),
            Err(AcceleratorDoctrineError::InvalidPilotBoundary)
        );

        let mut doctrine = ACCELERATOR_DOCTRINE;
        doctrine.min_pilot_kernel_wall_share_bps += 1;
        assert_eq!(
            validate_accelerator_doctrine(
                ACCELERATOR_DOCTRINE_SCHEMA_VERSION,
                &doctrine,
                &BACKEND_EVIDENCE_FIELDS,
                &ACCELERATOR_CANDIDATES,
            ),
            Err(AcceleratorDoctrineError::DoctrineDrift)
        );

        let mut evidence = BACKEND_EVIDENCE_FIELDS;
        evidence[4].requirement = "valid-looking but semantically different workload identity";
        assert_eq!(
            validate_accelerator_doctrine(
                ACCELERATOR_DOCTRINE_SCHEMA_VERSION,
                &ACCELERATOR_DOCTRINE,
                &evidence,
                &ACCELERATOR_CANDIDATES,
            ),
            Err(AcceleratorDoctrineError::EvidenceRowDrift { row: "AE-05" })
        );
    }

    #[test]
    fn g0_named_existing_sources_and_candidate_sources_exist() {
        let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
        let repository = manifest_dir
            .parent()
            .and_then(Path::parent)
            .expect("fs-govern lives under repository/crates");
        for row in BACKEND_EVIDENCE_FIELDS {
            assert!(
                repository.join("crates").join(row.record_owner).is_dir(),
                "{} names missing workspace crate {}",
                row.id,
                row.record_owner
            );
            if let Some(locator) = row.source_locator {
                assert!(repository.join(locator).is_file());
            } else {
                assert_eq!(row.status, EvidenceRecordStatus::ExplicitlyNew);
            }
        }
        for row in ACCELERATOR_CANDIDATES {
            assert!(repository.join("crates").join(row.crate_name).is_dir());
            assert!(
                repository.join(row.source_locator).is_file(),
                "{} names missing CPU source {}",
                row.id,
                row.source_locator
            );
        }
    }

    #[test]
    fn g5_machine_catalog_is_deterministic_and_complete() {
        let first = accelerator_doctrine_json().expect("canonical doctrine renders");
        assert_eq!(
            first,
            accelerator_doctrine_json().expect("canonical doctrine rerenders")
        );
        assert!(first.starts_with("{\"schema_version\":1,\"doctrine_bead\":"));
        assert!(first.contains("\"ambition\":\"[M]\""));
        assert!(first.contains("\"permanent_cpu_reference\":true"));
        assert!(first.contains("\"top_three_wall\":5000"));
        for row in BACKEND_EVIDENCE_FIELDS {
            assert!(first.contains(&format!("\"id\":\"{}\"", row.id)));
            assert!(first.contains(&format!("\"record_type\":\"{}\"", row.record_type)));
        }
        for row in ACCELERATOR_CANDIDATES {
            assert!(first.contains(&format!("\"id\":\"{}\"", row.id)));
        }
    }

    #[test]
    fn g0_documented_doctrine_is_exactly_code_derived() {
        let start = DOCTRINE_DOC
            .find(BLOCK_BEGIN)
            .expect("doctrine block begin marker")
            + BLOCK_BEGIN.len();
        let tail = &DOCTRINE_DOC[start..];
        let end = tail.find(BLOCK_END).expect("doctrine block end marker");
        let documented = &tail[..end];
        let generated = accelerator_doctrine_markdown().expect("doctrine renders");
        assert_eq!(documented, format!("\n{generated}"));
        assert!(DOCTRINE_DOC.contains("does not ship an accelerator backend"));
        assert!(DOCTRINE_DOC.contains("refused-with-evidence"));
        assert!(DOCTRINE_DOC.contains("workflow-level"));
    }
}
