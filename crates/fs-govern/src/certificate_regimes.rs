//! Certificate-regime routing doctrine (bead `frankensim-extreal-program-f85xj.9.1`).
//!
//! This module fixes a closed, machine-readable claim-to-evidence table. It is
//! descriptive governance data: selecting a row does not mint evidence,
//! authenticate an artifact, or admit a runtime claim. The downstream router
//! is tracked separately by `frankensim-extreal-program-f85xj.9.3`.

use core::fmt;

pub use fs_evidence::{ClaimClass, EvidenceRegime};

/// Schema version for the closed certificate-regime table.
pub const CERTIFICATE_REGIME_SCHEMA_VERSION: u16 = 1;

/// Bead that owns the future executable claim router.
pub const CERTIFICATE_REGIME_ROUTER_BEAD: &str = "frankensim-extreal-program-f85xj.9.3";

/// Stable no-authority boundary carried by every rendered doctrine artifact.
pub const CERTIFICATE_REGIME_NO_CLAIM: &str = "descriptive routing doctrine only; no evidence, scientific truth, artifact authenticity, or runtime admission is minted";

/// Marketing and report-language boundary for wide but honest enclosures.
pub const CERTIFICATE_REPORT_BOUNDARY: &str = "A mathematically valid enclosure remains a valid result even when it is too wide to support the engineering decision. It must be reported as inconclusive or NoUsefulBound for that claim, never promoted to an engineering certificate by changing the adjective.";

/// Current maturity of one serving capability.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CapabilityStatus {
    /// A narrow implementation exists at the named source locator.
    Available,
    /// Some necessary machinery exists, but not the complete claim route.
    Thin,
    /// The capability is deliberately staged and has no implementation locator.
    Staged,
}

impl CapabilityStatus {
    /// Stable machine code.
    #[must_use]
    pub const fn code(self) -> &'static str {
        match self {
            Self::Available => "available",
            Self::Thin => "thin",
            Self::Staged => "staged",
        }
    }
}

/// One crate/capability mapping retained by a doctrine row.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CapabilityRef {
    /// Workspace crate expected to own or stage the capability.
    pub crate_name: &'static str,
    /// Stable capability name.
    pub capability: &'static str,
    /// Honest current maturity.
    pub status: CapabilityStatus,
    /// Repository-relative implementation evidence; absent only when staged.
    pub source_locator: Option<&'static str>,
}

/// One closed claim-to-evidence routing row.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct CertificateRegimeRow {
    /// Canonical row id.
    pub id: &'static str,
    /// Claim class.
    pub claim: ClaimClass,
    /// Required evidence object or explicit no-useful-bound outcome.
    pub evidence: EvidenceRegime,
    /// Exact applicability boundary.
    pub scope: &'static str,
    /// Current crate/capability mapping.
    pub capabilities: &'static [CapabilityRef],
    /// What the row explicitly does not establish.
    pub no_claim: &'static str,
}

const ROOT_CAPABILITIES: &[CapabilityRef] = &[
    CapabilityRef {
        crate_name: "fs-ivl",
        capability: "interval-root-isolation",
        status: CapabilityStatus::Available,
        source_locator: Some("crates/fs-ivl/src/newton.rs"),
    },
    CapabilityRef {
        crate_name: "fs-ivl",
        capability: "univariate-taylor-enclosure",
        status: CapabilityStatus::Available,
        source_locator: Some("crates/fs-ivl/src/taylor.rs"),
    },
];

const REACHABILITY_CAPABILITIES: &[CapabilityRef] = &[CapabilityRef {
    crate_name: "fs-ivl",
    capability: "validated-reachability-tube",
    status: CapabilityStatus::Staged,
    source_locator: None,
}];

const BALANCE_CAPABILITIES: &[CapabilityRef] = &[CapabilityRef {
    crate_name: "fs-evidence",
    capability: "discrete-balance-defect-receipt",
    status: CapabilityStatus::Available,
    source_locator: Some("crates/fs-evidence/src/balance.rs"),
}];

const STABILITY_CAPABILITIES: &[CapabilityRef] = &[
    CapabilityRef {
        crate_name: "fs-spectral",
        capability: "residual-enclosed-spectral-service",
        status: CapabilityStatus::Available,
        source_locator: Some("crates/fs-spectral/src/service.rs"),
    },
    CapabilityRef {
        crate_name: "fs-sos",
        capability: "quadratic-lyapunov-check",
        status: CapabilityStatus::Available,
        source_locator: Some("crates/fs-sos/src/lib.rs"),
    },
];

const MEAN_LOAD_CAPABILITIES: &[CapabilityRef] = &[
    CapabilityRef {
        crate_name: "fs-eproc",
        capability: "time-uniform-mean-evidence",
        status: CapabilityStatus::Available,
        source_locator: Some("crates/fs-eproc/src/lib.rs"),
    },
    CapabilityRef {
        crate_name: "fs-uq",
        capability: "model-uncertainty-context",
        status: CapabilityStatus::Thin,
        source_locator: Some("crates/fs-uq/src/lib.rs"),
    },
];

const SPECTRUM_CAPABILITIES: &[CapabilityRef] = &[
    CapabilityRef {
        crate_name: "fs-spectral",
        capability: "spectral-computation-evidence",
        status: CapabilityStatus::Thin,
        source_locator: Some("crates/fs-spectral/src/service.rs"),
    },
    CapabilityRef {
        crate_name: "fs-uq",
        capability: "distributional-spectrum-validation",
        status: CapabilityStatus::Staged,
        source_locator: None,
    },
];

const RELIABILITY_CAPABILITIES: &[CapabilityRef] = &[
    CapabilityRef {
        crate_name: "fs-eproc",
        capability: "anytime-sequential-evidence",
        status: CapabilityStatus::Available,
        source_locator: Some("crates/fs-eproc/src/lib.rs"),
    },
    CapabilityRef {
        crate_name: "fs-uq",
        capability: "rare-event-duty-cycle-model",
        status: CapabilityStatus::Staged,
        source_locator: None,
    },
];

const NO_USEFUL_BOUND_CAPABILITIES: &[CapabilityRef] = &[CapabilityRef {
    crate_name: "fs-govern",
    capability: "certificate-regime-no-useful-bound-policy",
    status: CapabilityStatus::Available,
    source_locator: Some("crates/fs-govern/src/certificate_regimes.rs"),
}];

/// Closed v1 doctrine table in canonical schema order.
pub const CERTIFICATE_REGIMES: [CertificateRegimeRow; 8] = [
    CertificateRegimeRow {
        id: "CR-01",
        claim: ClaimClass::RootOrEventTime,
        evidence: EvidenceRegime::IntervalRootOrTaylorEnclosure,
        scope: "one declared finite parameter/time domain with explicit function and derivative enclosure assumptions",
        capabilities: ROOT_CAPABILITIES,
        no_claim: "does not establish behavior after the isolated event or outside the admitted box",
    },
    CertificateRegimeRow {
        id: "CR-02",
        claim: ClaimClass::ShortHorizonReachability,
        evidence: EvidenceRegime::ValidatedReachabilityTube,
        scope: "one declared finite horizon, dynamics version, initial set, disturbance set, and event geometry",
        capabilities: REACHABILITY_CAPABILITIES,
        no_claim: "validated reachability tubes are staged; ordinary interval propagation is not a tube certificate",
    },
    CertificateRegimeRow {
        id: "CR-03",
        claim: ClaimClass::ConservedQuantity,
        evidence: EvidenceRegime::DiscreteBalanceCertificate,
        scope: "one discrete operator, mesh/complex, boundary partition, source convention, units, and time slab",
        capabilities: BALANCE_CAPABILITIES,
        no_claim: "a balance receipt does not validate the constitutive model or measurement chain",
    },
    CertificateRegimeRow {
        id: "CR-04",
        claim: ClaimClass::LocalStability,
        evidence: EvidenceRegime::SpectralOrLyapunovCertificate,
        scope: "one stated equilibrium/orbit, model and linearization version, parameter domain, norm, and spectral or Lyapunov assumptions",
        capabilities: STABILITY_CAPABILITIES,
        no_claim: "local stability evidence does not imply global attraction, nonlinear robustness, or long-time predictive accuracy",
    },
    CertificateRegimeRow {
        id: "CR-05",
        claim: ClaimClass::LongHorizonMeanLoad,
        evidence: EvidenceRegime::StatisticalObservableWithModelEvidence,
        scope: "one declared observable, population/regime, sampling design, dependence model, stopping rule, and model-form evidence",
        capabilities: MEAN_LOAD_CAPABILITIES,
        no_claim: "sampling evidence does not prove the simulator model, and a trajectory enclosure is not required or implied",
    },
    CertificateRegimeRow {
        id: "CR-06",
        claim: ClaimClass::BroadbandSpectrum,
        evidence: EvidenceRegime::DistributionalSpectralValidation,
        scope: "one declared spectrum/statistic, frequency band, windowing convention, operating regime, observation process, and comparison metric",
        capabilities: SPECTRUM_CAPABILITIES,
        no_claim: "computed eigenvalue or FFT evidence alone does not validate a turbulent distribution or broadband field",
    },
    CertificateRegimeRow {
        id: "CR-07",
        claim: ClaimClass::DutyCycleReliability,
        evidence: EvidenceRegime::SequentialRareEventStatistics,
        scope: "one declared duty-cycle population, failure event, censoring/dependence model, stopping policy, and model/field evidence",
        capabilities: RELIABILITY_CAPABILITIES,
        no_claim: "anytime validity under the stated sampling law does not validate the failure model or transfer to an undeclared population",
    },
    CertificateRegimeRow {
        id: "CR-08",
        claim: ClaimClass::ExactLongChaoticTrajectory,
        evidence: EvidenceRegime::NoUsefulBound,
        scope: "a request for one exact trajectory materially beyond the admitted predictability horizon of the stated chaotic model",
        capabilities: NO_USEFUL_BOUND_CAPABILITIES,
        no_claim: "NoUsefulBound is an honest routing result, not proof that all statistical, local, event, or finite-horizon claims are impossible",
    },
];

/// A positive, bounded role for interval/Taylor evidence in chaotic systems.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct IntervalRole {
    /// Stable role code.
    pub id: &'static str,
    /// What interval/Taylor evidence can support.
    pub supports: &'static str,
    /// Boundary that prevents overclaiming.
    pub boundary: &'static str,
}

/// Positive interval roles retained even when exact long-trajectory bounds are useless.
pub const INTERVAL_ROLES: [IntervalRole; 5] = [
    IntervalRole {
        id: "IR-01",
        supports: "short-time event and root isolation",
        boundary: "only the declared finite search domain and enclosure assumptions",
    },
    IntervalRole {
        id: "IR-02",
        supports: "local exclusion and collision-free boxes",
        boundary: "only the admitted box, geometry version, and finite horizon",
    },
    IntervalRole {
        id: "IR-03",
        supports: "parameter and constitutive-law bounds",
        boundary: "only the stated parameter domain and law version",
    },
    IntervalRole {
        id: "IR-04",
        supports: "discrete conservation and balance audits",
        boundary: "only the declared operator, boundary, sources, units, and time slab",
    },
    IntervalRole {
        id: "IR-05",
        supports: "finite-horizon tube ingredients",
        boundary: "not a validated reachability tube until the staged tube construction and its proof obligations exist",
    },
];

/// One worked thermal routing example.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ThermalExample {
    /// Stable example id.
    pub id: &'static str,
    /// Requested thermal claim.
    pub claim: &'static str,
    /// Required route.
    pub route: &'static str,
    /// Explicit refusal boundary.
    pub refusal: &'static str,
}

/// Worked steady-temperature and duty-cycle reliability examples.
pub const THERMAL_EXAMPLES: [ThermalExample; 2] = [
    ThermalExample {
        id: "THERMAL-STEADY-01",
        claim: "maximum steady-state component temperature over declared material, load, boundary, and discretization uncertainty",
        route: "finite interval/Taylor or residual evidence for the stated steady problem, with every model and discretization assumption retained",
        refusal: "does not establish reliability over a population of duty cycles",
    },
    ThermalExample {
        id: "THERMAL-DUTY-01",
        claim: "probability of thermal-limit exceedance over a declared duty-cycle population",
        route: "sequential or rare-event statistics plus sampling, dependence, model-form, calibration, and field-transfer evidence",
        refusal: "a steady temperature interval or one long simulated trajectory cannot mint this reliability claim",
    },
];

/// Fail-closed doctrine validation error.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CertificateRegimeError {
    /// The caller presented an unsupported schema version.
    SchemaVersion {
        /// Presented schema version.
        found: u16,
    },
    /// The closed v1 row count changed.
    RowCount {
        /// Presented row count.
        found: usize,
    },
    /// A row is not in canonical schema order.
    RowOrder {
        /// Zero-based row position.
        index: usize,
        /// Claim required at this position.
        expected: ClaimClass,
        /// Claim found at this position.
        found: ClaimClass,
    },
    /// A canonical row id changed.
    RowId {
        /// Zero-based row position.
        index: usize,
        /// Required canonical row id.
        expected: &'static str,
        /// Presented row id.
        found: &'static str,
    },
    /// The same claim appears more than once.
    DuplicateClaim {
        /// Repeated claim class.
        claim: ClaimClass,
    },
    /// A claim is paired with the wrong evidence regime.
    EvidenceMismatch {
        /// Claim whose route changed.
        claim: ClaimClass,
        /// Required evidence regime.
        expected: EvidenceRegime,
        /// Presented evidence regime.
        found: EvidenceRegime,
    },
    /// Required row text or capability metadata is empty.
    EmptyField {
        /// Canonical row id.
        row: &'static str,
        /// Empty field name.
        field: &'static str,
    },
    /// A row has no serving or policy capability.
    MissingCapability {
        /// Canonical row id.
        row: &'static str,
    },
    /// A capability names a non-workspace crate convention.
    InvalidCrateName {
        /// Canonical row id.
        row: &'static str,
        /// Presented crate name.
        crate_name: &'static str,
    },
    /// A staged capability falsely names implementation evidence.
    StagedCapabilityHasLocator {
        /// Canonical row id.
        row: &'static str,
        /// Capability that falsely named live source.
        capability: &'static str,
    },
    /// A non-staged capability has no source evidence.
    LiveCapabilityMissingLocator {
        /// Canonical row id.
        row: &'static str,
        /// Capability missing live source evidence.
        capability: &'static str,
    },
    /// A source locator is absolute or escapes the repository.
    InvalidSourceLocator {
        /// Canonical row id.
        row: &'static str,
        /// Rejected source locator.
        locator: &'static str,
    },
}

impl fmt::Display for CertificateRegimeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}

impl std::error::Error for CertificateRegimeError {}

const CANONICAL_ROW_IDS: [&str; 8] = [
    "CR-01", "CR-02", "CR-03", "CR-04", "CR-05", "CR-06", "CR-07", "CR-08",
];

/// Validate a candidate table against the closed v1 doctrine.
///
/// This is public so mutation and migration tooling can prove that a proposed
/// table is rejected before it reaches a future router.
pub fn validate_certificate_regimes(
    schema_version: u16,
    rows: &[CertificateRegimeRow],
) -> Result<(), CertificateRegimeError> {
    if schema_version != CERTIFICATE_REGIME_SCHEMA_VERSION {
        return Err(CertificateRegimeError::SchemaVersion {
            found: schema_version,
        });
    }
    if rows.len() != CERTIFICATE_REGIMES.len() {
        return Err(CertificateRegimeError::RowCount { found: rows.len() });
    }

    let mut seen_claims = [false; ClaimClass::ALL.len()];
    for (index, row) in rows.iter().enumerate() {
        if seen_claims[row.claim.canonical_index()] {
            return Err(CertificateRegimeError::DuplicateClaim { claim: row.claim });
        }
        seen_claims[row.claim.canonical_index()] = true;

        let expected_claim = ClaimClass::ALL[index];
        if row.claim != expected_claim {
            return Err(CertificateRegimeError::RowOrder {
                index,
                expected: expected_claim,
                found: row.claim,
            });
        }
        if row.id != CANONICAL_ROW_IDS[index] {
            return Err(CertificateRegimeError::RowId {
                index,
                expected: CANONICAL_ROW_IDS[index],
                found: row.id,
            });
        }

        let expected_evidence = row.claim.required_evidence();
        if row.evidence != expected_evidence {
            return Err(CertificateRegimeError::EvidenceMismatch {
                claim: row.claim,
                expected: expected_evidence,
                found: row.evidence,
            });
        }
        for (field, value) in [("scope", row.scope), ("no_claim", row.no_claim)] {
            if value.trim().is_empty() {
                return Err(CertificateRegimeError::EmptyField { row: row.id, field });
            }
        }
        if row.capabilities.is_empty() {
            return Err(CertificateRegimeError::MissingCapability { row: row.id });
        }
        for capability in row.capabilities {
            if capability.crate_name.trim().is_empty() {
                return Err(CertificateRegimeError::EmptyField {
                    row: row.id,
                    field: "capability.crate_name",
                });
            }
            if !capability.crate_name.starts_with("fs-") {
                return Err(CertificateRegimeError::InvalidCrateName {
                    row: row.id,
                    crate_name: capability.crate_name,
                });
            }
            if capability.capability.trim().is_empty() {
                return Err(CertificateRegimeError::EmptyField {
                    row: row.id,
                    field: "capability.capability",
                });
            }
            match (capability.status, capability.source_locator) {
                (CapabilityStatus::Staged, Some(_)) => {
                    return Err(CertificateRegimeError::StagedCapabilityHasLocator {
                        row: row.id,
                        capability: capability.capability,
                    });
                }
                (CapabilityStatus::Available | CapabilityStatus::Thin, None) => {
                    return Err(CertificateRegimeError::LiveCapabilityMissingLocator {
                        row: row.id,
                        capability: capability.capability,
                    });
                }
                (_, Some(locator))
                    if locator.starts_with('/')
                        || locator.split('/').any(|component| component == "..") =>
                {
                    return Err(CertificateRegimeError::InvalidSourceLocator {
                        row: row.id,
                        locator,
                    });
                }
                _ => {}
            }
        }
    }
    Ok(())
}

/// Validate and return the canonical doctrine.
pub fn certificate_regimes() -> Result<&'static [CertificateRegimeRow], CertificateRegimeError> {
    validate_certificate_regimes(CERTIFICATE_REGIME_SCHEMA_VERSION, &CERTIFICATE_REGIMES)?;
    Ok(&CERTIFICATE_REGIMES)
}

/// Total typed lookup in the canonical doctrine.
#[must_use]
pub fn certificate_regime(claim: ClaimClass) -> &'static CertificateRegimeRow {
    &CERTIFICATE_REGIMES[claim.canonical_index()]
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

/// Deterministic machine-readable doctrine table.
pub fn certificate_regime_json() -> Result<String, CertificateRegimeError> {
    let rows = certificate_regimes()?;
    let mut output = String::from("{\"schema_version\":");
    output.push_str(&CERTIFICATE_REGIME_SCHEMA_VERSION.to_string());
    output.push_str(",\"authority\":");
    push_json_string(&mut output, CERTIFICATE_REGIME_NO_CLAIM);
    output.push_str(",\"router_bead\":");
    push_json_string(&mut output, CERTIFICATE_REGIME_ROUTER_BEAD);
    output.push_str(",\"rows\":[");
    for (row_index, row) in rows.iter().enumerate() {
        if row_index != 0 {
            output.push(',');
        }
        output.push_str("{\"id\":");
        push_json_string(&mut output, row.id);
        output.push_str(",\"claim\":");
        push_json_string(&mut output, row.claim.code());
        output.push_str(",\"evidence\":");
        push_json_string(&mut output, row.evidence.code());
        output.push_str(",\"scope\":");
        push_json_string(&mut output, row.scope);
        output.push_str(",\"capabilities\":[");
        for (capability_index, capability) in row.capabilities.iter().enumerate() {
            if capability_index != 0 {
                output.push(',');
            }
            output.push_str("{\"crate\":");
            push_json_string(&mut output, capability.crate_name);
            output.push_str(",\"capability\":");
            push_json_string(&mut output, capability.capability);
            output.push_str(",\"status\":");
            push_json_string(&mut output, capability.status.code());
            output.push_str(",\"source_locator\":");
            if let Some(locator) = capability.source_locator {
                push_json_string(&mut output, locator);
            } else {
                output.push_str("null");
            }
            output.push('}');
        }
        output.push_str("],\"no_claim\":");
        push_json_string(&mut output, row.no_claim);
        output.push('}');
    }
    output.push_str("]}");
    Ok(output)
}

/// Code-derived Markdown table embedded verbatim in `docs/CERTIFICATE_REGIMES.md`.
pub fn certificate_regime_markdown_table() -> Result<String, CertificateRegimeError> {
    let rows = certificate_regimes()?;
    let mut output = String::from(
        "| ID | Claim class | Required evidence object | Current capability map | Applicability and no-claim boundary |\n\
         | --- | --- | --- | --- | --- |\n",
    );
    for row in rows {
        let capabilities = row
            .capabilities
            .iter()
            .map(|capability| {
                format!(
                    "`{}` / `{}` ({})",
                    capability.crate_name,
                    capability.capability,
                    capability.status.code()
                )
            })
            .collect::<Vec<_>>()
            .join("<br>");
        output.push_str(&format!(
            "| `{}` | {} | `{}` — {} | {} | Scope: {}<br>No claim: {} |\n",
            row.id,
            row.claim.title(),
            row.evidence.code(),
            row.evidence.title(),
            capabilities,
            row.scope,
            row.no_claim
        ));
    }
    Ok(output)
}
