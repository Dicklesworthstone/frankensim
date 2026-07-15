//! Typed verification-and-validation artifacts.
//!
//! This module owns the in-memory schemas and the structural admission rules.
//! Canonical transport lives in the sibling codec module; scientific authority
//! remains external and must be supplied by an authenticated package policy.

use core::fmt;
use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::ContentHash;

/// Version of the seven V&V artifact schemas.
pub const VV_SCHEMA_VERSION: u32 = 1;
/// Version of the structural rule matrix enforced by [`VvCase::admit`].
pub const VV_RULESET_VERSION: u32 = 1;
/// Stable family identity for canonical V&V payloads.
pub const VV_ARTIFACT_FAMILY: &str = "org.frankensim.fs-evidence.vv-artifact.v1";
/// Maximum UTF-8 bytes in a machine identity.
pub const MAX_VV_ID_BYTES: usize = 256;
/// Maximum UTF-8 bytes in one descriptive field.
pub const MAX_VV_TEXT_BYTES: usize = 64 * 1024;
/// Maximum rows accepted in any one collection.
pub const MAX_VV_ITEMS: usize = 4_096;
/// Maximum admitted dense/covariance dimension; bounds cubic validation work.
pub const MAX_VV_MATRIX_DIMENSION: usize = 128;

/// Stable, machine-readable structural rules.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum VvRule {
    SchemaIdentity,
    SchemaCardinality,
    SplitPartitionsDisjoint,
    SplitBlindHoldoutSealed,
    ColorCategoricalOnly,
    ValidationRequiresPhysicalReferent,
    QoiDependencyClosed,
    QoiDependencyIsolated,
    WaterfallModeDeclared,
    WaterfallArithmetic,
    WaterfallDependenceDeclared,
    ExperimentInstrumentCalibration,
    ExperimentClockSynchronization,
    ExperimentRepeatabilityCovariance,
    ExperimentDataAuthenticity,
    DiagnosticObservability,
    DiagnosticIdentifiability,
    DiagnosticConfounding,
    DiagnosticInverseCrime,
    ValidationMetricUncertainty,
    SolutionVerificationComplete,
    ApplicabilityDecision,
    ApplicabilityPolicy,
    ProcessConformanceSeparate,
    AssumptionRowComplete,
    AssumptionDomainEnforced,
    AssumptionA001,
    AssumptionA002,
    AssumptionA003,
    AssumptionA004,
    AssumptionA005,
    AssumptionA006,
    AssumptionA007,
    AssumptionA008,
    ReceiptBinding,
}

impl VvRule {
    /// Stable lowercase slug emitted in every refusal.
    #[must_use]
    pub const fn slug(self) -> &'static str {
        match self {
            Self::SchemaIdentity => "vv-schema-identity",
            Self::SchemaCardinality => "vv-schema-cardinality",
            Self::SplitPartitionsDisjoint => "vv-split-partitions-disjoint",
            Self::SplitBlindHoldoutSealed => "vv-split-blind-holdout-sealed",
            Self::ColorCategoricalOnly => "vv-color-categorical-only",
            Self::ValidationRequiresPhysicalReferent => "vv-validation-requires-physical-referent",
            Self::QoiDependencyClosed => "vv-qoi-dependency-closed",
            Self::QoiDependencyIsolated => "vv-qoi-dependency-isolated",
            Self::WaterfallModeDeclared => "vv-waterfall-mode-declared",
            Self::WaterfallArithmetic => "vv-waterfall-arithmetic",
            Self::WaterfallDependenceDeclared => "vv-waterfall-dependence-declared",
            Self::ExperimentInstrumentCalibration => "vv-experiment-instrument-calibration",
            Self::ExperimentClockSynchronization => "vv-experiment-clock-synchronization",
            Self::ExperimentRepeatabilityCovariance => "vv-experiment-repeatability-covariance",
            Self::ExperimentDataAuthenticity => "vv-experiment-data-authenticity",
            Self::DiagnosticObservability => "vv-diagnostic-observability",
            Self::DiagnosticIdentifiability => "vv-diagnostic-identifiability",
            Self::DiagnosticConfounding => "vv-diagnostic-confounding",
            Self::DiagnosticInverseCrime => "vv-diagnostic-inverse-crime",
            Self::ValidationMetricUncertainty => "vv-validation-metric-uncertainty",
            Self::SolutionVerificationComplete => "vv-solution-verification-complete",
            Self::ApplicabilityDecision => "vv-applicability-decision",
            Self::ApplicabilityPolicy => "vv-applicability-policy",
            Self::ProcessConformanceSeparate => "vv-process-conformance-separate",
            Self::AssumptionRowComplete => "vv-assumption-row-complete",
            Self::AssumptionDomainEnforced => "vv-assumption-domain-enforced",
            Self::AssumptionA001 => "vv-assumption-a001",
            Self::AssumptionA002 => "vv-assumption-a002",
            Self::AssumptionA003 => "vv-assumption-a003",
            Self::AssumptionA004 => "vv-assumption-a004",
            Self::AssumptionA005 => "vv-assumption-a005",
            Self::AssumptionA006 => "vv-assumption-a006",
            Self::AssumptionA007 => "vv-assumption-a007",
            Self::AssumptionA008 => "vv-assumption-a008",
            Self::ReceiptBinding => "vv-receipt-binding",
        }
    }
}

/// One deterministic, actionable schema refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VvViolation {
    rule: VvRule,
    artifact_id: Option<String>,
    qoi_id: Option<String>,
    field: &'static str,
    detail: String,
}

impl VvViolation {
    pub(crate) fn new(
        rule: VvRule,
        artifact_id: Option<String>,
        qoi_id: Option<String>,
        field: &'static str,
        detail: impl Into<String>,
    ) -> Self {
        Self {
            rule,
            artifact_id,
            qoi_id,
            field,
            detail: detail.into(),
        }
    }

    #[must_use]
    pub const fn rule(&self) -> VvRule {
        self.rule
    }

    #[must_use]
    pub const fn rule_slug(&self) -> &'static str {
        self.rule.slug()
    }

    #[must_use]
    pub fn artifact_id(&self) -> Option<&str> {
        self.artifact_id.as_deref()
    }

    #[must_use]
    pub fn qoi_id(&self) -> Option<&str> {
        self.qoi_id.as_deref()
    }

    #[must_use]
    pub const fn field(&self) -> &'static str {
        self.field
    }

    #[must_use]
    pub fn detail(&self) -> &str {
        &self.detail
    }
}

impl fmt::Display for VvViolation {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.rule.slug(), self.detail)?;
        if let Some(artifact) = &self.artifact_id {
            write!(f, " [artifact={artifact}]")?;
        }
        if let Some(qoi) = &self.qoi_id {
            write!(f, " [qoi={qoi}]")?;
        }
        write!(f, " [field={}]", self.field)
    }
}

impl std::error::Error for VvViolation {}

/// Deterministically ordered collection of schema refusals.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VvErrors {
    violations: Vec<VvViolation>,
}

impl VvErrors {
    pub(crate) fn one(violation: VvViolation) -> Self {
        Self {
            violations: vec![violation],
        }
    }

    pub(crate) fn from_vec(mut violations: Vec<VvViolation>) -> Self {
        violations.sort_by(|a, b| {
            (
                a.artifact_id.as_deref(),
                a.qoi_id.as_deref(),
                a.rule,
                a.field,
                a.detail.as_str(),
            )
                .cmp(&(
                    b.artifact_id.as_deref(),
                    b.qoi_id.as_deref(),
                    b.rule,
                    b.field,
                    b.detail.as_str(),
                ))
        });
        violations.dedup();
        Self { violations }
    }

    #[must_use]
    pub fn violations(&self) -> &[VvViolation] {
        &self.violations
    }
}

impl fmt::Display for VvErrors {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{} V&V rule violation(s)", self.violations.len())?;
        for violation in &self.violations {
            write!(f, "; {violation}")?;
        }
        Ok(())
    }
}

impl std::error::Error for VvErrors {}

fn invalid(
    rule: VvRule,
    artifact: Option<&str>,
    qoi: Option<&str>,
    field: &'static str,
    detail: impl Into<String>,
) -> VvErrors {
    VvErrors::one(VvViolation::new(
        rule,
        artifact.map(str::to_owned),
        qoi.map(str::to_owned),
        field,
        detail,
    ))
}

fn validate_id(value: &str, field: &'static str) -> Result<(), VvErrors> {
    let valid_byte = |byte: u8| {
        byte.is_ascii_alphanumeric()
            || matches!(
                byte,
                b'-' | b'_' | b'.' | b'/' | b':' | b'@' | b'+' | b'=' | b'^' | b'*'
            )
    };
    let placeholder = [
        "-",
        "?",
        "n/a",
        "na",
        "none",
        "not run",
        "pending",
        "placeholder",
        "tbd",
        "todo",
        "unknown",
    ]
    .iter()
    .any(|candidate| value.eq_ignore_ascii_case(candidate));
    if value.is_empty()
        || value.len() > MAX_VV_ID_BYTES
        || value.trim() != value
        || placeholder
        || !value.bytes().all(valid_byte)
    {
        return Err(invalid(
            VvRule::SchemaIdentity,
            None,
            None,
            field,
            format!("invalid bounded machine identity {value:?}"),
        ));
    }
    Ok(())
}

fn validate_text(value: &str, field: &'static str) -> Result<(), VvErrors> {
    if value.trim().is_empty()
        || value.len() > MAX_VV_TEXT_BYTES
        || value.chars().any(|ch| ch == '\0')
    {
        return Err(invalid(
            VvRule::SchemaIdentity,
            None,
            None,
            field,
            "text must be non-blank, NUL-free, and within the transport bound",
        ));
    }
    Ok(())
}

macro_rules! vv_id {
    ($name:ident, $field:literal) => {
        #[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(String);

        impl $name {
            pub fn try_new(value: impl Into<String>) -> Result<Self, VvErrors> {
                let value = value.into();
                validate_id(&value, $field)?;
                Ok(Self(value))
            }

            #[must_use]
            pub fn as_str(&self) -> &str {
                &self.0
            }

            pub(crate) fn from_canonical(value: String) -> Result<Self, VvErrors> {
                Self::try_new(value)
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str(&self.0)
            }
        }
    };
}

vv_id!(ArtifactId, "artifact_id");
vv_id!(QoiId, "qoi_id");
vv_id!(ObservationId, "observation_id");
vv_id!(AssumptionId, "assumption_id");
vv_id!(AxisId, "axis_id");
vv_id!(UnitId, "unit_id");

/// The seven top-level artifact families.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ArtifactKind {
    ContextOfUse,
    ValidationPlan,
    ExperimentArtifact,
    CalibrationSplit,
    SolutionVerificationReceipt,
    PredictionAssessment,
    AssumptionsLedger,
}

impl ArtifactKind {
    #[must_use]
    pub const fn slug(self) -> &'static str {
        match self {
            Self::ContextOfUse => "context-of-use",
            Self::ValidationPlan => "validation-plan",
            Self::ExperimentArtifact => "experiment-artifact",
            Self::CalibrationSplit => "calibration-split",
            Self::SolutionVerificationReceipt => "solution-verification-receipt",
            Self::PredictionAssessment => "prediction-assessment",
            Self::AssumptionsLedger => "assumptions-ledger",
        }
    }
}

/// Content-bound reference to one V&V artifact.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ArtifactRef {
    kind: ArtifactKind,
    id: ArtifactId,
    hash: ContentHash,
}

impl ArtifactRef {
    #[must_use]
    pub fn new(kind: ArtifactKind, id: ArtifactId, hash: ContentHash) -> Self {
        Self { kind, id, hash }
    }

    #[must_use]
    pub const fn kind(&self) -> ArtifactKind {
        self.kind
    }

    #[must_use]
    pub const fn id(&self) -> &ArtifactId {
        &self.id
    }

    #[must_use]
    pub const fn hash(&self) -> ContentHash {
        self.hash
    }
}

/// Explicit random-seed declaration required by every top-level artifact.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SeedDeclaration {
    Fixed(u64),
    NotApplicable { reason: String },
}

/// Explicit bounded budget or an equally explicit not-applicable reason.
#[derive(Debug, Clone, PartialEq)]
pub enum DeclaredBudget<T> {
    Limit(T),
    NotApplicable { reason: String },
}

/// Agent-native Five-Explicits header shared by all seven artifacts.
#[derive(Debug, Clone, PartialEq)]
pub struct ArtifactHeader {
    id: ArtifactId,
    units: Vec<UnitId>,
    seed: SeedDeclaration,
    accuracy: DeclaredBudget<f64>,
    time_ms: DeclaredBudget<u64>,
    memory_bytes: DeclaredBudget<u64>,
    versions: BTreeMap<String, String>,
    capabilities: BTreeSet<String>,
}

impl ArtifactHeader {
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        id: ArtifactId,
        units: Vec<UnitId>,
        seed: SeedDeclaration,
        accuracy: DeclaredBudget<f64>,
        time_ms: DeclaredBudget<u64>,
        memory_bytes: DeclaredBudget<u64>,
        versions: Vec<(String, String)>,
        capabilities: Vec<String>,
    ) -> Result<Self, VvErrors> {
        if units.is_empty() || units.len() > MAX_VV_ITEMS {
            return Err(invalid(
                VvRule::SchemaIdentity,
                Some(id.as_str()),
                None,
                "header.units",
                "units must be explicit and bounded",
            ));
        }
        let mut units = units;
        units.sort();
        units.dedup();
        if let SeedDeclaration::NotApplicable { reason } = &seed {
            validate_text(reason, "header.seed.not_applicable")?;
        }
        match &accuracy {
            DeclaredBudget::Limit(value) if value.is_finite() && *value >= 0.0 => {}
            DeclaredBudget::NotApplicable { reason } => {
                validate_text(reason, "header.accuracy.not_applicable")?;
            }
            DeclaredBudget::Limit(_) => {
                return Err(invalid(
                    VvRule::SchemaIdentity,
                    Some(id.as_str()),
                    None,
                    "header.accuracy",
                    "accuracy budget must be finite and non-negative",
                ));
            }
        }
        for (field, budget) in [
            ("header.time_ms", &time_ms),
            ("header.memory_bytes", &memory_bytes),
        ] {
            match budget {
                DeclaredBudget::Limit(value) if *value > 0 => {}
                DeclaredBudget::NotApplicable { reason } => validate_text(reason, field)?,
                DeclaredBudget::Limit(_) => {
                    return Err(invalid(
                        VvRule::SchemaIdentity,
                        Some(id.as_str()),
                        None,
                        field,
                        "resource budgets must be positive or explicitly not applicable",
                    ));
                }
            }
        }
        if versions.is_empty() || versions.len() > MAX_VV_ITEMS {
            return Err(invalid(
                VvRule::SchemaIdentity,
                Some(id.as_str()),
                None,
                "header.versions",
                "at least one exact version is required",
            ));
        }
        let mut version_map = BTreeMap::new();
        for (component, version) in versions {
            validate_id(&component, "header.version.component")?;
            validate_id(&version, "header.version.value")?;
            if version_map.insert(component, version).is_some() {
                return Err(invalid(
                    VvRule::SchemaIdentity,
                    Some(id.as_str()),
                    None,
                    "header.versions",
                    "version component identities must be unique",
                ));
            }
        }
        if capabilities.is_empty() || capabilities.len() > MAX_VV_ITEMS {
            return Err(invalid(
                VvRule::SchemaIdentity,
                Some(id.as_str()),
                None,
                "header.capabilities",
                "capabilities must be explicit and bounded",
            ));
        }
        let mut capability_set = BTreeSet::new();
        for capability in capabilities {
            validate_id(&capability, "header.capability")?;
            if !capability_set.insert(capability) {
                return Err(invalid(
                    VvRule::SchemaIdentity,
                    Some(id.as_str()),
                    None,
                    "header.capabilities",
                    "capability identities must be unique",
                ));
            }
        }
        Ok(Self {
            id,
            units,
            seed,
            accuracy,
            time_ms,
            memory_bytes,
            versions: version_map,
            capabilities: capability_set,
        })
    }

    #[must_use]
    pub const fn id(&self) -> &ArtifactId {
        &self.id
    }

    #[must_use]
    pub fn units(&self) -> &[UnitId] {
        &self.units
    }

    #[must_use]
    pub const fn seed(&self) -> &SeedDeclaration {
        &self.seed
    }

    #[must_use]
    pub const fn accuracy(&self) -> &DeclaredBudget<f64> {
        &self.accuracy
    }

    #[must_use]
    pub const fn time_ms(&self) -> &DeclaredBudget<u64> {
        &self.time_ms
    }

    #[must_use]
    pub const fn memory_bytes(&self) -> &DeclaredBudget<u64> {
        &self.memory_bytes
    }

    #[must_use]
    pub const fn versions(&self) -> &BTreeMap<String, String> {
        &self.versions
    }

    #[must_use]
    pub const fn capabilities(&self) -> &BTreeSet<String> {
        &self.capabilities
    }
}

/// One strict continuous applicability constraint.
#[derive(Debug, Clone, PartialEq)]
pub struct NumericDomainAxis {
    axis: AxisId,
    unit: UnitId,
    lo: f64,
    hi: f64,
}

impl NumericDomainAxis {
    pub fn try_new(axis: AxisId, unit: UnitId, lo: f64, hi: f64) -> Result<Self, VvErrors> {
        if !lo.is_finite() || !hi.is_finite() || lo > hi {
            return Err(invalid(
                VvRule::AssumptionDomainEnforced,
                None,
                None,
                "numeric_domain",
                "numeric applicability bounds must be finite and ordered",
            ));
        }
        Ok(Self { axis, unit, lo, hi })
    }

    #[must_use]
    pub const fn axis(&self) -> &AxisId {
        &self.axis
    }

    #[must_use]
    pub const fn unit(&self) -> &UnitId {
        &self.unit
    }

    #[must_use]
    pub const fn bounds(&self) -> (f64, f64) {
        (self.lo, self.hi)
    }
}

/// One strict categorical applicability constraint.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CategoricalDomainAxis {
    axis: AxisId,
    allowed: BTreeSet<String>,
}

impl CategoricalDomainAxis {
    pub fn try_new(axis: AxisId, allowed: Vec<String>) -> Result<Self, VvErrors> {
        if allowed.is_empty() || allowed.len() > MAX_VV_ITEMS {
            return Err(invalid(
                VvRule::AssumptionDomainEnforced,
                None,
                None,
                "categorical_domain",
                "a categorical axis needs a bounded non-empty allowed set",
            ));
        }
        let mut canonical = BTreeSet::new();
        for value in allowed {
            validate_id(&value, "categorical_value")?;
            if !canonical.insert(value) {
                return Err(invalid(
                    VvRule::AssumptionDomainEnforced,
                    None,
                    None,
                    "categorical_domain",
                    "categorical values must be unique",
                ));
            }
        }
        Ok(Self {
            axis,
            allowed: canonical,
        })
    }

    #[must_use]
    pub const fn axis(&self) -> &AxisId {
        &self.axis
    }

    #[must_use]
    pub const fn allowed(&self) -> &BTreeSet<String> {
        &self.allowed
    }
}

/// Mixed numeric/categorical domain used by context and assumption rows.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ApplicabilityDomain {
    numeric: BTreeMap<AxisId, NumericDomainAxis>,
    categorical: BTreeMap<AxisId, CategoricalDomainAxis>,
}

impl ApplicabilityDomain {
    #[must_use]
    pub fn unconstrained() -> Self {
        Self::default()
    }

    pub fn try_new(
        numeric: Vec<NumericDomainAxis>,
        categorical: Vec<CategoricalDomainAxis>,
    ) -> Result<Self, VvErrors> {
        if numeric.len().saturating_add(categorical.len()) > MAX_VV_ITEMS {
            return Err(invalid(
                VvRule::SchemaCardinality,
                None,
                None,
                "applicability_domain",
                "too many applicability axes",
            ));
        }
        let mut numeric_map = BTreeMap::new();
        let mut categorical_map = BTreeMap::new();
        for row in numeric {
            if categorical_map.contains_key(row.axis())
                || numeric_map.insert(row.axis.clone(), row).is_some()
            {
                return Err(invalid(
                    VvRule::AssumptionDomainEnforced,
                    None,
                    None,
                    "applicability_domain",
                    "an applicability axis may be declared exactly once",
                ));
            }
        }
        for row in categorical {
            if numeric_map.contains_key(row.axis())
                || categorical_map.insert(row.axis.clone(), row).is_some()
            {
                return Err(invalid(
                    VvRule::AssumptionDomainEnforced,
                    None,
                    None,
                    "applicability_domain",
                    "an applicability axis may be declared exactly once",
                ));
            }
        }
        Ok(Self {
            numeric: numeric_map,
            categorical: categorical_map,
        })
    }

    #[must_use]
    pub const fn numeric(&self) -> &BTreeMap<AxisId, NumericDomainAxis> {
        &self.numeric
    }

    #[must_use]
    pub const fn categorical(&self) -> &BTreeMap<AxisId, CategoricalDomainAxis> {
        &self.categorical
    }

    fn violations(&self, point: &ApplicabilityPoint) -> Vec<DomainViolation> {
        let mut violations = Vec::new();
        for (axis, range) in &self.numeric {
            match point.numeric.get(axis) {
                Some(value) if value.is_finite() && *value >= range.lo && *value <= range.hi => {}
                Some(value) => violations.push(DomainViolation::Numeric {
                    axis: axis.clone(),
                    value: *value,
                    lo: range.lo,
                    hi: range.hi,
                }),
                None => violations.push(DomainViolation::Missing { axis: axis.clone() }),
            }
        }
        for (axis, range) in &self.categorical {
            match point.categorical.get(axis) {
                Some(value) if range.allowed.contains(value) => {}
                Some(value) => violations.push(DomainViolation::Categorical {
                    axis: axis.clone(),
                    value: value.clone(),
                }),
                None => violations.push(DomainViolation::Missing { axis: axis.clone() }),
            }
        }
        violations
    }
}

/// Runtime point checked against an applicability domain.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct ApplicabilityPoint {
    numeric: BTreeMap<AxisId, f64>,
    categorical: BTreeMap<AxisId, String>,
}

impl ApplicabilityPoint {
    pub fn try_new(
        numeric: Vec<(AxisId, f64)>,
        categorical: Vec<(AxisId, String)>,
    ) -> Result<Self, VvErrors> {
        if numeric.len().saturating_add(categorical.len()) > MAX_VV_ITEMS {
            return Err(invalid(
                VvRule::SchemaCardinality,
                None,
                None,
                "applicability_point",
                "too many applicability values",
            ));
        }
        let mut numeric_map = BTreeMap::new();
        let mut categorical_map = BTreeMap::new();
        for (axis, value) in numeric {
            if !value.is_finite()
                || categorical_map.contains_key(&axis)
                || numeric_map.insert(axis, value).is_some()
            {
                return Err(invalid(
                    VvRule::ApplicabilityDecision,
                    None,
                    None,
                    "applicability_point",
                    "runtime numeric axes must be finite and unique",
                ));
            }
        }
        for (axis, value) in categorical {
            validate_id(&value, "categorical_value")?;
            if numeric_map.contains_key(&axis) || categorical_map.insert(axis, value).is_some() {
                return Err(invalid(
                    VvRule::ApplicabilityDecision,
                    None,
                    None,
                    "applicability_point",
                    "runtime categorical axes must be unique",
                ));
            }
        }
        Ok(Self {
            numeric: numeric_map,
            categorical: categorical_map,
        })
    }

    #[must_use]
    pub const fn numeric(&self) -> &BTreeMap<AxisId, f64> {
        &self.numeric
    }

    #[must_use]
    pub const fn categorical(&self) -> &BTreeMap<AxisId, String> {
        &self.categorical
    }
}

/// A concrete reason a point is outside its declared domain.
#[derive(Debug, Clone, PartialEq)]
pub enum DomainViolation {
    Missing {
        axis: AxisId,
    },
    Numeric {
        axis: AxisId,
        value: f64,
        lo: f64,
        hi: f64,
    },
    Categorical {
        axis: AxisId,
        value: String,
    },
    Assumption {
        id: AssumptionId,
    },
}

impl DomainViolation {
    fn has_same_canonical_bits(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Missing { axis: left }, Self::Missing { axis: right }) => left == right,
            (
                Self::Numeric {
                    axis: left_axis,
                    value: left_value,
                    lo: left_lo,
                    hi: left_hi,
                },
                Self::Numeric {
                    axis: right_axis,
                    value: right_value,
                    lo: right_lo,
                    hi: right_hi,
                },
            ) => {
                left_axis == right_axis
                    && left_value.to_bits() == right_value.to_bits()
                    && left_lo.to_bits() == right_lo.to_bits()
                    && left_hi.to_bits() == right_hi.to_bits()
            }
            (
                Self::Categorical {
                    axis: left_axis,
                    value: left_value,
                },
                Self::Categorical {
                    axis: right_axis,
                    value: right_value,
                },
            ) => left_axis == right_axis && left_value == right_value,
            (Self::Assumption { id: left }, Self::Assumption { id: right }) => left == right,
            _ => false,
        }
    }
}

/// Required treatment when applicability fails.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ApplicabilityPolicy {
    Demote,
    Refuse,
}

/// Derived applicability result; there is no silent extrapolation variant.
#[derive(Debug, Clone, PartialEq)]
pub enum ApplicabilityDecision {
    InDomain,
    Demoted { violations: Vec<DomainViolation> },
    Refused { violations: Vec<DomainViolation> },
}

impl ApplicabilityDecision {
    fn derive(policy: ApplicabilityPolicy, violations: Vec<DomainViolation>) -> Self {
        if violations.is_empty() {
            Self::InDomain
        } else {
            match policy {
                ApplicabilityPolicy::Demote => Self::Demoted { violations },
                ApplicabilityPolicy::Refuse => Self::Refused { violations },
            }
        }
    }

    fn has_same_canonical_bits(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::InDomain, Self::InDomain) => true,
            (Self::Demoted { violations: left }, Self::Demoted { violations: right })
            | (Self::Refused { violations: left }, Self::Refused { violations: right }) => {
                left.len() == right.len()
                    && left
                        .iter()
                        .zip(right)
                        .all(|(left, right)| left.has_same_canonical_bits(right))
            }
            _ => false,
        }
    }
}

/// Acceptance predicate for one QoI.
#[derive(Debug, Clone, PartialEq)]
pub enum AcceptanceCriterion {
    ClosedRange { lo: f64, hi: f64 },
    AbsoluteErrorAtMost { limit: f64 },
    RelativeErrorAtMost { limit: f64 },
    CategoryEquals { expected: String },
}

impl AcceptanceCriterion {
    fn validate(&self) -> Result<(), VvErrors> {
        match self {
            Self::ClosedRange { lo, hi } if lo.is_finite() && hi.is_finite() && lo <= hi => Ok(()),
            Self::AbsoluteErrorAtMost { limit } | Self::RelativeErrorAtMost { limit }
                if limit.is_finite() && *limit >= 0.0 =>
            {
                Ok(())
            }
            Self::CategoryEquals { expected } => validate_id(expected, "acceptance.category"),
            _ => Err(invalid(
                VvRule::SchemaIdentity,
                None,
                None,
                "acceptance",
                "acceptance bounds must be finite, ordered, and non-negative",
            )),
        }
    }
}

/// One decision-relevant quantity of interest.
#[derive(Debug, Clone, PartialEq)]
pub struct QoiSpec {
    id: QoiId,
    name: String,
    unit: UnitId,
    acceptance: AcceptanceCriterion,
}

impl QoiSpec {
    pub fn try_new(
        id: QoiId,
        name: impl Into<String>,
        unit: UnitId,
        acceptance: AcceptanceCriterion,
    ) -> Result<Self, VvErrors> {
        let name = name.into();
        validate_text(&name, "qoi.name")?;
        acceptance.validate()?;
        Ok(Self {
            id,
            name,
            unit,
            acceptance,
        })
    }

    #[must_use]
    pub const fn id(&self) -> &QoiId {
        &self.id
    }

    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    #[must_use]
    pub const fn unit(&self) -> &UnitId {
        &self.unit
    }

    #[must_use]
    pub const fn acceptance(&self) -> &AcceptanceCriterion {
        &self.acceptance
    }
}

/// What decision the model serves, its QoIs, criteria, and applicability.
#[derive(Debug, Clone, PartialEq)]
pub struct ContextOfUse {
    header: ArtifactHeader,
    decision: String,
    qois: BTreeMap<QoiId, QoiSpec>,
    applicability: ApplicabilityDomain,
    applicability_policy: ApplicabilityPolicy,
}

impl ContextOfUse {
    pub fn try_new(
        header: ArtifactHeader,
        decision: impl Into<String>,
        qois: Vec<QoiSpec>,
        applicability: ApplicabilityDomain,
        applicability_policy: ApplicabilityPolicy,
    ) -> Result<Self, VvErrors> {
        let decision = decision.into();
        validate_text(&decision, "context.decision")?;
        if qois.is_empty() || qois.len() > MAX_VV_ITEMS {
            return Err(invalid(
                VvRule::QoiDependencyClosed,
                Some(header.id().as_str()),
                None,
                "context.qois",
                "a context needs a bounded non-empty QoI set",
            ));
        }
        let mut by_id = BTreeMap::new();
        for qoi in qois {
            let qoi_id = qoi.id.clone();
            if by_id.insert(qoi_id.clone(), qoi).is_some() {
                return Err(invalid(
                    VvRule::QoiDependencyClosed,
                    Some(header.id().as_str()),
                    Some(qoi_id.as_str()),
                    "context.qois",
                    "QoI identities must be unique",
                ));
            }
        }
        Ok(Self {
            header,
            decision,
            qois: by_id,
            applicability,
            applicability_policy,
        })
    }

    #[must_use]
    pub const fn id(&self) -> &ArtifactId {
        self.header.id()
    }

    #[must_use]
    pub const fn header(&self) -> &ArtifactHeader {
        &self.header
    }

    #[must_use]
    pub fn decision(&self) -> &str {
        &self.decision
    }

    #[must_use]
    pub const fn qois(&self) -> &BTreeMap<QoiId, QoiSpec> {
        &self.qois
    }

    #[must_use]
    pub const fn applicability(&self) -> &ApplicabilityDomain {
        &self.applicability
    }

    #[must_use]
    pub const fn applicability_policy(&self) -> ApplicabilityPolicy {
        self.applicability_policy
    }
}

/// One retained diagnostic decision.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DiagnosticRecord {
    passed: bool,
    artifact_hash: ContentHash,
    detail: String,
}

impl DiagnosticRecord {
    pub fn try_new(
        passed: bool,
        artifact_hash: ContentHash,
        detail: impl Into<String>,
    ) -> Result<Self, VvErrors> {
        let detail = detail.into();
        validate_text(&detail, "diagnostic.detail")?;
        Ok(Self {
            passed,
            artifact_hash,
            detail,
        })
    }

    #[must_use]
    pub const fn passed(&self) -> bool {
        self.passed
    }

    #[must_use]
    pub const fn artifact_hash(&self) -> ContentHash {
        self.artifact_hash
    }

    #[must_use]
    pub fn detail(&self) -> &str {
        &self.detail
    }
}

/// Required observability, identifiability, confounding, and inverse-crime checks.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DiagnosticPlan {
    observability: DiagnosticRecord,
    identifiability: DiagnosticRecord,
    confounding: DiagnosticRecord,
    inverse_crime: DiagnosticRecord,
}

impl DiagnosticPlan {
    #[must_use]
    pub fn new(
        observability: DiagnosticRecord,
        identifiability: DiagnosticRecord,
        confounding: DiagnosticRecord,
        inverse_crime: DiagnosticRecord,
    ) -> Self {
        Self {
            observability,
            identifiability,
            confounding,
            inverse_crime,
        }
    }

    #[must_use]
    pub const fn observability(&self) -> &DiagnosticRecord {
        &self.observability
    }

    #[must_use]
    pub const fn identifiability(&self) -> &DiagnosticRecord {
        &self.identifiability
    }

    #[must_use]
    pub const fn confounding(&self) -> &DiagnosticRecord {
        &self.confounding
    }

    #[must_use]
    pub const fn inverse_crime(&self) -> &DiagnosticRecord {
        &self.inverse_crime
    }

    fn violations(&self, artifact: &ArtifactId, qoi: &QoiId) -> Vec<VvViolation> {
        let checks = [
            (
                VvRule::DiagnosticObservability,
                "diagnostics.observability",
                &self.observability,
            ),
            (
                VvRule::DiagnosticIdentifiability,
                "diagnostics.identifiability",
                &self.identifiability,
            ),
            (
                VvRule::DiagnosticConfounding,
                "diagnostics.confounding",
                &self.confounding,
            ),
            (
                VvRule::DiagnosticInverseCrime,
                "diagnostics.inverse_crime",
                &self.inverse_crime,
            ),
        ];
        checks
            .into_iter()
            .filter(|(_, _, record)| !record.passed)
            .map(|(rule, field, record)| {
                VvViolation::new(
                    rule,
                    Some(artifact.as_str().to_owned()),
                    Some(qoi.as_str().to_owned()),
                    field,
                    record.detail.clone(),
                )
            })
            .collect()
    }
}

/// Metric requested by a validation plan.
#[derive(Debug, Clone, PartialEq)]
pub enum ValidationMetricSpec {
    IntervalAgreement,
    NormalizedDiscrepancy { maximum: f64 },
    PosteriorPredictive { minimum_tail_probability: f64 },
}

impl ValidationMetricSpec {
    pub(super) fn canonical_key(&self) -> (u8, u64) {
        match self {
            Self::IntervalAgreement => (0, 0),
            Self::NormalizedDiscrepancy { maximum } => (1, canonical_float_key(*maximum)),
            Self::PosteriorPredictive {
                minimum_tail_probability,
            } => (2, canonical_float_key(*minimum_tail_probability)),
        }
    }

    fn validate(&self) -> Result<(), VvErrors> {
        match self {
            Self::IntervalAgreement => Ok(()),
            Self::NormalizedDiscrepancy { maximum } if maximum.is_finite() && *maximum >= 0.0 => {
                Ok(())
            }
            Self::PosteriorPredictive {
                minimum_tail_probability,
            } if minimum_tail_probability.is_finite()
                && *minimum_tail_probability > 0.0
                && *minimum_tail_probability < 1.0 =>
            {
                Ok(())
            }
            _ => Err(invalid(
                VvRule::ValidationMetricUncertainty,
                None,
                None,
                "validation_metric_spec",
                "metric thresholds must be finite and inside their declared domain",
            )),
        }
    }
}

fn canonical_float_key(value: f64) -> u64 {
    if value == 0.0 { 0 } else { value.to_bits() }
}

/// Per-QoI validation plan row.
#[derive(Debug, Clone, PartialEq)]
pub struct QoiValidationPlan {
    qoi: QoiId,
    experiments: Vec<ArtifactRef>,
    split: ArtifactRef,
    metrics: Vec<ValidationMetricSpec>,
    diagnostics: DiagnosticPlan,
}

impl QoiValidationPlan {
    pub fn try_new(
        qoi: QoiId,
        experiments: Vec<ArtifactRef>,
        split: ArtifactRef,
        metrics: Vec<ValidationMetricSpec>,
        diagnostics: DiagnosticPlan,
    ) -> Result<Self, VvErrors> {
        if experiments.is_empty() || experiments.len() > MAX_VV_ITEMS {
            return Err(invalid(
                VvRule::ValidationRequiresPhysicalReferent,
                None,
                Some(qoi.as_str()),
                "validation_plan.experiments",
                "each QoI needs a bounded non-empty experiment set",
            ));
        }
        if experiments
            .iter()
            .any(|reference| reference.kind != ArtifactKind::ExperimentArtifact)
        {
            return Err(invalid(
                VvRule::ValidationRequiresPhysicalReferent,
                None,
                Some(qoi.as_str()),
                "validation_plan.experiments",
                "validation experiment references must name ExperimentArtifact values",
            ));
        }
        if split.kind != ArtifactKind::CalibrationSplit {
            return Err(invalid(
                VvRule::SplitPartitionsDisjoint,
                None,
                Some(qoi.as_str()),
                "validation_plan.split",
                "the partition reference must name a CalibrationSplit",
            ));
        }
        if metrics.is_empty() || metrics.len() > MAX_VV_ITEMS {
            return Err(invalid(
                VvRule::ValidationMetricUncertainty,
                None,
                Some(qoi.as_str()),
                "validation_plan.metrics",
                "each QoI needs a bounded non-empty metric set",
            ));
        }
        for metric in &metrics {
            metric.validate()?;
        }
        let metric_count = metrics.len();
        let mut metrics = metrics;
        metrics.sort_by_key(ValidationMetricSpec::canonical_key);
        metrics.dedup();
        if metrics.len() != metric_count {
            return Err(invalid(
                VvRule::ValidationMetricUncertainty,
                None,
                Some(qoi.as_str()),
                "validation_plan.metrics",
                "validation metric specifications must be unique",
            ));
        }
        let experiment_count = experiments.len();
        let mut experiments = experiments;
        experiments.sort();
        experiments.dedup();
        if experiments.len() != experiment_count {
            return Err(invalid(
                VvRule::SchemaIdentity,
                None,
                Some(qoi.as_str()),
                "validation_plan.experiments",
                "experiment references must be unique",
            ));
        }
        Ok(Self {
            qoi,
            experiments,
            split,
            metrics,
            diagnostics,
        })
    }

    #[must_use]
    pub const fn qoi(&self) -> &QoiId {
        &self.qoi
    }

    #[must_use]
    pub fn experiments(&self) -> &[ArtifactRef] {
        &self.experiments
    }

    #[must_use]
    pub const fn split(&self) -> &ArtifactRef {
        &self.split
    }

    #[must_use]
    pub fn metrics(&self) -> &[ValidationMetricSpec] {
        &self.metrics
    }

    #[must_use]
    pub const fn diagnostics(&self) -> &DiagnosticPlan {
        &self.diagnostics
    }
}

/// Which experiments and diagnostics validate each Context-of-Use QoI.
#[derive(Debug, Clone, PartialEq)]
pub struct ValidationPlan {
    header: ArtifactHeader,
    context: ArtifactRef,
    by_qoi: BTreeMap<QoiId, QoiValidationPlan>,
}

impl ValidationPlan {
    pub fn try_new(
        header: ArtifactHeader,
        context: ArtifactRef,
        rows: Vec<QoiValidationPlan>,
    ) -> Result<Self, VvErrors> {
        if context.kind != ArtifactKind::ContextOfUse {
            return Err(invalid(
                VvRule::QoiDependencyClosed,
                Some(header.id().as_str()),
                None,
                "validation_plan.context",
                "validation plan context must reference ContextOfUse",
            ));
        }
        if rows.is_empty() || rows.len() > MAX_VV_ITEMS {
            return Err(invalid(
                VvRule::QoiDependencyClosed,
                Some(header.id().as_str()),
                None,
                "validation_plan.by_qoi",
                "validation plan needs a bounded non-empty QoI map",
            ));
        }
        let mut by_qoi = BTreeMap::new();
        for row in rows {
            let qoi = row.qoi.clone();
            if by_qoi.insert(qoi.clone(), row).is_some() {
                return Err(invalid(
                    VvRule::QoiDependencyClosed,
                    Some(header.id().as_str()),
                    Some(qoi.as_str()),
                    "validation_plan.by_qoi",
                    "a QoI may have exactly one validation-plan row",
                ));
            }
        }
        Ok(Self {
            header,
            context,
            by_qoi,
        })
    }

    #[must_use]
    pub const fn id(&self) -> &ArtifactId {
        self.header.id()
    }

    #[must_use]
    pub const fn header(&self) -> &ArtifactHeader {
        &self.header
    }

    #[must_use]
    pub const fn context(&self) -> &ArtifactRef {
        &self.context
    }

    #[must_use]
    pub const fn by_qoi(&self) -> &BTreeMap<QoiId, QoiValidationPlan> {
        &self.by_qoi
    }
}

/// Provenance class of observations. Only `Physical` can validate physics.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ExperimentOrigin {
    Physical {
        apparatus_id: ArtifactId,
        facility_id: ArtifactId,
    },
    SyntheticHighFidelity {
        producer: ArtifactId,
    },
    SecondImplementation {
        producer: ArtifactId,
    },
}

impl ExperimentOrigin {
    #[must_use]
    pub const fn is_physical(&self) -> bool {
        matches!(self, Self::Physical { .. })
    }
}

/// Calibration evidence for one instrument.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InstrumentCalibration {
    instrument_id: ArtifactId,
    certificate_hash: ContentHash,
    current: bool,
}

impl InstrumentCalibration {
    #[must_use]
    pub fn new(instrument_id: ArtifactId, certificate_hash: ContentHash, current: bool) -> Self {
        Self {
            instrument_id,
            certificate_hash,
            current,
        }
    }

    #[must_use]
    pub const fn instrument_id(&self) -> &ArtifactId {
        &self.instrument_id
    }

    #[must_use]
    pub const fn certificate_hash(&self) -> ContentHash {
        self.certificate_hash
    }

    #[must_use]
    pub const fn current(&self) -> bool {
        self.current
    }
}

/// Explicit clock topology for an experiment.
#[derive(Debug, Clone, PartialEq)]
pub enum ClockSynchronization {
    SingleClock {
        clock_id: ArtifactId,
    },
    Synchronized {
        clock_ids: Vec<ArtifactId>,
        method: String,
        max_skew_seconds: f64,
        evidence_hash: ContentHash,
    },
}

impl ClockSynchronization {
    pub fn synchronized(
        mut clock_ids: Vec<ArtifactId>,
        method: impl Into<String>,
        max_skew_seconds: f64,
        evidence_hash: ContentHash,
    ) -> Result<Self, VvErrors> {
        let method = method.into();
        validate_text(&method, "experiment.clocks.method")?;
        let clock_count = clock_ids.len();
        clock_ids.sort();
        clock_ids.dedup();
        if clock_ids.len() < 2
            || clock_ids.len() > MAX_VV_ITEMS
            || clock_ids.len() != clock_count
            || !max_skew_seconds.is_finite()
            || max_skew_seconds < 0.0
            || !evidence_hash.as_bytes().iter().any(|byte| *byte != 0)
        {
            return Err(invalid(
                VvRule::ExperimentClockSynchronization,
                None,
                None,
                "experiment.clocks",
                "synchronized clocks need at least two clocks and a finite non-negative skew",
            ));
        }
        Ok(Self::Synchronized {
            clock_ids,
            method,
            max_skew_seconds,
            evidence_hash,
        })
    }

    fn validated_canonical(self) -> Result<Self, VvErrors> {
        match self {
            Self::SingleClock { clock_id } => Ok(Self::SingleClock { clock_id }),
            Self::Synchronized {
                clock_ids,
                method,
                max_skew_seconds,
                evidence_hash,
            } => Self::synchronized(clock_ids, method, max_skew_seconds, evidence_hash),
        }
    }
}

/// Symmetric covariance stored as a lower triangle in declared QoI order.
#[derive(Debug, Clone, PartialEq)]
pub struct CovarianceMatrix {
    dimension: usize,
    lower_triangle: Vec<f64>,
}

impl CovarianceMatrix {
    pub fn try_new(dimension: usize, lower_triangle: Vec<f64>) -> Result<Self, VvErrors> {
        let expected = dimension
            .checked_mul(dimension.saturating_add(1))
            .and_then(|value| value.checked_div(2));
        if dimension == 0
            || dimension > MAX_VV_MATRIX_DIMENSION
            || expected != Some(lower_triangle.len())
            || lower_triangle.iter().any(|value| !value.is_finite())
        {
            return Err(invalid(
                VvRule::ExperimentRepeatabilityCovariance,
                None,
                None,
                "experiment.covariance",
                "covariance must be finite and have exactly n(n+1)/2 lower-triangle entries",
            ));
        }
        let candidate = Self {
            dimension,
            lower_triangle,
        };
        if (0..dimension).any(|index| candidate.get(index, index) < 0.0) {
            return Err(invalid(
                VvRule::ExperimentRepeatabilityCovariance,
                None,
                None,
                "experiment.covariance",
                "covariance diagonal entries must be non-negative",
            ));
        }
        if !candidate.is_positive_semidefinite() {
            return Err(invalid(
                VvRule::ExperimentRepeatabilityCovariance,
                None,
                None,
                "experiment.covariance",
                "covariance must be positive semidefinite",
            ));
        }
        Ok(candidate)
    }

    #[must_use]
    pub const fn dimension(&self) -> usize {
        self.dimension
    }

    #[must_use]
    pub fn lower_triangle(&self) -> &[f64] {
        &self.lower_triangle
    }

    fn get(&self, row: usize, column: usize) -> f64 {
        let (row, column) = if row >= column {
            (row, column)
        } else {
            (column, row)
        };
        self.lower_triangle[row * (row + 1) / 2 + column]
    }

    fn is_positive_semidefinite(&self) -> bool {
        // Deterministic LDL^T. A negative pivot always refuses; clamping even a
        // tiny negative pivot could turn a false scientific certificate green.
        let n = self.dimension;
        let mut lower = vec![0.0; n * n];
        let mut diagonal = vec![0.0; n];
        for row in 0..n {
            for column in 0..row {
                let mut value = self.get(row, column);
                for previous in 0..column {
                    value -= lower[row * n + previous]
                        * diagonal[previous]
                        * lower[column * n + previous];
                }
                if diagonal[column] == 0.0 {
                    if value != 0.0 {
                        return false;
                    }
                    lower[row * n + column] = 0.0;
                } else {
                    lower[row * n + column] = value / diagonal[column];
                }
            }
            let mut pivot = self.get(row, row);
            for previous in 0..row {
                let value = lower[row * n + previous];
                pivot -= value * value * diagonal[previous];
            }
            if pivot < 0.0 || !pivot.is_finite() {
                return false;
            }
            diagonal[row] = pivot;
            lower[row * n + row] = 1.0;
        }
        true
    }
}

/// Repeatability sample count and covariance.
#[derive(Debug, Clone, PartialEq)]
pub struct RepeatabilitySummary {
    replicates: u32,
    covariance: CovarianceMatrix,
}

impl RepeatabilitySummary {
    pub fn try_new(replicates: u32, covariance: CovarianceMatrix) -> Result<Self, VvErrors> {
        if replicates < 2 {
            return Err(invalid(
                VvRule::ExperimentRepeatabilityCovariance,
                None,
                None,
                "experiment.repeatability",
                "repeatability requires at least two replicates",
            ));
        }
        Ok(Self {
            replicates,
            covariance,
        })
    }

    #[must_use]
    pub const fn replicates(&self) -> u32 {
        self.replicates
    }

    #[must_use]
    pub const fn covariance(&self) -> &CovarianceMatrix {
        &self.covariance
    }
}

/// Exact source bytes and custody evidence for an experimental dataset.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DataAuthenticity {
    source_bytes_hash: ContentHash,
    custody_receipt_hash: ContentHash,
    authenticated: bool,
}

impl DataAuthenticity {
    #[must_use]
    pub fn new(
        source_bytes_hash: ContentHash,
        custody_receipt_hash: ContentHash,
        authenticated: bool,
    ) -> Self {
        Self {
            source_bytes_hash,
            custody_receipt_hash,
            authenticated,
        }
    }

    #[must_use]
    pub const fn source_bytes_hash(&self) -> ContentHash {
        self.source_bytes_hash
    }

    #[must_use]
    pub const fn custody_receipt_hash(&self) -> ContentHash {
        self.custody_receipt_hash
    }

    #[must_use]
    pub const fn authenticated(&self) -> bool {
        self.authenticated
    }
}

/// Physical or synthetic observation artifact with metrology and authenticity.
#[derive(Debug, Clone, PartialEq)]
pub struct ExperimentArtifact {
    header: ArtifactHeader,
    dataset_id: ArtifactId,
    origin: ExperimentOrigin,
    qois: BTreeSet<QoiId>,
    observation_ids: BTreeSet<ObservationId>,
    observations_hash: ContentHash,
    instruments: Vec<InstrumentCalibration>,
    clocks: ClockSynchronization,
    repeatability: RepeatabilitySummary,
    authenticity: DataAuthenticity,
}

impl ExperimentArtifact {
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        header: ArtifactHeader,
        dataset_id: ArtifactId,
        origin: ExperimentOrigin,
        qois: Vec<QoiId>,
        observation_ids: Vec<ObservationId>,
        observations_hash: ContentHash,
        instruments: Vec<InstrumentCalibration>,
        clocks: ClockSynchronization,
        repeatability: RepeatabilitySummary,
        authenticity: DataAuthenticity,
    ) -> Result<Self, VvErrors> {
        if qois.is_empty() || qois.len() > MAX_VV_ITEMS {
            return Err(invalid(
                VvRule::QoiDependencyClosed,
                Some(header.id().as_str()),
                None,
                "experiment.qois",
                "experiment QoIs must be explicit and bounded",
            ));
        }
        let qoi_count = qois.len();
        let qois = qois.into_iter().collect::<BTreeSet<_>>();
        if qois.len() != qoi_count {
            return Err(invalid(
                VvRule::QoiDependencyClosed,
                Some(header.id().as_str()),
                None,
                "experiment.qois",
                "experiment QoI identities must be unique",
            ));
        }
        if repeatability.covariance.dimension != qois.len() {
            return Err(invalid(
                VvRule::ExperimentRepeatabilityCovariance,
                Some(header.id().as_str()),
                None,
                "experiment.covariance",
                "covariance dimension must equal the canonical QoI count",
            ));
        }
        if observation_ids.is_empty() || observation_ids.len() > MAX_VV_ITEMS {
            return Err(invalid(
                VvRule::SchemaCardinality,
                Some(header.id().as_str()),
                None,
                "experiment.observations",
                "experiment observation identities must be explicit and bounded",
            ));
        }
        let observation_count = observation_ids.len();
        let observation_ids = observation_ids.into_iter().collect::<BTreeSet<_>>();
        if observation_ids.len() != observation_count {
            return Err(invalid(
                VvRule::SchemaIdentity,
                Some(header.id().as_str()),
                None,
                "experiment.observations",
                "observation identities must be unique",
            ));
        }
        if instruments.is_empty()
            || instruments.len() > MAX_VV_ITEMS
            || instruments.iter().any(|instrument| !instrument.current)
        {
            return Err(invalid(
                VvRule::ExperimentInstrumentCalibration,
                Some(header.id().as_str()),
                None,
                "experiment.instruments",
                "every experiment needs current calibration evidence for every instrument",
            ));
        }
        let unique_instruments = instruments
            .iter()
            .map(|instrument| instrument.instrument_id())
            .collect::<BTreeSet<_>>();
        if unique_instruments.len() != instruments.len() {
            return Err(invalid(
                VvRule::ExperimentInstrumentCalibration,
                Some(header.id().as_str()),
                None,
                "experiment.instruments",
                "instrument calibration rows must be unique by instrument identity",
            ));
        }
        let mut instruments = instruments;
        instruments.sort_by(|left, right| left.instrument_id.cmp(&right.instrument_id));
        let clocks = clocks.validated_canonical().map_err(|_| {
            invalid(
                VvRule::ExperimentClockSynchronization,
                Some(header.id().as_str()),
                None,
                "experiment.clocks",
                "clock synchronization must be structurally valid and canonicalizable",
            )
        })?;
        if !authenticity.authenticated {
            return Err(invalid(
                VvRule::ExperimentDataAuthenticity,
                Some(header.id().as_str()),
                None,
                "experiment.authenticity",
                "dataset authenticity must be admitted by the configured policy",
            ));
        }
        Ok(Self {
            header,
            dataset_id,
            origin,
            qois,
            observation_ids,
            observations_hash,
            instruments,
            clocks,
            repeatability,
            authenticity,
        })
    }

    #[must_use]
    pub const fn id(&self) -> &ArtifactId {
        self.header.id()
    }

    #[must_use]
    pub const fn header(&self) -> &ArtifactHeader {
        &self.header
    }

    #[must_use]
    pub const fn dataset_id(&self) -> &ArtifactId {
        &self.dataset_id
    }

    #[must_use]
    pub const fn origin(&self) -> &ExperimentOrigin {
        &self.origin
    }

    #[must_use]
    pub const fn qois(&self) -> &BTreeSet<QoiId> {
        &self.qois
    }

    #[must_use]
    pub const fn observation_ids(&self) -> &BTreeSet<ObservationId> {
        &self.observation_ids
    }

    #[must_use]
    pub const fn observations_hash(&self) -> ContentHash {
        self.observations_hash
    }

    #[must_use]
    pub fn instruments(&self) -> &[InstrumentCalibration] {
        &self.instruments
    }

    #[must_use]
    pub const fn clocks(&self) -> &ClockSynchronization {
        &self.clocks
    }

    #[must_use]
    pub const fn repeatability(&self) -> &RepeatabilitySummary {
        &self.repeatability
    }

    #[must_use]
    pub const fn authenticity(&self) -> &DataAuthenticity {
        &self.authenticity
    }
}

fn commitment_for_blind_rows(
    preregistration_hash: ContentHash,
    rows: &BTreeSet<ObservationId>,
) -> ContentHash {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(preregistration_hash.as_bytes());
    bytes.extend_from_slice(&(rows.len() as u64).to_le_bytes());
    for row in rows {
        bytes.extend_from_slice(&(row.as_str().len() as u64).to_le_bytes());
        bytes.extend_from_slice(row.as_str().as_bytes());
    }
    fs_blake3::hash_domain("org.frankensim.fs-evidence.vv-blind-holdout.v1", &bytes)
}

/// Authority record required before blind holdout rows become validation input.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct BlindReleaseReceipt {
    split: ArtifactRef,
    blind_commitment: ContentHash,
    authority_receipt_hash: ContentHash,
}

impl BlindReleaseReceipt {
    pub fn new(
        split: ArtifactRef,
        blind_commitment: ContentHash,
        authority_receipt_hash: ContentHash,
    ) -> Result<Self, VvErrors> {
        if split.kind != ArtifactKind::CalibrationSplit
            || !authority_receipt_hash
                .as_bytes()
                .iter()
                .any(|byte| *byte != 0)
        {
            return Err(invalid(
                VvRule::SplitBlindHoldoutSealed,
                Some(split.id().as_str()),
                None,
                "blind_release",
                "blind release must bind an exact split and non-empty authority evidence",
            ));
        }
        Ok(Self {
            split,
            blind_commitment,
            authority_receipt_hash,
        })
    }

    #[must_use]
    pub const fn split(&self) -> &ArtifactRef {
        &self.split
    }

    #[must_use]
    pub const fn blind_commitment(&self) -> ContentHash {
        self.blind_commitment
    }

    #[must_use]
    pub const fn authority_receipt_hash(&self) -> ContentHash {
        self.authority_receipt_hash
    }
}

/// Evidence-bearing split partition. Calibration is deliberately absent.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum EvidencePartition {
    Validation,
    BlindHoldout { release: BlindReleaseReceipt },
}

/// Sealed observation subset that can be consumed by validation metrics.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct ObservationSelection {
    split: ArtifactRef,
    ids: BTreeSet<ObservationId>,
    partition: EvidencePartition,
}

impl ObservationSelection {
    #[must_use]
    pub const fn split(&self) -> &ArtifactRef {
        &self.split
    }

    #[must_use]
    pub const fn ids(&self) -> &BTreeSet<ObservationId> {
        &self.ids
    }

    #[must_use]
    pub const fn partition(&self) -> &EvidencePartition {
        &self.partition
    }

    pub(crate) fn from_canonical(
        split: ArtifactRef,
        ids: Vec<ObservationId>,
        partition: EvidencePartition,
    ) -> Result<Self, VvErrors> {
        if split.kind != ArtifactKind::CalibrationSplit
            || ids.is_empty()
            || ids.len() > MAX_VV_ITEMS
        {
            return Err(invalid(
                VvRule::SplitPartitionsDisjoint,
                Some(split.id().as_str()),
                None,
                "selection",
                "canonical selection must name a split and a bounded non-empty row set",
            ));
        }
        let count = ids.len();
        let ids = ids.into_iter().collect::<BTreeSet<_>>();
        if ids.len() != count {
            return Err(invalid(
                VvRule::SplitPartitionsDisjoint,
                Some(split.id().as_str()),
                None,
                "selection.ids",
                "canonical selection rows must be unique",
            ));
        }
        Ok(Self {
            split,
            ids,
            partition,
        })
    }
}

/// Pre-registered calibration, validation, and blind-holdout partition.
#[derive(Debug, Clone, PartialEq)]
pub struct CalibrationSplit {
    header: ArtifactHeader,
    experiment: ArtifactRef,
    preregistration_hash: ContentHash,
    calibration: BTreeSet<ObservationId>,
    validation: BTreeSet<ObservationId>,
    blind_holdout: BTreeSet<ObservationId>,
    blind_commitment: ContentHash,
}

impl CalibrationSplit {
    pub fn try_new(
        header: ArtifactHeader,
        experiment: ArtifactRef,
        preregistration_hash: ContentHash,
        calibration: Vec<ObservationId>,
        validation: Vec<ObservationId>,
        blind_holdout: Vec<ObservationId>,
    ) -> Result<Self, VvErrors> {
        if experiment.kind != ArtifactKind::ExperimentArtifact {
            return Err(invalid(
                VvRule::SplitPartitionsDisjoint,
                Some(header.id().as_str()),
                None,
                "split.experiment",
                "a split must reference one ExperimentArtifact",
            ));
        }
        if calibration.is_empty()
            || validation.is_empty()
            || blind_holdout.is_empty()
            || calibration
                .len()
                .saturating_add(validation.len())
                .saturating_add(blind_holdout.len())
                > MAX_VV_ITEMS
        {
            return Err(invalid(
                VvRule::SplitPartitionsDisjoint,
                Some(header.id().as_str()),
                None,
                "split.partitions",
                "all three partitions must be non-empty and bounded",
            ));
        }
        let calibration_count = calibration.len();
        let validation_count = validation.len();
        let blind_count = blind_holdout.len();
        let calibration = calibration.into_iter().collect::<BTreeSet<_>>();
        let validation = validation.into_iter().collect::<BTreeSet<_>>();
        let blind_holdout = blind_holdout.into_iter().collect::<BTreeSet<_>>();
        if calibration.len() != calibration_count
            || validation.len() != validation_count
            || blind_holdout.len() != blind_count
        {
            return Err(invalid(
                VvRule::SplitPartitionsDisjoint,
                Some(header.id().as_str()),
                None,
                "split.partitions",
                "observation identities must be unique within each partition",
            ));
        }
        if !calibration.is_disjoint(&validation)
            || !calibration.is_disjoint(&blind_holdout)
            || !validation.is_disjoint(&blind_holdout)
        {
            return Err(invalid(
                VvRule::SplitPartitionsDisjoint,
                Some(header.id().as_str()),
                None,
                "split.partitions",
                "calibration, validation, and blind holdout must be pairwise disjoint",
            ));
        }
        let blind_commitment = commitment_for_blind_rows(preregistration_hash, &blind_holdout);
        Ok(Self {
            header,
            experiment,
            preregistration_hash,
            calibration,
            validation,
            blind_holdout,
            blind_commitment,
        })
    }

    #[must_use]
    pub const fn id(&self) -> &ArtifactId {
        self.header.id()
    }

    #[must_use]
    pub const fn header(&self) -> &ArtifactHeader {
        &self.header
    }

    #[must_use]
    pub const fn experiment(&self) -> &ArtifactRef {
        &self.experiment
    }

    #[must_use]
    pub const fn preregistration_hash(&self) -> ContentHash {
        self.preregistration_hash
    }

    #[must_use]
    pub const fn calibration_ids(&self) -> &BTreeSet<ObservationId> {
        &self.calibration
    }

    #[must_use]
    pub const fn validation_ids(&self) -> &BTreeSet<ObservationId> {
        &self.validation
    }

    #[must_use]
    pub fn blind_holdout_len(&self) -> usize {
        self.blind_holdout.len()
    }

    #[must_use]
    pub const fn blind_commitment(&self) -> ContentHash {
        self.blind_commitment
    }

    pub fn validation_selection(
        &self,
        split: ArtifactRef,
        ids: Vec<ObservationId>,
    ) -> Result<ObservationSelection, VvErrors> {
        self.selection(split, ids, EvidencePartition::Validation)
    }

    pub fn blind_selection(
        &self,
        split: ArtifactRef,
        ids: Vec<ObservationId>,
        release: BlindReleaseReceipt,
    ) -> Result<ObservationSelection, VvErrors> {
        self.selection(split, ids, EvidencePartition::BlindHoldout { release })
    }

    fn selection(
        &self,
        split: ArtifactRef,
        ids: Vec<ObservationId>,
        partition: EvidencePartition,
    ) -> Result<ObservationSelection, VvErrors> {
        if split.kind != ArtifactKind::CalibrationSplit || split.id != *self.id() {
            return Err(invalid(
                VvRule::SplitPartitionsDisjoint,
                Some(self.id().as_str()),
                None,
                "selection.split",
                "selection must reference the split that minted it",
            ));
        }
        if ids.is_empty() || ids.len() > MAX_VV_ITEMS {
            return Err(invalid(
                VvRule::SplitPartitionsDisjoint,
                Some(self.id().as_str()),
                None,
                "selection.ids",
                "a validation selection must be bounded and non-empty",
            ));
        }
        let id_count = ids.len();
        let ids = ids.into_iter().collect::<BTreeSet<_>>();
        if ids.len() != id_count {
            return Err(invalid(
                VvRule::SplitPartitionsDisjoint,
                Some(self.id().as_str()),
                None,
                "selection.ids",
                "observation identities must be unique within a selection",
            ));
        }
        let allowed = match &partition {
            EvidencePartition::Validation => &self.validation,
            EvidencePartition::BlindHoldout { release }
                if release.split.id == *self.id()
                    && release.split.hash == split.hash
                    && release.blind_commitment == self.blind_commitment =>
            {
                &self.blind_holdout
            }
            EvidencePartition::BlindHoldout { .. } => {
                return Err(invalid(
                    VvRule::SplitBlindHoldoutSealed,
                    Some(self.id().as_str()),
                    None,
                    "selection.blind_release",
                    "blind rows require a release bound to this split and commitment",
                ));
            }
        };
        if !ids.is_subset(allowed) {
            let calibration_reuse = ids.iter().any(|id| self.calibration.contains(id));
            return Err(invalid(
                if calibration_reuse {
                    VvRule::ValidationRequiresPhysicalReferent
                } else {
                    VvRule::SplitPartitionsDisjoint
                },
                Some(self.id().as_str()),
                None,
                "selection.ids",
                if calibration_reuse {
                    "calibration observations cannot be reused as validation evidence"
                } else {
                    "selection contains observations outside its declared partition"
                },
            ));
        }
        Ok(ObservationSelection {
            split,
            ids,
            partition,
        })
    }

    fn all_ids(&self) -> BTreeSet<ObservationId> {
        self.calibration
            .union(&self.validation)
            .cloned()
            .chain(self.blind_holdout.iter().cloned())
            .collect()
    }

    pub(crate) fn blind_holdout_ids_for_codec(&self) -> &BTreeSet<ObservationId> {
        &self.blind_holdout
    }
}

fn next_up(value: f64) -> f64 {
    if value.is_nan() || value == f64::INFINITY {
        value
    } else if value == 0.0 {
        f64::from_bits(1)
    } else {
        let bits = value.to_bits();
        f64::from_bits(if value > 0.0 { bits + 1 } else { bits - 1 })
    }
}

/// One numerical uncertainty component and its retained evidence.
#[derive(Debug, Clone, PartialEq)]
pub struct NumericalUncertainty {
    half_width: f64,
    evidence_hash: ContentHash,
}

impl NumericalUncertainty {
    pub fn try_new(half_width: f64, evidence_hash: ContentHash) -> Result<Self, VvErrors> {
        if !half_width.is_finite() || half_width < 0.0 {
            return Err(invalid(
                VvRule::SolutionVerificationComplete,
                None,
                None,
                "solution_verification.uncertainty",
                "numerical uncertainty must be finite and non-negative",
            ));
        }
        Ok(Self {
            half_width,
            evidence_hash,
        })
    }

    #[must_use]
    pub const fn half_width(&self) -> f64 {
        self.half_width
    }

    #[must_use]
    pub const fn evidence_hash(&self) -> ContentHash {
        self.evidence_hash
    }
}

/// Mesh/time/nonlinear/iterative numerical uncertainty for one exact solve/QoI.
#[derive(Debug, Clone, PartialEq)]
pub struct SolutionVerificationReceipt {
    header: ArtifactHeader,
    solve_id: ArtifactId,
    qoi: QoiId,
    unit: UnitId,
    mesh: NumericalUncertainty,
    time: NumericalUncertainty,
    nonlinear: NumericalUncertainty,
    iterative: NumericalUncertainty,
    combined_half_width: f64,
}

impl SolutionVerificationReceipt {
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        header: ArtifactHeader,
        solve_id: ArtifactId,
        qoi: QoiId,
        unit: UnitId,
        mesh: NumericalUncertainty,
        time: NumericalUncertainty,
        nonlinear: NumericalUncertainty,
        iterative: NumericalUncertainty,
    ) -> Result<Self, VvErrors> {
        let mut combined = 0.0;
        for component in [&mesh, &time, &nonlinear, &iterative] {
            combined = next_up(combined + component.half_width);
        }
        if !combined.is_finite() {
            return Err(invalid(
                VvRule::SolutionVerificationComplete,
                Some(header.id().as_str()),
                Some(qoi.as_str()),
                "solution_verification.combined",
                "numerical uncertainty total overflowed",
            ));
        }
        Ok(Self {
            header,
            solve_id,
            qoi,
            unit,
            mesh,
            time,
            nonlinear,
            iterative,
            combined_half_width: combined,
        })
    }

    #[must_use]
    pub const fn id(&self) -> &ArtifactId {
        self.header.id()
    }

    #[must_use]
    pub const fn header(&self) -> &ArtifactHeader {
        &self.header
    }

    #[must_use]
    pub const fn solve_id(&self) -> &ArtifactId {
        &self.solve_id
    }

    #[must_use]
    pub const fn qoi(&self) -> &QoiId {
        &self.qoi
    }

    #[must_use]
    pub const fn unit(&self) -> &UnitId {
        &self.unit
    }

    #[must_use]
    pub const fn mesh(&self) -> &NumericalUncertainty {
        &self.mesh
    }

    #[must_use]
    pub const fn time(&self) -> &NumericalUncertainty {
        &self.time
    }

    #[must_use]
    pub const fn nonlinear(&self) -> &NumericalUncertainty {
        &self.nonlinear
    }

    #[must_use]
    pub const fn iterative(&self) -> &NumericalUncertainty {
        &self.iterative
    }

    #[must_use]
    pub const fn combined_half_width(&self) -> f64 {
        self.combined_half_width
    }
}

/// Exact target of a QoI-specific dependency.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EvidenceTarget {
    VvArtifact(ArtifactRef),
    External {
        family: ArtifactId,
        id: ArtifactId,
        hash: ContentHash,
    },
}

impl EvidenceTarget {
    #[must_use]
    pub fn hash(&self) -> ContentHash {
        match self {
            Self::VvArtifact(reference) => reference.hash,
            Self::External { hash, .. } => *hash,
        }
    }
}

/// Semantic role of an exact dependency.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DependencyRole {
    CodeVerification,
    SolutionVerification,
    PhysicalValidation,
    ModelDiscrepancy,
    ParameterData,
    PosteriorPredictive,
    ProcessConformance,
}

/// One load-bearing QoI-specific dependency edge.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct EvidenceDependency {
    qoi: QoiId,
    role: DependencyRole,
    target: EvidenceTarget,
    observations: Option<ObservationSelection>,
}

impl EvidenceDependency {
    #[must_use]
    pub fn new(qoi: QoiId, role: DependencyRole, target: EvidenceTarget) -> Self {
        Self {
            qoi,
            role,
            target,
            observations: None,
        }
    }

    #[must_use]
    pub fn physical_validation(
        qoi: QoiId,
        experiment: ArtifactRef,
        observations: ObservationSelection,
    ) -> Self {
        Self {
            qoi,
            role: DependencyRole::PhysicalValidation,
            target: EvidenceTarget::VvArtifact(experiment),
            observations: Some(observations),
        }
    }

    #[must_use]
    pub const fn qoi(&self) -> &QoiId {
        &self.qoi
    }

    #[must_use]
    pub const fn role(&self) -> DependencyRole {
        self.role
    }

    #[must_use]
    pub const fn target(&self) -> &EvidenceTarget {
        &self.target
    }

    #[must_use]
    pub const fn observations(&self) -> Option<&ObservationSelection> {
        self.observations.as_ref()
    }
}

/// Six required prediction-uncertainty categories.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PredictionUncertaintyKind {
    ModelForm,
    Parameter,
    Numerical,
    Data,
    Aleatory,
    Epistemic,
}

impl PredictionUncertaintyKind {
    pub const ALL: [Self; 6] = [
        Self::ModelForm,
        Self::Parameter,
        Self::Numerical,
        Self::Data,
        Self::Aleatory,
        Self::Epistemic,
    ];
}

/// One uncertainty term, interpreted according to the waterfall mode.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UncertaintyTerm {
    kind: PredictionUncertaintyKind,
    magnitude_bits: u64,
    source: EvidenceTarget,
}

impl UncertaintyTerm {
    pub fn try_new(
        kind: PredictionUncertaintyKind,
        magnitude: f64,
        source: EvidenceTarget,
    ) -> Result<Self, VvErrors> {
        if !magnitude.is_finite() || magnitude < 0.0 {
            return Err(invalid(
                VvRule::WaterfallArithmetic,
                None,
                None,
                "waterfall.term",
                "uncertainty magnitudes must be finite and non-negative",
            ));
        }
        Ok(Self {
            kind,
            magnitude_bits: magnitude.to_bits(),
            source,
        })
    }

    #[must_use]
    pub const fn kind(&self) -> PredictionUncertaintyKind {
        self.kind
    }

    #[must_use]
    pub fn magnitude(&self) -> f64 {
        f64::from_bits(self.magnitude_bits)
    }

    #[must_use]
    pub const fn source(&self) -> &EvidenceTarget {
        &self.source
    }
}

/// Dense, declared correlation matrix for probabilistic waterfalls.
#[derive(Debug, Clone, PartialEq)]
pub struct CorrelationMatrix {
    dimension: usize,
    values: Vec<f64>,
}

impl CorrelationMatrix {
    pub fn try_new(dimension: usize, values: Vec<f64>) -> Result<Self, VvErrors> {
        if dimension == 0
            || dimension > MAX_VV_MATRIX_DIMENSION
            || dimension.checked_mul(dimension) != Some(values.len())
            || values.iter().any(|value| !value.is_finite())
        {
            return Err(invalid(
                VvRule::WaterfallDependenceDeclared,
                None,
                None,
                "waterfall.correlation",
                "correlation matrix must be finite, square, and bounded",
            ));
        }
        let tolerance = 64.0 * f64::EPSILON * dimension as f64;
        for row in 0..dimension {
            for column in 0..dimension {
                let value = values[row * dimension + column];
                if value < -1.0 || value > 1.0 {
                    return Err(invalid(
                        VvRule::WaterfallDependenceDeclared,
                        None,
                        None,
                        "waterfall.correlation",
                        "correlations must lie in [-1, 1]",
                    ));
                }
                if (value - values[column * dimension + row]).abs() > tolerance {
                    return Err(invalid(
                        VvRule::WaterfallDependenceDeclared,
                        None,
                        None,
                        "waterfall.correlation",
                        "correlation matrix must be symmetric",
                    ));
                }
            }
            if (values[row * dimension + row] - 1.0).abs() > tolerance {
                return Err(invalid(
                    VvRule::WaterfallDependenceDeclared,
                    None,
                    None,
                    "waterfall.correlation",
                    "correlation diagonal must equal one",
                ));
            }
        }
        // Reuse covariance admission for the positive-semidefinite gate.
        let mut lower = Vec::with_capacity(dimension * (dimension + 1) / 2);
        for row in 0..dimension {
            for column in 0..=row {
                lower.push(values[row * dimension + column]);
            }
        }
        CovarianceMatrix::try_new(dimension, lower).map_err(|_| {
            invalid(
                VvRule::WaterfallDependenceDeclared,
                None,
                None,
                "waterfall.correlation",
                "correlation matrix must be positive semidefinite",
            )
        })?;
        Ok(Self { dimension, values })
    }

    #[must_use]
    pub const fn dimension(&self) -> usize {
        self.dimension
    }

    #[must_use]
    pub fn values(&self) -> &[f64] {
        &self.values
    }
}

/// Explicit interpretation of waterfall magnitudes.
#[derive(Debug, Clone, PartialEq)]
pub enum WaterfallMode {
    GuaranteedBound,
    Probabilistic {
        confidence: f64,
        dependence: CorrelationMatrix,
    },
}

/// Six-source uncertainty budget with a derived total.
#[derive(Debug, Clone, PartialEq)]
pub struct UncertaintyWaterfall {
    qoi: QoiId,
    unit: UnitId,
    mode: WaterfallMode,
    terms: Vec<UncertaintyTerm>,
    total: f64,
}

impl UncertaintyWaterfall {
    pub fn try_new(
        qoi: QoiId,
        unit: UnitId,
        mode: WaterfallMode,
        terms: Vec<UncertaintyTerm>,
    ) -> Result<Self, VvErrors> {
        if terms.len() != PredictionUncertaintyKind::ALL.len() {
            return Err(invalid(
                VvRule::WaterfallModeDeclared,
                None,
                Some(qoi.as_str()),
                "waterfall.terms",
                "waterfall needs exactly the six declared uncertainty categories",
            ));
        }
        let mut terms = terms;
        terms.sort_by_key(UncertaintyTerm::kind);
        if !terms
            .iter()
            .map(UncertaintyTerm::kind)
            .eq(PredictionUncertaintyKind::ALL)
        {
            return Err(invalid(
                VvRule::WaterfallModeDeclared,
                None,
                Some(qoi.as_str()),
                "waterfall.terms",
                "waterfall categories must be complete and unique",
            ));
        }
        let total = match &mode {
            WaterfallMode::GuaranteedBound => {
                let mut total = 0.0;
                for term in &terms {
                    total = next_up(total + term.magnitude());
                }
                total
            }
            WaterfallMode::Probabilistic {
                confidence,
                dependence,
            } => {
                if !confidence.is_finite()
                    || *confidence <= 0.0
                    || *confidence >= 1.0
                    || dependence.dimension != terms.len()
                {
                    return Err(invalid(
                        VvRule::WaterfallDependenceDeclared,
                        None,
                        Some(qoi.as_str()),
                        "waterfall.mode",
                        "probabilistic mode needs a confidence in (0,1) and a 6x6 dependence matrix",
                    ));
                }
                let mut variance = 0.0;
                for row in 0..terms.len() {
                    for column in 0..terms.len() {
                        variance += terms[row].magnitude()
                            * dependence.values[row * terms.len() + column]
                            * terms[column].magnitude();
                    }
                }
                if variance < -64.0 * f64::EPSILON || !variance.is_finite() {
                    return Err(invalid(
                        VvRule::WaterfallArithmetic,
                        None,
                        Some(qoi.as_str()),
                        "waterfall.total",
                        "probabilistic variance must be finite and non-negative",
                    ));
                }
                next_up(variance.max(0.0).sqrt())
            }
        };
        if !total.is_finite() {
            return Err(invalid(
                VvRule::WaterfallArithmetic,
                None,
                Some(qoi.as_str()),
                "waterfall.total",
                "uncertainty total overflowed",
            ));
        }
        Ok(Self {
            qoi,
            unit,
            mode,
            terms,
            total,
        })
    }

    #[must_use]
    pub const fn qoi(&self) -> &QoiId {
        &self.qoi
    }

    #[must_use]
    pub const fn unit(&self) -> &UnitId {
        &self.unit
    }

    #[must_use]
    pub const fn mode(&self) -> &WaterfallMode {
        &self.mode
    }

    #[must_use]
    pub fn terms(&self) -> &[UncertaintyTerm] {
        &self.terms
    }

    #[must_use]
    pub const fn total(&self) -> f64 {
        self.total
    }
}

/// Comparison to observations with both experimental and numerical uncertainty.
#[derive(Debug, Clone, PartialEq)]
pub struct ValidationMetric {
    name: ArtifactId,
    qoi: QoiId,
    observations: ObservationSelection,
    observed: f64,
    predicted: f64,
    experimental_uncertainty: f64,
    numerical_uncertainty: f64,
    combined_uncertainty: f64,
}

impl ValidationMetric {
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        name: ArtifactId,
        qoi: QoiId,
        observations: ObservationSelection,
        observed: f64,
        predicted: f64,
        experimental_uncertainty: f64,
        numerical_uncertainty: f64,
    ) -> Result<Self, VvErrors> {
        if !observed.is_finite()
            || !predicted.is_finite()
            || !experimental_uncertainty.is_finite()
            || experimental_uncertainty < 0.0
            || !numerical_uncertainty.is_finite()
            || numerical_uncertainty < 0.0
        {
            return Err(invalid(
                VvRule::ValidationMetricUncertainty,
                None,
                Some(qoi.as_str()),
                "validation_metric",
                "metric values and both uncertainty contributions must be finite and non-negative",
            ));
        }
        let combined_uncertainty = next_up(experimental_uncertainty + numerical_uncertainty);
        if !combined_uncertainty.is_finite() {
            return Err(invalid(
                VvRule::ValidationMetricUncertainty,
                None,
                Some(qoi.as_str()),
                "validation_metric.combined_uncertainty",
                "combined metric uncertainty overflowed",
            ));
        }
        Ok(Self {
            name,
            qoi,
            observations,
            observed,
            predicted,
            experimental_uncertainty,
            numerical_uncertainty,
            combined_uncertainty,
        })
    }

    #[must_use]
    pub const fn name(&self) -> &ArtifactId {
        &self.name
    }

    #[must_use]
    pub const fn qoi(&self) -> &QoiId {
        &self.qoi
    }

    #[must_use]
    pub const fn observations(&self) -> &ObservationSelection {
        &self.observations
    }

    #[must_use]
    pub const fn observed(&self) -> f64 {
        self.observed
    }

    #[must_use]
    pub const fn predicted(&self) -> f64 {
        self.predicted
    }

    #[must_use]
    pub const fn experimental_uncertainty(&self) -> f64 {
        self.experimental_uncertainty
    }

    #[must_use]
    pub const fn numerical_uncertainty(&self) -> f64 {
        self.numerical_uncertainty
    }

    #[must_use]
    pub const fn combined_uncertainty(&self) -> f64 {
        self.combined_uncertainty
    }
}

/// One held-out posterior-predictive diagnostic.
#[derive(Debug, Clone, PartialEq)]
pub struct PosteriorPredictiveCheck {
    name: ArtifactId,
    qoi: QoiId,
    observations: ObservationSelection,
    tail_probability: f64,
    minimum_tail_probability: f64,
    artifact_hash: ContentHash,
}

impl PosteriorPredictiveCheck {
    pub fn try_new(
        name: ArtifactId,
        qoi: QoiId,
        observations: ObservationSelection,
        tail_probability: f64,
        minimum_tail_probability: f64,
        artifact_hash: ContentHash,
    ) -> Result<Self, VvErrors> {
        if !tail_probability.is_finite()
            || !(0.0..=1.0).contains(&tail_probability)
            || !minimum_tail_probability.is_finite()
            || minimum_tail_probability <= 0.0
            || minimum_tail_probability >= 1.0
        {
            return Err(invalid(
                VvRule::ValidationMetricUncertainty,
                None,
                Some(qoi.as_str()),
                "posterior_predictive",
                "posterior-predictive probabilities must lie in their declared domains",
            ));
        }
        Ok(Self {
            name,
            qoi,
            observations,
            tail_probability,
            minimum_tail_probability,
            artifact_hash,
        })
    }

    #[must_use]
    pub const fn name(&self) -> &ArtifactId {
        &self.name
    }

    #[must_use]
    pub const fn qoi(&self) -> &QoiId {
        &self.qoi
    }

    #[must_use]
    pub const fn observations(&self) -> &ObservationSelection {
        &self.observations
    }

    #[must_use]
    pub const fn tail_probability(&self) -> f64 {
        self.tail_probability
    }

    #[must_use]
    pub const fn minimum_tail_probability(&self) -> f64 {
        self.minimum_tail_probability
    }

    #[must_use]
    pub const fn artifact_hash(&self) -> ContentHash {
        self.artifact_hash
    }

    #[must_use]
    pub fn passed(&self) -> bool {
        self.tail_probability >= self.minimum_tail_probability
    }
}

/// Independent report axes; these are categories, never numbers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EvidenceAxis {
    CodeVerification,
    SolutionVerification,
    NumericalUncertainty,
    ParameterDataUncertainty,
    ModelFormValidation,
    PredictionDomainRelevance,
    ComparisonToExperiment,
}

impl EvidenceAxis {
    pub const ALL: [Self; 7] = [
        Self::CodeVerification,
        Self::SolutionVerification,
        Self::NumericalUncertainty,
        Self::ParameterDataUncertainty,
        Self::ModelFormValidation,
        Self::PredictionDomainRelevance,
        Self::ComparisonToExperiment,
    ];
}

/// Categorical status of one evidence axis.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EvidenceAxisStatus {
    Present { artifacts: Vec<ContentHash> },
    Missing { reason: String },
    Refused { rule: VvRule, reason: String },
}

/// Complete categorical evidence-axis report with no numeric score API.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EvidenceAxes {
    axes: BTreeMap<EvidenceAxis, EvidenceAxisStatus>,
}

impl EvidenceAxes {
    pub fn try_new(rows: Vec<(EvidenceAxis, EvidenceAxisStatus)>) -> Result<Self, VvErrors> {
        if rows.len() != EvidenceAxis::ALL.len() {
            return Err(invalid(
                VvRule::ColorCategoricalOnly,
                None,
                None,
                "evidence_axes",
                "all seven categorical evidence axes are required",
            ));
        }
        let mut axes = BTreeMap::new();
        for (axis, mut status) in rows {
            match &mut status {
                EvidenceAxisStatus::Present { artifacts } => {
                    if artifacts.is_empty() || artifacts.len() > MAX_VV_ITEMS {
                        return Err(invalid(
                            VvRule::ColorCategoricalOnly,
                            None,
                            None,
                            "evidence_axes.present",
                            "a present category needs bounded artifact evidence",
                        ));
                    }
                    artifacts.sort();
                    artifacts.dedup();
                }
                EvidenceAxisStatus::Missing { reason }
                | EvidenceAxisStatus::Refused { reason, .. } => {
                    validate_text(reason, "evidence_axes.reason")?;
                }
            }
            if axes.insert(axis, status).is_some() {
                return Err(invalid(
                    VvRule::ColorCategoricalOnly,
                    None,
                    None,
                    "evidence_axes",
                    "evidence axes must be unique",
                ));
            }
        }
        Ok(Self { axes })
    }

    #[must_use]
    pub const fn axes(&self) -> &BTreeMap<EvidenceAxis, EvidenceAxisStatus> {
        &self.axes
    }
}

/// QoI-specific physical prediction assessment.
#[derive(Debug, Clone, PartialEq)]
pub struct PredictionAssessment {
    header: ArtifactHeader,
    context: ArtifactRef,
    validation_plan: ArtifactRef,
    qoi: QoiId,
    dependencies: Vec<EvidenceDependency>,
    waterfall: UncertaintyWaterfall,
    validation_metrics: Vec<ValidationMetric>,
    posterior_checks: Vec<PosteriorPredictiveCheck>,
    applicability_point: ApplicabilityPoint,
    applicability: ApplicabilityDecision,
    evidence_axes: EvidenceAxes,
    assumption_checks: BTreeMap<AssumptionId, bool>,
}

impl PredictionAssessment {
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        header: ArtifactHeader,
        context: ArtifactRef,
        validation_plan: ArtifactRef,
        qoi: QoiId,
        dependencies: Vec<EvidenceDependency>,
        waterfall: UncertaintyWaterfall,
        validation_metrics: Vec<ValidationMetric>,
        posterior_checks: Vec<PosteriorPredictiveCheck>,
        applicability_point: ApplicabilityPoint,
        applicability: ApplicabilityDecision,
        evidence_axes: EvidenceAxes,
        assumption_checks: Vec<(AssumptionId, bool)>,
    ) -> Result<Self, VvErrors> {
        if context.kind != ArtifactKind::ContextOfUse
            || validation_plan.kind != ArtifactKind::ValidationPlan
        {
            return Err(invalid(
                VvRule::QoiDependencyClosed,
                Some(header.id().as_str()),
                Some(qoi.as_str()),
                "prediction.context",
                "prediction must reference ContextOfUse and ValidationPlan",
            ));
        }
        if dependencies.is_empty() || dependencies.len() > MAX_VV_ITEMS {
            return Err(invalid(
                VvRule::QoiDependencyClosed,
                Some(header.id().as_str()),
                Some(qoi.as_str()),
                "prediction.dependencies",
                "prediction needs a bounded non-empty dependency closure",
            ));
        }
        if waterfall.qoi != qoi
            || validation_metrics.iter().any(|metric| metric.qoi != qoi)
            || posterior_checks.iter().any(|check| check.qoi != qoi)
            || dependencies.iter().any(|dependency| dependency.qoi != qoi)
        {
            return Err(invalid(
                VvRule::QoiDependencyIsolated,
                Some(header.id().as_str()),
                Some(qoi.as_str()),
                "prediction.qoi",
                "all prediction evidence must name exactly the enclosing QoI",
            ));
        }
        let dependency_count = dependencies.len();
        let mut dependencies = dependencies;
        dependencies.sort();
        dependencies.dedup();
        if dependencies.len() != dependency_count {
            return Err(invalid(
                VvRule::QoiDependencyClosed,
                Some(header.id().as_str()),
                Some(qoi.as_str()),
                "prediction.dependencies",
                "exact dependency edges must be unique",
            ));
        }
        let metric_count = validation_metrics.len();
        let mut validation_metrics = validation_metrics;
        validation_metrics.sort_by(|left, right| left.name.cmp(&right.name));
        validation_metrics.dedup_by(|left, right| left.name == right.name);
        if validation_metrics.len() != metric_count {
            return Err(invalid(
                VvRule::ValidationMetricUncertainty,
                Some(header.id().as_str()),
                Some(qoi.as_str()),
                "prediction.validation_metrics",
                "validation metric identities must be unique",
            ));
        }
        let posterior_count = posterior_checks.len();
        let mut posterior_checks = posterior_checks;
        posterior_checks.sort_by(|left, right| left.name.cmp(&right.name));
        posterior_checks.dedup_by(|left, right| left.name == right.name);
        if posterior_checks.len() != posterior_count {
            return Err(invalid(
                VvRule::ValidationMetricUncertainty,
                Some(header.id().as_str()),
                Some(qoi.as_str()),
                "prediction.posterior_checks",
                "posterior-predictive check identities must be unique",
            ));
        }
        let mut checks = BTreeMap::new();
        for (id, passed) in assumption_checks {
            if checks.insert(id, passed).is_some() {
                return Err(invalid(
                    VvRule::AssumptionRowComplete,
                    Some(header.id().as_str()),
                    Some(qoi.as_str()),
                    "prediction.assumption_checks",
                    "assumption checks must be unique",
                ));
            }
        }
        Ok(Self {
            header,
            context,
            validation_plan,
            qoi,
            dependencies,
            waterfall,
            validation_metrics,
            posterior_checks,
            applicability_point,
            applicability,
            evidence_axes,
            assumption_checks: checks,
        })
    }

    #[must_use]
    pub const fn id(&self) -> &ArtifactId {
        self.header.id()
    }

    #[must_use]
    pub const fn header(&self) -> &ArtifactHeader {
        &self.header
    }

    #[must_use]
    pub const fn context(&self) -> &ArtifactRef {
        &self.context
    }

    #[must_use]
    pub const fn validation_plan(&self) -> &ArtifactRef {
        &self.validation_plan
    }

    #[must_use]
    pub const fn qoi(&self) -> &QoiId {
        &self.qoi
    }

    #[must_use]
    pub fn dependencies(&self) -> &[EvidenceDependency] {
        &self.dependencies
    }

    #[must_use]
    pub const fn waterfall(&self) -> &UncertaintyWaterfall {
        &self.waterfall
    }

    #[must_use]
    pub fn validation_metrics(&self) -> &[ValidationMetric] {
        &self.validation_metrics
    }

    #[must_use]
    pub fn posterior_checks(&self) -> &[PosteriorPredictiveCheck] {
        &self.posterior_checks
    }

    #[must_use]
    pub const fn applicability_point(&self) -> &ApplicabilityPoint {
        &self.applicability_point
    }

    #[must_use]
    pub const fn applicability(&self) -> &ApplicabilityDecision {
        &self.applicability
    }

    #[must_use]
    pub const fn evidence_axes(&self) -> &EvidenceAxes {
        &self.evidence_axes
    }

    #[must_use]
    pub const fn assumption_checks(&self) -> &BTreeMap<AssumptionId, bool> {
        &self.assumption_checks
    }
}

/// Evidence required to support one runtime assumption.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AssumptionEvidence {
    requirement: String,
    artifact: Option<EvidenceTarget>,
}

impl AssumptionEvidence {
    pub fn try_new(
        requirement: impl Into<String>,
        artifact: Option<EvidenceTarget>,
    ) -> Result<Self, VvErrors> {
        let requirement = requirement.into();
        validate_text(&requirement, "assumption.evidence.requirement")?;
        Ok(Self {
            requirement,
            artifact,
        })
    }

    #[must_use]
    pub fn requirement(&self) -> &str {
        &self.requirement
    }

    #[must_use]
    pub const fn artifact(&self) -> Option<&EvidenceTarget> {
        self.artifact.as_ref()
    }

    #[must_use]
    pub fn with_artifact(mut self, artifact: EvidenceTarget) -> Self {
        self.artifact = Some(artifact);
        self
    }
}

/// Runtime signal that makes one assumption falsifiable during execution.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RuntimeMonitorSpec {
    signal: String,
    evidence_hash: Option<ContentHash>,
}

impl RuntimeMonitorSpec {
    pub fn try_new(
        signal: impl Into<String>,
        evidence_hash: Option<ContentHash>,
    ) -> Result<Self, VvErrors> {
        let signal = signal.into();
        validate_text(&signal, "assumption.monitor.signal")?;
        Ok(Self {
            signal,
            evidence_hash,
        })
    }

    #[must_use]
    pub fn signal(&self) -> &str {
        &self.signal
    }

    #[must_use]
    pub const fn evidence_hash(&self) -> Option<ContentHash> {
        self.evidence_hash
    }

    #[must_use]
    pub fn with_evidence(mut self, evidence_hash: ContentHash) -> Self {
        self.evidence_hash = Some(evidence_hash);
        self
    }
}

/// Required response when an assumption is false or cannot be monitored.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ViolationEffect {
    Demote { reason: String },
    EscalateOrRefuse { target_lane: ArtifactId },
    Refuse { reason: String },
}

impl ViolationEffect {
    fn validate(&self) -> Result<(), VvErrors> {
        match self {
            Self::Demote { reason } | Self::Refuse { reason } => {
                validate_text(reason, "assumption.violation_effect")
            }
            Self::EscalateOrRefuse { .. } => Ok(()),
        }
    }
}

/// Expiry or mandatory review cadence for an assumption row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ReviewGate {
    Phase { gate: ArtifactId },
    EverySolve,
    EveryQuery,
    EveryUpdate,
}

/// Complete operational row in the runtime-physics assumptions ledger.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AssumptionRow {
    id: AssumptionId,
    predicate: String,
    scope: String,
    evidence: AssumptionEvidence,
    monitor: RuntimeMonitorSpec,
    violation_effect: ViolationEffect,
    owner: ArtifactId,
    review_gate: ReviewGate,
}

impl AssumptionRow {
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        id: AssumptionId,
        predicate: impl Into<String>,
        scope: impl Into<String>,
        evidence: AssumptionEvidence,
        monitor: RuntimeMonitorSpec,
        violation_effect: ViolationEffect,
        owner: ArtifactId,
        review_gate: ReviewGate,
    ) -> Result<Self, VvErrors> {
        let predicate = predicate.into();
        let scope = scope.into();
        validate_text(&predicate, "assumption.predicate")?;
        validate_text(&scope, "assumption.scope")?;
        violation_effect.validate()?;
        Ok(Self {
            id,
            predicate,
            scope,
            evidence,
            monitor,
            violation_effect,
            owner,
            review_gate,
        })
    }

    #[must_use]
    pub const fn id(&self) -> &AssumptionId {
        &self.id
    }

    #[must_use]
    pub fn predicate(&self) -> &str {
        &self.predicate
    }

    #[must_use]
    pub fn scope(&self) -> &str {
        &self.scope
    }

    #[must_use]
    pub const fn evidence(&self) -> &AssumptionEvidence {
        &self.evidence
    }

    #[must_use]
    pub const fn monitor(&self) -> &RuntimeMonitorSpec {
        &self.monitor
    }

    #[must_use]
    pub const fn violation_effect(&self) -> &ViolationEffect {
        &self.violation_effect
    }

    #[must_use]
    pub const fn owner(&self) -> &ArtifactId {
        &self.owner
    }

    #[must_use]
    pub const fn review_gate(&self) -> &ReviewGate {
        &self.review_gate
    }

    #[must_use]
    pub fn with_evidence(mut self, artifact: EvidenceTarget) -> Self {
        self.evidence = self.evidence.with_artifact(artifact);
        self
    }

    #[must_use]
    pub fn with_monitor_evidence(mut self, evidence_hash: ContentHash) -> Self {
        self.monitor = self.monitor.with_evidence(evidence_hash);
        self
    }

    fn has_same_seed_semantics(&self, expected: &Self) -> bool {
        self.id == expected.id
            && self.predicate == expected.predicate
            && self.scope == expected.scope
            && self.evidence.requirement == expected.evidence.requirement
            && self.monitor.signal == expected.monitor.signal
            && self.violation_effect == expected.violation_effect
            && self.owner == expected.owner
            && self.review_gate == expected.review_gate
    }
}

fn seed_assumption_row(
    id: &str,
    predicate: &str,
    scope: &str,
    evidence: &str,
    monitor: &str,
    violation_effect: ViolationEffect,
    owner: &str,
    review_gate: ReviewGate,
) -> Result<AssumptionRow, VvErrors> {
    AssumptionRow::try_new(
        AssumptionId::try_new(id)?,
        predicate,
        scope,
        AssumptionEvidence::try_new(evidence, None)?,
        RuntimeMonitorSpec::try_new(monitor, None)?,
        violation_effect,
        ArtifactId::try_new(owner)?,
        review_gate,
    )
}

fn program_seed_rows() -> Result<Vec<AssumptionRow>, VvErrors> {
    Ok(vec![
        seed_assumption_row(
            "A-001",
            "The rigid or reduced-body rung is adequate for the named quantities of interest.",
            "Rigid-body and reduced-order multibody dynamics.",
            "Mode-truncation and interface-compliance evidence.",
            "Mode-truncation and interface-compliance indicators.",
            ViolationEffect::EscalateOrRefuse {
                target_lane: ArtifactId::try_new("fs-mbd:flexible-rung")?,
            },
            "fs-mbd",
            ReviewGate::Phase {
                gate: ArtifactId::try_new("E3")?,
            },
        )?,
        seed_assumption_row(
            "A-002",
            "The magnetoquasistatic wavelength and displacement-current regime is valid.",
            "Electromagnetic solves that select the magnetoquasistatic rung.",
            "Dimensionless electromagnetic-regime receipt.",
            "Frequency, wavelength, and displacement-current regime indicators.",
            ViolationEffect::EscalateOrRefuse {
                target_lane: ArtifactId::try_new("fs-em:electroquasistatic-or-full-wave")?,
            },
            "fs-em",
            ReviewGate::Phase {
                gate: ArtifactId::try_new("E4-E6")?,
            },
        )?,
        seed_assumption_row(
            "A-003",
            "Zero-dimensional cylinders are spatially mixed and one-dimensional ducts are section averaged.",
            "Gas-network cylinder and duct representations.",
            "Stratification, knock, wave, and secondary-flow adequacy evidence.",
            "Stratification, knock, wave, and secondary-flow indicators.",
            ViolationEffect::EscalateOrRefuse {
                target_lane: ArtifactId::try_new("fs-gas:multi-zone-or-multi-dimensional")?,
            },
            "fs-gas",
            ReviewGate::Phase {
                gate: ArtifactId::try_new("E5")?,
            },
        )?,
        seed_assumption_row(
            "A-004",
            "The smooth-contact law and continuum-scale representation are adequate.",
            "Contact, thin-feature, roughness, rarefaction, and short-duration impact regimes.",
            "Curvature, thickness, roughness, Knudsen-number, and contact-duration evidence.",
            "Curvature, thickness, roughness, Knudsen-number, and contact-duration indicators.",
            ViolationEffect::EscalateOrRefuse {
                target_lane: ArtifactId::try_new(
                    "fs-contact+rep-router:alternate-law-or-representation",
                )?,
            },
            "fs-contact+rep-router",
            ReviewGate::Phase {
                gate: ArtifactId::try_new("E2-E6")?,
            },
        )?,
        seed_assumption_row(
            "A-005",
            "The declared symmetry is preserved by geometry, loads, boundary conditions, and the evolving state.",
            "Any solve using symmetry reduction.",
            "Symmetry-residual evidence for the full reduced boundary.",
            "Symmetry residual at every state update.",
            ViolationEffect::EscalateOrRefuse {
                target_lane: ArtifactId::try_new("full-model")?,
            },
            "domain-owner",
            ReviewGate::EverySolve,
        )?,
        seed_assumption_row(
            "A-006",
            "Material, process, and query validity predicates hold for the consuming calculation.",
            "Every property query and downstream material use.",
            "A property-usage receipt with source, process, and validity-domain lineage.",
            "PropertyUsageReceipt monitor at every query.",
            ViolationEffect::Demote {
                reason: "Demote, recalibrate, or refuse according to the property policy."
                    .to_owned(),
            },
            "fs-matdb-consumer",
            ReviewGate::EveryQuery,
        )?,
        seed_assumption_row(
            "A-007",
            "The selected closure, correlation, turbulence, or lubrication law is applicable.",
            "Every empirical or reduced constitutive-law invocation.",
            "Dimensionless-group and held-out-discrepancy evidence.",
            "Dimensionless groups and held-out discrepancy at every solve.",
            ViolationEffect::EscalateOrRefuse {
                target_lane: ArtifactId::try_new("resolved-or-high-fidelity-law")?,
            },
            "law-owner",
            ReviewGate::Phase {
                gate: ArtifactId::try_new("phase-gate")?,
            },
        )?,
        seed_assumption_row(
            "A-008",
            "The probability model, dependence structure, and sampled population are representative for the decision.",
            "Uncertainty, reliability, posterior, and population-level claims.",
            "Sampling lineage, posterior-predictive, dependence, and population evidence.",
            "Population drift, posterior-predictive, and dependence diagnostics.",
            ViolationEffect::Demote {
                reason: "Demote the probabilistic claim or use a robust worst-case decision."
                    .to_owned(),
            },
            "fs-uq",
            ReviewGate::Phase {
                gate: ArtifactId::try_new("E7")?,
            },
        )?,
    ])
}

fn assumption_rule(id: &AssumptionId) -> VvRule {
    match id.as_str() {
        "A-001" => VvRule::AssumptionA001,
        "A-002" => VvRule::AssumptionA002,
        "A-003" => VvRule::AssumptionA003,
        "A-004" => VvRule::AssumptionA004,
        "A-005" => VvRule::AssumptionA005,
        "A-006" => VvRule::AssumptionA006,
        "A-007" => VvRule::AssumptionA007,
        "A-008" => VvRule::AssumptionA008,
        _ => VvRule::AssumptionRowComplete,
    }
}

/// Runtime-physics assumptions and their falsification/escalation policies.
#[derive(Debug, Clone, PartialEq)]
pub struct AssumptionsLedger {
    header: ArtifactHeader,
    rows: BTreeMap<AssumptionId, AssumptionRow>,
}

impl AssumptionsLedger {
    pub fn try_new(header: ArtifactHeader, rows: Vec<AssumptionRow>) -> Result<Self, VvErrors> {
        if rows.is_empty() || rows.len() > MAX_VV_ITEMS {
            return Err(invalid(
                VvRule::AssumptionRowComplete,
                Some(header.id().as_str()),
                None,
                "assumptions.rows",
                "the assumptions ledger needs a bounded non-empty row set",
            ));
        }
        let mut by_id = BTreeMap::new();
        for row in rows {
            let id = row.id.clone();
            if by_id.insert(id.clone(), row).is_some() {
                return Err(invalid(
                    VvRule::AssumptionRowComplete,
                    Some(header.id().as_str()),
                    None,
                    "assumptions.rows",
                    format!("duplicate assumption row {id}"),
                ));
            }
        }
        Ok(Self {
            header,
            rows: by_id,
        })
    }

    pub fn try_program_seed(header: ArtifactHeader) -> Result<Self, VvErrors> {
        Self::try_new(header, program_seed_rows()?)
    }

    #[must_use]
    pub const fn id(&self) -> &ArtifactId {
        self.header.id()
    }

    #[must_use]
    pub const fn header(&self) -> &ArtifactHeader {
        &self.header
    }

    #[must_use]
    pub const fn rows(&self) -> &BTreeMap<AssumptionId, AssumptionRow> {
        &self.rows
    }

    pub fn replace_row(&mut self, row: AssumptionRow) -> Result<(), VvErrors> {
        if !self.rows.contains_key(row.id()) {
            return Err(invalid(
                VvRule::AssumptionRowComplete,
                Some(self.id().as_str()),
                None,
                "assumptions.rows",
                "new assumption identities require a new ledger artifact",
            ));
        }
        self.rows.insert(row.id.clone(), row);
        Ok(())
    }

    fn seed_violations(&self) -> Vec<VvViolation> {
        let Ok(expected_rows) = program_seed_rows() else {
            return vec![VvViolation::new(
                VvRule::AssumptionRowComplete,
                Some(self.id().as_str().to_owned()),
                None,
                "assumptions.program_seed",
                "the built-in seed rows could not be constructed",
            )];
        };
        expected_rows
            .into_iter()
            .filter_map(|expected| match self.rows.get(expected.id()) {
                Some(actual) if actual.has_same_seed_semantics(&expected) => None,
                Some(_) => Some(VvViolation::new(
                    assumption_rule(expected.id()),
                    Some(self.id().as_str().to_owned()),
                    None,
                    "assumptions.program_seed",
                    format!(
                        "{} does not preserve the required seed semantics",
                        expected.id()
                    ),
                )),
                None => Some(VvViolation::new(
                    assumption_rule(expected.id()),
                    Some(self.id().as_str().to_owned()),
                    None,
                    "assumptions.program_seed",
                    format!("required program assumption {} is missing", expected.id()),
                )),
            })
            .collect()
    }
}

/// Any one of the seven top-level artifact schemas.
#[derive(Debug, Clone, PartialEq)]
pub enum VvArtifact {
    ContextOfUse(ContextOfUse),
    ValidationPlan(ValidationPlan),
    ExperimentArtifact(ExperimentArtifact),
    CalibrationSplit(CalibrationSplit),
    SolutionVerificationReceipt(SolutionVerificationReceipt),
    PredictionAssessment(PredictionAssessment),
    AssumptionsLedger(AssumptionsLedger),
}

impl VvArtifact {
    #[must_use]
    pub const fn kind(&self) -> ArtifactKind {
        match self {
            Self::ContextOfUse(_) => ArtifactKind::ContextOfUse,
            Self::ValidationPlan(_) => ArtifactKind::ValidationPlan,
            Self::ExperimentArtifact(_) => ArtifactKind::ExperimentArtifact,
            Self::CalibrationSplit(_) => ArtifactKind::CalibrationSplit,
            Self::SolutionVerificationReceipt(_) => ArtifactKind::SolutionVerificationReceipt,
            Self::PredictionAssessment(_) => ArtifactKind::PredictionAssessment,
            Self::AssumptionsLedger(_) => ArtifactKind::AssumptionsLedger,
        }
    }

    #[must_use]
    pub const fn header(&self) -> &ArtifactHeader {
        match self {
            Self::ContextOfUse(value) => value.header(),
            Self::ValidationPlan(value) => value.header(),
            Self::ExperimentArtifact(value) => value.header(),
            Self::CalibrationSplit(value) => value.header(),
            Self::SolutionVerificationReceipt(value) => value.header(),
            Self::PredictionAssessment(value) => value.header(),
            Self::AssumptionsLedger(value) => value.header(),
        }
    }

    #[must_use]
    pub const fn id(&self) -> &ArtifactId {
        self.header().id()
    }
}

macro_rules! vv_artifact_from {
    ($type:ty, $variant:ident) => {
        impl From<$type> for VvArtifact {
            fn from(value: $type) -> Self {
                Self::$variant(value)
            }
        }
    };
}

vv_artifact_from!(ContextOfUse, ContextOfUse);
vv_artifact_from!(ValidationPlan, ValidationPlan);
vv_artifact_from!(ExperimentArtifact, ExperimentArtifact);
vv_artifact_from!(CalibrationSplit, CalibrationSplit);
vv_artifact_from!(SolutionVerificationReceipt, SolutionVerificationReceipt);
vv_artifact_from!(PredictionAssessment, PredictionAssessment);
vv_artifact_from!(AssumptionsLedger, AssumptionsLedger);

/// A complete, closed V&V case before structural admission.
#[derive(Debug, Clone, PartialEq)]
pub struct VvCase {
    context: ContextOfUse,
    validation_plan: ValidationPlan,
    experiments: BTreeMap<ArtifactId, ExperimentArtifact>,
    splits: BTreeMap<ArtifactId, CalibrationSplit>,
    solution_verification: BTreeMap<ArtifactId, SolutionVerificationReceipt>,
    predictions: BTreeMap<ArtifactId, PredictionAssessment>,
    assumptions: AssumptionsLedger,
}

impl VvCase {
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        context: ContextOfUse,
        validation_plan: ValidationPlan,
        experiments: Vec<ExperimentArtifact>,
        splits: Vec<CalibrationSplit>,
        solution_verification: Vec<SolutionVerificationReceipt>,
        predictions: Vec<PredictionAssessment>,
        assumptions: AssumptionsLedger,
    ) -> Result<Self, VvErrors> {
        let total = 3usize
            .saturating_add(experiments.len())
            .saturating_add(splits.len())
            .saturating_add(solution_verification.len())
            .saturating_add(predictions.len());
        if total > MAX_VV_ITEMS {
            return Err(invalid(
                VvRule::SchemaCardinality,
                Some(context.id().as_str()),
                None,
                "case.artifacts",
                "the complete V&V case exceeds the artifact bound",
            ));
        }

        let mut experiment_map = BTreeMap::new();
        for artifact in experiments {
            let id = artifact.id().clone();
            if experiment_map.insert(id.clone(), artifact).is_some() {
                return Err(duplicate_artifact(ArtifactKind::ExperimentArtifact, &id));
            }
        }
        let mut split_map = BTreeMap::new();
        for artifact in splits {
            let id = artifact.id().clone();
            if split_map.insert(id.clone(), artifact).is_some() {
                return Err(duplicate_artifact(ArtifactKind::CalibrationSplit, &id));
            }
        }
        let mut solution_map = BTreeMap::new();
        for artifact in solution_verification {
            let id = artifact.id().clone();
            if solution_map.insert(id.clone(), artifact).is_some() {
                return Err(duplicate_artifact(
                    ArtifactKind::SolutionVerificationReceipt,
                    &id,
                ));
            }
        }
        let mut prediction_map = BTreeMap::new();
        for artifact in predictions {
            let id = artifact.id().clone();
            if prediction_map.insert(id.clone(), artifact).is_some() {
                return Err(duplicate_artifact(ArtifactKind::PredictionAssessment, &id));
            }
        }
        Ok(Self {
            context,
            validation_plan,
            experiments: experiment_map,
            splits: split_map,
            solution_verification: solution_map,
            predictions: prediction_map,
            assumptions,
        })
    }

    #[must_use]
    pub const fn context(&self) -> &ContextOfUse {
        &self.context
    }

    #[must_use]
    pub const fn validation_plan(&self) -> &ValidationPlan {
        &self.validation_plan
    }

    #[must_use]
    pub const fn experiments(&self) -> &BTreeMap<ArtifactId, ExperimentArtifact> {
        &self.experiments
    }

    #[must_use]
    pub const fn splits(&self) -> &BTreeMap<ArtifactId, CalibrationSplit> {
        &self.splits
    }

    #[must_use]
    pub const fn solution_verification(
        &self,
    ) -> &BTreeMap<ArtifactId, SolutionVerificationReceipt> {
        &self.solution_verification
    }

    #[must_use]
    pub const fn predictions(&self) -> &BTreeMap<ArtifactId, PredictionAssessment> {
        &self.predictions
    }

    #[must_use]
    pub const fn assumptions(&self) -> &AssumptionsLedger {
        &self.assumptions
    }

    #[must_use]
    pub fn artifacts(&self) -> Vec<VvArtifact> {
        let mut artifacts = Vec::with_capacity(
            3 + self.experiments.len()
                + self.splits.len()
                + self.solution_verification.len()
                + self.predictions.len(),
        );
        artifacts.push(VvArtifact::ContextOfUse(self.context.clone()));
        artifacts.push(VvArtifact::ValidationPlan(self.validation_plan.clone()));
        artifacts.extend(
            self.experiments
                .values()
                .cloned()
                .map(VvArtifact::ExperimentArtifact),
        );
        artifacts.extend(
            self.splits
                .values()
                .cloned()
                .map(VvArtifact::CalibrationSplit),
        );
        artifacts.extend(
            self.solution_verification
                .values()
                .cloned()
                .map(VvArtifact::SolutionVerificationReceipt),
        );
        artifacts.extend(
            self.predictions
                .values()
                .cloned()
                .map(VvArtifact::PredictionAssessment),
        );
        artifacts.push(VvArtifact::AssumptionsLedger(self.assumptions.clone()));
        artifacts
    }
}

fn push_receipt_string(bytes: &mut Vec<u8>, value: &str) {
    bytes.extend_from_slice(&(value.len() as u64).to_le_bytes());
    bytes.extend_from_slice(value.as_bytes());
}

fn receipt_identity(
    schema_version: u32,
    ruleset_version: u32,
    case_hash: ContentHash,
    context_id: &ArtifactId,
    qois: &BTreeSet<QoiId>,
    artifact_hashes: &ArtifactHashMap,
) -> ContentHash {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(&schema_version.to_le_bytes());
    bytes.extend_from_slice(&ruleset_version.to_le_bytes());
    bytes.extend_from_slice(case_hash.as_bytes());
    push_receipt_string(&mut bytes, context_id.as_str());
    bytes.extend_from_slice(&(qois.len() as u64).to_le_bytes());
    for qoi in qois {
        push_receipt_string(&mut bytes, qoi.as_str());
    }
    bytes.extend_from_slice(&(artifact_hashes.len() as u64).to_le_bytes());
    for ((kind, id), hash) in artifact_hashes {
        push_receipt_string(&mut bytes, kind.slug());
        push_receipt_string(&mut bytes, id.as_str());
        bytes.extend_from_slice(hash.as_bytes());
    }
    fs_blake3::hash_domain(
        "org.frankensim.fs-evidence.vv-schema-admission-receipt.v1",
        &bytes,
    )
}

/// Content-bound proof that the current V&V structural rules admitted a case.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SchemaAdmissionReceipt {
    schema_version: u32,
    ruleset_version: u32,
    case_hash: ContentHash,
    context_id: ArtifactId,
    qois: BTreeSet<QoiId>,
    artifact_hashes: ArtifactHashMap,
    receipt_hash: ContentHash,
}

impl SchemaAdmissionReceipt {
    fn new(
        case_hash: ContentHash,
        context_id: ArtifactId,
        qois: BTreeSet<QoiId>,
        artifact_hashes: ArtifactHashMap,
    ) -> Self {
        let receipt_hash = receipt_identity(
            VV_SCHEMA_VERSION,
            VV_RULESET_VERSION,
            case_hash,
            &context_id,
            &qois,
            &artifact_hashes,
        );
        Self {
            schema_version: VV_SCHEMA_VERSION,
            ruleset_version: VV_RULESET_VERSION,
            case_hash,
            context_id,
            qois,
            artifact_hashes,
            receipt_hash,
        }
    }

    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    #[must_use]
    pub const fn ruleset_version(&self) -> u32 {
        self.ruleset_version
    }

    #[must_use]
    pub const fn case_hash(&self) -> ContentHash {
        self.case_hash
    }

    #[must_use]
    pub const fn context_id(&self) -> &ArtifactId {
        &self.context_id
    }

    #[must_use]
    pub const fn qois(&self) -> &BTreeSet<QoiId> {
        &self.qois
    }

    #[must_use]
    pub const fn artifact_hashes(&self) -> &ArtifactHashMap {
        &self.artifact_hashes
    }

    #[must_use]
    pub const fn receipt_hash(&self) -> ContentHash {
        self.receipt_hash
    }

    #[must_use]
    pub fn has_valid_binding(&self) -> bool {
        self.schema_version == VV_SCHEMA_VERSION
            && self.ruleset_version == VV_RULESET_VERSION
            && self.receipt_hash
                == receipt_identity(
                    self.schema_version,
                    self.ruleset_version,
                    self.case_hash,
                    &self.context_id,
                    &self.qois,
                    &self.artifact_hashes,
                )
    }

    /// Re-run admission and verify that this receipt binds the exact case.
    pub fn verify_case(&self, case: &VvCase) -> Result<(), VvErrors> {
        case.validate()?;
        let case_hash = case.content_hash().map_err(|error| {
            invalid(
                VvRule::ReceiptBinding,
                Some(case.context.id().as_str()),
                None,
                "receipt.case_hash",
                format!("canonical case identity failed: {error}"),
            )
        })?;
        let artifact_hashes = case.artifact_hashes()?;
        let qois = case.context.qois.keys().cloned().collect::<BTreeSet<_>>();
        if !self.has_valid_binding()
            || self.case_hash != case_hash
            || self.context_id != *case.context.id()
            || self.qois != qois
            || self.artifact_hashes != artifact_hashes
        {
            return Err(invalid(
                VvRule::ReceiptBinding,
                Some(case.context.id().as_str()),
                None,
                "receipt",
                "schema admission receipt does not bind the exact canonical case",
            ));
        }
        Ok(())
    }
}

/// Opaque positive result of exact schema admission.
#[derive(Debug, Clone, PartialEq)]
pub struct AdmittedVvCase {
    case: VvCase,
    receipt: SchemaAdmissionReceipt,
}

impl AdmittedVvCase {
    #[must_use]
    pub const fn case(&self) -> &VvCase {
        &self.case
    }

    #[must_use]
    pub const fn receipt(&self) -> &SchemaAdmissionReceipt {
        &self.receipt
    }

    #[must_use]
    pub fn into_parts(self) -> (VvCase, SchemaAdmissionReceipt) {
        (self.case, self.receipt)
    }
}

impl VvCase {
    /// Admit an exact closed case and mint a content-bound structural receipt.
    pub fn admit(self) -> Result<AdmittedVvCase, VvErrors> {
        self.validate()?;
        let artifact_hashes = self.artifact_hashes()?;
        let case_hash = self.content_hash().map_err(|error| {
            invalid(
                VvRule::ReceiptBinding,
                Some(self.context.id().as_str()),
                None,
                "receipt.case_hash",
                format!("canonical case identity failed: {error}"),
            )
        })?;
        let receipt = SchemaAdmissionReceipt::new(
            case_hash,
            self.context.id().clone(),
            self.context.qois.keys().cloned().collect(),
            artifact_hashes,
        );
        Ok(AdmittedVvCase {
            case: self,
            receipt,
        })
    }
}

fn duplicate_artifact(kind: ArtifactKind, id: &ArtifactId) -> VvErrors {
    invalid(
        VvRule::SchemaIdentity,
        Some(id.as_str()),
        None,
        "case.artifacts",
        format!("duplicate {} artifact identity", kind.slug()),
    )
}

/// Canonical content identities keyed by artifact family and identity.
pub type ArtifactHashMap = BTreeMap<(ArtifactKind, ArtifactId), ContentHash>;

fn reference_violation(
    reference: &ArtifactRef,
    hashes: &ArtifactHashMap,
    artifact: &ArtifactId,
    qoi: Option<&QoiId>,
    field: &'static str,
    rule: VvRule,
) -> Option<VvViolation> {
    match hashes.get(&(reference.kind, reference.id.clone())) {
        Some(hash) if *hash == reference.hash => None,
        Some(_) => Some(VvViolation::new(
            rule,
            Some(artifact.as_str().to_owned()),
            qoi.map(|value| value.as_str().to_owned()),
            field,
            format!(
                "{} reference {} has a stale or mismatched content hash",
                reference.kind.slug(),
                reference.id
            ),
        )),
        None => Some(VvViolation::new(
            rule,
            Some(artifact.as_str().to_owned()),
            qoi.map(|value| value.as_str().to_owned()),
            field,
            format!(
                "{} reference {} is absent from the case",
                reference.kind.slug(),
                reference.id
            ),
        )),
    }
}

fn target_reference_violation(
    target: &EvidenceTarget,
    hashes: &ArtifactHashMap,
    artifact: &ArtifactId,
    qoi: &QoiId,
    field: &'static str,
) -> Option<VvViolation> {
    match target {
        EvidenceTarget::VvArtifact(reference) => reference_violation(
            reference,
            hashes,
            artifact,
            Some(qoi),
            field,
            VvRule::QoiDependencyClosed,
        ),
        EvidenceTarget::External { .. } => None,
    }
}

fn selection_violations(
    selection: &ObservationSelection,
    splits: &BTreeMap<ArtifactId, CalibrationSplit>,
    hashes: &ArtifactHashMap,
    artifact: &ArtifactId,
    qoi: &QoiId,
    field: &'static str,
) -> Vec<VvViolation> {
    let mut violations = Vec::new();
    if let Some(violation) = reference_violation(
        selection.split(),
        hashes,
        artifact,
        Some(qoi),
        field,
        VvRule::SplitPartitionsDisjoint,
    ) {
        violations.push(violation);
        return violations;
    }
    let Some(split) = splits.get(selection.split().id()) else {
        return violations;
    };
    let allowed = match selection.partition() {
        EvidencePartition::Validation => Some(&split.validation),
        EvidencePartition::BlindHoldout { release }
            if release.split == selection.split
                && release.blind_commitment == split.blind_commitment
                && release
                    .authority_receipt_hash
                    .as_bytes()
                    .iter()
                    .any(|byte| *byte != 0) =>
        {
            Some(&split.blind_holdout)
        }
        EvidencePartition::BlindHoldout { .. } => {
            violations.push(VvViolation::new(
                VvRule::SplitBlindHoldoutSealed,
                Some(artifact.as_str().to_owned()),
                Some(qoi.as_str().to_owned()),
                field,
                "blind selection is not bound to this split and preregistered commitment",
            ));
            None
        }
    };
    if let Some(allowed) = allowed {
        if !selection.ids.is_subset(allowed) {
            let calibration_reuse = selection
                .ids
                .iter()
                .any(|id| split.calibration.contains(id));
            violations.push(VvViolation::new(
                if calibration_reuse {
                    VvRule::ValidationRequiresPhysicalReferent
                } else {
                    VvRule::SplitPartitionsDisjoint
                },
                Some(artifact.as_str().to_owned()),
                Some(qoi.as_str().to_owned()),
                field,
                if calibration_reuse {
                    "calibration observations cannot support validation"
                } else {
                    "selection contains observations outside its declared partition"
                },
            ));
        }
    }
    violations
}

impl VvCase {
    fn artifact_hashes(&self) -> Result<ArtifactHashMap, VvErrors> {
        let mut hashes = BTreeMap::new();
        let artifacts = self.artifacts();
        let mut ids = BTreeMap::<ArtifactId, ArtifactKind>::new();
        for artifact in &artifacts {
            if let Some(existing) = ids.insert(artifact.id().clone(), artifact.kind()) {
                return Err(invalid(
                    VvRule::SchemaIdentity,
                    Some(artifact.id().as_str()),
                    None,
                    "case.artifacts",
                    format!(
                        "artifact identity is shared by {} and {}",
                        existing.slug(),
                        artifact.kind().slug()
                    ),
                ));
            }
            let hash = artifact.content_hash().map_err(|error| {
                invalid(
                    VvRule::SchemaIdentity,
                    Some(artifact.id().as_str()),
                    None,
                    "case.artifact_hash",
                    format!("canonical artifact identity failed: {error}"),
                )
            })?;
            hashes.insert((artifact.kind(), artifact.id().clone()), hash);
        }
        Ok(hashes)
    }

    /// Run all cross-artifact structural rules without minting a receipt.
    pub fn validate(&self) -> Result<(), VvErrors> {
        let hashes = self.artifact_hashes()?;
        let mut violations = self.assumptions.seed_violations();
        self.validate_context_and_plan(&hashes, &mut violations);
        self.validate_experiments_and_splits(&hashes, &mut violations);
        self.validate_solutions(&mut violations);
        self.validate_predictions(&hashes, &mut violations);
        self.validate_assumption_evidence(&hashes, &mut violations);
        if violations.is_empty() {
            Ok(())
        } else {
            Err(VvErrors::from_vec(violations))
        }
    }

    fn validate_context_and_plan(
        &self,
        hashes: &ArtifactHashMap,
        violations: &mut Vec<VvViolation>,
    ) {
        let context_qois = self.context.qois.keys().cloned().collect::<BTreeSet<_>>();
        let plan_qois = self
            .validation_plan
            .by_qoi
            .keys()
            .cloned()
            .collect::<BTreeSet<_>>();
        if context_qois != plan_qois {
            violations.push(VvViolation::new(
                VvRule::QoiDependencyClosed,
                Some(self.validation_plan.id().as_str().to_owned()),
                None,
                "validation_plan.by_qoi",
                "validation-plan QoIs must exactly equal the ContextOfUse QoIs",
            ));
        }
        if let Some(violation) = reference_violation(
            self.validation_plan.context(),
            hashes,
            self.validation_plan.id(),
            None,
            "validation_plan.context",
            VvRule::QoiDependencyClosed,
        ) {
            violations.push(violation);
        }
        if self.validation_plan.context.id != *self.context.id() {
            violations.push(VvViolation::new(
                VvRule::QoiDependencyClosed,
                Some(self.validation_plan.id().as_str().to_owned()),
                None,
                "validation_plan.context",
                "validation plan must reference the case ContextOfUse",
            ));
        }
        for qoi in self.context.qois.values() {
            if !self.context.header.units.contains(qoi.unit()) {
                violations.push(VvViolation::new(
                    VvRule::SchemaIdentity,
                    Some(self.context.id().as_str().to_owned()),
                    Some(qoi.id().as_str().to_owned()),
                    "context.header.units",
                    "every QoI unit must be declared in the ContextOfUse header",
                ));
            }
            if !self.validation_plan.header.units.contains(qoi.unit()) {
                violations.push(VvViolation::new(
                    VvRule::SchemaIdentity,
                    Some(self.validation_plan.id().as_str().to_owned()),
                    Some(qoi.id().as_str().to_owned()),
                    "validation_plan.header.units",
                    "every planned QoI unit must be declared in the ValidationPlan header",
                ));
            }
        }
        for row in self.validation_plan.by_qoi.values() {
            violations.extend(
                row.diagnostics
                    .violations(self.validation_plan.id(), row.qoi()),
            );
            for experiment_ref in &row.experiments {
                if let Some(violation) = reference_violation(
                    experiment_ref,
                    hashes,
                    self.validation_plan.id(),
                    Some(row.qoi()),
                    "validation_plan.experiments",
                    VvRule::ValidationRequiresPhysicalReferent,
                ) {
                    violations.push(violation);
                    continue;
                }
                if let Some(experiment) = self.experiments.get(experiment_ref.id()) {
                    if !experiment.origin.is_physical() || !experiment.qois.contains(row.qoi()) {
                        violations.push(VvViolation::new(
                            VvRule::ValidationRequiresPhysicalReferent,
                            Some(self.validation_plan.id().as_str().to_owned()),
                            Some(row.qoi().as_str().to_owned()),
                            "validation_plan.experiments",
                            "each validation referent must be a physical experiment for this exact QoI",
                        ));
                    }
                }
            }
            if let Some(violation) = reference_violation(
                row.split(),
                hashes,
                self.validation_plan.id(),
                Some(row.qoi()),
                "validation_plan.split",
                VvRule::SplitPartitionsDisjoint,
            ) {
                violations.push(violation);
            } else if let Some(split) = self.splits.get(row.split().id()) {
                if !row.experiments.contains(split.experiment()) {
                    violations.push(VvViolation::new(
                        VvRule::ValidationRequiresPhysicalReferent,
                        Some(self.validation_plan.id().as_str().to_owned()),
                        Some(row.qoi().as_str().to_owned()),
                        "validation_plan.split",
                        "the validation split must belong to a declared physical experiment",
                    ));
                }
            }
        }
    }

    fn validate_experiments_and_splits(
        &self,
        hashes: &ArtifactHashMap,
        violations: &mut Vec<VvViolation>,
    ) {
        for experiment in self.experiments.values() {
            let instrument_ids = experiment
                .instruments
                .iter()
                .map(|instrument| instrument.instrument_id())
                .collect::<BTreeSet<_>>();
            if instrument_ids.len() != experiment.instruments.len()
                || experiment.instruments.iter().any(|row| !row.current)
            {
                violations.push(VvViolation::new(
                    VvRule::ExperimentInstrumentCalibration,
                    Some(experiment.id().as_str().to_owned()),
                    None,
                    "experiment.instruments",
                    "instrument identities must be unique and every calibration current",
                ));
            }
            if experiment.clocks.clone().validated_canonical().is_err() {
                violations.push(VvViolation::new(
                    VvRule::ExperimentClockSynchronization,
                    Some(experiment.id().as_str().to_owned()),
                    None,
                    "experiment.clocks",
                    "clock topology or synchronization evidence is invalid",
                ));
            }
            for qoi in &experiment.qois {
                match self.context.qois.get(qoi) {
                    Some(spec) if experiment.header.units.contains(spec.unit()) => {}
                    Some(_) => violations.push(VvViolation::new(
                        VvRule::SchemaIdentity,
                        Some(experiment.id().as_str().to_owned()),
                        Some(qoi.as_str().to_owned()),
                        "experiment.header.units",
                        "experiment header must declare every observed QoI unit",
                    )),
                    None => violations.push(VvViolation::new(
                        VvRule::QoiDependencyClosed,
                        Some(experiment.id().as_str().to_owned()),
                        Some(qoi.as_str().to_owned()),
                        "experiment.qois",
                        "experiment QoI is absent from the ContextOfUse",
                    )),
                }
            }
        }
        for split in self.splits.values() {
            if let Some(violation) = reference_violation(
                split.experiment(),
                hashes,
                split.id(),
                None,
                "split.experiment",
                VvRule::SplitPartitionsDisjoint,
            ) {
                violations.push(violation);
                continue;
            }
            if let Some(experiment) = self.experiments.get(split.experiment().id()) {
                if split.all_ids() != experiment.observation_ids {
                    violations.push(VvViolation::new(
                        VvRule::SplitPartitionsDisjoint,
                        Some(split.id().as_str().to_owned()),
                        None,
                        "split.partitions",
                        "the three partitions must cover the experiment observation identities exactly",
                    ));
                }
            }
        }
    }

    fn validate_solutions(&self, violations: &mut Vec<VvViolation>) {
        for receipt in self.solution_verification.values() {
            match self.context.qois.get(receipt.qoi()) {
                Some(qoi)
                    if qoi.unit() == receipt.unit()
                        && receipt.header.units.contains(receipt.unit()) => {}
                Some(_) => violations.push(VvViolation::new(
                    VvRule::SolutionVerificationComplete,
                    Some(receipt.id().as_str().to_owned()),
                    Some(receipt.qoi().as_str().to_owned()),
                    "solution_verification.unit",
                    "solution-verification units must equal the ContextOfUse QoI unit",
                )),
                None => violations.push(VvViolation::new(
                    VvRule::QoiDependencyClosed,
                    Some(receipt.id().as_str().to_owned()),
                    Some(receipt.qoi().as_str().to_owned()),
                    "solution_verification.qoi",
                    "solution-verification QoI is absent from the ContextOfUse",
                )),
            }
        }
    }

    fn validate_assumption_evidence(
        &self,
        hashes: &ArtifactHashMap,
        violations: &mut Vec<VvViolation>,
    ) {
        for row in self.assumptions.rows.values() {
            if let Some(EvidenceTarget::VvArtifact(reference)) = row.evidence.artifact() {
                if let Some(violation) = reference_violation(
                    reference,
                    hashes,
                    self.assumptions.id(),
                    None,
                    "assumptions.evidence",
                    assumption_rule(row.id()),
                ) {
                    violations.push(violation);
                }
            }
        }
    }
}

impl VvCase {
    fn validate_predictions(&self, hashes: &ArtifactHashMap, violations: &mut Vec<VvViolation>) {
        let mut by_qoi = BTreeMap::<QoiId, &PredictionAssessment>::new();
        for prediction in self.predictions.values() {
            if by_qoi.insert(prediction.qoi.clone(), prediction).is_some() {
                violations.push(VvViolation::new(
                    VvRule::QoiDependencyIsolated,
                    Some(prediction.id().as_str().to_owned()),
                    Some(prediction.qoi().as_str().to_owned()),
                    "prediction.qoi",
                    "each ContextOfUse QoI must have exactly one prediction assessment",
                ));
            }
        }
        for qoi in self.context.qois.keys() {
            if !by_qoi.contains_key(qoi) {
                violations.push(VvViolation::new(
                    VvRule::QoiDependencyClosed,
                    Some(self.context.id().as_str().to_owned()),
                    Some(qoi.as_str().to_owned()),
                    "case.predictions",
                    "ContextOfUse QoI has no prediction assessment",
                ));
            }
        }
        for prediction in self.predictions.values() {
            self.validate_prediction(prediction, hashes, violations);
        }
    }

    fn validate_prediction(
        &self,
        prediction: &PredictionAssessment,
        hashes: &ArtifactHashMap,
        violations: &mut Vec<VvViolation>,
    ) {
        let qoi_id = prediction.qoi();
        let Some(qoi) = self.context.qois.get(qoi_id) else {
            violations.push(VvViolation::new(
                VvRule::QoiDependencyClosed,
                Some(prediction.id().as_str().to_owned()),
                Some(qoi_id.as_str().to_owned()),
                "prediction.qoi",
                "prediction QoI is absent from the ContextOfUse",
            ));
            return;
        };
        for (reference, field) in [
            (prediction.context(), "prediction.context"),
            (prediction.validation_plan(), "prediction.validation_plan"),
        ] {
            if let Some(violation) = reference_violation(
                reference,
                hashes,
                prediction.id(),
                Some(qoi_id),
                field,
                VvRule::QoiDependencyClosed,
            ) {
                violations.push(violation);
            }
        }
        if prediction.context.id != *self.context.id()
            || prediction.validation_plan.id != *self.validation_plan.id()
        {
            violations.push(VvViolation::new(
                VvRule::QoiDependencyClosed,
                Some(prediction.id().as_str().to_owned()),
                Some(qoi_id.as_str().to_owned()),
                "prediction.context",
                "prediction must bind this case's exact context and validation plan",
            ));
        }
        if prediction.waterfall.unit != *qoi.unit() || !prediction.header.units.contains(qoi.unit())
        {
            violations.push(VvViolation::new(
                VvRule::SchemaIdentity,
                Some(prediction.id().as_str().to_owned()),
                Some(qoi_id.as_str().to_owned()),
                "prediction.unit",
                "prediction waterfall and header must declare the exact QoI unit",
            ));
        }
        let Some(plan_row) = self.validation_plan.by_qoi.get(qoi_id) else {
            return;
        };

        let mut physical_targets = BTreeSet::new();
        let mut solution_targets = BTreeSet::new();
        let mut process_targets = BTreeSet::new();
        let mut role_targets = BTreeMap::<DependencyRole, BTreeSet<ContentHash>>::new();
        let mut referenced_solutions = Vec::new();
        for dependency in &prediction.dependencies {
            if dependency.qoi != *qoi_id {
                violations.push(VvViolation::new(
                    VvRule::QoiDependencyIsolated,
                    Some(prediction.id().as_str().to_owned()),
                    Some(qoi_id.as_str().to_owned()),
                    "prediction.dependencies",
                    "dependency QoI differs from the enclosing prediction QoI",
                ));
            }
            if let Some(violation) = target_reference_violation(
                dependency.target(),
                hashes,
                prediction.id(),
                qoi_id,
                "prediction.dependencies",
            ) {
                violations.push(violation);
                continue;
            }
            let target_hash = dependency.target.hash();
            role_targets
                .entry(dependency.role)
                .or_default()
                .insert(target_hash);
            match dependency.role {
                DependencyRole::PhysicalValidation => {
                    physical_targets.insert(target_hash);
                    self.validate_physical_dependency(
                        prediction, dependency, plan_row, hashes, violations,
                    );
                }
                DependencyRole::SolutionVerification => {
                    solution_targets.insert(target_hash);
                    match dependency.target() {
                        EvidenceTarget::VvArtifact(reference)
                            if reference.kind == ArtifactKind::SolutionVerificationReceipt =>
                        {
                            if let Some(receipt) = self.solution_verification.get(reference.id()) {
                                if receipt.qoi() == qoi_id && receipt.unit() == qoi.unit() {
                                    referenced_solutions.push(receipt);
                                } else {
                                    violations.push(VvViolation::new(
                                        VvRule::QoiDependencyIsolated,
                                        Some(prediction.id().as_str().to_owned()),
                                        Some(qoi_id.as_str().to_owned()),
                                        "prediction.solution_verification",
                                        "solution receipt must match this exact QoI and unit",
                                    ));
                                }
                            }
                        }
                        _ => violations.push(VvViolation::new(
                            VvRule::SolutionVerificationComplete,
                            Some(prediction.id().as_str().to_owned()),
                            Some(qoi_id.as_str().to_owned()),
                            "prediction.solution_verification",
                            "solution-verification dependency must target a V&V solution receipt",
                        )),
                    }
                    if dependency.observations.is_some() {
                        violations.push(VvViolation::new(
                            VvRule::QoiDependencyIsolated,
                            Some(prediction.id().as_str().to_owned()),
                            Some(qoi_id.as_str().to_owned()),
                            "prediction.solution_verification",
                            "solution-verification edges cannot carry observation selections",
                        ));
                    }
                }
                DependencyRole::ProcessConformance => {
                    process_targets.insert(target_hash);
                    if !matches!(dependency.target(), EvidenceTarget::External { .. })
                        || dependency.observations.is_some()
                    {
                        violations.push(VvViolation::new(
                            VvRule::ProcessConformanceSeparate,
                            Some(prediction.id().as_str().to_owned()),
                            Some(qoi_id.as_str().to_owned()),
                            "prediction.process_conformance",
                            "process conformance must be a separate external receipt without observations",
                        ));
                    }
                }
                _ if dependency.observations.is_some() => violations.push(VvViolation::new(
                    VvRule::QoiDependencyIsolated,
                    Some(prediction.id().as_str().to_owned()),
                    Some(qoi_id.as_str().to_owned()),
                    "prediction.dependencies.observations",
                    "only physical-validation edges may carry held-out observations",
                )),
                _ => {}
            }
        }
        if physical_targets.is_empty() || solution_targets.is_empty() {
            violations.push(VvViolation::new(
                VvRule::QoiDependencyClosed,
                Some(prediction.id().as_str().to_owned()),
                Some(qoi_id.as_str().to_owned()),
                "prediction.dependencies",
                "each prediction needs exact physical-validation and solution-verification edges",
            ));
        }
        self.validate_prediction_metrics(
            prediction,
            plan_row,
            &referenced_solutions,
            hashes,
            violations,
        );
        self.validate_waterfall(
            prediction,
            &referenced_solutions,
            &process_targets,
            violations,
        );
        self.validate_prediction_applicability(prediction, violations);
        self.validate_evidence_axes(
            prediction,
            hashes,
            &role_targets,
            &process_targets,
            violations,
        );
    }

    fn validate_physical_dependency(
        &self,
        prediction: &PredictionAssessment,
        dependency: &EvidenceDependency,
        plan_row: &QoiValidationPlan,
        hashes: &ArtifactHashMap,
        violations: &mut Vec<VvViolation>,
    ) {
        let qoi = prediction.qoi();
        let (reference, experiment) = match dependency.target() {
            EvidenceTarget::VvArtifact(reference)
                if reference.kind == ArtifactKind::ExperimentArtifact =>
            {
                let Some(experiment) = self.experiments.get(reference.id()) else {
                    return;
                };
                (reference, experiment)
            }
            _ => {
                violations.push(VvViolation::new(
                    VvRule::ValidationRequiresPhysicalReferent,
                    Some(prediction.id().as_str().to_owned()),
                    Some(qoi.as_str().to_owned()),
                    "prediction.physical_validation",
                    "physical validation must target an admitted ExperimentArtifact",
                ));
                return;
            }
        };
        if !experiment.origin.is_physical()
            || !experiment.qois.contains(qoi)
            || !plan_row.experiments.contains(reference)
        {
            violations.push(VvViolation::new(
                VvRule::ValidationRequiresPhysicalReferent,
                Some(prediction.id().as_str().to_owned()),
                Some(qoi.as_str().to_owned()),
                "prediction.physical_validation",
                "only a plan-declared physical experiment for this QoI can support validation",
            ));
        }
        let Some(selection) = dependency.observations() else {
            violations.push(VvViolation::new(
                VvRule::ValidationRequiresPhysicalReferent,
                Some(prediction.id().as_str().to_owned()),
                Some(qoi.as_str().to_owned()),
                "prediction.physical_validation",
                "physical validation requires a held-out or released-blind observation selection",
            ));
            return;
        };
        violations.extend(selection_violations(
            selection,
            &self.splits,
            hashes,
            prediction.id(),
            qoi,
            "prediction.physical_validation.observations",
        ));
        if selection.split != plan_row.split {
            violations.push(VvViolation::new(
                VvRule::SplitPartitionsDisjoint,
                Some(prediction.id().as_str().to_owned()),
                Some(qoi.as_str().to_owned()),
                "prediction.physical_validation.observations",
                "physical-validation observations must come from the QoI's preregistered split",
            ));
        }
        if let Some(split) = self.splits.get(selection.split.id()) {
            if split.experiment != *reference {
                violations.push(VvViolation::new(
                    VvRule::ValidationRequiresPhysicalReferent,
                    Some(prediction.id().as_str().to_owned()),
                    Some(qoi.as_str().to_owned()),
                    "prediction.physical_validation.observations",
                    "observation split and physical experiment dependency do not match",
                ));
            }
        }
    }

    fn validate_prediction_metrics(
        &self,
        prediction: &PredictionAssessment,
        plan_row: &QoiValidationPlan,
        referenced_solutions: &[&SolutionVerificationReceipt],
        hashes: &ArtifactHashMap,
        violations: &mut Vec<VvViolation>,
    ) {
        let qoi = prediction.qoi();
        if prediction.validation_metrics.is_empty() {
            violations.push(VvViolation::new(
                VvRule::ValidationMetricUncertainty,
                Some(prediction.id().as_str().to_owned()),
                Some(qoi.as_str().to_owned()),
                "prediction.validation_metrics",
                "physical prediction assessment requires actual validation metrics",
            ));
        }
        let required_numerical = referenced_solutions
            .iter()
            .map(|receipt| receipt.combined_half_width())
            .fold(0.0_f64, f64::max);
        for metric in &prediction.validation_metrics {
            violations.extend(selection_violations(
                metric.observations(),
                &self.splits,
                hashes,
                prediction.id(),
                qoi,
                "prediction.validation_metrics.observations",
            ));
            if metric.observations.split != plan_row.split {
                violations.push(VvViolation::new(
                    VvRule::SplitPartitionsDisjoint,
                    Some(prediction.id().as_str().to_owned()),
                    Some(qoi.as_str().to_owned()),
                    "prediction.validation_metrics.observations",
                    "validation metric must use the QoI's preregistered held-out split",
                ));
            }
            if metric.numerical_uncertainty < required_numerical {
                violations.push(VvViolation::new(
                    VvRule::ValidationMetricUncertainty,
                    Some(prediction.id().as_str().to_owned()),
                    Some(qoi.as_str().to_owned()),
                    "prediction.validation_metrics.numerical_uncertainty",
                    "validation metric understates the admitted solution-verification envelope",
                ));
            }
        }
        for check in &prediction.posterior_checks {
            violations.extend(selection_violations(
                check.observations(),
                &self.splits,
                hashes,
                prediction.id(),
                qoi,
                "prediction.posterior_checks.observations",
            ));
            if check.observations.split != plan_row.split {
                violations.push(VvViolation::new(
                    VvRule::SplitPartitionsDisjoint,
                    Some(prediction.id().as_str().to_owned()),
                    Some(qoi.as_str().to_owned()),
                    "prediction.posterior_checks.observations",
                    "posterior check must use the QoI's preregistered held-out split",
                ));
            }
        }
        let posterior_required = plan_row
            .metrics
            .iter()
            .any(|metric| matches!(metric, ValidationMetricSpec::PosteriorPredictive { .. }));
        if posterior_required && prediction.posterior_checks.is_empty() {
            violations.push(VvViolation::new(
                VvRule::ValidationMetricUncertainty,
                Some(prediction.id().as_str().to_owned()),
                Some(qoi.as_str().to_owned()),
                "prediction.posterior_checks",
                "the validation plan requires a posterior-predictive check",
            ));
        }
    }

    fn validate_waterfall(
        &self,
        prediction: &PredictionAssessment,
        referenced_solutions: &[&SolutionVerificationReceipt],
        process_targets: &BTreeSet<ContentHash>,
        violations: &mut Vec<VvViolation>,
    ) {
        let dependency_targets = prediction
            .dependencies
            .iter()
            .map(|dependency| dependency.target.clone())
            .collect::<BTreeSet<_>>();
        let solution_targets = prediction
            .dependencies
            .iter()
            .filter(|dependency| dependency.role == DependencyRole::SolutionVerification)
            .map(|dependency| dependency.target.clone())
            .collect::<BTreeSet<_>>();
        let numerical_floor = referenced_solutions
            .iter()
            .map(|receipt| receipt.combined_half_width())
            .fold(0.0_f64, f64::max);
        for term in &prediction.waterfall.terms {
            let hash = term.source.hash();
            if !dependency_targets.contains(&term.source) {
                violations.push(VvViolation::new(
                    VvRule::QoiDependencyClosed,
                    Some(prediction.id().as_str().to_owned()),
                    Some(prediction.qoi().as_str().to_owned()),
                    "prediction.waterfall.source",
                    "every uncertainty source must be an exact QoI dependency",
                ));
            }
            if process_targets.contains(&hash) {
                violations.push(VvViolation::new(
                    VvRule::ProcessConformanceSeparate,
                    Some(prediction.id().as_str().to_owned()),
                    Some(prediction.qoi().as_str().to_owned()),
                    "prediction.waterfall.source",
                    "process-conformance evidence cannot substitute for uncertainty evidence",
                ));
            }
            if term.kind == PredictionUncertaintyKind::Numerical
                && (!solution_targets.contains(&term.source) || term.magnitude() < numerical_floor)
            {
                violations.push(VvViolation::new(
                    VvRule::SolutionVerificationComplete,
                    Some(prediction.id().as_str().to_owned()),
                    Some(prediction.qoi().as_str().to_owned()),
                    "prediction.waterfall.numerical",
                    "numerical waterfall term must link a solution receipt without understating it",
                ));
            }
        }
    }

    fn validate_prediction_applicability(
        &self,
        prediction: &PredictionAssessment,
        violations: &mut Vec<VvViolation>,
    ) {
        let expected_checks = self.assumptions.rows.keys().collect::<BTreeSet<_>>();
        let actual_checks = prediction.assumption_checks.keys().collect::<BTreeSet<_>>();
        if expected_checks != actual_checks {
            violations.push(VvViolation::new(
                VvRule::AssumptionRowComplete,
                Some(prediction.id().as_str().to_owned()),
                Some(prediction.qoi().as_str().to_owned()),
                "prediction.assumption_checks",
                "prediction must report every ledger assumption exactly once",
            ));
        }
        let mut domain_violations = self
            .context
            .applicability
            .violations(&prediction.applicability_point);
        let has_domain_violation = !domain_violations.is_empty();
        let mut assumption_forces_refusal = false;
        for (id, row) in &self.assumptions.rows {
            let declared_pass = prediction
                .assumption_checks
                .get(id)
                .copied()
                .unwrap_or(false);
            // `RuntimeMonitorSpec` and `AssumptionEvidence` are mandatory,
            // non-blank schema fields. Their optional hashes refine retained
            // evidence; the per-assessment Boolean is the runtime predicate
            // result and cannot be inferred from hash presence alone.
            if !declared_pass {
                domain_violations.push(DomainViolation::Assumption { id: id.clone() });
                assumption_forces_refusal |= matches!(
                    &row.violation_effect,
                    ViolationEffect::EscalateOrRefuse { .. } | ViolationEffect::Refuse { .. }
                );
            }
        }
        let expected = if domain_violations.is_empty() {
            ApplicabilityDecision::InDomain
        } else if assumption_forces_refusal
            || (has_domain_violation
                && self.context.applicability_policy == ApplicabilityPolicy::Refuse)
        {
            ApplicabilityDecision::Refused {
                violations: domain_violations,
            }
        } else {
            ApplicabilityDecision::Demoted {
                violations: domain_violations,
            }
        };
        if !prediction.applicability.has_same_canonical_bits(&expected) {
            violations.push(VvViolation::new(
                VvRule::ApplicabilityDecision,
                Some(prediction.id().as_str().to_owned()),
                Some(prediction.qoi().as_str().to_owned()),
                "prediction.applicability",
                "applicability decision does not equal the domain, policy, monitor, and assumption result",
            ));
        }
    }

    fn validate_evidence_axes(
        &self,
        prediction: &PredictionAssessment,
        hashes: &ArtifactHashMap,
        role_targets: &BTreeMap<DependencyRole, BTreeSet<ContentHash>>,
        process_targets: &BTreeSet<ContentHash>,
        violations: &mut Vec<VvViolation>,
    ) {
        let targets = |role| role_targets.get(&role).cloned().unwrap_or_default();
        let mut allowed = BTreeMap::<EvidenceAxis, BTreeSet<ContentHash>>::new();
        allowed.insert(
            EvidenceAxis::CodeVerification,
            targets(DependencyRole::CodeVerification),
        );
        allowed.insert(
            EvidenceAxis::SolutionVerification,
            targets(DependencyRole::SolutionVerification),
        );
        let mut numerical = targets(DependencyRole::SolutionVerification);
        for receipt in self
            .solution_verification
            .values()
            .filter(|receipt| receipt.qoi() == prediction.qoi())
        {
            numerical.extend([
                receipt.mesh.evidence_hash,
                receipt.time.evidence_hash,
                receipt.nonlinear.evidence_hash,
                receipt.iterative.evidence_hash,
            ]);
        }
        allowed.insert(EvidenceAxis::NumericalUncertainty, numerical);
        let mut parameter_data = targets(DependencyRole::ParameterData);
        for term in &prediction.waterfall.terms {
            if matches!(
                term.kind,
                PredictionUncertaintyKind::Parameter | PredictionUncertaintyKind::Data
            ) {
                parameter_data.insert(term.source.hash());
            }
        }
        allowed.insert(EvidenceAxis::ParameterDataUncertainty, parameter_data);
        let mut model_form = targets(DependencyRole::ModelDiscrepancy);
        model_form.extend(targets(DependencyRole::PhysicalValidation));
        model_form.extend(targets(DependencyRole::PosteriorPredictive));
        model_form.extend(
            prediction
                .posterior_checks
                .iter()
                .map(|check| check.artifact_hash),
        );
        allowed.insert(EvidenceAxis::ModelFormValidation, model_form);
        let mut domain = BTreeSet::new();
        if let Some(hash) = hashes.get(&(ArtifactKind::ContextOfUse, self.context.id().clone())) {
            domain.insert(*hash);
        }
        for row in self.assumptions.rows.values() {
            if let Some(target) = row.evidence.artifact() {
                domain.insert(target.hash());
            }
            if let Some(hash) = row.monitor.evidence_hash {
                domain.insert(hash);
            }
        }
        allowed.insert(EvidenceAxis::PredictionDomainRelevance, domain);
        allowed.insert(
            EvidenceAxis::ComparisonToExperiment,
            targets(DependencyRole::PhysicalValidation),
        );

        for (axis, status) in &prediction.evidence_axes.axes {
            let EvidenceAxisStatus::Present { artifacts } = status else {
                continue;
            };
            let permitted = allowed.get(axis).cloned().unwrap_or_default();
            if artifacts
                .iter()
                .any(|hash| !permitted.contains(hash) || process_targets.contains(hash))
            {
                violations.push(VvViolation::new(
                    if artifacts.iter().any(|hash| process_targets.contains(hash)) {
                        VvRule::ProcessConformanceSeparate
                    } else {
                        VvRule::QoiDependencyClosed
                    },
                    Some(prediction.id().as_str().to_owned()),
                    Some(prediction.qoi().as_str().to_owned()),
                    "prediction.evidence_axes",
                    "present evidence-axis hashes must be exact, role-compatible dependency closure",
                ));
            }
        }
    }
}
