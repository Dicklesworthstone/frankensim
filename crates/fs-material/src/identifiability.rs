//! Authority-separated law/experiment identifiability schemas (I10.1).
//!
//! An identifiability statement is meaningful only for one exact constitutive
//! law, material card, state schema, specimen/process, protocol/refinement,
//! observation model, covariance, nuisance/discrepancy policy, and data split.
//! The current public API closes those inputs in four deliberately distinct
//! stages:
//!
//! - [`IdentifiabilityProblemDocument`] is an unresolved, coordinate-free
//!   physical and statistical question;
//! - [`AdmittedIdentifiabilityProblem`] binds exact source content and retains
//!   separate [`ProblemId`] and [`SourceAdmissionId`] identities;
//! - [`IdentifiabilityExecutionPlan`] binds coordinates, algorithms, numerical
//!   policy, budgets, seeds, builds, and replay authority in an [`ExecutionId`];
//! - [`IdentifiabilityAssessment`] binds product-typed claims and evidence in an
//!   [`AssessmentId`] without silently promoting a content receipt to a theorem.
//!
//! The older single-case schema lives only in repository history and the
//! explicitly compile-disabled archaeological fixture. Its former identities
//! cannot mint authority in the current multi-case chain.
//! No identity in this module authenticates a laboratory or proves scientific
//! correctness by itself; those remain explicit evidence and trust-policy
//! obligations.

mod authoritative;

pub use authoritative::*;

use core::fmt;
use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::vv::{
    ArtifactHeader, ArtifactId, ArtifactKind, ArtifactRef, BlindReleaseReceipt, CalibrationSplit,
    ContextOfUse, CovarianceMatrix, DeclaredBudget, ExperimentArtifact, ObservationId,
    ObservationManifest, ObservationManifestRow, QoiId, SeedDeclaration, UnitId, VV_SCHEMA_VERSION,
};
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, MATDB_SCHEMA_VERSION, MaterialCard,
};
use fs_qty::{Dims, QUANTITY_SPEC_ENCODED_LEN, QuantitySpec};

/// Maximum identifier or short-reason byte length.
pub const MAX_IDENTIFIABILITY_ID_BYTES: usize = 256;
/// Maximum long diagnostic/reason byte length.
pub const MAX_IDENTIFIABILITY_TEXT_BYTES: usize = 16 * 1024;
/// Maximum rows in any parameter/observation/path/gauge collection.
pub const MAX_IDENTIFIABILITY_ITEMS: usize = 4096;
/// Maximum canonical study bytes accepted or emitted.
pub const MAX_IDENTIFIABILITY_CANONICAL_BYTES: usize = 4 * 1024 * 1024;

fn hash_is_nonzero(hash: ContentHash) -> bool {
    hash.as_bytes().iter().any(|byte| *byte != 0)
}

fn canonical_f64(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}

fn same_f64(left: f64, right: f64) -> bool {
    canonical_f64(left).to_bits() == canonical_f64(right).to_bits()
}

fn validate_token(value: &str, field: &'static str) -> Result<(), IdentifiabilityError> {
    if value.is_empty()
        || value.len() > MAX_IDENTIFIABILITY_ID_BYTES
        || !value.bytes().all(|byte| byte.is_ascii_graphic())
    {
        return Err(IdentifiabilityError::InvalidText {
            field,
            detail: format!(
                "expected a nonempty ASCII machine token without whitespace of at most {MAX_IDENTIFIABILITY_ID_BYTES} bytes"
            ),
        });
    }
    Ok(())
}

fn validate_reason(value: &str, field: &'static str) -> Result<(), IdentifiabilityError> {
    if value.is_empty()
        || value.len() > MAX_IDENTIFIABILITY_TEXT_BYTES
        || value.trim() != value
        || value.chars().any(char::is_control)
    {
        return Err(IdentifiabilityError::InvalidText {
            field,
            detail: format!(
                "expected nonempty trimmed text of at most {MAX_IDENTIFIABILITY_TEXT_BYTES} bytes"
            ),
        });
    }
    Ok(())
}

macro_rules! typed_id {
    ($name:ident, $field:literal) => {
        #[doc = concat!("Typed ", $field, " identifier.")]
        #[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(String);

        impl $name {
            #[doc = concat!("Construct a bounded ", $field, " identifier.")]
            pub fn try_new(value: impl Into<String>) -> Result<Self, IdentifiabilityError> {
                let value = value.into();
                validate_token(&value, $field)?;
                Ok(Self(value))
            }

            /// Inspect the canonical identifier text.
            #[must_use]
            pub fn as_str(&self) -> &str {
                &self.0
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str(&self.0)
            }
        }
    };
}

typed_id!(ParameterRoleId, "parameter role");
typed_id!(CoordinateId, "estimation coordinate");
typed_id!(ObservationChannelId, "observation channel");
typed_id!(GaugeClassId, "gauge class");

/// Deterministic refusal at the identifiability-study boundary.
#[derive(Debug, Clone, PartialEq)]
pub enum IdentifiabilityError {
    /// Retained schema is stale or from the future.
    UnsupportedSchemaVersion {
        /// Version found in the supplied artifact.
        declared: u32,
        /// Sole version accepted by this implementation.
        supported: u32,
    },
    /// A bounded identifier/reason was malformed.
    InvalidText {
        /// Semantic field that failed validation.
        field: &'static str,
        /// Bounded refusal detail.
        detail: String,
    },
    /// A required content identity was the all-zero sentinel.
    ZeroIdentity {
        /// Semantic field carrying the zero identity.
        field: &'static str,
    },
    /// A count was empty, oversized, or inconsistent.
    Cardinality {
        /// Collection or structure whose cardinality was invalid.
        field: &'static str,
        /// Bounded refusal detail.
        detail: String,
    },
    /// Two rows claimed the same identity.
    Duplicate {
        /// Identity class in which a duplicate occurred.
        field: &'static str,
        /// Duplicated identity text.
        id: String,
    },
    /// A numeric interval/prior/transform was malformed.
    InvalidNumeric {
        /// Numeric field or relation that was invalid.
        field: &'static str,
        /// Bounded refusal detail.
        detail: String,
    },
    /// Exact version pins disagree.
    VersionMismatch {
        /// Versioned component that disagreed.
        field: &'static str,
        /// Required version.
        expected: u32,
        /// Supplied version.
        actual: u32,
    },
    /// One reference points outside the admitted closed graph.
    UnknownReference {
        /// Reference class that failed resolution.
        field: &'static str,
        /// Unresolved identity text.
        id: String,
    },
    /// A model card is not an exact member of the material card.
    ModelNotInMaterialCard,
    /// A model/state initialization policy was violated.
    InitialStatePolicy {
        /// Bounded policy mismatch detail.
        detail: String,
    },
    /// An estimated parameter has no observation path and no honest refusal.
    DisconnectedEstimatedParameter {
        /// Estimated parameter lacking an observation path.
        parameter: ParameterRoleId,
    },
    /// A nuisance parameter is not calibrated by this exact split.
    NuisanceCalibration {
        /// Nuisance parameter lacking calibration authorization.
        parameter: ParameterRoleId,
    },
    /// A declared gauge quotient is structurally invalid.
    InvalidGauge {
        /// Gauge class that failed validation.
        gauge: GaugeClassId,
        /// Bounded refusal detail.
        detail: String,
    },
    /// Covariance order/dimension disagrees with the observation schema.
    Covariance {
        /// Bounded covariance refusal detail.
        detail: String,
    },
    /// V&V artifact construction or canonicalization failed.
    Vv {
        /// Bounded V&V refusal detail.
        detail: String,
    },
    /// Material-card validation or identity failed.
    Material {
        /// Bounded material-card refusal detail.
        detail: String,
    },
    /// Canonical claims disagree with caller-resolved source artifacts.
    SourceMismatch {
        /// Canonical claim that disagreed with its source.
        field: &'static str,
    },
    /// Canonical transport is malformed or exceeds its public cap.
    Canonical {
        /// Byte offset at which transport validation refused.
        at: usize,
        /// Bounded transport refusal detail.
        detail: String,
    },
}

impl fmt::Display for IdentifiabilityError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedSchemaVersion {
                declared,
                supported,
            } => write!(
                f,
                "identifiability schema v{declared} is unsupported; expected exactly v{supported}"
            ),
            Self::InvalidText { field, detail } | Self::InvalidNumeric { field, detail } => {
                write!(f, "invalid {field}: {detail}")
            }
            Self::ZeroIdentity { field } => write!(f, "{field} uses the all-zero identity"),
            Self::Cardinality { field, detail } => {
                write!(f, "invalid {field} cardinality: {detail}")
            }
            Self::Duplicate { field, id } => write!(f, "duplicate {field} identity {id:?}"),
            Self::VersionMismatch {
                field,
                expected,
                actual,
            } => write!(
                f,
                "{field} version mismatch: expected {expected}, found {actual}"
            ),
            Self::UnknownReference { field, id } => {
                write!(f, "{field} references unknown identity {id:?}")
            }
            Self::ModelNotInMaterialCard => {
                f.write_str("constitutive model card is not an exact member of the material card")
            }
            Self::InitialStatePolicy { detail } => {
                write!(f, "initial-state policy mismatch: {detail}")
            }
            Self::DisconnectedEstimatedParameter { parameter } => write!(
                f,
                "estimated parameter {parameter} has no declared observation path and is not explicitly unidentifiable"
            ),
            Self::NuisanceCalibration { parameter } => write!(
                f,
                "nuisance parameter {parameter} is not calibrated by the study's exact split"
            ),
            Self::InvalidGauge { gauge, detail } => {
                write!(f, "invalid gauge class {gauge}: {detail}")
            }
            Self::Covariance { detail } => write!(f, "invalid observation covariance: {detail}"),
            Self::Vv { detail } => write!(f, "V&V artifact refusal: {detail}"),
            Self::Material { detail } => write!(f, "material-card refusal: {detail}"),
            Self::SourceMismatch { field } => {
                write!(
                    f,
                    "canonical {field} claim disagrees with the resolved source artifact"
                )
            }
            Self::Canonical { at, detail } => {
                write!(
                    f,
                    "canonical study transport refused at byte {at}: {detail}"
                )
            }
        }
    }
}

impl std::error::Error for IdentifiabilityError {}

/// Closed finite interval in coherent SI coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ParameterDomain {
    lo: f64,
    hi: f64,
}

impl ParameterDomain {
    /// Construct a finite ordered domain. Degenerate domains are permitted for
    /// fixed parameters but not for estimated/nuisance parameters.
    pub fn try_new(lo: f64, hi: f64) -> Result<Self, IdentifiabilityError> {
        if !lo.is_finite() || !hi.is_finite() || lo > hi {
            return Err(IdentifiabilityError::InvalidNumeric {
                field: "parameter domain",
                detail: format!("expected finite lo <= hi, got [{lo:?}, {hi:?}]"),
            });
        }
        Ok(Self {
            lo: canonical_f64(lo),
            hi: canonical_f64(hi),
        })
    }

    /// Inclusive bounds.
    #[must_use]
    pub const fn bounds(self) -> (f64, f64) {
        (self.lo, self.hi)
    }

    fn is_degenerate(self) -> bool {
        same_f64(self.lo, self.hi)
    }
}

/// Prior on the canonical physical parameter, never on a transient optimizer
/// coordinate. This placement is what lets reparameterized studies share a
/// sound quotient identity.
#[derive(Debug, Clone, PartialEq)]
pub enum ParameterPrior {
    /// Prior deliberately absent; the reason and policy version stay bound.
    None {
        /// Bounded reason no prior is declared.
        reason: String,
        /// Prior-policy semantics version.
        version: u32,
    },
    /// Uniform prior over a canonical physical interval.
    Uniform {
        /// Closed support interval in coherent SI coordinates.
        domain: ParameterDomain,
        /// Prior-family semantics version.
        version: u32,
    },
    /// Gaussian density in coherent SI units, normalized after conditioning
    /// on the enclosing [`StudyParameter`] physical domain.
    Gaussian {
        /// Mean in coherent SI coordinates.
        mean: f64,
        /// Positive standard deviation in coherent SI coordinates.
        standard_deviation: f64,
        /// Prior-family semantics version.
        version: u32,
    },
    /// Log-normal density over a positive dimensional parameter, normalized
    /// after conditioning on the enclosing physical domain. `reference`
    /// supplies the coherent-SI scale used by the logarithm.
    LogNormal {
        /// Mean of the dimensionless logarithmic coordinate.
        log_mean: f64,
        /// Positive standard deviation of the logarithmic coordinate.
        log_standard_deviation: f64,
        /// Positive coherent-SI reference scale used by the logarithm.
        reference: f64,
        /// Prior-family semantics version.
        version: u32,
    },
}

impl ParameterPrior {
    fn validate_against(
        &mut self,
        parameter_domain: ParameterDomain,
    ) -> Result<(), IdentifiabilityError> {
        let version = match self {
            Self::None { reason, version } => {
                validate_reason(reason, "prior absence reason")?;
                *version
            }
            Self::Uniform { domain, version } => {
                if domain.is_degenerate()
                    || domain.lo < parameter_domain.lo
                    || domain.hi > parameter_domain.hi
                {
                    return Err(IdentifiabilityError::InvalidNumeric {
                        field: "uniform prior support",
                        detail: "uniform-prior support must have positive width and lie inside the physical parameter domain; atomic mass requires an explicit discrete/Dirac prior semantics"
                            .to_string(),
                    });
                }
                *version
            }
            Self::Gaussian {
                mean,
                standard_deviation,
                version,
            } => {
                if !mean.is_finite()
                    || !standard_deviation.is_finite()
                    || *standard_deviation <= 0.0
                {
                    return Err(IdentifiabilityError::InvalidNumeric {
                        field: "Gaussian prior",
                        detail: "mean must be finite and standard deviation positive".to_string(),
                    });
                }
                *mean = canonical_f64(*mean);
                *version
            }
            Self::LogNormal {
                log_mean,
                log_standard_deviation,
                reference,
                version,
            } => {
                if parameter_domain.lo <= 0.0 {
                    return Err(IdentifiabilityError::InvalidNumeric {
                        field: "log-normal prior support",
                        detail: "log-normal prior requires an entirely positive physical domain"
                            .to_string(),
                    });
                }
                if !log_mean.is_finite()
                    || !log_standard_deviation.is_finite()
                    || *log_standard_deviation <= 0.0
                    || !reference.is_finite()
                    || *reference <= 0.0
                {
                    return Err(IdentifiabilityError::InvalidNumeric {
                        field: "log-normal prior",
                        detail:
                            "log moments must be finite, spreads positive, and reference positive"
                                .to_string(),
                    });
                }
                *log_mean = canonical_f64(*log_mean);
                *version
            }
        };
        if version == 0 {
            return Err(IdentifiabilityError::InvalidNumeric {
                field: "prior version",
                detail: "version zero is not a published prior semantics".to_string(),
            });
        }
        Ok(())
    }

    /// Declared prior-family semantics version.
    #[must_use]
    pub const fn version(&self) -> u32 {
        match self {
            Self::None { version, .. }
            | Self::Uniform { version, .. }
            | Self::Gaussian { version, .. }
            | Self::LogNormal { version, .. } => *version,
        }
    }
}

/// Bijective coordinate chart mapping an optimizer coordinate to the
/// canonical physical parameter.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CoordinateTransform {
    /// Canonical value equals coordinate value.
    Identity,
    /// `physical = scale * coordinate + offset`, with a typed nonzero scale
    /// and a coherent-SI physical offset.
    Affine {
        /// Nonzero scale multiplying the optimizer coordinate.
        scale: f64,
        /// Quantity carried by the affine scale.
        scale_quantity: QuantitySpec,
        /// Coherent-SI physical offset.
        offset: f64,
    },
    /// `physical = reference * exp(coordinate)`, using deterministic fs-math.
    LogPositive {
        /// Positive coherent-SI reference scale.
        reference: f64,
    },
}

impl CoordinateTransform {
    fn validate(self) -> Result<Self, IdentifiabilityError> {
        match self {
            Self::Identity => Ok(self),
            Self::Affine {
                scale,
                scale_quantity,
                offset,
            } if scale.is_finite() && scale != 0.0 && offset.is_finite() => Ok(Self::Affine {
                scale: canonical_f64(scale),
                scale_quantity,
                offset: canonical_f64(offset),
            }),
            Self::LogPositive { reference } if reference.is_finite() && reference > 0.0 => Ok(self),
            Self::Affine { .. } => Err(IdentifiabilityError::InvalidNumeric {
                field: "affine coordinate transform",
                detail: "scale must be finite/nonzero and offset finite".to_string(),
            }),
            Self::LogPositive { .. } => Err(IdentifiabilityError::InvalidNumeric {
                field: "log coordinate transform",
                detail: "reference must be finite and positive".to_string(),
            }),
        }
    }

    fn map(self, value: f64) -> f64 {
        match self {
            Self::Identity => value,
            Self::Affine {
                scale,
                scale_quantity: _,
                offset,
            } => scale.mul_add(value, offset),
            Self::LogPositive { reference } => reference * fs_math::det::exp(value),
        }
    }

    fn mapped_domain(
        self,
        domain: ParameterDomain,
    ) -> Result<ParameterDomain, IdentifiabilityError> {
        let (lo, hi) = domain.bounds();
        let left = self.map(lo);
        let right = self.map(hi);
        ParameterDomain::try_new(left.min(right), left.max(right))
    }
}

/// Exact optimizer coordinate used to represent one canonical parameter.
#[derive(Debug, Clone, PartialEq)]
pub struct ParameterCoordinate {
    id: CoordinateId,
    quantity: QuantitySpec,
    domain: ParameterDomain,
    transform: CoordinateTransform,
}

impl ParameterCoordinate {
    /// Construct one coordinate chart and validate the transform's local
    /// numeric invariants. The owning schema checks quantity compatibility,
    /// full-domain bijectivity, and physical-domain coverage when the current
    /// authority-separated execution stage admits it through
    /// [`IdentifiabilityExecutionPlan::try_new`].
    pub fn try_new(
        id: CoordinateId,
        quantity: QuantitySpec,
        domain: ParameterDomain,
        transform: CoordinateTransform,
    ) -> Result<Self, IdentifiabilityError> {
        Ok(Self {
            id,
            quantity,
            domain,
            transform: transform.validate()?,
        })
    }

    /// Coordinate identity.
    #[must_use]
    pub const fn id(&self) -> &CoordinateId {
        &self.id
    }

    /// Quantity descriptor of the optimizer coordinate.
    #[must_use]
    pub const fn quantity(&self) -> QuantitySpec {
        self.quantity
    }

    /// Coordinate-domain bounds.
    #[must_use]
    pub const fn domain(&self) -> ParameterDomain {
        self.domain
    }

    /// Coordinate-to-physical map.
    #[must_use]
    pub const fn transform(&self) -> CoordinateTransform {
        self.transform
    }
}

/// Quantity and exact nominal coherent-SI value of one model-card parameter.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ModelParameterBinding {
    quantity: QuantitySpec,
    nominal_bits: u64,
}

impl ModelParameterBinding {
    /// Full semantic quantity descriptor inherited from the model card.
    #[must_use]
    pub const fn quantity(&self) -> QuantitySpec {
        self.quantity
    }

    /// Exact coherent-SI nominal value inherited from the model card.
    #[must_use]
    pub fn nominal(&self) -> f64 {
        f64::from_bits(self.nominal_bits)
    }
}

/// Exact immutable material/law/parameter/state binding consumed by a study.
///
/// The caller-supplied graph digest is content binding only: `ConstitutiveGraph`
/// does not yet own a semantic identity, so this type does not authenticate the
/// digest or pretend its current FNV state-layout version is a graph identity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MaterialModelBinding {
    material_card: ContentHash,
    model_card: ContentHash,
    parameter_block: ContentHash,
    graph: ContentHash,
    law: LawId,
    law_version: u32,
    state_schema_version: u32,
    initial_state_policy: InitialStatePolicy,
    matdb_schema_version: u32,
    parameter_roster: BTreeMap<ParameterRoleId, ModelParameterBinding>,
}

impl MaterialModelBinding {
    /// Bind an exact model-card member of an immutable material card.
    ///
    /// # Errors
    /// Refuses invalid cards, a model not present byte-for-byte in `material`,
    /// malformed parameter role names, duplicate roles, or an all-zero graph
    /// content binding.
    pub fn from_cards(
        material: &MaterialCard,
        model: &ConstitutiveModelCard,
        graph: ContentHash,
    ) -> Result<Self, IdentifiabilityError> {
        model
            .validate()
            .map_err(|error| IdentifiabilityError::Material {
                detail: error.to_string(),
            })?;
        if !hash_is_nonzero(graph) {
            return Err(IdentifiabilityError::ZeroIdentity {
                field: "constitutive graph binding",
            });
        }
        let model_card = model.content_hash();
        if !material
            .models()
            .iter()
            .any(|candidate| candidate.content_hash() == model_card)
        {
            return Err(IdentifiabilityError::ModelNotInMaterialCard);
        }
        validate_token(&model.law.0, "constitutive law id")?;
        let mut parameter_roster = BTreeMap::new();
        for (name, parameter) in &model.parameters {
            let role = ParameterRoleId::try_new(name.clone())?;
            if parameter_roster
                .insert(
                    role.clone(),
                    ModelParameterBinding {
                        quantity: QuantitySpec::dimensional(parameter.dims),
                        nominal_bits: canonical_f64(parameter.value).to_bits(),
                    },
                )
                .is_some()
            {
                return Err(IdentifiabilityError::Duplicate {
                    field: "model parameter",
                    id: role.to_string(),
                });
            }
        }
        let parameter_block =
            model
                .canonical_parameters_hash()
                .map_err(|error| IdentifiabilityError::Material {
                    detail: error.to_string(),
                })?;
        Ok(Self {
            material_card: material.content_hash(),
            model_card,
            parameter_block,
            graph,
            law: model.law.clone(),
            law_version: model.law_version,
            state_schema_version: model.state_schema_version,
            initial_state_policy: model.initial_state,
            matdb_schema_version: material.schema_version(),
            parameter_roster,
        })
    }

    /// Exact material-card identity.
    #[must_use]
    pub const fn material_card(&self) -> ContentHash {
        self.material_card
    }

    /// Exact constitutive-model-card identity.
    #[must_use]
    pub const fn model_card(&self) -> ContentHash {
        self.model_card
    }

    /// Narrow canonical parameter-block identity.
    #[must_use]
    pub const fn parameter_block(&self) -> ContentHash {
        self.parameter_block
    }

    /// Caller-supplied constitutive-graph content binding.
    #[must_use]
    pub const fn graph(&self) -> ContentHash {
        self.graph
    }

    /// Constitutive law identity.
    #[must_use]
    pub const fn law(&self) -> &LawId {
        &self.law
    }

    /// Exact law semantics version.
    #[must_use]
    pub const fn law_version(&self) -> u32 {
        self.law_version
    }

    /// Exact internal-state schema version.
    #[must_use]
    pub const fn state_schema_version(&self) -> u32 {
        self.state_schema_version
    }

    /// Initial-state policy declared by the exact model card.
    #[must_use]
    pub const fn initial_state_policy(&self) -> InitialStatePolicy {
        self.initial_state_policy
    }

    /// Material-database schema version under which both cards were admitted.
    #[must_use]
    pub const fn matdb_schema_version(&self) -> u32 {
        self.matdb_schema_version
    }

    /// Parameter roles/dimensions from the exact model card.
    #[must_use]
    pub const fn parameter_roster(&self) -> &BTreeMap<ParameterRoleId, ModelParameterBinding> {
        &self.parameter_roster
    }
}

/// Exact initial-state binding for the constitutive state schema.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InitialStateBinding {
    /// The model card authorizes the canonical all-zero internal state.
    Zero {
        /// Exact internal-state schema version.
        schema_version: u32,
    },
    /// An explicit state artifact is required and retained.
    Explicit {
        /// Exact internal-state schema version.
        schema_version: u32,
        /// Content identity of the explicit state artifact.
        artifact: ContentHash,
    },
}

impl InitialStateBinding {
    /// Exact state-schema version used by this binding.
    #[must_use]
    pub const fn schema_version(self) -> u32 {
        match self {
            Self::Zero { schema_version } | Self::Explicit { schema_version, .. } => schema_version,
        }
    }

    /// Explicit state-artifact content identity, when the model does not use
    /// its canonical zero state.
    #[must_use]
    pub const fn artifact(self) -> Option<ContentHash> {
        match self {
            Self::Zero { .. } => None,
            Self::Explicit { artifact, .. } => Some(artifact),
        }
    }

    fn validate_against(self, model: &MaterialModelBinding) -> Result<(), IdentifiabilityError> {
        let schema_version = match self {
            Self::Zero { schema_version } | Self::Explicit { schema_version, .. } => schema_version,
        };
        if schema_version != model.state_schema_version {
            return Err(IdentifiabilityError::VersionMismatch {
                field: "initial state schema",
                expected: model.state_schema_version,
                actual: schema_version,
            });
        }
        match (model.initial_state_policy, self) {
            (InitialStatePolicy::ZeroInternalState, Self::Zero { .. }) => Ok(()),
            (InitialStatePolicy::RequiresDeclaredState, Self::Explicit { artifact, .. })
                if hash_is_nonzero(artifact) =>
            {
                Ok(())
            }
            (InitialStatePolicy::RequiresDeclaredState, Self::Explicit { .. }) => {
                Err(IdentifiabilityError::ZeroIdentity {
                    field: "explicit initial-state artifact",
                })
            }
            (InitialStatePolicy::ZeroInternalState, Self::Explicit { .. }) => {
                Err(IdentifiabilityError::InitialStatePolicy {
                    detail: "model card declares the canonical zero state, but an unrelated explicit state was supplied"
                        .to_string(),
                })
            }
            (InitialStatePolicy::RequiresDeclaredState, Self::Zero { .. }) => {
                Err(IdentifiabilityError::InitialStatePolicy {
                    detail: "model card requires a declared state artifact; zero/default cannot substitute"
                        .to_string(),
                })
            }
        }
    }
}

/// Content-bound spatial frame and orientation convention.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FrameBinding {
    id: ArtifactId,
    transform: ContentHash,
    convention: String,
}

impl FrameBinding {
    /// Construct a spatial frame binding.
    pub fn try_new(
        id: ArtifactId,
        transform: ContentHash,
        convention: impl Into<String>,
    ) -> Result<Self, IdentifiabilityError> {
        let convention = convention.into();
        validate_token(&convention, "frame convention")?;
        if !hash_is_nonzero(transform) {
            return Err(IdentifiabilityError::ZeroIdentity {
                field: "frame transform",
            });
        }
        Ok(Self {
            id,
            transform,
            convention,
        })
    }

    /// Frame identity.
    #[must_use]
    pub const fn id(&self) -> &ArtifactId {
        &self.id
    }

    /// Exact transform/orientation artifact.
    #[must_use]
    pub const fn transform(&self) -> ContentHash {
        self.transform
    }

    /// Orientation/handedness convention.
    #[must_use]
    pub fn convention(&self) -> &str {
        &self.convention
    }
}

/// Exact specimen, process, geometry, and frame binding.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SpecimenBinding {
    id: ArtifactId,
    geometry: ContentHash,
    process: ContentHash,
    preparation: ContentHash,
    frame: FrameBinding,
}

impl SpecimenBinding {
    /// Construct one specimen binding.
    pub fn try_new(
        id: ArtifactId,
        geometry: ContentHash,
        process: ContentHash,
        preparation: ContentHash,
        frame: FrameBinding,
    ) -> Result<Self, IdentifiabilityError> {
        for (field, hash) in [
            ("specimen geometry", geometry),
            ("specimen process", process),
            ("specimen preparation", preparation),
        ] {
            if !hash_is_nonzero(hash) {
                return Err(IdentifiabilityError::ZeroIdentity { field });
            }
        }
        Ok(Self {
            id,
            geometry,
            process,
            preparation,
            frame,
        })
    }

    /// Specimen identity.
    #[must_use]
    pub const fn id(&self) -> &ArtifactId {
        &self.id
    }

    /// Exact specimen-geometry content identity.
    #[must_use]
    pub const fn geometry(&self) -> ContentHash {
        self.geometry
    }

    /// Exact manufacturing/process content identity.
    #[must_use]
    pub const fn process(&self) -> ContentHash {
        self.process
    }

    /// Exact preparation/conditioning content identity.
    #[must_use]
    pub const fn preparation(&self) -> ContentHash {
        self.preparation
    }

    /// Specimen frame.
    #[must_use]
    pub const fn frame(&self) -> &FrameBinding {
        &self.frame
    }
}

/// Exact load/environment/time/refinement protocol.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProtocolBinding {
    id: ArtifactId,
    version: u32,
    state_schema_version: u32,
    refinement_version: u32,
    load_path: ContentHash,
    environment_path: ContentHash,
    time_grid: ContentHash,
    clock: ArtifactId,
}

impl ProtocolBinding {
    /// Construct one exact protocol binding.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        id: ArtifactId,
        version: u32,
        state_schema_version: u32,
        refinement_version: u32,
        load_path: ContentHash,
        environment_path: ContentHash,
        time_grid: ContentHash,
        clock: ArtifactId,
    ) -> Result<Self, IdentifiabilityError> {
        if version == 0 || refinement_version == 0 {
            return Err(IdentifiabilityError::InvalidNumeric {
                field: "protocol/refinement version",
                detail: "version zero is not a published semantics".to_string(),
            });
        }
        for (field, hash) in [
            ("protocol load path", load_path),
            ("protocol environment path", environment_path),
            ("protocol time grid", time_grid),
        ] {
            if !hash_is_nonzero(hash) {
                return Err(IdentifiabilityError::ZeroIdentity { field });
            }
        }
        Ok(Self {
            id,
            version,
            state_schema_version,
            refinement_version,
            load_path,
            environment_path,
            time_grid,
            clock,
        })
    }

    /// Protocol identity.
    #[must_use]
    pub const fn id(&self) -> &ArtifactId {
        &self.id
    }

    /// Protocol semantics version.
    #[must_use]
    pub const fn version(&self) -> u32 {
        self.version
    }

    /// Constitutive state schema expected by this protocol.
    #[must_use]
    pub const fn state_schema_version(&self) -> u32 {
        self.state_schema_version
    }

    /// Mesh/time/refinement policy version.
    #[must_use]
    pub const fn refinement_version(&self) -> u32 {
        self.refinement_version
    }

    /// Exact prescribed-load path identity.
    #[must_use]
    pub const fn load_path(&self) -> ContentHash {
        self.load_path
    }

    /// Exact environmental path identity.
    #[must_use]
    pub const fn environment_path(&self) -> ContentHash {
        self.environment_path
    }

    /// Exact time-grid identity.
    #[must_use]
    pub const fn time_grid(&self) -> ContentHash {
        self.time_grid
    }

    /// Experiment clock identity.
    #[must_use]
    pub const fn clock(&self) -> &ArtifactId {
        &self.clock
    }
}

/// Exact Context-of-Use identity plus the QoI/unit index needed for local
/// law/experiment closure. Construction requires the concrete V&V artifact;
/// callers cannot manufacture the derived index independently.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContextBinding {
    reference: ArtifactRef,
    qoi_units: BTreeMap<QoiId, UnitId>,
}

impl ContextBinding {
    /// Derive a binding from one concrete canonical V&V ContextOfUse.
    pub fn from_vv(context: &ContextOfUse) -> Result<Self, IdentifiabilityError> {
        let hash = context
            .content_hash()
            .map_err(|error| IdentifiabilityError::Vv {
                detail: error.to_string(),
            })?;
        if !hash_is_nonzero(hash) {
            return Err(IdentifiabilityError::ZeroIdentity {
                field: "context-of-use reference",
            });
        }
        let qoi_units = context
            .qois()
            .iter()
            .map(|(qoi, spec)| (qoi.clone(), spec.unit().clone()))
            .collect();
        Ok(Self {
            reference: ArtifactRef::new(ArtifactKind::ContextOfUse, context.id().clone(), hash),
            qoi_units,
        })
    }

    /// Exact context artifact reference.
    #[must_use]
    pub const fn reference(&self) -> &ArtifactRef {
        &self.reference
    }

    /// Exact context QoI-to-declared-unit index.
    #[must_use]
    pub const fn qoi_units(&self) -> &BTreeMap<QoiId, UnitId> {
        &self.qoi_units
    }
}

/// Immutable raw-data, custody, calibration/validation, and blind-holdout
/// lineage derived from concrete fs-evidence V&V artifacts.
#[derive(Clone, PartialEq, Eq)]
pub struct DataLineage {
    experiment: ArtifactRef,
    split: ArtifactRef,
    raw_manifest: ContentHash,
    source_bytes: ContentHash,
    custody_receipt: ContentHash,
    preregistration: ContentHash,
    blind_commitment: ContentHash,
    qois: BTreeSet<QoiId>,
    observation_ids: BTreeSet<ObservationId>,
    row_bindings: BTreeMap<ObservationId, ObservationManifestRow>,
    calibration_ids: BTreeSet<ObservationId>,
    validation_ids: BTreeSet<ObservationId>,
    blind_sources: BTreeMap<ObservationId, ContentHash>,
    parser: ContentHash,
    parser_version: u32,
    preprocessing: ContentHash,
    split_grouping: ArtifactId,
    vv_schema_version: u32,
}

impl fmt::Debug for DataLineage {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter
            .debug_struct("DataLineage")
            .field("experiment_binding", &"<redacted>")
            .field("split_binding", &"<redacted>")
            .field("qoi_count", &self.qois.len())
            .field("observation_count", &self.observation_ids.len())
            .field("partition_counts", &self.partition_counts())
            .field("blind_commitment", &"<redacted>")
            .field("parser_version", &self.parser_version)
            .field("vv_schema_version", &self.vv_schema_version)
            .field("row_bindings", &"<redacted>")
            .field("source_and_custody_hashes", &"<redacted>")
            .finish_non_exhaustive()
    }
}

impl DataLineage {
    fn validate_structural(&self) -> Result<(), IdentifiabilityError> {
        if self.experiment.kind() != ArtifactKind::ExperimentArtifact
            || self.split.kind() != ArtifactKind::CalibrationSplit
        {
            return Err(IdentifiabilityError::Vv {
                detail: "data lineage references have the wrong artifact kinds".to_string(),
            });
        }
        for (field, hash) in [
            ("experiment reference", self.experiment.hash()),
            ("calibration split reference", self.split.hash()),
            ("raw observation manifest", self.raw_manifest),
            ("raw source bytes", self.source_bytes),
            ("data custody receipt", self.custody_receipt),
            ("split preregistration", self.preregistration),
            ("blind-holdout commitment", self.blind_commitment),
            ("observation parser", self.parser),
            ("observation preprocessing pipeline", self.preprocessing),
        ] {
            if !hash_is_nonzero(hash) {
                return Err(IdentifiabilityError::ZeroIdentity { field });
            }
        }
        let blind_ids = self.blind_sources.keys().cloned().collect::<BTreeSet<_>>();
        let row_ids = self.row_bindings.keys().cloned().collect::<BTreeSet<_>>();
        let unique_row_locators = self
            .row_bindings
            .values()
            .map(|row| row.source_ref().locator_identity())
            .collect::<BTreeSet<_>>();
        let row_qois = self
            .row_bindings
            .values()
            .map(|row| row.qoi().clone())
            .collect::<BTreeSet<_>>();
        let canonical_manifest = ObservationManifest::try_new(
            self.row_bindings
                .iter()
                .map(|(id, row)| (id.clone(), row.clone()))
                .collect(),
        )
        .map_err(|error| IdentifiabilityError::Vv {
            detail: format!("data-lineage observation manifest is invalid: {error}"),
        })?;
        let mut partition_union = self.calibration_ids.clone();
        partition_union.extend(self.validation_ids.iter().cloned());
        partition_union.extend(blind_ids.iter().cloned());
        if self.parser_version == 0
            || self.qois.is_empty()
            || self.observation_ids.is_empty()
            || row_ids != self.observation_ids
            || unique_row_locators.len() != self.row_bindings.len()
            || self
                .row_bindings
                .values()
                .any(|row| !hash_is_nonzero(row.locator_hash()))
            || self
                .row_bindings
                .values()
                .any(|row| row.source_ref().dataset_source_bytes_hash() != self.source_bytes)
            || row_qois != self.qois
            || self.calibration_ids.is_empty()
            || self.validation_ids.is_empty()
            || self.blind_sources.is_empty()
            || !self.calibration_ids.is_disjoint(&self.validation_ids)
            || !self.calibration_ids.is_disjoint(&blind_ids)
            || !self.validation_ids.is_disjoint(&blind_ids)
            || partition_union != self.observation_ids
            || self.blind_sources.iter().any(|(id, source)| {
                !hash_is_nonzero(*source)
                    || self
                        .row_bindings
                        .get(id)
                        .is_none_or(|row| row.locator_hash() != *source)
            })
        {
            return Err(IdentifiabilityError::Cardinality {
                field: "data lineage",
                detail: "QoIs/typed row sources/partitions/parser version are inconsistent"
                    .to_string(),
            });
        }
        if canonical_manifest.canonical_hash() != self.raw_manifest {
            return Err(IdentifiabilityError::SourceMismatch {
                field: "raw observation manifest",
            });
        }
        Ok(())
    }

    /// Derive a closed lineage from the concrete canonical artifacts. An
    /// arbitrary caller-built `ArtifactRef` is not accepted as existence or
    /// split-consistency proof.
    pub fn from_vv(
        experiment: &ExperimentArtifact,
        split: &CalibrationSplit,
        parser: ContentHash,
        parser_version: u32,
        preprocessing: ContentHash,
        split_grouping: ArtifactId,
    ) -> Result<Self, IdentifiabilityError> {
        if parser_version == 0 {
            return Err(IdentifiabilityError::InvalidNumeric {
                field: "observation parser version",
                detail: "version zero is not a published parser semantics".to_string(),
            });
        }
        for (field, hash) in [
            ("observation parser", parser),
            ("observation preprocessing pipeline", preprocessing),
        ] {
            if !hash_is_nonzero(hash) {
                return Err(IdentifiabilityError::ZeroIdentity { field });
            }
        }
        let experiment_hash =
            experiment
                .content_hash()
                .map_err(|error| IdentifiabilityError::Vv {
                    detail: error.to_string(),
                })?;
        let experiment_ref = ArtifactRef::new(
            ArtifactKind::ExperimentArtifact,
            experiment.id().clone(),
            experiment_hash,
        );
        if split.experiment() != &experiment_ref {
            return Err(IdentifiabilityError::Vv {
                detail: "calibration split does not bind this exact experiment kind/id/hash"
                    .to_string(),
            });
        }
        let split_hash = split
            .content_hash()
            .map_err(|error| IdentifiabilityError::Vv {
                detail: error.to_string(),
            })?;
        let split_ref = ArtifactRef::new(
            ArtifactKind::CalibrationSplit,
            split.id().clone(),
            split_hash,
        );
        let authenticity = experiment.authenticity();
        for (field, hash) in [
            (
                "raw observation manifest",
                experiment.manifest().canonical_hash(),
            ),
            ("raw source bytes", authenticity.source_bytes_hash()),
            ("data custody receipt", authenticity.custody_receipt_hash()),
            ("split preregistration", split.preregistration_hash()),
            ("blind-holdout commitment", split.blind_commitment()),
        ] {
            if !hash_is_nonzero(hash) {
                return Err(IdentifiabilityError::ZeroIdentity { field });
            }
        }
        let calibration_ids = split.calibration_ids().clone();
        let validation_ids = split.validation_ids().clone();
        let blind_sources = split.blind_sources().clone();
        let mut partition_union = calibration_ids.clone();
        partition_union.extend(validation_ids.iter().cloned());
        partition_union.extend(blind_sources.keys().cloned());
        if partition_union != *experiment.observation_ids() {
            return Err(IdentifiabilityError::Vv {
                detail: "calibration/validation/blind partitions are not the exact experiment manifest row set"
                    .to_string(),
            });
        }
        for (row, source) in &blind_sources {
            if experiment.manifest().locator_hash_of(row) != Some(*source) {
                return Err(IdentifiabilityError::Vv {
                    detail: format!(
                        "blind row {} is not bound to its exact experiment-manifest source",
                        row.as_str()
                    ),
                });
            }
        }
        let lineage = Self {
            experiment: experiment_ref,
            split: split_ref,
            raw_manifest: experiment.manifest().canonical_hash(),
            source_bytes: authenticity.source_bytes_hash(),
            custody_receipt: authenticity.custody_receipt_hash(),
            preregistration: split.preregistration_hash(),
            blind_commitment: split.blind_commitment(),
            qois: experiment.qois().clone(),
            observation_ids: experiment.observation_ids().clone(),
            row_bindings: experiment.manifest().rows().clone(),
            calibration_ids,
            validation_ids,
            blind_sources,
            parser,
            parser_version,
            preprocessing,
            split_grouping,
            vv_schema_version: VV_SCHEMA_VERSION,
        };
        lineage.validate_structural()?;
        Ok(lineage)
    }

    /// Exact experiment artifact reference.
    #[must_use]
    pub const fn experiment(&self) -> &ArtifactRef {
        &self.experiment
    }

    /// Exact calibration/validation/blind split artifact reference.
    #[must_use]
    pub const fn split(&self) -> &ArtifactRef {
        &self.split
    }

    /// Derived raw observation-manifest identity.
    #[must_use]
    pub const fn raw_manifest(&self) -> ContentHash {
        self.raw_manifest
    }

    /// Sealed blind-holdout commitment.
    #[must_use]
    pub const fn blind_commitment(&self) -> ContentHash {
        self.blind_commitment
    }

    /// QoIs present in the exact experiment.
    #[must_use]
    pub const fn qois(&self) -> &BTreeSet<QoiId> {
        &self.qois
    }

    /// Exact row-level semantic bindings used by crate-internal admission and
    /// data-reuse checks. Validation and blind rows deliberately remain
    /// outside the public estimation capability.
    #[must_use]
    pub(crate) const fn row_bindings(&self) -> &BTreeMap<ObservationId, ObservationManifestRow> {
        &self.row_bindings
    }

    /// Digest of the retained raw source bytes.
    #[must_use]
    pub const fn source_bytes(&self) -> ContentHash {
        self.source_bytes
    }

    /// Exact custody-receipt content identity inherited from the experiment.
    #[must_use]
    pub const fn custody_receipt(&self) -> ContentHash {
        self.custody_receipt
    }

    /// Exact split-preregistration content identity.
    #[must_use]
    pub const fn preregistration(&self) -> ContentHash {
        self.preregistration
    }

    /// Exact parser implementation/configuration content identity.
    #[must_use]
    pub const fn parser(&self) -> ContentHash {
        self.parser
    }

    /// Published parser semantics version.
    #[must_use]
    pub const fn parser_version(&self) -> u32 {
        self.parser_version
    }

    /// Exact preprocessing-pipeline content identity.
    #[must_use]
    pub const fn preprocessing(&self) -> ContentHash {
        self.preprocessing
    }

    /// Grouping identity used to construct the calibration/validation/blind
    /// partition.
    #[must_use]
    pub const fn split_grouping(&self) -> &ArtifactId {
        &self.split_grouping
    }

    /// V&V artifact schema under which this lineage was derived.
    #[must_use]
    pub const fn vv_schema_version(&self) -> u32 {
        self.vv_schema_version
    }

    /// Inspect the full semantic binding for an estimation-authorized row.
    /// Validation, blind, and unknown rows return `None`.
    #[must_use]
    pub fn row_binding(&self, row: &ObservationId) -> Option<&ObservationManifestRow> {
        if !self.calibration_ids.contains(row) {
            return None;
        }
        self.row_bindings.get(row)
    }

    /// Rows preregistered for calibration/estimation. Validation and blind
    /// rows are intentionally unavailable through this accessor.
    #[must_use]
    pub const fn calibration_ids(&self) -> &BTreeSet<ObservationId> {
        &self.calibration_ids
    }

    /// Counts of `(calibration, validation, blind holdout)` rows.
    #[must_use]
    pub fn partition_counts(&self) -> (usize, usize, usize) {
        (
            self.calibration_ids.len(),
            self.validation_ids.len(),
            self.blind_sources.len(),
        )
    }
}

fn matrix_get(matrix: &CovarianceMatrix, row: usize, column: usize) -> f64 {
    let (row, column) = if row >= column {
        (row, column)
    } else {
        (column, row)
    };
    matrix.lower_triangle()[row * (row + 1) / 2 + column]
}

fn checked_derivative_dims(output: Dims, input: Dims) -> Option<Dims> {
    let mut result = [0i8; 6];
    for (index, value) in result.iter_mut().enumerate() {
        *value = output.0[index].checked_sub(input.0[index])?;
    }
    Some(Dims(result))
}

fn checked_add_dims(left: Dims, right: Dims) -> Option<Dims> {
    let mut result = [0i8; 6];
    for (index, value) in result.iter_mut().enumerate() {
        *value = left.0[index].checked_add(right.0[index])?;
    }
    Some(Dims(result))
}

fn validate_header_profile(header: &ArtifactHeader) -> Result<(), IdentifiabilityError> {
    validate_token(header.id().as_str(), "study artifact id")?;
    for unit in header.units() {
        validate_token(unit.as_str(), "study unit")?;
    }
    if let SeedDeclaration::NotApplicable { reason } = header.seed() {
        validate_reason(reason, "seed no-claim reason")?;
    }
    if let DeclaredBudget::NotApplicable { reason } = header.accuracy() {
        validate_reason(reason, "accuracy no-claim reason")?;
    }
    for budget in [header.time_ms(), header.memory_bytes()] {
        if let DeclaredBudget::NotApplicable { reason } = budget {
            validate_reason(reason, "resource no-claim reason")?;
        }
    }
    for (component, version) in header.versions() {
        validate_token(component, "study version component")?;
        validate_token(version, "study version value")?;
    }
    for capability in header.capabilities() {
        validate_token(capability, "study capability")?;
    }
    Ok(())
}

struct CanonicalWriter {
    bytes: Vec<u8>,
    error: Option<IdentifiabilityError>,
}

impl CanonicalWriter {
    fn new() -> Self {
        Self {
            bytes: Vec::new(),
            error: None,
        }
    }

    fn raw(&mut self, value: &[u8]) {
        if self.error.is_some() {
            return;
        }
        let Some(requested) = self.bytes.len().checked_add(value.len()) else {
            self.error = Some(IdentifiabilityError::Canonical {
                at: self.bytes.len(),
                detail: "canonical study length overflow".to_string(),
            });
            return;
        };
        if requested > MAX_IDENTIFIABILITY_CANONICAL_BYTES {
            self.error = Some(IdentifiabilityError::Canonical {
                at: self.bytes.len(),
                detail: format!(
                    "canonical study would require {requested} bytes; limit is {MAX_IDENTIFIABILITY_CANONICAL_BYTES}"
                ),
            });
            return;
        }
        // The logical byte ceiling is checked before asking Vec to grow. Vec
        // and the allocator may retain or round additional capacity, so this
        // is deliberately not an exact resident-memory claim; geometric growth
        // keeps repeated small-field encoding amortized rather than quadratic.
        if self.bytes.try_reserve(value.len()).is_err() {
            self.error = Some(IdentifiabilityError::Canonical {
                at: self.bytes.len(),
                detail: "canonical study allocation refused".to_string(),
            });
            return;
        }
        self.bytes.extend_from_slice(value);
    }

    fn byte(&mut self, value: u8) {
        self.raw(&[value]);
    }

    fn u32(&mut self, value: u32) {
        self.raw(&value.to_le_bytes());
    }

    fn u64(&mut self, value: u64) {
        self.raw(&value.to_le_bytes());
    }

    fn f64(&mut self, value: f64) {
        self.u64(canonical_f64(value).to_bits());
    }

    fn count(&mut self, count: usize, field: &'static str) -> Result<(), IdentifiabilityError> {
        if self.error.is_some() {
            return Ok(());
        }
        self.u32(
            u32::try_from(count).map_err(|_| IdentifiabilityError::Cardinality {
                field,
                detail: "count exceeds u32 canonical framing".to_string(),
            })?,
        );
        Ok(())
    }

    fn text(&mut self, value: &str, field: &'static str) -> Result<(), IdentifiabilityError> {
        self.count(value.len(), field)?;
        self.raw(value.as_bytes());
        Ok(())
    }

    fn hash(&mut self, value: ContentHash) {
        self.raw(value.as_bytes());
    }

    fn quantity(&mut self, value: QuantitySpec) {
        self.raw(&value.canonical_bytes());
    }

    fn finish(self) -> Result<Vec<u8>, IdentifiabilityError> {
        if let Some(error) = self.error {
            return Err(error);
        }
        Ok(self.bytes)
    }
}

fn encode_artifact_id(
    writer: &mut CanonicalWriter,
    id: &ArtifactId,
) -> Result<(), IdentifiabilityError> {
    writer.text(id.as_str(), "artifact id")
}

fn encode_qoi_id(writer: &mut CanonicalWriter, id: &QoiId) -> Result<(), IdentifiabilityError> {
    writer.text(id.as_str(), "QoI id")
}

fn encode_observation_row_id(
    writer: &mut CanonicalWriter,
    id: &ObservationId,
) -> Result<(), IdentifiabilityError> {
    writer.text(id.as_str(), "observation row id")
}

fn encode_header(
    writer: &mut CanonicalWriter,
    header: &ArtifactHeader,
    exact: bool,
) -> Result<(), IdentifiabilityError> {
    writer.byte(u8::from(exact));
    if exact {
        encode_artifact_id(writer, header.id())?;
    }
    writer.count(header.units().len(), "header units")?;
    for unit in header.units() {
        writer.text(unit.as_str(), "header unit")?;
    }
    match header.seed() {
        SeedDeclaration::Fixed(seed) => {
            writer.byte(0);
            writer.u64(*seed);
        }
        SeedDeclaration::NotApplicable { reason } => {
            writer.byte(1);
            writer.text(reason, "seed no-claim reason")?;
        }
    }
    match header.accuracy() {
        DeclaredBudget::Limit(value) => {
            writer.byte(0);
            writer.f64(*value);
        }
        DeclaredBudget::NotApplicable { reason } => {
            writer.byte(1);
            writer.text(reason, "accuracy no-claim reason")?;
        }
    }
    for budget in [header.time_ms(), header.memory_bytes()] {
        match budget {
            DeclaredBudget::Limit(value) => {
                writer.byte(0);
                writer.u64(*value);
            }
            DeclaredBudget::NotApplicable { reason } => {
                writer.byte(1);
                writer.text(reason, "resource no-claim reason")?;
            }
        }
    }
    writer.count(header.versions().len(), "header versions")?;
    for (component, version) in header.versions() {
        writer.text(component, "header version component")?;
        writer.text(version, "header version value")?;
    }
    writer.count(header.capabilities().len(), "header capabilities")?;
    for capability in header.capabilities() {
        writer.text(capability, "header capability")?;
    }
    Ok(())
}

fn encode_parameter_domain(writer: &mut CanonicalWriter, domain: ParameterDomain) {
    writer.f64(domain.lo);
    writer.f64(domain.hi);
}

fn encode_prior(
    writer: &mut CanonicalWriter,
    prior: &ParameterPrior,
) -> Result<(), IdentifiabilityError> {
    match prior {
        ParameterPrior::None { reason, version } => {
            writer.byte(0);
            writer.u32(*version);
            writer.text(reason, "prior absence reason")?;
        }
        ParameterPrior::Uniform { domain, version } => {
            writer.byte(1);
            writer.u32(*version);
            encode_parameter_domain(writer, *domain);
        }
        ParameterPrior::Gaussian {
            mean,
            standard_deviation,
            version,
        } => {
            writer.byte(2);
            writer.u32(*version);
            writer.f64(*mean);
            writer.f64(*standard_deviation);
        }
        ParameterPrior::LogNormal {
            log_mean,
            log_standard_deviation,
            reference,
            version,
        } => {
            writer.byte(3);
            writer.u32(*version);
            writer.f64(*log_mean);
            writer.f64(*log_standard_deviation);
            writer.f64(*reference);
        }
    }
    Ok(())
}

fn encode_coordinate(
    writer: &mut CanonicalWriter,
    coordinate: &ParameterCoordinate,
) -> Result<(), IdentifiabilityError> {
    writer.text(coordinate.id.as_str(), "coordinate id")?;
    writer.quantity(coordinate.quantity);
    encode_parameter_domain(writer, coordinate.domain);
    match coordinate.transform {
        CoordinateTransform::Identity => writer.byte(0),
        CoordinateTransform::Affine {
            scale,
            scale_quantity,
            offset,
        } => {
            writer.byte(1);
            writer.f64(scale);
            writer.quantity(scale_quantity);
            writer.f64(offset);
        }
        CoordinateTransform::LogPositive { reference } => {
            writer.byte(2);
            writer.f64(reference);
        }
    }
    Ok(())
}

fn encode_initial_state(writer: &mut CanonicalWriter, state: InitialStateBinding) {
    match state {
        InitialStateBinding::Zero { schema_version } => {
            writer.byte(0);
            writer.u32(schema_version);
        }
        InitialStateBinding::Explicit {
            schema_version,
            artifact,
        } => {
            writer.byte(1);
            writer.u32(schema_version);
            writer.hash(artifact);
        }
    }
}

fn encode_frame(
    writer: &mut CanonicalWriter,
    frame: &FrameBinding,
) -> Result<(), IdentifiabilityError> {
    encode_artifact_id(writer, &frame.id)?;
    writer.hash(frame.transform);
    writer.text(&frame.convention, "frame convention")?;
    Ok(())
}

fn encode_specimen(
    writer: &mut CanonicalWriter,
    specimen: &SpecimenBinding,
) -> Result<(), IdentifiabilityError> {
    encode_artifact_id(writer, &specimen.id)?;
    writer.hash(specimen.geometry);
    writer.hash(specimen.process);
    writer.hash(specimen.preparation);
    encode_frame(writer, &specimen.frame)
}

fn encode_protocol(
    writer: &mut CanonicalWriter,
    protocol: &ProtocolBinding,
) -> Result<(), IdentifiabilityError> {
    encode_artifact_id(writer, &protocol.id)?;
    writer.u32(protocol.version);
    writer.u32(protocol.state_schema_version);
    writer.u32(protocol.refinement_version);
    writer.hash(protocol.load_path);
    writer.hash(protocol.environment_path);
    writer.hash(protocol.time_grid);
    encode_artifact_id(writer, &protocol.clock)
}

struct CanonicalReader<'a> {
    bytes: &'a [u8],
    at: usize,
}

impl<'a> CanonicalReader<'a> {
    fn new(bytes: &'a [u8]) -> Result<Self, IdentifiabilityError> {
        if bytes.len() > MAX_IDENTIFIABILITY_CANONICAL_BYTES {
            return Err(IdentifiabilityError::Canonical {
                at: bytes.len(),
                detail: format!(
                    "input exceeds {MAX_IDENTIFIABILITY_CANONICAL_BYTES} canonical bytes"
                ),
            });
        }
        Ok(Self { bytes, at: 0 })
    }

    fn take(
        &mut self,
        count: usize,
        field: &'static str,
    ) -> Result<&'a [u8], IdentifiabilityError> {
        let end = self
            .at
            .checked_add(count)
            .ok_or_else(|| IdentifiabilityError::Canonical {
                at: self.at,
                detail: format!("{field} length overflows address space"),
            })?;
        let value =
            self.bytes
                .get(self.at..end)
                .ok_or_else(|| IdentifiabilityError::Canonical {
                    at: self.at,
                    detail: format!("truncated {field}"),
                })?;
        self.at = end;
        Ok(value)
    }

    fn byte(&mut self, field: &'static str) -> Result<u8, IdentifiabilityError> {
        Ok(self.take(1, field)?[0])
    }

    fn u32(&mut self, field: &'static str) -> Result<u32, IdentifiabilityError> {
        let bytes: [u8; 4] =
            self.take(4, field)?
                .try_into()
                .map_err(|_| IdentifiabilityError::Canonical {
                    at: self.at,
                    detail: format!("invalid {field} width"),
                })?;
        Ok(u32::from_le_bytes(bytes))
    }

    fn u64(&mut self, field: &'static str) -> Result<u64, IdentifiabilityError> {
        let bytes: [u8; 8] =
            self.take(8, field)?
                .try_into()
                .map_err(|_| IdentifiabilityError::Canonical {
                    at: self.at,
                    detail: format!("invalid {field} width"),
                })?;
        Ok(u64::from_le_bytes(bytes))
    }

    fn f64(&mut self, field: &'static str) -> Result<f64, IdentifiabilityError> {
        Ok(f64::from_bits(self.u64(field)?))
    }

    fn length(
        &mut self,
        maximum: usize,
        field: &'static str,
    ) -> Result<usize, IdentifiabilityError> {
        let value =
            usize::try_from(self.u32(field)?).map_err(|_| IdentifiabilityError::Canonical {
                at: self.at,
                detail: format!("{field} length is not representable"),
            })?;
        if value > maximum {
            return Err(IdentifiabilityError::Canonical {
                at: self.at,
                detail: format!("{field} length {value} exceeds {maximum}"),
            });
        }
        Ok(value)
    }

    fn count(&mut self, field: &'static str) -> Result<usize, IdentifiabilityError> {
        self.length(MAX_IDENTIFIABILITY_ITEMS, field)
    }

    fn text(
        &mut self,
        maximum: usize,
        field: &'static str,
    ) -> Result<String, IdentifiabilityError> {
        let length = self.length(maximum, field)?;
        let at = self.at;
        let bytes = self.take(length, field)?;
        String::from_utf8(bytes.to_vec()).map_err(|_| IdentifiabilityError::Canonical {
            at,
            detail: format!("{field} is not UTF-8"),
        })
    }

    fn token(&mut self, field: &'static str) -> Result<String, IdentifiabilityError> {
        let value = self.text(MAX_IDENTIFIABILITY_ID_BYTES, field)?;
        validate_token(&value, field)?;
        Ok(value)
    }

    fn reason(&mut self, field: &'static str) -> Result<String, IdentifiabilityError> {
        let value = self.text(MAX_IDENTIFIABILITY_TEXT_BYTES, field)?;
        validate_reason(&value, field)?;
        Ok(value)
    }

    fn hash(&mut self, field: &'static str) -> Result<ContentHash, IdentifiabilityError> {
        let bytes: [u8; 32] =
            self.take(32, field)?
                .try_into()
                .map_err(|_| IdentifiabilityError::Canonical {
                    at: self.at,
                    detail: format!("invalid {field} hash width"),
                })?;
        Ok(ContentHash(bytes))
    }

    fn quantity(&mut self, field: &'static str) -> Result<QuantitySpec, IdentifiabilityError> {
        let at = self.at;
        QuantitySpec::from_canonical_bytes(self.take(QUANTITY_SPEC_ENCODED_LEN, field)?).map_err(
            |error| IdentifiabilityError::Canonical {
                at,
                detail: format!("invalid {field} quantity token: {error}"),
            },
        )
    }

    fn expect_byte(
        &mut self,
        expected: u8,
        field: &'static str,
    ) -> Result<(), IdentifiabilityError> {
        let actual = self.byte(field)?;
        if actual == expected {
            Ok(())
        } else {
            Err(IdentifiabilityError::Canonical {
                at: self.at.saturating_sub(1),
                detail: format!("{field} expected tag {expected}, found {actual}"),
            })
        }
    }

    fn finish(self) -> Result<(), IdentifiabilityError> {
        if self.at == self.bytes.len() {
            Ok(())
        } else {
            Err(IdentifiabilityError::Canonical {
                at: self.at,
                detail: format!("{} trailing byte(s)", self.bytes.len() - self.at),
            })
        }
    }
}

fn decode_artifact_id(
    reader: &mut CanonicalReader<'_>,
) -> Result<ArtifactId, IdentifiabilityError> {
    ArtifactId::try_new(reader.token("artifact id")?).map_err(|error| IdentifiabilityError::Vv {
        detail: error.to_string(),
    })
}

fn decode_qoi_id(reader: &mut CanonicalReader<'_>) -> Result<QoiId, IdentifiabilityError> {
    QoiId::try_new(reader.token("QoI id")?).map_err(|error| IdentifiabilityError::Vv {
        detail: error.to_string(),
    })
}

fn decode_observation_row_id(
    reader: &mut CanonicalReader<'_>,
) -> Result<ObservationId, IdentifiabilityError> {
    ObservationId::try_new(reader.token("observation row id")?).map_err(|error| {
        IdentifiabilityError::Vv {
            detail: error.to_string(),
        }
    })
}

fn decode_header(reader: &mut CanonicalReader<'_>) -> Result<ArtifactHeader, IdentifiabilityError> {
    reader.expect_byte(1, "exact header marker")?;
    let id = decode_artifact_id(reader)?;
    let unit_count = reader.count("header units")?;
    let mut units = Vec::with_capacity(unit_count);
    for _ in 0..unit_count {
        units.push(
            UnitId::try_new(reader.token("header unit")?).map_err(|error| {
                IdentifiabilityError::Vv {
                    detail: error.to_string(),
                }
            })?,
        );
    }
    let seed = match reader.byte("seed tag")? {
        0 => SeedDeclaration::Fixed(reader.u64("seed")?),
        1 => SeedDeclaration::NotApplicable {
            reason: reader.reason("seed no-claim reason")?,
        },
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown seed tag {tag}"),
            });
        }
    };
    let accuracy = match reader.byte("accuracy budget tag")? {
        0 => DeclaredBudget::Limit(reader.f64("accuracy budget")?),
        1 => DeclaredBudget::NotApplicable {
            reason: reader.reason("accuracy no-claim reason")?,
        },
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown accuracy budget tag {tag}"),
            });
        }
    };
    let mut resources = Vec::with_capacity(2);
    for _ in 0..2 {
        resources.push(match reader.byte("resource budget tag")? {
            0 => DeclaredBudget::Limit(reader.u64("resource budget")?),
            1 => DeclaredBudget::NotApplicable {
                reason: reader.reason("resource no-claim reason")?,
            },
            tag => {
                return Err(IdentifiabilityError::Canonical {
                    at: reader.at.saturating_sub(1),
                    detail: format!("unknown resource budget tag {tag}"),
                });
            }
        });
    }
    let version_count = reader.count("header versions")?;
    let mut versions = Vec::with_capacity(version_count);
    for _ in 0..version_count {
        versions.push((
            reader.token("header version component")?,
            reader.token("header version value")?,
        ));
    }
    let capability_count = reader.count("header capabilities")?;
    let mut capabilities = Vec::with_capacity(capability_count);
    for _ in 0..capability_count {
        capabilities.push(reader.token("header capability")?);
    }
    ArtifactHeader::try_new(
        id,
        units,
        seed,
        accuracy,
        resources.remove(0),
        resources.remove(0),
        versions,
        capabilities,
    )
    .map_err(|error| IdentifiabilityError::Vv {
        detail: error.to_string(),
    })
}

fn decode_parameter_domain(
    reader: &mut CanonicalReader<'_>,
) -> Result<ParameterDomain, IdentifiabilityError> {
    ParameterDomain::try_new(
        reader.f64("parameter-domain lower bound")?,
        reader.f64("parameter-domain upper bound")?,
    )
}

fn decode_prior(reader: &mut CanonicalReader<'_>) -> Result<ParameterPrior, IdentifiabilityError> {
    Ok(match reader.byte("prior tag")? {
        0 => ParameterPrior::None {
            version: reader.u32("prior version")?,
            reason: reader.reason("prior absence reason")?,
        },
        1 => ParameterPrior::Uniform {
            version: reader.u32("prior version")?,
            domain: decode_parameter_domain(reader)?,
        },
        2 => ParameterPrior::Gaussian {
            version: reader.u32("prior version")?,
            mean: reader.f64("Gaussian prior mean")?,
            standard_deviation: reader.f64("Gaussian prior standard deviation")?,
        },
        3 => ParameterPrior::LogNormal {
            version: reader.u32("prior version")?,
            log_mean: reader.f64("log-normal prior mean")?,
            log_standard_deviation: reader.f64("log-normal prior standard deviation")?,
            reference: reader.f64("log-normal prior reference")?,
        },
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown prior tag {tag}"),
            });
        }
    })
}

fn decode_coordinate(
    reader: &mut CanonicalReader<'_>,
) -> Result<ParameterCoordinate, IdentifiabilityError> {
    let id = CoordinateId::try_new(reader.token("coordinate id")?)?;
    let quantity = reader.quantity("coordinate")?;
    let domain = decode_parameter_domain(reader)?;
    let transform = match reader.byte("coordinate transform tag")? {
        0 => CoordinateTransform::Identity,
        1 => CoordinateTransform::Affine {
            scale: reader.f64("affine scale")?,
            scale_quantity: reader.quantity("affine scale")?,
            offset: reader.f64("affine offset")?,
        },
        2 => CoordinateTransform::LogPositive {
            reference: reader.f64("log reference")?,
        },
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown coordinate transform tag {tag}"),
            });
        }
    };
    ParameterCoordinate::try_new(id, quantity, domain, transform)
}

fn decode_initial_state(
    reader: &mut CanonicalReader<'_>,
) -> Result<InitialStateBinding, IdentifiabilityError> {
    match reader.byte("initial state tag")? {
        0 => Ok(InitialStateBinding::Zero {
            schema_version: reader.u32("initial state schema version")?,
        }),
        1 => Ok(InitialStateBinding::Explicit {
            schema_version: reader.u32("initial state schema version")?,
            artifact: reader.hash("initial state artifact")?,
        }),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown initial state tag {tag}"),
        }),
    }
}

fn decode_frame(reader: &mut CanonicalReader<'_>) -> Result<FrameBinding, IdentifiabilityError> {
    FrameBinding::try_new(
        decode_artifact_id(reader)?,
        reader.hash("frame transform")?,
        reader.token("frame convention")?,
    )
}

fn decode_specimen(
    reader: &mut CanonicalReader<'_>,
) -> Result<SpecimenBinding, IdentifiabilityError> {
    SpecimenBinding::try_new(
        decode_artifact_id(reader)?,
        reader.hash("specimen geometry")?,
        reader.hash("specimen process")?,
        reader.hash("specimen preparation")?,
        decode_frame(reader)?,
    )
}

fn decode_protocol(
    reader: &mut CanonicalReader<'_>,
) -> Result<ProtocolBinding, IdentifiabilityError> {
    ProtocolBinding::try_new(
        decode_artifact_id(reader)?,
        reader.u32("protocol version")?,
        reader.u32("protocol state schema version")?,
        reader.u32("refinement version")?,
        reader.hash("load path")?,
        reader.hash("environment path")?,
        reader.hash("time grid")?,
        decode_artifact_id(reader)?,
    )
}

#[cfg(test)]
mod canonical_writer_resource_tests {
    use super::*;

    #[test]
    fn canonical_writer_bounds_logical_transport_and_refuses_before_reserve() {
        const BLOCK_BYTES: usize = 4096;
        let block = [0xa5; BLOCK_BYTES];
        let block_count = MAX_IDENTIFIABILITY_CANONICAL_BYTES / BLOCK_BYTES;

        let mut exact = CanonicalWriter::new();
        for _ in 0..block_count {
            exact.raw(&block);
        }
        assert_eq!(
            exact
                .finish()
                .expect("the exact transport limit admits")
                .len(),
            MAX_IDENTIFIABILITY_CANONICAL_BYTES,
        );

        let mut oversized = CanonicalWriter::new();
        for _ in 0..block_count {
            oversized.raw(&block);
        }
        let capacity_at_limit = oversized.bytes.capacity();
        oversized.byte(0xff);
        assert_eq!(
            oversized.bytes.len(),
            MAX_IDENTIFIABILITY_CANONICAL_BYTES,
            "the refusing append must not grow the logical transport",
        );
        assert_eq!(
            oversized.bytes.capacity(),
            capacity_at_limit,
            "the over-limit append must refuse before asking Vec for more capacity",
        );
        oversized.raw(&block);
        assert_eq!(
            oversized.bytes.len(),
            MAX_IDENTIFIABILITY_CANONICAL_BYTES,
            "writes after the first refusal must be no-ops",
        );
        assert!(matches!(
            oversized.finish(),
            Err(IdentifiabilityError::Canonical {
                at: MAX_IDENTIFIABILITY_CANONICAL_BYTES,
                ..
            })
        ));
    }
}
