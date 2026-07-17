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
//! The older single-case schema remains crate-private below for design
//! archaeology. Its `StudySpecId` and `PhysicalStudyId` cannot mint authority
//! in the current multi-case chain.
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
    ObservationManifest, ObservationManifestRow, ObservationSourceRef, QoiId, SeedDeclaration,
    UnitId, VV_SCHEMA_VERSION,
};
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, MATDB_SCHEMA_VERSION, MaterialCard,
};
use fs_qty::{Dims, QUANTITY_SPEC_ENCODED_LEN, QuantitySpec};

/// Retained binary/canonical semantics for the crate-private, single-case
/// [`IdentifiabilityStudySpec`] prototype.
pub(crate) const LEGACY_IDENTIFIABILITY_SCHEMA_VERSION: u32 = 2;
/// Maximum identifier or short-reason byte length.
pub const MAX_IDENTIFIABILITY_ID_BYTES: usize = 256;
/// Maximum long diagnostic/reason byte length.
pub const MAX_IDENTIFIABILITY_TEXT_BYTES: usize = 16 * 1024;
/// Maximum rows in any parameter/observation/path/gauge collection.
pub const MAX_IDENTIFIABILITY_ITEMS: usize = 4096;
/// Maximum canonical study bytes accepted or emitted.
pub const MAX_IDENTIFIABILITY_CANONICAL_BYTES: usize = 4 * 1024 * 1024;

const SPEC_DOMAIN: &str = "org.frankensim.fs-material.identifiability-spec.v2";
const PHYSICAL_DOMAIN: &str = "org.frankensim.fs-material.identifiability-physical.v2";
const CANONICAL_MAGIC: &[u8] = b"fs-material-identifiability-study\0";

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

/// Exact replay identity of one admitted study specification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct StudySpecId(ContentHash);

impl StudySpecId {
    /// Underlying domain-separated digest.
    #[must_use]
    pub const fn digest(self) -> ContentHash {
        self.0
    }
}

/// Reparameterization-quotient identity of one admitted physical study.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub(crate) struct PhysicalStudyId(ContentHash);

impl PhysicalStudyId {
    /// Underlying domain-separated digest.
    #[must_use]
    pub const fn digest(self) -> ContentHash {
        self.0
    }
}

/// Deterministic refusal at the identifiability-study boundary.
#[derive(Debug, Clone, PartialEq)]
pub enum IdentifiabilityError {
    /// Retained schema is stale or from the future.
    UnsupportedSchemaVersion { declared: u32, supported: u32 },
    /// A bounded identifier/reason was malformed.
    InvalidText { field: &'static str, detail: String },
    /// A required content identity was the all-zero sentinel.
    ZeroIdentity { field: &'static str },
    /// A count was empty, oversized, or inconsistent.
    Cardinality { field: &'static str, detail: String },
    /// Two rows claimed the same identity.
    Duplicate { field: &'static str, id: String },
    /// A numeric interval/prior/transform was malformed.
    InvalidNumeric { field: &'static str, detail: String },
    /// Exact version pins disagree.
    VersionMismatch {
        field: &'static str,
        expected: u32,
        actual: u32,
    },
    /// One reference points outside the admitted closed graph.
    UnknownReference { field: &'static str, id: String },
    /// A model card is not an exact member of the material card.
    ModelNotInMaterialCard,
    /// A model/state initialization policy was violated.
    InitialStatePolicy { detail: String },
    /// An estimated parameter has no observation path and no honest refusal.
    DisconnectedEstimatedParameter { parameter: ParameterRoleId },
    /// A nuisance parameter is not calibrated by this exact split.
    NuisanceCalibration { parameter: ParameterRoleId },
    /// A declared gauge quotient is structurally invalid.
    InvalidGauge { gauge: GaugeClassId, detail: String },
    /// Covariance order/dimension disagrees with the observation schema.
    Covariance { detail: String },
    /// V&V artifact construction or canonicalization failed.
    Vv { detail: String },
    /// Material-card validation or identity failed.
    Material { detail: String },
    /// Canonical claims disagree with caller-resolved source artifacts.
    SourceMismatch { field: &'static str },
    /// Canonical transport is malformed or exceeds its public cap.
    Canonical { at: usize, detail: String },
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
            Self::InvalidText { field, detail } => write!(f, "invalid {field}: {detail}"),
            Self::ZeroIdentity { field } => write!(f, "{field} uses the all-zero identity"),
            Self::Cardinality { field, detail } => {
                write!(f, "invalid {field} cardinality: {detail}")
            }
            Self::Duplicate { field, id } => write!(f, "duplicate {field} identity {id:?}"),
            Self::InvalidNumeric { field, detail } => write!(f, "invalid {field}: {detail}"),
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
    None { reason: String, version: u32 },
    /// Uniform prior over a canonical physical interval.
    Uniform {
        domain: ParameterDomain,
        version: u32,
    },
    /// Gaussian density in coherent SI units, normalized after conditioning
    /// on the enclosing [`ParameterSpec`] physical domain.
    Gaussian {
        mean: f64,
        standard_deviation: f64,
        version: u32,
    },
    /// Log-normal density over a positive dimensional parameter, normalized
    /// after conditioning on the enclosing physical domain. `reference`
    /// supplies the coherent-SI scale used by the logarithm.
    LogNormal {
        log_mean: f64,
        log_standard_deviation: f64,
        reference: f64,
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
        scale: f64,
        scale_quantity: QuantitySpec,
        offset: f64,
    },
    /// `physical = reference * exp(coordinate)`, using deterministic fs-math.
    LogPositive { reference: f64 },
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
    /// full-domain bijectivity, and physical-domain coverage when the chart is
    /// admitted: [`ParameterSpec::try_new`] does so for the retained prototype,
    /// and [`IdentifiabilityExecutionPlan::try_new`] does so for the current
    /// authority-separated execution stage.
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

/// How one canonical law parameter participates in inference.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParameterClass {
    /// Decision-facing parameter to estimate.
    Target,
    /// Joint nuisance parameter calibrated by the exact study split.
    Nuisance { calibration: ArtifactRef },
    /// Fixed parameter with an exact source artifact.
    Fixed { source: ContentHash },
}

/// Semantic owner of a parameter. Only constitutive-model parameters are
/// required to appear in the exact fs-matdb parameter roster; the remaining
/// owners make state, instrument, discrepancy, and controlled-input nuisance
/// variables explicit instead of smuggling them into a flat vector.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParameterOwner {
    /// Parameter belongs to the constitutive model card.
    ConstitutiveModel,
    /// Parameter describes an initial/internal state coordinate.
    InitialState,
    /// Parameter belongs to an instrument or observation operator.
    Instrument,
    /// Parameter belongs to the model-discrepancy family.
    Discrepancy,
    /// Parameter is a controlled protocol input with uncertainty.
    ControlledInput,
}

/// Population/realization scope of one parameter.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParameterScope {
    /// Shared across the complete study population.
    Global,
    /// Shared only within one material lot/batch.
    MaterialLot { lot: ArtifactId },
    /// Specific to the bound specimen.
    Specimen { specimen: ArtifactId },
    /// Spatially varying field described by a content-addressed basis/mesh.
    Field { support: ContentHash },
    /// One level of a declared hierarchical population model.
    Hierarchical { population: ArtifactId, level: u32 },
}

/// Honest structural observability declaration at schema time.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ParameterObservability {
    /// At least one declared observation path exists; identifiability itself
    /// remains to be assessed.
    Candidate,
    /// The parameter is deliberately retained as unidentifiable rather than
    /// silently dropped or falsely estimated.
    ExplicitlyUnidentifiable {
        reason: String,
        witness: ContentHash,
    },
    /// Parameter is fixed and is not an inference degree of freedom.
    NotEstimated { reason: String },
}

/// One canonical physical parameter plus its exact estimation coordinate.
#[derive(Debug, Clone, PartialEq)]
pub struct ParameterSpec {
    role: ParameterRoleId,
    quantity: QuantitySpec,
    domain: ParameterDomain,
    prior: ParameterPrior,
    coordinate: ParameterCoordinate,
    owner: ParameterOwner,
    scope: ParameterScope,
    class: ParameterClass,
    observability: ParameterObservability,
}

impl ParameterSpec {
    /// Build a parameter row. Cross-row/model/path obligations are checked by
    /// [`IdentifiabilityStudySpec::try_new`].
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        role: ParameterRoleId,
        quantity: QuantitySpec,
        domain: ParameterDomain,
        mut prior: ParameterPrior,
        coordinate: ParameterCoordinate,
        owner: ParameterOwner,
        scope: ParameterScope,
        class: ParameterClass,
        observability: ParameterObservability,
    ) -> Result<Self, IdentifiabilityError> {
        prior.validate_against(domain)?;
        if !matches!(&class, ParameterClass::Fixed { .. }) && domain.is_degenerate() {
            return Err(IdentifiabilityError::InvalidNumeric {
                field: "estimated parameter domain",
                detail: format!("parameter {role} needs a non-degenerate domain"),
            });
        }
        if let ParameterClass::Fixed { source } = &class {
            if !hash_is_nonzero(*source) {
                return Err(IdentifiabilityError::ZeroIdentity {
                    field: "fixed-parameter source",
                });
            }
            if !domain.is_degenerate()
                || !matches!(&prior, ParameterPrior::None { .. })
                || !matches!(&observability, ParameterObservability::NotEstimated { .. })
            {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "fixed parameter",
                    detail: "fixed parameters require a degenerate exact-value domain, no prior, and NotEstimated semantics"
                        .to_string(),
                });
            }
        }
        match &scope {
            ParameterScope::Field { support } if !hash_is_nonzero(*support) => {
                return Err(IdentifiabilityError::ZeroIdentity {
                    field: "parameter field support",
                });
            }
            ParameterScope::Hierarchical { level: 0, .. } => {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "hierarchical parameter level",
                    detail: "level zero is reserved for the global population".to_string(),
                });
            }
            _ => {}
        }
        if let ParameterObservability::ExplicitlyUnidentifiable { reason, witness } = &observability
        {
            validate_reason(reason, "unidentifiable reason")?;
            if !hash_is_nonzero(*witness) {
                return Err(IdentifiabilityError::ZeroIdentity {
                    field: "unidentifiable witness",
                });
            }
        }
        if let ParameterObservability::NotEstimated { reason } = &observability {
            validate_reason(reason, "not-estimated reason")?;
            if !matches!(&class, ParameterClass::Fixed { .. }) {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "parameter observability",
                    detail: "only fixed parameters may use NotEstimated in legacy schema v2"
                        .to_string(),
                });
            }
        }
        match coordinate.transform {
            CoordinateTransform::Identity => {
                if coordinate.quantity != quantity {
                    return Err(IdentifiabilityError::InvalidNumeric {
                        field: "coordinate quantity",
                        detail: format!(
                            "identity coordinate for {role} must retain the exact physical QuantitySpec"
                        ),
                    });
                }
            }
            CoordinateTransform::Affine { scale_quantity, .. } => {
                let mapped_dims =
                    checked_add_dims(coordinate.quantity.dims(), scale_quantity.dims())
                        .ok_or_else(|| IdentifiabilityError::InvalidNumeric {
                            field: "affine coordinate scale quantity",
                            detail: format!("dimension exponents overflow for {role}"),
                        })?;
                if mapped_dims != quantity.dims()
                    || (scale_quantity.dims() == Dims([0; 6])
                        && coordinate.quantity.dims() == quantity.dims()
                        && coordinate.quantity != quantity)
                {
                    return Err(IdentifiabilityError::InvalidNumeric {
                        field: "affine coordinate quantity",
                        detail: format!(
                            "typed scale times coordinate for {role} must produce the physical quantity without a semantic-kind alias"
                        ),
                    });
                }
            }
            CoordinateTransform::LogPositive { .. } => {
                if coordinate.quantity.dims() != Dims([0; 6]) || domain.lo <= 0.0 {
                    return Err(IdentifiabilityError::InvalidNumeric {
                        field: "log coordinate quantity",
                        detail: format!(
                            "log coordinate for {role} must be dimensionless and map to a positive domain"
                        ),
                    });
                }
            }
        }
        let mapped = coordinate.transform.mapped_domain(coordinate.domain)?;
        if !same_f64(mapped.lo, domain.lo) || !same_f64(mapped.hi, domain.hi) {
            return Err(IdentifiabilityError::InvalidNumeric {
                field: "coordinate-domain image",
                detail: format!(
                    "coordinate {} does not map exactly onto canonical domain [{:?}, {:?}]",
                    coordinate.id, domain.lo, domain.hi
                ),
            });
        }
        Ok(Self {
            role,
            quantity,
            domain,
            prior,
            coordinate,
            owner,
            scope,
            class,
            observability,
        })
    }

    /// Canonical law-parameter role.
    #[must_use]
    pub const fn role(&self) -> &ParameterRoleId {
        &self.role
    }

    /// Canonical physical quantity descriptor.
    #[must_use]
    pub const fn quantity(&self) -> QuantitySpec {
        self.quantity
    }

    /// Canonical physical domain.
    #[must_use]
    pub const fn domain(&self) -> ParameterDomain {
        self.domain
    }

    /// Canonical physical prior.
    #[must_use]
    pub const fn prior(&self) -> &ParameterPrior {
        &self.prior
    }

    /// Exact optimizer coordinate.
    #[must_use]
    pub const fn coordinate(&self) -> &ParameterCoordinate {
        &self.coordinate
    }

    /// Semantic owner of the parameter.
    #[must_use]
    pub const fn owner(&self) -> ParameterOwner {
        self.owner
    }

    /// Population/realization scope.
    #[must_use]
    pub const fn scope(&self) -> &ParameterScope {
        &self.scope
    }

    /// Inference role.
    #[must_use]
    pub const fn class(&self) -> &ParameterClass {
        &self.class
    }

    /// Structural observability declaration.
    #[must_use]
    pub const fn observability(&self) -> &ParameterObservability {
        &self.observability
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
    fn validate_structural(&self) -> Result<(), IdentifiabilityError> {
        for (field, hash) in [
            ("material card", self.material_card),
            ("constitutive model card", self.model_card),
            ("canonical parameter block", self.parameter_block),
            ("constitutive graph binding", self.graph),
        ] {
            if !hash_is_nonzero(hash) {
                return Err(IdentifiabilityError::ZeroIdentity { field });
            }
        }
        validate_token(&self.law.0, "constitutive law id")?;
        if self.law_version == 0
            || self.matdb_schema_version == 0
            || self.parameter_roster.is_empty()
            || self.parameter_roster.len() > MAX_IDENTIFIABILITY_ITEMS
        {
            return Err(IdentifiabilityError::Cardinality {
                field: "material/model binding",
                detail:
                    "law/card versions must be published and the parameter roster bounded/nonempty"
                        .to_string(),
            });
        }
        Ok(())
    }

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
    Zero { schema_version: u32 },
    /// An explicit state artifact is required and retained.
    Explicit {
        schema_version: u32,
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

    fn validate_structural(&self) -> Result<(), IdentifiabilityError> {
        if self.reference.kind() != ArtifactKind::ContextOfUse {
            return Err(IdentifiabilityError::Vv {
                detail: "context binding references the wrong artifact kind".to_string(),
            });
        }
        if !hash_is_nonzero(self.reference.hash()) {
            return Err(IdentifiabilityError::ZeroIdentity {
                field: "context-of-use reference",
            });
        }
        if self.qoi_units.is_empty() || self.qoi_units.len() > MAX_IDENTIFIABILITY_ITEMS {
            return Err(IdentifiabilityError::Cardinality {
                field: "context QoIs",
                detail: "context needs a bounded nonempty QoI/unit index".to_string(),
            });
        }
        Ok(())
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

    /// All raw observation row identities retained by the experiment manifest.
    ///
    /// This is crate-private because the set includes validation and sealed
    /// blind-holdout identities. Public estimation APIs expose only
    /// [`Self::calibration_ids`] and calibration-authorized row bindings.
    #[must_use]
    pub(crate) const fn observation_ids(&self) -> &BTreeSet<ObservationId> {
        &self.observation_ids
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

/// Exact sensor/observation-operator semantics for one channel.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SensorBinding {
    device: ArtifactId,
    channel: ArtifactId,
    model: ContentHash,
    model_version: u32,
    calibration_certificate: ContentHash,
    transfer_function: ContentHash,
    filter: ContentHash,
    spatial_support: ContentHash,
    clock: ArtifactId,
    delay_nanoseconds: i64,
    anti_aliasing: ContentHash,
}

impl SensorBinding {
    /// Construct a fully versioned sensor channel.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        device: ArtifactId,
        channel: ArtifactId,
        model: ContentHash,
        model_version: u32,
        calibration_certificate: ContentHash,
        transfer_function: ContentHash,
        filter: ContentHash,
        spatial_support: ContentHash,
        clock: ArtifactId,
        delay_nanoseconds: i64,
        anti_aliasing: ContentHash,
    ) -> Result<Self, IdentifiabilityError> {
        if model_version == 0 {
            return Err(IdentifiabilityError::InvalidNumeric {
                field: "sensor model version",
                detail: "version zero is not a published sensor semantics".to_string(),
            });
        }
        for (field, hash) in [
            ("sensor model", model),
            ("sensor calibration certificate", calibration_certificate),
            ("sensor transfer function", transfer_function),
            ("sensor filter", filter),
            ("sensor spatial support", spatial_support),
            ("sensor anti-aliasing policy", anti_aliasing),
        ] {
            if !hash_is_nonzero(hash) {
                return Err(IdentifiabilityError::ZeroIdentity { field });
            }
        }
        Ok(Self {
            device,
            channel,
            model,
            model_version,
            calibration_certificate,
            transfer_function,
            filter,
            spatial_support,
            clock,
            delay_nanoseconds,
            anti_aliasing,
        })
    }

    /// Sensor clock identity.
    #[must_use]
    pub const fn clock(&self) -> &ArtifactId {
        &self.clock
    }

    /// Observation-operator semantics version.
    #[must_use]
    pub const fn model_version(&self) -> u32 {
        self.model_version
    }
}

/// Missingness/censoring/dropout semantics. Absence is explicit rather than
/// silently interpreted as complete data.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MissingnessModel {
    /// Independently supported claim that the retained rows are complete.
    Complete { evidence: ContentHash },
    /// Explicit modeled missingness/censoring mechanism.
    Modeled { artifact: ContentHash, version: u32 },
    /// No missingness claim is available.
    Unknown { reason: String },
}

impl MissingnessModel {
    fn validate(&self) -> Result<(), IdentifiabilityError> {
        match self {
            Self::Complete { evidence } if hash_is_nonzero(*evidence) => Ok(()),
            Self::Modeled { artifact, version } if hash_is_nonzero(*artifact) && *version > 0 => {
                Ok(())
            }
            Self::Unknown { reason } => validate_reason(reason, "missingness no-claim reason"),
            Self::Complete { .. } => Err(IdentifiabilityError::ZeroIdentity {
                field: "data-completeness evidence",
            }),
            Self::Modeled { .. } => Err(IdentifiabilityError::InvalidNumeric {
                field: "missingness model",
                detail: "model needs a nonzero artifact and positive version".to_string(),
            }),
        }
    }
}

/// Marginal sensor/noise semantics in the channel's coherent SI units.
#[derive(Debug, Clone, PartialEq)]
pub enum NoiseModel {
    /// Deterministic but non-rigorous absolute error floor.
    Bounded { half_width: f64 },
    /// Gaussian marginal standard deviation.
    Gaussian { standard_deviation: f64 },
    /// Student-t marginal with finite degrees of freedom.
    StudentT { scale: f64, degrees_of_freedom: f64 },
    /// Empirical marginal distribution artifact.
    Empirical {
        artifact: ContentHash,
        version: u32,
        reference_scale: f64,
    },
    /// Noise is not characterized; practical identifiability cannot be
    /// upgraded through this channel.
    Unknown { reason: String },
}

impl NoiseModel {
    fn validate(&mut self) -> Result<(), IdentifiabilityError> {
        match self {
            Self::Bounded { half_width } if half_width.is_finite() && *half_width >= 0.0 => {
                *half_width = canonical_f64(*half_width);
                Ok(())
            }
            Self::Gaussian { standard_deviation }
                if standard_deviation.is_finite() && *standard_deviation > 0.0 =>
            {
                Ok(())
            }
            Self::StudentT {
                scale,
                degrees_of_freedom,
            } if scale.is_finite()
                && *scale > 0.0
                && degrees_of_freedom.is_finite()
                && *degrees_of_freedom > 0.0 =>
            {
                Ok(())
            }
            Self::Empirical {
                artifact,
                version,
                reference_scale,
            } if hash_is_nonzero(*artifact)
                && *version > 0
                && reference_scale.is_finite()
                && *reference_scale > 0.0 =>
            {
                Ok(())
            }
            Self::Unknown { reason } => validate_reason(reason, "noise no-claim reason"),
            _ => Err(IdentifiabilityError::InvalidNumeric {
                field: "noise model",
                detail: "noise scales/parameters must be finite and physically admissible"
                    .to_string(),
            }),
        }
    }
}

/// One measured observable with exact physical, sensor, protocol, and raw-row
/// semantics.
#[derive(Debug, Clone, PartialEq)]
pub struct ObservationSpec {
    id: ObservationChannelId,
    qoi: QoiId,
    quantity: QuantitySpec,
    frame: FrameBinding,
    graph_node: String,
    graph_port: String,
    operator: ContentHash,
    operator_version: u32,
    aggregation: ContentHash,
    sensor: SensorBinding,
    noise: NoiseModel,
    missingness: MissingnessModel,
    saturation: Option<ParameterDomain>,
    protocol_version: u32,
    refinement_version: u32,
    source_rows: BTreeSet<ObservationId>,
}

impl ObservationSpec {
    /// Construct one measured-observable channel.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        id: ObservationChannelId,
        qoi: QoiId,
        quantity: QuantitySpec,
        frame: FrameBinding,
        graph_node: impl Into<String>,
        graph_port: impl Into<String>,
        operator: ContentHash,
        operator_version: u32,
        aggregation: ContentHash,
        sensor: SensorBinding,
        mut noise: NoiseModel,
        missingness: MissingnessModel,
        saturation: Option<ParameterDomain>,
        protocol_version: u32,
        refinement_version: u32,
        source_rows: Vec<ObservationId>,
    ) -> Result<Self, IdentifiabilityError> {
        let graph_node = graph_node.into();
        let graph_port = graph_port.into();
        validate_token(&graph_node, "observation graph node")?;
        validate_token(&graph_port, "observation graph port")?;
        if !hash_is_nonzero(operator) || !hash_is_nonzero(aggregation) {
            return Err(IdentifiabilityError::ZeroIdentity {
                field: "observation operator/aggregation",
            });
        }
        if operator_version == 0 || protocol_version == 0 || refinement_version == 0 {
            return Err(IdentifiabilityError::InvalidNumeric {
                field: "observation version",
                detail: "operator, protocol, and refinement versions must be positive".to_string(),
            });
        }
        noise.validate()?;
        missingness.validate()?;
        if source_rows.is_empty() || source_rows.len() > MAX_IDENTIFIABILITY_ITEMS {
            return Err(IdentifiabilityError::Cardinality {
                field: "observation source rows",
                detail: "each channel needs a bounded nonempty raw-row set".to_string(),
            });
        }
        let source_count = source_rows.len();
        let source_rows = source_rows.into_iter().collect::<BTreeSet<_>>();
        if source_rows.len() != source_count {
            return Err(IdentifiabilityError::Duplicate {
                field: "observation source row",
                id: id.to_string(),
            });
        }
        Ok(Self {
            id,
            qoi,
            quantity,
            frame,
            graph_node,
            graph_port,
            operator,
            operator_version,
            aggregation,
            sensor,
            noise,
            missingness,
            saturation,
            protocol_version,
            refinement_version,
            source_rows,
        })
    }

    /// Observation-channel identity.
    #[must_use]
    pub const fn id(&self) -> &ObservationChannelId {
        &self.id
    }

    /// V&V quantity of interest observed by the channel.
    #[must_use]
    pub const fn qoi(&self) -> &QoiId {
        &self.qoi
    }

    /// Exact semantic quantity descriptor.
    #[must_use]
    pub const fn quantity(&self) -> QuantitySpec {
        self.quantity
    }

    /// Raw source-row identities consumed by this channel.
    #[must_use]
    pub const fn source_rows(&self) -> &BTreeSet<ObservationId> {
        &self.source_rows
    }
}

/// Exact observation-distribution functional influenced by a parameter.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum InfluenceFunctional {
    /// Modeled mean of the primary observation channel.
    Mean,
    /// Variance of the primary observation channel.
    Variance,
    /// Covariance of the primary channel with a second named channel.
    Covariance { other: ObservationChannelId },
    /// Dimensionless censoring/saturation probability.
    CensoringProbability,
    /// Dimensionless dropout/missingness probability.
    MissingnessProbability,
}

/// Whether influence reaches the functional directly or through a retained
/// hidden/internal-state trajectory.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum InfluenceRoute {
    /// Direct parameter-to-functional route.
    Direct,
    /// Route mediated by the bound constitutive state evolution.
    StateMediated,
}

/// Evidence state of one parameter-to-observation path. Connectivity is kept
/// distinct from a nonzero sensitivity or observability theorem.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum InfluenceStatus {
    /// Schema-level reachability only; no nonzero claim.
    DeclaredConnectivity,
    /// Symbolic nonzero result from an external checker.
    SymbolicallyNonzero { receipt: ContentHash },
    /// Numerical nonzero witness at a named realization.
    NumericallyWitnessed { receipt: ContentHash },
    /// Constructive proof/witness that this path is zero.
    ProvenZero { witness: ContentHash },
    /// Path semantics are unresolved.
    Unknown { reason: String },
}

impl InfluenceStatus {
    fn validate(&self) -> Result<(), IdentifiabilityError> {
        match self {
            Self::DeclaredConnectivity => Ok(()),
            Self::SymbolicallyNonzero { receipt } | Self::NumericallyWitnessed { receipt }
                if hash_is_nonzero(*receipt) =>
            {
                Ok(())
            }
            Self::ProvenZero { witness } if hash_is_nonzero(*witness) => Ok(()),
            Self::Unknown { reason } => validate_reason(reason, "observation-path unknown reason"),
            _ => Err(IdentifiabilityError::ZeroIdentity {
                field: "observation-path evidence",
            }),
        }
    }

    fn can_support_candidate(&self) -> bool {
        !matches!(self, Self::ProvenZero { .. })
    }
}

/// Complete declared route from one canonical parameter role to one measured
/// observation distribution component.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ObservationPath {
    parameter: ParameterRoleId,
    observation: ObservationChannelId,
    functional: InfluenceFunctional,
    route: InfluenceRoute,
    graph_path: ContentHash,
    derivative_quantity: QuantitySpec,
    status: InfluenceStatus,
}

impl ObservationPath {
    /// Construct one path row. Endpoint and derivative-unit closure are checked
    /// at whole-study admission.
    pub fn try_new(
        parameter: ParameterRoleId,
        observation: ObservationChannelId,
        functional: InfluenceFunctional,
        route: InfluenceRoute,
        graph_path: ContentHash,
        derivative_quantity: QuantitySpec,
        status: InfluenceStatus,
    ) -> Result<Self, IdentifiabilityError> {
        if !hash_is_nonzero(graph_path) {
            return Err(IdentifiabilityError::ZeroIdentity {
                field: "parameter-observation graph path",
            });
        }
        status.validate()?;
        Ok(Self {
            parameter,
            observation,
            functional,
            route,
            graph_path,
            derivative_quantity,
            status,
        })
    }

    /// Canonical parameter endpoint.
    #[must_use]
    pub const fn parameter(&self) -> &ParameterRoleId {
        &self.parameter
    }

    /// Observation endpoint.
    #[must_use]
    pub const fn observation(&self) -> &ObservationChannelId {
        &self.observation
    }

    /// Distribution functional affected by this route.
    #[must_use]
    pub const fn functional(&self) -> &InfluenceFunctional {
        &self.functional
    }

    /// Direct or hidden-state-mediated route class.
    #[must_use]
    pub const fn route(&self) -> InfluenceRoute {
        self.route
    }

    /// Declared derivative quantity for the functional with respect to the
    /// canonical physical parameter.
    #[must_use]
    pub const fn derivative_quantity(&self) -> QuantitySpec {
        self.derivative_quantity
    }

    /// Current evidence status of the route.
    #[must_use]
    pub const fn status(&self) -> &InfluenceStatus {
        &self.status
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

/// Unit-safe cross-channel dependence: marginal scales live in each
/// [`NoiseModel`], while this matrix is a dimensionless correlation matrix in
/// an explicit channel order (`Sigma = D R D^T`).
#[derive(Debug, Clone, PartialEq)]
pub struct NoiseDependence {
    order: Vec<ObservationChannelId>,
    correlation: CovarianceMatrix,
    evidence: ContentHash,
}

impl NoiseDependence {
    /// Construct a normalized correlation declaration.
    pub fn try_new(
        order: Vec<ObservationChannelId>,
        correlation: CovarianceMatrix,
        evidence: ContentHash,
    ) -> Result<Self, IdentifiabilityError> {
        if order.is_empty()
            || order.len() > MAX_IDENTIFIABILITY_ITEMS
            || order.len() != correlation.dimension()
        {
            return Err(IdentifiabilityError::Covariance {
                detail: "channel order must be bounded and match matrix dimension".to_string(),
            });
        }
        let unique = order.iter().cloned().collect::<BTreeSet<_>>();
        if unique.len() != order.len() {
            return Err(IdentifiabilityError::Covariance {
                detail: "correlation channel order contains duplicates".to_string(),
            });
        }
        if !hash_is_nonzero(evidence) {
            return Err(IdentifiabilityError::ZeroIdentity {
                field: "noise-correlation evidence",
            });
        }
        for index in 0..order.len() {
            if !same_f64(matrix_get(&correlation, index, index), 1.0) {
                return Err(IdentifiabilityError::Covariance {
                    detail: format!(
                        "correlation diagonal for channel {} is not exactly one",
                        order[index]
                    ),
                });
            }
        }
        Ok(Self {
            order,
            correlation,
            evidence,
        })
    }

    fn canonicalized(
        &self,
        canonical_order: &[ObservationChannelId],
    ) -> Result<Self, IdentifiabilityError> {
        if canonical_order.len() != self.order.len()
            || canonical_order.iter().cloned().collect::<BTreeSet<_>>()
                != self.order.iter().cloned().collect::<BTreeSet<_>>()
        {
            return Err(IdentifiabilityError::Covariance {
                detail: "correlation order is not the exact observation-channel set".to_string(),
            });
        }
        let positions = self
            .order
            .iter()
            .enumerate()
            .map(|(index, id)| (id.clone(), index))
            .collect::<BTreeMap<_, _>>();
        let mut lower = Vec::with_capacity(canonical_order.len() * (canonical_order.len() + 1) / 2);
        for (row, row_id) in canonical_order.iter().enumerate() {
            for column_id in &canonical_order[..=row] {
                lower.push(matrix_get(
                    &self.correlation,
                    positions[row_id],
                    positions[column_id],
                ));
            }
        }
        let correlation =
            CovarianceMatrix::try_new(canonical_order.len(), lower).map_err(|error| {
                IdentifiabilityError::Vv {
                    detail: error.to_string(),
                }
            })?;
        Self::try_new(canonical_order.to_vec(), correlation, self.evidence)
    }

    /// Canonical channel order.
    #[must_use]
    pub fn order(&self) -> &[ObservationChannelId] {
        &self.order
    }

    /// Dimensionless normalized correlation matrix.
    #[must_use]
    pub const fn correlation(&self) -> &CovarianceMatrix {
        &self.correlation
    }
}

/// External evidence verdict. The five identifiability axes each carry one of
/// these independently; no strength ordering is implemented.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum EvidenceStatus {
    /// Axis not run/evaluated.
    NotAssessed { reason: String },
    /// Result is unresolved despite an attempted analysis.
    Unknown { reason: String },
    /// Constructive negative result.
    Refuted { witness: ContentHash },
    /// Positive result from a named external method and receipt.
    Supported {
        method: String,
        receipt: ContentHash,
    },
}

impl EvidenceStatus {
    fn validate(&self, field: &'static str) -> Result<(), IdentifiabilityError> {
        match self {
            Self::NotAssessed { reason } | Self::Unknown { reason } => {
                validate_reason(reason, field)
            }
            Self::Refuted { witness } if hash_is_nonzero(*witness) => Ok(()),
            Self::Supported { method, receipt } if hash_is_nonzero(*receipt) => {
                validate_token(method, "identifiability method")
            }
            _ => Err(IdentifiabilityError::ZeroIdentity { field }),
        }
    }
}

/// Orthogonal identifiability evidence axes. Local, generic, global,
/// structural, and practical claims cannot be compared or promoted by enum
/// ordering.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct IdentifiabilityEvidence {
    structural: EvidenceStatus,
    local: EvidenceStatus,
    generic: EvidenceStatus,
    global: EvidenceStatus,
    practical: EvidenceStatus,
}

impl IdentifiabilityEvidence {
    /// Construct the five independent evidence axes.
    pub fn try_new(
        structural: EvidenceStatus,
        local: EvidenceStatus,
        generic: EvidenceStatus,
        global: EvidenceStatus,
        practical: EvidenceStatus,
    ) -> Result<Self, IdentifiabilityError> {
        for (field, state) in [
            ("structural evidence", &structural),
            ("local evidence", &local),
            ("generic evidence", &generic),
            ("global evidence", &global),
            ("practical evidence", &practical),
        ] {
            state.validate(field)?;
        }
        Ok(Self {
            structural,
            local,
            generic,
            global,
            practical,
        })
    }

    /// Structural-identifiability state.
    #[must_use]
    pub const fn structural(&self) -> &EvidenceStatus {
        &self.structural
    }

    /// Local-identifiability state.
    #[must_use]
    pub const fn local(&self) -> &EvidenceStatus {
        &self.local
    }

    /// Generic-identifiability state.
    #[must_use]
    pub const fn generic(&self) -> &EvidenceStatus {
        &self.generic
    }

    /// Global-identifiability state.
    #[must_use]
    pub const fn global(&self) -> &EvidenceStatus {
        &self.global
    }

    /// Practical-identifiability state.
    #[must_use]
    pub const fn practical(&self) -> &EvidenceStatus {
        &self.practical
    }
}

/// Explicit group-action/quotient declaration for gauge-equivalent physical
/// parameters. Hashes bind supplied artifacts; they are not theorem receipts
/// unless the evidence state says so.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GaugeClass {
    id: GaugeClassId,
    members: BTreeSet<ParameterRoleId>,
    continuous_dimension: u32,
    group_action: ContentHash,
    quotient_map: ContentHash,
    inverse_or_slice: ContentHash,
    stabilizer_strata: ContentHash,
    evidence: EvidenceStatus,
}

impl GaugeClass {
    /// Construct an explicit gauge quotient declaration.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        id: GaugeClassId,
        members: Vec<ParameterRoleId>,
        continuous_dimension: u32,
        group_action: ContentHash,
        quotient_map: ContentHash,
        inverse_or_slice: ContentHash,
        stabilizer_strata: ContentHash,
        evidence: EvidenceStatus,
    ) -> Result<Self, IdentifiabilityError> {
        let count = members.len();
        let members = members.into_iter().collect::<BTreeSet<_>>();
        if count < 2 || count > MAX_IDENTIFIABILITY_ITEMS || members.len() != count {
            return Err(IdentifiabilityError::InvalidGauge {
                gauge: id,
                detail: "a gauge needs at least two distinct bounded members".to_string(),
            });
        }
        if continuous_dimension == 0 || continuous_dimension as usize > members.len() {
            return Err(IdentifiabilityError::InvalidGauge {
                gauge: id,
                detail: "continuous gauge dimension must lie in 1..=member count".to_string(),
            });
        }
        for (field, hash) in [
            ("gauge group action", group_action),
            ("gauge quotient map", quotient_map),
            ("gauge inverse/slice", inverse_or_slice),
            ("gauge stabilizer strata", stabilizer_strata),
        ] {
            if !hash_is_nonzero(hash) {
                return Err(IdentifiabilityError::ZeroIdentity { field });
            }
        }
        evidence.validate("gauge evidence")?;
        Ok(Self {
            id,
            members,
            continuous_dimension,
            group_action,
            quotient_map,
            inverse_or_slice,
            stabilizer_strata,
            evidence,
        })
    }

    /// Gauge identity.
    #[must_use]
    pub const fn id(&self) -> &GaugeClassId {
        &self.id
    }

    /// Canonical parameter members.
    #[must_use]
    pub const fn members(&self) -> &BTreeSet<ParameterRoleId> {
        &self.members
    }
}

/// Per-channel discrepancy semantics. Missing and zero are deliberately
/// different variants.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DiscrepancyModel {
    /// No discrepancy model is available; this is a no-claim, not zero.
    NoModel { reason: String },
    /// Zero discrepancy is asserted by an explicit external receipt.
    Zero { evidence: ContentHash },
    /// Content-addressed discrepancy family with explicit confounding policy.
    Modeled {
        artifact: ContentHash,
        version: u32,
        confounding_constraint: ContentHash,
        evidence: EvidenceStatus,
    },
}

impl DiscrepancyModel {
    fn validate(&self) -> Result<(), IdentifiabilityError> {
        match self {
            Self::NoModel { reason } => validate_reason(reason, "discrepancy no-model reason"),
            Self::Zero { evidence } if hash_is_nonzero(*evidence) => Ok(()),
            Self::Modeled {
                artifact,
                version,
                confounding_constraint,
                evidence,
            } if hash_is_nonzero(*artifact)
                && *version > 0
                && hash_is_nonzero(*confounding_constraint) =>
            {
                evidence.validate("discrepancy evidence")
            }
            Self::Zero { .. } => Err(IdentifiabilityError::ZeroIdentity {
                field: "zero-discrepancy evidence",
            }),
            Self::Modeled { .. } => Err(IdentifiabilityError::InvalidNumeric {
                field: "discrepancy model",
                detail: "modeled discrepancy needs nonzero artifacts and a positive version"
                    .to_string(),
            }),
        }
    }
}

/// Discrepancy semantics for one observation channel.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DiscrepancySpec {
    observation: ObservationChannelId,
    model: DiscrepancyModel,
}

impl DiscrepancySpec {
    /// Construct one explicit discrepancy row.
    pub fn try_new(
        observation: ObservationChannelId,
        model: DiscrepancyModel,
    ) -> Result<Self, IdentifiabilityError> {
        model.validate()?;
        Ok(Self { observation, model })
    }

    /// Observation channel governed by this row.
    #[must_use]
    pub const fn observation(&self) -> &ObservationChannelId {
        &self.observation
    }

    /// Explicit discrepancy semantics.
    #[must_use]
    pub const fn model(&self) -> &DiscrepancyModel {
        &self.model
    }
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

fn require_header_version(
    header: &ArtifactHeader,
    component: &'static str,
    expected: u32,
) -> Result<(), IdentifiabilityError> {
    let expected_text = expected.to_string();
    match header.versions().get(component) {
        Some(actual) if actual == &expected_text => Ok(()),
        Some(actual) => Err(IdentifiabilityError::InvalidText {
            field: "study header version",
            detail: format!("component {component:?} must be {expected_text:?}, found {actual:?}"),
        }),
        None => Err(IdentifiabilityError::InvalidText {
            field: "study header version",
            detail: format!("component {component:?} is missing"),
        }),
    }
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

/// Canonical, structurally admitted law/experiment inverse-problem schema.
///
/// Fields stay private so callers cannot mutate a canonicalized map or splice
/// paths/discrepancy rows after admission. Scientific/laboratory authority is
/// still external: content hashes bind claimed artifacts but do not authenticate
/// their issuers.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct IdentifiabilityStudySpec {
    header: ArtifactHeader,
    context_of_use: ContextBinding,
    model: MaterialModelBinding,
    initial_state: InitialStateBinding,
    specimen: SpecimenBinding,
    protocol: ProtocolBinding,
    data: DataLineage,
    parameters: BTreeMap<ParameterRoleId, ParameterSpec>,
    gauges: BTreeMap<GaugeClassId, GaugeClass>,
    observations: BTreeMap<ObservationChannelId, ObservationSpec>,
    paths: Vec<ObservationPath>,
    noise_dependence: NoiseDependence,
    discrepancies: BTreeMap<ObservationChannelId, DiscrepancySpec>,
    evidence: IdentifiabilityEvidence,
}

/// Exact source artifacts required to turn canonical claims back into an
/// admitted study. Supplying these references proves byte/derivation
/// consistency only; issuer authentication remains an external policy.
#[derive(Debug, Clone, Copy)]
pub(crate) struct IdentifiabilitySourceArtifacts<'a> {
    context: &'a ContextOfUse,
    material: &'a MaterialCard,
    model: &'a ConstitutiveModelCard,
    graph: ContentHash,
    experiment: &'a ExperimentArtifact,
    split: &'a CalibrationSplit,
}

impl<'a> IdentifiabilitySourceArtifacts<'a> {
    /// Bind the exact immutable artifacts expected by canonical decoding.
    #[must_use]
    pub const fn new(
        context: &'a ContextOfUse,
        material: &'a MaterialCard,
        model: &'a ConstitutiveModelCard,
        graph: ContentHash,
        experiment: &'a ExperimentArtifact,
        split: &'a CalibrationSplit,
    ) -> Self {
        Self {
            context,
            material,
            model,
            graph,
            experiment,
            split,
        }
    }
}

impl IdentifiabilityStudySpec {
    /// Validate and canonicalize a complete study schema.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        header: ArtifactHeader,
        context_of_use: ContextBinding,
        model: MaterialModelBinding,
        initial_state: InitialStateBinding,
        specimen: SpecimenBinding,
        protocol: ProtocolBinding,
        data: DataLineage,
        parameters: Vec<ParameterSpec>,
        gauges: Vec<GaugeClass>,
        observations: Vec<ObservationSpec>,
        mut paths: Vec<ObservationPath>,
        noise_dependence: NoiseDependence,
        discrepancies: Vec<DiscrepancySpec>,
        evidence: IdentifiabilityEvidence,
    ) -> Result<Self, IdentifiabilityError> {
        validate_header_profile(&header)?;
        context_of_use.validate_structural()?;
        model.validate_structural()?;
        data.validate_structural()?;
        if data.vv_schema_version != VV_SCHEMA_VERSION {
            return Err(IdentifiabilityError::VersionMismatch {
                field: "V&V artifact schema",
                expected: VV_SCHEMA_VERSION,
                actual: data.vv_schema_version,
            });
        }
        if model.matdb_schema_version != MATDB_SCHEMA_VERSION {
            return Err(IdentifiabilityError::VersionMismatch {
                field: "material-card schema",
                expected: MATDB_SCHEMA_VERSION,
                actual: model.matdb_schema_version,
            });
        }
        initial_state.validate_against(&model)?;
        if protocol.state_schema_version != model.state_schema_version {
            return Err(IdentifiabilityError::VersionMismatch {
                field: "protocol state schema",
                expected: model.state_schema_version,
                actual: protocol.state_schema_version,
            });
        }
        for (component, expected) in [
            (
                "fs-material-identifiability",
                LEGACY_IDENTIFIABILITY_SCHEMA_VERSION,
            ),
            ("fs-evidence-vv", VV_SCHEMA_VERSION),
            ("fs-matdb", MATDB_SCHEMA_VERSION),
            ("constitutive-law", model.law_version),
            ("constitutive-state", model.state_schema_version),
            ("experiment-protocol", protocol.version),
            ("refinement-policy", protocol.refinement_version),
        ] {
            require_header_version(&header, component, expected)?;
        }
        if !header.capabilities().contains("identifiability.study") {
            return Err(IdentifiabilityError::InvalidText {
                field: "study capability",
                detail: "missing explicit identifiability.study capability".to_string(),
            });
        }

        if parameters.is_empty() || parameters.len() > MAX_IDENTIFIABILITY_ITEMS {
            return Err(IdentifiabilityError::Cardinality {
                field: "parameters",
                detail: "the study needs a bounded nonempty parameter schema".to_string(),
            });
        }
        let mut parameter_map = BTreeMap::new();
        for parameter in parameters {
            let role = parameter.role.clone();
            if parameter_map.insert(role.clone(), parameter).is_some() {
                return Err(IdentifiabilityError::Duplicate {
                    field: "parameter role",
                    id: role.to_string(),
                });
            }
        }
        for (role, roster_parameter) in &model.parameter_roster {
            let Some(parameter) = parameter_map.get(role) else {
                return Err(IdentifiabilityError::UnknownReference {
                    field: "model parameter declaration",
                    id: role.to_string(),
                });
            };
            if parameter.owner != ParameterOwner::ConstitutiveModel
                || parameter.quantity != roster_parameter.quantity
            {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "model parameter quantity/owner",
                    detail: format!(
                        "parameter {role} must be constitutive-owned with model-card dimensions"
                    ),
                });
            }
            let nominal = roster_parameter.nominal();
            if nominal < parameter.domain.lo || nominal > parameter.domain.hi {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "model parameter nominal domain",
                    detail: format!(
                        "model-card nominal value for {role} lies outside the declared physical domain"
                    ),
                });
            }
        }
        for parameter in parameter_map.values() {
            if parameter.owner == ParameterOwner::ConstitutiveModel
                && !model.parameter_roster.contains_key(&parameter.role)
            {
                return Err(IdentifiabilityError::UnknownReference {
                    field: "constitutive model parameter",
                    id: parameter.role.to_string(),
                });
            }
            if let ParameterScope::Specimen { specimen: scoped } = &parameter.scope {
                if scoped != specimen.id() {
                    return Err(IdentifiabilityError::UnknownReference {
                        field: "parameter specimen scope",
                        id: scoped.as_str().to_string(),
                    });
                }
            }
            if let ParameterClass::Nuisance { calibration } = &parameter.class {
                if calibration != data.split() {
                    return Err(IdentifiabilityError::NuisanceCalibration {
                        parameter: parameter.role.clone(),
                    });
                }
            }
        }

        if observations.is_empty() || observations.len() > MAX_IDENTIFIABILITY_ITEMS {
            return Err(IdentifiabilityError::Cardinality {
                field: "observations",
                detail: "the study needs a bounded nonempty observation schema".to_string(),
            });
        }
        let mut observation_map = BTreeMap::new();
        for observation in observations {
            let id = observation.id.clone();
            if observation.protocol_version != protocol.version {
                return Err(IdentifiabilityError::VersionMismatch {
                    field: "observation protocol",
                    expected: protocol.version,
                    actual: observation.protocol_version,
                });
            }
            if observation.refinement_version != protocol.refinement_version {
                return Err(IdentifiabilityError::VersionMismatch {
                    field: "observation refinement",
                    expected: protocol.refinement_version,
                    actual: observation.refinement_version,
                });
            }
            if observation.sensor.clock() != protocol.clock() {
                return Err(IdentifiabilityError::UnknownReference {
                    field: "observation clock",
                    id: observation.sensor.clock().as_str().to_string(),
                });
            }
            if !data.qois.contains(&observation.qoi)
                || !context_of_use.qoi_units.contains_key(&observation.qoi)
            {
                return Err(IdentifiabilityError::UnknownReference {
                    field: "experiment/context QoI",
                    id: observation.qoi.as_str().to_string(),
                });
            }
            if !observation.source_rows.is_subset(&data.calibration_ids) {
                return Err(IdentifiabilityError::Vv {
                    detail: format!(
                        "observation {} uses validation/blind/unknown rows; identifiability fitting may consume only preregistered calibration rows",
                        observation.id
                    ),
                });
            }
            if observation_map.insert(id.clone(), observation).is_some() {
                return Err(IdentifiabilityError::Duplicate {
                    field: "observation channel",
                    id: id.to_string(),
                });
            }
        }

        if paths.is_empty() || paths.len() > MAX_IDENTIFIABILITY_ITEMS {
            return Err(IdentifiabilityError::Cardinality {
                field: "observation paths",
                detail: "the study needs a bounded nonempty path schema".to_string(),
            });
        }
        paths.sort_by(|left, right| {
            (
                &left.parameter,
                &left.observation,
                &left.functional,
                left.route,
                left.graph_path,
            )
                .cmp(&(
                    &right.parameter,
                    &right.observation,
                    &right.functional,
                    right.route,
                    right.graph_path,
                ))
        });
        for pair in paths.windows(2) {
            if pair[0].parameter == pair[1].parameter
                && pair[0].observation == pair[1].observation
                && pair[0].functional == pair[1].functional
                && pair[0].route == pair[1].route
                && pair[0].graph_path == pair[1].graph_path
            {
                return Err(IdentifiabilityError::Duplicate {
                    field: "observation path route",
                    id: format!("{}->{}", pair[0].parameter, pair[0].observation),
                });
            }
        }
        let mut candidate_paths = BTreeSet::new();
        for path in &paths {
            let Some(parameter) = parameter_map.get(&path.parameter) else {
                return Err(IdentifiabilityError::UnknownReference {
                    field: "path parameter",
                    id: path.parameter.to_string(),
                });
            };
            let Some(observation) = observation_map.get(&path.observation) else {
                return Err(IdentifiabilityError::UnknownReference {
                    field: "path observation",
                    id: path.observation.to_string(),
                });
            };
            let functional_dims = match &path.functional {
                InfluenceFunctional::Mean => observation.quantity.dims(),
                InfluenceFunctional::Variance => {
                    checked_add_dims(observation.quantity.dims(), observation.quantity.dims())
                        .ok_or_else(|| IdentifiabilityError::InvalidNumeric {
                            field: "variance functional quantity",
                            detail: format!("dimension exponents overflow for {}", observation.id),
                        })?
                }
                InfluenceFunctional::Covariance { other } => {
                    if other == &path.observation {
                        return Err(IdentifiabilityError::InvalidNumeric {
                            field: "covariance functional",
                            detail: "self-covariance must use the variance functional".to_string(),
                        });
                    }
                    let Some(other_observation) = observation_map.get(other) else {
                        return Err(IdentifiabilityError::UnknownReference {
                            field: "covariance companion observation",
                            id: other.to_string(),
                        });
                    };
                    checked_add_dims(
                        observation.quantity.dims(),
                        other_observation.quantity.dims(),
                    )
                    .ok_or_else(|| IdentifiabilityError::InvalidNumeric {
                        field: "covariance functional quantity",
                        detail: format!(
                            "dimension exponents overflow for {} x {}",
                            observation.id, other_observation.id
                        ),
                    })?
                }
                InfluenceFunctional::CensoringProbability => {
                    if observation.saturation.is_none() {
                        return Err(IdentifiabilityError::InvalidNumeric {
                            field: "censoring influence",
                            detail: format!(
                                "observation {} has no saturation/censoring model",
                                observation.id
                            ),
                        });
                    }
                    Dims([0; 6])
                }
                InfluenceFunctional::MissingnessProbability => {
                    if matches!(&observation.missingness, MissingnessModel::Complete { .. }) {
                        return Err(IdentifiabilityError::InvalidNumeric {
                            field: "missingness influence",
                            detail: format!(
                                "observation {} claims complete data and cannot also have a missingness path",
                                observation.id
                            ),
                        });
                    }
                    Dims([0; 6])
                }
            };
            let expected_dims = checked_derivative_dims(functional_dims, parameter.quantity.dims())
                .ok_or_else(|| IdentifiabilityError::InvalidNumeric {
                    field: "path derivative quantity",
                    detail: format!(
                        "dimension exponents overflow for {} -> {}",
                        parameter.role, observation.id
                    ),
                })?;
            if path.derivative_quantity.semantic_type().is_some()
                || path.derivative_quantity.dims() != expected_dims
            {
                return Err(IdentifiabilityError::InvalidNumeric {
                    field: "path derivative quantity",
                    detail: format!(
                        "{} -> {} derivative must be an explicit dimension-only functional/input quotient",
                        parameter.role, observation.id
                    ),
                });
            }
            if path.status.can_support_candidate() {
                candidate_paths.insert(path.parameter.clone());
            }
        }
        for parameter in parameter_map.values() {
            if matches!(
                &parameter.class,
                ParameterClass::Target | ParameterClass::Nuisance { .. }
            ) && matches!(&parameter.observability, ParameterObservability::Candidate)
                && !candidate_paths.contains(&parameter.role)
            {
                return Err(IdentifiabilityError::DisconnectedEstimatedParameter {
                    parameter: parameter.role.clone(),
                });
            }
        }

        if gauges.len() > MAX_IDENTIFIABILITY_ITEMS {
            return Err(IdentifiabilityError::Cardinality {
                field: "gauge classes",
                detail: "too many gauge declarations".to_string(),
            });
        }
        let mut gauge_map = BTreeMap::new();
        for gauge in gauges {
            for member in &gauge.members {
                let Some(parameter) = parameter_map.get(member) else {
                    return Err(IdentifiabilityError::UnknownReference {
                        field: "gauge member",
                        id: member.to_string(),
                    });
                };
                if matches!(&parameter.class, ParameterClass::Fixed { .. }) {
                    return Err(IdentifiabilityError::InvalidGauge {
                        gauge: gauge.id.clone(),
                        detail: format!("fixed parameter {member} cannot be a free gauge member"),
                    });
                }
            }
            let id = gauge.id.clone();
            if gauge_map.insert(id.clone(), gauge).is_some() {
                return Err(IdentifiabilityError::Duplicate {
                    field: "gauge class",
                    id: id.to_string(),
                });
            }
        }

        let canonical_observations = observation_map.keys().cloned().collect::<Vec<_>>();
        let noise_dependence = noise_dependence.canonicalized(&canonical_observations)?;
        if discrepancies.len() != observation_map.len() {
            return Err(IdentifiabilityError::Cardinality {
                field: "discrepancy rows",
                detail: "every observation needs exactly one explicit discrepancy row".to_string(),
            });
        }
        let mut discrepancy_map = BTreeMap::new();
        for discrepancy in discrepancies {
            let id = discrepancy.observation.clone();
            if !observation_map.contains_key(&id) {
                return Err(IdentifiabilityError::UnknownReference {
                    field: "discrepancy observation",
                    id: id.to_string(),
                });
            }
            if discrepancy_map.insert(id.clone(), discrepancy).is_some() {
                return Err(IdentifiabilityError::Duplicate {
                    field: "discrepancy observation",
                    id: id.to_string(),
                });
            }
        }

        Ok(Self {
            header,
            context_of_use,
            model,
            initial_state,
            specimen,
            protocol,
            data,
            parameters: parameter_map,
            gauges: gauge_map,
            observations: observation_map,
            paths,
            noise_dependence,
            discrepancies: discrepancy_map,
            evidence,
        })
    }

    /// Retained schema version of this crate-private historical prototype.
    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        LEGACY_IDENTIFIABILITY_SCHEMA_VERSION
    }

    /// Five-Explicits V&V header.
    #[must_use]
    pub const fn header(&self) -> &ArtifactHeader {
        &self.header
    }

    /// Exact material/model binding.
    #[must_use]
    pub const fn model(&self) -> &MaterialModelBinding {
        &self.model
    }

    /// Canonical parameter map.
    #[must_use]
    pub const fn parameters(&self) -> &BTreeMap<ParameterRoleId, ParameterSpec> {
        &self.parameters
    }

    /// Canonical observation map.
    #[must_use]
    pub const fn observations(&self) -> &BTreeMap<ObservationChannelId, ObservationSpec> {
        &self.observations
    }

    /// Canonical observation-path rows.
    #[must_use]
    pub fn paths(&self) -> &[ObservationPath] {
        &self.paths
    }

    /// Raw/split/blind data lineage.
    #[must_use]
    pub const fn data(&self) -> &DataLineage {
        &self.data
    }

    /// Five independent identifiability evidence axes.
    #[must_use]
    pub const fn evidence(&self) -> &IdentifiabilityEvidence {
        &self.evidence
    }

    /// Exact canonical replay bytes.
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, IdentifiabilityError> {
        encode_study(self, true)
    }

    /// Decode exact canonical bytes, recompute every source-derived binding,
    /// and re-run whole-study admission. Raw bytes alone can never mint an
    /// admitted source claim.
    pub fn from_canonical_bytes(
        bytes: &[u8],
        sources: IdentifiabilitySourceArtifacts<'_>,
    ) -> Result<Self, IdentifiabilityError> {
        let study = decode_study(bytes)?;
        let resolved_context = ContextBinding::from_vv(sources.context)?;
        if study.context_of_use != resolved_context {
            return Err(IdentifiabilityError::SourceMismatch {
                field: "context of use",
            });
        }
        if study.model.graph != sources.graph {
            return Err(IdentifiabilityError::SourceMismatch {
                field: "constitutive graph",
            });
        }
        let resolved_model =
            MaterialModelBinding::from_cards(sources.material, sources.model, sources.graph)?;
        if study.model != resolved_model {
            return Err(IdentifiabilityError::SourceMismatch {
                field: "material/model",
            });
        }
        let resolved_data = DataLineage::from_vv(
            sources.experiment,
            sources.split,
            study.data.parser,
            study.data.parser_version,
            study.data.preprocessing,
            study.data.split_grouping.clone(),
        )?;
        if study.data != resolved_data {
            return Err(IdentifiabilityError::SourceMismatch {
                field: "experiment/split/data lineage",
            });
        }
        Ok(study)
    }
}

/// Complete retained identity preimage and bounded item count.
#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct StudyIdentityReceipt<I> {
    id: I,
    schema_version: u32,
    canonical_bytes: Vec<u8>,
    item_count: u64,
}

impl<I: Copy> StudyIdentityReceipt<I> {
    /// Typed identity.
    #[must_use]
    pub const fn id(&self) -> I {
        self.id
    }

    /// Canonical schema version.
    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Complete canonical preimage retained for independent adjudication.
    #[must_use]
    pub fn canonical_bytes(&self) -> &[u8] {
        &self.canonical_bytes
    }

    /// Total bounded row count across parameters, gauges, observations,
    /// paths, and discrepancies.
    #[must_use]
    pub const fn item_count(&self) -> u64 {
        self.item_count
    }
}

/// Opaque admitted study with exact and reparameterization-quotient receipts.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct AdmittedIdentifiabilityStudy {
    spec: IdentifiabilityStudySpec,
    exact: StudyIdentityReceipt<StudySpecId>,
    physical: StudyIdentityReceipt<PhysicalStudyId>,
}

impl AdmittedIdentifiabilityStudy {
    /// Mint receipts for a structurally admitted/canonicalized specification.
    pub fn admit(spec: IdentifiabilityStudySpec) -> Result<Self, IdentifiabilityError> {
        let exact_bytes = encode_study(&spec, true)?;
        let physical_bytes = encode_study(&spec, false)?;
        let item_count = u64::try_from(
            spec.parameters
                .len()
                .saturating_add(spec.gauges.len())
                .saturating_add(spec.observations.len())
                .saturating_add(spec.paths.len())
                .saturating_add(spec.discrepancies.len()),
        )
        .map_err(|_| IdentifiabilityError::Cardinality {
            field: "identity receipt item count",
            detail: "row count exceeds u64".to_string(),
        })?;
        let exact = StudyIdentityReceipt {
            id: StudySpecId(hash_domain(SPEC_DOMAIN, &exact_bytes)),
            schema_version: LEGACY_IDENTIFIABILITY_SCHEMA_VERSION,
            canonical_bytes: exact_bytes,
            item_count,
        };
        let physical = StudyIdentityReceipt {
            id: PhysicalStudyId(hash_domain(PHYSICAL_DOMAIN, &physical_bytes)),
            schema_version: LEGACY_IDENTIFIABILITY_SCHEMA_VERSION,
            canonical_bytes: physical_bytes,
            item_count,
        };
        Ok(Self {
            spec,
            exact,
            physical,
        })
    }

    /// Admitted immutable specification.
    #[must_use]
    pub const fn spec(&self) -> &IdentifiabilityStudySpec {
        &self.spec
    }

    /// Exact replay identity receipt.
    #[must_use]
    pub const fn exact_receipt(&self) -> &StudyIdentityReceipt<StudySpecId> {
        &self.exact
    }

    /// Validated reparameterization-quotient identity receipt.
    #[must_use]
    pub const fn physical_receipt(&self) -> &StudyIdentityReceipt<PhysicalStudyId> {
        &self.physical
    }
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

    fn i64(&mut self, value: i64) {
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

fn encode_artifact_kind(writer: &mut CanonicalWriter, kind: ArtifactKind) {
    writer.byte(match kind {
        ArtifactKind::ContextOfUse => 0,
        ArtifactKind::ValidationPlan => 1,
        ArtifactKind::ExperimentArtifact => 2,
        ArtifactKind::CalibrationSplit => 3,
        ArtifactKind::SolutionVerificationReceipt => 4,
        ArtifactKind::PredictionAssessment => 5,
        ArtifactKind::AssumptionsLedger => 6,
    });
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

fn encode_artifact_ref(
    writer: &mut CanonicalWriter,
    reference: &ArtifactRef,
) -> Result<(), IdentifiabilityError> {
    encode_artifact_kind(writer, reference.kind());
    encode_artifact_id(writer, reference.id())?;
    writer.hash(reference.hash());
    Ok(())
}

fn encode_context(
    writer: &mut CanonicalWriter,
    context: &ContextBinding,
) -> Result<(), IdentifiabilityError> {
    encode_artifact_ref(writer, &context.reference)?;
    writer.count(context.qoi_units.len(), "context QoIs")?;
    for (qoi, unit) in &context.qoi_units {
        encode_qoi_id(writer, qoi)?;
        writer.text(unit.as_str(), "context QoI unit")?;
    }
    Ok(())
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

fn encode_parameter_class(
    writer: &mut CanonicalWriter,
    class: &ParameterClass,
) -> Result<(), IdentifiabilityError> {
    match class {
        ParameterClass::Target => writer.byte(0),
        ParameterClass::Nuisance { calibration } => {
            writer.byte(1);
            encode_artifact_ref(writer, calibration)?;
        }
        ParameterClass::Fixed { source } => {
            writer.byte(2);
            writer.hash(*source);
        }
    }
    Ok(())
}

fn encode_parameter_scope(
    writer: &mut CanonicalWriter,
    scope: &ParameterScope,
) -> Result<(), IdentifiabilityError> {
    match scope {
        ParameterScope::Global => writer.byte(0),
        ParameterScope::MaterialLot { lot } => {
            writer.byte(1);
            encode_artifact_id(writer, lot)?;
        }
        ParameterScope::Specimen { specimen } => {
            writer.byte(2);
            encode_artifact_id(writer, specimen)?;
        }
        ParameterScope::Field { support } => {
            writer.byte(3);
            writer.hash(*support);
        }
        ParameterScope::Hierarchical { population, level } => {
            writer.byte(4);
            encode_artifact_id(writer, population)?;
            writer.u32(*level);
        }
    }
    Ok(())
}

fn encode_parameter(
    writer: &mut CanonicalWriter,
    parameter: &ParameterSpec,
    exact: bool,
) -> Result<(), IdentifiabilityError> {
    writer.text(parameter.role.as_str(), "parameter role")?;
    writer.quantity(parameter.quantity);
    encode_parameter_domain(writer, parameter.domain);
    encode_prior(writer, &parameter.prior)?;
    writer.byte(match parameter.owner {
        ParameterOwner::ConstitutiveModel => 0,
        ParameterOwner::InitialState => 1,
        ParameterOwner::Instrument => 2,
        ParameterOwner::Discrepancy => 3,
        ParameterOwner::ControlledInput => 4,
    });
    encode_parameter_scope(writer, &parameter.scope)?;
    encode_parameter_class(writer, &parameter.class)?;
    match &parameter.observability {
        ParameterObservability::Candidate => writer.byte(0),
        ParameterObservability::ExplicitlyUnidentifiable { reason, witness } => {
            writer.byte(1);
            writer.text(reason, "unidentifiable reason")?;
            writer.hash(*witness);
        }
        ParameterObservability::NotEstimated { reason } => {
            writer.byte(2);
            writer.text(reason, "not-estimated reason")?;
        }
    }
    writer.byte(u8::from(exact));
    if exact {
        encode_coordinate(writer, &parameter.coordinate)?;
    }
    Ok(())
}

fn encode_model(
    writer: &mut CanonicalWriter,
    model: &MaterialModelBinding,
) -> Result<(), IdentifiabilityError> {
    writer.hash(model.material_card);
    writer.hash(model.model_card);
    writer.hash(model.parameter_block);
    writer.hash(model.graph);
    writer.text(&model.law.0, "constitutive law id")?;
    writer.u32(model.law_version);
    writer.u32(model.state_schema_version);
    writer.byte(match model.initial_state_policy {
        InitialStatePolicy::ZeroInternalState => 0,
        InitialStatePolicy::RequiresDeclaredState => 1,
    });
    writer.u32(model.matdb_schema_version);
    writer.count(model.parameter_roster.len(), "model parameter roster")?;
    for (role, parameter) in &model.parameter_roster {
        writer.text(role.as_str(), "model parameter role")?;
        writer.quantity(parameter.quantity);
        writer.u64(parameter.nominal_bits);
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

fn encode_data(
    writer: &mut CanonicalWriter,
    data: &DataLineage,
) -> Result<(), IdentifiabilityError> {
    encode_artifact_ref(writer, &data.experiment)?;
    encode_artifact_ref(writer, &data.split)?;
    for hash in [
        data.raw_manifest,
        data.source_bytes,
        data.custody_receipt,
        data.preregistration,
        data.blind_commitment,
    ] {
        writer.hash(hash);
    }
    writer.count(data.qois.len(), "data QoIs")?;
    for qoi in &data.qois {
        encode_qoi_id(writer, qoi)?;
    }
    writer.count(data.observation_ids.len(), "raw observation ids")?;
    for id in &data.observation_ids {
        encode_observation_row_id(writer, id)?;
    }
    writer.count(data.row_bindings.len(), "raw observation bindings")?;
    for (id, binding) in &data.row_bindings {
        encode_observation_row_id(writer, id)?;
        let source = binding.source_ref();
        writer.hash(source.dataset_source_bytes_hash());
        writer.text(source.locator_domain(), "raw observation locator domain")?;
        writer.u32(source.locator_contract_version());
        writer.hash(source.locator_hash());
        writer.hash(source.extraction_receipt_hash());
        encode_qoi_id(writer, binding.qoi())?;
        encode_artifact_id(writer, binding.instrument())?;
        encode_artifact_id(writer, binding.acquisition_channel())?;
        encode_artifact_id(writer, binding.clock())?;
    }
    writer.count(data.calibration_ids.len(), "calibration observation ids")?;
    for id in &data.calibration_ids {
        encode_observation_row_id(writer, id)?;
    }
    writer.count(data.validation_ids.len(), "validation observation ids")?;
    for id in &data.validation_ids {
        encode_observation_row_id(writer, id)?;
    }
    writer.count(data.blind_sources.len(), "blind observation sources")?;
    for (id, source) in &data.blind_sources {
        encode_observation_row_id(writer, id)?;
        writer.hash(*source);
    }
    writer.hash(data.parser);
    writer.u32(data.parser_version);
    writer.hash(data.preprocessing);
    encode_artifact_id(writer, &data.split_grouping)?;
    writer.u32(data.vv_schema_version);
    Ok(())
}

fn encode_sensor(
    writer: &mut CanonicalWriter,
    sensor: &SensorBinding,
) -> Result<(), IdentifiabilityError> {
    encode_artifact_id(writer, &sensor.device)?;
    encode_artifact_id(writer, &sensor.channel)?;
    writer.hash(sensor.model);
    writer.u32(sensor.model_version);
    writer.hash(sensor.calibration_certificate);
    writer.hash(sensor.transfer_function);
    writer.hash(sensor.filter);
    writer.hash(sensor.spatial_support);
    encode_artifact_id(writer, &sensor.clock)?;
    writer.i64(sensor.delay_nanoseconds);
    writer.hash(sensor.anti_aliasing);
    Ok(())
}

fn encode_noise(
    writer: &mut CanonicalWriter,
    noise: &NoiseModel,
) -> Result<(), IdentifiabilityError> {
    match noise {
        NoiseModel::Bounded { half_width } => {
            writer.byte(0);
            writer.f64(*half_width);
        }
        NoiseModel::Gaussian { standard_deviation } => {
            writer.byte(1);
            writer.f64(*standard_deviation);
        }
        NoiseModel::StudentT {
            scale,
            degrees_of_freedom,
        } => {
            writer.byte(2);
            writer.f64(*scale);
            writer.f64(*degrees_of_freedom);
        }
        NoiseModel::Empirical {
            artifact,
            version,
            reference_scale,
        } => {
            writer.byte(3);
            writer.hash(*artifact);
            writer.u32(*version);
            writer.f64(*reference_scale);
        }
        NoiseModel::Unknown { reason } => {
            writer.byte(4);
            writer.text(reason, "noise no-claim reason")?;
        }
    }
    Ok(())
}

fn encode_missingness(
    writer: &mut CanonicalWriter,
    missingness: &MissingnessModel,
) -> Result<(), IdentifiabilityError> {
    match missingness {
        MissingnessModel::Complete { evidence } => {
            writer.byte(0);
            writer.hash(*evidence);
        }
        MissingnessModel::Modeled { artifact, version } => {
            writer.byte(1);
            writer.hash(*artifact);
            writer.u32(*version);
        }
        MissingnessModel::Unknown { reason } => {
            writer.byte(2);
            writer.text(reason, "missingness no-claim reason")?;
        }
    }
    Ok(())
}

fn encode_observation(
    writer: &mut CanonicalWriter,
    observation: &ObservationSpec,
) -> Result<(), IdentifiabilityError> {
    writer.text(observation.id.as_str(), "observation channel id")?;
    encode_qoi_id(writer, &observation.qoi)?;
    writer.quantity(observation.quantity);
    encode_frame(writer, &observation.frame)?;
    writer.text(&observation.graph_node, "observation graph node")?;
    writer.text(&observation.graph_port, "observation graph port")?;
    writer.hash(observation.operator);
    writer.u32(observation.operator_version);
    writer.hash(observation.aggregation);
    encode_sensor(writer, &observation.sensor)?;
    encode_noise(writer, &observation.noise)?;
    encode_missingness(writer, &observation.missingness)?;
    match observation.saturation {
        Some(domain) => {
            writer.byte(1);
            encode_parameter_domain(writer, domain);
        }
        None => writer.byte(0),
    }
    writer.u32(observation.protocol_version);
    writer.u32(observation.refinement_version);
    writer.count(observation.source_rows.len(), "observation source rows")?;
    for row in &observation.source_rows {
        encode_observation_row_id(writer, row)?;
    }
    Ok(())
}

fn encode_influence_status(
    writer: &mut CanonicalWriter,
    status: &InfluenceStatus,
) -> Result<(), IdentifiabilityError> {
    match status {
        InfluenceStatus::DeclaredConnectivity => writer.byte(0),
        InfluenceStatus::SymbolicallyNonzero { receipt } => {
            writer.byte(1);
            writer.hash(*receipt);
        }
        InfluenceStatus::NumericallyWitnessed { receipt } => {
            writer.byte(2);
            writer.hash(*receipt);
        }
        InfluenceStatus::ProvenZero { witness } => {
            writer.byte(3);
            writer.hash(*witness);
        }
        InfluenceStatus::Unknown { reason } => {
            writer.byte(4);
            writer.text(reason, "observation path unknown reason")?;
        }
    }
    Ok(())
}

fn encode_path(
    writer: &mut CanonicalWriter,
    path: &ObservationPath,
) -> Result<(), IdentifiabilityError> {
    writer.text(path.parameter.as_str(), "path parameter")?;
    writer.text(path.observation.as_str(), "path observation")?;
    match &path.functional {
        InfluenceFunctional::Mean => writer.byte(0),
        InfluenceFunctional::Variance => writer.byte(1),
        InfluenceFunctional::Covariance { other } => {
            writer.byte(2);
            writer.text(other.as_str(), "covariance companion observation")?;
        }
        InfluenceFunctional::CensoringProbability => writer.byte(3),
        InfluenceFunctional::MissingnessProbability => writer.byte(4),
    }
    writer.byte(match path.route {
        InfluenceRoute::Direct => 0,
        InfluenceRoute::StateMediated => 1,
    });
    writer.hash(path.graph_path);
    writer.quantity(path.derivative_quantity);
    encode_influence_status(writer, &path.status)
}

fn encode_evidence_status(
    writer: &mut CanonicalWriter,
    status: &EvidenceStatus,
) -> Result<(), IdentifiabilityError> {
    match status {
        EvidenceStatus::NotAssessed { reason } => {
            writer.byte(0);
            writer.text(reason, "evidence not-assessed reason")?;
        }
        EvidenceStatus::Unknown { reason } => {
            writer.byte(1);
            writer.text(reason, "evidence unknown reason")?;
        }
        EvidenceStatus::Refuted { witness } => {
            writer.byte(2);
            writer.hash(*witness);
        }
        EvidenceStatus::Supported { method, receipt } => {
            writer.byte(3);
            writer.text(method, "evidence method")?;
            writer.hash(*receipt);
        }
    }
    Ok(())
}

fn encode_evidence(
    writer: &mut CanonicalWriter,
    evidence: &IdentifiabilityEvidence,
) -> Result<(), IdentifiabilityError> {
    for status in [
        &evidence.structural,
        &evidence.local,
        &evidence.generic,
        &evidence.global,
        &evidence.practical,
    ] {
        encode_evidence_status(writer, status)?;
    }
    Ok(())
}

fn encode_gauge(
    writer: &mut CanonicalWriter,
    gauge: &GaugeClass,
) -> Result<(), IdentifiabilityError> {
    writer.text(gauge.id.as_str(), "gauge id")?;
    writer.count(gauge.members.len(), "gauge members")?;
    for member in &gauge.members {
        writer.text(member.as_str(), "gauge member")?;
    }
    writer.u32(gauge.continuous_dimension);
    writer.hash(gauge.group_action);
    writer.hash(gauge.quotient_map);
    writer.hash(gauge.inverse_or_slice);
    writer.hash(gauge.stabilizer_strata);
    encode_evidence_status(writer, &gauge.evidence)
}

fn encode_noise_dependence(
    writer: &mut CanonicalWriter,
    dependence: &NoiseDependence,
) -> Result<(), IdentifiabilityError> {
    writer.count(dependence.order.len(), "correlation order")?;
    for id in &dependence.order {
        writer.text(id.as_str(), "correlation channel")?;
    }
    writer.count(
        dependence.correlation.lower_triangle().len(),
        "correlation lower triangle",
    )?;
    for value in dependence.correlation.lower_triangle() {
        writer.f64(*value);
    }
    writer.hash(dependence.evidence);
    Ok(())
}

fn encode_discrepancy(
    writer: &mut CanonicalWriter,
    discrepancy: &DiscrepancySpec,
) -> Result<(), IdentifiabilityError> {
    writer.text(discrepancy.observation.as_str(), "discrepancy observation")?;
    match &discrepancy.model {
        DiscrepancyModel::NoModel { reason } => {
            writer.byte(0);
            writer.text(reason, "discrepancy no-model reason")?;
        }
        DiscrepancyModel::Zero { evidence } => {
            writer.byte(1);
            writer.hash(*evidence);
        }
        DiscrepancyModel::Modeled {
            artifact,
            version,
            confounding_constraint,
            evidence,
        } => {
            writer.byte(2);
            writer.hash(*artifact);
            writer.u32(*version);
            writer.hash(*confounding_constraint);
            encode_evidence_status(writer, evidence)?;
        }
    }
    Ok(())
}

fn encode_study(
    study: &IdentifiabilityStudySpec,
    exact: bool,
) -> Result<Vec<u8>, IdentifiabilityError> {
    let mut writer = CanonicalWriter::new();
    writer.raw(CANONICAL_MAGIC);
    writer.u32(LEGACY_IDENTIFIABILITY_SCHEMA_VERSION);
    writer.byte(u8::from(exact));
    if exact {
        encode_header(&mut writer, &study.header, true)?;
    }
    encode_context(&mut writer, &study.context_of_use)?;
    encode_model(&mut writer, &study.model)?;
    encode_initial_state(&mut writer, study.initial_state);
    encode_specimen(&mut writer, &study.specimen)?;
    encode_protocol(&mut writer, &study.protocol)?;
    encode_data(&mut writer, &study.data)?;
    writer.count(study.parameters.len(), "parameters")?;
    for parameter in study.parameters.values() {
        encode_parameter(&mut writer, parameter, exact)?;
    }
    writer.count(study.gauges.len(), "gauge classes")?;
    for gauge in study.gauges.values() {
        encode_gauge(&mut writer, gauge)?;
    }
    writer.count(study.observations.len(), "observations")?;
    for observation in study.observations.values() {
        encode_observation(&mut writer, observation)?;
    }
    writer.count(study.paths.len(), "observation paths")?;
    for path in &study.paths {
        encode_path(&mut writer, path)?;
    }
    encode_noise_dependence(&mut writer, &study.noise_dependence)?;
    writer.count(study.discrepancies.len(), "discrepancy rows")?;
    for discrepancy in study.discrepancies.values() {
        encode_discrepancy(&mut writer, discrepancy)?;
    }
    encode_evidence(&mut writer, &study.evidence)?;
    writer.finish()
}

/// Fail closed on stale/future retained study schemas.
pub(crate) fn check_legacy_identifiability_schema_version(
    declared: u32,
) -> Result<(), IdentifiabilityError> {
    if declared == LEGACY_IDENTIFIABILITY_SCHEMA_VERSION {
        Ok(())
    } else {
        Err(IdentifiabilityError::UnsupportedSchemaVersion {
            declared,
            supported: LEGACY_IDENTIFIABILITY_SCHEMA_VERSION,
        })
    }
}

/// Historical single-case prototype description. It is deliberately not an
/// active identity-governance declaration and cannot mint current authority.
const LEGACY_IDENTIFIABILITY_SPEC_NONAUTHORITATIVE_DESCRIPTION: &[&str] = &[
    "legacy-nonauthoritative-identifiability-description-v0",
    "id=fs-material:identifiability-spec",
    "version_const=LEGACY_IDENTIFIABILITY_SCHEMA_VERSION",
    "version=2",
    "domain=org.frankensim.fs-material.identifiability-spec.v2",
    "domain_const=SPEC_DOMAIN",
    "encoder=encode_study(exact=true)",
    "digest=fs-blake3-domain-separated",
    "encoding=bounded-typed-binary",
    "sources=IdentifiabilityStudySpec",
    "semantic_fields=all-fields-including-header-id-and-parameter-coordinate",
    "excluded_fields=none",
    "consumers=AdmittedIdentifiabilityStudy::admit,StudySpecId",
    "version_guard=check_legacy_identifiability_schema_version",
    "coupling_surface=fs-material:identifiability-spec",
];

/// Historical coordinate-quotient prototype description. Current physical
/// authority is represented by `authoritative::ProblemId` after source
/// resolution.
const LEGACY_IDENTIFIABILITY_PHYSICAL_NONAUTHORITATIVE_DESCRIPTION: &[&str] = &[
    "legacy-nonauthoritative-identifiability-description-v0",
    "id=fs-material:identifiability-physical",
    "version_const=LEGACY_IDENTIFIABILITY_SCHEMA_VERSION",
    "version=2",
    "domain=org.frankensim.fs-material.identifiability-physical.v2",
    "domain_const=PHYSICAL_DOMAIN",
    "encoder=encode_study(exact=false)",
    "digest=fs-blake3-domain-separated",
    "encoding=bounded-typed-binary",
    "sources=IdentifiabilityStudySpec",
    "semantic_fields=all-physical-study-fields-except-header-artifact-id-and-validated-coordinate-chart",
    "excluded_fields=ArtifactHeader.id:wire-only,ParameterCoordinate:validated-bijective-chart-only",
    "consumers=AdmittedIdentifiabilityStudy::admit,PhysicalStudyId",
    "version_guard=check_legacy_identifiability_schema_version",
    "coupling_surface=fs-material:identifiability-physical",
];

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

    fn i64(&mut self, field: &'static str) -> Result<i64, IdentifiabilityError> {
        let bytes: [u8; 8] =
            self.take(8, field)?
                .try_into()
                .map_err(|_| IdentifiabilityError::Canonical {
                    at: self.at,
                    detail: format!("invalid {field} width"),
                })?;
        Ok(i64::from_le_bytes(bytes))
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

fn decode_artifact_kind(
    reader: &mut CanonicalReader<'_>,
) -> Result<ArtifactKind, IdentifiabilityError> {
    match reader.byte("artifact kind")? {
        0 => Ok(ArtifactKind::ContextOfUse),
        1 => Ok(ArtifactKind::ValidationPlan),
        2 => Ok(ArtifactKind::ExperimentArtifact),
        3 => Ok(ArtifactKind::CalibrationSplit),
        4 => Ok(ArtifactKind::SolutionVerificationReceipt),
        5 => Ok(ArtifactKind::PredictionAssessment),
        6 => Ok(ArtifactKind::AssumptionsLedger),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown artifact kind tag {tag}"),
        }),
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

fn decode_artifact_ref(
    reader: &mut CanonicalReader<'_>,
) -> Result<ArtifactRef, IdentifiabilityError> {
    Ok(ArtifactRef::new(
        decode_artifact_kind(reader)?,
        decode_artifact_id(reader)?,
        reader.hash("artifact reference")?,
    ))
}

fn decode_context(
    reader: &mut CanonicalReader<'_>,
) -> Result<ContextBinding, IdentifiabilityError> {
    let reference = decode_artifact_ref(reader)?;
    let qoi_count = reader.count("context QoIs")?;
    let mut qoi_units = BTreeMap::new();
    for _ in 0..qoi_count {
        let qoi = decode_qoi_id(reader)?;
        let unit = UnitId::try_new(reader.token("context QoI unit")?).map_err(|error| {
            IdentifiabilityError::Vv {
                detail: error.to_string(),
            }
        })?;
        if qoi_units.insert(qoi.clone(), unit).is_some() {
            return Err(IdentifiabilityError::Duplicate {
                field: "context QoI",
                id: qoi.as_str().to_string(),
            });
        }
    }
    Ok(ContextBinding {
        reference,
        qoi_units,
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

fn decode_parameter_scope(
    reader: &mut CanonicalReader<'_>,
) -> Result<ParameterScope, IdentifiabilityError> {
    match reader.byte("parameter scope tag")? {
        0 => Ok(ParameterScope::Global),
        1 => Ok(ParameterScope::MaterialLot {
            lot: decode_artifact_id(reader)?,
        }),
        2 => Ok(ParameterScope::Specimen {
            specimen: decode_artifact_id(reader)?,
        }),
        3 => Ok(ParameterScope::Field {
            support: reader.hash("parameter field support")?,
        }),
        4 => Ok(ParameterScope::Hierarchical {
            population: decode_artifact_id(reader)?,
            level: reader.u32("hierarchical level")?,
        }),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown parameter scope tag {tag}"),
        }),
    }
}

fn decode_parameter_class(
    reader: &mut CanonicalReader<'_>,
) -> Result<ParameterClass, IdentifiabilityError> {
    match reader.byte("parameter class tag")? {
        0 => Ok(ParameterClass::Target),
        1 => Ok(ParameterClass::Nuisance {
            calibration: decode_artifact_ref(reader)?,
        }),
        2 => Ok(ParameterClass::Fixed {
            source: reader.hash("fixed parameter source")?,
        }),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown parameter class tag {tag}"),
        }),
    }
}

fn decode_parameter(
    reader: &mut CanonicalReader<'_>,
) -> Result<ParameterSpec, IdentifiabilityError> {
    let role = ParameterRoleId::try_new(reader.token("parameter role")?)?;
    let quantity = reader.quantity("parameter")?;
    let domain = decode_parameter_domain(reader)?;
    let prior = decode_prior(reader)?;
    let owner = match reader.byte("parameter owner tag")? {
        0 => ParameterOwner::ConstitutiveModel,
        1 => ParameterOwner::InitialState,
        2 => ParameterOwner::Instrument,
        3 => ParameterOwner::Discrepancy,
        4 => ParameterOwner::ControlledInput,
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown parameter owner tag {tag}"),
            });
        }
    };
    let scope = decode_parameter_scope(reader)?;
    let class = decode_parameter_class(reader)?;
    let observability = match reader.byte("parameter observability tag")? {
        0 => ParameterObservability::Candidate,
        1 => ParameterObservability::ExplicitlyUnidentifiable {
            reason: reader.reason("unidentifiable reason")?,
            witness: reader.hash("unidentifiable witness")?,
        },
        2 => ParameterObservability::NotEstimated {
            reason: reader.reason("not-estimated reason")?,
        },
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown parameter observability tag {tag}"),
            });
        }
    };
    reader.expect_byte(1, "exact coordinate marker")?;
    let coordinate = decode_coordinate(reader)?;
    ParameterSpec::try_new(
        role,
        quantity,
        domain,
        prior,
        coordinate,
        owner,
        scope,
        class,
        observability,
    )
}

fn decode_model(
    reader: &mut CanonicalReader<'_>,
) -> Result<MaterialModelBinding, IdentifiabilityError> {
    let material_card = reader.hash("material card")?;
    let model_card = reader.hash("constitutive model card")?;
    let parameter_block = reader.hash("canonical parameter block")?;
    let graph = reader.hash("constitutive graph binding")?;
    let law = LawId(reader.token("constitutive law id")?);
    let law_version = reader.u32("law version")?;
    let state_schema_version = reader.u32("state schema version")?;
    let initial_state_policy = match reader.byte("initial-state policy tag")? {
        0 => InitialStatePolicy::ZeroInternalState,
        1 => InitialStatePolicy::RequiresDeclaredState,
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown initial-state policy tag {tag}"),
            });
        }
    };
    let matdb_schema_version = reader.u32("material-card schema version")?;
    let roster_count = reader.count("model parameter roster")?;
    let mut parameter_roster = BTreeMap::new();
    for _ in 0..roster_count {
        let role = ParameterRoleId::try_new(reader.token("model parameter role")?)?;
        let quantity = reader.quantity("model parameter")?;
        let nominal_bits = reader.u64("model parameter nominal bits")?;
        let nominal = f64::from_bits(nominal_bits);
        if !nominal.is_finite() || canonical_f64(nominal).to_bits() != nominal_bits {
            return Err(IdentifiabilityError::InvalidNumeric {
                field: "model parameter nominal",
                detail: format!("parameter {role} has nonfinite or noncanonical nominal bits"),
            });
        }
        if parameter_roster
            .insert(
                role.clone(),
                ModelParameterBinding {
                    quantity,
                    nominal_bits,
                },
            )
            .is_some()
        {
            return Err(IdentifiabilityError::Duplicate {
                field: "model parameter roster",
                id: role.to_string(),
            });
        }
    }
    Ok(MaterialModelBinding {
        material_card,
        model_card,
        parameter_block,
        graph,
        law,
        law_version,
        state_schema_version,
        initial_state_policy,
        matdb_schema_version,
        parameter_roster,
    })
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

fn decode_data(reader: &mut CanonicalReader<'_>) -> Result<DataLineage, IdentifiabilityError> {
    let experiment = decode_artifact_ref(reader)?;
    let split = decode_artifact_ref(reader)?;
    let raw_manifest = reader.hash("raw observation manifest")?;
    let source_bytes = reader.hash("raw source bytes")?;
    let custody_receipt = reader.hash("custody receipt")?;
    let preregistration = reader.hash("split preregistration")?;
    let blind_commitment = reader.hash("blind commitment")?;
    let qoi_count = reader.count("data QoIs")?;
    let mut qois = BTreeSet::new();
    for _ in 0..qoi_count {
        let qoi = decode_qoi_id(reader)?;
        if !qois.insert(qoi.clone()) {
            return Err(IdentifiabilityError::Duplicate {
                field: "data QoI",
                id: qoi.as_str().to_string(),
            });
        }
    }
    let observation_count = reader.count("raw observation ids")?;
    let mut observation_ids = BTreeSet::new();
    for _ in 0..observation_count {
        let id = decode_observation_row_id(reader)?;
        if !observation_ids.insert(id.clone()) {
            return Err(IdentifiabilityError::Duplicate {
                field: "raw observation row",
                id: id.as_str().to_string(),
            });
        }
    }
    let row_binding_count = reader.count("raw observation bindings")?;
    let mut row_bindings = BTreeMap::new();
    for _ in 0..row_binding_count {
        let id = decode_observation_row_id(reader)?;
        let source = ObservationSourceRef::try_new(
            reader.hash("raw observation dataset source bytes")?,
            reader.token("raw observation locator domain")?,
            reader.u32("raw observation locator contract version")?,
            reader.hash("raw observation locator")?,
            reader.hash("raw observation extraction receipt")?,
        )
        .map_err(|error| IdentifiabilityError::Vv {
            detail: error.to_string(),
        })?;
        let binding = ObservationManifestRow::try_new(
            source,
            decode_qoi_id(reader)?,
            decode_artifact_id(reader)?,
            decode_artifact_id(reader)?,
            decode_artifact_id(reader)?,
        )
        .map_err(|error| IdentifiabilityError::Vv {
            detail: error.to_string(),
        })?;
        if row_bindings.insert(id.clone(), binding).is_some() {
            return Err(IdentifiabilityError::Duplicate {
                field: "raw observation binding",
                id: id.as_str().to_string(),
            });
        }
    }
    let calibration_id_count = reader.count("calibration observation ids")?;
    let mut calibration_ids = BTreeSet::new();
    for _ in 0..calibration_id_count {
        let id = decode_observation_row_id(reader)?;
        if !calibration_ids.insert(id.clone()) {
            return Err(IdentifiabilityError::Duplicate {
                field: "calibration observation row",
                id: id.as_str().to_string(),
            });
        }
    }
    let validation_id_count = reader.count("validation observation ids")?;
    let mut validation_ids = BTreeSet::new();
    for _ in 0..validation_id_count {
        let id = decode_observation_row_id(reader)?;
        if !validation_ids.insert(id.clone()) {
            return Err(IdentifiabilityError::Duplicate {
                field: "validation observation row",
                id: id.as_str().to_string(),
            });
        }
    }
    let blind_source_count = reader.count("blind observation sources")?;
    let mut blind_sources = BTreeMap::new();
    for _ in 0..blind_source_count {
        let id = decode_observation_row_id(reader)?;
        let source = reader.hash("blind observation source")?;
        if blind_sources.insert(id.clone(), source).is_some() {
            return Err(IdentifiabilityError::Duplicate {
                field: "blind observation row",
                id: id.as_str().to_string(),
            });
        }
    }
    Ok(DataLineage {
        experiment,
        split,
        raw_manifest,
        source_bytes,
        custody_receipt,
        preregistration,
        blind_commitment,
        qois,
        observation_ids,
        row_bindings,
        calibration_ids,
        validation_ids,
        blind_sources,
        parser: reader.hash("observation parser")?,
        parser_version: reader.u32("observation parser version")?,
        preprocessing: reader.hash("preprocessing pipeline")?,
        split_grouping: decode_artifact_id(reader)?,
        vv_schema_version: reader.u32("V&V schema version")?,
    })
}

fn decode_sensor(reader: &mut CanonicalReader<'_>) -> Result<SensorBinding, IdentifiabilityError> {
    SensorBinding::try_new(
        decode_artifact_id(reader)?,
        decode_artifact_id(reader)?,
        reader.hash("sensor model")?,
        reader.u32("sensor model version")?,
        reader.hash("sensor calibration certificate")?,
        reader.hash("sensor transfer function")?,
        reader.hash("sensor filter")?,
        reader.hash("sensor spatial support")?,
        decode_artifact_id(reader)?,
        reader.i64("sensor delay nanoseconds")?,
        reader.hash("sensor anti-aliasing policy")?,
    )
}

fn decode_noise(reader: &mut CanonicalReader<'_>) -> Result<NoiseModel, IdentifiabilityError> {
    match reader.byte("noise model tag")? {
        0 => Ok(NoiseModel::Bounded {
            half_width: reader.f64("bounded-noise half width")?,
        }),
        1 => Ok(NoiseModel::Gaussian {
            standard_deviation: reader.f64("Gaussian noise standard deviation")?,
        }),
        2 => Ok(NoiseModel::StudentT {
            scale: reader.f64("Student-t noise scale")?,
            degrees_of_freedom: reader.f64("Student-t degrees of freedom")?,
        }),
        3 => Ok(NoiseModel::Empirical {
            artifact: reader.hash("empirical noise artifact")?,
            version: reader.u32("empirical noise version")?,
            reference_scale: reader.f64("empirical noise reference scale")?,
        }),
        4 => Ok(NoiseModel::Unknown {
            reason: reader.reason("noise no-claim reason")?,
        }),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown noise model tag {tag}"),
        }),
    }
}

fn decode_missingness(
    reader: &mut CanonicalReader<'_>,
) -> Result<MissingnessModel, IdentifiabilityError> {
    match reader.byte("missingness model tag")? {
        0 => Ok(MissingnessModel::Complete {
            evidence: reader.hash("data-completeness evidence")?,
        }),
        1 => Ok(MissingnessModel::Modeled {
            artifact: reader.hash("missingness model artifact")?,
            version: reader.u32("missingness model version")?,
        }),
        2 => Ok(MissingnessModel::Unknown {
            reason: reader.reason("missingness no-claim reason")?,
        }),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown missingness model tag {tag}"),
        }),
    }
}

fn decode_observation(
    reader: &mut CanonicalReader<'_>,
) -> Result<ObservationSpec, IdentifiabilityError> {
    let id = ObservationChannelId::try_new(reader.token("observation channel id")?)?;
    let qoi = decode_qoi_id(reader)?;
    let quantity = reader.quantity("observation")?;
    let frame = decode_frame(reader)?;
    let graph_node = reader.token("observation graph node")?;
    let graph_port = reader.token("observation graph port")?;
    let operator = reader.hash("observation operator")?;
    let operator_version = reader.u32("observation operator version")?;
    let aggregation = reader.hash("observation aggregation")?;
    let sensor = decode_sensor(reader)?;
    let noise = decode_noise(reader)?;
    let missingness = decode_missingness(reader)?;
    let saturation = match reader.byte("observation saturation tag")? {
        0 => None,
        1 => Some(decode_parameter_domain(reader)?),
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown observation saturation tag {tag}"),
            });
        }
    };
    let protocol_version = reader.u32("observation protocol version")?;
    let refinement_version = reader.u32("observation refinement version")?;
    let source_count = reader.count("observation source rows")?;
    let mut source_rows = Vec::with_capacity(source_count);
    for _ in 0..source_count {
        source_rows.push(decode_observation_row_id(reader)?);
    }
    ObservationSpec::try_new(
        id,
        qoi,
        quantity,
        frame,
        graph_node,
        graph_port,
        operator,
        operator_version,
        aggregation,
        sensor,
        noise,
        missingness,
        saturation,
        protocol_version,
        refinement_version,
        source_rows,
    )
}

fn decode_influence_status(
    reader: &mut CanonicalReader<'_>,
) -> Result<InfluenceStatus, IdentifiabilityError> {
    match reader.byte("observation path status tag")? {
        0 => Ok(InfluenceStatus::DeclaredConnectivity),
        1 => Ok(InfluenceStatus::SymbolicallyNonzero {
            receipt: reader.hash("symbolic path receipt")?,
        }),
        2 => Ok(InfluenceStatus::NumericallyWitnessed {
            receipt: reader.hash("numerical path receipt")?,
        }),
        3 => Ok(InfluenceStatus::ProvenZero {
            witness: reader.hash("zero-path witness")?,
        }),
        4 => Ok(InfluenceStatus::Unknown {
            reason: reader.reason("observation path unknown reason")?,
        }),
        tag => Err(IdentifiabilityError::Canonical {
            at: reader.at.saturating_sub(1),
            detail: format!("unknown observation path status tag {tag}"),
        }),
    }
}

fn decode_path(reader: &mut CanonicalReader<'_>) -> Result<ObservationPath, IdentifiabilityError> {
    let parameter = ParameterRoleId::try_new(reader.token("path parameter")?)?;
    let observation = ObservationChannelId::try_new(reader.token("path observation")?)?;
    let functional = match reader.byte("influence functional tag")? {
        0 => InfluenceFunctional::Mean,
        1 => InfluenceFunctional::Variance,
        2 => InfluenceFunctional::Covariance {
            other: ObservationChannelId::try_new(
                reader.token("covariance companion observation")?,
            )?,
        },
        3 => InfluenceFunctional::CensoringProbability,
        4 => InfluenceFunctional::MissingnessProbability,
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown influence functional tag {tag}"),
            });
        }
    };
    let route = match reader.byte("influence route tag")? {
        0 => InfluenceRoute::Direct,
        1 => InfluenceRoute::StateMediated,
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown influence route tag {tag}"),
            });
        }
    };
    ObservationPath::try_new(
        parameter,
        observation,
        functional,
        route,
        reader.hash("parameter-observation graph path")?,
        reader.quantity("path derivative")?,
        decode_influence_status(reader)?,
    )
}

fn decode_evidence_status(
    reader: &mut CanonicalReader<'_>,
    field: &'static str,
) -> Result<EvidenceStatus, IdentifiabilityError> {
    let status = match reader.byte("identifiability evidence tag")? {
        0 => EvidenceStatus::NotAssessed {
            reason: reader.reason("evidence not-assessed reason")?,
        },
        1 => EvidenceStatus::Unknown {
            reason: reader.reason("evidence unknown reason")?,
        },
        2 => EvidenceStatus::Refuted {
            witness: reader.hash("refuted-evidence witness")?,
        },
        3 => EvidenceStatus::Supported {
            method: reader.token("evidence method")?,
            receipt: reader.hash("supported-evidence receipt")?,
        },
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown identifiability evidence tag {tag}"),
            });
        }
    };
    status.validate(field)?;
    Ok(status)
}

fn decode_evidence(
    reader: &mut CanonicalReader<'_>,
) -> Result<IdentifiabilityEvidence, IdentifiabilityError> {
    IdentifiabilityEvidence::try_new(
        decode_evidence_status(reader, "structural evidence")?,
        decode_evidence_status(reader, "local evidence")?,
        decode_evidence_status(reader, "generic evidence")?,
        decode_evidence_status(reader, "global evidence")?,
        decode_evidence_status(reader, "practical evidence")?,
    )
}

fn decode_gauge(reader: &mut CanonicalReader<'_>) -> Result<GaugeClass, IdentifiabilityError> {
    let id = GaugeClassId::try_new(reader.token("gauge id")?)?;
    let member_count = reader.count("gauge members")?;
    let mut members = Vec::with_capacity(member_count);
    for _ in 0..member_count {
        members.push(ParameterRoleId::try_new(reader.token("gauge member")?)?);
    }
    GaugeClass::try_new(
        id,
        members,
        reader.u32("continuous gauge dimension")?,
        reader.hash("gauge group action")?,
        reader.hash("gauge quotient map")?,
        reader.hash("gauge inverse or slice")?,
        reader.hash("gauge stabilizer strata")?,
        decode_evidence_status(reader, "gauge evidence")?,
    )
}

fn decode_noise_dependence(
    reader: &mut CanonicalReader<'_>,
) -> Result<NoiseDependence, IdentifiabilityError> {
    let order_count = reader.count("correlation order")?;
    let mut order = Vec::with_capacity(order_count);
    for _ in 0..order_count {
        order.push(ObservationChannelId::try_new(
            reader.token("correlation channel")?,
        )?);
    }
    let expected_entries = order_count
        .checked_mul(order_count.saturating_add(1))
        .and_then(|value| value.checked_div(2))
        .ok_or_else(|| IdentifiabilityError::Covariance {
            detail: "correlation lower-triangle size overflows".to_string(),
        })?;
    let entry_count = reader.length(
        MAX_IDENTIFIABILITY_CANONICAL_BYTES / core::mem::size_of::<f64>(),
        "correlation lower triangle",
    )?;
    if entry_count != expected_entries {
        return Err(IdentifiabilityError::Covariance {
            detail: format!(
                "correlation dimension {order_count} needs {expected_entries} lower-triangle entries, found {entry_count}"
            ),
        });
    }
    let mut lower_triangle = Vec::with_capacity(entry_count);
    for _ in 0..entry_count {
        lower_triangle.push(reader.f64("correlation entry")?);
    }
    let correlation = CovarianceMatrix::try_new(order_count, lower_triangle).map_err(|error| {
        IdentifiabilityError::Vv {
            detail: error.to_string(),
        }
    })?;
    NoiseDependence::try_new(
        order,
        correlation,
        reader.hash("noise-correlation evidence")?,
    )
}

fn decode_discrepancy(
    reader: &mut CanonicalReader<'_>,
) -> Result<DiscrepancySpec, IdentifiabilityError> {
    let observation = ObservationChannelId::try_new(reader.token("discrepancy observation")?)?;
    let model = match reader.byte("discrepancy model tag")? {
        0 => DiscrepancyModel::NoModel {
            reason: reader.reason("discrepancy no-model reason")?,
        },
        1 => DiscrepancyModel::Zero {
            evidence: reader.hash("zero-discrepancy evidence")?,
        },
        2 => DiscrepancyModel::Modeled {
            artifact: reader.hash("discrepancy model artifact")?,
            version: reader.u32("discrepancy model version")?,
            confounding_constraint: reader.hash("discrepancy confounding constraint")?,
            evidence: decode_evidence_status(reader, "discrepancy evidence")?,
        },
        tag => {
            return Err(IdentifiabilityError::Canonical {
                at: reader.at.saturating_sub(1),
                detail: format!("unknown discrepancy model tag {tag}"),
            });
        }
    };
    DiscrepancySpec::try_new(observation, model)
}

fn decode_study(bytes: &[u8]) -> Result<IdentifiabilityStudySpec, IdentifiabilityError> {
    let mut reader = CanonicalReader::new(bytes)?;
    let magic_at = reader.at;
    if reader.take(CANONICAL_MAGIC.len(), "study magic")? != CANONICAL_MAGIC {
        return Err(IdentifiabilityError::Canonical {
            at: magic_at,
            detail: "wrong identifiability-study magic".to_string(),
        });
    }
    check_legacy_identifiability_schema_version(reader.u32("study schema version")?)?;
    reader.expect_byte(1, "exact study marker")?;
    let header = decode_header(&mut reader)?;
    let context_of_use = decode_context(&mut reader)?;
    let model = decode_model(&mut reader)?;
    let initial_state = decode_initial_state(&mut reader)?;
    let specimen = decode_specimen(&mut reader)?;
    let protocol = decode_protocol(&mut reader)?;
    let data = decode_data(&mut reader)?;

    let parameter_count = reader.count("parameters")?;
    let mut parameters = Vec::with_capacity(parameter_count);
    for _ in 0..parameter_count {
        parameters.push(decode_parameter(&mut reader)?);
    }
    let gauge_count = reader.count("gauge classes")?;
    let mut gauges = Vec::with_capacity(gauge_count);
    for _ in 0..gauge_count {
        gauges.push(decode_gauge(&mut reader)?);
    }
    let observation_count = reader.count("observations")?;
    let mut observations = Vec::with_capacity(observation_count);
    for _ in 0..observation_count {
        observations.push(decode_observation(&mut reader)?);
    }
    let path_count = reader.count("observation paths")?;
    let mut paths = Vec::with_capacity(path_count);
    for _ in 0..path_count {
        paths.push(decode_path(&mut reader)?);
    }
    let noise_dependence = decode_noise_dependence(&mut reader)?;
    let discrepancy_count = reader.count("discrepancy rows")?;
    let mut discrepancies = Vec::with_capacity(discrepancy_count);
    for _ in 0..discrepancy_count {
        discrepancies.push(decode_discrepancy(&mut reader)?);
    }
    let evidence = decode_evidence(&mut reader)?;
    reader.finish()?;

    let study = IdentifiabilityStudySpec::try_new(
        header,
        context_of_use,
        model,
        initial_state,
        specimen,
        protocol,
        data,
        parameters,
        gauges,
        observations,
        paths,
        noise_dependence,
        discrepancies,
        evidence,
    )?;
    let canonical = study.canonical_bytes()?;
    if canonical != bytes {
        let at = canonical
            .iter()
            .zip(bytes)
            .position(|(left, right)| left != right)
            .unwrap_or(canonical.len().min(bytes.len()));
        return Err(IdentifiabilityError::Canonical {
            at,
            detail: "input is valid data but not the unique exact canonical encoding".to_string(),
        });
    }
    Ok(study)
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
