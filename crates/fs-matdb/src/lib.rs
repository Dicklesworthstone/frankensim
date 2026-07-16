//! fs-matdb: the L1 immutable cross-physics materials data layer (bead
//! frankensim-ext-matdb-core-5hmy, PR-1 of 5).
//!
//! "Real material properties" means a NAMED CONDITION with uncertainty,
//! not a marketing-grade label. Every downstream physics claim inherits
//! its weakest load-bearing material datum through the evidence system,
//! which only works if material data is TYPED (dims-checked at
//! insertion), IMMUTABLE (claims are never overwritten — conflicting
//! observations stay separate), PROVENANCE-COMPLETE (source, license,
//! and content addresses are load-bearing fields, not comments), and
//! queried through receipts (PR-4).
//!
//! PR-1 scope: [`ObservationDataset`], [`PropertyClaim`] with
//! `fs_evidence::ValidityDomain` integration (THE single validity type —
//! a competing type is forbidden), the [`PropertyKey`] dimension
//! registry, and the dims-checked, license-gated, fail-closed
//! [`ClaimSet`] insertion path. MaterialCard/ConstitutiveModelCard
//! (PR-2), InterfaceSystemCard (PR-3), the query/receipt path (PR-4) and
//! the receipt mutation battery (PR-5) follow.
//!
//! Layer: L1. Deps: fs-qty, fs-evidence, fs-blake3 ONLY. This crate owns
//! NO executable closures and NO per-run state; those belong to L3
//! adapters. It never imports L2 transforms, L3 state types, or L6
//! persistence.

use std::collections::BTreeMap;
use std::fmt;

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::ValidityDomain;
use fs_qty::Dims;

mod cards;
mod interface;
mod pack;
mod query;

pub use cards::{
    CANONICAL_PARAMETER_BLOCK_IDENTITY_DOMAIN, CANONICAL_PARAMETER_BLOCK_IDENTITY_VERSION,
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, MATDB_SCHEMA_VERSION,
    MaterialCard, MaterialStateId,
};
pub use interface::{InterfaceSystemCard, SurfaceSpec, SystemContext};
pub use pack::{
    JointStatistics, MATDB_PACK_SCHEMA_VERSION, MATDB_PACK_TARGET_BASIS, NormalizationReceipt,
    NormalizationTarget, NormalizedPack, PackError, StatisticComponent, StatisticMember,
    ValidityBoundSide,
};
pub use query::{
    EvaluationDecision, MATDB_EVALUATOR_VERSION, MaterialAnswer, PropertySample,
    PropertyUsageReceipt, QueryPoint, SelectionPolicy,
};

/// Hash domain for property-claim canonical identity.
const CLAIM_HASH_DOMAIN: &str = "org.frankensim.fs-matdb.property-claim.v1";
/// Hash domain for observation-dataset canonical identity.
const OBSERVATION_HASH_DOMAIN: &str = "org.frankensim.fs-matdb.observation-dataset.v1";

/// Everything that can go wrong at the immutable boundary. Total and
/// typed: refusals teach, and nothing inserts partially.
#[derive(Debug, Clone, PartialEq)]
pub enum MatDbError {
    /// A claim's value dimensions disagree with its key's registered
    /// dimensions.
    DimsMismatch {
        /// The property key.
        key: PropertyKey,
        /// The key's registered dimensions.
        expected: Dims,
        /// The offered value's dimensions.
        found: Dims,
    },
    /// The provenance record has an empty license field. License is
    /// LOAD-BEARING: unlicensed data cannot enter the store.
    MissingLicense {
        /// The claim/observation source citation (may itself be empty).
        source: String,
    },
    /// The provenance record has an empty source citation.
    MissingSource,
    /// A numeric payload field is non-finite where finiteness is
    /// structural (values, uncertainty widths, curve knots).
    NonFinite {
        /// Which field refused.
        field: &'static str,
        /// The offending bits (exact, for the receipt).
        bits: u64,
    },
    /// A validity-domain axis carries a NaN endpoint: the domain is
    /// unusable and the claim would be dead on arrival.
    UnusableValidity {
        /// The axis whose bounds are unusable.
        axis: String,
    },
    /// An uncertainty model parameter is structurally invalid (negative
    /// half-width, confidence outside (0,1), or empty sample count).
    InvalidUncertainty {
        /// What is wrong.
        reason: &'static str,
    },
    /// A curve payload has fewer than two knots, unordered abscissae, or
    /// duplicated abscissae.
    MalformedCurve {
        /// What is wrong.
        reason: &'static str,
    },
    /// A referenced observation id is not present in the claim set's
    /// observation registry.
    UnknownObservation {
        /// The dangling reference.
        observation: ObservationId,
    },
    /// A constitutive model card has no parameters: an empty block is a
    /// name, not a model.
    EmptyParameterBlock {
        /// The law whose card is empty.
        law: LawId,
    },
    /// A law parameter is non-finite.
    NonFiniteParameter {
        /// The law whose parameter refused.
        law: LawId,
        /// The parameter name.
        parameter: String,
        /// The offending bits (exact, for the receipt).
        bits: u64,
    },
    /// A genesis material card claimed a nonzero revision: lineage
    /// starts at 0 and only supersession advances it.
    RevisionNotZero {
        /// The offered revision.
        offered: u32,
    },
    /// A supersession is structurally impossible (named-state mismatch
    /// or exhausted revision counter).
    SupersedesMismatch {
        /// What is wrong.
        reason: &'static str,
    },
    /// An interface surface has no texture-frame id: a claim against an
    /// unnamed surface state is unreproducible.
    MissingTextureFrame {
        /// The surface's bulk material state.
        material: MaterialStateId,
    },
    /// A required interface-system field (medium, environment, history,
    /// or a present-but-blank third body) is empty.
    MissingSystemField {
        /// Which field.
        field: &'static str,
    },
    /// A query coordinate is non-finite.
    NonFiniteQueryPoint {
        /// The axis.
        axis: String,
        /// The offending bits.
        bits: u64,
    },
    /// No claim carries the queried property name.
    UnknownProperty {
        /// The name.
        property: String,
    },
    /// The point lies outside every claim's validity: THE extrapolation
    /// refusal. Nothing is evaluated.
    NoClaimInDomain {
        /// The property.
        property: String,
        /// How many claims were considered (all out of domain).
        considered: usize,
    },
    /// The selection policy cannot narrow the in-domain candidates to
    /// one claim; fusion must be explicit, so ambiguity refuses.
    AmbiguousSelection {
        /// The property.
        property: String,
        /// The surviving candidates.
        candidates: Vec<ClaimId>,
    },
    /// A curve claim's abscissa axis is missing from the query point.
    MissingQueryAxis {
        /// The required axis.
        axis: String,
    },
    /// The requested abscissa lies outside the curve's knot span
    /// (extrapolation beyond data refuses even inside validity).
    OutsideKnotSpan {
        /// The abscissa axis.
        axis: String,
        /// The requested coordinate.
        requested: f64,
        /// First knot.
        lo: f64,
        /// Last knot.
        hi: f64,
    },
    /// The claim's payload/policy combination cannot answer this query.
    UnsupportedEvaluation {
        /// What is wrong.
        reason: &'static str,
    },
    /// A receipt names a selection-policy tag no policy owns.
    UnknownPolicyTag {
        /// The foreign tag.
        tag: String,
    },
    /// A receipt was minted by a different evaluator version; it cannot
    /// be adjudicated by this one.
    EvaluatorVersionDrift {
        /// The receipt's evaluator version.
        receipt: u32,
        /// This evaluator's version.
        current: u32,
    },
    /// A replayed receipt field failed to reproduce: the receipt is
    /// incomplete, tampered, or stale against this claim set.
    ReceiptMismatch {
        /// The first divergent field.
        field: &'static str,
    },
}

impl fmt::Display for MatDbError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            MatDbError::DimsMismatch {
                key,
                expected,
                found,
            } => write!(
                f,
                "property '{}': value dims {found:?} disagree with the registered dims \
                 {expected:?}",
                key.name()
            ),
            MatDbError::MissingLicense { source } => write!(
                f,
                "provenance for '{source}' has no license: unlicensed data cannot enter the store"
            ),
            MatDbError::MissingSource => {
                f.write_str("provenance has no source citation; a bare value is not a datum")
            }
            MatDbError::NonFinite { field, bits } => {
                write!(f, "field '{field}' is non-finite (bits {bits:#018x})")
            }
            MatDbError::UnusableValidity { axis } => {
                write!(
                    f,
                    "validity axis '{axis}' has NaN bounds: the domain is unusable"
                )
            }
            MatDbError::InvalidUncertainty { reason } => {
                write!(f, "uncertainty model invalid: {reason}")
            }
            MatDbError::MalformedCurve { reason } => write!(f, "curve payload invalid: {reason}"),
            MatDbError::UnknownObservation { observation } => write!(
                f,
                "claim references unknown observation dataset {observation:?}"
            ),
            MatDbError::EmptyParameterBlock { law } => write!(
                f,
                "constitutive card for law '{}' has an empty parameter block",
                law.0
            ),
            MatDbError::NonFiniteParameter {
                law,
                parameter,
                bits,
            } => write!(
                f,
                "law '{}' parameter '{parameter}' is non-finite (bits {bits:#018x})",
                law.0
            ),
            MatDbError::RevisionNotZero { offered } => write!(
                f,
                "a genesis material card must be revision 0, not {offered}; use supersede"
            ),
            MatDbError::SupersedesMismatch { reason } => {
                write!(f, "supersession impossible: {reason}")
            }
            MatDbError::MissingTextureFrame { material } => write!(
                f,
                "surface of {material} has no texture-frame id; interface claims against an \
                 unnamed surface state are unreproducible"
            ),
            MatDbError::MissingSystemField { field } => write!(
                f,
                "interface system field '{field}' is blank; the system identity is incomplete"
            ),
            MatDbError::NonFiniteQueryPoint { axis, bits } => {
                write!(f, "query axis '{axis}' is non-finite (bits {bits:#018x})")
            }
            MatDbError::UnknownProperty { property } => {
                write!(f, "no claim carries property '{property}'")
            }
            MatDbError::NoClaimInDomain {
                property,
                considered,
            } => write!(
                f,
                "extrapolation refused: the query point is outside the validity of every one of \
                 the {considered} claim(s) for '{property}'"
            ),
            MatDbError::AmbiguousSelection {
                property,
                candidates,
            } => write!(
                f,
                "{} in-domain claims for '{property}' survive the policy; fusion must be an \
                 explicit policy, so ambiguity refuses",
                candidates.len()
            ),
            MatDbError::MissingQueryAxis { axis } => {
                write!(f, "the query point lacks required abscissa axis '{axis}'")
            }
            MatDbError::OutsideKnotSpan {
                axis,
                requested,
                lo,
                hi,
            } => write!(
                f,
                "abscissa '{axis}' = {requested} is outside the knot span [{lo}, {hi}]; \
                 extrapolation beyond data refuses"
            ),
            MatDbError::UnsupportedEvaluation { reason } => {
                write!(f, "unsupported evaluation: {reason}")
            }
            MatDbError::UnknownPolicyTag { tag } => {
                write!(f, "no selection policy owns receipt tag '{tag}'")
            }
            MatDbError::EvaluatorVersionDrift { receipt, current } => write!(
                f,
                "receipt evaluator version {receipt} cannot be adjudicated by evaluator {current}"
            ),
            MatDbError::ReceiptMismatch { field } => write!(
                f,
                "receipt field '{field}' failed to reproduce under replay: incomplete, tampered, \
                 or stale"
            ),
        }
    }
}

impl std::error::Error for MatDbError {}

/// A property key: the NAME of a material property plus its registered
/// dimensions. The registry is the dims-check authority at insertion —
/// a claim whose value dims disagree with its key refuses.
///
/// The name is free-form (the initial vocabulary in the bead is
/// expressible without a closed enum), but the (name, dims) pair is the
/// identity: `density` with mass/volume dims and `density` with any
/// other dims are DIFFERENT keys, and a [`ClaimSet`] refuses to register
/// the same name twice with different dims.
#[derive(Debug, Clone, PartialEq)]
pub struct PropertyKey {
    name: String,
    dims: Dims,
}

impl PropertyKey {
    /// A key with its registered dimensions.
    #[must_use]
    pub fn new(name: impl Into<String>, dims: Dims) -> PropertyKey {
        PropertyKey {
            name: name.into(),
            dims,
        }
    }

    /// The property name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// The registered dimensions.
    #[must_use]
    pub fn dims(&self) -> Dims {
        self.dims
    }
}

/// Provenance: WHERE a datum came from and under WHAT license. Both
/// fields are load-bearing; the content hash pins the exact acquired
/// artifact when one exists.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Provenance {
    /// Source citation (standard, paper, datasheet, lab notebook id).
    pub source: String,
    /// License / usage terms. EMPTY REFUSES: unlicensed data cannot
    /// enter the store, per the crate contract.
    pub license: String,
    /// Content address of the acquired artifact (report PDF, CSV,
    /// instrument export), when retained.
    pub artifact: Option<ContentHash>,
}

impl Provenance {
    pub(crate) fn validate(&self) -> Result<(), MatDbError> {
        if self.source.trim().is_empty() {
            return Err(MatDbError::MissingSource);
        }
        if self.license.trim().is_empty() {
            return Err(MatDbError::MissingLicense {
                source: self.source.clone(),
            });
        }
        Ok(())
    }
}

/// Identifier of a registered [`ObservationDataset`] (its content hash).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct ObservationId(pub ContentHash);

/// Identifier of an inserted [`PropertyClaim`] (its content hash).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct ClaimId(pub ContentHash);

/// One observed dataset: the specimen/process context, the method, the
/// raw artifact, and its statistical caveats. Observations are the
/// ground truth CLAIMS point at — a claim without observations is
/// citation-only and can never be more than Estimated-class evidence
/// downstream (the query path enforces the color discipline in PR-4;
/// PR-1 stores the linkage).
#[derive(Debug, Clone, PartialEq)]
pub struct ObservationDataset {
    /// Specimen and process record: alloy/temper, cure state, print
    /// orientation — the named condition the data is true OF.
    pub specimen: String,
    /// Method and instrument (standard designation where one exists,
    /// e.g. "ASTM E8 / frame X").
    pub method: String,
    /// Content address of the observation artifact (raw table, curve).
    pub artifact: ContentHash,
    /// Covariance/censoring notes: what was measured jointly, what was
    /// truncated or detection-limited. Free text in PR-1; typed in the
    /// PR-4 receipt path.
    pub caveats: String,
    /// Where the dataset came from and under what license.
    pub provenance: Provenance,
}

impl ObservationDataset {
    /// Canonical content identity over every field.
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        let mut payload = Vec::new();
        for part in [
            self.specimen.as_bytes(),
            self.method.as_bytes(),
            &self.artifact.0,
            self.caveats.as_bytes(),
            self.provenance.source.as_bytes(),
            self.provenance.license.as_bytes(),
        ] {
            payload.extend_from_slice(&u64::try_from(part.len()).unwrap_or(u64::MAX).to_le_bytes());
            payload.extend_from_slice(part);
        }
        if let Some(artifact) = &self.provenance.artifact {
            payload.extend_from_slice(&artifact.0);
        }
        hash_domain(OBSERVATION_HASH_DOMAIN, &payload)
    }
}

/// The value payload of a property claim. PR-1 carries the scalar and
/// curve forms (enough for the density/moduli/conductivity class of the
/// initial vocabulary); tensor, distribution, and model-parameter
/// payloads land with the constitutive card work (PR-2) so their frames
/// and state schemas arrive together.
#[derive(Debug, Clone, PartialEq)]
pub enum PropertyValue {
    /// A single SI value with its dimensions.
    Scalar {
        /// SI base-unit value.
        value: f64,
        /// The value's dimensions.
        dims: Dims,
    },
    /// A 1-D curve `property(abscissa)` — e.g. conductivity(T), BH
    /// curves, S-N data — with strictly increasing abscissae.
    Curve {
        /// The abscissa axis name (must match a validity axis name when
        /// the claim constrains it, e.g. "T").
        abscissa: String,
        /// The abscissa's dimensions.
        abscissa_dims: Dims,
        /// `(x, y)` knots, strictly increasing in `x`, all finite.
        knots: Vec<(f64, f64)>,
        /// The ordinate's dimensions.
        dims: Dims,
    },
}

impl PropertyValue {
    /// The payload's ordinate dimensions (what the dims check compares
    /// against the key's registration).
    #[must_use]
    pub fn dims(&self) -> Dims {
        match self {
            PropertyValue::Scalar { dims, .. } | PropertyValue::Curve { dims, .. } => *dims,
        }
    }

    fn validate(&self) -> Result<(), MatDbError> {
        match self {
            PropertyValue::Scalar { value, .. } => {
                if !value.is_finite() {
                    return Err(MatDbError::NonFinite {
                        field: "scalar value",
                        bits: value.to_bits(),
                    });
                }
                Ok(())
            }
            PropertyValue::Curve { knots, .. } => {
                if knots.len() < 2 {
                    return Err(MatDbError::MalformedCurve {
                        reason: "a curve needs at least two knots",
                    });
                }
                for &(x, y) in knots {
                    if !x.is_finite() {
                        return Err(MatDbError::NonFinite {
                            field: "curve abscissa",
                            bits: x.to_bits(),
                        });
                    }
                    if !y.is_finite() {
                        return Err(MatDbError::NonFinite {
                            field: "curve ordinate",
                            bits: y.to_bits(),
                        });
                    }
                }
                if !knots.windows(2).all(|w| w[0].0 < w[1].0) {
                    return Err(MatDbError::MalformedCurve {
                        reason: "abscissae must be strictly increasing",
                    });
                }
                Ok(())
            }
        }
    }
}

/// The uncertainty model attached to a claim. Uncertainty is not
/// optional decoration: a claim with `Unstated` uncertainty is admitted
/// (honesty over refusal — the datum exists) but is marked, and the
/// PR-4 query path will never let it launder into a certified band.
#[derive(Debug, Clone, PartialEq)]
pub enum UncertaintyModel {
    /// The source states no uncertainty. Admitted and marked.
    Unstated,
    /// Symmetric absolute half-width at a stated confidence.
    HalfWidth {
        /// Half-width in the value's SI units (finite, `>= 0`).
        half_width: f64,
        /// Confidence level, strictly inside `(0, 1)`.
        confidence: f64,
    },
    /// Relative (fractional) half-width at a stated confidence.
    RelativeHalfWidth {
        /// Fractional half-width (finite, `>= 0`).
        fraction: f64,
        /// Confidence level, strictly inside `(0, 1)`.
        confidence: f64,
    },
}

impl UncertaintyModel {
    fn validate(&self) -> Result<(), MatDbError> {
        match *self {
            UncertaintyModel::Unstated => Ok(()),
            UncertaintyModel::HalfWidth {
                half_width,
                confidence,
            } => {
                if !half_width.is_finite() || half_width < 0.0 {
                    return Err(MatDbError::InvalidUncertainty {
                        reason: "half-width must be finite and non-negative",
                    });
                }
                validate_confidence(confidence)
            }
            UncertaintyModel::RelativeHalfWidth {
                fraction,
                confidence,
            } => {
                if !fraction.is_finite() || fraction < 0.0 {
                    return Err(MatDbError::InvalidUncertainty {
                        reason: "relative half-width must be finite and non-negative",
                    });
                }
                validate_confidence(confidence)
            }
        }
    }
}

fn validate_confidence(confidence: f64) -> Result<(), MatDbError> {
    if confidence.is_finite() && confidence > 0.0 && confidence < 1.0 {
        Ok(())
    } else {
        Err(MatDbError::InvalidUncertainty {
            reason: "confidence must be strictly between 0 and 1",
        })
    }
}

/// How a consumer may evaluate the claim BETWEEN its knots or ACROSS its
/// validity box. Extrapolation is never implicit: the PR-4 query path
/// turns an out-of-domain evaluation into a typed refusal or an
/// explicit, receipt-recorded extrapolation decision.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InterpolationPolicy {
    /// Piecewise-linear between knots; refuse outside the knot span.
    LinearInside,
    /// The value holds across the whole validity box (a plateau claim).
    ConstantWithinValidity,
    /// No interpolation claim at all: exact tabulated points only.
    TabulatedOnly,
}

/// One immutable property claim: a value, where it is valid, how
/// uncertain it is, and where it came from. Conflicting claims for the
/// same key COEXIST — fusion is an explicit query-time policy (PR-4),
/// never a map overwrite that invents a canonical value.
#[derive(Debug, Clone, PartialEq)]
pub struct PropertyClaim {
    /// The property being claimed.
    pub key: PropertyKey,
    /// The value payload (dims must match the key's registration).
    pub value: PropertyValue,
    /// The region of condition space the claim is good for, over typed
    /// axes (T, p, frequency, field, strain rate, composition, history,
    /// named dimensionless groups). THE single validity type, reused
    /// from fs-evidence.
    pub validity: ValidityDomain,
    /// The claim's uncertainty model.
    pub uncertainty: UncertaintyModel,
    /// How the claim may be evaluated between/around its data.
    pub interpolation: InterpolationPolicy,
    /// Observation datasets backing the claim (may be empty: a
    /// citation-only claim is admitted but can never be Validated-class
    /// downstream — specimen/process match requires observations).
    pub observations: Vec<ObservationId>,
    /// Where the claim came from and under what license.
    pub provenance: Provenance,
}

impl PropertyClaim {
    /// Canonical content identity over every semantic field.
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        let mut payload = Vec::new();
        let mut push = |part: &[u8]| {
            payload.extend_from_slice(&u64::try_from(part.len()).unwrap_or(u64::MAX).to_le_bytes());
            payload.extend_from_slice(part);
        };
        push(self.key.name.as_bytes());
        push(&dims_bytes(self.key.dims));
        match &self.value {
            PropertyValue::Scalar { value, dims } => {
                push(b"scalar");
                push(&value.to_bits().to_le_bytes());
                push(&dims_bytes(*dims));
            }
            PropertyValue::Curve {
                abscissa,
                abscissa_dims,
                knots,
                dims,
            } => {
                push(b"curve");
                push(abscissa.as_bytes());
                push(&dims_bytes(*abscissa_dims));
                for &(x, y) in knots {
                    push(&x.to_bits().to_le_bytes());
                    push(&y.to_bits().to_le_bytes());
                }
                push(&dims_bytes(*dims));
            }
        }
        for (axis, &(lo, hi)) in self.validity.bounds() {
            push(axis.as_bytes());
            push(&lo.to_bits().to_le_bytes());
            push(&hi.to_bits().to_le_bytes());
        }
        match self.uncertainty {
            UncertaintyModel::Unstated => push(b"unstated"),
            UncertaintyModel::HalfWidth {
                half_width,
                confidence,
            } => {
                push(b"half-width");
                push(&half_width.to_bits().to_le_bytes());
                push(&confidence.to_bits().to_le_bytes());
            }
            UncertaintyModel::RelativeHalfWidth {
                fraction,
                confidence,
            } => {
                push(b"relative-half-width");
                push(&fraction.to_bits().to_le_bytes());
                push(&confidence.to_bits().to_le_bytes());
            }
        }
        push(match self.interpolation {
            InterpolationPolicy::LinearInside => b"linear-inside".as_slice(),
            InterpolationPolicy::ConstantWithinValidity => b"constant-within-validity",
            InterpolationPolicy::TabulatedOnly => b"tabulated-only",
        });
        for observation in &self.observations {
            push(&observation.0.0);
        }
        push(self.provenance.source.as_bytes());
        push(self.provenance.license.as_bytes());
        if let Some(artifact) = &self.provenance.artifact {
            push(&artifact.0);
        }
        hash_domain(CLAIM_HASH_DOMAIN, &payload)
    }
}

pub(crate) fn dims_bytes(dims: Dims) -> Vec<u8> {
    dims.0.iter().map(|&e| e.cast_unsigned()).collect()
}

/// The PR-1 immutable container: registered observations plus property
/// claims, keyed by content identity. Insertion is the ONLY mutation,
/// and it is fail-closed and append-only:
///
/// - a key registers its dims on first use; a later claim reusing the
///   name with different dims refuses ([`MatDbError::DimsMismatch`]);
/// - claims referencing unregistered observations refuse;
/// - unlicensed or source-less provenance refuses;
/// - conflicting claims for one key ALL coexist under distinct
///   [`ClaimId`]s — nothing is overwritten, ever.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ClaimSet {
    observations: BTreeMap<ObservationId, ObservationDataset>,
    claims: BTreeMap<ClaimId, PropertyClaim>,
    key_dims: BTreeMap<String, Dims>,
    by_key: BTreeMap<String, Vec<ClaimId>>,
}

impl ClaimSet {
    /// An empty set.
    #[must_use]
    pub fn new() -> ClaimSet {
        ClaimSet::default()
    }

    /// Register an observation dataset. Idempotent by content identity.
    ///
    /// # Errors
    /// Provenance refusals ([`MatDbError::MissingLicense`],
    /// [`MatDbError::MissingSource`]).
    pub fn register_observation(
        &mut self,
        dataset: ObservationDataset,
    ) -> Result<ObservationId, MatDbError> {
        dataset.provenance.validate()?;
        let id = ObservationId(dataset.content_hash());
        self.observations.entry(id).or_insert(dataset);
        Ok(id)
    }

    /// Insert a property claim. Append-only: the same content inserts
    /// idempotently, DIFFERENT content for the same key coexists.
    ///
    /// # Errors
    /// Every gate is typed: dims mismatch against the key registry,
    /// provenance refusals, non-finite payloads, malformed curves,
    /// invalid uncertainty, NaN validity axes, dangling observation
    /// references.
    pub fn insert_claim(&mut self, claim: PropertyClaim) -> Result<ClaimId, MatDbError> {
        claim.provenance.validate()?;
        claim.value.validate()?;
        claim.uncertainty.validate()?;
        for (axis, &(lo, hi)) in claim.validity.bounds() {
            if lo.is_nan() || hi.is_nan() {
                return Err(MatDbError::UnusableValidity { axis: axis.clone() });
            }
        }
        let found = claim.value.dims();
        if found != claim.key.dims() {
            return Err(MatDbError::DimsMismatch {
                key: claim.key.clone(),
                expected: claim.key.dims(),
                found,
            });
        }
        if let Some(&registered) = self.key_dims.get(claim.key.name())
            && registered != claim.key.dims()
        {
            return Err(MatDbError::DimsMismatch {
                key: claim.key.clone(),
                expected: registered,
                found,
            });
        }
        for observation in &claim.observations {
            if !self.observations.contains_key(observation) {
                return Err(MatDbError::UnknownObservation {
                    observation: *observation,
                });
            }
        }
        let id = ClaimId(claim.content_hash());
        self.key_dims
            .insert(claim.key.name().to_string(), claim.key.dims());
        let name = claim.key.name().to_string();
        if !self.claims.contains_key(&id) {
            self.claims.insert(id, claim);
            self.by_key.entry(name).or_default().push(id);
        }
        Ok(id)
    }

    /// A registered observation by id.
    #[must_use]
    pub fn observation(&self, id: ObservationId) -> Option<&ObservationDataset> {
        self.observations.get(&id)
    }

    /// An inserted claim by id.
    #[must_use]
    pub fn claim(&self, id: ClaimId) -> Option<&PropertyClaim> {
        self.claims.get(&id)
    }

    /// EVERY claim for a property name, in insertion order. Conflicting
    /// observations stay separate — this is the anti-laundering surface
    /// the PR-4 fusion policies will consume.
    #[must_use]
    pub fn claims_for(&self, name: &str) -> Vec<(ClaimId, &PropertyClaim)> {
        self.by_key
            .get(name)
            .into_iter()
            .flatten()
            .filter_map(|id| self.claims.get(id).map(|claim| (*id, claim)))
            .collect()
    }

    /// The dims registered for a property name, if any claim used it.
    #[must_use]
    pub fn registered_dims(&self, name: &str) -> Option<Dims> {
        self.key_dims.get(name).copied()
    }

    /// Number of stored claims.
    #[must_use]
    pub fn claim_count(&self) -> usize {
        self.claims.len()
    }

    /// All claims in content-id order (the canonical order card hashes
    /// bind).
    pub fn claims_ordered(&self) -> impl Iterator<Item = (ClaimId, &PropertyClaim)> {
        self.claims.iter().map(|(id, claim)| (*id, claim))
    }

    /// All registered observation ids in content-id order.
    pub fn observation_ids(&self) -> impl Iterator<Item = ObservationId> + '_ {
        self.observations.keys().copied()
    }
}
