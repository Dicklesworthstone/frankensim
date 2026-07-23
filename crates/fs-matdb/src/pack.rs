//! Canonical, runtime-loadable material-pack artifacts.
//!
//! Raw handbook/CSV/NASA parsing and redistribution policy deliberately do
//! **not** live here.  This L1 module owns only the bounded normalized wire
//! format that an offline compiler may emit after those policy decisions.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::ValidityDomain;
use fs_qty::Dims;

use crate::{
    ClaimId, ClaimSet, InterpolationPolicy, MatDbError, ObservationDataset, ObservationId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, UncertaintyModel,
};

/// Current normalized material-pack wire schema.
pub const MATDB_PACK_SCHEMA_VERSION: u32 = 1;
/// Canonical coherent-SI target basis and its explicit six-base order.
pub const MATDB_PACK_TARGET_BASIS: &str = "SI-six-base[m,kg,s,K,A,mol]";

const MAGIC: &[u8; 8] = b"FSMATPK\0";
const PACK_HASH_DOMAIN: &str = "org.frankensim.fs-matdb.normalized-pack.v1";
const MAX_PACK_BYTES: usize = 256 * 1024 * 1024;
const MAX_STRING_BYTES: usize = 1_048_576;
const MAX_OBSERVATIONS: usize = 100_000;
const MAX_CLAIMS: usize = 100_000;
const MAX_VALIDITY_AXES: usize = 4_096;
const MAX_CURVE_KNOTS: usize = 4_000_000;
const MAX_OBSERVATIONS_PER_CLAIM: usize = 100_000;
const MAX_JOINT_BLOCKS: usize = 100_000;
// PSD admission is cubic.  Keep the public block cap small enough that a
// hostile pack cannot turn a bounded decode into an unbounded CPU event.
const MAX_JOINT_MEMBERS: usize = 256;
const MAX_NORMALIZATIONS: usize = 100_000;
const MAX_PSD_CUBIC_WORK: u64 = 134_217_728;
const STATISTIC_MEMBER_BYTES: usize = 33;
const OBSERVATION_ID_BYTES: usize = 32;
const MIN_JOINT_BLOCK_BYTES: usize = 41;
const MIN_NORMALIZATION_BYTES: usize = 97;

/// Numeric component of a property claim participating in joint statistics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum StatisticComponent {
    /// The value of a scalar claim.
    Scalar,
    /// One curve knot's abscissa coordinate.
    CurveAbscissa {
        /// Zero-based knot index.
        knot: u32,
    },
    /// One curve knot's ordinate value.
    CurveOrdinate {
        /// Zero-based knot index.
        knot: u32,
    },
}

/// Exact scalar/curve component named by a covariance matrix row or column.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct StatisticMember {
    claim: ClaimId,
    component: StatisticComponent,
}

impl StatisticMember {
    /// The value of a scalar claim.
    #[must_use]
    pub const fn scalar(claim: ClaimId) -> Self {
        Self {
            claim,
            component: StatisticComponent::Scalar,
        }
    }

    /// One curve knot's abscissa coordinate.
    #[must_use]
    pub const fn curve_abscissa(claim: ClaimId, knot: u32) -> Self {
        Self {
            claim,
            component: StatisticComponent::CurveAbscissa { knot },
        }
    }

    /// One curve knot's ordinate value.
    #[must_use]
    pub const fn curve_ordinate(claim: ClaimId, knot: u32) -> Self {
        Self {
            claim,
            component: StatisticComponent::CurveOrdinate { knot },
        }
    }

    /// Owning property claim.
    #[must_use]
    pub const fn claim(self) -> ClaimId {
        self.claim
    }

    /// Numeric component within the claim.
    #[must_use]
    pub const fn component(self) -> StatisticComponent {
        self.component
    }
}

/// Which endpoint of a normalized validity interval is addressed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ValidityBoundSide {
    /// Inclusive lower endpoint.
    Lower,
    /// Inclusive upper endpoint.
    Upper,
}

/// Exact normalized numeric field addressed by a transform receipt.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum NormalizationTarget {
    /// Scalar value or one curve-knot coordinate/value.
    ClaimValue(StatisticMember),
    /// Absolute/relative uncertainty parameter on a claim.
    ClaimUncertainty {
        /// Owning claim.
        claim: ClaimId,
    },
    /// One endpoint of a named validity axis.
    ValidityBound {
        /// Owning claim.
        claim: ClaimId,
        /// Existing validity-axis name.
        axis: String,
        /// Lower or upper endpoint.
        side: ValidityBoundSide,
    },
    /// One packed covariance entry in a named joint block.
    JointCovariance {
        /// Owning observation.
        observation: ObservationId,
        /// Stable block name within the observation.
        block_id: String,
        /// Matrix row in explicit member order.
        row: u32,
        /// Matrix column in explicit member order (`column <= row`).
        column: u32,
    },
}

/// A typed refusal at the normalized-pack boundary.
#[derive(Debug, Clone, PartialEq)]
pub enum PackError {
    /// A semantic field is absent or violates the canonical pack profile.
    InvalidField {
        /// Stable field/rule identity.
        field: &'static str,
        /// Teaching detail.
        detail: String,
    },
    /// A collection or byte input exceeded its public processing budget.
    ResourceLimit {
        /// Stable resource identity.
        resource: &'static str,
        /// Configured maximum.
        limit: usize,
        /// Exact observation or proven lower bound.
        observed: usize,
    },
    /// The binary envelope is malformed or non-canonical.
    Malformed {
        /// Byte offset at which decoding refused.
        at: usize,
        /// Stable diagnostic detail.
        detail: String,
    },
    /// A decoded semantic object failed the ordinary fs-matdb admission gate.
    MatDb(MatDbError),
    /// A serialized semantic id did not reproduce after decoding.
    IdentityMismatch {
        /// Object class whose identity moved.
        kind: &'static str,
        /// Serialized identity.
        expected: ContentHash,
        /// Reconstructed identity.
        actual: ContentHash,
    },
}

impl fmt::Display for PackError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidField { field, detail } => {
                write!(f, "normalized pack field '{field}' refused: {detail}")
            }
            Self::ResourceLimit {
                resource,
                limit,
                observed,
            } => write!(
                f,
                "normalized pack resource '{resource}' exceeds {limit} (observed {observed})"
            ),
            Self::Malformed { at, detail } => {
                write!(f, "malformed normalized pack at byte {at}: {detail}")
            }
            Self::MatDb(error) => write!(f, "normalized pack semantic admission failed: {error}"),
            Self::IdentityMismatch {
                kind,
                expected,
                actual,
            } => write!(
                f,
                "normalized pack {kind} identity mismatch: encoded {expected}, reconstructed {actual}"
            ),
        }
    }
}

impl std::error::Error for PackError {}

impl From<MatDbError> for PackError {
    fn from(value: MatDbError) -> Self {
        Self::MatDb(value)
    }
}

/// Joint covariance/correlation for one observation dataset.
///
/// `members` is an explicit, strictly increasing list of scalar values or
/// curve-knot components. Both matrices use packed lower-triangle order
/// `(0,0), (1,0), (1,1), ...`.
/// Covariance entries are coherent-SI values whose dimensions are implied by
/// the product of the corresponding component dimensions. Correlation, when
/// present, is dimensionless and has an exact unit diagonal.
#[derive(Debug, Clone, PartialEq)]
pub struct JointStatistics {
    observation: ObservationId,
    block_id: String,
    members: Vec<StatisticMember>,
    covariance: Vec<f64>,
    correlation: Option<Vec<f64>>,
}

impl JointStatistics {
    /// Construct a typed joint-statistics block.  Full validation against a
    /// claim set occurs when the block enters [`NormalizedPack`].
    #[must_use]
    pub fn new(
        observation: ObservationId,
        block_id: impl Into<String>,
        members: Vec<StatisticMember>,
        covariance: Vec<f64>,
        correlation: Option<Vec<f64>>,
    ) -> Self {
        Self {
            observation,
            block_id: block_id.into(),
            members,
            covariance,
            correlation,
        }
    }

    /// Observation dataset whose joint samples supplied the matrices.
    #[must_use]
    pub fn observation(&self) -> ObservationId {
        self.observation
    }

    /// Stable source-defined joint block name within the observation.
    #[must_use]
    pub fn block_id(&self) -> &str {
        &self.block_id
    }

    /// Ordered matrix member claims.
    #[must_use]
    pub fn members(&self) -> &[StatisticMember] {
        &self.members
    }

    /// Packed lower-triangle covariance in coherent SI.
    #[must_use]
    pub fn covariance(&self) -> &[f64] {
        &self.covariance
    }

    /// Packed lower-triangle correlation, when the source supplied it.
    #[must_use]
    pub fn correlation(&self) -> Option<&[f64]> {
        self.correlation.as_deref()
    }
}

/// Auditable unit/basis normalization applied by the offline compiler.
///
/// The source literal itself need not be redistributable, so only its content
/// hash crosses the pack boundary.  `si = source * scale + offset` records the
/// exact affine transform used for this field.  A frame pair is metadata only:
/// this scalar/curve schema makes no tensor-rotation claim.
#[derive(Debug, Clone, PartialEq)]
pub struct NormalizationReceipt {
    target: NormalizationTarget,
    source_literal: ContentHash,
    dims: Dims,
    scale: f64,
    offset: f64,
    source_basis: String,
    target_basis: String,
    source_frame: Option<String>,
    target_frame: Option<String>,
}

impl NormalizationReceipt {
    /// Construct an immutable transform receipt.
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub fn new(
        target: NormalizationTarget,
        source_literal: ContentHash,
        dims: Dims,
        scale: f64,
        offset: f64,
        source_basis: impl Into<String>,
        target_basis: impl Into<String>,
        source_frame: Option<String>,
        target_frame: Option<String>,
    ) -> Self {
        Self {
            target,
            source_literal,
            dims,
            scale,
            offset,
            source_basis: source_basis.into(),
            target_basis: target_basis.into(),
            source_frame,
            target_frame,
        }
    }

    /// Structurally linked normalized field.
    #[must_use]
    pub fn target(&self) -> &NormalizationTarget {
        &self.target
    }

    /// Hash of the exact source literal.
    #[must_use]
    pub fn source_literal(&self) -> ContentHash {
        self.source_literal
    }

    /// Six-base dimensions of the normalized field.
    #[must_use]
    pub fn dims(&self) -> Dims {
        self.dims
    }

    /// Multiplicative term in `si = source * scale + offset`.
    #[must_use]
    pub fn scale(&self) -> f64 {
        self.scale
    }

    /// Additive term in `si = source * scale + offset`.
    #[must_use]
    pub fn offset(&self) -> f64 {
        self.offset
    }

    /// Declared source basis.
    #[must_use]
    pub fn source_basis(&self) -> &str {
        &self.source_basis
    }

    /// Declared normalized basis.
    #[must_use]
    pub fn target_basis(&self) -> &str {
        &self.target_basis
    }

    /// Optional source frame identifier.
    #[must_use]
    pub fn source_frame(&self) -> Option<&str> {
        self.source_frame.as_deref()
    }

    /// Optional target frame identifier.
    #[must_use]
    pub fn target_frame(&self) -> Option<&str> {
        self.target_frame.as_deref()
    }
}

/// Runtime-loadable result of an admitted offline material-pack compilation.
#[derive(Debug, Clone, PartialEq)]
pub struct NormalizedPack {
    pack_id: String,
    compiler: String,
    source_artifact: ContentHash,
    redistribution_terms: String,
    claims: ClaimSet,
    joint_statistics: Vec<JointStatistics>,
    normalizations: Vec<NormalizationReceipt>,
}

impl NormalizedPack {
    /// Admit a pack into the canonical runtime profile.
    ///
    /// The constructor validates references, matrix shape/PSD, canonical
    /// ordering, finite portable values, signed-zero policy, and transform
    /// receipts.  License interpretation remains the offline compiler's job;
    /// this boundary requires its nonempty retained redistribution decision.
    #[allow(clippy::needless_pass_by_value)] // by-value admission: the caller surrenders its set so no non-canonical alias survives
    pub fn new(
        pack_id: impl Into<String>,
        compiler: impl Into<String>,
        source_artifact: ContentHash,
        redistribution_terms: impl Into<String>,
        claims: ClaimSet,
        mut joint_statistics: Vec<JointStatistics>,
        mut normalizations: Vec<NormalizationReceipt>,
    ) -> Result<Self, PackError> {
        let pack_id = pack_id.into();
        let compiler = compiler.into();
        let redistribution_terms = redistribution_terms.into();
        require_text("pack_id", &pack_id)?;
        require_text("compiler", &compiler)?;
        require_text("redistribution_terms", &redistribution_terms)?;
        validate_claim_set(&claims)?;
        let claims = canonical_claim_set(&claims)?;

        if joint_statistics.len() > MAX_JOINT_BLOCKS {
            return Err(limit(
                "joint_statistics",
                MAX_JOINT_BLOCKS,
                joint_statistics.len(),
            ));
        }
        joint_statistics.sort_by(|left, right| {
            (left.observation, left.block_id.as_str())
                .cmp(&(right.observation, right.block_id.as_str()))
        });
        for pair in joint_statistics.windows(2) {
            if pair[0].observation == pair[1].observation && pair[0].block_id == pair[1].block_id {
                return Err(invalid(
                    "joint_statistics",
                    format!(
                        "duplicate block {:?} for observation {}",
                        pair[0].block_id, pair[0].observation.0
                    ),
                ));
            }
        }
        validate_psd_work_budget(&joint_statistics)?;
        for block in &joint_statistics {
            validate_joint_statistics(&claims, block)?;
        }
        validate_disjoint_joint_statistics(&joint_statistics)?;

        if normalizations.len() > MAX_NORMALIZATIONS {
            return Err(limit(
                "normalizations",
                MAX_NORMALIZATIONS,
                normalizations.len(),
            ));
        }
        normalizations.sort_by(|left, right| left.target.cmp(&right.target));
        for pair in normalizations.windows(2) {
            if pair[0].target == pair[1].target {
                return Err(invalid(
                    "normalizations",
                    format!("duplicate target receipt {:?}", pair[0].target),
                ));
            }
        }
        for receipt in &normalizations {
            validate_normalization(&claims, &joint_statistics, receipt)?;
        }
        validate_normalization_coherence(&normalizations)?;

        let pack = Self {
            pack_id,
            compiler,
            source_artifact,
            redistribution_terms,
            claims,
            joint_statistics,
            normalizations,
        };
        let encoded_bytes = pack.to_bytes().len();
        if encoded_bytes > MAX_PACK_BYTES {
            return Err(limit("pack_bytes", MAX_PACK_BYTES, encoded_bytes));
        }
        Ok(pack)
    }

    /// Stable pack name supplied by the source manifest.
    #[must_use]
    pub fn pack_id(&self) -> &str {
        &self.pack_id
    }

    /// Compiler/version identity that made the admission decisions.
    #[must_use]
    pub fn compiler(&self) -> &str {
        &self.compiler
    }

    /// Hash of the exact raw source envelope.
    #[must_use]
    pub fn source_artifact(&self) -> ContentHash {
        self.source_artifact
    }

    /// Retained redistribution decision/terms.
    #[must_use]
    pub fn redistribution_terms(&self) -> &str {
        &self.redistribution_terms
    }

    /// Immutable fs-matdb claims reconstructed at runtime.
    #[must_use]
    pub fn claims(&self) -> &ClaimSet {
        &self.claims
    }

    /// Canonically ordered joint-statistics blocks.
    #[must_use]
    pub fn joint_statistics(&self) -> &[JointStatistics] {
        &self.joint_statistics
    }

    /// Canonically ordered unit/basis transform receipts.
    #[must_use]
    pub fn normalizations(&self) -> &[NormalizationReceipt] {
        &self.normalizations
    }

    /// Canonical binary representation consumed by L1.
    #[must_use]
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut writer = Writer::default();
        writer.bytes.extend_from_slice(MAGIC);
        writer.u32(MATDB_PACK_SCHEMA_VERSION);
        writer.string(&self.pack_id);
        writer.string(&self.compiler);
        writer.hash(self.source_artifact);
        writer.string(&self.redistribution_terms);

        let observations: Vec<_> = self
            .claims
            .observation_ids()
            .filter_map(|id| self.claims.observation(id).map(|dataset| (id, dataset)))
            .collect();
        writer.count(observations.len());
        for (id, dataset) in observations {
            writer.hash(id.0);
            encode_observation(&mut writer, dataset);
        }

        let claims: Vec<_> = self.claims.claims_ordered().collect();
        writer.count(claims.len());
        for (id, claim) in claims {
            writer.hash(id.0);
            encode_claim(&mut writer, claim);
        }

        writer.count(self.joint_statistics.len());
        for block in &self.joint_statistics {
            writer.hash(block.observation.0);
            writer.string(&block.block_id);
            writer.count(block.members.len());
            for member in &block.members {
                encode_statistic_member(&mut writer, *member);
            }
            writer.f64s(&block.covariance);
            match &block.correlation {
                None => writer.u8(0),
                Some(values) => {
                    writer.u8(1);
                    writer.f64s(values);
                }
            }
        }

        writer.count(self.normalizations.len());
        for receipt in &self.normalizations {
            encode_normalization_target(&mut writer, &receipt.target);
            writer.hash(receipt.source_literal);
            writer.dims(receipt.dims);
            writer.f64(receipt.scale);
            writer.f64(receipt.offset);
            writer.string(&receipt.source_basis);
            writer.string(&receipt.target_basis);
            writer.optional_string(receipt.source_frame.as_deref());
            writer.optional_string(receipt.target_frame.as_deref());
        }
        writer.bytes
    }

    /// Domain-separated identity of the canonical pack bytes.
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        hash_domain(PACK_HASH_DOMAIN, &self.to_bytes())
    }

    /// Verify an externally pinned artifact identity before decoding.
    ///
    /// Nested observation/claim ids protect semantic reconstruction, while
    /// this expected whole-pack hash is the authority for top-level metadata,
    /// matrices, and normalization receipts.
    pub fn from_bytes_verified(expected: ContentHash, bytes: &[u8]) -> Result<Self, PackError> {
        if bytes.len() > MAX_PACK_BYTES {
            return Err(limit("pack_bytes", MAX_PACK_BYTES, bytes.len()));
        }
        let actual = hash_domain(PACK_HASH_DOMAIN, bytes);
        if actual != expected {
            return Err(PackError::IdentityMismatch {
                kind: "pack",
                expected,
                actual,
            });
        }
        Self::from_bytes(bytes)
    }

    /// Decode and semantically re-admit a canonical normalized pack.
    #[allow(clippy::too_many_lines)] // one canonical decode sequence: splitting would hide the wire order
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, PackError> {
        if bytes.len() > MAX_PACK_BYTES {
            return Err(limit("pack_bytes", MAX_PACK_BYTES, bytes.len()));
        }
        let mut reader = Reader::new(bytes);
        reader.expect(MAGIC, "normalized pack magic")?;
        let version = reader.u32()?;
        if version != MATDB_PACK_SCHEMA_VERSION {
            return Err(reader.malformed(format!(
                "unsupported schema version {version}; expected {MATDB_PACK_SCHEMA_VERSION}"
            )));
        }
        let pack_id = reader.string()?;
        let compiler = reader.string()?;
        let source_artifact = reader.hash()?;
        let redistribution_terms = reader.string()?;

        let observation_count = reader.count("observations", MAX_OBSERVATIONS)?;
        let mut claims = ClaimSet::new();
        for _ in 0..observation_count {
            let expected = reader.hash()?;
            let dataset = decode_observation(&mut reader)?;
            let actual = claims.register_observation(dataset)?.0;
            if actual != expected {
                return Err(PackError::IdentityMismatch {
                    kind: "observation",
                    expected,
                    actual,
                });
            }
        }

        let claim_count = reader.count("claims", MAX_CLAIMS)?;
        for _ in 0..claim_count {
            let expected = reader.hash()?;
            let claim = decode_claim(&mut reader)?;
            let actual = claims.insert_claim(claim)?.0;
            if actual != expected {
                return Err(PackError::IdentityMismatch {
                    kind: "claim",
                    expected,
                    actual,
                });
            }
        }

        let block_count = reader.count("joint_statistics", MAX_JOINT_BLOCKS)?;
        reader.require_items(
            block_count,
            MIN_JOINT_BLOCK_BYTES,
            "joint-statistics blocks",
        )?;
        let mut joint_statistics = Vec::with_capacity(block_count);
        for _ in 0..block_count {
            let observation = ObservationId(reader.hash()?);
            let block_id = reader.string()?;
            let member_count = reader.count("joint_members", MAX_JOINT_MEMBERS)?;
            reader.require_items(
                member_count,
                STATISTIC_MEMBER_BYTES,
                "joint-statistics members",
            )?;
            let mut members = Vec::with_capacity(member_count);
            for _ in 0..member_count {
                members.push(decode_statistic_member(&mut reader)?);
            }
            let triangle = triangle_len(member_count)?;
            let covariance = reader.fixed_f64s(triangle)?;
            let correlation = match reader.u8()? {
                0 => None,
                1 => Some(reader.fixed_f64s(triangle)?),
                tag => {
                    return Err(reader.malformed(format!("unknown correlation presence tag {tag}")));
                }
            };
            joint_statistics.push(JointStatistics::new(
                observation,
                block_id,
                members,
                covariance,
                correlation,
            ));
        }

        let normalization_count = reader.count("normalizations", MAX_NORMALIZATIONS)?;
        reader.require_items(
            normalization_count,
            MIN_NORMALIZATION_BYTES,
            "normalization receipts",
        )?;
        let mut normalizations = Vec::with_capacity(normalization_count);
        for _ in 0..normalization_count {
            normalizations.push(NormalizationReceipt::new(
                decode_normalization_target(&mut reader)?,
                reader.hash()?,
                reader.dims()?,
                reader.f64()?,
                reader.f64()?,
                reader.string()?,
                reader.string()?,
                reader.optional_string()?,
                reader.optional_string()?,
            ));
        }
        reader.finish()?;
        let pack = Self::new(
            pack_id,
            compiler,
            source_artifact,
            redistribution_terms,
            claims,
            joint_statistics,
            normalizations,
        )?;
        if pack.to_bytes() != bytes {
            return Err(PackError::Malformed {
                at: bytes.len(),
                detail: "decoded fields do not reproduce the canonical byte stream".to_string(),
            });
        }
        Ok(pack)
    }
}

fn invalid(field: &'static str, detail: impl Into<String>) -> PackError {
    PackError::InvalidField {
        field,
        detail: detail.into(),
    }
}

fn limit(resource: &'static str, limit: usize, observed: usize) -> PackError {
    PackError::ResourceLimit {
        resource,
        limit,
        observed,
    }
}

fn require_text(field: &'static str, value: &str) -> Result<(), PackError> {
    if value.trim().is_empty() {
        return Err(invalid(field, "must not be blank"));
    }
    if value.len() > MAX_STRING_BYTES {
        return Err(limit(field, MAX_STRING_BYTES, value.len()));
    }
    Ok(())
}

fn bounded_text(field: &'static str, value: &str) -> Result<(), PackError> {
    if value.len() > MAX_STRING_BYTES {
        return Err(limit(field, MAX_STRING_BYTES, value.len()));
    }
    Ok(())
}

fn portable(value: f64, field: &'static str) -> Result<(), PackError> {
    if !value.is_finite() {
        return Err(invalid(
            field,
            format!("must be finite (bits {:#018x})", value.to_bits()),
        ));
    }
    if value.to_bits() == (-0.0f64).to_bits() {
        return Err(invalid(
            field,
            "negative zero is non-canonical; normalize it to positive zero",
        ));
    }
    Ok(())
}

fn canonical_claim_set(source: &ClaimSet) -> Result<ClaimSet, PackError> {
    let mut canonical = ClaimSet::new();
    for expected in source.observation_ids() {
        let dataset = source
            .observation(expected)
            .ok_or_else(|| invalid("observations", "observation index/lookup mismatch"))?
            .clone();
        let actual = canonical.register_observation(dataset)?;
        if actual != expected {
            return Err(PackError::IdentityMismatch {
                kind: "observation",
                expected: expected.0,
                actual: actual.0,
            });
        }
    }
    for (expected, claim) in source.claims_ordered() {
        let actual = canonical.insert_claim(claim.clone())?;
        if actual != expected {
            return Err(PackError::IdentityMismatch {
                kind: "claim",
                expected: expected.0,
                actual: actual.0,
            });
        }
    }
    Ok(canonical)
}

#[allow(clippy::too_many_lines)] // one admission checklist, checked in its stated order
fn validate_claim_set(claims: &ClaimSet) -> Result<(), PackError> {
    let observation_ids: Vec<_> = claims
        .observation_ids()
        .take(MAX_OBSERVATIONS + 1)
        .collect();
    if observation_ids.is_empty() {
        return Err(invalid(
            "observations",
            "a compiled pack must retain at least one raw observation artifact",
        ));
    }
    if observation_ids.len() > MAX_OBSERVATIONS {
        return Err(limit(
            "observations",
            MAX_OBSERVATIONS,
            observation_ids.len(),
        ));
    }
    let ordered_claims: Vec<_> = claims.claims_ordered().take(MAX_CLAIMS + 1).collect();
    if ordered_claims.is_empty() {
        return Err(invalid(
            "claims",
            "a compiled pack must contain at least one admitted property claim",
        ));
    }
    if ordered_claims.len() > MAX_CLAIMS {
        return Err(limit("claims", MAX_CLAIMS, ordered_claims.len()));
    }
    for (_, claim) in ordered_claims {
        require_text("claim.key.name", claim.key.name())?;
        match &claim.value {
            PropertyValue::Scalar { value, .. } => portable(*value, "claim.scalar")?,
            PropertyValue::Curve {
                abscissa, knots, ..
            } => {
                require_text("claim.curve.abscissa", abscissa)?;
                if knots.len() > MAX_CURVE_KNOTS {
                    return Err(limit("curve_knots", MAX_CURVE_KNOTS, knots.len()));
                }
                for &(x, y) in knots {
                    portable(x, "claim.curve.abscissa_value")?;
                    portable(y, "claim.curve.ordinate_value")?;
                }
            }
        }
        if claim.validity.bounds().len() > MAX_VALIDITY_AXES {
            return Err(limit(
                "validity_axes",
                MAX_VALIDITY_AXES,
                claim.validity.bounds().len(),
            ));
        }
        for (axis, &(lo, hi)) in claim.validity.bounds() {
            require_text("claim.validity.axis", axis)?;
            portable(lo, "claim.validity.lo")?;
            portable(hi, "claim.validity.hi")?;
            if lo > hi {
                return Err(invalid(
                    "claim.validity",
                    format!("axis {axis:?} has empty bounds [{lo}, {hi}]"),
                ));
            }
        }
        match claim.uncertainty {
            UncertaintyModel::Unstated => {}
            UncertaintyModel::HalfWidth {
                half_width,
                confidence,
            } => {
                portable(half_width, "claim.uncertainty.half_width")?;
                portable(confidence, "claim.uncertainty.confidence")?;
            }
            UncertaintyModel::RelativeHalfWidth {
                fraction,
                confidence,
            } => {
                portable(fraction, "claim.uncertainty.fraction")?;
                portable(confidence, "claim.uncertainty.confidence")?;
            }
        }
        if claim.observations.len() > MAX_OBSERVATIONS_PER_CLAIM {
            return Err(limit(
                "claim_observations",
                MAX_OBSERVATIONS_PER_CLAIM,
                claim.observations.len(),
            ));
        }
        if !claim.observations.windows(2).all(|pair| pair[0] < pair[1]) {
            return Err(invalid(
                "claim.observations",
                "observation ids must be strictly increasing and deduplicated",
            ));
        }
        require_text("claim.provenance.source", &claim.provenance.source)?;
        require_text("claim.provenance.license", &claim.provenance.license)?;
    }
    for observation in observation_ids {
        let dataset = claims
            .observation(observation)
            .ok_or_else(|| invalid("observations", "observation index/lookup mismatch"))?;
        require_text("observation.specimen", &dataset.specimen)?;
        require_text("observation.method", &dataset.method)?;
        bounded_text("observation.caveats", &dataset.caveats)?;
        require_text("observation.provenance.source", &dataset.provenance.source)?;
        require_text(
            "observation.provenance.license",
            &dataset.provenance.license,
        )?;
    }
    Ok(())
}

fn validate_joint_statistics(claims: &ClaimSet, block: &JointStatistics) -> Result<(), PackError> {
    require_text("joint_statistics.block_id", &block.block_id)?;
    if claims.observation(block.observation).is_none() {
        return Err(invalid(
            "joint_statistics.observation",
            format!("unknown observation {}", block.observation.0),
        ));
    }
    if block.members.is_empty() {
        return Err(invalid(
            "joint_statistics.members",
            "a covariance block needs at least one member",
        ));
    }
    if block.members.len() > MAX_JOINT_MEMBERS {
        return Err(limit(
            "joint_members",
            MAX_JOINT_MEMBERS,
            block.members.len(),
        ));
    }
    if !block.members.windows(2).all(|pair| pair[0] < pair[1]) {
        return Err(invalid(
            "joint_statistics.members",
            "statistic members must be strictly increasing and deduplicated",
        ));
    }
    for member in &block.members {
        let claim = claims.claim(member.claim).ok_or_else(|| {
            invalid(
                "joint_statistics.members",
                format!("unknown claim {}", member.claim.0),
            )
        })?;
        statistic_member_dims(claims, *member)?;
        if claim
            .observations
            .binary_search(&block.observation)
            .is_err()
        {
            return Err(invalid(
                "joint_statistics.members",
                format!(
                    "claim {} does not cite joint observation {}",
                    member.claim.0, block.observation.0
                ),
            ));
        }
    }
    for (row_index, row) in block.members.iter().enumerate() {
        let row_dims = statistic_member_dims(claims, *row)?;
        for column in &block.members[..=row_index] {
            let column_dims = statistic_member_dims(claims, *column)?;
            if row_dims.checked_plus(column_dims).is_none() {
                return Err(invalid(
                    "joint_statistics.covariance",
                    format!(
                        "covariance dimensions overflow for claims {} and {}",
                        row.claim.0, column.claim.0
                    ),
                ));
            }
        }
    }
    let expected = triangle_len(block.members.len())?;
    if block.covariance.len() != expected {
        return Err(invalid(
            "joint_statistics.covariance",
            format!(
                "packed lower triangle needs {expected} entries for {} members, found {}",
                block.members.len(),
                block.covariance.len()
            ),
        ));
    }
    let implied_correlation = validate_covariance(&block.covariance, block.members.len())?;
    if let Some(correlation) = &block.correlation {
        if correlation.len() != expected {
            return Err(invalid(
                "joint_statistics.correlation",
                format!(
                    "packed lower triangle needs {expected} entries for {} members, found {}",
                    block.members.len(),
                    correlation.len()
                ),
            ));
        }
        if (0..block.members.len()).any(|row| block.covariance[packed_index(row, row)] == 0.0) {
            return Err(invalid(
                "joint_statistics.correlation",
                "correlation is undefined when a covariance member has zero variance",
            ));
        }
        validate_correlation(correlation, block.members.len())?;
        validate_correlation_consistency(&implied_correlation, correlation)?;
    }
    Ok(())
}

fn statistic_member_dims(claims: &ClaimSet, member: StatisticMember) -> Result<Dims, PackError> {
    statistic_member_dims_at(claims, member, "joint_statistics.members")
}

fn statistic_member_dims_at(
    claims: &ClaimSet,
    member: StatisticMember,
    field: &'static str,
) -> Result<Dims, PackError> {
    let claim = claims
        .claim(member.claim)
        .ok_or_else(|| invalid(field, format!("unknown claim {}", member.claim.0)))?;
    match (&claim.value, member.component) {
        (PropertyValue::Scalar { dims, .. }, StatisticComponent::Scalar) => Ok(*dims),
        (
            PropertyValue::Curve {
                abscissa_dims,
                knots,
                ..
            },
            StatisticComponent::CurveAbscissa { knot },
        ) => {
            validate_knot_index(knot, knots.len(), field)?;
            Ok(*abscissa_dims)
        }
        (PropertyValue::Curve { knots, dims, .. }, StatisticComponent::CurveOrdinate { knot }) => {
            validate_knot_index(knot, knots.len(), field)?;
            Ok(*dims)
        }
        (PropertyValue::Scalar { .. }, component) => Err(invalid(
            field,
            format!("scalar claim cannot supply component {component:?}"),
        )),
        (PropertyValue::Curve { .. }, StatisticComponent::Scalar) => Err(invalid(
            field,
            "curve claim requires an explicit abscissa/ordinate knot component",
        )),
    }
}

fn validate_knot_index(knot: u32, knot_count: usize, field: &'static str) -> Result<(), PackError> {
    let index =
        usize::try_from(knot).map_err(|_| invalid(field, "knot index does not fit usize"))?;
    if index < knot_count {
        Ok(())
    } else {
        Err(invalid(
            field,
            format!("curve knot {knot} is out of range for {knot_count} knots"),
        ))
    }
}

fn triangle_len(size: usize) -> Result<usize, PackError> {
    size.checked_add(1)
        .and_then(|next| size.checked_mul(next))
        .map(|product| product / 2)
        .ok_or_else(|| invalid("joint_statistics.members", "triangle size overflow"))
}

fn packed_index(row: usize, column: usize) -> usize {
    row * (row + 1) / 2 + column
}

fn validate_psd_work_budget(blocks: &[JointStatistics]) -> Result<(), PackError> {
    let mut work = 0u64;
    let public_limit = usize::try_from(MAX_PSD_CUBIC_WORK).unwrap_or(usize::MAX);
    for block in blocks {
        let size = u64::try_from(block.members.len()).map_err(|_| {
            invalid(
                "joint_statistics.members",
                "member count does not fit the PSD work counter",
            )
        })?;
        let matrices = if block.correlation.is_some() { 2 } else { 1 };
        let block_work = size
            .checked_mul(size)
            .and_then(|value| value.checked_mul(size))
            .and_then(|value| value.checked_mul(matrices))
            .ok_or_else(|| limit("psd_cubic_work", public_limit, usize::MAX))?;
        work = work
            .checked_add(block_work)
            .ok_or_else(|| limit("psd_cubic_work", public_limit, usize::MAX))?;
        if work > MAX_PSD_CUBIC_WORK {
            return Err(PackError::ResourceLimit {
                resource: "psd_cubic_work",
                limit: public_limit,
                observed: usize::try_from(work).unwrap_or(usize::MAX),
            });
        }
    }
    Ok(())
}

fn validate_disjoint_joint_statistics(blocks: &[JointStatistics]) -> Result<(), PackError> {
    let mut seen = BTreeSet::new();
    for block in blocks {
        for &member in &block.members {
            if !seen.insert((block.observation, member)) {
                return Err(invalid(
                    "joint_statistics.members",
                    format!(
                        "member {member:?} occurs in more than one block for observation {}",
                        block.observation.0
                    ),
                ));
            }
        }
    }
    Ok(())
}

fn validate_covariance(packed: &[f64], size: usize) -> Result<Vec<f64>, PackError> {
    for &value in packed {
        portable(value, "covariance")?;
    }
    let mut variances = Vec::with_capacity(size);
    for row in 0..size {
        let diagonal = packed[packed_index(row, row)];
        if diagonal < 0.0 {
            return Err(invalid(
                "joint_statistics.covariance",
                format!("diagonal {row} is negative ({diagonal})"),
            ));
        }
        variances.push(diagonal);
    }

    // Derive the exact canonical correlation representation used for
    // pairwise/source agreement. The raw covariance itself is separately
    // gated by outward-rounded interval LDLT below.
    let mut equilibrated = vec![0.0; packed.len()];
    for row in 0..size {
        equilibrated[packed_index(row, row)] = if variances[row] == 0.0 { 0.0 } else { 1.0 };
        for column in 0..row {
            let value = packed[packed_index(row, column)];
            let row_variance = variances[row];
            let column_variance = variances[column];
            if row_variance == 0.0 || column_variance == 0.0 {
                if value != 0.0 {
                    return Err(invalid(
                        "joint_statistics.covariance",
                        format!(
                            "zero variance requires exact zero covariance at ({row},{column}), found {value}"
                        ),
                    ));
                }
                continue;
            }
            let denominator = row_variance.sqrt() * column_variance.sqrt();
            if !denominator.is_finite() || denominator == 0.0 {
                return Err(invalid(
                    "joint_statistics.covariance",
                    format!("cannot equilibrate entry ({row},{column})"),
                ));
            }
            let correlation = value / denominator;
            portable(correlation, "covariance.equilibrated")?;
            if correlation.abs() > 1.0 {
                return Err(invalid(
                    "joint_statistics.covariance",
                    format!(
                        "entry ({row},{column}) implies correlation outside [-1,1]: {correlation}"
                    ),
                ));
            }
            equilibrated[packed_index(row, column)] = correlation;
        }
    }
    validate_psd(packed, size, "covariance")?;
    Ok(equilibrated)
}

fn validate_correlation(packed: &[f64], size: usize) -> Result<(), PackError> {
    for &value in packed {
        portable(value, "correlation")?;
    }
    for row in 0..size {
        let diagonal = packed[packed_index(row, row)];
        if diagonal.to_bits() != 1.0f64.to_bits() {
            return Err(invalid(
                "joint_statistics.correlation",
                format!("diagonal {row} must be exactly 1.0, found {diagonal}"),
            ));
        }
        for column in 0..row {
            let value = packed[packed_index(row, column)];
            if value.abs() > 1.0 {
                return Err(invalid(
                    "joint_statistics.correlation",
                    format!("entry ({row},{column}) lies outside [-1,1]: {value}"),
                ));
            }
        }
    }
    validate_psd(packed, size, "correlation")
}

fn validate_correlation_consistency(implied: &[f64], supplied: &[f64]) -> Result<(), PackError> {
    for (index, (&expected, &actual)) in implied.iter().zip(supplied).enumerate() {
        if expected.to_bits() != actual.to_bits() {
            return Err(invalid(
                "joint_statistics.correlation",
                format!(
                    "entry {index} is inconsistent with covariance: expected {expected}, found {actual}"
                ),
            ));
        }
    }
    Ok(())
}

#[derive(Clone, Copy)]
struct Enclosure {
    lo: f64,
    hi: f64,
}

impl Enclosure {
    const ZERO: Self = Self { lo: 0.0, hi: 0.0 };

    fn point(value: f64) -> Self {
        Self {
            lo: value,
            hi: value,
        }
    }

    fn is_exact(self, value: f64) -> bool {
        self.lo.to_bits() == value.to_bits() && self.hi.to_bits() == value.to_bits()
    }

    fn negate(self) -> Self {
        Self {
            lo: -self.hi,
            hi: -self.lo,
        }
    }

    fn outward(
        raw_lo: f64,
        raw_hi: f64,
        field: &'static str,
        operation: &str,
    ) -> Result<Self, PackError> {
        if !raw_lo.is_finite() || !raw_hi.is_finite() {
            return Err(invalid(
                field,
                format!("PSD interval {operation} overflowed"),
            ));
        }
        let lo = next_down(raw_lo);
        let hi = next_up(raw_hi);
        if !lo.is_finite() || !hi.is_finite() {
            return Err(invalid(
                field,
                format!("PSD interval {operation} exceeded finite range"),
            ));
        }
        Ok(Self { lo, hi })
    }

    fn subtract(self, rhs: Self, field: &'static str) -> Result<Self, PackError> {
        if rhs.is_exact(0.0) {
            return Ok(self);
        }
        Self::outward(self.lo - rhs.hi, self.hi - rhs.lo, field, "subtraction")
    }

    fn multiply(self, rhs: Self, field: &'static str) -> Result<Self, PackError> {
        if self.is_exact(0.0) || rhs.is_exact(0.0) {
            return Ok(Self::ZERO);
        }
        if self.is_exact(1.0) {
            return Ok(rhs);
        }
        if rhs.is_exact(1.0) {
            return Ok(self);
        }
        if self.is_exact(-1.0) {
            return Ok(rhs.negate());
        }
        if rhs.is_exact(-1.0) {
            return Ok(self.negate());
        }
        let products = [
            self.lo * rhs.lo,
            self.lo * rhs.hi,
            self.hi * rhs.lo,
            self.hi * rhs.hi,
        ];
        if products.iter().any(|value| !value.is_finite()) {
            return Err(invalid(field, "PSD interval multiplication overflowed"));
        }
        let raw_lo = products.iter().copied().fold(f64::INFINITY, f64::min);
        let raw_hi = products.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        Self::outward(raw_lo, raw_hi, field, "multiplication")
    }

    fn square(self, field: &'static str) -> Result<Self, PackError> {
        if self.is_exact(0.0) {
            return Ok(Self::ZERO);
        }
        if self.is_exact(1.0) || self.is_exact(-1.0) {
            return Ok(Self::point(1.0));
        }
        let (raw_lo, raw_hi) = if self.lo >= 0.0 {
            (self.lo * self.lo, self.hi * self.hi)
        } else if self.hi <= 0.0 {
            (self.hi * self.hi, self.lo * self.lo)
        } else {
            (0.0, (self.lo * self.lo).max(self.hi * self.hi))
        };
        Self::outward(raw_lo, raw_hi, field, "squaring")
    }

    fn divide(self, rhs: Self, field: &'static str) -> Result<Self, PackError> {
        if rhs.lo <= 0.0 {
            return Err(invalid(
                field,
                "PSD interval division lacks a positive denominator proof",
            ));
        }
        if self.is_exact(0.0) {
            return Ok(Self::ZERO);
        }
        if rhs.is_exact(1.0) {
            return Ok(self);
        }
        let quotients = [
            self.lo / rhs.lo,
            self.lo / rhs.hi,
            self.hi / rhs.lo,
            self.hi / rhs.hi,
        ];
        if quotients.iter().any(|value| !value.is_finite()) {
            return Err(invalid(field, "PSD interval division overflowed"));
        }
        let raw_lo = quotients.iter().copied().fold(f64::INFINITY, f64::min);
        let raw_hi = quotients.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        Self::outward(raw_lo, raw_hi, field, "division")
    }
}

fn next_up(value: f64) -> f64 {
    if value.is_nan() || value == f64::INFINITY {
        return value;
    }
    if value == 0.0 {
        return f64::from_bits(1);
    }
    let bits = value.to_bits();
    if value > 0.0 {
        f64::from_bits(bits + 1)
    } else {
        f64::from_bits(bits - 1)
    }
}

fn next_down(value: f64) -> f64 {
    if value.is_nan() || value == f64::NEG_INFINITY {
        return value;
    }
    if value == 0.0 {
        return f64::from_bits((1u64 << 63) | 1);
    }
    let bits = value.to_bits();
    if value > 0.0 {
        f64::from_bits(bits - 1)
    } else {
        f64::from_bits(bits + 1)
    }
}

fn validate_psd(packed: &[f64], size: usize, field: &'static str) -> Result<(), PackError> {
    for &value in packed {
        portable(value, field)?;
    }

    // Every arithmetic operation carries an outward-rounded interval. A
    // positive pivot is admitted only when its lower bound is positive; an
    // exact structural zero is admitted only with exact-zero residuals. This
    // deliberately refuses ill-conditioned valid matrices rather than
    // certifying a rounded zero whose exact pivot may be negative.
    let mut lower = vec![Enclosure::ZERO; packed.len()];
    let mut pivots = vec![Enclosure::ZERO; size];
    for column in 0..size {
        let mut pivot = Enclosure::point(packed[packed_index(column, column)]);
        for prior in 0..column {
            let value = lower[packed_index(column, prior)];
            let contribution = value.square(field)?.multiply(pivots[prior], field)?;
            pivot = pivot.subtract(contribution, field)?;
        }
        if pivot.hi < 0.0 {
            return Err(invalid(
                field,
                format!(
                    "matrix is not positive semidefinite at pivot {column}: upper bound {}",
                    pivot.hi
                ),
            ));
        }
        let zero_pivot = pivot.is_exact(0.0);
        if !zero_pivot && pivot.lo <= 0.0 {
            return Err(invalid(
                field,
                format!(
                    "rounding-ambiguous PSD pivot {column} is enclosed by [{}, {}]",
                    pivot.lo, pivot.hi
                ),
            ));
        }
        if zero_pivot {
            pivots[column] = Enclosure::ZERO;
            for row in column + 1..size {
                let mut residual = Enclosure::point(packed[packed_index(row, column)]);
                for prior in 0..column {
                    let contribution = lower[packed_index(row, prior)]
                        .multiply(lower[packed_index(column, prior)], field)?
                        .multiply(pivots[prior], field)?;
                    residual = residual.subtract(contribution, field)?;
                }
                if !residual.is_exact(0.0) {
                    return Err(invalid(
                        field,
                        format!(
                            "zero pivot {column} lacks an exact-zero residual proof at row {row}: [{}, {}]",
                            residual.lo, residual.hi
                        ),
                    ));
                }
            }
        } else {
            pivots[column] = pivot;
            for row in column + 1..size {
                let mut residual = Enclosure::point(packed[packed_index(row, column)]);
                for prior in 0..column {
                    let contribution = lower[packed_index(row, prior)]
                        .multiply(lower[packed_index(column, prior)], field)?
                        .multiply(pivots[prior], field)?;
                    residual = residual.subtract(contribution, field)?;
                }
                lower[packed_index(row, column)] = residual.divide(pivot, field)?;
            }
        }
    }
    Ok(())
}

fn validate_normalization(
    claims: &ClaimSet,
    blocks: &[JointStatistics],
    receipt: &NormalizationReceipt,
) -> Result<(), PackError> {
    if let Some(expected_dims) = normalization_target_dims(claims, blocks, &receipt.target)?
        && expected_dims != receipt.dims
    {
        return Err(invalid(
            "normalization.dims",
            format!(
                "target {:?} requires {expected_dims:?}, found {:?}",
                receipt.target, receipt.dims
            ),
        ));
    }
    require_text("normalization.source_basis", &receipt.source_basis)?;
    require_text("normalization.target_basis", &receipt.target_basis)?;
    if receipt.target_basis != MATDB_PACK_TARGET_BASIS {
        return Err(invalid(
            "normalization.target_basis",
            format!(
                "must be the canonical target {MATDB_PACK_TARGET_BASIS:?}, found {:?}",
                receipt.target_basis
            ),
        ));
    }
    portable(receipt.scale, "normalization.scale")?;
    portable(receipt.offset, "normalization.offset")?;
    if receipt.scale == 0.0 {
        return Err(invalid(
            "normalization.scale",
            "a zero scale is not an invertible unit/basis transform",
        ));
    }
    if matches!(
        &receipt.target,
        NormalizationTarget::ClaimUncertainty { .. } | NormalizationTarget::JointCovariance { .. }
    ) && receipt.offset.to_bits() != 0.0f64.to_bits()
    {
        return Err(invalid(
            "normalization.offset",
            "uncertainty and covariance transforms must not translate",
        ));
    }
    let positive_scale_required = match &receipt.target {
        NormalizationTarget::ClaimUncertainty { .. } => true,
        NormalizationTarget::JointCovariance { row, column, .. } => row == column,
        _ => false,
    };
    if positive_scale_required && receipt.scale < 0.0 {
        return Err(invalid(
            "normalization.scale",
            "uncertainty magnitudes and variances require a positive scale",
        ));
    }
    match (&receipt.source_frame, &receipt.target_frame) {
        (None, None) => {}
        (Some(source), Some(target)) => {
            require_text("normalization.source_frame", source)?;
            require_text("normalization.target_frame", target)?;
        }
        _ => {
            return Err(invalid(
                "normalization.frames",
                "source and target frame must either both be present or both be absent",
            ));
        }
    }
    Ok(())
}

fn validate_normalization_coherence(receipts: &[NormalizationReceipt]) -> Result<(), PackError> {
    let mut validity_dims = BTreeMap::new();
    for receipt in receipts {
        if let NormalizationTarget::ValidityBound { claim, axis, .. } = &receipt.target {
            let key = (*claim, axis.as_str());
            if let Some(previous) = validity_dims.insert(key, receipt.dims)
                && previous != receipt.dims
            {
                return Err(invalid(
                    "normalization.dims",
                    format!(
                        "validity axis {axis:?} on claim {} has contradictory dimensions {previous:?} and {:?}",
                        claim.0, receipt.dims
                    ),
                ));
            }
        }
    }
    Ok(())
}

fn normalization_target_dims(
    claims: &ClaimSet,
    blocks: &[JointStatistics],
    target: &NormalizationTarget,
) -> Result<Option<Dims>, PackError> {
    match target {
        NormalizationTarget::ClaimValue(member) => {
            statistic_member_dims_at(claims, *member, "normalization.target").map(Some)
        }
        NormalizationTarget::ClaimUncertainty { claim } => {
            let claim = claims.claim(*claim).ok_or_else(|| {
                invalid(
                    "normalization.target",
                    format!("unknown uncertainty claim {}", claim.0),
                )
            })?;
            match claim.uncertainty {
                UncertaintyModel::Unstated => Err(invalid(
                    "normalization.target",
                    "Unstated uncertainty has no numeric field to normalize",
                )),
                UncertaintyModel::HalfWidth { .. } => Ok(Some(claim.value.dims())),
                UncertaintyModel::RelativeHalfWidth { .. } => Ok(Some(Dims::NONE)),
            }
        }
        NormalizationTarget::ValidityBound { claim, axis, .. } => {
            require_text("normalization.target.validity_axis", axis)?;
            let claim = claims.claim(*claim).ok_or_else(|| {
                invalid(
                    "normalization.target",
                    format!("unknown validity claim {}", claim.0),
                )
            })?;
            if claim.validity.bound(axis).is_none() {
                return Err(invalid(
                    "normalization.target",
                    format!(
                        "claim {} has no validity axis {axis:?}",
                        claim.content_hash()
                    ),
                ));
            }
            // ValidityDomain currently stores values but not axis dimensions;
            // the compiler receipt is the sole dimension authority here.
            Ok(None)
        }
        NormalizationTarget::JointCovariance {
            observation,
            block_id,
            row,
            column,
        } => {
            require_text("normalization.target.block_id", block_id)?;
            let block = blocks
                .binary_search_by(|block| {
                    (block.observation, block.block_id.as_str())
                        .cmp(&(*observation, block_id.as_str()))
                })
                .map(|index| &blocks[index])
                .map_err(|_| {
                    invalid(
                        "normalization.target",
                        format!(
                            "unknown covariance block {block_id:?} for observation {}",
                            observation.0
                        ),
                    )
                })?;
            let row = usize::try_from(*row)
                .map_err(|_| invalid("normalization.target", "row does not fit usize"))?;
            let column = usize::try_from(*column)
                .map_err(|_| invalid("normalization.target", "column does not fit usize"))?;
            if column > row || row >= block.members.len() {
                return Err(invalid(
                    "normalization.target",
                    format!(
                        "covariance coordinate ({row},{column}) is outside lower triangle of {} members",
                        block.members.len()
                    ),
                ));
            }
            let row_dims =
                statistic_member_dims_at(claims, block.members[row], "normalization.target")?;
            let column_dims =
                statistic_member_dims_at(claims, block.members[column], "normalization.target")?;
            row_dims.checked_plus(column_dims).map(Some).ok_or_else(|| {
                invalid(
                    "normalization.target",
                    "covariance target dimensions overflow six-base exponents",
                )
            })
        }
    }
}

#[derive(Default)]
struct Writer {
    bytes: Vec<u8>,
}

impl Writer {
    fn u8(&mut self, value: u8) {
        self.bytes.push(value);
    }

    fn u32(&mut self, value: u32) {
        self.bytes.extend_from_slice(&value.to_le_bytes());
    }

    fn u64(&mut self, value: u64) {
        self.bytes.extend_from_slice(&value.to_le_bytes());
    }

    fn count(&mut self, value: usize) {
        self.u32(u32::try_from(value).unwrap_or(u32::MAX));
    }

    fn string(&mut self, value: &str) {
        self.count(value.len());
        self.bytes.extend_from_slice(value.as_bytes());
    }

    fn optional_string(&mut self, value: Option<&str>) {
        match value {
            None => self.u8(0),
            Some(value) => {
                self.u8(1);
                self.string(value);
            }
        }
    }

    fn hash(&mut self, value: ContentHash) {
        self.bytes.extend_from_slice(&value.0);
    }

    fn dims(&mut self, value: Dims) {
        for exponent in value.0 {
            self.u8(exponent.cast_unsigned());
        }
    }

    fn f64(&mut self, value: f64) {
        self.u64(value.to_bits());
    }

    fn f64s(&mut self, values: &[f64]) {
        for &value in values {
            self.f64(value);
        }
    }
}

struct Reader<'a> {
    bytes: &'a [u8],
    cursor: usize,
}

impl<'a> Reader<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, cursor: 0 }
    }

    fn malformed(&self, detail: impl Into<String>) -> PackError {
        PackError::Malformed {
            at: self.cursor,
            detail: detail.into(),
        }
    }

    fn take(&mut self, length: usize) -> Result<&'a [u8], PackError> {
        let end = self
            .cursor
            .checked_add(length)
            .ok_or_else(|| self.malformed("byte offset overflow"))?;
        let slice = self
            .bytes
            .get(self.cursor..end)
            .ok_or_else(|| self.malformed(format!("truncated field needs {length} bytes")))?;
        self.cursor = end;
        Ok(slice)
    }

    fn require_remaining(&self, length: usize, field: &str) -> Result<(), PackError> {
        if self.bytes.len().saturating_sub(self.cursor) < length {
            Err(self.malformed(format!("truncated {field} needs {length} bytes")))
        } else {
            Ok(())
        }
    }

    fn require_items(
        &self,
        count: usize,
        minimum_width: usize,
        field: &str,
    ) -> Result<(), PackError> {
        let minimum_bytes = count
            .checked_mul(minimum_width)
            .ok_or_else(|| self.malformed(format!("{field} byte length overflow")))?;
        self.require_remaining(minimum_bytes, field)
    }

    fn expect(&mut self, expected: &[u8], name: &str) -> Result<(), PackError> {
        let actual = self.take(expected.len())?;
        if actual == expected {
            Ok(())
        } else {
            Err(self.malformed(format!("invalid {name}")))
        }
    }

    fn u8(&mut self) -> Result<u8, PackError> {
        Ok(self.take(1)?[0])
    }

    fn u32(&mut self) -> Result<u32, PackError> {
        let bytes: [u8; 4] = self
            .take(4)?
            .try_into()
            .map_err(|_| self.malformed("u32 width"))?;
        Ok(u32::from_le_bytes(bytes))
    }

    fn u64(&mut self) -> Result<u64, PackError> {
        let bytes: [u8; 8] = self
            .take(8)?
            .try_into()
            .map_err(|_| self.malformed("u64 width"))?;
        Ok(u64::from_le_bytes(bytes))
    }

    fn count(&mut self, resource: &'static str, maximum: usize) -> Result<usize, PackError> {
        let raw = self.u32()?;
        let count = usize::try_from(raw)
            .map_err(|_| self.malformed(format!("{resource} count does not fit usize")))?;
        if count > maximum {
            Err(limit(resource, maximum, count))
        } else {
            Ok(count)
        }
    }

    fn string(&mut self) -> Result<String, PackError> {
        let length = self.count("string_bytes", MAX_STRING_BYTES)?;
        let start = self.cursor;
        let bytes = self.take(length)?;
        std::str::from_utf8(bytes)
            .map(str::to_string)
            .map_err(|error| PackError::Malformed {
                at: start + error.valid_up_to(),
                detail: "string field is not UTF-8".to_string(),
            })
    }

    fn optional_string(&mut self) -> Result<Option<String>, PackError> {
        match self.u8()? {
            0 => Ok(None),
            1 => self.string().map(Some),
            tag => Err(self.malformed(format!("unknown optional-string tag {tag}"))),
        }
    }

    fn hash(&mut self) -> Result<ContentHash, PackError> {
        let bytes: [u8; 32] = self
            .take(32)?
            .try_into()
            .map_err(|_| self.malformed("content-hash width"))?;
        Ok(ContentHash(bytes))
    }

    fn dims(&mut self) -> Result<Dims, PackError> {
        let mut dims = [0i8; 6];
        for exponent in &mut dims {
            *exponent = i8::from_ne_bytes([self.u8()?]);
        }
        Ok(Dims(dims))
    }

    fn f64(&mut self) -> Result<f64, PackError> {
        Ok(f64::from_bits(self.u64()?))
    }

    fn fixed_f64s(&mut self, count: usize) -> Result<Vec<f64>, PackError> {
        let byte_count = count
            .checked_mul(8)
            .ok_or_else(|| self.malformed("float vector byte length overflow"))?;
        if self.bytes.len().saturating_sub(self.cursor) < byte_count {
            return Err(self.malformed(format!("truncated float vector needs {byte_count} bytes")));
        }
        let mut values = Vec::with_capacity(count);
        for _ in 0..count {
            values.push(self.f64()?);
        }
        Ok(values)
    }

    fn finish(self) -> Result<(), PackError> {
        if self.cursor == self.bytes.len() {
            Ok(())
        } else {
            Err(PackError::Malformed {
                at: self.cursor,
                detail: format!(
                    "{} trailing bytes after canonical pack",
                    self.bytes.len() - self.cursor
                ),
            })
        }
    }
}

fn encode_statistic_member(writer: &mut Writer, member: StatisticMember) {
    writer.hash(member.claim.0);
    match member.component {
        StatisticComponent::Scalar => writer.u8(0),
        StatisticComponent::CurveAbscissa { knot } => {
            writer.u8(1);
            writer.u32(knot);
        }
        StatisticComponent::CurveOrdinate { knot } => {
            writer.u8(2);
            writer.u32(knot);
        }
    }
}

fn decode_statistic_member(reader: &mut Reader<'_>) -> Result<StatisticMember, PackError> {
    let claim = ClaimId(reader.hash()?);
    match reader.u8()? {
        0 => Ok(StatisticMember::scalar(claim)),
        1 => Ok(StatisticMember::curve_abscissa(claim, reader.u32()?)),
        2 => Ok(StatisticMember::curve_ordinate(claim, reader.u32()?)),
        tag => Err(reader.malformed(format!("unknown statistic-component tag {tag}"))),
    }
}

fn encode_normalization_target(writer: &mut Writer, target: &NormalizationTarget) {
    match target {
        NormalizationTarget::ClaimValue(member) => {
            writer.u8(0);
            encode_statistic_member(writer, *member);
        }
        NormalizationTarget::ClaimUncertainty { claim } => {
            writer.u8(1);
            writer.hash(claim.0);
        }
        NormalizationTarget::ValidityBound { claim, axis, side } => {
            writer.u8(2);
            writer.hash(claim.0);
            writer.string(axis);
            writer.u8(match side {
                ValidityBoundSide::Lower => 0,
                ValidityBoundSide::Upper => 1,
            });
        }
        NormalizationTarget::JointCovariance {
            observation,
            block_id,
            row,
            column,
        } => {
            writer.u8(3);
            writer.hash(observation.0);
            writer.string(block_id);
            writer.u32(*row);
            writer.u32(*column);
        }
    }
}

fn decode_normalization_target(reader: &mut Reader<'_>) -> Result<NormalizationTarget, PackError> {
    match reader.u8()? {
        0 => Ok(NormalizationTarget::ClaimValue(decode_statistic_member(
            reader,
        )?)),
        1 => Ok(NormalizationTarget::ClaimUncertainty {
            claim: ClaimId(reader.hash()?),
        }),
        2 => {
            let claim = ClaimId(reader.hash()?);
            let axis = reader.string()?;
            let side = match reader.u8()? {
                0 => ValidityBoundSide::Lower,
                1 => ValidityBoundSide::Upper,
                tag => {
                    return Err(reader.malformed(format!("unknown validity-bound-side tag {tag}")));
                }
            };
            Ok(NormalizationTarget::ValidityBound { claim, axis, side })
        }
        3 => Ok(NormalizationTarget::JointCovariance {
            observation: ObservationId(reader.hash()?),
            block_id: reader.string()?,
            row: reader.u32()?,
            column: reader.u32()?,
        }),
        tag => Err(reader.malformed(format!("unknown normalization-target tag {tag}"))),
    }
}

fn encode_provenance(writer: &mut Writer, provenance: &Provenance) {
    writer.string(&provenance.source);
    writer.string(&provenance.license);
    match provenance.artifact {
        None => writer.u8(0),
        Some(hash) => {
            writer.u8(1);
            writer.hash(hash);
        }
    }
}

fn decode_provenance(reader: &mut Reader<'_>) -> Result<Provenance, PackError> {
    let source = reader.string()?;
    let license = reader.string()?;
    let artifact = match reader.u8()? {
        0 => None,
        1 => Some(reader.hash()?),
        tag => return Err(reader.malformed(format!("unknown provenance-artifact tag {tag}"))),
    };
    Ok(Provenance {
        source,
        license,
        artifact,
    })
}

fn encode_observation(writer: &mut Writer, dataset: &ObservationDataset) {
    writer.string(&dataset.specimen);
    writer.string(&dataset.method);
    writer.hash(dataset.artifact);
    writer.string(&dataset.caveats);
    encode_provenance(writer, &dataset.provenance);
}

fn decode_observation(reader: &mut Reader<'_>) -> Result<ObservationDataset, PackError> {
    Ok(ObservationDataset {
        specimen: reader.string()?,
        method: reader.string()?,
        artifact: reader.hash()?,
        caveats: reader.string()?,
        provenance: decode_provenance(reader)?,
    })
}

fn encode_claim(writer: &mut Writer, claim: &PropertyClaim) {
    writer.string(claim.key.name());
    writer.dims(claim.key.dims());
    match &claim.value {
        PropertyValue::Scalar { value, dims } => {
            writer.u8(0);
            writer.f64(*value);
            writer.dims(*dims);
        }
        PropertyValue::Curve {
            abscissa,
            abscissa_dims,
            knots,
            dims,
        } => {
            writer.u8(1);
            writer.string(abscissa);
            writer.dims(*abscissa_dims);
            writer.count(knots.len());
            for &(x, y) in knots {
                writer.f64(x);
                writer.f64(y);
            }
            writer.dims(*dims);
        }
    }
    writer.count(claim.validity.bounds().len());
    for (axis, &(lo, hi)) in claim.validity.bounds() {
        writer.string(axis);
        writer.f64(lo);
        writer.f64(hi);
    }
    match claim.uncertainty {
        UncertaintyModel::Unstated => writer.u8(0),
        UncertaintyModel::HalfWidth {
            half_width,
            confidence,
        } => {
            writer.u8(1);
            writer.f64(half_width);
            writer.f64(confidence);
        }
        UncertaintyModel::RelativeHalfWidth {
            fraction,
            confidence,
        } => {
            writer.u8(2);
            writer.f64(fraction);
            writer.f64(confidence);
        }
    }
    writer.u8(match claim.interpolation {
        InterpolationPolicy::LinearInside => 0,
        InterpolationPolicy::ConstantWithinValidity => 1,
        InterpolationPolicy::TabulatedOnly => 2,
    });
    writer.count(claim.observations.len());
    for observation in &claim.observations {
        writer.hash(observation.0);
    }
    encode_provenance(writer, &claim.provenance);
}

fn decode_claim(reader: &mut Reader<'_>) -> Result<PropertyClaim, PackError> {
    let name = reader.string()?;
    let key_dims = reader.dims()?;
    let value = match reader.u8()? {
        0 => PropertyValue::Scalar {
            value: reader.f64()?,
            dims: reader.dims()?,
        },
        1 => {
            let abscissa = reader.string()?;
            let abscissa_dims = reader.dims()?;
            let knot_count = reader.count("curve_knots", MAX_CURVE_KNOTS)?;
            let knot_bytes = knot_count
                .checked_mul(16)
                .ok_or_else(|| reader.malformed("curve-knot byte length overflow"))?;
            reader.require_remaining(knot_bytes, "curve knots")?;
            let mut knots = Vec::with_capacity(knot_count);
            for _ in 0..knot_count {
                knots.push((reader.f64()?, reader.f64()?));
            }
            PropertyValue::Curve {
                abscissa,
                abscissa_dims,
                knots,
                dims: reader.dims()?,
            }
        }
        tag => return Err(reader.malformed(format!("unknown property-value tag {tag}"))),
    };
    let validity_count = reader.count("validity_axes", MAX_VALIDITY_AXES)?;
    let mut validity = ValidityDomain::unconstrained();
    let mut previous_axis: Option<String> = None;
    for _ in 0..validity_count {
        let axis = reader.string()?;
        if previous_axis
            .as_ref()
            .is_some_and(|previous| previous >= &axis)
        {
            return Err(
                reader.malformed("validity axes are not strictly increasing and deduplicated")
            );
        }
        let lo = reader.f64()?;
        let hi = reader.f64()?;
        previous_axis = Some(axis.clone());
        validity = validity.with(axis, lo, hi);
    }
    let uncertainty = match reader.u8()? {
        0 => UncertaintyModel::Unstated,
        1 => UncertaintyModel::HalfWidth {
            half_width: reader.f64()?,
            confidence: reader.f64()?,
        },
        2 => UncertaintyModel::RelativeHalfWidth {
            fraction: reader.f64()?,
            confidence: reader.f64()?,
        },
        tag => return Err(reader.malformed(format!("unknown uncertainty tag {tag}"))),
    };
    let interpolation = match reader.u8()? {
        0 => InterpolationPolicy::LinearInside,
        1 => InterpolationPolicy::ConstantWithinValidity,
        2 => InterpolationPolicy::TabulatedOnly,
        tag => return Err(reader.malformed(format!("unknown interpolation tag {tag}"))),
    };
    let observation_count = reader.count("claim_observations", MAX_OBSERVATIONS_PER_CLAIM)?;
    reader.require_items(
        observation_count,
        OBSERVATION_ID_BYTES,
        "claim observation ids",
    )?;
    let mut observations = Vec::with_capacity(observation_count);
    for _ in 0..observation_count {
        observations.push(ObservationId(reader.hash()?));
    }
    Ok(PropertyClaim {
        key: PropertyKey::new(name, key_dims),
        value,
        validity,
        uncertainty,
        interpolation,
        observations,
        provenance: decode_provenance(reader)?,
    })
}
