//! Bounded canonical transport for admitted Euler render trajectories.
//!
//! The codec stores the accepted trajectory inputs once, together with the
//! declared timeline-composition seams needed by render consumers. Derived
//! visualization and audio controls are deliberately not duplicated: their
//! schema versions are pinned in the envelope and they are regenerated from
//! the decoded [`RenderTrajectory`].

use core::fmt;
use std::io::{self, Cursor, Read, Seek, SeekFrom, Write};

use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_mbd::{MassProperties, Pose, RigidBodyState, UnitQuaternion, Vec3};

use crate::{
    EULER_CONTROL_STREAM_SCHEMA_VERSION,
    coupled_runner::{ChannelOwnership, ChannelWrench, ContactTransitionKind},
    render_trajectory::{
        DerivedEulerQois, EULER_RENDER_TRAJECTORY_SCHEMA_VERSION,
        MAX_RENDER_TRAJECTORY_NO_CLAIMS, MAX_RENDER_TRAJECTORY_SAMPLES,
        MAX_RENDER_TRANSITIONS_PER_SAMPLE, RenderBaseFrame, RenderBaseModeState,
        RenderChannelAvailability, RenderContactBranch, RenderContactGeometry,
        RenderContactTransition, RenderMassProperties, RenderNumericalRefusalReason,
        RenderSampleDisposition, RenderSupportFeature, RenderTerminalEvent, RenderTrajectory,
        RenderTrajectoryAuthority, RenderTrajectoryError, RenderTrajectoryMetadata,
        RenderTrajectorySampleInput, RenderUnitSystem, RenderWorldFrame,
    },
    timeline_resampling::{
        DeclaredDiscontinuityKind, DeclaredTimelineDiscontinuity,
        EULER_TIMELINE_RESAMPLER_VERSION,
    },
};

/// Canonical wire-schema version.
pub const EULER_RENDER_TRAJECTORY_CODEC_VERSION: u16 = 1;
/// Frozen chunking policy. Every non-final chunk has exactly this many samples.
pub const EULER_RENDER_TRAJECTORY_SAMPLES_PER_CHUNK: usize = 1_024;
/// Hard transport ceiling. Streaming paths do not allocate this amount.
pub const MAX_RENDER_TRAJECTORY_ARTIFACT_BYTES: u64 = 32 * 1_024 * 1_024 * 1_024;
/// Hard aggregate text ceiling implied by two metadata strings and 64 no-claims.
pub const MAX_RENDER_TRAJECTORY_TEXT_BYTES: usize =
    (MAX_RENDER_TRAJECTORY_NO_CLAIMS + 2) * 1_024;
/// Hard aggregate transition ceiling implied by the semantic sample limits.
pub const MAX_RENDER_TRAJECTORY_TOTAL_TRANSITIONS: usize =
    MAX_RENDER_TRAJECTORY_SAMPLES * MAX_RENDER_TRANSITIONS_PER_SAMPLE;

/// Domain for the complete canonical artifact identity.
pub const EULER_RENDER_TRAJECTORY_ARTIFACT_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.render-trajectory-artifact.v1";
/// Domain for the embedded pre-trailer integrity fingerprint.
pub const EULER_RENDER_TRAJECTORY_PAYLOAD_FINGERPRINT_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.render-trajectory-payload.v1";
/// Domain for independently checked sample-chunk integrity.
pub const EULER_RENDER_TRAJECTORY_CHUNK_FINGERPRINT_DOMAIN: &str =
    "org.frankensim.fs-euler-disc-e2e.render-trajectory-chunk.v1";

const MAGIC: &[u8; 8] = b"FSEULTRJ";
const FLOAT_POLICY_RAW_IEEE754_LE: u8 = 1;
const CHUNKING_VERSION: u16 = 1;
const INTERPOLATION_CUBIC_HERMITE_SLERP_V1: u8 = 1;
const HEADER_RESERVED: u32 = 0;
const HEADER_LEN: u64 = 116;
const PAYLOAD_FINGERPRINT_LEN: u64 = 32;
const CHUNK_HEADER_LEN: u64 = 56;
const DECLARED_DISCONTINUITY_RECORD_LEN: u64 = 9;
const MAX_METADATA_BYTES: usize = 128 * 1_024;
const MAX_SAMPLE_RECORD_BYTES: usize = 4_096;
const MAX_CHUNK_PAYLOAD_BYTES: usize =
    EULER_RENDER_TRAJECTORY_SAMPLES_PER_CHUNK * (4 + MAX_SAMPLE_RECORD_BYTES);

/// Caller-controlled decode/encode resource budget under hard schema ceilings.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RenderTrajectoryCodecBudget {
    /// Maximum complete artifact bytes accepted or emitted.
    pub max_artifact_bytes: u64,
    /// Maximum retained samples.
    pub max_samples: usize,
    /// Maximum aggregate localized contact transitions.
    pub max_total_transitions: usize,
    /// Maximum aggregate UTF-8 metadata text bytes.
    pub max_total_text_bytes: usize,
}

impl RenderTrajectoryCodecBudget {
    /// Production default. It preserves the semantic sample ceiling while
    /// bounding unusually event-dense or text-heavy inputs independently.
    pub const DEFAULT: Self = Self {
        max_artifact_bytes: MAX_RENDER_TRAJECTORY_ARTIFACT_BYTES,
        max_samples: MAX_RENDER_TRAJECTORY_SAMPLES,
        max_total_transitions: 1_000_000,
        max_total_text_bytes: MAX_RENDER_TRAJECTORY_TEXT_BYTES,
    };
}

impl Default for RenderTrajectoryCodecBudget {
    fn default() -> Self {
        Self::DEFAULT
    }
}

/// Integrity and identity facts for one complete canonical artifact.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct RenderTrajectoryCodecReceipt {
    artifact_identity: ContentHash,
    payload_fingerprint: ContentHash,
    source_campaign_identity: ContentHash,
    byte_len: u64,
    sample_count: u32,
    transition_count: u32,
    chunk_count: u32,
}

impl RenderTrajectoryCodecReceipt {
    /// Domain-separated identity of the complete bytes, including the trailer.
    #[must_use]
    pub const fn artifact_identity(self) -> ContentHash {
        self.artifact_identity
    }

    /// Embedded fingerprint of every canonical byte preceding the trailer.
    #[must_use]
    pub const fn payload_fingerprint(self) -> ContentHash {
        self.payload_fingerprint
    }

    /// Declared upstream campaign/operation identity bound by the envelope.
    #[must_use]
    pub const fn source_campaign_identity(self) -> ContentHash {
        self.source_campaign_identity
    }

    /// Exact complete wire length.
    #[must_use]
    pub const fn byte_len(self) -> u64 {
        self.byte_len
    }

    /// Exact retained sample count.
    #[must_use]
    pub const fn sample_count(self) -> u32 {
        self.sample_count
    }

    /// Exact aggregate contact-transition count.
    #[must_use]
    pub const fn transition_count(self) -> u32 {
        self.transition_count
    }

    /// Exact canonical chunk count.
    #[must_use]
    pub const fn chunk_count(self) -> u32 {
        self.chunk_count
    }
}

/// Owned admitted artifact. It is safe to move; derived controls borrow its
/// trajectory only after construction and are never stored self-referentially.
#[derive(Clone, Debug, PartialEq)]
pub struct EulerRenderTrajectoryArtifact {
    trajectory: RenderTrajectory,
    source_campaign_identity: ContentHash,
    declared_discontinuities: Vec<DeclaredTimelineDiscontinuity>,
    receipt: RenderTrajectoryCodecReceipt,
}

impl EulerRenderTrajectoryArtifact {
    /// Bind an already-admitted trajectory to a declared source campaign and
    /// timeline composition, computing its canonical receipt without retaining
    /// a second monolithic byte buffer.
    pub fn try_from_trajectory(
        source_campaign_identity: ContentHash,
        trajectory: RenderTrajectory,
        declared_discontinuities: Vec<DeclaredTimelineDiscontinuity>,
        budget: RenderTrajectoryCodecBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, RenderTrajectoryCodecError> {
        validate_context(
            source_campaign_identity,
            &trajectory,
            &declared_discontinuities,
        )?;
        let mut sink = io::sink();
        let receipt = encode_to_writer(
            &trajectory,
            source_campaign_identity,
            &declared_discontinuities,
            budget,
            &mut sink,
            &mut || checkpoint(cx),
        )?;
        checkpoint(cx)?;
        Ok(Self {
            trajectory,
            source_campaign_identity,
            declared_discontinuities,
            receipt,
        })
    }

    /// Decode, integrity-check, semantically re-admit, and byte-compare a
    /// complete seekable artifact. The reader is left at artifact end.
    pub fn read_from<R: Read + Seek>(
        reader: &mut R,
        budget: RenderTrajectoryCodecBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, RenderTrajectoryCodecError> {
        decode_from_reader(reader, budget, &mut || checkpoint(cx))
    }

    /// Decode a complete in-memory artifact without requiring a second retained
    /// copy in the resulting object.
    pub fn from_canonical_bytes(
        bytes: &[u8],
        budget: RenderTrajectoryCodecBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, RenderTrajectoryCodecError> {
        let mut cursor = Cursor::new(bytes);
        Self::read_from(&mut cursor, budget, cx)
    }

    /// Decode and additionally require an expected out-of-band artifact root.
    pub fn from_canonical_bytes_verified(
        bytes: &[u8],
        expected_identity: ContentHash,
        budget: RenderTrajectoryCodecBudget,
        cx: &Cx<'_>,
    ) -> Result<Self, RenderTrajectoryCodecError> {
        let artifact = Self::from_canonical_bytes(bytes, budget, cx)?;
        if artifact.receipt.artifact_identity != expected_identity {
            return Err(RenderTrajectoryCodecError::ArtifactIdentityMismatch {
                expected: expected_identity,
                actual: artifact.receipt.artifact_identity,
            });
        }
        Ok(artifact)
    }

    /// Stream the exact canonical representation to a writer.
    pub fn write_to<W: Write>(
        &self,
        writer: &mut W,
        budget: RenderTrajectoryCodecBudget,
        cx: &Cx<'_>,
    ) -> Result<RenderTrajectoryCodecReceipt, RenderTrajectoryCodecError> {
        let receipt = encode_to_writer(
            &self.trajectory,
            self.source_campaign_identity,
            &self.declared_discontinuities,
            budget,
            writer,
            &mut || checkpoint(cx),
        )?;
        if receipt != self.receipt {
            return Err(RenderTrajectoryCodecError::ReceiptMismatch);
        }
        Ok(receipt)
    }

    /// Convenience in-memory encoding for bounded callers and tests. Large
    /// artifacts should use [`Self::write_to`].
    pub fn canonical_bytes(
        &self,
        budget: RenderTrajectoryCodecBudget,
        cx: &Cx<'_>,
    ) -> Result<Vec<u8>, RenderTrajectoryCodecError> {
        let requested = usize::try_from(self.receipt.byte_len).map_err(|_| {
            RenderTrajectoryCodecError::Capacity {
                artifact: "canonical artifact bytes",
                requested: self.receipt.byte_len,
            }
        })?;
        let mut bytes = Vec::new();
        bytes.try_reserve_exact(requested).map_err(|_| {
            RenderTrajectoryCodecError::Capacity {
                artifact: "canonical artifact bytes",
                requested: self.receipt.byte_len,
            }
        })?;
        self.write_to(&mut bytes, budget, cx)?;
        if bytes.len() != requested {
            return Err(RenderTrajectoryCodecError::NonCanonical);
        }
        Ok(bytes)
    }

    /// Decoded, semantically admitted trajectory used by render/control clients.
    #[must_use]
    pub const fn trajectory(&self) -> &RenderTrajectory {
        &self.trajectory
    }

    /// Consume the envelope and return its admitted trajectory.
    #[must_use]
    pub fn into_trajectory(self) -> RenderTrajectory {
        self.trajectory
    }

    /// Declared upstream campaign/operation identity (not a path).
    #[must_use]
    pub const fn source_campaign_identity(&self) -> ContentHash {
        self.source_campaign_identity
    }

    /// Ordered producer-declared timeline seams bound by this artifact.
    #[must_use]
    pub fn declared_discontinuities(&self) -> &[DeclaredTimelineDiscontinuity] {
        &self.declared_discontinuities
    }

    /// Complete integrity and identity receipt.
    #[must_use]
    pub const fn receipt(&self) -> RenderTrajectoryCodecReceipt {
        self.receipt
    }
}

/// Typed refusal from the trajectory transport boundary.
#[derive(Clone, Debug, PartialEq)]
pub enum RenderTrajectoryCodecError {
    /// The execution scope requested cancellation.
    Cancelled,
    /// A caller budget was zero or exceeded a hard schema ceiling.
    InvalidBudget(&'static str),
    /// The declared source campaign identity was all zeroes.
    ZeroSourceCampaignIdentity,
    /// The complete transport exceeded its admitted byte budget.
    ArtifactTooLarge {
        /// Observed or projected bytes.
        bytes: u64,
        /// Active caller limit.
        maximum: u64,
    },
    /// A bounded allocation could not be reserved.
    Capacity {
        /// Stable artifact/component name.
        artifact: &'static str,
        /// Requested items or bytes.
        requested: u64,
    },
    /// Seek/read/write failure at the durable transport boundary.
    Io {
        /// Stable operation label.
        operation: &'static str,
        /// Operating-system-independent error class.
        kind: io::ErrorKind,
    },
    /// Input ended before a complete canonical field was available.
    Truncated {
        /// Stable field/record label.
        field: &'static str,
        /// Byte offset at which the read began.
        offset: u64,
    },
    /// Fixed magic did not identify this artifact family.
    InvalidMagic,
    /// The wire schema is unsupported.
    UnsupportedCodecVersion(u16),
    /// A pinned semantic/consumer version or policy did not match v1.
    ContractMismatch(&'static str),
    /// A closed enum, boolean, option, or reserved field used another value.
    InvalidTag {
        /// Stable field name.
        field: &'static str,
        /// Refused numeric tag.
        tag: u64,
    },
    /// A text field was not valid UTF-8.
    InvalidUtf8(&'static str),
    /// A primitive decoded value could not construct its validated domain type.
    InvalidValue(&'static str),
    /// A length/count was impossible, over budget, or inconsistent.
    InvalidLength {
        /// Stable field name.
        field: &'static str,
        /// Refused value.
        value: u64,
        /// Active maximum or exact expected value.
        maximum: u64,
    },
    /// Chunks were missing, duplicated, reordered, or not canonically sized.
    InvalidChunk {
        /// Observed chunk ordinal.
        chunk: u32,
        /// Stable violated descriptor field.
        field: &'static str,
    },
    /// A per-chunk digest did not match its descriptor and payload.
    ChunkFingerprintMismatch(u32),
    /// The embedded whole-prefix fingerprint was stale or corrupt.
    PayloadFingerprintMismatch,
    /// Expected and decoded out-of-band artifact identities differed.
    ArtifactIdentityMismatch {
        /// Caller-supplied expected root.
        expected: ContentHash,
        /// Root of the decoded complete bytes.
        actual: ContentHash,
    },
    /// A decoded representation was not the unique canonical byte fixed point.
    NonCanonical,
    /// Re-encoding no longer matched the retained receipt.
    ReceiptMismatch,
    /// A declared timeline seam was invalid for the admitted trajectory.
    InvalidDeclaredDiscontinuity(usize),
    /// Semantic trajectory admission refused the decoded payload.
    Trajectory(RenderTrajectoryError),
}

impl fmt::Display for RenderTrajectoryCodecError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for RenderTrajectoryCodecError {}

impl From<RenderTrajectoryError> for RenderTrajectoryCodecError {
    fn from(error: RenderTrajectoryError) -> Self {
        Self::Trajectory(error)
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), RenderTrajectoryCodecError> {
    cx.checkpoint()
        .map_err(|_| RenderTrajectoryCodecError::Cancelled)
}

trait CanonicalSink {
    fn put(&mut self, bytes: &[u8]) -> Result<(), RenderTrajectoryCodecError>;

    fn u8(&mut self, value: u8) -> Result<(), RenderTrajectoryCodecError> {
        self.put(&[value])
    }

    fn u16(&mut self, value: u16) -> Result<(), RenderTrajectoryCodecError> {
        self.put(&value.to_le_bytes())
    }

    fn u32(&mut self, value: u32) -> Result<(), RenderTrajectoryCodecError> {
        self.put(&value.to_le_bytes())
    }

    fn u64(&mut self, value: u64) -> Result<(), RenderTrajectoryCodecError> {
        self.put(&value.to_le_bytes())
    }

    fn f64(&mut self, value: f64) -> Result<(), RenderTrajectoryCodecError> {
        self.u64(value.to_bits())
    }

    fn hash(&mut self, value: ContentHash) -> Result<(), RenderTrajectoryCodecError> {
        self.put(value.as_bytes())
    }

    fn vec3(&mut self, value: Vec3) -> Result<(), RenderTrajectoryCodecError> {
        self.f64(value.x)?;
        self.f64(value.y)?;
        self.f64(value.z)
    }

    fn quaternion(
        &mut self,
        value: UnitQuaternion,
    ) -> Result<(), RenderTrajectoryCodecError> {
        for component in value.components() {
            self.f64(component)?;
        }
        Ok(())
    }

    fn string(&mut self, value: &str) -> Result<(), RenderTrajectoryCodecError> {
        let len = u32::try_from(value.len()).map_err(|_| {
            RenderTrajectoryCodecError::InvalidLength {
                field: "string",
                value: u64::try_from(value.len()).unwrap_or(u64::MAX),
                maximum: u64::from(u32::MAX),
            }
        })?;
        self.u32(len)?;
        self.put(value.as_bytes())
    }
}

#[derive(Default)]
struct SizeSink {
    len: u64,
}

impl CanonicalSink for SizeSink {
    fn put(&mut self, bytes: &[u8]) -> Result<(), RenderTrajectoryCodecError> {
        let added = u64::try_from(bytes.len()).map_err(|_| {
            RenderTrajectoryCodecError::InvalidLength {
                field: "canonical size",
                value: u64::MAX,
                maximum: u64::MAX,
            }
        })?;
        self.len = self.len.checked_add(added).ok_or(
            RenderTrajectoryCodecError::InvalidLength {
                field: "canonical size",
                value: u64::MAX,
                maximum: u64::MAX,
            },
        )?;
        Ok(())
    }
}

struct VecSink {
    bytes: Vec<u8>,
    maximum: usize,
    artifact: &'static str,
}

impl VecSink {
    fn with_exact_capacity(
        capacity: usize,
        maximum: usize,
        artifact: &'static str,
    ) -> Result<Self, RenderTrajectoryCodecError> {
        if capacity > maximum {
            return Err(RenderTrajectoryCodecError::InvalidLength {
                field: artifact,
                value: u64::try_from(capacity).unwrap_or(u64::MAX),
                maximum: u64::try_from(maximum).unwrap_or(u64::MAX),
            });
        }
        let mut bytes = Vec::new();
        bytes.try_reserve_exact(capacity).map_err(|_| {
            RenderTrajectoryCodecError::Capacity {
                artifact,
                requested: u64::try_from(capacity).unwrap_or(u64::MAX),
            }
        })?;
        Ok(Self {
            bytes,
            maximum,
            artifact,
        })
    }

    fn finish(self, expected: usize) -> Result<Vec<u8>, RenderTrajectoryCodecError> {
        if self.bytes.len() != expected {
            return Err(RenderTrajectoryCodecError::NonCanonical);
        }
        Ok(self.bytes)
    }
}

impl CanonicalSink for VecSink {
    fn put(&mut self, bytes: &[u8]) -> Result<(), RenderTrajectoryCodecError> {
        let new_len = self.bytes.len().checked_add(bytes.len()).ok_or(
            RenderTrajectoryCodecError::InvalidLength {
                field: self.artifact,
                value: u64::MAX,
                maximum: u64::try_from(self.maximum).unwrap_or(u64::MAX),
            },
        )?;
        if new_len > self.maximum {
            return Err(RenderTrajectoryCodecError::InvalidLength {
                field: self.artifact,
                value: u64::try_from(new_len).unwrap_or(u64::MAX),
                maximum: u64::try_from(self.maximum).unwrap_or(u64::MAX),
            });
        }
        if new_len > self.bytes.capacity() {
            self.bytes.try_reserve(new_len - self.bytes.len()).map_err(|_| {
                RenderTrajectoryCodecError::Capacity {
                    artifact: self.artifact,
                    requested: u64::try_from(new_len).unwrap_or(u64::MAX),
                }
            })?;
        }
        self.bytes.extend_from_slice(bytes);
        Ok(())
    }
}

fn encode_metadata<S: CanonicalSink>(
    metadata: &RenderTrajectoryMetadata,
    sink: &mut S,
) -> Result<(), RenderTrajectoryCodecError> {
    sink.u16(metadata.schema_version)?;
    sink.u8(world_frame_tag(metadata.world_frame))?;
    sink.u8(unit_system_tag(metadata.units))?;
    sink.u8(authority_tag(metadata.authority))?;
    sink.u8(availability_bits(metadata.channel_availability))?;
    sink.hash(metadata.specimen_profile_identity)?;
    sink.hash(metadata.specimen_chart_identity)?;
    sink.hash(metadata.mass_properties.identity)?;
    let mass = metadata.mass_properties.properties;
    sink.f64(mass.mass())?;
    sink.vec3(mass.center_of_mass_body())?;
    sink.vec3(mass.principal_inertia_body())?;
    encode_rigid_state(metadata.initial_state, sink)?;
    sink.f64(metadata.initial_base_mode.displacement_m)?;
    sink.f64(metadata.initial_base_mode.velocity_m_per_s)?;
    sink.hash(metadata.base_model_identity)?;
    sink.vec3(metadata.base_frame.origin_world_m)?;
    sink.quaternion(metadata.base_frame.orientation_base_to_world)?;
    sink.hash(metadata.model_identity)?;
    sink.hash(metadata.configuration_identity)?;
    sink.u64(metadata.configuration_fingerprint)?;
    sink.f64(metadata.timestep_s)?;
    sink.string(&metadata.producer_version)?;
    sink.string(&metadata.applicability)?;
    sink.u32(u32_from_usize("metadata.no_claims", metadata.no_claims.len())?)?;
    for no_claim in &metadata.no_claims {
        sink.string(no_claim)?;
    }
    Ok(())
}

fn encode_rigid_state<S: CanonicalSink>(
    state: RigidBodyState,
    sink: &mut S,
) -> Result<(), RenderTrajectoryCodecError> {
    sink.vec3(state.pose().position_world())?;
    sink.quaternion(state.pose().orientation())?;
    sink.vec3(state.linear_momentum_world())?;
    sink.vec3(state.angular_momentum_body())
}

fn encode_sample<S: CanonicalSink>(
    input: &RenderTrajectorySampleInput,
    sink: &mut S,
) -> Result<(), RenderTrajectoryCodecError> {
    sink.f64(input.interval_start_time_s)?;
    sink.f64(input.time_s)?;
    sink.u8(world_frame_tag(input.world_frame))?;
    sink.u8(unit_system_tag(input.units))?;
    sink.vec3(input.center_of_mass_world_m)?;
    for component in input.orientation_body_to_world {
        sink.f64(component)?;
    }
    sink.vec3(input.linear_momentum_world_kg_m_per_s)?;
    sink.vec3(input.angular_momentum_body_kg_m2_per_s)?;
    sink.vec3(input.symmetry_axis_world)?;
    sink.u8(contact_branch_tag(input.contact_branch))?;
    match input.contact_geometry {
        None => sink.u8(0)?,
        Some(geometry) => {
            sink.u8(1)?;
            sink.vec3(geometry.point_world_m)?;
            sink.vec3(geometry.normal_world)?;
            match geometry.support_feature {
                RenderSupportFeature::CylinderRim => sink.u8(1)?,
                RenderSupportFeature::ProfileFeature(index) => {
                    sink.u8(2)?;
                    sink.u64(u64::try_from(index).map_err(|_| {
                        RenderTrajectoryCodecError::InvalidLength {
                            field: "support feature index",
                            value: u64::MAX,
                            maximum: u64::MAX,
                        }
                    })?)?;
                }
            }
        }
    }
    sink.f64(input.signed_gap_m)?;
    sink.u8(u8::from(input.interval_contact_active))?;
    sink.f64(input.interval_normal_force_n)?;
    sink.u32(u32_from_usize(
        "sample.contact_transitions",
        input.contact_transitions.len(),
    )?)?;
    for transition in &input.contact_transitions {
        sink.u8(transition_kind_tag(transition.kind))?;
        sink.f64(transition.time_s)?;
        sink.f64(transition.bracket_start_s)?;
        sink.f64(transition.bracket_end_s)?;
    }
    match input.base_mode {
        None => sink.u8(0)?,
        Some(base) => {
            sink.u8(1)?;
            sink.f64(base.displacement_m)?;
            sink.f64(base.velocity_m_per_s)?;
        }
    }
    encode_channels(input.channels, sink)?;
    sink.f64(input.mechanical_energy_j)?;
    sink.f64(input.energy_defect_j)?;
    sink.f64(input.qois.inclination_rad)?;
    sink.f64(input.qois.precession_rad_per_s)?;
    sink.f64(input.qois.spin_rad_per_s)?;
    sink.f64(input.qois.precession_acceleration_rad_per_s2)?;
    let (disposition, backend_code) = disposition_encoding(input.disposition);
    sink.u8(disposition)?;
    sink.u32(backend_code)?;
    match input.terminal_event {
        None => sink.u8(0)?,
        Some(event) => {
            sink.u8(1)?;
            sink.f64(event.time_s)?;
            sink.f64(event.bracket_start_s)?;
            sink.f64(event.bracket_end_s)?;
        }
    }
    Ok(())
}

fn encode_channels<S: CanonicalSink>(
    channels: ChannelOwnership,
    sink: &mut S,
) -> Result<(), RenderTrajectoryCodecError> {
    for channel in [
        channels.gravity,
        channels.contact,
        channels.rolling,
        channels.base,
        channels.gas,
    ] {
        sink.vec3(channel.force_world_n)?;
        sink.vec3(channel.torque_world_nm)?;
        sink.f64(channel.work_j)?;
    }
    Ok(())
}

fn world_frame_tag(value: RenderWorldFrame) -> u8 {
    match value {
        RenderWorldFrame::RightHandedZUp => 1,
        RenderWorldFrame::RightHandedYUp => 2,
    }
}

fn unit_system_tag(value: RenderUnitSystem) -> u8 {
    match value {
        RenderUnitSystem::SiRadians => 1,
        RenderUnitSystem::SiDegrees => 2,
    }
}

const fn authority_tag(value: RenderTrajectoryAuthority) -> u8 {
    match value {
        RenderTrajectoryAuthority::SimulationEvidence => 1,
    }
}

fn availability_bits(value: RenderChannelAvailability) -> u8 {
    u8::from(value.gravity)
        | (u8::from(value.contact) << 1)
        | (u8::from(value.rolling) << 2)
        | (u8::from(value.base) << 3)
        | (u8::from(value.gas) << 4)
}

const fn contact_branch_tag(value: RenderContactBranch) -> u8 {
    match value {
        RenderContactBranch::Open => 1,
        RenderContactBranch::Closed => 2,
    }
}

const fn transition_kind_tag(value: ContactTransitionKind) -> u8 {
    match value {
        ContactTransitionKind::Opening => 1,
        ContactTransitionKind::Reimpact => 2,
    }
}

const fn declared_discontinuity_kind_tag(value: DeclaredDiscontinuityKind) -> u8 {
    match value {
        DeclaredDiscontinuityKind::ContinuationSeam => 1,
        DeclaredDiscontinuityKind::ProducerDeclared => 2,
    }
}

const fn disposition_encoding(value: RenderSampleDisposition) -> (u8, u32) {
    match value {
        RenderSampleDisposition::Continue => (0, 0),
        RenderSampleDisposition::TerminalInclination => (1, 0),
        RenderSampleDisposition::HorizonCensored => (2, 0),
        RenderSampleDisposition::NumericalRefusal(
            RenderNumericalRefusalReason::ReimpactLimitExceeded,
        ) => (3, 0),
        RenderSampleDisposition::NumericalRefusal(
            RenderNumericalRefusalReason::ContactEventLocalizationFailed,
        ) => (4, 0),
        RenderSampleDisposition::NumericalRefusal(
            RenderNumericalRefusalReason::NonFiniteEnergyOrBaseState,
        ) => (5, 0),
        RenderSampleDisposition::NumericalRefusal(
            RenderNumericalRefusalReason::BackendSpecific(code),
        ) => (6, code),
    }
}

fn u32_from_usize(
    field: &'static str,
    value: usize,
) -> Result<u32, RenderTrajectoryCodecError> {
    u32::try_from(value).map_err(|_| RenderTrajectoryCodecError::InvalidLength {
        field,
        value: u64::try_from(value).unwrap_or(u64::MAX),
        maximum: u64::from(u32::MAX),
    })
}

#[derive(Clone, Copy, Debug)]
struct WirePlan {
    total_len: u64,
    metadata_len: u32,
    discontinuity_count: u32,
    discontinuity_len: u64,
    sample_count: u32,
    transition_count: u32,
    chunk_count: u32,
    first_time_s: f64,
    last_time_s: f64,
    terminal_tag: u8,
}

#[derive(Clone, Copy, Debug)]
struct Header {
    plan: WirePlan,
    source_campaign_identity: ContentHash,
    availability: u8,
    world_frame: u8,
    units: u8,
}

fn validate_budget(budget: RenderTrajectoryCodecBudget) -> Result<(), RenderTrajectoryCodecError> {
    if budget.max_artifact_bytes == 0
        || budget.max_artifact_bytes > MAX_RENDER_TRAJECTORY_ARTIFACT_BYTES
    {
        return Err(RenderTrajectoryCodecError::InvalidBudget(
            "max_artifact_bytes",
        ));
    }
    if budget.max_samples == 0 || budget.max_samples > MAX_RENDER_TRAJECTORY_SAMPLES {
        return Err(RenderTrajectoryCodecError::InvalidBudget("max_samples"));
    }
    if budget.max_total_transitions > MAX_RENDER_TRAJECTORY_TOTAL_TRANSITIONS {
        return Err(RenderTrajectoryCodecError::InvalidBudget(
            "max_total_transitions",
        ));
    }
    if budget.max_total_text_bytes == 0
        || budget.max_total_text_bytes > MAX_RENDER_TRAJECTORY_TEXT_BYTES
    {
        return Err(RenderTrajectoryCodecError::InvalidBudget(
            "max_total_text_bytes",
        ));
    }
    Ok(())
}

fn validate_context(
    source_campaign_identity: ContentHash,
    trajectory: &RenderTrajectory,
    declared_discontinuities: &[DeclaredTimelineDiscontinuity],
) -> Result<(), RenderTrajectoryCodecError> {
    if source_campaign_identity
        .as_bytes()
        .iter()
        .all(|byte| *byte == 0)
    {
        return Err(RenderTrajectoryCodecError::ZeroSourceCampaignIdentity);
    }
    let samples = trajectory.samples();
    let first_time_s = samples[0].input().time_s;
    let last_time_s = samples[samples.len() - 1].input().time_s;
    let mut previous_time_s = None;
    for (index, discontinuity) in declared_discontinuities.iter().enumerate() {
        let matches_sample = samples
            .iter()
            .any(|sample| same_time(sample.input().time_s, discontinuity.time_s));
        if !discontinuity.time_s.is_finite()
            || discontinuity.time_s < first_time_s
            || discontinuity.time_s > last_time_s
            || previous_time_s.is_some_and(|time| discontinuity.time_s <= time)
            || !matches_sample
        {
            return Err(RenderTrajectoryCodecError::InvalidDeclaredDiscontinuity(
                index,
            ));
        }
        previous_time_s = Some(discontinuity.time_s);
    }
    Ok(())
}

fn same_time(first: f64, second: f64) -> bool {
    let first_bits = first.to_bits();
    let second_bits = second.to_bits();
    first_bits == second_bits || (first_bits << 1 == 0 && second_bits << 1 == 0)
}

fn text_byte_count(metadata: &RenderTrajectoryMetadata) -> Result<usize, RenderTrajectoryCodecError> {
    let mut total = metadata
        .producer_version
        .len()
        .checked_add(metadata.applicability.len())
        .ok_or(RenderTrajectoryCodecError::InvalidLength {
            field: "metadata text bytes",
            value: u64::MAX,
            maximum: u64::try_from(MAX_RENDER_TRAJECTORY_TEXT_BYTES).unwrap_or(u64::MAX),
        })?;
    for no_claim in &metadata.no_claims {
        total = total.checked_add(no_claim.len()).ok_or(
            RenderTrajectoryCodecError::InvalidLength {
                field: "metadata text bytes",
                value: u64::MAX,
                maximum: u64::try_from(MAX_RENDER_TRAJECTORY_TEXT_BYTES).unwrap_or(u64::MAX),
            },
        )?;
    }
    Ok(total)
}

fn checked_len_add(
    current: u64,
    added: u64,
    field: &'static str,
) -> Result<u64, RenderTrajectoryCodecError> {
    current
        .checked_add(added)
        .ok_or(RenderTrajectoryCodecError::InvalidLength {
            field,
            value: u64::MAX,
            maximum: MAX_RENDER_TRAJECTORY_ARTIFACT_BYTES,
        })
}

fn measure_wire(
    trajectory: &RenderTrajectory,
    declared_discontinuities: &[DeclaredTimelineDiscontinuity],
    budget: RenderTrajectoryCodecBudget,
    checkpoint: &mut impl FnMut() -> Result<(), RenderTrajectoryCodecError>,
) -> Result<WirePlan, RenderTrajectoryCodecError> {
    validate_budget(budget)?;
    let samples = trajectory.samples();
    if samples.len() > budget.max_samples {
        return Err(RenderTrajectoryCodecError::InvalidLength {
            field: "sample_count",
            value: u64::try_from(samples.len()).unwrap_or(u64::MAX),
            maximum: u64::try_from(budget.max_samples).unwrap_or(u64::MAX),
        });
    }
    let text_bytes = text_byte_count(trajectory.metadata())?;
    if text_bytes > budget.max_total_text_bytes {
        return Err(RenderTrajectoryCodecError::InvalidLength {
            field: "metadata text bytes",
            value: u64::try_from(text_bytes).unwrap_or(u64::MAX),
            maximum: u64::try_from(budget.max_total_text_bytes).unwrap_or(u64::MAX),
        });
    }

    let mut metadata_size = SizeSink::default();
    encode_metadata(trajectory.metadata(), &mut metadata_size)?;
    let metadata_len = u32::try_from(metadata_size.len).map_err(|_| {
        RenderTrajectoryCodecError::InvalidLength {
            field: "metadata_len",
            value: metadata_size.len,
            maximum: u64::from(u32::MAX),
        }
    })?;
    if metadata_size.len > u64::try_from(MAX_METADATA_BYTES).unwrap_or(u64::MAX) {
        return Err(RenderTrajectoryCodecError::InvalidLength {
            field: "metadata_len",
            value: metadata_size.len,
            maximum: u64::try_from(MAX_METADATA_BYTES).unwrap_or(u64::MAX),
        });
    }

    let discontinuity_count = u32_from_usize(
        "declared_discontinuity_count",
        declared_discontinuities.len(),
    )?;
    let discontinuity_len = u64::from(discontinuity_count)
        .checked_mul(DECLARED_DISCONTINUITY_RECORD_LEN)
        .ok_or(RenderTrajectoryCodecError::InvalidLength {
            field: "declared_discontinuity_len",
            value: u64::MAX,
            maximum: MAX_RENDER_TRAJECTORY_ARTIFACT_BYTES,
        })?;
    let sample_count = u32_from_usize("sample_count", samples.len())?;
    let chunk_count_usize = samples.len().div_ceil(EULER_RENDER_TRAJECTORY_SAMPLES_PER_CHUNK);
    let chunk_count = u32_from_usize("chunk_count", chunk_count_usize)?;
    let mut transition_count = 0usize;
    let mut total_len = checked_len_add(HEADER_LEN, metadata_size.len, "artifact length")?;
    total_len = checked_len_add(total_len, discontinuity_len, "artifact length")?;

    for chunk in samples.chunks(EULER_RENDER_TRAJECTORY_SAMPLES_PER_CHUNK) {
        checkpoint()?;
        let mut chunk_payload_len = 0u64;
        for sample in chunk {
            let input = sample.input();
            transition_count = transition_count
                .checked_add(input.contact_transitions.len())
                .ok_or(RenderTrajectoryCodecError::InvalidLength {
                    field: "transition_count",
                    value: u64::MAX,
                    maximum: u64::try_from(budget.max_total_transitions).unwrap_or(u64::MAX),
                })?;
            if transition_count > budget.max_total_transitions {
                return Err(RenderTrajectoryCodecError::InvalidLength {
                    field: "transition_count",
                    value: u64::try_from(transition_count).unwrap_or(u64::MAX),
                    maximum: u64::try_from(budget.max_total_transitions).unwrap_or(u64::MAX),
                });
            }
            let mut sample_size = SizeSink::default();
            encode_sample(input, &mut sample_size)?;
            if sample_size.len > u64::try_from(MAX_SAMPLE_RECORD_BYTES).unwrap_or(u64::MAX) {
                return Err(RenderTrajectoryCodecError::InvalidLength {
                    field: "sample_record_len",
                    value: sample_size.len,
                    maximum: u64::try_from(MAX_SAMPLE_RECORD_BYTES).unwrap_or(u64::MAX),
                });
            }
            chunk_payload_len = checked_len_add(chunk_payload_len, 4, "chunk payload length")?;
            chunk_payload_len = checked_len_add(
                chunk_payload_len,
                sample_size.len,
                "chunk payload length",
            )?;
        }
        if chunk_payload_len > u64::try_from(MAX_CHUNK_PAYLOAD_BYTES).unwrap_or(u64::MAX) {
            return Err(RenderTrajectoryCodecError::InvalidLength {
                field: "chunk_payload_len",
                value: chunk_payload_len,
                maximum: u64::try_from(MAX_CHUNK_PAYLOAD_BYTES).unwrap_or(u64::MAX),
            });
        }
        total_len = checked_len_add(total_len, CHUNK_HEADER_LEN, "artifact length")?;
        total_len = checked_len_add(total_len, chunk_payload_len, "artifact length")?;
    }
    total_len = checked_len_add(total_len, PAYLOAD_FINGERPRINT_LEN, "artifact length")?;
    if total_len > budget.max_artifact_bytes {
        return Err(RenderTrajectoryCodecError::ArtifactTooLarge {
            bytes: total_len,
            maximum: budget.max_artifact_bytes,
        });
    }
    let transition_count = u32_from_usize("transition_count", transition_count)?;
    let first_time_s = samples[0].input().time_s;
    let last = samples[samples.len() - 1].input();
    Ok(WirePlan {
        total_len,
        metadata_len,
        discontinuity_count,
        discontinuity_len,
        sample_count,
        transition_count,
        chunk_count,
        first_time_s,
        last_time_s: last.time_s,
        terminal_tag: disposition_encoding(last.disposition).0,
    })
}

fn encode_header(
    plan: WirePlan,
    source_campaign_identity: ContentHash,
    metadata: &RenderTrajectoryMetadata,
) -> Result<Vec<u8>, RenderTrajectoryCodecError> {
    let capacity = usize::try_from(HEADER_LEN).map_err(|_| {
        RenderTrajectoryCodecError::Capacity {
            artifact: "header",
            requested: HEADER_LEN,
        }
    })?;
    let mut sink = VecSink::with_exact_capacity(capacity, capacity, "header")?;
    sink.put(MAGIC)?;
    sink.u16(EULER_RENDER_TRAJECTORY_CODEC_VERSION)?;
    sink.u16(EULER_RENDER_TRAJECTORY_SCHEMA_VERSION)?;
    sink.u16(EULER_CONTROL_STREAM_SCHEMA_VERSION)?;
    sink.u16(EULER_TIMELINE_RESAMPLER_VERSION)?;
    sink.u8(FLOAT_POLICY_RAW_IEEE754_LE)?;
    sink.u8(INTERPOLATION_CUBIC_HERMITE_SLERP_V1)?;
    sink.u16(CHUNKING_VERSION)?;
    sink.u32(u32_from_usize(
        "samples_per_chunk",
        EULER_RENDER_TRAJECTORY_SAMPLES_PER_CHUNK,
    )?)?;
    sink.u32(HEADER_RESERVED)?;
    sink.u64(plan.total_len)?;
    sink.u32(plan.metadata_len)?;
    sink.u32(plan.discontinuity_count)?;
    sink.u64(plan.discontinuity_len)?;
    sink.u32(plan.sample_count)?;
    sink.u32(plan.transition_count)?;
    sink.u32(plan.chunk_count)?;
    sink.f64(plan.first_time_s)?;
    sink.f64(plan.last_time_s)?;
    sink.hash(source_campaign_identity)?;
    sink.u8(plan.terminal_tag)?;
    sink.u8(availability_bits(metadata.channel_availability))?;
    sink.u8(world_frame_tag(metadata.world_frame))?;
    sink.u8(unit_system_tag(metadata.units))?;
    sink.finish(capacity)
}

fn chunk_descriptor(
    chunk_index: u32,
    first_sample_index: u32,
    sample_count: u32,
    transition_count: u32,
    payload_len: u64,
) -> Result<Vec<u8>, RenderTrajectoryCodecError> {
    let mut sink = VecSink::with_exact_capacity(24, 24, "chunk descriptor")?;
    sink.u32(chunk_index)?;
    sink.u32(first_sample_index)?;
    sink.u32(sample_count)?;
    sink.u32(transition_count)?;
    sink.u64(payload_len)?;
    sink.finish(24)
}

struct ArtifactWriter<'writer, W> {
    writer: &'writer mut W,
    prefix_hasher: DomainHasher,
    artifact_hasher: DomainHasher,
    written: u64,
}

impl<'writer, W: Write> ArtifactWriter<'writer, W> {
    fn new(writer: &'writer mut W) -> Self {
        Self {
            writer,
            prefix_hasher: DomainHasher::new(EULER_RENDER_TRAJECTORY_PAYLOAD_FINGERPRINT_DOMAIN),
            artifact_hasher: DomainHasher::new(EULER_RENDER_TRAJECTORY_ARTIFACT_IDENTITY_DOMAIN),
            written: 0,
        }
    }

    fn prefix(&mut self, bytes: &[u8]) -> Result<(), RenderTrajectoryCodecError> {
        write_all(self.writer, bytes, "write artifact")?;
        self.prefix_hasher.update(bytes);
        self.artifact_hasher.update(bytes);
        self.written = checked_len_add(
            self.written,
            u64::try_from(bytes.len()).unwrap_or(u64::MAX),
            "written artifact length",
        )?;
        Ok(())
    }

    fn trailer(&mut self, bytes: &[u8]) -> Result<(), RenderTrajectoryCodecError> {
        write_all(self.writer, bytes, "write fingerprint trailer")?;
        self.artifact_hasher.update(bytes);
        self.written = checked_len_add(
            self.written,
            u64::try_from(bytes.len()).unwrap_or(u64::MAX),
            "written artifact length",
        )?;
        Ok(())
    }
}

fn write_all(
    writer: &mut impl Write,
    bytes: &[u8],
    operation: &'static str,
) -> Result<(), RenderTrajectoryCodecError> {
    writer
        .write_all(bytes)
        .map_err(|error| RenderTrajectoryCodecError::Io {
            operation,
            kind: error.kind(),
        })
}

#[allow(
    clippy::too_many_lines,
    reason = "the canonical envelope is emitted in one visibly ordered protocol path"
)]
fn encode_to_writer<W: Write>(
    trajectory: &RenderTrajectory,
    source_campaign_identity: ContentHash,
    declared_discontinuities: &[DeclaredTimelineDiscontinuity],
    budget: RenderTrajectoryCodecBudget,
    writer: &mut W,
    checkpoint: &mut impl FnMut() -> Result<(), RenderTrajectoryCodecError>,
) -> Result<RenderTrajectoryCodecReceipt, RenderTrajectoryCodecError> {
    validate_context(
        source_campaign_identity,
        trajectory,
        declared_discontinuities,
    )?;
    let plan = measure_wire(trajectory, declared_discontinuities, budget, checkpoint)?;
    checkpoint()?;
    let header = encode_header(plan, source_campaign_identity, trajectory.metadata())?;
    let mut output = ArtifactWriter::new(writer);
    output.prefix(&header)?;

    let metadata_len = usize::try_from(plan.metadata_len).map_err(|_| {
        RenderTrajectoryCodecError::Capacity {
            artifact: "metadata",
            requested: u64::from(plan.metadata_len),
        }
    })?;
    let mut metadata = VecSink::with_exact_capacity(
        metadata_len,
        MAX_METADATA_BYTES,
        "metadata",
    )?;
    encode_metadata(trajectory.metadata(), &mut metadata)?;
    output.prefix(&metadata.finish(metadata_len)?)?;

    for discontinuity in declared_discontinuities {
        let mut record = VecSink::with_exact_capacity(9, 9, "declared discontinuity")?;
        record.f64(discontinuity.time_s)?;
        record.u8(declared_discontinuity_kind_tag(discontinuity.kind))?;
        output.prefix(&record.finish(9)?)?;
    }

    let samples = trajectory.samples();
    let mut first_sample_index = 0usize;
    for (chunk_index, chunk) in samples
        .chunks(EULER_RENDER_TRAJECTORY_SAMPLES_PER_CHUNK)
        .enumerate()
    {
        checkpoint()?;
        let mut payload_size = SizeSink::default();
        let mut chunk_transition_count = 0usize;
        for sample in chunk {
            let mut sample_size = SizeSink::default();
            encode_sample(sample.input(), &mut sample_size)?;
            payload_size.u32(u32::try_from(sample_size.len).map_err(|_| {
                RenderTrajectoryCodecError::InvalidLength {
                    field: "sample_record_len",
                    value: sample_size.len,
                    maximum: u64::from(u32::MAX),
                }
            })?)?;
            payload_size.len = checked_len_add(
                payload_size.len,
                sample_size.len,
                "chunk payload length",
            )?;
            chunk_transition_count = chunk_transition_count
                .checked_add(sample.input().contact_transitions.len())
                .ok_or(RenderTrajectoryCodecError::InvalidLength {
                    field: "chunk_transition_count",
                    value: u64::MAX,
                    maximum: u64::from(u32::MAX),
                })?;
        }
        let payload_len = usize::try_from(payload_size.len).map_err(|_| {
            RenderTrajectoryCodecError::Capacity {
                artifact: "chunk payload",
                requested: payload_size.len,
            }
        })?;
        let mut payload = VecSink::with_exact_capacity(
            payload_len,
            MAX_CHUNK_PAYLOAD_BYTES,
            "chunk payload",
        )?;
        for sample in chunk {
            let mut sample_size = SizeSink::default();
            encode_sample(sample.input(), &mut sample_size)?;
            payload.u32(u32::try_from(sample_size.len).map_err(|_| {
                RenderTrajectoryCodecError::InvalidLength {
                    field: "sample_record_len",
                    value: sample_size.len,
                    maximum: u64::from(u32::MAX),
                }
            })?)?;
            encode_sample(sample.input(), &mut payload)?;
        }
        let payload = payload.finish(payload_len)?;
        let descriptor = chunk_descriptor(
            u32_from_usize("chunk_index", chunk_index)?,
            u32_from_usize("first_sample_index", first_sample_index)?,
            u32_from_usize("chunk_sample_count", chunk.len())?,
            u32_from_usize("chunk_transition_count", chunk_transition_count)?,
            u64::try_from(payload.len()).unwrap_or(u64::MAX),
        )?;
        let mut chunk_hasher = DomainHasher::new(EULER_RENDER_TRAJECTORY_CHUNK_FINGERPRINT_DOMAIN);
        chunk_hasher.update(&descriptor);
        chunk_hasher.update(&payload);
        output.prefix(&descriptor)?;
        output.prefix(chunk_hasher.finalize().as_bytes())?;
        output.prefix(&payload)?;
        first_sample_index += chunk.len();
    }

    let payload_fingerprint = output.prefix_hasher.finalize();
    output.trailer(payload_fingerprint.as_bytes())?;
    if output.written != plan.total_len {
        return Err(RenderTrajectoryCodecError::NonCanonical);
    }
    let artifact_identity = output.artifact_hasher.finalize();
    Ok(RenderTrajectoryCodecReceipt {
        artifact_identity,
        payload_fingerprint,
        source_campaign_identity,
        byte_len: plan.total_len,
        sample_count: plan.sample_count,
        transition_count: plan.transition_count,
        chunk_count: plan.chunk_count,
    })
}

struct SliceDecoder<'bytes> {
    bytes: &'bytes [u8],
    position: usize,
    base_offset: u64,
    text_bytes: usize,
    max_text_bytes: usize,
}

impl<'bytes> SliceDecoder<'bytes> {
    const fn new(bytes: &'bytes [u8], base_offset: u64, max_text_bytes: usize) -> Self {
        Self {
            bytes,
            position: 0,
            base_offset,
            text_bytes: 0,
            max_text_bytes,
        }
    }

    fn take(
        &mut self,
        length: usize,
        field: &'static str,
    ) -> Result<&'bytes [u8], RenderTrajectoryCodecError> {
        let offset = self
            .base_offset
            .checked_add(u64::try_from(self.position).unwrap_or(u64::MAX))
            .unwrap_or(u64::MAX);
        let end = self.position.checked_add(length).ok_or(
            RenderTrajectoryCodecError::InvalidLength {
                field,
                value: u64::MAX,
                maximum: u64::try_from(self.bytes.len()).unwrap_or(u64::MAX),
            },
        )?;
        let value = self
            .bytes
            .get(self.position..end)
            .ok_or(RenderTrajectoryCodecError::Truncated { field, offset })?;
        self.position = end;
        Ok(value)
    }

    fn u8(&mut self, field: &'static str) -> Result<u8, RenderTrajectoryCodecError> {
        Ok(self.take(1, field)?[0])
    }

    fn u16(&mut self, field: &'static str) -> Result<u16, RenderTrajectoryCodecError> {
        let bytes: [u8; 2] = self.take(2, field)?.try_into().map_err(|_| {
            RenderTrajectoryCodecError::Truncated {
                field,
                offset: self.base_offset,
            }
        })?;
        Ok(u16::from_le_bytes(bytes))
    }

    fn u32(&mut self, field: &'static str) -> Result<u32, RenderTrajectoryCodecError> {
        let bytes: [u8; 4] = self.take(4, field)?.try_into().map_err(|_| {
            RenderTrajectoryCodecError::Truncated {
                field,
                offset: self.base_offset,
            }
        })?;
        Ok(u32::from_le_bytes(bytes))
    }

    fn u64(&mut self, field: &'static str) -> Result<u64, RenderTrajectoryCodecError> {
        let bytes: [u8; 8] = self.take(8, field)?.try_into().map_err(|_| {
            RenderTrajectoryCodecError::Truncated {
                field,
                offset: self.base_offset,
            }
        })?;
        Ok(u64::from_le_bytes(bytes))
    }

    fn f64(&mut self, field: &'static str) -> Result<f64, RenderTrajectoryCodecError> {
        Ok(f64::from_bits(self.u64(field)?))
    }

    fn hash(&mut self, field: &'static str) -> Result<ContentHash, RenderTrajectoryCodecError> {
        let bytes: [u8; 32] = self.take(32, field)?.try_into().map_err(|_| {
            RenderTrajectoryCodecError::Truncated {
                field,
                offset: self.base_offset,
            }
        })?;
        Ok(ContentHash(bytes))
    }

    fn vec3(&mut self, field: &'static str) -> Result<Vec3, RenderTrajectoryCodecError> {
        Ok(Vec3::new(
            self.f64(field)?,
            self.f64(field)?,
            self.f64(field)?,
        ))
    }

    fn quaternion(
        &mut self,
        field: &'static str,
    ) -> Result<UnitQuaternion, RenderTrajectoryCodecError> {
        UnitQuaternion::new(
            self.f64(field)?,
            self.f64(field)?,
            self.f64(field)?,
            self.f64(field)?,
        )
        .map_err(|_| RenderTrajectoryCodecError::InvalidValue(field))
    }

    fn boolean(&mut self, field: &'static str) -> Result<bool, RenderTrajectoryCodecError> {
        match self.u8(field)? {
            0 => Ok(false),
            1 => Ok(true),
            tag => Err(RenderTrajectoryCodecError::InvalidTag {
                field,
                tag: u64::from(tag),
            }),
        }
    }

    fn string(&mut self, field: &'static str) -> Result<String, RenderTrajectoryCodecError> {
        let length = usize::try_from(self.u32(field)?).map_err(|_| {
            RenderTrajectoryCodecError::InvalidLength {
                field,
                value: u64::MAX,
                maximum: u64::try_from(self.max_text_bytes).unwrap_or(u64::MAX),
            }
        })?;
        self.text_bytes = self.text_bytes.checked_add(length).ok_or(
            RenderTrajectoryCodecError::InvalidLength {
                field: "metadata text bytes",
                value: u64::MAX,
                maximum: u64::try_from(self.max_text_bytes).unwrap_or(u64::MAX),
            },
        )?;
        if self.text_bytes > self.max_text_bytes {
            return Err(RenderTrajectoryCodecError::InvalidLength {
                field: "metadata text bytes",
                value: u64::try_from(self.text_bytes).unwrap_or(u64::MAX),
                maximum: u64::try_from(self.max_text_bytes).unwrap_or(u64::MAX),
            });
        }
        let bytes = self.take(length, field)?;
        let value = core::str::from_utf8(bytes)
            .map_err(|_| RenderTrajectoryCodecError::InvalidUtf8(field))?;
        let mut owned = String::new();
        owned.try_reserve_exact(value.len()).map_err(|_| {
            RenderTrajectoryCodecError::Capacity {
                artifact: "metadata string",
                requested: u64::try_from(value.len()).unwrap_or(u64::MAX),
            }
        })?;
        owned.push_str(value);
        Ok(owned)
    }

    fn finish(self, field: &'static str) -> Result<(), RenderTrajectoryCodecError> {
        if self.position != self.bytes.len() {
            return Err(RenderTrajectoryCodecError::InvalidLength {
                field,
                value: u64::try_from(self.bytes.len() - self.position).unwrap_or(u64::MAX),
                maximum: 0,
            });
        }
        Ok(())
    }
}

fn decode_world_frame(
    decoder: &mut SliceDecoder<'_>,
    field: &'static str,
) -> Result<RenderWorldFrame, RenderTrajectoryCodecError> {
    match decoder.u8(field)? {
        1 => Ok(RenderWorldFrame::RightHandedZUp),
        2 => Ok(RenderWorldFrame::RightHandedYUp),
        tag => Err(RenderTrajectoryCodecError::InvalidTag {
            field,
            tag: u64::from(tag),
        }),
    }
}

fn decode_unit_system(
    decoder: &mut SliceDecoder<'_>,
    field: &'static str,
) -> Result<RenderUnitSystem, RenderTrajectoryCodecError> {
    match decoder.u8(field)? {
        1 => Ok(RenderUnitSystem::SiRadians),
        2 => Ok(RenderUnitSystem::SiDegrees),
        tag => Err(RenderTrajectoryCodecError::InvalidTag {
            field,
            tag: u64::from(tag),
        }),
    }
}

fn decode_authority(
    decoder: &mut SliceDecoder<'_>,
) -> Result<RenderTrajectoryAuthority, RenderTrajectoryCodecError> {
    match decoder.u8("metadata.authority")? {
        1 => Ok(RenderTrajectoryAuthority::SimulationEvidence),
        tag => Err(RenderTrajectoryCodecError::InvalidTag {
            field: "metadata.authority",
            tag: u64::from(tag),
        }),
    }
}

fn decode_availability_bits(
    bits: u8,
    field: &'static str,
) -> Result<RenderChannelAvailability, RenderTrajectoryCodecError> {
    if bits & !0x1f != 0 {
        return Err(RenderTrajectoryCodecError::InvalidTag {
            field,
            tag: u64::from(bits),
        });
    }
    Ok(RenderChannelAvailability {
        gravity: bits & 1 != 0,
        contact: bits & 2 != 0,
        rolling: bits & 4 != 0,
        base: bits & 8 != 0,
        gas: bits & 16 != 0,
    })
}

fn decode_rigid_state(
    decoder: &mut SliceDecoder<'_>,
    field: &'static str,
) -> Result<RigidBodyState, RenderTrajectoryCodecError> {
    let position = decoder.vec3(field)?;
    let orientation = decoder.quaternion(field)?;
    let linear_momentum = decoder.vec3(field)?;
    let angular_momentum = decoder.vec3(field)?;
    let pose = Pose::new(position, orientation)
        .map_err(|_| RenderTrajectoryCodecError::InvalidValue(field))?;
    RigidBodyState::new(pose, linear_momentum, angular_momentum)
        .map_err(|_| RenderTrajectoryCodecError::InvalidValue(field))
}

#[allow(
    clippy::too_many_lines,
    reason = "metadata decoding mirrors the frozen canonical field order"
)]
fn decode_metadata(
    bytes: &[u8],
    base_offset: u64,
    max_text_bytes: usize,
) -> Result<RenderTrajectoryMetadata, RenderTrajectoryCodecError> {
    let mut decoder = SliceDecoder::new(bytes, base_offset, max_text_bytes);
    let schema_version = decoder.u16("metadata.schema_version")?;
    let world_frame = decode_world_frame(&mut decoder, "metadata.world_frame")?;
    let units = decode_unit_system(&mut decoder, "metadata.units")?;
    let authority = decode_authority(&mut decoder)?;
    let channel_availability = decode_availability_bits(
        decoder.u8("metadata.channel_availability")?,
        "metadata.channel_availability",
    )?;
    let specimen_profile_identity = decoder.hash("metadata.specimen_profile_identity")?;
    let specimen_chart_identity = decoder.hash("metadata.specimen_chart_identity")?;
    let mass_identity = decoder.hash("metadata.mass_properties.identity")?;
    let mass = decoder.f64("metadata.mass_properties.mass")?;
    let center_of_mass_body = decoder.vec3("metadata.mass_properties.center_of_mass_body")?;
    let principal_inertia_body =
        decoder.vec3("metadata.mass_properties.principal_inertia_body")?;
    let properties = MassProperties::new(mass, center_of_mass_body, principal_inertia_body)
        .map_err(|_| RenderTrajectoryCodecError::InvalidValue("metadata.mass_properties"))?;
    let initial_state = decode_rigid_state(&mut decoder, "metadata.initial_state")?;
    let initial_base_mode = RenderBaseModeState {
        displacement_m: decoder.f64("metadata.initial_base_mode.displacement_m")?,
        velocity_m_per_s: decoder.f64("metadata.initial_base_mode.velocity_m_per_s")?,
    };
    let base_model_identity = decoder.hash("metadata.base_model_identity")?;
    let base_frame = RenderBaseFrame {
        origin_world_m: decoder.vec3("metadata.base_frame.origin_world_m")?,
        orientation_base_to_world: decoder
            .quaternion("metadata.base_frame.orientation_base_to_world")?,
    };
    let model_identity = decoder.hash("metadata.model_identity")?;
    let configuration_identity = decoder.hash("metadata.configuration_identity")?;
    let configuration_fingerprint = decoder.u64("metadata.configuration_fingerprint")?;
    let timestep_s = decoder.f64("metadata.timestep_s")?;
    let producer_version = decoder.string("metadata.producer_version")?;
    let applicability = decoder.string("metadata.applicability")?;
    let no_claim_count = usize::try_from(decoder.u32("metadata.no_claim_count")?).map_err(
        |_| RenderTrajectoryCodecError::InvalidLength {
            field: "metadata.no_claim_count",
            value: u64::MAX,
            maximum: u64::try_from(MAX_RENDER_TRAJECTORY_NO_CLAIMS).unwrap_or(u64::MAX),
        },
    )?;
    if no_claim_count > MAX_RENDER_TRAJECTORY_NO_CLAIMS {
        return Err(RenderTrajectoryCodecError::InvalidLength {
            field: "metadata.no_claim_count",
            value: u64::try_from(no_claim_count).unwrap_or(u64::MAX),
            maximum: u64::try_from(MAX_RENDER_TRAJECTORY_NO_CLAIMS).unwrap_or(u64::MAX),
        });
    }
    let mut no_claims = Vec::new();
    no_claims.try_reserve_exact(no_claim_count).map_err(|_| {
        RenderTrajectoryCodecError::Capacity {
            artifact: "metadata no-claims",
            requested: u64::try_from(no_claim_count).unwrap_or(u64::MAX),
        }
    })?;
    for _ in 0..no_claim_count {
        no_claims.push(decoder.string("metadata.no_claim")?);
    }
    decoder.finish("metadata trailing bytes")?;
    Ok(RenderTrajectoryMetadata {
        schema_version,
        world_frame,
        units,
        specimen_profile_identity,
        specimen_chart_identity,
        mass_properties: RenderMassProperties {
            identity: mass_identity,
            properties,
        },
        initial_state,
        initial_base_mode,
        base_model_identity,
        base_frame,
        model_identity,
        channel_availability,
        configuration_identity,
        configuration_fingerprint,
        timestep_s,
        producer_version,
        applicability,
        no_claims,
        authority,
    })
}

fn decode_contact_branch(
    decoder: &mut SliceDecoder<'_>,
) -> Result<RenderContactBranch, RenderTrajectoryCodecError> {
    match decoder.u8("sample.contact_branch")? {
        1 => Ok(RenderContactBranch::Open),
        2 => Ok(RenderContactBranch::Closed),
        tag => Err(RenderTrajectoryCodecError::InvalidTag {
            field: "sample.contact_branch",
            tag: u64::from(tag),
        }),
    }
}

fn decode_transition_kind(
    decoder: &mut SliceDecoder<'_>,
) -> Result<ContactTransitionKind, RenderTrajectoryCodecError> {
    match decoder.u8("sample.transition.kind")? {
        1 => Ok(ContactTransitionKind::Opening),
        2 => Ok(ContactTransitionKind::Reimpact),
        tag => Err(RenderTrajectoryCodecError::InvalidTag {
            field: "sample.transition.kind",
            tag: u64::from(tag),
        }),
    }
}

fn decode_option_tag(
    decoder: &mut SliceDecoder<'_>,
    field: &'static str,
) -> Result<bool, RenderTrajectoryCodecError> {
    match decoder.u8(field)? {
        0 => Ok(false),
        1 => Ok(true),
        tag => Err(RenderTrajectoryCodecError::InvalidTag {
            field,
            tag: u64::from(tag),
        }),
    }
}

fn decode_channels(
    decoder: &mut SliceDecoder<'_>,
) -> Result<ChannelOwnership, RenderTrajectoryCodecError> {
    fn channel(
        decoder: &mut SliceDecoder<'_>,
    ) -> Result<ChannelWrench, RenderTrajectoryCodecError> {
        Ok(ChannelWrench {
            force_world_n: decoder.vec3("sample.channel.force_world_n")?,
            torque_world_nm: decoder.vec3("sample.channel.torque_world_nm")?,
            work_j: decoder.f64("sample.channel.work_j")?,
        })
    }
    Ok(ChannelOwnership {
        gravity: channel(decoder)?,
        contact: channel(decoder)?,
        rolling: channel(decoder)?,
        base: channel(decoder)?,
        gas: channel(decoder)?,
    })
}

fn decode_disposition(
    decoder: &mut SliceDecoder<'_>,
) -> Result<RenderSampleDisposition, RenderTrajectoryCodecError> {
    let tag = decoder.u8("sample.disposition")?;
    let backend_code = decoder.u32("sample.disposition.backend_code")?;
    match (tag, backend_code) {
        (0, 0) => Ok(RenderSampleDisposition::Continue),
        (1, 0) => Ok(RenderSampleDisposition::TerminalInclination),
        (2, 0) => Ok(RenderSampleDisposition::HorizonCensored),
        (3, 0) => Ok(RenderSampleDisposition::NumericalRefusal(
            RenderNumericalRefusalReason::ReimpactLimitExceeded,
        )),
        (4, 0) => Ok(RenderSampleDisposition::NumericalRefusal(
            RenderNumericalRefusalReason::ContactEventLocalizationFailed,
        )),
        (5, 0) => Ok(RenderSampleDisposition::NumericalRefusal(
            RenderNumericalRefusalReason::NonFiniteEnergyOrBaseState,
        )),
        (6, code) => Ok(RenderSampleDisposition::NumericalRefusal(
            RenderNumericalRefusalReason::BackendSpecific(code),
        )),
        _ => Err(RenderTrajectoryCodecError::InvalidTag {
            field: "sample.disposition",
            tag: u64::from(tag),
        }),
    }
}

#[allow(
    clippy::too_many_lines,
    reason = "sample decoding mirrors the frozen canonical field order"
)]
fn decode_sample(
    bytes: &[u8],
    base_offset: u64,
) -> Result<RenderTrajectorySampleInput, RenderTrajectoryCodecError> {
    let mut decoder = SliceDecoder::new(bytes, base_offset, 0);
    let interval_start_time_s = decoder.f64("sample.interval_start_time_s")?;
    let time_s = decoder.f64("sample.time_s")?;
    let world_frame = decode_world_frame(&mut decoder, "sample.world_frame")?;
    let units = decode_unit_system(&mut decoder, "sample.units")?;
    let center_of_mass_world_m = decoder.vec3("sample.center_of_mass_world_m")?;
    let orientation_body_to_world = [
        decoder.f64("sample.orientation_body_to_world")?,
        decoder.f64("sample.orientation_body_to_world")?,
        decoder.f64("sample.orientation_body_to_world")?,
        decoder.f64("sample.orientation_body_to_world")?,
    ];
    let linear_momentum_world_kg_m_per_s =
        decoder.vec3("sample.linear_momentum_world_kg_m_per_s")?;
    let angular_momentum_body_kg_m2_per_s =
        decoder.vec3("sample.angular_momentum_body_kg_m2_per_s")?;
    let symmetry_axis_world = decoder.vec3("sample.symmetry_axis_world")?;
    let contact_branch = decode_contact_branch(&mut decoder)?;
    let contact_geometry = if decode_option_tag(&mut decoder, "sample.contact_geometry")? {
        let point_world_m = decoder.vec3("sample.contact_geometry.point_world_m")?;
        let normal_world = decoder.vec3("sample.contact_geometry.normal_world")?;
        let support_feature = match decoder.u8("sample.contact_geometry.support_feature")? {
            1 => RenderSupportFeature::CylinderRim,
            2 => {
                let raw = decoder.u64("sample.contact_geometry.profile_feature")?;
                RenderSupportFeature::ProfileFeature(usize::try_from(raw).map_err(|_| {
                    RenderTrajectoryCodecError::InvalidLength {
                        field: "sample.contact_geometry.profile_feature",
                        value: raw,
                        maximum: u64::try_from(usize::MAX).unwrap_or(u64::MAX),
                    }
                })?)
            }
            tag => {
                return Err(RenderTrajectoryCodecError::InvalidTag {
                    field: "sample.contact_geometry.support_feature",
                    tag: u64::from(tag),
                });
            }
        };
        Some(RenderContactGeometry {
            point_world_m,
            normal_world,
            support_feature,
        })
    } else {
        None
    };
    let signed_gap_m = decoder.f64("sample.signed_gap_m")?;
    let interval_contact_active = decoder.boolean("sample.interval_contact_active")?;
    let interval_normal_force_n = decoder.f64("sample.interval_normal_force_n")?;
    let transition_count = usize::try_from(decoder.u32("sample.transition_count")?).map_err(
        |_| RenderTrajectoryCodecError::InvalidLength {
            field: "sample.transition_count",
            value: u64::MAX,
            maximum: u64::try_from(MAX_RENDER_TRANSITIONS_PER_SAMPLE).unwrap_or(u64::MAX),
        },
    )?;
    if transition_count > MAX_RENDER_TRANSITIONS_PER_SAMPLE {
        return Err(RenderTrajectoryCodecError::InvalidLength {
            field: "sample.transition_count",
            value: u64::try_from(transition_count).unwrap_or(u64::MAX),
            maximum: u64::try_from(MAX_RENDER_TRANSITIONS_PER_SAMPLE).unwrap_or(u64::MAX),
        });
    }
    let mut contact_transitions = Vec::new();
    contact_transitions
        .try_reserve_exact(transition_count)
        .map_err(|_| RenderTrajectoryCodecError::Capacity {
            artifact: "sample contact transitions",
            requested: u64::try_from(transition_count).unwrap_or(u64::MAX),
        })?;
    for _ in 0..transition_count {
        contact_transitions.push(RenderContactTransition {
            kind: decode_transition_kind(&mut decoder)?,
            time_s: decoder.f64("sample.transition.time_s")?,
            bracket_start_s: decoder.f64("sample.transition.bracket_start_s")?,
            bracket_end_s: decoder.f64("sample.transition.bracket_end_s")?,
        });
    }
    let base_mode = if decode_option_tag(&mut decoder, "sample.base_mode")? {
        Some(RenderBaseModeState {
            displacement_m: decoder.f64("sample.base_mode.displacement_m")?,
            velocity_m_per_s: decoder.f64("sample.base_mode.velocity_m_per_s")?,
        })
    } else {
        None
    };
    let channels = decode_channels(&mut decoder)?;
    let mechanical_energy_j = decoder.f64("sample.mechanical_energy_j")?;
    let energy_defect_j = decoder.f64("sample.energy_defect_j")?;
    let qois = DerivedEulerQois {
        inclination_rad: decoder.f64("sample.qois.inclination_rad")?,
        precession_rad_per_s: decoder.f64("sample.qois.precession_rad_per_s")?,
        spin_rad_per_s: decoder.f64("sample.qois.spin_rad_per_s")?,
        precession_acceleration_rad_per_s2: decoder
            .f64("sample.qois.precession_acceleration_rad_per_s2")?,
    };
    let disposition = decode_disposition(&mut decoder)?;
    let terminal_event = if decode_option_tag(&mut decoder, "sample.terminal_event")? {
        Some(RenderTerminalEvent {
            time_s: decoder.f64("sample.terminal_event.time_s")?,
            bracket_start_s: decoder.f64("sample.terminal_event.bracket_start_s")?,
            bracket_end_s: decoder.f64("sample.terminal_event.bracket_end_s")?,
        })
    } else {
        None
    };
    decoder.finish("sample record trailing bytes")?;
    Ok(RenderTrajectorySampleInput {
        interval_start_time_s,
        time_s,
        world_frame,
        units,
        center_of_mass_world_m,
        orientation_body_to_world,
        linear_momentum_world_kg_m_per_s,
        angular_momentum_body_kg_m2_per_s,
        symmetry_axis_world,
        contact_branch,
        contact_geometry,
        signed_gap_m,
        interval_contact_active,
        interval_normal_force_n,
        contact_transitions,
        base_mode,
        channels,
        mechanical_energy_j,
        energy_defect_j,
        qois,
        disposition,
        terminal_event,
    })
}

fn decode_header(bytes: &[u8], base_offset: u64) -> Result<Header, RenderTrajectoryCodecError> {
    let mut decoder = SliceDecoder::new(bytes, base_offset, 0);
    if decoder.take(MAGIC.len(), "header.magic")? != MAGIC {
        return Err(RenderTrajectoryCodecError::InvalidMagic);
    }
    let codec_version = decoder.u16("header.codec_version")?;
    if codec_version != EULER_RENDER_TRAJECTORY_CODEC_VERSION {
        return Err(RenderTrajectoryCodecError::UnsupportedCodecVersion(
            codec_version,
        ));
    }
    if decoder.u16("header.trajectory_schema_version")?
        != EULER_RENDER_TRAJECTORY_SCHEMA_VERSION
    {
        return Err(RenderTrajectoryCodecError::ContractMismatch(
            "trajectory schema version",
        ));
    }
    if decoder.u16("header.control_schema_version")? != EULER_CONTROL_STREAM_SCHEMA_VERSION {
        return Err(RenderTrajectoryCodecError::ContractMismatch(
            "control schema version",
        ));
    }
    if decoder.u16("header.timeline_resampler_version")? != EULER_TIMELINE_RESAMPLER_VERSION {
        return Err(RenderTrajectoryCodecError::ContractMismatch(
            "timeline resampler version",
        ));
    }
    if decoder.u8("header.float_policy")? != FLOAT_POLICY_RAW_IEEE754_LE {
        return Err(RenderTrajectoryCodecError::ContractMismatch(
            "floating-point policy",
        ));
    }
    if decoder.u8("header.interpolation_policy")? != INTERPOLATION_CUBIC_HERMITE_SLERP_V1 {
        return Err(RenderTrajectoryCodecError::ContractMismatch(
            "interpolation policy",
        ));
    }
    if decoder.u16("header.chunking_version")? != CHUNKING_VERSION {
        return Err(RenderTrajectoryCodecError::ContractMismatch(
            "chunking version",
        ));
    }
    if usize::try_from(decoder.u32("header.samples_per_chunk")?).ok()
        != Some(EULER_RENDER_TRAJECTORY_SAMPLES_PER_CHUNK)
    {
        return Err(RenderTrajectoryCodecError::ContractMismatch(
            "samples per chunk",
        ));
    }
    let reserved = decoder.u32("header.reserved")?;
    if reserved != HEADER_RESERVED {
        return Err(RenderTrajectoryCodecError::InvalidTag {
            field: "header.reserved",
            tag: u64::from(reserved),
        });
    }
    let total_len = decoder.u64("header.total_len")?;
    let metadata_len = decoder.u32("header.metadata_len")?;
    let discontinuity_count = decoder.u32("header.discontinuity_count")?;
    let discontinuity_len = decoder.u64("header.discontinuity_len")?;
    let sample_count = decoder.u32("header.sample_count")?;
    let transition_count = decoder.u32("header.transition_count")?;
    let chunk_count = decoder.u32("header.chunk_count")?;
    let first_time_s = decoder.f64("header.first_time_s")?;
    let last_time_s = decoder.f64("header.last_time_s")?;
    let source_campaign_identity = decoder.hash("header.source_campaign_identity")?;
    let terminal_tag = decoder.u8("header.terminal_tag")?;
    let availability = decoder.u8("header.channel_availability")?;
    let world_frame = decoder.u8("header.world_frame")?;
    let units = decoder.u8("header.units")?;
    decoder.finish("header trailing bytes")?;
    Ok(Header {
        plan: WirePlan {
            total_len,
            metadata_len,
            discontinuity_count,
            discontinuity_len,
            sample_count,
            transition_count,
            chunk_count,
            first_time_s,
            last_time_s,
            terminal_tag,
        },
        source_campaign_identity,
        availability,
        world_frame,
        units,
    })
}

fn validate_header(
    header: Header,
    available: u64,
    budget: RenderTrajectoryCodecBudget,
) -> Result<(), RenderTrajectoryCodecError> {
    validate_budget(budget)?;
    let plan = header.plan;
    if plan.total_len > budget.max_artifact_bytes {
        return Err(RenderTrajectoryCodecError::ArtifactTooLarge {
            bytes: plan.total_len,
            maximum: budget.max_artifact_bytes,
        });
    }
    if plan.total_len != available {
        return Err(RenderTrajectoryCodecError::InvalidLength {
            field: "header.total_len",
            value: plan.total_len,
            maximum: available,
        });
    }
    let metadata_len = usize::try_from(plan.metadata_len).map_err(|_| {
        RenderTrajectoryCodecError::InvalidLength {
            field: "header.metadata_len",
            value: u64::from(plan.metadata_len),
            maximum: u64::try_from(MAX_METADATA_BYTES).unwrap_or(u64::MAX),
        }
    })?;
    if metadata_len > MAX_METADATA_BYTES {
        return Err(RenderTrajectoryCodecError::InvalidLength {
            field: "header.metadata_len",
            value: u64::from(plan.metadata_len),
            maximum: u64::try_from(MAX_METADATA_BYTES).unwrap_or(u64::MAX),
        });
    }
    let sample_count = usize::try_from(plan.sample_count).map_err(|_| {
        RenderTrajectoryCodecError::InvalidLength {
            field: "header.sample_count",
            value: u64::from(plan.sample_count),
            maximum: u64::try_from(budget.max_samples).unwrap_or(u64::MAX),
        }
    })?;
    if sample_count == 0 || sample_count > budget.max_samples {
        return Err(RenderTrajectoryCodecError::InvalidLength {
            field: "header.sample_count",
            value: u64::from(plan.sample_count),
            maximum: u64::try_from(budget.max_samples).unwrap_or(u64::MAX),
        });
    }
    let transition_count = usize::try_from(plan.transition_count).map_err(|_| {
        RenderTrajectoryCodecError::InvalidLength {
            field: "header.transition_count",
            value: u64::from(plan.transition_count),
            maximum: u64::try_from(budget.max_total_transitions).unwrap_or(u64::MAX),
        }
    })?;
    if transition_count > budget.max_total_transitions {
        return Err(RenderTrajectoryCodecError::InvalidLength {
            field: "header.transition_count",
            value: u64::from(plan.transition_count),
            maximum: u64::try_from(budget.max_total_transitions).unwrap_or(u64::MAX),
        });
    }
    let expected_chunks = sample_count.div_ceil(EULER_RENDER_TRAJECTORY_SAMPLES_PER_CHUNK);
    if usize::try_from(plan.chunk_count).ok() != Some(expected_chunks) {
        return Err(RenderTrajectoryCodecError::InvalidLength {
            field: "header.chunk_count",
            value: u64::from(plan.chunk_count),
            maximum: u64::try_from(expected_chunks).unwrap_or(u64::MAX),
        });
    }
    if usize::try_from(plan.discontinuity_count).ok().is_none_or(|count| count > sample_count) {
        return Err(RenderTrajectoryCodecError::InvalidLength {
            field: "header.discontinuity_count",
            value: u64::from(plan.discontinuity_count),
            maximum: u64::from(plan.sample_count),
        });
    }
    let expected_discontinuity_len = u64::from(plan.discontinuity_count)
        .checked_mul(DECLARED_DISCONTINUITY_RECORD_LEN)
        .ok_or(RenderTrajectoryCodecError::InvalidLength {
            field: "header.discontinuity_len",
            value: plan.discontinuity_len,
            maximum: u64::MAX,
        })?;
    if plan.discontinuity_len != expected_discontinuity_len {
        return Err(RenderTrajectoryCodecError::InvalidLength {
            field: "header.discontinuity_len",
            value: plan.discontinuity_len,
            maximum: expected_discontinuity_len,
        });
    }
    if !plan.first_time_s.is_finite() || !plan.last_time_s.is_finite() {
        return Err(RenderTrajectoryCodecError::InvalidValue(
            "header trajectory times",
        ));
    }
    if header
        .source_campaign_identity
        .as_bytes()
        .iter()
        .all(|byte| *byte == 0)
    {
        return Err(RenderTrajectoryCodecError::ZeroSourceCampaignIdentity);
    }
    if plan.terminal_tag > 6 {
        return Err(RenderTrajectoryCodecError::InvalidTag {
            field: "header.terminal_tag",
            tag: u64::from(plan.terminal_tag),
        });
    }
    decode_availability_bits(header.availability, "header.channel_availability")?;
    if !matches!(header.world_frame, 1 | 2) {
        return Err(RenderTrajectoryCodecError::InvalidTag {
            field: "header.world_frame",
            tag: u64::from(header.world_frame),
        });
    }
    if !matches!(header.units, 1 | 2) {
        return Err(RenderTrajectoryCodecError::InvalidTag {
            field: "header.units",
            tag: u64::from(header.units),
        });
    }

    let minimum = HEADER_LEN
        .checked_add(u64::from(plan.metadata_len))
        .and_then(|length| length.checked_add(plan.discontinuity_len))
        .and_then(|length| {
            length.checked_add(u64::from(plan.chunk_count).saturating_mul(CHUNK_HEADER_LEN))
        })
        .and_then(|length| length.checked_add(PAYLOAD_FINGERPRINT_LEN))
        .ok_or(RenderTrajectoryCodecError::InvalidLength {
            field: "minimum artifact length",
            value: u64::MAX,
            maximum: plan.total_len,
        })?;
    if plan.total_len < minimum {
        return Err(RenderTrajectoryCodecError::InvalidLength {
            field: "header.total_len",
            value: plan.total_len,
            maximum: minimum,
        });
    }
    Ok(())
}

struct IntegrityReader<'reader, R> {
    reader: &'reader mut R,
    offset: u64,
    prefix_hasher: DomainHasher,
    artifact_hasher: DomainHasher,
}

impl<'reader, R: Read> IntegrityReader<'reader, R> {
    fn new(reader: &'reader mut R, offset: u64) -> Self {
        Self {
            reader,
            offset,
            prefix_hasher: DomainHasher::new(EULER_RENDER_TRAJECTORY_PAYLOAD_FINGERPRINT_DOMAIN),
            artifact_hasher: DomainHasher::new(EULER_RENDER_TRAJECTORY_ARTIFACT_IDENTITY_DOMAIN),
        }
    }

    fn prefix(
        &mut self,
        bytes: &mut [u8],
        field: &'static str,
    ) -> Result<(), RenderTrajectoryCodecError> {
        read_exact(self.reader, bytes, field, self.offset)?;
        self.prefix_hasher.update(bytes);
        self.artifact_hasher.update(bytes);
        self.offset = checked_len_add(
            self.offset,
            u64::try_from(bytes.len()).unwrap_or(u64::MAX),
            "reader offset",
        )?;
        Ok(())
    }

    fn trailer(
        &mut self,
        bytes: &mut [u8],
        field: &'static str,
    ) -> Result<(), RenderTrajectoryCodecError> {
        read_exact(self.reader, bytes, field, self.offset)?;
        self.artifact_hasher.update(bytes);
        self.offset = checked_len_add(
            self.offset,
            u64::try_from(bytes.len()).unwrap_or(u64::MAX),
            "reader offset",
        )?;
        Ok(())
    }
}

fn read_exact(
    reader: &mut impl Read,
    bytes: &mut [u8],
    field: &'static str,
    offset: u64,
) -> Result<(), RenderTrajectoryCodecError> {
    reader.read_exact(bytes).map_err(|error| {
        if error.kind() == io::ErrorKind::UnexpectedEof {
            RenderTrajectoryCodecError::Truncated { field, offset }
        } else {
            RenderTrajectoryCodecError::Io {
                operation: "read artifact",
                kind: error.kind(),
            }
        }
    })
}

fn seek<S: Seek + ?Sized>(
    reader: &mut S,
    position: SeekFrom,
    operation: &'static str,
) -> Result<u64, RenderTrajectoryCodecError> {
    reader
        .seek(position)
        .map_err(|error| RenderTrajectoryCodecError::Io {
            operation,
            kind: error.kind(),
        })
}

fn allocate_bytes(
    length: usize,
    maximum: usize,
    artifact: &'static str,
) -> Result<Vec<u8>, RenderTrajectoryCodecError> {
    if length > maximum {
        return Err(RenderTrajectoryCodecError::InvalidLength {
            field: artifact,
            value: u64::try_from(length).unwrap_or(u64::MAX),
            maximum: u64::try_from(maximum).unwrap_or(u64::MAX),
        });
    }
    let mut bytes = Vec::new();
    bytes.try_reserve_exact(length).map_err(|_| {
        RenderTrajectoryCodecError::Capacity {
            artifact,
            requested: u64::try_from(length).unwrap_or(u64::MAX),
        }
    })?;
    bytes.resize(length, 0);
    Ok(bytes)
}

fn decode_declared_discontinuities(
    bytes: &[u8],
    base_offset: u64,
    count: usize,
) -> Result<Vec<DeclaredTimelineDiscontinuity>, RenderTrajectoryCodecError> {
    let mut decoder = SliceDecoder::new(bytes, base_offset, 0);
    let mut discontinuities = Vec::new();
    discontinuities
        .try_reserve_exact(count)
        .map_err(|_| RenderTrajectoryCodecError::Capacity {
            artifact: "declared discontinuities",
            requested: u64::try_from(count).unwrap_or(u64::MAX),
        })?;
    for _ in 0..count {
        let time_s = decoder.f64("declared_discontinuity.time_s")?;
        let kind = match decoder.u8("declared_discontinuity.kind")? {
            1 => DeclaredDiscontinuityKind::ContinuationSeam,
            2 => DeclaredDiscontinuityKind::ProducerDeclared,
            tag => {
                return Err(RenderTrajectoryCodecError::InvalidTag {
                    field: "declared_discontinuity.kind",
                    tag: u64::from(tag),
                });
            }
        };
        discontinuities.push(DeclaredTimelineDiscontinuity { time_s, kind });
    }
    decoder.finish("declared discontinuity trailing bytes")?;
    Ok(discontinuities)
}

#[derive(Debug)]
struct Preflight {
    header: Header,
    metadata: RenderTrajectoryMetadata,
    declared_discontinuities: Vec<DeclaredTimelineDiscontinuity>,
    sample_section_offset: u64,
    receipt: RenderTrajectoryCodecReceipt,
}

#[allow(
    clippy::too_many_lines,
    reason = "the bounded first pass validates the complete chunk protocol before sample retention"
)]
fn preflight<R: Read + Seek>(
    reader: &mut R,
    budget: RenderTrajectoryCodecBudget,
    start: u64,
    available: u64,
    checkpoint: &mut impl FnMut() -> Result<(), RenderTrajectoryCodecError>,
) -> Result<Preflight, RenderTrajectoryCodecError> {
    seek(reader, SeekFrom::Start(start), "seek artifact start")?;
    let mut input = IntegrityReader::new(reader, 0);
    let mut header_bytes = [0u8; HEADER_LEN as usize];
    input.prefix(&mut header_bytes, "header")?;
    let header = decode_header(&header_bytes, 0)?;
    validate_header(header, available, budget)?;
    checkpoint()?;

    let metadata_len = usize::try_from(header.plan.metadata_len).map_err(|_| {
        RenderTrajectoryCodecError::Capacity {
            artifact: "metadata",
            requested: u64::from(header.plan.metadata_len),
        }
    })?;
    let mut metadata_bytes = allocate_bytes(metadata_len, MAX_METADATA_BYTES, "metadata")?;
    let metadata_offset = input.offset;
    input.prefix(&mut metadata_bytes, "metadata")?;
    let metadata = decode_metadata(
        &metadata_bytes,
        metadata_offset,
        budget.max_total_text_bytes,
    )?;
    if availability_bits(metadata.channel_availability) != header.availability
        || world_frame_tag(metadata.world_frame) != header.world_frame
        || unit_system_tag(metadata.units) != header.units
    {
        return Err(RenderTrajectoryCodecError::ContractMismatch(
            "header metadata summary",
        ));
    }

    let discontinuity_len = usize::try_from(header.plan.discontinuity_len).map_err(|_| {
        RenderTrajectoryCodecError::Capacity {
            artifact: "declared discontinuities",
            requested: header.plan.discontinuity_len,
        }
    })?;
    let mut discontinuity_bytes = allocate_bytes(
        discontinuity_len,
        budget.max_samples.saturating_mul(9),
        "declared discontinuities",
    )?;
    let discontinuity_offset = input.offset;
    input.prefix(&mut discontinuity_bytes, "declared discontinuities")?;
    let declared_discontinuities = decode_declared_discontinuities(
        &discontinuity_bytes,
        discontinuity_offset,
        usize::try_from(header.plan.discontinuity_count).map_err(|_| {
            RenderTrajectoryCodecError::InvalidLength {
                field: "header.discontinuity_count",
                value: u64::from(header.plan.discontinuity_count),
                maximum: u64::from(header.plan.sample_count),
            }
        })?,
    )?;
    let sample_section_offset = input.offset;

    let mut observed_samples = 0usize;
    let mut observed_transitions = 0usize;
    let mut observed_first_time = None;
    let mut observed_last_time = None;
    let mut observed_terminal_tag = None;
    for expected_chunk_index in 0..header.plan.chunk_count {
        checkpoint()?;
        let descriptor_offset = input.offset;
        let mut descriptor = [0u8; 24];
        input.prefix(&mut descriptor, "chunk descriptor")?;
        let mut descriptor_decoder = SliceDecoder::new(&descriptor, descriptor_offset, 0);
        let chunk_index = descriptor_decoder.u32("chunk.index")?;
        let first_sample_index = descriptor_decoder.u32("chunk.first_sample_index")?;
        let chunk_sample_count = descriptor_decoder.u32("chunk.sample_count")?;
        let declared_chunk_transitions = descriptor_decoder.u32("chunk.transition_count")?;
        let payload_len_u64 = descriptor_decoder.u64("chunk.payload_len")?;
        descriptor_decoder.finish("chunk descriptor trailing bytes")?;
        if chunk_index != expected_chunk_index {
            return Err(RenderTrajectoryCodecError::InvalidChunk {
                chunk: chunk_index,
                field: "index",
            });
        }
        if usize::try_from(first_sample_index).ok() != Some(observed_samples) {
            return Err(RenderTrajectoryCodecError::InvalidChunk {
                chunk: chunk_index,
                field: "first_sample_index",
            });
        }
        let remaining = usize::try_from(header.plan.sample_count)
            .unwrap_or(usize::MAX)
            .saturating_sub(observed_samples);
        let expected_chunk_samples = remaining.min(EULER_RENDER_TRAJECTORY_SAMPLES_PER_CHUNK);
        if usize::try_from(chunk_sample_count).ok() != Some(expected_chunk_samples) {
            return Err(RenderTrajectoryCodecError::InvalidChunk {
                chunk: chunk_index,
                field: "sample_count",
            });
        }
        let payload_len = usize::try_from(payload_len_u64).map_err(|_| {
            RenderTrajectoryCodecError::InvalidLength {
                field: "chunk.payload_len",
                value: payload_len_u64,
                maximum: u64::try_from(MAX_CHUNK_PAYLOAD_BYTES).unwrap_or(u64::MAX),
            }
        })?;
        let mut expected_chunk_fingerprint_bytes = [0u8; 32];
        input.prefix(
            &mut expected_chunk_fingerprint_bytes,
            "chunk fingerprint",
        )?;
        let expected_chunk_fingerprint = ContentHash(expected_chunk_fingerprint_bytes);
        let payload_offset = input.offset;
        let mut payload = allocate_bytes(
            payload_len,
            MAX_CHUNK_PAYLOAD_BYTES,
            "chunk payload",
        )?;
        input.prefix(&mut payload, "chunk payload")?;
        let mut chunk_hasher = DomainHasher::new(EULER_RENDER_TRAJECTORY_CHUNK_FINGERPRINT_DOMAIN);
        chunk_hasher.update(&descriptor);
        chunk_hasher.update(&payload);
        if chunk_hasher.finalize() != expected_chunk_fingerprint {
            return Err(RenderTrajectoryCodecError::ChunkFingerprintMismatch(
                chunk_index,
            ));
        }

        let mut payload_decoder = SliceDecoder::new(&payload, payload_offset, 0);
        let mut chunk_transitions = 0usize;
        for _ in 0..expected_chunk_samples {
            let record_len_u32 = payload_decoder.u32("sample_record_len")?;
            let record_len = usize::try_from(record_len_u32).map_err(|_| {
                RenderTrajectoryCodecError::InvalidLength {
                    field: "sample_record_len",
                    value: u64::from(record_len_u32),
                    maximum: u64::try_from(MAX_SAMPLE_RECORD_BYTES).unwrap_or(u64::MAX),
                }
            })?;
            if record_len == 0 || record_len > MAX_SAMPLE_RECORD_BYTES {
                return Err(RenderTrajectoryCodecError::InvalidLength {
                    field: "sample_record_len",
                    value: u64::from(record_len_u32),
                    maximum: u64::try_from(MAX_SAMPLE_RECORD_BYTES).unwrap_or(u64::MAX),
                });
            }
            let record_offset = payload_decoder
                .base_offset
                .checked_add(u64::try_from(payload_decoder.position).unwrap_or(u64::MAX))
                .unwrap_or(u64::MAX);
            let record = payload_decoder.take(record_len, "sample record")?;
            let sample = decode_sample(record, record_offset)?;
            chunk_transitions = chunk_transitions
                .checked_add(sample.contact_transitions.len())
                .ok_or(RenderTrajectoryCodecError::InvalidLength {
                    field: "chunk transition count",
                    value: u64::MAX,
                    maximum: u64::from(declared_chunk_transitions),
                })?;
            observed_transitions = observed_transitions
                .checked_add(sample.contact_transitions.len())
                .ok_or(RenderTrajectoryCodecError::InvalidLength {
                    field: "transition count",
                    value: u64::MAX,
                    maximum: u64::from(header.plan.transition_count),
                })?;
            observed_first_time.get_or_insert(sample.time_s);
            observed_last_time = Some(sample.time_s);
            observed_terminal_tag = Some(disposition_encoding(sample.disposition).0);
            observed_samples += 1;
        }
        payload_decoder.finish("chunk payload trailing bytes")?;
        if u32::try_from(chunk_transitions).ok() != Some(declared_chunk_transitions) {
            return Err(RenderTrajectoryCodecError::InvalidChunk {
                chunk: chunk_index,
                field: "transition_count",
            });
        }
    }
    if u32::try_from(observed_samples).ok() != Some(header.plan.sample_count) {
        return Err(RenderTrajectoryCodecError::InvalidLength {
            field: "observed sample_count",
            value: u64::try_from(observed_samples).unwrap_or(u64::MAX),
            maximum: u64::from(header.plan.sample_count),
        });
    }
    if u32::try_from(observed_transitions).ok() != Some(header.plan.transition_count) {
        return Err(RenderTrajectoryCodecError::InvalidLength {
            field: "observed transition_count",
            value: u64::try_from(observed_transitions).unwrap_or(u64::MAX),
            maximum: u64::from(header.plan.transition_count),
        });
    }
    if observed_first_time.is_none_or(|time| time.to_bits() != header.plan.first_time_s.to_bits())
        || observed_last_time.is_none_or(|time| time.to_bits() != header.plan.last_time_s.to_bits())
        || observed_terminal_tag != Some(header.plan.terminal_tag)
    {
        return Err(RenderTrajectoryCodecError::ContractMismatch(
            "header trajectory summary",
        ));
    }

    let actual_payload_fingerprint = input.prefix_hasher.finalize();
    let mut trailer = [0u8; 32];
    input.trailer(&mut trailer, "payload fingerprint trailer")?;
    let expected_payload_fingerprint = ContentHash(trailer);
    if actual_payload_fingerprint != expected_payload_fingerprint {
        return Err(RenderTrajectoryCodecError::PayloadFingerprintMismatch);
    }
    if input.offset != header.plan.total_len {
        return Err(RenderTrajectoryCodecError::InvalidLength {
            field: "observed artifact length",
            value: input.offset,
            maximum: header.plan.total_len,
        });
    }
    let artifact_identity = input.artifact_hasher.finalize();
    let receipt = RenderTrajectoryCodecReceipt {
        artifact_identity,
        payload_fingerprint: actual_payload_fingerprint,
        source_campaign_identity: header.source_campaign_identity,
        byte_len: header.plan.total_len,
        sample_count: header.plan.sample_count,
        transition_count: header.plan.transition_count,
        chunk_count: header.plan.chunk_count,
    };
    Ok(Preflight {
        header,
        metadata,
        declared_discontinuities,
        sample_section_offset,
        receipt,
    })
}

fn checked_absolute(start: u64, relative: u64) -> Result<u64, RenderTrajectoryCodecError> {
    start
        .checked_add(relative)
        .ok_or(RenderTrajectoryCodecError::InvalidLength {
            field: "artifact absolute offset",
            value: u64::MAX,
            maximum: u64::MAX,
        })
}

fn decode_sample_inputs<R: Read + Seek>(
    reader: &mut R,
    start: u64,
    preflight: &Preflight,
    budget: RenderTrajectoryCodecBudget,
    checkpoint: &mut impl FnMut() -> Result<(), RenderTrajectoryCodecError>,
) -> Result<Vec<RenderTrajectorySampleInput>, RenderTrajectoryCodecError> {
    let section_start = checked_absolute(start, preflight.sample_section_offset)?;
    seek(
        reader,
        SeekFrom::Start(section_start),
        "seek sample section",
    )?;
    let sample_count = usize::try_from(preflight.header.plan.sample_count).map_err(|_| {
        RenderTrajectoryCodecError::InvalidLength {
            field: "header.sample_count",
            value: u64::from(preflight.header.plan.sample_count),
            maximum: u64::try_from(budget.max_samples).unwrap_or(u64::MAX),
        }
    })?;
    let mut inputs = Vec::new();
    inputs.try_reserve_exact(sample_count).map_err(|_| {
        RenderTrajectoryCodecError::Capacity {
            artifact: "render trajectory sample inputs",
            requested: u64::try_from(sample_count).unwrap_or(u64::MAX),
        }
    })?;
    let mut relative_offset = preflight.sample_section_offset;
    for expected_chunk_index in 0..preflight.header.plan.chunk_count {
        checkpoint()?;
        let mut descriptor = [0u8; 24];
        read_exact(reader, &mut descriptor, "chunk descriptor", relative_offset)?;
        relative_offset = checked_len_add(relative_offset, 24, "reader offset")?;
        let mut decoder = SliceDecoder::new(&descriptor, relative_offset - 24, 0);
        let chunk_index = decoder.u32("chunk.index")?;
        let first_sample_index = decoder.u32("chunk.first_sample_index")?;
        let chunk_sample_count = decoder.u32("chunk.sample_count")?;
        let _chunk_transition_count = decoder.u32("chunk.transition_count")?;
        let payload_len_u64 = decoder.u64("chunk.payload_len")?;
        decoder.finish("chunk descriptor trailing bytes")?;
        if chunk_index != expected_chunk_index
            || usize::try_from(first_sample_index).ok() != Some(inputs.len())
        {
            return Err(RenderTrajectoryCodecError::InvalidChunk {
                chunk: chunk_index,
                field: "second-pass descriptor",
            });
        }
        let expected_samples = sample_count
            .saturating_sub(inputs.len())
            .min(EULER_RENDER_TRAJECTORY_SAMPLES_PER_CHUNK);
        if usize::try_from(chunk_sample_count).ok() != Some(expected_samples) {
            return Err(RenderTrajectoryCodecError::InvalidChunk {
                chunk: chunk_index,
                field: "second-pass sample_count",
            });
        }
        let mut fingerprint = [0u8; 32];
        read_exact(
            reader,
            &mut fingerprint,
            "chunk fingerprint",
            relative_offset,
        )?;
        relative_offset = checked_len_add(relative_offset, 32, "reader offset")?;
        let payload_len = usize::try_from(payload_len_u64).map_err(|_| {
            RenderTrajectoryCodecError::InvalidLength {
                field: "chunk.payload_len",
                value: payload_len_u64,
                maximum: u64::try_from(MAX_CHUNK_PAYLOAD_BYTES).unwrap_or(u64::MAX),
            }
        })?;
        let mut payload = allocate_bytes(
            payload_len,
            MAX_CHUNK_PAYLOAD_BYTES,
            "chunk payload",
        )?;
        let payload_offset = relative_offset;
        read_exact(reader, &mut payload, "chunk payload", payload_offset)?;
        relative_offset = checked_len_add(
            relative_offset,
            u64::try_from(payload_len).unwrap_or(u64::MAX),
            "reader offset",
        )?;
        let mut payload_decoder = SliceDecoder::new(&payload, payload_offset, 0);
        for _ in 0..expected_samples {
            let record_len_u32 = payload_decoder.u32("sample_record_len")?;
            let record_len = usize::try_from(record_len_u32).map_err(|_| {
                RenderTrajectoryCodecError::InvalidLength {
                    field: "sample_record_len",
                    value: u64::from(record_len_u32),
                    maximum: u64::try_from(MAX_SAMPLE_RECORD_BYTES).unwrap_or(u64::MAX),
                }
            })?;
            if record_len == 0 || record_len > MAX_SAMPLE_RECORD_BYTES {
                return Err(RenderTrajectoryCodecError::InvalidLength {
                    field: "sample_record_len",
                    value: u64::from(record_len_u32),
                    maximum: u64::try_from(MAX_SAMPLE_RECORD_BYTES).unwrap_or(u64::MAX),
                });
            }
            let record_offset = payload_decoder
                .base_offset
                .checked_add(u64::try_from(payload_decoder.position).unwrap_or(u64::MAX))
                .unwrap_or(u64::MAX);
            let record = payload_decoder.take(record_len, "sample record")?;
            inputs.push(decode_sample(record, record_offset)?);
        }
        payload_decoder.finish("chunk payload trailing bytes")?;
    }
    if inputs.len() != sample_count {
        return Err(RenderTrajectoryCodecError::InvalidLength {
            field: "decoded sample_count",
            value: u64::try_from(inputs.len()).unwrap_or(u64::MAX),
            maximum: u64::try_from(sample_count).unwrap_or(u64::MAX),
        });
    }
    Ok(inputs)
}

struct CompareWriter<'reader, R> {
    reader: &'reader mut R,
    mismatch: bool,
    compared: u64,
}

impl<'reader, R> CompareWriter<'reader, R> {
    fn new(reader: &'reader mut R) -> Self {
        Self {
            reader,
            mismatch: false,
            compared: 0,
        }
    }
}

impl<R: Read> Write for CompareWriter<'_, R> {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        let mut compared = 0usize;
        let mut observed = [0u8; 8_192];
        while compared < bytes.len() {
            let block_len = (bytes.len() - compared).min(observed.len());
            self.reader.read_exact(&mut observed[..block_len])?;
            if observed[..block_len] != bytes[compared..compared + block_len] {
                self.mismatch = true;
            }
            compared += block_len;
        }
        self.compared = self
            .compared
            .checked_add(u64::try_from(bytes.len()).unwrap_or(u64::MAX))
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "comparison overflow"))?;
        Ok(bytes.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

fn decode_from_reader<R: Read + Seek>(
    reader: &mut R,
    budget: RenderTrajectoryCodecBudget,
    checkpoint: &mut impl FnMut() -> Result<(), RenderTrajectoryCodecError>,
) -> Result<EulerRenderTrajectoryArtifact, RenderTrajectoryCodecError> {
    validate_budget(budget)?;
    checkpoint()?;
    let start = seek(reader, SeekFrom::Current(0), "query artifact start")?;
    let end = seek(reader, SeekFrom::End(0), "query artifact end")?;
    if end < start {
        return Err(RenderTrajectoryCodecError::InvalidLength {
            field: "seekable artifact extent",
            value: end,
            maximum: start,
        });
    }
    let available = end - start;
    if available > budget.max_artifact_bytes {
        return Err(RenderTrajectoryCodecError::ArtifactTooLarge {
            bytes: available,
            maximum: budget.max_artifact_bytes,
        });
    }
    let preflight = preflight(reader, budget, start, available, checkpoint)?;
    checkpoint()?;
    let inputs = decode_sample_inputs(reader, start, &preflight, budget, checkpoint)?;
    let trajectory = RenderTrajectory::try_new(preflight.metadata.clone(), inputs)?;
    validate_context(
        preflight.header.source_campaign_identity,
        &trajectory,
        &preflight.declared_discontinuities,
    )?;

    seek(reader, SeekFrom::Start(start), "seek canonical comparison")?;
    let (canonical_receipt, mismatch, compared) = {
        let mut comparison = CompareWriter::new(reader);
        let receipt = encode_to_writer(
            &trajectory,
            preflight.header.source_campaign_identity,
            &preflight.declared_discontinuities,
            budget,
            &mut comparison,
            checkpoint,
        )?;
        (receipt, comparison.mismatch, comparison.compared)
    };
    if mismatch || compared != preflight.header.plan.total_len {
        return Err(RenderTrajectoryCodecError::NonCanonical);
    }
    if canonical_receipt != preflight.receipt {
        return Err(RenderTrajectoryCodecError::ReceiptMismatch);
    }
    checkpoint()?;
    Ok(EulerRenderTrajectoryArtifact {
        trajectory,
        source_campaign_identity: preflight.header.source_campaign_identity,
        declared_discontinuities: preflight.declared_discontinuities,
        receipt: preflight.receipt,
    })
}
