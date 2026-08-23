//! Euler-disc flagship integration boundary.
//!
//! The crate freezes the scientific Context of Use, claim taxonomy, evidence
//! ceilings, and binding no-claims. Its bounded executable rungs compose
//! geometry-derived specimen properties, ideal rolling, and a reduced flexible
//! base response. They do not establish compliant contact, gas-film, experiment,
//! or target-outcome prediction.

#![forbid(unsafe_code)]

pub mod air;
pub mod audio_artifact;
pub mod audio_excitation;
pub mod audio_resampling;
pub mod base_response;
pub mod baseline;
#[cfg(feature = "cinematic-finalization")]
pub mod cinematic_finalize;
#[cfg(feature = "cinematic-render")]
pub mod cinematic_fixture;
#[cfg(feature = "cinematic-orchestration")]
pub mod cinematic_job;
pub mod contact_dynamics;
#[cfg(feature = "scientific-contract")]
pub mod contract;
pub mod control_stream;
pub mod convergence;
pub mod coupled_runner;
pub mod external_air;
pub mod mechanics;
pub mod mechanism_registry;
pub mod modal_base_response;
pub mod modal_synthesis;
pub mod normal_contact;
pub mod patch_kinematics;
pub mod physical_audio_artifact;
pub mod ports;
pub mod production_coupling;
#[cfg(feature = "scientific-contract")]
pub mod protocol;
pub mod reduced_decay;
#[cfg(feature = "render-checkpoint-ledger")]
pub mod render_checkpoint;
pub mod render_motion_bridge;
#[cfg(feature = "cinematic-render")]
pub mod render_scene_bridge;
#[cfg(feature = "render-sharding-ledger")]
pub mod render_sharding;
pub mod render_trajectory;
pub mod render_trajectory_codec;
pub mod rolling_contact;
pub mod spatial_audio;
pub mod specimen;
pub mod structural_acoustics;
pub mod tangential_contact;
pub mod terminal_sound;
pub mod timeline_resampling;

pub use audio_artifact::{
    AUDIO_ARTIFACT_CANCELLATION_POLL_BYTES, AUDIO_ARTIFACT_CANCELLATION_POLL_FRAMES,
    AUDIO_ARTIFACT_SCHEMA_VERSION, AUDIO_TRUE_PEAK_OVERSAMPLE_FACTOR, AudioArtifactBudget,
    AudioArtifactError, AudioArtifactManifest, AudioArtifactRole, AudioChannelLayoutReceipt,
    AudioDryMixSpec, AudioMasterSource, AudioMeters, AudioSignalPath, DecodedStereoWav,
    EULER_WAV_CODEC_VERSION, MAX_AUDIO_MASTER_GAIN_DB, MAX_WAV_COMMENT_BYTES, SoundWavArtifact,
    StemGainPan, StereoSample, WavCodecReceipt, WavMetadata, WavSampleEncoding, decode_stereo_wav,
    encode_stereo_wav, measure_audio, mix_dry_modal_stems, verify_wav_against_manifest,
};
pub use audio_excitation::{
    AUDIO_EXCITATION_ALGORITHM_VERSION, AUDIO_EXCITATION_CANCELLATION_POLL_INTERVALS,
    ArtisticEventExcitation, ArtisticTextureConfig, ArtisticTextureEnvelope,
    AudioExcitationAvailability, AudioExcitationBudget, AudioExcitationCheckpoint,
    AudioExcitationChunk, AudioExcitationDiagnostics, AudioExcitationError, AudioExcitationEvent,
    AudioExcitationGrid, AudioExcitationInterval, AudioExcitationMapper, AudioExcitationModelInput,
    AudioExcitationReconstructionStatus, AudioExcitationReduction, AudioExcitationStems,
    ContactModeShape, ContactParticipationPolicy, ExcitationSourceAvailability,
    MAX_AUDIO_EXCITATION_AZIMUTHAL_HARMONIC, MAX_AUDIO_EXCITATION_CHUNK_INTERVALS,
    ModalSpatialEnvelope, ModeContactParticipationRule, SpatialEnvelopeSource,
    procedural_texture_unit_sample,
};
pub use audio_resampling::{
    AUDIO_RECONSTRUCTION_FILTER_VERSION, AUDIO_RESAMPLING_ALGORITHM_VERSION,
    AUDIO_RESAMPLING_CANCELLATION_POLL_FRAMES, AUDIO_RESAMPLING_CANCELLATION_POLL_MODES,
    AudioEventFractionalDelay, AudioReconstructionFilterDiagnostics, AudioReconstructionFilterSpec,
    AudioResampler, AudioResamplingBoundaryPolicy, AudioResamplingBudget,
    AudioResamplingCheckpoint, AudioResamplingChunk, AudioResamplingCrop,
    AudioResamplingDiagnostics, AudioResamplingError, AudioResamplingModelInput,
    AudioVideoAlignment, AudioVideoSyncMarker, EVENT_SAMPLE_SNAP_TOLERANCE_FRAMES,
    GeneralizedForceMeasureInterval, GeneralizedForceReconstructionInput,
    MAX_AUDIO_FILTER_PASSBAND_RIPPLE_DB, MAX_AUDIO_RECONSTRUCTION_FILTER_TAPS,
    MAX_SOURCE_CLOCK_ALIGNMENT_ERROR_FRAMES, MIN_AUDIO_FILTER_STOPBAND_ATTENUATION_DB,
    ReconstructedGeneralizedForce, ResampledAudioEvent, reconstruct_generalized_force_measures,
};
pub use control_stream::{
    AudioControlFilter, AudioControlInterval, AudioVisualCoverage, AudioVisualHorizon,
    AvailableChannelControl, ChannelControl, ChannelControlSet, ChannelWorkIntegralChecks,
    CoarsenedAudioBin, CoarsenedAudioControls, ContactEventMeasure, ContactFrameCoordinates,
    ControlContactEvent, ControlStreamError, EULER_CONTROL_STREAM_SCHEMA_VERSION,
    EulerControlStream, VisualizationControlPoint, WorkIntegralCheck,
};
pub use physical_audio_artifact::{
    PhysicalListeningMasterError, PhysicalPressureListeningMaster, PressureListeningMasterPolicy,
};
pub use render_trajectory::{
    DerivedEulerQois, EULER_RENDER_TRAJECTORY_SCHEMA_VERSION, MAX_RENDER_TRAJECTORY_NO_CLAIMS,
    MAX_RENDER_TRAJECTORY_SAMPLES, MAX_RENDER_TRANSITIONS_PER_SAMPLE,
    REDUCED_DECAY_RENDER_TAIL_HORIZON_S, RenderBaseFrame, RenderBaseModeState,
    RenderChannelAvailability, RenderContactBranch, RenderContactGeometry, RenderContactTransition,
    RenderMassProperties, RenderNormalForceSampling, RenderNumericalRefusalReason,
    RenderSampleDisposition, RenderSupportFeature, RenderTerminalEvent, RenderTrajectory,
    RenderTrajectoryAuthority, RenderTrajectoryError, RenderTrajectoryMetadata,
    RenderTrajectorySample, RenderTrajectorySampleInput, RenderUnitSystem, RenderWorldFrame,
};
pub use render_trajectory_codec::{
    EULER_RENDER_TRAJECTORY_ARTIFACT_IDENTITY_DOMAIN,
    EULER_RENDER_TRAJECTORY_CHUNK_FINGERPRINT_DOMAIN, EULER_RENDER_TRAJECTORY_CODEC_VERSION,
    EULER_RENDER_TRAJECTORY_PAYLOAD_FINGERPRINT_DOMAIN, EULER_RENDER_TRAJECTORY_SAMPLES_PER_CHUNK,
    EulerRenderTrajectoryArtifact, MAX_RENDER_TRAJECTORY_ARTIFACT_BYTES,
    MAX_RENDER_TRAJECTORY_TEXT_BYTES, MAX_RENDER_TRAJECTORY_TOTAL_TRANSITIONS,
    RenderTrajectoryCodecBudget, RenderTrajectoryCodecError, RenderTrajectoryCodecReceipt,
};
pub use spatial_audio::{
    ListenerPose, ListenerPoseTrack, MAX_SPATIAL_AUDIO_ROOM_IR_TAPS,
    MAX_SPATIAL_AUDIO_SAMPLE_RATE_HZ, MAX_SPATIAL_AUDIO_SOURCES, MicrophoneDirectivity,
    OfflineSpatializer, SPATIAL_AUDIO_ALGORITHM_VERSION, SPATIAL_AUDIO_CANCELLATION_POLL_FRAMES,
    SourcePositionTrack, SpatialAudioAuthority, SpatialAudioBudget, SpatialAudioConfig,
    SpatialAudioDiagnostics, SpatialAudioError, SpatialAudioOutput, SpatialAudioRenderInput,
    SpatialAudioSource, SpatialDelayPolicy, SpatialMonoSignal, SpatialOutputHorizon,
    SpatialStemComponent, StereoRoomImpulseResponse, bypass_dry_stereo,
};
pub use structural_acoustics::{
    AcousticModeRadiation, AcousticObserver, FarFieldSourceCoefficientPaM, ModalAcousticRadiation,
    ModalLossSpectrum, PhysicalContactForceSampling, PhysicalModalAudioModel,
    PhysicalModalInitialState, PhysicalModalPressureFrame, PhysicalPressureSignal,
    PointForceProjection, RETARDED_FAR_FIELD_OBSERVER_NO_CLAIM, ResolvedAcousticMedium,
    RetardedFarFieldObserverControls, STRUCTURAL_BROADBAND_SOURCE_NO_CLAIM,
    STRUCTURAL_MODAL_BASIS_SCHEMA_VERSION, STRUCTURAL_RESIDUAL_FLEXIBILITY_NO_CLAIM,
    STRUCTURAL_RESIDUAL_FLEXIBILITY_SCHEMA_VERSION, StructuralBroadbandRadiationArtifact,
    StructuralBroadbandRadiationRequest, StructuralBroadbandSourceRuntime,
    StructuralBroadbandSourceStem, StructuralMeshControls, StructuralModalBasis,
    StructuralModalBasisError, StructuralMode, StructuralModeRequest,
    StructuralResidualFlexibilityAuthority, StructuralResidualFlexibilityControls,
    StructuralResidualFlexibilityEstimateBasis, StructuralResidualFlexibilityEstimateComparison,
    StructuralResidualFlexibilityEstimateResponse, StructuralResidualModalLossSpectrum,
    build_structural_broadband_radiation_artifact, build_structural_modal_basis,
    build_structural_residual_flexibility_estimate_basis,
    compare_structural_residual_flexibility_estimates, modal_loss_spectrum_from_prony,
    modal_loss_spectrum_from_rayleigh, superpose_pressure_signals,
    synthesize_retarded_far_field_world_observers,
};
pub use timeline_resampling::{
    DeclaredDiscontinuityKind, DeclaredTimelineDiscontinuity, EULER_TIMELINE_RESAMPLER_VERSION,
    EventEvaluationSide, ExposureEventPolicy, ExposurePartition, ExposureSegment,
    ResampledTimelineSample, TimelineEvent, TimelineInterpolationMethod, TimelineResampler,
    TimelineResamplingError, TimelineSampleSource,
};

pub use base_response::{
    BaseGeometryScope, BaseResponseDiagnostics, BaseResponseError, BaseResponseInput,
    BaseResponseRefinement, BaseResponseRun, BaseResponseSample, ContactLoadScope,
    LevelSupportInput, MAX_BASE_RESPONSE_STEPS, MovingContactLoad, refine_reduced_base_response,
    run_reduced_base_response,
};
pub use baseline::{
    BaselineDynamicsClass, BaselineEnergyLedger, BaselineEquilibriumReceipt, BaselineRefusal,
    BaselineRefusalReason, BaselineRunOutput, BaselineSample, BaselineState,
    BaselineSupportDiagnostic, BaselineTerminal, BaselineTrajectory, SquatDiscInput,
    run_ideal_conservative_baseline,
};
pub use contact_dynamics::{
    ContactDynamicsError, ContactDynamicsInput, ContactDynamicsRun, ContactGeometry,
    ContactStepReceipt, ContactTermination, DiscGeometry as ContactDiscGeometry, EnergyLedger,
    NO_CLAIM_BOUNDARY as CONTACT_NO_CLAIM_BOUNDARY, ProfileContactDynamicsInput,
    ProfileContactGeometry, ProfileContactPatchGeometry, ProfileRollingInitializer,
    StickFeasibility, TimestepRefinement, contact_geometry, declared_profile_rolling_initializer,
    profile_contact_geometry, profile_contact_patch_geometry, profile_state_at_ground_contact,
    refine_profile_timestep_by_two, refine_timestep_by_two, run_contact_dynamics,
    run_profile_contact_dynamics, small_angle_rolling_profile_initializer, state_at_ground_contact,
    state_at_profile_ground_contact,
};
pub use modal_synthesis::{
    MAX_MODAL_SPATIAL_PARTICIPATION, MAX_MODAL_SYNTHESIS_CHUNK_FRAMES,
    MAX_MODAL_SYNTHESIS_TOTAL_FRAMES, MODAL_CANCELLATION_POLL_FRAMES,
    MODAL_SYNTHESIS_ALGORITHM_VERSION, ModalComponentValues, ModalCouplingClass, ModalDriveFrame,
    ModalModeEnergy, ModalModeState, ModalPresetAuthority, ModalSpatialParticipation,
    ModalStemFrame, ModalSynthesisBudget, ModalSynthesisCheckpoint, ModalSynthesisChunk,
    ModalSynthesisDiagnostics, ModalSynthesisError, ModalSynthesisModel, ModalSynthesisModelInput,
    RepresentativeDiscMaterial, RepresentativeModalPreset, representative_modal_preset,
};
pub use ports::{
    ChannelActivity, ContributionDomain, ContributionOwnership, DecompositionReceipt,
    EULER_COMPOSITION_PORT_SCHEMA_VERSION, EnergyClosureDisposition, EnergyContribution,
    EnergyLedgerCheckpoint, EnergyTerms, EulerChannel, EulerEnergyLedger, EulerPortError,
    EulerPortRegistry, GeneralizedVelocityCoordinate, MAX_EULER_PORT_DECLARATIONS, PatchRegion,
    PortDeclaration, PortInterval, SurfacePair,
};

#[cfg(feature = "scientific-contract")]
pub use contract::{
    AuthorityCeiling, CLAIM_POLICY_ASSESSMENT_LOG_DOMAIN, CONTRACT_CHECK_RECEIPT_DOMAIN,
    CORE_NO_CLAIMS, ContractError, ContractIdentity, EULER_ASSESSMENT_IDENTITY_DOMAIN,
    EULER_CLAIM_GRAPH_IDENTITY_DOMAIN, EULER_CLAIM_POLICY_SCHEMA_VERSION,
    EULER_CONTRACT_IDENTITY_DOMAIN, EULER_CONTRACT_SCHEMA_VERSION,
    EULER_EVIDENCE_PACKET_IDENTITY_DOMAIN, EULER_OWNER_MATRIX_IDENTITY_DOMAIN,
    EULER_OWNER_MATRIX_SCHEMA_VERSION, EULER_PREREQUISITE_RECEIPT_IDENTITY_DOMAIN,
    EulerAcceptanceFamily, EulerClaimGraph, EulerClaimKind, EulerClaimSpec, EulerContextExtension,
    EulerScientificContract, EvidenceRequirement, FS_GOVERN_AUTHORITY_SOURCE_SCHEMA,
    FS_IR_CAMPAIGN_SOURCE_SCHEMA, HYPOTHESIS_SOURCE_DECLARATION_DOMAIN, HypothesisSource,
    MAX_EULER_CLAIMS, MAX_EULER_NO_CLAIMS, MAX_OWNER_MATRIX_BYTES, OwnerMatrix,
    OwnerMatrixIdentity, OwnerRole, OwnerRow, ScientificRisk,
};
#[cfg(feature = "scientific-contract")]
pub use protocol::{
    AssessmentDisposition, ClaimEvidencePacket, ClaimPolicyAssessment, ClaimPolicyAssessmentLog,
    ContractCheckReceipt, DeclaredEvidenceAccessClass, EULER_PROTOCOL_SCHEMA_VERSION,
    EvidenceAuthorityClass, EvidenceAuthorityDeclaration, EvidenceRecord,
    FROZEN_CLAIM_GRAPH_HASH_HEX, FROZEN_CONTEXT_HASH_HEX, FROZEN_CONTRACT_IDENTITY_HEX,
    MAX_PREREQUISITE_RECEIPTS, MAX_PROTOCOL_ID_BYTES, MAX_VALIDITY_DOMAIN_AXES,
    MAX_VALIDITY_DOMAIN_CANONICAL_BYTES, PrerequisiteAssessmentReceipt, ProtocolBudget,
    ProtocolSeed, ReportedScientificDisposition, StructurallyAdmittedEulerContract,
    admit_frozen_contract, build_frozen_contract, check_frozen_contract,
};
pub use reduced_decay::{
    BILDSTEN_PUBLISHED_POWER_COEFFICIENT, BildstenBoundaryLayerChannel,
    CHANNEL_CROSSOVER_BISECTION_STEPS, ChannelCrossoverDiagnostic, ChannelCrossoverNotComparable,
    ChannelPowers, ChannelWork, DryContourChannel, MAX_REDUCED_DECAY_STEPS,
    MAX_SMALL_ANGLE_THETA_RAD, REDUCED_DECAY_MODEL_ID, ReducedDecayError, ReducedDecayInput,
    ReducedDecayProvenance, ReducedDecayRun, ReducedDecaySample, ReducedDecayTerminal,
    RefinementEvidence, STANDARD_GRAVITY_M_PER_S2, channel_crossover_diagnostic,
    refinement_evidence, run_reduced_decay, structured_runner_output,
};
