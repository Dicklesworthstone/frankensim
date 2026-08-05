//! Euler-disc flagship integration boundary.
//!
//! The crate freezes the scientific Context of Use, claim taxonomy, evidence
//! ceilings, and binding no-claims. Its bounded executable rungs compose
//! geometry-derived specimen properties, ideal rolling, and a reduced flexible
//! base response. They do not establish compliant contact, gas-film, experiment,
//! or target-outcome prediction.

#![forbid(unsafe_code)]

pub mod air;
pub mod base_response;
pub mod baseline;
pub mod contact_dynamics;
#[cfg(feature = "scientific-contract")]
pub mod contract;
pub mod control_stream;
pub mod convergence;
pub mod coupled_runner;
pub mod external_air;
pub mod mechanics;
pub mod normal_contact;
pub mod patch_kinematics;
pub mod ports;
pub mod production_coupling;
#[cfg(feature = "scientific-contract")]
pub mod protocol;
pub mod reduced_decay;
pub mod render_motion_bridge;
pub mod render_trajectory;
pub mod rolling_contact;
pub mod specimen;
pub mod tangential_contact;
pub mod timeline_resampling;

pub use control_stream::{
    AudioControlFilter, AudioControlInterval, AvailableChannelControl, ChannelControl,
    ChannelControlSet, ChannelWorkIntegralChecks, CoarsenedAudioBin, CoarsenedAudioControls,
    ContactEventMeasure, ContactFrameCoordinates, ControlContactEvent, ControlStreamError,
    EULER_CONTROL_STREAM_SCHEMA_VERSION, EulerControlStream, VisualizationControlPoint,
    WorkIntegralCheck,
};
pub use render_trajectory::{
    DerivedEulerQois, EULER_RENDER_TRAJECTORY_SCHEMA_VERSION, MAX_RENDER_TRAJECTORY_NO_CLAIMS,
    MAX_RENDER_TRAJECTORY_SAMPLES, MAX_RENDER_TRANSITIONS_PER_SAMPLE, RenderBaseFrame,
    RenderBaseModeState, RenderChannelAvailability, RenderContactBranch, RenderContactGeometry,
    RenderContactTransition, RenderMassProperties, RenderNumericalRefusalReason,
    RenderSampleDisposition, RenderSupportFeature, RenderTerminalEvent, RenderTrajectory,
    RenderTrajectoryAuthority, RenderTrajectoryError, RenderTrajectoryMetadata,
    RenderTrajectorySample, RenderTrajectorySampleInput, RenderUnitSystem, RenderWorldFrame,
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
    ProfileContactGeometry, ProfileRollingInitializer, StickFeasibility, TimestepRefinement,
    contact_geometry, profile_contact_geometry, profile_state_at_ground_contact,
    refine_profile_timestep_by_two, refine_timestep_by_two, run_contact_dynamics,
    run_profile_contact_dynamics, small_angle_rolling_profile_initializer, state_at_ground_contact,
    state_at_profile_ground_contact,
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
