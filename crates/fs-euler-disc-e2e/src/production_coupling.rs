//! One atomic, mechanically homogeneous Euler-disc accepted substep.
//!
//! This is an orchestration boundary over existing adapters.  It is deliberately
//! Estimate-only: it does not calibrate a disc, rank air correlations, resolve
//! impact, supply thin-gap pressure, or claim a resolved/as-built base.

use core::fmt::{self, Write as _};

use fs_blake3::{ContentHash, DomainHasher};
use fs_contact::normal_patch::{NormalPatchEmbedState, NormalPatchPort, NormalPatchReceipt};
use fs_couple::{StableId, modal_acoustic_time::ModalAcousticState};
use fs_exec::Cx;
use fs_mbd::{
    Gravity, MassProperties, Pose, RigidBodyIntegrator, RigidBodyState, StepReceipt, Vec3, Wrench,
};
use fs_rep_frep::{
    AxisymmetricCurvatureAuthority, AxisymmetricError, AxisymmetricIdentity,
    AxisymmetricMassProperties, AxisymmetricQuerySession, AxisymmetricSupportError,
};
use fs_tribo::{
    InputAuthority, InterfaceSystemRef,
    partial_slip::{NormalPatchAuthority, NormalPatchView},
    rolling_loss::{PatchCurvature, RollingPatchReceipt},
    surface_excitation::{
        AdmittedSurfaceTracePair, FilteredSurfacePairReceipt, HertzRoughnessExcitationInput,
        HertzRoughnessExcitationReceipt, ProjectedHertzFootprint, SurfaceExcitationError,
        SurfaceTraceMotion, UniformSurfaceTrace, evaluate_hertz_filtered_surface_pair,
        evaluate_hertz_roughness_excitation, evaluate_point_surface_pair,
    },
};

use crate::{
    air::{
        AirFilmError, AirFilmProposal, AirFilmTransactionState, AirVec3, TiltedDiscAirFilmInput,
        TiltedDiscKinematics,
    },
    base_response::{
        BaseResponseError, ReducedBaseCheckpoint, ReducedBasePort, ReducedBaseStepInput,
        ReducedBaseStepProposal,
    },
    contact_dynamics::{
        ContactDynamicsError, ProfileContactGeometry, ProfileContactPatchGeometry,
        profile_contact_geometry_from_query_session, profile_contact_patch_geometry_from_mass,
        profile_contact_patch_geometry_from_query_session, profile_mass_to_mbd,
    },
    external_air::{
        EulerDiscBodyFrame, EulerDiscExteriorState, EulerExternalAirCandidate,
        EulerExternalAirInput, EulerExternalAirWorkProposal, EulerExternalAirWorkState,
        ExternalAirError, evaluate_euler_disc_external_air,
    },
    modal_base_response::{
        RectangularModalBaseCheckpoint, RectangularModalBaseError, RectangularModalBasePort,
        RectangularModalBaseProposal, RectangularModalBaseStepInput,
    },
    normal_contact::{
        ActiveNormalContact, EulerNormalContactInput, EulerNormalContactOutcome,
        NormalContactError, NormalContactIntegrationRegime, evaluate_normal_contact,
    },
    patch_kinematics::{
        CurvatureMetadata, MovingOneModePatchKinematicsInput, PatchContactStatus, PatchKinematics,
        PatchKinematicsError, ProfileSupportKinematics, compute_moving_one_mode_patch_kinematics,
    },
    rolling_contact::{
        RollingContactError, RollingContactInput, RollingContactProposal, RollingContactState,
        commit_rolling_contact, prepare_rolling_contact,
    },
    specimen::ResolvedDiscProfile,
    tangential_contact::{
        EulerTangentialContactAdapter, TangentialContactError, TangentialContactReceipt,
        TangentialContactRequest, TangentialContactState,
    },
};

/// The one gas mechanism that owns a production substep's aerodynamic work.
///
/// This sum type prevents exterior free-gas and thin-gap film wrenches from
/// being supplied together for the same interval.
#[derive(Debug, Clone, PartialEq)]
pub enum GasChannelStepInput {
    /// Free exterior gas; one caller-named correlation is selected without ranking.
    ExteriorFreeGas {
        /// Exterior-air cards and candidate alternatives.
        input: EulerExternalAirInput,
        /// Exact selected correlation identity.
        selected_correlation_id: String,
        /// Exactly-once exterior work key.
        exchange_key: u64,
    },
    /// Thin-gap gas-film card, whose transactional state owns restart and work.
    ThinGap {
        /// Film cards; disc pose/rates and duration are replaced from fs-mbd.
        input: TiltedDiscAirFilmInput,
        /// Exactly-once wall-to-gas work key.
        exchange_key: u64,
    },
}

/// Accepted state for exactly one mutually exclusive gas mechanism.
#[derive(Debug, Clone, PartialEq)]
pub enum GasChannelState {
    /// Exterior-air exact-once work ledger.
    ExteriorFreeGas(EulerExternalAirWorkState),
    /// Thin-gap gas-film checkpoint and exact-once work ledger.
    ThinGap(AirFilmTransactionState),
}

/// Accepted receipt for the gas mechanism selected by a substep.
#[derive(Debug, Clone, PartialEq)]
pub enum GasChannelReceipt {
    /// One explicitly named exterior-air candidate and its staged work.
    ExteriorFreeGas {
        /// Selected candidate; alternatives were not ranked or averaged.
        candidate: EulerExternalAirCandidate,
        /// Exactly-once exterior work proposal.
        work: EulerExternalAirWorkProposal,
    },
    /// One staged thin-gap candidate with its restart/work transaction.
    ThinGap {
        /// Gas-film proposal; its receipt wrench is gas-on-disc in world axes.
        proposal: AirFilmProposal,
    },
}

/// Stable owner identity for an accepted production-coupling lineage.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ProductionCouplingIdentity {
    /// Case/run identity.
    pub case_id: String,
    /// Configuration identity binding every checkpoint in this composition.
    pub configuration_id: String,
    /// Common inertial frame identity.
    pub world_frame_id: String,
}

/// Selected structural support backend for the coupled mechanics owner.
#[derive(Debug, Clone)]
pub enum ProductionBasePort {
    /// Legacy load-shaped one-mode flat-plate reduction.
    ReducedOneMode(ReducedBasePort),
    /// Resolved multi-mode rectangular plate with moving point contact.
    RectangularModal(RectangularModalBasePort),
}

impl From<ReducedBasePort> for ProductionBasePort {
    fn from(value: ReducedBasePort) -> Self {
        Self::ReducedOneMode(value)
    }
}

impl From<RectangularModalBasePort> for ProductionBasePort {
    fn from(value: RectangularModalBasePort) -> Self {
        Self::RectangularModal(value)
    }
}

/// Accepted support state paired with the selected production backend.
#[derive(Debug, Clone, PartialEq)]
pub enum ProductionBaseCheckpoint {
    /// Legacy one-mode checkpoint.
    ReducedOneMode(ReducedBaseCheckpoint),
    /// Resolved rectangular modal checkpoint.
    RectangularModal(RectangularModalBaseCheckpoint),
}

impl ProductionBaseCheckpoint {
    fn accepted_version(&self) -> u64 {
        match self {
            Self::ReducedOneMode(state) => state.accepted_version(),
            Self::RectangularModal(state) => state.accepted_version(),
        }
    }

    fn elapsed_time_s(&self) -> f64 {
        match self {
            Self::ReducedOneMode(state) => state.elapsed_time_s(),
            Self::RectangularModal(state) => state.elapsed_time_s(),
        }
    }

    fn last_surface_state(&self) -> (f64, f64) {
        match self {
            Self::ReducedOneMode(state) => {
                (state.modal_displacement_m(), state.modal_velocity_m_per_s())
            }
            Self::RectangularModal(state) => {
                let surface = state.last_surface_state();
                (surface.displacement_m, surface.velocity_m_per_s)
            }
        }
    }

    fn last_contact_point_world_m(&self) -> Vec3 {
        match self {
            Self::ReducedOneMode(_) => Vec3::ZERO,
            Self::RectangularModal(state) => {
                let [x, y, z] = state.last_contact_point_base_m();
                Vec3::new(x, y, z)
            }
        }
    }
}

/// Backend-independent support accounting consumed by render and sound bridges.
#[derive(Debug, Clone, PartialEq)]
pub struct ProductionBaseStepReceipt {
    /// Consumed support version.
    pub parent_version: u64,
    /// Produced support version.
    pub next_version: u64,
    /// Accepted subinterval [s].
    pub timestep_s: f64,
    /// Compressive force applied into the support [N].
    pub compressive_normal_force_on_base_n: f64,
    /// Equal-and-opposite reaction on the disc [N].
    pub normal_reaction_on_disc_world_n: [f64; 3],
    /// Local contact displacement at interval start [m].
    pub modal_displacement_start_m: f64,
    /// Local contact displacement at interval end [m].
    pub modal_displacement_end_m: f64,
    /// Local contact velocity at interval start [m/s].
    pub modal_velocity_start_m_per_s: f64,
    /// Local contact velocity at interval end [m/s].
    pub modal_velocity_end_m_per_s: f64,
    /// Retained structural-energy change [J].
    pub stored_energy_change_j: f64,
    /// Viscous damping work [J].
    pub damping_work_j: f64,
    /// External contact work on the support [J].
    pub external_contact_work_j: f64,
    /// Structural energy-closure residual [J].
    pub energy_closure_residual_j: f64,
    /// Available norm of the support reaction [N].
    pub end_support_reaction_norm_n: f64,
}

#[derive(Debug, Clone, PartialEq)]
enum ProductionBaseProposalBackend {
    ReducedOneMode(ReducedBaseStepProposal),
    RectangularModal(RectangularModalBaseProposal),
}

/// Prepared support transition whose accepted state remains private until commit.
#[derive(Debug, Clone, PartialEq)]
pub struct ProductionBaseStepProposal {
    backend: ProductionBaseProposalBackend,
    receipt: ProductionBaseStepReceipt,
}

impl ProductionBaseStepProposal {
    /// Backend-independent interval accounting.
    #[must_use]
    pub const fn receipt(&self) -> &ProductionBaseStepReceipt {
        &self.receipt
    }

    fn rectangular_modal_audio_parts(&self) -> Option<(&[f64], &[ModalAcousticState])> {
        match &self.backend {
            ProductionBaseProposalBackend::RectangularModal(proposal) => {
                Some((proposal.modal_force_n_per_sqrt_kg(), proposal.next_states()))
            }
            ProductionBaseProposalBackend::ReducedOneMode(_) => None,
        }
    }
}

/// Immutable adapters and rigid-body properties used by a production substep.
#[derive(Debug, Clone)]
pub struct ProductionCouplingModel {
    /// Identity binding checkpoints to this model.
    pub identity: ProductionCouplingIdentity,
    /// Disc mass and principal inertia.
    pub disc_mass_properties: MassProperties,
    /// Uniform world-frame gravity.
    pub gravity: Gravity,
    /// Already assembled structural support port.
    pub base_port: ProductionBasePort,
    /// Explicitly selected partial-slip adapter/lane.
    pub tangential_adapter: EulerTangentialContactAdapter,
}

/// Run-scoped proof that one immutable profile and mechanics model agree.
///
/// Construction revalidates the public `ResolvedDiscProfile` and retains the
/// freshly integrated properties plus one admitted chart-query session.
/// Dynamic support and curvature remain live, cancellable local queries, but
/// they do not repeat immutable whole-chart construction validation.
pub(crate) struct AdmittedAxisymmetricProfile<'run> {
    model: &'run ProductionCouplingModel,
    profile: &'run ResolvedDiscProfile,
    mass_properties: AxisymmetricMassProperties,
    query_session: AxisymmetricQuerySession<'run>,
}

#[derive(Clone, Copy)]
enum CheckpointValidation {
    Required,
    /// The enclosing serial driver validated its start and exclusively owns
    /// every successor created by the commit paths below.
    TrustedInternalSuccessor,
}

#[derive(Clone, Copy)]
enum SurfaceTraceEvaluation<'a> {
    Checked,
    Admitted(&'a AdmittedSurfaceTracePair),
}

impl SurfaceTraceEvaluation<'_> {
    fn evaluate_point_surface_pair(
        self,
        interface: &InterfaceSystemRef,
        surface_a: SurfaceTraceMotion<'_>,
        surface_b: SurfaceTraceMotion<'_>,
    ) -> Result<FilteredSurfacePairReceipt, SurfaceExcitationError> {
        match self {
            Self::Checked => evaluate_point_surface_pair(interface, surface_a, surface_b),
            Self::Admitted(traces) => {
                traces.evaluate_point_surface_pair(interface, surface_a, surface_b)
            }
        }
    }

    fn evaluate_hertz_filtered_surface_pair(
        self,
        interface: &InterfaceSystemRef,
        surface_a: SurfaceTraceMotion<'_>,
        surface_b: SurfaceTraceMotion<'_>,
        footprint: ProjectedHertzFootprint,
    ) -> Result<FilteredSurfacePairReceipt, SurfaceExcitationError> {
        match self {
            Self::Checked => {
                evaluate_hertz_filtered_surface_pair(interface, surface_a, surface_b, footprint)
            }
            Self::Admitted(traces) => traces
                .evaluate_hertz_filtered_surface_pair(interface, surface_a, surface_b, footprint),
        }
    }
}

/// Cloneable accepted state across all participating adapters.
#[derive(Debug, Clone, PartialEq)]
pub struct ProductionCouplingCheckpoint {
    identity: ProductionCouplingIdentity,
    /// Monotone outer accepted-substep version.
    pub committed_version: u64,
    /// Disc rigid-body state.
    pub disc_state: RigidBodyState,
    /// Deterministic integrity binding for the complete accepted snapshot.
    checkpoint_fingerprint: ContentHash,
    normal_state: NormalPatchEmbedState,
    tangential_state: TangentialContactState,
    rolling_state: RollingContactState,
    gas_channel_state: GasChannelState,
    base_state: ProductionBaseCheckpoint,
}

impl ProductionCouplingCheckpoint {
    /// Stable owner identity for this accepted lineage.
    #[must_use]
    pub const fn identity(&self) -> &ProductionCouplingIdentity {
        &self.identity
    }

    /// Content hash binding identity, version, rigid state, and every channel snapshot.
    #[must_use]
    pub const fn fingerprint(&self) -> ContentHash {
        self.checkpoint_fingerprint
    }

    /// Accepted physical time retained by the base port [s].
    #[must_use]
    pub fn elapsed_time_s(&self) -> f64 {
        self.base_state.elapsed_time_s()
    }

    /// Accepted one-mode base displacement [m].
    #[must_use]
    pub fn base_displacement_m(&self) -> f64 {
        self.base_state.last_surface_state().0
    }

    /// Accepted one-mode base velocity [m/s].
    #[must_use]
    pub fn base_velocity_m_per_s(&self) -> f64 {
        self.base_state.last_surface_state().1
    }

    /// Number of accepted finite-contact intervals in this eventful lineage.
    ///
    /// Open-flight steps advance the outer checkpoint and gas/base channels,
    /// but intentionally retain every dry-contact state. Contact work keys
    /// must therefore be indexed by this counter rather than by
    /// [`Self::committed_version`].
    #[must_use]
    pub const fn committed_contact_intervals(&self) -> u64 {
        self.tangential_state.committed_version()
    }

    /// Complete accepted tangential checkpoint, including reversible history.
    #[must_use]
    pub const fn tangential_state(&self) -> &TangentialContactState {
        &self.tangential_state
    }
}

/// All caller-owned cards and one explicit selection for one attempted substep.
#[derive(Debug, Clone, PartialEq)]
pub struct ProductionCouplingStepInput {
    /// Outer checkpoint version observed by the caller before preparing cards.
    pub expected_checkpoint_version: u64,
    /// Positive event-localized duration [s].
    pub duration_s: f64,
    /// Sample time used by the normal smooth fixed-branch state [s].
    pub time_s: f64,
    /// Actual support direction/profile feature/relative-gap curvature plus the moving-base bridge.
    pub patch: MovingOneModePatchKinematicsInput,
    /// Normal material/law card; kinematics, state, time, and duration are replaced by this owner.
    pub normal: EulerNormalContactInput,
    /// Optional measured or explicitly declared small-amplitude topography.
    ///
    /// When present, the actual accepted Hertz footprint and normal tangent
    /// filter these traces before the perturbation is applied to both bodies.
    /// This is a first-order contact linearization, not a sound preset.
    pub surface_excitation: Option<ProductionSurfaceExcitationStepInput>,
    /// Optional nonlinear surface geometry resolved inside the contact law.
    ///
    /// This and `surface_excitation` are mutually exclusive: applying both
    /// would count the same height once in the gap and again as a force.
    pub surface_geometry: Option<ProductionSurfaceGeometryStepInput>,
    /// Tangential card/ownership; kinematics, normal patch, version, and duration are replaced.
    pub tangential: TangentialContactRequest,
    /// Rolling card/ownership; state/checkpoint/patch/interval are replaced.
    pub rolling: RollingContactInput,
    /// Exactly one selected gas mechanism; the enum forbids exterior/film double counting.
    pub gas_channel: GasChannelStepInput,
    /// Base-port replay identity for this accepted interval.
    pub base_step_id: String,
    /// Moving-load location at interval start.
    pub base_load_progress_start: f64,
    /// Moving-load location at interval end.
    pub base_load_progress_end: f64,
}

/// Caller-owned inputs for one mechanically open interval.
///
/// Contact, rolling, and tangential laws are intentionally absent. The disc
/// advances under gravity plus exactly one selected gas mechanism while the
/// support advances with zero applied contact load. This is the reusable
/// counterpart to [`ProductionCouplingStepInput`] for separation intervals;
/// it does not infer an impact impulse or coefficient of restitution.
#[derive(Debug, Clone, PartialEq)]
pub struct ProductionOpenFlightStepInput {
    /// Outer checkpoint version observed by the caller.
    pub expected_checkpoint_version: u64,
    /// Positive event-localized duration [s].
    pub duration_s: f64,
    /// Exactly one selected gas mechanism.
    pub gas_channel: GasChannelStepInput,
    /// Base-port replay identity for this zero-contact-load interval.
    pub base_step_id: String,
    /// Moving-load coordinate retained for deterministic support state lineage.
    pub base_load_progress_start: f64,
    /// Moving-load coordinate retained for deterministic support state lineage.
    pub base_load_progress_end: f64,
}

/// One owned surface trace and its current material-frame contact coordinate.
#[derive(Debug, Clone, PartialEq)]
pub struct ProductionSurfaceTraceStepInput {
    /// Measured or explicitly declared spatial height trace.
    pub trace: UniformSurfaceTrace,
    /// Current contact-path coordinate in the trace's material frame [m].
    pub path_coordinate_m: f64,
    /// Coordinate rate in that same material frame [m/s].
    pub path_speed_m_per_s: f64,
}

impl ProductionSurfaceTraceStepInput {
    fn as_motion(&self) -> SurfaceTraceMotion<'_> {
        SurfaceTraceMotion {
            trace: &self.trace,
            path_coordinate_m: self.path_coordinate_m,
            path_speed_m_per_s: self.path_speed_m_per_s,
        }
    }
}

/// Small-amplitude topography channel for one production contact substep.
///
/// The two surfaces retain the exact ordering of `interface`. Material names,
/// audible frequencies, and renderer choices are absent: this channel contains
/// only interface identity, spatial geometry, contact-path kinematics, and the
/// declared linearization domain.
#[derive(Debug, Clone, PartialEq)]
pub struct ProductionSurfaceExcitationStepInput {
    /// Ordered dry-interface identity used by the normal-contact law.
    pub interface: InterfaceSystemRef,
    /// First ordered surface trace and motion.
    pub surface_a: ProductionSurfaceTraceStepInput,
    /// Second ordered surface trace and motion.
    pub surface_b: ProductionSurfaceTraceStepInput,
    /// Travel direction measured from the accepted patch major axis [rad].
    pub travel_angle_from_patch_major_rad: f64,
    /// Maximum admitted absolute filtered-height/nominal-approach ratio.
    pub maximum_linearized_height_fraction: f64,
}

/// Nonlinear surface-geometry coupling for one production contact substep.
///
/// Unlike [`ProductionSurfaceExcitationStepInput`], this mode changes the
/// unilateral gap before the normal law is evaluated. The accepted Hertz
/// footprint then filters the traces and the contact is re-resolved until the
/// filtered height is self-consistent. It is the appropriate rung at first
/// touch, re-contact, or whenever a tangent perturbation about a smooth-contact
/// approach would be singular.
#[derive(Debug, Clone, PartialEq)]
pub struct ProductionSurfaceGeometryStepInput {
    /// Ordered dry-interface identity used by the normal-contact law.
    pub interface: InterfaceSystemRef,
    /// First ordered surface trace and motion.
    pub surface_a: ProductionSurfaceTraceStepInput,
    /// Second ordered surface trace and motion.
    pub surface_b: ProductionSurfaceTraceStepInput,
    /// Travel direction measured from the accepted patch major axis [rad].
    pub travel_angle_from_patch_major_rad: f64,
    /// Positive fixed-point iteration ceiling.
    pub maximum_iterations: usize,
    /// Absolute convergence tolerance for filtered height [m].
    pub absolute_height_tolerance_m: f64,
    /// Absolute convergence tolerance for filtered height rate [m/s].
    pub absolute_height_rate_tolerance_m_per_s: f64,
    /// Relative convergence tolerance for filtered height and height rate.
    pub relative_tolerance: f64,
}

/// Accepted self-consistent surface geometry used by one normal-contact solve.
#[derive(Debug, Clone, PartialEq)]
pub struct ProductionSurfaceGeometryReceipt {
    /// Final Hertz-filtered ordered surface-pair geometry.
    pub filtered_pair: FilteredSurfacePairReceipt,
    /// Number of normal/footprint fixed-point evaluations performed.
    pub iterations: usize,
    /// Final absolute change in combined filtered height [m].
    pub height_residual_m: f64,
    /// Final absolute change in combined height material derivative [m/s].
    pub height_rate_residual_m_per_s: f64,
    /// Smooth nominal signed gap before surface geometry [m].
    pub nominal_signed_gap_m: f64,
    /// Signed gap actually supplied to the accepted normal law [m].
    pub resolved_signed_gap_m: f64,
    /// Smooth-contact force at the same nominal bulk geometry [N].
    pub nominal_smooth_force_n: f64,
    /// Surface-resolved normal force sent to mechanics and the support [N].
    pub resolved_force_n: f64,
}

struct ResolvedProductionContact {
    patch_kinematics: PatchKinematics,
    normal: ActiveNormalContact,
    surface_geometry: Option<ProductionSurfaceGeometryReceipt>,
}

struct PreparedProductionContact {
    patch_kinematics: PatchKinematics,
    normal: ActiveNormalContact,
    surface_excitation: Option<HertzRoughnessExcitationReceipt>,
    surface_geometry: Option<ProductionSurfaceGeometryReceipt>,
    tangential: TangentialContactReceipt,
    rolling: RollingContactProposal,
    gas_channel: GasChannelReceipt,
    normal_force_n: f64,
    applied_tangential_force_world_n: Vec3,
    applied_tangential_free_torsional_torque_world_nm: Vec3,
    total_force_world_n: Vec3,
    total_moment_about_com_world_n_m: Vec3,
}

enum ResolvedProductionEvent {
    CompliantContact(ResolvedProductionContact),
    OpenFlight(PatchKinematics),
}

enum ProductionMidpointInput<'a> {
    Borrowed(&'a ProductionCouplingStepInput),
    Reusable(&'a mut ProductionCouplingStepInput),
}

impl ProductionMidpointInput<'_> {
    fn as_ref(&self) -> &ProductionCouplingStepInput {
        match self {
            Self::Borrowed(input) => input,
            Self::Reusable(input) => input,
        }
    }
}

impl ResolvedProductionEvent {
    const fn patch_kinematics(&self) -> &PatchKinematics {
        match self {
            Self::CompliantContact(contact) => &contact.patch_kinematics,
            Self::OpenFlight(patch) => patch,
        }
    }

    const fn branch(&self) -> ProductionTrajectoryBranch {
        match self {
            Self::CompliantContact(_) => ProductionTrajectoryBranch::CompliantContact,
            Self::OpenFlight(_) => ProductionTrajectoryBranch::OpenFlight,
        }
    }
}

/// Binds one accepted rigid/profile state to the horizontal-plane production patch.
///
/// This is the profile-native bridge into [`ProductionCouplingStepInput`]. It
/// resolves the support arm and smooth principal curvatures from the same
/// body-frame ground direction, replaces the stale rigid-body/mass fields, and
/// aligns the one-mode base counterpart beneath that support point. The caller
/// still owns material cards and must select an [`crate::normal_contact::EulerNormalGeometry`]
/// compatible with the returned principal-curvature pair. No smoothing radius
/// or desired experimental ranking is inferred here.
///
/// `base_vertical_displacement_m` and `base_vertical_velocity_m_per_s` must be
/// the values from the accepted base checkpoint that will be passed to
/// [`ProductionCouplingModel::step`]. The step owner replaces them from its
/// private checkpoint before evaluation; passing the same accepted values here
/// keeps classification and the reconstructed moving counterpart coherent.
fn bind_horizontal_plane_axisymmetric_profile_contact_input(
    input: &mut ProductionCouplingStepInput,
    resolved: ProfileContactPatchGeometry,
    disc_mass_properties: MassProperties,
    disc_state: RigidBodyState,
    base_vertical_displacement_m: f64,
    base_vertical_velocity_m_per_s: f64,
) -> Result<ProfileContactPatchGeometry, ProductionCouplingError> {
    if !base_vertical_displacement_m.is_finite() || !base_vertical_velocity_m_per_s.is_finite() {
        return Err(ProductionCouplingError::InvalidInput {
            field: "profile base mode",
        });
    }
    let profile_support = &resolved.support;
    let mut support_kinematics =
        ProfileSupportKinematics::from_profile_contact_geometry(*profile_support);
    support_kinematics.gap_m -= base_vertical_displacement_m;
    if !support_kinematics.gap_m.is_finite() {
        return Err(ProductionCouplingError::InvalidInput {
            field: "profile moving-base gap",
        });
    }
    let curvature_identity = StableId::new(format!(
        "axisymmetric/{:016x}/feature-{}/principal-curvature-v1/{:016x}-{:016x}-{:016x}",
        resolved.profile_identity.0,
        resolved.curvature.source_feature,
        resolved.curvature.meridional_m_inverse.to_bits(),
        resolved.curvature.azimuthal_m_inverse.to_bits(),
        resolved.curvature.uncertainty_m_inverse.to_bits(),
    ))
    .map_err(|_| ProductionCouplingError::InvalidInput {
        field: "profile curvature identity",
    })?;
    input.patch.bridge.profile_support = support_kinematics;
    input.patch.bridge.disc_state = disc_state;
    input.patch.bridge.disc_mass_properties = disc_mass_properties;
    input
        .patch
        .bridge
        .base_mode
        .undeformed_contact_point_world_m = Vec3::new(
        profile_support.contact.point_world_m.x,
        profile_support.contact.point_world_m.y,
        0.0,
    );
    input.patch.bridge.base_mode.vertical_displacement_m = base_vertical_displacement_m;
    input.patch.bridge.base_mode.vertical_velocity_m_per_s = base_vertical_velocity_m_per_s;
    input.patch.bridge.normal_world = Vec3::new(0.0, 0.0, 1.0);
    input.patch.patch.source_feature = profile_support.support_source_feature;
    input.patch.patch.curvature = CurvatureMetadata::Known {
        curvature_identity,
        authority: match resolved.curvature.authority {
            AxisymmetricCurvatureAuthority::Estimate => InputAuthority::Estimated,
        },
        first_principal_m_inverse: resolved.curvature.meridional_m_inverse,
        second_principal_m_inverse: resolved.curvature.azimuthal_m_inverse,
        uncertainty_m_inverse: resolved.curvature.uncertainty_m_inverse,
    };
    input.rolling.contact_arm_world_m = profile_support.contact.radius_world_m;
    Ok(resolved)
}

/// Published mechanical result and retained Estimate-only boundary.
#[derive(Debug, Clone, PartialEq)]
pub struct ProductionCouplingReceipt {
    /// Full constructed patch kinematics.
    pub patch_kinematics: PatchKinematics,
    /// Active finite-patch normal response.
    pub normal: ActiveNormalContact,
    /// Optional finite-patch-filtered topography force actually included in
    /// the accepted rigid-body and support wrenches.
    pub surface_excitation: Option<HertzRoughnessExcitationReceipt>,
    /// Optional self-consistent surface geometry already included in the
    /// accepted normal response, rigid-body wrench, and support load.
    pub surface_geometry: Option<ProductionSurfaceGeometryReceipt>,
    /// Prepared and accepted tangential response.
    pub tangential: TangentialContactReceipt,
    /// Tangential force actually included in the rigid-body midpoint wrench [N].
    pub applied_tangential_force_world_n: Vec3,
    /// Free tangential torque actually included in the midpoint wrench [N m].
    pub applied_tangential_free_torsional_torque_world_nm: Vec3,
    /// Prepared and accepted rolling response.
    pub rolling: RollingContactProposal,
    /// The one gas-channel receipt contributing its real wrench to fs-mbd.
    pub gas_channel: GasChannelReceipt,
    /// Accepted selected-base transition accounting.
    pub base: ProductionBaseStepProposal,
    /// Total real world-frame force sent to fs-mbd [N].
    pub total_force_world_n: Vec3,
    /// Total real world-frame moment about disc COM sent to fs-mbd [N m].
    pub total_moment_about_com_world_n_m: Vec3,
    /// fs-mbd accepted disc state.
    pub next_disc_state: RigidBodyState,
    /// Complete rigid-body before/after state and energy diagnostics from the
    /// exact step that produced `next_disc_state`.
    pub rigid_step: StepReceipt,
    /// This composition retains only source-adapter Estimate authority.
    pub estimate_only: bool,
}

/// Accepted open-flight transition over the same shared checkpoint as contact.
#[derive(Debug, Clone, PartialEq)]
pub struct ProductionOpenFlightReceipt {
    /// The one gas-channel receipt contributing its real wrench to fs-mbd.
    pub gas_channel: GasChannelReceipt,
    /// Zero-contact-load support transition; the support may continue ringing.
    pub base: ProductionBaseStepProposal,
    /// Gas force sent to fs-mbd [N]; gravity remains integrator-owned.
    pub total_force_world_n: Vec3,
    /// Gas moment about disc COM sent to fs-mbd [N m].
    pub total_moment_about_com_world_n_m: Vec3,
    /// fs-mbd accepted disc state.
    pub next_disc_state: RigidBodyState,
    /// Complete rigid-body before/after and energy diagnostics.
    pub rigid_step: StepReceipt,
    /// This composition retains only source-adapter Estimate authority.
    pub estimate_only: bool,
}

impl ProductionCouplingReceipt {
    /// Filter two named surface traces through this accepted Hertz patch.
    ///
    /// The normal law supplies the approach, force, consistent tangent, and
    /// circular/elliptic footprint. The caller still supplies material-frame
    /// path coordinates and speeds because the contact solver cannot infer a
    /// texture frame from bulk kinematics. The ordered interface identity must
    /// exactly match the one used by the accepted normal-contact law.
    pub fn evaluate_surface_excitation(
        &self,
        interface: &InterfaceSystemRef,
        surface_a: SurfaceTraceMotion<'_>,
        surface_b: SurfaceTraceMotion<'_>,
        travel_angle_from_patch_major_rad: f64,
        maximum_linearized_height_fraction: f64,
    ) -> Result<HertzRoughnessExcitationReceipt, ProductionSurfaceExcitationError> {
        evaluate_surface_excitation_for_normal(
            &self.normal,
            interface,
            surface_a,
            surface_b,
            travel_angle_from_patch_major_rad,
            maximum_linearized_height_fraction,
        )
    }
}

fn evaluate_surface_excitation_for_normal(
    active: &ActiveNormalContact,
    interface: &InterfaceSystemRef,
    surface_a: SurfaceTraceMotion<'_>,
    surface_b: SurfaceTraceMotion<'_>,
    travel_angle_from_patch_major_rad: f64,
    maximum_linearized_height_fraction: f64,
) -> Result<HertzRoughnessExcitationReceipt, ProductionSurfaceExcitationError> {
    let NormalPatchReceipt::Point(normal) = &active.generic.receipt else {
        return Err(ProductionSurfaceExcitationError::UnsupportedLineContact);
    };
    if interface.ordered_system_id() != normal.interface_system_id {
        return Err(
            ProductionSurfaceExcitationError::InterfaceIdentityMismatch {
                accepted_normal_interface: normal.interface_system_id.clone(),
                supplied_interface: interface.ordered_system_id().to_owned(),
            },
        );
    }
    let (semi_major_axis_m, semi_minor_axis_m) = normal
        .elliptic_patch_axes
        .map_or((normal.patch_radius_m, normal.patch_radius_m), |axes| {
            (axes.semi_major_axis_m, axes.semi_minor_axis_m)
        });
    evaluate_hertz_roughness_excitation(HertzRoughnessExcitationInput {
        interface,
        surface_a,
        surface_b,
        footprint: ProjectedHertzFootprint {
            semi_major_axis_m,
            semi_minor_axis_m,
            travel_angle_from_major_rad: travel_angle_from_patch_major_rad,
        },
        nominal_approach_m: normal.approach_m,
        nominal_normal_force_n: normal.normal_force_n,
        normal_tangent_n_per_m: normal.tangent_n_per_m,
        maximum_linearized_height_fraction,
    })
    .map_err(ProductionSurfaceExcitationError::Surface)
}

/// Refusal while adapting one accepted production contact to surface excitation.
#[derive(Debug, Clone, PartialEq)]
pub enum ProductionSurfaceExcitationError {
    /// The production composition does not admit line-contact audio forcing.
    UnsupportedLineContact,
    /// A trace/interface card was not the ordered system used by normal contact.
    InterfaceIdentityMismatch {
        /// Ordered interface retained by the accepted normal response.
        accepted_normal_interface: String,
        /// Ordered interface offered with the surface traces.
        supplied_interface: String,
    },
    /// The reusable fs-tribo filtering/linearization leaf refused the query.
    Surface(SurfaceExcitationError),
    /// The self-consistent surface-height/contact-footprint solve exhausted
    /// its caller-declared iteration budget without meeting both tolerances.
    SurfaceGeometryDidNotConverge {
        /// Number of fixed-point evaluations attempted.
        iterations: usize,
        /// Final absolute filtered-height change [m].
        height_residual_m: f64,
        /// Final absolute filtered-height-rate change [m/s].
        height_rate_residual_m_per_s: f64,
    },
}

impl fmt::Display for ProductionSurfaceExcitationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for ProductionSurfaceExcitationError {}

fn validate_surface_geometry_identity(
    normal: &EulerNormalContactInput,
    surface: &ProductionSurfaceGeometryStepInput,
) -> Result<(), ProductionCouplingError> {
    if surface.interface.ordered_system_id() != normal.material.interface.ordered_system_id() {
        return Err(ProductionCouplingError::SurfaceExcitation(
            ProductionSurfaceExcitationError::InterfaceIdentityMismatch {
                accepted_normal_interface: normal.material.interface.ordered_system_id().to_owned(),
                supplied_interface: surface.interface.ordered_system_id().to_owned(),
            },
        ));
    }
    if surface.maximum_iterations == 0
        || surface.maximum_iterations > 128
        || !(surface.absolute_height_tolerance_m.is_finite()
            && surface.absolute_height_tolerance_m > 0.0)
        || !(surface.absolute_height_rate_tolerance_m_per_s.is_finite()
            && surface.absolute_height_rate_tolerance_m_per_s > 0.0)
        || !(surface.relative_tolerance.is_finite()
            && surface.relative_tolerance > 0.0
            && surface.relative_tolerance <= 1.0)
    {
        return Err(ProductionCouplingError::InvalidInput {
            field: "surface geometry solve controls",
        });
    }
    Ok(())
}

fn compute_surface_adjusted_patch(
    mut patch_input: MovingOneModePatchKinematicsInput,
    surface: &FilteredSurfacePairReceipt,
) -> Result<PatchKinematics, ProductionCouplingError> {
    patch_input.bridge.base_mode.vertical_displacement_m += surface.combined_effective_height_m;
    patch_input.bridge.base_mode.vertical_velocity_m_per_s +=
        surface.combined_effective_height_rate_m_per_s;
    if !(patch_input
        .bridge
        .base_mode
        .vertical_displacement_m
        .is_finite()
        && patch_input
            .bridge
            .base_mode
            .vertical_velocity_m_per_s
            .is_finite())
    {
        return Err(ProductionCouplingError::InvalidInput {
            field: "surface-adjusted base geometry",
        });
    }
    compute_moving_one_mode_patch_kinematics(patch_input).map_err(ProductionCouplingError::Patch)
}

fn signed_patch_gap_m(patch: &PatchKinematics) -> Result<f64, ProductionCouplingError> {
    let gap_m = patch
        .disc_point
        .point_world
        .sub(patch.base_point.point_world)
        .dot(patch.tangent_basis.normal_world);
    if gap_m.is_finite() {
        Ok(gap_m)
    } else {
        Err(ProductionCouplingError::InvalidInput {
            field: "surface geometry signed gap",
        })
    }
}

fn surface_footprint(
    normal: &ActiveNormalContact,
    travel_angle_from_patch_major_rad: f64,
) -> Result<ProjectedHertzFootprint, ProductionCouplingError> {
    let NormalPatchReceipt::Point(point) = &normal.generic.receipt else {
        return Err(ProductionCouplingError::UnsupportedLineNormalContact);
    };
    let (semi_major_axis_m, semi_minor_axis_m) = point
        .elliptic_patch_axes
        .map_or((point.patch_radius_m, point.patch_radius_m), |axes| {
            (axes.semi_major_axis_m, axes.semi_minor_axis_m)
        });
    Ok(ProjectedHertzFootprint {
        semi_major_axis_m,
        semi_minor_axis_m,
        travel_angle_from_major_rad: travel_angle_from_patch_major_rad,
    })
}

fn resolve_surface_geometry_contact(
    patch_input: MovingOneModePatchKinematicsInput,
    normal_input: &EulerNormalContactInput,
    surface: &ProductionSurfaceGeometryStepInput,
    surface_traces: SurfaceTraceEvaluation<'_>,
) -> Result<
    (
        PatchKinematics,
        ActiveNormalContact,
        ProductionSurfaceGeometryReceipt,
    ),
    ProductionCouplingError,
> {
    validate_surface_geometry_identity(normal_input, surface)?;
    let nominal_patch = compute_moving_one_mode_patch_kinematics(patch_input.clone())
        .map_err(ProductionCouplingError::Patch)?;
    let nominal_signed_gap_m = signed_patch_gap_m(&nominal_patch)?;
    let mut nominal_normal = normal_input.clone();
    nominal_normal.kinematics = nominal_patch;
    let (nominal_smooth_force_n, nominal_footprint) =
        match evaluate_normal_contact(&nominal_normal).map_err(ProductionCouplingError::Normal)? {
            EulerNormalContactOutcome::Active(active) => (
                point_normal_force(&active)?,
                Some(surface_footprint(
                    &active,
                    surface.travel_angle_from_patch_major_rad,
                )?),
            ),
            EulerNormalContactOutcome::InactiveSeparated { .. } => (0.0, None),
        };
    let filtered_for_footprint = |footprint: Option<ProjectedHertzFootprint>| {
        match footprint {
            Some(footprint)
                if footprint.semi_major_axis_m > 0.0 && footprint.semi_minor_axis_m > 0.0 =>
            {
                surface_traces.evaluate_hertz_filtered_surface_pair(
                    &surface.interface,
                    surface.surface_a.as_motion(),
                    surface.surface_b.as_motion(),
                    footprint,
                )
            }
            _ => surface_traces.evaluate_point_surface_pair(
                &surface.interface,
                surface.surface_a.as_motion(),
                surface.surface_b.as_motion(),
            ),
        }
        .map_err(ProductionSurfaceExcitationError::Surface)
        .map_err(ProductionCouplingError::SurfaceExcitation)
    };
    // Start on the finite-contact branch already selected by the nominal bulk
    // geometry.  A point-height seed is still the correct grazing limit, but
    // using it inside an established Hertz patch can converge to a different
    // roughness/load branch at shallow Euler-disc inclinations.
    let mut filtered = filtered_for_footprint(nominal_footprint)?;
    let mut last_height_residual_m = f64::INFINITY;
    let mut last_height_rate_residual_m_per_s = f64::INFINITY;

    for iteration in 1..=surface.maximum_iterations {
        let patch = compute_surface_adjusted_patch(patch_input.clone(), &filtered)?;
        let mut trial_normal = normal_input.clone();
        trial_normal.kinematics = patch.clone();
        let active = match evaluate_normal_contact(&trial_normal)
            .map_err(ProductionCouplingError::Normal)?
        {
            EulerNormalContactOutcome::Active(active) => active,
            EulerNormalContactOutcome::InactiveSeparated { .. } => {
                return Err(ProductionCouplingError::UnsupportedMechanism {
                    status: PatchContactStatus::Separated,
                });
            }
        };
        let footprint = surface_footprint(&active, surface.travel_angle_from_patch_major_rad)?;
        // Static-load bracketing legitimately evaluates the exact grazing
        // state before a finite Hertz patch exists.  Preserve the point-limit
        // geometry at that endpoint; clamping in an artificial patch radius
        // would change both the load root and the filtered topography.
        let next_filtered = filtered_for_footprint(Some(footprint))?;
        last_height_residual_m = (next_filtered.combined_effective_height_m
            - filtered.combined_effective_height_m)
            .abs();
        last_height_rate_residual_m_per_s = (next_filtered.combined_effective_height_rate_m_per_s
            - filtered.combined_effective_height_rate_m_per_s)
            .abs();
        let height_limit_m = surface.absolute_height_tolerance_m
            + surface.relative_tolerance
                * next_filtered
                    .combined_effective_height_m
                    .abs()
                    .max(filtered.combined_effective_height_m.abs());
        let rate_limit_m_per_s = surface.absolute_height_rate_tolerance_m_per_s
            + surface.relative_tolerance
                * next_filtered
                    .combined_effective_height_rate_m_per_s
                    .abs()
                    .max(filtered.combined_effective_height_rate_m_per_s.abs());
        if last_height_residual_m <= height_limit_m
            && last_height_rate_residual_m_per_s <= rate_limit_m_per_s
        {
            // The receipt and the returned mechanical state must describe the
            // same accepted filtered geometry.  Returning `active` here would
            // expose the preceding fixed-point iterate while naming
            // `next_filtered` in the receipt.
            let resolved_patch =
                compute_surface_adjusted_patch(patch_input.clone(), &next_filtered)?;
            let mut resolved_normal = normal_input.clone();
            resolved_normal.kinematics = resolved_patch.clone();
            let resolved_active = match evaluate_normal_contact(&resolved_normal)
                .map_err(ProductionCouplingError::Normal)?
            {
                EulerNormalContactOutcome::Active(active) => active,
                EulerNormalContactOutcome::InactiveSeparated { .. } => {
                    return Err(ProductionCouplingError::UnsupportedMechanism {
                        status: PatchContactStatus::Separated,
                    });
                }
            };
            let resolved_force_n = point_normal_force(&resolved_active)?;
            return Ok((
                resolved_patch.clone(),
                resolved_active,
                ProductionSurfaceGeometryReceipt {
                    filtered_pair: next_filtered,
                    iterations: iteration,
                    height_residual_m: last_height_residual_m,
                    height_rate_residual_m_per_s: last_height_rate_residual_m_per_s,
                    nominal_signed_gap_m,
                    resolved_signed_gap_m: signed_patch_gap_m(&resolved_patch)?,
                    nominal_smooth_force_n,
                    resolved_force_n,
                },
            ));
        }
        filtered = next_filtered;
    }
    Err(ProductionCouplingError::SurfaceExcitation(
        ProductionSurfaceExcitationError::SurfaceGeometryDidNotConverge {
            iterations: surface.maximum_iterations,
            height_residual_m: last_height_residual_m,
            height_rate_residual_m_per_s: last_height_rate_residual_m_per_s,
        },
    ))
}

/// Terminal state of a bounded smooth-contact trajectory attempt.
///
/// A refusal is an expected, typed handoff boundary: the failed proposed
/// substep is not included in [`SmoothContactTrajectory::accepted_steps`] and
/// [`SmoothContactTrajectory::last_accepted_checkpoint`] remains restartable.
#[derive(Debug, Clone, PartialEq)]
pub enum SmoothContactTrajectoryTermination {
    /// The caller-selected accepted-step budget was consumed.
    StepLimitReached {
        /// Maximum accepted steps requested for this invocation.
        maximum_accepted_steps: usize,
    },
    /// The next homogeneous smooth-contact step was refused without commit.
    Refused {
        /// Outer checkpoint version from which the failed step was proposed.
        attempted_checkpoint_version: u64,
        /// Exact lower-level refusal; impact/event/gas-film owners remain separate.
        error: ProductionCouplingError,
    },
}

/// Replayable prefix produced by [`ProductionCouplingModel::run_smooth_contact_trajectory`].
///
/// The input factory is caller-supplied and must be deterministic for an exact
/// replay claim. It is called from each accepted checkpoint so support
/// direction, resolved patch metadata, and all cards can be recomputed after
/// the preceding mechanics update rather than frozen at the initial state.
#[derive(Debug, Clone, PartialEq)]
pub struct SmoothContactTrajectory {
    /// Checkpoint supplied to this invocation; also the prefix-resume boundary.
    pub start_checkpoint: ProductionCouplingCheckpoint,
    /// The last checkpoint actually accepted by all channels.
    pub last_accepted_checkpoint: ProductionCouplingCheckpoint,
    /// One receipt for every accepted mechanics substep, in order.
    pub accepted_steps: Vec<ProductionCouplingReceipt>,
    /// Bounded completion or the first uncommitted typed refusal.
    pub termination: SmoothContactTrajectoryTermination,
}

/// Mechanically homogeneous branch selected from pre-constitutive gap/rate
/// kinematics at an accepted checkpoint.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProductionTrajectoryBranch {
    /// Finite-patch normal/tangential/rolling contact is active.
    CompliantContact,
    /// No contact port exists; only gravity, gas, and unforced support dynamics advance.
    OpenFlight,
}

/// One accepted eventful trajectory interval.
#[derive(Debug, Clone, PartialEq)]
pub enum ProductionTrajectoryStepReceipt {
    /// Contact interval evaluated by the complete production composition.
    CompliantContact(ProductionCouplingReceipt),
    /// Open interval evaluated without fabricated contact channels.
    OpenFlight(ProductionOpenFlightReceipt),
}

/// A branch transition bracketed by the fixed accepted time grid.
///
/// No exact event time is claimed: the transition occurred after
/// `bracket_start_s` and no later than `bracket_end_s`. Refinement of the
/// caller-selected timestep must shrink this bracket for a timing claim.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProductionBranchTransition {
    /// Branch used by the preceding accepted interval.
    pub from: ProductionTrajectoryBranch,
    /// Branch selected at the current accepted checkpoint.
    pub to: ProductionTrajectoryBranch,
    /// Start of the unresolved transition bracket [s].
    pub bracket_start_s: f64,
    /// End of the unresolved transition bracket [s].
    pub bracket_end_s: f64,
}

/// One accepted interval with its shared before/after clock and branch.
#[derive(Debug, Clone, PartialEq)]
pub struct ProductionTrajectoryStep {
    /// Accepted interval start [s].
    pub start_time_s: f64,
    /// Accepted interval end [s].
    pub end_time_s: f64,
    /// Homogeneous mechanics branch for the interval.
    pub branch: ProductionTrajectoryBranch,
    /// Surface-resolved signed gap that selected `branch` at the interval
    /// start [m]. This includes nonlinear footprint-filtered topography.
    pub resolved_signed_gap_start_m: f64,
    /// Exact branch-specific receipt.
    pub receipt: ProductionTrajectoryStepReceipt,
}

/// Borrowed plate-acoustic data from one actually accepted mechanics step.
#[derive(Debug, Clone, Copy)]
pub struct ProductionModalAudioStep<'a> {
    /// Accepted interval start [s].
    pub start_time_s: f64,
    /// Accepted interval end [s].
    pub end_time_s: f64,
    /// Held generalized force for every retained plate mode.
    pub modal_force_n_per_sqrt_kg: &'a [f64],
    /// Actual accepted modal state at the closing boundary.
    pub accepted_states: &'a [ModalAcousticState],
    /// Exact rigid-body before/after states for the same accepted interval.
    pub rigid_step: &'a StepReceipt,
}

impl ProductionTrajectoryStep {
    fn modal_audio_step(&self) -> Option<ProductionModalAudioStep<'_>> {
        let (base, rigid_step) = match &self.receipt {
            ProductionTrajectoryStepReceipt::CompliantContact(receipt) => {
                (&receipt.base, &receipt.rigid_step)
            }
            ProductionTrajectoryStepReceipt::OpenFlight(receipt) => {
                (&receipt.base, &receipt.rigid_step)
            }
        };
        let (modal_force_n_per_sqrt_kg, accepted_states) = base.rectangular_modal_audio_parts()?;
        Some(ProductionModalAudioStep {
            start_time_s: self.start_time_s,
            end_time_s: self.end_time_s,
            modal_force_n_per_sqrt_kg,
            accepted_states,
            rigid_step,
        })
    }
}

/// Terminal state of a bounded event-aware compliant trajectory attempt.
#[derive(Debug, Clone, PartialEq)]
pub enum ProductionEventTrajectoryTermination {
    /// Caller-selected accepted-step budget was consumed.
    StepLimitReached {
        /// Maximum accepted intervals requested by the caller.
        maximum_accepted_steps: usize,
    },
    /// The next candidate refused atomically; the prior checkpoint remains restartable.
    Refused {
        /// Version of the checkpoint at which classification or stepping refused.
        attempted_checkpoint_version: u64,
        /// Exact typed refusal from the constituent or branch admission.
        error: ProductionCouplingError,
    },
}

/// Replayable open/contact trajectory over one shared coupled checkpoint.
#[derive(Debug, Clone, PartialEq)]
pub struct ProductionEventTrajectory {
    /// Checkpoint supplied to this invocation.
    pub start_checkpoint: ProductionCouplingCheckpoint,
    /// Last state accepted by every active channel.
    pub last_accepted_checkpoint: ProductionCouplingCheckpoint,
    /// Ordered accepted homogeneous intervals.
    pub accepted_steps: Vec<ProductionTrajectoryStep>,
    /// Fixed-grid brackets for every observed branch change.
    pub transitions: Vec<ProductionBranchTransition>,
    /// Bounded completion or first uncommitted refusal.
    pub termination: ProductionEventTrajectoryTermination,
}

/// One bounded-memory mechanics control interval reduced from homogeneous
/// accepted substeps.
///
/// Force is reduced as an exact impulse and reported as its interval mean;
/// endpoint state and base motion remain the actual last accepted values. A
/// branch change always flushes the preceding interval, so contact and open
/// work are never averaged together.
#[derive(Debug, Clone, PartialEq)]
pub struct ProductionControlInterval {
    /// First accepted mechanics time in this homogeneous group [s].
    pub start_time_s: f64,
    /// Last accepted mechanics time in this homogeneous group [s].
    pub end_time_s: f64,
    /// Number of accepted mechanics substeps reduced into this record.
    pub mechanics_substeps: usize,
    /// Homogeneous contact/open branch for the complete interval.
    pub branch: ProductionTrajectoryBranch,
    /// Surface-resolved signed gap that selected the branch at interval start
    /// [m]. This is not the nominal smooth-profile gap.
    pub resolved_signed_gap_start_m: f64,
    /// Surface-resolved signed gap at the accepted endpoint [m], when endpoint
    /// classification succeeded. A refused final prefix may legitimately lack
    /// this value and cannot mint a resolved-gap render/audio sample.
    pub resolved_signed_gap_end_m: Option<f64>,
    /// Exact rigid state before the first reduced substep.
    pub state_before: RigidBodyState,
    /// Exact rigid state after the last reduced substep.
    pub state_after: RigidBodyState,
    /// Exact normal-force impulse over the reduced interval [N s].
    pub normal_impulse_n_s: f64,
    /// Impulse divided by exact interval duration [N].
    pub mean_normal_force_n: f64,
    /// Base displacement at interval start [m].
    pub base_displacement_start_m: f64,
    /// Base displacement at interval end [m].
    pub base_displacement_end_m: f64,
    /// Base velocity at interval start [m/s].
    pub base_velocity_start_m_per_s: f64,
    /// Base velocity at interval end [m/s].
    pub base_velocity_end_m_per_s: f64,
    /// Accepted base version before the first reduced substep.
    pub base_parent_version: u64,
    /// Accepted base version after the last reduced substep.
    pub base_next_version: u64,
    /// Sum of exact discrete disc wrench-work residuals [J].
    pub disc_work_residual_j: f64,
    /// Accepted disc mechanical energy at the interval endpoint [J].
    pub mechanical_energy_end_j: f64,
}

/// Bounded-memory event trajectory suitable for audio/render control streams.
#[derive(Debug, Clone, PartialEq)]
pub struct ProductionControlTrajectory {
    /// Initial complete coupled checkpoint.
    pub start_checkpoint: ProductionCouplingCheckpoint,
    /// Last complete checkpoint accepted by every constituent.
    pub last_accepted_checkpoint: ProductionCouplingCheckpoint,
    /// Homogeneous reduced control intervals.
    pub intervals: Vec<ProductionControlInterval>,
    /// Exact number of underlying accepted mechanics substeps.
    pub accepted_mechanics_steps: usize,
    /// Fixed-grid branch-transition brackets at mechanics resolution.
    pub transitions: Vec<ProductionBranchTransition>,
    /// Bounded completion or first atomic refusal.
    pub termination: ProductionEventTrajectoryTermination,
}

impl ProductionControlTrajectory {
    /// Concatenate two exactly adjacent prefixes without replaying mechanics.
    pub fn concatenate(mut self, next: Self) -> Result<Self, ProductionCouplingError> {
        if self.last_accepted_checkpoint != next.start_checkpoint
            || !matches!(
                self.termination,
                ProductionEventTrajectoryTermination::StepLimitReached { .. }
            )
            || self
                .intervals
                .last()
                .zip(next.intervals.first())
                .is_some_and(|(left, right)| {
                    left.end_time_s.to_bits() != right.start_time_s.to_bits()
                        || left.state_after != right.state_before
                        || left.base_next_version != right.base_parent_version
                        || left.resolved_signed_gap_end_m.is_none_or(|gap| {
                            gap.to_bits() != right.resolved_signed_gap_start_m.to_bits()
                        })
                })
        {
            return Err(ProductionCouplingError::InvalidInput {
                field: "control trajectory concatenation",
            });
        }
        self.intervals.extend(next.intervals);
        self.accepted_mechanics_steps = self
            .accepted_mechanics_steps
            .checked_add(next.accepted_mechanics_steps)
            .ok_or(ProductionCouplingError::InvalidInput {
                field: "concatenated mechanics step count",
            })?;
        self.transitions.extend(next.transitions);
        self.last_accepted_checkpoint = next.last_accepted_checkpoint;
        self.termination = next.termination;
        Ok(self)
    }

    /// Coarsen already reduced controls while preserving exact impulse,
    /// endpoints, branch boundaries, and summed disc-work residual.
    pub fn coarsened(
        &self,
        maximum_source_intervals_per_output: usize,
    ) -> Result<Self, ProductionCouplingError> {
        if maximum_source_intervals_per_output == 0 {
            return Err(ProductionCouplingError::InvalidInput {
                field: "control coarsening factor",
            });
        }
        let mut output = Vec::new();
        output
            .try_reserve_exact(
                self.intervals
                    .len()
                    .div_ceil(maximum_source_intervals_per_output),
            )
            .map_err(|_| ProductionCouplingError::InvalidInput {
                field: "coarsened control capacity",
            })?;
        let mut active: Option<ProductionControlInterval> = None;
        let mut source_interval_count = 0_usize;
        for interval in &self.intervals {
            let must_flush = active.as_ref().is_some_and(|current| {
                current.branch != interval.branch
                    || source_interval_count == maximum_source_intervals_per_output
            });
            if must_flush {
                output.push(active.take().expect("active coarsening interval"));
                source_interval_count = 0;
            }
            if let Some(current) = &mut active {
                if current.end_time_s.to_bits() != interval.start_time_s.to_bits()
                    || current.state_after != interval.state_before
                    || current.base_next_version != interval.base_parent_version
                    || current.resolved_signed_gap_end_m.is_none_or(|gap| {
                        gap.to_bits() != interval.resolved_signed_gap_start_m.to_bits()
                    })
                {
                    return Err(ProductionCouplingError::InvalidInput {
                        field: "coarsened control lineage",
                    });
                }
                current.end_time_s = interval.end_time_s;
                current.mechanics_substeps = current
                    .mechanics_substeps
                    .checked_add(interval.mechanics_substeps)
                    .ok_or(ProductionCouplingError::InvalidInput {
                        field: "coarsened mechanics substeps",
                    })?;
                current.state_after = interval.state_after;
                current.normal_impulse_n_s += interval.normal_impulse_n_s;
                current.base_displacement_end_m = interval.base_displacement_end_m;
                current.base_velocity_end_m_per_s = interval.base_velocity_end_m_per_s;
                current.base_next_version = interval.base_next_version;
                current.disc_work_residual_j += interval.disc_work_residual_j;
                current.mechanical_energy_end_j = interval.mechanical_energy_end_j;
                current.resolved_signed_gap_end_m = interval.resolved_signed_gap_end_m;
                let duration_s = current.end_time_s - current.start_time_s;
                current.mean_normal_force_n = current.normal_impulse_n_s / duration_s;
            } else {
                active = Some(interval.clone());
            }
            source_interval_count += 1;
        }
        if let Some(active) = active {
            output.push(active);
        }
        Ok(Self {
            start_checkpoint: self.start_checkpoint.clone(),
            last_accepted_checkpoint: self.last_accepted_checkpoint.clone(),
            intervals: output,
            accepted_mechanics_steps: self.accepted_mechanics_steps,
            transitions: self.transitions.clone(),
            termination: self.termination.clone(),
        })
    }
}

/// Typed refusal; no caller checkpoint is changed on every error path.
#[derive(Debug, Clone, PartialEq)]
pub enum ProductionCouplingError {
    /// Model/checkpoint identity or version mismatch.
    CheckpointMismatch,
    /// Caller prepared cards from a different accepted outer checkpoint.
    CheckpointVersionMismatch { expected: u64, observed: u64 },
    /// Public checkpoint fields were altered without rebuilding its state binding.
    CheckpointIntegrityMismatch,
    /// A caller-owned adapter card is not bound to this case or world frame.
    InputIdentityMismatch { field: &'static str },
    /// Invalid outer scalar or identity.
    InvalidInput { field: &'static str },
    /// Profile support/curvature or profile-derived mass refused publication.
    ProfileContact(ContactDynamicsError),
    /// Profile-derived mass/inertia disagrees with the immutable mechanics model.
    ProfileModelMassMismatch,
    /// Cached resolved mass properties disagree with a fresh evaluation of the
    /// retained chart at the retained density.
    ResolvedProfileMassMismatch,
    /// A supposedly resolved specimen no longer matches its retained chart identity.
    ResolvedProfileIdentityMismatch {
        /// Identity retained by the resolved specimen.
        profile_identity: AxisymmetricIdentity,
        /// Identity retained by the specimen's actual chart.
        chart_identity: AxisymmetricIdentity,
    },
    /// The profile/moving-base bridge could not form one patch.
    Patch(PatchKinematicsError),
    /// Nonsmooth/unavailable curvature cannot be silently approximated.
    CurvatureUnavailable,
    /// Separation, unknown, or impact requires a separate event owner.
    UnsupportedMechanism { status: PatchContactStatus },
    /// Point-resultant mechanics cannot mix line-normal units into fs-mbd.
    UnsupportedLineNormalContact,
    /// Generic normal admission refused the sample.
    Normal(NormalContactError),
    /// Finite-patch topography filtering or its declared linearization refused.
    SurfaceExcitation(ProductionSurfaceExcitationError),
    /// Tangential generic admission/refusal.
    Tangential(TangentialContactError),
    /// Rolling generic admission/refusal, including work ownership overlap.
    Rolling(RollingContactError),
    /// Exterior-air admission/refusal.
    ExteriorAir(ExternalAirError),
    /// Thin-gap gas-film admission, topology, or exact-once-work refusal.
    AirFilm(AirFilmError),
    /// The input gas mechanism does not match the accepted gas checkpoint mode.
    GasChannelMismatch,
    /// The named exterior candidate was absent; alternatives are never auto-selected.
    ExteriorCandidateUnavailable,
    /// Reduced base explicitly refuses resolved/as-built/unsupported requests.
    Base(BaseResponseError),
    /// Resolved moving-contact modal base refused projection or time advance.
    ModalBase(RectangularModalBaseError),
    /// Checkpoint/proposal backend differs from the immutable selected base port.
    BaseBackendMismatch,
    /// The selected support backend cannot initialize a static contact preload.
    BaseStaticPreloadUnsupported,
    /// Caller-supplied initial pose and static support load do not satisfy the selected normal law.
    StaticPreloadMismatch {
        /// Requested equilibrium normal load [N].
        target_force_n: f64,
        /// Normal-law load at the supplied pose and preloaded support [N].
        observed_force_n: f64,
        /// Caller-declared relative acceptance tolerance.
        maximum_relative_mismatch: f64,
    },
    /// fs-mbd refused the actual summed wrench.
    Dynamics(fs_mbd::DynamicsError),
}

impl fmt::Display for ProductionCouplingError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for ProductionCouplingError {}

impl ProductionBasePort {
    /// Stable physical-model and selected-configuration identities.
    #[must_use]
    pub fn identity_parts(&self) -> (&str, &str) {
        match self {
            Self::ReducedOneMode(port) => {
                let identity = port.identity();
                (&identity.model_id, &identity.configuration_id)
            }
            Self::RectangularModal(port) => {
                let identity = port.identity();
                (&identity.model_id, &identity.configuration_id)
            }
        }
    }

    fn initial_checkpoint(&self) -> ProductionBaseCheckpoint {
        match self {
            Self::ReducedOneMode(port) => {
                ProductionBaseCheckpoint::ReducedOneMode(port.initial_checkpoint())
            }
            Self::RectangularModal(port) => {
                ProductionBaseCheckpoint::RectangularModal(port.initial_checkpoint())
            }
        }
    }

    fn initial_static_contact_checkpoint(
        &self,
        contact_point_world_m: Vec3,
        compressive_normal_force_on_base_n: f64,
    ) -> Result<ProductionBaseCheckpoint, ProductionCouplingError> {
        match self {
            Self::RectangularModal(port) => port
                .initial_static_contact_checkpoint(
                    [contact_point_world_m.x, contact_point_world_m.y, 0.0],
                    compressive_normal_force_on_base_n,
                )
                .map(ProductionBaseCheckpoint::RectangularModal)
                .map_err(ProductionCouplingError::ModalBase),
            Self::ReducedOneMode(_) => Err(ProductionCouplingError::BaseStaticPreloadUnsupported),
        }
    }

    fn surface_state(
        &self,
        checkpoint: &ProductionBaseCheckpoint,
        contact_point_world_m: Vec3,
    ) -> Result<(f64, f64), ProductionCouplingError> {
        match (self, checkpoint) {
            (Self::ReducedOneMode(_), ProductionBaseCheckpoint::ReducedOneMode(state)) => {
                Ok((state.modal_displacement_m(), state.modal_velocity_m_per_s()))
            }
            (Self::RectangularModal(port), ProductionBaseCheckpoint::RectangularModal(state)) => {
                let surface = port
                    .surface_state(
                        state,
                        [contact_point_world_m.x, contact_point_world_m.y, 0.0],
                    )
                    .map_err(ProductionCouplingError::ModalBase)?;
                Ok((surface.displacement_m, surface.velocity_m_per_s))
            }
            _ => Err(ProductionCouplingError::BaseBackendMismatch),
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn propose(
        &self,
        checkpoint: &ProductionBaseCheckpoint,
        step_id: String,
        duration_s: f64,
        compressive_normal_force_on_base_n: f64,
        contact_point_start_world_m: Vec3,
        contact_point_force_world_m: Vec3,
        contact_point_end_world_m: Vec3,
        legacy_load_progress_start: f64,
        legacy_load_progress_end: f64,
    ) -> Result<ProductionBaseStepProposal, ProductionCouplingError> {
        match (self, checkpoint) {
            (Self::ReducedOneMode(port), ProductionBaseCheckpoint::ReducedOneMode(state)) => {
                let proposal = port
                    .propose(
                        state,
                        &ReducedBaseStepInput {
                            step_id,
                            expected_version: state.accepted_version(),
                            duration_s,
                            compressive_normal_force_on_base_n,
                            load_progress_start: legacy_load_progress_start,
                            load_progress_end: legacy_load_progress_end,
                        },
                    )
                    .map_err(ProductionCouplingError::Base)?;
                let base = proposal.receipt();
                let receipt = ProductionBaseStepReceipt {
                    parent_version: base.parent_version,
                    next_version: base.next_version,
                    timestep_s: base.timestep_s,
                    compressive_normal_force_on_base_n: base.compressive_normal_force_on_base_n,
                    normal_reaction_on_disc_world_n: base.normal_reaction_on_disc_world_n,
                    modal_displacement_start_m: base.modal_displacement_start_m,
                    modal_displacement_end_m: base.modal_displacement_end_m,
                    modal_velocity_start_m_per_s: base.modal_velocity_start_m_per_s,
                    modal_velocity_end_m_per_s: base.modal_velocity_end_m_per_s,
                    stored_energy_change_j: base.stored_energy_change_j,
                    damping_work_j: base.damping_work_j,
                    external_contact_work_j: base.external_contact_work_j,
                    energy_closure_residual_j: base.energy_closure_residual_j,
                    end_support_reaction_norm_n: base.end_support_reaction_norm_n,
                };
                Ok(ProductionBaseStepProposal {
                    backend: ProductionBaseProposalBackend::ReducedOneMode(proposal),
                    receipt,
                })
            }
            (Self::RectangularModal(port), ProductionBaseCheckpoint::RectangularModal(state)) => {
                let point = |world: Vec3| [world.x, world.y, 0.0];
                let proposal = port
                    .propose(
                        state,
                        &RectangularModalBaseStepInput {
                            step_id,
                            expected_version: state.accepted_version(),
                            duration_s,
                            contact_point_start_base_m: point(contact_point_start_world_m),
                            contact_point_force_base_m: point(contact_point_force_world_m),
                            contact_point_end_base_m: point(contact_point_end_world_m),
                            compressive_normal_force_on_base_n,
                        },
                    )
                    .map_err(ProductionCouplingError::ModalBase)?;
                let base = proposal.receipt();
                let receipt = ProductionBaseStepReceipt {
                    parent_version: base.parent_version,
                    next_version: base.next_version,
                    timestep_s: base.duration_s,
                    compressive_normal_force_on_base_n: base.compressive_normal_force_on_base_n,
                    normal_reaction_on_disc_world_n: base.normal_reaction_on_disc_base_n,
                    modal_displacement_start_m: base.surface_start.displacement_m,
                    modal_displacement_end_m: base.surface_end.displacement_m,
                    modal_velocity_start_m_per_s: base.surface_start.velocity_m_per_s,
                    modal_velocity_end_m_per_s: base.surface_end.velocity_m_per_s,
                    stored_energy_change_j: base.stored_energy_change_j,
                    damping_work_j: base.viscous_dissipation_j,
                    external_contact_work_j: base.external_work_j,
                    energy_closure_residual_j: base.energy_closure_residual_j,
                    end_support_reaction_norm_n: compressive_normal_force_on_base_n,
                };
                Ok(ProductionBaseStepProposal {
                    backend: ProductionBaseProposalBackend::RectangularModal(proposal),
                    receipt,
                })
            }
            _ => Err(ProductionCouplingError::BaseBackendMismatch),
        }
    }

    fn accept(
        &self,
        checkpoint: &ProductionBaseCheckpoint,
        proposal: ProductionBaseStepProposal,
    ) -> Result<ProductionBaseCheckpoint, ProductionCouplingError> {
        match (self, checkpoint, proposal.backend) {
            (
                Self::ReducedOneMode(port),
                ProductionBaseCheckpoint::ReducedOneMode(state),
                ProductionBaseProposalBackend::ReducedOneMode(proposal),
            ) => port
                .accept(state, proposal)
                .map(ProductionBaseCheckpoint::ReducedOneMode)
                .map_err(ProductionCouplingError::Base),
            (
                Self::RectangularModal(port),
                ProductionBaseCheckpoint::RectangularModal(state),
                ProductionBaseProposalBackend::RectangularModal(proposal),
            ) => port
                .accept(state, proposal)
                .map(ProductionBaseCheckpoint::RectangularModal)
                .map_err(ProductionCouplingError::ModalBase),
            _ => Err(ProductionCouplingError::BaseBackendMismatch),
        }
    }
}

impl<'run> AdmittedAxisymmetricProfile<'run> {
    fn resolve_contact(
        &self,
        disc_state: RigidBodyState,
        cx: &Cx<'_>,
    ) -> Result<ProfileContactPatchGeometry, ProductionCouplingError> {
        profile_contact_patch_geometry_from_query_session(
            &self.query_session,
            self.profile.identity,
            self.mass_properties,
            disc_state.pose(),
            cx,
        )
        .map_err(ProductionCouplingError::ProfileContact)
    }

    fn bind_contact_at_states(
        &self,
        input: &mut ProductionCouplingStepInput,
        disc_state: RigidBodyState,
        base_state: &ProductionBaseCheckpoint,
        cx: &Cx<'_>,
    ) -> Result<ProfileContactPatchGeometry, ProductionCouplingError> {
        let resolved = self.resolve_contact(disc_state, cx)?;
        let (base_displacement_m, base_velocity_m_per_s) = self
            .model
            .base_port
            .surface_state(base_state, resolved.support.contact.point_world_m)?;
        bind_horizontal_plane_axisymmetric_profile_contact_input(
            input,
            resolved,
            self.model.disc_mass_properties,
            disc_state,
            base_displacement_m,
            base_velocity_m_per_s,
        )
    }

    pub(crate) fn bind_contact(
        &self,
        input: &mut ProductionCouplingStepInput,
        checkpoint: &ProductionCouplingCheckpoint,
        cx: &Cx<'_>,
    ) -> Result<ProfileContactPatchGeometry, ProductionCouplingError> {
        self.model.validate_checkpoint(checkpoint)?;
        self.bind_contact_after_checkpoint_validation(input, checkpoint, cx)
    }

    pub(crate) fn bind_contact_after_checkpoint_validation(
        &self,
        input: &mut ProductionCouplingStepInput,
        checkpoint: &ProductionCouplingCheckpoint,
        cx: &Cx<'_>,
    ) -> Result<ProfileContactPatchGeometry, ProductionCouplingError> {
        let resolved =
            self.bind_contact_at_states(input, checkpoint.disc_state, &checkpoint.base_state, cx)?;
        input.expected_checkpoint_version = checkpoint.committed_version;
        Ok(resolved)
    }

    pub(crate) fn contact_geometry(
        &self,
        pose: Pose,
        cx: &Cx<'_>,
    ) -> Result<ProfileContactGeometry, ContactDynamicsError> {
        profile_contact_geometry_from_query_session(
            &self.query_session,
            self.mass_properties,
            pose,
            cx,
        )
    }

    pub(crate) fn run_eventful_profile_midpoint_control_trajectory_observed(
        &self,
        start_checkpoint: ProductionCouplingCheckpoint,
        maximum_accepted_steps: usize,
        mechanics_steps_per_control_interval: usize,
        surface_traces: &AdmittedSurfaceTracePair,
        cx: &Cx<'_>,
        rebind_input_for_checkpoint: impl FnMut(
            &ProductionCouplingCheckpoint,
            &mut Option<ProductionCouplingStepInput>,
        ) -> Result<(), ProductionCouplingError>,
        mut refresh_midpoint_input: impl FnMut(
            &mut ProductionCouplingStepInput,
            RigidBodyState,
        ) -> Result<(), ProductionCouplingError>,
        observe_modal_audio_step: impl for<'step> FnMut(ProductionModalAudioStep<'step>),
    ) -> Result<ProductionControlTrajectory, ProductionCouplingError> {
        self.model.run_eventful_control_trajectory_observed_with(
            start_checkpoint,
            maximum_accepted_steps,
            mechanics_steps_per_control_interval,
            CheckpointValidation::TrustedInternalSuccessor,
            SurfaceTraceEvaluation::Admitted(surface_traces),
            rebind_input_for_checkpoint,
            |checkpoint, input| {
                self.model.step_eventful_profile_midpoint_with(
                    checkpoint,
                    ProductionMidpointInput::Reusable(input),
                    self.profile,
                    CheckpointValidation::TrustedInternalSuccessor,
                    Some(self),
                    SurfaceTraceEvaluation::Admitted(surface_traces),
                    cx,
                    |input, state| refresh_midpoint_input(input, state),
                )
            },
            observe_modal_audio_step,
        )
    }
}

impl ProductionCouplingModel {
    /// Revalidate immutable public profile state once for a private run.
    pub(crate) fn admit_axisymmetric_profile<'run>(
        &'run self,
        profile: &'run ResolvedDiscProfile,
        cx: &Cx<'_>,
    ) -> Result<AdmittedAxisymmetricProfile<'run>, ProductionCouplingError> {
        validate_identity(&self.identity)?;
        let chart_identity = profile.chart.construction_certificate().identity;
        if profile.identity != chart_identity {
            return Err(ProductionCouplingError::ResolvedProfileIdentityMismatch {
                profile_identity: profile.identity,
                chart_identity,
            });
        }
        let mass_properties = profile
            .chart
            .mass_properties(profile.density_kg_per_m3, cx)
            .map_err(|detail| {
                ProductionCouplingError::ProfileContact(ContactDynamicsError::ProfileMassRefusal {
                    detail,
                })
            })?;
        if profile.mass_properties != mass_properties {
            return Err(ProductionCouplingError::ResolvedProfileMassMismatch);
        }
        let disc_mass_properties = profile_mass_to_mbd(mass_properties)
            .map_err(ProductionCouplingError::ProfileContact)?;
        if disc_mass_properties != self.disc_mass_properties {
            return Err(ProductionCouplingError::ProfileModelMassMismatch);
        }
        let query_session: AxisymmetricQuerySession<'run> =
            profile.chart.admit_query_session(cx).map_err(|detail| {
                let detail = match detail {
                    AxisymmetricError::Cancelled => AxisymmetricSupportError::Cancelled,
                    detail => AxisymmetricSupportError::InvalidChart(detail),
                };
                ProductionCouplingError::ProfileContact(
                    ContactDynamicsError::ProfileSupportRefusal { detail },
                )
            })?;
        Ok(AdmittedAxisymmetricProfile {
            model: self,
            profile,
            mass_properties,
            query_session,
        })
    }

    /// Verify that a checkpoint belongs to this exact model lineage and that
    /// none of its public or private accepted state has been altered.
    pub fn validate_checkpoint(
        &self,
        checkpoint: &ProductionCouplingCheckpoint,
    ) -> Result<(), ProductionCouplingError> {
        validate_identity(&self.identity)?;
        if checkpoint.identity != self.identity {
            return Err(ProductionCouplingError::CheckpointMismatch);
        }
        if checkpoint.checkpoint_fingerprint
            != production_checkpoint_fingerprint(
                &checkpoint.identity,
                checkpoint.committed_version,
                checkpoint.disc_state,
                &checkpoint.normal_state,
                &checkpoint.tangential_state,
                &checkpoint.rolling_state,
                &checkpoint.gas_channel_state,
                &checkpoint.base_state,
            )
        {
            return Err(ProductionCouplingError::CheckpointIntegrityMismatch);
        }
        Ok(())
    }

    /// Resolves the initial smooth profile patch against this model's zero-version base state.
    ///
    /// The immutable mechanics mass must exactly match the mass derived from
    /// `profile`; this prevents a real profile patch from
    /// being advanced by a stale synthetic rigid-body model. This initializer
    /// exists because constructing the first tangential checkpoint requires the
    /// same profile-native normal patch that the first coupled step will use.
    pub fn bind_initial_horizontal_plane_axisymmetric_profile_contact(
        &self,
        input: &mut ProductionCouplingStepInput,
        profile: &ResolvedDiscProfile,
        disc_state: RigidBodyState,
        cx: &Cx<'_>,
    ) -> Result<ProfileContactPatchGeometry, ProductionCouplingError> {
        validate_identity(&self.identity)?;
        let base_state = self.base_port.initial_checkpoint();
        self.bind_initial_horizontal_plane_axisymmetric_profile_contact_with_base_state(
            input,
            profile,
            disc_state,
            &base_state,
            cx,
        )
    }

    fn bind_initial_horizontal_plane_axisymmetric_profile_contact_with_base_state(
        &self,
        input: &mut ProductionCouplingStepInput,
        profile: &ResolvedDiscProfile,
        disc_state: RigidBodyState,
        base_state: &ProductionBaseCheckpoint,
        cx: &Cx<'_>,
    ) -> Result<ProfileContactPatchGeometry, ProductionCouplingError> {
        validate_identity(&self.identity)?;
        let (resolved, disc_mass_properties) =
            self.resolve_axisymmetric_profile_contact(profile, disc_state, cx)?;
        let (base_displacement_m, base_velocity_m_per_s) = self
            .base_port
            .surface_state(base_state, resolved.support.contact.point_world_m)?;
        let resolved = bind_horizontal_plane_axisymmetric_profile_contact_input(
            input,
            resolved,
            disc_mass_properties,
            disc_state,
            base_displacement_m,
            base_velocity_m_per_s,
        )?;
        input.expected_checkpoint_version = 0;
        Ok(resolved)
    }

    /// Build the first complete coupled checkpoint from one resolved profile.
    ///
    /// This is the profile-native admission boundary used by product drivers.
    /// It binds the actual support point and principal curvatures, evaluates
    /// the selected normal law at the supplied state, derives the finite patch
    /// seen by the selected tangential law, starts rolling history at zero,
    /// and seals all channel states into one checkpoint. No material name,
    /// contact radius, normal force, or coefficient is inferred here: every
    /// constitutive choice remains in `input`, while `gas_channel_state`
    /// remains an explicitly selected mutually exclusive gas mechanism.
    ///
    /// The method mutates only the caller-owned template fields whose values
    /// are derived from the admitted profile and initial checkpoint. A refusal
    /// returns no partially accepted coupled state.
    pub fn initialize_horizontal_plane_axisymmetric_profile_trajectory(
        &self,
        input: &mut ProductionCouplingStepInput,
        profile: &ResolvedDiscProfile,
        disc_state: RigidBodyState,
        normal_state: NormalPatchEmbedState,
        gas_channel_state: GasChannelState,
        maximum_tangential_work_keys: usize,
        cx: &Cx<'_>,
    ) -> Result<ProductionCouplingCheckpoint, ProductionCouplingError> {
        self.bind_initial_horizontal_plane_axisymmetric_profile_contact(
            input, profile, disc_state, cx,
        )?;
        let mut normal_input = input.normal.clone();
        normal_input.state = normal_state.clone();
        normal_input.time_s = 0.0;
        normal_input.step_s = input.duration_s;
        let (patch_kinematics, normal) = if let Some(surface) = &input.surface_geometry {
            let (patch, normal, _) = resolve_surface_geometry_contact(
                input.patch.clone(),
                &normal_input,
                surface,
                SurfaceTraceEvaluation::Checked,
            )?;
            (patch, normal)
        } else {
            let patch = compute_moving_one_mode_patch_kinematics(input.patch.clone())
                .map_err(ProductionCouplingError::Patch)?;
            normal_input.kinematics = patch.clone();
            let normal = match evaluate_normal_contact(&normal_input)
                .map_err(ProductionCouplingError::Normal)?
            {
                EulerNormalContactOutcome::Active(active) => active,
                EulerNormalContactOutcome::InactiveSeparated { .. } => {
                    return Err(ProductionCouplingError::UnsupportedMechanism {
                        status: PatchContactStatus::Separated,
                    });
                }
            };
            (patch, normal)
        };
        let normal_patch = normal_patch_view(&normal, &patch_kinematics)?;
        let tangential_state = self
            .tangential_adapter
            .initial_state(
                &normal_patch,
                &input.tangential.interface,
                maximum_tangential_work_keys,
            )
            .map_err(ProductionCouplingError::Tangential)?;
        self.initial_checkpoint(
            disc_state,
            normal_state,
            tangential_state,
            RollingContactState::zero(),
            gas_channel_state,
        )
    }

    fn resolve_horizontal_plane_axisymmetric_static_contact_at_state(
        &self,
        input: &mut ProductionCouplingStepInput,
        profile: &ResolvedDiscProfile,
        disc_state: RigidBodyState,
        normal_state: &NormalPatchEmbedState,
        static_contact_force_n: f64,
        cx: &Cx<'_>,
    ) -> Result<
        (
            ProductionBaseCheckpoint,
            Option<(PatchKinematics, ActiveNormalContact)>,
        ),
        ProductionCouplingError,
    > {
        let (resolved, _) = self.resolve_axisymmetric_profile_contact(profile, disc_state, cx)?;
        let base_state = self.base_port.initial_static_contact_checkpoint(
            resolved.support.contact.point_world_m,
            static_contact_force_n,
        )?;
        self.bind_initial_horizontal_plane_axisymmetric_profile_contact_with_base_state(
            input,
            profile,
            disc_state,
            &base_state,
            cx,
        )?;
        let mut normal_input = input.normal.clone();
        normal_input.state = normal_state.clone();
        normal_input.time_s = 0.0;
        normal_input.step_s = input.duration_s;
        let resolved_contact = if let Some(surface) = &input.surface_geometry {
            match resolve_surface_geometry_contact(
                input.patch.clone(),
                &normal_input,
                surface,
                SurfaceTraceEvaluation::Checked,
            ) {
                Ok((patch, normal, _)) => Some((patch, normal)),
                Err(ProductionCouplingError::UnsupportedMechanism {
                    status: PatchContactStatus::Separated,
                }) => None,
                Err(error) => return Err(error),
            }
        } else {
            let patch = compute_moving_one_mode_patch_kinematics(input.patch.clone())
                .map_err(ProductionCouplingError::Patch)?;
            normal_input.kinematics = patch.clone();
            match evaluate_normal_contact(&normal_input).map_err(ProductionCouplingError::Normal)? {
                EulerNormalContactOutcome::Active(active) => Some((patch, active)),
                EulerNormalContactOutcome::InactiveSeparated { .. } => None,
            }
        };
        Ok((base_state, resolved_contact))
    }

    /// Build a coupled checkpoint with the resolved modal support already in static equilibrium.
    ///
    /// `disc_state` must already encode the contact approach consistent with
    /// `static_contact_force_n`. The method independently initializes the
    /// selected support backend under that force, rebinds the true profile
    /// against its resulting local displacement, evaluates the selected normal
    /// law, and refuses unless the two forces agree within the declared relative
    /// tolerance. This prevents a product driver from injecting an artificial
    /// start-up ring by pairing a gravity-loaded disc with an unloaded plate.
    #[allow(clippy::too_many_arguments)]
    pub fn initialize_horizontal_plane_axisymmetric_profile_static_contact(
        &self,
        input: &mut ProductionCouplingStepInput,
        profile: &ResolvedDiscProfile,
        disc_state: RigidBodyState,
        normal_state: NormalPatchEmbedState,
        gas_channel_state: GasChannelState,
        maximum_tangential_work_keys: usize,
        strict_sequential_tangential_work: bool,
        static_contact_force_n: f64,
        maximum_relative_force_mismatch: f64,
        cx: &Cx<'_>,
    ) -> Result<ProductionCouplingCheckpoint, ProductionCouplingError> {
        if !(static_contact_force_n.is_finite()
            && static_contact_force_n > 0.0
            && maximum_relative_force_mismatch.is_finite()
            && maximum_relative_force_mismatch >= 0.0)
        {
            return Err(ProductionCouplingError::InvalidInput {
                field: "static contact initialization",
            });
        }
        let (base_state, resolved_contact) = self
            .resolve_horizontal_plane_axisymmetric_static_contact_at_state(
                input,
                profile,
                disc_state,
                &normal_state,
                static_contact_force_n,
                cx,
            )?;
        let Some((patch_kinematics, normal)) = resolved_contact else {
            return Err(ProductionCouplingError::StaticPreloadMismatch {
                target_force_n: static_contact_force_n,
                observed_force_n: 0.0,
                maximum_relative_mismatch: maximum_relative_force_mismatch,
            });
        };
        let observed_force_n = observed_static_normal_force(input, &normal)?;
        let relative_mismatch = (observed_force_n - static_contact_force_n).abs()
            / static_contact_force_n.max(f64::MIN_POSITIVE);
        if relative_mismatch > maximum_relative_force_mismatch {
            return Err(ProductionCouplingError::StaticPreloadMismatch {
                target_force_n: static_contact_force_n,
                observed_force_n,
                maximum_relative_mismatch: maximum_relative_force_mismatch,
            });
        }
        let normal_patch = normal_patch_view(&normal, &patch_kinematics)?;
        let tangential_state = if strict_sequential_tangential_work {
            self.tangential_adapter.initial_state_strict_sequence(
                &normal_patch,
                &input.tangential.interface,
                maximum_tangential_work_keys,
                input.tangential.work_ownership.clone(),
            )
        } else {
            self.tangential_adapter.initial_state(
                &normal_patch,
                &input.tangential.interface,
                maximum_tangential_work_keys,
            )
        }
        .map_err(ProductionCouplingError::Tangential)?;
        self.initial_checkpoint_with_base_state(
            disc_state,
            normal_state,
            tangential_state,
            RollingContactState::zero(),
            gas_channel_state,
            base_state,
        )
    }

    /// Solve the vertical pose coordinate for a declared static normal load.
    ///
    /// The orientation and both momenta from `seed_state` are preserved. The
    /// profile is first placed at geometric contact with the undeformed plane;
    /// the selected modal base is then initialized under `target_force_n`, and
    /// a deterministic signed bisection finds the vertical contact offset whose
    /// selected finite-patch law returns that same force. Positive offset moves
    /// inward; negative offset permits a realized surface summit to touch before
    /// the smooth profile does. Every trial goes through the ordinary
    /// profile/patch/normal adapters and therefore retains their applicability
    /// refusals rather than extrapolating a Hertz law.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_horizontal_plane_axisymmetric_static_contact_state(
        &self,
        input: &mut ProductionCouplingStepInput,
        profile: &ResolvedDiscProfile,
        seed_state: RigidBodyState,
        normal_state: NormalPatchEmbedState,
        target_force_n: f64,
        maximum_relative_force_mismatch: f64,
        maximum_iterations: usize,
        cx: &Cx<'_>,
    ) -> Result<RigidBodyState, ProductionCouplingError> {
        if !(target_force_n.is_finite()
            && target_force_n > 0.0
            && maximum_relative_force_mismatch.is_finite()
            && maximum_relative_force_mismatch > 0.0
            && maximum_iterations > 0)
        {
            return Err(ProductionCouplingError::InvalidInput {
                field: "static contact pose solve",
            });
        }
        let (seed_contact, _) =
            self.resolve_axisymmetric_profile_contact(profile, seed_state, cx)?;
        let seed_position = seed_state.pose().position_world();
        let grounded_position = Vec3::new(
            seed_position.x,
            seed_position.y,
            seed_position.z - seed_contact.support.contact.gap_m,
        );
        let grounded_pose = Pose::new(grounded_position, seed_state.pose().orientation())
            .map_err(ProductionCouplingError::Dynamics)?;
        let grounded_state = RigidBodyState::new(
            grounded_pose,
            seed_state.linear_momentum_world(),
            seed_state.angular_momentum_body(),
        )
        .map_err(ProductionCouplingError::Dynamics)?;
        let (grounded_contact, _) =
            self.resolve_axisymmetric_profile_contact(profile, grounded_state, cx)?;
        let base_state = self.base_port.initial_static_contact_checkpoint(
            grounded_contact.support.contact.point_world_m,
            target_force_n,
        )?;
        let (base_displacement_m, _) = self
            .base_port
            .surface_state(&base_state, grounded_contact.support.contact.point_world_m)?;
        let state_at_approach =
            |approach_m: f64| -> Result<RigidBodyState, ProductionCouplingError> {
                let pose = Pose::new(
                    Vec3::new(
                        grounded_position.x,
                        grounded_position.y,
                        grounded_position.z + base_displacement_m - approach_m,
                    ),
                    grounded_state.pose().orientation(),
                )
                .map_err(ProductionCouplingError::Dynamics)?;
                RigidBodyState::new(
                    pose,
                    grounded_state.linear_momentum_world(),
                    grounded_state.angular_momentum_body(),
                )
                .map_err(ProductionCouplingError::Dynamics)
            };
        let canonical_template = input.clone();
        let force_at_approach = |approach_m: f64| -> Result<
            (RigidBodyState, f64, ProductionCouplingStepInput),
            ProductionCouplingError,
        > {
            let state = state_at_approach(approach_m)?;
            // Root evaluations must not inherit derived patch fields from a
            // preceding bracket candidate.  Retain the fully rebound template
            // alongside its state so the selected root has one coherent pair.
            let mut template = canonical_template.clone();
            let (_, resolved_contact) = self
                .resolve_horizontal_plane_axisymmetric_static_contact_at_state(
                    &mut template,
                    profile,
                    state,
                    &normal_state,
                    target_force_n,
                    cx,
                )?;
            let force_n = resolved_contact
                .as_ref()
                .map(|(_, normal)| observed_static_normal_force(&template, normal))
                .transpose()?
                .unwrap_or(0.0);
            Ok((state, force_n, template))
        };

        let curvature_scale_m = [
            grounded_contact.curvature.meridional_m_inverse,
            grounded_contact.curvature.azimuthal_m_inverse,
        ]
        .into_iter()
        .filter(|value| value.is_finite() && *value > 0.0)
        .map(|value| value.recip())
        .fold(f64::INFINITY, f64::min);
        if !curvature_scale_m.is_finite() {
            return Err(ProductionCouplingError::CurvatureUnavailable);
        }
        let bracket_step_m = (curvature_scale_m * 1.0e-9).max(1.0e-15);
        let (zero_state, zero_force_n, zero_template) = force_at_approach(0.0)?;
        let (
            mut lower_m,
            lower_state,
            lower_force_n,
            lower_template,
            mut upper_m,
            upper_state,
            upper_force_n,
            upper_template,
        ) = if zero_force_n <= target_force_n {
            let mut trial_m = bracket_step_m;
            let (upper_state, upper_force_n, upper_template) = loop {
                let trial = force_at_approach(trial_m)?;
                if trial.1 >= target_force_n {
                    break trial;
                }
                trial_m *= 2.0;
                if !trial_m.is_finite() || trial_m > curvature_scale_m {
                    return Err(ProductionCouplingError::StaticPreloadMismatch {
                        target_force_n,
                        observed_force_n: trial.1,
                        maximum_relative_mismatch: maximum_relative_force_mismatch,
                    });
                }
            };
            (
                0.0,
                zero_state,
                zero_force_n,
                zero_template,
                trial_m,
                upper_state,
                upper_force_n,
                upper_template,
            )
        } else {
            // Positive topography can already compress the interface at the
            // smooth-profile touching pose. Search outward until the selected
            // surface/contact model falls below the declared preload; zero is
            // not a valid lower bracket in that case.
            let mut trial_m = -bracket_step_m;
            let (lower_state, lower_force_n, lower_template) = loop {
                let trial = force_at_approach(trial_m)?;
                if trial.1 <= target_force_n {
                    break trial;
                }
                trial_m *= 2.0;
                if !trial_m.is_finite() || trial_m.abs() > curvature_scale_m {
                    return Err(ProductionCouplingError::StaticPreloadMismatch {
                        target_force_n,
                        observed_force_n: trial.1,
                        maximum_relative_mismatch: maximum_relative_force_mismatch,
                    });
                }
            };
            (
                trial_m,
                lower_state,
                lower_force_n,
                lower_template,
                0.0,
                zero_state,
                zero_force_n,
                zero_template,
            )
        };
        let (mut best_state, mut best_force_n, mut best_template) =
            if (lower_force_n - target_force_n).abs() < (upper_force_n - target_force_n).abs() {
                (lower_state, lower_force_n, lower_template)
            } else {
                (upper_state, upper_force_n, upper_template)
            };
        for _ in 0..maximum_iterations {
            let midpoint_m = 0.5 * (lower_m + upper_m);
            let (candidate, force_n, candidate_template) = force_at_approach(midpoint_m)?;
            if (force_n - target_force_n).abs() / target_force_n.max(f64::MIN_POSITIVE)
                <= maximum_relative_force_mismatch
            {
                *input = candidate_template;
                return Ok(candidate);
            }
            if (force_n - target_force_n).abs() < (best_force_n - target_force_n).abs() {
                best_state = candidate;
                best_force_n = force_n;
                best_template = candidate_template;
            }
            if force_n < target_force_n {
                lower_m = midpoint_m;
            } else {
                upper_m = midpoint_m;
            }
        }
        let relative_mismatch =
            (best_force_n - target_force_n).abs() / target_force_n.max(f64::MIN_POSITIVE);
        if relative_mismatch <= maximum_relative_force_mismatch {
            *input = best_template;
            Ok(best_state)
        } else {
            Err(ProductionCouplingError::StaticPreloadMismatch {
                target_force_n,
                observed_force_n: best_force_n,
                maximum_relative_mismatch: maximum_relative_force_mismatch,
            })
        }
    }

    /// Rebuilds one smooth profile patch directly from an accepted coupled checkpoint.
    ///
    /// Base displacement/velocity and rigid-body state come from the private,
    /// integrity-checked checkpoint, so callers cannot accidentally classify
    /// contact with stale base motion. The horizontal base contributes zero
    /// relative-gap curvature.
    pub fn bind_horizontal_plane_axisymmetric_profile_contact(
        &self,
        input: &mut ProductionCouplingStepInput,
        profile: &ResolvedDiscProfile,
        checkpoint: &ProductionCouplingCheckpoint,
        cx: &Cx<'_>,
    ) -> Result<ProfileContactPatchGeometry, ProductionCouplingError> {
        self.validate_checkpoint(checkpoint)?;
        let resolved = self.bind_horizontal_plane_axisymmetric_profile_contact_at_states(
            input,
            profile,
            checkpoint.disc_state,
            &checkpoint.base_state,
            cx,
        )?;
        input.expected_checkpoint_version = checkpoint.committed_version;
        Ok(resolved)
    }

    fn bind_horizontal_plane_axisymmetric_profile_contact_at_states(
        &self,
        input: &mut ProductionCouplingStepInput,
        profile: &ResolvedDiscProfile,
        disc_state: RigidBodyState,
        base_state: &ProductionBaseCheckpoint,
        cx: &Cx<'_>,
    ) -> Result<ProfileContactPatchGeometry, ProductionCouplingError> {
        let (resolved, disc_mass_properties) =
            self.resolve_axisymmetric_profile_contact(profile, disc_state, cx)?;
        let (base_displacement_m, base_velocity_m_per_s) = self
            .base_port
            .surface_state(base_state, resolved.support.contact.point_world_m)?;
        bind_horizontal_plane_axisymmetric_profile_contact_input(
            input,
            resolved,
            disc_mass_properties,
            disc_state,
            base_displacement_m,
            base_velocity_m_per_s,
        )
    }

    fn resolve_axisymmetric_profile_contact(
        &self,
        profile: &ResolvedDiscProfile,
        disc_state: RigidBodyState,
        cx: &Cx<'_>,
    ) -> Result<(ProfileContactPatchGeometry, MassProperties), ProductionCouplingError> {
        let chart_identity = profile.chart.construction_certificate().identity;
        if profile.identity != chart_identity {
            return Err(ProductionCouplingError::ResolvedProfileIdentityMismatch {
                profile_identity: profile.identity,
                chart_identity,
            });
        }
        // `ResolvedDiscProfile` is publicly constructible, so its cached mass
        // cannot be treated as an admission token. Re-evaluate the same chart
        // and density at every public binding boundary before using either the
        // cached properties or a mechanics model that may have copied them.
        // This is deterministic exact-formula work over the profile segments;
        // a future sealed admitted-profile token may cache it without weakening
        // this trust boundary.
        let recomputed_mass_properties = profile
            .chart
            .mass_properties(profile.density_kg_per_m3, cx)
            .map_err(|detail| {
                ProductionCouplingError::ProfileContact(ContactDynamicsError::ProfileMassRefusal {
                    detail,
                })
            })?;
        if profile.mass_properties != recomputed_mass_properties {
            return Err(ProductionCouplingError::ResolvedProfileMassMismatch);
        }
        let resolved = profile_contact_patch_geometry_from_mass(
            &profile.chart,
            recomputed_mass_properties,
            disc_state.pose(),
            cx,
        )
        .map_err(ProductionCouplingError::ProfileContact)?;
        let disc_mass_properties = profile_mass_to_mbd(resolved.support.mass_properties)
            .map_err(ProductionCouplingError::ProfileContact)?;
        if disc_mass_properties != self.disc_mass_properties {
            return Err(ProductionCouplingError::ProfileModelMassMismatch);
        }
        Ok((resolved, disc_mass_properties))
    }

    /// Validates immutable identities and creates one complete initial checkpoint.
    pub fn initial_checkpoint(
        &self,
        disc_state: RigidBodyState,
        normal_state: NormalPatchEmbedState,
        tangential_state: TangentialContactState,
        rolling_state: RollingContactState,
        gas_channel_state: GasChannelState,
    ) -> Result<ProductionCouplingCheckpoint, ProductionCouplingError> {
        let base_state = self.base_port.initial_checkpoint();
        self.initial_checkpoint_with_base_state(
            disc_state,
            normal_state,
            tangential_state,
            rolling_state,
            gas_channel_state,
            base_state,
        )
    }

    fn initial_checkpoint_with_base_state(
        &self,
        disc_state: RigidBodyState,
        normal_state: NormalPatchEmbedState,
        tangential_state: TangentialContactState,
        rolling_state: RollingContactState,
        gas_channel_state: GasChannelState,
        base_state: ProductionBaseCheckpoint,
    ) -> Result<ProductionCouplingCheckpoint, ProductionCouplingError> {
        validate_identity(&self.identity)?;
        if !disc_state.pose().position_world().is_finite() {
            return Err(ProductionCouplingError::InvalidInput {
                field: "disc_state",
            });
        }
        validate_gas_state_identity(&self.identity, &gas_channel_state)?;
        let checkpoint_fingerprint = production_checkpoint_fingerprint(
            &self.identity,
            0,
            disc_state,
            &normal_state,
            &tangential_state,
            &rolling_state,
            &gas_channel_state,
            &base_state,
        );
        Ok(ProductionCouplingCheckpoint {
            identity: self.identity.clone(),
            committed_version: 0,
            disc_state,
            checkpoint_fingerprint,
            normal_state,
            tangential_state,
            rolling_state,
            gas_channel_state,
            base_state,
        })
    }

    /// Attempts one homogeneous smooth substep and atomically advances every channel only after fs-mbd accepts.
    pub fn step(
        &self,
        checkpoint: &ProductionCouplingCheckpoint,
        input: &ProductionCouplingStepInput,
    ) -> Result<(ProductionCouplingCheckpoint, ProductionCouplingReceipt), ProductionCouplingError>
    {
        self.step_with_resolved_contact(checkpoint, input, None)
    }

    fn step_with_resolved_contact(
        &self,
        checkpoint: &ProductionCouplingCheckpoint,
        input: &ProductionCouplingStepInput,
        resolved_contact: Option<ResolvedProductionContact>,
    ) -> Result<(ProductionCouplingCheckpoint, ProductionCouplingReceipt), ProductionCouplingError>
    {
        self.validate_checkpoint(checkpoint)?;
        if input.expected_checkpoint_version != checkpoint.committed_version {
            return Err(ProductionCouplingError::CheckpointVersionMismatch {
                expected: input.expected_checkpoint_version,
                observed: checkpoint.committed_version,
            });
        }
        validate_step_identity(&self.identity, input)?;
        if !(input.duration_s.is_finite()
            && input.duration_s > 0.0
            && input.time_s.is_finite()
            && input.time_s >= 0.0)
        {
            return Err(ProductionCouplingError::InvalidInput {
                field: "duration_s or time_s",
            });
        }
        if input.surface_excitation.is_some() && input.surface_geometry.is_some() {
            return Err(ProductionCouplingError::InvalidInput {
                field: "mutually exclusive surface coupling modes",
            });
        }
        let resolved_contact = match resolved_contact {
            Some(contact) => contact,
            None => self.resolve_active_contact(checkpoint, input)?,
        };
        let prepared = self.prepare_contact_channels(
            checkpoint,
            input,
            checkpoint.disc_state,
            resolved_contact,
            false,
        )?;
        let contact_point = prepared.patch_kinematics.base_point.point_world;
        let base = self.base_port.propose(
            &checkpoint.base_state,
            input.base_step_id.clone(),
            input.duration_s,
            prepared.normal_force_n,
            contact_point,
            contact_point,
            contact_point,
            input.base_load_progress_start,
            input.base_load_progress_end,
        )?;
        let rigid_step = RigidBodyIntegrator::new(self.gravity)
            .step(
                checkpoint.disc_state,
                self.disc_mass_properties,
                Wrench {
                    force_world: prepared.total_force_world_n,
                    torque_body: checkpoint
                        .disc_state
                        .pose()
                        .orientation()
                        .rotate_world_to_body(prepared.total_moment_about_com_world_n_m),
                },
                input.duration_s,
            )
            .map_err(ProductionCouplingError::Dynamics)?;
        self.commit_contact_step(checkpoint, prepared, base, rigid_step)
    }

    fn prepare_contact_channels(
        &self,
        checkpoint: &ProductionCouplingCheckpoint,
        input: &ProductionCouplingStepInput,
        disc_state: RigidBodyState,
        resolved_contact: ResolvedProductionContact,
        apply_tangential_midpoint_probe: bool,
    ) -> Result<PreparedProductionContact, ProductionCouplingError> {
        let ResolvedProductionContact {
            patch_kinematics,
            normal,
            surface_geometry,
        } = resolved_contact;
        if matches!(
            patch_kinematics.patch.curvature,
            CurvatureMetadata::Unavailable { .. }
        ) {
            return Err(ProductionCouplingError::CurvatureUnavailable);
        }
        let admitted_status = match input.normal.integration_regime {
            NormalContactIntegrationRegime::SmoothQuasistatic => matches!(
                patch_kinematics.status,
                PatchContactStatus::Touching | PatchContactStatus::Grazing
            ),
            NormalContactIntegrationRegime::CompliantTransient => matches!(
                patch_kinematics.status,
                PatchContactStatus::Approaching
                    | PatchContactStatus::Receding
                    | PatchContactStatus::Touching
                    | PatchContactStatus::Grazing
                    | PatchContactStatus::ImpactCandidate
            ),
        };
        if !admitted_status {
            return Err(ProductionCouplingError::UnsupportedMechanism {
                status: patch_kinematics.status,
            });
        }

        let surface_excitation = input
            .surface_excitation
            .as_ref()
            .map(|surface| {
                evaluate_surface_excitation_for_normal(
                    &normal,
                    &surface.interface,
                    surface.surface_a.as_motion(),
                    surface.surface_b.as_motion(),
                    surface.travel_angle_from_patch_major_rad,
                    surface.maximum_linearized_height_fraction,
                )
            })
            .transpose()
            .map_err(ProductionCouplingError::SurfaceExcitation)?;
        let normal_patch = normal_patch_view(&normal, &patch_kinematics)?;
        let rolling_patch = rolling_patch_receipt(&normal, &patch_kinematics)?;

        let mut tangential_request = input.tangential.clone();
        tangential_request.expected_state_version = checkpoint.tangential_state.committed_version();
        tangential_request.patch_kinematics = patch_kinematics.clone();
        tangential_request.normal_patch = normal_patch;
        tangential_request.dt_s = input.duration_s;
        let tangential = self
            .tangential_adapter
            .prepare(&checkpoint.tangential_state, &tangential_request)
            .map_err(ProductionCouplingError::Tangential)?;
        let (applied_tangential_force_world_n, applied_tangential_free_torsional_torque_world_nm) =
            if apply_tangential_midpoint_probe {
                let mut probe_request = tangential_request.clone();
                probe_request.dt_s *= 0.5;
                let probe = self
                    .tangential_adapter
                    .prepare(&checkpoint.tangential_state, &probe_request)
                    .map_err(ProductionCouplingError::Tangential)?;
                (
                    probe.force_on_disc_world_n,
                    probe.free_torsional_torque_on_disc_world_nm,
                )
            } else {
                (
                    tangential.force_on_disc_world_n,
                    tangential.free_torsional_torque_on_disc_world_nm,
                )
            };

        let mut rolling_input = input.rolling.clone();
        rolling_input.patch = rolling_patch;
        rolling_input.state = checkpoint.rolling_state.generic_state.clone();
        rolling_input.checkpoint = checkpoint.rolling_state.checkpoint.clone();
        rolling_input.partial_slip_ownership = Some(tangential_request.work_ownership.clone());
        rolling_input.interval_s = input.duration_s;
        bind_rolling_kinematics_from_patch(&mut rolling_input, &patch_kinematics)?;
        let rolling = prepare_rolling_contact(&checkpoint.rolling_state, &rolling_input)
            .map_err(ProductionCouplingError::Rolling)?;

        let (gas_channel, gas_force, gas_moment) = prepare_gas_channel(
            &input.gas_channel,
            &checkpoint.gas_channel_state,
            disc_state,
            self.disc_mass_properties,
            input.duration_s,
        )?;

        let nominal_normal_force_n = point_normal_force(&normal)?;
        let normal_force_n = nominal_normal_force_n
            + surface_excitation
                .as_ref()
                .map_or(0.0, |surface| surface.normal_force_perturbation_n);
        if !(normal_force_n.is_finite() && normal_force_n >= 0.0) {
            return Err(ProductionCouplingError::InvalidInput {
                field: "topography-perturbed normal force",
            });
        }
        let (mut normal_force, mut normal_moment) = point_normal_wrench(&normal)?;
        if surface_excitation.is_some() {
            if !(nominal_normal_force_n.is_finite() && nominal_normal_force_n > 0.0) {
                return Err(ProductionCouplingError::InvalidInput {
                    field: "surface-excitation nominal normal force",
                });
            }
            // The roughness leaf is explicitly a tangent perturbation about
            // this accepted contact. Preserve the same action line while
            // scaling its point-resultant force and moment; recomputing the
            // Hertz footprint here would double-apply the linearization.
            let scale = normal_force_n / nominal_normal_force_n;
            normal_force = normal_force.scale(scale);
            normal_moment = normal_moment.scale(scale);
        }
        let tangential_moment = tangential
            .application_arm_world_m
            .cross(applied_tangential_force_world_n)
            .add(applied_tangential_free_torsional_torque_world_nm);
        let total_force_world_n = normal_force
            .add(applied_tangential_force_world_n)
            .add(rolling.step.body_wrench.contour_force_world_n)
            .add(gas_force);
        let total_moment_about_com_world_n_m = normal_moment
            .add(tangential_moment)
            .add(rolling.step.body_wrench.total_moment_about_com_world_nm)
            .add(gas_moment);
        if !(total_force_world_n.is_finite() && total_moment_about_com_world_n_m.is_finite()) {
            return Err(ProductionCouplingError::InvalidInput {
                field: "summed wrench",
            });
        }
        Ok(PreparedProductionContact {
            patch_kinematics,
            normal,
            surface_excitation,
            surface_geometry,
            tangential,
            rolling,
            gas_channel,
            normal_force_n,
            applied_tangential_force_world_n,
            applied_tangential_free_torsional_torque_world_nm,
            total_force_world_n,
            total_moment_about_com_world_n_m,
        })
    }

    fn commit_contact_step(
        &self,
        checkpoint: &ProductionCouplingCheckpoint,
        prepared: PreparedProductionContact,
        base: ProductionBaseStepProposal,
        rigid_step: StepReceipt,
    ) -> Result<(ProductionCouplingCheckpoint, ProductionCouplingReceipt), ProductionCouplingError>
    {
        let PreparedProductionContact {
            patch_kinematics,
            normal,
            surface_excitation,
            surface_geometry,
            tangential,
            rolling,
            gas_channel,
            normal_force_n: _,
            applied_tangential_force_world_n,
            applied_tangential_free_torsional_torque_world_nm,
            total_force_world_n,
            total_moment_about_com_world_n_m,
        } = prepared;
        let next_disc_state = rigid_step.state_after;
        let committed_version = checkpoint.committed_version.checked_add(1).ok_or(
            ProductionCouplingError::InvalidInput {
                field: "committed_version",
            },
        )?;
        let normal_state = normal.generic.next_state.clone();
        let tangential_state = self
            .tangential_adapter
            .commit(&checkpoint.tangential_state, &tangential)
            .map_err(ProductionCouplingError::Tangential)?;
        let rolling_state = commit_rolling_contact(&checkpoint.rolling_state, &rolling)
            .map_err(ProductionCouplingError::Rolling)?;
        let gas_channel_state = commit_gas_channel(&checkpoint.gas_channel_state, &gas_channel)?;
        let base_state = self
            .base_port
            .accept(&checkpoint.base_state, base.clone())?;
        let checkpoint_fingerprint = production_checkpoint_fingerprint(
            &self.identity,
            committed_version,
            next_disc_state,
            &normal_state,
            &tangential_state,
            &rolling_state,
            &gas_channel_state,
            &base_state,
        );
        let next = ProductionCouplingCheckpoint {
            identity: self.identity.clone(),
            committed_version,
            disc_state: next_disc_state,
            checkpoint_fingerprint,
            normal_state,
            tangential_state,
            rolling_state,
            gas_channel_state,
            base_state,
        };
        Ok((
            next,
            ProductionCouplingReceipt {
                patch_kinematics,
                normal,
                surface_excitation,
                surface_geometry,
                tangential,
                applied_tangential_force_world_n,
                applied_tangential_free_torsional_torque_world_nm,
                rolling,
                gas_channel,
                base,
                total_force_world_n,
                total_moment_about_com_world_n_m,
                next_disc_state,
                rigid_step,
                estimate_only: true,
            },
        ))
    }

    /// Advances one homogeneous open-flight substep atomically.
    ///
    /// The normal, tangential, and rolling checkpoints are retained exactly;
    /// no zero-force contact receipt is fabricated. The support is advanced by
    /// its real unforced dynamics, and the gas channel is evaluated from the
    /// same accepted rigid-body state used by fs-mbd.
    pub fn step_open_flight(
        &self,
        checkpoint: &ProductionCouplingCheckpoint,
        input: &ProductionOpenFlightStepInput,
    ) -> Result<(ProductionCouplingCheckpoint, ProductionOpenFlightReceipt), ProductionCouplingError>
    {
        self.validate_checkpoint(checkpoint)?;
        self.step_open_flight_after_checkpoint_validation(checkpoint, input)
    }

    fn step_open_flight_after_checkpoint_validation(
        &self,
        checkpoint: &ProductionCouplingCheckpoint,
        input: &ProductionOpenFlightStepInput,
    ) -> Result<(ProductionCouplingCheckpoint, ProductionOpenFlightReceipt), ProductionCouplingError>
    {
        if input.expected_checkpoint_version != checkpoint.committed_version {
            return Err(ProductionCouplingError::CheckpointVersionMismatch {
                expected: input.expected_checkpoint_version,
                observed: checkpoint.committed_version,
            });
        }
        validate_gas_step_identity(&self.identity, &input.gas_channel)?;
        if !(input.duration_s.is_finite() && input.duration_s > 0.0) {
            return Err(ProductionCouplingError::InvalidInput {
                field: "open-flight duration_s",
            });
        }
        let (gas_channel, gas_force, gas_moment) = prepare_gas_channel(
            &input.gas_channel,
            &checkpoint.gas_channel_state,
            checkpoint.disc_state,
            self.disc_mass_properties,
            input.duration_s,
        )?;
        let contact_point = checkpoint.base_state.last_contact_point_world_m();
        let base = self.base_port.propose(
            &checkpoint.base_state,
            input.base_step_id.clone(),
            input.duration_s,
            0.0,
            contact_point,
            contact_point,
            contact_point,
            input.base_load_progress_start,
            input.base_load_progress_end,
        )?;
        if !(gas_force.is_finite() && gas_moment.is_finite()) {
            return Err(ProductionCouplingError::InvalidInput {
                field: "open-flight gas wrench",
            });
        }
        let rigid_step = RigidBodyIntegrator::new(self.gravity)
            .step(
                checkpoint.disc_state,
                self.disc_mass_properties,
                Wrench {
                    force_world: gas_force,
                    torque_body: checkpoint
                        .disc_state
                        .pose()
                        .orientation()
                        .rotate_world_to_body(gas_moment),
                },
                input.duration_s,
            )
            .map_err(ProductionCouplingError::Dynamics)?;
        let next_disc_state = rigid_step.state_after;
        let committed_version = checkpoint.committed_version.checked_add(1).ok_or(
            ProductionCouplingError::InvalidInput {
                field: "committed_version",
            },
        )?;
        let gas_channel_state = commit_gas_channel(&checkpoint.gas_channel_state, &gas_channel)?;
        let base_state = self
            .base_port
            .accept(&checkpoint.base_state, base.clone())?;
        let checkpoint_fingerprint = production_checkpoint_fingerprint(
            &self.identity,
            committed_version,
            next_disc_state,
            &checkpoint.normal_state,
            &checkpoint.tangential_state,
            &checkpoint.rolling_state,
            &gas_channel_state,
            &base_state,
        );
        let next = ProductionCouplingCheckpoint {
            identity: self.identity.clone(),
            committed_version,
            disc_state: next_disc_state,
            checkpoint_fingerprint,
            normal_state: checkpoint.normal_state.clone(),
            tangential_state: checkpoint.tangential_state.clone(),
            rolling_state: checkpoint.rolling_state.clone(),
            gas_channel_state,
            base_state,
        };
        Ok((
            next,
            ProductionOpenFlightReceipt {
                gas_channel,
                base,
                total_force_world_n: gas_force,
                total_moment_about_com_world_n_m: gas_moment,
                next_disc_state,
                rigid_step,
                estimate_only: true,
            },
        ))
    }

    /// Advances one event-selected interval without retaining a trajectory.
    ///
    /// This is the streaming primitive for long simulations. It performs the
    /// same finite-footprint surface/contact resolution as
    /// [`Self::run_eventful_compliant_trajectory`], then returns exactly one
    /// branch-specific receipt. Callers may immediately reduce that receipt
    /// into force measures, acoustic cells, checkpoints, or visualization
    /// controls instead of retaining millions of heavyweight transactions.
    pub fn step_eventful_compliant(
        &self,
        checkpoint: &ProductionCouplingCheckpoint,
        input: &ProductionCouplingStepInput,
    ) -> Result<(ProductionCouplingCheckpoint, ProductionTrajectoryStep), ProductionCouplingError>
    {
        let resolved = self.resolve_production_event(checkpoint, input)?;
        let resolved_signed_gap_start_m = signed_patch_gap_m(resolved.patch_kinematics())?;
        let branch = resolved.branch();
        let start_time_s = checkpoint.elapsed_time_s();
        let (next, receipt) = match resolved {
            ResolvedProductionEvent::CompliantContact(contact) => {
                let (next, receipt) =
                    self.step_with_resolved_contact(checkpoint, input, Some(contact))?;
                (
                    next,
                    ProductionTrajectoryStepReceipt::CompliantContact(receipt),
                )
            }
            ResolvedProductionEvent::OpenFlight(_) => {
                let open = ProductionOpenFlightStepInput {
                    expected_checkpoint_version: input.expected_checkpoint_version,
                    duration_s: input.duration_s,
                    gas_channel: input.gas_channel.clone(),
                    base_step_id: input.base_step_id.clone(),
                    base_load_progress_start: input.base_load_progress_start,
                    base_load_progress_end: input.base_load_progress_end,
                };
                let (next, receipt) = self.step_open_flight(checkpoint, &open)?;
                (next, ProductionTrajectoryStepReceipt::OpenFlight(receipt))
            }
        };
        let end_time_s = next.elapsed_time_s();
        Ok((
            next,
            ProductionTrajectoryStep {
                start_time_s,
                end_time_s,
                branch,
                resolved_signed_gap_start_m,
                receipt,
            },
        ))
    }

    /// Advances one profile-bound interval from a shared predicted midpoint.
    ///
    /// The start evaluation and both half-step predictors are discard-only. All
    /// constitutive candidates are prepared again from `checkpoint` at the
    /// coupled midpoint and only those full-step candidates are committed.
    pub fn step_eventful_profile_midpoint<R>(
        &self,
        checkpoint: &ProductionCouplingCheckpoint,
        input: &ProductionCouplingStepInput,
        profile: &ResolvedDiscProfile,
        cx: &Cx<'_>,
        refresh_midpoint_input: R,
    ) -> Result<(ProductionCouplingCheckpoint, ProductionTrajectoryStep), ProductionCouplingError>
    where
        R: FnMut(
            &mut ProductionCouplingStepInput,
            RigidBodyState,
        ) -> Result<(), ProductionCouplingError>,
    {
        self.step_eventful_profile_midpoint_with(
            checkpoint,
            ProductionMidpointInput::Borrowed(input),
            profile,
            CheckpointValidation::Required,
            None,
            SurfaceTraceEvaluation::Checked,
            cx,
            refresh_midpoint_input,
        )
    }

    fn step_eventful_profile_midpoint_with<R>(
        &self,
        checkpoint: &ProductionCouplingCheckpoint,
        input: ProductionMidpointInput<'_>,
        profile: &ResolvedDiscProfile,
        checkpoint_validation: CheckpointValidation,
        admitted: Option<&AdmittedAxisymmetricProfile<'_>>,
        surface_traces: SurfaceTraceEvaluation<'_>,
        cx: &Cx<'_>,
        mut refresh_midpoint_input: R,
    ) -> Result<(ProductionCouplingCheckpoint, ProductionTrajectoryStep), ProductionCouplingError>
    where
        R: FnMut(
            &mut ProductionCouplingStepInput,
            RigidBodyState,
        ) -> Result<(), ProductionCouplingError>,
    {
        let start_input = input.as_ref();
        let resolved = match checkpoint_validation {
            CheckpointValidation::Required => {
                self.resolve_production_event(checkpoint, start_input)
            }
            CheckpointValidation::TrustedInternalSuccessor => self
                .resolve_production_event_after_checkpoint_validation(
                    checkpoint,
                    start_input,
                    surface_traces,
                ),
        }?;
        let resolved_signed_gap_start_m = signed_patch_gap_m(resolved.patch_kinematics())?;
        let start_time_s = checkpoint.elapsed_time_s();
        let ResolvedProductionEvent::CompliantContact(contact) = resolved else {
            let open = ProductionOpenFlightStepInput {
                expected_checkpoint_version: start_input.expected_checkpoint_version,
                duration_s: start_input.duration_s,
                gas_channel: start_input.gas_channel.clone(),
                base_step_id: start_input.base_step_id.clone(),
                base_load_progress_start: start_input.base_load_progress_start,
                base_load_progress_end: start_input.base_load_progress_end,
            };
            let (next, receipt) = match checkpoint_validation {
                CheckpointValidation::Required => self.step_open_flight(checkpoint, &open),
                CheckpointValidation::TrustedInternalSuccessor => {
                    self.step_open_flight_after_checkpoint_validation(checkpoint, &open)
                }
            }?;
            let end_time_s = next.elapsed_time_s();
            return Ok((
                next,
                ProductionTrajectoryStep {
                    start_time_s,
                    end_time_s,
                    branch: ProductionTrajectoryBranch::OpenFlight,
                    resolved_signed_gap_start_m,
                    receipt: ProductionTrajectoryStepReceipt::OpenFlight(receipt),
                },
            ));
        };

        let start = self.prepare_contact_channels(
            checkpoint,
            start_input,
            checkpoint.disc_state,
            contact,
            false,
        )?;
        let duration_s = start_input.duration_s;
        let base_step_id = start_input.base_step_id.clone();
        let base_load_progress_start = start_input.base_load_progress_start;
        let base_load_progress_end = start_input.base_load_progress_end;
        let half_duration_s = 0.5 * duration_s;
        let predicted_disc_state = RigidBodyIntegrator::new(self.gravity)
            .step(
                checkpoint.disc_state,
                self.disc_mass_properties,
                Wrench {
                    force_world: start.total_force_world_n,
                    torque_body: checkpoint
                        .disc_state
                        .pose()
                        .orientation()
                        .rotate_world_to_body(start.total_moment_about_com_world_n_m),
                },
                half_duration_s,
            )
            .map_err(ProductionCouplingError::Dynamics)?
            .state_after;
        let start_point = start.patch_kinematics.base_point.point_world;
        let midpoint_progress =
            base_load_progress_start + 0.5 * (base_load_progress_end - base_load_progress_start);
        let predicted_base = self.base_port.propose(
            &checkpoint.base_state,
            base_step_id.clone(),
            half_duration_s,
            start.normal_force_n,
            start_point,
            start_point,
            start_point,
            base_load_progress_start,
            midpoint_progress,
        )?;
        let predicted_base_state = self
            .base_port
            .accept(&checkpoint.base_state, predicted_base)?;

        let mut owned_midpoint_input;
        let midpoint_input = match input {
            ProductionMidpointInput::Borrowed(input) => {
                owned_midpoint_input = input.clone();
                &mut owned_midpoint_input
            }
            ProductionMidpointInput::Reusable(input) => input,
        };
        midpoint_input.time_s += half_duration_s;
        if let Some(admitted) = admitted {
            admitted.bind_contact_at_states(
                &mut *midpoint_input,
                predicted_disc_state,
                &predicted_base_state,
                cx,
            )?;
        } else {
            self.bind_horizontal_plane_axisymmetric_profile_contact_at_states(
                &mut *midpoint_input,
                profile,
                predicted_disc_state,
                &predicted_base_state,
                cx,
            )?;
        }
        refresh_midpoint_input(&mut *midpoint_input, predicted_disc_state)?;
        validate_step_identity(&self.identity, &*midpoint_input)?;
        let midpoint_contact = self.resolve_active_contact_from_patch(
            checkpoint,
            &*midpoint_input,
            midpoint_input.patch.clone(),
            surface_traces,
        )?;
        if select_event_branch(
            &midpoint_contact.patch_kinematics,
            midpoint_input.normal.integration_regime,
        )? != ProductionTrajectoryBranch::CompliantContact
        {
            return Err(ProductionCouplingError::UnsupportedMechanism {
                status: midpoint_contact.patch_kinematics.status,
            });
        }
        let midpoint = self.prepare_contact_channels(
            checkpoint,
            &*midpoint_input,
            predicted_disc_state,
            midpoint_contact,
            true,
        )?;
        let rigid_step = RigidBodyIntegrator::new(self.gravity)
            .step(
                checkpoint.disc_state,
                self.disc_mass_properties,
                Wrench {
                    force_world: midpoint.total_force_world_n,
                    torque_body: predicted_disc_state
                        .pose()
                        .orientation()
                        .rotate_world_to_body(midpoint.total_moment_about_com_world_n_m),
                },
                duration_s,
            )
            .map_err(ProductionCouplingError::Dynamics)?;
        let force_point = midpoint.patch_kinematics.base_point.point_world;
        let end_point = force_point.scale(2.0).sub(start_point);
        let base = self.base_port.propose(
            &checkpoint.base_state,
            base_step_id,
            duration_s,
            midpoint.normal_force_n,
            start_point,
            force_point,
            end_point,
            base_load_progress_start,
            base_load_progress_end,
        )?;
        let (next, receipt) = self.commit_contact_step(checkpoint, midpoint, base, rigid_step)?;
        let end_time_s = next.elapsed_time_s();
        Ok((
            next,
            ProductionTrajectoryStep {
                start_time_s,
                end_time_s,
                branch: ProductionTrajectoryBranch::CompliantContact,
                resolved_signed_gap_start_m,
                receipt: ProductionTrajectoryStepReceipt::CompliantContact(receipt),
            },
        ))
    }

    /// Runs a bounded fixed-grid open/contact compliant trajectory.
    ///
    /// The deterministic factory must rebuild geometry, material/interface
    /// cards, texture coordinates, and ownership keys from each accepted
    /// checkpoint. `Separated` selects true open flight. `Touching` and
    /// `Grazing` select finite-patch contact. `Approaching` and
    /// `ImpactCandidate` additionally require
    /// [`NormalContactIntegrationRegime::CompliantTransient`].
    /// `Unknown` refuses. Branch changes are bracketed by the accepted step,
    /// never promoted to exact impact times, and no restitution impulse is
    /// invented.
    pub fn run_eventful_compliant_trajectory<F>(
        &self,
        start_checkpoint: ProductionCouplingCheckpoint,
        maximum_accepted_steps: usize,
        mut input_for_checkpoint: F,
    ) -> ProductionEventTrajectory
    where
        F: FnMut(
            &ProductionCouplingCheckpoint,
        ) -> Result<ProductionCouplingStepInput, ProductionCouplingError>,
    {
        let mut last_accepted_checkpoint = start_checkpoint.clone();
        let mut accepted_steps = Vec::new();
        let mut transitions = Vec::new();
        let mut preceding_branch = None;
        let mut preceding_start_time_s = last_accepted_checkpoint.elapsed_time_s();

        for _ in 0..maximum_accepted_steps {
            let attempted_checkpoint_version = last_accepted_checkpoint.committed_version;
            let input = match input_for_checkpoint(&last_accepted_checkpoint) {
                Ok(input) => input,
                Err(error) => {
                    return production_event_refusal(
                        start_checkpoint,
                        last_accepted_checkpoint,
                        accepted_steps,
                        transitions,
                        attempted_checkpoint_version,
                        error,
                    );
                }
            };
            let (next_checkpoint, step) =
                match self.step_eventful_compliant(&last_accepted_checkpoint, &input) {
                    Ok(accepted) => accepted,
                    Err(error) => {
                        return production_event_refusal(
                            start_checkpoint,
                            last_accepted_checkpoint,
                            accepted_steps,
                            transitions,
                            attempted_checkpoint_version,
                            error,
                        );
                    }
                };
            if let Some(previous) = preceding_branch
                && previous != step.branch
            {
                transitions.push(ProductionBranchTransition {
                    from: previous,
                    to: step.branch,
                    bracket_start_s: preceding_start_time_s,
                    bracket_end_s: step.start_time_s,
                });
            }
            preceding_branch = Some(step.branch);
            preceding_start_time_s = step.start_time_s;
            accepted_steps.push(step);
            last_accepted_checkpoint = next_checkpoint;
        }

        // Classify the final accepted endpoint as well as every interval
        // start. Otherwise a branch crossing during the last permitted step
        // would disappear merely because no subsequent step was requested.
        // This calls the same deterministic factory once more but advances no
        // constituent and commits no ownership key. The dry finite-footprint
        // evaluation still requires capacity for its next logical normal key.
        if let Some(previous) = preceding_branch {
            let attempted_checkpoint_version = last_accepted_checkpoint.committed_version;
            let input = match input_for_checkpoint(&last_accepted_checkpoint) {
                Ok(input) => input,
                Err(error) => {
                    return production_event_refusal(
                        start_checkpoint,
                        last_accepted_checkpoint,
                        accepted_steps,
                        transitions,
                        attempted_checkpoint_version,
                        error,
                    );
                }
            };
            let terminal_branch =
                match self.resolve_production_event(&last_accepted_checkpoint, &input) {
                    Ok(resolved) => resolved.branch(),
                    Err(error) => {
                        return production_event_refusal(
                            start_checkpoint,
                            last_accepted_checkpoint,
                            accepted_steps,
                            transitions,
                            attempted_checkpoint_version,
                            error,
                        );
                    }
                };
            if previous != terminal_branch {
                transitions.push(ProductionBranchTransition {
                    from: previous,
                    to: terminal_branch,
                    bracket_start_s: preceding_start_time_s,
                    bracket_end_s: last_accepted_checkpoint.elapsed_time_s(),
                });
            }
        }
        ProductionEventTrajectory {
            start_checkpoint,
            last_accepted_checkpoint,
            accepted_steps,
            transitions,
            termination: ProductionEventTrajectoryTermination::StepLimitReached {
                maximum_accepted_steps,
            },
        }
    }

    /// Runs an eventful trajectory while retaining only homogeneous reduced
    /// control intervals and the final coupled checkpoint.
    ///
    /// `mechanics_steps_per_control_interval` is a reduction factor, not a
    /// change to the physics timestep. Normal force is accumulated as impulse,
    /// discrete disc-work residuals are summed, and every branch change flushes
    /// the current interval before the new branch is admitted. Consequently an
    /// 8x reduction from 384 kHz mechanics to 48 kHz controls preserves force
    /// measure exactly without retaining 3.07 million full receipts.
    pub fn run_eventful_control_trajectory<F>(
        &self,
        start_checkpoint: ProductionCouplingCheckpoint,
        maximum_accepted_steps: usize,
        mechanics_steps_per_control_interval: usize,
        input_for_checkpoint: F,
    ) -> Result<ProductionControlTrajectory, ProductionCouplingError>
    where
        F: FnMut(
            &ProductionCouplingCheckpoint,
        ) -> Result<ProductionCouplingStepInput, ProductionCouplingError>,
    {
        self.run_eventful_control_trajectory_observed(
            start_checkpoint,
            maximum_accepted_steps,
            mechanics_steps_per_control_interval,
            input_for_checkpoint,
            |_| {},
        )
    }

    /// Runs the reduced-control trajectory while observing each accepted
    /// rectangular-modal mechanics step before its information is reduced.
    ///
    /// The observer is never called for a refused proposal or for the dry
    /// terminal branch classification. It receives only borrowed data and must
    /// reduce it during the call rather than archive every mechanics receipt.
    pub fn run_eventful_control_trajectory_observed<F, O>(
        &self,
        start_checkpoint: ProductionCouplingCheckpoint,
        maximum_accepted_steps: usize,
        mechanics_steps_per_control_interval: usize,
        mut input_for_checkpoint: F,
        observe_modal_audio_step: O,
    ) -> Result<ProductionControlTrajectory, ProductionCouplingError>
    where
        F: FnMut(
            &ProductionCouplingCheckpoint,
        ) -> Result<ProductionCouplingStepInput, ProductionCouplingError>,
        O: for<'step> FnMut(ProductionModalAudioStep<'step>),
    {
        self.run_eventful_control_trajectory_observed_with(
            start_checkpoint,
            maximum_accepted_steps,
            mechanics_steps_per_control_interval,
            CheckpointValidation::Required,
            SurfaceTraceEvaluation::Checked,
            |checkpoint, reusable_input| {
                *reusable_input = Some(input_for_checkpoint(checkpoint)?);
                Ok(())
            },
            |checkpoint, input| self.step_eventful_compliant(checkpoint, input),
            observe_modal_audio_step,
        )
    }

    /// Runs profile-native controls with one shared contact/rigid/base midpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn run_eventful_profile_midpoint_control_trajectory_observed<F, R, O>(
        &self,
        start_checkpoint: ProductionCouplingCheckpoint,
        maximum_accepted_steps: usize,
        mechanics_steps_per_control_interval: usize,
        profile: &ResolvedDiscProfile,
        cx: &Cx<'_>,
        mut input_for_checkpoint: F,
        mut refresh_midpoint_input: R,
        observe_modal_audio_step: O,
    ) -> Result<ProductionControlTrajectory, ProductionCouplingError>
    where
        F: FnMut(
            &ProductionCouplingCheckpoint,
        ) -> Result<ProductionCouplingStepInput, ProductionCouplingError>,
        R: FnMut(
            &mut ProductionCouplingStepInput,
            RigidBodyState,
        ) -> Result<(), ProductionCouplingError>,
        O: for<'step> FnMut(ProductionModalAudioStep<'step>),
    {
        self.run_eventful_control_trajectory_observed_with(
            start_checkpoint,
            maximum_accepted_steps,
            mechanics_steps_per_control_interval,
            CheckpointValidation::Required,
            SurfaceTraceEvaluation::Checked,
            |checkpoint, reusable_input| {
                *reusable_input = Some(input_for_checkpoint(checkpoint)?);
                Ok(())
            },
            |checkpoint, input| {
                self.step_eventful_profile_midpoint(
                    checkpoint,
                    input,
                    profile,
                    cx,
                    |input, state| refresh_midpoint_input(input, state),
                )
            },
            observe_modal_audio_step,
        )
    }

    fn run_eventful_control_trajectory_observed_with<F, S, O>(
        &self,
        start_checkpoint: ProductionCouplingCheckpoint,
        maximum_accepted_steps: usize,
        mechanics_steps_per_control_interval: usize,
        checkpoint_validation: CheckpointValidation,
        surface_traces: SurfaceTraceEvaluation<'_>,
        mut input_for_checkpoint: F,
        mut advance: S,
        mut observe_modal_audio_step: O,
    ) -> Result<ProductionControlTrajectory, ProductionCouplingError>
    where
        F: FnMut(
            &ProductionCouplingCheckpoint,
            &mut Option<ProductionCouplingStepInput>,
        ) -> Result<(), ProductionCouplingError>,
        S: FnMut(
            &ProductionCouplingCheckpoint,
            &mut ProductionCouplingStepInput,
        ) -> Result<
            (ProductionCouplingCheckpoint, ProductionTrajectoryStep),
            ProductionCouplingError,
        >,
        O: for<'step> FnMut(ProductionModalAudioStep<'step>),
    {
        if mechanics_steps_per_control_interval == 0 {
            return Err(ProductionCouplingError::InvalidInput {
                field: "mechanics_steps_per_control_interval",
            });
        }
        self.validate_checkpoint(&start_checkpoint)?;
        let mut intervals = Vec::new();
        let requested_intervals = maximum_accepted_steps
            .checked_add(mechanics_steps_per_control_interval - 1)
            .map(|value| value / mechanics_steps_per_control_interval)
            .ok_or(ProductionCouplingError::InvalidInput {
                field: "control trajectory capacity",
            })?;
        intervals
            .try_reserve_exact(requested_intervals)
            .map_err(|_| ProductionCouplingError::InvalidInput {
                field: "control trajectory capacity",
            })?;
        let mut transitions = Vec::new();
        let mut accumulator = None;
        let mut last_accepted_checkpoint = start_checkpoint.clone();
        let mut accepted_mechanics_steps = 0_usize;
        let mut preceding_branch = None;
        let mut preceding_start_time_s = start_checkpoint.elapsed_time_s();
        let mut reusable_input = None;
        let termination = loop {
            if accepted_mechanics_steps == maximum_accepted_steps {
                break ProductionEventTrajectoryTermination::StepLimitReached {
                    maximum_accepted_steps,
                };
            }
            let attempted_checkpoint_version = last_accepted_checkpoint.committed_version;
            if let Err(error) = input_for_checkpoint(&last_accepted_checkpoint, &mut reusable_input)
            {
                break ProductionEventTrajectoryTermination::Refused {
                    attempted_checkpoint_version,
                    error,
                };
            }
            let Some(input) = reusable_input.as_mut() else {
                break ProductionEventTrajectoryTermination::Refused {
                    attempted_checkpoint_version,
                    error: ProductionCouplingError::InvalidInput {
                        field: "control trajectory input provider",
                    },
                };
            };
            let (next_checkpoint, step) = match advance(&last_accepted_checkpoint, input) {
                Ok(accepted) => accepted,
                Err(error) => {
                    break ProductionEventTrajectoryTermination::Refused {
                        attempted_checkpoint_version,
                        error,
                    };
                }
            };
            if let Some(modal_audio_step) = step.modal_audio_step() {
                observe_modal_audio_step(modal_audio_step);
            }
            if let Some(active) = &mut accumulator {
                set_control_accumulator_endpoint_gap(active, step.resolved_signed_gap_start_m)?;
            }
            let reduction_boundary = accumulator.as_ref().is_some_and(|active| {
                active.mechanics_substeps == mechanics_steps_per_control_interval
            });
            if let Some(previous) = preceding_branch
                && previous != step.branch
            {
                flush_control_accumulator(&mut intervals, &mut accumulator)?;
                transitions.push(ProductionBranchTransition {
                    from: previous,
                    to: step.branch,
                    bracket_start_s: preceding_start_time_s,
                    bracket_end_s: step.start_time_s,
                });
            } else if reduction_boundary {
                flush_control_accumulator(&mut intervals, &mut accumulator)?;
            }
            match &mut accumulator {
                Some(active) => extend_control_accumulator(self, active, &step)?,
                None => accumulator = Some(start_control_accumulator(self, &step)?),
            }
            accepted_mechanics_steps = accepted_mechanics_steps.checked_add(1).ok_or(
                ProductionCouplingError::InvalidInput {
                    field: "accepted mechanics step count",
                },
            )?;
            preceding_branch = Some(step.branch);
            preceding_start_time_s = step.start_time_s;
            last_accepted_checkpoint = next_checkpoint;
        };
        // Preserve the ordinary event driver's endpoint classification without
        // advancing any constituent or committing its dry-evaluation work key.
        // The supplied normal ledger must nevertheless admit the next logical
        // key so the finite-footprint classifier can evaluate that endpoint.
        if let Some(previous) = preceding_branch {
            let attempted_checkpoint_version = last_accepted_checkpoint.committed_version;
            let endpoint = input_for_checkpoint(&last_accepted_checkpoint, &mut reusable_input)
                .and_then(|()| {
                    let input =
                        reusable_input
                            .as_ref()
                            .ok_or(ProductionCouplingError::InvalidInput {
                                field: "control trajectory input provider",
                            })?;
                    let resolved = match checkpoint_validation {
                        CheckpointValidation::Required => {
                            self.resolve_production_event(&last_accepted_checkpoint, input)
                        }
                        CheckpointValidation::TrustedInternalSuccessor => self
                            .resolve_production_event_after_checkpoint_validation(
                                &last_accepted_checkpoint,
                                input,
                                surface_traces,
                            ),
                    };
                    resolved.and_then(|resolved| {
                        Ok((
                            resolved.branch(),
                            signed_patch_gap_m(resolved.patch_kinematics())?,
                        ))
                    })
                });
            match endpoint {
                Ok((terminal_branch, resolved_signed_gap_end_m)) => {
                    if let Some(active) = &mut accumulator {
                        set_control_accumulator_endpoint_gap(active, resolved_signed_gap_end_m)?;
                    }
                    if previous != terminal_branch {
                        transitions.push(ProductionBranchTransition {
                            from: previous,
                            to: terminal_branch,
                            bracket_start_s: preceding_start_time_s,
                            bracket_end_s: last_accepted_checkpoint.elapsed_time_s(),
                        });
                    }
                }
                Err(error)
                    if matches!(
                        termination,
                        ProductionEventTrajectoryTermination::StepLimitReached { .. }
                    ) =>
                {
                    // Endpoint-classification refusals are part of the
                    // physical prefix contract and must not be hidden on a
                    // nominally completed horizon. A prior in-loop refusal is
                    // retained verbatim because it is the first failed state.
                    flush_control_accumulator(&mut intervals, &mut accumulator)?;
                    return Ok(ProductionControlTrajectory {
                        start_checkpoint,
                        last_accepted_checkpoint,
                        intervals,
                        accepted_mechanics_steps,
                        transitions,
                        termination: ProductionEventTrajectoryTermination::Refused {
                            attempted_checkpoint_version,
                            error,
                        },
                    });
                }
                Err(_) => {}
            }
        }
        flush_control_accumulator(&mut intervals, &mut accumulator)?;
        Ok(ProductionControlTrajectory {
            start_checkpoint,
            last_accepted_checkpoint,
            intervals,
            accepted_mechanics_steps,
            transitions,
            termination,
        })
    }

    fn resolve_patch_kinematics(
        &self,
        checkpoint: &ProductionCouplingCheckpoint,
        input: &ProductionCouplingStepInput,
    ) -> Result<PatchKinematics, ProductionCouplingError> {
        compute_moving_one_mode_patch_kinematics(self.prepared_patch_input(checkpoint, input)?)
            .map_err(ProductionCouplingError::Patch)
    }

    fn resolve_active_contact(
        &self,
        checkpoint: &ProductionCouplingCheckpoint,
        input: &ProductionCouplingStepInput,
    ) -> Result<ResolvedProductionContact, ProductionCouplingError> {
        self.validate_checkpoint(checkpoint)?;
        self.resolve_active_contact_after_checkpoint_validation(
            checkpoint,
            input,
            SurfaceTraceEvaluation::Checked,
        )
    }

    fn resolve_active_contact_after_checkpoint_validation(
        &self,
        checkpoint: &ProductionCouplingCheckpoint,
        input: &ProductionCouplingStepInput,
        surface_traces: SurfaceTraceEvaluation<'_>,
    ) -> Result<ResolvedProductionContact, ProductionCouplingError> {
        let patch_input =
            self.prepared_patch_input_after_checkpoint_validation(checkpoint, input)?;
        self.resolve_active_contact_from_patch(checkpoint, input, patch_input, surface_traces)
    }

    fn resolve_active_contact_from_patch(
        &self,
        checkpoint: &ProductionCouplingCheckpoint,
        input: &ProductionCouplingStepInput,
        patch_input: MovingOneModePatchKinematicsInput,
        surface_traces: SurfaceTraceEvaluation<'_>,
    ) -> Result<ResolvedProductionContact, ProductionCouplingError> {
        if input.surface_excitation.is_some() && input.surface_geometry.is_some() {
            return Err(ProductionCouplingError::InvalidInput {
                field: "mutually exclusive surface coupling modes",
            });
        }
        let mut normal_input = input.normal.clone();
        normal_input.state = checkpoint.normal_state.clone();
        normal_input.time_s = input.time_s;
        normal_input.step_s = input.duration_s;
        if let Some(surface) = &input.surface_geometry {
            let (patch_kinematics, normal, surface_geometry) = resolve_surface_geometry_contact(
                patch_input,
                &normal_input,
                surface,
                surface_traces,
            )?;
            Ok(ResolvedProductionContact {
                patch_kinematics,
                normal,
                surface_geometry: Some(surface_geometry),
            })
        } else {
            let patch_kinematics = compute_moving_one_mode_patch_kinematics(patch_input)
                .map_err(ProductionCouplingError::Patch)?;
            normal_input.kinematics = patch_kinematics.clone();
            let normal = match evaluate_normal_contact(&normal_input)
                .map_err(ProductionCouplingError::Normal)?
            {
                EulerNormalContactOutcome::Active(active) => active,
                EulerNormalContactOutcome::InactiveSeparated { .. } => {
                    return Err(ProductionCouplingError::UnsupportedMechanism {
                        status: PatchContactStatus::Separated,
                    });
                }
            };
            Ok(ResolvedProductionContact {
                patch_kinematics,
                normal,
                surface_geometry: None,
            })
        }
    }

    fn resolve_production_event(
        &self,
        checkpoint: &ProductionCouplingCheckpoint,
        input: &ProductionCouplingStepInput,
    ) -> Result<ResolvedProductionEvent, ProductionCouplingError> {
        self.validate_checkpoint(checkpoint)?;
        self.resolve_production_event_after_checkpoint_validation(
            checkpoint,
            input,
            SurfaceTraceEvaluation::Checked,
        )
    }

    fn resolve_production_event_after_checkpoint_validation(
        &self,
        checkpoint: &ProductionCouplingCheckpoint,
        input: &ProductionCouplingStepInput,
        surface_traces: SurfaceTraceEvaluation<'_>,
    ) -> Result<ResolvedProductionEvent, ProductionCouplingError> {
        match self.resolve_active_contact_after_checkpoint_validation(
            checkpoint,
            input,
            surface_traces,
        ) {
            Ok(contact) => {
                let branch = select_event_branch(
                    &contact.patch_kinematics,
                    input.normal.integration_regime,
                )?;
                match branch {
                    ProductionTrajectoryBranch::CompliantContact => {
                        Ok(ResolvedProductionEvent::CompliantContact(contact))
                    }
                    // A dry unilateral patch may still have a well-defined
                    // zero-load footprint at the exact separation boundary.
                    // If that same resolved patch is receding, carry its
                    // kinematics into the existing open-flight step rather
                    // than rejecting a physically valid contact-to-flight
                    // handoff or resolving the rough surface a second time.
                    ProductionTrajectoryBranch::OpenFlight => Ok(
                        ResolvedProductionEvent::OpenFlight(contact.patch_kinematics),
                    ),
                }
            }
            Err(
                separated @ ProductionCouplingError::UnsupportedMechanism {
                    status: PatchContactStatus::Separated,
                },
            ) => {
                // Without an active footprint, the exact unilateral event
                // boundary is the point-support limit. It may select true open
                // flight, but it cannot override a contact classification when
                // the finite-footprint solve itself found no admissible patch.
                let patch = self.resolve_event_patch_kinematics_after_checkpoint_validation(
                    checkpoint,
                    input,
                    surface_traces,
                )?;
                match select_event_branch(&patch, input.normal.integration_regime)? {
                    ProductionTrajectoryBranch::OpenFlight => {
                        Ok(ResolvedProductionEvent::OpenFlight(patch))
                    }
                    ProductionTrajectoryBranch::CompliantContact => Err(separated),
                }
            }
            Err(error) => Err(error),
        }
    }

    fn resolve_event_patch_kinematics_after_checkpoint_validation(
        &self,
        checkpoint: &ProductionCouplingCheckpoint,
        input: &ProductionCouplingStepInput,
        surface_traces: SurfaceTraceEvaluation<'_>,
    ) -> Result<PatchKinematics, ProductionCouplingError> {
        let patch_input =
            self.prepared_patch_input_after_checkpoint_validation(checkpoint, input)?;
        let Some(surface) = &input.surface_geometry else {
            return compute_moving_one_mode_patch_kinematics(patch_input)
                .map_err(ProductionCouplingError::Patch);
        };
        if input.surface_excitation.is_some() {
            return Err(ProductionCouplingError::InvalidInput {
                field: "mutually exclusive surface coupling modes",
            });
        }
        validate_surface_geometry_identity(&input.normal, surface)?;
        let point = surface_traces
            .evaluate_point_surface_pair(
                &surface.interface,
                surface.surface_a.as_motion(),
                surface.surface_b.as_motion(),
            )
            .map_err(ProductionSurfaceExcitationError::Surface)
            .map_err(ProductionCouplingError::SurfaceExcitation)?;
        compute_surface_adjusted_patch(patch_input, &point)
    }

    fn prepared_patch_input(
        &self,
        checkpoint: &ProductionCouplingCheckpoint,
        input: &ProductionCouplingStepInput,
    ) -> Result<MovingOneModePatchKinematicsInput, ProductionCouplingError> {
        self.validate_checkpoint(checkpoint)?;
        self.prepared_patch_input_after_checkpoint_validation(checkpoint, input)
    }

    fn prepared_patch_input_after_checkpoint_validation(
        &self,
        checkpoint: &ProductionCouplingCheckpoint,
        input: &ProductionCouplingStepInput,
    ) -> Result<MovingOneModePatchKinematicsInput, ProductionCouplingError> {
        if input.expected_checkpoint_version != checkpoint.committed_version {
            return Err(ProductionCouplingError::CheckpointVersionMismatch {
                expected: input.expected_checkpoint_version,
                observed: checkpoint.committed_version,
            });
        }
        validate_step_identity(&self.identity, input)?;
        if !(input.duration_s.is_finite()
            && input.duration_s > 0.0
            && input.time_s.is_finite()
            && input.time_s >= 0.0)
        {
            return Err(ProductionCouplingError::InvalidInput {
                field: "duration_s or time_s",
            });
        }
        let mut patch_input = input.patch.clone();
        patch_input.bridge.disc_state = checkpoint.disc_state;
        let (base_displacement_m, base_velocity_m_per_s) = self.base_port.surface_state(
            &checkpoint.base_state,
            patch_input.bridge.profile_support.disc_point_world_m,
        )?;
        patch_input.bridge.base_mode.vertical_displacement_m = base_displacement_m;
        patch_input.bridge.base_mode.vertical_velocity_m_per_s = base_velocity_m_per_s;
        Ok(patch_input)
    }

    /// Runs at most `maximum_accepted_steps` smooth homogeneous substeps.
    ///
    /// This is deliberately a trajectory *prefix* helper, not an event owner:
    /// when a rebuilt patch leaves the touching/grazing envelope, [`Self::step`]
    /// returns its typed refusal and this method preserves the preceding
    /// checkpoint for an impact, separation, thin-gap, or terminal-event owner.
    /// The caller must supply a deterministic factory to make replay or a
    /// prefix/resume comparison meaningful.
    pub fn run_smooth_contact_trajectory<F>(
        &self,
        start_checkpoint: ProductionCouplingCheckpoint,
        maximum_accepted_steps: usize,
        mut input_for_checkpoint: F,
    ) -> SmoothContactTrajectory
    where
        F: FnMut(
            &ProductionCouplingCheckpoint,
        ) -> Result<ProductionCouplingStepInput, ProductionCouplingError>,
    {
        let mut last_accepted_checkpoint = start_checkpoint.clone();
        // Do not reserve an untrusted caller limit up front: bounded execution
        // must not turn a large limit into an immediate allocation failure.
        let mut accepted_steps = Vec::new();
        for _ in 0..maximum_accepted_steps {
            let attempted_checkpoint_version = last_accepted_checkpoint.committed_version;
            let input = match input_for_checkpoint(&last_accepted_checkpoint) {
                Ok(input) => input,
                Err(error) => {
                    return SmoothContactTrajectory {
                        start_checkpoint,
                        last_accepted_checkpoint,
                        accepted_steps,
                        termination: SmoothContactTrajectoryTermination::Refused {
                            attempted_checkpoint_version,
                            error,
                        },
                    };
                }
            };
            match self.step(&last_accepted_checkpoint, &input) {
                Ok((next_checkpoint, receipt)) => {
                    last_accepted_checkpoint = next_checkpoint;
                    accepted_steps.push(receipt);
                }
                Err(error) => {
                    return SmoothContactTrajectory {
                        start_checkpoint,
                        last_accepted_checkpoint,
                        accepted_steps,
                        termination: SmoothContactTrajectoryTermination::Refused {
                            attempted_checkpoint_version,
                            error,
                        },
                    };
                }
            }
        }
        SmoothContactTrajectory {
            start_checkpoint,
            last_accepted_checkpoint,
            accepted_steps,
            termination: SmoothContactTrajectoryTermination::StepLimitReached {
                maximum_accepted_steps,
            },
        }
    }
}

fn bind_rolling_kinematics_from_patch(
    rolling: &mut RollingContactInput,
    patch: &PatchKinematics,
) -> Result<(), ProductionCouplingError> {
    let normal = patch.tangent_basis.normal_world;
    let rotational_surface_velocity = patch
        .disc_point
        .angular_velocity_world
        .cross(patch.disc_point.arm_world);
    let projected =
        rotational_surface_velocity.sub(normal.scale(rotational_surface_velocity.dot(normal)));
    let projected_norm = projected.norm_squared().sqrt();
    let contour_axis = if projected_norm > 0.0 && projected_norm.is_finite() {
        projected.scale(projected_norm.recip())
    } else {
        patch.tangent_basis.first_world
    };
    let rolling_axis_raw = normal.cross(contour_axis);
    let rolling_axis_norm = rolling_axis_raw.norm_squared().sqrt();
    if !(contour_axis.is_finite() && rolling_axis_norm.is_finite() && rolling_axis_norm > 0.0) {
        return Err(ProductionCouplingError::InvalidInput {
            field: "derived rolling frame",
        });
    }
    let rolling_axis = rolling_axis_raw.scale(rolling_axis_norm.recip());
    rolling.contour_tangent_axis_world = contour_axis;
    rolling.rolling_axis_world = rolling_axis;
    rolling.contour_speed_mps = rotational_surface_velocity.dot(contour_axis);
    rolling.rolling_rate_rad_s = patch.disc_point.angular_velocity_world.dot(rolling_axis);
    rolling.spin_rate_rad_s = patch.normal_spin_rad_per_s;
    if !(rolling.contour_speed_mps.is_finite()
        && rolling.rolling_rate_rad_s.is_finite()
        && rolling.spin_rate_rad_s.is_finite())
    {
        return Err(ProductionCouplingError::InvalidInput {
            field: "derived rolling rates",
        });
    }
    Ok(())
}

fn select_event_branch(
    patch: &PatchKinematics,
    integration_regime: NormalContactIntegrationRegime,
) -> Result<ProductionTrajectoryBranch, ProductionCouplingError> {
    let status = patch.status;
    let signed_gap_m = patch
        .disc_point
        .point_world
        .sub(patch.base_point.point_world)
        .dot(patch.tangent_basis.normal_world);
    if !signed_gap_m.is_finite() {
        return Err(ProductionCouplingError::InvalidInput {
            field: "event-branch signed gap",
        });
    }
    // This is a dry, non-adhesive unilateral interface. A certainly positive
    // gap has no contact port, even while the bodies are approaching, and an
    // exactly closed interface that is receding opens without evaluating a
    // zero-indentation Hertz tangent. Negative gap remains on the compliant
    // branch until the accepted state actually reaches separation.
    if signed_gap_m > 0.0 || (signed_gap_m == 0.0 && status == PatchContactStatus::Receding) {
        return Ok(ProductionTrajectoryBranch::OpenFlight);
    }
    match status {
        PatchContactStatus::Separated => Ok(ProductionTrajectoryBranch::OpenFlight),
        PatchContactStatus::Touching | PatchContactStatus::Grazing => {
            Ok(ProductionTrajectoryBranch::CompliantContact)
        }
        PatchContactStatus::Approaching
        | PatchContactStatus::Receding
        | PatchContactStatus::ImpactCandidate
            if integration_regime == NormalContactIntegrationRegime::CompliantTransient =>
        {
            Ok(ProductionTrajectoryBranch::CompliantContact)
        }
        status => Err(ProductionCouplingError::UnsupportedMechanism { status }),
    }
}

struct ProductionControlAccumulator {
    start_time_s: f64,
    end_time_s: f64,
    mechanics_substeps: usize,
    branch: ProductionTrajectoryBranch,
    resolved_signed_gap_start_m: f64,
    resolved_signed_gap_end_m: Option<f64>,
    state_before: RigidBodyState,
    state_after: RigidBodyState,
    normal_impulse_n_s: f64,
    base_displacement_start_m: f64,
    base_displacement_end_m: f64,
    base_velocity_start_m_per_s: f64,
    base_velocity_end_m_per_s: f64,
    base_parent_version: u64,
    base_next_version: u64,
    disc_work_residual_j: f64,
    mechanical_energy_end_j: f64,
}

fn production_step_observables<'a>(
    step: &'a ProductionTrajectoryStep,
) -> Result<
    (
        &'a StepReceipt,
        &'a ProductionBaseStepReceipt,
        f64,
        Vec3,
        Vec3,
    ),
    ProductionCouplingError,
> {
    match (&step.branch, &step.receipt) {
        (
            ProductionTrajectoryBranch::CompliantContact,
            ProductionTrajectoryStepReceipt::CompliantContact(receipt),
        ) => Ok((
            &receipt.rigid_step,
            receipt.base.receipt(),
            receipt.base.receipt().compressive_normal_force_on_base_n,
            receipt.total_force_world_n,
            receipt.total_moment_about_com_world_n_m,
        )),
        (
            ProductionTrajectoryBranch::OpenFlight,
            ProductionTrajectoryStepReceipt::OpenFlight(receipt),
        ) => Ok((
            &receipt.rigid_step,
            receipt.base.receipt(),
            0.0,
            receipt.total_force_world_n,
            receipt.total_moment_about_com_world_n_m,
        )),
        _ => Err(ProductionCouplingError::InvalidInput {
            field: "event branch/receipt mismatch",
        }),
    }
}

fn production_disc_work_residual_j(
    model: &ProductionCouplingModel,
    rigid_step: &StepReceipt,
    force_world_n: Vec3,
    moment_about_com_world_n_m: Vec3,
) -> Result<f64, ProductionCouplingError> {
    let before = rigid_step.state_before;
    let after = rigid_step.state_after;
    let velocity_before = before
        .center_of_mass_velocity_world(model.disc_mass_properties)
        .map_err(ProductionCouplingError::Dynamics)?;
    let velocity_after = after
        .center_of_mass_velocity_world(model.disc_mass_properties)
        .map_err(ProductionCouplingError::Dynamics)?;
    let omega_before_body = model
        .disc_mass_properties
        .angular_velocity_body_checked(before.angular_momentum_body())
        .map_err(ProductionCouplingError::Dynamics)?;
    let omega_after_body = model
        .disc_mass_properties
        .angular_velocity_body_checked(after.angular_momentum_body())
        .map_err(ProductionCouplingError::Dynamics)?;
    let velocity_mid = velocity_before.add(velocity_after).scale(0.5);
    let omega_mid_body = omega_before_body.add(omega_after_body).scale(0.5);
    let torque_body = before
        .pose()
        .orientation()
        .rotate_world_to_body(moment_about_com_world_n_m);
    let wrench_work_j = rigid_step.duration_seconds
        * (force_world_n.dot(velocity_mid) + torque_body.dot(omega_mid_body));
    let residual = rigid_step.diagnostics_after.mechanical_energy
        - rigid_step.diagnostics_before.mechanical_energy
        - wrench_work_j;
    if residual.is_finite() {
        Ok(residual)
    } else {
        Err(ProductionCouplingError::InvalidInput {
            field: "production disc work residual",
        })
    }
}

fn start_control_accumulator(
    model: &ProductionCouplingModel,
    step: &ProductionTrajectoryStep,
) -> Result<ProductionControlAccumulator, ProductionCouplingError> {
    let (rigid, base, force_n, force, moment) = production_step_observables(step)?;
    Ok(ProductionControlAccumulator {
        start_time_s: step.start_time_s,
        end_time_s: step.end_time_s,
        mechanics_substeps: 1,
        branch: step.branch,
        resolved_signed_gap_start_m: step.resolved_signed_gap_start_m,
        resolved_signed_gap_end_m: None,
        state_before: rigid.state_before,
        state_after: rigid.state_after,
        normal_impulse_n_s: force_n * rigid.duration_seconds,
        base_displacement_start_m: base.modal_displacement_start_m,
        base_displacement_end_m: base.modal_displacement_end_m,
        base_velocity_start_m_per_s: base.modal_velocity_start_m_per_s,
        base_velocity_end_m_per_s: base.modal_velocity_end_m_per_s,
        base_parent_version: base.parent_version,
        base_next_version: base.next_version,
        disc_work_residual_j: production_disc_work_residual_j(model, rigid, force, moment)?,
        mechanical_energy_end_j: rigid.diagnostics_after.mechanical_energy,
    })
}

fn extend_control_accumulator(
    model: &ProductionCouplingModel,
    accumulator: &mut ProductionControlAccumulator,
    step: &ProductionTrajectoryStep,
) -> Result<(), ProductionCouplingError> {
    let (rigid, base, force_n, force, moment) = production_step_observables(step)?;
    if accumulator.branch != step.branch
        || accumulator.end_time_s.to_bits() != step.start_time_s.to_bits()
        || accumulator.state_after != rigid.state_before
        || accumulator.base_next_version != base.parent_version
        || accumulator
            .resolved_signed_gap_end_m
            .is_none_or(|gap_m| gap_m.to_bits() != step.resolved_signed_gap_start_m.to_bits())
    {
        return Err(ProductionCouplingError::InvalidInput {
            field: "control-interval lineage",
        });
    }
    accumulator.end_time_s = step.end_time_s;
    accumulator.mechanics_substeps = accumulator.mechanics_substeps.checked_add(1).ok_or(
        ProductionCouplingError::InvalidInput {
            field: "control-interval substep count",
        },
    )?;
    accumulator.state_after = rigid.state_after;
    accumulator.normal_impulse_n_s += force_n * rigid.duration_seconds;
    accumulator.base_displacement_end_m = base.modal_displacement_end_m;
    accumulator.base_velocity_end_m_per_s = base.modal_velocity_end_m_per_s;
    accumulator.base_next_version = base.next_version;
    accumulator.disc_work_residual_j +=
        production_disc_work_residual_j(model, rigid, force, moment)?;
    accumulator.mechanical_energy_end_j = rigid.diagnostics_after.mechanical_energy;
    accumulator.resolved_signed_gap_end_m = None;
    Ok(())
}

fn set_control_accumulator_endpoint_gap(
    accumulator: &mut ProductionControlAccumulator,
    resolved_signed_gap_end_m: f64,
) -> Result<(), ProductionCouplingError> {
    if !resolved_signed_gap_end_m.is_finite() || accumulator.resolved_signed_gap_end_m.is_some() {
        return Err(ProductionCouplingError::InvalidInput {
            field: "control-interval resolved endpoint gap",
        });
    }
    accumulator.resolved_signed_gap_end_m = Some(resolved_signed_gap_end_m);
    Ok(())
}

fn finish_control_accumulator(
    accumulator: ProductionControlAccumulator,
) -> Result<ProductionControlInterval, ProductionCouplingError> {
    let duration_s = accumulator.end_time_s - accumulator.start_time_s;
    let mean_normal_force_n = accumulator.normal_impulse_n_s / duration_s;
    if !(duration_s.is_finite()
        && duration_s > 0.0
        && mean_normal_force_n.is_finite()
        && mean_normal_force_n >= 0.0
        && accumulator.resolved_signed_gap_start_m.is_finite()
        && accumulator
            .resolved_signed_gap_end_m
            .is_none_or(f64::is_finite)
        && accumulator.disc_work_residual_j.is_finite())
    {
        return Err(ProductionCouplingError::InvalidInput {
            field: "reduced control interval",
        });
    }
    Ok(ProductionControlInterval {
        start_time_s: accumulator.start_time_s,
        end_time_s: accumulator.end_time_s,
        mechanics_substeps: accumulator.mechanics_substeps,
        branch: accumulator.branch,
        resolved_signed_gap_start_m: accumulator.resolved_signed_gap_start_m,
        resolved_signed_gap_end_m: accumulator.resolved_signed_gap_end_m,
        state_before: accumulator.state_before,
        state_after: accumulator.state_after,
        normal_impulse_n_s: accumulator.normal_impulse_n_s,
        mean_normal_force_n,
        base_displacement_start_m: accumulator.base_displacement_start_m,
        base_displacement_end_m: accumulator.base_displacement_end_m,
        base_velocity_start_m_per_s: accumulator.base_velocity_start_m_per_s,
        base_velocity_end_m_per_s: accumulator.base_velocity_end_m_per_s,
        base_parent_version: accumulator.base_parent_version,
        base_next_version: accumulator.base_next_version,
        disc_work_residual_j: accumulator.disc_work_residual_j,
        mechanical_energy_end_j: accumulator.mechanical_energy_end_j,
    })
}

fn flush_control_accumulator(
    intervals: &mut Vec<ProductionControlInterval>,
    accumulator: &mut Option<ProductionControlAccumulator>,
) -> Result<(), ProductionCouplingError> {
    if let Some(active) = accumulator.take() {
        intervals.push(finish_control_accumulator(active)?);
    }
    Ok(())
}

fn production_event_refusal(
    start_checkpoint: ProductionCouplingCheckpoint,
    last_accepted_checkpoint: ProductionCouplingCheckpoint,
    accepted_steps: Vec<ProductionTrajectoryStep>,
    transitions: Vec<ProductionBranchTransition>,
    attempted_checkpoint_version: u64,
    error: ProductionCouplingError,
) -> ProductionEventTrajectory {
    ProductionEventTrajectory {
        start_checkpoint,
        last_accepted_checkpoint,
        accepted_steps,
        transitions,
        termination: ProductionEventTrajectoryTermination::Refused {
            attempted_checkpoint_version,
            error,
        },
    }
}

fn normal_patch_view(
    normal: &ActiveNormalContact,
    patch: &PatchKinematics,
) -> Result<NormalPatchView, ProductionCouplingError> {
    let NormalPatchReceipt::Point(receipt) = &normal.generic.receipt else {
        return Err(ProductionCouplingError::UnsupportedLineNormalContact);
    };
    let (longitudinal, lateral) = receipt
        .elliptic_patch_axes
        .map_or((receipt.patch_radius_m, receipt.patch_radius_m), |axes| {
            (axes.semi_major_axis_m, axes.semi_minor_axis_m)
        });
    NormalPatchView::new(
        patch.patch.patch_identity.as_str(),
        normal.material_card_id.clone(),
        normal.material_source_id.clone(),
        normal_authority(receipt.authority),
        receipt.normal_force_n,
        longitudinal,
        lateral,
        receipt.pressure.second_moment_m2,
    )
    .map_err(|_| ProductionCouplingError::InvalidInput {
        field: "normal patch view",
    })
}

fn rolling_patch_receipt(
    normal: &ActiveNormalContact,
    patch: &PatchKinematics,
) -> Result<RollingPatchReceipt, ProductionCouplingError> {
    let NormalPatchReceipt::Point(receipt) = &normal.generic.receipt else {
        return Err(ProductionCouplingError::UnsupportedLineNormalContact);
    };
    let CurvatureMetadata::Known {
        authority,
        first_principal_m_inverse,
        second_principal_m_inverse,
        ..
    } = &patch.patch.curvature
    else {
        return Err(ProductionCouplingError::CurvatureUnavailable);
    };
    RollingPatchReceipt::new(
        patch.patch.patch_identity.as_str(),
        normal.material_card_id.clone(),
        normal.material_source_id.clone(),
        *authority,
        receipt.normal_force_n,
        receipt.patch_radius_m,
        PatchCurvature::Principal {
            first_per_m: *first_principal_m_inverse,
            second_per_m: *second_principal_m_inverse,
        },
    )
    .map_err(|_| ProductionCouplingError::InvalidInput {
        field: "rolling patch receipt",
    })
}

fn normal_authority(authority: InputAuthority) -> NormalPatchAuthority {
    match authority {
        InputAuthority::CallerDeclared => NormalPatchAuthority::CallerDeclared,
        InputAuthority::SyntheticFixture => NormalPatchAuthority::SyntheticFixture,
        InputAuthority::Estimated => NormalPatchAuthority::Estimated,
    }
}

fn point_normal_force(normal: &ActiveNormalContact) -> Result<f64, ProductionCouplingError> {
    let NormalPatchReceipt::Point(receipt) = &normal.generic.receipt else {
        return Err(ProductionCouplingError::UnsupportedLineNormalContact);
    };
    Ok(receipt.normal_force_n)
}

fn observed_static_normal_force(
    input: &ProductionCouplingStepInput,
    normal: &ActiveNormalContact,
) -> Result<f64, ProductionCouplingError> {
    let surface_perturbation_n = input
        .surface_excitation
        .as_ref()
        .map(|surface| {
            evaluate_surface_excitation_for_normal(
                normal,
                &surface.interface,
                surface.surface_a.as_motion(),
                surface.surface_b.as_motion(),
                surface.travel_angle_from_patch_major_rad,
                surface.maximum_linearized_height_fraction,
            )
            .map(|receipt| receipt.normal_force_perturbation_n)
        })
        .transpose()
        .map_err(ProductionCouplingError::SurfaceExcitation)?
        .unwrap_or(0.0);
    Ok(point_normal_force(normal)? + surface_perturbation_n)
}

fn point_normal_wrench(
    normal: &ActiveNormalContact,
) -> Result<(Vec3, Vec3), ProductionCouplingError> {
    let NormalPatchPort::Point(port) = &normal.generic.port else {
        return Err(ProductionCouplingError::UnsupportedLineNormalContact);
    };
    Ok((
        Vec3::new(
            port.action_force_n[0],
            port.action_force_n[1],
            port.action_force_n[2],
        ),
        Vec3::new(
            port.action_moment_n_m[0],
            port.action_moment_n_m[1],
            port.action_moment_n_m[2],
        ),
    ))
}

fn validate_identity(identity: &ProductionCouplingIdentity) -> Result<(), ProductionCouplingError> {
    for value in [
        &identity.case_id,
        &identity.configuration_id,
        &identity.world_frame_id,
    ] {
        if value.trim().is_empty() || value.len() > 256 || !value.is_ascii() {
            return Err(ProductionCouplingError::InvalidInput { field: "identity" });
        }
    }
    Ok(())
}

/// Canonical enough state binding for this crate's wholly in-memory checkpoint.
///
/// Every nested checkpoint type supplies a deterministic derived `Debug`
/// representation (including ordered/BTree retained work keys). The field is
/// intentionally recomputed before every step, so a caller cannot modify the
/// public mechanics state or version while retaining private channel snapshots.
fn production_checkpoint_fingerprint(
    identity: &ProductionCouplingIdentity,
    committed_version: u64,
    disc_state: RigidBodyState,
    normal_state: &NormalPatchEmbedState,
    tangential_state: &TangentialContactState,
    rolling_state: &RollingContactState,
    gas_channel_state: &GasChannelState,
    base_state: &ProductionBaseCheckpoint,
) -> ContentHash {
    let mut hasher = DomainHasher::new("fs-euler-disc-e2e/production-coupling-checkpoint/v1");
    write!(
        hasher,
        "{identity:?}|{committed_version}|{disc_state:?}|{normal_state:?}|{tangential_state:?}|{rolling_state:?}|{gas_channel_state:?}|{base_state:?}"
    )
    .expect("DomainHasher's fmt::Write implementation is infallible");
    hasher.finalize()
}

fn validate_step_identity(
    identity: &ProductionCouplingIdentity,
    input: &ProductionCouplingStepInput,
) -> Result<(), ProductionCouplingError> {
    for (actual, expected, field) in [
        (
            input.normal.identity.case_id.as_str(),
            identity.case_id.as_str(),
            "normal.case_id",
        ),
        (
            input.rolling.identity.case_id.as_str(),
            identity.case_id.as_str(),
            "rolling.case_id",
        ),
        (
            input.rolling.identity.world_frame_id.as_str(),
            identity.world_frame_id.as_str(),
            "rolling.world_frame_id",
        ),
    ] {
        if actual != expected {
            return Err(ProductionCouplingError::InputIdentityMismatch { field });
        }
    }
    validate_gas_step_identity(identity, &input.gas_channel)
}

fn validate_gas_step_identity(
    identity: &ProductionCouplingIdentity,
    gas_channel: &GasChannelStepInput,
) -> Result<(), ProductionCouplingError> {
    match gas_channel {
        GasChannelStepInput::ExteriorFreeGas { input, .. } => {
            for (actual, expected, field) in [
                (
                    input.identity.case_id.as_str(),
                    identity.case_id.as_str(),
                    "exterior_air.case_id",
                ),
                (
                    input.identity.world_frame_id.as_str(),
                    identity.world_frame_id.as_str(),
                    "exterior_air.world_frame_id",
                ),
            ] {
                if actual != expected {
                    return Err(ProductionCouplingError::InputIdentityMismatch { field });
                }
            }
        }
        GasChannelStepInput::ThinGap { input, .. } => {
            for (actual, expected, field) in [
                (
                    input.identity.case_id.as_str(),
                    identity.case_id.as_str(),
                    "air_film.case_id",
                ),
                (
                    input.identity.configuration_id.as_str(),
                    identity.configuration_id.as_str(),
                    "air_film.configuration_id",
                ),
                (
                    input.identity.frame_id.as_str(),
                    identity.world_frame_id.as_str(),
                    "air_film.frame_id",
                ),
            ] {
                if actual != expected {
                    return Err(ProductionCouplingError::InputIdentityMismatch { field });
                }
            }
        }
    }
    Ok(())
}

fn validate_gas_state_identity(
    identity: &ProductionCouplingIdentity,
    state: &GasChannelState,
) -> Result<(), ProductionCouplingError> {
    let GasChannelState::ThinGap(state) = state else {
        return Ok(());
    };
    let gas_identity = state.identity();
    for (actual, expected, field) in [
        (
            gas_identity.case_id.as_str(),
            identity.case_id.as_str(),
            "air_film_state.case_id",
        ),
        (
            gas_identity.configuration_id.as_str(),
            identity.configuration_id.as_str(),
            "air_film_state.configuration_id",
        ),
        (
            gas_identity.frame_id.as_str(),
            identity.world_frame_id.as_str(),
            "air_film_state.frame_id",
        ),
    ] {
        if actual != expected {
            return Err(ProductionCouplingError::InputIdentityMismatch { field });
        }
    }
    Ok(())
}

fn prepare_gas_channel(
    input: &GasChannelStepInput,
    state: &GasChannelState,
    disc_state: RigidBodyState,
    properties: MassProperties,
    duration_s: f64,
) -> Result<(GasChannelReceipt, Vec3, Vec3), ProductionCouplingError> {
    match (input, state) {
        (
            GasChannelStepInput::ExteriorFreeGas {
                input,
                selected_correlation_id,
                exchange_key,
            },
            GasChannelState::ExteriorFreeGas(state),
        ) => {
            let mut exterior_air_input = input.clone();
            exterior_air_input.state = exterior_state_from_disc(disc_state, properties)
                .map_err(ProductionCouplingError::Dynamics)?;
            let exterior_set = evaluate_euler_disc_external_air(&exterior_air_input)
                .map_err(ProductionCouplingError::ExteriorAir)?;
            let candidate = exterior_set
                .candidates
                .into_iter()
                .find(|candidate| candidate.world_wrench.correlation.id == *selected_correlation_id)
                .ok_or(ProductionCouplingError::ExteriorCandidateUnavailable)?;
            let force = flux_vec3(candidate.world_wrench.force_world_n);
            let moment = flux_vec3(candidate.world_wrench.torque_world_n_m);
            let work = state
                .prepare(*exchange_key, duration_s, &candidate)
                .map_err(ProductionCouplingError::ExteriorAir)?;
            Ok((
                GasChannelReceipt::ExteriorFreeGas { candidate, work },
                force,
                moment,
            ))
        }
        (
            GasChannelStepInput::ThinGap {
                input,
                exchange_key,
            },
            GasChannelState::ThinGap(state),
        ) => {
            let mut air_film_input = input.clone();
            air_film_input.disc = air_film_disc_state_from_disc(disc_state, properties)
                .map_err(ProductionCouplingError::Dynamics)?;
            air_film_input.timestep_s = duration_s;
            let proposal = state
                .prepare(*exchange_key, &air_film_input)
                .map_err(ProductionCouplingError::AirFilm)?;
            let wrench = proposal.step.receipt.wrench;
            Ok((
                GasChannelReceipt::ThinGap { proposal },
                air_vec3(wrench.force_world_n),
                air_vec3(wrench.moment_about_com_world_n_m),
            ))
        }
        _ => Err(ProductionCouplingError::GasChannelMismatch),
    }
}

fn commit_gas_channel(
    state: &GasChannelState,
    receipt: &GasChannelReceipt,
) -> Result<GasChannelState, ProductionCouplingError> {
    match (state, receipt) {
        (
            GasChannelState::ExteriorFreeGas(state),
            GasChannelReceipt::ExteriorFreeGas { work, .. },
        ) => state
            .commit(work)
            .map(GasChannelState::ExteriorFreeGas)
            .map_err(ProductionCouplingError::ExteriorAir),
        (GasChannelState::ThinGap(state), GasChannelReceipt::ThinGap { proposal }) => state
            .commit(proposal)
            .map(GasChannelState::ThinGap)
            .map_err(ProductionCouplingError::AirFilm),
        _ => Err(ProductionCouplingError::GasChannelMismatch),
    }
}

/// Derives the exterior-air state from the checkpoint that receives its wrench.
///
/// The exterior card supplies geometry, gas, and correlation alternatives, but
/// cannot smuggle a different pose or rate into a mechanically homogeneous
/// substep.
fn exterior_state_from_disc(
    state: RigidBodyState,
    properties: MassProperties,
) -> Result<EulerDiscExteriorState, fs_mbd::DynamicsError> {
    let pose = state.pose();
    let orientation = pose.orientation();
    let angular_velocity_body =
        properties.angular_velocity_body_checked(state.angular_momentum_body())?;
    Ok(EulerDiscExteriorState {
        center_world_m: mbd_to_flux_vec3(pose.position_world()),
        center_velocity_world_m_per_s: mbd_to_flux_vec3(
            state.center_of_mass_velocity_world(properties)?,
        ),
        angular_velocity_world_rad_per_s: mbd_to_flux_vec3(
            orientation.rotate_body_to_world(angular_velocity_body),
        ),
        body_frame: EulerDiscBodyFrame {
            x_world: mbd_to_flux_vec3(orientation.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0))),
            z_world: mbd_to_flux_vec3(orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0))),
        },
    })
}

/// Derives the gas-film disc kinematics from the same fs-mbd checkpoint that
/// receives the gas-on-body wrench. `AirFilmWrench` is documented in world
/// axes about COM, so no sign or frame conversion is inferred at composition.
fn air_film_disc_state_from_disc(
    state: RigidBodyState,
    properties: MassProperties,
) -> Result<TiltedDiscKinematics, fs_mbd::DynamicsError> {
    let pose = state.pose();
    let orientation = pose.orientation();
    let angular_velocity_body =
        properties.angular_velocity_body_checked(state.angular_momentum_body())?;
    Ok(TiltedDiscKinematics {
        center_world_m: mbd_to_air_vec3(pose.position_world()),
        normal_away_from_base_world: mbd_to_air_vec3(
            orientation.rotate_body_to_world(Vec3::new(0.0, 0.0, 1.0)),
        ),
        center_velocity_world_m_per_s: mbd_to_air_vec3(
            state.center_of_mass_velocity_world(properties)?,
        ),
        angular_velocity_world_rad_per_s: mbd_to_air_vec3(
            orientation.rotate_body_to_world(angular_velocity_body),
        ),
    })
}

fn flux_vec3(value: fs_flux::Vec3) -> Vec3 {
    Vec3::new(value.x, value.y, value.z)
}

fn mbd_to_flux_vec3(value: Vec3) -> fs_flux::Vec3 {
    fs_flux::Vec3::new(value.x, value.y, value.z)
}

fn air_vec3(value: AirVec3) -> Vec3 {
    Vec3::new(value.x, value.y, value.z)
}

fn mbd_to_air_vec3(value: Vec3) -> AirVec3 {
    AirVec3::new(value.x, value.y, value.z)
}
