//! One atomic, mechanically homogeneous Euler-disc accepted substep.
//!
//! This is an orchestration boundary over existing adapters.  It is deliberately
//! Estimate-only: it does not calibrate a disc, rank air correlations, resolve
//! impact, supply thin-gap pressure, or claim a resolved/as-built base.

use core::fmt;

use fs_blake3::{ContentHash, hash_domain};
use fs_contact::normal_patch::{NormalPatchEmbedState, NormalPatchPort, NormalPatchReceipt};
use fs_couple::StableId;
use fs_exec::Cx;
use fs_mbd::{
    Gravity, MassProperties, RigidBodyIntegrator, RigidBodyState, StepReceipt, Vec3, Wrench,
};
use fs_rep_frep::{AxisymmetricCurvatureAuthority, AxisymmetricIdentity};
use fs_tribo::{
    InputAuthority, InterfaceSystemRef,
    partial_slip::{NormalPatchAuthority, NormalPatchView},
    rolling_loss::{PatchCurvature, RollingPatchReceipt},
    surface_excitation::{
        HertzRoughnessExcitationInput, HertzRoughnessExcitationReceipt, ProjectedHertzFootprint,
        SurfaceExcitationError, SurfaceTraceMotion, evaluate_hertz_roughness_excitation,
    },
};

use crate::{
    air::{
        AirFilmError, AirFilmProposal, AirFilmTransactionState, AirVec3, TiltedDiscAirFilmInput,
        TiltedDiscKinematics,
    },
    base_response::{
        ReducedBaseCheckpoint, ReducedBasePort, ReducedBaseStepInput, ReducedBaseStepProposal,
    },
    contact_dynamics::{
        ContactDynamicsError, ProfileContactPatchGeometry,
        profile_contact_patch_geometry_from_mass, profile_mass_to_mbd,
    },
    external_air::{
        EulerDiscBodyFrame, EulerDiscExteriorState, EulerExternalAirCandidate,
        EulerExternalAirInput, EulerExternalAirWorkProposal, EulerExternalAirWorkState,
        ExternalAirError, evaluate_euler_disc_external_air,
    },
    normal_contact::{
        ActiveNormalContact, EulerNormalContactInput, EulerNormalContactOutcome,
        NormalContactError, evaluate_normal_contact,
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

/// Immutable adapters and rigid-body properties used by a production substep.
#[derive(Debug, Clone)]
pub struct ProductionCouplingModel {
    /// Identity binding checkpoints to this model.
    pub identity: ProductionCouplingIdentity,
    /// Disc mass and principal inertia.
    pub disc_mass_properties: MassProperties,
    /// Uniform world-frame gravity.
    pub gravity: Gravity,
    /// Already assembled moving-one-mode flexible-base port.
    pub base_port: ReducedBasePort,
    /// Explicitly selected partial-slip adapter/lane.
    pub tangential_adapter: EulerTangentialContactAdapter,
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
    base_state: ReducedBaseCheckpoint,
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
        self.base_state.modal_displacement_m()
    }

    /// Accepted one-mode base velocity [m/s].
    #[must_use]
    pub fn base_velocity_m_per_s(&self) -> f64 {
        self.base_state.modal_velocity_m_per_s()
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
    /// Prepared and accepted tangential response.
    pub tangential: TangentialContactReceipt,
    /// Prepared and accepted rolling response.
    pub rolling: RollingContactProposal,
    /// The one gas-channel receipt contributing its real wrench to fs-mbd.
    pub gas_channel: GasChannelReceipt,
    /// Accepted moving-one-mode base transition accounting.
    pub base: ReducedBaseStepProposal,
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
        let NormalPatchReceipt::Point(normal) = &self.normal.generic.receipt else {
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
}

impl fmt::Display for ProductionSurfaceExcitationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for ProductionSurfaceExcitationError {}

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
    Base(crate::base_response::BaseResponseError),
    /// fs-mbd refused the actual summed wrench.
    Dynamics(fs_mbd::DynamicsError),
}

impl fmt::Display for ProductionCouplingError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for ProductionCouplingError {}

impl ProductionCouplingModel {
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
        let (resolved, disc_mass_properties) =
            self.resolve_axisymmetric_profile_contact(profile, disc_state, cx)?;
        let resolved = bind_horizontal_plane_axisymmetric_profile_contact_input(
            input,
            resolved,
            disc_mass_properties,
            disc_state,
            base_state.modal_displacement_m(),
            base_state.modal_velocity_m_per_s(),
        )?;
        input.expected_checkpoint_version = 0;
        Ok(resolved)
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
        let (resolved, disc_mass_properties) =
            self.resolve_axisymmetric_profile_contact(profile, checkpoint.disc_state, cx)?;
        let resolved = bind_horizontal_plane_axisymmetric_profile_contact_input(
            input,
            resolved,
            disc_mass_properties,
            checkpoint.disc_state,
            checkpoint.base_state.modal_displacement_m(),
            checkpoint.base_state.modal_velocity_m_per_s(),
        )?;
        input.expected_checkpoint_version = checkpoint.committed_version;
        Ok(resolved)
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
        validate_identity(&self.identity)?;
        if !disc_state.pose().position_world().is_finite() {
            return Err(ProductionCouplingError::InvalidInput {
                field: "disc_state",
            });
        }
        validate_gas_state_identity(&self.identity, &gas_channel_state)?;
        let base_state = self.base_port.initial_checkpoint();
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
        let mut patch_input = input.patch.clone();
        patch_input.bridge.disc_state = checkpoint.disc_state;
        patch_input.bridge.base_mode.vertical_displacement_m =
            checkpoint.base_state.modal_displacement_m();
        patch_input.bridge.base_mode.vertical_velocity_m_per_s =
            checkpoint.base_state.modal_velocity_m_per_s();
        let patch_kinematics = compute_moving_one_mode_patch_kinematics(patch_input)
            .map_err(ProductionCouplingError::Patch)?;
        if matches!(
            patch_kinematics.patch.curvature,
            CurvatureMetadata::Unavailable { .. }
        ) {
            return Err(ProductionCouplingError::CurvatureUnavailable);
        }
        if !matches!(
            patch_kinematics.status,
            PatchContactStatus::Touching | PatchContactStatus::Grazing
        ) {
            return Err(ProductionCouplingError::UnsupportedMechanism {
                status: patch_kinematics.status,
            });
        }

        let mut normal_input = input.normal.clone();
        normal_input.kinematics = patch_kinematics.clone();
        normal_input.state = checkpoint.normal_state.clone();
        normal_input.time_s = input.time_s;
        normal_input.step_s = input.duration_s;
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

        let mut rolling_input = input.rolling.clone();
        rolling_input.patch = rolling_patch;
        rolling_input.state = checkpoint.rolling_state.generic_state.clone();
        rolling_input.checkpoint = checkpoint.rolling_state.checkpoint.clone();
        rolling_input.partial_slip_ownership = Some(tangential_request.work_ownership.clone());
        rolling_input.interval_s = input.duration_s;
        let rolling = prepare_rolling_contact(&checkpoint.rolling_state, &rolling_input)
            .map_err(ProductionCouplingError::Rolling)?;

        let (gas_channel, gas_force, gas_moment) = prepare_gas_channel(
            &input.gas_channel,
            &checkpoint.gas_channel_state,
            checkpoint.disc_state,
            self.disc_mass_properties,
            input.duration_s,
        )?;

        let normal_force_n = point_normal_force(&normal)?;
        let base = self
            .base_port
            .propose(
                &checkpoint.base_state,
                &ReducedBaseStepInput {
                    step_id: input.base_step_id.clone(),
                    expected_version: checkpoint.base_state.accepted_version(),
                    duration_s: input.duration_s,
                    compressive_normal_force_on_base_n: normal_force_n,
                    load_progress_start: input.base_load_progress_start,
                    load_progress_end: input.base_load_progress_end,
                },
            )
            .map_err(ProductionCouplingError::Base)?;

        let (normal_force, normal_moment) = point_normal_wrench(&normal)?;
        let tangential_moment = tangential
            .application_arm_world_m
            .cross(tangential.force_on_disc_world_n)
            .add(tangential.free_torsional_torque_on_disc_world_nm);
        let total_force_world_n = normal_force
            .add(tangential.force_on_disc_world_n)
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
        let rigid_step = RigidBodyIntegrator::new(self.gravity)
            .step(
                checkpoint.disc_state,
                self.disc_mass_properties,
                Wrench {
                    force_world: total_force_world_n,
                    torque_body: checkpoint
                        .disc_state
                        .pose()
                        .orientation()
                        .rotate_world_to_body(total_moment_about_com_world_n_m),
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
            .accept(&checkpoint.base_state, base.clone())
            .map_err(ProductionCouplingError::Base)?;
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
                tangential,
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
    base_state: &ReducedBaseCheckpoint,
) -> ContentHash {
    hash_domain(
        "fs-euler-disc-e2e/production-coupling-checkpoint/v1",
        format!(
            "{identity:?}|{committed_version}|{disc_state:?}|{normal_state:?}|{tangential_state:?}|{rolling_state:?}|{gas_channel_state:?}|{base_state:?}"
        )
        .as_bytes(),
    )
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
    match &input.gas_channel {
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
