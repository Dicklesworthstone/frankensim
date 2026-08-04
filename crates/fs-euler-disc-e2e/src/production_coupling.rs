//! One atomic, mechanically homogeneous Euler-disc accepted substep.
//!
//! This is an orchestration boundary over existing adapters.  It is deliberately
//! Estimate-only: it does not calibrate a disc, rank air correlations, resolve
//! impact, supply thin-gap pressure, or claim a resolved/as-built base.

use core::fmt;

use fs_contact::normal_patch::{NormalPatchEmbedState, NormalPatchPort, NormalPatchReceipt};
use fs_mbd::{Gravity, MassProperties, RigidBodyIntegrator, RigidBodyState, Vec3, Wrench};
use fs_tribo::{
    InputAuthority,
    partial_slip::{NormalPatchAuthority, NormalPatchView},
    rolling_loss::{PatchCurvature, RollingPatchReceipt},
};

use crate::{
    base_response::{
        ReducedBaseCheckpoint, ReducedBasePort, ReducedBaseStepInput, ReducedBaseStepProposal,
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
        PatchKinematicsError, compute_moving_one_mode_patch_kinematics,
    },
    rolling_contact::{
        RollingContactError, RollingContactInput, RollingContactProposal, RollingContactState,
        commit_rolling_contact, prepare_rolling_contact,
    },
    tangential_contact::{
        EulerTangentialContactAdapter, TangentialContactError, TangentialContactReceipt,
        TangentialContactRequest, TangentialContactState,
    },
};

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
    normal_state: NormalPatchEmbedState,
    tangential_state: TangentialContactState,
    rolling_state: RollingContactState,
    exterior_air_state: EulerExternalAirWorkState,
    base_state: ReducedBaseCheckpoint,
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
    /// Exterior free-gas alternatives. The named candidate is selected exactly, never ranked.
    pub exterior_air: EulerExternalAirInput,
    /// Explicit correlation identity selected by the caller from `exterior_air.alternatives`.
    pub selected_exterior_correlation_id: String,
    /// Exactly-once exterior work key for this accepted interval.
    pub exterior_exchange_key: u64,
    /// Base-port replay identity for this accepted interval.
    pub base_step_id: String,
    /// Moving-load location at interval start.
    pub base_load_progress_start: f64,
    /// Moving-load location at interval end.
    pub base_load_progress_end: f64,
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
    /// Explicitly selected exterior candidate, retained without ranking alternatives.
    pub exterior_air: EulerExternalAirCandidate,
    /// Exactly-once exterior work staged/accepted with this mechanics step.
    pub exterior_air_work: EulerExternalAirWorkProposal,
    /// Accepted moving-one-mode base transition accounting.
    pub base: ReducedBaseStepProposal,
    /// Total real world-frame force sent to fs-mbd [N].
    pub total_force_world_n: Vec3,
    /// Total real world-frame moment about disc COM sent to fs-mbd [N m].
    pub total_moment_about_com_world_n_m: Vec3,
    /// fs-mbd accepted disc state.
    pub next_disc_state: RigidBodyState,
    /// This composition retains only source-adapter Estimate authority.
    pub estimate_only: bool,
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

/// Typed refusal; no caller checkpoint is changed on every error path.
#[derive(Debug, Clone, PartialEq)]
pub enum ProductionCouplingError {
    /// Model/checkpoint identity or version mismatch.
    CheckpointMismatch,
    /// Caller prepared cards from a different accepted outer checkpoint.
    CheckpointVersionMismatch { expected: u64, observed: u64 },
    /// A caller-owned adapter card is not bound to this case or world frame.
    InputIdentityMismatch { field: &'static str },
    /// Invalid outer scalar or identity.
    InvalidInput { field: &'static str },
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
    /// Exterior-air admission/refusal; thin-gap input is explicitly rejected downstream.
    ExteriorAir(ExternalAirError),
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
    /// Validates immutable identities and creates one complete initial checkpoint.
    pub fn initial_checkpoint(
        &self,
        disc_state: RigidBodyState,
        normal_state: NormalPatchEmbedState,
        tangential_state: TangentialContactState,
        rolling_state: RollingContactState,
        exterior_air_state: EulerExternalAirWorkState,
    ) -> Result<ProductionCouplingCheckpoint, ProductionCouplingError> {
        validate_identity(&self.identity)?;
        if !disc_state.pose().position_world().is_finite() {
            return Err(ProductionCouplingError::InvalidInput {
                field: "disc_state",
            });
        }
        Ok(ProductionCouplingCheckpoint {
            identity: self.identity.clone(),
            committed_version: 0,
            disc_state,
            normal_state,
            tangential_state,
            rolling_state,
            exterior_air_state,
            base_state: self.base_port.initial_checkpoint(),
        })
    }

    /// Attempts one homogeneous smooth substep and atomically advances every channel only after fs-mbd accepts.
    pub fn step(
        &self,
        checkpoint: &ProductionCouplingCheckpoint,
        input: &ProductionCouplingStepInput,
    ) -> Result<(ProductionCouplingCheckpoint, ProductionCouplingReceipt), ProductionCouplingError>
    {
        validate_identity(&self.identity)?;
        if checkpoint.identity != self.identity {
            return Err(ProductionCouplingError::CheckpointMismatch);
        }
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

        let mut exterior_air_input = input.exterior_air.clone();
        exterior_air_input.state =
            exterior_state_from_disc(checkpoint.disc_state, self.disc_mass_properties)
                .map_err(ProductionCouplingError::Dynamics)?;
        let exterior_set = evaluate_euler_disc_external_air(&exterior_air_input)
            .map_err(ProductionCouplingError::ExteriorAir)?;
        let exterior_air = exterior_set
            .candidates
            .into_iter()
            .find(|candidate| {
                candidate.world_wrench.correlation.id == input.selected_exterior_correlation_id
            })
            .ok_or(ProductionCouplingError::ExteriorCandidateUnavailable)?;
        let exterior_air_work = checkpoint
            .exterior_air_state
            .prepare(input.exterior_exchange_key, input.duration_s, &exterior_air)
            .map_err(ProductionCouplingError::ExteriorAir)?;

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
            .add(flux_vec3(exterior_air.world_wrench.force_world_n));
        let total_moment_about_com_world_n_m = normal_moment
            .add(tangential_moment)
            .add(rolling.step.body_wrench.total_moment_about_com_world_nm)
            .add(flux_vec3(exterior_air.world_wrench.torque_world_n_m));
        if !(total_force_world_n.is_finite() && total_moment_about_com_world_n_m.is_finite()) {
            return Err(ProductionCouplingError::InvalidInput {
                field: "summed wrench",
            });
        }
        let next_disc_state = RigidBodyIntegrator::new(self.gravity)
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
            .map_err(ProductionCouplingError::Dynamics)?
            .state_after;

        let next = ProductionCouplingCheckpoint {
            identity: self.identity.clone(),
            committed_version: checkpoint.committed_version.checked_add(1).ok_or(
                ProductionCouplingError::InvalidInput {
                    field: "committed_version",
                },
            )?,
            disc_state: next_disc_state,
            normal_state: normal.generic.next_state.clone(),
            tangential_state: self
                .tangential_adapter
                .commit(&checkpoint.tangential_state, &tangential)
                .map_err(ProductionCouplingError::Tangential)?,
            rolling_state: commit_rolling_contact(&checkpoint.rolling_state, &rolling)
                .map_err(ProductionCouplingError::Rolling)?,
            exterior_air_state: checkpoint
                .exterior_air_state
                .commit(&exterior_air_work)
                .map_err(ProductionCouplingError::ExteriorAir)?,
            base_state: self
                .base_port
                .accept(&checkpoint.base_state, base.clone())
                .map_err(ProductionCouplingError::Base)?,
        };
        Ok((
            next,
            ProductionCouplingReceipt {
                patch_kinematics,
                normal,
                tangential,
                rolling,
                exterior_air,
                exterior_air_work,
                base,
                total_force_world_n,
                total_moment_about_com_world_n_m,
                next_disc_state,
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
        InputAuthority::Estimated,
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
        (
            input.exterior_air.identity.case_id.as_str(),
            identity.case_id.as_str(),
            "exterior_air.case_id",
        ),
        (
            input.exterior_air.identity.world_frame_id.as_str(),
            identity.world_frame_id.as_str(),
            "exterior_air.world_frame_id",
        ),
    ] {
        if actual != expected {
            return Err(ProductionCouplingError::InputIdentityMismatch { field });
        }
    }
    Ok(())
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

fn flux_vec3(value: fs_flux::Vec3) -> Vec3 {
    Vec3::new(value.x, value.y, value.z)
}

fn mbd_to_flux_vec3(value: Vec3) -> fs_flux::Vec3 {
    fs_flux::Vec3::new(value.x, value.y, value.z)
}
