//! Euler-disc translation of finite-patch kinematics into tangential contact.
//!
//! This adapter owns no friction coefficients or contact-state classification.
//! It refuses unavailable creepage and pre-constitutive patch states, then
//! delegates the actual return map to `fs-tribo` (or the `fs-contact` smooth
//! transaction wrapper).  The emitted receipt keeps the generic discrete work
//! closure distinct from the disc endpoint-power diagnostic.

use core::fmt;

use fs_contact::tangential::smooth::{
    SmoothTangentialAdapter, SmoothTangentialError, SmoothTangentialQuery, SmoothTangentialState,
};
use fs_mbd::Vec3;
use fs_tribo::partial_slip::{
    GeneralizedWorkOwnership, NormalPatchView, PartialSlipCheckpoint, PartialSlipError,
    PartialSlipInterface, PartialSlipKinematics, PartialSlipLaw, PartialSlipState,
    PartialSlipStateKind, TangentFrame, TangentialWrench,
};

use crate::patch_kinematics::{Creepage, PatchContactStatus, PatchKinematics, SurfaceOrder};

/// Stable model identity for the Euler translation layer.
pub const EULER_TANGENTIAL_CONTACT_ADAPTER_MODEL_ID: &str =
    "fs-euler-disc-e2e/tangential-contact-adapter-v1";

/// The explicitly selected generic tangential lane.
#[derive(Clone, Debug, PartialEq)]
pub enum TangentialContactLane {
    /// Direct nonsmooth partial-slip return map.
    PartialSlip { law: PartialSlipLaw },
    /// The generic smooth transaction wrapper over the same partial-slip law.
    Smooth {
        law: PartialSlipLaw,
        adapter: SmoothTangentialAdapter,
    },
}

/// Immutable Euler adapter configuration.
#[derive(Clone, Debug, PartialEq)]
pub struct EulerTangentialContactAdapter {
    adapter_id: String,
    source_id: String,
    lane: TangentialContactLane,
}

impl EulerTangentialContactAdapter {
    /// Creates a named adapter with an explicitly selected generic lane.
    pub fn new(
        adapter_id: impl Into<String>,
        source_id: impl Into<String>,
        lane: TangentialContactLane,
    ) -> Result<Self, TangentialContactError> {
        let value = Self {
            adapter_id: adapter_id.into(),
            source_id: source_id.into(),
            lane,
        };
        value.validate_identity()?;
        Ok(value)
    }

    /// Stable adapter identity retained by every checkpoint.
    #[must_use]
    pub fn adapter_id(&self) -> &str {
        &self.adapter_id
    }

    /// Source/configuration identity retained by every checkpoint.
    #[must_use]
    pub fn source_id(&self) -> &str {
        &self.source_id
    }

    /// Selected direct or smooth generic lane.
    #[must_use]
    pub const fn lane(&self) -> &TangentialContactLane {
        &self.lane
    }

    /// Creates a zero-history state bound to a particular finite normal patch
    /// and ordered interface.  A zero or absent work-key capacity refuses.
    pub fn initial_state(
        &self,
        normal_patch: &NormalPatchView,
        interface: &PartialSlipInterface,
        max_committed_work_keys: usize,
    ) -> Result<TangentialContactState, TangentialContactError> {
        self.validate_identity()?;
        if max_committed_work_keys == 0 {
            return Err(TangentialContactError::InvalidInput {
                field: "max_committed_work_keys",
            });
        }
        match &self.lane {
            TangentialContactLane::PartialSlip { law } => Ok(TangentialContactState::PartialSlip {
                adapter_id: self.adapter_id.clone(),
                source_id: self.source_id.clone(),
                checkpoint: PartialSlipCheckpoint::new(
                    normal_patch.clone(),
                    interface.clone(),
                    law.clone(),
                    PartialSlipState::zero(),
                )?,
                committed_version: 0,
                committed_work_keys: Vec::new(),
                max_committed_work_keys,
            }),
            TangentialContactLane::Smooth { law, adapter } => Ok(TangentialContactState::Smooth {
                adapter_id: self.adapter_id.clone(),
                source_id: self.source_id.clone(),
                state: adapter.initial_state(
                    law,
                    normal_patch,
                    interface,
                    max_committed_work_keys,
                )?,
            }),
        }
    }

    /// Restores a checkpoint only when the adapter, selected lane, normal
    /// patch, interface, and generic upstream checkpoint all still agree.
    pub fn restore_state(
        &self,
        normal_patch: &NormalPatchView,
        interface: &PartialSlipInterface,
        state: TangentialContactState,
    ) -> Result<TangentialContactState, TangentialContactError> {
        self.validate_identity()?;
        match (&self.lane, &state) {
            (
                TangentialContactLane::PartialSlip { law },
                TangentialContactState::PartialSlip {
                    adapter_id,
                    source_id,
                    checkpoint,
                    committed_work_keys,
                    max_committed_work_keys,
                    ..
                },
            ) => {
                self.validate_state_identity(adapter_id, source_id)?;
                if *max_committed_work_keys == 0
                    || committed_work_keys.len() > *max_committed_work_keys
                {
                    return Err(TangentialContactError::InvalidCheckpoint {
                        field: "partial-slip work-key budget",
                    });
                }
                law.restore_checkpoint(normal_patch, interface, checkpoint)?;
                Ok(state)
            }
            (
                TangentialContactLane::Smooth { law, adapter },
                TangentialContactState::Smooth {
                    adapter_id,
                    source_id,
                    state: smooth_state,
                },
            ) => {
                self.validate_state_identity(adapter_id, source_id)?;
                adapter.restore_checkpoint(
                    law,
                    normal_patch,
                    interface,
                    smooth_state.checkpoint(),
                )?;
                Ok(state)
            }
            _ => Err(TangentialContactError::CheckpointLaneMismatch),
        }
    }

    /// Prepares exactly one tangential contact candidate without mutating
    /// `state`.  Call [`Self::commit`] only after the consuming mechanics step
    /// accepts the candidate; [`Self::rollback`] returns the exact prior state.
    pub fn prepare(
        &self,
        state: &TangentialContactState,
        request: &TangentialContactRequest,
    ) -> Result<TangentialContactReceipt, TangentialContactError> {
        self.validate_request(request)?;
        let prior_state =
            self.restore_state(&request.normal_patch, &request.interface, state.clone())?;
        if request.expected_state_version != prior_state.committed_version() {
            return Err(TangentialContactError::StateVersionMismatch {
                expected: prior_state.committed_version(),
                observed: request.expected_state_version,
            });
        }
        let frame = TangentFrame::new(
            vec3(request.patch_kinematics.tangent_basis.normal_world),
            vec3(request.patch_kinematics.tangent_basis.first_world),
        )?;
        let Creepage::Available {
            longitudinal,
            lateral,
            reference_rolling_speed_m_per_s,
        } = request.patch_kinematics.creepage
        else {
            return Err(TangentialContactError::CreepageUnavailable);
        };
        let kinematics = PartialSlipKinematics {
            creepage: [longitudinal, lateral],
            rolling_speed_mps: reference_rolling_speed_m_per_s,
            torsional_spin_rad_per_s: request.patch_kinematics.normal_spin_rad_per_s,
            dt_s: request.dt_s,
        };
        let (step, candidate) = match (&self.lane, &prior_state) {
            (
                TangentialContactLane::PartialSlip { law },
                TangentialContactState::PartialSlip {
                    checkpoint,
                    committed_version,
                    committed_work_keys,
                    max_committed_work_keys,
                    ..
                },
            ) => {
                if committed_work_keys
                    .iter()
                    .any(|key| key == &request.work_ownership)
                {
                    return Err(TangentialContactError::DuplicateWorkOwnership);
                }
                if committed_work_keys.len() >= *max_committed_work_keys {
                    return Err(TangentialContactError::WorkOwnershipCapacityExceeded {
                        max: *max_committed_work_keys,
                    });
                }
                let law_state =
                    law.restore_checkpoint(&request.normal_patch, &request.interface, checkpoint)?;
                let step = law.advance(
                    &request.normal_patch,
                    &request.interface,
                    frame,
                    kinematics,
                    &request.work_ownership,
                    &law_state,
                )?;
                let next_version = committed_version.checked_add(1).ok_or(
                    TangentialContactError::InvalidDerived {
                        field: "partial-slip committed version",
                    },
                )?;
                let mut next_keys = committed_work_keys.clone();
                next_keys.push(request.work_ownership.clone());
                let next_state = TangentialContactState::PartialSlip {
                    adapter_id: self.adapter_id.clone(),
                    source_id: self.source_id.clone(),
                    checkpoint: step.checkpoint.clone(),
                    committed_version: next_version,
                    committed_work_keys: next_keys,
                    max_committed_work_keys: *max_committed_work_keys,
                };
                (step, TangentialCandidate::PartialSlip { next_state })
            }
            (
                TangentialContactLane::Smooth { law, adapter },
                TangentialContactState::Smooth { state, .. },
            ) => {
                let smooth = adapter.prepare(
                    law,
                    &request.normal_patch,
                    &request.interface,
                    state,
                    &SmoothTangentialQuery {
                        query_id: request.request_id.clone(),
                        expected_state_version: request.expected_state_version,
                        frame,
                        kinematics,
                        work_ownership: request.work_ownership.clone(),
                    },
                )?;
                (
                    smooth.partial_slip_step.clone(),
                    TangentialCandidate::Smooth { receipt: smooth },
                )
            }
            _ => return Err(TangentialContactError::CheckpointLaneMismatch),
        };
        let (force_on_disc_world_n, free_torsional_torque_on_disc_world_nm) =
            wrench_on_disc(step.wrench, request.patch_kinematics.surfaces.order());
        let application_arm_world_m = request.patch_kinematics.disc_point.arm_world;
        let disc_endpoint_power_w = force_on_disc_world_n
            .dot(request.patch_kinematics.disc_point.point_velocity_world)
            + free_torsional_torque_on_disc_world_nm
                .dot(request.patch_kinematics.disc_point.angular_velocity_world);
        finite(disc_endpoint_power_w, "disc_endpoint_power_w")?;
        let disc_endpoint_work_j = disc_endpoint_power_w * request.dt_s;
        finite(disc_endpoint_work_j, "disc_endpoint_work_j")?;
        nonnegative(
            step.dissipation.tangential_and_torsional_microslip_j,
            "irreversible_loss_j",
        )?;
        nonnegative(step.dissipation.heat_j, "heat_j")?;
        Ok(TangentialContactReceipt {
            request_id: request.request_id.clone(),
            parent_version: prior_state.committed_version(),
            mode: step.state,
            force_on_disc_world_n,
            free_torsional_torque_on_disc_world_nm,
            application_arm_world_m,
            microslip_fraction: step.slip_partition.microslip_fraction,
            signed_storage_change_j: step.storage_change_j,
            stored_energy_j: step.stored_energy_j,
            irreversible_loss_j: step.dissipation.tangential_and_torsional_microslip_j,
            heat_j: step.dissipation.heat_j,
            exact_relative_body_power_w: step.generalized_work.reconstructed_body_power_w,
            exact_relative_body_work_j: -step.generalized_work.work_into_interface_j,
            endpoint_relative_power_w: step.generalized_work.endpoint_body_power_w,
            disc_endpoint_power_w,
            disc_endpoint_work_j,
            checkpoint: candidate.next_state(self, &prior_state)?,
            prior_state,
            candidate,
        })
    }

    /// Commits a prepared candidate exactly once.
    pub fn commit(
        &self,
        state: &TangentialContactState,
        receipt: &TangentialContactReceipt,
    ) -> Result<TangentialContactState, TangentialContactError> {
        if *state != receipt.prior_state {
            return Err(TangentialContactError::ReceiptDoesNotMatchState);
        }
        match (&self.lane, &receipt.candidate) {
            (
                TangentialContactLane::PartialSlip { .. },
                TangentialCandidate::PartialSlip { next_state },
            ) => Ok(next_state.clone()),
            (
                TangentialContactLane::Smooth { adapter, .. },
                TangentialCandidate::Smooth { receipt },
            ) => {
                let TangentialContactState::Smooth {
                    state: smooth_state,
                    ..
                } = state
                else {
                    return Err(TangentialContactError::CheckpointLaneMismatch);
                };
                Ok(TangentialContactState::Smooth {
                    adapter_id: self.adapter_id.clone(),
                    source_id: self.source_id.clone(),
                    state: adapter.commit(smooth_state, receipt)?,
                })
            }
            _ => Err(TangentialContactError::CheckpointLaneMismatch),
        }
    }

    /// Rolls back an uncommitted candidate to its exact parent state.
    pub fn rollback(
        &self,
        receipt: &TangentialContactReceipt,
    ) -> Result<TangentialContactState, TangentialContactError> {
        match (&self.lane, &receipt.candidate) {
            (
                TangentialContactLane::PartialSlip { .. },
                TangentialCandidate::PartialSlip { .. },
            ) => Ok(receipt.prior_state.clone()),
            (
                TangentialContactLane::Smooth { adapter, .. },
                TangentialCandidate::Smooth { receipt },
            ) => {
                let restored = adapter.rollback(receipt)?;
                Ok(TangentialContactState::Smooth {
                    adapter_id: self.adapter_id.clone(),
                    source_id: self.source_id.clone(),
                    state: restored,
                })
            }
            _ => Err(TangentialContactError::CheckpointLaneMismatch),
        }
    }

    fn validate_request(
        &self,
        request: &TangentialContactRequest,
    ) -> Result<(), TangentialContactError> {
        self.validate_identity()?;
        nonblank(&request.request_id, "request_id")?;
        if !request.dt_s.is_finite() || request.dt_s <= 0.0 {
            return Err(TangentialContactError::InvalidInput { field: "dt_s" });
        }
        if !request.patch_kinematics.disc_point.arm_world.is_finite() {
            return Err(TangentialContactError::InvalidInput {
                field: "patch_kinematics.disc_point.arm_world",
            });
        }
        if request.normal_patch.patch_id() != request.patch_kinematics.patch.patch_identity.as_str()
        {
            return Err(TangentialContactError::NormalPatchIdentityMismatch);
        }
        match request.patch_kinematics.status {
            PatchContactStatus::Touching | PatchContactStatus::Grazing => Ok(()),
            status => Err(TangentialContactError::UnsupportedPatchStatus { status }),
        }
    }

    fn validate_identity(&self) -> Result<(), TangentialContactError> {
        nonblank(&self.adapter_id, "adapter_id")?;
        nonblank(&self.source_id, "source_id")
    }

    fn validate_state_identity(
        &self,
        adapter_id: &str,
        source_id: &str,
    ) -> Result<(), TangentialContactError> {
        if adapter_id != self.adapter_id {
            return Err(TangentialContactError::CheckpointIdentityMismatch {
                field: "adapter_id",
            });
        }
        if source_id != self.source_id {
            return Err(TangentialContactError::CheckpointIdentityMismatch { field: "source_id" });
        }
        Ok(())
    }
}

/// One explicit Euler translation request.  No coefficient is accepted here:
/// the named generic law owns friction data and refuses invalid coefficients.
#[derive(Clone, Debug, PartialEq)]
pub struct TangentialContactRequest {
    /// Stable attempt identity.
    pub request_id: String,
    /// State version observed by the mechanics transaction.
    pub expected_state_version: u64,
    /// Pre-constitutive finite-patch geometry and kinematics.
    pub patch_kinematics: PatchKinematics,
    /// Finite-patch normal response from its owning normal-contact rung.
    pub normal_patch: NormalPatchView,
    /// Ordered dry interface/history identity.
    pub interface: PartialSlipInterface,
    /// Exact generic work ownership for this patch/time interval.
    pub work_ownership: GeneralizedWorkOwnership,
    /// Positive constitutive interval duration in seconds.
    pub dt_s: f64,
}

/// Replay-bound state.  It is the checkpoint surface for this adapter.
#[derive(Clone, Debug, PartialEq)]
pub enum TangentialContactState {
    /// Direct partial-slip lane with an explicit exactly-once work-key ledger.
    PartialSlip {
        adapter_id: String,
        source_id: String,
        checkpoint: PartialSlipCheckpoint,
        committed_version: u64,
        committed_work_keys: Vec<GeneralizedWorkOwnership>,
        max_committed_work_keys: usize,
    },
    /// Smooth lane, whose generic adapter owns the full transaction checkpoint.
    Smooth {
        adapter_id: String,
        source_id: String,
        state: SmoothTangentialState,
    },
}

impl TangentialContactState {
    /// Version of the last committed contact interval.
    #[must_use]
    pub const fn committed_version(&self) -> u64 {
        match self {
            Self::PartialSlip {
                committed_version, ..
            } => *committed_version,
            Self::Smooth { state, .. } => state.committed_version(),
        }
    }
}

/// Candidate receipt retaining all physical accounting without a defaulted
/// material coefficient, omitted torsion, or inferred work owner.
#[derive(Clone, Debug, PartialEq)]
pub struct TangentialContactReceipt {
    /// Stable caller attempt identity.
    pub request_id: String,
    /// Committed state version used for this candidate.
    pub parent_version: u64,
    /// Generic returned contact mode.
    pub mode: PartialSlipStateKind,
    /// Action on the Euler disc in world newtons.
    pub force_on_disc_world_n: Vec3,
    /// Free torsional torque on the Euler disc in world newton metres.
    pub free_torsional_torque_on_disc_world_nm: Vec3,
    /// World arm from disc centre of mass to the application point.
    pub application_arm_world_m: Vec3,
    /// Generic lumped microslip partition, not a resolved slip-area fraction.
    pub microslip_fraction: f64,
    /// Signed reversible storage change, J.
    pub signed_storage_change_j: f64,
    /// Reversible storage after the accepted constitutive interval, J.
    pub stored_energy_j: f64,
    /// Non-negative irreversible microslip loss, J.
    pub irreversible_loss_j: f64,
    /// Non-negative heat assigned by the generic rung, J.
    pub heat_j: f64,
    /// Exact discrete relative-body power reconstructed by the generic law, W.
    pub exact_relative_body_power_w: f64,
    /// Exact discrete relative-body work reconstructed by the generic law, J.
    pub exact_relative_body_work_j: f64,
    /// Endpoint relative-wrench power kept separately from discrete closure, W.
    pub endpoint_relative_power_w: f64,
    /// Disc endpoint-power diagnostic from force, arm, and free torque, W.
    pub disc_endpoint_power_w: f64,
    /// `disc_endpoint_power_w * dt_s`, retained as an endpoint diagnostic, J.
    pub disc_endpoint_work_j: f64,
    /// Exact next checkpoint identity/state for deterministic recontact replay.
    pub checkpoint: TangentialContactState,
    prior_state: TangentialContactState,
    candidate: TangentialCandidate,
}

#[derive(Clone, Debug, PartialEq)]
enum TangentialCandidate {
    PartialSlip {
        next_state: TangentialContactState,
    },
    Smooth {
        receipt: fs_contact::tangential::smooth::SmoothTangentialReceipt,
    },
}

impl TangentialCandidate {
    fn next_state(
        &self,
        adapter: &EulerTangentialContactAdapter,
        prior: &TangentialContactState,
    ) -> Result<TangentialContactState, TangentialContactError> {
        match self {
            Self::PartialSlip { next_state } => Ok(next_state.clone()),
            Self::Smooth { receipt } => {
                let TangentialContactLane::Smooth {
                    adapter: smooth, ..
                } = &adapter.lane
                else {
                    return Err(TangentialContactError::CheckpointLaneMismatch);
                };
                let TangentialContactState::Smooth { state, .. } = prior else {
                    return Err(TangentialContactError::CheckpointLaneMismatch);
                };
                Ok(TangentialContactState::Smooth {
                    adapter_id: adapter.adapter_id.clone(),
                    source_id: adapter.source_id.clone(),
                    state: smooth.commit(state, receipt)?,
                })
            }
        }
    }
}

/// Typed refusal surface.  No branch supplies target coefficients, ranks a
/// configuration, or upgrades the generic law's authority.
#[derive(Clone, Debug, PartialEq)]
pub enum TangentialContactError {
    /// A required identity was blank.
    MissingIdentity { field: &'static str },
    /// A scalar/vector interval input was invalid.
    InvalidInput { field: &'static str },
    /// A finite input generated an unrepresentable result.
    InvalidDerived { field: &'static str },
    /// Creepage cannot be normalized at the supplied entrainment speed.
    CreepageUnavailable,
    /// A pre-constitutive status cannot be treated as a tangential-law mode.
    UnsupportedPatchStatus { status: PatchContactStatus },
    /// The supplied normal receipt belongs to a different finite patch.
    NormalPatchIdentityMismatch,
    /// Caller state is stale or future for this adapter.
    StateVersionMismatch { expected: u64, observed: u64 },
    /// A checkpoint was created by a different wrapper configuration.
    CheckpointIdentityMismatch { field: &'static str },
    /// A direct and smooth checkpoint cannot be exchanged.
    CheckpointLaneMismatch,
    /// The direct lane's exactly-once work key was reused.
    DuplicateWorkOwnership,
    /// The direct lane's caller-declared work-key budget is exhausted.
    WorkOwnershipCapacityExceeded { max: usize },
    /// Checkpoint internals violate the explicit bounded state contract.
    InvalidCheckpoint { field: &'static str },
    /// A receipt was committed against a different same-version parent.
    ReceiptDoesNotMatchState,
    /// The delegated generic partial-slip rung refused the request.
    PartialSlip(PartialSlipError),
    /// The delegated generic smooth adapter refused the transaction.
    Smooth(SmoothTangentialError),
}

impl fmt::Display for TangentialContactError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self:?}")
    }
}

impl std::error::Error for TangentialContactError {}

impl From<PartialSlipError> for TangentialContactError {
    fn from(value: PartialSlipError) -> Self {
        Self::PartialSlip(value)
    }
}

impl From<SmoothTangentialError> for TangentialContactError {
    fn from(value: SmoothTangentialError) -> Self {
        Self::Smooth(value)
    }
}

fn wrench_on_disc(wrench: TangentialWrench, order: SurfaceOrder) -> (Vec3, Vec3) {
    let force = Vec3::new(wrench.force_n[0], wrench.force_n[1], wrench.force_n[2]);
    let torque = Vec3::new(
        wrench.torque_nm[0],
        wrench.torque_nm[1],
        wrench.torque_nm[2],
    );
    match order {
        SurfaceOrder::DiscThenBase => (force, torque),
        SurfaceOrder::BaseThenDisc => (force.scale(-1.0), torque.scale(-1.0)),
    }
}

fn vec3(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}

fn nonblank(value: &str, field: &'static str) -> Result<(), TangentialContactError> {
    if value.trim().is_empty() {
        Err(TangentialContactError::MissingIdentity { field })
    } else {
        Ok(())
    }
}

fn finite(value: f64, field: &'static str) -> Result<(), TangentialContactError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(TangentialContactError::InvalidDerived { field })
    }
}

fn nonnegative(value: f64, field: &'static str) -> Result<(), TangentialContactError> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(TangentialContactError::InvalidDerived { field })
    }
}
