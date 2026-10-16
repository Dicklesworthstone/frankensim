//! Thin world-frame adapter for generic rolling and contour-loss candidates.
//!
//! This module maps one already-declared `fs_tribo::rolling_loss` candidate
//! into a body wrench. It neither changes a constitutive coefficient nor
//! selects, blends, ranks, or calibrates the alternatives. Normal response,
//! gas, base, sliding, and microslip mechanisms are outside this adapter.

use core::fmt;

use fs_mbd::Vec3;
use fs_tribo::{
    InterfaceSystemRef,
    partial_slip::GeneralizedWorkOwnership,
    rolling_loss::{
        RollingKinematics, RollingLossCheckpoint, RollingLossError, RollingLossLaw,
        RollingLossState, RollingLossStep, RollingPatchReceipt, RollingWorkOwnership,
    },
};

/// Stable identity of this coordinate-only adapter.
pub const ROLLING_CONTACT_ADAPTER_ID: &str = "euler-disc/rolling-contact-adapter-v1";

const AXIS_TOLERANCE: f64 = 256.0 * f64::EPSILON;

/// Typed refusal from the rolling-contact adapter.
#[derive(Debug, Clone, PartialEq)]
pub enum RollingContactError {
    /// A frame, port, or domain identity was blank.
    MissingIdentity { field: &'static str },
    /// A supplied scalar/vector is not representable in this adapter's domain.
    InvalidInput { field: &'static str },
    /// The declared world axes are not an orthonormal pair.
    NonOrthonormalAxes,
    /// A supplied checkpoint does not restore to the supplied candidate state.
    CheckpointStateMismatch,
    /// The caller passed raw generic state/checkpoint different from the
    /// transaction state it asked this adapter to advance.
    StateInputMismatch,
    /// A proposal was accepted against a stale or different parent state.
    ProposalDoesNotMatchState,
    /// The generic rolling-loss leaf refused the retained caller inputs.
    GenericRefusal(RollingLossError),
}

impl fmt::Display for RollingContactError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for RollingContactError {}

impl From<RollingLossError> for RollingContactError {
    fn from(value: RollingLossError) -> Self {
        Self::GenericRefusal(value)
    }
}

/// Caller identities for one world-frame adapter invocation.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RollingContactIdentity {
    /// Caller case or run identity.
    pub case_id: String,
    /// Must equal [`ROLLING_CONTACT_ADAPTER_ID`].
    pub adapter_id: String,
    /// World-frame identity shared by the arm and both axes.
    pub world_frame_id: String,
    /// Caller composition-port identity.
    pub port_id: String,
    /// Caller contribution-domain identity.
    pub domain_id: String,
}

impl RollingContactIdentity {
    fn validate(&self) -> Result<(), RollingContactError> {
        for (value, field) in [
            (self.case_id.as_str(), "case_id"),
            (self.world_frame_id.as_str(), "world_frame_id"),
            (self.port_id.as_str(), "port_id"),
            (self.domain_id.as_str(), "domain_id"),
        ] {
            if value.trim().is_empty() {
                return Err(RollingContactError::MissingIdentity { field });
            }
        }
        if self.adapter_id != ROLLING_CONTACT_ADAPTER_ID {
            return Err(RollingContactError::InvalidInput {
                field: "adapter_id",
            });
        }
        Ok(())
    }
}

/// Explicit input for exactly one generic rolling-loss candidate.
///
/// The scalar speeds are signed against their corresponding unit world axes.
/// The adapter forwards them unchanged to `fs_tribo`; it does not infer either
/// speed from a body pose, contact radius, or another loss channel.
#[derive(Debug, Clone, PartialEq)]
pub struct RollingContactInput {
    /// Frame, port, and contribution-domain identities.
    pub identity: RollingContactIdentity,
    /// Caller-retained normal-patch receipt.
    pub patch: RollingPatchReceipt,
    /// Ordered dry interface and history receipt.
    pub interface: InterfaceSystemRef,
    /// One distinct generic candidate; alternatives are evaluated separately.
    pub law: RollingLossLaw,
    /// Candidate state to advance only after outer-solver acceptance.
    pub state: RollingLossState,
    /// Optional identity-bound predecessor receipt.
    pub checkpoint: Option<RollingLossCheckpoint>,
    /// Dedicated rolling-loss ownership key.
    pub ownership: RollingWorkOwnership,
    /// Optional partial-slip owner to reject conservatively before composition.
    pub partial_slip_ownership: Option<GeneralizedWorkOwnership>,
    /// Contact-point arm from body center of mass in the declared world frame [m].
    pub contact_arm_world_m: Vec3,
    /// Unit tangent axis signed with material contour speed.
    pub contour_tangent_axis_world: Vec3,
    /// Unit rolling axis signed with material rolling rate.
    pub rolling_axis_world: Vec3,
    /// Signed material contact-point contour speed [m/s].
    pub contour_speed_mps: f64,
    /// Signed rolling rate [rad/s].
    pub rolling_rate_rad_s: f64,
    /// Relative spin rate [rad/s], retained but unavailable in these rungs.
    pub spin_rate_rad_s: f64,
    /// Absolute temperature [K].
    pub temperature_kelvin: f64,
    /// Caller-declared excitation frequency [Hz].
    pub excitation_frequency_hz: f64,
    /// Accepted candidate interval duration [s].
    pub interval_s: f64,
}

/// Gas/base-independent body-wrench components in the declared world frame.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RollingContactBodyWrench {
    /// Contour force applied at the retained contact point [N].
    pub contour_force_world_n: Vec3,
    /// Moment about COM induced by `contact_arm_world_m × contour_force_world_n` [N m].
    pub contour_force_moment_about_com_world_nm: Vec3,
    /// Free rolling couple, independent of contact-arm recentering [N m].
    pub rolling_free_couple_world_nm: Vec3,
    /// Sum of the contact-force moment and free rolling couple [N m].
    pub total_moment_about_com_world_nm: Vec3,
    /// Wrench power reconstructed from the retained scalar speeds [W].
    pub body_power_w: f64,
}

/// Explicit availability of mechanisms deliberately absent from these rungs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpinMicroslipAvailability {
    /// Spin/microslip is not supplied by the selected rolling/contour candidate.
    Unavailable,
}

/// One mapped candidate result, retaining the exact generic response.
#[derive(Debug, Clone, PartialEq)]
pub struct RollingContactStep {
    /// Exact caller identities used for this mapping.
    pub identity: RollingContactIdentity,
    /// Immutable generic candidate model identity, retained without selection.
    pub law_model_id: String,
    /// Immutable generic candidate source-card identity, retained without promotion.
    pub law_source_id: String,
    /// Unmodified generic candidate response, including work, heat, storage,
    /// applicability, uncertainty, next state, and checkpoint.
    pub generic: RollingLossStep,
    /// World-frame rolling/contour components only; no gas/base contribution.
    pub body_wrench: RollingContactBodyWrench,
    /// Explicitly unavailable spin/microslip contribution.
    pub spin_microslip: SpinMicroslipAvailability,
}

/// Caller-owned accepted-step state for one rolling-loss channel.
///
/// The generic rolling law already produces a candidate `next_state` and
/// identity-bound checkpoint. This wrapper makes the outer accept/refuse
/// decision explicit without adding a second constitutive state machine. Work
/// ownership remains the caller's cross-channel responsibility: the adapter
/// rejects overlap with a supplied partial-slip owner, while the composition
/// ledger must record accepted rolling work exactly once for its interval.
#[derive(Debug, Clone, PartialEq)]
pub struct RollingContactState {
    /// Generic cumulative irreversible-loss state.
    pub generic_state: RollingLossState,
    /// Optional checkpoint for identity-bound deterministic replay.
    pub checkpoint: Option<RollingLossCheckpoint>,
    /// Number of outer mechanics steps that accepted a candidate.
    pub committed_version: u64,
}

impl RollingContactState {
    /// Starts a fresh rolling channel with no accepted interval.
    #[must_use]
    pub const fn zero() -> Self {
        Self {
            generic_state: RollingLossState::zero(),
            checkpoint: None,
            committed_version: 0,
        }
    }
}

impl Default for RollingContactState {
    fn default() -> Self {
        Self::zero()
    }
}

/// An uncommitted rolling candidate paired with the exact parent and successor
/// state snapshots. It contains no mutable global state.
#[derive(Debug, Clone, PartialEq)]
pub struct RollingContactProposal {
    /// Mapped physical candidate and generic work/heat receipt.
    pub step: RollingContactStep,
    /// Parent version observed by the proposal.
    pub parent_version: u64,
    prior_state: RollingContactState,
    next_state: RollingContactState,
}

/// Evaluates one generic rolling-loss candidate and maps it into a world wrench.
///
/// The caller decides whether to commit `generic.next_state`; no external state
/// is mutated by this adapter.
pub fn evaluate_rolling_contact(
    input: &RollingContactInput,
) -> Result<RollingContactStep, RollingContactError> {
    input.identity.validate()?;
    validate_axes(input.contour_tangent_axis_world, input.rolling_axis_world)?;
    finite_vec(input.contact_arm_world_m, "contact_arm_world_m")?;
    for (value, field) in [
        (input.contour_speed_mps, "contour_speed_mps"),
        (input.rolling_rate_rad_s, "rolling_rate_rad_s"),
        (input.spin_rate_rad_s, "spin_rate_rad_s"),
        (input.temperature_kelvin, "temperature_kelvin"),
        (input.excitation_frequency_hz, "excitation_frequency_hz"),
        (input.interval_s, "interval_s"),
    ] {
        finite(value, field)?;
    }
    if let Some(partial_slip) = &input.partial_slip_ownership {
        input
            .ownership
            .require_disjoint_from_partial_slip(partial_slip)?;
    }
    if let Some(checkpoint) = &input.checkpoint {
        let restored = input.law.restore_checkpoint(
            &input.patch,
            &input.interface,
            &input.ownership,
            checkpoint,
        )?;
        if restored != input.state {
            return Err(RollingContactError::CheckpointStateMismatch);
        }
    }
    let kinematics = RollingKinematics::new(
        input.contour_speed_mps,
        input.rolling_rate_rad_s,
        input.spin_rate_rad_s,
        input.temperature_kelvin,
        input.excitation_frequency_hz,
        input.interval_s,
    )?;
    let generic = input.law.advance(
        &input.patch,
        &input.interface,
        kinematics,
        &input.ownership,
        &input.state,
    )?;
    let contour_force_world_n = checked_scale(
        input.contour_tangent_axis_world,
        generic.wrench.contour_force_n,
        "contour_force_world_n",
    )?;
    let contour_force_moment_about_com_world_nm = checked_cross(
        input.contact_arm_world_m,
        contour_force_world_n,
        "contour_force_moment_about_com_world_nm",
    )?;
    let rolling_free_couple_world_nm = checked_scale(
        input.rolling_axis_world,
        generic.wrench.rolling_moment_nm,
        "rolling_free_couple_world_nm",
    )?;
    let total_moment_about_com_world_nm = checked_add(
        contour_force_moment_about_com_world_nm,
        rolling_free_couple_world_nm,
        "total_moment_about_com_world_nm",
    )?;
    let body_power_w = checked_scalar_add(
        checked_mul(
            generic.wrench.contour_force_n,
            input.contour_speed_mps,
            "contour_body_power_w",
        )?,
        checked_mul(
            generic.wrench.rolling_moment_nm,
            input.rolling_rate_rad_s,
            "rolling_body_power_w",
        )?,
        "body_power_w",
    )?;
    if !approximately_equal(body_power_w, generic.generalized_work.endpoint_body_power_w) {
        return Err(RollingContactError::InvalidInput {
            field: "generic_power_mapping",
        });
    }
    Ok(RollingContactStep {
        identity: input.identity.clone(),
        law_model_id: input.law.model_id().to_owned(),
        law_source_id: input.law.source_id().to_owned(),
        generic,
        body_wrench: RollingContactBodyWrench {
            contour_force_world_n,
            contour_force_moment_about_com_world_nm,
            rolling_free_couple_world_nm,
            total_moment_about_com_world_nm,
            body_power_w,
        },
        spin_microslip: SpinMicroslipAvailability::Unavailable,
    })
}

/// Stages one rolling candidate against an explicit accepted-step state.
///
/// `input.state` must be an exact copy of `state`. If the caller supplies a
/// checkpoint, it must be the retained checkpoint from that state; omitting it
/// is intentional for a time-evolving patch, interface, or ownership identity
/// because a generic checkpoint is bound to those exact inputs. The returned
/// proposal changes no state. Call [`commit_rolling_contact`] only after the
/// outer rigid-body step accepts the complete multi-channel candidate, or
/// [`rollback_rolling_contact`] to recover the parent snapshot.
pub fn prepare_rolling_contact(
    state: &RollingContactState,
    input: &RollingContactInput,
) -> Result<RollingContactProposal, RollingContactError> {
    if input.state != state.generic_state
        || input
            .checkpoint
            .as_ref()
            .is_some_and(|checkpoint| state.checkpoint.as_ref() != Some(checkpoint))
    {
        return Err(RollingContactError::StateInputMismatch);
    }
    let step = evaluate_rolling_contact(input)?;
    let committed_version =
        state
            .committed_version
            .checked_add(1)
            .ok_or(RollingContactError::InvalidInput {
                field: "committed_version",
            })?;
    Ok(RollingContactProposal {
        parent_version: state.committed_version,
        next_state: RollingContactState {
            generic_state: step.generic.next_state.clone(),
            checkpoint: Some(step.generic.checkpoint.clone()),
            committed_version,
        },
        prior_state: state.clone(),
        step,
    })
}

/// Accepts a staged rolling candidate exactly against its recorded parent.
pub fn commit_rolling_contact(
    state: &RollingContactState,
    proposal: &RollingContactProposal,
) -> Result<RollingContactState, RollingContactError> {
    if *state != proposal.prior_state || state.committed_version != proposal.parent_version {
        return Err(RollingContactError::ProposalDoesNotMatchState);
    }
    Ok(proposal.next_state.clone())
}

/// Refuses a staged rolling candidate and returns its exact parent snapshot.
#[must_use]
pub fn rollback_rolling_contact(proposal: &RollingContactProposal) -> RollingContactState {
    proposal.prior_state.clone()
}

fn validate_axes(contour: Vec3, rolling: Vec3) -> Result<(), RollingContactError> {
    let contour_norm = stable_norm(contour, "contour_tangent_axis_world")?;
    let rolling_norm = stable_norm(rolling, "rolling_axis_world")?;
    if (contour_norm - 1.0).abs() > AXIS_TOLERANCE
        || (rolling_norm - 1.0).abs() > AXIS_TOLERANCE
        || contour.dot(rolling).abs() > AXIS_TOLERANCE
    {
        return Err(RollingContactError::NonOrthonormalAxes);
    }
    Ok(())
}

fn finite(value: f64, field: &'static str) -> Result<(), RollingContactError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(RollingContactError::InvalidInput { field })
    }
}

fn finite_vec(value: Vec3, field: &'static str) -> Result<(), RollingContactError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(RollingContactError::InvalidInput { field })
    }
}

fn stable_norm(value: Vec3, field: &'static str) -> Result<f64, RollingContactError> {
    finite_vec(value, field)?;
    let scale = value.x.abs().max(value.y.abs()).max(value.z.abs());
    if scale == 0.0 {
        return Ok(0.0);
    }
    let scaled = Vec3::new(value.x / scale, value.y / scale, value.z / scale);
    let norm = scale * scaled.dot(scaled).sqrt();
    if norm.is_finite() {
        Ok(norm)
    } else {
        Err(RollingContactError::InvalidInput { field })
    }
}

fn checked_mul(left: f64, right: f64, field: &'static str) -> Result<f64, RollingContactError> {
    let value = left * right;
    finite(value, field)?;
    Ok(value)
}

fn checked_scale(
    vector: Vec3,
    scalar: f64,
    field: &'static str,
) -> Result<Vec3, RollingContactError> {
    let value = vector.scale(scalar);
    finite_vec(value, field)?;
    Ok(value)
}

fn checked_cross(
    left: Vec3,
    right: Vec3,
    field: &'static str,
) -> Result<Vec3, RollingContactError> {
    let value = left.cross(right);
    finite_vec(value, field)?;
    Ok(value)
}

fn checked_add(left: Vec3, right: Vec3, field: &'static str) -> Result<Vec3, RollingContactError> {
    let value = left.add(right);
    finite_vec(value, field)?;
    Ok(value)
}

fn checked_scalar_add(
    left: f64,
    right: f64,
    field: &'static str,
) -> Result<f64, RollingContactError> {
    let value = left + right;
    finite(value, field)?;
    Ok(value)
}

fn approximately_equal(left: f64, right: f64) -> bool {
    (left - right).abs() <= 512.0 * f64::EPSILON * left.abs().max(right.abs()).max(1.0)
}
