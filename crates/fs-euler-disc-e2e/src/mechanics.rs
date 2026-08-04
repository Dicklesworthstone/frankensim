//! Atomic, fixed-step composition of declared Euler-disc mechanical wrenches.
//!
//! This is a composition boundary, not a contact, impact, base, gas, film, or
//! constitutive model.  It admits caller-supplied, fixed-branch normal,
//! tangential, rolling, and impact contributions only.  The underlying
//! [`RigidBodyIntegrator`] is the smooth constant-wrench midpoint rung; an
//! eventful impact law must be resolved before it reaches this boundary.

use core::fmt;
use std::collections::BTreeSet;

use fs_couple::{PortKind, PortOrientation, StableId};
use fs_mbd::{
    DynamicsError, Gravity, MassProperties, RigidBodyIntegrator, RigidBodyState, StepReceipt, Vec3,
    Wrench,
};

use crate::ports::{
    ChannelActivity, ContributionOwnership, EulerChannel, EulerPortError, EulerPortRegistry,
    PortDeclaration,
};

/// Upper bound on mechanical contributions in one atomic macro-step.
pub const MAX_MECHANICS_CONTRIBUTIONS: usize = 32;

/// The reference point assumed by this composition boundary.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MechanicsReference {
    /// The application arm is relative to the integrated body's centre of mass.
    CenterOfMass,
}

/// The role of a supplied action/reaction vector relative to the integrated body.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BodyAction {
    /// The supplied world force and free torque act on the integrated body.
    ActionOnIntegratedBody,
    /// A reaction on a different body must not be summed into this body's wrench.
    ReactionOnOtherBody,
}

/// Explicit status of a channel this rung does not execute.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InactiveMechanicsChannel {
    /// The declared Euler channel.
    pub channel: EulerChannel,
    /// The channel must be inactive or unavailable, never active.
    pub activity: ChannelActivity,
}

/// Required no-model declarations for base and gas lanes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct InactiveMechanicsChannels {
    /// Flexible or rigid base channel declaration.
    pub base: InactiveMechanicsChannel,
    /// Exterior-gas channel declaration.
    pub external_gas: InactiveMechanicsChannel,
    /// Thin-film gas channel declaration.
    pub gas_film: InactiveMechanicsChannel,
}

/// A caller-declared constant wrench contribution for one macro-step.
#[derive(Clone, Debug, PartialEq)]
pub struct MechanicsContribution {
    /// Complete port identity, channel, law/source identities, and ownership domain.
    pub port: PortDeclaration,
    /// This contribution must be the action on the integrated body.
    pub body_action: BodyAction,
    /// Declared location of the arm reference.
    pub reference: MechanicsReference,
    /// Application-point arm from the centre of mass, expressed in body metres.
    pub application_arm_body_m: Vec3,
    /// Force at the application point, expressed in world newtons.
    pub force_world_n: Vec3,
    /// Torque about the application point, expressed in world newton metres.
    pub free_torque_world_nm: Vec3,
    /// Caller-claimed discrete work on the body over this macro-step, in joules.
    pub claimed_discrete_work_j: f64,
    /// Caller-claimed signed recoverable-storage change, in joules.
    ///
    /// A negative value is an unloading/release of previously recoverable
    /// storage into the integrated body.
    pub claimed_storage_j: f64,
    /// Caller-claimed non-negative irreversible work, in joules.
    pub claimed_dissipation_j: f64,
    /// Caller-claimed heat partition/diagnostic of irreversible work, in joules.
    ///
    /// This must lie in `[0, claimed_dissipation_j]`; it is diagnostic only and
    /// is never added as another loss in either energy identity.
    pub claimed_heat_j: f64,
    /// Non-negative uncertainty allowance attached to this contribution, in joules.
    pub uncertainty_j: f64,
}

/// Declared absolute-plus-relative energy acceptance bound.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EnergyAcceptanceBound {
    /// Absolute allowance in joules.
    pub absolute_j: f64,
    /// Dimensionless allowance multiplied by the compared energy scale.
    pub relative: f64,
}

impl EnergyAcceptanceBound {
    /// Tests whether an energy residual is admitted at the supplied energy scale.
    #[must_use]
    pub fn admits(self, residual_j: f64, scale_j: f64) -> bool {
        self.absolute_j.is_finite()
            && self.absolute_j >= 0.0
            && self.relative.is_finite()
            && self.relative >= 0.0
            && residual_j.is_finite()
            && scale_j.is_finite()
            && scale_j >= 0.0
            && residual_j.abs() <= self.absolute_j + self.relative * scale_j
    }
}

/// Complete input to one atomic smooth mechanics macro-step.
#[derive(Clone, Debug, PartialEq)]
pub struct MechanicsMacroStepInput {
    /// State before the attempted step.
    pub state: RigidBodyState,
    /// Mass and principal inertia of the integrated body.
    pub mass_properties: MassProperties,
    /// Uniform world-frame gravity.
    pub gravity: Gravity,
    /// Fixed macro-step duration in seconds.
    pub duration_seconds: f64,
    /// Exact port identities that must be present once each.
    pub expected_contribution_keys: Vec<StableId>,
    /// Constant-wrench contributions admitted for this step.
    pub contributions: Vec<MechanicsContribution>,
    /// Explicit no-model declaration for base/external-gas/gas-film channels.
    pub inactive_channels: InactiveMechanicsChannels,
    /// Identity of the common world frame used by all mechanical port bindings.
    pub world_frame: StableId,
    /// Acceptance bound for caller claims and the mechanical energy closure.
    pub energy_bound: EnergyAcceptanceBound,
    /// A boundary cancellation flag.  Cancellation is checked before integration.
    pub cancelled_before_step: bool,
}

/// Sum of admitted contribution wrenches, centered at the body's centre of mass.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RecenteredWrenchReceipt {
    /// Sum of application-point forces in world newtons.
    pub force_world_n: Vec3,
    /// Sum of application-arm plus free torques about the COM in world newton metres.
    pub torque_center_of_mass_world_nm: Vec3,
    /// Same COM torque rotated into fs-mbd's principal body frame.
    pub torque_center_of_mass_body_nm: Vec3,
}

/// Independently reconstructed discrete midpoint data for one contribution.
#[derive(Clone, Debug, PartialEq)]
pub struct ContributionMechanicsReceipt {
    /// Contribution/port identity.
    pub identity: StableId,
    /// Typed channel family.
    pub channel: EulerChannel,
    /// Reconstructed midpoint application point in world metres.
    pub midpoint_point_world_m: Vec3,
    /// Reconstructed midpoint application arm in world metres.
    pub midpoint_arm_world_m: Vec3,
    /// Endpoint-centred midpoint point velocity in world m/s.
    pub midpoint_point_velocity_world_m_per_s: Vec3,
    /// Endpoint-centred midpoint angular velocity in world rad/s.
    pub midpoint_angular_velocity_world_rad_per_s: Vec3,
    /// Recentered torque about the COM at the reconstructed midpoint, in N m.
    pub midpoint_torque_center_of_mass_world_nm: Vec3,
    /// Independently reconstructed work, in joules.
    pub reconstructed_discrete_work_j: f64,
    /// Claimed work minus reconstructed work, in joules.
    pub claimed_work_residual_j: f64,
}

/// Energy accounting published with an admitted macro-step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MechanicsEnergyReceipt {
    /// Mechanical energy before the step (kinetic plus gravity potential), in J.
    pub mechanical_energy_before_j: f64,
    /// Mechanical energy after the step (kinetic plus gravity potential), in J.
    pub mechanical_energy_after_j: f64,
    /// Change in mechanical energy, in J.
    pub mechanical_energy_change_j: f64,
    /// Sum of independently reconstructed body work, in J.
    pub reconstructed_body_work_j: f64,
    /// Sum of caller-claimed body work, in J.
    pub claimed_body_work_j: f64,
    /// Signed change in recoverable contact storage, in J.
    pub recoverable_storage_change_j: f64,
    /// Non-negative irreversible contact/impact work, in J.
    pub irreversible_work_j: f64,
    /// Heat partition/diagnostic of irreversible work, in J.
    pub heat_j: f64,
    /// Sum of caller-provided uncertainty allowances, in J.
    pub uncertainty_j: f64,
    /// `delta(K + Ug) - reconstructed_body_work`, in J.
    pub body_energy_residual_j: f64,
    /// `-claimed_body_work - (recoverable_storage_change + irreversible_work)`, in J.
    pub interface_energy_residual_j: f64,
    /// Sum of the body and interface residuals, in J.
    pub combined_energy_residual_j: f64,
    /// Declared acceptance bound used for this receipt.
    pub acceptance_bound: EnergyAcceptanceBound,
}

/// Atomically published result of one accepted macro-step.
#[derive(Clone, Debug, PartialEq)]
pub struct MechanicsMacroStepReceipt {
    /// The underlying fixed-step receipt.
    pub rigid_body_step: StepReceipt,
    /// Recentered total wrench passed to the rigid-body integrator.
    pub resultant: RecenteredWrenchReceipt,
    /// One exact-once reconstruction receipt per declared contribution.
    pub contributions: Vec<ContributionMechanicsReceipt>,
    /// Energy accounting and closure diagnostic.
    pub energy: MechanicsEnergyReceipt,
}

/// Typed refusal from this composition boundary; no state is published on error.
#[derive(Clone, Debug, PartialEq)]
pub enum MechanicsMacroStepError {
    /// The caller cancelled at the whole-step boundary before any integration.
    CancelledBeforeStep,
    /// A scalar, vector, or duration was non-finite or out of range.
    InvalidInput(&'static str),
    /// The explicit expected key list repeated an identity.
    DuplicateExpectedKey { identity: StableId },
    /// A supplied port identity was not expected.
    UnexpectedContribution { identity: StableId },
    /// An expected port identity was absent.
    MissingContribution { identity: StableId },
    /// A supplied contribution repeated a port identity.
    DuplicateContribution { identity: StableId },
    /// This macro-step only accepts normal/tangential/rolling/impact channels.
    UnsupportedActiveChannel { channel: EulerChannel },
    /// An inactive or unavailable declaration attempted to supply a wrench.
    InactiveContribution { identity: StableId },
    /// The contribution is labelled as a reaction on another body.
    ReactionCannotActOnIntegratedBody { identity: StableId },
    /// The supplied application arm was not COM-relative.
    UnsupportedReference { identity: StableId },
    /// The declared coordinate frame differs from the input world-frame identity.
    WrongWorldFrame { identity: StableId },
    /// The declared coordinate orientation is not the common-frame convention.
    WrongCoordinateOrientation { identity: StableId },
    /// The port kind does not match the supplied force/free-torque representation.
    WrongPortKind { identity: StableId },
    /// Additive/shared ownership is deliberately not admitted by this atomic rung.
    OverlappingOwnership { identity: StableId },
    /// The port registry found duplicate or overlapping ownership.
    PortOwnership(EulerPortError),
    /// fs-mbd rejected state, mass properties, gravity, or the resultant wrench.
    Dynamics(DynamicsError),
    /// A caller work claim disagreed with the independent reconstruction.
    ContributionWorkMismatch { identity: StableId, residual_j: f64 },
    /// The rigid-body energy change disagreed with reconstructed body work.
    BodyEnergyMismatch { residual_j: f64 },
    /// Claimed body work did not balance recoverable storage plus irreversible work.
    InterfaceEnergyMismatch { residual_j: f64 },
    /// The composed body-plus-interface accounting did not close.
    CombinedEnergyMismatch { residual_j: f64 },
}

impl fmt::Display for MechanicsMacroStepError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::CancelledBeforeStep => formatter.write_str("cancelled before mechanics step"),
            Self::InvalidInput(field) => write!(formatter, "invalid mechanics input: {field}"),
            Self::DuplicateExpectedKey { identity } => {
                write!(
                    formatter,
                    "duplicate expected contribution key: {}",
                    identity.as_str()
                )
            }
            Self::UnexpectedContribution { identity } => {
                write!(formatter, "unexpected contribution: {}", identity.as_str())
            }
            Self::MissingContribution { identity } => {
                write!(formatter, "missing contribution: {}", identity.as_str())
            }
            Self::DuplicateContribution { identity } => {
                write!(formatter, "duplicate contribution: {}", identity.as_str())
            }
            Self::UnsupportedActiveChannel { channel } => {
                write!(
                    formatter,
                    "unsupported active mechanics channel: {channel:?}"
                )
            }
            Self::InactiveContribution { identity } => {
                write!(
                    formatter,
                    "inactive contribution supplied: {}",
                    identity.as_str()
                )
            }
            Self::ReactionCannotActOnIntegratedBody { identity } => write!(
                formatter,
                "reaction cannot be summed as an action on the integrated body: {}",
                identity.as_str()
            ),
            Self::UnsupportedReference { identity } => {
                write!(
                    formatter,
                    "unsupported contribution reference: {}",
                    identity.as_str()
                )
            }
            Self::WrongWorldFrame { identity } => {
                write!(
                    formatter,
                    "wrong contribution world frame: {}",
                    identity.as_str()
                )
            }
            Self::WrongCoordinateOrientation { identity } => write!(
                formatter,
                "wrong contribution coordinate orientation: {}",
                identity.as_str()
            ),
            Self::WrongPortKind { identity } => {
                write!(
                    formatter,
                    "wrong contribution port kind: {}",
                    identity.as_str()
                )
            }
            Self::OverlappingOwnership { identity } => {
                write!(
                    formatter,
                    "shared/additive ownership is not admitted: {}",
                    identity.as_str()
                )
            }
            Self::PortOwnership(error) => write!(formatter, "port ownership refusal: {error}"),
            Self::Dynamics(error) => write!(formatter, "rigid-body refusal: {error}"),
            Self::ContributionWorkMismatch {
                identity,
                residual_j,
            } => write!(
                formatter,
                "contribution work claim mismatches reconstruction for {}: {residual_j} J",
                identity.as_str()
            ),
            Self::BodyEnergyMismatch { residual_j } => {
                write!(formatter, "mechanics body-energy mismatch: {residual_j} J")
            }
            Self::InterfaceEnergyMismatch { residual_j } => {
                write!(
                    formatter,
                    "mechanics interface-energy mismatch: {residual_j} J"
                )
            }
            Self::CombinedEnergyMismatch { residual_j } => {
                write!(
                    formatter,
                    "mechanics combined-energy mismatch: {residual_j} J"
                )
            }
        }
    }
}

impl std::error::Error for MechanicsMacroStepError {}

/// Advances one admitted fixed mechanics macro-step atomically.
///
/// The returned [`MechanicsMacroStepReceipt`] is the only publication path. A
/// refusal occurs before the output exists, so callers retain their input
/// state.  Base, exterior-gas, and gas-film lanes are status declarations only
/// and cannot introduce hidden work through this rung.
pub fn run_mechanics_macro_step(
    input: MechanicsMacroStepInput,
) -> Result<MechanicsMacroStepReceipt, MechanicsMacroStepError> {
    if input.cancelled_before_step {
        return Err(MechanicsMacroStepError::CancelledBeforeStep);
    }
    validate_input_shape(&input)?;
    validate_inactive_channels(&input.inactive_channels)?;
    validate_expected_keys(&input.expected_contribution_keys, &input.contributions)?;

    let mut ports = Vec::with_capacity(input.contributions.len());
    let mut total_force_world_n = Vec3::ZERO;
    let mut total_torque_center_of_mass_world_nm = Vec3::ZERO;
    for contribution in &input.contributions {
        validate_contribution(contribution, &input.world_frame)?;
        let arm_world_m = input
            .state
            .pose()
            .orientation()
            .rotate_body_to_world(contribution.application_arm_body_m);
        let torque_world_nm = arm_world_m
            .cross(contribution.force_world_n)
            .add(contribution.free_torque_world_nm);
        if !arm_world_m.is_finite() || !torque_world_nm.is_finite() {
            return Err(MechanicsMacroStepError::InvalidInput(
                "derived contribution wrench",
            ));
        }
        total_force_world_n = total_force_world_n.add(contribution.force_world_n);
        total_torque_center_of_mass_world_nm =
            total_torque_center_of_mass_world_nm.add(torque_world_nm);
        ports.push(contribution.port.clone());
    }
    if !total_force_world_n.is_finite() || !total_torque_center_of_mass_world_nm.is_finite() {
        return Err(MechanicsMacroStepError::InvalidInput("resultant wrench"));
    }
    let registry_identity = StableId::new("euler-mechanics-step-registry")
        .map_err(|_| MechanicsMacroStepError::InvalidInput("internal registry identity"))?;
    EulerPortRegistry::try_new(registry_identity, ports)
        .map_err(MechanicsMacroStepError::PortOwnership)?;

    let torque_center_of_mass_body_nm = input
        .state
        .pose()
        .orientation()
        .rotate_world_to_body(total_torque_center_of_mass_world_nm);
    if !torque_center_of_mass_body_nm.is_finite() {
        return Err(MechanicsMacroStepError::InvalidInput(
            "body-frame resultant torque",
        ));
    }
    let resultant = RecenteredWrenchReceipt {
        force_world_n: total_force_world_n,
        torque_center_of_mass_world_nm: total_torque_center_of_mass_world_nm,
        torque_center_of_mass_body_nm,
    };
    let integrator = RigidBodyIntegrator::new(input.gravity);
    let rigid_body_step = integrator
        .step(
            input.state,
            input.mass_properties,
            Wrench {
                force_world: resultant.force_world_n,
                torque_body: resultant.torque_center_of_mass_body_nm,
            },
            input.duration_seconds,
        )
        .map_err(MechanicsMacroStepError::Dynamics)?;

    let mut contribution_receipts = Vec::with_capacity(input.contributions.len());
    let mut reconstructed_body_work_j = 0.0;
    let mut claimed_body_work_j = 0.0;
    let mut recoverable_storage_change_j = 0.0;
    let mut irreversible_work_j = 0.0;
    let mut heat_j = 0.0;
    let mut uncertainty_j = 0.0;
    for contribution in &input.contributions {
        let receipt = reconstruct_contribution(
            contribution,
            &rigid_body_step,
            input.mass_properties,
            input.duration_seconds,
        )?;
        let scale_j = energy_scale(
            contribution.claimed_discrete_work_j,
            receipt.reconstructed_discrete_work_j,
            contribution.uncertainty_j,
        );
        if !input
            .energy_bound
            .admits(receipt.claimed_work_residual_j, scale_j)
        {
            return Err(MechanicsMacroStepError::ContributionWorkMismatch {
                identity: receipt.identity,
                residual_j: receipt.claimed_work_residual_j,
            });
        }
        reconstructed_body_work_j += receipt.reconstructed_discrete_work_j;
        claimed_body_work_j += contribution.claimed_discrete_work_j;
        recoverable_storage_change_j += contribution.claimed_storage_j;
        irreversible_work_j += contribution.claimed_dissipation_j;
        heat_j += contribution.claimed_heat_j;
        uncertainty_j += contribution.uncertainty_j;
        contribution_receipts.push(receipt);
    }
    let mechanical_energy_before_j = rigid_body_step.diagnostics_before.mechanical_energy;
    let mechanical_energy_after_j = rigid_body_step.diagnostics_after.mechanical_energy;
    let mechanical_energy_change_j = mechanical_energy_after_j - mechanical_energy_before_j;
    let body_energy_residual_j = mechanical_energy_change_j - reconstructed_body_work_j;
    let interface_energy_residual_j =
        -claimed_body_work_j - (recoverable_storage_change_j + irreversible_work_j);
    let combined_energy_residual_j =
        combined_energy_residual_j(body_energy_residual_j, interface_energy_residual_j);
    if !all_finite(&[
        mechanical_energy_before_j,
        mechanical_energy_after_j,
        mechanical_energy_change_j,
        reconstructed_body_work_j,
        claimed_body_work_j,
        recoverable_storage_change_j,
        irreversible_work_j,
        heat_j,
        uncertainty_j,
        body_energy_residual_j,
        interface_energy_residual_j,
        combined_energy_residual_j,
    ]) {
        return Err(MechanicsMacroStepError::InvalidInput(
            "derived energy accounting",
        ));
    }
    let body_scale_j = energy_scale(
        mechanical_energy_change_j,
        reconstructed_body_work_j,
        uncertainty_j,
    );
    if !input
        .energy_bound
        .admits(body_energy_residual_j, body_scale_j)
    {
        return Err(MechanicsMacroStepError::BodyEnergyMismatch {
            residual_j: body_energy_residual_j,
        });
    }
    let interface_scale_j = energy_scale(
        -claimed_body_work_j,
        recoverable_storage_change_j + irreversible_work_j,
        uncertainty_j,
    );
    if !input
        .energy_bound
        .admits(interface_energy_residual_j, interface_scale_j)
    {
        return Err(MechanicsMacroStepError::InterfaceEnergyMismatch {
            residual_j: interface_energy_residual_j,
        });
    }
    let combined_scale_j = energy_scale(
        mechanical_energy_change_j,
        -(recoverable_storage_change_j + irreversible_work_j),
        uncertainty_j,
    );
    if !input
        .energy_bound
        .admits(combined_energy_residual_j, combined_scale_j)
    {
        return Err(MechanicsMacroStepError::CombinedEnergyMismatch {
            residual_j: combined_energy_residual_j,
        });
    }
    Ok(MechanicsMacroStepReceipt {
        rigid_body_step,
        resultant,
        contributions: contribution_receipts,
        energy: MechanicsEnergyReceipt {
            mechanical_energy_before_j,
            mechanical_energy_after_j,
            mechanical_energy_change_j,
            reconstructed_body_work_j,
            claimed_body_work_j,
            recoverable_storage_change_j,
            irreversible_work_j,
            heat_j,
            uncertainty_j,
            body_energy_residual_j,
            interface_energy_residual_j,
            combined_energy_residual_j,
            acceptance_bound: input.energy_bound,
        },
    })
}

fn validate_input_shape(input: &MechanicsMacroStepInput) -> Result<(), MechanicsMacroStepError> {
    if input.contributions.len() > MAX_MECHANICS_CONTRIBUTIONS {
        return Err(MechanicsMacroStepError::InvalidInput("contribution count"));
    }
    if !input.duration_seconds.is_finite() || input.duration_seconds <= 0.0 {
        return Err(MechanicsMacroStepError::InvalidInput("duration_seconds"));
    }
    if !input.energy_bound.absolute_j.is_finite()
        || input.energy_bound.absolute_j < 0.0
        || !input.energy_bound.relative.is_finite()
        || input.energy_bound.relative < 0.0
    {
        return Err(MechanicsMacroStepError::InvalidInput("energy_bound"));
    }
    Ok(())
}

fn validate_inactive_channels(
    channels: &InactiveMechanicsChannels,
) -> Result<(), MechanicsMacroStepError> {
    for (expected, declaration) in [
        (EulerChannel::Base, &channels.base),
        (EulerChannel::ExternalGas, &channels.external_gas),
        (EulerChannel::GasFilm, &channels.gas_film),
    ] {
        if declaration.channel != expected || declaration.activity.is_active() {
            return Err(MechanicsMacroStepError::UnsupportedActiveChannel {
                channel: declaration.channel,
            });
        }
    }
    Ok(())
}

fn validate_expected_keys(
    expected: &[StableId],
    contributions: &[MechanicsContribution],
) -> Result<(), MechanicsMacroStepError> {
    let mut expected_set = BTreeSet::new();
    for key in expected {
        if !expected_set.insert(key.clone()) {
            return Err(MechanicsMacroStepError::DuplicateExpectedKey {
                identity: key.clone(),
            });
        }
    }
    let mut actual_set = BTreeSet::new();
    for contribution in contributions {
        let identity = contribution.port.identity().clone();
        if !actual_set.insert(identity.clone()) {
            return Err(MechanicsMacroStepError::DuplicateContribution { identity });
        }
        if !expected_set.contains(&identity) {
            return Err(MechanicsMacroStepError::UnexpectedContribution { identity });
        }
    }
    for identity in expected_set {
        if !actual_set.contains(&identity) {
            return Err(MechanicsMacroStepError::MissingContribution { identity });
        }
    }
    Ok(())
}

fn validate_contribution(
    contribution: &MechanicsContribution,
    world_frame: &StableId,
) -> Result<(), MechanicsMacroStepError> {
    let identity = contribution.port.identity().clone();
    if !matches!(
        contribution.port.channel(),
        EulerChannel::NormalContact
            | EulerChannel::TangentialContact
            | EulerChannel::RollingContourSpin
            | EulerChannel::Impact
    ) {
        return Err(MechanicsMacroStepError::UnsupportedActiveChannel {
            channel: contribution.port.channel(),
        });
    }
    if !contribution.port.activity().is_active() {
        return Err(MechanicsMacroStepError::InactiveContribution { identity });
    }
    if contribution.body_action != BodyAction::ActionOnIntegratedBody {
        return Err(MechanicsMacroStepError::ReactionCannotActOnIntegratedBody { identity });
    }
    if contribution.reference != MechanicsReference::CenterOfMass {
        return Err(MechanicsMacroStepError::UnsupportedReference { identity });
    }
    if contribution.port.domain().coordinate().binding().frame() != world_frame {
        return Err(MechanicsMacroStepError::WrongWorldFrame { identity });
    }
    if contribution
        .port
        .domain()
        .coordinate()
        .binding()
        .orientation()
        != PortOrientation::AlongFrame
    {
        return Err(MechanicsMacroStepError::WrongCoordinateOrientation { identity });
    }
    let expected_port_kind = if contribution.port.channel() == EulerChannel::RollingContourSpin {
        PortKind::RotationalTorqueAngularVelocity
    } else {
        PortKind::MechanicalForceVelocity
    };
    if contribution.port.port_kind() != expected_port_kind {
        return Err(MechanicsMacroStepError::WrongPortKind { identity });
    }
    if !matches!(
        contribution.port.ownership(),
        ContributionOwnership::Exclusive
    ) {
        return Err(MechanicsMacroStepError::OverlappingOwnership { identity });
    }
    if !contribution.application_arm_body_m.is_finite()
        || !contribution.force_world_n.is_finite()
        || !contribution.free_torque_world_nm.is_finite()
        || !all_finite(&[
            contribution.claimed_discrete_work_j,
            contribution.claimed_storage_j,
            contribution.claimed_dissipation_j,
            contribution.claimed_heat_j,
            contribution.uncertainty_j,
        ])
        || contribution.claimed_dissipation_j < 0.0
        || contribution.claimed_heat_j < 0.0
        || contribution.claimed_heat_j > contribution.claimed_dissipation_j
        || contribution.uncertainty_j < 0.0
    {
        return Err(MechanicsMacroStepError::InvalidInput(
            "mechanics contribution",
        ));
    }
    Ok(())
}

fn reconstruct_contribution(
    contribution: &MechanicsContribution,
    step: &StepReceipt,
    properties: MassProperties,
    duration_seconds: f64,
) -> Result<ContributionMechanicsReceipt, MechanicsMacroStepError> {
    let before = step
        .state_before
        .point_kinematics(properties, contribution.application_arm_body_m)
        .map_err(MechanicsMacroStepError::Dynamics)?;
    let after = step
        .state_after
        .point_kinematics(properties, contribution.application_arm_body_m)
        .map_err(MechanicsMacroStepError::Dynamics)?;
    // This endpoint-centred reconstruction is the independently evaluated
    // discrete midpoint used by this composition layer.  It intentionally
    // does not claim a separate collocation state or a variational integrator.
    let midpoint_point_world_m = before.point_world.add(after.point_world).scale(0.5);
    let midpoint_arm_world_m = before.arm_world.add(after.arm_world).scale(0.5);
    let midpoint_point_velocity_world_m_per_s = before
        .point_velocity_world
        .add(after.point_velocity_world)
        .scale(0.5);
    let midpoint_angular_velocity_world_rad_per_s = before
        .angular_velocity_world
        .add(after.angular_velocity_world)
        .scale(0.5);
    let midpoint_torque_center_of_mass_world_nm = midpoint_arm_world_m
        .cross(contribution.force_world_n)
        .add(contribution.free_torque_world_nm);
    let power_w = contribution
        .force_world_n
        .dot(midpoint_point_velocity_world_m_per_s)
        + contribution
            .free_torque_world_nm
            .dot(midpoint_angular_velocity_world_rad_per_s);
    let reconstructed_discrete_work_j = power_w * duration_seconds;
    if !midpoint_point_world_m.is_finite()
        || !midpoint_arm_world_m.is_finite()
        || !midpoint_point_velocity_world_m_per_s.is_finite()
        || !midpoint_angular_velocity_world_rad_per_s.is_finite()
        || !midpoint_torque_center_of_mass_world_nm.is_finite()
        || !reconstructed_discrete_work_j.is_finite()
    {
        return Err(MechanicsMacroStepError::InvalidInput(
            "reconstructed contribution mechanics",
        ));
    }
    Ok(ContributionMechanicsReceipt {
        identity: contribution.port.identity().clone(),
        channel: contribution.port.channel(),
        midpoint_point_world_m,
        midpoint_arm_world_m,
        midpoint_point_velocity_world_m_per_s,
        midpoint_angular_velocity_world_rad_per_s,
        midpoint_torque_center_of_mass_world_nm,
        reconstructed_discrete_work_j,
        claimed_work_residual_j: contribution.claimed_discrete_work_j
            - reconstructed_discrete_work_j,
    })
}

fn energy_scale(first_j: f64, second_j: f64, uncertainty_j: f64) -> f64 {
    first_j.abs().max(second_j.abs()).max(uncertainty_j)
}

fn combined_energy_residual_j(
    body_energy_residual_j: f64,
    interface_energy_residual_j: f64,
) -> f64 {
    body_energy_residual_j - interface_energy_residual_j
}

fn all_finite(values: &[f64]) -> bool {
    values.iter().all(|value| value.is_finite())
}

#[cfg(test)]
mod tests {
    use super::{combined_energy_residual_j, energy_scale};

    #[test]
    fn combined_residual_subtracts_nonzero_interface_residual() {
        let mechanical_energy_change_j = -0.12;
        let reconstructed_body_work_j = -0.15;
        let claimed_body_work_j = -0.15;
        let recoverable_storage_change_j = 0.10;
        let irreversible_work_j = 0.07;
        let body_energy_residual_j = mechanical_energy_change_j - reconstructed_body_work_j;
        let interface_energy_residual_j =
            -claimed_body_work_j - (recoverable_storage_change_j + irreversible_work_j);

        assert_ne!(body_energy_residual_j, 0.0);
        assert_ne!(interface_energy_residual_j, 0.0);
        let aggregate_residual_j =
            combined_energy_residual_j(body_energy_residual_j, interface_energy_residual_j);
        let physical_full_system_residual_j =
            mechanical_energy_change_j + recoverable_storage_change_j + irreversible_work_j;
        assert!(
            (aggregate_residual_j - physical_full_system_residual_j).abs() <= 1.0e-15,
            "aggregate={aggregate_residual_j:?}, physical={physical_full_system_residual_j:?}"
        );
    }

    #[test]
    fn energy_scale_keeps_sub_joule_relative_tolerances_sub_joule() {
        assert_eq!(energy_scale(0.001, 0.0002, 0.0), 0.001);
    }
}
