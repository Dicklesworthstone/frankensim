//! Simultaneous normal impacts with free and prescribed articulated boundaries.
//!
//! Every tree uses its original native inverse-mass action. Prescribed supports
//! retain their full reaction impulse and signed work in the shared physical
//! endpoint audit; they are never represented by a large fictitious mass.

use super::{
    ArticulatedContactError, ArticulatedContactResult, ArticulatedContactSystem,
    ArticulatedNormalContact, resolve_contacts,
};
use crate::articulated::{
    ArticulatedError, ArticulatedModel, ArticulatedMomentum, BaseState, FreeFloatingBaseState,
    FreeFloatingImpulseResponse, ImpulseVelocityChange, PrescribedBaseImpulseResponse,
    apply_free_floating_impulses, apply_prescribed_base_impulses, free_floating_impulse_delta,
    prescribed_base_impulse_delta,
};
use fs_contact::multicontact::MultiImpactConfig;
use fs_ga::{Twist, Vec3, Wrench};

/// Physical boundary on one original articulated tree during the impulse.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ArticulatedContactBoundary {
    /// Six base velocities recoil with every available joint.
    Free(FreeFloatingBaseState),
    /// Base pose and twist remain prescribed; the support receives the reaction.
    Prescribed(BaseState),
}

/// One tree in a simultaneous mixed-boundary normal contact event.
#[derive(Debug, Clone, Copy)]
pub struct MixedArticulatedContactSystem<'a> {
    /// Original validated model, including full inertia tensors and root joints.
    pub model: &'a ArticulatedModel,
    /// Actual boundary state at the supplied event configuration.
    pub boundary: ArticulatedContactBoundary,
    /// Joint configuration, unchanged by the event.
    pub joint_position: &'a [f64],
    /// Incoming joint velocities in compact model order.
    pub joint_velocity: &'a [f64],
}

/// Actual native response of one tree, in input-system order.
#[derive(Debug, Clone, PartialEq)]
pub enum MixedArticulatedImpulseResponse {
    /// Original free-tree response with all base and joint recoil.
    Free(FreeFloatingImpulseResponse),
    /// Original supported-tree response, including support impulse and work.
    Prescribed(PrescribedBaseImpulseResponse),
}

/// Resolve simultaneous frictionless contacts across free and supported trees.
///
/// The existing normal-contact geometry, coupled restitution law, accumulated
/// impulses and physical endpoint audit are shared with the free-tree API.
/// Each probe selects the native free or prescribed ABA response, preserving
/// cross-contact coupling, full tensors, offset COMs and scalar root joints.
/// All point pairs use equal and opposite impulses. Support reactions are
/// external impulses; their signed work is separate from moving-plane work.
///
/// The same bounds of 64 contacts and 128 trees apply. Singular or immobile
/// active contact directions refuse through the original positive-definite
/// contact solver. No arbitrary reactions, support compliance, position repair,
/// friction, collision detection or continuous trajectory are inferred.
/// All borrowed input states remain unchanged on refusal or cancellation.
pub fn resolve_mixed_articulated_contacts(
    systems: &[MixedArticulatedContactSystem<'_>],
    contacts: &[ArticulatedNormalContact],
    config: MultiImpactConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedContactResult<MixedArticulatedImpulseResponse>, ArticulatedContactError> {
    resolve_contacts(systems, contacts, config, cancelled)
}

// Views share the physical assembly without changing the original public free
// state/response types or cloning their joint vectors.
pub(super) struct SystemView<'a> {
    pub model: &'a ArticulatedModel,
    pub base: BaseState,
    pub joint_position: &'a [f64],
    pub joint_velocity: &'a [f64],
}

pub(super) struct ImpulseView<'a> {
    pub base_after: BaseState,
    pub joint_velocity_after: &'a [f64],
    pub momentum_before: ArticulatedMomentum,
    pub momentum_after: ArticulatedMomentum,
    pub kinetic_energy_change: f64,
    pub linear_momentum_change_world: Vec3,
    pub angular_momentum_change_world: Vec3,
    pub support: Option<(Wrench, f64)>,
}

pub(super) trait ContactResponse {
    fn impulse_view(&self) -> ImpulseView<'_>;
}

pub(super) trait ContactSystem {
    type Response: ContactResponse;
    fn system_view(&self) -> SystemView<'_>;
    fn inverse_mass(
        &self,
        joint_impulses: &[f64],
        body_impulses: &[Wrench],
    ) -> Result<ImpulseVelocityChange, ArticulatedError>;
    fn apply(
        &self,
        joint_impulses: &[f64],
        body_impulses: &[Wrench],
    ) -> Result<Self::Response, ArticulatedError>;
}

impl ContactResponse for FreeFloatingImpulseResponse {
    fn impulse_view(&self) -> ImpulseView<'_> {
        ImpulseView {
            base_after: BaseState::prescribed(
                self.base_after.world_from_base,
                self.base_after.twist_body,
                Twist::zero(),
            ),
            joint_velocity_after: &self.joint_velocity_after,
            momentum_before: self.momentum_before,
            momentum_after: self.momentum_after,
            kinetic_energy_change: self.kinetic_energy_change,
            linear_momentum_change_world: self.linear_momentum_change_world,
            angular_momentum_change_world: self.angular_momentum_change_world,
            support: None,
        }
    }
}

impl ContactResponse for PrescribedBaseImpulseResponse {
    fn impulse_view(&self) -> ImpulseView<'_> {
        ImpulseView {
            base_after: self.base_after,
            joint_velocity_after: &self.joint_velocity_after,
            momentum_before: self.momentum_before,
            momentum_after: self.momentum_after,
            kinetic_energy_change: self.kinetic_energy_change,
            linear_momentum_change_world: self.linear_momentum_change_world,
            angular_momentum_change_world: self.angular_momentum_change_world,
            support: Some((self.support_impulse_world, self.support_work)),
        }
    }
}

impl ContactResponse for MixedArticulatedImpulseResponse {
    fn impulse_view(&self) -> ImpulseView<'_> {
        match self {
            Self::Free(response) => response.impulse_view(),
            Self::Prescribed(response) => response.impulse_view(),
        }
    }
}

impl ContactSystem for ArticulatedContactSystem<'_> {
    type Response = FreeFloatingImpulseResponse;
    fn system_view(&self) -> SystemView<'_> {
        SystemView {
            model: self.model,
            base: BaseState::prescribed(
                self.base.world_from_base,
                self.base.twist_body,
                Twist::zero(),
            ),
            joint_position: self.joint_position,
            joint_velocity: self.joint_velocity,
        }
    }
    fn inverse_mass(
        &self,
        joint_impulses: &[f64],
        body_impulses: &[Wrench],
    ) -> Result<ImpulseVelocityChange, ArticulatedError> {
        free_floating_impulse_delta(
            self.model,
            self.joint_position,
            joint_impulses,
            body_impulses,
        )
    }
    fn apply(
        &self,
        joint_impulses: &[f64],
        body_impulses: &[Wrench],
    ) -> Result<Self::Response, ArticulatedError> {
        apply_free_floating_impulses(
            self.model,
            self.base,
            self.joint_position,
            self.joint_velocity,
            joint_impulses,
            body_impulses,
        )
    }
}

impl ContactSystem for MixedArticulatedContactSystem<'_> {
    type Response = MixedArticulatedImpulseResponse;
    fn system_view(&self) -> SystemView<'_> {
        SystemView {
            model: self.model,
            base: match self.boundary {
                ArticulatedContactBoundary::Free(base) => {
                    BaseState::prescribed(base.world_from_base, base.twist_body, Twist::zero())
                }
                ArticulatedContactBoundary::Prescribed(base) => base,
            },
            joint_position: self.joint_position,
            joint_velocity: self.joint_velocity,
        }
    }
    fn inverse_mass(
        &self,
        joint_impulses: &[f64],
        body_impulses: &[Wrench],
    ) -> Result<ImpulseVelocityChange, ArticulatedError> {
        match self.boundary {
            ArticulatedContactBoundary::Free(_) => free_floating_impulse_delta(
                self.model,
                self.joint_position,
                joint_impulses,
                body_impulses,
            ),
            ArticulatedContactBoundary::Prescribed(_) => prescribed_base_impulse_delta(
                self.model,
                self.joint_position,
                joint_impulses,
                body_impulses,
            )
            .map(|delta| delta.velocity),
        }
    }
    fn apply(
        &self,
        joint_impulses: &[f64],
        body_impulses: &[Wrench],
    ) -> Result<Self::Response, ArticulatedError> {
        match self.boundary {
            ArticulatedContactBoundary::Free(base) => apply_free_floating_impulses(
                self.model,
                base,
                self.joint_position,
                self.joint_velocity,
                joint_impulses,
                body_impulses,
            )
            .map(MixedArticulatedImpulseResponse::Free),
            ArticulatedContactBoundary::Prescribed(base) => apply_prescribed_base_impulses(
                self.model,
                base,
                self.joint_position,
                self.joint_velocity,
                joint_impulses,
                body_impulses,
            )
            .map(MixedArticulatedImpulseResponse::Prescribed),
        }
    }
}
