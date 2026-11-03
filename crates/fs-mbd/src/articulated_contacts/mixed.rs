//! Simultaneous normal and Coulomb impacts with native articulated boundaries.
//!
//! Every tree uses its original native inverse-mass action. Prescribed supports
//! retain their full reaction impulse and signed work in the shared physical
//! endpoint audit; they are never represented by a large fictitious mass.

use super::friction::elliptic::{
    ArticulatedEllipticContact, ArticulatedEllipticContactResult, ArticulatedEllipticSourceError,
    ArticulatedEllipticSourceResult, ArticulatedImpactGeometry, resolve_elliptic_event,
    resolve_elliptic_event_with,
};
use super::friction::{
    ArticulatedFrictionalContact, ArticulatedFrictionalContactResult, ArticulatedImpactKinematics,
    resolve_frictional_event,
};
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
use fs_contact::frictional_multicontact::elliptic::EllipticFriction;
use fs_contact::frictional_multicontact::{FrictionalContactConfig, FrictionalContactResult};
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
///
/// A unique-velocity law determines motion and total boundary work without
/// identifying individual contact/support loads or source-work contributions.
/// Prescribed reactions retain the selected impulse representative's full six
/// components and signed work, in their original native frames.
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
/// The same bounds of 64 contacts and 128 trees apply. Redundant rows can
/// determine physical motion with representative impulses; the original law
/// separately admits total work. Locked directions or unresolved motion/work
/// refuse. No support compliance, position repair, friction, collision
/// detection or continuous trajectory is inferred by this normal entry.
/// All borrowed input states remain unchanged on refusal or cancellation.
pub fn resolve_mixed_articulated_contacts(
    systems: &[MixedArticulatedContactSystem<'_>],
    contacts: &[ArticulatedNormalContact],
    config: MultiImpactConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<ArticulatedContactResult<MixedArticulatedImpulseResponse>, ArticulatedContactError> {
    resolve_contacts(systems, contacts, config, cancelled)
}

impl MixedArticulatedImpulseResponse {
    /// Unchanged configuration with the actually returned base twist.
    #[must_use]
    pub fn base_after(&self) -> ArticulatedContactBoundary {
        match self {
            Self::Free(response) => ArticulatedContactBoundary::Free(response.base_after),
            Self::Prescribed(response) => {
                ArticulatedContactBoundary::Prescribed(response.base_after)
            }
        }
    }

    /// Actual returned generalized velocities in the original compact order.
    #[must_use]
    pub fn joint_velocity_after(&self) -> &[f64] {
        self.impulse_view().joint_velocity_after
    }
}

/// Resolve simultaneous Coulomb impacts over mixed free and prescribed bases.
///
/// Each tree uses its original constrained or free ABA mobility, retaining all
/// cross-contact and normal/tangent coupling. Prescribed supports transmit full
/// force/couple impulses and signed work; free trees recoil in the same global
/// solve. Native endpoint kinematics, momentum, energy and limits are admitted
/// before any complete candidate publishes.
///
/// The original law distinguishes unique impulses from unique motion and total
/// work. In the latter case individual contact/support loads and source-work
/// contributions remain representative. This adds no support compliance,
/// contact detection, force duration, or frictional heat partition.
pub fn resolve_mixed_articulated_frictional_contacts(
    systems: &[MixedArticulatedContactSystem<'_>],
    contacts: &[ArticulatedFrictionalContact],
    config: FrictionalContactConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<
    ArticulatedFrictionalContactResult<FrictionalContactResult, MixedArticulatedImpulseResponse>,
    ArticulatedContactError,
> {
    resolve_frictional_event(systems, contacts, config, cancelled)
}

/// Resolve directional Coulomb impacts over free and prescribed-base trees.
///
/// Principal friction axes remain in each original material contact frame.
/// The same native constrained mobility, full support reactions, signed total
/// boundary work and actual endpoint audits as the circular mixed-base entry
/// point apply. A unique-velocity witness does not identify individual support
/// loads or their source-work allocation.
pub fn resolve_mixed_articulated_elliptic_contacts(
    systems: &[MixedArticulatedContactSystem<'_>],
    contacts: &[ArticulatedEllipticContact],
    config: FrictionalContactConfig,
    cancelled: impl FnMut() -> bool,
) -> Result<
    ArticulatedEllipticContactResult<MixedArticulatedImpulseResponse>,
    ArticulatedContactError,
> {
    resolve_elliptic_event(systems, contacts, config, cancelled)
}

/// Resolve source-owned directional friction over free and prescribed mechanisms.
///
/// Each geometrically active contact samples its borrowed accepted source history
/// once, in canonical contact-ID order, at the original material-point locations,
/// contact axes and relative velocity. Inactive rows retain native endpoint
/// kinematics but neither query nor advance source history. The returned source
/// candidates retain their original input indices; histories need not be Clone.
///
/// The existing source transaction and mixed mechanical owner remain responsible
/// for admission: source failure, joint-speed refusal, unresolved momentum/work
/// or cancellation returns neither mechanical nor history candidates. The caller
/// publishes all accepted owners together; callback-owned external side effects
/// must be staged by their owner. Native total-boundary-work determinacy is
/// retained without identifying redundant local loads or individual work shares.
///
/// Coefficients are held at their preimpact query. This adds no implicit
/// velocity-dependent constitutive iteration, contact discovery or thermal solve.
pub fn resolve_mixed_articulated_elliptic_contacts_with<S, E>(
    systems: &[MixedArticulatedContactSystem<'_>],
    contacts: &[ArticulatedImpactGeometry],
    histories: &[S],
    config: FrictionalContactConfig,
    resolve: impl FnMut(&S, &ArticulatedImpactKinematics) -> Result<(EllipticFriction, S), E>,
    cancelled: impl FnMut() -> bool,
) -> Result<
    ArticulatedEllipticSourceResult<S, MixedArticulatedImpulseResponse>,
    ArticulatedEllipticSourceError<E>,
> {
    resolve_elliptic_event_with(systems, contacts, histories, config, resolve, cancelled)
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
