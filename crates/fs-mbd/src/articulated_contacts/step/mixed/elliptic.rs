//! Directional source friction through the original mixed force/contact/drift.
//!
//! Native load kicks, free/prescribed drift, whole-path geometry and every
//! work/momentum gate remain in one owner. Original source histories are
//! sampled at the actual post-kick contact state and stay private until the
//! entire physical step, including its drift and path, is accepted.

use super::{
    ArticulatedContactStepConfig, ArticulatedContactStepError, MixedArticulatedContactPath,
    MixedArticulatedContactStepResult, MixedArticulatedTrajectorySystem, PrescribedBaseMotionError,
    PrescribedBaseMotionState, PrescribedBasePath, contact_error, step_with,
};
use crate::articulated_contacts::friction::ArticulatedImpactKinematics;
use crate::articulated_contacts::friction::elliptic::{
    ArticulatedEllipticContact, ArticulatedEllipticSourceError, ArticulatedImpactGeometry,
    ArticulatedSourceCandidate,
};
use crate::articulated_contacts::mixed::{
    resolve_mixed_articulated_elliptic_contacts, resolve_mixed_articulated_elliptic_contacts_with,
};
use fs_contact::frictional_multicontact::elliptic::{
    EllipticFriction, EllipticFrictionalContactResult,
};

/// Complete native step with physical directional contact impulses.
pub type MixedArticulatedEllipticContactStepResult<P> =
    MixedArticulatedContactStepResult<P, EllipticFrictionalContactResult>;

/// Complete mechanical and source candidates, published together by the caller.
pub struct MixedArticulatedEllipticSourceStepResult<S, P> {
    pub mechanical: MixedArticulatedEllipticContactStepResult<P>,
    /// Active source states in canonical contact-ID order, with original indices.
    pub candidates: Vec<ArticulatedSourceCandidate<S>>,
}

/// No mechanical endpoint or source candidate is returned on either refusal.
#[derive(Debug)]
pub enum MixedArticulatedEllipticSourceStepError<E, SE> {
    Step(ArticulatedContactStepError<E>),
    Source {
        contact_id: u64,
        original_index: usize,
        source: SE,
    },
}
impl<E, SE> From<ArticulatedContactStepError<E>>
    for MixedArticulatedEllipticSourceStepError<E, SE>
{
    fn from(error: ArticulatedContactStepError<E>) -> Self {
        Self::Step(error)
    }
}
impl<E: core::fmt::Debug, SE: core::fmt::Debug> core::fmt::Display
    for MixedArticulatedEllipticSourceStepError<E, SE>
{
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "source directional mixed contact step: {self:?}")
    }
}
impl<E: core::fmt::Debug, SE: core::fmt::Debug> std::error::Error
    for MixedArticulatedEllipticSourceStepError<E, SE>
{
}

/// Advance loaded mixed trees with an explicit oriented Coulomb ellipse.
///
/// The original coupled contact law retains normal/tangent and joint/base
/// response. Equal principal coefficients recover the circular law. Native
/// force work, all prescribed-base/plane work, unloaded drift, whole-path
/// geometry and global momentum/work admission retain their existing meaning.
///
/// Coefficients are held for this first-order step. Its force/contact splitting
/// loss is not a continuum heat allocation. Normal response, geometry coverage,
/// support paths and timestep refinement remain explicit caller inputs.
#[allow(clippy::too_many_arguments)]
pub fn step_mixed_articulated_elliptic_contacts<P, E>(
    systems: &[MixedArticulatedTrajectorySystem<'_>],
    contacts: &[ArticulatedEllipticContact],
    start_time_s: f64,
    config: ArticulatedContactStepConfig,
    support_path: impl FnMut(
        usize,
        &PrescribedBaseMotionState,
        f64,
    ) -> Result<PrescribedBasePath, PrescribedBaseMotionError>,
    admit_path: impl FnOnce(&MixedArticulatedContactPath<'_, EllipticFriction>) -> Result<P, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<MixedArticulatedEllipticContactStepResult<P>, ArticulatedContactStepError<E>> {
    step_with(
        systems,
        contacts,
        start_time_s,
        config,
        |inputs, contacts, controls, cancelled| {
            let contact =
                resolve_mixed_articulated_elliptic_contacts(inputs, contacts, controls, cancelled)
                    .map_err(contact_error::<E>)?;
            let loss = contact.contact.dissipated_energy;
            Ok((contact, loss))
        },
        support_path,
        admit_path,
        cancelled,
    )
}

/// Resolve original source histories at the actual native post-load-kick state.
///
/// Histories and geometry retain original input order. The native event owner
/// samples each active source once in stable contact-ID order, using its actual
/// incoming point velocities, assembled frame and physical witness positions.
/// Inactive histories remain unchanged. The callback owns source conditions and
/// objective texture transport; it receives no invented pressure or slip path.
///
/// Candidates remain private through native contact, prescribed/free drift,
/// endpoint gaps, complete path admission, physical balances and final
/// cancellation. Histories need no Clone bound. A late refusal returns neither
/// a mechanical endpoint nor source candidates. Callback-owned external state
/// must also be staged by its owner.
///
/// This is the same held-coefficient first-order discretization as the supplied
/// law entry, not a continuous source iteration, finite-slip history evolution
/// or thermal partition. The caller publishes every returned state together.
#[allow(clippy::too_many_arguments)]
pub fn step_mixed_articulated_elliptic_contacts_with<S, P, E, SE>(
    systems: &[MixedArticulatedTrajectorySystem<'_>],
    contacts: &[ArticulatedImpactGeometry],
    histories: &[S],
    start_time_s: f64,
    config: ArticulatedContactStepConfig,
    mut resolve: impl FnMut(&S, &ArticulatedImpactKinematics) -> Result<(EllipticFriction, S), SE>,
    support_path: impl FnMut(
        usize,
        &PrescribedBaseMotionState,
        f64,
    ) -> Result<PrescribedBasePath, PrescribedBaseMotionError>,
    admit_path: impl FnOnce(&MixedArticulatedContactPath<'_, ()>) -> Result<P, E>,
    cancelled: impl FnMut() -> bool,
) -> Result<
    MixedArticulatedEllipticSourceStepResult<S, P>,
    MixedArticulatedEllipticSourceStepError<E, SE>,
> {
    if histories.len() != contacts.len() {
        return Err(ArticulatedContactStepError::<E>::Invalid(
            "one accepted source history per original contact required",
        )
        .into());
    }
    let mut candidates = Vec::new();
    let mechanical = step_with(
        systems,
        contacts,
        start_time_s,
        config,
        |inputs,
         contacts,
         controls,
         cancelled|
         -> Result<_, MixedArticulatedEllipticSourceStepError<E, SE>> {
            let source = resolve_mixed_articulated_elliptic_contacts_with(
                inputs,
                contacts,
                histories,
                controls,
                &mut resolve,
                cancelled,
            )
            .map_err(source_error::<E, SE>)?;
            let loss = source.mechanical.contact.dissipated_energy;
            candidates = source.candidates;
            Ok((source.mechanical, loss))
        },
        support_path,
        admit_path,
        cancelled,
    )?;
    Ok(MixedArticulatedEllipticSourceStepResult {
        mechanical,
        candidates,
    })
}

fn source_error<E, SE>(
    error: ArticulatedEllipticSourceError<SE>,
) -> MixedArticulatedEllipticSourceStepError<E, SE> {
    match error {
        ArticulatedEllipticSourceError::Mechanics(error) => {
            MixedArticulatedEllipticSourceStepError::Step(contact_error(error))
        }
        ArticulatedEllipticSourceError::Source {
            contact_id,
            original_index,
            source,
        } => MixedArticulatedEllipticSourceStepError::Source {
            contact_id,
            original_index,
            source,
        },
    }
}
