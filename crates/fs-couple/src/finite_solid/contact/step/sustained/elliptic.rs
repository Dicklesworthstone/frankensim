//! Directional coefficients through the shared native sustained-contact path.
//!
//! The original contact owner transforms the ellipse by work conjugacy, solves
//! the complete mixed system and returns PHYSICAL impulses to native bodies.
//! The complete transformed determinacy receipt is retained, while force work,
//! angular transport, consistent solid response, geometry, refinement and
//! publication remain in the same trajectory core as the circular entry point.

use super::{
    FiniteSolidSustainedContactContext, FiniteSolidSustainedTrajectoryError,
    FiniteSolidSustainedTrajectoryReport,
};
use crate::finite_solid::FiniteSolidBody;
use crate::finite_solid::contact::mixed::{
    FiniteSolidMixedContact, FiniteSolidMixedEllipticContactError,
};
use crate::finite_solid::contact::scene::{FiniteSolidSphereScene, FiniteSolidSphereSceneKey};
use crate::finite_solid::contact::step::FiniteSolidSphereTrajectoryControls;
use crate::finite_solid::contact::surface::FiniteSolidContactSurface;
use fs_contact::frictional_multicontact::elliptic::{
    EllipticFriction, EllipticFrictionalContactResult,
};
use fs_exec::Cx;
use fs_mbd::{Gravity, MassProperties, RigidBodyState, Wrench};

pub type FiniteSolidSustainedEllipticTrajectoryReport =
    FiniteSolidSustainedTrajectoryReport<EllipticFrictionalContactResult>;
pub type FiniteSolidSustainedEllipticTrajectoryError<E> = FiniteSolidSustainedTrajectoryError<
    EllipticFrictionalContactResult,
    FiniteSolidSustainedEllipticFailure<E>,
>;

/// Retains original scene identity and the complete original mechanical refusal.
#[derive(Debug)]
pub enum FiniteSolidSustainedEllipticFailure<E> {
    Coefficients {
        key: FiniteSolidSphereSceneKey,
        original_index: usize,
        error: E,
    },
    Contact(FiniteSolidMixedEllipticContactError),
}

impl<'a> FiniteSolidBody<'a> {
    /// Advance actual mixed sphere contact with explicitly oriented ellipse laws.
    ///
    /// Resolve each current scene row in canonical contact-ID order using the
    /// ACTUAL post-force-kick states. The original contact frame is retained;
    /// principal coefficients and axis are interpreted in that frame. Newton
    /// targets remain those already derived from pre-force incoming velocity.
    /// Equal principal coefficients recover circular friction; a singular
    /// one-axis ellipse retains the original contact owner's typed refusal.
    ///
    /// Coefficients are held over each proposed interval, including its shared
    /// original impulse solve and native kinetic drift. A rejected interval may
    /// call the resolver again with a smaller duration and different kicked
    /// velocity. The resolver must not publish external constitutive history;
    /// source-history evolution uses the separately staged source owner.
    /// No coefficient averaging, surrogate mass/inertia, or heat allocation occurs.
    #[allow(clippy::too_many_arguments)]
    pub fn sustained_elliptic_sphere_trajectory_candidate<E>(
        &self,
        surface: &FiniteSolidContactSurface,
        scene: &FiniteSolidSphereScene,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        wrenches: &[Wrench],
        gravity: Gravity,
        duration_s: f64,
        origin_m: [f64; 3],
        controls: FiniteSolidSphereTrajectoryControls,
        mut friction: impl FnMut(
            &FiniteSolidSustainedContactContext<'_, 'a, '_>,
            &FiniteSolidSphereSceneKey,
            &FiniteSolidMixedContact,
        ) -> Result<EllipticFriction, E>,
        cx: &Cx<'_>,
    ) -> Result<
        (
            Self,
            Vec<RigidBodyState>,
            FiniteSolidSustainedEllipticTrajectoryReport,
        ),
        FiniteSolidSustainedEllipticTrajectoryError<E>,
    > {
        let (solid, bodies, (), report) = self.sustained_sphere_trajectory_with(
            surface,
            scene,
            bodies,
            properties,
            nodal_forces_n,
            body_force_n_m3,
            wrenches,
            gravity,
            duration_s,
            origin_m,
            controls,
            &(),
            |context, _| {
                let mut contacts = Vec::with_capacity(context.contacts.len());
                for (original_index, (&contact, key)) in context
                    .contacts
                    .iter()
                    .zip(context.contact_keys)
                    .enumerate()
                {
                    let coefficient = friction(&context, key, &contact).map_err(|error| {
                        FiniteSolidSustainedEllipticFailure::Coefficients {
                            key: *key,
                            original_index,
                            error,
                        }
                    })?;
                    contacts.push(contact.with_friction(coefficient));
                }
                let (solid, bodies, report) = context
                    .solid
                    .mixed_elliptic_contact_candidate(
                        context.bodies,
                        context.properties,
                        &contacts,
                        context.origin_m,
                        context.controls,
                        context.cx,
                    )
                    .map_err(FiniteSolidSustainedEllipticFailure::Contact)?;
                Ok((solid, bodies, report, ()))
            },
            |_, _, _| Ok(None),
            cx,
        )?;
        Ok((solid, bodies, report))
    }
}
