//! Native geometry dispatch; the motion and impact owners remain shared.

use super::{
    FiniteSolidSphere, FiniteSolidSphereScene, FiniteSolidSphereSceneKey,
    MixedViscoSpherePathRefusal, MixedViscoSphereTrajectoryError, QuadraticSolidContactSurface,
    QuadraticSolidSurfaceError, QuadraticSphereSceneError,
};
use crate::inelastic_impact::{
    multiple::InelasticMixedContact, surface::quadratic::QuadraticSolidSurfaceControls,
};
use fs_exec::Cx;
use fs_mbd::RigidBodyState;
use fs_solid::hyper3::mixed::visco::dynamics::Tet10ImpulsePoint;

#[derive(Clone, Copy)]
pub(super) enum TrajectoryScene<'s> {
    /// No rigid-pair law was supplied; possible encounters must refuse.
    Surface(&'s [FiniteSolidSphere]),
    /// Every unordered rigid pair has an explicit, independently resolved law.
    Complete(&'s FiniteSolidSphereScene),
}

pub(super) struct TrajectoryGeometry {
    pub contacts: Vec<InelasticMixedContact<Tet10ImpulsePoint>>,
    pub keys: Vec<FiniteSolidSphereSceneKey>,
    pub maximum_contact_gap_m: f64,
}

impl<'s> TrajectoryScene<'s> {
    pub(super) fn spheres(self) -> &'s [FiniteSolidSphere] {
        match self {
            Self::Surface(spheres) => spheres,
            Self::Complete(scene) => scene.spheres(),
        }
    }

    pub(super) fn contacts(
        self,
        surface: &QuadraticSolidContactSurface,
        positions: &[[f64; 3]],
        bodies: &[RigidBodyState],
        controls: QuadraticSolidSurfaceControls,
        cx: &Cx<'_>,
    ) -> Result<TrajectoryGeometry, MixedViscoSphereTrajectoryError> {
        use MixedViscoSphereTrajectoryError as E;
        match self {
            Self::Complete(scene) => {
                let query = scene
                    .quadratic_contacts_at_positions(surface, positions, bodies, controls, cx)
                    .map_err(E::Scene)?;
                Ok(TrajectoryGeometry {
                    contacts: query.contacts,
                    keys: query.keys,
                    maximum_contact_gap_m: query.maximum_contact_gap_m,
                })
            }
            Self::Surface(spheres) => {
                let query = surface
                    .sphere_separations_at_positions(positions, bodies, spheres, controls, cx)
                    .map_err(E::Surface)?;
                let contacts = query.admitted_contacts().map_err(E::Surface)?;
                let keys = query
                    .contacts
                    .iter()
                    .map(|w| FiniteSolidSphereSceneKey::Solid {
                        sphere: w.sphere_id,
                        contact_id: w.contact.id,
                        face: w.face,
                    })
                    .collect();
                let maximum_contact_gap_m =
                    query.contacts.iter().map(|w| w.gap_m).fold(0.0, f64::max);
                Ok(TrajectoryGeometry {
                    contacts,
                    keys,
                    maximum_contact_gap_m,
                })
            }
        }
    }
}

/// Only a bounded endpoint localization ambiguity can shorten a trial here.
/// Penetration, failed source admission and exhausted query work remain errors.
pub(super) fn refinable_endpoint(
    error: MixedViscoSphereTrajectoryError,
) -> Result<MixedViscoSpherePathRefusal, MixedViscoSphereTrajectoryError> {
    use MixedViscoSphereTrajectoryError as E;
    match error {
        E::Surface(
            error @ (QuadraticSolidSurfaceError::UnresolvedGap { .. }
            | QuadraticSolidSurfaceError::UnresolvedActivation { .. }),
        )
        | E::Scene(QuadraticSphereSceneError::Surface(
            error @ (QuadraticSolidSurfaceError::UnresolvedGap { .. }
            | QuadraticSolidSurfaceError::UnresolvedActivation { .. }),
        )) => Ok(MixedViscoSpherePathRefusal::Endpoint(error)),
        error => Err(error),
    }
}
