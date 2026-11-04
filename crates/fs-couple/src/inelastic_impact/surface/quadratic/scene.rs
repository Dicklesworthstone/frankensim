//! Original curved solid traces and rigid pairs in one canonical contact scene.
//!
//! The original prepared scene owns sphere identities, explicit pair laws,
//! disjoint contact IDs and rigid-pair witnesses. This adapter contributes the
//! admitted six-node boundary query and its exact native TET10 material point.
//! No linear surface receipt, corner projection or replacement point is made.

use super::{
    QuadraticSolidContactSurface, QuadraticSolidSurfaceControls, QuadraticSolidSurfaceError,
    QuadraticSolidSurfaceQuery,
};
use fs_exec::Cx;
use fs_mbd::RigidBodyState;
use fs_query::triangle::TriangleQueryError;
use fs_solid::hyper3::mixed::{
    conservative::ConservativeTetSolid,
    visco::dynamics::{MixedViscoTetMotion, Tet10ImpulsePoint},
};

/// Same prepared sphere topology, declared pair laws and canonical contact keys.
pub use crate::finite_solid::contact::scene::{
    FiniteSolidSphereScene, FiniteSolidSphereSceneError, FiniteSolidSphereSceneKey,
    FiniteSolidSphereSceneQuery, SpherePairLaw,
};

mod path;
pub use path::{
    QuadraticSphereScenePathKey, QuadraticSphereScenePathReport, QuadraticSphereScenePathVerdict,
};

/// Complete admitted contact rows with their original curved geometry receipt.
pub type QuadraticSphereSceneQuery =
    FiniteSolidSphereSceneQuery<Tet10ImpulsePoint, QuadraticSolidSurfaceQuery>;

/// Original surface, pair and path refusals retain their physical evidence.
#[derive(Debug)]
pub enum QuadraticSphereSceneError {
    Surface(QuadraticSolidSurfaceError),
    Scene(FiniteSolidSphereSceneError),
    Query(TriangleQueryError),
    PathBudget,
    Cancelled,
}
impl From<QuadraticSolidSurfaceError> for QuadraticSphereSceneError {
    fn from(error: QuadraticSolidSurfaceError) -> Self {
        Self::Surface(error)
    }
}
impl From<FiniteSolidSphereSceneError> for QuadraticSphereSceneError {
    fn from(error: FiniteSolidSphereSceneError) -> Self {
        Self::Scene(error)
    }
}
impl core::fmt::Display for QuadraticSphereSceneError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "quadratic sphere scene refused: {self:?}")
    }
}
impl std::error::Error for QuadraticSphereSceneError {}

impl FiniteSolidSphereScene {
    /// Detect all curved solid and rigid-pair rows at the same original indexed
    /// coordinates. Native callers first obtain this complete image from their
    /// exact reference-bound surface adapter; force kicks can reuse it unchanged.
    pub fn quadratic_contacts_at_positions(
        &self,
        surface: &QuadraticSolidContactSurface,
        positions: &[[f64; 3]],
        bodies: &[RigidBodyState],
        controls: QuadraticSolidSurfaceControls,
        cx: &Cx<'_>,
    ) -> Result<QuadraticSphereSceneQuery, QuadraticSphereSceneError> {
        validate_bodies(self, bodies, cx)?;
        let query = surface.sphere_separations_at_positions(
            positions,
            bodies,
            self.spheres(),
            controls,
            cx,
        )?;
        self.contacts_from_quadratic_query(query, bodies, controls, cx)
    }

    fn contacts_from_quadratic_query(
        &self,
        query: QuadraticSolidSurfaceQuery,
        bodies: &[RigidBodyState],
        controls: QuadraticSolidSurfaceControls,
        cx: &Cx<'_>,
    ) -> Result<QuadraticSphereSceneQuery, QuadraticSphereSceneError> {
        // Admission covers actual gaps, conservative lower bounds and every
        // unresolved activation candidate. The original signed P2 traces are
        // passed through unchanged to the common simultaneous response owner.
        let admitted = query.admitted_contacts()?;
        let contacts = admitted
            .into_iter()
            .zip(&query.contacts)
            .map(|(contact, witness)| (contact, witness.sphere_id, witness.face, witness.gap_m))
            .collect();
        Ok(self.contacts_from_admitted_surface(query, contacts, bodies, controls.surface, cx)?)
    }
}

macro_rules! native_scene {
    ($method:ident, $motion:ty, $query:ident) => {
        impl FiniteSolidSphereScene {
            /// Reuse the original native geometry and complete prepared scene.
            /// Every native reference coordinate/incidence must match before
            /// original material traces join the simultaneous pair response.
            pub fn $method(
                &self,
                surface: &QuadraticSolidContactSurface,
                solid: &$motion,
                bodies: &[RigidBodyState],
                controls: QuadraticSolidSurfaceControls,
                cx: &Cx<'_>,
            ) -> Result<QuadraticSphereSceneQuery, QuadraticSphereSceneError> {
                validate_bodies(self, bodies, cx)?;
                let query = surface.$query(solid, bodies, self.spheres(), controls, cx)?;
                self.contacts_from_quadratic_query(query, bodies, controls, cx)
            }
        }
    };
}

native_scene!(
    mixed_visco_contacts,
    MixedViscoTetMotion<'_>,
    mixed_visco_sphere_separations
);
native_scene!(
    conservative_contacts,
    ConservativeTetSolid<'_>,
    conservative_sphere_separations
);

fn validate_bodies(
    scene: &FiniteSolidSphereScene,
    bodies: &[RigidBodyState],
    cx: &Cx<'_>,
) -> Result<(), QuadraticSphereSceneError> {
    cx.checkpoint()
        .map_err(|_| QuadraticSphereSceneError::Cancelled)?;
    if bodies.len() != scene.spheres().len() {
        return Err(
            FiniteSolidSphereSceneError::Invalid("scene/native body count mismatch").into(),
        );
    }
    Ok(())
}
