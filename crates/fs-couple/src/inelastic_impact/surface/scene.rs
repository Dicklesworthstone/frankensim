//! Complete canonical sphere scenes on original native Tet4 material bodies.
//!
//! Surface queries borrow the accepted plastic, Maxwell or thermal geometry.
//! The original scene composes those witnesses with its own rigid-pair broad
//! phase, separation tests, explicit pair laws and disjoint contact-ID scheme.
//! Only solid point coordinates are expressed as the native Tet4 descriptor.
//! No geometry, material history, source identity or normal law is rebound.

use super::{FiniteSolidContactSurface, FiniteSolidSurfaceControls};
use fs_exec::Cx;
use fs_material::plastic::finite::{FiniteJ2Material, thermal::FiniteThermalJ2Material};
use fs_mbd::RigidBodyState;
use fs_solid::hyper3::{
    plastic::{
        dynamics::FinitePlasticTetMotion,
        thermal::{FiniteThermoplasticTetSolid, TetImpulsePoint},
    },
    visco::dynamics::FiniteViscoTetMotion,
};

/// Same prepared topology, explicit pair laws, original keys and refusals.
pub use crate::finite_solid::contact::scene::{
    FiniteSolidSphereScene, FiniteSolidSphereSceneError, FiniteSolidSphereSceneKey,
    FiniteSolidSphereSceneQuery, SpherePairLaw,
};
/// Complete query in the original scene-ID order, ready for one native mixed
/// solve. It also remains valid across a force kick at unchanged coordinates.
pub type InelasticSphereSceneQuery = FiniteSolidSphereSceneQuery<TetImpulsePoint>;

macro_rules! native_scene {
    ($method:ident, [$($generics:tt)*], $motion:ty, $surface_query:ident) => {
        impl FiniteSolidSphereScene {
            /// Detect and admit all current sphere/solid and sphere/sphere rows.
            /// Original exact-reference checks, penetration refusals, pair laws
            /// and contact/work budgets remain with the existing owners.
            pub fn $method<$($generics)*>(
                &self,
                surface: &FiniteSolidContactSurface,
                solid: &$motion,
                bodies: &[RigidBodyState],
                controls: FiniteSolidSurfaceControls,
                cx: &Cx<'_>,
            ) -> Result<InelasticSphereSceneQuery, FiniteSolidSphereSceneError> {
                cx.checkpoint().map_err(|_| FiniteSolidSphereSceneError::Cancelled)?;
                if bodies.len() != self.spheres().len() {
                    return Err(FiniteSolidSphereSceneError::Invalid(
                        "scene/native body count mismatch",
                    ));
                }
                let query = surface.$surface_query(solid, bodies, self.spheres(), controls, cx)?;
                self.contacts_from_surface_query(query, bodies, controls, cx, |point| {
                    TetImpulsePoint {
                        tetrahedron: point.tetrahedron(),
                        weights: point.weights(),
                    }
                })
            }
        }
    };
}

native_scene!(plastic_contacts, ['a, L: FiniteJ2Material],
    FinitePlasticTetMotion<'a, L>, plastic_sphere_separations);
native_scene!(visco_contacts, ['a], FiniteViscoTetMotion<'a>, visco_sphere_separations);
native_scene!(thermoplastic_tet_contacts, ['a, L: FiniteThermalJ2Material],
    FiniteThermoplasticTetSolid<'a, L>, thermoplastic_tet_sphere_separations);
