//! Current Tet4 exterior detection feeding retained-history native impacts.
//!
//! Preparation, exact reference matching, BVH refits, witness deduplication and
//! geometry refusals remain with the original finite-solid surface owner. Only
//! the accepted body's geometry and material-point payload differ here. Queries
//! borrow all plastic/Maxwell/thermal history; the original simultaneous contact
//! owner stages a complete candidate after detection and momentum/work audits.
//! This discrete event supplies no swept-path claim, penetration repair or heat
//! partition. Hex8 and curved P2 boundaries require their own exact geometry.

/// Curved quadratic material boundaries through native momentum response.
pub mod quadratic;

/// Complete canonical rigid/solid and rigid/rigid sphere scenes on native bodies.
pub mod scene;

use super::multiple::{
    FiniteSolidContactCounterpart, InelasticMixedContact, InelasticMixedContactCandidate,
    InelasticMixedContactControls, InelasticMixedContactError, RigidImpactPoint,
    plastic_mixed_contact_candidate, thermoplastic_tet_mixed_contact_candidate,
    visco_mixed_contact_candidate,
};
pub use crate::finite_solid::contact::surface::{
    FiniteSolidContactSurface, FiniteSolidSphere, FiniteSolidSphereWitness,
    FiniteSolidSurfaceControls, FiniteSolidSurfaceError, FiniteSolidSurfaceQuery,
};
use fs_exec::Cx;
use fs_material::plastic::finite::{FiniteJ2Material, thermal::FiniteThermalJ2Material};
use fs_mbd::{MassProperties, RigidBodyState};
use fs_solid::hyper3::{
    plastic::{
        FinitePlasticTetError,
        dynamics::FinitePlasticTetMotion,
        thermal::{FiniteThermoplasticTetError, FiniteThermoplasticTetSolid, TetImpulsePoint},
    },
    visco::{FiniteViscoTetError, dynamics::FiniteViscoTetMotion},
};

/// Original geometry receipt and complete simultaneous mechanical candidate.
#[derive(Debug)]
pub struct InelasticSphereContactCandidate<M> {
    pub mechanical: InelasticMixedContactCandidate<M>,
    pub geometry: FiniteSolidSurfaceQuery,
}

/// Preserve original geometric admission and native contact/source failures.
#[derive(Debug)]
pub enum InelasticSphereContactError<E> {
    Geometry(FiniteSolidSurfaceError),
    Contact(InelasticMixedContactError<E>),
}
impl<E: core::fmt::Debug> core::fmt::Display for InelasticSphereContactError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native inelastic sphere contact refused: {self:?}")
    }
}
impl<E: core::fmt::Debug> std::error::Error for InelasticSphereContactError<E> {}

impl FiniteSolidSurfaceQuery {
    /// Admit the original actual witnesses and express their unchanged material
    /// coordinates in the native Tet4 payload. Stable contact IDs are retained.
    pub fn admitted_inelastic_contacts(
        &self,
    ) -> Result<Vec<InelasticMixedContact>, FiniteSolidSurfaceError> {
        Ok(self
            .admitted_contacts()?
            .into_iter()
            .map(|contact| InelasticMixedContact {
                id: contact.id,
                point: RigidImpactPoint {
                    body: contact.rigid_body,
                    arm_body: contact.arm_body_m,
                },
                counterpart: FiniteSolidContactCounterpart::Solid(TetImpulsePoint {
                    tetrahedron: contact.solid_point.tetrahedron(),
                    weights: contact.solid_point.weights(),
                }),
                normal_world: contact.normal_world,
                tangent_world: contact.tangent_world,
                activation_tolerance_m: contact.activation_tolerance_m,
                friction_coefficient: contact.friction_coefficient,
                normal_law: contact.normal_law,
            })
            .collect())
    }
}

// Only native geometry access and the original candidate type vary. Detection
// and impulse assembly stay with their existing owners, without history clones
// during preparation or queries and without an elastic proxy body.
macro_rules! native_surface {
    ($prepare:ident, $query:ident, $candidate:ident,
     [$($generics:tt)*], $solid:ident: $motion:ty, $geometry:expr,
     $error:ty, $impact:ident) => {
        impl FiniteSolidContactSurface {
            /// Prepare the original closed reference exterior and cached BVH.
            pub fn $prepare<$($generics)*>(
                $solid: &$motion,
                maximum_face_entries: usize,
                cx: &Cx<'_>,
            ) -> Result<Self, FiniteSolidSurfaceError> {
                let geometry = $geometry;
                Self::from_reference(
                    geometry.reference_nodes_m(),
                    geometry.reference_tetrahedra(),
                    geometry.assembly_budget(),
                    maximum_face_entries,
                    cx,
                )
            }

            /// Refit CURRENT native coordinates while preserving exact reference
            /// identity and original material incidence. Raw gaps may be negative.
            pub fn $query<$($generics)*>(
                &self,
                $solid: &$motion,
                rigid_bodies: &[RigidBodyState],
                spheres: &[FiniteSolidSphere],
                controls: FiniteSolidSurfaceControls,
                cx: &Cx<'_>,
            ) -> Result<FiniteSolidSurfaceQuery, FiniteSolidSurfaceError> {
                let geometry = $geometry;
                let positions = self.current_positions_from_reference(
                    geometry.reference_nodes_m(),
                    geometry.reference_tetrahedra(),
                    geometry.displacement_m(),
                    cx,
                )?;
                self.sphere_separations_at_positions(&positions, rigid_bodies, spheres, controls, cx)
            }

            /// Detect and admit all sphere/exterior contacts, then resolve them
            /// together with the original mass, material histories and supports.
            /// Sphere/sphere detection and time integration remain separate owners.
            #[allow(clippy::too_many_arguments)]
            pub fn $candidate<$($generics)*>(
                &self,
                $solid: &$motion,
                held: &[usize],
                rigid_bodies: &[RigidBodyState],
                properties: &[MassProperties],
                spheres: &[FiniteSolidSphere],
                origin: [f64; 3],
                geometry_controls: FiniteSolidSurfaceControls,
                impact_controls: InelasticMixedContactControls,
                cx: &Cx<'_>,
            ) -> Result<InelasticSphereContactCandidate<$motion>, InelasticSphereContactError<$error>> {
                let geometry = self.$query($solid, rigid_bodies, spheres, geometry_controls, cx)
                    .map_err(InelasticSphereContactError::Geometry)?;
                let contacts = geometry.admitted_inelastic_contacts()
                    .map_err(InelasticSphereContactError::Geometry)?;
                let mechanical = $impact(
                    $solid, held, rigid_bodies, properties, &contacts, origin, impact_controls, cx,
                ).map_err(InelasticSphereContactError::Contact)?;
                Ok(InelasticSphereContactCandidate { mechanical, geometry })
            }
        }
    };
}

native_surface!(
    for_plastic, plastic_sphere_separations, plastic_sphere_contact_candidate,
    ['a, L: FiniteJ2Material], solid: FinitePlasticTetMotion<'a, L>, solid.solid(),
    FinitePlasticTetError, plastic_mixed_contact_candidate
);
native_surface!(
    for_visco, visco_sphere_separations, visco_sphere_contact_candidate,
    ['a], solid: FiniteViscoTetMotion<'a>, solid.solid(),
    FiniteViscoTetError, visco_mixed_contact_candidate
);
native_surface!(
    for_thermoplastic_tet, thermoplastic_tet_sphere_separations,
    thermoplastic_tet_sphere_contact_candidate,
    ['a, L: FiniteThermalJ2Material], solid: FiniteThermoplasticTetSolid<'a, L>, solid,
    FiniteThermoplasticTetError, thermoplastic_tet_mixed_contact_candidate
);
