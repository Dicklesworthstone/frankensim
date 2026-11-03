//! Native curved witnesses through the original complete-state impulse owners.

use super::{
    QuadraticSolidContactSurface, QuadraticSolidSurfaceControls, QuadraticSolidSurfaceError,
    QuadraticSolidSurfaceQuery,
};
use crate::finite_solid::contact::{FiniteSolidRigidContactError, surface::FiniteSolidSphere};
use crate::inelastic_impact::multiple::{
    InelasticMixedContactCandidate, InelasticMixedContactControls,
    mixed_visco_mixed_contact_candidate,
    volume::{
        VolumeContactControls, conservative_contact_candidate, mixed_visco_contact_candidate,
    },
};
use fs_contact::frictional_multicontact::{FrictionalContactError, FrictionalContactResult};
use fs_exec::Cx;
use fs_mbd::{MassProperties, RigidBodyState};
use fs_solid::hyper3::mixed::{
    conservative::{ConservativeTetError, ConservativeTetSolid},
    impulse::MixedImpulseReport,
    visco::{
        MixedViscoTetError,
        dynamics::{MixedViscoTetMotion, SolidImpulseReport},
    },
};

/// Current curved geometry and the native complete mechanical transaction.
#[derive(Debug)]
pub struct QuadraticSphereImpactCandidate<M, I = SolidImpulseReport> {
    pub mechanical: InelasticMixedContactCandidate<M, FrictionalContactResult, I>,
    pub geometry: QuadraticSolidSurfaceQuery,
}

/// Retain the original geometry, material, contact and physical-audit refusals.
#[derive(Debug)]
pub enum QuadraticSphereImpactError<S, I = SolidImpulseReport> {
    Geometry(QuadraticSolidSurfaceError),
    Contact(FiniteSolidRigidContactError<FrictionalContactResult, FrictionalContactError, S, I>),
}
impl<S: core::fmt::Debug, I: core::fmt::Debug> core::fmt::Display
    for QuadraticSphereImpactError<S, I>
{
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "quadratic sphere impact refused: {self:?}")
    }
}
impl<S: core::fmt::Debug, I: core::fmt::Debug> std::error::Error
    for QuadraticSphereImpactError<S, I>
{
}

// Both native owners expose the same immutable geometry. Query preparation and
// refits borrow their material histories; only a successful impulse stages one
// complete successor through the existing material/momentum owner.
macro_rules! native_geometry {
    ($prepare:ident, $query:ident, $solid:ident: $motion:ty, $geometry:expr) => {
        impl QuadraticSolidContactSurface {
            /// Bind the full original TET10 incidence and curved reference nodes.
            pub fn $prepare(
                $solid: &$motion,
                maximum_face_entries: usize,
                cx: &Cx<'_>,
            ) -> Result<Self, QuadraticSolidSurfaceError> {
                let geometry = $geometry;
                Self::from_reference(
                    geometry.reference_nodes_m(),
                    geometry.reference_tetrahedra(),
                    geometry.assembly_budget(),
                    maximum_face_entries,
                    cx,
                )
            }

            /// Refit actual current quadratic coordinates and retain the original
            /// material barycentrics. A raw query may report penetration; only
            /// its separate admission can authorize a mechanical candidate.
            pub fn $query(
                &self,
                $solid: &$motion,
                bodies: &[RigidBodyState],
                spheres: &[FiniteSolidSphere],
                controls: QuadraticSolidSurfaceControls,
                cx: &Cx<'_>,
            ) -> Result<QuadraticSolidSurfaceQuery, QuadraticSolidSurfaceError> {
                let geometry = $geometry;
                let positions = self.current_positions_from_reference(
                    geometry.reference_nodes_m(),
                    geometry.reference_tetrahedra(),
                    geometry.displacement_m(),
                    cx,
                )?;
                self.sphere_separations_at_positions(&positions, bodies, spheres, controls, cx)
            }
        }
    };
}

native_geometry!(for_mixed_visco, mixed_visco_sphere_separations,
    solid: MixedViscoTetMotion<'_>, solid.solid());
native_geometry!(for_conservative, conservative_sphere_separations,
    solid: ConservativeTetSolid<'_>, solid);

impl QuadraticSolidContactSurface {
    /// Detect actual curved-boundary sphere contacts and resolve their complete
    /// finite-compliance Maxwell response in one original simultaneous solve.
    /// No witness is moved to hide its gap, torque, or geometric uncertainty.
    #[allow(clippy::too_many_arguments)]
    pub fn mixed_visco_sphere_contact_candidate<'a>(
        &self,
        solid: &MixedViscoTetMotion<'a>,
        held: &[usize],
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        spheres: &[FiniteSolidSphere],
        origin_m: [f64; 3],
        geometry_controls: QuadraticSolidSurfaceControls,
        impact_controls: InelasticMixedContactControls,
        cx: &Cx<'_>,
    ) -> Result<
        QuadraticSphereImpactCandidate<MixedViscoTetMotion<'a>>,
        QuadraticSphereImpactError<MixedViscoTetError>,
    > {
        let geometry = self
            .mixed_visco_sphere_separations(solid, bodies, spheres, geometry_controls, cx)
            .map_err(QuadraticSphereImpactError::Geometry)?;
        let contacts = geometry
            .admitted_contacts()
            .map_err(QuadraticSphereImpactError::Geometry)?;
        let mechanical = mixed_visco_mixed_contact_candidate(
            solid,
            held,
            bodies,
            properties,
            &contacts,
            origin_m,
            impact_controls,
            cx,
        )
        .map_err(QuadraticSphereImpactError::Contact)?;
        Ok(QuadraticSphereImpactCandidate {
            mechanical,
            geometry,
        })
    }

    /// Same admitted curved witnesses with the native exact-pressure mobility.
    /// One factor supplies all contact directions and the final impulse; pressure
    /// impulses remain separate from retained physical pressure and history.
    #[allow(clippy::too_many_arguments)]
    pub fn mixed_visco_volume_sphere_contact_candidate<'a>(
        &self,
        solid: &MixedViscoTetMotion<'a>,
        held: &[usize],
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        spheres: &[FiniteSolidSphere],
        pressure_gauge_node: Option<usize>,
        origin_m: [f64; 3],
        geometry_controls: QuadraticSolidSurfaceControls,
        impact_controls: VolumeContactControls,
        cx: &Cx<'_>,
    ) -> Result<
        QuadraticSphereImpactCandidate<MixedViscoTetMotion<'a>, MixedImpulseReport>,
        QuadraticSphereImpactError<MixedViscoTetError, MixedImpulseReport>,
    > {
        let geometry = self
            .mixed_visco_sphere_separations(solid, bodies, spheres, geometry_controls, cx)
            .map_err(QuadraticSphereImpactError::Geometry)?;
        let contacts = geometry
            .admitted_contacts()
            .map_err(QuadraticSphereImpactError::Geometry)?;
        let mechanical = mixed_visco_contact_candidate(
            solid,
            held,
            bodies,
            properties,
            &contacts,
            pressure_gauge_node,
            origin_m,
            impact_controls,
            cx,
        )
        .map_err(QuadraticSphereImpactError::Contact)?;
        Ok(QuadraticSphereImpactCandidate {
            mechanical,
            geometry,
        })
    }

    /// Curved-boundary detection with the original conservative solid response.
    /// Native formulation dispatch, cumulative work and retained material history
    /// remain with the same conservative owner used by authored point impulses.
    #[allow(clippy::too_many_arguments)]
    pub fn conservative_sphere_contact_candidate<'a>(
        &self,
        solid: &ConservativeTetSolid<'a>,
        held: &[usize],
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        spheres: &[FiniteSolidSphere],
        pressure_gauge_node: Option<usize>,
        origin_m: [f64; 3],
        geometry_controls: QuadraticSolidSurfaceControls,
        impact_controls: VolumeContactControls,
        cx: &Cx<'_>,
    ) -> Result<
        QuadraticSphereImpactCandidate<ConservativeTetSolid<'a>, MixedImpulseReport>,
        QuadraticSphereImpactError<ConservativeTetError, MixedImpulseReport>,
    > {
        let geometry = self
            .conservative_sphere_separations(solid, bodies, spheres, geometry_controls, cx)
            .map_err(QuadraticSphereImpactError::Geometry)?;
        let contacts = geometry
            .admitted_contacts()
            .map_err(QuadraticSphereImpactError::Geometry)?;
        let mechanical = conservative_contact_candidate(
            solid,
            held,
            bodies,
            properties,
            &contacts,
            pressure_gauge_node,
            origin_m,
            impact_controls,
            cx,
        )
        .map_err(QuadraticSphereImpactError::Contact)?;
        Ok(QuadraticSphereImpactCandidate {
            mechanical,
            geometry,
        })
    }
}
