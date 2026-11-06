//! Native thermal midpoint motion with sustained static friction after fracture.
//!
//! Both temperature-held and explicitly resolved thermal transitions use the
//! same physical contact law and complete original midpoint response as the
//! mechanical owner. Thermal parameter/eigenstrain work remains explicit.

use fs_exec::Cx;
use fs_material::brittle::source::expansion::ResolvedThermalBrittle;

use crate::brittle::solid::checkpoint;
use crate::brittle::solid::dynamics::fracture::contact::friction::BrittleFaceFriction;
use crate::brittle::solid::dynamics::fracture::contact::friction::source::BrittleFaceFrictionSource;
pub use crate::brittle::solid::dynamics::fracture::contact::friction::sustained::source::BrittleFaceSourcedSustainedFrictionError;
use crate::brittle::solid::dynamics::fracture::contact::friction::sustained::source::{
    BrittleFaceSourcedSustainedFrictionResult, solve_released_face_sustained_sourced_frictional,
};
pub use crate::brittle::solid::dynamics::fracture::contact::friction::sustained::{
    BrittleFaceSustainedFrictionControls, BrittleFaceSustainedFrictionError,
};
use crate::brittle::solid::dynamics::fracture::contact::friction::sustained::{
    BrittleFaceSustainedFrictionResult, solve_released_face_sustained_frictional,
};
use crate::brittle::solid::dynamics::fracture::contact::sustained::BrittleFaceSustainedError;
use crate::brittle::solid::expansion::dynamics::ThermalBrittleTetMotion;
use crate::brittle::solid::expansion::dynamics::midpoint::ThermalBrittleTetMidpointEnergy;
use crate::plastic3::IncrementSettings;

pub type ThermalBrittleFaceSustainedFrictionStep<'a> = BrittleFaceSustainedFrictionResult<
    ThermalBrittleTetMotion<'a>,
    ThermalBrittleTetMidpointEnergy,
>;

pub type ThermalBrittleFaceSourcedSustainedFrictionStep<'a> =
    BrittleFaceSourcedSustainedFrictionResult<
        ThermalBrittleTetMotion<'a>,
        ThermalBrittleTetMidpointEnergy,
    >;

impl<'a> ThermalBrittleTetMotion<'a> {
    /// Carry normal and shear load through matching, sticking fracture faces at
    /// the current thermal source state. Forces are physical interval means.
    /// The shared native force solve retains all material, mass, work, geometry
    /// and cancellation checks. Sliding and unlocalized arrivals refuse.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_released_face_sustained_frictional(
        &self,
        dt_s: f64,
        friction: &[BrittleFaceFriction],
        prescribed_m: &[(usize, f64)],
        mean_nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: BrittleFaceSustainedFrictionControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleFaceSustainedFrictionStep<'a>, BrittleFaceSustainedFrictionError>
    {
        self.advance_released_face_sustained_frictional_at_materials(
            self.solid.materials,
            dt_s,
            friction,
            prescribed_m,
            mean_nodal_forces_n,
            settings,
            controls,
            cx,
        )
    }

    /// Couple an explicit thermal source transition to normal and static shear
    /// forces. Every trial replays the same accepted mechanical history and
    /// source successor through the original thermal midpoint. Total strain,
    /// eigenstrain and mechanical strain retain their native meanings.
    ///
    /// Parameter/eigenstrain work, elastic storage, fracture loss, kinetic
    /// storage, external work and support work remain separate. Original
    /// reference mass/density remains fixed when thermal source density changes.
    /// Static constraint work is independently bounded; it supplies no heating.
    ///
    /// This retains the original coefficient-first straight thermal/strain
    /// increment and the declared fixed contact ellipse. Resolve source knots,
    /// reversals and contact-law parameter changes with explicit subintervals.
    /// Any refusal leaves the original thermal body available for retry.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_released_face_sustained_frictional_at_materials<'next>(
        &self,
        next_materials: &'next [ResolvedThermalBrittle],
        dt_s: f64,
        friction: &[BrittleFaceFriction],
        prescribed_m: &[(usize, f64)],
        mean_nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: BrittleFaceSustainedFrictionControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleFaceSustainedFrictionStep<'next>, BrittleFaceSustainedFrictionError>
    where
        'a: 'next,
    {
        checkpoint(cx).map_err(BrittleFaceSustainedFrictionError::from)?;
        if friction.len() != self.solid.released_faces.len() {
            return Err(BrittleFaceSustainedFrictionError::Sustained(
                BrittleFaceSustainedError::Invalid(
                    "complete retained-face friction inputs required",
                ),
            ));
        }
        let directions: Vec<_> = friction
            .iter()
            .map(|source| source.tangent_direction)
            .collect();
        solve_released_face_sustained_frictional(
            self.released_face_motion_view(),
            &directions,
            prescribed_m,
            mean_nodal_forces_n,
            controls,
            cx,
            |face, row, _, _, _| {
                Ok::<_, BrittleFaceSustainedFrictionError>(friction[face].vertices[row.vertex])
            },
            |total| {
                self.midpoint_candidate(
                    next_materials,
                    dt_s,
                    prescribed_m,
                    total,
                    settings,
                    controls.contact.work,
                    cx,
                )
            },
        )
    }

    /// Resolve original interface sources for a temperature-held sticking
    /// interval, returning source candidates only with the accepted native body.
    /// Static roles and actual ordered material identities are required.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_released_face_sustained_sourced_frictional(
        &self,
        dt_s: f64,
        sources: &[BrittleFaceFrictionSource<'_>],
        prescribed_m: &[(usize, f64)],
        mean_nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: BrittleFaceSustainedFrictionControls,
        cx: &Cx<'_>,
    ) -> Result<
        ThermalBrittleFaceSourcedSustainedFrictionStep<'a>,
        BrittleFaceSourcedSustainedFrictionError,
    > {
        self.advance_released_face_sustained_sourced_frictional_at_materials(
            self.solid.materials,
            dt_s,
            sources,
            prescribed_m,
            mean_nodal_forces_n,
            settings,
            controls,
            cx,
        )
    }

    /// Couple explicit thermal motion to a source-held static contact ellipse.
    /// Contact conditions are caller supplied and held through this interval;
    /// the contact source refresh does not advance a finite-step thermal store,
    /// travel or native slip history. Bulk thermal evolution follows the actual
    /// material successor with original parameter and eigenstrain work.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_released_face_sustained_sourced_frictional_at_materials<'next>(
        &self,
        next_materials: &'next [ResolvedThermalBrittle],
        dt_s: f64,
        sources: &[BrittleFaceFrictionSource<'_>],
        prescribed_m: &[(usize, f64)],
        mean_nodal_forces_n: &[(usize, f64)],
        settings: IncrementSettings,
        controls: BrittleFaceSustainedFrictionControls,
        cx: &Cx<'_>,
    ) -> Result<
        ThermalBrittleFaceSourcedSustainedFrictionStep<'next>,
        BrittleFaceSourcedSustainedFrictionError,
    >
    where
        'a: 'next,
    {
        solve_released_face_sustained_sourced_frictional(
            self.released_face_motion_view(),
            |element| {
                next_materials
                    .get(element)
                    .map(|source| source.band().material().material())
            },
            sources,
            prescribed_m,
            mean_nodal_forces_n,
            controls,
            cx,
            |total| {
                self.midpoint_candidate(
                    next_materials,
                    dt_s,
                    prescribed_m,
                    total,
                    settings,
                    controls.contact.work,
                    cx,
                )
            },
        )
    }
}
