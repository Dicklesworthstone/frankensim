//! Coupled directional-friction impulses on native thermal brittle fragments.
//!
//! The mechanical and thermal owners share the same sparse contact assembly,
//! consistent mass, simultaneous elliptic law, actual endpoint and work audits.
//! This adapter retains the complete accepted thermal/constitutive state and
//! reference matter while publishing only the admitted velocity successor.

/// Native thermal midpoint motion with sustained static friction forces.
pub mod sustained;

/// Complete physical durations with directional impacts, static contact and release.
pub mod time;

use std::sync::Arc;

use fs_contact::frictional_multicontact::elliptic::{
    EllipticContactResiduals, EllipticFrictionalContactResult,
};
use fs_exec::Cx;

use crate::brittle::solid::checkpoint;
use crate::brittle::solid::dynamics::fracture::contact::ReleasedFaceImpactBody;
use crate::brittle::solid::dynamics::fracture::contact::friction::source::released_face_sourced_frictional_impact;
pub use crate::brittle::solid::dynamics::fracture::contact::friction::source::{
    BrittleFaceFrictionSample, BrittleFaceFrictionSource, BrittleFaceSourcedFrictionError,
    BrittleVertexFrictionSource,
};
pub use crate::brittle::solid::dynamics::fracture::contact::friction::{
    BrittleFaceFriction, BrittleFaceFrictionControls, BrittleFaceFrictionError,
    BrittleFaceFrictionFrame, EllipticFriction,
};
use crate::brittle::solid::dynamics::fracture::contact::friction::{
    ReleasedFaceFrictionResponse, released_face_frictional_impact_with,
};
pub use crate::brittle::solid::dynamics::fracture::contact::{
    BrittleFaceImpactRow, SolidImpulseReport,
};
use crate::brittle::solid::expansion::ThermalBrittleTetSolid;
use crate::brittle::solid::expansion::dynamics::ThermalBrittleTetMotion;

/// Accepted thermal velocity successor, with the shared original contact law.
#[derive(Debug)]
pub struct ThermalBrittleFaceFrictionStep<'a> {
    pub motion: ThermalBrittleTetMotion<'a>,
    pub contact: EllipticFrictionalContactResult,
    pub impulse: SolidImpulseReport,
    pub rows: Vec<BrittleFaceImpactRow>,
    pub frames: Vec<BrittleFaceFrictionFrame>,
    pub measured_residuals: EllipticContactResiduals,
    pub maximum_velocity_residual_mps: f64,
    /// Actual total midpoint contact work; no thermal destination is inferred.
    pub relative_work_j: f64,
    pub contact_work_error_j: f64,
    pub relative_work_error_j: f64,
    pub energy_balance_error_j: f64,
    pub absolute_contact_couple_kg_m2_s: f64,
}

/// One accepted thermal event and exactly its active source candidates.
#[derive(Debug)]
pub struct ThermalBrittleFaceSourcedFrictionStep<'a> {
    pub impact: ThermalBrittleFaceFrictionStep<'a>,
    pub samples: Vec<BrittleFaceFrictionSample>,
}

impl<'a> ThermalBrittleTetMotion<'a> {
    /// Apply the shared simultaneous directional-friction law at the current
    /// thermal configuration. Inputs follow retained face and vertex order.
    ///
    /// The full native consistent mass, moving support reactions, actual
    /// endpoint ellipse and total impact work retain the same gates as the
    /// mechanical owner. Coefficients are explicit model inputs. Active paired
    /// points must satisfy the same strict opening and tangential geometry;
    /// this instantaneous event does not supply finite-sliding geometry.
    ///
    /// Source selections, temperatures, eigenstrain reference, mechanical
    /// strain/history, stored energy, cumulative fracture loss, reference
    /// density, positions, time and step allowance remain exactly retained.
    /// Total contact loss does not become heat or a second fracture debit.
    /// Any refusal leaves the complete accepted thermal motion available.
    pub fn impact_released_faces_frictional(
        &self,
        restitution: &[f64],
        friction: &[BrittleFaceFriction],
        held_velocity_dofs: &[usize],
        origin_m: [f64; 3],
        controls: BrittleFaceFrictionControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleFaceFrictionStep<'a>, BrittleFaceFrictionError> {
        checkpoint(cx)?;
        if friction.len() != self.solid.released_faces.len() {
            return Err(BrittleFaceFrictionError::Invalid(
                "complete retained-face friction inputs required",
            ));
        }
        let directions: Vec<_> = friction.iter().map(|f| f.tangent_direction).collect();
        let response = released_face_frictional_impact_with::<BrittleFaceFrictionError>(
            self.frictional_impact_body(),
            restitution,
            &directions,
            held_velocity_dofs,
            origin_m,
            controls,
            cx,
            |face, row, _, _| Ok(friction[face].vertices[row.vertex]),
        )?;
        ThermalBrittleFaceFrictionStep::from_response(self, response, cx)
    }

    /// Sample each active original source at its actual incoming slip and
    /// return the admitted thermal body/source candidates together.
    ///
    /// Both ordered interface material IDs must match the actual thermal
    /// element bindings. The caller supplies complete source conditions and
    /// an explicit coefficient role; the original source owner retains its
    /// domains, selected claims, transported texture and accepted history.
    /// Inactive/tied vertices are not sampled. No source, thermal or mechanical
    /// input is mutated, including after a late contact or native work refusal.
    pub fn impact_released_faces_sourced_frictional(
        &self,
        restitution: &[f64],
        sources: &[BrittleFaceFrictionSource<'_>],
        held_velocity_dofs: &[usize],
        origin_m: [f64; 3],
        controls: BrittleFaceFrictionControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleFaceSourcedFrictionStep<'a>, BrittleFaceSourcedFrictionError> {
        let (response, samples) = released_face_sourced_frictional_impact(
            self.frictional_impact_body(),
            |element| {
                self.solid
                    .materials
                    .get(element)
                    .map(|source| source.band().material().material())
            },
            restitution,
            sources,
            held_velocity_dofs,
            origin_m,
            controls,
            cx,
        )?;
        let impact = ThermalBrittleFaceFrictionStep::from_response(self, response, cx)?;
        Ok(ThermalBrittleFaceSourcedFrictionStep { impact, samples })
    }

    fn frictional_impact_body(&self) -> ReleasedFaceImpactBody<'_> {
        ReleasedFaceImpactBody {
            reference_nodes_m: self.solid.geometry.nodes_m(),
            displacement_m: &self.solid.solution.displacement_m,
            released_faces: &self.solid.released_faces,
            mass: &self.mass,
            velocity_mps: &self.velocity_mps,
            kinetic_energy_j: self.kinetic_energy_j,
        }
    }

    /// Preserve the complete native thermal owner while accepting only the
    /// already audited contact velocity and kinetic-energy successor.
    pub(crate) fn friction_velocity_successor(
        &self,
        velocity_mps: Vec<[f64; 3]>,
        kinetic_energy_j: f64,
        cx: &Cx<'_>,
    ) -> Result<Self, BrittleFaceFrictionError> {
        let mut solid = ThermalBrittleTetSolid {
            geometry: self.solid.geometry,
            materials: self.solid.materials,
            bands: Arc::clone(&self.solid.bands),
            released_faces: Arc::clone(&self.solid.released_faces),
            budget: self.solid.budget,
            solution: self.solid.solution.clone(),
            stored_energy_j: self.solid.stored_energy_j,
            fracture_dissipation_j: self.solid.fracture_dissipation_j,
        };
        solid
            .refresh_rate_stresses(&velocity_mps, &self.free_linear_strain_rates_per_s, cx)
            .map_err(BrittleFaceFrictionError::Motion)?;
        checkpoint(cx)?;
        Ok(Self {
            solid,
            velocity_mps,
            free_linear_strain_rates_per_s: Arc::clone(&self.free_linear_strain_rates_per_s),
            reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
            mass: Arc::clone(&self.mass),
            kinetic_energy_j,
            time_s: self.time_s,
            maximum_steps: self.maximum_steps,
            accepted_steps: self.accepted_steps,
        })
    }
}

impl<'a> ThermalBrittleFaceFrictionStep<'a> {
    fn from_response(
        motion: &ThermalBrittleTetMotion<'a>,
        response: ReleasedFaceFrictionResponse,
        cx: &Cx<'_>,
    ) -> Result<Self, BrittleFaceFrictionError> {
        Ok(Self {
            motion: motion.friction_velocity_successor(
                response.velocity_mps,
                response.impulse.kinetic_after_j,
                cx,
            )?,
            contact: response.contact,
            impulse: response.impulse,
            rows: response.rows,
            frames: response.frames,
            measured_residuals: response.measured_residuals,
            maximum_velocity_residual_mps: response.maximum_velocity_residual_mps,
            relative_work_j: response.relative_work_j,
            contact_work_error_j: response.contact_work_error_j,
            relative_work_error_j: response.relative_work_error_j,
            energy_balance_error_j: response.energy_balance_error_j,
            absolute_contact_couple_kg_m2_s: response.absolute_contact_couple_kg_m2_s,
        })
    }
}
