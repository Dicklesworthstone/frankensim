//! Simultaneous physical impacts between released thermal-fracture faces.
//!
//! The native contact and consistent-mass owners are shared with isothermal
//! fracture. Thermal source bindings, eigenstrain reference, native mechanical
//! history, total displacement, stored energy and cumulative fracture loss are
//! retained exactly. A fixed-configuration impact changes only velocity and
//! kinetic storage; it does not advance temperature, physical time or the
//! physical-step allowance, and it does not infer impact heat.
//!
//! These are paired, fixed-reference, frictionless normal events. Collision
//! location, persistent contact forces and finite sliding remain separate
//! physical problems.

/// Bounded native midpoint prefixes at the admitted thermal state.
pub mod approach;

use std::sync::Arc;

use fs_contact::multicontact::MultiImpactResult;
use fs_exec::Cx;

use super::super::ThermalBrittleTetMotion;
use crate::brittle::solid::checkpoint;
pub use crate::brittle::solid::dynamics::fracture::contact::{
    BrittleFaceImpactControls, BrittleFaceImpactError, BrittleFaceImpactRow, SolidImpulseReport,
};
use crate::brittle::solid::dynamics::fracture::contact::{
    ReleasedFaceImpactBody, ReleasedFaceImpactResponse, released_face_impact,
};
use crate::brittle::solid::expansion::ThermalBrittleTetSolid;

/// Complete native velocity successor with all thermal history unchanged.
#[derive(Debug)]
pub struct ThermalBrittleFaceImpactStep<'a> {
    pub motion: ThermalBrittleTetMotion<'a>,
    pub contact: MultiImpactResult,
    pub impulse: SolidImpulseReport,
    pub rows: Vec<BrittleFaceImpactRow>,
    pub maximum_velocity_residual_mps: f64,
    /// Native kinetic change minus support work minus contact relative work.
    pub relative_work_error_j: f64,
    /// Native kinetic change minus support work plus contact dissipation.
    pub energy_balance_error_j: f64,
    /// Sum of absolute paired contact couples, before any cancellation.
    pub absolute_contact_couple_kg_m2_s: f64,
}

impl<'a> ThermalBrittleTetMotion<'a> {
    /// Apply the shared simultaneous normal event to retained thermal fragments.
    ///
    /// Restitution follows `solid().released_faces()` order. Active, non-tied
    /// paired vertices enter one coupled contact problem, including initially
    /// separating pairs. Penetration and active tangential separation refuse;
    /// this operation never repairs positions. The original native consistent
    /// mass supplies every mobility response and the combined impulse.
    ///
    /// Held velocity DOFs retain their current speeds, with their external
    /// reaction work exposed by `impulse`. Actual velocities, contact work,
    /// passivity and internal impulse moments must all pass the same gates as
    /// the isothermal owner before a thermal successor is returned. Refusal
    /// leaves the original accepted thermal motion available for retry.
    pub fn impact_released_faces(
        &self,
        restitution: &[f64],
        held_velocity_dofs: &[usize],
        origin_m: [f64; 3],
        controls: BrittleFaceImpactControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleFaceImpactStep<'a>, BrittleFaceImpactError> {
        let ReleasedFaceImpactResponse {
            velocity_mps,
            contact,
            impulse,
            rows,
            maximum_velocity_residual_mps,
            relative_work_error_j,
            energy_balance_error_j,
            absolute_contact_couple_kg_m2_s,
        } = released_face_impact(
            ReleasedFaceImpactBody {
                reference_nodes_m: self.solid.geometry.nodes_m(),
                displacement_m: &self.solid.solution.displacement_m,
                released_faces: &self.solid.released_faces,
                mass: &self.mass,
                velocity_mps: &self.velocity_mps,
                kinetic_energy_j: self.kinetic_energy_j,
            },
            restitution,
            held_velocity_dofs,
            origin_m,
            controls,
            cx,
        )?;
        let solid = ThermalBrittleTetSolid {
            geometry: self.solid.geometry,
            materials: self.solid.materials,
            bands: Arc::clone(&self.solid.bands),
            released_faces: Arc::clone(&self.solid.released_faces),
            budget: self.solid.budget,
            solution: self.solid.solution.clone(),
            stored_energy_j: self.solid.stored_energy_j,
            fracture_dissipation_j: self.solid.fracture_dissipation_j,
        };
        checkpoint(cx)?;
        Ok(ThermalBrittleFaceImpactStep {
            motion: ThermalBrittleTetMotion {
                solid,
                velocity_mps,
                reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
                mass: Arc::clone(&self.mass),
                kinetic_energy_j: impulse.kinetic_after_j,
                time_s: self.time_s,
                maximum_steps: self.maximum_steps,
                accepted_steps: self.accepted_steps,
            },
            contact,
            impulse,
            rows,
            maximum_velocity_residual_mps,
            relative_work_error_j,
            energy_balance_error_j,
            absolute_contact_couple_kg_m2_s,
        })
    }
}
