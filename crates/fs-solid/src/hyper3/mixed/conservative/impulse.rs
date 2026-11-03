//! Same-position pressure impulses preserve the complete conservative history.

use super::{
    ConservativeFormulation, ConservativeTetError, ConservativeTetSolid, MixedUpError, checkpoint,
    gate, invalid,
};
use crate::hyper3::mixed::impulse::{
    MixedImpulseControls, MixedImpulseReport, constrained_impulse,
};
use fs_exec::Cx;
use std::sync::Arc;

#[derive(Debug)]
pub struct ConservativeTetImpulse<'a> {
    pub solid: ConservativeTetSolid<'a>,
    pub impulse: MixedImpulseReport,
}

impl<'a> ConservativeTetSolid<'a> {
    /// Apply an instantaneous load through original mass AND exact weak volume
    /// constraints at the accepted configuration. Geometry, material points,
    /// reference density, stored pressure and time are retained unchanged.
    ///
    /// Held components keep their current velocity. A pressure gauge is zero
    /// impulse on one original corner and is admitted only for a true connected
    /// constant null mode. Other rank deficiencies refuse in the native solver.
    ///
    /// Pre-event velocity must satisfy the stated weak volume-rate allowance.
    /// The time integrator enforces endpoint volume, so its displacement-chord
    /// velocity is not automatically an admissible tangent velocity. This event
    /// does not conceal a large pre-existing defect by a dissipative projection.
    ///
    /// # Errors
    /// Non-exact formulation, invalid input, incompatible velocity/gauge,
    /// insufficient sparse/workspace budgets, native factorization or physical
    /// balance failure, or cancellation. No accepted state is mutated.
    #[allow(clippy::too_many_arguments)]
    pub fn impulse_candidate(
        &self,
        held_velocity_dofs: &[usize],
        nodal_impulses_kg_m_s: &[(usize, f64)],
        pressure_gauge_node: Option<usize>,
        origin_m: [f64; 3],
        controls: MixedImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<ConservativeTetImpulse<'a>, ConservativeTetError> {
        checkpoint(cx)?;
        let ConservativeFormulation::ExactIncompressible {
            maximum_point_volume_work_defect,
            ..
        } = self.definition.formulation
        else {
            return Err(invalid(
                "pressure impulse requires exact-incompressible conservative motion",
            ));
        };
        let positions: Vec<[f64; 3]> = self
            .definition
            .nodes_m
            .iter()
            .zip(&self.displacement_m)
            .map(|(x, u)| core::array::from_fn(|c| x[c] + u[c]))
            .collect();
        let (velocity_mps, impulse) = constrained_impulse(
            &self.mesh,
            &self.mass,
            &positions,
            &self.velocity_mps,
            held_velocity_dofs,
            nodal_impulses_kg_m_s,
            pressure_gauge_node,
            origin_m,
            controls,
            |e, q| {
                let point = &self.points[e * self.mesh.points_per_element + q];
                point
                    .volume_state(
                        point.accepted_deformation(),
                        maximum_point_volume_work_defect,
                    )
                    .map(|volume| volume.endpoint_cofactor)
                    .map_err(|error| MixedUpError::Material {
                        element: e,
                        point: q,
                        error,
                    })
            },
            cx,
        )?;
        // Only external/support work enters the conservative trajectory ledger.
        // Pressure work is independently zero-gated; it is not supplied energy.
        let work = impulse.external_work_j + impulse.momentum.support_work_j;
        let increment = work - self.work_correction_j;
        let cumulative_work_j = self.cumulative_work_j + increment;
        let work_correction_j = (cumulative_work_j - self.cumulative_work_j) - increment;
        gate(
            "impulsive accumulated total work",
            cumulative_work_j
                - (self.stored_energy_j + impulse.momentum.kinetic_after_j
                    - self.initial_total_energy_j),
            controls.balance.work_tolerance_j,
        )?;
        let solid = Self {
            definition: self.definition,
            mesh: Arc::clone(&self.mesh),
            mass: Arc::clone(&self.mass),
            corner_reference_density: self.corner_reference_density.clone(),
            points: self.points.clone(),
            bulk_points: self.bulk_points.clone(),
            bulk_volume_ratios: self.bulk_volume_ratios.clone(),
            displacement_m: self.displacement_m.clone(),
            velocity_mps,
            algorithmic_pressure_pa: self.algorithmic_pressure_pa.clone(),
            pressure_residual_m3: self.pressure_residual_m3.clone(),
            stored_energy_j: self.stored_energy_j,
            bulk_stored_energy_j: self.bulk_stored_energy_j,
            kinetic_energy_j: impulse.momentum.kinetic_after_j,
            initial_total_energy_j: self.initial_total_energy_j,
            cumulative_work_j,
            work_correction_j,
            time_s: self.time_s,
        };
        checkpoint(cx)?;
        Ok(ConservativeTetImpulse { solid, impulse })
    }
}
