//! Zero-duration momentum transfer with the original mixed viscoelastic state.

use std::sync::Arc;

use fs_exec::Cx;
use fs_material::hyper::mixed::VolumetricTreatment;

use super::{
    MixedViscoTetError, MixedViscoTetMotion, MixedViscoTetSolid, SolidImpulseControls,
    SolidImpulseReport, checkpoint, inertia_error,
};

/// Complete candidate retaining continuous pressure and every Maxwell history.
#[derive(Debug)]
pub struct MixedViscoMotionImpulse<'a> {
    pub motion: MixedViscoTetMotion<'a>,
    /// Actual momentum, kinetic change, applied impulses and signed support work.
    pub impulse: SolidImpulseReport,
}

impl<'a> MixedViscoTetMotion<'a> {
    /// Apply instantaneous nodal impulses through the original consistent P2 mass.
    ///
    /// Only the finite-compliance `OriginalMaterial` formulation is admitted.
    /// Exact incompressibility needs a pressure-constrained impulse solve and
    /// returns `ExactIncompressibleImpulse`; no pressure mass is invented.
    ///
    /// Held velocity DOFs keep their accepted values, including moving supports.
    /// The native mass owner retains signed corner/edge entries, repeated-load
    /// accumulation, bounded CG and actual free-row/momentum/work gates. Angular
    /// momentum uses current deformed positions about the supplied origin.
    ///
    /// Only velocity changes. Pressure, deformation, all branch metrics and
    /// dissipated histories, stresses, storage, time and reference density are
    /// unchanged. Subsequent `advance` calls continue those same histories.
    /// The caller owns the contact/restitution law and opposite body's impulse;
    /// this event creates no material relaxation loss or heat.
    ///
    /// # Errors
    /// Exact incompressibility, invalid inputs, native mass-solve/work/momentum
    /// refusal, or cancellation. No accepted field changes on any refusal.
    pub fn impulse_candidate(
        &self,
        held_velocity_dofs: &[usize],
        nodal_impulses_kg_m_s: &[(usize, f64)],
        origin_m: [f64; 3],
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<MixedViscoMotionImpulse<'a>, MixedViscoTetError> {
        checkpoint(cx)?;
        if self.solid.volumetric == VolumetricTreatment::ExactIncompressible {
            return Err(MixedViscoTetError::ExactIncompressibleImpulse);
        }
        let mut positions_m = Vec::with_capacity(self.velocity_mps.len());
        for (node, (reference, displacement)) in self
            .solid
            .nodes_m
            .iter()
            .zip(&self.solid.displacement_m)
            .enumerate()
        {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            positions_m.push(core::array::from_fn(|axis| {
                reference[axis] + displacement[axis]
            }));
        }
        let (velocity_mps, impulse) = self
            .mass
            .impulse_candidate(
                &positions_m,
                &self.velocity_mps,
                held_velocity_dofs,
                nodal_impulses_kg_m_s,
                origin_m,
                controls,
                cx,
            )
            .map_err(inertia_error)?;
        let solid = MixedViscoTetSolid {
            nodes_m: self.solid.nodes_m,
            elements: self.solid.elements,
            material: self.solid.material,
            volumetric: self.solid.volumetric,
            mesh: Arc::clone(&self.solid.mesh),
            budget: self.solid.budget,
            quadrature_points_per_dir: self.solid.quadrature_points_per_dir,
            maximum_quadrature_points: self.solid.maximum_quadrature_points,
            displacement_m: self.solid.displacement_m.clone(),
            pressure_pa: self.solid.pressure_pa.clone(),
            states: self.solid.states.clone(),
            first_piola_pa: self.solid.first_piola_pa.clone(),
            remainder_energy_j: self.solid.remainder_energy_j,
            bulk_energy_j: self.solid.bulk_energy_j,
            time_s: self.solid.time_s,
        };
        checkpoint(cx)?;
        Ok(MixedViscoMotionImpulse {
            motion: Self {
                solid,
                velocity_mps,
                reference_density_kg_m3: Arc::clone(&self.reference_density_kg_m3),
                mass: Arc::clone(&self.mass),
            },
            impulse,
        })
    }
}
