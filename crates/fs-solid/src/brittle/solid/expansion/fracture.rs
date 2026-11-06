//! Completed tensile cracks retain their thermal source and mechanical history.
//!
//! Connectivity release uses the original native admission and remapping
//! algorithm. The total FE strain and mechanical constitutive strain stay
//! distinct: no held-temperature state is replayed at total displacement strain.
//! Expansion references, coefficients, damage, stress, geometry and accumulated
//! fracture loss survive without another surface energy or material rebind.

use fs_exec::Cx;

use super::ThermalBrittleTetSolid;
use crate::brittle::solid::BrittleTetError;
use crate::brittle::solid::fracture::{
    BrittleFractureView, BrittleTetFractureControls, BrittleTetFractureEnergy,
    BrittleTetFractureError, ReleasedBrittleFace, require_released_faces, split_completed_fracture,
};
use crate::brittle::topology::TetFractureTopology;
use crate::loads3::TetLoadGeometry3;

/// Complete thermal solid on the released mesh, without a temperature change.
#[derive(Debug)]
pub struct ThermalBrittleTetFractureStep<'a> {
    pub solid: ThermalBrittleTetSolid<'a>,
    pub energy: BrittleTetFractureEnergy,
}

impl<'a> ThermalBrittleTetSolid<'a> {
    /// Retained matching traces for subsequent motion and explicit contact.
    #[must_use]
    pub fn released_faces(&self) -> &[ReleasedBrittleFace] {
        &self.released_faces
    }

    /// Release selected conforming faces after both native bands fully fail.
    ///
    /// The topology must originate from this body's current reference mesh.
    /// Both incident mechanical histories must be traction-free and aligned
    /// with the selected face. The original shared release gates preserve
    /// source ownership, total and mechanical strains, physical stress, band
    /// widths, storage and cumulative fracture. A cut does no physical work
    /// and neither advances temperature nor accepts another fracture step.
    ///
    /// Previous face guards are remapped through element-local corners. Future
    /// static and inertial endpoints retain their paired opening checks.
    /// Old equilibrium reactions are discarded because nodal numbering changes;
    /// callers must supply supports and forces in the new numbering. No finite
    /// sliding, automatic crack selection or partial-damage insertion is added.
    pub fn split_fracture<'next>(
        &self,
        topology: &TetFractureTopology<'_>,
        target_geometry: &'next TetLoadGeometry3<'next>,
        controls: BrittleTetFractureControls,
        cx: &Cx<'_>,
    ) -> Result<ThermalBrittleTetFractureStep<'next>, BrittleTetFractureError>
    where
        'a: 'next,
    {
        let transfer = split_completed_fracture(
            BrittleFractureView {
                geometry: self.geometry,
                bands: &self.bands,
                released_faces: &self.released_faces,
                budget: self.budget,
                solution: &self.solution,
                stored_energy_j: self.stored_energy_j,
                fracture_dissipation_j: self.fracture_dissipation_j,
            },
            |element| self.materials[element].band().law(),
            topology,
            target_geometry,
            controls,
            cx,
        )?;
        Ok(ThermalBrittleTetFractureStep {
            solid: ThermalBrittleTetSolid {
                geometry: target_geometry,
                materials: self.materials,
                bands: transfer.bands,
                released_faces: transfer.released_faces,
                budget: self.budget,
                solution: transfer.solution,
                stored_energy_j: transfer.energy.stored_after_j,
                fracture_dissipation_j: transfer.energy.fracture_dissipation_after_j,
            },
            energy: transfer.energy,
        })
    }

    /// Shared geometric guard independent of total/mechanical strain convention.
    pub(crate) fn require_released_faces(
        &self,
        displacement_m: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<(), BrittleTetError> {
        require_released_faces(
            self.geometry.nodes_m().len(),
            &self.released_faces,
            displacement_m,
            cx,
        )
    }
}
