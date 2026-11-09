//! Adaptive callbacks retain the exact native reservoir/face-temperature binding.

use super::super::super::enthalpy::adaptive::SpeciesEnthalpyDiffusionModel;
use super::super::super::{
    CancelGate, ChtError, MaterialAdvectionStep, SpeciesDiffusionBoundary, VoxelDomain,
};
use super::{GasSpeciesEnthalpyBoundarySource, GasSpeciesEnthalpyField};

impl<'s, 'a, const N: usize> SpeciesEnthalpyDiffusionModel<N>
    for GasSpeciesEnthalpyBoundarySource<'s, 'a, N>
{
    type Prepared = GasSpeciesEnthalpyField<'a, N>;

    fn prepare(
        &self,
        domain: &VoxelDomain,
        carrier: &MaterialAdvectionStep,
        species_mass_kg: &[[f64; N]],
        specific_enthalpy_j_kg: &[f64],
        gate: &CancelGate,
    ) -> Result<Self::Prepared, ChtError> {
        GasSpeciesEnthalpyBoundarySource::prepare(
            self,
            domain,
            carrier,
            species_mass_kg,
            specific_enthalpy_j_kg,
            gate,
        )
    }

    fn temperature_k<'p>(&self, prepared: &'p Self::Prepared) -> &'p [f64] {
        prepared.temperatures_k()
    }

    fn diffusivity_m2_s<'p>(&self, prepared: &'p Self::Prepared) -> &'p [f64] {
        prepared.diffusivity_m2_s()
    }

    fn face_partial_enthalpy_j_kg(
        &self,
        prepared: &Self::Prepared,
        domain: &VoxelDomain,
        cells: [usize; 2],
        axis: usize,
        gate: &CancelGate,
    ) -> Result<[f64; N], ChtError> {
        self.face_partials(prepared, domain, cells, axis, gate)
    }

    fn boundary_partial_enthalpy_j_kg(
        &self,
        prepared: &Self::Prepared,
        domain: &VoxelDomain,
        boundary_index: usize,
        boundary: &SpeciesDiffusionBoundary<N>,
        gate: &CancelGate,
    ) -> Result<[f64; N], ChtError> {
        super::poll(gate)?;
        if self.boundaries().get(boundary_index) != Some(boundary) {
            return Err(super::invalid(
                "adaptive gas reservoir differs from its bound face-temperature descriptor",
            ));
        }
        self.boundary_partials(prepared, domain, boundary_index, gate)
    }

    fn maximum_prepared_bytes(&self, cells: usize) -> Option<usize> {
        GasSpeciesEnthalpyBoundarySource::maximum_prepared_bytes(self, cells).ok()
    }

    fn caloric_defect_j(&self, prepared: &Self::Prepared) -> Result<f64, ChtError> {
        Ok(prepared.absolute_caloric_residual_j())
    }
}
