//! Original gas closure for the state-refreshed native diffusion driver.

use super::super::enthalpy::adaptive::SpeciesEnthalpyDiffusionModel;
use super::{
    CancelGate, ChtError, GasSpeciesEnthalpyField, GasSpeciesEnthalpySource, MaterialAdvectionStep,
    VoxelDomain,
};

impl<'a, const N: usize> SpeciesEnthalpyDiffusionModel<N> for GasSpeciesEnthalpySource<'a, N> {
    type Prepared = GasSpeciesEnthalpyField<'a, N>;

    fn prepare(
        &self,
        domain: &VoxelDomain,
        carrier: &MaterialAdvectionStep,
        species_mass_kg: &[[f64; N]],
        specific_enthalpy_j_kg: &[f64],
        gate: &CancelGate,
    ) -> Result<Self::Prepared, ChtError> {
        GasSpeciesEnthalpySource::prepare(
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

    fn maximum_prepared_bytes(&self, cells: usize) -> Option<usize> {
        GasSpeciesEnthalpySource::maximum_prepared_bytes(self, cells).ok()
    }

    fn caloric_defect_j(&self, prepared: &Self::Prepared) -> Result<f64, ChtError> {
        Ok(prepared.absolute_caloric_residual_j())
    }
}
