//! Source-preserving temporal fluid changes in the original native DG owner.
//!
//! These adapters select complete material responses. The existing shared
//! mass, Bessel-arm, wall-history and parameter-work transfer does the update.

use super::super::gas_update::poll;
use super::*;

impl<M: FluidAcoustics + Clone> DgDuct<DgFluidMaterial<M>> {
    /// Rebind a continuing uniform fluid while retaining its original source.
    ///
    /// The prescribed, instantaneous material intervention preserves acoustic
    /// compression `p/K_s` and axial momentum `rho*q/A`. Geometry, ports,
    /// accepted time, prior work/dissipation, and every elastic/Kelvin/Maxwell
    /// wall coordinate survive. The resulting storage change is signed external
    /// parameter work, charged exactly once through the existing native ledger.
    ///
    /// Source/query/phase metadata are retained even when the new coefficients
    /// equal the old ones. A changed operator invalidates its old propagator;
    /// prepare the native propagator again before continuing. All admission,
    /// work-budget and cancellation checks precede publishing either owner.
    ///
    /// This is not a heat, mass-injection, mean-flow or mean-pressure wall-load
    /// solve. Viscothermal memory requires the dedicated method below; spatial
    /// and molecular history transformations are not inferred for a generic
    /// fluid. There is no change of geometry, topology or physical time here.
    pub fn rebind_uniform_material(
        &mut self,
        state: &mut DgState,
        material: &M,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgMaterialUpdate, FluidDuctError> {
        poll(cancelled)?;
        let acoustic = material.acoustic_state()?;
        let medium = DgFluidMaterial {
            source: material.clone(),
            acoustic,
            thermoviscous: None,
        };
        Ok(self.rebind_uniform_medium_data(
            state,
            medium,
            None,
            maximum_absolute_parameter_work_j,
            cancelled,
        )?)
    }
}

impl<M: ThermoviscousFluid + Clone> DgDuct<DgFluidMaterial<M>> {
    /// Rebind the same live fluid WITH its original circular diffusion history.
    ///
    /// The existing physical Bessel-arm transfer preserves series momentum and
    /// thermal-arm compression volume, using the actually frozen old Cp/Cv.
    /// Pressure, density, Cp/Cv, viscosity and conductivity may all change.
    /// Changing transport alone leaves instantaneous stored energy unchanged
    /// and changes subsequent diffusion. Structural histories remain intact.
    ///
    /// Original radius, pole count, angular band, mesh and ports are retained;
    /// the replacement state must independently resolve that band at the given
    /// error bound. Creating or removing thermal-arm support (Cp/Cv crossing
    /// exactly one), spatial memory, or a molecular model refuses. Failed
    /// source evaluation, cancellation or excess parameter work changes neither
    /// the original source/operator nor the continuing state.
    ///
    /// This uses the same native exponential after its propagator is rebuilt.
    /// It assigns no background thermodynamic heat, species exchange or static
    /// wall traction; the owner prescribing the change owes its parameter work.
    pub fn rebind_uniform_viscothermal_material(
        &mut self,
        state: &mut DgState,
        material: &M,
        maximum_relative_error: f64,
        maximum_absolute_parameter_work_j: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<DgMaterialUpdate, FluidDuctError> {
        poll(cancelled)?;
        let old_material = self.gas.thermoviscous.ok_or(DuctError::BadParameter {
            what: "live generic fluid diffusion requires its original frozen caloric and transport state",
        })?;
        let replacement = material.thermoviscous_state()?;
        let acoustic = material.acoustic_state()?;
        if replacement.acoustic() != acoustic {
            return Err(DuctError::BadParameter {
                what: "live DG thermoviscous and acoustic views must be the same original material state",
            }.into());
        }
        let medium = DgFluidMaterial {
            source: material.clone(),
            acoustic,
            thermoviscous: Some(replacement),
        };
        Ok(self.rebind_uniform_viscothermal_medium_data(
            state,
            medium,
            None,
            old_material,
            replacement,
            replacement.transport().prandtl(),
            maximum_relative_error,
            maximum_absolute_parameter_work_j,
            cancelled,
        )?)
    }
}
