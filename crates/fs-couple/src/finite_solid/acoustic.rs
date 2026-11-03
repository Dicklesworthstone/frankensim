//! Local finite-strain waves retain the original material and reference mass.
//!
//! Preparation uses the same total-F tangent as finite solid forces. Held
//! thermal references therefore retain their complete Fth pullback instead of
//! querying the extracted bare elastic law. Current density is rho0/det(F).
//! A body's separately bound reference-density field is outside this per-material
//! view; it cannot substitute the cell material's density for that body's mass.

use super::FiniteSolidMaterial;
use fs_material::MaterialError;
use fs_material::hyper::acoustic::{DirectionalHyperelasticWaves, HyperelasticAcousticState};

/// Original source binding and native acoustic snapshot at one total deformation.
/// Temperature and material history stay fixed; no adiabatic modulus, heat
/// evolution, attenuation or global stability certificate is introduced.
#[derive(Debug)]
pub struct FiniteSolidAcousticState<'a> {
    material: &'a FiniteSolidMaterial,
    native: HyperelasticAcousticState,
}

impl FiniteSolidMaterial {
    /// Prepare local waves using this exact source material's consistent
    /// tangent and original reference density. The thermal adapter, when
    /// present, owns elastic-domain admission through Fe=F Fth^-1.
    ///
    /// This is the per-material density binding. A body created with a separate
    /// source density field requires that field's local inertial density.
    pub fn acoustic_state(
        &self,
        total_deformation_gradient: &[f64; 9],
    ) -> Result<FiniteSolidAcousticState<'_>, MaterialError> {
        let native = HyperelasticAcousticState::from_material(
            self,
            total_deformation_gradient,
            self.reference_density_kg_m3(),
        )?;
        Ok(FiniteSolidAcousticState {
            material: self,
            native,
        })
    }
}

impl<'a> FiniteSolidAcousticState<'a> {
    /// Exact original model, density receipts, thermal path and named frames.
    pub const fn material(&self) -> &'a FiniteSolidMaterial {
        self.material
    }

    /// Native prepared tangent and conserved reference/current densities.
    pub const fn native(&self) -> &HyperelasticAcousticState {
        &self.native
    }

    /// Signed modes and energy velocities in the current spatial frame.
    /// Native direction, arithmetic and instability semantics remain unchanged.
    pub fn directional_waves(
        &self,
        direction: [f64; 3],
    ) -> Result<DirectionalHyperelasticWaves, MaterialError> {
        self.native.directional_waves(direction)
    }
}
