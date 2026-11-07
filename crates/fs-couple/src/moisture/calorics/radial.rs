//! Physical isothermal energy closure over actual curved-shell water storage.

use super::*;
use crate::moisture::{MoistureSdirk2Report, radial::MoistureRadial};

/// Original radial transport, finite vapor and sourced thermostat duty.
pub type RadialSorptionEnergyCandidate = SorptionEnergyCandidate<MoistureRadial>;
/// Both original radial SDIRK stages with one original-to-final energy account.
pub type RadialSdirk2SorptionEnergyCandidate =
    SorptionEnergyCandidate<MoistureRadial, MoistureSdirk2Report>;

impl FiniteMoistureAtmosphere {
    /// Resolve finite-vapor radial diffusion and its actual isothermal heat duty.
    ///
    /// Each shell's original dry mass weights the source-integrated differential
    /// isosteric heat. Nonuniform loading and internal redistribution therefore
    /// remain visible even when the net vapor exchange is small or zero. The
    /// gas and bound-water energy increments share the original vapor datum.
    /// Source, work or closed-energy refusal publishes neither candidate.
    ///
    /// This is the external heat needed to hold prescribed temperature. Local
    /// bound-water storage changes are not unique heat deposition; no thermal
    /// diffusion, latent-heat default or temperature feedback is inferred.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_radial_with_sorption_energy(
        &self,
        field: &MoistureRadial,
        dt_s: f64,
        faces: &[Option<DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        atmosphere: AtmosphereMoistureSolve,
        law: &MoistureSorptionEnergy,
        controls: SorptionEnergyControls,
        cx: &Cx<'_>,
    ) -> Result<RadialSorptionEnergyCandidate> {
        admit(field.field(), law, controls, cx)?;
        let native =
            self.advance_radial(field, dt_s, faces, transport, nonlinear, atmosphere, cx)?;
        finish(self, field, native, law, controls, cx)
    }

    /// The same actual storage-energy closure after conservative radial SDIRK2.
    ///
    /// The full accepted endpoint determines one physical isothermal duty. The
    /// intermediate stage is diagnostic and is not added as another energy or
    /// water transfer. The source integral retains its exact piecewise-linear
    /// loading dependence; temporal accuracy of the transport remains separate.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_radial_sdirk2_with_sorption_energy(
        &self,
        field: &MoistureRadial,
        dt_s: f64,
        faces: &[Option<DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        atmosphere: AtmosphereMoistureSolve,
        law: &MoistureSorptionEnergy,
        controls: SorptionEnergyControls,
        cx: &Cx<'_>,
    ) -> Result<RadialSdirk2SorptionEnergyCandidate> {
        admit(field.field(), law, controls, cx)?;
        let native =
            self.advance_radial_sdirk2(field, dt_s, faces, transport, nonlinear, atmosphere, cx)?;
        finish(self, field, native, law, controls, cx)
    }
}
