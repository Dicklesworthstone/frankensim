//! Endpoint caloric closure of the original two-stage finite-air diffusion.

use super::*;

/// Native layered SDIRK stages and one complete isothermal heat account.
pub type LayeredSdirk2SorptionEnergyCandidate =
    LayeredSorptionEnergyCandidate<MoistureLayeredSlab, MoistureLayeredSdirk2Report>;
/// Native coated sphere/cylinder stages and one complete heat account.
pub type LayeredRadialSdirk2SorptionEnergyCandidate =
    LayeredSorptionEnergyCandidate<MoistureLayeredRadial, MoistureLayeredSdirk2Report>;

impl FiniteMoistureAtmosphere {
    /// Solve both native material/gas SDIRK stages, then account the actual
    /// original-to-final bound storage and finite-gas energy once. Every layer
    /// keeps its own original loading-dependent isosteric heat source.
    ///
    /// The first stage is a transport diagnostic, not an additional heat or
    /// water debit. Heat integration uses the accepted final cell inventory and
    /// one source-segment budget over all layers. Internal redistribution is
    /// retained even when the gas receives no water. Source, work, energy or
    /// cancellation refusal discards the complete private native candidate.
    ///
    /// Temperature remains prescribed by the same external thermostat. This
    /// adds no spatial heat equation or inferred local heat deposition. Smooth
    /// native source branches retain the transport solver's second-order scope.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_layered_sdirk2_with_sorption_energy(
        &self,
        body: &MoistureLayeredSlab,
        dt_s: f64,
        exposed: [bool; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        atmosphere: AtmosphereMoistureSolve,
        laws: &[MoistureSorptionEnergy],
        controls: SorptionEnergyControls,
        cx: &Cx<'_>,
    ) -> Result<LayeredSdirk2SorptionEnergyCandidate> {
        admit_layers(
            body.layers().iter().map(|layer| layer.slab().field()),
            laws,
            controls,
            cx,
        )?;
        let native =
            self.advance_layered_sdirk2(body, dt_s, exposed, transport, nonlinear, atmosphere, cx)?;
        let energy = account_layers(
            self,
            body.layers().iter().map(|layer| layer.slab().field()),
            native
                .moisture()
                .layers()
                .iter()
                .map(|layer| layer.slab().field()),
            laws,
            native.step().gas_internal_energy_change_j,
            controls,
            cx,
        )?;
        Ok(LayeredSorptionEnergyCandidate { native, energy })
    }

    /// The same closed caloric account for native two-stage coated-sphere or
    /// coated-cylinder diffusion. Each original shell contributes its actual
    /// dry mass, not a slab or equal-volume surrogate. Gas energy enters once.
    ///
    /// Held-temperature scope, source work and atomic failure semantics follow
    /// [`Self::advance_layered_sdirk2_with_sorption_energy`].
    #[allow(clippy::too_many_arguments)]
    pub fn advance_layered_radial_sdirk2_with_sorption_energy(
        &self,
        body: &MoistureLayeredRadial,
        dt_s: f64,
        exposed: [bool; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        atmosphere: AtmosphereMoistureSolve,
        laws: &[MoistureSorptionEnergy],
        controls: SorptionEnergyControls,
        cx: &Cx<'_>,
    ) -> Result<LayeredRadialSdirk2SorptionEnergyCandidate> {
        admit_layers(
            body.layers().iter().map(|layer| layer.radial().field()),
            laws,
            controls,
            cx,
        )?;
        let native = self.advance_layered_radial_sdirk2(
            body, dt_s, exposed, transport, nonlinear, atmosphere, cx,
        )?;
        let energy = account_layers(
            self,
            body.layers().iter().map(|layer| layer.radial().field()),
            native
                .moisture()
                .layers()
                .iter()
                .map(|layer| layer.radial().field()),
            laws,
            native.step().gas_internal_energy_change_j,
            controls,
            cx,
        )?;
        Ok(LayeredSorptionEnergyCandidate { native, energy })
    }
}
