//! Independently fitted pathline heat and common species mixing on the
//! original accepted molecular--Euler--molecular physical interval.

use super::*;

impl<L: TransportLaw> MolecularEuler3d<L> {
    /// Additional numerical storage for all three histories: 448 bytes/cell.
    /// Add this to the original coupled or stabilized workspace requirement,
    /// plus supplied forcing scratch. Storage is independent of species count.
    pub fn required_lagrangian_dynamic_species_sgs_workspace_bytes(
        flow: &MixtureEuler3d,
    ) -> Result<usize, DissipationError> {
        sgs::SubgridTransport::required_lagrangian_species_bytes(flow)
    }

    /// Bind independent Lagrangian stress, heat and common species diffusion.
    ///
    /// The accepted physical START state supplies density, velocity,
    /// temperature and every mass fraction to one shared history owner.
    /// Separate signed contractions and memory times determine the three
    /// coefficients. No prescribed turbulent Prandtl/Schmidt ratio is inferred.
    /// The first-order pathline update uses the previous actual accepted
    /// duration, including shortening; all three histories publish together.
    ///
    /// Both molecular halves, every SSP/RKL stage and each whole-split retry
    /// hold the same coefficients. Existing conservative faces evaluate current
    /// strain, density, composition, source Cp and referenced partial enthalpies,
    /// enforce viscosity/heat/species bounds, and include all added diffusion
    /// in the original stability limit. There is no separate heating source.
    ///
    /// The supplied molecular operator must have no other SGS closure.
    /// The cap covers the complete pair and all history/staging capacities;
    /// preparation and stepping allocate no history fields or species matrix.
    /// Fully periodic Cartesian flow only. This starts fresh history and adds
    /// no wall, turbulent reaction-rate, SGS energy or shock-interaction model.
    pub fn with_lagrangian_dynamic_species_sgs(
        flow: MixtureEuler3d,
        molecular: MolecularTransport3d<L>,
        config: sgs::MixtureLagrangianSpeciesSubgridConfig,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        poll(&mut cancelled)?;
        config.validate()?;
        if molecular.sgs.is_some() {
            return Err(DissipationError::Invalid(
                "Lagrangian mixture attachment requires a molecular operator without SGS",
            ));
        }
        let required = Self::required_workspace_bytes(&flow, &molecular)?
            .checked_add(Self::required_lagrangian_dynamic_species_sgs_workspace_bytes(&flow)?)
            .ok_or(DissipationError::Invalid(
                "coupled Lagrangian mixture species workspace size overflow",
            ))?;
        admit_workspace(required, maximum_workspace_bytes)?;
        let sgs = sgs::SubgridTransport::new_lagrangian_species(&flow, config)?;
        poll(&mut cancelled)?;
        let mut result = Self::new(flow, molecular, maximum_workspace_bytes, &mut cancelled)?;
        result.molecular.sgs = Some(sgs);
        admit_workspace(result.workspace_bytes(), maximum_workspace_bytes)?;
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Independent memory controls, or None for fixed-ratio and other closures.
    pub fn lagrangian_dynamic_species_subgrid_config(
        &self,
    ) -> Option<sgs::MixtureLagrangianSpeciesSubgridConfig> {
        self.molecular
            .sgs
            .as_ref()
            .and_then(sgs::SubgridTransport::lagrangian_species_config)
    }

    /// Last committed accepted-START temperature fit, independently of momentum.
    /// Uses the same sample time as `lagrangian_subgrid_sample_time`.
    /// None for fixed-Prandtl controls or until a complete split accepts.
    pub fn lagrangian_subgrid_heat_report(&self) -> Option<sgs::DynamicScalarReport> {
        self.molecular
            .sgs
            .as_ref()
            .and_then(sgs::SubgridTransport::lagrangian_heat_report)
    }

    /// Last committed common-species fit, independently of momentum and heat.
    /// Cell diffusivity diagnostics do not replace current-face admission.
    /// None for fixed-Schmidt controls or until a complete split accepts.
    pub fn lagrangian_subgrid_species_report(&self) -> Option<sgs::DynamicScalarReport> {
        self.molecular
            .sgs
            .as_ref()
            .and_then(sgs::SubgridTransport::lagrangian_species_report)
    }
}
