//! Independent pathline memory for mixture heat and common species mixing.
//!
//! Every mass fraction participates in the common species contraction. The
//! callback reads each original cell's species row directly; no species-by-cell
//! transpose or separately drifting species histories are stored. All three
//! fitted coefficients share one accepted trajectory and remain held during
//! a complete physical split, while the conservative faces use current strain,
//! density, source Cp and referenced partial enthalpies.

use super::*;
pub use fs_lbm::sgs3::dynamic::lagrangian::scalar::LagrangianScalarConfig;
pub use fs_lbm::sgs3::dynamic::lagrangian::scalar::components::LagrangianComponentsConfig;

/// Independent physical-time memories for stress, heat and common species flux.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixtureLagrangianSpeciesSubgridConfig {
    pub momentum: LagrangianSgsConfig,
    pub heat: LagrangianScalarConfig,
    pub species: LagrangianComponentsConfig,
    /// Finite nonnegative bound on cell and conservative-face D [m²/s].
    /// The original Maxwell--Stefan matrix has no single molecular diffusivity;
    /// its coupled stability bound is not used as a surrogate Schmidt ratio.
    pub maximum_species_diffusivity_m2_s: f64,
}

impl MixtureLagrangianSpeciesSubgridConfig {
    pub fn validate(self) -> Result<(), DissipationError> {
        self.momentum.validate().map_err(dynamic::dynamic_error)?;
        self.heat.validate().map_err(dynamic::dynamic_error)?;
        self.species.validate().map_err(dynamic::dynamic_error)?;
        if !self.maximum_species_diffusivity_m2_s.is_finite()
            || self.maximum_species_diffusivity_m2_s < 0.0
        {
            return Err(DissipationError::Invalid(
                "finite nonnegative Lagrangian species diffusivity bound required",
            ));
        }
        Ok(())
    }
}

impl SubgridTransport {
    /// Shared histories/filter: 376 bytes/cell; rho/u/T staging, three held
    /// coefficients and the existing source rho*Cv cache add 72 bytes/cell.
    pub(in crate::compressible::spatial::mixture::molecular) fn required_lagrangian_species_bytes(
        flow: &MixtureEuler3d,
    ) -> Result<usize, DissipationError> {
        LagrangianDynamicComponents::required_bytes(flow.cell_count())
            .map_err(dynamic::dynamic_error)?
            .checked_add(flow.cell_count().checked_mul(9 * size_of::<f64>()).ok_or(
                DissipationError::Invalid("Lagrangian mixture species workspace size overflow"),
            )?)
            .ok_or(DissipationError::Invalid(
                "Lagrangian mixture species workspace size overflow",
            ))
    }

    pub(in crate::compressible::spatial::mixture::molecular) fn new_lagrangian_species(
        flow: &MixtureEuler3d,
        config: MixtureLagrangianSpeciesSubgridConfig,
    ) -> Result<Self, DissipationError> {
        config.validate()?;
        Self::new_lagrangian_owner(flow, LagrangianConfig::Species(config), None)
    }

    pub(in crate::compressible::spatial::mixture::molecular) fn lagrangian_species_config(
        &self,
    ) -> Option<MixtureLagrangianSpeciesSubgridConfig> {
        match &self.closure {
            ViscosityClosure::Lagrangian(dynamic) => dynamic.species_config(),
            _ => None,
        }
    }

    pub(in crate::compressible::spatial::mixture::molecular) fn lagrangian_heat_report(
        &self,
    ) -> Option<DynamicScalarReport> {
        match &self.closure {
            ViscosityClosure::Lagrangian(dynamic) => dynamic.last_heat_report,
            _ => None,
        }
    }

    pub(in crate::compressible::spatial::mixture::molecular) fn lagrangian_species_report(
        &self,
    ) -> Option<DynamicScalarReport> {
        match &self.closure {
            ViscosityClosure::Lagrangian(dynamic) => dynamic.last_species_report,
            _ => None,
        }
    }
}

impl LagrangianTransport {
    fn species_config(&self) -> Option<MixtureLagrangianSpeciesSubgridConfig> {
        match self.config {
            LagrangianConfig::FixedTransport(_) => None,
            LagrangianConfig::Species(config) => Some(config),
        }
    }

    pub(in crate::compressible::spatial::mixture::molecular::sgs) fn face_heat_coefficient(
        &self,
        left: usize,
        right: usize,
    ) -> Option<(f64, f64)> {
        self.species_config().map(|config| {
            (
                0.5 * self.heat_coefficient[left] + 0.5 * self.heat_coefficient[right],
                config.heat.maximum_diffusivity_to_molecular,
            )
        })
    }

    pub(in crate::compressible::spatial::mixture::molecular::sgs) fn face_species_coefficient(
        &self,
        left: usize,
        right: usize,
    ) -> Option<(f64, f64)> {
        self.species_config().map(|config| {
            (
                0.5 * self.species_coefficient[left] + 0.5 * self.species_coefficient[right],
                config.maximum_species_diffusivity_m2_s,
            )
        })
    }
}
