//! Common species diffusion fitted from every resolved mass fraction.
//!
//! The shared Germano fit minimizes the sum of all species' scalar residuals,
//! with one normalization for the complete species bank. Its common positive
//! diffusivity preserves the mass-frame closure `sum(j_i) = 0` without fitting
//! or clipping species independently. Heat, species and momentum remain three
//! independent fits: species mixing can survive a zero momentum coefficient.
//! This is a gradient-diffusion exploration closure, not a turbulent reaction
//! rate, flame-wrinkling, wall or shock interaction model.

use super::thermal::ThermalWorkspace;
use super::*;
pub use fs_lbm::sgs3::dynamic::scalar::components::DynamicComponentsConfig;
use fs_lbm::sgs3::dynamic::scalar::components::PeriodicDynamicComponents;

/// Independent momentum, heat and common species diffusivity controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixtureDynamicSpeciesSubgridConfig {
    pub transport: DynamicThermalSubgridConfig,
    pub species: DynamicComponentsConfig,
    /// Finite nonnegative admission bound on the fitted species D [m²/s].
    /// The original coupled Maxwell--Stefan law has no single molecular D;
    /// its total stability bound is deliberately not used as a fake ratio.
    pub maximum_species_diffusivity_m2_s: f64,
}

impl MixtureDynamicSpeciesSubgridConfig {
    pub fn validate(self) -> Result<(), DissipationError> {
        self.transport.validate().map_err(|_| {
            DissipationError::Invalid("invalid dynamic mixture momentum/heat controls")
        })?;
        self.species.validate().map_err(dynamic_error)?;
        if !self.maximum_species_diffusivity_m2_s.is_finite()
            || self.maximum_species_diffusivity_m2_s < 0.0
        {
            return Err(DissipationError::Invalid(
                "finite nonnegative dynamic species diffusivity bound required",
            ));
        }
        Ok(())
    }
}

impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Combined flow/operator storage. Fitting all species adds 88 bytes/cell
    /// over independently fitted momentum/heat, for nominally 352 bytes/cell
    /// over molecular-only flow. Storage does not grow with the species count.
    pub fn required_workspace_bytes_with_dynamic_species_sgs(
        flow: &MixtureEuler3d,
        source: &SpatialTransport<L>,
    ) -> Result<usize, DissipationError> {
        let components =
            PeriodicDynamicComponents::required_bytes(flow.cell_count()).map_err(dynamic_error)?;
        let weights =
            flow.cell_count()
                .checked_mul(size_of::<f64>())
                .ok_or(DissipationError::Invalid(
                    "dynamic species workspace overflow",
                ))?;
        Self::required_workspace_bytes_with_dynamic_thermal_sgs(flow, source)?
            .checked_add(components)
            .and_then(|bytes| bytes.checked_add(weights))
            .ok_or(DissipationError::Invalid(
                "dynamic species workspace overflow",
            ))
    }

    /// Fit all three transports from the same original SSP-stage fields.
    /// Every species participates with its actual mass fraction and the same
    /// normalization; no marker species or turbulent Schmidt number is chosen.
    /// The original conservative faces and reacting split remain authoritative.
    pub fn with_dynamic_species_sgs(
        flow: &MixtureEuler3d,
        source: SpatialTransport<L>,
        config: MixtureDynamicSpeciesSubgridConfig,
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        poll(&mut cancelled)?;
        config.validate()?;
        admit_domain(flow)?;
        admit_workspace(
            Self::required_workspace_bytes_with_dynamic_species_sgs(flow, &source)?,
            maximum_workspace_bytes,
        )?;
        let mut result = Self::new(flow, source, maximum_workspace_bytes, &mut cancelled)?;
        result.sgs = Some(SubgridTransport {
            closure: ViscosityClosure::DynamicSpecies(DynamicSpeciesTransport {
                config,
                work: ThermalWorkspace::new(&flow.grid)?,
                species: PeriodicDynamicComponents::new(
                    flow.grid.shape,
                    flow.gas.species_count(),
                    PeriodicDynamicComponents::required_bytes(flow.cell_count())
                        .map_err(dynamic_error)?,
                )
                .map_err(dynamic_error)?,
                rho_diffusivity: allocated(flow.cell_count(), 0.0)?,
            }),
            dynamic_viscosity: allocated(flow.cell_count(), 0.0)?,
            heat_capacity_v_density: allocated(flow.cell_count(), 0.0)?,
        });
        let bytes = result
            .workspace_bytes()
            .checked_add(flow.workspace_bytes())
            .ok_or(DissipationError::Invalid(
                "dynamic species workspace overflow",
            ))?;
        admit_workspace(bytes, maximum_workspace_bytes)?;
        poll(&mut cancelled)?;
        Ok(result)
    }

    /// Active independent species controls, without a prescribed Schmidt ratio.
    pub fn dynamic_species_subgrid_config(&self) -> Option<MixtureDynamicSpeciesSubgridConfig> {
        match &self.sgs.as_ref()?.closure {
            ViscosityClosure::DynamicSpecies(dynamic) => Some(dynamic.config),
            _ => None,
        }
    }
}

pub(in crate::compressible::spatial::mixture::molecular::sgs) struct DynamicSpeciesTransport {
    pub(in crate::compressible::spatial::mixture::molecular::sgs) config:
        MixtureDynamicSpeciesSubgridConfig,
    work: ThermalWorkspace,
    species: PeriodicDynamicComponents,
    rho_diffusivity: Vec<f64>,
}

impl DynamicSpeciesTransport {
    pub(in crate::compressible::spatial::mixture::molecular::sgs) fn workspace_bytes(
        &self,
    ) -> usize {
        self.work.workspace_bytes()
            + self.species.workspace_bytes()
            + self.rho_diffusivity.capacity() * size_of::<f64>()
    }

    pub(in crate::compressible::spatial::mixture::molecular::sgs) fn prepare(
        &mut self,
        states: &[f64],
        primitive: &[Primitive],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<bool, DissipationError> {
        self.work.load(primitive, cancelled)?;
        let width = self.species.component_count() + 4;
        let (momentum, heat, species) = self
            .work
            .momentum
            .prepare_with_scalar_and_components(
                self.config.transport.momentum,
                self.config.transport.heat,
                self.config.species,
                &self.work.density,
                &self.work.velocity,
                &self.work.temperature,
                &mut self.work.heat,
                |species, cell| states[cell * width + species] / primitive[cell].density,
                &mut self.species,
                &mut *cancelled,
            )
            .map_err(dynamic_error)?;
        if species.maximum_diffusivity > self.config.maximum_species_diffusivity_m2_s {
            return Err(DissipationError::Invalid(
                "dynamic species diffusivity exceeds its physical admission bound",
            ));
        }
        self.work.collect_heat(cancelled)?;
        for (cell, (&diffusivity, weight)) in self
            .species
            .diffusivities()
            .iter()
            .zip(&mut self.rho_diffusivity)
            .enumerate()
        {
            if cell % 64 == 0 {
                poll(cancelled)?;
            }
            *weight = primitive[cell].density * diffusivity;
            if !weight.is_finite() {
                return Err(DissipationError::Invalid(
                    "unrepresentable dynamic species mass diffusion",
                ));
            }
        }
        poll(cancelled)?;
        Ok(momentum.maximum_eddy_viscosity != 0.0
            || heat.maximum_diffusivity != 0.0
            || species.maximum_diffusivity != 0.0)
    }

    pub(in crate::compressible::spatial::mixture::molecular::sgs) fn viscosities(&self) -> &[f64] {
        self.work.momentum.viscosities()
    }

    pub(in crate::compressible::spatial::mixture::molecular::sgs) fn face_heat(
        &self,
        left: usize,
        right: usize,
    ) -> (f64, f64) {
        (
            self.work.face_heat(left, right),
            self.config.transport.heat.maximum_diffusivity_to_molecular,
        )
    }

    pub(in crate::compressible::spatial::mixture::molecular::sgs) fn face_species(
        &self,
        left: usize,
        right: usize,
    ) -> f64 {
        0.5 * self.rho_diffusivity[left] + 0.5 * self.rho_diffusivity[right]
    }
}
