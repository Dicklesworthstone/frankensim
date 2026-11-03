//! Independent momentum, heat and species fits on an original wall channel.
//!
//! Only homogeneous tangential planes are filtered. The source owns molecular
//! wall stress, heat and exactly zero species wall flux. Impermeability alone
//! does not prescribe a mass-fraction gradient in a coupled Maxwell--Stefan /
//! Soret law: both physical and filtered species derivatives therefore use
//! three interior fluid centres near each wall, without invented wall values.

use super::thermal::ThermalFields;
use super::*;
use crate::compressible::spatial::molecular::sgs::{channel_axis, channel_wall_values};
use fs_lbm::sgs3::dynamic::channel::components::ChannelDynamicComponents;
pub use fs_lbm::sgs3::dynamic::channel::components::ChannelDynamicComponentsConfig;
use fs_lbm::sgs3::dynamic::channel::{ChannelDynamicScalar, ChannelDynamicSgs};

/// Independent plane-averaged momentum, thermal and common species transport.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MixtureChannelDynamicSpeciesSubgridConfig {
    pub transport: ChannelDynamicThermalSubgridConfig,
    pub species: ChannelDynamicComponentsConfig,
    /// Finite nonnegative physical bound on the fitted common D [m²/s].
    /// Coupled molecular diffusion is not reduced to an invented scalar ratio.
    pub maximum_species_diffusivity_m2_s: f64,
}

impl MixtureChannelDynamicSpeciesSubgridConfig {
    pub fn validate(self) -> Result<(), DissipationError> {
        self.transport.validate().map_err(|_| {
            DissipationError::Invalid("invalid channel mixture momentum/heat controls")
        })?;
        self.species.validate().map_err(dynamic_error)?;
        if !self.maximum_species_diffusivity_m2_s.is_finite()
            || self.maximum_species_diffusivity_m2_s < 0.0
        {
            return Err(DissipationError::Invalid(
                "finite nonnegative channel species diffusivity bound required",
            ));
        }
        Ok(())
    }
}

impl<L: TransportLaw> MolecularTransport3d<L> {
    /// Complete original flow/operator storage plus 352 nominal numerical
    /// bytes/cell for independent channel fits, without species-count growth.
    pub fn required_workspace_bytes_with_channel_dynamic_species_sgs(
        flow: &MixtureEuler3d,
        source: &SpatialTransport<L>,
    ) -> Result<usize, DissipationError> {
        Self::required_workspace_bytes_with_dynamic_species_sgs(flow, source)
    }

    /// Attach the shared channel fit to the original molecular/mixture owner.
    /// Exactly two axes must be periodic; the remaining lower/upper Euler sides
    /// must reflect and receive the explicitly supplied no-slip wall pair.
    /// Every dimension needs at least three cells. Each wall may translate
    /// tangentially and use its original fixed-temperature or adiabatic law.
    ///
    /// Every SSP stage and retry uses current rho/u/T and ALL mass fractions.
    /// Maxwell--Stefan, Soret/Dufour, mass-frame closure, source enthalpies and
    /// the conservative energy/stability owner are unchanged. Modeled wall
    /// fluxes remain zero. This is resolved-wall exploratory gradient diffusion,
    /// without wall functions, catalytic chemistry or turbulent reaction rates.
    #[allow(clippy::too_many_arguments)]
    pub fn with_channel_dynamic_species_sgs(
        flow: &MixtureEuler3d,
        source: SpatialTransport<L>,
        config: MixtureChannelDynamicSpeciesSubgridConfig,
        walls: [MolecularWall3d; 2],
        maximum_workspace_bytes: usize,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, DissipationError> {
        poll(&mut cancelled)?;
        config.validate()?;
        let normal = channel_axis(&flow.grid, flow.boundaries).map_err(dynamic_error)?;
        if walls
            .iter()
            .any(|wall| matches!(wall, MolecularWall3d::SlipAdiabatic))
        {
            return Err(DissipationError::Invalid(
                "channel SGS requires explicit no-slip walls",
            ));
        }
        admit_workspace(
            Self::required_workspace_bytes_with_channel_dynamic_species_sgs(flow, &source)?,
            maximum_workspace_bytes,
        )?;
        let mut result = Self::new(flow, source, maximum_workspace_bytes, &mut cancelled)?;
        for (upper, wall) in [false, true].into_iter().zip(walls) {
            result = result.with_wall(side_for(normal, upper), wall)?;
        }
        result.sgs = Some(SubgridTransport {
            closure: ViscosityClosure::ChannelDynamicSpecies(ChannelDynamicSpeciesTransport {
                config,
                momentum: ChannelDynamicSgs::new_cartesian(
                    flow.grid.shape,
                    flow.grid.spacing,
                    normal,
                    ChannelDynamicSgs::required_bytes(flow.cell_count()).map_err(dynamic_error)?,
                )
                .map_err(dynamic_error)?,
                heat: ChannelDynamicScalar::new(
                    flow.grid.shape,
                    normal,
                    ChannelDynamicScalar::required_bytes(flow.cell_count())
                        .map_err(dynamic_error)?,
                )
                .map_err(dynamic_error)?,
                species: ChannelDynamicComponents::new(
                    flow.grid.shape,
                    normal,
                    flow.gas.species_count(),
                    ChannelDynamicComponents::required_bytes(flow.cell_count())
                        .map_err(dynamic_error)?,
                )
                .map_err(dynamic_error)?,
                fields: ThermalFields::new(flow.cell_count())?,
                rho_diffusivity: allocated(flow.cell_count(), 0.0)?,
            }),
            dynamic_viscosity: allocated(flow.cell_count(), 0.0)?,
            heat_capacity_v_density: allocated(flow.cell_count(), 0.0)?,
        });
        let bytes = result
            .workspace_bytes()
            .checked_add(flow.workspace_bytes())
            .ok_or(DissipationError::Invalid(
                "channel species workspace overflow",
            ))?;
        admit_workspace(bytes, maximum_workspace_bytes)?;
        poll(&mut cancelled)?;
        Ok(result)
    }

    pub fn channel_dynamic_species_subgrid_config(
        &self,
    ) -> Option<MixtureChannelDynamicSpeciesSubgridConfig> {
        match &self.sgs.as_ref()?.closure {
            ViscosityClosure::ChannelDynamicSpecies(dynamic) => Some(dynamic.config),
            _ => None,
        }
    }
}

pub(in crate::compressible::spatial::mixture::molecular::sgs) struct ChannelDynamicSpeciesTransport
{
    pub(in crate::compressible::spatial::mixture::molecular::sgs) config:
        MixtureChannelDynamicSpeciesSubgridConfig,
    momentum: ChannelDynamicSgs,
    heat: ChannelDynamicScalar,
    species: ChannelDynamicComponents,
    fields: ThermalFields,
    rho_diffusivity: Vec<f64>,
}

impl ChannelDynamicSpeciesTransport {
    pub(in crate::compressible::spatial::mixture::molecular::sgs) fn workspace_bytes(
        &self,
    ) -> usize {
        self.momentum.workspace_bytes()
            + self.heat.workspace_bytes()
            + self.species.workspace_bytes()
            + self.fields.workspace_bytes()
            + self.rho_diffusivity.capacity() * size_of::<f64>()
    }

    pub(in crate::compressible::spatial::mixture::molecular::sgs) fn prepare(
        &mut self,
        states: &[f64],
        primitive: &[Primitive],
        walls: &[MolecularWall3d; 6],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<bool, DissipationError> {
        let (wall_velocity, heat_walls) =
            channel_wall_values(self.momentum.normal_axis(), walls).map_err(dynamic_error)?;
        self.fields.load(primitive, cancelled)?;
        let width = self.species.component_count() + 4;
        let (momentum, heat, species) = self
            .momentum
            .prepare_with_scalar_and_components(
                self.config.transport.momentum,
                self.config.transport.heat,
                self.config.species,
                &self.fields.density,
                &self.fields.velocity,
                wall_velocity,
                &self.fields.temperature,
                heat_walls,
                &mut self.heat,
                |species, cell| states[cell * width + species] / primitive[cell].density,
                &mut self.species,
                &mut *cancelled,
            )
            .map_err(dynamic_error)?;
        if species.maximum_diffusivity > self.config.maximum_species_diffusivity_m2_s {
            return Err(DissipationError::Invalid(
                "channel species diffusivity exceeds its physical admission bound",
            ));
        }
        self.fields
            .collect_heat(self.heat.diffusivities(), cancelled)?;
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
                    "unrepresentable channel species mass diffusion",
                ));
            }
        }
        poll(cancelled)?;
        Ok(momentum.maximum_eddy_viscosity != 0.0
            || heat.maximum_diffusivity != 0.0
            || species.maximum_diffusivity != 0.0)
    }

    pub(in crate::compressible::spatial::mixture::molecular::sgs) fn viscosities(&self) -> &[f64] {
        self.momentum.viscosities()
    }

    pub(in crate::compressible::spatial::mixture::molecular::sgs) fn face_heat(
        &self,
        left: usize,
        right: usize,
    ) -> (f64, f64) {
        (
            self.fields.face_heat(left, right),
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
