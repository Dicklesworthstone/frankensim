//! Different accepted particle materials in the original physical grid cells.

use super::{MaterialParticleGrid, MaterialParticleSpectrum, ParticleScatteringError};
use crate::volumes::VolumeGrid;
use crate::volumes::transport::scattering::{ScatteringGrid, ScatteringOptics};
use crate::volumes::transport::{
    AbsorbingGrid, AbsorbingGridBudget, CellPalette, CellPaletteBudget, VolumeTransportError,
};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;

/// Original prepared material spectra and one compact material index per cell.
/// Density and optional LTE temperature remain borrowed simulation buffers.
pub struct MaterialParticleCellGrid<'a> {
    materials: CellPalette<MaterialParticleSpectrum>,
    grid: MaterialParticleGrid<'a>,
}

impl<'a> MaterialParticleCellGrid<'a> {
    /// Bind spatially varying accepted particle materials, temperatures or size
    /// distributions to index-one ambient transport. `materials` is the original
    /// prepared spectrum palette; its cell map uses the density grid's x-major
    /// order, and must have exactly one entry per cell. Palette preparation owns
    /// its explicit model-count and map budgets before this call.
    ///
    /// Each cell uses its mapped material's absorption, scattering and angular
    /// law. The two native optical palettes share the SAME immutable index map;
    /// no per-ray material search or substitute transport loop is introduced.
    /// All palette entries must have an index-one host.
    ///
    /// Density is dispersed particle mass in kg/m³. Each cell's volume fraction
    /// is checked against its own accepted particle bulk density and the caller's
    /// `maximum_volume_fraction` in (0,1]. This admission limit does not establish
    /// independent-particle accuracy or account for particle correlations.
    /// Optional LTE temperature must equal the mapped original material endpoint
    /// in every cell; omitting the field declares a nonemitting model.
    pub fn try_new(
        cx: &Cx<'_>,
        grid: VolumeGrid<'a>,
        materials: CellPalette<MaterialParticleSpectrum>,
        temperature_k: Option<&'a [f64]>,
        maximum_volume_fraction: f64,
        budget: AbsorbingGridBudget,
    ) -> Result<Self, ParticleScatteringError> {
        checkpoint(cx)?;
        if grid.data.len() > budget.max_grid_cells {
            return Err(VolumeTransportError::Budget.into());
        }
        if materials.cell_indices().len() != grid.data.len() {
            return Err(VolumeTransportError::Palette.into());
        }
        if !maximum_volume_fraction.is_finite()
            || maximum_volume_fraction <= 0.0
            || maximum_volume_fraction > 1.0
        {
            return Err(ParticleScatteringError::VolumeFractionLimit);
        }
        if temperature_k.is_some_and(|values| values.len() != grid.data.len()) {
            return Err(VolumeTransportError::Temperature.into());
        }
        let mut maximum_observed_volume_fraction = 0.0_f64;
        for (index, &density) in grid.data.iter().enumerate() {
            if index % 1024 == 0 {
                checkpoint(cx)?;
            }
            let material = materials
                .at_cell(index)
                .ok_or(VolumeTransportError::Palette)?;
            if !density.is_finite() || density < 0.0 {
                return Err(VolumeTransportError::Density { index }.into());
            }
            let fraction = density / material.density_kg_m3;
            if fraction > maximum_volume_fraction {
                return Err(ParticleScatteringError::VolumeFraction { index });
            }
            maximum_observed_volume_fraction = maximum_observed_volume_fraction.max(fraction);
            if let Some(temperatures) = temperature_k {
                let accepted = material
                    .temperature_k
                    .ok_or(ParticleScatteringError::ThermalState)?;
                if temperatures[index] != accepted {
                    return Err(ParticleScatteringError::TemperatureMismatch { index });
                }
            }
        }
        let palette_budget = CellPaletteBudget {
            max_models: materials.values().len(),
            max_cells: materials.cell_indices().len(),
        };
        let mut absorption = Vec::with_capacity(palette_budget.max_models);
        let mut scattering = Vec::with_capacity(palette_budget.max_models);
        let mut hash = DomainHasher::new("org.frankensim.render.material-particle-cell-grid.v1");
        hash.update(&(materials.values().len() as u64).to_le_bytes());
        for material in materials.values() {
            checkpoint(cx)?;
            if material.spectrum.host_index_n != 1.0 {
                return Err(ParticleScatteringError::HostIndex);
            }
            absorption.push(material.spectrum.absorption);
            scattering.push(ScatteringOptics {
                scattering: material.spectrum.scattering,
                phase: material.spectrum.phase.clone(),
            });
            hash.update(material.binding_identity.as_bytes());
        }
        let absorption = materials.with_values(cx, absorption, palette_budget)?;
        let scattering = materials.with_values(cx, scattering, palette_budget)?;
        let absorption =
            AbsorbingGrid::try_new_indexed(cx, grid, absorption, temperature_k, budget)?;
        let medium = ScatteringGrid::try_new_indexed(cx, absorption, scattering)?;
        hash.update(medium.content_identity().as_bytes());
        hash.update(&maximum_volume_fraction.to_bits().to_le_bytes());
        Ok(Self {
            materials,
            grid: MaterialParticleGrid {
                medium,
                identity: hash.finalize(),
                maximum_observed_volume_fraction,
            },
        })
    }

    /// Original prepared spectra, source identities and shared material map.
    #[must_use]
    pub const fn materials(&self) -> &CellPalette<MaterialParticleSpectrum> {
        &self.materials
    }

    /// Native transport input for rendering, photometry and beam deposition.
    #[must_use]
    pub const fn medium(&self) -> &ScatteringGrid<'a> {
        self.grid.medium()
    }

    /// Original source palette, numerical medium and declared dilution limit.
    #[must_use]
    pub const fn content_identity(&self) -> ContentHash {
        self.grid.content_identity()
    }

    /// Largest cell fraction using each mapped material's own bulk density.
    #[must_use]
    pub const fn maximum_volume_fraction(&self) -> f64 {
        self.grid.maximum_volume_fraction()
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), ParticleScatteringError> {
    cx.checkpoint()
        .map_err(|_| VolumeTransportError::Cancelled.into())
}
