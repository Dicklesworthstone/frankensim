//! Spatially resolved material absorption on the original density cells.
//!
//! A common wavelength grid and a borrowed material map supply each cell's
//! physical mass absorption. Extinction is integrated through the actual
//! layers before the shared finite-pixel footprint averages transmission.

use super::{
    ColumnAbsorption, DensityAbsorptionSpectrum, DvrError, MAX_OPACITY_WAVELENGTHS,
    SpectralAbsorptionImageSettings, SpectralOpacityImage, Sum, TransmissionAverage,
    TransportAverage, VolumeGrid, optical_product, render,
};

/// Borrowed spectral material palette and original per-cell material indices.
///
/// `mass_absorption_m2_kg[material * wavelengths_nm.len() + wavelength]` is
/// intensity mass absorption in m²/kg. Every row uses the same one to sixteen
/// positive, strictly increasing vacuum wavelengths. Coefficients are finite
/// and nonnegative; no interpolation, measured-data authority, material mixing
/// law or dependence on temperature is inferred.
///
/// `cell_materials` has exactly the density grid's original x-fastest layout
/// and length. Each index selects one material/state row, including in vacuum
/// cells. Density remains `grid.data[cell] * density_unit_kg_m3`; a material
/// index does not replace or rescale the physical density. All slices are
/// borrowed, with no converted extinction field or expanded spectral grid.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IndexedDensityAbsorptionSpectrum<'a> {
    /// Common positive, strictly increasing vacuum wavelength samples [nm].
    pub wavelengths_nm: &'a [f64],
    /// Material-major, wavelength-fastest mass absorption coefficients [m²/kg].
    pub mass_absorption_m2_kg: &'a [f64],
    /// One admitted palette-row index per original density cell.
    pub cell_materials: &'a [u32],
    /// Positive density represented by one stored field unit [kg/m³].
    pub density_unit_kg_m3: f64,
    /// Positive metres represented by one grid coordinate unit.
    pub length_unit_m: f64,
}

/// Explicit image, material-palette and wavelength evaluation limits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct IndexedSpectralAbsorptionImageSettings {
    /// Shared image limits. Spectral work counts intersected original cells
    /// times wavelengths, including zero density and zero absorption.
    pub spectral: SpectralAbsorptionImageSettings,
    /// Maximum complete material rows, including rows unused by the map.
    pub max_materials: usize,
}

/// Exact spectral opacity of spatially varying absorbing materials.
///
/// Each ray column integrates `sum(rho_cell * kappa_material * cell_length)`
/// separately at every supplied wavelength. Its transmitted intensity is then
/// averaged over the exact pixel footprint. Density, material coefficients and
/// optical depth are never averaged across distinct rays before attenuation.
///
/// The geometry, averaging, output layout and range-stable optical products
/// are the same owners used by uniform spectral opacity. The material map and
/// complete palette are admitted before output allocation; invalid indices are
/// refused even when their cells have zero density. Scattering, interfaces,
/// emission and subcell material arrangements are outside this model.
pub fn render_density_spectral_indexed_box(
    grid: &VolumeGrid<'_>,
    settings: IndexedSpectralAbsorptionImageSettings,
    spectrum: IndexedDensityAbsorptionSpectrum<'_>,
) -> Result<SpectralOpacityImage, DvrError> {
    let samples = render::<TransportAverage, _>(
        grid,
        settings.spectral,
        IndexedColumns {
            spectrum,
            max_materials: settings.max_materials,
        },
    )?;
    Ok(SpectralOpacityImage::from_samples(
        settings.spectral.image.resolution,
        spectrum.wavelengths_nm,
        samples,
    ))
}

/// The same indexed-material transport, retaining only direct transmission.
///
/// Writes pixel-major, wavelength-fastest values directly, without a full
/// opacity/depth image. Every cell uses its own original material index; the
/// per-column exponential and footprint filtering are shared with the uniform
/// transmission path.
pub fn render_density_spectral_indexed_transmittance_box(
    grid: &VolumeGrid<'_>,
    settings: IndexedSpectralAbsorptionImageSettings,
    spectrum: IndexedDensityAbsorptionSpectrum<'_>,
) -> Result<Vec<f64>, DvrError> {
    render::<TransmissionAverage, _>(
        grid,
        settings.spectral,
        IndexedColumns {
            spectrum,
            max_materials: settings.max_materials,
        },
    )
}

pub(super) struct IndexedColumns<'a> {
    pub(super) spectrum: IndexedDensityAbsorptionSpectrum<'a>,
    pub(super) max_materials: usize,
}

impl ColumnAbsorption for IndexedColumns<'_> {
    fn wavelengths(&self) -> &[f64] {
        self.spectrum.wavelengths_nm
    }

    fn cell_coefficients(&self, cell: usize) -> &[f64] {
        let lanes = self.spectrum.wavelengths_nm.len();
        let first = self.spectrum.cell_materials[cell] as usize * lanes;
        &self.spectrum.mass_absorption_m2_kg[first..first + lanes]
    }

    fn physical_units(&self) -> (f64, f64) {
        (
            self.spectrum.density_unit_kg_m3,
            self.spectrum.length_unit_m,
        )
    }

    fn admit(&self, grid: &VolumeGrid<'_>) -> Result<(), DvrError> {
        let lanes = self.spectrum.wavelengths_nm.len();
        // Reuse the bounded uniform spectral/unit admission on the first row.
        // The remaining palette values and cell indices are scanned only after
        // the shared image and spectral work admission has succeeded.
        DensityAbsorptionSpectrum {
            wavelengths_nm: self.spectrum.wavelengths_nm,
            mass_absorption_m2_kg: self
                .spectrum
                .mass_absorption_m2_kg
                .get(..lanes)
                .ok_or(DvrError::InvalidAbsorptionSpectrum)?,
            density_unit_kg_m3: self.spectrum.density_unit_kg_m3,
            length_unit_m: self.spectrum.length_unit_m,
        }
        .validate()?;
        if self.spectrum.mass_absorption_m2_kg.len() % lanes != 0 {
            return Err(DvrError::InvalidAbsorptionSpectrum);
        }
        let materials = self.spectrum.mass_absorption_m2_kg.len() / lanes;
        if materials > self.max_materials {
            return Err(DvrError::AbsorptionMaterialBudgetExceeded {
                requested: materials,
                limit: self.max_materials,
            });
        }
        if self.spectrum.cell_materials.len() != grid.data.len() {
            return Err(DvrError::InvalidAbsorptionMaterialMap);
        }
        Ok(())
    }

    fn validate_values(&self) -> Result<(), DvrError> {
        if self
            .spectrum
            .mass_absorption_m2_kg
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(DvrError::InvalidAbsorptionSpectrum);
        }
        let materials =
            self.spectrum.mass_absorption_m2_kg.len() / self.spectrum.wavelengths_nm.len();
        if let Some((cell, &material)) = self
            .spectrum
            .cell_materials
            .iter()
            .enumerate()
            .find(|(_, material)| **material as usize >= materials)
        {
            return Err(DvrError::AbsorptionMaterialIndexOutOfBounds {
                cell,
                material,
                materials,
            });
        }
        Ok(())
    }

    fn spectral_steps(&self, depth_cells: usize) -> usize {
        depth_cells
    }

    fn visit_column(
        &self,
        grid: &VolumeGrid<'_>,
        column: usize,
        plane: usize,
        mut visit: impl FnMut(usize, f64) -> Result<(), DvrError>,
    ) -> Result<(), DvrError> {
        let lanes = self.spectrum.wavelengths_nm.len();
        let mut depths = [Sum::default(); MAX_OPACITY_WAVELENGTHS];
        for iz in (0..grid.dims[2]).rev() {
            let cell = column + plane * iz;
            let first = self.spectrum.cell_materials[cell] as usize * lanes;
            let coefficients = &self.spectrum.mass_absorption_m2_kg[first..first + lanes];
            for (lane, &coefficient) in coefficients.iter().enumerate() {
                depths[lane].add(optical_product(&[
                    grid.data[cell],
                    coefficient,
                    self.spectrum.density_unit_kg_m3,
                    grid.cell[2],
                    self.spectrum.length_unit_m,
                ]));
            }
        }
        for (lane, depth) in depths.into_iter().take(lanes).enumerate() {
            visit(lane, depth.value)?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests;
