//! Conservative incident-beam absorption on the original density cells.
//!
//! Irradiance, wavelength integration and physical column area determine
//! incident watts. Absorption is evaluated before any spatial reduction.
//! This is incident absorbed power, without scattering, thermal emission,
//! temperature advancement or an inferred net radiative heat balance.

use super::{
    ColumnAbsorption, DensityAbsorptionSpectrum, DvrError, IndexedDensityAbsorptionSpectrum,
    MAX_OPACITY_WAVELENGTHS, Sum, VolumeGrid, indexed::IndexedColumns, optical_product,
    validate_transfer_grid,
};
use fs_math::det;

/// Explicit normal irradiance at the original vacuum wavelength nodes.
///
/// Values are finite nonnegative W/(m² nm). They are not camera radiance;
/// converting radiance to irradiance requires a separate angular integral.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum SpectralBeamIrradiance<'a> {
    /// One spectrum over the entire illuminated XY face.
    Uniform(&'a [f64]),
    /// One spectrum per original XY column, with X fastest and wavelength
    /// innermost. Each spectrum is constant across its column's physical area.
    Columns(&'a [f64]),
}

impl SpectralBeamIrradiance<'_> {
    fn validate(self, columns: usize, lanes: usize) -> Result<(), DvrError> {
        let (values, expected) = match self {
            Self::Uniform(values) => (values, lanes),
            Self::Columns(values) => (
                values,
                columns
                    .checked_mul(lanes)
                    .ok_or(DvrError::InvalidBeamIrradianceSize)?,
            ),
        };
        if values.len() != expected {
            return Err(DvrError::InvalidBeamIrradianceSize);
        }
        if values
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
        {
            return Err(DvrError::InvalidBeamIrradiance);
        }
        Ok(())
    }

    fn at(self, column: usize, lane: usize, lanes: usize) -> f64 {
        match self {
            Self::Uniform(values) => values[lane],
            Self::Columns(values) => values[column * lanes + lane],
        }
    }
}

/// Physical propagation through the complete original Z extent.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BeamPropagation {
    /// Enter at the minimum-Z face, leave at the maximum-Z face.
    PositiveZ,
    /// Enter at the maximum-Z face, leave at the minimum-Z face.
    NegativeZ,
}

/// Admission limits for a complete dense physical power report.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SpectralBeamDepositionSettings {
    /// Maximum original density cells, including vacuum cells.
    pub max_grid_cells: usize,
    /// Maximum original cells times wavelengths, including zero-power lanes.
    pub max_spectral_evaluations: u64,
    /// Maximum allocated map values: original cells plus original XY columns.
    /// Three fixed total-power slots and bounded wavelength metadata are extra.
    pub max_output_values: usize,
}

/// Complete beam-work limits plus an explicit material-palette bound.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct IndexedSpectralBeamDepositionSettings {
    pub deposition: SpectralBeamDepositionSettings,
    /// Maximum supplied material rows, including unused rows.
    pub max_materials: usize,
}

/// Spatial absorbed watts and the complete transmitted incident-beam balance.
///
/// Every wavelength obeys incident = absorbed + transmitted to floating-point
/// roundoff. Maps retain the original X-fastest grid order. No camera sampling,
/// exposure, display response, time step or heat capacity enters this report.
/// The spectral integral is only the caller's supplied quadrature; no missing
/// wavelength response or broadband-source coverage is inferred.
#[derive(Clone, Debug, PartialEq)]
pub struct SpectralBeamDeposition {
    dimensions: [usize; 3],
    wavelengths_nm: [f64; MAX_OPACITY_WAVELENGTHS],
    quadrature_weights_nm: [f64; MAX_OPACITY_WAVELENGTHS],
    lanes: usize,
    propagation: BeamPropagation,
    incident: [f64; MAX_OPACITY_WAVELENGTHS],
    absorbed: [f64; MAX_OPACITY_WAVELENGTHS],
    transmitted: [f64; MAX_OPACITY_WAVELENGTHS],
    // One allocation serves the typed maps and the browser/export buffer.
    power_buffer: Vec<f64>,
}

impl SpectralBeamDeposition {
    #[must_use]
    pub const fn dimensions(&self) -> [usize; 3] {
        self.dimensions
    }

    #[must_use]
    pub fn wavelengths_nm(&self) -> &[f64] {
        &self.wavelengths_nm[..self.lanes]
    }

    #[must_use]
    pub fn quadrature_weights_nm(&self) -> &[f64] {
        &self.quadrature_weights_nm[..self.lanes]
    }

    #[must_use]
    pub const fn propagation(&self) -> BeamPropagation {
        self.propagation
    }

    /// Original per-cell absorbed power, summed over supplied spectral lanes [W].
    #[must_use]
    pub fn cell_absorbed_power_w(&self) -> &[f64] {
        &self.power_buffer[3..3 + self.cell_count()]
    }

    /// Power leaving each original XY column, summed over spectral lanes [W].
    #[must_use]
    pub fn column_transmitted_power_w(&self) -> &[f64] {
        &self.power_buffer[3 + self.cell_count()..]
    }

    /// Whole-face incident power per already integrated spectral lane [W].
    #[must_use]
    pub fn incident_power_w(&self) -> &[f64] {
        &self.incident[..self.lanes]
    }

    /// Whole-field absorbed power per already integrated spectral lane [W].
    #[must_use]
    pub fn absorbed_power_w(&self) -> &[f64] {
        &self.absorbed[..self.lanes]
    }

    /// Whole-face transmitted power per already integrated spectral lane [W].
    #[must_use]
    pub fn transmitted_power_w(&self) -> &[f64] {
        &self.transmitted[..self.lanes]
    }

    #[must_use]
    pub fn total_incident_power_w(&self) -> f64 {
        self.power_buffer[0]
    }

    #[must_use]
    pub fn total_absorbed_power_w(&self) -> f64 {
        self.power_buffer[1]
    }

    #[must_use]
    pub fn total_transmitted_power_w(&self) -> f64 {
        self.power_buffer[2]
    }

    /// Transfer the single allocation for native/binary/browser consumption.
    /// Layout is [total incident W, total absorbed W, total transmitted W,
    /// cell absorbed W..., column transmitted W...]. No map copy is made.
    #[must_use]
    pub fn into_power_buffer(self) -> Vec<f64> {
        self.power_buffer
    }

    fn cell_count(&self) -> usize {
        self.dimensions[0] * self.dimensions[1] * self.dimensions[2]
    }
}

/// Deposit a normal collimated spectral beam in a uniformly composed volume.
///
/// Explicit irradiance-per-nanometre, nonnegative quadrature widths and the
/// original cell's SI XY area determine incident watts. Each cell absorbs
/// P_incident*exp(-prefix_depth)*(1-exp(-cell_depth)); compensated cumulative
/// depth retains transmission through many individually sub-ULP thin cells.
/// The same physical density/length scales as spectral imaging are used.
///
/// Returned cell watts can supply an admitted thermal power port. They do not
/// advance temperature or assume that absorbed optical energy immediately
/// thermalizes. Emission and scattered radiation are separate physical owners.
/// All grid, material, source and work admission precedes output allocation;
/// arithmetic refusal returns no partial report.
pub fn deposit_density_spectral_beam(
    grid: &VolumeGrid<'_>,
    settings: SpectralBeamDepositionSettings,
    spectrum: DensityAbsorptionSpectrum<'_>,
    irradiance: SpectralBeamIrradiance<'_>,
    quadrature_weights_nm: &[f64],
    propagation: BeamPropagation,
) -> Result<SpectralBeamDeposition, DvrError> {
    deposit(
        grid,
        settings,
        spectrum,
        irradiance,
        quadrature_weights_nm,
        propagation,
    )
}

/// The same conservative beam deposition with each original cell's material.
/// Coefficients may be prepared mixtures; their density mass basis stays explicit.
pub fn deposit_density_spectral_indexed_beam(
    grid: &VolumeGrid<'_>,
    settings: IndexedSpectralBeamDepositionSettings,
    spectrum: IndexedDensityAbsorptionSpectrum<'_>,
    irradiance: SpectralBeamIrradiance<'_>,
    quadrature_weights_nm: &[f64],
    propagation: BeamPropagation,
) -> Result<SpectralBeamDeposition, DvrError> {
    deposit(
        grid,
        settings.deposition,
        IndexedColumns {
            spectrum,
            max_materials: settings.max_materials,
        },
        irradiance,
        quadrature_weights_nm,
        propagation,
    )
}

fn deposit<S: ColumnAbsorption>(
    grid: &VolumeGrid<'_>,
    settings: SpectralBeamDepositionSettings,
    spectrum: S,
    irradiance: SpectralBeamIrradiance<'_>,
    quadrature_weights_nm: &[f64],
    propagation: BeamPropagation,
) -> Result<SpectralBeamDeposition, DvrError> {
    spectrum.admit(grid)?;
    let lanes = spectrum.wavelengths().len();
    let cells = grid.data.len();
    if cells > settings.max_grid_cells {
        return Err(DvrError::GridCellBudgetExceeded {
            requested: cells,
            limit: settings.max_grid_cells,
        });
    }
    let [nx, ny, nz] = grid.dims;
    let columns = nx
        .checked_mul(ny)
        .filter(|_| !grid.dims.contains(&0))
        .ok_or(DvrError::InvalidGrid)?;
    if columns.checked_mul(nz) != Some(cells) {
        return Err(DvrError::InvalidGrid);
    }
    let evaluations = (cells as u128)
        .checked_mul(lanes as u128)
        .and_then(|value| u64::try_from(value).ok())
        .ok_or(DvrError::SpectralEvaluationCountOverflow)?;
    if evaluations > settings.max_spectral_evaluations {
        return Err(DvrError::SpectralEvaluationBudgetExceeded {
            requested: evaluations,
            limit: settings.max_spectral_evaluations,
        });
    }
    let map_values = cells
        .checked_add(columns)
        .ok_or(DvrError::SpectralEvaluationCountOverflow)?;
    if map_values > settings.max_output_values {
        return Err(DvrError::BeamOutputBudgetExceeded {
            requested: map_values,
            limit: settings.max_output_values,
        });
    }
    if quadrature_weights_nm.len() != lanes
        || quadrature_weights_nm
            .iter()
            .any(|value| !value.is_finite() || *value < 0.0)
    {
        return Err(DvrError::InvalidSpectralWeights);
    }
    spectrum.validate_values()?;
    validate_transfer_grid(grid)?;
    if let Some(index) = grid.data.iter().position(|value| *value < 0.0) {
        return Err(DvrError::NegativeGridAbsorption { index });
    }
    irradiance.validate(columns, lanes)?;
    let values = map_values
        .checked_add(3)
        .ok_or(DvrError::SpectralEvaluationCountOverflow)?;
    let mut power_buffer = Vec::new();
    power_buffer
        .try_reserve_exact(values)
        .map_err(|_| DvrError::BeamAllocationRefused { values })?;
    power_buffer.resize(values, 0.0);
    let mut incident = [Sum::default(); MAX_OPACITY_WAVELENGTHS];
    let mut absorbed = [Sum::default(); MAX_OPACITY_WAVELENGTHS];
    let mut transmitted = [Sum::default(); MAX_OPACITY_WAVELENGTHS];
    let (density_unit, length_unit) = spectrum.physical_units();
    for column in 0..columns {
        let mut entering = [0.0; MAX_OPACITY_WAVELENGTHS];
        let mut depths = [Sum::default(); MAX_OPACITY_WAVELENGTHS];
        for lane in 0..lanes {
            // Keep every original scale until multiplication is range-safe.
            entering[lane] = optical_product(&[
                irradiance.at(column, lane, lanes),
                quadrature_weights_nm[lane],
                grid.cell[0],
                grid.cell[1],
                length_unit,
                length_unit,
            ]);
            add_power(&mut incident[lane], entering[lane])?;
        }
        for step in 0..nz {
            let z = match propagation {
                BeamPropagation::PositiveZ => step,
                BeamPropagation::NegativeZ => nz - 1 - step,
            };
            let cell = column + columns * z;
            let coefficients = spectrum.cell_coefficients(cell);
            let mut cell_power = Sum::default();
            for lane in 0..lanes {
                let factors = [
                    grid.data[cell],
                    coefficients[lane],
                    density_unit,
                    grid.cell[2],
                    length_unit,
                ];
                let tau = optical_product(&factors);
                let power = absorbed_power(entering[lane], depths[lane].value, tau, &factors);
                add_power(&mut cell_power, power)?;
                add_power(&mut absorbed[lane], power)?;
                depths[lane].add(tau);
            }
            power_buffer[3 + cell] = cell_power.value;
        }
        let mut column_power = Sum::default();
        for lane in 0..lanes {
            let power = attenuated_power(entering[lane], depths[lane].value);
            add_power(&mut column_power, power)?;
            add_power(&mut transmitted[lane], power)?;
        }
        power_buffer[3 + cells + column] = column_power.value;
    }
    for (slot, totals) in [incident, absorbed, transmitted].into_iter().enumerate() {
        let mut total = Sum::default();
        for value in &totals[..lanes] {
            add_power(&mut total, value.value)?;
        }
        power_buffer[slot] = total.value;
    }
    let mut wavelengths_nm = [0.0; MAX_OPACITY_WAVELENGTHS];
    wavelengths_nm[..lanes].copy_from_slice(spectrum.wavelengths());
    let mut weights = [0.0; MAX_OPACITY_WAVELENGTHS];
    weights[..lanes].copy_from_slice(quadrature_weights_nm);
    Ok(SpectralBeamDeposition {
        dimensions: grid.dims,
        wavelengths_nm,
        quadrature_weights_nm: weights,
        lanes,
        propagation,
        incident: incident.map(|sum| sum.value),
        absorbed: absorbed.map(|sum| sum.value),
        transmitted: transmitted.map(|sum| sum.value),
        power_buffer,
    })
}

fn add_power(sum: &mut Sum, value: f64) -> Result<(), DvrError> {
    if !value.is_finite() || value < 0.0 {
        return Err(DvrError::UnrepresentableBeamPower);
    }
    sum.add(value);
    if sum.value.is_finite() {
        Ok(())
    } else {
        Err(DvrError::UnrepresentableBeamPower)
    }
}

fn attenuated_power(incident: f64, depth: f64) -> f64 {
    if incident == 0.0 {
        return 0.0;
    }
    if depth == 0.0 {
        return incident;
    }
    let transmission = det::exp(-depth);
    if transmission >= f64::MIN_POSITIVE {
        incident * transmission
    } else {
        det::exp(det::ln(incident) - depth)
    }
}

fn absorbed_power(incident: f64, prefix: f64, tau: f64, factors: &[f64; 5]) -> f64 {
    if incident == 0.0 || factors.contains(&0.0) {
        return 0.0;
    }
    if tau < f64::MIN_POSITIVE {
        // In this range 1-exp(-tau)=tau to binary64 precision. Preserve its
        // original positive factors: a bright beam can expose a rounded-zero tau.
        return det::exp(
            det::ln(incident) - prefix + factors.iter().map(|value| det::ln(*value)).sum::<f64>(),
        );
    }
    attenuated_power(incident, prefix) * -det::expm1(-tau)
}

#[cfg(test)]
mod tests;
