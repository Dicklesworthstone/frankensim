//! Absorption and local-equilibrium emission in borrowed simulation cells.
//!
//! Density is kg/m³, mass absorption is m²/kg, lengths are metres and the
//! source radiance is W/(m² sr nm). Each cell solves dL/ds = -rho*kappa*L +
//! rho*kappa*B_lambda(T) analytically. DDA visits only intersected cells;
//! there is no stochastic transmittance noise or fixed marching step size.
//! This is a piecewise-constant, non-scattering, index-one medium model.
//! An LBM lattice density requires an explicit physical scaling before entry.

use super::VolumeGrid;
use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN};
use crate::thermal_surface::planck_radiance_w_m2_sr_nm;
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_math::det;

mod cells;
/// Exact per-cell beam absorption and the separate scattered-power balance.
pub mod deposition;
/// Physical emitting-cell direction proposals and their solid-angle density.
pub mod emission;
mod palette;
pub use palette::{CellPalette, CellPaletteBudget};
/// Spectral free flights and phase functions for physical participating media.
pub mod scattering;
/// Immutable accepted medium snapshots on one physical observation clock.
pub mod timeline;

const MAX_KNOTS: usize = 16;

/// Validated mass absorption coefficient, in m²/kg, on an explicit wavelength
/// domain. Legacy visible constructors retain 380–780 nm; broader source tables
/// retain their own endpoints and never infer missing UV or infrared response.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MassAbsorptionSpectrum {
    knots: [(f64, f64); MAX_KNOTS],
    count: usize,
    source_identity: Option<ContentHash>,
}

impl MassAbsorptionSpectrum {
    /// Wavelength-independent absorption over 380–780 nm.
    pub fn constant(value_m2_kg: f64) -> Result<Self, VolumeTransportError> {
        Self::table(&[(LAMBDA_MIN, value_m2_kg), (LAMBDA_MAX, value_m2_kg)])
    }

    /// Linear interpolation of 2–16 `(vacuum_nm, m²/kg)` knots.
    /// The strictly increasing table must cover exactly 380–780 nm.
    pub fn table(knots: &[(f64, f64)]) -> Result<Self, VolumeTransportError> {
        if !(2..=MAX_KNOTS).contains(&knots.len())
            || knots[0].0 != LAMBDA_MIN
            || knots[knots.len() - 1].0 != LAMBDA_MAX
        {
            return Err(VolumeTransportError::AbsorptionSpectrum);
        }
        Self::admit(knots, None)
    }

    /// Admit an explicitly identified source over exactly its supplied support.
    ///
    /// The 2–16 knots are `(positive vacuum_nm, nonnegative m²/kg)`, strictly
    /// increasing in wavelength. Between knots the declared model is linear;
    /// callers must resolve spectral features adequately. The nonzero identity
    /// binds source/condition data but does not confer measurement authority.
    /// No value outside the first and last source knots is extrapolated.
    pub fn source_table(
        knots: &[(f64, f64)],
        source_identity: ContentHash,
    ) -> Result<Self, VolumeTransportError> {
        if source_identity == ContentHash([0; 32]) {
            return Err(VolumeTransportError::SourceIdentity);
        }
        Self::admit(knots, Some(source_identity))
    }

    fn admit(
        knots: &[(f64, f64)],
        source_identity: Option<ContentHash>,
    ) -> Result<Self, VolumeTransportError> {
        if !(2..=MAX_KNOTS).contains(&knots.len())
            || knots
                .iter()
                .any(|(l, k)| !l.is_finite() || *l <= 0.0 || !k.is_finite() || *k < 0.0)
            || knots.windows(2).any(|pair| pair[0].0 >= pair[1].0)
        {
            return Err(VolumeTransportError::AbsorptionSpectrum);
        }
        let mut result = Self {
            knots: [(0.0, 0.0); MAX_KNOTS],
            count: knots.len(),
            source_identity,
        };
        result.knots[..knots.len()].copy_from_slice(knots);
        Ok(result)
    }

    /// Exact inclusive support of the retained source/model, in vacuum nm.
    #[must_use]
    pub const fn wavelength_bounds_nm(&self) -> (f64, f64) {
        (self.knots[0].0, self.knots[self.count - 1].0)
    }

    /// Original explicitly supplied source identity, when present.
    #[must_use]
    pub const fn source_identity(&self) -> Option<ContentHash> {
        self.source_identity
    }

    /// Evaluate without extrapolating outside the admitted source domain.
    pub fn at(&self, wavelength_nm: f64) -> Result<f64, VolumeTransportError> {
        let (lower, upper) = self.wavelength_bounds_nm();
        if !wavelength_nm.is_finite() || !(lower..=upper).contains(&wavelength_nm) {
            return Err(VolumeTransportError::Wavelength);
        }
        let knots = &self.knots[..self.count];
        let upper = knots
            .partition_point(|(l, _)| *l < wavelength_nm)
            .clamp(1, self.count - 1);
        let (a, ka) = knots[upper - 1];
        let (b, kb) = knots[upper];
        if ka == kb {
            return Ok(ka);
        }
        let t = (wavelength_nm - a) / (b - a);
        let value = (1.0 - t) * ka + t * kb;
        if value.is_finite() {
            Ok(value)
        } else {
            Err(VolumeTransportError::Unrepresentable)
        }
    }
}

/// Caller-owned limits for admission and per-ray traversal.
#[derive(Clone, Copy, Debug)]
pub struct AbsorbingGridBudget {
    /// Maximum number of borrowed cells scanned at admission. Emitting fields
    /// may also retain at most one cumulative f64 weight per cell plus one.
    pub max_grid_cells: usize,
    /// Maximum visited cells per ray segment, including zero-length crossings.
    pub max_steps_per_segment: u32,
}

/// No partial transport result is returned on refusal or cancellation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VolumeTransportError {
    /// The caller cancelled this operation.
    Cancelled,
    /// Dimensions, bounds or cell sizes are invalid or not representable.
    Grid,
    /// Density must be finite and nonnegative, in kg/m³.
    Density {
        /// Invalid field index.
        index: usize,
    },
    /// Temperature length or a finite, strictly positive kelvin value is invalid.
    Temperature,
    /// Invalid spectral knots or negative/non-finite mass absorption.
    AbsorptionSpectrum,
    /// Invalid spectral knots or negative/non-finite mass scattering.
    ScatteringSpectrum,
    /// Empty optical palette, invalid model index, or cell-map length mismatch.
    Palette,
    /// An explicitly sourced optical table requires a nonzero source identity.
    SourceIdentity,
    /// Invalid phase asymmetry, direction, or scattering angle.
    Phase,
    /// An incident beam-power lane is negative or non-finite.
    BeamPower {
        /// Invalid spectral lane index.
        lane: usize,
    },
    /// A spectral packet is empty or a random variate is outside [0, 1).
    Sample,
    /// Wavelength is nonfinite, nonpositive or outside an admitted source domain.
    Wavelength,
    /// Ray direction is not finite and unit, or distance is negative/NaN.
    Ray,
    /// Admission storage/scan or per-segment work exceeds the caller's limit.
    Budget,
    /// Finite arithmetic cannot resolve a cell crossing or physical radiance.
    Unrepresentable,
}

impl core::fmt::Display for VolumeTransportError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "volume transport refused: {self:?}")
    }
}
impl core::error::Error for VolumeTransportError {}

/// Integrated foreground transport for a fixed wavelength packet.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SegmentTransport<const N: usize> {
    /// Multiplicative attenuation of radiance arriving from the far endpoint.
    pub transmittance: [f64; N],
    /// Foreground local-equilibrium emission, W/(m² sr nm).
    pub emission_w_m2_sr_nm: [f64; N],
}

/// Immutable, borrowed accepted density/temperature snapshot.
/// Its content identity binds numerical state and the declared optical model,
/// without asserting that the coefficients are measured specimen properties.
pub struct AbsorbingGrid<'a> {
    grid: VolumeGrid<'a>,
    absorption: MassAbsorptionSpectrum,
    absorption_palette: Option<CellPalette<MassAbsorptionSpectrum>>,
    temperature_k: Option<&'a [f64]>,
    max_steps: u32,
    identity: ContentHash,
    vacuum: bool,
    emission: emission::EmissionDistribution,
}

impl<'a> AbsorbingGrid<'a> {
    /// Number of original X-major physical density cells.
    #[must_use]
    pub fn cell_count(&self) -> usize {
        self.grid.data.len()
    }

    /// Borrowed physical mass density in kg/m³ at an original cell index.
    #[must_use]
    pub fn cell_density_kg_m3(&self, index: usize) -> Option<f64> {
        self.grid.data.get(index).copied()
    }

    /// Validate, hash and borrow the original simulation buffers; cache bounded
    /// emission proposal weights when a visible thermal source is present.
    pub fn try_new(
        cx: &Cx<'_>,
        grid: VolumeGrid<'a>,
        absorption: MassAbsorptionSpectrum,
        temperature_k: Option<&'a [f64]>,
        budget: AbsorbingGridBudget,
    ) -> Result<Self, VolumeTransportError> {
        Self::admit(cx, grid, absorption, None, temperature_k, budget)
    }

    /// Bind cell-resolved absorption to the original borrowed density and
    /// optional temperature buffers. The map must contain exactly one model
    /// index per native cell, and EVERY supplied model must have identical
    /// wavelength endpoints. Unreferenced palette entries never create haze
    /// or thermal emission, but remain part of the declaration's identity.
    pub fn try_new_indexed(
        cx: &Cx<'_>,
        grid: VolumeGrid<'a>,
        absorption: CellPalette<MassAbsorptionSpectrum>,
        temperature_k: Option<&'a [f64]>,
        budget: AbsorbingGridBudget,
    ) -> Result<Self, VolumeTransportError> {
        checkpoint(cx)?;
        if grid.data.len() > budget.max_grid_cells || budget.max_steps_per_segment == 0 {
            return Err(VolumeTransportError::Budget);
        }
        if absorption.cell_indices().len() != grid.data.len() {
            return Err(VolumeTransportError::Palette);
        }
        let representative = absorption.values()[0];
        for chunk in absorption.values().chunks(1024) {
            checkpoint(cx)?;
            if chunk
                .iter()
                .any(|value| value.wavelength_bounds_nm() != representative.wavelength_bounds_nm())
            {
                return Err(VolumeTransportError::AbsorptionSpectrum);
            }
        }
        Self::admit(
            cx,
            grid,
            representative,
            Some(absorption),
            temperature_k,
            budget,
        )
    }

    fn admit(
        cx: &Cx<'_>,
        grid: VolumeGrid<'a>,
        absorption: MassAbsorptionSpectrum,
        absorption_palette: Option<CellPalette<MassAbsorptionSpectrum>>,
        temperature_k: Option<&'a [f64]>,
        budget: AbsorbingGridBudget,
    ) -> Result<Self, VolumeTransportError> {
        checkpoint(cx)?;
        if grid.data.len() > budget.max_grid_cells || budget.max_steps_per_segment == 0 {
            return Err(VolumeTransportError::Budget);
        }
        let (lo, hi) = grid.bounds();
        if grid.dims.contains(&0)
            || grid.dims.iter().try_fold(1usize, |n, d| n.checked_mul(*d)) != Some(grid.data.len())
            || (0..3).any(|a| {
                !lo[a].is_finite()
                    || !hi[a].is_finite()
                    || !grid.cell[a].is_finite()
                    || grid.cell[a] <= 0.0
                    || !((hi[a] - lo[a]).is_finite() && hi[a] > lo[a])
                    || lo[a] + grid.cell[a] <= lo[a]
                    || hi[a] - grid.cell[a] >= hi[a]
            })
        {
            return Err(VolumeTransportError::Grid);
        }
        if temperature_k.is_some_and(|t| t.len() != grid.data.len()) {
            return Err(VolumeTransportError::Temperature);
        }
        let mut hash = DomainHasher::new("org.frankensim.render.absorbing-grid.v1");
        for axis in 0..3 {
            hash.update(&(grid.dims[axis] as u64).to_le_bytes());
            hash.update(&grid.origin[axis].to_bits().to_le_bytes());
            hash.update(&grid.cell[axis].to_bits().to_le_bytes());
        }
        if let Some(palette) = &absorption_palette {
            hash.update(b"indexed-absorption-v1");
            hash.update(&(palette.values().len() as u64).to_le_bytes());
            for chunk in palette.values().chunks(1024) {
                checkpoint(cx)?;
                for value in chunk {
                    update_absorption_identity(&mut hash, value);
                }
            }
            hash.update(&(palette.cell_indices().len() as u64).to_le_bytes());
            for chunk in palette.cell_indices().chunks(1024) {
                checkpoint(cx)?;
                for index in chunk {
                    hash.update(&index.to_le_bytes());
                }
            }
        } else {
            // Preserve the original uniform optical identity exactly.
            update_absorption_identity(&mut hash, &absorption);
        }
        hash.update(&[u8::from(temperature_k.is_some())]);
        let mut positive_density = false;
        let mut mapped_absorption = false;
        let mut previous_model = None;
        let mut model_absorbs = false;
        for (index, &density) in grid.data.iter().enumerate() {
            if index % 1024 == 0 {
                checkpoint(cx)?;
            }
            if !density.is_finite() || density < 0.0 {
                return Err(VolumeTransportError::Density { index });
            }
            positive_density |= density > 0.0;
            if !mapped_absorption
                && density > 0.0
                && let Some(palette) = &absorption_palette
            {
                let model = palette.cell_indices()[index];
                if previous_model != Some(model) {
                    let spectrum = &palette.values()[model as usize];
                    model_absorbs = spectrum.knots[..spectrum.count]
                        .iter()
                        .any(|(_, k)| *k > 0.0);
                    previous_model = Some(model);
                }
                mapped_absorption = model_absorbs;
            }
            hash.update(&density.to_bits().to_le_bytes());
            if let Some(temperatures) = temperature_k {
                let t = temperatures[index];
                if !t.is_finite() || t <= 0.0 {
                    return Err(VolumeTransportError::Temperature);
                }
                hash.update(&t.to_bits().to_le_bytes());
            }
        }
        checkpoint(cx)?;
        let vacuum = if absorption_palette.is_some() {
            !mapped_absorption
        } else {
            !positive_density
                || absorption.knots[..absorption.count]
                    .iter()
                    .all(|(_, k)| *k == 0.0)
        };
        let emission = if !vacuum && let Some(temperatures) = temperature_k {
            emission::EmissionDistribution::try_new(
                cx,
                &grid,
                temperatures,
                absorption.wavelength_bounds_nm(),
            )?
        } else {
            emission::EmissionDistribution::default()
        };
        Ok(Self {
            grid,
            absorption,
            absorption_palette,
            temperature_k,
            max_steps: budget.max_steps_per_segment,
            identity: hash.finalize(),
            vacuum,
            emission,
        })
    }

    /// Numerical snapshot identity for progressive-render binding.
    #[must_use]
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }

    /// Inclusive wavelength domain supported by the accepted absorption law.
    #[must_use]
    pub const fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.absorption.wavelength_bounds_nm()
    }

    /// Cell-resolved optical metadata, when present. Density and temperature
    /// remain the original borrowed fields regardless of palette ownership.
    #[must_use]
    pub const fn absorption_palette(&self) -> Option<&CellPalette<MassAbsorptionSpectrum>> {
        self.absorption_palette.as_ref()
    }

    /// Original mass-absorption law for one native x-major cell index.
    pub fn absorption_at_cell(
        &self,
        cell: usize,
    ) -> Result<&MassAbsorptionSpectrum, VolumeTransportError> {
        if cell >= self.grid.data.len() {
            return Err(VolumeTransportError::Palette);
        }
        Ok(self
            .absorption_palette
            .as_ref()
            .map_or(&self.absorption, |palette| {
                // Palette length and every index were admitted before publication.
                &palette.values()[palette.cell_indices()[cell] as usize]
            }))
    }

    /// Evaluate the cell's mass absorption in m²/kg without ray allocation.
    pub fn coefficients_at_cell<const N: usize>(
        &self,
        cell: usize,
        wavelengths_nm: &[f64; N],
    ) -> Result<[f64; N], VolumeTransportError> {
        let spectrum = self.absorption_at_cell(cell)?;
        let mut result = [0.0; N];
        for (value, &wavelength) in result.iter_mut().zip(wavelengths_nm) {
            *value = spectrum.at(wavelength)?;
        }
        Ok(result)
    }

    /// World-space bounds of the borrowed physical cell support, in metres.
    #[must_use]
    pub fn bounds_m(&self) -> ([f64; 3], [f64; 3]) {
        self.grid.bounds()
    }

    /// Every wavelength has identically zero absorption and emission.
    #[must_use]
    pub const fn is_vacuum(&self) -> bool {
        self.vacuum
    }

    /// A non-vacuum snapshot declares a positive-temperature LTE source.
    /// Very cold visible spectra may still underflow to zero.
    #[must_use]
    pub fn has_thermal_emission(&self) -> bool {
        !self.vacuum && self.temperature_k.is_some()
    }

    /// Integrate a unit-direction segment, clipped to the finite grid bounds.
    /// Positive infinity is accepted for an environment-bound ray.
    pub fn transport<const N: usize>(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        direction: Vec3,
        max_distance_m: f64,
        wavelengths_nm: &[f64; N],
    ) -> Result<SegmentTransport<N>, VolumeTransportError> {
        self.integrate::<N, true>(cx, origin, direction, max_distance_m, wavelengths_nm)
    }

    /// Shadow transmittance; accumulates optical depth before one exp per lane
    /// and never evaluates thermal spectra.
    pub fn transmittance<const N: usize>(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        direction: Vec3,
        max_distance_m: f64,
        wavelengths_nm: &[f64; N],
    ) -> Result<[f64; N], VolumeTransportError> {
        Ok(self
            .integrate::<N, false>(cx, origin, direction, max_distance_m, wavelengths_nm)?
            .transmittance)
    }

    fn integrate<const N: usize, const EMIT: bool>(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        direction: Vec3,
        distance: f64,
        wavelengths: &[f64; N],
    ) -> Result<SegmentTransport<N>, VolumeTransportError> {
        checkpoint(cx)?;
        cells::validate_ray(origin, direction, distance)?;
        let mut kappa = [0.0; N];
        for lane in 0..N {
            kappa[lane] = self.absorption.at(wavelengths[lane])?;
        }
        let mut result = SegmentTransport {
            transmittance: [1.0; N],
            emission_w_m2_sr_nm: [0.0; N],
        };
        if self.vacuum || distance == 0.0 {
            return Ok(result);
        }
        let Some(mut ray) =
            cells::CellRay::clip(&self.grid, origin, direction, distance, self.max_steps)?
        else {
            return Ok(result);
        };
        let mut optical_depth = [0.0; N];
        let mut source_cache = [(f64::NAN, 0.0); N];
        let mut previous_model = None;
        while let Some(segment) = ray.next(cx)? {
            let length = segment.length_m;
            let cell = segment.cell;
            let density = self.grid.data[cell];
            if length > 0.0 && density > 0.0 {
                if let Some(palette) = &self.absorption_palette {
                    let model = palette.cell_indices()[cell];
                    if previous_model != Some(model) {
                        let spectrum = &palette.values()[model as usize];
                        for lane in 0..N {
                            kappa[lane] = spectrum.at(wavelengths[lane])?;
                        }
                        previous_model = Some(model);
                    }
                }
                for lane in 0..N {
                    if kappa[lane] == 0.0 {
                        continue;
                    }
                    let tau = positive_product(density, kappa[lane], length);
                    if EMIT && let Some(temperatures) = self.temperature_k {
                        let weight = det::exp(-optical_depth[lane]) * -det::expm1(-tau);
                        if weight > 0.0 {
                            // Uniform-temperature runs share the expensive
                            // spectral source evaluation without a field copy.
                            if source_cache[lane].0 != temperatures[cell] {
                                source_cache[lane] = (
                                    temperatures[cell],
                                    planck_radiance_w_m2_sr_nm(
                                        temperatures[cell],
                                        wavelengths[lane],
                                    )
                                    .map_err(|_| VolumeTransportError::Unrepresentable)?,
                                );
                            }
                            result.emission_w_m2_sr_nm[lane] += weight * source_cache[lane].1;
                            if !result.emission_w_m2_sr_nm[lane].is_finite() {
                                return Err(VolumeTransportError::Unrepresentable);
                            }
                        }
                    }
                    optical_depth[lane] += tau;
                }
            }
        }
        for lane in 0..N {
            result.transmittance[lane] = det::exp(-optical_depth[lane]);
        }
        Ok(result)
    }
}

fn update_absorption_identity(hash: &mut DomainHasher, spectrum: &MassAbsorptionSpectrum) {
    hash.update(&(spectrum.count as u64).to_le_bytes());
    for &(wavelength, coefficient) in &spectrum.knots[..spectrum.count] {
        hash.update(&wavelength.to_bits().to_le_bytes());
        hash.update(&coefficient.to_bits().to_le_bytes());
    }
    if let Some(source) = spectrum.source_identity() {
        hash.update(b"source-table");
        hash.update(source.as_bytes());
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), VolumeTransportError> {
    cx.checkpoint().map_err(|_| VolumeTransportError::Cancelled)
}

fn positive_product(a: f64, b: f64, c: f64) -> f64 {
    let ab = a * b;
    let product = ab * c;
    if ab >= f64::MIN_POSITIVE && product > 0.0 && product.is_finite() {
        product
    } else {
        det::exp(det::ln(a) + det::ln(b) + det::ln(c))
    }
}

#[cfg(test)]
mod indexed;
#[cfg(test)]
mod tests;
