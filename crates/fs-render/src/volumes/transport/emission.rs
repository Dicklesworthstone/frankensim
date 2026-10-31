//! Direction proposals for physical, piecewise-constant thermal volume sources.
//!
//! A source point proposes only a direction. Its radial coordinate is then
//! marginalized: `q(omega) = integral p_volume(x + r*omega) r² dr`. Evaluating
//! the complete source line integral on that direction avoids the singular
//! point-source `1/r²` estimator for receivers inside an emitting cell.
//! Occluders clip radiance, never this proposal density.

use super::{
    AbsorbingGrid, CellPalette, ThermalEmissionModel, VolumeGrid, VolumeTransportError, cells,
    checkpoint, planck_radiance_w_m2_sr_nm,
};
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_math::det;

/// Source-domain weighting, support mixture and point-to-direction map.
/// Version 2 retains the declared wavelength bounds; visible models keep their
/// original proposal arithmetic and sampling distribution.
pub const EMISSION_PROPOSAL_SEMANTICS_VERSION: u32 = 2;

const UNIFORM_FRACTION: f64 = 1.0 / 16.0;
// The wavelength maximum of spectral radiance per wavelength, in nm K.
const WIEN_NM_K: f64 = 2.897_771_955_185_172_7e6;

#[derive(Default)]
pub(super) struct EmissionDistribution {
    // Actual intervals are authoritative for BOTH inversion and density. The
    // field itself remains borrowed; these are only sampling probabilities.
    cdf: Vec<f64>,
    inverse_cell_volume: f64,
    log_cell_volume: f64,
}

impl EmissionDistribution {
    pub(super) fn try_new(
        cx: &Cx<'_>,
        grid: &VolumeGrid<'_>,
        temperatures: &[f64],
        wavelength_bounds_nm: (f64, f64),
    ) -> Result<Self, VolumeTransportError> {
        let mut cache = SourceCache::default();
        let mut count = 0usize;
        let mut maximum = f64::NEG_INFINITY;
        for (index, (&density, &temperature)) in grid.data.iter().zip(temperatures).enumerate() {
            if index % 1024 == 0 {
                checkpoint(cx)?;
            }
            if density > 0.0 {
                count += 1;
                maximum = maximum
                    .max(det::ln(density) + cache.log_source(temperature, wavelength_bounds_nm)?);
            }
        }
        checkpoint(cx)?;
        // The clamped Wien maximum covers the COMPLETE admitted source band.
        // An empty proposal cannot discard a representable source in that band.
        if maximum == f64::NEG_INFINITY {
            return Ok(Self::default());
        }
        let length = grid
            .data
            .len()
            .checked_add(1)
            .ok_or(VolumeTransportError::Budget)?;
        let mut cdf = Vec::new();
        cdf.try_reserve_exact(length)
            .map_err(|_| VolumeTransportError::Budget)?;
        cdf.push(0.0);
        let mut total = 0.0;
        for (index, (&density, &temperature)) in grid.data.iter().zip(temperatures).enumerate() {
            if index % 1024 == 0 {
                checkpoint(cx)?;
            }
            let weight = if density > 0.0 {
                det::exp(
                    det::ln(density) + cache.log_source(temperature, wavelength_bounds_nm)?
                        - maximum,
                )
            } else {
                0.0
            };
            total += weight;
            cdf.push(weight);
        }
        let uniform = UNIFORM_FRACTION / count as f64;
        let mut sum = 0.0;
        for (index, &density) in grid.data.iter().enumerate() {
            if index % 1024 == 0 {
                checkpoint(cx)?;
            }
            let probability = if density > 0.0 {
                uniform + (1.0 - UNIFORM_FRACTION) * (cdf[index + 1] / total)
            } else {
                0.0
            };
            let next = sum + probability;
            if density > 0.0 && next <= sum {
                return Err(VolumeTransportError::Unrepresentable);
            }
            sum = next;
            cdf[index + 1] = sum;
        }
        // Normalize accumulated roundoff, then verify that no required support
        // interval vanished. PDF evaluation uses these resulting differences.
        for (index, &density) in grid.data.iter().enumerate() {
            if index % 1024 == 0 {
                checkpoint(cx)?;
            }
            cdf[index + 1] /= sum;
            if density > 0.0 && cdf[index + 1] <= cdf[index] {
                return Err(VolumeTransportError::Unrepresentable);
            }
        }
        let volume = grid.cell[0] * grid.cell[1] * grid.cell[2];
        let inverse = 1.0 / volume;
        checkpoint(cx)?;
        Ok(Self {
            cdf,
            inverse_cell_volume: if inverse.is_finite() && inverse > 0.0 {
                inverse
            } else {
                0.0
            },
            log_cell_volume: grid.cell.iter().map(|&length| det::ln(length)).sum(),
        })
    }

    /// Prepare the same directional proposal for constituent-temperature
    /// sources, without copying a temperature field into physical cells.
    /// The largest active temperature bounds Planck radiance at every admitted
    /// wavelength; it affects importance weights only, never source evaluation.
    pub(super) fn try_new_sources(
        cx: &Cx<'_>,
        grid: &VolumeGrid<'_>,
        sources: &CellPalette<ThermalEmissionModel>,
    ) -> Result<Self, VolumeTransportError> {
        checkpoint(cx)?;
        if sources.cell_indices().len() != grid.data.len() {
            return Err(VolumeTransportError::Palette);
        }
        let bounds = sources.values()[0].wavelength_bounds_nm();
        let mut cache = SourceCache::default();
        let mut models = Vec::new();
        models
            .try_reserve_exact(sources.values().len())
            .map_err(|_| VolumeTransportError::Budget)?;
        for model in sources.values() {
            checkpoint(cx)?;
            if model.wavelength_bounds_nm() != bounds {
                return Err(VolumeTransportError::Wavelength);
            }
            models.push(None);
        }
        let mut count = 0usize;
        let mut maximum = f64::NEG_INFINITY;
        for (cell, &density) in grid.data.iter().enumerate() {
            if cell % 1024 == 0 {
                checkpoint(cx)?;
            }
            if density == 0.0 {
                continue;
            }
            let index = sources.cell_indices()[cell] as usize;
            // Planck maxima are evaluated once per OCCUPIED optical model.
            // Unreferenced models and zero-density cells cannot create a source.
            // Retain an activity bit even when a cold importance weight vanishes.
            let (active, log_weight) = if let Some(weight) = models[index] {
                weight
            } else {
                let model = &sources.values()[index];
                let absorption = model.maximum_absorption_m2_kg();
                let weight = if let Some(temperature) = model.maximum_temperature_k()
                    && absorption > 0.0
                {
                    (
                        true,
                        det::ln(absorption) + cache.log_source(temperature, bounds)?,
                    )
                } else {
                    (false, f64::NEG_INFINITY)
                };
                models[index] = Some(weight);
                weight
            };
            if active {
                count += 1;
                maximum = maximum.max(det::ln(density) + log_weight);
            }
        }
        checkpoint(cx)?;
        if maximum == f64::NEG_INFINITY {
            return Ok(Self::default());
        }
        let active = |cell: usize, density: f64| {
            density > 0.0
                && models[sources.cell_indices()[cell] as usize].is_some_and(|weight| weight.0)
        };
        let length = grid
            .data
            .len()
            .checked_add(1)
            .ok_or(VolumeTransportError::Budget)?;
        let mut cdf = Vec::new();
        cdf.try_reserve_exact(length)
            .map_err(|_| VolumeTransportError::Budget)?;
        cdf.push(0.0);
        let mut total = 0.0;
        for (cell, &density) in grid.data.iter().enumerate() {
            if cell % 1024 == 0 {
                checkpoint(cx)?;
            }
            let weight = if active(cell, density) {
                det::exp(
                    det::ln(density)
                        + models[sources.cell_indices()[cell] as usize]
                            .expect("occupied source model")
                            .1
                        - maximum,
                )
            } else {
                0.0
            };
            total += weight;
            cdf.push(weight);
        }
        let uniform = UNIFORM_FRACTION / count as f64;
        let mut sum = 0.0;
        for (cell, &density) in grid.data.iter().enumerate() {
            if cell % 1024 == 0 {
                checkpoint(cx)?;
            }
            let probability = if active(cell, density) {
                uniform + (1.0 - UNIFORM_FRACTION) * (cdf[cell + 1] / total)
            } else {
                0.0
            };
            let next = sum + probability;
            if active(cell, density) && next <= sum {
                return Err(VolumeTransportError::Unrepresentable);
            }
            sum = next;
            cdf[cell + 1] = sum;
        }
        for (cell, &density) in grid.data.iter().enumerate() {
            if cell % 1024 == 0 {
                checkpoint(cx)?;
            }
            cdf[cell + 1] /= sum;
            if active(cell, density) && cdf[cell + 1] <= cdf[cell] {
                return Err(VolumeTransportError::Unrepresentable);
            }
        }
        let volume = grid.cell[0] * grid.cell[1] * grid.cell[2];
        let inverse = 1.0 / volume;
        checkpoint(cx)?;
        Ok(Self {
            cdf,
            inverse_cell_volume: if inverse.is_finite() && inverse > 0.0 {
                inverse
            } else {
                0.0
            },
            log_cell_volume: grid.cell.iter().map(|&length| det::ln(length)).sum(),
        })
    }

    fn probability(&self, cell: usize) -> f64 {
        self.cdf[cell + 1] - self.cdf[cell]
    }

    fn radial_density(
        &self,
        probability: f64,
        start: f64,
        length: f64,
    ) -> Result<f64, VolumeTransportError> {
        // Factored cube difference avoids distant-camera cancellation.
        let moment = length * (start * start + start * length + length * length / 3.0);
        let geometric = moment * self.inverse_cell_volume;
        let ordinary = probability * geometric;
        if moment >= f64::MIN_POSITIVE
            && geometric >= f64::MIN_POSITIVE
            && ordinary >= f64::MIN_POSITIVE
            && ordinary.is_finite()
        {
            return Ok(ordinary);
        }
        // Extreme metric scales can overflow r² or cell volume even when the
        // final angular density is ordinary. Scale before squaring and divide
        // in log space only on that exceptional path.
        let scale = start.max(length);
        let a = start / scale;
        let b = length / scale;
        let log_moment =
            det::ln(length) + 2.0 * det::ln(scale) + det::ln(a * a + a * b + b * b / 3.0);
        let value = det::exp(det::ln(probability) + log_moment - self.log_cell_volume);
        if value.is_finite() && value > 0.0 {
            Ok(value)
        } else {
            Err(VolumeTransportError::Unrepresentable)
        }
    }
}

struct SourceCache {
    values: [(f64, f64); 8],
    next: usize,
}

impl Default for SourceCache {
    fn default() -> Self {
        Self {
            values: [(f64::NAN, 0.0); 8],
            next: 0,
        }
    }
}

impl SourceCache {
    fn log_source(
        &mut self,
        temperature: f64,
        wavelength_bounds_nm: (f64, f64),
    ) -> Result<f64, VolumeTransportError> {
        if let Some((_, value)) = self.values.iter().find(|(t, _)| *t == temperature) {
            return Ok(*value);
        }
        let wavelength =
            (WIEN_NM_K / temperature).clamp(wavelength_bounds_nm.0, wavelength_bounds_nm.1);
        let source = planck_radiance_w_m2_sr_nm(temperature, wavelength)
            .map_err(|_| VolumeTransportError::Unrepresentable)?;
        let value = if source > 0.0 {
            det::ln(source)
        } else {
            f64::NEG_INFINITY
        };
        self.values[self.next] = (temperature, value);
        self.next = (self.next + 1) % self.values.len();
        Ok(value)
    }
}

impl AbsorbingGrid<'_> {
    /// A cached, support-complete proposal exists for source-domain thermal emission.
    /// Vacuum, absent sources, or complete source-band underflow retains
    /// the ordinary transport path without proposal storage.
    #[must_use]
    pub fn has_emission_proposal(&self) -> bool {
        !self.emission.cdf.is_empty()
    }

    /// Select a source cell and a uniform point using four independent `[0,1)`
    /// variates, then retain only its direction from the physical receiver.
    /// This allocates nothing and uses logarithmic-time cell selection.
    /// A source point coinciding exactly with the receiver returns `None`.
    pub fn sample_emission_direction(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        sample: [f64; 4],
    ) -> Result<Option<Vec3>, VolumeTransportError> {
        checkpoint(cx)?;
        if sample
            .iter()
            .any(|value| !value.is_finite() || !(0.0..1.0).contains(value))
        {
            return Err(VolumeTransportError::Sample);
        }
        if [origin.x, origin.y, origin.z]
            .iter()
            .any(|value| !value.is_finite())
        {
            return Err(VolumeTransportError::Ray);
        }
        if !self.has_emission_proposal() {
            return Ok(None);
        }
        let cell = self.emission.cdf.partition_point(|&cdf| cdf <= sample[0]) - 1;
        let [nx, ny, _] = self.grid.dims;
        let index = [cell % nx, (cell / nx) % ny, cell / (nx * ny)];
        let position = std::array::from_fn::<_, 3, _>(|axis| {
            (index[axis] as f64 + sample[axis + 1])
                .mul_add(self.grid.cell[axis], self.grid.origin[axis])
        });
        let displacement = Point3::new(position[0], position[1], position[2]).delta_from(origin);
        let scale = displacement
            .x
            .abs()
            .max(displacement.y.abs())
            .max(displacement.z.abs());
        if !scale.is_finite() {
            return Err(VolumeTransportError::Unrepresentable);
        }
        if scale == 0.0 {
            return Ok(None);
        }
        let scaled = Vec3::new(
            displacement.x / scale,
            displacement.y / scale,
            displacement.z / scale,
        );
        Ok(Some(scaled.scale(1.0 / det::sqrt(scaled.dot(scaled)))))
    }

    /// Full proposal density in sr⁻¹, marginalizing every positive ray/cell
    /// span. It intentionally includes source cells behind geometric blockers:
    /// visibility changes the sampled integrand, never the distribution.
    /// Uses the same exact DDA, cancellation and segment-step limit as transport.
    pub fn emission_direction_pdf(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        direction: Vec3,
    ) -> Result<f64, VolumeTransportError> {
        checkpoint(cx)?;
        cells::validate_ray(origin, direction, f64::INFINITY)?;
        if !self.has_emission_proposal() {
            return Ok(0.0);
        }
        let Some(mut ray) =
            cells::CellRay::clip(&self.grid, origin, direction, f64::INFINITY, self.max_steps)?
        else {
            return Ok(0.0);
        };
        let mut pdf = 0.0;
        while let Some(segment) = ray.next(cx)? {
            let probability = self.emission.probability(segment.cell);
            if probability > 0.0 && segment.length_m > 0.0 {
                pdf += self.emission.radial_density(
                    probability,
                    ray.distance_from_origin(segment.start_m),
                    segment.length_m,
                )?;
                if !pdf.is_finite() {
                    return Err(VolumeTransportError::Unrepresentable);
                }
            }
        }
        Ok(pdf)
    }
}

#[cfg(test)]
mod tests;
