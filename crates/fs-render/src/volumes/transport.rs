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

const MAX_KNOTS: usize = 16;

/// Validated visible-band mass absorption coefficient, in m²/kg.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MassAbsorptionSpectrum {
    knots: [(f64, f64); MAX_KNOTS],
    count: usize,
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
            || knots
                .iter()
                .any(|(l, k)| !l.is_finite() || !k.is_finite() || *k < 0.0)
            || knots.windows(2).any(|pair| pair[0].0 >= pair[1].0)
        {
            return Err(VolumeTransportError::AbsorptionSpectrum);
        }
        let mut result = Self {
            knots: [(0.0, 0.0); MAX_KNOTS],
            count: knots.len(),
        };
        result.knots[..knots.len()].copy_from_slice(knots);
        Ok(result)
    }

    /// Evaluate without extrapolating outside the admitted visible band.
    pub fn at(&self, wavelength_nm: f64) -> Result<f64, VolumeTransportError> {
        if !wavelength_nm.is_finite() || !(LAMBDA_MIN..=LAMBDA_MAX).contains(&wavelength_nm) {
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
    /// Maximum number of borrowed cells scanned at admission.
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
    /// Wavelength lies outside 380–780 nm.
    Wavelength,
    /// Ray direction is not finite and unit, or distance is negative/NaN.
    Ray,
    /// Admission scan or per-segment work exceeds the caller's limit.
    Budget,
    /// Finite arithmetic cannot resolve a cell crossing or physical radiance.
    Unrepresentable,
}

impl core::fmt::Display for VolumeTransportError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "absorbing-volume transport refused: {self:?}")
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
    temperature_k: Option<&'a [f64]>,
    max_steps: u32,
    identity: ContentHash,
    vacuum: bool,
}

impl<'a> AbsorbingGrid<'a> {
    /// Validate, hash and borrow the original simulation buffers once.
    pub fn try_new(
        cx: &Cx<'_>,
        grid: VolumeGrid<'a>,
        absorption: MassAbsorptionSpectrum,
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
        hash.update(&(absorption.count as u64).to_le_bytes());
        for &(l, k) in &absorption.knots[..absorption.count] {
            hash.update(&l.to_bits().to_le_bytes());
            hash.update(&k.to_bits().to_le_bytes());
        }
        hash.update(&[u8::from(temperature_k.is_some())]);
        let mut positive_density = false;
        for (index, &density) in grid.data.iter().enumerate() {
            if index % 1024 == 0 {
                checkpoint(cx)?;
            }
            if !density.is_finite() || density < 0.0 {
                return Err(VolumeTransportError::Density { index });
            }
            positive_density |= density > 0.0;
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
        let vacuum = !positive_density
            || absorption.knots[..absorption.count]
                .iter()
                .all(|(_, k)| *k == 0.0);
        Ok(Self {
            grid,
            absorption,
            temperature_k,
            max_steps: budget.max_steps_per_segment,
            identity: hash.finalize(),
            vacuum,
        })
    }

    /// Numerical snapshot identity for progressive-render binding.
    #[must_use]
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
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
        let o = [origin.x, origin.y, origin.z];
        let d = [direction.x, direction.y, direction.z];
        if o.iter().chain(d.iter()).any(|v| !v.is_finite())
            || (direction.dot(direction) - 1.0).abs() > 2.0e-10
            || distance.is_nan()
            || distance < 0.0
        {
            return Err(VolumeTransportError::Ray);
        }
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
        let (lo, hi) = self.grid.bounds();
        let (mut enter, mut exit) = (0.0_f64, distance);
        let mut collapsed_interval = false;
        for a in 0..3 {
            if d[a] == 0.0 {
                if o[a] < lo[a] || o[a] >= hi[a] {
                    return Ok(result);
                }
            } else {
                let t1 = (lo[a] - o[a]) / d[a];
                let t2 = (hi[a] - o[a]) / d[a];
                collapsed_interval |= t1 == t2;
                enter = enter.max(t1.min(t2));
                exit = exit.min(t1.max(t2));
            }
        }
        if exit < enter {
            return Ok(result);
        }
        if exit == enter {
            return if collapsed_interval {
                Err(VolumeTransportError::Unrepresentable)
            } else {
                Ok(result)
            };
        }
        if !enter.is_finite() || !exit.is_finite() {
            return Err(VolumeTransportError::Unrepresentable);
        }
        // Rebase at entry so far-away camera parameters do not accumulate in
        // each cell length. Recompute crossings from integer cell boundaries.
        let start =
            std::array::from_fn::<_, 3, _>(|a| d[a].mul_add(enter, o[a]).clamp(lo[a], hi[a]));
        let mut end = exit - enter;
        for a in 0..3 {
            if d[a] != 0.0 {
                let far = if d[a] > 0.0 { hi[a] } else { lo[a] };
                end = end.min((far - start[a]) / d[a]);
            }
        }
        let mut index = [0usize; 3];
        for a in 0..3 {
            let q = (start[a] - lo[a]) / self.grid.cell[a];
            let floor = q.floor();
            let cell = if d[a] < 0.0 && q == floor {
                floor - 1.0
            } else {
                floor
            };
            index[a] = cell.clamp(0.0, (self.grid.dims[a] - 1) as f64) as usize;
        }
        let crossing = |a: usize, i: usize| {
            if d[a] == 0.0 {
                return f64::INFINITY;
            }
            let boundary = if d[a] > 0.0 { i + 1 } else { i };
            let coordinate = (boundary as f64).mul_add(self.grid.cell[a], lo[a]);
            (coordinate - start[a]) / d[a]
        };
        let mut next = std::array::from_fn::<_, 3, _>(|a| crossing(a, index[a]));
        let mut current = 0.0;
        let mut optical_depth = [0.0; N];
        let mut source_cache = [(f64::NAN, 0.0); N];
        for _ in 0..self.max_steps {
            checkpoint(cx)?;
            let stop = next.into_iter().fold(end, f64::min).max(current);
            let length = stop - current;
            let cell = index[0] + self.grid.dims[0] * (index[1] + self.grid.dims[1] * index[2]);
            let density = self.grid.data[cell];
            if length > 0.0 && density > 0.0 {
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
            current = stop;
            if current >= end {
                for lane in 0..N {
                    result.transmittance[lane] = det::exp(-optical_depth[lane]);
                }
                return Ok(result);
            }
            for a in 0..3 {
                if next[a] <= stop {
                    if (d[a] < 0.0 && index[a] == 0)
                        || (d[a] > 0.0 && index[a] + 1 == self.grid.dims[a])
                    {
                        return Err(VolumeTransportError::Unrepresentable);
                    }
                    if d[a] > 0.0 {
                        index[a] += 1;
                    } else {
                        index[a] -= 1;
                    }
                    next[a] = crossing(a, index[a]);
                }
            }
        }
        Err(VolumeTransportError::Budget)
    }
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), VolumeTransportError> {
    cx.checkpoint().map_err(|_| VolumeTransportError::Cancelled)
}

fn positive_product(a: f64, b: f64, c: f64) -> f64 {
    let product = a * b * c;
    if product > 0.0 && product.is_finite() {
        product
    } else {
        det::exp(det::ln(a) + det::ln(b) + det::ln(c))
    }
}

#[cfg(test)]
mod tests;
