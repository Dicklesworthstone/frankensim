//! Spatially resolved absorption of a collimated spectral beam.
//!
//! The powers supplied here are watts per spectral lane, already integrated
//! over the beam cross-section and the lane's spectral quadrature weight.
//! They are not radiance or irradiance. A ray represents a narrow beam whose
//! cross-section sees the same piecewise-constant material along the ray.
//!
//! Each cell removes `P * (1 - exp(-rho * (kappa_a + kappa_s) * length))`
//! from the directed beam. Only the absorption share becomes local heat;
//! the scattering share remains radiation. Tracking that scattered radiation
//! requires the separate multiple-scattering transport API. Thermal emission
//! is also a separate source, so this report is deposited incident power,
//! not the material's net radiative heating.

use super::cells::{CellRay, validate_ray};
use super::scattering::ScatteringGrid;
use super::{AbsorbingGrid, VolumeTransportError, checkpoint, positive_product};
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_math::det;

/// Incident-beam power removed in one original simulation cell.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CellBeamPower<const N: usize> {
    /// X-major index in the density buffer borrowed by the medium.
    pub cell_index: usize,
    /// Local heat input, watts per spectral lane.
    pub absorbed_power_w: [f64; N],
    /// Radiation scattered out of the directed beam, watts per spectral lane.
    /// This term must not be deposited as heat.
    pub scattered_power_w: [f64; N],
}

/// Complete power balance and spatial heat input for one finite beam segment.
///
/// Every lane obeys incident = transmitted + absorbed + scattered, to floating
/// point roundoff. Cells follow ray traversal order; a cell with no removed
/// power is omitted. No result is returned on cancellation or budget refusal.
#[derive(Clone, Debug, PartialEq)]
pub struct BeamPowerDeposition<const N: usize> {
    /// Identity of the exact accepted medium and its optical model.
    pub medium_identity: ContentHash,
    /// Inclusive vacuum-wavelength support of the source model, in nm.
    /// This remains the admitted source domain, not the packet's extrema.
    pub wavelength_bounds_nm: (f64, f64),
    /// Vacuum wavelength of each already weighted spectral power lane.
    pub wavelengths_nm: [f64; N],
    /// Power entering the requested segment, watts per spectral lane.
    pub incident_power_w: [f64; N],
    /// Power remaining in the directed beam at the segment endpoint.
    pub transmitted_power_w: [f64; N],
    /// Sum of the cells' local heat input, watts per spectral lane.
    pub absorbed_power_w: [f64; N],
    /// Sum scattered out of the directed beam, watts per spectral lane.
    pub scattered_power_w: [f64; N],
    /// Per-cell power, bounded by the caller's `max_deposition_cells`.
    pub cells: Vec<CellBeamPower<N>>,
}

/// Geometry/density retained by temperature-dependent incident-beam owners.
/// The shared DDA remains the only authority for physical cell crossings.
#[derive(Clone, Copy, Debug)]
#[cfg(feature = "optical-heating")]
pub(crate) struct BeamCellSegment {
    pub cell_index: usize,
    pub density_kg_m3: f64,
    pub length_m: f64,
}

impl AbsorbingGrid<'_> {
    #[cfg(feature = "optical-heating")]
    pub(crate) fn beam_cell_path(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        unit_direction: Vec3,
        max_distance_m: f64,
        max_cells: usize,
    ) -> Result<Vec<BeamCellSegment>, VolumeTransportError> {
        checkpoint(cx)?;
        validate_ray(origin, unit_direction, max_distance_m)?;
        let mut path = Vec::new();
        if let Some(mut cells) = CellRay::clip(
            &self.grid,
            origin,
            unit_direction,
            max_distance_m,
            self.max_steps,
        )? {
            while let Some(cell) = cells.next(cx)? {
                let density = self.grid.data[cell.cell];
                if density == 0.0 || cell.length_m == 0.0 {
                    continue;
                }
                if path.len() == max_cells {
                    return Err(VolumeTransportError::Budget);
                }
                if path.len() == path.capacity() {
                    path.try_reserve_exact(path.capacity().max(16).min(max_cells - path.len()))
                        .map_err(|_| VolumeTransportError::Budget)?;
                }
                path.push(BeamCellSegment {
                    cell_index: cell.cell,
                    density_kg_m3: density,
                    length_m: cell.length_m,
                });
            }
        }
        checkpoint(cx)?;
        Ok(path)
    }

    /// Integrate spatial heat input from an incident beam by exact cell
    /// crossings and Beer–Lambert attenuation.
    ///
    /// `unit_direction` must be unit length; distances are metres. Each
    /// `incident_power_w` lane must be finite and nonnegative, and already
    /// includes its spectral integration weight and beam area. The returned
    /// absorbed watts may be multiplied by an exposure time for deposited
    /// joules, or supplied directly to a thermal solver's cell power sources.
    ///
    /// `max_deposition_cells` bounds the output allocation independently of
    /// the medium's traversal budget. Zero permits an empty deposition. The
    /// borrowed density buffer and the execution context are never mutated.
    #[allow(clippy::too_many_arguments)]
    pub fn deposit_beam_power<const N: usize>(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        unit_direction: Vec3,
        max_distance_m: f64,
        wavelengths_nm: &[f64; N],
        incident_power_w: &[f64; N],
        max_deposition_cells: usize,
    ) -> Result<BeamPowerDeposition<N>, VolumeTransportError> {
        deposit(
            cx,
            self,
            (origin, unit_direction, max_distance_m),
            wavelengths_nm,
            incident_power_w,
            max_deposition_cells,
            None,
        )
    }
}

/// The same stable Beer-law power partition used by cell deposition, with a
/// coefficient supplied by the temperature-dependent constitutive owner.
#[cfg(feature = "optical-heating")]
pub(crate) fn absorb_beam(power: f64, density: f64, length: f64, coefficient: f64) -> (f64, f64) {
    if power == 0.0 || coefficient == 0.0 {
        return (0.0, power);
    }
    let tau = positive_product(density, coefficient, length);
    let absorbed = removed_power(
        power,
        -det::expm1(-tau),
        density,
        length,
        coefficient,
        coefficient,
    );
    let transmission = det::exp(-tau);
    let transmitted = if transmission >= f64::MIN_POSITIVE {
        power * transmission
    } else {
        det::exp(det::ln(power) - tau)
    };
    (absorbed, transmitted)
}

pub(super) fn deposit<const N: usize>(
    cx: &Cx<'_>,
    medium: &AbsorbingGrid<'_>,
    ray: (Point3, Vec3, f64),
    wavelengths_nm: &[f64; N],
    incident_power_w: &[f64; N],
    max_deposition_cells: usize,
    scattering: Option<&ScatteringGrid<'_>>,
) -> Result<BeamPowerDeposition<N>, VolumeTransportError> {
    checkpoint(cx)?;
    validate_ray(ray.0, ray.1, ray.2)?;
    if N == 0 {
        return Err(VolumeTransportError::Sample);
    }
    let indexed = medium.absorption_palette().is_some()
        || scattering.is_some_and(|grid| grid.scattering_palette().is_some());
    let key_at_cell = |cell| {
        if let Some(grid) = scattering {
            grid.coefficient_key(cell)
        } else {
            Ok((
                medium
                    .absorption_palette()
                    .map_or(0, |palette| palette.cell_indices()[cell]),
                0,
            ))
        }
    };
    let coefficients_at_cell = |cell| {
        if let Some(grid) = scattering {
            grid.coefficients_at_cell(cell, wavelengths_nm)
        } else {
            let absorption = medium.coefficients_at_cell(cell, wavelengths_nm)?;
            Ok((absorption, [0.0; N], absorption))
        }
    };
    let mut coefficient_key = key_at_cell(0)?;
    let (mut absorption, mut scatter, mut extinction) = coefficients_at_cell(0)?;
    for lane in 0..N {
        if !incident_power_w[lane].is_finite() || incident_power_w[lane] < 0.0 {
            return Err(VolumeTransportError::BeamPower { lane });
        }
    }
    let ray = if (!indexed && extinction.iter().all(|k| *k == 0.0))
        || scattering.map_or(medium.is_vacuum(), ScatteringGrid::is_vacuum)
    {
        (ray.0, ray.1, 0.0)
    } else {
        ray
    };
    deposit_prepared(
        &mut || checkpoint(cx),
        medium,
        ray,
        wavelengths_nm,
        incident_power_w,
        max_deposition_cells,
        (
            scattering.map_or(medium.identity, ScatteringGrid::content_identity),
            medium.wavelength_bounds_nm(),
        ),
        |cell| {
            if indexed {
                let key = key_at_cell(cell)?;
                if key != coefficient_key {
                    (absorption, scatter, extinction) = coefficients_at_cell(cell)?;
                    coefficient_key = key;
                }
            }
            Ok((absorption, scatter, extinction))
        },
    )
}

/// Native exact spatial beam integration for either uniform or borrowed cell
/// coefficients. The original DDA and stable optical-depth/power arithmetic
/// are the sole implementation for both owners.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
pub(super) fn deposit_prepared<const N: usize>(
    check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
    medium: &AbsorbingGrid<'_>,
    ray: (Point3, Vec3, f64),
    wavelengths_nm: &[f64; N],
    incident_power_w: &[f64; N],
    max_deposition_cells: usize,
    metadata: (ContentHash, (f64, f64)),
    mut coefficients: impl FnMut(usize) -> Result<([f64; N], [f64; N], [f64; N]), VolumeTransportError>,
) -> Result<BeamPowerDeposition<N>, VolumeTransportError> {
    check()?;
    validate_ray(ray.0, ray.1, ray.2)?;
    for (lane, power) in incident_power_w.iter().enumerate() {
        if !power.is_finite() || *power < 0.0 {
            return Err(VolumeTransportError::BeamPower { lane });
        }
    }
    let (identity, bounds) = metadata;
    let mut result = BeamPowerDeposition {
        medium_identity: identity,
        wavelength_bounds_nm: bounds,
        wavelengths_nm: *wavelengths_nm,
        incident_power_w: *incident_power_w,
        transmitted_power_w: *incident_power_w,
        absorbed_power_w: [0.0; N],
        scattered_power_w: [0.0; N],
        cells: Vec::new(),
    };
    if incident_power_w.iter().all(|p| *p == 0.0) {
        return Ok(result);
    }
    let Some(mut cells) = CellRay::clip(&medium.grid, ray.0, ray.1, ray.2, medium.max_steps)?
    else {
        return Ok(result);
    };
    let mut absorbed_correction = [0.0; N];
    let mut scattered_correction = [0.0; N];
    while let Some(cell) = cells.next_with_checkpoint(check)? {
        let rho = medium.grid.data[cell.cell];
        if rho == 0.0 || cell.length_m == 0.0 {
            continue;
        }
        let (absorption, scatter, extinction) = coefficients(cell.cell)?;
        let mut deposit = CellBeamPower {
            cell_index: cell.cell,
            absorbed_power_w: [0.0; N],
            scattered_power_w: [0.0; N],
        };
        for lane in 0..N {
            let power = result.transmitted_power_w[lane];
            let total = extinction[lane];
            if power == 0.0 || total == 0.0 {
                continue;
            }
            let tau = positive_product(rho, total, cell.length_m);
            let opacity = -det::expm1(-tau);
            deposit.absorbed_power_w[lane] =
                removed_power(power, opacity, rho, cell.length_m, absorption[lane], total);
            deposit.scattered_power_w[lane] =
                removed_power(power, opacity, rho, cell.length_m, scatter[lane], total);
            let transmittance = det::exp(-tau);
            // A tiny transmission fraction can still carry representable
            // watts when the incoming beam is large.
            result.transmitted_power_w[lane] = if transmittance >= f64::MIN_POSITIVE {
                power * transmittance
            } else {
                det::exp(det::ln(power) - tau)
            };
        }
        if deposit.absorbed_power_w.iter().all(|p| *p == 0.0)
            && deposit.scattered_power_w.iter().all(|p| *p == 0.0)
        {
            continue;
        }
        if result.cells.len() == max_deposition_cells {
            return Err(VolumeTransportError::Budget);
        }
        if result.cells.len() == result.cells.capacity() {
            let additional = result
                .cells
                .capacity()
                .max(16)
                .min(max_deposition_cells - result.cells.len());
            result
                .cells
                .try_reserve_exact(additional)
                .map_err(|_| VolumeTransportError::Budget)?;
        }
        for lane in 0..N {
            add_compensated(
                &mut result.absorbed_power_w[lane],
                &mut absorbed_correction[lane],
                deposit.absorbed_power_w[lane],
            )?;
            add_compensated(
                &mut result.scattered_power_w[lane],
                &mut scattered_correction[lane],
                deposit.scattered_power_w[lane],
            )?;
        }
        result.cells.push(deposit);
    }
    check()?;
    Ok(result)
}

/// Absorption or scattering share of removed power. The ordinary path is
/// short; the log fallback preserves a representable heat input when either
/// the opacity or the coefficient ratio underflows before power is applied.
pub(crate) fn removed_power(
    power: f64,
    opacity: f64,
    density: f64,
    length: f64,
    coefficient: f64,
    extinction: f64,
) -> f64 {
    if coefficient == 0.0 {
        return 0.0;
    }
    if opacity == 0.0 {
        // Optical depth is below the representable range, where 1-exp(-tau)
        // equals tau to far better than f64 precision.
        return det::exp(
            det::ln(power) + det::ln(density) + det::ln(length) + det::ln(coefficient),
        );
    }
    let removed = power * opacity;
    let share = coefficient / extinction;
    let result = removed * share;
    if removed >= f64::MIN_POSITIVE && share >= f64::MIN_POSITIVE && result > 0.0 {
        result
    } else {
        det::exp(det::ln(power) + det::ln(opacity) + det::ln(coefficient) - det::ln(extinction))
    }
}

fn add_compensated(
    sum: &mut f64,
    correction: &mut f64,
    value: f64,
) -> Result<(), VolumeTransportError> {
    let increment = value - *correction;
    let next = *sum + increment;
    if !next.is_finite() {
        return Err(VolumeTransportError::Unrepresentable);
    }
    *correction = (next - *sum) - increment;
    *sum = next;
    Ok(())
}

#[cfg(test)]
mod tests;
