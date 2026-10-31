//! Shared, monomorphized transport over prepared optical coefficients. Field
//! and palette packet views avoid spectrum evaluation at cell crossings; the
//! allocation-free renderer view caches the last native material pair.

use super::super::{AbsorbingGrid, SegmentTransport, cells::CellRay, positive_product};
use super::{
    FreeFlightSample, PhaseFunction, ScatteringEvent, ScatteringGrid, VolumeTransportError,
    collision_weights, divide_product,
};
use crate::thermal_surface::planck_radiance_w_m2_sr_nm;
use fs_blake3::ContentHash;
use fs_geom::{Point3, Vec3};
use fs_math::det;

pub(super) trait PreparedOptics<const N: usize> {
    fn geometry(&self) -> &AbsorbingGrid<'_>;
    fn identity(&self) -> ContentHash;
    fn bounds(&self) -> (f64, f64);
    fn wavelengths(&self) -> &[f64; N];
    fn coefficient(
        &self,
        cell: usize,
        lane: usize,
    ) -> Result<(f64, f64, f64), VolumeTransportError>;
    fn coefficients(
        &self,
        cell: usize,
    ) -> Result<([f64; N], [f64; N], [f64; N]), VolumeTransportError> {
        let mut a = [0.0; N];
        let mut s = [0.0; N];
        let mut t = [0.0; N];
        for lane in 0..N {
            (a[lane], s[lane], t[lane]) = self.coefficient(cell, lane)?;
        }
        Ok((a, s, t))
    }
    fn phase(&self, cell: usize) -> &PhaseFunction;
    fn has_scattering(&self) -> bool;
    fn lane_has_scattering(&self, lane: usize) -> bool;
    fn lane_has_extinction(&self, lane: usize) -> bool;
    fn has_emission(&self) -> bool;
}

pub(super) struct UniformOptics<'g, 'a, const N: usize> {
    pub medium: &'g ScatteringGrid<'a>,
    pub wavelengths: [f64; N],
    pub absorption: [f64; N],
    pub scattering: [f64; N],
    pub extinction: [f64; N],
}
impl<const N: usize> PreparedOptics<N> for UniformOptics<'_, '_, N> {
    fn geometry(&self) -> &AbsorbingGrid<'_> {
        &self.medium.absorption
    }
    fn identity(&self) -> ContentHash {
        self.medium.identity
    }
    fn bounds(&self) -> (f64, f64) {
        self.medium.wavelength_bounds_nm()
    }
    fn wavelengths(&self) -> &[f64; N] {
        &self.wavelengths
    }
    fn coefficient(&self, _: usize, lane: usize) -> Result<(f64, f64, f64), VolumeTransportError> {
        Ok((
            self.absorption[lane],
            self.scattering[lane],
            self.extinction[lane],
        ))
    }
    fn coefficients(
        &self,
        _: usize,
    ) -> Result<([f64; N], [f64; N], [f64; N]), VolumeTransportError> {
        Ok((self.absorption, self.scattering, self.extinction))
    }
    fn phase(&self, _: usize) -> &PhaseFunction {
        &self.medium.phase
    }
    fn has_scattering(&self) -> bool {
        self.medium.has_scattering
    }
    fn lane_has_scattering(&self, lane: usize) -> bool {
        self.medium.has_scattering && self.scattering[lane] > 0.0
    }
    fn lane_has_extinction(&self, lane: usize) -> bool {
        self.extinction[lane] > 0.0
            && (!self.medium.absorption.is_vacuum() || self.medium.has_scattering)
    }
    fn has_emission(&self) -> bool {
        self.medium.has_thermal_emission()
    }
}

/// A lane view reads only that lane's coefficients, without copying the packet.
pub(super) struct Monochromatic<'m, M, const N: usize> {
    pub medium: &'m M,
    pub lane: usize,
    pub wavelength: [f64; 1],
}
impl<M: PreparedOptics<N>, const N: usize> PreparedOptics<1> for Monochromatic<'_, M, N> {
    fn geometry(&self) -> &AbsorbingGrid<'_> {
        self.medium.geometry()
    }
    fn identity(&self) -> ContentHash {
        self.medium.identity()
    }
    fn bounds(&self) -> (f64, f64) {
        self.medium.bounds()
    }
    fn wavelengths(&self) -> &[f64; 1] {
        &self.wavelength
    }
    fn coefficient(&self, cell: usize, _: usize) -> Result<(f64, f64, f64), VolumeTransportError> {
        self.medium.coefficient(cell, self.lane)
    }
    fn phase(&self, cell: usize) -> &PhaseFunction {
        self.medium.phase(cell)
    }
    fn has_scattering(&self) -> bool {
        self.medium.lane_has_scattering(self.lane)
    }
    fn lane_has_scattering(&self, _: usize) -> bool {
        self.has_scattering()
    }
    fn lane_has_extinction(&self, _: usize) -> bool {
        self.medium.lane_has_extinction(self.lane)
    }
    fn has_emission(&self) -> bool {
        self.medium.has_emission()
    }
}

#[allow(clippy::too_many_arguments)]
pub(super) fn sample<const N: usize>(
    medium: &impl PreparedOptics<N>,
    check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
    origin: Point3,
    direction: Vec3,
    max_distance_m: f64,
    u_lane: f64,
    u_distance: f64,
    include_emission: bool,
) -> Result<FreeFlightSample<N>, VolumeTransportError> {
    check()?;
    let wavelengths_nm = medium.wavelengths();
    if !medium.has_scattering() {
        let result = transport(
            medium,
            check,
            origin,
            direction,
            max_distance_m,
            include_emission,
        )?;
        return Ok(FreeFlightSample {
            event: None,
            weight: result.transmittance,
            emission_w_m2_sr_nm: result.emission_w_m2_sr_nm,
            escape_probability: 1.0,
        });
    }
    let mut result = FreeFlightSample {
        event: None,
        weight: [1.0; N],
        emission_w_m2_sr_nm: [0.0; N],
        escape_probability: 1.0,
    };
    let Some(mut ray) = CellRay::clip(
        &medium.geometry().grid,
        origin,
        direction,
        max_distance_m,
        medium.geometry().max_steps,
    )?
    else {
        return Ok(result);
    };
    let hero = ((u_lane * N as f64) as usize).min(N - 1);
    let mut remaining_tau = -det::ln(u_distance.max(f64::MIN_POSITIVE));
    let mut depth = [0.0; N];
    let mut source_cache = [(f64::NAN, 0.0); N];
    while let Some(segment) = ray.next_with_checkpoint(check)? {
        let rho = medium.geometry().grid.data[segment.cell];
        if rho == 0.0 || segment.length_m == 0.0 {
            continue;
        }
        let (absorption, scattering, extinction) = medium.coefficients(segment.cell)?;
        let tau = extinction.map(|k| {
            if k == 0.0 {
                0.0
            } else {
                positive_product(rho, k, segment.length_m)
            }
        });
        if result.event.is_none() && extinction[hero] > 0.0 && remaining_tau < tau[hero] {
            let local_m = divide_product(remaining_tau, rho, extinction[hero]);
            let parameter = segment.start_m + local_m;
            if !parameter.is_finite()
                || parameter < segment.start_m
                || parameter >= segment.start_m + segment.length_m
                || (local_m > 0.0 && parameter == segment.start_m)
                || (remaining_tau > 0.0 && local_m == 0.0)
            {
                return Err(VolumeTransportError::Unrepresentable);
            }
            let event_depth = std::array::from_fn(|lane| {
                depth[lane]
                    + if extinction[lane] == 0.0 || local_m == 0.0 {
                        0.0
                    } else {
                        positive_product(rho, extinction[lane], local_m)
                    }
            });
            result.weight = collision_weights(&event_depth, &extinction, &scattering)?;
            let position = ray.position(parameter);
            let distance_m = ray.distance_from_origin(parameter);
            if !distance_m.is_finite() || (distance_m > 0.0 && position == origin) {
                return Err(VolumeTransportError::Unrepresentable);
            }
            result.event = Some(ScatteringEvent {
                cell_index: segment.cell,
                position,
                distance_m,
            });
            result.escape_probability = 0.0;
            if !include_emission || !medium.has_emission() {
                return Ok(result);
            }
        }
        if result.event.is_none() {
            remaining_tau -= tau[hero];
        }
        if include_emission {
            integrate_source_cell(
                medium.geometry(),
                segment.cell,
                &tau,
                &absorption,
                &extinction,
                wavelengths_nm,
                &mut depth,
                &mut source_cache,
                &mut result.emission_w_m2_sr_nm,
            )?;
        } else {
            for lane in 0..N {
                depth[lane] += tau[lane];
            }
        }
    }
    if result.event.is_none() {
        let minimum = depth.into_iter().fold(f64::INFINITY, f64::min);
        if !minimum.is_finite() {
            return Err(VolumeTransportError::Unrepresentable);
        }
        let relative = depth.map(|tau| det::exp(minimum - tau));
        let sum: f64 = relative.iter().sum();
        result.weight = relative.map(|t| t * (N as f64 / sum));
        result.escape_probability = depth.iter().map(|tau| det::exp(-tau)).sum::<f64>() / N as f64;
    }
    Ok(result)
}

/// Deterministic attenuation/LTE shares the same prepared coefficients and DDA.
pub(super) fn transport<const N: usize>(
    medium: &impl PreparedOptics<N>,
    check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
    origin: Point3,
    direction: Vec3,
    max_distance_m: f64,
    include_emission: bool,
) -> Result<SegmentTransport<N>, VolumeTransportError> {
    check()?;
    let mut result = SegmentTransport {
        transmittance: [1.0; N],
        emission_w_m2_sr_nm: [0.0; N],
    };
    let Some(mut ray) = CellRay::clip(
        &medium.geometry().grid,
        origin,
        direction,
        max_distance_m,
        medium.geometry().max_steps,
    )?
    else {
        return Ok(result);
    };
    let mut depth = [0.0; N];
    let mut source_cache = [(f64::NAN, 0.0); N];
    while let Some(segment) = ray.next_with_checkpoint(check)? {
        let rho = medium.geometry().grid.data[segment.cell];
        if rho == 0.0 || segment.length_m == 0.0 {
            continue;
        }
        let (absorption, _, extinction) = medium.coefficients(segment.cell)?;
        let tau = extinction.map(|k| {
            if k == 0.0 {
                0.0
            } else {
                positive_product(rho, k, segment.length_m)
            }
        });
        if include_emission {
            integrate_source_cell(
                medium.geometry(),
                segment.cell,
                &tau,
                &absorption,
                &extinction,
                medium.wavelengths(),
                &mut depth,
                &mut source_cache,
                &mut result.emission_w_m2_sr_nm,
            )?;
        } else {
            for lane in 0..N {
                depth[lane] += tau[lane];
            }
        }
    }
    result.transmittance = depth.map(|tau| det::exp(-tau));
    check()?;
    Ok(result)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn integrate_source_cell<const N: usize>(
    geometry: &AbsorbingGrid<'_>,
    cell: usize,
    tau: &[f64; N],
    absorption: &[f64; N],
    extinction: &[f64; N],
    wavelengths_nm: &[f64; N],
    depth: &mut [f64; N],
    source_cache: &mut [(f64, f64); N],
    emission: &mut [f64; N],
) -> Result<(), VolumeTransportError> {
    for lane in 0..N {
        if absorption[lane] > 0.0
            && let Some(temperatures) = geometry.temperature_k
        {
            let opacity = -det::expm1(-tau[lane]);
            if opacity > 0.0 && depth[lane].is_finite() {
                if source_cache[lane].0 != temperatures[cell] {
                    source_cache[lane] = (
                        temperatures[cell],
                        planck_radiance_w_m2_sr_nm(temperatures[cell], wavelengths_nm[lane])
                            .map_err(|_| VolumeTransportError::Unrepresentable)?,
                    );
                }
                let source = source_cache[lane].1;
                let fraction = absorption[lane] / extinction[lane];
                let mut contribution = det::exp(-depth[lane]) * opacity * fraction * source;
                if contribution == 0.0 && source > 0.0 {
                    contribution = det::exp(
                        -depth[lane] + det::ln(opacity) + det::ln(absorption[lane])
                            - det::ln(extinction[lane])
                            + det::ln(source),
                    );
                }
                emission[lane] += contribution;
                if !emission[lane].is_finite() {
                    return Err(VolumeTransportError::Unrepresentable);
                }
            }
        }
        depth[lane] += tau[lane];
    }
    Ok(())
}

/// Allocation-free current-main palette resolution for renderer ray segments.
/// The last material pair is cached on the stack; spectra are only reevaluated
/// when a physical crossing changes that pair.
pub(super) struct LiveOptics<'g, 'a, const N: usize> {
    medium: &'g ScatteringGrid<'a>,
    wavelengths: [f64; N],
    cache: std::cell::Cell<Option<((u32, u32), ([f64; N], [f64; N], [f64; N]))>>,
}
impl<'g, 'a, const N: usize> LiveOptics<'g, 'a, N> {
    pub fn new(
        medium: &'g ScatteringGrid<'a>,
        wavelengths: [f64; N],
        coefficients: (&[f64; N], &[f64; N], &[f64; N]),
    ) -> Self {
        let key = (
            medium
                .absorption
                .absorption_palette()
                .map_or(0, |p| p.cell_indices()[0]),
            medium
                .scattering_palette()
                .map_or(0, |p| p.cell_indices()[0]),
        );
        Self {
            medium,
            wavelengths,
            cache: std::cell::Cell::new(Some((
                key,
                (*coefficients.0, *coefficients.1, *coefficients.2),
            ))),
        }
    }
}
impl<const N: usize> PreparedOptics<N> for LiveOptics<'_, '_, N> {
    fn geometry(&self) -> &AbsorbingGrid<'_> {
        &self.medium.absorption
    }
    fn identity(&self) -> ContentHash {
        self.medium.identity
    }
    fn bounds(&self) -> (f64, f64) {
        self.medium.wavelength_bounds_nm()
    }
    fn wavelengths(&self) -> &[f64; N] {
        &self.wavelengths
    }
    fn coefficient(
        &self,
        cell: usize,
        lane: usize,
    ) -> Result<(f64, f64, f64), VolumeTransportError> {
        let (a, s, t) = self.coefficients(cell)?;
        Ok((a[lane], s[lane], t[lane]))
    }
    fn coefficients(
        &self,
        cell: usize,
    ) -> Result<([f64; N], [f64; N], [f64; N]), VolumeTransportError> {
        let key = self.medium.coefficient_key(cell)?;
        if let Some((last, values)) = self.cache.get()
            && key == last
        {
            return Ok(values);
        }
        let values = self.medium.coefficients_at_cell(cell, &self.wavelengths)?;
        self.cache.set(Some((key, values)));
        Ok(values)
    }
    fn phase(&self, cell: usize) -> &PhaseFunction {
        self.medium
            .phase_at_cell(cell)
            .expect("admitted original DDA cell")
    }
    fn has_scattering(&self) -> bool {
        self.medium.has_scattering
    }
    // This renderer-only view conservatively retains all physical paths. The
    // deposition view below prepares exact occupied-cell lane flags once.
    fn lane_has_scattering(&self, _: usize) -> bool {
        self.medium.has_scattering
    }
    fn lane_has_extinction(&self, _: usize) -> bool {
        !self.medium.is_vacuum()
    }
    fn has_emission(&self) -> bool {
        self.medium.has_thermal_emission()
    }
}

/// Photon histories amortize each distinct native spectrum's preparation over
/// the complete operation. The original cell maps and phase palette stay
/// borrowed; storage is O(distinct models × packet width), not O(cells × knots).
pub(super) struct PaletteOptics<'g, 'a, const N: usize> {
    medium: &'g ScatteringGrid<'a>,
    wavelengths: [f64; N],
    absorption: Vec<[f64; N]>,
    scattering: Vec<[f64; N]>,
    scattering_lanes: [bool; N],
    extinction_lanes: [bool; N],
}
impl<'g, 'a, const N: usize> PaletteOptics<'g, 'a, N> {
    pub fn try_new(
        medium: &'g ScatteringGrid<'a>,
        wavelengths: [f64; N],
        check: &mut (impl FnMut() -> Result<(), VolumeTransportError> + ?Sized),
    ) -> Result<Self, VolumeTransportError> {
        check()?;
        if N == 0 {
            return Err(VolumeTransportError::Sample);
        }
        let a_models = medium
            .absorption
            .absorption_palette()
            .map_or(1, |p| p.values().len());
        let s_models = medium.scattering_palette().map_or(1, |p| p.values().len());
        let mut absorption = Vec::new();
        let mut scattering = Vec::new();
        absorption
            .try_reserve_exact(a_models)
            .map_err(|_| VolumeTransportError::Budget)?;
        scattering
            .try_reserve_exact(s_models)
            .map_err(|_| VolumeTransportError::Budget)?;
        for index in 0..a_models {
            check()?;
            let spectrum = medium
                .absorption
                .absorption_palette()
                .map_or(&medium.absorption.absorption, |p| &p.values()[index]);
            let mut values = [0.0; N];
            for lane in 0..N {
                values[lane] = spectrum.at(wavelengths[lane])?;
            }
            absorption.push(values);
        }
        for index in 0..s_models {
            check()?;
            let spectrum = medium
                .scattering_palette()
                .map_or(&medium.scattering, |p| &p.values()[index].scattering);
            let mut values = [0.0; N];
            for lane in 0..N {
                values[lane] = spectrum.at(wavelengths[lane])?;
            }
            scattering.push(values);
        }
        let mut result = Self {
            medium,
            wavelengths,
            absorption,
            scattering,
            scattering_lanes: [false; N],
            extinction_lanes: [false; N],
        };
        for (cell, &rho) in medium.absorption.grid.data.iter().enumerate() {
            if cell % 1024 == 0 {
                check()?;
            }
            if rho == 0.0 {
                continue;
            }
            let (a, s) = result.packets(cell);
            for lane in 0..N {
                let t = a[lane] + s[lane];
                if !t.is_finite() {
                    return Err(VolumeTransportError::Unrepresentable);
                }
                result.scattering_lanes[lane] |= s[lane] > 0.0;
                result.extinction_lanes[lane] |= t > 0.0;
            }
        }
        check()?;
        Ok(result)
    }
    fn packets(&self, cell: usize) -> ([f64; N], [f64; N]) {
        let a = self
            .medium
            .absorption
            .absorption_palette()
            .map_or(0, |p| p.cell_indices()[cell] as usize);
        let s = self
            .medium
            .scattering_palette()
            .map_or(0, |p| p.cell_indices()[cell] as usize);
        (self.absorption[a], self.scattering[s])
    }
}
impl<const N: usize> PreparedOptics<N> for PaletteOptics<'_, '_, N> {
    fn geometry(&self) -> &AbsorbingGrid<'_> {
        &self.medium.absorption
    }
    fn identity(&self) -> ContentHash {
        self.medium.identity
    }
    fn bounds(&self) -> (f64, f64) {
        self.medium.wavelength_bounds_nm()
    }
    fn wavelengths(&self) -> &[f64; N] {
        &self.wavelengths
    }
    fn coefficient(
        &self,
        cell: usize,
        lane: usize,
    ) -> Result<(f64, f64, f64), VolumeTransportError> {
        let a = self
            .medium
            .absorption
            .absorption_palette()
            .map_or(0, |p| p.cell_indices()[cell] as usize);
        let s = self
            .medium
            .scattering_palette()
            .map_or(0, |p| p.cell_indices()[cell] as usize);
        let a = self.absorption[a][lane];
        let s = self.scattering[s][lane];
        Ok((a, s, a + s))
    }
    fn coefficients(
        &self,
        cell: usize,
    ) -> Result<([f64; N], [f64; N], [f64; N]), VolumeTransportError> {
        let (a, s) = self.packets(cell);
        Ok((a, s, std::array::from_fn(|lane| a[lane] + s[lane])))
    }
    fn phase(&self, cell: usize) -> &PhaseFunction {
        self.medium
            .phase_at_cell(cell)
            .expect("admitted original DDA cell")
    }
    fn has_scattering(&self) -> bool {
        self.medium.has_scattering
    }
    fn lane_has_scattering(&self, lane: usize) -> bool {
        self.scattering_lanes[lane]
    }
    fn lane_has_extinction(&self, lane: usize) -> bool {
        self.extinction_lanes[lane]
    }
    fn has_emission(&self) -> bool {
        self.medium.has_thermal_emission()
    }
}
