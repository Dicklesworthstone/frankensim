//! Spectral multiple-scattering building blocks on borrowed physical cells.
//!
//! A cell has extinction `rho * (kappa_a + kappa_s)` in m⁻¹. Exact DDA
//! crossings integrate its piecewise-constant optical depth and invert one
//! exponential free flight; there is no step-size bias or global majorant.
//! One uniformly chosen spectral lane proposes the distance, but every lane
//! uses the complete mixture density. This is the spectral-channel estimator
//! in PBRT, *Sampling Volume Scattering*, generalized to exact cell integrals:
//! <https://pbr-book.org/3ed-2018/Light_Transport_II_Volume_Rendering/Sampling_Volume_Scattering>.
//!
//! LTE emission is integrated over the complete clipped segment independently
//! of the sampled collision. Its source is `rho * kappa_a * B_lambda(T)`;
//! scattering never creates thermal emission. Add that source with the incoming
//! path throughput **before** applying the returned stochastic continuation
//! weight. This splits the source and scattering terms of the transport
//! equation without double counting or terminating emission at a random event.

use super::cells::{CellRay, validate_ray};
use super::{
    AbsorbingGrid, CellPalette, MassAbsorptionSpectrum, SegmentTransport, VolumeTransportError,
    checkpoint, positive_product,
};
use crate::thermal_surface::planck_radiance_w_m2_sr_nm;
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_math::det;

#[cfg(feature = "photometry")]
mod contraction;
mod indexed;
mod phase;
pub use indexed::ScatteringOptics;
pub use phase::{
    MAX_PHASE_KNOTS, MAX_PHASE_WAVELENGTHS, PhaseFunction, PhasePacketEvaluation,
    PhasePacketSample, PhaseSample, SpectralPhaseTable, TabulatedPhase,
};
/// Incident power absorbed after any number of native scattering events.
pub mod deposition;
/// Physical homogeneous-sphere cross sections and resolved angular scattering.
pub mod mie;

/// Mass scattering coefficient in m²/kg, independently specified
/// from absorption. The interpolation and admitted wavelength band are the
/// same as [`MassAbsorptionSpectrum`].
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MassScatteringSpectrum {
    // Tabulated native coefficients, or the exact power-law value at its
    // reference wavelength (retained in both slots for canonical hashing).
    spectrum: MassAbsorptionSpectrum,
    rayleigh_reference_nm: Option<f64>,
}

impl MassScatteringSpectrum {
    /// Wavelength-independent scattering over 380–780 nm.
    pub fn constant(value_m2_kg: f64) -> Result<Self, VolumeTransportError> {
        MassAbsorptionSpectrum::constant(value_m2_kg)
            .map(|spectrum| Self {
                spectrum,
                rayleigh_reference_nm: None,
            })
            .map_err(|_| VolumeTransportError::ScatteringSpectrum)
    }

    /// Linear interpolation of 2–16 strictly increasing `(vacuum_nm, m²/kg)`
    /// knots with endpoints exactly 380 and 780 nm.
    pub fn table(knots: &[(f64, f64)]) -> Result<Self, VolumeTransportError> {
        MassAbsorptionSpectrum::table(knots)
            .map(|spectrum| Self {
                spectrum,
                rayleigh_reference_nm: None,
            })
            .map_err(|_| VolumeTransportError::ScatteringSpectrum)
    }

    /// Retain a source table over its exact positive vacuum-wavelength support.
    /// Coefficients are in m²/kg; no UV/IR continuation is inferred. Source and
    /// interpolation semantics are those of [`MassAbsorptionSpectrum::source_table`].
    pub fn source_table(
        knots: &[(f64, f64)],
        source_identity: ContentHash,
    ) -> Result<Self, VolumeTransportError> {
        MassAbsorptionSpectrum::source_table(knots, source_identity)
            .map(|spectrum| Self {
                spectrum,
                rayleigh_reference_nm: None,
            })
            .map_err(|error| match error {
                VolumeTransportError::AbsorptionSpectrum => {
                    VolumeTransportError::ScatteringSpectrum
                }
                other => other,
            })
    }

    /// Exact inclusive source/model support, in vacuum nm.
    #[must_use]
    pub const fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.spectrum.wavelength_bounds_nm()
    }

    /// Original explicitly supplied source identity, when present.
    #[must_use]
    pub const fn source_identity(&self) -> Option<ContentHash> {
        self.spectrum.source_identity()
    }

    /// Evaluate the admitted coefficient without wavelength extrapolation.
    pub fn at(&self, wavelength_nm: f64) -> Result<f64, VolumeTransportError> {
        let value = self.spectrum.at(wavelength_nm)?;
        let Some(reference) = self.rayleigh_reference_nm else {
            return Ok(value);
        };
        let ratio = reference / wavelength_nm;
        let square = ratio * ratio;
        let result = value * square * square;
        if !result.is_finite() || (value > 0.0 && result == 0.0) {
            return Err(VolumeTransportError::Unrepresentable);
        }
        Ok(result)
    }

    /// Exact `kappa(reference) * (reference / lambda)^4` molecular Rayleigh
    /// spectrum, without piecewise-linear wavelength bias. The reference must
    /// be in the visible band. Refractivity and depolarization are assumed
    /// wavelength independent; this is not a dispersion or resonance model.
    pub fn rayleigh(reference_nm: f64, value_m2_kg: f64) -> Result<Self, VolumeTransportError> {
        if !reference_nm.is_finite()
            || !(crate::spectral::LAMBDA_MIN..=crate::spectral::LAMBDA_MAX).contains(&reference_nm)
        {
            return Err(VolumeTransportError::Wavelength);
        }
        let result = Self {
            spectrum: MassAbsorptionSpectrum::constant(value_m2_kg)
                .map_err(|_| VolumeTransportError::ScatteringSpectrum)?,
            rayleigh_reference_nm: Some(reference_nm),
        };
        result.at(crate::spectral::LAMBDA_MIN)?;
        result.at(crate::spectral::LAMBDA_MAX)?;
        Ok(result)
    }

    fn knots(&self) -> &[(f64, f64)] {
        &self.spectrum.knots[..self.spectrum.count]
    }
}

/// A sampled interaction in the physical medium, before any surface endpoint.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ScatteringEvent {
    /// X-major index of the original density cell containing this collision.
    pub cell_index: usize,
    /// World-space point in metres, computed from the rebased grid entry.
    pub position: Point3,
    /// Physical distance from the requested segment origin, in metres.
    pub distance_m: f64,
}

/// One unbiased spectral free-flight estimate and an independent source term.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeFlightSample<const N: usize> {
    /// `Some` continues through a phase-function vertex; `None` reaches the
    /// original surface/environment endpoint without a scattering vertex.
    pub event: Option<ScatteringEvent>,
    /// Collision: `T_lambda * kappa_s_lambda / mean(T * kappa_t)`.
    /// Escape: `T_lambda / mean(T)`. Snapshots with no scattering use
    /// deterministic attenuation instead. Never attenuate this weight again.
    pub weight: [f64; N],
    /// Complete foreground LTE source integral, W/(m² sr nm), added using
    /// incoming throughput before the stochastic continuation weight.
    pub emission_w_m2_sr_nm: [f64; N],
    /// Probability of sampled escape, `mean(T_lambda)`, on an escape result;
    /// zero on a collision. Deterministic absorption-only transport has value
    /// one. This supplies the survival factor for complete directional MIS.
    pub escape_probability: f64,
}

/// An immutable, index-one participating medium. Density and temperature are
/// borrowed from the accepted absorption snapshot. Scattering and angular
/// response may be global or selected by a compact immutable cell palette.
pub struct ScatteringGrid<'a> {
    absorption: AbsorbingGrid<'a>,
    scattering: MassScatteringSpectrum,
    phase: PhaseFunction,
    scattering_palette: Option<CellPalette<ScatteringOptics>>,
    identity: ContentHash,
    has_scattering: bool,
}

impl<'a> ScatteringGrid<'a> {
    /// Extend an already validated physical snapshot. Spatial density multiplies
    /// both spectra; temperature controls only the absorption-linked LTE source.
    /// The absorption snapshot's cell and traversal budgets remain binding.
    /// Coefficient domains must match; a spectral phase table must cover that
    /// entire domain. Missing wavelengths are never treated as transparent.
    pub fn try_new(
        cx: &Cx<'_>,
        absorption: AbsorbingGrid<'a>,
        scattering: MassScatteringSpectrum,
        phase: PhaseFunction,
    ) -> Result<Self, VolumeTransportError> {
        Self::prepare(cx, absorption, scattering, phase, None)
    }

    /// Borrow the original absorption owner for the exact zero-scattering path.
    #[must_use]
    pub const fn absorption_grid(&self) -> &AbsorbingGrid<'a> {
        &self.absorption
    }

    /// Numerical snapshot, scattering spectrum and phase-law identity.
    #[must_use]
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }

    /// Shared inclusive support of absorption, scattering and the phase law.
    #[must_use]
    pub const fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.absorption.wavelength_bounds_nm()
    }

    /// Whether at least one wavelength and cell can scatter.
    #[must_use]
    pub const fn has_scattering(&self) -> bool {
        self.has_scattering
    }

    /// Every wavelength and cell has zero extinction and zero emission.
    #[must_use]
    pub fn is_vacuum(&self) -> bool {
        self.absorption.is_vacuum() && !self.has_scattering
    }

    /// Positive absorption and an accepted temperature field supply LTE emission.
    #[must_use]
    pub fn has_thermal_emission(&self) -> bool {
        self.absorption.has_thermal_emission()
    }

    /// The global angular law, if one was declared. Indexed snapshots return
    /// `None`; their physical response must be queried with `phase_at_cell`.
    #[must_use]
    pub fn phase(&self) -> Option<PhaseFunction> {
        self.phase_ref().cloned()
    }

    /// Borrow a declared global phase without cloning its angular tables.
    /// An indexed medium has no implicit representative angular law.
    #[must_use]
    pub const fn phase_ref(&self) -> Option<&PhaseFunction> {
        if self.scattering_palette.is_some() {
            None
        } else {
            Some(&self.phase)
        }
    }

    /// Resolve a collimated beam's attenuation into spatial absorption and
    /// scattering, using the exact physical cell crossings.
    ///
    /// Inputs are already spectrally weighted watts, with a unit ray direction
    /// and metre distances. Only `absorbed_power_w` is heat. The scattered
    /// power leaves this directed beam and must be transported separately;
    /// this operation does not approximate multiple-scattering heating or
    /// subtract the medium's independent thermal emission.
    ///
    /// `max_deposition_cells` bounds retained cell records in addition to the
    /// accepted medium's traversal budget. No partial report is returned.
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
    ) -> Result<super::deposition::BeamPowerDeposition<N>, VolumeTransportError> {
        super::deposition::deposit(
            cx,
            &self.absorption,
            (origin, unit_direction, max_distance_m),
            wavelengths_nm,
            incident_power_w,
            max_deposition_cells,
            Some(self),
        )
    }

    fn coefficients<const N: usize>(
        &self,
        wavelengths: &[f64; N],
    ) -> Result<([f64; N], [f64; N], [f64; N]), VolumeTransportError> {
        self.coefficients_at_cell(0, wavelengths)
    }

    /// Exact extinction on a physical shadow ray, without emission or random
    /// sampling. Positive infinity is valid for environment-bound directions.
    pub fn transmittance<const N: usize>(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        direction: Vec3,
        max_distance_m: f64,
        wavelengths_nm: &[f64; N],
    ) -> Result<[f64; N], VolumeTransportError> {
        checkpoint(cx)?;
        validate_ray(origin, direction, max_distance_m)?;
        let (_, _, mut extinction) = self.coefficients(wavelengths_nm)?;
        let mut coefficient_key = self.coefficient_key(0)?;
        if !self.has_scattering {
            return self.absorption.transmittance(
                cx,
                origin,
                direction,
                max_distance_m,
                wavelengths_nm,
            );
        }
        let Some(mut ray) = CellRay::clip(
            &self.absorption.grid,
            origin,
            direction,
            max_distance_m,
            self.absorption.max_steps,
        )?
        else {
            return Ok([1.0; N]);
        };
        let mut depth = [0.0; N];
        while let Some(segment) = ray.next(cx)? {
            let rho = self.absorption.grid.data[segment.cell];
            if rho > 0.0 && segment.length_m > 0.0 {
                let key = self.coefficient_key(segment.cell)?;
                if key != coefficient_key {
                    extinction = self.coefficients_at_cell(segment.cell, wavelengths_nm)?.2;
                    coefficient_key = key;
                }
                for lane in 0..N {
                    if extinction[lane] > 0.0 {
                        depth[lane] += positive_product(rho, extinction[lane], segment.length_m);
                    }
                }
            }
        }
        Ok(depth.map(|tau| det::exp(-tau)))
    }

    /// Exact attenuation and LTE source on a segment, without sampling a collision.
    ///
    /// Direct source-light sampling needs the deterministic emission term of the
    /// transfer equation. Scattering contributes extinction here, while only
    /// absorption creates thermal radiation. The stochastic in-scattering term
    /// remains the responsibility of [`Self::sample_segment`].
    pub fn transport<const N: usize>(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        direction: Vec3,
        max_distance_m: f64,
        wavelengths_nm: &[f64; N],
    ) -> Result<SegmentTransport<N>, VolumeTransportError> {
        checkpoint(cx)?;
        validate_ray(origin, direction, max_distance_m)?;
        let (mut absorption, _, mut extinction) = self.coefficients(wavelengths_nm)?;
        let mut coefficient_key = self.coefficient_key(0)?;
        if !self.has_scattering {
            return self.absorption.transport(
                cx,
                origin,
                direction,
                max_distance_m,
                wavelengths_nm,
            );
        }
        let mut result = SegmentTransport {
            transmittance: [1.0; N],
            emission_w_m2_sr_nm: [0.0; N],
        };
        let Some(mut ray) = CellRay::clip(
            &self.absorption.grid,
            origin,
            direction,
            max_distance_m,
            self.absorption.max_steps,
        )?
        else {
            return Ok(result);
        };
        let mut depth = [0.0; N];
        let mut source_cache = [(f64::NAN, 0.0); N];
        while let Some(segment) = ray.next(cx)? {
            let rho = self.absorption.grid.data[segment.cell];
            if rho == 0.0 || segment.length_m == 0.0 {
                continue;
            }
            let key = self.coefficient_key(segment.cell)?;
            if key != coefficient_key {
                (absorption, _, extinction) =
                    self.coefficients_at_cell(segment.cell, wavelengths_nm)?;
                coefficient_key = key;
            }
            let tau = extinction.map(|k| {
                if k == 0.0 {
                    0.0
                } else {
                    positive_product(rho, k, segment.length_m)
                }
            });
            self.integrate_source_cell(
                segment.cell,
                &tau,
                &absorption,
                &extinction,
                wavelengths_nm,
                &mut depth,
                &mut source_cache,
                &mut result.emission_w_m2_sr_nm,
            )?;
        }
        result.transmittance = depth.map(|tau| det::exp(-tau));
        Ok(result)
    }

    /// Sample a complete spectral mixture free flight, without allocating.
    ///
    /// Both variates lie in `[0,1)`: one selects a uniform packet lane and the
    /// other is a survival variate for its exponential optical depth. A zero
    /// variate uses `f64::MIN_POSITIVE`, avoiding an infinite sampled distance.
    /// Zero-extinction lanes remain valid proposals that always escape. The
    /// bounded DDA stops at an event if there is no thermal source; otherwise it
    /// finishes the complete source integral before returning the event.
    #[allow(clippy::too_many_arguments)]
    pub fn sample_segment<const N: usize>(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        direction: Vec3,
        max_distance_m: f64,
        wavelengths_nm: &[f64; N],
        u_lane: f64,
        u_distance: f64,
    ) -> Result<FreeFlightSample<N>, VolumeTransportError> {
        checkpoint(cx)?;
        validate_ray(origin, direction, max_distance_m)?;
        if !unit_variate(u_lane) || !unit_variate(u_distance) {
            return Err(VolumeTransportError::Sample);
        }
        let (absorption, scattering, extinction) = self.coefficients(wavelengths_nm)?;
        self.sample_segment_prepared(
            cx,
            origin,
            direction,
            max_distance_m,
            wavelengths_nm,
            u_lane,
            u_distance,
            (&absorption, &scattering, &extinction),
            true,
        )
    }

    // Reuse the exact free-flight kernel for source-power histories. Their
    // coefficients are prepared once, and incident heating must not include
    // the medium's independent LTE source or spend work integrating it.
    #[allow(clippy::too_many_arguments, clippy::type_complexity)]
    fn sample_segment_prepared<const N: usize>(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        direction: Vec3,
        max_distance_m: f64,
        wavelengths_nm: &[f64; N],
        u_lane: f64,
        u_distance: f64,
        coefficients: (&[f64; N], &[f64; N], &[f64; N]),
        include_emission: bool,
    ) -> Result<FreeFlightSample<N>, VolumeTransportError> {
        checkpoint(cx)?;
        let (mut absorption, mut scattering, mut extinction) =
            (*coefficients.0, *coefficients.1, *coefficients.2);
        let mut coefficient_key = self.coefficient_key(0)?;
        if !self.has_scattering {
            let result = if include_emission {
                self.absorption
                    .transport(cx, origin, direction, max_distance_m, wavelengths_nm)?
            } else {
                SegmentTransport {
                    transmittance: self.absorption.transmittance(
                        cx,
                        origin,
                        direction,
                        max_distance_m,
                        wavelengths_nm,
                    )?,
                    emission_w_m2_sr_nm: [0.0; N],
                }
            };
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
            &self.absorption.grid,
            origin,
            direction,
            max_distance_m,
            self.absorption.max_steps,
        )?
        else {
            return Ok(result);
        };
        let hero = ((u_lane * N as f64) as usize).min(N - 1);
        let mut remaining_tau = -det::ln(u_distance.max(f64::MIN_POSITIVE));
        let mut depth = [0.0; N];
        let mut source_cache = [(f64::NAN, 0.0); N];
        while let Some(segment) = ray.next(cx)? {
            let rho = self.absorption.grid.data[segment.cell];
            if rho == 0.0 || segment.length_m == 0.0 {
                continue;
            }
            let key = self.coefficient_key(segment.cell)?;
            if key != coefficient_key {
                (absorption, scattering, extinction) =
                    self.coefficients_at_cell(segment.cell, wavelengths_nm)?;
                coefficient_key = key;
            }
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
                if !include_emission || !self.has_thermal_emission() {
                    return Ok(result);
                }
            }
            if result.event.is_none() {
                remaining_tau -= tau[hero];
            }
            if include_emission {
                self.integrate_source_cell(
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
            result.escape_probability =
                depth.iter().map(|tau| det::exp(-tau)).sum::<f64>() / N as f64;
        }
        Ok(result)
    }

    #[allow(clippy::too_many_arguments)]
    fn integrate_source_cell<const N: usize>(
        &self,
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
                && let Some(temperatures) = self.absorption.temperature_k
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
}

fn collision_weights<const N: usize>(
    depth: &[f64; N],
    extinction: &[f64; N],
    scattering: &[f64; N],
) -> Result<[f64; N], VolumeTransportError> {
    let log_density =
        std::array::from_fn::<_, N, _>(|lane| det::ln(extinction[lane]) - depth[lane]);
    let peak = log_density.into_iter().fold(f64::NEG_INFINITY, f64::max);
    if !peak.is_finite() {
        return Err(VolumeTransportError::Unrepresentable);
    }
    let sum: f64 = log_density
        .iter()
        .map(|value| det::exp(*value - peak))
        .sum();
    let weights = std::array::from_fn(|lane| {
        if scattering[lane] == 0.0 {
            0.0
        } else {
            det::exp(det::ln(scattering[lane]) - depth[lane] - peak) * (N as f64 / sum)
        }
    });
    if weights.iter().all(|weight| weight.is_finite()) {
        Ok(weights)
    } else {
        Err(VolumeTransportError::Unrepresentable)
    }
}

fn divide_product(numerator: f64, a: f64, b: f64) -> f64 {
    if numerator == 0.0 {
        return 0.0;
    }
    let denominator = a * b;
    let result = numerator / denominator;
    if denominator >= f64::MIN_POSITIVE && result > 0.0 && result.is_finite() {
        result
    } else {
        det::exp(det::ln(numerator) - det::ln(a) - det::ln(b))
    }
}

fn unit_variate(u: f64) -> bool {
    u.is_finite() && (0.0..1.0).contains(&u)
}

#[cfg(test)]
mod tests;
