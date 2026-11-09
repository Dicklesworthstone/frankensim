//! Observe original luminescent cell powers without assigning a temperature.
//!
//! Each supplied power is integrated over one discrete spectral node/bin and
//! all directions, in watts. It is spread uniformly through its original cell
//! and emitted isotropically: `j = P / (4 pi V)`, in W/(m³ sr). The native DDA
//! integrates `dL/ds = -(rho*kappa + alpha_active)*L + j` exactly for that cellwise model,
//! including emission when absorption is zero. No spectral width is inferred.
//!
//! This is a stationary, non-scattering, index-one line-of-sight observation.
//! Absorption attenuates the observed light; this view does not feed that energy
//! into thermal cells or create secondary fluorescence. Native thermal emission
//! is a separate additive component. Supplied powers must be accepted endpoint
//! values, or interval averages for fixed geometry and optics. There is no
//! time-of-flight, directional dipole emission, subcell pump profile, surface
//! reflection/refraction or calibrated source data inferred by this adapter.
//! [`ScatteringLuminescenceGrid`] explicitly extends the same accepted source
//! with elastic multiple scattering; the direct observer retains its semantics.

use super::cells::{CellRay, validate_ray};
use super::mixture::PositiveSum;
use super::scattering::extinction::log_total_linear_extinction;
use super::{AbsorbingGrid, VolumeTransportError, checkpoint, positive_product};
use fs_blake3::{ContentHash, DomainHasher};
use fs_exec::Cx;
use fs_geom::{Point3, Vec3};
use fs_math::det;

mod camera;
mod exposure;
mod scattered;
pub use camera::{
    LuminescenceCameraError, LuminescenceCameraImage, LuminescenceCameraPixel,
    LuminescenceCameraRadianceImage, LuminescenceCameraSettings, VolumeCameraBackground,
    render_luminescence_camera, render_luminescence_camera_with_background,
    render_scattered_luminescence_camera,
};
pub use exposure::{
    LuminescenceCameraExposure, LuminescenceExposureBudget, LuminescenceExposureError,
    LuminescenceExposureField, LuminescenceExposureNode, render_luminescence_exposure,
};
pub use scattered::{
    LuminescenceScatteringBudget, ScatteredLuminescenceSample, ScatteringLuminescenceGrid,
};

const SOLID_ANGLE_SR: f64 = 4.0 * core::f64::consts::PI;

/// Nonthermal foreground radiance and attenuation for the declared nodes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LuminescenceTransport<const M: usize> {
    /// Multiplicative attenuation of line/bin radiance at the far endpoint.
    pub transmittance: [f64; M],
    /// Observed emitted line/bin radiance, W/(m² sr), NOT per nanometre.
    pub emitted_radiance_w_m2_sr: [f64; M],
}

// One native traversal supplies the public ray result and camera work/depth.
// Finite depth is retained independently of underflowed unit transmission.
struct LuminescenceRaySample<const M: usize> {
    transport: LuminescenceTransport<M>,
    optical_depth: [f64; M],
    cell_visits: u32,
}

/// Sparse immutable emissivity bound to the original optical grid.
///
/// Preparation is O(K log K) for K supplied cells; each ray allocates nothing,
/// uses the native cell-step/cancellation budget and looks up only crossed
/// emitting cells. The original density field and optical metadata are borrowed.
pub struct LuminescenceGrid<'m, 'g, const M: usize> {
    medium: &'m AbsorbingGrid<'g>,
    source_identity: ContentHash,
    wavelengths_nm: [f64; M],
    emissivity: Vec<(usize, [f64; M])>,
    additional_absorption: Vec<(usize, [f64; M])>,
    identity: ContentHash,
}

impl<'m, 'g, const M: usize> LuminescenceGrid<'m, 'g, M> {
    /// Admit 1–16 increasing source-supported vacuum wavelengths and sparse
    /// original-cell powers. Input cell order is arbitrary; duplicate or unknown
    /// indices refuse. Zero-power entries are allowed and count toward the
    /// explicit preparation bound. `source_identity` binds the photophysical
    /// source/condition; it does not confer measurement authority.
    ///
    /// The declared spatial model is uniform isotropic emission in each grid
    /// cell, even when its incident pump was a collimated ray. The caller must
    /// resolve the pumped support accordingly. Cell volume and positive
    /// emissivities must be representable. Invalid cell mappings use `Palette`;
    /// invalid emitted powers use the existing spectral `BeamPower` diagnostic.
    pub fn try_new(
        cx: &Cx<'_>,
        medium: &'m AbsorbingGrid<'g>,
        source_identity: ContentHash,
        wavelengths_nm: [f64; M],
        cell_powers_w: &[(usize, [f64; M])],
        max_emitting_cells: usize,
    ) -> Result<Self, VolumeTransportError> {
        checkpoint(cx)?;
        if !(1..=16).contains(&M) {
            return Err(VolumeTransportError::Sample);
        }
        if source_identity == ContentHash([0; 32]) {
            return Err(VolumeTransportError::SourceIdentity);
        }
        if cell_powers_w.len() > max_emitting_cells {
            return Err(VolumeTransportError::Budget);
        }
        if wavelengths_nm.windows(2).any(|w| w[0] >= w[1]) {
            return Err(VolumeTransportError::Wavelength);
        }
        for wavelength in wavelengths_nm {
            medium.absorption.at(wavelength)?;
        }
        let volume = positive_product(
            medium.grid.cell[0],
            medium.grid.cell[1],
            medium.grid.cell[2],
        );
        if !volume.is_finite() || volume <= 0.0 {
            return Err(VolumeTransportError::Unrepresentable);
        }
        let mut emissivity = Vec::new();
        emissivity
            .try_reserve_exact(cell_powers_w.len())
            .map_err(|_| VolumeTransportError::Budget)?;
        emissivity.extend_from_slice(cell_powers_w);
        emissivity.sort_unstable_by_key(|(cell, _)| *cell);
        let mut hash = DomainHasher::new("org.frankensim.render.luminescence-grid.v1");
        hash.update(medium.content_identity().as_bytes());
        hash.update(source_identity.as_bytes());
        hash.update(&(M as u64).to_le_bytes());
        for wavelength in wavelengths_nm {
            hash.update(&wavelength.to_bits().to_le_bytes());
        }
        hash.update(&(emissivity.len() as u64).to_le_bytes());
        let mut previous = None;
        for (cell, powers) in &mut emissivity {
            checkpoint(cx)?;
            if *cell >= medium.cell_count() || previous == Some(*cell) {
                return Err(VolumeTransportError::Palette);
            }
            previous = Some(*cell);
            hash.update(&(*cell as u64).to_le_bytes());
            for (lane, power) in powers.iter_mut().enumerate() {
                if !power.is_finite() || *power < 0.0 {
                    return Err(VolumeTransportError::BeamPower { lane });
                }
                hash.update(&power.to_bits().to_le_bytes());
                if *power > 0.0 {
                    let per_sr = *power / SOLID_ANGLE_SR;
                    let coefficient = per_sr / volume;
                    *power = if per_sr >= f64::MIN_POSITIVE
                        && coefficient.is_finite()
                        && coefficient >= f64::MIN_POSITIVE
                    {
                        coefficient
                    } else {
                        det::exp(det::ln(*power) - det::ln(SOLID_ANGLE_SR) - det::ln(volume))
                    };
                    if !power.is_finite() || *power <= 0.0 {
                        return Err(VolumeTransportError::Unrepresentable);
                    }
                }
            }
        }
        checkpoint(cx)?;
        Ok(Self {
            medium,
            source_identity,
            wavelengths_nm,
            emissivity,
            additional_absorption: Vec::new(),
            identity: hash.finalize(),
        })
    }

    /// Add accepted, cellwise absorption in inverse metres at the original
    /// spectral nodes. This supports finite molecular populations whose active
    /// opacity is separate from the host's passive mass absorption. Absorbing
    /// cells need not emit. The caller must exclude this same contribution from
    /// the host and supply the opacity from the same accepted state as powers.
    ///
    /// The sparse immutable binding is admitted once. Duplicate/unknown cells,
    /// invalid coefficients, cancellation or a preparation cap refuse before
    /// publishing the snapshot. Changing opacity during an exposure requires
    /// time-resolved observations; mean opacity and mean power are insufficient.
    /// The multiple-scattering adapter uses this same opacity in its flights,
    /// spectral continuation weights and complete source integrals.
    pub fn with_additional_absorption(
        mut self,
        cx: &Cx<'_>,
        absorption_identity: ContentHash,
        cell_absorption_m_inv: &[(usize, [f64; M])],
        max_absorbing_cells: usize,
    ) -> Result<Self, VolumeTransportError> {
        checkpoint(cx)?;
        if absorption_identity == ContentHash([0; 32]) || !self.additional_absorption.is_empty() {
            return Err(VolumeTransportError::SourceIdentity);
        }
        if cell_absorption_m_inv.len() > max_absorbing_cells {
            return Err(VolumeTransportError::Budget);
        }
        let mut absorption = Vec::new();
        absorption
            .try_reserve_exact(cell_absorption_m_inv.len())
            .map_err(|_| VolumeTransportError::Budget)?;
        absorption.extend_from_slice(cell_absorption_m_inv);
        absorption.sort_unstable_by_key(|(cell, _)| *cell);
        let mut hash = DomainHasher::new("org.frankensim.render.luminescence-active-opacity.v1");
        hash.update(self.identity.as_bytes());
        hash.update(absorption_identity.as_bytes());
        hash.update(&(absorption.len() as u64).to_le_bytes());
        let mut previous = None;
        for (cell, coefficients) in &absorption {
            checkpoint(cx)?;
            if *cell >= self.medium.cell_count() || previous == Some(*cell) {
                return Err(VolumeTransportError::Palette);
            }
            previous = Some(*cell);
            hash.update(&(*cell as u64).to_le_bytes());
            for coefficient in coefficients {
                if !coefficient.is_finite() || *coefficient < 0.0 {
                    return Err(VolumeTransportError::AbsorptionSpectrum);
                }
                hash.update(&coefficient.to_bits().to_le_bytes());
            }
        }
        checkpoint(cx)?;
        self.additional_absorption = absorption;
        self.identity = hash.finalize();
        Ok(self)
    }

    /// Snapshot identity including the medium, source, spectrum and cell powers.
    #[must_use]
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }

    /// Supplied photophysical source and condition identity.
    #[must_use]
    pub const fn source_identity(&self) -> ContentHash {
        self.source_identity
    }

    /// Original discrete emission-node/bin wavelengths, in vacuum nm.
    #[must_use]
    pub const fn wavelengths_nm(&self) -> &[f64; M] {
        &self.wavelengths_nm
    }

    /// Observe nonthermal emission along a camera/receiver ray. The direction
    /// points from the observer into the scene; foreground absorption therefore
    /// attenuates more distant emitters. Positive infinite distance is admitted.
    /// No partial result is published on cancellation or work-budget refusal.
    pub fn transport(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        unit_direction: Vec3,
        max_distance_m: f64,
    ) -> Result<LuminescenceTransport<M>, VolumeTransportError> {
        Ok(self
            .transport_counted(
                cx,
                origin,
                unit_direction,
                max_distance_m,
                u64::from(self.medium.max_steps),
            )?
            .transport)
    }

    fn transport_counted(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        unit_direction: Vec3,
        max_distance_m: f64,
        max_cell_visits: u64,
    ) -> Result<LuminescenceRaySample<M>, VolumeTransportError> {
        checkpoint(cx)?;
        validate_ray(origin, unit_direction, max_distance_m)?;
        let max_steps = self
            .medium
            .max_steps
            .min(u32::try_from(max_cell_visits).unwrap_or(u32::MAX));
        let mut result = LuminescenceTransport {
            transmittance: [1.0; M],
            emitted_radiance_w_m2_sr: [0.0; M],
        };
        let Some(mut ray) = CellRay::clip(
            &self.medium.grid,
            origin,
            unit_direction,
            max_distance_m,
            max_steps,
        )?
        else {
            return Ok(LuminescenceRaySample {
                transport: result,
                optical_depth: [0.0; M],
                cell_visits: 0,
            });
        };
        let mut kappa = self.medium.coefficients_at_cell(0, &self.wavelengths_nm)?;
        let mut previous_model = None;
        let mut depth = [0.0; M];
        let mut radiance: [PositiveSum; M] = std::array::from_fn(|_| PositiveSum::default());
        let mut cell_visits = 0;
        while let Some(segment) = ray.next(cx)? {
            // CellRay already charged this bounded u32 visit, including
            // zero-length face crossings and transparent emissive cells.
            cell_visits += 1;
            if segment.length_m == 0.0 {
                continue;
            }
            let density = self.medium.grid.data[segment.cell];
            if let Some(palette) = &self.medium.absorption_palette {
                let model = palette.cell_indices()[segment.cell];
                if previous_model != Some(model) {
                    kappa = self
                        .medium
                        .coefficients_at_cell(segment.cell, &self.wavelengths_nm)?;
                    previous_model = Some(model);
                }
            }
            let source = self
                .emissivity
                .binary_search_by_key(&segment.cell, |(cell, _)| *cell)
                .ok()
                .map(|index| &self.emissivity[index].1);
            let additional = self
                .additional_absorption
                .binary_search_by_key(&segment.cell, |(cell, _)| *cell)
                .ok()
                .map(|index| &self.additional_absorption[index].1);
            for lane in 0..M {
                let host_tau = if density == 0.0 || kappa[lane] == 0.0 {
                    0.0
                } else {
                    positive_product(density, kappa[lane], segment.length_m)
                };
                let extra = additional.map_or(0.0, |values| values[lane]);
                let tau = if extra == 0.0 {
                    host_tau
                } else {
                    host_tau + positive_product(extra, segment.length_m, 1.0)
                };
                if let Some(source) = source {
                    let value = if extra == 0.0 {
                        segment_emission(
                            source[lane],
                            density,
                            kappa[lane],
                            segment.length_m,
                            tau,
                            depth[lane],
                        )
                    } else {
                        segment_emission_integral(
                            source[lane],
                            segment.length_m,
                            tau,
                            depth[lane],
                            || -log_total_linear_extinction(density, kappa[lane], extra),
                        )
                    };
                    radiance[lane].add(value)?;
                }
                depth[lane] += tau;
            }
        }
        for (lane, value) in radiance.into_iter().enumerate() {
            result.transmittance[lane] = det::exp(-depth[lane]);
            result.emitted_radiance_w_m2_sr[lane] = value.finish()?;
        }
        checkpoint(cx)?;
        Ok(LuminescenceRaySample {
            transport: result,
            optical_depth: depth,
            cell_visits,
        })
    }
}

fn segment_emission(
    j: f64,
    density: f64,
    kappa: f64,
    length: f64,
    tau: f64,
    foreground: f64,
) -> f64 {
    segment_emission_integral(j, length, tau, foreground, || {
        -det::ln(density) - det::ln(kappa)
    })
}

fn segment_emission_integral(
    j: f64,
    length: f64,
    tau: f64,
    foreground: f64,
    opaque_log_length: impl FnOnce() -> f64,
) -> f64 {
    if j == 0.0 || foreground == f64::INFINITY {
        return 0.0;
    }
    let fraction = if tau == 0.0 {
        1.0
    } else {
        -det::expm1(-tau) / tau
    };
    let unattenuated = j * length;
    let transmitted = det::exp(-foreground);
    let value = unattenuated * fraction * transmitted;
    if unattenuated.is_finite()
        && unattenuated >= f64::MIN_POSITIVE
        && value.is_finite()
        && value >= f64::MIN_POSITIVE
    {
        return value;
    }
    // A very bright source can remain visible behind a transmittance that
    // underflows on its own. Opaque-cell optical depth may overflow while j/alpha
    // remains finite. Evaluate those rare cases before rounding their factors.
    let log_length = if tau == 0.0 {
        det::ln(length)
    } else if tau.is_infinite() {
        opaque_log_length()
    } else {
        det::ln(length) + det::ln(-det::expm1(-tau)) - det::ln(tau)
    };
    det::exp(det::ln(j) + log_length - foreground)
}
