//! Degenerate two-photon absorption of one collimated monochromatic beam.
//!
//! Each native cell solves `dI/ds = -alpha I - beta I²`, with irradiance
//! `I = P / area`, `alpha = rho * kappa_a` in m⁻¹ and bulk `beta` in m/W.
//! The exact positive solution is `P_out = P exp(-alpha L)/(1 + q)`, where
//! `q = beta I L_eff` and `L_eff = -expm1(-alpha L)/alpha`. The zero-alpha
//! limit is `L_eff = L`. See Sheik-Bahae et al., IEEE JQE 26 (1990),
//! equations (23)–(25), <https://doi.org/10.1109/3.53394>.
//!
//! Area is constant and irradiance uniform across the beam. The caller binds
//! a passive channel that thermalizes absorbed energy on the thermal clock;
//! no carrier/excitation storage, fluorescence, saturation, nonlinear
//! refraction, diffraction or cross-frequency process is inferred. A pulse
//! must be resolved at its instantaneous power, not its mean power. Overlapping
//! beams cannot generally be deposited separately and added: this law is
//! nonlinear. Linear transport APIs and their superposition remain unchanged.

use super::{
    AbsorbingGrid, BeamPowerDeposition, CellBeamPower, CellRay, ContentHash, Cx, Point3, Vec3,
    VolumeTransportError, absorb_beam, add_compensated, checkpoint, det, positive_product,
    validate_ray,
};
use crate::volumes::transport::CellPalette;
use fs_blake3::DomainHasher;

/// A source-declared constant bulk beta over an explicit vacuum-wavelength
/// interval and irradiance ceiling. Beta already includes both absorbed
/// photons; multiplying the loss by two would double count energy. It is not
/// a molecular cross section or a mass coefficient. Its source must cover the
/// original material density, state and polarization; no density scaling is
/// inferred. Zero-density cells remain vacuum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TwoPhotonAbsorption {
    beta_m_w: f64,
    wavelength_bounds_nm: (f64, f64),
    maximum_irradiance_w_m2: f64,
    source_identity: ContentHash,
}

/// No partial beam report is returned for any refused source, ray or budget.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum TwoPhotonError {
    /// Invalid beta, wavelength interval, irradiance ceiling or source identity.
    Source,
    /// Beam cross-sectional area must be positive and finite, in square metres.
    BeamArea,
    /// Local incident irradiance exceeds this original cell's source support.
    Irradiance {
        /// Original native cell whose incident irradiance exceeds source support.
        cell_index: usize,
    },
    /// Original geometry, spectrum, cancellation or output-budget refusal.
    Transport(VolumeTransportError),
}

impl From<VolumeTransportError> for TwoPhotonError {
    fn from(value: VolumeTransportError) -> Self {
        Self::Transport(value)
    }
}
impl core::fmt::Display for TwoPhotonError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "two-photon absorption refused: {self:?}")
    }
}
impl core::error::Error for TwoPhotonError {}

impl TwoPhotonAbsorption {
    /// Declare a nonnegative bulk coefficient and its measured source support.
    pub fn try_new(
        beta_m_w: f64,
        wavelength_bounds_nm: (f64, f64),
        maximum_irradiance_w_m2: f64,
        source_identity: ContentHash,
    ) -> Result<Self, TwoPhotonError> {
        let (lo, hi) = wavelength_bounds_nm;
        if !beta_m_w.is_finite()
            || beta_m_w < 0.0
            || !lo.is_finite()
            || !hi.is_finite()
            || lo <= 0.0
            || hi <= lo
            || !maximum_irradiance_w_m2.is_finite()
            || maximum_irradiance_w_m2 <= 0.0
            || source_identity == ContentHash([0; 32])
        {
            return Err(TwoPhotonError::Source);
        }
        Ok(Self {
            beta_m_w,
            wavelength_bounds_nm,
            maximum_irradiance_w_m2,
            source_identity,
        })
    }
}

/// Prepared nonlinear absorption over the original borrowed physical grid.
/// Uniform and indexed laws share the same DDA and closed-form cell solve.
/// Preparation hashes the optical models once; each query visits only crossed
/// cells and allocates only its bounded deposition output.
pub struct TwoPhotonAbsorbingGrid<'m, 'g> {
    medium: &'m AbsorbingGrid<'g>,
    uniform: TwoPhotonAbsorption,
    palette: Option<CellPalette<TwoPhotonAbsorption>>,
    identity: ContentHash,
    bounds: (f64, f64),
}

impl<'m, 'g> TwoPhotonAbsorbingGrid<'m, 'g> {
    /// Bind a uniform bulk beta to every positive-density native cell.
    pub fn try_new(
        cx: &Cx<'_>,
        medium: &'m AbsorbingGrid<'g>,
        law: TwoPhotonAbsorption,
    ) -> Result<Self, TwoPhotonError> {
        Self::prepare(cx, medium, law, None)
    }

    /// Bind independently sourced nonlinear laws in native x-major cell order.
    /// The existing palette admission bounds model storage and index scanning.
    /// All supplied source bands must have a nonempty common interval with
    /// the original linear absorption. No source is extrapolated.
    pub fn try_new_indexed(
        cx: &Cx<'_>,
        medium: &'m AbsorbingGrid<'g>,
        laws: CellPalette<TwoPhotonAbsorption>,
    ) -> Result<Self, TwoPhotonError> {
        checkpoint(cx)?;
        if laws.cell_indices().len() != medium.grid.data.len() {
            return Err(VolumeTransportError::Palette.into());
        }
        Self::prepare(cx, medium, laws.values()[0], Some(laws))
    }

    fn prepare(
        cx: &Cx<'_>,
        medium: &'m AbsorbingGrid<'g>,
        uniform: TwoPhotonAbsorption,
        palette: Option<CellPalette<TwoPhotonAbsorption>>,
    ) -> Result<Self, TwoPhotonError> {
        checkpoint(cx)?;
        let mut bounds = medium.wavelength_bounds_nm();
        let mut hash = DomainHasher::new("org.frankensim.render.two-photon-grid.v1");
        hash.update(&medium.identity.0);
        hash.update(&[u8::from(palette.is_some())]);
        let values = palette
            .as_ref()
            .map_or(core::slice::from_ref(&uniform), |p| p.values());
        hash.update(&(values.len() as u64).to_le_bytes());
        for chunk in values.chunks(256) {
            checkpoint(cx)?;
            for law in chunk {
                bounds.0 = bounds.0.max(law.wavelength_bounds_nm.0);
                bounds.1 = bounds.1.min(law.wavelength_bounds_nm.1);
                for value in [
                    law.beta_m_w,
                    law.wavelength_bounds_nm.0,
                    law.wavelength_bounds_nm.1,
                    law.maximum_irradiance_w_m2,
                ] {
                    hash.update(&value.to_bits().to_le_bytes());
                }
                hash.update(&law.source_identity.0);
            }
        }
        if bounds.0 >= bounds.1 {
            return Err(VolumeTransportError::Wavelength.into());
        }
        if let Some(palette) = &palette {
            for chunk in palette.cell_indices().chunks(1024) {
                checkpoint(cx)?;
                for index in chunk {
                    hash.update(&index.to_le_bytes());
                }
            }
        }
        checkpoint(cx)?;
        Ok(Self {
            medium,
            uniform,
            palette,
            identity: hash.finalize(),
            bounds,
        })
    }

    /// Identity of the borrowed medium, nonlinear laws and original cell map.
    #[must_use]
    pub const fn content_identity(&self) -> ContentHash {
        self.identity
    }

    /// Common supported wavelength interval of the linear and nonlinear sources.
    #[must_use]
    pub const fn wavelength_bounds_nm(&self) -> (f64, f64) {
        self.bounds
    }

    /// Deposit one instantaneous monochromatic beam, in watts, over an explicit
    /// constant cross-sectional area in m². Geometry must be uniform over that
    /// cross section. This is power already inside the medium; no entrance
    /// Fresnel reflection is invented. The result feeds the existing
    /// `OpticalHeating` observer directly when prompt thermalization applies.
    ///
    /// Only absorption is present, so scattered power is exactly zero. A local
    /// irradiance-domain failure, cancellation or exhausted cell/output budget
    /// refuses the entire result. Per-cell alpha comes from the original linear
    /// palette. Each bulk beta is applied as authored, not multiplied by rho.
    #[allow(clippy::too_many_arguments)]
    pub fn deposit_beam_power(
        &self,
        cx: &Cx<'_>,
        origin: Point3,
        unit_direction: Vec3,
        max_distance_m: f64,
        wavelength_nm: f64,
        incident_power_w: f64,
        beam_area_m2: f64,
        max_deposition_cells: usize,
    ) -> Result<BeamPowerDeposition<1>, TwoPhotonError> {
        checkpoint(cx)?;
        validate_ray(origin, unit_direction, max_distance_m)?;
        if !wavelength_nm.is_finite() || !(self.bounds.0..=self.bounds.1).contains(&wavelength_nm) {
            return Err(VolumeTransportError::Wavelength.into());
        }
        if !incident_power_w.is_finite() || incident_power_w < 0.0 {
            return Err(VolumeTransportError::BeamPower { lane: 0 }.into());
        }
        if !beam_area_m2.is_finite() || beam_area_m2 <= 0.0 {
            return Err(TwoPhotonError::BeamArea);
        }
        let mut result = BeamPowerDeposition {
            medium_identity: self.identity,
            wavelength_bounds_nm: self.bounds,
            wavelengths_nm: [wavelength_nm],
            incident_power_w: [incident_power_w],
            transmitted_power_w: [incident_power_w],
            absorbed_power_w: [0.0],
            scattered_power_w: [0.0],
            cells: Vec::new(),
        };
        if incident_power_w == 0.0 {
            return Ok(result);
        }
        let Some(mut cells) = CellRay::clip(
            &self.medium.grid,
            origin,
            unit_direction,
            max_distance_m,
            self.medium.max_steps,
        )?
        else {
            return Ok(result);
        };
        let mut correction = 0.0;
        let mut last_linear = None;
        let mut kappa = 0.0;
        while let Some(cell) = cells.next(cx)? {
            let rho = self.medium.grid.data[cell.cell];
            let power = result.transmitted_power_w[0];
            if rho == 0.0 || cell.length_m == 0.0 || power == 0.0 {
                continue;
            }
            let law = self.palette.as_ref().map_or(&self.uniform, |p| {
                &p.values()[p.cell_indices()[cell.cell] as usize]
            });
            if power / beam_area_m2 > law.maximum_irradiance_w_m2 {
                return Err(TwoPhotonError::Irradiance {
                    cell_index: cell.cell,
                });
            }
            let linear = self
                .medium
                .absorption_palette()
                .map_or(0, |p| p.cell_indices()[cell.cell]);
            if last_linear != Some(linear) {
                kappa = self
                    .medium
                    .absorption_at_cell(cell.cell)?
                    .at(wavelength_nm)?;
                last_linear = Some(linear);
            }
            let (absorbed, transmitted) =
                absorb_two_photon(power, beam_area_m2, rho, cell.length_m, kappa, law.beta_m_w);
            if !absorbed.is_finite() || !transmitted.is_finite() {
                return Err(VolumeTransportError::Unrepresentable.into());
            }
            result.transmitted_power_w[0] = transmitted;
            if absorbed == 0.0 {
                continue;
            }
            if result.cells.len() == max_deposition_cells {
                return Err(VolumeTransportError::Budget.into());
            }
            if result.cells.len() == result.cells.capacity() {
                result
                    .cells
                    .try_reserve_exact(
                        result
                            .cells
                            .capacity()
                            .max(16)
                            .min(max_deposition_cells - result.cells.len()),
                    )
                    .map_err(|_| VolumeTransportError::Budget)?;
            }
            add_compensated(&mut result.absorbed_power_w[0], &mut correction, absorbed)?;
            result.cells.push(CellBeamPower {
                cell_index: cell.cell,
                absorbed_power_w: [absorbed],
                scattered_power_w: [0.0],
            });
        }
        checkpoint(cx)?;
        Ok(result)
    }
}

pub(crate) fn absorb_two_photon(
    power: f64,
    area: f64,
    rho: f64,
    length: f64,
    kappa: f64,
    beta: f64,
) -> (f64, f64) {
    if power == 0.0 {
        return (0.0, 0.0);
    }
    if beta == 0.0 {
        return absorb_beam(power, rho, length, kappa);
    }
    let tau = if kappa == 0.0 {
        0.0
    } else {
        positive_product(rho, kappa, length)
    };
    let opacity = -det::expm1(-tau);
    let phi = if tau == 0.0 { 1.0 } else { opacity / tau };
    let effective_length = length * phi;
    let irradiance = power / area;
    let rate = beta * irradiance;
    let mut q = rate * effective_length;
    // Avoid losing finite watts to an underflowed optical fraction or to an
    // overflowing intermediate product. Logarithms are only the range fallback.
    let mut log_q = None;
    if !irradiance.is_normal()
        || !rate.is_normal()
        || !phi.is_normal()
        || !effective_length.is_normal()
        || !q.is_normal()
    {
        let log_length = if tau == 0.0 {
            det::ln(length)
        } else if tau.is_finite() {
            det::ln(length) + det::ln(opacity) - det::ln(tau)
        } else {
            -det::ln(rho) - det::ln(kappa)
        };
        let value = det::ln(beta) + det::ln(power) - det::ln(area) + log_length;
        q = det::exp(value);
        log_q = Some(value);
    }
    let inverse = 1.0 / (1.0 + q);
    let fraction = if q <= 1.0 {
        (opacity + q) * inverse
    } else {
        (1.0 + opacity / q) / (1.0 + 1.0 / q)
    };
    let log_denominator = || {
        if q.is_finite() {
            det::ln(1.0 + q)
        } else {
            log_q.unwrap_or_else(|| det::ln(q))
        }
    };
    let absorbed = if fraction >= f64::MIN_POSITIVE {
        power * fraction
    } else {
        let log_linear = if opacity > 0.0 {
            det::ln(opacity)
        } else if kappa == 0.0 {
            f64::NEG_INFINITY
        } else {
            det::ln(rho) + det::ln(kappa) + det::ln(length)
        };
        let log_nonlinear = log_q.unwrap_or_else(|| det::ln(q));
        let hi = log_linear.max(log_nonlinear);
        let log_removed = hi + det::ln(det::exp(log_linear - hi) + det::exp(log_nonlinear - hi));
        det::exp(det::ln(power) + log_removed - log_denominator())
    };
    let transmittance = det::exp(-tau);
    let transmitted = if transmittance >= f64::MIN_POSITIVE && inverse >= f64::MIN_POSITIVE {
        (power * transmittance) * inverse
    } else {
        det::exp(det::ln(power) - tau - log_denominator())
    };
    (absorbed, transmitted)
}

#[cfg(test)]
mod tests;
