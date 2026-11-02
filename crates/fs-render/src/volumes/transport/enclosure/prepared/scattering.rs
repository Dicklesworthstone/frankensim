//! Column-integrated scattering response of the original cells and phase laws.
//! This is NOT a radiance solution: callers still supply geometric source
//! irradiance and both transmission legs, retaining each wavelength until the
//! final observable integration. Thermal emission and higher orders are absent.

use super::*;
use crate::volumes::transport::scattering::PhaseFunction;

mod remainder;
pub use remainder::{VolumeContractionEnvelope, VOLUME_CONTRACTION_ENCLOSURE_VERSION};

/// Prepared original extinction plus separately retained scattering response.
/// No density, cell map, angular table or Cartesian palette product is copied.
/// One extra inline N-lane coefficient packet is needed only for uniform optics;
/// indexed optics reuse the scattering table already in PreparedGridExtinction.
pub struct PreparedGridScattering<'g, 'field, const N: usize> {
    extinction: PreparedGridExtinction<'g, 'field, N>,
    uniform_scattering: Option<[Interval; N]>,
    active: bool,
}

impl<'g, 'field, const N: usize> PreparedGridExtinction<'g, 'field, N> {
    /// Add native scattering/phase response without repeating extinction setup.
    /// The cap covers retained extinction and any extra uniform scattering
    /// intervals. Wavelengths remain the SAME packet supplied at preparation;
    /// they cannot be rebound to inconsistent phase wavelengths afterwards.
    /// Absorption-only fields are admitted as exactly zero single scattering.
    pub fn into_single_scattering(self, cx: &Cx<'_>, max_coefficient_values: usize)
        -> Result<PreparedGridScattering<'g, 'field, N>, VolumeTransportError> {
        checkpoint(cx)?;
        let count = match &self.coefficients {
            Coefficients::Uniform(_) => N,
            Coefficients::Indexed { absorption, scattering } => absorption.len()
                .checked_add(scattering.len()).and_then(|n| n.checked_mul(N))
                .ok_or(VolumeTransportError::Budget)?,
        };
        let extra = if self.scattering.is_some() && matches!(&self.coefficients, Coefficients::Uniform(_)) { N } else { 0 };
        if count.checked_add(extra).is_none_or(|n| n > max_coefficient_values) {
            return Err(VolumeTransportError::Budget);
        }
        let active = self.scattering.is_some_and(ScatteringGrid::has_scattering);
        let uniform_scattering = if extra > 0 {
            let medium = self.scattering.ok_or(VolumeTransportError::ScatteringSpectrum)?;
            Some(scattering_packet(medium.scattering_at_cell(0)?, &self.wavelengths_nm)?)
        } else { None };
        // Validate all declared angular domains once. No wavelength extrapolation
        // is allowed, even for a model currently unused by occupied cells.
        if let Some(medium) = self.scattering {
            if let Some(palette) = medium.scattering_palette() {
                for model in palette.values() {
                    checkpoint(cx)?;
                    for &wavelength in &self.wavelengths_nm {
                        model.phase.validate_wavelength(wavelength)?;
                    }
                }
            } else {
                let phase = medium.phase_at_cell(0)?;
                for &wavelength in &self.wavelengths_nm { phase.validate_wavelength(wavelength)?; }
            }
        }
        checkpoint(cx)?;
        Ok(PreparedGridScattering { extinction: self, uniform_scattering, active })
    }
}

impl<'g, 'field, const N: usize> PreparedGridScattering<'g, 'field, N> {
    /// Same original extinction owner; point and complete-ray-family queries
    /// retain their original cell/node budgets and numerical operation ordering.
    #[must_use]
    pub const fn extinction(&self) -> &PreparedGridExtinction<'g, 'field, N> { &self.extinction }

    /// Whether any occupied native cell has positive scattering support.
    #[must_use]
    pub const fn has_scattering(&self) -> bool { self.active }

    /// Enclose integral rho*kappa_s(lambda)*phase(lambda,mu) ds, in sr^-1,
    /// along every segment joining the endpoint boxes. `cosine` must contain
    /// every actual propagation angle on those segments. All possible cell
    /// ownership and column lengths come from the existing grid traversal;
    /// nonuniform density is integrated, never sampled at a panel midpoint.
    ///
    /// Per-cell phases remain coupled to their own scattering coefficient.
    /// A one-model cache avoids reevaluating phase packets across consecutive
    /// cells using the same palette entry. Queries allocate nothing. Failure
    /// or cancellation returns no partial response. This does not bound the
    /// incident field, transmittance, physical coefficients, or discretization.
    pub fn scattering_integral_box(
        &self, cx: &Cx<'_>, start: [Interval; 3], end: [Interval; 3],
        cosine: Interval, budget: GridEnclosureBudget,
    ) -> Result<[Interval; N], VolumeTransportError> {
        checkpoint(cx)?;
        if !finite(cosine) || cosine.lo() < -1.0 || cosine.hi() > 1.0 {
            return Err(VolumeTransportError::Phase);
        }
        // Even an analytically zero field retains input and traversal admission.
        let mut result = [Interval::point(0.0); N];
        let mut cached_index = None;
        let mut response = [Interval::point(0.0); N];
        let medium = self.extinction.scattering;
        self.extinction.grid.visit_box_column(cx, start, end, budget, |cell, column| {
            let Some(medium) = medium else { return Ok(()); };
            if !self.active || column.hi() == 0.0 { return Ok(()); }
            let index = medium.scattering_palette().map_or(0, |p| p.cell_indices()[cell] as usize);
            if cached_index != Some(index) {
                let scattering = match &self.extinction.coefficients {
                    Coefficients::Uniform(_) => self.uniform_scattering.as_ref()
                        .ok_or(VolumeTransportError::ScatteringSpectrum)?,
                    Coefficients::Indexed { scattering, .. } => &scattering[index],
                };
                let phase = medium.phase_at_cell(cell)?;
                let spectral = matches!(phase, PhaseFunction::SpectralTable(_) | PhaseFunction::SpectralMixture(_));
                let gray = if spectral { None } else { Some(phase.enclosure_cosine(cx, cosine)?) };
                for lane in 0..N {
                    checkpoint(cx)?;
                    response[lane] = if scattering[lane].hi() == 0.0 { Interval::point(0.0) } else {
                        let phase = match gray {
                            Some(value) => value,
                            None => phase.enclosure_at(cx, self.extinction.wavelengths_nm[lane], cosine)?,
                        };
                        nonnegative(scattering[lane] * phase)
                    };
                    if !finite(response[lane]) { return Err(VolumeTransportError::Unrepresentable); }
                }
                cached_index = Some(index);
            }
            for lane in 0..N {
                if lane % 64 == 0 { checkpoint(cx)?; }
                if response[lane].hi() == 0.0 { continue; }
                result[lane] = nonnegative(result[lane] + column * response[lane]);
                if !finite(result[lane]) { return Err(VolumeTransportError::Unrepresentable); }
            }
            Ok(())
        })?;
        checkpoint(cx)?;
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
