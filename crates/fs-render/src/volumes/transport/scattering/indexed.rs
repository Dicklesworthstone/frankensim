//! Cell-selected spectra and angular laws on the original physical grid.

use super::*;

/// One prepared scattering material. The mass coefficient and normalized
/// angular response belong to the same declared optical state.
#[derive(Clone, Debug, PartialEq)]
pub struct ScatteringOptics {
    /// Elastic scattering per kilogram of the density field's material.
    pub scattering: MassScatteringSpectrum,
    /// Reciprocal scalar angular response, including spectral dependence.
    pub phase: PhaseFunction,
}

impl<'a> ScatteringGrid<'a> {
    /// Bind distinct prepared scattering materials to physical cells.
    ///
    /// Absorption may independently be global or indexed. Both cell maps use
    /// the original grid's x-major order; their palette numbers need not agree.
    /// Every model must cover the absorption owner's complete wavelength band.
    /// Density multiplies each cell's own coefficients. No constitutive solve,
    /// phase-table copy, or palette allocation occurs during ray traversal.
    pub fn try_new_indexed(
        cx: &Cx<'_>,
        absorption: AbsorbingGrid<'a>,
        palette: CellPalette<ScatteringOptics>,
    ) -> Result<Self, VolumeTransportError> {
        checkpoint(cx)?;
        if palette.cell_indices().len() != absorption.grid.data.len() {
            return Err(VolumeTransportError::Palette);
        }
        let first = &palette.values()[0];
        Self::prepare(
            cx,
            absorption,
            first.scattering,
            first.phase.clone(),
            Some(palette),
        )
    }

    /// Complete immutable scattering model table and native-order cell map.
    #[must_use]
    pub const fn scattering_palette(&self) -> Option<&CellPalette<ScatteringOptics>> {
        self.scattering_palette.as_ref()
    }

    /// Angular law of the ORIGINAL collision cell. Use the index returned by
    /// `sample_segment`, avoiding a second position-to-cell classification at
    /// floating-point faces. Out-of-range cell access refuses.
    pub fn phase_at_cell(&self, cell: usize) -> Result<&PhaseFunction, VolumeTransportError> {
        self.check_cell(cell)?;
        Ok(self
            .scattering_palette
            .as_ref()
            .map_or(&self.phase, |palette| {
                &palette.at_cell(cell).expect("admitted cell map").phase
            }))
    }

    /// Exact mass spectrum selected by a physical cell, without interpolation
    /// between material identities or copied optical tables.
    pub fn scattering_at_cell(
        &self,
        cell: usize,
    ) -> Result<&MassScatteringSpectrum, VolumeTransportError> {
        self.check_cell(cell)?;
        Ok(self
            .scattering_palette
            .as_ref()
            .map_or(&self.scattering, |palette| {
                &palette.at_cell(cell).expect("admitted cell map").scattering
            }))
    }

    pub(in crate::volumes::transport) fn coefficient_key(
        &self,
        cell: usize,
    ) -> Result<(u32, u32), VolumeTransportError> {
        self.check_cell(cell)?;
        Ok((
            self.absorption.absorption_palette().map_or(0, |palette| {
                palette
                    .index_at_cell(cell)
                    .expect("admitted absorption map")
            }),
            self.scattering_palette.as_ref().map_or(0, |palette| {
                palette
                    .index_at_cell(cell)
                    .expect("admitted scattering map")
            }),
        ))
    }

    pub(in crate::volumes::transport) fn coefficients_at_cell<const N: usize>(
        &self,
        cell: usize,
        wavelengths: &[f64; N],
    ) -> Result<([f64; N], [f64; N], [f64; N]), VolumeTransportError> {
        if N == 0 {
            return Err(VolumeTransportError::Sample);
        }
        let a = self.absorption.absorption_at_cell(cell)?;
        let s = self.scattering_at_cell(cell)?;
        let (mut absorption, mut scattering, mut extinction) = ([0.0; N], [0.0; N], [0.0; N]);
        for lane in 0..N {
            absorption[lane] = a.at(wavelengths[lane])?;
            scattering[lane] = s.at(wavelengths[lane])?;
            extinction[lane] = absorption[lane] + scattering[lane];
            if !extinction[lane].is_finite() {
                return Err(VolumeTransportError::Unrepresentable);
            }
        }
        Ok((absorption, scattering, extinction))
    }

    fn check_cell(&self, cell: usize) -> Result<(), VolumeTransportError> {
        if cell >= self.absorption.grid.data.len() {
            Err(VolumeTransportError::Palette)
        } else {
            Ok(())
        }
    }

    pub(super) fn prepare(
        cx: &Cx<'_>,
        absorption: AbsorbingGrid<'a>,
        scattering: MassScatteringSpectrum,
        phase: PhaseFunction,
        scattering_palette: Option<CellPalette<ScatteringOptics>>,
    ) -> Result<Self, VolumeTransportError> {
        checkpoint(cx)?;
        let bounds = absorption.wavelength_bounds_nm();
        if let Some(palette) = &scattering_palette {
            for optical in palette.values() {
                checkpoint(cx)?;
                validate_model(&optical.scattering, &optical.phase, bounds)?;
            }
        } else {
            validate_model(&scattering, &phase, bounds)?;
        }
        let indexed = scattering_palette.is_some() || absorption.absorption_palette().is_some();
        let mut has_scattering = false;
        if indexed {
            let mut last_pair = None;
            let mut active = false;
            for (cell, &rho) in absorption.grid.data.iter().enumerate() {
                if cell % 1024 == 0 {
                    checkpoint(cx)?;
                }
                let pair = (
                    absorption
                        .absorption_palette()
                        .map_or(0, |p| p.cell_indices()[cell]),
                    scattering_palette
                        .as_ref()
                        .map_or(0, |p| p.cell_indices()[cell]),
                );
                if last_pair != Some(pair) {
                    let a = absorption.absorption_at_cell(cell)?;
                    let s = scattering_palette.as_ref().map_or(&scattering, |p| {
                        &p.at_cell(cell).expect("admitted cell map").scattering
                    });
                    validate_extinction(a, s)?;
                    active = !s.is_zero();
                    last_pair = Some(pair);
                }
                has_scattering |= rho > 0.0 && active;
            }
        } else {
            // Retain the existing uniform admission arithmetic and hash bytes.
            validate_extinction(&absorption.absorption, &scattering)?;
            if !scattering.is_zero() {
                for chunk in absorption.grid.data.chunks(1024) {
                    checkpoint(cx)?;
                    if chunk.iter().any(|density| *density > 0.0) {
                        has_scattering = true;
                        break;
                    }
                }
            }
        }
        let mut hash = DomainHasher::new(if scattering_palette.is_some() {
            "org.frankensim.render.indexed-scattering-grid.v1"
        } else {
            "org.frankensim.render.scattering-grid.v1"
        });
        hash.update(absorption.content_identity().as_bytes());
        if let Some(palette) = &scattering_palette {
            hash.update(&(palette.values().len() as u64).to_le_bytes());
            for optical in palette.values() {
                checkpoint(cx)?;
                hash_model(&optical.scattering, &optical.phase, &mut hash);
            }
            hash.update(&(palette.cell_indices().len() as u64).to_le_bytes());
            for chunk in palette.cell_indices().chunks(1024) {
                checkpoint(cx)?;
                for index in chunk {
                    hash.update(&index.to_le_bytes());
                }
            }
        } else {
            hash_model(&scattering, &phase, &mut hash);
        }
        checkpoint(cx)?;
        Ok(Self {
            absorption,
            scattering,
            phase,
            scattering_palette,
            identity: hash.finalize(),
            has_scattering,
        })
    }
}

fn validate_model(
    scattering: &MassScatteringSpectrum,
    phase: &PhaseFunction,
    bounds: (f64, f64),
) -> Result<(), VolumeTransportError> {
    phase.validate()?;
    if scattering.wavelength_bounds_nm() != bounds {
        return Err(VolumeTransportError::Wavelength);
    }
    phase.validate_wavelength(bounds.0)?;
    phase.validate_wavelength(bounds.1)
}

fn validate_extinction(
    absorption: &MassAbsorptionSpectrum,
    scattering: &MassScatteringSpectrum,
) -> Result<(), VolumeTransportError> {
    // The sum is convex on every combined affine / positive lambda^-4
    // interval. Endpoints therefore suffice for finite-extinction admission.
    for &(wavelength, _) in absorption.knots[..absorption.count]
        .iter()
        .chain(scattering.knots())
    {
        if !(absorption.at(wavelength)? + scattering.at(wavelength)?).is_finite() {
            return Err(VolumeTransportError::Unrepresentable);
        }
    }
    Ok(())
}

fn hash_model(scattering: &MassScatteringSpectrum, phase: &PhaseFunction, hash: &mut DomainHasher) {
    scattering.hash_into(hash);
    phase.hash_into(hash);
}

#[cfg(test)]
mod tests;
