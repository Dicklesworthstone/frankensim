//! Prepared physical-grid extinction for the existing direct-lux enclosure.
//!
//! Coefficients are prepared once at the shared photometric wavelengths. Each
//! uniform-optics ray computes one column mass; indexed optics integrate each
//! original cell spectrum with the same bounded traversal. Original
//! spectra and retained observer weights are integrated BEFORE conversion to
//! lux; no RGB fit or wavelength-independent haze substitute is introduced.

use super::{finite, hash_float, nonnegative};
use crate::photometry::source::{
    LUMINOUS_EFFICACY_LM_W, SPECTRAL_BINS, PhotometricSource,
    photopic_weight_nm, wavelength_nm,
};
use crate::photometry::transport::PhotometricMedium;
use crate::volumes::transport::VolumeTransportError;
use crate::volumes::transport::enclosure::{
    GRID_EXTINCTION_ENCLOSURE_SEMANTICS_VERSION, GridEnclosureBudget,
    MAX_PREPARED_EXTINCTION_VALUES, PreparedGridExtinction,
};
use fs_blake3::DomainHasher;
use fs_exec::Cx;
use fs_geom::Point3;
use fs_ivl::{Interval, decay::exp_neg};

pub(super) enum PreparedMedium<'field> {
    Gray { extinction: Interval, vacuum: bool },
    Grid { packet: PreparedGridExtinction<'field, 'field, SPECTRAL_BINS>, budget: GridEnclosureBudget },
}

impl<'field> PreparedMedium<'field> {
    pub(super) fn prepare(
        medium: &PhotometricMedium<'field>, hash: &mut DomainHasher, cx: &Cx<'_>,
    ) -> Result<Self, VolumeTransportError> {
        let wavelengths: [f64; SPECTRAL_BINS] = std::array::from_fn(wavelength_nm);
        let packet = match medium {
            PhotometricMedium::Homogeneous(medium) => {
                hash.update(&[0]);
                hash_float(hash, medium.absorption_per_m());
                hash_float(hash, medium.scattering_per_m());
                return Ok(Self::Gray {
                    extinction: nonnegative(Interval::point(medium.absorption_per_m())
                        + Interval::point(medium.scattering_per_m())),
                    vacuum: medium.extinction_per_m() == 0.0,
                });
            }
            PhotometricMedium::Absorbing(grid) => {
                hash.update(&[1]);
                hash.update(grid.content_identity().as_bytes());
                (*grid).prepare_extinction_enclosure(cx, &wavelengths, MAX_PREPARED_EXTINCTION_VALUES)?
            }
            PhotometricMedium::Scattering(grid) => {
                hash.update(&[2]);
                hash.update(grid.content_identity().as_bytes());
                (*grid).prepare_extinction_enclosure(cx, &wavelengths, MAX_PREPARED_EXTINCTION_VALUES)?
            }
        };
        hash.update(&GRID_EXTINCTION_ENCLOSURE_SEMANTICS_VERSION.to_le_bytes());
        let observer = crate::photometry::observer::OBSERVER_ID.as_bytes();
        hash.update(&(observer.len() as u64).to_le_bytes());
        hash.update(observer);
        hash_float(hash, LUMINOUS_EFFICACY_LM_W);
        hash.update(&(SPECTRAL_BINS as u64).to_le_bytes());
        for (bin, wavelength) in wavelengths.into_iter().enumerate() {
            hash_float(hash, wavelength);
            hash_float(hash, photopic_weight_nm(bin));
        }
        Ok(Self::Grid { packet, budget: GridEnclosureBudget::default() })
    }

    pub(super) fn needs_spectrum(&self) -> bool { matches!(self, Self::Grid { .. }) }

    pub(super) fn set_grid_budget(&mut self, new_budget: GridEnclosureBudget) {
        if let Self::Grid { budget, .. } = self { *budget = new_budget; }
    }

    pub(super) fn photopic_transfer(
        &self,
        cx: &Cx<'_>,
        source: &PhotometricSource,
        receiver: Point3,
        distance: Interval,
    ) -> Result<Interval, VolumeTransportError> {
        let (packet, budget) = match self {
            Self::Gray { extinction, vacuum } => {
                return if *vacuum { Ok(Interval::point(1.0)) } else {
                    exp_neg(nonnegative(*extinction * distance))
                        .ok_or(VolumeTransportError::Unrepresentable)
                };
            }
            Self::Grid { packet, budget } => (packet, *budget),
        };
        let transmission = packet.transmittance(cx, source.position_m(), receiver, budget)?;
        // Retained normalization need not equal one in real arithmetic. Bound
        // its actual weighted sum rather than inventing perfect normalization.
        // The enclosing caller also retains the unchanged nominal gray fast
        // path, so any discrepancy with that scalar convention widens the result.
        let mut sum = Interval::point(0.0);
        for (bin, (&power, transmission)) in source.spectrum().bins_w_per_lumen_nm()
            .iter().zip(transmission).enumerate() {
            let weight = photopic_weight_nm(bin);
            if power == 0.0 || weight == 0.0 || transmission.hi() == 0.0 { continue; }
            sum = nonnegative(sum + Interval::point(power) * Interval::point(weight) * transmission);
        }
        let result = nonnegative(Interval::point(LUMINOUS_EFFICACY_LM_W) * sum);
        if !finite(result) { return Err(VolumeTransportError::Unrepresentable); }
        Ok(result)
    }
}

#[cfg(test)]
mod tests;
