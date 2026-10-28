//! Entropy observations of the same source-bound gray/diffuse emitting surface.
use super::GrayDiffuseThermalSurface;
use fs_conduction::radiation_entropy::{
    GrayPhotonEmission, PhotonEmissionFlux, PhotonEntropyQuadrature, RadiationEntropyError,
};

impl GrayDiffuseThermalSurface {
    /// Shared conduction-layer photon closure at this immutable surface's actual
    /// temperature and emissivity. Reuse this value for a spectral sweep.
    /// The surface retains all source identity; this does not authenticate data.
    ///
    /// The closure describes the emitted component into vacuum. A general
    /// illuminated surface requires the outgoing mode occupations including
    /// reflection; adding emitted and reflected entropy separately is not valid.
    pub fn photon_emission(self) -> Result<GrayPhotonEmission, RadiationEntropyError> {
        GrayPhotonEmission::new(self.temperature_k(), self.emissivity())
    }

    /// Physical emitted entropy radiance, W/(m² sr K nm), before display mapping.
    /// This performs no tracer sampling, material query, or physical state update.
    pub fn entropy_radiance_w_m2_sr_k_nm(self, wavelength_nm: f64)
        -> Result<f64, RadiationEntropyError>
    {
        Ok(self.photon_emission()?.spectrum_nm(wavelength_nm)?.entropy_w_m2_sr_k_nm)
    }

    /// Complete-spectrum emitted energy and photon entropy per area of one side.
    /// The existing visible-band image and bolometric cooling debit are unchanged.
    /// No extra heat is removed. Budget/cancellation errors retain their original
    /// type; the caller may cache this immutable result for repeated observations.
    pub fn photon_emission_flux(self, controls: PhotonEntropyQuadrature,
        cancelled: impl FnMut() -> bool) -> Result<PhotonEmissionFlux, RadiationEntropyError>
    {
        self.photon_emission()?.bolometric(controls, cancelled)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_blake3::ContentHash;

    fn surface(temperature: f64, emissivity: f64) -> GrayDiffuseThermalSurface {
        GrayDiffuseThermalSurface::try_new(temperature, emissivity, ContentHash([7; 32])).unwrap()
    }
    fn close(a: f64, b: f64) {
        assert!((a-b).abs() <= 2e-12*b.abs(), "{a:e} != {b:e}");
    }

    #[test]
    fn original_surface_and_shared_photon_closure_observe_the_same_state() {
        for t in [300.0, 1000.0, 5778.0] {
            for e in [0.0, 0.1, 0.9, 1.0] {
                let s = surface(t, e);
                let saved = s;
                let identity = s.identity();
                let flux = s.photon_emission_flux(PhotonEntropyQuadrature::default(), || false).unwrap();
                close(flux.energy_w_m2, s.bolometric_exitance_w_m2());
                for wavelength in [400.0, 1000.0, 10_000.0] {
                    let spectrum = s.photon_emission().unwrap().spectrum_nm(wavelength).unwrap();
                    close(spectrum.energy_w_m2_sr_nm, s.radiance_w_m2_sr_nm(wavelength).unwrap());
                    assert_eq!(s.entropy_radiance_w_m2_sr_k_nm(wavelength).unwrap(), spectrum.entropy_w_m2_sr_k_nm);
                }
                assert_eq!(s, saved);
                assert_eq!(s.identity(), identity);
                assert_eq!(s.reflectance(), 1.0-e);
            }
        }
    }

    #[test]
    fn gray_entropy_is_not_a_color_band_or_emissivity_scaled_blackbody() {
        let black = surface(400.0, 1.0).photon_emission_flux(PhotonEntropyQuadrature::default(), || false).unwrap();
        let gray = surface(400.0, 0.1).photon_emission_flux(PhotonEntropyQuadrature::default(), || false).unwrap();
        close(black.entropy_w_m2_k, (4.0/3.0)*black.energy_w_m2/400.0);
        close(gray.energy_w_m2, 0.1*black.energy_w_m2);
        assert!(gray.entropy_w_m2_k > 0.1*black.entropy_w_m2_k);
        assert!(gray.entropy_w_m2_k < black.entropy_w_m2_k);
        let dark = surface(400.0, 0.0).photon_emission_flux(PhotonEntropyQuadrature::default(), || false).unwrap();
        assert_eq!(dark.entropy_w_m2_k, 0.0);
        assert_eq!(dark.evaluations, 0);
    }

    #[test]
    fn observation_refusal_does_not_change_radiance_or_source_identity() {
        let s = surface(500.0, 0.2);
        let saved = (s.identity(), s.radiance_w_m2_sr_nm(500.0).unwrap());
        assert_eq!(s.photon_emission_flux(PhotonEntropyQuadrature::default(), || true), Err(RadiationEntropyError::Cancelled));
        let tiny = PhotonEntropyQuadrature { relative_tolerance: 1e-12, maximum_evaluations: 33 };
        assert!(matches!(s.photon_emission_flux(tiny, || false), Err(RadiationEntropyError::BudgetExceeded { .. })));
        assert_eq!(s.entropy_radiance_w_m2_sr_k_nm(f64::NAN), Err(RadiationEntropyError::Wavelength));
        assert_eq!(saved, (s.identity(), s.radiance_w_m2_sr_nm(500.0).unwrap()));
        assert!(s.photon_emission_flux(PhotonEntropyQuadrature::default(), || false).is_ok());
    }
}
