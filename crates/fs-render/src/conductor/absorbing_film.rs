//! A dispersive passive film with coherent absorption and reflection.
//!
//! The boundary-field characteristic matrix sums all internal amplitudes.
//! Its complex propagation entries are scaled by exp(-Im(delta)), which
//! cancels from reflection and avoids exponentially growing intermediates.
//! The exp(-i omega t), eta+i*k and polarization conventions match
//! `thin_film` and Byrnes, <https://arxiv.org/abs/1603.02720>.

use super::{
    CONDUCTOR_IOR_SAMPLE_COUNT, Complex, ComplexIor, ConductorDataStatus, ConductorError,
    ConductorIorSample, ConductorSource, MAX_COHERENT_COATING_THICKNESS_M,
    MAX_CONDUCTOR_KNOT_SPAN_NM, MAX_IOR_COMPONENT,
};
use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN};
use fs_blake3::DomainHasher;
use fs_math::det;

/// One visible-band absolute complex-index sample of a passive film.
/// Unlike an opaque substrate, a film can have exactly zero extinction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoatingIorSample(ConductorIorSample);

impl CoatingIorSample {
    /// Vacuum wavelength in nm, positive phase index and nonnegative
    /// dimensionless extinction coefficient k. Both index components are
    /// bounded by 100; negative k (gain) is refused.
    pub fn try_new(wavelength_nm: f64, eta: f64, k: f64) -> Result<Self, ConductorError> {
        if !wavelength_nm.is_finite()
            || !(LAMBDA_MIN..=LAMBDA_MAX).contains(&wavelength_nm)
            || !eta.is_finite()
            || eta <= 0.0
            || eta > MAX_IOR_COMPONENT
            || !k.is_finite()
            || !(0.0..=MAX_IOR_COMPONENT).contains(&k)
        {
            return Err(ConductorError::InvalidSample);
        }
        Ok(Self(ConductorIorSample {
            wavelength_nm,
            eta,
            k: super::canonical_zero(k),
        }))
    }

    /// Vacuum wavelength in nanometres.
    #[must_use]
    pub const fn wavelength_nm(self) -> f64 {
        self.0.wavelength_nm
    }

    /// Absolute real refractive index.
    #[must_use]
    pub const fn eta(self) -> f64 {
        self.0.eta
    }

    /// Dimensionless imaginary refractive index; zero is transparent.
    #[must_use]
    pub const fn k(self) -> f64 {
        self.0.k
    }
}

/// One uniform isotropic passive coating on an opaque substrate.
///
/// Nine explicit vacuum-wavelength knots span 380--780 nm. Both eta and k
/// interpolate linearly, retaining dispersion and absorption inside the
/// coherent phase calculation rather than applying an intensity tint after
/// interference. The sample/source declaration does not assert measurement
/// fidelity or enforce a broadband causal dispersion fit.
///
/// Thickness is normal to the local microfacet. As for the transparent film,
/// the caller declares locally parallel boundaries and sufficient coherence;
/// this model does not include arbitrary multilayers, scattering in the film,
/// thickness distributions, or coherence between distinct scene paths.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentAbsorbingCoating {
    samples: [ConductorIorSample; CONDUCTOR_IOR_SAMPLE_COUNT],
    thickness_m: f64,
    source: ConductorSource,
}

impl CoherentAbsorbingCoating {
    /// Admit a source-bound film. Knots must be strictly increasing, cover
    /// exactly 380--780 nm, and have no gap above 100 nm. Thickness is in
    /// metres, including exact zero, and is bounded by 100 micrometres.
    pub fn try_new(
        samples: [CoatingIorSample; CONDUCTOR_IOR_SAMPLE_COUNT],
        thickness_m: f64,
        source: ConductorSource,
    ) -> Result<Self, ConductorError> {
        if !thickness_m.is_finite()
            || !(0.0..=MAX_COHERENT_COATING_THICKNESS_M).contains(&thickness_m)
        {
            return Err(ConductorError::InvalidCoating);
        }
        let samples = samples.map(|sample| sample.0);
        if samples[0].wavelength_nm.to_bits() != LAMBDA_MIN.to_bits()
            || samples[CONDUCTOR_IOR_SAMPLE_COUNT - 1]
                .wavelength_nm
                .to_bits()
                != LAMBDA_MAX.to_bits()
            || samples.windows(2).any(|pair| {
                let span = pair[1].wavelength_nm - pair[0].wavelength_nm;
                span <= 0.0 || span > MAX_CONDUCTOR_KNOT_SPAN_NM
            })
        {
            return Err(ConductorError::InvalidTable);
        }
        Ok(Self {
            samples,
            thickness_m: super::canonical_zero(thickness_m),
            source,
        })
    }

    /// The film's absolute eta+i*k at the supplied vacuum wavelength.
    pub fn complex_ior(&self, wavelength_nm: f64) -> Result<ComplexIor, ConductorError> {
        if !wavelength_nm.is_finite() || !(LAMBDA_MIN..=LAMBDA_MAX).contains(&wavelength_nm) {
            return Err(ConductorError::InvalidWavelength);
        }
        Ok(super::interpolate_ior(&self.samples, wavelength_nm))
    }

    /// Exact admitted spectral knots in ascending vacuum wavelength.
    #[must_use]
    pub fn samples(&self) -> [CoatingIorSample; CONDUCTOR_IOR_SAMPLE_COUNT] {
        self.samples.map(CoatingIorSample)
    }

    /// Uniform physical thickness normal to the microfacet, in metres.
    #[must_use]
    pub const fn thickness_m(&self) -> f64 {
        self.thickness_m
    }

    /// Caller-declared optical/thickness source, without upgraded authority.
    #[must_use]
    pub const fn source(&self) -> ConductorSource {
        self.source
    }

    pub(super) fn update_identity(&self, hasher: &mut DomainHasher) {
        hasher.update(b"coherent-absorbing-coating.v1");
        for sample in self.samples {
            hasher.update(&sample.wavelength_nm.to_bits().to_le_bytes());
            hasher.update(&sample.eta.to_bits().to_le_bytes());
            hasher.update(&sample.k.to_bits().to_le_bytes());
        }
        hasher.update(&self.thickness_m.to_bits().to_le_bytes());
        hasher.update(self.source.identity().as_bytes());
        hasher.update(&[match self.source.status() {
            ConductorDataStatus::Representative => 0,
            ConductorDataStatus::CallerAssertedMeasured => 1,
            ConductorDataStatus::MaterialStateClaims => 2,
        }]);
    }

    // The common conductor owner already validated wavelength and substrate.
    pub(super) fn fresnel(
        &self,
        wavelength_nm: f64,
        incident_eta: f64,
        substrate: ComplexIor,
        cosine: f64,
    ) -> Result<f64, ConductorError> {
        if !incident_eta.is_finite() || incident_eta <= 0.0 || incident_eta > MAX_IOR_COMPONENT {
            return Err(ConductorError::InvalidIncidentMedium);
        }
        if !cosine.is_finite() || !(0.0..=1.0).contains(&cosine) {
            return Err(ConductorError::InvalidIncidence);
        }
        let film = super::interpolate_ior(&self.samples, wavelength_nm);
        if self.thickness_m == 0.0
            || (film.k == 0.0 && film.eta == incident_eta)
            || (film.eta == substrate.eta && film.k == substrate.k)
            || cosine == 0.0
        {
            return super::fresnel_conductor(incident_eta, substrate.eta, substrate.k, cosine);
        }

        let [r_s, r_p] = reflection_amplitudes(
            wavelength_nm, incident_eta, substrate, cosine, film, self.thickness_m,
        );
        let reflectance = 0.5 * (r_s.norm_squared() + r_p.norm_squared());
        if !reflectance.is_finite() || !(-1.0e-12..=1.0 + 1.0e-12).contains(&reflectance) {
            return Err(ConductorError::NumericalFailure);
        }
        Ok(reflectance.clamp(0.0, 1.0))
    }
}

// Shared by scalar power and polarization. Keep the native complex
// propagation and interface operation order; only the final readout differs.
pub(super) fn reflection_amplitudes(
    wavelength_nm: f64,
    incident_eta: f64,
    substrate: ComplexIor,
    cosine: f64,
    film: ComplexIor,
    thickness_m: f64,
) -> [Complex; 2] {
    let n0_squared = incident_eta * incident_eta;
    let transverse_squared = n0_squared * (1.0 - cosine * cosine);
    let n1_squared = Complex::new(film.eta, film.k).square();
    let q1_squared = n1_squared - Complex::new(transverse_squared, 0.0);
    let n2_squared = Complex::new(substrate.eta, substrate.k).square();
    let q2 = super::thin_film::forward_sqrt(n2_squared - Complex::new(transverse_squared, 0.0));
    let phase_scale = core::f64::consts::TAU * (thickness_m * 1.0e9) / wavelength_nm;
    let (cos_phase, sin_phase_over_q) = propagation(q1_squared, phase_scale);
    let q0 = incident_eta * cosine;
    let r_s = reflection(q0, q2, cos_phase, sin_phase_over_q, q1_squared * sin_phase_over_q);
    // The p wave uses magnetic-field admittance q/n² at both interfaces.
    let r_p = reflection(q0 / n0_squared, q2 / n2_squared, cos_phase,
        n1_squared * sin_phase_over_q, q1_squared * sin_phase_over_q / n1_squared);
    [r_s, r_p]
}

fn propagation(q_squared: Complex, phase_scale: f64) -> (Complex, Complex) {
    let q = super::thin_film::forward_sqrt(q_squared);
    let delta = q * phase_scale;
    if det::hypot(delta.re, delta.im) < 1.0e-4 {
        // Analytic critical-angle limit, also preserving arbitrarily weak
        // positive extinction: no 0/0 or subtraction of nearly equal roots.
        let square = delta.square();
        let fourth = square.square();
        return (
            Complex::new(1.0, 0.0) - square * 0.5 + fourth * (1.0 / 24.0),
            (Complex::new(1.0, 0.0) - square * (1.0 / 6.0) + fourth * (1.0 / 120.0)) * phase_scale,
        );
    }
    let (sine, cosine) = det::sin_cos(delta.re);
    let half_sum = 0.5 * (1.0 + det::exp(-2.0 * delta.im));
    let half_difference = -0.5 * det::expm1(-2.0 * delta.im);
    (
        Complex::new(cosine * half_sum, -sine * half_difference),
        Complex::new(sine * half_sum, cosine * half_difference) / q,
    )
}

fn reflection(
    incident_admittance: f64,
    substrate_admittance: Complex,
    cosine: Complex,
    sine_over_admittance: Complex,
    admittance_times_sine: Complex,
) -> Complex {
    let minus_i = Complex::new(0.0, -1.0);
    let field = cosine + minus_i * sine_over_admittance * substrate_admittance;
    let flux = substrate_admittance * cosine + minus_i * admittance_times_sine;
    (field * incident_admittance - flux) / (field * incident_admittance + flux)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conductor::{CoherentDielectricCoating, ConductorOptics};
    use crate::dielectric::CauchyIor;

    fn source() -> ConductorSource {
        ConductorSource::try_new(
            fs_blake3::hash_domain("org.frankensim.test.absorbing-film", b"synthetic"),
            ConductorDataStatus::Representative,
        )
        .unwrap()
    }

    fn film(eta: f64, k: f64, thickness_m: f64) -> CoherentAbsorbingCoating {
        CoherentAbsorbingCoating::try_new(
            core::array::from_fn(|i| {
                CoatingIorSample::try_new(380.0 + 50.0 * i as f64, eta, k).unwrap()
            }),
            thickness_m,
            source(),
        )
        .unwrap()
    }

    fn substrate(eta: f64, k: f64) -> ConductorOptics {
        ConductorOptics::try_new(
            core::array::from_fn(|i| {
                ConductorIorSample::try_new(380.0 + 50.0 * i as f64, eta, k).unwrap()
            }),
            source(),
        )
        .unwrap()
    }

    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
        );
    }

    #[test]
    fn g1_absorbing_film_matches_independent_complex_airy_answers() {
        // Python cmath oracle: r=(r01+r12*exp(2i*delta)) /
        // (1+r01*r12*exp(2i*delta)), independently for s and p. This differs
        // from the production boundary-field matrix. Fixtures include weak
        // loss at an interference null, evanescent incidence, a critical
        // film and an optically opaque film with more than 10,000 skin depths.
        for (n0, eta1, k1, eta2, k2, d_nm, wavelength, cosine, expected) in [
            (
                1.0,
                1.45,
                0.07,
                0.2,
                3.0,
                100.0,
                550.0,
                1.0,
                0.515_088_324_939_222_6,
            ),
            (
                1.33,
                2.2,
                0.35,
                3.0,
                4.0,
                77.0,
                430.0,
                0.6,
                0.134_219_925_249_046_88,
            ),
            (
                1.0,
                1.5,
                0.001,
                2.25,
                1.0e-8,
                100.0,
                600.0,
                1.0,
                6.335_829_231_524_68e-7,
            ),
            (
                1.5,
                1.0,
                0.01,
                0.2,
                3.0,
                80.0,
                580.0,
                0.55,
                0.744_887_377_519_035_4,
            ),
            (
                1.5,
                1.2,
                1.0e-10,
                0.2,
                3.0,
                120.0,
                600.0,
                0.6,
                0.862_195_634_034_291_3,
            ),
            (
                1.0,
                2.0,
                4.0,
                3.0,
                2.0,
                100_000.0,
                500.0,
                0.37,
                0.642_872_936_071_509_1,
            ),
            (
                1.0,
                1.8,
                0.6,
                0.17,
                3.1,
                35.0,
                780.0,
                0.81,
                0.483_650_748_473_711_8,
            ),
        ] {
            let actual = substrate(eta2, k2)
                .with_absorbing_coating(film(eta1, k1, d_nm * 1.0e-9))
                .fresnel(wavelength, n0, cosine)
                .unwrap();
            // The critical-angle Airy oracle loses several digits itself.
            let tolerance = if k1 == 1.0e-10 { 2.0e-11 } else { 5.0e-14 };
            close(actual, expected, tolerance);
        }
    }

    #[test]
    fn g3_transparent_limit_matches_existing_film_including_critical_angle() {
        let bare = substrate(0.2, 3.0);
        for (n0, n1) in [(1.0, 1.5), (1.5, 1.2)] {
            for thickness_m in [1.0e-14, 100.0e-9, 120.0e-9, 10.0e-6] {
                let transparent = bare.with_coating(
                    CoherentDielectricCoating::try_new(
                        CauchyIor::try_constant(n1).unwrap(),
                        thickness_m,
                        source(),
                    )
                    .unwrap(),
                );
                let passive = bare.with_absorbing_coating(film(n1, 0.0, thickness_m));
                for wavelength in [380.0, 470.0, 600.0, 780.0] {
                    for cosine in [0.01, 0.3, 0.6 - 1.0e-9, 0.6, 0.6 + 1.0e-9, 1.0] {
                        close(
                            passive.fresnel(wavelength, n0, cosine).unwrap(),
                            transparent.fresnel(wavelength, n0, cosine).unwrap(),
                            5.0e-13,
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn g3_absorbing_zero_and_matched_interfaces_keep_bare_bits() {
        let bare = substrate(0.2, 3.0);
        for coating in [
            film(1.5, 0.3, 0.0),
            film(1.33, 0.0, 80.0e-9),
            film(0.2, 3.0, 80.0e-9),
        ] {
            let coated = bare.with_absorbing_coating(coating);
            assert_ne!(bare.content_identity(), coated.content_identity());
            assert_eq!(coated.absorbing_coating(), Some(coating));
            assert_eq!(coated.coating(), None);
            for wavelength in [380.0, 555.0, 780.0] {
                for cosine in [0.0, 0.1, 0.5, 1.0] {
                    assert_eq!(
                        coated.fresnel(wavelength, 1.33, cosine).unwrap().to_bits(),
                        bare.fresnel(wavelength, 1.33, cosine).unwrap().to_bits(),
                    );
                }
            }
        }
    }

    #[test]
    fn g1_opaque_film_screens_substrate_and_preserves_its_front_interface() {
        let coating = film(1.8, 0.6, MAX_COHERENT_COATING_THICKNESS_M);
        for wavelength in [380.0, 530.0, 780.0] {
            for cosine in [0.01, 0.3, 0.8, 1.0] {
                let expected = super::super::fresnel_conductor(1.33, 1.8, 0.6, cosine).unwrap();
                for (eta, k) in [(0.2, 3.0), (4.0, 1.0)] {
                    close(
                        substrate(eta, k)
                            .with_absorbing_coating(coating)
                            .fresnel(wavelength, 1.33, cosine)
                            .unwrap(),
                        expected,
                        5.0e-14,
                    );
                }
            }
        }
    }

    #[test]
    fn g0_absorbing_spectrum_is_retained_and_bound_to_identity() {
        let samples = core::array::from_fn(|i| {
            CoatingIorSample::try_new(
                380.0 + 50.0 * i as f64,
                1.2 + 0.1 * i as f64,
                0.02 * i as f64,
            )
            .unwrap()
        });
        let coating = CoherentAbsorbingCoating::try_new(samples, 100.0e-9, source()).unwrap();
        assert_eq!(coating.samples(), samples);
        close(coating.complex_ior(505.0).unwrap().eta(), 1.45, 1.0e-15);
        close(coating.complex_ior(505.0).unwrap().k(), 0.05, 1.0e-15);
        assert_eq!(
            coating.complex_ior(379.0),
            Err(ConductorError::InvalidWavelength)
        );
        assert_eq!(
            coating.complex_ior(f64::NAN),
            Err(ConductorError::InvalidWavelength)
        );
        let bare = substrate(0.2, 3.0);
        let identity = bare.with_absorbing_coating(coating).content_identity();
        for i in 0..CONDUCTOR_IOR_SAMPLE_COUNT {
            let mut changed = samples;
            changed[i] = CoatingIorSample::try_new(
                samples[i].wavelength_nm(),
                samples[i].eta(),
                samples[i].k() + 0.001,
            )
            .unwrap();
            assert_ne!(
                identity,
                bare.with_absorbing_coating(
                    CoherentAbsorbingCoating::try_new(changed, 100.0e-9, source()).unwrap()
                )
                .content_identity(),
            );
        }
        for bad in [f64::NAN, f64::INFINITY, -1.0, 100.1] {
            assert_eq!(
                CoatingIorSample::try_new(500.0, 1.5, bad),
                Err(ConductorError::InvalidSample)
            );
        }
        for bad in [f64::NAN, -1.0, 2.0e-4] {
            assert_eq!(
                CoherentAbsorbingCoating::try_new(samples, bad, source()),
                Err(ConductorError::InvalidCoating),
            );
        }
        let mut unordered = samples;
        unordered.swap(2, 3);
        assert_eq!(
            CoherentAbsorbingCoating::try_new(unordered, 100.0e-9, source()),
            Err(ConductorError::InvalidTable),
        );
    }

    #[test]
    fn g0_passive_absorbing_film_stays_bounded_over_angle_and_optical_depth() {
        let bare = substrate(0.2, 3.0);
        for eta in [0.2, 1.0, 1.5, 3.0] {
            for k in [0.0, 1.0e-12, 0.1, 4.0, 100.0] {
                for thickness_m in [1.0e-14, 1.0e-9, 100.0e-9, 1.0e-4] {
                    let coated = bare.with_absorbing_coating(film(eta, k, thickness_m));
                    for wavelength in [380.0, 550.0, 780.0] {
                        for cosine in [0.0, 1.0e-8, 0.1, 0.6, 0.9, 1.0] {
                            let reflectance = coated.fresnel(wavelength, 1.5, cosine).unwrap();
                            assert!((0.0..=1.0).contains(&reflectance));
                        }
                    }
                }
            }
        }
    }
}