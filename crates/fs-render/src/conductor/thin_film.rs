//! Coherent single-layer optics for a transparent film on an opaque substrate.
//!
//! Maxwell boundary conditions give a 2x2 characteristic matrix for each
//! polarization. This is the closed-form sum of every internal reflection,
//! including the substrate's complex phase, rather than a sum of intensities.
//! Conventions follow S. J. Byrnes, *Multilayer optical calculations*, sections
//! 2--5 and appendix D: <https://arxiv.org/abs/1603.02720>.
//! Here fields use exp(-i omega t), passive substrates eta + i k, and the
//! forward longitudinal wave number has nonnegative imaginary part.

use super::{Complex, ComplexIor, ConductorDataStatus, ConductorError, ConductorSource};
use crate::dielectric::DielectricIor;
use fs_blake3::DomainHasher;
use fs_math::det;

/// Upper thickness of this explicitly coherent, uniform local coating model.
/// This 100 micrometre bound limits phase arguments and does not establish
/// that a particular source has sufficient temporal or spatial coherence.
pub const MAX_COHERENT_COATING_THICKNESS_M: f64 = 1.0e-4;

/// One transparent isotropic film, conformal to each conductor microfacet.
///
/// The source must describe its visible phase-index law and thickness.
/// The model sums coherent s/p amplitudes within this film and then returns
/// unpolarized reflected power. Thickness is normal to the local microfacet.
/// It assumes locally parallel boundaries and a laterally uniform film; the
/// caller must ensure coherence across its optical path difference. It adds
/// no absorption in the film, polarization tracking between bounces, arbitrary
/// multilayer stacks, thickness distributions, diffraction, or scene-wide
/// coherent propagation. The substrate retains its complex-index absorption.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoherentDielectricCoating {
    ior: DielectricIor,
    thickness_m: f64,
    source: ConductorSource,
}

impl CoherentDielectricCoating {
    /// Declare a source-bound film. Wavelength evaluation uses the existing
    /// shared visible-band Cauchy or Sellmeier law; thickness is in metres,
    /// including exact zero. The phase index controls both Snell's law and
    /// coherent optical thickness at each wavelength; group index is not used.
    pub fn try_new(
        ior: impl Into<DielectricIor>,
        thickness_m: f64,
        source: ConductorSource,
    ) -> Result<Self, ConductorError> {
        if !thickness_m.is_finite()
            || !(0.0..=MAX_COHERENT_COATING_THICKNESS_M).contains(&thickness_m)
        {
            return Err(ConductorError::InvalidCoating);
        }
        Ok(Self {
            ior: ior.into(),
            thickness_m: super::canonical_zero(thickness_m),
            source,
        })
    }

    /// Visible-band phase-index law, with no inferred material preset.
    #[must_use]
    pub const fn ior(self) -> DielectricIor {
        self.ior
    }

    /// Uniform physical thickness normal to the microfacet, in metres.
    #[must_use]
    pub const fn thickness_m(self) -> f64 {
        self.thickness_m
    }

    /// Caller-declared optical/thickness source, without upgraded authority.
    #[must_use]
    pub const fn source(self) -> ConductorSource {
        self.source
    }

    pub(super) fn update_identity(self, hasher: &mut DomainHasher) {
        match self.ior {
            DielectricIor::Cauchy(law) => {
                // Retain the original Cauchy coating's complete identity bytes.
                hasher.update(b"coherent-dielectric-coating.v1");
                for coefficient in law.coefficients() {
                    hasher.update(&coefficient.to_bits().to_le_bytes());
                }
            }
            DielectricIor::Sellmeier(law) => {
                hasher.update(b"coherent-dielectric-coating.sellmeier.v1");
                let (strengths, resonances_um2) = law.coefficients();
                for coefficient in strengths.into_iter().chain(resonances_um2) {
                    hasher.update(&coefficient.to_bits().to_le_bytes());
                }
            }
        }
        hasher.update(&self.thickness_m.to_bits().to_le_bytes());
        hasher.update(self.source.identity().as_bytes());
        hasher.update(&[match self.source.status() {
            ConductorDataStatus::Representative => 0,
            ConductorDataStatus::CallerAssertedMeasured => 1,
            ConductorDataStatus::MaterialStateClaims => 2,
        }]);
    }

    // Called after ConductorOptics has admitted wavelength and substrate IOR.
    pub(super) fn fresnel(
        self,
        wavelength_nm: f64,
        incident_eta: f64,
        substrate: ComplexIor,
        cosine: f64,
    ) -> Result<f64, ConductorError> {
        if !incident_eta.is_finite()
            || incident_eta <= 0.0
            || incident_eta > super::MAX_IOR_COMPONENT
        {
            return Err(ConductorError::InvalidIncidentMedium);
        }
        if !cosine.is_finite() || !(0.0..=1.0).contains(&cosine) {
            return Err(ConductorError::InvalidIncidence);
        }
        let film_eta = self
            .ior
            .eval(wavelength_nm)
            .map_err(|_| ConductorError::NumericalFailure)?;
        // With no film or no first interface, the intensity is exactly the
        // bare-substrate result. Keep its established bits in these limits.
        if self.thickness_m == 0.0 || film_eta == incident_eta || cosine == 0.0 {
            return super::fresnel_conductor(incident_eta, substrate.eta, substrate.k, cosine);
        }

        let n0_squared = incident_eta * incident_eta;
        let n1_squared = film_eta * film_eta;
        let transverse_squared = n0_squared * (1.0 - cosine * cosine);
        let q0 = incident_eta * cosine;
        let q1_squared = n1_squared - transverse_squared;
        let n2_squared = Complex::new(substrate.eta, substrate.k).square();
        let q2 = forward_sqrt(n2_squared - Complex::new(transverse_squared, 0.0));

        let phase_scale = core::f64::consts::TAU * (self.thickness_m * 1.0e9) / wavelength_nm;
        let (cos_phase, sin_phase_over_q) = propagation(q1_squared, phase_scale);

        // Using q/n^2 for p is the magnetic-field convention; at normal
        // incidence its reflected amplitude has the opposite sign from s.
        // Powers, and the coherent relative signs at BOTH interfaces, agree.
        let r_s = reflection(
            q0,
            q2,
            cos_phase,
            sin_phase_over_q,
            q1_squared * sin_phase_over_q,
        );
        let r_p = reflection(
            q0 / n0_squared,
            q2 / n2_squared,
            cos_phase,
            n1_squared * sin_phase_over_q,
            q1_squared * sin_phase_over_q / n1_squared,
        );
        let reflectance = 0.5 * (r_s.norm_squared() + r_p.norm_squared());
        if !reflectance.is_finite() || !(-1.0e-12..=1.0 + 1.0e-12).contains(&reflectance) {
            return Err(ConductorError::NumericalFailure);
        }
        Ok(reflectance.clamp(0.0, 1.0))
    }
}

// The passive substrate puts n^2-k_parallel^2 in the upper half-plane.
// Recover the small component from 2*re*im=value.im instead of subtracting
// nearly equal magnitudes; weak absorption must survive interference nulls.
fn forward_sqrt(value: Complex) -> Complex {
    let magnitude = det::hypot(value.re, value.im);
    if magnitude == 0.0 {
        return Complex::new(0.0, 0.0);
    }
    let large = det::sqrt(magnitude) * det::sqrt(0.5 * (1.0 + value.re.abs() / magnitude));
    let small = value.im / (2.0 * large);
    if value.re >= 0.0 {
        Complex::new(large, small)
    } else {
        Complex::new(small, large)
    }
}

// Matrix entries cos(a q) and sin(a q)/q, up to a common real scale that
// cancels from reflection. The sinc limit removes the Airy formula's 0/0
// at the film's critical angle. Evanescent layers use exp(-2 a |q|), never
// exponentially growing transfer entries, so frustrated TIR stays finite.
fn propagation(q_squared: f64, phase_scale: f64) -> (f64, f64) {
    if q_squared >= 0.0 {
        let q = det::sqrt(q_squared);
        let phase = phase_scale * q;
        let (sine, cosine) = det::sin_cos(phase);
        let sine_over_q = if phase.abs() < 1.0e-4 {
            let phase_squared = phase * phase;
            phase_scale * (1.0 - phase_squared / 6.0 + phase_squared * phase_squared / 120.0)
        } else {
            sine / q
        };
        (cosine, sine_over_q)
    } else {
        let kappa = det::sqrt(-q_squared);
        let decay = 2.0 * phase_scale * kappa;
        (
            0.5 * (1.0 + det::exp(-decay)),
            -0.5 * det::expm1(-decay) / kappa,
        )
    }
}

fn reflection(
    incident_admittance: f64,
    substrate_admittance: Complex,
    cosine: f64,
    sine_over_admittance: f64,
    admittance_times_sine: f64,
) -> Complex {
    // M = [ cos(delta), -i sin(delta)/y; -i y sin(delta), cos(delta) ].
    let field = Complex::new(cosine, 0.0)
        - Complex::new(-substrate_admittance.im, substrate_admittance.re) * sine_over_admittance;
    let flux = substrate_admittance * cosine - Complex::new(0.0, admittance_times_sine);
    (field * incident_admittance - flux) / (field * incident_admittance + flux)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conductor::{CONDUCTOR_IOR_SAMPLE_COUNT, ConductorIorSample, ConductorOptics};
    use crate::dielectric::{CauchyIor, SellmeierIor};

    fn source() -> ConductorSource {
        ConductorSource::try_new(
            fs_blake3::hash_domain("org.frankensim.test.coherent-film", b"synthetic"),
            ConductorDataStatus::Representative,
        )
        .unwrap()
    }

    fn coating(index: f64, thickness_m: f64) -> CoherentDielectricCoating {
        CoherentDielectricCoating::try_new(
            CauchyIor::try_constant(index).unwrap(),
            thickness_m,
            source(),
        )
        .unwrap()
    }

    fn substrate(eta: f64, k: f64) -> ConductorOptics {
        let samples = core::array::from_fn::<_, CONDUCTOR_IOR_SAMPLE_COUNT, _>(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, eta, k).unwrap()
        });
        ConductorOptics::try_new(samples, source()).unwrap()
    }

    fn close(actual: f64, expected: f64, tolerance: f64) {
        assert!(
            (actual - expected).abs() <= tolerance,
            "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
        );
    }

    #[test]
    fn g1_coherent_film_matches_independent_complex_airy_answers() {
        // Independent oracle: r=(r01+r12 exp(2i delta))/(1+r01 r12 exp(2i delta)),
        // evaluated separately for s/p with Python's cmath. Production uses
        // boundary-field characteristic matrices, including the critical limit.
        for (n0, n1, eta, k, d_nm, lambda, cosine, expected) in [
            (
                1.0,
                1.45,
                0.2,
                3.0,
                125.0,
                550.0,
                1.0,
                0.916_288_333_911_421_1,
            ),
            (
                1.0,
                1.45,
                0.2,
                3.0,
                125.0,
                550.0,
                0.5,
                0.892_500_332_685_915_1,
            ),
            (
                1.33,
                2.1,
                2.7,
                4.1,
                82.0,
                430.0,
                0.73,
                0.536_716_136_691_466_2,
            ),
            (
                1.5,
                1.0,
                0.2,
                3.0,
                200.0,
                550.0,
                0.5,
                0.992_093_574_833_073_9,
            ),
            (1.0, 2.0, 3.0, 4.0, 137.0, 610.0, 0.2, 0.551_252_088_805_326),
        ] {
            let optics = substrate(eta, k).with_coating(coating(n1, d_nm * 1.0e-9));
            close(
                optics.fresnel(lambda, n0, cosine).unwrap(),
                expected,
                4.0e-14,
            );
        }
    }

    #[test]
    fn g1_coherent_film_quarter_wave_cancels_and_half_wave_recovers_bare() {
        let (n0, n1, eta, k) = (1.0_f64, 1.5_f64, 2.25_f64, 1.0e-8_f64);
        let lambda = 600.0;
        let bare = substrate(eta, k);
        let quarter = bare.with_coating(coating(n1, lambda * 1.0e-9 / (4.0 * n1)));
        let half = bare.with_coating(coating(n1, lambda * 1.0e-9 / (2.0 * n1)));
        // The quarter-wave input admittance is n1^2/n_substrate. For
        // n1^2=n0*eta, its real reflection cancels; weak substrate loss remains.
        let expected = ((n0 * eta - n1 * n1).powi(2) + (n0 * k).powi(2))
            / ((n0 * eta + n1 * n1).powi(2) + (n0 * k).powi(2));
        let quarter_reflection = quarter.fresnel(lambda, n0, 1.0).unwrap();
        close(quarter_reflection, expected, 1.0e-23);
        assert!(quarter_reflection < 1.0e-16);
        assert!(bare.fresnel(lambda, n0, 1.0).unwrap() > 0.14);
        close(
            half.fresnel(lambda, n0, 1.0).unwrap(),
            bare.fresnel(lambda, n0, 1.0).unwrap(),
            2.0e-14,
        );
        assert!(quarter.fresnel(430.0, n0, 1.0).unwrap() > 0.04);
    }

    #[test]
    fn g3_coherent_film_zero_and_index_matched_limits_keep_bare_bits() {
        let bare = ConductorOptics::representative_tungsten();
        for n0 in [1.0, 1.33, 1.5] {
            for film in [coating(2.1, 0.0), coating(n0, 123.0e-9)] {
                let coated = bare.with_coating(film);
                for lambda in [380.0, 510.0, 780.0] {
                    for cosine in [0.0, 0.01, 0.6, 1.0] {
                        assert_eq!(
                            coated.fresnel(lambda, n0, cosine).unwrap().to_bits(),
                            bare.fresnel(lambda, n0, cosine).unwrap().to_bits(),
                        );
                    }
                }
            }
        }
        let scaled = substrate(0.2, 3.0);
        close(
            scaled
                .with_coating(coating(1.45, 75.0e-9))
                .fresnel(500.0, 1.0, 0.7)
                .unwrap(),
            scaled
                .with_coating(coating(1.45, 112.5e-9))
                .fresnel(750.0, 1.0, 0.7)
                .unwrap(),
            2.0e-14,
        );
    }

    #[test]
    fn g1_coherent_film_critical_and_evanescent_limits_are_continuous() {
        let bare = substrate(0.2, 3.0);
        let thin = bare.with_coating(coating(1.0, 200.0e-9));
        let critical = (1.0_f64 - (1.0_f64 / 1.5).powi(2)).sqrt();
        let at_critical = thin.fresnel(550.0, 1.5, critical).unwrap();
        assert!((0.0..1.0).contains(&at_critical));
        for offset in [-1.0e-9, -1.0e-12, 1.0e-12, 1.0e-9] {
            close(
                thin.fresnel(550.0, 1.5, critical + offset).unwrap(),
                at_critical,
                3.0e-9,
            );
        }
        // Below the critical cosine an evanescent film still couples energy
        // into the absorbing substrate. Returning TIR=1 would erase tunnelling.
        assert!(thin.fresnel(550.0, 1.5, 0.5).unwrap() < 0.999);
        let thick = bare.with_coating(coating(1.0, MAX_COHERENT_COATING_THICKNESS_M));
        close(thick.fresnel(550.0, 1.5, 0.5).unwrap(), 1.0, 2.0e-14);
    }

    #[test]
    fn g1_coherent_film_critical_substrate_retains_tiny_positive_loss() {
        let zero = forward_sqrt(Complex::new(0.0, 0.0));
        assert_eq!([zero.re, zero.im], [0.0, 0.0]);
        let tiny = forward_sqrt(Complex::new(0.0, 4.0e-300));
        let expected = (2.0e-300_f64).sqrt();
        assert!((tiny.re / expected - 1.0).abs() < 1.0e-14);
        assert!((tiny.im / expected - 1.0).abs() < 1.0e-14);
        // n0 sin(theta)=2 exactly; substrate k>0 is admitted even where
        // its square and an unscaled norm would underflow. The lossless
        // critical-substrate limit is unit reflection through this film.
        let optics = substrate(2.0, 1.0e-300).with_coating(coating(1.5, 20.0e-9));
        close(optics.fresnel(550.0, 2.5, 0.6).unwrap(), 1.0, 2.0e-14);
    }

    #[test]
    fn g0_coherent_film_refuses_bad_thickness_and_binds_physical_inputs() {
        for thickness in [
            -1.0e-9,
            f64::NAN,
            f64::INFINITY,
            2.0 * MAX_COHERENT_COATING_THICKNESS_M,
        ] {
            assert_eq!(
                CoherentDielectricCoating::try_new(
                    CauchyIor::try_constant(1.5).unwrap(),
                    thickness,
                    source()
                ),
                Err(ConductorError::InvalidCoating),
            );
        }
        let bare = ConductorOptics::representative_stainless_steel();
        let film = coating(1.45, 100.0e-9);
        let coated = bare.with_coating(film);
        assert_eq!(coated.coating(), Some(film));
        assert_eq!(film.thickness_m(), 100.0e-9);
        assert_eq!(
            film.ior(),
            DielectricIor::Cauchy(CauchyIor::try_constant(1.45).unwrap()),
        );
        assert_eq!(film.source(), source());
        assert_ne!(coated.content_identity(), bare.content_identity());
        for other in [coating(1.45, 101.0e-9), coating(1.46, 100.0e-9)] {
            assert_ne!(
                coated.content_identity(),
                bare.with_coating(other).content_identity()
            );
        }
        let dispersive = CoherentDielectricCoating::try_new(
            CauchyIor::try_new(1.45, 0.01, 0.0).unwrap(),
            100.0e-9,
            source(),
        )
        .unwrap();
        assert!(
            (coated.fresnel(400.0, 1.0, 0.65).unwrap()
                - bare
                    .with_coating(dispersive)
                    .fresnel(400.0, 1.0, 0.65)
                    .unwrap())
            .abs()
                > 1.0e-3
        );
        assert_eq!(
            coated.fresnel(379.0, 1.0, 1.0),
            Err(ConductorError::InvalidWavelength)
        );
        assert_eq!(
            coated.fresnel(500.0, 0.0, 1.0),
            Err(ConductorError::InvalidIncidentMedium)
        );
        assert_eq!(
            coated.fresnel(500.0, 1.0, f64::NAN),
            Err(ConductorError::InvalidIncidence)
        );
    }

    #[test]
    fn g5_cauchy_coating_identity_stays_frozen() {
        let source = ConductorSource::try_new(
            fs_blake3::ContentHash([77; 32]),
            ConductorDataStatus::Representative,
        )
        .unwrap();
        let samples = core::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, 0.2, 3.0).unwrap()
        });
        let optics = ConductorOptics::try_new(samples, source)
            .unwrap()
            .with_coating(
                CoherentDielectricCoating::try_new(
                    CauchyIor::try_new(1.45, 0.01, 0.0).unwrap(),
                    100.0e-9,
                    source,
                )
                .unwrap(),
            );
        // Captured by executing the pre-Sellmeier production coating/hash owner.
        assert_eq!(
            *optics.content_identity().as_bytes(),
            [
                249, 199, 58, 11, 211, 199, 239, 151, 74, 39, 38, 39, 85, 123, 191, 224, 131, 116,
                98, 82, 67, 86, 162, 56, 163, 172, 212, 1, 234, 59, 237, 23
            ],
        );
    }

    fn sellmeier_coating(thickness_m: f64) -> CoherentDielectricCoating {
        CoherentDielectricCoating::try_new(
            SellmeierIor::try_new([1.1, 0.2, 0.8], [0.01, 0.04, 100.0]).unwrap(),
            thickness_m,
            source(),
        )
        .unwrap()
    }

    #[test]
    fn g1_sellmeier_film_matches_independent_complex_airy_spectrum() {
        // Synthetic Sellmeier strengths/resonances, evaluated independently in
        // Python with the complex Airy formula from Byrnes, sections 2--5.
        // These answers include both propagating and evanescent films.
        for (n0, eta, k, thickness_nm, wavelength_nm, cosine, expected) in [
            (1.0, 0.2, 3.0, 125.0, 400.0, 1.0, 0.9261497100473655),
            (1.0, 0.2, 3.0, 125.0, 550.0, 0.5, 0.9024456746212396),
            (1.33, 2.7, 4.1, 82.0, 430.0, 0.73, 0.46954541554138274),
            (2.4, 0.2, 3.0, 200.0, 550.0, 0.5, 0.9975089344662884),
            (1.0, 3.0, 4.0, 137.0, 750.0, 0.2, 0.36497674703594385),
            (1.0, 0.2, 3.0, 125.0, 780.0, 1.0, 0.8756898552426906),
            (1.33, 2.7, 4.1, 82.0, 650.0, 0.73, 0.49067376008162433),
        ] {
            let optics = substrate(eta, k).with_coating(sellmeier_coating(thickness_nm * 1.0e-9));
            close(
                optics.fresnel(wavelength_nm, n0, cosine).unwrap(),
                expected,
                4.0e-14,
            );
        }
    }

    #[test]
    fn g3_sellmeier_film_constant_and_zero_thickness_limits_keep_fresnel_bits() {
        let bare = substrate(0.2, 3.0);
        let constant = SellmeierIor::try_new([1.25, 0.0, 0.0], [0.0; 3]).unwrap();
        let sellmeier = bare.with_coating(
            CoherentDielectricCoating::try_new(constant, 125.0e-9, source()).unwrap(),
        );
        let cauchy = bare.with_coating(coating(1.5, 125.0e-9));
        assert_ne!(sellmeier.content_identity(), cauchy.content_identity());
        for wavelength in [380.0, 430.0, 550.0, 710.0, 780.0] {
            for incident in [1.0, 1.33, 2.4] {
                for cosine in [0.0, 1.0e-6, 0.05, 0.5, 1.0] {
                    assert_eq!(
                        sellmeier
                            .fresnel(wavelength, incident, cosine)
                            .unwrap()
                            .to_bits(),
                        cauchy
                            .fresnel(wavelength, incident, cosine)
                            .unwrap()
                            .to_bits(),
                    );
                    assert_eq!(
                        bare.with_coating(sellmeier_coating(0.0))
                            .fresnel(wavelength, incident, cosine)
                            .unwrap()
                            .to_bits(),
                        bare.fresnel(wavelength, incident, cosine)
                            .unwrap()
                            .to_bits(),
                    );
                }
            }
        }
    }

    #[test]
    fn g3_sellmeier_film_dispersion_and_every_coefficient_reach_transport() {
        let bare = substrate(0.2, 3.0);
        let film = sellmeier_coating(125.0e-9);
        let optics = bare.with_coating(film);
        let reference =
            bare.with_coating(coating(film.ior().eval(550.0).unwrap(), film.thickness_m()));
        close(
            optics.fresnel(550.0, 1.0, 0.6).unwrap(),
            reference.fresnel(550.0, 1.0, 0.6).unwrap(),
            1.0e-14,
        );
        assert!(
            (optics.fresnel(400.0, 1.0, 0.6).unwrap()
                - reference.fresnel(400.0, 1.0, 0.6).unwrap())
            .abs()
                > 1.0e-4
        );
        let DielectricIor::Sellmeier(law) = film.ior() else {
            unreachable!()
        };
        for index in 0..6 {
            let (mut strengths, mut resonances) = law.coefficients();
            if index < 3 {
                strengths[index] *= 1.01;
            } else {
                resonances[index - 3] *= 1.01;
            }
            let changed = bare.with_coating(
                CoherentDielectricCoating::try_new(
                    SellmeierIor::try_new(strengths, resonances).unwrap(),
                    film.thickness_m(),
                    source(),
                )
                .unwrap(),
            );
            assert_ne!(changed.content_identity(), optics.content_identity());
            assert_ne!(
                changed.fresnel(550.0, 1.0, 0.6).unwrap().to_bits(),
                optics.fresnel(550.0, 1.0, 0.6).unwrap().to_bits(),
            );
        }
        for thickness in [1.0e-12, 100.0e-9, MAX_COHERENT_COATING_THICKNESS_M] {
            let coated = bare.with_coating(sellmeier_coating(thickness));
            for wavelength in [380.0, 430.0, 550.0, 710.0, 780.0] {
                for incident in [1.0, 1.33, 2.4] {
                    for cosine in [0.0, 1.0e-6, 0.05, 0.2, 0.5, 0.8, 1.0] {
                        let response = coated.fresnel(wavelength, incident, cosine).unwrap();
                        assert!((0.0..=1.0).contains(&response));
                    }
                }
            }
        }
    }

    #[test]
    fn g0_coherent_film_visible_angular_grid_is_passive() {
        for bare in [
            ConductorOptics::representative_tungsten(),
            ConductorOptics::representative_stainless_steel(),
        ] {
            for n0 in [1.0, 1.33, 2.4] {
                for n1 in [1.0, 1.45, 2.0, 3.0] {
                    for thickness in [1.0e-12, 100.0e-9, MAX_COHERENT_COATING_THICKNESS_M] {
                        let optics = bare.with_coating(coating(n1, thickness));
                        for lambda in [380.0, 430.0, 550.0, 710.0, 780.0] {
                            for cosine in [0.0, 1.0e-6, 0.05, 0.2, 0.5, 0.8, 0.999, 1.0] {
                                let r = optics.fresnel(lambda, n0, cosine).unwrap();
                                assert!((0.0..=1.0).contains(&r));
                            }
                        }
                    }
                }
            }
        }
    }
}
