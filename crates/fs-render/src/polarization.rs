//! Allocation-free Stokes/Mueller polarization algebra for spectral light transport.
//!
//! Stokes components use the convention S=(I,Q,U,V), with Q positive for
//! horizontal linear polarization, U positive at +45 degrees and V positive
//! for the chosen right-handed circular basis. All operations assume the same
//! propagation direction and transverse reference frame. Rotate frames before
//! composing optical elements; this module does not silently rotate them.
//! Incoherent light is represented by Stokes vectors, not Jones amplitudes.

/// Polarized Fresnel, spectral materials and native smooth-sampler integration.
#[cfg(feature = "tracer")]
pub mod dielectric;

/// A physical (possibly partially polarized) incoherent light state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Stokes {
    pub i: f64,
    pub q: f64,
    pub u: f64,
    pub v: f64,
}

impl Stokes {
    /// Reject nonfinite, negative-intensity or super-physical polarization.
    pub fn try_new(i: f64, q: f64, u: f64, v: f64) -> Option<Self> {
        if ![i, q, u, v].iter().all(|x| x.is_finite()) || i < 0.0 {
            return None;
        }
        // Division by I avoids squaring large radiances and overflow.
        if i == 0.0 {
            return (q == 0.0 && u == 0.0 && v == 0.0).then_some(Self { i, q, u, v });
        }
        let norm = (q / i).hypot(u / i).hypot(v / i);
        (norm <= 1.0 + 16.0 * f64::EPSILON).then_some(Self { i, q, u, v })
    }

    pub const fn unpolarized(i: f64) -> Self {
        Self { i, q: 0.0, u: 0.0, v: 0.0 }
    }

    /// Rotate the reference axes counterclockwise by angle (radians).
    /// The Stokes Q/U coordinates rotate by minus twice the angle.
    pub fn rotate_axes(self, angle: f64) -> Self {
        let (s, c) = (2.0 * angle).sin_cos();
        Self { i: self.i, q: c * self.q + s * self.u,
            u: -s * self.q + c * self.u, v: self.v }
    }

    /// Degree of polarization; zero for a dark state.
    pub fn degree(self) -> f64 {
        if self.i == 0.0 { 0.0 } else {
            (self.q / self.i).hypot(self.u / self.i).hypot(self.v / self.i)
        }
    }
}

/// Mueller matrix acting on a Stokes column vector. The matrix is row-major.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Mueller(pub [[f64; 4]; 4]);

impl Mueller {
    pub const IDENTITY: Self = Self([
        [1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0],
    ]);

    /// Passive ideal linear polarizer with intensity throughput 1/2
    /// for unpolarized light, oriented at angle from the Q-positive axis.
    pub fn linear_polarizer(angle: f64) -> Self {
        let (s, c) = (2.0 * angle).sin_cos();
        let a = [1.0, c, s, 0.0];
        Self(core::array::from_fn(|r| core::array::from_fn(|k| 0.5 * a[r] * a[k])))
    }

    /// Ideal lossless linear retarder with phase difference in radians.
    /// At zero orientation it leaves Q unchanged and rotates the U/V plane.
    pub fn linear_retarder(angle: f64, phase: f64) -> Self {
        let (s, c) = (2.0 * angle).sin_cos();
        let (sp, cp) = phase.sin_cos();
        let n = [c, s];
        let a = [[cp + (1.0-cp)*n[0]*n[0], (1.0-cp)*n[0]*n[1]],
                 [(1.0-cp)*n[1]*n[0], cp + (1.0-cp)*n[1]*n[1]]];
        Self([
            [1.0, 0.0, 0.0, 0.0],
            [0.0, a[0][0], a[0][1], -s*sp],
            [0.0, a[1][0], a[1][1], c*sp],
            [0.0, s*sp, -c*sp, cp],
        ])
    }

    pub fn apply(self, input: Stokes) -> Stokes {
        let x = [input.i, input.q, input.u, input.v];
        let y = self.0.map(|row| row.iter().zip(x).map(|(a,b)| a*b).sum::<f64>());
        Stokes { i: y[0], q: y[1], u: y[2], v: y[3] }
    }

    /// Composition: self acts after rhs.
    pub fn compose(self, rhs: Self) -> Self {
        Self(core::array::from_fn(|r| core::array::from_fn(|c|
            (0..4).map(|k| self.0[r][k] * rhs.0[k][c]).sum()
        )))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn near(a: f64, b: f64) { assert!((a-b).abs() < 1e-12, "{a} != {b}"); }

    #[test]
    fn malus_law_and_crossed_polarizers() {
        let source = Stokes::unpolarized(10.0);
        let horizontal = Mueller::linear_polarizer(0.0).apply(source);
        near(horizontal.i, 5.0);
        near(horizontal.degree(), 1.0);
        near(Mueller::linear_polarizer(core::f64::consts::FRAC_PI_2).apply(horizontal).i, 0.0);
        near(Mueller::linear_polarizer(core::f64::consts::FRAC_PI_4).apply(horizontal).i, 2.5);
    }

    #[test]
    fn retarder_preserves_intensity_and_polarization_degree() {
        let source = Stokes::try_new(2.0, 1.0, 1.0, 0.0).unwrap();
        let out = Mueller::linear_retarder(0.0, core::f64::consts::FRAC_PI_2).apply(source);
        near(out.i, source.i);
        near(out.q, 1.0);
        near(out.u, 0.0);
        near(out.v, -1.0);
        near(out.degree(), source.degree());
    }

    #[test]
    fn rotations_and_composition() {
        let x = Stokes::try_new(1.0, 1.0, 0.0, 0.0).unwrap();
        let rotated = x.rotate_axes(core::f64::consts::FRAC_PI_4);
        near(rotated.q, 0.0);
        near(rotated.u, -1.0);
        let a = Mueller::linear_polarizer(0.3);
        let b = Mueller::linear_retarder(0.7, 0.8);
        let lhs = a.compose(b).apply(x);
        let rhs = a.apply(b.apply(x));
        near(lhs.i, rhs.i); near(lhs.q, rhs.q); near(lhs.u, rhs.u); near(lhs.v, rhs.v);
    }

    #[test]
    fn invalid_stokes_states_are_rejected() {
        assert!(Stokes::try_new(1.0, 1.1, 0.0, 0.0).is_none());
        assert!(Stokes::try_new(-1.0, 0.0, 0.0, 0.0).is_none());
        assert!(Stokes::try_new(0.0, 0.0, 0.0, 1.0).is_none());
        assert!(Stokes::try_new(f64::NAN, 0.0, 0.0, 0.0).is_none());
    }
}

/// Polarization-preserving optics connected to calibrated fixture output.
#[cfg(feature = "photometry")]
pub mod photometry {
    //! Polarization optics connected to calibrated fixture spectra and lumens.
    //!
    //! A chain retains the full Stokes state in each existing 380--780 nm bin.
    //! Only its final intensity is projected into a native `PhotometricSource`,
    //! so crossed polarizers extinguish light instead of independently halving it.
    //! The original source spectrum/flux remains the normalization reference.
    //!
    //! This is a uniform, paraxial source-optics model in one transverse frame,
    //! not a spatial ray-optics solve. Optics do not change placement or angular
    //! distribution. The caller places each outgoing Fresnel channel. Birefringent
    //! waveplates model phase delay only, without walk-off or interface reflections.
    //! The native intensity source preserves spectral power, but does NOT make the
    //! existing scalar path tracer propagate Q/U/V after emission.

    use core::f64::consts::{PI, TAU};

    use crate::dielectric::{DielectricEvent, DielectricIor};
    use crate::photometry::source::{PhotometricSource, SPECTRAL_BINS, SourceError, wavelength_nm};
    use crate::polarization::dielectric::spectral_fresnel_dielectric_polarized_borrowed;
    use crate::polarization::{Mueller, Stokes};

    /// Precomputed passive optics at the existing photometric quadrature bins.
    /// Construction is restricted to physical optical elements, not arbitrary
    /// caller-authored Mueller matrices. No per-ray matrix work is added.
    #[derive(Clone, Debug)]
    pub struct SpectralOptics {
        matrices: [Mueller; SPECTRAL_BINS],
    }

    impl SpectralOptics {
        /// Linear polarizer with distinct power transmissions along its principal
        /// axis and perpendicular axis. Both arrays must lie in `[0,1]`.
        /// This includes finite extinction ratio and wavelength-dependent loss.
        pub fn linear_diattenuator(
            angle_rad: f64,
            parallel: [f64; SPECTRAL_BINS],
            perpendicular: [f64; SPECTRAL_BINS],
        ) -> Result<Self, SourceError> {
            if !angle_rad.is_finite() || parallel.iter().chain(&perpendicular)
                .any(|t| !t.is_finite() || !(0.0..=1.0).contains(t)) {
                return Err(SourceError::InvalidFilter);
            }
            let (s, c) = fs_math::det::sin_cos(2.0 * angle_rad.rem_euclid(PI));
            let matrices = core::array::from_fn(|bin| {
                let (p, t) = (parallel[bin], perpendicular[bin]);
                if p == t {
                    return Mueller(Mueller::IDENTITY.0.map(|row| row.map(|value| value * p)));
                }
                let a = 0.5 * p + 0.5 * t;
                let b = 0.5 * p - 0.5 * t;
                let g = p.sqrt() * t.sqrt();
                Mueller([
                    [a, b * c, b * s, 0.0],
                    [b * c, a * c * c + g * s * s, (a - g) * c * s, 0.0],
                    [b * s, (a - g) * c * s, a * s * s + g * c * c, 0.0],
                    [0.0, 0.0, 0.0, g],
                ])
            });
            Ok(Self { matrices })
        }

        /// Ideal achromatic polarizer. Unpolarized light loses half its power.
        pub fn linear_polarizer(angle_rad: f64) -> Result<Self, SourceError> {
            Self::linear_diattenuator(angle_rad, [1.0; SPECTRAL_BINS], [0.0; SPECTRAL_BINS])
        }

        /// A neutral spectral filter attenuates all polarization components equally.
        pub fn neutral(transmission: [f64; SPECTRAL_BINS]) -> Result<Self, SourceError> {
            Self::linear_diattenuator(0.0, transmission, transmission)
        }

        /// Lossless linear retarder. Positive phase delays the perpendicular
        /// component relative to the principal axis, matching `Mueller` conventions.
        pub fn linear_retarder(
            angle_rad: f64,
            phase_rad: [f64; SPECTRAL_BINS],
        ) -> Result<Self, SourceError> {
            if !angle_rad.is_finite() || phase_rad.iter().any(|phase| !phase.is_finite()) {
                return Err(SourceError::InvalidFilter);
            }
            // Reduce before doubling to avoid overflow for finite authored angles.
            let angle = angle_rad.rem_euclid(PI);
            let matrices = phase_rad.map(|phase| Mueller::linear_retarder(angle, phase.rem_euclid(TAU)));
            Ok(Self { matrices })
        }

        /// Phase-only, normal-incidence waveplate using two native spectral index
        /// laws. `thickness_m` is metres; the indices are the two principal modes.
        /// Interface losses and walk-off are deliberately not inferred.
        pub fn waveplate(
            angle_rad: f64,
            principal: impl Into<DielectricIor>,
            perpendicular: impl Into<DielectricIor>,
            thickness_m: f64,
        ) -> Result<Self, SourceError> {
            if !thickness_m.is_finite() || thickness_m < 0.0 {
                return Err(SourceError::InvalidFilter);
            }
            let principal = principal.into();
            let perpendicular = perpendicular.into();
            let mut phase = [0.0; SPECTRAL_BINS];
            for (bin, value) in phase.iter_mut().enumerate() {
                let wavelength = wavelength_nm(bin);
                let a = principal
                    .eval(wavelength)
                    .map_err(|_| SourceError::InvalidFilter)?;
                let b = perpendicular
                    .eval(wavelength)
                    .map_err(|_| SourceError::InvalidFilter)?;
                *value = TAU * ((b - a) * thickness_m / (wavelength * 1e-9));
                if !value.is_finite() {
                    return Err(SourceError::NumericalRange);
                }
            }
            Self::linear_retarder(angle_rad, phase)
        }

        /// One Fresnel power channel at a fixed incidence angle, in s/p frames.
        /// Reflection and transmission can be applied to separate copies of the
        /// original source; the caller owns the corresponding outgoing placement.
        /// Both indices use the native Cauchy, Sellmeier or tabulated phase law.
        pub fn dielectric_channel(
            incident_cosine: f64,
            incident: impl Into<DielectricIor>,
            transmitted: impl Into<DielectricIor>,
            event: DielectricEvent,
        ) -> Result<Self, SourceError> {
            let incident = incident.into();
            let transmitted = transmitted.into();
            let mut matrices = [Mueller::IDENTITY; SPECTRAL_BINS];
            for (bin, matrix) in matrices.iter_mut().enumerate() {
                let f = spectral_fresnel_dielectric_polarized_borrowed(
                    incident_cosine,
                    &incident,
                    &transmitted,
                    wavelength_nm(bin),
                )
                .map_err(|_| SourceError::InvalidFilter)?;
                *matrix = match event {
                    DielectricEvent::Reflection => f.reflection,
                    DielectricEvent::Transmission => f.transmission,
                };
            }
            Ok(Self { matrices })
        }
    }

    /// A calibrated source plus its retained wavelength-dependent polarization.
    #[derive(Clone, Debug)]
    pub struct PolarizedSource {
        reference: PhotometricSource,
        intensity: PhotometricSource,
        relative: [Stokes; SPECTRAL_BINS],
    }

    impl PolarizedSource {
        /// Wrap an existing calibrated source without changing its emitted flux.
        #[must_use]
        pub fn unpolarized(source: PhotometricSource) -> Self {
            Self { intensity: source.clone(), reference: source,
                relative: [Stokes::unpolarized(1.0); SPECTRAL_BINS] }
        }

        /// Set normalized per-bin Stokes states. I must be exactly one so authored
        /// polarization cannot silently replace the source's calibrated brightness.
        pub fn try_new(
            source: PhotometricSource,
            normalized: [Stokes; SPECTRAL_BINS],
        ) -> Result<Self, SourceError> {
            if normalized.iter().any(|s| s.i != 1.0
                || Stokes::try_new(s.i, s.q, s.u, s.v).is_none()) {
                return Err(SourceError::InvalidSpectrum);
            }
            Ok(Self { intensity: source.clone(), reference: source, relative: normalized })
        }

        /// Apply the next optic without discarding Q/U/V or restoring lost flux.
        /// Uses bounded fixed-size work; failures leave the original unchanged.
        pub fn through(&self, optic: &SpectralOptics) -> Result<Self, SourceError> {
            let mut relative = self.relative;
            for (bin, state) in relative.iter_mut().enumerate() {
                *state = passive_roundoff(optic.matrices[bin].apply(*state), state.i)?;
            }
            // Project absolute surviving power only once relative to the ORIGINAL
            // spectrum. Normalizing Stokes after each element would erase extinction.
            let intensity = self.reference.filtered(relative.map(|s| s.i))?;
            Ok(Self { reference: self.reference.clone(), intensity, relative })
        }

        /// Per-bin Stokes power relative to the original source's spectral power.
        /// At dark bins the state is exactly zero, not an arbitrary normalized axis.
        #[must_use]
        pub const fn relative_stokes(&self) -> &[Stokes; SPECTRAL_BINS] { &self.relative }

        /// Native calibrated intensity view, usable by existing illuminance,
        /// haze and spectral rendering consumers. Q/U/V is retained in this wrapper,
        /// not propagated by those scalar consumers.
        #[must_use]
        pub const fn intensity_source(&self) -> &PhotometricSource { &self.intensity }

        /// Explicitly finish source optics and hand calibrated intensity to a scene.
        #[must_use]
        pub fn into_intensity_source(self) -> PhotometricSource { self.intensity }
    }

    impl PhotometricSource {
        /// Start a polarization-preserving optical chain on an unpolarized fixture.
        #[must_use]
        pub fn with_unpolarized_optics(self) -> PolarizedSource { PolarizedSource::unpolarized(self) }
    }

    // Only private, physically constructed matrices reach this boundary. Project
    // floating-point roundoff at zero intensity/the Stokes cone; never admit a
    // general unphysical matrix or a meaningful power gain by clamping it away.
    fn passive_roundoff(mut s: Stokes, incident: f64) -> Result<Stokes, SourceError> {
        let tolerance = 64.0 * f64::EPSILON * incident;
        let polarized = s.q.hypot(s.u).hypot(s.v);
        if ![s.i, s.q, s.u, s.v].iter().all(|v| v.is_finite())
            || s.i < -tolerance || s.i > incident + tolerance
            || polarized > s.i.max(0.0) + tolerance {
            return Err(SourceError::NumericalRange);
        }
        s.i = s.i.clamp(0.0, incident);
        if s.i == 0.0 { return Ok(Stokes::unpolarized(0.0)); }
        if polarized > s.i {
            let scale = s.i / polarized;
            s.q *= scale;
            s.u *= scale;
            s.v *= scale;
        }
        Stokes::try_new(s.i, s.q, s.u, s.v).ok_or(SourceError::NumericalRange)
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use crate::dielectric::{CauchyIor, DielectricIorSample, SellmeierIor, TabulatedIor};
        use crate::instances::RigidTransform;
        use crate::photometry::source::{AngularDistribution, SourceSpectrum};
        use fs_geom::Vec3;

        fn source() -> PhotometricSource {
            PhotometricSource::try_new(
                RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0; 3]).unwrap(),
                1000.0, 1.0, AngularDistribution::isotropic(),
                SourceSpectrum::try_from_bins([1.0; SPECTRAL_BINS]).unwrap(),
            ).unwrap()
        }

        fn close(a: f64, b: f64) { assert!((a - b).abs() < 2e-10, "{a} != {b}"); }

        #[test]
        fn crossed_and_rotated_polarizers_change_native_flux_and_candela() {
            let original = source();
            let horizontal = original.clone().with_unpolarized_optics()
                .through(&SpectralOptics::linear_polarizer(0.0).unwrap()).unwrap();
            close(horizontal.intensity_source().flux_lm(), 500.0);
            close(horizontal.intensity_source().intensity_cd(Vec3::new(0.0, 0.0, 1.0)).unwrap(),
                original.intensity_cd(Vec3::new(0.0, 0.0, 1.0)).unwrap() * 0.5);
            let crossed = horizontal.through(&SpectralOptics::linear_polarizer(PI / 2.0).unwrap()).unwrap();
            close(crossed.intensity_source().flux_lm(), 0.0);
            let angled = horizontal.through(&SpectralOptics::linear_polarizer(PI / 4.0).unwrap()).unwrap();
            close(angled.into_intensity_source().flux_lm(), 250.0);
            close(original.flux_lm(), 1000.0);
        }

        #[test]
        fn retarder_retains_circular_state_for_the_next_analyzer() {
            let linear = [Stokes::try_new(1.0, 0.0, 1.0, 0.0).unwrap(); SPECTRAL_BINS];
            let quarter = SpectralOptics::linear_retarder(0.0, [PI / 2.0; SPECTRAL_BINS]).unwrap();
            let circular = PolarizedSource::try_new(source(), linear).unwrap().through(&quarter).unwrap();
            close(circular.intensity_source().flux_lm(), 1000.0);
            close(circular.relative_stokes()[0].v, -1.0);
            for angle in [0.0, 0.37, PI / 4.0] {
                let out = circular.through(&SpectralOptics::linear_polarizer(angle).unwrap()).unwrap();
                close(out.intensity_source().flux_lm(), 500.0);
            }
        }

        #[test]
        fn finite_extinction_is_not_an_ideal_polarizer_or_depolarizer() {
            let optic = SpectralOptics::linear_diattenuator(0.0,
                [0.8; SPECTRAL_BINS], [0.2; SPECTRAL_BINS]).unwrap();
            let first = source().with_unpolarized_optics().through(&optic).unwrap();
            close(first.relative_stokes()[0].i, 0.5);
            close(first.relative_stokes()[0].q, 0.3);
            let crossed = SpectralOptics::linear_diattenuator(PI / 2.0,
                [0.8; SPECTRAL_BINS], [0.2; SPECTRAL_BINS]).unwrap();
            close(first.through(&crossed).unwrap().intensity_source().flux_lm(), 160.0);
            let linear = [Stokes::try_new(1.0, 0.0, 1.0, 0.0).unwrap(); SPECTRAL_BINS];
            let out = PolarizedSource::try_new(source(), linear).unwrap().through(&optic).unwrap();
            close(out.relative_stokes()[0].u, 0.4);
            close(out.relative_stokes()[0].degree(), 1.0);
        }

        #[test]
        fn neutral_filter_matches_native_spectral_calibration() {
            let transmission = core::array::from_fn(|bin| (bin + 1) as f64 / SPECTRAL_BINS as f64);
            let original = source();
            let expected = original.filtered(transmission).unwrap();
            let out = original.with_unpolarized_optics()
                .through(&SpectralOptics::neutral(transmission).unwrap()).unwrap();
            assert_eq!(out.intensity_source().flux_lm(), expected.flux_lm());
            for bin in 0..SPECTRAL_BINS {
                let wavelength = wavelength_nm(bin);
                assert_eq!(out.intensity_source().spectrum().value_w_per_lumen_nm(wavelength),
                    expected.spectrum().value_w_per_lumen_nm(wavelength));
            }
        }

        #[test]
        fn dispersive_fresnel_channels_close_calibrated_luminous_power() {
            let air = CauchyIor::try_constant(1.0).unwrap();
            let material = crate::dielectric::DielectricGlass::representative_crown();
            let sampled = DielectricIor::from(
                TabulatedIor::try_new(core::array::from_fn(|bin| {
                    DielectricIorSample::try_new(
                        380.0 + 50.0 * bin as f64,
                        1.6 - 0.005 * bin as f64,
                    )
                    .unwrap()
                }))
                .unwrap(),
            );
            let source = source().with_unpolarized_optics();
            for glass in [material.ior(), &sampled] {
                let r = source
                    .through(
                        &SpectralOptics::dielectric_channel(
                            0.6,
                            air,
                            glass,
                            DielectricEvent::Reflection,
                        )
                        .unwrap(),
                    )
                    .unwrap();
                let t = source
                    .through(
                        &SpectralOptics::dielectric_channel(
                            0.6,
                            air,
                            glass,
                            DielectricEvent::Transmission,
                        )
                        .unwrap(),
                    )
                    .unwrap();
                close(
                    r.intensity_source().flux_lm() + t.intensity_source().flux_lm(),
                    1000.0,
                );
                assert!(r.relative_stokes()[0].q > 0.0);
                assert_ne!(
                    r.relative_stokes()[0].i,
                    r.relative_stokes()[SPECTRAL_BINS - 1].i
                );
            }
        }

        #[test]
        fn waveplate_uses_principal_spectral_indices_and_metric_thickness() {
            let a = CauchyIor::try_constant(1.5).unwrap();
            let b = CauchyIor::try_constant(1.6).unwrap();
            let bin = 30;
            let thickness = wavelength_nm(bin) * 1e-9 / (4.0 * (1.6 - 1.5));
            let plate = SpectralOptics::waveplate(0.0, a, b, thickness).unwrap();
            let a_sellmeier = SellmeierIor::try_new([1.25, 0.0, 0.0], [0.0; 3]).unwrap();
            let b_sampled = DielectricIor::from(
                TabulatedIor::try_new(core::array::from_fn(|i| {
                    DielectricIorSample::try_new(380.0 + 50.0 * i as f64, 1.6).unwrap()
                }))
                .unwrap(),
            );
            assert_eq!(
                plate.matrices,
                SpectralOptics::waveplate(0.0, a_sellmeier, &b_sampled, thickness)
                    .unwrap()
                    .matrices
            );
            let linear = [Stokes::try_new(1.0, 0.0, 1.0, 0.0).unwrap(); SPECTRAL_BINS];
            let out = PolarizedSource::try_new(source(), linear)
                .unwrap()
                .through(&plate)
                .unwrap();
            close(out.relative_stokes()[bin].v, -1.0);
            close(out.intensity_source().flux_lm(), 1000.0);
            assert!(out.relative_stokes()[0].u.abs() > 0.01);
        }

        #[test]
        fn invalid_filters_and_brightness_redefinition_are_refused() {
            assert!(SpectralOptics::linear_polarizer(f64::NAN).is_err());
            assert!(SpectralOptics::neutral([1.01; SPECTRAL_BINS]).is_err());
            assert!(SpectralOptics::linear_retarder(0.0, [f64::INFINITY; SPECTRAL_BINS]).is_err());
            assert!(PolarizedSource::try_new(source(), [Stokes::unpolarized(2.0); SPECTRAL_BINS]).is_err());
            let invalid = Stokes { i: 1.0, q: 2.0, u: 0.0, v: 0.0 };
            assert!(PolarizedSource::try_new(source(), [invalid; SPECTRAL_BINS]).is_err());
            let dark = source().with_unpolarized_optics()
                .through(&SpectralOptics::neutral([0.0; SPECTRAL_BINS]).unwrap()).unwrap()
                .through(&SpectralOptics::linear_polarizer(0.3).unwrap()).unwrap();
            assert_eq!(dark.intensity_source().flux_lm(), 0.0);
            assert!(dark.relative_stokes().iter().all(|s| *s == Stokes::unpolarized(0.0)));
        }
    }
}
