//! Polarized dielectric optics and frame-aware ray-power transport.
//!
//! Interface matrices use +Q along s, p chosen so (s,p,propagation) is
//! right-handed, U=2 Re(E_s conj(E_p)), V=2 Im(E_s conj(E_p)), and
//! exp(-i omega t). Reflection reverses p at normal incidence. Transmission
//! matrices transfer power; the legacy sampling adapter includes its radiance
//! eta-squared factor. Only lossless isotropic interfaces are admitted, with
//! native Beer-Lambert absorption as a separate bulk effect.
//!
//! [`transport`] carries the frame and current material across successive
//! physical ray interactions. It is forward power transport, not a replacement
//! for the scalar camera-path integrator.

mod interface;
#[cfg(feature = "photometry")]
pub(super) use interface::spectral_fresnel_dielectric_polarized_borrowed;
pub use interface::{
    PolarizedDielectricSample, PolarizedFresnel, attenuate_stokes,
    fresnel_dielectric_polarized, sample_smooth_dielectric_polarized,
    spectral_fresnel_dielectric_polarized,
};
pub mod transport;
/// Native geometry traversal, explicit boundary media and polarized receivers.
pub mod scene;

/// Phase-preserving reflection from explicitly supplied smooth metal optics.
pub mod mirrors {
    //! Smooth opaque reflection in the existing forward polarized beam path.
    //!
    //! Optical tables and coherent coatings belong to `crate::conductor`.
    //! This adapter owns only ray frames, reflected direction and power weight.
    //! No transmitted metal ray, roughness, shading-normal correction, random
    //! branch decision or spatial deposition of absorbed power is inferred.

    use fs_geom::Vec3;
    use fs_math::det;
    use crate::conductor::{ConductorError, ConductorOptics};
    use crate::dielectric::DielectricEvent;
    use crate::polarization::{Mueller, Stokes};
    use super::transport::{BeamError, BeamMedium, BeamSample, PolarizationFrame, PolarizedBeam};

    /// An optical-data or beam failure, never a dark sample substituted for an error.
    #[derive(Clone, Copy, Debug, PartialEq, Eq)]
    pub enum ConductorBeamError {
        Optics(ConductorError),
        Beam(BeamError),
    }
    impl From<ConductorError> for ConductorBeamError {
        fn from(error: ConductorError) -> Self { Self::Optics(error) }
    }
    impl From<BeamError> for ConductorBeamError {
        fn from(error: BeamError) -> Self { Self::Beam(error) }
    }
    impl core::fmt::Display for ConductorBeamError {
        fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
            write!(f, "polarized conductor: {self:?}")
        }
    }
    impl core::error::Error for ConductorBeamError {}

    impl PolarizedBeam {
        /// Reflect from a smooth opaque metal at the current physical ray origin.
        ///
        /// Uses the CURRENT medium's spectral phase index (or its retained local
        /// gradient-index sample) and the supplied metal's complex-index law with
        /// any transparent or absorbing coherent coating. The normal is normalized
        /// and face-forwarded. Frames rotate automatically; normal incidence keeps
        /// the original horizontal axis. Reflection leaves the medium unchanged.
        ///
        /// All paths follow the reflected direction with true reflected power,
        /// so this delta event has probability one. Do not apply Fresnel twice
        /// or add a radiance eta factor. No thermal deposition is performed.
        pub fn reflect_conductor(
            self,
            normal: Vec3,
            optics: &ConductorOptics,
        ) -> Result<BeamSample, ConductorBeamError> {
            let k = self.frame().direction();
            let mut normal = unit(normal)?;
            if normal.dot(k) > 0.0 { normal = normal.scale(-1.0); }
            let cosine = (-normal.dot(k)).clamp(0.0, 1.0);
            if cosine == 0.0 { return Err(BeamError::InvalidFrame.into()); }
            let eta = match self.medium() {
                BeamMedium::Homogeneous(material) => {
                    material.ior().eval(self.wavelength_nm()).map_err(BeamError::from)?
                }
                BeamMedium::GradientIndex(local) => {
                    if local.position() != self.ray().origin
                        || local.wavelength_nm() != self.wavelength_nm() {
                        return Err(BeamError::MediumMismatch.into());
                    }
                    local.sample().phase_index()
                }
            };
            let matrix = optics.fresnel_mueller(self.wavelength_nm(), eta, cosine)?;
            let transverse = cross(k, normal);
            let s = if transverse.x == 0.0 && transverse.y == 0.0 && transverse.z == 0.0 {
                self.frame().horizontal()
            } else { unit(transverse)? };
            let incident = self.with_reference_axis(s)?;
            let power = apply_passive(matrix, incident.power())?;
            let direction = Vec3::new(
                k.x + 2.0 * cosine * normal.x,
                k.y + 2.0 * cosine * normal.y,
                k.z + 2.0 * cosine * normal.z,
            );
            let frame = PolarizationFrame::try_new(direction, s)?;
            let beam = PolarizedBeam::try_new(
                self.ray().origin, frame, self.wavelength_nm(), self.medium(), power,
            )?;
            Ok(BeamSample { beam, event: DielectricEvent::Reflection, pdf: 1.0, delta: true })
        }
    }

    fn unit(v: Vec3) -> Result<Vec3, BeamError> {
        if ![v.x, v.y, v.z].iter().all(|x| x.is_finite()) {
            return Err(BeamError::InvalidFrame);
        }
        let scale = v.x.abs().max(v.y.abs()).max(v.z.abs());
        if scale == 0.0 { return Err(BeamError::InvalidFrame); }
        let v = Vec3::new(v.x / scale, v.y / scale, v.z / scale);
        Ok(v.scale(1.0 / det::sqrt(v.dot(v))))
    }
    fn cross(a: Vec3, b: Vec3) -> Vec3 {
        Vec3::new(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x)
    }

    fn apply_passive(matrix: Mueller, input: Stokes) -> Result<Stokes, BeamError> {
        if input.i == 0.0 { return Ok(Stokes::unpolarized(0.0)); }
        // Normalize before mixing signed terms to avoid overflow at large power.
        let mut out = matrix.apply(Stokes {
            i: 1.0, q: input.q / input.i, u: input.u / input.i, v: input.v / input.i,
        });
        let tolerance = 128.0 * f64::EPSILON;
        let norm = det::hypot(det::hypot(out.q, out.u), out.v);
        if ![out.i, out.q, out.u, out.v, norm].iter().all(|x| x.is_finite())
            || out.i < -tolerance || out.i > 1.0 + tolerance
            || norm > out.i.max(0.0) + tolerance {
            return Err(BeamError::NumericalRange);
        }
        out.i = out.i.clamp(0.0, 1.0);
        if out.i == 0.0 { return Ok(Stokes::unpolarized(0.0)); }
        if norm > out.i {
            let ratio = out.i / norm;
            out.q *= ratio; out.u *= ratio; out.v *= ratio;
        }
        Ok(Stokes { i: out.i * input.i, q: out.q * input.i, u: out.u * input.i, v: out.v * input.i })
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use fs_geom::Point3;
        use crate::conductor::{ConductorDataStatus, ConductorIorSample, ConductorSource};
        use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};

        fn glass(index: f64) -> DielectricGlass {
            DielectricGlass::new(CauchyIor::try_constant(index).unwrap(),
                BeerLambertAbsorption::CLEAR, GlassProvenance::Custom)
        }
        fn metal() -> ConductorOptics {
            let source = ConductorSource::try_new(
                fs_blake3::hash_domain("org.frankensim.test.polarized-mirror", b"synthetic"),
                ConductorDataStatus::Representative,
            ).unwrap();
            ConductorOptics::try_new(core::array::from_fn(|i| {
                ConductorIorSample::try_new(380.0 + 50.0 * i as f64, 0.2, 3.0).unwrap()
            }), source).unwrap()
        }
        fn beam(k: Vec3, axis: Vec3, state: Stokes) -> PolarizedBeam {
            PolarizedBeam::try_new(Point3::new(1.0, 2.0, 3.0),
                PolarizationFrame::try_new(k, axis).unwrap(), 550.0, glass(1.0), state).unwrap()
        }
        fn close(a: f64, b: f64) { assert!((a - b).abs() < 5e-13, "{a} != {b}"); }
        fn same(a: Stokes, b: Stokes) {
            for (a,b) in [a.i,a.q,a.u,a.v].into_iter().zip([b.i,b.q,b.u,b.v]) { close(a,b); }
        }

        #[test]
        fn g1_metal_reflection_matches_complex_field_phase_and_retains_medium() {
            let b = beam(Vec3::new(0.75_f64.sqrt(), 0.0, -0.5), Vec3::new(0.0, -1.0, 0.0),
                Stokes::try_new(1.0, 0.0, 1.0, 0.0).unwrap());
            let out = b.reflect_conductor(Vec3::new(0.0, 0.0, 1.0), &metal()).unwrap();
            same(out.beam.power(), Stokes { i: 0.9184110846593692, q: 0.04400178561613738,
                u: -0.5749762549092700, v: 0.7148042176556463 });
            assert_eq!(out.event, DielectricEvent::Reflection);
            assert!(out.delta); assert_eq!(out.pdf, 1.0);
            assert_eq!(out.beam.medium(), b.medium());
            assert_eq!(out.beam.ray().origin, b.ray().origin);
            close(out.beam.ray().dir.x, b.ray().dir.x);
            close(out.beam.ray().dir.z, -b.ray().dir.z);
            close(out.beam.power().degree(), 1.0);
        }

        #[test]
        fn g3_mirror_is_covariant_under_reference_axis_and_normal_sign_changes() {
            let b = beam(Vec3::new(0.6, 0.0, -0.8), Vec3::new(0.0, 1.0, 0.0),
                Stokes::try_new(2.0, 0.4, -0.7, 0.2).unwrap());
            let changed = b.with_reference_axis(Vec3::new(1.0, 2.0, 3.0)).unwrap();
            let n = Vec3::new(0.0, 0.0, 1.0);
            let a = b.reflect_conductor(n, &metal()).unwrap();
            let c = changed.reflect_conductor(n.scale(-10.0), &metal()).unwrap();
            same(a.beam.power(), c.beam.power());
            assert_eq!(a.beam.ray(), c.beam.ray());
            for scale in [f64::MAX, f64::from_bits(1)] {
                same(b.reflect_conductor(n.scale(scale), &metal()).unwrap().beam.power(), a.beam.power());
            }
        }

        #[test]
        fn g1_normal_mirror_flips_circular_handedness_without_power_gain() {
            let m = metal();
            let b = beam(Vec3::new(0.0,0.0,1.0), Vec3::new(1.0,0.0,0.0),
                Stokes::try_new(1.0, 0.0, 0.0, 1.0).unwrap());
            let out = b.reflect_conductor(Vec3::new(0.0,0.0,-1.0), &m).unwrap().beam;
            let r = m.fresnel(550.0, 1.0, 1.0).unwrap();
            same(out.power(), Stokes { i:r, q:0.0, u:0.0, v:-r });
            let twice = out.reflect_conductor(Vec3::new(0.0,0.0,1.0), &m).unwrap().beam;
            same(twice.power(), Stokes { i:r*r, q:0.0, u:0.0, v:r*r });
            let large = beam(Vec3::new(0.0,0.0,1.0), Vec3::new(1.0,0.0,0.0),
                Stokes::try_new(f64::MAX, f64::MAX, 0.0, 0.0).unwrap());
            let large_out = large.reflect_conductor(Vec3::new(0.0,0.0,1.0), &m).unwrap().beam;
            assert!(large_out.power().i.is_finite());
            close(large_out.power().i / f64::MAX, r);
        }

        #[test]
        fn g0_invalid_mirror_operations_do_not_mutate_the_beam() {
            let b = beam(Vec3::new(0.0,0.0,1.0), Vec3::new(1.0,0.0,0.0), Stokes::unpolarized(1.0));
            let original = b;
            for n in [Vec3::new(0.0,0.0,0.0), Vec3::new(1.0,0.0,0.0), Vec3::new(f64::NAN,0.0,1.0)] {
                assert!(b.reflect_conductor(n, &metal()).is_err());
                assert_eq!(b, original);
            }
        }

        #[test]
        fn g1_absorbing_coating_reaches_reflected_polarization_without_double_loss() {
            use crate::conductor::{CoatingIorSample, CoherentAbsorbingCoating};
            let bare = metal();
            let film = CoherentAbsorbingCoating::try_new(
                core::array::from_fn(|i| {
                    CoatingIorSample::try_new(380.0 + 50.0 * i as f64, 1.45, 0.07).unwrap()
                }), 100.0e-9, bare.source(),
            ).unwrap();
            let optics = bare.with_absorbing_coating(film);
            let b = beam(Vec3::new(0.0, 0.0, 1.0), Vec3::new(1.0, 0.0, 0.0),
                Stokes::try_new(1.0, 0.0, 0.6, 0.8).unwrap());
            let out = b.reflect_conductor(Vec3::new(0.0, 0.0, -1.0), &optics).unwrap();
            // Independent complex Airy result, not scalar-Fresnel agreement alone.
            let r = 0.515_088_324_939_222_6;
            same(out.beam.power(), Stokes { i:r, q:0.0, u:-0.6*r, v:-0.8*r });
            assert_eq!(out.beam.medium(), b.medium());
            assert!(bare.fresnel(550.0, 1.0, 1.0).unwrap() - r > 0.4);
        }

        #[test]
        fn g1_mirror_uses_local_phase_index_without_flattening_gradient_medium() {
            use crate::gradient_index::HomogeneousIndexField;
            let material = glass(1.33);
            let position = Point3::new(1.0, 2.0, 3.0);
            let local = BeamMedium::try_gradient_at(
                &HomogeneousIndexField::new(&material), position, 550.0,
            ).unwrap();
            let b = PolarizedBeam::try_new(position,
                PolarizationFrame::try_new(Vec3::new(0.0, 0.0, 1.0), Vec3::new(1.0, 0.0, 0.0)).unwrap(),
                550.0, local, Stokes::try_new(1.0, 0.0, 0.0, 1.0).unwrap(),
            ).unwrap();
            let metal = metal();
            let out = b.reflect_conductor(Vec3::new(0.0, 0.0, -1.0), &metal).unwrap().beam;
            let r = metal.fresnel(550.0, 1.33, 1.0).unwrap();
            same(out.power(), Stokes { i:r, q:0.0, u:0.0, v:-r });
            assert_eq!(out.medium(), local);
            assert_eq!(out.advance(0.1), Err(BeamError::GradientFieldRequired));
            assert!((r - metal.fresnel(550.0, 1.0, 1.0).unwrap()).abs() > 1e-3);
        }
    }
}
