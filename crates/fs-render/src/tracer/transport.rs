//! Reciprocal transport bookkeeping shared by camera and light subpaths.
//!
//! The legacy path tracer only needed the density of the direction it had
//! just sampled. Bidirectional techniques additionally need the density of the
//! same edge when the path is generated in reverse, and both densities must be
//! expressed in a common measure before MIS. This module owns those semantics;
//! it deliberately contains no scene- or Euler-specific policy.

mod polarized_brushed_conductor;

/// Direction in which a scattering function transports throughput.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum TransportMode {
    /// Sensor-to-scene paths carry radiance and include the refractive
    /// `eta_i^2 / eta_t^2` Jacobian on transmission.
    Radiance,
    /// Emitter-to-scene paths carry importance. The corresponding refractive
    /// Jacobian is supplied when the light and camera techniques are compared.
    Importance,
}

/// Forward and reverse densities of one sampled scattering event.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct DirectionalPdfPair {
    forward_solid_angle: f64,
    reverse_solid_angle: f64,
    delta: bool,
}

impl DirectionalPdfPair {
    /// Admit one pair of finite solid-angle densities.
    pub(super) fn continuous(forward_solid_angle: f64, reverse_solid_angle: f64) -> Option<Self> {
        if valid_pdf(forward_solid_angle) && valid_pdf(reverse_solid_angle) {
            Some(Self {
                forward_solid_angle,
                reverse_solid_angle,
                delta: false,
            })
        } else {
            None
        }
    }

    /// A delta event has no finite solid-angle density. Its discrete branch
    /// probability belongs to path throughput, not an area-measure MIS PDF.
    pub(super) const fn delta() -> Self {
        Self {
            forward_solid_angle: 0.0,
            reverse_solid_angle: 0.0,
            delta: true,
        }
    }

    pub(super) const fn forward_solid_angle(self) -> f64 {
        self.forward_solid_angle
    }

    pub(super) const fn reverse_solid_angle(self) -> f64 {
        self.reverse_solid_angle
    }

    pub(super) const fn is_delta(self) -> bool {
        self.delta
    }

    /// Convert both directions to the area measures at their respective edge
    /// targets. For an edge `x -> y`, the forward density uses the cosine at
    /// `y`; the reverse density uses the cosine at `x`.
    pub(super) fn to_area(
        self,
        distance_squared: f64,
        target_abs_cosine: f64,
        source_abs_cosine: f64,
    ) -> Option<AreaPdfPair> {
        if self.delta {
            return Some(AreaPdfPair::delta());
        }
        if !distance_squared.is_finite()
            || distance_squared <= 0.0
            || !valid_cosine(target_abs_cosine)
            || !valid_cosine(source_abs_cosine)
        {
            return None;
        }
        let forward_area = self.forward_solid_angle * target_abs_cosine / distance_squared;
        let reverse_area = self.reverse_solid_angle * source_abs_cosine / distance_squared;
        if valid_pdf(forward_area) && valid_pdf(reverse_area) {
            Some(AreaPdfPair {
                forward_area,
                reverse_area,
                delta: false,
            })
        } else {
            None
        }
    }
}

/// Forward and reverse densities after solid-angle-to-area conversion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(super) struct AreaPdfPair {
    forward_area: f64,
    reverse_area: f64,
    delta: bool,
}

impl AreaPdfPair {
    const fn delta() -> Self {
        Self {
            forward_area: 0.0,
            reverse_area: 0.0,
            delta: true,
        }
    }

    pub(super) const fn forward_area(self) -> f64 {
        self.forward_area
    }

    pub(super) const fn reverse_area(self) -> f64 {
        self.reverse_area
    }

    pub(super) const fn is_delta(self) -> bool {
        self.delta
    }
}

/// Refractive throughput factor for the declared transport direction.
pub(super) fn refractive_transport_factor(
    mode: TransportMode,
    eta_incident: f64,
    eta_transmitted: f64,
) -> Option<f64> {
    if !eta_incident.is_finite()
        || !eta_transmitted.is_finite()
        || eta_incident <= 0.0
        || eta_transmitted <= 0.0
    {
        return None;
    }
    let factor = match mode {
        TransportMode::Radiance => {
            let ratio = eta_incident / eta_transmitted;
            ratio * ratio
        }
        TransportMode::Importance => 1.0,
    };
    factor.is_finite().then_some(factor)
}

fn valid_pdf(pdf: f64) -> bool {
    pdf.is_finite() && pdf >= 0.0
}

fn valid_cosine(cosine: f64) -> bool {
    cosine.is_finite() && (0.0..=1.0).contains(&cosine)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn g0_area_measure_conversion_is_reciprocal_under_edge_reversal() {
        let pair = DirectionalPdfPair::continuous(0.25, 0.75).unwrap();
        let forward = pair.to_area(4.0, 0.5, 0.8).unwrap();
        let reverse = DirectionalPdfPair::continuous(0.75, 0.25)
            .unwrap()
            .to_area(4.0, 0.8, 0.5)
            .unwrap();
        assert_eq!(
            forward.forward_area().to_bits(),
            reverse.reverse_area().to_bits()
        );
        assert_eq!(
            forward.reverse_area().to_bits(),
            reverse.forward_area().to_bits()
        );
        assert!(!forward.is_delta());
    }

    #[test]
    fn g0_delta_event_never_masquerades_as_finite_area_density() {
        let area = DirectionalPdfPair::delta().to_area(2.0, 0.4, 0.6).unwrap();
        assert!(area.is_delta());
        assert_eq!(area.forward_area().to_bits(), 0.0_f64.to_bits());
        assert_eq!(area.reverse_area().to_bits(), 0.0_f64.to_bits());
    }

    #[test]
    fn g0_radiance_eta_factors_cancel_across_a_lossless_slab() {
        let entry = refractive_transport_factor(TransportMode::Radiance, 1.0, 1.5).unwrap();
        let exit = refractive_transport_factor(TransportMode::Radiance, 1.5, 1.0).unwrap();
        assert!((entry * exit - 1.0).abs() <= 2.0 * f64::EPSILON);
        assert_eq!(
            refractive_transport_factor(TransportMode::Importance, 1.0, 1.5)
                .unwrap()
                .to_bits(),
            1.0_f64.to_bits()
        );
    }
}

// This adapter lives with the native proposal owner to reuse its private
// kernels, not to introduce another GGX distribution or change scalar paths.
mod polarized_conductor {
    use fs_geom::Vec3;
    use fs_math::det;
    use crate::conductor::{ConductorError, ConductorOptics, ConductorSurface};
    use crate::dielectric::DielectricEvent;
    use crate::polarization::Stokes;
    use crate::polarization::dielectric::mirrors::ConductorBeamError;
    use crate::polarization::dielectric::transport::{
        BeamError, BeamSample, PolarizationFrame, PolarizedBeam,
    };
    use super::super::{sample_ggx_vndf_reflection, smith_g1};

    impl PolarizedBeam {
        /// Sample an isotropic rough metal with the native visible-normal GGX
        /// proposal. The two explicit uniform samples must lie in `[0,1)`.
        /// Invalid samples return `ConductorError::InvalidSample`, not a miss.
        ///
        /// The geometric normal is normalized and face-forwarded. No separate
        /// tilted shading normal is inferred. The existing conductor optical
        /// owner evaluates bare, transparent-coated or absorbing-coated Fresnel
        /// phase at the sampled facet, using the beam's actual phase index.
        ///
        /// For this proposal, `BSDF * cos / pdf = G1(outgoing) * Fresnel`.
        /// Cancel the distribution and incoming masking analytically rather
        /// than divide by a near-zero scalar reflectance or subtract absorbed
        /// power. The returned PDF is the unchanged native solid-angle proposal;
        /// only the direction is normalized into the physical beam frame.
        ///
        /// `Ok(None)` retains a native below-surface facet rejection. Never
        /// resample it or normalize the PDF by the accepted fraction. This is
        /// the original single-scattering GGX model, not inter-facet multiple
        /// scattering, a thermal deposition law or a brushed-surface model.
        pub fn sample_rough_conductor(
            self,
            normal: Vec3,
            optics: &ConductorOptics,
            surface: ConductorSurface,
            samples: [f64; 2],
        ) -> Result<Option<BeamSample>, ConductorBeamError> {
            if samples.iter().any(|u| !u.is_finite() || !(0.0..1.0).contains(u)) {
                return Err(ConductorError::InvalidSample.into());
            }
            let mut normal = normalized(normal)?;
            let wo = self.frame().direction().scale(-1.0);
            if normal.dot(wo) < 0.0 { normal = normal.scale(-1.0); }
            if normal.dot(wo) <= 0.0 { return Err(BeamError::InvalidFrame.into()); }
            let alpha = surface.roughness_alpha();
            let Some((direction, pdf)) = sample_ggx_vndf_reflection(
                normal, wo, alpha, samples[0], samples[1],
            ) else { return Ok(None); };
            if !pdf.is_finite() || pdf <= 0.0 {
                return Err(BeamError::NumericalRange.into());
            }
            let facet = normalized(Vec3::new(
                wo.x + direction.x, wo.y + direction.y, wo.z + direction.z,
            ))?;
            // Reuse the one frame-aware conductor response; do not evaluate
            // its scalar Fresnel law again or divide by its mean reflectance.
            let reflected = self.reflect_conductor(facet, optics)?.beam;
            let mask = smith_g1(alpha, normal.dot(direction).clamp(0.0, 1.0));
            if !mask.is_finite() || !(0.0..=1.0).contains(&mask) {
                return Err(BeamError::NumericalRange.into());
            }
            let s = reflected.power();
            let power = Stokes { i: s.i * mask, q: s.q * mask, u: s.u * mask, v: s.v * mask };
            let frame = PolarizationFrame::try_new(direction, reflected.frame().horizontal())?;
            let beam = Self::try_new(
                self.ray().origin, frame, self.wavelength_nm(), self.medium(), power,
            )?;
            Ok(Some(BeamSample { beam, event: DielectricEvent::Reflection, pdf, delta: false }))
        }
    }

    fn normalized(v: Vec3) -> Result<Vec3, BeamError> {
        if ![v.x, v.y, v.z].iter().all(|x| x.is_finite()) {
            return Err(BeamError::InvalidFrame);
        }
        let scale = v.x.abs().max(v.y.abs()).max(v.z.abs());
        if scale == 0.0 { return Err(BeamError::InvalidFrame); }
        let v = Vec3::new(v.x / scale, v.y / scale, v.z / scale);
        Ok(v.scale(1.0 / det::sqrt(v.dot(v))))
    }

    #[cfg(test)]
    mod tests {
        use super::*;
        use fs_geom::Point3;
        use crate::conductor::{CoatingIorSample, CoherentAbsorbingCoating, CoherentDielectricCoating,
            ConductorDataStatus, ConductorIorSample, ConductorSource};
        use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
        use crate::tracer::{Material, basis_all_sphere, opaque_bsdf_eval};

        fn medium() -> DielectricGlass {
            DielectricGlass::new(CauchyIor::try_constant(1.33).unwrap(),
                BeerLambertAbsorption::CLEAR, GlassProvenance::Custom)
        }
        fn bare(eta: f64, k: f64) -> ConductorOptics {
            let source = ConductorSource::try_new(fs_blake3::ContentHash([91; 32]),
                ConductorDataStatus::Representative).unwrap();
            ConductorOptics::try_new(core::array::from_fn(|i|
                ConductorIorSample::try_new(380.0 + 50.0 * i as f64, eta, k).unwrap()), source).unwrap()
        }
        fn optics() -> [ConductorOptics; 3] {
            let metal = bare(0.2, 3.0);
            let transparent = metal.with_coating(CoherentDielectricCoating::try_new(
                CauchyIor::try_constant(1.45).unwrap(), 125e-9, metal.source()).unwrap());
            let absorbing = metal.with_absorbing_coating(CoherentAbsorbingCoating::try_new(
                core::array::from_fn(|i| CoatingIorSample::try_new(380.0+50.0*i as f64,1.45,0.07).unwrap()),
                100e-9, metal.source()).unwrap());
            [metal, transparent, absorbing]
        }
        fn beam(wo: Vec3, axis: Vec3, state: Stokes) -> PolarizedBeam {
            PolarizedBeam::try_new(Point3::new(0.0,0.0,0.0),
                PolarizationFrame::try_new(wo.scale(-1.0),axis).unwrap(),550.0,medium(),state).unwrap()
        }
        fn near(a: f64, b: f64) {
            assert!((a-b).abs() <= 2e-10*a.abs().max(b.abs()).max(1.0), "{a:.17e} != {b:.17e}");
        }

        #[test]
        fn g1_rough_metals_retain_native_proposals_and_scalar_power() {
            let mut accepted=0; let mut rejected=0;
            for raw in [Vec3::new(0.0,0.0,1.0), Vec3::new(0.0,0.0,-1.0), Vec3::new(2e-5,-3e-5,-1.0)] {
                let n=normalized(raw).unwrap(); let (t,b)=basis_all_sphere(n);
                for cosine in [0.1_f64,0.6,1.0] {
                    let sine=(1.0-cosine*cosine).sqrt();
                    let wo=normalized(Vec3::new(n.x*cosine+t.x*sine,n.y*cosine+t.y*sine,n.z*cosine+t.z*sine)).unwrap();
                    let source=beam(wo,b,Stokes::unpolarized(2.0));
                    let wo=source.frame().direction().scale(-1.0);
                    for alpha in [0.02,0.2,1.0] {
                        let surface=ConductorSurface::try_rough(alpha).unwrap();
                        for law in optics() {
                            let material=Material::Conductor { optics:law, surface };
                            for draw in [[0.0,0.0],[0.17,0.63],[0.9,0.2],[0.999,0.7]] {
                                let native=sample_ggx_vndf_reflection(n,wo,alpha,draw[0],draw[1]);
                                let actual=source.sample_rough_conductor(n,&law,surface,draw).unwrap();
                                match (native,actual) {
                                    (None,None)=>rejected+=1,
                                    (Some((wi,pdf)),Some(out))=> {
                                        accepted+=1;
                                        assert_eq!(out.pdf,pdf); assert!(!out.delta);
                                        assert_eq!(out.event,DielectricEvent::Reflection);
                                        let expected=opaque_bsdf_eval(&material,n,wo,wi,550.0,Some(medium()),None).unwrap()
                                            * n.dot(wi)/pdf*2.0;
                                        near(out.beam.power().i,expected);
                                        for (a,b) in [wi.x,wi.y,wi.z].into_iter().zip([
                                            out.beam.ray().dir.x,out.beam.ray().dir.y,out.beam.ray().dir.z]) { near(a,b); }
                                        assert_eq!(out.beam.medium(),source.medium());
                                        assert!(out.beam.power().i <= source.power().i+1e-13);
                                    }
                                    other=>panic!("native proposal support changed: {other:?}"),
                                }
                            }
                        }
                    }
                }
            }
            assert!(accepted>0 && rejected>0);
        }

        #[test]
        fn g3_rough_metal_reference_covariance_preserves_partial_and_circular_states() {
            let n=Vec3::new(0.0,0.0,1.0);
            let wo=Vec3::new(0.6,0.0,0.8);
            let surface=ConductorSurface::try_rough(0.2).unwrap();
            for state in [Stokes::try_new(2.0,0.6,-0.8,0.4).unwrap(),Stokes::try_new(1.0,0.0,0.0,1.0).unwrap()] {
                let source=beam(wo,Vec3::new(0.0,1.0,0.0),state);
                let changed=source.with_reference_axis(Vec3::new(1.0,2.0,3.0)).unwrap();
                for law in optics() {
                    let a=source.sample_rough_conductor(n,&law,surface,[0.17,0.63]).unwrap().unwrap();
                    let b=changed.sample_rough_conductor(n.scale(-2.0),&law,surface,[0.17,0.63]).unwrap().unwrap();
                    for (a,b) in [a.beam.power().i,a.beam.power().q,a.beam.power().u,a.beam.power().v]
                        .into_iter().zip([b.beam.power().i,b.beam.power().q,b.beam.power().u,b.beam.power().v]) { near(a,b); }
                    assert_eq!(a.pdf,b.pdf);
                    assert!(Stokes::try_new(a.beam.power().i,a.beam.power().q,a.beam.power().u,a.beam.power().v).is_some());
                }
            }
        }

        #[test]
        fn g0_rough_metal_rejections_are_not_resampled_or_admission_failures() {
            let source=beam(Vec3::new(0.0,0.0,1.0),Vec3::new(1.0,0.0,0.0),Stokes::unpolarized(1.0));
            let original=source; let law=bare(0.2,3.0); let surface=ConductorSurface::try_rough(1.0).unwrap();
            let n=Vec3::new(0.0,0.0,1.0);
            assert_eq!(source.sample_rough_conductor(n,&law,surface,[0.9,0.2]).unwrap(),None);
            for draw in [[-0.1,0.0],[1.0,0.0],[0.0,f64::NAN],[0.0,f64::INFINITY]] {
                assert_eq!(source.sample_rough_conductor(n,&law,surface,draw),
                    Err(ConductorBeamError::Optics(ConductorError::InvalidSample)));
            }
            for normal in [Vec3::new(0.0,0.0,0.0),Vec3::new(1.0,0.0,0.0),Vec3::new(f64::NAN,0.0,1.0)] {
                assert!(source.sample_rough_conductor(normal,&law,surface,[0.0,0.0]).is_err());
            }
            assert_eq!(source,original);
        }

        #[test]
        fn g1_rough_coating_null_does_not_divide_by_reflectance() {
            let metal=bare(2.25,1e-8);
            let law=metal.with_coating(CoherentDielectricCoating::try_new(
                CauchyIor::try_constant(1.5).unwrap(),100e-9,metal.source()).unwrap());
            let vacuum=DielectricGlass::new(CauchyIor::try_constant(1.0).unwrap(),
                BeerLambertAbsorption::CLEAR,GlassProvenance::Custom);
            let source=PolarizedBeam::try_new(Point3::new(0.0,0.0,0.0),
                PolarizationFrame::try_new(Vec3::new(0.0,0.0,-1.0),Vec3::new(1.0,0.0,0.0)).unwrap(),
                600.0,vacuum,Stokes::try_new(1.0,0.0,0.0,1.0).unwrap()).unwrap();
            for alpha in [1e-4,0.1,1.0] {
                let out=source.sample_rough_conductor(Vec3::new(0.0,0.0,1.0),&law,
                    ConductorSurface::try_rough(alpha).unwrap(),[0.0,0.0]).unwrap().unwrap();
                assert!(out.beam.power().i>0.0 && out.beam.power().i<1e-16);
                assert!((out.beam.power().v+out.beam.power().i).abs()<1e-28);
                assert!(out.pdf.is_finite() && out.pdf>0.0);
            }
        }
    }
}
