use super::*;
use crate::thermal_surface::{ThermalSurfaceError, planck_radiance_w_m2_sr_nm};

const TEMPERATURE_K: f64 = 2400.0;

fn thermal_source(medium: usize) -> IsothermalVolumeEmitter {
    IsothermalVolumeEmitter::try_new(medium, TEMPERATURE_K, ContentHash([91; 32])).unwrap()
}

#[test]
fn g1_isothermal_slab_uses_true_absorption_with_stable_optically_thin_emission() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(9, 1.0, 0.0, 10.0)]).unwrap();
        let blackbody = planck_radiance_w_m2_sr_nm(TEMPERATURE_K, 550.0).unwrap();
        for index in [1.0, 1.5] {
            for absorption in [0.0_f64, 1e-12, 0.4, 30.0] {
                let material = glass(index, absorption);
                let media = [SceneMedium::Homogeneous(&material)];
                let native = CurvedPolarizedScene::try_new(
                    cx, &geometry, &media, vec![opaque(9)], 100_000,
                ).unwrap();
                let cold = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap()
                    .sample_ray(cx, ray(0.0, 1.0), 550.0, 0, settings(), 31, 2).unwrap();
                let camera = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap()
                    .with_isothermal_volume_emitters(cx, &[thermal_source(0)]).unwrap();
                let hot = camera.sample_ray(cx, ray(0.0, 1.0), 550.0, 0, settings(), 31, 2).unwrap();
                assert_eq!(cold.path(), hot.path(), "source must not change the traced path");
                assert_eq!(cold.transport_weight(), hot.transport_weight());
                let expected = -(-absorption).exp_m1();
                near(hot.radiance_w_m2_sr_nm() / (index * index * blackbody), expected,
                    5e-7 * expected);
                if absorption == 0.0 { assert_eq!(cold, hot); }
                else { assert!(hot.radiance_w_m2_sr_nm() > 0.0); }
            }
        }
    });
}

struct WarmGradient { group_offset: f64 }
impl GradientIndexField for WarmGradient {
    fn sample(&self, p: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
        let n = 1.2 + 0.2 * p.x;
        GradientIndexSample::try_new(n, Vec3::new(0.2, 0.0, 0.0), 0.3 + p.x,
            Some(n + self.group_offset))
    }
}

#[test]
fn g1_grin_lte_integrates_spatial_absorption_on_the_curved_arc_without_group_index_scaling() {
    with_cx(|cx, _| {
        // n(x)=n0+g*x: s(z)=n0/g*sinh(g*z/n0). Integrate alpha(x)=.3+x
        // on that analytic curve, independently of the native RK integral.
        let n0 = 1.2_f64;
        let g = 0.2_f64;
        let u = g / n0;
        let arc = n0 / g * u.sinh();
        let x = n0 / g * (u.cosh() - 1.0);
        let depth = 0.3 * arc + (n0 / g).powi(2)
            * (0.5 * u.sinh() * u.cosh() + 0.5 * u - u.sinh());
        assert!((depth - 0.3).abs() > 0.02, "a straight chord has the wrong opacity");
        let geometry = InstanceScene::try_new(vec![plane(9, 1.0, x, 0.01)]).unwrap();
        for group_offset in [0.0, 1.3] {
            let field = WarmGradient { group_offset };
            let media = [SceneMedium::GradientIndex {
                field: &field, identity: ContentHash([92; 32]),
            }];
            let native = CurvedPolarizedScene::try_new(
                cx, &geometry, &media, vec![opaque(9)], 100_000,
            ).unwrap();
            let camera = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap()
                .with_isothermal_volume_emitters(cx, &[thermal_source(0)]).unwrap();
            let sample = camera.sample_ray(cx, ray(0.0, 1.0), 550.0, 0, settings(), 31, 2).unwrap();
            let blackbody = planck_radiance_w_m2_sr_nm(TEMPERATURE_K, 550.0).unwrap();
            near(sample.path().traveled_m(), arc, 2e-7);
            near(sample.path().optical_depth(), depth, 2e-7);
            near(sample.radiance_w_m2_sr_nm() / blackbody,
                n0 * n0 * -(-depth).exp_m1(), 3e-7);
            let mut limited = settings();
            limited.max_distance_m = 0.5;
            assert!(matches!(camera.sample_ray(cx, ray(0.0, 1.0), 550.0, 0, limited, 31, 2),
                Err(CameraRadianceError::Unresolved(PathTermination::DistanceLimit))));
        }
    });
}

#[test]
fn g1_refractive_lte_boundary_has_one_camera_phase_factor_and_one_fresnel_probability() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![
            plane(1, 0.4, 0.0, 10.0), plane(9, 1.0, 0.0, 10.0),
        ]).unwrap();
        let outside = glass(1.2, 0.1);
        let inside = glass(1.7, 0.4);
        let media = [SceneMedium::Homogeneous(&outside), SceneMedium::Homogeneous(&inside)];
        let native = CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![
            SurfaceBinding { object_id: 1, interaction: Interaction::Dielectric {
                inside: 0, outside: 1, surface: DielectricSurface::SMOOTH,
            }}, opaque(9),
        ], 100_000).unwrap();
        let camera = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap()
            .with_isothermal_volume_emitters(cx, &[thermal_source(1)]).unwrap();
        let trace = |branch| camera.sample_ray_with_samples(cx, ray(0.0, 1.0), 550.0,
            0, settings(), |_| [0.5, 0.5, branch]).unwrap();
        let reflection = trace(0.0);
        let transmission = trace(0.99);
        let blackbody = planck_radiance_w_m2_sr_nm(TEMPERATURE_K, 550.0).unwrap();
        let conditional = 1.2_f64.powi(2) * (-0.04_f64).exp() * -(-0.24_f64).exp_m1();
        assert_eq!(reflection.radiance_w_m2_sr_nm(), 0.0);
        near(transmission.radiance_w_m2_sr_nm() / blackbody, conditional, 4e-7);
        let r = ((1.7_f64 - 1.2) / (1.7 + 1.2)).powi(2);
        near((r * reflection.radiance_w_m2_sr_nm()
            + (1.0 - r) * transmission.radiance_w_m2_sr_nm()) / blackbody,
            (1.0 - r) * conditional, 4e-7);
    });
}

#[test]
fn g1_infinite_absorbing_medium_has_finite_refractive_blackbody_equilibrium() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(9, -1.0, 0.0, 10.0)]).unwrap();
        for absorption in [0.0, 1e-12, 0.7] {
            let material = glass(1.5, absorption);
            let media = [SceneMedium::Homogeneous(&material)];
            let native = CurvedPolarizedScene::try_new(cx, &geometry, &media,
                vec![opaque(9)], 100_000).unwrap();
            let camera = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap()
                .with_isothermal_volume_emitters(cx, &[thermal_source(0)]).unwrap();
            let sample = camera.sample_ray(cx, ray(0.0, 1.0), 550.0, 0, settings(), 31, 2).unwrap();
            assert_eq!(sample.path().termination(), PathTermination::Escaped);
            let blackbody = planck_radiance_w_m2_sr_nm(TEMPERATURE_K, 550.0).unwrap();
            near(sample.radiance_w_m2_sr_nm() / blackbody,
                if absorption == 0.0 { 0.0 } else { 2.25 }, 1e-14);
        }
    });
}

#[test]
fn g0_lte_bindings_retain_sources_refuse_ambiguity_and_preserve_disabled_bits() {
    assert_eq!(IsothermalVolumeEmitter::try_new(0, 0.0, ContentHash([1; 32])),
        Err(ThermalSurfaceError::Temperature));
    assert_eq!(IsothermalVolumeEmitter::try_new(0, TEMPERATURE_K, ContentHash([0; 32])),
        Err(ThermalSurfaceError::SourceIdentity));
    assert_eq!(thermal_source(0).source_identity(), ContentHash([91; 32]));
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(9, 1.0, 0.0, 10.0)]).unwrap();
        let material = glass(1.4, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let native = CurvedPolarizedScene::try_new(cx, &geometry, &media,
            vec![opaque(9)], 100_000).unwrap();
        let sources = [emitter(9, 3.0, EmissionSides::Both)];
        let camera = || CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        assert!(matches!(camera().with_isothermal_volume_emitters(cx, &[thermal_source(1)]),
            Err(CameraRadianceError::InvalidVolumeEmitter { medium_index: 1 })));
        assert!(camera().with_isothermal_volume_emitters(cx,
            &[thermal_source(0), thermal_source(0)]).is_err());
        let disabled = camera().with_isothermal_volume_emitters(cx, &[thermal_source(0)]).unwrap()
            .with_isothermal_volume_emitters(cx, &[]).unwrap();
        assert_eq!(camera().sample_ray(cx, ray(0.0, 1.0), 550.0, 0, settings(), 31, 2),
            disabled.sample_ray(cx, ray(0.0, 1.0), 550.0, 0, settings(), 31, 2));
    });
}

#[test]
fn g1_thermal_pixels_use_native_spectral_lanes_and_match_zero_scattering() {
    use crate::camera::{Aperture, CameraProjection, PhysicalCamera};
    use super::super::raster::CurvedCameraRaster;
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(9, 1.0, 0.0, 10.0)]).unwrap();
        let material = glass(1.4, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let native = CurvedPolarizedScene::try_new(cx, &geometry, &media,
            vec![opaque(9)], 100_000).unwrap();
        let camera = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap()
            .with_isothermal_volume_emitters(cx, &[thermal_source(0)]).unwrap();
        let physical = PhysicalCamera::try_look_at(Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 1.0), Vec3::new(0.0, 1.0, 0.0),
            CameraProjection::try_half_tangent(0.15).unwrap(), 1.0,
            Aperture::try_circular(0.0).unwrap()).unwrap();
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        let pixel = raster.sample_pixel(cx, &camera, 0, settings(), 17, 0, 4).unwrap();
        let depth = 0.2 / pixel.camera_ray().dir.z;
        for (&wavelength, lane) in pixel.wavelengths_nm().iter().zip(pixel.spectral_samples()) {
            let blackbody = planck_radiance_w_m2_sr_nm(TEMPERATURE_K, wavelength).unwrap();
            near(lane.radiance_w_m2_sr_nm() / blackbody, 1.4_f64.powi(2) * -(-depth).exp_m1(), 2e-7);
        }
        #[cfg(feature = "volumes")]
        {
            use super::super::participating::CurvedRayleighCameraScene;
            use crate::polarization::dielectric::scene::curved::participating::{
                CurvedRayleighScene, CurvedRayleighTraceSettings,
            };
            let laws = [None];
            let scattering = CurvedRayleighScene::try_new(cx, &native, &laws).unwrap();
            let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &scattering).unwrap();
            let limits = CurvedRayleighTraceSettings { surfaces: settings(), max_volume_events: 8 };
            let sample = raster.sample_pixel_rayleigh(cx, &molecular, 0, limits, 17, 0, 4).unwrap();
            assert_eq!(pixel, *sample.pixel_sample());
            assert_eq!(*sample.free_flight_draws(), [0; 4]);
            assert_eq!(camera.sample_ray(cx, ray(0.0, 1.0), 550.0, 0, settings(), 17, 4).unwrap(),
                *molecular.sample_ray(cx, ray(0.0, 1.0), 550.0, 0, limits, 17, 4).unwrap().radiance_sample());
        }
    });
}
