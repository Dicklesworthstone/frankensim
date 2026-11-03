use super::*;
use crate::thermal_surface::planck_radiance_w_m2_sr_nm;

fn source<'a>(medium: usize, field: &'a (dyn TemperatureField + Sync)) -> SpatialVolumeEmitter<'a> {
    SpatialVolumeEmitter::try_new(medium, field, ContentHash([107; 32]), 1e-11).unwrap()
}

// Test-only inverse of SI Planck: author a physical T field with analytically
// integrable B(T). Production always evaluates the shared original Planck law.
fn temperature_for_radiance(radiance: f64) -> f64 {
    let h = 6.626_070_15e-34;
    let c = 299_792_458.0_f64;
    let k = 1.380_649e-23;
    let wavelength = 550e-9_f64;
    (h * c / (wavelength * k)) / (2.0 * h * c * c / wavelength.powi(5) * 1e-9 / radiance).ln_1p()
}
fn linear_source_integral(absorption: f64, distance: f64, slope: f64) -> f64 {
    let survival = (-absorption * distance).exp();
    let i0 = -(-absorption * distance).exp_m1();
    i0 + slope * (i0 / absorption - distance * survival)
}

#[test]
fn g1_constant_spatial_temperature_matches_isothermal_and_empty_bindings_preserve_bits() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(9, 1.0, 0.0, 10.0)]).unwrap();
        let temperature = |_: Point3| Ok(2400.0);
        let blackbody = planck_radiance_w_m2_sr_nm(2400.0, 550.0).unwrap();
        for absorption in [1e-12, 0.4, 30.0] {
            let material = glass(1.5, absorption);
            let media = [SceneMedium::Homogeneous(&material)];
            let native =
                CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(9)], 100_000)
                    .unwrap();
            let original = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
            let cold = original
                .sample_ray(cx, ray(0.0, 1.0), 550.0, 0, settings(), 41, 3)
                .unwrap();
            let disabled = original
                .with_spatial_volume_emitters(cx, &[])
                .unwrap()
                .sample_ray(cx, ray(0.0, 1.0), 550.0, 0, settings(), 41, 3)
                .unwrap();
            assert_eq!(disabled, cold);
            let isothermal = CurvedCameraScene::try_new(cx, &native, &[], None)
                .unwrap()
                .with_isothermal_volume_emitters(
                    cx,
                    &[
                        IsothermalVolumeEmitter::try_new(0, 2400.0, ContentHash([106; 32]))
                            .unwrap(),
                    ],
                )
                .unwrap()
                .sample_ray(cx, ray(0.0, 1.0), 550.0, 0, settings(), 41, 3)
                .unwrap();
            let spatial = CurvedCameraScene::try_new(cx, &native, &[], None)
                .unwrap()
                .with_spatial_volume_emitters(cx, &[source(0, &temperature)])
                .unwrap()
                .sample_ray(cx, ray(0.0, 1.0), 550.0, 0, settings(), 41, 3)
                .unwrap();
            let expected = -(-absorption).exp_m1();
            near(
                spatial.radiance_w_m2_sr_nm() / (2.25 * blackbody),
                expected,
                5e-7 * expected,
            );
            near(
                spatial.radiance_w_m2_sr_nm() / blackbody,
                isothermal.radiance_w_m2_sr_nm() / blackbody,
                2e-6 * expected,
            );
        }
    });
}

struct WarmCurve {
    gradient: f64,
    group_offset: f64,
}
impl GradientIndexField for WarmCurve {
    fn sample(&self, p: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
        let n = 1.2 + self.gradient * p.x;
        GradientIndexSample::try_new(
            n,
            Vec3::new(self.gradient, 0.0, 0.0),
            0.4,
            Some(n + self.group_offset),
        )
    }
}

#[test]
fn g1_temperature_gradient_follows_true_refracted_arc_with_phase_not_group_normalization() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(9, 1.0, 0.0, 10.0)]).unwrap();
        for gradient in [0.0_f64, 0.7] {
            let temperature = |p: Point3| {
                let arc = if gradient == 0.0 {
                    p.z
                } else {
                    1.2 / gradient * (gradient * p.z / 1.2).sinh()
                };
                Ok(temperature_for_radiance(1.0 + 2.0 * arc))
            };
            let arc = if gradient == 0.0 {
                1.0
            } else {
                1.2 / gradient * (gradient / 1.2).sinh()
            };
            let expected = 1.2_f64.powi(2) * linear_source_integral(0.4, arc, 2.0);
            let mut radiance = None;
            for group_offset in [0.0, 0.8] {
                let field = WarmCurve {
                    gradient,
                    group_offset,
                };
                let media = [SceneMedium::GradientIndex {
                    field: &field,
                    identity: ContentHash([108; 32]),
                }];
                let native =
                    CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(9)], 100_000)
                        .unwrap();
                let camera = CurvedCameraScene::try_new(cx, &native, &[], None)
                    .unwrap()
                    .with_spatial_volume_emitters(cx, &[source(0, &temperature)])
                    .unwrap();
                let sample = camera
                    .sample_ray(cx, ray(0.0, 1.0), 550.0, 0, settings(), 41, 3)
                    .unwrap();
                near(sample.path().traveled_m(), arc, 3e-7);
                near(sample.radiance_w_m2_sr_nm(), expected, 2e-6);
                if let Some(previous) = radiance {
                    near(sample.radiance_w_m2_sr_nm(), previous, 2e-7);
                }
                radiance = Some(sample.radiance_w_m2_sr_nm());
            }
            if gradient != 0.0 {
                assert!(
                    (expected - 1.2_f64.powi(2) * linear_source_integral(0.4, 1.0, 2.0)).abs()
                        > 0.05,
                    "a straight optical/source march must be measurably wrong"
                );
            }
        }
    });
}

#[test]
fn g1_spatial_emission_after_dielectric_boundary_uses_original_camera_prefix_and_index() {
    with_cx(|cx, _| {
        let geometry =
            InstanceScene::try_new(vec![plane(1, 0.4, 0.0, 10.0), plane(9, 1.0, 0.0, 10.0)])
                .unwrap();
        let outside = glass(1.2, 0.1);
        let inside = glass(1.7, 0.4);
        let media = [
            SceneMedium::Homogeneous(&outside),
            SceneMedium::Homogeneous(&inside),
        ];
        let native = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::Dielectric {
                        inside: 0,
                        outside: 1,
                        surface: DielectricSurface::SMOOTH,
                    },
                },
                opaque(9),
            ],
            100_000,
        )
        .unwrap();
        let temperature = |p: Point3| Ok(temperature_for_radiance(1.0 + 2.0 * (p.z - 0.4)));
        let camera = CurvedCameraScene::try_new(cx, &native, &[], None)
            .unwrap()
            .with_spatial_volume_emitters(cx, &[source(1, &temperature)])
            .unwrap();
        let reflected = camera
            .sample_ray_with_samples(cx, ray(0.0, 1.0), 550.0, 0, settings(), |_| [0.5, 0.5, 0.0])
            .unwrap();
        let transmitted = camera
            .sample_ray_with_samples(cx, ray(0.0, 1.0), 550.0, 0, settings(), |_| {
                [0.5, 0.5, 0.99]
            })
            .unwrap();
        assert_eq!(reflected.radiance_w_m2_sr_nm(), 0.0);
        let expected = 1.2_f64.powi(2) * (-0.04_f64).exp() * linear_source_integral(0.4, 0.6, 2.0);
        near(transmitted.radiance_w_m2_sr_nm(), expected, 2e-6);
    });
}

#[test]
fn g0_spatial_sources_retain_domains_refuse_unknown_tails_and_skip_transparent_temperature() {
    with_cx(|cx, _| {
        let bounded = |p: Point3| {
            if p.z > 0.3 {
                Err(GradientIndexError::OutsideFieldDomain)
            } else {
                Ok(2400.0)
            }
        };
        for z in [-1.0, 1.0] {
            let geometry = InstanceScene::try_new(vec![plane(9, z, 0.0, 10.0)]).unwrap();
            let material = glass(1.5, 0.4);
            let media = [SceneMedium::Homogeneous(&material)];
            let native =
                CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(9)], 100_000)
                    .unwrap();
            let camera = CurvedCameraScene::try_new(cx, &native, &[], None)
                .unwrap()
                .with_spatial_volume_emitters(cx, &[source(0, &bounded)])
                .unwrap();
            let error = camera
                .sample_ray(cx, ray(0.0, 1.0), 550.0, 0, settings(), 41, 3)
                .unwrap_err();
            if z < 0.0 {
                assert_eq!(
                    error,
                    CameraRadianceError::Unresolved(PathTermination::Escaped)
                );
            } else {
                assert_eq!(
                    error,
                    CameraRadianceError::Scene(CurvedSceneError::Propagation(
                        GradientIndexError::OutsideFieldDomain
                    ))
                );
            }
        }
        let geometry = InstanceScene::try_new(vec![plane(9, 1.0, 0.0, 10.0)]).unwrap();
        let clear = glass(1.5, 0.0);
        let media = [SceneMedium::Homogeneous(&clear)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![opaque(9)], 100_000).unwrap();
        let forbidden = |_: Point3| -> Result<f64, GradientIndexError> {
            panic!("transparent source queried");
        };
        let source = source(0, &forbidden);
        let original = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
        let expected = original
            .sample_ray(cx, ray(0.0, 1.0), 550.0, 0, settings(), 41, 3)
            .unwrap();
        let spatial = original
            .with_spatial_volume_emitters(cx, &[source])
            .unwrap();
        assert_eq!(
            spatial
                .sample_ray(cx, ray(0.0, 1.0), 550.0, 0, settings(), 41, 3)
                .unwrap(),
            expected
        );
        assert!(matches!(
            CurvedCameraScene::try_new(cx, &native, &[], None)
                .unwrap()
                .with_spatial_volume_emitters(cx, &[source, source]),
            Err(CameraRadianceError::InvalidVolumeEmitter { .. })
        ));
        let isothermal =
            IsothermalVolumeEmitter::try_new(0, 2400.0, ContentHash([109; 32])).unwrap();
        assert!(matches!(
            CurvedCameraScene::try_new(cx, &native, &[], None)
                .unwrap()
                .with_isothermal_volume_emitters(cx, &[isothermal])
                .unwrap()
                .with_spatial_volume_emitters(cx, &[source]),
            Err(CameraRadianceError::InvalidVolumeEmitter { .. })
        ));
        assert!(matches!(
            spatial.with_isothermal_volume_emitters(cx, &[isothermal]),
            Err(CameraRadianceError::InvalidVolumeEmitter { .. })
        ));
    });
}
