use super::*;
use crate::camera::{Aperture, CameraProjection, PhysicalCamera};
use crate::conductor::{
    CoherentDielectricCoating, ConductorDataStatus, ConductorOptics, ConductorSource,
    ConductorSurface,
};
use crate::polarization::dielectric::scene::curved::camera::raster::CurvedCameraRaster;
use crate::spectral::SpectralReflectance;
use crate::surface_detail::{BrushingPattern, SurfaceDetail};

fn coated_metal() -> ConductorOptics {
    let source =
        ConductorSource::try_new(ContentHash([91; 32]), ConductorDataStatus::Representative)
            .unwrap();
    ConductorOptics::representative_tungsten().with_coating(
        CoherentDielectricCoating::try_new(CauchyIor::try_constant(1.45).unwrap(), 125e-9, source)
            .unwrap(),
    )
}

fn diffuse() -> SurfaceBinding {
    SurfaceBinding {
        object_id: 1,
        interaction: Interaction::DepolarizingDiffuse {
            medium: 0,
            reflectance: 0,
        },
    }
}

#[test]
fn g1_reflective_emission_adds_le_before_diffuse_scattering_and_counts_black_once() {
    with_cx(|cx, _| {
        let geometry =
            InstanceScene::try_new(vec![plane(1, 0.4, 0.0, 10.0), plane(9, -0.6, 0.0, 10.0)])
                .unwrap();
        let host = glass(1.33, 0.2);
        let media = [SceneMedium::Homogeneous(&host)];
        let sources = [
            emitter(1, 2.0, EmissionSides::Both),
            emitter(9, 5.0, EmissionSides::Both),
        ];
        for rho in [0.0, 0.37, 1.0] {
            let reflectances = [SpectralReflectance::constant(rho).unwrap()];
            let native = CurvedPolarizedScene::try_new_with_materials(
                cx,
                &geometry,
                &media,
                &[],
                &[],
                &reflectances,
                vec![diffuse(), opaque(9)],
                100_000,
            )
            .unwrap();
            let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
            let mut calls = 0;
            let sample = camera
                .sample_ray_with_samples(cx, ray(0.0, 1.0), 550.0, 0, settings(), |_| {
                    calls += 1;
                    [0.0, 0.25, 0.5]
                })
                .unwrap();
            // The current surface's rho and sampled PDF cannot multiply its Le.
            let expected = 2.0 * (-0.08_f64).exp() + 5.0 * rho * (-0.28_f64).exp();
            near(sample.radiance_w_m2_sr_nm(), expected, 3e-6);
            assert_eq!(calls, u32::from(rho > 0.0));
            assert_eq!(sample.path().surface_events(), 1 + u32::from(rho > 0.0));
            // Camera source bindings have no effect on the original power API.
            let input = native
                .launch_beam(
                    Point3::new(0.0, 0.0, 0.0),
                    PolarizationFrame::try_new(Vec3::new(0.0, 0.0, 1.0), Vec3::new(1.0, 0.0, 0.0))
                        .unwrap(),
                    550.0,
                    0,
                    Stokes::unpolarized(1.0),
                )
                .unwrap();
            let power = native
                .trace_with_samples(cx, input, 0, settings(), |_| [0.0, 0.25, 0.5])
                .unwrap();
            assert_eq!(power, *sample.path());
        }
    });
}

#[test]
fn g1_reflective_emission_uses_each_incident_index_before_fresnel_branch_selection() {
    with_cx(|cx, _| {
        let geometry =
            InstanceScene::try_new(vec![plane(1, 0.4, 0.0, 10.0), plane(9, 1.0, 0.0, 10.0)])
                .unwrap();
        let initial = glass(1.0, 0.1);
        let transmitted = glass(1.5, 0.4);
        let media = [
            SceneMedium::Homogeneous(&initial),
            SceneMedium::Homogeneous(&transmitted),
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
        let sources = [
            emitter(1, 2.0, EmissionSides::Both),
            emitter(9, 5.0, EmissionSides::Both),
        ];
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let first = 2.0 * (-0.04_f64).exp();
        let last = 5.0 * (-0.28_f64).exp() / (1.5 * 1.5);
        let reflected = camera
            .sample_ray_with_samples(cx, ray(0.0, 1.0), 550.0, 0, settings(), |_| [0.5, 0.5, 0.0])
            .unwrap();
        let transmitted = camera
            .sample_ray_with_samples(cx, ray(0.0, 1.0), 550.0, 0, settings(), |_| {
                [0.5, 0.5, 0.99]
            })
            .unwrap();
        near(reflected.radiance_w_m2_sr_nm(), first, 2e-6);
        // Infinite absorbing escape kills only the environment contribution.
        assert_eq!(reflected.transport_weight(), 0.0);
        near(transmitted.radiance_w_m2_sr_nm(), first + last, 3e-6);
        near(
            0.04 * reflected.radiance_w_m2_sr_nm() + 0.96 * transmitted.radiance_w_m2_sr_nm(),
            first + 0.96 * last,
            3e-6,
        );
        let mut limited = settings();
        limited.max_surface_events = 1;
        assert!(matches!(
            camera.sample_ray_with_samples(cx, ray(0.0, 1.0), 550.0, 0, limited, |_| [
                0.5, 0.5, 0.99
            ]),
            Err(CameraRadianceError::Unresolved(
                PathTermination::SurfaceLimit
            ))
        ));
    });
}

#[test]
fn g1_reflective_emission_in_grin_has_separate_source_indices_and_beer_prefixes() {
    struct VerticalIndex;
    impl GradientIndexField for VerticalIndex {
        fn sample(&self, p: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
            let n = 1.2 + 0.02 * p.z;
            GradientIndexSample::try_new(n, Vec3::new(0.0, 0.0, 0.02), 0.1, Some(n))
        }
    }
    with_cx(|cx, _| {
        let geometry =
            InstanceScene::try_new(vec![plane(1, 0.4, 0.0, 10.0), plane(9, -0.6, 0.0, 10.0)])
                .unwrap();
        let media = [SceneMedium::GradientIndex {
            field: &VerticalIndex,
            identity: ContentHash([93; 32]),
        }];
        let metals = [coated_metal()];
        let native = CurvedPolarizedScene::try_new_with_conductors(
            cx,
            &geometry,
            &media,
            &metals,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::SmoothConductor {
                        optics: 0,
                        medium: 0,
                    },
                },
                opaque(9),
            ],
            100_000,
        )
        .unwrap();
        let sources = [
            emitter(1, 2.0, EmissionSides::Both),
            emitter(9, 3.0, EmissionSides::Both),
        ];
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let sample = camera
            .sample_ray_with_samples(cx, ray(0.0, 1.0), 550.0, 0, settings(), |_| {
                panic!("delta mirror and terminal emitter need no draws")
            })
            .unwrap();
        let fresnel = metals[0].fresnel(550.0, 1.208, 1.0).unwrap();
        let first = 2.0 * (-0.04_f64).exp() * (1.2_f64 / 1.208).powi(2);
        let last = 3.0 * fresnel * (-0.14_f64).exp() * (1.2_f64 / 1.188).powi(2);
        near(sample.radiance_w_m2_sr_nm(), first + last, 3e-6);
        near(sample.transport_weight(), last / 3.0, 2e-6);
        assert_eq!(sample.path().surface_events(), 2);
    });
}

#[test]
fn g0_reflective_emission_survives_null_proposals_and_uses_original_geometric_sidedness() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(1, 0.4, 0.0, 10.0)]).unwrap();
        let host = glass(1.0, 0.1);
        let media = [SceneMedium::Homogeneous(&host)];
        let metals = [coated_metal()];
        for interaction in [
            Interaction::RoughConductor {
                optics: 0,
                medium: 0,
                surface: ConductorSurface::try_rough(1.0).unwrap(),
            },
            Interaction::BrushedConductor {
                optics: 0,
                medium: 0,
                detail: SurfaceDetail::try_new(0.5, 2.0, BrushingPattern::Radial, 0.2, None)
                    .unwrap(),
            },
        ] {
            let native = CurvedPolarizedScene::try_new_with_conductors(
                cx,
                &geometry,
                &media,
                &metals,
                vec![SurfaceBinding {
                    object_id: 1,
                    interaction,
                }],
                100_000,
            )
            .unwrap();
            for (sides, z, direction, expected) in [
                (EmissionSides::Both, 0.0, 1.0, 2.0 * (-0.04_f64).exp()),
                (EmissionSides::Front, 0.0, 1.0, 0.0),
                (EmissionSides::Front, 1.0, -1.0, 2.0 * (-0.06_f64).exp()),
            ] {
                let sources = [emitter(1, 2.0, sides)];
                let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
                let input = Ray {
                    origin: Point3::new(0.4, 0.0, z),
                    dir: Vec3::new(0.0, 0.0, direction),
                };
                let mut calls = 0;
                let sample = camera
                    .sample_ray_with_samples(cx, input, 550.0, 0, settings(), |_| {
                        calls += 1;
                        [0.999_999, 0.3, 0.0]
                    })
                    .unwrap();
                assert_eq!(calls, 1);
                assert_eq!(
                    sample.endpoint(),
                    RadianceEndpoint::RejectedMicrofacet { object_id: 1 }
                );
                assert_eq!(sample.transport_weight(), 0.0);
                near(sample.radiance_w_m2_sr_nm(), expected, 2e-6);
            }
        }
    });
}

#[test]
fn g1_reflective_emission_reaches_spectral_raster_and_all_molecular_camera_routes() {
    with_cx(|cx, _| {
        let geometry =
            InstanceScene::try_new(vec![plane(1, 0.4, 0.0, 10.0), plane(9, -0.6, 0.0, 10.0)])
                .unwrap();
        let host = glass(1.33, 0.0);
        let media = [SceneMedium::Homogeneous(&host)];
        let metals = [coated_metal()];
        let native = CurvedPolarizedScene::try_new_with_conductors(
            cx,
            &geometry,
            &media,
            &metals,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::SmoothConductor {
                        optics: 0,
                        medium: 0,
                    },
                },
                opaque(9),
            ],
            100_000,
        )
        .unwrap();
        let sources = [
            emitter(1, 3.0, EmissionSides::Both),
            emitter(9, 5.0, EmissionSides::Both),
        ];
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let physical = PhysicalCamera::try_look_at(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 0.4),
            Vec3::new(0.0, 1.0, 0.0),
            CameraProjection::try_half_tangent(0.15).unwrap(),
            0.4,
            Aperture::try_circular(0.0).unwrap(),
        )
        .unwrap();
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        let pixel = raster
            .sample_pixel(cx, &camera, 0, settings(), 17, 0, 4)
            .unwrap();
        assert_eq!(
            pixel,
            raster
                .sample_pixel(cx, &camera, 0, settings(), 17, 0, 4)
                .unwrap()
        );
        for (sample, &wavelength) in pixel.spectral_samples().iter().zip(pixel.wavelengths_nm()) {
            let fresnel = metals[0]
                .fresnel(wavelength, 1.33, pixel.camera_ray().dir.z)
                .unwrap();
            near(sample.radiance_w_m2_sr_nm(), 3.0 + 5.0 * fresnel, 2e-6);
        }
        #[cfg(feature = "volumes")]
        {
            use crate::polarization::dielectric::scene::curved::camera::participating::CurvedRayleighCameraScene;
            use crate::polarization::dielectric::scene::curved::participating::{
                CurvedRayleighMedium, CurvedRayleighScene, CurvedRayleighTraceSettings,
            };
            use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};
            let limits = CurvedRayleighTraceSettings {
                surfaces: settings(),
                max_volume_events: 8,
            };
            let laws = [None];
            let scattering = CurvedRayleighScene::try_new(cx, &native, &laws).unwrap();
            let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &scattering).unwrap();
            let molecular_pixel = raster
                .sample_pixel_rayleigh(cx, &molecular, 0, limits, 17, 0, 4)
                .unwrap();
            assert_eq!(pixel, *molecular_pixel.pixel_sample());
            assert_eq!(
                camera
                    .sample_ray(cx, ray(0.0, 1.0), 550.0, 0, settings(), 17, 4)
                    .unwrap(),
                *molecular
                    .sample_ray(cx, ray(0.0, 1.0), 550.0, 0, limits, 17, 4)
                    .unwrap()
                    .radiance_sample()
            );

            // One real Rayleigh collision AFTER the emitting mirror. At normal
            // incidence its beam is unpolarized; M00/pdf=1, so the independent
            // scalar reference is still Le_mirror + F*Le_terminal. No duplicate
            // survival factor belongs in this conditional path contribution.
            let spectrum = MassScatteringSpectrum::constant(1.0).unwrap();
            let phase = PhaseFunction::Rayleigh;
            let laws = [Some(
                CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
            )];
            let scattering = CurvedRayleighScene::try_new(cx, &native, &laws).unwrap();
            let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &scattering).unwrap();
            let flights = |ordinal| match ordinal {
                0 => [0.9, 0.9, 0.3],
                1 => [1.0 - (-0.2_f64).exp(), 0.9, 0.3],
                2 => [1.0 - 1e-12, 0.9, 0.3],
                _ => panic!("finite terminal source precedes the next collision"),
            };
            let sample = molecular
                .sample_ray_with_samples(
                    cx,
                    ray(0.0, 1.0),
                    550.0,
                    0,
                    limits,
                    |_| panic!("delta mirror"),
                    flights,
                )
                .unwrap();
            assert_eq!(sample.volume_events(), 1);
            near(
                sample.radiance_sample().radiance_w_m2_sr_nm(),
                3.0 + 5.0 * metals[0].fresnel(550.0, 1.33, 1.0).unwrap(),
                3e-6,
            );
            let tracked = molecular
                .sample_ray_with_tracking_samples(
                    cx,
                    ray(0.0, 1.0),
                    550.0,
                    0,
                    limits,
                    |_| panic!("delta mirror"),
                    |ordinal| {
                        let [a, b, c] = flights(ordinal);
                        [a, b, c, 0.25]
                    },
                )
                .unwrap();
            assert_eq!(sample, tracked);
        }
    });
}

#[test]
fn g0_reflective_emission_overflow_and_late_cancellation_return_no_partial_sample() {
    with_cx(|cx, gate| {
        let geometry =
            InstanceScene::try_new(vec![plane(1, 0.4, 0.0, 10.0), plane(9, -0.6, 0.0, 10.0)])
                .unwrap();
        let host = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&host)];
        let reflectances = [SpectralReflectance::constant(1.0).unwrap()];
        let native = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![diffuse(), opaque(9)],
            100_000,
        )
        .unwrap();
        let sources = [
            emitter(1, f64::MAX * 0.75, EmissionSides::Both),
            emitter(9, f64::MAX * 0.75, EmissionSides::Both),
        ];
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        assert!(matches!(
            camera.sample_ray_with_samples(cx, ray(0.0, 1.0), 550.0, 0, settings(), |_| [
                0.0, 0.25, 0.5
            ]),
            Err(CameraRadianceError::NumericalRange)
        ));
        let sources = [
            emitter(1, 2.0, EmissionSides::Both),
            emitter(9, 5.0, EmissionSides::Both),
        ];
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        assert!(matches!(
            camera.sample_ray_with_samples(cx, ray(0.0, 1.0), 550.0, 0, settings(), |_| {
                gate.request();
                [0.0, 0.25, 0.5]
            }),
            Err(CameraRadianceError::Scene(CurvedSceneError::Scene(
                SceneError::Cancelled
            )))
        ));
    });
}
