use super::*;
use crate::camera::{Aperture, CameraProjection, PhysicalCamera};
use crate::conductor::{
    CauchyCoatingLayer, CoherentDielectricCoating, CoherentPeriodicCoating, ConductorDataStatus,
    ConductorOptics, ConductorSource, ConductorSurface,
};
use crate::dielectric::{DielectricCoating, DielectricEvent};
use crate::polarization::dielectric::scene::curved::camera::raster::CurvedCameraRaster;
use crate::surface_detail::{BrushingPattern, SurfaceDetail, shading_frame};

fn source() -> ConductorSource {
    ConductorSource::try_new(ContentHash([86; 32]), ConductorDataStatus::Representative).unwrap()
}

fn metal() -> ConductorOptics {
    ConductorOptics::representative_tungsten().with_coating(
        CoherentDielectricCoating::try_new(
            CauchyIor::try_constant(1.45).unwrap(),
            125e-9,
            source(),
        )
        .unwrap(),
    )
}

fn conductor_binding(surface: ConductorSurface) -> SurfaceBinding {
    SurfaceBinding {
        object_id: 1,
        interaction: Interaction::RoughConductor {
            optics: 0,
            medium: 0,
            surface,
        },
    }
}

fn smith_lambda(alpha: f64, cosine: f64) -> f64 {
    0.5 * ((1.0 + alpha * alpha * (1.0 - cosine * cosine) / (cosine * cosine)).sqrt() - 1.0)
}

#[test]
fn g1_rough_camera_coated_metal_matches_fresnel_masking_and_physical_beer_path() {
    with_cx(|cx, _| {
        let geometry =
            InstanceScene::try_new(vec![plane(1, 0.4, 0.0, 100.0), plane(9, -0.6, 0.0, 100.0)])
                .unwrap();
        let metals = [metal()];
        let emitters = [emitter(9, 3.0, EmissionSides::Both)];
        for index in [1.0, 1.33] {
            let host = glass(index, 0.12);
            let media = [SceneMedium::Homogeneous(&host)];
            for alpha in [0.05, 0.3] {
                let scene = CurvedPolarizedScene::try_new_with_conductors(
                    cx,
                    &geometry,
                    &media,
                    &metals,
                    vec![
                        conductor_binding(ConductorSurface::try_rough(alpha).unwrap()),
                        opaque(9),
                    ],
                    100_000,
                )
                .unwrap();
                let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
                for cosine in [0.6_f64, 1.0] {
                    let ray = Ray {
                        origin: Point3::new(0.0, 0.0, 0.0),
                        dir: Vec3::new((1.0 - cosine * cosine).sqrt(), 0.0, cosine),
                    };
                    for draw in [[0.0, 0.0, 0.8], [0.13, 0.27, 0.8]] {
                        let mut limits = settings();
                        limits.max_distance_m = 16.0;
                        let mut calls = 0;
                        let result = camera
                            .sample_ray_with_samples(cx, ray, 550.0, 0, limits, |event| {
                                assert_eq!(event, 0);
                                calls += 1;
                                draw
                            })
                            .unwrap();
                        assert_eq!(calls, 1);
                        assert_eq!(
                            result.endpoint(),
                            RadianceEndpoint::Surface { object_id: 9 }
                        );
                        let wo = ray.dir.scale(-1.0);
                        let wi = result.path().beam().ray().dir;
                        let half = Vec3::new(wo.x + wi.x, wo.y + wi.y, wo.z + wi.z);
                        let half = half.scale(1.0 / half.norm());
                        let fresnel = metals[0].fresnel(550.0, index, wo.dot(half)).unwrap();
                        // The unchanged VNDF proposal cancels to F*G2/G1(wo).
                        let lambda_o = smith_lambda(alpha, cosine);
                        let lambda_i = smith_lambda(alpha, -wi.z);
                        let expected = fresnel
                            * (1.0 + lambda_o) / (1.0 + lambda_o + lambda_i)
                            * (-0.12 * result.path().traveled_m()).exp();
                        near(result.transport_weight(), expected, 3e-6);
                        near(result.radiance_w_m2_sr_nm(), 3.0 * expected, 9e-6);
                        assert_eq!(result.path().medium_index(), 0);
                        assert_eq!(result.path().surface_events(), 2);
                    }
                }
            }
        }
    });
}

#[test]
fn g1_rough_camera_glass_uses_native_radiance_weights_in_both_directions_and_branches() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![
            plane(1, 0.4, 0.0, 100.0),
            plane(8, -0.6, 0.0, 100.0),
            plane(9, 1.4, 0.0, 100.0),
        ])
        .unwrap();
        let inside = glass(1.0, 0.0);
        let outside = glass(1.5, 0.0);
        let media = [
            SceneMedium::Homogeneous(&inside),
            SceneMedium::Homogeneous(&outside),
        ];
        let emitters = [
            emitter(8, 3.0, EmissionSides::Both),
            emitter(9, 3.0, EmissionSides::Both),
        ];
        let coating = DielectricCoating::try_new(
            CoherentPeriodicCoating::try_new(
                &[
                    CauchyCoatingLayer::try_new(CauchyIor::try_constant(1.7).unwrap(), 0.05, 73e-9)
                        .unwrap(),
                ],
                1,
                source(),
            )
            .unwrap(),
        )
        .unwrap();
        let bare = DielectricSurface::try_rough(0.2).unwrap();
        for surface in [bare, bare.with_coating(coating)] {
            let scene = CurvedPolarizedScene::try_new(
                cx,
                &geometry,
                &media,
                vec![
                    SurfaceBinding {
                        object_id: 1,
                        interaction: Interaction::Dielectric {
                            outside: 1,
                            inside: 0,
                            surface,
                        },
                    },
                    opaque(8),
                    opaque(9),
                ],
                100_000,
            )
            .unwrap();
            let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
            for (origin, direction, initial, eta_i, eta_t) in
                [(0.0, 1.0, 0, 1.0, 1.5), (0.8, -1.0, 1, 1.5, 1.0)]
            {
                let input = ray(origin, direction);
                let normal = input.dir.scale(-1.0);
                for branch in [0.0, 0.999] {
                    let draw = [0.13, 0.27, branch];
                    let reference = surface
                        .interface(550.0, direction > 0.0)
                        .unwrap()
                        .sample_rough(normal, normal, eta_i, eta_t, 0.2, draw[0], draw[1], draw[2])
                        .unwrap()
                        .unwrap();
                    let result = camera
                        .sample_ray_with_samples(cx, input, 550.0, initial, settings(), |_| draw)
                        .unwrap();
                    let transmission = reference.event == DielectricEvent::Transmission;
                    assert_eq!(transmission, branch > 0.5);
                    let target = if transmission { 1 - initial } else { initial };
                    let emitter_id = if reference.direction.z > 0.0 { 9 } else { 8 };
                    assert_eq!(
                        result.endpoint(),
                        RadianceEndpoint::Surface {
                            object_id: emitter_id
                        }
                    );
                    assert_eq!(result.path().medium_index(), target);
                    near(result.transport_weight(), reference.radiance_weight, 3e-6);
                    near(
                        result.radiance_w_m2_sr_nm(),
                        3.0 * reference.radiance_weight,
                        9e-6,
                    );
                    near(
                        result.path().beam().ray().dir.dot(reference.direction),
                        1.0,
                        2e-10,
                    );
                    let power_factor = if transmission {
                        (eta_t / eta_i) * (eta_t / eta_i)
                    } else {
                        1.0
                    };
                    near(
                        result.path().beam().power().i,
                        reference.radiance_weight * power_factor,
                        3e-6,
                    );
                }
            }
        }
    });
}

struct VerticalIndex;
impl GradientIndexField for VerticalIndex {
    fn sample(&self, point: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
        let n = 1.2 + 0.02 * point.z;
        GradientIndexSample::try_new(n, Vec3::new(0.0, 0.0, 0.02), 0.1, Some(n))
    }
}

#[test]
fn g1_rough_camera_grin_reflection_uses_contact_index_and_one_endpoint_radiance_ratio() {
    with_cx(|cx, _| {
        let geometry =
            InstanceScene::try_new(vec![plane(1, 0.4, 0.0, 10.0), plane(9, -0.6, 0.0, 10.0)])
                .unwrap();
        let field = VerticalIndex;
        let identity = ContentHash([87; 32]);
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity,
        }];
        let metals = [metal()];
        let scene = CurvedPolarizedScene::try_new_with_conductors(
            cx,
            &geometry,
            &media,
            &metals,
            vec![
                conductor_binding(ConductorSurface::try_rough(0.3).unwrap()),
                opaque(9),
            ],
            100_000,
        )
        .unwrap();
        let emitters = [emitter(9, 3.0, EmissionSides::Both)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let result = camera
            .sample_ray_with_samples(cx, ray(0.0, 1.0), 550.0, 0, settings(), |_| [0.0, 0.0, 0.9])
            .unwrap();
        // The gradient is parallel to this chosen facet path: z=0 -> 0.4 -> -0.6.
        let fresnel = metals[0].fresnel(550.0, 1.208, 1.0).unwrap();
        let ratio = 1.2_f64 / 1.188;
        let expected = fresnel * (-0.1_f64 * 1.4).exp() * ratio * ratio;
        near(result.transport_weight(), expected, 3e-6);
        near(result.radiance_w_m2_sr_nm(), 3.0 * expected, 9e-6);
        near(result.path().traveled_m(), 1.4, 2e-6);
        assert_eq!(result.path().medium_snapshot_identity(), Some(identity));
        assert_eq!(
            result.endpoint(),
            RadianceEndpoint::Surface { object_id: 9 }
        );
    });
}

#[test]
fn g0_tilted_brush_camera_null_samples_preserve_both_horizons_and_sample_count() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(1, 0.4, 0.0, 10.0)]).unwrap();
        let host = glass(1.0, 0.0);
        let target = glass(1.5, 0.0);
        let media = [
            SceneMedium::Homogeneous(&host),
            SceneMedium::Homogeneous(&target),
        ];
        let metals = [metal()];
        for pattern in [BrushingPattern::Radial, BrushingPattern::Circular] {
            let binding = SurfaceBinding {
                object_id: 1,
                interaction: Interaction::BrushedConductor {
                    optics: 0,
                    medium: 0,
                    detail: SurfaceDetail::try_new(0.5, 2.0, pattern, 0.2, None).unwrap(),
                },
            };
            let scene = CurvedPolarizedScene::try_new_with_conductors(
                cx,
                &geometry,
                &media,
                &metals,
                vec![binding],
                100_000,
            )
            .unwrap();
            let camera = CurvedCameraScene::try_new(cx, &scene, &[], None).unwrap();
            let input = Ray {
                origin: Point3::new(0.4, 0.0, 0.0),
                dir: Vec3::new(0.0, 0.0, 1.0),
            };
            let mut calls = 0;
            let result = camera
                .sample_ray_with_samples(cx, input, 550.0, 0, settings(), |_| {
                    calls += 1;
                    [0.999_999, 0.3, 0.0]
                })
                .unwrap();
            assert_eq!(calls, 1);
            assert_eq!(result.radiance_w_m2_sr_nm(), 0.0);
            assert_eq!(result.transport_weight(), 0.0);
            assert_eq!(
                result.endpoint(),
                RadianceEndpoint::RejectedMicrofacet { object_id: 1 }
            );
            assert_eq!(
                result.path().termination(),
                PathTermination::RejectedMicrofacet
            );
            near(result.path().beam().power().i, 1.0, 1e-14);
            assert!(
                camera
                    .sample_ray_with_samples(cx, input, 550.0, 0, settings(), |_| [
                        f64::NAN,
                        0.3,
                        0.0
                    ])
                    .is_err()
            );
        }
    });
}

#[test]
fn g1_tilted_brush_camera_raster_matches_radiance_weights_and_keeps_null_lanes() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(1, 0.4, 0.0, 100.0)]).unwrap();
        let host = glass(1.33, 0.0);
        let media = [SceneMedium::Homogeneous(&host)];
        let metals = [metal()];
        let detail = SurfaceDetail::try_new(0.5, 2.0, BrushingPattern::Radial, 0.2, None).unwrap();
        let native = CurvedPolarizedScene::try_new_with_conductors(
            cx,
            &geometry,
            &media,
            &metals,
            vec![SurfaceBinding {
                object_id: 1,
                interaction: Interaction::BrushedConductor {
                    optics: 0,
                    medium: 0,
                    detail,
                },
            }],
            100_000,
        )
        .unwrap();
        let environment = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap();
        let scene = CurvedCameraScene::try_new(
            cx,
            &native,
            &[],
            Some(EnvironmentRadiance::try_new(&environment, 2.0).unwrap()),
        )
        .unwrap();
        let camera = PhysicalCamera::try_look_at(
            Point3::new(0.4, 0.0, 0.0),
            Point3::new(0.4, 0.0, 0.4),
            Vec3::new(0.0, 1.0, 0.0),
            CameraProjection::try_half_tangent(0.15).unwrap(),
            0.4,
            Aperture::try_circular(0.0).unwrap(),
        )
        .unwrap();
        let raster = CurvedCameraRaster::try_new(&camera, 1, 1).unwrap();
        let mut rejected = 0;
        let mut positive = 0;
        for sample in 0..16 {
            let result = raster
                .sample_pixel(cx, &scene, 0, settings(), 81, 0, sample)
                .unwrap();
            assert_eq!(
                result,
                raster
                    .sample_pixel(cx, &scene, 0, settings(), 81, 0, sample)
                    .unwrap()
            );
            assert!(result.xyz().iter().all(|x| x.is_finite() && *x >= 0.0));
            for (lane, wavelength) in result
                .spectral_samples()
                .iter()
                .zip(result.wavelengths_nm())
            {
                if lane.path().termination() == PathTermination::RejectedMicrofacet {
                    rejected += 1;
                    assert_eq!(lane.radiance_w_m2_sr_nm(), 0.0);
                } else {
                    positive += 1;
                    assert_eq!(lane.endpoint(), RadianceEndpoint::Environment);
                    assert!(lane.radiance_w_m2_sr_nm() > 0.0);
                    let contact = lane.path().last_contact().unwrap();
                    let frame = shading_frame(
                        &detail,
                        contact.witness.surface_point,
                        Vec3::new(0.0, 0.0, -1.0),
                    )
                    .unwrap();
                    let wo = result.camera_ray().dir.scale(-1.0);
                    let wi = lane.path().beam().ray().dir;
                    let half = Vec3::new(wo.x + wi.x, wo.y + wi.y, wo.z + wi.z);
                    let half = half.scale(1.0 / half.norm());
                    let (ax, ay) = detail.alpha_pair();
                    // The weight is Fresnel times two closed-form factors, both evaluated
                    // in the shading frame:
                    // - 4d5123920d: height-correlated Smith visibility G2(wo, wi) / G1(wo)
                    //   = (1 + lambda_o) / (1 + lambda_o + lambda_i). It replaced the
                    //   separable G1(wi) = 1 / (1 + lambda_i).
                    // - 6a13950eae: the bounded spherical-cap VNDF proposal. Each weight
                    //   carries the cap mass 1 - (1 - k) * G1(wo) / 2, with
                    //   k = (1 - a^2) s^2 / (s^2 + a^2 c^2), a = min(ax, ay, 1),
                    //   c = cos_o and s = 1 + sin_o (tracer.rs ggx_reflection_cap_scale
                    //   and ggx_reflection_cap_mass).
                    // R-c measured observed/expected = 0.826579323586119186 on sample 0;
                    // cap 0.826189551574754710 times G2/G1 over G1(wi) 1.00047177068582127
                    // matches it to 2 ulp.
                    let lambda = |v: [f64; 3]| {
                        0.5 * ((1.0
                            + (ax * ax * v[0] * v[0] + ay * ay * v[1] * v[1]) / (v[2] * v[2]))
                            .sqrt()
                            - 1.0)
                    };
                    let wo_local = [
                        frame.tangent.dot(wo),
                        frame.bitangent.dot(wo),
                        frame.shading_normal.dot(wo),
                    ];
                    let wi_local = [
                        frame.tangent.dot(wi),
                        frame.bitangent.dot(wi),
                        frame.shading_normal.dot(wi),
                    ];
                    let (lambda_o, lambda_i) = (lambda(wo_local), lambda(wi_local));
                    let correlated_visibility = (1.0 + lambda_o) / (1.0 + lambda_o + lambda_i);
                    let a = ax.min(ay).min(1.0);
                    let c = wo_local[2].clamp(0.0, 1.0);
                    let s = 1.0 + (1.0 - c * c).max(0.0).sqrt();
                    let k = (1.0 - a * a) * s * s / (s * s + a * a * c * c);
                    let cap_mass = 1.0 - 0.5 * (1.0 - k) / (1.0 + lambda_o);
                    let expected = metals[0].fresnel(*wavelength, 1.33, wo.dot(half)).unwrap()
                        * correlated_visibility
                        * cap_mass;
                    near(lane.transport_weight(), expected, 3e-6);
                }
            }
        }
        assert!(rejected > 0 && positive > 0);
        assert_eq!(rejected + positive, 16 * 4);
    });
}

#[test]
fn g1_tilted_brush_camera_uses_radiance_while_native_and_molecular_power_stay_adjoint() {
    with_cx(|cx, _| {
        let geometry =
            InstanceScene::try_new(vec![plane(1, 0.4, 0.0, 10.0), plane(9, -0.6, 0.0, 10.0)])
                .unwrap();
        let host = glass(1.33, 0.0);
        let media = [SceneMedium::Homogeneous(&host)];
        let metals = [metal()];
        let detail = SurfaceDetail::try_new(0.1, 3.0, BrushingPattern::Radial, 0.25, None).unwrap();
        let native = CurvedPolarizedScene::try_new_with_conductors(
            cx,
            &geometry,
            &media,
            &metals,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::BrushedConductor {
                        optics: 0,
                        medium: 0,
                        detail,
                    },
                },
                opaque(9),
            ],
            100_000,
        )
        .unwrap();
        let emitters = [emitter(9, 3.0, EmissionSides::Both)];
        let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let input = Ray {
            origin: Point3::new(0.4, 0.0, 0.0),
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        let draw = [0.13, 0.27, 0.5];
        let sample = camera
            .sample_ray_with_samples(cx, input, 550.0, 0, settings(), |_| draw)
            .unwrap();
        assert_eq!(
            sample.endpoint(),
            RadianceEndpoint::Surface { object_id: 9 }
        );
        let beam = native
            .launch_beam(
                input.origin,
                PolarizationFrame::try_new(input.dir, Vec3::new(1.0, 0.0, 0.0)).unwrap(),
                550.0,
                0,
                Stokes::unpolarized(1.0),
            )
            .unwrap();
        let power = native
            .trace_with_samples(cx, beam, 0, settings(), |_| draw)
            .unwrap();
        assert_eq!(power.beam().ray(), sample.path().beam().ray());
        let frame = shading_frame(
            &detail,
            Point3::new(0.4, 0.0, 0.4),
            Vec3::new(0.0, 0.0, -1.0),
        )
        .unwrap();
        let wo = input.dir.scale(-1.0);
        let wi = power.beam().ray().dir;
        let correction = frame.shading_normal.dot(wo) / frame.geometric_normal.dot(wo)
            * frame.geometric_normal.dot(wi)
            / frame.shading_normal.dot(wi);
        near(
            power.beam().power().i,
            sample.transport_weight() * correction,
            3e-6,
        );
        assert!((correction - 1.0).abs() > 1e-3);

        #[cfg(feature = "volumes")]
        {
            use crate::polarization::dielectric::scene::curved::camera::participating::CurvedRayleighCameraScene;
            use crate::polarization::dielectric::scene::curved::participating::{
                CurvedRayleighMedium, CurvedRayleighScene, CurvedRayleighTraceSettings,
            };
            use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};
            let spectrum = MassScatteringSpectrum::constant(0.001).unwrap();
            let phase = PhaseFunction::Rayleigh;
            let law = CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap();
            let limits = CurvedRayleighTraceSettings {
                surfaces: settings(),
                max_volume_events: 8,
            };
            for law in [None, Some(law)] {
                let laws = [law];
                let scattering = CurvedRayleighScene::try_new(cx, &native, &laws).unwrap();
                let molecular =
                    CurvedRayleighCameraScene::try_new(cx, &camera, &scattering).unwrap();
                // For the positive law these flights exceed each finite segment;
                // the path weight must not acquire a second survival factor.
                let result = molecular
                    .sample_ray_with_samples(
                        cx,
                        input,
                        550.0,
                        0,
                        limits,
                        |_| draw,
                        |_| [0.5, 0.3, 0.7],
                    )
                    .unwrap();
                assert_eq!(*result.radiance_sample(), sample);
                let tracked = molecular
                    .sample_ray_with_tracking_samples(
                        cx,
                        input,
                        550.0,
                        0,
                        limits,
                        |_| draw,
                        |_| [0.5, 0.3, 0.7, 0.2],
                    )
                    .unwrap();
                assert_eq!(*tracked.radiance_sample(), sample);
                let old_power = scattering
                    .trace_with_samples(cx, beam, 0, limits, |_| draw, |_| [0.5, 0.3, 0.7])
                    .unwrap();
                assert_eq!(*old_power.path(), power);
                if law.is_none() {
                    assert_eq!(
                        camera
                            .sample_ray(cx, input, 550.0, 0, settings(), 51, 9)
                            .unwrap(),
                        *molecular
                            .sample_ray(cx, input, 550.0, 0, limits, 51, 9)
                            .unwrap()
                            .radiance_sample()
                    );
                }
            }
        }
    });
}
