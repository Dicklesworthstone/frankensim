use super::*;
use crate::charts::TriMesh;
use crate::conductor::{
    CoherentDielectricCoating, ConductorDataStatus, ConductorIorSample, ConductorOptics,
    ConductorSource, ConductorSurface,
};
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricEvent, DielectricGlass, DielectricSurface,
    GlassProvenance, sample_rough_dielectric,
};
use crate::gradient_index::GradientIndexOptions;
use crate::gradient_index::event::GradientEventOptions;
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::polarization::dielectric::scene::SurfaceBinding;
use crate::polarization::dielectric::scene::curved::SceneMedium;
use crate::surface_detail::{BrushingPattern, RadialMark, SurfaceDetail, shading_frame};
use fs_blake3::{ContentHash, hash_domain};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::Point3;

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 31,
                kernel_id: 7,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}
fn glass(n: f64, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(n).unwrap(),
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}
fn plane(id: u64, z: f64, placement: RigidTransform) -> GeometryInstance {
    GeometryInstance::try_new(
        id,
        hash_domain("org.frankensim.test.rough-camera", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(
            vec![[-20.0, -20.0, z], [20.0, -20.0, z], [0.0, 20.0, z]],
            vec![[0, 1, 2]],
        )),
        placement,
    )
    .unwrap()
}
fn settings() -> CurvedTraceSettings {
    CurvedTraceSettings {
        max_surface_events: 8,
        max_distance_m: 8.0,
        max_total_steps: 16_384,
        max_query_work: 1_000_000,
        event: GradientEventOptions {
            integration: GradientIndexOptions {
                position_tolerance_m: 1e-11,
                direction_tolerance: 1e-12,
                ..GradientIndexOptions::default()
            },
            surface_tolerance_m: 1e-7,
            ..GradientEventOptions::default()
        },
    }
}
fn sources(ids: &[u64]) -> Vec<SurfaceEmitter> {
    ids.iter()
        .map(|&id| {
            SurfaceEmitter::new(
                id,
                RadianceSpectrum::try_constant(4.0).unwrap(),
                EmissionSides::Both,
            )
        })
        .collect()
}
fn opaque(id: u64) -> SurfaceBinding {
    SurfaceBinding {
        object_id: id,
        interaction: Interaction::Opaque,
    }
}
fn unit(v: Vec3) -> Vec3 {
    v.scale(1.0 / v.norm())
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}, tolerance {tolerance}"
    );
}
fn metal() -> ConductorOptics {
    let source =
        ConductorSource::try_new(ContentHash([91; 32]), ConductorDataStatus::Representative)
            .unwrap();
    ConductorOptics::try_new(
        core::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, 0.2, 3.0).unwrap()
        }),
        source,
    )
    .unwrap()
}

#[test]
fn g1_oblique_rough_dielectric_camera_matches_radiance_bsdf_in_both_media_directions() {
    with_cx(|cx| {
        let identity = RigidTransform::identity();
        let geometry = InstanceScene::try_new(vec![
            plane(1, 0.0, identity),
            plane(9, 1.0, identity),
            plane(10, -1.0, identity),
        ])
        .unwrap();
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let cosine = 0.9_f64;
        let wo = Vec3::new((1.0 - cosine * cosine).sqrt(), 0.0, cosine);
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 0.5),
            dir: wo.scale(-1.0),
        };
        let emitters = sources(&[9, 10]);
        for (eta_i, eta_t) in [(1.0, 1.5), (1.5, 1.0)] {
            let incident = glass(eta_i, 0.1);
            let transmitted = glass(eta_t, 0.3);
            let media = [
                SceneMedium::Homogeneous(&incident),
                SceneMedium::Homogeneous(&transmitted),
            ];
            for surface in [
                DielectricSurface::POLISHED,
                DielectricSurface::try_rough(0.25).unwrap(),
            ] {
                let native = CurvedPolarizedScene::try_new(
                    cx,
                    &geometry,
                    &media,
                    vec![
                        SurfaceBinding {
                            object_id: 1,
                            interaction: Interaction::Dielectric {
                                outside: 0,
                                inside: 1,
                                surface,
                            },
                        },
                        opaque(9),
                        opaque(10),
                    ],
                    100_000,
                )
                .unwrap();
                let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
                for (event_draw, event) in [
                    (0.0, DielectricEvent::Reflection),
                    (0.999, DielectricEvent::Transmission),
                ] {
                    let draw = [0.17, 0.63, event_draw];
                    let scalar = sample_rough_dielectric(
                        normal,
                        wo,
                        eta_i,
                        eta_t,
                        surface.roughness_alpha().unwrap(),
                        draw[0],
                        draw[1],
                        draw[2],
                    )
                    .unwrap()
                    .unwrap();
                    assert_eq!(scalar.event, event);
                    let mut draws = 0;
                    let sample = camera
                        .sample_ray_with_samples(cx, ray, 550.0, 0, settings(), |_| {
                            draws += 1;
                            draw
                        })
                        .unwrap();
                    assert_eq!(draws, 1);
                    let transmitted = event == DielectricEvent::Transmission;
                    let distance_in = 0.5 / cosine;
                    let distance_out = 1.0 / scalar.direction.z.abs();
                    let absorption_out = if transmitted { 0.3 } else { 0.1 };
                    let beer = (-0.1 * distance_in - absorption_out * distance_out).exp();
                    // This scalar radiance-mode oracle ALREADY contains the
                    // projected cosine, VNDF PDF, branch proposal and eta².
                    close(
                        sample.radiance_w_m2_sr_nm(),
                        4.0 * scalar.radiance_weight * beer,
                        3e-6,
                    );
                    close(sample.path().traveled_m(), distance_in + distance_out, 2e-6);
                    assert_eq!(
                        sample.endpoint(),
                        RadianceEndpoint::Surface {
                            object_id: if transmitted { 10 } else { 9 }
                        }
                    );
                }
            }
        }
    });
}

#[test]
fn g1_rough_metal_camera_has_one_fresnel_and_one_smith_mask_with_physical_beer() {
    with_cx(|cx| {
        let identity = RigidTransform::identity();
        let geometry =
            InstanceScene::try_new(vec![plane(1, 0.0, identity), plane(9, 1.0, identity)]).unwrap();
        let material = glass(1.33, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let emitters = sources(&[9]);
        let bare = metal();
        let coated = bare.with_coating(
            CoherentDielectricCoating::try_new(
                CauchyIor::try_constant(1.45).unwrap(),
                125e-9,
                bare.source(),
            )
            .unwrap(),
        );
        let wo = Vec3::new(0.6, 0.0, 0.8);
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 0.5),
            dir: wo.scale(-1.0),
        };
        let alpha = 0.2_f64;
        for law in [bare, coated] {
            let optics = [law];
            let native = CurvedPolarizedScene::try_new_with_conductors(
                cx,
                &geometry,
                &media,
                &optics,
                vec![
                    SurfaceBinding {
                        object_id: 1,
                        interaction: Interaction::RoughConductor {
                            optics: 0,
                            medium: 0,
                            surface: ConductorSurface::try_rough(alpha).unwrap(),
                        },
                    },
                    opaque(9),
                ],
                100_000,
            )
            .unwrap();
            let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
            for wavelength in [450.0, 650.0] {
                let sample = camera
                    .sample_ray_with_samples(cx, ray, wavelength, 0, settings(), |_| {
                        [0.17, 0.63, 0.5]
                    })
                    .unwrap();
                assert_eq!(
                    sample.endpoint(),
                    RadianceEndpoint::Surface { object_id: 9 }
                );
                let wi = sample.path().beam().ray().dir;
                let facet = unit(Vec3::new(wo.x + wi.x, wo.y + wi.y, wo.z + wi.z));
                let fresnel = law.fresnel(wavelength, 1.33, wo.dot(facet)).unwrap();
                // Independently cancel D in the correlated GGX f*cos / PDF.
                // The current reflection proposal is bounded by a spherical
                // cap, so the old unbounded outgoing G1 is not its weight.
                let lambda = |c: f64| {
                    0.5 * ((1.0 + alpha * alpha * (1.0 - c * c) / (c * c)).sqrt() - 1.0)
                };
                let s2 = (1.0 + (1.0 - wo.z * wo.z).sqrt()).powi(2);
                let k = (1.0 - alpha * alpha) * s2 / (s2 + alpha * alpha * wo.z * wo.z);
                let stretched = (alpha * alpha + (1.0 - alpha * alpha) * wo.z * wo.z).sqrt();
                let masking = (k * wo.z + stretched)
                    / (2.0 * wo.z * (1.0 + lambda(wo.z) + lambda(wi.z)));
                let distance = 0.5 / wo.z + 1.0 / wi.z;
                close(
                    sample.radiance_w_m2_sr_nm(),
                    4.0 * fresnel * masking * (-0.2 * distance).exp(),
                    3e-6,
                );
                close(sample.path().traveled_m(), distance, 2e-6);
            }
        }
    });
}

#[test]
fn g0_native_rejected_facets_are_exact_zero_camera_samples_without_resampling() {
    with_cx(|cx| {
        let geometry =
            InstanceScene::try_new(vec![plane(1, 0.0, RigidTransform::identity())]).unwrap();
        let outside = glass(1.0, 0.2);
        let inside = glass(1.5, 0.0);
        let media = [
            SceneMedium::Homogeneous(&outside),
            SceneMedium::Homogeneous(&inside),
        ];
        let optics = [metal()];
        for interaction in [
            Interaction::RoughConductor {
                optics: 0,
                medium: 0,
                surface: ConductorSurface::try_rough(1.0).unwrap(),
            },
            Interaction::Dielectric {
                outside: 0,
                inside: 1,
                surface: DielectricSurface::try_rough(1.0).unwrap(),
            },
        ] {
            let native = CurvedPolarizedScene::try_new_with_conductors(
                cx,
                &geometry,
                &media,
                &optics,
                vec![SurfaceBinding {
                    object_id: 1,
                    interaction,
                }],
                100_000,
            )
            .unwrap();
            let camera = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
            let ray = Ray {
                origin: Point3::new(0.0, 0.0, 1.0),
                dir: Vec3::new(0.0, 0.0, -1.0),
            };
            let mut draws = 0;
            let sample = camera
                .sample_ray_with_samples(cx, ray, 550.0, 0, settings(), |_| {
                    draws += 1;
                    [0.9, 0.2, 0.0]
                })
                .unwrap();
            assert_eq!(draws, 1);
            assert_eq!(sample.radiance_w_m2_sr_nm(), 0.0);
            assert_eq!(sample.transport_weight(), 0.0);
            assert_eq!(
                sample.endpoint(),
                RadianceEndpoint::RejectedMicrofacet { object_id: 1 }
            );
            assert_eq!(
                sample.path().termination(),
                PathTermination::RejectedMicrofacet
            );
            assert_eq!(sample.path().surface_events(), 1);
            assert!(
                sample.path().beam().power().i > 0.0,
                "retained state is before the rejected proposal"
            );
            assert!(
                camera
                    .sample_ray_with_samples(cx, ray, 550.0, 0, settings(), |_| [
                        f64::NAN,
                        0.2,
                        0.0
                    ])
                    .is_err()
            );
        }
    });
}

#[test]
fn g3_untilted_brushing_moves_with_its_body_and_tilted_shading_is_admitted() {
    with_cx(|cx| {
        let detail = SurfaceDetail::try_new(
            0.1,
            4.0,
            BrushingPattern::Radial,
            0.0,
            Some(RadialMark {
                azimuth_rad: 0.0,
                half_width_rad: 0.1,
            }),
        )
        .unwrap();
        assert_eq!(detail.normal_perturbation_rad(), 0.0);
        let material = glass(1.33, 0.2);
        let media = [SceneMedium::Homogeneous(&material)];
        let optics = [metal()];
        let emitters = sources(&[9]);
        let angle = 0.7_f64;
        let moved = RigidTransform::try_new(
            [0.0, (angle / 2.0).sin(), 0.0, (angle / 2.0).cos()],
            [2.0, -3.0, 4.0],
        )
        .unwrap();
        let normal = Vec3::new(0.0, 0.0, 1.0);
        let local_frame = shading_frame(&detail, Point3::new(0.4, 0.0, 0.0), normal).unwrap();
        assert!(local_frame.mark_covered);
        let (ax, ay) = detail.alpha_pair_for_mark(true);
        // The bounded spherical-cap VNDF proposal (6a13950eae) leaves its mass
        // 1 - (1 - k) G1(wo) / 2 in the weight, with k = (1 - a^2) s^2 / (s^2 + a^2 c^2),
        // a = min(ax, ay, 1), c = cos_o and s = 1 + sin_o (tracer.rs ggx_reflection_cap_scale and
        // ggx_reflection_cap_mass). The camera looks straight down the untilted normal, so
        // c = 1, s = 1, G1(wo) = 1 and the mass is 1 / (1 + a^2).
        let a = ax.min(ay).min(1.0);
        let k = (1.0 - a * a) / (1.0 + a * a);
        let cap_mass = 1.0 - 0.5 * (1.0 - k);
        let mut baseline = None;
        for placement in [RigidTransform::identity(), moved] {
            let geometry =
                InstanceScene::try_new(vec![plane(1, 0.0, placement), plane(9, 1.0, placement)])
                    .unwrap();
            let native = CurvedPolarizedScene::try_new_with_conductors(
                cx,
                &geometry,
                &media,
                &optics,
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
            let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
            let ray = Ray {
                origin: placement.transform_point(Point3::new(0.4, 0.0, 0.5)),
                dir: placement.transform_vector(normal.scale(-1.0)),
            };
            let sample = camera
                .sample_ray_with_samples(cx, ray, 550.0, 0, settings(), |_| [0.13, 0.27, 0.5])
                .unwrap();
            assert_eq!(
                sample.endpoint(),
                RadianceEndpoint::Surface { object_id: 9 }
            );
            let wi = placement
                .inverse()
                .unwrap()
                .transform_vector(sample.path().beam().ray().dir);
            let x = wi.dot(local_frame.tangent);
            let y = wi.dot(local_frame.bitangent);
            let z = wi.dot(local_frame.shading_normal);
            let masking =
                2.0 / (1.0 + (1.0 + (ax * ax * x * x + ay * ay * y * y) / (z * z)).sqrt());
            let facet = unit(Vec3::new(wi.x, wi.y, 1.0 + wi.z));
            let fresnel = optics[0].fresnel(550.0, 1.33, facet.z).unwrap();
            let distance = 0.5 + 1.0 / wi.z;
            close(
                sample.radiance_w_m2_sr_nm(),
                4.0 * fresnel * cap_mass * masking * (-0.2 * distance).exp(),
                3e-6,
            );
            if let Some(value) = baseline {
                close(sample.radiance_w_m2_sr_nm(), value, 3e-6);
            } else {
                baseline = Some(sample.radiance_w_m2_sr_nm());
            }
        }
        let tilted = SurfaceDetail::try_new(0.1, 4.0, BrushingPattern::Radial, 0.1, None).unwrap();
        assert_eq!(tilted.normal_perturbation_rad(), 0.1);
        let geometry =
            InstanceScene::try_new(vec![plane(1, 0.0, RigidTransform::identity())]).unwrap();
        let native = CurvedPolarizedScene::try_new_with_conductors(
            cx,
            &geometry,
            &media,
            &optics,
            vec![SurfaceBinding {
                object_id: 1,
                interaction: Interaction::BrushedConductor {
                    optics: 0,
                    medium: 0,
                    detail: tilted,
                },
            }],
            100_000,
        )
        .unwrap();
        assert!(CurvedCameraScene::try_new(cx, &native, &[], None).is_ok());
    });
}
