//! Calibrated native point fixtures at continuous non-Lambertian receivers.

use super::*;
use crate::conductor::{
    ConductorDataStatus, ConductorIorSample, ConductorOptics, ConductorSource, ConductorSurface,
};
use crate::dielectric::{DielectricEvent, DielectricSurface};
use crate::polarization::Stokes;
use crate::polarization::dielectric::transport::PolarizedBeam;
use crate::surface_detail::{BrushingPattern, SurfaceDetail};
use core::f64::consts::{FRAC_PI_2, PI};

fn metal() -> ConductorOptics {
    let provenance = ConductorSource::try_new(
        ContentHash([91; 32]),
        ConductorDataStatus::Representative,
    )
    .unwrap();
    ConductorOptics::try_new(
        core::array::from_fn(|i| {
            ConductorIorSample::try_new(380.0 + 50.0 * i as f64, 0.2, 3.0).unwrap()
        }),
        provenance,
    )
    .unwrap()
}

fn oriented_panel(id: u64, center: Point3, normal: Vec3, radius: f64) -> GeometryInstance {
    let normal = unit(normal);
    let reference = if normal.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let tangent = unit(Vec3::new(
        reference.x - normal.x * reference.dot(normal),
        reference.y - normal.y * reference.dot(normal),
        reference.z - normal.z * reference.dot(normal),
    ));
    let bitangent = Vec3::new(
        normal.y * tangent.z - normal.z * tangent.y,
        normal.z * tangent.x - normal.x * tangent.z,
        normal.x * tangent.y - normal.y * tangent.x,
    );
    let point = |u, v| {
        [
            center.x + radius * (u * tangent.x + v * bitangent.x),
            center.y + radius * (u * tangent.y + v * bitangent.y),
            center.z + radius * (u * tangent.z + v * bitangent.z),
        ]
    };
    GeometryInstance::try_new(
        id,
        hash_domain("rough-point-fixture-panel", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(
            vec![point(-1.0, -1.0), point(1.0, -1.0), point(0.0, 2.0)],
            vec![[0, 1, 2]],
        )),
        RigidTransform::identity(),
    )
    .unwrap()
}

#[test]
fn g1_point_fixture_rough_glass_has_both_media_and_radiance_reciprocity() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![
            plane(1, 0.0, 0.0, 2.0),
            plane(8, -2.0, 0.0, 2.0),
            plane(9, 2.0, 0.0, 2.0),
        ])
        .unwrap();
        let indices = [1.0_f64, 1.5];
        let absorption = [0.13_f64, 0.29];
        let materials = [
            glass(indices[0], absorption[0]),
            glass(indices[1], absorption[1]),
        ];
        let media = [
            SceneMedium::Homogeneous(&materials[0]),
            SceneMedium::Homogeneous(&materials[1]),
        ];
        for alpha in [0.12_f64, 0.3] {
            let scene = CurvedPolarizedScene::try_new_with_materials(
                cx, &geometry, &media, &[], &[], &[],
                vec![
                    SurfaceBinding {
                        object_id: 1,
                        interaction: Interaction::Dielectric {
                            outside: 0,
                            inside: 1,
                            surface: DielectricSurface::try_rough(alpha).unwrap(),
                        },
                    },
                    opaque(8), opaque(9),
                ],
                100_000,
            )
            .unwrap();
            let camera = CurvedCameraScene::try_new(cx, &scene, &[], None).unwrap();
            let mut transmission = [0.0; 2];
            for camera_medium in 0..2 {
                let camera_side = if camera_medium == 0 { 1.0 } else { -1.0 };
                let camera_ray = Ray {
                    origin: Point3::new(0.0, 0.0, 0.5 * camera_side),
                    dir: Vec3::new(0.0, 0.0, -camera_side),
                };
                let ordinary = camera.sample_ray_with_samples(
                    cx, camera_ray, 550.0, camera_medium, trace_settings(),
                    |_| [0.0, 0.0, 0.99],
                ).unwrap();
                for source_medium in 0..2 {
                    let source_side = if source_medium == 0 { 1.0 } else { -1.0 };
                    let lamp = source(Point3::new(0.0, 0.0, source_side));
                    let mut binding = branch(&lamp, 1, Vec3::new(0.0, 0.0, -source_side));
                    binding.initial_medium = source_medium;
                    let branches = [binding];
                    let lit = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
                    let sample = lit.sample_ray_with_samples(
                        cx, camera_ray, 550.0, camera_medium, settings(),
                        |_| [0.0, 0.0, 0.99],
                    ).unwrap();
                    // Independent normal-incidence Fresnel and GGX laws.
                    // The transmission law is the radiance BSDF; adding a
                    // second local eta-squared conversion would break reciprocity.
                    let fresnel = ((indices[0] - indices[1])
                        / (indices[0] + indices[1])).powi(2);
                    let response = if camera_medium == source_medium {
                        fresnel / (4.0 * PI * alpha * alpha)
                    } else {
                        (1.0 - fresnel) / (PI * alpha * alpha
                            * (1.0 - indices[source_medium] / indices[camera_medium]).powi(2))
                    };
                    let incoming = intensity(&lamp)
                        * (-absorption[source_medium] - 0.5 * absorption[camera_medium]).exp();
                    let expected = incoming * response;
                    near(sample.fixture_radiance_w_m2_sr_nm(), expected);
                    near(sample.radiance_sample().radiance_w_m2_sr_nm(), expected);
                    assert_eq!(sample.radiance_sample().path(), ordinary.path());
                    assert_eq!(sample.connections(), 1);
                    if camera_medium != source_medium {
                        transmission[camera_medium] = sample.fixture_radiance_w_m2_sr_nm() / incoming;
                    }
                }
            }
            near(
                transmission[0] / indices[0].powi(2),
                transmission[1] / indices[1].powi(2),
            );
        }
    });
}

#[test]
fn g1_point_fixture_rough_and_brushed_metal_use_original_complex_index_and_finish() {
    with_cx(|cx, _| {
        // A nonzero local radius exercises the radial anisotropic brush rather
        // than its documented isotropic fallback at the origin.
        let x = 0.2;
        let camera_ray = Ray {
            origin: Point3::new(x, 0.0, 0.0),
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        let geometry = InstanceScene::try_new(vec![
            plane(1, 1.0, x, 2.0), plane(9, -1.0, x, 2.0),
        ]).unwrap();
        let conductors = [metal()];
        let detail = SurfaceDetail::try_new(0.12, 3.0, BrushingPattern::Radial, 0.0, None).unwrap();
        let laws = [
            (Interaction::RoughConductor {
                optics: 0, medium: 0, surface: ConductorSurface::try_rough(0.2).unwrap(),
            }, (0.2, 0.2)),
            (Interaction::BrushedConductor {
                optics: 0, medium: 0, detail,
            }, detail.alpha_pair()),
        ];
        for n in [1.0_f64, 1.33] {
            let material = glass(n, 0.1);
            let media = [SceneMedium::Homogeneous(&material)];
            for (interaction, (alpha_x, alpha_y)) in laws {
                let scene = CurvedPolarizedScene::try_new_with_materials(
                    cx, &geometry, &media, &conductors, &[], &[],
                    vec![SurfaceBinding { object_id: 1, interaction }, opaque(9)],
                    100_000,
                ).unwrap();
                let camera = CurvedCameraScene::try_new(cx, &scene, &[], None).unwrap();
                let lamp = source(camera_ray.origin);
                let branches = [branch(&lamp, 1, camera_ray.dir)];
                let lit = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
                let ordinary = camera.sample_ray_with_samples(
                    cx, camera_ray, 550.0, 0, trace_settings(), |_| [0.0; 3],
                ).unwrap();
                let sample = lit.sample_ray_with_samples(
                    cx, camera_ray, 550.0, 0, settings(), |_| [0.0; 3],
                ).unwrap();
                let fresnel = ((n - 0.2).powi(2) + 9.0) / ((n + 0.2).powi(2) + 9.0);
                let expected = intensity(&lamp) * (-0.2_f64).exp()
                    * fresnel / (4.0 * PI * alpha_x * alpha_y);
                near(sample.fixture_radiance_w_m2_sr_nm(), expected);
                near(sample.radiance_sample().radiance_w_m2_sr_nm(), expected);
                assert_eq!(sample.radiance_sample().path(), ordinary.path());
                assert_eq!(sample.connections(), 1);
            }
        }
    });
}

#[test]
fn g1_point_fixture_through_glass_uses_bent_bundle_area_and_resolved_angular_intensity() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![
            oriented_panel(2, Point3::new(0.0, 0.0, 0.25), Vec3::new(0.0, 0.0, -1.0), 2.0),
            plane(3, 0.5, 0.0, 2.0),
            plane(5, 1.0, 0.0, 2.0),
            plane(9, -0.5, 0.0, 2.0),
        ]).unwrap();
        let air = glass(1.0, 0.1);
        let glass = glass(1.5, 0.3);
        let media = [SceneMedium::Homogeneous(&air), SceneMedium::Homogeneous(&glass)];
        let conductors = [metal()];
        let alpha = 0.2;
        let slab = |object_id| SurfaceBinding {
            object_id,
            interaction: Interaction::Dielectric {
                outside: 0, inside: 1, surface: DielectricSurface::SMOOTH,
            },
        };
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx, &geometry, &media, &conductors, &[], &[],
            vec![
                slab(2), slab(3),
                SurfaceBinding { object_id: 5, interaction: Interaction::RoughConductor {
                    optics: 0, medium: 0, surface: ConductorSurface::try_rough(alpha).unwrap(),
                }},
                opaque(9),
            ],
            100_000,
        ).unwrap();
        let camera = CurvedCameraScene::try_new(cx, &scene, &[], None).unwrap();
        let lamp = PhotometricSource::try_new(
            RigidTransform::identity(), 100.0, 1.0,
            AngularDistribution::cosine_power(100.0, FRAC_PI_2).unwrap(),
            SourceSpectrum::try_from_bins([1.0; SPECTRAL_BINS]).unwrap(),
        ).unwrap();
        let events = [DielectricEvent::Transmission; 2];
        let mut binding = branch(&lamp, 5, unit(Vec3::new(0.1, 0.0, 1.0)));
        binding.policy = BundleBranchPolicy::Sequence(&events);
        let branches = [binding];
        let lit = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
        let sample = lit.sample_ray_with_samples(
            cx, ray(0.75), 550.0, 0, settings(), |_| [0.0, 0.0, 0.99],
        ).unwrap();
        // At the axial solution the two refractions change the source angular
        // footprint to (air_length + glass_length/n)^2. Geometric distance^2
        // is wrong even though the central ray itself is straight.
        // Calibration must use the resolved axial launch, not the biased seed.
        let area = (0.75_f64 + 0.25 / 1.5).powi(2);
        let axial_intensity = lamp.flux_lm() * lamp.spectrum().bins_w_per_lumen_nm()[0]
            * 101.0 / (2.0 * PI);
        let normal_metal = ((1.0_f64 - 0.2).powi(2) + 9.0)
            / ((1.0_f64 + 0.2).powi(2) + 9.0);
        let expected = axial_intensity / area * 0.96_f64.powi(2)
            * (-0.1_f64 * (0.75 + 0.25) - 0.3 * 0.25).exp()
            * normal_metal / (4.0 * PI * alpha * alpha);
        near(sample.fixture_radiance_w_m2_sr_nm(), expected);
        near(sample.radiance_sample().radiance_w_m2_sr_nm(), expected);
        assert_eq!(sample.connections(), 1);
        assert!(sample.trace_evaluations() > 1);
    });
}

#[test]
fn g1_point_fixture_rough_metal_retains_non_coplanar_camera_polarization() {
    with_cx(|cx, _| {
        let mirror_normal = unit(Vec3::new(-1.0, 0.0, 1.0));
        let receiver = Point3::new(1.0, 0.0, 1.0);
        let source_point = Point3::new(-1.0, 1.2, 1.8);
        let incoming = unit(receiver.delta_from(source_point));
        let outgoing = incoming.scale(-1.0);
        let geometry = InstanceScene::try_new(vec![
            oriented_panel(1, Point3::new(0.0, 0.0, 1.0), mirror_normal, 0.15),
            oriented_panel(2, receiver, Vec3::new(1.0, 0.0, 0.0), 0.5),
        ]).unwrap();
        let material = glass(1.0, 0.1);
        let media = [SceneMedium::Homogeneous(&material)];
        let conductors = [metal()];
        let alpha = 0.3;
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx, &geometry, &media, &conductors, &[], &[],
            vec![
                SurfaceBinding { object_id: 1, interaction: Interaction::SmoothConductor {
                    optics: 0, medium: 0,
                }},
                SurfaceBinding { object_id: 2, interaction: Interaction::RoughConductor {
                    optics: 0, medium: 0, surface: ConductorSurface::try_rough(alpha).unwrap(),
                }},
            ],
            100_000,
        ).unwrap();
        let camera = CurvedCameraScene::try_new(cx, &scene, &[], None).unwrap();
        let lamp = source(source_point);
        let branches = [branch(&lamp, 2, incoming)];
        let lit = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
        let ordinary = camera.sample_ray_with_samples(
            cx, ray(0.0), 550.0, 0, trace_settings(), |_| [0.0; 3],
        ).unwrap();
        let sample = lit.sample_ray_with_samples(
            cx, ray(0.0), 550.0, 0, settings(), |_| [0.0; 3],
        ).unwrap();
        // Use the native smooth optical operators as an independent Mueller
        // reference, then analytic GGX D, Smith G1 and calibrated source geometry.
        // This does not call the continuous receiver or its connection evaluator.
        let prefix = camera.launch(ray(0.0), 550.0, 0).unwrap().advance(1.0).unwrap()
            .reflect_conductor(mirror_normal, &conductors[0]).unwrap().beam
            .advance(1.0).unwrap();
        let facet = unit(Vec3::new(-1.0 + outgoing.x, outgoing.y, outgoing.z));
        let reflected = prefix.reflect_conductor(facet, &conductors[0]).unwrap().beam;
        let cosine = -outgoing.x;
        let half_cosine = -facet.x;
        let denominator = half_cosine * half_cosine * (alpha * alpha - 1.0) + 1.0;
        let d = alpha * alpha / (PI * denominator * denominator);
        let masking = 2.0 * cosine
            / (cosine + (alpha * alpha + (1.0 - alpha * alpha) * cosine * cosine).sqrt());
        let distance2 = receiver.delta_from(source_point).dot(receiver.delta_from(source_point));
        let irradiance = intensity(&lamp) * cosine / distance2
            * (-0.1 * distance2.sqrt()).exp();
        let expected = reflected.power().i * d * masking / (4.0 * cosine) * irradiance;
        near(sample.fixture_radiance_w_m2_sr_nm(), expected);
        near(sample.radiance_sample().radiance_w_m2_sr_nm(), expected);
        assert_eq!(sample.radiance_sample().path(), ordinary.path());
        let depolarized = PolarizedBeam::try_new(
            prefix.ray().origin, prefix.frame(), 550.0, prefix.medium(),
            Stokes::unpolarized(prefix.power().i),
        ).unwrap().reflect_conductor(facet, &conductors[0]).unwrap().beam;
        assert!((reflected.power().i - depolarized.power().i).abs()
            > 1e-5 * reflected.power().i,
            "depolarizing this non-coplanar metallic path must change its response");
        assert_eq!(sample.connections(), 1);
    });
}

#[test]
fn g0_point_fixture_continuous_receiver_admission_keeps_delta_surfaces_explicit() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(1, 1.0, 0.0, 2.0)]).unwrap();
        let first = glass(1.0, 0.0);
        let second = glass(1.5, 0.0);
        let media = [SceneMedium::Homogeneous(&first), SceneMedium::Homogeneous(&second)];
        let conductors = [metal()];
        for interaction in [
            Interaction::Dielectric { outside: 0, inside: 1, surface: DielectricSurface::SMOOTH },
            Interaction::RoughConductor { optics: 0, medium: 0, surface: ConductorSurface::SMOOTH },
            Interaction::SmoothConductor { optics: 0, medium: 0 },
            Interaction::Opaque,
        ] {
            let scene = CurvedPolarizedScene::try_new_with_materials(
                cx, &geometry, &media, &conductors, &[], &[],
                vec![SurfaceBinding { object_id: 1, interaction }], 100_000,
            ).unwrap();
            let camera = CurvedCameraScene::try_new(cx, &scene, &[], None).unwrap();
            let lamp = source(Point3::new(0.0, 0.0, 0.0));
            let branches = [branch(&lamp, 1, Vec3::new(0.0, 0.0, 1.0))];
            assert!(matches!(
                CurvedPointLightCameraScene::try_new(cx, &camera, &branches),
                Err(PointLightCameraError::InvalidBinding { index: 0 }),
            ));
        }
    });
}
