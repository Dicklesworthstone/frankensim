//! Distinct slab orders illuminating the original continuous surface laws.

use super::*;
use crate::camera::{Aperture, CameraProjection, PhysicalCamera};
use crate::conductor::{ConductorOptics, ConductorSurface};
use crate::dielectric::DielectricEvent;
use crate::polarization::dielectric::scene::curved::camera::raster::point_lights::PointLightRasterError;
use crate::polarization::dielectric::scene::curved::camera::raster::{
    CAMERA_WAVELENGTHS, CurvedCameraRaster,
};
use crate::surface_detail::{BrushingPattern, SurfaceDetail};
use core::f64::consts::{FRAC_PI_2, PI};

const X: f64 = 0.2;
const ABSORPTIONS: [f64; 3] = [0.11, 0.29, 0.07];
const SCATTERING: [f64; 3] = [0.17, 0.23, 0.0];
const DIRECT: [DielectricEvent; 2] = [DielectricEvent::Transmission; 2];
const ROUND_TRIP: [DielectricEvent; 4] = [
    DielectricEvent::Transmission,
    DielectricEvent::Reflection,
    DielectricEvent::Reflection,
    DielectricEvent::Transmission,
];

fn slab_geometry() -> InstanceScene {
    InstanceScene::try_new(
        [
            (2, 0.25, false),
            (3, 0.5, true),
            (5, 1.0, true),
            (9, -0.5, true),
        ]
        .map(|(id, z, upward)| {
            GeometryInstance::try_new(
                id,
                hash_domain("rough-image-slab", &id.to_le_bytes()),
                SharedGeometry::mesh(TriMesh::new(
                    vec![[X - 3.0, -3.0, z], [X + 3.0, -3.0, z], [X, 3.0, z]],
                    vec![if upward { [0, 1, 2] } else { [0, 2, 1] }],
                )),
                RigidTransform::identity(),
            )
            .unwrap()
        })
        .to_vec(),
    )
    .unwrap()
}

fn fluid(index: f64, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(index).unwrap(),
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}

fn slab_bindings(receiver: Interaction) -> Vec<SurfaceBinding> {
    vec![
        SurfaceBinding {
            object_id: 2,
            interaction: Interaction::Dielectric {
                outside: 0,
                inside: 1,
                surface: DielectricSurface::SMOOTH,
            },
        },
        SurfaceBinding {
            object_id: 3,
            interaction: Interaction::Dielectric {
                outside: 0,
                inside: 1,
                surface: DielectricSurface::SMOOTH,
            },
        },
        SurfaceBinding {
            object_id: 5,
            interaction: receiver,
        },
        SurfaceBinding {
            object_id: 9,
            interaction: Interaction::Opaque,
        },
    ]
}

fn slab_source(x: f64) -> PhotometricSource {
    PhotometricSource::try_new(
        RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [x, 0.0, 0.0]).unwrap(),
        100.0,
        1.0,
        AngularDistribution::cosine_power(8.0, FRAC_PI_2).unwrap(),
        SourceSpectrum::try_from_bins([1.0; SPECTRAL_BINS]).unwrap(),
    )
    .unwrap()
}

fn seeds() -> [BundleImageSeed<'static>; 3] {
    [
        BundleBranchPolicy::Sequence(&DIRECT),
        BundleBranchPolicy::Sequence(&ROUND_TRIP),
        BundleBranchPolicy::Sequence(&DIRECT),
    ]
    .map(|policy| BundleImageSeed {
        initial_direction: Vec3::new(0.0, 0.0, 1.0),
        policy,
    })
}

fn image_group<'a>(
    source: &'a PhotometricSource,
    seeds: &'a [BundleImageSeed<'a>],
) -> PointLightImageGroup<'a> {
    PointLightImageGroup {
        source,
        receiver: PointLightImageReceiver::SurfaceObject(5),
        initial_medium: 0,
        source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
        seeds,
        image_options: BundleImageOptions {
            max_images: 2,
            max_trace_evaluations: 256,
            ..BundleImageOptions::default()
        },
        position_tolerance_m: 2e-7,
    }
}

fn image_settings() -> RayleighPointLightCameraSettings {
    let mut limits = settings();
    limits.point_lights.max_connections = 2;
    limits.point_lights.max_trace_evaluations = 256;
    limits.max_source_free_flight_draws = 0;
    limits.max_source_candidates = 0;
    limits
}

#[test]
fn g1_incoherent_slab_orders_keep_rough_material_radiometry_and_source_survival() {
    with_cx(|cx| {
        let geometry = slab_geometry();
        let air = fluid(1.0, ABSORPTIONS[0]);
        let slab = fluid(1.5, ABSORPTIONS[1]);
        let outside = fluid(1.33, ABSORPTIONS[2]);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&slab),
            SceneMedium::Homogeneous(&outside),
        ];
        let conductors = [ConductorOptics::representative_tungsten()];
        let metal = conductors[0].complex_ior(550.0).unwrap();
        let metal_r = ((1.0 - metal.eta()).powi(2) + metal.k().powi(2))
            / ((1.0 + metal.eta()).powi(2) + metal.k().powi(2));
        let alpha = 0.2;
        let glass_r = ((1.0_f64 - 1.33) / (1.0 + 1.33)).powi(2);
        let glass_law = Interaction::Dielectric {
            outside: 2,
            inside: 0,
            surface: DielectricSurface::try_rough(alpha).unwrap(),
        };
        // X is away from the radial brush's isotropic pole.
        let detail = SurfaceDetail::try_new(0.12, 3.0, BrushingPattern::Radial, 0.0, None).unwrap();
        let (ax, ay) = detail.alpha_pair();
        let cases = [
            (
                glass_law,
                0,
                0.75,
                1.0,
                [0.0, 0.0, 0.99],
                glass_r / (4.0 * PI * alpha * alpha),
            ),
            (
                glass_law,
                2,
                1.25,
                -1.0,
                [0.0; 3],
                (1.0 - glass_r) / (PI * alpha * alpha * (1.0_f64 - 1.0 / 1.33).powi(2)),
            ),
            (
                Interaction::RoughConductor {
                    optics: 0,
                    medium: 0,
                    surface: ConductorSurface::try_rough(alpha).unwrap(),
                },
                0,
                0.75,
                1.0,
                [0.0, 0.0, 0.99],
                metal_r / (4.0 * PI * alpha * alpha),
            ),
            (
                Interaction::BrushedConductor {
                    optics: 0,
                    medium: 0,
                    detail,
                },
                0,
                0.75,
                1.0,
                [0.0, 0.0, 0.99],
                metal_r / (4.0 * PI * ax * ay),
            ),
        ];
        let source = slab_source(X);
        let seeds = seeds();
        let groups = [image_group(&source, &seeds)];
        let spectrum = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [
            Some(CurvedRayleighMedium::try_new(&spectrum, SCATTERING[0], &phase).unwrap()),
            Some(CurvedRayleighMedium::try_new(&spectrum, SCATTERING[1], &phase).unwrap()),
            None,
        ];
        for (case, (interaction, medium, z, direction, surface_draw, bsdf)) in
            cases.into_iter().enumerate()
        {
            let scene = CurvedPolarizedScene::try_new_with_materials(
                cx,
                &geometry,
                &media,
                &conductors,
                &[],
                &[],
                slab_bindings(interaction),
                100_000,
            )
            .unwrap();
            let camera = CurvedCameraScene::try_new(cx, &scene, &[], None).unwrap();
            let transport = CurvedRayleighScene::try_new(cx, &scene, &volumes).unwrap();
            let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
            let fixtures = CurvedPointLightCameraScene::try_new(cx, &camera, &[]).unwrap();
            let owner =
                CurvedRayleighPointLightCameraScene::try_new(cx, &molecular, &fixtures).unwrap();
            let images = CurvedRayleighImageLightCameraScene::try_new(cx, &owner, &groups).unwrap();
            let ray = Ray {
                origin: Point3::new(X, 0.0, z),
                dir: Vec3::new(0.0, 0.0, direction),
            };
            let ordinary = molecular
                .sample_ray_with_tracking_samples(
                    cx,
                    ray,
                    550.0,
                    medium,
                    camera_settings(),
                    |_| surface_draw,
                    |_| [0.99, 0.5, 0.2, 0.9],
                )
                .unwrap();
            let mut workspace = CurvedTransmittanceImageWorkspace::try_new(2, 3).unwrap();
            let sample = |limits, workspace: &mut CurvedTransmittanceImageWorkspace| {
                images
                    .sample_ray_with_tracking_samples(
                        cx,
                        ray,
                        550.0,
                        medium,
                        limits,
                        workspace,
                        |_| surface_draw,
                        |_| [0.99, 0.5, 0.2, 0.9],
                        |_, _| panic!("uniform survival consumes no source variates"),
                    )
                    .unwrap()
            };
            let combined = sample(image_settings(), &mut workspace);
            // Independent unfolded-slab image geometry and normal Fresnel.
            // Each round trip adds two POSITIVE glass lever arms and R^2.
            // Camera survival is already sampled; only source survival is
            // multiplied here. This is an incoherent sum, not film interference.
            let intensity =
                source.flux_lm() * source.spectrum().bins_w_per_lumen_nm()[0] * 9.0 / (2.0 * PI);
            let source_sum = (0..2)
                .map(|q| {
                    let glass_length = (2 * q + 1) as f64 * 0.25;
                    let area = (0.75 + glass_length / 1.5).powi(2);
                    intensity * 0.96_f64.powi(2) * 0.04_f64.powi(2 * q) / area
                        * (-(ABSORPTIONS[0] + SCATTERING[0]) * 0.75
                            - (ABSORPTIONS[1] + SCATTERING[1]) * glass_length)
                            .exp()
                })
                .sum::<f64>();
            let expected = bsdf * source_sum * (-ABSORPTIONS[medium] * 0.25).exp();
            near(combined.fixture_radiance_w_m2_sr_nm(), expected);
            near(combined.radiance_sample().radiance_w_m2_sr_nm(), expected);
            assert_eq!(
                combined.radiance_sample().path(),
                ordinary.radiance_sample().path()
            );
            assert_eq!(combined.connections(), 2);
            assert_eq!(combined.volume_events(), 0);
            assert_eq!(combined.free_flight_draws(), ordinary.free_flight_draws());
            assert_eq!(combined.source_free_flight_draws(), 0);
            assert_eq!(combined.source_candidates(), 0);
            let mut exact = image_settings();
            exact.point_lights.max_trace_evaluations = combined.trace_evaluations();
            assert_eq!(sample(exact, &mut workspace), combined);

            if case == 0 {
                let alias = PointLightImageGroup {
                    receiver: PointLightImageReceiver::DiffuseObject(5),
                    ..groups[0]
                };
                assert_eq!(alias.receiver, groups[0].receiver);
                assert!(matches!(
                    CurvedRayleighImageLightCameraScene::try_new(cx, &owner, &[groups[0], alias]),
                    Err(PointLightCameraError::DuplicateBinding {
                        first: 0,
                        second: 1
                    })
                ));
                let branches = [PointLightBranch {
                    source: &source,
                    receiver_object_id: 5,
                    initial_medium: 0,
                    initial_direction: seeds[0].initial_direction,
                    source_transverse_axis: groups[0].source_transverse_axis,
                    policy: seeds[0].policy,
                    position_tolerance_m: groups[0].position_tolerance_m,
                }];
                let bound = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
                let bound =
                    CurvedRayleighPointLightCameraScene::try_new(cx, &molecular, &bound).unwrap();
                for group in [groups[0], alias] {
                    assert!(matches!(
                        CurvedRayleighImageLightCameraScene::try_new(cx, &bound, &[group]),
                        Err(PointLightCameraError::DuplicateBinding {
                            first: 0,
                            second: 1
                        })
                    ));
                }
                let aliases = [alias];
                let alias_owner =
                    CurvedRayleighImageLightCameraScene::try_new(cx, &owner, &aliases).unwrap();
                let alias_sample = alias_owner
                    .sample_ray_with_tracking_samples(
                        cx,
                        ray,
                        550.0,
                        medium,
                        image_settings(),
                        &mut workspace,
                        |_| surface_draw,
                        |_| [0.99, 0.5, 0.2, 0.9],
                        |_, _| panic!("uniform survival consumes no source variates"),
                    )
                    .unwrap();
                assert_eq!(alias_sample, combined);
                let mut short = image_settings();
                short.point_lights.max_connections = 1;
                assert!(
                    images
                        .sample_ray_with_tracking_samples(
                            cx,
                            ray,
                            550.0,
                            medium,
                            short,
                            &mut workspace,
                            |_| surface_draw,
                            |_| [0.99, 0.5, 0.2, 0.9],
                            |_, _| 0.5
                        )
                        .is_err()
                );
            }
        }
    });
}

#[test]
fn g3_rough_transmitting_image_pixels_sum_original_spectral_branches() {
    with_cx(|cx| {
        let geometry = slab_geometry();
        let air = fluid(1.0, ABSORPTIONS[0]);
        let slab = DielectricGlass::new(
            CauchyIor::try_new(1.47, 0.004, 0.0).unwrap(),
            BeerLambertAbsorption::try_constant(ABSORPTIONS[1]).unwrap(),
            GlassProvenance::Custom,
        );
        let outside = fluid(1.33, ABSORPTIONS[2]);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&slab),
            SceneMedium::Homogeneous(&outside),
        ];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &[],
            slab_bindings(Interaction::Dielectric {
                outside: 2,
                inside: 0,
                surface: DielectricSurface::try_rough(0.2).unwrap(),
            }),
            100_000,
        )
        .unwrap();
        let camera = CurvedCameraScene::try_new(cx, &scene, &[], None).unwrap();
        let volumes = [None, None, None];
        let transport = CurvedRayleighScene::try_new(cx, &scene, &volumes).unwrap();
        let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let fixtures = CurvedPointLightCameraScene::try_new(cx, &camera, &[]).unwrap();
        let owner =
            CurvedRayleighPointLightCameraScene::try_new(cx, &molecular, &fixtures).unwrap();
        // The displaced fixture makes TT and TRRT arrive with different
        // directions, calibrated angular intensities and Fresnel polarization.
        let source = slab_source(X - 0.5);
        let seeds = seeds();
        let groups = [image_group(&source, &seeds)];
        let images = CurvedRayleighImageLightCameraScene::try_new(cx, &owner, &groups).unwrap();
        let physical = PhysicalCamera::try_look_at(
            Point3::new(X, 0.0, 1.25),
            Point3::new(X, 0.0, 1.0),
            Vec3::new(0.0, 1.0, 0.0),
            CameraProjection::try_half_tangent(0.04).unwrap(),
            0.25,
            Aperture::try_circular(0.0).unwrap(),
        )
        .unwrap();
        let raster = CurvedCameraRaster::try_new(&physical, 3, 2).unwrap();
        let mut workspace = CurvedTransmittanceImageWorkspace::try_new(2, 3).unwrap();
        let limits = image_settings();
        // This seed chooses moderate receiver facets and transmits the two
        // smooth slab contacts in every lane, so each camera visits the
        // receiving surface once while both source image orders still count.
        let (seed, pixel, sample) = (91, 4, 9);
        let combined = raster
            .sample_pixel_rayleigh_image_lights(
                cx,
                &images,
                2,
                limits,
                &mut workspace,
                seed,
                pixel,
                sample,
            )
            .unwrap();
        let ordinary = raster
            .sample_pixel_rayleigh(cx, &molecular, 2, camera_settings(), seed, pixel, sample)
            .unwrap();
        let singles = [&DIRECT[..], &ROUND_TRIP[..]].map(|events| {
            let branches = [PointLightBranch {
                source: &source,
                receiver_object_id: 5,
                initial_medium: 0,
                initial_direction: Vec3::new(0.0, 0.0, 1.0),
                source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
                policy: BundleBranchPolicy::Sequence(events),
                position_tolerance_m: 2e-7,
            }];
            let fixtures = CurvedPointLightCameraScene::try_new(cx, &camera, &branches).unwrap();
            let single =
                CurvedRayleighPointLightCameraScene::try_new(cx, &molecular, &fixtures).unwrap();
            raster
                .sample_pixel_rayleigh_point_lights(cx, &single, 2, limits, seed, pixel, sample)
                .unwrap()
        });
        for lane in 0..CAMERA_WAVELENGTHS {
            let expected = singles
                .iter()
                .map(|single| single.point_light_sample().fixture_radiance_w_m2_sr_nm()[lane])
                .sum();
            assert!(singles.iter().all(|single| {
                single.point_light_sample().fixture_radiance_w_m2_sr_nm()[lane] > 0.0
            }));
            near(
                combined.point_light_sample().fixture_radiance_w_m2_sr_nm()[lane],
                expected,
            );
            near(
                combined.pixel_sample().spectral_samples()[lane].radiance_w_m2_sr_nm(),
                expected,
            );
            assert_eq!(
                combined.pixel_sample().spectral_samples()[lane].path(),
                ordinary.pixel_sample().spectral_samples()[lane].path()
            );
            assert_eq!(combined.point_light_sample().connections()[lane], 2);
        }
        let expected_xyz = core::array::from_fn::<_, 3, _>(|axis| {
            singles
                .iter()
                .map(|single| single.pixel_sample().xyz()[axis])
                .sum::<f64>()
        });
        for (actual, expected) in combined.pixel_sample().xyz().into_iter().zip(expected_xyz) {
            near(actual, expected);
        }
        assert_eq!(
            combined,
            raster
                .sample_pixel_rayleigh_image_lights(
                    cx,
                    &images,
                    2,
                    limits,
                    &mut workspace,
                    seed,
                    pixel,
                    sample
                )
                .unwrap()
        );
        let mut too_small = CurvedTransmittanceImageWorkspace::try_new(1, 1).unwrap();
        assert!(matches!(
            raster.sample_pixel_rayleigh_image_lights(
                cx,
                &images,
                2,
                limits,
                &mut too_small,
                seed,
                pixel,
                sample
            ),
            Err(PointLightRasterError::Transport { lane: 0, .. })
        ));
    });
}
