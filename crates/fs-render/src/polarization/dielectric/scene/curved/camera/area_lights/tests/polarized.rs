//! G1/G3: the original source law reaches both area bundles and rough receivers.

use super::*;
use crate::dielectric::{DielectricEvent, DielectricSurface};
use crate::polarization::dielectric::scene::curved::camera::area_lights::camera::{
    AreaLightBinding, AreaLightCameraSettings, AreaLightCameraWorkspace, CurvedAreaLightCameraScene,
};
use crate::polarization::dielectric::scene::curved::camera::{
    CurvedCameraScene, RadianceSpectrum, SourcePolarization,
};

fn placed_source(placement: RigidTransform) -> GeometryInstance {
    let original = source(0.0);
    GeometryInstance::try_new(
        2,
        original.geometry_identity(),
        original.geometry().clone(),
        placement,
    )
    .unwrap()
}

fn polarized_source(q: f64) -> SurfaceEmitter {
    SurfaceEmitter::new(
        2,
        RadianceSpectrum::try_constant(4.0).unwrap(),
        EmissionSides::Front,
    )
    .with_polarization(
        SourcePolarization::try_constant(
            Vec3::new(0.0, 1.0, 0.0),
            Stokes {
                i: 1.0,
                q,
                u: 0.0,
                v: 0.2,
            },
        )
        .unwrap(),
    )
    .unwrap()
}

#[test]
fn g1_brewster_area_power_uses_source_polarization_and_body_frame_not_caller_basis() {
    with_cx(|cx| {
        // Source +z hits a plane at z=1 at the 1 -> 1.5 Brewster angle.
        // The reflected direction is (12,0,5)/13 and reaches a z-plane
        // detector after another metre. Rs=25/169 and Rp=0 independently
        // of the bundle implementation; the unfolded distance is 2 m.
        let point = Point3::new(12.0 / 13.0, 0.0, 18.0 / 13.0);
        let tangent = Vec3::new(2.0, 0.0, 3.0).scale(1.0 / 13.0_f64.sqrt());
        let air = glass(1.0, 0.0);
        let substrate = glass(1.5, 0.0);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&substrate),
        ];
        let decisions = [DielectricEvent::Reflection];
        let seeds = [BundleImageSeed {
            initial_direction: Vec3::new(0.0, 0.0, 1.0),
            policy: BundleBranchPolicy::Sequence(&decisions),
        }];
        let target =
            BundleConnectionTarget::try_new(9, point, Vec3::new(0.0, 0.0, -1.0), 2e-7).unwrap();
        for angle in [0.0, core::f64::consts::FRAC_PI_2] {
            let placement = RigidTransform::try_new(
                [0.0, 0.0, (0.5 * angle).sin(), (0.5 * angle).cos()],
                [0.0; 3],
            )
            .unwrap();
            let geometry = InstanceScene::try_new(vec![
                mesh(
                    1,
                    vec![
                        [-0.3 * tangent.x, -0.3, 1.0 - 0.3 * tangent.z],
                        [-0.3 * tangent.x, 0.3, 1.0 - 0.3 * tangent.z],
                        [0.6 * tangent.x, 0.0, 1.0 + 0.6 * tangent.z],
                    ],
                    vec![[0, 1, 2]],
                    RigidTransform::identity(),
                ),
                placed_source(placement),
                detector(point),
            ])
            .unwrap();
            let scene = CurvedPolarizedScene::try_new(
                cx,
                &geometry,
                &media,
                vec![
                    SurfaceBinding {
                        object_id: 1,
                        interaction: Interaction::Dielectric {
                            outside: 0,
                            inside: 1,
                            surface: DielectricSurface::SMOOTH,
                        },
                    },
                    SurfaceBinding {
                        object_id: 2,
                        interaction: Interaction::Opaque,
                    },
                    detector_binding(),
                ],
                100_000,
            )
            .unwrap();
            for q in [-0.6, 0.6] {
                let emitter = polarized_source(q);
                let panel = CurvedAreaEmitter::try_new(cx, &scene, &emitter, 1).unwrap();
                let sample = panel.sample_area([4.0 / 9.0, 0.5]).unwrap();
                let mut workspace = BundleImageWorkspace::try_new(1, 1).unwrap();
                let q_in_incidence_frame = q * (2.0 * angle).cos();
                let reflected = (25.0 / 169.0) * 0.5 * (1.0 + q_in_incidence_frame);
                let expected = 4.0 * 0.27 / 4.0 * (5.0 / 13.0) * reflected;
                let mut original = None;
                for fallback in [
                    Vec3::new(1.0, 0.0, 0.0),
                    Vec3::new(0.0, 1.0, 0.0),
                    // The old unpolarized fallback is singular for this ray;
                    // it is irrelevant once the source owns a valid axis.
                    Vec3::new(0.0, 0.0, 1.0),
                ] {
                    let observed = panel
                        .incident_irradiance(
                            cx,
                            sample,
                            BoundarySide::Positive,
                            0,
                            fallback,
                            550.0,
                            target,
                            options(),
                            &seeds,
                            BundleImageOptions {
                                max_images: 1,
                                ..BundleImageOptions::default()
                            },
                            &mut workspace,
                        )
                        .unwrap();
                    let power = observed.spectral_irradiance_w_m2_nm();
                    near(power, expected, 4e-6);
                    if let Some(previous) = original {
                        assert_eq!(power, previous);
                    }
                    original = Some(power);
                    assert_eq!(observed.image_count(), 1);
                }
            }
        }
    });
}

#[test]
fn g1_area_nee_keeps_polarized_fresnel_at_the_native_rough_dielectric_receiver() {
    with_cx(|cx| {
        // The camera sees the receiver normally; the area-source centre is
        // 60 degrees off its normal. Their half-vector has 30-degree facet
        // incidence. A scalar source average cannot reproduce the +/-Q sweep.
        let cosine = 0.5_f64;
        let sine = 3.0_f64.sqrt() / 2.0;
        let alpha = 0.3_f64;
        let placement =
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [-2.0 * sine / cosine, 0.0, -1.0])
                .unwrap();
        let geometry = InstanceScene::try_new(vec![
            placed_source(placement),
            detector(Point3::new(0.0, 0.0, 1.0)),
        ])
        .unwrap();
        let air = glass(1.0, 0.0);
        let substrate = glass(1.5, 0.0);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&substrate),
        ];
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![
                SurfaceBinding {
                    object_id: 2,
                    interaction: Interaction::Opaque,
                },
                SurfaceBinding {
                    object_id: 9,
                    interaction: Interaction::Dielectric {
                        outside: 1,
                        inside: 0,
                        surface: DielectricSurface::try_rough(alpha).unwrap(),
                    },
                },
            ],
            100_000,
        )
        .unwrap();
        let camera_ray = Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        let seeds = [BundleImageSeed {
            initial_direction: Vec3::new(sine, 0.0, cosine),
            policy: BundleBranchPolicy::Sequence(&[]),
        }];
        let config = AreaLightCameraSettings {
            surfaces: options().scene.trace,
            connection: options(),
            max_area_samples: 4,
            max_trace_evaluations: 128,
        };

        let facet_cosine = ((1.0 + cosine) / 2.0).sqrt();
        let transmitted_cosine = (1.0 - (1.0 - facet_cosine * facet_cosine) / 2.25).sqrt();
        let rs = ((facet_cosine - 1.5 * transmitted_cosine)
            / (facet_cosine + 1.5 * transmitted_cosine))
            .powi(2);
        let rp = ((1.5 * facet_cosine - transmitted_cosine)
            / (1.5 * facet_cosine + transmitted_cosine))
            .powi(2);
        let d = alpha * alpha
            / (core::f64::consts::PI
                * (facet_cosine * facet_cosine * (alpha * alpha - 1.0) + 1.0).powi(2));
        let g = 2.0 * cosine
            / (cosine + (alpha * alpha + (1.0 - alpha * alpha) * cosine * cosine).sqrt());
        let distance_squared = 4.0 / (cosine * cosine);
        let p_light = distance_squared / (0.27 * cosine);
        // VNDF at the normal view gives D/4; the native dielectric proposal
        // also includes its unpolarized reflection probability, regardless Q.
        let p_bsdf = 0.5 * (rs + rp) * d / 4.0;
        let weight = p_light * p_light / (p_light * p_light + p_bsdf * p_bsdf);
        let irradiance = 4.0 * 0.27 * cosine * cosine / distance_squared;

        let emitters = [SurfaceEmitter::new(
            2,
            RadianceSpectrum::try_constant(4.0).unwrap(),
            EmissionSides::Front,
        )
        .with_polarization(
            SourcePolarization::try_spectral(
                Vec3::new(0.0, 1.0, 0.0),
                vec![
                    (
                        380.0,
                        Stokes {
                            i: 1.0,
                            q: 0.6,
                            u: 0.0,
                            v: 0.0,
                        },
                    ),
                    (
                        780.0,
                        Stokes {
                            i: 1.0,
                            q: -0.6,
                            u: 0.0,
                            v: 0.0,
                        },
                    ),
                ],
            )
            .unwrap(),
        )
        .unwrap()];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let panel = CurvedAreaEmitter::from_camera(cx, &camera, 2, 1).unwrap();
        // The middle wavelength is exactly unpolarized, but it still belongs
        // to the same declared source frame. A singular caller fallback must
        // neither reject preparation nor replace that frame during a query.
        for (wavelength, q) in [(380.0, 0.6), (580.0, 0.0), (780.0, -0.6)] {
            let baseline = camera
                .sample_ray_with_samples(cx, camera_ray, wavelength, 0, config.surfaces, |_| {
                    [0.0; 3]
                })
                .unwrap();
            assert_eq!(baseline.radiance_w_m2_sr_nm(), 0.0);
            let fresnel = 0.5 * ((1.0 + q) * rs + (1.0 - q) * rp);
            let expected = fresnel * d * g / (4.0 * cosine) * irradiance * weight;
            for fallback in [
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
                seeds[0].initial_direction,
            ] {
                let bindings = [AreaLightBinding {
                    source: &panel,
                    receiver_object_id: 9,
                    source_side: BoundarySide::Positive,
                    initial_medium: 0,
                    source_transverse_axis: fallback,
                    seeds: &seeds,
                    image_options: BundleImageOptions {
                        max_images: 1,
                        ..BundleImageOptions::default()
                    },
                    position_tolerance_m: 2e-7,
                }];
                let owner = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
                let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
                let observed = owner
                    .sample_ray_with_samples(
                        cx,
                        camera_ray,
                        wavelength,
                        0,
                        config,
                        &mut workspace,
                        |_| [0.0; 3],
                        |_, _| [4.0 / 9.0, 0.5],
                    )
                    .unwrap();
                near(observed.area_radiance_w_m2_sr_nm(), expected, 1e-5);
                near(
                    observed.radiance_sample().radiance_w_m2_sr_nm(),
                    expected,
                    1e-5,
                );
                assert_eq!(observed.radiance_sample().path(), baseline.path());
                assert_eq!((observed.area_samples(), observed.images()), (1, 1));
            }
        }
    });
}
