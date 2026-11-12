//! Native directional sources, receiving optics, visibility and MIS.

use super::super::super::{EmissionSides, RadianceSpectrum, SourcePolarization, SurfaceEmitter};
use super::*;
use crate::charts::TriMesh;
use crate::conductor::{
    ConductorDataStatus, ConductorIorSample, ConductorOptics, ConductorSource, ConductorSurface,
};
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricSurface, GlassProvenance,
};
use crate::gradient_index::GradientIndexOptions;
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::lighting::EnvironmentMap;
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, SceneMedium};
use crate::polarization::dielectric::scene::{Interaction, SurfaceBinding};
use crate::polarization::dielectric::transport::PolarizationFrame;
use crate::spectral::SpectralReflectance;
use fs_blake3::{ContentHash, hash_domain};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};

mod microsurface;

#[cfg(feature = "volumes")]
mod participating;

const X: Vec3 = Vec3::new(1.0, 0.0, 0.0);
const Y: Vec3 = Vec3::new(0.0, 1.0, 0.0);
const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);

fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 73,
                kernel_id: 91,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}

fn settings() -> EnvironmentCameraSettings {
    let event = GradientEventOptions {
        integration: GradientIndexOptions {
            position_tolerance_m: 1e-11,
            direction_tolerance: 1e-12,
            optical_path_tolerance_m: 1e-11,
            optical_depth_tolerance: 1e-12,
            ..GradientIndexOptions::default()
        },
        surface_tolerance_m: 1e-8,
        ..GradientEventOptions::default()
    };
    EnvironmentCameraSettings {
        surfaces: CurvedTraceSettings {
            max_surface_events: 8,
            max_distance_m: 4.0,
            event,
            ..CurvedTraceSettings::default()
        },
        visibility: event,
        visibility_distance_m: 4.0,
        max_environment_samples: 8,
        max_query_work: 1_000_000,
    }
}

fn material(index: f64, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(index).unwrap(),
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}

fn patch(id: u64, center_x: f64, z: f64, radius: f64) -> GeometryInstance {
    GeometryInstance::try_new(
        id,
        hash_domain("environment-nee-native", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(
            vec![
                [center_x - radius, -radius, z],
                [center_x, 2.0 * radius, z],
                [center_x + radius, -radius, z],
            ],
            vec![[0, 1, 2]],
        )),
        RigidTransform::identity(),
    )
    .unwrap()
}

fn camera_ray() -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, 0.0),
        dir: Z,
    }
}
fn state(q: f64, u: f64, v: f64) -> Stokes {
    Stokes { i: 1.0, q, u, v }
}
fn near(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * b.abs().max(1e-14),
        "{a:.17e} != {b:.17e}"
    );
}

#[test]
fn original_texel_pdf_spectrum_and_rotated_source_frame_remain_one_proposal() {
    with_cx(|cx| {
        let map = EnvironmentMap::try_from_linear_srgb(
            2,
            2,
            vec![[1.0; 3], [3.0; 3], [2.0; 3], [4.0; 3]],
            0.37,
        )
        .unwrap();
        let laws: Vec<_> = (0..4)
            .map(|i| {
                SourcePolarization::try_spectral(
                    X,
                    vec![
                        (380.0, state(0.1 * i as f64, 0.1, 0.4)),
                        (780.0, state(-0.1 * i as f64, -0.2, -0.3)),
                    ],
                )
                .unwrap()
            })
            .collect();
        let source = EnvironmentRadiance::try_new(&map, 7.0)
            .unwrap()
            .with_polarization_map(cx, &laws)
            .unwrap();
        for samples in [[0.0, 0.3], [0.17, 0.21], [0.5, 0.3], [0.95, 0.91]] {
            let (original, texel) = map.sample_with_texel(samples[0], samples[1]).unwrap();
            for wavelength in [380.0, 580.0, 780.0] {
                let sample = source
                    .sample_direction(wavelength, samples)
                    .unwrap()
                    .unwrap();
                assert_eq!(sample.pdf, original.pdf_solid_angle);
                assert_eq!(
                    sample.stokes,
                    laws[texel].normalized_stokes(wavelength).unwrap()
                );
                assert_eq!(
                    sample.radiance,
                    original.emission.0.eval(wavelength) * original.emission.1 * 7.0
                );
                let axis = Vec3::new(0.37_f64.cos(), 0.0, 0.37_f64.sin());
                let expected =
                    PolarizationFrame::try_new(sample.direction.scale(-1.0), axis).unwrap();
                near(
                    sample.frame.horizontal().dot(expected.horizontal()),
                    1.0,
                    2e-13,
                );
            }
        }
        assert!(matches!(
            source.sample_direction(550.0, [1.0, 0.2]),
            Err(EnvironmentCameraError::InvalidSample)
        ));
        assert!(source.sample_direction(850.0, [0.2, 0.3]).is_err());
        for (light, continuation) in [(1.0, 4.0), (1e200, 1e200), (1e-200, 2e-200)] {
            near(
                power_weight(light, continuation).unwrap()
                    + power_weight(continuation, light).unwrap(),
                1.0,
                2e-15,
            );
        }
        near(
            positive_product(&[1e200, 1e200, 1e-200], 1e200).unwrap(),
            1.0,
            2e-12,
        );
        assert!(positive_product(&[f64::from_bits(1)], 4.0).is_err());
    });
}

#[test]
fn native_diffuse_environment_mis_preserves_path_and_unweighted_surface_emission() {
    with_cx(|cx| {
        let map = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap();
        let polarization = SourcePolarization::try_constant(X, state(0.6, 0.2, 0.3)).unwrap();
        let source = EnvironmentRadiance::try_new(&map, 3.0)
            .unwrap()
            .with_polarization(&polarization)
            .unwrap();
        let geometry = InstanceScene::try_new(vec![patch(9, 0.0, 1.0, 0.3)]).unwrap();
        let air = material(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let diffuse = [SpectralReflectance::constant(0.6).unwrap()];
        let native = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &diffuse,
            vec![SurfaceBinding {
                object_id: 9,
                interaction: Interaction::DepolarizingDiffuse {
                    medium: 0,
                    reflectance: 0,
                },
            }],
            100_000,
        )
        .unwrap();
        let emitters = [SurfaceEmitter::new(
            9,
            RadianceSpectrum::try_constant(5.0).unwrap(),
            EmissionSides::Front,
        )];
        let camera = CurvedCameraScene::try_new(cx, &native, &emitters, Some(source)).unwrap();
        let owner = CurvedEnvironmentCameraScene::try_new(cx, &camera).unwrap();
        let observed = owner
            .sample_ray_with_samples(
                cx,
                camera_ray(),
                550.0,
                0,
                settings(),
                |event| {
                    assert_eq!(event, 0);
                    [0.0, 0.0, 0.0]
                },
                |event| {
                    assert_eq!(event, 0);
                    [0.5, 0.75]
                },
            )
            .unwrap();
        let baseline = camera
            .sample_ray_with_samples(cx, camera_ray(), 550.0, 0, settings().surfaces, |_| {
                [0.0, 0.0, 0.0]
            })
            .unwrap();
        let le = source.evaluate(Z.scale(-1.0), 550.0).unwrap();
        near(baseline.radiance_w_m2_sr_nm(), 5.0 + 0.6 * le, 4e-6);
        near(
            observed.environment_radiance_w_m2_sr_nm(),
            0.6 * le * 4.0 / 17.0,
            4e-6,
        );
        near(observed.environment_endpoint_weight(), 16.0 / 17.0, 2e-13);
        near(
            observed.radiance_sample().radiance_w_m2_sr_nm(),
            5.0 + 0.6 * le * 20.0 / 17.0,
            4e-6,
        );
        assert_eq!(observed.radiance_sample().path(), baseline.path());
        assert_eq!(
            (
                observed.environment_samples(),
                observed.visibility_queries()
            ),
            (1, 1)
        );
        assert!(observed.field_evaluations() > 0 && observed.query_work() > 0);
        let analyzed = owner
            .sample_ray_analyzed_with_samples(
                cx,
                camera_ray(),
                550.0,
                0,
                CameraAnalyzer::try_linear(X, 0.0).unwrap(),
                settings(),
                |_| [0.0; 3],
                |_| [0.5, 0.75],
            )
            .unwrap();
        near(
            analyzed.radiance_sample().radiance_w_m2_sr_nm(),
            0.5 * observed.radiance_sample().radiance_w_m2_sr_nm(),
            4e-6,
        );
        let insufficient = EnvironmentCameraSettings {
            visibility_distance_m: 1e-10,
            ..settings()
        };
        assert!(
            owner
                .sample_ray_with_samples(
                    cx,
                    camera_ray(),
                    550.0,
                    0,
                    insufficient,
                    |_| [0.0; 3],
                    |_| [0.5, 0.75]
                )
                .is_err()
        );
    });
}

#[test]
fn rough_metal_environment_nee_retains_linear_extinction_and_circular_handedness() {
    with_cx(|cx| {
        let map = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap();
        let geometry = InstanceScene::try_new(vec![patch(9, 0.0, 1.0, 0.3)]).unwrap();
        let air = material(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let metals = [ConductorOptics::try_new(
            core::array::from_fn(|i| {
                ConductorIorSample::try_new(380.0 + 50.0 * i as f64, 0.3, 3.0).unwrap()
            }),
            ConductorSource::try_new(ContentHash([135; 32]), ConductorDataStatus::Representative)
                .unwrap(),
        )
        .unwrap()];
        let alpha = 0.4_f64;
        let native = CurvedPolarizedScene::try_new_with_conductors(
            cx,
            &geometry,
            &media,
            &metals,
            vec![SurfaceBinding {
                object_id: 9,
                interaction: Interaction::RoughConductor {
                    medium: 0,
                    optics: 0,
                    surface: ConductorSurface::try_rough(alpha).unwrap(),
                },
            }],
            100_000,
        )
        .unwrap();
        let reflection = ((0.3_f64 - 1.0).powi(2) + 9.0) / ((0.3_f64 + 1.0).powi(2) + 9.0);
        for (stokes, analyzer, selected) in [
            (
                state(1.0, 0.0, 0.0),
                CameraAnalyzer::try_linear(X, 0.0).unwrap(),
                true,
            ),
            (
                state(-1.0, 0.0, 0.0),
                CameraAnalyzer::try_linear(X, 0.0).unwrap(),
                false,
            ),
            (
                state(0.0, 0.0, -1.0),
                CameraAnalyzer::try_circular(X, true).unwrap(),
                true,
            ),
            (
                state(0.0, 0.0, 1.0),
                CameraAnalyzer::try_circular(X, true).unwrap(),
                false,
            ),
        ] {
            let law = SourcePolarization::try_constant(X, stokes).unwrap();
            let source = EnvironmentRadiance::try_new(&map, 2.0)
                .unwrap()
                .with_polarization(&law)
                .unwrap();
            let camera = CurvedCameraScene::try_new(cx, &native, &[], Some(source)).unwrap();
            let owner = CurvedEnvironmentCameraScene::try_new(cx, &camera).unwrap();
            let observed = owner
                .sample_ray_analyzed_with_samples(
                    cx,
                    camera_ray(),
                    550.0,
                    0,
                    analyzer,
                    settings(),
                    |_| [0.0; 3],
                    |_| [0.5, 0.75],
                )
                .unwrap();
            let le = source.evaluate(Z.scale(-1.0), 550.0).unwrap();
            // NEE at the normal view: f cos / p_env = F D(n) / (4 p_env) = F Le / a^2, with
            // D(n) = 1 / (pi a^2) and the uniform map's p_env = 1 / (4 pi), times the power
            // heuristic against the reflection proposal. The bounded spherical-cap proposal
            // (6a13950eae; tracer.rs ggx_vndf_reflection_pdf) has density D(n) / 4 divided by its
            // cap mass 1 - (1 - k) G1(wo) / 2, which at this view (c = 1, s = 1, G1(wo) = 1,
            // k = (1 - a^2) / (1 + a^2)) is 1 / (1 + a^2).
            let a2 = alpha * alpha;
            let k = (1.0 - a2) / (1.0 + a2);
            let cap_mass = 1.0 - 0.5 * (1.0 - k);
            let p_env = 1.0 / (4.0 * core::f64::consts::PI);
            let p_proposal = 1.0 / (4.0 * core::f64::consts::PI * a2) / cap_mass;
            let weight = p_env * p_env / (p_env * p_env + p_proposal * p_proposal);
            let expected = if selected {
                reflection * le / a2 * weight
            } else {
                0.0
            };
            near(observed.environment_radiance_w_m2_sr_nm(), expected, 4e-6);
            if !selected {
                assert_eq!(observed.environment_radiance_w_m2_sr_nm(), 0.0);
            }
        }
    });
}

#[test]
fn native_first_boundary_blocks_shadow_and_delta_departure_removes_endpoint_competitor() {
    with_cx(|cx| {
        let map = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap();
        let source = EnvironmentRadiance::try_new(&map, 3.0).unwrap();
        let air = material(1.0, 0.0);
        let glass = material(1.5, 0.0);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&glass),
        ];
        let diffuse = [SpectralReflectance::constant(0.6).unwrap()];
        for mirror in [false, true] {
            let geometry =
                InstanceScene::try_new(vec![patch(1, 0.8, 0.4, 0.15), patch(9, 0.0, 1.0, 0.3)])
                    .unwrap();
            let first = if mirror {
                Interaction::Dielectric {
                    outside: 1,
                    inside: 0,
                    surface: DielectricSurface::SMOOTH,
                }
            } else {
                Interaction::Opaque
            };
            let native = CurvedPolarizedScene::try_new_with_materials(
                cx,
                &geometry,
                &media,
                &[],
                &[],
                &diffuse,
                vec![
                    SurfaceBinding {
                        object_id: 1,
                        interaction: first,
                    },
                    SurfaceBinding {
                        object_id: 9,
                        interaction: Interaction::DepolarizingDiffuse {
                            medium: 0,
                            reflectance: 0,
                        },
                    },
                ],
                100_000,
            )
            .unwrap();
            let camera = CurvedCameraScene::try_new(cx, &native, &[], Some(source)).unwrap();
            let owner = CurvedEnvironmentCameraScene::try_new(cx, &camera).unwrap();
            // The narrow blocker misses the primary camera ray. A sampled
            // direction (0.8,0,-0.6) reaches it after leaving the receiver.
            let phi = (-0.6_f64).atan2(0.8).rem_euclid(core::f64::consts::TAU);
            let environment = [0.5, phi / core::f64::consts::TAU];
            let continuation = |event| {
                if mirror && event == 0 {
                    [0.64, 0.25, 0.0]
                } else {
                    [0.0; 3]
                }
            };
            let observed = owner
                .sample_ray_with_samples(
                    cx,
                    camera_ray(),
                    550.0,
                    0,
                    settings(),
                    continuation,
                    |_| environment,
                )
                .unwrap();
            let baseline = camera
                .sample_ray_with_samples(
                    cx,
                    camera_ray(),
                    550.0,
                    0,
                    settings().surfaces,
                    continuation,
                )
                .unwrap();
            assert_eq!(observed.environment_radiance_w_m2_sr_nm(), 0.0);
            assert_eq!(observed.visibility_queries(), 1);
            assert_eq!(observed.radiance_sample().path(), baseline.path());
            if mirror {
                assert_eq!(baseline.path().surface_events(), 2);
                assert_eq!(observed.environment_endpoint_weight(), 1.0);
                assert_eq!(observed.radiance_sample(), &baseline);
            } else {
                near(observed.environment_endpoint_weight(), 16.0 / 17.0, 2e-13);
                near(
                    observed.radiance_sample().radiance_w_m2_sr_nm(),
                    baseline.radiance_w_m2_sr_nm() * 16.0 / 17.0,
                    4e-6,
                );
            }
        }
    });
}
