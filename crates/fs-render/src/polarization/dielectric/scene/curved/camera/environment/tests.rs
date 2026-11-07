//! G0/G1/G3: native escaped sources, physical frames and spectral support.

use super::super::{CameraAnalyzer, CurvedCameraScene, RadianceEndpoint};
use super::*;
use crate::charts::{Ray, TriMesh};
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricSurface, GlassProvenance,
};
use crate::gradient_index::{GradientIndexOptions, event::GradientEventOptions};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::lighting::EnvironmentMap;
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::curved::{
    CurvedPolarizedScene, CurvedTraceSettings, SceneMedium,
};
use crate::polarization::dielectric::scene::{Interaction, SurfaceBinding};
use fs_blake3::hash_domain;
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_geom::Point3;

const X: Vec3 = Vec3::new(1.0, 0.0, 0.0);
const Y: Vec3 = Vec3::new(0.0, 1.0, 0.0);
const Z: Vec3 = Vec3::new(0.0, 0.0, 1.0);

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 67,
                    kernel_id: 91,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        );
    });
}

fn material(index: f64, absorption: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(index).unwrap(),
        BeerLambertAbsorption::try_constant(absorption).unwrap(),
        GlassProvenance::Custom,
    )
}

fn settings() -> CurvedTraceSettings {
    CurvedTraceSettings {
        max_surface_events: 8,
        max_distance_m: 4.0,
        max_total_steps: 16_384,
        max_query_work: 1_000_000,
        event: GradientEventOptions {
            integration: GradientIndexOptions {
                position_tolerance_m: 1e-11,
                direction_tolerance: 1e-12,
                optical_path_tolerance_m: 1e-11,
                optical_depth_tolerance: 1e-12,
                ..GradientIndexOptions::default()
            },
            surface_tolerance_m: 1e-8,
            ..GradientEventOptions::default()
        },
    }
}

fn triangle(id: u64, vertices: [[f64; 3]; 3]) -> GeometryInstance {
    GeometryInstance::try_new(
        id,
        hash_domain("polarized-environment-test", &id.to_le_bytes()),
        SharedGeometry::mesh(TriMesh::new(vertices.to_vec(), vec![[0, 1, 2]])),
        RigidTransform::identity(),
    )
    .unwrap()
}

fn with_camera(
    cx: &Cx<'_>,
    source: EnvironmentRadiance<'_>,
    absorption: f64,
    run: impl FnOnce(&CurvedCameraScene<'_, '_>),
) {
    let geometry = InstanceScene::try_new(vec![triangle(
        9,
        [[10.0, 10.0, 10.0], [11.0, 10.0, 10.0], [10.0, 11.0, 10.0]],
    )])
    .unwrap();
    let air = material(1.0, absorption);
    let media = [SceneMedium::Homogeneous(&air)];
    let scene = CurvedPolarizedScene::try_new(
        cx,
        &geometry,
        &media,
        vec![SurfaceBinding {
            object_id: 9,
            interaction: Interaction::Opaque,
        }],
        100_000,
    )
    .unwrap();
    let camera = CurvedCameraScene::try_new(cx, &scene, &[], Some(source)).unwrap();
    run(&camera);
}

fn ray(direction: Vec3) -> Ray {
    Ray {
        origin: Point3::new(0.0, 0.0, 0.0),
        dir: direction,
    }
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 4e-12 * b.abs().max(1.0), "{a:e} != {b:e}");
}
fn state(q: f64, u: f64, v: f64) -> Stokes {
    Stokes { i: 1.0, q, u, v }
}

#[test]
fn g1_escaped_source_retains_full_stokes_and_exact_extinction() {
    with_cx(|cx, _| {
        let map = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[2.0; 3]], 0.0).unwrap();
        let scalar = EnvironmentRadiance::try_new(&map, 3.0).unwrap();
        let emitted = scalar.evaluate(Z, 550.0).unwrap();
        let analyzers = [
            CameraAnalyzer::try_linear(X, 0.0).unwrap(),
            CameraAnalyzer::try_new(
                X,
                Stokes {
                    i: 0.5,
                    q: -0.5,
                    u: 0.0,
                    v: 0.0,
                },
            )
            .unwrap(),
            CameraAnalyzer::try_new(
                X,
                Stokes {
                    i: 0.5,
                    q: 0.0,
                    u: 0.5,
                    v: 0.0,
                },
            )
            .unwrap(),
            CameraAnalyzer::try_circular(X, true).unwrap(),
            CameraAnalyzer::try_circular(X, false).unwrap(),
        ];
        for input in [
            state(0.0, 0.0, 0.0),
            state(1.0, 0.0, 0.0),
            state(-1.0, 0.0, 0.0),
            state(0.0, 1.0, 0.0),
            state(0.0, 0.0, 1.0),
            state(0.0, 0.0, -1.0),
            state(0.3, -0.2, 0.4),
        ] {
            let law = SourcePolarization::try_constant(X, input).unwrap();
            with_camera(cx, scalar.with_polarization(&law).unwrap(), 0.0, |camera| {
                near(
                    camera
                        .sample_ray(cx, ray(Z), 550.0, 0, settings(), 1, 1)
                        .unwrap()
                        .radiance_w_m2_sr_nm(),
                    emitted,
                );
                for analyzer in analyzers {
                    let a = analyzer.row();
                    let expected = emitted * (a.i + a.q * input.q + a.u * input.u + a.v * input.v);
                    let sample = camera
                        .sample_ray_analyzed(cx, ray(Z), 550.0, 0, analyzer, settings(), 1, 1)
                        .unwrap();
                    near(sample.radiance_w_m2_sr_nm(), expected);
                    if expected == 0.0 {
                        assert_eq!(sample.radiance_w_m2_sr_nm(), 0.0);
                    }
                    assert_eq!(sample.endpoint(), RadianceEndpoint::Environment);
                    assert_eq!(sample.transport_weight(), a.i);
                }
            });
        }
    });
}

#[test]
fn g1_directional_spectral_maps_share_original_texels_and_wavelengths() {
    with_cx(|cx, _| {
        let rotation = 0.37;
        let map = EnvironmentMap::try_from_linear_srgb(
            4,
            2,
            (0..8).map(|index| [1.0 + index as f32; 3]).collect(),
            rotation,
        )
        .unwrap();
        let scalar = EnvironmentRadiance::try_new(&map, 2.0).unwrap();
        let texels: Vec<_> = (0..8)
            .map(|index| {
                let q = if index % 2 == 0 { 1.0 } else { -1.0 };
                SourcePolarization::try_spectral(
                    Y,
                    vec![
                        (LAMBDA_MIN, state(q, 0.0, 0.0)),
                        (LAMBDA_MAX, state(0.0, 0.0, q)),
                    ],
                )
                .unwrap()
            })
            .collect();
        let source = scalar.with_polarization_map(cx, &texels).unwrap();
        assert!(core::ptr::eq(
            source.polarization_at_texel(6).unwrap(),
            &texels[6]
        ));
        assert!(source.polarization_at_texel(8).is_none());
        with_camera(cx, source, 0.0, |camera| {
            for index in 0..8 {
                let theta = core::f64::consts::PI * (0.5 + (index / 4) as f64) / 2.0;
                let phi = core::f64::consts::TAU * (0.5 + (index % 4) as f64) / 4.0 + rotation;
                let direction = Vec3::new(
                    theta.sin() * phi.cos(),
                    theta.cos(),
                    theta.sin() * phi.sin(),
                );
                assert_eq!(map.evaluate(direction).unwrap().texel_index, index);
                for wavelength in [LAMBDA_MIN, 0.5 * (LAMBDA_MIN + LAMBDA_MAX), LAMBDA_MAX] {
                    let input = texels[index].normalized_stokes(wavelength).unwrap();
                    let emitted = scalar.evaluate(direction, wavelength).unwrap();
                    for (analyzer, factor) in [
                        (
                            CameraAnalyzer::try_linear(Y, 0.0).unwrap(),
                            0.5 * (1.0 + input.q),
                        ),
                        (
                            CameraAnalyzer::try_circular(Y, true).unwrap(),
                            0.5 * (1.0 + input.v),
                        ),
                    ] {
                        let sample = camera
                            .sample_ray_analyzed(
                                cx,
                                ray(direction),
                                wavelength,
                                0,
                                analyzer,
                                settings(),
                                1,
                                1,
                            )
                            .unwrap();
                        near(sample.radiance_w_m2_sr_nm(), emitted * factor);
                    }
                }
            }
        });
    });
}

#[test]
fn g3_map_rotation_moves_the_reference_axis_and_circular_poles_are_regular() {
    with_cx(|cx, _| {
        let linear = SourcePolarization::try_constant(X, state(1.0, 0.0, 0.0)).unwrap();
        for rotation in [0.0, 0.13, 0.71, core::f64::consts::FRAC_PI_2] {
            let map = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], rotation).unwrap();
            let scalar = EnvironmentRadiance::try_new(&map, 1.0).unwrap();
            let emitted = scalar.evaluate(Y, 550.0).unwrap();
            with_camera(
                cx,
                scalar.with_polarization(&linear).unwrap(),
                0.0,
                |camera| {
                    let sample = camera
                        .sample_ray_analyzed(
                            cx,
                            ray(Y),
                            550.0,
                            0,
                            CameraAnalyzer::try_linear(X, 0.0).unwrap(),
                            settings(),
                            1,
                            1,
                        )
                        .unwrap();
                    near(
                        sample.radiance_w_m2_sr_nm(),
                        emitted * rotation.cos().powi(2),
                    );
                },
            );
        }
        let map = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap();
        let scalar = EnvironmentRadiance::try_new(&map, 1.0).unwrap();
        let circular = SourcePolarization::try_constant(Z, state(0.0, 0.0, 1.0)).unwrap();
        with_camera(
            cx,
            scalar.with_polarization(&circular).unwrap(),
            0.0,
            |camera| {
                let sample = camera
                    .sample_ray_analyzed(
                        cx,
                        ray(Z),
                        550.0,
                        0,
                        CameraAnalyzer::try_circular(X, true).unwrap(),
                        settings(),
                        1,
                        1,
                    )
                    .unwrap();
                near(
                    sample.radiance_w_m2_sr_nm(),
                    scalar.evaluate(Z, 550.0).unwrap(),
                );
            },
        );
    });
}

#[test]
fn g1_environment_observation_keeps_native_total_internal_reflection_phase() {
    with_cx(|cx, _| {
        // Plane x=z turns the reverse camera ray -X into -Z. At 45 degrees
        // and n=3/2 -> 1, rs=4/5-3i/5 and rp=7/25-24i/25. Independently,
        // rs*conj(rp)=4/5+3i/5: U'=0.8U-0.6V, V'=0.6U+0.8V.
        let geometry = InstanceScene::try_new(vec![triangle(
            1,
            [[-1.0, -1.0, -1.0], [-1.0, 1.0, -1.0], [1.0, 0.0, 1.0]],
        )])
        .unwrap();
        let glass = material(1.5, 0.0);
        let air = material(1.0, 0.0);
        let media = [
            SceneMedium::Homogeneous(&glass),
            SceneMedium::Homogeneous(&air),
        ];
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![SurfaceBinding {
                object_id: 1,
                interaction: Interaction::Dielectric {
                    outside: 0,
                    inside: 1,
                    surface: DielectricSurface::SMOOTH,
                },
            }],
            100_000,
        )
        .unwrap();
        let map = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap();
        let scalar = EnvironmentRadiance::try_new(&map, 5.0).unwrap();
        for input in [
            state(0.2, 0.4, 0.6),
            state(0.0, 0.0, 1.0),
            state(0.0, 1.0, 0.0),
        ] {
            let law = SourcePolarization::try_constant(Y, input).unwrap();
            let camera = CurvedCameraScene::try_new(
                cx,
                &scene,
                &[],
                Some(scalar.with_polarization(&law).unwrap()),
            )
            .unwrap();
            for (analyzer, fraction) in [
                (
                    CameraAnalyzer::try_linear(Y, core::f64::consts::FRAC_PI_4).unwrap(),
                    0.5 * (1.0 + 0.8 * input.u - 0.6 * input.v),
                ),
                (
                    CameraAnalyzer::try_circular(Y, true).unwrap(),
                    0.5 * (1.0 + 0.6 * input.u + 0.8 * input.v),
                ),
            ] {
                let sample = camera
                    .sample_ray_analyzed_with_samples(
                        cx,
                        Ray {
                            origin: Point3::new(0.4, 0.0, 0.0),
                            dir: X.scale(-1.0),
                        },
                        550.0,
                        0,
                        analyzer,
                        settings(),
                        |_| [0.5, 0.5, 0.0],
                    )
                    .unwrap();
                near(
                    sample.radiance_w_m2_sr_nm(),
                    scalar.evaluate(Z.scale(-1.0), 550.0).unwrap() * fraction,
                );
                assert_eq!(sample.path().surface_events(), 1);
                assert_eq!(sample.endpoint(), RadianceEndpoint::Environment);
            }
        }
    });
}

#[test]
fn g0_support_identity_limits_cancellation_and_infinite_absorption_remain_explicit() {
    with_cx(|cx, gate| {
        let map = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap();
        let scalar = EnvironmentRadiance::try_new(&map, 1.0).unwrap();
        let narrow = SourcePolarization::try_spectral(
            X,
            vec![(400.0, state(1.0, 0.0, 0.0)), (700.0, state(1.0, 0.0, 0.0))],
        )
        .unwrap();
        assert!(scalar.with_polarization(&narrow).is_err());
        assert!(scalar.with_polarization_map(cx, &[]).is_err());
        assert!(scalar.with_polarization_map(cx, &[narrow]).is_err());
        let unpolarized = SourcePolarization::try_constant(Z, state(0.0, 0.0, 0.0)).unwrap();
        let law = scalar.with_polarization(&unpolarized).unwrap();
        let mut reference = None;
        with_camera(cx, scalar, 0.0, |camera| {
            reference = Some(
                camera
                    .sample_ray(cx, ray(Z), 550.0, 0, settings(), 1, 1)
                    .unwrap(),
            );
        });
        with_camera(cx, law, 0.0, |camera| {
            assert_eq!(
                camera
                    .sample_ray(cx, ray(Z), 550.0, 0, settings(), 1, 1)
                    .unwrap(),
                reference.unwrap()
            );
            assert!(
                camera
                    .sample_ray(cx, ray(Z), 1000.0, 0, settings(), 1, 1)
                    .is_err()
            );
        });
        let singular = SourcePolarization::try_constant(Z, state(1.0, 0.0, 0.0)).unwrap();
        with_camera(
            cx,
            scalar.with_polarization(&singular).unwrap(),
            0.0,
            |camera| {
                assert!(
                    camera
                        .sample_ray(cx, ray(Z), 550.0, 0, settings(), 1, 1)
                        .is_err()
                );
            },
        );
        // Infinite true absorption is exact zero before source/frame lookup,
        // even when the authored linear reference would be singular.
        with_camera(
            cx,
            scalar.with_polarization(&singular).unwrap(),
            0.2,
            |camera| {
                assert_eq!(
                    camera
                        .sample_ray(cx, ray(Z), 550.0, 0, settings(), 1, 1)
                        .unwrap()
                        .radiance_w_m2_sr_nm(),
                    0.0
                );
            },
        );
        let texels = [unpolarized];
        gate.request();
        assert!(scalar.with_polarization_map(cx, &texels).is_err());
    });
}
