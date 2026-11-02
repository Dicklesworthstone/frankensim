use super::*;
use crate::camera::Aperture;
use crate::conductor::{
    CoherentDielectricCoating, ConductorDataStatus, ConductorIorSample, ConductorSource,
};
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, GlassProvenance};
use crate::motion::{ShotTimeBounds, ShutterConvention, ShutterDistribution};
use crate::polarization::Stokes;
use crate::polarization::dielectric::transport::{PolarizationFrame, PolarizedBeam};
use crate::spectral::lift_rgb;
use crate::surface_detail::{
    BrushingPattern, RadialMark, ShadingFrame, SurfaceDetail, shading_frame,
};
use fs_blake3::ContentHash;
use fs_exec::{CancelGate, StreamKey};

fn near(left: f64, right: f64) {
    assert!(
        (left - right).abs() <= 3.0e-9 * left.abs().max(right.abs()).max(1.0),
        "{left:.17e} != {right:.17e}"
    );
}

fn near_direction(left: Vec3, right: Vec3) {
    near(left.x, right.x);
    near(left.y, right.y);
    near(left.z, right.z);
}

fn optics() -> ConductorOptics {
    let source =
        ConductorSource::try_new(ContentHash([43; 32]), ConductorDataStatus::Representative)
            .unwrap();
    ConductorOptics::try_new(
        core::array::from_fn(|index| {
            ConductorIorSample::try_new(380.0 + 50.0 * index as f64, 0.2, 3.0).unwrap()
        }),
        source,
    )
    .unwrap()
    .with_coating(
        CoherentDielectricCoating::try_new(
            CauchyIor::try_constant(1.45).unwrap(),
            125.0e-9,
            source,
        )
        .unwrap(),
    )
}

fn detail(alpha: f64, ratio: f64, tilt: f64) -> SurfaceDetail {
    SurfaceDetail::try_new(
        alpha,
        ratio,
        BrushingPattern::Radial,
        tilt,
        Some(RadialMark {
            azimuth_rad: 0.0,
            half_width_rad: 0.1,
        }),
    )
    .unwrap()
}

fn quad(z: f64) -> TriMesh {
    TriMesh::new(
        vec![
            [-2.0, -2.0, z],
            [2.0, -2.0, z],
            [2.0, 2.0, z],
            [-2.0, 2.0, z],
        ],
        vec![[0, 1, 2], [0, 2, 3]],
    )
}

fn scene(detail: SurfaceDetail) -> Scene {
    Scene {
        primitives: vec![Primitive {
            shape: Shape::Mesh(quad(0.0)),
            material: Material::BrushedConductor {
                optics: optics(),
                detail,
            },
            emission: None,
        }],
        lights: Vec::new(),
        environment: None,
        camera: Camera {
            eye: Point3::new(0.4, 0.0, 1.0),
            forward: Vec3::new(0.0, 0.0, -1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 1.0,
        },
    }
}

fn vertex<'scene>(
    frame: ShadingFrame,
    wo: Vec3,
    medium: &'scene DielectricGlass,
) -> Vertex<'scene> {
    let host = Some(MediumEntry {
        boundary_primitive: 1,
        glass: medium,
    });
    Vertex {
        kind: VertexKind::Surface { primitive_index: 0 },
        geometric_normal: Some(frame.geometric_normal),
        brushed_frame: Some(frame),
        wo: Some(wo),
        positive_medium: host,
        negative_medium: host,
        ..Vertex::camera(Point3::new(0.4, 0.0, 0.0))
    }
}

fn with_test_cx<R>(operation: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 43,
                kernel_id: 1,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        operation(&cx)
    })
}

fn hit_vertex<'scene>(scene: &'scene Scene, cx: &Cx<'_>, ray: Ray) -> Vertex<'scene> {
    let intersection = intersect(scene, cx, &ray, None).unwrap().unwrap();
    let frame = surface_frame(&intersection.hit, &ray).unwrap();
    Vertex::surface(
        scene,
        intersection,
        ray,
        frame,
        &MediumStack::new(),
        1.0,
        1.0,
        None,
    )
    .unwrap()
}

#[test]
fn g1_brushed_bdpt_sample_pdf_and_adjoint_match_native_owners() {
    let medium = DielectricGlass::new(
        CauchyIor::try_constant(1.33).unwrap(),
        BeerLambertAbsorption::CLEAR,
        GlassProvenance::Custom,
    );
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let mut accepted = 0;
    let mut rejected = 0;
    let mut importance_above_one = false;
    for alpha in [0.05, 0.3] {
        for ratio in [1.0, 3.0] {
            for tilt in [0.0, 0.2, 0.35] {
                let detail = detail(alpha, ratio, tilt);
                let scene = scene(detail);
                let material = &scene.primitives[0].material;
                let frame = shading_frame(&detail, Point3::new(0.4, 0.0, 0.0), normal).unwrap();
                assert!(frame.mark_covered);
                for cosine in [0.005_f64, 0.4, 1.0] {
                    let wo = Vec3::new(0.0, (1.0 - cosine * cosine).sqrt(), cosine);
                    let vertex = vertex(frame, wo, &medium);
                    let beam = PolarizedBeam::try_new(
                        vertex.point,
                        PolarizationFrame::try_new(wo.scale(-1.0), frame.tangent).unwrap(),
                        550.0,
                        medium,
                        Stokes::unpolarized(1.0),
                    )
                    .unwrap();
                    for draw in [[0.0, 0.0], [0.13, 0.27], [0.6, 0.9], [0.99, 0.4]] {
                        let native =
                            bsdf_sample(material, normal, wo, draw[0], draw[1], Some(&frame));
                        let polarized = beam
                            .sample_brushed_conductor(&optics(), &detail, &frame, draw)
                            .unwrap();
                        for mode in [TransportMode::Radiance, TransportMode::Importance] {
                            let sampled = sample_surface(
                                &scene,
                                &vertex,
                                550.0,
                                mode,
                                [draw[0], draw[1], 0.4, 0.7],
                            )
                            .unwrap();
                            let Some((wi, pdf)) = native else {
                                assert!(sampled.is_none());
                                assert!(polarized.is_none());
                                rejected += 1;
                                continue;
                            };
                            let sampled = sampled.expect("native visible-normal sample survives");
                            accepted += 1;
                            assert_eq!(sampled.direction, wi);
                            assert_eq!(sampled.pdf_fwd_solid_angle, pdf);
                            near(
                                pdf,
                                surface_pdf_solid_angle(&scene, &vertex, wo, wi, 550.0).unwrap(),
                            );
                            near(
                                sampled.pdf_rev_solid_angle,
                                surface_pdf_solid_angle(&scene, &vertex, wi, wo, 550.0).unwrap(),
                            );
                            assert!(!sampled.delta);
                            assert_eq!(sampled.event, None);
                            let value =
                                surface_scattering(&scene, &vertex, wo, wi, 550.0, mode).unwrap();
                            near(sampled.weight, value * normal.dot(wi) / pdf);
                            if mode == TransportMode::Importance {
                                let adjoint = surface_scattering(
                                    &scene,
                                    &vertex,
                                    wi,
                                    wo,
                                    550.0,
                                    TransportMode::Radiance,
                                )
                                .unwrap();
                                assert_eq!(value.to_bits(), adjoint.to_bits());
                                // Independent forward-power owner includes the same coated
                                // Fresnel and host index, with the normal correction applied
                                // after Mueller reflection instead of by reversing this BSDF.
                                let polarized = polarized.as_ref().unwrap();
                                near(sampled.weight, polarized.beam.power().i);
                                near(pdf, polarized.pdf);
                                importance_above_one |= sampled.weight > 1.0;
                            } else {
                                let native_value = opaque_bsdf_eval(
                                    material,
                                    normal,
                                    wo,
                                    wi,
                                    550.0,
                                    Some(medium),
                                    Some(&frame),
                                )
                                .unwrap();
                                assert_eq!(value.to_bits(), native_value.to_bits());
                            }
                        }
                    }
                }
            }
        }
    }
    assert!(accepted > 0 && rejected > 0 && importance_above_one);
}

#[test]
fn g3_brushed_bdpt_hit_frame_and_both_transport_modes_rotate_with_the_body() {
    with_test_cx(|cx| {
        let detail = SurfaceDetail::try_new(
            0.1,
            4.0,
            BrushingPattern::Circular,
            0.15,
            Some(RadialMark {
                azimuth_rad: 0.0,
                half_width_rad: 0.1,
            }),
        )
        .unwrap();
        let angle = 0.7_f64;
        let placement = RigidTransform::try_new(
            [0.0, (angle / 2.0).sin(), 0.0, (angle / 2.0).cos()],
            [2.0, -3.0, 4.0],
        )
        .unwrap();
        let geometry = SharedGeometry::mesh(quad(0.0));
        let make_scene = |transform| {
            let mut scene = scene(detail);
            scene.primitives[0].shape = Shape::Instance(
                GeometryInstance::try_new(1, ContentHash([44; 32]), geometry.clone(), transform)
                    .unwrap(),
            );
            scene
        };
        let original_scene = make_scene(RigidTransform::identity());
        let moved_scene = make_scene(placement);
        for side in [1.0, -1.0] {
            let wo = Vec3::new(0.0, 0.6, side * 0.8);
            let ray = Ray {
                origin: Point3::new(0.4, 0.0, 0.0).offset(wo.scale(2.0)),
                dir: wo.scale(-1.0),
            };
            let moved_ray = Ray {
                origin: placement.transform_point(ray.origin),
                dir: placement.transform_vector(ray.dir),
            };
            let original = hit_vertex(&original_scene, cx, ray);
            let moved = hit_vertex(&moved_scene, cx, moved_ray);
            let frame = original.brushed_frame.unwrap();
            let moved_frame = moved.brushed_frame.unwrap();
            assert!(frame.mark_covered && moved_frame.mark_covered);
            for (a, b) in [
                (frame.geometric_normal, moved_frame.geometric_normal),
                (frame.shading_normal, moved_frame.shading_normal),
                (frame.tangent, moved_frame.tangent),
                (frame.bitangent, moved_frame.bitangent),
            ] {
                near_direction(placement.transform_vector(a), b);
            }
            near_direction(original.geometric_normal.unwrap(), Vec3::new(0.0, 0.0, 1.0));
            assert!(frame.geometric_normal.dot(wo) > 0.0);
            for mode in [TransportMode::Radiance, TransportMode::Importance] {
                let a = sample_surface(
                    &original_scene,
                    &original,
                    550.0,
                    mode,
                    [0.13, 0.27, 0.0, 0.0],
                )
                .unwrap()
                .unwrap();
                let b = sample_surface(&moved_scene, &moved, 550.0, mode, [0.13, 0.27, 0.0, 0.0])
                    .unwrap()
                    .unwrap();
                near_direction(placement.transform_vector(a.direction), b.direction);
                near(a.weight, b.weight);
                near(a.pdf_fwd_solid_angle, b.pdf_fwd_solid_angle);
                near(a.pdf_rev_solid_angle, b.pdf_rev_solid_angle);
            }
        }
    });
}

#[test]
fn g3_brushed_bdpt_undefined_radial_pole_retains_isotropic_reflection() {
    with_test_cx(|cx| {
        let mut scene = scene(detail(0.2, 3.0, 0.2));
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 1.0),
            dir: Vec3::new(0.0, 0.0, -1.0),
        };
        let at_pole = hit_vertex(&scene, cx, ray);
        assert!(at_pole.brushed_frame.is_none());
        let samples = [0.13, 0.27, 0.0, 0.0];
        let brushed = [TransportMode::Radiance, TransportMode::Importance].map(|mode| {
            sample_surface(&scene, &at_pole, 550.0, mode, samples)
                .unwrap()
                .unwrap()
        });
        // Drop the hit's scene borrow before installing the independently
        // admitted isotropic material at exactly the same surface.
        scene.primitives[0].material = Material::Conductor {
            optics: optics(),
            surface: ConductorSurface::try_rough(0.2).unwrap(),
        };
        let isotropic_vertex = hit_vertex(&scene, cx, ray);
        for (mode, actual) in [TransportMode::Radiance, TransportMode::Importance]
            .into_iter()
            .zip(brushed)
        {
            let expected = sample_surface(&scene, &isotropic_vertex, 550.0, mode, samples)
                .unwrap()
                .unwrap();
            assert_eq!(actual.direction, expected.direction);
            assert_eq!(actual.pdf_fwd_solid_angle, expected.pdf_fwd_solid_angle);
            assert_eq!(actual.pdf_rev_solid_angle, expected.pdf_rev_solid_angle);
            near(actual.weight, expected.weight);
        }
    });
}

#[test]
fn g0_brushed_bdpt_rejects_both_horizons_without_reorienting_the_finish() {
    let detail = detail(0.2, 3.0, 0.2);
    let scene = scene(detail);
    let normal = Vec3::new(0.0, 0.0, 1.0);
    let frame = shading_frame(&detail, Point3::new(0.4, 0.0, 0.0), normal).unwrap();
    let medium = DielectricGlass::new(
        CauchyIor::try_constant(1.0).unwrap(),
        BeerLambertAbsorption::CLEAR,
        GlassProvenance::Custom,
    );
    let shading_below = Vec3::new(0.0, -(1.0_f64 - 0.01 * 0.01).sqrt(), 0.01);
    let vertex = vertex(frame, normal, &medium);
    for rejected in [shading_below, normal.scale(-1.0)] {
        assert_eq!(
            surface_pdf_solid_angle(&scene, &vertex, normal, rejected, 550.0).unwrap(),
            0.0
        );
        for mode in [TransportMode::Radiance, TransportMode::Importance] {
            assert_eq!(
                surface_scattering(&scene, &vertex, normal, rejected, 550.0, mode).unwrap(),
                0.0
            );
        }
    }
    let below = Vertex {
        wo: Some(shading_below),
        ..vertex
    };
    for mode in [TransportMode::Radiance, TransportMode::Importance] {
        assert!(
            sample_surface(&scene, &below, 550.0, mode, [0.13; 4])
                .unwrap()
                .is_none()
        );
    }
    assert_eq!(below.brushed_frame, Some(frame));
}

#[test]
fn g1_brushed_bdpt_finite_emitter_reaches_camera_connections_and_light_splats() {
    with_test_cx(|cx| {
        let mut scene = scene(detail(0.2, 3.0, 0.2));
        let white = lift_rgb([1.0, 1.0, 1.0]);
        scene.primitives.push(Primitive {
            shape: Shape::Mesh(quad(2.0)),
            material: Material::Lambertian { reflectance: white },
            emission: Some((white, 2.0)),
        });
        scene.lights.push(RectLight {
            corner: Point3::new(-2.0, -2.0, 2.0),
            edge_u: Vec3::new(4.0, 0.0, 0.0),
            edge_v: Vec3::new(0.0, 4.0, 0.0),
            prim: 1,
            emission: (white, 2.0),
        });
        let physical = PhysicalCamera::try_legacy_compatible(
            scene.camera.eye,
            scene.camera.forward,
            scene.camera.up,
            scene.camera.half_tan,
            1.0,
            Aperture::try_circular(0.0).unwrap(),
        )
        .unwrap();
        let camera = AnimatedCamera::try_static(1, 0.0, 1.0, physical).unwrap();
        let shutter = ShutterInterval::resolve(
            0.5,
            0.0,
            ShutterConvention::Centered,
            ShutterDistribution::StratifiedCounterV1 { strata: 4 },
            ShotTimeBounds::try_new(0.0, 1.0).unwrap(),
        )
        .unwrap();
        let settings = Settings {
            width: 1,
            height: 1,
            spp: 64,
            max_depth: 1,
            sampler: Sampler::Iid,
            strategy: DirectStrategy::Mis,
            seed: 0x4252_5553,
        };
        for strategy_set in [
            BidirectionalStrategySet::Complete,
            BidirectionalStrategySet::CameraConnected,
        ] {
            let output = render_cinematic_bidirectional_with_strategy_set(
                &scene,
                &camera,
                CutSide::After,
                cx,
                &settings,
                shutter,
                strategy_set,
            )
            .unwrap();
            assert_eq!(output.film.spp_done, settings.spp);
            assert!(
                output.film.xyz[0]
                    .iter()
                    .all(|value| value.is_finite() && *value > 0.0)
            );
            assert!(output.strategies.nonzero > 0);
            if strategy_set == BidirectionalStrategySet::Complete {
                assert!(output.strategies.camera_splats > 0);
            } else {
                assert_eq!(output.strategies.camera_splats, 0);
            }
        }
    });
}
