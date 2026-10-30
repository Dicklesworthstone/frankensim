//! Exact native-traversal and whole-transport qualification of prepared geometry.

use super::*;
use crate::instances::{GeometryInstance, RigidTransform, SharedGeometry};
use crate::photometry::area::{EmissionSidedness, PhotometricAreaEmitter};
use crate::photometry::path::{IlluminanceSampling, PathIlluminanceJob, PathLuminanceJob};
use crate::photometry::source::{AngularDistribution, SourceSpectrum};
use fs_exec::{CancelGate, RunId, StreamKey};

fn scope<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 71,
                kernel_id: 99,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}

fn panel(center: [f64; 3], half: f64, id: Option<u64>) -> Primitive {
    let corners = [
        [-half, -half, 0.0],
        [half, -half, 0.0],
        [half, half, 0.0],
        [-half, half, 0.0],
    ];
    let shape = if let Some(id) = id {
        let transform =
            RigidTransform::try_new([0.0, 0.0, (PI / 12.0).sin(), (PI / 12.0).cos()], center)
                .unwrap();
        Shape::Instance(
            GeometryInstance::try_new(
                id,
                fs_blake3::hash_domain(
                    "org.frankensim.test.prepared-geometry",
                    &half.to_bits().to_le_bytes(),
                ),
                SharedGeometry::mesh(TriMesh::new(corners.to_vec(), vec![[0, 1, 2], [0, 2, 3]])),
                transform,
            )
            .unwrap(),
        )
    } else {
        Shape::Mesh(TriMesh::new(
            corners
                .map(|p| [p[0] + center[0], p[1] + center[1], center[2]])
                .to_vec(),
            vec![[0, 1, 2], [0, 2, 3]],
        ))
    };
    Primitive {
        shape,
        material: Material::Lambertian {
            reflectance: LiftedSpectrum { c: [0.0; 3] },
        },
        emission: None,
    }
}

fn decoys(primitives: &mut Vec<Primitive>) {
    // At least eight primitives engage the spatial index instead of the small
    // scene fallback. Their bounds are separate from the tested central rays.
    for i in 0..16 {
        primitives.push(panel(
            [20.0 + f64::from(i), 10.0, 3.0],
            0.3,
            Some(100 + i as u64),
        ));
    }
}

fn scene(primitives: Vec<Primitive>) -> Scene {
    Scene {
        primitives,
        lights: Vec::new(),
        environment: None,
        camera: Camera {
            eye: Point3::new(0.0, 0.0, 2.0),
            forward: Vec3::new(0.0, 0.0, -1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.2,
        },
    }
}

fn ray(origin: [f64; 3], direction: [f64; 3]) -> Ray {
    let d = Vec3::new(direction[0], direction[1], direction[2]);
    Ray {
        origin: Point3::new(origin[0], origin[1], origin[2]),
        dir: d.scale(1.0 / d.norm()),
    }
}

fn compare_hit(
    index: &geometry::PreparedGeometry,
    primitives: &[Primitive],
    cx: &Cx<'_>,
    ray: &Ray,
    limit: f64,
) -> Option<usize> {
    let accelerated = index.intersect_scene(primitives, cx, ray, limit).unwrap();
    let linear = intersect_primitives(primitives, cx, ray, None, Some(limit)).unwrap();
    match (accelerated, linear) {
        (None, None) => None,
        (Some(a), Some(b)) => {
            assert_eq!(a.primitive_index, b.primitive_index);
            assert_eq!(a.hit, b.hit, "native hit fields and work counts");
            assert_eq!(
                a.instance_hit, b.instance_hit,
                "instance identity and local witness"
            );
            assert_eq!(a.hit.t.to_bits(), b.hit.t.to_bits());
            assert_eq!(
                [a.hit.point.x, a.hit.point.y, a.hit.point.z].map(f64::to_bits),
                [b.hit.point.x, b.hit.point.y, b.hit.point.z].map(f64::to_bits)
            );
            Some(a.primitive_index)
        }
        _ => panic!("prepared and original traversal disagree on hit existence"),
    }
}

#[test]
fn g3_prepared_mixed_mesh_and_rigid_hits_equal_native_linear_queries() {
    scope(|cx, _| {
        let primitives: Vec<_> = (0..40)
            .map(|i| {
                panel(
                    [
                        2.0 * f64::from(i % 8) - 7.0,
                        2.0 * f64::from(i / 8) - 4.0,
                        3.0 + f64::from(i % 3),
                    ],
                    0.45,
                    (i % 2 == 0).then_some(i as u64 + 1),
                )
            })
            .collect();
        let index = geometry::PreparedGeometry::try_new(&primitives, cx).unwrap();
        assert!(index.is_accelerated_for_test());
        let mut hits = 0;
        let mut misses = 0;
        for x in -12..=12 {
            for y in -7..=7 {
                let ray = ray([0.0; 3], [0.1 * f64::from(x), 0.1 * f64::from(y), 1.0]);
                for limit in [2.0, 4.0, 10.0] {
                    if compare_hit(&index, &primitives, cx, &ray, limit).is_some() {
                        hits += 1;
                    } else {
                        misses += 1;
                    }
                }
            }
        }
        assert!(hits > 20 && misses > 20);
    });
}

#[test]
fn g3_prepared_bounds_preserve_mixed_ties_endpoints_edges_and_grazing() {
    scope(|cx, _| {
        for (order, expected) in [([Some(9), None, Some(3)], 2), ([None, Some(9), Some(3)], 0)] {
            let mut primitives: Vec<_> = order
                .into_iter()
                .map(|id| panel([0.0, 0.0, 2.0], 1.0, id))
                .collect();
            decoys(&mut primitives);
            let index = geometry::PreparedGeometry::try_new(&primitives, cx).unwrap();
            assert!(index.is_accelerated_for_test());
            let axial = ray([0.0; 3], [-0.0, 0.0, 1.0]);
            assert_eq!(
                compare_hit(&index, &primitives, cx, &axial, 2.0),
                Some(expected)
            );
            assert_eq!(
                compare_hit(&index, &primitives, cx, &axial, 2.0_f64.next_down()),
                None
            );
            // The native rigid leaf normalizes directions and returns world
            // distance. Public near-unit inputs must retain the original
            // distinction from mesh ray parameters at a finite endpoint.
            for magnitude in [1.0 - 1.0e-9, 1.0 + 1.0e-9] {
                let nonunit = Ray {
                    origin: Point3::new(0.0, 0.0, 0.0),
                    dir: Vec3::new(0.0, 0.0, magnitude),
                };
                for limit in [2.0_f64.next_down(), 2.0, 2.0_f64.next_up()] {
                    compare_hit(&index, &primitives, cx, &nonunit, limit);
                }
            }
            for x in [1.0_f64.next_down(), 1.0, 1.0_f64.next_up()] {
                compare_hit(&index, &primitives, cx, &ray([0.0; 3], [x, 0.0, 2.0]), 10.0);
            }
            for direction in [[1.0, 0.0, 0.0], [1.0, -0.0, -1.0e-14]] {
                compare_hit(
                    &index,
                    &primitives,
                    cx,
                    &ray([-10.0, 0.0, 2.0 + 1.0e-13], direction),
                    100.0,
                );
            }
        }
    });
}

#[test]
fn g3_prepared_bounds_keep_near_far_and_large_translated_hits() {
    scope(|cx, _| {
        for (origin, distance) in [
            ([0.0; 3], 1.0e-4),
            ([0.0; 3], 20_001.0),
            ([1.0e12, -1.0e12, 3.0e12], 2.0),
        ] {
            for id in [None, Some(7)] {
                let mut primitives =
                    vec![panel([origin[0], origin[1], origin[2] + distance], 1.0, id)];
                decoys(&mut primitives);
                let index = geometry::PreparedGeometry::try_new(&primitives, cx).unwrap();
                assert!(index.is_accelerated_for_test());
                assert_eq!(
                    compare_hit(
                        &index,
                        &primitives,
                        cx,
                        &ray(origin, [0.0, 0.0, 1.0]),
                        distance + 1.0
                    ),
                    Some(0)
                );
            }
        }
    });
}

fn source(position: [f64; 3], candela: f64) -> PhotometricSource {
    PhotometricSource::try_new(
        RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], position).unwrap(),
        4.0 * PI * candela,
        1.0,
        AngularDistribution::isotropic(),
        SourceSpectrum::equal_energy(),
    )
    .unwrap()
}

fn probe(x: f64, y: f64) -> IlluminanceProbe {
    IlluminanceProbe::try_new(Point3::new(x, y, 0.0), Vec3::new(0.0, 0.0, 1.0)).unwrap()
}

#[test]
fn g1_prepared_direct_lux_preserves_receiver_endpoint_and_real_shadow() {
    scope(|cx, _| {
        let mut primitives = vec![
            panel([0.0; 3], 8.0, None),
            panel([1.0, 0.0, 2.0], 0.25, Some(11)),
        ];
        decoys(&mut primitives);
        let sources = [
            source([0.0, 0.0, 4.0], 20.0),
            source([-2.0, 0.0, 4.0], 40.0),
            source([2.0, 0.0, 4.0], 80.0),
        ];
        let make = || {
            PhotometricScene::try_new(
                &primitives,
                &sources,
                HomogeneousMedium::vacuum(),
                TransportOptions::default(),
                cx,
            )
            .unwrap()
        };
        let indexed = make();
        let mut linear = make();
        linear.disable_acceleration_for_test();
        for x in [0.0, 0.5, 2.0] {
            let p = probe(x, 0.0);
            let a = indexed.illuminance(cx, &p).unwrap();
            let b = linear.illuminance(cx, &p).unwrap();
            assert_eq!(a.lux.to_bits(), b.lux.to_bits());
            assert_eq!(a, b);
            if x == 0.0 {
                let expected = 20.0 / 16.0 + 160.0 / (20.0 * 20.0_f64.sqrt());
                assert!(
                    (a.lux - expected).abs() < 1.0e-12 * expected,
                    "actual={} expected={expected}",
                    a.lux
                );
            }
        }
    });
}

fn image_bits(job: &PathLuminanceJob<'_, '_>) -> Vec<Option<([u64; 3], u32, Option<u64>)>> {
    job.estimates()
        .map(|v| {
            v.map(|v| {
                (
                    v.xyz.map(f64::to_bits),
                    v.samples,
                    v.estimated_standard_error_cd_m2.map(f64::to_bits),
                )
            })
        })
        .collect()
}

fn probe_bits(job: &PathIlluminanceJob<'_, '_>) -> Vec<Option<([u64; 3], u32, Option<u64>)>> {
    job.estimates()
        .map(|v| {
            v.map(|v| {
                (
                    [v.direct_lux, v.continuous_lux, v.illuminance_lux].map(f64::to_bits),
                    v.samples,
                    v.estimated_standard_error_lux.map(f64::to_bits),
                )
            })
        })
        .collect()
}

#[test]
fn g4_g5_prepared_primary_area_shadow_and_continuation_paths_resume_identically() {
    for strategy in [DirectStrategy::Mis, DirectStrategy::BsdfOnly] {
        scope(|cx, gate| {
            let mut primitives = vec![
                panel([0.0; 3], 6.0, None),
                panel([0.0, 0.0, 4.0], 1.5, Some(7)),
                panel([0.8, 0.0, 2.0], 0.25, None),
            ];
            decoys(&mut primitives);
            let scene = scene(primitives);
            let points = [source([-1.0, 1.0, 3.0], 50.0)];
            let surfaces = [PhotometricAreaEmitter::try_new(
                1,
                200.0,
                1.0,
                SourceSpectrum::equal_energy(),
                EmissionSidedness::Both,
            )
            .unwrap()];
            let emitters = PhotometricEmitters {
                points: if strategy == DirectStrategy::BsdfOnly {
                    &[]
                } else {
                    &points
                },
                surfaces: &surfaces,
            };
            let receivers = [probe(0.0, 0.0), probe(0.5, 0.25)];
            let settings = Settings {
                width: 2,
                height: 2,
                spp: 32,
                max_depth: 3,
                sampler: Sampler::Iid,
                strategy,
                seed: 71,
            };
            let image = || {
                PathLuminanceJob::new_with_emitters(
                    &scene,
                    emitters,
                    FixtureMedium::Vacuum,
                    TransportOptions::default(),
                    settings,
                    cx,
                )
                .unwrap()
            };
            let probes = || {
                PathIlluminanceJob::new_with_emitters(
                    &scene,
                    emitters,
                    FixtureMedium::Vacuum,
                    &receivers,
                    TransportOptions::default(),
                    IlluminanceSampling {
                        samples_per_probe: 32,
                        max_depth: 3,
                        sampler: Sampler::Iid,
                        seed: 71,
                    },
                    cx,
                )
                .unwrap()
            };
            let mut reference = image();
            let mut reference_probes = probes();
            reference.disable_acceleration_for_test();
            reference_probes.disable_acceleration_for_test();
            reference.advance(cx, usize::MAX).unwrap();
            reference_probes.advance(cx, usize::MAX).unwrap();
            assert!(
                reference
                    .estimates()
                    .flatten()
                    .any(|v| v.luminance_cd_m2 > 0.0)
            );
            assert!(
                reference_probes
                    .estimates()
                    .flatten()
                    .any(|v| v.continuous_lux > 0.0)
            );
            let exec = RenderExecutionConfig::try_new(2, 2, 4, 64 << 20, RunId(2701)).unwrap();
            let mut actual = image();
            let mut actual_probes = probes();
            actual.advance(cx, 3).unwrap();
            actual_probes.advance(cx, 2).unwrap();
            actual.advance_with_execution(cx, 5, &exec).unwrap();
            actual_probes.advance_with_execution(cx, 7, &exec).unwrap();
            let prefix = image_bits(&actual);
            let probe_prefix = probe_bits(&actual_probes);
            gate.request();
            assert!(matches!(
                actual.advance_with_execution(cx, 20, &exec),
                Err(PhotometryError::Cancelled)
            ));
            assert!(matches!(
                actual_probes.advance_with_execution(cx, 20, &exec),
                Err(PhotometryError::Cancelled)
            ));
            assert_eq!(image_bits(&actual), prefix);
            assert_eq!(probe_bits(&actual_probes), probe_prefix);
            scope(|resumed, _| {
                actual.advance(resumed, 11).unwrap();
                actual_probes.advance(resumed, 5).unwrap();
                actual
                    .advance_with_execution(resumed, usize::MAX, &exec)
                    .unwrap();
                actual_probes
                    .advance_with_execution(resumed, usize::MAX, &exec)
                    .unwrap();
            });
            assert_eq!(image_bits(&actual), image_bits(&reference));
            assert_eq!(probe_bits(&actual_probes), probe_bits(&reference_probes));
        });
    }
}

/// Original linear visibility and prepared visibility execute identical
/// calibrated queries. Geometry construction is outside the measured index
/// preparation; preparation and query times are reported separately.
#[test]
#[ignore = "host-specific dense 300-fixture prepared geometry qualification"]
fn photometric_prepared_geometry_300_fixture_throughput() {
    scope(|cx, _| {
        let mut primitives = vec![panel([0.0; 3], 20.0, None)];
        for i in 0..1024 {
            primitives.push(panel(
                [
                    3.0 * f64::from(i % 32) - 48.0,
                    3.0 * f64::from(i / 32) - 48.0,
                    2.0 + f64::from(i % 3),
                ],
                0.4,
                (i % 2 == 0).then_some(i as u64 + 1),
            ));
        }
        let sources: Vec<_> = (0..300)
            .map(|i| {
                source(
                    [f64::from(i % 20) - 9.5, f64::from(i / 20) - 7.0, 8.0],
                    100.0,
                )
            })
            .collect();
        let probes: Vec<_> = (0..16)
            .map(|i| probe(2.0 * f64::from(i % 4) - 3.0, 2.0 * f64::from(i / 4) - 3.0))
            .collect();
        let start = std::time::Instant::now();
        let indexed = PhotometricScene::try_new(
            &primitives,
            &sources,
            HomogeneousMedium::vacuum(),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        let build_seconds = start.elapsed().as_secs_f64();
        let mut linear = PhotometricScene::try_new(
            &primitives,
            &sources,
            HomogeneousMedium::vacuum(),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        linear.disable_acceleration_for_test();
        let mut expected = None;
        for (name, view) in [("linear", &linear), ("prepared", &indexed)] {
            let mut bits = Vec::with_capacity(128);
            let start = std::time::Instant::now();
            for _ in 0..8 {
                for probe in &probes {
                    bits.push(view.illuminance(cx, probe).unwrap().lux.to_bits());
                }
            }
            let seconds = start.elapsed().as_secs_f64();
            assert!(bits.iter().any(|&value| f64::from_bits(value) > 0.0));
            if let Some(expected) = &expected {
                assert_eq!(&bits, expected);
            } else {
                expected = Some(bits);
            }
            println!(
                "photometric_geometry traversal={name} primitives={} fixtures=300 probes=16 repeats=8 source_connections=38400 build_seconds={build_seconds:.6} query_seconds={seconds:.6} source_connections_per_second={:.3}",
                primitives.len(),
                38_400.0 / seconds
            );
        }
    });
}
