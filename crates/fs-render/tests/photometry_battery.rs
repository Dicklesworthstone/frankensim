//! G0/G1/G3/G4/G5 qualification of calibrated point-source transport.
//! References here are closed-form photometry, independent of the source and
//! transport evaluators. Distances are metres and reported illuminance is lux.

#![cfg(feature = "photometry")]

use core::f64::consts::{FRAC_PI_2, PI};

use asupersync::types::Budget;
use fs_blake3::hash_domain;
use fs_evidence::NumericalKind;
use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::fixtures::SphereChart;
use fs_geom::{Point3, Vec3};
use fs_render::charts::TriMesh;
use fs_render::instances::{GeometryInstance, RigidTransform, SharedGeometry};
use fs_render::photometry::medium::{HomogeneousMedium, PhaseFunction};
use fs_render::photometry::source::{AngularDistribution, PhotometricSource, SourceSpectrum};
use fs_render::photometry::{
    CalculationGrid, IlluminanceJob, IlluminanceProbe, PhotometricScene, PhotometryError,
    TransportOptions,
};
use fs_render::spectral::LiftedSpectrum;
use fs_render::tracer::{Material, Primitive, Shape};

fn with_cx<R>(operation: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 0x5048_4f54,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        operation(&cx, &gate)
    })
}

fn placement(position: [f64; 3]) -> RigidTransform {
    RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], position).expect("finite translation")
}

fn isotropic(position: [f64; 3], intensity_cd: f64, dimmer: f64) -> PhotometricSource {
    PhotometricSource::try_new(
        placement(position),
        4.0 * PI * intensity_cd,
        dimmer,
        AngularDistribution::isotropic(),
        SourceSpectrum::equal_energy(),
    )
    .expect("isotropic fixture")
}

fn probe(position: [f64; 3], normal: [f64; 3]) -> IlluminanceProbe {
    IlluminanceProbe::try_new(
        Point3::new(position[0], position[1], position[2]),
        Vec3::new(normal[0], normal[1], normal[2]),
    )
    .expect("finite probe")
}

fn assert_close(observed: f64, expected: f64, tolerance: f64, context: &str) {
    let residual = (observed - expected).abs();
    assert!(
        observed.is_finite() && residual <= tolerance * expected.abs().max(1.0),
        "{context}: observed={observed:.17e} expected={expected:.17e} residual={residual:.6e}",
    );
}

fn square(z: f64, half_width: f64) -> Primitive {
    Primitive {
        shape: Shape::Mesh(TriMesh::new(
            vec![
                [-half_width, -half_width, z],
                [half_width, -half_width, z],
                [half_width, half_width, z],
                [-half_width, half_width, z],
            ],
            vec![[0, 1, 2], [0, 2, 3]],
        )),
        material: Material::Lambertian {
            reflectance: LiftedSpectrum { c: [0.0; 3] },
        },
        emission: None,
    }
}

#[test]
fn g1_inverse_square_and_receiver_cosine_laws_in_lux() {
    with_cx(|cx, _| {
        let sources = [isotropic([0.0, 0.0, 6.0], 120.0, 1.0)];
        let scene = PhotometricScene::try_new(
            &[],
            &sources,
            HomogeneousMedium::vacuum(),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        let near = scene
            .illuminance(cx, &probe([0.0, 0.0, 3.0], [0.0, 0.0, 1.0]))
            .unwrap();
        let far = scene
            .illuminance(cx, &probe([0.0; 3], [0.0, 0.0, 1.0]))
            .unwrap();
        let tilted = scene
            .illuminance(cx, &probe([0.0; 3], [3.0_f64.sqrt() * 0.5, 0.0, 0.5]))
            .unwrap();
        let rear = scene
            .illuminance(cx, &probe([0.0; 3], [0.0, 0.0, -1.0]))
            .unwrap();
        assert_close(near.lux, 120.0 / 9.0, 2.0e-13, "3 m illuminance");
        assert_close(far.lux, 120.0 / 36.0, 2.0e-13, "6 m illuminance");
        assert_close(near.lux / far.lux, 4.0, 2.0e-13, "doubling distance");
        assert_close(tilted.lux, 60.0 / 36.0, 2.0e-13, "60 degree receiver");
        assert_eq!(near.numerical.kind, NumericalKind::NoClaim);
        assert!(
            near.total_error_bound_lux.is_none(),
            "analytic agreement cannot certify unprovided fixture-calibration/model uncertainty"
        );
        assert_eq!(
            rear.lux, 0.0,
            "back-facing probe receives no direct illuminance"
        );
    });
}

#[test]
fn g1_flux_normalized_lambertian_source_has_both_source_and_receiver_cosines() {
    with_cx(|cx, _| {
        // I(theta)=100*cos(theta) cd integrates to 100*pi lm over a hemisphere.
        let sources = [PhotometricSource::try_new(
            RigidTransform::identity(),
            100.0 * PI,
            1.0,
            AngularDistribution::cosine_power(1.0, FRAC_PI_2).unwrap(),
            SourceSpectrum::equal_energy(),
        )
        .unwrap()];
        let scene = PhotometricScene::try_new(
            &[],
            &sources,
            HomogeneousMedium::vacuum(),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        let axial = scene
            .illuminance(cx, &probe([0.0, 0.0, 5.0], [0.0, 0.0, -1.0]))
            .unwrap();
        let off_axis = scene
            .illuminance(cx, &probe([3.0, 0.0, 4.0], [0.0, 0.0, -1.0]))
            .unwrap();
        let backward = scene
            .illuminance(cx, &probe([0.0, 0.0, -5.0], [0.0, 0.0, 1.0]))
            .unwrap();
        assert_close(
            axial.lux,
            100.0 / 25.0,
            2.0e-13,
            "normalized on-axis source",
        );
        assert_close(
            off_axis.lux,
            100.0 * 0.8 * 0.8 / 25.0,
            2.0e-13,
            "two independent cosines",
        );
        assert_eq!(backward.lux, 0.0);
    });
}

#[test]
fn g3_fixture_pan_tilt_and_common_rigid_motion_preserve_photometry() {
    with_cx(|cx, _| {
        let half = 0.5_f64.sqrt();
        let transform = RigidTransform::try_new([0.0, half, 0.0, half], [4.0, -3.0, 2.0]).unwrap();
        let distribution = AngularDistribution::uniform_cone(PI / 6.0).unwrap();
        let flux = 750.0;
        let original_sources = [PhotometricSource::try_new(
            RigidTransform::identity(),
            flux,
            0.7,
            distribution.clone(),
            SourceSpectrum::equal_energy(),
        )
        .unwrap()];
        let moved_sources = [PhotometricSource::try_new(
            transform,
            flux,
            0.7,
            distribution,
            SourceSpectrum::equal_energy(),
        )
        .unwrap()];
        let original = PhotometricScene::try_new(
            &[],
            &original_sources,
            HomogeneousMedium::vacuum(),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        let moved = PhotometricScene::try_new(
            &[],
            &moved_sources,
            HomogeneousMedium::vacuum(),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        for position in [[0.0, 0.0, 4.0], [0.8, 0.0, 4.0], [4.0, 0.0, 4.0]] {
            let before = probe(position, [0.0, 0.0, -1.0]);
            let after = IlluminanceProbe::try_new(
                transform.transform_point(Point3::new(position[0], position[1], position[2])),
                transform.transform_vector(Vec3::new(0.0, 0.0, -1.0)),
            )
            .unwrap();
            let a = original.illuminance(cx, &before).unwrap().lux;
            let b = moved.illuminance(cx, &after).unwrap().lux;
            assert_close(b, a, 1.0e-12, "common rigid transform");
        }
        let on_axis = moved
            .illuminance(cx, &probe([8.0, -3.0, 2.0], [-1.0, 0.0, 0.0]))
            .unwrap();
        let solid_angle = 2.0 * PI * (1.0 - (PI / 6.0).cos());
        assert_close(
            on_axis.lux,
            flux * 0.7 / (solid_angle * 16.0),
            1.0e-12,
            "moving-head 90 degree tilt",
        );
    });
}

#[test]
fn g1_extended_lambertian_disk_emitter_list_converges_to_analytic_illuminance() {
    with_cx(|cx, _| {
        // Uniform-radiance disk: E = pi * L * a² / (h² + a²).
        // Each point is a differential Lambertian patch with flux pi*L*dA.
        // Equal-area radial midpoint quadrature must converge at second order.
        let (radius, height, radiance) = (1.0_f64, 2.0_f64, 20.0_f64);
        let reference = PI * radiance * radius * radius / (height * height + radius * radius);
        let evaluate = |rings: u32| {
            let mut sources = Vec::new();
            let patch_flux = PI * radiance * (PI * radius * radius / f64::from(4 * rings));
            for ring in 0..rings {
                let r = radius * ((f64::from(ring) + 0.5) / f64::from(rings)).sqrt();
                for (x, y) in [(r, 0.0), (0.0, r), (-r, 0.0), (0.0, -r)] {
                    sources.push(
                        PhotometricSource::try_new(
                            placement([x, y, 0.0]),
                            patch_flux,
                            1.0,
                            AngularDistribution::cosine_power(1.0, FRAC_PI_2).unwrap(),
                            SourceSpectrum::equal_energy(),
                        )
                        .unwrap(),
                    );
                }
            }
            let scene = PhotometricScene::try_new(
                &[],
                &sources,
                HomogeneousMedium::vacuum(),
                TransportOptions::default(),
                cx,
            )
            .unwrap();
            scene
                .illuminance(cx, &probe([0.0, 0.0, height], [0.0, 0.0, -1.0]))
                .unwrap()
                .lux
        };
        let coarse_error = (evaluate(8) - reference).abs();
        let fine_error = (evaluate(32) - reference).abs();
        assert!((15.5..16.5).contains(&(coarse_error / fine_error)));
        assert_close(evaluate(128), reference, 1.0e-6, "uniform-radiance disk");
    });
}

#[test]
fn g3_emitters_and_dimmers_add_in_physical_units() {
    with_cx(|cx, _| {
        let sources = [
            isotropic([0.0, 0.0, 2.0], 100.0, 0.5),
            isotropic([0.0, 0.0, 4.0], 80.0, 0.25),
        ];
        let scene = PhotometricScene::try_new(
            &[],
            &sources,
            HomogeneousMedium::vacuum(),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        let result = scene
            .illuminance(cx, &probe([0.0; 3], [0.0, 0.0, 1.0]))
            .unwrap();
        assert_close(
            result.lux,
            50.0 / 4.0 + 20.0 / 16.0,
            2.0e-13,
            "two dimmed sources",
        );
    });
}

#[test]
fn g1_native_mesh_visibility_clips_only_between_receiver_and_fixture() {
    with_cx(|cx, _| {
        let sources = [isotropic([0.0, 0.0, 4.0], 160.0, 1.0)];
        let receiver = probe([0.0; 3], [0.0, 0.0, 1.0]);
        for (z, expected) in [(0.0, 10.0), (2.0, 0.0), (6.0, 10.0)] {
            let geometry = [square(z, 10.0)];
            let scene = PhotometricScene::try_new(
                &geometry,
                &sources,
                HomogeneousMedium::vacuum(),
                TransportOptions::default(),
                cx,
            )
            .unwrap();
            let result = scene.illuminance(cx, &receiver).unwrap();
            assert_close(
                result.lux,
                expected,
                2.0e-13,
                "receiver, blocker, or behind-source mesh",
            );
        }
    });
}

#[test]
fn g1_rigid_mesh_instances_preserve_inverse_square_and_segment_visibility() {
    with_cx(|cx, _| {
        let sources = [isotropic([0.0, 0.0, 4.0], 160.0, 1.0)];
        let receiver = probe([0.0; 3], [0.0, 0.0, 1.0]);
        let Shape::Mesh(mesh) = square(0.0, 10.0).shape else {
            unreachable!("square supplies native mesh geometry");
        };
        let shared = SharedGeometry::mesh(mesh);
        let identity = hash_domain("org.frankensim.test.photometry", b"placed-square");
        for (z, expected) in [(0.0, 10.0), (2.0, 0.0), (6.0, 10.0)] {
            let instance =
                GeometryInstance::try_new(1, identity, shared.clone(), placement([0.0, 0.0, z]))
                    .unwrap();
            let geometry = [Primitive {
                shape: Shape::Instance(instance),
                material: Material::Lambertian {
                    reflectance: LiftedSpectrum { c: [0.0; 3] },
                },
                emission: None,
            }];
            let scene = PhotometricScene::try_new(
                &geometry,
                &sources,
                HomogeneousMedium::vacuum(),
                TransportOptions::default(),
                cx,
            )
            .unwrap();
            assert_close(
                scene.illuminance(cx, &receiver).unwrap().lux,
                expected,
                2.0e-13,
                "placed receiver, blocker, or behind-source mesh",
            );
        }
    });
}

#[test]
fn g1_homogeneous_extinction_includes_absorption_and_out_scattering() {
    with_cx(|cx, _| {
        let sources = [isotropic([0.0, 0.0, 5.0], 250.0, 1.0)];
        let receiver = probe([0.0; 3], [0.0, 0.0, 1.0]);
        for (absorption, scattering) in [(0.2, 0.0), (0.0, 0.2), (0.07, 0.13)] {
            let medium = HomogeneousMedium::try_new(
                absorption,
                scattering,
                PhaseFunction::HenyeyGreenstein { g: 0.7 },
            )
            .unwrap();
            let scene =
                PhotometricScene::try_new(&[], &sources, medium, TransportOptions::default(), cx)
                    .unwrap();
            let result = scene.illuminance(cx, &receiver).unwrap();
            assert_close(
                result.lux,
                10.0 * (-1.0_f64).exp(),
                2.0e-13,
                "Beer-Lambert direct beam",
            );
        }
    });
}

#[test]
fn g0_singular_source_and_unresolved_distance_refuse() {
    with_cx(|cx, _| {
        let sources = [isotropic([0.0; 3], 100.0, 1.0)];
        let scene = PhotometricScene::try_new(
            &[],
            &sources,
            HomogeneousMedium::vacuum(),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        assert!(matches!(
            scene.illuminance(cx, &probe([0.0; 3], [0.0, 0.0, 1.0])),
            Err(PhotometryError::SourceTooClose { source: 0 })
        ));
        let options = TransportOptions {
            max_distance_m: 1.0,
            ..TransportOptions::default()
        };
        let bounded =
            PhotometricScene::try_new(&[], &sources, HomogeneousMedium::vacuum(), options, cx)
                .unwrap();
        assert!(
            matches!(
                bounded.illuminance(cx, &probe([0.0, 0.0, 2.0], [0.0, 0.0, -1.0])),
                Err(PhotometryError::DistanceLimitExceeded)
            ),
            "an unresolved visibility distance must not silently become zero or unoccluded light"
        );
    });
}

#[test]
fn g0_invalid_mesh_and_material_refuse_before_transport() {
    with_cx(|cx, _| {
        let sources = [isotropic([0.0, 0.0, 5.0], 100.0, 1.0)];
        let mut bad_mesh = square(0.0, 1.0);
        if let Shape::Mesh(mesh) = &mut bad_mesh.shape {
            mesh.vertices[0][0] = f64::NAN;
        }
        let geometry = [bad_mesh];
        assert!(
            PhotometricScene::try_new(
                &geometry,
                &sources,
                HomogeneousMedium::vacuum(),
                TransportOptions::default(),
                cx
            )
            .is_err()
        );
        let mut bad_material = square(0.0, 1.0);
        bad_material.material = Material::Lambertian {
            reflectance: LiftedSpectrum {
                c: [f64::NAN, 0.0, 0.0],
            },
        };
        let geometry = [bad_material];
        assert!(
            PhotometricScene::try_new(
                &geometry,
                &sources,
                HomogeneousMedium::vacuum(),
                TransportOptions::default(),
                cx
            )
            .is_err()
        );
    });
}

#[test]
fn g0_chart_geometry_refuses_before_uncertain_endpoint_visibility() {
    with_cx(|cx, _| {
        let geometry = [Primitive {
            shape: Shape::Chart(Box::new(SphereChart {
                center: Point3::new(0.0, 0.0, 0.0),
                radius: 1.0,
            })),
            material: Material::Lambertian {
                reflectance: LiftedSpectrum { c: [0.0; 3] },
            },
            emission: None,
        }];
        // Near a grazing receiver, chart proximity can be reached before the
        // excluded endpoint even with the default geometric tolerance.
        let sources = [isotropic([1.0, 0.0, 1.05], 100.0, 1.0)];
        for ray_epsilon_m in [TransportOptions::default().ray_epsilon_m, 1.0e-9] {
            let result = PhotometricScene::try_new(
                &geometry,
                &sources,
                HomogeneousMedium::vacuum(),
                TransportOptions {
                    ray_epsilon_m,
                    ..TransportOptions::default()
                },
                cx,
            );
            assert!(matches!(
                result,
                Err(PhotometryError::UnsupportedPrimitive {
                    primitive: 0,
                    reason,
                }) if reason.contains("chart")
            ));
        }
    });
}

#[test]
fn g0_rigid_chart_instances_refuse_before_uncertain_endpoint_visibility() {
    with_cx(|cx, _| {
        let instance = GeometryInstance::try_new(
            1,
            hash_domain("org.frankensim.test.photometry", b"placed-sphere"),
            SharedGeometry::chart(SphereChart {
                center: Point3::new(0.0, 0.0, 0.0),
                radius: 1.0,
            }),
            placement([2.0, 0.0, 0.0]),
        )
        .unwrap();
        let geometry = [Primitive {
            shape: Shape::Instance(instance),
            material: Material::Lambertian {
                reflectance: LiftedSpectrum { c: [0.0; 3] },
            },
            emission: None,
        }];
        let sources = [isotropic([3.0, 0.0, 1.05], 100.0, 1.0)];
        for ray_epsilon_m in [TransportOptions::default().ray_epsilon_m, 1.0e-9] {
            let result = PhotometricScene::try_new(
                &geometry,
                &sources,
                HomogeneousMedium::vacuum(),
                TransportOptions {
                    ray_epsilon_m,
                    ..TransportOptions::default()
                },
                cx,
            );
            assert!(matches!(
                result,
                Err(PhotometryError::UnsupportedPrimitive {
                    primitive: 0,
                    reason,
                }) if reason.contains("chart")
            ));
        }
    });
}

#[test]
fn g4_g5_progressive_grid_cancel_and_resume_preserve_accepted_results_bitwise() {
    with_cx(|cx, gate| {
        let sources = [
            isotropic([0.0, 0.0, 3.0], 100.0, 0.8),
            isotropic([2.0, -1.0, 4.0], 70.0, 1.0),
        ];
        let grid = CalculationGrid::try_new(
            Point3::new(-1.0, -1.0, 0.0),
            Vec3::new(0.5, 0.0, 0.0),
            Vec3::new(0.0, 0.5, 0.0),
            5,
            5,
        )
        .unwrap();
        let scene = PhotometricScene::try_new(
            &[],
            &sources,
            HomogeneousMedium::vacuum(),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        let reference = scene.illuminance_grid(&grid, cx).unwrap();
        let probes: Vec<_> = (0..25).map(|i| grid.probe(i).unwrap()).collect();
        let mut job = IlluminanceJob::new(&scene, &probes).unwrap();
        assert_eq!(job.advance(cx, 3).unwrap(), 3);
        let accepted: Vec<_> = job.results().iter().map(|e| e.lux.to_bits()).collect();
        gate.request();
        assert!(matches!(
            job.advance(cx, 5),
            Err(PhotometryError::Cancelled)
        ));
        assert_eq!(
            job.results()
                .iter()
                .map(|e| e.lux.to_bits())
                .collect::<Vec<_>>(),
            accepted
        );
        with_cx(|resume, _| {
            while !job.is_complete() {
                job.advance(resume, 7).unwrap();
            }
        });
        assert_eq!(job.results().len(), reference.len());
        for (index, (got, want)) in job.results().iter().zip(&reference).enumerate() {
            assert_eq!(
                got.lux.to_bits(),
                want.lux.to_bits(),
                "progressive probe={index}"
            );
            assert_eq!(
                got.numerical, want.numerical,
                "progressive uncertainty probe={index}"
            );
        }
    });
}

#[test]
fn g0_invalid_probe_and_degenerate_grid_refuse() {
    assert!(
        IlluminanceProbe::try_new(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 0.0)).is_err()
    );
    assert!(
        IlluminanceProbe::try_new(Point3::new(f64::NAN, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0))
            .is_err()
    );
    assert!(
        CalculationGrid::try_new(
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(2.0, 0.0, 0.0),
            4,
            4,
        )
        .is_err()
    );
}

/// Run with --release --test photometry_battery -- --ignored --nocapture.
/// Measures this host; it makes no claim about the requested 1440p image target.
#[test]
#[ignore = "release throughput qualification; reports host timing without a portable speed claim"]
fn photometry_300_fixture_grid_throughput() {
    with_cx(|cx, _| {
        let sources: Vec<_> = (0..300)
            .map(|i| {
                isotropic(
                    [(i % 20) as f64 - 9.5, (i / 20) as f64 - 7.0, 8.0],
                    100.0,
                    1.0,
                )
            })
            .collect();
        let geometry = [square(0.0, 20.0)];
        let grid = CalculationGrid::try_new(
            Point3::new(-6.3, -3.1, 0.0),
            Vec3::new(0.2, 0.0, 0.0),
            Vec3::new(0.0, 0.2, 0.0),
            64,
            32,
        )
        .unwrap();
        let scene = PhotometricScene::try_new(
            &geometry,
            &sources,
            HomogeneousMedium::vacuum(),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        let start = std::time::Instant::now();
        let results = scene.illuminance_grid(&grid, cx).unwrap();
        let elapsed = start.elapsed().as_secs_f64();
        for (index, observed) in results.iter().enumerate() {
            let x = -6.3 + 0.2 * (index % 64) as f64;
            let y = -3.1 + 0.2 * (index / 64) as f64;
            let expected: f64 = (0..300)
                .map(|i| {
                    let dx = (i % 20) as f64 - 9.5 - x;
                    let dy = (i / 20) as f64 - 7.0 - y;
                    let radius_squared = dx * dx + dy * dy + 64.0;
                    800.0 / (radius_squared * radius_squared.sqrt())
                })
                .sum();
            assert_close(observed.lux, expected, 2.0e-12, "300-source analytic sum");
        }
        println!(
            "photometry grid: fixtures=300 probes={} elapsed_s={elapsed:.6} probes_per_s={:.3} source_probe_pairs_per_s={:.3}",
            results.len(),
            results.len() as f64 / elapsed,
            (results.len() * sources.len()) as f64 / elapsed
        );
    });
}
