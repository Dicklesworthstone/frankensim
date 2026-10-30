//! Calibrated camera references: projection, pupil averages, medium and replay.

use super::*;
use crate::camera::{Aperture, CameraProjection, ExposureMetadata, PhysicalCamera};
use crate::charts::TriMesh;
use crate::instances::RigidTransform;
use crate::photometry::area::{
    EmissionSidedness, GoboTransmission, PhotometricAreaEmitter, ThinLensProjector,
};
use crate::photometry::source::{AngularDistribution, SourceSpectrum};
use crate::spectral::LiftedSpectrum;
use crate::tracer::{Camera, DirectStrategy, Material, Primitive, Shape};
use crate::volumes::VolumeGrid;
use crate::volumes::transport::scattering::{
    MassScatteringSpectrum, PhaseFunction, ScatteringGrid,
};
use crate::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};
use core::f64::consts::PI;
use fs_exec::{Budget, CancelGate, ExecMode, RunId, StreamKey};
use fs_geom::{Point3, Vec3};

fn scope<R>(f: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 173,
                kernel_id: 0x5048_4341,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}

fn scene() -> Scene {
    Scene {
        primitives: Vec::new(),
        lights: Vec::new(),
        environment: None,
        camera: Camera {
            eye: Point3::new(0.0, 0.0, 2.0),
            forward: Vec3::new(0.0, 0.0, -1.0),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.0,
        },
    }
}

fn panel(x: [f64; 2], y: [f64; 2], z: f64, front: bool) -> Primitive {
    Primitive {
        shape: Shape::Mesh(TriMesh::new(
            vec![
                [x[0], y[0], z],
                [x[1], y[0], z],
                [x[1], y[1], z],
                [x[0], y[1], z],
            ],
            if front {
                vec![[0, 1, 2], [0, 2, 3]]
            } else {
                vec![[0, 2, 1], [0, 3, 2]]
            },
        )),
        material: Material::Lambertian {
            reflectance: LiftedSpectrum { c: [0.0; 3] },
        },
        emission: None,
    }
}

fn camera(scene: &Scene, focus: f64, aperture: Aperture) -> PhysicalCamera {
    PhysicalCamera::try_legacy_compatible(
        scene.camera.eye,
        scene.camera.forward,
        scene.camera.up,
        scene.camera.half_tan,
        focus,
        aperture,
    )
    .unwrap()
}

fn emitter(index: usize, area: f64, luminance: f64) -> PhotometricAreaEmitter {
    PhotometricAreaEmitter::try_new(
        index,
        PI * area * luminance,
        1.0,
        SourceSpectrum::equal_energy(),
        EmissionSidedness::Front,
    )
    .unwrap()
}

fn settings(spp: u32) -> Settings {
    Settings {
        width: 1,
        height: 1,
        spp,
        max_depth: 1,
        sampler: Sampler::OwenSobol,
        strategy: DirectStrategy::BsdfOnly,
        seed: 173,
    }
}

fn job<'a, 'field>(
    scene: &'a Scene,
    camera: &'a PhysicalCamera,
    emitters: PhotometricEmitters<'a>,
    medium: PathMedium<'a, 'field>,
    settings: Settings,
    cx: &Cx<'_>,
) -> PathLuminanceJob<'a, 'field> {
    PathLuminanceJob::new_with_camera(
        scene,
        camera,
        emitters,
        medium,
        TransportOptions::default(),
        settings,
        cx,
    )
    .unwrap()
}

fn pixel(
    scene: &Scene,
    camera: &PhysicalCamera,
    emitters: &[PhotometricAreaEmitter],
    medium: PathMedium<'_, '_>,
    settings: Settings,
    cx: &Cx<'_>,
) -> PathLuminanceEstimate {
    let mut job = job(
        scene,
        camera,
        PhotometricEmitters {
            points: &[],
            surfaces: emitters,
        },
        medium,
        settings,
        cx,
    );
    job.advance(cx, usize::MAX).unwrap();
    assert!(job.is_complete());
    let value = job.pixel_estimate(0).unwrap();
    assert_eq!(value.samples, settings.spp);
    assert!(value.total_error_bound_cd_m2.is_none());
    value
}

type PixelBits = ([u64; 3], u32, Option<u64>, Option<u64>);
fn bits(job: &PathLuminanceJob<'_, '_>) -> Vec<Option<PixelBits>> {
    job.estimates()
        .map(|v| {
            v.map(|v| {
                (
                    v.xyz.map(f64::to_bits),
                    v.samples,
                    v.estimated_standard_error_cd_m2.map(f64::to_bits),
                    v.total_error_bound_cd_m2.map(f64::to_bits),
                )
            })
        })
        .collect()
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs(),
        "observed={actual:.12e}, reference={expected:.12e}, tolerance={relative}"
    );
}

#[test]
fn g5_physical_pinhole_preserves_legacy_bits_and_replaces_unused_camera() {
    scope(|cx, _| {
        let mut scene = scene();
        scene.camera.half_tan = 0.2;
        scene
            .primitives
            .push(panel([-2.0, 2.0], [-2.0, 2.0], 0.0, true));
        let surfaces = [emitter(0, 16.0, 25.0)];
        // Position-dependent reflected intensity exposes changed camera-ray
        // arithmetic; uniform emitted radiance alone would hide that defect.
        let points = [PhotometricSource::try_new(
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.3, 0.4, 1.3]).unwrap(),
            80.0 * PI,
            1.0,
            AngularDistribution::isotropic(),
            SourceSpectrum::equal_energy(),
        )
        .unwrap()];
        let camera = camera(&scene, 2.0, Aperture::try_circular(0.0).unwrap());
        for sampler in [Sampler::Iid, Sampler::OwenSobol] {
            let policy = Settings {
                width: 3,
                height: 2,
                sampler,
                strategy: DirectStrategy::Mis,
                ..settings(64)
            };
            let emitters = PhotometricEmitters {
                points: &points,
                surfaces: &surfaces,
            };
            let mut legacy = PathLuminanceJob::new_with_emitters(
                &scene,
                emitters,
                PathMedium::Vacuum,
                TransportOptions::default(),
                policy,
                cx,
            )
            .unwrap();
            let mut explicit = job(&scene, &camera, emitters, PathMedium::Vacuum, policy, cx);
            legacy.advance(cx, usize::MAX).unwrap();
            explicit.advance(cx, usize::MAX).unwrap();
            assert_eq!(bits(&explicit), bits(&legacy));
            assert!(
                explicit
                    .estimates()
                    .flatten()
                    .all(|v| v.luminance_cd_m2 > 0.0)
            );
        }
        let valid = pixel(
            &scene,
            &camera,
            &surfaces,
            PathMedium::Vacuum,
            settings(64),
            cx,
        );
        scene.camera.eye.x = f64::NAN;
        scene.camera.forward = Vec3::new(0.0, 0.0, 0.0);
        assert!(
            PathLuminanceJob::new_with_emitters(
                &scene,
                PhotometricEmitters {
                    points: &[],
                    surfaces: &surfaces
                },
                PathMedium::Vacuum,
                TransportOptions::default(),
                settings(64),
                cx,
            )
            .is_err()
        );
        let explicit = pixel(
            &scene,
            &camera,
            &surfaces,
            PathMedium::Vacuum,
            settings(64),
            cx,
        );
        assert_eq!(explicit.xyz.map(f64::to_bits), valid.xyz.map(f64::to_bits));
    });
}

#[test]
fn g1_physical_sensor_and_focal_length_set_calibrated_panel_coverage() {
    scope(|cx, _| {
        let mut scene = scene();
        scene
            .primitives
            .push(panel([-0.5, 0.5], [-0.25, 0.25], 0.0, true));
        let surfaces = [emitter(0, 0.5, 40.0)];
        for focal_length in [0.025, 0.050] {
            let camera = PhysicalCamera::try_look_at(
                scene.camera.eye,
                Point3::new(0.0, 0.0, 0.0),
                scene.camera.up,
                CameraProjection::try_focal_sensor(focal_length, 0.024).unwrap(),
                2.0,
                Aperture::try_circular(0.0).unwrap(),
            )
            .unwrap();
            // Uniform raster coordinates cover a square of half-width z*h/(2*f).
            // The panel occupies the product of its clipped X and Y fractions.
            let half_footprint = 2.0 * 0.024 / (2.0 * focal_length);
            let coverage =
                (0.5_f64 / half_footprint).min(1.0) * (0.25_f64 / half_footprint).min(1.0);
            let value = pixel(
                &scene,
                &camera,
                &surfaces,
                PathMedium::Vacuum,
                settings(32_768),
                cx,
            );
            close(value.luminance_cd_m2, 40.0 * coverage, 0.01);
            assert!(value.estimated_standard_error_cd_m2.is_none());
        }
    });
}

#[test]
fn g1_camera_defocus_matches_circular_segment_coverage_and_sharp_focus() {
    scope(|cx, _| {
        let mut scene = scene();
        scene
            .primitives
            .push(panel([0.1, 2.0], [-2.0, 2.0], 0.0, true));
        let surfaces = [emitter(0, 7.6, 25.0)];
        let policy = Settings {
            sampler: Sampler::Iid,
            ..settings(65_536)
        };
        for focus in [2.0, 4.0, 1.0] {
            let camera = camera(&scene, focus, Aperture::try_circular(0.4).unwrap());
            let value = pixel(&scene, &camera, &surfaces, PathMedium::Vacuum, policy, cx);
            // At depth d the centre pixel's circular footprint has radius
            // R*|1-d/focus|. The illuminated half-plane starts at x=0.1 m.
            let blur_radius = 0.4_f64 * (1.0_f64 - 2.0 / focus).abs();
            if blur_radius == 0.0 {
                assert_eq!(value.luminance_cd_m2, 0.0, "focus plane must be sharp");
            } else {
                let q = 0.1 / blur_radius;
                let fraction = (q.acos() - q * (1.0 - q * q).sqrt()) / PI;
                let expected = 25.0 * fraction;
                let tolerance = 5.0 * value.estimated_standard_error_cd_m2.unwrap() + 0.002 * 25.0;
                assert!(
                    (value.luminance_cd_m2 - expected).abs() < tolerance,
                    "focus={focus}, observed={}, reference={expected}, tolerance={tolerance}",
                    value.luminance_cd_m2
                );
            }
        }
    });
}

#[test]
fn g1_camera_pupil_and_exposure_do_not_scale_uniform_luminance() {
    scope(|cx, _| {
        let mut scene = scene();
        scene.camera.half_tan = 0.15;
        scene
            .primitives
            .push(panel([-10.0, 10.0], [-10.0, 10.0], 0.0, true));
        let surfaces = [emitter(0, 400.0, 25.0)];
        let mut reference = None;
        for aperture in [
            Aperture::try_circular(0.0).unwrap(),
            Aperture::try_circular(0.03).unwrap(),
            Aperture::try_circular(0.5).unwrap(),
            Aperture::try_regular_polygon(0.5, 5, 0.37).unwrap(),
        ] {
            let camera = camera(&scene, 3.0, aperture);
            let value = pixel(
                &scene,
                &camera,
                &surfaces,
                PathMedium::Vacuum,
                settings(8192),
                cx,
            );
            close(value.luminance_cd_m2, 25.0, 0.002);
            let expected = reference.get_or_insert(value.xyz.map(f64::to_bits));
            assert_eq!(value.xyz.map(f64::to_bits), *expected);
            let exposed =
                camera.with_exposure_metadata(ExposureMetadata::try_new(3200.0, 7.0).unwrap());
            let exposed_value = pixel(
                &scene,
                &exposed,
                &surfaces,
                PathMedium::Vacuum,
                settings(8192),
                cx,
            );
            assert_eq!(exposed_value.xyz.map(f64::to_bits), *expected);
        }
    });
}

#[test]
fn g1_camera_pupil_paths_match_homogeneous_beer_lambert_integral() {
    scope(|cx, _| {
        let mut scene = scene();
        scene
            .primitives
            .push(panel([-3.0, 3.0], [-3.0, 3.0], 0.0, true));
        let surfaces = [emitter(0, 36.0, 25.0)];
        let density = [1.0];
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [-3.0, -3.0, 0.0], [6.0, 6.0, 2.0]),
            MassAbsorptionSpectrum::constant(0.7).unwrap(),
            None,
            AbsorbingGridBudget {
                max_grid_cells: 1,
                max_steps_per_segment: 16,
            },
        )
        .unwrap();
        let camera = camera(&scene, 2.0, Aperture::try_circular(0.8).unwrap());
        let vacuum = pixel(
            &scene,
            &camera,
            &surfaces,
            PathMedium::Vacuum,
            settings(32_768),
            cx,
        );
        let attenuated = pixel(
            &scene,
            &camera,
            &surfaces,
            PathMedium::Absorbing(&medium),
            settings(32_768),
            cx,
        );
        // Squared pupil radius is uniform. Independently integrate the segment
        // length sqrt(d²+r²) over r² in [0,R²], using its elementary primitive.
        let sigma = 0.7_f64;
        let distance = 2.0_f64;
        let radius = 0.8_f64;
        let edge_distance = (distance * distance + radius * radius).sqrt();
        let reference = 2.0 / (radius * radius * sigma * sigma)
            * ((1.0 + sigma * distance) * (-sigma * distance).exp()
                - (1.0 + sigma * edge_distance) * (-sigma * edge_distance).exp());
        close(
            attenuated.luminance_cd_m2 / vacuum.luminance_cd_m2,
            reference,
            0.003,
        );
        assert!(
            reference < (-sigma * distance).exp(),
            "off-axis pupil rays travel farther"
        );
    });
}

#[test]
fn g4_g5_physical_camera_projector_paths_resume_with_identical_parallel_bits() {
    for sampler in [Sampler::Iid, Sampler::OwenSobol] {
        scope(|cx, gate| {
            let mut scene = scene();
            scene.camera.half_tan = 0.2;
            scene
                .primitives
                .push(panel([-3.0, 3.0], [-3.0, 3.0], 0.0, true));
            scene
                .primitives
                .push(panel([-0.05, 0.05], [-0.05, 0.05], 3.0, false));
            let frame = RigidTransform::try_new([1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 3.0]).unwrap();
            let projector = ThinLensProjector::try_new(
                frame,
                0.1,
                0.3 / 3.1,
                GoboTransmission::try_new(0.04, 0.04, 2, 2, vec![1.0, 0.25, 0.5, 0.0]).unwrap(),
            )
            .unwrap();
            let surfaces = [emitter(1, 0.01, 500.0).with_projector(projector).unwrap()];
            let points: Vec<_> = [-0.6, 0.6]
                .map(|x| {
                    PhotometricSource::try_new(
                        RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [x, 0.5, 1.5]).unwrap(),
                        80.0 * PI,
                        1.0,
                        AngularDistribution::isotropic(),
                        SourceSpectrum::equal_energy(),
                    )
                    .unwrap()
                })
                .into_iter()
                .collect();
            let camera = camera(
                &scene,
                2.5,
                Aperture::try_regular_polygon(0.1, 6, 0.2).unwrap(),
            );
            let policy = Settings {
                width: 3,
                height: 2,
                max_depth: 3,
                sampler,
                strategy: DirectStrategy::Mis,
                ..settings(32)
            };
            let emitters = PhotometricEmitters {
                points: &points,
                surfaces: &surfaces,
            };
            let sampling = PointFixtureSampling::Importance {
                connections_per_vertex: 1,
            };
            let mut serial = job(&scene, &camera, emitters, PathMedium::Vacuum, policy, cx)
                .with_point_sampling(sampling)
                .unwrap();
            let mut resumed = job(&scene, &camera, emitters, PathMedium::Vacuum, policy, cx)
                .with_point_sampling(sampling)
                .unwrap();
            serial.advance(cx, usize::MAX).unwrap();
            assert!(serial.estimates().flatten().all(|v| {
                v.luminance_cd_m2 > 0.0
                    && v.estimated_standard_error_cd_m2.is_some() == (sampler == Sampler::Iid)
                    && v.total_error_bound_cd_m2.is_none()
            }));
            let first_execution =
                RenderExecutionConfig::try_new(2, 2, 2, 64 << 20, RunId(901)).unwrap();
            resumed.advance(cx, 3).unwrap();
            resumed
                .advance_with_execution(cx, 7, &first_execution)
                .unwrap();
            let prefix = bits(&resumed);
            gate.request();
            assert!(matches!(
                resumed.advance_with_execution(cx, 20, &first_execution),
                Err(PhotometryError::Cancelled)
            ));
            assert_eq!(resumed.completed_path_samples(), 10);
            assert_eq!(bits(&resumed), prefix);
            scope(|fresh, _| {
                let execution =
                    RenderExecutionConfig::try_new(3, 1, 4, 64 << 20, RunId(902)).unwrap();
                while !resumed.is_complete() {
                    resumed
                        .advance_with_execution(fresh, 17, &execution)
                        .unwrap();
                    resumed.advance(fresh, 3).unwrap();
                }
            });
            assert_eq!(bits(&resumed), bits(&serial));
        });
    }
}

#[test]
#[ignore = "host-specific camera and scattering throughput qualification"]
fn photometric_physical_camera_throughput() {
    scope(|cx, _| {
        let mut scene = scene();
        scene.camera.half_tan = 0.2;
        scene
            .primitives
            .push(panel([-3.0, 3.0], [-3.0, 3.0], 0.0, true));
        scene
            .primitives
            .push(panel([-3.0, 3.0], [-3.0, 3.0], 3.0, false));
        let surfaces = [emitter(0, 36.0, 25.0)];
        let density = [0.2, 1.0, 0.6, 0.4, 1.4, 0.3, 0.7, 1.1];
        let haze = ScatteringGrid::try_new(
            cx,
            AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([2; 3], &density, [-3.0, -3.0, 0.0], [3.0, 3.0, 1.5]),
                MassAbsorptionSpectrum::constant(0.02).unwrap(),
                None,
                AbsorbingGridBudget {
                    max_grid_cells: 8,
                    max_steps_per_segment: 64,
                },
            )
            .unwrap(),
            MassScatteringSpectrum::constant(0.08).unwrap(),
            PhaseFunction::henyey_greenstein(0.35).unwrap(),
        )
        .unwrap();
        let policy = Settings {
            width: 16,
            height: 12,
            max_depth: 3,
            sampler: Sampler::Iid,
            strategy: DirectStrategy::Mis,
            ..settings(256)
        };
        for repeat in 0..3 {
            for radius in [0.0, 0.1] {
                let camera = camera(&scene, 2.5, Aperture::try_circular(radius).unwrap());
                let mut job = job(
                    &scene,
                    &camera,
                    PhotometricEmitters {
                        points: &[],
                        surfaces: &surfaces,
                    },
                    PathMedium::Scattering(&haze),
                    policy,
                    cx,
                );
                let start = std::time::Instant::now();
                job.advance(cx, usize::MAX).unwrap();
                let seconds = start.elapsed().as_secs_f64();
                let mean = job
                    .estimates()
                    .flatten()
                    .map(|v| v.luminance_cd_m2)
                    .sum::<f64>()
                    / 192.0;
                assert!(job.is_complete() && mean > 0.0 && mean.is_finite());
                println!(
                    "camera_throughput repeat={repeat} radius_m={radius} paths={} seconds={seconds:.6} paths_per_second={:.1} mean_cd_m2={mean:.9}",
                    job.completed_path_samples(),
                    job.completed_path_samples() as f64 / seconds
                );
            }
        }
    });
}
