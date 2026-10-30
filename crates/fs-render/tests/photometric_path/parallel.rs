//! G3/G4/G5: calibrated native tiles preserve the original path estimator.

use super::*;
use fs_exec::RunId;
use fs_render::photometry::path::{IlluminanceSampling, PathIlluminanceJob};
use fs_render::photometry::{CalculationGrid, IlluminanceProbe};
use fs_render::tracer::{RenderExecutionConfig, RenderExecutionError};

type PixelBits = ([u64; 3], u64, u32, Option<u64>, Option<u64>);
type ProbeBits = ([u64; 3], u32, Option<u64>, Option<u64>);

fn image_bits(job: &PathLuminanceJob<'_, '_>) -> Vec<Option<PixelBits>> {
    job.estimates()
        .map(|value| {
            value.map(|v| {
                (
                    v.xyz.map(f64::to_bits),
                    v.luminance_cd_m2.to_bits(),
                    v.samples,
                    v.estimated_standard_error_cd_m2.map(f64::to_bits),
                    v.total_error_bound_cd_m2.map(f64::to_bits),
                )
            })
        })
        .collect()
}

fn probe_bits(job: &PathIlluminanceJob<'_, '_>) -> Vec<Option<ProbeBits>> {
    job.estimates()
        .map(|value| {
            value.map(|v| {
                (
                    [v.direct_lux, v.continuous_lux, v.illuminance_lux].map(f64::to_bits),
                    v.samples,
                    v.estimated_standard_error_lux.map(f64::to_bits),
                    v.total_error_bound_lux.map(f64::to_bits),
                )
            })
        })
        .collect()
}

fn execution(width: u32, height: u32, workers: usize, run: u64) -> RenderExecutionConfig {
    RenderExecutionConfig::try_new(width, height, workers, 64 << 20, RunId(run))
        .unwrap()
        .with_quantum_weights((0..workers).map(|index| (index % 3 + 1) as u32).collect())
        .unwrap()
}

fn room() -> Scene {
    let mut scene = receiver_scene();
    scene.camera.half_tan = 0.4;
    scene.primitives[0].material = Material::Lambertian {
        reflectance: LiftedSpectrum {
            c: [0.0, 0.6, -0.2],
        },
    };
    scene
        .primitives
        .push(floor(4.0, 6.0, LiftedSpectrum { c: [0.0; 3] }));
    scene
}

fn sources() -> [PhotometricSource; 2] {
    [
        isotropic([0.75, 0.25, 3.0], 90.0, SourceSpectrum::equal_energy()),
        isotropic([-1.0, -0.5, 2.5], 35.0, SourceSpectrum::equal_energy()),
    ]
}

fn haze<'a>(cx: &Cx<'_>, density: &'a [f64; 8]) -> ScatteringGrid<'a> {
    ScatteringGrid::try_new(
        cx,
        AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([2; 3], density, [-6.0, -6.0, 0.0], [6.0, 6.0, 2.0]),
            MassAbsorptionSpectrum::table(&[(380.0, 0.02), (780.0, 0.09)]).unwrap(),
            None,
            medium_budget(),
        )
        .unwrap(),
        MassScatteringSpectrum::constant(0.08).unwrap(),
        PhaseFunction::henyey_greenstein(0.35).unwrap(),
    )
    .unwrap()
}

fn image_settings(spp: u32) -> Settings {
    Settings {
        width: 4,
        height: 3,
        sampler: Sampler::Iid,
        ..settings(3, spp)
    }
}

fn sampling(spp: u32) -> IlluminanceSampling {
    IlluminanceSampling {
        samples_per_probe: spp,
        max_depth: 3,
        sampler: Sampler::Iid,
        seed: 173,
    }
}

fn grid() -> CalculationGrid {
    CalculationGrid::try_new(
        Point3::new(-0.6, -0.25, 0.0),
        Vec3::new(0.6, 0.0, 0.0),
        Vec3::new(0.0, 0.5, 0.0),
        3,
        2,
    )
    .unwrap()
}

fn image_job<'scene, 'field>(
    scene: &'scene Scene,
    sources: &'scene [PhotometricSource],
    medium: PathMedium<'scene, 'field>,
    settings: Settings,
    cx: &Cx<'_>,
) -> PathLuminanceJob<'scene, 'field> {
    PathLuminanceJob::new(
        scene,
        sources,
        medium,
        TransportOptions::default(),
        settings,
        cx,
    )
    .unwrap()
}

fn probe_job<'scene, 'field>(
    scene: &'scene Scene,
    sources: &'scene [PhotometricSource],
    medium: PathMedium<'scene, 'field>,
    sampling: IlluminanceSampling,
    cx: &Cx<'_>,
) -> PathIlluminanceJob<'scene, 'field> {
    PathIlluminanceJob::for_grid(
        scene,
        sources,
        medium,
        grid(),
        TransportOptions::default(),
        sampling,
        cx,
    )
    .unwrap()
}

#[test]
fn g5_parallel_luminance_matches_serial_xyz_and_iid_error_through_heterogeneous_haze() {
    with_cx(|cx, _| {
        let scene = room();
        let sources = sources();
        let density = [0.2, 1.0, 0.6, 0.4, 1.4, 0.3, 0.7, 1.1];
        let medium = haze(cx, &density);
        let settings = image_settings(16);
        let mut serial = image_job(
            &scene,
            &sources,
            PathMedium::Scattering(&medium),
            settings,
            cx,
        );
        serial.advance(cx, usize::MAX).unwrap();
        let expected = image_bits(&serial);
        assert!(serial.estimates().flatten().all(|value| {
            value.luminance_cd_m2 > 0.0
                && value.estimated_standard_error_cd_m2.is_some()
                && value.total_error_bound_cd_m2.is_none()
        }));

        // Twenty-four workers establish deterministic support, not a scaling
        // claim. Unequal scheduling weights and tile shapes change assignment.
        for (workers, width, height, chunks) in [
            (1, 2, 3, &[1, 13, 39, usize::MAX][..]),
            (2, 3, 2, &[5, 17, 47, usize::MAX][..]),
            (8, 1, 2, &[3, 19, 57, usize::MAX][..]),
            (24, 1, 1, &[usize::MAX][..]),
        ] {
            let policy = execution(width, height, workers, 100 + workers as u64);
            let mut job = image_job(
                &scene,
                &sources,
                PathMedium::Scattering(&medium),
                settings,
                cx,
            );
            for &chunk in chunks {
                let before = job.completed_path_samples();
                let accepted = job.advance_with_execution(cx, chunk, &policy).unwrap();
                assert_eq!(accepted, chunk.min((192 - before) as usize));
            }
            assert!(job.is_complete());
            assert_eq!(job.completed_path_samples(), 192);
            assert_eq!(image_bits(&job), expected, "workers={workers}");
        }
    });
}

#[test]
fn g5_parallel_probe_grid_matches_serial_lux_and_iid_error_through_heterogeneous_haze() {
    with_cx(|cx, _| {
        let scene = room();
        let sources = sources();
        let density = [0.2, 1.0, 0.6, 0.4, 1.4, 0.3, 0.7, 1.1];
        let medium = haze(cx, &density);
        let sampling = sampling(24);
        let mut serial = probe_job(
            &scene,
            &sources,
            PathMedium::Scattering(&medium),
            sampling,
            cx,
        );
        serial.advance(cx, usize::MAX).unwrap();
        let expected = probe_bits(&serial);
        assert!(serial.estimates().flatten().all(|value| {
            value.direct_lux > 0.0
                && value.continuous_lux > 0.0
                && value.estimated_standard_error_lux.is_some()
                && value.total_error_bound_lux.is_none()
        }));

        for (workers, width, height, chunks) in [
            (1, 3, 1, &[1, 11, 23, usize::MAX][..]),
            (2, 2, 3, &[5, 17, 29, usize::MAX][..]),
            (8, 2, 1, &[3, 19, 31, usize::MAX][..]),
            (24, 1, 1, &[usize::MAX][..]),
        ] {
            let policy = execution(width, height, workers, 200 + workers as u64);
            let mut job = probe_job(
                &scene,
                &sources,
                PathMedium::Scattering(&medium),
                sampling,
                cx,
            );
            for &chunk in chunks {
                let before = job.completed_path_samples();
                let accepted = job.advance_with_execution(cx, chunk, &policy).unwrap();
                assert_eq!(accepted, chunk.min((144 - before) as usize));
            }
            assert!(job.is_complete());
            assert_eq!(job.completed_path_samples(), 144);
            assert_eq!(probe_bits(&job), expected, "workers={workers}");
        }
    });
}

#[test]
fn g4_g5_cancelled_parallel_calls_preserve_prefix_and_mixed_resume_bits() {
    with_cx(|cx, gate| {
        let scene = room();
        let sources = sources();
        let mut image_reference =
            image_job(&scene, &sources, PathMedium::Vacuum, image_settings(9), cx);
        let mut probe_reference = probe_job(&scene, &sources, PathMedium::Vacuum, sampling(9), cx);
        image_reference.advance(cx, usize::MAX).unwrap();
        probe_reference.advance(cx, usize::MAX).unwrap();
        let mut image = image_job(&scene, &sources, PathMedium::Vacuum, image_settings(9), cx);
        let mut probes = probe_job(&scene, &sources, PathMedium::Vacuum, sampling(9), cx);
        let policy = execution(2, 2, 2, 301);
        image.advance(cx, 3).unwrap();
        probes.advance(cx, 2).unwrap();
        image.advance_with_execution(cx, 5, &policy).unwrap();
        probes.advance_with_execution(cx, 7, &policy).unwrap();
        let image_prefix = image_bits(&image);
        let probe_prefix = probe_bits(&probes);

        gate.request();
        assert!(matches!(
            image.advance_with_execution(cx, 20, &policy),
            Err(PhotometryError::Cancelled)
        ));
        assert!(matches!(
            probes.advance_with_execution(cx, 20, &policy),
            Err(PhotometryError::Cancelled)
        ));
        assert_eq!(image.completed_path_samples(), 8);
        assert_eq!(probes.completed_path_samples(), 9);
        assert_eq!(image_bits(&image), image_prefix);
        assert_eq!(probe_bits(&probes), probe_prefix);

        with_cx(|resumed, _| {
            let changed_policy = execution(3, 1, 4, 302);
            while !image.is_complete() {
                image
                    .advance_with_execution(resumed, 17, &changed_policy)
                    .unwrap();
                image.advance(resumed, 3).unwrap();
            }
            while !probes.is_complete() {
                probes.advance(resumed, 5).unwrap();
                probes
                    .advance_with_execution(resumed, 11, &changed_policy)
                    .unwrap();
            }
        });
        assert_eq!(image_bits(&image), image_bits(&image_reference));
        assert_eq!(probe_bits(&probes), probe_bits(&probe_reference));
    });
}

#[test]
fn g3_parallel_physical_range_refusal_discards_the_incomplete_wave() {
    with_cx(|cx, _| {
        let mut scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        scene
            .primitives
            .push(floor(20.0, 1.0e6, LiftedSpectrum { c: [0.0; 3] }));
        let sources = [];
        let probes = [-1.0, -1.0, 1.0].map(|z| {
            IlluminanceProbe::try_new(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, z)).unwrap()
        });
        let options = TransportOptions {
            max_distance_m: 10.0,
            ..TransportOptions::default()
        };
        let mut job = PathIlluminanceJob::new(
            &scene,
            &sources,
            PathMedium::Vacuum,
            &probes,
            options,
            sampling(4),
            cx,
        )
        .unwrap();
        job.advance(cx, 1).unwrap();
        let prefix = probe_bits(&job);
        // One safe downward path and one out-of-range upward path share this
        // two-address wave. Neither may publish when the latter refuses.
        let policy = execution(2, 1, 2, 401);
        assert!(matches!(
            job.advance_with_execution(cx, 2, &policy),
            Err(PhotometryError::Geometry(TracerError::FixturePathRange))
        ));
        assert_eq!(job.completed_path_samples(), 1);
        assert_eq!(probe_bits(&job), prefix);

        with_cx(|resumed, _| {
            // The image entry uses the same refusal, including an empty
            // accepted prefix; a failed wave must not become a black image.
            let mut image = PathLuminanceJob::new(
                &scene,
                &sources,
                PathMedium::Vacuum,
                options,
                image_settings(4),
                resumed,
            )
            .unwrap();
            assert!(matches!(
                image.advance_with_execution(resumed, 4, &policy),
                Err(PhotometryError::Geometry(TracerError::FixturePathRange))
            ));
            assert_eq!(image.completed_path_samples(), 0);
            assert!(image.estimates().all(|value| value.is_none()));
        });
    });
}

#[test]
fn g3_parallel_memory_refusal_preserves_both_jobs_and_allows_retry() {
    with_cx(|cx, _| {
        let scene = room();
        let sources = sources();
        let mut image = image_job(&scene, &sources, PathMedium::Vacuum, image_settings(3), cx);
        let mut probes = probe_job(&scene, &sources, PathMedium::Vacuum, sampling(3), cx);
        image.advance(cx, 2).unwrap();
        probes.advance(cx, 2).unwrap();
        let image_prefix = image_bits(&image);
        let probe_prefix = probe_bits(&probes);
        let refused = RenderExecutionConfig::try_new(2, 2, 2, 1, RunId(501)).unwrap();
        assert!(matches!(
            image.advance_with_execution(cx, 12, &refused),
            Err(PhotometryError::Execution(RenderExecutionError::Memory(_)))
        ));
        assert!(matches!(
            probes.advance_with_execution(cx, 12, &refused),
            Err(PhotometryError::Execution(RenderExecutionError::Memory(_)))
        ));
        assert_eq!(image.completed_path_samples(), 2);
        assert_eq!(probes.completed_path_samples(), 2);
        assert_eq!(image_bits(&image), image_prefix);
        assert_eq!(probe_bits(&probes), probe_prefix);

        let admitted = execution(2, 2, 2, 502);
        image
            .advance_with_execution(cx, usize::MAX, &admitted)
            .unwrap();
        probes
            .advance_with_execution(cx, usize::MAX, &admitted)
            .unwrap();
        let mut image_reference =
            image_job(&scene, &sources, PathMedium::Vacuum, image_settings(3), cx);
        let mut probe_reference = probe_job(&scene, &sources, PathMedium::Vacuum, sampling(3), cx);
        image_reference.advance(cx, usize::MAX).unwrap();
        probe_reference.advance(cx, usize::MAX).unwrap();
        assert_eq!(image_bits(&image), image_bits(&image_reference));
        assert_eq!(probe_bits(&probes), probe_bits(&probe_reference));
    });
}

/// Fixed analytic receiver workloads; worker timings include dispatch and
/// collection. These are not general image-frame or 24-core scaling claims.
#[test]
#[ignore = "host-specific serial versus native-tile photometric throughput"]
fn photometric_parallel_300_fixture_throughput() {
    with_cx(|cx, _| {
        let scene = receiver_scene();
        let sources: Vec<_> = (0..300)
            .map(|index| {
                isotropic(
                    [
                        f64::from(index % 20) - 9.5,
                        f64::from(index / 20) - 7.0,
                        8.0,
                    ],
                    100.0,
                    SourceSpectrum::equal_energy(),
                )
            })
            .collect();
        let settings = Settings {
            width: 4,
            height: 4,
            ..settings(1, 512)
        };
        let expected_lux: f64 = (0..300)
            .map(|index| {
                let x = f64::from(index % 20) - 9.5;
                let y = f64::from(index / 20) - 7.0;
                let r2 = x * x + y * y + 64.0;
                800.0 / (r2 * r2.sqrt())
            })
            .sum();
        let mut expected_bits = None;
        for workers in [0, 1, 2, 4, 8] {
            let policy = RenderExecutionConfig::try_new(
                32,
                32,
                workers.max(1),
                64 << 20,
                RunId(600 + workers as u64),
            )
            .unwrap();
            let mut job = image_job(&scene, &sources, PathMedium::Vacuum, settings, cx);
            let start = std::time::Instant::now();
            let accepted = if workers == 0 {
                job.advance(cx, usize::MAX)
            } else {
                job.advance_with_execution(cx, usize::MAX, &policy)
            }
            .unwrap();
            let seconds = start.elapsed().as_secs_f64();
            assert!(job.is_complete());
            assert_eq!(accepted, 8192);
            for value in job.estimates().flatten() {
                close(
                    value.luminance_cd_m2,
                    0.5 * expected_lux / PI,
                    0.01,
                    "parallel image analytic calibration",
                );
            }
            let bits = image_bits(&job);
            if let Some(reference) = &expected_bits {
                assert_eq!(&bits, reference);
            } else {
                expected_bits = Some(bits);
            }
            println!(
                "photometric parallel image: workers={workers} fixtures=300 paths={accepted} elapsed_s={seconds:.6} paths_per_s={:.3}",
                accepted as f64 / seconds
            );
        }

        let mut scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        scene
            .primitives
            .push(floor(2.0, 100.0, LiftedSpectrum { c: [0.0; 3] }));
        let sources: Vec<_> = (0..300)
            .map(|_| isotropic([0.0, 0.0, 1.0], 1.0, SourceSpectrum::equal_energy()))
            .collect();
        let probes =
            [
                IlluminanceProbe::try_new(Point3::new(0.0, 0.0, 0.0), Vec3::new(0.0, 0.0, 1.0))
                    .unwrap(),
            ];
        let sampling = IlluminanceSampling {
            samples_per_probe: 8192,
            max_depth: 1,
            sampler: Sampler::OwenSobol,
            seed: 173,
        };
        let expected_reflected = 150.0 * (1.0 - 4.0 * PI / (9.0 * 3.0_f64.sqrt()));
        let mut expected_bits = None;
        for workers in [0, 1, 2, 4, 8] {
            let policy = RenderExecutionConfig::try_new(
                32,
                32,
                workers.max(1),
                64 << 20,
                RunId(700 + workers as u64),
            )
            .unwrap();
            let mut job = PathIlluminanceJob::new(
                &scene,
                &sources,
                PathMedium::Vacuum,
                &probes,
                TransportOptions::default(),
                sampling,
                cx,
            )
            .unwrap();
            let start = std::time::Instant::now();
            let accepted = if workers == 0 {
                job.advance(cx, usize::MAX)
            } else {
                job.advance_with_execution(cx, usize::MAX, &policy)
            }
            .unwrap();
            let seconds = start.elapsed().as_secs_f64();
            assert!(job.is_complete());
            assert_eq!(accepted, 8192);
            let value = job.probe_estimate(0).unwrap();
            close(
                value.direct_lux,
                300.0,
                1.0e-12,
                "parallel probe exact direct illuminance",
            );
            close(
                value.continuous_lux,
                expected_reflected,
                0.02,
                "parallel probe reflected illuminance",
            );
            let bits = probe_bits(&job);
            if let Some(reference) = &expected_bits {
                assert_eq!(&bits, reference);
            } else {
                expected_bits = Some(bits);
            }
            println!(
                "photometric parallel probe: workers={workers} fixtures=300 paths={accepted} elapsed_s={seconds:.6} paths_per_s={:.3} approximate_fixture_vertex_connections_per_s={:.3}",
                accepted as f64 / seconds,
                (accepted * sources.len()) as f64 / seconds
            );
        }
    });
}
