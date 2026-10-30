//! G1/G3/G4/G5 qualification of calibrated sources in the shared path tracer.
//! References use declared candela, metres and cd/m², with no display exposure.

#![cfg(feature = "photometry")]

#[path = "photometric_path/probe_throughput.rs"]
mod probe_throughput;

#[path = "photometric_path/illuminance.rs"]
mod illuminance;

#[path = "photometric_path/parallel.rs"]
mod parallel;

#[path = "photometric_path/area.rs"]
mod area;

#[path = "photometric_path/dielectric.rs"]
mod dielectric;

#[path = "photometric_path/brushed.rs"]
mod brushed;

use core::f64::consts::{FRAC_PI_2, PI};

use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_geom::{Point3, Vec3};
use fs_render::charts::TriMesh;
use fs_render::instances::RigidTransform;
use fs_render::lighting::EnvironmentMap;
use fs_render::photometry::path::{PathLuminanceEstimate, PathLuminanceJob, PathMedium};
use fs_render::photometry::source::{
    AngularDistribution, PhotometricSource, SPECTRAL_BINS, SourceSpectrum, wavelength_nm,
};
use fs_render::photometry::{PhotometryError, TransportOptions};
use fs_render::spectral::LiftedSpectrum;
use fs_render::tracer::{
    Camera, DirectStrategy, Material, Primitive, Sampler, Scene, Settings, Shape, TracerError,
};
use fs_render::volumes::VolumeGrid;
use fs_render::volumes::transport::scattering::{
    MassScatteringSpectrum, PhaseFunction, ScatteringGrid,
};
use fs_render::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};

fn with_cx<R>(operation: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 173,
                kernel_id: 0x5048_5054,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        operation(&cx, &gate)
    })
}

fn settings(depth: u32, spp: u32) -> Settings {
    Settings {
        width: 1,
        height: 1,
        spp,
        max_depth: depth,
        sampler: Sampler::OwenSobol,
        strategy: DirectStrategy::Mis,
        seed: 173,
    }
}

fn empty_scene(eye: [f64; 3], direction: [f64; 3]) -> Scene {
    Scene {
        primitives: Vec::new(),
        lights: Vec::new(),
        environment: None,
        camera: Camera {
            eye: Point3::new(eye[0], eye[1], eye[2]),
            forward: Vec3::new(direction[0], direction[1], direction[2]),
            up: Vec3::new(0.0, 1.0, 0.0),
            half_tan: 0.0,
        },
    }
}

fn floor(z: f64, half_width: f64, reflectance: LiftedSpectrum) -> Primitive {
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
        material: Material::Lambertian { reflectance },
        emission: None,
    }
}

fn receiver_scene() -> Scene {
    let mut scene = empty_scene([0.0, 0.0, 2.0], [0.0, 0.0, -1.0]);
    scene
        .primitives
        .push(floor(0.0, 10.0, LiftedSpectrum { c: [0.0; 3] }));
    scene
}

fn isotropic(position: [f64; 3], candela: f64, spectrum: SourceSpectrum) -> PhotometricSource {
    PhotometricSource::try_new(
        RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], position).unwrap(),
        4.0 * PI * candela,
        1.0,
        AngularDistribution::isotropic(),
        spectrum,
    )
    .unwrap()
}

fn medium_budget() -> AbsorbingGridBudget {
    AbsorbingGridBudget {
        max_grid_cells: 16,
        max_steps_per_segment: 64,
    }
}

fn render_pixel(
    scene: &Scene,
    sources: &[PhotometricSource],
    medium: PathMedium<'_, '_>,
    settings: Settings,
    cx: &Cx<'_>,
) -> PathLuminanceEstimate {
    let mut job = PathLuminanceJob::new(
        scene,
        sources,
        medium,
        TransportOptions::default(),
        settings,
        cx,
    )
    .unwrap();
    assert!(job.pixel_estimate(0).is_none());
    assert_eq!(job.advance(cx, usize::MAX).unwrap(), settings.spp as usize);
    assert!(job.is_complete());
    let estimate = job.pixel_estimate(0).unwrap();
    assert_eq!(estimate.samples, settings.spp);
    assert_eq!(
        estimate.xyz[1].to_bits(),
        estimate.luminance_cd_m2.to_bits()
    );
    assert!(estimate.xyz.iter().all(|value| value.is_finite()));
    assert!(estimate.total_error_bound_cd_m2.is_none());
    estimate
}

fn close(actual: f64, expected: f64, relative: f64, context: &str) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= relative * expected.abs().max(1.0e-12),
        "{context}: observed={actual:.12e}, expected={expected:.12e}, relative tolerance={relative}"
    );
}

#[test]
fn g1_shared_path_reports_absolute_lambertian_luminance() {
    with_cx(|cx, _| {
        let scene = receiver_scene();
        let sources = [isotropic(
            [0.0, 0.0, 3.0],
            90.0,
            SourceSpectrum::equal_energy(),
        )];
        let estimate = render_pixel(&scene, &sources, PathMedium::Vacuum, settings(1, 8192), cx);
        // 90 cd / (3 m)^2 = 10 lx; a rho=0.5 Lambertian returns 5/pi cd/m².
        close(
            estimate.luminance_cd_m2,
            5.0 / PI,
            0.002,
            "absolute Lambertian calibration",
        );
        assert!(
            estimate.estimated_standard_error_cd_m2.is_none(),
            "QMC does not supply an IID standard error"
        );
    });
}

#[test]
fn g1_shared_path_retains_monochromatic_fixture_and_material_spectra() {
    with_cx(|cx, _| {
        let mut scene = receiver_scene();
        scene.primitives[0].material = Material::Lambertian {
            reflectance: LiftedSpectrum {
                c: [0.0, 2.0, -0.5],
            },
        };
        let mut values = Vec::new();
        for bin in [12, 60] {
            let mut bins = [0.0; SPECTRAL_BINS];
            bins[bin] = 1.0;
            let sources = [isotropic(
                [0.0, 0.0, 3.0],
                90.0,
                SourceSpectrum::try_from_bins(bins).unwrap(),
            )];
            let estimate = render_pixel(
                &scene,
                &sources,
                PathMedium::Vacuum,
                // Increase narrow-band packet occupancy while retaining the
                // 0.2% physical tolerance; this is not a certified QMC bound.
                settings(1, 65_536),
                cx,
            );
            // Independent closed-form sigmoid reflectance at the admitted
            // spectral bin, rather than the production BRDF evaluator.
            let x = (wavelength_nm(bin) - 380.0) / 400.0;
            let polynomial = 2.0 * x - 0.5;
            let rho = 0.5 + 0.5 * polynomial / (1.0 + polynomial * polynomial).sqrt();
            close(
                estimate.luminance_cd_m2,
                10.0 * rho / PI,
                0.002,
                "monochromatic reflected luminance",
            );
            values.push(estimate.luminance_cd_m2);
        }
        assert!(
            values[1] > 1.8 * values[0],
            "equal lux must preserve the surface's spectral response"
        );
    });
}

#[test]
fn g1_absorbing_snapshot_integrates_both_legs_and_distinct_density_cells() {
    with_cx(|cx, _| {
        let scene = receiver_scene();
        let sources = [isotropic(
            [0.0, 0.0, 3.0],
            90.0,
            SourceSpectrum::equal_energy(),
        )];
        let uniform = [1.0];
        let layered = [1.0, 3.0];
        for (density, dimensions, cell, optical_depth) in [
            (&uniform[..], [1, 1, 1], [10.0, 10.0, 4.0], 1.0_f64),
            (&layered[..], [1, 1, 2], [10.0, 10.0, 2.0], 1.4_f64),
        ] {
            let grid = AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new(dimensions, density, [-5.0, -5.0, 0.0], cell),
                MassAbsorptionSpectrum::constant(0.2).unwrap(),
                None,
                medium_budget(),
            )
            .unwrap();
            let estimate = render_pixel(
                &scene,
                &sources,
                PathMedium::Absorbing(&grid),
                settings(1, 8192),
                cx,
            );
            // Uniform: (camera 2 + light 3)*0.2 = 1.0. Layered: camera
            // 2*1 plus light (2*1 + 1*3), all multiplied by 0.2 = 1.4.
            close(
                estimate.luminance_cd_m2,
                (5.0 / PI) * (-optical_depth).exp(),
                0.002,
                "physical density and both extinction legs",
            );
        }
    });
}

#[test]
fn g1_point_fixture_single_scattering_matches_collinear_slab_integral() {
    with_cx(|cx, _| {
        let scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        let sources = [isotropic(
            [0.0, 0.0, 2.0],
            100.0,
            SourceSpectrum::equal_energy(),
        )];
        let density = [1.0];
        let grid = ScatteringGrid::try_new(
            cx,
            AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &density, [-1.0, -1.0, 0.0], [2.0, 2.0, 1.0]),
                MassAbsorptionSpectrum::constant(0.2).unwrap(),
                None,
                medium_budget(),
            )
            .unwrap(),
            MassScatteringSpectrum::constant(0.3).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        let estimate = render_pixel(
            &scene,
            &sources,
            PathMedium::Scattering(&grid),
            settings(1, 32_768),
            cx,
        );
        // Camera/light medium lengths sum to 1 m for every collision z.
        // Integral_0^1 [I*sigma_s/(4*pi*(2-z)^2)] exp(-sigma_t) dz
        // = I*sigma_s*exp(-sigma_t)/(8*pi), with sigma_t=0.5/m.
        let expected = 100.0 * 0.3 * (-0.5_f64).exp() / (8.0 * PI);
        close(
            estimate.luminance_cd_m2,
            expected,
            0.025,
            "independent single-scattering integral",
        );
    });
}

#[test]
fn g1_fixture_light_reaches_a_shadowed_receiver_through_a_second_surface() {
    with_cx(|cx, _| {
        let mut scene = receiver_scene();
        scene.camera.eye.z = 0.5;
        scene
            .primitives
            .push(floor(2.0, 4.0, LiftedSpectrum { c: [0.0; 3] }));
        // The upward hemisphere illuminates the ceiling. It emits exactly
        // zero toward the floor, so a direct-only engine must remain black.
        let sources = [PhotometricSource::try_new(
            RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, 1.0]).unwrap(),
            2.0 * PI * 100.0,
            1.0,
            AngularDistribution::uniform_cone(FRAC_PI_2).unwrap(),
            SourceSpectrum::equal_energy(),
        )
        .unwrap()];
        let direct = render_pixel(&scene, &sources, PathMedium::Vacuum, settings(1, 256), cx);
        assert_eq!(direct.luminance_cd_m2, 0.0);
        let reflected = render_pixel(
            &scene,
            &sources,
            PathMedium::Vacuum,
            settings(2, 32_768),
            cx,
        );

        // Independent area quadrature over the actual square ceiling. The
        // two Lambertian factors are 0.5/pi; source gap is 1 m; floor/ceiling
        // separation is 2 m. This does not call any renderer transport helper.
        let subdivisions = 256;
        let dx = 8.0 / f64::from(subdivisions);
        let mut integral = 0.0;
        for ix in 0..subdivisions {
            for iy in 0..subdivisions {
                let x = -4.0 + (f64::from(ix) + 0.5) * dx;
                let y = -4.0 + (f64::from(iy) + 0.5) * dx;
                let radial_squared = x * x + y * y;
                let source_squared = radial_squared + 1.0;
                let connection_squared = radial_squared + 4.0;
                integral += dx * dx
                    / (source_squared
                        * source_squared.sqrt()
                        * connection_squared
                        * connection_squared);
            }
        }
        let expected = 0.5 * 0.5 * 100.0 * 4.0 * integral / (PI * PI);
        assert!(expected > 0.1);
        close(
            reflected.luminance_cd_m2,
            expected,
            0.04,
            "two-surface reflected fixture light",
        );
    });
}

fn estimate_bits(estimate: PathLuminanceEstimate) -> ([u64; 3], u32, Option<u64>) {
    (
        estimate.xyz.map(f64::to_bits),
        estimate.samples,
        estimate.estimated_standard_error_cd_m2.map(f64::to_bits),
    )
}

#[test]
fn g4_g5_progressive_paths_publish_whole_samples_and_resume_exactly() {
    with_cx(|cx, gate| {
        let scene = receiver_scene();
        let sources = [isotropic(
            [0.0, 0.0, 3.0],
            90.0,
            SourceSpectrum::equal_energy(),
        )];
        let mut settings = settings(1, 64);
        settings.width = 2;
        settings.sampler = Sampler::Iid;
        let mut reference = PathLuminanceJob::new(
            &scene,
            &sources,
            PathMedium::Vacuum,
            TransportOptions::default(),
            settings,
            cx,
        )
        .unwrap();
        assert_eq!(reference.advance(cx, usize::MAX).unwrap(), 128);
        let reference_bits: Vec<_> = (0..2)
            .map(|i| estimate_bits(reference.pixel_estimate(i).unwrap()))
            .collect();
        let mut job = PathLuminanceJob::new(
            &scene,
            &sources,
            PathMedium::Vacuum,
            TransportOptions::default(),
            settings,
            cx,
        )
        .unwrap();
        assert!(job.pixel_estimate(0).is_none());
        assert!(job.pixel_estimate(2).is_none());
        assert_eq!(job.advance(cx, 1).unwrap(), 1);
        let first = job.pixel_estimate(0).unwrap();
        assert_eq!(first.samples, 1);
        assert!(first.estimated_standard_error_cd_m2.is_none());
        assert!(job.pixel_estimate(1).is_none());
        assert_eq!(job.advance(cx, 2).unwrap(), 2);
        assert_eq!(job.completed_path_samples(), 3);
        let second = job.pixel_estimate(0).unwrap();
        assert_eq!(second.samples, 2);
        assert_eq!(job.pixel_estimate(1).unwrap().samples, 1);
        // For two IID observations, the estimated standard error equals
        // |x2-x1|/2 = |mean(x1,x2)-x1|, independently of the accumulator.
        close(
            second.estimated_standard_error_cd_m2.unwrap(),
            (second.luminance_cd_m2 - first.luminance_cd_m2).abs(),
            1.0e-12,
            "IID standard-error units",
        );
        let prefix: Vec<_> = (0..2)
            .map(|i| estimate_bits(job.pixel_estimate(i).unwrap()))
            .collect();
        gate.request();
        assert!(matches!(
            job.advance(cx, 8),
            Err(PhotometryError::Cancelled)
        ));
        assert_eq!(job.completed_path_samples(), 3);
        assert_eq!(
            (0..2)
                .map(|i| estimate_bits(job.pixel_estimate(i).unwrap()))
                .collect::<Vec<_>>(),
            prefix
        );
        with_cx(|resumed, _| {
            while !job.is_complete() {
                job.advance(resumed, 7).unwrap();
            }
            assert_eq!(job.advance(resumed, 1).unwrap(), 0);
        });
        assert_eq!(
            (0..2)
                .map(|i| estimate_bits(job.pixel_estimate(i).unwrap()))
                .collect::<Vec<_>>(),
            reference_bits
        );
    });
}

#[test]
fn g0_unqualified_sampler_and_uncalibrated_scene_lighting_refuse() {
    with_cx(|cx, _| {
        let mut scene = receiver_scene();
        let sources = [isotropic(
            [0.0, 0.0, 3.0],
            90.0,
            SourceSpectrum::equal_energy(),
        )];
        let mut settings = settings(1, 8);
        assert!(
            PathLuminanceJob::new(
                &scene,
                &sources,
                PathMedium::Vacuum,
                TransportOptions::default(),
                settings,
                cx
            )
            .is_ok()
        );
        settings.sampler = Sampler::OwenSobolFullPath;
        assert!(
            PathLuminanceJob::new(
                &scene,
                &sources,
                PathMedium::Vacuum,
                TransportOptions::default(),
                settings,
                cx
            )
            .is_err()
        );
        settings.sampler = Sampler::OwenSobol;
        scene.environment =
            Some(EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap());
        assert!(
            PathLuminanceJob::new(
                &scene,
                &sources,
                PathMedium::Vacuum,
                TransportOptions::default(),
                settings,
                cx
            )
            .is_err()
        );
        scene.environment = None;
        scene.primitives[0].emission = Some((LiftedSpectrum { c: [0.0; 3] }, 1.0));
        assert!(
            PathLuminanceJob::new(
                &scene,
                &sources,
                PathMedium::Vacuum,
                TransportOptions::default(),
                settings,
                cx
            )
            .is_err()
        );
    });
}

#[test]
fn g3_scattering_path_refuses_an_unsampled_fixture_singularity() {
    with_cx(|cx, _| {
        let scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        let sources = [isotropic(
            [0.0, 0.0, 0.3],
            100.0,
            SourceSpectrum::equal_energy(),
        )];
        let density = [1.0];
        let grid = ScatteringGrid::try_new(
            cx,
            AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &density, [-1.0, -1.0, 0.0], [2.0, 2.0, 1.0]),
                MassAbsorptionSpectrum::constant(0.0).unwrap(),
                None,
                medium_budget(),
            )
            .unwrap(),
            MassScatteringSpectrum::constant(0.1).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        let mut job = PathLuminanceJob::new(
            &scene,
            &sources,
            PathMedium::Scattering(&grid),
            TransportOptions::default(),
            settings(1, 8),
            cx,
        )
        .unwrap();
        assert!(matches!(
            job.advance(cx, 1),
            Err(PhotometryError::Geometry(TracerError::FixtureConnection {
                source: 0,
                ..
            }))
        ));
        assert_eq!(job.completed_path_samples(), 0);
        assert!(job.pixel_estimate(0).is_none());
    });
}

/// Execute explicitly with `--ignored --nocapture`. This is a fixed-ray,
/// one-surface spectral path workload, not a general image-frame benchmark.
#[test]
#[ignore = "host-specific calibrated path throughput qualification"]
fn photometric_path_300_fixture_throughput() {
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
        let mut settings = settings(1, 512);
        settings.width = 4;
        settings.height = 4;
        let mut job = PathLuminanceJob::new(
            &scene,
            &sources,
            PathMedium::Vacuum,
            TransportOptions::default(),
            settings,
            cx,
        )
        .unwrap();
        let start = std::time::Instant::now();
        let accepted = job.advance(cx, usize::MAX).unwrap();
        let elapsed = start.elapsed().as_secs_f64();
        assert!(job.is_complete());
        assert_eq!(accepted, 16 * 512);
        let expected_lux: f64 = (0..300)
            .map(|index| {
                let x = f64::from(index % 20) - 9.5;
                let y = f64::from(index / 20) - 7.0;
                let r_squared = x * x + y * y + 64.0;
                800.0 / (r_squared * r_squared.sqrt())
            })
            .sum();
        let expected_luminance = 0.5 * expected_lux / PI;
        for pixel in 0..16 {
            let estimate = job.pixel_estimate(pixel).unwrap();
            assert_eq!(estimate.samples, 512);
            close(
                estimate.luminance_cd_m2,
                expected_luminance,
                0.01,
                "300-fixture spectral path analytic calibration",
            );
        }
        println!(
            "photometric path: fixtures=300 pixels=16 spp=512 depth=1 paths={accepted} elapsed_s={elapsed:.6} paths_per_s={:.3} fixture_vertex_connections_per_s={:.3}",
            accepted as f64 / elapsed,
            (accepted * sources.len()) as f64 / elapsed,
        );
    });
}
