//! Manufactured angular laws through gray, spectral and finite-source haze.
//! Independent integrals qualify the represented single-scattering model;
//! neither input measurements nor interpolation/model error are certified.

use super::*;
use fs_exec::RunId;
use fs_render::charts::Ray;
use fs_render::photometry::area::{EmissionSidedness, PhotometricAreaEmitter, PhotometricEmitters};
use fs_render::photometry::medium::{HomogeneousMedium, ScatterSettings};
use fs_render::photometry::path::{IlluminanceSampling, PathIlluminanceJob};
use fs_render::photometry::{IlluminanceProbe, PhotometricScene};
use fs_render::tracer::RenderExecutionConfig;

const HAZE_FIXTURES: [([f64; 3], f64); 2] = [([0.8, 0.0, 0.8], 100.0), ([-1.2, 0.0, 0.2], 35.0)];

fn linear_phase(a: f64) -> PhaseFunction {
    // Deliberately nonuniform knots; linear interpolation is exactly the
    // manufactured p(mu)=(1+a*mu)/(4*pi), whose mean cosine is a/3.
    PhaseFunction::tabulated(&[-1.0, -0.6, 0.05, 0.7, 1.0].map(|mu| (mu, 1.0 + a * mu))).unwrap()
}

fn slab<'a>(
    cx: &Cx<'_>,
    density: &'a [f64; 2],
    absorption: MassAbsorptionSpectrum,
    scattering: MassScatteringSpectrum,
    phase: PhaseFunction,
) -> ScatteringGrid<'a> {
    ScatteringGrid::try_new(
        cx,
        AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 2], density, [-3.0, -3.0, 0.0], [6.0, 6.0, 0.5]),
            absorption,
            None,
            medium_budget(),
        )
        .unwrap(),
        scattering,
        phase,
    )
    .unwrap()
}

fn sources(spectrum: SourceSpectrum) -> Vec<PhotometricSource> {
    HAZE_FIXTURES
        .iter()
        .map(|&(position, candela)| isotropic(position, candela, spectrum.clone()))
        .collect()
}

fn iid_settings(depth: u32, spp: u32) -> Settings {
    Settings {
        sampler: Sampler::Iid,
        ..settings(depth, spp)
    }
}

fn statistical_match(value: PathLuminanceEstimate, expected: f64, oracle_error: f64) {
    let se = value
        .estimated_standard_error_cd_m2
        .expect("IID sampling standard error");
    // A fixed regression threshold using estimated variance; this is not a
    // confidence interval or a total physical/numerical error certificate.
    let tolerance = 6.0 * se + oracle_error + 2.0e-10 * expected.abs().max(1.0);
    assert!(se.is_finite() && se >= 0.0);
    assert!(
        (value.luminance_cd_m2 - expected).abs() <= tolerance,
        "observed={} expected={expected} estimated_se={se} tolerance={tolerance}",
        value.luminance_cd_m2
    );
    assert!(value.total_error_bound_cd_m2.is_none());
    assert_eq!(value.numerical.kind, fs_evidence::NumericalKind::NoClaim);
}

// Independent Simpson integration split at the actual density interface.
// For the straight source segment, ds/dz=r/(z_source-z). Its column mass is
// therefore r*|C(z_source)-C(z)|/|z_source-z|, with C the elementary piecewise
// linear density primitive. This does not call production DDA, phase, Beer,
// source-spectrum or observer functions.
fn haze_reference(density: [f64; 2], a: f64, absorption: f64, scattering: f64, panels: u32) -> f64 {
    let column = |z: f64| density[0] * z.min(0.5) + density[1] * (z - 0.5).max(0.0);
    let mut total = 0.0;
    for (layer, rho) in density.into_iter().enumerate() {
        let low = layer as f64 * 0.5;
        let step = 0.5 / f64::from(panels);
        let integrand = |z: f64| {
            HAZE_FIXTURES
                .iter()
                .map(|&(position, candela)| {
                    let dz = position[2] - z;
                    let r2 = position[0].powi(2) + position[1].powi(2) + dz * dz;
                    let r = r2.sqrt();
                    let source_column = if dz == 0.0 {
                        rho * r
                    } else {
                        r * ((column(position[2]) - column(z)) / dz)
                    };
                    let phase = (1.0 + a * dz / r) / (4.0 * PI);
                    candela / r2
                        * rho
                        * scattering
                        * phase
                        * (-(absorption + scattering) * (column(z) + source_column)).exp()
                })
                .sum::<f64>()
        };
        let mut sum = integrand(low) + integrand(low + 0.5);
        for panel in 1..panels {
            sum +=
                if panel % 2 == 0 { 2.0 } else { 4.0 } * integrand(low + f64::from(panel) * step);
        }
        total += sum * step / 3.0;
    }
    total
}

#[test]
fn g1_tabulated_linear_phase_matches_gray_and_native_haze_for_both_signs() {
    with_cx(|cx, _| {
        let sources = sources(SourceSpectrum::equal_energy());
        let scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        let density = [1.0; 2];
        let mut values = Vec::new();
        for a in [-0.75, 0.75] {
            let phase = linear_phase(a);
            let gray = HomogeneousMedium::try_new(0.1, 0.2, phase.clone()).unwrap();
            let gray_scene =
                PhotometricScene::try_new(&[], &sources, gray, TransportOptions::default(), cx)
                    .unwrap();
            let gray_value = gray_scene
                .luminance_along_ray(
                    cx,
                    &Ray {
                        origin: Point3::new(0.0, 0.0, 0.0),
                        dir: Vec3::new(0.0, 0.0, 1.0),
                    },
                    1.0,
                    ScatterSettings {
                        initial_panels: 64,
                        max_panels: 8192,
                        absolute_tolerance_cd_m2: 0.0,
                        relative_tolerance: 1.0e-8,
                    },
                )
                .unwrap();
            let expected = haze_reference(density, a, 0.1, 0.2, 4096);
            close(
                haze_reference(density, a, 0.1, 0.2, 2048),
                expected,
                2.0e-10,
                "independent layer-Simpson refinement",
            );
            close(
                gray_value.luminance_cd_m2,
                expected,
                4.0e-8,
                "gray table phase and propagation sign",
            );
            assert_eq!(gray_value.surface_cd_m2, 0.0);
            assert!(gray_value.single_scatter.tolerance_met);
            assert!(gray_value.total_error_bound_cd_m2.is_none());
            assert!(
                gray_value
                    .single_scatter
                    .omitted_multiple_scattering_bound_cd_m2
                    .is_none()
            );
            let medium = slab(
                cx,
                &density,
                MassAbsorptionSpectrum::constant(0.1).unwrap(),
                MassScatteringSpectrum::constant(0.2).unwrap(),
                phase,
            );
            let value = render_pixel(
                &scene,
                &sources,
                PathMedium::Scattering(&medium),
                iid_settings(1, 32_768),
                cx,
            );
            statistical_match(value, expected, 2.0e-10 * expected);
            values.push(value.luminance_cd_m2);
        }
        // Reversing the same angular law changes this asymmetric rig's
        // radiance. A silently substituted isotropic law cannot pass both.
        assert!(
            values[1] > 1.1 * values[0],
            "forward/backward values: {values:?}"
        );
    });
}

#[test]
fn g1_tabulated_two_layer_spectral_haze_matches_independent_cie_xyz() {
    with_cx(|cx, _| {
        let scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        let density = [0.4, 1.6];
        let medium = slab(
            cx,
            &density,
            MassAbsorptionSpectrum::table(&[(380.0, 0.04), (780.0, 0.3)]).unwrap(),
            MassScatteringSpectrum::table(&[(380.0, 0.25), (780.0, 0.07)]).unwrap(),
            linear_phase(-0.6),
        );
        let mut measured = Vec::new();
        for (bin, wavelength) in [(13, 447.5), (54, 652.5)] {
            let mut bins = [0.0; SPECTRAL_BINS];
            bins[bin] = 1.0;
            let sources = sources(SourceSpectrum::try_from_bins(bins).unwrap());
            let t = (wavelength - 380.0) / 400.0;
            let absorption = 0.04 + 0.26 * t;
            let scattering = 0.25 - 0.18 * t;
            let expected = haze_reference(density, -0.6, absorption, scattering, 4096);
            close(
                haze_reference(density, -0.6, absorption, scattering, 2048),
                expected,
                2.0e-8,
                "independent heterogeneous spectral Beer legs",
            );
            let value = render_pixel(
                &scene,
                &sources,
                PathMedium::Scattering(&medium),
                iid_settings(1, 65_536),
                cx,
            );
            statistical_match(value, expected, 2.0e-8 * expected);
            // For a monochromatic represented bin, all three channels have
            // exactly the original CIE table's X:Y:Z ratios. No production
            // observer helper supplies this oracle.
            let cie = super::spectral_reference::cie_reference_xyz(wavelength);
            for (actual, component) in value.xyz.into_iter().zip(cie) {
                close(
                    actual,
                    value.luminance_cd_m2 * component / cie[1],
                    3.0e-12,
                    "independent monochromatic CIE ratio",
                );
            }
            measured.push(value.luminance_cd_m2);
        }
        assert!(
            measured[0] > 1.3 * measured[1],
            "spectral haze: {measured:?}"
        );
    });
}

// The unit-area panel spans x=[1.5,2.5], y=[-0.5,0.5] at z=2 m and
// emits L=100 cd/m² downward. The camera misses it; its only one-vertex
// contribution is L*cos_source/r² * phase * sigma_s with both Beer legs.
fn finite_panel_reference(panels: u32) -> f64 {
    let step = 1.0 / f64::from(panels);
    let mut sum = 0.0;
    for iz in 0..panels {
        let z = (f64::from(iz) + 0.5) * step;
        let dz = 2.0 - z;
        for iy in 0..panels {
            let y = -0.5 + (f64::from(iy) + 0.5) * step;
            for ix in 0..panels {
                let x = 1.5 + (f64::from(ix) + 0.5) * step;
                let r2 = x * x + y * y + dz * dz;
                let r = r2.sqrt();
                let mu = dz / r;
                let light_distance_in_slab = r * (1.0 - z) / dz;
                sum += 100.0 * mu / r2 * (1.0 + 0.6 * mu) / (4.0 * PI)
                    * 0.2
                    * (-0.3 * (z + light_distance_in_slab)).exp();
            }
        }
    }
    sum * step.powi(3)
}

#[test]
fn g1_tabulated_finite_emitter_haze_matches_independent_solid_angle_integral() {
    with_cx(|cx, _| {
        let mut scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        scene.primitives.push(Primitive {
            shape: Shape::Mesh(TriMesh::new(
                vec![
                    [1.5, -0.5, 2.0],
                    [2.5, -0.5, 2.0],
                    [2.5, 0.5, 2.0],
                    [1.5, 0.5, 2.0],
                ],
                vec![[0, 2, 1], [0, 3, 2]],
            )),
            material: Material::Lambertian {
                reflectance: LiftedSpectrum { c: [0.0; 3] },
            },
            emission: None,
        });
        let surfaces = [PhotometricAreaEmitter::try_new(
            0,
            100.0 * PI,
            1.0,
            SourceSpectrum::equal_energy(),
            EmissionSidedness::Front,
        )
        .unwrap()];
        let density = [1.0; 2];
        let medium = slab(
            cx,
            &density,
            MassAbsorptionSpectrum::constant(0.1).unwrap(),
            MassScatteringSpectrum::constant(0.2).unwrap(),
            linear_phase(0.6),
        );
        let expected = finite_panel_reference(64);
        close(
            finite_panel_reference(32),
            expected,
            0.0002,
            "independent panel solid-angle quadrature refinement",
        );
        for strategy in [DirectStrategy::NeeOnly, DirectStrategy::Mis] {
            let mut job = PathLuminanceJob::new_with_emitters(
                &scene,
                PhotometricEmitters {
                    points: &[],
                    surfaces: &surfaces,
                },
                PathMedium::Scattering(&medium),
                TransportOptions::default(),
                Settings {
                    strategy,
                    ..iid_settings(1, 32_768)
                },
                cx,
            )
            .unwrap();
            job.advance(cx, usize::MAX).unwrap();
            statistical_match(job.pixel_estimate(0).unwrap(), expected, 0.0002 * expected);
        }
    });
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

fn room_haze<'a>(cx: &Cx<'_>, density: &'a [f64; 8], phase: PhaseFunction) -> ScatteringGrid<'a> {
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
        MassScatteringSpectrum::constant(0.25).unwrap(),
        phase,
    )
    .unwrap()
}

fn room_sources() -> [PhotometricSource; 2] {
    [
        isotropic([0.75, 0.25, 3.0], 90.0, SourceSpectrum::equal_energy()),
        isotropic([-1.0, -0.5, 2.5], 35.0, SourceSpectrum::equal_energy()),
    ]
}

type ImageBits = ([u64; 3], u32, Option<u64>, Option<u64>);
type ProbeBits = ([u64; 3], u32, Option<u64>, Option<u64>);

fn image_bits(job: &PathLuminanceJob<'_, '_>) -> Vec<Option<ImageBits>> {
    job.estimates()
        .map(|value| {
            value.map(|v| {
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

#[test]
fn g4_g5_tabulated_multiple_scattering_images_and_probes_resume_with_identical_bits() {
    with_cx(|cx, gate| {
        let scene = room();
        let sources = room_sources();
        let density = [0.2, 1.0, 0.6, 0.4, 1.4, 0.3, 0.7, 1.1];
        let medium = room_haze(cx, &density, linear_phase(0.75));
        let policy = Settings {
            width: 4,
            height: 3,
            ..iid_settings(4, 32)
        };
        let new_image = |depth| {
            PathLuminanceJob::new(
                &scene,
                &sources,
                PathMedium::Scattering(&medium),
                TransportOptions::default(),
                Settings {
                    max_depth: depth,
                    ..policy
                },
                cx,
            )
            .unwrap()
        };
        let probes = [
            IlluminanceProbe::try_new(Point3::new(0.0, 0.0, 0.1), Vec3::new(0.0, 0.0, 1.0))
                .unwrap(),
            IlluminanceProbe::try_new(Point3::new(0.5, 0.2, 0.1), Vec3::new(0.0, 0.0, 1.0))
                .unwrap(),
        ];
        let new_probes = || {
            PathIlluminanceJob::new(
                &scene,
                &sources,
                PathMedium::Scattering(&medium),
                &probes,
                TransportOptions::default(),
                IlluminanceSampling {
                    samples_per_probe: 32,
                    max_depth: 4,
                    sampler: Sampler::Iid,
                    seed: 173,
                },
                cx,
            )
            .unwrap()
        };
        let mut first_order = new_image(1);
        first_order.advance(cx, usize::MAX).unwrap();
        let mut reference = new_image(4);
        reference.advance(cx, usize::MAX).unwrap();
        let sum = |job: &PathLuminanceJob<'_, '_>| {
            job.estimates()
                .flatten()
                .map(|v| v.luminance_cd_m2)
                .sum::<f64>()
        };
        assert!(
            sum(&reference) > 1.05 * sum(&first_order),
            "higher transport orders must contribute"
        );
        let mut probe_reference = new_probes();
        probe_reference.advance(cx, usize::MAX).unwrap();
        assert!(
            probe_reference
                .estimates()
                .flatten()
                .all(|v| v.direct_lux > 0.0 && v.continuous_lux > 0.0)
        );

        let mut image = new_image(4);
        let mut probe_job = new_probes();
        let execution =
            RenderExecutionConfig::try_new(2, 2, 2, 64 << 20, RunId(0x70686131)).unwrap();
        image.advance(cx, 7).unwrap();
        image.advance_with_execution(cx, 13, &execution).unwrap();
        probe_job.advance_with_execution(cx, 5, &execution).unwrap();
        let image_prefix = image_bits(&image);
        let probe_prefix = probe_bits(&probe_job);
        gate.request();
        assert!(matches!(
            image.advance_with_execution(cx, 19, &execution),
            Err(PhotometryError::Cancelled)
        ));
        assert!(matches!(
            probe_job.advance(cx, 17),
            Err(PhotometryError::Cancelled)
        ));
        assert_eq!(image.completed_path_samples(), 20);
        assert_eq!(probe_job.completed_path_samples(), 5);
        assert_eq!(image_bits(&image), image_prefix);
        assert_eq!(probe_bits(&probe_job), probe_prefix);
        with_cx(|next, _| {
            let execution =
                RenderExecutionConfig::try_new(3, 1, 4, 64 << 20, RunId(0x70686132)).unwrap();
            image.advance_with_execution(next, 37, &execution).unwrap();
            image.advance(next, usize::MAX).unwrap();
            probe_job.advance(next, 7).unwrap();
            probe_job
                .advance_with_execution(next, usize::MAX, &execution)
                .unwrap();
        });
        assert!(image.is_complete() && probe_job.is_complete());
        assert_eq!(image_bits(&image), image_bits(&reference));
        assert_eq!(probe_bits(&probe_job), probe_bits(&probe_reference));
        assert!(
            image
                .estimates()
                .flatten()
                .all(|v| v.total_error_bound_cd_m2.is_none())
        );
        assert!(
            probe_job
                .estimates()
                .flatten()
                .all(|v| v.total_error_bound_lux.is_none())
        );
    });
}

fn hg_table(count: usize, g: f64) -> PhaseFunction {
    PhaseFunction::tabulated(
        &(0..count)
            .map(|index| {
                let mu = 2.0 * index as f64 / (count - 1) as f64 - 1.0;
                (mu, (1.0 - g * g) / (1.0 + g * g - 2.0 * g * mu).powf(1.5))
            })
            .collect::<Vec<_>>(),
    )
    .unwrap()
}

#[test]
#[ignore = "host-specific analytic and tabulated HG multiple-scattering throughput qualification"]
fn photometric_tabulated_phase_throughput() {
    with_cx(|cx, _| {
        let scene = room();
        let sources = room_sources();
        let density = [0.2, 1.0, 0.6, 0.4, 1.4, 0.3, 0.7, 1.1];
        let execution =
            RenderExecutionConfig::try_new(4, 4, 4, 64 << 20, RunId(0x70686133)).unwrap();
        for (label, phase) in [
            (
                "analytic_hg",
                PhaseFunction::henyey_greenstein(0.65).unwrap(),
            ),
            ("table_257", hg_table(257, 0.65)),
            ("table_1025", hg_table(1025, 0.65)),
        ] {
            let medium = room_haze(cx, &density, phase);
            let mut times = [0.0; 3];
            let mut replay = None;
            let mut mean = 0.0;
            for (repeat, seconds) in times.iter_mut().enumerate() {
                let mut job = PathLuminanceJob::new(
                    &scene,
                    &sources,
                    PathMedium::Scattering(&medium),
                    TransportOptions::default(),
                    Settings {
                        width: 16,
                        height: 16,
                        ..iid_settings(4, 64)
                    },
                    cx,
                )
                .unwrap();
                let start = std::time::Instant::now();
                job.advance_with_execution(cx, usize::MAX, &execution)
                    .unwrap();
                *seconds = start.elapsed().as_secs_f64();
                assert!(job.is_complete());
                assert_eq!(job.completed_path_samples(), 16_384);
                let bits = image_bits(&job);
                if let Some(expected) = &replay {
                    assert_eq!(&bits, expected, "{label}, repeat {repeat}");
                } else {
                    mean = job
                        .estimates()
                        .flatten()
                        .map(|v| v.luminance_cd_m2)
                        .sum::<f64>()
                        / 256.0;
                    assert!(mean > 0.0);
                    replay = Some(bits);
                }
            }
            times.sort_by(f64::total_cmp);
            eprintln!(
                "photometric_tabulated_phase_throughput phase={label} workers=4 paths=16384 repeats=3 seconds_median={:.6} paths_per_second={:.1} mean_cd_m2={mean:.9} seconds_sorted={times:?}",
                times[1],
                16_384.0 / times[1],
            );
        }
    });
}
