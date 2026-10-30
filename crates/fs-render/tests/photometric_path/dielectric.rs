//! G1/G3/G4/G5: calibrated finite sources use the native dielectric transport.
//! Analytic expectations describe the admitted finite-depth path families,
//! not a bound on omitted paths or on the declared optical input data.

use super::*;
use fs_exec::RunId;
use fs_render::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricSurface, GlassProvenance,
};
use fs_render::photometry::area::{EmissionSidedness, PhotometricAreaEmitter, PhotometricEmitters};
use fs_render::photometry::medium::HomogeneousMedium;
use fs_render::photometry::path::{IlluminanceSampling, PathIlluminanceJob};
use fs_render::photometry::source::{LUMINOUS_EFFICACY_LM_W, photopic_weight_nm};
use fs_render::photometry::{IlluminanceProbe, PhotometricScene};
use fs_render::spectral::lift_rgb;
use fs_render::tracer::RenderExecutionConfig;

const THICKNESS_M: f64 = 0.4;
const SOURCE_HALF_M: f64 = 2.0;
const LUMINANCE: f64 = 100.0;

fn glass(index: f64, extinction_per_m: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(index).unwrap(),
        BeerLambertAbsorption::try_constant(extinction_per_m).unwrap(),
        GlassProvenance::Custom,
    )
}

fn slab_scene(glass: DielectricGlass, surface: DielectricSurface, source_z: f64) -> Scene {
    let h = 4.0;
    let mesh = TriMesh::new(
        vec![
            [-h, -h, -THICKNESS_M],
            [h, -h, -THICKNESS_M],
            [h, h, -THICKNESS_M],
            [-h, h, -THICKNESS_M],
            [-h, -h, 0.0],
            [h, -h, 0.0],
            [h, h, 0.0],
            [-h, h, 0.0],
        ],
        vec![
            [0, 2, 1],
            [0, 3, 2],
            [4, 5, 6],
            [4, 6, 7],
            [0, 1, 5],
            [0, 5, 4],
            [1, 2, 6],
            [1, 6, 5],
            [2, 3, 7],
            [2, 7, 6],
            [3, 0, 4],
            [3, 4, 7],
        ],
    );
    let mut source = floor(source_z, SOURCE_HALF_M, LiftedSpectrum { c: [0.0; 3] });
    if source_z > 1.0 {
        let Shape::Mesh(mesh) = &mut source.shape else {
            unreachable!()
        };
        // A source above the camera emits downward toward the front reflection.
        *mesh = TriMesh::new(mesh.vertices.clone(), vec![[0, 2, 1], [0, 3, 2]]);
    }
    let mut scene = empty_scene([0.0, 0.0, 1.0], [0.0, 0.0, -1.0]);
    scene.primitives.push(Primitive {
        shape: Shape::Mesh(mesh),
        material: Material::Dielectric { glass, surface },
        emission: None,
    });
    scene.primitives.push(source);
    scene
}

fn emitter(spectrum: SourceSpectrum) -> PhotometricAreaEmitter {
    PhotometricAreaEmitter::try_new(
        1,
        LUMINANCE * PI * 4.0 * SOURCE_HALF_M * SOURCE_HALF_M,
        1.0,
        spectrum,
        EmissionSidedness::Front,
    )
    .unwrap()
}

fn pixel(
    scene: &Scene,
    sources: &[PhotometricAreaEmitter],
    policy: Settings,
    cx: &Cx<'_>,
) -> PathLuminanceEstimate {
    let mut job = PathLuminanceJob::new_with_emitters(
        scene,
        PhotometricEmitters {
            points: &[],
            surfaces: sources,
        },
        PathMedium::Vacuum,
        TransportOptions::default(),
        policy,
        cx,
    )
    .unwrap();
    job.advance(cx, usize::MAX).unwrap();
    let value = job.pixel_estimate(0).unwrap();
    assert_eq!(value.samples, policy.spp);
    assert!(value.total_error_bound_cd_m2.is_none());
    value
}

#[test]
fn g1_calibrated_fresnel_reflection_and_transmission_retain_radiance_eta_factors() {
    with_cx(|cx, _| {
        let index: f64 = 1.5;
        let reflected = ((index - 1.0) / (index + 1.0)).powi(2);
        let sources = [emitter(SourceSpectrum::equal_energy())];
        for (source_z, extinction, depth, fraction, label) in [
            (2.0, 0.0, 2, reflected, "front Fresnel reflection"),
            (
                -0.2,
                0.0,
                2,
                (1.0 - reflected) / (index * index),
                "one interface retains its radiance eta factor",
            ),
            (
                -1.0,
                0.7,
                3,
                (1.0 - reflected).powi(2) * (-0.7 * THICKNESS_M).exp(),
                "two interfaces cancel eta factors and retain internal Beer loss",
            ),
        ] {
            let scene = slab_scene(
                glass(index, extinction),
                DielectricSurface::SMOOTH,
                source_z,
            );
            // The native first two smooth-interface splits enumerate reflection
            // and transmission. Depth 3 stops before an internal round trip can
            // return to the source, so the analytic term contains exactly T².
            let value = pixel(&scene, &sources, settings(depth, 8192), cx);
            close(value.luminance_cd_m2, LUMINANCE * fraction, 0.002, label);
        }
    });
}

#[test]
fn g1_calibrated_cauchy_and_colored_beer_slab_match_spectral_quadrature() {
    with_cx(|cx, _| {
        let reference_rgb = [0.98, 0.72, 0.22];
        let reference_distance_m = 0.2;
        let tint = lift_rgb(reference_rgb);
        let spectrum = SourceSpectrum::try_from_bins(std::array::from_fn(|bin| {
            1.0 + 3.0 * bin as f64 / SPECTRAL_BINS as f64
        }))
        .unwrap();
        let optic = DielectricGlass::new(
            CauchyIor::try_new(1.38, 0.018, 0.0).unwrap(),
            BeerLambertAbsorption::try_from_rgb_transmittance(reference_rgb, reference_distance_m)
                .unwrap(),
            GlassProvenance::Custom,
        );
        let mut scene = slab_scene(optic, DielectricSurface::SMOOTH, -1.0);
        // Independently evaluate Fresnel and Beer from the declared spectrum
        // coefficients. The RGB lift defines the input model, not a measurement
        // oracle; no production glass or transport evaluator supplies this sum.
        let sources = [emitter(spectrum.clone())];
        for incident_cosine in [1.0_f64, 0.8] {
            let sine = (1.0 - incident_cosine * incident_cosine).sqrt();
            scene.camera.forward = Vec3::new(sine, 0.0, -incident_cosine);
            let expected_fraction: f64 = (0..SPECTRAL_BINS)
                .map(|bin| {
                    let lambda = wavelength_nm(bin);
                    let um = lambda * 1.0e-3;
                    let index = 1.38 + 0.018 / (um * um);
                    let internal_cosine = (1.0 - (sine / index).powi(2)).sqrt();
                    let rs = (incident_cosine - index * internal_cosine)
                        / (incident_cosine + index * internal_cosine);
                    let rp = (index * incident_cosine - internal_cosine)
                        / (index * incident_cosine + internal_cosine);
                    let reflection = 0.5 * (rs * rs + rp * rp);
                    let x = (lambda - 380.0) / 400.0;
                    let polynomial = (tint.c[0] * x + tint.c[1]) * x + tint.c[2];
                    let reference_transmission =
                        0.5 + 0.5 * polynomial / (1.0 + polynomial * polynomial).sqrt();
                    spectrum.bins_w_per_lumen_nm()[bin]
                        * LUMINOUS_EFFICACY_LM_W
                        * photopic_weight_nm(bin)
                        * (1.0 - reflection).powi(2)
                        * reference_transmission
                            .powf(THICKNESS_M / internal_cosine / reference_distance_m)
                })
                .sum();
            let value = pixel(&scene, &sources, settings(3, 8192), cx);
            close(
                value.luminance_cd_m2,
                LUMINANCE * expected_fraction,
                0.004,
                "full-spectrum Fresnel and dispersive internal Beer distance",
            );
        }
    });
}

#[test]
fn g1_calibrated_glass_displaces_ambient_absorption_without_double_attenuation() {
    with_cx(|cx, _| {
        let density = [1.0];
        let ambient = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [-4.0, -4.0, -2.0], [8.0, 8.0, 4.0]),
            MassAbsorptionSpectrum::constant(0.6).unwrap(),
            None,
            medium_budget(),
        )
        .unwrap();
        assert_eq!(ambient.bounds_m(), ([-4.0, -4.0, -2.0], [4.0, 4.0, 2.0]));
        let sources = [emitter(SourceSpectrum::equal_energy())];
        for (source_z, ambient_length, glass_length, depth) in [
            (-0.2, 1.0, 0.2, 2),
            (-1.0, 2.0 - THICKNESS_M, THICKNESS_M, 3),
        ] {
            let scene = slab_scene(glass(1.0, 1.3), DielectricSurface::SMOOTH, source_z);
            let mut job = PathLuminanceJob::new_with_emitters(
                &scene,
                PhotometricEmitters {
                    points: &[],
                    surfaces: &sources,
                },
                PathMedium::Absorbing(&ambient),
                TransportOptions::default(),
                settings(depth, 8192),
                cx,
            )
            .unwrap();
            job.advance(cx, usize::MAX).unwrap();
            let optical_depth: f64 = 0.6 * ambient_length + 1.3 * glass_length;
            close(
                job.pixel_estimate(0).unwrap().luminance_cd_m2,
                LUMINANCE * (-optical_depth).exp(),
                0.002,
                "only the active material attenuates each physical segment",
            );
        }
    });
}

fn slab_probe_reference(index: f64, extinction: f64, panels: u32) -> f64 {
    // Integrate incident directions at an ambient probe 1 m above the slab.
    // Source is 1 m below its top. Snell changes the internal footprint while
    // the two radiance eta factors cancel. For the square source, azimuth
    // support has a closed form; only a one-dimensional mu integral remains.
    let air_distance = 2.0 - THICKNESS_M;
    let mut sum = 0.0;
    for panel in 0..panels {
        let mu = (f64::from(panel) + 0.5) / f64::from(panels);
        let sin_incident = (1.0 - mu * mu).sqrt();
        let sin_internal = sin_incident / index;
        let mu_internal = (1.0 - sin_internal * sin_internal).sqrt();
        let radius = air_distance * sin_incident / mu + THICKNESS_M * sin_internal / mu_internal;
        let azimuth = if radius <= SOURCE_HALF_M {
            2.0 * PI
        } else if radius < SOURCE_HALF_M * 2.0_f64.sqrt() {
            8.0 * (SOURCE_HALF_M / radius).asin() - 2.0 * PI
        } else {
            0.0
        };
        let rs = (mu - index * mu_internal) / (mu + index * mu_internal);
        let rp = (index * mu - mu_internal) / (index * mu + mu_internal);
        let transmission = 1.0 - 0.5 * (rs * rs + rp * rp);
        sum += mu
            * azimuth
            * transmission
            * transmission
            * (-extinction * THICKNESS_M / mu_internal).exp();
    }
    LUMINANCE * sum / f64::from(panels)
}

#[test]
fn g1_calibrated_glass_probe_matches_refracted_square_source_integral() {
    with_cx(|cx, _| {
        let scene = slab_scene(glass(1.5, 0.7), DielectricSurface::SMOOTH, -1.0);
        let sources = [emitter(SourceSpectrum::equal_energy())];
        let probes =
            [
                IlluminanceProbe::try_new(Point3::new(0.0, 0.0, 1.0), Vec3::new(0.0, 0.0, -1.0))
                    .unwrap(),
            ];
        let mut job = PathIlluminanceJob::new_with_emitters(
            &scene,
            PhotometricEmitters {
                points: &[],
                surfaces: &sources,
            },
            PathMedium::Vacuum,
            &probes,
            TransportOptions::default(),
            IlluminanceSampling {
                samples_per_probe: 16_384,
                max_depth: 3,
                sampler: Sampler::OwenSobol,
                seed: 173,
            },
            cx,
        )
        .unwrap();
        job.advance(cx, usize::MAX).unwrap();
        let estimate = job.probe_estimate(0).unwrap();
        let reference = slab_probe_reference(1.5, 0.7, 65_536);
        close(
            slab_probe_reference(1.5, 0.7, 32_768),
            reference,
            2.0e-5,
            "independent angular reference refinement",
        );
        close(
            estimate.illuminance_lux,
            reference,
            0.02,
            "incident lux through a refracting and absorbing slab",
        );
        assert_eq!(estimate.direct_lux, 0.0);
        assert!(estimate.estimated_standard_error_lux.is_none());
        assert!(estimate.total_error_bound_lux.is_none());
    });
}

#[test]
fn g3_calibrated_rough_glass_nee_mis_and_bsdf_agree_and_absorb() {
    with_cx(|cx, _| {
        let surface = DielectricSurface::try_rough(0.25).unwrap();
        let scene = slab_scene(glass(1.5, 0.0), surface, -0.2);
        let sources = [emitter(SourceSpectrum::equal_energy())];
        let estimates = [
            DirectStrategy::BsdfOnly,
            DirectStrategy::NeeOnly,
            DirectStrategy::Mis,
            DirectStrategy::PowerMis,
        ]
        .map(|strategy| {
            pixel(
                &scene,
                &sources,
                Settings {
                    strategy,
                    sampler: Sampler::Iid,
                    ..settings(6, 8192)
                },
                cx,
            )
        });
        let reference = estimates[0];
        assert!(reference.luminance_cd_m2 > 0.0);
        for estimate in &estimates[1..] {
            let allowance = 6.0
                * (reference.estimated_standard_error_cd_m2.unwrap().powi(2)
                    + estimate.estimated_standard_error_cd_m2.unwrap().powi(2))
                .sqrt()
                + 0.005 * reference.luminance_cd_m2;
            // A deterministic qualification band based on estimated IID noise;
            // it is not advertised as a guaranteed confidence interval.
            assert!(
                (estimate.luminance_cd_m2 - reference.luminance_cd_m2).abs() <= allowance,
                "rough-glass strategies disagree: reference={reference:?}, candidate={estimate:?}, allowance={allowance}"
            );
        }
        let absorbed = pixel(
            &slab_scene(glass(1.5, 4.0), surface, -0.2),
            &sources,
            Settings {
                sampler: Sampler::Iid,
                ..settings(6, 8192)
            },
            cx,
        );
        assert!(absorbed.luminance_cd_m2 > 0.0);
        assert!(absorbed.luminance_cd_m2 < 0.7 * estimates[2].luminance_cd_m2);
    });
}

fn image_bits(job: &PathLuminanceJob<'_, '_>) -> Vec<Option<([u64; 3], u32, Option<u64>)>> {
    job.estimates()
        .map(|value| {
            value.map(|estimate| {
                assert!(estimate.total_error_bound_cd_m2.is_none());
                (
                    estimate.xyz.map(f64::to_bits),
                    estimate.samples,
                    estimate.estimated_standard_error_cd_m2.map(f64::to_bits),
                )
            })
        })
        .collect()
}

fn probe_bits(job: &PathIlluminanceJob<'_, '_>) -> Vec<Option<([u64; 3], u32, Option<u64>)>> {
    job.estimates()
        .map(|value| {
            value.map(|estimate| {
                assert!(estimate.total_error_bound_lux.is_none());
                (
                    [
                        estimate.direct_lux,
                        estimate.continuous_lux,
                        estimate.illuminance_lux,
                    ]
                    .map(f64::to_bits),
                    estimate.samples,
                    estimate.estimated_standard_error_lux.map(f64::to_bits),
                )
            })
        })
        .collect()
}

#[test]
fn g4_g5_calibrated_glass_parallel_cancellation_and_mixed_resume_preserve_bits() {
    with_cx(|cx, gate| {
        let mut scene = slab_scene(glass(1.5, 0.7), DielectricSurface::SMOOTH, -1.0);
        scene.camera.half_tan = 0.25;
        let sources = [emitter(SourceSpectrum::equal_energy())];
        let emitters = PhotometricEmitters {
            points: &[],
            surfaces: &sources,
        };
        let probes = [-0.5, 0.0, 0.5].map(|x| {
            IlluminanceProbe::try_new(Point3::new(x, 0.0, 1.0), Vec3::new(0.0, 0.0, -1.0)).unwrap()
        });
        let new_image = || {
            PathLuminanceJob::new_with_emitters(
                &scene,
                emitters,
                PathMedium::Vacuum,
                TransportOptions::default(),
                Settings {
                    width: 3,
                    height: 2,
                    sampler: Sampler::Iid,
                    ..settings(5, 16)
                },
                cx,
            )
            .unwrap()
        };
        let new_probes = || {
            PathIlluminanceJob::new_with_emitters(
                &scene,
                emitters,
                PathMedium::Vacuum,
                &probes,
                TransportOptions::default(),
                IlluminanceSampling {
                    samples_per_probe: 16,
                    max_depth: 5,
                    sampler: Sampler::Iid,
                    seed: 173,
                },
                cx,
            )
            .unwrap()
        };
        let mut image_reference = new_image();
        let mut probe_reference = new_probes();
        image_reference.advance(cx, usize::MAX).unwrap();
        probe_reference.advance(cx, usize::MAX).unwrap();
        let mut image = new_image();
        let mut probes = new_probes();
        let execution =
            RenderExecutionConfig::try_new(2, 2, 2, 64 << 20, RunId(0x676c6173)).unwrap();
        image.advance(cx, 3).unwrap();
        probes.advance_with_execution(cx, 7, &execution).unwrap();
        let image_prefix = image_bits(&image);
        let probe_prefix = probe_bits(&probes);
        gate.request();
        assert!(matches!(
            image.advance_with_execution(cx, 16, &execution),
            Err(PhotometryError::Cancelled)
        ));
        assert!(matches!(
            probes.advance_with_execution(cx, 16, &execution),
            Err(PhotometryError::Cancelled)
        ));
        assert_eq!(image_bits(&image), image_prefix);
        assert_eq!(probe_bits(&probes), probe_prefix);
        with_cx(|resumed, _| {
            let execution =
                RenderExecutionConfig::try_new(3, 1, 4, 64 << 20, RunId(0x676c6174)).unwrap();
            image
                .advance_with_execution(resumed, 23, &execution)
                .unwrap();
            image.advance(resumed, usize::MAX).unwrap();
            probes.advance(resumed, 5).unwrap();
            probes
                .advance_with_execution(resumed, usize::MAX, &execution)
                .unwrap();
        });
        assert_eq!(image_bits(&image), image_bits(&image_reference));
        assert_eq!(probe_bits(&probes), probe_bits(&probe_reference));
    });
}

#[test]
fn g0_direct_adapter_and_unsupported_point_glass_connections_still_refuse() {
    with_cx(|cx, _| {
        let scene = slab_scene(glass(1.5, 0.0), DielectricSurface::SMOOTH, -1.0);
        let surfaces = [emitter(SourceSpectrum::equal_energy())];
        let points = [isotropic(
            [0.0, 0.0, 2.0],
            10.0,
            SourceSpectrum::equal_energy(),
        )];
        assert!(matches!(
            PhotometricScene::try_new(
                &scene.primitives,
                &[],
                HomogeneousMedium::vacuum(),
                TransportOptions::default(),
                cx
            ),
            Err(PhotometryError::UnsupportedPrimitive { primitive: 0, .. })
        ));
        assert!(matches!(
            PathLuminanceJob::new_with_emitters(
                &scene,
                PhotometricEmitters {
                    points: &points,
                    surfaces: &surfaces
                },
                PathMedium::Vacuum,
                TransportOptions::default(),
                settings(3, 8),
                cx,
            ),
            Err(PhotometryError::UnsupportedPrimitive { primitive: 0, .. })
        ));
        let invalid_emitter = [PhotometricAreaEmitter::try_new(
            0,
            100.0,
            1.0,
            SourceSpectrum::equal_energy(),
            EmissionSidedness::Both,
        )
        .unwrap()];
        assert!(matches!(
            PathLuminanceJob::new_with_emitters(
                &scene,
                PhotometricEmitters {
                    points: &[],
                    surfaces: &invalid_emitter
                },
                PathMedium::Vacuum,
                TransportOptions::default(),
                settings(3, 8),
                cx,
            ),
            Err(PhotometryError::UnsupportedPrimitive { primitive: 0, .. })
        ));
    });
}

#[test]
#[ignore = "run explicitly to qualify calibrated native glass CPU throughput"]
fn photometric_glass_parallel_throughput() {
    use std::time::Instant;

    with_cx(|cx, _| {
        let sources = [emitter(SourceSpectrum::equal_energy())];
        for (label, surface) in [
            ("smooth", DielectricSurface::SMOOTH),
            ("rough", DielectricSurface::try_rough(0.2).unwrap()),
        ] {
            let mut scene = slab_scene(glass(1.5, 0.7), surface, -1.0);
            scene.camera.half_tan = 0.4;
            let new_job = || {
                PathLuminanceJob::new_with_emitters(
                    &scene,
                    PhotometricEmitters {
                        points: &[],
                        surfaces: &sources,
                    },
                    PathMedium::Vacuum,
                    TransportOptions::default(),
                    Settings {
                        width: 32,
                        height: 32,
                        sampler: Sampler::Iid,
                        ..settings(5, 128)
                    },
                    cx,
                )
                .unwrap()
            };
            let mut serial = new_job();
            let start = Instant::now();
            serial.advance(cx, usize::MAX).unwrap();
            let seconds = start.elapsed().as_secs_f64();
            let paths = serial.completed_path_samples();
            assert_eq!(paths, 131_072);
            assert!(
                serial
                    .estimates()
                    .flatten()
                    .all(|pixel| pixel.luminance_cd_m2 > 0.0)
            );
            let expected = image_bits(&serial);
            eprintln!(
                "photometric_glass_throughput surface={label} execution=serial paths={paths} seconds={seconds:.6} paths_per_second={:.1}",
                paths as f64 / seconds,
            );
            for workers in [1, 4, 8] {
                let mut job = new_job();
                let execution = RenderExecutionConfig::try_new(
                    8,
                    8,
                    workers,
                    64 << 20,
                    RunId(0x676c_6200 + workers as u64),
                )
                .unwrap();
                let start = Instant::now();
                job.advance_with_execution(cx, usize::MAX, &execution)
                    .unwrap();
                let seconds = start.elapsed().as_secs_f64();
                assert!(job.is_complete());
                assert_eq!(
                    image_bits(&job),
                    expected,
                    "surface={label}, workers={workers}"
                );
                eprintln!(
                    "photometric_glass_throughput surface={label} workers={workers} paths={paths} seconds={seconds:.6} paths_per_second={:.1}",
                    paths as f64 / seconds,
                );
            }
        }
    });
}
