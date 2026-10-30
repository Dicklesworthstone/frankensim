//! G1/G3/G4/G5: sampled point connections retain calibrated spectral transport.

use super::*;
use fs_exec::RunId;
use fs_render::photometry::IlluminanceProbe;
use fs_render::photometry::path::{
    IlluminanceSampling, MAX_POINT_CONNECTIONS_PER_VERTEX, PathIlluminanceJob, PointFixtureSampling,
};
use fs_render::tracer::RenderExecutionConfig;

fn importance(connections_per_vertex: u32) -> PointFixtureSampling {
    PointFixtureSampling::Importance {
        connections_per_vertex,
    }
}

fn image_job<'scene, 'field>(
    scene: &'scene Scene,
    sources: &'scene [PhotometricSource],
    medium: PathMedium<'scene, 'field>,
    settings: Settings,
    policy: PointFixtureSampling,
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
    .with_point_sampling(policy)
    .unwrap()
}

fn iid_settings(depth: u32, spp: u32) -> Settings {
    Settings {
        sampler: Sampler::Iid,
        ..settings(depth, spp)
    }
}

fn probe_at(x: f64, y: f64) -> IlluminanceProbe {
    IlluminanceProbe::try_new(Point3::new(x, y, 0.0), Vec3::new(0.0, 0.0, 1.0)).unwrap()
}

fn probe_job<'scene, 'field>(
    scene: &'scene Scene,
    sources: &'scene [PhotometricSource],
    medium: PathMedium<'scene, 'field>,
    probes: &'scene [IlluminanceProbe],
    settings: Settings,
    policy: PointFixtureSampling,
    cx: &Cx<'_>,
) -> PathIlluminanceJob<'scene, 'field> {
    PathIlluminanceJob::new(
        scene,
        sources,
        medium,
        probes,
        TransportOptions::default(),
        IlluminanceSampling {
            samples_per_probe: settings.spp,
            max_depth: settings.max_depth,
            sampler: settings.sampler,
            seed: settings.seed,
        },
        cx,
    )
    .unwrap()
    .with_point_sampling(policy)
    .unwrap()
}

fn assert_statistical_match(actual: f64, expected: f64, se: Option<f64>, reference_error: f64) {
    let se = se.expect("IID paths must report estimated sampling standard error");
    assert!(se.is_finite() && se >= 0.0);
    // A deterministic regression threshold using an estimated SE, not a
    // confidence certificate. The separate allowance belongs to the oracle.
    let tolerance = 6.0 * se + reference_error + 2.0e-10 * expected.abs().max(1.0);
    assert!(
        actual.is_finite() && (actual - expected).abs() <= tolerance,
        "observed={actual:.12e} expected={expected:.12e} estimated_se={se:.6e} tolerance={tolerance:.6e}"
    );
}

struct FixtureReference {
    position: [f64; 3],
    intensity_to_origin_cd: f64,
    spectral_shape: [f64; SPECTRAL_BINS],
}

fn colored_rig() -> (Vec<PhotometricSource>, Vec<FixtureReference>) {
    let positions = [
        [0.0, 0.0, 3.0],
        [-2.0, 1.0, 4.0],
        [1.0, -1.0, 3.0],
        [-1.0, 0.0, 3.0],
        [0.5, 0.0, 3.5],
    ];
    let mut sources = Vec::new();
    let mut references = Vec::new();
    for (index, position) in positions.into_iter().enumerate() {
        let spectral_shape = std::array::from_fn(|bin| {
            let x = (wavelength_nm(bin) - 380.0) / 400.0;
            match index {
                1 => 0.1 + 3.0 * x * x,
                2 => 0.2 + 2.0 * (1.0 - x).powi(2),
                3 => 0.1 + 4.0 * x * (1.0 - x),
                _ => 1.0,
            }
        });
        let r2: f64 = position.iter().map(|value| value * value).sum();
        let cosine = position[2] / r2.sqrt();
        let (rotation, flux, dimmer, distribution, intensity) = match index {
            0 => (
                [0.0, 0.0, 0.0, 1.0],
                4.0 * PI * 90.0,
                1.0,
                AngularDistribution::isotropic(),
                90.0,
            ),
            1 => (
                [0.0, 0.0, 0.0, 1.0],
                4.0 * PI * 25.0,
                0.6,
                AngularDistribution::isotropic(),
                15.0,
            ),
            2 => (
                [1.0, 0.0, 0.0, 0.0],
                200.0,
                0.7,
                AngularDistribution::cosine_power(2.0, FRAC_PI_2).unwrap(),
                140.0 * 3.0 * cosine * cosine / (2.0 * PI),
            ),
            3 => (
                [1.0, 0.0, 0.0, 0.0],
                300.0,
                0.2,
                AngularDistribution::uniform_cone(0.7).unwrap(),
                60.0 / (2.0 * PI * (1.0 - 0.7_f64.cos())),
            ),
            // Large total flux directed away from the measured surface must
            // neither illuminate it nor consume the proposal's useful weight.
            _ => (
                [0.0, 0.0, 0.0, 1.0],
                100_000.0,
                1.0,
                AngularDistribution::uniform_cone(0.2).unwrap(),
                0.0,
            ),
        };
        sources.push(
            PhotometricSource::try_new(
                RigidTransform::try_new(rotation, position).unwrap(),
                flux,
                dimmer,
                distribution,
                SourceSpectrum::try_from_bins(spectral_shape).unwrap(),
            )
            .unwrap(),
        );
        references.push(FixtureReference {
            position,
            intensity_to_origin_cd: intensity,
            spectral_shape,
        });
    }
    (sources, references)
}

fn spectral_lambertian_reference(
    fixtures: &[FixtureReference],
    blocked: Option<usize>,
    absorption: bool,
) -> f64 {
    // This uses the declared observer, but independently normalizes the raw
    // source bins and expands the sigmoid material rather than calling either
    // the source-intensity evaluator or a production BRDF/transport helper.
    fixtures
        .iter()
        .enumerate()
        .filter(|(index, _)| Some(*index) != blocked)
        .map(|(_, fixture)| {
            let r2: f64 = fixture.position.iter().map(|v| v * v).sum();
            let distance = r2.sqrt();
            let incident_lux =
                fixture.intensity_to_origin_cd * fixture.position[2] / (r2 * distance);
            let mut normalization = 0.0;
            let mut reflected = 0.0;
            for (bin, &shape) in fixture.spectral_shape.iter().enumerate() {
                let lambda = wavelength_nm(bin);
                let x = (lambda - 380.0) / 400.0;
                let polynomial = 2.0 * x - 0.5;
                let rho = 0.5 + 0.5 * polynomial / (1.0 + polynomial * polynomial).sqrt();
                let transmission = if absorption {
                    (-(0.05 + 0.15 * x) * (distance + 2.0)).exp()
                } else {
                    1.0
                };
                let photopic = shape * super::spectral_reference::cie_reference_xyz(lambda)[1];
                normalization += photopic;
                reflected += photopic * rho * transmission;
            }
            incident_lux * reflected / normalization / PI
        })
        .sum()
}

fn colored_receiver() -> Scene {
    let mut scene = receiver_scene();
    scene.primitives[0].material = Material::Lambertian {
        reflectance: LiftedSpectrum {
            c: [0.0, 2.0, -0.5],
        },
    };
    scene
}

fn rectangle(x0: f64, x1: f64, y0: f64, y1: f64, z: f64) -> Primitive {
    Primitive {
        shape: Shape::Mesh(TriMesh::new(
            vec![[x0, y0, z], [x1, y0, z], [x1, y1, z], [x0, y1, z]],
            vec![[0, 1, 2], [0, 2, 3]],
        )),
        material: Material::Lambertian {
            reflectance: LiftedSpectrum { c: [0.0; 3] },
        },
        emission: None,
    }
}

#[test]
fn g0_point_sampling_policy_is_bounded_and_cannot_change_accepted_paths() {
    with_cx(|cx, _| {
        let scene = receiver_scene();
        let sources = [isotropic(
            [0.0, 0.0, 3.0],
            90.0,
            SourceSpectrum::equal_energy(),
        )];
        let settings = iid_settings(1, 32);
        let probes = [probe_at(0.0, 0.0)];
        let mut original = PathLuminanceJob::new(
            &scene,
            &sources,
            PathMedium::Vacuum,
            TransportOptions::default(),
            settings,
            cx,
        )
        .unwrap();
        assert_eq!(original.point_sampling(), PointFixtureSampling::All);
        let mut explicit = image_job(
            &scene,
            &sources,
            PathMedium::Vacuum,
            settings,
            PointFixtureSampling::All,
            cx,
        );
        original.advance(cx, usize::MAX).unwrap();
        explicit.advance(cx, usize::MAX).unwrap();
        assert_eq!(
            estimate_bits(original.pixel_estimate(0).unwrap()),
            estimate_bits(explicit.pixel_estimate(0).unwrap())
        );
        let mut image = image_job(
            &scene,
            &sources,
            PathMedium::Vacuum,
            settings,
            importance(1),
            cx,
        );
        let mut probe = probe_job(
            &scene,
            &sources,
            PathMedium::Vacuum,
            &probes,
            settings,
            importance(1),
            cx,
        );
        for invalid in [0, MAX_POINT_CONNECTIONS_PER_VERTEX + 1] {
            assert!(image.set_point_sampling(importance(invalid)).is_err());
            assert!(probe.set_point_sampling(importance(invalid)).is_err());
            assert_eq!(image.point_sampling(), importance(1));
            assert_eq!(probe.point_sampling(), importance(1));
        }
        image
            .set_point_sampling(importance(MAX_POINT_CONNECTIONS_PER_VERTEX))
            .unwrap();
        probe
            .set_point_sampling(importance(MAX_POINT_CONNECTIONS_PER_VERTEX))
            .unwrap();
        image.advance(cx, 1).unwrap();
        probe.advance(cx, 1).unwrap();
        let accepted = estimate_bits(image.pixel_estimate(0).unwrap());
        let accepted_probe = probe.probe_estimate(0).unwrap();
        assert!(image.set_point_sampling(PointFixtureSampling::All).is_err());
        assert!(probe.set_point_sampling(PointFixtureSampling::All).is_err());
        assert_eq!(estimate_bits(image.pixel_estimate(0).unwrap()), accepted);
        assert_eq!(probe.probe_estimate(0).unwrap(), accepted_probe);
        assert_eq!(image.completed_path_samples(), 1);
        assert_eq!(probe.completed_path_samples(), 1);

        // A rig with no candidate in the receiving hemisphere remains dark.
        let dark = [isotropic(
            [0.0, 0.0, -3.0],
            90.0,
            SourceSpectrum::equal_energy(),
        )];
        let mut empty = image_job(&scene, &[], PathMedium::Vacuum, settings, importance(1), cx);
        let mut backwards = image_job(
            &scene,
            &dark,
            PathMedium::Vacuum,
            settings,
            importance(1),
            cx,
        );
        for job in [&mut empty, &mut backwards] {
            job.advance(cx, usize::MAX).unwrap();
            let estimate = job.pixel_estimate(0).unwrap();
            assert_eq!(estimate.luminance_cd_m2, 0.0);
            assert_eq!(estimate.estimated_standard_error_cd_m2, Some(0.0));
        }
    });
}

#[test]
fn g1_sampled_unequal_spectra_and_beams_match_independent_eighty_bin_luminance() {
    with_cx(|cx, _| {
        let scene = colored_receiver();
        let (sources, references) = colored_rig();
        let expected = spectral_lambertian_reference(&references, None, false);
        for policy in [PointFixtureSampling::All, importance(1), importance(4)] {
            let mut job = image_job(
                &scene,
                &sources,
                PathMedium::Vacuum,
                iid_settings(1, 16_384),
                policy,
                cx,
            );
            job.advance(cx, usize::MAX).unwrap();
            let estimate = job.pixel_estimate(0).unwrap();
            assert_statistical_match(
                estimate.luminance_cd_m2,
                expected,
                estimate.estimated_standard_error_cd_m2,
                0.0,
            );
            assert!(estimate.estimated_standard_error_cd_m2.unwrap() < 0.01 * expected);
            assert!(estimate.total_error_bound_cd_m2.is_none());
        }
    });
}

#[test]
fn g1_sampled_connections_retain_real_occlusion_and_both_spectral_beer_legs() {
    with_cx(|cx, _| {
        let mut scene = colored_receiver();
        // Only fixture 1 crosses (-1, 0.5, 2). The primary camera ray and
        // every other nonzero source connection miss this finite blocker.
        scene
            .primitives
            .push(rectangle(-1.15, -0.85, 0.35, 0.65, 2.0));
        let (sources, references) = colored_rig();
        let density = [1.0];
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [-5.0, -5.0, 0.0], [10.0, 10.0, 5.0]),
            MassAbsorptionSpectrum::table(&[(380.0, 0.05), (780.0, 0.20)]).unwrap(),
            None,
            medium_budget(),
        )
        .unwrap();
        let expected = spectral_lambertian_reference(&references, Some(1), true);
        assert!(expected < spectral_lambertian_reference(&references, None, true));
        for policy in [PointFixtureSampling::All, importance(1), importance(4)] {
            let mut job = image_job(
                &scene,
                &sources,
                PathMedium::Absorbing(&medium),
                iid_settings(1, 16_384),
                policy,
                cx,
            );
            job.advance(cx, usize::MAX).unwrap();
            let estimate = job.pixel_estimate(0).unwrap();
            assert_statistical_match(
                estimate.luminance_cd_m2,
                expected,
                estimate.estimated_standard_error_cd_m2,
                0.0,
            );
        }
    });
}

const HAZE_SOURCES: [([f64; 3], f64); 4] = [
    ([0.0, 0.0, 2.0], 100.0),
    ([0.0, 0.0, 4.0], 500.0),
    ([0.0, 0.0, -1.0], 50.0),
    ([1.0, 0.0, 2.0], 80.0),
];

fn haze_reference(hg: bool, panels: u32) -> f64 {
    let step = 1.0 / f64::from(panels);
    let mut integral = 0.0;
    for panel in 0..panels {
        let z = (f64::from(panel) + 0.5) * step;
        for (position, candela) in HAZE_SOURCES {
            let dz = position[2] - z;
            let r2 = position[0].powi(2) + position[1].powi(2) + dz * dz;
            let r = r2.sqrt();
            let cosine = dz / r;
            let phase = if hg {
                let g = 0.65;
                (1.0 - g * g) / (4.0 * PI * (1.0 + g * g - 2.0 * g * cosine).powf(1.5))
            } else {
                3.0 * (1.0 + cosine * cosine) / (16.0 * PI)
            };
            let source_medium_distance = if dz > 0.0 {
                r * (1.0 - z) / dz
            } else {
                r * z / -dz
            };
            integral +=
                step * candela / r2 * 0.3 * phase * (-0.5 * (z + source_medium_distance)).exp();
        }
    }
    integral
}

#[test]
fn g1_sampled_haze_connections_match_hg_and_rayleigh_single_scatter_integrals() {
    with_cx(|cx, _| {
        let scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        let sources: Vec<_> = HAZE_SOURCES
            .into_iter()
            .map(|(p, i)| isotropic(p, i, SourceSpectrum::equal_energy()))
            .collect();
        let density = [1.0];
        for hg in [true, false] {
            let phase = if hg {
                PhaseFunction::henyey_greenstein(0.65).unwrap()
            } else {
                PhaseFunction::Rayleigh
            };
            let medium = ScatteringGrid::try_new(
                cx,
                AbsorbingGrid::try_new(
                    cx,
                    VolumeGrid::new([1; 3], &density, [-2.0, -2.0, 0.0], [4.0, 4.0, 1.0]),
                    MassAbsorptionSpectrum::constant(0.2).unwrap(),
                    None,
                    medium_budget(),
                )
                .unwrap(),
                MassScatteringSpectrum::constant(0.3).unwrap(),
                phase,
            )
            .unwrap();
            let expected = haze_reference(hg, 4096);
            close(
                haze_reference(hg, 2048),
                expected,
                2.0e-7,
                "independent haze quadrature refinement",
            );
            for policy in [PointFixtureSampling::All, importance(1), importance(4)] {
                let mut job = image_job(
                    &scene,
                    &sources,
                    PathMedium::Scattering(&medium),
                    iid_settings(1, 32_768),
                    policy,
                    cx,
                );
                job.advance(cx, usize::MAX).unwrap();
                let estimate = job.pixel_estimate(0).unwrap();
                assert_statistical_match(
                    estimate.luminance_cd_m2,
                    expected,
                    estimate.estimated_standard_error_cd_m2,
                    2.0e-7 * expected,
                );
            }
        }
    });
}

#[test]
fn g1_probe_direct_lux_stays_exact_while_reflected_connections_are_sampled() {
    with_cx(|cx, _| {
        let mut scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        scene
            .primitives
            .push(floor(2.0, 100.0, LiftedSpectrum { c: [0.0; 3] }));
        let sources: Vec<_> = [1.0, 7.0, 20.0, 72.0]
            .into_iter()
            .map(|intensity| isotropic([0.0, 0.0, 1.0], intensity, SourceSpectrum::equal_energy()))
            .collect();
        let probes = [probe_at(0.0, 0.0)];
        // Independent infinite-ceiling integral; finite square missing tail
        // is below 8*rho*I/(5*100^5) = 8e-9 lx for this rho=0.5, I=100 cd.
        let reflected = 50.0 * (1.0 - 4.0 * PI / (9.0 * 3.0_f64.sqrt()));
        let mut exact_bits = None;
        for policy in [PointFixtureSampling::All, importance(1), importance(4)] {
            let mut job = probe_job(
                &scene,
                &sources,
                PathMedium::Vacuum,
                &probes,
                iid_settings(1, 32_768),
                policy,
                cx,
            );
            job.advance(cx, 1).unwrap();
            let cached = job.probe_estimate(0).unwrap().direct_lux.to_bits();
            job.advance(cx, usize::MAX).unwrap();
            let value = job.probe_estimate(0).unwrap();
            close(value.direct_lux, 100.0, 2.0e-13, "exact direct fixture lux");
            assert_eq!(value.direct_lux.to_bits(), cached);
            if let Some(bits) = exact_bits {
                assert_eq!(cached, bits);
            } else {
                exact_bits = Some(cached);
            }
            assert_statistical_match(
                value.continuous_lux,
                reflected,
                value.estimated_standard_error_lux,
                8.0e-9,
            );
            assert!(value.estimated_standard_error_lux.unwrap() > 0.0);
            assert_eq!(
                value.illuminance_lux.to_bits(),
                (value.direct_lux + value.continuous_lux).to_bits()
            );
            assert!(value.total_error_bound_lux.is_none());
        }
    });
}

#[test]
fn g3_unselected_sources_cannot_hide_range_singularity_or_support_refusals() {
    with_cx(|cx, _| {
        let scene = receiver_scene();
        let settings = iid_settings(1, 8);
        for far in [true, false] {
            let position = if far {
                [20.0, 0.0, 1.0]
            } else {
                [0.0, 0.0, 1.0e-8]
            };
            let sources = [
                isotropic([0.0, 0.0, 3.0], 90.0, SourceSpectrum::equal_energy()),
                PhotometricSource::try_new(
                    RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], position).unwrap(),
                    1.0e-20,
                    1.0,
                    AngularDistribution::uniform_cone(0.2).unwrap(),
                    SourceSpectrum::equal_energy(),
                )
                .unwrap(),
            ];
            // The second source has exactly zero intensity toward the floor,
            // hence cannot be selected, but its positive-flux state still
            // owes the same geometric admission as the all-source estimator.
            let mut job = PathLuminanceJob::new(
                &scene,
                &sources,
                PathMedium::Vacuum,
                TransportOptions {
                    max_distance_m: 10.0,
                    ..TransportOptions::default()
                },
                settings,
                cx,
            )
            .unwrap()
            .with_point_sampling(importance(1))
            .unwrap();
            assert!(matches!(
                job.advance(cx, 1),
                Err(PhotometryError::Geometry(TracerError::FixtureConnection {
                    source: 1,
                    ..
                }))
            ));
            assert_eq!(job.completed_path_samples(), 0);
            assert!(job.pixel_estimate(0).is_none());
        }
        let sources = [
            isotropic([0.0, 0.0, 3.0], 1.0, SourceSpectrum::equal_energy()),
            isotropic([0.0, 0.0, 3.0], 1.0e-30, SourceSpectrum::equal_energy()),
        ];
        let mut unsupported = image_job(
            &scene,
            &sources,
            PathMedium::Vacuum,
            settings,
            importance(1),
            cx,
        );
        assert!(matches!(
            unsupported.advance(cx, 1),
            Err(PhotometryError::Geometry(
                TracerError::FixtureConnection { .. }
            ))
        ));
        assert_eq!(unsupported.completed_path_samples(), 0);

        let scene = empty_scene([0.0; 3], [0.0, 0.0, 1.0]);
        let sources = [
            isotropic([0.0, 0.0, 2.0], 100.0, SourceSpectrum::equal_energy()),
            isotropic([0.0, 0.0, 0.3], 1.0e-20, SourceSpectrum::equal_energy()),
        ];
        let density = [1.0];
        let medium = ScatteringGrid::try_new(
            cx,
            AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &density, [-1.0, -1.0, 0.0], [2.0, 2.0, 1.0]),
                MassAbsorptionSpectrum::constant(0.0).unwrap(),
                None,
                medium_budget(),
            )
            .unwrap(),
            MassScatteringSpectrum::constant(0.001).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        let mut singular = image_job(
            &scene,
            &sources,
            PathMedium::Scattering(&medium),
            settings,
            importance(1),
            cx,
        );
        assert!(matches!(
            singular.advance(cx, 1),
            Err(PhotometryError::Geometry(TracerError::FixtureConnection {
                source: 1,
                ..
            }))
        ));
        assert_eq!(singular.completed_path_samples(), 0);
    });
}

fn execution(workers: usize, run: u64) -> RenderExecutionConfig {
    RenderExecutionConfig::try_new(2, 2, workers, 64 << 20, RunId(run)).unwrap()
}

fn image_bits(job: &PathLuminanceJob<'_, '_>) -> Vec<Option<([u64; 3], u32, Option<u64>)>> {
    job.estimates()
        .map(|value| value.map(estimate_bits))
        .collect()
}

fn probe_bits(job: &PathIlluminanceJob<'_, '_>) -> Vec<Option<([u64; 3], u32, Option<u64>)>> {
    job.estimates()
        .map(|value| {
            value.map(|v| {
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
fn g4_g5_sampled_image_and_probe_paths_replay_across_parallel_cancel_and_resume() {
    for sampler in [Sampler::Iid, Sampler::OwenSobol] {
        with_cx(|cx, gate| {
            let mut scene = colored_receiver();
            scene.camera.half_tan = 0.3;
            scene
                .primitives
                .push(floor(5.0, 6.0, LiftedSpectrum { c: [0.0; 3] }));
            let (sources, _) = colored_rig();
            let density = [0.2, 1.0, 0.6, 0.4, 1.4, 0.3, 0.7, 1.1];
            let medium = ScatteringGrid::try_new(
                cx,
                AbsorbingGrid::try_new(
                    cx,
                    VolumeGrid::new([2; 3], &density, [-6.0, -6.0, 0.0], [6.0, 6.0, 2.5]),
                    MassAbsorptionSpectrum::table(&[(380.0, 0.02), (780.0, 0.09)]).unwrap(),
                    None,
                    medium_budget(),
                )
                .unwrap(),
                MassScatteringSpectrum::constant(0.08).unwrap(),
                PhaseFunction::henyey_greenstein(0.35).unwrap(),
            )
            .unwrap();
            let settings = Settings {
                width: 4,
                height: 3,
                sampler,
                ..settings(3, 24)
            };
            let probes = [probe_at(0.0, 0.0), probe_at(0.5, 0.25)];
            let image = || {
                image_job(
                    &scene,
                    &sources,
                    PathMedium::Scattering(&medium),
                    settings,
                    importance(4),
                    cx,
                )
            };
            let probe = || {
                probe_job(
                    &scene,
                    &sources,
                    PathMedium::Scattering(&medium),
                    &probes,
                    settings,
                    importance(4),
                    cx,
                )
            };
            let mut image_reference = image();
            let mut probe_reference = probe();
            image_reference.advance(cx, usize::MAX).unwrap();
            probe_reference.advance(cx, usize::MAX).unwrap();
            for workers in [1, 4] {
                let mut parallel_image = image();
                let mut parallel_probe = probe();
                let exec = execution(workers, 2100 + workers as u64);
                parallel_image
                    .advance_with_execution(cx, usize::MAX, &exec)
                    .unwrap();
                parallel_probe
                    .advance_with_execution(cx, usize::MAX, &exec)
                    .unwrap();
                assert_eq!(image_bits(&parallel_image), image_bits(&image_reference));
                assert_eq!(probe_bits(&parallel_probe), probe_bits(&probe_reference));
            }
            let mut mixed_image = image();
            let mut mixed_probe = probe();
            mixed_image.advance(cx, 3).unwrap();
            mixed_probe.advance(cx, 2).unwrap();
            mixed_image
                .advance_with_execution(cx, 5, &execution(2, 2110))
                .unwrap();
            mixed_probe
                .advance_with_execution(cx, 7, &execution(2, 2111))
                .unwrap();
            let image_prefix = image_bits(&mixed_image);
            let probe_prefix = probe_bits(&mixed_probe);
            gate.request();
            assert!(matches!(
                mixed_image.advance_with_execution(cx, 20, &execution(4, 2112)),
                Err(PhotometryError::Cancelled)
            ));
            assert!(matches!(
                mixed_probe.advance_with_execution(cx, 20, &execution(4, 2113)),
                Err(PhotometryError::Cancelled)
            ));
            assert_eq!(image_bits(&mixed_image), image_prefix);
            assert_eq!(probe_bits(&mixed_probe), probe_prefix);
            with_cx(|resumed, _| {
                let exec = execution(4, 2114);
                while !mixed_image.is_complete() {
                    mixed_image
                        .advance_with_execution(resumed, 17, &exec)
                        .unwrap();
                    mixed_image.advance(resumed, 3).unwrap();
                }
                while !mixed_probe.is_complete() {
                    mixed_probe.advance(resumed, 5).unwrap();
                    mixed_probe
                        .advance_with_execution(resumed, 11, &exec)
                        .unwrap();
                }
            });
            assert_eq!(image_bits(&mixed_image), image_bits(&image_reference));
            assert_eq!(probe_bits(&mixed_probe), probe_bits(&probe_reference));
        });
    }
}

/// The same 20x15 rig and 4x4 fixed-ray image as the existing 300-fixture
/// workload, now with eight real opaque baffles and a 64-cell absorption grid.
/// It reports observed error per unit time; speed alone does not qualify a
/// stochastic connection policy. Repeat this ignored test for host timing.
#[test]
#[ignore = "host-specific 300-fixture connection sampling error/throughput qualification"]
fn photometric_sampled_300_fixture_throughput_and_error() {
    with_cx(|cx, _| {
        let mut scene = receiver_scene();
        for strip in 0..8 {
            let center = -4.25 + f64::from(strip);
            scene
                .primitives
                .push(rectangle(center - 0.15, center + 0.15, -4.0, 4.0, 4.0));
        }
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
        let density = [1.0; 64];
        let medium = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([4; 3], &density, [-12.0, -10.0, 0.0], [6.0, 5.0, 2.0]),
            MassAbsorptionSpectrum::constant(0.08).unwrap(),
            None,
            AbsorbingGridBudget {
                max_grid_cells: 64,
                max_steps_per_segment: 128,
            },
        )
        .unwrap();
        let expected: f64 = (0..300)
            .map(|index| {
                let x = f64::from(index % 20) - 9.5;
                let y = f64::from(index / 20) - 7.0;
                if (0..8).any(|strip| (0.5 * x - (-4.25 + f64::from(strip))).abs() < 0.15) {
                    return 0.0;
                }
                let r2 = x * x + y * y + 64.0;
                0.5 / PI * 800.0 / (r2 * r2.sqrt()) * (-0.08 * (2.0 + r2.sqrt())).exp()
            })
            .sum();
        assert!(expected > 0.0);
        for policy in [
            PointFixtureSampling::All,
            importance(1),
            importance(4),
            importance(8),
        ] {
            let mut squared_error = 0.0;
            let mut seconds_total = 0.0;
            let mut outputs = 0;
            for seed in [173, 941, 1877] {
                let settings = Settings {
                    width: 4,
                    height: 4,
                    seed,
                    ..iid_settings(1, 512)
                };
                let mut job = image_job(
                    &scene,
                    &sources,
                    PathMedium::Absorbing(&medium),
                    settings,
                    policy,
                    cx,
                );
                let start = std::time::Instant::now();
                let paths = job.advance(cx, usize::MAX).unwrap();
                let seconds = start.elapsed().as_secs_f64();
                assert_eq!(paths, 8192);
                assert!(job.is_complete());
                let mut seed_squared_error = 0.0;
                let mut se_squared = 0.0;
                for value in job.estimates().flatten() {
                    assert_statistical_match(
                        value.luminance_cd_m2,
                        expected,
                        value.estimated_standard_error_cd_m2,
                        0.0,
                    );
                    seed_squared_error += (value.luminance_cd_m2 - expected).powi(2);
                    se_squared += value.estimated_standard_error_cd_m2.unwrap().powi(2);
                    outputs += 1;
                }
                squared_error += seed_squared_error;
                seconds_total += seconds;
                println!(
                    "photometric_point_sampling policy={policy:?} seed={seed} fixtures=300 primitives={} cells=64 paths={paths} seconds={seconds:.6} paths_per_second={:.3} reference_cd_m2={expected:.9} rmse_cd_m2={:.9} rms_estimated_se_cd_m2={:.9}",
                    scene.primitives.len(),
                    paths as f64 / seconds,
                    (seed_squared_error / 16.0).sqrt(),
                    (se_squared / 16.0).sqrt(),
                );
            }
            let mse = squared_error / f64::from(outputs);
            let mean_seconds = seconds_total / 3.0;
            println!(
                "photometric_point_sampling_summary policy={policy:?} seeds=3 pixels_per_seed=16 spp=512 mean_seconds={mean_seconds:.6} rmse_cd_m2={:.9} relative_rmse={:.9} mse_times_seconds={:.9}",
                mse.sqrt(),
                mse.sqrt() / expected,
                mse * mean_seconds,
            );
        }
    });
}

/// Paired clear-scene qualification of the same 300-fixture workload. Equal
/// spectra and a gray Lambertian make each actual contribution proportional
/// to the geometric proposal; K=1 should add no source-selection variance in
/// exact arithmetic. Spectral-packet variance is still measured, not removed.
#[test]
#[ignore = "host-specific clear 300-fixture connection sampling error/throughput qualification"]
fn photometric_sampled_300_fixture_clear_throughput_and_error() {
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
        let expected: f64 = (0..300)
            .map(|index| {
                let x = f64::from(index % 20) - 9.5;
                let y = f64::from(index / 20) - 7.0;
                let r2 = x * x + y * y + 64.0;
                0.5 / PI * 800.0 / (r2 * r2.sqrt())
            })
            .sum();
        for policy in [PointFixtureSampling::All, importance(1)] {
            let mut squared_error = 0.0;
            let mut seconds_total = 0.0;
            for seed in [173, 941, 1877] {
                let settings = Settings {
                    width: 4,
                    height: 4,
                    seed,
                    ..iid_settings(1, 512)
                };
                let mut job = image_job(&scene, &sources, PathMedium::Vacuum, settings, policy, cx);
                let start = std::time::Instant::now();
                let paths = job.advance(cx, usize::MAX).unwrap();
                let seconds = start.elapsed().as_secs_f64();
                assert_eq!(paths, 8192);
                assert!(job.is_complete());
                let mut seed_squared_error = 0.0;
                let mut se_squared = 0.0;
                for value in job.estimates().flatten() {
                    assert_statistical_match(
                        value.luminance_cd_m2,
                        expected,
                        value.estimated_standard_error_cd_m2,
                        0.0,
                    );
                    seed_squared_error += (value.luminance_cd_m2 - expected).powi(2);
                    se_squared += value.estimated_standard_error_cd_m2.unwrap().powi(2);
                }
                squared_error += seed_squared_error;
                seconds_total += seconds;
                println!(
                    "photometric_point_sampling_clear policy={policy:?} seed={seed} fixtures=300 primitives=1 cells=0 paths={paths} seconds={seconds:.6} paths_per_second={:.3} reference_cd_m2={expected:.9} rmse_cd_m2={:.9} rms_estimated_se_cd_m2={:.9}",
                    paths as f64 / seconds,
                    (seed_squared_error / 16.0).sqrt(),
                    (se_squared / 16.0).sqrt(),
                );
            }
            let mse = squared_error / 48.0;
            let mean_seconds = seconds_total / 3.0;
            println!(
                "photometric_point_sampling_clear_summary policy={policy:?} seeds=3 pixels_per_seed=16 spp=512 mean_seconds={mean_seconds:.6} rmse_cd_m2={:.9} relative_rmse={:.9} mse_times_seconds={:.9}",
                mse.sqrt(),
                mse.sqrt() / expected,
                mse * mean_seconds,
            );
        }
    });
}
