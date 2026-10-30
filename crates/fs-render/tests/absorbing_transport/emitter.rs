//! Physical volume sources illuminate surfaces and other participating cells.

use super::*;
use fs_exec::RunId;
use fs_render::motion::{ShotTimeBounds, ShutterConvention, ShutterDistribution, ShutterInterval};
use fs_render::tracer::{
    MediumExposure, MediumExposureFilm, RenderExecutionConfig, ScatteringFilm,
    render_absorbing_with_execution, render_medium_exposure, render_medium_exposure_range,
    render_medium_exposure_with_execution, render_scattering, render_scattering_range,
    render_scattering_range_with_execution, render_scattering_with_execution,
};
use fs_render::volumes::transport::scattering::{
    MassScatteringSpectrum, PhaseFunction, ScatteringGrid,
};
use fs_render::volumes::transport::timeline::{
    MediumFrame, MediumTimeline, MediumTimelineBudget, PhysicalMedium,
};

const TEMPERATURE_K: f64 = 2400.0;
const ABSORPTION: f64 = 0.7;

fn hot_grid<'field>(
    cx: &Cx<'_>,
    grid: VolumeGrid<'field>,
    temperature: &'field [f64],
    absorption: f64,
) -> AbsorbingGrid<'field> {
    AbsorbingGrid::try_new(
        cx,
        grid,
        MassAbsorptionSpectrum::constant(absorption).unwrap(),
        Some(temperature),
        budget(),
    )
    .unwrap()
}

fn source_box<'field>(
    cx: &Cx<'_>,
    density: &'field [f64],
    temperature: &'field [f64],
) -> AbsorbingGrid<'field> {
    hot_grid(
        cx,
        VolumeGrid::new([1; 3], density, [-0.5, -0.5, 1.0], [1.0; 3]),
        temperature,
        ABSORPTION,
    )
}

fn receiver_scene() -> (Scene, LiftedSpectrum) {
    let reflectance = lift_rgb([0.6; 3]);
    let mut scene = empty_scene();
    scene.camera.eye.z = 0.5;
    scene
        .primitives
        .push(quad(0.0, 10.0, Material::Lambertian { reflectance }));
    (scene, reflectance)
}

// The receiver is at the origin below [-1/2,1/2]^2 × [1,2]. Every visible
// source direction enters through z=1. Integrate that square independently
// of the production sampler/DDA: the surface-to-solid-angle Jacobian and
// receiver cosine multiply to 1/r^4, and slab intersections give chord length.
fn box_irradiance_factor(absorption: f64) -> f64 {
    let n = 160_u32;
    let spacing = 1.0 / f64::from(n);
    let mut integral = 0.0;
    for ix in 0..n {
        for iy in 0..n {
            let x = -0.5 + (f64::from(ix) + 0.5) * spacing;
            let y = -0.5 + (f64::from(iy) + 0.5) * spacing;
            let r2 = x * x + y * y + 1.0;
            let r = r2.sqrt();
            let exit = (2.0 * r).min(0.5 * r / x.abs()).min(0.5 * r / y.abs());
            let emission_fraction = -(-absorption * (exit - r)).exp_m1();
            integral += emission_fraction * spacing * spacing / (r2 * r2);
        }
    }
    integral / core::f64::consts::PI
}

fn receiver_reference(reflectance: LiftedSpectrum) -> [f64; 3] {
    let geometric_factor = box_irradiance_factor(ABSORPTION);
    xyz_of_spectrum(|lambda| {
        reflectance.eval(lambda) * planck(TEMPERATURE_K, lambda) * geometric_factor
    })
}

#[test]
fn g1_hot_volume_is_a_light_source_at_the_final_diffuse_vertex() {
    with_cx(|cx, _| {
        let density = [1.0];
        let temperature = [TEMPERATURE_K];
        let medium = source_box(cx, &density, &temperature);
        let (scene, reflectance) = receiver_scene();
        let mut settings = settings(DirectStrategy::NeeOnly);
        settings.spp = 16_384;
        settings.max_depth = 1;
        let direct = render_absorbing(&scene, &medium, cx, &settings).unwrap();
        assert_mean(&direct, receiver_reference(reflectance), 0.025);
        for strategy in [DirectStrategy::Mis, DirectStrategy::PowerMis] {
            settings.strategy = strategy;
            assert_eq!(
                render_absorbing(&scene, &medium, cx, &settings).unwrap(),
                direct,
                "a terminal source sample has no competing continuation"
            );
        }
        settings.strategy = DirectStrategy::BsdfOnly;
        assert_eq!(
            render_absorbing(&scene, &medium, cx, &settings)
                .unwrap()
                .xyz,
            vec![[0.0; 3]],
            "the camera sees only the receiver; no source is directly visible"
        );
    });
}

#[test]
fn g1_volume_source_nee_and_continuation_mis_do_not_double_count_radiance() {
    with_cx(|cx, _| {
        let density = [1.0];
        let temperature = [TEMPERATURE_K];
        let medium = source_box(cx, &density, &temperature);
        let (scene, reflectance) = receiver_scene();
        let expected = receiver_reference(reflectance);
        for strategy in [
            DirectStrategy::BsdfOnly,
            DirectStrategy::NeeOnly,
            DirectStrategy::Mis,
            DirectStrategy::PowerMis,
        ] {
            let mut settings = settings(strategy);
            settings.spp = 32_768;
            settings.max_depth = 2;
            assert_mean(
                &render_absorbing(&scene, &medium, cx, &settings).unwrap(),
                expected,
                0.03,
            );
        }
    });
}

#[test]
fn g1_hidden_cells_remain_in_source_pdf_while_opaque_visibility_clips_radiance() {
    with_cx(|cx, _| {
        let density = [1.0, 1.0];
        let temperature = [TEMPERATURE_K; 2];
        let medium = hot_grid(
            cx,
            VolumeGrid::new([1, 1, 2], &density, [-0.5, -0.5, 1.0], [1.0; 3]),
            &temperature,
            ABSORPTION,
        );
        let black = Material::Lambertian {
            reflectance: lift_rgb([0.0; 3]),
        };
        let (mut scene, reflectance) = receiver_scene();
        // The z=2 plane hides the second cell but leaves the complete first
        // cell visible. Samples selected in the hidden cell still propose
        // directions through the visible source. Clipping the direction PDF
        // at this plane would spuriously brighten those directions.
        scene.primitives.push(quad(2.0, 10.0, black));
        let mut settings = settings(DirectStrategy::Mis);
        settings.spp = 32_768;
        settings.max_depth = 1;
        assert_mean(
            &render_absorbing(&scene, &medium, cx, &settings).unwrap(),
            receiver_reference(reflectance),
            0.03,
        );
        let (mut blocked, _) = receiver_scene();
        blocked.primitives.push(quad(0.75, 10.0, black));
        assert_eq!(
            render_absorbing(&blocked, &medium, cx, &settings)
                .unwrap()
                .xyz,
            vec![[0.0; 3]],
            "an opaque surface before the first emitting cell blocks its light"
        );
    });
}

// Independent Neumann terms for a homogeneous [-1,1]^3 source. The camera
// starts at the center and looks along -z. Average the analytic source along
// each direction at a collision, then integrate sigma_s exp(-sigma_t t) over
// the one-metre camera segment. This is the direct term plus one scattering.
fn single_scatter_cube_factor(absorption: f64, scattering: f64) -> f64 {
    let extinction = absorption + scattering;
    let source_ratio = absorption / extinction;
    let mut result = source_ratio * -(-extinction).exp_m1();
    let (nt, nmu, nphi) = (80_u32, 64_u32, 128_u32);
    let dt = 1.0 / f64::from(nt);
    for it in 0..nt {
        let t = (f64::from(it) + 0.5) * dt;
        let mut angular_mean = 0.0;
        for imu in 0..nmu {
            let mu = -1.0 + 2.0 * (f64::from(imu) + 0.5) / f64::from(nmu);
            let transverse = (1.0 - mu * mu).sqrt();
            for iphi in 0..nphi {
                let phi = core::f64::consts::TAU * (f64::from(iphi) + 0.5) / f64::from(nphi);
                let x = transverse * phi.cos();
                let y = transverse * phi.sin();
                let z_exit = if mu > 0.0 {
                    (1.0 + t) / mu
                } else {
                    (t - 1.0) / mu
                };
                let length = (1.0 / x.abs()).min(1.0 / y.abs()).min(z_exit);
                angular_mean += -(-extinction * length).exp_m1();
            }
        }
        angular_mean /= f64::from(nmu) * f64::from(nphi);
        result += dt * scattering * (-extinction * t).exp() * source_ratio * angular_mean;
    }
    result
}

#[test]
fn g1_phase_source_nee_matches_independent_single_scattering_integral() {
    with_cx(|cx, _| {
        let density = [1.0];
        let temperature = [TEMPERATURE_K];
        let (absorption, scattering) = (0.3, 0.7);
        let medium = ScatteringGrid::try_new(
            cx,
            hot_grid(
                cx,
                VolumeGrid::new([1; 3], &density, [-1.0; 3], [2.0; 3]),
                &temperature,
                absorption,
            ),
            MassScatteringSpectrum::constant(scattering).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        let mut scene = empty_scene();
        scene.camera.eye = Point3::new(0.0, 0.0, 0.0);
        let factor = single_scatter_cube_factor(absorption, scattering);
        let expected = xyz_of_spectrum(|lambda| planck(TEMPERATURE_K, lambda) * factor);
        for (strategy, max_depth) in [(DirectStrategy::NeeOnly, 1), (DirectStrategy::BsdfOnly, 2)] {
            let mut settings = settings(strategy);
            settings.max_depth = max_depth;
            settings.spp = 32_768;
            assert_mean(
                &render_scattering(&scene, &medium, cx, &settings).unwrap(),
                expected,
                0.03,
            );
        }
    });
}

fn execution(width: u32, height: u32, workers: usize) -> RenderExecutionConfig {
    RenderExecutionConfig::try_new(width, height, workers, 8 << 20, RunId(971)).unwrap()
}

#[test]
fn g3_volume_source_sampling_replays_across_tiles_workers_and_sample_ranges() {
    with_cx(|cx, _| {
        let density = [1.0];
        let temperature = [TEMPERATURE_K];
        let medium = ScatteringGrid::try_new(
            cx,
            source_box(cx, &density, &temperature),
            MassScatteringSpectrum::constant(0.5).unwrap(),
            PhaseFunction::Rayleigh,
        )
        .unwrap();
        let (mut scene, _) = receiver_scene();
        scene.camera.half_tan = 0.25;
        for sampler in [Sampler::Iid, Sampler::OwenSobolFullPath] {
            let mut settings = settings(DirectStrategy::PowerMis);
            settings.width = 3;
            settings.height = 2;
            settings.spp = 32;
            settings.max_depth = 5;
            settings.sampler = sampler;
            let serial = render_scattering(&scene, &medium, cx, &settings).unwrap();
            assert!(serial.xyz.iter().all(|xyz| xyz[1] > 0.0));
            for config in [execution(1, 1, 1), execution(2, 2, 3)] {
                assert_eq!(
                    render_scattering_with_execution(&scene, &medium, cx, &settings, &config)
                        .unwrap()
                        .film,
                    serial
                );
            }
            let absorption = medium.absorption_grid();
            assert_eq!(
                render_absorbing_with_execution(
                    &scene,
                    absorption,
                    cx,
                    &settings,
                    &execution(2, 1, 3)
                )
                .unwrap()
                .film,
                render_absorbing(&scene, absorption, cx, &settings).unwrap()
            );
            let mut progressive = ScatteringFilm::try_new(&medium, &settings).unwrap();
            render_scattering_range(&scene, &medium, cx, &settings, &mut progressive, 0, 11)
                .unwrap();
            render_scattering_range_with_execution(
                &scene,
                &medium,
                cx,
                &settings,
                &mut progressive,
                11,
                32,
                &execution(2, 1, 3),
            )
            .unwrap();
            assert_eq!(*progressive.film(), serial);
        }
    });
}

#[test]
fn g3_exposure_samples_the_source_distribution_of_its_selected_physical_state() {
    with_cx(|cx, _| {
        let density = [1.0];
        let temperatures = [[1800.0], [2600.0]];
        let first = source_box(cx, &density, &temperatures[0]);
        let second = ScatteringGrid::try_new(
            cx,
            source_box(cx, &density, &temperatures[1]),
            MassScatteringSpectrum::constant(0.4).unwrap(),
            PhaseFunction::Isotropic,
        )
        .unwrap();
        let frames = [
            MediumFrame {
                start_s: 0.0,
                end_s: 0.5,
                medium: PhysicalMedium::Absorbing(&first),
            },
            MediumFrame {
                start_s: 0.5,
                end_s: 1.0,
                medium: PhysicalMedium::Scattering(&second),
            },
        ];
        let history =
            MediumTimeline::try_new(cx, &frames, MediumTimelineBudget { max_frames: 2 }).unwrap();
        let shutter = ShutterInterval::resolve(
            0.0,
            1.0,
            ShutterConvention::FrontLoaded,
            ShutterDistribution::StratifiedCounterV1 { strata: 16 },
            ShotTimeBounds::try_new(0.0, 1.0).unwrap(),
        )
        .unwrap();
        let exposure = MediumExposure::try_new(cx, &history, shutter).unwrap();
        let (scene, _) = receiver_scene();
        let mut settings = settings(DirectStrategy::Mis);
        settings.spp = 64;
        settings.max_depth = 4;
        let mut first_film = AbsorbingFilm::try_new(&first, &settings).unwrap();
        let mut second_film = ScatteringFilm::try_new(&second, &settings).unwrap();
        let mut actual = MediumExposureFilm::try_new(&exposure, &settings).unwrap();
        let mut expected = [0.0; 3];
        let mut counts = [0; 2];
        for sample in 0..settings.spp {
            let before = [first_film.film().xyz[0], second_film.film().xyz[0]];
            render_absorbing_range(
                &scene,
                &first,
                cx,
                &settings,
                &mut first_film,
                sample,
                sample + 1,
            )
            .unwrap();
            render_scattering_range(
                &scene,
                &second,
                cx,
                &settings,
                &mut second_film,
                sample,
                sample + 1,
            )
            .unwrap();
            let after = [first_film.film().xyz[0], second_film.film().xyz[0]];
            let time =
                shutter.time_at(shutter.sample_for_stream(settings.seed, 0, u64::from(sample)));
            let state = usize::from(time >= 0.5);
            counts[state] += 1;
            render_medium_exposure_range(
                &scene,
                &exposure,
                cx,
                &settings,
                &mut actual,
                sample,
                sample + 1,
            )
            .unwrap();
            for channel in 0..3 {
                expected[channel] += after[state][channel] - before[state][channel];
                assert!(
                    (actual.film().xyz[0][channel] - expected[channel]).abs()
                        <= 2.0e-12 * expected[channel].abs().max(1.0),
                    "source snapshot {state}, sample {sample}, channel {channel}"
                );
            }
        }
        assert_eq!(counts, [32, 32]);
        assert!(expected[1] > 0.0);
        assert_eq!(
            *actual.film(),
            render_medium_exposure(&scene, &exposure, cx, &settings).unwrap()
        );
        assert_eq!(
            *actual.film(),
            render_medium_exposure_with_execution(
                &scene,
                &exposure,
                cx,
                &settings,
                &execution(1, 1, 3)
            )
            .unwrap()
            .film
        );
    });
}
