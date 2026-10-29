//! Physical fields share the existing executor's sample order and transaction.

use super::*;
use fs_exec::RunId;
use fs_render::tracer::{
    RenderExecutionConfig, RenderExecutionError, render_absorbing_range_with_execution,
    render_absorbing_with_execution, render_scattering_range_with_execution,
    render_scattering_with_execution, render_with_execution,
};

const MEMORY_BYTES: u64 = 8 << 20;

fn policy(width: u32, height: u32, workers: usize, run: u64) -> RenderExecutionConfig {
    RenderExecutionConfig::try_new(width, height, workers, MEMORY_BYTES, RunId(run)).unwrap()
}

fn assert_film_bits_equal(actual: &Film, expected: &Film) {
    assert_eq!(
        (
            actual.width,
            actual.height,
            actual.spp_done,
            actual.time_mode
        ),
        (
            expected.width,
            expected.height,
            expected.spp_done,
            expected.time_mode
        )
    );
    assert_eq!(actual.xyz.len(), expected.xyz.len());
    for (pixel, (actual, expected)) in actual.xyz.iter().zip(&expected.xyz).enumerate() {
        for channel in 0..3 {
            assert_eq!(
                actual[channel].to_bits(),
                expected[channel].to_bits(),
                "pixel {pixel}, channel {channel}"
            );
        }
    }
}

#[test]
fn g3_parallel_scattering_matches_serial_across_tiles_and_progressive_splits() {
    with_cx(|cx, _| {
        let density = [1.0];
        let medium = physical_grid(
            cx,
            &density,
            [-1.0; 3],
            [2.0; 3],
            0.2,
            0.7,
            PhaseFunction::HenyeyGreenstein { g: 0.4 },
        );
        let (mut scene, _) = environment_scene();
        scene.camera.eye = Point3::new(0.0, 0.0, 0.0);
        scene.camera.half_tan = 0.4;
        for sampler in [Sampler::Iid, Sampler::OwenSobolFullPath] {
            let mut settings = settings(DirectStrategy::PowerMis);
            settings.width = 5;
            settings.height = 3;
            settings.spp = 24;
            settings.max_depth = 8;
            settings.sampler = sampler;
            let serial = render_scattering(&scene, &medium, cx, &settings).unwrap();
            let film_bytes = u64::from(settings.width) * u64::from(settings.height) * 24;
            for (tw, th, workers, run) in [(1, 1, 1, 1), (3, 2, 3, 2)] {
                let execution = policy(tw, th, workers, run);
                let output =
                    render_scattering_with_execution(&scene, &medium, cx, &settings, &execution)
                        .unwrap();
                assert_film_bits_equal(&output.film, &serial);
                assert_eq!(output.report.retained_film_bytes, 0);
                assert_eq!(output.report.staging_film_bytes, film_bytes);
                assert_eq!(output.report.memory.used_bytes, 0);
                assert!(output.report.memory.peak_bytes <= MEMORY_BYTES);
            }
            let mut progressive = ScatteringFilm::try_new(&medium, &settings).unwrap();
            render_scattering_range(&scene, &medium, cx, &settings, &mut progressive, 0, 7)
                .unwrap();
            let report = render_scattering_range_with_execution(
                &scene,
                &medium,
                cx,
                &settings,
                &mut progressive,
                7,
                17,
                &policy(3, 2, 3, 3),
            )
            .unwrap();
            assert_eq!(report.retained_film_bytes, film_bytes);
            assert!(report.memory.peak_bytes >= 2 * film_bytes);
            render_scattering_range_with_execution(
                &scene,
                &medium,
                cx,
                &settings,
                &mut progressive,
                17,
                24,
                &policy(2, 1, 2, 4),
            )
            .unwrap();
            assert_film_bits_equal(progressive.film(), &serial);
        }
    });
}

#[test]
fn g3_parallel_absorption_vacuum_and_emission_share_serial_pixels() {
    with_cx(|cx, _| {
        let density = [1.0];
        let (scene, _) = environment_scene();
        let mut settings = settings(DirectStrategy::Mis);
        settings.width = 4;
        settings.height = 3;
        settings.spp = 16;
        let execution = policy(2, 2, 3, 5);
        for coefficient in [0.0, 0.6] {
            let medium = AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1; 3], &density, [-1.0; 3], [2.0; 3]),
                MassAbsorptionSpectrum::constant(coefficient).unwrap(),
                None,
                budget(),
            )
            .unwrap();
            let serial = render_absorbing(&scene, &medium, cx, &settings).unwrap();
            let output =
                render_absorbing_with_execution(&scene, &medium, cx, &settings, &execution)
                    .unwrap();
            assert_film_bits_equal(&output.film, &serial);
            if coefficient == 0.0 {
                assert_film_bits_equal(
                    &output.film,
                    &render_with_execution(&scene, cx, &settings, &execution)
                        .unwrap()
                        .film,
                );
            }
            let mut progressive = AbsorbingFilm::try_new(&medium, &settings).unwrap();
            render_absorbing_range(&scene, &medium, cx, &settings, &mut progressive, 0, 5).unwrap();
            render_absorbing_range_with_execution(
                &scene,
                &medium,
                cx,
                &settings,
                &mut progressive,
                5,
                16,
                &execution,
            )
            .unwrap();
            assert_film_bits_equal(progressive.film(), &serial);
        }
        // Medium emission alone must pass the same scene admission, without
        // inventing a surface or environment source for the parallel path.
        let hot = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1; 3], &density, [-1.0; 3], [2.0; 3]),
            MassAbsorptionSpectrum::constant(0.6).unwrap(),
            Some(&[1800.0]),
            budget(),
        )
        .unwrap();
        let scene = empty_scene();
        assert_film_bits_equal(
            &render_absorbing_with_execution(&scene, &hot, cx, &settings, &execution)
                .unwrap()
                .film,
            &render_absorbing(&scene, &hot, cx, &settings).unwrap(),
        );
    });
}

#[test]
fn g4_parallel_medium_memory_traversal_and_cancel_refusals_are_atomic() {
    with_cx(|cx, gate| {
        let density = [1.0; 2];
        let temperatures = [1800.0; 2];
        let make = |steps| {
            ScatteringGrid::try_new(
                cx,
                AbsorbingGrid::try_new(
                    cx,
                    VolumeGrid::new([1, 1, 2], &density, [-1.0, -1.0, 0.0], [2.0, 2.0, 1.0]),
                    MassAbsorptionSpectrum::constant(0.2).unwrap(),
                    Some(&temperatures),
                    AbsorbingGridBudget {
                        max_grid_cells: 2,
                        max_steps_per_segment: steps,
                    },
                )
                .unwrap(),
                MassScatteringSpectrum::constant(0.7).unwrap(),
                PhaseFunction::Isotropic,
            )
            .unwrap()
        };
        let medium = make(8);
        let restricted = make(1);
        assert_eq!(medium.content_identity(), restricted.content_identity());
        let (scene, _) = environment_scene();
        let mut settings = settings(DirectStrategy::PowerMis);
        settings.width = 3;
        settings.height = 2;
        settings.spp = 8;
        let mut film = ScatteringFilm::try_new(&medium, &settings).unwrap();
        render_scattering_range(&scene, &medium, cx, &settings, &mut film, 0, 2).unwrap();
        let before = film.clone();
        let tiny = RenderExecutionConfig::try_new(2, 1, 2, 1, RunId(6)).unwrap();
        assert!(matches!(
            render_scattering_range_with_execution(
                &scene, &medium, cx, &settings, &mut film, 2, 8, &tiny,
            ),
            Err(RenderExecutionError::Memory(_))
        ));
        assert_eq!(film, before);
        let execution = policy(2, 1, 2, 7);
        assert_eq!(
            render_scattering_range_with_execution(
                &scene,
                &restricted,
                cx,
                &settings,
                &mut film,
                2,
                8,
                &execution,
            ),
            Err(RenderExecutionError::Tracer(TracerError::Volume(
                VolumeTransportError::Budget
            )))
        );
        assert_eq!(film, before);
        gate.request();
        assert_eq!(
            render_scattering_range_with_execution(
                &scene, &medium, cx, &settings, &mut film, 2, 8, &execution,
            ),
            Err(RenderExecutionError::Tracer(TracerError::Cancelled))
        );
        assert_eq!(film, before);
    });
}
