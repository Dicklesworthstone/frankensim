//! Original finite cells, native flights and reciprocal thermal reception.

use super::*;
use crate::charts::Ray;
use crate::polarization::dielectric::scene::curved::camera::participating::CurvedRayleighCameraScene;
use crate::polarization::dielectric::scene::curved::camera::{
    CameraAnalyzer, CurvedCameraScene, EmissionSides, IsothermalVolumeEmitter, RadianceSpectrum,
    SurfaceEmitter,
};
use crate::thermal_surface::planck_radiance_w_m2_sr_nm;
use crate::volumes::VolumeGrid;
use crate::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};

fn grid_budget() -> AbsorbingGridBudget {
    AbsorbingGridBudget {
        max_grid_cells: 16,
        max_steps_per_segment: 16,
    }
}

fn compare_stokes(actual: Stokes, source: Stokes, factor: f64) {
    for (a, s) in [actual.i, actual.q, actual.u, actual.v]
        .into_iter()
        .zip([source.i, source.q, source.u, source.v])
    {
        close(a, s * factor, 2e-10);
    }
}

#[test]
fn finite_absorbing_cells_preserve_full_stokes_on_native_segments_and_escape_tails() {
    with_cx(|cx, _| {
        let air = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let density = [2.0, 3.0];
        let temperature = [1400.0, 2000.0];
        let input = Stokes {
            i: 2.0,
            q: 0.5,
            u: -0.4,
            v: 0.2,
        };
        for (surface_z, grid_z, escaped) in [(3.0, 0.5, false), (-3.0, 10.0, true)] {
            let geometry = InstanceScene::try_new(vec![plane(9, 2, surface_z)]).unwrap();
            let scene =
                CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector(9)], 100_000)
                    .unwrap();
            let absorber = AbsorbingGrid::try_new(
                cx,
                VolumeGrid::new([1, 1, 2], &density, [-1.0, -1.0, grid_z], [2.0, 2.0, 0.4]),
                MassAbsorptionSpectrum::constant(0.7).unwrap(),
                Some(&temperature),
                grid_budget(),
            )
            .unwrap();
            let grids = [Some(&absorber)];
            let laws = [None];
            let transport = CurvedRayleighScene::try_new(cx, &scene, &laws)
                .unwrap()
                .with_absorbing_grids(cx, &grids, 2)
                .unwrap();
            assert!(core::ptr::eq(
                transport.absorbing_grid(0).unwrap(),
                &absorber
            ));
            let beam = scene
                .launch_beam(Point3::new(0.0, 0.0, 0.0), frame(), 550.0, 0, input)
                .unwrap();
            let result = transport
                .trace_with_samples(
                    cx,
                    beam,
                    0,
                    settings(),
                    |_| panic!("detector has no continuation"),
                    |_| panic!("absorption never samples a scattering flight"),
                )
                .unwrap();
            let path = result.path();
            assert_eq!(
                path.termination(),
                if escaped {
                    PathTermination::Escaped
                } else {
                    PathTermination::Detected
                }
            );
            close(path.optical_depth(), 1.4, 2e-12);
            compare_stokes(path.beam().power(), input, (-1.4_f64).exp());
            assert_eq!(path.absorption_cell_visits(), 2);
            assert_eq!(
                path.reduced_volume_emission, 0.0,
                "source-power paths do not add radiance"
            );
            assert_eq!(result.free_flight_draws(), 0);
            if escaped {
                assert_eq!(path.traveled_m(), 0.0);
            }
        }
    });
}

#[test]
fn finite_absorbing_lte_uses_original_temperature_order_across_null_flights_once() {
    with_cx(|cx, _| {
        let air = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let geometry = InstanceScene::try_new(vec![plane(9, 2, 3.0)]).unwrap();
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![SurfaceBinding {
                object_id: 9,
                interaction: Interaction::Opaque,
            }],
            100_000,
        )
        .unwrap();
        let density = [1.0, 1.0];
        let temperature = [1400.0, 2000.0];
        let absorber = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 2], &density, [-1.0, -1.0, 0.0], [2.0, 2.0, 0.5]),
            MassAbsorptionSpectrum::constant(1.0).unwrap(),
            Some(&temperature),
            grid_budget(),
        )
        .unwrap();
        let scattering_density = [0.0, 1.0];
        let scattering_grid = VolumeGrid::new(
            [1, 1, 2],
            &scattering_density,
            [-1.0, -1.0, 0.0],
            [2.0, 2.0, 0.5],
        );
        let mass = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let laws = [Some(
            CurvedRayleighMedium::try_new_volume_grid(&mass, &scattering_grid, &phase, 2).unwrap(),
        )];
        let grids = [Some(&absorber)];
        let transport = CurvedRayleighScene::try_new(cx, &scene, &laws)
            .unwrap()
            .with_absorbing_grids(cx, &grids, 4)
            .unwrap();
        let emitters = [SurfaceEmitter::new(
            9,
            RadianceSpectrum::try_constant(5.0).unwrap(),
            EmissionSides::Both,
        )];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let owner = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let ray = Ray {
            origin: Point3::new(0.0, 0.0, 0.0),
            dir: Vec3::new(0.0, 0.0, 1.0),
        };
        let flights = |ordinal| match ordinal {
            0 => [1.0 - (-0.25_f64).exp(), 0.5, 0.0, 0.0],
            1 => [1.0 - 1e-12, 0.5, 0.0, 0.0],
            _ => panic!("one null collision and one original finite-grid exit"),
        };
        let native = absorber
            .transport(cx, ray.origin, ray.dir, 3.0, &[550.0])
            .unwrap();
        let t = (-0.5_f64).exp();
        let thermal = (1.0 - t)
            * (planck_radiance_w_m2_sr_nm(1400.0, 550.0).unwrap()
                + t * planck_radiance_w_m2_sr_nm(2000.0, 550.0).unwrap());
        close(native.emission_w_m2_sr_nm[0], thermal, 2e-12);
        let expected = 5.0 * t * t + thermal;
        let measured = owner
            .sample_ray_with_tracking_samples(
                cx,
                ray,
                550.0,
                0,
                settings(),
                |_| panic!("opaque source"),
                flights,
            )
            .unwrap();
        close(
            measured.radiance_sample().radiance_w_m2_sr_nm(),
            expected,
            2e-10,
        );
        assert_eq!(measured.volume_events(), 0);
        assert_eq!(measured.free_flight_draws(), 2);
        assert_eq!(
            measured.radiance_sample().path().absorption_cell_visits(),
            3
        );
        for analyzer in [
            CameraAnalyzer::try_linear(Vec3::new(0.0, 1.0, 0.0), 0.37).unwrap(),
            CameraAnalyzer::try_circular(Vec3::new(0.0, 1.0, 0.0), true).unwrap(),
        ] {
            let analyzed = owner
                .sample_ray_analyzed_with_tracking_samples(
                    cx,
                    ray,
                    550.0,
                    0,
                    analyzer,
                    settings(),
                    |_| panic!("opaque source"),
                    flights,
                )
                .unwrap();
            close(
                analyzed.radiance_sample().radiance_w_m2_sr_nm(),
                expected * 0.5,
                2e-10,
            );
        }
    });
}

#[test]
fn finite_absorbing_real_rayleigh_collision_keeps_native_mueller_and_only_true_absorption() {
    with_cx(|cx, _| {
        let air = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let geometry = InstanceScene::try_new(vec![plane(9, 2, 3.0)]).unwrap();
        let scene =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector(9)], 100_000)
                .unwrap();
        let density = [1.0];
        let absorber = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 1], &density, [-1.0, -1.0, 0.0], [2.0; 3]),
            MassAbsorptionSpectrum::constant(0.4).unwrap(),
            None,
            grid_budget(),
        )
        .unwrap();
        let mass = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let laws = [Some(
            CurvedRayleighMedium::try_new(&mass, 1.0, &phase).unwrap(),
        )];
        let grids = [Some(&absorber)];
        let owner = CurvedRayleighScene::try_new(cx, &scene, &laws)
            .unwrap()
            .with_absorbing_grids(cx, &grids, 1)
            .unwrap();
        let beam = scene
            .launch_beam(
                Point3::new(0.0, 0.0, 0.0),
                frame(),
                550.0,
                0,
                Stokes {
                    i: 1.0,
                    q: 0.2,
                    u: -0.1,
                    v: 0.3,
                },
            )
            .unwrap();
        let angles = [0.5 + f64::EPSILON, 0.0];
        let expected = phase.sample_polarized(beam, angles).unwrap();
        let observed = owner
            .trace_with_samples(
                cx,
                beam,
                0,
                CurvedRayleighTraceSettings {
                    max_volume_events: 1,
                    ..settings()
                },
                |_| panic!("collision precedes surface"),
                |_| [1.0 - (-0.25_f64).exp(), angles[0], angles[1]],
            )
            .unwrap();
        assert_eq!(observed.volume_events(), 1);
        assert_eq!(observed.path().termination(), PathTermination::VolumeLimit);
        close(observed.path().traveled_m(), 0.25, 2e-10);
        close(observed.path().optical_depth(), 0.1, 2e-10);
        compare_stokes(
            observed.path().beam().power(),
            expected.beam.power(),
            (-0.1_f64).exp(),
        );
    });
}

#[test]
fn finite_absorbing_admission_budget_spectrum_and_duplicate_temperature_refuse() {
    with_cx(|cx, gate| {
        let air = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let geometry = InstanceScene::try_new(vec![plane(9, 2, 3.0)]).unwrap();
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![SurfaceBinding {
                object_id: 9,
                interaction: Interaction::Opaque,
            }],
            100_000,
        )
        .unwrap();
        let density = [1.0, 1.0];
        let absorber = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 2], &density, [-1.0, -1.0, 0.0], [2.0, 2.0, 0.5]),
            MassAbsorptionSpectrum::source_table(
                &[(450.0, 0.4), (650.0, 0.4)],
                ContentHash([91; 32]),
            )
            .unwrap(),
            None,
            grid_budget(),
        )
        .unwrap();
        let grids = [Some(&absorber)];
        let laws = [None];
        let bind = |cap| {
            CurvedRayleighScene::try_new(cx, &scene, &laws)
                .unwrap()
                .with_absorbing_grids(cx, &grids, cap)
                .unwrap()
        };
        let input = |wavelength| {
            scene
                .launch_beam(
                    Point3::new(0.0, 0.0, 0.0),
                    frame(),
                    wavelength,
                    0,
                    Stokes::unpolarized(1.0),
                )
                .unwrap()
        };
        assert!(
            bind(1)
                .trace(cx, input(550.0), 0, settings(), 1, 2)
                .is_err()
        );
        assert!(
            bind(2)
                .trace(cx, input(700.0), 0, settings(), 1, 2)
                .is_err()
        );
        let thermal = IsothermalVolumeEmitter::try_new(0, 1700.0, ContentHash([90; 32])).unwrap();
        let camera = CurvedCameraScene::try_new(cx, &scene, &[], None)
            .unwrap()
            .with_isothermal_volume_emitters(cx, &[thermal])
            .unwrap();
        let transport = bind(2);
        assert!(CurvedRayleighCameraScene::try_new(cx, &camera, &transport).is_err());
        for material in [glass(1.2, 0.0), glass(1.0, 0.2)] {
            let hosts = [SceneMedium::Homogeneous(&material)];
            let original =
                CurvedPolarizedScene::try_new(cx, &geometry, &hosts, vec![detector(9)], 100_000)
                    .unwrap();
            let invalid = CurvedRayleighScene::try_new(cx, &original, &laws)
                .unwrap()
                .with_absorbing_grids(cx, &grids, 2)
                .unwrap();
            let beam = original
                .launch_beam(
                    Point3::new(0.0, 0.0, 0.0),
                    frame(),
                    550.0,
                    0,
                    Stokes::unpolarized(1.0),
                )
                .unwrap();
            assert!(invalid.trace(cx, beam, 0, settings(), 1, 2).is_err());
        }
        gate.request();
        assert!(
            transport
                .trace(cx, input(550.0), 0, settings(), 1, 2)
                .is_err()
        );
    });
}
