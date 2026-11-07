//! Original grid support, physical collision phases, and complementary MIS.

use super::*;
use super::super::participating::{
    CurvedParticipatingEnvironmentCameraScene, ParticipatingEnvironmentCameraSettings,
};
use super::super::weights::positive_product_with_log_survival;
use crate::polarization::dielectric::scene::curved::camera::participating::CurvedRayleighCameraScene;
use crate::polarization::dielectric::scene::curved::participating::{
    CurvedRayleighMedium, CurvedRayleighScene, CurvedRayleighTraceSettings,
};
use crate::volumes::VolumeGrid;
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};
use crate::volumes::transport::scattering::mie::{MieBudget, MieSphere};

fn participating_settings() -> ParticipatingEnvironmentCameraSettings {
    ParticipatingEnvironmentCameraSettings {
        environment: settings(),
        max_volume_events: 8,
        max_source_cells: 64,
    }
}

fn native_settings() -> CurvedRayleighTraceSettings {
    CurvedRayleighTraceSettings {
        surfaces: settings().surfaces,
        max_volume_events: 8,
    }
}

#[test]
fn absent_and_zero_scattering_preserve_every_surface_environment_sample_bit() {
    with_cx(|cx| {
        let map = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap();
        let source = EnvironmentRadiance::try_new(&map, 3.0).unwrap();
        let geometry = InstanceScene::try_new(vec![patch(9, 0.0, 1.0, 0.3)]).unwrap();
        let air = material(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let diffuse = [SpectralReflectance::constant(0.6).unwrap()];
        let native = CurvedPolarizedScene::try_new_with_materials(
            cx, &geometry, &media, &[], &[], &diffuse,
            vec![SurfaceBinding { object_id: 9, interaction: Interaction::DepolarizingDiffuse {
                medium: 0, reflectance: 0,
            }}],
            100_000,
        ).unwrap();
        let sources = [SurfaceEmitter::new(9, RadianceSpectrum::try_constant(5.0).unwrap(),
            EmissionSides::Front)];
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, Some(source)).unwrap();
        let ordinary = CurvedEnvironmentCameraScene::try_new(cx, &camera).unwrap()
            .sample_ray_with_samples(cx, camera_ray(), 550.0, 0, settings(),
                |_| [0.0; 3], |_| [0.5, 0.75]).unwrap();
        let zero = MassScatteringSpectrum::constant(0.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        for law in [None, Some(CurvedRayleighMedium::try_new(&zero, 2.0, &phase).unwrap())] {
            let volumes = [law];
            let transport = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
            let participating = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
            let owner = CurvedParticipatingEnvironmentCameraScene::try_new(cx, &participating).unwrap();
            let observed = owner.sample_ray_with_tracking_samples(
                cx, camera_ray(), 550.0, 0,
                ParticipatingEnvironmentCameraSettings { max_source_cells: 0, ..participating_settings() },
                |_| [0.0; 3],
                |_| panic!("zero scattering consumes no flight variate"),
                |_| [0.5, 0.75],
                |_| panic!("no physical collision"),
            ).unwrap();
            assert_eq!(observed.environment_sample(), &ordinary);
            assert_eq!(observed.volume_events(), 0);
            assert_eq!(observed.free_flight_draws(), 0);
            assert_eq!(observed.source_cell_visits(), 0);
        }
    });
}

#[test]
fn native_grid_surface_nee_integrates_column_once_and_null_candidates_keep_mis_anchor() {
    with_cx(|cx| {
        let map = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap();
        let source = EnvironmentRadiance::try_new(&map, 3.0).unwrap();
        let geometry = InstanceScene::try_new(vec![patch(9, 0.0, 1.0, 0.3)]).unwrap();
        let air = material(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let diffuse = [SpectralReflectance::constant(0.6).unwrap()];
        let native = CurvedPolarizedScene::try_new_with_materials(
            cx, &geometry, &media, &[], &[], &diffuse,
            vec![SurfaceBinding { object_id: 9, interaction: Interaction::DepolarizingDiffuse {
                medium: 0, reflectance: 0,
            }}],
            100_000,
        ).unwrap();
        let camera = CurvedCameraScene::try_new(cx, &native, &[], Some(source)).unwrap();
        let density = [0.0, 1.5];
        let grid = VolumeGrid::new([1, 1, 2], &density, [-1.0, -1.0, 0.1], [2.0, 2.0, 0.4]);
        let mass = MassScatteringSpectrum::constant(0.7).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [Some(CurvedRayleighMedium::try_new_volume_grid(
            &mass, &grid, &phase, 2,
        ).unwrap())];
        let transport = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        let participating = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let owner = CurvedParticipatingEnvironmentCameraScene::try_new(cx, &participating).unwrap();
        let flight = |ordinal| match ordinal {
            0 | 2 => [1.0 - 1e-12, 0.5, 0.0, 0.0],
            // After diffuse reflection the original majorant flight enters at
            // z=.9 and travels .6 m to the zero-density lower cell at z=.3.
            1 => [1.0 - (-1.05_f64 * 0.6).exp(), 0.5, 0.0, 0.0],
            _ => panic!("grid exit must not generate an infinite null tail"),
        };
        let observed = owner.sample_ray_with_tracking_samples(
            cx, camera_ray(), 550.0, 0, participating_settings(),
            |_| [0.0; 3], flight,
            |ordinal| { assert_eq!(ordinal, 0); [0.5, 0.75] },
            |_| panic!("null candidate is not a physical volume vertex"),
        ).unwrap();
        let baseline = participating.sample_ray_with_tracking_samples(
            cx, camera_ray(), 550.0, 0, native_settings(), |_| [0.0; 3], flight,
        ).unwrap();
        let le = source.evaluate(Z.scale(-1.0), 550.0).unwrap();
        let direct = 0.6 * le * (4.0 / 17.0) * (-0.7_f64 * 1.5 * 0.4).exp();
        near(observed.environment_radiance_w_m2_sr_nm(), direct, 4e-6);
        near(observed.environment_endpoint_weight(), 16.0 / 17.0, 2e-13);
        near(observed.radiance_sample().radiance_w_m2_sr_nm(),
            direct + baseline.radiance_sample().radiance_w_m2_sr_nm() * 16.0 / 17.0, 4e-6);
        assert_eq!(observed.radiance_sample().path(), baseline.radiance_sample().path());
        assert_eq!(observed.volume_events(), 0);
        assert_eq!(observed.free_flight_draws(), 3);
        assert_eq!(observed.source_cell_visits(), 2);
        assert_eq!(observed.environment_sample().environment_samples(), 1);
        assert!(owner.sample_ray_with_tracking_samples(
            cx, camera_ray(), 550.0, 0,
            ParticipatingEnvironmentCameraSettings { max_source_cells: 1, ..participating_settings() },
            |_| [0.0; 3], flight, |_| [0.5, 0.75], |_| [0.5, 0.0],
        ).is_err());
    });
}

#[test]
fn accepted_rayleigh_and_spectral_mie_vertices_use_physical_source_stokes_and_native_phase_pdf() {
    with_cx(|cx| {
        let sphere = MieSphere::try_new(80e-9, 1.5, 0.0, 1.0, 2000.0).unwrap();
        let phases: Vec<_> = [380.0, 780.0].into_iter().map(|wavelength| {
            let solution = sphere.at_wavelength(cx, wavelength, MieBudget::default()).unwrap();
            (wavelength, solution.polarized_phase_table(cx, 1e-3, 2048).unwrap().phase_function())
        }).collect();
        let mie = PhaseFunction::spectral_table(&phases).unwrap();
        assert!(mie.eval_cos(0.0).is_err(), "gray evaluation must not substitute for spectral Mie");
        let map = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap();
        let geometry = InstanceScene::try_new(vec![patch(9, 20.0, 1.0, 0.3)]).unwrap();
        let air = material(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let native = CurvedPolarizedScene::try_new(cx, &geometry, &media,
            vec![SurfaceBinding { object_id: 9, interaction: Interaction::Opaque }], 100_000).unwrap();
        let density = [1.0];
        let grid = VolumeGrid::new([1, 1, 1], &density, [-1.0, -1.0, 0.0], [2.0, 2.0, 1.0]);
        let mass = MassScatteringSpectrum::constant(1.0).unwrap();
        let flight = |ordinal| match ordinal {
            0 => [1.0 - (-0.25_f64).exp(), 0.5 + f64::EPSILON, 0.0, 0.0],
            1 => [1.0 - 1e-12, 0.5, 0.0, 0.0],
            _ => panic!("one real collision then the retained grid exit"),
        };
        for (index, phase) in [PhaseFunction::Rayleigh, mie].iter().enumerate() {
            let volumes = [Some(CurvedRayleighMedium::try_new_volume_grid(
                &mass, &grid, phase, 1,
            ).unwrap())];
            let transport = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
            for q in [1.0, -1.0] {
                let law = SourcePolarization::try_constant(Y, state(q, 0.0, 0.0)).unwrap();
                let source = EnvironmentRadiance::try_new(&map, 2.0).unwrap()
                    .with_polarization(&law).unwrap();
                let camera = CurvedCameraScene::try_new(cx, &native, &[], Some(source)).unwrap();
                let participating = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
                let owner = CurvedParticipatingEnvironmentCameraScene::try_new(cx, &participating).unwrap();
                let observed = owner.sample_ray_with_tracking_samples(
                    cx, camera_ray(), 550.0, 0, participating_settings(),
                    |_| panic!("no surface continuation"), flight,
                    |_| panic!("no surface source proposal"),
                    |ordinal| { assert_eq!(ordinal, 0); [0.5, 0.0] },
                ).unwrap();
                let baseline = participating.sample_ray_with_tracking_samples(
                    cx, camera_ray(), 550.0, 0, native_settings(),
                    |_| panic!("no surface continuation"), flight,
                ).unwrap();
                let source_sample = source.sample_direction(550.0, [0.5, 0.0]).unwrap().unwrap();
                let incoming = camera.launch(camera_ray(), 550.0, 0).unwrap();
                let factor = phase.polarized_connection_factor(
                    incoming, source_sample.frame, source_sample.stokes,
                ).unwrap();
                let pdf = phase.eval_at(550.0, 0.0).unwrap();
                let direct = factor * source_sample.radiance / source_sample.pdf
                    * power_weight(source_sample.pdf, pdf).unwrap() * (-1.0_f64).exp();
                near(observed.environment_radiance_w_m2_sr_nm(), direct, 4e-6);
                if index == 0 {
                    let expected = if q > 0.0 {
                        source_sample.radiance * (24.0 / 25.0) * (-1.0_f64).exp()
                    } else { 0.0 };
                    near(observed.environment_radiance_w_m2_sr_nm(), expected, 4e-6);
                    if q < 0.0 {
                        assert_eq!(observed.environment_radiance_w_m2_sr_nm(), 0.0);
                    }
                }
                let outgoing = observed.radiance_sample().path().beam().ray().dir;
                let continuation_pdf = phase.eval_at(550.0, Z.dot(outgoing)).unwrap();
                let environment_pdf = map.evaluate(outgoing).unwrap().pdf_solid_angle;
                let weight = power_weight(continuation_pdf, environment_pdf).unwrap();
                near(observed.environment_endpoint_weight(), weight, 2e-12);
                near(observed.radiance_sample().radiance_w_m2_sr_nm(),
                    direct + baseline.radiance_sample().radiance_w_m2_sr_nm() * weight, 4e-6);
                assert_eq!(observed.radiance_sample().path(), baseline.radiance_sample().path());
                assert_eq!(observed.volume_events(), 1);
                assert_eq!(observed.free_flight_draws(), 2);
                assert_eq!(observed.source_cell_visits(), if factor > 0.0 { 1 } else { 0 });
                assert_eq!(observed.environment_sample().environment_samples(), 1);
            }
        }
    });
}

#[test]
fn arbitrary_density_binding_does_not_infer_clear_infinite_tail_from_grid_values() {
    with_cx(|cx| {
        let map = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap();
        let source = EnvironmentRadiance::try_new(&map, 3.0).unwrap();
        let geometry = InstanceScene::try_new(vec![
            patch(9, 0.0, 1.0, 0.3), patch(10, 0.0, -1.0, 0.3),
        ]).unwrap();
        let air = material(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let diffuse = [SpectralReflectance::constant(0.6).unwrap()];
        let native = CurvedPolarizedScene::try_new_with_materials(
            cx, &geometry, &media, &[], &[], &diffuse,
            vec![
                SurfaceBinding { object_id: 9, interaction: Interaction::DepolarizingDiffuse {
                    medium: 0, reflectance: 0,
                }},
                SurfaceBinding { object_id: 10, interaction: Interaction::Opaque },
            ],
            100_000,
        ).unwrap();
        let emitters = [SurfaceEmitter::new(10, RadianceSpectrum::try_constant(5.0).unwrap(),
            EmissionSides::Both)];
        let camera = CurvedCameraScene::try_new(cx, &native, &emitters, Some(source)).unwrap();
        let density = [1.0];
        let grid = VolumeGrid::new([1, 1, 1], &density, [-1.0, -1.0, 0.1], [2.0, 2.0, 0.8]);
        let mass = MassScatteringSpectrum::constant(0.7).unwrap();
        let phase = PhaseFunction::Rayleigh;
        // Erasing the native support behind an arbitrary field is deliberate:
        // its identity and sampled zeros do not certify the infinite tail.
        let volumes = [Some(CurvedRayleighMedium::try_new_heterogeneous(
            &mass, &grid, 1.0, ContentHash([177; 32]), &phase,
        ).unwrap())];
        let transport = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        let participating = CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
        let owner = CurvedParticipatingEnvironmentCameraScene::try_new(cx, &participating).unwrap();
        let flight = |_| [1.0 - 1e-12, 0.5, 0.0, 0.0];
        let observed = owner.sample_ray_with_tracking_samples(
            cx, camera_ray(), 550.0, 0, participating_settings(),
            |_| [0.0; 3], flight, |_| [0.5, 0.75], |_| panic!("no collision"),
        ).unwrap();
        let baseline = participating.sample_ray_with_tracking_samples(
            cx, camera_ray(), 550.0, 0, native_settings(), |_| [0.0; 3], flight,
        ).unwrap();
        assert_eq!(observed.radiance_sample(), baseline.radiance_sample());
        assert_eq!(observed.environment_radiance_w_m2_sr_nm(), 0.0);
        assert_eq!(observed.environment_endpoint_weight(), 1.0);
        assert_eq!(observed.environment_sample().visibility_queries(), 0);
        assert_eq!(observed.source_cell_visits(), 0);
    });
}

#[test]
fn finite_log_survival_preserves_bright_sources_and_only_infinite_tail_is_exact_zero() {
    assert_eq!(positive_product_with_log_survival(&[3.0, 7.0], 2.0, 0.0).unwrap(),
        positive_product(&[3.0, 7.0], 2.0).unwrap());
    let rescued = positive_product_with_log_survival(&[1e300], 1.0, -800.0).unwrap();
    let expected = (1e300_f64.ln() - 800.0).exp();
    assert!((rescued / expected - 1.0).abs() < 2e-12);
    assert!(rescued > 0.0);
    assert_eq!(positive_product_with_log_survival(&[1e300], 1.0, f64::NEG_INFINITY).unwrap(), 0.0);
    assert_eq!(positive_product_with_log_survival(&[0.0], 1.0, -1e6).unwrap(), 0.0);
    assert!(positive_product_with_log_survival(&[1.0], 1.0, -1e6).is_err());
    assert!(positive_product_with_log_survival(&[1.0], 1.0, 0.1).is_err());
    assert!(positive_product_with_log_survival(&[1.0], 1.0, f64::NAN).is_err());
}
