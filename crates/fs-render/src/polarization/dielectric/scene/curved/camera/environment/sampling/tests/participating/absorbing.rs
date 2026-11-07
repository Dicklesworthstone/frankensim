//! Full-Stokes environment reception through finite scattering/absorbing cells.

use super::*;
use crate::volumes::transport::{AbsorbingGrid, AbsorbingGridBudget, MassAbsorptionSpectrum};

#[test]
fn finite_absorbing_environment_nee_and_camera_mis_share_true_extinction_once() {
    with_cx(|cx| {
        let map = EnvironmentMap::try_from_linear_srgb(1, 1, vec![[1.0; 3]], 0.0).unwrap();
        let geometry = InstanceScene::try_new(vec![patch(9, 20.0, 1.0, 0.3)]).unwrap();
        let air = material(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let native = CurvedPolarizedScene::try_new(
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
        let density = [1.0];
        let grid = VolumeGrid::new([1, 1, 1], &density, [-1.0, -1.0, 0.0], [2.0, 2.0, 1.0]);
        let absorber = AbsorbingGrid::try_new(
            cx,
            VolumeGrid::new([1, 1, 1], &density, [-1.0, -1.0, 0.0], [2.0, 2.0, 1.0]),
            MassAbsorptionSpectrum::constant(0.4).unwrap(),
            None,
            AbsorbingGridBudget {
                max_grid_cells: 1,
                max_steps_per_segment: 1,
            },
        )
        .unwrap();
        let mass = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let laws = [Some(
            CurvedRayleighMedium::try_new_volume_grid(&mass, &grid, &phase, 1).unwrap(),
        )];
        let grids = [Some(&absorber)];
        let transport = CurvedRayleighScene::try_new(cx, &native, &laws)
            .unwrap()
            .with_absorbing_grids(cx, &grids, 4)
            .unwrap();
        let flight = |ordinal| match ordinal {
            0 => [1.0 - (-0.25_f64).exp(), 0.5 + f64::EPSILON, 0.0, 0.0],
            1 => [1.0 - 1e-12, 0.5, 0.0, 0.0],
            _ => panic!("one physical collision then native finite-grid escape"),
        };
        for q in [1.0, -1.0] {
            let polarization = SourcePolarization::try_constant(Y, state(q, 0.0, 0.0)).unwrap();
            let source = EnvironmentRadiance::try_new(&map, 2.0)
                .unwrap()
                .with_polarization(&polarization)
                .unwrap();
            let camera = CurvedCameraScene::try_new(cx, &native, &[], Some(source)).unwrap();
            let participating =
                CurvedRayleighCameraScene::try_new(cx, &camera, &transport).unwrap();
            let owner =
                CurvedParticipatingEnvironmentCameraScene::try_new(cx, &participating).unwrap();
            let sample = owner
                .sample_ray_with_tracking_samples(
                    cx,
                    camera_ray(),
                    550.0,
                    0,
                    ParticipatingEnvironmentCameraSettings {
                        max_source_cells: 2,
                        ..participating_settings()
                    },
                    |_| panic!("no surface"),
                    flight,
                    |_| panic!("no surface proposal"),
                    |_| [0.5, 0.0],
                )
                .unwrap();
            let baseline = participating
                .sample_ray_with_tracking_samples(
                    cx,
                    camera_ray(),
                    550.0,
                    0,
                    native_settings(),
                    |_| panic!("no surface"),
                    flight,
                )
                .unwrap();
            let source_sample = source.sample_direction(550.0, [0.5, 0.0]).unwrap().unwrap();
            let expected_direct = if q > 0.0 {
                source_sample.radiance * (24.0 / 25.0) * (-0.1_f64 - 1.0 - 0.4).exp()
            } else {
                0.0
            };
            near(
                sample.environment_radiance_w_m2_sr_nm(),
                expected_direct,
                4e-6,
            );
            if q < 0.0 {
                assert_eq!(sample.environment_radiance_w_m2_sr_nm(), 0.0);
            }
            let outgoing = sample.radiance_sample().path().beam().ray().dir;
            let phase_pdf = phase.eval_at(550.0, Z.dot(outgoing)).unwrap();
            let env_pdf = map.evaluate(outgoing).unwrap().pdf_solid_angle;
            let weight = power_weight(phase_pdf, env_pdf).unwrap();
            near(sample.environment_endpoint_weight(), weight, 2e-12);
            near(
                sample.radiance_sample().radiance_w_m2_sr_nm(),
                expected_direct + baseline.radiance_sample().radiance_w_m2_sr_nm() * weight,
                4e-6,
            );
            assert_eq!(
                sample.radiance_sample().path(),
                baseline.radiance_sample().path()
            );
            assert_eq!(sample.radiance_sample().path().absorption_cell_visits(), 2);
            assert_eq!(sample.volume_events(), 1);
            assert_eq!(sample.free_flight_draws(), 2);
            assert_eq!(sample.source_cell_visits(), if q > 0.0 { 2 } else { 0 });
            if q > 0.0 {
                assert!(
                    owner
                        .sample_ray_with_tracking_samples(
                            cx,
                            camera_ray(),
                            550.0,
                            0,
                            ParticipatingEnvironmentCameraSettings {
                                max_source_cells: 1,
                                ..participating_settings()
                            },
                            |_| panic!("no surface"),
                            flight,
                            |_| panic!("no surface proposal"),
                            |_| [0.5, 0.0],
                        )
                        .is_err()
                );
            }
        }
    });
}
