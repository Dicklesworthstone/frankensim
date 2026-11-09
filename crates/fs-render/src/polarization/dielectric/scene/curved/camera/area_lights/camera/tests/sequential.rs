//! A steady coated real lens must retain native finite-area radiometry.

use super::*;
use crate::camera::{
    Aperture, CameraProjection, LensSample, LensSurface, PhysicalCamera, SequentialLens,
};
use crate::conductor::{CauchyCoatingLayer, CoherentPeriodicCoating};
use crate::dielectric::DielectricCoating;

#[test]
fn coated_steady_lens_preserves_area_radiometry_without_group_timing() {
    with_cx(|cx| {
        let geometry =
            InstanceScene::try_new(vec![panel(1, 0.0, 1.0, 1.0), panel(9, 0.0, -1.0, 0.3)])
                .unwrap();
        let vacuum = glass(1.0, 0.2);
        let media = [SceneMedium::Homogeneous(&vacuum)];
        let reflectances = [SpectralReflectance::constant(RHO).unwrap()];
        let scene = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![diffuse(1, 0), opaque(9)],
            100_000,
        )
        .unwrap();
        let emitters = [emitter(9)];
        let camera = CurvedCameraScene::try_new(cx, &scene, &emitters, None).unwrap();
        let light = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let seeds = [seed(Vec3::new(0.0, 0.0, 1.0))];
        let bindings = [AreaLightBinding {
            source: &light,
            receiver_object_id: 1,
            source_side: BoundarySide::Positive,
            initial_medium: 0,
            source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            seeds: &seeds,
            image_options: image_options(1),
            position_tolerance_m: 2e-7,
        }];
        let lights = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
        let film_index = 1.5_f64.sqrt();
        let coating = DielectricCoating::try_new(
            CoherentPeriodicCoating::try_new(
                &[CauchyCoatingLayer::try_new(
                    CauchyIor::try_constant(film_index).unwrap(),
                    0.0,
                    550e-9 / (4.0 * film_index),
                )
                .unwrap()],
                1,
                ConductorSource::try_new(
                    hash_domain(
                        "org.frankensim.test.steady-lens-area",
                        b"authored quarter-wave",
                    ),
                    ConductorDataStatus::Representative,
                )
                .unwrap(),
            )
            .unwrap(),
        )
        .unwrap();
        let lens = SequentialLens::try_new(
            550.0,
            &[
                LensSurface::try_new(0.08, 20.0, 0.01, Some(glass(1.5, 2.0)))
                    .unwrap()
                    .with_coating(coating),
                LensSurface::try_new(0.09, -20.0, 0.012, None)
                    .unwrap()
                    .with_coating(coating),
                LensSurface::try_new(0.1, 0.0, 0.02, None).unwrap(),
            ],
        )
        .unwrap();
        let physical = PhysicalCamera::try_look_at(
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 1.0, 0.0),
            CameraProjection::try_focal_sensor(0.05, 0.002).unwrap(),
            1.0,
            Aperture::try_circular(0.01).unwrap(),
        )
        .unwrap()
        .with_sequential_lens(lens)
        .unwrap();
        let beam = physical
            .generate_sequential_spectral_beams_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                &[550.0],
            )
            .unwrap()[0]
            .unwrap();
        assert_eq!(beam.group_delay_s(), None);
        assert!(
            physical
                .generate_sequential_spectral_timed_beams_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    LensSample::CENTER,
                    &[550.0],
                )
                .is_err()
        );
        // Two matched quarter-wave films cancel axial Fresnel loss. Only
        // the authored centimetre of absorbing bulk remains in this seed.
        let lens_transmission = (-2.0_f64 * 0.01).exp();
        near(beam.power().i, lens_transmission, 8e-13);
        let area_draw = [4.0 / 9.0, 0.5];
        let mut workspace = AreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        let ideal = lights
            .sample_ray_with_samples(
                cx,
                beam.ray(),
                550.0,
                0,
                settings(),
                &mut workspace,
                |_| [0.0; 3],
                |_, _| area_draw,
            )
            .unwrap();
        let mut surface_draws = Vec::new();
        let mut area_draws = Vec::new();
        let actual = lights
            .sample_sequential_ray_with_samples(
                cx,
                &beam,
                0,
                settings(),
                &mut workspace,
                |event| {
                    surface_draws.push(event);
                    [0.0; 3]
                },
                |event, binding| {
                    area_draws.push((event, binding));
                    area_draw
                },
            )
            .unwrap();
        assert_eq!(surface_draws, [0]);
        assert_eq!(area_draws, [(0, 0)]);
        assert_eq!(actual.images(), ideal.images());
        assert_eq!(actual.trace_evaluations(), ideal.trace_evaluations());
        assert!(actual.area_radiance_w_m2_sr_nm() > 0.0);
        assert!(actual.radiance_sample().radiance_w_m2_sr_nm() > actual.area_radiance_w_m2_sr_nm());
        // Independent axial finite-area geometry and power-heuristic MIS.
        let p_light = 4.0 / AREA;
        let p_bsdf = 1.0 / core::f64::consts::PI;
        let mis = p_light * p_light / (p_light * p_light + p_bsdf * p_bsdf);
        let expected =
            lens_transmission * (-0.2_f64 * 2.9).exp() * RHO / core::f64::consts::PI * LE * AREA
                / 4.0
                * mis;
        near(actual.area_radiance_w_m2_sr_nm(), expected, 3e-7);
        near(
            actual.radiance_sample().radiance_w_m2_sr_nm(),
            ideal.radiance_sample().radiance_w_m2_sr_nm() * lens_transmission,
            8e-13,
        );
    });
}
