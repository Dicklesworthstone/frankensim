//! Chromatic detector rows precede the native coated lens and rough glass.

use super::*;
use crate::camera::{LensSample, LensSurface, SequentialLens};
use crate::conductor::{
    CauchyCoatingLayer, CoherentPeriodicCoating, ConductorDataStatus, ConductorSource,
};
use crate::dielectric::DielectricCoating;

const SEED: u64 = 79;
const SAMPLE: u32 = 6;

fn coated_camera(stop_radius: f64) -> PhysicalCamera {
    let film = DielectricCoating::try_new(
        CoherentPeriodicCoating::try_new(
            &[
                CauchyCoatingLayer::try_new(CauchyIor::try_constant(1.3).unwrap(), 0.02, 110e-9)
                    .unwrap(),
            ],
            1,
            ConductorSource::try_new(
                hash_domain(
                    "org.frankensim.test.rough-glass-detector",
                    b"authored absorbing coating",
                ),
                ConductorDataStatus::Representative,
            )
            .unwrap(),
        )
        .unwrap(),
    )
    .unwrap();
    let glass = DielectricGlass::new(
        CauchyIor::try_new(1.48, 0.006, 0.0).unwrap(),
        BeerLambertAbsorption::try_constant(2.0).unwrap(),
        GlassProvenance::Custom,
    );
    let lens = SequentialLens::try_new(
        550.0,
        &[
            LensSurface::try_new(0.08, 20.0, 0.01, Some(glass))
                .unwrap()
                .with_coating(film),
            LensSurface::try_new(0.09, -20.0, 0.012, None)
                .unwrap()
                .with_coating(film),
            LensSurface::try_new(0.12, 0.0, stop_radius, None).unwrap(),
        ],
    )
    .unwrap();
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.5),
        Point3::new(0.0, 0.0, 0.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(0.05, 0.002).unwrap(),
        1.0,
        Aperture::try_circular(0.01).unwrap(),
    )
    .unwrap()
    .with_sequential_lens(lens)
    .unwrap()
}

fn packet(
    camera: &PhysicalCamera,
    cx: &Cx<'_>,
    wavelengths: &[f64; 4],
    response: Option<AnalyzerResponse<'_>>,
) -> [Option<SequentialCameraBeam>; 4] {
    let key = [SEED as u32, 0];
    let raster = philox4x32_10([0, SAMPLE, 0, 0x6372_7331], key);
    let pupil = philox4x32_10([0, SAMPLE, 0, 0x6372_7332], key);
    let half = camera.raster_vertical_half_tan();
    let x = (2.0 * unit(raster[0]) - 1.0) * half;
    let y = (1.0 - 2.0 * unit(raster[1])) * half;
    let lens = LensSample::try_new(unit(pupil[0]), unit(pupil[1])).unwrap();
    if let Some(response) = response {
        let rows = wavelengths.map(|wavelength| response.at(wavelength).unwrap().row());
        let axis = response.at(wavelengths[0]).unwrap().reference_axis();
        camera
            .generate_sequential_spectral_beams_with_detector_rows_from_tangent_offsets(
                cx,
                x,
                y,
                lens,
                wavelengths,
                axis,
                &rows,
            )
            .unwrap()
    } else {
        camera
            .generate_sequential_spectral_beams_from_tangent_offsets(cx, x, y, lens, wavelengths)
            .unwrap()
    }
}

#[test]
fn g1_full_glass_qe_lens_retains_active_coatings_polarization_and_clipped_proposal_mass() {
    with_cx(|cx, gate| {
        let geometry = geometry();
        let materials = [
            material(1.0, 0.0),
            DielectricGlass::new(
                CauchyIor::try_new(1.48, 0.006, 0.0).unwrap(),
                BeerLambertAbsorption::try_constant(0.24).unwrap(),
                GlassProvenance::Custom,
            ),
        ];
        let media = [
            SceneMedium::Homogeneous(&materials[0]),
            SceneMedium::Homogeneous(&materials[1]),
        ];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let camera = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
        let polarized = [lamp(1.0, 0.17), lamp(-1.0, 0.83)];
        let branches = [
            branch(polarized[0].intensity_source(), 0),
            branch(polarized[1].intensity_source(), 1),
        ];
        let polarization = polarized.each_ref().map(|source| {
            PointLightSourcePolarization::try_new(source.intensity_source(), source).unwrap()
        });
        let points = CurvedPointLightCameraScene::try_new(cx, &camera, &branches)
            .unwrap()
            .with_source_polarizations(cx, &polarization)
            .unwrap();
        let smith = SmithDielectricPointLightCameraScene::new(
            &points,
            SmithWalkSettings::try_new(256).unwrap(),
        );
        let analyzer = analyzer();
        let constant = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.3).unwrap();
        let sensor = sensor();
        let physical = coated_camera(0.02);
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        let render_cx = cx.with_stream_seed(SEED);
        for selected in [
            None,
            Some(AnalyzerResponse::Constant(constant)),
            Some(AnalyzerResponse::Spectral(&analyzer)),
        ] {
            let result = match selected {
                Some(AnalyzerResponse::Spectral(spectral)) =>
                    raster.sample_pixel_sensor_sequential_smith_dielectric_point_lights_spectrally_analyzed(
                        cx, &smith, &sensor, 0, spectral, settings(), SEED, 0, SAMPLE),
                response => raster.sample_pixel_sensor_sequential_smith_dielectric_point_lights(
                    cx, &smith, &sensor, 0, response.map(|r| r.at(550.0).unwrap()), settings(), SEED, 0, SAMPLE),
            }.unwrap();
            let pixel = result.pixel_sample();
            let wavelengths = pixel
                .detector_wavelengths()
                .map(|draw| draw.wavelength_nm());
            let beams = packet(&physical, cx, &wavelengths, selected);
            let mut energy = [0.0; 3];
            let mut electrons = [0.0; 3];
            for lane in 0..4 {
                let beam = beams[lane].unwrap();
                assert_eq!(
                    beam.group_delay_s(),
                    None,
                    "active steady optics need no invented clock"
                );
                assert_eq!(pixel.camera_rays()[lane], Some(beam.ray()));
                assert_eq!(pixel.geometry_weights()[lane], beam.geometry_weight());
                let manual = smith
                    .sample_sequential_ray(
                        &render_cx,
                        &beam,
                        0,
                        settings(),
                        SEED ^ (0x6c67_7300_u64 << 32) ^ ((lane as u64) << 30),
                        u64::from(SAMPLE),
                    )
                    .unwrap();
                assert_eq!(
                    pixel.spectral_samples()[lane],
                    Some(*manual.radiance_sample())
                );
                check_work(&result, lane, &manual);
                assert_eq!(result.connections()[lane], 2);
                let radiance = manual.radiance_sample().radiance_w_m2_sr_nm();
                assert!(radiance > 0.0);
                for channel in 0..3 {
                    let value = radiance
                        * beam.geometry_weight()
                        * qe(wavelengths[lane])[channel]
                        * (LAMBDA_MAX - LAMBDA_MIN)
                        / 4.0;
                    energy[channel] += value;
                    electrons[channel] +=
                        value * wavelengths[lane] * 1e-9 / (6.62607015e-34 * 299_792_458.0);
                }
            }
            for channel in 0..3 {
                near(
                    pixel.qe_weighted_radiance_w_m2_sr()[channel],
                    energy[channel],
                );
                near(
                    pixel.photoelectron_rate_per_m2_sr()[channel],
                    electrons[channel],
                );
            }
        }

        let visible = raster
            .sample_pixel_sequential_smith_dielectric_point_lights_spectrally_analyzed(
                cx,
                &smith,
                0,
                &analyzer,
                settings(),
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        let beams = packet(
            &physical,
            cx,
            visible.pixel_sample().wavelengths_nm(),
            Some(AnalyzerResponse::Spectral(&analyzer)),
        );
        let mut xyz = [0.0; 3];
        for lane in 0..4 {
            let beam = beams[lane].unwrap();
            let manual = smith
                .sample_sequential_ray(
                    &render_cx,
                    &beam,
                    0,
                    settings(),
                    SEED ^ (0x6c67_7300_u64 << 32) ^ ((lane as u64) << 30),
                    u64::from(SAMPLE),
                )
                .unwrap();
            assert_eq!(
                visible.pixel_sample().spectral_samples()[lane],
                Some(*manual.radiance_sample())
            );
            check_work(&visible, lane, &manual);
            for (sum, response) in xyz.iter_mut().zip([
                cie_x(beam.wavelength_nm()),
                cie_y(beam.wavelength_nm()),
                cie_z(beam.wavelength_nm()),
            ]) {
                *sum += manual.radiance_sample().radiance_w_m2_sr_nm()
                    * beam.geometry_weight()
                    * response
                    * (LAMBDA_MAX - LAMBDA_MIN)
                    / 4.0
                    / y_integral();
            }
        }
        for channel in 0..3 {
            near(visible.pixel_sample().xyz()[channel], xyz[channel]);
        }

        let full = raster
            .sample_pixel_sensor_sequential_smith_dielectric_point_lights_spectrally_analyzed(
                cx,
                &smith,
                &sensor,
                0,
                &analyzer,
                settings(),
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        let wavelengths = full
            .pixel_sample()
            .detector_wavelengths()
            .map(|draw| draw.wavelength_nm());
        let beams = packet(
            &physical,
            cx,
            &wavelengths,
            Some(AnalyzerResponse::Spectral(&analyzer)),
        );
        let radii = beams.map(|beam| {
            let ray = beam.unwrap().ray();
            let distance = (0.38 - ray.origin.z) / ray.dir.z;
            (ray.origin.x + distance * ray.dir.x).hypot(ray.origin.y + distance * ray.dir.y)
        });
        let min_radius = radii.iter().copied().fold(f64::INFINITY, f64::min);
        let max_radius = radii.iter().copied().fold(0.0, f64::max);
        assert!(max_radius - min_radius > 1e-10);
        let stop = 0.5 * (min_radius + max_radius);
        let clipped_camera = coated_camera(stop);
        let clipped = CurvedCameraRaster::try_new(&clipped_camera, 1, 1)
            .unwrap()
            .sample_pixel_sensor_sequential_smith_dielectric_point_lights_spectrally_analyzed(
                cx,
                &smith,
                &sensor,
                0,
                &analyzer,
                settings(),
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        let mut live = 0;
        let mut energy = [0.0; 3];
        let mut electrons = [0.0; 3];
        for lane in 0..4 {
            let present = radii[lane] < stop;
            assert_eq!(
                clipped.pixel_sample().spectral_samples()[lane].is_some(),
                present
            );
            if !present {
                assert_eq!(clipped.connections()[lane], 0);
                assert_eq!(clipped.trace_evaluations()[lane], 0);
                assert_eq!(clipped.field_evaluations()[lane], 0);
                assert_eq!(clipped.query_work()[lane], 0);
                continue;
            }
            live += 1;
            assert_eq!(
                clipped.pixel_sample().spectral_samples()[lane],
                full.pixel_sample().spectral_samples()[lane]
            );
            assert_eq!(
                clipped.trace_evaluations()[lane],
                full.trace_evaluations()[lane]
            );
            for channel in 0..3 {
                let value = full.pixel_sample().spectral_samples()[lane]
                    .unwrap()
                    .radiance_w_m2_sr_nm()
                    * beams[lane].unwrap().geometry_weight()
                    * qe(wavelengths[lane])[channel]
                    * (LAMBDA_MAX - LAMBDA_MIN)
                    / 4.0;
                energy[channel] += value;
                electrons[channel] +=
                    value * wavelengths[lane] * 1e-9 / (6.62607015e-34 * 299_792_458.0);
            }
        }
        assert!(live > 0 && live < 4);
        for channel in 0..3 {
            near(
                clipped.pixel_sample().qe_weighted_radiance_w_m2_sr()[channel],
                energy[channel],
            );
            near(
                clipped.pixel_sample().photoelectron_rate_per_m2_sr()[channel],
                electrons[channel],
            );
        }
        let closed_camera = coated_camera(1e-12);
        let closed = CurvedCameraRaster::try_new(&closed_camera, 1, 1).unwrap();
        let dark = closed
            .sample_pixel_sensor_sequential_smith_dielectric_point_lights_spectrally_analyzed(
                cx,
                &smith,
                &sensor,
                0,
                &analyzer,
                settings(),
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        assert_eq!(dark.pixel_sample().camera_rays(), &[None; 4]);
        assert_eq!(
            dark.pixel_sample().qe_weighted_radiance_w_m2_sr(),
            &[0.0; 3]
        );
        assert_eq!(dark.connections(), &[0; 4]);
        assert_eq!(dark.trace_evaluations(), &[0; 4]);
        let incomplete = SpectralCameraAnalyzer::try_new(
            Vec3::new(1.0, 0.0, 0.0),
            vec![
                (500.0, Stokes::unpolarized(0.5)),
                (600.0, Stokes::unpolarized(0.5)),
            ],
        )
        .unwrap();
        assert!(matches!(
            closed
                .sample_pixel_sensor_sequential_smith_dielectric_point_lights_spectrally_analyzed(
                    cx,
                    &smith,
                    &sensor,
                    0,
                    &incomplete,
                    settings(),
                    SEED,
                    0,
                    SAMPLE
                ),
            Err(PointLightRasterError::Raster(CurvedRasterError::Transport(
                _
            )))
        ));
        assert!(matches!(
            raster.sample_pixel_sensor_smith_dielectric_point_lights_spectrally_analyzed(
                cx,
                &smith,
                &sensor,
                0,
                &incomplete,
                settings(),
                SEED,
                0,
                SAMPLE
            ),
            Err(PointLightRasterError::Raster(CurvedRasterError::Transport(
                _
            )))
        ));
        gate.request();
        assert!(
            closed
                .sample_pixel_sensor_sequential_smith_dielectric_point_lights(
                    cx,
                    &smith,
                    &sensor,
                    0,
                    None,
                    settings(),
                    SEED,
                    0,
                    SAMPLE
                )
                .is_err()
        );
    });
}
