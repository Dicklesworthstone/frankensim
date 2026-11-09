//! Active coated molecular lenses reach steady CIE/QE without a group model.

use super::*;
use crate::camera::{
    Aperture, CameraProjection, LensSample, LensSurface, PhysicalCamera, SequentialCameraBeam,
    SequentialLens,
};
use crate::conductor::{
    CauchyCoatingLayer, CoherentPeriodicCoating, ConductorDataStatus, ConductorSource,
};
use crate::dielectric::DielectricCoating;
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::curved::camera::analyzer::AnalyzerResponse;
use crate::polarization::dielectric::scene::curved::camera::{
    CameraAnalyzer, SpectralCameraAnalyzer,
};
use crate::spectral::{cie_x, cie_y, cie_z, y_integral};

const SEED: u64 = 79;
const SAMPLE: u32 = 6;
const SCATTERING: f64 = 0.01;

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
                    "org.frankensim.test.steady-molecular-area-pixel",
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
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
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

fn check_area<P>(
    sample: &RayleighAreaLightPixelSample<P>,
    beams: &[Option<SequentialCameraBeam>; 4],
) {
    assert_eq!(sample.volume_events(), &[0; 4]);
    assert!(sample.free_flight_draws().iter().all(|draws| *draws > 0));
    assert_eq!(sample.source_free_flight_draws(), &[0; 4]);
    assert_eq!(sample.source_candidates(), &[0; 4]);
    assert_eq!(sample.area_pixel_sample().area_samples(), &[1; 4]);
    assert_eq!(sample.area_pixel_sample().images(), &[1; 4]);
    for lane in 0..4 {
        let beam = beams[lane].unwrap();
        assert_eq!(beam.group_delay_s(), None);
        let ray = beam.ray();
        let prefix = (1.0 - ray.origin.z) / ray.dir.z;
        let receiver = Point3::new(
            ray.origin.x + prefix * ray.dir.x,
            ray.origin.y + prefix * ray.dir.y,
            1.0,
        );
        let draw = philox4x32_10(
            [0, SAMPLE, 0, 0x6372_6100 + lane as u32 * 64],
            [SEED as u32, 0],
        );
        let root = unit(draw[0]).sqrt();
        let v = unit(draw[1]);
        let point = Point3::new(
            -0.3 + 0.6 * root * (1.0 - v) + 0.3 * root * v,
            -0.3 + 0.9 * root * v,
            -1.0,
        );
        let displacement = receiver.delta_from(point);
        let r2 = displacement.dot(displacement);
        let distance = r2.sqrt();
        let cosine = 2.0 / distance;
        let light_pdf = r2 / (0.27 * cosine);
        let continuation_pdf = cosine / core::f64::consts::PI;
        let mis =
            light_pdf * light_pdf / (light_pdf * light_pdf + continuation_pdf * continuation_pdf);
        // Original selected camera flights are uncollided. Molecular survival
        // belongs only to the independent source segment; lens Stokes intensity
        // already includes the detector row, coated Fresnel response and Beer.
        let expected = beam.power().i
            * (-0.2 * (prefix + distance)).exp()
            * (-SCATTERING * distance).exp()
            * 0.6
            / core::f64::consts::PI
            * 4.0
            * 0.27
            * cosine
            * cosine
            / r2
            * mis;
        near(
            sample.area_pixel_sample().area_radiance_w_m2_sr_nm()[lane],
            expected,
        );
        assert!(expected > 0.0);
    }
}

#[test]
fn steady_coated_molecular_pixels_preserve_area_radiometry_and_all_six_detector_modes() {
    fixtures::with_cx(|cx, gate| {
        let geometry = InstanceScene::try_new(vec![
            panel(1, 1.0, 10.0),
            panel(9, -1.0, 0.3),
            panel(10, -1.5, 50.0),
        ])
        .unwrap();
        let vacuum = DielectricGlass::new(
            CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.2).unwrap(),
            GlassProvenance::Custom,
        );
        let media = [SceneMedium::Homogeneous(&vacuum)];
        let reflectances = [SpectralReflectance::constant(0.6).unwrap()];
        let native = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::DepolarizingDiffuse {
                        medium: 0,
                        reflectance: 0,
                    },
                },
                SurfaceBinding {
                    object_id: 9,
                    interaction: Interaction::Opaque,
                },
                SurfaceBinding {
                    object_id: 10,
                    interaction: Interaction::Opaque,
                },
            ],
            100_000,
        )
        .unwrap();
        let emitters = [SurfaceEmitter::new(
            9,
            RadianceSpectrum::try_constant(4.0).unwrap(),
            EmissionSides::Front,
        )];
        let camera = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let source = CurvedAreaEmitter::from_camera(cx, &camera, 9, 1).unwrap();
        let seeds = [BundleImageSeed {
            initial_direction: Vec3::new(0.0, 0.0, 1.0),
            policy: BundleBranchPolicy::TransmitOrReflectTir,
        }];
        let bindings = [AreaLightBinding {
            source: &source,
            receiver_object_id: 1,
            source_side: BoundarySide::Positive,
            initial_medium: 0,
            source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            seeds: &seeds,
            image_options: BundleImageOptions {
                max_images: 1,
                ..BundleImageOptions::default()
            },
            position_tolerance_m: 2e-7,
        }];
        let area = CurvedAreaLightCameraScene::try_new(cx, &camera, &bindings).unwrap();
        let spectrum = MassScatteringSpectrum::constant(SCATTERING).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [Some(
            CurvedRayleighMedium::try_new(&spectrum, 1.0, &phase).unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        let molecular = CurvedRayleighCameraScene::try_new(cx, &camera, &scattering).unwrap();
        let scene = CurvedRayleighAreaLightCameraScene::try_new(cx, &molecular, &area).unwrap();
        let physical = coated_camera(0.02);
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        let mut workspace = RayleighAreaLightCameraWorkspace::try_new(1, 1, 8).unwrap();
        assert!(
            physical
                .generate_sequential_spectral_timed_beams_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    LensSample::CENTER,
                    &[550.0]
                )
                .is_err()
        );
        let analyzer = CameraAnalyzer::try_new(
            Vec3::new(1.0, 0.0, 0.0),
            Stokes {
                i: 0.55,
                q: 0.2,
                u: -0.1,
                v: 0.05,
            },
        )
        .unwrap();
        let spectral = SpectralCameraAnalyzer::try_new(
            Vec3::new(1.0, 0.0, 0.0),
            vec![
                (
                    LAMBDA_MIN,
                    Stokes {
                        i: 0.45,
                        q: 0.1,
                        u: 0.0,
                        v: 0.05,
                    },
                ),
                (
                    LAMBDA_MAX,
                    Stokes {
                        i: 0.65,
                        q: -0.15,
                        u: 0.08,
                        v: -0.03,
                    },
                ),
            ],
        )
        .unwrap();
        let sensor = QuantumEfficiencySensor::try_new(vec![
            (LAMBDA_MIN, [0.25, 0.75, 0.0]),
            (LAMBDA_MAX, [0.25, 0.75, 0.0]),
        ])
        .unwrap();
        let responses = [
            None,
            Some(AnalyzerResponse::Constant(analyzer)),
            Some(AnalyzerResponse::Spectral(&spectral)),
        ];
        let mut visible_work = None;
        let mut sensor_work = None;
        for (mode, response) in responses.into_iter().enumerate() {
            let visible = match mode {
                0 => raster.sample_pixel_sequential_rayleigh_area_lights(
                    cx,
                    &scene,
                    0,
                    settings(),
                    &mut workspace,
                    SEED,
                    0,
                    SAMPLE,
                ),
                1 => raster.sample_pixel_sequential_rayleigh_area_lights_analyzed(
                    cx,
                    &scene,
                    0,
                    analyzer,
                    settings(),
                    &mut workspace,
                    SEED,
                    0,
                    SAMPLE,
                ),
                _ => raster.sample_pixel_sequential_rayleigh_area_lights_spectrally_analyzed(
                    cx,
                    &scene,
                    0,
                    &spectral,
                    settings(),
                    &mut workspace,
                    SEED,
                    0,
                    SAMPLE,
                ),
            }
            .unwrap();
            let pixel = visible.pixel_sample();
            let beams = packet(&physical, cx, pixel.wavelengths_nm(), response);
            check_area(&visible, &beams);
            let mut xyz = [0.0; 3];
            for lane in 0..4 {
                let beam = beams[lane].unwrap();
                assert_eq!(pixel.camera_rays()[lane], Some(beam.ray()));
                assert_eq!(pixel.geometry_weights()[lane], beam.geometry_weight());
                let weight = pixel.spectral_samples()[lane]
                    .unwrap()
                    .radiance_w_m2_sr_nm()
                    * beam.geometry_weight()
                    * (LAMBDA_MAX - LAMBDA_MIN)
                    / 4.0
                    / y_integral();
                for (sum, observer) in xyz.iter_mut().zip([
                    cie_x(beam.wavelength_nm()),
                    cie_y(beam.wavelength_nm()),
                    cie_z(beam.wavelength_nm()),
                ]) {
                    *sum += weight * observer;
                }
            }
            for channel in 0..3 {
                near(pixel.xyz()[channel], xyz[channel]);
            }
            let work = (
                *visible.free_flight_draws(),
                *visible.area_pixel_sample().area_samples(),
                *visible.area_pixel_sample().trace_evaluations(),
            );
            if let Some(original) = visible_work {
                assert_eq!(work, original);
            }
            visible_work = Some(work);

            let detected = match mode {
                0 => raster.sample_pixel_sensor_sequential_rayleigh_area_lights(
                    cx,
                    &scene,
                    &sensor,
                    0,
                    settings(),
                    &mut workspace,
                    SEED,
                    0,
                    SAMPLE,
                ),
                1 => raster.sample_pixel_sensor_sequential_rayleigh_area_lights_analyzed(
                    cx,
                    &scene,
                    &sensor,
                    0,
                    analyzer,
                    settings(),
                    &mut workspace,
                    SEED,
                    0,
                    SAMPLE,
                ),
                _ => raster
                    .sample_pixel_sensor_sequential_rayleigh_area_lights_spectrally_analyzed(
                        cx,
                        &scene,
                        &sensor,
                        0,
                        &spectral,
                        settings(),
                        &mut workspace,
                        SEED,
                        0,
                        SAMPLE,
                    ),
            }
            .unwrap();
            let pixel = detected.pixel_sample();
            let wavelengths = pixel
                .detector_wavelengths()
                .map(|draw| draw.wavelength_nm());
            let beams = packet(&physical, cx, &wavelengths, response);
            check_area(&detected, &beams);
            let mut energy = [0.0; 3];
            let mut electrons = [0.0; 3];
            for lane in 0..4 {
                let beam = beams[lane].unwrap();
                assert_eq!(pixel.camera_rays()[lane], Some(beam.ray()));
                assert_eq!(pixel.geometry_weights()[lane], beam.geometry_weight());
                for channel in 0..3 {
                    let value = pixel.spectral_samples()[lane]
                        .unwrap()
                        .radiance_w_m2_sr_nm()
                        * beam.geometry_weight()
                        * [0.25, 0.75, 0.0][channel]
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
            let work = (
                *detected.free_flight_draws(),
                *detected.area_pixel_sample().area_samples(),
                *detected.area_pixel_sample().trace_evaluations(),
            );
            if let Some(original) = sensor_work {
                assert_eq!(work, original);
            }
            sensor_work = Some(work);
        }

        let closed_camera = coated_camera(1e-12);
        let closed = CurvedCameraRaster::try_new(&closed_camera, 1, 1).unwrap();
        let dark = closed
            .sample_pixel_sensor_sequential_rayleigh_area_lights(
                cx,
                &scene,
                &sensor,
                0,
                settings(),
                &mut workspace,
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
        assert_eq!(dark.area_pixel_sample().area_samples(), &[0; 4]);
        assert_eq!(dark.area_pixel_sample().trace_evaluations(), &[0; 4]);
        assert_eq!(dark.volume_events(), &[0; 4]);
        assert_eq!(dark.free_flight_draws(), &[0; 4]);
        assert_eq!(dark.source_free_flight_draws(), &[0; 4]);
        assert_eq!(dark.source_candidates(), &[0; 4]);
        let incomplete = SpectralCameraAnalyzer::try_new(
            Vec3::new(1.0, 0.0, 0.0),
            vec![
                (500.0, Stokes::unpolarized(0.5)),
                (600.0, Stokes::unpolarized(0.5)),
            ],
        )
        .unwrap();
        assert!(
            closed
                .sample_pixel_sequential_rayleigh_area_lights_spectrally_analyzed(
                    cx,
                    &scene,
                    0,
                    &incomplete,
                    settings(),
                    &mut workspace,
                    SEED,
                    0,
                    SAMPLE
                )
                .is_err()
        );
        assert!(
            closed
                .sample_pixel_sensor_sequential_rayleigh_area_lights_spectrally_analyzed(
                    cx,
                    &scene,
                    &sensor,
                    0,
                    &incomplete,
                    settings(),
                    &mut workspace,
                    SEED,
                    0,
                    SAMPLE
                )
                .is_err()
        );
        gate.request();
        assert!(
            closed
                .sample_pixel_sensor_sequential_rayleigh_area_lights(
                    cx,
                    &scene,
                    &sensor,
                    0,
                    settings(),
                    &mut workspace,
                    SEED,
                    0,
                    SAMPLE
                )
                .is_err()
        );
    });
}
