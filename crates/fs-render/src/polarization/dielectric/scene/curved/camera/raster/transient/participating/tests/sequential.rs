use super::*;
use crate::camera::{
    Aperture, CameraProjection, LensSample, LensSurface, PhysicalCamera, SequentialCameraBeam,
    SequentialLens,
};
use crate::dielectric::{DielectricIorSample, TabulatedIor};
use crate::polarization::dielectric::scene::curved::camera::CameraAnalyzer;
use fs_geom::{Point3, Vec3};
use std::sync::Mutex;

const SEED: u64 = 79;
const SAMPLE: u32 = 6;
const C0: f64 = 299_792_458.0;

fn visible_glass() -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_new(1.48, 0.006, 0.0).unwrap(),
        BeerLambertAbsorption::try_constant(2.0).unwrap(),
        GlassProvenance::Custom,
    )
}

fn lens_camera(material: DielectricGlass, reference: f64, stop: Option<f64>) -> PhysicalCamera {
    let mut caps = vec![
        LensSurface::try_new(0.08, 20.0, 0.01, Some(material)).unwrap(),
        LensSurface::try_new(0.09, -20.0, 0.012, None).unwrap(),
        LensSurface::try_new(0.1, 0.0, 0.02, None).unwrap(),
    ];
    if let Some(radius) = stop {
        caps.push(LensSurface::try_new(0.12, 0.0, radius, None).unwrap());
    }
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(0.05, 0.002).unwrap(),
        1.0,
        Aperture::try_circular(0.01).unwrap(),
    )
    .unwrap()
    .with_sequential_lens(SequentialLens::try_new(reference, &caps).unwrap())
    .unwrap()
}

fn packet(
    camera: &PhysicalCamera,
    cx: &Cx<'_>,
    wavelengths: &[f64; CAMERA_WAVELENGTHS],
) -> [Option<SequentialCameraBeam>; CAMERA_WAVELENGTHS] {
    let key = [SEED as u32, 0];
    let raster = philox4x32_10([0, SAMPLE, 0, 0x6372_7331], key);
    let pupil = philox4x32_10([0, SAMPLE, 0, 0x6372_7332], key);
    let half = camera.raster_vertical_half_tan();
    camera
        .generate_sequential_spectral_timed_beams_from_tangent_offsets(
            cx,
            (2.0 * unit(raster[0]) - 1.0) * half,
            (1.0 - 2.0 * unit(raster[1])) * half,
            LensSample::try_new(unit(pupil[0]), unit(pupil[1])).unwrap(),
            wavelengths,
        )
        .unwrap()
}

#[test]
fn analyzed_scattering_pixels_keep_clipped_measure_sensor_clock_and_work_identity() {
    with_cx(|cx, gate| {
        let geometry = geometry();
        let vacuum = DielectricGlass::new(
            CauchyIor::try_constant(1.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.12).unwrap(),
            GlassProvenance::Custom,
        );
        let media = [SceneMedium::Homogeneous(&vacuum)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let sources = sources();
        let surface = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let wide = lens_camera(visible_glass(), 550.0, None);
        let full = CurvedCameraRaster::try_new(&wide, 1, 1)
            .unwrap()
            .sample_pixel_sequential(cx, &surface, 0, settings(), SEED, 0, SAMPLE)
            .unwrap();
        let mut radii = packet(&wide, cx, full.wavelengths_nm()).map(|beam| {
            let ray = beam.unwrap().ray();
            let p = ray
                .origin
                .offset(ray.dir.scale((0.12 - ray.origin.z) / ray.dir.z));
            p.x.hypot(p.y)
        });
        radii.sort_by(f64::total_cmp);
        assert!(radii[2] - radii[1] > 1e-8);
        let physical = lens_camera(visible_glass(), 550.0, Some((radii[1] + radii[2]) * 0.5));
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        let analyzer = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.3).unwrap();
        let queries = Mutex::new(Vec::new());
        let history = History::new(|query, _| {
            queries.lock().unwrap().push(query);
            Ok(1.0)
        });
        let transient = TransientCurvedCameraScene::try_new(cx, &surface, &[history.binding()])
            .unwrap()
            .with_analyzer(analyzer);
        let clock = Dd::from_f64(2.0_f64.powi(30)) + Dd::from_f64(1e-8);
        let phase = PhaseFunction::Rayleigh;
        for coefficient in [0.0, 0.5] {
            let sigma = MassScatteringSpectrum::constant(coefficient).unwrap();
            let volumes = [Some(
                CurvedRayleighMedium::try_new(&sigma, 1.0, &phase).unwrap(),
            )];
            let scattering = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
            let steady = CurvedRayleighCameraScene::try_new(cx, &surface, &scattering).unwrap();
            let timed =
                TransientCurvedRayleighCameraScene::try_new(cx, &transient, &scattering).unwrap();
            let expected = raster
                .sample_pixel_sequential_rayleigh_analyzed(
                    cx,
                    &steady,
                    0,
                    analyzer,
                    trace_settings(),
                    SEED,
                    0,
                    SAMPLE,
                )
                .unwrap();
            let before = queries.lock().unwrap().len();
            let actual = raster
                .sample_pixel_sequential_rayleigh_transient_at_time(
                    cx,
                    &timed,
                    0,
                    clock,
                    trace_settings(),
                    SEED,
                    0,
                    SAMPLE,
                )
                .unwrap();
            assert_eq!(actual.pixel_sample(), expected.pixel_sample());
            assert_eq!(actual.volume_events(), expected.volume_events());
            assert_eq!(actual.free_flight_draws(), expected.free_flight_draws());
            assert_eq!(actual.observation_time_s(), clock);
            assert_eq!(
                actual.pixel_sample().camera_rays().iter().flatten().count(),
                2
            );
            let beams = packet(&physical, cx, actual.pixel_sample().wavelengths_nm());
            let seen = queries.lock().unwrap()[before..].to_vec();
            assert_eq!(seen.len(), 2);
            let mut source = 0;
            for lane in 0..CAMERA_WAVELENGTHS {
                if let Some(beam) = beams[lane] {
                    let path = actual.pixel_sample().spectral_samples()[lane].unwrap();
                    let delay = beam.group_delay_s().unwrap()
                        + Dd::from_f64(path.path().group_optical_path_length_m().unwrap())
                            / Dd::from_f64(C0);
                    assert_eq!(seen[source].observation_time_s, clock);
                    assert_eq!(seen[source].group_delay_s, delay);
                    assert_eq!(seen[source].emission_time_s, clock - delay);
                    source += 1;
                    // At 0.5/m the four first flights are 0.89..1.52 m,
                    // all inside the original 2 m box. No seeded luck is needed.
                    if coefficient > 0.0 {
                        assert!(actual.volume_events()[lane] >= 1);
                    } else {
                        assert_eq!(actual.volume_events()[lane], 0);
                    }
                } else {
                    assert_eq!(actual.volume_events()[lane], 0);
                    assert_eq!(actual.free_flight_draws()[lane], 0);
                }
            }
            assert_eq!(
                actual,
                raster
                    .sample_pixel_sequential_rayleigh_transient_at_time(
                        cx,
                        &timed,
                        0,
                        clock,
                        trace_settings(),
                        SEED,
                        0,
                        SAMPLE,
                    )
                    .unwrap()
            );
        }
        let volumes = [None];
        let scattering = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        let timed =
            TransientCurvedRayleighCameraScene::try_new(cx, &transient, &scattering).unwrap();
        let blocked_camera = lens_camera(visible_glass(), 550.0, Some(radii[0] * 0.5));
        let blocked = CurvedCameraRaster::try_new(&blocked_camera, 1, 1).unwrap();
        let before = history.calls.load(Ordering::Relaxed);
        let black = blocked
            .sample_pixel_sequential_rayleigh_transient_at_time(
                cx,
                &timed,
                0,
                clock,
                trace_settings(),
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        assert!(
            black
                .pixel_sample()
                .camera_rays()
                .iter()
                .all(Option::is_none)
        );
        assert_eq!(black.volume_events(), &[0; CAMERA_WAVELENGTHS]);
        assert_eq!(black.free_flight_draws(), &[0; CAMERA_WAVELENGTHS]);
        assert_eq!(history.calls.load(Ordering::Relaxed), before);
        assert_eq!(
            blocked
                .sample_pixel_sequential_rayleigh_transient_at_time(
                    cx,
                    &timed,
                    0,
                    Dd {
                        hi: f64::NAN,
                        lo: 0.0
                    },
                    trace_settings(),
                    SEED,
                    0,
                    SAMPLE,
                )
                .unwrap_err(),
            TransientRasterError::Source(TransientCameraError::InvalidObservationTime)
        );
        gate.request();
        assert!(
            blocked
                .sample_pixel_sequential_rayleigh_transient_at_time(
                    cx,
                    &timed,
                    0,
                    clock,
                    trace_settings(),
                    SEED,
                    0,
                    SAMPLE,
                )
                .is_err()
        );
    });
}

#[test]
fn infrared_molecular_histories_share_qe_paths_and_refuse_partial_source_failures() {
    with_cx(|cx, _| {
        let ior = |base: f64, slope: f64| {
            TabulatedIor::try_new_in_band(
                core::array::from_fn(|i| {
                    DielectricIorSample::try_new_in_band(
                        8000.0 + 750.0 * i as f64,
                        base + slope * i as f64,
                    )
                    .unwrap()
                }),
                750.0,
            )
            .unwrap()
        };
        let material = DielectricGlass::new(
            ior(1.5625, -1.0 / 128.0),
            BeerLambertAbsorption::try_constant(2.0).unwrap(),
            GlassProvenance::Custom,
        );
        let vacuum = DielectricGlass::new(
            ior(1.0, 0.0),
            BeerLambertAbsorption::try_constant(0.12).unwrap(),
            GlassProvenance::Custom,
        );
        let geometry = geometry();
        let media = [SceneMedium::Homogeneous(&vacuum)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let sources = [SurfaceEmitter::new(
            7,
            RadianceSpectrum::try_new_in_band(vec![(8000.0, 3.0), (14000.0, 3.0)]).unwrap(),
            EmissionSides::Both,
        )];
        let surface = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let sigma = MassScatteringSpectrum::source_table(
            &[(8000.0, 0.5), (14000.0, 0.5)],
            ContentHash([92; 32]),
        )
        .unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [Some(
            CurvedRayleighMedium::try_new(&sigma, 1.0, &phase).unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        let steady = CurvedRayleighCameraScene::try_new(cx, &surface, &scattering).unwrap();
        let physical = lens_camera(material, 11000.0, None);
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let clock = epoch + Dd::from_f64(1e-7);
        let history = History::new(|query: SurfaceEmissionQuery, _| {
            Ok((-((query.emission_time_s - epoch).to_f64() / 1e-7).powi(2)).exp())
        });
        let transient =
            TransientCurvedCameraScene::try_new(cx, &surface, &[history.binding()]).unwrap();
        let timed =
            TransientCurvedRayleighCameraScene::try_new(cx, &transient, &scattering).unwrap();
        let sensor = QuantumEfficiencySensor::try_new(vec![
            (8000.0, [0.0, 1.0, 0.0]),
            (14000.0, [1.0, 0.0, 0.0]),
        ])
        .unwrap();
        let base = raster
            .sample_pixel_sensor_sequential_rayleigh(
                cx,
                &steady,
                &sensor,
                0,
                trace_settings(),
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        let actual = raster
            .sample_pixel_sensor_sequential_rayleigh_transient_at_time(
                cx,
                &timed,
                &sensor,
                0,
                clock,
                trace_settings(),
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        assert_eq!(actual.observation_time_s(), clock);
        assert_eq!(actual.volume_events(), base.volume_events());
        assert_eq!(actual.free_flight_draws(), base.free_flight_draws());
        assert!(actual.volume_events().iter().all(|&n| n > 0));
        assert_eq!(history.calls.load(Ordering::Relaxed), 4);
        let output = actual.pixel_sample();
        let wavelengths = output
            .detector_wavelengths()
            .map(|draw| draw.wavelength_nm());
        let beams = packet(&physical, cx, &wavelengths);
        let mut energy = [0.0; 3];
        let mut electrons = [0.0; 3];
        for lane in 0..CAMERA_WAVELENGTHS {
            let beam = beams[lane].unwrap();
            let reference = base.pixel_sample().spectral_samples()[lane].unwrap();
            let path = reference.path();
            assert_eq!(output.spectral_samples()[lane].unwrap().path(), path);
            let delay = beam.group_delay_s().unwrap()
                + Dd::from_f64(path.group_optical_path_length_m().unwrap()) / Dd::from_f64(C0);
            let scale = (-(((clock - epoch) - delay).to_f64() / 1e-7).powi(2)).exp();
            assert!(scale > 0.0);
            let radiance = reference.radiance_w_m2_sr_nm() * scale;
            near(
                output.spectral_samples()[lane]
                    .unwrap()
                    .radiance_w_m2_sr_nm(),
                radiance,
            );
            let nm = wavelengths[lane];
            let qe = [(nm - 8000.0) / 6000.0, (14000.0 - nm) / 6000.0, 0.0];
            for c in 0..3 {
                let contribution = radiance * beam.geometry_weight() * qe[c] * 6000.0 / 4.0;
                energy[c] += contribution;
                electrons[c] += contribution * nm * 1e-9 / (6.62607015e-34 * C0);
            }
        }
        for c in 0..3 {
            near(output.qe_weighted_radiance_w_m2_sr()[c], energy[c]);
            near(output.photoelectron_rate_per_m2_sr()[c], electrons[c]);
        }
        assert_eq!(output.photoelectron_rate_per_m2_sr()[2], 0.0);
        let failure = History::new(|_, _| Err(SurfaceHistoryError::OutsideTimeDomain));
        let failed =
            TransientCurvedCameraScene::try_new(cx, &surface, &[failure.binding()]).unwrap();
        let failed_scene =
            TransientCurvedRayleighCameraScene::try_new(cx, &failed, &scattering).unwrap();
        assert_eq!(
            raster
                .sample_pixel_sensor_sequential_rayleigh_transient_at_time(
                    cx,
                    &failed_scene,
                    &sensor,
                    0,
                    clock,
                    trace_settings(),
                    SEED,
                    0,
                    SAMPLE,
                )
                .unwrap_err(),
            TransientRasterError::Source(TransientCameraError::History {
                object_id: 7,
                error: SurfaceHistoryError::OutsideTimeDomain,
            })
        );
        assert_eq!(failure.calls.load(Ordering::Relaxed), 1);
        let mut limited = trace_settings();
        limited.max_volume_events = 1;
        assert_eq!(
            raster
                .sample_pixel_sensor_sequential_rayleigh_transient_at_time(
                    cx,
                    &failed_scene,
                    &sensor,
                    0,
                    clock,
                    limited,
                    SEED,
                    0,
                    SAMPLE,
                )
                .unwrap_err(),
            TransientRasterError::Source(TransientCameraError::Camera(
                CameraRadianceError::Unresolved(PathTermination::VolumeLimit),
            ))
        );
        assert_eq!(failure.calls.load(Ordering::Relaxed), 1);
    });
}
