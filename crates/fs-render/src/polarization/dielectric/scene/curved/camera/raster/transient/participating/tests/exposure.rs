use super::*;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::piecewise_linear::PiecewiseLinearHistory;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::{
    ExposureWindow, SurfaceExposureQuery,
};

struct AverageOnly<'a> {
    history: PiecewiseLinearHistory<'a>,
    calls: AtomicUsize,
}
impl SurfaceEmissionHistory for AverageOnly<'_> {
    fn sample(&self, _: &Cx<'_>, _: SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError> {
        panic!("finite shutter must query each source mean exactly once per path")
    }
    fn average(
        &self,
        cx: &Cx<'_>,
        query: SurfaceExposureQuery,
    ) -> Result<f64, SurfaceHistoryError> {
        self.calls.fetch_add(1, Ordering::Relaxed);
        self.history.average(cx, query)
    }
}

#[test]
fn g3_scattered_exposure_keeps_native_rng_and_integrates_all_qe_channels_from_four_paths() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let glass = glass();
        let media = [SceneMedium::Homogeneous(&glass)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let sources = sources();
        let surface = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let sigma = MassScatteringSpectrum::constant(0.5).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [Some(
            CurvedRayleighMedium::try_new(&sigma, 1.0, &phase).unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        let steady = CurvedRayleighCameraScene::try_new(cx, &surface, &scattering).unwrap();
        let constant = [(0.0, 1.0), (1.0, 1.0)];
        let unity = AverageOnly {
            history: PiecewiseLinearHistory::try_new(cx, Dd::ZERO, &constant).unwrap(),
            calls: AtomicUsize::new(0),
        };
        let transient = TransientCurvedCameraScene::try_new(
            cx,
            &surface,
            &[SurfaceEmissionHistoryBinding {
                object_id: 7,
                source_identity: ContentHash([93; 32]),
                history: &unity,
            }],
        )
        .unwrap();
        let timed =
            TransientCurvedRayleighCameraScene::try_new(cx, &transient, &scattering).unwrap();
        let optical = camera(0.01);
        let raster = CurvedCameraRaster::try_new(&optical, 3, 2).unwrap();
        let window = ExposureWindow::try_new(Dd::from_f64(1e-6), Dd::from_f64(2e-6)).unwrap();
        let original = raster
            .sample_pixel_rayleigh(cx, &steady, 0, trace_settings(), 73, 4, 9)
            .unwrap();
        let exposed = raster
            .sample_pixel_rayleigh_during_exposure(
                cx,
                &timed,
                0,
                window,
                trace_settings(),
                73,
                4,
                9,
            )
            .unwrap();
        assert_eq!(exposed.mean_pixel_sample(), original.pixel_sample());
        assert_eq!(exposed.volume_events(), original.volume_events());
        assert_eq!(exposed.free_flight_draws(), original.free_flight_draws());
        assert!(exposed.volume_events()[3] >= 1);
        assert_eq!(unity.calls.load(Ordering::Relaxed), 4);

        // The exact source interpolant L(t)/L0=t/(1 microsecond) has mean
        // 1.5-delay/(1 microsecond) over this shutter, independently per path.
        let linear = [(0.0, 0.0), (1.0, 1e6)];
        let ramp = AverageOnly {
            history: PiecewiseLinearHistory::try_new(cx, Dd::ZERO, &linear).unwrap(),
            calls: AtomicUsize::new(0),
        };
        let transient = TransientCurvedCameraScene::try_new(
            cx,
            &surface,
            &[SurfaceEmissionHistoryBinding {
                object_id: 7,
                source_identity: ContentHash([94; 32]),
                history: &ramp,
            }],
        )
        .unwrap();
        let timed =
            TransientCurvedRayleighCameraScene::try_new(cx, &transient, &scattering).unwrap();
        let one = QuantumEfficiencySensor::try_new(vec![(420.0, [0.5]), (680.0, [0.5])]).unwrap();
        let three = QuantumEfficiencySensor::try_new(vec![
            (420.0, [0.2, 0.3, 0.0]),
            (680.0, [0.2, 0.3, 0.0]),
        ])
        .unwrap();
        let original = raster
            .sample_pixel_sensor_rayleigh(cx, &steady, &one, 0, trace_settings(), 73, 4, 9)
            .unwrap();
        let a = raster
            .sample_pixel_sensor_rayleigh_during_exposure(
                cx,
                &timed,
                &one,
                0,
                window,
                trace_settings(),
                73,
                4,
                9,
            )
            .unwrap();
        let b = raster
            .sample_pixel_sensor_rayleigh_during_exposure(
                cx,
                &timed,
                &three,
                0,
                window,
                trace_settings(),
                73,
                4,
                9,
            )
            .unwrap();
        assert_eq!(ramp.calls.load(Ordering::Relaxed), 8);
        assert_eq!(a.volume_events(), original.volume_events());
        assert_eq!(a.free_flight_draws(), original.free_flight_draws());
        assert_eq!(a.volume_events(), b.volume_events());
        assert_eq!(a.free_flight_draws(), b.free_flight_draws());
        assert_eq!(
            a.mean_pixel_sample().spectral_samples(),
            b.mean_pixel_sample().spectral_samples()
        );
        let mut energy = 0.0;
        let mut electrons = 0.0;
        for lane in 0..CAMERA_WAVELENGTHS {
            let pixel = a.mean_pixel_sample();
            let draw = pixel.detector_wavelengths()[lane];
            let wave_um = draw.wavelength_nm() * 1e-3;
            let steady = original.pixel_sample().spectral_samples()[lane];
            let delay = steady.path().traveled_m() * (1.4 + 0.06 / wave_um.powi(2)) / 299_792_458.0;
            let mean = steady.radiance_w_m2_sr_nm() * (1.5 - delay / 1e-6);
            near(pixel.spectral_samples()[lane].radiance_w_m2_sr_nm(), mean);
            let integral =
                mean * window.duration_s().to_f64() * pixel.diffraction_weights()[lane] / 4.0;
            energy += integral * draw.radiance_weights_nm()[0];
            electrons += integral * draw.electron_weights_nm_per_joule()[0];
        }
        near(a.qe_weighted_energy_j_m2_sr()[0], energy);
        near(a.photoelectrons_per_m2_sr()[0], electrons);
        near(
            a.photoelectrons_per_m2_sr()[0],
            b.photoelectrons_per_m2_sr()[0] + b.photoelectrons_per_m2_sr()[1],
        );
        assert_eq!(b.photoelectrons_per_m2_sr()[2], 0.0);

        let unsupported = History::new(|_, _| Ok(1.0));
        let transient =
            TransientCurvedCameraScene::try_new(cx, &surface, &[unsupported.binding()]).unwrap();
        let timed =
            TransientCurvedRayleighCameraScene::try_new(cx, &transient, &scattering).unwrap();
        assert!(matches!(
            raster.sample_pixel_sensor_rayleigh_during_exposure(
                cx,
                &timed,
                &one,
                0,
                window,
                trace_settings(),
                73,
                4,
                9
            ),
            Err(TransientRasterError::Source(
                TransientCameraError::History {
                    object_id: 7,
                    error: SurfaceHistoryError::ExposureUnsupported
                }
            ))
        ));
        assert_eq!(unsupported.calls.load(Ordering::Relaxed), 0);
    });
}
