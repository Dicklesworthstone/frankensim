use super::*;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::piecewise_linear::PiecewiseLinearHistory;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::{
    ExposureWindow, SurfaceExposureQuery,
};

struct ExposureOnly<'a> {
    history: PiecewiseLinearHistory<'a>,
    queries: Mutex<Vec<SurfaceExposureQuery>>,
}
impl SurfaceEmissionHistory for ExposureOnly<'_> {
    fn sample(&self, _: &Cx<'_>, _: SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError> {
        panic!("one source-owned mean must replace all temporal point samples")
    }
    fn average(
        &self,
        cx: &Cx<'_>,
        query: SurfaceExposureQuery,
    ) -> Result<f64, SurfaceHistoryError> {
        cx.checkpoint()?;
        self.queries.lock().unwrap().push(query);
        self.history.average(cx, query)
    }
}

#[test]
fn g1_narrow_source_pulse_integrates_after_two_molecular_turns_without_extra_paths() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let glass = DielectricGlass::new(
            CauchyIor::try_new(1.4, 0.02, 0.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.2).unwrap(),
            GlassProvenance::Custom,
        );
        let media = [SceneMedium::Homogeneous(&glass)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let sources = sources();
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let sigma = MassScatteringSpectrum::constant(1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [Some(
            CurvedRayleighMedium::try_new(&sigma, 1.0, &phase).unwrap(),
        )];
        let scattering = CurvedRayleighScene::try_new(cx, &native, &volumes).unwrap();
        let steady = CurvedRayleighCameraScene::try_new(cx, &camera, &scattering).unwrap();
        let original = steady
            .sample_ray_with_samples(
                cx,
                ray(),
                550.0,
                0,
                settings(),
                |_| panic!("opaque"),
                two_turns,
            )
            .unwrap();
        let epoch = Dd::from_f64(1_073_741_824.0);
        let delay = 1.5 * (1.4 + 0.06 / 0.55_f64.powi(2)) / LIGHT_M_S;
        let window = ExposureWindow::try_new(
            epoch + Dd::from_f64(delay),
            epoch + Dd::from_f64(delay + 2e-9),
        )
        .unwrap();
        let knots = [
            (-1.0, 0.0),
            (0.999e-9, 0.0),
            (1e-9, 1e3),
            (1.001e-9, 0.0),
            (1.0, 0.0),
        ];
        let pulse = PiecewiseLinearHistory::try_new(cx, epoch, &knots).unwrap();
        let history = ExposureOnly {
            history: pulse,
            queries: Mutex::new(Vec::new()),
        };
        let transient = TransientCurvedCameraScene::try_new(
            cx,
            &camera,
            &[SurfaceEmissionHistoryBinding {
                object_id: 9,
                source_identity: ContentHash([91; 32]),
                history: &history,
            }],
        )
        .unwrap();
        let timed =
            TransientCurvedRayleighCameraScene::try_new(cx, &transient, &scattering).unwrap();
        let mut flights = Vec::new();
        let exposed = timed
            .sample_ray_with_samples_during_exposure(
                cx,
                ray(),
                550.0,
                0,
                window,
                settings(),
                |_| panic!("opaque"),
                |ordinal| {
                    flights.push(ordinal);
                    two_turns(ordinal)
                },
            )
            .unwrap();
        let pulse_area = 0.5 * 1e3 * (knots[3].0 - knots[1].0);
        let steady_radiance = 8.0 * (-0.3_f64).exp();
        near(
            exposed.integrated_radiance_j_m2_sr_nm(),
            steady_radiance * pulse_area,
            3e-15,
        );
        near(
            exposed.mean_radiance_sample().radiance_w_m2_sr_nm(),
            steady_radiance * pulse_area / window.duration_s().to_f64(),
            3e-6,
        );
        assert_eq!(
            exposed.mean_radiance_sample().path(),
            original.radiance_sample().path()
        );
        assert_eq!(exposed.volume_events(), 2);
        assert_eq!(exposed.free_flight_draws(), 3);
        assert_eq!(flights, [0, 1, 2]);
        let queries = history.queries.lock().unwrap();
        assert_eq!(queries.len(), 1);
        near(queries[0].opening.group_delay_s.to_f64(), delay, 3e-16);
        near(
            (queries[0].opening.emission_time_s - epoch).to_f64(),
            0.0,
            3e-16,
        );
        drop(queries);

        // An exact unit mean must retain even the steady estimator's bits,
        // with the independent fourth tracking variate still unmodified.
        let constant = [(-1.0, 1.0), (1.0, 1.0)];
        let unity = ExposureOnly {
            history: PiecewiseLinearHistory::try_new(cx, epoch, &constant).unwrap(),
            queries: Mutex::new(Vec::new()),
        };
        let transient = TransientCurvedCameraScene::try_new(
            cx,
            &camera,
            &[SurfaceEmissionHistoryBinding {
                object_id: 9,
                source_identity: ContentHash([92; 32]),
                history: &unity,
            }],
        )
        .unwrap();
        let timed =
            TransientCurvedRayleighCameraScene::try_new(cx, &transient, &scattering).unwrap();
        let mut flights = Vec::new();
        let exposed = timed
            .sample_ray_with_tracking_samples_during_exposure(
                cx,
                ray(),
                550.0,
                0,
                window,
                settings(),
                |_| panic!("opaque"),
                |ordinal| {
                    flights.push(ordinal);
                    let [a, b, c] = two_turns(ordinal);
                    [a, b, c, 0.5]
                },
            )
            .unwrap();
        assert_eq!(exposed.mean_radiance_sample(), original.radiance_sample());
        assert_eq!(exposed.volume_events(), original.volume_events());
        assert_eq!(exposed.free_flight_draws(), original.free_flight_draws());
        assert_eq!(flights, [0, 1, 2]);
        assert_eq!(unity.queries.lock().unwrap().len(), 1);
    });
}
