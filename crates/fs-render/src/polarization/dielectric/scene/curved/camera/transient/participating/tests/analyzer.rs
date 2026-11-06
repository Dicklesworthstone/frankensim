use super::*;
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::curved::camera::analyzer::CameraAnalyzer;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::piecewise_linear::PiecewiseLinearHistory;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::{
    ExposureWindow, SurfaceExposureQuery,
};

struct RecordedPulse<'a> {
    source: PiecewiseLinearHistory<'a>,
    intervals: Mutex<Vec<SurfaceExposureQuery>>,
}
impl SurfaceEmissionHistory for RecordedPulse<'_> {
    fn sample(&self, cx: &Cx<'_>, query: SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError> {
        self.source.sample(cx, query)
    }
    fn average(
        &self,
        cx: &Cx<'_>,
        query: SurfaceExposureQuery,
    ) -> Result<f64, SurfaceHistoryError> {
        self.intervals.lock().unwrap().push(query);
        self.source.average(cx, query)
    }
}

#[test]
fn g1_analyzed_molecular_pulse_keeps_complete_stokes_and_its_original_interval() {
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
        let scalar = steady
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
        let delay = Dd::from_f64(
            scalar
                .radiance_sample()
                .path()
                .group_optical_path_length_m()
                .unwrap(),
        ) / Dd::from_f64(LIGHT_M_S);
        let window =
            ExposureWindow::try_new(epoch + delay, epoch + delay + Dd::from_f64(2e-9)).unwrap();
        // A two-picosecond pulse wholly inside the shutter, away from its
        // midpoint. Point-sampled exposure would incorrectly report zero.
        let knots = [
            (-1.0, 0.0),
            (1.2e-9, 0.0),
            (1.201e-9, 1e3),
            (1.202e-9, 0.0),
            (1.0, 0.0),
        ];
        let pulse = RecordedPulse {
            source: PiecewiseLinearHistory::try_new(cx, epoch, &knots).unwrap(),
            intervals: Mutex::new(Vec::new()),
        };
        let area = 0.5 * 1e3 * (knots[3].0 - knots[1].0);
        let mut total_energy = 0.0;
        for sign in [-1.0, 1.0] {
            let analyzer = CameraAnalyzer::try_new(
                Vec3::new(1.0, 0.0, 0.0),
                Stokes {
                    i: 0.5,
                    q: sign * 0.2,
                    u: sign * 0.15,
                    v: sign * 0.1,
                },
            )
            .unwrap();
            let expected = steady
                .sample_ray_analyzed_with_samples(
                    cx,
                    ray(),
                    550.0,
                    0,
                    analyzer,
                    settings(),
                    |_| panic!("opaque"),
                    two_turns,
                )
                .unwrap();
            let transient = TransientCurvedCameraScene::try_new(
                cx,
                &camera,
                &[SurfaceEmissionHistoryBinding {
                    object_id: 9,
                    source_identity: ContentHash([91; 32]),
                    history: &pulse,
                }],
            )
            .unwrap()
            .with_analyzer(analyzer);
            let timed =
                TransientCurvedRayleighCameraScene::try_new(cx, &transient, &scattering).unwrap();
            let time = epoch + delay + Dd::from_f64(knots[2].0);
            let instantaneous = timed
                .sample_ray_with_samples_at_time(
                    cx,
                    ray(),
                    550.0,
                    0,
                    time,
                    settings(),
                    |_| panic!("opaque"),
                    two_turns,
                )
                .unwrap();
            near(
                instantaneous.radiance_sample().radiance_w_m2_sr_nm(),
                expected.radiance_sample().radiance_w_m2_sr_nm() * 1e3,
                2e-6,
            );
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
            let tracked = timed
                .sample_ray_with_tracking_samples_during_exposure(
                    cx,
                    ray(),
                    550.0,
                    0,
                    window,
                    settings(),
                    |_| panic!("opaque"),
                    |ordinal| {
                        let [a, b, c] = two_turns(ordinal);
                        [a, b, c, 0.5]
                    },
                )
                .unwrap();
            assert_eq!(exposed, tracked);
            assert_eq!(
                exposed.mean_radiance_sample().path(),
                expected.radiance_sample().path()
            );
            assert_eq!(
                (exposed.volume_events(), exposed.free_flight_draws()),
                (2, 3)
            );
            assert_eq!(flights, [0, 1, 2]);
            near(
                exposed.integrated_radiance_j_m2_sr_nm(),
                expected.radiance_sample().radiance_w_m2_sr_nm() * area,
                2e-21,
            );
            total_energy += exposed.integrated_radiance_j_m2_sr_nm();
        }
        near(
            total_energy,
            scalar.radiance_sample().radiance_w_m2_sr_nm() * area,
            4e-21,
        );
        let intervals = pulse.intervals.lock().unwrap();
        assert_eq!(intervals.len(), 4);
        for query in intervals.iter() {
            assert_eq!(query.opening.group_delay_s, delay);
            near((query.opening.emission_time_s - epoch).to_f64(), 0.0, 1e-22);
            assert_eq!(query.duration_s, window.duration_s());
        }
    });
}
