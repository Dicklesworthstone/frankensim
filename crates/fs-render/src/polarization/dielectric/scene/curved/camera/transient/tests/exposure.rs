//! G1/G3/G4: original paths, exact temporal integrals and refused partial work.

use super::super::exposure::piecewise_linear::{
    PiecewiseLinearHistory, PiecewiseLinearHistoryError,
};
use super::super::exposure::{ExposureWindow, SurfaceExposureQuery};
use super::*;

struct Averaged<F> {
    evaluate: F,
    queries: Mutex<Vec<SurfaceExposureQuery>>,
}

impl<F> Averaged<F> {
    fn new(evaluate: F) -> Self {
        Self {
            evaluate,
            queries: Mutex::new(Vec::new()),
        }
    }
}

impl<F: Fn(SurfaceExposureQuery) -> Result<f64, SurfaceHistoryError> + Sync> SurfaceEmissionHistory
    for Averaged<F>
{
    fn sample(&self, _: &Cx<'_>, _: SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError> {
        panic!("finite exposure must never replace source integration with point samples")
    }
    fn average(
        &self,
        cx: &Cx<'_>,
        query: SurfaceExposureQuery,
    ) -> Result<f64, SurfaceHistoryError> {
        cx.checkpoint()?;
        self.queries.lock().unwrap().push(query);
        (self.evaluate)(query)
    }
}

#[test]
fn g1_piecewise_source_integrates_narrow_interior_pulse_on_original_epoch() {
    with_cx(|cx, _| {
        let knots = [
            (0.0, 0.0),
            (1e-9, 0.0),
            (1.001e-9, 5e5),
            (1.002e-9, 0.0),
            (1.0, 0.0),
        ];
        let area = 0.5 * 5e5 * (knots[3].0 - knots[1].0);
        for epoch in [Dd::ZERO, Dd::from_f64(2.0_f64.powi(30))] {
            let history = PiecewiseLinearHistory::try_new(cx, epoch, &knots).unwrap();
            assert!(core::ptr::eq(history.knots(), knots.as_slice()));
            near(
                history.mean_between(cx, epoch, epoch + Dd::ONE).unwrap(),
                area,
                1e-20,
            );
            assert_eq!(history.value_at(cx, epoch).unwrap(), 0.0);
            // Both clipped endpoints lie between original knots. The mean
            // of the clipped triangular cap is 3/4 of its peak.
            let clipped = history
                .mean_between(
                    cx,
                    epoch + Dd::from_f64(1.0005e-9),
                    epoch + Dd::from_f64(1.0015e-9),
                )
                .unwrap();
            near(clipped, 375_000.0, 0.01);
        }
        // The pulse duration fraction underflows binary64, but its weighted
        // mean remains representable. No cumulative-area subtraction occurs.
        let wide = [(0.0, 0.0), (1e-250, 1e250), (2e-250, 0.0), (1e100, 0.0)];
        let history = PiecewiseLinearHistory::try_new(cx, Dd::ZERO, &wide).unwrap();
        near(
            history
                .mean_between(cx, Dd::ZERO, Dd::from_f64(1e100))
                .unwrap(),
            1e-100,
            1e-114,
        );
        // The tiny duration fraction also occurs for an enormous window
        // whose triangle area overflows. Its complete mean is still three.
        let huge = [(0.0, 0.0), (3.0, f64::MAX), (6.0, 0.0), (f64::MAX, 0.0)];
        let history = PiecewiseLinearHistory::try_new(cx, Dd::ZERO, &huge).unwrap();
        near(
            history
                .mean_between(cx, Dd::ZERO, Dd::from_f64(f64::MAX))
                .unwrap(),
            3.0,
            1e-14,
        );
    });
}

#[test]
fn g3_unit_exposure_preserves_trace_variates_and_each_sources_retarded_window() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(1, 0.4), plane(9, 1.0)]).unwrap();
        let incident = material(CauchyIor::try_constant(1.0).unwrap(), 0.1);
        let transmitted = material(CauchyIor::try_constant(1.5).unwrap(), 0.4);
        let media = [
            SceneMedium::Homogeneous(&incident),
            SceneMedium::Homogeneous(&transmitted),
        ];
        let native = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![interface(), opaque(9)],
            100_000,
        )
        .unwrap();
        let sources = [source(1, 2.0), source(9, 5.0)];
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let unit = Averaged::new(|_| Ok(1.0));
        let timed = TransientCurvedCameraScene::try_new(
            cx,
            &camera,
            &[binding(1, &unit), binding(9, &unit)],
        )
        .unwrap();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let window =
            ExposureWindow::try_new(epoch + Dd::from_f64(1e-8), epoch + Dd::from_f64(2e-8))
                .unwrap();
        let mut original_events = Vec::new();
        let original = camera
            .sample_ray_with_samples(cx, ray(), 550.0, 0, settings(), |event| {
                original_events.push(event);
                [0.25, 0.5, 0.99]
            })
            .unwrap();
        let mut exposed_events = Vec::new();
        let exposed = timed
            .sample_ray_with_samples_during_exposure(
                cx,
                ray(),
                550.0,
                0,
                window,
                settings(),
                |event| {
                    exposed_events.push(event);
                    [0.25, 0.5, 0.99]
                },
            )
            .unwrap();
        assert_eq!(*exposed.mean_radiance_sample(), original);
        assert_eq!(original_events, exposed_events);
        near(
            exposed.integrated_radiance_j_m2_sr_nm(),
            original.radiance_w_m2_sr_nm() * 1e-8,
            1e-21,
        );
        let queries = unit.queries.lock().unwrap();
        assert_eq!(queries.len(), 2);
        assert_eq!(
            (queries[0].opening.object_id, queries[1].opening.object_id),
            (1, 9)
        );
        for (query, group_m) in queries.iter().zip([0.4, 0.4 + 0.6 * 1.5]) {
            near(
                query.opening.group_delay_s.to_f64(),
                group_m / 299_792_458.0,
                1e-16,
            );
            assert_eq!(query.opening.observation_time_s, window.open_s());
            assert_eq!(query.observation_close_s, window.close_s());
            near(
                (query.emission_close_s - query.opening.emission_time_s).to_f64(),
                1e-8,
                1e-22,
            );
        }
        drop(queries);
        for path in 0..4 {
            let original = camera
                .sample_ray(cx, ray(), 550.0, 0, settings(), 31, path)
                .unwrap();
            let exposed = timed
                .sample_ray_during_exposure(cx, ray(), 550.0, 0, window, settings(), 31, path)
                .unwrap();
            assert_eq!(*exposed.mean_radiance_sample(), original);
        }
    });
}

#[test]
fn g4_window_and_original_piecewise_domain_refuse_without_extrapolation() {
    with_cx(|cx, gate| {
        for (open, close) in [(1.0, 1.0), (2.0, 1.0), (0.0, f64::INFINITY)] {
            assert!(ExposureWindow::try_new(Dd::from_f64(open), Dd::from_f64(close)).is_err());
        }
        for knots in [
            &[(0.0, 1.0)][..],
            &[(0.0, 1.0), (0.0, 2.0)],
            &[(0.0, -1.0), (1.0, 1.0)],
        ] {
            assert_eq!(
                PiecewiseLinearHistory::try_new(cx, Dd::ZERO, knots).unwrap_err(),
                PiecewiseLinearHistoryError::InvalidKnots
            );
        }
        let knots = [(0.0, 0.0), (1.0, 2.0)];
        let history = PiecewiseLinearHistory::try_new(cx, Dd::ZERO, &knots).unwrap();
        assert_eq!(history.value_at(cx, Dd::ONE).unwrap(), 2.0);
        assert_eq!(history.mean_between(cx, Dd::ZERO, Dd::ONE).unwrap(), 1.0);
        assert_eq!(
            history.mean_between(cx, Dd::from_f64(-1e-6), Dd::ONE),
            Err(SurfaceHistoryError::OutsideTimeDomain)
        );
        assert_eq!(
            history.mean_between(cx, Dd::ZERO, Dd::ZERO),
            Err(SurfaceHistoryError::NumericalRange)
        );
        gate.request();
        assert_eq!(
            history.mean_between(cx, Dd::ZERO, Dd::ONE),
            Err(SurfaceHistoryError::Cancelled)
        );
    });
}

#[test]
fn g4_unqualified_source_and_late_failed_integral_publish_no_exposure() {
    with_cx(|cx, gate| {
        let geometry = InstanceScene::try_new(vec![plane(1, 0.4), plane(9, 1.0)]).unwrap();
        let glass = material(CauchyIor::try_constant(1.3).unwrap(), 0.0);
        let media = [
            SceneMedium::Homogeneous(&glass),
            SceneMedium::Homogeneous(&glass),
        ];
        let native = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![interface(), opaque(9)],
            100_000,
        )
        .unwrap();
        let sources = [source(1, 2.0), source(9, 5.0)];
        let camera = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let window = ExposureWindow::try_new(Dd::ZERO, Dd::from_f64(1e-8)).unwrap();
        let unit = Averaged::new(|_| Ok(1.0));
        let instantaneous = History::new(|_| Ok(1.0));
        let unqualified = TransientCurvedCameraScene::try_new(
            cx,
            &camera,
            &[binding(1, &unit), binding(9, &instantaneous)],
        )
        .unwrap();
        assert!(matches!(
            unqualified.sample_ray_during_exposure(cx, ray(), 550.0, 0, window, settings(), 3, 0),
            Err(TransientCameraError::History {
                object_id: 9,
                error: SurfaceHistoryError::ExposureUnsupported
            })
        ));
        assert!(instantaneous.queries.lock().unwrap().is_empty());
        for value in [-1.0, f64::NAN, f64::INFINITY] {
            let invalid = Averaged::new(|_| Ok(value));
            let bad = TransientCurvedCameraScene::try_new(
                cx,
                &camera,
                &[binding(1, &unit), binding(9, &invalid)],
            )
            .unwrap();
            assert!(matches!(
                bad.sample_ray_during_exposure(cx, ray(), 550.0, 0, window, settings(), 3, 0),
                Err(TransientCameraError::InvalidHistoryScale { object_id: 9 })
            ));
        }
        let cancel = Averaged::new(|_| {
            gate.request();
            Ok(1.0)
        });
        let interrupted = TransientCurvedCameraScene::try_new(
            cx,
            &camera,
            &[binding(1, &unit), binding(9, &cancel)],
        )
        .unwrap();
        assert!(
            interrupted
                .sample_ray_during_exposure(cx, ray(), 550.0, 0, window, settings(), 3, 0)
                .is_err()
        );
        assert_eq!(cancel.queries.lock().unwrap().len(), 1);
    });
}
