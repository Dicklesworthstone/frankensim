use super::*;
use crate::polarization::Stokes;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::ExposureWindow;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::piecewise_linear::PiecewiseLinearHistory;

#[test]
fn g1_analyzer_keeps_native_refraction_and_source_histories_for_both_clock_owners() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(vec![plane(1, 0.4), plane(9, 1.0)]).unwrap();
        let incident = material(CauchyIor::try_constant(1.0).unwrap(), 0.1);
        let transmitted = material(CauchyIor::try_new(1.5, 0.01, 0.0).unwrap(), 0.4);
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
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let knots = [(-1.0, 1.0), (1.0, 1.0)];
        let unity = PiecewiseLinearHistory::try_new(cx, epoch, &knots).unwrap();
        let bindings = [binding(1, &unity), binding(9, &unity)];
        let scalar = TransientCurvedCameraScene::try_new(cx, &camera, &bindings).unwrap();
        assert_eq!(scalar.analyzer(), None);
        let ray = Ray {
            dir: Vec3::new(0.8, 0.0, 0.6),
            ..ray()
        };
        let identity = CameraAnalyzer::try_new(
            ray.dir,
            Stokes {
                i: 1.0,
                q: 0.0,
                u: 0.0,
                v: 0.0,
            },
        )
        .unwrap();
        let selector = CameraAnalyzer::try_new(
            Vec3::new(0.0, 1.0, 0.0),
            Stokes {
                i: 0.5,
                q: 0.2,
                u: -0.15,
                v: 0.1,
            },
        )
        .unwrap();
        let time = epoch + Dd::from_f64(1e-8);
        let window = ExposureWindow::try_new(time, time + Dd::from_f64(1e-8)).unwrap();
        for analyzer in [identity, selector] {
            let timed = TransientCurvedCameraScene::try_new(cx, &camera, &bindings)
                .unwrap()
                .with_analyzer(analyzer);
            assert_eq!(timed.analyzer(), Some(analyzer));
            assert!(core::ptr::eq(timed.source_camera(), &camera));
            assert_eq!(timed.history_identity(9), scalar.history_identity(9));
            for wavelength in [450.0, 650.0] {
                for branch in [0.0, 0.99] {
                    let mut original_events = Vec::new();
                    let original = camera
                        .sample_ray_analyzed_with_samples(
                            cx,
                            ray,
                            wavelength,
                            0,
                            analyzer,
                            settings(),
                            |event| {
                                original_events.push(event);
                                [0.25, 0.5, branch]
                            },
                        )
                        .unwrap();
                    let mut observed_events = Vec::new();
                    let observed = timed
                        .sample_ray_with_samples_at_time(
                            cx,
                            ray,
                            wavelength,
                            0,
                            time,
                            settings(),
                            |event| {
                                observed_events.push(event);
                                [0.25, 0.5, branch]
                            },
                        )
                        .unwrap();
                    let exposed = timed
                        .sample_ray_with_samples_during_exposure(
                            cx,
                            ray,
                            wavelength,
                            0,
                            window,
                            settings(),
                            |_| [0.25, 0.5, branch],
                        )
                        .unwrap();
                    assert_eq!(observed.radiance_sample(), &original);
                    assert_eq!(exposed.mean_radiance_sample(), &original);
                    assert_eq!(original_events, observed_events);
                    near(
                        exposed.integrated_radiance_j_m2_sr_nm(),
                        original.radiance_w_m2_sr_nm() * window.duration_s().to_f64(),
                        1e-22,
                    );
                    let unselected = scalar
                        .sample_ray_with_samples_at_time(
                            cx,
                            ray,
                            wavelength,
                            0,
                            time,
                            settings(),
                            |_| [0.25, 0.5, branch],
                        )
                        .unwrap();
                    if analyzer == identity {
                        assert_eq!(observed, unselected);
                    } else if branch == 0.99 {
                        // Oblique Fresnel diattenuation couples to Q before
                        // the retarded sources. A scalar half-filter is wrong.
                        assert!(
                            (original.radiance_w_m2_sr_nm()
                                - 0.5 * unselected.radiance_sample().radiance_w_m2_sr_nm())
                            .abs()
                                > 1e-4
                        );
                    }
                }
                for path in 0..3 {
                    let original = camera
                        .sample_ray_analyzed(cx, ray, wavelength, 0, analyzer, settings(), 29, path)
                        .unwrap();
                    let observed = timed
                        .sample_ray_at_time(cx, ray, wavelength, 0, time, settings(), 29, path)
                        .unwrap();
                    let exposed = timed
                        .sample_ray_during_exposure(
                            cx,
                            ray,
                            wavelength,
                            0,
                            window,
                            settings(),
                            29,
                            path,
                        )
                        .unwrap();
                    assert_eq!(observed.radiance_sample(), &original);
                    assert_eq!(exposed.mean_radiance_sample(), &original);
                }
            }
        }
    });
}
