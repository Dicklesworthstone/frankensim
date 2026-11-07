//! G1/G3/G4: real-lens shutter clocks, original paths and source-owned integrals.

use super::*;
use crate::dielectric::DielectricSurface;
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::{
    ExposureWindow, SurfaceExposureQuery,
};

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
        panic!("an exposure must use the original source interval integral")
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
fn lens_exposure_retains_each_sources_interval_and_original_path_variates() {
    with_cx(|cx, _| {
        let geometry = InstanceScene::try_new(
            [(7, 0.4), (9, 1.0)]
                .map(|(id, z)| {
                    GeometryInstance::try_new(
                        id,
                        fs_blake3::hash_domain("lens-shutter-source", &id.to_le_bytes()),
                        SharedGeometry::mesh(TriMesh::new(
                            vec![[-10.0, -10.0, z], [10.0, -10.0, z], [0.0, 10.0, z]],
                            vec![[0, 1, 2]],
                        )),
                        RigidTransform::identity(),
                    )
                    .unwrap()
                })
                .to_vec(),
        )
        .unwrap();
        let vacuum = glass(1.0, 0.0, 0.12);
        let far_glass = glass(1.5, 0.0, 0.4);
        let media = [
            SceneMedium::Homogeneous(&vacuum),
            SceneMedium::Homogeneous(&far_glass),
        ];
        let native = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![
                SurfaceBinding {
                    object_id: 7,
                    interaction: Interaction::Dielectric {
                        inside: 0,
                        outside: 1,
                        surface: DielectricSurface::SMOOTH,
                    },
                },
                SurfaceBinding {
                    object_id: 9,
                    interaction: Interaction::Opaque,
                },
            ],
            100_000,
        )
        .unwrap();
        let sources = [7, 9].map(|id| {
            SurfaceEmitter::new(
                id,
                RadianceSpectrum::try_constant(3.0).unwrap(),
                EmissionSides::Both,
            )
        });
        let scene = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let unit = Averaged::new(|_| Ok(1.0));
        let histories = [
            binding(&unit),
            SurfaceEmissionHistoryBinding {
                object_id: 9,
                ..binding(&unit)
            },
        ];
        let transient = TransientCurvedCameraScene::try_new(cx, &scene, &histories).unwrap();
        let physical = camera();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let window =
            ExposureWindow::try_new(epoch + Dd::from_f64(1e-8), epoch + Dd::from_f64(2e-8))
                .unwrap();
        for nm in [450.0, 650.0] {
            let beam = physical
                .generate_sequential_spectral_timed_beams_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    LensSample::CENTER,
                    &[nm],
                )
                .unwrap()[0]
                .unwrap();
            let mut steady_events = Vec::new();
            let steady = scene
                .sample_sequential_ray_with_samples(cx, &beam, 0, settings(), |event| {
                    steady_events.push(event);
                    [0.25, 0.5, 0.99]
                })
                .unwrap();
            let mut exposed_events = Vec::new();
            let exposed = transient
                .sample_sequential_ray_with_samples_during_exposure(
                    cx,
                    &beam,
                    0,
                    window,
                    settings(),
                    |event| {
                        exposed_events.push(event);
                        [0.25, 0.5, 0.99]
                    },
                )
                .unwrap();
            assert_eq!(*exposed.mean_radiance_sample(), steady);
            assert_eq!(steady_events, exposed_events);
            assert_eq!(exposed.window(), window);
            assert_eq!(
                exposed.integrated_radiance_j_m2_sr_nm(),
                (Dd::from_f64(steady.radiance_w_m2_sr_nm()) * window.duration_s()).to_f64()
            );
            let queries = core::mem::take(&mut *unit.queries.lock().unwrap());
            assert_eq!(queries.len(), 2);
            let group = 1.48 + 0.018 / (nm * 1e-3).powi(2);
            let lens_group_m = 0.09 + 0.01 * group;
            for (query, scene_group_m) in queries.iter().zip([0.3, 1.2]) {
                assert_eq!(query.opening.observation_time_s, window.open_s());
                assert_eq!(query.observation_close_s, window.close_s());
                assert_eq!(query.duration_s, window.duration_s());
                assert!(
                    (query.opening.group_delay_s.to_f64()
                        - (lens_group_m + scene_group_m) / LIGHT_M_S)
                        .abs()
                        < 1e-16
                );
                assert_eq!(
                    query.opening.emission_time_s,
                    window.open_s() - query.opening.group_delay_s
                );
                assert_eq!(
                    query.emission_close_s,
                    window.close_s() - query.opening.group_delay_s
                );
            }
            assert_eq!(
                (queries[0].opening.object_id, queries[1].opening.object_id),
                (7, 9)
            );
            // Seeded work follows its same branch choices and performs one trace.
            let seeded = scene
                .sample_sequential_ray(cx, &beam, 0, settings(), 13, 7)
                .unwrap();
            assert_eq!(
                *transient
                    .sample_sequential_ray_during_exposure(cx, &beam, 0, window, settings(), 13, 7,)
                    .unwrap()
                    .mean_radiance_sample(),
                seeded
            );
            unit.queries.lock().unwrap().clear();
        }
    });
}

#[test]
fn exposure_requires_timed_beams_and_interval_histories_and_preserves_cancellation() {
    with_cx(|cx, gate| {
        let geometry = geometry();
        let vacuum = glass(1.0, 0.0, 0.0);
        let media = [SceneMedium::Homogeneous(&vacuum)];
        let native = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![SurfaceBinding {
                object_id: 7,
                interaction: Interaction::Opaque,
            }],
            100_000,
        )
        .unwrap();
        let sources = [SurfaceEmitter::new(
            7,
            RadianceSpectrum::try_constant(3.0).unwrap(),
            EmissionSides::Both,
        )];
        let scene = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let instant = History::new(|_| Ok(1.0));
        let transient =
            TransientCurvedCameraScene::try_new(cx, &scene, &[binding(&instant)]).unwrap();
        let physical = camera();
        let ordinary = physical
            .generate_sequential_spectral_beams_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                &[550.0],
            )
            .unwrap()[0]
            .unwrap();
        let window = ExposureWindow::try_new(Dd::ZERO, Dd::from_f64(1e-8)).unwrap();
        assert_eq!(
            transient
                .sample_sequential_ray_during_exposure(cx, &ordinary, 0, window, settings(), 13, 7,)
                .unwrap_err(),
            TransientCameraError::MissingLensGroupDelay
        );
        let beam = physical
            .generate_sequential_spectral_timed_beams_from_tangent_offsets(
                cx,
                0.0,
                0.0,
                LensSample::CENTER,
                &[550.0],
            )
            .unwrap()[0]
            .unwrap();
        assert_eq!(
            transient
                .sample_sequential_ray_during_exposure(cx, &beam, 0, window, settings(), 13, 7,)
                .unwrap_err(),
            TransientCameraError::History {
                object_id: 7,
                error: SurfaceHistoryError::ExposureUnsupported,
            }
        );
        assert!(instant.queries.lock().unwrap().is_empty());
        let analyzed = transient
            .with_analyzer(CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.3).unwrap());
        assert_eq!(
            analyzed
                .sample_sequential_ray_during_exposure(cx, &beam, 0, window, settings(), 13, 7,)
                .unwrap_err(),
            TransientCameraError::AnalyzerMustPrecedeLens
        );
        let cancel = Averaged::new(|_| {
            gate.request();
            Ok(1.0)
        });
        let interrupted =
            TransientCurvedCameraScene::try_new(cx, &scene, &[binding(&cancel)]).unwrap();
        assert!(
            interrupted
                .sample_sequential_ray_during_exposure(cx, &beam, 0, window, settings(), 13, 7,)
                .is_err()
        );
        assert_eq!(cancel.queries.lock().unwrap().len(), 1);
    });
}
