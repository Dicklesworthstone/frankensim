use super::*;
use crate::camera::{
    Aperture, CameraProjection, LensSample, LensSurface, PhysicalCamera, SequentialLens,
};
use crate::charts::TriMesh;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::gradient_index::GradientIndexOptions;
use crate::gradient_index::event::GradientEventOptions;
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::polarization::dielectric::scene::curved::camera::{RadianceSpectrum, SurfaceEmitter};
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, SceneMedium};
use crate::polarization::dielectric::scene::{Interaction, SurfaceBinding};
use fs_exec::{CancelGate, ExecMode, StreamKey};
use fs_geom::Vec3;
use std::sync::Mutex;

mod exposure;

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 17,
                kernel_id: 3,
                tile: 0,
                iteration: 0,
            },
            asupersync::types::Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx, &gate);
    });
}

fn glass(a: f64, b: f64, extinction: f64) -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_new(a, b, 0.0).unwrap(),
        BeerLambertAbsorption::try_constant(extinction).unwrap(),
        GlassProvenance::Custom,
    )
}

fn camera() -> PhysicalCamera {
    PhysicalCamera::try_look_at(
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(0.0, 0.0, 1.0),
        Vec3::new(0.0, 1.0, 0.0),
        CameraProjection::try_focal_sensor(0.05, 0.002).unwrap(),
        1.0,
        Aperture::try_circular(0.01).unwrap(),
    )
    .unwrap()
    .with_sequential_lens(
        SequentialLens::try_new(
            550.0,
            &[
                LensSurface::try_new(0.08, 20.0, 0.01, Some(glass(1.48, 0.006, 2.0))).unwrap(),
                LensSurface::try_new(0.09, -20.0, 0.012, None).unwrap(),
                LensSurface::try_new(0.1, 0.0, 0.02, None).unwrap(),
            ],
        )
        .unwrap(),
    )
    .unwrap()
}

fn geometry() -> InstanceScene {
    InstanceScene::try_new(vec![
        GeometryInstance::try_new(
            7,
            ContentHash([81; 32]),
            SharedGeometry::mesh(TriMesh::new(
                vec![[-10.0, -10.0, 1.0], [10.0, -10.0, 1.0], [0.0, 10.0, 1.0]],
                vec![[0, 1, 2]],
            )),
            RigidTransform::identity(),
        )
        .unwrap(),
    ])
    .unwrap()
}

fn settings() -> CurvedTraceSettings {
    CurvedTraceSettings {
        max_surface_events: 8,
        max_distance_m: 2.0,
        max_total_steps: 16384,
        max_query_work: 1000000,
        // GradientEventOptions::validate admits a surface tolerance only at eight or more
        // integration position tolerances (gradient_index/event.rs); the default 1e-9 m position
        // tolerance with this 1e-9 m surface tolerance was refused as InvalidOptions before any
        // trace. 1e-10 m keeps that ratio at ten.
        event: GradientEventOptions {
            integration: GradientIndexOptions {
                position_tolerance_m: 1e-10,
                ..GradientIndexOptions::default()
            },
            surface_tolerance_m: 1e-9,
            ..GradientEventOptions::default()
        },
    }
}

struct History<F> {
    sample: F,
    queries: Mutex<Vec<SurfaceEmissionQuery>>,
}
impl<F> History<F> {
    fn new(sample: F) -> Self {
        Self {
            sample,
            queries: Mutex::new(Vec::new()),
        }
    }
}
impl<F: Fn(SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError> + Sync> SurfaceEmissionHistory
    for History<F>
{
    fn sample(&self, cx: &Cx<'_>, query: SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError> {
        cx.checkpoint()?;
        self.queries.lock().unwrap().push(query);
        (self.sample)(query)
    }
}

fn binding(history: &dyn SurfaceEmissionHistory) -> SurfaceEmissionHistoryBinding<'_> {
    SurfaceEmissionHistoryBinding {
        object_id: 7,
        source_identity: ContentHash([82; 32]),
        history,
    }
}

#[test]
fn real_lens_delay_keeps_sensor_clock_and_original_steady_radiance_at_large_epoch() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let vacuum = glass(1.0, 0.0, 0.12);
        let media = [SceneMedium::Homogeneous(&vacuum)];
        let native = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![SurfaceBinding {
                object_id: 7,
                interaction: Interaction::Opaque,
            }],
            100000,
        )
        .unwrap();
        let sources = [SurfaceEmitter::new(
            7,
            RadianceSpectrum::try_constant(3.0).unwrap(),
            EmissionSides::Both,
        )];
        let scene = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let history = History::new(|_| Ok(1.0));
        let transient =
            TransientCurvedCameraScene::try_new(cx, &scene, &[binding(&history)]).unwrap();
        let camera = camera();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let clock = epoch + Dd::from_f64(1e-8);
        let mut delays = Vec::new();
        for wavelength in [450.0, 650.0] {
            let optical = camera
                .generate_sequential_spectral_timed_beams_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    LensSample::CENTER,
                    &[wavelength],
                )
                .unwrap()[0]
                .unwrap();
            // On axis: 9 cm in vacuum and 1 cm in Cauchy glass. Group
            // index is A+3B/lambda_um², not the phase A+B/lambda_um².
            let group = 1.48 + 0.018 / (wavelength * 1e-3).powi(2);
            let expected_lens = (0.09 + 0.01 * group) / LIGHT_M_S;
            assert!((optical.group_delay_s().unwrap().to_f64() - expected_lens).abs() < 1e-23);
            let steady = scene
                .sample_sequential_ray_with_samples(cx, &optical, 0, settings(), |_| {
                    [0.2, 0.5, 0.9]
                })
                .unwrap();
            let observed = transient
                .sample_sequential_ray_with_samples_at_time(
                    cx,
                    &optical,
                    0,
                    clock,
                    settings(),
                    |_| [0.2, 0.5, 0.9],
                )
                .unwrap();
            assert_eq!(*observed.radiance_sample(), steady);
            assert_eq!(observed.observation_time_s(), clock);
            let query = *history.queries.lock().unwrap().last().unwrap();
            let total = optical.group_delay_s().unwrap()
                + Dd::from_f64(steady.path().group_optical_path_length_m().unwrap())
                    / Dd::from_f64(LIGHT_M_S);
            assert_eq!(query.observation_time_s, clock);
            assert_eq!(query.group_delay_s, total);
            assert_eq!(query.emission_time_s, clock - total);
            assert_eq!(query.position_m, steady.path().beam().ray().origin);
            assert!(
                ((query.emission_time_s - epoch).to_f64()
                    - (1e-8 - (0.99 + 0.01 * group) / LIGHT_M_S))
                    .abs()
                    < 1e-17
            );
            delays.push(total.to_f64());
            assert_eq!(
                transient
                    .sample_sequential_ray_at_time(cx, &optical, 0, clock, settings(), 13, 7)
                    .unwrap(),
                observed
            );
        }
        assert!((delays[0] - delays[1]).abs() > 1e-13);
    });
}

#[test]
fn timing_and_analyzer_admission_preserve_original_source_errors_and_cancellation() {
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
            100000,
        )
        .unwrap();
        let sources = [SurfaceEmitter::new(
            7,
            RadianceSpectrum::try_constant(3.0).unwrap(),
            EmissionSides::Both,
        )];
        let scene = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let history = History::new(|_| Err(SurfaceHistoryError::OutsideTimeDomain));
        let transient =
            TransientCurvedCameraScene::try_new(cx, &scene, &[binding(&history)]).unwrap();
        let camera = camera();
        let ordinary = camera
            .generate_sequential_spectral_beams_from_tangent_offsets(
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
                .sample_sequential_ray_at_time(cx, &ordinary, 0, Dd::ZERO, settings(), 13, 7)
                .unwrap_err(),
            TransientCameraError::MissingLensGroupDelay
        );
        assert!(history.queries.lock().unwrap().is_empty());
        let timed = camera
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
                .sample_sequential_ray_at_time(cx, &timed, 0, Dd::ZERO, settings(), 13, 7)
                .unwrap_err(),
            TransientCameraError::History {
                object_id: 7,
                error: SurfaceHistoryError::OutsideTimeDomain
            }
        );
        let analyzed = transient
            .with_analyzer(CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.3).unwrap());
        assert_eq!(
            analyzed
                .sample_sequential_ray_at_time(cx, &timed, 0, Dd::ZERO, settings(), 13, 7)
                .unwrap_err(),
            TransientCameraError::AnalyzerMustPrecedeLens
        );
        assert_eq!(history.queries.lock().unwrap().len(), 1);
        gate.request();
        assert!(
            analyzed
                .sample_sequential_ray_at_time(cx, &timed, 0, Dd::ZERO, settings(), 13, 7)
                .is_err()
        );
    });
}
