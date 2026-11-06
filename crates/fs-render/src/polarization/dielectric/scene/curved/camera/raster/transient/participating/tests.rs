mod exposure;

use super::super::super::tests::{bindings, camera, settings, with_cx};
use super::*;
use crate::charts::TriMesh;
use crate::dielectric::{BeerLambertAbsorption, CauchyIor, DielectricGlass, GlassProvenance};
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::polarization::dielectric::scene::PathTermination;
use crate::polarization::dielectric::scene::curved::camera::participating::CurvedRayleighCameraScene;
use crate::polarization::dielectric::scene::curved::camera::transient::{
    SurfaceEmissionHistory, SurfaceEmissionHistoryBinding, SurfaceEmissionQuery,
    SurfaceHistoryError, TransientCameraError, TransientCurvedCameraScene,
};
use crate::polarization::dielectric::scene::curved::camera::{
    CurvedCameraScene, EmissionSides, RadianceSpectrum, SurfaceEmitter,
};
use crate::polarization::dielectric::scene::curved::participating::{
    CurvedRayleighMedium, CurvedRayleighScene,
};
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, SceneMedium};
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};
use fs_blake3::ContentHash;
use std::sync::atomic::{AtomicUsize, Ordering};

// The closed, emitting box receives paths turned in any molecular direction.
fn geometry() -> InstanceScene {
    InstanceScene::try_new(vec![
        GeometryInstance::try_new(
            7,
            ContentHash([89; 32]),
            SharedGeometry::mesh(TriMesh::new(
                vec![
                    [-2.0, -2.0, -2.0],
                    [2.0, -2.0, -2.0],
                    [2.0, 2.0, -2.0],
                    [-2.0, 2.0, -2.0],
                    [-2.0, -2.0, 2.0],
                    [2.0, -2.0, 2.0],
                    [2.0, 2.0, 2.0],
                    [-2.0, 2.0, 2.0],
                ],
                vec![
                    [0, 2, 1],
                    [0, 3, 2],
                    [4, 5, 6],
                    [4, 6, 7],
                    [0, 1, 5],
                    [0, 5, 4],
                    [3, 7, 6],
                    [3, 6, 2],
                    [0, 4, 7],
                    [0, 7, 3],
                    [1, 2, 6],
                    [1, 6, 5],
                ],
            )),
            RigidTransform::identity(),
        )
        .unwrap(),
    ])
    .unwrap()
}

fn glass() -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_new(1.4, 0.02, 0.0).unwrap(),
        BeerLambertAbsorption::try_constant(0.12).unwrap(),
        GlassProvenance::Custom,
    )
}

fn sources() -> [SurfaceEmitter; 1] {
    [SurfaceEmitter::new(
        7,
        RadianceSpectrum::try_constant(3.0).unwrap(),
        EmissionSides::Both,
    )]
}

fn trace_settings() -> CurvedRayleighTraceSettings {
    let mut surfaces = settings();
    surfaces.max_distance_m = 128.0;
    CurvedRayleighTraceSettings {
        surfaces,
        max_volume_events: 64,
    }
}

struct History<F> {
    evaluate: F,
    calls: AtomicUsize,
}
impl<F> History<F> {
    fn new(evaluate: F) -> Self {
        Self {
            evaluate,
            calls: AtomicUsize::new(0),
        }
    }
    fn binding(&self) -> SurfaceEmissionHistoryBinding<'_>
    where
        F: Fn(SurfaceEmissionQuery, usize) -> Result<f64, SurfaceHistoryError> + Sync,
    {
        SurfaceEmissionHistoryBinding {
            object_id: 7,
            source_identity: ContentHash([87; 32]),
            history: self,
        }
    }
}
impl<F> SurfaceEmissionHistory for History<F>
where
    F: Fn(SurfaceEmissionQuery, usize) -> Result<f64, SurfaceHistoryError> + Sync,
{
    fn sample(&self, cx: &Cx<'_>, query: SurfaceEmissionQuery) -> Result<f64, SurfaceHistoryError> {
        cx.checkpoint()?;
        (self.evaluate)(query, self.calls.fetch_add(1, Ordering::Relaxed))
    }
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-11 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g3_actual_scattering_preserves_static_pixels_rng_and_channel_shared_work() {
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
        let unity = History::new(|_, _| Ok(1.0));
        let transient =
            TransientCurvedCameraScene::try_new(cx, &surface, &[unity.binding()]).unwrap();
        let timed =
            TransientCurvedRayleighCameraScene::try_new(cx, &transient, &scattering).unwrap();
        let optical = camera(0.01);
        let raster = CurvedCameraRaster::try_new(&optical, 3, 2).unwrap();
        let clock = Dd::from_f64(1e-6);
        let expected = raster
            .sample_pixel_rayleigh(cx, &steady, 0, trace_settings(), 73, 4, 9)
            .unwrap();
        let actual = raster
            .sample_pixel_rayleigh_transient_at_time(
                cx,
                &timed,
                0,
                clock,
                trace_settings(),
                73,
                4,
                9,
            )
            .unwrap();
        assert_eq!(actual.pixel_sample(), expected.pixel_sample());
        assert_eq!(actual.volume_events(), expected.volume_events());
        assert_eq!(actual.free_flight_draws(), expected.free_flight_draws());
        assert_eq!(actual.observation_time_s(), clock);
        // Lane 3's native first flight is 1.473589... m, inside the 2 m box.
        assert!(actual.volume_events()[3] >= 1);
        assert_eq!(unity.calls.load(Ordering::Relaxed), 4);
        raster
            .sample_pixel_rayleigh_transient_at_time(
                cx,
                &timed,
                0,
                clock,
                trace_settings(),
                73,
                1,
                3,
            )
            .unwrap();
        assert_eq!(
            actual,
            raster
                .sample_pixel_rayleigh_transient_at_time(
                    cx,
                    &timed,
                    0,
                    clock,
                    trace_settings(),
                    73,
                    4,
                    9
                )
                .unwrap()
        );

        let one = QuantumEfficiencySensor::try_new(vec![(420.0, [0.5]), (680.0, [0.5])]).unwrap();
        let three = QuantumEfficiencySensor::try_new(vec![
            (420.0, [0.2, 0.3, 0.0]),
            (680.0, [0.2, 0.3, 0.0]),
        ])
        .unwrap();
        let before = unity.calls.load(Ordering::Relaxed);
        let a = raster
            .sample_pixel_sensor_rayleigh_transient_at_time(
                cx,
                &timed,
                &one,
                0,
                clock,
                trace_settings(),
                73,
                4,
                9,
            )
            .unwrap();
        let b = raster
            .sample_pixel_sensor_rayleigh_transient_at_time(
                cx,
                &timed,
                &three,
                0,
                clock,
                trace_settings(),
                73,
                4,
                9,
            )
            .unwrap();
        assert_eq!(unity.calls.load(Ordering::Relaxed) - before, 8);
        assert_eq!(a.volume_events(), b.volume_events());
        assert_eq!(a.free_flight_draws(), b.free_flight_draws());
        assert_eq!(
            a.pixel_sample().spectral_samples(),
            b.pixel_sample().spectral_samples()
        );
        let expected = raster
            .sample_pixel_sensor_rayleigh(cx, &steady, &three, 0, trace_settings(), 73, 4, 9)
            .unwrap();
        assert_eq!(b.pixel_sample(), expected.pixel_sample());
        near(
            a.pixel_sample().photoelectron_rate_per_m2_sr()[0],
            b.pixel_sample().photoelectron_rate_per_m2_sr()[0]
                + b.pixel_sample().photoelectron_rate_per_m2_sr()[1],
        );
        assert_eq!(b.pixel_sample().photoelectron_rate_per_m2_sr()[2], 0.0);
    });
}

#[test]
fn g1_qe_reads_chromatic_source_times_after_actual_molecular_travel() {
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
        let ramp =
            History::new(
                |query: SurfaceEmissionQuery, _| Ok(query.emission_time_s.to_f64() / 1e-6),
            );
        let transient =
            TransientCurvedCameraScene::try_new(cx, &surface, &[ramp.binding()]).unwrap();
        let timed =
            TransientCurvedRayleighCameraScene::try_new(cx, &transient, &scattering).unwrap();
        let optical = camera(0.01);
        let raster = CurvedCameraRaster::try_new(&optical, 3, 2).unwrap();
        let sensor =
            QuantumEfficiencySensor::try_new(vec![(420.0, [0.5]), (680.0, [0.5])]).unwrap();
        let original = raster
            .sample_pixel_sensor_rayleigh(cx, &steady, &sensor, 0, trace_settings(), 73, 4, 9)
            .unwrap();
        let result = raster
            .sample_pixel_sensor_rayleigh_transient_at_time(
                cx,
                &timed,
                &sensor,
                0,
                Dd::from_f64(1e-6),
                trace_settings(),
                73,
                4,
                9,
            )
            .unwrap();
        let pixel = result.pixel_sample();
        let mut energy = 0.0;
        let mut electrons = 0.0;
        for lane in 0..CAMERA_WAVELENGTHS {
            let draw = pixel.detector_wavelengths()[lane];
            let wave_um = draw.wavelength_nm() * 1e-3;
            let source = original.pixel_sample().spectral_samples()[lane];
            // Independent homogeneous-medium group index and the retained
            // scattered arclength determine this linear history exactly.
            let delay = source.path().traveled_m() * (1.4 + 0.06 / wave_um.powi(2)) / 299_792_458.0;
            let radiance = source.radiance_w_m2_sr_nm() * (1.0 - delay / 1e-6);
            near(
                pixel.spectral_samples()[lane].radiance_w_m2_sr_nm(),
                radiance,
            );
            let weighted = radiance * pixel.diffraction_weights()[lane] / 4.0;
            energy += weighted * draw.radiance_weights_nm()[0];
            electrons += weighted * draw.electron_weights_nm_per_joule()[0];
        }
        near(pixel.qe_weighted_radiance_w_m2_sr()[0], energy);
        near(pixel.photoelectron_rate_per_m2_sr()[0], electrons);
        assert!(result.volume_events()[3] >= 1);
        assert_eq!(ramp.calls.load(Ordering::Relaxed), 4);
    });
}

#[test]
fn g0_late_source_failure_volume_limit_and_cancellation_refuse_whole_pixels() {
    with_cx(|cx, gate| {
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
        let failing = History::new(|_, call| {
            if call == 3 {
                Err(SurfaceHistoryError::OutsideTimeDomain)
            } else {
                Ok(1.0)
            }
        });
        let transient =
            TransientCurvedCameraScene::try_new(cx, &surface, &[failing.binding()]).unwrap();
        let timed =
            TransientCurvedRayleighCameraScene::try_new(cx, &transient, &scattering).unwrap();
        let optical = camera(0.01);
        let raster = CurvedCameraRaster::try_new(&optical, 3, 2).unwrap();
        let sensor =
            QuantumEfficiencySensor::try_new(vec![(420.0, [0.5]), (680.0, [0.5])]).unwrap();
        let clock = Dd::from_f64(1e-6);
        assert!(matches!(
            raster.sample_pixel_sensor_rayleigh_transient_at_time(
                cx,
                &timed,
                &sensor,
                0,
                clock,
                trace_settings(),
                73,
                4,
                9
            ),
            Err(TransientRasterError::Source(
                TransientCameraError::History {
                    object_id: 7,
                    error: SurfaceHistoryError::OutsideTimeDomain
                }
            ))
        ));
        assert_eq!(failing.calls.load(Ordering::Relaxed), 4);
        let mut limited = trace_settings();
        limited.max_volume_events = 1;
        assert!(matches!(
            raster.sample_pixel_rayleigh_transient_at_time(cx, &timed, 0, clock, limited, 73, 4, 9),
            Err(TransientRasterError::Source(TransientCameraError::Camera(
                CameraRadianceError::Unresolved(PathTermination::VolumeLimit)
            )))
        ));
        gate.request();
        assert!(
            raster
                .sample_pixel_rayleigh_transient_at_time(
                    cx,
                    &timed,
                    0,
                    clock,
                    trace_settings(),
                    73,
                    4,
                    9
                )
                .is_err()
        );
    });
}
