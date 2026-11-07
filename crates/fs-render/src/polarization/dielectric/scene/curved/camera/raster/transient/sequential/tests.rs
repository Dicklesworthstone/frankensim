use super::super::super::tests::{bindings, geometry, settings, with_cx};
use super::*;
use crate::camera::{
    Aperture, CameraError, CameraProjection, LensSample, LensSurface, PhysicalCamera,
    SequentialLens,
};
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricIorSample, GlassProvenance,
    TabulatedIor,
};
use crate::polarization::dielectric::scene::curved::camera::transient::{
    SurfaceEmissionHistory, SurfaceEmissionHistoryBinding, SurfaceEmissionQuery,
    SurfaceHistoryError,
};
use crate::polarization::dielectric::scene::curved::camera::{
    CameraAnalyzer, CurvedCameraScene, EmissionSides, RadianceSpectrum, SurfaceEmitter,
};
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, SceneMedium};
use fs_blake3::ContentHash;
use fs_geom::{Point3, Vec3};
use std::sync::Mutex;

mod exposure;

const SEED: u64 = 79;
const SAMPLE: u32 = 6;
const C: f64 = 299792458.0;

fn glass() -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_new(1.48, 0.006, 0.0).unwrap(),
        BeerLambertAbsorption::try_constant(2.0).unwrap(),
        GlassProvenance::Custom,
    )
}
fn camera(material: DielectricGlass, reference: f64, stop: Option<f64>) -> PhysicalCamera {
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
    nm: &[f64; 4],
    analyzer: Option<CameraAnalyzer>,
) -> [Option<SequentialCameraBeam>; 4] {
    let key = [SEED as u32, 0];
    let raster = philox4x32_10([0, SAMPLE, 0, 0x6372_7331], key);
    let pupil = philox4x32_10([0, SAMPLE, 0, 0x6372_7332], key);
    let half = camera.raster_vertical_half_tan();
    let x = (2.0 * unit(raster[0]) - 1.0) * half;
    let y = (1.0 - 2.0 * unit(raster[1])) * half;
    let lens = LensSample::try_new(unit(pupil[0]), unit(pupil[1])).unwrap();
    if let Some(a) = analyzer {
        camera
            .generate_sequential_spectral_timed_beams_analyzed_from_tangent_offsets(
                cx,
                x,
                y,
                lens,
                nm,
                a.reference_axis(),
                a.row(),
            )
            .unwrap()
    } else {
        camera
            .generate_sequential_spectral_timed_beams_from_tangent_offsets(cx, x, y, lens, nm)
            .unwrap()
    }
}
fn vacuum() -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_constant(1.0).unwrap(),
        BeerLambertAbsorption::try_constant(0.12).unwrap(),
        GlassProvenance::Custom,
    )
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
        source_identity: ContentHash([84; 32]),
        history,
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}");
}

#[test]
fn sensor_analyzer_and_clipped_packets_keep_original_measure_and_dd_clock() {
    with_cx(|cx, gate| {
        let geometry = geometry();
        let vacuum = vacuum();
        let media = [SceneMedium::Homogeneous(&vacuum)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100000).unwrap();
        let sources = [SurfaceEmitter::new(
            7,
            RadianceSpectrum::try_constant(3.0).unwrap(),
            EmissionSides::Both,
        )];
        let scene = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let history = History::new(|_| Ok(1.0));
        let analyzer = CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), 0.3).unwrap();
        let timed = TransientCurvedCameraScene::try_new(cx, &scene, &[binding(&history)])
            .unwrap()
            .with_analyzer(analyzer);
        let wide = camera(glass(), 550.0, None);
        let full = CurvedCameraRaster::try_new(&wide, 1, 1)
            .unwrap()
            .sample_pixel_sequential(cx, &scene, 0, settings(), SEED, 0, SAMPLE)
            .unwrap();
        let mut radii = packet(&wide, cx, full.wavelengths_nm(), None).map(|beam| {
            let ray = beam.unwrap().ray();
            let p = ray
                .origin
                .offset(ray.dir.scale((0.12 - ray.origin.z) / ray.dir.z));
            p.x.hypot(p.y)
        });
        radii.sort_by(f64::total_cmp);
        assert!(radii[2] - radii[1] > 1e-8);
        let physical = camera(glass(), 550.0, Some((radii[1] + radii[2]) * 0.5));
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        let clock = Dd::from_f64(2.0_f64.powi(30)) + Dd::from_f64(1e-8);
        let actual = raster
            .sample_pixel_sequential_transient_at_time(
                cx,
                &timed,
                0,
                clock,
                settings(),
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        let steady = raster
            .sample_pixel_sequential_analyzed(cx, &scene, 0, analyzer, settings(), SEED, 0, SAMPLE)
            .unwrap();
        assert_eq!(*actual.pixel_sample(), steady);
        assert_eq!(actual.observation_time_s(), clock);
        assert_eq!(steady.camera_rays().iter().flatten().count(), 2);
        let queries = history.queries.lock().unwrap().clone();
        assert_eq!(queries.len(), 2);
        let beams = packet(&physical, cx, steady.wavelengths_nm(), Some(analyzer));
        for (lane, beam) in beams.iter().enumerate() {
            let Some(beam) = beam else {
                continue;
            };
            let query = queries
                .iter()
                .find(|query| query.wavelength_nm == beam.wavelength_nm())
                .unwrap();
            let group = steady.spectral_samples()[lane]
                .unwrap()
                .path()
                .group_optical_path_length_m()
                .unwrap();
            assert_eq!(query.observation_time_s, clock);
            assert_eq!(
                query.group_delay_s,
                beam.group_delay_s().unwrap() + Dd::from_f64(group) / Dd::from_f64(C)
            );
            assert_eq!(query.emission_time_s, clock - query.group_delay_s);
        }
        assert_eq!(
            actual,
            raster
                .sample_pixel_sequential_transient_at_time(
                    cx,
                    &timed,
                    0,
                    clock,
                    settings(),
                    SEED,
                    0,
                    SAMPLE
                )
                .unwrap()
        );
        let closed_camera = camera(glass(), 550.0, Some(radii[0] * 0.5));
        let closed = CurvedCameraRaster::try_new(&closed_camera, 1, 1).unwrap();
        let before = history.queries.lock().unwrap().len();
        let dark = closed
            .sample_pixel_sequential_transient_at_time(
                cx,
                &timed,
                0,
                clock,
                settings(),
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        assert_eq!(dark.pixel_sample().xyz(), [0.0; 3]);
        assert_eq!(dark.pixel_sample().camera_rays(), &[None; 4]);
        assert_eq!(dark.observation_time_s(), clock);
        assert_eq!(history.queries.lock().unwrap().len(), before);
        assert_eq!(
            closed
                .sample_pixel_sequential_transient_at_time(
                    cx,
                    &timed,
                    0,
                    Dd::from_f64(f64::NAN),
                    settings(),
                    SEED,
                    0,
                    SAMPLE
                )
                .unwrap_err(),
            TransientRasterError::Source(TransientCameraError::InvalidObservationTime)
        );
        gate.request();
        assert!(
            closed
                .sample_pixel_sequential_transient_at_time(
                    cx,
                    &timed,
                    0,
                    clock,
                    settings(),
                    SEED,
                    0,
                    SAMPLE
                )
                .is_err()
        );
    });
}

#[test]
fn infrared_sensor_integrates_retarded_lens_beer_radiance_and_preserves_source_failures() {
    with_cx(|cx, _| {
        // Exact binary index steps give one differentiable affine table law;
        // the group law is qualified even at its interior knots.
        let ior = TabulatedIor::try_new_in_band(
            core::array::from_fn(|i| {
                DielectricIorSample::try_new_in_band(
                    8000.0 + 750.0 * i as f64,
                    1.5625 - i as f64 / 128.0,
                )
                .unwrap()
            }),
            750.0,
        )
        .unwrap();
        let material = DielectricGlass::new(
            ior,
            BeerLambertAbsorption::try_constant(2.0).unwrap(),
            GlassProvenance::Custom,
        );
        let geometry = geometry();
        // An explicit IR table also admits the exterior's source band.
        let vacuum_ior = TabulatedIor::try_new_in_band(
            core::array::from_fn(|i| {
                DielectricIorSample::try_new_in_band(8000.0 + 750.0 * i as f64, 1.0).unwrap()
            }),
            750.0,
        )
        .unwrap();
        let vacuum = DielectricGlass::new(
            vacuum_ior,
            BeerLambertAbsorption::try_constant(0.12).unwrap(),
            GlassProvenance::Custom,
        );
        let media = [SceneMedium::Homogeneous(&vacuum)];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100000).unwrap();
        let sources = [SurfaceEmitter::new(
            7,
            RadianceSpectrum::try_new_in_band(vec![(8000.0, 3.0), (14000.0, 3.0)]).unwrap(),
            EmissionSides::Both,
        )];
        let scene = CurvedCameraScene::try_new(cx, &native, &sources, None).unwrap();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let clock = epoch + Dd::from_f64(7e-9);
        let pulse = History::new(|q: SurfaceEmissionQuery| {
            Ok((-((q.emission_time_s - epoch).to_f64() / 1e-9).powi(2)).exp())
        });
        let timed = TransientCurvedCameraScene::try_new(cx, &scene, &[binding(&pulse)]).unwrap();
        let physical = camera(material, 11000.0, None);
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        let mut trace_settings = settings();
        trace_settings.event.surface_tolerance_m = 1e-9;
        let sensor = QuantumEfficiencySensor::try_new(vec![
            (8000.0, [0.0, 1.0, 0.0]),
            (14000.0, [1.0, 0.0, 0.0]),
        ])
        .unwrap();
        let sampled = raster
            .sample_pixel_sensor_sequential_transient_at_time(
                cx,
                &timed,
                &sensor,
                0,
                clock,
                trace_settings,
                SEED,
                0,
                SAMPLE,
            )
            .unwrap();
        let output = sampled.pixel_sample();
        let wavelengths = output
            .detector_wavelengths()
            .map(|draw| draw.wavelength_nm());
        let beams = packet(&physical, cx, &wavelengths, None);
        let mut energy = [0.0; 3];
        let mut electrons = [0.0; 3];
        for (i, beam) in beams.iter().enumerate() {
            let beam = beam.unwrap();
            let ray = beam.ray();
            let length = (2.0 - ray.origin.z) / ray.dir.z;
            let delay = beam.group_delay_s().unwrap() + Dd::from_f64(length) / Dd::from_f64(C);
            let source_time = ((clock - epoch) - delay).to_f64();
            let radiance =
                3.0 * beam.power().i * (-0.12 * length - (source_time / 1e-9).powi(2)).exp();
            near(
                output.spectral_samples()[i].unwrap().radiance_w_m2_sr_nm(),
                radiance,
                2e-7,
            );
            let nm = beam.wavelength_nm();
            let qe = [(nm - 8000.0) / 6000.0, (14000.0 - nm) / 6000.0, 0.0];
            for c in 0..3 {
                let weight = radiance * beam.geometry_weight() * qe[c] * 6000.0 / 4.0;
                energy[c] += weight;
                electrons[c] += weight * nm * 1e-9 / (6.62607015e-34 * C);
            }
        }
        for c in 0..3 {
            near(
                output.qe_weighted_radiance_w_m2_sr()[c],
                energy[c],
                2e-7 * energy[c].max(1.0),
            );
            near(
                output.photoelectron_rate_per_m2_sr()[c],
                electrons[c],
                2e-7 * electrons[c].max(1.0),
            );
        }
        assert_eq!(pulse.queries.lock().unwrap().len(), 4);
        assert!(
            pulse
                .queries
                .lock()
                .unwrap()
                .iter()
                .all(|q| q.observation_time_s == clock)
        );
        assert_eq!(output.photoelectron_rate_per_m2_sr()[2], 0.0);
        let failure = History::new(|_| Err(SurfaceHistoryError::OutsideTimeDomain));
        let failed = TransientCurvedCameraScene::try_new(cx, &scene, &[binding(&failure)]).unwrap();
        assert_eq!(
            raster
                .sample_pixel_sensor_sequential_transient_at_time(
                    cx,
                    &failed,
                    &sensor,
                    0,
                    clock,
                    trace_settings,
                    SEED,
                    0,
                    SAMPLE
                )
                .unwrap_err(),
            TransientRasterError::Source(TransientCameraError::History {
                object_id: 7,
                error: SurfaceHistoryError::OutsideTimeDomain
            })
        );
        assert_eq!(failure.queries.lock().unwrap().len(), 1);
        assert!(matches!(
            raster.sample_pixel_sequential_transient_at_time(
                cx,
                &timed,
                0,
                clock,
                trace_settings,
                SEED,
                0,
                SAMPLE
            ),
            Err(TransientRasterError::Raster(CurvedRasterError::Camera(
                CameraError::SequentialLensOutsideModel
            )))
        ));
    });
}
