//! Retarded point radiometry through the original CIE and calibrated QE pixels.

use super::*;
use crate::polarization::dielectric::scene::curved::camera::raster::transient::point_lights::TransientPointLightRasterError;
use crate::polarization::dielectric::scene::curved::camera::raster::{LAMBDA_MAX, LAMBDA_MIN};
use crate::polarization::dielectric::scene::curved::camera::transient::exposure::ExposureWindow;
use crate::polarization::dielectric::scene::curved::camera::transient::point_lights::{
    PointEmissionHistory, PointEmissionHistoryBinding, PointEmissionQuery, PointExposureQuery,
    TransientCurvedPointLightCameraScene, TransientPointLightCameraError,
};
use crate::polarization::dielectric::scene::curved::camera::transient::{
    SurfaceHistoryError, TransientCameraError, TransientCurvedCameraScene,
};
use crate::spectral::detector::QuantumEfficiencySensor;
use fs_blake3::ContentHash;
use fs_geom::Point3;
use fs_math::dd::Dd;
use std::sync::atomic::{AtomicUsize, Ordering};

pub(in super::super) mod sequential;

pub(in super::super) const LIGHT_SPEED: f64 = 299_792_458.0;
pub(in super::super) const ID: ContentHash = ContentHash([127; 32]);

pub(in super::super) struct LinearHistory<'a> {
    source: &'a PhotometricSource,
    epoch: Dd,
    offset: f64,
    slope: f64,
    calls: AtomicUsize,
    fail_at: AtomicUsize,
}
impl<'a> LinearHistory<'a> {
    pub(in super::super) fn new(source: &'a PhotometricSource, epoch: Dd, slope: f64) -> Self {
        Self {
            source,
            epoch,
            offset: 1.0,
            slope,
            calls: AtomicUsize::new(0),
            fail_at: AtomicUsize::new(0),
        }
    }
    fn evaluate(
        &self,
        cx: &Cx<'_>,
        query: PointEmissionQuery<'_>,
        time: Dd,
    ) -> Result<f64, SurfaceHistoryError> {
        cx.checkpoint()?;
        assert!(core::ptr::eq(query.source, self.source));
        assert_eq!(query.source_identity, ID);
        assert_eq!(query.position_m, self.source.position_m());
        assert_eq!(
            query.emission_time_s,
            query.observation_time_s - query.group_delay_s
        );
        let count = self.calls.fetch_add(1, Ordering::Relaxed) + 1;
        if count == self.fail_at.load(Ordering::Relaxed) {
            return Err(SurfaceHistoryError::OutsideTimeDomain);
        }
        Ok(self.offset + self.slope * (time - self.epoch).to_f64())
    }
}
impl PointEmissionHistory for LinearHistory<'_> {
    fn sample(
        &self,
        cx: &Cx<'_>,
        query: PointEmissionQuery<'_>,
    ) -> Result<f64, SurfaceHistoryError> {
        self.evaluate(cx, query, query.emission_time_s)
    }
    fn average(
        &self,
        cx: &Cx<'_>,
        query: PointExposureQuery<'_>,
    ) -> Result<f64, SurfaceHistoryError> {
        // Exact integral of this affine fixture history over its own interval.
        self.evaluate(
            cx,
            query.opening,
            query.opening.emission_time_s + query.duration_s * Dd::from_f64(0.5),
        )
    }
}
pub(in super::super) fn history_binding<'a>(
    source: &'a PhotometricSource,
    history: &'a dyn PointEmissionHistory,
) -> PointEmissionHistoryBinding<'a> {
    PointEmissionHistoryBinding {
        source,
        source_identity: ID,
        history,
    }
}
pub(in super::super) fn fixture_source() -> PhotometricSource {
    PhotometricSource::try_new(
        RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.4, 0.0, 0.0]).unwrap(),
        1000.0,
        1.0,
        AngularDistribution::isotropic(),
        SourceSpectrum::try_from_bins(core::array::from_fn(|i| 1.0 + i as f64 / 80.0)).unwrap(),
    )
    .unwrap()
}
pub(in super::super) fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 8e-6 * expected.abs() + 1e-14,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_point_pixels_use_each_chromatic_clock_and_original_cie_qe_quadrature() {
    with_cx(|cx, gate| {
        let geometry = geometry();
        let material = DielectricGlass::new(
            CauchyIor::try_new(1.4, 0.006, 0.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.2).unwrap(),
            GlassProvenance::Custom,
        );
        let media = [SceneMedium::Homogeneous(&material)];
        let reflectances = [SpectralReflectance::constant(0.37).unwrap()];
        let native = CurvedPolarizedScene::try_new_with_materials(
            cx,
            &geometry,
            &media,
            &[],
            &[],
            &reflectances,
            bindings(),
            100_000,
        )
        .unwrap();
        let surface = CurvedCameraScene::try_new(cx, &native, &[], None).unwrap();
        let source = fixture_source();
        let branches = [PointLightBranch {
            source: &source,
            receiver_object_id: 7,
            initial_medium: 0,
            initial_direction: Vec3::new(0.0, 0.0, 1.0),
            source_transverse_axis: Vec3::new(0.0, 1.0, 0.0),
            policy: BundleBranchPolicy::TransmitOrReflectTir,
            position_tolerance_m: 1e-6,
        }];
        let lights = CurvedPointLightCameraScene::try_new(cx, &surface, &branches).unwrap();
        let epoch = Dd::from_f64(2.0_f64.powi(30));
        let clock = epoch + Dd::from_f64(30e-9);
        let timed = TransientCurvedCameraScene::try_new(cx, &surface, &[]).unwrap();
        let history = LinearHistory::new(&source, epoch, 1e7);
        let histories = [history_binding(&source, &history)];
        let observer =
            TransientCurvedPointLightCameraScene::try_new(cx, &timed, &lights, &histories).unwrap();
        let unit = LinearHistory::new(&source, epoch, 0.0);
        let unit_histories = [history_binding(&source, &unit)];
        let steady_observer =
            TransientCurvedPointLightCameraScene::try_new(cx, &timed, &lights, &unit_histories)
                .unwrap();
        let physical = camera(0.0);
        let raster = CurvedCameraRaster::try_new(&physical, 3, 2).unwrap();
        let (seed, pixel_id, sample_id) = (0x1234_5678_0000_0049_u64, 4, 9);
        let ordinary = raster
            .sample_pixel_point_lights(cx, &lights, 0, point_settings(), seed, pixel_id, sample_id)
            .unwrap();
        let steady = raster
            .sample_pixel_point_lights_at_time(
                cx,
                &steady_observer,
                0,
                clock,
                point_settings(),
                seed,
                pixel_id,
                sample_id,
            )
            .unwrap();
        assert_eq!(steady.point_light_sample(), &ordinary);
        let result = raster
            .sample_pixel_point_lights_at_time(
                cx,
                &observer,
                0,
                clock,
                point_settings(),
                seed,
                pixel_id,
                sample_id,
            )
            .unwrap();
        assert_eq!(result.observation_time_s(), clock);
        let pixel = result.pixel_sample();
        let ray = pixel.camera_ray();
        let camera_length = 2.0 / ray.dir.z;
        let receiver = ray.origin.offset(ray.dir.scale(camera_length));
        let delta = receiver.delta_from(source.position_m());
        let distance_squared = delta.dot(delta);
        let distance = distance_squared.sqrt();
        let direction = delta.scale(1.0 / distance);
        let expected = |wavelength: f64, observation: Dd| {
            let group_index = 1.4 + 0.018 / (wavelength * 1e-3).powi(2);
            let delay = group_index * (camera_length + distance) / LIGHT_SPEED;
            let scale = 1.0 + 1e7 * ((observation - epoch).to_f64() - delay);
            0.37 / core::f64::consts::PI
                * source
                    .spectral_intensity_w_sr_nm(direction, wavelength)
                    .unwrap()
                * direction.z
                / distance_squared
                * (-0.2 * (camera_length + distance)).exp()
                * scale
        };
        let mut xyz = [0.0; 3];
        for (lane, &wavelength) in pixel.wavelengths_nm().iter().enumerate() {
            let radiance = expected(wavelength, clock);
            near(
                pixel.spectral_samples()[lane].radiance_w_m2_sr_nm(),
                radiance,
            );
            assert_eq!(
                pixel.spectral_samples()[lane].path(),
                ordinary.pixel_sample().spectral_samples()[lane].path()
            );
            assert_eq!(result.point_light_sample().connections()[lane], 1);
            for (sum, response) in
                xyz.iter_mut()
                    .zip([cie_x(wavelength), cie_y(wavelength), cie_z(wavelength)])
            {
                *sum += radiance * response * (LAMBDA_MAX - LAMBDA_MIN) / 4.0 / y_integral();
            }
        }
        for (actual, expected) in pixel.xyz().into_iter().zip(xyz) {
            near(actual, expected);
        }
        assert_eq!(
            result.point_light_sample().trace_evaluations(),
            ordinary.trace_evaluations()
        );
        assert_eq!(
            result.point_light_sample().field_evaluations(),
            ordinary.field_evaluations()
        );
        assert_eq!(
            result.point_light_sample().query_work(),
            ordinary.query_work()
        );

        let one = QuantumEfficiencySensor::try_new(vec![(LAMBDA_MIN, [1.0]), (LAMBDA_MAX, [1.0])])
            .unwrap();
        let three = QuantumEfficiencySensor::try_new(vec![
            (LAMBDA_MIN, [0.25, 0.75, 0.0]),
            (LAMBDA_MAX, [0.25, 0.75, 0.0]),
        ])
        .unwrap();
        let a = raster
            .sample_pixel_sensor_point_lights_at_time(
                cx,
                &observer,
                &one,
                0,
                clock,
                point_settings(),
                seed,
                pixel_id,
                sample_id,
            )
            .unwrap();
        let b = raster
            .sample_pixel_sensor_point_lights_at_time(
                cx,
                &observer,
                &three,
                0,
                clock,
                point_settings(),
                seed,
                pixel_id,
                sample_id,
            )
            .unwrap();
        assert_eq!(
            a.pixel_sample().camera_rays(),
            b.pixel_sample().camera_rays()
        );
        assert_eq!(
            a.pixel_sample().spectral_samples(),
            b.pixel_sample().spectral_samples()
        );
        assert_eq!(
            a.point_light_sample().trace_evaluations(),
            b.point_light_sample().trace_evaluations()
        );
        let window = ExposureWindow::try_new(clock, clock + Dd::from_f64(4e-9)).unwrap();
        let exposed = raster
            .sample_pixel_sensor_point_lights_during_exposure(
                cx,
                &observer,
                &three,
                0,
                window,
                point_settings(),
                seed,
                pixel_id,
                sample_id,
            )
            .unwrap();
        let color = raster
            .sample_pixel_point_lights_during_exposure(
                cx,
                &observer,
                0,
                window,
                point_settings(),
                seed,
                pixel_id,
                sample_id,
            )
            .unwrap();
        let middle = raster
            .sample_pixel_point_lights_at_time(
                cx,
                &observer,
                0,
                clock + Dd::from_f64(2e-9),
                point_settings(),
                seed,
                pixel_id,
                sample_id,
            )
            .unwrap();
        for channel in 0..3 {
            near(
                color.pixel_sample().mean_pixel_sample().xyz()[channel],
                middle.pixel_sample().xyz()[channel],
            );
        }
        let mut energy = [0.0; 3];
        let mut photons = [0.0; 3];
        for (lane, draw) in b.pixel_sample().detector_wavelengths().iter().enumerate() {
            let wavelength = draw.wavelength_nm();
            let mean = expected(wavelength, clock + Dd::from_f64(2e-9));
            near(
                b.pixel_sample().spectral_samples()[lane].radiance_w_m2_sr_nm(),
                expected(wavelength, clock),
            );
            for channel in 0..3 {
                let contribution =
                    mean * [0.25, 0.75, 0.0][channel] * (LAMBDA_MAX - LAMBDA_MIN) / 4.0 * 4e-9;
                energy[channel] += contribution;
                photons[channel] +=
                    contribution * wavelength * 1e-9 / (6.62607015e-34 * LIGHT_SPEED);
            }
        }
        for channel in 0..3 {
            near(
                exposed.pixel_sample().qe_weighted_energy_j_m2_sr()[channel],
                energy[channel],
            );
            near(
                exposed.pixel_sample().photoelectrons_per_m2_sr()[channel],
                photons[channel],
            );
        }
        assert_eq!(
            exposed.trace_evaluations(),
            b.point_light_sample().trace_evaluations()
        );
        history.calls.store(0, Ordering::Relaxed);
        history.fail_at.store(3, Ordering::Relaxed);
        assert!(matches!(
            raster.sample_pixel_sensor_point_lights_at_time(
                cx,
                &observer,
                &three,
                0,
                clock,
                point_settings(),
                seed,
                pixel_id,
                sample_id
            ),
            Err(TransientPointLightRasterError::Transport {
                lane: 2,
                error: TransientPointLightCameraError::History {
                    source_identity: ID,
                    error: SurfaceHistoryError::OutsideTimeDomain
                }
            })
        ));
        history.fail_at.store(0, Ordering::Relaxed);
        let replay = raster
            .sample_pixel_sensor_point_lights_at_time(
                cx,
                &observer,
                &three,
                0,
                clock,
                point_settings(),
                seed,
                pixel_id,
                sample_id,
            )
            .unwrap();
        assert_eq!(replay, b);
        gate.request();
        assert!(
            raster
                .sample_pixel_point_lights_during_exposure(
                    cx,
                    &observer,
                    0,
                    window,
                    point_settings(),
                    seed,
                    pixel_id,
                    sample_id
                )
                .is_err()
        );
    });
}
