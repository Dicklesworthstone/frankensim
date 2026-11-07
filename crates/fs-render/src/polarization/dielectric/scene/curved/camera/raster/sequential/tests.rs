use super::super::tests::{bindings, geometry, settings, with_cx};
use super::*;
use crate::camera::{Aperture, CameraProjection, LensSurface, SequentialLens};
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricIorSample, GlassProvenance,
    TabulatedIor,
};
use crate::gradient_index::{GradientIndexError, GradientIndexField, GradientIndexSample};
use crate::polarization::dielectric::scene::curved::camera::{
    EmissionSides, RadianceSpectrum, SurfaceEmitter,
};
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, SceneMedium};
use crate::spectral::detector::QuantumEfficiencySensor;
use fs_blake3::ContentHash;
use fs_geom::{Point3, Vec3};
use std::sync::atomic::{AtomicUsize, Ordering};

#[cfg(feature = "volumes")]
mod participating;

const SEED: u64 = 79;
const SAMPLE: u32 = 6;

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

fn camera(glass: DielectricGlass, reference_nm: f64, stop: Option<f64>) -> PhysicalCamera {
    let mut caps = vec![
        LensSurface::try_new(0.08, 20.0, 0.01, Some(glass)).unwrap(),
        LensSurface::try_new(0.09, -20.0, 0.012, None).unwrap(),
        // A wide vacuum stop puts every GRIN launch on the same n=1 plane.
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
    .with_sequential_lens(SequentialLens::try_new(reference_nm, &caps).unwrap())
    .unwrap()
}

fn visible_glass() -> DielectricGlass {
    DielectricGlass::new(
        CauchyIor::try_new(1.48, 0.006, 0.0).unwrap(),
        BeerLambertAbsorption::try_constant(2.0).unwrap(),
        GlassProvenance::Custom,
    )
}

fn optical_packet(
    camera: &PhysicalCamera,
    cx: &Cx<'_>,
    wavelengths: &[f64; CAMERA_WAVELENGTHS],
) -> [SequentialCameraBeam; CAMERA_WAVELENGTHS] {
    let key = [SEED as u32, (SEED >> 32) as u32];
    let raster = philox4x32_10([0, SAMPLE, 0, 0x6372_7331], key);
    let lens = philox4x32_10([0, SAMPLE, 0, 0x6372_7332], key);
    let half = camera.raster_vertical_half_tan();
    camera
        .generate_sequential_spectral_beams_from_tangent_offsets(
            cx,
            (2.0 * unit(raster[0]) - 1.0) * half,
            (1.0 - 2.0 * unit(raster[1])) * half,
            LensSample::try_new(unit(lens[0]), unit(lens[1])).unwrap(),
            wavelengths,
        )
        .unwrap()
        .map(Option::unwrap)
}

struct Field {
    curved: bool,
    calls: AtomicUsize,
}
impl GradientIndexField for Field {
    fn sample(
        &self,
        p: Point3,
        wavelength: f64,
    ) -> Result<GradientIndexSample, GradientIndexError> {
        self.calls.fetch_add(1, Ordering::Relaxed);
        let color = (550.0 / wavelength).powi(2);
        let gradient = if self.curved {
            0.02 + 0.01 * color
        } else {
            0.0
        };
        let n = 1.0 + gradient * (p.z - 0.1);
        let group = n + if self.curved {
            0.02 * color * (p.z - 0.1)
        } else {
            0.0
        };
        GradientIndexSample::try_new(n, Vec3::new(0.0, 0.0, gradient), 0.12, Some(group))
    }
    fn is_spatially_uniform(&self) -> bool {
        !self.curved
    }
}

fn emitter(min: f64, max: f64) -> SurfaceEmitter {
    SurfaceEmitter::new(
        7,
        RadianceSpectrum::try_new_in_band(vec![(min, 3.0), (max, 3.0)]).unwrap(),
        EmissionSides::Both,
    )
}

fn stop_radii(beams: &[SequentialCameraBeam; CAMERA_WAVELENGTHS]) -> [f64; CAMERA_WAVELENGTHS] {
    let mut radii = beams.map(|beam| {
        let ray = beam.ray();
        let point = ray
            .origin
            .offset(ray.dir.scale((0.12 - ray.origin.z) / ray.dir.z));
        point.x.hypot(point.y)
    });
    radii.sort_by(f64::total_cmp);
    assert!(
        radii[2] - radii[1] > 1e-8,
        "fixture must spectrally separate the aperture footprints"
    );
    radii
}

#[test]
fn resolved_lens_continues_each_native_curved_path_and_analyzer_without_double_loss() {
    with_cx(|cx, _| {
        let geometry = geometry();
        let field = Field {
            curved: true,
            calls: AtomicUsize::new(0),
        };
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity: ContentHash([61; 32]),
        }];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let emitters = [emitter(380.0, 780.0)];
        let scene = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let physical = camera(visible_glass(), 550.0, None);
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        let sample = raster
            .sample_pixel_sequential(cx, &scene, 0, settings(), SEED, 0, SAMPLE)
            .unwrap();
        let beams = optical_packet(&physical, cx, sample.wavelengths_nm());
        let mut xyz = [0.0; 3];
        for (i, optical) in beams.iter().enumerate() {
            let wavelength = optical.wavelength_nm();
            let ray = optical.ray();
            let path = sample.spectral_samples()[i].unwrap();
            let n1 = 1.0 + (0.02 + 0.01 * (550.0 / wavelength).powi(2)) * 1.9;
            let p2 = ray.dir.x.powi(2) + ray.dir.y.powi(2);
            // Planar GRIN conserves n*sin(theta). Integrating ds/dz gives
            // this rationalized arc length, independently of its RK solver.
            let length = 1.9 * (n1 + 1.0) / ((n1 * n1 - p2).sqrt() + (1.0 - p2).sqrt());
            let radiance = 3.0 * optical.power().i / n1.powi(2) * (-0.12 * length).exp();
            assert_eq!(sample.camera_rays()[i], Some(ray));
            assert_eq!(sample.geometry_weights()[i], optical.geometry_weight());
            near(path.launch_phase_index(), 1.0, 0.0);
            near(path.radiance_w_m2_sr_nm(), radiance, 2e-7);
            near(path.path().traveled_m(), length, 2e-7);
            near(path.path().beam().ray().dir.x, ray.dir.x / n1, 2e-8);
            near(path.path().beam().ray().dir.y, ray.dir.y / n1, 2e-8);
            for (sum, cie) in
                xyz.iter_mut()
                    .zip([cie_x(wavelength), cie_y(wavelength), cie_z(wavelength)])
            {
                *sum += radiance * optical.geometry_weight() * 100.0 / y_integral() * cie;
            }
        }
        for (actual, expected) in sample.xyz().into_iter().zip(xyz) {
            near(actual, expected, 2e-7);
        }
        assert_ne!(
            sample.camera_rays()[0].unwrap().dir,
            sample.camera_rays()[1].unwrap().dir
        );
        assert_eq!(
            sample,
            raster
                .sample_pixel_sequential(cx, &scene, 0, settings(), SEED, 0, SAMPLE)
                .unwrap()
        );
        assert_eq!(
            raster.sample_pixel(cx, &scene, 0, settings(), SEED, 0, SAMPLE),
            Err(CameraError::SequentialLensRequiresWeightedRays.into())
        );
        let analyze = |angle| {
            raster
                .sample_pixel_sequential_analyzed(
                    cx,
                    &scene,
                    0,
                    CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), angle).unwrap(),
                    settings(),
                    SEED,
                    0,
                    SAMPLE,
                )
                .unwrap()
        };
        let horizontal = analyze(0.0);
        let vertical = analyze(core::f64::consts::FRAC_PI_2);
        for c in 0..3 {
            near(
                horizontal.xyz()[c] + vertical.xyz()[c],
                sample.xyz()[c],
                5e-13,
            );
        }
        assert!((horizontal.xyz()[1] - vertical.xyz()[1]).abs() > 1e-7);
    });
}

#[test]
fn wavelength_aperture_clipping_keeps_original_measure_and_never_launches_missing_rays() {
    with_cx(|cx, gate| {
        let geometry = geometry();
        let field = Field {
            curved: false,
            calls: AtomicUsize::new(0),
        };
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity: ContentHash([62; 32]),
        }];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let emitters = [emitter(380.0, 780.0)];
        let scene = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let wide = camera(visible_glass(), 550.0, None);
        let full = CurvedCameraRaster::try_new(&wide, 1, 1)
            .unwrap()
            .sample_pixel_sequential(cx, &scene, 0, settings(), SEED, 0, SAMPLE)
            .unwrap();
        let beams = optical_packet(&wide, cx, full.wavelengths_nm());
        let radii = stop_radii(&beams);
        let clipped_camera = camera(visible_glass(), 550.0, Some((radii[1] + radii[2]) * 0.5));
        let clipped = CurvedCameraRaster::try_new(&clipped_camera, 1, 1)
            .unwrap()
            .sample_pixel_sequential(cx, &scene, 0, settings(), SEED, 0, SAMPLE)
            .unwrap();
        assert_eq!(clipped.spectral_samples().iter().flatten().count(), 2);
        let mut xyz = [0.0; 3];
        for (i, optical) in beams.iter().enumerate() {
            assert_eq!(
                clipped.camera_rays()[i].is_some(),
                clipped.spectral_samples()[i].is_some()
            );
            if let Some(path) = clipped.spectral_samples()[i] {
                let ray = clipped.camera_rays()[i].unwrap();
                let length = (2.0 - ray.origin.z) / ray.dir.z;
                let radiance = 3.0 * optical.power().i * (-0.12 * length).exp();
                near(path.radiance_w_m2_sr_nm(), radiance, 2e-7);
                let nm = optical.wavelength_nm();
                for (sum, cie) in xyz.iter_mut().zip([cie_x(nm), cie_y(nm), cie_z(nm)]) {
                    *sum += radiance * optical.geometry_weight() * 100.0 / y_integral() * cie;
                }
            } else {
                assert_eq!(clipped.geometry_weights()[i], 0.0);
            }
        }
        for (actual, expected) in clipped.xyz().into_iter().zip(xyz) {
            near(actual, expected, 2e-7);
        }
        assert_eq!(clipped.wavelengths_nm(), full.wavelengths_nm());
        let closed_camera = camera(visible_glass(), 550.0, Some(radii[0] * 0.5));
        let closed = CurvedCameraRaster::try_new(&closed_camera, 1, 1).unwrap();
        let calls_before = field.calls.load(Ordering::Relaxed);
        let black = closed
            .sample_pixel_sequential(cx, &scene, 0, settings(), SEED, 0, SAMPLE)
            .unwrap();
        assert_eq!(black.xyz(), [0.0; 3]);
        assert_eq!(black.camera_rays(), &[None; CAMERA_WAVELENGTHS]);
        assert_eq!(black.spectral_samples(), &[None; CAMERA_WAVELENGTHS]);
        assert_eq!(field.calls.load(Ordering::Relaxed), calls_before);
        let unsupported =
            QuantumEfficiencySensor::try_new(vec![(8000.0, [1.0]), (14000.0, [1.0])]).unwrap();
        assert!(
            closed
                .sample_pixel_sensor_sequential(
                    cx,
                    &scene,
                    &unsupported,
                    0,
                    settings(),
                    SEED,
                    0,
                    SAMPLE
                )
                .is_err()
        );
        gate.request();
        assert!(
            closed
                .sample_pixel_sequential(cx, &scene, 0, settings(), SEED, 0, SAMPLE)
                .is_err()
        );
    });
}

#[test]
fn infrared_qe_retains_clipped_lane_weights_and_sensor_analyzers() {
    with_cx(|cx, _| {
        let ior = TabulatedIor::try_new_in_band(
            core::array::from_fn(|i| {
                let nm = 8000.0 + 750.0 * i as f64;
                DielectricIorSample::try_new_in_band(nm, 1.5 + (11000.0 - nm) * 2e-5).unwrap()
            }),
            750.0,
        )
        .unwrap();
        let glass = DielectricGlass::new(
            ior,
            BeerLambertAbsorption::try_constant(2.0).unwrap(),
            GlassProvenance::Custom,
        );
        let geometry = geometry();
        let field = Field {
            curved: false,
            calls: AtomicUsize::new(0),
        };
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity: ContentHash([63; 32]),
        }];
        let native =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, bindings(), 100_000).unwrap();
        let emitters = [emitter(8000.0, 14000.0)];
        let scene = CurvedCameraScene::try_new(cx, &native, &emitters, None).unwrap();
        let sensor = QuantumEfficiencySensor::try_new(vec![
            (8000.0, [0.0, 1.0, 0.0]),
            (14000.0, [1.0, 0.0, 0.0]),
        ])
        .unwrap();
        let wide = camera(glass, 11000.0, None);
        let full = CurvedCameraRaster::try_new(&wide, 1, 1)
            .unwrap()
            .sample_pixel_sensor_sequential(cx, &scene, &sensor, 0, settings(), SEED, 0, SAMPLE)
            .unwrap();
        let wavelengths = full.detector_wavelengths().map(|draw| draw.wavelength_nm());
        let beams = optical_packet(&wide, cx, &wavelengths);
        let radii = stop_radii(&beams);
        let physical = camera(glass, 11000.0, Some((radii[1] + radii[2]) * 0.5));
        let raster = CurvedCameraRaster::try_new(&physical, 1, 1).unwrap();
        let sample = raster
            .sample_pixel_sensor_sequential(cx, &scene, &sensor, 0, settings(), SEED, 0, SAMPLE)
            .unwrap();
        assert_eq!(sample.spectral_samples().iter().flatten().count(), 2);
        assert_eq!(sample.detector_wavelengths(), full.detector_wavelengths());
        let mut energy = [0.0; 3];
        let mut electrons = [0.0; 3];
        for (i, optical) in beams.iter().enumerate() {
            assert_eq!(
                sample.camera_rays()[i].is_some(),
                sample.spectral_samples()[i].is_some()
            );
            let Some(path) = sample.spectral_samples()[i] else {
                continue;
            };
            let ray = sample.camera_rays()[i].unwrap();
            let nm = optical.wavelength_nm();
            let radiance =
                3.0 * optical.power().i * (-0.12 * (2.0 - ray.origin.z) / ray.dir.z).exp();
            near(path.radiance_w_m2_sr_nm(), radiance, 2e-7);
            let qe = [(nm - 8000.0) / 6000.0, (14000.0 - nm) / 6000.0, 0.0];
            for c in 0..3 {
                // Complementary QE has uniform p(lambda)=1/6000 nm. Keep
                // all four draws in this denominator, including the two misses.
                let weight = radiance * optical.geometry_weight() * qe[c] * 6000.0 / 4.0;
                energy[c] += weight;
                electrons[c] += weight * nm * 1e-9 / (6.626_070_15e-34 * 299_792_458.0);
            }
        }
        for c in 0..3 {
            near(
                sample.qe_weighted_radiance_w_m2_sr()[c],
                energy[c],
                2e-7 * energy[c].max(1.0),
            );
            near(
                sample.photoelectron_rate_per_m2_sr()[c],
                electrons[c],
                2e-7 * electrons[c].max(1.0),
            );
        }
        assert_eq!(sample.photoelectron_rate_per_m2_sr()[2], 0.0);
        assert_eq!(
            sample,
            raster
                .sample_pixel_sensor_sequential(cx, &scene, &sensor, 0, settings(), SEED, 0, SAMPLE)
                .unwrap()
        );
        let analyze = |angle| {
            raster
                .sample_pixel_sensor_sequential_analyzed(
                    cx,
                    &scene,
                    &sensor,
                    0,
                    CameraAnalyzer::try_linear(Vec3::new(1.0, 0.0, 0.0), angle).unwrap(),
                    settings(),
                    SEED,
                    0,
                    SAMPLE,
                )
                .unwrap()
        };
        let a = analyze(0.0);
        let b = analyze(core::f64::consts::FRAC_PI_2);
        for c in 0..3 {
            near(
                a.qe_weighted_radiance_w_m2_sr()[c] + b.qe_weighted_radiance_w_m2_sr()[c],
                sample.qe_weighted_radiance_w_m2_sr()[c],
                2e-9,
            );
            near(
                a.photoelectron_rate_per_m2_sr()[c] + b.photoelectron_rate_per_m2_sr()[c],
                sample.photoelectron_rate_per_m2_sr()[c],
                2e-12 * sample.photoelectron_rate_per_m2_sr()[c].max(1.0),
            );
        }
        assert!(
            raster
                .sample_pixel_sequential(cx, &scene, 0, settings(), SEED, 0, SAMPLE)
                .is_err()
        );
    });
}
