use super::super::tests::{close, detector, glass, request, source, square, with_cx};
use super::*;
use crate::dielectric::DielectricSurface;
use crate::instances::{GeometryInstance, InstanceScene, RigidTransform, SharedGeometry};
use crate::photometry::source::{AngularDistribution, SPECTRAL_BINS};
use crate::polarization::dielectric::scene::curved::participating::CurvedRayleighMedium;
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, SceneMedium};
use crate::polarization::dielectric::scene::{Interaction, PathTermination, SurfaceBinding};
use crate::polarization::photometry::SpectralOptics;
use crate::volumes::transport::scattering::{MassScatteringSpectrum, PhaseFunction};
use fs_blake3::hash_domain;
use fs_geom::Point3;
use fs_geom::fixtures::SphereChart;

fn sphere(id: u64, radius: f64) -> GeometryInstance {
    GeometryInstance::try_new(
        id,
        hash_domain("curved-scattering-detector-sphere", &id.to_le_bytes()),
        SharedGeometry::chart(SphereChart {
            center: Point3::new(0.0, 0.0, -1.0),
            radius,
        }),
        RigidTransform::identity(),
    )
    .unwrap()
}

fn lens_binding() -> SurfaceBinding {
    SurfaceBinding {
        object_id: 1,
        interaction: Interaction::Dielectric {
            outside: 0,
            inside: 1,
            surface: DielectricSurface::SMOOTH,
        },
    }
}

fn options(samples: u32) -> RayleighDetectorRequest {
    let mut detector = request(samples);
    detector.initial_medium = 1;
    detector.trace.max_surface_events = 128;
    detector.trace.max_distance_m = 128.0;
    detector.trace.max_total_steps = 65_536;
    RayleighDetectorRequest {
        detector,
        max_volume_events: 128,
    }
}

fn raster(z: f64, width: f64) -> DetectorRaster {
    DetectorRaster::try_new(
        RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], [0.0, 0.0, z]).unwrap(),
        width,
        width,
        3,
        3,
    )
    .unwrap()
}

#[test]
fn g3_zero_scattering_preserves_complete_detector_and_polarized_sensor_measurements() {
    with_cx(|cx, gate| {
        let source = source(AngularDistribution::uniform_cone(0.12).unwrap());
        let filtered = source
            .clone()
            .with_unpolarized_optics()
            .through(
                &SpectralOptics::neutral(core::array::from_fn(|bin| {
                    if bin % 2 == 0 { 0.2 } else { 0.7 }
                }))
                .unwrap(),
            )
            .unwrap()
            .through(&SpectralOptics::linear_polarizer(0.3).unwrap())
            .unwrap();
        let geometry = InstanceScene::try_new(vec![square(9, 0.0, 2.0, true)]).unwrap();
        let carrier = glass(1.5, 0.3);
        let media = [SceneMedium::Homogeneous(&carrier)];
        let surfaces =
            CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector()], 100_000)
                .unwrap();
        let coefficient = MassScatteringSpectrum::rayleigh(550.0, 1.0).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [Some(
            CurvedRayleighMedium::try_new(&coefficient, 0.0, &phase).unwrap(),
        )];
        let scene = CurvedRayleighScene::try_new(cx, &surfaces, &volumes).unwrap();
        let request = RayleighDetectorRequest {
            detector: request(8),
            max_volume_events: 4,
        };
        assert_eq!(
            source
                .measure_rayleigh_detector(cx, &scene, request)
                .unwrap(),
            source
                .measure_detector(cx, &surfaces, request.detector)
                .unwrap()
        );
        assert_eq!(
            filtered
                .measure_rayleigh_detector(cx, &scene, request)
                .unwrap(),
            filtered
                .measure_detector(cx, &surfaces, request.detector)
                .unwrap()
        );
        assert_eq!(
            source
                .render_rayleigh_detector(cx, &scene, request, raster(0.0, 2.0))
                .unwrap(),
            source
                .render_detector(cx, &surfaces, request.detector, raster(0.0, 2.0))
                .unwrap()
        );
        assert_eq!(
            filtered
                .render_rayleigh_detector(cx, &scene, request, raster(0.0, 2.0))
                .unwrap(),
            filtered
                .render_detector(cx, &surfaces, request.detector, raster(0.0, 2.0))
                .unwrap()
        );
        assert_eq!(
            source
                .measure_rayleigh_detector(
                    cx,
                    &scene,
                    RayleighDetectorRequest {
                        max_volume_events: 0,
                        ..request
                    }
                )
                .unwrap_err(),
            FiniteDetectorError::InvalidRequest
        );
        gate.request();
        assert_eq!(
            source
                .measure_rayleigh_detector(cx, &scene, request)
                .unwrap_err(),
            FiniteDetectorError::Cancelled
        );
    });
}

#[test]
fn g1_lossless_scattering_glass_conserves_calibrated_spectral_power_in_expectation() {
    with_cx(|cx, _| {
        // The actual sphere detector captures every escaped direction. No
        // aperture factor or lens magnification belongs in this power balance.
        let source = source(AngularDistribution::isotropic());
        let geometry = InstanceScene::try_new(vec![sphere(1, 0.5), sphere(9, 2.0)]).unwrap();
        let air = glass(1.0, 0.0);
        let lens = glass(1.5, 0.0);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&lens),
        ];
        let surfaces = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![lens_binding(), detector()],
            100_000,
        )
        .unwrap();
        let coefficient = MassScatteringSpectrum::rayleigh(550.0, 0.7).unwrap();
        let phase = PhaseFunction::Rayleigh;
        let volumes = [
            None,
            Some(CurvedRayleighMedium::try_new(&coefficient, 1.0, &phase).unwrap()),
        ];
        let scene = CurvedRayleighScene::try_new(cx, &surfaces, &volumes).unwrap();
        let result = source
            .measure_rayleigh_detector(cx, &scene, options(96))
            .unwrap();
        assert_eq!(result.detected_paths(), result.launched_paths());
        let mut emitted_w = 0.0;
        for bin in 0..SPECTRAL_BINS {
            let emitted = source.flux_lm() * source.spectrum().bins_w_per_lumen_nm()[bin];
            close(
                result.spectral_power_w_nm()[bin].i,
                emitted,
                6.0 * result.spectral_standard_error_w_nm()[bin][0] + 2e-11 * emitted,
            );
            emitted_w += emitted * crate::photometry::source::BIN_WIDTH_NM;
        }
        close(
            result.radiant_stokes_w().i,
            emitted_w,
            6.0 * result.radiant_standard_error_w()[0] + 2e-11 * emitted_w,
        );
        close(
            result.photopic_stokes_lm().i,
            source.flux_lm(),
            6.0 * result.photopic_standard_error_lm()[0] + 2e-11 * source.flux_lm(),
        );
        // Rayleigh M/PDF may exceed one for a polarized packet. The changing
        // random power confirms actual polarized multiple scattering occurred;
        // clamping weights or applying scattering Beer again violates balance.
        assert!(
            result
                .spectral_standard_error_w_nm()
                .iter()
                .any(|error| error[0] > 1e-8)
        );
    });
}

#[test]
fn g1_scattering_in_glass_forms_a_sensor_halo_and_refuses_surviving_volume_caps() {
    with_cx(|cx, _| {
        let source = source(AngularDistribution::uniform_cone(0.02).unwrap());
        let geometry =
            InstanceScene::try_new(vec![sphere(1, 0.5), square(9, 2.0, 2.0, true)]).unwrap();
        let air = glass(1.0, 0.0);
        let lens = glass(1.5, 0.0);
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&lens),
        ];
        let surfaces = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![lens_binding(), detector()],
            100_000,
        )
        .unwrap();
        let coefficient = MassScatteringSpectrum::constant(1.5).unwrap();
        let phase = PhaseFunction::rayleigh_depolarized(0.0279).unwrap();
        let volumes = [
            None,
            Some(CurvedRayleighMedium::try_new(&coefficient, 1.0, &phase).unwrap()),
        ];
        let scene = CurvedRayleighScene::try_new(cx, &surfaces, &volumes).unwrap();
        let request = options(64);
        let layout = raster(2.0, 4.0);
        let clear = source
            .render_detector(cx, &surfaces, request.detector, layout)
            .unwrap();
        assert!(
            clear
                .pixels()
                .iter()
                .enumerate()
                .all(|(pixel, value)| pixel == 4 || value.radiant_stokes_w().i == 0.0)
        );
        let image = source
            .render_rayleigh_detector(cx, &scene, request, layout)
            .unwrap();
        assert!(
            image
                .pixels()
                .iter()
                .enumerate()
                .any(|(pixel, value)| pixel != 4 && value.radiant_stokes_w().i > 0.0)
        );
        assert!(image.pixels()[4].radiant_stokes_w().i < clear.pixels()[4].radiant_stokes_w().i);
        let total = source
            .measure_rayleigh_detector(cx, &scene, request)
            .unwrap();
        assert_eq!(image.total_detector(), &total);
        let pixels_w: f64 = image
            .pixels()
            .iter()
            .map(|pixel| pixel.radiant_stokes_w().i)
            .sum();
        close(
            pixels_w + image.cropped().radiant_stokes_w().i,
            total.radiant_stokes_w().i,
            2e-12 * total.radiant_stokes_w().i,
        );
        let pixels_lm: f64 = image
            .pixels()
            .iter()
            .map(|pixel| pixel.photopic_stokes_lm().i)
            .sum();
        close(
            pixels_lm + image.cropped().photopic_stokes_lm().i,
            total.photopic_stokes_lm().i,
            2e-12 * total.photopic_stokes_lm().i,
        );
        // All 80 bins share nondispersive geometry, scattering and random draws.
        // The broadband standard error must retain that perfect correlation.
        let first = total.spectral_power_w_nm()[0].i;
        close(
            total.photopic_standard_error_lm()[0] / total.photopic_stokes_lm().i,
            total.spectral_standard_error_w_nm()[0][0] / first,
            2e-12,
        );
        assert!(matches!(
            source.measure_rayleigh_detector(
                cx,
                &scene,
                RayleighDetectorRequest {
                    max_volume_events: 1,
                    ..request
                }
            ),
            Err(FiniteDetectorError::IncompletePath {
                termination: PathTermination::VolumeLimit,
                ..
            })
        ));
    });
}
