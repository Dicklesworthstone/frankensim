use super::super::tests::{close, detector, glass, request, source, square, with_cx};
use super::super::{FiniteDetectorError, FiniteDetectorRequest};
use crate::dielectric::{
    BeerLambertAbsorption, CauchyIor, DielectricGlass, DielectricSurface, GlassProvenance,
};
use crate::gradient_index::{GradientIndexError, GradientIndexField, GradientIndexSample};
use crate::instances::{InstanceScene, RigidTransform};
use crate::photometry::{
    observer::xyz_at_bin,
    source::{
        AngularDistribution, BIN_WIDTH_NM, LUMINOUS_EFFICACY_LM_W, PhotometricSource,
        SPECTRAL_BINS, SourceSpectrum, wavelength_nm,
    },
};
use crate::polarization::dielectric::scene::curved::{CurvedPolarizedScene, SceneMedium};
use crate::polarization::dielectric::scene::measurement::raster::DetectorRaster;
use crate::polarization::dielectric::scene::{Interaction, PathTermination, SurfaceBinding};
use crate::polarization::photometry::SpectralOptics;
use fs_blake3::hash_domain;
use fs_geom::{Point3, Vec3};

fn raster(center: [f64; 3], width: f64, height: f64, cols: u32, rows: u32) -> DetectorRaster {
    DetectorRaster::try_new(
        RigidTransform::try_new([0.0, 0.0, 0.0, 1.0], center).unwrap(),
        width,
        height,
        cols,
        rows,
    )
    .unwrap()
}

fn add_xyz(sum: &mut [f64; 3], bin: usize, spectral_power: f64) {
    for (channel, weight) in xyz_at_bin(bin).into_iter().enumerate() {
        sum[channel] += spectral_power * BIN_WIDTH_NM * LUMINOUS_EFFICACY_LM_W * weight;
    }
}

#[test]
fn g1_pixels_and_crop_partition_power_with_all_launch_variance() {
    with_cx(|cx, _| {
        let source = source(AngularDistribution::isotropic());
        let geometry = InstanceScene::try_new(vec![square(9, 0.0, 1.0, true)]).unwrap();
        let air = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector()], 100_000)
            .unwrap();
        // A1m² image inside an actual4m² detector makes crop accounting observable.
        let layout = raster([0.0; 3], 1.0, 1.0, 2, 2);
        let options = request(128);
        let image = source.render_detector(cx, &scene, options, layout).unwrap();
        assert_eq!(*image.raster(), layout);
        assert_eq!(image.pixels().len(), 4);
        assert!(image.cropped().detected_paths() > 0);
        let mut emitted_xyz = [0.0; 3];
        let mut emitted_w = 0.0;
        for bin in 0..SPECTRAL_BINS {
            let power = source.flux_lm() * source.spectrum().bins_w_per_lumen_nm()[bin];
            emitted_w += power * BIN_WIDTH_NM;
            add_xyz(&mut emitted_xyz, bin, power);
        }
        let mut power = 0.0;
        let mut xyz = [0.0; 3];
        let mut detected = 0;
        let mut received = 0;
        for pixel in image
            .pixels()
            .iter()
            .chain(core::iter::once(image.cropped()))
        {
            let hits = pixel.detected_paths() / SPECTRAL_BINS as u64;
            assert_eq!(pixel.detected_paths() % SPECTRAL_BINS as u64, 0);
            assert_eq!(u64::from(pixel.received_samples()), hits);
            let fraction = hits as f64 / f64::from(options.samples);
            let sem = (fraction * (1.0 - fraction) / f64::from(options.samples - 1)).sqrt();
            close(
                pixel.radiant_stokes_w().i,
                emitted_w * fraction,
                2e-12 * emitted_w,
            );
            close(
                pixel.radiant_standard_error_w()[0],
                emitted_w * sem,
                2e-12 * emitted_w,
            );
            close(
                pixel.photopic_stokes_lm().i,
                emitted_xyz[1] * fraction,
                2e-12 * emitted_xyz[1],
            );
            close(
                pixel.photopic_standard_error_lm()[0],
                emitted_xyz[1] * sem,
                2e-12 * emitted_xyz[1],
            );
            for channel in 0..3 {
                close(
                    pixel.xyz_flux()[channel],
                    emitted_xyz[channel] * fraction,
                    2e-12 * emitted_xyz[channel],
                );
                close(
                    pixel.xyz_standard_error()[channel],
                    emitted_xyz[channel] * sem,
                    2e-12 * emitted_xyz[channel],
                );
                xyz[channel] += pixel.xyz_flux()[channel];
            }
            power += pixel.radiant_stokes_w().i;
            detected += pixel.detected_paths();
            received += pixel.received_samples();
        }
        close(
            power,
            image.total_detector().radiant_stokes_w().i,
            2e-12 * emitted_w,
        );
        close(
            xyz[1],
            image.total_detector().photopic_stokes_lm().i,
            2e-12 * emitted_xyz[1],
        );
        assert_eq!(detected, image.total_detector().detected_paths());
        assert_eq!(received, image.total_detector().received_samples());
        // Each quadrant subtends atan(ab / [z sqrt(z²+a²+b²)]), with a=b=0.5m.
        let pixel_probability = (0.25 / 1.5_f64.sqrt()).atan() / (4.0 * core::f64::consts::PI);
        let expected_sem =
            (pixel_probability * (1.0 - pixel_probability) / f64::from(options.samples)).sqrt();
        for (index, pixel) in image.pixels().iter().enumerate() {
            let observed = f64::from(pixel.received_samples()) / f64::from(options.samples);
            close(observed, pixel_probability, 5.0 * expected_sem);
            let irradiance = image.pixel_irradiance(index).unwrap();
            close(
                irradiance.photopic_lux.i,
                pixel.photopic_stokes_lm().i / 0.25,
                2e-9,
            );
            close(
                irradiance.photopic_standard_error_lux[0],
                pixel.photopic_standard_error_lm()[0] / 0.25,
                2e-9,
            );
        }
        assert!(image.pixel_irradiance(4).is_err());
    });
}

#[test]
fn g1_exact_grin_focus_concentrates_finite_flux_in_the_central_pixel() {
    // The conserved longitudinal momentum of n(r)=1.5sech(r) gives
    // sinh(r)=tan(theta)sin(z-z_source), so all meridional rays meet atπ.
    struct FocusingField;
    impl GradientIndexField for FocusingField {
        fn sample(&self, p: Point3, _: f64) -> Result<GradientIndexSample, GradientIndexError> {
            let radius = p.x.hypot(p.y);
            let n = 1.5 / radius.cosh();
            let g = if radius == 0.0 {
                0.0
            } else {
                -n * radius.tanh() / radius
            };
            GradientIndexSample::try_new(n, Vec3::new(g * p.x, g * p.y, 0.0), 0.0, Some(n))
        }
    }
    with_cx(|cx, _| {
        let source = source(AngularDistribution::uniform_cone(0.2).unwrap());
        let z = core::f64::consts::PI - 1.0;
        let geometry = InstanceScene::try_new(vec![square(9, z, 0.02, true)]).unwrap();
        let field = FocusingField;
        let media = [SceneMedium::GradientIndex {
            field: &field,
            identity: hash_domain("detector-image-focus", b"n=1.5sech(r),g=1/m"),
        }];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector()], 100_000)
            .unwrap();
        let mut options = request(4);
        options.trace.event.integration.max_step_m = 0.05;
        let image = source
            .render_detector(cx, &scene, options, raster([0.0, 0.0, z], 0.03, 0.03, 3, 3))
            .unwrap();
        for (index, pixel) in image.pixels().iter().enumerate() {
            if index == 4 {
                assert_eq!(pixel.received_samples(), options.samples);
                assert_eq!(
                    pixel.detected_paths(),
                    u64::from(options.samples) * SPECTRAL_BINS as u64
                );
                close(pixel.photopic_stokes_lm().i, source.flux_lm(), 2e-9);
            } else {
                assert_eq!(pixel.received_samples(), 0);
                assert_eq!(pixel.photopic_stokes_lm().i, 0.0);
                assert_eq!(pixel.photopic_standard_error_lm(), [0.0; 4]);
            }
        }
        assert_eq!(image.cropped().detected_paths(), 0);
        close(
            image.pixel_irradiance(4).unwrap().photopic_lux.i,
            source.flux_lm() / 0.0001,
            1e-4,
        );
    });
}

#[test]
fn g1_chromatic_refraction_places_each_wavelength_in_its_snell_pixel() {
    with_cx(|cx, _| {
        let theta = 0.6_f64;
        let source = PhotometricSource::try_new(
            RigidTransform::try_new(
                [0.0, (theta / 2.0).sin(), 0.0, (theta / 2.0).cos()],
                [0.0, 0.0, -1.0],
            )
            .unwrap(),
            900.0,
            1.0,
            AngularDistribution::uniform_cone(1e-6).unwrap(),
            SourceSpectrum::equal_energy(),
        )
        .unwrap();
        let geometry =
            InstanceScene::try_new(vec![square(1, 0.0, 2.0, false), square(9, 1.0, 2.0, true)])
                .unwrap();
        let air = glass(1.0, 0.0);
        let material = DielectricGlass::new(
            CauchyIor::try_new(1.3, 0.08, 0.0).unwrap(),
            BeerLambertAbsorption::try_constant(0.0).unwrap(),
            GlassProvenance::Custom,
        );
        let media = [
            SceneMedium::Homogeneous(&air),
            SceneMedium::Homogeneous(&material),
        ];
        let scene = CurvedPolarizedScene::try_new(
            cx,
            &geometry,
            &media,
            vec![
                SurfaceBinding {
                    object_id: 1,
                    interaction: Interaction::Dielectric {
                        outside: 0,
                        inside: 1,
                        surface: DielectricSurface::SMOOTH,
                    },
                },
                detector(),
            ],
            100_000,
        )
        .unwrap();
        let image = source
            .render_detector(
                cx,
                &scene,
                request(4),
                raster([1.06, 0.0, 1.0], 0.14, 0.02, 7, 1),
            )
            .unwrap();
        let mut expected_w = [0.0; 7];
        let mut expected_xyz = [[0.0; 3]; 7];
        // Independent Snell construction, with1m in air then1m in glass:
        // x=tan(theta)+sin(theta)/sqrt(n(lambda)^2-sin(theta)^2).
        // The tiny source cone cannot cross a pixel edge: the closest spectral
        // midpoint is >5e-5m away, versus <4e-6m angular/contact displacement.
        for bin in 0..SPECTRAL_BINS {
            let lambda_um = wavelength_nm(bin) * 0.001;
            let n = 1.3 + 0.08 / (lambda_um * lambda_um);
            let x = theta.tan() + theta.sin() / (n * n - theta.sin().powi(2)).sqrt();
            let column = ((x - 0.99) / 0.02).floor() as usize;
            assert!(column < 7);
            let edge_distance = (0..=7)
                .map(|edge| (x - (0.99 + edge as f64 * 0.02)).abs())
                .fold(f64::INFINITY, f64::min);
            assert!(edge_distance > 5e-5);
            // Use actually measured transmitted flux, so the independent
            // spatial oracle does not assume which stochastic branch survived.
            let power = image.total_detector().spectral_power_w_nm()[bin].i;
            expected_w[column] += power * BIN_WIDTH_NM;
            add_xyz(&mut expected_xyz[column], bin, power);
        }
        assert_eq!(image.cropped().detected_paths(), 0);
        for (column, pixel) in image.pixels().iter().enumerate() {
            close(pixel.radiant_stokes_w().i, expected_w[column], 2e-11);
            for channel in 0..3 {
                close(
                    pixel.xyz_flux()[channel],
                    expected_xyz[column][channel],
                    2e-8,
                );
            }
            close(pixel.photopic_stokes_lm().i, expected_xyz[column][1], 2e-8);
        }
        assert!(image.pixels()[0].xyz_flux()[2] > 0.0);
        assert_eq!(image.pixels()[6].xyz_flux()[2], 0.0);
    });
}

#[test]
fn g1_image_preserves_source_filter_polarization_and_absolute_xyz() {
    with_cx(|cx, _| {
        let original = source(AngularDistribution::uniform_cone(0.1).unwrap());
        let filter = core::array::from_fn(|bin| if bin % 2 == 0 { 0.25 } else { 0.75 });
        let source = original
            .clone()
            .with_unpolarized_optics()
            .through(&SpectralOptics::neutral(filter).unwrap())
            .unwrap()
            .through(&SpectralOptics::linear_polarizer(0.0).unwrap())
            .unwrap();
        let geometry = InstanceScene::try_new(vec![square(9, 0.0, 2.0, true)]).unwrap();
        let air = glass(1.0, 0.0);
        let media = [SceneMedium::Homogeneous(&air)];
        let scene = CurvedPolarizedScene::try_new(cx, &geometry, &media, vec![detector()], 100_000)
            .unwrap();
        let layout = raster([0.0; 3], 2.0, 2.0, 1, 1);
        let image = source
            .render_detector(cx, &scene, request(4), layout)
            .unwrap();
        let mut expected_xyz = [0.0; 3];
        for bin in 0..SPECTRAL_BINS {
            let expected = original.flux_lm()
                * original.spectrum().bins_w_per_lumen_nm()[bin]
                * filter[bin]
                * 0.5;
            add_xyz(&mut expected_xyz, bin, expected);
        }
        let pixel = &image.pixels()[0];
        close(pixel.photopic_stokes_lm().i, expected_xyz[1], 2e-9);
        close(pixel.photopic_stokes_lm().q, expected_xyz[1], 2e-9);
        close(pixel.photopic_stokes_lm().u, 0.0, 2e-9);
        close(pixel.photopic_stokes_lm().v, 0.0, 2e-9);
        for channel in 0..3 {
            close(pixel.xyz_flux()[channel], expected_xyz[channel], 2e-9);
        }
        assert_eq!(pixel.received_samples(), 4);
        let limited = FiniteDetectorRequest {
            trace: crate::polarization::dielectric::scene::curved::CurvedTraceSettings {
                max_distance_m: 0.1,
                ..request(4).trace
            },
            ..request(4)
        };
        assert!(matches!(
            source.render_detector(cx, &scene, limited, layout),
            Err(FiniteDetectorError::IncompletePath {
                sample: 0,
                bin: 0,
                termination: PathTermination::DistanceLimit
            })
        ));
    });
}
