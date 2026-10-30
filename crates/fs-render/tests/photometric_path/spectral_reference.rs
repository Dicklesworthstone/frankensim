//! Independent CIE references through filters, materials, camera paths and NEE.
//! Tolerances qualify these represented spectra, not real-fixture calibration.

use super::*;
use std::sync::OnceLock;

use fs_render::photometry::area::{EmissionSidedness, PhotometricAreaEmitter, PhotometricEmitters};
use fs_render::photometry::medium::HomogeneousMedium;
use fs_render::photometry::path::{IlluminanceSampling, PathIlluminanceJob};
use fs_render::photometry::{CalculationGrid, IlluminanceJob, IlluminanceProbe, PhotometricScene};

/// Independent reference from the retained original one-nanometre CIE table.
/// This deliberately does not use any production observer or source helper.
pub(super) fn cie_reference_xyz(wavelength: f64) -> [f64; 3] {
    static ROWS: OnceLock<Vec<[f64; 4]>> = OnceLock::new();
    let rows = ROWS.get_or_init(|| {
        include_str!("../../data/cie/CIE_xyz_1931_2deg.csv")
            .lines()
            .map(|line| {
                let values: Vec<f64> = line
                    .split(',')
                    .map(|value| value.parse().expect("published CIE decimal"))
                    .collect();
                values
                    .try_into()
                    .expect("wavelength and three CIE functions")
            })
            .collect()
    });
    if !wavelength.is_finite() || !(360.0..=830.0).contains(&wavelength) {
        return [0.0; 3];
    }
    let upper = rows.partition_point(|row| row[0] <= wavelength);
    let lower = upper - 1;
    if upper == rows.len() {
        return [rows[lower][1], rows[lower][2], rows[lower][3]];
    }
    let t = (wavelength - rows[lower][0]) / (rows[upper][0] - rows[lower][0]);
    core::array::from_fn(|channel| {
        (1.0 - t) * rows[lower][channel + 1] + t * rows[upper][channel + 1]
    })
}

fn two_band_spectrum() -> SourceSpectrum {
    let mut spectrum = [0.0; SPECTRAL_BINS];
    spectrum[13] = 1.0; // 445--450 nm, evaluated at 447.5 nm.
    spectrum[54] = 1.0; // 650--655 nm, evaluated at 652.5 nm.
    SourceSpectrum::try_from_bins(spectrum).unwrap()
}

fn point_lux(sources: &[([f64; 3], f64)], position: [f64; 3], normal: [f64; 3]) -> f64 {
    sources
        .iter()
        .map(|(source, candela)| {
            let d: [f64; 3] = core::array::from_fn(|axis| source[axis] - position[axis]);
            let r_squared = d.iter().map(|v| v * v).sum::<f64>();
            let projected = d.iter().zip(normal).map(|(d, n)| d * n).sum::<f64>();
            candela * projected.max(0.0) / (r_squared * r_squared.sqrt())
        })
        .sum()
}

#[test]
fn g1_cie_selective_filters_preserve_multi_fixture_si_illuminance() {
    with_cx(|cx, _| {
        // Published 447/448 and 652/653 nm averages. Equal radiant power in
        // the two bins does not mean equal luminous power.
        let blue_fraction = 0.033_702_4 / (0.033_702_4 + 0.093_655_8);
        let mut blue_filter = [0.0; SPECTRAL_BINS];
        blue_filter[13] = 1.0;
        let mut red_filter = [0.0; SPECTRAL_BINS];
        red_filter[54] = 1.0;
        let sources = [
            isotropic([0.0, 0.0, 6.0], 120.0, two_band_spectrum()),
            isotropic([3.0, 0.0, 6.0], 75.0, two_band_spectrum())
                .filtered(blue_filter)
                .unwrap(),
            isotropic([-3.0, 4.0, 6.0], 60.0, two_band_spectrum())
                .filtered(red_filter)
                .unwrap(),
        ];
        close(
            sources[1].flux_lm(),
            4.0 * PI * 75.0 * blue_fraction,
            2.0e-13,
            "CIE-weighted blue filter flux",
        );
        close(
            sources[2].flux_lm(),
            4.0 * PI * 60.0 * (1.0 - blue_fraction),
            2.0e-13,
            "CIE-weighted red filter flux",
        );
        let reference_sources = [
            ([0.0, 0.0, 6.0], 120.0),
            ([3.0, 0.0, 6.0], 75.0 * blue_fraction),
            ([-3.0, 4.0, 6.0], 60.0 * (1.0 - blue_fraction)),
        ];
        let scene = PhotometricScene::try_new(
            &[],
            &sources,
            HomogeneousMedium::vacuum(),
            TransportOptions::default(),
            cx,
        )
        .unwrap();
        let grid = CalculationGrid::try_new(
            Point3::new(-1.0, -0.5, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            3,
            2,
        )
        .unwrap();
        let mut chunked = IlluminanceJob::for_grid(&scene, grid).unwrap();
        assert_eq!(chunked.advance(cx, 2).unwrap(), 2);
        assert_eq!(chunked.advance(cx, usize::MAX).unwrap(), 4);
        let mut whole = IlluminanceJob::for_grid(&scene, grid).unwrap();
        whole.advance(cx, usize::MAX).unwrap();
        for (index, (actual, replay)) in chunked.results().iter().zip(whole.results()).enumerate() {
            let position = [-1.0 + (index % 3) as f64, -0.5 + (index / 3) as f64, 0.0];
            close(
                actual.lux,
                point_lux(&reference_sources, position, [0.0, 0.0, 1.0]),
                2.0e-13,
                "independent filtered stage-grid lux",
            );
            assert_eq!(actual.lux.to_bits(), replay.lux.to_bits());
            assert!(actual.total_error_bound_lux.is_none());
        }
        // A second height and an inclined receiver exercise actual metre
        // distances and cosine projection independently of spectral weights.
        for normal in [[0.0, 0.0, 1.0], [0.6, 0.0, 0.8], [0.0, 0.0, -1.0]] {
            let receiver = IlluminanceProbe::try_new(
                Point3::new(0.0, 0.0, 3.0),
                Vec3::new(normal[0], normal[1], normal[2]),
            )
            .unwrap();
            close(
                scene.illuminance(cx, &receiver).unwrap().lux,
                point_lux(&reference_sources, [0.0, 0.0, 3.0], normal),
                2.0e-13,
                "inclined multi-fixture receiver",
            );
        }
    });
}

#[test]
fn g1_colored_lambertian_paths_match_independent_cie_xyz_and_reflectance() {
    with_cx(|cx, _| {
        let mut scene = receiver_scene();
        scene.primitives[0].material = Material::Lambertian {
            reflectance: LiftedSpectrum {
                c: [0.0, 2.0, -0.5],
            },
        };
        let sources = [isotropic([0.0, 0.0, 3.0], 90.0, two_band_spectrum())];
        let estimate = render_pixel(
            &scene,
            &sources,
            PathMedium::Vacuum,
            settings(1, 65_536),
            cx,
        );
        let wavelengths = [447.5, 652.5];
        let y_sum: f64 = wavelengths
            .map(|lambda| cie_reference_xyz(lambda)[1])
            .iter()
            .sum();
        let mut expected = [0.0; 3];
        for lambda in wavelengths {
            let p = 2.0 * (lambda - 380.0) / 400.0 - 0.5;
            let rho = 0.5 + 0.5 * p / (1.0 + p * p).sqrt();
            // 90 cd at 3 m gives 10 lx before spectral reflection.
            for (value, observer) in expected.iter_mut().zip(cie_reference_xyz(lambda)) {
                *value += (10.0 / PI) * rho * observer / y_sum;
            }
        }
        for (actual, reference) in estimate.xyz.into_iter().zip(expected) {
            close(actual, reference, 0.002, "absolute reflected CIE XYZ");
        }
        assert!(estimate.estimated_standard_error_cd_m2.is_none());
        assert!(estimate.total_error_bound_cd_m2.is_none());
    });
}

#[test]
fn g1_blue_finite_emitter_camera_and_probe_nee_use_the_same_cie_observer() {
    with_cx(|cx, _| {
        let mut scene = empty_scene([0.0, 0.0, 2.0], [0.0, 0.0, -1.0]);
        scene
            .primitives
            .push(floor(0.0, 1.0, LiftedSpectrum { c: [0.0; 3] }));
        let mut bins = [0.0; SPECTRAL_BINS];
        bins[13] = 1.0;
        let surfaces = [PhotometricAreaEmitter::try_new(
            0,
            100.0 * PI,
            1.0,
            SourceSpectrum::try_from_bins(bins).unwrap(),
            EmissionSidedness::Front,
        )
        .unwrap()];
        let emitters = PhotometricEmitters {
            points: &[],
            surfaces: &surfaces,
        };
        let mut image = PathLuminanceJob::new_with_emitters(
            &scene,
            emitters,
            PathMedium::Vacuum,
            TransportOptions::default(),
            settings(1, 65_536),
            cx,
        )
        .unwrap();
        image.advance(cx, usize::MAX).unwrap();
        let image_value = image.pixel_estimate(0).unwrap();
        // Actual area=4 m², one emitting hemisphere: L=Phi/(pi*A)=25 cd/m².
        let blue = cie_reference_xyz(447.5);
        for (actual, observer) in image_value.xyz.into_iter().zip(blue) {
            close(
                actual,
                25.0 * observer / blue[1],
                0.002,
                "finite blue emitter XYZ",
            );
        }
        let probes =
            [
                IlluminanceProbe::try_new(Point3::new(0.0, 0.0, 2.0), Vec3::new(0.0, 0.0, -1.0))
                    .unwrap(),
            ];
        let mut probe = PathIlluminanceJob::new_with_emitters(
            &scene,
            emitters,
            PathMedium::Vacuum,
            &probes,
            TransportOptions::default(),
            IlluminanceSampling {
                samples_per_probe: 65_536,
                max_depth: 1,
                sampler: Sampler::OwenSobol,
                seed: 173,
            },
            cx,
        )
        .unwrap();
        probe.advance(cx, usize::MAX).unwrap();
        let probe_value = probe.probe_estimate(0).unwrap();
        // Closed-form projected solid angle of a 2x2 m square viewed axially
        // from 2 m. This independently checks both the hit and area-NEE arms.
        let projected_solid_angle = 4.0 / 5.0_f64.sqrt() * (1.0 / 5.0_f64.sqrt()).atan();
        close(
            probe_value.continuous_lux,
            25.0 * projected_solid_angle,
            0.02,
            "finite blue emitter probe illuminance",
        );
        assert_eq!(probe_value.direct_lux, 0.0);
        assert!(image_value.total_error_bound_cd_m2.is_none());
        assert!(probe_value.total_error_bound_lux.is_none());
        assert!(probe_value.estimated_standard_error_lux.is_none());
    });
}
