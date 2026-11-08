//! G0/G1/G3: disconnected physical pupils, Young fringes and absolute Fresnel phase.
use super::*;
use fs_math::c64::C64;
use std::f64::consts::{PI, TAU};

fn rotate(q: [f64; 2], angle: f64) -> [f64; 2] {
    let (sin, cos) = angle.sin_cos();
    [cos * q[0] - sin * q[1], sin * q[0] + cos * q[1]]
}

fn sinc(q: f64) -> f64 {
    if q == 0.0 { 1.0 } else { q.sin() / q }
}

#[test]
fn rectangular_pair_samples_only_its_two_physical_openings() {
    let (width, height, separation, rotation) = (0.0005, 0.0006, 0.0015, 0.37);
    let aperture = Aperture::try_rectangular_pair(width, height, separation, rotation).unwrap();
    assert_eq!(
        aperture.rectangular_pair_dimensions_m(),
        Some([width, height, separation])
    );
    assert_eq!(aperture.rectangular_dimensions_m(), None);
    assert_eq!(aperture.rotation_rad(), Some(rotation));
    assert_eq!(aperture.blades(), None);
    assert_close(
        aperture.transmitting_area_m2().unwrap(),
        2.0 * width * height,
        1e-21,
        "open area",
    );
    assert_close(
        aperture.radius_m(),
        0.5 * (width + separation).hypot(height),
        1e-18,
        "outer radius",
    );
    let camera = axial_camera(aperture, 0.25);
    with_cx(|_, cx| {
        for (u, expected_x) in [
            (0.0, -0.5 * (width + separation)),
            (0.25, -0.5 * separation),
            (0.5, 0.5 * (separation - width)),
            (0.75, 0.5 * separation),
        ] {
            let ray = camera
                .generate_ray_from_tangent_offsets(
                    cx,
                    0.0,
                    0.0,
                    LensSample::try_new(u, 0.5).unwrap(),
                )
                .unwrap();
            let expected = rotate([expected_x, 0.0], rotation);
            assert_point_close(
                ray.origin,
                Point3::new(expected[0], expected[1], 0.0),
                "paired opening sample",
            );
        }
    });
    for separation in [0.0, 0.5 * width, f64::NAN, f64::INFINITY, 1024.0 * width] {
        assert!(Aperture::try_rectangular_pair(width, height, separation, rotation).is_err());
    }
    assert!(Aperture::try_rectangular_pair(width, height, width, rotation).is_ok());
    let narrow = axial_camera(
        Aperture::try_rectangular_pair(1e-6, 2e-6, 1e-4, 0.0).unwrap(),
        0.25,
    )
    .with_in_focus_diffraction(FocusedDiffraction::try_new(0.0).unwrap())
    .unwrap();
    assert!(
        narrow
            .in_focus_diffraction_density([0.0; 2], 550.0)
            .is_err()
    );
}

#[test]
fn rectangular_pair_retains_young_fringe_nulls_and_signed_amplitude() {
    let (width, height, separation, rotation): (f64, f64, f64, f64) =
        (0.0005, 0.0006, 0.0015, 0.37);
    let camera = axial_camera(
        Aperture::try_rectangular_pair(width, height, separation, rotation).unwrap(),
        0.25,
    )
    .with_in_focus_diffraction(FocusedDiffraction::try_new(0.0).unwrap())
    .unwrap();
    for wavelength in [400.0, 550.0, 700.0] {
        let lambda = wavelength * 1e-9;
        let peak = (2.0 * width * height).sqrt() / lambda;
        for [x, y] in [
            [0.0, 0.0],
            [0.5, 0.0],
            [1.0, 0.0],
            [1.5, 0.7],
            [2.0, -0.4],
            [4.25, 1.3],
        ] {
            // x is measured in fringe periods, y in the single-height zeros.
            let offset = rotate([lambda * x / separation, lambda * y / height], rotation);
            let expected = peak * sinc(PI * x * width / separation) * sinc(PI * y) * (PI * x).cos();
            let actual = camera
                .in_focus_diffraction_amplitude(offset, wavelength)
                .unwrap();
            assert_close(actual.re, expected, 3e-12 * peak, "signed pair transfer");
            assert_close(actual.im, 0.0, 1e-15 * peak, "symmetric pair phase");
            assert_close(
                camera
                    .in_focus_diffraction_density(offset, wavelength)
                    .unwrap(),
                expected * expected,
                8e-12 * peak * peak,
                "Young interference power",
            );
        }
        let null = camera
            .in_focus_diffraction_density(
                rotate([0.5 * lambda / separation, 0.0], rotation),
                wavelength,
            )
            .unwrap();
        assert!(null < 1e-25 * peak * peak);
    }
}

// Independent area-coordinate integration: the translation stays inside the
// squared physical coordinate, rather than using the implementation's shift theorem.
fn translated_average(center: f64, half_width: f64, q: f64, phase: f64) -> C64 {
    const INTERVALS: usize = 4096;
    let mut sum = C64::ZERO;
    for index in 0..=INTERVALS {
        let coordinate = center + half_width * (-1.0 + 2.0 * index as f64 / INTERVALS as f64);
        let (sin, cos) = (phase * coordinate * coordinate - q * coordinate).sin_cos();
        let weight = if index == 0 || index == INTERVALS {
            1.0
        } else if index % 2 == 0 {
            2.0
        } else {
            4.0
        };
        sum = sum + C64::new(cos, sin).scale(weight);
    }
    sum.scale(1.0 / (3.0 * INTERVALS as f64))
}

#[test]
fn rectangular_pair_fresnel_and_source_depth_paths_keep_translation_phase() {
    let (width, height, separation, rotation): (f64, f64, f64, f64) =
        (0.0005, 0.0006, 0.0015, 0.37);
    let aperture = Aperture::try_rectangular_pair(width, height, separation, rotation).unwrap();
    let radius = aperture.radius_m();
    let lambda = 550e-9;
    let scale = lambda / (TAU * radius);
    let peak = (2.0 * width * height).sqrt() / lambda;
    let focus = 0.25;
    for phase in [1e-6, 6.0, -7.0] {
        let depth = 1.0 / (1.0 / focus + phase * lambda / (PI * radius * radius));
        let camera = axial_camera(aperture.clone(), focus)
            .with_defocused_diffraction(DefocusedDiffraction::try_new(Some(depth), 0.0).unwrap())
            .unwrap();
        let depths = if depth < focus {
            [depth, focus]
        } else {
            [focus, depth]
        };
        let layered = axial_camera(aperture.clone(), focus)
            .with_layered_diffraction(LayeredDiffraction::try_new(&depths, 0.0).unwrap())
            .unwrap();
        for q in [[0.0; 2], [2e-7, -3e-7], [0.4, -1.1], [20.0, 1.0]] {
            let x = translated_average(
                -0.5 * separation / radius,
                0.5 * width / radius,
                q[0],
                phase,
            ) + translated_average(
                0.5 * separation / radius,
                0.5 * width / radius,
                q[0],
                phase,
            );
            let expected =
                (x * translated_average(0.0, 0.5 * height / radius, q[1], phase)).scale(0.5 * peak);
            let offset = rotate(q.map(|v| v * scale), rotation);
            let actual = camera
                .defocused_diffraction_amplitude(offset, 550.0)
                .unwrap();
            assert_close(actual.re, expected.re, 5e-10 * peak, "paired Fresnel real");
            assert_close(
                actual.im,
                expected.im,
                5e-10 * peak,
                "paired Fresnel imaginary",
            );
            assert_close(
                layered
                    .layered_diffraction_density(offset, depth, 550.0)
                    .unwrap(),
                actual.norm_sq(),
                4e-11 * peak * peak,
                "paired source-depth density",
            );
        }
    }
}
