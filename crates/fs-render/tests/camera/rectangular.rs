//! G0/G1/G3: finite slit geometry, separable Fraunhofer and Fresnel oracles,
//! and routing through the same coherent and source-depth camera owners.
use super::*;
use fs_math::c64::C64;
use std::f64::consts::{PI, TAU};

fn rotate(value: [f64; 2], angle: f64) -> [f64; 2] {
    let (sin, cos) = angle.sin_cos();
    [
        cos * value[0] - sin * value[1],
        sin * value[0] + cos * value[1],
    ]
}

fn sinc(value: f64) -> f64 {
    if value == 0.0 {
        1.0
    } else {
        value.sin() / value
    }
}

#[test]
fn rectangular_geometry_retains_physical_dimensions_rotation_and_admission() {
    let (width, height, rotation) = (0.002, 0.0005, 0.37);
    let aperture = Aperture::try_rectangular(width, height, rotation).unwrap();
    assert_eq!(aperture.rectangular_dimensions_m(), Some([width, height]));
    assert_eq!(aperture.rotation_rad(), Some(rotation));
    assert_eq!(aperture.blades(), None);
    assert_close(
        aperture.radius_m(),
        0.5 * width.hypot(height),
        1e-18,
        "circumradius",
    );
    assert_close(
        aperture.transmitting_area_m2().unwrap(),
        width * height,
        1e-20,
        "area",
    );
    let camera = axial_camera(aperture, 0.25);
    with_cx(|_, cx| {
        for [u, v] in [[0.0, 0.0], [0.5, 0.5], [0.2, 0.8], [0.99, 0.01]] {
            let ray = camera
                .generate_ray_from_tangent_offsets(cx, 0.0, 0.0, LensSample::try_new(u, v).unwrap())
                .unwrap();
            let expected = rotate([width * (u - 0.5), height * (v - 0.5)], rotation);
            assert_point_close(
                ray.origin,
                Point3::new(expected[0], expected[1], 0.0),
                "uniform rectangle sampling",
            );
        }
    });
    for [w, h, angle] in [
        [0.0, height, rotation],
        [-width, height, rotation],
        [width, f64::INFINITY, rotation],
        [width, height, f64::NAN],
        [1.0, 1.0 / 1025.0, 0.0],
        [f64::MIN_POSITIVE, f64::MIN_POSITIVE, 0.0],
    ] {
        assert!(Aperture::try_rectangular(w, h, angle).is_err());
    }
    assert!(Aperture::try_rectangular(1.0, 1.0 / 1024.0, 0.0).is_ok());
    // The short slit dimension controls the scalar-wave admission. Its long
    // diagonal must not disguise a subwavelength opening.
    let too_narrow = axial_camera(Aperture::try_rectangular(0.001, 2e-6, 0.0).unwrap(), 0.25)
        .with_in_focus_diffraction(FocusedDiffraction::try_new(0.0).unwrap())
        .unwrap();
    assert!(
        too_narrow
            .in_focus_diffraction_density([0.0; 2], 550.0)
            .is_err()
    );
}

#[test]
fn rectangular_fraunhofer_keeps_signed_sinc_lobes_and_anisotropic_zeros() {
    let width: f64 = 0.008;
    let rotation = 0.37;
    for aspect in [1.0, 4.0, 64.0] {
        let height = width / aspect;
        let camera = axial_camera(
            Aperture::try_rectangular(width, height, rotation).unwrap(),
            0.25,
        )
        .with_in_focus_diffraction(FocusedDiffraction::try_new(0.0).unwrap())
        .unwrap();
        for wavelength_nm in [400.0, 550.0, 700.0] {
            let lambda = wavelength_nm * 1e-9;
            let peak = (width * height).sqrt() / lambda;
            for [x, y] in [
                [0.0, 0.0],
                [1e-7, -2e-7],
                [0.2, -0.7],
                [1.0, 0.0],
                [0.0, 1.0],
                [1.5, 0.0],
                [2.5, -1.8],
            ] {
                let offset = rotate([x * lambda / width, y * lambda / height], rotation);
                let expected = peak * sinc(PI * x) * sinc(PI * y);
                let amplitude = camera
                    .in_focus_diffraction_amplitude(offset, wavelength_nm)
                    .unwrap();
                let density = camera
                    .in_focus_diffraction_density(offset, wavelength_nm)
                    .unwrap();
                assert_close(
                    amplitude.re,
                    expected,
                    2e-12 * peak,
                    "rectangular sinc transfer",
                );
                assert_eq!(amplitude.im, 0.0);
                assert_close(
                    density,
                    expected * expected,
                    4e-12 * peak * peak,
                    "rectangular power",
                );
            }
        }
    }
}

// Independent composite Simpson integration of the one-dimensional pupil.
// No native edge formula, quadrature rule or Fresnel implementation is used.
fn fresnel_average(curvature: f64, frequency: f64) -> C64 {
    const INTERVALS: usize = 4096;
    let mut sum = C64::ZERO;
    for index in 0..=INTERVALS {
        let t = -1.0 + 2.0 * index as f64 / INTERVALS as f64;
        let phase = curvature * t * t - frequency * t;
        let (sin, cos) = phase.sin_cos();
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
fn rectangular_defocus_and_layered_fields_match_separable_fresnel_integrals() {
    let (width, height, rotation): (f64, f64, f64) = (0.0016, 0.0006, 0.37);
    let aperture = Aperture::try_rectangular(width, height, rotation).unwrap();
    let radius = aperture.radius_m();
    let focus = 0.25;
    let lambda = 550e-9;
    let scale = lambda / (TAU * radius);
    let half = [0.5 * width / radius, 0.5 * height / radius];
    let peak = (width * height).sqrt() / lambda;
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
            let offset = rotate(q.map(|v| v * scale), rotation);
            let expected = (fresnel_average(phase * half[0] * half[0], q[0] * half[0])
                * fresnel_average(phase * half[1] * half[1], q[1] * half[1]))
            .scale(peak);
            let actual = camera
                .defocused_diffraction_amplitude(offset, 550.0)
                .unwrap();
            assert_close(
                actual.re,
                expected.re,
                4e-10 * peak,
                "Fresnel real transfer",
            );
            assert_close(
                actual.im,
                expected.im,
                4e-10 * peak,
                "Fresnel imaginary transfer",
            );
            assert_close(
                layered
                    .layered_diffraction_density(offset, depth, 550.0)
                    .unwrap(),
                actual.norm_sq(),
                4e-11 * peak * peak,
                "layered rectangular pupil",
            );
        }
    }
}
