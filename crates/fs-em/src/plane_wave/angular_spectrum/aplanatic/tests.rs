use super::super::modes::{fields, normal_flux};
use super::*;

fn area(grid: AngularSpectrumGrid) -> f64 {
    (grid.width as f64 * grid.spacing_m[0]) * (grid.height as f64 * grid.spacing_m[1])
}

fn sample(
    grid: AngularSpectrumGrid,
    pupil: AplanaticPupil,
    shape: impl FnMut(AplanaticPupilSample) -> [C64; 2],
) -> Vec<[C64; 2]> {
    pupil
        .sample_modes(grid, area(grid), 1.0, 1.0, shape, &mut || false)
        .unwrap()
}

fn field_at(
    grid: AngularSpectrumGrid,
    coefficients: &[[C64; 2]],
    position: [f64; 3],
) -> ([C64; 3], [C64; 3]) {
    let mut electric = [C64::ZERO; 3];
    let mut magnetic = [C64::ZERO; 3];
    for (index, &coefficient) in coefficients.iter().enumerate() {
        if coefficient == [C64::ZERO; 2] {
            continue;
        }
        let kx = transverse_number(index % grid.width, grid.width, grid.spacing_m[0]);
        let ky = transverse_number(index / grid.width, grid.height, grid.spacing_m[1]);
        let kz = incident_cosine(1.0, kx, ky).unwrap();
        let (e, h) = fields(kx, ky, C64::new(kz, 0.0), 1.0, coefficient).unwrap();
        let angle = kx * position[0] + ky * position[1] + kz * position[2];
        let phase = C64::new(det::cos(angle), -det::sin(angle));
        for axis in 0..3 {
            electric[axis] = electric[axis] + e[axis] * phase;
            magnetic[axis] = magnetic[axis] + h[axis] * phase;
        }
    }
    (electric, magnetic)
}

fn close(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:?} != {b:?}, error {}",
        (a - b).abs()
    );
}

#[test]
fn aplanatic_rotation_preserves_transversality_and_requested_vector_power() {
    // G0/G1: verify the lens map through all six original Maxwell components,
    // including each oblique longitudinal field and its actual Poynting flux.
    let grid = AngularSpectrumGrid {
        width: 31,
        height: 27,
        spacing_m: [2.5, 2.0],
    };
    let pupil = AplanaticPupil::try_new(0.95, [0.0; 3], 0.17).unwrap();
    let coefficients = sample(grid, pupil, |point| {
        let [x, y, _] = point.direction;
        [
            C64::new(1.0 + 0.3 * x, 0.4 * y),
            C64::new(-0.2 + 0.5 * y, 0.6 + 0.1 * x),
        ]
    });
    let mut power = 0.0;
    let mut longitudinal = 0.0_f64;
    for (index, &coefficient) in coefficients.iter().enumerate() {
        if coefficient == [C64::ZERO; 2] {
            continue;
        }
        let kx = transverse_number(index % grid.width, grid.width, grid.spacing_m[0]);
        let ky = transverse_number(index / grid.width, grid.height, grid.spacing_m[1]);
        let kz = incident_cosine(1.0, kx, ky).unwrap();
        let (e, h) = fields(kx, ky, C64::new(kz, 0.0), 1.0, coefficient).unwrap();
        close(
            e[0].scale(kx) + e[1].scale(ky) + e[2].scale(kz),
            C64::ZERO,
            2e-18,
        );
        close(
            h[0].scale(kx) + h[1].scale(ky) + h[2].scale(kz),
            C64::ZERO,
            2e-18,
        );
        longitudinal = longitudinal.max(e[2].abs());
        power += area(grid) * normal_flux(e, h);
    }
    assert!(longitudinal > 1e-4);
    assert!(
        (power / pupil.incident_power_w() - 1.0).abs() < 3e-14,
        "{power}"
    );
}

#[test]
fn radial_and_azimuthal_pupils_resolve_longitudinal_focus() {
    // G1/G3: rotational symmetry leaves radial illumination purely electric
    // longitudinal on axis and azimuthal illumination magnetic longitudinal.
    let grid = AngularSpectrumGrid {
        width: 63,
        height: 63,
        spacing_m: [2.0; 2],
    };
    let pupil = AplanaticPupil::try_new(0.95, [0.0; 3], 1.0).unwrap();
    let radial = sample(grid, pupil, |point| {
        let [x, y, _] = point.direction;
        [C64::new(x, 0.0), C64::new(y, 0.0)]
    });
    let azimuthal = sample(grid, pupil, |point| {
        let [x, y, _] = point.direction;
        [C64::new(-y, 0.0), C64::new(x, 0.0)]
    });
    let (er, hr) = field_at(grid, &radial, [0.0; 3]);
    let (ea, ha) = field_at(grid, &azimuthal, [0.0; 3]);
    assert!(er[2].re < -0.2 && ha[2].re > 0.2);
    for z in [
        er[0], er[1], hr[0], hr[1], hr[2], ea[0], ea[1], ea[2], ha[0], ha[1],
    ] {
        close(z, C64::ZERO, 3e-15);
    }
    close(er[2], -ha[2], 5e-15);
}

#[test]
fn focus_translation_and_authored_phase_remain_coherent() {
    let grid = AngularSpectrumGrid {
        width: 17,
        height: 19,
        spacing_m: [2.1, 2.2],
    };
    let focus = [1.3, -0.8, 7.25];
    let offset = [0.2, 0.35, -0.4];
    let shape = |point: AplanaticPupilSample| {
        let [x, y] = point.pupil_coordinate;
        [
            C64::new(1.0 - 0.2 * x * x, 0.3 * y),
            C64::new(-0.1 * y, 0.4 + 0.1 * x),
        ]
    };
    let reference = sample(
        grid,
        AplanaticPupil::try_new(0.9, [0.0; 3], 0.4).unwrap(),
        shape,
    );
    let shifted = sample(
        grid,
        AplanaticPupil::try_new(0.9, focus, 0.4).unwrap(),
        shape,
    );
    let phase = C64::new(det::cos(0.73), det::sin(0.73));
    let phased = sample(
        grid,
        AplanaticPupil::try_new(0.9, focus, 0.4).unwrap(),
        |point| shape(point).map(|z| z * phase),
    );
    let target = std::array::from_fn(|axis| focus[axis] + offset[axis]);
    let (a, b) = field_at(grid, &reference, offset);
    let (c, d) = field_at(grid, &shifted, target);
    let (e, f) = field_at(grid, &phased, target);
    for axis in 0..3 {
        close(a[axis], c[axis], 2e-15);
        close(b[axis], d[axis], 2e-15);
        close(a[axis] * phase, e[axis], 2e-15);
        close(b[axis] * phase, f[axis], 2e-15);
    }
}

#[test]
fn normalized_jones_shape_survives_extreme_authored_magnitudes() {
    let grid = AngularSpectrumGrid {
        width: 11,
        height: 13,
        spacing_m: [2.0; 2],
    };
    let pupil = AplanaticPupil::try_new(0.9, [0.0; 3], 0.8).unwrap();
    let reference = sample(grid, pupil, |_| [C64::ONE, C64::ZERO]);
    for magnitude in [f64::from_bits(1), 1e-300, 1e300] {
        let actual = sample(grid, pupil, |_| [C64::new(magnitude, 0.0), C64::ZERO]);
        assert_eq!(actual, reference);
    }
}

#[test]
fn uniform_high_na_peak_converges_to_the_richards_wolf_integral() {
    // G1/G2: independently integrate the continuous x-polarized pupil on axis.
    // The angular integral is pi*int_(cos(alpha))^1 (1+c)*sqrt(c) dc.
    // Requested watts normalize its continuum pupil area pi*sin(alpha)^2.
    let sine = 0.9;
    let cosine = det::sqrt(1.0 - sine * sine);
    let root = det::sqrt(cosine);
    let integral = core::f64::consts::PI
        * ((2.0 / 5.0) * (1.0 - cosine * cosine * root) + (2.0 / 3.0) * (1.0 - cosine * root));
    let expected =
        integral / core::f64::consts::TAU * det::sqrt(2.0 / (core::f64::consts::PI * sine * sine));
    let mut errors = Vec::new();
    for n in [16, 64, 256] {
        let grid = AngularSpectrumGrid {
            width: n,
            height: n,
            spacing_m: [2.0; 2],
        };
        let pupil = AplanaticPupil::try_new(sine, [0.0; 3], 1.0).unwrap();
        let coefficients = sample(grid, pupil, |_| [C64::ONE, C64::ZERO]);
        let (e, _) = field_at(grid, &coefficients, [0.0; 3]);
        close(e[1], C64::ZERO, 2e-15);
        close(e[2], C64::ZERO, 2e-15);
        errors.push((e[0].re / expected - 1.0).abs());
    }
    assert!(
        errors[2] < 3e-4 && errors[2] < errors[1] && errors[1] < errors[0],
        "{errors:?}"
    );
}

#[test]
fn angular_domain_and_cancellation_refuse_before_unbounded_work() {
    let grid = AngularSpectrumGrid {
        width: 17,
        height: 19,
        spacing_m: [2.0; 2],
    };
    for sine in [0.0, -0.1, 1.0, 1.1, f64::NAN] {
        assert!(AplanaticPupil::try_new(sine, [0.0; 3], 1.0).is_err());
    }
    assert!(AplanaticPupil::try_new(0.9, [f64::INFINITY, 0.0, 0.0], 1.0).is_err());
    assert!(AplanaticPupil::try_new(0.9, [0.0; 3], 0.0).is_err());
    let pupil = AplanaticPupil::try_new(0.9, [0.0; 3], 1.0).unwrap();
    let clipped = AngularSpectrumGrid {
        spacing_m: [4.0, 2.0],
        ..grid
    };
    let mut calls = 0;
    assert!(
        pupil
            .sample_modes(
                clipped,
                area(clipped),
                1.0,
                1.0,
                |_| {
                    calls += 1;
                    [C64::ONE, C64::ZERO]
                },
                &mut || false
            )
            .is_err()
    );
    assert_eq!(calls, 0);
    let mut polls = 0;
    let result = pupil.sample_modes(
        grid,
        area(grid),
        1.0,
        1.0,
        |_| {
            calls += 1;
            [C64::ONE, C64::ZERO]
        },
        &mut || {
            polls += 1;
            polls == 3
        },
    );
    assert!(matches!(result, Err(AngularSpectrumError::Cancelled)));
    assert!(calls <= 32);
    assert!(
        pupil
            .sample_modes(grid, area(grid), 1.0, 1.0, |_| [C64::ZERO; 2], &mut || {
                false
            })
            .is_err()
    );
    assert!(
        pupil
            .sample_modes(
                grid,
                area(grid),
                1.0,
                1.0,
                |_| [C64::new(f64::NAN, 0.0); 2],
                &mut || false
            )
            .is_err()
    );
}
