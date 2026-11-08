//! Independent physical-volume and signed-face controls for coherent heat.
use super::*;
use super::super::{
    Input, Incidence, OrientedCrystalFieldSide, coherent_illumination_field,
    coherent_illumination_poynting, solve, solve_in_plane,
};

const C0: f64 = 299_792_458.0;
const MU: f64 = 1.256_637_062_12e-6;
const LAMBDA: f64 = 550e-9;

fn input(axis: [f64; 3], indices: [C64; 2], p: f64) -> Input {
    Input {
        omega: core::f64::consts::TAU * C0 / LAMBDA,
        ni: 1.0,
        indices,
        mui: MU,
        mut_: MU,
        axis,
        incidence: Incidence::TangentialIndex(p),
        electric: [C64::new(0.7, 0.2), C64::new(-0.3, 0.5)],
        loss_f_m: indices.map(|n| -2.0 * n.re * n.im / (MU * C0 * C0)),
    }
}

fn limits() -> AbsorbingCrystalBoxLimits {
    AbsorbingCrystalBoxLimits {
        maximum_mode_pairs: 4096,
        maximum_workspace_bytes: 1 << 20,
    }
}

fn integrate(responses: &[AbsorbingCrystalResponse], region: AbsorbingCrystalBox) -> BoxIntegral {
    coherent_box_absorption(responses.iter(), region, limits(), 0, || false).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}; tolerance={tolerance:e}"
    );
}

#[test]
fn analytic_crossed_te_box_matches_independent_six_face_maxwell_fixture() {
    // Independently evaluated complex Fresnel TE waves in an isotropic
    // incident port, with ordinary-axis transmission in an optic-axis-z
    // crystal. The unequal signed directions retain lateral energy flow.
    let n = C64::new(1.4, -0.12);
    let mut data = input([0.0, 0.0, 1.0], [n, C64::new(1.9, -0.2)], 0.25);
    data.mui = 1e-6;
    data.mut_ = 1e-6;
    data.ni = C0 * det::sqrt(1e-11 * data.mui);
    data.loss_f_m = data.indices.map(|n| -2.0 * n.re * n.im / (data.mut_ * C0 * C0));
    data.electric = [C64::ZERO, C64::ONE];
    let first = solve(data).unwrap();
    data.incidence = Incidence::TangentialIndex(0.75);
    data.electric = [C64::ZERO, C64::new(0.4, 0.9)];
    let responses = [first, solve(data).unwrap()];
    let region = AbsorbingCrystalBox::try_new(
        [0.13, -0.17, 0.1].map(|v| v * LAMBDA),
        [0.31, 0.4, 0.9].map(|v| v * LAMBDA),
    )
    .unwrap();
    let result = integrate(&responses, region);
    let heat = 8.953_983_708_141_13e-17;
    let faces = [
        -4.737_108_596_954_948e-17,
        7.422_234_151_334_246e-17,
        0.0,
        0.0,
        -1.346_252_473_501_958_3e-16,
        1.823_415_472_499_167_8e-17,
    ];
    close(result.absorbed_power_w, heat, 3e-12 * heat);
    for (actual, expected) in result.outward_face_power_w.into_iter().zip(faces) {
        close(actual, expected, 3e-12 * heat);
    }
    close(result.residual_w, 0.0, 5e-12 * heat);
    let normal_drop = -result.outward_face_power_w[4] - result.outward_face_power_w[5];
    assert!((normal_drop - result.absorbed_power_w).abs() > 0.1 * heat);
    assert_eq!(result.mode_pair_count, 16);
}

#[test]
fn crossed_te_fringe_volume_and_phase_translation_match_closed_form() {
    let n = C64::new(1.7, -0.18);
    let p = 0.6;
    let mut data = input([0.0, 0.0, 1.0], [n; 2], p);
    data.electric = [C64::ZERO, C64::ONE];
    let first = solve(data).unwrap();
    let k0 = data.omega / C0;
    let qi = det::sqrt(1.0 - p * p);
    let q = (n * n - C64::from_re(p * p)).sqrt();
    let transmission = C64::from_re(2.0 * qi) / (C64::from_re(qi) + q);
    let alpha = -k0 * q.im;
    for phase in [0.0, 0.7, core::f64::consts::PI] {
        data.incidence = Incidence::TangentialIndex(-p);
        data.electric[1] = C64::new(det::cos(phase), det::sin(phase));
        let responses = [first, solve(data).unwrap()];
        for width in [0.015 * LAMBDA, 0.8 * LAMBDA, 3.0 * LAMBDA] {
            for depth in [0.02 * LAMBDA, 2.0 * LAMBDA, 30.0 * LAMBDA] {
                let origin = [-0.13 * LAMBDA, 0.2 * LAMBDA, 0.13 * LAMBDA];
                let height = 0.7 * LAMBDA;
                let region = AbsorbingCrystalBox::try_new(origin, [width, height, depth]).unwrap();
                let result = integrate(&responses, region);
                let delta = k0 * p * width;
                let intensity_x = 2.0 * width
                    * (1.0 + det::cos(2.0 * k0 * p * origin[0] + delta + phase)
                        * det::sin(delta) / delta);
                let depth_integral = det::exp(-2.0 * alpha * origin[2])
                    * (-det::expm1(-2.0 * alpha * depth)) / (2.0 * alpha);
                let expected = data.omega * (-n.re * n.im) / (MU * C0 * C0)
                    * transmission.norm_sq() * intensity_x * height * depth_integral;
                close(result.absorbed_power_w, expected, 3e-12 * expected);
                let entrance = q.re / (2.0 * C0 * MU) * transmission.norm_sq()
                    * intensity_x * height * det::exp(-2.0 * alpha * origin[2]);
                close(result.outward_face_power_w[4], -entrance, 3e-12 * entrance);
                close(
                    result.outward_face_power_w[5],
                    entrance * det::exp(-2.0 * alpha * depth),
                    3e-12 * entrance,
                );
                close(result.residual_w, 0.0, 5e-12 * entrance);
            }
        }
    }
}

#[test]
fn short_dark_fringe_keeps_cubic_heat_through_opaque_common_depth() {
    let n = C64::new(1.6, -0.12);
    let p = 0.6;
    let mut data = input([0.0, 0.0, 1.0], [n; 2], p);
    data.electric = [C64::ZERO, C64::ONE];
    let first = solve(data).unwrap();
    data.incidence = Incidence::TangentialIndex(-p);
    data.electric[1] = -C64::ONE;
    let responses = [first, solve(data).unwrap()];
    let k0 = data.omega / C0;
    let qi = det::sqrt(1.0 - p * p);
    let q = (n * n - C64::from_re(p * p)).sqrt();
    let transmission = C64::from_re(2.0 * qi) / (C64::from_re(qi) + q);
    let alpha = -k0 * q.im;
    let height = 0.7 * LAMBDA;
    let depth = 30.0 * LAMBDA;
    let depth_integral = -det::expm1(-2.0 * alpha * depth) / (2.0 * alpha);
    for delta in [1e-4, 1e-8, 1e-16] {
        let width = delta / (k0 * p);
        let region = AbsorbingCrystalBox::try_new(
            [-0.5 * width, 0.0, 0.0],
            [width, height, depth],
        )
        .unwrap();
        let result = integrate(&responses, region);
        // Integral_x |exp(-ikx)-exp(+ikx)|² =
        // width * delta²/3 * (1-delta²/20+...), with delta=k*width.
        let intensity_x = width * delta * delta / 3.0 * (1.0 - delta * delta / 20.0);
        let expected = data.omega * (-n.re * n.im) / (MU * C0 * C0)
            * transmission.norm_sq() * intensity_x * height * depth_integral;
        assert!(result.absorbed_power_w > 0.0);
        close(result.absorbed_power_w, expected, 2e-10 * expected);
    }
}

fn simpson_weight(index: usize, intervals: usize) -> f64 {
    if index == 0 || index == intervals { 1.0 } else if index % 2 == 0 { 2.0 } else { 4.0 }
}

// Independent direct physical-point quadrature. It does not use box pair
// kernels, window integrals or the production short-box Gauss rule.
fn point_quadrature(
    responses: &[AbsorbingCrystalResponse],
    region: AbsorbingCrystalBox,
) -> (f64, [f64; 6]) {
    let n = 16;
    let origin = region.origin_m();
    let step = region.extent_m().map(|v| v / n as f64);
    let side = OrientedCrystalFieldSide::Transmitted;
    let mut heat = 0.0;
    for x in 0..=n {
        for y in 0..=n {
            for z in 0..=n {
                let point = [
                    origin[0] + x as f64 * step[0],
                    origin[1] + y as f64 * step[1],
                    origin[2] + z as f64 * step[2],
                ];
                heat += simpson_weight(x, n) * simpson_weight(y, n) * simpson_weight(z, n)
                    * coherent_illumination_field(responses.iter(), side, point, || false)
                        .unwrap().absorbed_power_density_w_m3;
            }
        }
    }
    heat *= step[0] * step[1] * step[2] / 27.0;
    let mut faces = [0.0; 6];
    for axis in 0..3 {
        let a = (axis + 1) % 3;
        let b = (axis + 2) % 3;
        for upper in 0..2 {
            let mut sum = 0.0;
            for i in 0..=n {
                for j in 0..=n {
                    let mut point = origin;
                    point[axis] += upper as f64 * region.extent_m()[axis];
                    point[a] += i as f64 * step[a];
                    point[b] += j as f64 * step[b];
                    sum += simpson_weight(i, n) * simpson_weight(j, n)
                        * coherent_illumination_poynting(responses.iter(), side, point, || false)
                            .unwrap()[axis];
                }
            }
            faces[2 * axis + upper] = sum * step[a] * step[b] / 9.0
                * if upper == 0 { -1.0 } else { 1.0 };
        }
    }
    (heat, faces)
}

#[test]
fn rotated_anisotropic_box_matches_local_work_six_faces_and_xyz_partition() {
    let data = input(
        [0.3, -0.4, 1.0],
        [C64::new(1.5, -0.09), C64::new(1.9, -0.23)],
        0.45,
    );
    let responses = [
        (0.45, [1.0, 0.0]), (-0.3, [0.6, 0.8]), (0.7, [-0.8, 0.6]),
    ].map(|(p, plane)| {
        let mut data = data;
        data.incidence = Incidence::TangentialIndex(p);
        solve_in_plane(data, plane).unwrap()
    });
    let origin = [0.11, -0.27, 0.17].map(|v| v * LAMBDA);
    let extent = [0.23, 0.31, 0.43].map(|v| v * LAMBDA);
    let region = AbsorbingCrystalBox::try_new(origin, extent).unwrap();
    let result = integrate(&responses, region);
    let (heat, faces) = point_quadrature(&responses, region);
    let scale = result.outward_face_power_w.iter().map(|v| v.abs()).sum::<f64>();
    close(result.absorbed_power_w, heat, 3e-6 * scale);
    for (actual, expected) in result.outward_face_power_w.into_iter().zip(faces) {
        close(actual, expected, 3e-6 * scale);
    }
    close(result.residual_w, 0.0, 3e-12 * scale);
    let mut partition_heat = Dd::ZERO;
    let mut outer_faces = [Dd::ZERO; 6];
    for x in 0..2 {
        for y in 0..2 {
            for z in 0..2 {
                let index = [x, y, z];
                let half = extent.map(|v| v * 0.5);
                let start = core::array::from_fn(|axis| origin[axis] + index[axis] as f64 * half[axis]);
                let part = integrate(&responses, AbsorbingCrystalBox::try_new(start, half).unwrap());
                partition_heat = partition_heat + Dd::from_f64(part.absorbed_power_w);
                for axis in 0..3 {
                    let face = 2 * axis + index[axis];
                    outer_faces[face] = outer_faces[face] + Dd::from_f64(part.outward_face_power_w[face]);
                }
            }
        }
    }
    close(partition_heat.to_f64(), result.absorbed_power_w, 5e-12 * scale);
    for (sum, expected) in outer_faces.into_iter().zip(result.outward_face_power_w) {
        close(sum.to_f64(), expected, 5e-12 * scale);
    }
}

#[test]
fn one_beam_and_large_transverse_origins_retain_original_slab_absorption() {
    let response = solve(input(
        [0.3, 0.7, -0.4],
        [C64::new(1.5, -0.09), C64::new(1.9, -0.23)],
        0.0,
    )).unwrap();
    let extent = [2e-6, 3e-6, 0.5 * LAMBDA];
    let origin = [1e100, -1e100, 0.1 * LAMBDA];
    assert_eq!(origin[0] + extent[0], origin[0]);
    let region = AbsorbingCrystalBox::try_new(origin, extent).unwrap();
    let result = integrate(&[response], region);
    let original = response.slab_absorbed_power_w(
        extent[0] * extent[1], origin[2], origin[2] + extent[2],
    ).unwrap();
    close(result.absorbed_power_w, original, 3e-12 * original);
    close(result.residual_w, 0.0, 5e-12 * original);
}

#[test]
fn exact_coherent_darkness_lossless_limits_and_volume_work_admission() {
    let mut data = input([0.0, 0.0, 1.0], [C64::new(1.6, -0.12); 2], 0.6);
    data.electric = [C64::ZERO, C64::ONE];
    let first = solve(data).unwrap();
    data.electric[1] = -C64::ONE;
    let responses = [first, solve(data).unwrap()];
    let region = AbsorbingCrystalBox::try_new(
        [0.1 * LAMBDA, -0.2 * LAMBDA, 0.0],
        [2.0 * LAMBDA, LAMBDA, 30.0 * LAMBDA],
    ).unwrap();
    let dark = integrate(&responses, region);
    assert_eq!(dark.absorbed_power_w, 0.0);
    assert_eq!(dark.outward_face_power_w, [0.0; 6]);
    assert_eq!(dark.residual_w, 0.0);
    let empty: [AbsorbingCrystalResponse; 0] = [];
    assert_eq!(integrate(&empty, region).absorbed_power_w, 0.0);
    for invalid in [
        AbsorbingCrystalBox::try_new([0.0, 0.0, -1.0], [1.0; 3]),
        AbsorbingCrystalBox::try_new([0.0; 3], [1.0, 0.0, 1.0]),
        AbsorbingCrystalBox::try_new([f64::NAN, 0.0, 0.0], [1.0; 3]),
    ] {
        assert!(invalid.is_err());
    }
    let mut limited = limits();
    limited.maximum_mode_pairs = 15;
    assert!(coherent_box_absorption(responses.iter(), region, limited, 0, || false).is_err());
    limited = limits();
    limited.maximum_workspace_bytes = dark.workspace_bytes - 1;
    assert!(coherent_box_absorption(responses.iter(), region, limited, 0, || false).is_err());
    let mut polls = 0;
    let repeat = coherent_box_absorption([first].iter(), region, limits(), 0, || {
        polls += 1;
        false
    }).unwrap();
    for cancel_at in [1, 5, polls] {
        let mut seen = 0;
        let result = coherent_box_absorption([first].iter(), region, limits(), 0, || {
            seen += 1;
            seen == cancel_at
        });
        assert!(matches!(result, Err(IlluminationError::Cancelled)));
        assert_eq!(seen, cancel_at);
    }
    let again = integrate(&[first], region);
    assert_eq!(repeat.absorbed_power_w, again.absorbed_power_w);
    assert_eq!(repeat.outward_face_power_w, again.outward_face_power_w);
    let lossless = solve(input(
        [0.3, 0.7, -0.4], [C64::from_re(1.5), C64::from_re(1.9)], 0.6,
    )).unwrap();
    let result = integrate(&[lossless], region);
    assert_eq!(result.absorbed_power_w, 0.0);
    let scale = result.outward_face_power_w.iter().map(|v| v.abs()).sum::<f64>();
    close(result.residual_w, 0.0, 5e-12 * scale);
}
