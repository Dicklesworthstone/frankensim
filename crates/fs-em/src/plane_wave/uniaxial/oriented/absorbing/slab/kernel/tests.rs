use super::*;
use crate::plane_wave::uniaxial::oriented::{Incidence, slab::SlabRegion};
use fs_math::det;

const MU: f64 = 1.256_637_062_12e-6;
fn input(axis: [f64; 3], indices: [C64; 2], p: f64) -> Input {
    Input {
        omega: C0,
        ni: 1.2,
        mui: MU,
        mut_: MU,
        axis,
        indices,
        incidence: Incidence::TangentialIndex(p),
        electric: [C64::new(0.7, 0.2), C64::new(-0.3, 0.5)],
        loss_f_m: indices.map(|n| -2.0 * n.re * n.im / (MU * C0 * C0)),
    }
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; tolerance {tolerance:e}"
    );
}
fn complex_near(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:?} != {b:?}; tolerance {tolerance:e}"
    );
}
fn phase(q: C64, depth: f64) -> C64 {
    C64::new(det::cos(q.re * depth), -det::sin(q.re * depth)).scale(det::exp(q.im * depth))
}
// Independent scalar infinite Airy sum, in tangential electric admittances.
fn airy(yi: C64, yc: C64, yt: C64, phase: C64) -> (C64, C64) {
    let ri = (yi - yc) / (yi + yc);
    let rt = (yc - yt) / (yc + yt);
    let denominator = C64::ONE + ri * rt * phase * phase;
    (
        (ri + rt * phase * phase) / denominator,
        yi.scale(2.0) / (yi + yc) * yc.scale(2.0) / (yc + yt) * phase / denominator,
    )
}

#[test]
fn absorbing_isotropic_oblique_film_matches_independent_airy_with_magnetic_contrast() {
    for n in [
        C64::new(1.7, -0.15),
        C64::new(0.4, -2.3),
        C64::new(1.7, -1e-8),
    ] {
        for p in [0.0, 0.7, 1.15] {
            for mui in [0.7, 1.0, 1.4] {
                let mut d = input([1.0, -2.0, 3.0], [n; 2], p);
                d.mut_ = MU * mui;
                d.loss_f_m = d.loss_f_m.map(|v| v / mui);
                let (nt, mt, depth) = (1.4, MU * 1.1, 0.83);
                let r = solve(d, nt, mt, depth).unwrap();
                let qi = det::sqrt(d.ni * d.ni - p * p);
                let qt = det::sqrt(nt * nt - p * p);
                let mut qc = (n * n - C64::from_re(p * p)).sqrt();
                if qc.im > 0.0 {
                    qc = -qc;
                }
                let mut expected_r = 0.0;
                let mut expected_t = 0.0;
                let mut expected_i = 0.0;
                for j in 0..2 {
                    let (yi, yc, yt) = if j == 0 {
                        (
                            d.ni * d.ni / (d.mui * qi),
                            n * n / qc.scale(d.mut_),
                            nt * nt / (mt * qt),
                        )
                    } else {
                        (qi / d.mui, qc.scale(1.0 / d.mut_), qt / mt)
                    };
                    let (rr, tt) = airy(C64::from_re(yi), yc, C64::from_re(yt), phase(qc, depth));
                    complex_near(r.reflected_electric_v_m[j], rr * d.electric[j], 3e-13);
                    complex_near(r.transmitted_electric_v_m[j], tt * d.electric[j], 3e-13);
                    expected_i += yi * d.electric[j].norm_sq();
                    expected_r += yi * (rr * d.electric[j]).norm_sq();
                    expected_t += yt * (tt * d.electric[j]).norm_sq();
                }
                near(r.reflected_power_fraction, expected_r / expected_i, 4e-13);
                near(r.transmitted_power_fraction, expected_t / expected_i, 4e-13);
                assert!(r.absorbed_power_fraction > 0.0);
            }
        }
    }
}

#[test]
fn rotated_dichroic_film_retains_two_complex_etalon_polarizations() {
    let axis = [0.6, 0.8, 0.0];
    let d = input(axis, [C64::new(1.4, -0.03), C64::new(2.1, -0.35)], 0.0);
    let depth = 1.7;
    let r = solve(d, 1.0, MU, depth).unwrap();
    let principal = [
        d.electric[0].scale(-0.8) + d.electric[1].scale(0.6),
        d.electric[0].scale(0.6) + d.electric[1].scale(0.8),
    ];
    let transmitted = core::array::from_fn::<_, 2, _>(|j| {
        let (_, t) = airy(
            C64::from_re(d.ni),
            d.indices[j],
            C64::ONE,
            phase(d.indices[j], depth),
        );
        t * principal[j]
    });
    complex_near(
        r.transmitted_electric_v_m[0],
        transmitted[0].scale(-0.8) + transmitted[1].scale(0.6),
        3e-14,
    );
    complex_near(
        r.transmitted_electric_v_m[1],
        transmitted[0].scale(0.6) + transmitted[1].scale(0.8),
        3e-14,
    );
    assert!(
        r.backward
            .iter()
            .any(|w| w.electric_v_m.iter().any(|v| v.abs() > 0.01))
    );
}

#[test]
fn tilted_absorbing_film_heat_matches_local_and_integrated_poynting_theorem() {
    let d = input(
        [0.4, 0.7, 0.6],
        [C64::new(1.5, -0.08), C64::new(2.0, -0.23)],
        0.85,
    );
    let depth = 2.3;
    let r = solve(d, 1.0, MU, depth).unwrap();
    assert!(
        (r.forward[1].normal_wavenumber_rad_m + r.backward[1].normal_wavenumber_rad_m).abs() > 0.1
    );
    let heat = |z| {
        r.field_at(SlabRegion::Crystal, [0.0, 0.0, z])
            .unwrap()
            .absorbed_power_density_w_m3
    };
    let power = |z| {
        r.field_at(SlabRegion::Crystal, [0.0, 0.0, z])
            .unwrap()
            .field
            .poynting_w_m2[2]
    };
    for z in [0.2, 0.9, 1.7, 2.1] {
        let h = 1e-4;
        let derivative = (power(z - 2.0 * h) - 8.0 * power(z - h) + 8.0 * power(z + h)
            - power(z + 2.0 * h))
            / (12.0 * h);
        assert!(heat(z) > 0.0);
        near(-derivative, heat(z), 2e-12);
    }
    let count = 400usize;
    let step = depth / (count as f64);
    let mut integral = heat(0.0) + heat(depth);
    for j in 1..count {
        integral += (if j % 2 == 0 { 2.0 } else { 4.0 }) * heat(j as f64 * step);
    }
    integral *= step / 3.0;
    near(
        integral,
        r.incident_normal_power_w_m2 * r.absorbed_power_fraction,
        2e-13,
    );
    near(power(0.0) - power(depth), integral, 2e-13);
}

#[test]
fn opaque_crystal_has_finite_halfspace_limit_and_no_growing_exponential() {
    let d = input(
        [0.4, 0.7, 0.6],
        [C64::new(1.5, -0.1), C64::new(2.0, -0.3)],
        0.7,
    );
    let r = solve(d, 1.0, MU, 1e5).unwrap();
    let halfspace = crate::plane_wave::uniaxial::oriented::absorbing::kernel::solve(d).unwrap();
    near(
        r.reflected_power_fraction,
        halfspace.reflected_power_fraction,
        3e-14,
    );
    near(
        r.absorbed_power_fraction,
        halfspace.transmitted_power_fraction,
        3e-14,
    );
    assert_eq!(r.transmitted_power_fraction, 0.0);
    for z in [0.0, 1.0, 1e4, 1e5] {
        let f = r.field_at(SlabRegion::Crystal, [0.0, 0.0, z]).unwrap();
        assert!(f.absorbed_power_density_w_m3.is_finite());
        assert!(f.absorbed_power_density_w_m3 >= 0.0);
    }
}

#[test]
fn exact_zero_extinction_has_no_heat_and_continuous_lossless_limit() {
    for p in [0.0, 0.7, 1.05] {
        let d = input([0.3, 0.8, 0.5], [C64::from_re(1.5), C64::from_re(1.9)], p);
        let r = solve(d, 1.0, MU, 0.9).unwrap();
        let small_loss = input(d.axis, [C64::new(1.5, -1e-11), C64::new(1.9, -1e-11)], p);
        let limiting = solve(small_loss, 1.0, MU, 0.9).unwrap();
        near(
            r.reflected_power_fraction + r.transmitted_power_fraction,
            1.0,
            3e-14,
        );
        near(
            r.reflected_power_fraction,
            limiting.reflected_power_fraction,
            2e-10,
        );
        near(
            r.transmitted_power_fraction,
            limiting.transmitted_power_fraction,
            2e-10,
        );
        assert_eq!(r.absorbed_power_fraction, 0.0);
        for j in 0..3 {
            complex_near(
                r.transmitted_electric_v_m[j],
                limiting.transmitted_electric_v_m[j],
                2e-10,
            );
        }
        assert_eq!(
            r.field_at(SlabRegion::Crystal, [0.0, 0.0, 0.4])
                .unwrap()
                .absorbed_power_density_w_m3,
            0.0
        );
    }
}

#[test]
fn rotated_field_and_tensor_heat_are_covariant_and_regions_are_enforced() {
    let d = input(
        [0.3, 0.8, 0.5],
        [C64::new(1.5, -0.07), C64::new(1.9, -0.2)],
        0.7,
    );
    let r = solve(d, 1.0, MU, 0.9).unwrap();
    let plane = [0.6, 0.8];
    let rotate = |v: Vector| {
        [
            v[0].scale(0.6) - v[1].scale(0.8),
            v[0].scale(0.8) + v[1].scale(0.6),
            v[2],
        ]
    };
    let mut turned = r;
    let a = unit(d.axis).unwrap();
    turned
        .rotate(
            plane,
            [0.6 * a[0] - 0.8 * a[1], 0.8 * a[0] + 0.6 * a[1], a[2]],
        )
        .unwrap();
    let f = r.field_at(SlabRegion::Crystal, [0.2, 0.0, 0.4]).unwrap();
    let g = turned
        .field_at(SlabRegion::Crystal, [0.12, 0.16, 0.4])
        .unwrap();
    for j in 0..3 {
        complex_near(
            g.field.electric_v_m[j],
            rotate(f.field.electric_v_m)[j],
            5e-14,
        );
    }
    near(
        g.absorbed_power_density_w_m3,
        f.absorbed_power_density_w_m3,
        2e-17,
    );
    assert!(r.field_at(SlabRegion::Crystal, [0.0, 0.0, -0.1]).is_err());
    assert!(
        r.field_at(SlabRegion::Transmitted, [0.0, 0.0, 0.4])
            .is_err()
    );
    assert!(
        r.field_at(SlabRegion::TotalIncidentSide, [0.0, 0.0, 0.4])
            .is_err()
    );
}

#[test]
fn positive_grazing_cosine_preserves_incident_and_matched_exit_ports() {
    let n = C64::new(1.7, -0.1);
    let mut d = input([1.0, 2.0, 3.0], [n; 2], 0.0);
    let cosine = 1e-9;
    d.incidence = Incidence::Cosine(cosine);
    d.electric = [C64::ZERO, C64::ONE];
    let qi = d.ni * cosine;
    let p = d.ni * det::sqrt((1.0 - cosine) * (1.0 + cosine));
    let mut qc = (n * n - C64::from_re(p * p)).sqrt();
    if qc.im > 0.0 {
        qc = -qc;
    }
    for nt in [d.ni, 1.4] {
        let r = solve(d, nt, MU, 0.9).unwrap();
        let qt = if nt == d.ni {
            qi
        } else {
            det::sqrt(nt * nt - p * p)
        };
        let (rr, tt) = airy(C64::from_re(qi), qc, C64::from_re(qt), phase(qc, 0.9));
        near(r.incident_normal_wavenumber_rad_m, qi, 1e-23);
        if nt == d.ni {
            near(r.exit_normal_wavenumber_rad_m.re, qi, 1e-23);
            assert_eq!(r.exit_normal_wavenumber_rad_m.im, 0.0);
        }
        complex_near(r.reflected_electric_v_m[1], rr, 3e-14);
        complex_near(r.transmitted_electric_v_m[1], tt, 3e-14);
        near(r.reflected_power_fraction, rr.norm_sq(), 3e-14);
        near(
            r.transmitted_power_fraction,
            (qt / qi) * tt.norm_sq(),
            3e-14,
        );
        assert!(r.transmitted_power_fraction > 0.0);
    }
}

#[test]
fn two_face_isotropic_drive_matches_independent_bidirectional_airy_scattering() {
    let right = [C64::new(0.25, -0.4), C64::new(0.5, 0.1)];
    for n in [C64::new(1.7, -0.15), C64::new(0.4, -2.3)] {
        for p in [0.0, 0.7, 1.15] {
            let mut d = input([1.0, -2.0, 3.0], [n; 2], p);
            d.mut_ *= 0.7;
            d.loss_f_m = d.loss_f_m.map(|v| v / 0.7);
            let (nt, mt, depth) = (1.4, MU * 1.1, 0.83);
            let r = solve_from_both_sides(d, nt, mt, depth, right).unwrap();
            let qi = det::sqrt(d.ni * d.ni - p * p);
            let qt = det::sqrt(nt * nt - p * p);
            let mut qc = (n * n - C64::from_re(p * p)).sqrt();
            if qc.im > 0.0 {
                qc = -qc;
            }
            let mut left_input = 0.0;
            let mut right_input = 0.0;
            let mut left_output = 0.0;
            let mut right_output = 0.0;
            for j in 0..2 {
                let (yi, yc, yt) = if j == 0 {
                    (
                        d.ni * d.ni / (d.mui * qi),
                        n * n / qc.scale(d.mut_),
                        nt * nt / (mt * qt),
                    )
                } else {
                    (qi / d.mui, qc.scale(1.0 / d.mut_), qt / mt)
                };
                let (rl, tl) = airy(C64::from_re(yi), yc, C64::from_re(yt), phase(qc, depth));
                let (rr, tr) = airy(C64::from_re(yt), yc, C64::from_re(yi), phase(qc, depth));
                let left = rl * d.electric[j] + tr * right[j];
                let exit = tl * d.electric[j] + rr * right[j];
                complex_near(r.reflected_electric_v_m[j], left, 4e-13);
                complex_near(r.transmitted_electric_v_m[j], exit, 4e-13);
                complex_near(r.right_incident_electric_v_m[j], right[j], 4e-15);
                left_input += yi * d.electric[j].norm_sq();
                right_input += yt * right[j].norm_sq();
                left_output += yi * left.norm_sq();
                right_output += yt * exit.norm_sq();
            }
            near(r.incident_normal_power_w_m2, left_input / (2.0 * C0), 4e-16);
            near(
                r.right_incident_normal_power_w_m2,
                right_input / (2.0 * C0),
                4e-16,
            );
            let incoming = left_input + right_input;
            near(
                r.total_incident_normal_power_w_m2(),
                incoming / (2.0 * C0),
                8e-16,
            );
            near(r.reflected_power_fraction, left_output / incoming, 4e-13);
            near(r.transmitted_power_fraction, right_output / incoming, 4e-13);
            near(
                r.absorbed_power_fraction,
                1.0 - (left_output + right_output) / incoming,
                4e-13,
            );
        }
    }
}

#[test]
fn either_port_can_be_dark_with_unchanged_one_face_and_lossless_power_limits() {
    let d = input([1.0, 0.0, 0.0], [C64::from_re(1.5), C64::from_re(1.9)], 0.7);
    let depth = 0.9;
    let left = solve(d, d.ni, d.mui, depth).unwrap();
    assert_eq!(
        left,
        solve_from_both_sides(d, d.ni, d.mui, depth, [C64::ZERO; 2]).unwrap()
    );
    let dark = Input {
        electric: [C64::ZERO; 2],
        ..d
    };
    let right = solve_from_both_sides(dark, d.ni, d.mui, depth, d.electric).unwrap();
    assert_eq!(right.incident_normal_power_w_m2, 0.0);
    near(
        right.right_incident_normal_power_w_m2,
        left.incident_normal_power_w_m2,
        1e-17,
    );
    near(
        right.reflected_power_fraction,
        left.transmitted_power_fraction,
        3e-14,
    );
    near(
        right.transmitted_power_fraction,
        left.reflected_power_fraction,
        3e-14,
    );
    near(
        right.reflected_power_fraction + right.transmitted_power_fraction,
        1.0,
        3e-14,
    );
    assert_eq!(right.absorbed_power_fraction, 0.0);
    for j in 0..3 {
        complex_near(
            right.reflected_electric_v_m[j],
            left.transmitted_electric_v_m[j].scale(if j == 2 { -1.0 } else { 1.0 }),
            4e-14,
        );
        complex_near(
            right.transmitted_electric_v_m[j],
            left.reflected_electric_v_m[j].scale(if j == 2 { -1.0 } else { 1.0 }),
            4e-14,
        );
    }
    assert!(solve_from_both_sides(dark, d.ni, d.mui, depth, [C64::ZERO; 2]).is_err());
    for p in [1.0, 1.05] {
        let grazing_exit = Input {
            incidence: Incidence::TangentialIndex(p),
            ..d
        };
        assert!(solve_from_both_sides(grazing_exit, 1.0, MU, depth, d.electric).is_err());
    }
    let evanescent = solve(
        Input {
            incidence: Incidence::TangentialIndex(1.05),
            ..d
        },
        1.0,
        MU,
        depth,
    )
    .unwrap();
    let empty = evanescent
        .exit_incident_field_at([0.0, 0.0, depth + 1e4])
        .unwrap();
    assert_eq!(empty.field.electric_v_m, [C64::ZERO; 3]);
    assert_eq!(empty.field.magnetic_a_m, [C64::ZERO; 3]);
    assert_eq!(empty.absorbed_power_density_w_m3, 0.0);
    assert!(
        right
            .exit_incident_field_at([0.0, 0.0, depth - 0.1])
            .is_err()
    );
    assert!(right.total_exit_field_at([f64::NAN, 0.0, depth]).is_err());
}

#[test]
fn relative_port_phase_controls_absorption_and_retains_a_standing_wave_heat_node() {
    let n = C64::new(1.7, -0.15);
    let mut d = input([1.0, 0.0, 0.0], [n; 2], 0.0);
    d.ni = 1.0;
    d.electric = [C64::ZERO, C64::ONE];
    let depth = 0.8;
    let same = solve_from_both_sides(d, 1.0, MU, depth, d.electric).unwrap();
    let opposite = solve_from_both_sides(d, 1.0, MU, depth, d.electric.map(|v| -v)).unwrap();
    let (rr, tt) = airy(C64::ONE, n, C64::ONE, phase(n, depth));
    near(
        same.absorbed_power_fraction,
        1.0 - (rr + tt).norm_sq(),
        4e-14,
    );
    near(
        opposite.absorbed_power_fraction,
        1.0 - (rr - tt).norm_sq(),
        4e-14,
    );
    assert!((same.absorbed_power_fraction - opposite.absorbed_power_fraction).abs() > 0.1);
    for response in [&same, &opposite] {
        let heat = super::super::absorption::integrated_heat(response, 1.0, 0.0, depth).unwrap();
        near(
            heat,
            response.total_incident_normal_power_w_m2() * response.absorbed_power_fraction,
            3e-16,
        );
        let top = response
            .field_at(SlabRegion::Crystal, [0.0, 0.0, 0.0])
            .unwrap();
        let bottom = response
            .field_at(SlabRegion::Crystal, [0.0, 0.0, depth])
            .unwrap();
        assert!(top.field.poynting_w_m2[2] > 0.0);
        assert!(bottom.field.poynting_w_m2[2] < 0.0);
        near(
            heat,
            top.field.poynting_w_m2[2] - bottom.field.poynting_w_m2[2],
            3e-16,
        );
    }
    let middle = depth * 0.5;
    let heat_at = |z| {
        opposite
            .field_at(SlabRegion::Crystal, [0.0, 0.0, z])
            .unwrap()
            .absorbed_power_density_w_m3
    };
    assert!(heat_at(middle) < 1e-25);
    let h = 1e-5;
    let node_heat =
        super::super::absorption::integrated_heat(&opposite, 1.0, middle - h, middle + h).unwrap();
    let simpson = h / 3.0 * (heat_at(middle - h) + 4.0 * heat_at(middle) + heat_at(middle + h));
    assert!(node_heat > 0.0);
    near(node_heat, simpson, 2e-6 * simpson);
}

#[test]
fn tilted_two_face_fields_match_both_faces_and_signed_poynting_after_azimuth_rotation() {
    let d = input(
        [0.4, 0.7, 0.6],
        [C64::new(1.5, -0.08), C64::new(2.0, -0.23)],
        0.85,
    );
    let depth = 2.3;
    let right = [C64::new(-0.2, 0.4), C64::new(0.6, -0.1)];
    let response = solve_from_both_sides(d, 1.0, MU * 0.85, depth, right).unwrap();
    let mut turned = response;
    let rotate = |v: Vector| {
        [
            v[0].scale(0.6) - v[1].scale(0.8),
            v[0].scale(0.8) + v[1].scale(0.6),
            v[2],
        ]
    };
    let a = unit(d.axis).unwrap();
    turned
        .rotate(
            [0.6, 0.8],
            [0.6 * a[0] - 0.8 * a[1], 0.8 * a[0] + 0.6 * a[1], a[2]],
        )
        .unwrap();
    for r in [&response, &turned] {
        for z in [0.0, depth] {
            let point = [0.12, -0.16, z];
            let inside = r.field_at(SlabRegion::Crystal, point).unwrap();
            let outside = if z == 0.0 {
                r.field_at(SlabRegion::TotalIncidentSide, point).unwrap()
            } else {
                r.total_exit_field_at(point).unwrap()
            };
            for j in 0..2 {
                complex_near(
                    inside.field.electric_v_m[j],
                    outside.field.electric_v_m[j],
                    4e-13,
                );
                complex_near(
                    inside.field.magnetic_a_m[j],
                    outside.field.magnetic_a_m[j],
                    4e-15,
                );
            }
            near(
                inside.field.poynting_w_m2[2],
                outside.field.poynting_w_m2[2],
                4e-16,
            );
            assert_eq!(outside.absorbed_power_density_w_m3, 0.0);
        }
        let heat = super::super::absorption::integrated_heat(r, 0.7, 0.0, depth).unwrap();
        near(
            heat,
            0.7 * r.total_incident_normal_power_w_m2() * r.absorbed_power_fraction,
            5e-16,
        );
    }
    let original = response
        .total_exit_field_at([0.2, 0.0, depth + 0.4])
        .unwrap();
    let rotated = turned
        .total_exit_field_at([0.12, 0.16, depth + 0.4])
        .unwrap();
    for j in 0..3 {
        complex_near(
            rotated.field.electric_v_m[j],
            rotate(original.field.electric_v_m)[j],
            4e-13,
        );
        complex_near(
            rotated.field.magnetic_a_m[j],
            rotate(original.field.magnetic_a_m)[j],
            4e-15,
        );
    }
    let incoming = response.exit_incident_field_at([0.0, 0.0, depth]).unwrap();
    near(
        incoming.field.poynting_w_m2[2],
        -response.right_incident_normal_power_w_m2,
        3e-17,
    );
    let total = response.total_exit_field_at([0.0, 0.0, depth]).unwrap();
    let outgoing = response
        .field_at(SlabRegion::Transmitted, [0.0, 0.0, depth])
        .unwrap();
    near(
        total.field.poynting_w_m2[2],
        outgoing.field.poynting_w_m2[2] + incoming.field.poynting_w_m2[2],
        3e-17,
    );
}
