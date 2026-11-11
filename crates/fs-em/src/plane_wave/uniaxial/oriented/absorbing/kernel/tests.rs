use super::*;
use fs_math::det;

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

fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:.17e} != {b:.17e}; tolerance {tolerance:e}"
    );
}
fn complex_close(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:?} != {b:?}; tolerance {tolerance:e}"
    );
}

#[test]
fn isotropic_complex_limit_matches_independent_fresnel_with_magnetic_contrast() {
    for n in [
        C64::new(1.6, -0.15),
        C64::new(0.3, -3.2),
        C64::new(0.8, -1e-8),
    ] {
        for p in [0.0, 0.6, 0.99] {
            let mut d = input([0.3, -0.2, 0.7], [n; 2], p);
            d.mut_ = MU * 1.4;
            d.loss_f_m = d.loss_f_m.map(|v| v / 1.4);
            let r = solve(d).unwrap();
            let qi = det::sqrt(d.ni * d.ni - p * p);
            let qt = (n * n - C64::from_re(p * p)).sqrt();
            for j in 0..2 {
                let yi = if j == 0 {
                    d.ni * d.ni / (d.mui * qi)
                } else {
                    qi / d.mui
                };
                let yt = if j == 0 {
                    (n * n) / qt.scale(d.mut_)
                } else {
                    qt.scale(1.0 / d.mut_)
                };
                let reflection = (C64::from_re(yi) - yt) / (C64::from_re(yi) + yt);
                complex_close(
                    r.reflected_electric_v_m[j],
                    reflection * d.electric[j],
                    2e-14,
                );
                let sum = r.transmitted[0].electric_v_m[j] + r.transmitted[1].electric_v_m[j];
                complex_close(sum, (C64::ONE + reflection) * d.electric[j], 2e-14);
            }
            close(
                r.reflected_power_fraction + r.transmitted_power_fraction,
                1.0,
                2e-14,
            );
            for m in r.transmitted {
                complex_close(m.normal_wavenumber_rad_m, qt.scale(d.omega / C0), 3e-8);
                assert!(m.normal_wavenumber_rad_m.im < 0.0);
            }
        }
    }
}

#[test]
fn principal_axis_dichroism_has_original_complex_phases_and_analytic_heat() {
    let indices = [C64::new(1.5, -0.08), C64::new(1.9, -0.3)];
    let d = input([1.0, 0.0, 0.0], indices, 0.0);
    let r = solve(d).unwrap();
    for z in [0.0, 0.2 * LAMBDA, 2.0 * LAMBDA] {
        let f = r
            .field(OrientedCrystalFieldSide::Transmitted, [0.0, 0.0, z])
            .unwrap();
        let mut expected_heat = 0.0;
        for (j, n) in [indices[1], indices[0]].into_iter().enumerate() {
            let transmission = C64::from_re(2.0) / (C64::ONE + n);
            let theta = -d.omega / C0 * n.re * z;
            let phase =
                C64::new(det::cos(theta), det::sin(theta)).scale(det::exp(d.omega / C0 * n.im * z));
            let field = transmission * d.electric[j] * phase;
            complex_close(f.field.electric_v_m[j], field, 3e-15);
            expected_heat += d.omega * (-n.re * n.im) / (MU * C0 * C0) * field.norm_sq();
        }
        close(
            f.absorbed_power_density_w_m3,
            expected_heat,
            2e-13 * expected_heat,
        );
    }
    // Axis normal to the face admits a different extraordinary phase index.
    let p = 0.7;
    let d = input([0.0, 0.0, 1.0], indices, p);
    let r = solve(d).unwrap();
    let qe = indices[0] * (C64::ONE - C64::from_re(p * p) / (indices[1] * indices[1])).sqrt();
    complex_close(
        r.transmitted[1].normal_wavenumber_rad_m,
        qe.scale(d.omega / C0),
        3e-8,
    );
}

fn check_maxwell(data: Input, response: &AbsorbingCrystalResponse) {
    let axis = transparent::unit(data.axis).unwrap();
    let plane = response.incidence_plane_axis();
    let epsilon = data
        .indices
        .map(|n| (n * n).scale(1.0 / (data.mut_ * C0 * C0)));
    let mut te = [C64::ZERO; 3];
    let mut th = [C64::ZERO; 3];
    for mode in response.transmitted {
        let k = [
            C64::from_re(response.tangential_wavenumber_rad_m * plane[0]),
            C64::from_re(response.tangential_wavenumber_rad_m * plane[1]),
            mode.normal_wavenumber_rad_m,
        ];
        let dot = (0..3).fold(C64::ZERO, |s, j| s + mode.electric_v_m[j].scale(axis[j]));
        let displacement: Vector = core::array::from_fn(|j| {
            mode.electric_v_m[j] * epsilon[0] + dot * (epsilon[1] - epsilon[0]).scale(axis[j])
        });
        let faraday = cross(k, mode.electric_v_m);
        let ampere = cross(k, mode.magnetic_a_m);
        for j in 0..3 {
            complex_close(
                faraday[j],
                mode.magnetic_a_m[j].scale(data.omega * data.mut_),
                5e-8,
            );
            complex_close(ampere[j], displacement[j].scale(-data.omega), 1e-9);
            te[j] = te[j] + mode.electric_v_m[j];
            th[j] = th[j] + mode.magnetic_a_m[j];
        }
    }
    for j in 0..2 {
        complex_close(
            response.incident_electric_v_m[j] + response.reflected_electric_v_m[j],
            te[j],
            3e-14,
        );
        complex_close(
            response.incident_magnetic_a_m[j] + response.reflected_magnetic_a_m[j],
            th[j],
            3e-16,
        );
    }
    close(
        flux(te, th) / response.incident_normal_power_w_m2,
        response.transmitted_power_fraction,
        2e-14,
    );
}

#[test]
fn arbitrary_azimuth_retains_global_maxwell_tensor_and_authored_spatial_phase() {
    let d = input(
        [1.0, 2.0, -0.6],
        [C64::new(1.4, -0.12), C64::new(2.1, -0.42)],
        0.65,
    );
    for authored_plane in [[3.0, 4.0], [-2.0, 1.0], [0.0, -7.0]] {
        let r = solve_in_plane(d, authored_plane).unwrap();
        check_maxwell(d, &r);
        for j in 0..2 {
            complex_close(r.incident_electric_v_m[j], d.electric[j], 4e-16);
        }
        let plane = r.incidence_plane_axis();
        let point = [0.13 * LAMBDA, -0.41 * LAMBDA, 0.23 * LAMBDA];
        let mut expected_e = [C64::ZERO; 3];
        let mut expected_h = [C64::ZERO; 3];
        for mode in r.transmitted {
            let angle = r.tangential_wavenumber_rad_m * (plane[0] * point[0] + plane[1] * point[1])
                + mode.normal_wavenumber_rad_m.re * point[2];
            let phase = C64::new(det::cos(angle), -det::sin(angle))
                .scale(det::exp(mode.normal_wavenumber_rad_m.im * point[2]));
            for j in 0..3 {
                expected_e[j] = expected_e[j] + mode.electric_v_m[j] * phase;
                expected_h[j] = expected_h[j] + mode.magnetic_a_m[j] * phase;
            }
        }
        let actual = r
            .field(OrientedCrystalFieldSide::Transmitted, point)
            .unwrap();
        for j in 0..3 {
            complex_close(actual.field.electric_v_m[j], expected_e[j], 2e-15);
            complex_close(actual.field.magnetic_a_m[j], expected_h[j], 2e-17);
        }
        let at_face = [point[0], point[1], 0.0];
        let left = r
            .field(OrientedCrystalFieldSide::TotalIncidentSide, at_face)
            .unwrap();
        let right = r
            .field(OrientedCrystalFieldSide::Transmitted, at_face)
            .unwrap();
        for j in 0..2 {
            complex_close(
                left.field.electric_v_m[j],
                right.field.electric_v_m[j],
                3e-15,
            );
            complex_close(
                left.field.magnetic_a_m[j],
                right.field.magnetic_a_m[j],
                3e-17,
            );
        }
        close(
            left.field.poynting_w_m2[2],
            right.field.poynting_w_m2[2],
            3e-17,
        );
    }
}

#[test]
fn azimuth_rotational_covariance_preserves_original_heat_and_volume() {
    let d = input(
        [1.0, 2.0, -0.6],
        [C64::new(1.4, -0.12), C64::new(2.1, -0.42)],
        -0.65,
    );
    let reference = solve(d).unwrap();
    let plane = [0.6, 0.8];
    let rotate = |v: Vector| {
        [
            v[0].scale(plane[0]) - v[1].scale(plane[1]),
            v[0].scale(plane[1]) + v[1].scale(plane[0]),
            v[2],
        ]
    };
    let [ax, ay, az] = d.axis;
    let e = rotate([d.electric[0], d.electric[1], C64::ZERO]);
    let mut rotated = d;
    rotated.axis = [
        ax * plane[0] - ay * plane[1],
        ax * plane[1] + ay * plane[0],
        az,
    ];
    rotated.electric = [e[0], e[1]];
    let response = solve_in_plane(rotated, plane).unwrap();
    check_maxwell(rotated, &response);
    let original_point = [0.17 * LAMBDA, -0.11 * LAMBDA, 0.41 * LAMBDA];
    let rotated_point = [
        original_point[0] * plane[0] - original_point[1] * plane[1],
        original_point[0] * plane[1] + original_point[1] * plane[0],
        original_point[2],
    ];
    let side = OrientedCrystalFieldSide::Transmitted;
    let old = reference.field(side, original_point).unwrap();
    let actual = response.field(side, rotated_point).unwrap();
    let expected_e = rotate(old.field.electric_v_m);
    let expected_h = rotate(old.field.magnetic_a_m);
    let expected_flux = rotate(old.field.poynting_w_m2.map(C64::from_re));
    for j in 0..3 {
        complex_close(actual.field.electric_v_m[j], expected_e[j], 4e-15);
        complex_close(actual.field.magnetic_a_m[j], expected_h[j], 4e-17);
        close(actual.field.poynting_w_m2[j], expected_flux[j].re, 4e-17);
    }
    close(
        actual.absorbed_power_density_w_m3,
        old.absorbed_power_density_w_m3,
        2e-13 * old.absorbed_power_density_w_m3,
    );
    let area = 0.3;
    let (start, end) = (0.13 * LAMBDA, 2.7 * LAMBDA);
    let watts = response.slab_absorbed_power_w(area, start, end).unwrap();
    close(
        watts,
        reference.slab_absorbed_power_w(area, start, end).unwrap(),
        3e-13 * watts,
    );
    let front = response
        .electromagnetic_field(side, [0.0, 0.0, start])
        .unwrap();
    let back = response
        .electromagnetic_field(side, [0.0, 0.0, end])
        .unwrap();
    close(
        watts,
        area * (front.poynting_w_m2[2] - back.poynting_w_m2[2]),
        3e-13 * watts,
    );
    rotated.incidence = Incidence::TangentialIndex(0.65);
    let equivalent = solve_in_plane(rotated, plane.map(|v| -v)).unwrap();
    let equivalent_field = equivalent.field(side, rotated_point).unwrap();
    for j in 0..3 {
        complex_close(
            equivalent_field.field.electric_v_m[j],
            actual.field.electric_v_m[j],
            4e-15,
        );
        complex_close(
            equivalent_field.field.magnetic_a_m[j],
            actual.field.magnetic_a_m[j],
            4e-17,
        );
    }
    close(
        equivalent.slab_absorbed_power_w(area, start, end).unwrap(),
        watts,
        3e-13 * watts,
    );
}

#[test]
fn azimuth_admission_and_normal_incidence_keep_legacy_and_transparent_limits() {
    for indices in [
        [C64::new(1.4, -0.12), C64::new(2.1, -0.42)],
        [C64::from_re(1.4), C64::from_re(2.1)],
    ] {
        let d = input([1.0, 2.0, -0.6], indices, 0.0);
        let legacy = solve(d).unwrap();
        assert_eq!(solve_in_plane(d, [1.0, 0.0]).unwrap(), legacy);
        assert_eq!(solve_in_plane(d, [f64::MAX, 0.0]).unwrap(), legacy);
        let side = OrientedCrystalFieldSide::Transmitted;
        let point = [0.2 * LAMBDA, 0.8 * LAMBDA, 0.3 * LAMBDA];
        let expected = legacy.field(side, point).unwrap();
        for plane in [[3.0, 4.0], [0.0, -f64::MAX], [1e-300, -2e-300]] {
            let response = solve_in_plane(d, plane).unwrap();
            check_maxwell(d, &response);
            let actual = response.field(side, point).unwrap();
            for j in 0..3 {
                complex_close(
                    actual.field.electric_v_m[j],
                    expected.field.electric_v_m[j],
                    4e-15,
                );
                complex_close(
                    actual.field.magnetic_a_m[j],
                    expected.field.magnetic_a_m[j],
                    4e-17,
                );
            }
            close(
                actual.absorbed_power_density_w_m3,
                expected.absorbed_power_density_w_m3,
                2e-13 * expected.absorbed_power_density_w_m3,
            );
        }
        for plane in [[0.0, 0.0], [f64::NAN, 1.0], [1.0, f64::INFINITY]] {
            assert!(solve_in_plane(d, plane).is_err());
        }
    }
}

#[test]
fn rotated_absorption_retains_maxwell_coupling_and_coherent_cross_mode_power() {
    let mut d = input(
        [1.0, 2.0, -0.6],
        [C64::new(1.4, -0.12), C64::new(2.1, -0.42)],
        0.8,
    );
    for axis in [[1.0, 2.0, -0.6], [-0.2, 0.4, 0.8], [0.0, 1.0, 0.0]] {
        d.axis = axis;
        let r = solve(d).unwrap();
        check_maxwell(d, &r);
        assert!(
            r.transmitted
                .iter()
                .all(|m| m.normal_wavenumber_rad_m.im < 0.0 && m.normal_power_w_m2 > 0.0)
        );
        let mut reversed = d;
        reversed.axis = axis.map(|v| -v);
        let flipped = solve(reversed).unwrap();
        for j in 0..3 {
            complex_close(
                r.reflected_electric_v_m[j],
                flipped.reflected_electric_v_m[j],
                1e-14,
            );
        }
    }
    d.axis = [1.0, 2.0, -0.6];
    d.electric = [C64::ZERO, C64::ONE];
    let r = solve(d).unwrap();
    assert!(r.reflected_electric_v_m[0].abs() > 0.01);
    let point = [0.1 * LAMBDA, 0.0, 0.2 * LAMBDA];
    let total = r
        .field(OrientedCrystalFieldSide::Transmitted, point)
        .unwrap();
    let ordinary = r.field(OrientedCrystalFieldSide::Ordinary, point).unwrap();
    let extra = r
        .field(OrientedCrystalFieldSide::Extraordinary, point)
        .unwrap();
    assert!(
        (total.field.poynting_w_m2[2]
            - ordinary.field.poynting_w_m2[2]
            - extra.field.poynting_w_m2[2])
            .abs()
            > 1e-7
    );
    assert!(
        (total.absorbed_power_density_w_m3
            - ordinary.absorbed_power_density_w_m3
            - extra.absorbed_power_density_w_m3)
            .abs()
            > 1.0
    );
}

#[test]
fn local_anisotropic_heating_matches_poynting_loss_and_finite_depth_integral() {
    let d = input(
        [0.2, 0.7, -0.4],
        [C64::new(1.5, -0.09), C64::new(1.9, -0.23)],
        0.65,
    );
    let r = solve(d).unwrap();
    let at = |z| {
        r.field(OrientedCrystalFieldSide::Transmitted, [0.0, 0.0, z])
            .unwrap()
    };
    let step = LAMBDA * 1e-5;
    for z in [0.1 * LAMBDA, 0.7 * LAMBDA, 2.0 * LAMBDA] {
        let expected = -(at(z + step).field.poynting_w_m2[2] - at(z - step).field.poynting_w_m2[2])
            / (2.0 * step);
        let actual = at(z).absorbed_power_density_w_m3;
        assert!(actual > 0.0);
        close(actual, expected, 3e-8 * actual);
    }
    let depth = 6.0 * LAMBDA;
    let samples = 512;
    let dz = depth / samples as f64;
    let integral = (0..=samples)
        .map(|j| {
            let weight = if j == 0 || j == samples {
                1.0
            } else if j % 2 == 0 {
                2.0
            } else {
                4.0
            };
            weight * at(j as f64 * dz).absorbed_power_density_w_m3
        })
        .sum::<f64>()
        * dz
        / 3.0;
    let flux_loss = at(0.0).field.poynting_w_m2[2] - at(depth).field.poynting_w_m2[2];
    close(integral, flux_loss, 3e-8 * flux_loss);
    let incident = r
        .field(
            OrientedCrystalFieldSide::TotalIncidentSide,
            [0.0, 0.0, -LAMBDA],
        )
        .unwrap();
    assert_eq!(incident.absorbed_power_density_w_m3, 0.0);
}

#[test]
fn original_zero_loss_retains_transparent_fields_and_evanescence_is_not_heat() {
    let d = input(
        [0.3, 0.4, -0.2],
        [C64::from_re(1.5), C64::from_re(1.9)],
        0.7,
    );
    let r = solve(d).unwrap();
    let reference = transparent::solve(transparent::Input {
        omega: d.omega,
        ni: d.ni,
        no: 1.5,
        ne: 1.9,
        mui: d.mui,
        mut_: d.mut_,
        axis: d.axis,
        incidence: d.incidence,
        electric: d.electric,
    })
    .unwrap();
    assert_eq!(r.reflected_electric_v_m, reference.reflected_electric_v_m);
    assert_eq!(
        r.reflected_power_fraction,
        reference.reflected_power_fraction
    );
    assert_eq!(
        r.transmitted_power_fraction,
        reference.transmitted_power_fraction
    );
    assert_eq!(
        r.field(OrientedCrystalFieldSide::Transmitted, [0.0, 0.0, LAMBDA])
            .unwrap()
            .absorbed_power_density_w_m3,
        0.0
    );
    let mut d = input(
        [1.0, 0.0, 0.0],
        [C64::from_re(0.6), C64::new(1.9, -0.2)],
        0.8,
    );
    d.electric = [C64::ZERO, C64::ONE];
    let r = solve(d).unwrap();
    assert!(r.transmitted[0].normal_wavenumber_rad_m.im < 0.0);
    assert_eq!(r.transmitted[1].electric_v_m, [C64::ZERO; 3]);
    let field = r
        .field(OrientedCrystalFieldSide::Transmitted, [0.0, 0.0, LAMBDA])
        .unwrap();
    assert_eq!(field.absorbed_power_density_w_m3, 0.0);
    close(field.field.poynting_w_m2[2], 0.0, 1e-20);
}

#[test]
fn passive_source_and_physical_half_space_admission_are_explicit() {
    let d = input(
        [1.0, 1.0, 1.0],
        [C64::new(1.4, -0.1), C64::new(1.9, -0.2)],
        0.5,
    );
    let r = solve(d).unwrap();
    assert!(
        r.field(OrientedCrystalFieldSide::Transmitted, [0.0, 0.0, -1e-6])
            .is_err()
    );
    assert!(
        r.field(OrientedCrystalFieldSide::Incident, [0.0, 0.0, 1e-6])
            .is_err()
    );
    assert!(
        r.field(OrientedCrystalFieldSide::Transmitted, [f64::NAN, 0.0, 1e-6])
            .is_err()
    );
    let mut invalid = d;
    invalid.indices[1].im = 0.1;
    assert!(solve(invalid).is_err());
    invalid = d;
    invalid.electric = [C64::ZERO; 2];
    assert!(solve(invalid).is_err());
    invalid = d;
    invalid.incidence = Incidence::TangentialIndex(1.0);
    assert!(solve(invalid).is_err());
    invalid = d;
    invalid.loss_f_m[0] = 0.0;
    assert!(solve(invalid).is_err());
}

#[test]
fn analytic_crystal_volume_heat_matches_independent_depth_quadrature_and_flux() {
    let response = solve(input(
        [0.2, 0.7, -0.4],
        [C64::new(1.5, -0.09), C64::new(1.9, -0.23)],
        0.65,
    ))
    .unwrap();
    let area = 2.5e-8;
    let at = |z| {
        response
            .field(OrientedCrystalFieldSide::Transmitted, [0.0, 0.0, z])
            .unwrap()
    };
    for [start, end] in [[0.0, 0.01], [0.2, 1.7], [0.1, 20.0]].map(|v| v.map(|z| z * LAMBDA)) {
        let actual = response.slab_absorbed_power_w(area, start, end).unwrap();
        let samples = 4096;
        let dz = (end - start) / samples as f64;
        let quadrature = (0..=samples)
            .map(|j| {
                let weight = if j == 0 || j == samples {
                    1.0
                } else if j % 2 == 0 {
                    2.0
                } else {
                    4.0
                };
                weight * at(start + j as f64 * dz).absorbed_power_density_w_m3
            })
            .sum::<f64>()
            * area
            * dz
            / 3.0;
        close(actual, quadrature, 3e-10 * actual);
        let flux_loss = area * (at(start).field.poynting_w_m2[2] - at(end).field.poynting_w_m2[2]);
        close(actual, flux_loss, 3e-12 * actual);
    }
}

#[test]
fn integrated_heat_retains_ordinary_extraordinary_interference_and_interval_additivity() {
    let mut data = input(
        [1.0, 2.0, -0.6],
        [C64::new(1.4, -0.12), C64::new(2.1, -0.42)],
        0.8,
    );
    data.electric = [C64::ZERO, C64::ONE];
    let response = solve(data).unwrap();
    let [start, split, end] = [0.1, 0.6, 2.0].map(|z| z * LAMBDA);
    let area = 3e-9;
    let total = response.slab_absorbed_power_w(area, start, end).unwrap();
    let partial = response.slab_absorbed_power_w(area, start, split).unwrap()
        + response.slab_absorbed_power_w(area, split, end).unwrap();
    close(total, partial, 3e-14 * total);
    close(
        response
            .slab_absorbed_power_w(area * 7.0, start, end)
            .unwrap(),
        total * 7.0,
        3e-14 * total,
    );
    let modal_sum = (0..2)
        .map(|m| {
            let mut single = response;
            single.transmitted[m].electric_v_m = [C64::ZERO; 3];
            single.slab_absorbed_power_w(area, start, end).unwrap()
        })
        .sum::<f64>();
    assert!((total - modal_sum).abs() > 1e-3 * total);
}

#[test]
fn isotropic_slab_heat_matches_fresnel_beer_lambert_for_both_polarizations() {
    for extinction in [0.0, 1e-20, 0.2, 3.0] {
        let n = C64::new(1.6, -extinction);
        let response = solve(input([1.0, 2.0, 3.0], [n; 2], 0.0)).unwrap();
        let area = 4e-7;
        let start = 0.3 * LAMBDA;
        let width = 1.2 * LAMBDA;
        let transmitted = C64::from_re(2.0) / (C64::ONE + n);
        let entering = response
            .incident_electric_v_m
            .iter()
            .map(|&electric| (electric * transmitted).norm_sq())
            .sum::<f64>()
            * n.re
            / (2.0 * MU * C0);
        let alpha = response.omega_rad_s() / C0 * extinction;
        let expected =
            area * entering * (-2.0 * alpha * start).exp() * -(-2.0 * alpha * width).exp_m1();
        let actual = response
            .slab_absorbed_power_w(area, start, start + width)
            .unwrap();
        close(actual, expected, 8e-14 * expected);
    }
}

#[test]
fn weak_thin_opaque_and_deep_crystal_volumes_keep_physical_heat_limits() {
    let data = input(
        [0.2, 0.7, -0.4],
        [C64::new(1.5, -0.09), C64::new(1.9, -0.23)],
        0.65,
    );
    let response = solve(data).unwrap();
    let area = 7e-8;
    let entrance = response
        .field(OrientedCrystalFieldSide::Transmitted, [0.0; 3])
        .unwrap();
    let thin = LAMBDA * 1e-14;
    close(
        response.slab_absorbed_power_w(area, 0.0, thin).unwrap(),
        entrance.absorbed_power_density_w_m3 * area * thin,
        3e-13 * entrance.absorbed_power_density_w_m3 * area * thin,
    );
    let opaque = LAMBDA * 1e6;
    close(
        response.slab_absorbed_power_w(area, 0.0, opaque).unwrap(),
        entrance.field.poynting_w_m2[2] * area,
        3e-14 * entrance.field.poynting_w_m2[2] * area,
    );
    assert_eq!(
        response
            .slab_absorbed_power_w(area, opaque, opaque * 2.0)
            .unwrap(),
        0.0
    );
    let weak = solve(input([1.0, 0.0, 0.0], [C64::new(1.5, -1e-20); 2], 0.0)).unwrap();
    assert!(weak.slab_absorbed_power_w(area, 0.0, LAMBDA).unwrap() > 0.0);
    let zero_loss = solve(input(
        [0.2, 0.7, -0.4],
        data.indices.map(|n| C64::from_re(n.re)),
        0.65,
    ))
    .unwrap();
    assert_eq!(
        zero_loss.slab_absorbed_power_w(area, 0.0, opaque).unwrap(),
        0.0
    );
    let mut evanescent = input(
        [1.0, 0.0, 0.0],
        [C64::from_re(0.6), C64::new(1.9, -0.2)],
        0.8,
    );
    evanescent.electric = [C64::ZERO, C64::ONE];
    assert_eq!(
        solve(evanescent)
            .unwrap()
            .slab_absorbed_power_w(area, 0.0, opaque)
            .unwrap(),
        0.0
    );
}

#[test]
fn one_factor_statistics_preserve_complex_polarization_and_peak_rms_convention() {
    let mut data = input(
        [0.2, 0.7, -0.4],
        [C64::new(1.5, -0.09), C64::new(1.9, -0.23)],
        0.0,
    );
    data.electric = [C64::ONE, C64::new(0.0, 1.0)];
    let response = solve(data).unwrap();
    let incident = illumination_field_statistics(
        [&response],
        OrientedCrystalFieldSide::Incident,
        [0.0; 3],
        || false,
    )
    .unwrap();
    complex_close(
        incident.electric_coherency_v2_m2[0][1],
        C64::new(0.0, -0.5),
        2e-15,
    );
    close(incident.electric_coherency_v2_m2[0][0].re, 0.5, 2e-15);
    close(incident.electric_coherency_v2_m2[1][1].re, 0.5, 2e-15);
    for (side, depth) in [
        (OrientedCrystalFieldSide::TotalIncidentSide, -0.3 * LAMBDA),
        (OrientedCrystalFieldSide::Transmitted, 0.7 * LAMBDA),
    ] {
        let point = [0.2 * LAMBDA, 0.0, depth];
        let expected = response.field(side, point).unwrap();
        let actual = illumination_field_statistics([&response], side, point, || false).unwrap();
        assert_eq!(actual.poynting_w_m2, expected.field.poynting_w_m2);
        assert_eq!(
            actual.absorbed_power_density_w_m3,
            expected.absorbed_power_density_w_m3
        );
        for i in 0..3 {
            for j in 0..3 {
                complex_close(
                    actual.electric_coherency_v2_m2[i][j],
                    expected.field.electric_v_m[i].scale(0.5)
                        * expected.field.electric_v_m[j].conj(),
                    3e-15,
                );
                assert_eq!(
                    actual.electric_coherency_v2_m2[i][j],
                    actual.electric_coherency_v2_m2[j][i].conj()
                );
            }
        }
    }
    let area = 2e-8;
    let end = 1.7 * LAMBDA;
    assert_eq!(
        illumination_absorption([&response], area, 0.0, end, || false)
            .unwrap()
            .1,
        response.slab_absorbed_power_w(area, 0.0, end).unwrap()
    );
}

#[test]
fn unpolarized_oblique_factors_are_basis_invariant_and_do_not_interfere() {
    let p: f64 = 0.8;
    let c = (1.0 - p * p).sqrt();
    let data = input(
        [1.0, 2.0, -0.6],
        [C64::new(1.4, -0.12), C64::new(2.1, -0.42)],
        p,
    );
    let prepare = |electric| {
        let mut d = data;
        d.electric = electric;
        solve(d).unwrap()
    };
    let factors = [
        prepare([C64::from_re(c), C64::ZERO]),
        prepare([C64::ZERO, C64::ONE]),
    ];
    let h = core::f64::consts::FRAC_1_SQRT_2;
    let rotated = [
        prepare([C64::from_re(c * h), C64::from_re(h)]),
        prepare([C64::from_re(-c * h), C64::from_re(h)]),
    ];
    let incident = illumination_field_statistics(
        factors.iter(),
        OrientedCrystalFieldSide::Incident,
        [0.0; 3],
        || false,
    )
    .unwrap();
    close(
        (0..3)
            .map(|j| incident.electric_coherency_v2_m2[j][j].re)
            .sum(),
        1.0,
        3e-15,
    );
    close(incident.poynting_w_m2[2], c / (MU * C0), 3e-15);
    let point = [0.0, 0.0, 0.7 * LAMBDA];
    let stats = |bank: &[AbsorbingCrystalResponse]| {
        illumination_field_statistics(
            bank.iter(),
            OrientedCrystalFieldSide::Transmitted,
            point,
            || false,
        )
        .unwrap()
    };
    let a = stats(&factors);
    let b = stats(&rotated);
    for i in 0..3 {
        for j in 0..3 {
            complex_close(
                a.electric_coherency_v2_m2[i][j],
                b.electric_coherency_v2_m2[i][j],
                3e-14,
            );
        }
        close(a.poynting_w_m2[i], b.poynting_w_m2[i], 3e-16);
    }
    close(
        a.absorbed_power_density_w_m3,
        b.absorbed_power_density_w_m3,
        3e-13 * a.absorbed_power_density_w_m3,
    );
    let heat = |bank: &[AbsorbingCrystalResponse]| {
        illumination_absorption(bank.iter(), 2e-8, 0.1 * LAMBDA, 1.7 * LAMBDA, || false)
            .unwrap()
            .1
    };
    close(heat(&factors), heat(&rotated), 3e-13 * heat(&factors));
    let coherent = [prepare([C64::from_re(c), C64::ONE])];
    assert!((heat(&coherent) - heat(&factors)).abs() > 1e-3 * heat(&factors));
    let phased = [factors[0], prepare([C64::ZERO, C64::new(0.0, 1.0)])];
    close(heat(&phased), heat(&factors), 3e-14 * heat(&factors));
}

#[test]
fn independent_frequencies_and_directions_keep_original_loss_and_volume_balance() {
    let mut first = input(
        [0.2, 0.7, -0.4],
        [C64::new(1.5, -0.04), C64::new(1.9, -0.18)],
        0.3,
    );
    first.omega *= 0.7;
    let mut second = input(
        first.axis,
        [C64::new(1.8, -0.2), C64::new(2.3, -0.33)],
        -0.7,
    );
    second.omega *= 1.4;
    second.electric = [C64::new(0.4, -0.2), C64::new(-0.1, 0.6)];
    let factors = [solve(first).unwrap(), solve(second).unwrap()];
    let side = OrientedCrystalFieldSide::Transmitted;
    let point = [0.2 * LAMBDA, 0.0, 0.4 * LAMBDA];
    let statistics = illumination_field_statistics(factors.iter(), side, point, || false).unwrap();
    let expected_heat = factors
        .iter()
        .map(|r| r.field(side, point).unwrap().absorbed_power_density_w_m3)
        .sum::<f64>();
    close(
        statistics.absorbed_power_density_w_m3,
        expected_heat,
        3e-15 * expected_heat,
    );
    assert_eq!(
        illumination_poynting(factors.iter(), side, point, || false).unwrap(),
        statistics.poynting_w_m2
    );
    let area = 7e-8;
    let start = 0.1 * LAMBDA;
    let split = 0.7 * LAMBDA;
    let end = 2.3 * LAMBDA;
    let absorbed = illumination_absorption(factors.iter(), area, start, end, || false)
        .unwrap()
        .1;
    let split_heat = illumination_absorption(factors.iter(), area, start, split, || false)
        .unwrap()
        .1
        + illumination_absorption(factors.iter(), area, split, end, || false)
            .unwrap()
            .1;
    close(absorbed, split_heat, 3e-14 * absorbed);
    let flux = |z| illumination_poynting(factors.iter(), side, [0.0, 0.0, z], || false).unwrap()[2];
    close(absorbed, area * (flux(start) - flux(end)), 3e-13 * absorbed);
    let power = illumination_power(factors.iter(), || false).unwrap();
    close(
        power.incident_w_m2,
        power.reflected_w_m2 + power.transmitted_w_m2,
        3e-14 * power.incident_w_m2,
    );
}

#[test]
fn flux_observation_does_not_require_local_heat_density_to_be_representable() {
    let mut d = input([0.0, 1.0, 0.0], [C64::new(1.5, -0.2); 2], 0.0);
    d.electric = [C64::from_re(1e154), C64::ZERO];
    let response = solve(d).unwrap();
    let side = OrientedCrystalFieldSide::Transmitted;
    // The same original source has finite E/H and energy flux, while local
    // dielectric power density exceeds f64. A small physical
    // volume still has finite integrated watts and needs a usable flux trace.
    assert!(response.field(side, [0.0; 3]).is_err());
    assert!(illumination_field_statistics([&response], side, [0.0; 3], || false).is_err());
    let front = illumination_poynting([&response], side, [0.0; 3], || false).unwrap();
    let back = illumination_poynting([&response], side, [0.0, 0.0, 1e-6], || false).unwrap();
    assert!(front.iter().all(|v| v.is_finite()));
    assert!(front[2] > back[2] && back[2] > 0.0);
    let (_, heat) = illumination_absorption([&response], 1e-6, 0.0, 1e-6, || false).unwrap();
    close(heat, 1e-6 * (front[2] - back[2]), 3e-13 * heat);
}

#[test]
fn dark_statistics_still_admit_geometry_and_variable_work_has_typed_cancellation() {
    let empty: [AbsorbingCrystalResponse; 0] = [];
    let side = OrientedCrystalFieldSide::Transmitted;
    let dark = illumination_field_statistics(empty.iter(), side, [0.0; 3], || false).unwrap();
    assert_eq!(dark.electric_coherency_v2_m2, [[C64::ZERO; 3]; 3]);
    assert_eq!(dark.poynting_w_m2, [0.0; 3]);
    assert_eq!(dark.absorbed_power_density_w_m3, 0.0);
    assert_eq!(
        illumination_absorption(empty.iter(), 2.0, 1.0, 4.0, || false).unwrap(),
        (6.0, 0.0)
    );
    assert!(illumination_field_statistics(empty.iter(), side, [0.0, 0.0, -1.0], || false).is_err());
    assert!(illumination_poynting(empty.iter(), side, [f64::NAN, 0.0, 0.0], || false).is_err());
    assert!(illumination_absorption(empty.iter(), 0.0, 0.0, 1.0, || false).is_err());
    assert!(matches!(
        illumination_power(empty.iter(), || true),
        Err(IlluminationError::Cancelled)
    ));
    let response = solve(input(
        [1.0, 2.0, -0.6],
        [C64::new(1.4, -0.12), C64::new(2.1, -0.42)],
        0.8,
    ))
    .unwrap();
    let mut polls = 0;
    let result = illumination_absorption([&response, &response], 2e-8, 0.0, LAMBDA, || {
        polls += 1;
        polls == 4
    });
    assert!(matches!(result, Err(IlluminationError::Cancelled)));
    assert_eq!(polls, 4);
}

#[test]
fn coherent_beams_interfere_before_loss_and_preserve_weak_residuals() {
    let mut data = input([0.0, 0.0, 1.0], [C64::new(1.6, -0.12); 2], 0.6);
    data.electric = [C64::ZERO, C64::new(0.8, 0.3)];
    let wave = solve(data).unwrap();
    let mut reversed = data;
    reversed.electric = data.electric.map(|value| -value);
    let opposite = solve(reversed).unwrap();
    use OrientedCrystalFieldSide::*;
    for side in [
        Incident,
        Reflected,
        TotalIncidentSide,
        Ordinary,
        Extraordinary,
        Transmitted,
    ] {
        let z = if matches!(side, Ordinary | Extraordinary | Transmitted) {
            0.27 * LAMBDA
        } else {
            -0.27 * LAMBDA
        };
        let point = [0.13 * LAMBDA, -0.2 * LAMBDA, z];
        let single = wave.field(side, point).unwrap();
        let doubled =
            coherent_illumination_field([&wave, &wave], side, point, || false).unwrap();
        for j in 0..3 {
            complex_close(
                doubled.field.electric_v_m[j],
                single.field.electric_v_m[j].scale(2.0),
                5e-15,
            );
            complex_close(
                doubled.field.magnetic_a_m[j],
                single.field.magnetic_a_m[j].scale(2.0),
                5e-17,
            );
            close(
                doubled.field.poynting_w_m2[j],
                4.0 * single.field.poynting_w_m2[j],
                5e-17,
            );
        }
        close(
            doubled.absorbed_power_density_w_m3,
            4.0 * single.absorbed_power_density_w_m3,
            5e-14 * single.absorbed_power_density_w_m3.max(1.0),
        );
        let dark =
            coherent_illumination_field([&wave, &opposite], side, point, || false).unwrap();
        for value in dark.field.electric_v_m {
            assert!(value.abs() < 1e-24);
        }
        assert!(dark.absorbed_power_density_w_m3 < 1e-40);
    }
    // An unrelated faint coherent component must survive two large opposite
    // fields. Summing beam powers, or rounding an intermediate field sum,
    // loses either the dark port or its physically resolved weak signal.
    let mut strong_data = data;
    strong_data.electric = data.electric.map(|value| value.scale(1e16));
    let strong = solve(strong_data).unwrap();
    strong_data.electric = strong_data.electric.map(|value| -value);
    let strong_opposite = solve(strong_data).unwrap();
    let point = [0.0; 3];
    let expected = wave.field(Transmitted, point).unwrap();
    let residual = coherent_illumination_field(
        [&strong, &wave, &strong_opposite],
        Transmitted,
        point,
        || false,
    )
    .unwrap();
    for j in 0..3 {
        complex_close(
            residual.field.electric_v_m[j],
            expected.field.electric_v_m[j],
            5e-14,
        );
        complex_close(
            residual.field.magnetic_a_m[j],
            expected.field.magnetic_a_m[j],
            5e-16,
        );
    }
    close(
        residual.absorbed_power_density_w_m3,
        expected.absorbed_power_density_w_m3,
        2e-13 * expected.absorbed_power_density_w_m3,
    );
    let independent =
        illumination_field_statistics([&wave, &opposite], Transmitted, point, || false).unwrap();
    close(
        independent.absorbed_power_density_w_m3,
        2.0 * expected.absorbed_power_density_w_m3,
        2e-14 * expected.absorbed_power_density_w_m3,
    );
}

#[test]
fn crossed_absorbing_te_beams_match_fresnel_heat_fringe_and_phase_translation() {
    let n = C64::new(1.7, -0.18);
    let p = 0.6;
    let mut data = input([0.0, 0.0, 1.0], [n; 2], p);
    data.electric = [C64::ZERO, C64::ONE];
    let forward = solve(data).unwrap();
    let k0 = data.omega / C0;
    let qi = det::sqrt(1.0 - p * p);
    let q = (n * n - C64::from_re(p * p)).sqrt();
    let transmission = C64::from_re(2.0 * qi) / (C64::from_re(qi) + q);
    let side = OrientedCrystalFieldSide::Transmitted;
    for relative_phase in [0.0, 0.7, core::f64::consts::PI] {
        let mut other = data;
        other.incidence = Incidence::TangentialIndex(-p);
        other.electric[1] =
            C64::new(det::cos(relative_phase), det::sin(relative_phase));
        let backward = solve(other).unwrap();
        for z in [0.0, 0.2 * LAMBDA, 1.3 * LAMBDA] {
            let propagation = C64::new(det::cos(k0 * q.re * z), -det::sin(k0 * q.re * z))
                .scale(det::exp(k0 * q.im * z));
            for fringe_fraction in [0.0, 0.13, 0.25, 0.41, 0.5, 0.73] {
                let x = fringe_fraction * core::f64::consts::TAU / (k0 * p);
                let envelope = 2.0 * det::cos(k0 * p * x + 0.5 * relative_phase);
                let expected_e = transmission
                    * propagation
                    * C64::new(det::cos(0.5 * relative_phase), det::sin(0.5 * relative_phase))
                        .scale(envelope);
                let intensity = 4.0
                    * transmission.norm_sq()
                    * det::exp(2.0 * k0 * q.im * z)
                    * det::cos(k0 * p * x + 0.5 * relative_phase).powi(2);
                let expected_heat =
                    data.omega * (-n.re * n.im) / (MU * C0 * C0) * intensity;
                let actual = coherent_illumination_field(
                    [&forward, &backward],
                    side,
                    [x, 0.31 * LAMBDA, z],
                    || false,
                )
                .unwrap();
                complex_close(actual.field.electric_v_m[1], expected_e, 8e-15);
                close(
                    actual.absorbed_power_density_w_m3,
                    expected_heat,
                    3e-13 * (data.omega * data.loss_f_m[0]),
                );
                close(
                    actual.field.poynting_w_m2[2],
                    q.re / (2.0 * C0 * MU) * intensity,
                    1e-16,
                );
                close(actual.field.poynting_w_m2[0], 0.0, 1e-17);
            }
        }
        if relative_phase == 0.0 {
            let node = coherent_illumination_field(
                [&forward, &backward],
                side,
                [core::f64::consts::FRAC_PI_2 / (k0 * p), 0.0, 0.27 * LAMBDA],
                || false,
            )
            .unwrap();
            // A magnetic antinode has no dielectric work when E cancels.
            assert!(node.field.magnetic_a_m[2].abs() > 1e-5);
            assert!(node.absorbed_power_density_w_m3 < 1e-24);
        }
    }
}

#[test]
fn coherent_noncollinear_crystal_fields_obey_three_dimensional_energy_balance() {
    let data = input(
        [0.2, 0.7, -0.4],
        [C64::new(1.5, -0.09), C64::new(1.9, -0.23)],
        0.45,
    );
    let responses = [
        (0.45, [1.0, 0.0]),
        (-0.3, [0.6, 0.8]),
        (0.7, [-0.8, 0.6]),
    ]
    .map(|(p, axis)| {
        let mut data = data;
        data.incidence = Incidence::TangentialIndex(p);
        solve_in_plane(data, axis).unwrap()
    });
    let side = OrientedCrystalFieldSide::Transmitted;
    let at = |point| coherent_illumination_field(responses.iter(), side, point, || false).unwrap();
    let step = 1e-5 * LAMBDA;
    for position in [[0.0, 0.0, 0.1], [0.2, -0.3, 0.4], [-0.4, 0.7, 0.9]] {
        let point = position.map(|value| value * LAMBDA);
        let actual = at(point);
        let mut divergence = 0.0;
        for j in 0..3 {
            let mut plus = point;
            let mut minus = point;
            plus[j] += step;
            minus[j] -= step;
            divergence +=
                (at(plus).field.poynting_w_m2[j] - at(minus).field.poynting_w_m2[j])
                    / (2.0 * step);
        }
        close(
            actual.absorbed_power_density_w_m3,
            -divergence,
            5e-8 * actual.absorbed_power_density_w_m3,
        );
        assert_eq!(
            coherent_illumination_poynting(responses.iter(), side, point, || false).unwrap(),
            actual.field.poynting_w_m2
        );
        let face = [point[0], point[1], 0.0];
        let inside = at(face);
        let outside = coherent_illumination_field(
            responses.iter(),
            OrientedCrystalFieldSide::TotalIncidentSide,
            face,
            || false,
        )
        .unwrap();
        for j in 0..2 {
            complex_close(
                outside.field.electric_v_m[j],
                inside.field.electric_v_m[j],
                8e-14,
            );
            complex_close(
                outside.field.magnetic_a_m[j],
                inside.field.magnetic_a_m[j],
                8e-16,
            );
        }
        close(
            outside.field.poynting_w_m2[2],
            inside.field.poynting_w_m2[2],
            5e-16,
        );
        assert_eq!(outside.absorbed_power_density_w_m3, 0.0);
    }
    let point = [0.0, 0.0, 0.1 * LAMBDA];
    let coherent = at(point).absorbed_power_density_w_m3;
    let independent = responses
        .iter()
        .map(|response| response.field(side, point).unwrap().absorbed_power_density_w_m3)
        .sum::<f64>();
    assert!((coherent - independent).abs() > 0.1 * independent);
}

#[test]
fn coherent_dark_fields_and_cancellation_preserve_observation_admission() {
    let empty: [AbsorbingCrystalResponse; 0] = [];
    let side = OrientedCrystalFieldSide::Transmitted;
    let dark = coherent_illumination_field(empty.iter(), side, [0.0; 3], || false).unwrap();
    assert_eq!(dark.field.electric_v_m, [C64::ZERO; 3]);
    assert_eq!(dark.field.magnetic_a_m, [C64::ZERO; 3]);
    assert_eq!(dark.field.poynting_w_m2, [0.0; 3]);
    assert_eq!(dark.absorbed_power_density_w_m3, 0.0);
    assert!(coherent_illumination_field(empty.iter(), side, [0.0, 0.0, -1.0], || false).is_err());
    assert!(coherent_illumination_poynting(empty.iter(), side, [f64::NAN, 0.0, 0.0], || false).is_err());
    let response = solve(input(
        [1.0, 2.0, -0.6],
        [C64::new(1.4, -0.12), C64::new(2.1, -0.42)],
        0.8,
    ))
    .unwrap();
    for cancel_at in 1..=5 {
        let mut polls = 0;
        let result = coherent_illumination_field([&response, &response], side, [0.0; 3], || {
            polls += 1;
            polls == cancel_at
        });
        assert!(matches!(result, Err(IlluminationError::Cancelled)));
        assert_eq!(polls, cancel_at);
    }
    for cancel_at in 1..=4 {
        let mut polls = 0;
        let result = coherent_illumination_poynting([&response, &response], side, [0.0; 3], || {
            polls += 1;
            polls == cancel_at
        });
        assert!(matches!(result, Err(IlluminationError::Cancelled)));
        assert_eq!(polls, cancel_at);
    }
}

#[test]
fn coherent_poynting_does_not_require_representable_dielectric_heat() {
    let mut data = input([0.0, 0.0, 1.0], [C64::new(1.6, -0.12); 2], 0.6);
    data.electric = [C64::ZERO, C64::from_re(1e153)];
    let response = solve(data).unwrap();
    let side = OrientedCrystalFieldSide::Transmitted;
    let flux = coherent_illumination_poynting([&response, &response], side, [0.0; 3], || false)
        .unwrap();
    assert!(flux.iter().all(|value| value.is_finite()));
    assert!(flux[2] > 0.0);
    assert!(
        coherent_illumination_field([&response, &response], side, [0.0; 3], || false).is_err()
    );
}

#[test]
fn g1_the_extraordinary_mode_keeps_the_loss_the_double_precision_root_lost() {
    // The optic axis in the interface plane at normal incidence: the ordinary tensor element does not reach the
    // extraordinary mode. The tip refused both cases as ambiguous, with Im -5.17e-17 and -4.93e-17 on the two
    // roots; the exactly transparent one is a dichroic polarizer.
    for loss in [-1.0e-24, 0.0] {
        let indices = [C64::new(1.45, -0.12), C64::new(2.1, loss)];
        let basis = eigenbasis(input([1.0, 0.0, 0.0], indices, 0.0)).unwrap();
        let z = basis.modes[1].z;
        assert!(
            z.re > 0.0 && (z.re - 2.1 / 2.0).abs() <= f64::EPSILON,
            "{z:?}"
        );
        assert!((z.im - loss / 2.0).abs() <= 1.0e-28, "{z:?}");
        assert!(flux(basis.modes[1].e, basis.modes[1].h) > 0.0);
    }
}

#[test]
fn g1_a_nearly_in_plane_axis_is_treated_like_the_in_plane_one() {
    for tilt in [1.0e-9, 1.0e-12] {
        let indices = [C64::new(1.45, -0.12), C64::new(2.1, -1.0e-24)];
        let basis = eigenbasis(input([1.0, 0.0, tilt], indices, 0.0)).unwrap();
        let z = basis.modes[1].z;
        assert!(z.re > 0.0 && z.im < 0.0, "{tilt:e} {z:?}");
    }
}

#[test]
fn g0_the_certified_roots_match_the_double_precision_roots_to_their_rounding() {
    // The roots of the 27209 print, in scaled units, with the inputs that gave them.
    let cases = [
        (
            [1.0, 0.0, 0.0],
            -1.0e-3,
            0.0,
            C64::new(1.05, -5.000_000_000_000_017e-4),
        ),
        (
            [1.0, 0.0, 0.0],
            -1.0e-24,
            0.3,
            C64::new(1.027_754_329_547_13, -3.748_697_141_589_883e-3),
        ),
        (
            [1.0, 0.0, 0.0],
            -1.0e-24,
            0.7,
            C64::new(0.922_655_119_953_977_6, -2.273_441_819_896_722_4e-2),
        ),
        (
            [0.6, 0.0, 0.8],
            -1.0e-24,
            0.0,
            C64::new(0.806_124_793_378_505_7, -5.251_381_597_692_664e-2),
        ),
        (
            [0.6, 0.0, 0.8],
            -1.0e-24,
            0.3,
            C64::new(0.749_428_317_644_494_1, -6.084_845_472_171_465e-2),
        ),
        (
            [0.6, 0.0, 0.8],
            -1.0e-3,
            0.7,
            C64::new(0.640_333_370_224_624_9, -7.103_034_402_419_4e-2),
        ),
    ];
    for (axis, loss, p, old) in cases {
        let indices = [C64::new(1.45, -0.12), C64::new(2.1, loss)];
        let basis = eigenbasis(input(axis, indices, p)).unwrap();
        let z = basis.modes[1].z;
        assert!(
            (z - old).abs() <= 1.0e-15,
            "{axis:?} {loss:e} {p}: {z:?} against {old:?}"
        );
        assert!(z.im < 0.0, "{axis:?} {loss:e} {p}: {z:?}");
    }
    // With a loss of 1e-3 on the in-plane axis at normal incidence the root is n_e / 2 itself.
    let indices = [C64::new(1.45, -0.12), C64::new(2.1, -1.0e-3)];
    let z = eigenbasis(input([1.0, 0.0, 0.0], indices, 0.0))
        .unwrap()
        .modes[1]
        .z;
    assert_eq!((z.re, z.im), (2.1 / 2.0, -1.0e-3 / 2.0));
}
