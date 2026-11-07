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
    let epsilon = data
        .indices
        .map(|n| (n * n).scale(1.0 / (data.mut_ * C0 * C0)));
    let mut te = [C64::ZERO; 3];
    let mut th = [C64::ZERO; 3];
    for mode in response.transmitted {
        let k = [
            C64::from_re(response.tangential_wavenumber_rad_m),
            C64::ZERO,
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
