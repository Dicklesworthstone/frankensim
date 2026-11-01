//! Independent interior Maxwell and dielectric-loss regressions.
use super::finite_layers::debye;
use super::illumination::{act, matrix, product};
use super::*;
use fs_em::plane_wave::layers::{DielectricLayer, IlluminationBudget, PlanarDielectricLayers};
use fs_math::det;

const OMEGA: f64 = 1e9;
fn k() -> f64 {
    OMEGA * det::sqrt(E) * det::sqrt(M)
}
fn phase(x: f64) -> C64 {
    C64::new(det::cos(x), det::sin(x))
}
fn tangent(e: [C64; 3], h: [C64; 3], pol: Polarization) -> [C64; 2] {
    match pol {
        Polarization::Te => [e[1], -h[0]],
        Polarization::Tm => [e[0], h[1]],
    }
}

#[test]
fn homogeneous_counterpropagating_fields_retain_vector_and_standing_wave_structure() {
    let outside = material(E, M);
    let thickness = 0.7 / k();
    let layers = [DielectricLayer {
        response: DielectricResponse::Lossless(&outside),
        thickness_m: thickness,
    }];
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        OMEGA,
        &layers,
    )
    .unwrap();
    let input = [C64::new(0.8, 0.3), C64::new(-0.2, 0.7)];
    let base = det::sqrt(E / M);
    for c in [0.3, 0.75, 1.0] {
        for pol in [Polarization::Te, Polarization::Tm] {
            let bank = stack.prepare_illumination(c, pol).unwrap();
            let y = bank.port_admittance_s()[0];
            let boundary = bank.evaluate(input, &mut [0.0]).unwrap();
            for fraction in [0.0, 0.23, 0.5, 1.0] {
                let z = thickness * fraction;
                let point = bank.prepare_point(0, z).unwrap();
                assert!(std::ptr::eq(point.illumination(), &bank));
                assert_eq!(point.layer_index(), 0);
                assert_eq!(point.depth_m(), z);
                let field = point.evaluate(input).unwrap();
                let left = input[0] * phase(-k() * c * z);
                let right = input[1] * phase(-k() * c * (thickness - z));
                let expected_e = left + right;
                let expected_h = (left - right).scale(y);
                let actual = tangent(field.electric_v_m, field.magnetic_a_m, pol);
                complex_near(actual[0], expected_e, 3e-14);
                complex_near(actual[1], expected_h, 3e-16);
                let sine = det::sqrt((1.0 - c) * (1.0 + c));
                match pol {
                    Polarization::Te => {
                        complex_near(field.magnetic_a_m[2], expected_e.scale(base * sine), 3e-16)
                    }
                    Polarization::Tm => {
                        complex_near(field.electric_v_m[2], expected_h.scale(-sine / base), 3e-14)
                    }
                }
                near(
                    field.statistics.complex_poynting_w_m2[2].re,
                    0.5 * y * (input[0].norm_sq() - input[1].norm_sq()),
                    3e-16,
                );
                assert_eq!(field.statistics.absorbed_power_w_m3, 0.0);
                if fraction == 0.0 || fraction == 1.0 {
                    let face = usize::from(fraction == 1.0);
                    complex_near(actual[0], boundary.boundary_tangential_e_v_m[face], 3e-14);
                    complex_near(actual[1], boundary.boundary_conjugate_h_a_m[face], 3e-16);
                }
            }
            if c == 1.0 {
                let node = bank
                    .prepare_point(0, 0.5 * thickness)
                    .unwrap()
                    .evaluate([C64::ONE, -C64::ONE])
                    .unwrap();
                near(
                    node.electric_v_m.iter().map(|e| e.norm_sq()).sum(),
                    0.0,
                    1e-28,
                );
                assert!(node.magnetic_a_m.iter().map(|h| h.norm_sq()).sum::<f64>() > y * y);
                near(node.statistics.complex_poynting_w_m2[2].re, 0.0, 1e-18);
            }
        }
    }
}

/// Conventional characteristic matrices independently resolve local fields;
/// quadrature of positive local loss must recover the original integrated heat.
#[test]
fn lossy_local_fields_match_independent_matrices_and_integrated_heat() {
    let left = material(E, M);
    let right = material(4.0 * E, 1.2 * M);
    let passive = debye(3.0 * E, 4.0 * E, 0.03, 0.001);
    let transparent = material(2.0 * E, 0.7 * M);
    let layers = [
        DielectricLayer {
            response: DielectricResponse::Debye(&passive),
            thickness_m: 0.2 / k(),
        },
        DielectricLayer {
            response: DielectricResponse::Lossless(&transparent),
            thickness_m: 0.3 / k(),
        },
    ];
    let stack = PlanarDielectricLayers::try_new(
        &left,
        DielectricResponse::Lossless(&right),
        OMEGA,
        &layers,
    )
    .unwrap();
    for c in [1.0, 0.72] {
        for pol in [Polarization::Te, Polarization::Tm] {
            let bank = stack.prepare_illumination(c, pol).unwrap();
            let y = bank.port_admittance_s();
            let ratio = y[1] / y[0];
            let m0 = matrix(C64::new(5.0, -2.0), 1.0, 0.2, c, pol);
            let m1 = matrix(C64::from_re(2.0), 0.7, 0.3, c, pol);
            let all = product(m0, m1);
            let plus = act(all, [C64::ONE, C64::from_re(ratio)]);
            let minus = act(all, [C64::ONE, C64::from_re(-ratio)]);
            for angle in [0.0, 0.6, 2.1] {
                let input = [C64::new(0.4, 0.7), phase(angle).scale(0.8)];
                let transmitted =
                    (input[0].scale(2.0) - (minus[0] + minus[1]) * input[1]) / (plus[0] + plus[1]);
                let after = [
                    transmitted + input[1],
                    (transmitted - input[1]).scale(ratio),
                ];
                let between = act(m1, after);
                for (i, eps, mu, depth, end) in [
                    (0, C64::new(5.0, -2.0), 1.0, 0.2, between),
                    (1, C64::from_re(2.0), 0.7, 0.3, after),
                ] {
                    for fraction in [0.0, 0.2, 0.5, 0.85, 1.0] {
                        let expected = act(matrix(eps, mu, depth * (1.0 - fraction), c, pol), end);
                        let field = bank
                            .prepare_point(i, layers[i].thickness_m * fraction)
                            .unwrap()
                            .evaluate(input)
                            .unwrap();
                        let actual = tangent(field.electric_v_m, field.magnetic_a_m, pol);
                        complex_near(actual[0], expected[0], 2e-13);
                        complex_near(actual[1], expected[1].scale(y[0]), 2e-15);
                        let expected_heat = OMEGA
                            * (-eps.im * E)
                            * field
                                .statistics
                                .electric_rms_squared_v2_m2
                                .iter()
                                .sum::<f64>();
                        near(
                            field.statistics.absorbed_power_w_m3,
                            expected_heat,
                            2e-14 * expected_heat.max(1.0),
                        );
                    }
                }
                let mut heat = [0.0; 2];
                bank.evaluate(input, &mut heat).unwrap();
                let steps = 256;
                let mut integral = 0.0;
                for j in 0..=steps {
                    let weight = if j == 0 || j == steps {
                        1.0
                    } else if j % 2 == 0 {
                        2.0
                    } else {
                        4.0
                    };
                    let q = bank
                        .prepare_point(0, layers[0].thickness_m * j as f64 / steps as f64)
                        .unwrap()
                        .evaluate(input)
                        .unwrap()
                        .statistics
                        .absorbed_power_w_m3;
                    integral += weight * q;
                }
                integral *= layers[0].thickness_m / (3.0 * steps as f64);
                near(integral, heat[0], 2e-11 * heat[0]);
                assert_eq!(heat[1], 0.0);
                let z = 0.41 * layers[0].thickness_m;
                let dz = 1e-5 * layers[0].thickness_m;
                let flux = |depth| {
                    bank.prepare_point(0, depth)
                        .unwrap()
                        .evaluate(input)
                        .unwrap()
                        .statistics
                        .complex_poynting_w_m2[2]
                        .re
                };
                let local = bank
                    .prepare_point(0, z)
                    .unwrap()
                    .evaluate(input)
                    .unwrap()
                    .statistics
                    .absorbed_power_w_m3;
                near(
                    -(flux(z + dz) - flux(z - dz)) / (2.0 * dz),
                    local,
                    2e-8 * local,
                );
            }
        }
    }
}

#[test]
fn tangential_electric_node_can_still_heat_a_tm_film_through_its_longitudinal_field() {
    let outside = material(E, M);
    let passive = debye(3.0 * E, 4.0 * E, 0.1, 0.001);
    let layers = [DielectricLayer {
        response: DielectricResponse::Debye(&passive),
        thickness_m: 0.1 / k(),
    }];
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        OMEGA,
        &layers,
    )
    .unwrap();
    for pol in [Polarization::Te, Polarization::Tm] {
        let bank = stack.prepare_illumination(0.7, pol).unwrap();
        let field = bank
            .prepare_point(0, 0.5 * layers[0].thickness_m)
            .unwrap()
            .evaluate([C64::ONE, -C64::ONE])
            .unwrap();
        let e = tangent(field.electric_v_m, field.magnetic_a_m, pol)[0];
        assert!(e.abs() < 1e-13);
        match pol {
            Polarization::Te => assert!(field.statistics.absorbed_power_w_m3 < 1e-26),
            Polarization::Tm => {
                assert!(field.electric_v_m[2].abs() > 0.1);
                assert!(field.statistics.absorbed_power_w_m3 > 1e-4);
                near(
                    field.statistics.absorbed_power_w_m3,
                    0.5 * OMEGA * 2.0 * E * field.electric_v_m[2].norm_sq(),
                    1e-16,
                );
            }
        }
    }
}

#[test]
fn exactly_critical_and_evanescent_films_keep_finite_interior_fields() {
    let outside = material(E, M);
    let critical = material(0.75 * E, M);
    for depth in [0.4, 1e100] {
        let layers = [DielectricLayer {
            response: DielectricResponse::Lossless(&critical),
            thickness_m: depth / k(),
        }];
        let stack = PlanarDielectricLayers::try_new(
            &outside,
            DielectricResponse::Lossless(&outside),
            OMEGA,
            &layers,
        )
        .unwrap();
        for pol in [Polarization::Te, Polarization::Tm] {
            let bank = stack.prepare_illumination(0.5, pol).unwrap();
            for fraction in [0.0, 0.25, 0.5, 1.0] {
                let field = bank
                    .prepare_point(0, layers[0].thickness_m * fraction)
                    .unwrap()
                    .evaluate([C64::ONE, C64::ZERO])
                    .unwrap();
                let actual = tangent(field.electric_v_m, field.magnetic_a_m, pol);
                let coupling = match pol {
                    Polarization::Te => 0.5,
                    Polarization::Tm => 0.375,
                };
                let denominator = C64::new(2.0, coupling * depth);
                let constant = C64::from_re(2.0) / denominator;
                let varying = C64::ONE
                    + C64::new(0.0, coupling * depth * (1.0 - 2.0 * fraction)) / denominator;
                let expected = match pol {
                    Polarization::Te => [varying, constant],
                    Polarization::Tm => [constant, varying],
                };
                complex_near(actual[0], expected[0], 2e-13);
                complex_near(
                    actual[1],
                    expected[1].scale(bank.port_admittance_s()[0]),
                    2e-15,
                );
                assert_eq!(field.statistics.absorbed_power_w_m3, 0.0);
            }
        }
    }
    let low = material(0.16 * E, M);
    let layers = [DielectricLayer {
        response: DielectricResponse::Lossless(&low),
        thickness_m: 0.8 / k(),
    }];
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        OMEGA,
        &layers,
    )
    .unwrap();
    for pol in [Polarization::Te, Polarization::Tm] {
        let bank = stack.prepare_illumination(0.6, pol).unwrap();
        let input = [C64::new(0.5, 0.2), C64::new(-0.1, 0.7)];
        let boundary = bank.evaluate(input, &mut [0.0]).unwrap();
        let end = [
            boundary.boundary_tangential_e_v_m[1],
            boundary.boundary_conjugate_h_a_m[1].scale(1.0 / bank.port_admittance_s()[0]),
        ];
        let reference = act(matrix(C64::from_re(0.16), 1.0, 0.8 * 0.65, 0.6, pol), end);
        let field = bank
            .prepare_point(0, 0.35 * layers[0].thickness_m)
            .unwrap()
            .evaluate(input)
            .unwrap();
        let actual = tangent(field.electric_v_m, field.magnetic_a_m, pol);
        complex_near(actual[0], reference[0], 2e-13);
        complex_near(
            actual[1],
            reference[1].scale(bank.port_admittance_s()[0]),
            2e-15,
        );
    }
}

#[test]
fn opaque_film_points_retain_the_nearest_face_phase_even_for_enormous_thickness() {
    let outside = material(E, M);
    let passive = debye(3.0 * E, 4.0 * E, 0.1, 0.001);
    let layers = [DielectricLayer {
        response: DielectricResponse::Debye(&passive),
        thickness_m: 1e200 / k(),
    }];
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        OMEGA,
        &layers,
    )
    .unwrap();
    let interface =
        PlanarDielectricInterface::try_new(&outside, DielectricResponse::Debye(&passive)).unwrap();
    let shallow = 0.3 / k();
    for pol in [Polarization::Te, Polarization::Tm] {
        let bank = stack.prepare_illumination(0.8, pol).unwrap();
        let response = interface.response(OMEGA, 0.8, pol).unwrap();
        let q = response.transmitted_normal_wavenumber_rad_m;
        let propagation = phase(-q.re * shallow).scale(det::exp(q.im * shallow));
        let expected = response.transmitted_tangential_e * propagation;
        let left = bank
            .prepare_point(0, shallow)
            .unwrap()
            .evaluate([C64::ONE, C64::ZERO])
            .unwrap();
        let right_point = bank.prepare_point_from_right(0, shallow).unwrap();
        assert_eq!(right_point.distance_from_right_m(), shallow);
        let right = right_point.evaluate([C64::ZERO, C64::ONE]).unwrap();
        let l = tangent(left.electric_v_m, left.magnetic_a_m, pol);
        let r = tangent(right.electric_v_m, right.magnetic_a_m, pol);
        complex_near(l[0], expected, 2e-13);
        complex_near(r[0], expected, 2e-13);
        complex_near(l[1], -r[1], 2e-15);
        near(
            left.statistics.absorbed_power_w_m3,
            right.statistics.absorbed_power_w_m3,
            1e-16,
        );
        let middle = bank
            .prepare_point(0, 0.5 * layers[0].thickness_m)
            .unwrap()
            .evaluate([C64::ONE, C64::ONE])
            .unwrap();
        assert_eq!(middle.electric_v_m, [C64::ZERO; 3]);
        assert_eq!(middle.magnetic_a_m, [C64::ZERO; 3]);
        assert_eq!(middle.statistics.absorbed_power_w_m3, 0.0);
    }
}

#[test]
fn partial_coherence_matches_dense_field_covariance_and_explicit_independence() {
    let outside = material(E, M);
    let passive = debye(3.0 * E, 4.0 * E, 0.1, 0.001);
    let layers = [DielectricLayer {
        response: DielectricResponse::Debye(&passive),
        thickness_m: 0.4 / k(),
    }];
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        OMEGA,
        &layers,
    )
    .unwrap();
    let factors = [
        [C64::new(0.8, 0.2), C64::new(-0.1, 0.4)],
        [C64::new(0.2, -0.3), C64::new(0.7, 0.1)],
    ];
    let covariance: [[C64; 2]; 2] = std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            factors
                .iter()
                .fold(C64::ZERO, |v, b| v + b[i] * b[j].conj())
        })
    });
    let budget = IlluminationBudget {
        maximum_modes: 2,
        maximum_layer_mode_pairs: 2,
    };
    for pol in [Polarization::Te, Polarization::Tm] {
        let bank = stack.prepare_illumination(0.65, pol).unwrap();
        let point = bank.prepare_point(0, 0.37 * layers[0].thickness_m).unwrap();
        let basis = [
            point.evaluate([C64::ONE, C64::ZERO]).unwrap(),
            point.evaluate([C64::ZERO, C64::ONE]).unwrap(),
        ];
        let contract = |a: [C64; 2], b: [C64; 2]| {
            let mut sum = C64::ZERO;
            for i in 0..2 {
                for j in 0..2 {
                    sum = sum + a[i] * covariance[i][j] * b[j].conj();
                }
            }
            sum.scale(0.5)
        };
        let actual = point
            .evaluate_coherence_modes(&factors, budget, || false)
            .unwrap();
        for i in 0..3 {
            let e = basis.map(|f| f.electric_v_m[i]);
            let h = basis.map(|f| f.magnetic_a_m[i]);
            near(
                actual.electric_rms_squared_v2_m2[i],
                contract(e, e).re,
                2e-14,
            );
            near(
                actual.magnetic_rms_squared_a2_m2[i],
                contract(h, h).re,
                2e-18,
            );
            let j = (i + 1) % 3;
            let l = (i + 2) % 3;
            let p = contract(
                basis.map(|f| f.electric_v_m[j]),
                basis.map(|f| f.magnetic_a_m[l]),
            ) - contract(
                basis.map(|f| f.electric_v_m[l]),
                basis.map(|f| f.magnetic_a_m[j]),
            );
            complex_near(actual.complex_poynting_w_m2[i], p, 2e-16);
        }
        near(
            actual.absorbed_power_w_m3,
            OMEGA * 2.0 * E * actual.electric_rms_squared_v2_m2.iter().sum::<f64>(),
            2e-16,
        );
        let rotated = factors.map(|b| b.map(|z| z * phase(0.79)));
        let rotation = point
            .evaluate_coherence_modes(&rotated, budget, || false)
            .unwrap();
        near(
            rotation.absorbed_power_w_m3,
            actual.absorbed_power_w_m3,
            2e-16,
        );
        let independent = point.evaluate_independent_sources([C64::ONE; 2]).unwrap();
        near(
            independent.absorbed_power_w_m3,
            basis[0].statistics.absorbed_power_w_m3 + basis[1].statistics.absorbed_power_w_m3,
            2e-16,
        );
        assert_eq!(
            point
                .evaluate_coherence_modes(&[], budget, || false)
                .unwrap(),
            Default::default()
        );
        assert!(
            point
                .evaluate_coherence_modes(
                    &factors,
                    IlluminationBudget {
                        maximum_modes: 1,
                        maximum_layer_mode_pairs: 2
                    },
                    || false
                )
                .is_err()
        );
        assert!(
            point
                .evaluate_coherence_modes(&factors, budget, || true)
                .is_err()
        );
        assert!(
            point
                .evaluate([C64::new(f64::NAN, 0.0), C64::ZERO])
                .is_err()
        );
    }
}

#[test]
fn point_queries_refuse_invalid_geometry_and_zero_volume_without_changing_stack_identity() {
    let outside = material(E, M);
    let layers = [
        DielectricLayer {
            response: DielectricResponse::Lossless(&outside),
            thickness_m: 0.0,
        },
        DielectricLayer {
            response: DielectricResponse::Lossless(&outside),
            thickness_m: 0.1 / k(),
        },
    ];
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        OMEGA,
        &layers,
    )
    .unwrap();
    let bank = stack.prepare_illumination(0.8, Polarization::Te).unwrap();
    assert!(bank.prepare_point(0, 0.0).is_err());
    assert!(bank.prepare_point(2, 0.0).is_err());
    for z in [-1.0, f64::NAN, f64::INFINITY, 2.0 * layers[1].thickness_m] {
        assert!(bank.prepare_point(1, z).is_err());
        assert!(bank.prepare_point_from_right(1, z).is_err());
    }
    let from_left = bank
        .prepare_point(1, 0.2 * layers[1].thickness_m)
        .unwrap()
        .evaluate([C64::ONE, C64::ONE])
        .unwrap();
    let from_right = bank
        .prepare_point_from_right(1, 0.8 * layers[1].thickness_m)
        .unwrap()
        .evaluate([C64::ONE, C64::ONE])
        .unwrap();
    for i in 0..3 {
        complex_near(from_left.electric_v_m[i], from_right.electric_v_m[i], 1e-14);
    }
}
