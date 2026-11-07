//! Independent two-source transfer equations and volume work, G1/G3 only.
use super::super::tests::{
    admittance, boundary, complex_near, electric, field, field_near, films, input, matrix, near,
    power, product, propagated,
};
use super::*;

fn right_electric() -> [C64; 2] {
    [C64::new(0.2, -0.4), C64::new(0.6, 0.1)]
}

// Solve the dimensional characteristic matrix for both unknown outgoing
// waves at once. This oracle never uses either production scattering column.
fn outgoing(data: Input, films: &[FilmInput], left: [C64; 2], right: [C64; 2]) -> [[C64; 2]; 2] {
    let mut result = [[C64::ZERO; 2]; 2];
    for c in 0..2 {
        let mut m = [[C64::ONE, C64::ZERO], [C64::ZERO, C64::ONE]];
        for film in films {
            m = product(m, matrix(data, *film, c, film.thickness));
        }
        let yl = admittance(data, data.epsilon_i, data.mu_i, c);
        let yr = admittance(data, data.epsilon_t, data.mu_t, c);
        let f = m[0][0] + m[0][1] * yr;
        let g = m[0][0] - m[0][1] * yr;
        let h = m[1][0] + m[1][1] * yr;
        let k = m[1][0] - m[1][1] * yr;
        result[1][c] = (yl.scale(2.0) * left[c] - (yl * g + k) * right[c]) / (yl * f + h);
        result[0][c] = f * result[1][c] + g * right[c] - left[c];
    }
    result
}
fn reverse_input(data: Input) -> Input {
    Input {
        epsilon_i: data.epsilon_t,
        mu_i: data.mu_t,
        epsilon_t: data.epsilon_i,
        mu_t: data.mu_i,
        ..data
    }
}
fn reorient(value: InterfaceBoundaryField) -> InterfaceBoundaryField {
    InterfaceBoundaryField {
        electric_v_m: [
            value.electric_v_m[0],
            value.electric_v_m[1],
            -value.electric_v_m[2],
        ],
        magnetic_a_m: [
            -value.magnetic_a_m[0],
            -value.magnetic_a_m[1],
            value.magnetic_a_m[2],
        ],
    }
}
fn port_checks(result: AbsorbingTwoSidedReflection) {
    let mut sum = 0.0;
    for face in 0..2 {
        let sign = if face == 0 { 1.0 } else { -1.0 };
        let p = result.power[face];
        near(sign * power(result.incoming[face]), p.incoming_w_m2, 2e-13);
        near(-sign * power(result.outgoing[face]), p.outgoing_w_m2, 2e-13);
        near(
            sign * power(result.boundary[face]),
            p.net_into_stack_w_m2,
            2e-13,
        );
        near(
            p.incoming_w_m2 - p.outgoing_w_m2 + p.interference_into_stack_w_m2,
            p.net_into_stack_w_m2,
            2e-13,
        );
        sum += p.net_into_stack_w_m2;
    }
    near(sum, result.absorbed_power_w_m2, 4e-13);
}

#[test]
fn matched_absorbing_film_preserves_phase_control_and_both_exterior_cross_fluxes() {
    let data = Input {
        omega: 1.0,
        kx: 0.0,
        epsilon_i: C64::new(-3.0, -4.0),
        epsilon_t: C64::new(-3.0, -4.0),
        mu_i: 1.0,
        mu_t: 1.0,
    };
    let d = 0.4_f64;
    let source = Prepared::new(
        data,
        &[FilmInput {
            epsilon: data.epsilon_i,
            mu: 1.0,
            thickness: d,
        }],
    )
    .unwrap();
    let prepared = PreparedTwoSided::new(&source).unwrap();
    let p = C64::new((-2.0 * d).exp() * d.cos(), -(-2.0 * d).exp() * d.sin());
    for channel in 0..2 {
        let mut left = [C64::ZERO; 2];
        left[channel] = C64::ONE;
        let mut heat = Vec::new();
        for b in [C64::ONE, -C64::ONE, C64::new(0.0, 1.0), C64::new(0.3, -0.6)] {
            let mut right = [C64::ZERO; 2];
            right[channel] = b;
            let solved = prepared.illuminate(left, right).unwrap();
            let result = solved.response;
            complex_near(result.outgoing[0].electric_v_m[channel], p * b, 3e-14);
            complex_near(result.outgoing[1].electric_v_m[channel], p, 3e-14);
            let expected = 0.5 * (1.0 - p.norm_sq()) * (1.0 + b.norm_sq()) - 4.0 * p.im * b.re;
            near(result.absorbed_power_w_m2, expected, 8e-14);
            near(solved.absorbed[0], expected, 8e-14);
            for z in [0.0, 0.07, 0.2, 0.4] {
                let decay =
                    |x: f64| C64::new((-2.0 * x).exp() * x.cos(), -(-2.0 * x).exp() * x.sin());
                let e = decay(z) + b * decay(d - z);
                complex_near(
                    prepared
                        .field(left, right, 0, z)
                        .unwrap()
                        .field
                        .electric_v_m[channel],
                    e,
                    5e-14,
                );
            }
            port_checks(result);
            heat.push(result.absorbed_power_w_m2);
        }
        assert!(
            heat[0] > 10.0 * heat[1],
            "source phase must control physical absorption"
        );
        assert!(
            heat[0] > 1.0,
            "positive entrance interference can exceed separate source powers"
        );
    }
}

#[test]
fn oblique_two_source_stack_matches_independent_fields_boundaries_and_volume_loss() {
    for kx in [-1.1, 0.45] {
        let data = Input { kx, ..input() };
        let layers = films();
        let source = Prepared::new(data, &layers).unwrap();
        let prepared = PreparedTwoSided::new(&source).unwrap();
        let left = electric();
        let right = right_electric();
        let solved = prepared.illuminate(left, right).unwrap();
        let result = solved.response;
        let expected = outgoing(data, &layers, left, right);
        for face in 0..2 {
            for c in 0..2 {
                complex_near(
                    result.outgoing[face].electric_v_m[c],
                    expected[face][c],
                    2e-12,
                );
            }
        }
        let mut state = core::array::from_fn(|c| {
            [
                left[c] + expected[0][c],
                admittance(data, data.epsilon_i, data.mu_i, c) * (left[c] - expected[0][c]),
            ]
        });
        let mut preceding = result.boundary[0];
        let mut epsilon = data.epsilon_i;
        let mut mu = data.mu_i;
        for (i, film) in layers.iter().copied().enumerate() {
            let point = prepared.field(left, right, i, 0.0).unwrap();
            boundary(preceding, epsilon, mu, point.field, film.epsilon, film.mu);
            for fraction in [0.0, 0.23, 0.71, 1.0] {
                let z = fraction * film.thickness;
                field_near(
                    prepared.field(left, right, i, z).unwrap().field,
                    field(data, film, propagated(data, film, state, z)),
                    3e-12,
                );
            }
            let heat = |z: f64| {
                -0.5 * data.omega
                    * film.epsilon.im
                    * field(data, film, propagated(data, film, state, z))
                        .electric_v_m
                        .iter()
                        .map(|e| e.norm_sq())
                        .sum::<f64>()
            };
            let n = 512;
            let dz = film.thickness / (n as f64);
            let mut integral = heat(0.0) + heat(film.thickness);
            for j in 1..n {
                integral += (if j % 2 == 0 { 2.0 } else { 4.0 }) * heat((j as f64) * dz);
            }
            near(integral * dz / 3.0, solved.absorbed[i], 3e-11);
            state = propagated(data, film, state, film.thickness);
            preceding = prepared
                .field(left, right, i, film.thickness)
                .unwrap()
                .field;
            epsilon = film.epsilon;
            mu = film.mu;
        }
        boundary(
            preceding,
            epsilon,
            mu,
            result.boundary[1],
            data.epsilon_t,
            data.mu_t,
        );
        port_checks(result);
    }
}

#[test]
fn either_single_source_recovers_original_solution_and_lossless_ports_have_no_cross_flux() {
    for lossless in [false, true] {
        let mut data = input();
        if lossless {
            data.epsilon_i.im = 0.0;
            data.epsilon_t.im = 0.0;
        }
        let layers = films();
        let source = Prepared::new(data, &layers).unwrap();
        let reversed: Vec<_> = layers.iter().rev().copied().collect();
        let reverse = Prepared::new(reverse_input(data), &reversed).unwrap();
        let prepared = PreparedTwoSided::new(&source).unwrap();
        let a = electric();
        let b = right_electric();
        for right_on in [false, true] {
            let (left, right) = if right_on {
                ([C64::ZERO; 2], b)
            } else {
                (a, [C64::ZERO; 2])
            };
            let actual = prepared.illuminate(left, right).unwrap();
            let prior = if right_on {
                reverse.illuminate(b).unwrap()
            } else {
                source.illuminate(a).unwrap()
            };
            let from = usize::from(right_on);
            let to = 1 - from;
            near(
                actual.response.power[from].incoming_w_m2,
                prior.response.incident_normal_power_w_m2,
                2e-13,
            );
            near(
                actual.response.power[from].outgoing_w_m2,
                -prior.response.reflected_normal_power_w_m2,
                2e-13,
            );
            near(
                actual.response.power[from].interference_into_stack_w_m2,
                prior.response.interference_normal_power_w_m2,
                2e-13,
            );
            near(
                actual.response.power[to].outgoing_w_m2,
                prior.response.transmitted_normal_power_w_m2,
                2e-13,
            );
            assert_eq!(actual.response.power[to].interference_into_stack_w_m2, 0.0);
            for (i, film) in layers.iter().enumerate() {
                let index = if right_on { layers.len() - 1 - i } else { i };
                near(actual.absorbed[i], prior.absorbed[index], 2e-13);
                let z = 0.37 * film.thickness;
                let expected = if right_on {
                    reorient(reverse.field(b, index, film.thickness - z).unwrap().field)
                } else {
                    source.field(a, i, z).unwrap().field
                };
                field_near(
                    prepared.field(left, right, i, z).unwrap().field,
                    expected,
                    3e-12,
                );
            }
        }
        let result = prepared.illuminate(a, b).unwrap().response;
        if lossless {
            for p in result.power {
                assert_eq!(p.interference_into_stack_w_m2, 0.0);
            }
        }
        port_checks(result);
    }
}

#[test]
fn critical_and_opaque_films_retain_global_field_orientation_and_positive_heat() {
    let data = Input {
        omega: 1.0,
        kx: 1.0,
        ..input()
    };
    let film = FilmInput {
        epsilon: C64::ONE,
        mu: 1.0,
        thickness: 0.7,
    };
    let source = Prepared::new(data, &[film]).unwrap();
    let prepared = PreparedTwoSided::new(&source).unwrap();
    let left = electric();
    let right = right_electric();
    let result = prepared.illuminate(left, right).unwrap();
    assert_eq!(result.absorbed[0], 0.0);
    let lo = prepared.field(left, right, 0, 0.0).unwrap().field;
    let hi = prepared
        .field(left, right, 0, film.thickness)
        .unwrap()
        .field;
    let mid = prepared
        .field(left, right, 0, 0.5 * film.thickness)
        .unwrap();
    for c in 0..3 {
        complex_near(
            mid.field.electric_v_m[c],
            (lo.electric_v_m[c] + hi.electric_v_m[c]).scale(0.5),
            2e-13,
        );
        complex_near(
            mid.field.magnetic_a_m[c],
            (lo.magnetic_a_m[c] + hi.magnetic_a_m[c]).scale(0.5),
            2e-13,
        );
    }
    assert_eq!(mid.heating_w_m3, 0.0);
    port_checks(result.response);
    let opaque = FilmInput {
        epsilon: C64::new(2.0, -0.8),
        mu: 1.2,
        thickness: 1e6,
    };
    let source = Prepared::new(input(), &[opaque]).unwrap();
    let prepared = PreparedTwoSided::new(&source).unwrap();
    let both = prepared.illuminate(left, right).unwrap();
    let one = prepared.illuminate(left, [C64::ZERO; 2]).unwrap();
    let other = prepared.illuminate([C64::ZERO; 2], right).unwrap();
    near(both.absorbed[0], one.absorbed[0] + other.absorbed[0], 3e-13);
    for c in 0..2 {
        assert_eq!(prepared.scattering()[c][0][1], C64::ZERO);
        assert_eq!(prepared.scattering()[c][1][0], C64::ZERO);
    }
    assert_eq!(
        prepared.field(left, right, 0, 5e5).unwrap().heating_w_m3,
        0.0
    );
    port_checks(both.response);
}

#[test]
fn odd_coherent_illumination_retains_weak_heat_that_exterior_subtraction_erases() {
    let loss = 1e-20;
    let d = 0.001;
    let data = Input {
        omega: 1.0,
        kx: 0.0,
        epsilon_i: C64::ONE,
        epsilon_t: C64::ONE,
        mu_i: 1.0,
        mu_t: 1.0,
    };
    let source = Prepared::new(
        data,
        &[FilmInput {
            epsilon: C64::new(1.0, -loss),
            mu: 1.0,
            thickness: d,
        }],
    )
    .unwrap();
    let prepared = PreparedTwoSided::new(&source).unwrap();
    for c in 0..2 {
        let mut a = [C64::ZERO; 2];
        a[c] = C64::ONE;
        let b = a.map(|v| -v);
        let result = prepared.illuminate(a, b).unwrap();
        let x = d * d;
        let expected = loss * d * x * (1.0 / 6.0 - x * (1.0 / 120.0 - x / 5040.0));
        assert!(result.absorbed[0] > 0.0);
        assert!((result.absorbed[0] / expected - 1.0).abs() < 2e-10);
        near(
            prepared.field(a, b, 0, 0.5 * d).unwrap().heating_w_m3,
            0.0,
            1e-30,
        );
    }
}

#[test]
fn dark_fields_empty_stack_and_exterior_port_admission_remain_explicit() {
    let source = Prepared::new(input(), &[]).unwrap();
    let prepared = PreparedTwoSided::new(&source).unwrap();
    let result = prepared.illuminate(electric(), right_electric()).unwrap();
    boundary(
        result.response.boundary[0],
        input().epsilon_i,
        input().mu_i,
        result.response.boundary[1],
        input().epsilon_t,
        input().mu_t,
    );
    assert_eq!(result.response.absorbed_power_w_m2, 0.0);
    port_checks(result.response);
    let dark = prepared
        .illuminate([C64::ZERO; 2], [C64::ZERO; 2])
        .unwrap()
        .response;
    assert_eq!(dark.absorbed_power_w_m2, 0.0);
    for field in dark
        .incoming
        .into_iter()
        .chain(dark.outgoing)
        .chain(dark.boundary)
    {
        assert_eq!(field.electric_v_m, [C64::ZERO; 3]);
        assert_eq!(field.magnetic_a_m, [C64::ZERO; 3]);
    }
    assert!(
        prepared
            .illuminate([C64::new(f64::NAN, 0.0); 2], electric())
            .is_err()
    );
    assert!(
        prepared
            .field(electric(), right_electric(), 0, 0.0)
            .is_err()
    );
    for kx in [1.0, 1.2] {
        let data = Input {
            omega: 1.0,
            kx,
            epsilon_t: C64::ONE,
            mu_t: 1.0,
            ..input()
        };
        let source = Prepared::new(data, &[]).unwrap();
        source.illuminate(electric()).unwrap();
        assert!(PreparedTwoSided::new(&source).is_err());
    }
}
