//! G1/G3: dimensional transfer matrices, full boundary fields and volume work.
use super::*;

pub(super) fn input() -> Input {
    Input {
        omega: 1.7,
        kx: 0.45,
        epsilon_i: C64::new(2.5, -0.8),
        mu_i: 1.4,
        epsilon_t: C64::new(0.7, -0.2),
        mu_t: 0.9,
    }
}
pub(super) fn films() -> [FilmInput; 3] {
    [
        FilmInput {
            epsilon: C64::new(0.9, -0.3),
            mu: 0.8,
            thickness: 0.37,
        },
        FilmInput {
            epsilon: C64::new(-1.2, -0.5),
            mu: 1.3,
            thickness: 0.23,
        },
        FilmInput {
            epsilon: C64::new(3.1, -0.1),
            mu: 0.9,
            thickness: 0.41,
        },
    ]
}
pub(super) fn electric() -> [C64; 2] {
    [C64::new(0.7, 0.2), C64::new(-0.3, 0.5)]
}
pub(super) fn near(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + a.abs() + b.abs()),
        "{a:e} != {b:e}"
    );
}
pub(super) fn complex_near(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * (1.0 + a.abs() + b.abs()),
        "{a:?} != {b:?}"
    );
}
pub(super) fn root(data: Input, epsilon: C64, mu: f64) -> C64 {
    let q = (epsilon.scale(data.omega * data.omega * mu) - C64::from_re(data.kx * data.kx)).sqrt();
    if q.im > 0.0 { -q } else { q }
}
pub(super) fn admittance(data: Input, epsilon: C64, mu: f64, channel: usize) -> C64 {
    let q = root(data, epsilon, mu);
    if channel == 0 {
        epsilon.scale(data.omega) / q
    } else {
        q.scale(1.0 / (data.omega * mu))
    }
}

// Independent moderate-depth characteristic matrix. It propagates right-face
// [E,H_conjugate] to the left face, with exp(+i omega t), exp(-i q z).
// Direct cosh/sinh are deliberate: this test oracle is not the production
// scaled scattering/positive-loss recurrence and is not used for opaque films.
pub(super) fn matrix(data: Input, film: FilmInput, channel: usize, distance: f64) -> [[C64; 2]; 2] {
    let q = root(data, film.epsilon, film.mu);
    let q2 =
        film.epsilon.scale(data.omega * data.omega * film.mu) - C64::from_re(data.kx * data.kx);
    let (a, b) = if channel == 0 {
        (
            q2 / film.epsilon.scale(data.omega),
            film.epsilon.scale(data.omega),
        )
    } else {
        (
            C64::from_re(data.omega * film.mu),
            q2.scale(1.0 / (data.omega * film.mu)),
        )
    };
    let x = q.re * distance;
    let y = q.im * distance;
    let cosine = C64::new(x.cos() * y.cosh(), -x.sin() * y.sinh());
    let sine_over_q = if q == C64::ZERO {
        C64::from_re(distance)
    } else {
        C64::new(x.sin() * y.cosh(), x.cos() * y.sinh()) / q
    };
    let coupling = C64::new(0.0, 1.0) * sine_over_q;
    [[cosine, coupling * a], [coupling * b, cosine]]
}
pub(super) fn product(a: [[C64; 2]; 2], b: [[C64; 2]; 2]) -> [[C64; 2]; 2] {
    core::array::from_fn(|i| core::array::from_fn(|j| a[i][0] * b[0][j] + a[i][1] * b[1][j]))
}
fn gains(data: Input, films: &[FilmInput]) -> ([C64; 2], [C64; 2]) {
    let mut r = [C64::ZERO; 2];
    let mut t = [C64::ZERO; 2];
    for c in 0..2 {
        let mut m = [[C64::ONE, C64::ZERO], [C64::ZERO, C64::ONE]];
        for film in films {
            m = product(m, matrix(data, *film, c, film.thickness));
        }
        let yi = admittance(data, data.epsilon_i, data.mu_i, c);
        let yt = admittance(data, data.epsilon_t, data.mu_t, c);
        let f = m[0][0] + m[0][1] * yt;
        let g = m[1][0] + m[1][1] * yt;
        r[c] = (yi * f - g) / (yi * f + g);
        t[c] = yi.scale(2.0) / (yi * f + g);
    }
    (r, t)
}
fn initial_state(data: Input, r: [C64; 2]) -> [[C64; 2]; 2] {
    core::array::from_fn(|c| {
        [
            electric()[c] * (C64::ONE + r[c]),
            admittance(data, data.epsilon_i, data.mu_i, c) * electric()[c] * (C64::ONE - r[c]),
        ]
    })
}
pub(super) fn propagated(
    data: Input,
    film: FilmInput,
    left: [[C64; 2]; 2],
    z: f64,
) -> [[C64; 2]; 2] {
    core::array::from_fn(|c| {
        let m = matrix(data, film, c, -z);
        [
            m[0][0] * left[c][0] + m[0][1] * left[c][1],
            m[1][0] * left[c][0] + m[1][1] * left[c][1],
        ]
    })
}
pub(super) fn field(data: Input, film: FilmInput, state: [[C64; 2]; 2]) -> InterfaceBoundaryField {
    InterfaceBoundaryField {
        electric_v_m: [
            state[0][0],
            state[1][0],
            -state[0][1].scale(data.kx) / film.epsilon.scale(data.omega),
        ],
        magnetic_a_m: [
            -state[1][1],
            state[0][1],
            state[1][0].scale(data.kx / (data.omega * film.mu)),
        ],
    }
}
pub(super) fn field_near(a: InterfaceBoundaryField, b: InterfaceBoundaryField, tolerance: f64) {
    for c in 0..3 {
        complex_near(a.electric_v_m[c], b.electric_v_m[c], tolerance);
        complex_near(a.magnetic_a_m[c], b.magnetic_a_m[c], tolerance);
    }
}
pub(super) fn power(field: InterfaceBoundaryField) -> f64 {
    0.5 * (field.electric_v_m[0] * field.magnetic_a_m[1].conj()
        - field.electric_v_m[1] * field.magnetic_a_m[0].conj())
    .re
}
pub(super) fn combined(
    a: InterfaceBoundaryField,
    b: InterfaceBoundaryField,
) -> InterfaceBoundaryField {
    InterfaceBoundaryField {
        electric_v_m: core::array::from_fn(|c| a.electric_v_m[c] + b.electric_v_m[c]),
        magnetic_a_m: core::array::from_fn(|c| a.magnetic_a_m[c] + b.magnetic_a_m[c]),
    }
}
pub(super) fn boundary(
    left: InterfaceBoundaryField,
    epsilon_l: C64,
    mu_l: f64,
    right: InterfaceBoundaryField,
    epsilon_r: C64,
    mu_r: f64,
) {
    for c in 0..2 {
        complex_near(left.electric_v_m[c], right.electric_v_m[c], 2e-12);
        complex_near(left.magnetic_a_m[c], right.magnetic_a_m[c], 2e-12);
    }
    complex_near(
        epsilon_l * left.electric_v_m[2],
        epsilon_r * right.electric_v_m[2],
        2e-12,
    );
    complex_near(
        left.magnetic_a_m[2].scale(mu_l),
        right.magnetic_a_m[2].scale(mu_r),
        2e-12,
    );
}

#[test]
fn g1_absorbing_stack_matches_independent_dimensional_transfer_matrices() {
    for kx in [-1.4, -0.45, 0.0, 0.45, 1.4] {
        let data = Input { kx, ..input() };
        let layers = films();
        let prepared = Prepared::new(data, &layers).unwrap();
        let result = prepared.illuminate(electric()).unwrap().response;
        let (r, t) = gains(data, &layers);
        for c in 0..2 {
            complex_near(result.reflected_tangential_e_gain[c], r[c], 2e-12);
            complex_near(result.transmitted_tangential_e_gain[c], t[c], 2e-12);
        }
        let mut state = initial_state(data, r);
        for (index, film) in layers.into_iter().enumerate() {
            for fraction in [0.0, 0.19, 0.53, 1.0] {
                let z = fraction * film.thickness;
                let expected = field(data, film, propagated(data, film, state, z));
                field_near(
                    prepared.field(electric(), index, z).unwrap().field,
                    expected,
                    3e-12,
                );
            }
            state = propagated(data, film, state, film.thickness);
        }
    }
}

#[test]
fn g1_full_boundary_fields_and_independent_volume_work_close_with_incident_interference() {
    let data = input();
    let layers = films();
    let prepared = Prepared::new(data, &layers).unwrap();
    let illuminated = prepared.illuminate(electric()).unwrap();
    let result = illuminated.response;
    let mut preceding = combined(result.incident, result.reflected);
    let mut epsilon = data.epsilon_i;
    let mut mu = data.mu_i;
    let mut state = initial_state(data, gains(data, &layers).0);
    let mut total = 0.0;
    for (index, film) in layers.into_iter().enumerate() {
        let left = prepared.field(electric(), index, 0.0).unwrap();
        let right = prepared.field(electric(), index, film.thickness).unwrap();
        boundary(preceding, epsilon, mu, left.field, film.epsilon, film.mu);
        let heat = |z: f64| {
            let exact = field(data, film, propagated(data, film, state, z));
            -0.5 * data.omega
                * film.epsilon.im
                * exact.electric_v_m.iter().map(|e| e.norm_sq()).sum::<f64>()
        };
        // Simpson integration of the independent E(z), not boundary subtraction.
        let intervals = 512;
        let dz = film.thickness / intervals as f64;
        let mut integral = heat(0.0) + heat(film.thickness);
        for i in 1..intervals {
            integral += (if i % 2 == 0 { 2.0 } else { 4.0 }) * heat(i as f64 * dz);
        }
        integral *= dz / 3.0;
        assert!(integral > 0.0);
        near(illuminated.absorbed[index], integral, 2e-11);
        near(power(left.field) - power(right.field), integral, 2e-11);
        for z in [0.0, 0.37 * film.thickness, film.thickness] {
            let point = prepared.field(electric(), index, z).unwrap();
            near(point.heating_w_m3, heat(z), 3e-12);
            near(point.normal_power_w_m2, power(point.field), 2e-13);
        }
        total += integral;
        preceding = right.field;
        epsilon = film.epsilon;
        mu = film.mu;
        state = propagated(data, film, state, film.thickness);
    }
    boundary(
        preceding,
        epsilon,
        mu,
        result.transmitted,
        data.epsilon_t,
        data.mu_t,
    );
    near(result.absorbed_power_w_m2, total, 3e-11);
    near(
        power(result.incident),
        result.incident_normal_power_w_m2,
        3e-13,
    );
    near(
        power(result.reflected),
        result.reflected_normal_power_w_m2,
        3e-13,
    );
    near(
        power(result.transmitted),
        result.transmitted_normal_power_w_m2,
        3e-13,
    );
    let total_in = power(combined(result.incident, result.reflected));
    near(
        total_in - power(result.incident) - power(result.reflected),
        result.interference_normal_power_w_m2,
        3e-13,
    );
    near(
        total_in,
        result.transmitted_normal_power_w_m2 + total,
        3e-11,
    );
    assert!(result.interference_normal_power_w_m2.abs() > 1e-4);
}

#[test]
fn g1_empty_stack_retains_the_exact_absorbing_port_interference() {
    let data = Input {
        omega: 1.0,
        kx: 0.0,
        epsilon_i: C64::new(-3.0, -4.0),
        mu_i: 1.0,
        epsilon_t: C64::ONE,
        mu_t: 1.0,
    };
    for e in [[C64::ONE, C64::ZERO], [C64::ZERO, C64::ONE]] {
        let result = Prepared::new(data, &[])
            .unwrap()
            .illuminate(e)
            .unwrap()
            .response;
        for c in 0..2 {
            complex_near(
                result.reflected_tangential_e_gain[c],
                C64::new(0.5, -0.5),
                2e-14,
            );
            complex_near(
                result.transmitted_tangential_e_gain[c],
                C64::new(1.5, -0.5),
                2e-14,
            );
        }
        near(result.incident_normal_power_w_m2, 0.5, 2e-14);
        near(result.reflected_normal_power_w_m2, -0.25, 2e-14);
        near(result.interference_normal_power_w_m2, 1.0, 2e-14);
        near(result.transmitted_normal_power_w_m2, 1.25, 2e-14);
        assert_eq!(result.absorbed_power_w_m2, 0.0);
    }
}

#[test]
fn g3_splitting_a_lossy_film_preserves_fields_and_original_film_heat() {
    let data = input();
    let film = FilmInput {
        epsilon: C64::new(2.1, -0.7),
        mu: 0.9,
        thickness: 0.7,
    };
    let split = [
        FilmInput {
            thickness: 0.2,
            ..film
        },
        FilmInput {
            thickness: 0.5,
            ..film
        },
    ];
    let one = Prepared::new(data, &[film]).unwrap();
    let two = Prepared::new(data, &split).unwrap();
    let a = one.illuminate(electric()).unwrap();
    let b = two.illuminate(electric()).unwrap();
    for c in 0..2 {
        complex_near(
            a.response.reflected_tangential_e_gain[c],
            b.response.reflected_tangential_e_gain[c],
            2e-12,
        );
        complex_near(
            a.response.transmitted_tangential_e_gain[c],
            b.response.transmitted_tangential_e_gain[c],
            2e-12,
        );
    }
    near(a.absorbed[0], b.absorbed[0] + b.absorbed[1], 2e-12);
    for z in [0.0, 0.1, 0.2, 0.35, 0.7] {
        let (index, depth) = if z <= 0.2 { (0, z) } else { (1, z - 0.2) };
        field_near(
            one.field(electric(), 0, z).unwrap().field,
            two.field(electric(), index, depth).unwrap().field,
            2e-12,
        );
    }
}

#[test]
fn g1_exactly_critical_lossless_film_keeps_linear_maxwell_fields() {
    let data = Input {
        omega: 1.0,
        kx: 1.0,
        epsilon_i: C64::new(4.0, -0.2),
        mu_i: 1.0,
        epsilon_t: C64::new(2.0, -0.1),
        mu_t: 1.0,
    };
    let film = FilmInput {
        epsilon: C64::ONE,
        mu: 1.0,
        thickness: 1.25,
    };
    let prepared = Prepared::new(data, &[film]).unwrap();
    let illuminated = prepared.illuminate(electric()).unwrap();
    let (r, t) = gains(data, &[film]);
    for c in 0..2 {
        complex_near(
            illuminated.response.reflected_tangential_e_gain[c],
            r[c],
            3e-12,
        );
        complex_near(
            illuminated.response.transmitted_tangential_e_gain[c],
            t[c],
            3e-12,
        );
    }
    for z in [0.0, 0.4, film.thickness] {
        let point = prepared.field(electric(), 0, z).unwrap();
        field_near(
            point.field,
            field(
                data,
                film,
                propagated(data, film, initial_state(data, r), z),
            ),
            3e-12,
        );
        assert_eq!(point.heating_w_m3, 0.0);
    }
    assert_eq!(illuminated.absorbed[0], 0.0);
}

#[test]
fn g1_opaque_and_weakly_absorbing_films_keep_positive_volume_heat() {
    let data = Input {
        omega: 1.0,
        kx: 0.0,
        epsilon_i: C64::new(1.0, -0.1),
        mu_i: 1.0,
        epsilon_t: C64::ONE,
        mu_t: 1.0,
    };
    let e = [C64::ONE, C64::ZERO];
    let opaque = FilmInput {
        epsilon: C64::new(2.0, -0.6),
        mu: 1.0,
        thickness: 1e4,
    };
    let prepared = Prepared::new(data, &[opaque]).unwrap();
    let result = prepared.illuminate(e).unwrap();
    let yi = admittance(data, data.epsilon_i, data.mu_i, 0);
    let yf = admittance(data, opaque.epsilon, opaque.mu, 0);
    complex_near(
        result.response.reflected_tangential_e_gain[0],
        (yi - yf) / (yi + yf),
        2e-12,
    );
    assert_eq!(result.response.transmitted_normal_power_w_m2, 0.0);
    near(
        result.absorbed[0],
        power(combined(
            result.response.incident,
            result.response.reflected,
        )),
        2e-12,
    );
    assert!(result.absorbed[0] > 0.0);
    assert_eq!(
        prepared.field(e, 0, opaque.thickness).unwrap().heating_w_m3,
        0.0
    );

    let weak = FilmInput {
        epsilon: C64::new(1.0, -1e-20),
        mu: 1.0,
        thickness: 0.4,
    };
    let result = Prepared::new(data, &[weak]).unwrap().illuminate(e).unwrap();
    let t = yi.scale(2.0) / (yi + C64::ONE);
    let expected = 0.5 * 1e-20 * weak.thickness * t.norm_sq();
    assert!(result.absorbed[0] > 0.0);
    assert!(
        (result.absorbed[0] / expected - 1.0).abs() < 2e-12,
        "weak-film heat {:e} != {:e}",
        result.absorbed[0],
        expected
    );
}
