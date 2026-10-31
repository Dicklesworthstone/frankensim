use super::finite_layers::debye;
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
const BUDGET: IlluminationBudget = IlluminationBudget {
    maximum_modes: 16,
    maximum_layer_mode_pairs: 1024,
};

#[test]
fn simultaneous_bare_interface_fields_are_continuous_and_power_reciprocal() {
    let left = material(E, M);
    let right = material(4.0 * E, 1.3 * M);
    let stack =
        PlanarDielectricLayers::try_new(&left, DielectricResponse::Lossless(&right), OMEGA, &[])
            .unwrap();
    for pol in [Polarization::Te, Polarization::Tm] {
        let bank = stack.prepare_illumination(0.65, pol).unwrap();
        assert!(std::ptr::eq(bank.stack(), &stack));
        let s = bank.tangential_e_scattering();
        let y = bank.port_admittance_s();
        complex_near(s[0][1].scale(y[0]), s[1][0].scale(y[1]), 1e-16);
        let input = [C64::new(0.8, 0.3), C64::new(-0.1, 0.7)];
        let result = bank.evaluate(input, &mut []).unwrap();
        complex_near(
            result.boundary_tangential_e_v_m[0],
            result.boundary_tangential_e_v_m[1],
            1e-14,
        );
        complex_near(
            result.boundary_conjugate_h_a_m[0],
            result.boundary_conjugate_h_a_m[1],
            1e-16,
        );
        assert_eq!(result.power.absorbed_w_m2, 0.0);
        let single = stack.response(0.65, pol).unwrap();
        complex_near(s[0][0], single.reflected_tangential_e, 1e-14);
        complex_near(s[1][0], single.transmitted_tangential_e, 1e-14);
        for i in 0..2 {
            let flux = 0.5
                * (result.boundary_tangential_e_v_m[i] * result.boundary_conjugate_h_a_m[i].conj())
                    .re;
            let expected = if i == 0 {
                result.power.incoming_w_m2[i] - result.power.outgoing_w_m2[i]
            } else {
                result.power.outgoing_w_m2[i] - result.power.incoming_w_m2[i]
            };
            near(flux, expected, 1e-16);
        }
    }
}

// Independent conventional characteristic matrices and face continuity.
// These modest depths permit growing sin/cos; production must never need them.
fn matrix(
    epsilon_ratio: C64,
    mu_ratio: f64,
    depth: f64,
    c: f64,
    pol: Polarization,
) -> [[C64; 2]; 2] {
    let mut q = (epsilon_ratio.scale(mu_ratio) - C64::from_re((1.0 - c) * (1.0 + c))).sqrt();
    if q.im > 0.0 {
        q = -q;
    }
    let y = match pol {
        Polarization::Te => q.scale(1.0 / (mu_ratio * c)),
        Polarization::Tm => epsilon_ratio.scale(c) / q,
    };
    let z = q.scale(depth);
    let (s, co) = det::sin_cos(z.re);
    let ch = 0.5 * (det::exp(z.im) + det::exp(-z.im));
    let sh = 0.5 * (det::exp(z.im) - det::exp(-z.im));
    let cosine = C64::new(co * ch, -s * sh);
    let j_sine = C64::new(-co * sh, s * ch);
    [[cosine, j_sine / y], [j_sine * y, cosine]]
}
fn product(a: [[C64; 2]; 2], b: [[C64; 2]; 2]) -> [[C64; 2]; 2] {
    std::array::from_fn(|i| std::array::from_fn(|j| a[i][0] * b[0][j] + a[i][1] * b[1][j]))
}
fn act(a: [[C64; 2]; 2], x: [C64; 2]) -> [C64; 2] {
    a.map(|row| row[0] * x[0] + row[1] * x[1])
}

/// G1: independent boundary equations give all outgoing amplitudes; local
/// Poynting differences independently check each finite layer's actual heat.
#[test]
fn coherent_lossy_stack_matches_independent_boundary_and_layer_work() {
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
            for relative_phase in [0.0, 0.6, 2.1, 3.14] {
                let input = [C64::new(0.4, 0.7), phase(relative_phase).scale(0.8)];
                let transmitted =
                    (input[0].scale(2.0) - (minus[0] + minus[1]) * input[1]) / (plus[0] + plus[1]);
                let reflected = plus[0] * transmitted + minus[0] * input[1] - input[0];
                let mut heat = [0.0; 2];
                let result = bank.evaluate(input, &mut heat).unwrap();
                complex_near(result.outgoing_tangential_e_v_m[0], reflected, 2e-13);
                complex_near(result.outgoing_tangential_e_v_m[1], transmitted, 2e-13);
                let mut field = [
                    transmitted + input[1],
                    (transmitted - input[1]).scale(ratio),
                ];
                for (index, m) in [(1, m1), (0, m0)] {
                    let before = act(m, field);
                    let reference = 0.5
                        * y[0]
                        * ((before[0] * before[1].conj()).re - (field[0] * field[1].conj()).re);
                    near(heat[index], reference, 2e-15);
                    field = before;
                }
                assert_eq!(heat[1], 0.0, "transparent film cannot produce heat");
                near(result.power.absorbed_w_m2, heat[0], 1e-18);
            }
            let s = bank.tangential_e_scattering();
            complex_near(s[0][1].scale(y[0]), s[1][0].scale(y[1]), 1e-15);
        }
    }
}

/// G1/G3: a thin lossy layer can sit at an electric antinode or node depending
/// on source phase. Its weak odd loss remains positive instead of being lost
/// in subtraction of incoming/outgoing powers rounded to the same value.
#[test]
fn coherent_phase_controls_thin_film_heat_and_independent_sources_are_the_phase_average() {
    let outside = material(E, M);
    let film = debye(E, 2.0 * E, 0.03, 0.001);
    let layers = [DielectricLayer {
        response: DielectricResponse::Debye(&film),
        thickness_m: 1e-7 / k(),
    }];
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        OMEGA,
        &layers,
    )
    .unwrap();
    for pol in [Polarization::Te, Polarization::Tm] {
        let bank = stack.prepare_illumination(1.0, pol).unwrap();
        let mut heat = [0.0];
        let even = bank
            .evaluate([C64::ONE; 2], &mut heat)
            .unwrap()
            .power
            .absorbed_w_m2;
        let odd = bank
            .evaluate([C64::ONE, -C64::ONE], &mut heat)
            .unwrap()
            .power
            .absorbed_w_m2;
        assert!(odd > 0.0 && even > 1e10 * odd, "even={even}, odd={odd}");
        let independent = bank
            .evaluate_independent_sources([C64::ONE; 2], &mut heat)
            .unwrap();
        near(independent.absorbed_w_m2, 0.5 * (even + odd), 1e-13 * even);
        let one = bank
            .evaluate([C64::ONE, C64::ZERO], &mut heat)
            .unwrap()
            .power
            .absorbed_w_m2;
        near(independent.absorbed_w_m2, 2.0 * one, 1e-13 * even);
        let common = phase(1.123);
        near(
            bank.evaluate([common; 2], &mut heat)
                .unwrap()
                .power
                .absorbed_w_m2,
            even,
            1e-13 * even,
        );
    }
}

/// G0: complex source factors agree with dense covariance/Poynting algebra.
#[test]
fn partial_coherence_matches_full_source_covariance() {
    let outside = material(E, M);
    let film = debye(2.0 * E, 3.0 * E, 1.0, 0.01);
    let layers = [DielectricLayer {
        response: DielectricResponse::Debye(&film),
        thickness_m: 0.6 / k(),
    }];
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        OMEGA,
        &layers,
    )
    .unwrap();
    let bank = stack.prepare_illumination(0.73, Polarization::Tm).unwrap();
    let factors = [
        [C64::new(0.2, 0.5), C64::new(0.7, -0.1)],
        [C64::new(-0.3, 0.4), C64::new(0.2, 0.6)],
        [C64::from_re(0.8), C64::new(0.1, -0.2)],
    ];
    let mut c = [[C64::ZERO; 2]; 2];
    for i in 0..2 {
        for j in 0..2 {
            for b in factors {
                c[i][j] = c[i][j] + b[i] * b[j].conj();
            }
        }
    }
    let mut heat = [0.0];
    let result = bank
        .evaluate_coherence_modes(&factors, BUDGET, &mut heat, || false)
        .unwrap();
    let s = bank.tangential_e_scattering();
    let y = bank.port_admittance_s();
    for face in 0..2 {
        let mut out = C64::ZERO;
        for i in 0..2 {
            for j in 0..2 {
                out = out + s[face][i] * c[i][j] * s[face][j].conj();
            }
        }
        near(
            result.incoming_w_m2[face],
            0.5 * y[face] * c[face][face].re,
            1e-16,
        );
        near(result.outgoing_w_m2[face], 0.5 * y[face] * out.re, 1e-16);
    }
    near(
        heat[0],
        result.incoming_w_m2.iter().sum::<f64>() - result.outgoing_w_m2.iter().sum::<f64>(),
        1e-16,
    );
}

/// G1: critical/evanescent internal films and opaque absorption retain bounded
/// original slab limits for arbitrary simultaneous source phase.
#[test]
fn opaque_critical_and_evanescent_films_keep_physical_two_port_limits() {
    let outside = material(E, M);
    for (epsilon_ratio, cosine, depth) in [(0.75, 0.5, 0.4), (0.16, 0.6, 0.8), (0.16, 0.6, 1e6)] {
        let film = material(epsilon_ratio * E, M);
        let layers = [DielectricLayer {
            response: DielectricResponse::Lossless(&film),
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
            let bank = stack.prepare_illumination(cosine, pol).unwrap();
            let mut heat = [1.0];
            let result = bank
                .evaluate([C64::ONE, C64::new(0.2, 0.7)], &mut heat)
                .unwrap();
            assert_eq!(heat, [0.0]);
            near(
                result.power.incoming_w_m2.iter().sum(),
                result.power.outgoing_w_m2.iter().sum(),
                1e-15,
            );
        }
    }
    let film = debye(3.0 * E, 4.0 * E, 0.1, 0.003);
    let layers = [DielectricLayer {
        response: DielectricResponse::Debye(&film),
        thickness_m: 1e6 / k(),
    }];
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        OMEGA,
        &layers,
    )
    .unwrap();
    let bank = stack.prepare_illumination(0.8, Polarization::Te).unwrap();
    let s = bank.tangential_e_scattering();
    assert_eq!(s[0][1], C64::ZERO);
    assert_eq!(s[1][0], C64::ZERO);
    let mut heat = [0.0];
    let independent = bank
        .evaluate_independent_sources([C64::ONE; 2], &mut heat)
        .unwrap();
    let coherent = bank.evaluate([C64::ONE, phase(0.4)], &mut heat).unwrap();
    near(
        independent.absorbed_w_m2,
        coherent.power.absorbed_w_m2,
        1e-16,
    );
}

#[test]
fn unsupported_external_modes_shapes_nonfinite_and_work_refuse() {
    let outside = material(E, M);
    let loss = debye(E, E, 1.0, 0.01);
    let opaque_exit =
        PlanarDielectricLayers::try_new(&outside, DielectricResponse::Debye(&loss), OMEGA, &[])
            .unwrap();
    assert!(
        opaque_exit
            .prepare_illumination(1.0, Polarization::Te)
            .is_err()
    );
    let low = material(0.16 * E, M);
    let evanescent_exit =
        PlanarDielectricLayers::try_new(&outside, DielectricResponse::Lossless(&low), OMEGA, &[])
            .unwrap();
    assert!(
        evanescent_exit
            .prepare_illumination(0.5, Polarization::Tm)
            .is_err()
    );
    let layers = [DielectricLayer {
        response: DielectricResponse::Debye(&loss),
        thickness_m: 0.1 / k(),
    }];
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        OMEGA,
        &layers,
    )
    .unwrap();
    assert!(stack.prepare_illumination(0.0, Polarization::Te).is_err());
    let bank = stack.prepare_illumination(1.0, Polarization::Te).unwrap();
    assert!(bank.evaluate([C64::ONE; 2], &mut []).is_err());
    assert!(
        bank.evaluate([C64::new(f64::NAN, 0.0), C64::ZERO], &mut [0.0])
            .is_err()
    );
    assert!(
        bank.evaluate_coherence_modes(
            &[[C64::ONE; 2]],
            IlluminationBudget {
                maximum_modes: 1,
                maximum_layer_mode_pairs: 0
            },
            &mut [0.0],
            || false
        )
        .is_err()
    );
    assert!(
        bank.evaluate_coherence_modes(&[[C64::ONE; 2]], BUDGET, &mut [0.0], || true)
            .is_err()
    );
    let mut heat = [1.0];
    let silent = bank
        .evaluate_coherence_modes(&[], BUDGET, &mut heat, || false)
        .unwrap();
    assert_eq!(silent.absorbed_w_m2, 0.0);
    assert_eq!(heat, [0.0]);
}
