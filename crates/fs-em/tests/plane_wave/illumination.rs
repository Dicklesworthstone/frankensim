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

/// G1/G3: a near-dark outgoing channel retains the exact represented source
/// perturbation. The oracle reduces a*(b+delta)-b*a to a*delta, with no
/// cancellation or high-precision implementation shared with the observer.
#[test]
fn coherent_interface_retains_weak_dark_port_leakage() {
    let left = material(E, M);
    let right = material(4.0 * E, M);
    let stack =
        PlanarDielectricLayers::try_new(&left, DielectricResponse::Lossless(&right), OMEGA, &[])
            .unwrap();
    for pol in [Polarization::Te, Polarization::Tm] {
        let bank = stack.prepare_illumination(1.0, pol).unwrap();
        let [a, b] = bank.tangential_e_scattering()[0];
        assert_eq!(a.im, 0.0);
        assert_eq!(b.im, 0.0);
        let perturbed = b.re.next_up();
        let delta = perturbed - b.re;
        for exponent in [-450, 0, 450] {
            let scale = 2.0_f64.powi(exponent);
            let expected = a.scale(delta * scale);
            let response = bank
                .evaluate([C64::from_re(perturbed * scale), -a.scale(scale)], &mut [])
                .unwrap();
            let leakage = response.outgoing_tangential_e_v_m[0];
            complex_near(leakage, expected, 2e-15 * expected.abs());
            let expected_power = 0.5 * bank.port_admittance_s()[0] * expected.norm_sq();
            near(
                response.power.outgoing_w_m2[0],
                expected_power,
                2e-15 * expected_power,
            );
            let dark = bank
                .evaluate([b.scale(scale), -a.scale(scale)], &mut [])
                .unwrap();
            assert_eq!(dark.outgoing_tangential_e_v_m[0], C64::ZERO);
            assert_eq!(dark.power.outgoing_w_m2[0], 0.0);
        }
    }
}

/// G1/G3: preserve a small electric-field residue and its physical positive
/// loss inside a lossy film. Independent exact algebra on the retained source
/// columns supplies the near-null field; the original Debye law supplies heat.
#[test]
fn coherent_film_retains_weak_node_field_and_local_heat() {
    let outside = material(E, M);
    let film = debye(E, 2.0 * E, 0.03, 0.001);
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
    let bank = stack.prepare_illumination(1.0, Polarization::Te).unwrap();
    let point = bank.prepare_point(0, 0.23 * layers[0].thickness_m).unwrap();
    let a = point.evaluate([C64::ONE, C64::ZERO]).unwrap().electric_v_m[1];
    let b = point.evaluate([C64::ZERO, C64::ONE]).unwrap().electric_v_m[1];
    let perturbed = C64::new(b.re.next_up(), b.im);
    let delta = perturbed.re - b.re;
    for exponent in [-200, 0, 200] {
        let scale = 2.0_f64.powi(exponent);
        let expected = a.scale(delta * scale);
        let field = point
            .evaluate([perturbed.scale(scale), -a.scale(scale)])
            .unwrap();
        complex_near(field.electric_v_m[1], expected, 2e-14 * expected.abs());
        // omega*tau=1, so epsilon''=-delta_epsilon/2=-E.
        let expected_heat = 0.5 * OMEGA * E * expected.norm_sq();
        near(
            field.statistics.absorbed_power_w_m3,
            expected_heat,
            5e-14 * expected_heat,
        );
        assert!(field.statistics.absorbed_power_w_m3 > 0.0);
        let dark = point.evaluate([b.scale(scale), -a.scale(scale)]).unwrap();
        assert_eq!(dark.electric_v_m, [C64::ZERO; 3]);
        assert_eq!(dark.statistics.absorbed_power_w_m3, 0.0);
    }
}

// Independent conventional characteristic matrices and face continuity.
// These modest depths permit growing sin/cos; production must never need them.
pub(super) fn matrix(
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
pub(super) fn product(a: [[C64; 2]; 2], b: [[C64; 2]; 2]) -> [[C64; 2]; 2] {
    std::array::from_fn(|i| std::array::from_fn(|j| a[i][0] * b[0][j] + a[i][1] * b[1][j]))
}
pub(super) fn act(a: [[C64; 2]; 2], x: [C64; 2]) -> [C64; 2] {
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

/// G1/G3: independent characteristic matrices resolve the physical work in
/// each cell of a standing-wave film. Repartitioning its depth does not alter
/// whole-film heat, and the left-source restriction agrees with the existing
/// one-sided field owner for TE and longitudinally lossy TM illumination.
#[test]
fn depth_intervals_preserve_coherent_interference_and_original_layer_heat() {
    let left = material(E, M);
    let right = material(4.0 * E, 1.2 * M);
    let passive = debye(3.0 * E, 4.0 * E, 0.03, 0.001);
    let transparent = material(2.0 * E, 0.7 * M);
    let thickness = 0.7 / k();
    let stack = PlanarDielectricLayers::try_new(
        &left,
        DielectricResponse::Lossless(&right),
        OMEGA,
        &[
            DielectricLayer {
                response: DielectricResponse::Debye(&passive),
                thickness_m: thickness,
            },
            DielectricLayer {
                response: DielectricResponse::Lossless(&transparent),
                thickness_m: 0.3 / k(),
            },
        ],
    )
    .unwrap();
    for cosine in [1.0, 0.72] {
        for pol in [Polarization::Te, Polarization::Tm] {
            let bank = stack.prepare_illumination(cosine, pol).unwrap();
            let single = stack.prepare_layer_field(0, cosine, pol).unwrap();
            let y = bank.port_admittance_s();
            let ratio = y[1] / y[0];
            let epsilon = C64::new(5.0, -2.0);
            let m0 = matrix(epsilon, 1.0, 0.7, cosine, pol);
            let m1 = matrix(C64::from_re(2.0), 0.7, 0.3, cosine, pol);
            let all = product(m0, m1);
            let plus = act(all, [C64::ONE, C64::from_re(ratio)]);
            let minus = act(all, [C64::ONE, C64::from_re(-ratio)]);
            for angle in [0.0, 0.6, 2.1, 3.14] {
                let input = [C64::new(0.4, 0.7), phase(angle).scale(0.8)];
                let transmitted =
                    (input[0].scale(2.0) - (minus[0] + minus[1]) * input[1]) / (plus[0] + plus[1]);
                let end_field = act(
                    m1,
                    [
                        transmitted + input[1],
                        (transmitted - input[1]).scale(ratio),
                    ],
                );
                let mut heat = [0.0; 2];
                bank.evaluate(input, &mut heat).unwrap();
                let mut sum = 0.0;
                let cuts = [0.0, 0.12, 0.51, 0.93, 1.0];
                for cut in cuts.windows(2) {
                    let start = cut[0] * thickness;
                    let end = cut[1] * thickness;
                    let interval = bank.prepare_depth_interval(0, start, end).unwrap();
                    assert!(std::ptr::eq(interval.illumination(), &bank));
                    assert_eq!(interval.layer_index(), 0);
                    assert_eq!(interval.start_depth_m(), start);
                    assert_eq!(interval.end_depth_m(), end);
                    let before = act(
                        matrix(epsilon, 1.0, 0.7 * (1.0 - cut[0]), cosine, pol),
                        end_field,
                    );
                    let after = act(
                        matrix(epsilon, 1.0, 0.7 * (1.0 - cut[1]), cosine, pol),
                        end_field,
                    );
                    let expected = 0.5
                        * y[0]
                        * ((before[0] * before[1].conj()).re - (after[0] * after[1].conj()).re);
                    let actual = interval.evaluate(input).unwrap();
                    assert!(actual > 0.0);
                    near(actual, expected, 3e-12 * heat[0]);
                    sum += actual;
                    near(
                        interval.evaluate([C64::ONE, C64::ZERO]).unwrap(),
                        single
                            .absorbed_power_in_depth_interval_w_m2(start, end)
                            .unwrap(),
                        3e-12 * single.absorbed_power_w_m2(),
                    );
                }
                near(sum, heat[0], 3e-12 * heat[0]);
                assert_eq!(
                    bank.prepare_depth_interval(0, 0.0, thickness)
                        .unwrap()
                        .evaluate(input)
                        .unwrap(),
                    heat[0]
                );
                assert_eq!(
                    bank.prepare_depth_interval(1, 0.01 / k(), 0.2 / k())
                        .unwrap()
                        .evaluate(input)
                        .unwrap(),
                    0.0
                );
            }
        }
    }
}

/// G1: a coherent electric node has cubic, positive interval heat even when
/// boundary-flux subtraction cannot resolve it. Independent face sources are
/// the phase average; coherence factors add powers without losing cross terms
/// inside each factor. Bounds and work admission happen before evaluation.
#[test]
fn depth_intervals_retain_thin_node_heat_partial_coherence_and_refusals() {
    let outside = material(E, M);
    let passive = debye(E, 2.0 * E, 0.03, 0.001);
    let thickness = 1e-7 / k();
    let stack = PlanarDielectricLayers::try_new(
        &outside,
        DielectricResponse::Lossless(&outside),
        OMEGA,
        &[
            DielectricLayer {
                response: DielectricResponse::Debye(&passive),
                thickness_m: thickness,
            },
            DielectricLayer {
                response: DielectricResponse::Debye(&passive),
                thickness_m: 0.0,
            },
        ],
    )
    .unwrap();
    for pol in [Polarization::Te, Polarization::Tm] {
        let bank = stack.prepare_illumination(1.0, pol).unwrap();
        let mut whole = [0.0; 2];
        bank.evaluate([C64::ONE, -C64::ONE], &mut whole).unwrap();
        assert!(whole[0] > 0.0);
        let mut sum = 0.0;
        for (quarter, expected_fraction) in [7.0 / 16.0, 1.0 / 16.0, 1.0 / 16.0, 7.0 / 16.0]
            .into_iter()
            .enumerate()
        {
            let interval = bank
                .prepare_depth_interval(
                    0,
                    quarter as f64 * thickness / 4.0,
                    (quarter + 1) as f64 * thickness / 4.0,
                )
                .unwrap();
            let heat = interval.evaluate([C64::ONE, -C64::ONE]).unwrap();
            assert!(heat > 0.0);
            near(heat, whole[0] * expected_fraction, 2e-8 * whole[0]);
            sum += heat;
        }
        near(sum, whole[0], 2e-8 * whole[0]);
        let interval = bank
            .prepare_depth_interval(0, 0.0, 0.25 * thickness)
            .unwrap();
        let independent = interval
            .evaluate_independent_sources([C64::ONE; 2])
            .unwrap();
        let phase_average = [C64::ONE, C64::new(0.0, 1.0), -C64::ONE, C64::new(0.0, -1.0)]
            .into_iter()
            .map(|source| interval.evaluate([C64::ONE, source]).unwrap())
            .sum::<f64>()
            / 4.0;
        near(independent, phase_average, 3e-14 * independent);
        let factors = [
            [C64::new(0.7, -0.2), C64::new(0.3, 0.4)],
            [C64::new(-0.1, 0.5), C64::new(0.2, -0.6)],
        ];
        let expected =
            interval.evaluate(factors[0]).unwrap() + interval.evaluate(factors[1]).unwrap();
        near(
            interval
                .evaluate_coherence_modes(&factors, BUDGET, || false)
                .unwrap(),
            expected,
            1e-15 * expected,
        );
        assert_eq!(
            interval
                .evaluate_coherence_modes(&[], BUDGET, || false)
                .unwrap(),
            0.0
        );
        assert_eq!(
            bank.prepare_depth_interval(0, thickness, thickness)
                .unwrap()
                .evaluate([C64::ONE; 2])
                .unwrap(),
            0.0
        );
        assert_eq!(
            bank.prepare_depth_interval(1, 0.0, 0.0)
                .unwrap()
                .evaluate([C64::ONE; 2])
                .unwrap(),
            0.0
        );
        assert!(bank.prepare_depth_interval(2, 0.0, 0.0).is_err());
        for bounds in [
            [-1.0, 0.0],
            [0.0, 2.0 * thickness],
            [thickness, 0.0],
            [f64::NAN, thickness],
            [0.0, f64::INFINITY],
        ] {
            assert!(
                bank.prepare_depth_interval(0, bounds[0], bounds[1])
                    .is_err()
            );
        }
        let mut calls = 0;
        assert!(
            interval
                .evaluate_coherence_modes(
                    &factors,
                    IlluminationBudget {
                        maximum_modes: 2,
                        maximum_layer_mode_pairs: 1
                    },
                    || {
                        calls += 1;
                        false
                    }
                )
                .is_err()
        );
        assert_eq!(calls, 0);
        assert!(
            interval
                .evaluate_coherence_modes(&factors, BUDGET, || {
                    calls += 1;
                    calls == 2
                })
                .is_err()
        );
        assert!(
            interval
                .evaluate([C64::new(f64::NAN, 0.0), C64::ZERO])
                .is_err()
        );
        assert!(
            interval
                .evaluate_coherence_modes(
                    &[[C64::ZERO, C64::new(0.0, f64::INFINITY)]],
                    BUDGET,
                    || false
                )
                .is_err()
        );
    }
}
