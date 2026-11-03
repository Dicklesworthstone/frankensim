//! Cycle-mean momentum balance between two homogeneous lossless exteriors.
//!
//! For a planar film, tangential momentum transfer is k_parallel*A/omega;
//! normal transfer under left incidence is
//! [k_z_left*(2R+A) + (k_z_left-k_z_right)*T]/omega. R, T and A are actual
//! reflected/transmitted/absorbed powers, not dimensionless fractions. The
//! original positive physical loss and a separately resolved wavevector
//! difference retain weak absorption and refraction when direct subtraction
//! of incident/outgoing momentum would erase them. This is the exterior
//! Maxwell stress, without medium striction.
//! Reference: MIT CETA, radiation pressure on a slab at arbitrary incidence:
//! <https://web.mit.edu/ceta/obt/fund-media.html>.

use fs_math::{c64::C64, dd::Dd};

// Resolve (k_z_left-k_z_right)/omega from the original lossless material
// ratios. Rationalization and one fused product retain index contrasts that
// disappear from either epsilon_ratio*mu_ratio or the rounded normal roots.
// The original wave owner admits these inputs and the outward branch of q.
pub(super) fn normal_difference(
    incident_inverse_speed: f64,
    incidence_cosine: f64,
    exit_normal_ratio: C64,
    epsilon_ratio: f64,
    mu_ratio: f64,
) -> f64 {
    if exit_normal_ratio.im != 0.0 {
        // A lossless evanescent exit carries no normal transmitted power.
        // Its real normal wavevector is zero; never divide by its complex root.
        incident_inverse_speed * incidence_cosine
    } else {
        incident_inverse_speed
            * ((-epsilon_ratio).mul_add(mu_ratio, 1.0) / (incidence_cosine + exit_normal_ratio.re))
    }
}

pub(super) fn one_sided(
    wavevector_over_omega_s_m: [f64; 3],
    normal_difference: f64,
    reflected: f64,
    transmitted: f64,
    absorbed: f64,
) -> [f64; 3] {
    let [kx, ky, kz] = wavevector_over_omega_s_m;
    [
        kx * absorbed,
        ky * absorbed,
        sum_terms([
            kz * reflected * 2.0,
            kz * absorbed,
            normal_difference * transmitted,
        ]),
    ]
}

// Actual normal force, with each exterior's positive k_z/omega and admittance.
// Separate reflection, original physical absorption, refraction and coherent
// interference. normal_difference is resolved by the material/geometry owner;
// reconstructing it from the rounded normal_over_omega values can lose weak
// contrast. Absorption diagonals are powers per squared incident amplitude.
pub(super) fn two_sided_normal(
    input: [C64; 2],
    scattering: [[C64; 2]; 2],
    admittance: [f64; 2],
    normal_over_omega: [f64; 2],
    normal_difference: f64,
    absorption_diagonal: [f64; 2],
) -> f64 {
    let reflected = [scattering[0][0] * input[0], scattering[1][1] * input[1]];
    let transmitted = [scattering[0][1] * input[1], scattering[1][0] * input[0]];
    let cross_left = (reflected[0] * transmitted[0].conj()).re;
    let cross_right = (transmitted[1] * reflected[1].conj()).re;
    let [left, right] = normal_over_omega;
    let [y_left, y_right] = admittance;
    sum_terms([
        left * y_left * reflected[0].norm_sq(),
        -right * y_right * reflected[1].norm_sq(),
        left * absorption_diagonal[0] * input[0].norm_sq(),
        -right * absorption_diagonal[1] * input[1].norm_sq(),
        normal_difference * 0.5 * y_right * transmitted[1].norm_sq(),
        normal_difference * 0.5 * y_left * transmitted[0].norm_sq(),
        left * y_left * cross_left,
        -right * y_right * cross_right,
    ])
}

fn sum_terms<const N: usize>(terms: [f64; N]) -> f64 {
    terms
        .into_iter()
        .fold(Dd::ZERO, |s, x| s + Dd::from_f64(x))
        .to_f64()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normal_difference_retains_weak_contrast_and_admits_critical_and_tir() {
        let increment = f64::EPSILON;
        let epsilon = 1.0 + increment;
        let mu = 1.0 - increment;
        assert_eq!(epsilon * mu, 1.0);
        assert_eq!(
            normal_difference(1.0, 1.0, C64::ONE, epsilon, mu),
            0.5 * increment * increment,
        );
        for cosine in [0.1, 0.5, 1.0] {
            assert_eq!(
                normal_difference(2.0, cosine, C64::from_re(cosine), 1.0, 1.0),
                0.0
            );
            assert_eq!(
                normal_difference(2.0, cosine, C64::from_re(cosine), 2.0, 0.5),
                0.0
            );
        }
        assert_eq!(normal_difference(2.0, 0.5, C64::ZERO, 0.75, 1.0), 1.0);
        assert_eq!(
            normal_difference(2.0, 0.5, C64::new(0.0, 0.5_f64.sqrt()), 0.25, 1.0),
            1.0,
        );
    }

    #[test]
    fn impedance_matched_refraction_pulls_toward_lower_index_from_either_side() {
        // epsilon=mu=1 on the left and epsilon=mu=2 on the right: admittance
        // is one on both sides, with k_z/omega one and two. There is no reflected
        // field or absorption. Unit peak field carries power 1/2 on either side.
        let transparent = [[C64::ZERO, C64::ONE], [C64::ONE, C64::ZERO]];
        let difference = normal_difference(1.0, 1.0, C64::from_re(2.0), 2.0, 2.0);
        assert_eq!(difference, -1.0);
        assert_eq!(
            one_sided([0.0, 0.0, 1.0], difference, 0.0, 0.5, 0.0)[2],
            -0.5
        );
        for input in [[C64::ONE, C64::ZERO], [C64::ZERO, C64::ONE]] {
            assert_eq!(
                two_sided_normal(
                    input,
                    transparent,
                    [1.0; 2],
                    [1.0, 2.0],
                    difference,
                    [0.0; 2]
                ),
                -0.5,
            );
        }
        for phase in [C64::ONE, C64::new(0.0, 1.0), C64::from_re(-1.0)] {
            assert_eq!(
                two_sided_normal(
                    [C64::ONE, phase],
                    transparent,
                    [1.0; 2],
                    [1.0, 2.0],
                    difference,
                    [0.0; 2]
                ),
                -1.0,
            );
        }
    }

    #[test]
    fn weak_absorption_and_oblique_mirror_momentum_remain_distinct() {
        let k = [0.6, -0.3, (0.55_f64).sqrt()];
        assert_eq!(one_sided(k, 0.0, 0.0, 2.0, 0.0), [0.0; 3]);
        let mirror = one_sided(k, -0.4, 2.0, 0.0, 0.0);
        assert_eq!(mirror, [0.0, -0.0, 4.0 * k[2]]);
        let weak = one_sided(k, 0.0, 0.0, 1.0, 1e-28);
        assert_eq!(weak, k.map(|x| x * 1e-28));
        let input = [C64::ONE, C64::ZERO];
        let transparent = [[C64::ZERO, C64::ONE], [C64::ONE, C64::ZERO]];
        assert_eq!(
            two_sided_normal(input, transparent, [1.0; 2], [k[2]; 2], 0.0, [1e-28; 2]),
            k[2] * 1e-28,
        );
        // Reflection and refractive momentum cancel, leaving physical loss
        // far below either term's ulp. A left-to-right f64 sum loses this force.
        assert_eq!(one_sided([0.0, 0.0, 1.0], -1.0, 0.5, 1.0, 1e-28)[2], 1e-28);
    }

    #[test]
    fn weak_refractive_contrast_survives_rounded_exterior_wavevectors() {
        let transparent = [[C64::ZERO, C64::ONE], [C64::ONE, C64::ZERO]];
        for difference in [-1e-28, 1e-28] {
            // The caller has resolved this difference from its original
            // constitutive data even though both rounded k_z/omega equal one.
            assert_eq!(
                one_sided([0.2, -0.3, 1.0], difference, 0.0, 2.0, 0.0),
                [0.0, -0.0, 2.0 * difference],
            );
            for input in [[C64::ONE, C64::ZERO], [C64::ZERO, C64::ONE]] {
                assert_eq!(
                    two_sided_normal(input, transparent, [1.0; 2], [1.0; 2], difference, [0.0; 2]),
                    0.5 * difference,
                );
            }
        }
    }

    #[test]
    fn coherent_normal_force_matches_independent_unequal_exterior_flux() {
        // Build reciprocal passive scattering in power-normalized coordinates,
        // then convert to field amplitudes with the two original admittances.
        // The first pair is exactly lossless; the second has physical loss.
        for eigenvalues in [
            [C64::ONE, C64::new(0.0, 1.0)],
            [C64::new(0.2, 0.1), C64::new(-0.4, 0.3)],
        ] {
            let r = (eigenvalues[0] + eigenvalues[1]).scale(0.5);
            let t = (eigenvalues[0] - eigenvalues[1]).scale(0.5);
            for (y, k) in [([1.7_f64, 0.6_f64], [0.7, 1.3]), ([0.8, 2.3], [1.6, 0.4])] {
                let s = [
                    [r, t.scale((y[1] / y[0]).sqrt())],
                    [t.scale((y[0] / y[1]).sqrt()), r],
                ];
                let loss = y.map(|admittance| 0.5 * admittance * (1.0 - r.norm_sq() - t.norm_sq()));
                for input in [
                    [C64::ONE, C64::ZERO],
                    [C64::ZERO, C64::ONE],
                    [C64::new(0.2, 0.7), C64::new(-0.4, 0.8)],
                    [C64::ONE, C64::ONE],
                    [C64::ONE, C64::new(0.0, 1.0)],
                ] {
                    let output = s.map(|row| row[0] * input[0] + row[1] * input[1]);
                    // Independent exterior stress: incoming and outgoing waves
                    // both press on their own side, with distinct momentum/power.
                    let left = 0.5 * k[0] * y[0] * (input[0].norm_sq() + output[0].norm_sq());
                    let right = 0.5 * k[1] * y[1] * (input[1].norm_sq() + output[1].norm_sq());
                    let actual = two_sided_normal(input, s, y, k, k[0] - k[1], loss);
                    assert!((actual - (left - right)).abs() < 2e-15 * (left + right).max(1.0));
                    if input[1] == C64::ZERO {
                        let force = one_sided(
                            [0.2, -0.3, k[0]],
                            k[0] - k[1],
                            0.5 * y[0] * output[0].norm_sq(),
                            0.5 * y[1] * output[1].norm_sq(),
                            loss[0] * input[0].norm_sq(),
                        );
                        assert!((force[2] - actual).abs() < 2e-15);
                        assert_eq!(force[0], 0.2 * loss[0] * input[0].norm_sq());
                        assert_eq!(force[1], -0.3 * loss[0] * input[0].norm_sq());
                    }
                }
            }
        }
    }

    #[test]
    fn unequal_exterior_mirrors_use_each_incident_momentum() {
        let input = [C64::new(0.2, 0.7), C64::new(-0.4, 0.8)];
        let mirror = [[C64::from_re(-1.0), C64::ZERO], [C64::ZERO, C64::ONE]];
        let y = [1.7, 0.6];
        let k = [0.7, 1.3];
        let expected = k[0] * y[0] * input[0].norm_sq() - k[1] * y[1] * input[1].norm_sq();
        let actual = two_sided_normal(input, mirror, y, k, k[0] - k[1], [0.0; 2]);
        assert!((actual - expected).abs() < 2e-15);
    }

    #[test]
    fn matched_exteriors_preserve_coherent_force_and_reversal() {
        // Symmetric passive slab: the two scattering eigenvalues have modulus
        // <=1. This exercises nonzero phase-sensitive interference as well as
        // independent left/right illumination and unequal coherent amplitudes.
        for eigenvalues in [
            [C64::new(0.6, 0.8), C64::new(-0.8, 0.6)],
            [C64::new(0.2, 0.1), C64::new(-0.4, 0.3)],
        ] {
            let r = (eigenvalues[0] + eigenvalues[1]).scale(0.5);
            let t = (eigenvalues[0] - eigenvalues[1]).scale(0.5);
            let s = [[r, t], [t, r]];
            let y = 1.7;
            let normal = 0.73;
            let loss = 0.5 * y * (1.0 - r.norm_sq() - t.norm_sq());
            for input in [
                [C64::ONE, C64::ZERO],
                [C64::ZERO, C64::ONE],
                [C64::new(0.2, 0.7), C64::new(-0.4, 0.8)],
                [C64::ONE, C64::new(0.0, 1.0)],
            ] {
                let output = s.map(|row| row[0] * input[0] + row[1] * input[1]);
                let expected = 0.5
                    * normal
                    * y
                    * (input[0].norm_sq() + output[0].norm_sq()
                        - input[1].norm_sq()
                        - output[1].norm_sq());
                let actual = two_sided_normal(input, s, [y; 2], [normal; 2], 0.0, [loss; 2]);
                assert!((actual - expected).abs() < 2e-15);
                let reversed =
                    two_sided_normal([input[1], input[0]], s, [y; 2], [normal; 2], 0.0, [loss; 2]);
                assert!((actual + reversed).abs() < 2e-15);
            }
        }
    }
}
