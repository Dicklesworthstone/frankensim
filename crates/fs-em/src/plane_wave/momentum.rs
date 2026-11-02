//! Cycle-mean momentum balance in one homogeneous lossless exterior.
//!
//! For a planar film, tangential momentum transfer is k_parallel*A/omega;
//! normal transfer under left incidence is k_z*(2R+A)/omega. Here R and A
//! are actual reflected/absorbed powers, not dimensionless fractions. Using
//! the positive physical loss retains weak absorption even when 1-R-T rounds
//! to zero. This is the exterior Maxwell stress, without medium striction.
//! Reference: MIT CETA, radiation pressure on a slab at arbitrary incidence:
//! <https://web.mit.edu/ceta/obt/fund-media.html>.

use fs_math::{c64::C64, dd::Dd};

pub(super) fn one_sided(
    wavevector_over_omega_s_m: [f64; 3],
    reflected: f64,
    absorbed: f64,
) -> [f64; 3] {
    let [kx, ky, kz] = wavevector_over_omega_s_m;
    [
        kx * absorbed,
        ky * absorbed,
        kz * reflected * 2.0 + kz * absorbed,
    ]
}

// Normal momentum before multiplication by k_z/omega. Separate the two
// positive diagonal responses from their coherent interference: subtracting
// almost equal transmitted/incident powers would erase weak one-sided force.
pub(super) fn two_sided_normal(
    input: [C64; 2],
    scattering: [[C64; 2]; 2],
    admittance: f64,
    absorption_diagonal: [f64; 2],
) -> f64 {
    let reflected = [scattering[0][0] * input[0], scattering[1][1] * input[1]];
    let transmitted = [scattering[0][1] * input[1], scattering[1][0] * input[0]];
    let cross_left = (reflected[0] * transmitted[0].conj()).re;
    let cross_right = (transmitted[1] * reflected[1].conj()).re;
    let terms = [
        admittance * reflected[0].norm_sq(),
        -admittance * reflected[1].norm_sq(),
        absorption_diagonal[0] * input[0].norm_sq(),
        -absorption_diagonal[1] * input[1].norm_sq(),
        admittance * cross_left,
        -admittance * cross_right,
    ];
    terms
        .into_iter()
        .fold(Dd::ZERO, |s, x| s + Dd::from_f64(x))
        .to_f64()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn weak_absorption_and_oblique_mirror_momentum_remain_distinct() {
        let k = [0.6, -0.3, (0.55_f64).sqrt()];
        assert_eq!(one_sided(k, 0.0, 0.0), [0.0; 3]);
        let mirror = one_sided(k, 2.0, 0.0);
        assert_eq!(mirror, [0.0, -0.0, 4.0 * k[2]]);
        let weak = one_sided(k, 0.0, 1e-28);
        assert_eq!(weak, k.map(|x| x * 1e-28));
        let input = [C64::ONE, C64::ZERO];
        let transparent = [[C64::ZERO, C64::ONE], [C64::ONE, C64::ZERO]];
        assert_eq!(two_sided_normal(input, transparent, 1.0, [1e-28; 2]), 1e-28);
    }

    #[test]
    fn coherent_normal_force_matches_independent_exterior_flux() {
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
            let loss = 0.5 * y * (1.0 - r.norm_sq() - t.norm_sq());
            for input in [
                [C64::ONE, C64::ZERO],
                [C64::ZERO, C64::ONE],
                [C64::new(0.2, 0.7), C64::new(-0.4, 0.8)],
                [C64::ONE, C64::new(0.0, 1.0)],
            ] {
                let output = s.map(|row| row[0] * input[0] + row[1] * input[1]);
                let expected = 0.5
                    * y
                    * (input[0].norm_sq() + output[0].norm_sq()
                        - input[1].norm_sq()
                        - output[1].norm_sq());
                let actual = two_sided_normal(input, s, y, [loss; 2]);
                assert!((actual - expected).abs() < 2e-15);
                let reversed = two_sided_normal([input[1], input[0]], s, y, [loss; 2]);
                assert!((actual + reversed).abs() < 2e-15);
            }
        }
    }
}
