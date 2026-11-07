//! Factored complex-symmetric frequency contractions in original units.

use fs_math::{c64::C64, dd::Dd};

// Each reduced entry owns its product scale. A huge eliminated normal slope
// must not erase an unrelated weak in-plane slope. Keeping every product
// factored also permits cancellation before restoring physical units. DD is
// compensated accumulation, not a relative-accuracy theorem at cancellation.
#[derive(Clone, Copy, Default)]
pub(crate) struct Sum {
    fraction: Dd,
    exponent: i32,
}

fn split(value: f64) -> (f64, i32) {
    let (value, correction) = if value.abs() < f64::MIN_POSITIVE {
        (value * 18_014_398_509_481_984.0, -54)
    } else {
        (value, 0)
    };
    let bits = value.to_bits();
    let exponent = ((bits >> 52) & 0x7ff) as i32 - 1023 + correction;
    let fraction = f64::from_bits((bits & 0x800f_ffff_ffff_ffff) | (1023_u64 << 52));
    (fraction, exponent)
}

fn binary_scale(mut value: f64, mut exponent: i32) -> f64 {
    while exponent > 512 {
        value *= f64::from_bits((1535_u64) << 52);
        exponent -= 512;
    }
    while exponent < -512 {
        value *= f64::from_bits((511_u64) << 52);
        exponent += 512;
    }
    value * f64::from_bits(((exponent + 1023) as u64) << 52)
}

impl Sum {
    pub(crate) fn add(&mut self, factors: [f64; 4]) {
        self.add_ratio(&factors, &[]);
    }

    // Inputs are original finite factors; denominators must be nonzero.
    // Dividing their binary mantissas avoids overflowing 1/tau before the
    // physical spring strain has supplied its compensating scale.
    pub(crate) fn add_ratio(&mut self, numerator: &[f64], denominator: &[f64]) {
        self.add_binary_ratio(numerator, denominator, 0);
    }

    fn add_binary_ratio(&mut self, numerator: &[f64], denominator: &[f64], shift: i32) {
        if numerator.contains(&0.0) {
            return;
        }
        let (mut fraction, mut exponent) = (Dd::ONE, shift);
        for &factor in numerator {
            let (f, e) = split(factor);
            fraction = fraction * Dd::from_f64(f);
            exponent += e;
        }
        for &factor in denominator {
            let (f, e) = split(factor);
            fraction = fraction / Dd::from_f64(f);
            exponent -= e;
        }
        self.accumulate(fraction, exponent);
    }

    fn accumulate(&mut self, fraction: Dd, exponent: i32) {
        if self.fraction == Dd::ZERO {
            self.fraction = fraction;
            self.exponent = exponent;
            return;
        }
        let scale = self.exponent.max(exponent);
        let scaled = |value: Dd, shift| {
            Dd::from_pair(binary_scale(value.hi, shift), binary_scale(value.lo, shift))
        };
        self.fraction =
            scaled(self.fraction, self.exponent - scale) + scaled(fraction, exponent - scale);
        self.exponent = scale;
    }

    // Square only normalized channel mantissas. Both the source spectrum and
    // original tau remain binary-factored until this physical energy/heat sum.
    fn add_square(
        &mut self,
        channel: Self,
        spectrum: (f64, i32),
        weight: f64,
        denominator: &[f64],
    ) {
        if channel.fraction == Dd::ZERO || spectrum.0 == 0.0 {
            return;
        }
        let (_, shift) = split(if channel.fraction.hi != 0.0 {
            channel.fraction.hi
        } else {
            channel.fraction.lo
        });
        let normalized = Dd::from_pair(
            binary_scale(channel.fraction.hi, -shift),
            binary_scale(channel.fraction.lo, -shift),
        );
        let mut factor = Self::default();
        factor.add_ratio(&[weight, spectrum.0], denominator);
        self.accumulate(
            normalized * normalized * factor.fraction,
            2 * (channel.exponent + shift) + spectrum.1 + factor.exponent,
        );
    }

    pub(crate) fn value(self) -> f64 {
        binary_scale(self.fraction.to_f64(), self.exponent)
    }
}

pub(crate) fn apply<const N: usize, const M: usize>(
    t: &[[C64; M]; N],
    strain: &[C64; M],
) -> [C64; N] {
    core::array::from_fn(|i| {
        let (mut real, mut imaginary) = (Sum::default(), Sum::default());
        for (map, value) in t[i].iter().zip(strain) {
            real.add([map.re, value.re, 1.0, 1.0]);
            real.add([-map.im, value.im, 1.0, 1.0]);
            imaginary.add([map.re, value.im, 1.0, 1.0]);
            imaginary.add([map.im, value.re, 1.0, 1.0]);
        }
        C64::new(real.value(), imaginary.value())
    })
}

/// Apply the original complex branch transfer (storage - i loss) without
/// rounding its intermediate spring strain or any matrix entry first.
pub(crate) fn apply_spectrum<const N: usize, const M: usize>(
    matrix: &[[f64; M]; N],
    vector: &[C64; M],
    spectrum: [(f64, i32); 2],
) -> [C64; N] {
    let [a, b] = spectrum;
    core::array::from_fn(|i| {
        let (mut re, mut im) = (Sum::default(), Sum::default());
        for (&map, value) in matrix[i].iter().zip(vector) {
            re.add_binary_ratio(&[map, value.re, a.0], &[], a.1);
            re.add_binary_ratio(&[map, value.im, b.0], &[], b.1);
            im.add_binary_ratio(&[map, value.im, a.0], &[], a.1);
            im.add_binary_ratio(&[-map, value.re, b.0], &[], b.1);
        }
        C64::new(re.value(), im.value())
    })
}

/// Positive Gram norm with the complete original spectral/rate/heat scale.
/// No channel or squared strain is rounded into physical units prematurely.
pub(crate) fn weighted_norm<const N: usize, const M: usize>(
    matrix: &[[f64; M]; N],
    vector: &[C64; M],
    rate: f64,
    spectrum: (f64, i32),
    weight: f64,
    denominator: &[f64],
) -> f64 {
    weighted_norms(matrix, vector, rate, spectrum, [(weight, denominator)])[0]
}

/// Reuse each original Gram contraction for independent storage/heat scales.
pub(crate) fn weighted_norms<const N: usize, const M: usize, const K: usize>(
    matrix: &[[f64; M]; N],
    vector: &[C64; M],
    rate: f64,
    spectrum: (f64, i32),
    weights: [(f64, &[f64]); K],
) -> [f64; K] {
    let mut totals = [Sum::default(); K];
    for row in matrix {
        let (mut re, mut im) = (Sum::default(), Sum::default());
        for (&map, value) in row.iter().zip(vector) {
            re.add([map, value.re, rate, 1.0]);
            im.add([map, value.im, rate, 1.0]);
        }
        for (total, (weight, denominator)) in totals.iter_mut().zip(&weights) {
            total.add_square(re, spectrum, *weight, denominator);
            total.add_square(im, spectrum, *weight, denominator);
        }
    }
    totals.map(Sum::value)
}

pub(crate) fn project<const M: usize>(
    c: &[[C64; 6]; 6],
    t: &[[C64; M]; 6],
    i: usize,
    j: usize,
    scale: f64,
) -> C64 {
    let (mut real, mut imaginary) = (Sum::default(), Sum::default());
    for row in 0..6 {
        for column in 0..6 {
            let (a, b, d) = (t[row][i], c[row][column], t[column][j]);
            for factors in [
                [a.re, b.re, d.re, scale],
                [-a.im, b.im, d.re, scale],
                [-a.re, b.im, d.im, scale],
                [-a.im, b.re, d.im, scale],
            ] {
                real.add(factors);
            }
            for factors in [
                [a.re, b.re, d.im, scale],
                [a.im, b.re, d.re, scale],
                [a.re, b.im, d.re, scale],
                [-a.im, b.im, d.im, scale],
            ] {
                imaginary.add(factors);
            }
        }
    }
    C64::new(real.value(), imaginary.value())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn physical_stress_and_heat_survive_an_unrepresentable_spring_strain() {
        let branch = [(1.0, -2000), (1.0, -1000)];
        let strain = [C64::from_re(2_f64.powi(-200))];
        let stiffness = [[2_f64.powi(400)]];
        let gram = [[2_f64.powi(200)]];
        assert_eq!(apply_spectrum(&[[1.0]], &strain, branch)[0], C64::ZERO);
        assert_eq!(
            apply_spectrum(&stiffness, &strain, branch)[0].im,
            -2_f64.powi(-800)
        );
        assert_eq!(
            weighted_norm(&gram, &strain, 1.0, branch[0], 0.5, &[2_f64.powi(-1000)]),
            2_f64.powi(-1001)
        );
        // The opposite limit: B*epsilon overflows, but the complete
        // branch energy remains finite after original frequency scaling.
        assert_eq!(
            weighted_norm(
                &[[2_f64.powi(600)]],
                &[C64::from_re(2_f64.powi(500))],
                1.0,
                (1.0, -2000),
                0.25,
                &[]
            ),
            2_f64.powi(198)
        );
    }

    #[test]
    fn original_relaxation_time_divides_after_physical_spring_scale() {
        let mut heat = Sum::default();
        heat.add_ratio(&[0.5, 1e-200, 1e-200], &[1e-300]);
        assert!((heat.value() - 5e-101).abs() < 2e-115);
        let mut subnormal_tau = Sum::default();
        let spring = 2_f64.powi(-537);
        subnormal_tau.add_ratio(&[0.5, spring, spring], &[f64::from_bits(1)]);
        assert_eq!(subnormal_tau.value(), 0.5);
    }

    #[test]
    fn physical_rate_and_gram_are_applied_before_restoring_units() {
        let result = apply_spectrum(
            &[[1e200, 0.0], [0.0, -1e200]],
            &[C64::new(1e200, -1e200), C64::new(1e200, 1e200)],
            [(1e-300, 0), (0.0, 0)],
        );
        for (actual, expected) in result
            .into_iter()
            .zip([C64::new(1e100, -1e100), C64::new(-1e100, -1e100)])
        {
            assert!((actual.re - expected.re).abs() < 3e-16 * expected.re.abs());
            assert!((actual.im - expected.im).abs() < 3e-16 * expected.im.abs());
        }
    }

    #[test]
    fn full_strain_map_retains_signed_complex_cancellation() {
        let mut map = [[C64::ZERO; 3]; 6];
        map[0] = [
            C64::new(1e200, 1e200),
            C64::new(-1e200, -1e200),
            C64::new(3.0, -2.0),
        ];
        map[1][2] = C64::from_re(f64::from_bits(1));
        let strain = [
            C64::new(1e200, -1e200),
            C64::new(1e200, -1e200),
            C64::new(2.0, 1.0),
        ];
        let actual = apply(&map, &strain);
        assert_eq!(actual[0], C64::new(8.0, -1.0));
        assert_eq!(actual[1], C64::new(f64::from_bits(2), f64::from_bits(1)));
        assert!(actual[2..].iter().all(|v| *v == C64::ZERO));
    }

    #[test]
    fn complex_symmetric_projection_matches_direct_full_contraction() {
        let c = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                C64::new((i + j + 1) as f64 * 0.17, (i * j + 2) as f64 * -0.13)
            })
        });
        let t: [[C64; 3]; 6] = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                C64::new((i + 2 * j + 1) as f64 * 0.11, (2 * i + j + 1) as f64 * 0.07)
            })
        });
        for i in 0..3 {
            for j in 0..3 {
                let mut expected = C64::ZERO;
                for row in 0..6 {
                    for column in 0..6 {
                        expected =
                            expected + (t[row][i] * c[row][column] * t[column][j]).scale(0.5);
                    }
                }
                let actual = project(&c, &t, i, j, 0.5);
                assert!((actual - expected).abs() <= 3e-14 * expected.abs());
            }
        }
    }

    #[test]
    fn unrelated_huge_transverse_slope_preserves_weak_retained_response() {
        let mut c = [[C64::ZERO; 6]; 6];
        c[0][0] = C64::new(0.0, 1e-100);
        c[2][2] = C64::new(0.0, 1e300);
        let mut t = [[C64::ZERO; 3]; 6];
        t[0][0] = C64::ONE;
        t[1][1] = C64::ONE;
        t[3][2] = C64::ONE;
        assert_eq!(project(&c, &t, 0, 0, 1.0), c[0][0]);
        assert_eq!(project(&c, &t, 1, 1, 1.0), C64::ZERO);
    }

    #[test]
    fn unrepresentable_products_cancel_before_restoring_physical_units() {
        let mut c = [[C64::ZERO; 6]; 6];
        c[0][0] = C64::from_re(1e200);
        c[1][1] = C64::from_re(-1e200);
        c[2][2] = C64::new(3.0, -2.0);
        let mut t = [[C64::ZERO; 3]; 6];
        t[0][0] = C64::from_re(1e100);
        t[1][0] = C64::from_re(1e100);
        t[2][0] = C64::ONE;
        assert_eq!(project(&c, &t, 0, 0, 1.0), C64::new(3.0, -2.0));
    }

    #[test]
    fn minimum_subnormal_coefficient_survives_physical_constraint_scaling() {
        let mut c = [[C64::ZERO; 6]; 6];
        c[0][0] = C64::new(f64::from_bits(1), -f64::from_bits(1));
        let mut t = [[C64::ZERO; 3]; 6];
        t[0][0] = C64::from_re(2.0_f64.powi(500));
        assert_eq!(
            project(&c, &t, 0, 0, 0.5),
            C64::new(2.0_f64.powi(-75), -2.0_f64.powi(-75)),
        );
    }
}
