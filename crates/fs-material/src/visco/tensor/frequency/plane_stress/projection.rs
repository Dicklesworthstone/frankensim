//! Factored complex-symmetric frequency contractions in original units.

use fs_math::{c64::C64, dd::Dd};

// Each reduced entry owns its product scale. A huge eliminated normal slope
// must not erase an unrelated weak in-plane slope. Keeping every product
// factored also permits cancellation before restoring physical units. DD is
// compensated accumulation, not a relative-accuracy theorem at cancellation.
#[derive(Default)]
struct Sum {
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
    fn add(&mut self, factors: [f64; 4]) {
        if factors.contains(&0.0) {
            return;
        }
        let (mut fraction, mut exponent) = (Dd::ONE, 0);
        for factor in factors {
            let (f, e) = split(factor);
            fraction = fraction * Dd::from_f64(f);
            exponent += e;
        }
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

    fn value(self) -> f64 {
        binary_scale(self.fraction.to_f64(), self.exponent)
    }
}

pub(super) fn project(c: &[[C64; 6]; 6], t: &[[C64; 3]; 6], i: usize, j: usize, scale: f64) -> C64 {
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
    fn complex_symmetric_projection_matches_direct_full_contraction() {
        let c = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                C64::new((i + j + 1) as f64 * 0.17, (i * j + 2) as f64 * -0.13)
            })
        });
        let t = core::array::from_fn(|i| {
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
