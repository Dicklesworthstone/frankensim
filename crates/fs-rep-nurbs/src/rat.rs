//! Exact rational arithmetic over i128 — the scalar the EXACT spline
//! algebra runs in. Every operation is gcd-reduced and overflow-checked:
//! leaving the exactness domain is a structured event (checked panics
//! with a named message, bounded by construction in the conformance
//! fixtures), never silent wraparound.

use core::cmp::Ordering;

/// A reduced fraction (denominator > 0).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Rat {
    /// Numerator.
    num: i128,
    /// Denominator (always positive).
    den: i128,
}

fn gcd_unsigned(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a.max(1)
}

fn gcd(a: i128, b: i128) -> i128 {
    i128::try_from(gcd_unsigned(a.unsigned_abs(), b.unsigned_abs()))
        .expect("Rat gcd magnitude exceeds i128 after denominator-bounded reduction")
}

fn signed_magnitude(magnitude: u128, negative: bool, op: &str) -> i128 {
    if negative && magnitude == (1u128 << 127) {
        return i128::MIN;
    }
    let positive = i128::try_from(magnitude)
        .unwrap_or_else(|_| panic!("Rat {op}: i128 overflow — exactness domain exceeded"));
    if negative {
        positive
            .checked_neg()
            .unwrap_or_else(|| panic!("Rat {op}: i128 overflow — exactness domain exceeded"))
    } else {
        positive
    }
}

fn cmp_unsigned_fraction(
    mut left_num: u128,
    mut left_den: u128,
    mut right_num: u128,
    mut right_den: u128,
) -> Ordering {
    let mut reversed = false;
    loop {
        let (left_whole, left_rem) = (left_num / left_den, left_num % left_den);
        let (right_whole, right_rem) = (right_num / right_den, right_num % right_den);
        let whole_order = left_whole.cmp(&right_whole);
        if whole_order != Ordering::Equal {
            return if reversed {
                whole_order.reverse()
            } else {
                whole_order
            };
        }
        let terminal_order = match (left_rem == 0, right_rem == 0) {
            (true, true) => return Ordering::Equal,
            (true, false) => Some(Ordering::Less),
            (false, true) => Some(Ordering::Greater),
            (false, false) => None,
        };
        if let Some(order) = terminal_order {
            return if reversed { order.reverse() } else { order };
        }
        (left_num, left_den) = (left_den, left_rem);
        (right_num, right_den) = (right_den, right_rem);
        reversed = !reversed;
    }
}

impl Rat {
    /// Construct and reduce.
    ///
    /// # Panics
    /// On a zero denominator (structural misuse, not data).
    #[must_use]
    pub fn new(num: i128, den: i128) -> Rat {
        assert!(den != 0, "Rat with zero denominator");
        let negative = (num < 0) ^ (den < 0);
        let g = gcd_unsigned(num.unsigned_abs(), den.unsigned_abs());
        let reduced_num = num.unsigned_abs() / g;
        let reduced_den = den.unsigned_abs() / g;
        Rat {
            num: signed_magnitude(reduced_num, negative, "construct"),
            den: i128::try_from(reduced_den).unwrap_or_else(|_| {
                panic!("Rat construct: i128 overflow — exactness domain exceeded")
            }),
        }
    }

    /// From an integer.
    #[must_use]
    pub fn int(v: i64) -> Rat {
        Rat {
            num: i128::from(v),
            den: 1,
        }
    }

    /// Numerator (reduced).
    #[must_use]
    pub fn numerator(self) -> i128 {
        self.num
    }

    /// Denominator (reduced, positive).
    #[must_use]
    pub fn denominator(self) -> i128 {
        self.den
    }

    /// Nearest f64 (for reporting/plotting only — never for exactness).
    #[must_use]
    #[allow(clippy::cast_precision_loss)]
    pub fn to_f64(self) -> f64 {
        self.num as f64 / self.den as f64
    }

    fn checked(num: Option<i128>, den: Option<i128>, op: &str) -> Rat {
        let (Some(n), Some(d)) = (num, den) else {
            panic!("Rat {op}: i128 overflow — exactness domain exceeded");
        };
        Rat::new(n, d)
    }
}

impl core::ops::Add for Rat {
    type Output = Rat;
    fn add(self, o: Rat) -> Rat {
        // a/b + c/d over lcm. The provisional numerator can share a second
        // factor with the common-denominator gcd; cancel it before forming the
        // final denominator so a representable result does not falsely
        // overflow merely because the unreduced lcm does.
        let g = gcd(self.den, o.den);
        let (db, dd) = (self.den / g, o.den / g);
        let numerator = self
            .num
            .checked_mul(dd)
            .and_then(|left| {
                o.num
                    .checked_mul(db)
                    .and_then(|right| left.checked_add(right))
            })
            .unwrap_or_else(|| panic!("Rat add: i128 overflow — exactness domain exceeded"));
        let post_gcd = gcd(numerator, g);
        Rat::checked(
            Some(numerator / post_gcd),
            (self.den / post_gcd).checked_mul(dd),
            "add",
        )
    }
}

impl core::ops::Sub for Rat {
    type Output = Rat;
    fn sub(self, o: Rat) -> Rat {
        let g = gcd(self.den, o.den);
        let (db, dd) = (self.den / g, o.den / g);
        let numerator = self
            .num
            .checked_mul(dd)
            .and_then(|left| {
                o.num
                    .checked_mul(db)
                    .and_then(|right| left.checked_sub(right))
            })
            .unwrap_or_else(|| panic!("Rat sub: i128 overflow — exactness domain exceeded"));
        let post_gcd = gcd(numerator, g);
        Rat::checked(
            Some(numerator / post_gcd),
            (self.den / post_gcd).checked_mul(dd),
            "sub",
        )
    }
}

impl core::ops::Mul for Rat {
    type Output = Rat;
    fn mul(self, o: Rat) -> Rat {
        // Cross-reduce before multiplying.
        let g1 = gcd(self.num, o.den);
        let g2 = gcd(o.num, self.den);
        Rat::checked(
            (self.num / g1).checked_mul(o.num / g2),
            (self.den / g2).checked_mul(o.den / g1),
            "mul",
        )
    }
}

impl core::ops::Div for Rat {
    type Output = Rat;
    fn div(self, o: Rat) -> Rat {
        assert!(o.num != 0, "Rat division by zero");
        // (a/b) / (c/d) = (a*d)/(b*c). Cross-reduce magnitudes directly
        // rather than constructing d/c first: c may be i128::MIN even when
        // the final reduced quotient is perfectly representable.
        let g1 = gcd_unsigned(self.num.unsigned_abs(), o.num.unsigned_abs());
        let g2 = gcd_unsigned(o.den as u128, self.den as u128);
        let numerator_magnitude = (self.num.unsigned_abs() / g1)
            .checked_mul((o.den as u128) / g2)
            .unwrap_or_else(|| panic!("Rat div: i128 overflow — exactness domain exceeded"));
        let denominator_magnitude = ((self.den as u128) / g2)
            .checked_mul(o.num.unsigned_abs() / g1)
            .unwrap_or_else(|| panic!("Rat div: i128 overflow — exactness domain exceeded"));
        Rat {
            num: signed_magnitude(numerator_magnitude, (self.num < 0) ^ (o.num < 0), "div"),
            den: i128::try_from(denominator_magnitude)
                .unwrap_or_else(|_| panic!("Rat div: i128 overflow — exactness domain exceeded")),
        }
    }
}

impl core::ops::Neg for Rat {
    type Output = Rat;
    fn neg(self) -> Rat {
        Rat {
            num: self
                .num
                .checked_neg()
                .unwrap_or_else(|| panic!("Rat neg: i128 overflow — exactness domain exceeded")),
            den: self.den,
        }
    }
}

impl PartialOrd for Rat {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Rat {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self.num.cmp(&0), other.num.cmp(&0)) {
            (Ordering::Less, Ordering::Greater | Ordering::Equal)
            | (Ordering::Equal, Ordering::Greater) => Ordering::Less,
            (Ordering::Greater, Ordering::Less | Ordering::Equal)
            | (Ordering::Equal, Ordering::Less) => Ordering::Greater,
            (Ordering::Equal, Ordering::Equal) => Ordering::Equal,
            (Ordering::Greater, Ordering::Greater) => cmp_unsigned_fraction(
                self.num.unsigned_abs(),
                self.den as u128,
                other.num.unsigned_abs(),
                other.den as u128,
            ),
            (Ordering::Less, Ordering::Less) => cmp_unsigned_fraction(
                self.num.unsigned_abs(),
                self.den as u128,
                other.num.unsigned_abs(),
                other.den as u128,
            )
            .reverse(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn arithmetic_is_exact_and_reduced() {
        let a = Rat::new(1, 3);
        let b = Rat::new(1, 6);
        assert_eq!(a + b, Rat::new(1, 2));
        assert_eq!(a - b, Rat::new(1, 6));
        assert_eq!(a * b, Rat::new(1, 18));
        assert_eq!(a / b, Rat::int(2));
        assert_eq!(Rat::new(-2, -4), Rat::new(1, 2));
        assert_eq!(Rat::new(2, -4), Rat::new(-1, 2));
        assert!(Rat::new(1, 3) > Rat::new(1, 4));
        assert_eq!(Rat::new(i128::MIN, 1).numerator(), i128::MIN);
        assert_eq!(Rat::new(i128::MIN, i128::MIN), Rat::int(1));
        assert_eq!(
            Rat::new(i128::MIN, 1) / Rat::new(i128::MIN, 1),
            Rat::int(1),
            "division must cross-reduce before attempting an unrepresentable reciprocal"
        );
        assert_eq!(
            Rat::int(2) / Rat::new(i128::MIN, 1),
            Rat::new(-1, 1i128 << 126)
        );
        assert_eq!(
            Rat::new(i128::MIN, 1) / Rat::int(2),
            Rat::new(-(1i128 << 126), 1)
        );
        assert_eq!(Rat::new(0, i128::MIN), Rat::int(0));
        let shared = 56_713_727_820_156_410_577_229_101_238_628_035_241i128;
        let reducible_sum = Rat::new(1, 2 * shared) + Rat::new((shared - 3) / 2, 3 * shared);
        assert_eq!(
            reducible_sum,
            Rat::new(1, 6),
            "post-addition cancellation must happen before the final denominator product"
        );
        let reducible_difference = Rat::new(1, 2 * shared) - Rat::new((3 - shared) / 2, 3 * shared);
        assert_eq!(reducible_difference, Rat::new(1, 6));
        assert!(Rat::new(i128::MIN, 1) < Rat::new(i128::MIN + 1, 1));
        assert!(Rat::new(i128::MIN, 3) < Rat::new(i128::MIN + 1, 3));
        let named_overflow = std::panic::catch_unwind(|| -Rat::new(i128::MIN, 1));
        assert!(named_overflow.is_err(), "unrepresentable negation refuses");
    }
}
