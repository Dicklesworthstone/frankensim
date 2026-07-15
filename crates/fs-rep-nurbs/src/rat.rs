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

/// Minimal unsigned 256-bit scratch used to keep fallible rational midpoint
/// arithmetic exact without adding a runtime dependency to the crate.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct Wide {
    hi: u128,
    lo: u128,
}

impl Wide {
    const ZERO: Wide = Wide { hi: 0, lo: 0 };

    fn product(left: u128, right: u128) -> Wide {
        const MASK: u128 = u64::MAX as u128;
        let (left_lo, left_hi) = (left & MASK, left >> 64);
        let (right_lo, right_hi) = (right & MASK, right >> 64);
        let low_product = left_lo * right_lo;
        let cross_left = left_hi * right_lo;
        let cross_right = left_lo * right_hi;
        let high_product = left_hi * right_hi;
        let middle = (low_product >> 64) + (cross_left & MASK) + (cross_right & MASK);
        let lo = (low_product & MASK) | ((middle & MASK) << 64);
        let hi = high_product
            .checked_add(cross_left >> 64)
            .and_then(|value| value.checked_add(cross_right >> 64))
            .and_then(|value| value.checked_add(middle >> 64))
            .expect("u128 product decomposition must fit its high limb");
        Wide { hi, lo }
    }

    fn checked_add(self, other: Wide) -> Option<Wide> {
        let (lo, carry) = self.lo.overflowing_add(other.lo);
        let hi = self
            .hi
            .checked_add(other.hi)?
            .checked_add(u128::from(carry))?;
        Some(Wide { hi, lo })
    }

    fn subtract(self, other: Wide) -> Wide {
        debug_assert!(self >= other);
        let (lo, borrow) = self.lo.overflowing_sub(other.lo);
        let hi = self.hi - other.hi - u128::from(borrow);
        Wide { hi, lo }
    }

    fn bit(self, index: usize) -> u128 {
        if index >= 128 {
            (self.hi >> (index - 128)) & 1
        } else {
            (self.lo >> index) & 1
        }
    }

    fn div_rem_u128(self, divisor: u128) -> (Wide, u128) {
        debug_assert_ne!(divisor, 0);
        let mut quotient = Wide::ZERO;
        let mut remainder = 0u128;
        for index in (0..256).rev() {
            quotient.hi = (quotient.hi << 1) | (quotient.lo >> 127);
            quotient.lo <<= 1;

            let carry = remainder >> 127 != 0;
            let doubled = (remainder << 1) | self.bit(index);
            if carry || doubled >= divisor {
                // When `carry` is set, the true 129-bit value is
                // 2^128+doubled. Since the old remainder was below divisor,
                // subtracting once yields a value below divisor; wrapping
                // subtraction computes that low limb exactly.
                remainder = doubled.wrapping_sub(divisor);
                quotient.lo |= 1;
            } else {
                remainder = doubled;
            }
        }
        (quotient, remainder)
    }

    fn remainder_u128(self, divisor: u128) -> u128 {
        self.div_rem_u128(divisor).1
    }
}

fn signed_product_sum(
    left_num: i128,
    left_scale: u128,
    right_num: i128,
    right_scale: u128,
) -> Option<(bool, Wide)> {
    let left = Wide::product(left_num.unsigned_abs(), left_scale);
    let right = Wide::product(right_num.unsigned_abs(), right_scale);
    let left_negative = left_num < 0;
    let right_negative = right_num < 0;
    if left_negative == right_negative {
        let magnitude = left.checked_add(right)?;
        return Some((left_negative && magnitude != Wide::ZERO, magnitude));
    }
    Some(match left.cmp(&right) {
        Ordering::Greater => (left_negative, left.subtract(right)),
        Ordering::Less => (right_negative, right.subtract(left)),
        Ordering::Equal => (false, Wide::ZERO),
    })
}

fn checked_signed_wide(magnitude: Wide, negative: bool) -> Option<i128> {
    if magnitude.hi != 0 {
        return None;
    }
    if negative && magnitude.lo == (1u128 << 127) {
        return Some(i128::MIN);
    }
    let value = i128::try_from(magnitude.lo).ok()?;
    if negative {
        value.checked_neg()
    } else {
        Some(value)
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

    /// Exact midpoint when the fully reduced result stays inside the
    /// i128-backed representation domain. Wide scratch prevents provisional
    /// products or sums from falsely refusing a representable result.
    /// `None` is a fail-closed exactness refusal; no rounded midpoint is
    /// returned.
    #[must_use]
    pub(crate) fn checked_midpoint(self, other: Rat) -> Option<Rat> {
        let left_den = self.den as u128;
        let right_den = other.den as u128;
        let common = gcd_unsigned(left_den, right_den);
        let left_scale = right_den / common;
        let right_scale = left_den / common;
        let (negative, numerator) =
            signed_product_sum(self.num, left_scale, other.num, right_scale)?;
        if numerator == Wide::ZERO {
            return Some(Rat::int(0));
        }

        // For reduced inputs, N = a*(d/g) + c*(b/g) is coprime to the two
        // co-prime denominator wings b/g and d/g. Therefore only 2*g can
        // cancel from the midpoint denominator. Compute that gcd through a
        // wide remainder, divide before narrowing, then form the already
        // reduced denominator with checked u128 products.
        let doubled_common = common.checked_mul(2)?;
        let reduction = gcd_unsigned(numerator.remainder_u128(doubled_common), doubled_common);
        let (reduced_numerator, remainder) = numerator.div_rem_u128(reduction);
        debug_assert_eq!(remainder, 0);
        let reduced_denominator = (left_den / common)
            .checked_mul(right_den / common)?
            .checked_mul(doubled_common / reduction)?;
        let denominator = i128::try_from(reduced_denominator).ok()?;
        let numerator = checked_signed_wide(reduced_numerator, negative)?;
        Some(Rat::new(numerator, denominator))
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
    fn wide_product_and_division_cover_both_high_limbs() {
        let product = Wide::product(u128::MAX, u128::MAX);
        assert_eq!(
            product,
            Wide {
                hi: u128::MAX - 1,
                lo: 1,
            }
        );
        let (quotient, remainder) = product.div_rem_u128(u128::MAX);
        assert_eq!(
            quotient,
            Wide {
                hi: 0,
                lo: u128::MAX
            }
        );
        assert_eq!(remainder, 0);
    }

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

    #[test]
    fn checked_midpoint_avoids_representable_extreme_sum_overflow() {
        assert_eq!(
            Rat::new(i128::MAX, 1).checked_midpoint(Rat::new(i128::MAX - 2, 1)),
            Some(Rat::new(i128::MAX - 1, 1))
        );
        assert_eq!(
            Rat::new(i128::MIN, 1).checked_midpoint(Rat::new(i128::MIN + 2, 1)),
            Some(Rat::new(i128::MIN + 1, 1))
        );
        assert_eq!(
            Rat::new(i128::MIN, 1).checked_midpoint(Rat::new(i128::MAX, 1)),
            Some(Rat::new(-1, 2))
        );
        assert_eq!(
            Rat::new(1, 3).checked_midpoint(Rat::new(2, 5)),
            Some(Rat::new(11, 30))
        );
        let max_third = i128::MAX / 3;
        assert_eq!(
            Rat::new(i128::MAX, 3).checked_midpoint(Rat::new(i128::MAX - 5, 3)),
            Some(Rat::new(max_third * 2 - 1, 2)),
            "same-denominator cancellation must precede narrowing"
        );
        let denominator = i128::MAX - 4;
        assert_eq!(
            Rat::new(1, denominator).checked_midpoint(Rat::new(2, denominator)),
            Some(Rat::new(1, 2 * (denominator / 3))),
            "denominator cancellation must precede doubling"
        );
        assert_eq!(
            Rat::new(i128::MAX, 2).checked_midpoint(Rat::new(-i128::MAX, 3)),
            Some(Rat::new(i128::MAX, 12)),
            "opposite-signed wide products must cancel before narrowing"
        );
        assert_eq!(
            Rat::new(i128::MIN, 1).checked_midpoint(Rat::new(i128::MIN, 1)),
            Some(Rat::new(i128::MIN, 1))
        );
        assert_eq!(
            Rat::new(i128::MAX, 1).checked_midpoint(Rat::new(i128::MAX - 1, 1)),
            None,
            "an unrepresentable positive reduced numerator must refuse"
        );
        assert_eq!(
            Rat::new(1, i128::MAX).checked_midpoint(Rat::int(0)),
            None,
            "an unrepresentable reduced denominator must refuse"
        );
        assert_eq!(
            Rat::new(7, 11).checked_midpoint(Rat::new(-7, 11)),
            Some(Rat::int(0))
        );
    }
}
