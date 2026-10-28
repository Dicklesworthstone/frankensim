//! Shared source-path admission for kinetic reactors and compressible flow.

use super::{KineticsErrorV1, poll, positive};
use crate::Nasa9StandardStateModelV1;
use fs_qty::Temperature;

/// Admit source support, explicit caloric continuity and strictly positive Cv
/// over the entire original NASA-9 temperature interval.
///
/// Outward-rounded original monomial bounds are subdivided at most12levels and
///4095pieces per source region. An unresolved bound refuses; sampling positive
/// values is never treated as proof that an intervening polynomial stays positive.
/// No source coefficients or formation/entropy references are altered.
pub fn admit_caloric_interval(
    model: &Nasa9StandardStateModelV1,
    lower: Temperature,
    upper: Temperature,
    mut cancelled: impl FnMut() -> bool,
) -> Result<(), KineticsErrorV1> {
    poll(&mut cancelled)?;
    let lo = lower.value();
    let hi = upper.value();
    if !positive(lo) || !positive(hi) || lo > hi {
        return Err(KineticsErrorV1::Invalid("invalid caloric source interval"));
    }
    let first = model
        .region_for_temperature(lo)
        .ok_or(KineticsErrorV1::Invalid(
            "lower temperature lacks source support",
        ))?
        .0;
    let last = model
        .region_for_temperature(hi)
        .ok_or(KineticsErrorV1::Invalid(
            "upper temperature lacks source support",
        ))?
        .0;
    if first != last && model.caloric_reference().is_none() {
        return Err(KineticsErrorV1::Invalid(
            "cross-region path requires explicit caloric continuity",
        ));
    }
    for region in &model.regions()[first..=last] {
        positive_cv(
            *region.coefficients(),
            lo.max(region.temperature_min().value()),
            hi.min(region.temperature_max().value()),
            &mut cancelled,
        )?;
    }
    Ok(())
}

// Each original Cv/R monomial is monotone on positive T. Directed rounding
// bounds every multiply and sum, including the cp/R constant minus one.
fn cv_lower_bound(a: [f64; 9], lo: f64, hi: f64) -> f64 {
    let mut low = (a[2] - 1.0).next_down();
    for (i, (coefficient, l, h)) in [
        (a[0], (1.0 / hi).next_down(), (1.0 / lo).next_up()),
        (a[1], (1.0 / hi).next_down(), (1.0 / lo).next_up()),
        (a[3], lo, hi),
        (a[4], lo, hi),
        (a[5], lo, hi),
        (a[6], lo, hi),
    ]
    .into_iter()
    .enumerate()
    {
        if coefficient == 0.0 {
            continue;
        }
        let power = match i {
            0 => 2,
            1 | 2 => 1,
            3 => 2,
            4 => 3,
            _ => 4,
        };
        let (mut lower, mut upper) = (1.0_f64, 1.0_f64);
        for _ in 0..power {
            lower = (lower * l).next_down();
            upper = (upper * h).next_up();
        }
        let term = coefficient * if coefficient > 0.0 { lower } else { upper };
        low = (low + term.next_down()).next_down();
    }
    low
}
fn positive_cv(
    a: [f64; 9],
    lo: f64,
    hi: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), KineticsErrorV1> {
    let mut stack = [(0.0, 0.0, 0_u8); 14];
    stack[0] = (lo, hi, 0);
    let mut n = 1;
    let mut work = 0;
    while n > 0 {
        poll(cancelled)?;
        if work >= 4095 {
            return Err(KineticsErrorV1::Invalid(
                "positive Cv interval subdivision budget exhausted",
            ));
        }
        work += 1;
        n -= 1;
        let (l, h, depth) = stack[n];
        let bound = cv_lower_bound(a, l, h);
        if bound.is_finite() && bound > 0.0 {
            continue;
        }
        let mid = l + 0.5 * (h - l);
        if depth == 12 || mid <= l || mid >= h {
            return Err(KineticsErrorV1::Invalid(
                "positive component Cv could not be bounded on the source path",
            ));
        }
        stack[n] = (mid, h, depth + 1);
        stack[n + 1] = (l, mid, depth + 1);
        n += 2;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn g0_interval_cv_refuses_an_interior_negative_pocket_and_preserves_cancellation() {
        // Cv/R=(T-500)^2/10000-1: both endpoints are positive, midpoint negative.
        let mut a = [0.0; 9];
        a[2] = 25.0;
        a[3] = -0.1;
        a[4] = 0.0001;
        assert!(cv_lower_bound(a, 200.0, 800.0) < 0.0);
        assert!(positive_cv(a, 200.0, 800.0, &mut || false).is_err());
        a[2] = 27.0;
        assert!(positive_cv(a, 200.0, 800.0, &mut || false).is_ok());
        assert_eq!(
            positive_cv(a, 200.0, 800.0, &mut || true),
            Err(KineticsErrorV1::Cancelled)
        );
    }
}
