//! First-law residual on the independently represented heat and work accounts.

use fs_math::eft::two_sum;

/// Keep the error in Q - W until the represented stored increment cancels its
/// large part. These original EFT primitives need no allocation or ordering
/// precondition. Nonfinite/overflowing accounts remain nonfinite for refusal.
pub(super) fn residual(stored_j: f64, heat_j: f64, work_by_j: f64) -> f64 {
    let (target, target_tail) = two_sum(heat_j, -work_by_j);
    let (difference, difference_tail) = two_sum(stored_j, -target);
    difference + (difference_tail - target_tail)
}

#[cfg(test)]
mod tests {
    use super::residual;

    #[test]
    fn g0_small_heat_cannot_disappear_inside_the_mechanical_target() {
        let stored = 1e16;
        let heat = 1.1;
        let work = -9_999_999_999_999_998.0;
        assert_eq!(stored - (heat - work), 0.0);
        // The large represented mechanical terms cancel to exactly 2 J.
        // Their independently authored heat leaves this visible mismatch.
        assert_eq!(residual(stored, heat, work), 2.0 - heat);
        assert_eq!(residual(-stored, -heat, -work), heat - 2.0);
    }

    #[test]
    fn g0_closed_heat_work_and_insulated_mechanics_remain_exact() {
        assert_eq!(residual(1e16, 2.0, -9_999_999_999_999_998.0), 0.0);
        assert_eq!(residual(1e16, 0.0, -1e16), 0.0);
        assert_eq!(residual(0.0, 0.0, 0.0), 0.0);
    }

    #[test]
    fn g4_overflow_and_invalid_accounts_cannot_report_a_finite_balance() {
        assert!(!residual(0.0, f64::MAX, -f64::MAX).is_finite());
        assert!(!residual(f64::MAX, -f64::MAX, 0.0).is_finite());
        assert!(!residual(f64::NAN, 0.0, 0.0).is_finite());
    }
}
