//! Full native velocity-correction norms and explicit convergence tolerances.
//!
//! No pressure datum, implicit physical velocity floor, or relaxed step length
//! changes this decision. The caller separately requires BOTH momentum-system
//! residual criteria and admits its tolerances before iteration.

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(super) struct VelocityUpdate {
    /// Euclidean norm of the original BDM velocity-moment correction [m/s].
    pub norm_m_s: f64,
    /// Relative correction using the fixed initial/first-candidate reference.
    pub relative: f64,
}

impl VelocityUpdate {
    /// Measure the FULL Newton correction before line-search relaxation.
    pub fn measure(
        correction: impl IntoIterator<Item = f64>,
        reference_norm_m_s: f64,
    ) -> Result<Self, &'static str> {
        if !reference_norm_m_s.is_finite() || reference_norm_m_s < 0.0 {
            return Err("Forchheimer velocity reference is not representable");
        }
        let mut norm_m_s = 0.0_f64;
        for value in correction {
            if !value.is_finite() {
                return Err("Forchheimer full velocity correction is not representable");
            }
            norm_m_s = norm_m_s.hypot(value);
        }
        if !norm_m_s.is_finite() {
            return Err("Forchheimer full velocity-correction norm is not representable");
        }
        let relative = if norm_m_s == 0.0 {
            0.0
        } else {
            norm_m_s / reference_norm_m_s
        };
        if !relative.is_finite() {
            return Err("Forchheimer relative velocity correction is not representable");
        }
        Ok(Self { norm_m_s, relative })
    }

    /// Both tolerances have already been admitted by the nonlinear owner.
    /// An absolute tolerance of zero retains strict relative-only behavior,
    /// apart from the already-converged exactly-zero correction.
    pub fn converged(self, relative_tolerance: f64, absolute_tolerance_m_s: f64) -> bool {
        self.relative <= relative_tolerance || self.norm_m_s <= absolute_tolerance_m_s
    }
}

#[cfg(test)]
mod tests {
    use super::VelocityUpdate;

    #[test]
    fn exact_zero_needs_no_reference_or_absolute_allowance() {
        let update = VelocityUpdate::measure([0.0, -0.0], 0.0).unwrap();
        assert_eq!(update, VelocityUpdate::default());
        assert!(update.converged(1e-9, 0.0));
    }

    #[test]
    fn strict_default_does_not_ignore_a_small_but_relatively_large_correction() {
        let update = VelocityUpdate::measure([3e-16, 4e-16], 5e-16).unwrap();
        assert!((update.relative - 1.0).abs() < 2e-16);
        assert!(!update.converged(1e-9, 0.0));
        assert!(update.converged(1e-9, 1e-14));
    }

    #[test]
    fn absolute_threshold_uses_full_euclidean_moment_norm() {
        let update = VelocityUpdate::measure([3.0, -4.0], 5.0).unwrap();
        assert_eq!(update.norm_m_s, 5.0);
        assert_eq!(update.relative, 1.0);
        assert!(!update.converged(1e-9, 4.9));
        assert!(update.converged(1e-9, 5.0));
        assert!(update.converged(1.0, 0.0));
    }

    #[test]
    fn explicit_small_absolute_allowance_preserves_resolvable_weak_flow() {
        let update = VelocityUpdate::measure([1e-6, 0.0], 1e-6).unwrap();
        assert!(!update.converged(1e-9, 1e-12));
        // A large pressure datum never enters this physical update criterion.
        let later = VelocityUpdate::measure([3.3e-7, 0.0], 1e-6).unwrap();
        assert!(!later.converged(1e-9, 1e-12));
    }

    #[test]
    fn scale_changes_require_corresponding_explicit_absolute_units() {
        let small = VelocityUpdate::measure([3e-100, 4e-100], 5e-100).unwrap();
        let large = VelocityUpdate::measure([3e100, 4e100], 5e100).unwrap();
        assert!((small.relative - large.relative).abs() < 3e-16);
        assert!(small.converged(1e-9, 6e-100));
        assert!(!large.converged(1e-9, 6e-100));
        assert!(large.converged(1e-9, 6e100));
    }

    #[test]
    fn norms_avoid_naive_square_overflow_and_underflow() {
        for scale in [1e-200, 1e200] {
            let update = VelocityUpdate::measure([3.0 * scale, 4.0 * scale], 5.0 * scale).unwrap();
            assert!((update.norm_m_s / (5.0 * scale) - 1.0).abs() < 3e-16);
            assert!((update.relative - 1.0).abs() < 3e-16);
        }
    }

    #[test]
    fn malformed_or_unrepresentable_updates_refuse() {
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(VelocityUpdate::measure([value], 1.0).is_err());
            assert!(VelocityUpdate::measure([0.0], value).is_err());
        }
        assert!(VelocityUpdate::measure([0.0], -1.0).is_err());
        assert!(VelocityUpdate::measure([1.0], 0.0).is_err());
        assert!(VelocityUpdate::measure([f64::MAX, f64::MAX], f64::MAX).is_err());
        assert!(VelocityUpdate::measure([1.0], f64::from_bits(1)).is_err());
    }
}
