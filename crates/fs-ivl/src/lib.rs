//! fs-ivl — certified arithmetic (plan §6.4): outward-rounded intervals and
//! affine forms. This crate is WHAT THE WORD "CERTIFIED" MEANS everywhere
//! else in the system: every operation's postcondition is ENCLOSURE (the
//! result contains the true value set), tested as the G0 containment law.
//!
//! - [`Interval`]: directed rounding via fs-math's `next_up`/`next_down`
//!   nudging — no global rounding-mode state anywhere (grep-lintable:
//!   this workspace never touches the FPU control word). Elementary
//!   functions inherit fs-math's DECLARED ULP budgets.
//! - [`decay::exp_neg`]: nonnegative exponential attenuation with a proved
//!   series remainder, without an elementary-function ULP assumption.
//! - [`Affine`]/[`AffineCtx`]: noise-symbol forms that kill the dependency
//!   problem on correlated expressions (x − x, deep F-rep DAGs).
//! - The high-precision oracle rungs live in `fs_math::{eft, dd}` (L0;
//!   single implementation shared with fs-la's iterative refinement —
//!   recorded relocation, beads 6ys.8/6ys.12).
//! - [`TaylorModel1`]: bounded fixed-order univariate models with a rigorous
//!   interval remainder and fallible finite-input admission. Multivariate
//!   Taylor models and quad-double arithmetic remain follow-up scope.
//!
//! Determinism: everything here is straight-line IEEE arithmetic on
//! fs-math strict functions — cross-ISA bit-deterministic BY CONSTRUCTION
//! (golden-hashed in tests/conformance.rs, verified on both reference
//! ISAs).

pub mod affine;
pub mod angle;
pub mod expansion;
pub mod formal_binding;
pub mod formal_manifest;
pub mod formal_proofs;
pub mod interval;
pub mod multi_taylor;
pub mod newton;
pub mod predicates;
pub mod taylor;

/// Exponential attenuation without a libm error assumption.
///
/// Range reduction gives u = x / 2^k <= 1/16. After degree 12 the omitted
/// positive series for exp(u) is bounded by t_13 / (1 - u/14), since every
/// subsequent term ratio is at most u/14. Outward arithmetic encloses that
/// sum; reciprocation and k squarings enclose exp(-x).
pub mod decay {
    use crate::Interval;

    /// Version of the series, range-reduction and underflow policy.
    pub const DECAY_ENCLOSURE_SEMANTICS_VERSION: u32 = 1;

    /// Enclose exp(-optical_depth) for a nonnegative interval.
    ///
    /// Negative lower endpoints refuse rather than changing the domain.
    /// Positive infinity has limiting value zero. Finite positive underflow
    /// retains a positive upper bound. Work is at most two degree-12 series
    /// and thirty squarings, with no allocation or global rounding state.
    /// The ordinary fs-ivl IEEE arithmetic/subnormal requirements apply.
    #[must_use]
    pub fn exp_neg(optical_depth: Interval) -> Option<Interval> {
        if optical_depth.lo() < 0.0 {
            return None;
        }
        Some(Interval::new(
            endpoint(optical_depth.hi()).lo(),
            endpoint(optical_depth.lo()).hi(),
        ))
    }

    fn endpoint(x: f64) -> Interval {
        if x == 0.0 {
            return Interval::point(1.0);
        }
        if x == f64::INFINITY {
            return Interval::point(0.0);
        }
        if x >= 2048.0 {
            // e > 2, hence exp(-2048) < 2^-2048 < the least binary64.
            return Interval::new(0.0, f64::from_bits(1));
        }
        let mut u = x;
        let mut squares = 0;
        while u > 0.0625 {
            // Exact binary scaling in the normal range; at most 15 steps.
            u *= 0.5;
            squares += 1;
        }
        let one = Interval::point(1.0);
        let u = Interval::point(u);
        let mut term = one;
        let mut sum = one;
        for n in 1..=12 {
            term = nonnegative(term * u / Interval::point(f64::from(n)));
            sum = sum + term;
        }
        let next = nonnegative(term * u / Interval::point(13.0));
        let tail = next / (one - u / Interval::point(14.0));
        let exponential = Interval::new(sum.lo(), (sum + tail).hi());
        let mut result = unit(one / exponential);
        for _ in 0..squares {
            result = unit(result * result);
        }
        result
    }

    fn nonnegative(value: Interval) -> Interval {
        Interval::new(value.lo().max(0.0), value.hi())
    }

    fn unit(value: Interval) -> Interval {
        Interval::new(value.lo().max(0.0), value.hi().min(1.0))
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        #[test]
        fn known_values_and_range_reduction_boundaries() {
            // Independently rounded high-precision values, not fs-math exp.
            for (x, reference) in [
                (0.0625, 0.939_413_062_813_475_8_f64),
                (0.5, 0.606_530_659_712_633_4),
                (1.0, 0.367_879_441_171_442_33),
                (2.0, 0.135_335_283_236_612_7),
                (10.0, 0.000_045_399_929_762_484_854),
                (100.0, 3.720_075_976_020_836e-44),
                (700.0, 9.859_676_543_759_77e-305),
            ] {
                let bound = exp_neg(Interval::point(x)).unwrap();
                assert!(bound.lo() <= reference.next_down(), "x={x}: {bound:?}");
                assert!(bound.hi() >= reference.next_up(), "x={x}: {bound:?}");
                assert!(bound.width() <= reference * 1e-9);
            }
        }

        #[test]
        fn zero_infinity_and_invalid_domain_are_distinct() {
            assert_eq!(exp_neg(Interval::point(-0.0)), Some(Interval::point(1.0)));
            assert_eq!(
                exp_neg(Interval::point(f64::INFINITY)),
                Some(Interval::point(0.0))
            );
            assert_eq!(
                exp_neg(Interval::new(0.0, f64::INFINITY)),
                Some(Interval::new(0.0, 1.0))
            );
            assert!(exp_neg(Interval::new(-f64::from_bits(1), 1.0)).is_none());
            assert!(exp_neg(Interval::WHOLE).is_none());
        }

        #[test]
        fn finite_positive_underflow_does_not_claim_exact_zero() {
            for x in [745.0, 1000.0, 2048.0, f64::MAX] {
                let bound = exp_neg(Interval::point(x)).unwrap();
                assert!(bound.lo() >= 0.0);
                assert!(bound.hi() > 0.0 && bound.hi() < f64::MIN_POSITIVE);
            }
            let tiny = exp_neg(Interval::point(f64::from_bits(1))).unwrap();
            assert!(tiny.lo() < 1.0 && tiny.hi() == 1.0);
        }

        #[test]
        fn interval_endpoints_and_split_path_contain_a_common_truth() {
            for (a, b) in [(0.1, 0.2), (1.0, 2.0), (100.0, 300.0)] {
                let bound = exp_neg(Interval::new(a, b)).unwrap();
                assert_eq!(bound.lo(), exp_neg(Interval::point(b)).unwrap().lo());
                assert_eq!(bound.hi(), exp_neg(Interval::point(a)).unwrap().hi());
                let whole = exp_neg(Interval::point(a) + Interval::point(b)).unwrap();
                let parts = exp_neg(Interval::point(a)).unwrap()
                    * exp_neg(Interval::point(b)).unwrap();
                assert!(whole.intersect(parts).is_some());
            }
        }
    }
}

pub use affine::{Affine, AffineCtx};
pub use interval::Interval;
pub use multi_taylor::{
    MAX_MODEL_MEMORY_BYTES, MAX_MULTIVARIATE_DIM, MAX_MULTIVARIATE_ORDER, MAX_MULTIVARIATE_TERMS,
    TaylorModel, VariableInfo, binomial, generate_multi_indices, term_count,
};
pub use newton::{
    DEFAULT_MAX_ROOT_BOXES, RootBox, RootSearchConfig, RootSearchError, RootSearchOutcome,
    RootSearchReport, RootSearchTermination, krawczyk_step, lipschitz_bound, newton_roots,
    newton_roots_bounded, newton_roots_with_checkpoint,
};
pub use predicates::{
    Sign, Stage, incircle, incircle_with_stage, insphere, insphere_with_stage, orient2d,
    orient2d_sos, orient2d_with_stage, orient3d, orient3d_sos, orient3d_with_stage,
};
pub use taylor::{MAX_TAYLOR_ORDER, TaylorModel1, TaylorModelError};

/// Crate version, re-exported for provenance stamping (the Five Explicits'
/// "versions" pillar reaches down to individual crates).
pub const VERSION: &str = env!("CARGO_PKG_VERSION");
/// Semantic version of the outward interval operation and refusal policy.
///
/// This version covers basic-operation finite-overflow handling, elementary
/// endpoint nudging from fs-math's declared budgets, and extended-real
/// indeterminate-operation fallbacks. Certificate identities must also bind
/// fs-math's strict-core semantic version and golden fingerprint.
pub const INTERVAL_SEMANTICS_VERSION: u32 = 1;

#[cfg(test)]
mod tests {
    #[test]
    fn version_is_stamped() {
        assert!(!super::VERSION.is_empty());
        assert_eq!(super::INTERVAL_SEMANTICS_VERSION, 1);
    }
}
