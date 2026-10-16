//! The deterministic reduction library (plan §5.4, Decalogue P2): fixed-
//! shape pairwise trees keyed by tile/element index, Neumaier-compensated
//! accumulation, deterministic tie-breaking, and an order-sensitivity audit
//! that catches arrival-order bugs with localization.
//!
//! THE SHAPE RULE (part of the contract): a fold over `n` ordered items
//! splits at the largest power of two strictly below `n` and recurses —
//! a pure function of `n`, never of scheduling, arrival order, or thread
//! count. Every rounding decision in a deterministic-mode reduction is
//! therefore replayable from the plan alone.
//!
//! Cross-ISA stance (documented, not claimed): identical shapes make
//! results a function of scalar arithmetic only; remaining divergence
//! classes (FMA contraction, libm ULP) are fs-math's contract and the G5
//! cross-ISA report's subject.

use crate::kernel::Reduce;

/// Fixed accumulation block size for `det_sum`/`det_dot`: part of the
/// reduction-shape contract (changing it changes bits).
const BLOCK: usize = 256;

/// Fold ordered items up the fixed-shape pairwise tree. The tree depends
/// only on `items.len()` (see module docs), so the result is bit-identical
/// for a given input sequence regardless of how the inputs were produced.
/// Order-sensitive merges (e.g. concatenation) see items in ascending
/// index order.
#[must_use]
pub fn pairwise_fold<T: Reduce>(mut items: Vec<T>) -> T {
    match items.len() {
        0 => T::identity(),
        1 => items.pop().expect("len checked"),
        n => {
            let split = largest_pow2_below(n);
            let right = items.split_off(split);
            pairwise_fold(items).merge(pairwise_fold(right))
        }
    }
}

/// Largest power of two strictly below `n` (n >= 2).
fn largest_pow2_below(n: usize) -> usize {
    debug_assert!(n >= 2);
    let p = usize::BITS - (n - 1).leading_zeros() - 1;
    1 << p
}

/// A Neumaier-compensated partial sum: `Reduce`-composable, so compensated
/// GLOBAL sums ride the same fixed tree as everything else. Merging two
/// partials two-sums the sums and adds the carried compensations — a
/// deterministic, shape-stable combine.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Compensated {
    /// Running sum.
    pub sum: f64,
    /// Running compensation (lost low-order bits).
    pub comp: f64,
}

impl Compensated {
    /// The zero partial.
    #[must_use]
    pub const fn zero() -> Self {
        Compensated {
            sum: 0.0,
            comp: 0.0,
        }
    }

    /// Accumulate one term (Neumaier's variant: compensation also captures
    /// the case where the term dominates the running sum).
    #[must_use]
    pub fn accumulate(self, x: f64) -> Self {
        let t = self.sum + x;
        let comp = if self.sum.abs() >= x.abs() {
            self.comp + ((self.sum - t) + x)
        } else {
            self.comp + ((x - t) + self.sum)
        };
        Compensated { sum: t, comp }
    }

    /// The compensated total.
    #[must_use]
    pub fn value(self) -> f64 {
        self.sum + self.comp
    }
}

impl Reduce for Compensated {
    fn identity() -> Self {
        Compensated::zero()
    }

    fn merge(self, other: Self) -> Self {
        // Two-sum of the partial sums; compensations add exactly enough
        // that `value()` of the merge equals the compensated total.
        let s = self.sum + other.sum;
        let bb = s - self.sum;
        let err = (self.sum - (s - bb)) + (other.sum - bb);
        Compensated {
            sum: s,
            comp: self.comp + other.comp + err,
        }
    }
}

/// Deterministic compensated sum of a slice: per-block Neumaier
/// accumulation (block size fixed) merged up the pairwise tree. Shape is a
/// pure function of `xs.len()`.
#[must_use]
pub fn det_sum(xs: &[f64]) -> f64 {
    let partials: Vec<Compensated> = xs
        .chunks(BLOCK)
        .map(|c| c.iter().fold(Compensated::zero(), |a, &x| a.accumulate(x)))
        .collect();
    pairwise_fold(partials).value()
}

/// Deterministic compensated dot product (same blocking/tree as
/// [`det_sum`]; products are formed unfused — `a[i] * b[i]` — so the result
/// is identical whether or not the target fuses multiply-add).
///
/// # Panics
/// When lengths differ (programmer-error contract, checked first).
#[must_use]
pub fn det_dot(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len(), "dot operands must match");
    let partials: Vec<Compensated> = a
        .chunks(BLOCK)
        .zip(b.chunks(BLOCK))
        .map(|(ca, cb)| {
            ca.iter()
                .zip(cb)
                .fold(Compensated::zero(), |acc, (&x, &y)| acc.accumulate(x * y))
        })
        .collect();
    pairwise_fold(partials).value()
}

/// Deterministic Euclidean norm (sqrt of [`det_dot`] with itself).
#[must_use]
pub fn det_norm2(xs: &[f64]) -> f64 {
    det_dot(xs, xs).sqrt()
}

/// Deterministic minimum under IEEE total order (NaN sorts after +inf, so
/// data NaNs cannot poison comparisons nondeterministically).
#[must_use]
pub fn det_min(xs: &[f64]) -> Option<f64> {
    xs.iter().copied().min_by(f64::total_cmp)
}

/// Deterministic maximum under IEEE total order.
#[must_use]
pub fn det_max(xs: &[f64]) -> Option<f64> {
    xs.iter().copied().max_by(f64::total_cmp)
}

/// Deterministic argmin: ties break to the LOWEST index (the tie-breaking
/// law of plan §5.4 — "argmin ties → lowest logical index").
#[must_use]
pub fn det_argmin(xs: &[f64]) -> Option<usize> {
    let mut best: Option<(usize, f64)> = None;
    for (i, &x) in xs.iter().enumerate() {
        match best {
            None => best = Some((i, x)),
            Some((_, b)) if x.total_cmp(&b).is_lt() => best = Some((i, x)),
            _ => {}
        }
    }
    best.map(|(i, _)| i)
}

/// Deterministic argmax: ties break to the LOWEST index.
#[must_use]
pub fn det_argmax(xs: &[f64]) -> Option<usize> {
    let mut best: Option<(usize, f64)> = None;
    for (i, &x) in xs.iter().enumerate() {
        match best {
            None => best = Some((i, x)),
            Some((_, b)) if x.total_cmp(&b).is_gt() => best = Some((i, x)),
            _ => {}
        }
    }
    best.map(|(i, _)| i)
}

/// Deterministic inclusive prefix sum (compensated, sequential — scans are
/// prefix-dependent by definition, so the deterministic shape IS the
/// sequential order).
#[must_use]
pub fn det_prefix_sum(xs: &[f64]) -> Vec<f64> {
    let mut acc = Compensated::zero();
    xs.iter()
        .map(|&x| {
            acc = acc.accumulate(x);
            acc.value()
        })
        .collect()
}

/// An arrival-order bug report from [`audit_accumulator`]: feeding the
/// items in a permuted completion order changed the result, and `witness`
/// is the number of leading (index-ordered) items that already suffice to
/// expose the divergence — the localization handle for debugging.
#[derive(Debug, Clone, PartialEq)]
pub struct OrderSensitivity {
    /// Result under ascending index order.
    pub ordered: f64,
    /// Result under the seeded permuted order.
    pub permuted: f64,
    /// Smallest prefix length that already diverges.
    pub witness: usize,
}

impl core::fmt::Display for OrderSensitivity {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(
            f,
            "accumulator is arrival-order sensitive: ordered={} permuted={} (first {} items \
             suffice to expose it); key the reduction by logical index instead (P2)",
            self.ordered, self.permuted, self.witness
        )
    }
}

/// The G5 audit utility: drives an accumulator (`fold(acc, item)`) with
/// items in ascending order and in a seeded permuted "completion order".
/// A divergence is an arrival-order bug — reported with the smallest
/// exposing prefix (acceptance: "seeded arrival-order bug caught by the
/// audit with correct localization").
///
/// # Errors
/// [`OrderSensitivity`] when the accumulator's result depends on order.
pub fn audit_accumulator(
    items: &[f64],
    seed: u64,
    fold: impl Fn(f64, f64) -> f64,
) -> Result<(), OrderSensitivity> {
    let run = |subset: &[f64], permute: bool| -> f64 {
        let mut order: Vec<usize> = (0..subset.len()).collect();
        if permute {
            // Deterministic Fisher–Yates from the seed (high LCG bits).
            let mut state = seed | 1;
            for i in (1..order.len()).rev() {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                let j = ((state >> 32) as usize) % (i + 1);
                order.swap(i, j);
            }
        }
        order.into_iter().map(|i| subset[i]).fold(0.0, &fold)
    };
    let ordered = run(items, false);
    let permuted = run(items, true);
    if ordered.to_bits() == permuted.to_bits() {
        return Ok(());
    }
    // Localize: smallest prefix whose ordered/permuted results diverge.
    let mut witness = items.len();
    for k in 2..=items.len() {
        if run(&items[..k], false).to_bits() != run(&items[..k], true).to_bits() {
            witness = k;
            break;
        }
    }
    Err(OrderSensitivity {
        ordered,
        permuted,
        witness,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Test-only syntax tree whose merge records the exact reduction shape.
    /// Empty is a lawful two-sided identity, as required by `Reduce`.
    #[derive(Debug, PartialEq, Eq)]
    enum Shape {
        Empty,
        Leaf(u64),
        Node(Box<Shape>, Box<Shape>),
    }

    impl Reduce for Shape {
        fn identity() -> Self {
            Shape::Empty
        }

        fn merge(self, other: Self) -> Self {
            match (self, other) {
                (Shape::Empty, right) => right,
                (left, Shape::Empty) => left,
                (left, right) => Shape::Node(Box::new(left), Box::new(right)),
            }
        }
    }

    /// Independent statement of the contract: for n >= 2,
    /// next_power_of_two(n) / 2 is the largest power of two below n.
    fn expected_shape(start: u64, n: usize) -> Shape {
        match n {
            0 => Shape::Empty,
            1 => Shape::Leaf(start),
            _ => {
                let split = n.next_power_of_two() / 2;
                Shape::Node(
                    Box::new(expected_shape(start, split)),
                    Box::new(expected_shape(
                        start + u64::try_from(split).expect("bounded split fits u64"),
                        n - split,
                    )),
                )
            }
        }
    }

    #[test]
    fn generated_tree_shapes_match_the_independent_length_oracle() {
        fs_propcheck::check(
            "pairwise-fold-length-keyed-shape",
            0xE008_0001,
            512,
            |s| {
                if s.next_u64() & 1 == 0 {
                    let exponent =
                        u32::try_from(s.int_in(0, 12)).expect("exponent is non-negative");
                    let power = 1_i64 << exponent;
                    u64::try_from((power + s.int_in(-2, 2)).clamp(0, 4_096))
                        .expect("clamped length fits u64")
                } else {
                    u64::try_from(s.int_in(0, 4_096)).expect("length is non-negative")
                }
            },
            |n| {
                let n = usize::try_from(*n).expect("generated length fits usize");
                let leaves: Vec<Shape> = (0..n)
                    .map(|index| Shape::Leaf(u64::try_from(index).expect("bounded index fits u64")))
                    .collect();
                pairwise_fold(leaves) == expected_shape(0, n)
            },
        );
    }

    /// Double-double oracle for compensated-sum accuracy (test-only; the
    /// production dd ladder is fs-ivl's).
    fn dd_sum(xs: &[f64]) -> f64 {
        let (mut hi, mut lo) = (0.0f64, 0.0f64);
        for &x in xs {
            // two-sum(hi, x)
            let s = hi + x;
            let bb = s - hi;
            let err = (hi - (s - bb)) + (x - bb);
            hi = s;
            lo += err;
        }
        hi + lo
    }

    #[test]
    fn tree_shape_is_a_pure_function_of_length() {
        // Concatenation exposes the visit order: it must be ascending for
        // every length, i.e. the tree only regroups, never reorders.
        for n in 0..40usize {
            let items: Vec<Vec<u64>> = (0..n as u64).map(|i| vec![i]).collect();
            let folded = pairwise_fold(items);
            assert!(folded.iter().copied().eq(0..n as u64), "n={n}");
        }
        assert_eq!(largest_pow2_below(2), 1);
        assert_eq!(largest_pow2_below(3), 2);
        assert_eq!(largest_pow2_below(8), 4);
        assert_eq!(largest_pow2_below(9), 8);
    }

    #[test]
    fn compensated_sum_tracks_the_dd_oracle_on_ill_conditioned_series() {
        // Alternating large/small terms: naive summation loses the small
        // terms entirely; Neumaier + fixed tree must stay at dd accuracy.
        let mut xs = Vec::new();
        for i in 0..10_000 {
            xs.push(1e16);
            xs.push(0.371 + f64::from(i % 7));
            xs.push(-1e16);
        }
        let naive: f64 = xs.iter().sum();
        let det = det_sum(&xs);
        let oracle = dd_sum(&xs);
        let det_err = ((det - oracle) / oracle).abs();
        let naive_err = ((naive - oracle) / oracle).abs();
        assert!(det_err < 1e-12, "det_sum rel err {det_err}");
        assert!(
            naive_err > det_err * 1e3 || naive_err > 1e-10,
            "the series must actually be ill-conditioned (naive rel err {naive_err})"
        );
    }

    #[test]
    fn compensated_partials_merge_associatively_enough_for_trees() {
        let xs: Vec<f64> = (0..1000).map(|i| 1.0 / f64::from(i + 1)).collect();
        let whole = det_sum(&xs);
        // Same values, different BLOCK boundary simulation: accumulate as
        // one compensated stream, then as merged halves via Reduce.
        let a = xs[..500]
            .iter()
            .fold(Compensated::zero(), |c, &x| c.accumulate(x));
        let b = xs[500..]
            .iter()
            .fold(Compensated::zero(), |c, &x| c.accumulate(x));
        let merged = a.merge(b).value();
        assert!((whole - merged).abs() <= 1e-12 * whole.abs());
    }

    #[test]
    fn argmin_ties_break_to_lowest_index_and_nans_are_total_ordered() {
        let xs = [3.0, 1.0, 5.0, 1.0, 2.0];
        assert_eq!(det_argmin(&xs), Some(1), "tie -> lowest index");
        assert_eq!(det_argmax(&[2.0, 7.0, 7.0]), Some(1));
        let with_nan = [f64::NAN, 2.0, -1.0];
        assert_eq!(det_min(&with_nan), Some(-1.0));
        assert_eq!(
            det_max(&with_nan).map(f64::is_nan),
            Some(true),
            "NaN after +inf"
        );
        assert_eq!(det_argmin(&with_nan), Some(2));
        assert_eq!(det_argmin(&[]), None);
    }

    #[test]
    fn prefix_sum_is_inclusive_and_deterministic() {
        let p = det_prefix_sum(&[1.0, 2.0, 3.0]);
        assert_eq!(p, vec![1.0, 3.0, 6.0]);
        assert!(det_prefix_sum(&[]).is_empty());
    }

    #[test]
    fn audit_passes_exact_accumulators_and_localizes_order_bugs() {
        const SEED: u64 = 0xA0D1_2026;
        // Integer-valued sums are exact in f64 at this scale: order-free.
        let ints: Vec<f64> = (0..200).map(f64::from).collect();
        assert!(audit_accumulator(&ints, SEED, |a, x| a + x).is_ok());
        // The classic arrival-order bug: naive float accumulation of an
        // ill-conditioned series. The audit must catch AND localize it.
        let mut nasty = Vec::new();
        for _ in 0..50 {
            nasty.push(1e16);
            nasty.push(1.0);
            nasty.push(-1e16);
        }
        let bug = audit_accumulator(&nasty, SEED, |a, x| a + x)
            .expect_err("naive accumulation must be flagged");
        assert!(bug.witness >= 2 && bug.witness <= nasty.len(), "{bug}");
        assert!(bug.to_string().contains("logical index"), "teaching text");
    }
}
