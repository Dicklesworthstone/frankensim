//! Deterministic closed-AABB candidate generation with explicit work budgets.
//!
//! The caller owns the validity of each spatial or spacetime enclosure. This
//! module returns every overlapping pair of admitted boxes; it does not turn
//! sampled poses into a continuous-motion certificate. The certified-motion
//! wrapper and numerical drift paths share this spatial operation.
//!
//! Exact one-dimensional overlap counts select the least expensive sweep axis
//! before pair enumeration. Setup is `O(n log n)`, the scan is `O(n + k)` for
//! `k` overlapping intervals on that axis, and canonical output sorting is
//! `O(p log p)` for `p` retained candidates. These are operation-count bounds,
//! not hardware performance claims. Closed touching intervals are retained.

use fs_geom::{Aabb, Point3};

/// Hard ceiling on input boxes, including the certified-motion wrapper.
pub const MAX_SWEEP_ITEMS: usize = 1 << 16;

/// A conservative box carrying the caller's stable object identifier.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SweepItem {
    /// Unique identifier; independent of the item's position in the input.
    pub id: u64,
    /// Finite, non-inverted, closed bounds in a shared coordinate frame/unit.
    pub bounds: Aabb,
}

/// One canonical candidate, retaining both identity and input-array addressing.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SweepPair {
    /// Smaller object identifier.
    pub id_a: u64,
    /// Larger object identifier.
    pub id_b: u64,
    /// Input index of `id_a` (not necessarily smaller than `index_b`).
    pub index_a: usize,
    /// Input index of `id_b`.
    pub index_b: usize,
}

/// Per-call bounds on candidate output and pair-comparison work.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SweepConfig {
    /// Maximum retained overlapping pairs. Zero permits a provably empty set.
    pub max_candidates: usize,
    /// Maximum interval-overlapping pairs whose other coordinates are checked.
    /// Setup and output sorting are additional, bounded by input/output sizes.
    pub max_comparisons: usize,
}

/// A complete candidate set and the work performed to obtain it.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct SweepReport {
    /// All closed-box overlaps, sorted by `(id_a, id_b)` without duplicates.
    pub pairs: Vec<SweepPair>,
    /// Pairs compared after passing the selected-axis interval test.
    pub checked_pairs: usize,
    /// All other unordered input pairs, excluded by the selected axis.
    pub pruned_pairs: usize,
    /// Selected coordinate: 0 = x, 1 = y, 2 = z; ties favor x, then y.
    pub sweep_axis: usize,
}

/// Refusals never present a truncated candidate set as complete.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SweepError {
    /// The input exceeds the fixed setup-storage ceiling.
    TooManyItems {
        /// Supplied item count.
        items: usize,
        /// Admitted ceiling.
        max: usize,
    },
    /// Non-finite or inverted public AABB fields.
    InvalidBounds {
        /// Original input index.
        index: usize,
    },
    /// Stable identifiers must be unique within the call.
    DuplicateId {
        /// Repeated identifier.
        id: u64,
    },
    /// A swept-sphere input violates the named finite/nonnegative assumption.
    InvalidSphere(&'static str),
    /// Outward-rounded sphere bounds cannot be represented by finite endpoints.
    NumericalRange,
    /// Pair counts cannot be represented by the platform's `usize`.
    CountOverflow,
    /// Exact preflight count exceeds the caller's pair-work budget.
    ComparisonBudgetExhausted {
        /// Caller cap.
        max_comparisons: usize,
        /// Required checks on the least-overlapping coordinate axis.
        required_comparisons: usize,
        /// Selected coordinate axis.
        sweep_axis: usize,
    },
    /// At least one more candidate exists than the output budget permits.
    CandidateBudgetExhausted {
        /// Caller cap.
        max_candidates: usize,
        /// Canonically sorted retained scan prefix; explicitly incomplete.
        prefix: Vec<SweepPair>,
        /// Checks performed, including the first unretained candidate.
        checked_pairs: usize,
        /// Selected coordinate axis.
        sweep_axis: usize,
    },
    /// Cancellation was observed before publishing a complete report.
    Cancelled,
}

impl core::fmt::Display for SweepError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::TooManyItems { items, max } => {
                write!(f, "{items} boxes exceed the {max}-item broad-phase ceiling")
            }
            Self::InvalidBounds { index } => {
                write!(f, "box {index} has non-finite or inverted bounds")
            }
            Self::DuplicateId { id } => write!(f, "broad-phase object id {id} is repeated"),
            Self::InvalidSphere(reason) => write!(f, "invalid swept sphere: {reason}"),
            Self::NumericalRange => write!(f, "finite outward sphere bounds are unrepresentable"),
            Self::CountOverflow => write!(f, "broad-phase pair count exceeds usize"),
            Self::ComparisonBudgetExhausted {
                max_comparisons,
                required_comparisons,
                sweep_axis,
            } => write!(
                f,
                "axis {sweep_axis} requires {required_comparisons} pair checks, exceeding budget {max_comparisons}"
            ),
            Self::CandidateBudgetExhausted { max_candidates, .. } => write!(
                f,
                "candidate budget {max_candidates} exhausted; retained prefix is incomplete"
            ),
            Self::Cancelled => write!(f, "cancelled during broad-phase candidate generation"),
        }
    }
}

impl core::error::Error for SweepError {}

/// Enclose a sphere translated along the straight segment `start` to `end`.
///
/// `radius` and `padding` use the same coordinate unit as the centres and must
/// be finite and nonnegative. Each arithmetic bound is rounded outward,
/// including the sum of positive radius and padding. Finite extreme inputs
/// whose enclosure would require an infinite endpoint refuse.
///
/// This encloses straight centre motion only. Curved/accelerating trajectories,
/// state uncertainty, and geometry conversion errors need their own enclosure
/// or a justified padding supplied by the caller.
pub fn swept_sphere_aabb(
    start: [f64; 3],
    end: [f64; 3],
    radius: f64,
    padding: f64,
) -> Result<Aabb, SweepError> {
    if !start.into_iter().chain(end).all(f64::is_finite) {
        return Err(SweepError::InvalidSphere("centres must be finite"));
    }
    if !radius.is_finite() || radius < 0.0 {
        return Err(SweepError::InvalidSphere(
            "radius must be finite and nonnegative",
        ));
    }
    if !padding.is_finite() || padding < 0.0 {
        return Err(SweepError::InvalidSphere(
            "padding must be finite and nonnegative",
        ));
    }
    let extent = if radius == 0.0 {
        padding
    } else if padding == 0.0 {
        radius
    } else {
        (radius + padding).next_up()
    };
    let lo =
        core::array::from_fn::<_, 3, _>(|axis| (start[axis].min(end[axis]) - extent).next_down());
    let hi =
        core::array::from_fn::<_, 3, _>(|axis| (start[axis].max(end[axis]) + extent).next_up());
    if !extent.is_finite() || !lo.into_iter().chain(hi).all(f64::is_finite) {
        return Err(SweepError::NumericalRange);
    }
    Ok(Aabb {
        min: Point3::new(lo[0], lo[1], lo[2]),
        max: Point3::new(hi[0], hi[1], hi[2]),
    })
}

/// Return every overlap of the admitted closed boxes under explicit budgets.
///
/// Axis selection and traversal are deterministic under input permutation when
/// identifiers and boxes stay attached; only the returned array indices change.
/// Setup uses `O(n)` auxiliary indices; the scan needs no per-pair scratch.
/// Candidate storage grows only when a retained output pair is added.
/// Cancellation is checked during setup/scanning and before/after sorts, whose
/// lengths are bounded by the admitted input/output sizes.
pub fn sweep_candidates(
    items: &[SweepItem],
    config: SweepConfig,
    mut cancelled: impl FnMut() -> bool,
) -> Result<SweepReport, SweepError> {
    let mut pairs = Vec::new();
    let scan = visit_sweep(items, config.max_comparisons, &mut cancelled, |pair| {
        if pairs.len() == config.max_candidates {
            false
        } else {
            pairs.push(pair);
            true
        }
    })?;
    checkpoint(&mut cancelled)?;
    pairs.sort_unstable_by_key(|pair| (pair.id_a, pair.id_b));
    checkpoint(&mut cancelled)?;
    if !scan.completed {
        return Err(SweepError::CandidateBudgetExhausted {
            max_candidates: config.max_candidates,
            prefix: pairs,
            checked_pairs: scan.checked_pairs,
            sweep_axis: scan.sweep_axis,
        });
    }
    Ok(SweepReport {
        pairs,
        checked_pairs: scan.checked_pairs,
        pruned_pairs: scan.pruned_pairs,
        sweep_axis: scan.sweep_axis,
    })
}

/// The certified-motion wrapper needs every overflow pair, whereas the bounded
/// numerical API stops at the first overflow. Both use this same scan.
pub(crate) struct SweepVisitReport {
    pub(crate) checked_pairs: usize,
    pub(crate) pruned_pairs: usize,
    pub(crate) sweep_axis: usize,
    pub(crate) completed: bool,
}

pub(crate) fn visit_sweep(
    items: &[SweepItem],
    max_comparisons: usize,
    cancelled: &mut impl FnMut() -> bool,
    mut visit: impl FnMut(SweepPair) -> bool,
) -> Result<SweepVisitReport, SweepError> {
    let prepared = prepare(items, max_comparisons, cancelled)?;
    let mut checked_pairs = 0usize;
    for (rank, &i) in prepared.order.iter().enumerate() {
        checkpoint(cancelled)?;
        let end = coordinate(items[i].bounds.max, prepared.axis);
        for &j in &prepared.order[rank + 1..] {
            if coordinate(items[j].bounds.min, prepared.axis) > end {
                break;
            }
            if checked_pairs & 255 == 0 {
                checkpoint(cancelled)?;
            }
            checked_pairs += 1;
            if transverse_overlap(items[i].bounds, items[j].bounds, prepared.axis) {
                let (a, b) = if items[i].id < items[j].id {
                    (i, j)
                } else {
                    (j, i)
                };
                if !visit(SweepPair {
                    id_a: items[a].id,
                    id_b: items[b].id,
                    index_a: a,
                    index_b: b,
                }) {
                    checkpoint(cancelled)?;
                    return Ok(SweepVisitReport {
                        checked_pairs,
                        pruned_pairs: prepared.total_pairs - prepared.comparisons,
                        sweep_axis: prepared.axis,
                        completed: false,
                    });
                }
            }
        }
    }
    checkpoint(cancelled)?;
    debug_assert_eq!(checked_pairs, prepared.comparisons);
    Ok(SweepVisitReport {
        checked_pairs,
        pruned_pairs: prepared.total_pairs - checked_pairs,
        sweep_axis: prepared.axis,
        completed: true,
    })
}

struct PreparedSweep {
    order: Vec<usize>,
    axis: usize,
    comparisons: usize,
    total_pairs: usize,
}

fn prepare(
    items: &[SweepItem],
    max_comparisons: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<PreparedSweep, SweepError> {
    checkpoint(cancelled)?;
    if items.len() > MAX_SWEEP_ITEMS {
        return Err(SweepError::TooManyItems {
            items: items.len(),
            max: MAX_SWEEP_ITEMS,
        });
    }
    for (index, item) in items.iter().enumerate() {
        if index & 255 == 0 {
            checkpoint(cancelled)?;
        }
        if !item.bounds.is_finite() {
            return Err(SweepError::InvalidBounds { index });
        }
    }
    let n = items.len();
    let total_pairs = n
        .checked_mul(n.saturating_sub(1))
        .ok_or(SweepError::CountOverflow)?
        / 2;
    let mut minima: Vec<usize> = (0..n).collect();
    checkpoint(cancelled)?;
    minima.sort_unstable_by_key(|&index| items[index].id);
    checkpoint(cancelled)?;
    for (index, pair) in minima.windows(2).enumerate() {
        if index & 255 == 0 {
            checkpoint(cancelled)?;
        }
        if items[pair[0]].id == items[pair[1]].id {
            return Err(SweepError::DuplicateId {
                id: items[pair[0]].id,
            });
        }
    }
    let mut maxima = minima.clone();
    let mut best = PreparedSweep {
        order: Vec::new(),
        axis: 0,
        comparisons: usize::MAX,
        total_pairs,
    };
    for axis in 0..3 {
        checkpoint(cancelled)?;
        minima.sort_unstable_by(|&i, &j| {
            coordinate(items[i].bounds.min, axis)
                .total_cmp(&coordinate(items[j].bounds.min, axis))
                .then_with(|| items[i].id.cmp(&items[j].id))
        });
        checkpoint(cancelled)?;
        maxima.sort_unstable_by(|&i, &j| {
            coordinate(items[i].bounds.max, axis)
                .total_cmp(&coordinate(items[j].bounds.max, axis))
                .then_with(|| items[i].id.cmp(&items[j].id))
        });
        checkpoint(cancelled)?;
        let mut expired = 0usize;
        let mut comparisons = 0usize;
        for (rank, &index) in minima.iter().enumerate() {
            if rank & 255 == 0 {
                checkpoint(cancelled)?;
            }
            let start = coordinate(items[index].bounds.min, axis);
            // Strict expiry retains touching closed intervals, including
            // signed-zero endpoints. Every expired interval started earlier.
            while expired < n && coordinate(items[maxima[expired]].bounds.max, axis) < start {
                if expired & 255 == 0 {
                    checkpoint(cancelled)?;
                }
                expired += 1;
            }
            let active = rank.checked_sub(expired).ok_or(SweepError::CountOverflow)?;
            comparisons = comparisons
                .checked_add(active)
                .ok_or(SweepError::CountOverflow)?;
        }
        if comparisons < best.comparisons {
            best.order.clone_from(&minima);
            best.axis = axis;
            best.comparisons = comparisons;
        }
        if comparisons == 0 {
            break;
        }
    }
    checkpoint(cancelled)?;
    if best.comparisons > max_comparisons {
        return Err(SweepError::ComparisonBudgetExhausted {
            max_comparisons,
            required_comparisons: best.comparisons,
            sweep_axis: best.axis,
        });
    }
    Ok(best)
}

fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> Result<(), SweepError> {
    if cancelled() {
        Err(SweepError::Cancelled)
    } else {
        Ok(())
    }
}

fn coordinate(point: Point3, axis: usize) -> f64 {
    match axis {
        0 => point.x,
        1 => point.y,
        _ => point.z,
    }
}

fn transverse_overlap(a: Aabb, b: Aabb, axis: usize) -> bool {
    match axis {
        0 => a.min.y <= b.max.y && b.min.y <= a.max.y && a.min.z <= b.max.z && b.min.z <= a.max.z,
        1 => a.min.x <= b.max.x && b.min.x <= a.max.x && a.min.z <= b.max.z && b.min.z <= a.max.z,
        _ => a.min.x <= b.max.x && b.min.x <= a.max.x && a.min.y <= b.max.y && b.min.y <= a.max.y,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const UNLIMITED: SweepConfig = SweepConfig {
        max_candidates: usize::MAX,
        max_comparisons: usize::MAX,
    };

    fn item(id: u64, lo: [f64; 3], hi: [f64; 3]) -> SweepItem {
        SweepItem {
            id,
            bounds: Aabb {
                min: Point3::new(lo[0], lo[1], lo[2]),
                max: Point3::new(hi[0], hi[1], hi[2]),
            },
        }
    }

    fn identifiers(report: &SweepReport) -> Vec<(u64, u64)> {
        report
            .pairs
            .iter()
            .map(|pair| (pair.id_a, pair.id_b))
            .collect()
    }

    fn brute_force(items: &[SweepItem]) -> Vec<(u64, u64)> {
        let mut pairs = Vec::new();
        for i in 0..items.len() {
            for j in i + 1..items.len() {
                let a = items[i].bounds;
                let b = items[j].bounds;
                if a.min.x <= b.max.x
                    && b.min.x <= a.max.x
                    && a.min.y <= b.max.y
                    && b.min.y <= a.max.y
                    && a.min.z <= b.max.z
                    && b.min.z <= a.max.z
                {
                    pairs.push((items[i].id.min(items[j].id), items[i].id.max(items[j].id)));
                }
            }
        }
        pairs.sort_unstable();
        pairs
    }

    // G0/G3/G5: independent all-pairs geometry, permutation and replay checks.
    #[test]
    fn complete_candidates_match_brute_force_under_permutation() {
        let mut bits = 0x6d4d_4149_4e42_5048u64;
        let mut next = || {
            bits = bits
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (bits >> 32) as u32
        };
        for case in 0usize..240 {
            let n = case % 37;
            let mut items = Vec::new();
            for index in 0..n {
                let lo = core::array::from_fn(|_| (next() % 97) as f64 * 0.125 - 6.0);
                let hi = core::array::from_fn(|axis| lo[axis] + (next() % 41) as f64 * 0.125);
                items.push(item(u64::MAX - index as u64 * 17, lo, hi));
            }
            let expected = brute_force(&items);
            let report = sweep_candidates(&items, UNLIMITED, || false).unwrap();
            assert_eq!(identifiers(&report), expected, "case {case}");
            assert_eq!(
                report,
                sweep_candidates(&items, UNLIMITED, || false).unwrap()
            );
            assert_eq!(
                report.checked_pairs + report.pruned_pairs,
                n * n.saturating_sub(1) / 2
            );
            let counts: [usize; 3] = core::array::from_fn(|axis| {
                let mut count = 0;
                for i in 0..n {
                    for j in i + 1..n {
                        if coordinate(items[i].bounds.min, axis)
                            <= coordinate(items[j].bounds.max, axis)
                            && coordinate(items[j].bounds.min, axis)
                                <= coordinate(items[i].bounds.max, axis)
                        {
                            count += 1;
                        }
                    }
                }
                count
            });
            assert_eq!(report.checked_pairs, *counts.iter().min().unwrap());
            assert_eq!(
                report.sweep_axis,
                counts
                    .iter()
                    .position(|count| *count == report.checked_pairs)
                    .unwrap()
            );
            items.reverse();
            if n > 0 {
                items.rotate_left(case % n);
            }
            let reordered = sweep_candidates(&items, UNLIMITED, || false).unwrap();
            assert_eq!(identifiers(&reordered), expected);
            assert_eq!(reordered.sweep_axis, report.sweep_axis);
            for pair in reordered.pairs {
                assert_eq!(items[pair.index_a].id, pair.id_a);
                assert_eq!(items[pair.index_b].id, pair.id_b);
            }
        }
    }

    #[test]
    fn closed_touching_nested_boxes_and_signed_zero_are_retained() {
        let items = [
            item(70, [-0.0; 3], [0.0; 3]),
            item(10, [0.0; 3], [1.0; 3]),
            item(40, [-1.0; 3], [-0.0; 3]),
            item(30, [-2.0; 3], [2.0; 3]),
            item(20, [3.0; 3], [4.0; 3]),
        ];
        let report = sweep_candidates(&items, UNLIMITED, || false).unwrap();
        assert_eq!(
            identifiers(&report),
            vec![(10, 30), (10, 40), (10, 70), (30, 40), (30, 70), (40, 70)]
        );
        assert_eq!(report.sweep_axis, 0);
        assert_eq!(report.checked_pairs, 6);
        assert_eq!(report.pruned_pairs, 4);
    }

    #[test]
    fn separated_vertical_scene_uses_zero_pair_comparisons() {
        let items: Vec<_> = (0..4096)
            .map(|i| {
                item(
                    i as u64,
                    [-1.0, 4.0 * i as f64, -1.0],
                    [1.0, 4.0 * i as f64 + 1.0, 1.0],
                )
            })
            .collect();
        let report = sweep_candidates(
            &items,
            SweepConfig {
                max_candidates: 0,
                max_comparisons: 0,
            },
            || false,
        )
        .unwrap();
        assert!(report.pairs.is_empty());
        assert_eq!(report.sweep_axis, 1);
        assert_eq!(report.checked_pairs, 0);
        assert_eq!(report.pruned_pairs, 8_386_560);
    }

    #[test]
    fn budgets_refuse_before_work_or_at_the_first_unretained_pair() {
        let items: Vec<_> = (0..4).map(|i| item(9 - i, [0.0; 3], [1.0; 3])).collect();
        let mut visited = 0;
        let error = visit_sweep(&items, 5, &mut || false, |_| {
            visited += 1;
            true
        })
        .err()
        .unwrap();
        assert_eq!(visited, 0);
        assert_eq!(
            error,
            SweepError::ComparisonBudgetExhausted {
                max_comparisons: 5,
                required_comparisons: 6,
                sweep_axis: 0
            }
        );
        let error = sweep_candidates(
            &items,
            SweepConfig {
                max_candidates: 2,
                max_comparisons: 6,
            },
            || false,
        )
        .unwrap_err();
        let SweepError::CandidateBudgetExhausted {
            max_candidates,
            prefix,
            checked_pairs,
            sweep_axis,
        } = error
        else {
            panic!("expected candidate refusal")
        };
        assert_eq!(
            (max_candidates, prefix.len(), checked_pairs, sweep_axis),
            (2, 2, 3, 0)
        );
        assert!(
            prefix
                .windows(2)
                .all(|pair| (pair[0].id_a, pair[0].id_b) < (pair[1].id_a, pair[1].id_b))
        );
        let all = sweep_candidates(
            &items,
            SweepConfig {
                max_candidates: 6,
                max_comparisons: 6,
            },
            || false,
        )
        .unwrap();
        assert_eq!(all.pairs.len(), 6);
        // The certified wrapper's visitor can retain every overflow pair.
        let mut retained = Vec::new();
        let mut overflow = Vec::new();
        let scan = visit_sweep(&items, 6, &mut || false, |pair| {
            if retained.len() < 2 {
                retained.push(pair);
            } else {
                overflow.push(pair);
            }
            true
        })
        .unwrap();
        assert!(scan.completed);
        assert_eq!(
            (retained.len(), overflow.len(), scan.checked_pairs),
            (2, 4, 6)
        );
        retained.extend(overflow);
        retained.sort_unstable_by_key(|pair| (pair.id_a, pair.id_b));
        assert_eq!(retained, all.pairs);
    }

    #[test]
    fn outward_sphere_bounds_cover_the_entire_straight_path() {
        for scale in [1e-300, 1.0, 1e300] {
            let start = [-3.0 * scale, 2.0 * scale, 0.0];
            let end = [4.0 * scale, -2.0 * scale, scale];
            let radius = 0.25 * scale;
            let padding = 0.125 * scale;
            let bounds = swept_sphere_aabb(start, end, radius, padding).unwrap();
            for step in 0..=20 {
                let t = step as f64 / 20.0;
                let center: [f64; 3] =
                    core::array::from_fn(|axis| (1.0 - t) * start[axis] + t * end[axis]);
                for corner in 0..8 {
                    let point: [f64; 3] = core::array::from_fn(|axis| {
                        let sign = if corner & (1 << axis) == 0 { -1.0 } else { 1.0 };
                        center[axis] + sign * (radius + padding)
                    });
                    assert!(bounds.contains(Point3::new(point[0], point[1], point[2])));
                }
            }
        }
        let bounds = swept_sphere_aabb([0.0; 3], [0.0; 3], 0.0, 0.0).unwrap();
        assert_eq!(bounds.min.x, -f64::from_bits(1));
        assert_eq!(bounds.max.x, f64::from_bits(1));
        let bounds = swept_sphere_aabb([1.0; 3], [1.0; 3], 1.0, f64::EPSILON / 4.0).unwrap();
        assert!(bounds.max.x > 2.0);
        assert!(bounds.min.x < 0.0);
        assert_eq!(
            swept_sphere_aabb([f64::MAX; 3], [f64::MAX; 3], 0.0, 0.0),
            Err(SweepError::NumericalRange)
        );
    }

    #[test]
    fn malformed_inputs_and_body_ceiling_refuse_without_partial_candidates() {
        for bounds in [
            item(0, [f64::NAN, 0.0, 0.0], [1.0; 3]),
            item(0, [0.0; 3], [f64::INFINITY, 1.0, 1.0]),
            item(0, [2.0, 0.0, 0.0], [1.0; 3]),
        ] {
            assert_eq!(
                sweep_candidates(&[bounds], UNLIMITED, || false),
                Err(SweepError::InvalidBounds { index: 0 })
            );
        }
        let same = item(u64::MAX, [0.0; 3], [0.0; 3]);
        assert_eq!(
            sweep_candidates(&[same, same], UNLIMITED, || false),
            Err(SweepError::DuplicateId { id: u64::MAX })
        );
        for (radius, padding) in [
            (-1.0, 0.0),
            (0.0, -1.0),
            (f64::NAN, 0.0),
            (0.0, f64::INFINITY),
        ] {
            assert!(matches!(
                swept_sphere_aabb([0.0; 3], [0.0; 3], radius, padding),
                Err(SweepError::InvalidSphere(_))
            ));
        }
        assert!(matches!(
            swept_sphere_aabb([f64::NAN; 3], [0.0; 3], 1.0, 0.0),
            Err(SweepError::InvalidSphere(_))
        ));
        assert_eq!(
            swept_sphere_aabb([0.0; 3], [0.0; 3], f64::MAX, f64::MAX),
            Err(SweepError::NumericalRange)
        );
        let mut items: Vec<_> = (0..MAX_SWEEP_ITEMS)
            .map(|i| {
                item(
                    i as u64,
                    [2.0 * i as f64, 0.0, 0.0],
                    [2.0 * i as f64, 1.0, 1.0],
                )
            })
            .collect();
        let report = sweep_candidates(
            &items,
            SweepConfig {
                max_candidates: 0,
                max_comparisons: 0,
            },
            || false,
        )
        .unwrap();
        assert_eq!(report.checked_pairs, 0);
        items.push(item(MAX_SWEEP_ITEMS as u64, [0.0; 3], [1.0; 3]));
        assert_eq!(
            sweep_candidates(&items, UNLIMITED, || false),
            Err(SweepError::TooManyItems {
                items: MAX_SWEEP_ITEMS + 1,
                max: MAX_SWEEP_ITEMS
            })
        );
    }

    #[test]
    fn cancellation_never_publishes_an_incomplete_report() {
        assert_eq!(
            sweep_candidates(&[], UNLIMITED, || true),
            Err(SweepError::Cancelled)
        );
        let items: Vec<_> = (0..256).map(|i| item(i, [0.0; 3], [1.0; 3])).collect();
        let mut polls = 0;
        let result = sweep_candidates(&items, UNLIMITED, || {
            polls += 1;
            polls == 30
        });
        assert_eq!(result, Err(SweepError::Cancelled));
        assert_eq!(polls, 30);
        let empty = sweep_candidates(
            &[],
            SweepConfig {
                max_candidates: 0,
                max_comparisons: 0,
            },
            || false,
        )
        .unwrap();
        assert_eq!(empty.pairs.len(), 0);
        assert_eq!(
            (empty.checked_pairs, empty.pruned_pairs, empty.sweep_axis),
            (0, 0, 0)
        );
    }
}
