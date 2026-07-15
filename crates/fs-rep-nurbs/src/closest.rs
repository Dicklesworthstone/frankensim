//! Measured closest-point bracket estimates: BEST-FIRST branch-and-bound over
//! rational Bézier segments (in exact arithmetic the convex-hull property
//! supplies lower bounds; positive weights survive de Casteljau splitting).
//! The current Cartesian hull is evaluated in ordinary f64. Splits are LOCAL —
//! one segment per iteration — so heap selection costs O(log S) per split,
//! and the split junction point is a free upper-bound sample. Boxes are
//! heuristically expanded by one ULP against f64 rounding. Dense-oracle
//! conformance is useful evidence, but ordinary Cartesian division, distance
//! arithmetic, and evaluation are not outward-rounded; `[lower, upper]` is not
//! a rigorous enclosure until the interval/Taylor upgrade lands.

use crate::NurbsError;
use crate::basis::AdmittedKnotVector;
use crate::curve::{AdmittedNurbsCurve, BezierConversionPlan, NurbsCurve};
use crate::surface::{AdmittedNurbsSurface, NurbsSurface};
use fs_math::{det, next_down, next_up};
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::mem::size_of;

/// Defensive ceiling for the legacy allocation-bearing subdivision path.
/// Caller-owned work/memory budgets and cancellation belong to the successor
/// certifying API; this cap only prevents an unbounded `u32` request here.
pub(crate) const CLOSEST_MAX_SPLITS: u32 = 1_048_576;

/// Defensive ceiling for conservative initial validation and stage-faithful
/// knot-insertion, expanded-grid, run-scan, and queue-seeding work before the
/// split counter starts. This closes the `max_splits=0` admission gap on large
/// public spline structures.
pub(crate) const CLOSEST_MAX_BASE_WORK_UNITS: u128 = 16_777_216;

/// Defensive aggregate split-work ceiling for the legacy APIs.
const CLOSEST_MAX_SPLIT_WORK_UNITS: u128 = 1_073_741_824;

/// Defensive retained-payload ceiling. Curve admission composes source,
/// conversion, search frontier, and post-release polish phases; the legacy
/// surface path currently applies it to frontier/scratch payload.
const CLOSEST_MAX_RETAINED_BYTES: u128 = 256 * 1024 * 1024;

// Keep the source-snapshot charge aligned with the admitted curve/knot
// validation envelopes. This is an accounting coefficient, not a timing
// claim for any particular machine.
const CLOSEST_CURVE_SOURCE_SCAN_WORK_PER_ENTRY: u128 = 16;
const CLOSEST_SURFACE_SOURCE_SCAN_WORK_PER_CONTROL: u128 = 16;

/// Deterministic minimum-priority entry. `BinaryHeap` is a max-heap, so the
/// comparisons are reversed: lower bound first, then lower logical ID. IDs
/// are unique among resident entries and are reused only when the same popped
/// unsplittable leaf is reinserted.
struct MinEntry<T> {
    key: f64,
    logical_id: u64,
    value: T,
}

impl<T> PartialEq for MinEntry<T> {
    fn eq(&self, other: &Self) -> bool {
        self.key.to_bits() == other.key.to_bits() && self.logical_id == other.logical_id
    }
}

impl<T> Eq for MinEntry<T> {}

impl<T> PartialOrd for MinEntry<T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<T> Ord for MinEntry<T> {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .key
            .total_cmp(&self.key)
            .then_with(|| other.logical_id.cmp(&self.logical_id))
    }
}

/// A measured distance-bracket estimate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DistanceBracketEstimate {
    /// Convex-hull lower estimate with heuristic f64 inflation.
    pub lower: f64,
    /// Achieved evaluated-point distance estimate.
    pub upper: f64,
    /// Parameter of the best found point (curve: t; surface: (u, v)).
    pub param: [f64; 2],
    /// Branch-and-bound splits spent.
    pub iterations: u32,
}

pub(crate) fn norm3(value: [f64; 3]) -> f64 {
    if value.iter().any(|component| !component.is_finite()) {
        return f64::INFINITY;
    }
    let scale = value
        .iter()
        .fold(0.0f64, |largest, component| largest.max(component.abs()));
    if scale == 0.0 {
        return 0.0;
    }
    if !scale.is_finite() {
        return f64::INFINITY;
    }
    let normalized_square_sum: f64 = value
        .iter()
        .map(|component| (component / scale).powi(2))
        .sum();
    scale * det::sqrt(normalized_square_sum)
}

fn nonempty_span_count(knots: AdmittedKnotVector<'_, f64>) -> u128 {
    knots
        .knots()
        .windows(2)
        .filter(|pair| pair[1] > pair[0])
        .count() as u128
}

fn bezier_insertion_count(knots: AdmittedKnotVector<'_, f64>) -> Result<u128, NurbsError> {
    let (lo, hi) = knots.domain();
    let entries = knots.knots();
    let degree = knots.degree();
    let mut insertions = 0u128;
    let mut run_start = 0usize;
    while run_start < entries.len() {
        let knot = entries[run_start];
        let mut run_end = run_start + 1;
        while run_end < entries.len() && entries[run_end] == knot {
            run_end += 1;
        }
        let multiplicity = run_end - run_start;
        if knot > lo && knot < hi && multiplicity < degree {
            insertions = insertions
                .checked_add((degree - multiplicity) as u128)
                .ok_or_else(|| NurbsError::Domain {
                    what: "Bézier insertion-count accounting overflows u128".to_string(),
                })?;
        }
        run_start = run_end;
    }
    Ok(insertions)
}

fn checked_product(values: &[u128], what: &str) -> Result<u128, NurbsError> {
    values.iter().try_fold(1u128, |acc, value| {
        acc.checked_mul(*value).ok_or_else(|| NurbsError::Domain {
            what: format!("{what} overflows u128 work accounting"),
        })
    })
}

fn checked_sum(values: &[u128], what: &str) -> Result<u128, NurbsError> {
    values.iter().try_fold(0u128, |acc, value| {
        acc.checked_add(*value).ok_or_else(|| NurbsError::Domain {
            what: format!("{what} overflows u128 work accounting"),
        })
    })
}

fn curve_source_admission_work<const DIM: usize>(
    curve: &NurbsCurve<f64, DIM>,
) -> Result<u128, NurbsError> {
    let knot_work = (curve.knots.knots.len() as u128)
        .checked_mul(CLOSEST_CURVE_SOURCE_SCAN_WORK_PER_ENTRY)
        .and_then(|work| work.checked_add(curve.knots.degree as u128))
        .ok_or_else(|| NurbsError::Domain {
            what: "curve source-knot admission work overflows u128".to_string(),
        })?;
    let control_work = (curve.cpw.len() as u128)
        .checked_mul(CLOSEST_CURVE_SOURCE_SCAN_WORK_PER_ENTRY)
        .ok_or_else(|| NurbsError::Domain {
            what: "curve source-control admission work overflows u128".to_string(),
        })?;
    checked_sum(&[knot_work, control_work], "curve source admission work")
}

fn surface_source_admission_work(surface: &NurbsSurface<f64>) -> Result<u128, NurbsError> {
    let controls = (surface.knots_u.control_count() as u128)
        .checked_mul(surface.knots_v.control_count() as u128)
        .and_then(|count| count.checked_mul(CLOSEST_SURFACE_SOURCE_SCAN_WORK_PER_CONTROL))
        .ok_or_else(|| NurbsError::Domain {
            what: "surface source-control admission work overflows u128".to_string(),
        })?;
    checked_sum(
        &[
            surface.knots_u.validation_work()?,
            surface.knots_v.validation_work()?,
            controls,
        ],
        "surface source admission work",
    )
}

fn curve_base_work_units<const DIM: usize>(
    curve: AdmittedNurbsCurve<'_, f64, DIM>,
    source_admission_work: u128,
    conversion: BezierConversionPlan,
    seed_leaves: u128,
    order: u128,
    polish_work: u128,
    seed_heap_work: u128,
) -> Result<u128, NurbsError> {
    let seed_control_visits = checked_product(&[seed_leaves, order], "curve patch seeding")?;
    let seed_work = checked_product(
        &[seed_control_visits, 16],
        "curve seed copy, hull, and queue work",
    )?;
    checked_sum(
        &[
            source_admission_work,
            conversion.work_units,
            curve.knots().knots().len() as u128,
            seed_work,
            seed_heap_work,
            seed_leaves,
            polish_work,
        ],
        "curve base work",
    )
}

#[derive(Debug, Clone, Copy)]
struct SurfaceBasePlan {
    work_units: u128,
    seed_leaves: u128,
    order_u: u128,
    order_v: u128,
}

pub(crate) fn surface_base_work_units(surface: &NurbsSurface<f64>) -> Result<u128, NurbsError> {
    let source_admission_work = surface_source_admission_work(surface)?;
    enforce_base_work(source_admission_work, "surface source admission")?;
    let base_plan = surface_base_plan_from_admitted(surface.admit()?, source_admission_work)?;
    let (_, _, worst_heap_height) = surface_queue_shape(base_plan.seed_leaves, CLOSEST_MAX_SPLITS)?;
    let seed_heap_work = surface_seed_heap_work(base_plan.seed_leaves, worst_heap_height)?;
    checked_sum(
        &[base_plan.work_units, seed_heap_work],
        "surface base and worst-case seed-heap work",
    )
}

fn surface_base_plan_from_admitted(
    surface: AdmittedNurbsSurface<'_, f64>,
    source_admission_work: u128,
) -> Result<SurfaceBasePlan, NurbsError> {
    let knots_u = surface.knots_u();
    let knots_v = surface.knots_v();
    let expected_u = knots_u.control_count();
    let expected_v = knots_v.control_count();
    let spans_u = nonempty_span_count(knots_u);
    let spans_v = nonempty_span_count(knots_v);
    let degree_u = knots_u.degree() as u128;
    let degree_v = knots_v.degree() as u128;
    let order_u = degree_u.checked_add(1).ok_or_else(|| NurbsError::Domain {
        what: "closest-surface u-order accounting overflows u128".to_string(),
    })?;
    let order_v = degree_v.checked_add(1).ok_or_else(|| NurbsError::Domain {
        what: "closest-surface v-order accounting overflows u128".to_string(),
    })?;
    let seed_leaves = spans_u
        .checked_mul(spans_v)
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-surface seed-leaf accounting overflows u128".to_string(),
        })?;
    let insertions_u = bezier_insertion_count(knots_u)?;
    let insertions_v = bezier_insertion_count(knots_v)?;
    let expanded_u = (expected_u as u128)
        .checked_add(insertions_u)
        .ok_or_else(|| NurbsError::Domain {
            what: "surface u Bézier control-count accounting overflows u128".to_string(),
        })?;
    let expanded_v = (expected_v as u128)
        .checked_add(insertions_v)
        .ok_or_else(|| NurbsError::Domain {
            what: "surface v Bézier control-count accounting overflows u128".to_string(),
        })?;
    let expanded_knots_u = (knots_u.knots().len() as u128)
        .checked_add(insertions_u)
        .ok_or_else(|| NurbsError::Domain {
            what: "surface u Bézier knot-count accounting overflows u128".to_string(),
        })?;
    let expanded_knots_v = (knots_v.knots().len() as u128)
        .checked_add(insertions_v)
        .ok_or_else(|| NurbsError::Domain {
            what: "surface v Bézier knot-count accounting overflows u128".to_string(),
        })?;
    let expanded_grid = checked_product(&[expanded_u, expanded_v], "surface Bézier grid")?;
    let patch_controls =
        checked_product(&[seed_leaves, order_u, order_v], "surface patch seeding")?;
    let input_controls = checked_product(
        &[
            surface.homogeneous_control_net().len() as u128,
            expected_v as u128,
        ],
        "surface input grid",
    )?;
    let u_insertion_work = checked_product(
        &[
            insertions_u,
            expanded_v,
            checked_sum(
                &[
                    expanded_u
                        .checked_mul(2)
                        .ok_or_else(|| NurbsError::Domain {
                            what: "surface u insertion-grid work overflows u128".to_string(),
                        })?,
                    expanded_knots_u,
                    order_u.checked_mul(8).ok_or_else(|| NurbsError::Domain {
                        what: "surface u insertion-order work overflows u128".to_string(),
                    })?,
                ],
                "surface u insertion stage",
            )?,
        ],
        "surface u Bézier insertion",
    )?;
    let v_insertion_work = checked_product(
        &[
            insertions_v,
            expanded_u,
            checked_sum(
                &[
                    expanded_v
                        .checked_mul(2)
                        .ok_or_else(|| NurbsError::Domain {
                            what: "surface v insertion-grid work overflows u128".to_string(),
                        })?,
                    expanded_knots_v,
                    order_v.checked_mul(8).ok_or_else(|| NurbsError::Domain {
                        what: "surface v insertion-order work overflows u128".to_string(),
                    })?,
                ],
                "surface v insertion stage",
            )?,
        ],
        "surface v Bézier insertion",
    )?;
    let scan_work = checked_product(
        &[
            insertions_u.max(insertions_v) + 1,
            checked_sum(
                &[expanded_knots_u, expanded_knots_v],
                "surface knot-scan extent",
            )?,
        ],
        "surface Bézier run scanning",
    )?;
    let work_units = checked_sum(
        &[
            source_admission_work,
            input_controls,
            knots_u.knots().len() as u128,
            knots_v.knots().len() as u128,
            expanded_grid,
            patch_controls,
            u_insertion_work,
            v_insertion_work,
            scan_work,
        ],
        "surface base work",
    )?;
    Ok(SurfaceBasePlan {
        work_units,
        seed_leaves,
        order_u,
        order_v,
    })
}

fn enforce_base_work(units: u128, kind: &str) -> Result<(), NurbsError> {
    if units > CLOSEST_MAX_BASE_WORK_UNITS {
        return Err(NurbsError::Domain {
            what: format!(
                "closest-{kind} base work {units} exceeds defensive ceiling {CLOSEST_MAX_BASE_WORK_UNITS}"
            ),
        });
    }
    Ok(())
}

fn enforce_retained_bytes(retained_bytes: u128, kind: &str) -> Result<(), NurbsError> {
    if retained_bytes > CLOSEST_MAX_RETAINED_BYTES {
        return Err(NurbsError::Domain {
            what: format!(
                "closest-{kind} request can retain {retained_bytes} bytes above defensive ceiling {CLOSEST_MAX_RETAINED_BYTES}"
            ),
        });
    }
    Ok(())
}

fn subdivision_frontier_bytes(
    seed_leaves: u128,
    payload_per_leaf: u128,
    scratch_leaves: u128,
    work_per_split: u128,
    max_splits: u32,
    kind: &str,
) -> Result<u128, NurbsError> {
    let split_count = u128::from(max_splits);
    let retained_leaves = seed_leaves
        .checked_add(split_count)
        .and_then(|leaves| leaves.checked_add(if max_splits == 0 { 0 } else { scratch_leaves }))
        .ok_or_else(|| NurbsError::Domain {
            what: format!("closest-{kind} retained-leaf accounting overflows u128"),
        })?;
    let retained_bytes = retained_leaves
        .checked_mul(payload_per_leaf)
        .ok_or_else(|| NurbsError::Domain {
            what: format!("closest-{kind} retained-byte accounting overflows u128"),
        })?;
    let split_work = split_count
        .checked_mul(work_per_split)
        .ok_or_else(|| NurbsError::Domain {
            what: format!("closest-{kind} split-work accounting overflows u128"),
        })?;
    if split_work > CLOSEST_MAX_SPLIT_WORK_UNITS {
        return Err(NurbsError::Domain {
            what: format!(
                "closest-{kind} request needs {split_work} split-work units above defensive ceiling {CLOSEST_MAX_SPLIT_WORK_UNITS}"
            ),
        });
    }
    Ok(retained_bytes)
}

fn enforce_curve_retained_envelope(
    conversion_peak_bytes: u128,
    persistent_curve_bytes: u128,
    frontier_bytes: u128,
    polish_peak_bytes: u128,
) -> Result<(), NurbsError> {
    let traversal_peak = persistent_curve_bytes
        .checked_add(frontier_bytes)
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-curve aggregate retained-byte accounting overflows u128".to_string(),
        })?;
    enforce_retained_bytes(
        conversion_peak_bytes
            .max(traversal_peak)
            .max(polish_peak_bytes),
        "curve",
    )
}

#[derive(Debug, Clone, Copy)]
struct CurveClosestPlan {
    conversion: BezierConversionPlan,
    seed_leaves: usize,
    queue_capacity: usize,
}

fn binary_heap_height(capacity: usize) -> u128 {
    if capacity <= 1 {
        0
    } else {
        u128::from(usize::BITS - (capacity - 1).leading_zeros())
    }
}

fn curve_subdivision_work_per_split(order: u128, heap_height: u128) -> Result<u128, NurbsError> {
    let triangular = order
        .checked_mul(order.saturating_sub(1))
        .map(|product| product / 2)
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-curve split triangle overflows u128".to_string(),
        })?;
    let geometric_work = triangular
        .checked_mul(4)
        .and_then(|work| {
            order
                .checked_mul(16)
                .and_then(|linear| work.checked_add(linear))
        })
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-curve split work overflows u128".to_string(),
        })?;
    let heap_work = checked_product(
        &[heap_height, 3, 16],
        "closest-curve pop and child-push heap work",
    )?;
    checked_sum(
        &[geometric_work, heap_work],
        "closest-curve aggregate split work",
    )
}

fn preflight_curve_closest<const DIM: usize>(
    curve: AdmittedNurbsCurve<'_, f64, DIM>,
    max_splits: u32,
    source_admission_work: u128,
) -> Result<CurveClosestPlan, NurbsError> {
    // This count-only gate precedes the first knot-run scan performed by the
    // exact conversion planner.
    let pre_scan_work = curve.bezier_pre_scan_work()?;
    let work_before_plan_scan = source_admission_work
        .checked_add(pre_scan_work)
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-curve admission and plan-scan work overflows u128".to_string(),
        })?;
    enforce_base_work(work_before_plan_scan, "curve plan scan")?;
    let conversion = curve.bezier_conversion_plan()?;
    let knots = curve.knots();
    let order = (knots.degree() as u128)
        .checked_add(1)
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-curve order accounting overflows u128".to_string(),
        })?;
    let seed_leaves = knots
        .knots()
        .windows(2)
        .filter(|pair| pair[1] > pair[0])
        .count();
    if seed_leaves == 0 {
        return Err(NurbsError::Structure {
            what: "closest-curve admitted source has no nonempty knot span".to_string(),
        });
    }
    let split_capacity = usize::try_from(max_splits).map_err(|_| NurbsError::Domain {
        what: "closest-curve split count is not representable as usize".to_string(),
    })?;
    let queue_capacity =
        seed_leaves
            .checked_add(split_capacity)
            .ok_or_else(|| NurbsError::Domain {
                what: "closest-curve queue capacity overflows usize".to_string(),
            })?;
    let heap_height = binary_heap_height(queue_capacity);
    let seed_heap_work = checked_product(
        &[seed_leaves as u128, heap_height, 16],
        "closest-curve seed heap work",
    )?;
    let (derivative_work, derivative_bytes) = NurbsCurve::<f64, DIM>::derivative_envelope(
        curve.homogeneous_control_points().len(),
        knots.knots().len(),
        knots.degree(),
        2,
    )?;
    let final_eval_work = checked_sum(
        &[
            knots.knots().len() as u128,
            checked_product(&[order, order, 8], "curve polish evaluation")?,
        ],
        "curve final polish evaluation",
    )?;
    let polish_work = derivative_work
        .checked_mul(12)
        .and_then(|work| work.checked_add(final_eval_work))
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-curve polish-work accounting overflows u128".to_string(),
        })?;
    enforce_base_work(
        curve_base_work_units(
            curve,
            source_admission_work,
            conversion,
            seed_leaves as u128,
            order,
            polish_work,
            seed_heap_work,
        )?,
        "curve",
    )?;
    let work_per_split = curve_subdivision_work_per_split(order, heap_height)?;
    let payload_per_leaf = order
        .checked_mul(size_of::<[f64; 4]>() as u128)
        .and_then(|payload| payload.checked_add(size_of::<MinEntry<Seg>>() as u128))
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-curve leaf payload overflows u128".to_string(),
        })?;
    let frontier_bytes = subdivision_frontier_bytes(
        seed_leaves as u128,
        payload_per_leaf,
        3,
        work_per_split,
        max_splits,
        "curve",
    )?;
    let source_knot_bytes = (knots.knots().len() as u128)
        .checked_mul(size_of::<f64>() as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-curve source-knot byte accounting overflows u128".to_string(),
        })?;
    let source_control_bytes = (curve.homogeneous_control_points().len() as u128)
        .checked_mul(size_of::<[f64; 4]>() as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-curve source-control byte accounting overflows u128".to_string(),
        })?;
    let source_bytes = source_knot_bytes
        .checked_add(source_control_bytes)
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-curve source byte accounting overflows u128".to_string(),
        })?;
    let conversion_peak_bytes = source_bytes
        .checked_add(conversion.peak_allocated_bytes)
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-curve conversion peak accounting overflows u128".to_string(),
        })?;
    let persistent_curve_bytes = source_bytes
        .checked_add(conversion.converted_bytes)
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-curve persistent curve accounting overflows u128".to_string(),
        })?;
    let final_eval_bytes = checked_product(
        &[order, 3, size_of::<f64>() as u128],
        "closest-curve final-evaluation workspace",
    )?;
    let polish_peak_bytes = source_bytes
        .checked_add(derivative_bytes.max(final_eval_bytes))
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-curve polish peak accounting overflows u128".to_string(),
        })?;
    enforce_curve_retained_envelope(
        conversion_peak_bytes,
        persistent_curve_bytes,
        frontier_bytes,
        polish_peak_bytes,
    )?;
    Ok(CurveClosestPlan {
        conversion,
        seed_leaves,
        queue_capacity,
    })
}

pub(crate) fn surface_subdivision_work_per_split(
    surface: &NurbsSurface<f64>,
) -> Result<u128, NurbsError> {
    let order_u = (surface.knots_u.degree as u128)
        .checked_add(1)
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-surface u-order accounting overflows u128".to_string(),
        })?;
    let order_v = (surface.knots_v.degree as u128)
        .checked_add(1)
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-surface v-order accounting overflows u128".to_string(),
        })?;
    let seed_capacity = (surface.knots_u.control_count() as u128)
        .checked_mul(surface.knots_v.control_count() as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-surface seed-capacity accounting overflows u128".to_string(),
        })?;
    let (_, _, worst_heap_height) = surface_queue_shape(seed_capacity, CLOSEST_MAX_SPLITS)?;
    checked_sum(
        &[
            surface_geometric_work_per_split(order_u, order_v)?,
            surface_split_heap_work(worst_heap_height)?,
        ],
        "closest-surface worst-case split work",
    )
}

fn surface_seed_leaf_count(surface: &NurbsSurface<f64>) -> Result<u128, NurbsError> {
    let spans_u = surface
        .knots_u
        .knots
        .windows(2)
        .filter(|pair| pair[1] > pair[0])
        .count() as u128;
    let spans_v = surface
        .knots_v
        .knots
        .windows(2)
        .filter(|pair| pair[1] > pair[0])
        .count() as u128;
    spans_u
        .checked_mul(spans_v)
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-surface seed-leaf accounting overflows u128".to_string(),
        })
}

fn surface_geometric_work_per_split(order_u: u128, order_v: u128) -> Result<u128, NurbsError> {
    let controls = order_u
        .checked_mul(order_v)
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-surface patch payload overflows u128".to_string(),
        })?;
    let triangle_u = order_u
        .checked_mul(order_u.saturating_sub(1))
        .map(|product| product / 2)
        .and_then(|triangle| triangle.checked_mul(order_v))
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-surface u-split triangle overflows u128".to_string(),
        })?;
    let triangle_v = order_v
        .checked_mul(order_v.saturating_sub(1))
        .map(|product| product / 2)
        .and_then(|triangle| triangle.checked_mul(order_u))
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-surface v-split triangle overflows u128".to_string(),
        })?;
    triangle_u
        .max(triangle_v)
        .checked_mul(4)
        .and_then(|work| {
            controls
                .checked_mul(16)
                .and_then(|linear| work.checked_add(linear))
        })
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-surface split work overflows u128".to_string(),
        })
}

#[derive(Debug, Clone, Copy)]
struct SurfaceClosestPlan {
    seed_leaves: usize,
    queue_capacity: usize,
    seed_heap_work: u128,
}

fn surface_queue_shape(
    seed_leaves: u128,
    max_splits: u32,
) -> Result<(usize, usize, u128), NurbsError> {
    let seed_leaves = usize::try_from(seed_leaves).map_err(|_| NurbsError::Domain {
        what: "closest-surface seed-leaf count is not representable as usize".to_string(),
    })?;
    let split_capacity = usize::try_from(max_splits).map_err(|_| NurbsError::Domain {
        what: "closest-surface split count is not representable as usize".to_string(),
    })?;
    let queue_capacity =
        seed_leaves
            .checked_add(split_capacity)
            .ok_or_else(|| NurbsError::Domain {
                what: "closest-surface queue capacity overflows usize".to_string(),
            })?;
    Ok((
        seed_leaves,
        queue_capacity,
        binary_heap_height(queue_capacity),
    ))
}

fn surface_seed_heap_work(seed_leaves: u128, heap_height: u128) -> Result<u128, NurbsError> {
    checked_product(
        &[seed_leaves, heap_height, 16],
        "closest-surface seed heap work",
    )
}

fn surface_split_heap_work(heap_height: u128) -> Result<u128, NurbsError> {
    checked_product(
        &[heap_height, 3, 16],
        "closest-surface pop and child-push heap work",
    )
}

fn surface_closest_plan(
    order_u: u128,
    order_v: u128,
    seed_leaves: u128,
    max_splits: u32,
) -> Result<SurfaceClosestPlan, NurbsError> {
    let (seed_leaves, queue_capacity, heap_height) = surface_queue_shape(seed_leaves, max_splits)?;
    let geometric_work = surface_geometric_work_per_split(order_u, order_v)?;
    let heap_work = surface_split_heap_work(heap_height)?;
    let work_per_split = checked_sum(
        &[geometric_work, heap_work],
        "closest-surface aggregate split work",
    )?;
    let controls = order_u
        .checked_mul(order_v)
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-surface patch payload overflows u128".to_string(),
        })?;
    let payload_per_leaf = controls
        .checked_mul(size_of::<[f64; 4]>() as u128)
        .and_then(|payload| {
            order_u
                .checked_mul(size_of::<Vec<[f64; 4]>>() as u128)
                .and_then(|headers| payload.checked_add(headers))
        })
        .and_then(|payload| payload.checked_add(size_of::<MinEntry<Patch>>() as u128))
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-surface leaf payload overflows u128".to_string(),
        })?;
    let frontier_bytes = subdivision_frontier_bytes(
        seed_leaves as u128,
        payload_per_leaf,
        3,
        work_per_split,
        max_splits,
        "surface",
    )?;
    enforce_retained_bytes(frontier_bytes, "surface")?;
    let seed_heap_work = surface_seed_heap_work(seed_leaves as u128, heap_height)?;
    Ok(SurfaceClosestPlan {
        seed_leaves,
        queue_capacity,
        seed_heap_work,
    })
}

pub(crate) fn preflight_surface_subdivision(
    surface: &NurbsSurface<f64>,
    max_splits: u32,
) -> Result<(), NurbsError> {
    let order_u = (surface.knots_u.degree as u128)
        .checked_add(1)
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-surface u-order accounting overflows u128".to_string(),
        })?;
    let order_v = (surface.knots_v.degree as u128)
        .checked_add(1)
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-surface v-order accounting overflows u128".to_string(),
        })?;
    let seed_leaves = surface_seed_leaf_count(surface)?;
    surface_closest_plan(order_u, order_v, seed_leaves, max_splits).map(|_| ())
}

fn push_heap_within_admitted_capacity<T>(
    heap: &mut BinaryHeap<MinEntry<T>>,
    entry: MinEntry<T>,
    admitted_capacity: usize,
    stage: &str,
) -> Result<(), NurbsError> {
    if heap.len() >= admitted_capacity {
        return Err(NurbsError::Domain {
            what: format!("{stage} would exceed admitted queue capacity {admitted_capacity}"),
        });
    }
    heap.push(entry);
    Ok(())
}

fn dist3(a: [f64; 3], b: [f64; 3]) -> f64 {
    norm3([a[0] - b[0], a[1] - b[1], a[2] - b[2]])
}

fn cartesian(h: &[f64; 4]) -> [f64; 3] {
    [h[0] / h[3], h[1] / h[3], h[2] / h[3]]
}

fn validate_closest_request(q: [f64; 3], tol: f64, max_splits: u32) -> Result<(), NurbsError> {
    if q.iter().any(|coordinate| !coordinate.is_finite()) {
        return Err(NurbsError::Domain {
            what: "closest-point query coordinates must be finite".to_string(),
        });
    }
    if !tol.is_finite() || tol < 0.0 {
        return Err(NurbsError::Domain {
            what: "closest-point tolerance must be finite and non-negative".to_string(),
        });
    }
    if max_splits > CLOSEST_MAX_SPLITS {
        return Err(NurbsError::Domain {
            what: format!(
                "closest-point split request {max_splits} exceeds defensive ceiling {CLOSEST_MAX_SPLITS}"
            ),
        });
    }
    Ok(())
}

fn hull_lower_bound<'a>(q: [f64; 3], cps: impl Iterator<Item = &'a [f64; 4]>) -> f64 {
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    for h in cps {
        let c = cartesian(h);
        for k in 0..3 {
            min[k] = min[k].min(c[k]);
            max[k] = max[k].max(c[k]);
        }
    }
    let mut gaps = [0.0; 3];
    for k in 0..3 {
        // One-ULP expansion removes the old absolute-coordinate epsilon whose
        // width grew catastrophically under translation. This remains measured
        // rather than certified because upstream homogeneous arithmetic is not
        // interval-tracked.
        min[k] = next_down(min[k]);
        max[k] = next_up(max[k]);
        gaps[k] = if q[k] < min[k] {
            min[k] - q[k]
        } else if q[k] > max[k] {
            q[k] - max[k]
        } else {
            0.0
        };
    }
    norm3(gaps)
}

/// De Casteljau split of a homogeneous Bézier control net at 1/2.
fn split_bezier(cps: &[[f64; 4]]) -> Result<(Vec<[f64; 4]>, Vec<[f64; 4]>), NurbsError> {
    let n = cps.len();
    let mut tri = Vec::new();
    tri.try_reserve_exact(n).map_err(|_| NurbsError::Domain {
        what: "Bezier split triangle allocation was refused".to_string(),
    })?;
    tri.extend_from_slice(cps);
    let mut left = Vec::new();
    left.try_reserve_exact(n).map_err(|_| NurbsError::Domain {
        what: "Bezier split left-boundary allocation was refused".to_string(),
    })?;
    let mut right = Vec::new();
    right.try_reserve_exact(n).map_err(|_| NurbsError::Domain {
        what: "Bezier split right-boundary allocation was refused".to_string(),
    })?;
    right.resize(n, [0.0f64; 4]);
    left.push(tri[0]);
    right[n - 1] = tri[n - 1];
    for level in 1..n {
        for i in 0..n - level {
            let (head, tail) = tri.split_at_mut(i + 1);
            for (x, &y) in head[i].iter_mut().zip(&tail[0]) {
                *x = f64::midpoint(*x, y);
            }
        }
        left.push(tri[0]);
        right[n - 1 - level] = tri[n - 1 - level];
    }
    Ok((left, right))
}

struct Seg {
    cpw: Vec<[f64; 4]>,
    t0: f64,
    t1: f64,
}

/// Measured closest-point bracket estimate on a curve (best-first B&B +
/// Newton polish).
///
/// # Errors
/// Propagates evaluation/domain errors.
pub fn closest_point_curve(
    curve: &NurbsCurve<f64, 3>,
    q: [f64; 3],
    tol: f64,
    max_splits: u32,
) -> Result<DistanceBracketEstimate, NurbsError> {
    validate_closest_request(q, tol, max_splits)?;
    let source_admission_work = curve_source_admission_work(curve)?;
    enforce_base_work(source_admission_work, "curve source admission")?;
    closest_point_curve_after_request(curve.admit()?, q, tol, max_splits, source_admission_work)
}

impl AdmittedNurbsCurve<'_, f64, 3> {
    /// Measured closest-point estimate while reusing this admitted immutable
    /// source snapshot across planning, exact conversion, and Newton polish.
    ///
    /// # Errors
    /// Returns a structured refusal for an invalid request, excessive
    /// work/retained payload, allocation failure, or non-finite arithmetic.
    pub fn closest_point(
        &self,
        q: [f64; 3],
        tol: f64,
        max_splits: u32,
    ) -> Result<DistanceBracketEstimate, NurbsError> {
        validate_closest_request(q, tol, max_splits)?;
        closest_point_curve_after_request(*self, q, tol, max_splits, 0)
    }
}

fn closest_point_curve_after_request(
    curve: AdmittedNurbsCurve<'_, f64, 3>,
    q: [f64; 3],
    tol: f64,
    max_splits: u32,
    source_admission_work: u128,
) -> Result<DistanceBracketEstimate, NurbsError> {
    let plan = preflight_curve_closest(curve, max_splits, source_admission_work)?;
    let bezier_curve = curve.to_bezier_form()?;
    let bezier = bezier_curve.admitted_after_validation();
    let knots = bezier.knots();
    let entries = knots.knots();
    let controls = bezier.homogeneous_control_points();
    if entries.len() != plan.conversion.final_knot_count
        || controls.len() != plan.conversion.final_control_count
    {
        return Err(NurbsError::Structure {
            what: "closest-curve conversion did not match its admitted plan".to_string(),
        });
    }
    let p = knots.degree();
    let mut queue: BinaryHeap<MinEntry<Seg>> = BinaryHeap::new();
    queue
        .try_reserve_exact(plan.queue_capacity)
        .map_err(|_| NurbsError::Domain {
            what: format!(
                "closest-curve queue allocation was refused for admitted capacity {}",
                plan.queue_capacity
            ),
        })?;
    let mut next_logical_id = 0u64;
    let mut upper = f64::INFINITY;
    let mut best_t = knots.domain().0;
    for span in p..knots.control_count() {
        let (t0, t1) = (entries[span], entries[span + 1]);
        if t1 <= t0 {
            continue;
        }
        let mut cpw = Vec::new();
        cpw.try_reserve_exact(p + 1)
            .map_err(|_| NurbsError::Domain {
                what: "closest-curve seed control allocation was refused".to_string(),
            })?;
        cpw.extend_from_slice(&controls[span - p..=span]);
        for (h, tt) in [(&cpw[0], t0), (&cpw[p], t1)] {
            let d = dist3(cartesian(h), q);
            if d < upper {
                upper = d;
                best_t = tt;
            }
        }
        let lb = hull_lower_bound(q, cpw.iter());
        push_heap_within_admitted_capacity(
            &mut queue,
            MinEntry {
                key: lb,
                logical_id: next_logical_id,
                value: Seg { cpw, t0, t1 },
            },
            plan.queue_capacity,
            "closest-curve seed",
        )?;
        next_logical_id += 1;
    }
    if queue.len() != plan.seed_leaves {
        return Err(NurbsError::Structure {
            what: format!(
                "closest-curve conversion produced {} seed leaves after admitting {}",
                queue.len(),
                plan.seed_leaves
            ),
        });
    }
    if queue.is_empty() || !upper.is_finite() || queue.iter().any(|entry| !entry.key.is_finite()) {
        return Err(NurbsError::Domain {
            what: "closest-curve initial distance bounds are not finite".to_string(),
        });
    }
    let mut iterations = 0u32;
    while iterations < max_splits {
        let Some(entry) = queue.peek() else {
            break;
        };
        if upper - entry.key <= tol {
            break; // bracket closed
        }
        let Some(entry) = queue.pop() else {
            break;
        };
        let seg = entry.value;
        let tm = f64::midpoint(seg.t0, seg.t1);
        if tm == seg.t0 || tm == seg.t1 {
            // De Casteljau would still create a geometric half-split even
            // though the recorded parameter interval cannot shrink. Retain
            // the leaf so its lower estimate remains part of the bracket.
            push_heap_within_admitted_capacity(
                &mut queue,
                MinEntry {
                    key: entry.key,
                    logical_id: entry.logical_id,
                    value: seg,
                },
                plan.queue_capacity,
                "closest-curve unsplittable leaf",
            )?;
            break;
        }
        let (l, r) = split_bezier(&seg.cpw)?;
        // The split junction is C(mid): a free upper-bound sample.
        let d = dist3(cartesian(&l[l.len() - 1]), q);
        if d < upper {
            upper = d;
            best_t = tm;
        }
        for (cpw, t0, t1) in [(l, seg.t0, tm), (r, tm, seg.t1)] {
            let lb = hull_lower_bound(q, cpw.iter());
            if !lb.is_finite() {
                return Err(NurbsError::Domain {
                    what: "closest-curve child distance bound is not finite".to_string(),
                });
            }
            if lb < upper {
                push_heap_within_admitted_capacity(
                    &mut queue,
                    MinEntry {
                        key: lb,
                        logical_id: next_logical_id,
                        value: Seg { cpw, t0, t1 },
                    },
                    plan.queue_capacity,
                    "closest-curve child",
                )?;
                next_logical_id += 1;
            }
        }
        iterations += 1;
    }
    let lower = queue.peek().map_or(upper, |entry| entry.key);
    // The bracket no longer needs the converted generation or frontier.
    // Release both phases before optional derivative polish so the aggregate
    // peak is max(conversion, search, polish), not their sum.
    drop(queue);
    drop(bezier_curve);
    // Newton polish on g(t) = (C − q)·C' sharpens the upper bound.
    let (dlo, dhi) = curve.knots().domain();
    let mut t = best_t;
    for _ in 0..12 {
        // Polishing is an optional improvement to an already valid measured
        // B&B estimate. A derivative can be undefined at a repeated knot (or
        // unavailable in the legacy jet API); that must not erase the retained
        // geometric witness and bracket.
        let Ok(ders) = curve.derivatives(t, 2) else {
            break;
        };
        if ders.len() < 2 {
            break;
        }
        let second = ders.get(2).copied().unwrap_or([0.0; 3]);
        let diff = [ders[0][0] - q[0], ders[0][1] - q[1], ders[0][2] - q[2]];
        let g: f64 = (0..3).map(|k| diff[k] * ders[1][k]).sum();
        let gp: f64 = (0..3)
            .map(|k| ders[1][k] * ders[1][k] + diff[k] * second[k])
            .sum();
        if !g.is_finite() || !gp.is_finite() || gp.abs() < 1e-300 {
            break;
        }
        let next = (t - g / gp).clamp(dlo, dhi);
        if !next.is_finite() || next == t {
            break;
        }
        t = next;
    }
    let (upper, best_t) = if let Ok(point) = curve.eval(t) {
        let polished = dist3(point, q);
        if polished.is_finite() && polished < upper {
            (polished, t)
        } else {
            (upper, best_t)
        }
    } else {
        // Newton is optional polish. Its failure cannot erase the finite
        // branch-and-bound witness retained above.
        (upper, best_t)
    };
    Ok(DistanceBracketEstimate {
        lower: lower.min(upper),
        upper,
        param: [best_t, 0.0],
        iterations,
    })
}

/// A homogeneous Bézier control net (rows × cols).
type Net = Vec<Vec<[f64; 4]>>;

struct Patch {
    cpw: Net, // (pu+1) rows × (pv+1) cols
    u0: f64,
    u1: f64,
    v0: f64,
    v1: f64,
    depth_u: u32,
    depth_v: u32,
}

fn patch_lb(q: [f64; 3], net: &[Vec<[f64; 4]>]) -> f64 {
    hull_lower_bound(q, net.iter().flatten())
}

fn zero_net(rows: usize, cols: usize, stage: &str) -> Result<Net, NurbsError> {
    let mut net = Vec::new();
    net.try_reserve_exact(rows)
        .map_err(|_| NurbsError::Domain {
            what: format!("{stage} row-table allocation was refused"),
        })?;
    for _ in 0..rows {
        let mut row = Vec::new();
        row.try_reserve_exact(cols)
            .map_err(|_| NurbsError::Domain {
                what: format!("{stage} row allocation was refused"),
            })?;
        row.resize(cols, [0.0; 4]);
        net.push(row);
    }
    Ok(net)
}

fn split_patch_u(net: &[Vec<[f64; 4]>]) -> Result<(Net, Net), NurbsError> {
    // Split every v-column along u (rows are u direction).
    let rows = net.len();
    let cols = net[0].len();
    let mut left = zero_net(rows, cols, "u-split left patch")?;
    let mut right = zero_net(rows, cols, "u-split right patch")?;
    for j in 0..cols {
        let mut col = Vec::new();
        col.try_reserve_exact(rows)
            .map_err(|_| NurbsError::Domain {
                what: "u-split column allocation was refused".to_string(),
            })?;
        col.extend((0..rows).map(|i| net[i][j]));
        let (l, r) = split_bezier(&col)?;
        for i in 0..rows {
            left[i][j] = l[i];
            right[i][j] = r[i];
        }
    }
    Ok((left, right))
}

fn split_patch_v(net: &[Vec<[f64; 4]>]) -> Result<(Net, Net), NurbsError> {
    let (mut left, mut right) = (Vec::new(), Vec::new());
    left.try_reserve_exact(net.len())
        .map_err(|_| NurbsError::Domain {
            what: "v-split left row-table allocation was refused".to_string(),
        })?;
    right
        .try_reserve_exact(net.len())
        .map_err(|_| NurbsError::Domain {
            what: "v-split right row-table allocation was refused".to_string(),
        })?;
    for row in net {
        let (l, r) = split_bezier(row)?;
        left.push(l);
        right.push(r);
    }
    Ok((left, right))
}

/// Decompose a surface to Bézier patches via repeated knot insertion.
fn to_bezier_surface(
    surface: AdmittedNurbsSurface<'_, f64>,
) -> Result<NurbsSurface<f64>, NurbsError> {
    let mut work = surface.source().try_clone()?;
    loop {
        let mut inserted = false;
        for dir_u in [true, false] {
            let (kv, p) = if dir_u {
                (&work.knots_u, work.knots_u.degree)
            } else {
                (&work.knots_v, work.knots_v.degree)
            };
            let (lo, hi) = kv.domain()?;
            let mut target = None;
            let mut run_start = 0usize;
            while run_start < kv.knots.len() {
                let t = kv.knots[run_start];
                let mut run_end = run_start + 1;
                while run_end < kv.knots.len() && kv.knots[run_end] == t {
                    run_end += 1;
                }
                if t > lo && t < hi && run_end - run_start < p {
                    target = Some(t);
                    break;
                }
                run_start = run_end;
            }
            if let Some(t) = target {
                work = if dir_u {
                    work.insert_knot_u(t)?
                } else {
                    work.insert_knot_v(t)?
                };
                inserted = true;
            }
        }
        if !inserted {
            return Ok(work);
        }
    }
}

/// Seed the patch queue from a Bézier-form surface.
fn seed_patches(
    work: &NurbsSurface<f64>,
    q: [f64; 3],
    queue: &mut BinaryHeap<MinEntry<Patch>>,
    queue_capacity: usize,
    next_logical_id: &mut u64,
    upper: &mut f64,
    best: &mut [f64; 2],
) -> Result<(), NurbsError> {
    let (pu, pv) = (work.knots_u.degree, work.knots_v.degree);
    for su in pu..work.knots_u.control_count() {
        let (u0, u1) = (work.knots_u.knots[su], work.knots_u.knots[su + 1]);
        if u1 <= u0 {
            continue;
        }
        for sv in pv..work.knots_v.control_count() {
            let (v0, v1) = (work.knots_v.knots[sv], work.knots_v.knots[sv + 1]);
            if v1 <= v0 {
                continue;
            }
            let mut net = zero_net(pu + 1, pv + 1, "closest-surface seed patch")?;
            for (target_row, source_row) in net.iter_mut().zip(work.cpw[su - pu..=su].iter()) {
                target_row.copy_from_slice(&source_row[sv - pv..=sv]);
            }
            let d = dist3(cartesian(&net[0][0]), q);
            if d < *upper {
                *upper = d;
                *best = [u0, v0];
            }
            let lb = patch_lb(q, &net);
            push_heap_within_admitted_capacity(
                queue,
                MinEntry {
                    key: lb,
                    logical_id: *next_logical_id,
                    value: Patch {
                        cpw: net,
                        u0,
                        u1,
                        v0,
                        v1,
                        depth_u: 0,
                        depth_v: 0,
                    },
                },
                queue_capacity,
                "closest-surface seed",
            )?;
            *next_logical_id += 1;
        }
    }
    Ok(())
}

/// Measured closest point on a surface (best-first B&B over Bézier
/// patches with de Casteljau splits that balance normalized subdivision depth).
///
/// # Errors
/// Propagates evaluation/domain errors.
pub fn closest_point_surface(
    surface: &NurbsSurface<f64>,
    q: [f64; 3],
    tol: f64,
    max_splits: u32,
) -> Result<DistanceBracketEstimate, NurbsError> {
    validate_closest_request(q, tol, max_splits)?;
    let source_admission_work = surface_source_admission_work(surface)?;
    enforce_base_work(source_admission_work, "surface source admission")?;
    closest_point_surface_after_request(surface.admit()?, q, tol, max_splits, source_admission_work)
}

impl AdmittedNurbsSurface<'_, f64> {
    /// Measured closest-point estimate while reusing this admitted immutable
    /// source through conversion planning and final evaluation.
    ///
    /// # Errors
    /// Returns a structured refusal for an invalid request, excessive
    /// work/frontier payload, allocation failure, or non-finite arithmetic.
    pub fn closest_point(
        &self,
        q: [f64; 3],
        tol: f64,
        max_splits: u32,
    ) -> Result<DistanceBracketEstimate, NurbsError> {
        validate_closest_request(q, tol, max_splits)?;
        closest_point_surface_after_request(*self, q, tol, max_splits, 0)
    }
}

fn closest_point_surface_after_request(
    surface: AdmittedNurbsSurface<'_, f64>,
    q: [f64; 3],
    tol: f64,
    max_splits: u32,
    source_admission_work: u128,
) -> Result<DistanceBracketEstimate, NurbsError> {
    let base_plan = surface_base_plan_from_admitted(surface, source_admission_work)?;
    enforce_base_work(base_plan.work_units, "surface")?;
    let plan = surface_closest_plan(
        base_plan.order_u,
        base_plan.order_v,
        base_plan.seed_leaves,
        max_splits,
    )?;
    enforce_base_work(
        base_plan
            .work_units
            .checked_add(plan.seed_heap_work)
            .ok_or_else(|| NurbsError::Domain {
                what: "closest-surface base and seed-heap work overflows u128".to_string(),
            })?,
        "surface",
    )?;
    let work = to_bezier_surface(surface)?;
    let mut queue: BinaryHeap<MinEntry<Patch>> = BinaryHeap::new();
    queue
        .try_reserve_exact(plan.queue_capacity)
        .map_err(|_| NurbsError::Domain {
            what: format!(
                "closest-surface queue allocation was refused for admitted capacity {}",
                plan.queue_capacity
            ),
        })?;
    let mut next_logical_id = 0u64;
    let mut upper = f64::INFINITY;
    let mut best = [surface.knots_u().domain().0, surface.knots_v().domain().0];
    seed_patches(
        &work,
        q,
        &mut queue,
        plan.queue_capacity,
        &mut next_logical_id,
        &mut upper,
        &mut best,
    )?;
    if queue.len() != plan.seed_leaves {
        return Err(NurbsError::Structure {
            what: format!(
                "closest-surface conversion produced {} seed leaves after admitting {}",
                queue.len(),
                plan.seed_leaves
            ),
        });
    }
    if queue.is_empty() || !upper.is_finite() || queue.iter().any(|entry| !entry.key.is_finite()) {
        return Err(NurbsError::Domain {
            what: "closest-surface initial distance bounds are not finite".to_string(),
        });
    }
    let mut iterations = 0u32;
    while iterations < max_splits {
        let Some(entry) = queue.peek() else {
            break;
        };
        if upper - entry.key <= tol {
            break;
        }
        let Some(entry) = queue.pop() else {
            break;
        };
        let patch = entry.value;
        let midpoint_u = f64::midpoint(patch.u0, patch.u1);
        let midpoint_v = f64::midpoint(patch.v0, patch.v1);
        let can_split_u = midpoint_u != patch.u0 && midpoint_u != patch.u1;
        let can_split_v = midpoint_v != patch.v0 && midpoint_v != patch.v1;
        let preferred_u = patch.depth_u <= patch.depth_v;
        let split_u = match (can_split_u, can_split_v) {
            (true, true) => preferred_u,
            (true, false) => true,
            (false, true) => false,
            (false, false) => {
                push_heap_within_admitted_capacity(
                    &mut queue,
                    MinEntry {
                        key: entry.key,
                        logical_id: entry.logical_id,
                        value: patch,
                    },
                    plan.queue_capacity,
                    "closest-surface unsplittable leaf",
                )?;
                break;
            }
        };
        let midpoint = if split_u { midpoint_u } else { midpoint_v };
        if !midpoint.is_finite() {
            push_heap_within_admitted_capacity(
                &mut queue,
                MinEntry {
                    key: entry.key,
                    logical_id: entry.logical_id,
                    value: patch,
                },
                plan.queue_capacity,
                "closest-surface non-finite midpoint leaf",
            )?;
            break;
        }
        let (l, r) = if split_u {
            split_patch_u(&patch.cpw)?
        } else {
            split_patch_v(&patch.cpw)?
        };
        let halves = if split_u {
            [
                (l, patch.u0, midpoint, patch.v0, patch.v1),
                (r, midpoint, patch.u1, patch.v0, patch.v1),
            ]
        } else {
            [
                (l, patch.u0, patch.u1, patch.v0, midpoint),
                (r, patch.u0, patch.u1, midpoint, patch.v1),
            ]
        };
        for (net, u0, u1, v0, v1) in halves {
            // Corner sample improves the upper bound cheaply.
            let d = dist3(cartesian(&net[0][0]), q);
            if d < upper {
                upper = d;
                best = [u0, v0];
            }
            let lb = patch_lb(q, &net);
            if !lb.is_finite() {
                return Err(NurbsError::Domain {
                    what: "closest-surface child distance bound is not finite".to_string(),
                });
            }
            if lb < upper {
                push_heap_within_admitted_capacity(
                    &mut queue,
                    MinEntry {
                        key: lb,
                        logical_id: next_logical_id,
                        value: Patch {
                            cpw: net,
                            u0,
                            u1,
                            v0,
                            v1,
                            depth_u: patch.depth_u + u32::from(split_u),
                            depth_v: patch.depth_v + u32::from(!split_u),
                        },
                    },
                    plan.queue_capacity,
                    "closest-surface child",
                )?;
                next_logical_id += 1;
            }
        }
        iterations += 1;
    }
    let lower = queue.peek().map_or(upper, |entry| entry.key);
    // Sample the current best-lower-estimate patch center for a final
    // evaluated-point improvement. The former midpoint(best,best) expression
    // merely re-evaluated the already retained corner and never sampled a
    // patch center; it also failed to update `param` when the sample improved.
    if let Some(entry) = queue.peek() {
        let patch = &entry.value;
        let candidate = [
            f64::midpoint(patch.u0, patch.u1),
            f64::midpoint(patch.v0, patch.v1),
        ];
        if let Ok(point) = surface.eval(candidate[0], candidate[1]) {
            let distance = dist3(point, q);
            if distance.is_finite() && distance < upper {
                upper = distance;
                best = candidate;
            }
        }
    }
    Ok(DistanceBracketEstimate {
        lower: lower.min(upper),
        upper,
        param: best,
        iterations,
    })
}

#[cfg(test)]
mod tests {
    use super::{
        BinaryHeap, CLOSEST_MAX_BASE_WORK_UNITS, CLOSEST_MAX_RETAINED_BYTES, MinEntry,
        closest_point_curve, closest_point_surface, curve_subdivision_work_per_split,
        enforce_curve_retained_envelope, preflight_surface_subdivision, subdivision_frontier_bytes,
        surface_base_work_units,
    };
    use crate::{KnotVector, NurbsCurve, NurbsSurface};

    #[test]
    fn min_heap_order_is_key_then_logical_identity() {
        let mut heap = BinaryHeap::new();
        for (key, logical_id, value) in [(2.0, 8, 'd'), (1.0, 9, 'c'), (1.0, 3, 'a'), (1.0, 7, 'b')]
        {
            heap.push(MinEntry {
                key,
                logical_id,
                value,
            });
        }
        let popped: Vec<_> = core::iter::from_fn(|| heap.pop())
            .map(|entry| (entry.key, entry.logical_id, entry.value))
            .collect();
        assert_eq!(
            popped,
            vec![(1.0, 3, 'a'), (1.0, 7, 'b'), (1.0, 9, 'c'), (2.0, 8, 'd')]
        );
    }

    #[test]
    fn curve_closest_owning_and_admitted_paths_match_exactly() {
        let curve = NurbsCurve::<f64, 3>::new(
            KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2).expect("quadratic knots"),
            &[[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [2.0, 0.0, 0.0]],
            &[1.0; 3],
        )
        .expect("quadratic curve");
        let query = [1.0, 0.25, 0.0];
        let owning = closest_point_curve(&curve, query, 1e-9, 8).expect("owning closest");
        let admitted = curve.admit().expect("admitted curve");
        let first = admitted
            .closest_point(query, 1e-9, 8)
            .expect("admitted closest");
        let repeated = admitted
            .closest_point(query, 1e-9, 8)
            .expect("repeated admitted closest");
        assert_eq!(first, owning);
        assert_eq!(repeated, first);
    }

    #[test]
    fn surface_closest_owning_and_admitted_paths_match_exactly() {
        let knots = KnotVector::new(vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0], 2)
            .expect("two-span quadratic surface knots");
        let points: Vec<Vec<[f64; 3]>> = (0..4)
            .map(|u| {
                (0..4)
                    .map(|v| [f64::from(u) / 3.0, f64::from(v) / 3.0, 0.0])
                    .collect()
            })
            .collect();
        let surface = NurbsSurface::new(knots.clone(), knots, &points, &vec![vec![1.0; 4]; 4])
            .expect("quadratic multispan surface");
        let query = [0.25, 0.75, 0.5];
        let owning =
            closest_point_surface(&surface, query, 1e-9, 8).expect("owning surface closest");
        let admitted = surface.admit().expect("admitted surface");
        let first = admitted
            .closest_point(query, 1e-9, 8)
            .expect("admitted surface closest");
        let repeated = admitted
            .closest_point(query, 1e-9, 8)
            .expect("repeated admitted surface closest");
        assert_eq!(first, owning);
        assert_eq!(repeated, first);
    }

    #[test]
    fn curve_retained_envelope_composes_all_live_phases() {
        assert!(
            enforce_curve_retained_envelope(CLOSEST_MAX_RETAINED_BYTES, 0, 0, 0).is_ok(),
            "the exact conversion-phase cap is admissible"
        );
        assert!(
            enforce_curve_retained_envelope(0, CLOSEST_MAX_RETAINED_BYTES - 1, 1, 0).is_ok(),
            "persistent plus frontier bytes compose at the exact cap"
        );
        let search_error = enforce_curve_retained_envelope(0, CLOSEST_MAX_RETAINED_BYTES, 1, 0)
            .expect_err("one aggregate search byte above the cap must refuse");
        assert!(matches!(search_error, crate::NurbsError::Domain { .. }));
        let polish_error = enforce_curve_retained_envelope(0, 0, 0, CLOSEST_MAX_RETAINED_BYTES + 1)
            .expect_err("one polish byte above the cap must refuse");
        assert!(matches!(polish_error, crate::NurbsError::Domain { .. }));
        let overflow = enforce_curve_retained_envelope(0, u128::MAX, 1, 0)
            .expect_err("aggregate retained-byte overflow must refuse");
        assert!(matches!(overflow, crate::NurbsError::Domain { .. }));
    }

    #[test]
    fn stage_faithful_work_model_admits_ordinary_multispan_cubic_surface() {
        let knots = KnotVector::new(
            vec![0.0, 0.0, 0.0, 0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.0, 1.0, 1.0],
            3,
        )
        .expect("five-span cubic knots");
        let mut points = Vec::new();
        let mut weights = Vec::new();
        for i in 0..8 {
            let mut point_row = Vec::new();
            let mut weight_row = Vec::new();
            for j in 0..8 {
                point_row.push([f64::from(i), f64::from(j), 0.0]);
                weight_row.push(1.0);
            }
            points.push(point_row);
            weights.push(weight_row);
        }
        let surface = NurbsSurface::new(knots.clone(), knots, &points, &weights).expect("surface");
        let work = surface_base_work_units(&surface).expect("work estimate");
        assert!(
            work <= CLOSEST_MAX_BASE_WORK_UNITS,
            "a 5x5 cubic patch grid is ordinary work, not a cubic-in-total-grid refusal: {work}"
        );
        assert!(
            preflight_surface_subdivision(&surface, u32::MAX).is_err(),
            "the public surface path must price its worst-case retained frontier"
        );

        let high_order = 100_001u128;
        let split_work =
            curve_subdivision_work_per_split(high_order, 0).expect("representable split work");
        assert!(
            subdivision_frontier_bytes(1, 1, 3, split_work, 1, "curve").is_err(),
            "one quadratic de Casteljau split must not masquerade as one work unit"
        );
    }
}
