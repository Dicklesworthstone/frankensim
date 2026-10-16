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
use crate::curve::{
    AdmittedNurbsCurve, BezierConversionPlan, CurveAdmissionRun, CurveBezierRun,
    CurveDerivativesRun, CurveEvaluationRun, NurbsCurve,
};
use crate::surface::{
    AdmittedNurbsSurface, NurbsSurface, SurfaceAdmissionRun, SurfaceEvaluationRun,
};
use fs_exec::Cx;
use fs_math::{det, next_down, next_up};
use std::cmp::Ordering;
use std::collections::BinaryHeap;
use std::mem::size_of;

/// Defensive ceiling for the legacy allocation-bearing subdivision path.
/// Caller-owned work/memory budgets and certifying request-drain-finalize
/// cancellation belong to the successor API. The admitted measured path only
/// observes a borrowed `Cx`; this cap still prevents an unbounded `u32` request.
pub(crate) const CLOSEST_MAX_SPLITS: u32 = 1_048_576;

/// Defensive ceiling for conservative initial validation and stage-faithful
/// knot-insertion, expanded-grid, run-scan, and queue-seeding work before the
/// split counter starts. This closes the `max_splits=0` admission gap on large
/// public spline structures.
pub(crate) const CLOSEST_MAX_BASE_WORK_UNITS: u128 = 16_777_216;

/// Defensive aggregate split-work ceiling for the legacy APIs.
const CLOSEST_MAX_SPLIT_WORK_UNITS: u128 = 1_073_741_824;

/// Defensive retained-payload ceiling. Curve and surface admission compose
/// their borrowed source, exact-conversion allocations, search frontier, and
/// final evaluation or post-release polish phases.
const CLOSEST_MAX_RETAINED_BYTES: u128 = 256 * 1024 * 1024;

// Keep the source-snapshot charge aligned with the admitted curve/knot
// validation envelopes. This is an accounting coefficient, not a timing
// claim for any particular machine.
const CLOSEST_CURVE_SOURCE_SCAN_WORK_PER_ENTRY: u128 = 16;
const CLOSEST_SURFACE_SOURCE_SCAN_WORK_PER_CONTROL: u128 = 16;
const CLOSEST_CANCELLATION_STRIDE: usize = 64;

fn closest_poll_due(
    operations_since_poll: &mut usize,
    should_cancel: &mut impl FnMut() -> bool,
) -> bool {
    *operations_since_poll += 1;
    if *operations_since_poll < CLOSEST_CANCELLATION_STRIDE {
        return false;
    }
    *operations_since_poll = 0;
    should_cancel()
}

#[derive(Debug, Clone, PartialEq)]
enum ClosestWorkRun<T> {
    Complete(T),
    Cancelled,
}

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

/// Transactional terminal state of cancellation-aware measured closest-point
/// search.
#[must_use]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ClosestPointRun {
    /// A complete finite measured bracket estimate.
    Complete {
        /// The fully published closest-point estimate.
        estimate: DistanceBracketEstimate,
    },
    /// Cancellation was observed; no partial estimate was published.
    Cancelled,
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

fn nonempty_span_count_with_poll(
    knots: AdmittedKnotVector<'_, f64>,
    should_cancel: &mut impl FnMut() -> bool,
) -> ClosestWorkRun<u128> {
    if should_cancel() {
        return ClosestWorkRun::Cancelled;
    }
    let mut count = 0u128;
    let mut operations_since_poll = 0usize;
    for pair in knots.knots().windows(2) {
        if pair[1] > pair[0] {
            count += 1;
        }
        if closest_poll_due(&mut operations_since_poll, should_cancel) {
            return ClosestWorkRun::Cancelled;
        }
    }
    if should_cancel() {
        return ClosestWorkRun::Cancelled;
    }
    ClosestWorkRun::Complete(count)
}

fn bezier_insertion_count_with_poll(
    knots: AdmittedKnotVector<'_, f64>,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<ClosestWorkRun<u128>, NurbsError> {
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    let (lo, hi) = knots.domain();
    let entries = knots.knots();
    let degree = knots.degree();
    let mut insertions = 0u128;
    let mut run_start = 0usize;
    let mut operations_since_poll = 0usize;
    while run_start < entries.len() {
        let knot = entries[run_start];
        let mut run_end = run_start + 1;
        while run_end < entries.len() && entries[run_end] == knot {
            run_end += 1;
            if closest_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(ClosestWorkRun::Cancelled);
            }
        }
        let multiplicity = run_end - run_start;
        if knot > lo && knot < hi && multiplicity < degree {
            insertions = insertions
                .checked_add((degree - multiplicity) as u128)
                .ok_or_else(|| NurbsError::Domain {
                    what: "Bézier insertion-count accounting overflows u128".to_string(),
                })?;
        }
        if closest_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(ClosestWorkRun::Cancelled);
        }
        run_start = run_end;
    }
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    Ok(ClosestWorkRun::Complete(insertions))
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

// Requested heap payload for the two knot arrays, the outer control-row table,
// and each row's control allocation. Allocator metadata, rounding, and the
// inline `Vec` headers owned by the surface itself are intentionally excluded.
fn surface_storage_bytes(
    knot_count_u: u128,
    knot_count_v: u128,
    control_count_u: u128,
    control_count_v: u128,
) -> Result<u128, NurbsError> {
    let knot_bytes = checked_product(
        &[
            knot_count_u
                .checked_add(knot_count_v)
                .ok_or_else(|| NurbsError::Domain {
                    what: "surface knot-storage count overflows u128".to_string(),
                })?,
            size_of::<f64>() as u128,
        ],
        "surface knot storage",
    )?;
    let row_table_bytes = checked_product(
        &[control_count_u, size_of::<Vec<[f64; 4]>>() as u128],
        "surface row-table storage",
    )?;
    let control_bytes = checked_product(
        &[
            control_count_u,
            control_count_v,
            size_of::<[f64; 4]>() as u128,
        ],
        "surface control storage",
    )?;
    checked_sum(
        &[knot_bytes, row_table_bytes, control_bytes],
        "surface storage",
    )
}

// During direct tensor insertion, the borrowed source is accounted separately
// while the current derived surface and its successor overlap. The last
// insertion has the largest such pair, so its exact requested payload dominates
// every earlier generation without inventing one-dimensional curve scratch.
fn surface_conversion_peak_allocated_bytes(
    insertions_u: u128,
    insertions_v: u128,
    final_knot_count_u: u128,
    final_knot_count_v: u128,
    final_control_count_u: u128,
    final_control_count_v: u128,
) -> Result<(u128, u128), NurbsError> {
    let converted_bytes = surface_storage_bytes(
        final_knot_count_u,
        final_knot_count_v,
        final_control_count_u,
        final_control_count_v,
    )?;
    if insertions_u == 0 && insertions_v == 0 {
        return Ok((converted_bytes, converted_bytes));
    }

    // The conversion inserts U then V once per outer pass. Therefore V is the
    // final direction when its count is at least U's; otherwise U is final.
    let final_axis_is_v = insertions_v > 0 && insertions_v >= insertions_u;
    let (previous_knots_u, previous_knots_v, previous_controls_u, previous_controls_v) =
        if final_axis_is_v {
            (
                final_knot_count_u,
                final_knot_count_v
                    .checked_sub(1)
                    .ok_or_else(|| NurbsError::Domain {
                        what: "surface previous v-knot count underflows u128".to_string(),
                    })?,
                final_control_count_u,
                final_control_count_v
                    .checked_sub(1)
                    .ok_or_else(|| NurbsError::Domain {
                        what: "surface previous v-control count underflows u128".to_string(),
                    })?,
            )
        } else {
            (
                final_knot_count_u
                    .checked_sub(1)
                    .ok_or_else(|| NurbsError::Domain {
                        what: "surface previous u-knot count underflows u128".to_string(),
                    })?,
                final_knot_count_v,
                final_control_count_u
                    .checked_sub(1)
                    .ok_or_else(|| NurbsError::Domain {
                        what: "surface previous u-control count underflows u128".to_string(),
                    })?,
                final_control_count_v,
            )
        };
    let previous_bytes = surface_storage_bytes(
        previous_knots_u,
        previous_knots_v,
        previous_controls_u,
        previous_controls_v,
    )?;
    let peak_allocated_bytes =
        converted_bytes
            .checked_add(previous_bytes)
            .ok_or_else(|| NurbsError::Domain {
                what: "surface conversion peak allocated-byte accounting overflows u128"
                    .to_string(),
            })?;
    Ok((converted_bytes, peak_allocated_bytes))
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
    source_bytes: u128,
    converted_bytes: u128,
    conversion_peak_allocated_bytes: u128,
    final_knot_count_u: u128,
    final_knot_count_v: u128,
    final_control_count_u: u128,
    final_control_count_v: u128,
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
    let mut never_cancel = || false;
    match surface_base_plan_from_admitted_with_poll(
        surface,
        source_admission_work,
        &mut never_cancel,
    )? {
        ClosestWorkRun::Complete(plan) => Ok(plan),
        ClosestWorkRun::Cancelled => Err(NurbsError::Domain {
            what: "non-cancelling closest-surface planning observed cancellation".to_string(),
        }),
    }
}

fn surface_base_plan_from_admitted_with_poll(
    surface: AdmittedNurbsSurface<'_, f64>,
    source_admission_work: u128,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<ClosestWorkRun<SurfaceBasePlan>, NurbsError> {
    let knots_u = surface.knots_u();
    let knots_v = surface.knots_v();
    let expected_u = knots_u.control_count();
    let expected_v = knots_v.control_count();
    let input_controls = checked_product(
        &[
            surface.homogeneous_control_net().len() as u128,
            expected_v as u128,
        ],
        "surface input grid",
    )?;
    let work_before_plan_scans = checked_sum(
        &[
            source_admission_work,
            input_controls,
            knots_u.knots().len() as u128,
            knots_v.knots().len() as u128,
        ],
        "closest-surface pre-scan work",
    )?;
    enforce_base_work(work_before_plan_scans, "surface plan scan")?;
    let spans_u = match nonempty_span_count_with_poll(knots_u, should_cancel) {
        ClosestWorkRun::Complete(count) => count,
        ClosestWorkRun::Cancelled => return Ok(ClosestWorkRun::Cancelled),
    };
    let spans_v = match nonempty_span_count_with_poll(knots_v, should_cancel) {
        ClosestWorkRun::Complete(count) => count,
        ClosestWorkRun::Cancelled => return Ok(ClosestWorkRun::Cancelled),
    };
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
    let insertions_u = match bezier_insertion_count_with_poll(knots_u, should_cancel)? {
        ClosestWorkRun::Complete(count) => count,
        ClosestWorkRun::Cancelled => return Ok(ClosestWorkRun::Cancelled),
    };
    let insertions_v = match bezier_insertion_count_with_poll(knots_v, should_cancel)? {
        ClosestWorkRun::Complete(count) => count,
        ClosestWorkRun::Cancelled => return Ok(ClosestWorkRun::Cancelled),
    };
    let insertion_work = surface.projected_directional_insertion_work(
        usize::try_from(insertions_u).map_err(|_| NurbsError::Domain {
            what: "surface u Bézier insertion count exceeds usize".to_string(),
        })?,
        usize::try_from(insertions_v).map_err(|_| NurbsError::Domain {
            what: "surface v Bézier insertion count exceeds usize".to_string(),
        })?,
    )?;
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
    let source_bytes = surface_storage_bytes(
        knots_u.knots().len() as u128,
        knots_v.knots().len() as u128,
        expected_u as u128,
        expected_v as u128,
    )?;
    let (converted_bytes, conversion_peak_allocated_bytes) =
        surface_conversion_peak_allocated_bytes(
            insertions_u,
            insertions_v,
            expanded_knots_u,
            expanded_knots_v,
            expanded_u,
            expanded_v,
        )?;
    let expanded_grid = checked_product(&[expanded_u, expanded_v], "surface Bézier grid")?;
    let patch_controls =
        checked_product(&[seed_leaves, order_u, order_v], "surface patch seeding")?;
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
            insertion_work,
            scan_work,
        ],
        "surface base work",
    )?;
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    Ok(ClosestWorkRun::Complete(SurfaceBasePlan {
        work_units,
        seed_leaves,
        order_u,
        order_v,
        source_bytes,
        converted_bytes,
        conversion_peak_allocated_bytes,
        final_knot_count_u: expanded_knots_u,
        final_knot_count_v: expanded_knots_v,
        final_control_count_u: expanded_u,
        final_control_count_v: expanded_v,
    }))
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

fn surface_final_eval_workspace_bytes(order_u: u128, order_v: u128) -> Result<u128, NurbsError> {
    let u_basis_peak = order_u.checked_mul(3).ok_or_else(|| NurbsError::Domain {
        what: "closest-surface u-basis workspace overflows u128".to_string(),
    })?;
    let v_basis_peak = order_v
        .checked_mul(3)
        .and_then(|workspace| workspace.checked_add(order_u))
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-surface v-basis workspace overflows u128".to_string(),
        })?;
    u_basis_peak
        .max(v_basis_peak)
        .checked_mul(size_of::<f64>() as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-surface final-evaluation workspace overflows u128".to_string(),
        })
}

fn enforce_surface_retained_envelope(
    source_bytes: u128,
    conversion_peak_allocated_bytes: u128,
    converted_bytes: u128,
    frontier_bytes: u128,
    final_eval_workspace_bytes: u128,
) -> Result<(), NurbsError> {
    let conversion_peak = source_bytes
        .checked_add(conversion_peak_allocated_bytes)
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-surface conversion peak accounting overflows u128".to_string(),
        })?;
    let traversal_peak = checked_sum(
        &[source_bytes, converted_bytes, frontier_bytes],
        "closest-surface traversal retained bytes",
    )?;
    let final_eval_peak = traversal_peak
        .checked_add(final_eval_workspace_bytes)
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-surface final-evaluation peak accounting overflows u128".to_string(),
        })?;
    enforce_retained_bytes(
        conversion_peak.max(traversal_peak).max(final_eval_peak),
        "surface",
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
    let mut never_cancel = || false;
    match preflight_curve_closest_with_poll(
        curve,
        max_splits,
        source_admission_work,
        &mut never_cancel,
    )? {
        ClosestWorkRun::Complete(plan) => Ok(plan),
        ClosestWorkRun::Cancelled => Err(NurbsError::Domain {
            what: "non-cancelling closest-curve preflight observed cancellation".to_string(),
        }),
    }
}

fn preflight_curve_closest_with_poll<const DIM: usize>(
    curve: AdmittedNurbsCurve<'_, f64, DIM>,
    max_splits: u32,
    source_admission_work: u128,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<ClosestWorkRun<CurveClosestPlan>, NurbsError> {
    // This count-only gate precedes the first knot-run scan performed by the
    // exact conversion planner.
    let pre_scan_work = curve.bezier_pre_scan_work()?;
    let work_before_plan_scan = source_admission_work
        .checked_add(pre_scan_work)
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-curve admission and plan-scan work overflows u128".to_string(),
        })?;
    enforce_base_work(work_before_plan_scan, "curve plan scan")?;
    let Some(conversion) = curve.bezier_conversion_plan_with_poll(should_cancel)? else {
        return Ok(ClosestWorkRun::Cancelled);
    };
    let knots = curve.knots();
    let order = (knots.degree() as u128)
        .checked_add(1)
        .ok_or_else(|| NurbsError::Domain {
            what: "closest-curve order accounting overflows u128".to_string(),
        })?;
    let mut seed_leaves = 0usize;
    let mut operations_since_poll = 0usize;
    for pair in knots.knots().windows(2) {
        if pair[1] > pair[0] {
            seed_leaves = seed_leaves
                .checked_add(1)
                .ok_or_else(|| NurbsError::Domain {
                    what: "closest-curve seed-leaf count overflows usize".to_string(),
                })?;
        }
        if closest_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(ClosestWorkRun::Cancelled);
        }
    }
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
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    Ok(ClosestWorkRun::Complete(CurveClosestPlan {
        conversion,
        seed_leaves,
        queue_capacity,
    }))
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
    frontier_bytes: u128,
}

#[derive(Debug, Clone, Copy)]
struct SurfaceClosestExecutionPlan {
    base: SurfaceBasePlan,
    subdivision: SurfaceClosestPlan,
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
        frontier_bytes,
    })
}

fn preflight_surface_closest(
    surface: AdmittedNurbsSurface<'_, f64>,
    max_splits: u32,
    source_admission_work: u128,
) -> Result<SurfaceClosestExecutionPlan, NurbsError> {
    let mut never_cancel = || false;
    match preflight_surface_closest_with_poll(
        surface,
        max_splits,
        source_admission_work,
        &mut never_cancel,
    )? {
        ClosestWorkRun::Complete(plan) => Ok(plan),
        ClosestWorkRun::Cancelled => Err(NurbsError::Domain {
            what: "non-cancelling closest-surface preflight observed cancellation".to_string(),
        }),
    }
}

fn preflight_surface_closest_with_poll(
    surface: AdmittedNurbsSurface<'_, f64>,
    max_splits: u32,
    source_admission_work: u128,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<ClosestWorkRun<SurfaceClosestExecutionPlan>, NurbsError> {
    let base = match surface_base_plan_from_admitted_with_poll(
        surface,
        source_admission_work,
        should_cancel,
    )? {
        ClosestWorkRun::Complete(plan) => plan,
        ClosestWorkRun::Cancelled => return Ok(ClosestWorkRun::Cancelled),
    };
    enforce_base_work(base.work_units, "surface")?;
    let subdivision =
        surface_closest_plan(base.order_u, base.order_v, base.seed_leaves, max_splits)?;
    enforce_base_work(
        base.work_units
            .checked_add(subdivision.seed_heap_work)
            .ok_or_else(|| NurbsError::Domain {
                what: "closest-surface base and seed-heap work overflows u128".to_string(),
            })?,
        "surface",
    )?;
    let final_eval_workspace_bytes =
        surface_final_eval_workspace_bytes(base.order_u, base.order_v)?;
    enforce_surface_retained_envelope(
        base.source_bytes,
        base.conversion_peak_allocated_bytes,
        base.converted_bytes,
        subdivision.frontier_bytes,
        final_eval_workspace_bytes,
    )?;
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    Ok(ClosestWorkRun::Complete(SurfaceClosestExecutionPlan {
        base,
        subdivision,
    }))
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

fn hull_lower_bound_with_poll<'a>(
    q: [f64; 3],
    cps: impl Iterator<Item = &'a [f64; 4]>,
    should_cancel: &mut impl FnMut() -> bool,
) -> ClosestWorkRun<f64> {
    if should_cancel() {
        return ClosestWorkRun::Cancelled;
    }
    let mut min = [f64::INFINITY; 3];
    let mut max = [f64::NEG_INFINITY; 3];
    let mut operations_since_poll = 0usize;
    for h in cps {
        let c = cartesian(h);
        for k in 0..3 {
            min[k] = min[k].min(c[k]);
            max[k] = max[k].max(c[k]);
            if closest_poll_due(&mut operations_since_poll, should_cancel) {
                return ClosestWorkRun::Cancelled;
            }
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
        if closest_poll_due(&mut operations_since_poll, should_cancel) {
            return ClosestWorkRun::Cancelled;
        }
    }
    let lower = norm3(gaps);
    if should_cancel() {
        return ClosestWorkRun::Cancelled;
    }
    ClosestWorkRun::Complete(lower)
}

fn hull_lower_bound<'a>(q: [f64; 3], cps: impl Iterator<Item = &'a [f64; 4]>) -> f64 {
    let mut never_cancel = || false;
    match hull_lower_bound_with_poll(q, cps, &mut never_cancel) {
        ClosestWorkRun::Complete(lower) => lower,
        ClosestWorkRun::Cancelled => unreachable!("non-cancelling hull scan observed cancellation"),
    }
}

/// De Casteljau split of a homogeneous Bézier control net at 1/2.
fn split_bezier_with_poll(
    cps: &[[f64; 4]],
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<ClosestWorkRun<(Vec<[f64; 4]>, Vec<[f64; 4]>)>, NurbsError> {
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    let n = cps.len();
    let mut tri = Vec::new();
    tri.try_reserve_exact(n).map_err(|_| NurbsError::Domain {
        what: "Bezier split triangle allocation was refused".to_string(),
    })?;
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    let mut left = Vec::new();
    left.try_reserve_exact(n).map_err(|_| NurbsError::Domain {
        what: "Bezier split left-boundary allocation was refused".to_string(),
    })?;
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    let mut right = Vec::new();
    right.try_reserve_exact(n).map_err(|_| NurbsError::Domain {
        what: "Bezier split right-boundary allocation was refused".to_string(),
    })?;
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    let mut operations_since_poll = 0usize;
    for &control in cps {
        tri.push(control);
        if closest_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(ClosestWorkRun::Cancelled);
        }
    }
    for _ in 0..n {
        right.push([0.0f64; 4]);
        if closest_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(ClosestWorkRun::Cancelled);
        }
    }
    left.push(tri[0]);
    right[n - 1] = tri[n - 1];
    for level in 1..n {
        for i in 0..n - level {
            let (head, tail) = tri.split_at_mut(i + 1);
            for (x, &y) in head[i].iter_mut().zip(&tail[0]) {
                *x = f64::midpoint(*x, y);
                if closest_poll_due(&mut operations_since_poll, should_cancel) {
                    return Ok(ClosestWorkRun::Cancelled);
                }
            }
        }
        left.push(tri[0]);
        right[n - 1 - level] = tri[n - 1 - level];
        if closest_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(ClosestWorkRun::Cancelled);
        }
    }
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    Ok(ClosestWorkRun::Complete((left, right)))
}

fn split_bezier(cps: &[[f64; 4]]) -> Result<(Vec<[f64; 4]>, Vec<[f64; 4]>), NurbsError> {
    let mut never_cancel = || false;
    match split_bezier_with_poll(cps, &mut never_cancel)? {
        ClosestWorkRun::Complete(halves) => Ok(halves),
        ClosestWorkRun::Cancelled => Err(NurbsError::Domain {
            what: "non-cancelling Bezier split observed cancellation".to_string(),
        }),
    }
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

/// Measured closest-point search on an owning curve with bounded cancellation
/// polling and transactional estimate publication.
///
/// Request shape and count-derived source-work refusals retain their
/// synchronous precedence. The same `Cx` then spans structural admission and
/// the admitted measured search, while the aggregate preflight includes the
/// source-admission charge. Cancellation publishes neither admitted authority
/// nor a partial frontier-derived estimate.
///
/// # Errors
/// Returns the synchronous owning search's request, work, retained-memory,
/// allocation, structure, and finite-arithmetic refusals when they win before
/// an observed cancellation.
pub fn closest_point_curve_with_cx(
    curve: &NurbsCurve<f64, 3>,
    q: [f64; 3],
    tol: f64,
    max_splits: u32,
    cx: &Cx<'_>,
) -> Result<ClosestPointRun, NurbsError> {
    validate_closest_request(q, tol, max_splits)?;
    let source_admission_work = curve_source_admission_work(curve)?;
    enforce_base_work(source_admission_work, "curve source admission")?;
    let admitted = match curve.admit_with_cx(cx)? {
        CurveAdmissionRun::Complete { admitted } => admitted,
        CurveAdmissionRun::Cancelled => return Ok(ClosestPointRun::Cancelled),
    };
    admitted.closest_point_with_cx_after_request(q, tol, max_splits, source_admission_work, cx)
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

    /// Measured closest-point search with bounded cancellation polling and
    /// transactional estimate publication.
    ///
    /// Request shape and the constant-time pre-scan work gate retain their
    /// synchronous precedence. One `Cx` then spans Bezier planning/conversion,
    /// seed construction, branch-and-bound splitting, optional Newton polish,
    /// final evaluation, and publication. Cancellation during optional polish
    /// is not downgraded to an ordinary polish miss. This method does not
    /// consume the `Cx` budget or finalize its executor scope.
    ///
    /// # Errors
    /// Returns the synchronous search's request, work, retained-memory,
    /// allocation, structure, and finite-arithmetic refusals when they win
    /// before an observed cancellation.
    pub fn closest_point_with_cx(
        &self,
        q: [f64; 3],
        tol: f64,
        max_splits: u32,
        cx: &Cx<'_>,
    ) -> Result<ClosestPointRun, NurbsError> {
        validate_closest_request(q, tol, max_splits)?;
        self.closest_point_with_cx_after_request(q, tol, max_splits, 0, cx)
    }

    fn closest_point_with_cx_after_request(
        &self,
        q: [f64; 3],
        tol: f64,
        max_splits: u32,
        source_admission_work: u128,
        cx: &Cx<'_>,
    ) -> Result<ClosestPointRun, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        let plan = match preflight_curve_closest_with_poll(
            *self,
            max_splits,
            source_admission_work,
            &mut should_cancel,
        )? {
            ClosestWorkRun::Complete(plan) => plan,
            ClosestWorkRun::Cancelled => return Ok(ClosestPointRun::Cancelled),
        };
        let mut kernels = CurveClosestKernels {
            conversion: || match self.to_bezier_form_with_cx(cx)? {
                CurveBezierRun::Complete { curve } => Ok(ClosestWorkRun::Complete(curve)),
                CurveBezierRun::Cancelled => Ok(ClosestWorkRun::Cancelled),
            },
            derivatives: |t, order| match self.derivatives_with_cx(t, order, cx)? {
                CurveDerivativesRun::Complete { derivatives } => {
                    Ok(ClosestWorkRun::Complete(derivatives))
                }
                CurveDerivativesRun::Cancelled => Ok(ClosestWorkRun::Cancelled),
            },
            evaluation: |t| match self.eval_with_cx(t, cx)? {
                CurveEvaluationRun::Complete { point } => Ok(ClosestWorkRun::Complete(point)),
                CurveEvaluationRun::Cancelled => Ok(ClosestWorkRun::Cancelled),
            },
        };
        match closest_point_curve_after_plan_with_poll(
            *self,
            q,
            tol,
            max_splits,
            plan,
            &mut kernels,
            &mut should_cancel,
        )? {
            ClosestWorkRun::Complete(estimate) => Ok(ClosestPointRun::Complete { estimate }),
            ClosestWorkRun::Cancelled => Ok(ClosestPointRun::Cancelled),
        }
    }
}

struct CurveClosestKernels<C, D, E> {
    conversion: C,
    derivatives: D,
    evaluation: E,
}

fn closest_point_curve_after_request(
    curve: AdmittedNurbsCurve<'_, f64, 3>,
    q: [f64; 3],
    tol: f64,
    max_splits: u32,
    source_admission_work: u128,
) -> Result<DistanceBracketEstimate, NurbsError> {
    let plan = preflight_curve_closest(curve, max_splits, source_admission_work)?;
    let mut kernels = CurveClosestKernels {
        conversion: || curve.to_bezier_form().map(ClosestWorkRun::Complete),
        derivatives: |t, order| curve.derivatives(t, order).map(ClosestWorkRun::Complete),
        evaluation: |t| curve.eval(t).map(ClosestWorkRun::Complete),
    };
    let mut never_cancel = || false;
    match closest_point_curve_after_plan_with_poll(
        curve,
        q,
        tol,
        max_splits,
        plan,
        &mut kernels,
        &mut never_cancel,
    )? {
        ClosestWorkRun::Complete(estimate) => Ok(estimate),
        ClosestWorkRun::Cancelled => Err(NurbsError::Domain {
            what: "non-cancelling closest-curve search observed cancellation".to_string(),
        }),
    }
}

fn closest_point_curve_after_plan_with_poll<C, D, E>(
    curve: AdmittedNurbsCurve<'_, f64, 3>,
    q: [f64; 3],
    tol: f64,
    max_splits: u32,
    plan: CurveClosestPlan,
    kernels: &mut CurveClosestKernels<C, D, E>,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<ClosestWorkRun<DistanceBracketEstimate>, NurbsError>
where
    C: FnMut() -> Result<ClosestWorkRun<NurbsCurve<f64, 3>>, NurbsError>,
    D: FnMut(f64, usize) -> Result<ClosestWorkRun<Vec<[f64; 3]>>, NurbsError>,
    E: FnMut(f64) -> Result<ClosestWorkRun<[f64; 3]>, NurbsError>,
{
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    let bezier_curve = match (kernels.conversion)()? {
        ClosestWorkRun::Complete(curve) => curve,
        ClosestWorkRun::Cancelled => return Ok(ClosestWorkRun::Cancelled),
    };
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
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
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    queue
        .try_reserve_exact(plan.queue_capacity)
        .map_err(|_| NurbsError::Domain {
            what: format!(
                "closest-curve queue allocation was refused for admitted capacity {}",
                plan.queue_capacity
            ),
        })?;
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    let mut next_logical_id = 0u64;
    let mut upper = f64::INFINITY;
    let mut best_t = knots.domain().0;
    let mut operations_since_poll = 0usize;
    for span in p..knots.control_count() {
        let (t0, t1) = (entries[span], entries[span + 1]);
        if closest_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(ClosestWorkRun::Cancelled);
        }
        if t1 <= t0 {
            continue;
        }
        let mut cpw = Vec::new();
        if should_cancel() {
            return Ok(ClosestWorkRun::Cancelled);
        }
        cpw.try_reserve_exact(p + 1)
            .map_err(|_| NurbsError::Domain {
                what: "closest-curve seed control allocation was refused".to_string(),
            })?;
        if should_cancel() {
            return Ok(ClosestWorkRun::Cancelled);
        }
        for &control in &controls[span - p..=span] {
            cpw.push(control);
            if closest_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(ClosestWorkRun::Cancelled);
            }
        }
        for (h, tt) in [(&cpw[0], t0), (&cpw[p], t1)] {
            let d = dist3(cartesian(h), q);
            if d < upper {
                upper = d;
                best_t = tt;
            }
            if closest_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(ClosestWorkRun::Cancelled);
            }
        }
        let lb = match hull_lower_bound_with_poll(q, cpw.iter(), should_cancel) {
            ClosestWorkRun::Complete(lower) => lower,
            ClosestWorkRun::Cancelled => return Ok(ClosestWorkRun::Cancelled),
        };
        if should_cancel() {
            return Ok(ClosestWorkRun::Cancelled);
        }
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
        if should_cancel() {
            return Ok(ClosestWorkRun::Cancelled);
        }
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
    if queue.is_empty() || !upper.is_finite() {
        return Err(NurbsError::Domain {
            what: "closest-curve initial distance bounds are not finite".to_string(),
        });
    }
    for entry in &queue {
        if !entry.key.is_finite() {
            return Err(NurbsError::Domain {
                what: "closest-curve initial distance bounds are not finite".to_string(),
            });
        }
        if closest_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(ClosestWorkRun::Cancelled);
        }
    }
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    let mut iterations = 0u32;
    while iterations < max_splits {
        if should_cancel() {
            return Ok(ClosestWorkRun::Cancelled);
        }
        let Some(entry) = queue.peek() else {
            break;
        };
        if should_cancel() {
            return Ok(ClosestWorkRun::Cancelled);
        }
        if upper - entry.key <= tol {
            break; // bracket closed
        }
        let Some(entry) = queue.pop() else {
            break;
        };
        if should_cancel() {
            return Ok(ClosestWorkRun::Cancelled);
        }
        let seg = entry.value;
        let tm = f64::midpoint(seg.t0, seg.t1);
        if tm == seg.t0 || tm == seg.t1 {
            // De Casteljau would still create a geometric half-split even
            // though the recorded parameter interval cannot shrink. Retain
            // the leaf so its lower estimate remains part of the bracket.
            if should_cancel() {
                return Ok(ClosestWorkRun::Cancelled);
            }
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
            if should_cancel() {
                return Ok(ClosestWorkRun::Cancelled);
            }
            break;
        }
        let (l, r) = match split_bezier_with_poll(&seg.cpw, should_cancel)? {
            ClosestWorkRun::Complete(halves) => halves,
            ClosestWorkRun::Cancelled => return Ok(ClosestWorkRun::Cancelled),
        };
        if should_cancel() {
            return Ok(ClosestWorkRun::Cancelled);
        }
        // The split junction is C(mid): a free upper-bound sample.
        let d = dist3(cartesian(&l[l.len() - 1]), q);
        if d < upper {
            upper = d;
            best_t = tm;
        }
        for (cpw, t0, t1) in [(l, seg.t0, tm), (r, tm, seg.t1)] {
            let lb = match hull_lower_bound_with_poll(q, cpw.iter(), should_cancel) {
                ClosestWorkRun::Complete(lower) => lower,
                ClosestWorkRun::Cancelled => return Ok(ClosestWorkRun::Cancelled),
            };
            if !lb.is_finite() {
                return Err(NurbsError::Domain {
                    what: "closest-curve child distance bound is not finite".to_string(),
                });
            }
            if lb < upper {
                if should_cancel() {
                    return Ok(ClosestWorkRun::Cancelled);
                }
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
                if should_cancel() {
                    return Ok(ClosestWorkRun::Cancelled);
                }
                next_logical_id += 1;
            }
        }
        iterations += 1;
    }
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    let lower = queue.peek().map_or(upper, |entry| entry.key);
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    // The bracket no longer needs the converted generation or frontier.
    // Release both phases before optional derivative polish so the aggregate
    // peak is max(conversion, search, polish), not their sum.
    drop(queue);
    drop(bezier_curve);
    // Allocator and destructor latency is outside this measured poll bound;
    // gate the next phase immediately after both owned generations are gone.
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    // Newton polish on g(t) = (C − q)·C' sharpens the upper bound.
    let (dlo, dhi) = curve.knots().domain();
    let mut t = best_t;
    for _ in 0..12 {
        // Polishing is an optional improvement to an already valid measured
        // B&B estimate. A derivative can be undefined at a repeated knot (or
        // unavailable in the legacy jet API); that must not erase the retained
        // geometric witness and bracket.
        if should_cancel() {
            return Ok(ClosestWorkRun::Cancelled);
        }
        let ders = match (kernels.derivatives)(t, 2) {
            Ok(ClosestWorkRun::Complete(derivatives)) => derivatives,
            Ok(ClosestWorkRun::Cancelled) => return Ok(ClosestWorkRun::Cancelled),
            Err(_) => break,
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
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    let (upper, best_t) = match (kernels.evaluation)(t) {
        Ok(ClosestWorkRun::Complete(point)) => {
            let polished = dist3(point, q);
            if polished.is_finite() && polished < upper {
                (polished, t)
            } else {
                (upper, best_t)
            }
        }
        Ok(ClosestWorkRun::Cancelled) => return Ok(ClosestWorkRun::Cancelled),
        Err(_) => {
            // Newton is optional polish. Its failure cannot erase the finite
            // branch-and-bound witness retained above.
            (upper, best_t)
        }
    };
    let estimate = DistanceBracketEstimate {
        lower: lower.min(upper),
        upper,
        param: [best_t, 0.0],
        iterations,
    };
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    Ok(ClosestWorkRun::Complete(estimate))
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

fn patch_lb_with_poll(
    q: [f64; 3],
    net: &[Vec<[f64; 4]>],
    should_cancel: &mut impl FnMut() -> bool,
) -> ClosestWorkRun<f64> {
    hull_lower_bound_with_poll(q, net.iter().flatten(), should_cancel)
}

fn zero_net_with_poll(
    rows: usize,
    cols: usize,
    stage: &str,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<ClosestWorkRun<Net>, NurbsError> {
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    let mut net = Vec::new();
    net.try_reserve_exact(rows)
        .map_err(|_| NurbsError::Domain {
            what: format!("{stage} row-table allocation was refused"),
        })?;
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    let mut operations_since_poll = 0usize;
    for _ in 0..rows {
        if should_cancel() {
            return Ok(ClosestWorkRun::Cancelled);
        }
        let mut row = Vec::new();
        row.try_reserve_exact(cols)
            .map_err(|_| NurbsError::Domain {
                what: format!("{stage} row allocation was refused"),
            })?;
        if should_cancel() {
            return Ok(ClosestWorkRun::Cancelled);
        }
        for _ in 0..cols {
            row.push([0.0; 4]);
            if closest_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(ClosestWorkRun::Cancelled);
            }
        }
        net.push(row);
        if closest_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(ClosestWorkRun::Cancelled);
        }
    }
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    Ok(ClosestWorkRun::Complete(net))
}

fn split_patch_u_with_poll(
    net: &[Vec<[f64; 4]>],
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<ClosestWorkRun<(Net, Net)>, NurbsError> {
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    // Split every v-column along u (rows are u direction).
    let rows = net.len();
    let cols = net[0].len();
    let mut left = match zero_net_with_poll(rows, cols, "u-split left patch", should_cancel)? {
        ClosestWorkRun::Complete(net) => net,
        ClosestWorkRun::Cancelled => return Ok(ClosestWorkRun::Cancelled),
    };
    let mut right = match zero_net_with_poll(rows, cols, "u-split right patch", should_cancel)? {
        ClosestWorkRun::Complete(net) => net,
        ClosestWorkRun::Cancelled => return Ok(ClosestWorkRun::Cancelled),
    };
    let mut operations_since_poll = 0usize;
    for j in 0..cols {
        if closest_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(ClosestWorkRun::Cancelled);
        }
        let mut col = Vec::new();
        if should_cancel() {
            return Ok(ClosestWorkRun::Cancelled);
        }
        col.try_reserve_exact(rows)
            .map_err(|_| NurbsError::Domain {
                what: "u-split column allocation was refused".to_string(),
            })?;
        if should_cancel() {
            return Ok(ClosestWorkRun::Cancelled);
        }
        for row in net {
            col.push(row[j]);
            if closest_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(ClosestWorkRun::Cancelled);
            }
        }
        let (l, r) = match split_bezier_with_poll(&col, should_cancel)? {
            ClosestWorkRun::Complete(halves) => halves,
            ClosestWorkRun::Cancelled => return Ok(ClosestWorkRun::Cancelled),
        };
        for i in 0..rows {
            left[i][j] = l[i];
            if closest_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(ClosestWorkRun::Cancelled);
            }
            right[i][j] = r[i];
            if closest_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(ClosestWorkRun::Cancelled);
            }
        }
    }
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    Ok(ClosestWorkRun::Complete((left, right)))
}

fn split_patch_v_with_poll(
    net: &[Vec<[f64; 4]>],
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<ClosestWorkRun<(Net, Net)>, NurbsError> {
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    let (mut left, mut right) = (Vec::new(), Vec::new());
    left.try_reserve_exact(net.len())
        .map_err(|_| NurbsError::Domain {
            what: "v-split left row-table allocation was refused".to_string(),
        })?;
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    right
        .try_reserve_exact(net.len())
        .map_err(|_| NurbsError::Domain {
            what: "v-split right row-table allocation was refused".to_string(),
        })?;
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    let mut operations_since_poll = 0usize;
    for row in net {
        if closest_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(ClosestWorkRun::Cancelled);
        }
        let (l, r) = match split_bezier_with_poll(row, should_cancel)? {
            ClosestWorkRun::Complete(halves) => halves,
            ClosestWorkRun::Cancelled => return Ok(ClosestWorkRun::Cancelled),
        };
        left.push(l);
        if closest_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(ClosestWorkRun::Cancelled);
        }
        right.push(r);
        if closest_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(ClosestWorkRun::Cancelled);
        }
    }
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    Ok(ClosestWorkRun::Complete((left, right)))
}

/// Decompose a surface to Bézier patches via repeated knot insertion.
fn to_bezier_surface(
    surface: AdmittedNurbsSurface<'_, f64>,
) -> Result<NurbsSurface<f64>, NurbsError> {
    let mut never_cancel = || false;
    match to_bezier_surface_with_poll(surface, &mut never_cancel)? {
        ClosestWorkRun::Complete(surface) => Ok(surface),
        ClosestWorkRun::Cancelled => Err(NurbsError::Domain {
            what: "non-cancelling surface Bezier conversion observed cancellation".to_string(),
        }),
    }
}

fn to_bezier_surface_with_poll(
    surface: AdmittedNurbsSurface<'_, f64>,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<ClosestWorkRun<NurbsSurface<f64>>, NurbsError> {
    let mut work = match surface.try_clone_with_poll(should_cancel)? {
        Some(surface) => surface,
        None => return Ok(ClosestWorkRun::Cancelled),
    };
    let mut operations_since_poll = 0usize;
    loop {
        let mut inserted = false;
        for dir_u in [true, false] {
            if should_cancel() {
                return Ok(ClosestWorkRun::Cancelled);
            }
            let target = {
                let (kv, p) = if dir_u {
                    (&work.knots_u, work.knots_u.degree)
                } else {
                    (&work.knots_v, work.knots_v.degree)
                };
                // The cloned source and every derived publication are already
                // validated; keep that authority for domain lookup as well as
                // the subsequent insertion instead of rescanning structure.
                let (lo, hi) = kv.admitted_after_validation().domain();
                let mut target = None;
                let mut run_start = 0usize;
                while run_start < kv.knots.len() {
                    let t = kv.knots[run_start];
                    let mut run_end = run_start + 1;
                    while run_end < kv.knots.len() && kv.knots[run_end] == t {
                        run_end += 1;
                        if closest_poll_due(&mut operations_since_poll, should_cancel) {
                            return Ok(ClosestWorkRun::Cancelled);
                        }
                    }
                    if t > lo && t < hi && run_end - run_start < p {
                        target = Some(t);
                        break;
                    }
                    if closest_poll_due(&mut operations_since_poll, should_cancel) {
                        return Ok(ClosestWorkRun::Cancelled);
                    }
                    run_start = run_end;
                }
                target
            };
            if let Some(t) = target {
                let admitted_work = work.admitted_after_validation();
                let inserted_surface = if dir_u {
                    admitted_work.insert_knot_u_with_poll(t, should_cancel)?
                } else {
                    admitted_work.insert_knot_v_with_poll(t, should_cancel)?
                };
                work = match inserted_surface {
                    crate::surface::SurfaceInsertionRun::Complete { surface } => surface,
                    crate::surface::SurfaceInsertionRun::Cancelled => {
                        return Ok(ClosestWorkRun::Cancelled);
                    }
                };
                // Dropping the superseded generation is not preemptible; gate
                // the next phase immediately after replacement.
                if should_cancel() {
                    return Ok(ClosestWorkRun::Cancelled);
                }
                inserted = true;
            }
        }
        if !inserted {
            if should_cancel() {
                return Ok(ClosestWorkRun::Cancelled);
            }
            return Ok(ClosestWorkRun::Complete(work));
        }
    }
}

/// Seed the patch queue from a Bézier-form surface.
fn seed_patches_with_poll(
    work: &NurbsSurface<f64>,
    q: [f64; 3],
    queue: &mut BinaryHeap<MinEntry<Patch>>,
    queue_capacity: usize,
    next_logical_id: &mut u64,
    upper: &mut f64,
    best: &mut [f64; 2],
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<ClosestWorkRun<()>, NurbsError> {
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    let (pu, pv) = (work.knots_u.degree, work.knots_v.degree);
    let mut operations_since_poll = 0usize;
    for su in pu..work.knots_u.control_count() {
        let (u0, u1) = (work.knots_u.knots[su], work.knots_u.knots[su + 1]);
        if closest_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(ClosestWorkRun::Cancelled);
        }
        if u1 <= u0 {
            continue;
        }
        for sv in pv..work.knots_v.control_count() {
            let (v0, v1) = (work.knots_v.knots[sv], work.knots_v.knots[sv + 1]);
            if closest_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(ClosestWorkRun::Cancelled);
            }
            if v1 <= v0 {
                continue;
            }
            let mut net = match zero_net_with_poll(
                pu + 1,
                pv + 1,
                "closest-surface seed patch",
                should_cancel,
            )? {
                ClosestWorkRun::Complete(net) => net,
                ClosestWorkRun::Cancelled => return Ok(ClosestWorkRun::Cancelled),
            };
            for (target_row, source_row) in net.iter_mut().zip(work.cpw[su - pu..=su].iter()) {
                for (target, &source) in target_row.iter_mut().zip(&source_row[sv - pv..=sv]) {
                    *target = source;
                    if closest_poll_due(&mut operations_since_poll, should_cancel) {
                        return Ok(ClosestWorkRun::Cancelled);
                    }
                }
            }
            let d = dist3(cartesian(&net[0][0]), q);
            if d < *upper {
                *upper = d;
                *best = [u0, v0];
            }
            let lb = match patch_lb_with_poll(q, &net, should_cancel) {
                ClosestWorkRun::Complete(lower) => lower,
                ClosestWorkRun::Cancelled => return Ok(ClosestWorkRun::Cancelled),
            };
            if should_cancel() {
                return Ok(ClosestWorkRun::Cancelled);
            }
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
            if should_cancel() {
                return Ok(ClosestWorkRun::Cancelled);
            }
            *next_logical_id += 1;
        }
    }
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    Ok(ClosestWorkRun::Complete(()))
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

/// Measured closest-point search on an owning surface with bounded
/// cancellation polling and transactional estimate publication.
///
/// Request shape and count-derived source-work refusals retain their
/// synchronous precedence. The same `Cx` then spans structural admission and
/// the admitted measured search, while the aggregate preflight includes the
/// source-admission charge. Cancellation publishes neither admitted authority
/// nor a partial frontier-derived estimate.
///
/// # Errors
/// Returns the synchronous owning search's request, work, retained-memory,
/// allocation, structure, and finite-arithmetic refusals when they win before
/// an observed cancellation.
pub fn closest_point_surface_with_cx(
    surface: &NurbsSurface<f64>,
    q: [f64; 3],
    tol: f64,
    max_splits: u32,
    cx: &Cx<'_>,
) -> Result<ClosestPointRun, NurbsError> {
    validate_closest_request(q, tol, max_splits)?;
    let source_admission_work = surface_source_admission_work(surface)?;
    enforce_base_work(source_admission_work, "surface source admission")?;
    let admitted = match surface.admit_with_cx(cx)? {
        SurfaceAdmissionRun::Complete { admitted } => admitted,
        SurfaceAdmissionRun::Cancelled => return Ok(ClosestPointRun::Cancelled),
    };
    admitted.closest_point_with_cx_after_request(q, tol, max_splits, source_admission_work, cx)
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

    /// Measured surface closest-point search with bounded cancellation
    /// polling and transactional estimate publication.
    ///
    /// Request shape and the count-only pre-scan work gate retain synchronous
    /// precedence. One `Cx` then spans planning, exact tensor conversion,
    /// patch seeding/splitting, optional center evaluation, cleanup, and final
    /// publication. This method does not consume the `Cx` budget or finalize
    /// its executor scope.
    ///
    /// # Errors
    /// Returns the synchronous search's work, retained-memory, allocation,
    /// structure, and finite-arithmetic refusals when they win before an
    /// observed cancellation.
    pub fn closest_point_with_cx(
        &self,
        q: [f64; 3],
        tol: f64,
        max_splits: u32,
        cx: &Cx<'_>,
    ) -> Result<ClosestPointRun, NurbsError> {
        validate_closest_request(q, tol, max_splits)?;
        self.closest_point_with_cx_after_request(q, tol, max_splits, 0, cx)
    }

    fn closest_point_with_cx_after_request(
        &self,
        q: [f64; 3],
        tol: f64,
        max_splits: u32,
        source_admission_work: u128,
        cx: &Cx<'_>,
    ) -> Result<ClosestPointRun, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        let plan = match preflight_surface_closest_with_poll(
            *self,
            max_splits,
            source_admission_work,
            &mut should_cancel,
        )? {
            ClosestWorkRun::Complete(plan) => plan,
            ClosestWorkRun::Cancelled => return Ok(ClosestPointRun::Cancelled),
        };
        let mut kernels = SurfaceClosestKernels {
            conversion: || {
                let mut conversion_cancel = || cx.checkpoint().is_err();
                to_bezier_surface_with_poll(*self, &mut conversion_cancel)
            },
            evaluation: |u, v| match self.eval_with_cx(u, v, cx)? {
                SurfaceEvaluationRun::Complete { point } => Ok(ClosestWorkRun::Complete(point)),
                SurfaceEvaluationRun::Cancelled => Ok(ClosestWorkRun::Cancelled),
            },
        };
        match closest_point_surface_after_plan_with_poll(
            *self,
            q,
            tol,
            max_splits,
            plan,
            &mut kernels,
            &mut should_cancel,
        )? {
            ClosestWorkRun::Complete(estimate) => Ok(ClosestPointRun::Complete { estimate }),
            ClosestWorkRun::Cancelled => Ok(ClosestPointRun::Cancelled),
        }
    }
}

struct SurfaceClosestKernels<C, E> {
    conversion: C,
    evaluation: E,
}

fn closest_point_surface_after_request(
    surface: AdmittedNurbsSurface<'_, f64>,
    q: [f64; 3],
    tol: f64,
    max_splits: u32,
    source_admission_work: u128,
) -> Result<DistanceBracketEstimate, NurbsError> {
    let plan = preflight_surface_closest(surface, max_splits, source_admission_work)?;
    let mut kernels = SurfaceClosestKernels {
        conversion: || to_bezier_surface(surface).map(ClosestWorkRun::Complete),
        evaluation: |u, v| surface.eval(u, v).map(ClosestWorkRun::Complete),
    };
    let mut never_cancel = || false;
    match closest_point_surface_after_plan_with_poll(
        surface,
        q,
        tol,
        max_splits,
        plan,
        &mut kernels,
        &mut never_cancel,
    )? {
        ClosestWorkRun::Complete(estimate) => Ok(estimate),
        ClosestWorkRun::Cancelled => Err(NurbsError::Domain {
            what: "non-cancelling closest-surface search observed cancellation".to_string(),
        }),
    }
}

fn closest_point_surface_after_plan_with_poll<C, E>(
    surface: AdmittedNurbsSurface<'_, f64>,
    q: [f64; 3],
    tol: f64,
    max_splits: u32,
    execution_plan: SurfaceClosestExecutionPlan,
    kernels: &mut SurfaceClosestKernels<C, E>,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<ClosestWorkRun<DistanceBracketEstimate>, NurbsError>
where
    C: FnMut() -> Result<ClosestWorkRun<NurbsSurface<f64>>, NurbsError>,
    E: FnMut(f64, f64) -> Result<ClosestWorkRun<[f64; 3]>, NurbsError>,
{
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    let work = match (kernels.conversion)()? {
        ClosestWorkRun::Complete(surface) => surface,
        ClosestWorkRun::Cancelled => return Ok(ClosestWorkRun::Cancelled),
    };
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    let base_plan = execution_plan.base;
    let plan = execution_plan.subdivision;
    let actual_control_count_v = work.cpw.first().map_or(0, Vec::len) as u128;
    if work.knots_u.knots.len() as u128 != base_plan.final_knot_count_u
        || work.knots_v.knots.len() as u128 != base_plan.final_knot_count_v
        || work.cpw.len() as u128 != base_plan.final_control_count_u
        || actual_control_count_v != base_plan.final_control_count_v
    {
        return Err(NurbsError::Structure {
            what: "closest-surface conversion did not match its admitted plan".to_string(),
        });
    }
    let mut queue: BinaryHeap<MinEntry<Patch>> = BinaryHeap::new();
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    queue
        .try_reserve_exact(plan.queue_capacity)
        .map_err(|_| NurbsError::Domain {
            what: format!(
                "closest-surface queue allocation was refused for admitted capacity {}",
                plan.queue_capacity
            ),
        })?;
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    let mut next_logical_id = 0u64;
    let mut upper = f64::INFINITY;
    let mut best = [surface.knots_u().domain().0, surface.knots_v().domain().0];
    match seed_patches_with_poll(
        &work,
        q,
        &mut queue,
        plan.queue_capacity,
        &mut next_logical_id,
        &mut upper,
        &mut best,
        should_cancel,
    )? {
        ClosestWorkRun::Complete(()) => {}
        ClosestWorkRun::Cancelled => return Ok(ClosestWorkRun::Cancelled),
    }
    if queue.len() != plan.seed_leaves {
        return Err(NurbsError::Structure {
            what: format!(
                "closest-surface conversion produced {} seed leaves after admitting {}",
                queue.len(),
                plan.seed_leaves
            ),
        });
    }
    if queue.is_empty() || !upper.is_finite() {
        return Err(NurbsError::Domain {
            what: "closest-surface initial distance bounds are not finite".to_string(),
        });
    }
    let mut operations_since_poll = 0usize;
    for entry in &queue {
        if !entry.key.is_finite() {
            return Err(NurbsError::Domain {
                what: "closest-surface initial distance bounds are not finite".to_string(),
            });
        }
        if closest_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(ClosestWorkRun::Cancelled);
        }
    }
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    let mut iterations = 0u32;
    while iterations < max_splits {
        if should_cancel() {
            return Ok(ClosestWorkRun::Cancelled);
        }
        let Some(entry) = queue.peek() else {
            break;
        };
        if should_cancel() {
            return Ok(ClosestWorkRun::Cancelled);
        }
        if upper - entry.key <= tol {
            break;
        }
        let Some(entry) = queue.pop() else {
            break;
        };
        if should_cancel() {
            return Ok(ClosestWorkRun::Cancelled);
        }
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
                if should_cancel() {
                    return Ok(ClosestWorkRun::Cancelled);
                }
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
                if should_cancel() {
                    return Ok(ClosestWorkRun::Cancelled);
                }
                break;
            }
        };
        let midpoint = if split_u { midpoint_u } else { midpoint_v };
        if !midpoint.is_finite() {
            if should_cancel() {
                return Ok(ClosestWorkRun::Cancelled);
            }
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
            if should_cancel() {
                return Ok(ClosestWorkRun::Cancelled);
            }
            break;
        }
        let (l, r) = if split_u {
            match split_patch_u_with_poll(&patch.cpw, should_cancel)? {
                ClosestWorkRun::Complete(halves) => halves,
                ClosestWorkRun::Cancelled => return Ok(ClosestWorkRun::Cancelled),
            }
        } else {
            match split_patch_v_with_poll(&patch.cpw, should_cancel)? {
                ClosestWorkRun::Complete(halves) => halves,
                ClosestWorkRun::Cancelled => return Ok(ClosestWorkRun::Cancelled),
            }
        };
        if should_cancel() {
            return Ok(ClosestWorkRun::Cancelled);
        }
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
            let lb = match patch_lb_with_poll(q, &net, should_cancel) {
                ClosestWorkRun::Complete(lower) => lower,
                ClosestWorkRun::Cancelled => return Ok(ClosestWorkRun::Cancelled),
            };
            if !lb.is_finite() {
                return Err(NurbsError::Domain {
                    what: "closest-surface child distance bound is not finite".to_string(),
                });
            }
            if lb < upper {
                if should_cancel() {
                    return Ok(ClosestWorkRun::Cancelled);
                }
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
                if should_cancel() {
                    return Ok(ClosestWorkRun::Cancelled);
                }
                next_logical_id += 1;
            }
        }
        iterations += 1;
    }
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    let lower = queue.peek().map_or(upper, |entry| entry.key);
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    // Sample the current best-lower-estimate patch center for a final
    // evaluated-point improvement. The former midpoint(best,best) expression
    // merely re-evaluated the already retained corner and never sampled a
    // patch center; it also failed to update `param` when the sample improved.
    let candidate = queue.peek().map(|entry| {
        let patch = &entry.value;
        [
            f64::midpoint(patch.u0, patch.u1),
            f64::midpoint(patch.v0, patch.v1),
        ]
    });
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    if let Some(candidate) = candidate {
        match (kernels.evaluation)(candidate[0], candidate[1]) {
            Ok(ClosestWorkRun::Complete(point)) => {
                let distance = dist3(point, q);
                if distance.is_finite() && distance < upper {
                    upper = distance;
                    best = candidate;
                }
            }
            Ok(ClosestWorkRun::Cancelled) => return Ok(ClosestWorkRun::Cancelled),
            Err(_) => {}
        }
    }
    // The retained-envelope model keeps conversion and frontier live through
    // final evaluation. Their destructors are not preemptible; gate the final
    // publication immediately after cleanup.
    drop(queue);
    drop(work);
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    let estimate = DistanceBracketEstimate {
        lower: lower.min(upper),
        upper,
        param: best,
        iterations,
    };
    if should_cancel() {
        return Ok(ClosestWorkRun::Cancelled);
    }
    Ok(ClosestWorkRun::Complete(estimate))
}

#[cfg(test)]
mod tests {
    use super::{
        BinaryHeap, CLOSEST_MAX_BASE_WORK_UNITS, CLOSEST_MAX_RETAINED_BYTES, CLOSEST_MAX_SPLITS,
        ClosestPointRun, ClosestWorkRun, CurveClosestKernels, MinEntry, SurfaceClosestKernels,
        closest_point_curve, closest_point_curve_after_plan_with_poll, closest_point_curve_with_cx,
        closest_point_surface, closest_point_surface_after_plan_with_poll,
        closest_point_surface_with_cx, curve_subdivision_work_per_split,
        enforce_curve_retained_envelope, enforce_surface_retained_envelope,
        hull_lower_bound_with_poll, preflight_curve_closest, preflight_curve_closest_with_poll,
        preflight_surface_closest, preflight_surface_closest_with_poll,
        preflight_surface_subdivision, seed_patches_with_poll, split_bezier_with_poll,
        split_patch_u_with_poll, split_patch_v_with_poll, subdivision_frontier_bytes,
        surface_base_work_units, surface_conversion_peak_allocated_bytes,
        surface_final_eval_workspace_bytes, surface_storage_bytes, to_bezier_surface,
        to_bezier_surface_with_poll,
    };
    use crate::{KnotVector, NurbsCurve, NurbsError, NurbsSurface};
    use asupersync::types::Budget;
    use fs_exec::{CancelGate, Cx, ExecMode, StreamKey};
    use std::mem::size_of;

    fn with_closest_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new_clock_free();
        if cancelled {
            gate.request();
        }
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 0xC105_E570,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    fn quadratic_curve() -> NurbsCurve<f64, 3> {
        NurbsCurve::<f64, 3>::new(
            KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2).expect("quadratic knots"),
            &[[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [2.0, 0.0, 0.0]],
            &[1.0; 3],
        )
        .expect("quadratic curve")
    }

    fn long_linear_curve() -> NurbsCurve<f64, 3> {
        let mut knots = vec![0.0, 0.0];
        knots.extend((1..=128).map(|index| f64::from(index) / 129.0));
        knots.extend([1.0, 1.0]);
        let points: Vec<_> = (0..130)
            .map(|index| [f64::from(index) / 129.0, 0.0, 0.0])
            .collect();
        let weights = vec![1.0; points.len()];
        NurbsCurve::new(
            KnotVector::new(knots, 1).expect("long linear knots"),
            &points,
            &weights,
        )
        .expect("long linear curve")
    }

    fn high_degree_bezier_curve() -> NurbsCurve<f64, 3> {
        let mut knots = vec![0.0; 65];
        knots.extend([1.0; 65]);
        let controls = (0..65)
            .map(|index| [f64::from(index) / 64.0, 0.0, 0.0, 1.0])
            .collect();
        NurbsCurve::from_homogeneous(
            KnotVector::new(knots, 64).expect("degree-64 Bezier knots"),
            controls,
        )
        .expect("degree-64 Bezier curve")
    }

    fn multispan_quadratic_surface() -> NurbsSurface<f64> {
        let knots = KnotVector::new(vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0], 2)
            .expect("two-span quadratic surface knots");
        let points: Vec<Vec<[f64; 3]>> = (0..4)
            .map(|u| {
                (0..4)
                    .map(|v| [f64::from(u) / 3.0, f64::from(v) / 3.0, 0.0])
                    .collect()
            })
            .collect();
        let weights = vec![vec![1.0; 4]; 4];
        NurbsSurface::new(knots.clone(), knots, &points, &weights)
            .expect("quadratic multispan surface")
    }

    fn long_linear_surface() -> NurbsSurface<f64> {
        let mut knots_u = vec![0.0, 0.0];
        knots_u.extend((1..=128).map(|index| f64::from(index) / 129.0));
        knots_u.extend([1.0, 1.0]);
        let knots_v = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("linear V knots");
        let controls = (0..130)
            .map(|index| {
                let u = f64::from(index) / 129.0;
                vec![[u, 0.0, 0.0, 1.0], [u, 1.0, 0.0, 1.0]]
            })
            .collect();
        NurbsSurface::from_homogeneous(
            KnotVector::new(knots_u, 1).expect("long linear U knots"),
            knots_v,
            controls,
        )
        .expect("long linear surface")
    }

    fn high_degree_bezier_surface() -> NurbsSurface<f64> {
        let mut knots = vec![0.0; 17];
        knots.extend([1.0; 17]);
        let controls = (0..17)
            .map(|u| {
                (0..17)
                    .map(|v| [f64::from(u) / 16.0, f64::from(v) / 16.0, 0.0, 1.0])
                    .collect()
            })
            .collect();
        let knots = KnotVector::new(knots, 16).expect("degree-16 surface knots");
        NurbsSurface::from_homogeneous(knots.clone(), knots, controls)
            .expect("degree-16 Bezier surface")
    }

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
        let curve = quadratic_curve();
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
        for _ in 0..2 {
            let owning_cancellable = with_closest_cx(false, |cx| {
                closest_point_curve_with_cx(&curve, query, 1e-9, 8, cx)
                    .expect("active cancellable owning closest")
            });
            assert_eq!(
                owning_cancellable,
                ClosestPointRun::Complete { estimate: first }
            );
        }
        let cancellable = with_closest_cx(false, |cx| {
            admitted
                .closest_point_with_cx(query, 1e-9, 8, cx)
                .expect("active cancellable closest")
        });
        assert_eq!(cancellable, ClosestPointRun::Complete { estimate: first });
    }

    #[test]
    fn curve_closest_cancellation_preserves_request_refusal_precedence() {
        let curve = quadratic_curve();
        let admitted = curve.admit().expect("admitted curve");
        let mut pre_scan_polls = 0usize;
        let pre_scan_refusal = preflight_curve_closest_with_poll(
            admitted,
            0,
            CLOSEST_MAX_BASE_WORK_UNITS,
            &mut || {
                pre_scan_polls += 1;
                true
            },
        )
        .expect_err("count-only pre-scan work refusal must precede cancellation");
        assert!(matches!(pre_scan_refusal, NurbsError::Domain { .. }));
        assert_eq!(pre_scan_polls, 0);

        with_closest_cx(true, |cx| {
            assert_eq!(
                closest_point_curve_with_cx(&curve, [1.0, 0.25, 0.0], 1e-9, 8, cx)
                    .expect("valid pre-cancelled owning request"),
                ClosestPointRun::Cancelled
            );
            assert_eq!(
                admitted
                    .closest_point_with_cx([1.0, 0.25, 0.0], 1e-9, 8, cx)
                    .expect("valid pre-cancelled request"),
                ClosestPointRun::Cancelled
            );
            for (query, tolerance, splits) in [
                ([f64::NAN, 0.0, 0.0], 0.0, 0),
                ([0.0; 3], -1.0, 0),
                ([0.0; 3], f64::NAN, 0),
                ([0.0; 3], 0.0, CLOSEST_MAX_SPLITS + 1),
            ] {
                let synchronous = admitted
                    .closest_point(query, tolerance, splits)
                    .expect_err("malformed synchronous request");
                let owning_synchronous = closest_point_curve(&curve, query, tolerance, splits)
                    .expect_err("malformed synchronous owning request");
                let cancellable = admitted
                    .closest_point_with_cx(query, tolerance, splits, cx)
                    .expect_err("request refusal must precede cancellation");
                assert_eq!(cancellable, synchronous);
                let owning_cancellable =
                    closest_point_curve_with_cx(&curve, query, tolerance, splits, cx)
                        .expect_err("owning request refusal must precede cancellation");
                assert_eq!(owning_cancellable, owning_synchronous);
            }
        });
    }

    #[test]
    fn owning_curve_closest_forwards_source_work_before_polling() {
        let curve = quadratic_curve();
        let admitted = curve.admit().expect("admitted curve");
        let expected = preflight_curve_closest(admitted, 0, CLOSEST_MAX_BASE_WORK_UNITS)
            .expect_err("source plus plan work must exceed the base ceiling");
        with_closest_cx(true, |cx| {
            assert_eq!(
                admitted
                    .closest_point_with_cx_after_request(
                        [1.0, 0.25, 0.0],
                        1e-9,
                        0,
                        CLOSEST_MAX_BASE_WORK_UNITS,
                        cx,
                    )
                    .expect_err("aggregate source-work refusal must precede cancellation"),
                expected
            );
        });
    }

    #[test]
    fn curve_closest_planning_hull_and_split_have_replayable_ordinals() {
        let long = long_linear_curve();
        let admitted = long.admit().expect("admitted long curve");
        let high_degree = high_degree_bezier_curve();
        let high_degree_admitted = high_degree.admit().expect("admitted degree-64 curve");
        let controls = high_degree.homogeneous_control_points();
        let seed_plan =
            preflight_curve_closest(high_degree_admitted, 0, 0).expect("high-degree seed plan");
        for _ in 0..2 {
            let mut plan_polls = 0usize;
            let plan = preflight_curve_closest_with_poll(admitted, 0, 0, &mut || {
                plan_polls += 1;
                plan_polls == 2
            })
            .expect("cancellable planning");
            assert!(matches!(plan, ClosestWorkRun::Cancelled));
            assert_eq!(plan_polls, 2);

            let mut hull_polls = 0usize;
            let hull = hull_lower_bound_with_poll([0.5, 1.0, 0.0], controls.iter(), &mut || {
                hull_polls += 1;
                hull_polls == 2
            });
            assert_eq!(hull, ClosestWorkRun::Cancelled);
            assert_eq!(hull_polls, 2);

            let mut split_polls = 0usize;
            let split = split_bezier_with_poll(controls, &mut || {
                split_polls += 1;
                split_polls == 7
            })
            .expect("cancellable split");
            assert_eq!(split, ClosestWorkRun::Cancelled);
            assert_eq!(split_polls, 7);

            let mut seed_kernels = CurveClosestKernels {
                conversion: || {
                    high_degree_admitted
                        .to_bezier_form()
                        .map(ClosestWorkRun::Complete)
                },
                derivatives: |t, order| {
                    high_degree_admitted
                        .derivatives(t, order)
                        .map(ClosestWorkRun::Complete)
                },
                evaluation: |t| high_degree_admitted.eval(t).map(ClosestWorkRun::Complete),
            };
            let mut seed_polls = 0usize;
            let seed = closest_point_curve_after_plan_with_poll(
                high_degree_admitted,
                [0.5, 1.0, 0.0],
                0.0,
                0,
                seed_plan,
                &mut seed_kernels,
                &mut || {
                    seed_polls += 1;
                    seed_polls == 7
                },
            )
            .expect("cancellable seed copy");
            assert_eq!(seed, ClosestWorkRun::Cancelled);
            assert_eq!(seed_polls, 7);
        }
    }

    #[test]
    fn curve_closest_kernel_cancellation_is_not_an_optional_polish_miss() {
        let curve = quadratic_curve();
        let admitted = curve.admit().expect("admitted curve");
        let plan = preflight_curve_closest(admitted, 0, 0).expect("closest plan");
        let query = [1.0, 0.25, 0.0];

        let mut conversion_cancel = CurveClosestKernels {
            conversion: || Ok(ClosestWorkRun::Cancelled),
            derivatives: |t, order| admitted.derivatives(t, order).map(ClosestWorkRun::Complete),
            evaluation: |t| admitted.eval(t).map(ClosestWorkRun::Complete),
        };
        let result = closest_point_curve_after_plan_with_poll(
            admitted,
            query,
            1e-9,
            0,
            plan,
            &mut conversion_cancel,
            &mut || false,
        )
        .expect("conversion cancellation");
        assert_eq!(result, ClosestWorkRun::Cancelled);

        let mut derivative_cancel = CurveClosestKernels {
            conversion: || admitted.to_bezier_form().map(ClosestWorkRun::Complete),
            derivatives: |_t, _order| Ok(ClosestWorkRun::Cancelled),
            evaluation: |t| admitted.eval(t).map(ClosestWorkRun::Complete),
        };
        let result = closest_point_curve_after_plan_with_poll(
            admitted,
            query,
            1e-9,
            0,
            plan,
            &mut derivative_cancel,
            &mut || false,
        )
        .expect("derivative cancellation");
        assert_eq!(result, ClosestWorkRun::Cancelled);

        let mut evaluation_cancel = CurveClosestKernels {
            conversion: || admitted.to_bezier_form().map(ClosestWorkRun::Complete),
            derivatives: |_t, _order| {
                Err(NurbsError::Domain {
                    what: "ordinary optional polish miss".to_string(),
                })
            },
            evaluation: |_t| Ok(ClosestWorkRun::Cancelled),
        };
        let result = closest_point_curve_after_plan_with_poll(
            admitted,
            query,
            1e-9,
            0,
            plan,
            &mut evaluation_cancel,
            &mut || false,
        )
        .expect("evaluation cancellation");
        assert_eq!(result, ClosestWorkRun::Cancelled);
    }

    #[test]
    fn curve_closest_final_publication_checkpoint_replays_exactly() {
        let curve = quadratic_curve();
        let admitted = curve.admit().expect("admitted curve");
        let plan = preflight_curve_closest(admitted, 2, 0).expect("closest plan");
        let query = [1.0, 0.25, 0.0];

        let mut first_kernels = CurveClosestKernels {
            conversion: || admitted.to_bezier_form().map(ClosestWorkRun::Complete),
            derivatives: |t, order| admitted.derivatives(t, order).map(ClosestWorkRun::Complete),
            evaluation: |t| admitted.eval(t).map(ClosestWorkRun::Complete),
        };
        let mut complete_polls = 0usize;
        let complete = closest_point_curve_after_plan_with_poll(
            admitted,
            query,
            1e-9,
            2,
            plan,
            &mut first_kernels,
            &mut || {
                complete_polls += 1;
                false
            },
        )
        .expect("complete closest replay");
        assert!(matches!(complete, ClosestWorkRun::Complete(_)));
        assert!(complete_polls > 0);

        let mut replay_kernels = CurveClosestKernels {
            conversion: || admitted.to_bezier_form().map(ClosestWorkRun::Complete),
            derivatives: |t, order| admitted.derivatives(t, order).map(ClosestWorkRun::Complete),
            evaluation: |t| admitted.eval(t).map(ClosestWorkRun::Complete),
        };
        let mut replay_polls = 0usize;
        let replay = closest_point_curve_after_plan_with_poll(
            admitted,
            query,
            1e-9,
            2,
            plan,
            &mut replay_kernels,
            &mut || {
                replay_polls += 1;
                replay_polls == complete_polls
            },
        )
        .expect("cancelled closest replay");
        assert_eq!(replay, ClosestWorkRun::Cancelled);
        assert_eq!(replay_polls, complete_polls);
    }

    #[test]
    fn surface_closest_owning_admitted_and_cx_paths_match_exactly() {
        let surface = multispan_quadratic_surface();
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
        for _ in 0..2 {
            let owning_cancellable = with_closest_cx(false, |cx| {
                closest_point_surface_with_cx(&surface, query, 1e-9, 8, cx)
                    .expect("active cancellable owning surface closest")
            });
            assert_eq!(
                owning_cancellable,
                ClosestPointRun::Complete { estimate: first }
            );
        }
        let cancellable = with_closest_cx(false, |cx| {
            admitted
                .closest_point_with_cx(query, 1e-9, 8, cx)
                .expect("active cancellable surface closest")
        });
        assert_eq!(cancellable, ClosestPointRun::Complete { estimate: first });
    }

    #[test]
    fn surface_closest_cancellation_preserves_request_refusal_precedence() {
        let surface = multispan_quadratic_surface();
        let admitted = surface.admit().expect("admitted surface");
        let mut pre_scan_polls = 0usize;
        let pre_scan_refusal = preflight_surface_closest_with_poll(
            admitted,
            0,
            CLOSEST_MAX_BASE_WORK_UNITS,
            &mut || {
                pre_scan_polls += 1;
                true
            },
        )
        .expect_err("surface count-only work refusal must precede cancellation");
        assert!(matches!(pre_scan_refusal, NurbsError::Domain { .. }));
        assert_eq!(pre_scan_polls, 0);

        with_closest_cx(true, |cx| {
            assert_eq!(
                closest_point_surface_with_cx(&surface, [0.25, 0.75, 0.5], 1e-9, 8, cx)
                    .expect("valid pre-cancelled owning surface request"),
                ClosestPointRun::Cancelled
            );
            assert_eq!(
                admitted
                    .closest_point_with_cx([0.25, 0.75, 0.5], 1e-9, 8, cx)
                    .expect("valid pre-cancelled surface request"),
                ClosestPointRun::Cancelled
            );
            for (query, tolerance, splits) in [
                ([f64::NAN, 0.0, 0.0], 0.0, 0),
                ([0.0; 3], -1.0, 0),
                ([0.0; 3], f64::NAN, 0),
                ([0.0; 3], 0.0, CLOSEST_MAX_SPLITS + 1),
            ] {
                let synchronous = admitted
                    .closest_point(query, tolerance, splits)
                    .expect_err("malformed synchronous surface request");
                let owning_synchronous = closest_point_surface(&surface, query, tolerance, splits)
                    .expect_err("malformed synchronous owning surface request");
                let cancellable = admitted
                    .closest_point_with_cx(query, tolerance, splits, cx)
                    .expect_err("surface request refusal must precede cancellation");
                assert_eq!(cancellable, synchronous);
                let owning_cancellable =
                    closest_point_surface_with_cx(&surface, query, tolerance, splits, cx)
                        .expect_err("owning surface request refusal must precede cancellation");
                assert_eq!(owning_cancellable, owning_synchronous);
            }
        });
    }

    #[test]
    fn owning_surface_closest_forwards_source_work_before_polling() {
        let surface = multispan_quadratic_surface();
        let admitted = surface.admit().expect("admitted surface");
        let expected = preflight_surface_closest(admitted, 0, CLOSEST_MAX_BASE_WORK_UNITS)
            .expect_err("surface source plus plan work must exceed the base ceiling");
        with_closest_cx(true, |cx| {
            assert_eq!(
                admitted
                    .closest_point_with_cx_after_request(
                        [0.25, 0.75, 0.5],
                        1e-9,
                        0,
                        CLOSEST_MAX_BASE_WORK_UNITS,
                        cx,
                    )
                    .expect_err("surface aggregate work refusal must precede cancellation"),
                expected
            );
        });
    }

    #[test]
    fn surface_closest_planning_conversion_seed_and_splits_have_replayable_ordinals() {
        let long = long_linear_surface();
        let long_admitted = long.admit().expect("admitted long surface");
        let multispan = multispan_quadratic_surface();
        let multispan_admitted = multispan.admit().expect("admitted multispan surface");
        let high_degree = high_degree_bezier_surface();
        for _ in 0..2 {
            let mut plan_polls = 0usize;
            let plan = preflight_surface_closest_with_poll(long_admitted, 0, 0, &mut || {
                plan_polls += 1;
                plan_polls == 2
            })
            .expect("cancellable surface planning");
            assert!(matches!(plan, ClosestWorkRun::Cancelled));
            assert_eq!(plan_polls, 2);

            let mut conversion_polls = 0usize;
            let conversion = to_bezier_surface_with_poll(multispan_admitted, &mut || {
                conversion_polls += 1;
                conversion_polls == 7
            })
            .expect("cancellable surface conversion");
            assert!(matches!(conversion, ClosestWorkRun::Cancelled));
            assert_eq!(conversion_polls, 7);

            let mut queue = BinaryHeap::new();
            let mut next_logical_id = 0u64;
            let mut upper = f64::INFINITY;
            let mut best = [0.0; 2];
            let mut seed_polls = 0usize;
            let seed = seed_patches_with_poll(
                &high_degree,
                [0.5, 0.5, 1.0],
                &mut queue,
                1,
                &mut next_logical_id,
                &mut upper,
                &mut best,
                &mut || {
                    seed_polls += 1;
                    seed_polls == 43
                },
            )
            .expect("cancellable surface seed copy");
            assert!(matches!(seed, ClosestWorkRun::Cancelled));
            assert_eq!(seed_polls, 43);

            let mut split_u_polls = 0usize;
            let split_u =
                split_patch_u_with_poll(high_degree.homogeneous_control_net(), &mut || {
                    split_u_polls += 1;
                    split_u_polls == 90
                })
                .expect("cancellable U patch split");
            assert!(matches!(split_u, ClosestWorkRun::Cancelled));
            assert_eq!(split_u_polls, 90);

            let mut split_v_polls = 0usize;
            let split_v =
                split_patch_v_with_poll(high_degree.homogeneous_control_net(), &mut || {
                    split_v_polls += 1;
                    split_v_polls == 8
                })
                .expect("cancellable V patch split");
            assert!(matches!(split_v, ClosestWorkRun::Cancelled));
            assert_eq!(split_v_polls, 8);
        }
    }

    #[test]
    fn surface_closest_kernel_cancellation_is_not_an_optional_evaluation_miss() {
        let surface = multispan_quadratic_surface();
        let admitted = surface.admit().expect("admitted surface");
        let plan = preflight_surface_closest(admitted, 0, 0).expect("surface closest plan");
        let query = [0.25, 0.75, 0.5];

        let mut conversion_cancel = SurfaceClosestKernels {
            conversion: || {
                Ok::<ClosestWorkRun<NurbsSurface<f64>>, NurbsError>(ClosestWorkRun::Cancelled)
            },
            evaluation: |u, v| admitted.eval(u, v).map(ClosestWorkRun::Complete),
        };
        let result = closest_point_surface_after_plan_with_poll(
            admitted,
            query,
            1e-9,
            0,
            plan,
            &mut conversion_cancel,
            &mut || false,
        )
        .expect("surface conversion cancellation");
        assert!(matches!(result, ClosestWorkRun::Cancelled));

        let mut evaluation_cancel = SurfaceClosestKernels {
            conversion: || to_bezier_surface(admitted).map(ClosestWorkRun::Complete),
            evaluation: |_u, _v| {
                Ok::<ClosestWorkRun<[f64; 3]>, NurbsError>(ClosestWorkRun::Cancelled)
            },
        };
        let result = closest_point_surface_after_plan_with_poll(
            admitted,
            query,
            1e-9,
            0,
            plan,
            &mut evaluation_cancel,
            &mut || false,
        )
        .expect("surface evaluation cancellation");
        assert!(matches!(result, ClosestWorkRun::Cancelled));

        let mut optional_error = SurfaceClosestKernels {
            conversion: || to_bezier_surface(admitted).map(ClosestWorkRun::Complete),
            evaluation: |_u, _v| {
                Err(NurbsError::Domain {
                    what: "ordinary optional surface evaluation miss".to_string(),
                })
            },
        };
        let result = closest_point_surface_after_plan_with_poll(
            admitted,
            query,
            1e-9,
            0,
            plan,
            &mut optional_error,
            &mut || false,
        )
        .expect("optional surface evaluation error");
        assert!(matches!(result, ClosestWorkRun::Complete(_)));
    }

    #[test]
    fn surface_closest_final_publication_checkpoint_replays_exactly() {
        let surface = multispan_quadratic_surface();
        let admitted = surface.admit().expect("admitted surface");
        let plan = preflight_surface_closest(admitted, 2, 0).expect("surface closest plan");
        let query = [0.25, 0.75, 0.5];

        let mut first_kernels = SurfaceClosestKernels {
            conversion: || to_bezier_surface(admitted).map(ClosestWorkRun::Complete),
            evaluation: |u, v| admitted.eval(u, v).map(ClosestWorkRun::Complete),
        };
        let mut complete_polls = 0usize;
        let complete = closest_point_surface_after_plan_with_poll(
            admitted,
            query,
            1e-9,
            2,
            plan,
            &mut first_kernels,
            &mut || {
                complete_polls += 1;
                false
            },
        )
        .expect("complete surface closest replay");
        assert!(matches!(complete, ClosestWorkRun::Complete(_)));
        assert!(complete_polls > 0);

        let mut replay_kernels = SurfaceClosestKernels {
            conversion: || to_bezier_surface(admitted).map(ClosestWorkRun::Complete),
            evaluation: |u, v| admitted.eval(u, v).map(ClosestWorkRun::Complete),
        };
        let mut replay_polls = 0usize;
        let replay = closest_point_surface_after_plan_with_poll(
            admitted,
            query,
            1e-9,
            2,
            plan,
            &mut replay_kernels,
            &mut || {
                replay_polls += 1;
                replay_polls == complete_polls
            },
        )
        .expect("cancelled surface closest replay");
        assert!(matches!(replay, ClosestWorkRun::Cancelled));
        assert_eq!(replay_polls, complete_polls);
    }

    #[test]
    fn surface_closest_preflight_composes_nested_insertion_refusal() {
        // A 500 x 700 source is cheap enough to admit, but its one required U
        // refinement exceeds the insertion engine's aggregate work ceiling.
        // Keep every other interior U knot at full degree multiplicity so the
        // closest-point conversion has exactly that one nested insertion.
        let mut knots_u = vec![0.0; 3];
        for index in 1..=248 {
            let knot = f64::from(index) / 250.0;
            knots_u.extend([knot, knot]);
        }
        let target = 249.0 / 250.0;
        knots_u.push(target);
        knots_u.extend([1.0; 3]);
        let knots_u = KnotVector::new(knots_u, 2).expect("large quadratic u knots");

        let mut knots_v = vec![0.0; 2];
        knots_v.extend((1..=698).map(|index| f64::from(index) / 699.0));
        knots_v.extend([1.0; 2]);
        let knots_v = KnotVector::new(knots_v, 1).expect("large linear v knots");
        let controls = vec![vec![[0.0, 0.0, 0.0, 1.0]; 700]; 500];
        let surface = NurbsSurface::from_homogeneous(knots_u, knots_v, controls)
            .expect("large admitted-refusal surface");
        let admitted = surface.admit().expect("admitted large surface");

        let direct = admitted
            .insert_knot_u(target)
            .expect_err("the nested insertion must refuse its work envelope");
        let closest = admitted
            .closest_point([0.0; 3], 1e-6, 0)
            .expect_err("closest preflight must compose the nested refusal");
        assert_eq!(closest, direct);
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
    fn surface_payload_model_accounts_for_rows_conversion_and_eval_overlap() {
        let expected_storage = (7 + 8) * size_of::<f64>()
            + 4 * size_of::<Vec<[f64; 4]>>()
            + 4 * 5 * size_of::<[f64; 4]>();
        assert_eq!(
            surface_storage_bytes(7, 8, 4, 5).expect("representable surface storage"),
            expected_storage as u128
        );

        let no_insertion_bytes = (4 + 4) * size_of::<f64>() as u128
            + 2 * size_of::<Vec<[f64; 4]>>() as u128
            + 2 * 2 * size_of::<[f64; 4]>() as u128;
        assert_eq!(
            surface_conversion_peak_allocated_bytes(0, 0, 4, 4, 2, 2)
                .expect("no-insertion conversion"),
            (no_insertion_bytes, no_insertion_bytes)
        );
        let converted_bytes = (5 + 5) * size_of::<f64>() as u128
            + 3 * size_of::<Vec<[f64; 4]>>() as u128
            + 3 * 3 * size_of::<[f64; 4]>() as u128;
        let previous_bytes = (5 + 4) * size_of::<f64>() as u128
            + 3 * size_of::<Vec<[f64; 4]>>() as u128
            + 3 * 2 * size_of::<[f64; 4]>() as u128;
        assert_eq!(
            surface_conversion_peak_allocated_bytes(1, 1, 5, 5, 3, 3)
                .expect("two-direction conversion"),
            (converted_bytes, converted_bytes + previous_bytes)
        );

        assert_eq!(
            surface_final_eval_workspace_bytes(8, 1).expect("u-dominant workspace"),
            24 * size_of::<f64>() as u128
        );
        assert_eq!(
            surface_final_eval_workspace_bytes(1, 8).expect("v-dominant workspace"),
            25 * size_of::<f64>() as u128
        );
    }

    #[test]
    fn surface_retained_envelope_composes_every_live_phase() {
        assert!(
            enforce_surface_retained_envelope(1, CLOSEST_MAX_RETAINED_BYTES - 1, 0, 0, 0,).is_ok(),
            "source plus conversion allocation may reach the exact cap"
        );
        let conversion_error =
            enforce_surface_retained_envelope(1, CLOSEST_MAX_RETAINED_BYTES, 0, 0, 0)
                .expect_err("conversion residency one byte above the cap must refuse");
        assert!(matches!(conversion_error, crate::NurbsError::Domain { .. }));

        assert!(
            enforce_surface_retained_envelope(1, 0, 2, CLOSEST_MAX_RETAINED_BYTES - 4, 1,).is_ok(),
            "source, converted surface, frontier, and eval workspace compose at the exact cap"
        );
        let eval_error =
            enforce_surface_retained_envelope(1, 0, 2, CLOSEST_MAX_RETAINED_BYTES - 4, 2)
                .expect_err("aggregate final evaluation one byte above the cap must refuse");
        assert!(matches!(eval_error, crate::NurbsError::Domain { .. }));

        for overflow in [
            enforce_surface_retained_envelope(1, u128::MAX, 0, 0, 0),
            enforce_surface_retained_envelope(1, 0, u128::MAX, 0, 0),
            enforce_surface_retained_envelope(0, 0, 0, u128::MAX, 1),
        ] {
            assert!(
                matches!(overflow, Err(crate::NurbsError::Domain { .. })),
                "aggregate retained-byte overflow must refuse"
            );
        }
    }

    #[test]
    fn surface_payload_model_refuses_arithmetic_overflow() {
        for overflow in [
            surface_storage_bytes(u128::MAX, 1, 1, 1),
            surface_storage_bytes(1, 1, u128::MAX, 2),
        ] {
            assert!(
                matches!(overflow, Err(crate::NurbsError::Domain { .. })),
                "surface storage overflow must refuse"
            );
        }
        assert!(matches!(
            surface_conversion_peak_allocated_bytes(1, 1, u128::MAX, 1, 1, 1),
            Err(crate::NurbsError::Domain { .. })
        ));
        let late_peak_overflow_knot_count = u128::MAX / (2 * size_of::<f64>() as u128) + 2;
        assert!(matches!(
            surface_conversion_peak_allocated_bytes(1, 0, late_peak_overflow_knot_count, 2, 2, 1,),
            Err(crate::NurbsError::Domain { .. })
        ));
        assert!(matches!(
            surface_final_eval_workspace_bytes(u128::MAX, 1),
            Err(crate::NurbsError::Domain { .. })
        ));
        assert!(matches!(
            surface_final_eval_workspace_bytes(1, u128::MAX),
            Err(crate::NurbsError::Domain { .. })
        ));
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
