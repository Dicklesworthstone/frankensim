//! Rational B-spline curves over a generic scalar: homogeneous de Boor
//! evaluation, derivatives to arbitrary order (f64 path), EXACT Boehm
//! knot insertion, Bézier decomposition, and EXACT degree elevation via
//! per-segment Bézier elevation (the elevated curve carries a
//! full-multiplicity knot vector — valid, evaluation-identical; minimal
//! knot vectors are a documented follow-up).

use crate::NurbsError;
use crate::basis::{
    AdmittedKnotVector, BASIS_MAX_WORK_UNITS, BasisRun, KnotSpanRun, KnotValidationOutcome,
    KnotVector, Scalar,
};
use fs_exec::Cx;

// Conservative price for finite/weight/projection/canonical-lane validation of
// one homogeneous control. Structural admission must precede the full scan.
const CURVE_VALIDATION_WORK_PER_CONTROL: u128 = 16;
// Keep this conservative price aligned with KnotVector's private validation
// envelope: every derived insertion constructs and then revalidates its knots.
const CURVE_KNOT_VALIDATION_WORK_PER_ENTRY: u128 = 16;
const CURVE_BEZIER_SCAN_WORK_PER_KNOT: u128 = 4;
const CURVE_BEZIER_BLEND_WORK_PER_CONTROL: u128 = 32;
const CURVE_SPAN_BOX_WORK_PER_CONTROL: u128 = 16;
const CURVE_SPAN_BOX_MAX_RETAINED_BYTES: u128 = 64 * 1024 * 1024;
const CURVE_CONSTRUCTION_MAX_RETAINED_BYTES: u128 = 64 * 1024 * 1024;
const CURVE_COPY_MAX_RETAINED_BYTES: u128 = 64 * 1024 * 1024;
const CURVE_INSERTION_MAX_RETAINED_BYTES: u128 = 64 * 1024 * 1024;
const CURVE_REMOVAL_MAX_DERIVED_BYTES: u128 = 64 * 1024 * 1024;
const CURVE_BEZIER_MAX_RETAINED_BYTES: u128 = 64 * 1024 * 1024;
const CURVE_ELEVATION_MAX_RETAINED_BYTES: u128 = 64 * 1024 * 1024;
const CURVE_CANCELLATION_STRIDE: usize = 64;

fn curve_poll_due(
    operations_since_poll: &mut usize,
    should_cancel: &mut impl FnMut() -> bool,
) -> bool {
    *operations_since_poll += 1;
    if *operations_since_poll < CURVE_CANCELLATION_STRIDE {
        return false;
    }
    *operations_since_poll = 0;
    should_cancel()
}

#[derive(Debug)]
enum CurveWorkRun<T> {
    Complete(T),
    Cancelled,
}

/// Checked shape/work/retained-memory plan for exact Bezier conversion.
/// Fields are crate-visible so trim/closest primitives can compose this
/// conversion phase with their own simultaneously-live scratch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct BezierConversionPlan {
    pub(crate) insertions: usize,
    pub(crate) distinct_knot_count: usize,
    pub(crate) final_knot_count: usize,
    pub(crate) final_control_count: usize,
    pub(crate) work_units: u128,
    /// Bytes allocated by the conversion at its peak, excluding the borrowed
    /// source generation. With insertions this is the current derived curve
    /// plus the next derived curve under construction.
    pub(crate) peak_allocated_bytes: u128,
    /// Retained bytes in the returned converted curve.
    pub(crate) converted_bytes: u128,
}

/// Checked shape/work/retained-memory plan for exact degree elevation.
///
/// The borrowed source generation is excluded from retained-byte accounting.
/// The plan includes the simultaneously-live Bezier generation, elevation
/// metadata, and final knot/control payloads.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CurveElevationPlan {
    distinct_knot_count: usize,
    segment_count: usize,
    elevated_degree: usize,
    elevated_order: usize,
    elevated_degree_i64: i64,
    final_knot_count: usize,
    final_control_count: usize,
    work_units: u128,
    peak_retained_bytes: u128,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CurveInsertionPlan {
    new_knot_count: usize,
    new_control_count: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CurveRemovalPlan {
    new_knot_count: usize,
    new_control_count: usize,
    reconstruction_capacity: usize,
    verification_insertion: CurveInsertionPlan,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct CurveRemovalRange {
    left_control: usize,
    blend_start: usize,
    blend_end: usize,
    suffix_start: usize,
    forward_count: usize,
}

fn curve_storage_bytes<S: Scalar>(
    knot_count: usize,
    control_count: usize,
) -> Result<u128, NurbsError> {
    let knot_bytes = (knot_count as u128)
        .checked_mul(core::mem::size_of::<S>() as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "Bezier knot-storage accounting overflows u128".to_string(),
        })?;
    let control_bytes = (control_count as u128)
        .checked_mul(core::mem::size_of::<[S; 4]>() as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "Bezier control-storage accounting overflows u128".to_string(),
        })?;
    knot_bytes
        .checked_add(control_bytes)
        .ok_or_else(|| NurbsError::Domain {
            what: "Bezier curve-storage accounting overflows u128".to_string(),
        })
}

fn enforce_curve_elevation_envelope(
    work_units: u128,
    peak_retained_bytes: u128,
) -> Result<(), NurbsError> {
    if work_units > BASIS_MAX_WORK_UNITS {
        return Err(NurbsError::Domain {
            what: format!(
                "curve degree elevation requests {work_units} work units above defensive ceiling {BASIS_MAX_WORK_UNITS}"
            ),
        });
    }
    if peak_retained_bytes > CURVE_ELEVATION_MAX_RETAINED_BYTES {
        return Err(NurbsError::Domain {
            what: format!(
                "curve degree elevation can retain {peak_retained_bytes} derived payload bytes above defensive ceiling {CURVE_ELEVATION_MAX_RETAINED_BYTES}"
            ),
        });
    }
    Ok(())
}

fn plan_curve_elevation<S: Scalar>(
    degree: usize,
    bezier: BezierConversionPlan,
) -> Result<CurveElevationPlan, NurbsError> {
    if bezier.distinct_knot_count < 2 {
        return Err(NurbsError::Structure {
            what: "degree elevation requires at least two distinct knots".to_string(),
        });
    }
    let segment_count =
        bezier
            .distinct_knot_count
            .checked_sub(1)
            .ok_or_else(|| NurbsError::Structure {
                what: "degree-elevation segment count underflows usize".to_string(),
            })?;
    let elevated_degree = degree.checked_add(1).ok_or_else(|| NurbsError::Structure {
        what: "degree elevation overflows spline-degree arithmetic".to_string(),
    })?;
    let elevated_order = degree.checked_add(2).ok_or_else(|| NurbsError::Structure {
        what: "degree elevation overflows spline-order arithmetic".to_string(),
    })?;
    let elevated_degree_i64 =
        i64::try_from(elevated_degree).map_err(|_| NurbsError::Structure {
            what: "degree elevation exceeds the scalar integer-lift domain".to_string(),
        })?;
    let final_knot_count = bezier
        .final_knot_count
        .checked_add(bezier.distinct_knot_count)
        .ok_or_else(|| NurbsError::Domain {
            what: "degree-elevation final knot count overflows usize".to_string(),
        })?;
    let final_control_count = bezier
        .final_control_count
        .checked_add(segment_count)
        .ok_or_else(|| NurbsError::Domain {
            what: "degree-elevation final control count overflows usize".to_string(),
        })?;

    // Conservative aggregate envelope for the already-priced Bezier
    // conversion, post-conversion run/span scans, four-lane binomial blends,
    // output assembly, and both derived structural validation passes.
    let knot_run_scan = (bezier.final_knot_count as u128)
        .checked_mul(CURVE_BEZIER_SCAN_WORK_PER_KNOT)
        .ok_or_else(|| NurbsError::Domain {
            what: "degree-elevation knot-run scan work overflows u128".to_string(),
        })?;
    let span_scan = (bezier.final_control_count as u128)
        .checked_mul(4)
        .ok_or_else(|| NurbsError::Domain {
            what: "degree-elevation span scan work overflows u128".to_string(),
        })?;
    let blend_work = (degree as u128)
        .checked_mul(segment_count as u128)
        .and_then(|work| work.checked_mul(CURVE_BEZIER_BLEND_WORK_PER_CONTROL))
        .ok_or_else(|| NurbsError::Domain {
            what: "degree-elevation blend work overflows u128".to_string(),
        })?;
    let control_assembly = (final_control_count as u128)
        .checked_mul(4)
        .ok_or_else(|| NurbsError::Domain {
            what: "degree-elevation control assembly work overflows u128".to_string(),
        })?;
    let knot_assembly = final_knot_count as u128;
    let one_knot_validation = (final_knot_count as u128)
        .checked_mul(CURVE_KNOT_VALIDATION_WORK_PER_ENTRY)
        .and_then(|work| work.checked_add(elevated_degree as u128))
        .ok_or_else(|| NurbsError::Domain {
            what: "degree-elevation knot-validation work overflows u128".to_string(),
        })?;
    let derived_validation = one_knot_validation
        .checked_mul(2)
        .and_then(|work| {
            work.checked_add(
                (final_control_count as u128).checked_mul(CURVE_VALIDATION_WORK_PER_CONTROL)?,
            )
        })
        .ok_or_else(|| NurbsError::Domain {
            what: "degree-elevation derived-validation work overflows u128".to_string(),
        })?;
    let work_units = bezier
        .work_units
        .checked_add(knot_run_scan)
        .and_then(|work| work.checked_add(span_scan))
        .and_then(|work| work.checked_add(blend_work))
        .and_then(|work| work.checked_add(control_assembly))
        .and_then(|work| work.checked_add(knot_assembly))
        .and_then(|work| work.checked_add(derived_validation))
        .and_then(|work| work.checked_add(32))
        .ok_or_else(|| NurbsError::Domain {
            what: "degree-elevation aggregate work overflows u128".to_string(),
        })?;

    // Preserve deterministic refusal precedence: aggregate work is rejected
    // before retained-storage arithmetic or its ceiling is considered.
    if work_units > BASIS_MAX_WORK_UNITS {
        return Err(NurbsError::Domain {
            what: format!(
                "curve degree elevation requests {work_units} work units above defensive ceiling {BASIS_MAX_WORK_UNITS}"
            ),
        });
    }

    let break_bytes = (bezier.distinct_knot_count as u128)
        .checked_mul(core::mem::size_of::<S>() as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "degree-elevation break-storage accounting overflows u128".to_string(),
        })?;
    let multiplicity_bytes = (bezier.distinct_knot_count as u128)
        .checked_mul(core::mem::size_of::<usize>() as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "degree-elevation multiplicity-storage accounting overflows u128".to_string(),
        })?;
    let elevated_bytes = curve_storage_bytes::<S>(final_knot_count, final_control_count)?;
    let assembly_bytes = bezier
        .converted_bytes
        .checked_add(break_bytes)
        .and_then(|bytes| bytes.checked_add(multiplicity_bytes))
        .and_then(|bytes| bytes.checked_add(elevated_bytes))
        .ok_or_else(|| NurbsError::Domain {
            what: "degree-elevation retained-byte accounting overflows u128".to_string(),
        })?;
    let peak_retained_bytes = bezier.peak_allocated_bytes.max(assembly_bytes);
    let plan = CurveElevationPlan {
        distinct_knot_count: bezier.distinct_knot_count,
        segment_count,
        elevated_degree,
        elevated_order,
        elevated_degree_i64,
        final_knot_count,
        final_control_count,
        work_units,
        peak_retained_bytes,
    };
    enforce_curve_elevation_envelope(plan.work_units, plan.peak_retained_bytes)?;
    Ok(plan)
}

fn push_curve_elevation_value<T>(
    values: &mut Vec<T>,
    value: T,
    planned_len: usize,
    payload: &'static str,
) -> Result<(), NurbsError> {
    if values.len() >= planned_len {
        return Err(NurbsError::Structure {
            what: format!("degree-elevation {payload} exceeded its checked plan"),
        });
    }
    values.push(value);
    Ok(())
}

fn preflight_curve_copy<S: Scalar>(
    knot_count: usize,
    control_count: usize,
) -> Result<(), NurbsError> {
    let work_units = (control_count as u128)
        .checked_mul(4)
        .and_then(|work| work.checked_add(knot_count as u128))
        .and_then(|work| work.checked_add(2))
        .ok_or_else(|| NurbsError::Domain {
            what: "curve-copy work accounting overflows u128".to_string(),
        })?;
    if work_units > BASIS_MAX_WORK_UNITS {
        return Err(NurbsError::Domain {
            what: format!(
                "curve copy requests {work_units} work units above defensive ceiling {BASIS_MAX_WORK_UNITS}"
            ),
        });
    }
    let retained_bytes = (knot_count as u128)
        .checked_mul(core::mem::size_of::<S>() as u128)
        .and_then(|bytes| {
            bytes.checked_add(
                (control_count as u128).checked_mul(core::mem::size_of::<[S; 4]>() as u128)?,
            )
        })
        .ok_or_else(|| NurbsError::Domain {
            what: "curve-copy retained-byte accounting overflows u128".to_string(),
        })?;
    if retained_bytes > CURVE_COPY_MAX_RETAINED_BYTES {
        return Err(NurbsError::Domain {
            what: format!(
                "curve copy retains {retained_bytes} output bytes above defensive ceiling {CURVE_COPY_MAX_RETAINED_BYTES}"
            ),
        });
    }
    Ok(())
}

fn preflight_cartesian_curve_construction<S: Scalar>(
    control_count: usize,
) -> Result<(), NurbsError> {
    let retained_bytes = (control_count as u128)
        .checked_mul(core::mem::size_of::<[S; 4]>() as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "Cartesian curve control-storage accounting overflows u128".to_string(),
        })?;
    if retained_bytes > CURVE_CONSTRUCTION_MAX_RETAINED_BYTES {
        return Err(NurbsError::Domain {
            what: format!(
                "Cartesian curve retains {retained_bytes} homogeneous-control payload bytes above defensive ceiling {CURVE_CONSTRUCTION_MAX_RETAINED_BYTES}"
            ),
        });
    }
    Ok(())
}

fn validate_cartesian_curve_inputs_with_poll<S: Scalar, const DIM: usize>(
    points: &[[S; DIM]],
    weights: &[S],
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<CurveWorkRun<()>, NurbsError> {
    let mut operations_since_poll = 0usize;
    for &weight in weights {
        if !weight.is_admissible_weight() {
            return Err(NurbsError::Structure {
                what: "weights must be finite, positive, and numerically admissible".to_string(),
            });
        }
        if curve_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(CurveWorkRun::Cancelled);
        }
    }
    if should_cancel() {
        return Ok(CurveWorkRun::Cancelled);
    }

    operations_since_poll = 0;
    for point in points {
        if DIM == 0 && curve_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(CurveWorkRun::Cancelled);
        }
        for &coordinate in point {
            if !coordinate.is_finite() {
                return Err(NurbsError::Structure {
                    what: "control-point coordinates must be finite".to_string(),
                });
            }
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveWorkRun::Cancelled);
            }
        }
    }
    if should_cancel() {
        return Ok(CurveWorkRun::Cancelled);
    }
    Ok(CurveWorkRun::Complete(()))
}

fn build_cartesian_curve_controls_with_poll<S: Scalar, const DIM: usize>(
    points: &[[S; DIM]],
    weights: &[S],
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<CurveWorkRun<Vec<[S; 4]>>, NurbsError> {
    let mut cpw = Vec::new();
    cpw.try_reserve_exact(points.len())
        .map_err(|_| NurbsError::Domain {
            what: "curve homogeneous-control allocation was refused".to_string(),
        })?;
    if should_cancel() {
        return Ok(CurveWorkRun::Cancelled);
    }

    let mut operations_since_poll = 0usize;
    for (point, &weight) in points.iter().zip(weights) {
        let mut homogeneous = [S::zero(); 4];
        for (slot, &coordinate) in homogeneous.iter_mut().zip(point.iter()) {
            *slot = coordinate * weight;
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveWorkRun::Cancelled);
            }
        }
        homogeneous[3] = weight;
        cpw.push(homogeneous);
        if curve_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(CurveWorkRun::Cancelled);
        }
    }
    if should_cancel() {
        return Ok(CurveWorkRun::Cancelled);
    }
    Ok(CurveWorkRun::Complete(cpw))
}

fn validate_cartesian_curve_products_with_poll<S: Scalar, const DIM: usize>(
    points: &[[S; DIM]],
    cpw: &[[S; 4]],
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<CurveWorkRun<()>, NurbsError> {
    let mut operations_since_poll = 0usize;
    for (point, homogeneous) in points.iter().zip(cpw) {
        if DIM == 0 && curve_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(CurveWorkRun::Cancelled);
        }
        for (&coordinate, &weighted) in point.iter().zip(homogeneous) {
            if coordinate != S::zero() && weighted == S::zero() {
                return Err(NurbsError::Structure {
                    what: "Cartesian coordinate × weight underflowed a nonzero coordinate to zero"
                        .to_string(),
                });
            }
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveWorkRun::Cancelled);
            }
        }
    }
    if should_cancel() {
        return Ok(CurveWorkRun::Cancelled);
    }

    operations_since_poll = 0;
    for homogeneous in cpw {
        for &component in homogeneous {
            if !component.is_finite() {
                return Err(NurbsError::Structure {
                    what: "Cartesian coordinate × weight overflowed the homogeneous numeric domain"
                        .to_string(),
                });
            }
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveWorkRun::Cancelled);
            }
        }
    }
    if should_cancel() {
        return Ok(CurveWorkRun::Cancelled);
    }
    Ok(CurveWorkRun::Complete(()))
}

fn refinement_work_upper_bound(
    degree: usize,
    initial_knot_count: usize,
    initial_control_count: usize,
    direct_insertions: usize,
    bezier_insertions: usize,
) -> Result<u128, NurbsError> {
    let total_insertions = direct_insertions
        .checked_add(bezier_insertions)
        .ok_or_else(|| NurbsError::Domain {
            what: "curve refinement insertion count overflows usize".to_string(),
        })?;
    let final_knot_count = initial_knot_count
        .checked_add(total_insertions)
        .ok_or_else(|| NurbsError::Domain {
            what: "curve refinement knot count overflows usize".to_string(),
        })?;
    let final_control_count = initial_control_count
        .checked_add(total_insertions)
        .ok_or_else(|| NurbsError::Domain {
            what: "curve refinement control count overflows usize".to_string(),
        })?;

    // Conversion scans once to construct its plan and once per target-search
    // pass, including the final no-target pass. Reserve one additional scan so
    // callers such as trim can inspect the plan before invoking conversion
    // without leaving that preflight work uncharged.
    let scan_passes =
        (bezier_insertions as u128)
            .checked_add(3)
            .ok_or_else(|| NurbsError::Domain {
                what: "Bezier target-scan pass count overflows u128".to_string(),
            })?;
    let scan_work = (final_knot_count as u128)
        .checked_mul(CURVE_BEZIER_SCAN_WORK_PER_KNOT)
        .and_then(|work| work.checked_mul(scan_passes))
        .ok_or_else(|| NurbsError::Domain {
            what: "Bezier target-scan work overflows u128".to_string(),
        })?;
    let clone_work = (final_control_count as u128)
        .checked_mul(4)
        .and_then(|work| work.checked_add(final_knot_count as u128))
        .ok_or_else(|| NurbsError::Domain {
            what: "Bezier source-copy work overflows u128".to_string(),
        })?;

    // Each insertion performs a span search, copies both output arrays,
    // blends at most `degree` homogeneous controls, validates the constructed
    // KnotVector, then validates the published NurbsCurve (which scans those
    // knots a second time). Price every generation at the final, largest
    // shape so the aggregate bound is monotone and conservative.
    let knot_validation = (final_knot_count as u128)
        .checked_mul(CURVE_KNOT_VALIDATION_WORK_PER_ENTRY)
        .and_then(|work| work.checked_add(degree as u128))
        .ok_or_else(|| NurbsError::Domain {
            what: "Bezier knot-validation work overflows u128".to_string(),
        })?;
    let curve_validation = (final_control_count as u128)
        .checked_mul(CURVE_VALIDATION_WORK_PER_CONTROL)
        .and_then(|work| work.checked_add(knot_validation))
        .ok_or_else(|| NurbsError::Domain {
            what: "Bezier curve-validation work overflows u128".to_string(),
        })?;
    let copy_work = (final_control_count as u128)
        .checked_mul(4)
        .and_then(|work| work.checked_add(final_knot_count as u128))
        .ok_or_else(|| NurbsError::Domain {
            what: "Bezier insertion-copy work overflows u128".to_string(),
        })?;
    let blend_work = (degree as u128)
        .checked_mul(CURVE_BEZIER_BLEND_WORK_PER_CONTROL)
        .ok_or_else(|| NurbsError::Domain {
            what: "Bezier blend work overflows u128".to_string(),
        })?;
    let per_insertion = (final_control_count as u128)
        .checked_add(copy_work)
        .and_then(|work| work.checked_add(blend_work))
        .and_then(|work| work.checked_add(knot_validation))
        .and_then(|work| work.checked_add(curve_validation))
        .and_then(|work| work.checked_add(32))
        .ok_or_else(|| NurbsError::Domain {
            what: "Bezier per-insertion work overflows u128".to_string(),
        })?;
    let insertion_work = (total_insertions as u128)
        .checked_mul(per_insertion)
        .ok_or_else(|| NurbsError::Domain {
            what: "Bezier aggregate insertion work overflows u128".to_string(),
        })?;

    scan_work
        .checked_add(clone_work)
        .and_then(|work| work.checked_add(insertion_work))
        .ok_or_else(|| NurbsError::Domain {
            what: "Bezier aggregate refinement work overflows u128".to_string(),
        })
}

fn enforce_curve_insertion_envelope(
    work_units: u128,
    retained_bytes: u128,
) -> Result<(), NurbsError> {
    if work_units > BASIS_MAX_WORK_UNITS {
        return Err(NurbsError::Domain {
            what: format!(
                "curve knot insertion requests {work_units} work units above defensive ceiling {BASIS_MAX_WORK_UNITS}"
            ),
        });
    }
    if retained_bytes > CURVE_INSERTION_MAX_RETAINED_BYTES {
        return Err(NurbsError::Domain {
            what: format!(
                "curve knot insertion retains {retained_bytes} output bytes above defensive ceiling {CURVE_INSERTION_MAX_RETAINED_BYTES}"
            ),
        });
    }
    Ok(())
}

fn plan_curve_insertion<S: Scalar, const DIM: usize>(
    curve: AdmittedNurbsCurve<'_, S, DIM>,
) -> Result<CurveInsertionPlan, NurbsError> {
    let new_knot_count =
        curve
            .knots()
            .knots()
            .len()
            .checked_add(1)
            .ok_or_else(|| NurbsError::Domain {
                what: "inserted knot count overflows usize".to_string(),
            })?;
    let new_control_count = curve
        .homogeneous_control_points()
        .len()
        .checked_add(1)
        .ok_or_else(|| NurbsError::Domain {
            what: "inserted control count overflows usize".to_string(),
        })?;
    let work_units = refinement_work_upper_bound(
        curve.knots().degree(),
        curve.knots().knots().len(),
        curve.homogeneous_control_points().len(),
        1,
        0,
    )?;
    let retained_bytes = curve_storage_bytes::<S>(new_knot_count, new_control_count)?;
    enforce_curve_insertion_envelope(work_units, retained_bytes)?;
    Ok(CurveInsertionPlan {
        new_knot_count,
        new_control_count,
    })
}

fn curve_removal_work_units(
    degree: usize,
    knot_count: usize,
    control_count: usize,
    new_knot_count: usize,
    new_control_count: usize,
    reconstruction_capacity: usize,
) -> Result<u128, NurbsError> {
    // Cover occurrence discovery, knot/control copies, the exact forward
    // recurrence, both derived validation passes, full representation
    // comparison, and the restoring insertion at their largest shapes.
    let knot_work = (knot_count as u128)
        .checked_mul(16)
        .ok_or_else(|| NurbsError::Domain {
            what: "knot-removal knot work overflows u128".to_string(),
        })?;
    let control_work =
        (control_count as u128)
            .checked_mul(48)
            .ok_or_else(|| NurbsError::Domain {
                what: "knot-removal control work overflows u128".to_string(),
            })?;
    let reconstruction_work = (reconstruction_capacity as u128)
        .checked_mul(CURVE_BEZIER_BLEND_WORK_PER_CONTROL)
        .ok_or_else(|| NurbsError::Domain {
            what: "knot-removal reconstruction work overflows u128".to_string(),
        })?;
    let derived_validation_work = (new_knot_count as u128)
        .checked_mul(CURVE_KNOT_VALIDATION_WORK_PER_ENTRY)
        .and_then(|work| work.checked_mul(2))
        .and_then(|work| {
            work.checked_add(
                (new_control_count as u128).checked_mul(CURVE_VALIDATION_WORK_PER_CONTROL)?,
            )
        })
        .ok_or_else(|| NurbsError::Domain {
            what: "knot-removal derived-validation work overflows u128".to_string(),
        })?;
    let verification_work =
        refinement_work_upper_bound(degree, new_knot_count, new_control_count, 1, 0)?;
    knot_work
        .checked_add(control_work)
        .and_then(|work| work.checked_add(reconstruction_work))
        .and_then(|work| work.checked_add(derived_validation_work))
        .and_then(|work| work.checked_add(verification_work))
        .and_then(|work| work.checked_add(64))
        .ok_or_else(|| NurbsError::Domain {
            what: "knot-removal aggregate work overflows u128".to_string(),
        })
}

fn enforce_curve_removal_envelope<S: Scalar>(
    knot_count: usize,
    control_count: usize,
    new_knot_count: usize,
    new_control_count: usize,
    work_units: u128,
) -> Result<(), NurbsError> {
    if work_units > BASIS_MAX_WORK_UNITS {
        return Err(NurbsError::Domain {
            what: format!(
                "curve knot removal requests {work_units} work units above defensive ceiling {BASIS_MAX_WORK_UNITS}"
            ),
        });
    }

    // The removed candidate remains live while exact reinsertion constructs
    // the restored verifier generation. Source storage is borrowed.
    let derived_bytes = curve_storage_bytes::<S>(new_knot_count, new_control_count)?
        .checked_add(curve_storage_bytes::<S>(knot_count, control_count)?)
        .ok_or_else(|| NurbsError::Domain {
            what: "knot-removal aggregate derived-byte accounting overflows u128".to_string(),
        })?;
    if derived_bytes > CURVE_REMOVAL_MAX_DERIVED_BYTES {
        return Err(NurbsError::Domain {
            what: format!(
                "curve knot removal requires {derived_bytes} simultaneously-live derived bytes above defensive ceiling {CURVE_REMOVAL_MAX_DERIVED_BYTES}"
            ),
        });
    }
    Ok(())
}

fn plan_curve_removal<S: Scalar, const DIM: usize>(
    curve: AdmittedNurbsCurve<'_, S, DIM>,
) -> Result<CurveRemovalPlan, NurbsError> {
    let knot_count = curve.knots().knots().len();
    let control_count = curve.homogeneous_control_points().len();
    let new_knot_count = knot_count
        .checked_sub(1)
        .ok_or_else(|| NurbsError::Structure {
            what: "knot removal requires at least one source knot".to_string(),
        })?;
    let new_control_count = control_count
        .checked_sub(1)
        .ok_or_else(|| NurbsError::Structure {
            what: "knot removal requires at least one source control".to_string(),
        })?;
    let degree = curve.knots().degree();
    let reconstruction_capacity = degree.checked_add(1).ok_or_else(|| NurbsError::Structure {
        what: "knot-removal reconstruction capacity overflows usize".to_string(),
    })?;
    let work_units = curve_removal_work_units(
        degree,
        knot_count,
        control_count,
        new_knot_count,
        new_control_count,
        reconstruction_capacity,
    )?;
    enforce_curve_removal_envelope::<S>(
        knot_count,
        control_count,
        new_knot_count,
        new_control_count,
        work_units,
    )?;
    Ok(CurveRemovalPlan {
        new_knot_count,
        new_control_count,
        reconstruction_capacity,
        verification_insertion: CurveInsertionPlan {
            new_knot_count: knot_count,
            new_control_count: control_count,
        },
    })
}

fn curve_removal_range(
    degree: usize,
    removed_index: usize,
    prior_multiplicity: usize,
    reconstruction_capacity: usize,
) -> Result<CurveRemovalRange, NurbsError> {
    let span = removed_index
        .checked_sub(1)
        .ok_or_else(|| NurbsError::Structure {
            what: "interior knot removal has no left span".to_string(),
        })?;
    let left_control = span
        .checked_sub(degree)
        .ok_or_else(|| NurbsError::Structure {
            what: "knot-removal span is smaller than the degree".to_string(),
        })?;
    let blend_start = left_control
        .checked_add(1)
        .ok_or_else(|| NurbsError::Structure {
            what: "knot-removal blend start overflows usize".to_string(),
        })?;
    let blend_end = span
        .checked_sub(prior_multiplicity)
        .ok_or_else(|| NurbsError::Structure {
            what: "knot-removal multiplicity exceeds its span index".to_string(),
        })?;
    let suffix_start = blend_end
        .checked_add(1)
        .ok_or_else(|| NurbsError::Structure {
            what: "knot-removal suffix index overflows usize".to_string(),
        })?;
    let forward_count = blend_end
        .checked_sub(left_control)
        .and_then(|count| count.checked_add(1))
        .ok_or_else(|| NurbsError::Structure {
            what: "knot-removal reconstruction range is invalid".to_string(),
        })?;
    if forward_count > reconstruction_capacity {
        return Err(NurbsError::Structure {
            what: "knot-removal reconstruction exceeds its admitted capacity".to_string(),
        });
    }
    Ok(CurveRemovalRange {
        left_control,
        blend_start,
        blend_end,
        suffix_start,
        forward_count,
    })
}

fn bezier_pre_scan_work(knot_count: usize) -> Result<u128, NurbsError> {
    (knot_count as u128)
        .checked_mul(CURVE_BEZIER_SCAN_WORK_PER_KNOT)
        .ok_or_else(|| NurbsError::Domain {
            what: "Bezier pre-scan work overflows u128".to_string(),
        })
}

#[cfg(test)]
fn plan_bezier_conversion<S: Scalar, const DIM: usize>(
    curve: AdmittedNurbsCurve<'_, S, DIM>,
) -> Result<BezierConversionPlan, NurbsError> {
    let mut never_cancel = || false;
    match plan_bezier_conversion_with_poll(curve, &mut never_cancel)? {
        CurveWorkRun::Complete(plan) => Ok(plan),
        CurveWorkRun::Cancelled => Err(NurbsError::Domain {
            what: "non-cancelling Bezier planning observed cancellation".to_string(),
        }),
    }
}

fn plan_bezier_conversion_with_poll<S: Scalar, const DIM: usize>(
    curve: AdmittedNurbsCurve<'_, S, DIM>,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<CurveWorkRun<BezierConversionPlan>, NurbsError> {
    let knots = curve.knots();
    let degree = knots.degree();
    let (lo, hi) = knots.domain();
    let entries = knots.knots();
    let pre_scan_work = bezier_pre_scan_work(entries.len())?;
    if pre_scan_work > BASIS_MAX_WORK_UNITS {
        return Err(NurbsError::Domain {
            what: format!(
                "Bezier pre-scan requests {pre_scan_work} work units above defensive ceiling {BASIS_MAX_WORK_UNITS}"
            ),
        });
    }
    if should_cancel() {
        return Ok(CurveWorkRun::Cancelled);
    }

    let mut operations_since_poll = 0usize;
    let mut insertions = 0usize;
    let mut distinct_knot_count = 0usize;
    let mut run_start = 0usize;
    while run_start < entries.len() {
        distinct_knot_count =
            distinct_knot_count
                .checked_add(1)
                .ok_or_else(|| NurbsError::Domain {
                    what: "Bezier distinct-knot count overflows usize".to_string(),
                })?;
        let knot = entries[run_start];
        let mut run_end = run_start + 1;
        while run_end < entries.len() && entries[run_end] == knot {
            run_end += 1;
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveWorkRun::Cancelled);
            }
        }
        let multiplicity = run_end - run_start;
        if knot > lo && knot < hi && multiplicity < degree {
            insertions = insertions
                .checked_add(degree - multiplicity)
                .ok_or_else(|| NurbsError::Domain {
                    what: "Bezier insertion-count accounting overflows usize".to_string(),
                })?;
        }
        if curve_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(CurveWorkRun::Cancelled);
        }
        run_start = run_end;
    }

    let final_knot_count =
        entries
            .len()
            .checked_add(insertions)
            .ok_or_else(|| NurbsError::Domain {
                what: "Bezier final knot count overflows usize".to_string(),
            })?;
    let final_control_count = curve
        .homogeneous_control_points()
        .len()
        .checked_add(insertions)
        .ok_or_else(|| NurbsError::Domain {
            what: "Bezier final control count overflows usize".to_string(),
        })?;
    let work_units = refinement_work_upper_bound(
        degree,
        entries.len(),
        curve.homogeneous_control_points().len(),
        0,
        insertions,
    )?;
    if work_units > BASIS_MAX_WORK_UNITS {
        return Err(NurbsError::Domain {
            what: format!(
                "Bezier conversion requests {work_units} work units above defensive ceiling {BASIS_MAX_WORK_UNITS}"
            ),
        });
    }

    let converted_bytes = curve_storage_bytes::<S>(final_knot_count, final_control_count)?;
    let peak_allocated_bytes = if insertions == 0 {
        converted_bytes
    } else {
        converted_bytes
            .checked_mul(2)
            .ok_or_else(|| NurbsError::Domain {
                what: "Bezier conversion peak retained-byte accounting overflows u128".to_string(),
            })?
    };
    if peak_allocated_bytes > CURVE_BEZIER_MAX_RETAINED_BYTES {
        return Err(NurbsError::Domain {
            what: format!(
                "Bezier conversion can retain {peak_allocated_bytes} bytes above defensive ceiling {CURVE_BEZIER_MAX_RETAINED_BYTES}"
            ),
        });
    }

    if should_cancel() {
        return Ok(CurveWorkRun::Cancelled);
    }
    Ok(CurveWorkRun::Complete(BezierConversionPlan {
        insertions,
        distinct_knot_count,
        final_knot_count,
        final_control_count,
        work_units,
        peak_allocated_bytes,
        converted_bytes,
    }))
}

/// One span's Cartesian control box: (min, max, t0, t1).
pub type SpanBox<S, const DIM: usize> = ([S; DIM], [S; DIM], S, S);

fn enforce_span_box_envelope(work: u128, retained_bytes: u128) -> Result<(), NurbsError> {
    if work > BASIS_MAX_WORK_UNITS {
        return Err(NurbsError::Domain {
            what: format!(
                "curve span-box traversal requests {work} work units above defensive ceiling {BASIS_MAX_WORK_UNITS}"
            ),
        });
    }
    if retained_bytes > CURVE_SPAN_BOX_MAX_RETAINED_BYTES {
        return Err(NurbsError::Domain {
            what: format!(
                "curve span boxes retain {retained_bytes} bytes above defensive ceiling {CURVE_SPAN_BOX_MAX_RETAINED_BYTES}"
            ),
        });
    }
    Ok(())
}

fn preflight_span_boxes(
    control_count: usize,
    degree: usize,
    retained_bytes_per_box: usize,
) -> Result<usize, NurbsError> {
    let span_capacity = control_count
        .checked_sub(degree)
        .ok_or_else(|| NurbsError::Structure {
            what: "curve degree exceeds its admitted control count".to_string(),
        })?;
    let order = degree.checked_add(1).ok_or_else(|| NurbsError::Domain {
        what: "curve order overflows usize during span-box admission".to_string(),
    })?;
    let control_visits = (span_capacity as u128)
        .checked_mul(order as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "curve span-box control-scan work overflows u128".to_string(),
        })?;
    let traversal_work = (span_capacity as u128)
        .checked_mul(2)
        .and_then(|work| {
            control_visits
                .checked_mul(CURVE_SPAN_BOX_WORK_PER_CONTROL)
                .and_then(|control_work| work.checked_add(control_work))
        })
        .ok_or_else(|| NurbsError::Domain {
            what: "curve span-box traversal work overflows u128".to_string(),
        })?;
    let retained_bytes = (span_capacity as u128)
        .checked_mul(retained_bytes_per_box as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "curve span-box retained-byte accounting overflows u128".to_string(),
        })?;
    enforce_span_box_envelope(traversal_work, retained_bytes)?;
    Ok(span_capacity)
}

/// A rational curve in `DIM` dimensions: homogeneous control points
/// `(w·x…, w)` over a clamped knot vector.
///
/// ```compile_fail
/// use fs_rep_nurbs::{KnotVector, NurbsCurve};
/// let knots = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).unwrap();
/// let mut curve = NurbsCurve::new(knots, &[[0.0], [1.0]], &[1.0, 1.0]).unwrap();
/// curve.cpw.clear();
/// ```
#[derive(Debug, PartialEq)]
pub struct NurbsCurve<S: Scalar, const DIM: usize> {
    /// The knot vector.
    pub(crate) knots: KnotVector<S>,
    /// Homogeneous control points: `DIM` weighted coordinates + weight.
    pub(crate) cpw: Vec<[S; 4]>,
}

/// A validate-once borrow of one exact immutable curve snapshot.
#[derive(Debug, Clone, Copy)]
pub struct AdmittedNurbsCurve<'a, S: Scalar, const DIM: usize> {
    inner: &'a NurbsCurve<S, DIM>,
}

/// Transactional terminal state of cancellation-aware curve construction.
#[must_use]
#[derive(Debug, PartialEq)]
pub enum CurveConstructionRun<S: Scalar, const DIM: usize> {
    /// Validation completed and the sealed curve is safe to publish.
    Complete {
        /// Newly validated exact curve generation.
        curve: NurbsCurve<S, DIM>,
    },
    /// Cancellation was observed; the unpublished owned candidate was dropped.
    Cancelled,
}

/// Transactional terminal state of a cancellation-aware fallible curve copy.
#[must_use]
#[derive(Debug, PartialEq)]
pub enum CurveCloneRun<S: Scalar, const DIM: usize> {
    /// The complete sealed copy of the exact source representation.
    Complete {
        /// Copied curve generation.
        curve: NurbsCurve<S, DIM>,
    },
    /// Cancellation was observed; all partial copy storage was dropped.
    Cancelled,
}

/// Transactional terminal state of cancellation-aware curve admission.
#[must_use]
#[derive(Debug, Clone, Copy)]
pub enum CurveAdmissionRun<'a, S: Scalar, const DIM: usize> {
    /// The exact immutable source snapshot was fully validated.
    Complete {
        /// Lifetime-bound authority for the validated curve generation.
        admitted: AdmittedNurbsCurve<'a, S, DIM>,
    },
    /// Cancellation was observed; no admitted authority was published.
    Cancelled,
}

/// Transactional terminal state of cancellation-aware homogeneous evaluation.
#[must_use]
#[derive(Debug, Clone, PartialEq)]
pub enum CurveHomogeneousEvaluationRun<S: Scalar> {
    /// The complete finite homogeneous point `(w*x..., w)`.
    Complete {
        /// Evaluated four-lane homogeneous storage.
        homogeneous: [S; 4],
    },
    /// Cancellation was observed; no partial homogeneous point was published.
    Cancelled,
}

/// Transactional terminal state of cancellation-aware Cartesian evaluation.
#[must_use]
#[derive(Debug, Clone, PartialEq)]
pub enum CurveEvaluationRun<S: Scalar, const DIM: usize> {
    /// The complete finite Cartesian point.
    Complete {
        /// Evaluated point in the curve's declared dimension.
        point: [S; DIM],
    },
    /// Cancellation was observed; no partial point was published.
    Cancelled,
}

/// Transactional terminal state of cancellation-aware f64 derivatives.
#[must_use]
#[derive(Debug, Clone, PartialEq)]
pub enum CurveDerivativesRun<const DIM: usize> {
    /// The complete Cartesian jet from order zero through the requested order.
    Complete {
        /// Cartesian derivatives in ascending order.
        derivatives: Vec<[f64; DIM]>,
    },
    /// Cancellation was observed; no partial derivative jet was published.
    Cancelled,
}

/// Transactional terminal state of cancellation-aware exact knot insertion.
#[must_use]
#[derive(Debug, PartialEq)]
pub enum CurveInsertionRun<S: Scalar, const DIM: usize> {
    /// The complete sealed and validated derived generation.
    Complete {
        /// Exact refinement of the admitted source curve.
        curve: NurbsCurve<S, DIM>,
    },
    /// Cancellation was observed; all partial derived storage was dropped.
    Cancelled,
}

/// Transactional terminal state of cancellation-aware exact knot removal.
#[must_use]
#[derive(Debug, PartialEq)]
pub enum CurveRemovalRun<S: Scalar, const DIM: usize> {
    /// The complete sealed, validated, and reinsertion-verified generation.
    Complete {
        /// Exact coarsening of the admitted source curve.
        curve: NurbsCurve<S, DIM>,
    },
    /// Cancellation was observed; all partial derived storage was dropped.
    Cancelled,
}

/// Transactional terminal state of cancellation-aware curve span boxes.
#[must_use]
#[derive(Debug, Clone, PartialEq)]
pub enum CurveSpanBoxesRun<S: Scalar, const DIM: usize> {
    /// Every admitted nonempty span has a complete Cartesian control box.
    Complete {
        /// Boxes in deterministic source-span order.
        boxes: Vec<SpanBox<S, DIM>>,
    },
    /// Cancellation was observed; no partial box table was published.
    Cancelled,
}

/// Transactional terminal state of cancellation-aware exact Bezier conversion.
#[must_use]
#[derive(Debug, PartialEq)]
pub enum CurveBezierRun<S: Scalar, const DIM: usize> {
    /// The complete sealed and validated Bezier-form generation.
    Complete {
        /// Exact refinement of the admitted source curve.
        curve: NurbsCurve<S, DIM>,
    },
    /// Cancellation was observed; all partial derived generations were dropped.
    Cancelled,
}

/// Transactional terminal state of cancellation-aware exact degree elevation.
#[must_use]
#[derive(Debug, PartialEq)]
pub enum CurveElevationRun<S: Scalar, const DIM: usize> {
    /// The complete sealed and validated elevated generation.
    Complete {
        /// Exact degree elevation of the admitted source curve.
        curve: NurbsCurve<S, DIM>,
    },
    /// Cancellation was observed; all partial derived storage was dropped.
    Cancelled,
}

impl<S: Scalar, const DIM: usize> NurbsCurve<S, DIM> {
    fn validation_work_for(
        knots: &KnotVector<S>,
        control_count: usize,
    ) -> Result<u128, NurbsError> {
        knots
            .validation_work()?
            .checked_add(
                (control_count as u128)
                    .checked_mul(CURVE_VALIDATION_WORK_PER_CONTROL)
                    .ok_or_else(|| NurbsError::Domain {
                        what: "curve control-validation work overflows u128".to_string(),
                    })?,
            )
            .ok_or_else(|| NurbsError::Domain {
                what: "curve structure-validation work overflows u128".to_string(),
            })
    }

    fn enforce_validation_work(work: u128) -> Result<(), NurbsError> {
        if work > BASIS_MAX_WORK_UNITS {
            return Err(NurbsError::Domain {
                what: format!(
                    "curve structure validation requests {work} work units above defensive ceiling {BASIS_MAX_WORK_UNITS}"
                ),
            });
        }
        Ok(())
    }

    pub(crate) fn validate_live_structure(&self) -> Result<(), NurbsError> {
        let mut never_cancel = || false;
        match self.validate_live_structure_with_poll(&mut never_cancel)? {
            CurveWorkRun::Complete(()) => Ok(()),
            CurveWorkRun::Cancelled => Err(NurbsError::Domain {
                what: "non-cancelling curve validation observed cancellation".to_string(),
            }),
        }
    }

    fn validate_live_structure_with_poll(
        &self,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveWorkRun<()>, NurbsError> {
        if DIM > 3 {
            return Err(NurbsError::Structure {
                what: format!("curve dimension {DIM} exceeds the homogeneous storage limit 3"),
            });
        }
        Self::enforce_validation_work(Self::validation_work_for(&self.knots, self.cpw.len())?)?;
        match self.knots.validate_live_with_poll(|| should_cancel())? {
            KnotValidationOutcome::Complete => {}
            KnotValidationOutcome::Cancelled => return Ok(CurveWorkRun::Cancelled),
        }

        let invalid_control_net = || {
            NurbsError::Structure {
            what: "live curve control net must match its knots, retain finite homogeneous coordinates with admissible weights, and zero inactive coordinate lanes"
                .to_string(),
        }
        };
        if self.cpw.len() != self.knots.control_count() {
            return Err(invalid_control_net());
        }
        if should_cancel() {
            return Ok(CurveWorkRun::Cancelled);
        }

        let mut operations_since_poll = 0usize;
        for control in &self.cpw {
            if !control[3].is_admissible_weight() {
                return Err(invalid_control_net());
            }
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveWorkRun::Cancelled);
            }
            for &component in control {
                if !component.is_finite() {
                    return Err(invalid_control_net());
                }
                if curve_poll_due(&mut operations_since_poll, should_cancel) {
                    return Ok(CurveWorkRun::Cancelled);
                }
            }
            for &component in &control[..DIM] {
                if !component.quotient_is_finite(control[3]) {
                    return Err(invalid_control_net());
                }
                if curve_poll_due(&mut operations_since_poll, should_cancel) {
                    return Ok(CurveWorkRun::Cancelled);
                }
            }
            for &component in &control[DIM..3] {
                if component != S::zero() {
                    return Err(invalid_control_net());
                }
                if curve_poll_due(&mut operations_since_poll, should_cancel) {
                    return Ok(CurveWorkRun::Cancelled);
                }
            }
        }
        if should_cancel() {
            return Ok(CurveWorkRun::Cancelled);
        }
        Ok(CurveWorkRun::Complete(()))
    }

    /// Build from Cartesian control points + weights.
    ///
    /// # Errors
    /// [`NurbsError::Structure`] on count mismatch, non-finite coordinates, or
    /// non-positive/non-finite weights; [`NurbsError::Domain`] when validation
    /// work, the 64 MiB derived-control envelope, or homogeneous-control
    /// allocation is refused.
    pub fn new(
        knots: KnotVector<S>,
        points: &[[S; DIM]],
        weights: &[S],
    ) -> Result<Self, NurbsError> {
        let mut never_cancel = || false;
        match Self::new_with_poll(knots, points, weights, &mut never_cancel)? {
            CurveWorkRun::Complete(curve) => Ok(curve),
            CurveWorkRun::Cancelled => Err(NurbsError::Domain {
                what: "non-cancelling Cartesian curve construction observed cancellation"
                    .to_string(),
            }),
        }
    }

    /// Build from borrowed Cartesian controls and weights with bounded
    /// cancellation polling.
    ///
    /// Dimension, count, aggregate validation-work, and 64 MiB derived-output
    /// refusals precede cancellation. One `Cx` then spans knot validation,
    /// ordered weight and coordinate validation, fallible homogeneous-output
    /// allocation, Cartesian-to-homogeneous multiplication, underflow and
    /// overflow checks, and final owned publication. Cancellation drops the
    /// transferred knot vector and any partial derived output; the borrowed
    /// point and weight slices remain caller-owned. Individual allocator,
    /// scalar, and destructor operations are not preemptible. This primitive
    /// does not consume the `Cx` budget or own request -> drain -> finalize
    /// semantics.
    ///
    /// # Errors
    /// Returns the synchronous constructor's dimension, count, work,
    /// retained-memory, knot, weight, coordinate, arithmetic, or allocation
    /// refusal when it wins before an observed cancellation.
    pub fn new_with_cx(
        knots: KnotVector<S>,
        points: &[[S; DIM]],
        weights: &[S],
        cx: &Cx<'_>,
    ) -> Result<CurveConstructionRun<S, DIM>, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        match Self::new_with_poll(knots, points, weights, &mut should_cancel)? {
            CurveWorkRun::Complete(curve) => Ok(CurveConstructionRun::Complete { curve }),
            CurveWorkRun::Cancelled => Ok(CurveConstructionRun::Cancelled),
        }
    }

    fn new_with_poll(
        knots: KnotVector<S>,
        points: &[[S; DIM]],
        weights: &[S],
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveWorkRun<Self>, NurbsError> {
        if DIM > 3 {
            return Err(NurbsError::Structure {
                what: format!("curve dimension {DIM} exceeds the homogeneous storage limit 3"),
            });
        }
        if points.len() != knots.control_count() || weights.len() != points.len() {
            return Err(NurbsError::Structure {
                what: format!(
                    "knot vector wants {} control points, got {} points / {} weights",
                    knots.control_count(),
                    points.len(),
                    weights.len()
                ),
            });
        }
        Self::enforce_validation_work(Self::validation_work_for(&knots, points.len())?)?;
        preflight_cartesian_curve_construction::<S>(points.len())?;
        match knots.validate_live_with_poll(|| should_cancel())? {
            KnotValidationOutcome::Complete => {}
            KnotValidationOutcome::Cancelled => return Ok(CurveWorkRun::Cancelled),
        }
        if matches!(
            validate_cartesian_curve_inputs_with_poll(points, weights, should_cancel)?,
            CurveWorkRun::Cancelled
        ) {
            return Ok(CurveWorkRun::Cancelled);
        }
        let cpw = match build_cartesian_curve_controls_with_poll(points, weights, should_cancel)? {
            CurveWorkRun::Complete(cpw) => cpw,
            CurveWorkRun::Cancelled => return Ok(CurveWorkRun::Cancelled),
        };
        if matches!(
            validate_cartesian_curve_products_with_poll(points, &cpw, should_cancel)?,
            CurveWorkRun::Cancelled
        ) {
            return Ok(CurveWorkRun::Cancelled);
        }
        Ok(CurveWorkRun::Complete(NurbsCurve { knots, cpw }))
    }

    /// Build from a homogeneous control net, validating every coordinate and
    /// weight before the sealed representation is exposed.
    ///
    /// # Errors
    /// [`NurbsError::Structure`] on a knot/control count mismatch, an invalid
    /// knot vector, a non-finite coordinate, a noncanonical inactive lane, or
    /// an inadmissible weight; [`NurbsError::Domain`] when validation work is
    /// refused.
    pub fn from_homogeneous(knots: KnotVector<S>, cpw: Vec<[S; 4]>) -> Result<Self, NurbsError> {
        let mut never_cancel = || false;
        match Self::from_homogeneous_with_poll(knots, cpw, &mut never_cancel)? {
            CurveWorkRun::Complete(curve) => Ok(curve),
            CurveWorkRun::Cancelled => Err(NurbsError::Domain {
                what: "non-cancelling homogeneous curve construction observed cancellation"
                    .to_string(),
            }),
        }
    }

    /// Build from an owned homogeneous control net with bounded cancellation
    /// polling.
    ///
    /// Dimension and aggregate validation-work refusals precede cancellation.
    /// One `Cx` then spans knot validation, control-count, weight, finite,
    /// Cartesian-projection, and inactive-lane checks plus final owned
    /// publication. Cancellation drops both caller-transferred inputs without
    /// exposing a partially validated curve. Individual scalar operations and
    /// destruction are not preemptible. This primitive does not consume the
    /// `Cx` budget or own request -> drain -> finalize semantics.
    ///
    /// # Errors
    /// Returns the synchronous constructor's dimension, work, knot,
    /// control-count, weight, finite-arithmetic, or canonical-lane refusal when
    /// it wins before an observed cancellation.
    pub fn from_homogeneous_with_cx(
        knots: KnotVector<S>,
        cpw: Vec<[S; 4]>,
        cx: &Cx<'_>,
    ) -> Result<CurveConstructionRun<S, DIM>, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        match Self::from_homogeneous_with_poll(knots, cpw, &mut should_cancel)? {
            CurveWorkRun::Complete(curve) => Ok(CurveConstructionRun::Complete { curve }),
            CurveWorkRun::Cancelled => Ok(CurveConstructionRun::Cancelled),
        }
    }

    fn from_homogeneous_with_poll(
        knots: KnotVector<S>,
        cpw: Vec<[S; 4]>,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveWorkRun<Self>, NurbsError> {
        let candidate = NurbsCurve { knots, cpw };
        match candidate.validate_live_structure_with_poll(should_cancel)? {
            CurveWorkRun::Complete(()) => Ok(CurveWorkRun::Complete(candidate)),
            CurveWorkRun::Cancelled => Ok(CurveWorkRun::Cancelled),
        }
    }

    /// Borrow the curve's sealed knot vector.
    #[must_use]
    pub const fn knots(&self) -> &KnotVector<S> {
        &self.knots
    }

    /// Borrow the sealed homogeneous control points.
    #[must_use]
    pub fn homogeneous_control_points(&self) -> &[[S; 4]] {
        &self.cpw
    }

    /// Fallibly copy this sealed curve without revalidating unchanged data.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] when checked copy work/retained bytes or a
    /// destination allocation is refused.
    pub fn try_clone(&self) -> Result<Self, NurbsError> {
        let mut never_cancel = || false;
        match self.try_clone_with_poll(&mut never_cancel)? {
            CurveWorkRun::Complete(curve) => Ok(curve),
            CurveWorkRun::Cancelled => Err(NurbsError::Domain {
                what: "non-cancelling curve copy observed cancellation".to_string(),
            }),
        }
    }

    /// Fallibly copy this sealed curve with bounded cancellation polling.
    ///
    /// Count-derived work and a 64 MiB retained-output envelope precede
    /// cancellation. One gate then covers both fallible allocations, ordered
    /// knot and control copies at fixed logical-work strides, and final
    /// publication. The borrowed source is excluded from the output envelope.
    /// Individual allocator calls, scalar/array copies, and destructors are
    /// not preemptible. This primitive does not consume the `Cx` budget or own
    /// request -> drain -> finalize semantics.
    ///
    /// # Errors
    /// Returns the synchronous copy's work, retained-memory, or allocation
    /// refusal when it wins before an observed cancellation.
    pub fn try_clone_with_cx(&self, cx: &Cx<'_>) -> Result<CurveCloneRun<S, DIM>, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        match self.try_clone_with_poll(&mut should_cancel)? {
            CurveWorkRun::Complete(curve) => Ok(CurveCloneRun::Complete { curve }),
            CurveWorkRun::Cancelled => Ok(CurveCloneRun::Cancelled),
        }
    }

    fn try_clone_with_poll(
        &self,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveWorkRun<Self>, NurbsError> {
        preflight_curve_copy::<S>(self.knots.knots.len(), self.cpw.len())?;
        if should_cancel() {
            return Ok(CurveWorkRun::Cancelled);
        }
        let mut knot_entries = Vec::new();
        knot_entries
            .try_reserve_exact(self.knots.knots.len())
            .map_err(|_| NurbsError::Domain {
                what: "knot-vector copy allocation was refused".to_string(),
            })?;
        let mut operations_since_poll = 0usize;
        for &knot in &self.knots.knots {
            knot_entries.push(knot);
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveWorkRun::Cancelled);
            }
        }

        if should_cancel() {
            return Ok(CurveWorkRun::Cancelled);
        }
        let mut cpw = Vec::new();
        cpw.try_reserve_exact(self.cpw.len())
            .map_err(|_| NurbsError::Domain {
                what: "curve copy control-net allocation was refused".to_string(),
            })?;
        operations_since_poll = 0;
        for &control in &self.cpw {
            cpw.push(control);
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveWorkRun::Cancelled);
            }
        }
        if should_cancel() {
            return Ok(CurveWorkRun::Cancelled);
        }
        Ok(CurveWorkRun::Complete(NurbsCurve {
            knots: KnotVector {
                knots: knot_entries,
                degree: self.knots.degree,
            },
            cpw,
        }))
    }

    /// Validate this exact immutable curve snapshot once.
    ///
    /// # Errors
    /// [`NurbsError::Structure`] when the sealed source is internally invalid;
    /// [`NurbsError::Domain`] when validation work exceeds the defensive cap.
    pub fn admit(&self) -> Result<AdmittedNurbsCurve<'_, S, DIM>, NurbsError> {
        self.validate_live_structure()?;
        Ok(self.admitted_after_validation())
    }

    /// Validate this immutable curve with bounded cancellation polling.
    ///
    /// Dimension and checked validation-work refusal retain their synchronous
    /// precedence. The knot and homogeneous-control scans share one gate, and
    /// a final checkpoint gates publication of the lifetime-bound admitted
    /// view. Cancellation-aware homogeneous ownership transfer is provided
    /// separately by [`Self::from_homogeneous_with_cx`]. This method does not
    /// consume the `Cx` budget or finalize its executor scope.
    ///
    /// # Errors
    /// Returns the synchronous admission's dimension, work, knot, control-count,
    /// weight, finite-arithmetic, and canonical-lane refusals when they win
    /// before an observed cancellation.
    pub fn admit_with_cx<'a>(
        &'a self,
        cx: &Cx<'_>,
    ) -> Result<CurveAdmissionRun<'a, S, DIM>, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        match self.validate_live_structure_with_poll(&mut should_cancel)? {
            CurveWorkRun::Complete(()) => Ok(CurveAdmissionRun::Complete {
                admitted: self.admitted_after_validation(),
            }),
            CurveWorkRun::Cancelled => Ok(CurveAdmissionRun::Cancelled),
        }
    }

    pub(crate) const fn admitted_after_validation(&self) -> AdmittedNurbsCurve<'_, S, DIM> {
        AdmittedNurbsCurve { inner: self }
    }

    /// Homogeneous evaluation (the shared exact/fast core).
    ///
    /// # Errors
    /// [`NurbsError::Domain`] outside the domain.
    pub fn eval_homogeneous(&self, t: S) -> Result<[S; 4], NurbsError> {
        self.admit()?.eval_homogeneous(t)
    }

    /// Validate this owning curve and evaluate its homogeneous point with one
    /// cancellation gate.
    ///
    /// Dimension and checked structural-validation work refusals precede the
    /// first checkpoint. Cancellation then spans structural admission and the
    /// same admitted evaluation pipeline as
    /// [`AdmittedNurbsCurve::eval_homogeneous_with_cx`]. No partial admitted
    /// authority or homogeneous representation is published. This primitive
    /// does not consume the `Cx` budget or finalize its executor scope.
    ///
    /// # Errors
    /// Returns the synchronous owning evaluator's structure, parameter, work,
    /// allocation, and finite-arithmetic refusals when they win before an
    /// observed cancellation.
    pub fn eval_homogeneous_with_cx(
        &self,
        t: S,
        cx: &Cx<'_>,
    ) -> Result<CurveHomogeneousEvaluationRun<S>, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        self.eval_homogeneous_with_poll(t, &mut should_cancel)
    }

    fn eval_homogeneous_with_poll(
        &self,
        t: S,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveHomogeneousEvaluationRun<S>, NurbsError> {
        match self.validate_live_structure_with_poll(should_cancel)? {
            CurveWorkRun::Complete(()) => self
                .admitted_after_validation()
                .eval_homogeneous_with_poll(t, should_cancel),
            CurveWorkRun::Cancelled => Ok(CurveHomogeneousEvaluationRun::Cancelled),
        }
    }

    /// Cartesian evaluation.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] outside the domain.
    pub fn eval(&self, t: S) -> Result<[S; DIM], NurbsError> {
        self.admit()?.eval(t)
    }

    /// Validate this owning curve and evaluate its Cartesian point with one
    /// cancellation gate.
    ///
    /// Dimension and checked structural-validation work refusals precede the
    /// first checkpoint. Cancellation then spans structural admission and the
    /// same admitted evaluation pipeline as
    /// [`AdmittedNurbsCurve::eval_with_cx`]. No partial admitted authority or
    /// Cartesian point is published. This primitive does not consume the `Cx`
    /// budget or finalize its executor scope.
    ///
    /// # Errors
    /// Returns the synchronous owning evaluator's structure, parameter, work,
    /// allocation, weight, and finite-arithmetic refusals when they win before
    /// an observed cancellation.
    pub fn eval_with_cx(
        &self,
        t: S,
        cx: &Cx<'_>,
    ) -> Result<CurveEvaluationRun<S, DIM>, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        self.eval_with_poll(t, &mut should_cancel)
    }

    fn eval_with_poll(
        &self,
        t: S,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveEvaluationRun<S, DIM>, NurbsError> {
        match self.validate_live_structure_with_poll(should_cancel)? {
            CurveWorkRun::Complete(()) => self
                .admitted_after_validation()
                .eval_with_poll(t, should_cancel),
            CurveWorkRun::Cancelled => Ok(CurveEvaluationRun::Cancelled),
        }
    }
}

impl<'a, S: Scalar, const DIM: usize> AdmittedNurbsCurve<'a, S, DIM> {
    /// The exact immutable source bound to this view.
    #[must_use]
    pub const fn source(&self) -> &'a NurbsCurve<S, DIM> {
        self.inner
    }

    /// Borrow the already-validated knot-vector view.
    #[must_use]
    pub fn knots(&self) -> crate::basis::AdmittedKnotVector<'a, S> {
        self.inner.knots.admitted_after_validation()
    }

    /// Borrow the sealed homogeneous control points.
    #[must_use]
    pub fn homogeneous_control_points(&self) -> &'a [[S; 4]] {
        &self.inner.cpw
    }

    /// Homogeneous evaluation without rescanning curve or knot structure.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] outside the domain or when basis work/allocation
    /// is refused.
    pub fn eval_homogeneous(&self, t: S) -> Result<[S; 4], NurbsError> {
        let mut never_cancel = || false;
        match self.eval_homogeneous_with_poll(t, &mut never_cancel)? {
            CurveHomogeneousEvaluationRun::Complete { homogeneous } => Ok(homogeneous),
            CurveHomogeneousEvaluationRun::Cancelled => Err(NurbsError::Domain {
                what: "non-cancelling homogeneous curve evaluation observed cancellation"
                    .to_string(),
            }),
        }
    }

    /// Evaluate a homogeneous point with bounded cancellation polling.
    ///
    /// This admitted-only path reuses the exact immutable curve generation,
    /// delegates cancellable basis construction to its admitted knot view,
    /// polls the four-lane accumulation at fixed logical-work strides, and
    /// gates final homogeneous publication. It does not divide by the weight
    /// or make a Cartesian-finiteness claim. The caller remains responsible
    /// for owning admission, `Cx` budget consumption, and request -> drain ->
    /// finalize semantics around this primitive.
    ///
    /// # Errors
    /// Returns the same parameter, work, allocation, and homogeneous
    /// finite-arithmetic refusals as [`Self::eval_homogeneous`] when they win
    /// before an observed cancellation.
    pub fn eval_homogeneous_with_cx(
        &self,
        t: S,
        cx: &Cx<'_>,
    ) -> Result<CurveHomogeneousEvaluationRun<S>, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        self.eval_homogeneous_with_poll(t, &mut should_cancel)
    }

    /// Evaluate an admitted homogeneous point while sharing a compound
    /// caller's cancellation callback.
    pub(crate) fn eval_homogeneous_with_poll(
        &self,
        t: S,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveHomogeneousEvaluationRun<S>, NurbsError> {
        let (span, basis) = match self.knots().basis_with_poll(t, should_cancel)? {
            BasisRun::Complete { span, values } => (span, values),
            BasisRun::Cancelled => return Ok(CurveHomogeneousEvaluationRun::Cancelled),
        };
        self.eval_homogeneous_from_basis_with_poll(span, &basis, should_cancel)
    }

    fn eval_homogeneous_from_basis_with_poll(
        &self,
        span: usize,
        basis: &[S],
        mut should_cancel: impl FnMut() -> bool,
    ) -> Result<CurveHomogeneousEvaluationRun<S>, NurbsError> {
        let (homogeneous, _) = match self.accumulate_homogeneous_from_basis_with_poll(
            span,
            basis,
            &mut should_cancel,
        )? {
            CurveWorkRun::Complete(accumulation) => accumulation,
            CurveWorkRun::Cancelled => return Ok(CurveHomogeneousEvaluationRun::Cancelled),
        };
        if should_cancel() {
            return Ok(CurveHomogeneousEvaluationRun::Cancelled);
        }
        Ok(CurveHomogeneousEvaluationRun::Complete { homogeneous })
    }

    /// Cartesian evaluation.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] outside the domain.
    pub fn eval(&self, t: S) -> Result<[S; DIM], NurbsError> {
        if DIM > 3 {
            return Err(NurbsError::Structure {
                what: format!("curve dimension {DIM} exceeds the homogeneous storage limit 3"),
            });
        }
        let h = self.eval_homogeneous(t)?;
        if !h[3].is_admissible_weight() {
            return Err(NurbsError::Domain {
                what: "curve evaluation produced an inadmissible rational denominator".to_string(),
            });
        }
        let mut out = [S::zero(); DIM];
        for (o, &c) in out.iter_mut().zip(h.iter()) {
            *o = c / h[3];
        }
        if out.iter().copied().any(|component| !component.is_finite()) {
            return Err(NurbsError::Domain {
                what: "Cartesian curve evaluation left the finite numeric domain".to_string(),
            });
        }
        Ok(out)
    }

    /// Evaluate a Cartesian point with bounded cancellation polling.
    ///
    /// This entry point is deliberately admitted-only. It reuses the exact
    /// immutable curve generation, delegates cancellable basis construction
    /// to its admitted knot view, polls homogeneous accumulation at fixed
    /// logical-work strides, and gates final point publication. The caller
    /// remains responsible for owning admission, `Cx` budget consumption, and
    /// request -> drain -> finalize semantics around this primitive.
    ///
    /// # Errors
    /// Returns the same dimension, parameter, work, allocation, weight, and
    /// finite-arithmetic refusals as [`Self::eval`] when they win before an
    /// observed cancellation.
    pub fn eval_with_cx(
        &self,
        t: S,
        cx: &Cx<'_>,
    ) -> Result<CurveEvaluationRun<S, DIM>, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        self.eval_with_poll(t, &mut should_cancel)
    }

    /// Evaluate an admitted point while sharing a compound caller's
    /// cancellation callback.
    pub(crate) fn eval_with_poll(
        &self,
        t: S,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveEvaluationRun<S, DIM>, NurbsError> {
        if DIM > 3 {
            return Err(NurbsError::Structure {
                what: format!("curve dimension {DIM} exceeds the homogeneous storage limit 3"),
            });
        }
        let (span, basis) = match self.knots().basis_with_poll(t, should_cancel)? {
            BasisRun::Complete { span, values } => (span, values),
            BasisRun::Cancelled => return Ok(CurveEvaluationRun::Cancelled),
        };
        self.eval_from_basis_with_poll(span, &basis, should_cancel)
    }

    fn eval_from_basis_with_poll(
        &self,
        span: usize,
        basis: &[S],
        mut should_cancel: impl FnMut() -> bool,
    ) -> Result<CurveEvaluationRun<S, DIM>, NurbsError> {
        let (homogeneous, mut operations_since_poll) = match self
            .accumulate_homogeneous_from_basis_with_poll(span, basis, &mut should_cancel)?
        {
            CurveWorkRun::Complete(accumulation) => accumulation,
            CurveWorkRun::Cancelled => return Ok(CurveEvaluationRun::Cancelled),
        };
        if !homogeneous[3].is_admissible_weight() {
            return Err(NurbsError::Domain {
                what: "curve evaluation produced an inadmissible rational denominator".to_string(),
            });
        }

        let mut point = [S::zero(); DIM];
        for (slot, &component) in point.iter_mut().zip(&homogeneous) {
            *slot = component / homogeneous[3];
            if !slot.is_finite() {
                return Err(NurbsError::Domain {
                    what: "Cartesian curve evaluation left the finite numeric domain".to_string(),
                });
            }
            if curve_poll_due(&mut operations_since_poll, &mut should_cancel) {
                return Ok(CurveEvaluationRun::Cancelled);
            }
        }
        if should_cancel() {
            return Ok(CurveEvaluationRun::Cancelled);
        }
        Ok(CurveEvaluationRun::Complete { point })
    }

    fn accumulate_homogeneous_from_basis_with_poll(
        &self,
        span: usize,
        basis: &[S],
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveWorkRun<([S; 4], usize)>, NurbsError> {
        if should_cancel() {
            return Ok(CurveWorkRun::Cancelled);
        }

        let p = self.knots().degree();
        let mut operations_since_poll = 0usize;
        let mut homogeneous = [S::zero(); 4];
        for (r, &coefficient) in basis.iter().enumerate() {
            let control = &self.inner.cpw[span - p + r];
            for (accumulator, &component) in homogeneous.iter_mut().zip(control) {
                *accumulator = *accumulator + coefficient * component;
                if curve_poll_due(&mut operations_since_poll, should_cancel) {
                    return Ok(CurveWorkRun::Cancelled);
                }
            }
        }
        for &component in &homogeneous {
            if !component.is_finite() {
                return Err(NurbsError::Domain {
                    what: "homogeneous curve evaluation left the finite numeric domain".to_string(),
                });
            }
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveWorkRun::Cancelled);
            }
        }
        Ok(CurveWorkRun::Complete((homogeneous, operations_since_poll)))
    }

    /// Insert one knot while reusing this exact source admission.
    ///
    /// The returned curve is a new sealed generation and is validated before
    /// publication; only the unchanged source scan is elided.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] outside the open domain or when checked output
    /// sizing/allocation is refused; [`NurbsError::Structure`] if the derived
    /// generation is not a valid curve.
    pub fn insert_knot(&self, t: S) -> Result<NurbsCurve<S, DIM>, NurbsError> {
        self.inner.insert_knot_after_validation(t)
    }

    /// Insert one knot with bounded cancellation polling.
    ///
    /// Open-domain validation and the checked direct-insertion work/retained
    /// output envelope precede cancellation. One gate then covers span search,
    /// output allocation and copies, homogeneous blends, both derived
    /// structural validation passes, and final generation publication.
    /// Cancellation drops every partial derived allocation. The caller retains
    /// responsibility for owning admission, `Cx` budget consumption, and
    /// request -> drain -> finalize semantics.
    ///
    /// # Errors
    /// Returns the synchronous insertion's parameter, work, retained-memory,
    /// allocation, numeric-domain, or structural refusal when it wins before
    /// an observed cancellation.
    pub fn insert_knot_with_cx(
        &self,
        t: S,
        cx: &Cx<'_>,
    ) -> Result<CurveInsertionRun<S, DIM>, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        self.insert_knot_with_poll(t, &mut should_cancel)
    }

    /// Insert one knot while sharing a compound caller's cancellation
    /// callback across span lookup and derived-generation assembly.
    pub(crate) fn insert_knot_with_poll(
        &self,
        t: S,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveInsertionRun<S, DIM>, NurbsError> {
        let plan = self.inner.insertion_plan_after_parameter(t)?;
        let span = match self.knots().span_with_poll(t, should_cancel)? {
            KnotSpanRun::Complete { span } => span,
            KnotSpanRun::Cancelled => return Ok(CurveInsertionRun::Cancelled),
        };
        match self
            .inner
            .insert_knot_at_span_with_plan_and_poll(t, span, plan, should_cancel)?
        {
            CurveWorkRun::Complete(curve) => Ok(CurveInsertionRun::Complete { curve }),
            CurveWorkRun::Cancelled => Ok(CurveInsertionRun::Cancelled),
        }
    }

    /// Remove one exactly redundant knot while reusing this source admission.
    ///
    /// The returned curve is a new sealed generation. Exact reinsertion must
    /// reproduce every source knot and homogeneous control before publication.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] when `t` is outside the open domain, absent, or
    /// checked work/memory/allocation is refused; [`NurbsError::Structure`]
    /// when the knot is not exactly removable or a derived generation fails
    /// validation.
    pub fn remove_knot(&self, t: S) -> Result<NurbsCurve<S, DIM>, NurbsError> {
        let mut never_cancel = || false;
        match self.remove_knot_with_poll(t, &mut never_cancel)? {
            CurveRemovalRun::Complete { curve } => Ok(curve),
            CurveRemovalRun::Cancelled => Err(NurbsError::Domain {
                what: "non-cancelling knot removal observed cancellation".to_string(),
            }),
        }
    }

    /// Remove one exactly redundant knot with bounded cancellation polling.
    ///
    /// The open-domain check and count-derived aggregate work/64 MiB
    /// simultaneously-live derived-storage envelope precede cancellation.
    /// One gate then covers occurrence discovery, fallible knot/control
    /// allocation and copying, exact forward reconstruction, both candidate
    /// validation passes, exact reinsertion, representation comparison, and
    /// final publication. The restored verifier is dropped before the final
    /// checkpoint, but that individual destructor is not preemptible.
    /// Cancellation publishes no partial curve. The caller retains
    /// responsibility for owning admission, `Cx` budget consumption, and
    /// request -> drain -> finalize semantics.
    ///
    /// # Errors
    /// Returns the synchronous removal's parameter, work, memory, allocation,
    /// exact-reconstruction, or structural refusal when it wins before an
    /// observed cancellation.
    pub fn remove_knot_with_cx(
        &self,
        t: S,
        cx: &Cx<'_>,
    ) -> Result<CurveRemovalRun<S, DIM>, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        self.remove_knot_with_poll(t, &mut should_cancel)
    }

    /// Remove an admitted knot while sharing a compound caller's cancellation
    /// callback across reconstruction and exact reinsertion verification.
    pub(crate) fn remove_knot_with_poll(
        &self,
        t: S,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveRemovalRun<S, DIM>, NurbsError> {
        self.inner
            .remove_knot_after_validation_with_poll(t, should_cancel)
    }

    /// Convert this admitted generation to exact Bezier form without
    /// repeating structural validation of the unchanged source. Planning still
    /// scans knot runs once to determine the exact refinement envelope.
    ///
    /// Each newly derived generation is still validated before it becomes the
    /// input to the next insertion.
    ///
    /// # Errors
    /// Propagates checked sizing, allocation, numeric-domain, and structural
    /// refusals from exact knot insertion.
    pub fn to_bezier_form(&self) -> Result<NurbsCurve<S, DIM>, NurbsError> {
        self.inner.to_bezier_form_after_validation()
    }

    /// Convert this admitted generation to exact Bezier form with bounded
    /// cancellation polling.
    ///
    /// Planning, source copies, target scans, insertion copies/blends, every
    /// derived structural validation, and final publication observe the same
    /// gate. Cancellation drops all partial derived generations. The caller
    /// remains responsible for owning curve admission, `Cx` budget
    /// consumption, and request -> drain -> finalize semantics.
    ///
    /// # Errors
    /// Propagates the synchronous conversion's checked work, retained-memory,
    /// allocation, numeric-domain, and structural refusals when they win
    /// before an observed cancellation.
    pub fn to_bezier_form_with_cx(
        &self,
        cx: &Cx<'_>,
    ) -> Result<CurveBezierRun<S, DIM>, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        self.to_bezier_form_with_poll(&mut should_cancel)
    }

    /// Return the checked Bezier conversion envelope without allocating.
    #[cfg(test)]
    pub(crate) fn bezier_conversion_plan(&self) -> Result<BezierConversionPlan, NurbsError> {
        plan_bezier_conversion(*self)
    }

    /// Return the checked Bezier conversion envelope with bounded polling of
    /// its knot-run scan. `None` means cancellation won before publication.
    pub(crate) fn bezier_conversion_plan_with_poll(
        &self,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<Option<BezierConversionPlan>, NurbsError> {
        match plan_bezier_conversion_with_poll(*self, should_cancel)? {
            CurveWorkRun::Complete(plan) => Ok(Some(plan)),
            CurveWorkRun::Cancelled => Ok(None),
        }
    }

    /// Convert this admitted generation while sharing a compound caller's
    /// cancellation callback across planning and every derived generation.
    pub(crate) fn to_bezier_form_with_poll(
        &self,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveBezierRun<S, DIM>, NurbsError> {
        self.inner
            .to_bezier_form_after_validation_with_poll(should_cancel)
    }

    /// Elevate this admitted generation exactly by one degree without
    /// repeating structural validation of the unchanged source.
    ///
    /// # Errors
    /// Propagates checked work, retained-memory, allocation, numeric-domain,
    /// and structural refusals from exact Bezier conversion and elevation.
    pub fn elevate_degree(&self) -> Result<NurbsCurve<S, DIM>, NurbsError> {
        let mut never_cancel = || false;
        match self.elevate_degree_with_poll(&mut never_cancel)? {
            CurveElevationRun::Complete { curve } => Ok(curve),
            CurveElevationRun::Cancelled => Err(NurbsError::Domain {
                what: "non-cancelling degree elevation observed cancellation".to_string(),
            }),
        }
    }

    /// Elevate this admitted generation exactly by one degree with bounded
    /// cancellation polling and transactional publication.
    ///
    /// Checked conversion and elevation work/retained-memory refusals precede
    /// elevation allocation. One gate then spans exact Bezier conversion,
    /// metadata and output allocation, knot-run and span traversal, four-lane
    /// binomial blends, knot replication, both derived validation passes, and
    /// final publication. Cancellation drops all partial derived storage. The
    /// caller remains responsible for owning admission, `Cx` budget
    /// consumption, and request -> drain -> finalize semantics.
    ///
    /// # Errors
    /// Returns the synchronous elevation's checked work, retained-memory,
    /// allocation, numeric-domain, or structural refusal when it wins before
    /// an observed cancellation.
    pub fn elevate_degree_with_cx(
        &self,
        cx: &Cx<'_>,
    ) -> Result<CurveElevationRun<S, DIM>, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        self.elevate_degree_with_poll(&mut should_cancel)
    }

    /// Elevate an admitted curve while sharing a compound caller's
    /// cancellation callback across conversion and derived assembly.
    pub(crate) fn elevate_degree_with_poll(
        &self,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveElevationRun<S, DIM>, NurbsError> {
        self.inner
            .elevate_degree_after_validation_with_poll(should_cancel)
    }

    /// Constant-time charge required before scanning knot runs to construct a
    /// Bezier conversion plan.
    pub(crate) fn bezier_pre_scan_work(&self) -> Result<u128, NurbsError> {
        bezier_pre_scan_work(self.knots().knots().len())
    }

    /// Conservatively price direct insertions followed by Bezier conversion
    /// at the largest derived generation, without scanning or allocating.
    pub(crate) fn projected_refinement_work(
        &self,
        direct_insertions: usize,
        bezier_insertions: usize,
    ) -> Result<u128, NurbsError> {
        refinement_work_upper_bound(
            self.knots().degree(),
            self.knots().knots().len(),
            self.homogeneous_control_points().len(),
            direct_insertions,
            bezier_insertions,
        )
    }

    /// Build per-span Cartesian control boxes from this admitted generation
    /// without rescanning its sealed source structure.
    ///
    /// # Errors
    /// Returns a structured refusal when traversal work, retained bytes, or
    /// output allocation exceed the defensive envelope.
    pub fn span_boxes(&self) -> Result<Vec<SpanBox<S, DIM>>, NurbsError> {
        self.inner.span_boxes_after_validation()
    }

    /// Build the complete ordered span-box table with bounded cancellation
    /// polling and transactional publication.
    ///
    /// Checked traversal work and retained-output refusal precede cancellation.
    /// The gate then spans allocation, candidate-span traversal, Cartesian
    /// projection/bounds work, and final publication. This method does not
    /// consume the `Cx` budget or finalize its executor scope.
    ///
    /// # Errors
    /// Returns the synchronous span-box builder's work, memory, and allocation
    /// refusals when they win before an observed cancellation.
    pub fn span_boxes_with_cx(&self, cx: &Cx<'_>) -> Result<CurveSpanBoxesRun<S, DIM>, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        self.span_boxes_with_poll(&mut should_cancel)
    }

    /// Build admitted span boxes while sharing a compound caller's
    /// cancellation callback.
    pub(crate) fn span_boxes_with_poll(
        &self,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveSpanBoxesRun<S, DIM>, NurbsError> {
        self.inner
            .span_boxes_after_validation_with_poll(should_cancel)
    }
}

impl<S: Scalar, const DIM: usize> NurbsCurve<S, DIM> {
    /// EXACT Boehm knot insertion at `t` (multiplicity one per call).
    ///
    /// # Errors
    /// [`NurbsError::Domain`] outside the OPEN domain interior or when checked
    /// output sizing/allocation/validation is refused.
    pub fn insert_knot(&self, t: S) -> Result<Self, NurbsError> {
        self.admit()?.insert_knot(t)
    }

    fn insertion_plan_after_parameter(&self, t: S) -> Result<CurveInsertionPlan, NurbsError> {
        let admitted_knots = self.knots.admitted_after_validation();
        let (lo, hi) = admitted_knots.domain();
        if t <= lo || hi <= t {
            return Err(NurbsError::Domain {
                what: format!("insertion parameter {t:?} must be interior to {lo:?}..{hi:?}"),
            });
        }
        if !t.is_finite() {
            return Err(NurbsError::Domain {
                what: format!("parameter {t:?} outside {lo:?}..{hi:?}"),
            });
        }
        plan_curve_insertion(self.admitted_after_validation())
    }

    fn insert_knot_after_validation(&self, t: S) -> Result<Self, NurbsError> {
        let plan = self.insertion_plan_after_parameter(t)?;
        let admitted_knots = self.knots.admitted_after_validation();
        let k = admitted_knots.span(t)?;
        let mut never_cancel = || false;
        match self.insert_knot_at_span_with_plan_and_poll(t, k, plan, &mut never_cancel)? {
            CurveWorkRun::Complete(curve) => Ok(curve),
            CurveWorkRun::Cancelled => Err(NurbsError::Domain {
                what: "non-cancelling knot insertion observed cancellation".to_string(),
            }),
        }
    }

    fn insert_knot_at_span_with_poll(
        &self,
        t: S,
        k: usize,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveWorkRun<Self>, NurbsError> {
        let plan = self.insertion_plan_after_parameter(t)?;
        self.insert_knot_at_span_with_plan_and_poll(t, k, plan, should_cancel)
    }

    fn insert_knot_at_span_with_plan_and_poll(
        &self,
        t: S,
        k: usize,
        plan: CurveInsertionPlan,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveWorkRun<Self>, NurbsError> {
        let p = self.knots.degree;
        if should_cancel() {
            return Ok(CurveWorkRun::Cancelled);
        }
        let mut new_cpw = Vec::new();
        new_cpw
            .try_reserve_exact(plan.new_control_count)
            .map_err(|_| NurbsError::Domain {
                what: "inserted curve-control allocation was refused".to_string(),
            })?;
        let mut operations_since_poll = 0usize;
        for &control in &self.cpw[..=k - p] {
            new_cpw.push(control);
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveWorkRun::Cancelled);
            }
        }
        for i in (k - p + 1)..=k {
            let denom = self.knots.knots[i + p] - self.knots.knots[i];
            let alpha = (t - self.knots.knots[i]) / denom;
            let mut blended = [S::zero(); 4];
            for ((slot, &left), &right) in
                blended.iter_mut().zip(&self.cpw[i - 1]).zip(&self.cpw[i])
            {
                *slot = (S::one() - alpha) * left + alpha * right;
                if curve_poll_due(&mut operations_since_poll, should_cancel) {
                    return Ok(CurveWorkRun::Cancelled);
                }
            }
            new_cpw.push(blended);
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveWorkRun::Cancelled);
            }
        }
        for &control in &self.cpw[k..] {
            new_cpw.push(control);
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveWorkRun::Cancelled);
            }
        }

        if should_cancel() {
            return Ok(CurveWorkRun::Cancelled);
        }
        let mut new_knots = Vec::new();
        new_knots
            .try_reserve_exact(plan.new_knot_count)
            .map_err(|_| NurbsError::Domain {
                what: "inserted knot-vector allocation was refused".to_string(),
            })?;
        operations_since_poll = 0;
        for &knot in &self.knots.knots[..=k] {
            new_knots.push(knot);
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveWorkRun::Cancelled);
            }
        }
        new_knots.push(t);
        if curve_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(CurveWorkRun::Cancelled);
        }
        for &knot in &self.knots.knots[k + 1..] {
            new_knots.push(knot);
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveWorkRun::Cancelled);
            }
        }

        let knots = KnotVector {
            knots: new_knots,
            degree: p,
        };
        KnotVector::<S>::enforce_work(knots.validation_work()?, "knot-vector construction")?;
        match knots.validate_live_with_poll(|| should_cancel())? {
            KnotValidationOutcome::Complete => {}
            KnotValidationOutcome::Cancelled => return Ok(CurveWorkRun::Cancelled),
        }
        let candidate = NurbsCurve {
            knots,
            cpw: new_cpw,
        };
        match candidate.validate_live_structure_with_poll(should_cancel)? {
            CurveWorkRun::Complete(()) => {}
            CurveWorkRun::Cancelled => return Ok(CurveWorkRun::Cancelled),
        }
        if should_cancel() {
            return Ok(CurveWorkRun::Cancelled);
        }
        Ok(CurveWorkRun::Complete(candidate))
    }

    /// EXACT knot removal (inverse of [`Self::insert_knot`]) — succeeds
    /// only when the curve is exactly representable without the knot
    /// (e.g. a knot that was previously inserted); the reconstruction's
    /// consistency equation is checked with SCALAR EQUALITY, so in `Rat`
    /// this is a proof, not a tolerance.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] when `t` is not an interior knot;
    /// [`NurbsError::Structure`] when removal is not exact.
    pub fn remove_knot(&self, t: S) -> Result<Self, NurbsError> {
        self.admit()?.remove_knot(t)
    }

    fn removal_plan_after_parameter(&self, t: S) -> Result<CurveRemovalPlan, NurbsError> {
        let admitted = self.admitted_after_validation();
        let (lo, hi) = admitted.knots().domain();
        if !t.is_finite() || t <= lo || hi <= t {
            return Err(NurbsError::Domain {
                what: format!("{t:?} is not an interior knot"),
            });
        }
        plan_curve_removal(admitted)
    }

    fn find_removal_site_with_poll(
        &self,
        t: S,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveWorkRun<(usize, usize)>, NurbsError> {
        if should_cancel() {
            return Ok(CurveWorkRun::Cancelled);
        }
        let mut last_occurrence = None;
        let mut multiplicity = 0usize;
        let mut operations_since_poll = 0usize;
        for (index, &knot) in self.knots.knots.iter().enumerate() {
            if knot == t {
                last_occurrence = Some(index);
                multiplicity =
                    multiplicity
                        .checked_add(1)
                        .ok_or_else(|| NurbsError::Structure {
                            what: "knot-removal multiplicity overflows usize".to_string(),
                        })?;
            }
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveWorkRun::Cancelled);
            }
        }
        let Some(last_occurrence) = last_occurrence else {
            return Err(NurbsError::Domain {
                what: format!("{t:?} is not an interior knot"),
            });
        };
        let prior_multiplicity =
            multiplicity
                .checked_sub(1)
                .ok_or_else(|| NurbsError::Structure {
                    what: "knot-removal occurrence accounting underflowed".to_string(),
                })?;
        if should_cancel() {
            return Ok(CurveWorkRun::Cancelled);
        }
        Ok(CurveWorkRun::Complete((
            last_occurrence,
            prior_multiplicity,
        )))
    }

    fn copy_removed_knots_with_poll(
        &self,
        removed_index: usize,
        plan: CurveRemovalPlan,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveWorkRun<Vec<S>>, NurbsError> {
        if should_cancel() {
            return Ok(CurveWorkRun::Cancelled);
        }
        let mut new_knots = Vec::new();
        new_knots
            .try_reserve_exact(plan.new_knot_count)
            .map_err(|_| NurbsError::Domain {
                what: "knot-removal knot-vector allocation was refused".to_string(),
            })?;
        let mut operations_since_poll = 0usize;
        for (index, &knot) in self.knots.knots.iter().enumerate() {
            if index != removed_index {
                new_knots.push(knot);
            }
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveWorkRun::Cancelled);
            }
        }
        if new_knots.len() != plan.new_knot_count {
            return Err(NurbsError::Structure {
                what: "knot-removal knot copy produced the wrong length".to_string(),
            });
        }
        if should_cancel() {
            return Ok(CurveWorkRun::Cancelled);
        }
        Ok(CurveWorkRun::Complete(new_knots))
    }

    fn reconstruct_removed_controls_with_poll(
        &self,
        t: S,
        removed_index: usize,
        prior_multiplicity: usize,
        new_knots: &[S],
        plan: CurveRemovalPlan,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveWorkRun<Vec<[S; 4]>>, NurbsError> {
        let p = self.knots.degree;
        let range = curve_removal_range(
            p,
            removed_index,
            prior_multiplicity,
            plan.reconstruction_capacity,
        )?;
        if should_cancel() {
            return Ok(CurveWorkRun::Cancelled);
        }
        let mut forward = Vec::new();
        forward
            .try_reserve_exact(range.forward_count)
            .map_err(|_| NurbsError::Domain {
                what: "knot-removal reconstruction allocation was refused".to_string(),
            })?;
        let mut previous =
            *self
                .cpw
                .get(range.left_control)
                .ok_or_else(|| NurbsError::Structure {
                    what: "knot-removal left control is outside the source net".to_string(),
                })?;
        forward.push(previous);

        let mut operations_since_poll = 0usize;
        for index in range.blend_start..=range.blend_end {
            let right_knot = index.checked_add(p).ok_or_else(|| NurbsError::Structure {
                what: "knot-removal denominator index overflows usize".to_string(),
            })?;
            let denominator = *new_knots
                .get(right_knot)
                .ok_or_else(|| NurbsError::Structure {
                    what: "knot-removal denominator exceeds the derived knots".to_string(),
                })?
                - *new_knots.get(index).ok_or_else(|| NurbsError::Structure {
                    what: "knot-removal blend index exceeds the derived knots".to_string(),
                })?;
            let alpha = (t - new_knots[index]) / denominator;
            if alpha == S::zero() {
                return Err(NurbsError::Structure {
                    what: "degenerate removal alpha".to_string(),
                });
            }
            let source = self.cpw.get(index).ok_or_else(|| NurbsError::Structure {
                what: "knot-removal blend control exceeds the source net".to_string(),
            })?;
            let mut reconstructed = [S::zero(); 4];
            for ((slot, &value), &prior) in reconstructed.iter_mut().zip(source).zip(&previous) {
                *slot = (value - (S::one() - alpha) * prior) / alpha;
                if curve_poll_due(&mut operations_since_poll, should_cancel) {
                    return Ok(CurveWorkRun::Cancelled);
                }
            }
            forward.push(reconstructed);
            previous = reconstructed;
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveWorkRun::Cancelled);
            }
        }

        if forward.last() != self.cpw.get(range.suffix_start) {
            return Err(NurbsError::Structure {
                what: "knot is not exactly removable (curve genuinely uses it)".to_string(),
            });
        }
        self.assemble_removed_controls_with_poll(
            range.left_control,
            range.suffix_start,
            &forward,
            plan,
            should_cancel,
        )
    }

    fn assemble_removed_controls_with_poll(
        &self,
        left_control: usize,
        suffix_start: usize,
        forward: &[[S; 4]],
        plan: CurveRemovalPlan,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveWorkRun<Vec<[S; 4]>>, NurbsError> {
        if should_cancel() {
            return Ok(CurveWorkRun::Cancelled);
        }
        let mut new_cpw = Vec::new();
        new_cpw
            .try_reserve_exact(plan.new_control_count)
            .map_err(|_| NurbsError::Domain {
                what: "knot-removal control-net allocation was refused".to_string(),
            })?;
        let mut operations_since_poll = 0usize;
        for &control in &self.cpw[..left_control] {
            new_cpw.push(control);
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveWorkRun::Cancelled);
            }
        }
        for &control in &forward[..forward.len() - 1] {
            new_cpw.push(control);
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveWorkRun::Cancelled);
            }
        }
        for &control in &self.cpw[suffix_start..] {
            new_cpw.push(control);
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveWorkRun::Cancelled);
            }
        }
        if new_cpw.len() != plan.new_control_count {
            return Err(NurbsError::Structure {
                what: "knot-removal control reconstruction produced the wrong length".to_string(),
            });
        }
        if should_cancel() {
            return Ok(CurveWorkRun::Cancelled);
        }
        Ok(CurveWorkRun::Complete(new_cpw))
    }

    fn representation_matches_with_poll(
        &self,
        other: &Self,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> CurveWorkRun<bool> {
        if self.knots.degree != other.knots.degree
            || self.knots.knots.len() != other.knots.knots.len()
            || self.cpw.len() != other.cpw.len()
        {
            return CurveWorkRun::Complete(false);
        }
        if should_cancel() {
            return CurveWorkRun::Cancelled;
        }
        let mut operations_since_poll = 0usize;
        for (&left, &right) in self.knots.knots.iter().zip(&other.knots.knots) {
            if left != right {
                return CurveWorkRun::Complete(false);
            }
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return CurveWorkRun::Cancelled;
            }
        }
        for (left, right) in self.cpw.iter().zip(&other.cpw) {
            for (&left, &right) in left.iter().zip(right) {
                if left != right {
                    return CurveWorkRun::Complete(false);
                }
                if curve_poll_due(&mut operations_since_poll, should_cancel) {
                    return CurveWorkRun::Cancelled;
                }
            }
        }
        if should_cancel() {
            return CurveWorkRun::Cancelled;
        }
        CurveWorkRun::Complete(true)
    }

    fn remove_knot_after_validation_with_poll(
        &self,
        t: S,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveRemovalRun<S, DIM>, NurbsError> {
        let plan = self.removal_plan_after_parameter(t)?;
        let (removed_index, prior_multiplicity) =
            match self.find_removal_site_with_poll(t, should_cancel)? {
                CurveWorkRun::Complete(site) => site,
                CurveWorkRun::Cancelled => return Ok(CurveRemovalRun::Cancelled),
            };
        let new_knots =
            match self.copy_removed_knots_with_poll(removed_index, plan, should_cancel)? {
                CurveWorkRun::Complete(knots) => knots,
                CurveWorkRun::Cancelled => return Ok(CurveRemovalRun::Cancelled),
            };
        let new_cpw = match self.reconstruct_removed_controls_with_poll(
            t,
            removed_index,
            prior_multiplicity,
            &new_knots,
            plan,
            should_cancel,
        )? {
            CurveWorkRun::Complete(controls) => controls,
            CurveWorkRun::Cancelled => return Ok(CurveRemovalRun::Cancelled),
        };

        let knots = KnotVector {
            knots: new_knots,
            degree: self.knots.degree,
        };
        KnotVector::<S>::enforce_work(knots.validation_work()?, "knot-vector construction")?;
        match knots.validate_live_with_poll(|| should_cancel())? {
            KnotValidationOutcome::Complete => {}
            KnotValidationOutcome::Cancelled => return Ok(CurveRemovalRun::Cancelled),
        }
        let candidate = NurbsCurve {
            knots,
            cpw: new_cpw,
        };
        match candidate.validate_live_structure_with_poll(should_cancel)? {
            CurveWorkRun::Complete(()) => {}
            CurveWorkRun::Cancelled => return Ok(CurveRemovalRun::Cancelled),
        }

        let verification_plan = candidate.insertion_plan_after_parameter(t)?;
        if verification_plan != plan.verification_insertion {
            return Err(NurbsError::Structure {
                what: "knot-removal reinsertion shape disagrees with its admitted plan".to_string(),
            });
        }
        let span = match candidate
            .knots
            .admitted_after_validation()
            .span_with_poll(t, should_cancel)?
        {
            KnotSpanRun::Complete { span } => span,
            KnotSpanRun::Cancelled => return Ok(CurveRemovalRun::Cancelled),
        };
        let restored = match candidate.insert_knot_at_span_with_plan_and_poll(
            t,
            span,
            verification_plan,
            should_cancel,
        )? {
            CurveWorkRun::Complete(restored) => restored,
            CurveWorkRun::Cancelled => return Ok(CurveRemovalRun::Cancelled),
        };
        let matches = self.representation_matches_with_poll(&restored, should_cancel);
        drop(restored);
        match matches {
            CurveWorkRun::Complete(true) => {}
            CurveWorkRun::Complete(false) => {
                return Err(NurbsError::Structure {
                    what: "knot-removal candidate failed exact reinsertion verification"
                        .to_string(),
                });
            }
            CurveWorkRun::Cancelled => return Ok(CurveRemovalRun::Cancelled),
        }
        if should_cancel() {
            return Ok(CurveRemovalRun::Cancelled);
        }
        Ok(CurveRemovalRun::Complete { curve: candidate })
    }

    /// Decompose into Bézier segments by raising every interior knot to
    /// multiplicity `degree` (EXACT). Returns the refined curve.
    ///
    /// # Errors
    /// Propagates checked work, retained-memory, allocation, numeric-domain,
    /// and structural refusals from exact knot insertion.
    pub fn to_bezier_form(&self) -> Result<Self, NurbsError> {
        self.admit()?.to_bezier_form()
    }

    /// Validate this owning generation and convert it to exact Bezier form
    /// with one cancellation gate.
    ///
    /// Dimension and checked structural-validation work refusals precede the
    /// first checkpoint. Cancellation then spans structural admission and the
    /// same admitted conversion pipeline as
    /// [`AdmittedNurbsCurve::to_bezier_form_with_cx`]. No partial admitted
    /// authority or derived generation is published. This primitive does not
    /// consume the `Cx` budget or finalize its executor scope.
    ///
    /// # Errors
    /// Returns the synchronous owning conversion's structure, work,
    /// retained-memory, allocation, and finite-arithmetic refusals when they
    /// win before an observed cancellation.
    pub fn to_bezier_form_with_cx(
        &self,
        cx: &Cx<'_>,
    ) -> Result<CurveBezierRun<S, DIM>, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        self.to_bezier_form_with_poll(&mut should_cancel)
    }

    fn to_bezier_form_with_poll(
        &self,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveBezierRun<S, DIM>, NurbsError> {
        match self.validate_live_structure_with_poll(should_cancel)? {
            CurveWorkRun::Complete(()) => self
                .admitted_after_validation()
                .to_bezier_form_with_poll(should_cancel),
            CurveWorkRun::Cancelled => Ok(CurveBezierRun::Cancelled),
        }
    }

    fn to_bezier_form_after_validation(&self) -> Result<Self, NurbsError> {
        let mut never_cancel = || false;
        match self.to_bezier_form_after_validation_with_poll(&mut never_cancel)? {
            CurveBezierRun::Complete { curve } => Ok(curve),
            CurveBezierRun::Cancelled => Err(NurbsError::Domain {
                what: "non-cancelling Bezier conversion observed cancellation".to_string(),
            }),
        }
    }

    fn next_bezier_target_with_poll(
        &self,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> CurveWorkRun<Option<(S, usize)>> {
        if should_cancel() {
            return CurveWorkRun::Cancelled;
        }
        let p = self.knots.degree;
        let admitted_knots = self.knots.admitted_after_validation();
        let (lo, hi) = admitted_knots.domain();
        let entries = admitted_knots.knots();
        let mut operations_since_poll = 0usize;
        let mut run_start = 0usize;
        while run_start < entries.len() {
            let t = entries[run_start];
            let mut run_end = run_start + 1;
            while run_end < entries.len() && entries[run_end] == t {
                run_end += 1;
                if curve_poll_due(&mut operations_since_poll, should_cancel) {
                    return CurveWorkRun::Cancelled;
                }
            }
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return CurveWorkRun::Cancelled;
            }
            if t > lo && t < hi && run_end - run_start < p {
                if should_cancel() {
                    return CurveWorkRun::Cancelled;
                }
                return CurveWorkRun::Complete(Some((t, run_end - 1)));
            }
            run_start = run_end;
        }
        if should_cancel() {
            return CurveWorkRun::Cancelled;
        }
        CurveWorkRun::Complete(None)
    }

    fn to_bezier_form_after_validation_with_poll(
        &self,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveBezierRun<S, DIM>, NurbsError> {
        let plan = match plan_bezier_conversion_with_poll(
            self.admitted_after_validation(),
            should_cancel,
        )? {
            CurveWorkRun::Complete(plan) => plan,
            CurveWorkRun::Cancelled => return Ok(CurveBezierRun::Cancelled),
        };
        self.to_bezier_form_with_plan_and_poll(plan, should_cancel)
    }

    fn to_bezier_form_with_plan_and_poll(
        &self,
        plan: BezierConversionPlan,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveBezierRun<S, DIM>, NurbsError> {
        let mut current = match self.try_clone_with_poll(should_cancel)? {
            CurveWorkRun::Complete(curve) => curve,
            CurveWorkRun::Cancelled => return Ok(CurveBezierRun::Cancelled),
        };
        let mut completed_insertions = 0usize;
        loop {
            let target = match current.next_bezier_target_with_poll(should_cancel) {
                CurveWorkRun::Complete(target) => target,
                CurveWorkRun::Cancelled => return Ok(CurveBezierRun::Cancelled),
            };
            match target {
                Some((t, span)) => {
                    if completed_insertions >= plan.insertions {
                        return Err(NurbsError::Structure {
                            what: "Bezier conversion exceeded its checked insertion plan"
                                .to_string(),
                        });
                    }
                    current = match current.insert_knot_at_span_with_poll(t, span, should_cancel)? {
                        CurveWorkRun::Complete(curve) => curve,
                        CurveWorkRun::Cancelled => return Ok(CurveBezierRun::Cancelled),
                    };
                    completed_insertions =
                        completed_insertions.checked_add(1).ok_or_else(|| {
                            NurbsError::Structure {
                                what: "Bezier conversion insertion traversal overflowed usize"
                                    .to_string(),
                            }
                        })?;
                }
                None => {
                    if completed_insertions != plan.insertions
                        || current.knots.knots.len() != plan.final_knot_count
                        || current.cpw.len() != plan.final_control_count
                    {
                        return Err(NurbsError::Structure {
                            what: "Bezier conversion disagreed with its checked plan".to_string(),
                        });
                    }
                    if should_cancel() {
                        return Ok(CurveBezierRun::Cancelled);
                    }
                    return Ok(CurveBezierRun::Complete { curve: current });
                }
            }
        }
    }

    /// EXACT degree elevation by one: decompose to Bézier form, elevate
    /// each segment with the exact binomial rule, and reassemble with a
    /// full-multiplicity knot vector. Evaluation is IDENTICAL (the
    /// conformance suite proves it with rational equality).
    ///
    /// # Errors
    /// Propagates checked work, retained-memory, allocation, numeric-domain,
    /// and structural refusals. The borrowed source generation is excluded
    /// from the 64 MiB simultaneously-live derived-payload envelope.
    pub fn elevate_degree(&self) -> Result<Self, NurbsError> {
        self.admit()?.elevate_degree()
    }

    fn elevate_degree_after_validation_with_poll(
        &self,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveElevationRun<S, DIM>, NurbsError> {
        let admitted = self.admitted_after_validation();
        let p = admitted.knots().degree();
        let bezier_plan = match plan_bezier_conversion_with_poll(admitted, should_cancel)? {
            CurveWorkRun::Complete(plan) => plan,
            CurveWorkRun::Cancelled => return Ok(CurveElevationRun::Cancelled),
        };
        let plan = plan_curve_elevation::<S>(p, bezier_plan)?;
        let bez = match self.to_bezier_form_with_plan_and_poll(bezier_plan, should_cancel)? {
            CurveBezierRun::Complete { curve } => curve,
            CurveBezierRun::Cancelled => return Ok(CurveElevationRun::Cancelled),
        };
        if bez.knots.knots.len() != bezier_plan.final_knot_count
            || bez.cpw.len() != bezier_plan.final_control_count
        {
            return Err(NurbsError::Structure {
                what: "degree-elevation Bezier conversion disagreed with its checked plan"
                    .to_string(),
            });
        }

        // Collect distinct knots and their multiplicities in order. Ordinary
        // Bezier-form joins have multiplicity p and share one endpoint; a
        // legal full break has multiplicity p+1 and owns two independent
        // endpoints. Elevation must preserve that distinction.
        if should_cancel() {
            return Ok(CurveElevationRun::Cancelled);
        }
        let mut breaks: Vec<S> = Vec::new();
        breaks
            .try_reserve_exact(plan.distinct_knot_count)
            .map_err(|_| NurbsError::Domain {
                what: "degree-elevation break-table allocation was refused".to_string(),
            })?;
        if should_cancel() {
            return Ok(CurveElevationRun::Cancelled);
        }
        let mut multiplicities: Vec<usize> = Vec::new();
        multiplicities
            .try_reserve_exact(plan.distinct_knot_count)
            .map_err(|_| NurbsError::Domain {
                what: "degree-elevation multiplicity-table allocation was refused".to_string(),
            })?;
        let mut operations_since_poll = 0usize;
        for &u in &bez.knots.knots {
            if breaks.last() != Some(&u) {
                push_curve_elevation_value(
                    &mut breaks,
                    u,
                    plan.distinct_knot_count,
                    "break table",
                )?;
                push_curve_elevation_value(
                    &mut multiplicities,
                    1,
                    plan.distinct_knot_count,
                    "multiplicity table",
                )?;
            } else if let Some(multiplicity) = multiplicities.last_mut() {
                *multiplicity =
                    multiplicity
                        .checked_add(1)
                        .ok_or_else(|| NurbsError::Structure {
                            what: "degree-elevation knot multiplicity overflowed usize".to_string(),
                        })?;
            }
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveElevationRun::Cancelled);
            }
        }
        if breaks.len() != plan.distinct_knot_count
            || multiplicities.len() != plan.distinct_knot_count
        {
            return Err(NurbsError::Structure {
                what: "degree-elevation knot runs disagreed with their checked plan".to_string(),
            });
        }
        if should_cancel() {
            return Ok(CurveElevationRun::Cancelled);
        }

        // Elevate each Bézier segment: Q_0 = P_0; Q_{p+1} = P_p;
        // Q_i = (i/(p+1)) P_{i-1} + (1 - i/(p+1)) P_i.
        let mut new_cpw: Vec<[S; 4]> = Vec::new();
        new_cpw
            .try_reserve_exact(plan.final_control_count)
            .map_err(|_| NurbsError::Domain {
                what: "degree-elevation control-net allocation was refused".to_string(),
            })?;
        if should_cancel() {
            return Ok(CurveElevationRun::Cancelled);
        }
        let denominator = S::from_int(plan.elevated_degree_i64);
        let mut segment = 0usize;
        operations_since_poll = 0;
        for span in p..bez.knots.control_count() {
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveElevationRun::Cancelled);
            }
            if bez.knots.knots[span] >= bez.knots.knots[span + 1] {
                continue;
            }
            if segment >= plan.segment_count {
                return Err(NurbsError::Structure {
                    what: "degree elevation found more nonempty spans than its checked plan"
                        .to_string(),
                });
            }
            let pts = &bez.cpw[span - p..=span];
            let include_left_endpoint = if segment == 0 {
                true
            } else {
                let input_join_multiplicity =
                    *multiplicities
                        .get(segment)
                        .ok_or_else(|| NurbsError::Structure {
                            what: "degree elevation could not pair a span with its left knot run"
                                .to_string(),
                        })?;
                match input_join_multiplicity {
                    m if m == p => false,
                    m if m == plan.elevated_degree => true,
                    m => {
                        return Err(NurbsError::Structure {
                            what: format!(
                                "Bezier-form join multiplicity {m} is neither degree {p} nor full break {}",
                                plan.elevated_degree
                            ),
                        });
                    }
                }
            };
            if include_left_endpoint {
                push_curve_elevation_value(
                    &mut new_cpw,
                    pts[0],
                    plan.final_control_count,
                    "control net",
                )?;
                if curve_poll_due(&mut operations_since_poll, should_cancel) {
                    return Ok(CurveElevationRun::Cancelled);
                }
            }
            for i in 1..=p {
                let numerator = i64::try_from(i).map_err(|_| NurbsError::Structure {
                    what: "degree elevation exceeds the scalar integer-lift domain".to_string(),
                })?;
                let alpha = S::from_int(numerator) / denominator;
                let mut v = [S::zero(); 4];
                for ((slot, &a), &b) in v.iter_mut().zip(&pts[i - 1]).zip(&pts[i]) {
                    *slot = alpha * a + (S::one() - alpha) * b;
                    if curve_poll_due(&mut operations_since_poll, should_cancel) {
                        return Ok(CurveElevationRun::Cancelled);
                    }
                }
                push_curve_elevation_value(
                    &mut new_cpw,
                    v,
                    plan.final_control_count,
                    "control net",
                )?;
                if curve_poll_due(&mut operations_since_poll, should_cancel) {
                    return Ok(CurveElevationRun::Cancelled);
                }
            }
            push_curve_elevation_value(
                &mut new_cpw,
                pts[p],
                plan.final_control_count,
                "control net",
            )?;
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveElevationRun::Cancelled);
            }
            segment = segment
                .checked_add(1)
                .ok_or_else(|| NurbsError::Structure {
                    what: "degree-elevation segment traversal overflowed usize".to_string(),
                })?;
        }
        if segment != plan.segment_count || new_cpw.len() != plan.final_control_count {
            return Err(NurbsError::Structure {
                what: "degree elevation could not pair every distinct knot interval with one nonempty span"
                    .to_string(),
            });
        }
        if should_cancel() {
            return Ok(CurveElevationRun::Cancelled);
        }

        // Elevation raises every multiplicity by one, preserving continuity
        // order. Endpoints therefore have p+2 copies, C0 joins p+1, and full
        // discontinuities p+2.
        let mut new_knots = Vec::new();
        new_knots
            .try_reserve_exact(plan.final_knot_count)
            .map_err(|_| NurbsError::Domain {
                what: "degree-elevation knot allocation was refused".to_string(),
            })?;
        if should_cancel() {
            return Ok(CurveElevationRun::Cancelled);
        }
        operations_since_poll = 0;
        for (bi, (&b, &old_multiplicity)) in breaks.iter().zip(multiplicities.iter()).enumerate() {
            let mult = if bi == 0 || bi == breaks.len() - 1 {
                plan.elevated_order
            } else {
                old_multiplicity
                    .checked_add(1)
                    .ok_or_else(|| NurbsError::Structure {
                        what: "degree-elevation knot multiplicity overflowed usize".to_string(),
                    })?
            };
            for _ in 0..mult {
                push_curve_elevation_value(
                    &mut new_knots,
                    b,
                    plan.final_knot_count,
                    "knot vector",
                )?;
                if curve_poll_due(&mut operations_since_poll, should_cancel) {
                    return Ok(CurveElevationRun::Cancelled);
                }
            }
        }
        if new_knots.len() != plan.final_knot_count {
            return Err(NurbsError::Structure {
                what: "degree-elevation knot assembly disagreed with its checked plan".to_string(),
            });
        }
        if should_cancel() {
            return Ok(CurveElevationRun::Cancelled);
        }
        let knots = KnotVector {
            knots: new_knots,
            degree: plan.elevated_degree,
        };
        KnotVector::<S>::enforce_work(
            knots.validation_work()?,
            "degree-elevation knot construction",
        )?;
        match knots.validate_live_with_poll(|| should_cancel())? {
            KnotValidationOutcome::Complete => {}
            KnotValidationOutcome::Cancelled => return Ok(CurveElevationRun::Cancelled),
        }
        let elevated = NurbsCurve {
            knots,
            cpw: new_cpw,
        };
        match elevated.validate_live_structure_with_poll(should_cancel)? {
            CurveWorkRun::Complete(()) => {}
            CurveWorkRun::Cancelled => return Ok(CurveElevationRun::Cancelled),
        }
        if should_cancel() {
            return Ok(CurveElevationRun::Cancelled);
        }
        Ok(CurveElevationRun::Complete { curve: elevated })
    }

    /// Per-span control-point bounding boxes in Cartesian space (the
    /// convex-hull property: each span's curve lies inside its box).
    /// Requires Bézier form for the tight per-segment claim; on general
    /// knot vectors the box of the span's `p+1` control points still
    /// bounds that span.
    ///
    /// # Errors
    /// Propagates structural errors and checked work, retained-memory, or
    /// allocation refusals.
    pub fn span_boxes(&self) -> Result<Vec<SpanBox<S, DIM>>, NurbsError> {
        self.admit()?.span_boxes()
    }

    /// Validate this owning curve and build its complete ordered span-box
    /// table with one cancellation gate.
    ///
    /// Dimension and checked structural-validation work refusals precede the
    /// first checkpoint. Cancellation then spans structural admission and the
    /// same admitted span-box pipeline as
    /// [`AdmittedNurbsCurve::span_boxes_with_cx`]. No partial admitted
    /// authority or box table is published. This primitive does not consume
    /// the `Cx` budget or finalize its executor scope.
    ///
    /// # Errors
    /// Returns the synchronous owning builder's structure, traversal-work,
    /// retained-memory, finite-arithmetic, and allocation refusals when they
    /// win before an observed cancellation.
    pub fn span_boxes_with_cx(&self, cx: &Cx<'_>) -> Result<CurveSpanBoxesRun<S, DIM>, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        self.span_boxes_with_poll(&mut should_cancel)
    }

    fn span_boxes_with_poll(
        &self,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveSpanBoxesRun<S, DIM>, NurbsError> {
        match self.validate_live_structure_with_poll(should_cancel)? {
            CurveWorkRun::Complete(()) => self
                .admitted_after_validation()
                .span_boxes_with_poll(should_cancel),
            CurveWorkRun::Cancelled => Ok(CurveSpanBoxesRun::Cancelled),
        }
    }

    fn span_boxes_after_validation(&self) -> Result<Vec<SpanBox<S, DIM>>, NurbsError> {
        let mut never_cancel = || false;
        match self.span_boxes_after_validation_with_poll(&mut never_cancel)? {
            CurveSpanBoxesRun::Complete { boxes } => Ok(boxes),
            CurveSpanBoxesRun::Cancelled => Err(NurbsError::Domain {
                what: "non-cancelling curve span-box traversal observed cancellation".to_string(),
            }),
        }
    }

    fn span_boxes_after_validation_with_poll(
        &self,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveSpanBoxesRun<S, DIM>, NurbsError> {
        let p = self.knots.degree;
        let span_capacity = preflight_span_boxes(
            self.knots.control_count(),
            p,
            core::mem::size_of::<SpanBox<S, DIM>>(),
        )?;
        if should_cancel() {
            return Ok(CurveSpanBoxesRun::Cancelled);
        }
        let mut out = Vec::new();
        out.try_reserve_exact(span_capacity)
            .map_err(|_| NurbsError::Domain {
                what: "curve span-box allocation was refused".to_string(),
            })?;
        let mut operations_since_poll = 0usize;
        for span in p..self.knots.control_count() {
            let (t0, t1) = (self.knots.knots[span], self.knots.knots[span + 1]);
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveSpanBoxesRun::Cancelled);
            }
            if t1 <= t0 {
                continue;
            }
            let mut min = [S::zero(); DIM];
            let mut max = [S::zero(); DIM];
            let mut first = true;
            for cp in &self.cpw[span - p..=span] {
                let w = cp[3];
                for d in 0..DIM {
                    let c = cp[d] / w;
                    if first {
                        min[d] = c;
                        max[d] = c;
                    } else {
                        if c < min[d] {
                            min[d] = c;
                        }
                        if max[d] < c {
                            max[d] = c;
                        }
                    }
                    if curve_poll_due(&mut operations_since_poll, should_cancel) {
                        return Ok(CurveSpanBoxesRun::Cancelled);
                    }
                }
                first = false;
                if curve_poll_due(&mut operations_since_poll, should_cancel) {
                    return Ok(CurveSpanBoxesRun::Cancelled);
                }
            }
            out.push((min, max, t0, t1));
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveSpanBoxesRun::Cancelled);
            }
        }
        if should_cancel() {
            return Ok(CurveSpanBoxesRun::Cancelled);
        }
        Ok(CurveSpanBoxesRun::Complete { boxes: out })
    }
}

fn evaluate_homogeneous_derivative_net_with_poll(
    net: &[[f64; 4]],
    knots: &[f64],
    degree: usize,
    t: f64,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<CurveWorkRun<[f64; 4]>, NurbsError> {
    let expected_knots = net
        .len()
        .checked_add(degree)
        .and_then(|count| count.checked_add(1))
        .ok_or_else(|| NurbsError::Structure {
            what: "derivative-net knot-count arithmetic overflowed".to_string(),
        })?;
    let malformed_net = || NurbsError::Structure {
        what: "reduced homogeneous derivative net is malformed".to_string(),
    };
    if net.is_empty() || knots.len() != expected_knots {
        return Err(malformed_net());
    }
    if should_cancel() {
        return Ok(CurveWorkRun::Cancelled);
    }
    let mut operations_since_poll = 0usize;
    for &knot in knots {
        if !knot.is_finite() {
            return Err(malformed_net());
        }
        if curve_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(CurveWorkRun::Cancelled);
        }
    }
    for pair in knots.windows(2) {
        if pair[1] < pair[0] {
            return Err(malformed_net());
        }
        if curve_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(CurveWorkRun::Cancelled);
        }
    }
    let lo = knots[degree];
    let hi = knots[knots.len() - 1 - degree];
    if !t.is_finite() || t < lo || t > hi || lo >= hi {
        return Err(NurbsError::Domain {
            what: format!("derivative parameter {t} outside {lo}..{hi}"),
        });
    }
    let last_control = net.len() - 1;
    let span = if t == hi {
        let mut upper_span = None;
        for candidate in (0..=last_control).rev() {
            let nonempty = knots[candidate] < knots[candidate + 1];
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveWorkRun::Cancelled);
            }
            if nonempty {
                upper_span = Some(candidate);
                break;
            }
        }
        let Some(span) = upper_span else {
            return Err(NurbsError::Structure {
                what: "reduced derivative net has no nonempty upper span".to_string(),
            });
        };
        span
    } else {
        let mut span = degree;
        while span < last_control && knots[span + 1] <= t {
            span += 1;
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveWorkRun::Cancelled);
            }
        }
        span
    };
    if span < degree || span - degree + degree >= net.len() {
        return Err(NurbsError::Structure {
            what: "reduced derivative span does not index its control net".to_string(),
        });
    }

    let basis_len = degree.checked_add(1).ok_or_else(|| NurbsError::Domain {
        what: "derivative basis length overflowed".to_string(),
    })?;
    let mut basis = Vec::new();
    let mut left = Vec::new();
    let mut right = Vec::new();
    for (buffer, stage) in [
        (&mut basis, "basis"),
        (&mut left, "left basis workspace"),
        (&mut right, "right basis workspace"),
    ] {
        if should_cancel() {
            return Ok(CurveWorkRun::Cancelled);
        }
        buffer
            .try_reserve_exact(basis_len)
            .map_err(|_| NurbsError::Domain {
                what: format!("derivative {stage} allocation was refused"),
            })?;
        for _ in 0..basis_len {
            buffer.push(0.0);
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveWorkRun::Cancelled);
            }
        }
    }
    basis[0] = 1.0;
    for j in 1..=degree {
        left[j] = t - knots[span + 1 - j];
        if curve_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(CurveWorkRun::Cancelled);
        }
        right[j] = knots[span + j] - t;
        if curve_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(CurveWorkRun::Cancelled);
        }
        let mut saved = 0.0;
        for r in 0..j {
            let denominator = right[r + 1] + left[j - r];
            // Cox-de Boor's zero-width-span convention is 0/0 -> 0. This is
            // essential for reduced derivative nets whose inherited interior
            // multiplicity can exceed their reduced polynomial degree; away
            // from the break, the active nonzero span remains ordinary.
            let term = if denominator == 0.0 {
                0.0
            } else {
                basis[r] / denominator
            };
            basis[r] = saved + right[r + 1] * term;
            saved = left[j - r] * term;
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveWorkRun::Cancelled);
            }
        }
        basis[j] = saved;
        if curve_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(CurveWorkRun::Cancelled);
        }
    }
    for &basis_value in &basis {
        if !basis_value.is_finite() {
            return Err(NurbsError::Domain {
                what: "reduced derivative basis left the finite numeric domain".to_string(),
            });
        }
        if curve_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(CurveWorkRun::Cancelled);
        }
    }
    let mut value = [0.0; 4];
    for (offset, weight) in basis.into_iter().enumerate() {
        for (accumulator, control) in value.iter_mut().zip(net[span - degree + offset].iter()) {
            *accumulator += weight * control;
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveWorkRun::Cancelled);
            }
        }
    }
    for &component in &value {
        if !component.is_finite() {
            return Err(NurbsError::Domain {
                what: "reduced derivative evaluation left the finite numeric domain".to_string(),
            });
        }
        if curve_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(CurveWorkRun::Cancelled);
        }
    }
    Ok(CurveWorkRun::Complete(value))
}

impl<const DIM: usize> NurbsCurve<f64, DIM> {
    /// Defensive ceiling for the allocation- and quadratic-work-bearing legacy
    /// derivative API. Budgeted high-order jets belong to the typed successor.
    const MAX_DERIVATIVE_ORDER: usize = 64;
    /// Combined retained-net, quotient-recurrence, and allocation ceiling for
    /// the legacy whole-net derivative implementation.
    pub(crate) const MAX_DERIVATIVE_WORK_UNITS: u128 = 16_777_216;
    /// Hard retained payload bound for homogeneous nets and knot copies. Vec
    /// metadata and temporary basis arrays add a small bounded overhead.
    pub(crate) const MAX_DERIVATIVE_RETAINED_BYTES: u128 = 67_108_864;

    /// Derivatives up to `order` at `t` (rational quotient rule over the
    /// homogeneous derivative curves). Returns `[C(t), C'(t), …]`.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] outside the parameter domain or above the
    /// defensive legacy order ceiling.
    pub fn derivatives(&self, t: f64, order: usize) -> Result<Vec<[f64; DIM]>, NurbsError> {
        Self::preflight_derivative_shape(t, order)?;
        self.knots.preflight_parameter(t, "curve derivative")?;
        self.validate_live_structure()?;
        let knots = self.knots.admitted_after_validation();
        Self::preflight_derivative_request(knots, t, order)?;
        Self::derivatives_from_admitted_parts_after_preflight(knots, &self.cpw, t, order)
    }

    /// Validate the owning curve and evaluate its Cartesian jet with one
    /// cancellation gate.
    ///
    /// Constant-time shape, order, and parameter refusals precede the first
    /// checkpoint. Cancellation then spans structural admission and the same
    /// admitted derivative pipeline as [`AdmittedNurbsCurve::derivatives_with_cx`].
    /// No partial admitted authority or derivative jet is published. This
    /// primitive does not consume the `Cx` budget or finalize its executor
    /// scope.
    ///
    /// # Errors
    /// Returns the synchronous owning derivative evaluator's request,
    /// structure, work, memory, allocation, and finite-arithmetic refusals when
    /// they win before an observed cancellation.
    pub fn derivatives_with_cx(
        &self,
        t: f64,
        order: usize,
        cx: &Cx<'_>,
    ) -> Result<CurveDerivativesRun<DIM>, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        self.derivatives_with_poll(t, order, &mut should_cancel)
    }

    fn derivatives_with_poll(
        &self,
        t: f64,
        order: usize,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveDerivativesRun<DIM>, NurbsError> {
        Self::preflight_derivative_shape(t, order)?;
        self.knots.preflight_parameter(t, "curve derivative")?;
        match self.validate_live_structure_with_poll(should_cancel)? {
            CurveWorkRun::Complete(()) => {
                self.admitted_after_validation()
                    .derivatives_with_poll(t, order, should_cancel)
            }
            CurveWorkRun::Cancelled => Ok(CurveDerivativesRun::Cancelled),
        }
    }

    fn derivatives_after_validation(
        &self,
        t: f64,
        order: usize,
    ) -> Result<Vec<[f64; DIM]>, NurbsError> {
        Self::derivatives_from_admitted_parts(
            self.knots.admitted_after_validation(),
            &self.cpw,
            t,
            order,
        )
    }

    fn preflight_derivative_shape(t: f64, order: usize) -> Result<(), NurbsError> {
        if DIM > 3 {
            return Err(NurbsError::Structure {
                what: format!("curve dimension {DIM} exceeds the homogeneous storage limit 3"),
            });
        }
        if order > Self::MAX_DERIVATIVE_ORDER {
            return Err(NurbsError::Domain {
                what: format!(
                    "derivative order {order} exceeds defensive ceiling {}",
                    Self::MAX_DERIVATIVE_ORDER
                ),
            });
        }
        if !t.is_finite() {
            return Err(NurbsError::Domain {
                what: "derivative parameter must be finite".to_string(),
            });
        }
        Ok(())
    }

    pub(crate) fn derivative_envelope(
        control_count: usize,
        knot_count: usize,
        degree: usize,
        order: usize,
    ) -> Result<(u128, u128), NurbsError> {
        let homogeneous_order = order.min(degree);
        let retained_nets = (control_count as u128)
            .checked_add(knot_count as u128)
            .and_then(|extent| extent.checked_mul((homogeneous_order as u128).saturating_add(1)))
            .ok_or_else(|| NurbsError::Domain {
                what: "derivative retained-net accounting overflows u128".to_string(),
            })?;
        let quotient_extent = (order as u128)
            .checked_add(1)
            .and_then(|side| side.checked_mul(side))
            .and_then(|square| square.checked_mul((DIM as u128).saturating_add(4)))
            .ok_or_else(|| NurbsError::Domain {
                what: "derivative quotient-work accounting overflows u128".to_string(),
            })?;
        let basis_extent = (0..=homogeneous_order).try_fold(0u128, |total, derivative| {
            let reduced_degree = degree - derivative;
            let basis_order = reduced_degree
                .checked_add(1)
                .ok_or_else(|| NurbsError::Domain {
                    what: "derivative basis-order accounting overflows usize".to_string(),
                })?;
            let work = (reduced_degree as u128)
                .checked_mul(basis_order as u128)
                .map(|product| product / 2)
                .and_then(|triangular| triangular.checked_add(basis_order as u128))
                .ok_or_else(|| NurbsError::Domain {
                    what: "derivative basis-work accounting overflows u128".to_string(),
                })?;
            total.checked_add(work).ok_or_else(|| NurbsError::Domain {
                what: "derivative aggregate basis-work accounting overflows u128".to_string(),
            })
        })?;
        let requested_work = retained_nets
            .checked_add(quotient_extent)
            .and_then(|work| work.checked_add(basis_extent))
            .ok_or_else(|| NurbsError::Domain {
                what: "derivative total-work accounting overflows u128".to_string(),
            })?;

        let levels = (homogeneous_order as u128).saturating_add(1);
        let control_bytes = (control_count as u128)
            .checked_mul(levels)
            .and_then(|count| count.checked_mul(core::mem::size_of::<[f64; 4]>() as u128))
            .ok_or_else(|| NurbsError::Domain {
                what: "derivative retained-control accounting overflows u128".to_string(),
            })?;
        // Charge a full knot extent per derivative level even though the
        // implementation below borrows successively trimmed source slices.
        // The deliberate overestimate keeps the envelope stable if those
        // views later become owned again.
        let knot_bytes = (knot_count as u128)
            .checked_mul(levels)
            .and_then(|count| count.checked_mul(core::mem::size_of::<f64>() as u128))
            .ok_or_else(|| NurbsError::Domain {
                what: "derivative retained-knot accounting overflows u128".to_string(),
            })?;
        let table_bytes = levels
            .checked_mul(core::mem::size_of::<Vec<[f64; 4]>>() as u128)
            .ok_or_else(|| NurbsError::Domain {
                what: "derivative net-table accounting overflows u128".to_string(),
            })?;
        let jet_len = (order as u128)
            .checked_add(1)
            .ok_or_else(|| NurbsError::Domain {
                what: "derivative jet-length accounting overflows u128".to_string(),
            })?;
        let homogeneous_bytes = jet_len
            .checked_mul(core::mem::size_of::<[f64; 4]>() as u128)
            .ok_or_else(|| NurbsError::Domain {
                what: "derivative homogeneous-jet accounting overflows u128".to_string(),
            })?;
        let basis_bytes = (degree as u128)
            .checked_add(1)
            .and_then(|basis_len| basis_len.checked_mul(3))
            .and_then(|basis_values| basis_values.checked_mul(core::mem::size_of::<f64>() as u128))
            .ok_or_else(|| NurbsError::Domain {
                what: "derivative basis-workspace accounting overflows u128".to_string(),
            })?;
        let cartesian_bytes = jet_len
            .checked_mul(DIM as u128)
            .and_then(|count| count.checked_mul(core::mem::size_of::<f64>() as u128))
            .ok_or_else(|| NurbsError::Domain {
                what: "derivative Cartesian-jet accounting overflows u128".to_string(),
            })?;
        let retained_bytes = control_bytes
            .checked_add(knot_bytes)
            .and_then(|bytes| bytes.checked_add(table_bytes))
            .and_then(|bytes| bytes.checked_add(homogeneous_bytes))
            .and_then(|bytes| bytes.checked_add(basis_bytes.max(cartesian_bytes)))
            .ok_or_else(|| NurbsError::Domain {
                what: "derivative retained-byte accounting overflows u128".to_string(),
            })?;
        Ok((requested_work, retained_bytes))
    }

    pub(crate) fn preflight_derivative_request(
        knots: AdmittedKnotVector<'_, f64>,
        t: f64,
        order: usize,
    ) -> Result<(), NurbsError> {
        Self::preflight_derivative_shape(t, order)?;
        let (lo, hi) = knots.domain();
        if t < lo || t > hi {
            return Err(NurbsError::Domain {
                what: format!("derivative parameter {t} outside {lo}..{hi}"),
            });
        }
        let p = knots.degree();
        if t > lo && t < hi {
            let entries = knots.knots();
            let run_start = entries.partition_point(|&knot| knot < t);
            let run_end = entries.partition_point(|&knot| knot <= t);
            let multiplicity = run_end - run_start;
            if multiplicity > 0 && (multiplicity > p || order > p - multiplicity) {
                return Err(NurbsError::Domain {
                    what: format!(
                        "ordinary derivative order {order} is undefined at interior knot multiplicity {multiplicity} for degree {p}; request an explicit one-sided jet in the successor API"
                    ),
                });
            }
        }
        Ok(())
    }

    pub(crate) fn derivatives_from_admitted_parts(
        knots: AdmittedKnotVector<'_, f64>,
        cpw: &[[f64; 4]],
        t: f64,
        order: usize,
    ) -> Result<Vec<[f64; DIM]>, NurbsError> {
        Self::preflight_derivative_request(knots, t, order)?;
        Self::derivatives_from_admitted_parts_after_preflight(knots, cpw, t, order)
    }

    pub(crate) fn derivatives_from_admitted_parts_after_preflight(
        knots: AdmittedKnotVector<'_, f64>,
        cpw: &[[f64; 4]],
        t: f64,
        order: usize,
    ) -> Result<Vec<[f64; DIM]>, NurbsError> {
        let mut never_cancel = || false;
        match Self::derivatives_from_admitted_parts_after_preflight_with_poll(
            knots,
            cpw,
            t,
            order,
            &mut never_cancel,
        )? {
            CurveDerivativesRun::Complete { derivatives } => Ok(derivatives),
            CurveDerivativesRun::Cancelled => Err(NurbsError::Domain {
                what: "non-cancelling curve derivative evaluation observed cancellation"
                    .to_string(),
            }),
        }
    }

    pub(crate) fn derivatives_from_admitted_parts_after_preflight_with_poll(
        knots: AdmittedKnotVector<'_, f64>,
        cpw: &[[f64; 4]],
        t: f64,
        order: usize,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveDerivativesRun<DIM>, NurbsError> {
        if cpw.len() != knots.control_count() {
            return Err(NurbsError::Structure {
                what: "admitted derivative control count does not match its knot vector"
                    .to_string(),
            });
        }
        let p = knots.degree();
        let homogeneous_order = order.min(p);
        let (requested_work, retained_bytes) =
            Self::derivative_envelope(cpw.len(), knots.knots().len(), p, order)?;
        if retained_bytes > Self::MAX_DERIVATIVE_RETAINED_BYTES {
            return Err(NurbsError::Domain {
                what: format!(
                    "derivative request retains up to {retained_bytes} bytes, above ceiling {}",
                    Self::MAX_DERIVATIVE_RETAINED_BYTES
                ),
            });
        }
        if requested_work > Self::MAX_DERIVATIVE_WORK_UNITS {
            return Err(NurbsError::Domain {
                what: format!(
                    "derivative request needs {requested_work} defensive work units, above ceiling {}",
                    Self::MAX_DERIVATIVE_WORK_UNITS
                ),
            });
        }
        // Homogeneous derivative control nets by repeated differencing.
        if should_cancel() {
            return Ok(CurveDerivativesRun::Cancelled);
        }
        let mut nets: Vec<Vec<[f64; 4]>> = Vec::new();
        nets.try_reserve_exact(homogeneous_order + 1)
            .map_err(|_| NurbsError::Domain {
                what: "derivative net-table allocation was refused".to_string(),
            })?;
        if should_cancel() {
            return Ok(CurveDerivativesRun::Cancelled);
        }
        let mut initial_net = Vec::new();
        initial_net
            .try_reserve_exact(cpw.len())
            .map_err(|_| NurbsError::Domain {
                what: "derivative initial control-net allocation was refused".to_string(),
            })?;
        let mut operations_since_poll = 0usize;
        for &control in cpw {
            initial_net.push(control);
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveDerivativesRun::Cancelled);
            }
        }
        nets.push(initial_net);
        for k in 1..=homogeneous_order {
            let prev = &nets[k - 1];
            let degree = p - (k - 1);
            let trim = k - 1;
            let knot_end = knots.knots().len() - trim;
            let reduced_knots = &knots.knots()[trim..knot_end];
            if should_cancel() {
                return Ok(CurveDerivativesRun::Cancelled);
            }
            let mut next = Vec::new();
            next.try_reserve_exact(prev.len() - 1)
                .map_err(|_| NurbsError::Domain {
                    what: format!("derivative order {k} control-net allocation was refused"),
                })?;
            operations_since_poll = 0;
            #[allow(clippy::cast_precision_loss)]
            let degf = degree as f64;
            for i in 0..prev.len() - 1 {
                let denom = reduced_knots[i + degree + 1] - reduced_knots[i + 1];
                if curve_poll_due(&mut operations_since_poll, should_cancel) {
                    return Ok(CurveDerivativesRun::Cancelled);
                }
                let mut d = [0.0f64; 4];
                if denom != 0.0 {
                    for (slot, (a, b)) in d.iter_mut().zip(prev[i + 1].iter().zip(&prev[i])) {
                        *slot = degf * (a - b) / denom;
                        if curve_poll_due(&mut operations_since_poll, should_cancel) {
                            return Ok(CurveDerivativesRun::Cancelled);
                        }
                    }
                }
                next.push(d);
            }
            for control in &next {
                for &component in control {
                    if !component.is_finite() {
                        return Err(NurbsError::Domain {
                            what: format!(
                                "derivative order {k} control net left the finite numeric domain"
                            ),
                        });
                    }
                    if curve_poll_due(&mut operations_since_poll, should_cancel) {
                        return Ok(CurveDerivativesRun::Cancelled);
                    }
                }
            }
            nets.push(next);
        }
        // Evaluate each homogeneous derivative, then the quotient rule:
        // C^(k) = (A^(k) − Σ_{i=1..k} C(k−i) · w^(i) · binom(k,i)) / w.
        if should_cancel() {
            return Ok(CurveDerivativesRun::Cancelled);
        }
        let mut hom = Vec::new();
        hom.try_reserve_exact(order + 1)
            .map_err(|_| NurbsError::Domain {
                what: "derivative homogeneous-jet allocation was refused".to_string(),
            })?;
        for (derivative, net) in nets.iter().enumerate() {
            let knot_end = knots.knots().len() - derivative;
            let value = match evaluate_homogeneous_derivative_net_with_poll(
                net,
                &knots.knots()[derivative..knot_end],
                p - derivative,
                t,
                should_cancel,
            )? {
                CurveWorkRun::Complete(value) => value,
                CurveWorkRun::Cancelled => return Ok(CurveDerivativesRun::Cancelled),
            };
            hom.push(value);
        }
        // Polynomial homogeneous derivatives vanish above degree p, but a
        // rational quotient generally has nonzero derivatives of every order.
        // Retain those zero homogeneous jets so the quotient recurrence below
        // computes C^(k) correctly for k > p.
        operations_since_poll = 0;
        while hom.len() < order + 1 {
            hom.push([0.0; 4]);
            if curve_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(CurveDerivativesRun::Cancelled);
            }
        }
        let w0 = hom[0][3];
        if should_cancel() {
            return Ok(CurveDerivativesRun::Cancelled);
        }
        let mut out: Vec<[f64; DIM]> = Vec::new();
        out.try_reserve_exact(order + 1)
            .map_err(|_| NurbsError::Domain {
                what: "derivative Cartesian-jet allocation was refused".to_string(),
            })?;
        operations_since_poll = 0;
        for k in 0..=order {
            let mut num = [0.0f64; DIM];
            for (slot, &a) in num.iter_mut().zip(hom[k].iter()) {
                *slot = a;
                if curve_poll_due(&mut operations_since_poll, should_cancel) {
                    return Ok(CurveDerivativesRun::Cancelled);
                }
            }
            for i in 1..=k {
                let mut binomial = 1.0f64;
                for j in 0..i {
                    #[allow(clippy::cast_precision_loss)]
                    {
                        binomial = binomial * (k - j) as f64 / (j + 1) as f64;
                    }
                    if curve_poll_due(&mut operations_since_poll, should_cancel) {
                        return Ok(CurveDerivativesRun::Cancelled);
                    }
                }
                let c = binomial * hom[i][3];
                for (slot, prev) in num.iter_mut().zip(out[k - i].iter()) {
                    *slot -= c * prev;
                    if curve_poll_due(&mut operations_since_poll, should_cancel) {
                        return Ok(CurveDerivativesRun::Cancelled);
                    }
                }
            }
            let jet = num.map(|v| v / w0);
            for &component in &jet {
                if !component.is_finite() {
                    return Err(NurbsError::Domain {
                        what: format!(
                            "derivative order {k} left the finite floating-point numeric domain"
                        ),
                    });
                }
                if curve_poll_due(&mut operations_since_poll, should_cancel) {
                    return Ok(CurveDerivativesRun::Cancelled);
                }
            }
            out.push(jet);
        }
        if should_cancel() {
            return Ok(CurveDerivativesRun::Cancelled);
        }
        Ok(CurveDerivativesRun::Complete { derivatives: out })
    }
}

impl<const DIM: usize> AdmittedNurbsCurve<'_, f64, DIM> {
    /// Derivatives of the admitted curve without rescanning its source
    /// structure.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] outside the parameter domain or above the
    /// defensive legacy order/work/allocation ceilings.
    pub fn derivatives(&self, t: f64, order: usize) -> Result<Vec<[f64; DIM]>, NurbsError> {
        self.inner.derivatives_after_validation(t, order)
    }

    /// Evaluate the admitted curve's Cartesian jet with bounded cancellation
    /// polling and transactional publication.
    ///
    /// Shape, parameter, continuity, checked work, and retained-memory
    /// refusals retain their synchronous precedence. Cancellation then spans
    /// derivative-net construction, reduced basis evaluations, the rational
    /// quotient recurrence, and final publication. This method does not
    /// consume the `Cx` budget or finalize its executor scope.
    ///
    /// # Errors
    /// Returns the synchronous derivative evaluator's request, work, memory,
    /// allocation, and finite-arithmetic refusals when they win before an
    /// observed cancellation.
    pub fn derivatives_with_cx(
        &self,
        t: f64,
        order: usize,
        cx: &Cx<'_>,
    ) -> Result<CurveDerivativesRun<DIM>, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        self.derivatives_with_poll(t, order, &mut should_cancel)
    }

    /// Evaluate an admitted Cartesian jet while sharing a compound caller's
    /// cancellation callback.
    pub(crate) fn derivatives_with_poll(
        &self,
        t: f64,
        order: usize,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<CurveDerivativesRun<DIM>, NurbsError> {
        let knots = self.knots();
        NurbsCurve::<f64, DIM>::preflight_derivative_request(knots, t, order)?;
        NurbsCurve::<f64, DIM>::derivatives_from_admitted_parts_after_preflight_with_poll(
            knots,
            &self.inner.cpw,
            t,
            order,
            should_cancel,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::rat::Rat;
    use asupersync::types::Budget;
    use fs_exec::{CancelGate, ExecMode, StreamKey};

    fn with_curve_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
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
                    seed: 0xC0A7_E001,
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

    fn line_curve() -> NurbsCurve<f64, 1> {
        let knots = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("line knots");
        NurbsCurve::new(knots, &[[0.0], [1.0]], &[1.0, 1.0]).expect("line curve")
    }

    fn long_linear_curve() -> NurbsCurve<f64, 1> {
        let interior_count = 128usize;
        let mut knots = Vec::with_capacity(interior_count + 4);
        knots.extend([0.0, 0.0]);
        for index in 1..=interior_count {
            knots.push(index as f64 / (interior_count + 1) as f64);
        }
        knots.extend([1.0, 1.0]);
        let control_count = interior_count + 2;
        let points: Vec<[f64; 1]> = (0..control_count)
            .map(|index| [index as f64 / (control_count - 1) as f64])
            .collect();
        let weights = vec![1.0; control_count];
        NurbsCurve::new(
            KnotVector::new(knots, 1).expect("long linear knots"),
            &points,
            &weights,
        )
        .expect("long linear curve")
    }

    fn high_degree_insertion_curve() -> NurbsCurve<f64, 1> {
        let degree = 16usize;
        let mut knots = vec![0.0; degree + 1];
        knots.push(0.5);
        knots.extend(vec![1.0; degree + 1]);
        let control_count = degree + 2;
        let points: Vec<[f64; 1]> = (0..control_count)
            .map(|index| [index as f64 / (control_count - 1) as f64])
            .collect();
        let weights = vec![1.0; control_count];
        NurbsCurve::new(
            KnotVector::new(knots, degree).expect("high-degree insertion knots"),
            &points,
            &weights,
        )
        .expect("high-degree insertion curve")
    }

    fn quadratic_join_curve() -> NurbsCurve<f64, 1> {
        NurbsCurve::new(
            KnotVector::new(vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0], 2)
                .expect("quadratic join knots"),
            &[[0.0], [0.25], [0.75], [1.0]],
            &[1.0; 4],
        )
        .expect("quadratic join curve")
    }

    fn linear_full_break_curve() -> NurbsCurve<f64, 1> {
        NurbsCurve::new(
            KnotVector::new(vec![0.0, 0.0, 0.5, 0.5, 1.0, 1.0], 1)
                .expect("linear full-break knots"),
            &[[0.0], [0.25], [0.75], [1.0]],
            &[1.0; 4],
        )
        .expect("linear full-break curve")
    }

    #[test]
    fn cartesian_curve_construction_with_cx_is_transactional_and_exact() {
        let points = [[0.0], [1.0]];
        let weights = [1.0, 2.0];
        let expected = NurbsCurve::<f64, 1>::new(
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("line knots"),
            &points,
            &weights,
        )
        .expect("weighted line");
        with_curve_cx(true, |cx| {
            assert_eq!(
                NurbsCurve::<f64, 1>::new_with_cx(
                    KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("line knots"),
                    &points,
                    &weights,
                    cx,
                )
                .expect("valid pre-cancelled construction"),
                CurveConstructionRun::Cancelled
            );
        });
        with_curve_cx(false, |cx| {
            assert_eq!(
                NurbsCurve::<f64, 1>::new_with_cx(
                    KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("line knots"),
                    &points,
                    &weights,
                    cx,
                )
                .expect("active Cartesian construction"),
                CurveConstructionRun::Complete {
                    curve: expected.try_clone().expect("expected curve copy"),
                }
            );
        });

        let zero = Rat::int(0);
        let one = Rat::int(1);
        let two = Rat::int(2);
        let exact_points = [[zero], [one]];
        let exact_weights = [one, two];
        let exact = NurbsCurve::<Rat, 1>::new(
            KnotVector::new(vec![zero, zero, one, one], 1).expect("exact line knots"),
            &exact_points,
            &exact_weights,
        )
        .expect("exact weighted line");
        with_curve_cx(false, |cx| {
            assert_eq!(
                NurbsCurve::<Rat, 1>::new_with_cx(
                    KnotVector::new(vec![zero, zero, one, one], 1).expect("exact line knots"),
                    &exact_points,
                    &exact_weights,
                    cx,
                )
                .expect("active exact Cartesian construction"),
                CurveConstructionRun::Complete {
                    curve: exact.try_clone().expect("expected exact curve copy"),
                }
            );
        });

        with_curve_cx(true, |cx| {
            assert!(matches!(
                NurbsCurve::<f64, 4>::new_with_cx(
                    KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("line knots"),
                    &[[0.0; 4]; 2],
                    &[1.0; 2],
                    cx,
                ),
                Err(NurbsError::Structure { .. })
            ));
            assert!(matches!(
                NurbsCurve::<f64, 1>::new_with_cx(
                    KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("line knots"),
                    &[[0.0]],
                    &[1.0],
                    cx,
                ),
                Err(NurbsError::Structure { .. })
            ));
        });

        let control_bytes = core::mem::size_of::<[Rat; 4]>() as u128;
        let exact_control_count =
            usize::try_from(CURVE_CONSTRUCTION_MAX_RETAINED_BYTES / control_bytes)
                .expect("exact retained count fits usize");
        preflight_cartesian_curve_construction::<Rat>(exact_control_count)
            .expect("exact retained ceiling");
        assert!(matches!(
            preflight_cartesian_curve_construction::<Rat>(exact_control_count + 1),
            Err(NurbsError::Domain { .. })
        ));
    }

    #[test]
    fn cartesian_curve_construction_preserves_numeric_error_order() {
        let line_knots = || KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("line knots");
        let invalid_weight_points = [[0.0], [1.0]];
        let invalid_weights = [1.0, 0.0];
        let invalid_weight_error =
            NurbsCurve::<f64, 1>::new(line_knots(), &invalid_weight_points, &invalid_weights)
                .expect_err("legacy weight refusal");
        with_curve_cx(false, |cx| {
            assert_eq!(
                NurbsCurve::<f64, 1>::new_with_cx(
                    line_knots(),
                    &invalid_weight_points,
                    &invalid_weights,
                    cx,
                )
                .expect_err("cancellable weight refusal"),
                invalid_weight_error
            );
        });

        let nonfinite_points = [[0.0], [f64::INFINITY]];
        let finite_weights = [1.0, 1.0];
        let nonfinite_error =
            NurbsCurve::<f64, 1>::new(line_knots(), &nonfinite_points, &finite_weights)
                .expect_err("legacy coordinate refusal");
        with_curve_cx(false, |cx| {
            assert_eq!(
                NurbsCurve::<f64, 1>::new_with_cx(
                    line_knots(),
                    &nonfinite_points,
                    &finite_weights,
                    cx,
                )
                .expect_err("cancellable coordinate refusal"),
                nonfinite_error
            );
        });

        let underflow_points = [[f64::MIN_POSITIVE], [1.0]];
        let underflow_weights = [f64::MIN_POSITIVE, 1.0];
        let underflow_error =
            NurbsCurve::<f64, 1>::new(line_knots(), &underflow_points, &underflow_weights)
                .expect_err("legacy underflow refusal");
        with_curve_cx(false, |cx| {
            assert_eq!(
                NurbsCurve::<f64, 1>::new_with_cx(
                    line_knots(),
                    &underflow_points,
                    &underflow_weights,
                    cx,
                )
                .expect_err("cancellable underflow refusal"),
                underflow_error
            );
        });

        let overflow_points = [[f64::MAX], [1.0]];
        let overflow_weights = [2.0, 1.0];
        let overflow_error =
            NurbsCurve::<f64, 1>::new(line_knots(), &overflow_points, &overflow_weights)
                .expect_err("legacy overflow refusal");
        with_curve_cx(false, |cx| {
            assert_eq!(
                NurbsCurve::<f64, 1>::new_with_cx(
                    line_knots(),
                    &overflow_points,
                    &overflow_weights,
                    cx,
                )
                .expect_err("cancellable overflow refusal"),
                overflow_error
            );
        });
    }

    #[test]
    fn cartesian_curve_construction_cancels_in_scans_assembly_and_publication() {
        let inputs = || {
            let interior_count = 128usize;
            let mut knots = Vec::with_capacity(interior_count + 4);
            knots.extend([0.0, 0.0]);
            for index in 1..=interior_count {
                knots.push(index as f64 / (interior_count + 1) as f64);
            }
            knots.extend([1.0, 1.0]);
            let control_count = interior_count + 2;
            let points: Vec<[f64; 1]> = (0..control_count)
                .map(|index| [index as f64 / (control_count - 1) as f64])
                .collect();
            let weights = vec![1.0; control_count];
            (
                KnotVector::new(knots, 1).expect("long line knots"),
                points,
                weights,
            )
        };
        let run = |target| {
            let (knots, points, weights) = inputs();
            let mut polls = 0usize;
            let mut should_cancel = || {
                polls += 1;
                polls == target
            };
            let outcome =
                NurbsCurve::<f64, 1>::new_with_poll(knots, &points, &weights, &mut should_cancel)
                    .expect("valid cancellable construction");
            (matches!(outcome, CurveWorkRun::Cancelled), polls)
        };
        assert_eq!(run(12), run(12));
        assert_eq!(run(12), (true, 12));
        assert_eq!(run(20), run(20));
        assert_eq!(run(20), (true, 20));

        let (knots, points, weights) = inputs();
        let zero_dimensional_points = vec![[]; points.len()];
        let mut zero_dimensional_polls = 0usize;
        let mut cancel_inside_empty_points = || {
            zero_dimensional_polls += 1;
            zero_dimensional_polls == 15
        };
        assert!(matches!(
            NurbsCurve::<f64, 0>::new_with_poll(
                knots,
                &zero_dimensional_points,
                &weights,
                &mut cancel_inside_empty_points,
            )
            .expect("valid zero-dimensional cancellation"),
            CurveWorkRun::Cancelled
        ));
        assert_eq!(zero_dimensional_polls, 15);

        let (knots, points, weights) = inputs();
        let mut total_polls = 0usize;
        let mut never_cancel = || {
            total_polls += 1;
            false
        };
        assert!(matches!(
            NurbsCurve::<f64, 1>::new_with_poll(knots, &points, &weights, &mut never_cancel,)
                .expect("healthy Cartesian construction"),
            CurveWorkRun::Complete(_)
        ));
        let (knots, points, weights) = inputs();
        let mut replay_polls = 0usize;
        let mut cancel_at_publication = || {
            replay_polls += 1;
            replay_polls == total_polls
        };
        assert!(matches!(
            NurbsCurve::<f64, 1>::new_with_poll(
                knots,
                &points,
                &weights,
                &mut cancel_at_publication,
            )
            .expect("publication cancellation"),
            CurveWorkRun::Cancelled
        ));
        assert_eq!(replay_polls, total_polls);
    }

    #[test]
    fn homogeneous_curve_construction_with_cx_is_transactional_and_exact() {
        let expected = line_curve();
        with_curve_cx(true, |cx| {
            assert_eq!(
                NurbsCurve::<f64, 1>::from_homogeneous_with_cx(
                    expected.knots.clone(),
                    expected.cpw.clone(),
                    cx,
                )
                .expect("valid pre-cancelled construction"),
                CurveConstructionRun::Cancelled
            );
        });
        with_curve_cx(false, |cx| {
            assert_eq!(
                NurbsCurve::<f64, 1>::from_homogeneous_with_cx(
                    expected.knots.clone(),
                    expected.cpw.clone(),
                    cx,
                )
                .expect("active homogeneous construction"),
                CurveConstructionRun::Complete {
                    curve: expected.try_clone().expect("expected curve copy"),
                }
            );
        });

        let zero = Rat::int(0);
        let one = Rat::int(1);
        let exact = NurbsCurve::<Rat, 1>::from_homogeneous(
            KnotVector::new(vec![zero, zero, one, one], 1).expect("exact line knots"),
            vec![[zero, zero, zero, one], [one, zero, zero, one]],
        )
        .expect("exact homogeneous line");
        with_curve_cx(false, |cx| {
            assert_eq!(
                NurbsCurve::<Rat, 1>::from_homogeneous_with_cx(
                    exact.knots.clone(),
                    exact.cpw.clone(),
                    cx,
                )
                .expect("active exact construction"),
                CurveConstructionRun::Complete {
                    curve: exact.try_clone().expect("expected exact curve copy"),
                }
            );
        });

        with_curve_cx(true, |cx| {
            assert!(matches!(
                NurbsCurve::<f64, 4>::from_homogeneous_with_cx(
                    KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("line knots"),
                    vec![[0.0, 0.0, 0.0, 1.0]; 2],
                    cx,
                ),
                Err(NurbsError::Structure { .. })
            ));
        });

        let over_cap_control_count = 1_048_572usize;
        let line_knots = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("line knots");
        assert_eq!(
            NurbsCurve::<f64, 1>::validation_work_for(&line_knots, over_cap_control_count)
                .expect("cap-plus-one work"),
            BASIS_MAX_WORK_UNITS + 1
        );
        with_curve_cx(true, |cx| {
            assert!(matches!(
                NurbsCurve::<f64, 1>::from_homogeneous_with_cx(
                    line_knots,
                    vec![[0.0, 0.0, 0.0, 1.0]; over_cap_control_count],
                    cx,
                ),
                Err(NurbsError::Domain { .. })
            ));
        });

        let legacy_error = NurbsCurve::<f64, 1>::from_homogeneous(
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("line knots"),
            vec![[0.0, 0.0, 0.0, 1.0]],
        )
        .expect_err("legacy control-count mismatch");
        with_curve_cx(false, |cx| {
            assert_eq!(
                NurbsCurve::<f64, 1>::from_homogeneous_with_cx(
                    KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("line knots"),
                    vec![[0.0, 0.0, 0.0, 1.0]],
                    cx,
                )
                .expect_err("cancellable control-count mismatch"),
                legacy_error
            );
        });
    }

    #[test]
    fn homogeneous_curve_construction_cancels_inside_controls_and_at_publication() {
        let source = long_linear_curve();
        let run = |target| {
            let mut polls = 0usize;
            let mut should_cancel = || {
                polls += 1;
                polls == target
            };
            let outcome = NurbsCurve::<f64, 1>::from_homogeneous_with_poll(
                source.knots.clone(),
                source.cpw.clone(),
                &mut should_cancel,
            )
            .expect("valid cancellable construction");
            (matches!(outcome, CurveWorkRun::Cancelled), polls)
        };
        assert_eq!(run(13), run(13));
        assert_eq!(run(13), (true, 13));

        let source = line_curve();
        let mut total_polls = 0usize;
        let mut never_cancel = || {
            total_polls += 1;
            false
        };
        assert!(matches!(
            NurbsCurve::<f64, 1>::from_homogeneous_with_poll(
                source.knots.clone(),
                source.cpw.clone(),
                &mut never_cancel,
            )
            .expect("healthy construction"),
            CurveWorkRun::Complete(_)
        ));
        let mut replay_polls = 0usize;
        let mut cancel_at_publication = || {
            replay_polls += 1;
            replay_polls == total_polls
        };
        assert!(matches!(
            NurbsCurve::<f64, 1>::from_homogeneous_with_poll(
                source.knots.clone(),
                source.cpw.clone(),
                &mut cancel_at_publication,
            )
            .expect("publication cancellation"),
            CurveWorkRun::Cancelled
        ));
        assert_eq!(replay_polls, total_polls);
    }

    #[test]
    fn admitted_curve_evaluation_with_cx_is_transactional_and_exact() {
        let curve = line_curve();
        let admitted = curve.admit().expect("admitted line");
        with_curve_cx(true, |cx| {
            assert!(matches!(
                curve.eval_with_cx(0.25, cx).expect("valid owning request"),
                CurveEvaluationRun::Cancelled
            ));
            assert!(matches!(
                admitted.eval_with_cx(0.25, cx).expect("valid request"),
                CurveEvaluationRun::Cancelled
            ));
            assert!(matches!(
                admitted.eval_with_cx(-1.0, cx),
                Err(NurbsError::Domain { .. })
            ));
        });
        with_curve_cx(false, |cx| {
            assert_eq!(
                curve.eval_with_cx(0.25, cx).expect("active owning context"),
                CurveEvaluationRun::Complete {
                    point: curve.eval(0.25).expect("legacy owning evaluation"),
                }
            );
            assert_eq!(
                admitted.eval_with_cx(0.25, cx).expect("active context"),
                CurveEvaluationRun::Complete {
                    point: admitted.eval(0.25).expect("legacy evaluation"),
                }
            );
        });

        let mut admission_polls = 0usize;
        let mut observe_admission = || {
            admission_polls += 1;
            false
        };
        assert!(matches!(
            curve
                .validate_live_structure_with_poll(&mut observe_admission)
                .expect("healthy source admission"),
            CurveWorkRun::Complete(())
        ));
        let mut owning_polls = 0usize;
        let mut cancel_at_first_evaluation_poll = || {
            owning_polls += 1;
            owning_polls == admission_polls + 1
        };
        assert_eq!(
            curve
                .eval_with_poll(0.25, &mut cancel_at_first_evaluation_poll)
                .expect("owning evaluation cancellation"),
            CurveEvaluationRun::Cancelled
        );
        assert_eq!(owning_polls, admission_polls + 1);

        let invalid_dimension = NurbsCurve::<f64, 4> {
            knots: KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("line knots"),
            cpw: vec![[0.0, 0.0, 0.0, 1.0]; 2],
        };
        with_curve_cx(true, |cx| {
            assert!(matches!(
                invalid_dimension.eval_with_cx(0.25, cx),
                Err(NurbsError::Structure { .. })
            ));
            assert!(matches!(
                invalid_dimension.eval_homogeneous_with_cx(0.25, cx),
                Err(NurbsError::Structure { .. })
            ));
        });
    }

    #[test]
    fn admitted_homogeneous_curve_evaluation_with_cx_is_transactional_and_exact() {
        let curve = line_curve();
        let admitted = curve.admit().expect("admitted line");
        with_curve_cx(true, |cx| {
            assert_eq!(
                curve
                    .eval_homogeneous_with_cx(0.25, cx)
                    .expect("valid owning homogeneous request"),
                CurveHomogeneousEvaluationRun::Cancelled
            );
            assert_eq!(
                admitted
                    .eval_homogeneous_with_cx(0.25, cx)
                    .expect("valid homogeneous request"),
                CurveHomogeneousEvaluationRun::Cancelled
            );
            assert!(matches!(
                admitted.eval_homogeneous_with_cx(-1.0, cx),
                Err(NurbsError::Domain { .. })
            ));
        });
        with_curve_cx(false, |cx| {
            assert_eq!(
                curve
                    .eval_homogeneous_with_cx(0.25, cx)
                    .expect("active owning homogeneous context"),
                CurveHomogeneousEvaluationRun::Complete {
                    homogeneous: curve
                        .eval_homogeneous(0.25)
                        .expect("legacy owning homogeneous evaluation"),
                }
            );
            assert_eq!(
                admitted
                    .eval_homogeneous_with_cx(0.25, cx)
                    .expect("active homogeneous context"),
                CurveHomogeneousEvaluationRun::Complete {
                    homogeneous: admitted
                        .eval_homogeneous(0.25)
                        .expect("legacy homogeneous evaluation"),
                }
            );
        });

        let mut admission_polls = 0usize;
        let mut observe_admission = || {
            admission_polls += 1;
            false
        };
        assert!(matches!(
            curve
                .validate_live_structure_with_poll(&mut observe_admission)
                .expect("healthy source admission"),
            CurveWorkRun::Complete(())
        ));
        let mut owning_polls = 0usize;
        let mut cancel_at_first_evaluation_poll = || {
            owning_polls += 1;
            owning_polls == admission_polls + 1
        };
        assert_eq!(
            curve
                .eval_homogeneous_with_poll(0.25, &mut cancel_at_first_evaluation_poll)
                .expect("owning homogeneous evaluation cancellation"),
            CurveHomogeneousEvaluationRun::Cancelled
        );
        assert_eq!(owning_polls, admission_polls + 1);

        let exact = NurbsCurve::<Rat, 1>::new(
            KnotVector::new(vec![Rat::int(0), Rat::int(0), Rat::int(1), Rat::int(1)], 1)
                .expect("exact line knots"),
            &[[Rat::int(0)], [Rat::int(2)]],
            &[Rat::int(1), Rat::int(2)],
        )
        .expect("exact rational line");
        let exact_admitted = exact.admit().expect("admitted exact line");
        let parameter = Rat::new(1, 2);
        with_curve_cx(false, |cx| {
            assert_eq!(
                exact_admitted
                    .eval_homogeneous_with_cx(parameter, cx)
                    .expect("active exact homogeneous context"),
                CurveHomogeneousEvaluationRun::Complete {
                    homogeneous: exact_admitted
                        .eval_homogeneous(parameter)
                        .expect("legacy exact homogeneous evaluation"),
                }
            );
        });
    }

    #[test]
    fn curve_copy_with_cx_is_transactional_and_exact() {
        let curve = line_curve();
        with_curve_cx(true, |cx| {
            assert_eq!(
                curve.try_clone_with_cx(cx).expect("admitted copy request"),
                CurveCloneRun::Cancelled
            );
        });
        with_curve_cx(false, |cx| {
            assert_eq!(
                curve.try_clone_with_cx(cx).expect("active curve copy"),
                CurveCloneRun::Complete {
                    curve: curve.try_clone().expect("legacy curve copy"),
                }
            );
        });

        let exact = NurbsCurve::<Rat, 1>::new(
            KnotVector::new(vec![Rat::int(0), Rat::int(0), Rat::int(1), Rat::int(1)], 1)
                .expect("exact line knots"),
            &[[Rat::int(0)], [Rat::int(1)]],
            &[Rat::int(1), Rat::int(1)],
        )
        .expect("exact line");
        with_curve_cx(false, |cx| {
            assert_eq!(
                exact.try_clone_with_cx(cx).expect("active exact copy"),
                CurveCloneRun::Complete {
                    curve: exact.try_clone().expect("legacy exact copy"),
                }
            );
        });
    }

    #[test]
    fn curve_copy_cancels_inside_each_linear_copy_and_at_publication() {
        let curve = long_linear_curve();
        let run = |target| {
            let mut polls = 0usize;
            let mut should_cancel = || {
                polls += 1;
                polls == target
            };
            let outcome = curve
                .try_clone_with_poll(&mut should_cancel)
                .expect("bounded curve copy");
            (matches!(outcome, CurveWorkRun::Cancelled), polls)
        };
        assert_eq!(run(2), run(2));
        assert_eq!(run(2), (true, 2));
        assert_eq!(run(5), run(5));
        assert_eq!(run(5), (true, 5));

        let mut total_polls = 0usize;
        let mut never_cancel = || {
            total_polls += 1;
            false
        };
        assert!(matches!(
            curve
                .try_clone_with_poll(&mut never_cancel)
                .expect("healthy curve copy"),
            CurveWorkRun::Complete(_)
        ));
        let mut replay_polls = 0usize;
        let mut cancel_at_publication = || {
            replay_polls += 1;
            replay_polls == total_polls
        };
        assert!(matches!(
            curve
                .try_clone_with_poll(&mut cancel_at_publication)
                .expect("publication cancellation"),
            CurveWorkRun::Cancelled
        ));
        assert_eq!(replay_polls, total_polls);
    }

    #[test]
    fn curve_copy_preflight_refuses_work_before_retained_bytes() {
        let oversized = usize::MAX;
        let work_error = preflight_curve_copy::<f64>(oversized, oversized)
            .expect_err("work must refuse before retained-byte arithmetic");
        assert!(matches!(
            work_error,
            NurbsError::Domain { ref what } if what.contains("work units above defensive ceiling")
        ));

        let memory_error =
            preflight_curve_copy::<f64>(0, 2_100_000).expect_err("copy output must exceed 64 MiB");
        assert!(matches!(
            memory_error,
            NurbsError::Domain { ref what } if what.contains("retains")
        ));
    }

    #[test]
    fn admitted_knot_insertion_with_cx_is_transactional_and_exact() {
        let curve = line_curve();
        let admitted = curve.admit().expect("admitted line");
        let legacy = admitted.insert_knot(0.5).expect("legacy insertion");
        let endpoint_error = admitted
            .insert_knot(0.0)
            .expect_err("legacy endpoint refusal");
        let non_finite_error = admitted
            .insert_knot(f64::NAN)
            .expect_err("legacy non-finite refusal");

        with_curve_cx(true, |cx| {
            assert_eq!(
                admitted
                    .insert_knot_with_cx(0.5, cx)
                    .expect("valid request reaches cancellation"),
                CurveInsertionRun::Cancelled
            );
            assert_eq!(
                admitted
                    .insert_knot_with_cx(0.0, cx)
                    .expect_err("endpoint refusal must beat cancellation"),
                endpoint_error
            );
            assert_eq!(
                admitted
                    .insert_knot_with_cx(f64::NAN, cx)
                    .expect_err("non-finite refusal must beat cancellation"),
                non_finite_error
            );
        });
        with_curve_cx(false, |cx| {
            assert_eq!(
                admitted
                    .insert_knot_with_cx(0.5, cx)
                    .expect("active insertion"),
                CurveInsertionRun::Complete { curve: legacy }
            );
        });

        let exact_curve = NurbsCurve::<Rat, 1>::new(
            KnotVector::new(vec![Rat::int(0), Rat::int(0), Rat::int(1), Rat::int(1)], 1)
                .expect("exact line knots"),
            &[[Rat::int(0)], [Rat::int(1)]],
            &[Rat::int(1), Rat::int(1)],
        )
        .expect("exact line");
        let exact_admitted = exact_curve.admit().expect("admitted exact line");
        let half = Rat::new(1, 2);
        let exact_legacy = exact_admitted
            .insert_knot(half)
            .expect("legacy exact insertion");
        with_curve_cx(false, |cx| {
            assert_eq!(
                exact_admitted
                    .insert_knot_with_cx(half, cx)
                    .expect("active exact insertion"),
                CurveInsertionRun::Complete {
                    curve: exact_legacy,
                }
            );
        });
    }

    #[test]
    fn admitted_knot_removal_with_cx_is_transactional_and_exact() {
        let original = line_curve();
        let inserted = original.insert_knot(0.5).expect("insert midpoint");
        let admitted = inserted.admit().expect("admitted refined line");
        let legacy = admitted.remove_knot(0.5).expect("legacy removal");
        assert_eq!(legacy, original);

        let endpoint_error = admitted
            .remove_knot(0.0)
            .expect_err("legacy endpoint refusal");
        with_curve_cx(true, |cx| {
            assert_eq!(
                admitted
                    .remove_knot_with_cx(0.5, cx)
                    .expect("valid request reaches cancellation"),
                CurveRemovalRun::Cancelled
            );
            assert_eq!(
                admitted
                    .remove_knot_with_cx(0.0, cx)
                    .expect_err("endpoint refusal must beat cancellation"),
                endpoint_error
            );
            assert_eq!(
                admitted
                    .remove_knot_with_cx(0.25, cx)
                    .expect("cancellation wins before absent-knot scan refusal"),
                CurveRemovalRun::Cancelled
            );
        });
        with_curve_cx(false, |cx| {
            assert_eq!(
                admitted
                    .remove_knot_with_cx(0.5, cx)
                    .expect("active removal"),
                CurveRemovalRun::Complete {
                    curve: line_curve(),
                }
            );
            assert!(matches!(
                admitted.remove_knot_with_cx(0.25, cx),
                Err(NurbsError::Domain { .. })
            ));
        });

        let zero = Rat::int(0);
        let one = Rat::int(1);
        let half = Rat::new(1, 2);
        let exact_original = NurbsCurve::<Rat, 1>::new(
            KnotVector::new(vec![zero, zero, one, one], 1).expect("exact line knots"),
            &[[zero], [one]],
            &[one, one],
        )
        .expect("exact line");
        let exact_inserted = exact_original
            .insert_knot(half)
            .expect("exact midpoint insertion");
        let exact_admitted = exact_inserted.admit().expect("admitted exact refinement");
        with_curve_cx(false, |cx| {
            assert_eq!(
                exact_admitted
                    .remove_knot_with_cx(half, cx)
                    .expect("active exact removal"),
                CurveRemovalRun::Complete {
                    curve: exact_original,
                }
            );
        });
        let repeated = exact_inserted
            .insert_knot(half)
            .expect("raise midpoint multiplicity");
        assert_eq!(
            repeated
                .admit()
                .expect("admitted repeated knot")
                .remove_knot(half)
                .expect("remove inserted repeated knot"),
            exact_inserted
        );

        let discontinuous = NurbsCurve::<Rat, 1>::new(
            KnotVector::new(vec![zero, zero, half, half, one, one], 1).expect("full-break knots"),
            &[[zero], [one], [Rat::int(3)], [Rat::int(4)]],
            &[one; 4],
        )
        .expect("discontinuous curve");
        assert!(matches!(
            discontinuous
                .admit()
                .expect("admitted discontinuity")
                .remove_knot(half),
            Err(NurbsError::Structure { .. })
        ));
    }

    #[test]
    fn knot_removal_cancels_inside_scan_and_at_publication() {
        let shape = long_linear_curve();
        let points = vec![[0.0]; shape.cpw.len()];
        let weights = vec![1.0; shape.cpw.len()];
        let source = NurbsCurve::new(
            shape.knots.try_clone().expect("long knot copy"),
            &points,
            &weights,
        )
        .expect("long constant curve");
        let inserted = source.insert_knot(0.5).expect("insert absent midpoint");
        let admitted = inserted.admit().expect("admitted long refinement");

        let scan_run = || {
            let mut polls = 0usize;
            let mut should_cancel = || {
                polls += 1;
                polls == 2
            };
            let outcome = admitted
                .remove_knot_with_poll(0.5, &mut should_cancel)
                .expect("valid removal scan");
            (matches!(outcome, CurveRemovalRun::Cancelled), polls)
        };
        assert_eq!(scan_run(), scan_run());
        assert_eq!(scan_run(), (true, 2));

        let mut total_polls = 0usize;
        let mut never_cancel = || {
            total_polls += 1;
            false
        };
        assert!(matches!(
            admitted
                .remove_knot_with_poll(0.5, &mut never_cancel)
                .expect("healthy removal"),
            CurveRemovalRun::Complete { .. }
        ));
        let mut comparison_polls = 0usize;
        let comparison_target = total_polls
            .checked_sub(1)
            .expect("healthy removal has a pre-publication comparison checkpoint");
        let mut cancel_after_reinsertion = || {
            comparison_polls += 1;
            comparison_polls == comparison_target
        };
        assert_eq!(
            admitted
                .remove_knot_with_poll(0.5, &mut cancel_after_reinsertion)
                .expect("post-reinsertion comparison cancellation"),
            CurveRemovalRun::Cancelled
        );
        assert_eq!(comparison_polls, comparison_target);

        let mut replay_polls = 0usize;
        let mut cancel_at_publication = || {
            replay_polls += 1;
            replay_polls == total_polls
        };
        assert_eq!(
            admitted
                .remove_knot_with_poll(0.5, &mut cancel_at_publication)
                .expect("publication cancellation"),
            CurveRemovalRun::Cancelled
        );
        assert_eq!(replay_polls, total_polls);
    }

    #[test]
    fn knot_removal_envelope_refuses_work_before_combined_derived_storage() {
        let oversized = usize::MAX;
        let work_error = enforce_curve_removal_envelope::<f64>(
            oversized,
            oversized,
            oversized - 1,
            oversized - 1,
            BASIS_MAX_WORK_UNITS + 1,
        )
        .expect_err("work cap must refuse before derived-byte arithmetic");
        assert!(matches!(
            work_error,
            NurbsError::Domain { ref what } if what.contains("work units above defensive ceiling")
        ));

        let memory_error = enforce_curve_removal_envelope::<f64>(
            oversized,
            oversized,
            oversized - 1,
            oversized - 1,
            0,
        )
        .expect_err("combined derived storage must exceed 64 MiB");
        assert!(matches!(
            memory_error,
            NurbsError::Domain { ref what } if what.contains("simultaneously-live derived bytes")
        ));
    }

    #[test]
    fn knot_removal_cancels_inside_exact_reconstruction() {
        let degree = 16usize;
        let mut knots = vec![0.0; degree + 1];
        knots.extend(vec![1.0; degree + 1]);
        let points = vec![[0.0]; degree + 1];
        let weights = vec![1.0; degree + 1];
        let source = NurbsCurve::new(
            KnotVector::new(knots, degree).expect("high-degree knots"),
            &points,
            &weights,
        )
        .expect("high-degree constant curve");
        let inserted = source.insert_knot(0.5).expect("high-degree insertion");
        let plan = inserted
            .removal_plan_after_parameter(0.5)
            .expect("removal plan");
        let mut never_cancel = || false;
        let CurveWorkRun::Complete((removed_index, prior_multiplicity)) = inserted
            .find_removal_site_with_poll(0.5, &mut never_cancel)
            .expect("removal site")
        else {
            panic!("healthy site scan must complete");
        };
        let CurveWorkRun::Complete(new_knots) = inserted
            .copy_removed_knots_with_poll(removed_index, plan, &mut never_cancel)
            .expect("removed knot copy")
        else {
            panic!("healthy knot copy must complete");
        };

        let run = || {
            let mut polls = 0usize;
            let mut should_cancel = || {
                polls += 1;
                polls == 2
            };
            let outcome = inserted
                .reconstruct_removed_controls_with_poll(
                    0.5,
                    removed_index,
                    prior_multiplicity,
                    &new_knots,
                    plan,
                    &mut should_cancel,
                )
                .expect("valid exact reconstruction");
            (matches!(outcome, CurveWorkRun::Cancelled), polls)
        };
        assert_eq!(run(), run());
        assert_eq!(run(), (true, 2));
    }

    // G0/G4: the owning wrapper preserves refusal order and exact healthy output.
    #[test]
    fn owning_curve_derivatives_with_cx_are_transactional_and_exact() {
        let curve = line_curve();
        with_curve_cx(true, |cx| {
            assert!(matches!(
                curve
                    .derivatives_with_cx(0.25, 1, cx)
                    .expect("valid pre-cancelled derivative request"),
                CurveDerivativesRun::Cancelled
            ));

            let parameter_error = curve
                .derivatives(-1.0, 1)
                .expect_err("legacy parameter refusal");
            assert_eq!(
                curve
                    .derivatives_with_cx(-1.0, 1, cx)
                    .expect_err("parameter refusal must precede cancellation"),
                parameter_error
            );
            let order_error = curve
                .derivatives(0.25, 65)
                .expect_err("legacy order refusal");
            assert_eq!(
                curve
                    .derivatives_with_cx(0.25, 65, cx)
                    .expect_err("order refusal must precede cancellation"),
                order_error
            );

            let invalid_dimension = NurbsCurve::<f64, 4> {
                knots: KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("line knots"),
                cpw: vec![[0.0, 0.0, 0.0, 1.0]; 2],
            };
            let dimension_error = invalid_dimension
                .derivatives(0.25, 1)
                .expect_err("legacy dimension refusal");
            assert_eq!(
                invalid_dimension
                    .derivatives_with_cx(0.25, 1, cx)
                    .expect_err("dimension refusal must precede cancellation"),
                dimension_error
            );
        });
        with_curve_cx(false, |cx| {
            assert_eq!(
                curve
                    .derivatives_with_cx(0.25, 1, cx)
                    .expect("active owning derivative request"),
                CurveDerivativesRun::Complete {
                    derivatives: curve.derivatives(0.25, 1).expect("legacy derivatives"),
                }
            );
        });
    }

    // G4/G5: one callback spans owning admission and admitted derivative work.
    #[test]
    fn owning_curve_derivative_cancellation_replays_across_both_phases() {
        let long_curve = long_linear_curve();
        let run_admission = || {
            let mut polls = 0usize;
            let mut should_cancel = || {
                polls += 1;
                polls == 13
            };
            let outcome = long_curve
                .derivatives_with_poll(0.25, 1, &mut should_cancel)
                .expect("valid long-curve derivative request");
            (matches!(outcome, CurveDerivativesRun::Cancelled), polls)
        };
        assert_eq!(run_admission(), run_admission());
        assert_eq!(run_admission(), (true, 13));

        let curve = line_curve();
        let run_derivative = || {
            let mut polls = 0usize;
            let mut should_cancel = || {
                polls += 1;
                polls == 8
            };
            let outcome = curve
                .derivatives_with_poll(0.25, 1, &mut should_cancel)
                .expect("valid line derivative request");
            (matches!(outcome, CurveDerivativesRun::Cancelled), polls)
        };
        assert_eq!(run_derivative(), run_derivative());
        assert_eq!(run_derivative(), (true, 8));
    }

    #[test]
    fn admitted_curve_derivatives_with_cx_are_transactional_and_exact() {
        let curve = line_curve();
        let admitted = curve.admit().expect("admitted line");
        with_curve_cx(true, |cx| {
            assert!(matches!(
                admitted
                    .derivatives_with_cx(0.25, 1, cx)
                    .expect("valid derivative request"),
                CurveDerivativesRun::Cancelled
            ));

            let parameter_error = admitted
                .derivatives(-1.0, 1)
                .expect_err("legacy parameter refusal");
            assert_eq!(
                admitted
                    .derivatives_with_cx(-1.0, 1, cx)
                    .expect_err("parameter refusal must precede cancellation"),
                parameter_error
            );
            let order_error = admitted
                .derivatives(0.25, 65)
                .expect_err("legacy order refusal");
            assert_eq!(
                admitted
                    .derivatives_with_cx(0.25, 65, cx)
                    .expect_err("order refusal must precede cancellation"),
                order_error
            );
        });
        with_curve_cx(false, |cx| {
            for order in [0, 1] {
                assert_eq!(
                    admitted
                        .derivatives_with_cx(0.25, order, cx)
                        .expect("active derivative request"),
                    CurveDerivativesRun::Complete {
                        derivatives: admitted
                            .derivatives(0.25, order)
                            .expect("legacy derivative request"),
                    }
                );
            }
        });

        let weighted = NurbsCurve::new(
            KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("weighted line knots"),
            &[[0.0], [1.0]],
            &[1.0, 2.0],
        )
        .expect("weighted rational line");
        let admitted_weighted = weighted.admit().expect("admitted weighted line");
        with_curve_cx(false, |cx| {
            assert_eq!(
                admitted_weighted
                    .derivatives_with_cx(0.25, 2, cx)
                    .expect("active rational derivative request"),
                CurveDerivativesRun::Complete {
                    derivatives: admitted_weighted
                        .derivatives(0.25, 2)
                        .expect("legacy rational derivatives"),
                }
            );
        });

        let c0_knots = KnotVector::new(vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0], 2)
            .expect("C0 quadratic knots");
        let c0_curve = NurbsCurve::new(c0_knots, &[[0.0], [0.25], [0.5], [0.75], [1.0]], &[1.0; 5])
            .expect("C0 curve");
        let admitted_c0 = c0_curve.admit().expect("admitted C0 curve");
        let continuity_error = admitted_c0
            .derivatives(0.5, 1)
            .expect_err("legacy ordinary derivative refusal");
        with_curve_cx(true, |cx| {
            assert_eq!(
                admitted_c0
                    .derivatives_with_cx(0.5, 1, cx)
                    .expect_err("continuity refusal must precede cancellation"),
                continuity_error
            );
        });
    }

    #[test]
    fn curve_derivative_cancellation_replays_inside_work_and_at_publication() {
        let high_degree = high_degree_insertion_curve();
        let admitted_high_degree = high_degree.admit().expect("admitted high-degree curve");
        NurbsCurve::<f64, 1>::preflight_derivative_request(admitted_high_degree.knots(), 0.25, 1)
            .expect("valid high-degree derivative request");
        let run = || {
            let mut polls = 0usize;
            let mut should_cancel = || {
                polls += 1;
                polls == 4
            };
            let outcome =
                NurbsCurve::<f64, 1>::derivatives_from_admitted_parts_after_preflight_with_poll(
                    admitted_high_degree.knots(),
                    admitted_high_degree.homogeneous_control_points(),
                    0.25,
                    1,
                    &mut should_cancel,
                )
                .expect("valid high-degree derivative work");
            (matches!(outcome, CurveDerivativesRun::Cancelled), polls)
        };
        assert_eq!(run(), run());
        assert_eq!(run(), (true, 4));

        let curve = line_curve();
        let admitted = curve.admit().expect("admitted line");
        let mut total_polls = 0usize;
        let mut never_cancel = || {
            total_polls += 1;
            false
        };
        assert!(matches!(
            NurbsCurve::<f64, 1>::derivatives_from_admitted_parts_after_preflight_with_poll(
                admitted.knots(),
                admitted.homogeneous_control_points(),
                0.25,
                1,
                &mut never_cancel,
            )
            .expect("healthy derivative work"),
            CurveDerivativesRun::Complete { .. }
        ));
        assert_eq!(total_polls, 14);
        let mut replay_polls = 0usize;
        let mut cancel_at_publication = || {
            replay_polls += 1;
            replay_polls == 14
        };
        assert!(matches!(
            NurbsCurve::<f64, 1>::derivatives_from_admitted_parts_after_preflight_with_poll(
                admitted.knots(),
                admitted.homogeneous_control_points(),
                0.25,
                1,
                &mut cancel_at_publication,
            )
            .expect("publication cancellation"),
            CurveDerivativesRun::Cancelled
        ));
        assert_eq!(replay_polls, 14);

        let mut pre_cancelled = || true;
        assert!(matches!(
            NurbsCurve::<f64, 1>::derivatives_from_admitted_parts_after_preflight_with_poll(
                admitted.knots(),
                &[],
                0.25,
                1,
                &mut pre_cancelled,
            ),
            Err(NurbsError::Structure { .. })
        ));
    }

    #[test]
    fn admitted_curve_span_boxes_with_cx_are_transactional_and_exact() {
        let curve = line_curve();
        let admitted = curve.admit().expect("admitted line");
        with_curve_cx(true, |cx| {
            assert!(matches!(
                curve
                    .span_boxes_with_cx(cx)
                    .expect("valid owning span-box request"),
                CurveSpanBoxesRun::Cancelled
            ));
            assert!(matches!(
                admitted
                    .span_boxes_with_cx(cx)
                    .expect("valid span-box request"),
                CurveSpanBoxesRun::Cancelled
            ));
        });
        with_curve_cx(false, |cx| {
            assert_eq!(
                curve
                    .span_boxes_with_cx(cx)
                    .expect("active owning span-box request"),
                CurveSpanBoxesRun::Complete {
                    boxes: curve.span_boxes().expect("legacy owning span boxes"),
                }
            );
            assert_eq!(
                admitted
                    .span_boxes_with_cx(cx)
                    .expect("active span-box request"),
                CurveSpanBoxesRun::Complete {
                    boxes: admitted.span_boxes().expect("legacy span boxes"),
                }
            );
        });

        let mut admission_polls = 0usize;
        let mut observe_admission = || {
            admission_polls += 1;
            false
        };
        assert!(matches!(
            curve
                .validate_live_structure_with_poll(&mut observe_admission)
                .expect("healthy source admission"),
            CurveWorkRun::Complete(())
        ));
        let mut owning_polls = 0usize;
        let mut cancel_at_first_box_poll = || {
            owning_polls += 1;
            owning_polls == admission_polls + 1
        };
        assert_eq!(
            curve
                .span_boxes_with_poll(&mut cancel_at_first_box_poll)
                .expect("owning span-box cancellation"),
            CurveSpanBoxesRun::Cancelled
        );
        assert_eq!(owning_polls, admission_polls + 1);

        let invalid_dimension = NurbsCurve::<f64, 4> {
            knots: KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("line knots"),
            cpw: vec![[0.0, 0.0, 0.0, 1.0]; 2],
        };
        with_curve_cx(true, |cx| {
            assert!(matches!(
                invalid_dimension.span_boxes_with_cx(cx),
                Err(NurbsError::Structure { .. })
            ));
        });

        let degree = 1_024usize;
        let mut high_overlap_knots = vec![0.0; degree + 1];
        high_overlap_knots.extend(vec![0.5; degree]);
        high_overlap_knots.extend(vec![1.0; degree + 1]);
        let control_count = 2 * degree + 1;
        let high_overlap = NurbsCurve::new(
            KnotVector::new(high_overlap_knots, degree).expect("high-overlap knots"),
            &vec![[0.0]; control_count],
            &vec![1.0; control_count],
        )
        .expect("high-overlap curve");
        let admitted_high_overlap = high_overlap.admit().expect("admitted high-overlap curve");
        let work_error = admitted_high_overlap
            .span_boxes()
            .expect_err("legacy span-box work refusal");
        with_curve_cx(true, |cx| {
            assert_eq!(
                admitted_high_overlap
                    .span_boxes_with_cx(cx)
                    .expect_err("work refusal must precede cancellation"),
                work_error
            );
        });

        let zero = Rat::int(0);
        let one = Rat::int(1);
        let exact = NurbsCurve::new(
            KnotVector::new(vec![zero, zero, one, one], 1).expect("exact line knots"),
            &[[zero], [one]],
            &[one, one],
        )
        .expect("exact line");
        let admitted_exact = exact.admit().expect("admitted exact line");
        with_curve_cx(false, |cx| {
            assert_eq!(
                exact
                    .span_boxes_with_cx(cx)
                    .expect("active exact owning span-box request"),
                CurveSpanBoxesRun::Complete {
                    boxes: exact.span_boxes().expect("legacy exact owning span boxes"),
                }
            );
            assert_eq!(
                admitted_exact
                    .span_boxes_with_cx(cx)
                    .expect("active exact span-box request"),
                CurveSpanBoxesRun::Complete {
                    boxes: admitted_exact
                        .span_boxes()
                        .expect("legacy exact span boxes"),
                }
            );
        });
    }

    #[test]
    fn curve_span_box_cancellation_replays_inside_work_and_at_publication() {
        let high_degree = high_degree_insertion_curve();
        let run = || {
            let mut polls = 0usize;
            let mut should_cancel = || {
                polls += 1;
                polls == 2
            };
            let outcome = high_degree
                .span_boxes_after_validation_with_poll(&mut should_cancel)
                .expect("valid high-degree span-box work");
            (matches!(outcome, CurveSpanBoxesRun::Cancelled), polls)
        };
        assert_eq!(run(), run());
        assert_eq!(run(), (true, 2));

        let long = long_linear_curve();
        let zero_points: Vec<[f64; 0]> = vec![[]; long.cpw.len()];
        let zero_weights = vec![1.0; long.cpw.len()];
        let zero_dim = NurbsCurve::<f64, 0>::new(
            long.knots.try_clone().expect("zero-D knot copy"),
            &zero_points,
            &zero_weights,
        )
        .expect("long zero-D curve");
        let mut zero_dim_polls = 0usize;
        let mut cancel_zero_dim = || {
            zero_dim_polls += 1;
            zero_dim_polls == 2
        };
        assert!(matches!(
            zero_dim
                .span_boxes_after_validation_with_poll(&mut cancel_zero_dim)
                .expect("zero-D span-box cancellation"),
            CurveSpanBoxesRun::Cancelled
        ));
        assert_eq!(zero_dim_polls, 2);

        let curve = line_curve();
        let mut total_polls = 0usize;
        let mut never_cancel = || {
            total_polls += 1;
            false
        };
        assert!(matches!(
            curve
                .span_boxes_after_validation_with_poll(&mut never_cancel)
                .expect("healthy span-box work"),
            CurveSpanBoxesRun::Complete { .. }
        ));
        assert_eq!(total_polls, 2);
        let mut replay_polls = 0usize;
        let mut cancel_at_publication = || {
            replay_polls += 1;
            replay_polls == 2
        };
        assert!(matches!(
            curve
                .span_boxes_after_validation_with_poll(&mut cancel_at_publication)
                .expect("publication cancellation"),
            CurveSpanBoxesRun::Cancelled
        ));
        assert_eq!(replay_polls, 2);
    }

    #[test]
    fn curve_admission_with_cx_is_transactional_and_lifetime_bound() {
        let curve = line_curve();
        with_curve_cx(true, |cx| {
            assert!(matches!(
                curve.admit_with_cx(cx).expect("valid source"),
                CurveAdmissionRun::Cancelled
            ));

            let invalid_dimension = NurbsCurve::<f64, 4> {
                knots: KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("line knots"),
                cpw: vec![[0.0, 0.0, 0.0, 1.0]; 2],
            };
            assert!(matches!(
                invalid_dimension.admit_with_cx(cx),
                Err(NurbsError::Structure { .. })
            ));
        });
        with_curve_cx(false, |cx| {
            let CurveAdmissionRun::Complete { admitted } = curve
                .admit_with_cx(cx)
                .expect("healthy cancellable admission")
            else {
                panic!("active context must admit the valid curve");
            };
            assert!(core::ptr::eq(admitted.source(), &curve));
            assert!(matches!(
                admitted
                    .eval_with_cx(0.5, cx)
                    .expect("admitted cancellable evaluation"),
                CurveEvaluationRun::Complete { .. }
            ));
        });

        let mut malformed = line_curve();
        malformed.cpw.clear();
        let legacy_error = malformed.admit().expect_err("malformed legacy admission");
        with_curve_cx(false, |cx| {
            assert_eq!(
                malformed
                    .admit_with_cx(cx)
                    .expect_err("malformed cancellable admission"),
                legacy_error
            );
        });
    }

    #[test]
    fn curve_admission_replays_inside_controls_and_at_publication() {
        let long_curve = long_linear_curve();
        let run = || {
            let mut polls = 0usize;
            let mut should_cancel = || {
                polls += 1;
                polls == 13
            };
            let outcome = long_curve
                .validate_live_structure_with_poll(&mut should_cancel)
                .expect("valid long curve");
            (matches!(outcome, CurveWorkRun::Cancelled), polls)
        };
        assert_eq!(run(), run());
        assert_eq!(run(), (true, 13));

        let curve = line_curve();
        let mut total_polls = 0usize;
        let mut never_cancel = || {
            total_polls += 1;
            false
        };
        assert!(matches!(
            curve
                .validate_live_structure_with_poll(&mut never_cancel)
                .expect("healthy admission"),
            CurveWorkRun::Complete(())
        ));
        assert_eq!(total_polls, 7);
        let mut replay_polls = 0usize;
        let mut cancel_at_publication = || {
            replay_polls += 1;
            replay_polls == 7
        };
        assert!(matches!(
            curve
                .validate_live_structure_with_poll(&mut cancel_at_publication)
                .expect("publication cancellation"),
            CurveWorkRun::Cancelled
        ));
        assert_eq!(replay_polls, 7);
    }

    #[test]
    fn curve_evaluation_replays_inside_accumulation_and_at_publication() {
        let degree = 16usize;
        let mut knots = vec![0.0; degree + 1];
        knots.extend(vec![1.0; degree + 1]);
        #[allow(clippy::cast_precision_loss)]
        let controls: Vec<[f64; 1]> = (0..=degree)
            .map(|index| [index as f64 / degree as f64])
            .collect();
        let weights = vec![1.0; degree + 1];
        let curve = NurbsCurve::new(
            KnotVector::new(knots, degree).expect("high-degree knots"),
            &controls,
            &weights,
        )
        .expect("high-degree curve");
        let admitted = curve.admit().expect("admitted high-degree curve");
        let (span, basis) = admitted.knots().basis(0.5).expect("basis row");
        let run = || {
            let mut polls = 0usize;
            let outcome = admitted
                .eval_from_basis_with_poll(span, &basis, || {
                    polls += 1;
                    polls == 2
                })
                .expect("finite accumulation");
            (outcome, polls)
        };
        assert_eq!(run(), run());
        assert_eq!(run(), (CurveEvaluationRun::Cancelled, 2));

        let homogeneous_run = || {
            let mut polls = 0usize;
            let outcome = admitted
                .eval_homogeneous_from_basis_with_poll(span, &basis, || {
                    polls += 1;
                    polls == 2
                })
                .expect("finite homogeneous accumulation");
            (outcome, polls)
        };
        assert_eq!(homogeneous_run(), homogeneous_run());
        assert_eq!(
            homogeneous_run(),
            (CurveHomogeneousEvaluationRun::Cancelled, 2)
        );

        let line = line_curve();
        let admitted_line = line.admit().expect("admitted line");
        let (line_span, line_basis) = admitted_line.knots().basis(0.5).expect("line basis");
        let mut homogeneous_total_polls = 0usize;
        assert!(matches!(
            admitted_line
                .eval_homogeneous_from_basis_with_poll(line_span, &line_basis, || {
                    homogeneous_total_polls += 1;
                    false
                })
                .expect("healthy homogeneous line evaluation"),
            CurveHomogeneousEvaluationRun::Complete { .. }
        ));
        assert_eq!(homogeneous_total_polls, 2);
        let mut homogeneous_replay_polls = 0usize;
        assert_eq!(
            admitted_line
                .eval_homogeneous_from_basis_with_poll(line_span, &line_basis, || {
                    homogeneous_replay_polls += 1;
                    homogeneous_replay_polls == homogeneous_total_polls
                })
                .expect("homogeneous publication cancellation"),
            CurveHomogeneousEvaluationRun::Cancelled
        );
        assert_eq!(homogeneous_replay_polls, homogeneous_total_polls);

        let mut full_homogeneous_polls = 0usize;
        let mut never_cancel = || {
            full_homogeneous_polls += 1;
            false
        };
        assert!(matches!(
            admitted_line
                .eval_homogeneous_with_poll(0.5, &mut never_cancel)
                .expect("healthy full homogeneous evaluation"),
            CurveHomogeneousEvaluationRun::Complete { .. }
        ));
        assert_eq!(full_homogeneous_polls, 7);
        let mut full_homogeneous_replay = 0usize;
        let mut cancel_at_homogeneous_publication = || {
            full_homogeneous_replay += 1;
            full_homogeneous_replay == full_homogeneous_polls
        };
        assert_eq!(
            admitted_line
                .eval_homogeneous_with_poll(0.5, &mut cancel_at_homogeneous_publication)
                .expect("full homogeneous publication cancellation"),
            CurveHomogeneousEvaluationRun::Cancelled
        );
        assert_eq!(full_homogeneous_replay, full_homogeneous_polls);

        let mut total_polls = 0usize;
        assert!(matches!(
            admitted_line
                .eval_from_basis_with_poll(line_span, &line_basis, || {
                    total_polls += 1;
                    false
                })
                .expect("healthy line evaluation"),
            CurveEvaluationRun::Complete { .. }
        ));
        assert_eq!(total_polls, 2);
        let mut replay_polls = 0usize;
        assert_eq!(
            admitted_line
                .eval_from_basis_with_poll(line_span, &line_basis, || {
                    replay_polls += 1;
                    replay_polls == total_polls
                })
                .expect("publication cancellation"),
            CurveEvaluationRun::Cancelled
        );
        assert_eq!(replay_polls, total_polls);
    }

    #[test]
    fn bezier_planning_and_known_span_insertion_cancel_inside_linear_work() {
        let long_curve = long_linear_curve();
        let admitted_long = long_curve.admit().expect("admitted long curve");
        let plan_run = || {
            let mut polls = 0usize;
            let mut should_cancel = || {
                polls += 1;
                polls == 2
            };
            let outcome = plan_bezier_conversion_with_poll(admitted_long, &mut should_cancel)
                .expect("bounded Bezier plan");
            (matches!(outcome, CurveWorkRun::Cancelled), polls)
        };
        assert_eq!(plan_run(), plan_run());
        assert_eq!(plan_run(), (true, 2));

        let insertion_curve = high_degree_insertion_curve();
        let admitted_insertion = insertion_curve.admit().expect("admitted insertion curve");
        let insertion_run = || {
            let mut polls = 0usize;
            let mut should_cancel = || {
                polls += 1;
                polls == 2
            };
            let outcome = insertion_curve
                .insert_knot_at_span_with_poll(0.5, 17, &mut should_cancel)
                .expect("valid known-span insertion");
            (matches!(outcome, CurveWorkRun::Cancelled), polls)
        };
        assert_eq!(insertion_run(), insertion_run());
        assert_eq!(insertion_run(), (true, 2));

        let mut never_cancel = || false;
        let CurveWorkRun::Complete(cancellable) = insertion_curve
            .insert_knot_at_span_with_poll(0.5, 17, &mut never_cancel)
            .expect("healthy known-span insertion")
        else {
            panic!("healthy known-span insertion must complete");
        };
        assert_eq!(
            cancellable,
            admitted_insertion
                .insert_knot(0.5)
                .expect("legacy span-resolving insertion")
        );
    }

    #[test]
    fn knot_insertion_final_checkpoint_gates_curve_publication() {
        let curve = line_curve();
        let mut total_polls = 0usize;
        let mut never_cancel = || {
            total_polls += 1;
            false
        };
        let complete = curve
            .insert_knot_at_span_with_poll(0.5, 1, &mut never_cancel)
            .expect("healthy known-span insertion");
        assert!(matches!(complete, CurveWorkRun::Complete(_)));

        let mut replay_polls = 0usize;
        let mut cancel_at_publication = || {
            replay_polls += 1;
            replay_polls == total_polls
        };
        let cancelled = curve
            .insert_knot_at_span_with_poll(0.5, 1, &mut cancel_at_publication)
            .expect("publication cancellation");
        assert!(matches!(cancelled, CurveWorkRun::Cancelled));
        assert_eq!(replay_polls, total_polls);
    }

    #[test]
    fn bezier_conversion_with_cx_is_transactional_and_publication_gated() {
        let half = Rat::new(1, 2);
        let knots = KnotVector::new(
            vec![
                Rat::int(0),
                Rat::int(0),
                Rat::int(0),
                half,
                Rat::int(1),
                Rat::int(1),
                Rat::int(1),
            ],
            2,
        )
        .expect("quadratic knots");
        let curve = NurbsCurve::<Rat, 1>::new(
            knots,
            &[[Rat::int(0)], [Rat::int(1)], [Rat::int(2)], [Rat::int(3)]],
            &[Rat::int(1); 4],
        )
        .expect("quadratic curve");
        let admitted = curve.admit().expect("admitted curve");
        let synchronous = curve.to_bezier_form().expect("synchronous conversion");

        with_curve_cx(true, |cx| {
            assert!(matches!(
                curve
                    .to_bezier_form_with_cx(cx)
                    .expect("pre-cancelled owning conversion"),
                CurveBezierRun::Cancelled
            ));
            assert!(matches!(
                admitted
                    .to_bezier_form_with_cx(cx)
                    .expect("pre-cancelled conversion"),
                CurveBezierRun::Cancelled
            ));
        });
        with_curve_cx(false, |cx| {
            let CurveBezierRun::Complete {
                curve: owning_converted,
            } = curve
                .to_bezier_form_with_cx(cx)
                .expect("healthy cancellable owning conversion")
            else {
                panic!("active context must complete owning conversion");
            };
            assert_eq!(owning_converted, synchronous);
            let CurveBezierRun::Complete { curve: converted } = admitted
                .to_bezier_form_with_cx(cx)
                .expect("healthy cancellable conversion")
            else {
                panic!("active context must complete conversion");
            };
            assert_eq!(
                converted,
                admitted
                    .insert_knot(half)
                    .expect("exact reference insertion")
            );
        });

        let publication_curve = line_curve();
        let mut total_polls = 0usize;
        let mut never_cancel = || {
            total_polls += 1;
            false
        };
        let complete = publication_curve
            .to_bezier_form_after_validation_with_poll(&mut never_cancel)
            .expect("healthy conversion");
        assert!(matches!(complete, CurveBezierRun::Complete { .. }));
        assert_eq!(total_polls, 8);
        let mut replay_polls = 0usize;
        let mut cancel_at_publication = || {
            replay_polls += 1;
            replay_polls == total_polls
        };
        let cancelled = publication_curve
            .to_bezier_form_after_validation_with_poll(&mut cancel_at_publication)
            .expect("publication cancellation");
        assert!(matches!(cancelled, CurveBezierRun::Cancelled));
        assert_eq!(replay_polls, 8);

        let mut admission_polls = 0usize;
        let mut observe_admission = || {
            admission_polls += 1;
            false
        };
        assert!(matches!(
            publication_curve
                .validate_live_structure_with_poll(&mut observe_admission)
                .expect("healthy source admission"),
            CurveWorkRun::Complete(())
        ));
        let mut owning_polls = 0usize;
        let mut cancel_at_first_conversion_poll = || {
            owning_polls += 1;
            owning_polls == admission_polls + 1
        };
        assert!(matches!(
            publication_curve
                .to_bezier_form_with_poll(&mut cancel_at_first_conversion_poll)
                .expect("owning conversion cancellation"),
            CurveBezierRun::Cancelled
        ));
        assert_eq!(owning_polls, admission_polls + 1);

        let invalid_dimension = NurbsCurve::<f64, 4> {
            knots: KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("line knots"),
            cpw: vec![[0.0, 0.0, 0.0, 1.0]; 2],
        };
        with_curve_cx(true, |cx| {
            assert!(matches!(
                invalid_dimension.to_bezier_form_with_cx(cx),
                Err(NurbsError::Structure { .. })
            ));
        });
    }

    #[test]
    fn admitted_bezier_conversion_is_exact_and_preflighted() {
        let knots = KnotVector::new(
            vec![
                Rat::int(0),
                Rat::int(0),
                Rat::int(0),
                Rat::new(1, 2),
                Rat::int(1),
                Rat::int(1),
                Rat::int(1),
            ],
            2,
        )
        .expect("quadratic knots");
        let curve = NurbsCurve::<Rat, 1>::new(
            knots,
            &[[Rat::int(0)], [Rat::int(1)], [Rat::int(2)], [Rat::int(3)]],
            &[Rat::int(1); 4],
        )
        .expect("quadratic curve");
        let admitted = curve.admit().expect("admitted curve");
        let plan = admitted
            .bezier_conversion_plan()
            .expect("checked conversion plan");
        assert_eq!(plan.insertions, 1);
        assert_eq!(plan.distinct_knot_count, 3);
        assert_eq!(plan.final_knot_count, 8);
        assert_eq!(plan.final_control_count, 5);
        assert!(plan.work_units > 0);
        assert!(plan.converted_bytes > 0);
        assert!(plan.peak_allocated_bytes >= plan.converted_bytes);

        let bezier = admitted.to_bezier_form().expect("Bezier conversion");
        for parameter in [Rat::int(0), Rat::new(1, 4), Rat::new(3, 4), Rat::int(1)] {
            assert_eq!(
                admitted.eval(parameter).expect("source evaluation"),
                bezier.eval(parameter).expect("converted evaluation")
            );
        }
        assert_eq!(
            bezier
                .admit()
                .expect("admitted converted curve")
                .span_boxes()
                .expect("admitted span boxes")
                .len(),
            2
        );
    }

    #[test]
    fn degree_elevation_with_cx_is_transactional_and_exact() {
        let curve = quadratic_join_curve();
        let admitted = curve.admit().expect("admitted quadratic join");
        let expected = admitted.elevate_degree().expect("synchronous elevation");
        assert_eq!(expected, curve.elevate_degree().expect("owning elevation"));

        with_curve_cx(true, |cx| {
            assert_eq!(
                admitted
                    .elevate_degree_with_cx(cx)
                    .expect("pre-cancelled elevation"),
                CurveElevationRun::Cancelled
            );
        });
        with_curve_cx(false, |cx| {
            assert_eq!(
                admitted
                    .elevate_degree_with_cx(cx)
                    .expect("active elevation"),
                CurveElevationRun::Complete {
                    curve: expected.try_clone().expect("expected elevation copy"),
                }
            );
        });

        let zero = Rat::int(0);
        let one = Rat::int(1);
        let exact = NurbsCurve::<Rat, 1>::new(
            KnotVector::new(vec![zero, zero, zero, Rat::new(1, 2), one, one, one], 2)
                .expect("exact quadratic knots"),
            &[[zero], [Rat::int(1)], [Rat::int(2)], [Rat::int(3)]],
            &[one; 4],
        )
        .expect("exact quadratic curve");
        let exact_admitted = exact.admit().expect("admitted exact curve");
        let exact_expected = exact_admitted
            .elevate_degree()
            .expect("exact synchronous elevation");
        with_curve_cx(false, |cx| {
            assert_eq!(
                exact_admitted
                    .elevate_degree_with_cx(cx)
                    .expect("active exact elevation"),
                CurveElevationRun::Complete {
                    curve: exact_expected,
                }
            );
        });

        let full_break = linear_full_break_curve();
        let full_break_elevated = full_break
            .admit()
            .expect("admitted full break")
            .elevate_degree()
            .expect("full-break elevation");
        assert_eq!(full_break_elevated.knots.degree(), 2);
        assert_eq!(
            full_break_elevated
                .knots
                .knots()
                .iter()
                .filter(|&&knot| knot == 0.5)
                .count(),
            3
        );
        assert_eq!(full_break_elevated.cpw.len(), 6);
        assert_eq!(full_break_elevated.cpw[2], [0.25, 0.0, 0.0, 1.0]);
        assert_eq!(full_break_elevated.cpw[3], [0.75, 0.0, 0.0, 1.0]);
    }

    #[test]
    fn degree_elevation_cancellation_replays_every_checkpoint() {
        let curve = long_linear_curve();
        let admitted = curve.admit().expect("admitted long curve");
        let expected = admitted.elevate_degree().expect("reference elevation");

        let complete_run = || {
            let mut polls = 0usize;
            let mut never_cancel = || {
                polls += 1;
                false
            };
            let outcome = admitted
                .elevate_degree_with_poll(&mut never_cancel)
                .expect("healthy cancellable elevation");
            (outcome, polls)
        };
        let (complete, total_polls) = complete_run();
        assert_eq!(complete, CurveElevationRun::Complete { curve: expected });
        assert_eq!(complete_run().1, total_polls);
        assert!(
            total_polls > 16,
            "long elevation must expose phase checkpoints"
        );

        for cancel_at in 1..=total_polls {
            let mut polls = 0usize;
            let mut should_cancel = || {
                polls += 1;
                polls == cancel_at
            };
            assert_eq!(
                admitted
                    .elevate_degree_with_poll(&mut should_cancel)
                    .expect("valid cancellation replay"),
                CurveElevationRun::Cancelled,
                "checkpoint {cancel_at} published a partial or complete curve"
            );
            assert_eq!(polls, cancel_at);
        }
    }

    #[test]
    fn degree_elevation_plan_counts_single_join_and_full_break_shapes() {
        let cases = [
            (line_curve(), 2, 1, 2, 3, 6, 3),
            (quadratic_join_curve(), 3, 2, 3, 4, 11, 7),
            (linear_full_break_curve(), 3, 2, 2, 3, 9, 6),
        ];
        for (
            curve,
            distinct_knot_count,
            segment_count,
            elevated_degree,
            elevated_order,
            final_knot_count,
            final_control_count,
        ) in cases
        {
            let admitted = curve.admit().expect("admitted elevation source");
            let bezier = plan_bezier_conversion(admitted).expect("Bezier plan");
            let plan = plan_curve_elevation::<f64>(admitted.knots().degree(), bezier)
                .expect("degree-elevation plan");
            assert_eq!(plan.distinct_knot_count, distinct_knot_count);
            assert_eq!(plan.segment_count, segment_count);
            assert_eq!(plan.elevated_degree, elevated_degree);
            assert_eq!(plan.elevated_order, elevated_order);
            assert_eq!(plan.final_knot_count, final_knot_count);
            assert_eq!(plan.final_control_count, final_control_count);

            let elevated = curve.elevate_degree().expect("bounded degree elevation");
            assert_eq!(elevated.knots.knots().len(), plan.final_knot_count);
            assert_eq!(elevated.cpw.len(), plan.final_control_count);
        }
    }

    #[test]
    fn degree_elevation_plan_matches_work_and_peak_live_formulas() {
        let curve = quadratic_join_curve();
        let admitted = curve.admit().expect("admitted quadratic join");
        let bezier = plan_bezier_conversion(admitted).expect("Bezier plan");
        let plan = plan_curve_elevation::<f64>(admitted.knots().degree(), bezier)
            .expect("degree-elevation plan");

        let expected_work = bezier.work_units
            + 4 * bezier.final_knot_count as u128
            + 4 * bezier.final_control_count as u128
            + 32 * admitted.knots().degree() as u128 * plan.segment_count as u128
            + 4 * plan.final_control_count as u128
            + plan.final_knot_count as u128
            + 2 * (16 * plan.final_knot_count as u128 + plan.elevated_degree as u128)
            + 16 * plan.final_control_count as u128
            + 32;
        assert_eq!(plan.work_units, expected_work);

        let metadata_bytes = plan.distinct_knot_count as u128
            * (core::mem::size_of::<f64>() + core::mem::size_of::<usize>()) as u128;
        let assembly_bytes = bezier.converted_bytes
            + metadata_bytes
            + curve_storage_bytes::<f64>(plan.final_knot_count, plan.final_control_count)
                .expect("elevated storage bytes");
        assert_eq!(
            plan.peak_retained_bytes,
            bezier.peak_allocated_bytes.max(assembly_bytes)
        );
    }

    #[test]
    fn degree_elevation_envelope_is_boundary_exact_and_work_first() {
        assert!(
            enforce_curve_elevation_envelope(
                BASIS_MAX_WORK_UNITS,
                CURVE_ELEVATION_MAX_RETAINED_BYTES - 1,
            )
            .is_ok()
        );
        assert!(
            enforce_curve_elevation_envelope(
                BASIS_MAX_WORK_UNITS,
                CURVE_ELEVATION_MAX_RETAINED_BYTES,
            )
            .is_ok()
        );
        assert!(matches!(
            enforce_curve_elevation_envelope(
                BASIS_MAX_WORK_UNITS,
                CURVE_ELEVATION_MAX_RETAINED_BYTES + 1,
            ),
            Err(NurbsError::Domain { ref what }) if what.contains("retain")
        ));
        assert!(matches!(
            enforce_curve_elevation_envelope(
                BASIS_MAX_WORK_UNITS + 1,
                CURVE_ELEVATION_MAX_RETAINED_BYTES + 1,
            ),
            Err(NurbsError::Domain { ref what }) if what.contains("work")
        ));

        let pre_memory_refusal = BezierConversionPlan {
            insertions: 0,
            distinct_knot_count: 2,
            final_knot_count: 4,
            final_control_count: 2,
            work_units: BASIS_MAX_WORK_UNITS,
            peak_allocated_bytes: u128::MAX,
            converted_bytes: u128::MAX,
        };
        assert!(matches!(
            plan_curve_elevation::<f64>(1, pre_memory_refusal),
            Err(NurbsError::Domain { ref what }) if what.contains("work")
        ));
    }

    #[test]
    fn span_box_preflight_prices_scan_work_and_retained_output() {
        assert_eq!(
            preflight_span_boxes(2, 1, core::mem::size_of::<SpanBox<f64, 3>>())
                .expect("one linear span"),
            1
        );
        assert!(
            enforce_span_box_envelope(BASIS_MAX_WORK_UNITS, CURVE_SPAN_BOX_MAX_RETAINED_BYTES)
                .is_ok()
        );
        assert!(matches!(
            enforce_span_box_envelope(
                BASIS_MAX_WORK_UNITS + 1,
                CURVE_SPAN_BOX_MAX_RETAINED_BYTES
            ),
            Err(NurbsError::Domain { ref what }) if what.contains("work")
        ));
        assert!(matches!(
            enforce_span_box_envelope(
                BASIS_MAX_WORK_UNITS,
                CURVE_SPAN_BOX_MAX_RETAINED_BYTES + 1
            ),
            Err(NurbsError::Domain { ref what }) if what.contains("retain")
        ));
        assert!(matches!(
            preflight_span_boxes(600_000, 1, core::mem::size_of::<SpanBox<f64, 3>>()),
            Err(NurbsError::Domain { ref what }) if what.contains("traversal")
        ));
    }

    #[test]
    fn knot_insertion_preflight_prices_work_and_retained_output() {
        let curve = line_curve();
        let plan = plan_curve_insertion(curve.admit().expect("admitted line"))
            .expect("line insertion plan");
        assert_eq!(plan.new_knot_count, 5);
        assert_eq!(plan.new_control_count, 3);
        assert!(
            enforce_curve_insertion_envelope(
                BASIS_MAX_WORK_UNITS,
                CURVE_INSERTION_MAX_RETAINED_BYTES,
            )
            .is_ok()
        );
        assert!(matches!(
            enforce_curve_insertion_envelope(
                BASIS_MAX_WORK_UNITS + 1,
                CURVE_INSERTION_MAX_RETAINED_BYTES,
            ),
            Err(NurbsError::Domain { ref what }) if what.contains("work")
        ));
        assert!(matches!(
            enforce_curve_insertion_envelope(
                BASIS_MAX_WORK_UNITS,
                CURVE_INSERTION_MAX_RETAINED_BYTES + 1,
            ),
            Err(NurbsError::Domain { ref what }) if what.contains("retain")
        ));
    }

    #[test]
    fn derivative_parameter_refusal_precedes_live_structure_scan() {
        let knots = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("line knots");
        let malformed = NurbsCurve::<f64, 1> {
            knots,
            cpw: Vec::new(),
        };
        let error = malformed
            .derivatives(-1.0, 1)
            .expect_err("out-of-domain parameter must refuse before malformed controls");
        assert!(matches!(
            error,
            NurbsError::Domain { ref what } if what.contains("parameter")
        ));
    }

    #[test]
    fn derivative_envelope_and_fallible_copy_are_stable() {
        assert_eq!(
            NurbsCurve::<f64, 3>::derivative_envelope(2, 4, 1, 1)
                .expect("linear derivative envelope"),
            (44, 304 + 2 * core::mem::size_of::<Vec<[f64; 4]>>() as u128,)
        );
        let knots = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("line knots");
        let curve =
            NurbsCurve::<f64, 1>::new(knots, &[[0.0], [1.0]], &[1.0, 1.0]).expect("line curve");
        assert_eq!(curve.try_clone().expect("fallible curve copy"), curve);
    }
}
