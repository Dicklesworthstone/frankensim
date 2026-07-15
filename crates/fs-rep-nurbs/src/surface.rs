//! Tensor-product rational surfaces: two-stage de Boor evaluation
//! (generic scalar), EXACT directional knot insertion (Boehm along rows
//! or columns), first partial derivatives (f64), and per-span control
//! boxes (the convex-hull property in both directions).

use crate::NurbsError;
use crate::basis::{
    BASIS_MAX_WORK_UNITS, BasisRun, KnotSpanRun, KnotValidationOutcome, KnotVector, Scalar,
};
use crate::curve::{CurveDerivativesRun, NurbsCurve};
use fs_exec::Cx;

const SURFACE_VALIDATION_WORK_PER_CONTROL: u128 = 16;
// Keep this aligned with KnotVector's private conservative validation price.
const SURFACE_KNOT_VALIDATION_WORK_PER_ENTRY: u128 = 16;
// Match the curve refinement envelope's conservative four-lane blend price;
// copied controls are deliberately charged at the same worst-case rate.
const SURFACE_INSERTION_WORK_PER_CONTROL: u128 = 32;
const SURFACE_SPAN_BOX_WORK_PER_CONTROL: u128 = 16;
const SURFACE_SPAN_BOX_MAX_RETAINED_BYTES: u128 = 64 * 1024 * 1024;
const SURFACE_COPY_MAX_RETAINED_BYTES: u128 = 64 * 1024 * 1024;
const SURFACE_INSERTION_MAX_RETAINED_BYTES: u128 = 64 * 1024 * 1024;
const SURFACE_CANCELLATION_STRIDE: usize = 64;

fn surface_poll_due(
    operations_since_poll: &mut usize,
    should_cancel: &mut impl FnMut() -> bool,
) -> bool {
    *operations_since_poll += 1;
    if *operations_since_poll < SURFACE_CANCELLATION_STRIDE {
        return false;
    }
    *operations_since_poll = 0;
    should_cancel()
}

fn preflight_surface_copy<S: Scalar>(
    knot_count_u: usize,
    knot_count_v: usize,
    control_count_u: usize,
    control_count_v: usize,
) -> Result<(), NurbsError> {
    let control_count = (control_count_u as u128)
        .checked_mul(control_count_v as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "surface-copy control-count accounting overflows u128".to_string(),
        })?;
    let work_units = control_count
        .checked_mul(4)
        .and_then(|work| work.checked_add(knot_count_u as u128))
        .and_then(|work| work.checked_add(knot_count_v as u128))
        .and_then(|work| work.checked_add(control_count_u as u128))
        .and_then(|work| work.checked_add(4))
        .ok_or_else(|| NurbsError::Domain {
            what: "surface-copy work accounting overflows u128".to_string(),
        })?;
    if work_units > BASIS_MAX_WORK_UNITS {
        return Err(NurbsError::Domain {
            what: format!(
                "surface copy requests {work_units} work units above defensive ceiling {BASIS_MAX_WORK_UNITS}"
            ),
        });
    }

    let knot_bytes = (knot_count_u as u128)
        .checked_add(knot_count_v as u128)
        .and_then(|count| count.checked_mul(core::mem::size_of::<S>() as u128))
        .ok_or_else(|| NurbsError::Domain {
            what: "surface-copy knot-storage accounting overflows u128".to_string(),
        })?;
    let row_table_bytes = (control_count_u as u128)
        .checked_mul(core::mem::size_of::<Vec<[S; 4]>>() as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "surface-copy row-table accounting overflows u128".to_string(),
        })?;
    let control_bytes = control_count
        .checked_mul(core::mem::size_of::<[S; 4]>() as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "surface-copy control-storage accounting overflows u128".to_string(),
        })?;
    let retained_bytes = knot_bytes
        .checked_add(row_table_bytes)
        .and_then(|bytes| bytes.checked_add(control_bytes))
        .ok_or_else(|| NurbsError::Domain {
            what: "surface-copy retained-byte accounting overflows u128".to_string(),
        })?;
    if retained_bytes > SURFACE_COPY_MAX_RETAINED_BYTES {
        return Err(NurbsError::Domain {
            what: format!(
                "surface copy retains {retained_bytes} output bytes above defensive ceiling {SURFACE_COPY_MAX_RETAINED_BYTES}"
            ),
        });
    }
    Ok(())
}

/// One (u-span × v-span) control box: (min, max, (u0, u1), (v0, v1)).
pub type SurfaceSpanBox<S> = ([S; 3], [S; 3], (S, S), (S, S));

/// Value + first partials `(S, S_u, S_v)`.
pub type SurfacePartials = ([f64; 3], [f64; 3], [f64; 3]);

fn enforce_span_box_envelope(work: u128, retained_bytes: u128) -> Result<(), NurbsError> {
    if work > BASIS_MAX_WORK_UNITS {
        return Err(NurbsError::Domain {
            what: format!(
                "surface span-box traversal requests {work} work units above defensive ceiling {BASIS_MAX_WORK_UNITS}"
            ),
        });
    }
    if retained_bytes > SURFACE_SPAN_BOX_MAX_RETAINED_BYTES {
        return Err(NurbsError::Domain {
            what: format!(
                "surface span boxes retain {retained_bytes} bytes above defensive ceiling {SURFACE_SPAN_BOX_MAX_RETAINED_BYTES}"
            ),
        });
    }
    Ok(())
}

fn preflight_span_boxes(
    control_count_u: usize,
    control_count_v: usize,
    degree_u: usize,
    degree_v: usize,
    retained_bytes_per_box: usize,
) -> Result<usize, NurbsError> {
    let span_count_u =
        control_count_u
            .checked_sub(degree_u)
            .ok_or_else(|| NurbsError::Structure {
                what: "surface u degree exceeds its admitted control count".to_string(),
            })?;
    let span_count_v =
        control_count_v
            .checked_sub(degree_v)
            .ok_or_else(|| NurbsError::Structure {
                what: "surface v degree exceeds its admitted control count".to_string(),
            })?;
    let span_capacity =
        span_count_u
            .checked_mul(span_count_v)
            .ok_or_else(|| NurbsError::Domain {
                what: "surface span-box count overflows usize".to_string(),
            })?;
    let order_u = degree_u.checked_add(1).ok_or_else(|| NurbsError::Domain {
        what: "surface u order overflows usize during span-box admission".to_string(),
    })?;
    let order_v = degree_v.checked_add(1).ok_or_else(|| NurbsError::Domain {
        what: "surface v order overflows usize during span-box admission".to_string(),
    })?;
    let control_visits = (span_capacity as u128)
        .checked_mul(order_u as u128)
        .and_then(|work| work.checked_mul(order_v as u128))
        .ok_or_else(|| NurbsError::Domain {
            what: "surface span-box control-scan work overflows u128".to_string(),
        })?;
    // Worst case: Su outer span checks, two checks/write-accounting units per
    // candidate box, and a conservative 16 units for each overlapping control
    // visit (three Cartesian projections plus comparisons).
    let traversal_work =
        (span_count_u as u128)
            .checked_add((span_capacity as u128).checked_mul(2).ok_or_else(|| {
                NurbsError::Domain {
                    what: "surface span-box candidate work overflows u128".to_string(),
                }
            })?)
            .and_then(|work| {
                control_visits
                    .checked_mul(SURFACE_SPAN_BOX_WORK_PER_CONTROL)
                    .and_then(|control_work| work.checked_add(control_work))
            })
            .ok_or_else(|| NurbsError::Domain {
                what: "surface span-box traversal work overflows u128".to_string(),
            })?;
    let retained_bytes = (span_capacity as u128)
        .checked_mul(retained_bytes_per_box as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "surface span-box retained-byte accounting overflows u128".to_string(),
        })?;
    enforce_span_box_envelope(traversal_work, retained_bytes)?;
    Ok(span_capacity)
}

fn basis_operation_work(control_count: usize, degree: usize) -> Result<u128, NurbsError> {
    let order = degree.checked_add(1).ok_or_else(|| NurbsError::Domain {
        what: "surface-partial basis order overflows usize".to_string(),
    })?;
    (degree as u128)
        .checked_mul(order as u128)
        .map(|product| product / 2)
        .and_then(|work| work.checked_add(order as u128))
        .and_then(|work| work.checked_add(control_count as u128))
        .ok_or_else(|| NurbsError::Domain {
            what: "surface-partial basis work overflows u128".to_string(),
        })
}

fn enforce_partials_envelope(work: u128, retained_bytes: u128) -> Result<(), NurbsError> {
    if work > NurbsCurve::<f64, 3>::MAX_DERIVATIVE_WORK_UNITS {
        return Err(NurbsError::Domain {
            what: format!(
                "surface partials request {work} work units above defensive ceiling {}",
                NurbsCurve::<f64, 3>::MAX_DERIVATIVE_WORK_UNITS
            ),
        });
    }
    if retained_bytes > NurbsCurve::<f64, 3>::MAX_DERIVATIVE_RETAINED_BYTES {
        return Err(NurbsError::Domain {
            what: format!(
                "surface partials retain up to {retained_bytes} bytes above defensive ceiling {}",
                NurbsCurve::<f64, 3>::MAX_DERIVATIVE_RETAINED_BYTES
            ),
        });
    }
    Ok(())
}

fn preflight_partials_envelope(
    control_count_u: usize,
    control_count_v: usize,
    knot_count_u: usize,
    knot_count_v: usize,
    degree_u: usize,
    degree_v: usize,
) -> Result<(), NurbsError> {
    let order_u = degree_u.checked_add(1).ok_or_else(|| NurbsError::Domain {
        what: "surface-partial u order overflows usize".to_string(),
    })?;
    let order_v = degree_v.checked_add(1).ok_or_else(|| NurbsError::Domain {
        what: "surface-partial v order overflows usize".to_string(),
    })?;
    let contractions = (control_count_u as u128)
        .checked_mul(order_v as u128)
        .and_then(|u_visits| {
            (control_count_v as u128)
                .checked_mul(order_u as u128)
                .and_then(|v_visits| u_visits.checked_add(v_visits))
        })
        .and_then(|visits| visits.checked_mul(8))
        .ok_or_else(|| NurbsError::Domain {
            what: "surface-partial contraction work overflows u128".to_string(),
        })?;
    let (derivative_work_u, derivative_bytes_u) =
        NurbsCurve::<f64, 3>::derivative_envelope(control_count_u, knot_count_u, degree_u, 1)?;
    let (derivative_work_v, derivative_bytes_v) =
        NurbsCurve::<f64, 3>::derivative_envelope(control_count_v, knot_count_v, degree_v, 1)?;
    let work = basis_operation_work(control_count_u, degree_u)?
        .checked_add(basis_operation_work(control_count_v, degree_v)?)
        .and_then(|total| total.checked_add(contractions))
        // Conservatively charge one full knot extent per axis for the bounded
        // multiplicity lookups that prove ordinary first-derivative existence.
        .and_then(|total| total.checked_add(knot_count_u as u128))
        .and_then(|total| total.checked_add(knot_count_v as u128))
        .and_then(|total| total.checked_add(derivative_work_u))
        .and_then(|total| total.checked_add(derivative_work_v))
        .ok_or_else(|| NurbsError::Domain {
            what: "surface-partial aggregate work overflows u128".to_string(),
        })?;
    let scalar_bytes = core::mem::size_of::<f64>() as u128;
    let control_bytes = core::mem::size_of::<[f64; 4]>() as u128;
    let basis_bytes = (order_u as u128)
        .checked_add(order_v as u128)
        .and_then(|count| count.checked_mul(scalar_bytes))
        .ok_or_else(|| NurbsError::Domain {
            what: "surface-partial retained-basis accounting overflows u128".to_string(),
        })?;
    let u_peak = (control_count_u as u128)
        .checked_mul(control_bytes)
        .and_then(|bytes| bytes.checked_add(derivative_bytes_u))
        .and_then(|bytes| bytes.checked_add(basis_bytes))
        .ok_or_else(|| NurbsError::Domain {
            what: "surface-partial u-stage retained-byte accounting overflows u128".to_string(),
        })?;
    let retained_u_jet = 2u128
        .checked_mul(core::mem::size_of::<[f64; 3]>() as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "surface-partial retained u-jet accounting overflows u128".to_string(),
        })?;
    let v_peak = (control_count_v as u128)
        .checked_mul(control_bytes)
        .and_then(|bytes| bytes.checked_add(derivative_bytes_v))
        .and_then(|bytes| bytes.checked_add(basis_bytes))
        .and_then(|bytes| bytes.checked_add(retained_u_jet))
        .ok_or_else(|| NurbsError::Domain {
            what: "surface-partial v-stage retained-byte accounting overflows u128".to_string(),
        })?;
    let retained_bytes = u_peak.max(v_peak);
    enforce_partials_envelope(work, retained_bytes)
}

/// A rational tensor-product surface in 3D.
///
/// ```compile_fail
/// use fs_rep_nurbs::{KnotVector, NurbsSurface};
/// let knots = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).unwrap();
/// let mut surface = NurbsSurface::new(
///     knots.clone(), knots,
///     &vec![vec![[0.0, 0.0, 0.0]; 2]; 2],
///     &vec![vec![1.0; 2]; 2],
/// ).unwrap();
/// surface.cpw.clear();
/// ```
#[derive(Debug, PartialEq)]
pub struct NurbsSurface<S: Scalar> {
    /// Knots in u.
    pub(crate) knots_u: KnotVector<S>,
    /// Knots in v.
    pub(crate) knots_v: KnotVector<S>,
    /// Homogeneous control net `cpw[i][j]`, i along u, j along v.
    pub(crate) cpw: Vec<Vec<[S; 4]>>,
}

/// A validate-once borrow of one exact immutable surface snapshot.
#[derive(Debug, Clone, Copy)]
pub struct AdmittedNurbsSurface<'a, S: Scalar> {
    inner: &'a NurbsSurface<S>,
}

/// Transactional terminal state of a cancellation-aware fallible surface copy.
#[must_use]
#[derive(Debug, PartialEq)]
pub enum SurfaceCloneRun<S: Scalar> {
    /// The complete sealed copy of the exact source representation.
    Complete {
        /// Copied surface generation.
        surface: NurbsSurface<S>,
    },
    /// Cancellation was observed; all partial nested copy storage was dropped.
    Cancelled,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SurfaceValidationOutcome {
    Complete,
    Cancelled,
}

/// Transactional terminal state of cancellation-aware surface admission.
#[must_use]
#[derive(Debug, Clone, Copy)]
pub enum SurfaceAdmissionRun<'a, S: Scalar> {
    /// The exact immutable source snapshot was fully validated.
    Complete {
        /// Lifetime-bound authority for the validated surface generation.
        admitted: AdmittedNurbsSurface<'a, S>,
    },
    /// Cancellation was observed; no admitted authority was published.
    Cancelled,
}

/// Transactional terminal state of cancellation-aware Cartesian evaluation.
#[must_use]
#[derive(Debug, Clone, PartialEq)]
pub enum SurfaceEvaluationRun<S: Scalar> {
    /// The complete finite Cartesian point.
    Complete {
        /// Evaluated point on the admitted surface.
        point: [S; 3],
    },
    /// Cancellation was observed; no partial point was published.
    Cancelled,
}

/// Transactional terminal state of cancellation-aware f64 surface partials.
#[must_use]
#[derive(Debug, Clone, PartialEq)]
pub enum SurfacePartialsRun {
    /// The complete value and first directional partials.
    Complete {
        /// `(S, S_u, S_v)` at the requested parameter pair.
        partials: SurfacePartials,
    },
    /// Cancellation was observed; no value or directional jet was published.
    Cancelled,
}

/// Transactional terminal state of cancellation-aware directional insertion.
#[must_use]
#[derive(Debug, PartialEq)]
pub enum SurfaceInsertionRun<S: Scalar> {
    /// The complete sealed and validated derived surface.
    Complete {
        /// Exact directional refinement of the admitted source surface.
        surface: NurbsSurface<S>,
    },
    /// Cancellation was observed; all partial derived storage was dropped.
    Cancelled,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SurfaceInsertionAxis {
    U,
    V,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SurfaceInsertionPlan {
    axis: SurfaceInsertionAxis,
    new_control_count_u: usize,
    new_control_count_v: usize,
    new_knot_count_u: usize,
    new_knot_count_v: usize,
    #[cfg(test)]
    work_units: u128,
    #[cfg(test)]
    retained_bytes: u128,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct SurfaceInsertionEnvelope {
    work_units: u128,
    #[cfg(test)]
    retained_bytes: u128,
}

#[derive(Debug)]
enum SurfaceWorkRun<T> {
    Complete(T),
    Cancelled,
}

/// Transactional terminal state of cancellation-aware surface span boxes.
#[must_use]
#[derive(Debug, Clone, PartialEq)]
pub enum SurfaceSpanBoxesRun<S: Scalar> {
    /// Every admitted nonempty span pair has a complete Cartesian control box.
    Complete {
        /// Boxes in deterministic U-major, V-minor source-span order.
        boxes: Vec<SurfaceSpanBox<S>>,
    },
    /// Cancellation was observed; no partial box table was published.
    Cancelled,
}

fn enforce_surface_insertion_envelope(
    work_units: u128,
    retained_bytes: u128,
) -> Result<(), NurbsError> {
    if work_units > BASIS_MAX_WORK_UNITS {
        return Err(NurbsError::Domain {
            what: format!(
                "surface knot insertion requests {work_units} work units above defensive ceiling {BASIS_MAX_WORK_UNITS}"
            ),
        });
    }
    if retained_bytes > SURFACE_INSERTION_MAX_RETAINED_BYTES {
        return Err(NurbsError::Domain {
            what: format!(
                "surface knot insertion retains {retained_bytes} output bytes above defensive ceiling {SURFACE_INSERTION_MAX_RETAINED_BYTES}"
            ),
        });
    }
    Ok(())
}

fn surface_knot_validation_work(
    knot_count: usize,
    degree: usize,
    axis: &str,
) -> Result<u128, NurbsError> {
    (knot_count as u128)
        .checked_mul(SURFACE_KNOT_VALIDATION_WORK_PER_ENTRY)
        .and_then(|work| work.checked_add(degree as u128))
        .ok_or_else(|| NurbsError::Domain {
            what: format!("surface inserted {axis}-knot validation work overflows u128"),
        })
}

#[allow(clippy::too_many_arguments)]
fn surface_insertion_envelope_for_result<S: Scalar>(
    axis: SurfaceInsertionAxis,
    degree_u: usize,
    degree_v: usize,
    new_control_count_u: usize,
    new_control_count_v: usize,
    new_knot_count_u: usize,
    new_knot_count_v: usize,
) -> Result<SurfaceInsertionEnvelope, NurbsError> {
    let final_control_count = (new_control_count_u as u128)
        .checked_mul(new_control_count_v as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "surface inserted control count overflows u128".to_string(),
        })?;
    let validation_u = surface_knot_validation_work(new_knot_count_u, degree_u, "u")?;
    let validation_v = surface_knot_validation_work(new_knot_count_v, degree_v, "v")?;
    let final_validation_work = validation_u
        .checked_add(validation_v)
        .and_then(|work| {
            final_control_count
                .checked_mul(SURFACE_VALIDATION_WORK_PER_CONTROL)
                .and_then(|controls| work.checked_add(controls))
        })
        .ok_or_else(|| NurbsError::Domain {
            what: "surface inserted validation work overflows u128".to_string(),
        })?;
    let inserted_knot_validation_work = match axis {
        SurfaceInsertionAxis::U => validation_u,
        SurfaceInsertionAxis::V => validation_v,
    };
    let span_work = match axis {
        SurfaceInsertionAxis::U => new_control_count_u.checked_sub(1),
        SurfaceInsertionAxis::V => new_control_count_v.checked_sub(1),
    }
    .ok_or_else(|| NurbsError::Domain {
        what: "surface inserted direction has no source controls".to_string(),
    })? as u128;
    let assembly_work = final_control_count
        .checked_mul(SURFACE_INSERTION_WORK_PER_CONTROL)
        .ok_or_else(|| NurbsError::Domain {
            what: "surface insertion control-assembly work overflows u128".to_string(),
        })?;
    let knot_copy_work = (new_knot_count_u as u128)
        .checked_add(new_knot_count_v as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "surface insertion knot-copy work overflows u128".to_string(),
        })?;
    let work_units = span_work
        .checked_add(assembly_work)
        .and_then(|work| work.checked_add(knot_copy_work))
        .and_then(|work| work.checked_add(new_control_count_u as u128))
        .and_then(|work| work.checked_add(inserted_knot_validation_work))
        .and_then(|work| work.checked_add(final_validation_work))
        .and_then(|work| work.checked_add(32))
        .ok_or_else(|| NurbsError::Domain {
            what: "surface insertion aggregate work overflows u128".to_string(),
        })?;

    let knot_bytes = (new_knot_count_u as u128)
        .checked_add(new_knot_count_v as u128)
        .and_then(|count| count.checked_mul(core::mem::size_of::<S>() as u128))
        .ok_or_else(|| NurbsError::Domain {
            what: "surface insertion knot payload overflows u128".to_string(),
        })?;
    let row_table_bytes = (new_control_count_u as u128)
        .checked_mul(core::mem::size_of::<Vec<[S; 4]>>() as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "surface insertion row-table payload overflows u128".to_string(),
        })?;
    let control_bytes = final_control_count
        .checked_mul(core::mem::size_of::<[S; 4]>() as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "surface insertion control payload overflows u128".to_string(),
        })?;
    let retained_bytes = knot_bytes
        .checked_add(row_table_bytes)
        .and_then(|bytes| bytes.checked_add(control_bytes))
        .ok_or_else(|| NurbsError::Domain {
            what: "surface insertion retained payload overflows u128".to_string(),
        })?;
    enforce_surface_insertion_envelope(work_units, retained_bytes)?;
    Ok(SurfaceInsertionEnvelope {
        work_units,
        #[cfg(test)]
        retained_bytes,
    })
}

fn copy_surface_knot_vector_with_poll<S: Scalar>(
    source: &KnotVector<S>,
    insertion: Option<(usize, S)>,
    output_count: usize,
    stage: &str,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<SurfaceWorkRun<KnotVector<S>>, NurbsError> {
    if should_cancel() {
        return Ok(SurfaceWorkRun::Cancelled);
    }
    let mut entries = Vec::new();
    entries
        .try_reserve_exact(output_count)
        .map_err(|_| NurbsError::Domain {
            what: format!("surface {stage} knot allocation was refused"),
        })?;
    let mut operations_since_poll = 0usize;
    if let Some((span, t)) = insertion {
        for &knot in &source.knots()[..=span] {
            entries.push(knot);
            if surface_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(SurfaceWorkRun::Cancelled);
            }
        }
        entries.push(t);
        if surface_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(SurfaceWorkRun::Cancelled);
        }
        for &knot in &source.knots()[span + 1..] {
            entries.push(knot);
            if surface_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(SurfaceWorkRun::Cancelled);
            }
        }
    } else {
        for &knot in source.knots() {
            entries.push(knot);
            if surface_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(SurfaceWorkRun::Cancelled);
            }
        }
    }
    if should_cancel() {
        return Ok(SurfaceWorkRun::Cancelled);
    }
    Ok(SurfaceWorkRun::Complete(KnotVector {
        knots: entries,
        degree: source.degree(),
    }))
}

impl<S: Scalar> NurbsSurface<S> {
    fn validation_work_for(
        knots_u: &KnotVector<S>,
        knots_v: &KnotVector<S>,
    ) -> Result<u128, NurbsError> {
        let controls = (knots_u.control_count() as u128)
            .checked_mul(knots_v.control_count() as u128)
            .and_then(|count| count.checked_mul(SURFACE_VALIDATION_WORK_PER_CONTROL))
            .ok_or_else(|| NurbsError::Domain {
                what: "surface control-validation work overflows u128".to_string(),
            })?;
        knots_u
            .validation_work()?
            .checked_add(knots_v.validation_work()?)
            .and_then(|work| work.checked_add(controls))
            .ok_or_else(|| NurbsError::Domain {
                what: "surface structure-validation work overflows u128".to_string(),
            })
    }

    fn enforce_validation_work(work: u128) -> Result<(), NurbsError> {
        if work > BASIS_MAX_WORK_UNITS {
            return Err(NurbsError::Domain {
                what: format!(
                    "surface structure validation requests {work} work units above defensive ceiling {BASIS_MAX_WORK_UNITS}"
                ),
            });
        }
        Ok(())
    }

    fn insertion_plan_after_parameter(
        &self,
        t: S,
        axis: SurfaceInsertionAxis,
    ) -> Result<SurfaceInsertionPlan, NurbsError> {
        let direction_knots = match axis {
            SurfaceInsertionAxis::U => self.knots_u.admitted_after_validation(),
            SurfaceInsertionAxis::V => self.knots_v.admitted_after_validation(),
        };
        let (lo, hi) = direction_knots.domain();
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

        let old_control_count_u = self.knots_u.control_count();
        let old_control_count_v = self.knots_v.control_count();
        let old_knot_count_u = self.knots_u.knots().len();
        let old_knot_count_v = self.knots_v.knots().len();
        let (add_u, add_v) = match axis {
            SurfaceInsertionAxis::U => (1usize, 0usize),
            SurfaceInsertionAxis::V => (0usize, 1usize),
        };
        let new_control_count_u =
            old_control_count_u
                .checked_add(add_u)
                .ok_or_else(|| NurbsError::Domain {
                    what: "surface inserted u-control count overflows usize".to_string(),
                })?;
        let new_control_count_v =
            old_control_count_v
                .checked_add(add_v)
                .ok_or_else(|| NurbsError::Domain {
                    what: "surface inserted v-control count overflows usize".to_string(),
                })?;
        let new_knot_count_u =
            old_knot_count_u
                .checked_add(add_u)
                .ok_or_else(|| NurbsError::Domain {
                    what: "surface inserted u-knot count overflows usize".to_string(),
                })?;
        let new_knot_count_v =
            old_knot_count_v
                .checked_add(add_v)
                .ok_or_else(|| NurbsError::Domain {
                    what: "surface inserted v-knot count overflows usize".to_string(),
                })?;
        let envelope = surface_insertion_envelope_for_result::<S>(
            axis,
            self.knots_u.degree(),
            self.knots_v.degree(),
            new_control_count_u,
            new_control_count_v,
            new_knot_count_u,
            new_knot_count_v,
        )?;
        #[cfg(not(test))]
        let _ = envelope;
        Ok(SurfaceInsertionPlan {
            axis,
            new_control_count_u,
            new_control_count_v,
            new_knot_count_u,
            new_knot_count_v,
            #[cfg(test)]
            work_units: envelope.work_units,
            #[cfg(test)]
            retained_bytes: envelope.retained_bytes,
        })
    }

    fn validate_live_structure(&self) -> Result<(), NurbsError> {
        let mut never_cancel = || false;
        match self.validate_live_structure_with_poll(&mut never_cancel)? {
            SurfaceValidationOutcome::Complete => Ok(()),
            SurfaceValidationOutcome::Cancelled => Err(NurbsError::Domain {
                what: "non-cancelling surface validation observed cancellation".to_string(),
            }),
        }
    }

    fn validate_live_structure_with_poll(
        &self,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<SurfaceValidationOutcome, NurbsError> {
        Self::enforce_validation_work(Self::validation_work_for(&self.knots_u, &self.knots_v)?)?;
        match self.knots_u.validate_live_with_poll(|| should_cancel())? {
            KnotValidationOutcome::Complete => {}
            KnotValidationOutcome::Cancelled => return Ok(SurfaceValidationOutcome::Cancelled),
        }
        match self.knots_v.validate_live_with_poll(|| should_cancel())? {
            KnotValidationOutcome::Complete => {}
            KnotValidationOutcome::Cancelled => return Ok(SurfaceValidationOutcome::Cancelled),
        }

        let invalid_control_net = || NurbsError::Structure {
            what: concat!(
                "live surface control net must match both knot vectors and retain finite ",
                "homogeneous coordinates with admissible weights"
            )
            .to_string(),
        };
        let expected_u = self.knots_u.control_count();
        let expected_v = self.knots_v.control_count();
        if self.cpw.len() != expected_u {
            return Err(invalid_control_net());
        }
        if should_cancel() {
            return Ok(SurfaceValidationOutcome::Cancelled);
        }

        let mut operations_since_poll = 0usize;
        for row in &self.cpw {
            if row.len() != expected_v {
                return Err(invalid_control_net());
            }
            if surface_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(SurfaceValidationOutcome::Cancelled);
            }
            for control in row {
                if !control[3].is_admissible_weight() {
                    return Err(invalid_control_net());
                }
                if surface_poll_due(&mut operations_since_poll, should_cancel) {
                    return Ok(SurfaceValidationOutcome::Cancelled);
                }
                for &component in control {
                    if !component.is_finite() {
                        return Err(invalid_control_net());
                    }
                    if surface_poll_due(&mut operations_since_poll, should_cancel) {
                        return Ok(SurfaceValidationOutcome::Cancelled);
                    }
                }
                for &component in &control[..3] {
                    if !component.quotient_is_finite(control[3]) {
                        return Err(invalid_control_net());
                    }
                    if surface_poll_due(&mut operations_since_poll, should_cancel) {
                        return Ok(SurfaceValidationOutcome::Cancelled);
                    }
                }
            }
        }
        if should_cancel() {
            return Ok(SurfaceValidationOutcome::Cancelled);
        }
        Ok(SurfaceValidationOutcome::Complete)
    }

    /// Build from Cartesian control net + weights.
    ///
    /// # Errors
    /// [`NurbsError::Structure`] on grid/count mismatches, non-finite
    /// coordinates, or non-positive/non-finite weights; [`NurbsError::Domain`]
    /// when validation work or a control-net allocation is refused.
    pub fn new(
        knots_u: KnotVector<S>,
        knots_v: KnotVector<S>,
        points: &[Vec<[S; 3]>],
        weights: &[Vec<S>],
    ) -> Result<Self, NurbsError> {
        Self::enforce_validation_work(Self::validation_work_for(&knots_u, &knots_v)?)?;
        knots_u.validate_live()?;
        knots_v.validate_live()?;
        let (nu, nv) = (knots_u.control_count(), knots_v.control_count());
        if points.len() != nu || weights.len() != nu {
            return Err(NurbsError::Structure {
                what: format!("expected {nu} control rows, got {}", points.len()),
            });
        }
        // Scan the complete borrowed input before allocating any output row.
        // A malformed late row must not force a large, doomed prefix allocation.
        for (prow, wrow) in points.iter().zip(weights) {
            if prow.len() != nv || wrow.len() != nv {
                return Err(NurbsError::Structure {
                    what: format!("expected {nv} control columns"),
                });
            }
            if wrow.iter().copied().any(|w| !w.is_admissible_weight()) {
                return Err(NurbsError::Structure {
                    what: "weights must be finite, positive, and numerically admissible"
                        .to_string(),
                });
            }
            if prow
                .iter()
                .flat_map(|point| point.iter())
                .copied()
                .any(|coordinate| !coordinate.is_finite())
            {
                return Err(NurbsError::Structure {
                    what: "control-point coordinates must be finite".to_string(),
                });
            }
            for (point, &weight) in prow.iter().zip(wrow) {
                let homogeneous = [
                    point[0] * weight,
                    point[1] * weight,
                    point[2] * weight,
                    weight,
                ];
                if point
                    .iter()
                    .zip(homogeneous.iter())
                    .any(|(&coordinate, &weighted)| {
                        coordinate != S::zero() && weighted == S::zero()
                    })
                {
                    return Err(NurbsError::Structure {
                        what:
                            "Cartesian coordinate × weight underflowed a nonzero coordinate to zero"
                                .to_string(),
                    });
                }
                if homogeneous
                    .iter()
                    .copied()
                    .any(|component| !component.is_finite())
                {
                    return Err(NurbsError::Structure {
                        what: "Cartesian coordinate × weight overflowed the homogeneous numeric domain"
                            .to_string(),
                    });
                }
            }
        }
        let mut cpw = Vec::new();
        cpw.try_reserve_exact(nu).map_err(|_| NurbsError::Domain {
            what: "surface control-row allocation was refused".to_string(),
        })?;
        for (prow, wrow) in points.iter().zip(weights) {
            let mut row = Vec::new();
            row.try_reserve_exact(nv).map_err(|_| NurbsError::Domain {
                what: "surface homogeneous-control row allocation was refused".to_string(),
            })?;
            for (p, &w) in prow.iter().zip(wrow) {
                let homogeneous = [p[0] * w, p[1] * w, p[2] * w, w];
                row.push(homogeneous);
            }
            cpw.push(row);
        }
        Ok(NurbsSurface {
            knots_u,
            knots_v,
            cpw,
        })
    }

    /// Build from a homogeneous control net, validating the complete sealed
    /// representation before exposing it.
    ///
    /// # Errors
    /// [`NurbsError::Structure`] for invalid knots, grid shape, coordinates,
    /// Cartesian projections, or weights; [`NurbsError::Domain`] when
    /// validation work is refused.
    pub fn from_homogeneous(
        knots_u: KnotVector<S>,
        knots_v: KnotVector<S>,
        cpw: Vec<Vec<[S; 4]>>,
    ) -> Result<Self, NurbsError> {
        let candidate = NurbsSurface {
            knots_u,
            knots_v,
            cpw,
        };
        candidate.validate_live_structure()?;
        Ok(candidate)
    }

    /// Borrow the sealed u knot vector.
    #[must_use]
    pub const fn knots_u(&self) -> &KnotVector<S> {
        &self.knots_u
    }

    /// Borrow the sealed v knot vector.
    #[must_use]
    pub const fn knots_v(&self) -> &KnotVector<S> {
        &self.knots_v
    }

    /// Borrow the sealed homogeneous control net.
    #[must_use]
    pub fn homogeneous_control_net(&self) -> &[Vec<[S; 4]>] {
        &self.cpw
    }

    /// Fallibly copy this sealed surface without revalidating unchanged data.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] when checked copy work/retained bytes or a
    /// destination allocation is refused.
    pub fn try_clone(&self) -> Result<Self, NurbsError> {
        let mut never_cancel = || false;
        match self
            .admitted_after_validation()
            .try_clone_with_poll(&mut never_cancel)?
        {
            Some(surface) => Ok(surface),
            None => Err(NurbsError::Domain {
                what: "non-cancelling surface copy observed cancellation".to_string(),
            }),
        }
    }

    /// Fallibly copy this sealed surface with bounded cancellation polling.
    ///
    /// Count-derived work and a 64 MiB retained-output envelope precede
    /// cancellation. One gate then covers both knot copies, the outer row
    /// table and every inner control-row allocation/copy in deterministic
    /// row-major order, and final publication. The borrowed source is excluded
    /// from the output envelope. Individual allocator calls, scalar/array
    /// copies, and nested-vector destruction are not preemptible. This
    /// primitive does not consume the `Cx` budget or own request -> drain ->
    /// finalize semantics.
    ///
    /// # Errors
    /// Returns the synchronous copy's work, retained-memory, or allocation
    /// refusal when it wins before an observed cancellation.
    pub fn try_clone_with_cx(&self, cx: &Cx<'_>) -> Result<SurfaceCloneRun<S>, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        match self
            .admitted_after_validation()
            .try_clone_with_poll(&mut should_cancel)?
        {
            Some(surface) => Ok(SurfaceCloneRun::Complete { surface }),
            None => Ok(SurfaceCloneRun::Cancelled),
        }
    }

    /// Validate this exact immutable surface snapshot once.
    ///
    /// # Errors
    /// [`NurbsError::Structure`] when the sealed source is internally invalid;
    /// [`NurbsError::Domain`] when validation work exceeds the defensive cap.
    pub fn admit(&self) -> Result<AdmittedNurbsSurface<'_, S>, NurbsError> {
        self.validate_live_structure()?;
        Ok(AdmittedNurbsSurface { inner: self })
    }

    /// Bind an admitted view to a generation already validated by an internal
    /// transactional constructor.
    ///
    /// Compound crate-internal algorithms use this only while carrying the
    /// invariant from an admitted source through fully validated derived
    /// publications. Public callers must use [`Self::admit`].
    pub(crate) const fn admitted_after_validation(&self) -> AdmittedNurbsSurface<'_, S> {
        AdmittedNurbsSurface { inner: self }
    }

    /// Validate this immutable surface with bounded cancellation polling.
    ///
    /// Checked validation-work refusal retains synchronous precedence. The U
    /// knot, V knot, and row-major homogeneous-control scans share one gate,
    /// and a final checkpoint gates publication of the lifetime-bound admitted
    /// view. This method does not make construction cancellation-aware and
    /// does not consume the `Cx` budget or finalize its executor scope.
    ///
    /// # Errors
    /// Returns the synchronous admission's work, knot, control-grid, weight,
    /// and finite-arithmetic refusals when they win before an observed
    /// cancellation.
    pub fn admit_with_cx<'a>(
        &'a self,
        cx: &Cx<'_>,
    ) -> Result<SurfaceAdmissionRun<'a, S>, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        match self.validate_live_structure_with_poll(&mut should_cancel)? {
            SurfaceValidationOutcome::Complete => Ok(SurfaceAdmissionRun::Complete {
                admitted: AdmittedNurbsSurface { inner: self },
            }),
            SurfaceValidationOutcome::Cancelled => Ok(SurfaceAdmissionRun::Cancelled),
        }
    }

    /// Homogeneous evaluation.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] outside the domain.
    pub fn eval_homogeneous(&self, u: S, v: S) -> Result<[S; 4], NurbsError> {
        self.admit()?.eval_homogeneous(u, v)
    }

    /// Cartesian evaluation.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] outside the domain.
    pub fn eval(&self, u: S, v: S) -> Result<[S; 3], NurbsError> {
        self.admit()?.eval(u, v)
    }

    fn assemble_inserted_control_net_with_poll(
        &self,
        t: S,
        span: usize,
        plan: SurfaceInsertionPlan,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<SurfaceWorkRun<Vec<Vec<[S; 4]>>>, NurbsError> {
        if should_cancel() {
            return Ok(SurfaceWorkRun::Cancelled);
        }
        let mut output = Vec::new();
        output
            .try_reserve_exact(plan.new_control_count_u)
            .map_err(|_| NurbsError::Domain {
                what: "surface insertion row-table allocation was refused".to_string(),
            })?;
        let mut operations_since_poll = 0usize;

        match plan.axis {
            SurfaceInsertionAxis::U => {
                for _ in 0..plan.new_control_count_u {
                    if should_cancel() {
                        return Ok(SurfaceWorkRun::Cancelled);
                    }
                    let mut row = Vec::new();
                    row.try_reserve_exact(plan.new_control_count_v)
                        .map_err(|_| NurbsError::Domain {
                            what: "surface u-insertion output-row allocation was refused"
                                .to_string(),
                        })?;
                    output.push(row);
                }
                let degree = self.knots_u.degree();
                for column in 0..self.knots_v.control_count() {
                    for source_row in 0..=span - degree {
                        output[source_row].push(self.cpw[source_row][column]);
                        if surface_poll_due(&mut operations_since_poll, should_cancel) {
                            return Ok(SurfaceWorkRun::Cancelled);
                        }
                    }
                    for source_row in (span - degree + 1)..=span {
                        let denominator = self.knots_u.knots()[source_row + degree]
                            - self.knots_u.knots()[source_row];
                        let alpha = (t - self.knots_u.knots()[source_row]) / denominator;
                        let mut blended = [S::zero(); 4];
                        for ((slot, &left), &right) in blended
                            .iter_mut()
                            .zip(&self.cpw[source_row - 1][column])
                            .zip(&self.cpw[source_row][column])
                        {
                            *slot = (S::one() - alpha) * left + alpha * right;
                            if surface_poll_due(&mut operations_since_poll, should_cancel) {
                                return Ok(SurfaceWorkRun::Cancelled);
                            }
                        }
                        output[source_row].push(blended);
                        if surface_poll_due(&mut operations_since_poll, should_cancel) {
                            return Ok(SurfaceWorkRun::Cancelled);
                        }
                    }
                    for source_row in span..self.knots_u.control_count() {
                        output[source_row + 1].push(self.cpw[source_row][column]);
                        if surface_poll_due(&mut operations_since_poll, should_cancel) {
                            return Ok(SurfaceWorkRun::Cancelled);
                        }
                    }
                }
            }
            SurfaceInsertionAxis::V => {
                let degree = self.knots_v.degree();
                for source in &self.cpw {
                    if should_cancel() {
                        return Ok(SurfaceWorkRun::Cancelled);
                    }
                    let mut row = Vec::new();
                    row.try_reserve_exact(plan.new_control_count_v)
                        .map_err(|_| NurbsError::Domain {
                            what: "surface v-insertion output-row allocation was refused"
                                .to_string(),
                        })?;
                    for &control in &source[..=span - degree] {
                        row.push(control);
                        if surface_poll_due(&mut operations_since_poll, should_cancel) {
                            return Ok(SurfaceWorkRun::Cancelled);
                        }
                    }
                    for source_column in (span - degree + 1)..=span {
                        let denominator = self.knots_v.knots()[source_column + degree]
                            - self.knots_v.knots()[source_column];
                        let alpha = (t - self.knots_v.knots()[source_column]) / denominator;
                        let mut blended = [S::zero(); 4];
                        for ((slot, &left), &right) in blended
                            .iter_mut()
                            .zip(&source[source_column - 1])
                            .zip(&source[source_column])
                        {
                            *slot = (S::one() - alpha) * left + alpha * right;
                            if surface_poll_due(&mut operations_since_poll, should_cancel) {
                                return Ok(SurfaceWorkRun::Cancelled);
                            }
                        }
                        row.push(blended);
                        if surface_poll_due(&mut operations_since_poll, should_cancel) {
                            return Ok(SurfaceWorkRun::Cancelled);
                        }
                    }
                    for &control in &source[span..] {
                        row.push(control);
                        if surface_poll_due(&mut operations_since_poll, should_cancel) {
                            return Ok(SurfaceWorkRun::Cancelled);
                        }
                    }
                    output.push(row);
                    if surface_poll_due(&mut operations_since_poll, should_cancel) {
                        return Ok(SurfaceWorkRun::Cancelled);
                    }
                }
            }
        }
        if should_cancel() {
            return Ok(SurfaceWorkRun::Cancelled);
        }
        Ok(SurfaceWorkRun::Complete(output))
    }

    fn insert_knot_at_span_with_plan_and_poll(
        &self,
        t: S,
        span: usize,
        plan: SurfaceInsertionPlan,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<SurfaceInsertionRun<S>, NurbsError> {
        let cpw =
            match self.assemble_inserted_control_net_with_poll(t, span, plan, should_cancel)? {
                SurfaceWorkRun::Complete(cpw) => cpw,
                SurfaceWorkRun::Cancelled => return Ok(SurfaceInsertionRun::Cancelled),
            };

        let (inserted_source, unchanged_source, inserted_count, unchanged_count, stage) =
            match plan.axis {
                SurfaceInsertionAxis::U => (
                    &self.knots_u,
                    &self.knots_v,
                    plan.new_knot_count_u,
                    plan.new_knot_count_v,
                    "u-insertion",
                ),
                SurfaceInsertionAxis::V => (
                    &self.knots_v,
                    &self.knots_u,
                    plan.new_knot_count_v,
                    plan.new_knot_count_u,
                    "v-insertion",
                ),
            };
        let inserted = match copy_surface_knot_vector_with_poll(
            inserted_source,
            Some((span, t)),
            inserted_count,
            stage,
            should_cancel,
        )? {
            SurfaceWorkRun::Complete(knots) => knots,
            SurfaceWorkRun::Cancelled => return Ok(SurfaceInsertionRun::Cancelled),
        };
        KnotVector::<S>::enforce_work(
            inserted.validation_work()?,
            "surface inserted knot-vector validation",
        )?;
        match inserted.validate_live_with_poll(|| should_cancel())? {
            KnotValidationOutcome::Complete => {}
            KnotValidationOutcome::Cancelled => return Ok(SurfaceInsertionRun::Cancelled),
        }
        let unchanged = match copy_surface_knot_vector_with_poll(
            unchanged_source,
            None,
            unchanged_count,
            "unchanged-axis",
            should_cancel,
        )? {
            SurfaceWorkRun::Complete(knots) => knots,
            SurfaceWorkRun::Cancelled => return Ok(SurfaceInsertionRun::Cancelled),
        };
        let (knots_u, knots_v) = match plan.axis {
            SurfaceInsertionAxis::U => (inserted, unchanged),
            SurfaceInsertionAxis::V => (unchanged, inserted),
        };
        let candidate = NurbsSurface {
            knots_u,
            knots_v,
            cpw,
        };
        match candidate.validate_live_structure_with_poll(should_cancel)? {
            SurfaceValidationOutcome::Complete => {}
            SurfaceValidationOutcome::Cancelled => return Ok(SurfaceInsertionRun::Cancelled),
        }
        if should_cancel() {
            return Ok(SurfaceInsertionRun::Cancelled);
        }
        Ok(SurfaceInsertionRun::Complete { surface: candidate })
    }

    /// EXACT knot insertion in the u direction (Boehm on every column).
    ///
    /// # Errors
    /// [`NurbsError::Domain`] for a non-interior parameter.
    pub fn insert_knot_u(&self, t: S) -> Result<Self, NurbsError> {
        self.admit()?.insert_knot_u(t)
    }

    /// EXACT knot insertion in the v direction.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] for a non-interior parameter.
    pub fn insert_knot_v(&self, t: S) -> Result<Self, NurbsError> {
        self.admit()?.insert_knot_v(t)
    }

    /// Per-(u-span × v-span) Cartesian control boxes: each patch of the
    /// surface lies inside its sub-net's box (convex-hull property).
    ///
    /// # Errors
    /// Returns [`NurbsError::Structure`] when the sealed representation does
    /// not satisfy its knot/control-net invariants, or [`NurbsError::Domain`]
    /// when validation work or output allocation is refused.
    pub fn span_boxes(&self) -> Result<Vec<SurfaceSpanBox<S>>, NurbsError> {
        self.admit()?.span_boxes()
    }
}

impl<'a, S: Scalar> AdmittedNurbsSurface<'a, S> {
    /// The exact immutable source bound to this view.
    #[must_use]
    pub const fn source(&self) -> &'a NurbsSurface<S> {
        self.inner
    }

    /// Borrow the admitted u knot vector without rescanning it.
    #[must_use]
    pub fn knots_u(&self) -> crate::basis::AdmittedKnotVector<'a, S> {
        self.inner.knots_u.admitted_after_validation()
    }

    /// Borrow the admitted v knot vector without rescanning it.
    #[must_use]
    pub fn knots_v(&self) -> crate::basis::AdmittedKnotVector<'a, S> {
        self.inner.knots_v.admitted_after_validation()
    }

    /// Borrow the sealed homogeneous control net.
    #[must_use]
    pub fn homogeneous_control_net(&self) -> &'a [Vec<[S; 4]>] {
        &self.inner.cpw
    }

    /// Conservatively price a sequence of directional insertions at the
    /// largest derived generation, without scanning spans or allocating.
    ///
    /// The same per-insertion work and retained-output envelope used by the
    /// executable insertion path is enforced here so compound callers cannot
    /// admit a conversion that a nested insertion would later refuse.
    pub(crate) fn projected_directional_insertion_work(
        &self,
        insertions_u: usize,
        insertions_v: usize,
    ) -> Result<u128, NurbsError> {
        let final_control_count_u = self
            .knots_u()
            .control_count()
            .checked_add(insertions_u)
            .ok_or_else(|| NurbsError::Domain {
                what: "projected surface u-control count overflows usize".to_string(),
            })?;
        let final_control_count_v = self
            .knots_v()
            .control_count()
            .checked_add(insertions_v)
            .ok_or_else(|| NurbsError::Domain {
                what: "projected surface v-control count overflows usize".to_string(),
            })?;
        let final_knot_count_u = self
            .knots_u()
            .knots()
            .len()
            .checked_add(insertions_u)
            .ok_or_else(|| NurbsError::Domain {
                what: "projected surface u-knot count overflows usize".to_string(),
            })?;
        let final_knot_count_v = self
            .knots_v()
            .knots()
            .len()
            .checked_add(insertions_v)
            .ok_or_else(|| NurbsError::Domain {
                what: "projected surface v-knot count overflows usize".to_string(),
            })?;

        let mut work_units = 0u128;
        for (axis, insertions) in [
            (SurfaceInsertionAxis::U, insertions_u),
            (SurfaceInsertionAxis::V, insertions_v),
        ] {
            if insertions == 0 {
                continue;
            }
            let envelope = surface_insertion_envelope_for_result::<S>(
                axis,
                self.knots_u().degree(),
                self.knots_v().degree(),
                final_control_count_u,
                final_control_count_v,
                final_knot_count_u,
                final_knot_count_v,
            )?;
            work_units = (insertions as u128)
                .checked_mul(envelope.work_units)
                .and_then(|axis_work| work_units.checked_add(axis_work))
                .ok_or_else(|| NurbsError::Domain {
                    what: "projected surface insertion work overflows u128".to_string(),
                })?;
        }
        Ok(work_units)
    }

    /// Fallibly copy this admitted source while sharing a compound caller's
    /// cancellation callback. `None` means no partial copy was published.
    pub(crate) fn try_clone_with_poll(
        &self,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<Option<NurbsSurface<S>>, NurbsError> {
        preflight_surface_copy::<S>(
            self.inner.knots_u.knots.len(),
            self.inner.knots_v.knots.len(),
            self.inner.knots_u.control_count(),
            self.inner.knots_v.control_count(),
        )?;
        if should_cancel() {
            return Ok(None);
        }
        let mut knot_entries_u = Vec::new();
        knot_entries_u
            .try_reserve_exact(self.inner.knots_u.knots.len())
            .map_err(|_| NurbsError::Domain {
                what: "knot-vector copy allocation was refused".to_string(),
            })?;
        if should_cancel() {
            return Ok(None);
        }
        let mut operations_since_poll = 0usize;
        for &knot in &self.inner.knots_u.knots {
            knot_entries_u.push(knot);
            if surface_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(None);
            }
        }
        let knots_u = KnotVector {
            knots: knot_entries_u,
            degree: self.inner.knots_u.degree,
        };

        if should_cancel() {
            return Ok(None);
        }
        let mut knot_entries_v = Vec::new();
        knot_entries_v
            .try_reserve_exact(self.inner.knots_v.knots.len())
            .map_err(|_| NurbsError::Domain {
                what: "knot-vector copy allocation was refused".to_string(),
            })?;
        if should_cancel() {
            return Ok(None);
        }
        for &knot in &self.inner.knots_v.knots {
            knot_entries_v.push(knot);
            if surface_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(None);
            }
        }
        let knots_v = KnotVector {
            knots: knot_entries_v,
            degree: self.inner.knots_v.degree,
        };

        if should_cancel() {
            return Ok(None);
        }
        let mut cpw = Vec::new();
        cpw.try_reserve_exact(self.inner.cpw.len())
            .map_err(|_| NurbsError::Domain {
                what: "surface copy row-table allocation was refused".to_string(),
            })?;
        if should_cancel() {
            return Ok(None);
        }
        for source_row in &self.inner.cpw {
            if should_cancel() {
                return Ok(None);
            }
            let mut row = Vec::new();
            row.try_reserve_exact(source_row.len())
                .map_err(|_| NurbsError::Domain {
                    what: "surface copy control-row allocation was refused".to_string(),
                })?;
            if should_cancel() {
                return Ok(None);
            }
            for &control in source_row {
                row.push(control);
                if surface_poll_due(&mut operations_since_poll, should_cancel) {
                    return Ok(None);
                }
            }
            cpw.push(row);
            if surface_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(None);
            }
        }
        if should_cancel() {
            return Ok(None);
        }
        Ok(Some(NurbsSurface {
            knots_u,
            knots_v,
            cpw,
        }))
    }

    fn insert_knot(&self, t: S, axis: SurfaceInsertionAxis) -> Result<NurbsSurface<S>, NurbsError> {
        let plan = self.inner.insertion_plan_after_parameter(t, axis)?;
        let span = match axis {
            SurfaceInsertionAxis::U => self.knots_u().span(t)?,
            SurfaceInsertionAxis::V => self.knots_v().span(t)?,
        };
        let mut never_cancel = || false;
        match self
            .inner
            .insert_knot_at_span_with_plan_and_poll(t, span, plan, &mut never_cancel)?
        {
            SurfaceInsertionRun::Complete { surface } => Ok(surface),
            SurfaceInsertionRun::Cancelled => Err(NurbsError::Domain {
                what: "non-cancelling surface knot insertion observed cancellation".to_string(),
            }),
        }
    }

    fn insert_knot_with_poll(
        &self,
        t: S,
        axis: SurfaceInsertionAxis,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<SurfaceInsertionRun<S>, NurbsError> {
        let plan = self.inner.insertion_plan_after_parameter(t, axis)?;
        let span = match axis {
            SurfaceInsertionAxis::U => match self.knots_u().span_with_poll(t, should_cancel)? {
                KnotSpanRun::Complete { span } => span,
                KnotSpanRun::Cancelled => return Ok(SurfaceInsertionRun::Cancelled),
            },
            SurfaceInsertionAxis::V => match self.knots_v().span_with_poll(t, should_cancel)? {
                KnotSpanRun::Complete { span } => span,
                KnotSpanRun::Cancelled => return Ok(SurfaceInsertionRun::Cancelled),
            },
        };
        self.inner
            .insert_knot_at_span_with_plan_and_poll(t, span, plan, should_cancel)
    }

    /// Insert one exact U knot while sharing a compound caller's cancellation
    /// callback.
    pub(crate) fn insert_knot_u_with_poll(
        &self,
        t: S,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<SurfaceInsertionRun<S>, NurbsError> {
        self.insert_knot_with_poll(t, SurfaceInsertionAxis::U, should_cancel)
    }

    /// Insert one exact V knot while sharing a compound caller's cancellation
    /// callback.
    pub(crate) fn insert_knot_v_with_poll(
        &self,
        t: S,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<SurfaceInsertionRun<S>, NurbsError> {
        self.insert_knot_with_poll(t, SurfaceInsertionAxis::V, should_cancel)
    }

    /// Insert one exact knot in the U direction without rescanning the source.
    ///
    /// # Errors
    /// Returns a structured parameter, work, retained-memory, allocation,
    /// numeric-domain, or derived-structure refusal.
    pub fn insert_knot_u(&self, t: S) -> Result<NurbsSurface<S>, NurbsError> {
        self.insert_knot(t, SurfaceInsertionAxis::U)
    }

    /// Insert one exact knot in the V direction without rescanning the source.
    ///
    /// # Errors
    /// Returns a structured parameter, work, retained-memory, allocation,
    /// numeric-domain, or derived-structure refusal.
    pub fn insert_knot_v(&self, t: S) -> Result<NurbsSurface<S>, NurbsError> {
        self.insert_knot(t, SurfaceInsertionAxis::V)
    }

    /// Insert one exact U knot with bounded cancellation polling.
    ///
    /// Parameter and checked aggregate work/retained-output refusals precede
    /// cancellation. The gate then covers the directional span, complete
    /// tensor Boehm assembly, both knot copies, derived validation, and final
    /// publication. Cancellation exposes no partial surface. This method does
    /// not consume the `Cx` budget or finalize its executor scope.
    ///
    /// # Errors
    /// Returns the synchronous insertion's refusal when it wins before an
    /// observed cancellation.
    pub fn insert_knot_u_with_cx(
        &self,
        t: S,
        cx: &Cx<'_>,
    ) -> Result<SurfaceInsertionRun<S>, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        self.insert_knot_u_with_poll(t, &mut should_cancel)
    }

    /// Insert one exact V knot with bounded cancellation polling.
    ///
    /// Refusal order, transactionality, budget ownership, and executor-scope
    /// boundaries are identical to [`Self::insert_knot_u_with_cx`].
    ///
    /// # Errors
    /// Returns the synchronous insertion's refusal when it wins before an
    /// observed cancellation.
    pub fn insert_knot_v_with_cx(
        &self,
        t: S,
        cx: &Cx<'_>,
    ) -> Result<SurfaceInsertionRun<S>, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        self.insert_knot_v_with_poll(t, &mut should_cancel)
    }

    /// Homogeneous evaluation without rescanning surface or knot structure.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] outside either parameter domain or when basis
    /// work/allocation is refused.
    pub fn eval_homogeneous(&self, u: S, v: S) -> Result<[S; 4], NurbsError> {
        let knots_u = self.knots_u();
        let knots_v = self.knots_v();
        let (su, bu) = knots_u.basis(u)?;
        let (sv, bv) = knots_v.basis(v)?;
        let (pu, pv) = (knots_u.degree(), knots_v.degree());
        let mut acc = [S::zero(); 4];
        for (r, &wu) in bu.iter().enumerate() {
            for (c, &wv) in bv.iter().enumerate() {
                let cp = &self.inner.cpw[su - pu + r][sv - pv + c];
                let w = wu * wv;
                for (a, &x) in acc.iter_mut().zip(cp.iter()) {
                    *a = *a + w * x;
                }
            }
        }
        if acc.iter().copied().any(|component| !component.is_finite()) {
            return Err(NurbsError::Domain {
                what: "homogeneous surface evaluation left the finite numeric domain".to_string(),
            });
        }
        Ok(acc)
    }

    /// Cartesian evaluation without rescanning the sealed snapshot.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] outside either domain or for an inadmissible
    /// rational result.
    pub fn eval(&self, u: S, v: S) -> Result<[S; 3], NurbsError> {
        let h = self.eval_homogeneous(u, v)?;
        if !h[3].is_admissible_weight() {
            return Err(NurbsError::Domain {
                what: "surface evaluation produced an inadmissible rational denominator"
                    .to_string(),
            });
        }
        let point = [h[0] / h[3], h[1] / h[3], h[2] / h[3]];
        if point
            .iter()
            .copied()
            .any(|component| !component.is_finite())
        {
            return Err(NurbsError::Domain {
                what: "Cartesian surface evaluation left the finite numeric domain".to_string(),
            });
        }
        Ok(point)
    }

    /// Evaluate a Cartesian point with bounded cancellation polling.
    ///
    /// This admitted-only entry point carries one `Cx` through both basis
    /// rows, tensor contraction, rational projection, and final publication.
    /// The caller remains responsible for owning surface admission, `Cx`
    /// budget consumption, and request -> drain -> finalize semantics.
    ///
    /// # Errors
    /// Returns the synchronous evaluator's parameter, work, allocation,
    /// weight, and finite-arithmetic refusals when they win before an observed
    /// cancellation.
    pub fn eval_with_cx(
        &self,
        u: S,
        v: S,
        cx: &Cx<'_>,
    ) -> Result<SurfaceEvaluationRun<S>, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        self.eval_with_poll(u, v, &mut should_cancel)
    }

    /// Evaluate an admitted point while sharing a compound caller's
    /// cancellation callback.
    pub(crate) fn eval_with_poll(
        &self,
        u: S,
        v: S,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<SurfaceEvaluationRun<S>, NurbsError> {
        let (span_u, basis_u) = match self.knots_u().basis_with_poll(u, should_cancel)? {
            BasisRun::Complete { span, values } => (span, values),
            BasisRun::Cancelled => return Ok(SurfaceEvaluationRun::Cancelled),
        };
        let (span_v, basis_v) = match self.knots_v().basis_with_poll(v, should_cancel)? {
            BasisRun::Complete { span, values } => (span, values),
            BasisRun::Cancelled => return Ok(SurfaceEvaluationRun::Cancelled),
        };
        self.eval_from_basis_with_poll(span_u, &basis_u, span_v, &basis_v, should_cancel)
    }

    fn eval_from_basis_with_poll(
        &self,
        span_u: usize,
        basis_u: &[S],
        span_v: usize,
        basis_v: &[S],
        mut should_cancel: impl FnMut() -> bool,
    ) -> Result<SurfaceEvaluationRun<S>, NurbsError> {
        if should_cancel() {
            return Ok(SurfaceEvaluationRun::Cancelled);
        }
        let degree_u = self.knots_u().degree();
        let degree_v = self.knots_v().degree();
        let mut operations_since_poll = 0usize;
        let mut homogeneous = [S::zero(); 4];
        for (row_offset, &weight_u) in basis_u.iter().enumerate() {
            for (column_offset, &weight_v) in basis_v.iter().enumerate() {
                let control = &self.inner.cpw[span_u - degree_u + row_offset]
                    [span_v - degree_v + column_offset];
                let tensor_weight = weight_u * weight_v;
                if surface_poll_due(&mut operations_since_poll, &mut should_cancel) {
                    return Ok(SurfaceEvaluationRun::Cancelled);
                }
                for (accumulator, &component) in homogeneous.iter_mut().zip(control) {
                    *accumulator = *accumulator + tensor_weight * component;
                    if surface_poll_due(&mut operations_since_poll, &mut should_cancel) {
                        return Ok(SurfaceEvaluationRun::Cancelled);
                    }
                }
            }
        }
        for (index, &component) in homogeneous.iter().enumerate() {
            if !component.is_finite() {
                return Err(NurbsError::Domain {
                    what: "homogeneous surface evaluation left the finite numeric domain"
                        .to_string(),
                });
            }
            if index == 3 && !component.is_admissible_weight() {
                return Err(NurbsError::Domain {
                    what: "surface evaluation produced an inadmissible rational denominator"
                        .to_string(),
                });
            }
            if surface_poll_due(&mut operations_since_poll, &mut should_cancel) {
                return Ok(SurfaceEvaluationRun::Cancelled);
            }
        }

        let point = [
            homogeneous[0] / homogeneous[3],
            homogeneous[1] / homogeneous[3],
            homogeneous[2] / homogeneous[3],
        ];
        for &component in &point {
            if !component.is_finite() {
                return Err(NurbsError::Domain {
                    what: "Cartesian surface evaluation left the finite numeric domain".to_string(),
                });
            }
            if surface_poll_due(&mut operations_since_poll, &mut should_cancel) {
                return Ok(SurfaceEvaluationRun::Cancelled);
            }
        }
        if should_cancel() {
            return Ok(SurfaceEvaluationRun::Cancelled);
        }
        Ok(SurfaceEvaluationRun::Complete { point })
    }

    /// Per-span Cartesian control boxes without a second structural scan.
    ///
    /// # Errors
    /// Returns a structured refusal when nested control scans, retained output,
    /// or the output allocation exceed the defensive legacy envelope.
    pub fn span_boxes(&self) -> Result<Vec<SurfaceSpanBox<S>>, NurbsError> {
        let mut never_cancel = || false;
        match self.span_boxes_with_poll(&mut never_cancel)? {
            SurfaceSpanBoxesRun::Complete { boxes } => Ok(boxes),
            SurfaceSpanBoxesRun::Cancelled => Err(NurbsError::Domain {
                what: "non-cancelling surface span-box traversal observed cancellation".to_string(),
            }),
        }
    }

    /// Build the complete U-major, V-minor span-box table with bounded
    /// cancellation polling and transactional publication.
    ///
    /// Checked traversal work and retained-output refusal precede cancellation.
    /// The gate then spans allocation, candidate-span traversal, Cartesian
    /// projection/bounds work, and final publication. This method does not
    /// consume the `Cx` budget or finalize its executor scope.
    ///
    /// # Errors
    /// Returns the synchronous span-box builder's work, memory, and allocation
    /// refusals when they win before an observed cancellation.
    pub fn span_boxes_with_cx(&self, cx: &Cx<'_>) -> Result<SurfaceSpanBoxesRun<S>, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        self.span_boxes_with_poll(&mut should_cancel)
    }

    fn span_boxes_with_poll(
        &self,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<SurfaceSpanBoxesRun<S>, NurbsError> {
        let knots_u = self.knots_u();
        let knots_v = self.knots_v();
        let (pu, pv) = (knots_u.degree(), knots_v.degree());
        let span_capacity = preflight_span_boxes(
            knots_u.control_count(),
            knots_v.control_count(),
            pu,
            pv,
            core::mem::size_of::<SurfaceSpanBox<S>>(),
        )?;
        if should_cancel() {
            return Ok(SurfaceSpanBoxesRun::Cancelled);
        }
        let mut out = Vec::new();
        out.try_reserve_exact(span_capacity)
            .map_err(|_| NurbsError::Domain {
                what: "surface span-box allocation was refused".to_string(),
            })?;
        let mut operations_since_poll = 0usize;
        for su in pu..knots_u.control_count() {
            let (u0, u1) = (knots_u.knots()[su], knots_u.knots()[su + 1]);
            if surface_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(SurfaceSpanBoxesRun::Cancelled);
            }
            if u1 <= u0 {
                continue;
            }
            for sv in pv..knots_v.control_count() {
                let (v0, v1) = (knots_v.knots()[sv], knots_v.knots()[sv + 1]);
                if surface_poll_due(&mut operations_since_poll, should_cancel) {
                    return Ok(SurfaceSpanBoxesRun::Cancelled);
                }
                if v1 <= v0 {
                    continue;
                }
                let mut min = [S::zero(); 3];
                let mut max = [S::zero(); 3];
                let mut first = true;
                for row in &self.inner.cpw[su - pu..=su] {
                    if surface_poll_due(&mut operations_since_poll, should_cancel) {
                        return Ok(SurfaceSpanBoxesRun::Cancelled);
                    }
                    for cp in &row[sv - pv..=sv] {
                        let w = cp[3];
                        for d in 0..3 {
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
                            if surface_poll_due(&mut operations_since_poll, should_cancel) {
                                return Ok(SurfaceSpanBoxesRun::Cancelled);
                            }
                        }
                        first = false;
                        if surface_poll_due(&mut operations_since_poll, should_cancel) {
                            return Ok(SurfaceSpanBoxesRun::Cancelled);
                        }
                    }
                }
                out.push((min, max, (u0, u1), (v0, v1)));
                if surface_poll_due(&mut operations_since_poll, should_cancel) {
                    return Ok(SurfaceSpanBoxesRun::Cancelled);
                }
            }
        }
        if should_cancel() {
            return Ok(SurfaceSpanBoxesRun::Cancelled);
        }
        Ok(SurfaceSpanBoxesRun::Complete { boxes: out })
    }
}

impl NurbsSurface<f64> {
    /// Value and first partials `(S, S_u, S_v)` at `(u, v)` via extracted
    /// isocurve nets (the standard tensor-product route).
    ///
    /// # Errors
    /// [`NurbsError::Domain`] outside the domain.
    pub fn partials(&self, u: f64, v: f64) -> Result<SurfacePartials, NurbsError> {
        self.knots_u.preflight_parameter(u, "surface u-partial")?;
        self.knots_v.preflight_parameter(v, "surface v-partial")?;
        self.admit()?.partials(u, v)
    }
}

impl AdmittedNurbsSurface<'_, f64> {
    fn preflight_partials_request(&self, u: f64, v: f64) -> Result<(), NurbsError> {
        let knots_u = self.knots_u();
        let knots_v = self.knots_v();
        self.inner
            .knots_u
            .preflight_parameter(u, "surface u-partial")?;
        self.inner
            .knots_v
            .preflight_parameter(v, "surface v-partial")?;
        preflight_partials_envelope(
            knots_u.control_count(),
            knots_v.control_count(),
            knots_u.knots().len(),
            knots_v.knots().len(),
            knots_u.degree(),
            knots_v.degree(),
        )?;
        NurbsCurve::<f64, 3>::preflight_derivative_request(knots_u, u, 1)?;
        NurbsCurve::<f64, 3>::preflight_derivative_request(knots_v, v, 1)
    }

    /// Value and first partials without rescanning the source surface.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] outside the domain or when temporary isocurve
    /// construction/evaluation is refused.
    pub fn partials(&self, u: f64, v: f64) -> Result<SurfacePartials, NurbsError> {
        // Refuse the complete request before the first basis-workspace or
        // isocurve allocation. The U-first order is deterministic when both
        // parameters are invalid.
        self.preflight_partials_request(u, v)?;
        let knots_u = self.knots_u();
        let knots_v = self.knots_v();
        let (su, bu) = knots_u.basis(u)?;
        let (sv, bv) = knots_v.basis(v)?;
        let mut never_cancel = || false;
        match self.partials_from_basis_with_poll(u, v, su, &bu, sv, &bv, &mut never_cancel)? {
            SurfacePartialsRun::Complete { partials } => Ok(partials),
            SurfacePartialsRun::Cancelled => Err(NurbsError::Domain {
                what: "non-cancelling surface partial evaluation observed cancellation".to_string(),
            }),
        }
    }

    /// Evaluate the admitted surface's value and first partials with bounded
    /// cancellation polling and transactional publication.
    ///
    /// U then V request, aggregate-envelope, and ordinary-derivative refusals
    /// retain their synchronous precedence. One `Cx` then spans both basis
    /// rows, sequential isocurve construction and differentiation, and final
    /// publication. This method does not consume the `Cx` budget or finalize
    /// its executor scope.
    ///
    /// # Errors
    /// Returns the synchronous partial evaluator's parameter, continuity,
    /// work, memory, allocation, and finite-arithmetic refusals when they win
    /// before an observed cancellation.
    pub fn partials_with_cx(
        &self,
        u: f64,
        v: f64,
        cx: &Cx<'_>,
    ) -> Result<SurfacePartialsRun, NurbsError> {
        self.preflight_partials_request(u, v)?;
        let knots_u = self.knots_u();
        let knots_v = self.knots_v();
        let (span_u, basis_u) = match knots_u.basis_with_cx(u, cx)? {
            BasisRun::Complete { span, values } => (span, values),
            BasisRun::Cancelled => return Ok(SurfacePartialsRun::Cancelled),
        };
        let (span_v, basis_v) = match knots_v.basis_with_cx(v, cx)? {
            BasisRun::Complete { span, values } => (span, values),
            BasisRun::Cancelled => return Ok(SurfacePartialsRun::Cancelled),
        };
        let mut should_cancel = || cx.checkpoint().is_err();
        self.partials_from_basis_with_poll(
            u,
            v,
            span_u,
            &basis_u,
            span_v,
            &basis_v,
            &mut should_cancel,
        )
    }

    fn partials_from_basis_with_poll(
        &self,
        u: f64,
        v: f64,
        span_u: usize,
        basis_u: &[f64],
        span_v: usize,
        basis_v: &[f64],
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<SurfacePartialsRun, NurbsError> {
        let knots_u = self.knots_u();
        let knots_v = self.knots_v();
        let degree_u = knots_u.degree();
        let degree_v = knots_v.degree();

        // u-partial: build the v-evaluated control column, differentiate
        // as a u-curve; symmetrically for v. Each scope releases its temporary
        // control net before the next one is allocated.
        let du = {
            if should_cancel() {
                return Ok(SurfacePartialsRun::Cancelled);
            }
            let mut u_net = Vec::new();
            u_net
                .try_reserve_exact(self.inner.cpw.len())
                .map_err(|_| NurbsError::Domain {
                    what: "surface u-isocurve allocation was refused".to_string(),
                })?;
            let mut operations_since_poll = 0usize;
            for row in &self.inner.cpw {
                let mut acc = [0.0f64; 4];
                for (column_offset, &weight_v) in basis_v.iter().enumerate() {
                    let control = &row[span_v - degree_v + column_offset];
                    for (accumulator, &component) in acc.iter_mut().zip(control) {
                        *accumulator += weight_v * component;
                        if surface_poll_due(&mut operations_since_poll, should_cancel) {
                            return Ok(SurfacePartialsRun::Cancelled);
                        }
                    }
                }
                for &component in &acc {
                    if !component.is_finite() {
                        return Err(NurbsError::Domain {
                            what: "surface u-isocurve left the admissible homogeneous domain"
                                .to_string(),
                        });
                    }
                    if surface_poll_due(&mut operations_since_poll, should_cancel) {
                        return Ok(SurfacePartialsRun::Cancelled);
                    }
                }
                if !acc[3].is_admissible_weight() {
                    return Err(NurbsError::Domain {
                        what: "surface u-isocurve left the admissible homogeneous domain"
                            .to_string(),
                    });
                }
                if surface_poll_due(&mut operations_since_poll, should_cancel) {
                    return Ok(SurfacePartialsRun::Cancelled);
                }
                u_net.push(acc);
            }
            match NurbsCurve::<f64, 3>::derivatives_from_admitted_parts_after_preflight_with_poll(
                knots_u,
                &u_net,
                u,
                1,
                should_cancel,
            )? {
                CurveDerivativesRun::Complete { derivatives } => derivatives,
                CurveDerivativesRun::Cancelled => return Ok(SurfacePartialsRun::Cancelled),
            }
        };
        let dv = {
            if should_cancel() {
                return Ok(SurfacePartialsRun::Cancelled);
            }
            let mut v_net = Vec::new();
            v_net
                .try_reserve_exact(knots_v.control_count())
                .map_err(|_| NurbsError::Domain {
                    what: "surface v-isocurve allocation was refused".to_string(),
                })?;
            let mut operations_since_poll = 0usize;
            for j in 0..knots_v.control_count() {
                let mut acc = [0.0f64; 4];
                for (row_offset, &weight_u) in basis_u.iter().enumerate() {
                    let control = &self.inner.cpw[span_u - degree_u + row_offset][j];
                    for (accumulator, &component) in acc.iter_mut().zip(control) {
                        *accumulator += weight_u * component;
                        if surface_poll_due(&mut operations_since_poll, should_cancel) {
                            return Ok(SurfacePartialsRun::Cancelled);
                        }
                    }
                }
                for &component in &acc {
                    if !component.is_finite() {
                        return Err(NurbsError::Domain {
                            what: "surface v-isocurve left the admissible homogeneous domain"
                                .to_string(),
                        });
                    }
                    if surface_poll_due(&mut operations_since_poll, should_cancel) {
                        return Ok(SurfacePartialsRun::Cancelled);
                    }
                }
                if !acc[3].is_admissible_weight() {
                    return Err(NurbsError::Domain {
                        what: "surface v-isocurve left the admissible homogeneous domain"
                            .to_string(),
                    });
                }
                if surface_poll_due(&mut operations_since_poll, should_cancel) {
                    return Ok(SurfacePartialsRun::Cancelled);
                }
                v_net.push(acc);
            }
            match NurbsCurve::<f64, 3>::derivatives_from_admitted_parts_after_preflight_with_poll(
                knots_v,
                &v_net,
                v,
                1,
                should_cancel,
            )? {
                CurveDerivativesRun::Complete { derivatives } => derivatives,
                CurveDerivativesRun::Cancelled => return Ok(SurfacePartialsRun::Cancelled),
            }
        };
        if should_cancel() {
            return Ok(SurfacePartialsRun::Cancelled);
        }
        Ok(SurfacePartialsRun::Complete {
            partials: (du[0], du[1], dv[1]),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Rat;
    use asupersync::types::Budget;
    use fs_exec::{CancelGate, ExecMode, StreamKey};

    fn with_surface_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
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
                    seed: 0x5A2F_AC00,
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

    fn bilinear_surface() -> NurbsSurface<f64> {
        let knots_u = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("u knots");
        let knots_v = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("v knots");
        let points = vec![
            vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
        ];
        let weights = vec![vec![1.0; 2]; 2];
        NurbsSurface::new(knots_u, knots_v, &points, &weights).expect("bilinear plane")
    }

    fn high_degree_surface() -> NurbsSurface<f64> {
        let degree_u = 16usize;
        let mut knots_u = vec![0.0; degree_u + 1];
        knots_u.extend(vec![1.0; degree_u + 1]);
        let knots_v = vec![0.0, 0.0, 1.0, 1.0];
        let points: Vec<Vec<[f64; 3]>> = (0..=degree_u)
            .map(|row| {
                let x = row as f64 / degree_u as f64;
                vec![[x, 0.0, 0.0], [x, 1.0, 0.0]]
            })
            .collect();
        let weights = vec![vec![1.0; 2]; degree_u + 1];
        NurbsSurface::new(
            KnotVector::new(knots_u, degree_u).expect("high-degree u knots"),
            KnotVector::new(knots_v, 1).expect("linear v knots"),
            &points,
            &weights,
        )
        .expect("high-degree plane")
    }

    fn long_linear_surface() -> NurbsSurface<f64> {
        let control_count = 130usize;
        let make_knots = || {
            let mut entries = Vec::with_capacity(control_count + 2);
            entries.extend([0.0, 0.0]);
            for index in 1..control_count - 1 {
                entries.push(index as f64 / (control_count - 1) as f64);
            }
            entries.extend([1.0, 1.0]);
            KnotVector::new(entries, 1).expect("long linear knots")
        };
        let points: Vec<Vec<[f64; 3]>> = (0..control_count)
            .map(|row| {
                (0..control_count)
                    .map(|column| {
                        [
                            row as f64 / (control_count - 1) as f64,
                            column as f64 / (control_count - 1) as f64,
                            0.0,
                        ]
                    })
                    .collect()
            })
            .collect();
        let weights = vec![vec![1.0; control_count]; control_count];
        NurbsSurface::new(make_knots(), make_knots(), &points, &weights)
            .expect("long bilinear-grid surface")
    }

    fn asymmetric_surface() -> NurbsSurface<f64> {
        let knots_u =
            KnotVector::new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0], 2).expect("quadratic u knots");
        let knots_v = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("linear v knots");
        let points = vec![
            vec![[0.0, 0.0, 0.0], [0.0, 2.0, 1.0]],
            vec![[1.0, 0.0, 1.0], [1.0, 2.0, 0.0]],
            vec![[3.0, 0.0, 0.0], [3.0, 2.0, 2.0]],
        ];
        let weights = vec![vec![1.0, 2.0], vec![3.0, 1.0], vec![2.0, 3.0]];
        NurbsSurface::new(knots_u, knots_v, &points, &weights).expect("asymmetric rational surface")
    }

    fn exact_asymmetric_surface() -> NurbsSurface<Rat> {
        let zero = Rat::int(0);
        let one = Rat::int(1);
        let two = Rat::int(2);
        let three = Rat::int(3);
        let knots_u = KnotVector::new(vec![zero, zero, zero, one, one, one], 2)
            .expect("exact quadratic u knots");
        let knots_v = KnotVector::new(vec![zero, zero, one, one], 1).expect("exact linear v knots");
        let points = vec![
            vec![[zero, zero, zero], [zero, two, one]],
            vec![[one, zero, one], [one, two, zero]],
            vec![[three, zero, zero], [three, two, two]],
        ];
        let weights = vec![vec![one, two], vec![three, one], vec![two, three]];
        NurbsSurface::new(knots_u, knots_v, &points, &weights)
            .expect("exact asymmetric rational surface")
    }

    fn curve_oracle_insert_u<S: Scalar>(surface: &NurbsSurface<S>, t: S) -> NurbsSurface<S> {
        let control_count_v = surface.knots_v.control_count();
        let mut output: Option<Vec<Vec<[S; 4]>>> = None;
        let mut inserted_knots = None;
        for column in 0..control_count_v {
            let controls: Vec<[S; 4]> = surface.cpw.iter().map(|row| row[column]).collect();
            let curve = NurbsCurve::<S, 3>::from_homogeneous(
                surface.knots_u.try_clone().expect("oracle u knots"),
                controls,
            )
            .expect("oracle u isocurve");
            let refined = curve.insert_knot(t).expect("oracle u insertion");
            if output.is_none() {
                output = Some(vec![Vec::new(); refined.cpw.len()]);
            }
            for (row, &control) in output
                .as_mut()
                .expect("oracle rows initialized")
                .iter_mut()
                .zip(&refined.cpw)
            {
                row.push(control);
            }
            inserted_knots = Some(refined.knots);
        }
        NurbsSurface::from_homogeneous(
            inserted_knots.expect("surface has v controls"),
            surface
                .knots_v
                .try_clone()
                .expect("oracle unchanged v knots"),
            output.expect("surface has v controls"),
        )
        .expect("oracle u surface")
    }

    fn curve_oracle_insert_v<S: Scalar>(surface: &NurbsSurface<S>, t: S) -> NurbsSurface<S> {
        let mut output = Vec::new();
        let mut inserted_knots = None;
        for source in &surface.cpw {
            let curve = NurbsCurve::<S, 3>::from_homogeneous(
                surface.knots_v.try_clone().expect("oracle v knots"),
                source.to_vec(),
            )
            .expect("oracle v isocurve");
            let refined = curve.insert_knot(t).expect("oracle v insertion");
            inserted_knots = Some(refined.knots);
            output.push(refined.cpw);
        }
        NurbsSurface::from_homogeneous(
            surface
                .knots_u
                .try_clone()
                .expect("oracle unchanged u knots"),
            inserted_knots.expect("surface has u controls"),
            output,
        )
        .expect("oracle v surface")
    }

    #[test]
    fn surface_admission_with_cx_is_transactional_and_lifetime_bound() {
        let surface = bilinear_surface();
        with_surface_cx(true, |cx| {
            assert!(matches!(
                surface.admit_with_cx(cx).expect("valid source"),
                SurfaceAdmissionRun::Cancelled
            ));

            let mut invalid_u = bilinear_surface();
            invalid_u.knots_u.knots.clear();
            assert!(matches!(
                invalid_u.admit_with_cx(cx),
                Err(NurbsError::Structure { .. })
            ));
        });
        with_surface_cx(false, |cx| {
            let SurfaceAdmissionRun::Complete { admitted } = surface
                .admit_with_cx(cx)
                .expect("healthy cancellable admission")
            else {
                panic!("active context must admit the valid surface");
            };
            assert!(core::ptr::eq(admitted.source(), &surface));
            assert!(matches!(
                admitted
                    .eval_with_cx(0.5, 0.5, cx)
                    .expect("admitted cancellable evaluation"),
                SurfaceEvaluationRun::Complete { .. }
            ));
        });

        let mut malformed = bilinear_surface();
        malformed.cpw.clear();
        let legacy_error = malformed.admit().expect_err("malformed legacy admission");
        with_surface_cx(false, |cx| {
            assert_eq!(
                malformed
                    .admit_with_cx(cx)
                    .expect_err("malformed cancellable admission"),
                legacy_error
            );
        });
    }

    #[test]
    fn surface_admission_replays_inside_controls_and_at_publication() {
        let surface = high_degree_surface();
        let run = || {
            let mut polls = 0usize;
            let mut should_cancel = || {
                polls += 1;
                polls == 12
            };
            let outcome = surface
                .validate_live_structure_with_poll(&mut should_cancel)
                .expect("valid high-degree surface");
            (
                matches!(outcome, SurfaceValidationOutcome::Cancelled),
                polls,
            )
        };
        assert_eq!(run(), run());
        assert_eq!(run(), (true, 12));

        let bilinear = bilinear_surface();
        let mut total_polls = 0usize;
        let mut never_cancel = || {
            total_polls += 1;
            false
        };
        assert!(matches!(
            bilinear
                .validate_live_structure_with_poll(&mut never_cancel)
                .expect("healthy admission"),
            SurfaceValidationOutcome::Complete
        ));
        assert_eq!(total_polls, 12);
        let mut replay_polls = 0usize;
        let mut cancel_at_publication = || {
            replay_polls += 1;
            replay_polls == 12
        };
        assert!(matches!(
            bilinear
                .validate_live_structure_with_poll(&mut cancel_at_publication)
                .expect("publication cancellation"),
            SurfaceValidationOutcome::Cancelled
        ));
        assert_eq!(replay_polls, 12);
    }

    #[test]
    fn admitted_directional_insertion_with_cx_is_transactional_and_exact() {
        let surface = asymmetric_surface();
        let admitted = surface.admit().expect("admitted asymmetric surface");
        let inserted_u = admitted.insert_knot_u(0.5).expect("legacy u insertion");
        let inserted_v = admitted.insert_knot_v(0.5).expect("legacy v insertion");
        assert_eq!(inserted_u, curve_oracle_insert_u(&surface, 0.5));
        assert_eq!(inserted_v, curve_oracle_insert_v(&surface, 0.5));
        assert_eq!(inserted_u.knots_v, surface.knots_v);
        assert_eq!(inserted_v.knots_u, surface.knots_u);
        assert_eq!(inserted_u.cpw.len(), 4);
        assert!(inserted_u.cpw.iter().all(|row| row.len() == 2));
        assert_eq!(inserted_v.cpw.len(), 3);
        assert!(inserted_v.cpw.iter().all(|row| row.len() == 3));

        with_surface_cx(false, |cx| {
            assert_eq!(
                admitted
                    .insert_knot_u_with_cx(0.5, cx)
                    .expect("active u insertion"),
                SurfaceInsertionRun::Complete {
                    surface: inserted_u,
                }
            );
            assert_eq!(
                admitted
                    .insert_knot_v_with_cx(0.5, cx)
                    .expect("active v insertion"),
                SurfaceInsertionRun::Complete {
                    surface: inserted_v,
                }
            );
        });

        let exact_surface = exact_asymmetric_surface();
        let exact_admitted = exact_surface.admit().expect("admitted exact surface");
        let half = Rat::new(1, 2);
        let exact_u = exact_admitted
            .insert_knot_u(half)
            .expect("legacy exact u insertion");
        let exact_v = exact_admitted
            .insert_knot_v(half)
            .expect("legacy exact v insertion");
        assert_eq!(exact_u, curve_oracle_insert_u(&exact_surface, half));
        assert_eq!(exact_v, curve_oracle_insert_v(&exact_surface, half));
        for (u, v) in [
            (Rat::new(1, 4), Rat::new(1, 3)),
            (Rat::new(3, 4), Rat::new(2, 3)),
        ] {
            let source_point = exact_admitted.eval(u, v).expect("exact source evaluation");
            assert_eq!(
                exact_u.eval(u, v).expect("exact u-refined evaluation"),
                source_point
            );
            assert_eq!(
                exact_v.eval(u, v).expect("exact v-refined evaluation"),
                source_point
            );
        }
        with_surface_cx(false, |cx| {
            assert_eq!(
                exact_admitted
                    .insert_knot_u_with_cx(half, cx)
                    .expect("active exact u insertion"),
                SurfaceInsertionRun::Complete { surface: exact_u }
            );
            assert_eq!(
                exact_admitted
                    .insert_knot_v_with_cx(half, cx)
                    .expect("active exact v insertion"),
                SurfaceInsertionRun::Complete { surface: exact_v }
            );
        });
    }

    #[test]
    fn directional_insertion_refusals_precede_cancellation() {
        let surface = bilinear_surface();
        let admitted = surface.admit().expect("admitted plane");
        let u_endpoint = admitted
            .insert_knot_u(0.0)
            .expect_err("legacy u endpoint refusal");
        let v_endpoint = admitted
            .insert_knot_v(1.0)
            .expect_err("legacy v endpoint refusal");
        let u_non_finite = admitted
            .insert_knot_u(f64::NAN)
            .expect_err("legacy u non-finite refusal");
        let v_non_finite = admitted
            .insert_knot_v(f64::NAN)
            .expect_err("legacy v non-finite refusal");
        with_surface_cx(true, |cx| {
            assert_eq!(
                admitted
                    .insert_knot_u_with_cx(0.5, cx)
                    .expect("valid u request reaches cancellation"),
                SurfaceInsertionRun::Cancelled
            );
            assert_eq!(
                admitted
                    .insert_knot_v_with_cx(0.5, cx)
                    .expect("valid v request reaches cancellation"),
                SurfaceInsertionRun::Cancelled
            );
            assert_eq!(
                admitted
                    .insert_knot_u_with_cx(0.0, cx)
                    .expect_err("u endpoint refusal beats cancellation"),
                u_endpoint
            );
            assert_eq!(
                admitted
                    .insert_knot_v_with_cx(1.0, cx)
                    .expect_err("v endpoint refusal beats cancellation"),
                v_endpoint
            );
            assert_eq!(
                admitted
                    .insert_knot_u_with_cx(f64::NAN, cx)
                    .expect_err("u non-finite refusal beats cancellation"),
                u_non_finite
            );
            assert_eq!(
                admitted
                    .insert_knot_v_with_cx(f64::NAN, cx)
                    .expect_err("v non-finite refusal beats cancellation"),
                v_non_finite
            );
        });
    }

    #[test]
    fn directional_insertion_cancels_inside_assembly_and_at_publication() {
        let high_degree = high_degree_surface();
        let plan = high_degree
            .insertion_plan_after_parameter(0.5, SurfaceInsertionAxis::U)
            .expect("high-degree u plan");
        let span = high_degree
            .knots_u
            .admitted_after_validation()
            .span(0.5)
            .expect("high-degree u span");
        let run = || {
            let mut polls = 0usize;
            let mut should_cancel = || {
                polls += 1;
                // Entry + one checkpoint per output-row reservation precede
                // the first 64-operation assembly-stride checkpoint.
                polls == plan.new_control_count_u + 2
            };
            let outcome = high_degree
                .insert_knot_at_span_with_plan_and_poll(0.5, span, plan, &mut should_cancel)
                .expect("bounded high-degree insertion");
            (matches!(outcome, SurfaceInsertionRun::Cancelled), polls)
        };
        assert_eq!(run(), run());
        assert_eq!(run(), (true, plan.new_control_count_u + 2));

        let surface = bilinear_surface();
        for axis in [SurfaceInsertionAxis::U, SurfaceInsertionAxis::V] {
            let plan = surface
                .insertion_plan_after_parameter(0.5, axis)
                .expect("bilinear insertion plan");
            let span = match axis {
                SurfaceInsertionAxis::U => surface
                    .knots_u
                    .admitted_after_validation()
                    .span(0.5)
                    .expect("u span"),
                SurfaceInsertionAxis::V => surface
                    .knots_v
                    .admitted_after_validation()
                    .span(0.5)
                    .expect("v span"),
            };
            let mut total_polls = 0usize;
            let mut never_cancel = || {
                total_polls += 1;
                false
            };
            let complete = surface
                .insert_knot_at_span_with_plan_and_poll(0.5, span, plan, &mut never_cancel)
                .expect("healthy insertion");
            assert!(matches!(complete, SurfaceInsertionRun::Complete { .. }));

            let mut replay_polls = 0usize;
            let mut cancel_at_publication = || {
                replay_polls += 1;
                replay_polls == total_polls
            };
            let cancelled = surface
                .insert_knot_at_span_with_plan_and_poll(0.5, span, plan, &mut cancel_at_publication)
                .expect("publication cancellation");
            assert!(matches!(cancelled, SurfaceInsertionRun::Cancelled));
            assert_eq!(replay_polls, total_polls);
        }
    }

    #[test]
    fn admitted_surface_evaluation_with_cx_is_transactional_and_exact() {
        let surface = bilinear_surface();
        let admitted = surface.admit().expect("admitted plane");
        with_surface_cx(true, |cx| {
            assert!(matches!(
                admitted
                    .eval_with_cx(0.25, 0.75, cx)
                    .expect("valid request"),
                SurfaceEvaluationRun::Cancelled
            ));
            assert!(matches!(
                admitted.eval_with_cx(-1.0, 2.0, cx),
                Err(NurbsError::Domain { .. })
            ));
        });
        with_surface_cx(false, |cx| {
            assert_eq!(
                admitted
                    .eval_with_cx(0.25, 0.75, cx)
                    .expect("active context"),
                SurfaceEvaluationRun::Complete {
                    point: admitted.eval(0.25, 0.75).expect("legacy evaluation"),
                }
            );
        });

        let degree = 32usize;
        let mut high_overlap_knots = vec![0.0; degree + 1];
        high_overlap_knots.extend(vec![0.5; degree]);
        high_overlap_knots.extend(vec![1.0; degree + 1]);
        let control_count = 2 * degree + 1;
        let high_overlap = NurbsSurface::new(
            KnotVector::new(high_overlap_knots.clone(), degree).expect("high-overlap U knots"),
            KnotVector::new(high_overlap_knots, degree).expect("high-overlap V knots"),
            &vec![vec![[0.0; 3]; control_count]; control_count],
            &vec![vec![1.0; control_count]; control_count],
        )
        .expect("high-overlap surface");
        let admitted_high_overlap = high_overlap.admit().expect("admitted high-overlap surface");
        let work_error = admitted_high_overlap
            .span_boxes()
            .expect_err("legacy span-box work refusal");
        with_surface_cx(true, |cx| {
            assert_eq!(
                admitted_high_overlap
                    .span_boxes_with_cx(cx)
                    .expect_err("work refusal must precede cancellation"),
                work_error
            );
        });

        let zero = Rat::int(0);
        let one = Rat::int(1);
        let knots_u = KnotVector::new(vec![zero, zero, one, one], 1).expect("exact u knots");
        let knots_v = KnotVector::new(vec![zero, zero, one, one], 1).expect("exact v knots");
        let points = vec![
            vec![[zero, zero, zero], [zero, one, zero]],
            vec![[one, zero, zero], [one, one, zero]],
        ];
        let weights = vec![vec![one; 2]; 2];
        let exact_surface =
            NurbsSurface::new(knots_u, knots_v, &points, &weights).expect("exact bilinear plane");
        let admitted_exact = exact_surface.admit().expect("admitted exact plane");
        let quarter = Rat::new(1, 4);
        with_surface_cx(false, |cx| {
            assert_eq!(
                admitted_exact
                    .eval_with_cx(quarter, quarter, cx)
                    .expect("exact active context"),
                SurfaceEvaluationRun::Complete {
                    point: admitted_exact
                        .eval(quarter, quarter)
                        .expect("exact legacy evaluation"),
                }
            );
        });
    }

    #[test]
    fn admitted_surface_partials_with_cx_are_transactional_and_exact() {
        let surface = bilinear_surface();
        let admitted = surface.admit().expect("admitted bilinear plane");
        with_surface_cx(true, |cx| {
            assert!(matches!(
                admitted
                    .partials_with_cx(0.25, 0.75, cx)
                    .expect("valid partial request"),
                SurfacePartialsRun::Cancelled
            ));

            let u_error = admitted
                .partials(-1.0, 2.0)
                .expect_err("legacy U-parameter refusal");
            assert_eq!(
                admitted
                    .partials_with_cx(-1.0, 2.0, cx)
                    .expect_err("U refusal must precede cancellation"),
                u_error
            );
            let v_error = admitted
                .partials(0.5, 2.0)
                .expect_err("legacy V-parameter refusal");
            assert_eq!(
                admitted
                    .partials_with_cx(0.5, 2.0, cx)
                    .expect_err("V refusal must precede cancellation"),
                v_error
            );
        });
        with_surface_cx(false, |cx| {
            assert_eq!(
                admitted
                    .partials_with_cx(0.25, 0.75, cx)
                    .expect("active partial request"),
                SurfacePartialsRun::Complete {
                    partials: admitted
                        .partials(0.25, 0.75)
                        .expect("legacy partial request"),
                }
            );
        });
    }

    #[test]
    fn surface_partial_cancellation_replays_inside_work_and_at_publication() {
        let surface = high_degree_surface();
        let admitted = surface.admit().expect("admitted high-degree surface");
        let (span_u, basis_u) = admitted.knots_u().basis(0.25).expect("U basis");
        let (span_v, basis_v) = admitted.knots_v().basis(0.75).expect("V basis");
        let run = || {
            let mut polls = 0usize;
            let mut should_cancel = || {
                polls += 1;
                polls == 2
            };
            let outcome = admitted
                .partials_from_basis_with_poll(
                    0.25,
                    0.75,
                    span_u,
                    &basis_u,
                    span_v,
                    &basis_v,
                    &mut should_cancel,
                )
                .expect("valid high-degree partial work");
            (matches!(outcome, SurfacePartialsRun::Cancelled), polls)
        };
        assert_eq!(run(), run());
        assert_eq!(run(), (true, 2));

        let bilinear = bilinear_surface();
        let admitted_bilinear = bilinear.admit().expect("admitted bilinear surface");
        let (bilinear_span_u, bilinear_basis_u) =
            admitted_bilinear.knots_u().basis(0.25).expect("U basis");
        let (bilinear_span_v, bilinear_basis_v) =
            admitted_bilinear.knots_v().basis(0.75).expect("V basis");
        let mut total_polls = 0usize;
        let mut never_cancel = || {
            total_polls += 1;
            false
        };
        assert!(matches!(
            admitted_bilinear
                .partials_from_basis_with_poll(
                    0.25,
                    0.75,
                    bilinear_span_u,
                    &bilinear_basis_u,
                    bilinear_span_v,
                    &bilinear_basis_v,
                    &mut never_cancel,
                )
                .expect("healthy partial work"),
            SurfacePartialsRun::Complete { .. }
        ));
        assert_eq!(total_polls, 31);
        let mut replay_polls = 0usize;
        let mut cancel_at_publication = || {
            replay_polls += 1;
            replay_polls == 31
        };
        assert!(matches!(
            admitted_bilinear
                .partials_from_basis_with_poll(
                    0.25,
                    0.75,
                    bilinear_span_u,
                    &bilinear_basis_u,
                    bilinear_span_v,
                    &bilinear_basis_v,
                    &mut cancel_at_publication,
                )
                .expect("publication cancellation"),
            SurfacePartialsRun::Cancelled
        ));
        assert_eq!(replay_polls, 31);
    }

    #[test]
    fn admitted_surface_span_boxes_with_cx_are_transactional_and_exact() {
        let surface = bilinear_surface();
        let admitted = surface.admit().expect("admitted bilinear surface");
        with_surface_cx(true, |cx| {
            assert!(matches!(
                admitted
                    .span_boxes_with_cx(cx)
                    .expect("valid span-box request"),
                SurfaceSpanBoxesRun::Cancelled
            ));
        });
        with_surface_cx(false, |cx| {
            assert_eq!(
                admitted
                    .span_boxes_with_cx(cx)
                    .expect("active span-box request"),
                SurfaceSpanBoxesRun::Complete {
                    boxes: admitted.span_boxes().expect("legacy span boxes"),
                }
            );
        });

        let zero = Rat::int(0);
        let one = Rat::int(1);
        let exact = NurbsSurface::new(
            KnotVector::new(vec![zero, zero, one, one], 1).expect("exact U knots"),
            KnotVector::new(vec![zero, zero, one, one], 1).expect("exact V knots"),
            &vec![
                vec![[zero, zero, zero], [zero, one, zero]],
                vec![[one, zero, zero], [one, one, zero]],
            ],
            &vec![vec![one; 2]; 2],
        )
        .expect("exact bilinear surface");
        let admitted_exact = exact.admit().expect("admitted exact surface");
        with_surface_cx(false, |cx| {
            assert_eq!(
                admitted_exact
                    .span_boxes_with_cx(cx)
                    .expect("active exact span-box request"),
                SurfaceSpanBoxesRun::Complete {
                    boxes: admitted_exact
                        .span_boxes()
                        .expect("legacy exact span boxes"),
                }
            );
        });
    }

    #[test]
    fn surface_span_box_cancellation_replays_inside_work_and_at_publication() {
        let surface = high_degree_surface();
        let admitted = surface.admit().expect("admitted high-degree surface");
        let run = || {
            let mut polls = 0usize;
            let mut should_cancel = || {
                polls += 1;
                polls == 2
            };
            let outcome = admitted
                .span_boxes_with_poll(&mut should_cancel)
                .expect("valid high-degree span-box work");
            (matches!(outcome, SurfaceSpanBoxesRun::Cancelled), polls)
        };
        assert_eq!(run(), run());
        assert_eq!(run(), (true, 2));

        let bilinear = bilinear_surface();
        let admitted_bilinear = bilinear.admit().expect("admitted bilinear surface");
        let mut total_polls = 0usize;
        let mut never_cancel = || {
            total_polls += 1;
            false
        };
        assert!(matches!(
            admitted_bilinear
                .span_boxes_with_poll(&mut never_cancel)
                .expect("healthy span-box work"),
            SurfaceSpanBoxesRun::Complete { .. }
        ));
        assert_eq!(total_polls, 2);
        let mut replay_polls = 0usize;
        let mut cancel_at_publication = || {
            replay_polls += 1;
            replay_polls == 2
        };
        assert!(matches!(
            admitted_bilinear
                .span_boxes_with_poll(&mut cancel_at_publication)
                .expect("publication cancellation"),
            SurfaceSpanBoxesRun::Cancelled
        ));
        assert_eq!(replay_polls, 2);
    }

    #[test]
    fn surface_evaluation_replays_inside_tensor_work_and_at_publication() {
        let surface = high_degree_surface();
        let admitted = surface.admit().expect("admitted high-degree plane");
        let (span_u, basis_u) = admitted.knots_u().basis(0.5).expect("u basis");
        let (span_v, basis_v) = admitted.knots_v().basis(0.5).expect("v basis");
        let run = || {
            let mut polls = 0usize;
            let outcome = admitted
                .eval_from_basis_with_poll(span_u, &basis_u, span_v, &basis_v, || {
                    polls += 1;
                    polls == 2
                })
                .expect("finite tensor contraction");
            (outcome, polls)
        };
        assert_eq!(run(), run());
        assert_eq!(run(), (SurfaceEvaluationRun::Cancelled, 2));

        let bilinear = bilinear_surface();
        let admitted_bilinear = bilinear.admit().expect("admitted bilinear plane");
        let (bilinear_span_u, bilinear_basis_u) =
            admitted_bilinear.knots_u().basis(0.5).expect("u basis");
        let (bilinear_span_v, bilinear_basis_v) =
            admitted_bilinear.knots_v().basis(0.5).expect("v basis");
        let mut total_polls = 0usize;
        assert!(matches!(
            admitted_bilinear
                .eval_from_basis_with_poll(
                    bilinear_span_u,
                    &bilinear_basis_u,
                    bilinear_span_v,
                    &bilinear_basis_v,
                    || {
                        total_polls += 1;
                        false
                    },
                )
                .expect("healthy bilinear evaluation"),
            SurfaceEvaluationRun::Complete { .. }
        ));
        assert_eq!(total_polls, 2);
        let mut replay_polls = 0usize;
        assert_eq!(
            admitted_bilinear
                .eval_from_basis_with_poll(
                    bilinear_span_u,
                    &bilinear_basis_u,
                    bilinear_span_v,
                    &bilinear_basis_v,
                    || {
                        replay_polls += 1;
                        replay_polls == 2
                    },
                )
                .expect("publication cancellation"),
            SurfaceEvaluationRun::Cancelled
        );
        assert_eq!(replay_polls, 2);
    }

    #[test]
    fn directional_insertion_preflight_prices_work_and_complete_output() {
        let surface = asymmetric_surface();
        let plan_u = surface
            .insertion_plan_after_parameter(0.5, SurfaceInsertionAxis::U)
            .expect("u insertion plan");
        let plan_v = surface
            .insertion_plan_after_parameter(0.5, SurfaceInsertionAxis::V)
            .expect("v insertion plan");
        assert_eq!(
            (
                plan_u.new_control_count_u,
                plan_u.new_control_count_v,
                plan_u.new_knot_count_u,
                plan_u.new_knot_count_v,
            ),
            (4, 2, 7, 4)
        );
        assert_eq!(
            (
                plan_v.new_control_count_u,
                plan_v.new_control_count_v,
                plan_v.new_knot_count_u,
                plan_v.new_knot_count_v,
            ),
            (3, 3, 6, 5)
        );
        assert!(plan_u.work_units > 0 && plan_u.retained_bytes > 0);
        assert!(plan_v.work_units > 0 && plan_v.retained_bytes > 0);
        assert!(
            enforce_surface_insertion_envelope(
                BASIS_MAX_WORK_UNITS,
                SURFACE_INSERTION_MAX_RETAINED_BYTES,
            )
            .is_ok()
        );
        assert!(matches!(
            enforce_surface_insertion_envelope(
                BASIS_MAX_WORK_UNITS + 1,
                SURFACE_INSERTION_MAX_RETAINED_BYTES,
            ),
            Err(NurbsError::Domain { ref what }) if what.contains("work")
        ));
        assert!(matches!(
            enforce_surface_insertion_envelope(
                BASIS_MAX_WORK_UNITS,
                SURFACE_INSERTION_MAX_RETAINED_BYTES + 1,
            ),
            Err(NurbsError::Domain { ref what }) if what.contains("retain")
        ));
    }

    #[test]
    fn span_box_preflight_prices_nested_scans_and_retained_output() {
        let bytes_per_box = core::mem::size_of::<SurfaceSpanBox<f64>>();
        assert_eq!(
            preflight_span_boxes(2, 2, 1, 1, bytes_per_box).expect("one bilinear box"),
            1
        );
        assert!(
            enforce_span_box_envelope(BASIS_MAX_WORK_UNITS, SURFACE_SPAN_BOX_MAX_RETAINED_BYTES)
                .is_ok(),
            "both exact ceilings are admitted"
        );
        assert!(matches!(
            enforce_span_box_envelope(
                BASIS_MAX_WORK_UNITS + 1,
                SURFACE_SPAN_BOX_MAX_RETAINED_BYTES
            ),
            Err(NurbsError::Domain { ref what }) if what.contains("work")
        ));
        assert!(matches!(
            enforce_span_box_envelope(
                BASIS_MAX_WORK_UNITS,
                SURFACE_SPAN_BOX_MAX_RETAINED_BYTES + 1
            ),
            Err(NurbsError::Domain { ref what }) if what.contains("retain")
        ));

        let work_error = preflight_span_boxes(512, 512, 255, 255, bytes_per_box)
            .expect_err("high-degree overlap must be refused before allocation");
        assert!(matches!(
            work_error,
            NurbsError::Domain { ref what } if what.contains("traversal")
        ));

        let rat_box_bytes = core::mem::size_of::<SurfaceSpanBox<Rat>>();
        preflight_span_boxes(458, 458, 1, 1, rat_box_bytes)
            .expect("the Rat payload immediately below the retained-byte cap is admitted");
        let retained_error = preflight_span_boxes(459, 459, 1, 1, rat_box_bytes)
            .expect_err("the next Rat span grid exceeds retained bytes before allocation");
        assert!(matches!(
            retained_error,
            NurbsError::Domain { ref what } if what.contains("retain")
        ));
    }

    #[test]
    fn partials_preflight_prices_union_and_preserves_exact_cap_boundaries() {
        preflight_partials_envelope(2, 2, 4, 4, 1, 1)
            .expect("bilinear partial request is admitted");
        let work_error = preflight_partials_envelope(1_001, 1_001, 2_002, 2_002, 1_000, 1_000)
            .expect_err("high-order tensor work must refuse through the aggregate helper");
        assert!(matches!(
            work_error,
            NurbsError::Domain { ref what } if what.contains("work")
        ));
        let retained_error = preflight_partials_envelope(600_000, 2, 600_002, 4, 1, 1)
            .expect_err("large asymmetric isocurve retention must refuse through the helper");
        assert!(matches!(
            retained_error,
            NurbsError::Domain { ref what } if what.contains("retain")
        ));
        assert_eq!(
            NurbsCurve::<f64, 3>::derivative_envelope(2, 4, 1, 1)
                .expect("linear derivative envelope"),
            (44, 304 + 2 * core::mem::size_of::<Vec<[f64; 4]>>() as u128,)
        );
        assert!(
            enforce_partials_envelope(
                NurbsCurve::<f64, 3>::MAX_DERIVATIVE_WORK_UNITS,
                NurbsCurve::<f64, 3>::MAX_DERIVATIVE_RETAINED_BYTES,
            )
            .is_ok()
        );
        assert!(matches!(
            enforce_partials_envelope(
                NurbsCurve::<f64, 3>::MAX_DERIVATIVE_WORK_UNITS + 1,
                NurbsCurve::<f64, 3>::MAX_DERIVATIVE_RETAINED_BYTES,
            ),
            Err(NurbsError::Domain { ref what }) if what.contains("work")
        ));
        assert!(matches!(
            enforce_partials_envelope(
                NurbsCurve::<f64, 3>::MAX_DERIVATIVE_WORK_UNITS,
                NurbsCurve::<f64, 3>::MAX_DERIVATIVE_RETAINED_BYTES + 1,
            ),
            Err(NurbsError::Domain { ref what }) if what.contains("retain")
        ));
    }

    #[test]
    fn admitted_partials_refuse_parameters_u_first_and_accept_endpoints() {
        let line_u = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("u knots");
        let line_v = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("v knots");
        let points = vec![
            vec![[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            vec![[1.0, 0.0, 0.0], [1.0, 1.0, 0.0]],
        ];
        let weights = vec![vec![1.0; 2]; 2];
        let surface = NurbsSurface::new(line_u, line_v, &points, &weights).expect("bilinear plane");
        let admitted = surface.admit().expect("admitted plane");

        let u_error = admitted
            .partials(-1.0, 2.0)
            .expect_err("U must win when both parameters are invalid");
        assert!(matches!(
            u_error,
            NurbsError::Domain { ref what } if what.contains("u-partial")
        ));
        let v_error = admitted
            .partials(0.5, 2.0)
            .expect_err("invalid V must refuse before allocation");
        assert!(matches!(
            v_error,
            NurbsError::Domain { ref what } if what.contains("v-partial")
        ));
        assert!(matches!(
            admitted.partials(f64::NAN, 0.5),
            Err(NurbsError::Domain { ref what }) if what.contains("u-partial")
        ));
        assert!(matches!(
            admitted.partials(0.5, f64::NAN),
            Err(NurbsError::Domain { ref what }) if what.contains("v-partial")
        ));
        admitted
            .partials(0.0, 0.0)
            .expect("lower endpoints are admitted");
        admitted
            .partials(1.0, 1.0)
            .expect("upper endpoints are admitted");
    }

    #[test]
    fn admitted_partials_refuse_an_ordinary_derivative_at_a_c0_knot() {
        let knots_u = KnotVector::new(vec![0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0], 2)
            .expect("C0 quadratic u knots");
        let knots_v = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("linear v knots");
        let points: Vec<Vec<[f64; 3]>> = (0..5)
            .map(|u| vec![[f64::from(u), 0.0, 0.0], [f64::from(u), 1.0, 0.0]])
            .collect();
        let weights = vec![vec![1.0; 2]; 5];
        let surface = NurbsSurface::new(knots_u, knots_v, &points, &weights).expect("C0 surface");
        let admitted = surface.admit().expect("admitted C0 surface");
        admitted
            .partials(0.25, 0.5)
            .expect("ordinary partial inside a smooth span");
        let continuity_error = admitted
            .partials(0.5, 0.5)
            .expect_err("ordinary derivative at the C0 knot must refuse");
        assert!(matches!(
            &continuity_error,
            NurbsError::Domain { what } if what.contains("multiplicity 2")
        ));
        with_surface_cx(true, |cx| {
            assert_eq!(
                admitted
                    .partials_with_cx(0.5, 0.5, cx)
                    .expect_err("continuity refusal must precede cancellation"),
                continuity_error
            );
        });
    }

    #[test]
    fn partial_parameter_refusal_precedes_surface_structure_scan() {
        let line_u = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("u knots");
        let line_v = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("v knots");
        let malformed = NurbsSurface::<f64> {
            knots_u: line_u,
            knots_v: line_v,
            cpw: Vec::new(),
        };
        let u_error = malformed
            .partials(-1.0, 2.0)
            .expect_err("u is refused first when both parameters are invalid");
        assert!(matches!(
            u_error,
            NurbsError::Domain { ref what } if what.contains("u-partial")
        ));
        let v_error = malformed
            .partials(0.5, 2.0)
            .expect_err("v is refused before malformed controls are scanned");
        assert!(matches!(
            v_error,
            NurbsError::Domain { ref what } if what.contains("v-partial")
        ));
    }

    #[test]
    fn surface_copy_is_fallible_and_late_rows_are_validated_before_output() {
        let line_u = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("u knots");
        let line_v = KnotVector::new(vec![0.0, 0.0, 1.0, 1.0], 1).expect("v knots");
        let points = vec![vec![[0.0; 3]; 2], vec![[1.0; 3]]];
        let weights = vec![vec![1.0; 2], vec![1.0]];
        assert!(matches!(
            NurbsSurface::new(
                line_u.try_clone().expect("fallible u-knot copy"),
                line_v.try_clone().expect("fallible v-knot copy"),
                &points,
                &weights,
            ),
            Err(NurbsError::Structure { ref what }) if what.contains("control columns")
        ));

        let valid_points = vec![vec![[0.0; 3]; 2], vec![[1.0; 3]; 2]];
        let valid_weights = vec![vec![1.0; 2]; 2];
        let surface = NurbsSurface::new(line_u, line_v, &valid_points, &valid_weights)
            .expect("bilinear surface");
        assert_eq!(surface.try_clone().expect("fallible surface copy"), surface);
    }

    #[test]
    fn surface_copy_with_cx_is_transactional_and_exact() {
        let surface = asymmetric_surface();
        with_surface_cx(true, |cx| {
            assert_eq!(
                surface
                    .try_clone_with_cx(cx)
                    .expect("admitted copy request"),
                SurfaceCloneRun::Cancelled
            );
        });
        with_surface_cx(false, |cx| {
            assert_eq!(
                surface.try_clone_with_cx(cx).expect("active surface copy"),
                SurfaceCloneRun::Complete {
                    surface: surface.try_clone().expect("legacy surface copy"),
                }
            );
        });

        let exact = exact_asymmetric_surface();
        with_surface_cx(false, |cx| {
            assert_eq!(
                exact.try_clone_with_cx(cx).expect("active exact copy"),
                SurfaceCloneRun::Complete {
                    surface: exact.try_clone().expect("legacy exact copy"),
                }
            );
        });
    }

    #[test]
    fn surface_copy_cancels_inside_each_nested_copy_and_at_publication() {
        let surface = long_linear_surface();
        let admitted = surface.admitted_after_validation();
        let run = |target| {
            let mut polls = 0usize;
            let mut should_cancel = || {
                polls += 1;
                polls == target
            };
            let outcome = admitted
                .try_clone_with_poll(&mut should_cancel)
                .expect("bounded surface copy");
            (outcome.is_none(), polls)
        };
        assert_eq!(run(3), run(3));
        assert_eq!(run(3), (true, 3));
        assert_eq!(run(7), run(7));
        assert_eq!(run(7), (true, 7));
        assert_eq!(run(13), run(13));
        assert_eq!(run(13), (true, 13));

        let surface = bilinear_surface();
        let admitted = surface.admitted_after_validation();
        let mut total_polls = 0usize;
        let mut never_cancel = || {
            total_polls += 1;
            false
        };
        assert!(
            admitted
                .try_clone_with_poll(&mut never_cancel)
                .expect("healthy surface copy")
                .is_some()
        );
        let mut replay_polls = 0usize;
        let mut cancel_at_publication = || {
            replay_polls += 1;
            replay_polls == total_polls
        };
        assert!(
            admitted
                .try_clone_with_poll(&mut cancel_at_publication)
                .expect("publication cancellation")
                .is_none()
        );
        assert_eq!(replay_polls, total_polls);
    }

    #[test]
    fn surface_copy_preflight_refuses_work_before_nested_retained_bytes() {
        let excessive_columns = (BASIS_MAX_WORK_UNITS / 4 + 1) as usize;
        let work_error = preflight_surface_copy::<f64>(0, 0, 1, excessive_columns)
            .expect_err("work must refuse before retained-byte accounting");
        assert!(matches!(
            work_error,
            NurbsError::Domain { ref what } if what.contains("work units above defensive ceiling")
        ));

        let row_bytes = core::mem::size_of::<Vec<[f64; 4]>>() as u128;
        let control_bytes = core::mem::size_of::<[f64; 4]>() as u128;
        let row_count =
            (SURFACE_COPY_MAX_RETAINED_BYTES / (row_bytes + control_bytes) + 1) as usize;
        assert!(
            (row_count as u128) * control_bytes <= SURFACE_COPY_MAX_RETAINED_BYTES,
            "control payload alone stays within the envelope"
        );
        let memory_error = preflight_surface_copy::<f64>(0, 0, row_count, 1)
            .expect_err("nested row headers must count toward retained output");
        assert!(matches!(
            memory_error,
            NurbsError::Domain { ref what } if what.contains("retains")
        ));
    }
}
