//! Trimmed patches with CERTIFIED point classification. Trim loops are
//! held in EXACT RATIONAL form (2-D parameter-space NURBS over `Rat`) —
//! the dual representation the bead demands. Classification is proved,
//! not sampled: if the query point lies strictly outside every Bézier
//! span's control hull box, the curve and its control polygon are
//! homotopic in a region avoiding the point, so the EXACTLY-computed
//! control-polygon winding number IS the curve's winding number.
//! Ambiguous points (inside a hull box after bounded exact subdivision)
//! are honestly `Boundary`, never a guessed in/out.

use crate::NurbsError;
use crate::curve::{
    AdmittedNurbsCurve, BezierConversionPlan, CurveAdmissionRun, CurveEvaluationRun, NurbsCurve,
    SpanBox,
};
use crate::rat::Rat;
use fs_exec::Cx;

/// Defensive work ceiling for one exact trim classification across all loops.
/// This legacy cap bounds public allocation-bearing subdivision even when a
/// caller supplies `max_subdivision = u32::MAX`; explicit caller budgets belong
/// to the successor API.
pub(crate) const TRIM_CLASSIFY_MAX_WORK_UNITS: u128 = 1_048_576;

/// Aggregate retained-memory ceiling for the conversion, span-box, and
/// offending-interval phases of one exact trim classification.
const TRIM_CLASSIFY_MAX_RETAINED_BYTES: u128 = 64 * 1024 * 1024;
const TRIM_SPAN_BOX_WORK_PER_CONTROL: u128 = 16;
const TRIM_WINDING_WORK_PER_CONTROL: u128 = 128;
const TRIM_EXACT_MIDPOINT_WORK_UNITS: u128 = 1_024;
const TRIM_CANCELLATION_STRIDE: usize = 64;

fn trim_poll_due(
    operations_since_poll: &mut usize,
    should_cancel: &mut impl FnMut() -> bool,
) -> bool {
    *operations_since_poll += 1;
    if *operations_since_poll < TRIM_CANCELLATION_STRIDE {
        return false;
    }
    *operations_since_poll = 0;
    should_cancel()
}

/// Minimum charge for admitting one sealed loop before inspecting its
/// knot/control metadata. This makes a huge collection of individually tiny
/// loops reject in O(1), rather than spending unbounded time merely discovering
/// that the aggregate validation exceeds the legacy synchronous envelope.
const TRIM_MIN_LOOP_VALIDATION_WORK_UNITS: u128 = 64;

/// One closed trim loop: an exact rational curve in (u, v) parameter
/// space (closure is validated).
///
/// The exact curve is read-only after construction; callers use
/// [`TrimLoop::curve`] for inspection.
#[derive(Debug, PartialEq)]
pub struct TrimLoop {
    /// The exact 2-D curve.
    pub(crate) curve: NurbsCurve<Rat, 2>,
}

/// A validate-once borrow of one exact immutable trim-loop snapshot.
#[derive(Debug, Clone, Copy)]
pub struct AdmittedTrimLoop<'a> {
    inner: &'a TrimLoop,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TrimLoopValidationOutcome {
    Complete,
    Cancelled,
}

/// Transactional terminal state of cancellation-aware trim-loop admission.
#[must_use]
#[derive(Debug, Clone, Copy)]
pub enum TrimLoopAdmissionRun<'a> {
    /// The exact immutable loop snapshot was fully validated.
    Complete {
        /// Lifetime-bound authority for the validated trim-loop generation.
        admitted: AdmittedTrimLoop<'a>,
    },
    /// Cancellation was observed; no admitted authority was published.
    Cancelled,
}

fn validate_trim_loop_after_endpoints_with_poll(
    curve: AdmittedNurbsCurve<'_, Rat, 2>,
    start: [Rat; 2],
    end: [Rat; 2],
    mut should_cancel: impl FnMut() -> bool,
) -> Result<TrimLoopValidationOutcome, NurbsError> {
    if should_cancel() {
        return Ok(TrimLoopValidationOutcome::Cancelled);
    }
    if start != end {
        return Err(NurbsError::Structure {
            what: "trim loop must close exactly (rational endpoint equality)".to_string(),
        });
    }

    // A full interior knot break carries independent left and right limits.
    // Permit it only when those limits agree exactly in Cartesian space.
    let knots = curve.knots();
    let p = knots.degree();
    let knot_entries = knots.knots();
    let controls = curve.homogeneous_control_points();
    let mut operations_since_poll = 0usize;
    let mut run_start = 0usize;
    while run_start < knot_entries.len() {
        let mut run_end = run_start + 1;
        while run_end < knot_entries.len() && knot_entries[run_end] == knot_entries[run_start] {
            run_end += 1;
            if trim_poll_due(&mut operations_since_poll, &mut should_cancel) {
                return Ok(TrimLoopValidationOutcome::Cancelled);
            }
        }
        let is_interior = run_start != 0 && run_end != knot_entries.len();
        if is_interior && run_end - run_start == p + 1 {
            let left = controls[run_start - 1];
            let right = controls[run_start];
            for coordinate in 0..2 {
                if left[coordinate] * right[3] != right[coordinate] * left[3] {
                    return Err(NurbsError::Structure {
                        what: format!(
                            "trim loop is discontinuous at full knot break {:?}",
                            knot_entries[run_start]
                        ),
                    });
                }
                if trim_poll_due(&mut operations_since_poll, &mut should_cancel) {
                    return Ok(TrimLoopValidationOutcome::Cancelled);
                }
            }
        }
        run_start = run_end;
        if trim_poll_due(&mut operations_since_poll, &mut should_cancel) {
            return Ok(TrimLoopValidationOutcome::Cancelled);
        }
    }
    if should_cancel() {
        return Ok(TrimLoopValidationOutcome::Cancelled);
    }
    Ok(TrimLoopValidationOutcome::Complete)
}

impl TrimLoop {
    fn validate_live(&self) -> Result<(), NurbsError> {
        self.admit().map(|_| ())
    }

    /// Validate closure, continuity, knots, and controls once and bind the
    /// proof to this immutable borrow.
    ///
    /// # Errors
    /// [`NurbsError::Structure`] when the loop is not a valid closed continuous
    /// exact curve.
    pub fn admit(&self) -> Result<AdmittedTrimLoop<'_>, NurbsError> {
        let curve = self.curve.admit()?;
        let (lo, hi) = curve.knots().domain();
        let start = curve.eval(lo)?;
        let end = curve.eval(hi)?;
        match validate_trim_loop_after_endpoints_with_poll(curve, start, end, || false)? {
            TrimLoopValidationOutcome::Complete => Ok(AdmittedTrimLoop { inner: self }),
            TrimLoopValidationOutcome::Cancelled => Err(NurbsError::Domain {
                what: "non-cancelling trim-loop admission observed cancellation".to_string(),
            }),
        }
    }

    /// Validate this exact loop with bounded cancellation polling and publish
    /// only a lifetime-bound admitted view.
    ///
    /// The gate spans curve/knot admission, both exact endpoint evaluations,
    /// full-break continuity traversal, and final authority publication.
    /// Individual exact-rational operations are not preemptible. This method
    /// does not consume the `Cx` budget or finalize its executor scope.
    ///
    /// # Errors
    /// Returns the synchronous admission's work, allocation, structural,
    /// numeric-domain, closure, and continuity refusals when they win before
    /// an observed cancellation.
    pub fn admit_with_cx<'a>(
        &'a self,
        cx: &Cx<'_>,
    ) -> Result<TrimLoopAdmissionRun<'a>, NurbsError> {
        let curve = match self.curve.admit_with_cx(cx)? {
            CurveAdmissionRun::Complete { admitted } => admitted,
            CurveAdmissionRun::Cancelled => return Ok(TrimLoopAdmissionRun::Cancelled),
        };
        let (lo, hi) = curve.knots().domain();
        let start = match curve.eval_with_cx(lo, cx)? {
            CurveEvaluationRun::Complete { point } => point,
            CurveEvaluationRun::Cancelled => return Ok(TrimLoopAdmissionRun::Cancelled),
        };
        let end = match curve.eval_with_cx(hi, cx)? {
            CurveEvaluationRun::Complete { point } => point,
            CurveEvaluationRun::Cancelled => return Ok(TrimLoopAdmissionRun::Cancelled),
        };
        match validate_trim_loop_after_endpoints_with_poll(curve, start, end, || {
            cx.checkpoint().is_err()
        })? {
            TrimLoopValidationOutcome::Complete => Ok(TrimLoopAdmissionRun::Complete {
                admitted: AdmittedTrimLoop { inner: self },
            }),
            TrimLoopValidationOutcome::Cancelled => Ok(TrimLoopAdmissionRun::Cancelled),
        }
    }

    /// Validate closure and construct.
    ///
    /// # Errors
    /// [`NurbsError::Structure`] when the loop is not closed (exact
    /// endpoint equality — this is the rational representation).
    pub fn new(curve: NurbsCurve<Rat, 2>) -> Result<Self, NurbsError> {
        let candidate = TrimLoop { curve };
        candidate.validate_live()?;
        Ok(candidate)
    }

    /// Borrow the sealed exact curve.
    #[must_use]
    pub const fn curve(&self) -> &NurbsCurve<Rat, 2> {
        &self.curve
    }

    /// Fallibly copy this sealed loop without revalidating unchanged data.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] when a destination allocation is refused.
    pub fn try_clone(&self) -> Result<Self, NurbsError> {
        Ok(TrimLoop {
            curve: self.curve.try_clone()?,
        })
    }

    /// The same loop with reversed orientation (holes are wound opposite
    /// to outers under the nonzero rule): control points reversed, knot
    /// vector mirrored about the domain.
    ///
    /// # Errors
    /// [`NurbsError::Structure`] when closure, continuity, knots, or the control
    /// net are invalid.
    pub fn reversed_for_hole(&self) -> Result<TrimLoop, NurbsError> {
        let admitted = self.admit()?;
        let curve = admitted.curve();
        let admitted_knots = curve.knots();
        let (lo, hi) = admitted_knots.domain();
        let mut knots = Vec::new();
        knots
            .try_reserve_exact(admitted_knots.knots().len())
            .map_err(|_| NurbsError::Domain {
                what: "reversed trim-knot allocation was refused".to_string(),
            })?;
        for &knot in admitted_knots.knots().iter().rev() {
            knots.push(lo + (hi - knot));
        }
        let controls = curve.homogeneous_control_points();
        let mut cpw = Vec::new();
        cpw.try_reserve_exact(controls.len())
            .map_err(|_| NurbsError::Domain {
                what: "reversed trim-control allocation was refused".to_string(),
            })?;
        cpw.extend(controls.iter().rev().copied());
        let reversed_knots = crate::basis::KnotVector::new(knots, admitted_knots.degree())?;
        let reversed_curve = NurbsCurve::from_homogeneous(reversed_knots, cpw)?;
        // Reversal of an admitted curve preserves exact endpoint closure and
        // full-break continuity while only changing orientation.
        Ok(TrimLoop {
            curve: reversed_curve,
        })
    }
}

impl<'a> AdmittedTrimLoop<'a> {
    /// The exact immutable source bound to this view.
    #[must_use]
    pub const fn source(&self) -> &'a TrimLoop {
        self.inner
    }

    /// Borrow the admitted exact curve without rescanning it.
    #[must_use]
    pub fn curve(&self) -> crate::curve::AdmittedNurbsCurve<'a, Rat, 2> {
        self.inner.curve.admitted_after_validation()
    }
}

/// A certified classification verdict.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Classification {
    /// Certified inside the trimmed region (nonzero total winding).
    Inside,
    /// Certified outside.
    Outside,
    /// Within the certification band of some trim curve — no in/out
    /// claim is made (the honest verdict on tangent/sliver cases).
    Boundary,
}

/// A trimmed patch: parameter-space loops over any surface. (The surface
/// itself is not needed for classification, which happens in parameter
/// space; carrying it is the B-rep bookkeeping.)
///
/// ```compile_fail
/// use fs_rep_nurbs::TrimmedPatch;
/// let mut patch = TrimmedPatch::new(Vec::new());
/// patch.loops.clear();
/// ```
#[derive(Debug, PartialEq)]
pub struct TrimmedPatch {
    /// Outer boundary + hole loops (orientation encodes solidity via the
    /// nonzero-winding rule: outer CCW, holes CW).
    pub(crate) loops: Vec<TrimLoop>,
    /// Exact-subdivision depth before declaring `Boundary`.
    pub(crate) max_subdivision: u32,
}

/// A validate-once borrow of one exact immutable trimmed-patch snapshot.
#[derive(Debug, Clone, Copy)]
pub struct AdmittedTrimmedPatch<'a> {
    inner: &'a TrimmedPatch,
}

impl TrimmedPatch {
    pub(crate) fn validate_live_with_budget(
        &self,
        work_remaining: &mut u128,
    ) -> Result<(), NurbsError> {
        self.admit_with_budget(work_remaining).map(|_| ())
    }

    fn admit_with_budget<'a>(
        &'a self,
        work_remaining: &mut u128,
    ) -> Result<AdmittedTrimmedPatch<'a>, NurbsError> {
        let minimum_work = (self.loops.len() as u128)
            .checked_mul(TRIM_MIN_LOOP_VALIDATION_WORK_UNITS)
            .ok_or_else(|| NurbsError::Domain {
                what: "trim loop-count validation work overflows u128".to_string(),
            })?;
        if minimum_work > *work_remaining {
            return Err(NurbsError::Domain {
                what: format!(
                    "trim live validation needs at least {minimum_work} work units for {} loops, above the {work_remaining}-unit remaining budget",
                    self.loops.len()
                ),
            });
        }
        let validation_work = self.loops.iter().try_fold(0u128, |total, trim_loop| {
            total
                .checked_add(trim_loop_validation_work(&trim_loop.curve)?)
                .ok_or_else(|| NurbsError::Domain {
                    what: "trim live-validation accounting overflows u128".to_string(),
                })
        })?;
        spend_trim_work(work_remaining, validation_work, "live validation")?;
        for trim_loop in &self.loops {
            trim_loop.admit()?;
        }
        Ok(AdmittedTrimmedPatch { inner: self })
    }

    /// Construct with the default certification depth.
    #[must_use]
    pub fn new(loops: Vec<TrimLoop>) -> Self {
        TrimmedPatch {
            loops,
            max_subdivision: 12,
        }
    }

    /// Construct with an explicit exact-subdivision limit.
    #[must_use]
    pub fn with_max_subdivision(loops: Vec<TrimLoop>, max_subdivision: u32) -> Self {
        TrimmedPatch {
            loops,
            max_subdivision,
        }
    }

    /// Borrow the sealed loop collection.
    #[must_use]
    pub fn loops(&self) -> &[TrimLoop] {
        &self.loops
    }

    /// Exact-subdivision depth before an ambiguous query becomes `Boundary`.
    #[must_use]
    pub const fn max_subdivision(&self) -> u32 {
        self.max_subdivision
    }

    /// Fallibly copy this sealed patch without revalidating unchanged loops.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] when a destination allocation is refused.
    pub fn try_clone(&self) -> Result<Self, NurbsError> {
        let mut loops = Vec::new();
        loops
            .try_reserve_exact(self.loops.len())
            .map_err(|_| NurbsError::Domain {
                what: "trimmed-patch copy loop-table allocation was refused".to_string(),
            })?;
        for trim_loop in &self.loops {
            loops.push(trim_loop.try_clone()?);
        }
        Ok(TrimmedPatch {
            loops,
            max_subdivision: self.max_subdivision,
        })
    }

    /// Validate this exact immutable patch snapshot once under the defensive
    /// aggregate trim budget.
    ///
    /// # Errors
    /// Returns a structured refusal for excessive validation work or an
    /// invalid loop.
    pub fn admit(&self) -> Result<AdmittedTrimmedPatch<'_>, NurbsError> {
        let mut work_remaining = TRIM_CLASSIFY_MAX_WORK_UNITS;
        self.admit_with_budget(&mut work_remaining)
    }

    /// Certified classification of a parameter-space point.
    ///
    /// # Errors
    /// Propagates structural, defensive work/memory, allocation, and exact
    /// rational-domain refusals.
    pub fn classify(&self, q: [Rat; 2]) -> Result<Classification, NurbsError> {
        self.classify_box(q, q)
    }

    /// Certified classification of every point in a closed parameter-space
    /// box. A verdict is returned only after every trim-curve Bézier hull is
    /// separated from the entire box, which proves that winding is constant
    /// throughout the connected box. Otherwise bounded subdivision returns
    /// [`Classification::Boundary`] rather than guessing from its corners or
    /// centre.
    ///
    /// # Errors
    /// Returns [`NurbsError::Domain`] for an inverted box or defensive
    /// work/memory refusal, [`NurbsError::Exactness`] when an exact midpoint is
    /// not representable, and propagates structural subdivision errors.
    pub fn classify_box(&self, min: [Rat; 2], max: [Rat; 2]) -> Result<Classification, NurbsError> {
        validate_classification_box(min, max)?;
        let mut work_remaining = TRIM_CLASSIFY_MAX_WORK_UNITS;
        let admitted = self.admit_with_budget(&mut work_remaining)?;
        admitted.classify_box_with_budget(min, max, &mut work_remaining)
    }
}

impl<'a> AdmittedTrimmedPatch<'a> {
    /// The exact immutable source bound to this view.
    #[must_use]
    pub const fn source(&self) -> &'a TrimmedPatch {
        self.inner
    }

    /// Borrow the sealed, already-validated loops.
    #[must_use]
    pub fn loops(&self) -> &'a [TrimLoop] {
        &self.inner.loops
    }

    /// Iterate over already-validated loop views bound to this exact patch
    /// generation.
    pub fn admitted_loops(&self) -> impl ExactSizeIterator<Item = AdmittedTrimLoop<'a>> + 'a {
        let loops: &'a [TrimLoop] = &self.inner.loops;
        loops.iter().map(|inner| AdmittedTrimLoop { inner })
    }

    /// Exact-subdivision depth before an ambiguous query becomes `Boundary`.
    #[must_use]
    pub const fn max_subdivision(&self) -> u32 {
        self.inner.max_subdivision
    }

    /// Certified point classification reusing this exact patch admission.
    ///
    /// # Errors
    /// Propagates checked work, retained-memory, allocation, structural, and
    /// exact rational-domain refusals.
    pub fn classify(&self, q: [Rat; 2]) -> Result<Classification, NurbsError> {
        self.classify_box(q, q)
    }

    /// Certified connected-box classification reusing this exact patch
    /// admission.
    ///
    /// # Errors
    /// Returns [`NurbsError::Domain`] for an inverted box or defensive
    /// work/memory refusal and [`NurbsError::Exactness`] when an exact midpoint
    /// is not representable.
    pub fn classify_box(&self, min: [Rat; 2], max: [Rat; 2]) -> Result<Classification, NurbsError> {
        validate_classification_box(min, max)?;
        let mut work_remaining = TRIM_CLASSIFY_MAX_WORK_UNITS;
        self.classify_box_with_budget(min, max, &mut work_remaining)
    }

    fn classify_box_with_budget(
        &self,
        min: [Rat; 2],
        max: [Rat; 2],
        work_remaining: &mut u128,
    ) -> Result<Classification, NurbsError> {
        spend_trim_work(
            work_remaining,
            self.loops().len() as u128,
            "persistent trim-source retained-byte accounting",
        )?;
        let mut persistent_source_bytes = 0u128;
        for trim_loop in self.admitted_loops() {
            let curve = trim_loop.curve();
            let curve_bytes = trim_curve_storage_bytes(
                curve.knots().knots().len(),
                curve.homogeneous_control_points().len(),
            )?;
            persistent_source_bytes = persistent_source_bytes
                .checked_add(curve_bytes)
                .ok_or_else(|| NurbsError::Domain {
                    what: "aggregate trim-source retained bytes overflow u128".to_string(),
                })?;
        }
        enforce_trim_retained_bytes(persistent_source_bytes, "persistent source")?;
        spend_trim_work(
            work_remaining,
            TRIM_EXACT_MIDPOINT_WORK_UNITS * 2,
            "classification witness midpoints",
        )?;
        let witness = [
            exact_midpoint(min[0], max[0], "classification witness u")?,
            exact_midpoint(min[1], max[1], "classification witness v")?,
        ];
        let mut winding = 0i64;
        for trim_loop in self.admitted_loops() {
            match loop_winding_box(
                trim_loop.curve(),
                min,
                max,
                witness,
                self.max_subdivision(),
                persistent_source_bytes,
                work_remaining,
            )? {
                Some(loop_winding) => {
                    winding =
                        winding
                            .checked_add(loop_winding)
                            .ok_or_else(|| NurbsError::Domain {
                                what: "aggregate trim winding overflows i64".to_string(),
                            })?;
                }
                None => return Ok(Classification::Boundary),
            }
        }
        Ok(if winding != 0 {
            Classification::Inside
        } else {
            Classification::Outside
        })
    }
}

fn validate_classification_box(min: [Rat; 2], max: [Rat; 2]) -> Result<(), NurbsError> {
    if min[0] > max[0] || min[1] > max[1] {
        return Err(NurbsError::Domain {
            what: "trim classification box must be componentwise ordered".to_string(),
        });
    }
    Ok(())
}

fn trim_loop_validation_work(curve: &NurbsCurve<Rat, 2>) -> Result<u128, NurbsError> {
    let control_components =
        (curve.cpw.len() as u128)
            .checked_mul(4)
            .ok_or_else(|| NurbsError::Domain {
                what: "trim control-validation accounting overflows u128".to_string(),
            })?;
    let order = (curve.knots.degree as u128)
        .checked_add(1)
        .ok_or_else(|| NurbsError::Domain {
            what: "trim order-validation accounting overflows u128".to_string(),
        })?;
    let basis_triangle = order.checked_mul(order).ok_or_else(|| NurbsError::Domain {
        what: "trim basis-validation accounting overflows u128".to_string(),
    })?;
    let scanned_entries = (curve.knots.knots.len() as u128)
        .checked_add(control_components)
        .and_then(|work| work.checked_add(basis_triangle))
        .ok_or_else(|| NurbsError::Domain {
            what: "trim structure-validation accounting overflows u128".to_string(),
        })?;
    // Closure evaluates both endpoints through one admitted curve. Eight scans
    // remains a conservative legacy charge for closure, basis work, projection,
    // and the full-break continuity walk.
    scanned_entries
        .checked_mul(8)
        .map(|work| work.max(TRIM_MIN_LOOP_VALIDATION_WORK_UNITS))
        .ok_or_else(|| NurbsError::Domain {
            what: "trim repeated-validation accounting overflows u128".to_string(),
        })
}

fn exact_midpoint(left: Rat, right: Rat, stage: &str) -> Result<Rat, NurbsError> {
    left.checked_midpoint(right)
        .ok_or_else(|| NurbsError::Exactness {
            what: format!("trim {stage} midpoint exceeds the exact i128 rational domain"),
        })
}

fn trim_curve_storage_bytes(knot_count: usize, control_count: usize) -> Result<u128, NurbsError> {
    let knot_bytes = (knot_count as u128)
        .checked_mul(core::mem::size_of::<Rat>() as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "trim knot-storage accounting overflows u128".to_string(),
        })?;
    let control_bytes = (control_count as u128)
        .checked_mul(core::mem::size_of::<[Rat; 4]>() as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "trim control-storage accounting overflows u128".to_string(),
        })?;
    knot_bytes
        .checked_add(control_bytes)
        .ok_or_else(|| NurbsError::Domain {
            what: "trim curve-storage accounting overflows u128".to_string(),
        })
}

fn enforce_trim_retained_bytes(retained_bytes: u128, stage: &str) -> Result<(), NurbsError> {
    if retained_bytes > TRIM_CLASSIFY_MAX_RETAINED_BYTES {
        return Err(NurbsError::Domain {
            what: format!(
                "trim {stage} can retain {retained_bytes} bytes above defensive ceiling {TRIM_CLASSIFY_MAX_RETAINED_BYTES}"
            ),
        });
    }
    Ok(())
}

fn trim_bezier_conversion_plan(
    curve: AdmittedNurbsCurve<'_, Rat, 2>,
    persistent_source_bytes: u128,
    operation_source_is_persistent: bool,
) -> Result<BezierConversionPlan, NurbsError> {
    let plan = curve.bezier_conversion_plan()?;
    let operation_source_bytes = trim_curve_storage_bytes(
        curve.knots().knots().len(),
        curve.homogeneous_control_points().len(),
    )?;
    let additional_live_bytes = if operation_source_is_persistent {
        persistent_source_bytes
            .checked_sub(operation_source_bytes)
            .ok_or_else(|| NurbsError::Domain {
                what: "persistent trim-source bytes omit the active source curve".to_string(),
            })?
    } else {
        persistent_source_bytes
    };
    let conversion_peak = additional_live_bytes
        .checked_add(operation_source_bytes)
        .and_then(|bytes| bytes.checked_add(plan.peak_allocated_bytes))
        .ok_or_else(|| NurbsError::Domain {
            what: "trim Bezier conversion peak retained-byte accounting overflows u128".to_string(),
        })?;
    let span_capacity = plan
        .final_control_count
        .checked_sub(curve.knots().degree())
        .ok_or_else(|| NurbsError::Structure {
            what: "trim projected Bezier degree exceeds final control count".to_string(),
        })?;
    let classification_peak = (span_capacity as u128)
        .checked_mul(core::mem::size_of::<SpanBox<Rat, 2>>() as u128)
        .and_then(|box_bytes| {
            (span_capacity as u128)
                .checked_mul(core::mem::size_of::<(Rat, Rat)>() as u128)
                .and_then(|interval_bytes| box_bytes.checked_add(interval_bytes))
        })
        .and_then(|scratch| plan.converted_bytes.checked_add(scratch))
        .and_then(|bytes| bytes.checked_add(persistent_source_bytes))
        .ok_or_else(|| NurbsError::Domain {
            what: "trim converted classification peak accounting overflows u128".to_string(),
        })?;
    enforce_trim_retained_bytes(
        conversion_peak.max(classification_peak),
        "Bezier conversion/classification",
    )?;
    Ok(plan)
}

fn preflight_trim_bezier_conversion(
    curve: AdmittedNurbsCurve<'_, Rat, 2>,
    persistent_source_bytes: u128,
    operation_source_is_persistent: bool,
    work_remaining: &mut u128,
) -> Result<(BezierConversionPlan, u128), NurbsError> {
    let pre_scan_work = curve.bezier_pre_scan_work()?;
    spend_trim_work(
        work_remaining,
        pre_scan_work,
        "Bezier conversion-plan knot scan",
    )?;
    let plan = trim_bezier_conversion_plan(
        curve,
        persistent_source_bytes,
        operation_source_is_persistent,
    )?;
    let unspent_work =
        plan.work_units
            .checked_sub(pre_scan_work)
            .ok_or_else(|| NurbsError::Domain {
                what: "trim Bezier post-scan work accounting is inconsistent".to_string(),
            })?;
    Ok((plan, unspent_work))
}

fn trim_classification_pass_work(control_count: usize, degree: usize) -> Result<u128, NurbsError> {
    let span_count = control_count
        .checked_sub(degree)
        .ok_or_else(|| NurbsError::Structure {
            what: "trim Bezier degree exceeds its admitted control count".to_string(),
        })?;
    let order = degree.checked_add(1).ok_or_else(|| NurbsError::Domain {
        what: "trim Bezier order overflows usize".to_string(),
    })?;
    (span_count as u128)
        .checked_mul(order as u128)
        .and_then(|visits| visits.checked_mul(TRIM_SPAN_BOX_WORK_PER_CONTROL))
        .and_then(|work| {
            (span_count as u128)
                .checked_mul(2)
                .and_then(|traversal| work.checked_add(traversal))
        })
        .and_then(|work| {
            (control_count as u128)
                .checked_mul(TRIM_WINDING_WORK_PER_CONTROL)
                .and_then(|winding| work.checked_add(winding))
        })
        .ok_or_else(|| NurbsError::Domain {
            what: "trim span/winding work accounting overflows u128".to_string(),
        })
}

fn preflight_trim_span_scratch(
    curve: AdmittedNurbsCurve<'_, Rat, 2>,
    persistent_source_bytes: u128,
) -> Result<usize, NurbsError> {
    let knots = curve.knots();
    let control_count = curve.homogeneous_control_points().len();
    let span_capacity =
        control_count
            .checked_sub(knots.degree())
            .ok_or_else(|| NurbsError::Structure {
                what: "trim Bezier degree exceeds its admitted control count".to_string(),
            })?;
    let curve_bytes = trim_curve_storage_bytes(knots.knots().len(), control_count)?;
    let box_bytes = (span_capacity as u128)
        .checked_mul(core::mem::size_of::<SpanBox<Rat, 2>>() as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "trim span-box retained-byte accounting overflows u128".to_string(),
        })?;
    let interval_bytes = (span_capacity as u128)
        .checked_mul(core::mem::size_of::<(Rat, Rat)>() as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "trim offending-interval retained-byte accounting overflows u128".to_string(),
        })?;
    let peak = curve_bytes
        .checked_add(box_bytes)
        .and_then(|bytes| bytes.checked_add(interval_bytes))
        .and_then(|bytes| bytes.checked_add(persistent_source_bytes))
        .ok_or_else(|| NurbsError::Domain {
            what: "trim span-classification peak retained-byte accounting overflows u128"
                .to_string(),
        })?;
    enforce_trim_retained_bytes(peak, "span classification")?;
    Ok(span_capacity)
}

fn projected_subdivision_work(
    curve: AdmittedNurbsCurve<'_, Rat, 2>,
    offending_count: usize,
) -> Result<(u128, u128, usize), NurbsError> {
    let knots = curve.knots();
    let conversion_insertions = offending_count
        .checked_mul(
            knots
                .degree()
                .checked_sub(1)
                .ok_or_else(|| NurbsError::Structure {
                    what: "trim subdivision requires a positive spline degree".to_string(),
                })?,
        )
        .ok_or_else(|| NurbsError::Domain {
            what: "trim projected Bezier insertion count overflows usize".to_string(),
        })?;
    let midpoint_work = (offending_count as u128)
        .checked_mul(TRIM_EXACT_MIDPOINT_WORK_UNITS)
        .ok_or_else(|| NurbsError::Domain {
            what: "trim exact-midpoint work overflows u128".to_string(),
        })?;
    let refinement_work = curve
        .projected_refinement_work(offending_count, conversion_insertions)?
        .checked_add(midpoint_work)
        .ok_or_else(|| NurbsError::Domain {
            what: "trim midpoint/refinement work overflows u128".to_string(),
        })?;
    let total_growth = offending_count
        .checked_add(conversion_insertions)
        .ok_or_else(|| NurbsError::Domain {
            what: "trim projected subdivision growth overflows usize".to_string(),
        })?;
    let final_control_count = curve
        .homogeneous_control_points()
        .len()
        .checked_add(total_growth)
        .ok_or_else(|| NurbsError::Domain {
            what: "trim projected control count overflows usize".to_string(),
        })?;
    let final_span_work = trim_classification_pass_work(final_control_count, knots.degree())?;
    Ok((refinement_work, final_span_work, conversion_insertions))
}

fn preflight_trim_subdivision_retained(
    curve: AdmittedNurbsCurve<'_, Rat, 2>,
    offending_count: usize,
    interval_capacity: usize,
    persistent_source_bytes: u128,
) -> Result<(), NurbsError> {
    let knots = curve.knots();
    let growth_per_span = knots.degree().max(1);
    let final_growth = offending_count
        .checked_mul(growth_per_span)
        .ok_or_else(|| NurbsError::Domain {
            what: "trim retained subdivision growth overflows usize".to_string(),
        })?;
    let midpoint_knot_count = knots
        .knots()
        .len()
        .checked_add(offending_count)
        .ok_or_else(|| NurbsError::Domain {
            what: "trim midpoint knot count overflows usize".to_string(),
        })?;
    let midpoint_control_count = curve
        .homogeneous_control_points()
        .len()
        .checked_add(offending_count)
        .ok_or_else(|| NurbsError::Domain {
            what: "trim midpoint control count overflows usize".to_string(),
        })?;
    let final_knot_count = knots
        .knots()
        .len()
        .checked_add(final_growth)
        .ok_or_else(|| NurbsError::Domain {
            what: "trim converted knot count overflows usize".to_string(),
        })?;
    let final_control_count = curve
        .homogeneous_control_points()
        .len()
        .checked_add(final_growth)
        .ok_or_else(|| NurbsError::Domain {
            what: "trim converted control count overflows usize".to_string(),
        })?;
    let midpoint_bytes = trim_curve_storage_bytes(midpoint_knot_count, midpoint_control_count)?;
    let final_bytes = trim_curve_storage_bytes(final_knot_count, final_control_count)?;
    let interval_bytes = (interval_capacity as u128)
        .checked_mul(core::mem::size_of::<(Rat, Rat)>() as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "trim live interval retained-byte accounting overflows u128".to_string(),
        })?;
    let insertion_peak = midpoint_bytes
        .checked_mul(2)
        .and_then(|bytes| bytes.checked_add(interval_bytes))
        .and_then(|bytes| bytes.checked_add(persistent_source_bytes))
        .ok_or_else(|| NurbsError::Domain {
            what: "trim midpoint-insertion peak retained-byte accounting overflows u128"
                .to_string(),
        })?;
    let conversion_peak = final_bytes
        .checked_mul(2)
        .and_then(|allocated| midpoint_bytes.checked_add(allocated))
        .and_then(|bytes| bytes.checked_add(interval_bytes))
        .and_then(|bytes| bytes.checked_add(persistent_source_bytes))
        .ok_or_else(|| NurbsError::Domain {
            what: "trim subdivision-conversion peak retained-byte accounting overflows u128"
                .to_string(),
        })?;
    let span_capacity = final_control_count
        .checked_sub(knots.degree())
        .ok_or_else(|| NurbsError::Structure {
            what: "trim projected degree exceeds final control count".to_string(),
        })?;
    let box_bytes = (span_capacity as u128)
        .checked_mul(core::mem::size_of::<SpanBox<Rat, 2>>() as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "trim projected span-box bytes overflow u128".to_string(),
        })?;
    let final_interval_bytes = (span_capacity as u128)
        .checked_mul(core::mem::size_of::<(Rat, Rat)>() as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "trim projected interval bytes overflow u128".to_string(),
        })?;
    let classification_peak = final_bytes
        .checked_add(box_bytes)
        .and_then(|bytes| bytes.checked_add(final_interval_bytes))
        .and_then(|bytes| bytes.checked_add(persistent_source_bytes))
        .ok_or_else(|| NurbsError::Domain {
            what: "trim projected classification bytes overflow u128".to_string(),
        })?;
    enforce_trim_retained_bytes(
        insertion_peak.max(conversion_peak).max(classification_peak),
        "subdivision/conversion",
    )
}

/// Certified winding number of one closed rational curve about `q`, or
/// `None` when `q` cannot be separated from the curve within the
/// subdivision budget.
fn loop_winding_box(
    curve: AdmittedNurbsCurve<'_, Rat, 2>,
    query_min: [Rat; 2],
    query_max: [Rat; 2],
    witness: [Rat; 2],
    max_depth: u32,
    persistent_source_bytes: u128,
    work_remaining: &mut u128,
) -> Result<Option<i64>, NurbsError> {
    // Work in Bézier form so each span's control hull tightly bounds it.
    let (initial_plan, initial_conversion_work) =
        preflight_trim_bezier_conversion(curve, persistent_source_bytes, true, work_remaining)?;
    let initial_span_work =
        trim_classification_pass_work(initial_plan.final_control_count, curve.knots().degree())?;
    require_trim_work(
        work_remaining,
        initial_conversion_work
            .checked_add(initial_span_work)
            .ok_or_else(|| NurbsError::Domain {
                what: "initial trim conversion/span work overflows u128".to_string(),
            })?,
        "initial Bezier conversion and first span-box construction",
    )?;
    spend_trim_work(
        work_remaining,
        initial_conversion_work,
        "initial Bézier conversion",
    )?;
    let mut work = curve.to_bezier_form()?;
    let mut depth = 0u32;
    loop {
        let admitted_work = work.admitted_after_validation();
        spend_trim_work(
            work_remaining,
            trim_classification_pass_work(
                admitted_work.homogeneous_control_points().len(),
                admitted_work.knots().degree(),
            )?,
            "span-box construction",
        )?;
        let span_capacity = preflight_trim_span_scratch(admitted_work, persistent_source_bytes)?;
        let boxes = admitted_work.span_boxes()?;
        let mut offending = Vec::new();
        offending
            .try_reserve_exact(span_capacity)
            .map_err(|_| NurbsError::Domain {
                what: "trim offending-interval allocation was refused".to_string(),
            })?;
        for &(min, max, t0, t1) in &boxes {
            if max[0] >= query_min[0]
                && min[0] <= query_max[0]
                && max[1] >= query_min[1]
                && min[1] <= query_max[1]
            {
                offending.push((t0, t1));
            }
        }
        if offending.is_empty() {
            // Separated from the whole connected query box: winding is
            // constant throughout it, so one exact witness is sufficient.
            return Ok(Some(polygon_winding_homogeneous(
                admitted_work.homogeneous_control_points(),
                witness,
            )));
        }
        if depth >= max_depth {
            return Ok(None);
        }
        let (future_work, next_span_work, expected_conversion_insertions) =
            projected_subdivision_work(admitted_work, offending.len())?;
        require_trim_work(
            work_remaining,
            future_work
                .checked_add(next_span_work)
                .ok_or_else(|| NurbsError::Domain {
                    what: "trim subdivision/downstream work overflows u128".to_string(),
                })?,
            "midpoint subdivision, Bezier reconversion, and next span-box construction",
        )?;
        spend_trim_work(
            work_remaining,
            future_work,
            "midpoint subdivision and Bezier reconversion",
        )?;
        let interval_capacity = offending.capacity();
        preflight_trim_subdivision_retained(
            admitted_work,
            offending.len(),
            interval_capacity,
            persistent_source_bytes,
        )?;
        drop(boxes);
        for (t0, t1) in offending {
            let admitted_work = work.admitted_after_validation();
            let mid = exact_midpoint(t0, t1, "subdivision parameter")?;
            // Exact midpoint insertion splits the offending span.
            work = admitted_work.insert_knot(mid)?;
        }
        let admitted_work = work.admitted_after_validation();
        // The aggregate refinement charge was consumed before the first
        // midpoint insertion and includes this external plan scan.
        let conversion_plan =
            trim_bezier_conversion_plan(admitted_work, persistent_source_bytes, false)?;
        if conversion_plan.insertions != expected_conversion_insertions {
            return Err(NurbsError::Structure {
                what: format!(
                    "trim midpoint refinement projected {expected_conversion_insertions} Bezier insertions but derived generation requires {}",
                    conversion_plan.insertions
                ),
            });
        }
        work = admitted_work.to_bezier_form()?;
        depth = depth.checked_add(1).ok_or_else(|| NurbsError::Domain {
            what: "trim subdivision depth overflows u32".to_string(),
        })?;
    }
}

fn require_trim_work(remaining: &u128, requested: u128, stage: &str) -> Result<(), NurbsError> {
    if requested > *remaining {
        return Err(NurbsError::Domain {
            what: format!(
                "trim {stage} requests {requested} work units with only {remaining} remaining from the {TRIM_CLASSIFY_MAX_WORK_UNITS}-unit defensive budget"
            ),
        });
    }
    Ok(())
}

fn spend_trim_work(remaining: &mut u128, requested: u128, stage: &str) -> Result<(), NurbsError> {
    require_trim_work(remaining, requested, stage)?;
    *remaining -= requested;
    Ok(())
}

/// EXACT winding number of the Cartesian control polygon without allocating
/// a projected copy. Admission guarantees positive weights and at least one
/// control.
fn polygon_winding_homogeneous(cpw: &[[Rat; 4]], q: [Rat; 2]) -> i64 {
    let mut winding = 0i64;
    for index in 0..cpw.len() {
        let a_h = cpw[index];
        let b_h = cpw[(index + 1) % cpw.len()];
        let a = [a_h[0] / a_h[3], a_h[1] / a_h[3]];
        let b = [b_h[0] / b_h[3], b_h[1] / b_h[3]];
        // Upward crossing: a.y <= q.y < b.y and q strictly left of ab.
        let orient = (b[0] - a[0]) * (q[1] - a[1]) - (q[0] - a[0]) * (b[1] - a[1]);
        if a[1] <= q[1] && q[1] < b[1] && orient > Rat::int(0) {
            winding += 1;
        } else if b[1] <= q[1] && q[1] < a[1] && orient < Rat::int(0) {
            winding -= 1;
        }
    }
    winding
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basis::KnotVector;
    use asupersync::types::Budget;
    use fs_exec::{CancelGate, ExecMode, StreamKey};

    fn with_trim_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
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
                    seed: 0x7A1C_100F,
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

    fn point_trim_loop() -> TrimLoop {
        let knots = KnotVector::new(vec![Rat::int(0), Rat::int(0), Rat::int(1), Rat::int(1)], 1)
            .expect("point-loop knots");
        let curve = NurbsCurve::new(knots, &[[Rat::int(0), Rat::int(0)]; 2], &[Rat::int(1); 2])
            .expect("point-loop curve");
        TrimLoop::new(curve).expect("closed point loop")
    }

    fn long_trim_loop() -> TrimLoop {
        let mut knots = vec![Rat::int(0); 2];
        knots.extend((1..=128).map(|numerator| Rat::new(numerator, 129)));
        knots.extend([Rat::int(1); 2]);
        let knots = KnotVector::new(knots, 1).expect("long loop knots");
        let points = vec![[Rat::int(0), Rat::int(0)]; 130];
        let weights = vec![Rat::int(1); 130];
        let curve = NurbsCurve::new(knots, &points, &weights).expect("long point-loop curve");
        TrimLoop::new(curve).expect("closed long point loop")
    }

    #[test]
    fn trim_loop_admission_with_cx_is_transactional_and_lifetime_bound() {
        let trim_loop = point_trim_loop();
        with_trim_cx(true, |cx| {
            assert!(matches!(
                trim_loop
                    .admit_with_cx(cx)
                    .expect("valid pre-cancelled loop"),
                TrimLoopAdmissionRun::Cancelled
            ));
        });
        with_trim_cx(false, |cx| {
            let TrimLoopAdmissionRun::Complete { admitted } = trim_loop
                .admit_with_cx(cx)
                .expect("active trim-loop admission")
            else {
                panic!("active trim-loop admission must complete");
            };
            assert!(core::ptr::eq(admitted.source(), &trim_loop));
            assert_eq!(admitted.curve().homogeneous_control_points().len(), 2);
        });

        let open_curve = NurbsCurve::new(
            KnotVector::new(vec![Rat::int(0), Rat::int(0), Rat::int(1), Rat::int(1)], 1)
                .expect("open-loop knots"),
            &[[Rat::int(0), Rat::int(0)], [Rat::int(1), Rat::int(0)]],
            &[Rat::int(1); 2],
        )
        .expect("open curve");
        let open_loop = TrimLoop { curve: open_curve };
        with_trim_cx(false, |cx| {
            assert!(matches!(
                open_loop.admit_with_cx(cx),
                Err(NurbsError::Structure { ref what }) if what.contains("close exactly")
            ));
        });
    }

    #[test]
    fn trim_loop_continuity_scan_cancels_at_a_deterministic_stride() {
        let trim_loop = long_trim_loop();
        let curve = trim_loop.curve.admit().expect("admitted long trim curve");
        let (lo, hi) = curve.knots().domain();
        let start = curve.eval(lo).expect("long loop start");
        let end = curve.eval(hi).expect("long loop end");
        let run = || {
            let mut polls = 0usize;
            let outcome = validate_trim_loop_after_endpoints_with_poll(curve, start, end, || {
                polls += 1;
                polls == 2
            })
            .expect("valid cancellable continuity scan");
            (outcome, polls)
        };
        assert_eq!(run(), run());
        assert_eq!(run(), (TrimLoopValidationOutcome::Cancelled, 2));
    }

    #[test]
    fn trim_loop_admission_final_checkpoint_gates_publication() {
        let trim_loop = point_trim_loop();
        let curve = trim_loop.curve.admit().expect("admitted point-loop curve");
        let (lo, hi) = curve.knots().domain();
        let start = curve.eval(lo).expect("point-loop start");
        let end = curve.eval(hi).expect("point-loop end");
        let mut total_polls = 0usize;
        assert_eq!(
            validate_trim_loop_after_endpoints_with_poll(curve, start, end, || {
                total_polls += 1;
                false
            })
            .expect("healthy continuity scan"),
            TrimLoopValidationOutcome::Complete
        );
        assert!(total_polls >= 2);

        let mut replay_polls = 0usize;
        assert_eq!(
            validate_trim_loop_after_endpoints_with_poll(curve, start, end, || {
                replay_polls += 1;
                replay_polls == total_polls
            })
            .expect("publication cancellation"),
            TrimLoopValidationOutcome::Cancelled
        );
        assert_eq!(replay_polls, total_polls);
    }

    #[test]
    fn inverted_box_refusal_precedes_trim_admission() {
        let knots = KnotVector::new(vec![Rat::int(0), Rat::int(0), Rat::int(1), Rat::int(1)], 1)
            .expect("line knots");
        let malformed_loop = TrimLoop {
            curve: NurbsCurve {
                knots,
                cpw: Vec::new(),
            },
        };
        let patch = TrimmedPatch::new(vec![malformed_loop]);
        let error = patch
            .classify_box([Rat::int(1), Rat::int(0)], [Rat::int(0), Rat::int(1)])
            .expect_err("inverted box must refuse before malformed loop admission");
        assert!(matches!(
            error,
            NurbsError::Domain { ref what } if what.contains("componentwise ordered")
        ));
    }

    #[test]
    fn empty_patch_copy_preserves_sealed_configuration() {
        let patch = TrimmedPatch::with_max_subdivision(Vec::new(), 7);
        assert_eq!(patch.try_clone().expect("fallible patch copy"), patch);
    }

    #[test]
    fn classification_envelopes_refuse_before_runtime_allocation() {
        assert!(
            enforce_trim_retained_bytes(TRIM_CLASSIFY_MAX_RETAINED_BYTES, "test boundary").is_ok()
        );
        assert!(matches!(
            enforce_trim_retained_bytes(
                TRIM_CLASSIFY_MAX_RETAINED_BYTES + 1,
                "test boundary"
            ),
            Err(NurbsError::Domain { ref what }) if what.contains("retain")
        ));
        let synthetic_bytes = trim_curve_storage_bytes(usize::MAX, usize::MAX)
            .expect("usize-sized counts fit u128 accounting");
        assert!(matches!(
            enforce_trim_retained_bytes(synthetic_bytes, "synthetic counts"),
            Err(NurbsError::Domain { ref what }) if what.contains("retain")
        ));

        let mut no_work = 0;
        assert!(matches!(
            spend_trim_work(&mut no_work, 1, "test work precedence"),
            Err(NurbsError::Domain { ref what }) if what.contains("work")
        ));
    }

    #[test]
    fn conversion_plan_is_charged_to_the_shared_trim_budget() {
        let degree = 10usize;
        let mut knots = Vec::new();
        for _ in 0..=degree {
            knots.push(Rat::int(0));
        }
        for numerator in 1..20 {
            knots.push(Rat::new(numerator, 20));
        }
        for _ in 0..=degree {
            knots.push(Rat::int(1));
        }
        let knots = KnotVector::new(knots, degree).expect("high-degree trim knots");
        let points = vec![[Rat::int(0), Rat::int(0)]; 30];
        let weights = vec![Rat::int(1); 30];
        let curve = NurbsCurve::new(knots, &points, &weights).expect("high-degree trim curve");
        let persistent_source_bytes = trim_curve_storage_bytes(
            curve.knots().knots().len(),
            curve.homogeneous_control_points().len(),
        )
        .expect("source retained-byte accounting");
        let plan = trim_bezier_conversion_plan(
            curve.admit().expect("admitted high-degree trim curve"),
            persistent_source_bytes,
            true,
        )
        .expect("conversion plan remains inside the curve-local ceiling");
        assert!(
            plan.work_units > TRIM_CLASSIFY_MAX_WORK_UNITS,
            "fixture must exceed the smaller aggregate trim budget"
        );
        let patch = TrimmedPatch::new(vec![
            TrimLoop::new(curve).expect("closed high-degree trim loop"),
        ]);
        let error = patch
            .classify([Rat::int(2), Rat::int(2)])
            .expect_err("aggregate trim budget must refuse before conversion allocation");
        assert!(matches!(
            error,
            NurbsError::Domain { ref what }
                if what.contains("initial Bezier conversion") && what.contains("work")
        ));
    }

    #[test]
    fn extreme_representable_box_midpoint_does_not_overflow() {
        let patch = TrimmedPatch::new(Vec::new());
        let result = patch
            .classify_box(
                [Rat::new(i128::MAX - 2, 1), Rat::int(0)],
                [Rat::new(i128::MAX, 1), Rat::int(0)],
            )
            .expect("representable exact midpoint");
        assert_eq!(result, Classification::Outside);
    }

    #[test]
    fn unrepresentable_box_midpoint_is_a_typed_exactness_refusal() {
        let patch = TrimmedPatch::new(Vec::new());
        let error = patch
            .classify_box(
                [Rat::int(0), Rat::int(0)],
                [Rat::new(1, i128::MAX), Rat::int(0)],
            )
            .expect_err("reduced midpoint denominator exceeds i128");
        assert!(matches!(
            error,
            NurbsError::Exactness { ref what } if what.contains("midpoint")
        ));
    }
}
