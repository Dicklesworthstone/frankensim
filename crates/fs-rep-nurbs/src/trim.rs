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
    AdmittedNurbsCurve, BezierConversionPlan, CurveAdmissionRun, CurveCloneRun, CurveEvaluationRun,
    NurbsCurve, SpanBox,
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

#[derive(Debug, PartialEq, Eq)]
enum TrimWorkRun<T> {
    Complete(T),
    Cancelled,
}

#[derive(Debug, Clone, Copy)]
struct TrimClassificationQuery {
    min: [Rat; 2],
    max: [Rat; 2],
    witness: [Rat; 2],
    max_depth: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct TrimLoopReversalPlan {
    knot_count: usize,
    control_count: usize,
    degree: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct TrimmedPatchCopyPlan {
    loop_count: usize,
    #[cfg(test)]
    knot_count: usize,
    #[cfg(test)]
    control_count: usize,
    #[cfg(test)]
    work_units: u128,
    #[cfg(test)]
    retained_bytes: u128,
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

/// Transactional terminal state of cancellation-aware trim-loop construction.
#[must_use]
#[derive(Debug, PartialEq)]
pub enum TrimLoopConstructionRun {
    /// Validation completed and the sealed exact loop is safe to publish.
    Complete {
        /// Newly validated exact trim-loop generation.
        trim_loop: TrimLoop,
    },
    /// Cancellation was observed; the transferred curve was dropped without
    /// publishing a loop.
    Cancelled,
}

/// Transactional terminal state of a cancellation-aware fallible trim-loop
/// copy.
#[must_use]
#[derive(Debug, PartialEq)]
pub enum TrimLoopCloneRun {
    /// The complete sealed copy of the exact source representation.
    Complete {
        /// Copied exact trim-loop generation.
        trim_loop: TrimLoop,
    },
    /// Cancellation was observed; all partial copy storage was dropped.
    Cancelled,
}

/// Transactional terminal state of cancellation-aware trim-loop reversal.
#[must_use]
#[derive(Debug, PartialEq)]
pub enum TrimLoopReversalRun {
    /// The complete sealed and validated opposite-orientation loop.
    Complete {
        /// Reversed exact trim-loop generation.
        trim_loop: TrimLoop,
    },
    /// Cancellation was observed; all partial derived storage was dropped.
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

    /// Validate an owned exact curve as a closed trim loop with bounded
    /// cancellation polling.
    ///
    /// One `Cx` spans the existing curve/knot admission, both exact endpoint
    /// evaluations, full-break continuity traversal, and final owned
    /// publication. Construction allocates no derived loop payload; endpoint
    /// basis scratch retains its existing bounded allocation policy.
    /// Cancellation drops the caller-transferred curve without exposing a
    /// partially validated loop. Individual allocations, exact-rational
    /// operations, and destruction are not preemptible. This primitive does
    /// not consume the `Cx` budget or own request -> drain -> finalize
    /// semantics.
    ///
    /// # Errors
    /// Returns the synchronous constructor's work, allocation, structural,
    /// numeric-domain, closure, or continuity refusal when it wins before an
    /// observed cancellation.
    pub fn new_with_cx(
        curve: NurbsCurve<Rat, 2>,
        cx: &Cx<'_>,
    ) -> Result<TrimLoopConstructionRun, NurbsError> {
        let candidate = TrimLoop { curve };
        let validation = match candidate.admit_with_cx(cx)? {
            TrimLoopAdmissionRun::Complete { .. } => TrimLoopValidationOutcome::Complete,
            TrimLoopAdmissionRun::Cancelled => TrimLoopValidationOutcome::Cancelled,
        };
        let mut should_cancel = || cx.checkpoint().is_err();
        Ok(finish_trim_loop_construction_with_poll(
            candidate,
            validation,
            &mut should_cancel,
        ))
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

    /// Fallibly copy this sealed exact loop with bounded cancellation polling.
    ///
    /// The nested curve copy preserves its count-derived work and retained
    /// output refusals before cancellation, then one `Cx` spans both fallible
    /// allocations, ordered exact copies, and final loop publication. The
    /// immutable source is not revalidated. Individual allocator calls,
    /// exact-rational copies, and destructors are not preemptible. This
    /// primitive does not consume the `Cx` budget or own request -> drain ->
    /// finalize semantics.
    ///
    /// # Errors
    /// Returns the synchronous copy's work, retained-memory, or allocation
    /// refusal when it wins before an observed cancellation.
    pub fn try_clone_with_cx(&self, cx: &Cx<'_>) -> Result<TrimLoopCloneRun, NurbsError> {
        let curve_copy = self.curve.try_clone_with_cx(cx)?;
        let mut should_cancel = || cx.checkpoint().is_err();
        Ok(finish_trim_loop_clone_with_poll(
            curve_copy,
            &mut should_cancel,
        ))
    }

    /// The same loop with reversed orientation (holes are wound opposite
    /// to outers under the nonzero rule): control points reversed, knot
    /// vector mirrored about the domain.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] when checked work, retained storage, or a
    /// destination allocation is refused; [`NurbsError::Structure`] when the
    /// derived closure, continuity, knots, or control net are invalid.
    pub fn reversed_for_hole(&self) -> Result<TrimLoop, NurbsError> {
        self.admit()?.reversed_for_hole()
    }
}

fn finish_trim_loop_construction_with_poll(
    candidate: TrimLoop,
    validation: TrimLoopValidationOutcome,
    should_cancel: &mut impl FnMut() -> bool,
) -> TrimLoopConstructionRun {
    match validation {
        TrimLoopValidationOutcome::Complete => {
            if should_cancel() {
                TrimLoopConstructionRun::Cancelled
            } else {
                TrimLoopConstructionRun::Complete {
                    trim_loop: candidate,
                }
            }
        }
        TrimLoopValidationOutcome::Cancelled => TrimLoopConstructionRun::Cancelled,
    }
}

fn finish_trim_loop_clone_with_poll(
    curve_copy: CurveCloneRun<Rat, 2>,
    should_cancel: &mut impl FnMut() -> bool,
) -> TrimLoopCloneRun {
    match curve_copy {
        CurveCloneRun::Complete { curve } => {
            publish_trim_loop_clone_with_poll(curve, should_cancel)
        }
        CurveCloneRun::Cancelled => TrimLoopCloneRun::Cancelled,
    }
}

fn publish_trim_loop_clone_with_poll(
    curve: NurbsCurve<Rat, 2>,
    should_cancel: &mut impl FnMut() -> bool,
) -> TrimLoopCloneRun {
    if should_cancel() {
        return TrimLoopCloneRun::Cancelled;
    }
    TrimLoopCloneRun::Complete {
        trim_loop: TrimLoop { curve },
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

    /// Reverse this admitted loop without rescanning the immutable source.
    ///
    /// The derived loop is fully revalidated before publication.
    ///
    /// # Errors
    /// Returns checked work, retained-memory, allocation, or derived structural
    /// refusals.
    pub fn reversed_for_hole(&self) -> Result<TrimLoop, NurbsError> {
        let plan = preflight_trim_loop_reversal(*self)?;
        let mut never_cancel = || false;
        let candidate = match assemble_reversed_trim_loop_with_poll(*self, plan, &mut never_cancel)?
        {
            TrimWorkRun::Complete(candidate) => candidate,
            TrimWorkRun::Cancelled => {
                return Err(NurbsError::Domain {
                    what: "non-cancelling trim-loop reversal observed cancellation".to_string(),
                });
            }
        };
        candidate.admit()?;
        match publish_reversed_trim_loop_with_poll(candidate, &mut never_cancel) {
            TrimLoopReversalRun::Complete { trim_loop } => Ok(trim_loop),
            TrimLoopReversalRun::Cancelled => Err(NurbsError::Domain {
                what: "non-cancelling trim-loop reversal publication observed cancellation"
                    .to_string(),
            }),
        }
    }

    /// Reverse this admitted loop with bounded cancellation polling.
    ///
    /// Count-derived work and simultaneously-live retained storage are
    /// admitted before cancellation. One `Cx` then spans fallible knot/control
    /// allocation, exact knot mirroring, ordered copies, complete validation
    /// of the derived loop, and final owned publication. Cancellation exposes
    /// no partial loop. Individual allocator calls, exact-rational operations,
    /// and destructors are non-preemptible; this primitive does not consume
    /// the `Cx` budget or own request -> drain -> finalize semantics.
    ///
    /// # Errors
    /// Returns the synchronous reversal's work, memory, allocation, or
    /// structural refusal when it wins before an observed cancellation.
    pub fn reversed_for_hole_with_cx(
        &self,
        cx: &Cx<'_>,
    ) -> Result<TrimLoopReversalRun, NurbsError> {
        let plan = preflight_trim_loop_reversal(*self)?;
        let mut should_cancel = || cx.checkpoint().is_err();
        let candidate =
            match assemble_reversed_trim_loop_with_poll(*self, plan, &mut should_cancel)? {
                TrimWorkRun::Complete(candidate) => candidate,
                TrimWorkRun::Cancelled => return Ok(TrimLoopReversalRun::Cancelled),
            };
        match candidate.admit_with_cx(cx)? {
            TrimLoopAdmissionRun::Complete { .. } => {}
            TrimLoopAdmissionRun::Cancelled => return Ok(TrimLoopReversalRun::Cancelled),
        }
        Ok(publish_reversed_trim_loop_with_poll(
            candidate,
            &mut should_cancel,
        ))
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

/// Transactional terminal state of cancellation-aware certified trim
/// classification.
#[must_use]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrimClassificationRun {
    /// Classification completed and the entire terminal verdict is published.
    Complete {
        /// Certified in/out/boundary verdict.
        classification: Classification,
    },
    /// Cancellation was observed; no classification verdict was published.
    Cancelled,
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

/// Transactional terminal state of cancellation-aware trimmed-patch
/// admission.
#[must_use]
#[derive(Debug, Clone, Copy)]
pub enum TrimmedPatchAdmissionRun<'a> {
    /// Every exact loop in the immutable patch snapshot was fully validated.
    Complete {
        /// Lifetime-bound authority for the validated trimmed-patch
        /// generation.
        admitted: AdmittedTrimmedPatch<'a>,
    },
    /// Cancellation was observed; no admitted authority was published.
    Cancelled,
}

/// Transactional terminal state of a cancellation-aware fallible
/// trimmed-patch copy.
#[must_use]
#[derive(Debug, PartialEq)]
pub enum TrimmedPatchCloneRun {
    /// The complete sealed copy of the exact source representation.
    Complete {
        /// Copied trimmed-patch generation.
        trimmed_patch: TrimmedPatch,
    },
    /// Cancellation was observed; all partial nested copy storage was dropped.
    Cancelled,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TrimmedPatchValidationOutcome {
    Complete,
    Cancelled,
}

fn preflight_trimmed_patch_copy_counts(
    loop_count: usize,
    knot_count: usize,
    control_count: usize,
) -> Result<TrimmedPatchCopyPlan, NurbsError> {
    let work_units = (control_count as u128)
        .checked_mul(4)
        .and_then(|work| work.checked_add(knot_count as u128))
        .and_then(|work| work.checked_add((loop_count as u128).checked_mul(4)?))
        .and_then(|work| work.checked_add(2))
        .ok_or_else(|| NurbsError::Domain {
            what: "trimmed-patch copy-work accounting overflows u128".to_string(),
        })?;
    if work_units > TRIM_CLASSIFY_MAX_WORK_UNITS {
        return Err(NurbsError::Domain {
            what: format!(
                "trimmed-patch copy requests {work_units} work units above defensive ceiling {TRIM_CLASSIFY_MAX_WORK_UNITS}"
            ),
        });
    }

    let retained_bytes = (loop_count as u128)
        .checked_mul(core::mem::size_of::<TrimLoop>() as u128)
        .and_then(|bytes| {
            bytes
                .checked_add((knot_count as u128).checked_mul(core::mem::size_of::<Rat>() as u128)?)
        })
        .and_then(|bytes| {
            bytes.checked_add(
                (control_count as u128).checked_mul(core::mem::size_of::<[Rat; 4]>() as u128)?,
            )
        })
        .ok_or_else(|| NurbsError::Domain {
            what: "trimmed-patch copy retained-byte accounting overflows u128".to_string(),
        })?;
    if retained_bytes > TRIM_CLASSIFY_MAX_RETAINED_BYTES {
        return Err(NurbsError::Domain {
            what: format!(
                "trimmed-patch copy retains {retained_bytes} output bytes above defensive ceiling {TRIM_CLASSIFY_MAX_RETAINED_BYTES}"
            ),
        });
    }
    Ok(TrimmedPatchCopyPlan {
        loop_count,
        #[cfg(test)]
        knot_count,
        #[cfg(test)]
        control_count,
        #[cfg(test)]
        work_units,
        #[cfg(test)]
        retained_bytes,
    })
}

fn preflight_trimmed_patch_copy_with_poll(
    patch: &TrimmedPatch,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<TrimWorkRun<TrimmedPatchCopyPlan>, NurbsError> {
    // Refuse an impossible loop table in O(1), before a pre-existing
    // cancellation can mask the count-derived work/memory envelope.
    preflight_trimmed_patch_copy_counts(patch.loops.len(), 0, 0)?;
    if should_cancel() {
        return Ok(TrimWorkRun::Cancelled);
    }

    let mut knot_count = 0usize;
    let mut control_count = 0usize;
    let mut operations_since_poll = 0usize;
    for trim_loop in &patch.loops {
        knot_count = knot_count
            .checked_add(trim_loop.curve.knots.knots.len())
            .ok_or_else(|| NurbsError::Domain {
                what: "trimmed-patch copy knot-count accounting overflows usize".to_string(),
            })?;
        control_count = control_count
            .checked_add(trim_loop.curve.cpw.len())
            .ok_or_else(|| NurbsError::Domain {
                what: "trimmed-patch copy control-count accounting overflows usize".to_string(),
            })?;
        if trim_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(TrimWorkRun::Cancelled);
        }
    }
    let plan = preflight_trimmed_patch_copy_counts(patch.loops.len(), knot_count, control_count)?;
    if should_cancel() {
        return Ok(TrimWorkRun::Cancelled);
    }
    Ok(TrimWorkRun::Complete(plan))
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
        let mut never_cancel = || false;
        let mut admit_loop = |trim_loop: &TrimLoop| {
            trim_loop.admit()?;
            Ok(TrimmedPatchValidationOutcome::Complete)
        };
        match self.admit_with_budget_and_poll(work_remaining, &mut never_cancel, &mut admit_loop)? {
            TrimmedPatchAdmissionRun::Complete { admitted } => Ok(admitted),
            TrimmedPatchAdmissionRun::Cancelled => Err(NurbsError::Domain {
                what: "non-cancelling trimmed-patch admission observed cancellation".to_string(),
            }),
        }
    }

    fn admit_with_budget_and_poll<'a>(
        &'a self,
        work_remaining: &mut u128,
        should_cancel: &mut impl FnMut() -> bool,
        admit_loop: &mut impl FnMut(&TrimLoop) -> Result<TrimmedPatchValidationOutcome, NurbsError>,
    ) -> Result<TrimmedPatchAdmissionRun<'a>, NurbsError> {
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
        if should_cancel() {
            return Ok(TrimmedPatchAdmissionRun::Cancelled);
        }
        let mut validation_work = 0u128;
        let mut operations_since_poll = 0usize;
        for trim_loop in &self.loops {
            validation_work = validation_work
                .checked_add(trim_loop_validation_work(&trim_loop.curve)?)
                .ok_or_else(|| NurbsError::Domain {
                    what: "trim live-validation accounting overflows u128".to_string(),
                })?;
            if trim_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(TrimmedPatchAdmissionRun::Cancelled);
            }
        }
        spend_trim_work(work_remaining, validation_work, "live validation")?;
        if should_cancel() {
            return Ok(TrimmedPatchAdmissionRun::Cancelled);
        }
        operations_since_poll = 0;
        for trim_loop in &self.loops {
            match admit_loop(trim_loop)? {
                TrimmedPatchValidationOutcome::Complete => {}
                TrimmedPatchValidationOutcome::Cancelled => {
                    return Ok(TrimmedPatchAdmissionRun::Cancelled);
                }
            }
            if trim_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(TrimmedPatchAdmissionRun::Cancelled);
            }
        }
        if should_cancel() {
            return Ok(TrimmedPatchAdmissionRun::Cancelled);
        }
        Ok(TrimmedPatchAdmissionRun::Complete {
            admitted: AdmittedTrimmedPatch { inner: self },
        })
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
    /// [`NurbsError::Domain`] when checked aggregate work/retained bytes or a
    /// destination allocation is refused.
    pub fn try_clone(&self) -> Result<Self, NurbsError> {
        let mut never_cancel = || false;
        let mut clone_loop = |trim_loop: &TrimLoop| {
            Ok(TrimLoopCloneRun::Complete {
                trim_loop: trim_loop.try_clone()?,
            })
        };
        match self.try_clone_with_nested_and_poll(&mut never_cancel, &mut clone_loop)? {
            TrimmedPatchCloneRun::Complete { trimmed_patch } => Ok(trimmed_patch),
            TrimmedPatchCloneRun::Cancelled => Err(NurbsError::Domain {
                what: "non-cancelling trimmed-patch copy observed cancellation".to_string(),
            }),
        }
    }

    /// Fallibly copy this sealed exact patch with bounded cancellation polling.
    ///
    /// A count-only lower bound precedes cancellation. One fixed-stride
    /// metadata scan then admits aggregate nested copy work and a 64 MiB
    /// retained-output envelope covering the outer loop table plus every exact
    /// knot/control payload. The same `Cx` spans outer allocation, ordered
    /// nested loop copies, table moves, and final publication. The immutable
    /// source is not revalidated. Individual allocator calls, exact-rational
    /// copies, and nested destruction are not preemptible. This primitive does
    /// not consume the `Cx` budget or own request -> drain -> finalize
    /// semantics.
    ///
    /// # Errors
    /// Returns the synchronous copy's work, retained-memory, or allocation
    /// refusal when it wins before an observed cancellation.
    pub fn try_clone_with_cx(&self, cx: &Cx<'_>) -> Result<TrimmedPatchCloneRun, NurbsError> {
        let mut should_cancel = || cx.checkpoint().is_err();
        let mut clone_loop = |trim_loop: &TrimLoop| trim_loop.try_clone_with_cx(cx);
        self.try_clone_with_nested_and_poll(&mut should_cancel, &mut clone_loop)
    }

    fn try_clone_with_nested_and_poll(
        &self,
        should_cancel: &mut impl FnMut() -> bool,
        clone_loop: &mut impl FnMut(&TrimLoop) -> Result<TrimLoopCloneRun, NurbsError>,
    ) -> Result<TrimmedPatchCloneRun, NurbsError> {
        let plan = match preflight_trimmed_patch_copy_with_poll(self, should_cancel)? {
            TrimWorkRun::Complete(plan) => plan,
            TrimWorkRun::Cancelled => return Ok(TrimmedPatchCloneRun::Cancelled),
        };
        let mut loops = Vec::new();
        loops
            .try_reserve_exact(plan.loop_count)
            .map_err(|_| NurbsError::Domain {
                what: "trimmed-patch copy loop-table allocation was refused".to_string(),
            })?;
        if should_cancel() {
            return Ok(TrimmedPatchCloneRun::Cancelled);
        }
        let mut operations_since_poll = 0usize;
        for source_loop in &self.loops {
            let trim_loop = match clone_loop(source_loop)? {
                TrimLoopCloneRun::Complete { trim_loop } => trim_loop,
                TrimLoopCloneRun::Cancelled => return Ok(TrimmedPatchCloneRun::Cancelled),
            };
            loops.push(trim_loop);
            if trim_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(TrimmedPatchCloneRun::Cancelled);
            }
        }
        if should_cancel() {
            return Ok(TrimmedPatchCloneRun::Cancelled);
        }
        Ok(TrimmedPatchCloneRun::Complete {
            trimmed_patch: TrimmedPatch {
                loops,
                max_subdivision: self.max_subdivision,
            },
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

    /// Validate this exact immutable patch with bounded cancellation polling
    /// and publish only a lifetime-bound admitted view.
    ///
    /// The constant-time minimum loop-count work refusal precedes the first
    /// checkpoint. One `Cx` then spans the exact aggregate validation-work
    /// scan, every nested loop/curve admission, and final authority
    /// publication. Cancellation exposes no partially admitted loop table.
    /// This method does not consume the `Cx` budget or finalize its executor
    /// scope.
    ///
    /// # Errors
    /// Returns the synchronous admission's checked-work, knot, control,
    /// closure, continuity, and exact-arithmetic refusals when they win before
    /// an observed cancellation.
    pub fn admit_with_cx<'a>(
        &'a self,
        cx: &Cx<'_>,
    ) -> Result<TrimmedPatchAdmissionRun<'a>, NurbsError> {
        let mut work_remaining = TRIM_CLASSIFY_MAX_WORK_UNITS;
        let mut should_cancel = || cx.checkpoint().is_err();
        let mut admit_loop = |trim_loop: &TrimLoop| match trim_loop.admit_with_cx(cx)? {
            TrimLoopAdmissionRun::Complete { .. } => Ok(TrimmedPatchValidationOutcome::Complete),
            TrimLoopAdmissionRun::Cancelled => Ok(TrimmedPatchValidationOutcome::Cancelled),
        };
        self.admit_with_budget_and_poll(&mut work_remaining, &mut should_cancel, &mut admit_loop)
    }

    /// Certified classification of a parameter-space point.
    ///
    /// # Errors
    /// Propagates structural, defensive work/memory, allocation, and exact
    /// rational-domain refusals.
    pub fn classify(&self, q: [Rat; 2]) -> Result<Classification, NurbsError> {
        self.classify_box(q, q)
    }

    /// Certified point classification with bounded cancellation polling.
    ///
    /// # Errors
    /// Propagates the synchronous request, admission, work/memory,
    /// allocation, structural, and exact-arithmetic refusals when they win
    /// before an observed cancellation.
    pub fn classify_with_cx(
        &self,
        q: [Rat; 2],
        cx: &Cx<'_>,
    ) -> Result<TrimClassificationRun, NurbsError> {
        self.classify_box_with_cx(q, q, cx)
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

    /// Certified connected-box classification with bounded cancellation
    /// polling from source admission through terminal verdict publication.
    ///
    /// The componentwise box-order refusal and constant-time patch-admission
    /// work gate precede cancellation. One `Cx` then spans aggregate patch
    /// admission and the complete exact classification pipeline. Cancellation
    /// publishes neither partial admitted authority nor a partial verdict.
    /// This method does not consume the `Cx` budget or finalize its executor
    /// scope.
    ///
    /// # Errors
    /// Returns the synchronous path's request, work/memory, allocation,
    /// structural, and exact-arithmetic refusals when they win before an
    /// observed cancellation.
    pub fn classify_box_with_cx(
        &self,
        min: [Rat; 2],
        max: [Rat; 2],
        cx: &Cx<'_>,
    ) -> Result<TrimClassificationRun, NurbsError> {
        validate_classification_box(min, max)?;
        let mut work_remaining = TRIM_CLASSIFY_MAX_WORK_UNITS;
        let mut should_cancel = || cx.checkpoint().is_err();
        let mut admit_loop = |trim_loop: &TrimLoop| match trim_loop.admit_with_cx(cx)? {
            TrimLoopAdmissionRun::Complete { .. } => Ok(TrimmedPatchValidationOutcome::Complete),
            TrimLoopAdmissionRun::Cancelled => Ok(TrimmedPatchValidationOutcome::Cancelled),
        };
        let admitted = match self.admit_with_budget_and_poll(
            &mut work_remaining,
            &mut should_cancel,
            &mut admit_loop,
        )? {
            TrimmedPatchAdmissionRun::Complete { admitted } => admitted,
            TrimmedPatchAdmissionRun::Cancelled => return Ok(TrimClassificationRun::Cancelled),
        };
        admitted.classify_box_with_budget_and_poll(
            min,
            max,
            &mut work_remaining,
            &mut should_cancel,
        )
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

    /// Certified point classification with bounded cancellation polling while
    /// reusing this exact patch admission.
    ///
    /// # Errors
    /// Propagates the synchronous work/memory, allocation, structural, and
    /// exact-arithmetic refusals when they win before an observed
    /// cancellation.
    pub fn classify_with_cx(
        &self,
        q: [Rat; 2],
        cx: &Cx<'_>,
    ) -> Result<TrimClassificationRun, NurbsError> {
        self.classify_box_with_cx(q, q, cx)
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

    /// Certified connected-box classification with bounded cancellation
    /// polling while reusing this exact patch admission.
    ///
    /// The componentwise box-order refusal precedes the first checkpoint. The
    /// gate then spans retained-source accounting, exact witness construction,
    /// every loop conversion/subdivision/winding phase, cleanup, and terminal
    /// verdict publication. Cancellation publishes no partial verdict. This
    /// method does not consume the `Cx` budget or finalize its executor scope.
    ///
    /// # Errors
    /// Returns the synchronous path's request, work/memory, allocation,
    /// structural, and exact-arithmetic refusals when they win before an
    /// observed cancellation.
    pub fn classify_box_with_cx(
        &self,
        min: [Rat; 2],
        max: [Rat; 2],
        cx: &Cx<'_>,
    ) -> Result<TrimClassificationRun, NurbsError> {
        validate_classification_box(min, max)?;
        let mut work_remaining = TRIM_CLASSIFY_MAX_WORK_UNITS;
        let mut should_cancel = || cx.checkpoint().is_err();
        self.classify_box_with_budget_and_poll(min, max, &mut work_remaining, &mut should_cancel)
    }

    fn classify_box_with_budget(
        &self,
        min: [Rat; 2],
        max: [Rat; 2],
        work_remaining: &mut u128,
    ) -> Result<Classification, NurbsError> {
        let mut never_cancel = || false;
        match self.classify_box_with_budget_and_poll(min, max, work_remaining, &mut never_cancel)? {
            TrimClassificationRun::Complete { classification } => Ok(classification),
            TrimClassificationRun::Cancelled => Err(NurbsError::Domain {
                what: "non-cancelling trim classification observed cancellation".to_string(),
            }),
        }
    }

    fn classify_box_with_budget_and_poll(
        &self,
        min: [Rat; 2],
        max: [Rat; 2],
        work_remaining: &mut u128,
        should_cancel: &mut impl FnMut() -> bool,
    ) -> Result<TrimClassificationRun, NurbsError> {
        let (persistent_source_bytes, witness) = match prepare_trim_classification_with_poll(
            *self,
            min,
            max,
            work_remaining,
            should_cancel,
        )? {
            TrimWorkRun::Complete(prepared) => prepared,
            TrimWorkRun::Cancelled => return Ok(TrimClassificationRun::Cancelled),
        };
        let mut winding = 0i64;
        let mut operations_since_poll = 0usize;
        let query = TrimClassificationQuery {
            min,
            max,
            witness,
            max_depth: self.max_subdivision(),
        };
        for trim_loop in self.admitted_loops() {
            match loop_winding_box_with_poll(
                trim_loop.curve(),
                &query,
                persistent_source_bytes,
                work_remaining,
                should_cancel,
            )? {
                TrimWorkRun::Complete(Some(loop_winding)) => {
                    winding =
                        winding
                            .checked_add(loop_winding)
                            .ok_or_else(|| NurbsError::Domain {
                                what: "aggregate trim winding overflows i64".to_string(),
                            })?;
                }
                TrimWorkRun::Complete(None) => {
                    if should_cancel() {
                        return Ok(TrimClassificationRun::Cancelled);
                    }
                    return Ok(TrimClassificationRun::Complete {
                        classification: Classification::Boundary,
                    });
                }
                TrimWorkRun::Cancelled => return Ok(TrimClassificationRun::Cancelled),
            }
            if trim_poll_due(&mut operations_since_poll, should_cancel) {
                return Ok(TrimClassificationRun::Cancelled);
            }
        }
        let classification = if winding != 0 {
            Classification::Inside
        } else {
            Classification::Outside
        };
        if should_cancel() {
            return Ok(TrimClassificationRun::Cancelled);
        }
        Ok(TrimClassificationRun::Complete { classification })
    }
}

fn prepare_trim_classification_with_poll(
    patch: AdmittedTrimmedPatch<'_>,
    min: [Rat; 2],
    max: [Rat; 2],
    work_remaining: &mut u128,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<TrimWorkRun<(u128, [Rat; 2])>, NurbsError> {
    spend_trim_work(
        work_remaining,
        patch.loops().len() as u128,
        "persistent trim-source retained-byte accounting",
    )?;
    if should_cancel() {
        return Ok(TrimWorkRun::Cancelled);
    }
    let mut persistent_source_bytes = 0u128;
    let mut operations_since_poll = 0usize;
    for trim_loop in patch.admitted_loops() {
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
        if trim_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(TrimWorkRun::Cancelled);
        }
    }
    if should_cancel() {
        return Ok(TrimWorkRun::Cancelled);
    }
    enforce_trim_retained_bytes(persistent_source_bytes, "persistent source")?;
    spend_trim_work(
        work_remaining,
        TRIM_EXACT_MIDPOINT_WORK_UNITS * 2,
        "classification witness midpoints",
    )?;
    if should_cancel() {
        return Ok(TrimWorkRun::Cancelled);
    }
    let witness_u = exact_midpoint(min[0], max[0], "classification witness u")?;
    if should_cancel() {
        return Ok(TrimWorkRun::Cancelled);
    }
    let witness_v = exact_midpoint(min[1], max[1], "classification witness v")?;
    if should_cancel() {
        return Ok(TrimWorkRun::Cancelled);
    }
    Ok(TrimWorkRun::Complete((
        persistent_source_bytes,
        [witness_u, witness_v],
    )))
}

fn validate_classification_box(min: [Rat; 2], max: [Rat; 2]) -> Result<(), NurbsError> {
    if min[0] > max[0] || min[1] > max[1] {
        return Err(NurbsError::Domain {
            what: "trim classification box must be componentwise ordered".to_string(),
        });
    }
    Ok(())
}

fn trim_loop_validation_work_for(
    knot_count: usize,
    control_count: usize,
    degree: usize,
) -> Result<u128, NurbsError> {
    let control_components =
        (control_count as u128)
            .checked_mul(4)
            .ok_or_else(|| NurbsError::Domain {
                what: "trim control-validation accounting overflows u128".to_string(),
            })?;
    let order = (degree as u128)
        .checked_add(1)
        .ok_or_else(|| NurbsError::Domain {
            what: "trim order-validation accounting overflows u128".to_string(),
        })?;
    let basis_triangle = order.checked_mul(order).ok_or_else(|| NurbsError::Domain {
        what: "trim basis-validation accounting overflows u128".to_string(),
    })?;
    let scanned_entries = (knot_count as u128)
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

fn trim_loop_validation_work(curve: &NurbsCurve<Rat, 2>) -> Result<u128, NurbsError> {
    trim_loop_validation_work_for(curve.knots.knots.len(), curve.cpw.len(), curve.knots.degree)
}

fn preflight_trim_loop_reversal(
    trim_loop: AdmittedTrimLoop<'_>,
) -> Result<TrimLoopReversalPlan, NurbsError> {
    let curve = trim_loop.curve();
    preflight_trim_loop_reversal_counts(
        curve.knots().knots().len(),
        curve.homogeneous_control_points().len(),
        curve.knots().degree(),
    )
}

fn preflight_trim_loop_reversal_counts(
    knot_count: usize,
    control_count: usize,
    degree: usize,
) -> Result<TrimLoopReversalPlan, NurbsError> {
    let copy_work = (knot_count as u128)
        .checked_mul(16)
        .and_then(|work| work.checked_add((control_count as u128).checked_mul(4)?))
        .ok_or_else(|| NurbsError::Domain {
            what: "trim-loop reversal copy-work accounting overflows u128".to_string(),
        })?;
    let work_units = copy_work
        .checked_add(trim_loop_validation_work_for(
            knot_count,
            control_count,
            degree,
        )?)
        .and_then(|work| work.checked_add(64))
        .ok_or_else(|| NurbsError::Domain {
            what: "trim-loop reversal aggregate work overflows u128".to_string(),
        })?;
    if work_units > TRIM_CLASSIFY_MAX_WORK_UNITS {
        return Err(NurbsError::Domain {
            what: format!(
                "trim-loop reversal requests {work_units} work units above defensive ceiling {TRIM_CLASSIFY_MAX_WORK_UNITS}"
            ),
        });
    }

    let order = degree.checked_add(1).ok_or_else(|| NurbsError::Domain {
        what: "trim-loop reversal basis order overflows usize".to_string(),
    })?;
    let basis_workspace = (order as u128)
        .checked_mul(3)
        .and_then(|count| count.checked_mul(core::mem::size_of::<Rat>() as u128))
        .ok_or_else(|| NurbsError::Domain {
            what: "trim-loop reversal basis-workspace accounting overflows u128".to_string(),
        })?;
    let peak_retained_bytes = trim_curve_storage_bytes(knot_count, control_count)?
        .checked_add(basis_workspace)
        .ok_or_else(|| NurbsError::Domain {
            what: "trim-loop reversal retained-byte accounting overflows u128".to_string(),
        })?;
    enforce_trim_retained_bytes(peak_retained_bytes, "loop reversal")?;
    Ok(TrimLoopReversalPlan {
        knot_count,
        control_count,
        degree,
    })
}

fn assemble_reversed_trim_loop_with_poll(
    trim_loop: AdmittedTrimLoop<'_>,
    plan: TrimLoopReversalPlan,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<TrimWorkRun<TrimLoop>, NurbsError> {
    let curve = trim_loop.curve();
    let admitted_knots = curve.knots();
    let (lo, hi) = admitted_knots.domain();
    if should_cancel() {
        return Ok(TrimWorkRun::Cancelled);
    }
    let mut knots = Vec::new();
    knots
        .try_reserve_exact(plan.knot_count)
        .map_err(|_| NurbsError::Domain {
            what: "reversed trim-knot allocation was refused".to_string(),
        })?;
    let mut operations_since_poll = 0usize;
    for &knot in admitted_knots.knots().iter().rev() {
        knots.push(lo + (hi - knot));
        if trim_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(TrimWorkRun::Cancelled);
        }
    }
    if should_cancel() {
        return Ok(TrimWorkRun::Cancelled);
    }

    let controls = curve.homogeneous_control_points();
    let mut cpw = Vec::new();
    cpw.try_reserve_exact(plan.control_count)
        .map_err(|_| NurbsError::Domain {
            what: "reversed trim-control allocation was refused".to_string(),
        })?;
    operations_since_poll = 0;
    for &control in controls.iter().rev() {
        cpw.push(control);
        if trim_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(TrimWorkRun::Cancelled);
        }
    }
    if knots.len() != plan.knot_count || cpw.len() != plan.control_count {
        return Err(NurbsError::Structure {
            what: "trim-loop reversal assembly disagrees with its admitted shape".to_string(),
        });
    }
    if should_cancel() {
        return Ok(TrimWorkRun::Cancelled);
    }
    Ok(TrimWorkRun::Complete(TrimLoop {
        curve: NurbsCurve {
            knots: crate::basis::KnotVector {
                knots,
                degree: plan.degree,
            },
            cpw,
        },
    }))
}

fn publish_reversed_trim_loop_with_poll(
    trim_loop: TrimLoop,
    should_cancel: &mut impl FnMut() -> bool,
) -> TrimLoopReversalRun {
    if should_cancel() {
        return TrimLoopReversalRun::Cancelled;
    }
    TrimLoopReversalRun::Complete { trim_loop }
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

#[cfg(test)]
fn trim_bezier_conversion_plan(
    curve: AdmittedNurbsCurve<'_, Rat, 2>,
    persistent_source_bytes: u128,
    operation_source_is_persistent: bool,
) -> Result<BezierConversionPlan, NurbsError> {
    let plan = curve.bezier_conversion_plan()?;
    trim_bezier_conversion_plan_from_plan(
        curve,
        plan,
        persistent_source_bytes,
        operation_source_is_persistent,
    )
}

fn trim_bezier_conversion_plan_from_plan(
    curve: AdmittedNurbsCurve<'_, Rat, 2>,
    plan: BezierConversionPlan,
    persistent_source_bytes: u128,
    operation_source_is_persistent: bool,
) -> Result<BezierConversionPlan, NurbsError> {
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

fn trim_bezier_conversion_plan_with_poll(
    curve: AdmittedNurbsCurve<'_, Rat, 2>,
    persistent_source_bytes: u128,
    operation_source_is_persistent: bool,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<TrimWorkRun<BezierConversionPlan>, NurbsError> {
    let Some(plan) = curve.bezier_conversion_plan_with_poll(should_cancel)? else {
        return Ok(TrimWorkRun::Cancelled);
    };
    let plan = trim_bezier_conversion_plan_from_plan(
        curve,
        plan,
        persistent_source_bytes,
        operation_source_is_persistent,
    )?;
    if should_cancel() {
        return Ok(TrimWorkRun::Cancelled);
    }
    Ok(TrimWorkRun::Complete(plan))
}

fn preflight_trim_bezier_conversion_with_poll(
    curve: AdmittedNurbsCurve<'_, Rat, 2>,
    persistent_source_bytes: u128,
    operation_source_is_persistent: bool,
    work_remaining: &mut u128,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<TrimWorkRun<(BezierConversionPlan, u128)>, NurbsError> {
    let pre_scan_work = curve.bezier_pre_scan_work()?;
    spend_trim_work(
        work_remaining,
        pre_scan_work,
        "Bezier conversion-plan knot scan",
    )?;
    let Some(plan) = curve.bezier_conversion_plan_with_poll(should_cancel)? else {
        return Ok(TrimWorkRun::Cancelled);
    };
    let plan = trim_bezier_conversion_plan_from_plan(
        curve,
        plan,
        persistent_source_bytes,
        operation_source_is_persistent,
    )?;
    let unspent_work =
        plan.work_units
            .checked_sub(pre_scan_work)
            .ok_or_else(|| NurbsError::Domain {
                what: "trim Bezier post-scan work accounting is inconsistent".to_string(),
            })?;
    if should_cancel() {
        return Ok(TrimWorkRun::Cancelled);
    }
    Ok(TrimWorkRun::Complete((plan, unspent_work)))
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

fn collect_offending_intervals_with_poll(
    boxes: &[SpanBox<Rat, 2>],
    query_min: [Rat; 2],
    query_max: [Rat; 2],
    span_capacity: usize,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<TrimWorkRun<Vec<(Rat, Rat)>>, NurbsError> {
    if should_cancel() {
        return Ok(TrimWorkRun::Cancelled);
    }
    let mut offending = Vec::new();
    offending
        .try_reserve_exact(span_capacity)
        .map_err(|_| NurbsError::Domain {
            what: "trim offending-interval allocation was refused".to_string(),
        })?;
    if should_cancel() {
        return Ok(TrimWorkRun::Cancelled);
    }
    let mut operations_since_poll = 0usize;
    for &(min, max, t0, t1) in boxes {
        if max[0] >= query_min[0]
            && min[0] <= query_max[0]
            && max[1] >= query_min[1]
            && min[1] <= query_max[1]
        {
            offending.push((t0, t1));
        }
        if trim_poll_due(&mut operations_since_poll, should_cancel) {
            return Ok(TrimWorkRun::Cancelled);
        }
    }
    if should_cancel() {
        return Ok(TrimWorkRun::Cancelled);
    }
    Ok(TrimWorkRun::Complete(offending))
}

fn subdivide_offending_intervals_with_poll(
    mut work: NurbsCurve<Rat, 2>,
    boxes: Vec<SpanBox<Rat, 2>>,
    offending: Vec<(Rat, Rat)>,
    persistent_source_bytes: u128,
    work_remaining: &mut u128,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<TrimWorkRun<NurbsCurve<Rat, 2>>, NurbsError> {
    let admitted_work = work.admitted_after_validation();
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
    if should_cancel() {
        return Ok(TrimWorkRun::Cancelled);
    }
    for (t0, t1) in offending {
        let admitted_work = work.admitted_after_validation();
        if should_cancel() {
            return Ok(TrimWorkRun::Cancelled);
        }
        let mid = exact_midpoint(t0, t1, "subdivision parameter")?;
        if should_cancel() {
            return Ok(TrimWorkRun::Cancelled);
        }
        // Exact midpoint insertion splits the offending span.
        work = match admitted_work.insert_knot_with_poll(mid, should_cancel)? {
            crate::curve::CurveInsertionRun::Complete { curve } => curve,
            crate::curve::CurveInsertionRun::Cancelled => return Ok(TrimWorkRun::Cancelled),
        };
        // Dropping the superseded generation is not preemptible; gate the
        // next interval immediately after replacement.
        if should_cancel() {
            return Ok(TrimWorkRun::Cancelled);
        }
    }
    let admitted_work = work.admitted_after_validation();
    // The aggregate refinement charge was consumed before the first midpoint
    // insertion and includes this external plan scan.
    let conversion_plan = match trim_bezier_conversion_plan_with_poll(
        admitted_work,
        persistent_source_bytes,
        false,
        should_cancel,
    )? {
        TrimWorkRun::Complete(plan) => plan,
        TrimWorkRun::Cancelled => return Ok(TrimWorkRun::Cancelled),
    };
    if conversion_plan.insertions != expected_conversion_insertions {
        return Err(NurbsError::Structure {
            what: format!(
                "trim midpoint refinement projected {expected_conversion_insertions} Bezier insertions but derived generation requires {}",
                conversion_plan.insertions
            ),
        });
    }
    work = match admitted_work.to_bezier_form_with_poll(should_cancel)? {
        crate::curve::CurveBezierRun::Complete { curve } => curve,
        crate::curve::CurveBezierRun::Cancelled => return Ok(TrimWorkRun::Cancelled),
    };
    // As above, derived-generation destruction is nonpreemptible.
    if should_cancel() {
        return Ok(TrimWorkRun::Cancelled);
    }
    Ok(TrimWorkRun::Complete(work))
}

fn initialize_trim_bezier_work_with_poll(
    curve: AdmittedNurbsCurve<'_, Rat, 2>,
    persistent_source_bytes: u128,
    work_remaining: &mut u128,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<TrimWorkRun<NurbsCurve<Rat, 2>>, NurbsError> {
    let (initial_plan, initial_conversion_work) = match preflight_trim_bezier_conversion_with_poll(
        curve,
        persistent_source_bytes,
        true,
        work_remaining,
        should_cancel,
    )? {
        TrimWorkRun::Complete(plan) => plan,
        TrimWorkRun::Cancelled => return Ok(TrimWorkRun::Cancelled),
    };
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
    if should_cancel() {
        return Ok(TrimWorkRun::Cancelled);
    }
    let work = match curve.to_bezier_form_with_poll(should_cancel)? {
        crate::curve::CurveBezierRun::Complete { curve } => curve,
        crate::curve::CurveBezierRun::Cancelled => return Ok(TrimWorkRun::Cancelled),
    };
    if should_cancel() {
        return Ok(TrimWorkRun::Cancelled);
    }
    Ok(TrimWorkRun::Complete(work))
}

/// Certified winding number of one closed rational curve about `q`, or
/// `None` when `q` cannot be separated from the curve within the
/// subdivision budget. Cancellation drops every partial derived generation,
/// span table, and offending-interval table before returning.
fn loop_winding_box_with_poll(
    curve: AdmittedNurbsCurve<'_, Rat, 2>,
    query: &TrimClassificationQuery,
    persistent_source_bytes: u128,
    work_remaining: &mut u128,
    should_cancel: &mut impl FnMut() -> bool,
) -> Result<TrimWorkRun<Option<i64>>, NurbsError> {
    // Work in Bézier form so each span's control hull tightly bounds it.
    let mut work = match initialize_trim_bezier_work_with_poll(
        curve,
        persistent_source_bytes,
        work_remaining,
        should_cancel,
    )? {
        TrimWorkRun::Complete(work) => work,
        TrimWorkRun::Cancelled => return Ok(TrimWorkRun::Cancelled),
    };
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
        let boxes = match admitted_work.span_boxes_with_poll(should_cancel)? {
            crate::curve::CurveSpanBoxesRun::Complete { boxes } => boxes,
            crate::curve::CurveSpanBoxesRun::Cancelled => return Ok(TrimWorkRun::Cancelled),
        };
        let offending = match collect_offending_intervals_with_poll(
            &boxes,
            query.min,
            query.max,
            span_capacity,
            should_cancel,
        )? {
            TrimWorkRun::Complete(offending) => offending,
            TrimWorkRun::Cancelled => return Ok(TrimWorkRun::Cancelled),
        };
        if offending.is_empty() {
            // Separated from the whole connected query box: winding is
            // constant throughout it, so one exact witness is sufficient.
            let winding = match polygon_winding_homogeneous_with_poll(
                admitted_work.homogeneous_control_points(),
                query.witness,
                should_cancel,
            ) {
                TrimWorkRun::Complete(winding) => winding,
                TrimWorkRun::Cancelled => return Ok(TrimWorkRun::Cancelled),
            };
            drop(offending);
            drop(boxes);
            drop(work);
            if should_cancel() {
                return Ok(TrimWorkRun::Cancelled);
            }
            return Ok(TrimWorkRun::Complete(Some(winding)));
        }
        if depth >= query.max_depth {
            drop(offending);
            drop(boxes);
            drop(work);
            if should_cancel() {
                return Ok(TrimWorkRun::Cancelled);
            }
            return Ok(TrimWorkRun::Complete(None));
        }
        work = match subdivide_offending_intervals_with_poll(
            work,
            boxes,
            offending,
            persistent_source_bytes,
            work_remaining,
            should_cancel,
        )? {
            TrimWorkRun::Complete(work) => work,
            TrimWorkRun::Cancelled => return Ok(TrimWorkRun::Cancelled),
        };
        depth = depth.checked_add(1).ok_or_else(|| NurbsError::Domain {
            what: "trim subdivision depth overflows u32".to_string(),
        })?;
        if should_cancel() {
            return Ok(TrimWorkRun::Cancelled);
        }
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
/// control. Individual exact edge arithmetic is nonpreemptible.
fn polygon_winding_homogeneous_with_poll(
    cpw: &[[Rat; 4]],
    q: [Rat; 2],
    should_cancel: &mut impl FnMut() -> bool,
) -> TrimWorkRun<i64> {
    if should_cancel() {
        return TrimWorkRun::Cancelled;
    }
    let mut winding = 0i64;
    let mut operations_since_poll = 0usize;
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
        if trim_poll_due(&mut operations_since_poll, should_cancel) {
            return TrimWorkRun::Cancelled;
        }
    }
    if should_cancel() {
        return TrimWorkRun::Cancelled;
    }
    TrimWorkRun::Complete(winding)
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

    fn poly_trim_loop(points: &[[i64; 2]]) -> TrimLoop {
        let segment_count = points.len();
        let mut knots = vec![Rat::int(0), Rat::int(0)];
        let denominator = i128::try_from(segment_count).expect("test loop length fits i128");
        knots.extend((1..segment_count).map(|index| {
            Rat::new(
                i128::try_from(index).expect("test knot index fits i128"),
                denominator,
            )
        }));
        knots.extend([Rat::int(1), Rat::int(1)]);
        let knots = KnotVector::new(knots, 1).expect("polyline loop knots");
        let mut controls: Vec<_> = points
            .iter()
            .map(|point| [Rat::int(point[0]), Rat::int(point[1])])
            .collect();
        controls.push(controls[0]);
        let weights = vec![Rat::int(1); controls.len()];
        let curve = NurbsCurve::new(knots, &controls, &weights).expect("polyline loop curve");
        TrimLoop::new(curve).expect("closed polyline loop")
    }

    fn assemble_reversed_for_test(trim_loop: AdmittedTrimLoop<'_>) -> TrimLoop {
        let plan = preflight_trim_loop_reversal(trim_loop).expect("reversal preflight");
        let mut never_cancel = || false;
        let TrimWorkRun::Complete(candidate) =
            assemble_reversed_trim_loop_with_poll(trim_loop, plan, &mut never_cancel)
                .expect("reversal assembly")
        else {
            panic!("healthy reversal assembly must complete");
        };
        candidate
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
    fn trim_loop_construction_with_cx_is_transactional_and_exact() {
        let expected = point_trim_loop();
        with_trim_cx(true, |cx| {
            assert_eq!(
                TrimLoop::new_with_cx(expected.curve.try_clone().expect("curve copy"), cx,)
                    .expect("valid pre-cancelled loop construction"),
                TrimLoopConstructionRun::Cancelled
            );
        });
        with_trim_cx(false, |cx| {
            assert_eq!(
                TrimLoop::new_with_cx(expected.curve.try_clone().expect("curve copy"), cx,)
                    .expect("active loop construction"),
                TrimLoopConstructionRun::Complete {
                    trim_loop: expected.try_clone().expect("expected loop copy"),
                }
            );
        });

        let open_curve = || {
            NurbsCurve::new(
                KnotVector::new(vec![Rat::int(0), Rat::int(0), Rat::int(1), Rat::int(1)], 1)
                    .expect("open-loop knots"),
                &[[Rat::int(0), Rat::int(0)], [Rat::int(1), Rat::int(0)]],
                &[Rat::int(1); 2],
            )
            .expect("open curve")
        };
        let legacy_error = TrimLoop::new(open_curve()).expect_err("legacy open-loop refusal");
        with_trim_cx(false, |cx| {
            assert_eq!(
                TrimLoop::new_with_cx(open_curve(), cx).expect_err("cancellable open-loop refusal"),
                legacy_error
            );
        });
    }

    #[test]
    fn trim_loop_construction_propagates_validation_and_publication_cancellation() {
        let candidate = point_trim_loop();
        let mut unexpected_publication_polls = 0usize;
        assert_eq!(
            finish_trim_loop_construction_with_poll(
                candidate,
                TrimLoopValidationOutcome::Cancelled,
                &mut || {
                    unexpected_publication_polls += 1;
                    false
                },
            ),
            TrimLoopConstructionRun::Cancelled
        );
        assert_eq!(unexpected_publication_polls, 0);

        let candidate = point_trim_loop();
        let mut publication_polls = 0usize;
        assert_eq!(
            finish_trim_loop_construction_with_poll(
                candidate,
                TrimLoopValidationOutcome::Complete,
                &mut || {
                    publication_polls += 1;
                    true
                },
            ),
            TrimLoopConstructionRun::Cancelled
        );
        assert_eq!(publication_polls, 1);
    }

    #[test]
    fn trim_loop_copy_with_cx_is_transactional_and_exact() {
        let trim_loop = long_trim_loop();
        with_trim_cx(true, |cx| {
            assert_eq!(
                trim_loop
                    .try_clone_with_cx(cx)
                    .expect("valid pre-cancelled trim-loop copy"),
                TrimLoopCloneRun::Cancelled
            );
        });
        with_trim_cx(false, |cx| {
            assert_eq!(
                trim_loop
                    .try_clone_with_cx(cx)
                    .expect("active exact trim-loop copy"),
                TrimLoopCloneRun::Complete {
                    trim_loop: trim_loop.try_clone().expect("legacy trim-loop copy"),
                }
            );
        });
    }

    #[test]
    fn trim_loop_copy_propagates_nested_and_publication_cancellation() {
        let mut unexpected_publication_polls = 0usize;
        assert_eq!(
            finish_trim_loop_clone_with_poll(CurveCloneRun::Cancelled, &mut || {
                unexpected_publication_polls += 1;
                false
            },),
            TrimLoopCloneRun::Cancelled
        );
        assert_eq!(unexpected_publication_polls, 0);

        let trim_loop = point_trim_loop();
        let curve = trim_loop.curve.try_clone().expect("exact curve copy");
        let mut publication_polls = 0usize;
        assert_eq!(
            publish_trim_loop_clone_with_poll(curve, &mut || {
                publication_polls += 1;
                true
            }),
            TrimLoopCloneRun::Cancelled
        );
        assert_eq!(publication_polls, 1);
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
    fn trim_loop_reversal_with_cx_is_transactional_and_exact() {
        let trim_loop = poly_trim_loop(&[[0, 0], [2, 0], [2, 2], [0, 2]]);
        let admitted = trim_loop.admit().expect("admitted square loop");
        let legacy = admitted.reversed_for_hole().expect("legacy reversal");
        assert_eq!(
            legacy
                .admit()
                .expect("admitted reversed loop")
                .reversed_for_hole()
                .expect("double reversal"),
            trim_loop
        );

        with_trim_cx(true, |cx| {
            assert_eq!(
                admitted
                    .reversed_for_hole_with_cx(cx)
                    .expect("valid pre-cancelled reversal"),
                TrimLoopReversalRun::Cancelled
            );
        });
        with_trim_cx(false, |cx| {
            assert_eq!(
                admitted
                    .reversed_for_hole_with_cx(cx)
                    .expect("active reversal"),
                TrimLoopReversalRun::Complete {
                    trim_loop: legacy.try_clone().expect("expected reversal copy"),
                }
            );
        });

        let query = [Rat::int(1), Rat::int(1)];
        let mut never_cancel = || false;
        let TrimWorkRun::Complete(source_winding) = polygon_winding_homogeneous_with_poll(
            admitted.curve().homogeneous_control_points(),
            query,
            &mut never_cancel,
        ) else {
            panic!("healthy source winding must complete");
        };
        let reversed = legacy.admit().expect("admitted reversed loop");
        let TrimWorkRun::Complete(reversed_winding) = polygon_winding_homogeneous_with_poll(
            reversed.curve().homogeneous_control_points(),
            query,
            &mut never_cancel,
        ) else {
            panic!("healthy reversed winding must complete");
        };
        assert_eq!(reversed_winding, -source_winding);
        assert_ne!(source_winding, 0);
    }

    #[test]
    fn trim_loop_reversal_preserves_large_same_sign_domain_with_cx() {
        let lo = Rat::new(i128::MAX - 10, 1);
        let hi = Rat::new(i128::MAX - 1, 1);
        let trim_loop = TrimLoop::new(
            NurbsCurve::new(
                KnotVector::new(vec![lo, lo, hi, hi], 1).expect("large same-sign knots"),
                &[[Rat::int(0), Rat::int(0)]; 2],
                &[Rat::int(1); 2],
            )
            .expect("large-domain curve"),
        )
        .expect("large-domain loop");
        let admitted = trim_loop.admit().expect("admitted large-domain loop");
        with_trim_cx(false, |cx| {
            assert!(matches!(
                admitted
                    .reversed_for_hole_with_cx(cx)
                    .expect("active large-domain reversal"),
                TrimLoopReversalRun::Complete { .. }
            ));
        });
    }

    #[test]
    fn trim_loop_reversal_cancels_inside_both_copy_phases() {
        let trim_loop = long_trim_loop();
        let admitted = trim_loop.admit().expect("admitted long loop");
        let plan = preflight_trim_loop_reversal(admitted).expect("reversal preflight");
        let run = |target| {
            let mut polls = 0usize;
            let mut should_cancel = || {
                polls += 1;
                polls == target
            };
            let outcome = assemble_reversed_trim_loop_with_poll(admitted, plan, &mut should_cancel)
                .expect("bounded reversal assembly");
            (matches!(outcome, TrimWorkRun::Cancelled), polls)
        };
        assert_eq!(run(2), run(2));
        assert_eq!(run(2), (true, 2));
        assert_eq!(run(5), run(5));
        assert_eq!(run(5), (true, 5));
    }

    #[test]
    fn trim_loop_reversal_propagates_nested_and_publication_cancellation() {
        let trim_loop = point_trim_loop();
        let admitted = trim_loop.admit().expect("admitted point loop");
        let candidate = assemble_reversed_for_test(admitted);
        with_trim_cx(true, |cx| {
            assert!(matches!(
                candidate
                    .admit_with_cx(cx)
                    .expect("pre-cancelled derived admission"),
                TrimLoopAdmissionRun::Cancelled
            ));
        });

        let candidate = assemble_reversed_for_test(admitted);
        candidate.admit().expect("validated reversal candidate");
        let mut polls = 0usize;
        let mut cancel_at_publication = || {
            polls += 1;
            polls == 1
        };
        assert_eq!(
            publish_reversed_trim_loop_with_poll(candidate, &mut cancel_at_publication),
            TrimLoopReversalRun::Cancelled
        );
        assert_eq!(polls, 1);
    }

    #[test]
    fn trim_loop_reversal_refuses_work_before_retained_bytes() {
        let error = preflight_trim_loop_reversal_counts(usize::MAX, usize::MAX, 1)
            .expect_err("work must refuse before retained-byte accounting");
        assert!(matches!(
            error,
            NurbsError::Domain { ref what } if what.contains("work units above defensive ceiling")
        ));
    }

    #[test]
    fn trimmed_patch_admission_with_cx_is_transactional_and_lifetime_bound() {
        let patch = TrimmedPatch::with_max_subdivision(vec![point_trim_loop()], 7);
        with_trim_cx(true, |cx| {
            assert!(matches!(
                patch
                    .admit_with_cx(cx)
                    .expect("valid pre-cancelled trimmed patch"),
                TrimmedPatchAdmissionRun::Cancelled
            ));
        });
        with_trim_cx(false, |cx| {
            let TrimmedPatchAdmissionRun::Complete { admitted } = patch
                .admit_with_cx(cx)
                .expect("active trimmed-patch admission")
            else {
                panic!("active trimmed-patch admission must complete");
            };
            assert!(core::ptr::eq(admitted.source(), &patch));
            assert_eq!(admitted.loops().len(), 1);
            assert_eq!(admitted.max_subdivision(), 7);
        });
    }

    #[test]
    fn trimmed_patch_minimum_work_refusal_precedes_cancellation() {
        let patch = TrimmedPatch::new(vec![point_trim_loop()]);
        let mut work_remaining = 0u128;
        let mut polls = 0usize;
        let error = patch
            .admit_with_budget_and_poll(
                &mut work_remaining,
                &mut || {
                    polls += 1;
                    true
                },
                &mut |_trim_loop| -> Result<TrimmedPatchValidationOutcome, NurbsError> {
                    panic!("static minimum-work refusal must precede loop admission")
                },
            )
            .expect_err("static minimum-work refusal must precede cancellation");
        assert!(matches!(
            error,
            NurbsError::Domain { ref what } if what.contains("at least")
        ));
        assert_eq!(polls, 0);
    }

    #[test]
    fn trimmed_patch_plan_scan_cancels_at_a_replayable_stride() {
        let patch = TrimmedPatch::new((0..130).map(|_| point_trim_loop()).collect());
        let run = || {
            let mut work_remaining = TRIM_CLASSIFY_MAX_WORK_UNITS;
            let mut polls = 0usize;
            let mut admitted_loops = 0usize;
            let outcome = patch
                .admit_with_budget_and_poll(
                    &mut work_remaining,
                    &mut || {
                        polls += 1;
                        polls == 2
                    },
                    &mut |_trim_loop| {
                        admitted_loops += 1;
                        Ok(TrimmedPatchValidationOutcome::Complete)
                    },
                )
                .expect("cancellable trimmed-patch plan scan");
            (
                matches!(outcome, TrimmedPatchAdmissionRun::Cancelled),
                polls,
                admitted_loops,
                work_remaining,
            )
        };
        assert_eq!(run(), run());
        assert_eq!(run(), (true, 2, 0, TRIM_CLASSIFY_MAX_WORK_UNITS));
    }

    #[test]
    fn trimmed_patch_nested_cancellation_is_not_published() {
        let patch = TrimmedPatch::new(vec![point_trim_loop()]);
        let mut work_remaining = TRIM_CLASSIFY_MAX_WORK_UNITS;
        let mut admitted_loops = 0usize;
        let outcome = patch
            .admit_with_budget_and_poll(&mut work_remaining, &mut || false, &mut |_trim_loop| {
                admitted_loops += 1;
                Ok(TrimmedPatchValidationOutcome::Cancelled)
            })
            .expect("nested trim-loop cancellation");
        assert!(matches!(outcome, TrimmedPatchAdmissionRun::Cancelled));
        assert_eq!(admitted_loops, 1);
    }

    #[test]
    fn trimmed_patch_final_checkpoint_gates_authority_publication() {
        let patch = TrimmedPatch::new(Vec::new());
        let mut healthy_work = TRIM_CLASSIFY_MAX_WORK_UNITS;
        let mut total_polls = 0usize;
        let healthy = patch
            .admit_with_budget_and_poll(
                &mut healthy_work,
                &mut || {
                    total_polls += 1;
                    false
                },
                &mut |_trim_loop| -> Result<TrimmedPatchValidationOutcome, NurbsError> {
                    panic!("empty patch has no loop admission")
                },
            )
            .expect("healthy empty-patch admission");
        assert!(matches!(healthy, TrimmedPatchAdmissionRun::Complete { .. }));
        assert!(total_polls > 0);

        let mut replay_work = TRIM_CLASSIFY_MAX_WORK_UNITS;
        let mut replay_polls = 0usize;
        let replay = patch
            .admit_with_budget_and_poll(
                &mut replay_work,
                &mut || {
                    replay_polls += 1;
                    replay_polls == total_polls
                },
                &mut |_trim_loop| -> Result<TrimmedPatchValidationOutcome, NurbsError> {
                    panic!("empty patch has no loop admission")
                },
            )
            .expect("cancelled empty-patch admission replay");
        assert!(matches!(replay, TrimmedPatchAdmissionRun::Cancelled));
        assert_eq!(replay_polls, total_polls);
    }

    #[test]
    fn trim_classification_owning_admitted_and_cx_paths_match_exactly() {
        let patch = TrimmedPatch::with_max_subdivision(
            vec![poly_trim_loop(&[[0, 0], [10, 0], [10, 10], [0, 10]])],
            0,
        );
        let admitted = patch.admit().expect("admitted square trim patch");
        for (query, expected) in [
            ([Rat::int(5), Rat::int(5)], Classification::Inside),
            ([Rat::int(20), Rat::int(20)], Classification::Outside),
            ([Rat::int(0), Rat::int(5)], Classification::Boundary),
        ] {
            assert_eq!(
                patch.classify(query).expect("owning classification"),
                expected
            );
            assert_eq!(
                admitted.classify(query).expect("admitted classification"),
                expected
            );
            with_trim_cx(false, |cx| {
                assert_eq!(
                    patch
                        .classify_with_cx(query, cx)
                        .expect("owning cancellable classification"),
                    TrimClassificationRun::Complete {
                        classification: expected,
                    }
                );
                assert_eq!(
                    admitted
                        .classify_with_cx(query, cx)
                        .expect("admitted cancellable classification"),
                    TrimClassificationRun::Complete {
                        classification: expected,
                    }
                );
            });
        }
    }

    #[test]
    fn trim_classification_request_and_static_work_refusals_precede_cancellation() {
        let patch = TrimmedPatch::new(vec![point_trim_loop()]);
        let admitted = patch.admit().expect("admitted point trim patch");
        with_trim_cx(true, |cx| {
            assert_eq!(
                patch
                    .classify_with_cx([Rat::int(2), Rat::int(2)], cx)
                    .expect("valid pre-cancelled owning classification"),
                TrimClassificationRun::Cancelled
            );
            assert_eq!(
                admitted
                    .classify_with_cx([Rat::int(2), Rat::int(2)], cx)
                    .expect("valid pre-cancelled classification"),
                TrimClassificationRun::Cancelled
            );
            let synchronous = admitted
                .classify_box([Rat::int(1), Rat::int(0)], [Rat::int(0), Rat::int(1)])
                .expect_err("inverted synchronous box");
            let cancellable = admitted
                .classify_box_with_cx([Rat::int(1), Rat::int(0)], [Rat::int(0), Rat::int(1)], cx)
                .expect_err("inverted box must precede cancellation");
            assert_eq!(cancellable, synchronous);
            let owning = patch
                .classify_box_with_cx([Rat::int(1), Rat::int(0)], [Rat::int(0), Rat::int(1)], cx)
                .expect_err("owning inverted box must precede cancellation");
            assert_eq!(owning, synchronous);
        });

        let mut work_remaining = 0u128;
        let mut polls = 0usize;
        let error = admitted
            .classify_box_with_budget_and_poll(
                [Rat::int(2), Rat::int(2)],
                [Rat::int(2), Rat::int(2)],
                &mut work_remaining,
                &mut || {
                    polls += 1;
                    true
                },
            )
            .expect_err("constant loop-count work refusal must precede cancellation");
        assert!(matches!(
            error,
            NurbsError::Domain { ref what } if what.contains("persistent trim-source")
        ));
        assert_eq!(polls, 0);
    }

    #[test]
    fn trim_offending_and_winding_scans_cancel_at_replayable_strides() {
        let boxes: Vec<SpanBox<Rat, 2>> = (0..130)
            .map(|index| {
                let t0 = Rat::new(index, 130);
                let t1 = Rat::new(index + 1, 130);
                ([Rat::int(0); 2], [Rat::int(1); 2], t0, t1)
            })
            .collect();
        let run_boxes = || {
            let mut polls = 0usize;
            let outcome = collect_offending_intervals_with_poll(
                &boxes,
                [Rat::int(2); 2],
                [Rat::int(2); 2],
                boxes.len(),
                &mut || {
                    polls += 1;
                    polls == 3
                },
            )
            .expect("cancellable offending scan");
            (matches!(outcome, TrimWorkRun::Cancelled), polls)
        };
        assert_eq!(run_boxes(), run_boxes());
        assert_eq!(run_boxes(), (true, 3));

        let controls = vec![[Rat::int(0), Rat::int(0), Rat::int(0), Rat::int(1)]; 130];
        let run_winding = || {
            let mut polls = 0usize;
            let outcome = polygon_winding_homogeneous_with_poll(
                &controls,
                [Rat::int(2), Rat::int(2)],
                &mut || {
                    polls += 1;
                    polls == 2
                },
            );
            (matches!(outcome, TrimWorkRun::Cancelled), polls)
        };
        assert_eq!(run_winding(), run_winding());
        assert_eq!(run_winding(), (true, 2));
    }

    #[test]
    fn trim_conversion_plan_cancellation_preserves_consumed_pre_scan_charge() {
        let trim_loop = long_trim_loop();
        let curve = trim_loop.curve.admit().expect("admitted long trim curve");
        let persistent_source_bytes = trim_curve_storage_bytes(
            curve.knots().knots().len(),
            curve.homogeneous_control_points().len(),
        )
        .expect("persistent trim bytes");
        let pre_scan_work = curve.bezier_pre_scan_work().expect("pre-scan work");
        let run = || {
            let mut work_remaining = TRIM_CLASSIFY_MAX_WORK_UNITS;
            let mut polls = 0usize;
            let outcome = preflight_trim_bezier_conversion_with_poll(
                curve,
                persistent_source_bytes,
                true,
                &mut work_remaining,
                &mut || {
                    polls += 1;
                    polls == 2
                },
            )
            .expect("cancellable trim conversion plan");
            (
                matches!(outcome, TrimWorkRun::Cancelled),
                polls,
                work_remaining,
            )
        };
        assert_eq!(run(), run());
        assert_eq!(
            run(),
            (true, 2, TRIM_CLASSIFY_MAX_WORK_UNITS - pre_scan_work)
        );
    }

    #[test]
    fn trim_loop_terminal_cleanup_and_outer_publication_replay_exactly() {
        let trim_loop = point_trim_loop();
        let curve = trim_loop.curve.admit().expect("admitted point trim curve");
        let persistent_source_bytes = trim_curve_storage_bytes(
            curve.knots().knots().len(),
            curve.homogeneous_control_points().len(),
        )
        .expect("point trim bytes");
        for (query, expected) in [
            ([Rat::int(2), Rat::int(2)], Some(0i64)),
            ([Rat::int(0), Rat::int(0)], None),
        ] {
            let request = TrimClassificationQuery {
                min: query,
                max: query,
                witness: query,
                max_depth: 0,
            };
            let mut healthy_work = TRIM_CLASSIFY_MAX_WORK_UNITS;
            let mut total_polls = 0usize;
            let healthy = loop_winding_box_with_poll(
                curve,
                &request,
                persistent_source_bytes,
                &mut healthy_work,
                &mut || {
                    total_polls += 1;
                    false
                },
            )
            .expect("healthy trim loop terminal path");
            assert_eq!(healthy, TrimWorkRun::Complete(expected));
            assert!(total_polls > 0);

            let mut replay_work = TRIM_CLASSIFY_MAX_WORK_UNITS;
            let mut replay_polls = 0usize;
            let replay = loop_winding_box_with_poll(
                curve,
                &request,
                persistent_source_bytes,
                &mut replay_work,
                &mut || {
                    replay_polls += 1;
                    replay_polls == total_polls
                },
            )
            .expect("cancelled trim loop terminal replay");
            assert_eq!(replay, TrimWorkRun::Cancelled);
            assert_eq!(replay_polls, total_polls);
        }

        let patch = TrimmedPatch::new(Vec::new());
        let admitted = patch.admit().expect("admitted empty patch");
        let mut healthy_work = TRIM_CLASSIFY_MAX_WORK_UNITS;
        let mut total_polls = 0usize;
        let healthy = admitted
            .classify_box_with_budget_and_poll(
                [Rat::int(0); 2],
                [Rat::int(0); 2],
                &mut healthy_work,
                &mut || {
                    total_polls += 1;
                    false
                },
            )
            .expect("healthy empty classification");
        assert_eq!(
            healthy,
            TrimClassificationRun::Complete {
                classification: Classification::Outside,
            }
        );
        let mut replay_work = TRIM_CLASSIFY_MAX_WORK_UNITS;
        let mut replay_polls = 0usize;
        let replay = admitted
            .classify_box_with_budget_and_poll(
                [Rat::int(0); 2],
                [Rat::int(0); 2],
                &mut replay_work,
                &mut || {
                    replay_polls += 1;
                    replay_polls == total_polls
                },
            )
            .expect("cancelled empty classification replay");
        assert_eq!(replay, TrimClassificationRun::Cancelled);
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
    fn trimmed_patch_copy_with_cx_is_transactional_and_exact() {
        let patch = TrimmedPatch::with_max_subdivision(vec![point_trim_loop()], 7);
        with_trim_cx(true, |cx| {
            assert_eq!(
                patch
                    .try_clone_with_cx(cx)
                    .expect("valid pre-cancelled patch copy"),
                TrimmedPatchCloneRun::Cancelled
            );
        });
        with_trim_cx(false, |cx| {
            assert_eq!(
                patch
                    .try_clone_with_cx(cx)
                    .expect("active exact patch copy"),
                TrimmedPatchCloneRun::Complete {
                    trimmed_patch: patch.try_clone().expect("legacy patch copy"),
                }
            );
        });
    }

    #[test]
    fn trimmed_patch_copy_plan_is_aggregate_and_cancellable() {
        let patch = TrimmedPatch::new(vec![point_trim_loop()]);
        let mut never_cancel = || false;
        let TrimWorkRun::Complete(plan) =
            preflight_trimmed_patch_copy_with_poll(&patch, &mut never_cancel)
                .expect("healthy patch-copy plan")
        else {
            panic!("active plan must complete");
        };
        assert_eq!(plan.loop_count, 1);
        assert_eq!(plan.knot_count, 4);
        assert_eq!(plan.control_count, 2);
        assert_eq!(plan.work_units, 4 + 4 * 2 + 4 + 2);
        assert_eq!(
            plan.retained_bytes,
            core::mem::size_of::<TrimLoop>() as u128
                + 4 * core::mem::size_of::<Rat>() as u128
                + 2 * core::mem::size_of::<[Rat; 4]>() as u128
        );

        let many = TrimmedPatch::new((0..130).map(|_| point_trim_loop()).collect());
        let run = || {
            let mut polls = 0usize;
            let mut should_cancel = || {
                polls += 1;
                polls == 2
            };
            let outcome = preflight_trimmed_patch_copy_with_poll(&many, &mut should_cancel)
                .expect("cancellable metadata plan");
            (matches!(outcome, TrimWorkRun::Cancelled), polls)
        };
        assert_eq!(run(), run());
        assert_eq!(run(), (true, 2));

        let mut polls = 0usize;
        let mut copied_loops = 0usize;
        let move_outcome = many
            .try_clone_with_nested_and_poll(
                &mut || {
                    polls += 1;
                    polls == 6
                },
                &mut |trim_loop| {
                    copied_loops += 1;
                    Ok(TrimLoopCloneRun::Complete {
                        trim_loop: trim_loop.try_clone()?,
                    })
                },
            )
            .expect("cancellable outer table moves");
        assert_eq!(move_outcome, TrimmedPatchCloneRun::Cancelled);
        assert_eq!(polls, 6);
        assert_eq!(copied_loops, 64);
    }

    #[test]
    fn trimmed_patch_copy_propagates_nested_and_publication_cancellation() {
        let patch = TrimmedPatch::new(vec![point_trim_loop(), point_trim_loop()]);
        let mut nested_calls = 0usize;
        let outcome = patch
            .try_clone_with_nested_and_poll(&mut || false, &mut |_trim_loop| {
                nested_calls += 1;
                if nested_calls == 1 {
                    Ok(TrimLoopCloneRun::Cancelled)
                } else {
                    panic!("nested cancellation must stop later loop copies")
                }
            })
            .expect("nested patch-copy cancellation");
        assert_eq!(outcome, TrimmedPatchCloneRun::Cancelled);
        assert_eq!(nested_calls, 1);

        let empty = TrimmedPatch::with_max_subdivision(Vec::new(), 7);
        let mut total_polls = 0usize;
        let mut never_cancel = || {
            total_polls += 1;
            false
        };
        let mut clone_no_loop = |_trim_loop: &TrimLoop| -> Result<TrimLoopCloneRun, NurbsError> {
            panic!("empty patch has no nested copy")
        };
        assert!(matches!(
            empty
                .try_clone_with_nested_and_poll(&mut never_cancel, &mut clone_no_loop)
                .expect("healthy empty-patch copy"),
            TrimmedPatchCloneRun::Complete { .. }
        ));
        let mut replay_polls = 0usize;
        let mut cancel_at_publication = || {
            replay_polls += 1;
            replay_polls == total_polls
        };
        assert_eq!(
            empty
                .try_clone_with_nested_and_poll(&mut cancel_at_publication, &mut clone_no_loop,)
                .expect("empty-patch publication cancellation"),
            TrimmedPatchCloneRun::Cancelled
        );
        assert_eq!(replay_polls, total_polls);
    }

    #[test]
    fn trimmed_patch_copy_refuses_work_before_retained_bytes() {
        let error = preflight_trimmed_patch_copy_counts(usize::MAX, usize::MAX, usize::MAX)
            .expect_err("work must refuse before retained-byte accounting");
        assert!(matches!(
            error,
            NurbsError::Domain { ref what } if what.contains("work units above defensive ceiling")
        ));
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
