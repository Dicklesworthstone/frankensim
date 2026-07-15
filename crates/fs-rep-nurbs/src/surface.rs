//! Tensor-product rational surfaces: two-stage de Boor evaluation
//! (generic scalar), EXACT directional knot insertion (Boehm along rows
//! or columns), first partial derivatives (f64), and per-span control
//! boxes (the convex-hull property in both directions).

use crate::NurbsError;
use crate::basis::{BASIS_MAX_WORK_UNITS, BasisRun, KnotValidationOutcome, KnotVector, Scalar};
use crate::curve::NurbsCurve;
use fs_exec::Cx;

const SURFACE_VALIDATION_WORK_PER_CONTROL: u128 = 16;
const SURFACE_SPAN_BOX_WORK_PER_CONTROL: u128 = 16;
const SURFACE_SPAN_BOX_MAX_RETAINED_BYTES: u128 = 64 * 1024 * 1024;
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
    /// [`NurbsError::Domain`] when a destination allocation is refused.
    pub fn try_clone(&self) -> Result<Self, NurbsError> {
        let knots_u = self.knots_u.try_clone()?;
        let knots_v = self.knots_v.try_clone()?;
        let mut cpw = Vec::new();
        cpw.try_reserve_exact(self.cpw.len())
            .map_err(|_| NurbsError::Domain {
                what: "surface copy row-table allocation was refused".to_string(),
            })?;
        for source_row in &self.cpw {
            let mut row = Vec::new();
            row.try_reserve_exact(source_row.len())
                .map_err(|_| NurbsError::Domain {
                    what: "surface copy control-row allocation was refused".to_string(),
                })?;
            row.extend_from_slice(source_row);
            cpw.push(row);
        }
        Ok(NurbsSurface {
            knots_u,
            knots_v,
            cpw,
        })
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

    /// EXACT knot insertion in the u direction (Boehm on every column).
    ///
    /// # Errors
    /// [`NurbsError::Domain`] for a non-interior parameter.
    pub fn insert_knot_u(&self, t: S) -> Result<Self, NurbsError> {
        self.validate_live_structure()?;
        // Reuse the curve algorithm column-wise via a 1-D homogeneous
        // "curve" whose control points are rows of the net.
        let nv = self.knots_v.control_count();
        let mut new_rows: Option<Vec<Vec<[S; 4]>>> = None;
        let mut new_knots: Option<KnotVector<S>> = None;
        for j in 0..nv {
            let mut column = Vec::new();
            column
                .try_reserve_exact(self.cpw.len())
                .map_err(|_| NurbsError::Domain {
                    what: "surface u-insertion column allocation was refused".to_string(),
                })?;
            column.extend(self.cpw.iter().map(|row| row[j]));
            let curve = NurbsCurve::<S, 3>::from_homogeneous(self.knots_u.try_clone()?, column)?;
            let refined = curve.insert_knot(t)?;
            if new_rows.is_none() {
                let mut rows = Vec::new();
                rows.try_reserve_exact(refined.cpw.len())
                    .map_err(|_| NurbsError::Domain {
                        what: "surface u-insertion row-table allocation was refused".to_string(),
                    })?;
                for _ in 0..refined.cpw.len() {
                    let mut row = Vec::new();
                    row.try_reserve_exact(nv).map_err(|_| NurbsError::Domain {
                        what: "surface u-insertion output-row allocation was refused".to_string(),
                    })?;
                    rows.push(row);
                }
                new_rows = Some(rows);
            }
            let rows = new_rows.as_mut().expect("initialized above");
            for (i, cp) in refined.cpw.iter().enumerate() {
                rows[i].push(*cp);
            }
            new_knots = Some(refined.knots);
        }
        NurbsSurface::from_homogeneous(
            new_knots.expect("nv >= 1"),
            self.knots_v.try_clone()?,
            new_rows.expect("nv >= 1"),
        )
    }

    /// EXACT knot insertion in the v direction.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] for a non-interior parameter.
    pub fn insert_knot_v(&self, t: S) -> Result<Self, NurbsError> {
        self.validate_live_structure()?;
        let mut new_cpw = Vec::new();
        new_cpw
            .try_reserve_exact(self.cpw.len())
            .map_err(|_| NurbsError::Domain {
                what: "surface v-insertion row-table allocation was refused".to_string(),
            })?;
        let mut new_knots = None;
        for row in &self.cpw {
            let mut controls = Vec::new();
            controls
                .try_reserve_exact(row.len())
                .map_err(|_| NurbsError::Domain {
                    what: "surface v-insertion control-row allocation was refused".to_string(),
                })?;
            controls.extend_from_slice(row);
            let curve = NurbsCurve::<S, 3>::from_homogeneous(self.knots_v.try_clone()?, controls)?;
            let refined = curve.insert_knot(t)?;
            new_knots = Some(refined.knots);
            new_cpw.push(refined.cpw);
        }
        NurbsSurface::from_homogeneous(
            self.knots_u.try_clone()?,
            new_knots.expect("nu >= 1"),
            new_cpw,
        )
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
        let (span_u, basis_u) = match self.knots_u().basis_with_cx(u, cx)? {
            BasisRun::Complete { span, values } => (span, values),
            BasisRun::Cancelled => return Ok(SurfaceEvaluationRun::Cancelled),
        };
        let (span_v, basis_v) = match self.knots_v().basis_with_cx(v, cx)? {
            BasisRun::Complete { span, values } => (span, values),
            BasisRun::Cancelled => return Ok(SurfaceEvaluationRun::Cancelled),
        };
        self.eval_from_basis_with_poll(span_u, &basis_u, span_v, &basis_v, || {
            cx.checkpoint().is_err()
        })
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
        let mut out = Vec::new();
        out.try_reserve_exact(span_capacity)
            .map_err(|_| NurbsError::Domain {
                what: "surface span-box allocation was refused".to_string(),
            })?;
        for su in pu..knots_u.control_count() {
            let (u0, u1) = (knots_u.knots()[su], knots_u.knots()[su + 1]);
            if u1 <= u0 {
                continue;
            }
            for sv in pv..knots_v.control_count() {
                let (v0, v1) = (knots_v.knots()[sv], knots_v.knots()[sv + 1]);
                if v1 <= v0 {
                    continue;
                }
                let mut min = [S::zero(); 3];
                let mut max = [S::zero(); 3];
                let mut first = true;
                for row in &self.inner.cpw[su - pu..=su] {
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
                        }
                        first = false;
                    }
                }
                out.push((min, max, (u0, u1), (v0, v1)));
            }
        }
        Ok(out)
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
    /// Value and first partials without rescanning the source surface.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] outside the domain or when temporary isocurve
    /// construction/evaluation is refused.
    pub fn partials(&self, u: f64, v: f64) -> Result<SurfacePartials, NurbsError> {
        // Refuse the complete request before the first basis-workspace or
        // isocurve allocation. The U-first order is deterministic when both
        // parameters are invalid.
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
        NurbsCurve::<f64, 3>::preflight_derivative_request(knots_v, v, 1)?;

        let (su, bu) = knots_u.basis(u)?;
        let (sv, bv) = knots_v.basis(v)?;
        let pu = knots_u.degree();
        let pv = knots_v.degree();

        // u-partial: build the v-evaluated control column, differentiate
        // as a u-curve; symmetrically for v. Each scope releases its temporary
        // control net before the next one is allocated.
        let du = {
            let mut u_net = Vec::new();
            u_net
                .try_reserve_exact(self.inner.cpw.len())
                .map_err(|_| NurbsError::Domain {
                    what: "surface u-isocurve allocation was refused".to_string(),
                })?;
            for row in &self.inner.cpw {
                let mut acc = [0.0f64; 4];
                for (c, &wv) in bv.iter().enumerate() {
                    let cp = &row[sv - pv + c];
                    for (a, &x) in acc.iter_mut().zip(cp.iter()) {
                        *a += wv * x;
                    }
                }
                if acc.iter().any(|component| !component.is_finite())
                    || !acc[3].is_admissible_weight()
                {
                    return Err(NurbsError::Domain {
                        what: "surface u-isocurve left the admissible homogeneous domain"
                            .to_string(),
                    });
                }
                u_net.push(acc);
            }
            NurbsCurve::<f64, 3>::derivatives_from_admitted_parts_after_preflight(
                knots_u, &u_net, u, 1,
            )?
        };
        let dv = {
            let mut v_net = Vec::new();
            v_net
                .try_reserve_exact(knots_v.control_count())
                .map_err(|_| NurbsError::Domain {
                    what: "surface v-isocurve allocation was refused".to_string(),
                })?;
            for j in 0..knots_v.control_count() {
                let mut acc = [0.0f64; 4];
                for (r, &wu) in bu.iter().enumerate() {
                    let cp = &self.inner.cpw[su - pu + r][j];
                    for (a, &x) in acc.iter_mut().zip(cp.iter()) {
                        *a += wu * x;
                    }
                }
                if acc.iter().any(|component| !component.is_finite())
                    || !acc[3].is_admissible_weight()
                {
                    return Err(NurbsError::Domain {
                        what: "surface v-isocurve left the admissible homogeneous domain"
                            .to_string(),
                    });
                }
                v_net.push(acc);
            }
            NurbsCurve::<f64, 3>::derivatives_from_admitted_parts_after_preflight(
                knots_v, &v_net, v, 1,
            )?
        };
        Ok((du[0], du[1], dv[1]))
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
        assert!(matches!(
            admitted.partials(0.5, 0.5),
            Err(NurbsError::Domain { ref what }) if what.contains("multiplicity 2")
        ));
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
}
