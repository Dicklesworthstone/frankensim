//! Tensor-product rational surfaces: two-stage de Boor evaluation
//! (generic scalar), EXACT directional knot insertion (Boehm along rows
//! or columns), first partial derivatives (f64), and per-span control
//! boxes (the convex-hull property in both directions).

use crate::NurbsError;
use crate::basis::{BASIS_MAX_WORK_UNITS, KnotVector, Scalar};
use crate::curve::NurbsCurve;

const SURFACE_VALIDATION_WORK_PER_CONTROL: u128 = 16;
const SURFACE_SPAN_BOX_WORK_PER_CONTROL: u128 = 16;
const SURFACE_SPAN_BOX_MAX_RETAINED_BYTES: u128 = 64 * 1024 * 1024;

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
#[derive(Debug, Clone, PartialEq)]
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
        Self::enforce_validation_work(Self::validation_work_for(&self.knots_u, &self.knots_v)?)?;
        self.knots_u.validate_live()?;
        self.knots_v.validate_live()?;
        let expected_u = self.knots_u.control_count();
        let expected_v = self.knots_v.control_count();
        if self.cpw.len() != expected_u
            || self.cpw.iter().any(|row| {
                row.len() != expected_v
                    || row.iter().any(|control| {
                        !control[3].is_admissible_weight()
                            || control
                                .iter()
                                .copied()
                                .any(|component| !component.is_finite())
                            || control[..3]
                                .iter()
                                .copied()
                                .any(|component| !component.quotient_is_finite(control[3]))
                    })
            })
        {
            return Err(NurbsError::Structure {
                what: "live surface control net must match both knot vectors and retain finite homogeneous coordinates with admissible weights"
                    .to_string(),
            });
        }
        Ok(())
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
        let mut cpw = Vec::new();
        cpw.try_reserve_exact(nu).map_err(|_| NurbsError::Domain {
            what: "surface control-row allocation was refused".to_string(),
        })?;
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
            let mut row = Vec::new();
            row.try_reserve_exact(nv).map_err(|_| NurbsError::Domain {
                what: "surface homogeneous-control row allocation was refused".to_string(),
            })?;
            for (p, &w) in prow.iter().zip(wrow) {
                let homogeneous = [p[0] * w, p[1] * w, p[2] * w, w];
                if p.iter()
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

    /// Validate this exact immutable surface snapshot once.
    ///
    /// # Errors
    /// [`NurbsError::Structure`] when the sealed source is internally invalid;
    /// [`NurbsError::Domain`] when validation work exceeds the defensive cap.
    pub fn admit(&self) -> Result<AdmittedNurbsSurface<'_, S>, NurbsError> {
        self.validate_live_structure()?;
        Ok(AdmittedNurbsSurface { inner: self })
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
        // u-partial: build the v-evaluated control column, differentiate
        // as a u-curve; symmetrically for v.
        let knots_v = self.knots_v();
        let (sv, bv) = knots_v.basis(v)?;
        let pv = knots_v.degree();
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
            u_net.push(acc);
        }
        let u_curve =
            NurbsCurve::<f64, 3>::from_homogeneous(self.inner.knots_u.try_clone()?, u_net)?;
        let du = u_curve.admitted_after_validation().derivatives(u, 1)?;
        let knots_u = self.knots_u();
        let (su, bu) = knots_u.basis(u)?;
        let pu = knots_u.degree();
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
            v_net.push(acc);
        }
        let v_curve =
            NurbsCurve::<f64, 3>::from_homogeneous(self.inner.knots_v.try_clone()?, v_net)?;
        let dv = v_curve.admitted_after_validation().derivatives(v, 1)?;
        Ok((du[0], du[1], dv[1]))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Rat;

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
}
