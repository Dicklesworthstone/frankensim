//! Rational B-spline curves over a generic scalar: homogeneous de Boor
//! evaluation, derivatives to arbitrary order (f64 path), EXACT Boehm
//! knot insertion, Bézier decomposition, and EXACT degree elevation via
//! per-segment Bézier elevation (the elevated curve carries a
//! full-multiplicity knot vector — valid, evaluation-identical; minimal
//! knot vectors are a documented follow-up).

use crate::NurbsError;
use crate::basis::{KnotVector, Scalar};

/// One span's Cartesian control box: (min, max, t0, t1).
pub type SpanBox<S, const DIM: usize> = ([S; DIM], [S; DIM], S, S);

/// A rational curve in `DIM` dimensions: homogeneous control points
/// `(w·x…, w)` over a clamped knot vector.
#[derive(Debug, Clone, PartialEq)]
pub struct NurbsCurve<S: Scalar, const DIM: usize> {
    /// The knot vector.
    pub knots: KnotVector<S>,
    /// Homogeneous control points: `DIM` weighted coordinates + weight.
    pub cpw: Vec<[S; 4]>,
}

impl<S: Scalar, const DIM: usize> NurbsCurve<S, DIM> {
    pub(crate) fn validate_live_structure(&self) -> Result<(), NurbsError> {
        if DIM > 3 {
            return Err(NurbsError::Structure {
                what: format!("curve dimension {DIM} exceeds the homogeneous storage limit 3"),
            });
        }
        self.knots.validate_live()?;
        if self.cpw.len() != self.knots.control_count()
            || self.cpw.iter().any(|control| {
                !control[3].is_admissible_weight()
                    || control
                        .iter()
                        .copied()
                        .any(|component| !component.is_finite())
                    || control[..DIM]
                        .iter()
                        .copied()
                        .any(|component| !component.quotient_is_finite(control[3]))
            })
        {
            return Err(NurbsError::Structure {
                what: "live curve control net must match its knots and retain finite homogeneous coordinates with admissible weights"
                    .to_string(),
            });
        }
        Ok(())
    }

    /// Build from Cartesian control points + weights.
    ///
    /// # Errors
    /// [`NurbsError::Structure`] on count mismatch, non-finite coordinates, or
    /// non-positive/non-finite weights.
    pub fn new(
        knots: KnotVector<S>,
        points: &[[S; DIM]],
        weights: &[S],
    ) -> Result<Self, NurbsError> {
        knots.validate_live()?;
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
        if weights.iter().copied().any(|w| !w.is_admissible_weight()) {
            return Err(NurbsError::Structure {
                what: "weights must be finite, positive, and numerically admissible".to_string(),
            });
        }
        if points
            .iter()
            .flat_map(|point| point.iter())
            .copied()
            .any(|coordinate| !coordinate.is_finite())
        {
            return Err(NurbsError::Structure {
                what: "control-point coordinates must be finite".to_string(),
            });
        }
        let cpw: Vec<[S; 4]> = points
            .iter()
            .zip(weights)
            .map(|(p, &w)| {
                let mut h = [S::zero(); 4];
                for (slot, &c) in h.iter_mut().zip(p.iter()) {
                    *slot = c * w;
                }
                h[3] = w;
                h
            })
            .collect();
        if points.iter().zip(&cpw).any(|(point, homogeneous)| {
            point
                .iter()
                .zip(homogeneous.iter())
                .any(|(&coordinate, &weighted)| coordinate != S::zero() && weighted == S::zero())
        }) {
            return Err(NurbsError::Structure {
                what: "Cartesian coordinate × weight underflowed a nonzero coordinate to zero"
                    .to_string(),
            });
        }
        if cpw
            .iter()
            .flatten()
            .copied()
            .any(|component| !component.is_finite())
        {
            return Err(NurbsError::Structure {
                what: "Cartesian coordinate × weight overflowed the homogeneous numeric domain"
                    .to_string(),
            });
        }
        Ok(NurbsCurve { knots, cpw })
    }

    /// Homogeneous evaluation (the shared exact/fast core).
    ///
    /// # Errors
    /// [`NurbsError::Domain`] outside the domain.
    pub fn eval_homogeneous(&self, t: S) -> Result<[S; 4], NurbsError> {
        self.validate_live_structure()?;
        let (span, basis) = self.knots.basis(t)?;
        let p = self.knots.degree;
        let mut acc = [S::zero(); 4];
        for (r, &b) in basis.iter().enumerate() {
            let cp = &self.cpw[span - p + r];
            for (a, &c) in acc.iter_mut().zip(cp.iter()) {
                *a = *a + b * c;
            }
        }
        if acc.iter().copied().any(|component| !component.is_finite()) {
            return Err(NurbsError::Domain {
                what: "homogeneous curve evaluation left the finite numeric domain".to_string(),
            });
        }
        Ok(acc)
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

    /// EXACT Boehm knot insertion at `t` (multiplicity one per call).
    ///
    /// # Errors
    /// [`NurbsError::Domain`] outside the OPEN domain interior.
    pub fn insert_knot(&self, t: S) -> Result<Self, NurbsError> {
        self.validate_live_structure()?;
        let (lo, hi) = self.knots.domain()?;
        if t <= lo || hi <= t {
            return Err(NurbsError::Domain {
                what: format!("insertion parameter {t:?} must be interior to {lo:?}..{hi:?}"),
            });
        }
        let p = self.knots.degree;
        let k = self.knots.span(t)?;
        let mut new_cpw = Vec::with_capacity(self.cpw.len() + 1);
        new_cpw.extend_from_slice(&self.cpw[..=k - p]);
        for i in (k - p + 1)..=k {
            let denom = self.knots.knots[i + p] - self.knots.knots[i];
            let alpha = (t - self.knots.knots[i]) / denom;
            let mut q = [S::zero(); 4];
            for ((slot, &a), &b) in q.iter_mut().zip(&self.cpw[i - 1]).zip(&self.cpw[i]) {
                *slot = (S::one() - alpha) * a + alpha * b;
            }
            new_cpw.push(q);
        }
        new_cpw.extend_from_slice(&self.cpw[k..]);
        let mut new_knots = self.knots.knots.clone();
        new_knots.insert(k + 1, t);
        Ok(NurbsCurve {
            knots: KnotVector::new(new_knots, p)?,
            cpw: new_cpw,
        })
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
        self.validate_live_structure()?;
        let p = self.knots.degree;
        let (lo, hi) = self.knots.domain()?;
        if t <= lo || hi <= t || !self.knots.knots.contains(&t) {
            return Err(NurbsError::Domain {
                what: format!("{t:?} is not an interior knot"),
            });
        }
        // Index of the LAST occurrence of t.
        let r = self
            .knots
            .knots
            .iter()
            .rposition(|&u| u == t)
            .expect("contains checked");
        let mut new_knots = self.knots.knots.clone();
        new_knots.remove(r);
        let prior_multiplicity = new_knots.iter().filter(|&&knot| knot == t).count();
        // Insertion produced: Q_i = (1−α_i) P_{i−1} + α_i P_i over the
        // positive-alpha band i = k-p+1..=k-s, with α from the REMOVED
        // knot vector and prior multiplicity s. Rows after that band are exact
        // copies of the right suffix. Reconstruct the blended band forward;
        // the first suffix copy is an independent exact meet check.
        let k = r - 1; // span index of t in the removed vector
        let q = &self.cpw;
        let mut fwd: Vec<[S; 4]> = Vec::new(); // P_{k-p} .. computed forward
        let mut prev = q[k - p]; // P_{k-p} = Q_{k-p}
        fwd.push(prev);
        let blend_start = k - p + 1;
        let blend_end = k
            .checked_sub(prior_multiplicity)
            .ok_or_else(|| NurbsError::Structure {
                what: "knot-removal multiplicity exceeds its span index".to_string(),
            })?;
        for i in blend_start..=blend_end {
            let denom = new_knots[i + p] - new_knots[i];
            let alpha = (t - new_knots[i]) / denom;
            if alpha == S::zero() {
                return Err(NurbsError::Structure {
                    what: "degenerate removal alpha".to_string(),
                });
            }
            let mut pi = [S::zero(); 4];
            for ((slot, &qi), &pm) in pi.iter_mut().zip(&q[i]).zip(&prev) {
                *slot = (qi - (S::one() - alpha) * pm) / alpha;
            }
            fwd.push(pi);
            prev = pi;
        }
        let suffix_start = blend_end + 1;
        // Consistency: reconstructed P_{k-s} must equal the first untouched
        // suffix copy Q_{k-s+1} (= P_{k-s}).
        if fwd.last() != Some(&q[suffix_start]) {
            return Err(NurbsError::Structure {
                what: "knot is not exactly removable (curve genuinely uses it)".to_string(),
            });
        }
        let mut new_cpw: Vec<[S; 4]> = Vec::with_capacity(q.len() - 1);
        new_cpw.extend_from_slice(&q[..k - p]);
        new_cpw.extend_from_slice(&fwd[..fwd.len() - 1]);
        new_cpw.extend_from_slice(&q[suffix_start..]);
        let candidate = NurbsCurve {
            knots: KnotVector::new(new_knots, p)?,
            cpw: new_cpw,
        };
        // Exact end-to-end verifier: a successful removal must reproduce the
        // entire source representation under the public insertion algorithm,
        // not merely satisfy one local recurrence equation.
        if candidate.insert_knot(t)? != *self {
            return Err(NurbsError::Structure {
                what: "knot-removal candidate failed exact reinsertion verification".to_string(),
            });
        }
        Ok(candidate)
    }

    /// Decompose into Bézier segments by raising every interior knot to
    /// multiplicity `degree` (EXACT). Returns the refined curve.
    ///
    /// # Errors
    /// Propagates structural errors (none for valid inputs).
    pub fn to_bezier_form(&self) -> Result<Self, NurbsError> {
        self.validate_live_structure()?;
        let p = self.knots.degree;
        let mut cur = self.clone();
        loop {
            // Find an interior knot with multiplicity < p.
            let (lo, hi) = cur.knots.domain()?;
            let mut target = None;
            let mut i = 0;
            while i < cur.knots.knots.len() {
                let t = cur.knots.knots[i];
                let mut run_end = i + 1;
                while run_end < cur.knots.knots.len() && cur.knots.knots[run_end] == t {
                    run_end += 1;
                }
                if t > lo && t < hi && run_end - i < p {
                    target = Some(t);
                    break;
                }
                i = run_end;
            }
            match target {
                Some(t) => cur = cur.insert_knot(t)?,
                None => return Ok(cur),
            }
        }
    }

    /// EXACT degree elevation by one: decompose to Bézier form, elevate
    /// each segment with the exact binomial rule, and reassemble with a
    /// full-multiplicity knot vector. Evaluation is IDENTICAL (the
    /// conformance suite proves it with rational equality).
    ///
    /// # Errors
    /// Propagates structural/domain errors.
    pub fn elevate_degree(&self) -> Result<Self, NurbsError> {
        self.validate_live_structure()?;
        let p = self.knots.degree;
        let bez = self.to_bezier_form()?;
        // Collect distinct knots and their multiplicities in order. Ordinary
        // Bezier-form joins have multiplicity p and share one endpoint; a
        // legal full break has multiplicity p+1 and owns two independent
        // endpoints. Elevation must preserve that distinction.
        let mut breaks: Vec<S> = Vec::new();
        let mut multiplicities: Vec<usize> = Vec::new();
        for &u in &bez.knots.knots {
            if breaks.last() != Some(&u) {
                breaks.push(u);
                multiplicities.push(1);
            } else if let Some(multiplicity) = multiplicities.last_mut() {
                *multiplicity =
                    multiplicity
                        .checked_add(1)
                        .ok_or_else(|| NurbsError::Structure {
                            what: "degree-elevation knot multiplicity overflowed usize".to_string(),
                        })?;
            }
        }
        // Elevate each Bézier segment: Q_0 = P_0; Q_{p+1} = P_p;
        // Q_i = (i/(p+1)) P_{i-1} + (1 - i/(p+1)) P_i.
        let segment_spans: Vec<usize> = (p..bez.knots.control_count())
            .filter(|&span| bez.knots.knots[span] < bez.knots.knots[span + 1])
            .collect();
        let seg_count = breaks.len() - 1;
        if segment_spans.len() != seg_count {
            return Err(NurbsError::Structure {
                what: "degree elevation could not pair every distinct knot interval with one nonempty span"
                    .to_string(),
            });
        }
        let mut new_cpw: Vec<[S; 4]> = Vec::new();
        new_cpw
            .try_reserve(bez.cpw.len().saturating_add(seg_count))
            .map_err(|_| NurbsError::Domain {
                what: "degree-elevation control-net allocation was refused".to_string(),
            })?;
        let elevated_order = p.checked_add(2).ok_or_else(|| NurbsError::Structure {
            what: "degree elevation overflows spline-order arithmetic".to_string(),
        })?;
        let elevated_order_i64 =
            i64::try_from(elevated_order).map_err(|_| NurbsError::Structure {
                what: "degree elevation exceeds the scalar integer-lift domain".to_string(),
            })?;
        for (seg, &span) in segment_spans.iter().enumerate() {
            let pts = &bez.cpw[span - p..=span];
            let mut q = Vec::with_capacity(p + 2);
            q.push(pts[0]);
            for i in 1..=p {
                let numerator = i64::try_from(i).map_err(|_| NurbsError::Structure {
                    what: "degree elevation exceeds the scalar integer-lift domain".to_string(),
                })?;
                let alpha = S::from_int(numerator) / S::from_int(elevated_order_i64 - 1);
                let mut v = [S::zero(); 4];
                for ((slot, &a), &b) in v.iter_mut().zip(&pts[i - 1]).zip(&pts[i]) {
                    *slot = alpha * a + (S::one() - alpha) * b;
                }
                q.push(v);
            }
            q.push(pts[p]);
            if seg == 0 {
                new_cpw.extend_from_slice(&q);
            } else {
                let input_join_multiplicity = multiplicities[seg];
                match input_join_multiplicity {
                    m if m == p => {
                        // A Bezier-form C0 join shares its endpoint.
                        new_cpw.extend_from_slice(&q[1..]);
                    }
                    m if m == p + 1 => {
                        // A full break is discontinuous and owns both limiting
                        // endpoints. Do not manufacture continuity by dropping
                        // the right segment's first control point.
                        new_cpw.extend_from_slice(&q);
                    }
                    m => {
                        return Err(NurbsError::Structure {
                            what: format!(
                                "Bezier-form join multiplicity {m} is neither degree {p} nor full break {}",
                                p + 1
                            ),
                        });
                    }
                }
            }
        }
        // Elevation raises every multiplicity by one, preserving continuity
        // order. Endpoints therefore have p+2 copies, C0 joins p+1, and full
        // discontinuities p+2.
        let mut new_knots = Vec::new();
        new_knots
            .try_reserve(bez.knots.knots.len().saturating_add(breaks.len()))
            .map_err(|_| NurbsError::Domain {
                what: "degree-elevation knot allocation was refused".to_string(),
            })?;
        for (bi, (&b, &old_multiplicity)) in breaks.iter().zip(multiplicities.iter()).enumerate() {
            let mult = if bi == 0 || bi == breaks.len() - 1 {
                p + 2
            } else {
                old_multiplicity
                    .checked_add(1)
                    .ok_or_else(|| NurbsError::Structure {
                        what: "degree-elevation knot multiplicity overflowed usize".to_string(),
                    })?
            };
            for _ in 0..mult {
                new_knots.push(b);
            }
        }
        let elevated = NurbsCurve {
            knots: KnotVector::new(new_knots, p + 1)?,
            cpw: new_cpw,
        };
        elevated.validate_live_structure()?;
        Ok(elevated)
    }

    /// Per-span control-point bounding boxes in Cartesian space (the
    /// convex-hull property: each span's curve lies inside its box).
    /// Requires Bézier form for the tight per-segment claim; on general
    /// knot vectors the box of the span's `p+1` control points still
    /// bounds that span.
    ///
    /// # Errors
    /// Propagates domain errors (none for valid curves).
    pub fn span_boxes(&self) -> Result<Vec<SpanBox<S, DIM>>, NurbsError> {
        self.validate_live_structure()?;
        let p = self.knots.degree;
        let mut out = Vec::new();
        for span in p..self.knots.control_count() {
            let (t0, t1) = (self.knots.knots[span], self.knots.knots[span + 1]);
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
                }
                first = false;
            }
            out.push((min, max, t0, t1));
        }
        Ok(out)
    }
}

fn evaluate_homogeneous_derivative_net(
    net: &[[f64; 4]],
    knots: &[f64],
    degree: usize,
    t: f64,
) -> Result<[f64; 4], NurbsError> {
    let expected_knots = net
        .len()
        .checked_add(degree)
        .and_then(|count| count.checked_add(1))
        .ok_or_else(|| NurbsError::Structure {
            what: "derivative-net knot-count arithmetic overflowed".to_string(),
        })?;
    if net.is_empty()
        || knots.len() != expected_knots
        || knots.iter().any(|knot| !knot.is_finite())
        || knots.windows(2).any(|pair| pair[1] < pair[0])
    {
        return Err(NurbsError::Structure {
            what: "reduced homogeneous derivative net is malformed".to_string(),
        });
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
        let Some(span) = (0..=last_control)
            .rev()
            .find(|&candidate| knots[candidate] < knots[candidate + 1])
        else {
            return Err(NurbsError::Structure {
                what: "reduced derivative net has no nonempty upper span".to_string(),
            });
        };
        span
    } else {
        let mut span = degree;
        while span < last_control && knots[span + 1] <= t {
            span += 1;
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
        buffer
            .try_reserve_exact(basis_len)
            .map_err(|_| NurbsError::Domain {
                what: format!("derivative {stage} allocation was refused"),
            })?;
        buffer.resize(basis_len, 0.0);
    }
    basis[0] = 1.0;
    for j in 1..=degree {
        left[j] = t - knots[span + 1 - j];
        right[j] = knots[span + j] - t;
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
        }
        basis[j] = saved;
    }
    if basis.iter().any(|value| !value.is_finite()) {
        return Err(NurbsError::Domain {
            what: "reduced derivative basis left the finite numeric domain".to_string(),
        });
    }
    let mut value = [0.0; 4];
    for (offset, weight) in basis.into_iter().enumerate() {
        for (accumulator, control) in value.iter_mut().zip(net[span - degree + offset].iter()) {
            *accumulator += weight * control;
        }
    }
    if value.iter().any(|component| !component.is_finite()) {
        return Err(NurbsError::Domain {
            what: "reduced derivative evaluation left the finite numeric domain".to_string(),
        });
    }
    Ok(value)
}

impl<const DIM: usize> NurbsCurve<f64, DIM> {
    /// Defensive ceiling for the allocation- and quadratic-work-bearing legacy
    /// derivative API. Budgeted high-order jets belong to the typed successor.
    const MAX_DERIVATIVE_ORDER: usize = 64;
    /// Combined retained-net, quotient-recurrence, and allocation ceiling for
    /// the legacy whole-net derivative implementation.
    const MAX_DERIVATIVE_WORK_UNITS: u128 = 16_777_216;
    /// Hard retained payload bound for homogeneous nets and knot copies. Vec
    /// metadata and temporary basis arrays add a small bounded overhead.
    const MAX_DERIVATIVE_RETAINED_BYTES: u128 = 67_108_864;

    /// Derivatives up to `order` at `t` (rational quotient rule over the
    /// homogeneous derivative curves). Returns `[C(t), C'(t), …]`.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] outside the parameter domain or above the
    /// defensive legacy order ceiling.
    pub fn derivatives(&self, t: f64, order: usize) -> Result<Vec<[f64; DIM]>, NurbsError> {
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
        self.validate_live_structure()?;
        let p = self.knots.degree;
        let homogeneous_order = order.min(p);
        let structure_extent = (self.cpw.len() as u128)
            .checked_add(self.knots.knots.len() as u128)
            .ok_or_else(|| NurbsError::Domain {
                what: "derivative structure-size accounting overflows u128".to_string(),
            })?;
        let retained_nets = structure_extent
            .checked_mul((homogeneous_order as u128).saturating_add(1))
            .ok_or_else(|| NurbsError::Domain {
                what: "derivative retained-net accounting overflows u128".to_string(),
            })?;
        let retained_bytes = (self.cpw.len() as u128)
            .checked_mul((homogeneous_order as u128).saturating_add(1))
            .and_then(|count| count.checked_mul(core::mem::size_of::<[f64; 4]>() as u128))
            .and_then(|bytes| {
                (self.knots.knots.len() as u128)
                    .checked_mul((homogeneous_order as u128).saturating_add(1))
                    .and_then(|count| count.checked_mul(core::mem::size_of::<f64>() as u128))
                    .and_then(|knot_bytes| bytes.checked_add(knot_bytes))
            })
            .and_then(|bytes| {
                (p as u128)
                    .checked_add(1)
                    .and_then(|basis_len| basis_len.checked_mul(3))
                    .and_then(|basis_values| {
                        basis_values.checked_mul(core::mem::size_of::<f64>() as u128)
                    })
                    .and_then(|basis_bytes| bytes.checked_add(basis_bytes))
            })
            .ok_or_else(|| NurbsError::Domain {
                what: "derivative retained-byte accounting overflows u128".to_string(),
            })?;
        if retained_bytes > Self::MAX_DERIVATIVE_RETAINED_BYTES {
            return Err(NurbsError::Domain {
                what: format!(
                    "derivative request retains up to {retained_bytes} bytes, above ceiling {}",
                    Self::MAX_DERIVATIVE_RETAINED_BYTES
                ),
            });
        }
        let quotient_extent = (order as u128)
            .checked_add(1)
            .and_then(|side| side.checked_mul(side))
            .and_then(|square| square.checked_mul((DIM as u128).saturating_add(4)))
            .ok_or_else(|| NurbsError::Domain {
                what: "derivative quotient-work accounting overflows u128".to_string(),
            })?;
        let basis_extent = (0..=homogeneous_order).try_fold(0u128, |total, derivative| {
            let degree = p - derivative;
            let basis_order = degree.checked_add(1).ok_or_else(|| NurbsError::Domain {
                what: "derivative basis-order accounting overflows usize".to_string(),
            })?;
            let work = (degree as u128)
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
        if requested_work > Self::MAX_DERIVATIVE_WORK_UNITS {
            return Err(NurbsError::Domain {
                what: format!(
                    "derivative request needs {requested_work} defensive work units, above ceiling {}",
                    Self::MAX_DERIVATIVE_WORK_UNITS
                ),
            });
        }
        let (lo, hi) = self.knots.domain()?;
        if t > lo && t < hi {
            let multiplicity = self.knots.knots.iter().filter(|&&knot| knot == t).count();
            if multiplicity > 0 && (multiplicity > p || order > p - multiplicity) {
                return Err(NurbsError::Domain {
                    what: format!(
                        "ordinary derivative order {order} is undefined at interior knot multiplicity {multiplicity} for degree {p}; request an explicit one-sided jet in the successor API"
                    ),
                });
            }
        }
        // Homogeneous derivative control nets by repeated differencing.
        let mut nets: Vec<(Vec<[f64; 4]>, Vec<f64>, usize)> = Vec::new();
        nets.try_reserve_exact(homogeneous_order + 1)
            .map_err(|_| NurbsError::Domain {
                what: "derivative net-table allocation was refused".to_string(),
            })?;
        let mut initial_net = Vec::new();
        initial_net
            .try_reserve_exact(self.cpw.len())
            .map_err(|_| NurbsError::Domain {
                what: "derivative initial control-net allocation was refused".to_string(),
            })?;
        initial_net.extend_from_slice(&self.cpw);
        let mut initial_knots = Vec::new();
        initial_knots
            .try_reserve_exact(self.knots.knots.len())
            .map_err(|_| NurbsError::Domain {
                what: "derivative initial knot allocation was refused".to_string(),
            })?;
        initial_knots.extend_from_slice(&self.knots.knots);
        nets.push((initial_net, initial_knots, p));
        for k in 1..=homogeneous_order {
            let (prev, knots, deg) = &nets[k - 1];
            let mut next = Vec::new();
            next.try_reserve_exact(prev.len() - 1)
                .map_err(|_| NurbsError::Domain {
                    what: format!("derivative order {k} control-net allocation was refused"),
                })?;
            #[allow(clippy::cast_precision_loss)]
            let degf = *deg as f64;
            for i in 0..prev.len() - 1 {
                let denom = knots[i + deg + 1] - knots[i + 1];
                let mut d = [0.0f64; 4];
                if denom != 0.0 {
                    for (slot, (a, b)) in d.iter_mut().zip(prev[i + 1].iter().zip(&prev[i])) {
                        *slot = degf * (a - b) / denom;
                    }
                }
                next.push(d);
            }
            let mut new_knots = Vec::new();
            new_knots
                .try_reserve_exact(knots.len() - 2)
                .map_err(|_| NurbsError::Domain {
                    what: format!("derivative order {k} knot allocation was refused"),
                })?;
            new_knots.extend_from_slice(&knots[1..knots.len() - 1]);
            nets.push((next, new_knots, deg - 1));
        }
        // Evaluate each homogeneous derivative, then the quotient rule:
        // C^(k) = (A^(k) − Σ_{i=1..k} C(k−i) · w^(i) · binom(k,i)) / w.
        let mut hom = Vec::new();
        hom.try_reserve_exact(order + 1)
            .map_err(|_| NurbsError::Domain {
                what: "derivative homogeneous-jet allocation was refused".to_string(),
            })?;
        for (net, knots, degree) in &nets {
            hom.push(evaluate_homogeneous_derivative_net(net, knots, *degree, t)?);
        }
        // Polynomial homogeneous derivatives vanish above degree p, but a
        // rational quotient generally has nonzero derivatives of every order.
        // Retain those zero homogeneous jets so the quotient recurrence below
        // computes C^(k) correctly for k > p.
        hom.resize(order + 1, [0.0; 4]);
        let binom = |n: usize, k: usize| -> f64 {
            let mut b = 1.0f64;
            for j in 0..k {
                #[allow(clippy::cast_precision_loss)]
                {
                    b = b * (n - j) as f64 / (j + 1) as f64;
                }
            }
            b
        };
        let w0 = hom[0][3];
        let mut out: Vec<[f64; DIM]> = Vec::new();
        out.try_reserve_exact(order + 1)
            .map_err(|_| NurbsError::Domain {
                what: "derivative Cartesian-jet allocation was refused".to_string(),
            })?;
        for k in 0..=order {
            let mut num = [0.0f64; DIM];
            for (slot, &a) in num.iter_mut().zip(hom[k].iter()) {
                *slot = a;
            }
            for i in 1..=k {
                let c = binom(k, i) * hom[i][3];
                for (slot, prev) in num.iter_mut().zip(out[k - i].iter()) {
                    *slot -= c * prev;
                }
            }
            let jet = num.map(|v| v / w0);
            if jet.iter().any(|component| !component.is_finite()) {
                return Err(NurbsError::Domain {
                    what: format!(
                        "derivative order {k} left the finite floating-point numeric domain"
                    ),
                });
            }
            out.push(jet);
        }
        Ok(out)
    }
}
