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
use crate::curve::NurbsCurve;
use crate::rat::Rat;

/// Defensive work ceiling for one exact trim classification across all loops.
/// This legacy cap bounds public allocation-bearing subdivision even when a
/// caller supplies `max_subdivision = u32::MAX`; explicit caller budgets belong
/// to the successor API.
pub(crate) const TRIM_CLASSIFY_MAX_WORK_UNITS: u128 = 1_048_576;

/// Minimum charge for admitting one caller-mutable loop before inspecting its
/// knot/control metadata. This makes a huge collection of individually tiny
/// loops reject in O(1), rather than spending unbounded time merely discovering
/// that the aggregate validation exceeds the legacy synchronous envelope.
const TRIM_MIN_LOOP_VALIDATION_WORK_UNITS: u128 = 64;

/// One closed trim loop: an exact rational curve in (u, v) parameter
/// space (closure is validated).
#[derive(Debug, Clone, PartialEq)]
pub struct TrimLoop {
    /// The exact 2-D curve.
    pub curve: NurbsCurve<Rat, 2>,
}

impl TrimLoop {
    fn validate_live(&self) -> Result<(), NurbsError> {
        self.curve.validate_live_structure()?;
        let (lo, hi) = self.curve.knots.domain()?;
        let start = self.curve.eval(lo)?;
        let end = self.curve.eval(hi)?;
        if start != end {
            return Err(NurbsError::Structure {
                what: "trim loop must close exactly (rational endpoint equality)".to_string(),
            });
        }

        // A full interior knot break carries independent left and right
        // limits. It is admissible in a general piecewise curve, but a trim
        // loop must be continuous for the control-polygon homotopy/winding
        // proof. Permit the expressive full-break representation only when
        // those limits agree exactly in Cartesian space.
        let p = self.curve.knots.degree;
        let mut run_start = 0usize;
        while run_start < self.curve.knots.knots.len() {
            let mut run_end = run_start + 1;
            while run_end < self.curve.knots.knots.len()
                && self.curve.knots.knots[run_end] == self.curve.knots.knots[run_start]
            {
                run_end += 1;
            }
            let is_interior = run_start != 0 && run_end != self.curve.knots.knots.len();
            if is_interior && run_end - run_start == p + 1 {
                let left = self.curve.cpw[run_start - 1];
                let right = self.curve.cpw[run_start];
                for coordinate in 0..2 {
                    if left[coordinate] * right[3] != right[coordinate] * left[3] {
                        return Err(NurbsError::Structure {
                            what: format!(
                                "trim loop is discontinuous at full knot break {:?}",
                                self.curve.knots.knots[run_start]
                            ),
                        });
                    }
                }
            }
            run_start = run_end;
        }
        Ok(())
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

    /// The same loop with reversed orientation (holes are wound opposite
    /// to outers under the nonzero rule): control points reversed, knot
    /// vector mirrored about the domain.
    ///
    /// # Errors
    /// [`NurbsError::Structure`] when caller mutation has invalidated closure,
    /// continuity, knots, or the control net.
    pub fn reversed_for_hole(&self) -> Result<TrimLoop, NurbsError> {
        self.validate_live()?;
        let (lo, hi) = self.curve.knots.domain()?;
        let knots: Vec<Rat> = self
            .curve
            .knots
            .knots
            .iter()
            .rev()
            .map(|&u| lo + hi - u)
            .collect();
        let cpw: Vec<[Rat; 4]> = self.curve.cpw.iter().rev().copied().collect();
        TrimLoop::new(NurbsCurve {
            knots: crate::basis::KnotVector {
                knots,
                degree: self.curve.knots.degree,
            },
            cpw,
        })
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
#[derive(Debug, Clone, PartialEq)]
pub struct TrimmedPatch {
    /// Outer boundary + hole loops (orientation encodes solidity via the
    /// nonzero-winding rule: outer CCW, holes CW).
    pub loops: Vec<TrimLoop>,
    /// Exact-subdivision depth before declaring `Boundary`.
    pub max_subdivision: u32,
}

impl TrimmedPatch {
    pub(crate) fn validate_live_with_budget(
        &self,
        work_remaining: &mut u128,
    ) -> Result<(), NurbsError> {
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
            trim_loop.validate_live()?;
        }
        Ok(())
    }

    /// Construct with the default certification depth.
    #[must_use]
    pub fn new(loops: Vec<TrimLoop>) -> Self {
        TrimmedPatch {
            loops,
            max_subdivision: 12,
        }
    }

    /// Certified classification of a parameter-space point.
    ///
    /// # Errors
    /// Propagates structural errors from exact subdivision.
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
    /// Returns [`NurbsError::Domain`] for an inverted box and propagates
    /// structural errors from exact subdivision.
    pub fn classify_box(&self, min: [Rat; 2], max: [Rat; 2]) -> Result<Classification, NurbsError> {
        let mut work_remaining = TRIM_CLASSIFY_MAX_WORK_UNITS;
        self.validate_live_with_budget(&mut work_remaining)?;
        if min[0] > max[0] || min[1] > max[1] {
            return Err(NurbsError::Domain {
                what: "trim classification box must be componentwise ordered".to_string(),
            });
        }
        let two = Rat::int(2);
        let witness = [(min[0] + max[0]) / two, (min[1] + max[1]) / two];
        let mut winding = 0i64;
        for l in &self.loops {
            match loop_winding_box(
                &l.curve,
                min,
                max,
                witness,
                self.max_subdivision,
                &mut work_remaining,
            )? {
                Some(w) => winding += w,
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
    // Closure evaluates both endpoints and the live public representation is
    // revalidated by the generic curve/basis APIs. Eight scans is conservative
    // for that present implementation and includes the full-break continuity
    // walk.
    scanned_entries
        .checked_mul(8)
        .map(|work| work.max(TRIM_MIN_LOOP_VALIDATION_WORK_UNITS))
        .ok_or_else(|| NurbsError::Domain {
            what: "trim repeated-validation accounting overflows u128".to_string(),
        })
}

/// Certified winding number of one closed rational curve about `q`, or
/// `None` when `q` cannot be separated from the curve within the
/// subdivision budget.
fn loop_winding_box(
    curve: &NurbsCurve<Rat, 2>,
    query_min: [Rat; 2],
    query_max: [Rat; 2],
    witness: [Rat; 2],
    max_depth: u32,
    work_remaining: &mut u128,
) -> Result<Option<i64>, NurbsError> {
    // Work in Bézier form so each span's control hull tightly bounds it.
    spend_trim_work(
        work_remaining,
        bezier_conversion_work(curve)?,
        "initial Bézier conversion",
    )?;
    let mut work = curve.to_bezier_form()?;
    let mut depth = 0u32;
    loop {
        spend_trim_work(
            work_remaining,
            bezier_conversion_work(&work)?,
            "span-box construction",
        )?;
        let boxes = work.span_boxes()?;
        let offending: Vec<(Rat, Rat)> = boxes
            .iter()
            .filter(|(min, max, _, _)| {
                max[0] >= query_min[0]
                    && min[0] <= query_max[0]
                    && max[1] >= query_min[1]
                    && min[1] <= query_max[1]
            })
            .map(|&(_, _, t0, t1)| (t0, t1))
            .collect();
        if offending.is_empty() {
            // Separated from the whole connected query box: winding is
            // constant throughout it, so one exact witness is sufficient.
            return Ok(Some(polygon_winding(&control_polygon(&work), witness)));
        }
        if depth >= max_depth {
            return Ok(None);
        }
        spend_trim_work(
            work_remaining,
            offending.len() as u128,
            "exact midpoint subdivision",
        )?;
        for (t0, t1) in offending {
            let mid = (t0 + t1) / Rat::int(2);
            // Exact midpoint insertion splits the offending span.
            work = work.insert_knot(mid)?;
        }
        spend_trim_work(
            work_remaining,
            bezier_conversion_work(&work)?,
            "post-subdivision Bézier conversion",
        )?;
        work = work.to_bezier_form()?;
        depth = depth.checked_add(1).ok_or_else(|| NurbsError::Domain {
            what: "trim subdivision depth overflows u32".to_string(),
        })?;
    }
}

fn bezier_conversion_work(curve: &NurbsCurve<Rat, 2>) -> Result<u128, NurbsError> {
    let size = (curve.cpw.len() as u128)
        .checked_add(curve.knots.knots.len() as u128)
        .ok_or_else(|| NurbsError::Domain {
            what: "trim structure size overflows u128".to_string(),
        })?;
    size.checked_mul(size)
        .and_then(|square| square.checked_mul(size))
        .ok_or_else(|| NurbsError::Domain {
            what: "trim Bézier work estimate overflows u128".to_string(),
        })
}

fn spend_trim_work(remaining: &mut u128, requested: u128, stage: &str) -> Result<(), NurbsError> {
    if requested > *remaining {
        return Err(NurbsError::Domain {
            what: format!(
                "trim {stage} requests {requested} work units with only {remaining} remaining from the {TRIM_CLASSIFY_MAX_WORK_UNITS}-unit defensive budget"
            ),
        });
    }
    *remaining -= requested;
    Ok(())
}

/// The Cartesian control polygon (weights divided out — the hull
/// property holds for rational Bézier segments in Cartesian space).
fn control_polygon(curve: &NurbsCurve<Rat, 2>) -> Vec<[Rat; 2]> {
    curve
        .cpw
        .iter()
        .map(|cp| [cp[0] / cp[3], cp[1] / cp[3]])
        .collect()
}

/// EXACT winding number of a closed polygon about `q` (crossing rule
/// with exact rational orientation tests — no epsilons anywhere).
fn polygon_winding(poly: &[[Rat; 2]], q: [Rat; 2]) -> i64 {
    let mut winding = 0i64;
    let n = poly.len();
    for i in 0..n {
        let a = poly[i];
        let b = poly[(i + 1) % n];
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
