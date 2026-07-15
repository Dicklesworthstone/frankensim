//! Tensor-product rational surfaces: two-stage de Boor evaluation
//! (generic scalar), EXACT directional knot insertion (Boehm along rows
//! or columns), first partial derivatives (f64), and per-span control
//! boxes (the convex-hull property in both directions).

use crate::NurbsError;
use crate::basis::{KnotVector, Scalar};
use crate::curve::NurbsCurve;

/// One (u-span × v-span) control box: (min, max, (u0, u1), (v0, v1)).
pub type SurfaceSpanBox<S> = ([S; 3], [S; 3], (S, S), (S, S));

/// Value + first partials `(S, S_u, S_v)`.
pub type SurfacePartials = ([f64; 3], [f64; 3], [f64; 3]);

/// A rational tensor-product surface in 3D.
#[derive(Debug, Clone, PartialEq)]
pub struct NurbsSurface<S: Scalar> {
    /// Knots in u.
    pub knots_u: KnotVector<S>,
    /// Knots in v.
    pub knots_v: KnotVector<S>,
    /// Homogeneous control net `cpw[i][j]`, i along u, j along v.
    pub cpw: Vec<Vec<[S; 4]>>,
}

impl<S: Scalar> NurbsSurface<S> {
    fn validate_live_structure(&self) -> Result<(), NurbsError> {
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
    /// coordinates, or non-positive/non-finite weights.
    pub fn new(
        knots_u: KnotVector<S>,
        knots_v: KnotVector<S>,
        points: &[Vec<[S; 3]>],
        weights: &[Vec<S>],
    ) -> Result<Self, NurbsError> {
        knots_u.validate_live()?;
        knots_v.validate_live()?;
        let (nu, nv) = (knots_u.control_count(), knots_v.control_count());
        if points.len() != nu || weights.len() != nu {
            return Err(NurbsError::Structure {
                what: format!("expected {nu} control rows, got {}", points.len()),
            });
        }
        let mut cpw = Vec::with_capacity(nu);
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
            let mut row = Vec::with_capacity(nv);
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

    /// Homogeneous evaluation.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] outside the domain.
    pub fn eval_homogeneous(&self, u: S, v: S) -> Result<[S; 4], NurbsError> {
        self.validate_live_structure()?;
        let (su, bu) = self.knots_u.basis(u)?;
        let (sv, bv) = self.knots_v.basis(v)?;
        let (pu, pv) = (self.knots_u.degree, self.knots_v.degree);
        let mut acc = [S::zero(); 4];
        for (r, &wu) in bu.iter().enumerate() {
            for (c, &wv) in bv.iter().enumerate() {
                let cp = &self.cpw[su - pu + r][sv - pv + c];
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

    /// Cartesian evaluation.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] outside the domain.
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
            let column: Vec<[S; 4]> = self.cpw.iter().map(|row| row[j]).collect();
            let curve = NurbsCurve::<S, 3> {
                knots: self.knots_u.clone(),
                cpw: column,
            };
            let refined = curve.insert_knot(t)?;
            let rows =
                new_rows.get_or_insert_with(|| vec![Vec::with_capacity(nv); refined.cpw.len()]);
            for (i, cp) in refined.cpw.iter().enumerate() {
                rows[i].push(*cp);
            }
            new_knots = Some(refined.knots);
        }
        Ok(NurbsSurface {
            knots_u: new_knots.expect("nv >= 1"),
            knots_v: self.knots_v.clone(),
            cpw: new_rows.expect("nv >= 1"),
        })
    }

    /// EXACT knot insertion in the v direction.
    ///
    /// # Errors
    /// [`NurbsError::Domain`] for a non-interior parameter.
    pub fn insert_knot_v(&self, t: S) -> Result<Self, NurbsError> {
        self.validate_live_structure()?;
        let mut new_cpw = Vec::with_capacity(self.cpw.len());
        let mut new_knots = None;
        for row in &self.cpw {
            let curve = NurbsCurve::<S, 3> {
                knots: self.knots_v.clone(),
                cpw: row.clone(),
            };
            let refined = curve.insert_knot(t)?;
            new_knots = Some(refined.knots.clone());
            new_cpw.push(refined.cpw);
        }
        Ok(NurbsSurface {
            knots_u: self.knots_u.clone(),
            knots_v: new_knots.expect("nu >= 1"),
            cpw: new_cpw,
        })
    }

    /// Per-(u-span × v-span) Cartesian control boxes: each patch of the
    /// surface lies inside its sub-net's box (convex-hull property).
    ///
    /// # Errors
    /// Returns [`NurbsError::Structure`] when the mutable public representation
    /// no longer satisfies its knot/control-net invariants.
    pub fn span_boxes(&self) -> Result<Vec<SurfaceSpanBox<S>>, NurbsError> {
        self.validate_live_structure()?;
        let (pu, pv) = (self.knots_u.degree, self.knots_v.degree);
        let mut out = Vec::new();
        for su in pu..self.knots_u.control_count() {
            let (u0, u1) = (self.knots_u.knots[su], self.knots_u.knots[su + 1]);
            if u1 <= u0 {
                continue;
            }
            for sv in pv..self.knots_v.control_count() {
                let (v0, v1) = (self.knots_v.knots[sv], self.knots_v.knots[sv + 1]);
                if v1 <= v0 {
                    continue;
                }
                let mut min = [S::zero(); 3];
                let mut max = [S::zero(); 3];
                let mut first = true;
                for row in &self.cpw[su - pu..=su] {
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
        self.validate_live_structure()?;
        // u-partial: build the v-evaluated control column, differentiate
        // as a u-curve; symmetrically for v.
        let (sv, bv) = self.knots_v.basis(v)?;
        let pv = self.knots_v.degree;
        let u_net: Vec<[f64; 4]> = self
            .cpw
            .iter()
            .map(|row| {
                let mut acc = [0.0f64; 4];
                for (c, &wv) in bv.iter().enumerate() {
                    let cp = &row[sv - pv + c];
                    for (a, &x) in acc.iter_mut().zip(cp.iter()) {
                        *a += wv * x;
                    }
                }
                acc
            })
            .collect();
        let u_curve = NurbsCurve::<f64, 3> {
            knots: self.knots_u.clone(),
            cpw: u_net,
        };
        let du = u_curve.derivatives(u, 1)?;
        let (su, bu) = self.knots_u.basis(u)?;
        let pu = self.knots_u.degree;
        let v_net: Vec<[f64; 4]> = (0..self.knots_v.control_count())
            .map(|j| {
                let mut acc = [0.0f64; 4];
                for (r, &wu) in bu.iter().enumerate() {
                    let cp = &self.cpw[su - pu + r][j];
                    for (a, &x) in acc.iter_mut().zip(cp.iter()) {
                        *a += wu * x;
                    }
                }
                acc
            })
            .collect();
        let v_curve = NurbsCurve::<f64, 3> {
            knots: self.knots_v.clone(),
            cpw: v_net,
        };
        let dv = v_curve.derivatives(v, 1)?;
        Ok((du[0], du[1], dv[1]))
    }
}
