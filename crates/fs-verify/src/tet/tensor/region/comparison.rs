//! Pointwise comparison barriers, independent of the supplied numerical field.
//!
//! This deliberately narrow producer recognizes an axis-aligned slab whose
//! exterior triangles lie on its bounding-box planes. In direction s from the
//! insulated end, K*e_s = k*e_s is identical in every cell. Tangential tensor
//! blocks may differ. Nonzero Neumann data and contact are not admitted.
//!
//! Let F+ = max(0, sup f), F- = max(0, -inf f), L be the slab length, and
//! B+/B- the extrema of ALL Dirichlet and Robin reference data. With h* the
//! least Robin coefficient on the cooling end (1/h*=0 if that end is all D),
//! U(s)=B+ + F+*[L/h*+(L^2-s^2)/(2*k)] is a supersolution, and
//! V(s)=B- - F-*[L/h*+(L^2-s^2)/(2*k)] is a subsolution.
//! Their conormal fluxes are continuous across interior coefficient jumps.
//! On the other box planes the normal derivative vanishes, except where
//! Dirichlet data may be imposed. At the cooling end h*(U-reference)>=F+*L.
//! Testing the weak equation with (T-U)+ or (V-T)+ proves V<=T<=U a.e.
//! The anchor check excludes floating pure-Neumann components. There is no
//! assumption of a discrete maximum principle or a converged candidate.
use super::*;

impl AffineSourceTetProblem<'_> {
    /// Enclose EVERY continuum temperature of an admitted linear slab.
    ///
    /// `axis` is 0, 1 or 2. With `reverse=false`, the low-coordinate end is
    /// the insulated end of the comparison barrier and the high end must
    /// contain only Dirichlet/positive Robin faces; `reverse=true` swaps them.
    /// Actual Dirichlet/Robin data on other sides are retained, not removed.
    /// Source extrema include all four original P1 values of every element.
    ///
    /// Geometry and original boundary admission reuse the flux owner. The
    /// candidate is used ONLY for those checks, never for the continuum range.
    /// Domain conformity/non-overlap remains the original caller obligation.
    /// This is a nominal PDE bound, not CAD, material or model uncertainty.
    ///
    /// # Errors
    /// Unsupported geometry, directional conductivity, contact or nonzero
    /// Neumann data, invalid inputs, unanchored components, budget exhaustion,
    /// nonfinite interval arithmetic, or cancellation return no range.
    pub fn slab_temperature_range(
        &self, candidate: &[f64], axis: usize, reverse: bool,
        budget: FluxBudget, mut keep_going: impl FnMut() -> bool,
    ) -> Result<Iv, TetError> {
        poll(&mut keep_going)?;
        if axis >= 3 { return Err(TetError::Invalid("comparison axis must be 0, 1 or 2")); }
        let tensors = self.prepare(budget, &mut keep_going)?;
        let problem = self.problem(&tensors);
        let (cells, faces) = build(&problem, candidate, budget, &mut keep_going)?;
        let k = self.conductivity[0][axis][axis];
        for tensor in self.conductivity {
            poll(&mut keep_going)?;
            if tensor[axis][axis] != k || (0..3).any(|j| j != axis && tensor[axis][j] != 0.0) {
                return Err(TetError::Unsupported("slab comparison requires one constant uncoupled axial conductivity"));
            }
        }
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in self.vertices {
            poll(&mut keep_going)?;
            for d in 0..3 { lo[d] = lo[d].min(p[d]); hi[d] = hi[d].max(p[d]); }
        }
        let end = if reverse { lo[axis] } else { hi[axis] };
        let mut b_lo = f64::INFINITY;
        let mut b_hi = f64::NEG_INFINITY;
        let mut h_end = f64::INFINITY;
        let mut end_faces = 0usize;
        let mut anchored = vec![false; cells.len()];
        let mut queue = VecDeque::new();
        for face in &faces {
            poll(&mut keep_going)?;
            if matches!(face.condition, Some(BoundaryCondition::Contact { .. })) {
                return Err(TetError::Unsupported("slab comparison does not replace finite contact by continuity"));
            }
            if face.sides.len() == 2 { continue; }
            let on_plane = |d, value| face.vertices.iter().all(|&v| self.vertices[v][d] == value);
            if !(0..3).any(|d| on_plane(d, lo[d]) || on_plane(d, hi[d])) {
                return Err(TetError::Unsupported("slab exterior must lie on bounding-box planes"));
            }
            let cooling = on_plane(axis, end);
            if cooling { end_faces += 1; }
            let values = match face.condition {
                Some(BoundaryCondition::Dirichlet(values)) => Some(values),
                Some(BoundaryCondition::Robin { h, reference }) => {
                    if cooling { h_end = h_end.min(h); }
                    Some(reference)
                }
                Some(BoundaryCondition::Neumann(q)) if q == 0.0 && !cooling => None,
                _ => return Err(TetError::Unsupported("slab cooling end needs Dirichlet/Robin data; other fluxes must be zero")),
            };
            if let Some(values) = values {
                for value in values { b_lo = b_lo.min(value); b_hi = b_hi.max(value); }
                let e = face.sides[0].0;
                if !anchored[e] { anchored[e] = true; queue.push_back(e); }
            }
        }
        if end_faces == 0 || !b_lo.is_finite() || !b_hi.is_finite() {
            return Err(TetError::Unsupported("slab comparison needs an explicit cooling boundary"));
        }
        while let Some(e) = queue.pop_front() {
            poll(&mut keep_going)?;
            for &f in &cells[e].faces {
                for &(other, _) in &faces[f].sides {
                    if !anchored[other] { anchored[other] = true; queue.push_back(other); }
                }
            }
        }
        if anchored.iter().any(|a| !a) {
            return Err(TetError::Unsupported("unanchored component has no pointwise comparison bound"));
        }
        let mut positive = 0.0_f64;
        let mut negative = 0.0_f64;
        for row in self.source {
            poll(&mut keep_going)?;
            for &f in row { positive = positive.max(f); negative = negative.max(-f); }
        }
        let length = Iv::point(hi[axis]).sub(Iv::point(lo[axis]));
        if length.is_unbounded() || length.lo <= 0.0 { return Err(TetError::Unbounded); }
        let mut allowance = length.sq().div_pos(Iv::point(k).scale_pos(2.0));
        if h_end.is_finite() { allowance = allowance.add(length.div_pos(Iv::point(h_end))); }
        if allowance.is_unbounded() { return Err(TetError::Unbounded); }
        // Exact zero source needs no multiplication, rounding slack or invented
        // uncertainty: the constant comparison barriers are the boundary extrema.
        let lower = if negative == 0.0 { b_lo }
            else { Iv::point(b_lo).sub(Iv::point(negative).mul(allowance)).lo };
        let upper = if positive == 0.0 { b_hi }
            else { Iv::point(b_hi).add(Iv::point(positive).mul(allowance)).hi };
        let range = Iv { lo: lower, hi: upper };
        poll(&mut keep_going)?;
        if range.is_unbounded() { return Err(TetError::Unbounded); }
        Ok(range)
    }
}

#[cfg(test)]
mod tests;

mod two_sided;
#[cfg(test)]
mod two_sided_tests;
