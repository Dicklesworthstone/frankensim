//! Two cooled ends of the same admitted slab, without assuming symmetry.
//!
//! For unit source, set r_i=1/min(h_i) (zero for an all-Dirichlet end),
//! d=r_0+r_1+L/k and q_0=L*(r_1+L/(2k))/d, q_1=L-q_0.
//! w(s)=r_0*q_0+(q_0*s-s^2/2)/k solves -k*w''=1 with these
//! homogeneous Robin/Dirichlet end conditions. Its nonnegative maximum is
//! r_0*q_0+q_0^2/(2k)=r_1*q_1+q_1^2/(2k), attained at s=q_0.
//! Consequently B+ + F+*w and B- - F-*w are comparison barriers, including
//! mixed end patches with larger h and any side-face Dirichlet/Robin data.
//! Their conormal flux is continuous because K*e_s=k*e_s everywhere.
//! This is the weak comparison argument of the one-ended owner, with both
//! end conditions retained; it does not infer symmetry of the physical field.
use super::*;

// All arithmetic defining the maximum allowance is outward rounded. Two
// algebraically equal expressions avoid an unnecessarily wide result when
// one end is very poorly cooled. Taking the lesser upper bound is sound;
// neither interval's midpoint is used as a proof coefficient.
fn allowance(length: Iv, k: f64, h: [f64; 2]) -> Result<Iv, TetError> {
    let resistance = |h: f64| if h == f64::INFINITY { Iv::point(0.0) }
        else { Iv::point(1.0).div_pos(Iv::point(h)) };
    let r = h.map(resistance);
    let lk = length.div_pos(Iv::point(k));
    let denominator = r[0].add(r[1]).add(lk);
    let q0 = length.mul(r[1].add(lk.scale_pos(0.5))).div_pos(denominator);
    let q1 = length.mul(r[0].add(lk.scale_pos(0.5))).div_pos(denominator);
    let two_k = Iv::point(k).scale_pos(2.0);
    let left = r[0].mul(q0).add(q0.sq().div_pos(two_k));
    let right = r[1].mul(q1).add(q1.sq().div_pos(two_k));
    if [r[0], r[1], lk, denominator, q0, q1, left, right].iter().any(Iv::is_unbounded)
        || denominator.lo <= 0.0
    { return Err(TetError::Unbounded); }
    let upper = left.hi.min(right.hi);
    if upper <= 0.0 { return Err(TetError::Unbounded); }
    Ok(Iv { lo: 0.0, hi: upper })
}

impl AffineSourceTetProblem<'_> {
    /// Enclose all continuum temperatures using BOTH cooled ends of a slab.
    ///
    /// Both axis ends must consist entirely of prescribed-temperature or
    /// positive-Robin faces. End coefficients need not agree; each end uses
    /// its own minimum h among Robin patches. Dirichlet patches never weaken
    /// that minimum. All other boundary/source/tensor/geometry restrictions
    /// are those of `slab_temperature_range`. Actual reference temperatures
    /// on every face enter the extrema, not only the end temperatures.
    ///
    /// This can establish a finite material domain that a one-ended barrier
    /// cannot, without adding a primal solve or using nodal extrema. It is a
    /// nominal PDE comparison, not a mesh-error or physical-uncertainty bound.
    /// Non-overlapping conforming volume geometry remains a caller obligation.
    ///
    /// # Errors
    /// Incomplete end cooling, contact, nonzero side flux, unsupported tensors
    /// or geometry, invalid fields, unanchored components, exhausted budgets,
    /// nonfinite interval arithmetic or cancellation yield no partial range.
    pub fn slab_temperature_range_two_sided(
        &self, candidate: &[f64], axis: usize, budget: FluxBudget,
        mut keep_going: impl FnMut() -> bool,
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
        for point in self.vertices {
            poll(&mut keep_going)?;
            for d in 0..3 { lo[d] = lo[d].min(point[d]); hi[d] = hi[d].max(point[d]); }
        }
        let length = Iv::point(hi[axis]).sub(Iv::point(lo[axis]));
        if length.is_unbounded() || length.lo <= 0.0 { return Err(TetError::Unbounded); }
        let mut ends = [0usize; 2];
        let mut h = [f64::INFINITY; 2];
        let mut b_lo = f64::INFINITY;
        let mut b_hi = f64::NEG_INFINITY;
        let mut anchored = vec![false; cells.len()];
        let mut queue = VecDeque::new();
        for face in &faces {
            poll(&mut keep_going)?;
            if matches!(face.condition, Some(BoundaryCondition::Contact { .. })) {
                return Err(TetError::Unsupported("slab comparison does not replace finite contact by continuity"));
            }
            if face.sides.len() == 2 { continue; }
            let on = |d, value| face.vertices.iter().all(|&v| self.vertices[v][d] == value);
            if !(0..3).any(|d| on(d, lo[d]) || on(d, hi[d])) {
                return Err(TetError::Unsupported("slab exterior must lie on bounding-box planes"));
            }
            let end = if on(axis, lo[axis]) { Some(0) }
                else if on(axis, hi[axis]) { Some(1) } else { None };
            if let Some(i) = end { ends[i] += 1; }
            let values = match face.condition {
                Some(BoundaryCondition::Dirichlet(values)) => Some(values),
                Some(BoundaryCondition::Robin { h: coefficient, reference }) => {
                    if let Some(i) = end { h[i] = h[i].min(coefficient); }
                    Some(reference)
                }
                Some(BoundaryCondition::Neumann(q)) if q == 0.0 && end.is_none() => None,
                _ => return Err(TetError::Unsupported("both slab ends need complete Dirichlet/Robin coverage; other fluxes must be zero")),
            };
            if let Some(values) = values {
                for value in values { b_lo = b_lo.min(value); b_hi = b_hi.max(value); }
                let e = face.sides[0].0;
                if !anchored[e] { anchored[e] = true; queue.push_back(e); }
            }
        }
        if ends.contains(&0) || !b_lo.is_finite() || !b_hi.is_finite() {
            return Err(TetError::Unsupported("both slab ends need explicit cooling boundaries"));
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
        // No arithmetic slack or unnecessary resistance overflow for zero heat.
        let range = if positive == 0.0 && negative == 0.0 { Iv { lo: b_lo, hi: b_hi } }
        else {
            let a = allowance(length, k, h)?;
            Iv {
                lo: if negative == 0.0 { b_lo }
                    else { Iv::point(b_lo).sub(Iv::point(negative).mul(a)).lo },
                hi: if positive == 0.0 { b_hi }
                    else { Iv::point(b_hi).add(Iv::point(positive).mul(a)).hi },
            }
        };
        poll(&mut keep_going)?;
        if range.is_unbounded() { return Err(TetError::Unbounded); }
        Ok(range)
    }
}
