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

mod layered {
//! Comparison for a stack of separately meshed solids, with finite contact.
//!
//! Let s run from the insulated end (0) to the cooled end (L). For unit
//! heating, q(s)=s. In layer j, w'=-s/k_j; at a contact s_i the trace drop
//! toward the cooled end is R_i*s_i. Set w(L)=L/h_min (zero for all D).
//! Thus max w = L/h_min + sum_j (b_j^2-a_j^2)/(2*k_j) + sum_i R_i*s_i.
//! The contact jump is essential: replacing it by continuity solves another
//! problem. B+ + F+*w and B- - F-*w are weak super/subsolutions of the
//! ORIGINAL volume, Robin and contact forms. The positive-part comparison
//! uses monotonicity of (a-b)(a_+-b_+), not a discrete maximum principle.
use super::*;

#[derive(Clone, Copy)]
struct Layer { low: f64, high: f64, k: f64 }

impl AffineSourceTetProblem<'_> {
    /// Whole-continuum temperature range for an axis-aligned layered stack.
    ///
    /// Each face-connected solid must have constant uncoupled axial
    /// conductivity; different solids may have different values. Tangential
    /// tensor blocks may vary. Solids occupy successive non-overlapping axis
    /// intervals, joined by exact matching contacts with one resistance per
    /// joining plane. All remaining exterior faces lie on the global bounding
    /// box. The cooled end is completely Dirichlet/positive Robin; other
    /// Neumann fluxes are zero. Other exterior Dirichlet/Robin laws remain.
    ///
    /// Every contact is admitted by the original geometry owner; neither
    /// vertex identities nor temperatures are welded. The supplied field is
    /// used only for original geometry/Dirichlet admission, never for extrema.
    /// Conforming, non-overlapping volume geometry remains a caller obligation.
    /// A refusal is no range, not evidence the physical solution is invalid.
    /// This is a nominal PDE comparison, not material/geometry uncertainty.
    pub fn layered_slab_temperature_range(
        &self, candidate: &[f64], axis: usize, reverse: bool,
        budget: FluxBudget, mut keep_going: impl FnMut() -> bool,
    ) -> Result<Iv, TetError> {
        poll(&mut keep_going)?;
        if axis >= 3 { return Err(TetError::Invalid("comparison axis must be 0, 1 or 2")); }
        let tensors = self.prepare(budget, &mut keep_going)?;
        let problem = self.problem(&tensors);
        let (cells, faces) = build(&problem, candidate, budget, &mut keep_going)?;
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for p in self.vertices {
            poll(&mut keep_going)?;
            for d in 0..3 { lo[d] = lo[d].min(p[d]); hi[d] = hi[d].max(p[d]); }
        }
        // Components are defined by actual shared faces, NEVER contact links.
        let mut component = vec![usize::MAX; cells.len()];
        let mut layers = Vec::new();
        let mut queue = VecDeque::new();
        for root in 0..cells.len() {
            poll(&mut keep_going)?;
            if component[root] != usize::MAX { continue; }
            let id = layers.len();
            let mut layer = Layer { low: f64::INFINITY, high: f64::NEG_INFINITY,
                k: self.conductivity[root][axis][axis] };
            component[root] = id;
            queue.push_back(root);
            while let Some(e) = queue.pop_front() {
                poll(&mut keep_going)?;
                let tensor = &self.conductivity[e];
                if tensor[axis][axis] != layer.k
                    || (0..3).any(|j| j != axis && tensor[axis][j] != 0.0)
                { return Err(TetError::Unsupported("each layer needs one constant uncoupled axial conductivity")); }
                for &v in &self.tets[e] {
                    layer.low = layer.low.min(self.vertices[v][axis]);
                    layer.high = layer.high.max(self.vertices[v][axis]);
                }
                for &f in &cells[e].faces {
                    if matches!(faces[f].condition, Some(BoundaryCondition::Contact { .. })) { continue; }
                    for &(other, _) in &faces[f].sides {
                        if component[other] == usize::MAX {
                            component[other] = id;
                            queue.push_back(other);
                        }
                    }
                }
            }
            layers.push(layer);
        }
        let mut order: Vec<_> = (0..layers.len()).collect();
        order.sort_unstable_by(|&a, &b| layers[a].low.total_cmp(&layers[b].low));
        poll(&mut keep_going)?;
        let mut rank = vec![0; layers.len()];
        for (i, &id) in order.iter().enumerate() {
            poll(&mut keep_going)?;
            rank[id] = i;
            if layers[id].low >= layers[id].high
                || (i > 0 && layers[order[i-1]].high != layers[id].low)
            { return Err(TetError::Unsupported("layers must be consecutive non-overlapping axis intervals")); }
        }
        if layers[order[0]].low != lo[axis] || layers[*order.last().unwrap()].high != hi[axis] {
            return Err(TetError::Unsupported("layer extents do not cover the bounding axis"));
        }
        let mut resistance = vec![None::<f64>; layers.len()-1];
        let end = if reverse { lo[axis] } else { hi[axis] };
        let mut b_lo = f64::INFINITY;
        let mut b_hi = f64::NEG_INFINITY;
        let mut h_end = f64::INFINITY;
        let mut end_faces = 0usize;
        for face in &faces {
            poll(&mut keep_going)?;
            if let Some(BoundaryCondition::Contact { resistance: r, .. }) = face.condition {
                if face.sides.len() != 2 { return Err(TetError::Invalid("contact was not linked")); }
                let a = rank[component[face.sides[0].0]];
                let b = rank[component[face.sides[1].0]];
                let left = a.min(b);
                if a.max(b) != left+1 { return Err(TetError::Unsupported("contact must join consecutive layers")); }
                let plane = layers[order[left]].high;
                if !face.vertices.iter().all(|&v| self.vertices[v][axis] == plane) {
                    return Err(TetError::Unsupported("contact must lie on the exact layer joining plane"));
                }
                if resistance[left].is_some_and(|old| old != r) {
                    return Err(TetError::Unsupported("layer comparison requires one resistance per contact plane"));
                }
                resistance[left] = Some(r);
                continue;
            }
            if face.sides.len() == 2 { continue; }
            let on = |d, value| face.vertices.iter().all(|&v| self.vertices[v][d] == value);
            if !(0..3).any(|d| on(d, lo[d]) || on(d, hi[d])) {
                return Err(TetError::Unsupported("noncontact exterior must lie on global bounding-box planes"));
            }
            let cooling = on(axis, end);
            if cooling { end_faces += 1; }
            let values = match face.condition {
                Some(BoundaryCondition::Dirichlet(values)) => Some(values),
                Some(BoundaryCondition::Robin { h, reference }) => {
                    if cooling { h_end = h_end.min(h); }
                    Some(reference)
                }
                Some(BoundaryCondition::Neumann(q)) if q == 0.0 && !cooling => None,
                _ => return Err(TetError::Unsupported("layered cooling end needs Dirichlet/Robin data; other fluxes must be zero")),
            };
            if let Some(values) = values {
                for value in values { b_lo = b_lo.min(value); b_hi = b_hi.max(value); }
            }
        }
        // Each component is internally connected; every successive pair is
        // linked, and the stack has a positive boundary anchor. No floating
        // solid may inherit another solid's cooling by coordinate coincidence.
        if resistance.iter().any(Option::is_none) || end_faces == 0
            || !b_lo.is_finite() || !b_hi.is_finite()
        { return Err(TetError::Unsupported("every layer needs a contact path to the cooling end")); }
        let mut positive = 0.0_f64;
        let mut negative = 0.0_f64;
        for values in self.source {
            poll(&mut keep_going)?;
            for &f in values { positive = positive.max(f); negative = negative.max(-f); }
        }
        if positive == 0.0 && negative == 0.0 {
            poll(&mut keep_going)?;
            return Ok(Iv { lo: b_lo, hi: b_hi });
        }
        let distance = |x: f64| if reverse { Iv::point(hi[axis]).sub(Iv::point(x)) }
            else { Iv::point(x).sub(Iv::point(lo[axis])) };
        let length = Iv::point(hi[axis]).sub(Iv::point(lo[axis]));
        if length.is_unbounded() || length.lo <= 0.0 { return Err(TetError::Unbounded); }
        let mut allowance = Iv::zero();
        for &id in &order {
            poll(&mut keep_going)?;
            let layer = layers[id];
            let width = Iv::point(layer.high).sub(Iv::point(layer.low));
            // (b^2-a^2) = (b-a)*(a+b), avoiding subtractive cancellation.
            let integral = width.mul(distance(layer.low).add(distance(layer.high)))
                .div_pos(Iv::point(layer.k).scale_pos(2.0));
            if integral.is_unbounded() { return Err(TetError::Unbounded); }
            allowance = allowance.add(integral);
        }
        for (i, r) in resistance.iter().enumerate() {
            poll(&mut keep_going)?;
            allowance = allowance.add(Iv::point(r.unwrap()).mul(distance(layers[order[i]].high)));
        }
        if h_end.is_finite() { allowance = allowance.add(length.div_pos(Iv::point(h_end))); }
        if allowance.is_unbounded() || allowance.hi <= 0.0 { return Err(TetError::Unbounded); }
        let range = Iv {
            lo: if negative == 0.0 { b_lo }
                else { Iv::point(b_lo).sub(Iv::point(negative).mul(allowance)).lo },
            hi: if positive == 0.0 { b_hi }
                else { Iv::point(b_hi).add(Iv::point(positive).mul(allowance)).hi },
        };
        poll(&mut keep_going)?;
        if range.is_unbounded() { return Err(TetError::Unbounded); }
        Ok(range)
    }
}

#[cfg(test)]
mod tests {
use super::*;

struct Stack {
    vertices: Vec<[f64;3]>, tets: Vec<[usize;4]>,
    conductivity: Vec<ConductivityTensor>, source: Vec<[f64;4]>, boundary: Vec<BoundaryFace>,
}
impl Stack {
    fn new(ks: &[f64], rs: &[f64], heat: f64) -> Self {
        assert_eq!(ks.len(), rs.len()+1);
        let mut vertices = Vec::new(); let mut tets = Vec::new();
        let mut conductivity = Vec::new();
        for (layer, &k) in ks.iter().enumerate() {
            let offset = vertices.len();
            vertices.extend((0..8).map(|i| [layer as f64+(i&1) as f64,
                ((i>>1)&1) as f64, ((i>>2)&1) as f64]));
            for axes in [[0,1,2],[0,2,1],[1,0,2],[1,2,0],[2,0,1],[2,1,0]] {
                let mut tet = [0;4];
                for (i,axis) in axes.into_iter().enumerate() { tet[i+1]=tet[i]|(1<<axis); }
                tets.push(tet.map(|v|v+offset));
                conductivity.push([[k,0.0,0.0],[0.0,3.0,0.25],[0.0,0.25,4.0]]);
            }
        }
        let mut incidence = BTreeMap::<[usize;3],usize>::new();
        for tet in &tets { for opposite in 0..4 {
            let mut f=[0;3]; let mut slot=0;
            for (i,&v) in tet.iter().enumerate() {
                if i!=opposite { f[slot]=v; slot+=1; }
            }
            f.sort_unstable(); *incidence.entry(f).or_default()+=1;
        }}
        let exterior: Vec<_> = incidence.into_iter().filter(|(_,n)|*n==1).map(|(f,_)|f).collect();
        let mut boundary=Vec::new();
        for &face in &exterior {
            let x=vertices[face[0]][0];
            let planar=face.iter().all(|&v|vertices[v][0]==x);
            let condition=if planar && x>0.0 && x<(ks.len() as f64) {
                let partner=*exterior.iter().find(|&&other| other!=face && face.iter().all(|&v|
                    other.iter().any(|&w|vertices[v]==vertices[w]))).unwrap();
                BoundaryCondition::Contact {partner,resistance:rs[x as usize-1]}
            } else if planar && x==ks.len() as f64 {
                BoundaryCondition::Robin {h:4.0,reference:[300.0;3]}
            } else { BoundaryCondition::Neumann(0.0) };
            boundary.push(BoundaryFace {vertices:face,condition});
        }
        let source=vec![[heat;4];tets.len()];
        Self {vertices,tets,conductivity,source,boundary}
    }
    fn problem(&self)->AffineSourceTetProblem<'_> {
        AffineSourceTetProblem {vertices:&self.vertices,tets:&self.tets,
            conductivity:&self.conductivity,source:&self.source,boundary:&self.boundary}
    }
    fn range(&self)->Result<Iv,TetError> {
        self.problem().layered_slab_temperature_range(&vec![300.0;self.vertices.len()],0,false,
            FluxBudget::default(),||true)
    }
}

#[test]
fn unequal_layer_conductivities_and_contact_jumps_bound_both_heat_directions() {
    for (ks,rs) in [(vec![2.0,8.0],vec![0.125]),(vec![2.0,8.0,4.0],vec![0.125,0.25])] {
        let mut rise=ks.len() as f64/4.0;
        for (i,&k) in ks.iter().enumerate() { rise+=((2*i+1) as f64)/(2.0*k); }
        for (i,&r) in rs.iter().enumerate() { rise+=r*(i+1) as f64; }
        for heat in [-4.0,0.0,4.0] {
            let stack=Stack::new(&ks,&rs,heat);
            let range=stack.range().unwrap(); let exact=300.0+heat*rise;
            assert!(range.lo<=exact && exact<=range.hi,"{range:?} excludes {exact}");
            let endpoint=if heat<0.0 {range.lo} else {range.hi};
            assert!((endpoint-exact).abs()<2e-10);
            if heat==0.0 { assert_eq!(range,Iv::point(300.0)); }
            assert!(stack.problem().slab_temperature_range(&vec![300.0;stack.vertices.len()],0,false,
                FluxBudget::default(),||true).is_err(),"legacy no-contact admission must not widen");
        }
    }
    let stack=Stack::new(&[2.0,8.0],&[0.125],4.0);
    assert!(stack.range().unwrap().hi>=304.25);
    let nearly_perfect=Stack::new(&[2.0,8.0],&[1e-12],4.0);
    assert!(stack.range().unwrap().hi-nearly_perfect.range().unwrap().hi>0.49,
        "dropping the contact jump must be detectable");
}

#[test]
fn contact_proof_is_independent_of_candidate_and_preserves_axis_reversal() {
    let mut stack=Stack::new(&[2.0,8.0],&[0.125],4.0);
    let expected=stack.range().unwrap();
    let changed=stack.problem().layered_slab_temperature_range(&vec![123.0;16],0,false,
        FluxBudget::default(),||true).unwrap();
    assert_eq!(expected,changed);
    for p in &mut stack.vertices { p.swap(0,2); p[2]=2.0-p[2]; }
    for k in &mut stack.conductivity { *k=[2,1,0].map(|i|[2,1,0].map(|j|k[i][j])); }
    let rotated=stack.problem().layered_slab_temperature_range(&vec![300.0;16],2,true,
        FluxBudget::default(),||true).unwrap();
    assert!((rotated.hi-expected.hi).abs()<1e-11);
    stack.boundary.reverse();
    let reordered=stack.problem().layered_slab_temperature_range(&vec![300.0;16],2,true,
        FluxBudget::default(),||true).unwrap();
    assert_eq!(rotated,reordered);
}

#[test]
fn incomplete_or_nonlayered_contacts_cannot_inherit_a_cooling_path() {
    let mut stack=Stack::new(&[2.0,8.0],&[0.125],4.0);
    let row=stack.boundary.iter().position(|f|matches!(f.condition,BoundaryCondition::Contact {..})).unwrap();
    let BoundaryCondition::Contact {partner,..}=stack.boundary[row].condition else {unreachable!()};
    let other=stack.boundary.iter().position(|f|f.vertices==partner).unwrap();
    stack.boundary[row].condition=BoundaryCondition::Neumann(0.0);
    stack.boundary[other].condition=BoundaryCondition::Neumann(0.0);
    assert!(stack.range().is_err(),"a half-open interior interface is not a complete slab");
    let mut stack=Stack::new(&[2.0,8.0],&[0.125],4.0);
    for i in [row,other] {
        if let BoundaryCondition::Contact {resistance,..}=&mut stack.boundary[i].condition { *resistance=0.25; }
    }
    assert!(matches!(stack.range(),Err(TetError::Unsupported(_))),"one face cannot represent every resistance");
    let mut stack=Stack::new(&[2.0,8.0],&[0.125],4.0);
    stack.conductivity[0][0][0]=2.5;
    assert!(matches!(stack.range(),Err(TetError::Unsupported(_))));
    let mut stack=Stack::new(&[2.0,8.0],&[0.125],4.0);
    for p in &mut stack.vertices {p[0]+=0.125*p[1];}
    assert!(matches!(stack.range(),Err(TetError::Unsupported(_))));
}

#[test]
fn signed_source_extrema_and_original_budgets_remain_mandatory() {
    let mut stack=Stack::new(&[2.0,8.0],&[0.125],4.0);
    for (e,values) in stack.source.iter_mut().enumerate() {
        *values=stack.tets[e].map(|v|4.0*stack.vertices[v][0]-4.0);
    }
    let range=stack.range().unwrap();
    assert!(range.lo<=295.75 && range.hi>=304.25,"zero net heat cannot erase local heating/cooling");
    let p=stack.problem(); let candidate=vec![300.0;stack.vertices.len()];
    assert_eq!(p.layered_slab_temperature_range(&candidate,0,false,
        FluxBudget {max_cells:11,max_iterations:128},||true).unwrap_err(),TetError::Budget);
    let mut polls=0;
    p.layered_slab_temperature_range(&candidate,0,false,FluxBudget::default(),||{polls+=1;true}).unwrap();
    for stop in [1,polls/2,polls-1] {
        let mut count=0;
        assert_eq!(p.layered_slab_temperature_range(&candidate,0,false,FluxBudget::default(),||{
            count+=1;count<stop
        }).unwrap_err(),TetError::Cancelled);
    }
    assert!(p.layered_slab_temperature_range(&candidate,3,false,FluxBudget::default(),||true).is_err());
}
}
}
