//! Sparse Galerkin assembly from the actual constrained CutFEM cell blocks.
//! Geometry is retained; only density-dependent numerical data are rebuilt.
use super::*;
use fs_solver::op::multilevel::{MultilevelBudget, MultilevelControl, MultilevelError, SparseMultilevel, SymmetricGalerkin};

/// Geometry and numerical setup limits for recursive correction spaces.
#[derive(Debug, Clone, Copy)]
pub struct AdaptiveMultilevelOptions3 {
    /// Total scalar coefficients in all geometric transfer constructions and,
    /// separately, in the composed physical-node interpolation during setup.
    pub max_transfer_terms: usize,
    /// Exact constrained-diagonal accumulation allowance per density.
    pub max_diagonal_contributions: usize,
    /// Whole sparse hierarchy limits, including the small bottom factor.
    pub hierarchy: MultilevelBudget,
}
impl Default for AdaptiveMultilevelOptions3 {
    fn default() -> Self {
        Self { max_transfer_terms: 4_000_000, max_diagonal_contributions: 100_000_000,
            hierarchy: MultilevelBudget::default() }
    }
}
/// Two distinct setup operations, not interchangeable cost units. Partial
/// Galerkin construction and cancelled fine-application setup remain observable.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct AdaptiveSetupWork3 {
    pub operator_applications: usize,
    pub galerkin_products: usize,
}

/// Private geometry payload, always owned beside the fine operator.
pub(super) struct AdaptiveHierarchy3 {
    pub(super) transfers: Vec<Csr>,
    pub(super) options: AdaptiveMultilevelOptions3,
}
impl AdaptiveHierarchy3 {
    pub(super) fn new(fine: &AdaptiveElasticity3, coarser: &[&AdaptiveElasticity3],
        options: AdaptiveMultilevelOptions3, mut checkpoint: impl FnMut() -> ControlFlow<()>)
        -> Result<Self, AdaptivePreconditionError3> {
        let fail = AdaptivePreconditionError3::Hierarchy;
        poll(&mut checkpoint).map_err(AdaptivePreconditionError3::Physics)?;
        let budget = options.hierarchy;
        if !(2..=20).contains(&budget.max_levels) || !(1..=512).contains(&budget.max_coarsest_dofs)
            || coarser.is_empty() || coarser.len() >= budget.max_levels || fine.n() > budget.max_fine_dofs {
            return Err(fail(MultilevelError::Budget("hierarchy levels/fine/bottom limits")));
        }
        let mut transfers = Vec::with_capacity(coarser.len());
        let mut scalar_terms = 0usize; let mut vector_terms = 0usize; let mut expected = fine.n();
        for (level, &coarse) in coarser.iter().enumerate() {
            let source = if level == 0 { fine } else { coarser[level-1] };
            let transfer = AdaptiveTransfer3::new(coarse, source,
                options.max_transfer_terms.saturating_sub(scalar_terms), &mut checkpoint)
                .map_err(AdaptivePreconditionError3::Physics)?;
            let terms = transfer.rows.iter().try_fold(0usize, |n, r| n.checked_add(r.len()))
                .ok_or_else(|| fail(MultilevelError::Budget("scalar hierarchy transfer")))?;
            scalar_terms = scalar_terms.checked_add(terms).filter(|n| *n <= options.max_transfer_terms)
                .ok_or_else(|| fail(MultilevelError::Budget("scalar hierarchy transfer")))?;
            // Only the first transfer has the full fine clamp rows. Every
            // intermediate operator lives in COMPACT unconstrained coordinates.
            let p = transfer.vector_matrix(level > 0, budget.max_fine_dofs,
                budget.max_fine_dofs, budget.max_transfer_entries.saturating_sub(vector_terms), &mut checkpoint)?;
            if p.nrows() != expected || p.ncols() >= expected {
                return Err(fail(MultilevelError::Invalid("coarse geometries must strictly decrease in order")));
            }
            vector_terms = vector_terms.checked_add(p.nnz()).filter(|n| *n <= budget.max_transfer_entries)
                .ok_or_else(|| fail(MultilevelError::Budget("vector hierarchy transfer")))?;
            expected = p.ncols(); transfers.push(p);
        }
        if expected > budget.max_coarsest_dofs {
            return Err(fail(MultilevelError::Budget("bottom geometry too large; supply another coarse level")));
        }
        poll(&mut checkpoint).map_err(AdaptivePreconditionError3::Physics)?;
        Ok(Self { transfers, options })
    }

    pub(super) fn prepare<'a>(&self, op: &'a AdaptiveElasticity3, inverse: &[f64],
        mut checkpoint: impl FnMut(AdaptiveSetupWork3) -> ControlFlow<()>)
        -> Result<SparseMultilevel<'a, AdaptiveElasticity3>, AdaptivePreconditionError3> {
        let mut progress = |w: fs_solver::op::multilevel::MultilevelWork| checkpoint(AdaptiveSetupWork3 {
            operator_applications: 0, galerkin_products: w.galerkin_products,
        });
        let mut control = MultilevelControl::new(self.options.hierarchy, &mut progress)
            .map_err(AdaptivePreconditionError3::Hierarchy)?;
        let outcome = (|| {
            let first = first_galerkin(op, &self.transfers[0], self.options.max_transfer_terms, &mut control)?;
            SparseMultilevel::new(op, inverse, self.transfers.clone(), first, &mut control)
        })();
        // A budget may end between polling batches. Publish the final spent
        // count even when the in-progress matrix is discarded.
        let final_poll = control.poll();
        match outcome {
            Err(error) => Err(AdaptivePreconditionError3::Hierarchy(error)),
            Ok(prepared) => { final_poll.map_err(AdaptivePreconditionError3::Hierarchy)?; Ok(prepared) }
        }
    }
}

// X = T P, where T is the ACTUAL hanging-node reconstruction with homogeneous
// clamps removed. Scalar component separation is preserved by vector_matrix.
// Assemble X^T K_bulk X plus each density-scaled ghost outer product. This is
// the same first Galerkin form as fine probing, but without n_coarse fine applies
// or any dense n_coarse^2 storage. No separately integrated coarse K is used.
fn first_galerkin(op: &AdaptiveElasticity3, p: &Csr, max_terms: usize,
    control: &mut MultilevelControl<'_>) -> Result<Csr, MultilevelError> {
    control.poll()?;
    if p.nrows() != op.n() || p.ncols()%3 != 0 { return Err(MultilevelError::Invalid("vector Galerkin interpolation shape")); }
    let mut physical = Vec::with_capacity(op.rows.len()); let mut terms = 0usize;
    for row in &op.rows {
        control.poll()?;
        let mut combined = BTreeMap::new();
        for &(m, weight) in row {
            let (cols, values) = p.row(3*m);
            for (&column, &value) in cols.iter().zip(values) {
                if column%3 != 0 { return Err(MultilevelError::Invalid("interpolation mixes vector components")); }
                let key = column/3;
                if !combined.contains_key(&key) {
                    if terms >= max_terms { return Err(MultilevelError::Budget("composed physical interpolation")); }
                    terms += 1;
                }
                let entry = combined.entry(key).or_insert(0.0); *entry += weight*value;
                if !entry.is_finite() { return Err(MultilevelError::Invalid("composed interpolation overflow")); }
            }
        }
        physical.push(combined.into_iter().collect::<Vec<(usize,f64)>>());
    }
    let mut out = SymmetricGalerkin::new(p.ncols(), control)?;
    for (cell_id, cell) in op.raw.cells.iter().enumerate() {
        control.poll()?;
        for a in 0..8 { for b in 0..8 {
            control.poll()?;
            for &(m, wa) in &physical[cell.nodes[a]] { for &(n, wb) in &physical[cell.nodes[b]] {
                for c in 0..3 { for d in 0..3 {
                    let (i,j) = (3*m+c,3*n+d);
                    if i <= j {
                        let value = ((op.scales()[cell_id]*cell.stiffness[3*a+c][3*b+d])*wa)*wb;
                        out.add(i,j,value,control)?;
                    }
                } }
            } }
        } }
    }
    for face in &op.raw.ghosts {
        control.poll()?;
        let mut jump = BTreeMap::new();
        for (&node, &derivative) in face.nodes.iter().zip(&face.jump) {
            for &(m, weight) in &physical[node] {
                *jump.entry(m).or_insert(0.0) += derivative*weight;
            }
        }
        let weight = face.weight*0.5*(op.scales()[face.cells[0]]+op.scales()[face.cells[1]]);
        for (&m, &jm) in &jump { for (&n, &jn) in &jump {
            if m <= n { for c in 0..3 { out.add(3*m+c,3*n+c,(weight*jm)*jn,control)?; } }
        } }
    }
    out.finish(control)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{CutSdf3,HeightAxis};
    struct Slab;
    impl CutSdf3 for Slab {
        fn value(&self,p:[f64;3])->f64 {p[2]-0.73}
        fn enclose(&self,lo:[f64;3],hi:[f64;3])->Interval {Interval::new(lo[2],hi[2])-Interval::new(0.73,0.73)}
        fn derivative_enclose(&self,_:[f64;3],_:[f64;3],a:HeightAxis)->Interval {
            let d=if a==HeightAxis::Z{1.0}else{0.0};Interval::new(d,d)
        }
    }
    fn build(tree:&crate::octree3::Octree3)->AdaptiveElasticity3 {
        let mut cp=|_|ControlFlow::Continue(());
        let mut q=QuadratureControl3::new(crate::quad3::QuadratureOptions3{depth:1,..Default::default()},&mut cp).unwrap();
        AdaptiveElasticity3::build(HexCell::try_new([0.0;3],[1.0;3]).unwrap(),tree,&Slab,
            &IsotropicElastic::new(1.0,0.3,1.0).unwrap(),&|p|p[0]==0.0,ElasticityOptions3::default(),&mut q).unwrap()
    }
    #[test]
    fn g0_local_galerkin_matches_independent_full_operator_probing_on_hanging_cuts() {
        let tree=crate::octree3::Octree3::uniform(1,4,4096).unwrap();
        let refined=tree.refined(&[*tree.leaves().iter().next().unwrap()],||ControlFlow::Continue(())).unwrap();
        let coarse=build(&tree);let mut fine=build(&refined);
        let scales:Vec<f64>=(0..fine.cells()).map(|i|0.2+0.1*(i%7) as f64).collect();fine.set_scales(&scales).unwrap();
        let h=AdaptiveHierarchy3::new(&fine,&[&coarse],AdaptiveMultilevelOptions3::default(),||ControlFlow::Continue(())).unwrap();
        let mut callback=|_|ControlFlow::Continue(());let mut c=MultilevelControl::new(MultilevelBudget::default(),&mut callback).unwrap();
        let a=first_galerkin(&fine,&h.transfers[0],4_000_000,&mut c).unwrap();let p=&h.transfers[0];let pt=fs_sparse::ops::transpose(p);
        for j in 0..a.ncols() {
            let mut unit=vec![0.0;a.ncols()];unit[j]=1.0;let mut x=vec![0.0;fine.n()];p.spmv(&unit,&mut x);
            let mut ax=vec![0.0;fine.n()];fine.apply(&x,&mut ax);let mut expected=vec![0.0;a.ncols()];pt.spmv(&ax,&mut expected);
            let scale=expected.iter().map(|v|v.abs()).fold(1e-30_f64,f64::max);
            for (i,&v) in expected.iter().enumerate(){assert!((a.get(i,j)-v).abs()<1e-12*scale,"({i},{j})");}
        }
        assert!(c.work().galerkin_products>0);
    }
}
