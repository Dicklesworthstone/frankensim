//! Retained nonzero dependency rows for the existing pHS residual.
//!
//! This is a traversal plan, not a replacement matrix or an approximate sparse
//! operator. Read the current J-R before each held-input step. Never discard a
//! small but nonzero coupling; Newton and ledger arithmetic still read the owner.
use crate::{PhsError, PortHamiltonian};
use super::{dimensions, PreparedStepError};

#[path = "port_load.rs"]
mod port_load;

#[derive(Debug)]
struct ColumnPattern {
    starts: Vec<usize>,
    rows: Vec<usize>,
    next: Vec<usize>,
    ready: bool,
}

#[derive(Debug)]
pub(super) struct FlowPattern {
    n: usize,
    rows: Vec<usize>,
    columns: Vec<usize>,
    transpose: Option<ColumnPattern>,
}
impl FlowPattern {
    pub(super) fn new(n: usize) -> Result<Self, PhsError> {
        let square = n.checked_mul(n).ok_or_else(|| dimensions("flow pattern extent"))?;
        let row_count = n.checked_add(1).ok_or_else(|| dimensions("flow row extent"))?;
        let mut rows = Vec::new();
        rows.try_reserve_exact(row_count).map_err(|_| dimensions("flow row capacity"))?;
        rows.resize(row_count, 0);
        let mut columns = Vec::new();
        // A caller may switch to a dense operator of the same dimension. Reserve
        // for that at preparation, not while audio is being advanced.
        columns.try_reserve_exact(square).map_err(|_| dimensions("flow pattern capacity"))?;
        Ok(Self { n, rows, columns, transpose: None })
    }
    /// Allocate the alternate traversal only when a condensed analytic image
    /// selects it. No coefficients are copied, and no hot allocation is needed
    /// if a later same-size operator becomes completely dense.
    pub(super) fn set_column_traversal(&mut self, enabled: bool) -> Result<(), PhsError> {
        if enabled && self.transpose.is_none() {
            let mut starts=Vec::new();
            starts.try_reserve_exact(self.n+1).map_err(|_|dimensions("column pattern starts capacity"))?;
            starts.resize(self.n+1,0);
            let mut rows=Vec::new();
            rows.try_reserve_exact(self.n*self.n).map_err(|_|dimensions("column pattern rows capacity"))?;
            let mut next=Vec::new();
            next.try_reserve_exact(self.n).map_err(|_|dimensions("column pattern cursor capacity"))?;
            next.resize(self.n,0);
            self.transpose=Some(ColumnPattern {starts,rows,next,ready:false});
        } else if !enabled { self.transpose=None; }
        Ok(())
    }
    #[cfg(test)]
    pub(super) fn refresh(&mut self, sys: &PortHamiltonian) -> Result<(), PhsError> {
        self.refresh_for(sys,true)
    }
    pub(super) fn refresh_for(&mut self, sys:&PortHamiltonian, columns:bool) -> Result<(),PhsError> {
        if let Some(pattern)=&mut self.transpose {pattern.ready=false;}
        if sys.n != self.n || sys.j.len() != self.n*self.n || sys.r.len() != sys.j.len() {
            return Err(dimensions("flow pattern structure dimensions"));
        }
        self.columns.clear();
        for row in 0..self.n {
            self.rows[row] = self.columns.len();
            for col in 0..self.n {
                let value = sys.j[row*self.n+col] - sys.r[row*self.n+col];
                if !value.is_finite() { return Err(dimensions("nonfinite J-R flow coefficient")); }
                if value != 0.0 { self.columns.push(col); }
            }
        }
        self.rows[self.n] = self.columns.len();
        if let Some(pattern)=self.transpose.as_mut().filter(|_|columns) {
            pattern.starts.fill(0);
            for &col in &self.columns { pattern.starts[col+1]+=1; }
            for i in 0..self.n { pattern.starts[i+1]+=pattern.starts[i]; }
            pattern.next.copy_from_slice(&pattern.starts[..self.n]);
            pattern.rows.resize(self.columns.len(),0); // capacity reserved cold
            for row in 0..self.n {
                for &col in &self.columns[self.rows[row]..self.rows[row+1]] {
                    pattern.rows[pattern.next[col]]=row; pattern.next[col]+=1;
                }
            }
            pattern.ready=true;
        }
        Ok(())
    }
    pub(super) fn row(&self, row: usize) -> &[usize] {
        &self.columns[self.rows[row]..self.rows[row+1]]
    }
    /// Apply J-R to a possibly sparse direction. Visit source indices in the
    /// SAME increasing order as the row traversal, but skip only EXACT zero
    /// effort components. No magnitude threshold, stale coefficients or rank
    /// truncation. Returns the number of actual scalar products evaluated.
    /// Both input and output are solver scratch; the caller owns publication.
    pub(super) fn apply_direction<F>(&self, sys:&PortHamiltonian, direction:&[f64],
        out:&mut[f64], poll:&mut F) -> Result<usize,PreparedStepError>
    where F:FnMut()->Result<(),PreparedStepError> {
        let pattern=self.transpose.as_ref().ok_or_else(||dimensions("column traversal is not prepared"))?;
        if !pattern.ready || direction.len()!=self.n || out.len()!=self.n || sys.n!=self.n {
            return Err(dimensions("column traversal dimensions").into());
        }
        super::norm(direction)?;
        out.fill(0.0); let mut products=0usize;
        for (col,&effort) in direction.iter().enumerate() {
            if col%64==0 { poll()?; }
            if effort==0.0 { continue; }
            for &row in &pattern.rows[pattern.starts[col]..pattern.starts[col+1]] {
                out[row]+=(sys.j[row*self.n+col]-sys.r[row*self.n+col])*effort;
                products+=1;
            }
        }
        super::norm(out)?; Ok(products)
    }

}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{QuadraticStorage, RelaxationBranch, modal_bank};
    fn raw(n: usize, j: Vec<f64>, r: Vec<f64>) -> PortHamiltonian {
        PortHamiltonian::from_raw_parts(n, 0, j, r, vec![],
            Box::new(QuadraticStorage::new(vec![0.0;n*n],n).unwrap()))
    }
    #[test]
    fn canonical_modal_and_memory_rows_do_not_visit_dense_structural_zeros() {
        let count=32;
        let mut sys=modal_bank(&vec![1000.0;count],&vec![0.01;count],&vec![1.0;count]).unwrap();
        let branches=(0..6).map(|i| {let mut projection=vec![0.0;2*count];projection[2*i]=1.0;
            RelaxationBranch {projection,stiffness:1.0,relaxation_time_s:0.01}}).collect();
        sys=sys.with_relaxation_branches(branches).unwrap();
        let mut pattern=FlowPattern::new(sys.n).unwrap();pattern.refresh(&sys).unwrap();
        assert_eq!(sys.n,70);assert_eq!(pattern.columns.len(),102);
        assert_eq!(sys.n*sys.n,4900);
        for i in 0..count {
            assert_eq!(pattern.row(2*i),&[2*i+1]);
            assert_eq!(pattern.row(2*i+1),&[2*i,2*i+1]);
        }
        for i in 2*count..sys.n {assert_eq!(pattern.row(i),&[i]);}
    }
    #[test]
    fn finite_dense_and_retained_flow_have_identical_accumulation_bits() {
        let n=7;let mut j=vec![0.0;n*n];let mut r=j.clone();
        for row in 0..n {for col in 0..n {
            if (row+2*col)%3==0 {j[row*n+col]=(row as f64-col as f64)*0.125;}
            if (row+col)%4==0 {r[row*n+col]=0.25;}
        }}
        // Exact cancellation and subnormal couplings are both intentional.
        j[0]=r[0];j[3]=1e-310;r[3]=0.0;j[4]=-0.0;r[4]=0.0;
        let sys=raw(n,j,r);let mut pattern=FlowPattern::new(n).unwrap();pattern.refresh(&sys).unwrap();
        assert!(!pattern.row(0).contains(&0));assert!(pattern.row(0).contains(&3));
        for seed in 0..32 {let effort:Vec<_>=(0..n).map(|i| {
            if (i+seed)%5==0 {-0.0} else {((i+1)*(seed+1)) as f64*if i%2==0 {1e-80} else {-1e80}}
        }).collect();
            for row in 0..n {
                let mut dense=0.0;for col in 0..n {dense+=(sys.j[row*n+col]-sys.r[row*n+col])*effort[col];}
                let mut retained=0.0;for &col in pattern.row(row) {retained+=(sys.j[row*n+col]-sys.r[row*n+col])*effort[col];}
                assert_eq!(dense.to_bits(),retained.to_bits());
            }
        }
    }
    #[test]
    fn changing_operator_topology_reuses_capacity_without_stale_couplings() {
        let mut pattern=FlowPattern::new(3).unwrap();let capacity=pattern.columns.capacity();
        pattern.refresh(&raw(3,vec![0.0;9],vec![0.0;9])).unwrap();
        assert!(pattern.columns.is_empty());
        let dense=raw(3,vec![1.0;9],vec![0.0;9]);pattern.refresh(&dense).unwrap();
        assert_eq!(pattern.columns.len(),9);assert_eq!(pattern.columns.capacity(),capacity);
        let sparse=raw(3,vec![0.0,1.0,0.0, -1.0,0.0,0.0, 0.0,0.0,0.0],vec![0.0;9]);
        pattern.refresh(&sparse).unwrap();assert_eq!(pattern.row(0),&[1]);
        assert_eq!(pattern.row(1),&[0]);assert!(pattern.row(2).is_empty());
        assert_eq!(pattern.columns.capacity(),capacity);
    }
    #[test]
    fn empty_extent_and_overflow_refusals_remain_reusable() {
        let mut empty=FlowPattern::new(0).unwrap();empty.refresh(&raw(0,vec![],vec![])).unwrap();
        assert_eq!(empty.rows.as_slice(), &[0]);assert!(FlowPattern::new(usize::MAX).is_err());
        let mut pattern=FlowPattern::new(1).unwrap();
        for j in [f64::NAN,f64::INFINITY] {assert!(pattern.refresh(&raw(1,vec![j],vec![0.0])).is_err());}
        assert!(pattern.refresh(&raw(1,vec![f64::MAX],vec![-f64::MAX])).is_err());
        pattern.refresh(&raw(1,vec![1.0],vec![0.0])).unwrap();assert_eq!(pattern.row(0),&[0]);
    }
    #[test]
    fn sparse_directions_preserve_order_bits_subnormals_and_refreshed_topology() {
        let n=9; let mut j=vec![0.0;n*n]; let mut r=j.clone();
        for row in 0..n {for col in 0..n {
            if (row+col)%3==0 {j[row*n+col]=(row as f64-col as f64)*0.125;}
            if (row+2*col)%5==0 {r[row*n+col]=0.25;}
        }}
        j[3]=f64::from_bits(1); r[3]=0.0;
        let mut sys=raw(n,j,r); let mut pattern=FlowPattern::new(n).unwrap();
        assert!(pattern.transpose.is_none()); pattern.set_column_traversal(true).unwrap();
        let capacity=pattern.transpose.as_ref().unwrap().rows.capacity();
        for round in 0..3 {
            if round==1 {sys.j.fill(0.125);sys.r.fill(0.0);}
            if round==2 {sys.j.fill(0.0);sys.r.fill(0.0);sys.j[3]=f64::from_bits(1);}
            pattern.refresh(&sys).unwrap();
            assert_eq!(pattern.transpose.as_ref().unwrap().rows.capacity(),capacity);
            for seed in 0..16 {
                let direction:Vec<_>=(0..n).map(|i|if (i+seed)%3==0 {-0.0} else {((i+seed)%5) as f64-2.0}).collect();
                let mut actual=vec![0.0;n];
                let products=pattern.apply_direction(&sys,&direction,&mut actual,&mut ||Ok(())).unwrap();
                let mut expected_products=0;
                for row in 0..n {
                    let mut expected=0.0;
                    for &col in pattern.row(row) {
                        expected+=(sys.j[row*n+col]-sys.r[row*n+col])*direction[col];
                        expected_products+=usize::from(direction[col]!=0.0);
                    }
                    assert_eq!(actual[row].to_bits(),expected.to_bits());
                }
                assert_eq!(products,expected_products);
            }
        }
        let mut direction=vec![0.0;n];direction[3]=1.0;let mut out=vec![0.0;n];
        pattern.apply_direction(&sys,&direction,&mut out,&mut ||Ok(())).unwrap();
        assert_eq!(out[0].to_bits(),1,"the smallest coefficient is not numerical sparsity");
        assert_eq!(pattern.apply_direction(&sys,&direction,&mut out,&mut ||Err(PreparedStepError::Cancelled)).unwrap_err(),
            PreparedStepError::Cancelled);
        pattern.set_column_traversal(false).unwrap();assert!(pattern.transpose.is_none());
        assert!(pattern.apply_direction(&sys,&direction,&mut out,&mut ||Ok(())).is_err());
    }

}
