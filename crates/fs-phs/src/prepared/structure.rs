//! Retained nonzero dependency rows for the existing pHS residual.
//!
//! This is a traversal plan, not a replacement matrix or an approximate sparse
//! operator. Read the current J-R before each held-input step. Never discard a
//! small but nonzero coupling; Newton and ledger arithmetic still read the owner.
use crate::{PhsError, PortHamiltonian};
use super::dimensions;

#[path = "port_load.rs"]
mod port_load;

#[derive(Debug)]
pub(super) struct FlowPattern {
    n: usize,
    rows: Vec<usize>,
    columns: Vec<usize>,
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
        Ok(Self { n, rows, columns })
    }
    pub(super) fn refresh(&mut self, sys: &PortHamiltonian) -> Result<(), PhsError> {
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
        Ok(())
    }
    pub(super) fn row(&self, row: usize) -> &[usize] {
        &self.columns[self.rows[row]..self.rows[row+1]]
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
}
