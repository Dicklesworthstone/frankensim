//! Exact block elimination for a bordered Jacobian plus one rank-one update.
//!
//! Gonzalez's energy correction couples otherwise independent linear states.
//! Keep it as an extra scalar border, rather than discarding it or relying on
//! an invertible full uncorrected Jacobian (as Sherman--Morrison would require).
//! Every eliminated block is solved with the existing partial-pivot LU owner.
use super::{PreparedStepError, dimensions, norm, zeroed};
use crate::PhsError;
use fs_la::LuWorkspace;

#[derive(Debug)]
pub(super) struct Condensation {
    n: usize,
    pairs: Vec<[usize; 2]>,
    retained: Vec<usize>,
    owner: Vec<usize>,
    pub(super) left: Vec<f64>,
    pub(super) right: Vec<f64>,
    matrix: Vec<f64>,
    rhs: Vec<f64>,
    solution: Vec<f64>,
    responses: Vec<f64>,
    candidate: Vec<f64>,
    leaf_lu: LuWorkspace,
    border_lu: LuWorkspace,
}

impl Condensation {
    pub(super) fn new(n: usize, pairs: &[[usize; 2]]) -> Result<Self, PhsError> {
        if pairs.is_empty() || pairs.len().checked_mul(2).is_none_or(|size| size > n) {
            return Err(dimensions("condensed Newton pair count"));
        }
        // Reserve through the same fallible owner as the rest of StepWorkspace.
        let mut owner = Vec::new();
        owner.try_reserve_exact(n).map_err(|_| dimensions("condensed Newton index capacity"))?;
        owner.resize(n, usize::MAX);
        for (p, pair) in pairs.iter().enumerate() {
            for &index in pair {
                if index >= n || owner[index] != usize::MAX {
                    return Err(dimensions("condensed Newton needs disjoint in-range pairs"));
                }
                owner[index] = p;
            }
        }
        let retained: Vec<_> = (0..n).filter(|&i| owner[i] == usize::MAX).collect();
        let k = retained.len().checked_add(1).ok_or_else(|| dimensions("condensed border extent"))?;
        let entries = k.checked_mul(k).ok_or_else(|| dimensions("condensed matrix extent"))?;
        let response_size = (k+1).checked_mul(2*pairs.len())
            .ok_or_else(|| dimensions("condensed responses extent"))?;
        Ok(Self {
            n, pairs: pairs.to_vec(), retained, owner,
            left: zeroed(n)?, right: zeroed(n)?, matrix: zeroed(entries)?,
            rhs: zeroed(k)?, solution: zeroed(k)?, responses: zeroed(response_size)?,
            candidate: zeroed(n)?,
            leaf_lu: LuWorkspace::new(2).map_err(|_| dimensions("condensed leaf LU capacity"))?,
            border_lu: LuWorkspace::new(k).map_err(|_| dimensions("condensed border LU capacity"))?,
        })
    }
    pub(super) fn dimension(&self) -> usize { self.retained.len()+1 }

    /// Solve (base + left * right^T) x = rhs. Return false, without publication,
    /// when exact pair sparsity, finite elimination or the backward-error check
    /// fails. The caller then uses its original full dense LU on the SAME matrix.
    pub(super) fn solve<F>(&mut self, base: &[f64], rhs: &[f64], out: &mut [f64],
        poll: &mut F) -> Result<bool, PreparedStepError>
    where F: FnMut() -> Result<(), PreparedStepError> {
        let n=self.n; let k=self.dimension(); let last=k-1;
        if base.len()!=n*n || rhs.len()!=n || out.len()!=n {
            return Err(dimensions("condensed Newton solve dimensions").into());
        }
        norm(base)?; norm(rhs)?; norm(&self.left)?; norm(&self.right)?;
        // No epsilon sparsification: even a subnormal off-block term triggers
        // dense fallback. Layout reuse with another model cannot drop physics.
        for pair in &self.pairs {
            poll()?;
            for &row in pair { for col in 0..n {
                if self.owner[col]!=usize::MAX && self.owner[col]!=self.owner[row]
                    && base[row*n+col]!=0.0 { return Ok(false); }
            }}
        }
        self.matrix.fill(0.0); self.rhs.fill(0.0);
        for (i,&row) in self.retained.iter().enumerate() {
            self.rhs[i]=rhs[row];
            for (j,&col) in self.retained.iter().enumerate() { self.matrix[i*k+j]=base[row*n+col]; }
            self.matrix[i*k+last]=self.left[row]; self.matrix[last*k+i]=-self.right[row];
        }
        self.matrix[last*k+last]=1.0;
        // Augment with y=right^T*x; eliminate independent 2x2 leaf blocks.
        // Store D^-1[F | b] for back-substitution, not an explicit inverse.
        for (p,&[a,b]) in self.pairs.iter().enumerate() {
            poll()?;
            let d=[base[a*n+a],base[a*n+b],base[b*n+a],base[b*n+b]];
            for j in 0..=k {
                let input=if j==k { [rhs[a],rhs[b]] }
                    else if j==last { [self.left[a],self.left[b]] }
                    else { let col=self.retained[j]; [base[a*n+col],base[b*n+col]] };
                let mut response=[0.0;2];
                if self.leaf_lu.solve_into(&d,&input,&mut response).is_err() { return Ok(false); }
                for side in 0..2 { self.responses[(2*p+side)*(k+1)+j]=response[side]; }
            }
            for i in 0..k {
                let e=if i==last { [-self.right[a],-self.right[b]] }
                    else { let row=self.retained[i]; [base[row*n+a],base[row*n+b]] };
                for j in 0..k {
                    self.matrix[i*k+j]-=e[0]*self.responses[2*p*(k+1)+j]
                        +e[1]*self.responses[(2*p+1)*(k+1)+j];
                }
                self.rhs[i]-=e[0]*self.responses[2*p*(k+1)+k]
                    +e[1]*self.responses[(2*p+1)*(k+1)+k];
            }
        }
        poll()?;
        if self.border_lu.solve_into(&self.matrix,&self.rhs,&mut self.solution).is_err() { return Ok(false); }
        for (i,&row) in self.retained.iter().enumerate() { self.candidate[row]=self.solution[i]; }
        for (p,pair) in self.pairs.iter().enumerate() {
            poll()?;
            for (side,&row) in pair.iter().enumerate() {
                let start=(2*p+side)*(k+1);
                let mut value=self.responses[start+k];
                for j in 0..k { value-=self.responses[start+j]*self.solution[j]; }
                self.candidate[row]=value;
            }
        }
        if self.candidate.iter().any(|x|!x.is_finite()) { return Ok(false); }
        // Test the ORIGINAL full equations, not only the Schur system. Failure
        // does not relax Newton or energy admission: retry full partial-pivot LU.
        let dot=self.right.iter().zip(&self.candidate).map(|(v,x)|v*x).sum::<f64>();
        let dot_abs=self.right.iter().zip(&self.candidate).map(|(v,x)|(v*x).abs()).sum::<f64>();
        if !dot.is_finite() || !dot_abs.is_finite() { return Ok(false); }
        for row in 0..n {
            poll()?;
            let mut value=self.left[row]*dot;
            let mut magnitude=rhs[row].abs()+self.left[row].abs()*dot_abs;
            for col in 0..n {
                let term=base[row*n+col]*self.candidate[col]; value+=term; magnitude+=term.abs();
            }
            let tolerance=128.0*f64::EPSILON*(n+1) as f64*magnitude;
            if !value.is_finite() || !tolerance.is_finite() || (value-rhs[row]).abs()>tolerance {
                return Ok(false);
            }
        }
        poll()?;
        out.copy_from_slice(&self.candidate);
        Ok(true)
    }
}

#[cfg(test)]
#[path="condensed_tests.rs"]
mod tests;
