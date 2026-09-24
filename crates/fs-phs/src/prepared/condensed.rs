//! Exact block elimination for a bordered Jacobian plus one rank-one update.
//!
//! Gonzalez's energy correction stays in an extra scalar border. Factor each
//! independent pair and the equilibrated border once, then reuse those factors
//! for bounded iterative refinement against the ORIGINAL full equation. Neither
//! physical coordinates nor Newton/energy tolerances are rescaled or relaxed.
use super::{PreparedStepError, dimensions, norm, zeroed};
use crate::PhsError;
use fs_la::LuWorkspace;

const MAX_REFINEMENTS: usize = 2;

#[derive(Debug)]
pub(super) struct Condensation {
    n: usize,
    pairs: Vec<[usize; 2]>,
    retained: Vec<usize>,
    owner: Vec<usize>,
    pub(super) left: Vec<f64>,
    pub(super) right: Vec<f64>,
    scaled_left: Vec<f64>,
    scaled_right: Vec<f64>,
    matrix: Vec<f64>,
    rhs: Vec<f64>,
    solution: Vec<f64>,
    responses: Vec<f64>,
    candidate: Vec<f64>,
    correction: Vec<f64>,
    error: Vec<f64>,
    row_scale: Vec<f64>,
    col_scale: Vec<f64>,
    leaf_lu: Vec<LuWorkspace>,
    border_lu: LuWorkspace,
    refinements: usize,
}

// Exact binary scaling, including values whose exponent is subnormal. The
// reciprocal scaling changes only the auxiliary unknown; it never rescales
// physical state, storage, force, or the caller's original rank-one factors.
fn exponent(x: f64) -> i32 {
    let bits=x.abs().to_bits(); let field=(bits>>52) as i32;
    if field==0 { -1074+63-bits.leading_zeros() as i32 } else { field-1023 }
}
fn scale_binary(mut x: f64, mut shift: i32) -> f64 {
    // At most two factors are needed: finite f64 exponents span 2098 values,
    // and the chosen shift is half the difference of two such exponents.
    while shift>1023 { x*=f64::from_bits(2046_u64<<52); shift-=1023; }
    while shift< -1022 { x*=f64::MIN_POSITIVE; shift+=1022; }
    x*f64::from_bits(((shift+1023) as u64)<<52)
}

impl Condensation {
    pub(super) fn new(n: usize, pairs: &[[usize; 2]]) -> Result<Self, PhsError> {
        if pairs.is_empty() || pairs.len().checked_mul(2).is_none_or(|size| size > n) {
            return Err(dimensions("condensed Newton pair count"));
        }
        let mut owner = Vec::new();
        owner.try_reserve_exact(n).map_err(|_| dimensions("condensed Newton index capacity"))?;
        owner.resize(n, usize::MAX);
        let mut leaf_lu = Vec::new();
        leaf_lu.try_reserve_exact(pairs.len()).map_err(|_| dimensions("condensed leaf capacity"))?;
        for (p, pair) in pairs.iter().enumerate() {
            for &index in pair {
                if index >= n || owner[index] != usize::MAX {
                    return Err(dimensions("condensed Newton needs disjoint in-range pairs"));
                }
                owner[index] = p;
            }
            leaf_lu.push(LuWorkspace::new(2).map_err(|_| dimensions("condensed leaf LU capacity"))?);
        }
        let retained: Vec<_> = (0..n).filter(|&i| owner[i] == usize::MAX).collect();
        let k = retained.len().checked_add(1).ok_or_else(|| dimensions("condensed border extent"))?;
        let entries = k.checked_mul(k).ok_or_else(|| dimensions("condensed matrix extent"))?;
        let response_size = k.checked_add(1).and_then(|v| v.checked_mul(2*pairs.len()))
            .ok_or_else(|| dimensions("condensed responses extent"))?;
        Ok(Self {
            n, pairs: pairs.to_vec(), retained, owner,
            left: zeroed(n)?, right: zeroed(n)?, scaled_left: zeroed(n)?, scaled_right: zeroed(n)?,
            matrix: zeroed(entries)?,
            rhs: zeroed(k)?, solution: zeroed(k)?, responses: zeroed(response_size)?,
            candidate: zeroed(n)?, correction: zeroed(n)?, error: zeroed(n)?,
            row_scale: zeroed(k)?, col_scale: zeroed(k)?, leaf_lu,
            border_lu: LuWorkspace::new(k).map_err(|_| dimensions("condensed border LU capacity"))?,
            refinements: 0,
        })
    }
    pub(super) fn dimension(&self) -> usize { self.retained.len()+1 }

    fn balance(&mut self) {
        let left=self.left.iter().fold(0.0_f64,|s,x|s.max(x.abs()));
        let right=self.right.iter().fold(0.0_f64,|s,x|s.max(x.abs()));
        if left==0.0 || right==0.0 {
            // The original outer product is exactly zero. A dummy y=0 avoids
            // imposing an arbitrarily scaled, irrelevant auxiliary equation.
            self.scaled_left.fill(0.0); self.scaled_right.fill(0.0); return;
        }
        let shift=(exponent(right)-exponent(left))/2;
        let mut reversible=true;
        for i in 0..self.n {
            self.scaled_left[i]=scale_binary(self.left[i],shift);
            self.scaled_right[i]=scale_binary(self.right[i],-shift);
            reversible &= self.scaled_left[i].is_finite() && self.scaled_right[i].is_finite()
                && scale_binary(self.scaled_left[i],-shift).to_bits()==self.left[i].to_bits()
                && scale_binary(self.scaled_right[i],shift).to_bits()==self.right[i].to_bits();
        }
        // Do not discard an underflowed coefficient for better conditioning.
        if !reversible {
            self.scaled_left.copy_from_slice(&self.left); self.scaled_right.copy_from_slice(&self.right);
        }
    }

    /// Solve (base + left * right^T) x = rhs. Return false without publication
    /// on unsuitable structure, failed elimination or failed ORIGINAL backward
    /// error after bounded refinement. The caller then factors the same full
    /// equation. No finite coupling is dropped and no output is partly updated.
    pub(super) fn solve<F>(&mut self, base: &[f64], rhs: &[f64], out: &mut [f64],
        poll: &mut F) -> Result<bool, PreparedStepError>
    where F: FnMut() -> Result<(), PreparedStepError> {
        self.refinements=0;
        let n=self.n;
        if base.len()!=n*n || rhs.len()!=n || out.len()!=n {
            return Err(dimensions("condensed Newton solve dimensions").into());
        }
        norm(base)?; norm(rhs)?; norm(&self.left)?; norm(&self.right)?;
        if !self.factor(base,poll)? { return Ok(false); }
        self.error.copy_from_slice(rhs);
        if !self.solve_rhs(base,poll)? { return Ok(false); }
        self.candidate.copy_from_slice(&self.correction);
        loop {
            match self.check(base,rhs,poll)? {
                None => return Ok(false),
                Some(true) => {
                    poll()?;
                    out.copy_from_slice(&self.candidate);
                    return Ok(true);
                }
                Some(false) => {}
            }
            if self.refinements==MAX_REFINEMENTS { return Ok(false); }
            if !self.solve_rhs(base,poll)? { return Ok(false); }
            for (x,dx) in self.candidate.iter_mut().zip(&self.correction) { *x+=dx; }
            self.refinements+=1;
        }
    }

    fn factor<F>(&mut self, base: &[f64], poll: &mut F) -> Result<bool, PreparedStepError>
    where F: FnMut() -> Result<(), PreparedStepError> {
        let n=self.n; let k=self.dimension(); let last=k-1;
        // Exact sparsity: a subnormal off-block coefficient still needs full LU.
        for pair in &self.pairs {
            poll()?;
            for &row in pair { for col in 0..n {
                if self.owner[col]!=usize::MAX && self.owner[col]!=self.owner[row]
                    && base[row*n+col]!=0.0 { return Ok(false); }
            }}
        }
        self.balance(); self.matrix.fill(0.0);
        for (i,&row) in self.retained.iter().enumerate() {
            for (j,&col) in self.retained.iter().enumerate() { self.matrix[i*k+j]=base[row*n+col]; }
            self.matrix[i*k+last]=self.scaled_left[row]; self.matrix[last*k+i]=-self.scaled_right[row];
        }
        self.matrix[last*k+last]=1.0;
        // Augment y=scaled_right^T*x. Keep D^-1 F, NOT a formed inverse. Reuse each
        // original pivoted factor for all border columns and every later RHS.
        for (p,&[a,b]) in self.pairs.iter().enumerate() {
            poll()?;
            let d=[base[a*n+a],base[a*n+b],base[b*n+a],base[b*n+b]];
            if self.leaf_lu[p].factor(&d).is_err() { return Ok(false); }
            for j in 0..k {
                let input=if j==last { [self.scaled_left[a],self.scaled_left[b]] }
                    else { let col=self.retained[j]; [base[a*n+col],base[b*n+col]] };
                let mut response=[0.0;2];
                if self.leaf_lu[p].solve_factored_into(&input,&mut response).is_err() { return Ok(false); }
                for side in 0..2 { self.responses[(2*p+side)*(k+1)+j]=response[side]; }
            }
            for i in 0..k {
                let e=if i==last { [-self.scaled_right[a],-self.scaled_right[b]] }
                    else { let row=self.retained[i]; [base[row*n+a],base[row*n+b]] };
                for j in 0..k {
                    let at=i*k+j;
                    self.matrix[at]=(-e[0]).mul_add(self.responses[2*p*(k+1)+j],self.matrix[at]);
                    self.matrix[at]=(-e[1]).mul_add(self.responses[(2*p+1)*(k+1)+j],self.matrix[at]);
                }
            }
        }
        // The auxiliary y can have vastly different units/magnitude from x.
        // Equilibrate ONLY this algebraic solve. Store divisors instead of their
        // inverses so tiny, representable scales do not overflow on inversion.
        for i in 0..k {
            poll()?;
            let row=&mut self.matrix[i*k..(i+1)*k];
            let scale=row.iter().fold(0.0_f64,|s,x|s.max(x.abs()));
            if scale==0.0 || !scale.is_finite() || row.iter().any(|x|!x.is_finite()) { return Ok(false); }
            self.row_scale[i]=scale;
            for x in row {
                let value=*x/scale;
                if *x!=0.0 && value==0.0 { return Ok(false); }
                *x=value;
            }
        }
        for j in 0..k {
            let scale=(0..k).fold(0.0_f64,|s,i|s.max(self.matrix[i*k+j].abs()));
            if scale==0.0 || !scale.is_finite() { return Ok(false); }
            self.col_scale[j]=scale;
            for i in 0..k { self.matrix[i*k+j]/=scale; }
        }
        poll()?;
        Ok(self.border_lu.factor(&self.matrix).is_ok())
    }

    /// Solve the current ORIGINAL-space residual in `error` into `correction`.
    /// Matrix factors and D^-1 F are unchanged. RHS scratch can be overwritten
    /// on any refusal; caller state and the accepted Newton candidate cannot.
    fn solve_rhs<F>(&mut self, base: &[f64], poll: &mut F) -> Result<bool, PreparedStepError>
    where F: FnMut() -> Result<(), PreparedStepError> {
        let n=self.n; let k=self.dimension(); let last=k-1;
        self.rhs.fill(0.0);
        for (i,&row) in self.retained.iter().enumerate() { self.rhs[i]=self.error[row]; }
        for (p,&[a,b]) in self.pairs.iter().enumerate() {
            poll()?;
            let mut response=[0.0;2];
            if self.leaf_lu[p].solve_factored_into(&[self.error[a],self.error[b]],&mut response).is_err() { return Ok(false); }
            for side in 0..2 { self.responses[(2*p+side)*(k+1)+k]=response[side]; }
            for i in 0..k {
                let e=if i==last { [-self.scaled_right[a],-self.scaled_right[b]] }
                    else { let row=self.retained[i]; [base[row*n+a],base[row*n+b]] };
                self.rhs[i]=(-e[0]).mul_add(response[0],self.rhs[i]);
                self.rhs[i]=(-e[1]).mul_add(response[1],self.rhs[i]);
            }
        }
        for (b,scale) in self.rhs.iter_mut().zip(&self.row_scale) { *b/=scale; }
        poll()?;
        if self.border_lu.solve_factored_into(&self.rhs,&mut self.solution).is_err() { return Ok(false); }
        for (x,scale) in self.solution.iter_mut().zip(&self.col_scale) { *x/=scale; }
        for (i,&row) in self.retained.iter().enumerate() { self.correction[row]=self.solution[i]; }
        for (p,pair) in self.pairs.iter().enumerate() {
            poll()?;
            for (side,&row) in pair.iter().enumerate() {
                let start=(2*p+side)*(k+1);
                let mut value=self.responses[start+k];
                for j in 0..k { value=(-self.responses[start+j]).mul_add(self.solution[j],value); }
                self.correction[row]=value;
            }
        }
        Ok(self.correction.iter().all(|x|x.is_finite()))
    }

    /// Original componentwise backward-error gate, with an FMA residual for
    /// refinement. No normwise replacement, tolerance enlargement or diagonal
    /// regularization can hide a failed physical-coordinate equation.
    fn check<F>(&mut self, base: &[f64], rhs: &[f64], poll: &mut F)
        -> Result<Option<bool>, PreparedStepError>
    where F: FnMut() -> Result<(), PreparedStepError> {
        let n=self.n;
        if self.candidate.iter().any(|x|!x.is_finite()) { return Ok(None); }
        // Check the ORIGINAL floating-point matrix used by dense fallback.
        // Form each outer-product entry before applying it: right^T*x may
        // overflow even when every matrix coefficient and residual is finite.
        let mut accepted=true;
        for row in 0..n {
            poll()?;
            let mut residual=rhs[row]; let mut magnitude=rhs[row].abs();
            for col in 0..n {
                let a=base[row*n+col]; let x=self.candidate[col];
                let outer=self.left[row]*self.right[col];
                residual=(-(a+outer)).mul_add(x,residual);
                magnitude+=(a*x).abs()+(outer*x).abs();
            }
            let tolerance=128.0*f64::EPSILON*(n+1) as f64*magnitude;
            if !residual.is_finite() || !tolerance.is_finite() { return Ok(None); }
            self.error[row]=residual;
            accepted &= residual.abs()<=tolerance;
        }
        Ok(Some(accepted))
    }
}

#[cfg(test)]
#[path="condensed_tests.rs"]
mod tests;
