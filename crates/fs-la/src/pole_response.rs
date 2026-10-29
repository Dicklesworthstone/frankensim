//! Force residues and parameter tangents of an isolated simple nonlinear pole.
//!
//! Factor [T, sigma*T_s*v; e_p^T, 0] ONCE, with v_p=1. For RHS [f,0],
//! the last solution component times sigma is a=(w^H f)/(w^H T_s v),
//! and v*a is the residue of T(s)^-1 f at the pole. For RHS [-T_theta*v,0]
//! the same solve gives (v_theta,s_theta), with (v_theta)_p=0. No symmetry,
//! invertible descriptor mass, finite differences or separate adjoint solve.
//! See Guettel & Tisseur (2017), Acta Numerica 26, sections 2.2, 2.6 and 4.2,
//! doi:10.1017/S0962492917000034. The constrained complementary solution is
//! NOT, in general, the constant Laurent coefficient of a nonlinear resolvent.
//!
//! Preparation checks the supplied eigenpair against its original coefficient
//! scale AND its Newton correction. Each solve checks the original bordered
//! rows. These are floating-point backward checks, not condition certificates,
//! mode tracking across crossings, or completeness/passivity claims. Repeated
//! and defective poles have a singular border and may refuse. One ISA replay
//! uses the existing deterministic, cancellable complex LU owner.
use crate::eigen_complex::{LuComplex, LuComplexError, lu_complex_with_cancel};
use crate::nonlinear_eig::{MAX_NONLINEAR_POLE_DIMENSION, PoleEvaluation};
use fs_math::{c64::C64, det};

/// Limits apply independently to preparation and to EACH immutable batch call.
#[derive(Debug, Clone, Copy)]
pub struct PoleResponseOptions {
    /// Positive scale in the same units as the spectral variable (also at zero).
    pub spectral_scale: f64,
    pub pole_residual_tolerance: f64,
    pub pole_step_tolerance: f64,
    /// Maximum original-row residual / (|rhs_i| + ||row_i||_1 ||solution||_inf).
    /// The homogeneous normalization constraint uses this same scale.
    pub linear_residual_tolerance: f64,
    pub max_dimension: usize,
    /// Complete RHS count per batch, in 1..=256; no column is dropped.
    pub max_vectors: usize,
    /// Conservative numerical payload including the largest returned batch,
    /// excluding caller inputs, allocator overhead and process RSS.
    pub max_dense_bytes: usize,
    /// Conservative units: 32*(n+1)^3 + 128*(max_vectors+1)*(n+1)^2.
    /// This includes factorization, correction, substitutions and row audits;
    /// it is not a measured FLOP count or wall-clock guarantee.
    pub max_work: u64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum PoleResponseError {
    Invalid { what: &'static str },
    NotPole { residual: f64, relative_correction: f64 },
    Singular { index: usize },
    LinearResidual { residual: f64 },
    Cancelled,
}
impl core::fmt::Display for PoleResponseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "FS-LA-POLE-RESPONSE: {self:?}")
    }
}
impl std::error::Error for PoleResponseError {}
impl From<LuComplexError> for PoleResponseError {
    fn from(e: LuComplexError) -> Self {
        match e {
            LuComplexError::Singular { index } => Self::Singular { index },
            LuComplexError::Cancelled => Self::Cancelled,
        }
    }
}

/// R f = lim_(s->pole) (s-pole) T(s)^-1 f, for a fixed analytic source.
/// This is one pole's contribution, not a replacement for the full response.
#[derive(Debug, Clone)]
pub struct PoleResidue {
    pub vector: Vec<C64>,
    /// Residue = prepared.right_vector() * coefficient, in its fixed gauge.
    pub coefficient: C64,
    pub linear_residual: f64,
}

/// Directional derivative in the caller's explicitly declared parameter units.
#[derive(Debug, Clone)]
pub struct PoleDerivative {
    pub value: C64,
    /// Derivative of v with the prepared pivot held at one, NOT a moving
    /// unit-norm/phase convention. Physical observers must use that same gauge.
    pub vector: Vec<C64>,
    pub linear_residual: f64,
}

/// Immutable one-pole factorization, reusable for arbitrary admitted RHS batches.
#[derive(Debug)]
pub struct PreparedPole {
    value: C64,
    vector: Vec<C64>,
    pivot: usize,
    border: Vec<C64>,
    rows: Vec<f64>,
    columns: Vec<f64>,
    lu: LuComplex,
    coefficient_scale: f64,
    options: PoleResponseOptions,
    residual: f64,
    relative_correction: f64,
}

fn invalid(what: &'static str) -> PoleResponseError { PoleResponseError::Invalid { what } }
fn finite(z: C64) -> bool { z.re.is_finite() && z.im.is_finite() }
fn divide(z: C64, s: f64) -> C64 { C64::new(z.re/s, z.im/s) }
fn poll(c: &mut impl FnMut() -> bool) -> Result<(), PoleResponseError> {
    if c() { Err(PoleResponseError::Cancelled) } else { Ok(()) }
}
fn norm(x: &[C64]) -> f64 {
    let scale = x.iter().map(|z| z.abs()).fold(0.0, f64::max);
    if scale == 0.0 { return 0.0; }
    scale * det::sqrt(x.iter().map(|&z| divide(z, scale).norm_sq()).sum())
}
fn add(sum: &mut C64, correction: &mut C64, term: C64) {
    let y = term-*correction;
    let next = *sum+y;
    *correction = (next-*sum)-y;
    *sum = next;
}
fn admit(n: usize, o: &PoleResponseOptions) -> Result<usize, PoleResponseError> {
    let d = n.checked_add(1).ok_or_else(|| invalid("border dimension overflow"))?;
    if n == 0 || n > MAX_NONLINEAR_POLE_DIMENSION || n > o.max_dimension
        || !(1..=256).contains(&o.max_vectors)
        || !o.spectral_scale.is_finite() || o.spectral_scale <= 0.0
        || [o.pole_residual_tolerance,o.pole_step_tolerance,o.linear_residual_tolerance]
            .iter().any(|t| !t.is_finite() || *t <= 0.0 || *t >= 1.0)
    { return Err(invalid("finite positive scales, tolerances and bounded dimensions required")); }
    let d = d as u128;
    let bytes = 128*d*d + 64*d*(o.max_vectors as u128+4);
    let work = 32*d*d*d + 128*(o.max_vectors as u128+1)*d*d;
    if bytes > o.max_dense_bytes as u128 || work > o.max_work as u128 {
        return Err(invalid("pole factorization/batch exceeds its explicit payload or work budget"));
    }
    Ok(d as usize)
}

impl PreparedPole {
    /// Prepare at a supplied, already resolved pole. The operator is evaluated
    /// by its domain owner, with analytic T_s and a physical coefficient scale.
    /// The vector is checked again here; trusting a stored solver receipt is
    /// insufficient. Normalization and LU scaling do not change T or T_s.
    pub fn new(
        value: C64, vector: &[C64], evaluation: &PoleEvaluation,
        options: PoleResponseOptions, mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, PoleResponseError> {
        poll(&mut cancelled)?;
        let n = vector.len();
        let d = admit(n, &options)?;
        let scale = evaluation.residual_scale;
        if !finite(value) || !value.abs().is_finite()
            || evaluation.matrix.len() != n*n || evaluation.derivative.len() != n*n
            || evaluation.matrix.iter().chain(&evaluation.derivative).chain(vector).any(|&z| !finite(z))
            || !scale.is_finite() || scale <= 0.0
        { return Err(invalid("complete finite pole, matrices and positive coefficient scale required")); }
        let mut v = vector.to_vec();
        let mut pivot = 0;
        for i in 1..n { if v[i].abs() > v[pivot].abs() { pivot=i; } }
        let amplitude = v[pivot].abs();
        if !amplitude.is_finite() || amplitude == 0.0 { return Err(invalid("nonzero eigenvector required")); }
        for x in &mut v { *x=divide(*x, amplitude); }
        let phase = v[pivot];
        for x in &mut v { *x=*x/phase; }
        v[pivot]=C64::ONE;
        let mut border = vec![C64::ZERO; d*d];
        let mut residual_rhs = vec![C64::ZERO; d];
        for i in 0..n {
            poll(&mut cancelled)?;
            let (mut r,mut rc,mut dv,mut dc)=(C64::ZERO,C64::ZERO,C64::ZERO,C64::ZERO);
            for j in 0..n {
                let a=divide(evaluation.matrix[i*n+j],scale);
                let b=divide(evaluation.derivative[i*n+j],scale).scale(options.spectral_scale);
                border[i*d+j]=a;
                add(&mut r,&mut rc,a*v[j]);
                add(&mut dv,&mut dc,b*v[j]);
            }
            residual_rhs[i]=-r;
            border[i*d+n]=dv;
        }
        border[n*d+pivot]=C64::ONE;
        if border.iter().chain(&residual_rhs).any(|&z| !finite(z)) {
            return Err(invalid("normalized pole linearization is unrepresentable"));
        }
        let residual=norm(&residual_rhs[..n])/norm(&v);
        if !residual.is_finite() || residual>options.pole_residual_tolerance {
            return Err(PoleResponseError::NotPole { residual, relative_correction:f64::INFINITY });
        }
        let mut equilibrated=border.clone();
        let mut rows=vec![0.0;d];
        let mut columns=vec![0.0;d];
        for i in 0..d {
            poll(&mut cancelled)?;
            rows[i]=equilibrated[i*d..(i+1)*d].iter().map(|z|z.abs()).fold(0.0,f64::max);
            if rows[i]==0.0 { return Err(PoleResponseError::Singular { index:i }); }
            for x in &mut equilibrated[i*d..(i+1)*d] { *x=divide(*x,rows[i]); }
        }
        for j in 0..d {
            poll(&mut cancelled)?;
            columns[j]=(0..d).map(|i|equilibrated[i*d+j].abs()).fold(0.0,f64::max);
            if columns[j]==0.0 { return Err(PoleResponseError::Singular { index:j }); }
            for i in 0..d { equilibrated[i*d+j]=divide(equilibrated[i*d+j],columns[j]); }
        }
        let lu=lu_complex_with_cancel(&equilibrated,d,&mut cancelled)?;
        let mut result=Self { value,vector:v,pivot,border,rows,columns,lu,
            coefficient_scale:scale,options,residual,relative_correction:0.0 };
        let (correction,_)=result.solve_normalized(&residual_rhs,&mut cancelled)?;
        result.relative_correction=correction[n].abs()
            *(options.spectral_scale/options.spectral_scale.max(value.abs()));
        if !result.relative_correction.is_finite() || result.relative_correction>options.pole_step_tolerance {
            return Err(PoleResponseError::NotPole { residual,relative_correction:result.relative_correction });
        }
        poll(&mut cancelled)?;
        Ok(result)
    }

    pub fn value(&self) -> C64 { self.value }
    pub fn right_vector(&self) -> &[C64] { &self.vector }
    pub fn normalization_pivot(&self) -> usize { self.pivot }
    pub fn pole_residual(&self) -> f64 { self.residual }
    pub fn relative_correction(&self) -> f64 { self.relative_correction }

    fn solve_normalized(&self, rhs: &[C64], cancelled: &mut impl FnMut()->bool)
        -> Result<(Vec<C64>,f64),PoleResponseError>
    {
        let d=self.vector.len()+1;
        let mut x:Vec<_>=rhs.iter().zip(&self.rows).map(|(&b,&s)|divide(b,s)).collect();
        if x.iter().any(|&z|!finite(z)) { return Err(invalid("border RHS overflow")); }
        self.lu.solve_with_cancel(&mut x,&mut *cancelled)?;
        for (z,&s) in x.iter_mut().zip(&self.columns) { *z=divide(*z,s); }
        if x.iter().any(|&z|!finite(z)) { return Err(invalid("border solution overflow")); }
        let solution_norm=x.iter().map(|z|z.abs()).fold(0.0,f64::max);
        let mut worst=0.0_f64;
        for i in 0..d {
            poll(cancelled)?;
            let (mut sum,mut correction)=(-rhs[i],C64::ZERO);
            let mut denominator=rhs[i].abs();
            for j in 0..d {
                add(&mut sum,&mut correction,self.border[i*d+j]*x[j]);
                denominator+=self.border[i*d+j].abs()*solution_norm;
            }
            if !finite(sum) || !denominator.is_finite() { return Err(invalid("border residual overflow")); }
            let r=if denominator==0.0 { sum.abs() } else { sum.abs()/denominator };
            if !r.is_finite() { return Err(invalid("border residual is not representable")); }
            worst=worst.max(r);
        }
        if worst>self.options.linear_residual_tolerance {
            return Err(PoleResponseError::LinearResidual { residual:worst });
        }
        poll(cancelled)?;
        Ok((x,worst))
    }

    fn admit_batch(&self, inputs:&[&[C64]]) -> Result<(),PoleResponseError> {
        if inputs.is_empty() || inputs.len()>self.options.max_vectors
            || inputs.iter().any(|x|x.len()!=self.vector.len() || x.iter().any(|&z|!finite(z)))
        { return Err(invalid("complete nonempty finite RHS batch within its declared ceiling required")); }
        Ok(())
    }

    fn action(&self, input:&[C64], sign:f64, cancelled:&mut impl FnMut()->bool)
        -> Result<(Vec<C64>,C64,f64),PoleResponseError>
    {
        poll(cancelled)?;
        let n=self.vector.len();
        let amplitude=input.iter().map(|z|z.abs()).fold(0.0,f64::max);
        if !amplitude.is_finite() { return Err(invalid("RHS amplitude overflow")); }
        if amplitude==0.0 { return Ok((vec![C64::ZERO;n],C64::ZERO,0.0)); }
        let mut rhs=vec![C64::ZERO;n+1];
        for i in 0..n { rhs[i]=divide(divide(input[i],amplitude),self.coefficient_scale).scale(sign); }
        let (mut x,residual)=self.solve_normalized(&rhs,cancelled)?;
        let coefficient=x[n].scale(self.options.spectral_scale).scale(amplitude);
        x.truncate(n);
        for value in &mut x { *value=value.scale(amplitude); }
        if !finite(coefficient) || x.iter().any(|&z|!finite(z)) {
            return Err(invalid("physical pole action is unrepresentable"));
        }
        poll(cancelled)?;
        Ok((x,coefficient,residual))
    }

    /// Laurent residues for analytic forcing columns evaluated AT the pole.
    /// Each column uses two triangular solves and a residual audit, O(n^2).
    /// Complex source phase and nonnormal left/right participation are retained.
    /// A cancelled/failed batch publishes no partial output and changes no state.
    pub fn forcing_residues(&self, forces:&[&[C64]], mut cancelled:impl FnMut()->bool)
        -> Result<Vec<PoleResidue>,PoleResponseError>
    {
        poll(&mut cancelled)?;
        self.admit_batch(forces)?;
        let mut out=Vec::with_capacity(forces.len());
        for &f in forces {
            let (_,coefficient,linear_residual)=self.action(f,1.0,&mut cancelled)?;
            let vector:Vec<_>=self.vector.iter().map(|&v|v*coefficient).collect();
            if vector.iter().any(|&z|!finite(z)) { return Err(invalid("physical pole residue overflow")); }
            out.push(PoleResidue { vector,coefficient,linear_residual });
        }
        poll(&mut cancelled)?;
        Ok(out)
    }

    /// Inputs are the analytic actions (partial T / partial theta)*v, using
    /// EXACTLY right_vector()'s fixed gauge and holding the spectral value fixed.
    /// No finite-difference step or material/geometry law is inferred here.
    pub fn parameter_derivatives(&self, actions:&[&[C64]], mut cancelled:impl FnMut()->bool)
        -> Result<Vec<PoleDerivative>,PoleResponseError>
    {
        poll(&mut cancelled)?;
        self.admit_batch(actions)?;
        let mut out=Vec::with_capacity(actions.len());
        for &action in actions {
            let (vector,value,linear_residual)=self.action(action,-1.0,&mut cancelled)?;
            out.push(PoleDerivative { value,vector,linear_residual });
        }
        poll(&mut cancelled)?;
        Ok(out)
    }
}

#[cfg(test)]
mod tests;
